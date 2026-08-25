# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# This file is included inside module Sparlectra. Do not add a module wrapper
# here (same convention as powerflow_rectangular/ and powerflow_dc/).
#
# file: src/shortcircuit/short_circuit.jl
# purpose: balanced 3-phase initial symmetrical short-circuit current
#          (Ik'' max/min) per IEC 60909-0, positive sequence only. Consumes a
#          solved-or-unsolved `Net` plus the CGMES short-circuit harvest
#          (`CGMESImporter.CGMESShortCircuitData`); the power-flow solver and
#          its Ybus assembly stay untouched (issue #277).

"""
    ShortCircuitResult

Result of one [`runShortCircuit!`](@ref) evaluation (one `case`, one or many
fault buses). Deliberately a separate type — it neither extends nor embeds
into `ACPFlowReport` (precedent: `rundcpf!` with its own report type).

# Fields
- `case::Symbol`: `:max` or `:min`.
- `base_mva::Float64`: network power base used for the per-unit assembly.
- `c_factor_override::Float64`: `0.0` = IEC 60909-0 Table 1 by voltage
  level; positive = the scalar override used for every bus.
- `rows::Vector{NamedTuple}`: one row per requested fault bus with
  `bus` (name), `bus_idx`, `vn_kV`, `island`, `status` (`:ok` /
  `:no_source` / `:isolated`), `c`, `zk_ohm`, `rx_ratio`, `ik_kA`,
  `sk_MVA`, `kappa`, `ip_kA`, `contains_defaulted_data`, `reasons`.
  Safety-flag contract: whenever a source value was substituted with a
  documented default or a contribution was skipped, the affected rows carry
  `contains_defaulted_data = true` plus the reason list — for the `:max`
  case a skipped contribution makes `ik_kA` a **lower bound**, which is the
  non-conservative direction for equipment rating.
- `messages::Vector{String}`: the run's warnings/notices (same texts that
  were logged).
"""
struct ShortCircuitResult
  case::Symbol
  base_mva::Float64
  c_factor_override::Float64
  rows::Vector{NamedTuple}
  messages::Vector{String}
end

# IEC 60909-0 Table 1 voltage factors, keyed by nominal voltage level.
# LV (Un <= 1 kV): c_max = 1.05 (the 1.10 variant for +10 % tolerance bands
# is not implemented, nor is a configurable table), c_min = 0.95.
# MV/HV (Un > 1 kV): c_max = 1.10, c_min = 1.00.
function _sc_c_factor(vn_kV::Float64, case::Symbol, override::Float64)::Float64
  override > 0.0 && return override
  if vn_kV <= 1.0
    return case === :max ? 1.05 : 0.95
  end
  return case === :max ? 1.10 : 1.00
end

# Fictitious generator resistance R_Gf per IEC 60909-0 §6.6.3 (needed so the
# R/X ratio at the fault — and with it kappa/i_p — stays meaningful even in
# machine-dominated networks): 0.05·x''d for machines >= 100 MVA, 0.07·x''d
# below that (HV machines), 0.15·x''d for LV machines.
_sc_generator_resistance(x_pu::Float64, ratedS_MVA::Float64, ratedU_kV::Float64)::Float64 = (ratedU_kV > 1.0 ? (ratedS_MVA >= 100.0 ? 0.05 : 0.07) : 0.15) * x_pu

"""
    _takahashi_diag(F) -> (diagZ::Vector{ComplexF64}, ok::Bool, info::String)

Compute the full diagonal of `inv(A)` from the UMFPACK factorization `F`
of `A` via the Takahashi/Erisman-Tinney selected-inverse recurrences on
the filled factor pattern: one backward pass over `nnz(L) + nnz(U)`
entries instead of one triangular solve per bus, the algorithmic
alternative for all-bus short-circuit sweeps (measured 34x to 264x over
the serial sweep between n = 2000 and n = 16000, agreeing with direct
solves to about 1e-15 relative).

Applicability guard: requires the symmetric pivot ordering (`F.p == F.q`)
that UMFPACK's symmetric strategy picks on the structurally symmetric
`Ysc`; with an unsymmetric ordering the original-diagonal positions leave
the factor pattern. Returns `ok = false` with the reason in `info` in that
case (or on a pattern-closure violation, counted defensively); callers
must fall back to per-bus solves.
"""
function _takahashi_diag(F)
  p = F.p
  q = F.q
  p == q || return (ComplexF64[], false, "unsymmetric pivot ordering (p != q)")
  L = F.L
  U = F.U
  Rs = F.Rs
  n = size(L, 1)
  dU = Vector{ComplexF64}(undef, n)
  for j in 1:n
    dU[j] = U[j, j]
  end
  any(iszero, dU) && return (ComplexF64[], false, "zero pivot in U")
  all(isone, diag(L)) || return (ComplexF64[], false, "L is not unit-lower")

  # row access of U and L via their CSC transposes (column i = row i)
  Ut = SparseMatrixCSC(transpose(U))
  Lt = SparseMatrixCSC(transpose(L))
  Utc = Ut.colptr
  Utr = Ut.rowval
  Utv = Ut.nzval
  Lp = L.colptr
  Lr = L.rowval
  Lv = L.nzval
  Ltc = Lt.colptr
  Ltr = Lt.rowval

  # Z entries live on the pattern of (L + U)^T; per Erisman-Tinney the
  # recurrences restricted to that (filled) pattern are self-contained
  Z = Dict{Tuple{Int,Int},ComplexF64}()
  sizehint!(Z, nnz(L) + nnz(U))
  misses = 0
  uppers = Int[]
  for j in n:-1:1
    # lower entries of column j (positions i > j with U[j, i] != 0):
    # Z[i, j] = -sum_{k > j} Z[i, k] * L[k, j]
    for t in Utc[j]:(Utc[j+1]-1)
      i = Utr[t]
      i > j || continue
      acc = zero(ComplexF64)
      for s2 in Lp[j]:(Lp[j+1]-1)
        k = Lr[s2]
        k > j || continue
        z = get(Z, (i, k), nothing)
        z === nothing ? (misses += 1) : (acc += z * Lv[s2])
      end
      Z[(i, j)] = -acc
    end
    # diagonal: Z[j, j] = 1/d_j - sum_{k > j} (U[j, k]/d_j) * Z[k, j]
    accd = zero(ComplexF64)
    for t in Utc[j]:(Utc[j+1]-1)
      k = Utr[t]
      k > j || continue
      z = get(Z, (k, j), nothing)
      z === nothing ? (misses += 1) : (accd += (Utv[t] / dU[j]) * z)
    end
    Z[(j, j)] = inv(dU[j]) - accd
    # upper entries of column j (positions i < j with L[j, i] != 0),
    # descending i: Z[i, j] = -sum_{k > i} (U[i, k]/d_i) * Z[k, j]
    empty!(uppers)
    for t in Ltc[j]:(Ltc[j+1]-1)
      Ltr[t] < j && push!(uppers, Ltr[t])
    end
    sort!(uppers; rev = true)
    for i in uppers
      acc = zero(ComplexF64)
      for t in Utc[i]:(Utc[i+1]-1)
        k = Utr[t]
        k > i || continue
        z = get(Z, (k, j), nothing)
        z === nothing ? (misses += 1) : (acc += (Utv[t] / dU[i]) * z)
      end
      Z[(i, j)] = -acc
    end
  end
  misses == 0 || return (ComplexF64[], false, "pattern closure violated at $(misses) lookups")

  # map back: with L*U == (Rs .* A)[p, q] and p == q,
  # inv(A)[i, i] = Rs[i] * Z[r, r] at r = pinv[i]
  pinv = invperm(p)
  diagZ = Vector{ComplexF64}(undef, n)
  for i in 1:n
    z = get(Z, (pinv[i], pinv[i]), nothing)
    z === nothing && return (ComplexF64[], false, "diagonal position missing from pattern")
    diagZ[i] = Rs[i] * z
  end
  return (diagZ, true, "ok")
end

# One source-admittance table for the whole net: bus index → added shunt
# admittance (pu on `base_mva`), plus per-bus and per-island safety-flag
# reasons. Islands are flagged as a whole where a skipped contribution
# affects every fault location fed from that island (motor/feeder skips).
function _sc_source_admittances(net::Net, sc, case::Symbol, c_override::Float64, bus_island::Dict{Int,Int}, messages::Vector{String})
  y_add = Dict{Int,ComplexF64}()
  bus_reasons = Dict{Int,Vector{String}}()
  island_reasons = Dict{Int,Vector{String}}()
  flag_bus = (idx, why) -> push!(get!(Vector{String}, bus_reasons, idx), why)
  flag_island = (idx, why) -> begin
    island = get(bus_island, idx, 0)
    island == 0 && return
    push!(get!(Vector{String}, island_reasons, island), why)
  end
  note = why -> (push!(messages, why); @warn why)

  sbase = net.baseMVA

  # Synchronous machines: x''d on machine base → network base. Missing data
  # is substituted with documented defaults AND flagged (safety contract) —
  # a silent default on a safety-relevant quantity is not acceptable.
  for m in sc.synchronous_machines
    m.bus === nothing && continue
    busidx = get(net.busDict, String(m.bus), nothing)
    busidx === nothing && continue
    vn = getNodeVn(net.nodeVec[busidx])
    name = something(m.name, m.mrid)
    xdpp = m.satDirectSubtransX_pu
    if xdpp === nothing || !isfinite(Float64(xdpp)) || Float64(xdpp) <= 0.0
      xdpp = 0.2
      why = "short-circuit: SynchronousMachine $(name) has no usable x''_d — default 0.2 pu (machine base) substituted"
      note(why)
      flag_bus(busidx, why)
      flag_island(busidx, why)
    end
    srated = m.ratedS_MVA
    if srated === nothing || !isfinite(Float64(srated)) || Float64(srated) <= 0.0
      srated = sbase
      why = "short-circuit: SynchronousMachine $(name) has no usable ratedS — network base $(sbase) MVA substituted"
      note(why)
      flag_bus(busidx, why)
      flag_island(busidx, why)
    end
    urated = m.ratedU_kV
    if urated === nothing || !isfinite(Float64(urated)) || Float64(urated) <= 0.0
      urated = vn
      why = "short-circuit: SynchronousMachine $(name) has no usable ratedU — bus nominal $(vn) kV substituted"
      note(why)
      flag_bus(busidx, why)
      flag_island(busidx, why)
    end
    x_pu = Float64(xdpp) * (sbase / Float64(srated)) * (Float64(urated) / vn)^2
    r_pu = _sc_generator_resistance(x_pu, Float64(srated), Float64(urated))
    y_add[busidx] = get(y_add, busidx, 0.0 + 0.0im) + 1.0 / (r_pu + im * x_pu)
  end

  # Network-feeder equivalents (ExternalNetworkInjection): Z_Q from the
  # declared initial symmetrical short-circuit current at the connection
  # point, R/X split from the declared ratio (default X = |Z|, R = 0.1·X per
  # IEC 60909-0 §6.2 when unknown — flagged).
  for f in sc.external_network_injections
    f.bus === nothing && continue
    busidx = get(net.busDict, String(f.bus), nothing)
    busidx === nothing && continue
    vn = getNodeVn(net.nodeVec[busidx])
    name = something(f.name, f.mrid)
    ik_A = case === :max ? f.maxInitialSymShCCurrent_A : f.minInitialSymShCCurrent_A
    if ik_A === nothing || !isfinite(Float64(ik_A)) || Float64(ik_A) <= 0.0
      # A feeder without a declared current cannot be modeled; in the :max
      # case the missing infeed makes every result in the island a lower
      # bound — exactly what the safety flag exists for.
      why = "short-circuit: feeder $(name) has no usable $(case === :max ? "maxInitialSymShCCurrent" : "minInitialSymShCCurrent") — contribution skipped ($(case === :max ? "Ik''max is a lower bound" : "Ik''min may be overestimated"))"
      note(why)
      flag_island(busidx, why)
      continue
    end
    c = _sc_c_factor(vn, case, c_override)
    z_ohm = c * vn * 1000.0 / (sqrt(3.0) * Float64(ik_A))
    rx = case === :max ? f.maxR1ToX1Ratio : f.minR1ToX1Ratio
    if rx === nothing || !isfinite(Float64(rx)) || Float64(rx) < 0.0
      rx = 0.1
      why = "short-circuit: feeder $(name) has no usable R/X ratio — default R = 0.1·X substituted"
      note(why)
      flag_bus(busidx, why)
      flag_island(busidx, why)
    end
    x_ohm = z_ohm / sqrt(1.0 + Float64(rx)^2)
    r_ohm = Float64(rx) * x_ohm
    z_base = vn^2 / sbase
    y_add[busidx] = get(y_add, busidx, 0.0 + 0.0im) + 1.0 / ((r_ohm + im * x_ohm) / z_base)
  end

  # Boundary equivalents (EquivalentInjection) with declared positive-sequence
  # impedance: stamped like feeders. A missing resistance defaults to 0 with a
  # flag; records without any r/x never reach the harvest.
  for e in sc.equivalent_injections
    e.bus === nothing && continue
    busidx = get(net.busDict, String(e.bus), nothing)
    busidx === nothing && continue
    x_ohm = e.x_ohm
    (x_ohm === nothing || !isfinite(Float64(x_ohm)) || Float64(x_ohm) <= 0.0) && continue
    name = something(e.name, e.mrid)
    r_ohm = e.r_ohm
    if r_ohm === nothing || !isfinite(Float64(r_ohm)) || Float64(r_ohm) < 0.0
      r_ohm = 0.0
      why = "short-circuit: EquivalentInjection $(name) has no usable r — 0 Ω substituted"
      note(why)
      flag_bus(busidx, why)
      flag_island(busidx, why)
    end
    vn = getNodeVn(net.nodeVec[busidx])
    z_base = vn^2 / sbase
    y_add[busidx] = get(y_add, busidx, 0.0 + 0.0im) + 1.0 / ((Float64(r_ohm) + im * Float64(x_ohm)) / z_base)
  end

  # Asynchronous machines: included in Ik''max iff the
  # harvested attributes suffice for the IEC 60909-0 §6.7 motor impedance
  #   Z_M = (1/(I_LR/I_rM)) · U_rM² / S_rM,   S_rM = ratedS or P_rM/(η·cosφ).
  # A missing locked-rotor R/X ratio is substituted with the §6.7.2
  # guidance values (0.10 for MV motors ≥ 1 MW per pole pair, 0.15 below,
  # 0.42 for LV motors) and flagged; a motor whose |Z_M| is not determinable
  # is skipped with the island-wide lower-bound flag. Motors never
  # contribute to Ik''min (per standard, no flag needed there).
  harvested_asm_buses = Set{Int}()
  if case === :max
    for m in sc.asynchronous_machines
      m.bus === nothing && continue
      busidx = get(net.busDict, String(m.bus), nothing)
      busidx === nothing && continue
      push!(harvested_asm_buses, busidx)
      name = something(m.name, m.mrid)
      ilr = m.iaIrRatio
      urated = m.ratedU_kV
      srated = m.ratedS_MVA
      if (srated === nothing || !isfinite(Float64(srated)) || Float64(srated) <= 0.0) &&
         m.ratedMechanicalPower_MW !== nothing && m.efficiency_percent !== nothing && m.ratedPowerFactor !== nothing &&
         Float64(m.efficiency_percent) > 0.0 && Float64(m.ratedPowerFactor) > 0.0
        srated = Float64(m.ratedMechanicalPower_MW) / ((Float64(m.efficiency_percent) / 100.0) * Float64(m.ratedPowerFactor))
      end
      missing_attr = ilr === nothing || !isfinite(Float64(ilr)) || Float64(ilr) <= 0.0 ? "iaIrRatio" :
        urated === nothing || !isfinite(Float64(urated)) || Float64(urated) <= 0.0 ? "ratedU" :
        srated === nothing || !isfinite(Float64(srated)) || Float64(srated) <= 0.0 ? "ratedS (nor ratedMechanicalPower/efficiency/ratedPowerFactor)" : nothing
      if missing_attr !== nothing
        why = "short-circuit: AsynchronousMachine $(name) has no usable $(missing_attr) — contribution skipped; Ik''max in its island is a lower bound"
        note(why)
        flag_island(busidx, why)
        continue
      end
      z_ohm_total = (1.0 / Float64(ilr)) * Float64(urated)^2 / Float64(srated)
      rx = m.rxLockedRotorRatio
      if rx === nothing || !isfinite(Float64(rx)) || Float64(rx) < 0.0
        # §6.7.2 guidance split, keyed by motor class; MW per pole pair only
        # decides the MV band and falls back to the conservative 0.15.
        mw_per_pole = m.ratedMechanicalPower_MW !== nothing && m.polePairNumber !== nothing && Float64(m.polePairNumber) > 0.0 ? Float64(m.ratedMechanicalPower_MW) / Float64(m.polePairNumber) : nothing
        rx = Float64(urated) <= 1.0 ? 0.42 : (mw_per_pole !== nothing && mw_per_pole >= 1.0 ? 0.10 : 0.15)
        why = "short-circuit: AsynchronousMachine $(name) has no usable rxLockedRotorRatio — IEC 60909-0 §6.7.2 default R/X = $(rx) substituted"
        note(why)
        flag_bus(busidx, why)
        flag_island(busidx, why)
      end
      x_ohm = z_ohm_total / sqrt(1.0 + Float64(rx)^2)
      r_ohm = Float64(rx) * x_ohm
      vn = getNodeVn(net.nodeVec[busidx])
      z_base = vn^2 / sbase
      y_add[busidx] = get(y_add, busidx, 0.0 + 0.0im) + 1.0 / ((r_ohm + im * x_ohm) / z_base)
    end
    # Motors present in the net but absent from the harvest (delivery without
    # the EquipmentShortCircuit motor attributes):
    # skip with the mandatory island-wide lower-bound flag.
    for ps in net.prosumpsVec
      ps.comp.cTyp == AsynchronousMachine || continue
      busidx = Int(ps.comp.cFrom_bus)
      busidx in eachindex(net.nodeVec) || continue
      busidx in harvested_asm_buses && continue
      why = "short-circuit: AsynchronousMachine $(ps.comp.cName) — short-circuit attributes not harvested, contribution skipped; Ik''max in its island is a lower bound"
      note(why)
      flag_island(busidx, why)
    end
  end

  return y_add, bus_reasons, island_reasons
end

# Assemble the positive-sequence short-circuit admittance matrix of one
# island: series branch impedances only (IEC 60909-0 — loads, line charging
# and shunt compensators are dropped), source admittances on the diagonal.
# Reuses the equivalent-circuit helpers (calcBranchYser/calcBranchRatio)
# instead of re-implementing the branch model; deliberately NOT a patched
# createYBUS — the PF assembly stays untouched (task scope) and the
# series-only variant needs no sparsity bookkeeping beyond this COO loop.
function _sc_island_matrix(net::Net, island_buses::Vector{Int}, y_add::Dict{Int,ComplexF64})
  local_index = Dict{Int,Int}(bus => i for (i, bus) in enumerate(island_buses))
  n = length(island_buses)
  I = Int[]
  J = Int[]
  V = ComplexF64[]
  for br in net.branchVec
    # only fully closed branches carry a series path. A one-sided open
    # branch (r0.9.10) has no through connection, and its Schur input
    # admittance is charging, which this series-only matrix deliberately
    # drops for closed branches as well, so it contributes nothing here.
    _branch_terminal_state(br) == :closed || continue
    f = get(local_index, Int(br.fromBus), 0)
    t = get(local_index, Int(br.toBus), 0)
    (f == 0 || t == 0) && continue
    # A physical line/transformer has r >= 0. A NEGATIVE series resistance can
    # Short circuit uses the PHYSICAL branch impedance (#329): under a fault the
    # series-FACTS converter is bypassed, so the equipment resistance/reactance
    # governs the fault current, not a compensated operating point a power-flow
    # control run stamped onto r_pu/x_pu. calcBranchYserBase reads the base
    # (r_base_pu/x_base_pu). The base is physical, so this defensive check never
    # fires in the normal FACTS workflow; it guards only a corrupted base model.
    if br.r_base_pu < 0.0
      error("runShortCircuit!: branch $(getCompName(br.comp)) has a negative base series resistance (r_base = $(br.r_base_pu) pu). The physical equipment impedance must be non-negative.")
    end
    ys = calcBranchYserBase(br)
    ratio = calcBranchRatio(br)
    # Same stamping convention as calcAdmittance, minus the shunt halves.
    push!(I, f); push!(J, f); push!(V, ys / abs2(ratio))
    push!(I, f); push!(J, t); push!(V, -ys / conj(ratio))
    push!(I, t); push!(J, f); push!(V, -ys / ratio)
    push!(I, t); push!(J, t); push!(V, ys)
  end
  has_source = false
  for (bus, y) in y_add
    li = get(local_index, bus, 0)
    li == 0 && continue
    push!(I, li); push!(J, li); push!(V, y)
    has_source = true
  end
  return SparseArrays.sparse(I, J, V, n, n), local_index, has_source
end

"""
    runShortCircuit!(net, sc_data; buses = :all, case = :max, c_factor = 0.0) -> ShortCircuitResult
    runShortCircuit!(result::CGMESImporter.CGMESImportResult; kwargs...) -> ShortCircuitResult

Balanced 3-phase initial symmetrical short-circuit current per IEC 60909-0
(issue #277): equivalent voltage source at the fault location, positive
sequence only, series impedances only.

- `sc_data`: the CGMES short-circuit harvest (`CGMESImportResult.shortcircuit`).
- `buses`: `:all` (per-bus sweep as a loop over the per-bus path) or a vector
  of bus names.
- `case`: `:max` (equipment rating, `c_max`) or `:min` (protection
  sensitivity, `c_min`). Asynchronous machines never contribute to `:min`.
- `c_factor`: `0.0` uses the hardcoded IEC 60909-0 Table-1 values per voltage
  level (`short_circuit.c_factor` carries the same semantics in YAML); a
  positive scalar overrides the table for every bus.

Per fault bus the short-circuit impedance is obtained by a Z-bus column
solve: one sparse LU factorization per island and case, one triangular solve
per fault bus — the full inverse is never formed. Outputs per bus: `Ik''`
(kA), `Sk''` (MVA), and `kappa`/`i_p` from the R/X ratio at the fault
location (IEC 60909-0 method b: `κ = min(1.15·(1.02 + 0.98·e^{-3R/X}),
cap)` with cap 1.8 below 1 kV and 2.0 above).

Substituted defaults and skipped contributions are flagged on the affected
result rows (`contains_defaulted_data` + reasons), not only logged — for
`:max` a skipped contribution makes the result a lower bound. Buses in
islands without any short-circuit source report `status = :no_source` with
`NaN` currents; isolated buses report `status = :isolated`.

The function does not modify `net`; the `!` marks it as the acting entry
point of the module family (`runpf!`, `rundcpf!`). Failure behavior: throws
`ArgumentError` for an unknown bus name or an invalid `case`.

All-bus sweeps fan out over Julia threads (Phase 3 of the multi-core work):
the fault-bus list is split into `runtime.parallel.max_tasks` chunks, each
task solving on its own `copy` of the island factorization with reusable
RHS/solution buffers. Results are row-identical to the serial sweep. The
keywords `parallel_enabled`, `parallel_max_tasks`, and
`parallel_min_work_items` override the active `runtime.parallel.*`
configuration (`nothing` = use the configured values); with one Julia
thread or below `min_work_items` fault buses the sweep runs serially.

`sweep_method = :takahashi` (opt-in; default `:solves`) computes all
Thevenin impedances of an island in ONE selected-inverse pass over the
factorization instead of one triangular solve per bus
([`_takahashi_diag`](@ref); measured 34x to 264x over the serial sweep,
growing with island size). Values agree with `:solves` to machine
precision but not bitwise; islands where the method does not apply
(unsymmetric pivot ordering) automatically fall back to solves.
"""
function runShortCircuit!(
  net::Net,
  sc_data;
  buses = :all,
  case::Symbol = :max,
  c_factor::Real = 0.0,
  parallel_enabled::Union{Nothing,Bool} = nothing,
  parallel_max_tasks::Union{Nothing,Int} = nothing,
  parallel_min_work_items::Union{Nothing,Int} = nothing,
  sweep_method::Symbol = :solves,
)::ShortCircuitResult
  case in (:max, :min) || throw(ArgumentError("runShortCircuit!: case must be :max or :min; got $(case)."))
  sweep_method in (:solves, :takahashi) || throw(ArgumentError("runShortCircuit!: sweep_method must be :solves or :takahashi; got $(sweep_method)."))
  c_override = Float64(c_factor)
  (c_override == 0.0 || (isfinite(c_override) && 0.5 <= c_override <= 1.2)) || throw(ArgumentError("runShortCircuit!: c_factor must be 0 (automatic Table 1) or within [0.5, 1.2]; got $(c_factor)."))

  fault_buses = if buses === :all
    collect(eachindex(net.nodeVec))
  else
    [begin
      idx = get(net.busDict, String(b), nothing)
      idx === nothing && throw(ArgumentError("runShortCircuit!: unknown bus $(b)."))
      idx
    end for b in buses]
  end

  # Island decomposition reused from the AC island diagnostics (BFS over
  # in-service branches); isolated buses form single-bus islands there.
  components = _ac_island_components(net)
  bus_island = Dict{Int,Int}()
  for comp in components
    for bus in comp.buses
      bus_island[bus] = comp.island_id
    end
  end

  messages = String[]
  y_add, bus_reasons, island_reasons = _sc_source_admittances(net, sc_data, case, c_override, bus_island, messages)

  # Result rows carry the user-facing bus names (busDict keys), not the
  # generated component names ("Bus_<idx>_<vn>").
  bus_names = Dict{Int,String}(idx => name for (name, idx) in net.busDict)

  iso = Set(net.isoNodes)
  needed_islands = Set(get(bus_island, bus, 0) for bus in fault_buses if !(bus in iso))
  delete!(needed_islands, 0)

  # One factorization per island (and case), one triangular solve per fault
  # bus in it.
  island_solutions = Dict{Int,Any}()
  for comp in components
    comp.island_id in needed_islands || continue
    island_buses = [Int(b) for b in comp.buses if !(b in iso)]
    isempty(island_buses) && continue
    Ysc, local_index, has_source = _sc_island_matrix(net, island_buses, y_add)
    if !has_source
      island_solutions[comp.island_id] = (factor = nothing, local_index = local_index)
      continue
    end
    island_solutions[comp.island_id] = (factor = LinearAlgebra.lu(Ysc), local_index = local_index)
  end

  # sweep_method = :takahashi: ONE selected-inverse pass per island replaces
  # the per-bus triangular solves (Erisman-Tinney on the UMFPACK factors;
  # requires the symmetric pivot ordering UMFPACK picks on the structurally
  # symmetric Ysc). Islands where it does not apply fall back to solves;
  # values agree with :solves to machine precision, not bitwise.
  island_diag = Dict{Int,Vector{ComplexF64}}()
  if sweep_method === :takahashi
    for (island, sol) in island_solutions
      sol.factor === nothing && continue
      diagZ, ok, _ = _takahashi_diag(sol.factor)
      ok && (island_diag[island] = diagZ)
    end
  end

  # Phase 3: the per-bus row is a pure function of its inputs (no shared
  # container is touched), so the all-bus sweep can fan out over task
  # chunks. `solution_for(island)` supplies the factorization (the shared
  # one on the serial path, one `copy` per chunk on the parallel path;
  # a shared UMFPACK factor would serialize on its internal lock) and
  # `workspace_for(island, n)` supplies caller-owned RHS/solution buffers
  # (the former per-bus `zeros` allocation was the sweep's hot spot).
  compute_row = function (bus::Int, solution_for, workspace_for)
    node = net.nodeVec[bus]
    vn = getNodeVn(node)
    name = get(bus_names, bus, getCompName(node.comp))
    island = get(bus_island, bus, 0)
    reasons = String[]
    append!(reasons, get(bus_reasons, bus, String[]))
    append!(reasons, get(island_reasons, island, String[]))
    unique!(reasons)
    c = _sc_c_factor(vn, case, c_override)
    if bus in iso
      return (bus = name, bus_idx = bus, vn_kV = vn, island = island, status = :isolated, c = c, zk_ohm = NaN, rx_ratio = NaN, ik_kA = NaN, sk_MVA = NaN, kappa = NaN, ip_kA = NaN, contains_defaulted_data = !isempty(reasons), reasons = reasons)
    end
    sol = solution_for(island)
    if sol === nothing || sol.factor === nothing
      push!(reasons, "no short-circuit source in island $(island) — fault current undefined")
      return (bus = name, bus_idx = bus, vn_kV = vn, island = island, status = :no_source, c = c, zk_ohm = NaN, rx_ratio = NaN, ik_kA = NaN, sk_MVA = NaN, kappa = NaN, ip_kA = NaN, contains_defaulted_data = true, reasons = reasons)
    end
    li = sol.local_index[bus]
    dz = get(island_diag, island, nothing)
    z_pu = if dz === nothing
      rhs, solbuf = workspace_for(island, length(sol.local_index))
      fill!(rhs, 0.0 + 0.0im)
      rhs[li] = 1.0
      LinearAlgebra.ldiv!(solbuf, sol.factor, rhs)
      solbuf[li]
    else
      dz[li]
    end
    z_base = vn^2 / net.baseMVA
    z_ohm = z_pu * z_base
    zk = abs(z_ohm)
    ik = c * vn / (sqrt(3.0) * zk)
    sk = sqrt(3.0) * vn * ik
    # R/X at the fault location (method b); a non-positive reactance would
    # mean a badly conditioned source model — guard with a flag instead of
    # producing a silent kappa.
    rx = imag(z_ohm) > 0.0 ? max(real(z_ohm), 0.0) / imag(z_ohm) : NaN
    kappa = if isfinite(rx)
      min(1.15 * (1.02 + 0.98 * exp(-3.0 * rx)), vn <= 1.0 ? 1.8 : 2.0)
    else
      push!(reasons, "short-circuit: non-positive reactance at the fault point — kappa/i_p unavailable")
      NaN
    end
    ip = isfinite(kappa) ? kappa * sqrt(2.0) * ik : NaN
    return (bus = name, bus_idx = bus, vn_kV = vn, island = island, status = :ok, c = c, zk_ohm = zk, rx_ratio = rx, ik_kA = ik, sk_MVA = sk, kappa = kappa, ip_kA = ip, contains_defaulted_data = !isempty(reasons), reasons = reasons)
  end

  # lazily created per-owner workspaces: one (rhs, sol) buffer pair per
  # island, sized to the island's reduced system
  make_workspace_for = function ()
    ws = Dict{Int,Tuple{Vector{ComplexF64},Vector{ComplexF64}}}()
    return (island, n) -> get!(ws, island) do
      (zeros(ComplexF64, n), zeros(ComplexF64, n))
    end
  end

  parallel_on, parallel_cap, parallel_min_items = _resolve_parallel_runtime(parallel_enabled, parallel_max_tasks, parallel_min_work_items)
  use_parallel = parallel_on && Threads.nthreads() > 1 && parallel_cap > 1 && length(fault_buses) >= parallel_min_items

  if use_parallel
    # one task per chunk (number of chunks = max_tasks), one factor COPY and
    # one workspace pair per island PER CHUNK, never per bus; row slots are
    # written at distinct indices, then collected in input order. The factor
    # copies are created HERE, serially, before any task spawns: copy() on
    # one shared UmfpackLU is not safe against concurrent copies of the same
    # object (measured: silent garbage), only the SOLVES on distinct copies
    # are concurrency-safe.
    row_slots = Vector{Any}(undef, length(fault_buses))
    chunks = collect(Iterators.partition(eachindex(fault_buses), cld(length(fault_buses), parallel_cap)))
    chunk_solutions = [Dict{Int,Any}(island => (sol.factor === nothing ? sol : (factor = copy(sol.factor), local_index = sol.local_index)) for (island, sol) in island_solutions) for _ in chunks]
    tasks = [
      Threads.@spawn begin
        # `local` is load-bearing: without it these assignments would bind
        # the ENCLOSING function's variables (the serial branch defines the
        # same names), and every task would share one boxed slot — a data
        # race that silently corrupts results (measured before the fix)
        local task_solution_for = island -> get(chunk_solutions[ci], island, nothing)
        local task_workspace_for = make_workspace_for()
        for idx in chunks[ci]
          row_slots[idx] = compute_row(fault_buses[idx], task_solution_for, task_workspace_for)
        end
      end for ci in eachindex(chunks)
    ]
    foreach(wait, tasks)
    rows = NamedTuple[row_slots[i] for i in eachindex(fault_buses)]
  else
    shared_solution_for = island -> get(island_solutions, island, nothing)
    workspace_for = make_workspace_for()
    rows = NamedTuple[compute_row(bus, shared_solution_for, workspace_for) for bus in fault_buses]
  end

  return ShortCircuitResult(case, net.baseMVA, c_override, rows, messages)
end

runShortCircuit!(result::CGMESImporter.CGMESImportResult; kwargs...) = runShortCircuit!(result.net, result.shortcircuit; kwargs...)

# Native-data convenience (issue #299): run on the net's own sc_sources
# (filled by addExternalGrid!). An entirely empty container already yields
# :no_source rows through the engine's safety-flag machinery — deliberately
# no pre-check here.
runShortCircuit!(net::Net; kwargs...) = runShortCircuit!(net, net.sc_sources; kwargs...)

"""
    printShortCircuitResult(io::IO, result::ShortCircuitResult; max_rows = 50)
    printShortCircuitResult(result::ShortCircuitResult; max_rows = 50)

Pretty-print a [`ShortCircuitResult`](@ref): one table row per fault bus
(`Ik''`, `Sk''`, `kappa`, `i_p`, flag column), followed by the reason list of
every flagged row. The Web UI "Short circuit" run writes the same rows as
`short_circuit_max.csv` / `short_circuit_min.csv`.
"""
function printShortCircuitResult(io::IO, result::ShortCircuitResult; max_rows::Int = 50)
  println(io, "Balanced 3-phase short circuit (IEC 60909-0) — case: ", result.case, result.c_factor_override > 0.0 ? ", c-factor override $(result.c_factor_override)" : ", c per IEC Table 1")
  println(io, rpad("bus", 24), rpad("Un[kV]", 9), rpad("Ik''[kA]", 11), rpad("Sk''[MVA]", 12), rpad("kappa", 8), rpad("ip[kA]", 11), rpad("status", 11), "flagged")
  shown = 0
  for row in result.rows
    shown >= max_rows && (println(io, "… ", length(result.rows) - shown, " more row(s)"); break)
    fmt = x -> isfinite(x) ? string(round(x; sigdigits = 5)) : "-"
    println(io, rpad(String(row.bus), 24), rpad(fmt(row.vn_kV), 9), rpad(fmt(row.ik_kA), 11), rpad(fmt(row.sk_MVA), 12), rpad(fmt(row.kappa), 8), rpad(fmt(row.ip_kA), 11), rpad(String(row.status), 11), row.contains_defaulted_data ? "yes" : "no")
    shown += 1
  end
  flagged = [row for row in result.rows if row.contains_defaulted_data]
  if !isempty(flagged)
    println(io, "Flagged rows (defaulted data / skipped contributions — Ik''max lower bound where noted):")
    seen = Set{String}()
    for row in flagged
      for reason in row.reasons
        reason in seen && continue
        push!(seen, reason)
        println(io, "  - ", reason)
      end
    end
  end
  return nothing
end

printShortCircuitResult(result::ShortCircuitResult; max_rows::Int = 50) = printShortCircuitResult(stdout, result; max_rows = max_rows)

# Copyright 2023-2026 Udo Schmitz
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
fault buses). Deliberately a separate type: it neither extends nor embeds
into `ACPFlowReport` (precedent: `rundcpf!` with its own report type).

# Fields
- `case::Symbol`: `:max` or `:min`.
- `base_mva::Float64`: network power base used for the per-unit assembly.
- `c_factor_override::Float64`: `0.0` = IEC 60909-0 Table 1 by voltage
  level; positive = the scalar override used for every bus.
- `rows::Vector{NamedTuple}`: one row per requested fault bus with
  `bus` (name), `bus_idx`, `vn_kV`, `island`, `status` (`:ok` /
  `:no_source` / `:isolated`), `c`, `zk_ohm`, `rx_ratio`, `ik_kA`,
  `sk_MVA`, `kappa`, `ip_kA`, `contains_defaulted_data`, `lower_bound`,
  `reasons`.
  Safety-flag contract: whenever a source value was substituted with a
  documented default or a contribution was skipped, the affected rows carry
  `contains_defaulted_data = true` plus the reason list. In the `:max` case
  a skipped contribution makes `ik_kA` a **lower bound**, which is the
  non-conservative direction for equipment rating: every row of that
  island carries `lower_bound = true`. A row on a default or a substitute
  only is flagged but no lower bound.
- `messages::Vector{String}`: the run's warnings/notices (same texts that
  were logged).
- `status::Symbol`: the run status, `:succeeded`, `:warning` (some or all
  sources without their data, entered with defaults) or `:failed` (no
  source carries anything an impedance can come from, or no bus could be
  evaluated). The service reports the same status (one rule, see
  [Short circuit](@ref short_circuit_source_data)).
- `status_reason::Union{Nothing,String}`: the reason code that goes with
  the status (`short_circuit_partial_defaults`,
  `short_circuit_defaults_only`, `short_circuit_flagged_lower_bound`,
  `short_circuit_flagged_defaults`, `short_circuit_data_missing`), or
  `nothing` for a run without any flag.
"""
struct ShortCircuitResult
  case::Symbol
  base_mva::Float64
  c_factor_override::Float64
  rows::Vector{NamedTuple}
  messages::Vector{String}
  status::Symbol
  status_reason::Union{Nothing,String}
end

# prefixes of the messages that count the equivalent branches and the
# three-winding star legs with a negative r; printShortCircuitResult prints
# them as report lines
const _SC_EB_NEGATIVE_R = "short-circuit: equivalent branches with negative r: "
const _SC_3W_NEGATIVE_R = "short-circuit: three-winding legs with negative r: "

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
# R/X ratio at the fault (and with it kappa/i_p) stays meaningful even in
# machine-dominated networks): 0.05·x''d for machines >= 100 MVA, 0.07·x''d
# below that (HV machines), 0.15·x''d for LV machines.
_sc_generator_resistance(x_pu::Float64, ratedS_MVA::Float64, ratedU_kV::Float64)::Float64 = (ratedU_kV > 1.0 ? (ratedS_MVA >= 100.0 ? 0.05 : 0.07) : 0.15) * x_pu

_sc_usable(x) = x !== nothing && isfinite(Float64(x)) && Float64(x) > 0.0

# A maximum active power of 9999 MW and above is the way case files say
# "no limit stated" (IIDM and MATPOWER cases carry it on every unit); no
# machine is that large. Taken as a size it puts 0.2 pu on 9999 MVA and
# the fault current of a 14-bus case at 2667 kA (measured on ieee14).
const SC_PMAX_PLACEHOLDER_MW = 9999.0
_sc_usable_pmax(x) = _sc_usable(x) && Float64(x) < SC_PMAX_PLACEHOLDER_MW

# The same for the rated power of a machine: a CGMES delivery that states
# +-99999.9 MW as the limits of a unit writes ratedS = sqrt(2)*99999.9 MVA
# next to them. Taken as a size, 0.2 pu on 141421 MVA made a wind unit an
# ideal source with about 1180 kA at its 380 kV bus (measured on a real
# delivery). A rated power of 9999 MVA and above is no size.
const SC_RATEDS_PLACEHOLDER_MVA = 9999.0
_sc_rated_placeholder(x) = _sc_usable(x) && Float64(x) >= SC_RATEDS_PLACEHOLDER_MVA
_sc_usable_rated(x) = _sc_usable(x) && !_sc_rated_placeholder(x)

"""
    sc_machine_size(net, machine) -> (value, source)

The size a default reactance of a synchronous machine refers to, in MVA:
the rated power of the source record (`source = :rated_s`), else the
maximum active power of the unit (`:pmax`), else `nothing` (`:none`, the
caller falls back to the network base and flags it). The maximum active
power comes from the record (`pmax_MW`) where the adapter carries it, and
otherwise from the generating unit of the network at the machine's bus:
the unit of the same name, or the only unit with a maximum at that bus.
Only units of component type `Generator` or `SynchronousMachine` count; a
converter-connected unit (`PowerElectronicsConnection`) at the same bus is
not the machine. A maximum of
`SC_PMAX_PLACEHOLDER_MW` and above is a placeholder and no size, and so
is a rated power of `SC_RATEDS_PLACEHOLDER_MVA` and above: the search goes
on with the maximum active power as if the record had no rated power.
Several unnamed units at one bus give no answer rather than a guessed one.
"""
function sc_machine_size(net::Net, m)
  rated = hasproperty(m, :ratedS_MVA) ? m.ratedS_MVA : nothing
  _sc_usable_rated(rated) && return (value = Float64(rated), source = :rated_s)
  pmax = hasproperty(m, :pmax_MW) ? m.pmax_MW : nothing
  _sc_usable_pmax(pmax) && return (value = Float64(pmax), source = :pmax)
  (hasproperty(m, :bus) && m.bus !== nothing) || return (value = nothing, source = :none)
  busidx = get(net.busDict, String(m.bus), nothing)
  busidx === nothing && return (value = nothing, source = :none)
  units = [ps for ps in net.prosumpsVec if Int(ps.comp.cFrom_bus) == busidx && (ps.comp.cTyp == Generator || ps.comp.cTyp == SynchronousMachine) && _sc_usable_pmax(ps.maxP)]
  name = (hasproperty(m, :name) && m.name !== nothing) ? String(m.name) : ""
  named = [ps for ps in units if !isempty(name) && String(ps.comp.cName) == name]
  length(named) == 1 && return (value = Float64(named[1].maxP), source = :pmax)
  length(units) == 1 && return (value = Float64(units[1].maxP), source = :pmax)
  return (value = nothing, source = :none)
end

# The state of the source data, the same rule for every case format and for
# the library and the service. A source carries data when the quantity its
# impedance comes from is in the case: x''d of a synchronous machine, the
# short-circuit current of a feeder, the reactance of an equivalent
# injection, the locked-rotor ratio of a motor. A source carries NOTHING when
# the base of a default is missing as well (a machine without rated power; a
# feeder, an equivalent or a motor without its quantity cannot be defaulted
# at all). A rated power of SC_RATEDS_PLACEHOLDER_MVA and above is a
# placeholder and counts as none.
#   every source without anything -> refused (data missing)
#   no source with data           -> warning, the table rests on defaults
#   some sources without data     -> warning, the rows they feed are flagged
#   every source with data        -> success
_sc_field(record, name::Symbol) = hasproperty(record, name) ? getproperty(record, name) : nothing

# A synchronous machine on a wind or solar generating unit is a converter
# unit in fact, and the short circuit treats it as a current source like a
# PowerElectronicsConnection. The adapters
# mark its record with `converter_unit`, the class name of the unit
# (`WindGeneratingUnit` or `SolarGeneratingUnit`; the IIDM import maps
# energySource WIND and SOLAR onto the same names). Such a record is no
# voltage source and does not count in the state of the source data, as a
# PowerElectronicsConnection does not.
_sc_converter_unit(m) = (u = _sc_field(m, :converter_unit); (u isa AbstractString && !isempty(u)) ? String(u) : nothing)

# The rated power behind the current of such a unit: ratedS unless a
# placeholder, else the maximum active power of the unit (`pmax_MW`, the
# CGMES GeneratingUnit.maxOperatingP) unless a placeholder, else none. The
# maximum active power stands for the rated power at power factor 1.
function _sc_converter_unit_size(m)
  rated = _sc_field(m, :ratedS_MVA)
  _sc_usable_rated(rated) && return (value = Float64(rated), source = :rated_s)
  pmax = _sc_field(m, :pmax_MW)
  _sc_usable_pmax(pmax) && return (value = Float64(pmax), source = :pmax)
  return (value = nothing, source = :none)
end

function _sc_source_data_state(sc)
  total = 0
  with_data = 0
  without_any = 0
  for m in sc.synchronous_machines
    _sc_converter_unit(m) === nothing || continue
    total += 1
    has = _sc_usable(_sc_field(m, :satDirectSubtransX_pu))
    has && (with_data += 1)
    (!has && !_sc_usable_rated(_sc_field(m, :ratedS_MVA))) && (without_any += 1)
  end
  for f in sc.external_network_injections
    total += 1
    has = _sc_usable(_sc_field(f, :maxInitialSymShCCurrent_A)) || _sc_usable(_sc_field(f, :minInitialSymShCCurrent_A))
    has ? (with_data += 1) : (without_any += 1)
  end
  for e in sc.equivalent_injections
    total += 1
    _sc_usable(_sc_field(e, :x_ohm)) ? (with_data += 1) : (without_any += 1)
  end
  for m in sc.asynchronous_machines
    total += 1
    # a motor without its locked-rotor ratio is skipped by the engine,
    # there is no default for it
    _sc_usable(_sc_field(m, :iaIrRatio)) ? (with_data += 1) : (without_any += 1)
  end
  refused = total > 0 && without_any == total
  return (total = total, with_data = with_data, without_any = without_any, refused = refused, defaults_only = total > 0 && with_data == 0 && !refused, partial = with_data > 0 && with_data < total, warning = !refused && with_data < total)
end

# The run status and its reason code from the state of the source data and
# the rows. A flagged row is a lower bound only where the engine set its
# typed `lower_bound` (a skipped source in its island); a row on a default
# or on a substitute such as short_circuit.converter_current_factor is not,
# so it gets the neutral reason. The service takes status and reason from
# the result (_sc_outcome), so both say the same.
function _sc_run_status(state, rows::AbstractVector)
  (state.refused || !any(row -> row.status === :ok, rows)) && return (status = :failed, reason = "short_circuit_data_missing")
  state.defaults_only && return (status = :warning, reason = "short_circuit_defaults_only")
  state.partial && return (status = :warning, reason = "short_circuit_partial_defaults")
  any(row -> row.lower_bound, rows) && return (status = :succeeded, reason = "short_circuit_flagged_lower_bound")
  any(row -> row.contains_defaulted_data, rows) && return (status = :succeeded, reason = "short_circuit_flagged_defaults")
  return (status = :succeeded, reason = nothing)
end

# The reason on every row of an island fed by converter current sources
# alone.
const _SC_CONVERTER_ONLY_ISLAND = "short-circuit: converter-only island, arithmetic sum, no voltage source"

# The Takahashi/Erisman-Tinney recursion moved to the shared module
# src/takahashi.jl (the SE diagnostics reuse the identical backward pass
# for the selected inverse). This alias keeps every short-circuit call
# site and its behavior bitwise unchanged; see `takahashi_diag` for the
# full docstring including the applicability guards and the measured
# 34x to 264x speedup over the serial sweep.
const _takahashi_diag = takahashi_diag

# One source-admittance table for the whole net: bus index → added shunt
# admittance (pu on `base_mva`), the converter current sources (bus index →
# current in pu of the bus current base, `:max` only), plus per-bus and
# per-island safety-flag reasons. Islands are flagged as a whole where a
# skipped contribution affects every fault location fed from that island
# (motor/feeder/converter skips). `lower_islands` holds the islands in which
# a skipped source makes Ik''max a lower bound; the rows read their typed
# `lower_bound` from it, not from the reason text.
function _sc_source_admittances(net::Net, sc, case::Symbol, c_override::Float64, converter_factor::Float64, bus_island::Dict{Int,Int}, messages::Vector{String})
  y_add = Dict{Int,ComplexF64}()
  i_src = Dict{Int,Float64}()
  bus_reasons = Dict{Int,Vector{String}}()
  island_reasons = Dict{Int,Vector{String}}()
  lower_islands = Set{Int}()
  flag_bus = (idx, why) -> push!(get!(Vector{String}, bus_reasons, idx), why)
  flag_island = (idx, why) -> begin
    island = get(bus_island, idx, 0)
    island == 0 && return
    push!(get!(Vector{String}, island_reasons, island), why)
  end
  # a skipped source in the maximum case: Ik''max of the island is a lower bound
  flag_lower = idx -> begin
    island = get(bus_island, idx, 0)
    island == 0 || push!(lower_islands, island)
  end
  note = why -> (push!(messages, why); @warn why)

  sbase = net.baseMVA

  # Synchronous machines: x''d on machine base → network base. Missing data
  # is substituted with documented defaults AND flagged (safety contract):
  # a silent default on a safety-relevant quantity is not acceptable.
  for m in sc.synchronous_machines
    # a machine on a wind or solar unit is a current source, below with
    # the converters
    _sc_converter_unit(m) === nothing || continue
    m.bus === nothing && continue
    busidx = get(net.busDict, String(m.bus), nothing)
    busidx === nothing && continue
    vn = getNodeVn(net.nodeVec[busidx])
    name = something(m.name, m.mrid)
    xdpp = m.satDirectSubtransX_pu
    no_reactance = xdpp === nothing || !isfinite(Float64(xdpp)) || Float64(xdpp) <= 0.0
    no_rating = !_sc_usable_rated(m.ratedS_MVA)
    # a placeholder rating is named in the reason, so the row says why the
    # stated value was not used
    placeholder = _sc_rated_placeholder(m.ratedS_MVA) ? "a placeholder ratedS of $(m.ratedS_MVA) MVA ($(SC_RATEDS_PLACEHOLDER_MVA) MVA and above is no size)" : nothing
    srated = m.ratedS_MVA
    if no_reactance
      # The default is a per-unit value on the size of the machine, so the
      # size decides the current: rated power, else the maximum active
      # power of the unit, else the network base. The last step says
      # nothing about the machine and carries its own reason.
      xdpp = 0.2
      size = sc_machine_size(net, m)
      srated = something(size.value, sbase)
      why = if size.source === :rated_s
        "short-circuit: SynchronousMachine $(name) has no usable x''_d, default 0.2 pu (machine base) substituted"
      elseif size.source === :pmax && placeholder !== nothing
        "short-circuit: SynchronousMachine $(name) has no usable x''_d and $(placeholder): default 0.2 pu on its maximum active power $(srated) MW substituted"
      elseif size.source === :pmax
        "short-circuit: SynchronousMachine $(name) has no usable x''_d and no ratedS: default 0.2 pu on its maximum active power $(srated) MW substituted"
      elseif placeholder !== nothing
        "short-circuit: SynchronousMachine $(name) has no x''_d, $(placeholder) and no usable maximum active power: its contribution rests on defaults alone (0.2 pu on the network base $(sbase) MVA)"
      else
        "short-circuit: SynchronousMachine $(name) has neither x''_d nor ratedS nor a usable maximum active power: its contribution rests on defaults alone (0.2 pu on the network base $(sbase) MVA)"
      end
      note(why)
      flag_bus(busidx, why)
      flag_island(busidx, why)
    elseif no_rating
      # a stated reactance is a per-unit value on a base the record does not
      # name; the maximum active power is no substitute here, because an
      # adapter that converts from ohm (IIDM) has formed the value on the
      # network base already
      srated = sbase
      why = placeholder === nothing ? "short-circuit: SynchronousMachine $(name) has no usable ratedS, network base $(sbase) MVA substituted" : "short-circuit: SynchronousMachine $(name) has $(placeholder): network base $(sbase) MVA substituted"
      note(why)
      flag_bus(busidx, why)
      flag_island(busidx, why)
    end
    urated = m.ratedU_kV
    if urated === nothing || !isfinite(Float64(urated)) || Float64(urated) <= 0.0
      urated = vn
      why = "short-circuit: SynchronousMachine $(name) has no usable ratedU, bus nominal $(vn) kV substituted"
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
  # IEC 60909-0 §6.2 when unknown, flagged).
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
      # bound, exactly what the safety flag exists for.
      why = "short-circuit: feeder $(name) has no usable $(case === :max ? "maxInitialSymShCCurrent" : "minInitialSymShCCurrent"), contribution skipped ($(case === :max ? "Ik''max is a lower bound" : "Ik''min may be overestimated"))"
      note(why)
      flag_island(busidx, why)
      # in the minimum case the missing infeed overestimates, no lower bound
      case === :max && flag_lower(busidx)
      continue
    end
    c = _sc_c_factor(vn, case, c_override)
    z_ohm = c * vn * 1000.0 / (sqrt(3.0) * Float64(ik_A))
    rx = case === :max ? f.maxR1ToX1Ratio : f.minR1ToX1Ratio
    if rx === nothing || !isfinite(Float64(rx)) || Float64(rx) < 0.0
      rx = 0.1
      why = "short-circuit: feeder $(name) has no usable R/X ratio, default R = 0.1·X substituted"
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
      why = "short-circuit: EquivalentInjection $(name) has no usable r, 0 Ω substituted"
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
        why = "short-circuit: AsynchronousMachine $(name) has no usable $(missing_attr), contribution skipped; Ik''max in its island is a lower bound"
        note(why)
        flag_island(busidx, why)
        flag_lower(busidx)
        continue
      end
      z_ohm_total = (1.0 / Float64(ilr)) * Float64(urated)^2 / Float64(srated)
      rx = m.rxLockedRotorRatio
      if rx === nothing || !isfinite(Float64(rx)) || Float64(rx) < 0.0
        # §6.7.2 guidance split, keyed by motor class; MW per pole pair only
        # decides the MV band and falls back to the conservative 0.15.
        mw_per_pole = m.ratedMechanicalPower_MW !== nothing && m.polePairNumber !== nothing && Float64(m.polePairNumber) > 0.0 ? Float64(m.ratedMechanicalPower_MW) / Float64(m.polePairNumber) : nothing
        rx = Float64(urated) <= 1.0 ? 0.42 : (mw_per_pole !== nothing && mw_per_pole >= 1.0 ? 0.10 : 0.15)
        why = "short-circuit: AsynchronousMachine $(name) has no usable rxLockedRotorRatio, IEC 60909-0 §6.7.2 default R/X = $(rx) substituted"
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
      why = "short-circuit: AsynchronousMachine $(ps.comp.cName): short-circuit attributes not harvested, contribution skipped; Ik''max in its island is a lower bound"
      note(why)
      flag_island(busidx, why)
      flag_lower(busidx)
    end
    # Converter-connected units (CGMES 3.0 PowerElectronicsConnection: wind,
    # solar and battery parks behind a full-size converter) are current
    # sources (IEC 60909-0:2016, clause 6.9): no admittance on the diagonal,
    # an injected current instead, which runShortCircuit! superposes on the
    # equivalent voltage source. The current is maxIFault times the rated
    # current ratedS/(sqrt(3)*ratedU), here in per unit of the current base
    # of the bus, baseMVA/(sqrt(3)*Vn). A unit without maxIFault or ratedS
    # feeds nothing and flags its island; so does a converter prosumer whose
    # bus carries no record at all (a network whose source data has none).
    # maxIFault is CIM17 and in no CGMES 3.0 profile, so on a conforming
    # delivery every converter takes the flag branch, unless
    # converter_factor (short_circuit.converter_current_factor) is positive:
    # then a unit without maxIFault feeds that multiple of its rated current
    # instead, and the rows of its island name the substitute. A stated
    # maxIFault always wins over the factor: only an absent attribute
    # (nothing) takes it, a stated value that is not a positive current
    # (0, negative, NaN) feeds nothing and flags the island with its value.
    # Ik''min needs no flag: leaving a source out lowers the current, the
    # safe side for the minimum, and converters never enter it (like motors).
    pec_bus_names = nothing
    pec_name = busidx -> begin
      pec_bus_names === nothing && (pec_bus_names = _bus_name_by_idx(net))
      get(pec_bus_names, busidx, string(busidx))
    end
    no_model = busidx -> "short-circuit: PowerElectronicsConnection at bus $(pec_name(busidx)): no converter short-circuit model, contribution skipped; Ik''max in its island is a lower bound"
    recorded_buses = Set{Int}()
    for rec in sc.power_electronics_connections
      rec.bus === nothing && continue
      busidx = get(net.busDict, String(rec.bus), nothing)
      busidx === nothing && continue
      push!(recorded_buses, busidx)
      name = something(rec.name, rec.mrid)
      from_factor = rec.maxIFault_pu === nothing && converter_factor > 0.0
      if !from_factor && !_sc_usable(rec.maxIFault_pu)
        why = rec.maxIFault_pu === nothing ? no_model(busidx) : "short-circuit: PowerElectronicsConnection $(name) has maxIFault $(rec.maxIFault_pu), not a usable current: contribution skipped; Ik''max in its island is a lower bound"
        note(why)
        flag_island(busidx, why)
        flag_lower(busidx)
        continue
      end
      if !_sc_usable(rec.ratedS_MVA)
        why = "short-circuit: PowerElectronicsConnection $(name) has $(from_factor ? "no maxIFault (short_circuit.converter_current_factor set)" : "maxIFault") but no usable ratedS, so its rated current is unknown: contribution skipped; Ik''max in its island is a lower bound"
        note(why)
        flag_island(busidx, why)
        flag_lower(busidx)
        continue
      end
      ifault_pu = if from_factor
        why = "short-circuit: PowerElectronicsConnection $(name) has no maxIFault: short_circuit.converter_current_factor $(converter_factor) times the rated current substituted"
        note(why)
        flag_bus(busidx, why)
        flag_island(busidx, why)
        converter_factor
      else
        Float64(rec.maxIFault_pu)
      end
      vn = getNodeVn(net.nodeVec[busidx])
      urated = rec.ratedU_kV
      if !_sc_usable(urated)
        urated = vn
        why = "short-circuit: PowerElectronicsConnection $(name) has no usable ratedU: bus nominal $(vn) kV substituted"
        note(why)
        flag_bus(busidx, why)
        flag_island(busidx, why)
      end
      i_pu = ifault_pu * (Float64(rec.ratedS_MVA) / sbase) * (vn / Float64(urated))
      i_src[busidx] = get(i_src, busidx, 0.0) + i_pu
    end
    # Synchronous machines on a wind or solar unit (_sc_converter_unit): the
    # same current source with another origin. They state no maxIFault, so
    # the current is converter_factor times the rated current, the rated
    # power from ratedS, else from the maximum active power of the unit
    # (_sc_converter_unit_size). With the factor off or without a size the
    # unit feeds nothing and its island is a lower bound, as a
    # PowerElectronicsConnection without maxIFault or ratedS.
    for m in sc.synchronous_machines
      origin = _sc_converter_unit(m)
      origin === nothing && continue
      m.bus === nothing && continue
      busidx = get(net.busDict, String(m.bus), nothing)
      busidx === nothing && continue
      label = "SynchronousMachine $(something(m.name, m.mrid)) on a $(origin)"
      if converter_factor <= 0.0
        why = "short-circuit: $(label) is a converter unit without a converter short-circuit model (short_circuit.converter_current_factor off): contribution skipped; Ik''max in its island is a lower bound"
        note(why)
        flag_island(busidx, why)
        flag_lower(busidx)
        continue
      end
      size = _sc_converter_unit_size(m)
      placeholder = _sc_rated_placeholder(_sc_field(m, :ratedS_MVA)) ? " (ratedS $(m.ratedS_MVA) MVA is a placeholder)" : ""
      if size.value === nothing
        why = "short-circuit: $(label) is a converter unit without a usable ratedS or maximum active power$(placeholder), so its rated current is unknown: contribution skipped; Ik''max in its island is a lower bound"
        note(why)
        flag_island(busidx, why)
        flag_lower(busidx)
        continue
      end
      basis = size.source === :rated_s ? "ratedS $(size.value) MVA" : "maximum active power $(size.value) MW$(placeholder)"
      why = "short-circuit: $(label) is a converter unit: short_circuit.converter_current_factor $(converter_factor) times the rated current from its $(basis) substituted"
      note(why)
      flag_bus(busidx, why)
      flag_island(busidx, why)
      vn = getNodeVn(net.nodeVec[busidx])
      urated = _sc_field(m, :ratedU_kV)
      if !_sc_usable(urated)
        urated = vn
        why = "short-circuit: $(label) has no usable ratedU: bus nominal $(vn) kV substituted"
        note(why)
        flag_bus(busidx, why)
        flag_island(busidx, why)
      end
      i_src[busidx] = get(i_src, busidx, 0.0) + converter_factor * (size.value / sbase) * (vn / Float64(urated))
    end
    for ps in net.prosumpsVec
      ps.comp.cTyp == PowerElectronicsConnection || continue
      busidx = Int(ps.comp.cFrom_bus)
      (busidx in eachindex(net.nodeVec) && !(busidx in recorded_buses)) || continue
      why = no_model(busidx)
      note(why)
      flag_island(busidx, why)
      flag_lower(busidx)
    end
  end

  return y_add, i_src, bus_reasons, island_reasons, lower_islands
end

# Assemble the positive-sequence short-circuit admittance matrix of one
# island: series branch impedances only (IEC 60909-0: loads, line charging
# and shunt compensators are dropped), source admittances on the diagonal.
# Reuses the equivalent-circuit helpers (calcBranchYser/calcBranchRatio)
# instead of re-implementing the branch model; deliberately NOT a patched
# createYBUS: the PF assembly stays untouched by design and the
# series-only variant needs no sparsity bookkeeping beyond this COO loop.
# `island_buses` are group representatives (`reps`, the closed-link map of
# the power flow): every bus enters the matrix through the representative
# of its group, so a closed switch is one node, as in _merged_pf_net.
function _sc_island_matrix(net::Net, island_buses::Vector{Int}, y_add::Dict{Int,ComplexF64}, reps::Vector{Int}, neg_r_allowed::Set{Int})
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
    rf = reps[Int(br.fromBus)]
    rt = reps[Int(br.toBus)]
    # a branch between two buses of one group is shorted by the closed link
    # next to it; _merged_pf_net takes it out of service as well
    (rf == rt && Int(br.fromBus) != Int(br.toBus)) && continue
    f = get(local_index, rf, 0)
    t = get(local_index, rt, 0)
    (f == 0 || t == 0) && continue
    # Short circuit uses the PHYSICAL branch impedance (#329): under a fault the
    # series-FACTS converter is bypassed, so the equipment resistance/reactance
    # governs the fault current, not a compensated operating point a power-flow
    # control run stamped onto r_pu/x_pu. calcBranchYserBase reads the base
    # (r_base_pu/x_base_pu). A physical line or transformer has r >= 0, and the
    # base is physical, so this defensive check never fires in the normal FACTS
    # workflow; it guards only a corrupted base model. An equivalent branch
    # (network reduction) and a three-winding star leg are no single equipment
    # impedance: their negative r is taken over unchanged
    # (_negative_r_allowed_branches).
    if br.r_base_pu < 0.0 && !(Int(br.branchIdx) in neg_r_allowed)
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
    # a source on any member of a group feeds the group node
    li = get(local_index, reps[bus], 0)
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
  sensitivity, `c_min`). Asynchronous machines and converter-connected
  units never contribute to `:min`.
- `c_factor`: `0.0` uses the hardcoded IEC 60909-0 Table-1 values per voltage
  level (`short_circuit.c_factor` carries the same semantics in YAML); a
  positive scalar overrides the table for every bus.
- `converter_current_factor`: `0.0` (off) or the fault current of a
  converter without `maxIFault` as a multiple of its rated current (below).

Per fault bus the short-circuit impedance is obtained by a Z-bus column
solve: one sparse LU factorization per island and case, one triangular solve
per fault bus: the full inverse is never formed. Outputs per bus: `Ik''`
(kA), `Sk''` (MVA), and `kappa`/`i_p` from the R/X ratio at the fault
location (IEC 60909-0 method b: `κ = min(1.15·(1.02 + 0.98·e^{-3R/X}),
cap)` with cap 1.8 below 1 kV and 2.0 above).

Converter-connected units (`sc_data.power_electronics_connections`, one
record per mapped CGMES 3.0 `PowerElectronicsConnection`) feed `:max` as
current sources, the model IEC 60909-0:2016 gives power station units with
full size converter (clause 6.9). The current of a unit is
`I_sk = maxIFault·ratedS/(√3·ratedU)` (`maxIFault` is in per unit of the
rated current; a missing `ratedU` is replaced by the bus nominal voltage
and flagged). `maxIFault` is a CIM17 attribute that no CGMES 3.0 profile
defines: a conforming delivery has none, and its converters take the
lower-bound flag below. `converter_current_factor` (config
`short_circuit.converter_current_factor`, default `0.0` = off) replaces a
missing `maxIFault` by that multiple of the rated current; the rows of the
island are flagged with the substitute instead of the lower bound. A
stated `maxIFault` wins, also one that is not positive. The converter
share at fault bus `i` is superposed on the equivalent voltage source:
`Ik'' = c·Un/(√3·|Z_ii|) + (1/|Z_ii|)·Σ_j |Z_ij|·I_sk,j`, with `Z` the
impedance matrix of the island in which the voltage sources are
short-circuited (the matrix the first term uses) and `j` over the
converter buses. Magnitudes are added, which bounds the complex sum from
above for any phase of the converter currents. `κ` belongs to the first
term: `i_p = κ·√2·I''_vs + √2·I''_conv`. `zk_ohm` and `rx_ratio` stay those
of `Z_ii`, so with converters `ik_kA` is more than `c·Un/(√3·zk_ohm)`. A
unit without `maxIFault` (factor off), with a `maxIFault` that is not
positive or without `ratedS`, and a converter prosumer whose bus
has no record, feeds nothing and flags its island as a lower bound. An
island fed by converters alone has no voltage source and no `Z_ii`: every
bus of it reports `:ok` with the arithmetic sum of the converter currents
of the island, `i_p = √2·Ik''`, `NaN` for `zk_ohm`, `rx_ratio` and `κ`,
and the flagged reason `converter-only island, arithmetic sum, no voltage
source` (an upper bound, no lower bound of its own). Without any feeding
converter (factor off, no `maxIFault`) such an island stays `:no_source`.

A synchronous machine record with `converter_unit` set (a CGMES
`SynchronousMachine` on a `WindGeneratingUnit` or `SolarGeneratingUnit`,
an IIDM generator with `energySource` `WIND` or `SOLAR`) takes the same
path instead of the machine path: `converter_current_factor` times the
rated current, the rated power from `ratedS_MVA`, else from `pmax_MW` (the
maximum active power of the unit, at power factor 1), each unless a
placeholder (`SC_RATEDS_PLACEHOLDER_MVA`, `SC_PMAX_PLACEHOLDER_MW`). With
the factor off or without a size it feeds nothing and flags its island as
a lower bound. Such a machine never enters `:min`.

Substituted defaults and skipped contributions are flagged on the affected
result rows (`contains_defaulted_data` + reasons), not only logged: for
`:max` a skipped contribution makes the result a lower bound. Buses in
islands without any short-circuit source report `status = :no_source` with
`NaN` currents; isolated buses report `status = :isolated`.

Closed links (busbar couplers, retained switches) are contracted as in the
power flow: the buses of one closed-link group are one node of the
short-circuit matrix. Every member of a fed group reports the `Ik''` of
that node and the same island, sources on any member feed it, and the
reasons of every member appear on each member's row. A branch between two
buses of one group is shorted by the link and left out. A group counts as
isolated only when all its members are in `net.isoNodes`.

The function does not modify `net`; the `!` marks it as the acting entry
point of the module family (`runpf!`, `rundcpf!`). Failure behavior: throws
`ArgumentError` for an unknown bus name or an invalid `case`.

All-bus sweeps fan out over Julia threads:
the fault-bus list is split into `runtime.parallel.max_tasks` chunks, each
task solving on its own `copy` of the island factorization with reusable
RHS/solution buffers. Results are row-identical to the serial sweep. The
keywords `parallel_enabled`, `parallel_max_tasks`, and
`parallel_min_work_items` override the active `runtime.parallel.*`
configuration (`nothing` = use the configured values); with one Julia
thread or below `min_work_items` fault buses the sweep runs serially.

`sweep_method` selects the all-bus Thevenin sweep strategy. `:takahashi`
computes all Thevenin impedances of an island in ONE selected-inverse
pass over the factorization instead of one triangular solve per bus
([`takahashi_diag`](@ref); measured 34x to 264x over the serial sweep,
growing with island size). Values agree with `:solves` to machine
precision but not bitwise; islands where the method does not apply
(unsymmetric pivot ordering) automatically fall back to solves. The
default `:auto` uses Takahashi for islands with at least
`takahashi_min_buses` buses and plain solves below (the pass overhead
only pays off from a certain island size); `:solves` forces the per-bus
solves everywhere. The service and Web UI paths take both settings from
`short_circuit.sweep_method` / `short_circuit.takahashi_min_buses`, and
`converter_current_factor` from `short_circuit.converter_current_factor`.
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
  sweep_method::Symbol = :auto,
  takahashi_min_buses::Int = 50,
  converter_current_factor::Real = 0.0,
)::ShortCircuitResult
  case in (:max, :min) || throw(ArgumentError("runShortCircuit!: case must be :max or :min; got $(case)."))
  sweep_method in (:auto, :solves, :takahashi) || throw(ArgumentError("runShortCircuit!: sweep_method must be :auto, :solves, or :takahashi; got $(sweep_method)."))
  takahashi_min_buses >= 1 || throw(ArgumentError("runShortCircuit!: takahashi_min_buses must be at least 1; got $(takahashi_min_buses)."))
  c_override = Float64(c_factor)
  (c_override == 0.0 || (isfinite(c_override) && 0.5 <= c_override <= 1.2)) || throw(ArgumentError("runShortCircuit!: c_factor must be 0 (automatic Table 1) or within [0.5, 1.2]; got $(c_factor)."))
  conv_factor = Float64(converter_current_factor)
  (conv_factor == 0.0 || (isfinite(conv_factor) && conv_factor > 0.0)) || throw(ArgumentError("runShortCircuit!: converter_current_factor must be 0 (off) or positive; got $(converter_current_factor)."))

  fault_buses = if buses === :all
    collect(eachindex(net.nodeVec))
  else
    [begin
      idx = get(net.busDict, String(b), nothing)
      idx === nothing && throw(ArgumentError("runShortCircuit!: unknown bus $(b)."))
      idx
    end for b in buses]
  end

  # Closed-link groups are contracted as in the power flow (_merged_pf_net):
  # `reps` is the solver's group map, every bus counts through the
  # representative of its group. Without a closed link reps[b] == b and
  # everything below is the per-bus engine as before.
  reps = _active_link_representative_map(net)
  # members of the groups with more than one bus, keyed by representative;
  # a bus missing here is a group of its own
  group_members = Dict{Int,Vector{Int}}()
  for bus in eachindex(reps)
    r = reps[bus]
    r == bus && continue
    push!(get!(() -> [r], group_members, r), bus)
  end
  members_of(bus::Int) = get(group_members, reps[bus], nothing)

  # Island decomposition reused from the AC island diagnostics (BFS over
  # in-service branches plus the closed-link groups); isolated buses form
  # single-bus islands there, an isolated group one island of its members.
  components = _ac_island_components(net; links = true)
  bus_island = Dict{Int,Int}()
  for comp in components
    for bus in comp.buses
      bus_island[bus] = comp.island_id
    end
  end

  messages = String[]
  y_add, i_src, bus_reasons, island_reasons, lower_islands = _sc_source_admittances(net, sc_data, case, c_override, conv_factor, bus_island, messages)

  # Result rows carry the user-facing bus names (busDict keys), not the
  # generated component names ("Bus_<idx>_<vn>").
  bus_names = Dict{Int,String}(idx => name for (name, idx) in net.busDict)

  # A group is isolated when all its members are (net.isoNodes, the
  # branch-only marking): a bus on a closed switch to a fed bus is part of
  # the fed node, a group without any branch stays isolated.
  iso = Set(net.isoNodes)
  for grp in values(group_members)
    all(in(iso), grp) && continue
    foreach(b -> delete!(iso, b), grp)
  end
  needed_islands = Set(get(bus_island, bus, 0) for bus in fault_buses if !(bus in iso))
  delete!(needed_islands, 0)

  # One factorization per island (and case), one triangular solve per fault
  # bus in it. The matrix has one node per group (its representative).
  island_solutions = Dict{Int,Any}()
  neg_r_allowed = _negative_r_allowed_branches(net)
  for comp in components
    comp.island_id in needed_islands || continue
    island_buses = unique!([reps[Int(b)] for b in comp.buses if !(b in iso)])
    isempty(island_buses) && continue
    Ysc, local_index, has_source = _sc_island_matrix(net, island_buses, y_add, reps, neg_r_allowed)
    if !has_source
      island_solutions[comp.island_id] = (factor = nothing, local_index = local_index)
      continue
    end
    island_solutions[comp.island_id] = (factor = LinearAlgebra.lu(Ysc), local_index = local_index)
  end

  # equivalent branches and three-winding star legs with a negative series
  # resistance in the islands of the fault buses: taken over unchanged,
  # counted in one message per kind. A branch inside one closed-link group is
  # shorted and never enters the matrix (_sc_island_matrix), so it is not
  # counted either.
  neg_r_used = filter(br -> Int(br.branchIdx) in neg_r_allowed && br.r_base_pu < 0.0 && _branch_terminal_state(br) == :closed && reps[Int(br.fromBus)] != reps[Int(br.toBus)] && get(bus_island, Int(br.fromBus), 0) in needed_islands, net.branchVec)
  n_eb_neg = count(br -> haskey(net.equivalent_branch_metadata, Int(br.branchIdx)), neg_r_used)
  n_3w_neg = length(neg_r_used) - n_eb_neg
  n_eb_neg == 0 || push!(messages, string(_SC_EB_NEGATIVE_R, n_eb_neg, " in the islands of the fault buses, taken over unchanged"))
  n_3w_neg == 0 || push!(messages, string(_SC_3W_NEGATIVE_R, n_3w_neg, " in the islands of the fault buses, taken over unchanged"))

  # Converter current sources: per island the vector s_i = sum_j |Z_ij|·I_j
  # over the converter buses j, one column solve per converter bus on the
  # island factorization (Y·z = e_j gives column j of Z). The sweep divides
  # by |Z_ii|. Formed here, serially and before any task spawns, so the
  # parallel sweep only reads it. An island with converters but no voltage
  # source has no factorization: it gets the arithmetic sum of its converter
  # currents in pu (conv_only_sum, see the row below).
  island_conv = Dict{Int,Vector{Float64}}()
  conv_only_sum = Dict{Int,Float64}()
  for (bus, i_pu) in i_src
    island = get(bus_island, bus, 0)
    sol = get(island_solutions, island, nothing)
    sol === nothing && continue
    if sol.factor === nothing
      conv_only_sum[island] = get(conv_only_sum, island, 0.0) + i_pu
      continue
    end
    lj = get(sol.local_index, reps[bus], 0)
    lj == 0 && continue
    n = length(sol.local_index)
    rhs = zeros(ComplexF64, n)
    rhs[lj] = 1.0
    col = sol.factor \ rhs
    s = get!(() -> zeros(Float64, n), island_conv, island)
    s .+= abs.(col) .* i_pu
  end

  # sweep_method = :takahashi: ONE selected-inverse pass per island replaces
  # the per-bus triangular solves (Erisman-Tinney on the UMFPACK factors;
  # requires the symmetric pivot ordering UMFPACK picks on the structurally
  # symmetric Ysc). Islands where it does not apply fall back to solves;
  # values agree with :solves to machine precision, not bitwise.
  island_diag = Dict{Int,Vector{ComplexF64}}()
  if sweep_method !== :solves
    # :auto gates the pass per island: below takahashi_min_buses the
    # per-bus solves are cheaper than the selected-inverse setup
    island_sizes = Dict{Int,Int}(comp.island_id => length(comp.buses) for comp in components)
    for (island, sol) in island_solutions
      sol.factor === nothing && continue
      sweep_method === :auto && get(island_sizes, island, 0) < takahashi_min_buses && continue
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
    # a fault on a bus is a fault on its whole group: the reasons of every
    # member apply. `local`: the row closure runs in parallel tasks and must
    # not share a variable with the enclosing function
    local row_group = members_of(bus)
    if row_group === nothing
      append!(reasons, get(bus_reasons, bus, String[]))
    else
      for member in row_group
        append!(reasons, get(bus_reasons, member, String[]))
      end
    end
    append!(reasons, get(island_reasons, island, String[]))
    unique!(reasons)
    # `local`: see row_group above
    local lower = island in lower_islands
    c = _sc_c_factor(vn, case, c_override)
    if bus in iso
      return (bus = name, bus_idx = bus, vn_kV = vn, island = island, status = :isolated, c = c, zk_ohm = NaN, rx_ratio = NaN, ik_kA = NaN, sk_MVA = NaN, kappa = NaN, ip_kA = NaN, contains_defaulted_data = !isempty(reasons), lower_bound = lower, reasons = reasons)
    end
    sol = solution_for(island)
    i_conv_only = sol === nothing ? 0.0 : get(conv_only_sum, island, 0.0)
    if i_conv_only > 0.0
      # An island fed by converter current sources alone: no voltage
      # source, so no Thevenin impedance,
      # but the current sources still feed. IEC 60909 neglects loads and
      # shunt branches at a bolted fault, so their whole current flows into
      # it. The arithmetic sum is the upper bound (all currents in phase):
      # the row is flagged but no lower bound of its own, and the peak has
      # no kappa (no DC component), as the converter share above.
      ik = i_conv_only * net.baseMVA / (sqrt(3.0) * vn)
      push!(reasons, _SC_CONVERTER_ONLY_ISLAND)
      return (bus = name, bus_idx = bus, vn_kV = vn, island = island, status = :ok, c = c, zk_ohm = NaN, rx_ratio = NaN, ik_kA = ik, sk_MVA = sqrt(3.0) * vn * ik, kappa = NaN, ip_kA = sqrt(2.0) * ik, contains_defaulted_data = true, lower_bound = lower, reasons = reasons)
    end
    if sol === nothing || sol.factor === nothing
      push!(reasons, "no short-circuit source in island $(island), fault current undefined")
      return (bus = name, bus_idx = bus, vn_kV = vn, island = island, status = :no_source, c = c, zk_ohm = NaN, rx_ratio = NaN, ik_kA = NaN, sk_MVA = NaN, kappa = NaN, ip_kA = NaN, contains_defaulted_data = true, lower_bound = lower, reasons = reasons)
    end
    li = sol.local_index[reps[bus]]
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
    ik_vs = c * vn / (sqrt(3.0) * zk)
    # converter share (1/|Z_ii|)·sum_j |Z_ij|·I_j, from pu to kA with the
    # current base baseMVA/(sqrt(3)·Vn) of the fault bus
    s_conv = get(island_conv, island, nothing)
    ik_pf = s_conv === nothing ? 0.0 : s_conv[li] / abs(z_pu) * net.baseMVA / (sqrt(3.0) * vn)
    ik = ik_vs + ik_pf
    sk = sqrt(3.0) * vn * ik
    # R/X at the fault location (method b); a non-positive reactance would
    # mean a badly conditioned source model. Guard with a flag instead of
    # producing a silent kappa.
    rx = imag(z_ohm) > 0.0 ? max(real(z_ohm), 0.0) / imag(z_ohm) : NaN
    kappa = if isfinite(rx)
      min(1.15 * (1.02 + 0.98 * exp(-3.0 * rx)), vn <= 1.0 ? 1.8 : 2.0)
    else
      push!(reasons, "short-circuit: non-positive reactance at the fault point, kappa/i_p unavailable")
      NaN
    end
    # kappa belongs to the voltage-source part; the converter part enters
    # the peak with sqrt(2) alone
    ip = isfinite(kappa) ? kappa * sqrt(2.0) * ik_vs + sqrt(2.0) * ik_pf : NaN
    return (bus = name, bus_idx = bus, vn_kV = vn, island = island, status = :ok, c = c, zk_ohm = zk, rx_ratio = rx, ik_kA = ik, sk_MVA = sk, kappa = kappa, ip_kA = ip, contains_defaulted_data = !isempty(reasons), lower_bound = lower, reasons = reasons)
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
        # same names), and every task would share one boxed slot: a data
        # race that silently corrupts results (measured before the fix)
        local task_solution_for = island -> get(chunk_solutions[ci], island, nothing)
        local task_workspace_for = make_workspace_for()
        for idx in chunks[ci]
          # an all-bus sweep is the longest short-circuit run there is, so it
          # has to be interruptible; the token is process-wide, so it is
          # visible inside this spawned task
          sparlectra_check_abort()
          row_slots[idx] = compute_row(fault_buses[idx], task_solution_for, task_workspace_for)
        end
      end for ci in eachindex(chunks)
    ]
    foreach(wait, tasks)
    rows = NamedTuple[row_slots[i] for i in eachindex(fault_buses)]
  else
    shared_solution_for = island -> get(island_solutions, island, nothing)
    workspace_for = make_workspace_for()
    rows = NamedTuple[(sparlectra_check_abort(); compute_row(bus, shared_solution_for, workspace_for)) for bus in fault_buses]
  end

  run_status = _sc_run_status(_sc_source_data_state(sc_data), rows)
  return ShortCircuitResult(case, net.baseMVA, c_override, rows, messages, run_status.status, run_status.reason)
end

# the CGMES-result convenience overload lives in the adapters block
# (src/adapters/cgmes_glue.jl): it dispatches on an adapter type, and the
# core block must not reference adapter types at definition time

# Native-data convenience (issue #299): run on the net's own sc_sources
# (filled by addExternalGrid!). An entirely empty container already yields
# :no_source rows through the engine's safety-flag machinery, deliberately
# no pre-check here.
runShortCircuit!(net::Net; kwargs...) = runShortCircuit!(net, net.sc_sources; kwargs...)

# The head line of a result table: the run status with its reason, the
# number of flagged rows and the number of rows whose Ik''max is a lower
# bound. printShortCircuitResult starts with it, and the service puts the
# same line into the run metadata, from where the Web UI shows it.
function _sc_status_line(result::ShortCircuitResult)::String
  flagged = count(row -> row.contains_defaulted_data, result.rows)
  lower = count(row -> row.lower_bound, result.rows)
  status = result.status_reason === nothing ? String(result.status) : string(result.status, " (", result.status_reason, ")")
  return string("run status: ", status, ", ", flagged, " of ", length(result.rows), " rows flagged, ", lower, " rows with a lower bound")
end

"""
    printShortCircuitResult(io::IO, result::ShortCircuitResult; max_rows = 50)
    printShortCircuitResult(result::ShortCircuitResult; max_rows = 50)

Pretty-print a [`ShortCircuitResult`](@ref). The first line names the run
status with its reason, the number of flagged rows (defaulted or substituted
data, skipped sources) and the number of rows whose `Ik''max` is a lower
bound. Then one table row per fault bus (`Ik''`, `Sk''`, `kappa`, `i_p`,
flag column), followed by the reason list of every flagged row. The Web UI
"Short circuit" run writes the same rows as `short_circuit_max.csv` /
`short_circuit_min.csv` and shows the first line above the table.
"""
function printShortCircuitResult(io::IO, result::ShortCircuitResult; max_rows::Int = 50)
  println(io, _sc_status_line(result))
  println(io, "Balanced 3-phase short circuit (IEC 60909-0), case: ", result.case, result.c_factor_override > 0.0 ? ", c-factor override $(result.c_factor_override)" : ", c per IEC Table 1")
  foreach(m -> println(io, m), filter(m -> startswith(m, _SC_EB_NEGATIVE_R) || startswith(m, _SC_3W_NEGATIVE_R), result.messages))
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
    println(io, "Flagged rows (defaulted data / skipped contributions, Ik''max lower bound where noted):")
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

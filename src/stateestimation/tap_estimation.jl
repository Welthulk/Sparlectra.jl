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

# file: src/stateestimation/tap_estimation.jl
# purpose: transformer tap estimation (0.10.0): release API
#          (setTapEstimation!), the cascade tap model with fixed regulator
#          directions, and the r0 split from the current branch state. The
#          SE overlay that consumes these lives in state_estimation.jl.

"""
    _resolve_trafo_branch(net, trafo) -> Int

Resolve a transformer reference to its branch index: an `Int` branch index,
or a `String` matching the branch component name, component id, or a
numeric index. Errors when nothing matches or the branch is not a
transformer (`ratio == 0` marks a line).
"""
function _resolve_trafo_branch(net::Net, trafo)::Int
  k = 0
  if trafo isa Integer
    k = Int(trafo)
  else
    s = String(trafo)
    ki = tryparse(Int, s)
    if ki !== nothing
      k = ki
    else
      hits = [i for (i, br) in enumerate(net.branchVec) if getCompName(br.comp) == s || br.comp.cID == s]
      length(hits) == 1 || error("setTapEstimation!: transformer '$(s)' " * (isempty(hits) ? "not found" : "is ambiguous ($(length(hits)) matches)"))
      k = hits[1]
    end
  end
  1 <= k <= length(net.branchVec) || error("setTapEstimation!: branch index $(k) out of range")
  net.branchVec[k].ratio != 0.0 || error("setTapEstimation!: branch $(k) is not a transformer")
  return k
end

"""
    _cascade_tap(br, r1, r2, alpha_deg) -> (tap_ratio, phase_shift_deg)

Cascade tap model with FIXED regulator directions:

    t(r1, r2) = ratio_base / ((1 + r1) * (1 + r2 * cis(deg2rad(alpha))))

built strictly from the established convention: one `calcSkewAngleTap`
(`:reciprocal_from_side`) factor per regulator (r1 longitudinal, direction
0; r2 with the nameplate direction alpha), composed onto the neutral base
`(ratio, angle)`. No re-derived admittance terms anywhere: the caller
evaluates `calcAdmittance` on a scratch branch carrying these two fields.
"""
function _cascade_tap(br::Branch, r1::Float64, r2::Float64, alpha_deg::Float64)
  s1 = calcSkewAngleTap(tap_fraction = r1, skew_angle_deg = 0.0)
  s2 = calcSkewAngleTap(tap_fraction = r2, skew_angle_deg = alpha_deg)
  return (tap_ratio = br.ratio * s1.effective_ratio * s2.effective_ratio, phase_shift_deg = br.angle + s1.effective_shift_deg + s2.effective_shift_deg)
end

"""
    _tap_r0_split(br, mode, alpha_deg; tol = 1e-12) -> (r1_0, r2_0)

Initialize the regulator states from the CURRENT branch position: with
`t_base = calcComplexRatio(ratio, angle)` and `t_cur = calcBranchRatio(br)`
the cascade demands `v = t_base / t_cur = (1 + r1)(1 + r2 cis(alpha))`.
`r2_0` comes from the component along `cis(alpha)`, `r1_0` from the
remaining real factor. For single-regulator modes the frozen component
must vanish to `tol`; a violation means the device metadata contradicts
the branch state and errors instead of being silently absorbed.
"""
function _tap_r0_split(br::Branch, mode::Symbol, alpha_deg::Float64; tol::Float64 = 1e-12)
  t_base = calcComplexRatio(tapRatio = br.ratio, angleInDegrees = br.angle)
  t_cur = calcBranchRatio(br)
  v = t_base / t_cur
  a = real(v)
  b = imag(v)
  if mode == :ratio
    abs(b) <= tol * max(abs(v), 1.0) || error("setTapEstimation!: the branch state carries a phase component (imag = $(b)) but mode :ratio declares a longitudinal regulator only; the device metadata contradicts the branch state")
    return (a - 1.0, 0.0)
  elseif mode == :pst
    r2c = (v - 1.0) * cis(-deg2rad(alpha_deg))
    abs(imag(r2c)) <= tol * max(abs(v), 1.0) || error("setTapEstimation!: the branch state is inconsistent with a pure regulator of direction alpha = $(alpha_deg) degrees; the device metadata contradicts the branch state")
    return (0.0, real(r2c))
  elseif mode == :both
    sa = sin(deg2rad(alpha_deg))
    abs(sa) > 1e-9 || error("setTapEstimation!: alpha = $(alpha_deg) makes the two regulators collinear; use mode = :ratio")
    u = b / sa                     # u = (1 + r1) * r2
    onePlusR1 = a - u * cos(deg2rad(alpha_deg))
    onePlusR1 > 0.0 || error("setTapEstimation!: cascade split failed (non-physical longitudinal factor $(onePlusR1) from the branch state)")
    return (onePlusR1 - 1.0, u / onePlusR1)
  end
  error("setTapEstimation!: mode must be :ratio, :pst, or :both")
end

"""
    TapStateMap

State-map of the released transformer taps inside one SE (island) run:
per released trafo the snet branch position, mode, fixed direction alpha,
the state column of each regulator (0 = frozen at its r0), the r0 values,
and a dedicated scratch `Branch` whose `(tap_ratio, phase_shift_deg)` are
set from the cascade before every prediction (directive: the four
admittance terms always come from `calcAdmittance` on that scratch branch,
never re-derived). `col0` is the first tap state column, `ncols` the
number of tap states.
"""
struct TapStateMap
  branch_idxs::Vector{Int}
  modes::Vector{Symbol}
  alphas::Vector{Float64}
  r1cols::Vector{Int}
  r2cols::Vector{Int}
  r1_0::Vector{Float64}
  r2_0::Vector{Float64}
  scratch::Vector{Branch}
  col0::Int
  ncols::Int
end

## build the map over the (contracted, island) net; nothing when no closed
## transformer is released. Column layout: r1 before r2 per trafo, trafos in
## branch order, starting at col0 (after the shunt B states, before alpha).
function _build_tap_state_map(net::Net, col0::Int)::Union{Nothing,TapStateMap}
  idxs = [k for (k, br) in enumerate(net.branchVec) if br isa Branch && br.tap_est_mode != :none && _branch_terminal_state(br) == :closed]
  isempty(idxs) && return nothing
  modes = Symbol[]
  alphas = Float64[]
  r1c = Int[]
  r2c = Int[]
  r10 = Float64[]
  r20 = Float64[]
  scratch = Branch[]
  c = col0
  for k in idxs
    br = net.branchVec[k]
    a, b = _tap_r0_split(br, br.tap_est_mode, br.tap_est_alpha_deg)
    push!(modes, br.tap_est_mode)
    push!(alphas, br.tap_est_alpha_deg)
    push!(r10, a)
    push!(r20, b)
    if br.tap_est_mode in (:ratio, :both)
      push!(r1c, c)
      c += 1
    else
      push!(r1c, 0)
    end
    if br.tap_est_mode in (:pst, :both)
      push!(r2c, c)
      c += 1
    else
      push!(r2c, 0)
    end
    push!(scratch, deepcopy(br))
  end
  return TapStateMap(idxs, modes, alphas, r1c, r2c, r10, r20, scratch, col0, c - col0)
end

"""
Per-iteration STEP LIMIT of each released regulator state, in the unit of
the state itself (r1 as a tap FRACTION, r2 as the additional-voltage
amplitude or the regulating-vector ratio, matching `_tap_r0_split`).

The limit is derived from the changer's declared mechanical band: one
iteration may move a regulator by at most a quarter of its full travel.
For regulator 1 the multiplier band `[tap_min, tap_max]` maps through the
reciprocal `1/(1+r1)` onto `[1/tap_max - 1, 1/tap_min - 1]`.

Why a step limit and NOT a hard band on the state: without any bound the
Gauss-Newton step drives `r1` toward -1, where the cascade
`t = t_base/((1+r1)(1+r2 e^{jalpha}))` is singular and the stamped
admittance blows up (measured 2026-09-06: one weakly determined 50 kV
transformer ran to -329 electrical steps on a band of about -14 to +18 and
took the whole estimation down). But clamping the STATE into the band was
measured to be worse than no bound at all: five transformers that had
converged stuck to the lower bound and diverged. The transient excursion is
part of how a weakly determined regulator finds its position, so what must
be bounded is the distance ONE iteration may travel, not the region the
path may visit.
"""
function _tap_step_limits(map::TapStateMap, net::Net)::Dict{Int,Float64}
  # a quarter of the full mechanical travel per iteration
  const_share = 0.25
  limits = Dict{Int,Float64}()
  for j in eachindex(map.branch_idxs)
    br = net.branchVec[map.branch_idxs[j]]
    if map.r1cols[j] != 0 && br.tap_min > 0.0 && br.tap_max > 0.0
      lo = 1.0 / br.tap_max - 1.0
      hi = 1.0 / br.tap_min - 1.0
      limits[map.r1cols[j]] = const_share * abs(hi - lo)
    end
    if map.r2cols[j] != 0
      if br.phase_du_step > 0.0
        # additional-voltage stepper: the state IS the amplitude, the band is
        # linear in it
        du = br.phase_du_step
        lo = br.phase_du_min_step * du
        hi = br.phase_du_max_step * du
        limits[map.r2cols[j]] = const_share * abs(hi - lo)
      else
        # phase-angle stepper: the band lives in degrees and r2 follows
        # nonlinearly, so evaluate r2 at both band ends and take the span
        α = map.alphas[j]
        r2_at = function (shift_deg::Float64)
          ϕ = -deg2rad(shift_deg)
          d = sin(deg2rad(α) - ϕ)
          abs(d) < 1e-12 ? 0.0 : sin(ϕ) / d
        end
        a = r2_at(br.phase_min_deg)
        b = r2_at(br.phase_max_deg)
        (isfinite(a) && isfinite(b)) || continue
        limits[map.r2cols[j]] = const_share * abs(b - a)
      end
    end
  end
  return limits
end

## initial tap state values in column order (r1/r2 per trafo as released)
function _tap_r_init(map::TapStateMap)::Vector{Float64}
  out = Float64[]
  for j in eachindex(map.branch_idxs)
    map.r1cols[j] != 0 && push!(out, map.r1_0[j])
    map.r2cols[j] != 0 && push!(out, map.r2_0[j])
  end
  return out
end

## per-evaluation overlay: snet branchIdx => (scratch branch at the cascade
## position for the current x, plus its four admittance terms). The scratch
## branches are mutated in place (one per released trafo, serial FD).
function _tap_overlay(x::Vector{Float64}, map::TapStateMap, net::Net)
  d = Dict{Int,NamedTuple{(:branch, :Y11, :Y12, :Y21, :Y22),Tuple{Branch,ComplexF64,ComplexF64,ComplexF64,ComplexF64}}}()
  for j in eachindex(map.branch_idxs)
    r1 = map.r1cols[j] == 0 ? map.r1_0[j] : x[map.r1cols[j]]
    r2 = map.r2cols[j] == 0 ? map.r2_0[j] : x[map.r2cols[j]]
    sb = map.scratch[j]
    ct = _cascade_tap(sb, r1, r2, map.alphas[j])
    sb.tap_ratio = ct.tap_ratio
    sb.phase_shift_deg = ct.phase_shift_deg
    y11, y12, y21, y22 = calcAdmittance(sb, sb.comp.cVN, net.baseMVA)
    d[map.branch_idxs[j]] = (branch = sb, Y11 = y11, Y12 = y12, Y21 = y21, Y22 = y22)
  end
  return d
end

## remove the stamped admittance terms of the released transformers from the
## Ybus (once per run): the predictions re-add them analytically at the
## cascade position, so the Ybus stays constant across FD perturbations.
## Subtracts exactly what createYBUS stamped (calcAdmittance on the LIVE
## branch fields).
function _tap_unstamp!(Ybus::AbstractMatrix{ComplexF64}, net::Net, map::TapStateMap)
  for k in map.branch_idxs
    br = net.branchVec[k]
    y11, y12, y21, y22 = calcAdmittance(br, br.comp.cVN, net.baseMVA)
    f = Int(br.fromBus)
    t = Int(br.toBus)
    Ybus[f, f] -= y11
    Ybus[t, t] -= y22
    Ybus[f, t] -= y12
    Ybus[t, f] -= y21
  end
  return Ybus
end

"""
    _is_machine_transformer(net, k) -> Bool

Machine (generator step-up) transformer detection for the tap-release
guards: one terminal bus of branch `k` hangs on this branch alone (no
other closed branch, no closed link) and carries a generator prosumer.
Estimating such a tap is pointless: the generator terminal voltage behind
it is set by the machine, not observed independently, so the tap state
absorbs whatever the AVR does.
"""
function _is_machine_transformer(net::Net, k::Int)::Bool
  br = net.branchVec[k]
  for side in (Int(br.fromBus), Int(br.toBus))
    others = 0
    for (i, b) in enumerate(net.branchVec)
      i == k && continue
      (Int(b.fromBus) == side || Int(b.toBus) == side) || continue
      _branch_terminal_state(b) == :closed || continue
      others += 1
    end
    for l in net.linkVec
      (Int(l.fromBus) == side || Int(l.toBus) == side) || continue
      l.status == 1 || continue
      others += 1
    end
    others == 0 || continue
    any(ps -> getPosumerBusIndex(ps) == side && isGenerator(ps), net.prosumpsVec) && return true
  end
  return false
end

"""
    _tap_cut_component(net, k) -> Union{Nothing,Set{Int}}

Bridge test for the tap-release guards: connectivity of the two terminal
buses of branch `k` WITHOUT the branch itself, over closed branches and
closed links. Returns `nothing` when the terminals stay connected (the
transformer sits in a loop), otherwise the set of buses in the component
containing `toBus` (the side the transformer alone ties to the rest).
"""
function _tap_cut_component(net::Net, k::Int)::Union{Nothing,Set{Int}}
  br = net.branchVec[k]
  n = length(net.nodeVec)
  adj = [Int[] for _ in 1:n]
  for (i, b) in enumerate(net.branchVec)
    i == k && continue
    _branch_terminal_state(b) == :closed || continue
    push!(adj[Int(b.fromBus)], Int(b.toBus))
    push!(adj[Int(b.toBus)], Int(b.fromBus))
  end
  for l in net.linkVec
    l.status == 1 || continue
    push!(adj[Int(l.fromBus)], Int(l.toBus))
    push!(adj[Int(l.toBus)], Int(l.fromBus))
  end
  target = Int(br.fromBus)
  seen = falses(n)
  stack = [Int(br.toBus)]
  seen[Int(br.toBus)] = true
  while !isempty(stack)
    v = pop!(stack)
    v == target && return nothing
    for w in adj[v]
      seen[w] && continue
      seen[w] = true
      push!(stack, w)
    end
  end
  return Set{Int}(i for i in 1:n if seen[i])
end

"""
    _tap_guard_freezes(net, map, activeMeas) -> (freeze1, freeze2, reasons)

Structural guards on the released taps BEFORE the numerical observability
test: a transformer that is a BRIDGE (its removal cuts the net) whose
cut-off side carries no active voltage-magnitude measurement makes the tap
and the downstream voltage indistinguishable; both regulators of that
transformer are frozen (`:radial_no_voltage_pin`). `reasons[j]` is `:ok`
for untouched transformers.
"""
function _tap_guard_freezes(net::Net, map::TapStateMap, activeMeas::Vector{Measurement})
  nrel = length(map.branch_idxs)
  freeze1 = falses(nrel)
  freeze2 = falses(nrel)
  reasons = fill(:ok, nrel)
  vmBuses = Set{Int}(m.busIdx for m in activeMeas if m.typ == VmMeas && m.busIdx !== nothing)
  for j in eachindex(map.branch_idxs)
    cut = _tap_cut_component(net, map.branch_idxs[j])
    cut === nothing && continue
    if isempty(intersect(vmBuses, cut))
      freeze1[j] = true
      freeze2[j] = true
      reasons[j] = :radial_no_voltage_pin
    end
  end
  return freeze1, freeze2, reasons
end

"""
    _tap_apply_freezes(map, freeze1, freeze2, net, Ybus) -> Union{Nothing,TapStateMap}

Rebuild the tap state map after guard or observability freezes: frozen
regulator columns are removed (the regulator stays at its r0), and a
transformer with NO free column left leaves the map entirely and gets its
model stamp restored in the Ybus (bitwise identical to an unreleased
transformer). Returns `nothing` when nothing remains released. Column
numbers are re-packed from `map.col0`.
"""
function _tap_apply_freezes(map::TapStateMap, freeze1::AbstractVector{Bool}, freeze2::AbstractVector{Bool}, net::Net, Ybus::AbstractMatrix{ComplexF64})::Union{Nothing,TapStateMap}
  idxs = Int[]
  modes = Symbol[]
  alphas = Float64[]
  r1c = Int[]
  r2c = Int[]
  r10 = Float64[]
  r20 = Float64[]
  scratch = Branch[]
  c = map.col0
  for j in eachindex(map.branch_idxs)
    keep1 = map.r1cols[j] != 0 && !freeze1[j]
    keep2 = map.r2cols[j] != 0 && !freeze2[j]
    if !keep1 && !keep2
      # fully frozen: restore the exact stamped terms _tap_unstamp! removed
      # (same calcAdmittance call on the same branch, bitwise inverse)
      br = net.branchVec[map.branch_idxs[j]]
      y11, y12, y21, y22 = calcAdmittance(br, br.comp.cVN, net.baseMVA)
      fb = Int(br.fromBus)
      tb = Int(br.toBus)
      Ybus[fb, fb] += y11
      Ybus[tb, tb] += y22
      Ybus[fb, tb] += y12
      Ybus[tb, fb] += y21
      continue
    end
    push!(idxs, map.branch_idxs[j])
    push!(modes, map.modes[j])
    push!(alphas, map.alphas[j])
    push!(r10, map.r1_0[j])
    push!(r20, map.r2_0[j])
    push!(scratch, map.scratch[j])
    if keep1
      push!(r1c, c)
      c += 1
    else
      push!(r1c, 0)
    end
    if keep2
      push!(r2c, c)
      c += 1
    else
      push!(r2c, 0)
    end
  end
  isempty(idxs) && return nothing
  return TapStateMap(idxs, modes, alphas, r1c, r2c, r10, r20, scratch, map.col0, c - map.col0)
end

## machine terminal of a machine transformer: the degree-1 generator bus
## (0 when branch k is no machine transformer); shared by the guard and the
## post-processing back-calculation
function _machine_side(net::Net, k::Int)::Int
  br = net.branchVec[k]
  for side in (Int(br.fromBus), Int(br.toBus))
    others = 0
    for (i, b) in enumerate(net.branchVec)
      i == k && continue
      (Int(b.fromBus) == side || Int(b.toBus) == side) || continue
      _branch_terminal_state(b) == :closed || continue
      others += 1
    end
    for l in net.linkVec
      (Int(l.fromBus) == side || Int(l.toBus) == side) || continue
      l.status == 1 || continue
      others += 1
    end
    others == 0 || continue
    any(ps -> getPosumerBusIndex(ps) == side && isGenerator(ps), net.prosumpsVec) && return side
  end
  return 0
end

"""
    calcMachineTrafoTapFromSE(net; trafo, v_machine_pu = nothing, p_mw = nothing, q_mvar = nothing) -> NamedTuple

Back-calculate the tap position of a MACHINE (generator step-up)
transformer after a state estimation. Mutually exclusive with a tap
release on the same transformer: a released tap is estimated, a machine
tap is reconstructed here, never both.

A machine terminal is invisible to the estimator (the AVR sets its
voltage, nothing measures it independently), so taking both the terminal
voltage AND the terminal injection from the SE would circularly reproduce
the model tap. The non-circular information is: the estimated NETWORK-side
voltage (from `runse!(...; updateNet = true)`), the known machine terminal
voltage magnitude (the AVR setpoint `vm_pu` of the generator prosumer, or
`v_machine_pu`), the active-power dispatch (prosumer `p` schedule, or
`p_mw`), and the machine reactive power. CAREFUL with Q: under AVR voltage
control the machine's Q is NOT its schedule value; pass the MEASURED
machine reactive power (`q_mvar`, the SCADA telemetry every unit has). The
prosumer `q` value is only a fallback for uncontrolled machines. The
function solves the transformer terminal equation
for the machine angle (active-power balance) nested inside a scan over the
tap fraction (reactive-power balance) and reports the continuous
electrical step plus the nearest mechanical step on the `tap_step`
fraction grid. Reporting only, no write-back.

Errors: transformer released for estimation, not a machine transformer, no
registered SE result, or no machine voltage/injection available.
"""
function calcMachineTrafoTapFromSE(net::Net; trafo, v_machine_pu::Union{Nothing,Float64} = nothing, p_mw::Union{Nothing,Float64} = nothing, q_mvar::Union{Nothing,Float64} = nothing)
  k = _resolve_trafo_branch(net, trafo)
  br = net.branchVec[k]
  br.tap_est_mode == :none || error("calcMachineTrafoTapFromSE: transformer branch $(k) is released for tap estimation; back-calculation and release are mutually exclusive (disable the release first: setTapEstimation!(net; trafo = $(k), enabled = false))")
  m = _machine_side(net, k)
  m != 0 || error("calcMachineTrafoTapFromSE: branch $(k) is not a machine (generator step-up) transformer")
  st = _se_start_state(net)
  st === nothing && error("calcMachineTrafoTapFromSE: no state-estimation result is registered for this net; run runse!(...; updateNet = true) first")
  f = Int(br.fromBus)
  t = Int(br.toBus)
  netBus = m == f ? t : f
  # machine data: AVR setpoint and dispatch, never SE quantities
  gens = [ps for ps in net.prosumpsVec if getPosumerBusIndex(ps) == m && isGenerator(ps)]
  vset = v_machine_pu !== nothing ? v_machine_pu : something(isempty(gens) ? nothing : gens[1].vm_pu, nothing)
  vset === nothing && error("calcMachineTrafoTapFromSE: no machine voltage available; the generator carries no vm_pu setpoint, pass v_machine_pu")
  pset = p_mw !== nothing ? p_mw : sum(something(g.pVal, 0.0) for g in gens; init = 0.0)
  qset = q_mvar !== nothing ? q_mvar : sum(something(g.qVal, 0.0) for g in gens; init = 0.0)
  (p_mw !== nothing || any(g.pVal !== nothing for g in gens)) || error("calcMachineTrafoTapFromSE: no machine injection available; the generator carries no p schedule, pass p_mw")
  Sset = complex(pset, qset) / net.baseMVA
  Vnet = st.vm[netBus] * cis(deg2rad(st.va[netBus]))
  scratch = deepcopy(br)
  # terminal equation at the machine end for a candidate tap fraction r and
  # machine angle theta: S_pred = V_m * conj(I_m), V_m = vset * cis(theta)
  spred = function (r::Float64, θ::Float64)
    ct = _cascade_tap(br, r, 0.0, 0.0)
    scratch.tap_ratio = ct.tap_ratio
    scratch.phase_shift_deg = ct.phase_shift_deg
    y11, y12, y21, y22 = calcAdmittance(scratch, br.comp.cVN, net.baseMVA)
    Vm = vset * cis(θ)
    Im = m == f ? y11 * Vm + y12 * Vnet : y21 * Vnet + y22 * Vm
    return Vm * conj(Im)
  end
  θnet = deg2rad(st.va[netBus])
  # nested solve: for a given r the active power is monotone in the machine
  # angle (power-angle relation across one reactance); bisect theta on the
  # P balance, then scan r on the Q residual
  θ_for_r = function (r::Float64)
    # local names on purpose: the enclosing scope re-uses lo/hi for the
    # ternary search below, and an assignment here would share the boxed
    # variable and destroy that search's window (the serial cousin of the
    # @spawn boxing trap in the project's concurrency rules)
    local θlo = θnet - deg2rad(60.0)
    local θhi = θnet + deg2rad(60.0)
    plo = real(spred(r, θlo)) - real(Sset)
    phi = real(spred(r, θhi)) - real(Sset)
    plo * phi > 0.0 && return nothing   # dispatch outside the angle window
    for _ = 1:60
      mid = (θlo + θhi) / 2
      pm = real(spred(r, mid)) - real(Sset)
      if plo * pm <= 0.0
        θhi = mid
      else
        θlo = mid
        plo = pm
      end
    end
    return (θlo + θhi) / 2
  end
  step1 = br.tap_step > 0.0 ? br.tap_step : 0.00625
  rlo = 1.0 / br.tap_max - 1.0
  rhi = 1.0 / br.tap_min - 1.0
  # coarse scan plus refinement on |Q residual| (the Q balance carries the
  # voltage-ratio information; P is consumed by the angle solve)
  qres = function (r::Float64)
    θ = θ_for_r(r)
    θ === nothing && return Inf
    return abs(imag(spred(r, θ)) - imag(Sset))
  end
  best_r = rlo
  best = Inf
  n = 64
  for i = 0:n
    r = rlo + (rhi - rlo) * i / n
    q = qres(r)
    q < best && (best = q; best_r = r)
  end
  lo = max(rlo, best_r - (rhi - rlo) / n)
  hi = min(rhi, best_r + (rhi - rlo) / n)
  for _ = 1:60
    m1 = lo + (hi - lo) / 3
    m2 = hi - (hi - lo) / 3
    if qres(m1) <= qres(m2)
      hi = m2
    else
      lo = m1
    end
  end
  r = (lo + hi) / 2
  p = r / step1
  pfix = clamp(round(p, RoundNearestTiesAway), ceil(rlo / step1), floor(rhi / step1))
  isfinite(best) || error("calcMachineTrafoTapFromSE: no machine angle satisfies the active-power schedule across the transformer; check p_mw/v_machine_pu against the SE state")
  return (branch = k, name = getCompName(br.comp), machine_bus = m, network_bus = netBus, r_est = r, electrical_step = p, fixed_step = Int(pfix), q_residual_mvar = qres(r) * net.baseMVA, v_machine_pu = vset, converged_se = st.converged)
end

"""
    _fixate_taps(x, map, net) -> (rows, fixedMap, xFixed)

After the estimation the taps are always fixed
to the nearest MECHANICAL step and one more run is solved in which the tap
is no state any more. This helper computes, per released transformer, the
continuous (electrical) and rounded (fixed) step of each regulator on its
OWN grid (regulator 1: the ratio grid `tap_step` with the `tap_min`/`tap_max`
range; regulator 2: the phase grid `phase_step_deg` with its range),
out-of-range flags (the fixed step is clamped into the range), the r values
of the fixed position, a `TapStateMap` with zero state columns (all
regulators frozen at the fixed position), and the state vector without the
tap columns.
"""
function _fixate_taps(x::Vector{Float64}, map::TapStateMap, net::Net)
  rows = NamedTuple[]
  r1f = Float64[]
  r2f = Float64[]
  for j in eachindex(map.branch_idxs)
    br = net.branchVec[map.branch_idxs[j]]
    α = map.alphas[j]
    r1 = map.r1cols[j] == 0 ? map.r1_0[j] : x[map.r1cols[j]]
    r2 = map.r2cols[j] == 0 ? map.r2_0[j] : x[map.r2cols[j]]
    # regulator 1 mechanical grid: the tap FRACTION itself, r1 = n * tap_step
    # (the calcSkewAngleTap convention; CGMES stepVoltageIncrement/100 per
    # step). NOT the multiplier grid ratio * (1 + n * step): the two differ
    # at second order and only the fraction grid reproduces a cascade truth
    # exactly. The [tap_min, tap_max] multiplier band maps through the
    # reciprocal 1/(1+r1) onto the position range.
    step1 = br.tap_step > 0.0 ? br.tap_step : 0.00625
    p1 = r1 / step1
    p1min = (1.0 / br.tap_max - 1.0) / step1
    p1max = (1.0 / br.tap_min - 1.0) / step1
    oor1 = map.r1cols[j] != 0 && (p1 < p1min - 0.5 || p1 > p1max + 0.5)
    # ties away from zero: banker's rounding on a half-step estimate would
    # pick the even neighbour, which has no mechanical meaning
    p1fix = clamp(round(p1, RoundNearestTiesAway), ceil(p1min), floor(p1max))
    # only an ESTIMATED regulator snaps to the grid; a frozen one (no state
    # column) keeps its exact model position, off-grid or not
    r1fix = map.r1cols[j] == 0 ? r1 : p1fix * step1
    if br.phase_du_step > 0.0
      # additional-voltage stepper (Delta-u PST): the mechanical grid is
      # the additional-voltage AMPLITUDE itself, r2 = n * phase_du_step
      # along the nameplate direction; rounding is linear in r2 exactly
      # like regulator 1, and the shift angle merely FOLLOWS (atan). The
      # band lives in Delta-u steps.
      du = br.phase_du_step
      p2 = r2 / du
      p2min = br.phase_du_min_step
      p2max = br.phase_du_max_step
      oor2 = map.r2cols[j] != 0 && (p2 < p2min - 0.5 || p2 > p2max + 0.5)
      p2fix = clamp(round(p2, RoundNearestTiesAway), ceil(p2min), floor(p2max))
      r2fix = map.r2cols[j] == 0 ? r2 : p2fix * du
    else
      # regulator 2 on the phase grid: fixed shift, then the r2 that produces
      # exactly that regulating-vector angle: r2 = sin(phi)/sin(alpha - phi)
      # with phi = angle of the regulating vector = -shift. The
      # [phase_min_deg, phase_max_deg] band is RELATIVE to neutral (the
      # regulating-vector angle passes 0 at the neutral position).
      step2 = br.phase_step_deg > 0.0 ? br.phase_step_deg : 1.25
      s2 = calcSkewAngleTap(tap_fraction = r2, skew_angle_deg = α)
      p2 = s2.effective_shift_deg / step2
      p2min = br.phase_min_deg / step2
      p2max = br.phase_max_deg / step2
      oor2 = map.r2cols[j] != 0 && (p2 < p2min - 0.5 || p2 > p2max + 0.5)
      p2fix = clamp(round(p2, RoundNearestTiesAway), ceil(p2min), floor(p2max))
      r2fix = r2   # frozen regulator: exact model position
      if map.r2cols[j] != 0 && map.modes[j] in (:pst, :both)
        ϕ = -deg2rad(p2fix * step2)
        r2fix = abs(sin(deg2rad(α) - ϕ)) < 1e-12 ? 0.0 : sin(ϕ) / sin(deg2rad(α) - ϕ)
      end
    end
    push!(r1f, r1fix)
    push!(r2f, r2fix)
    push!(rows, (
      branch = map.branch_idxs[j],
      name = getCompName(br.comp),
      mode = map.modes[j],
      alpha_deg = α,
      r1_est = r1,
      r2_est = r2,
      electrical_step_1 = p1,
      fixed_step_1 = Int(p1fix),
      electrical_step_2 = p2,
      fixed_step_2 = Int(p2fix),
      out_of_range = oor1 || oor2,
    ))
  end
  fixedMap = TapStateMap(copy(map.branch_idxs), copy(map.modes), copy(map.alphas), zeros(Int, length(map.branch_idxs)), zeros(Int, length(map.branch_idxs)), r1f, r2f, map.scratch, map.col0, 0)
  # drop the tap columns from the state (they sit before the trailing offset)
  keepCols = [k for k in eachindex(x) if !(map.col0 <= k < map.col0 + map.ncols)]
  return rows, fixedMap, x[keepCols]
end

"""
    setTapEstimation!(net; trafo, mode = :ratio, alpha_deg = nothing, enabled = true)

Release a transformer's tap position(s) as state-estimation states
(0.10.0). `mode = :ratio` estimates the longitudinal regulator r1,
`mode = :pst` the angle regulator r2 along the FIXED nameplate direction
`alpha_deg` (required, typically 30/60/90; never estimated), and
`mode = :both` the two-regulator cascade. `enabled = false` clears the
release. The r0 split (see `_tap_r0_split`) runs at release time as a
consistency assert between the device metadata and the branch state.
Returns `(branch, mode, alpha_deg, r1_0, r2_0)`.
"""
function setTapEstimation!(net::Net; trafo, mode::Symbol = :ratio, alpha_deg::Union{Nothing,Real} = nothing, enabled::Bool = true)
  k = _resolve_trafo_branch(net, trafo)
  br = net.branchVec[k]
  if !enabled
    br.tap_est_mode = :none
    br.tap_est_alpha_deg = 0.0
    return (branch = k, mode = :none, alpha_deg = 0.0, r1_0 = 0.0, r2_0 = 0.0)
  end
  mode in (:ratio, :pst, :both) || error("setTapEstimation!: mode must be :ratio, :pst, or :both (got :$(mode))")
  (br.has_ratio_tap || br.has_phase_tap) || error("setTapEstimation!: branch $(k) carries no tap-changer data")
  α = 0.0
  if mode in (:pst, :both)
    alpha_deg === nothing && error("setTapEstimation!: alpha_deg is required for mode :$(mode) (nameplate regulator direction, e.g. 30/60/90 degrees; it is a Vorgabe and never estimated)")
    α = Float64(alpha_deg)
  end
  r1_0, r2_0 = _tap_r0_split(br, mode, α)
  br.tap_est_mode = mode
  br.tap_est_alpha_deg = α
  return (branch = k, mode = mode, alpha_deg = α, r1_0 = r1_0, r2_0 = r2_0)
end

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
# purpose: transformer tap estimation (0.10.0, step states since 0.60.0):
#          release API (setTapEstimation!), the state map of the released
#          regulators (one step state s per tap changer of the changer set,
#          tap_changer_data.jl), the overlay that stamps a released
#          transformer at its state, the determinability guard, the
#          fixation on the physical steps and the write-back. The WLS loop
#          that consumes these lives in state_estimation.jl.

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
    TapBranchState

One released transformer inside one SE (island) run: the snet branch
index, its changer set (`_effective_tap_changers`), the start position `s0`
of every changer (from the live branch, `_tap_s_from_branch`), the position
`spos` a changer without a state column sits at (`s0` in the estimation,
the fixed step in the fixation run), which changers the release selects
(`released`), the state column per changer (`cols`, 0 = no state), the
freeze reason and note per changer, a scratch `Branch` the overlay rewrites
per evaluation, the live impedance `z0 = (r, x, g_from, b_from, g_to, b_to)`
with its factors `F0` at `s0`, the step size `sigma` per changer
(`_tap_scale`), whether the run applies an impedance source on this branch
(`impedance`: `_tap_reactance_source` is not `:nominal`), the branch
`TapStepImpedance` for `:tsi` sources, and
`tree_phase`, the phase regulators on a tree edge without an angle anchor
that the vertex notice watches, `weak`, the factor of the weak-column
hint per changer (0.0 = no hint, see `_tap_weak_hints!`), and
`rel_ref`, per changer the names of the frozen phase regulators a phase
regulator is only relative to ("" = none, set with note
`:relative_to_frozen` by `_tap_relative_to_frozen!`); such a regulator is
never fixed, it stays a state in the fixed run. `cause` holds per changer
the structural cause the guard found (`_tap_guard_causes!`: `:none`,
`:radial_no_voltage_pin`, `:not_separable`, `:no_angle_anchor`,
`:symmetric_start`; `:zero_column` marks a column the noise test froze).
A cause freezes nothing by itself: it becomes the reason only when the
rank test drops the column (`_tap_rank_freezes!`).
"""
struct TapBranchState
  branch::Int
  set::SETapChangerSet
  s0::Vector{Float64}
  spos::Vector{Float64}
  released::Vector{Bool}
  cols::Vector{Int}
  frozen::Vector{Symbol}
  notes::Vector{Symbol}
  scratch::Branch
  z0::NTuple{6,Float64}
  F0::NTuple{4,Float64}
  sigma::Vector{Float64}
  impedance::Bool
  tsi::Union{Nothing,TapStepImpedance}
  tree_phase::Vector{Bool}
  weak::Vector{Float64}
  rel_ref::Vector{String}
  cause::Vector{Symbol}
end

"""
    TapStateMap

State map of the released tap regulators inside one SE (island) run: the
released transformers (`TapBranchState`, branch order), the first tap state
column `col0` (after the shunt B states, before the PMU offset), the number
of tap states `ncols`, and the `impedance_correction` mode of the run
(`:auto` or `:off`). One state per released, unfrozen regulator, the
electrical step s relative to the changer's neutral step; columns run in
branch order, then in the changer order of the set.
"""
struct TapStateMap
  trafos::Vector{TapBranchState}
  col0::Int
  ncols::Int
  impedance::Symbol
end

# whether a changer of `role` is moved by the release `mode`
_tap_role_released(role::Symbol, mode::Symbol)::Bool = mode === :both || (mode === :ratio && role === :ratio) || (mode === :pst && role === :phase)

# the same state with other columns (and an other overlay position)
_tap_with_cols(ts::TapBranchState, cols::Vector{Int}, spos::Vector{Float64}) = TapBranchState(ts.branch, ts.set, ts.s0, spos, ts.released, cols, ts.frozen, ts.notes, ts.scratch, ts.z0, ts.F0, ts.sigma, ts.impedance, ts.tsi, ts.tree_phase, ts.weak, ts.rel_ref, ts.cause)

# live impedance of a branch in the order the overlay scales it
_tap_live_impedance(br::Branch)::NTuple{6,Float64} = (br.r_pu, br.x_pu, br.g_from_pu, br.b_from_pu, br.g_to_pu, br.b_to_pu)

## Build the map over the (contracted, island) net. Returns `(map, out,
## nochanger)`: the map (`nothing` when no regulator got a state), the
## released transformers that never enter it (every released changer
## frozen with :model_mismatch because the branch is off its changer
## curves; their stamp stays), and the released branches without a changer
## of the released role (one report row each, :no_changer). Columns from
## `col0`; the guards freeze further regulators and re-pack (`_tap_repack`).
function _build_tap_state_map(net::Net, col0::Int, impedance::Symbol)
  states = TapBranchState[]
  out = TapBranchState[]
  nochanger = Int[]
  c = col0
  for (k, br) in enumerate(net.branchVec)
    (br isa Branch && br.tap_est_mode !== :none && _branch_terminal_state(br) == :closed) || continue
    set = _effective_tap_changers(br)
    set === nothing && continue
    # the set as this run takes its impedance (a typed set without the
    # analytic flag under off)
    set = _tap_set_under(set, impedance)
    n = length(set.changers)
    released = Bool[_tap_role_released(ch.role, br.tap_est_mode) for ch in set.changers]
    if !any(released)
      push!(nochanger, k)
      continue
    end
    frozen = fill(:none, n)
    notes = fill(:none, n)
    s0 = _tap_s_from_branch(set, br)
    mismatch = s0 === nothing
    mismatch && (s0 = Float64[ch.step - ch.neutral for ch in set.changers])
    tsi = br.tap_step_impedance
    # per transformer: per-step data always, the analytic correction only
    # without them and under auto
    active = _tap_reactance_source(set, impedance) !== :nominal
    F0 = active ? _tap_impedance_factors(set, s0; s0 = s0, tsi = tsi) : (1.0, 1.0, 1.0, 1.0)
    sigma = Float64[_tap_scale(ch; ratio_n = set.ratio_n) for ch in set.changers]
    cols = zeros(Int, n)
    if mismatch
      for j in 1:n
        released[j] && (frozen[j] = :model_mismatch)
      end
    else
      for j in 1:n
        released[j] || continue
        cols[j] = c
        c += 1
      end
    end
    ts = TapBranchState(k, set, s0, copy(s0), released, cols, frozen, notes, deepcopy(br), _tap_live_impedance(br), F0, sigma, active, tsi, falses(n), zeros(n), fill("", n), fill(:none, n))
    push!(mismatch ? out : states, ts)
  end
  isempty(states) && return (nothing, out, nochanger)
  return (TapStateMap(states, col0, c - col0, impedance), out, nochanger)
end

## Re-pack the columns after freezes: a regulator with a freeze reason
## loses its column (it stays at its start position), a transformer without
## a column left leaves the map and gets its model stamp back in the Ybus
## (the same calcAdmittance call _tap_unstamp! subtracted). Returns the new
## map (`nothing` when nothing is left) and the transformers that left it.
function _tap_repack(map::TapStateMap, net::Net, Ybus::AbstractMatrix{ComplexF64})
  keep = TapBranchState[]
  moved = TapBranchState[]
  c = map.col0
  for ts in map.trafos
    cols = zeros(Int, length(ts.cols))
    for j in eachindex(cols)
      (ts.cols[j] != 0 && ts.frozen[j] === :none) || continue
      cols[j] = c
      c += 1
    end
    if all(==(0), cols)
      _tap_restamp!(Ybus, net, ts.branch)
      push!(moved, _tap_with_cols(ts, cols, ts.spos))
    else
      push!(keep, _tap_with_cols(ts, cols, ts.spos))
    end
  end
  isempty(keep) && return (nothing, moved)
  return (TapStateMap(keep, map.col0, c - map.col0, map.impedance), moved)
end

## (trafo index, changer index, column) of every tap state, column order
function _tap_state_columns(map::TapStateMap)::Vector{Tuple{Int,Int,Int}}
  out = Tuple{Int,Int,Int}[]
  for (j, ts) in enumerate(map.trafos), c in eachindex(ts.cols)
    ts.cols[c] == 0 || push!(out, (j, c, ts.cols[c]))
  end
  return out
end

## start values of the tap states in column order
_tap_s_init(map::TapStateMap)::Vector{Float64} = Float64[map.trafos[j].s0[c] for (j, c, _) in _tap_state_columns(map)]

## step positions of one transformer at the state vector x
_tap_state_s(ts::TapBranchState, x::AbstractVector{Float64})::Vector{Float64} = Float64[ts.cols[c] == 0 ? ts.spos[c] : x[ts.cols[c]] for c in eachindex(ts.cols)]

"""
Per-iteration step limit of each tap state, in steps: a quarter of the
changer's band (`_tap_band`). Measured on a CGMES delivery: without a bound
one weakly determined transformer ran to -329 electrical steps on a band of
about -14 to +18 and took the whole estimation down, but clamping the
STATE into the band was worse than no bound at all (five transformers that
had converged stuck to the lower bound and diverged). The transient
excursion is part of how a weakly determined regulator finds its position,
so what is bounded is the distance one iteration may travel. A band of
zero width gives no limit.
"""
function _tap_step_limits(map::TapStateMap)::Dict{Int,Float64}
  limits = Dict{Int,Float64}()
  for ts in map.trafos, c in eachindex(ts.cols)
    ts.cols[c] == 0 && continue
    smin, smax = _tap_band(ts.set.changers[c])
    smax > smin && (limits[ts.cols[c]] = 0.25 * (smax - smin))
  end
  return limits
end

## scale of every state column for the convergence test: 1 for the voltage,
## shunt and offset states, the step size sigma for a tap state, so `tol`
## keeps its pu meaning (the old r state was s * du on the fraction grid)
function _tap_step_scale(map::Union{Nothing,TapStateMap}, n::Int)::Vector{Float64}
  sc = ones(Float64, n)
  map === nothing && return sc
  for ts in map.trafos, c in eachindex(ts.cols)
    ts.cols[c] == 0 || (sc[ts.cols[c]] = ts.sigma[c])
  end
  return sc
end

## r, x and the four shunt arms of the scratch branch from the live values
## z0 and the impedance factor ratio q = F(s) / F(s0) (g on both g arms,
## b on both b arms, the rule of _apply_tap_step_impedance!); q == 1 gives
## the live impedance bitwise
function _tap_scratch_impedance!(sb::Branch, z::NTuple{6,Float64})
  sb.r_pu = z[1]
  sb.x_pu = z[2]
  set_branch_shunt!(sb; g_from_pu = z[3], b_from_pu = z[4], g_to_pu = z[5], b_to_pu = z[6])
  return sb
end

_tap_z_at(ts::TapBranchState, F::NTuple{4,Float64})::NTuple{6,Float64} = (ts.z0[1] * (F[1] / ts.F0[1]), ts.z0[2] * (F[2] / ts.F0[2]), ts.z0[3] * (F[3] / ts.F0[3]), ts.z0[4] * (F[4] / ts.F0[4]), ts.z0[5] * (F[3] / ts.F0[3]), ts.z0[6] * (F[4] / ts.F0[4]))

## exact stamp of one transformer at the positions s
function _tap_stamp_exact!(ts::TapBranchState, s::Vector{Float64})
  sb = ts.scratch
  R, Phi = _tap_set_ratio_shift(ts.set, s)
  sb.tap_ratio = R
  sb.phase_shift_deg = Phi
  ts.impedance && _tap_scratch_impedance!(sb, _tap_z_at(ts, _tap_impedance_factors(ts.set, s; s0 = ts.s0, tsi = ts.tsi)))
  return sb
end

const _TapOverlayEntry = NamedTuple{(:branch, :Y11, :Y12, :Y21, :Y22),Tuple{Branch,ComplexF64,ComplexF64,ComplexF64,ComplexF64}}

## per-evaluation overlay, exact: snet branch index => (scratch branch at
## the positions of x, its four admittance terms). The scratch branches are
## rewritten in place (one per released transformer, serial FD).
function _tap_overlay(x::Vector{Float64}, map::TapStateMap, net::Net)
  d = Dict{Int,_TapOverlayEntry}()
  for ts in map.trafos
    sb = _tap_stamp_exact!(ts, _tap_state_s(ts, x))
    y11, y12, y21, y22 = calcAdmittance(sb, sb.comp.cVN, net.baseMVA)
    d[ts.branch] = (branch = sb, Y11 = y11, Y12 = y12, Y21 = y21, Y22 = y22)
  end
  return d
end

## Overlay of the FD Jacobian at the base point xbase: a transformer whose
## tap states equal those of xbase is stamped exactly, so every non-tap
## column sees bitwise the stamp of h0; a perturbed one is stamped on the
## analytic tangent, t = t(sbar) + sum_c
## (s_c - sbar_c) dt/ds_c and z = z(sbar) + sum_c (s_c - sbar_c) dz/ds_c. The
## forward difference along it is the chain rule dh/dt (by FD) times dt/ds
## (analytic), the impedance part included when the run applies one.
function _tap_overlay_tangent(x::Vector{Float64}, xbase::Vector{Float64}, map::TapStateMap, net::Net)
  d = Dict{Int,_TapOverlayEntry}()
  for ts in map.trafos
    s = _tap_state_s(ts, x)
    sbar = _tap_state_s(ts, xbase)
    if s == sbar
      sb = _tap_stamp_exact!(ts, s)
    else
      sb = ts.scratch
      R, Phi = _tap_set_ratio_shift(ts.set, sbar)
      t = R * cis(deg2rad(Phi))
      Fb = ts.impedance ? _tap_impedance_factors(ts.set, sbar; s0 = ts.s0, tsi = ts.tsi) : (1.0, 1.0, 1.0, 1.0)
      z = _tap_z_at(ts, Fb)
      for c in eachindex(s)
        ds = s[c] - sbar[c]
        ds == 0.0 && continue
        t += ds * _tap_set_tangent(ts.set, sbar, c)
        if ts.impedance
          dF = _tap_impedance_slope(ts.set, sbar, c; s0 = ts.s0, tsi = ts.tsi)
          dz = (ts.z0[1] * dF[1] / ts.F0[1], ts.z0[2] * dF[2] / ts.F0[2], ts.z0[3] * dF[3] / ts.F0[3], ts.z0[4] * dF[4] / ts.F0[4], ts.z0[5] * dF[3] / ts.F0[3], ts.z0[6] * dF[4] / ts.F0[4])
          z = z .+ ds .* dz
        end
      end
      sb.tap_ratio = abs(t)
      sb.phase_shift_deg = rad2deg(angle(t))
      ts.impedance && _tap_scratch_impedance!(sb, z)
    end
    y11, y12, y21, y22 = calcAdmittance(sb, sb.comp.cVN, net.baseMVA)
    d[ts.branch] = (branch = sb, Y11 = y11, Y12 = y12, Y21 = y21, Y22 = y22)
  end
  return d
end

## remove the stamped admittance terms of the released transformers from the
## Ybus (once per run): the predictions re-add them at the state position,
## so the Ybus stays constant across FD perturbations. Subtracts exactly what
## createYBUS stamped (calcAdmittance on the LIVE branch fields).
function _tap_unstamp!(Ybus::AbstractMatrix{ComplexF64}, net::Net, map::TapStateMap)
  for ts in map.trafos
    br = net.branchVec[ts.branch]
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

## the inverse of _tap_unstamp! for one branch (same call on the same
## branch), so a transformer whose regulators are all frozen is stamped as
## if it had never been released
function _tap_restamp!(Ybus::AbstractMatrix{ComplexF64}, net::Net, k::Int)
  br = net.branchVec[k]
  y11, y12, y21, y22 = calcAdmittance(br, br.comp.cVN, net.baseMVA)
  fb = Int(br.fromBus)
  tb = Int(br.toBus)
  Ybus[fb, fb] += y11
  Ybus[tb, tb] += y22
  Ybus[fb, tb] += y12
  Ybus[tb, fb] += y21
  return Ybus
end

"""
    _is_machine_transformer(net, k) -> Bool

Machine (generator step-up) transformer detection for the tap-release
guards: one terminal of branch `k` hangs on this branch alone and carries
a generator prosumer. The terminal is its closed-link cluster (the buses
closed links tie together count as one node, `_link_clusters`): a
generator behind a closed coupler is at the terminal, and a closed link
is no further connection. Estimating such a tap is pointless: the
generator terminal voltage behind it is set by the machine, not observed
independently, so the tap state absorbs whatever the AVR does. The test
of `_machine_sides`, so the guards, the exclusion and the calculation
after the estimation agree.
"""
_is_machine_transformer(net::Net, k::Int)::Bool = _machine_side(net, k) != 0

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

# whether a changer moves the angle of the branch ratio
function _tap_changes_angle(ch::SETapChanger)::Bool
  k = ch.kind
  k === :asymmetrical && return abs(sind(ch.alpha_deg)) > 1e-9
  (k === :symmetrical || k === :linear) && return true
  k === :tabular && return any(p -> p.angle_deg != 0.0, ch.table)
  return false
end

# whether a changer moves the magnitude of the branch ratio on its band
_tap_changes_magnitude(ch::SETapChanger)::Bool = ch.kind === :ratio || ch.kind === :asymmetrical || (ch.kind === :tabular && any(p -> p.ratio != ch.table[1].ratio, ch.table))

# (Re(f'/f), |f'/f|) of the changer factor f at s: the relative magnitude
# slope dm/m and its length together with the angle slope (radians)
function _tap_factor_log_slope(ch::SETapChanger, s::Float64)::Tuple{Float64,Float64}
  if ch.additive
    # R = R_n (1 + s du) with du = step_abs / R_n
    re = ch.du / (1.0 + s * ch.du)
    return (re, abs(re))
  end
  m, _ = _tap_changer_factor(ch, s)
  dm, da = _tap_changer_slope(ch, s)
  re = dm / m
  return (re, hypot(re, deg2rad(da)))
end

# Symmetric start of a phase regulator at s: the magnitude of
# its factor does not move to first order. Formula kinds: |Re(f'/f)| at
# most 1e-9 |f'/f| (always for :symmetrical and :linear, :asymmetrical at
# alpha 90 and s = 0). Tables: the one-sided slopes of the ratio left and
# right of s have opposite signs or one is zero (a vertex or a flat piece).
function _tap_symmetric_at(ch::SETapChanger, s::Float64)::Bool
  if ch.kind === :tabular
    _, _, dr_right, _ = _tap_table_eval(ch, s)
    _, _, dr_left, _ = _tap_table_eval(ch, s - 1e-9)
    return dr_left * dr_right <= 0.0
  end
  re, ab = _tap_factor_log_slope(ch, s)
  return abs(re) <= 1e-9 * ab
end

# turning point of |t| of an asymmetrical changer, relative to neutral and
# not limited to its band: d|1 + s du e^(j alpha)|/ds = 0 at
# s* = -cos(alpha) / du, and |t|^2 = (s du + cos(alpha))^2 + sin(alpha)^2 is
# symmetric about it, so s and 2 s* - s give the same |t|
_tap_asym_turning_point(ch::SETapChanger)::Float64 = -cosd(ch.alpha_deg) / ch.du

# whether |t| of a changer is not monotonic on its band: then only the start
# side decides which of two positions with the same magnitude the estimate
# takes (note :sign_from_start)
function _tap_magnitude_not_monotonic(ch::SETapChanger)::Bool
  smin, smax = _tap_band(ch)
  if ch.kind === :asymmetrical
    sstar = _tap_asym_turning_point(ch)
    return smin < sstar < smax
  elseif ch.kind === :tabular
    d = [ch.table[k + 1].ratio - ch.table[k].ratio for k in 1:length(ch.table) - 1]
    return any(>(0.0), d) && any(<(0.0), d)
  end
  return false
end

# positions s (relative to neutral) where |t| of a changer turns on its
# band: -cos(alpha)/du for an asymmetrical changer, the table points where
# the ratio changes from rising to falling or back (a flat stretch between
# them counts at its first point)
function _tap_magnitude_vertices(ch::SETapChanger)::Vector{Float64}
  smin, smax = _tap_band(ch)
  if ch.kind === :asymmetrical
    sstar = _tap_asym_turning_point(ch)
    return smin < sstar < smax ? [sstar] : Float64[]
  elseif ch.kind === :tabular
    out = Float64[]
    prev = 0.0
    for k in 1:length(ch.table) - 1
      d = ch.table[k + 1].ratio - ch.table[k].ratio
      d == 0.0 && continue
      (prev != 0.0 && sign(d) != sign(prev)) && push!(out, Float64(ch.table[k].step - ch.neutral))
      prev = d
    end
    return out
  end
  return Float64[]
end

## Columns with note :sign_from_start, the turning points of |t| on their
## band and the changer. The note promises that the estimate stays on the
## side of the start. With noise the best |t| can lie below the turning
## point, where Gauss-Newton jumps across it every iteration (measured on
## the alpha 90 case with s0 = 2 and truth 0, generator sigmas, seeds 2, 5,
## 7 and 9: no convergence in 40 iterations; seeds 3 and 8 converged on the
## other side). The WLS step keeps these states on their side
## (`_tap_keep_side_pins!`).
function _tap_side_vertices(map::Union{Nothing,TapStateMap})::Vector{Tuple{Int,Vector{Float64},SETapChanger}}
  out = Tuple{Int,Vector{Float64},SETapChanger}[]
  map === nothing && return out
  for ts in map.trafos, c in eachindex(ts.cols)
    (ts.cols[c] != 0 && ts.notes[c] === :sign_from_start) || continue
    ch = ts.set.changers[c]
    v = _tap_magnitude_vertices(ch)
    isempty(v) || push!(out, (ts.cols[c], v, ch))
  end
  return out
end

## Steps to hold in this WLS iteration (pairs column => step). A side column
## whose step `dx` reaches or crosses a turning point goes half the way to
## the nearest one (fraction to the boundary: the best |t| on the start side
## lies at the turning point, and the state approaches it without jumping
## across). When the magnitude share of the column at that half-way point
## (`_tap_magnitude_share`, the quantity of the vertex notice) squared is at
## most `jacEps`, the state goes onto the turning point and stays there for
## the rest of the solve (`held`, column => turning point): the information
## that separates the two sides is then below the forward-difference
## tolerance, and halving would not end (measured, seed 5: raw step -22,
## -40, -63 at s 0.5, 0.25, 0.125, then -61 at 0.0625 and +126 at 0.031).
function _tap_keep_side_pins!(held::Dict{Int,Float64}, side::Vector{Tuple{Int,Vector{Float64},SETapChanger}}, x::AbstractVector{Float64}, dx::AbstractVector{Float64}, jacEps::Float64)::Vector{Tuple{Int,Float64}}
  pins = Tuple{Int,Float64}[]
  for (col, verts, ch) in side
    s = x[col]
    if haskey(held, col)
      push!(pins, (col, held[col] - s))
      continue
    end
    snew = s + dx[col]
    best = Inf
    for v in verts
      v == s && continue
      ((v - s) * (v - snew) <= 0.0 && abs(v - s) < abs(best)) && (best = v - s)
    end
    isfinite(best) || continue
    if _tap_magnitude_share(ch, s + 0.5 * best)^2 <= jacEps
      held[col] = s + best
      push!(pins, (col, best))
    else
      push!(pins, (col, 0.5 * best))
    end
  end
  return pins
end

"""
    _tap_guard_causes!(net, map, activeMeas, slackIdx, withVaOffset)

Structural part of the determinability guard. It freezes nothing: it
writes per changer the structural cause (`ts.cause`) that names a rank
deficit when the rank test (`_tap_rank_freezes!`) finds one, plus the
notes of `map`. The first cause found wins, in this order:

1. Bridge rule: a transformer that is a bridge (its removal cuts the net)
   without an active voltage-magnitude measurement on its cut-off side,
   its regulators that move the magnitude of the ratio on their band
   (`_tap_changes_magnitude`, `:radial_no_voltage_pin`): the tap and the
   downstream voltage cannot be told apart. A pure angle changer there
   gets its cause from the angle rule (3).
2. Separability: the released regulators of one transformer, in set order,
   each have to raise the rank of the real 2 x k matrix of
   `(dt/ds_c)/t` at the start; one that does not gets `:not_separable`.
   The ratio plane has two real directions, so a third regulator always
   gets it, and so does the second of two along one direction (two ratio
   changers, a ratio changer and an asymmetrical one at alpha 0 or 180,
   two pure angle changers).
3. Angle anchors: the transformers with a released regulator that moves the
   angle are the edges of a graph over the classes of buses that the other
   closed branches and links tie together; buses with a Va row or at the
   measured end of an Ia row of `activeMeas` are anchors (united, and with
   the slack bus when the run has no PMU offset state; the caller passes
   the rows of the rank test, so an Ia row whose current is zero at the
   start does not count). Walked from the slack class: an edge inside one
   class is pinned by its loop; on a tree edge to
   a new class the phase regulators get `:no_angle_anchor` when the
   transformer has a released ratio regulator, else `:symmetric_start` when
   the magnitude of their factor does not move at the start
   (`_tap_symmetric_at`). A tree-edge phase regulator is watched by the
   vertex notice (`tree_phase`) and gets note `:sign_from_start` when the
   magnitude is not monotonic on its band; the noise test and the rank
   test clear both on the regulators they freeze (`_tap_freeze!`).
   A further edge between two reached classes (a mesh edge) gives its
   phase regulators note `:relative_only`. Which of them are only relative
   to a frozen phase regulator is decided after the freezes
   (`_tap_relative_to_frozen!`).

Until 0.60.0 these rules froze the regulators themselves. The rank test
replaced that: a regulator whose column the other
columns explain is not estimated, and the cause found here only names why.
"""
function _tap_guard_causes!(net::Net, map::TapStateMap, activeMeas::Vector{Measurement}, slackIdx::Int, withVaOffset::Bool)
  live(ts, c) = ts.cols[c] != 0 && ts.frozen[c] === :none
  vmBuses = Set{Int}(m.busIdx for m in activeMeas if m.typ == VmMeas && m.busIdx !== nothing)
  for ts in map.trafos
    # one bridge search over the whole net per transformer: seconds when
    # every branch is a transformer (case6495rte), so an abort point each
    sparlectra_check_abort()
    cut = _tap_cut_component(net, ts.branch)
    (cut === nothing || !isempty(intersect(vmBuses, cut))) && continue
    for c in eachindex(ts.cols)
      # only a regulator that moves |t|: a pure angle changer moves no
      # voltage magnitude, its deficit on a bridge is the missing angle
      # anchor, which the angle pass below names
      (live(ts, c) && _tap_changes_magnitude(ts.set.changers[c])) && (ts.cause[c] = :radial_no_voltage_pin)
    end
  end
  for ts in map.trafos
    basis = nothing
    for c in eachindex(ts.cols)
      (live(ts, c) && ts.cause[c] === :none) || continue
      R, _, dR, dPhi = _tap_set_partials(ts.set, ts.s0, c)
      col = (dR / R, deg2rad(dPhi))
      ncol = hypot(col[1], col[2])
      if basis === nothing
        if ncol > 1e-12
          basis = col
        else
          ts.cause[c] = :not_separable
        end
      elseif basis isa Tuple{Float64,Float64} && abs(basis[1] * col[2] - basis[2] * col[1]) > 1e-9 * hypot(basis[1], basis[2]) * ncol
        # the second direction: the plane is full, every later one gets the cause
        basis = :full
      else
        ts.cause[c] = :not_separable
      end
    end
  end
  _tap_angle_causes!(net, map, activeMeas, slackIdx, withVaOffset)
  return map
end

## Classes of the angle graph as the root bus of every bus: the
## closed branches outside `skip` and the closed links tie their buses
## together; the angle anchors (Va rows, the measured end of Ia rows, same
## PMU time base) are united, and with the slack bus when the run has no
## PMU offset state.
function _tap_angle_classes(net::Net, skip::AbstractSet{Int}, activeMeas::Vector{Measurement}, slackIdx::Int, withVaOffset::Bool)::Vector{Int}
  n = length(net.nodeVec)
  parent = collect(1:n)
  function root(i)
    while parent[i] != i
      parent[i] = parent[parent[i]]
      i = parent[i]
    end
    return i
  end
  unite(a, b) = (ra = root(a); rb = root(b); ra != rb && (parent[ra] = rb); nothing)
  for (i, b) in enumerate(net.branchVec)
    i in skip && continue
    _branch_terminal_state(b) == :closed || continue
    unite(Int(b.fromBus), Int(b.toBus))
  end
  for l in net.linkVec
    l.status == 1 && unite(Int(l.fromBus), Int(l.toBus))
  end
  anchors = Int[]
  for m in activeMeas
    if m.typ == VaMeas && m.busIdx !== nothing
      push!(anchors, Int(m.busIdx))
    elseif m.typ == IaMeas && m.branchIdx !== nothing && 1 <= m.branchIdx <= length(net.branchVec)
      b = net.branchVec[m.branchIdx]
      push!(anchors, m.direction == :to ? Int(b.toBus) : Int(b.fromBus))
    end
  end
  for a in anchors
    unite(a, anchors[1])
    withVaOffset || unite(a, slackIdx)
  end
  return Int[root(i) for i in 1:n]
end

function _tap_angle_causes!(net::Net, map::TapStateMap, activeMeas::Vector{Measurement}, slackIdx::Int, withVaOffset::Bool)
  live(ts, c) = ts.cols[c] != 0 && ts.frozen[c] === :none
  E = Dict{Int,Int}()
  for (j, ts) in enumerate(map.trafos)
    any(c -> live(ts, c) && _tap_changes_angle(ts.set.changers[c]), eachindex(ts.cols)) && (E[ts.branch] = j)
  end
  isempty(E) && return map
  cls = _tap_angle_classes(net, Set{Int}(keys(E)), activeMeas, slackIdx, withVaOffset)
  edges = sort!(collect(keys(E)))
  ends(k) = (cls[Int(net.branchVec[k].fromBus)], cls[Int(net.branchVec[k].toBus)])
  start = cls[slackIdx]
  visited = Set{Int}([start])
  queue = [start]
  tree = Set{Int}()
  while !isempty(queue)
    cl = popfirst!(queue)
    for k in edges
      ca, cb = ends(k)
      (ca == cb || (ca != cl && cb != cl)) && continue
      other = ca == cl ? cb : ca
      other in visited && continue
      push!(visited, other)
      push!(queue, other)
      push!(tree, k)
    end
  end
  for k in edges
    ts = map.trafos[E[k]]
    ca, cb = ends(k)
    ca == cb && continue
    phase = [c for c in eachindex(ts.cols) if live(ts, c) && _tap_changes_angle(ts.set.changers[c])]
    if !(k in tree)
      (ca in visited && cb in visited) || continue
      for c in phase
        ts.notes[c] === :none && (ts.notes[c] = :relative_only)
      end
      continue
    end
    ratio = [c for c in eachindex(ts.cols) if live(ts, c) && ts.set.changers[c].role === :ratio && !(c in phase)]
    for c in phase
      ch = ts.set.changers[c]
      if !isempty(ratio)
        ts.cause[c] === :none && (ts.cause[c] = :no_angle_anchor)
      elseif _tap_symmetric_at(ch, ts.s0[c])
        ts.cause[c] === :none && (ts.cause[c] = :symmetric_start)
      end
      # pinned only by its magnitude effect when the rank test keeps it
      # (also at a table vertex, where the column carries the one-sided
      # slope of the tangent): watched by the vertex notice; the noise
      # test and the rank test clear both on the regulators they freeze
      # (`_tap_freeze!`)
      ts.tree_phase[c] = true
      _tap_magnitude_not_monotonic(ch) && ts.notes[c] === :none && (ts.notes[c] = :sign_from_start)
    end
  end
  return map
end

"""
    _tap_relative_to_frozen!(net, map, frozenStates, activeMeas, slackIdx, withVaOffset)

Last part of the determinability guard, after the rank test, so every
freeze reason counts. A phase regulator that is estimated (no freeze
reason) and not pinned by its own magnitude effect (`tree_phase`) gets its
angle from the
angle difference of the two classes it joins. That difference comes from
the phase regulators whose angle is known: the frozen ones (any freeze
reason, also the transformers of `frozenStates`, the ones that left the
map, :model_mismatch included) and the estimated ones pinned by their
magnitude effect. When the two classes are tied only over a path that runs
through a frozen phase regulator, the measurements fix the regulator only
against the start position of that frozen one: it gets note
`:relative_to_frozen`, `rel_ref` names the frozen transformers on the
shortest such path, and it is not fixed (`_tap_fixable`), its continuous
value stays in the fit. The classes are those of `_tap_angle_classes` with
every transformer that has an estimated or a frozen phase regulator taken
out. Writes the notes of `map` in place; no column changes.
"""
function _tap_relative_to_frozen!(net::Net, map::TapStateMap, frozenStates::Vector{TapBranchState}, activeMeas::Vector{Measurement}, slackIdx::Int, withVaOffset::Bool)
  livePhase(ts) = [c for c in eachindex(ts.cols) if ts.cols[c] != 0 && ts.frozen[c] === :none && _tap_changes_angle(ts.set.changers[c])]
  frozenPhase(ts) = any(c -> ts.released[c] && ts.frozen[c] !== :none && _tap_changes_angle(ts.set.changers[c]), eachindex(ts.cols))
  live = Dict{Int,TapBranchState}()
  frozenSet = Set{Int}()
  for ts in map.trafos
    if !isempty(livePhase(ts))
      live[ts.branch] = ts
    elseif frozenPhase(ts)
      push!(frozenSet, ts.branch)
    end
  end
  for ts in frozenStates
    frozenPhase(ts) && push!(frozenSet, ts.branch)
  end
  (isempty(live) || isempty(frozenSet)) && return map
  cls = _tap_angle_classes(net, union(Set{Int}(keys(live)), frozenSet), activeMeas, slackIdx, withVaOffset)
  ends(k) = (cls[Int(net.branchVec[k].fromBus)], cls[Int(net.branchVec[k].toBus)])
  # estimated phase regulators pinned by their own magnitude effect
  pinned = sort!([k for (k, ts) in live if any(c -> ts.tree_phase[c], livePhase(ts))])
  frozenEdges = sort!(collect(frozenSet))
  # shortest path (in edges) from class a to class b over `edges` without
  # edge `skip`, `nothing` when there is none
  function classPath(a, b, edges, skip)
    prev = Dict{Int,Tuple{Int,Int}}(a => (0, 0))
    queue = [a]
    while !isempty(queue)
      cl = popfirst!(queue)
      cl == b && break
      for k in edges
        k == skip && continue
        ca, cb = ends(k)
        (ca == cb || (ca != cl && cb != cl)) && continue
        other = ca == cl ? cb : ca
        haskey(prev, other) && continue
        prev[other] = (cl, k)
        push!(queue, other)
      end
    end
    haskey(prev, b) || return nothing
    path = Int[]
    cl = b
    while cl != a
      cl, k = prev[cl]
      push!(path, k)
    end
    return path
  end
  for k in sort!(collect(keys(live)))
    ts = live[k]
    phase = [c for c in livePhase(ts) if !ts.tree_phase[c]]
    isempty(phase) && continue
    ca, cb = ends(k)
    ca == cb && continue
    # tied without a frozen one: estimated in its own right
    classPath(ca, cb, pinned, k) === nothing || continue
    path = classPath(ca, cb, vcat(frozenEdges, pinned), k)
    path === nothing && continue
    names = join([getCompName(net.branchVec[j].comp) for j in sort!([j for j in path if j in frozenSet])], ", ")
    for c in phase
      ts.notes[c] = :relative_to_frozen
      ts.rel_ref[c] = names
    end
  end
  return map
end

## Freeze changer `c` of `ts` with `reason` in the guard (noise test or
## rank test). The structural pass ran before both and wrote notes
## (`:sign_from_start`, `:relative_only`) and `tree_phase` for a regulator
## that stays estimated; on a frozen one they would contradict the freeze
## text of its notice, so they go here.
function _tap_freeze!(ts::TapBranchState, c::Int, reason::Symbol)
  ts.frozen[c] = reason
  ts.notes[c] = :none
  ts.tree_phase[c] = false
  return nothing
end

# rounding units of the FD noise bound of the tap column test (below)
const _TAP_FD_NOISE_ULPS = 10.0
# Singular-value bound of the rank test of the tap columns
# (`_tap_rank_freezes!`). It applies AFTER every column of H is scaled to
# norm 1. A column the other columns explain exactly (the
# bridge null vector of a phase regulator: phi + 1 with + 1 on every
# far-side angle; two regulators along one direction) leaves a singular
# value at rounding level, about 1e-16 on the central-difference H, while
# a column that adds rank leaves one far above 1e-8.
# Without the scaling the bound would compare columns of different size
# (MW rows against pu rows, one heavy pseudo measurement), the error of
# the former column rule under another name. The bound does not
# depend on the weights.
const _TAP_RANK_TOL = 1.0e-8
# Work budget of the rank test at the start point (`_tap_rank_freezes!`),
# in units of m*k^2 for m rows and k candidate columns. `_tap_rank_kept`
# keeps one dense m-vector per kept column (Gram-Schmidt, m*k^2) and takes
# the singular values of a growing k x k triangle per candidate (k^4/4,
# measured at about 1/11 of the m*k^2 weight). Measured on case6495rte
# (m = 55561): k = 400 in 6.6 s, k = 800 in 42.7 s. The service releases
# every transformer of a case, and that case stores every branch as a
# transformer with ratio 1: k = 8617, which would run for days and hold
# 3.8 GB of vectors, without an abort point. 1e11 is about one minute on
# that machine; above it no changer is estimated (reason :rank_budget) and
# the run is the run without tap states.
const _TAP_RANK_MAX_WORK = 1.0e11
_tap_rank_work(m::Integer, k::Integer)::Float64 = Float64(m) * Float64(k)^2 + Float64(k)^4 / 11.0

## Noise part of the determinability guard (a regulator no measurement
## sees): at the start point a tap column that is zero within the rounding
## noise of its own forward difference on every row is frozen
## (`:not_observable`, cause `:zero_column` for its notice). Bound: a row
## value h_i is a sum of terms and is computed to about eps(Float64) times
## the size of its terms, not of its value (a zero injection sums terms of
## thousands of MW to zero). The size comes from the voltage columns of the
## same row, in the unit of the row (MW, Mvar, pu):
##   T_i = |z_i| + sum_b (Vm_b |H[i, Vm_b]| + |H[i, Va_b]|)
## For a P or Q injection or flow the term of bus j has the size
## Vm_i Vm_j |Y_ij|, and Vm_j |dh/dVm_j| + |dh/dVa_j| is at least that size
## (cos and sin part of one term); the own-bus term is bounded by |h_i| plus
## the others. The slack has no Va column, so its term counts only with its
## cos part; T_i is an estimate within a small factor, not a strict bound.
## The difference of two evaluations divided by the FD step jac_eps /
## sigma_c gives the entry noise
##   n_i = _TAP_FD_NOISE_ULPS * eps(Float64) * T_i * sigma_c / jac_eps
## on every row where the column is non-zero; the column is dead when
## |H[i, col]| <= n_i on every one of those rows (an all-zero column too).
## The test is per row and does not read the weights: a heavy row elsewhere
## or a heavy noise-only row of the column itself (a zero injection where
## the regulator moves no flow) does not freeze a regulator that one light
## row sees above its noise. `H` is the forward-difference Jacobian at the
## start.
function _tap_column_freezes!(map::TapStateMap, H::AbstractMatrix{<:Real}, z::Vector{Float64}, x::Vector{Float64}, nbus::Int, jacEps::Float64)
  # T_i of the noise bound: Va columns 1:(nbus-1), Vm column of bus b at
  # (nbus-1)+b, Vm_b at the same index of the state x
  termSize = abs.(z)
  # one dense m-vector per voltage column: an abort point every 256 columns
  for col in 1:(nbus - 1)
    col % 256 == 0 && sparlectra_check_abort()
    termSize .+= abs.(H[:, col])
  end
  for b in 1:nbus
    b % 256 == 0 && sparlectra_check_abort()
    termSize .+= abs(x[(nbus - 1) + b]) .* abs.(H[:, (nbus - 1) + b])
  end
  for ts in map.trafos, c in eachindex(ts.cols)
    col = ts.cols[c]
    (col != 0 && ts.frozen[c] === :none) || continue
    sparlectra_check_abort()
    hc = Vector{Float64}(H[:, col])
    all(isfinite, hc) || error("tap column test: Jacobian column $(col) is not finite")
    scale = _TAP_FD_NOISE_ULPS * eps(Float64) * ts.sigma[c] / jacEps
    if all(i -> iszero(hc[i]) || abs(hc[i]) <= scale * termSize[i], eachindex(hc))
      _tap_freeze!(ts, c, :not_observable)
      ts.cause[c] = :zero_column
    end
  end
  return map
end

## Rank core of the tap columns (`_tap_rank_freezes!` at the start point,
## the recheck at the estimate): H with every column scaled to norm 1, the
## columns `stateCols` taken out by sparse least squares (QR), then the
## candidates `cand` in their order; returns per candidate whether it adds
## rank (one more singular value of [P_kept, p] above `_TAP_RANK_TOL`).
## Columns outside both lists do not enter.
function _tap_rank_kept(H::AbstractMatrix{<:Real}, stateCols::Vector{Int}, cand::Vector{Int})::Vector{Bool}
  colScaled(j) = (v = Vector{Float64}(H[:, j]); n = norm(v); n > 0.0 ? v ./ n : v)
  A = sparse(H[:, stateCols])
  cn = [norm(view(A, :, j)) for j in axes(A, 2)]
  A = A * Diagonal([n > 0.0 ? 1.0 / n : 0.0 for n in cn])
  F = qr(A)
  Q = Vector{Float64}[]
  T = zeros(0, 0)
  rankT = 0
  kept = falses(length(cand))
  for (i, j) in enumerate(cand)
    # one abort point per candidate: the loop is quadratic in the number of
    # candidates and was the step a Web UI abort could not stop
    sparlectra_check_abort()
    b = colScaled(j)
    p = b .- A * (F \ b)
    # twice Gram-Schmidt against the kept parts (orthonormal), for rounding
    coef = zeros(length(Q))
    for _ in 1:2
      for (k, q) in enumerate(Q)
        a = dot(q, p)
        coef[k] += a
        p .-= a .* q
      end
    end
    d = norm(p)
    na = length(Q)
    M = zeros(na + 1, na + 1)
    M[1:na, 1:na] .= T
    M[1:na, na + 1] .= coef
    M[na + 1, na + 1] = d
    rankM = count(>(_TAP_RANK_TOL), svdvals(M))
    if rankM > rankT && d > 0.0
      push!(Q, p ./ d)
      T = M
      rankT = rankM
      kept[i] = true
    end
  end
  return collect(kept)
end

## Rank test of the tap columns at the start point: H with every column scaled to norm 1; a tap column that
## adds no rank to the state columns (Va, Vm, shunt B, PMU offset) and to
## the tap columns kept before it is not estimated. Computed in two parts so
## that it scales to large nets: the state columns are taken out by sparse
## least squares (QR), which leaves the part p of each tap column outside
## their span (rank [A B] = rank A + rank P for P the parts of B); a
## candidate is kept when the singular values of [P_kept, p] count one more
## above `_TAP_RANK_TOL` than those of P_kept. They are computed on the
## small triangular factor [T c; 0 d] of [P_kept, p] (T the factor of the
## kept parts, c and d the coefficients and the length of p against them),
## which has the same singular values. Columns run in map order (branch
## order, then the changer order of the set), the ones without a
## structural cause before the ones with: of columns that share one null
## vector the last one is dropped, which is the one a structural rule
## names when one applies (else the last parallel member). A dropped regulator
## is frozen with its structural cause (`ts.cause`, `_tap_guard_causes!`),
## else with the cause of a phase regulator of another transformer that
## joins the same two angle classes (a parallel member), else with
## `:not_observable`. `H` is the central-difference Jacobian at the start
## (the caller says why).
function _tap_rank_freezes!(map::TapStateMap, H::AbstractMatrix{<:Real}, net::Net, activeMeas::Vector{Measurement}, slackIdx::Int, withVaOffset::Bool)
  live(ts, c) = ts.cols[c] != 0 && ts.frozen[c] === :none
  tapCols = Set{Int}(ts.cols[c] for ts in map.trafos for c in eachindex(ts.cols) if ts.cols[c] != 0)
  cand = [(ts, c) for ts in map.trafos for c in eachindex(ts.cols) if live(ts, c)]
  isempty(cand) && return map
  # the columns without a structural cause first: of two columns that
  # share one dependence the later one is dropped, so the drop lands on
  # the changer the structural rule names (the phase changer of a ratio and
  # phase pair on a bridge, the tree edge of parallel phase changers)
  cand = vcat([x for x in cand if x[1].cause[x[2]] === :none], [x for x in cand if x[1].cause[x[2]] !== :none])
  # above the work budget the test is not run and no changer is estimated
  # (_TAP_RANK_MAX_WORK): the run continues without tap states instead of
  # running for hours
  work = _tap_rank_work(size(H, 1), length(cand))
  if work > _TAP_RANK_MAX_WORK
    @warn "runse!: tap estimation: $(length(cand)) released regulator(s) on $(size(H, 1)) measurement row(s) exceed the work budget of the rank test ($(round(work; sigdigits = 3)) against $(_TAP_RANK_MAX_WORK)); no regulator is estimated, every one keeps its transferred step (reason rank_budget)"
    for (ts, c) in cand
      _tap_freeze!(ts, c, :rank_budget)
    end
    return map
  end
  kept = _tap_rank_kept(H, [j for j in axes(H, 2) if !(j in tapCols)], Int[ts.cols[c] for (ts, c) in cand])
  dropped = Tuple{TapBranchState,Int}[cand[i] for i in eachindex(cand) if !kept[i]]
  isempty(dropped) && return map
  # angle classes with every transformer of the map that carries a phase
  # regulator taken out: parallel members join the same pair of classes
  phaseTrafos = Set{Int}(ts.branch for ts in map.trafos if any(c -> ts.cols[c] != 0 && _tap_changes_angle(ts.set.changers[c]), eachindex(ts.cols)))
  cls = _tap_angle_classes(net, phaseTrafos, activeMeas, slackIdx, withVaOffset)
  pairOf(ts) = (br = net.branchVec[ts.branch]; minmax(cls[Int(br.fromBus)], cls[Int(br.toBus)]))
  for (ts, c) in dropped
    reason = ts.cause[c]
    if reason === :none && _tap_changes_angle(ts.set.changers[c])
      for other in map.trafos
        (other.branch != ts.branch && other.branch in phaseTrafos && pairOf(other) == pairOf(ts)) || continue
        k = findfirst(j -> other.cols[j] != 0 && other.cause[j] !== :none && _tap_changes_angle(other.set.changers[j]), eachindex(other.cols))
        k === nothing && continue
        reason = other.cause[k]
        break
      end
    end
    _tap_freeze!(ts, c, reason === :none ? :not_observable : reason)
  end
  return map
end

## note `:conditional` on the estimated ratio regulators of a transformer
## whose phase regulator the rank test dropped with `:no_angle_anchor` while
## it moves |t| on its band (kind b, or a table with varying ratio): the
## ratio step holds only together with the frozen phase position
function _tap_conditional_notes!(map::TapStateMap)
  for ts in map.trafos
    any(c -> ts.frozen[c] === :no_angle_anchor && _tap_changes_magnitude(ts.set.changers[c]), eachindex(ts.cols)) || continue
    for c in eachindex(ts.cols)
      (ts.cols[c] != 0 && ts.frozen[c] === :none && ts.set.changers[c].role === :ratio) && (ts.notes[c] = :conditional)
    end
  end
  return map
end

## Note `:conditional` of the ratio regulators of one transformer after the
## recheck at the estimate: set while a phase regulator of it that moves
## |t| is kept out with `:no_angle_anchor` or withdrawn by the recheck (the
## ratio step then holds only together with that frozen phase position),
## cleared when no such regulator is left (a phase regulator estimated again)
function _tap_recheck_conditional!(ts::TapBranchState)
  held = any(c -> ts.frozen[c] in (:no_angle_anchor, :not_observable_at_estimate, :mirror_step_fits, :weak_at_estimate, :not_converged_at_estimate) && ts.set.changers[c].role === :phase && _tap_changes_magnitude(ts.set.changers[c]), eachindex(ts.cols))
  for c in eachindex(ts.cols)
    (ts.released[c] && ts.frozen[c] === :none && ts.set.changers[c].role === :ratio) || continue
    if held
      ts.notes[c] === :none && (ts.notes[c] = :conditional)
    elseif ts.notes[c] === :conditional
      ts.notes[c] = :none
    end
  end
  return nothing
end

## The rows of the rank test at the state V (overlay `tapOv`, shunt B
## `sb`): a current row (Imag, Ia) counts only when its predicted current
## passes 3 sigma of its magnitude, the gate the WLS loop applies to the Ia
## rows (sigma of the paired Imag row, else `floorA`; an Imag row its own
## sigma). Below, the magnitude has a kink at zero current and the angle no
## value, so no difference gives a derivative: at a flat start an Imag row
## of size 1e-4 made a bridge tap look observable
function _tap_rank_rows(net::Net, activeMeas::Vector{Measurement}, V::Vector{ComplexF64}, sb, tapOv, floorA::Float64)::Vector{Int}
  function gate(m)
    m.typ in (ImagMeas, IaMeas) || return true
    pair = m.typ == ImagMeas ? m : findfirst(q -> q.typ == ImagMeas && q.branchIdx == m.branchIdx && q.busIdx == m.busIdx && q.direction == m.direction, activeMeas)
    thr = 3.0 * (pair === nothing ? floorA : (pair isa Measurement ? pair : activeMeas[pair]).sigma)
    return abs(_end_current_A(net, m, V, sb, tapOv)) >= thr
  end
  return Int[i for (i, m) in enumerate(activeMeas) if gate(m)]
end

# J difference below which the mirror step fits as well as the fixed step
# (recheck at the estimate): 9 is the
# 3 sigma threshold of chi2 with one degree of freedom
const _TAP_MIRROR_DJ = 9.0

# sigma_s bound (steps) of a changer the recheck at the estimate estimates
# again: at or above it the estimate is not taken over (:weak_at_estimate).
# The rank test there has no weight bound; a bridge transformer without a
# voltage on its cut-off side gains rank from the losses of its load
# alone, and its step is then known only to several steps
const _TAP_RECHECK_SIGMA_MAX = 0.5

## The mirror step of changer `ch` at its fixed step `s` (relative to
## neutral) for the mirror check of the recheck at the estimate: a
## symmetrical phase changer mirrors about its neutral step (-s, the angle
## changes sign, |t| stays 1), an asymmetrical one with alpha not 0 or 180
## about the turning point s* of |t| (2 s* - s, `_tap_asym_turning_point`).
## Without an angle anchor the measurements see such a changer only through
## quantities that are the same on both steps. Rounded to a physical step;
## `nothing` for the other kinds, when the mirror lies outside the band or
## rounds onto `s` itself (no second step).
function _tap_mirror_step(ch::SETapChanger, s::Float64)::Union{Nothing,Float64}
  ch.additive && return nothing
  if ch.kind === :symmetrical
    m = -s
  elseif ch.kind === :asymmetrical && ch.du !== nothing && ch.du != 0.0 && abs(sind(ch.alpha_deg)) > 1e-9
    m = 2.0 * _tap_asym_turning_point(ch) - s
  else
    return nothing
  end
  smin, smax = _tap_band(ch)
  (smin - 0.5 <= m <= smax + 0.5) || return nothing
  mf = _tap_fixed_position(ch, m)
  return mf == s ? nothing : mf
end

## A copy of `ts` for the recheck at the estimate with the columns `cols`
## and the positions `spos`: its own freeze, note, cause and hint vectors
## (the noise test writes into them), so the verdicts do not touch the
## states of the run
_tap_recheck_copy(ts::TapBranchState, cols::Vector{Int}, spos::Vector{Float64}) = TapBranchState(ts.branch, ts.set, ts.s0, spos, ts.released, cols, fill(:none, length(cols)), copy(ts.notes), ts.scratch, ts.z0, ts.F0, ts.sigma, ts.impedance, ts.tsi, copy(ts.tree_phase), zeros(length(cols)), copy(ts.rel_ref), copy(ts.cause))

## The map of the recheck at the estimate: over the fixed run (`fin` at the state `x`, in the
## transformer order of the estimation map `est`) and the transformers that
## left the map at the start (`moved`). Columns from `fin.col0`, in
## transformer and changer order: an estimated and fixed changer (kind
## `:est`, at its fixed step), a released changer the start test kept out
## with a rank reason (`:frz`, at its transferred step; of `moved` only the
## transformers with such a changer), and a state of the fixed run (a
## `:relative_to_frozen` regulator, at its value; a state column, no
## candidate). Returns `(rmap, xr, cand, movedIdx)`: the map of
## `_tap_recheck_copy` states, the state vector with these columns in
## place of the tap block of `x`, the candidates `(j, c, col, kind)` (`j`
## the transformer in `rmap`) and the indices into `moved` of the
## transformers appended after those of `fin`.
function _tap_recheck_map(fin::TapStateMap, est::TapStateMap, x::Vector{Float64}, moved::Vector{TapBranchState})
  trafos = TapBranchState[]
  cand = NamedTuple{(:j, :c, :col, :kind),Tuple{Int,Int,Int,Symbol}}[]
  vals = Float64[]
  movedIdx = Int[]
  col = fin.col0
  frz(ts, c) = ts.released[c] && ts.cols[c] == 0 && ts.frozen[c] in _TAP_START_RANK_REASONS
  function add!(ts, kinds)
    cols = zeros(Int, length(ts.cols))
    for c in eachindex(cols)
      kinds[c] === :none && continue
      cols[c] = col
      push!(vals, kinds[c] === :state ? x[ts.cols[c]] : ts.spos[c])
      kinds[c] === :state || push!(cand, (j = length(trafos) + 1, c = c, col = col, kind = kinds[c]))
      col += 1
    end
    push!(trafos, _tap_recheck_copy(ts, cols, copy(ts.spos)))
  end
  for (j, ts) in enumerate(fin.trafos)
    e = est.trafos[j]
    add!(ts, Symbol[ts.cols[c] != 0 ? :state : (e.cols[c] != 0 && ts.frozen[c] === :none && _tap_fixable(e, c)) ? :est : frz(ts, c) ? :frz : :none for c in eachindex(ts.cols)])
  end
  for (i, ts) in enumerate(moved)
    any(c -> frz(ts, c), eachindex(ts.cols)) || continue
    push!(movedIdx, i)
    add!(ts, Symbol[frz(ts, c) ? :frz : :none for c in eachindex(ts.cols)])
  end
  xr = vcat(x[1:(fin.col0 - 1)], vals, x[(fin.col0 + fin.ncols):end])
  return TapStateMap(trafos, fin.col0, col - fin.col0, fin.impedance), xr, cand, movedIdx
end

## Hint without effect on the freeze: the
## comparison with the strongest voltage column on the UNWEIGHTED
## sensitivity H. A column still estimated with norm(H[:, col]) / sigma_c
## <= jac_eps * max_b norm(H[:, Vm_b]) gets the factor max / (norm /
## sigma_c) in `ts.weak` (the row notice says "weak column, factor x below
## the strongest"). Before 0.60.0 the weighted form of this comparison froze
## the regulator; with weights one heavy pseudo measurement set the
## reference. Returns one (branch, changer, factor, bus of the strongest Vm
## column) per hint for the log line of the caller.
function _tap_weak_hints!(map::TapStateMap, H::AbstractMatrix{<:Real}, nbus::Int, jacEps::Float64)::Vector{NTuple{4,Any}}
  vmMax = 0.0
  vmBus = 0
  for b in 1:nbus
    v = norm(H[:, (nbus - 1) + b])
    v > vmMax && (vmMax = v; vmBus = b)
  end
  hints = NTuple{4,Any}[]
  for ts in map.trafos, c in eachindex(ts.cols)
    col = ts.cols[c]
    (col != 0 && ts.frozen[c] === :none) || continue
    strength = norm(H[:, col]) / ts.sigma[c]
    if strength <= jacEps * vmMax
      ts.weak[c] = vmMax / strength
      push!(hints, (ts.branch, c, ts.weak[c], vmBus))
    end
  end
  return hints
end

## Columns of the parallel-pair prior and their notes: released
## transformers between the same two buses whose difference no flow or
## current row (Pflow, Qflow, Imag, Ia) on a member can tell apart share an
## unobservable direction; their columns get a weak prior toward the start
## in the WLS loop and the regulators note `:relative_only`.
function _tap_prior_columns!(map::TapStateMap, net::Net, activeMeas::Vector{Measurement})::Vector{Int}
  key(ts) = (br = net.branchVec[ts.branch]; minmax(Int(br.fromBus), Int(br.toBus)))
  count = Dict{Tuple{Int,Int},Int}()
  for ts in map.trafos
    count[key(ts)] = get(count, key(ts), 0) + 1
  end
  measured = Set{Tuple{Int,Int}}()
  for ts in map.trafos
    any(m -> m.branchIdx == ts.branch && m.typ in (PflowMeas, QflowMeas, ImagMeas, IaMeas), activeMeas) && push!(measured, key(ts))
  end
  cols = Int[]
  for ts in map.trafos
    (count[key(ts)] > 1 && !(key(ts) in measured)) || continue
    for c in eachindex(ts.cols)
      ts.cols[c] == 0 && continue
      push!(cols, ts.cols[c])
      ts.notes[c] === :none && (ts.notes[c] = :relative_only)
    end
  end
  return cols
end

## closed-link clusters of the net: the representative (smallest bus index)
## of every bus, over the links with status 1; a bus without a closed link
## is its own cluster. The links have no impedance, so the buses of a
## cluster share one voltage, the node the SE link contraction fuses them
## into (_se_prepare).
function _link_clusters(net::Net)::Vector{Int}
  root = collect(1:length(net.nodeVec))
  find(i) = (while root[i] != i; root[i] = root[root[i]]; i = root[i]; end; i)
  for l in net.linkVec
    l.status == 1 || continue
    a = find(Int(l.fromBus))
    b = find(Int(l.toBus))
    a == b && continue
    a < b ? (root[b] = a) : (root[a] = b)
  end
  return [find(i) for i in eachindex(root)]
end

## machine terminal of a machine transformer: the terminal bus of branch k
## whose closed-link cluster hangs on this branch alone and carries a
## generator (0 when branch k is no machine transformer); shared by the
## guard and the post-processing back-calculation
_machine_side(net::Net, k::Int)::Int = get(_machine_sides(net, (k,)), k, 0)

## the machine terminal of every transformer in `ks` in one pass over the
## net: a terminal bus whose closed-link cluster (`_link_clusters`, one
## node) has no closed branch but this one and carries a generator on any
## of its buses. A generator behind a closed coupler is at the terminal,
## a closed link is no further connection. Transformers without such a
## terminal are left out.
function _machine_sides(net::Net, ks)::Dict{Int,Int}
  root = _link_clusters(net)
  deg = zeros(Int, length(net.nodeVec))
  for b in net.branchVec
    _branch_terminal_state(b) == :closed || continue
    deg[root[Int(b.fromBus)]] += 1
    deg[root[Int(b.toBus)]] += 1
  end
  gens = Set{Int}(root[getPosumerBusIndex(ps)] for ps in net.prosumpsVec if isGenerator(ps))
  out = Dict{Int,Int}()
  for k in ks
    br = net.branchVec[k]
    own = _branch_terminal_state(br) == :closed ? 1 : 0
    for side in (Int(br.fromBus), Int(br.toBus))
      if deg[root[side]] - own == 0 && root[side] in gens
        out[k] = side
        break
      end
    end
  end
  return out
end

## the buses of the closed-link cluster of bus `m` (m alone without links)
function _link_cluster_buses(net::Net, m::Int)::Vector{Int}
  root = _link_clusters(net)
  return [i for i in eachindex(root) if root[i] == root[m]]
end

## Why transformers stay out of the tap estimation by the configuration
## `te` (`state_estimation.tap_estimation`): `reasons` maps a branch index to
## `:excluded` (named in `te.exclude` by name or CGMES mRID; the user's
## choice wins over the machine rule) or `:machine_transformer` (generator
## step-up, unless `te.include_machine_transformers`), for the transformers
## in `ks` (default: every transformer); `unknown` lists the entries of
## `te.exclude` that name no transformer of the net.
function _tap_exclusions(net::Net, te::TapEstimationConfig; ks = (k for (k, br) in enumerate(net.branchVec) if br.ratio != 0.0))
  reasons = Dict{Int,Symbol}()
  ks = collect(Int, ks)
  unknown = String[]
  if !isempty(te.exclude)
    mrids = _transformer_mrids(net)
    names = Set{String}()
    for (k, br) in enumerate(net.branchVec)
      br.ratio != 0.0 || continue
      push!(names, getCompName(br.comp))
      haskey(mrids, k) && push!(names, mrids[k])
    end
    wanted = Set{String}(te.exclude)
    unknown = String[e for e in unique(te.exclude) if !(e in names)]
    for k in ks
      (getCompName(net.branchVec[k].comp) in wanted || (haskey(mrids, k) && mrids[k] in wanted)) && (reasons[k] = :excluded)
    end
  end
  if !te.include_machine_transformers
    for k in keys(_machine_sides(net, ks))
      haskey(reasons, k) || (reasons[k] = :machine_transformer)
    end
  end
  return (reasons = reasons, unknown = unknown)
end

## The transformers `te` keeps out of a run on `net`, by component id of the
## branch (reasons of `_tap_exclusions` for the released changers); empty
## with the switch off. Resolved once on the caller's net and handed to
## every solve of it: the link contraction of `_se_prepare` folds closed
## links into their representative bus and renames it, which changes the
## mRID match (keyed on bus names), so a solve on the contracted copy (the
## diagnostics) would otherwise keep out other transformers than `runse!`
## on the same net. The machine test sees the same node on both, it treats
## a closed-link cluster as one (`_machine_sides`).
function _tap_exclusion_map(net::Net, te::TapEstimationConfig)::Dict{String,Symbol}
  te.enabled || return Dict{String,Symbol}()
  ks = [k for (k, b) in enumerate(net.branchVec) if b.ratio != 0.0 && b.tap_est_mode !== :none]
  return Dict{String,Symbol}(net.branchVec[k].comp.cID => r for (k, r) in _tap_exclusions(net, te; ks = ks).reasons)
end

## the warning on exclude entries that name no transformer of `net`, once
## per call of the public entries (`who` names it)
function _tap_warn_unknown_exclusions(net::Net, te::TapEstimationConfig, who::AbstractString)
  isempty(te.exclude) && return nothing
  unknown = _tap_exclusions(net, te; ks = Int[]).unknown
  isempty(unknown) || @warn "$(who): state_estimation.tap_estimation.exclude names $(join(repr.(unknown), ", ")), which is no transformer of this net; the entry has no effect"
  return nothing
end

## Freeze the released changers of the transformers the configuration keeps
## out (`reasons` by component id of the branch, see `_tap_exclusions`):
## the reason replaces any other, they keep their start position.
function _tap_mark_exclusions!(states, net::Net, reasons::Dict{String,Symbol})
  isempty(reasons) && return nothing
  for ts in states
    r = get(reasons, net.branchVec[ts.branch].comp.cID, :none)
    r === :none && continue
    for c in eachindex(ts.frozen)
      ts.released[c] && (ts.frozen[c] = r)
    end
  end
  return nothing
end

## The rows of machine transformers (reason `:machine_transformer`) after a
## run that wrote its state into `net` (update_net): the ratio regulator
## row gets the step calculated from the AVR setpoint, the dispatch and the
## measured machine P and Q of the machine node (`_machine_tap_back_calc`),
## `status = :calculated`; nothing is written back. Without a measured Q on
## every bus of the machine node the row stays `:not_estimated` with reason
## `:q_not_measured_in_group` and its transferred step, and one info line
## names the transformer and the buses; without a measured P the scheduled
## P is used and the row gets note `:p_from_schedule`. A row whose
## calculation fails stays `:not_estimated`, its notice says why.
function _tap_machine_rows(rows::Vector{NamedTuple}, net::Net, measurements::Vector{Measurement})::Vector{NamedTuple}
  any(r -> r.reason === :machine_transformer, rows) || return rows
  vm = Float64[something(nd._vm_pu, 1.0) for nd in net.nodeVec]
  va = Float64[something(nd._va_deg, 0.0) for nd in net.nodeVec]
  out = NamedTuple[]
  for r in rows
    if r.reason !== :machine_transformer || r.role !== :ratio
      push!(out, r)
      continue
    end
    # one back-calculation per machine transformer (402 on case6495rte)
    sparlectra_check_abort()
    k = r.branch
    # the machine node is the closed-link cluster of the machine terminal:
    # the GSU is its only branch, so the sum of the measured injections of
    # its buses (the first active row per bus) is what the machine feeds
    # into the transformer. Only a fully measured cluster gives that sum,
    # the rule of the SE link contraction (_se_prepare): a partial sum
    # misses the unmeasured buses (a zero-injection row at the coupler bus
    # alone would read as a machine at zero)
    cl = _link_cluster_buses(net, _machine_side(net, k))
    hasrow(typ, b) = any(mm -> mm.active && mm.typ == typ && mm.busIdx == b, measurements)
    injsum = function (typ)
      s = 0.0
      for b in cl
        i = findfirst(mm -> mm.active && mm.typ == typ && mm.busIdx == b, measurements)
        i === nothing && return nothing
        s += measurements[i].value
      end
      return s
    end
    qm = injsum(QinjMeas)
    if qm === nothing
      # under voltage control the machine Q is not its schedule, a step
      # from the scheduled Q would be invented: the transferred step stays
      buses = join((net.nodeVec[b].comp.cName for b in cl if !hasrow(QinjMeas, b)), ", ")
      @info "tap estimation: machine transformer $(r.name) not calculated, Q is not measured at $(buses) of its machine node; the transferred step stays"
      push!(out, merge(r, (reason = :q_not_measured_in_group, notice = _tap_assumed_notice(r, string(_tap_notice_text(:q_not_measured_in_group), " (", buses, ")")))))
      continue
    end
    pm = injsum(PinjMeas)
    src = pm === nothing ? "the machine telemetry for Q and the generator schedule for P (P not measured on every bus of the machine node)" : "the machine telemetry"
    bt = try
      _machine_tap_back_calc(net, k, vm, va, true; p_mw = pm, q_mvar = qm)
    catch err
      # the named failures of the back-calculation (no ratio changer, no
      # machine voltage or dispatch, no machine angle meets the dispatch)
      # are ErrorExceptions; anything else is a bug and propagates
      err isa ErrorException || rethrow()
      push!(out, merge(r, (notice = string(r.notice, "; not calculated: ", sprint(showerror, err)),)))
      continue
    end
    # only the first ratio changer is calculated; another ratio row of the
    # same transformer keeps its start position
    if bt.regulator != r.regulator
      push!(out, r)
      continue
    end
    s_fix = Float64(bt.fixed_step - r.neutral_step)
    push!(out, merge(r, (
      s_est = bt.s_est,
      du_est = r.du === nothing ? nothing : bt.s_est * r.du,
      ratio_est = bt.ratio_est,
      angle_est_deg = bt.angle_est_deg,
      s_fixed = s_fix,
      step_fixed = Float64(bt.fixed_step),
      deviation_steps = Float64(bt.fixed_step) - r.step_model,
      status = :calculated,
      note = pm === nothing ? :p_from_schedule : r.note,
      notice = _tap_assumed_notice(r, string("calculated after the estimation from the AVR setpoint and ", src, ", not estimated (Q residual ", round(bt.q_residual_mvar; sigdigits = 3), " MVar)")),
    )))
  end
  return out
end

## the reasons of a regulator the determinability guard kept out at the
## start point (rank test with its structural causes, the noise test), and
## with the four of the recheck at the estimate (`_tap_recheck_map`)
const _TAP_START_RANK_REASONS = (:radial_no_voltage_pin, :not_separable, :no_angle_anchor, :symmetric_start, :not_observable)
const _TAP_RANK_REASONS = (_TAP_START_RANK_REASONS..., :not_observable_at_estimate, :mirror_step_fits, :weak_at_estimate, :not_converged_at_estimate)

## the reasons of a regulator the run kept out automatically (rank test,
## noise test, fallback without tap estimation): its row gets the step
## calculated from the control target as a comparison (`_tap_target_rows`)
const _TAP_TARGET_REASONS = (_TAP_RANK_REASONS..., :fallback_not_converged, :fallback_worse_fit)

## The rows after a run that wrote its state into `net` (update_net) with
## the comparison column `step_from_target` (for the automatically
## excluded transformers): for a regulator kept out
## automatically (`_TAP_TARGET_REASONS`) the step that meets the control
## target of its enabled tap controller (`_tap_step_from_target`), rounded
## to the nearest physical step. It is not applied: the regulator keeps
## its transferred step, the target is an assumption (a controller at a
## limit, in hand, with another target). When it differs from the
## transferred step by more than one step the notice says so, and the SE
## service writes the notice to run.log. Without a control target the
## column stays empty. Excluded and machine transformers are not touched.
function _tap_target_rows(rows::Vector{NamedTuple}, net::Net)::Vector{NamedTuple}
  any(r -> r.regulator > 0 && r.status === :not_estimated && r.reason in _TAP_TARGET_REASONS, rows) || return rows
  ctrls = PowerTransformerControl[c for c in _tap_controllers(net) if c.enabled]
  isempty(ctrls) && return rows
  V = ComplexF64[something(nd._vm_pu, 1.0) * cis(deg2rad(something(nd._va_deg, 0.0))) for nd in net.nodeVec]
  out = NamedTuple[]
  for r in rows
    if !(r.regulator > 0 && r.status === :not_estimated && r.reason in _TAP_TARGET_REASONS)
      push!(out, r)
      continue
    end
    st = _tap_step_from_target(net, r.branch, r.regulator, ctrls, V)
    if st === nothing
      push!(out, r)
      continue
    end
    far = abs(st - r.step_model) > 1.0 + 1e-9
    text = string("step from the control target ", round(Int, st), " differs from the transferred step ", round(Int, r.step_model), " by more than one step (not applied)")
    push!(out, merge(r, (step_from_target = st, notice = far ? join(filter(!isempty, [r.notice, text]), "; ") : r.notice)))
  end
  return out
end

## The absolute step of changer `c` of branch `k` that meets the control
## target of its tap controller in `ctrls` (matched by the transformer
## name, id or branch index as `addPowerTransformerControl!` resolves it,
## the controller of the changer's role; a ratio follower of a master/slave
## group takes the master's voltage target) at the estimated voltages `V` (complex pu, net order), the other
## changers at their transferred steps, the impedance at its live value;
## `nothing` when there is no target for the changer's role or it cannot
## be evaluated locally. Ratio changer, voltage control: the voltage at the
## other terminal and the complex power through the regulated terminal stay
## at their estimate, the terminal voltage follows from the branch equation
## (fixed point), and the step is the one whose |V| meets the target; a
## target bus that is no terminal counts through its estimated difference to
## the terminal on its side, which needs the transformer to be a bridge.
## Phase changer, active-power control over this transformer: both terminal
## voltages stay at their estimate and the step is the one whose flow meets
## the target; the network would share part of the change, so the step is
## a lower bound of the move. The search scans the band (65 points) and
## refines on a third of the interval as `_machine_tap_back_calc` does; the
## result is the nearest physical step (`_tap_fixed_position`).
function _tap_step_from_target(net::Net, k::Int, c::Int, ctrls::Vector{PowerTransformerControl}, V::Vector{ComplexF64})::Union{Nothing,Float64}
  br = net.branchVec[k]
  set = _effective_tap_changers(br)
  (set === nothing || !(1 <= c <= length(set.changers))) && return nothing
  ch = set.changers[c]
  # the controller is chosen per role: a transformer may carry a ratio
  # controller and a phase controller (one per TapChangerControl in the
  # CGMES import), and the first one by name need not be the one of this
  # changer
  ownsTrafo(x) = br.comp.cName == x.trafo || br.comp.cID == x.trafo || string(br.branchIdx) == x.trafo
  fitsRole(x) = ch.role === :ratio ? (x.control_ratio && x.mode in (:voltage, :voltage_and_branch_active_power) && x.target_bus !== nothing && x.target_vm_pu !== nothing) : (x.control_phase && x.mode in (:branch_active_power, :voltage_and_branch_active_power) && x.target_branch !== nothing && x.p_target_mw !== nothing)
  ci = findfirst(x -> ownsTrafo(x) && fitsRole(x), ctrls)
  # a ratio follower of a master/slave group has no controller of its own;
  # it regulates the master's target bus
  if ci === nothing && ch.role === :ratio
    isFollower(x) = any(fn -> br.comp.cName == fn || br.comp.cID == fn || string(br.branchIdx) == fn, x.followers)
    ci = findfirst(x -> fitsRole(x) && isFollower(x), ctrls)
  end
  ci === nothing && return nothing
  ctrl = ctrls[ci]
  s0 = _tap_s_from_branch(set, br)
  s0 === nothing && return nothing
  f = Int(br.fromBus)
  t = Int(br.toBus)
  busOf(name) = findfirst(nd -> getCompName(nd.comp) == name, net.nodeVec)
  scratch = deepcopy(br)
  function yAt(s::Float64)
    sx = copy(s0)
    sx[c] = s
    R, Phi = _tap_set_ratio_shift(set, sx)
    scratch.tap_ratio = R
    scratch.phase_shift_deg = Phi
    return calcAdmittance(scratch, br.comp.cVN, net.baseMVA)
  end
  objective = nothing
  if ch.role === :ratio
    (ctrl.control_ratio && ctrl.mode in (:voltage, :voltage_and_branch_active_power) && ctrl.target_bus !== nothing && ctrl.target_vm_pu !== nothing) || return nothing
    tb = busOf(ctrl.target_bus)
    tb === nothing && return nothing
    g = if tb == f || tb == t
      tb
    else
      cut = _tap_cut_component(net, k)
      cut === nothing ? 0 : (tb in cut ? t : f)
    end
    g == 0 && return nothing
    o = g == f ? t : f
    y0 = yAt(s0[c])
    Ig = g == f ? y0[1] * V[f] + y0[2] * V[t] : y0[3] * V[f] + y0[4] * V[t]
    Sg = V[g] * conj(Ig)
    target = ctrl.target_vm_pu - (abs(V[tb]) - abs(V[g]))
    objective = function (s::Float64)
      y = yAt(s)
      ygg, ygo = g == f ? (y[1], y[2]) : (y[4], y[3])
      v = V[g]
      for _ in 1:200
        vn = (conj(Sg / v) - ygo * V[o]) / ygg
        abs(vn - v) < 1e-13 && return abs(abs(vn) - target)
        v = vn
      end
      # no fixed point at this step: the load does not fit through the
      # branch there
      return Inf
    end
  else
    (ctrl.control_phase && ctrl.mode in (:branch_active_power, :voltage_and_branch_active_power) && ctrl.target_branch !== nothing && ctrl.p_target_mw !== nothing) || return nothing
    a = busOf(ctrl.target_branch[1])
    b = busOf(ctrl.target_branch[2])
    (a !== nothing && b !== nothing && minmax(a, b) == minmax(f, t)) || return nothing
    objective = function (s::Float64)
      y = yAt(s)
      Ia = a == f ? y[1] * V[f] + y[2] * V[t] : y[3] * V[f] + y[4] * V[t]
      return abs(real(V[a] * conj(Ia)) * net.baseMVA - ctrl.p_target_mw)
    end
  end
  slo, shi = _tap_band(ch)
  best_s = slo
  best = Inf
  n = 64
  for i = 0:n
    s = slo + (shi - slo) * i / n
    q = objective(s)
    q < best && (best = q; best_s = s)
  end
  isfinite(best) || return nothing
  lo = max(slo, best_s - (shi - slo) / n)
  hi = min(shi, best_s + (shi - slo) / n)
  for _ = 1:60
    m1 = lo + (hi - lo) / 3
    m2 = hi - (hi - lo) / 3
    objective(m1) <= objective(m2) ? (hi = m2) : (lo = m1)
  end
  return Float64(ch.neutral) + _tap_fixed_position(ch, (lo + hi) / 2)
end

"""
    calcMachineTrafoTapFromSE(net; trafo, v_machine_pu = nothing, p_mw = nothing, q_mvar = nothing) -> NamedTuple

Back-calculate the tap position of a MACHINE (generator step-up)
transformer after a state estimation. Mutually exclusive with a tap
release on the same transformer: a released tap is estimated, a machine
tap is reconstructed here, never both. A machine transformer has a
terminal that connects to nothing but this transformer and carries a
generator; the buses closed links tie together count as one terminal, so
a generator behind a closed coupler qualifies, and its setpoint and
schedule are read from the generators on those buses.

A machine terminal is invisible to the estimator (the AVR sets its
voltage, nothing measures it independently), so taking both the terminal
voltage AND the terminal injection from the SE would circularly reproduce
the model tap. The non-circular information is: the estimated NETWORK-side
voltage (from `runse!` with `state_estimation.update_net = true`), the known machine terminal
voltage magnitude (the AVR setpoint `vm_pu` of the generator prosumer, or
`v_machine_pu`), the active-power dispatch (prosumer `p` schedule, or
`p_mw`), and the machine reactive power. CAREFUL with Q: under AVR voltage
control the machine's Q is NOT its schedule value; pass the MEASURED
machine reactive power (`q_mvar`, the SCADA telemetry every unit has). The
prosumer `q` value is only a fallback for uncontrolled machines. The
function solves the transformer terminal equation for the machine angle
(active-power balance) nested inside a scan over the step `s` of the
branch's first ratio tap changer (reactive-power balance; the other
changers stay at their start position, the impedance at the live value).
Returns `s_est` (continuous step relative to neutral), `electrical_step =
neutral + s_est` and `fixed_step`, the nearest physical step on the
changer's own band, the changer index `regulator` and the branch ratio and
angle at `s_est`. Reporting only, no write-back. `runse!` runs the same
calculation by itself for every released machine transformer that
`state_estimation.tap_estimation.include_machine_transformers = false` keeps
out of the estimation (`tapEstimates` row `status = :calculated`).

Errors: transformer released for estimation, not a machine transformer, no
ratio tap changer, no registered SE result, or no machine voltage/injection
available.
"""
function calcMachineTrafoTapFromSE(net::Net; trafo, v_machine_pu::Union{Nothing,Float64} = nothing, p_mw::Union{Nothing,Float64} = nothing, q_mvar::Union{Nothing,Float64} = nothing)
  k = _resolve_trafo_branch(net, trafo)
  br = net.branchVec[k]
  br.tap_est_mode == :none || error("calcMachineTrafoTapFromSE: transformer branch $(k) is released for tap estimation; back-calculation and release are mutually exclusive (disable the release first: setTapEstimation!(net; trafo = $(k), enabled = false))")
  st = _se_start_state(net)
  st === nothing && error("calcMachineTrafoTapFromSE: no state-estimation result is registered for this net; run runse! with state_estimation.update_net = true first")
  return _machine_tap_back_calc(net, k, st.vm, st.va, st.converged; v_machine_pu = v_machine_pu, p_mw = p_mw, q_mvar = q_mvar)
end

## The back-calculation of `calcMachineTrafoTapFromSE` on the bus voltages
## `vm` (pu) and `va` (degrees) of an estimate, for runse! (the rows of
## machine transformers) and the public function. `converged` is passed
## through as `converged_se`. Throws an ErrorException for each named failure
## (not a machine transformer, no ratio changer, no machine voltage or
## injection, no machine angle meets the dispatch).
function _machine_tap_back_calc(net::Net, k::Int, vm::AbstractVector{Float64}, va::AbstractVector{Float64}, converged::Bool; v_machine_pu::Union{Nothing,Float64} = nothing, p_mw::Union{Nothing,Float64} = nothing, q_mvar::Union{Nothing,Float64} = nothing)
  br = net.branchVec[k]
  m = _machine_side(net, k)
  m != 0 || error("calcMachineTrafoTapFromSE: branch $(k) is not a machine (generator step-up) transformer")
  set = _effective_tap_changers(br)
  cr = set === nothing ? nothing : findfirst(ch -> ch.role === :ratio, set.changers)
  cr === nothing && error("calcMachineTrafoTapFromSE: branch $(k) has no ratio tap changer")
  ch = set.changers[cr]
  s0 = _tap_s_from_branch(set, br)
  s0 === nothing && (s0 = Float64[c.step - c.neutral for c in set.changers])
  f = Int(br.fromBus)
  t = Int(br.toBus)
  netBus = m == f ? t : f
  # machine data: AVR setpoint and dispatch, never SE quantities; the
  # generators on any bus of the machine terminal's closed-link cluster
  cl = Set{Int}(_link_cluster_buses(net, m))
  gens = [ps for ps in net.prosumpsVec if getPosumerBusIndex(ps) in cl && isGenerator(ps)]
  vset = v_machine_pu !== nothing ? v_machine_pu : something(isempty(gens) ? nothing : gens[1].vm_pu, nothing)
  vset === nothing && error("calcMachineTrafoTapFromSE: no machine voltage available; the generator carries no vm_pu setpoint, pass v_machine_pu")
  pset = p_mw !== nothing ? p_mw : sum(something(g.pVal, 0.0) for g in gens; init = 0.0)
  qset = q_mvar !== nothing ? q_mvar : sum(something(g.qVal, 0.0) for g in gens; init = 0.0)
  (p_mw !== nothing || any(g.pVal !== nothing for g in gens)) || error("calcMachineTrafoTapFromSE: no machine injection available; the generator carries no p schedule, pass p_mw")
  Sset = complex(pset, qset) / net.baseMVA
  Vnet = vm[netBus] * cis(deg2rad(va[netBus]))
  scratch = deepcopy(br)
  # terminal equation at the machine end for a candidate step s and machine
  # angle theta: S_pred = V_m * conj(I_m), V_m = vset * cis(theta)
  spred = function (s::Float64, θ::Float64)
    sx = copy(s0)
    sx[cr] = s
    R, Phi = _tap_set_ratio_shift(set, sx)
    scratch.tap_ratio = R
    scratch.phase_shift_deg = Phi
    y11, y12, y21, y22 = calcAdmittance(scratch, br.comp.cVN, net.baseMVA)
    Vm = vset * cis(θ)
    Im = m == f ? y11 * Vm + y12 * Vnet : y21 * Vnet + y22 * Vm
    return Vm * conj(Im)
  end
  θnet = deg2rad(va[netBus])
  # nested solve: for a given s the active power is monotone in the machine
  # angle (power-angle relation across one reactance); bisect theta on the
  # P balance, then scan s on the Q residual
  θ_for_s = function (s::Float64)
    # local names on purpose: the enclosing scope re-uses lo/hi for the
    # ternary search below, and an assignment here would share the boxed
    # variable and destroy that search's window (the serial cousin of the
    # @spawn boxing trap in the project's concurrency rules)
    local θlo = θnet - deg2rad(60.0)
    local θhi = θnet + deg2rad(60.0)
    plo = real(spred(s, θlo)) - real(Sset)
    phi = real(spred(s, θhi)) - real(Sset)
    plo * phi > 0.0 && return nothing   # dispatch outside the angle window
    for _ = 1:60
      mid = (θlo + θhi) / 2
      pm = real(spred(s, mid)) - real(Sset)
      if plo * pm <= 0.0
        θhi = mid
      else
        θlo = mid
        plo = pm
      end
    end
    return (θlo + θhi) / 2
  end
  slo, shi = _tap_band(ch)
  # coarse scan plus refinement on |Q residual| (the Q balance carries the
  # voltage-ratio information; P is consumed by the angle solve)
  qres = function (s::Float64)
    θ = θ_for_s(s)
    θ === nothing && return Inf
    return abs(imag(spred(s, θ)) - imag(Sset))
  end
  best_s = slo
  best = Inf
  n = 64
  for i = 0:n
    s = slo + (shi - slo) * i / n
    q = qres(s)
    q < best && (best = q; best_s = s)
  end
  lo = max(slo, best_s - (shi - slo) / n)
  hi = min(shi, best_s + (shi - slo) / n)
  for _ = 1:60
    m1 = lo + (hi - lo) / 3
    m2 = hi - (hi - lo) / 3
    if qres(m1) <= qres(m2)
      hi = m2
    else
      lo = m1
    end
  end
  s = (lo + hi) / 2
  isfinite(best) || error("calcMachineTrafoTapFromSE: no machine angle satisfies the active-power schedule across the transformer; check p_mw/v_machine_pu against the SE state")
  sx = copy(s0)
  sx[cr] = s
  R, Phi = _tap_set_ratio_shift(set, sx)
  return (branch = k, name = getCompName(br.comp), machine_bus = m, network_bus = netBus, regulator = cr, s_est = s, electrical_step = ch.neutral + s, fixed_step = round(Int, ch.neutral + _tap_fixed_position(ch, s)), ratio_est = R, angle_est_deg = Phi, q_residual_mvar = qres(s) * net.baseMVA, v_machine_pu = vset, converged_se = converged)
end

# Ratio grid positions of the branch fields, counted from `ratio`: the
# MATPOWER and SCF writers, the SCF reader and the scenario tap_pos op use
# them; the tap estimation does not since 0.60.0 (it works on the changer
# sets of tap_changer_data.jl).
# - a typed ratio model (taps_derived): the model steps, evaluated by the
#   resolver (_model_step_value), as the controllers move them;
# - `tap_ratio_grid` :linear (CGMES RatioTapChanger): the ratio n steps
#   away on the grid the tap controller steps on (_ratio_grid_shift:
#   ratio + n * tap_step, or the n-th table point);
# - `tap_ratio_grid` :fraction (applyTapNameplate!, constructor default):
#   the fraction r1 = n * tap_step of t = ratio / (1 + r1).
_tap_step1(br::Branch) = br.tap_step > 0.0 ? br.tap_step : 0.00625

# the typed ratio model of the branch when it has a step increment,
# `nothing` on the legacy grids
function _tap_r1_model(br::Branch)
  m = _derived_tap_model(br, :ratio)
  return (m === nothing || m.tapStepPercent == 0.0) ? nothing : m
end

function _tap_r1_at(br::Branch, n::Real)::Float64
  m = _tap_r1_model(br)
  m === nothing || return br.ratio / _model_step_value(br, :ratio, m.step + Int(n)) - 1.0
  br.tap_ratio_grid === :linear || return n * _tap_step1(br)
  return br.ratio / _ratio_grid_shift(br, br.ratio, n; window = (-Inf, Inf)) - 1.0
end

# (position, ratio) of every table point a controller reaches from `ratio`,
# position 0 being `ratio` itself, ordered by position
function _tap_table_positions(br::Branch)::Vector{Tuple{Int,Float64}}
  out = [(0, br.ratio)]
  for dir in (1, -1)
    r = br.ratio
    n = 0
    while true
      nxt = _ratio_grid_shift(br, r, dir; window = (-Inf, Inf))
      nxt == r && break
      n += dir
      r = nxt
      push!(out, (n, r))
    end
  end
  return sort!(out; by = first)
end

# continuous position of the fraction r1 on the grid of _tap_r1_at. A typed
# model is linear in its step; between two table points the position
# interpolates linearly, beyond the last one it continues with the last
# interval.
function _tap_r1_position(br::Branch, r1::Real)::Float64
  m = _tap_r1_model(br)
  br.tap_ratio_grid === :linear || m !== nothing || return r1 / _tap_step1(br)
  t = br.ratio / (1.0 + r1)
  if m !== nothing
    t0 = _model_step_value(br, :ratio, m.step)
    return (t - t0) / (_model_step_value(br, :ratio, m.step + 1) - t0)
  end
  br.tap_ratio_points === nothing && return (t - br.ratio) / br.tap_step
  pos = _tap_table_positions(br)
  length(pos) == 1 && return (t - br.ratio) / _tap_step1(br)
  k = findfirst(p -> p[2] >= t, pos)
  k = k === nothing ? length(pos) : max(k, 2)
  (n1, t1), (n2, t2) = pos[k - 1], pos[k]
  return n1 + (n2 - n1) * (t - t1) / (t2 - t1)
end

# position band, lower end first: the model range of a typed model, else
# the positions of the range ends tap_min/tap_max (absolute ratios on every
# legacy grid)
function _tap_r1_band(br::Branch)::Tuple{Float64,Float64}
  m = _tap_r1_model(br)
  m === nothing || return (Float64(min(m.lowStep, m.highStep) - m.step), Float64(max(m.lowStep, m.highStep) - m.step))
  return minmax(_tap_r1_position(br, br.ratio / br.tap_max - 1.0), _tap_r1_position(br, br.ratio / br.tap_min - 1.0))
end

# nearest whole position to `p` inside the band; ties away from zero, the
# band ends with a tolerance for the rounding of their own computation
function _tap_r1_fixed_position(br::Branch, p::Real)::Float64
  pmin, pmax = _tap_r1_band(br)
  return clamp(round(p, RoundNearestTiesAway), ceil(pmin - 1e-9), floor(pmax + 1e-9))
end

"""
    _tap_fixed_position(ch, s) -> Float64

Physical step nearest to the continuous estimate `s` (relative to the
changer's neutral step) on the changer's own absolute band (rounded to the
nearest physical step, up or down): the formula
kinds round with ties away from zero and clamp into `low - neutral ..
high - neutral`; a phase table takes its nearest point; a ratio table takes
its nearest point inside the table range and a whole step on the formula
outside it. A tie between two points goes to the one farther from neutral.
"""
function _tap_fixed_position(ch::SETapChanger, s::Real)::Float64
  smin, smax = _tap_band(ch)
  if ch.kind === :tabular
    cands = Float64[p.step - ch.neutral for p in ch.table]
    if ch.outside_du !== nothing
      lo, hi = cands[1], cands[end]
      for q in ceil(Int, smin - 1e-9):floor(Int, smax + 1e-9)
        (q < lo || q > hi) && push!(cands, Float64(q))
      end
    end
    best = cands[1]
    for q in cands
      d, db = abs(q - s), abs(best - s)
      (d < db || (d == db && abs(q) > abs(best))) && (best = q)
    end
    return best
  end
  return clamp(round(Float64(s), RoundNearestTiesAway), ceil(smin - 1e-9), floor(smax + 1e-9))
end

## whether tap state column c of `ts` is fixed to a physical step: every
## estimated regulator except one with note `:relative_to_frozen`, which
## stays a state
_tap_fixable(ts::TapBranchState, c::Int)::Bool = ts.cols[c] != 0 && ts.notes[c] !== :relative_to_frozen

## the tap state columns of `map` that the fixation removes
_tap_fixable_columns(map::TapStateMap)::Vector{Int} = Int[ts.cols[c] for ts in map.trafos for c in eachindex(ts.cols) if _tap_fixable(ts, c)]

## distance of a state to its physical step (`_tap_fixed_position`); the
## sequential fixation fixes the column with the smallest distance first
_tap_step_distance(ch::SETapChanger, s::Float64)::Float64 = abs(s - _tap_fixed_position(ch, s))

## The fixable column of `map` at `x` closest to its physical step (the
## lowest column on a tie), 0 when none is left
function _tap_next_fix(x::Vector{Float64}, map::TapStateMap)::Int
  best = 0
  bestD = Inf
  for ts in map.trafos, c in eachindex(ts.cols)
    _tap_fixable(ts, c) || continue
    d = _tap_step_distance(ts.set.changers[c], x[ts.cols[c]])
    d < bestD && (best = ts.cols[c]; bestD = d)
  end
  return best
end

"""
    _tap_fix_columns(x, map, fix) -> (fixedMap, xFixed, colmap)

Fixation of the tap states in the columns `fix`: each goes to
its nearest physical step (`_tap_fixed_position`), one frozen with
`:out_of_range` to its transferred step, and leaves the state; the
other tap columns stay states, renumbered from `map.col0` in their order.
Returns the map of the next run (every transformer of `map` in the same
order, a fixed regulator stamped at its step through `spos`), the state
vector without the fixed columns and the old-to-new number of every kept
tap column. With all columns in `fix` the next run has no tap state (the
fixed run); a single column per call is one step of the sequential
fixation.
"""
function _tap_fix_columns(x::Vector{Float64}, map::TapStateMap, fix::AbstractSet{Int})
  trafos = TapBranchState[]
  colmap = Dict{Int,Int}()
  next = map.col0
  for ts in map.trafos
    spos = copy(ts.spos)
    cols = zeros(Int, length(ts.cols))
    for c in eachindex(ts.cols)
      col = ts.cols[c]
      col == 0 && continue
      if col in fix
        # an estimate too far outside the band (_tap_mark_out_of_range!)
        # keeps the transferred step
        spos[c] = ts.frozen[c] === :out_of_range ? ts.s0[c] : _tap_fixed_position(ts.set.changers[c], x[col])
      else
        cols[c] = next
        colmap[col] = next
        next += 1
      end
    end
    push!(trafos, _tap_with_cols(ts, cols, spos))
  end
  # the tap columns sit before the trailing offset; drop the fixed ones
  keepCols = [k for k in eachindex(x) if !(k in fix)]
  return TapStateMap(trafos, map.col0, next - map.col0, map.impedance), x[keepCols], colmap
end

## Per transformer of the estimation map `mapEst` (the maps of later runs
## keep its order): the continuous estimate at `xEst`, the position of the
## final run at `(xFin, mapFin)` (a fixed step, the continuous value of a
## regulator that stayed a state, the start of a frozen one) and the
## out-of-range flag of the estimate, `(s_est, s_fixed, oor)`
function _tap_positions(xEst::Vector{Float64}, mapEst::TapStateMap, xFin::Vector{Float64}, mapFin::TapStateMap)
  out = Tuple{Vector{Float64},Vector{Float64},Vector{Bool}}[]
  for (j, ts) in enumerate(mapEst.trafos)
    s_est = _tap_state_s(ts, xEst)
    oor = falses(length(s_est))
    for c in eachindex(s_est)
      ts.cols[c] == 0 && continue
      smin, smax = _tap_band(ts.set.changers[c])
      oor[c] = s_est[c] < smin - 0.5 || s_est[c] > smax + 0.5
    end
    push!(out, (s_est, _tap_state_s(mapFin.trafos[j], xFin), collect(oor)))
  end
  return out
end

## Report rows of one released transformer, one per released
## changer: the estimate `s_est` (continuous, relative to neutral; the start
## position when frozen), the delta U s_est * du (`nothing` without du), the
## branch ratio and angle at the estimate, the position of the fixed run
## `s_fixed`, model step, fixed step and deviation in steps, the
## out-of-range flag, whether the regulator was fixed, its status
## (`:estimated`, `:out_of_range` or `:not_estimated`) with the reason,
## note and notice. An estimated changer with a partner the rank test or
## the recheck kept out, or one with status `:out_of_range`, gets the
## partner notice (`_tap_partner_text`). A regulator with note `:relative_to_frozen` is not
## fixed: `s_fixed` is its continuous value in the fixed run. `converged` enables the notices of the estimate (the
## vertex notice of a phase regulator on a tree edge without an angle
## anchor, the re-import guard of a CGMES changer). `sigma_s` holds per
## changer the standard deviation of the step from the covariance of the
## estimate (`nothing` = none, see `_se_state_sigmas`); the row carries it
## and the distance of `s_est` from its nearest physical step in units of
## it (`round_distance_sigma`, a hint only, the rounding does not read it).
## `step_from_target` is filled later for the regulators the run kept out
## automatically (`_tap_target_rows`). `voltage_hint` (`_tap_voltage_hint`)
## goes into the notice of a changer the rank test kept out.
## `extra_notice` holds per changer a text the notice ends with (the J of
## the mirror check of the recheck at the estimate, the estimate and sigma_s
## of a weak re-estimate, "" = none). `weak` holds per changer the
## continuous estimate (relative to neutral) and sigma_s of a weak
## re-estimate (`nothing` = none): a changer frozen with :weak_at_estimate
## shows them in s_est, du_est and sigma_s, its fixed step stays the
## transferred one. ratio_est and angle_est_deg of every row stay at the
## steps the state was solved with, so a partner row is not shifted by it.
function _tap_regulator_rows(ts::TapBranchState, br::Branch, ob::Int, mrid::AbstractString, s_est::Vector{Float64}, s_fix::Vector{Float64}, oor::AbstractVector{Bool}, fixed::Bool, converged::Bool; sigma_s::AbstractVector = fill(nothing, length(ts.cols)), voltage_hint::AbstractString = "", extra_notice::AbstractVector{<:AbstractString} = fill("", length(ts.cols)), weak::AbstractVector = fill(nothing, length(ts.cols)))::Vector{NamedTuple}
  rows = NamedTuple[]
  R, Phi = _tap_set_ratio_shift(ts.set, s_est)
  isWeak = [ts.frozen[c] === :weak_at_estimate && weak[c] !== nothing for c in eachindex(ts.cols)]
  # the value a row shows in s_est: the weak re-estimate replaces only its
  # own entry, after the composite R, Phi above
  s_row = [isWeak[c] ? first(weak[c]) : s_est[c] for c in eachindex(s_est)]
  for c in eachindex(ts.set.changers)
    ts.released[c] || continue
    ch = ts.set.changers[c]
    estimated = ts.cols[c] != 0 && ts.frozen[c] === :none
    relative = ts.notes[c] === :relative_to_frozen
    vertex = (converged && estimated && ts.tree_phase[c] && _tap_magnitude_share(ch, s_est[c]) < 0.1) ? :vertex : :none
    guard = :none
    if fixed && estimated && ts.set.origin === :cgmes && !ch.additive
      m, _ = _tap_changer_factor(ch, s_fix[c])
      (0.5 <= m <= 2.0) || (guard = :reimport_guard)
    end
    # an estimate too far outside the band, not taken over
    # (_tap_mark_out_of_range!): its sigma_s stays in the row
    outr = ts.frozen[c] === :out_of_range
    sg = (estimated || outr) ? sigma_s[c] : isWeak[c] ? last(weak[c]) : nothing
    fixedRow = fixed && estimated && !relative
    oorText = ""
    if oor[c]
      v, d = _tap_band_verdict(ch, s_est[c], sg)
      oorText = _tap_out_of_range_text(ch, s_est[c]; verdict = outr ? :out_of_range : (fixedRow && v === :band_end) ? :band_end : :inside, dist = d, s0 = ts.s0[c])
    end
    # an estimated changer whose partner on the same transformer is not
    # determinable or not taken over (out_of_range): its value holds for the
    # partner's transferred step only
    partner = estimated ? join([_tap_partner_text(ts, k) for k in eachindex(ts.set.changers) if k != c && ts.released[k] && (ts.frozen[k] in _TAP_RANK_REASONS || ts.frozen[k] === :out_of_range)], "; ") : ""
    push!(rows, (
      branch = ob,
      regulator = c,
      name = getCompName(br.comp),
      changer = ch.name,
      mrid = String(mrid),
      changer_mrid = ch.mrid,
      kind = ch.kind,
      role = ch.role,
      side = ch.side,
      mode = br.tap_est_mode,
      alpha_deg = ch.alpha_deg,
      du = ch.du,
      neutral_step = ch.neutral,
      s_est = s_row[c],
      du_est = ch.du === nothing ? nothing : s_row[c] * ch.du,
      ratio_est = R,
      angle_est_deg = Phi,
      s_fixed = s_fix[c],
      step_model = ch.neutral + ts.s0[c],
      step_fixed = ch.neutral + s_fix[c],
      deviation_steps = s_fix[c] - ts.s0[c],
      sigma_s = sg,
      round_distance_sigma = (sg === nothing || !(sg > 0.0)) ? nothing : _tap_step_distance(ch, s_row[c]) / sg,
      step_from_target = nothing,
      out_of_range = oor[c],
      fixed = fixedRow,
      # where the impedance of the transformer comes from in this run
      # (_tap_reactance_source): the per-step data, the analytic correction
      # or none (the start impedance)
      reactance_source = !ts.impedance ? :nominal : _tap_set_has_table(ts.set) ? :table : :analytic,
      # :estimated (a state of the run), :out_of_range (estimated too far
      # outside the band, the transferred step stays, reason :none) or
      # :not_estimated with the reason (a guard freeze, :disabled,
      # :fallback_not_converged, :excluded, :machine_transformer);
      # _tap_machine_rows turns a machine row into :calculated after a run
      # with update_net, or into reason :q_not_measured_in_group
      status = estimated ? :estimated : outr ? :out_of_range : :not_estimated,
      reason = (estimated || outr) ? :none : ts.frozen[c],
      note = ts.notes[c],
      # the import carried no step, step_model is the middle position;
      # every later rewrite of the notice keeps the text (_tap_assumed_notice)
      step_assumed = ch.step_assumed,
      notice = join(filter(!isempty, [ch.step_assumed ? _tap_notice_text(:step_assumed) : "", oorText, _tap_freeze_text(ts, c), _tap_regulator_notice(relative ? :none : ts.notes[c], vertex, guard), relative ? _tap_relative_text(ts.rel_ref[c]) : "", (estimated && ts.weak[c] > 0.0) ? _tap_weak_text(ts.weak[c]) : "", ts.frozen[c] in _TAP_RANK_REASONS ? String(voltage_hint) : "", partner, String(extra_notice[c])]), "; "),
    ))
  end
  return rows
end

## Notice on an estimated changer of `ts` whose partner `k` on the same
## transformer the rank test or the recheck kept out or the band verdict
## did not take over (`:out_of_range`): the estimate holds only together
## with the transferred step of `k`
function _tap_partner_text(ts::TapBranchState, k::Int)::String
  ch = ts.set.changers[k]
  name = isempty(ch.name) ? string("regulator ", k) : string("changer ", ch.name)
  return string(name, " of this transformer is not determinable (", ts.frozen[k], "), this estimate depends on its transferred step ", _tap_notice_step(ch.neutral + ts.s0[k]))
end

## The voltage hint of a transformer (snet branch `k`) whose regulator the
## rank test kept out: measured minus
## estimated |V| of the active Vm rows at its two terminals (final
## residuals `r`, pu), the larger of the two. A wrong transferred step
## shows there when a measurement sees that bus. Empty when no terminal
## carries an active Vm row.
function _tap_voltage_hint(net::Net, k::Int, activeMeas::Vector{Measurement}, r::AbstractVector{Float64})::String
  br = net.branchVec[k]
  ends = (Int(br.fromBus), Int(br.toBus))
  best = 0
  for (i, m) in enumerate(activeMeas)
    (m.typ == VmMeas && m.busIdx !== nothing && Int(m.busIdx) in ends) || continue
    (best == 0 || abs(r[i]) > abs(r[best])) && (best = i)
  end
  best == 0 && return ""
  return string("voltage hint: measured minus estimated |V| at ", getCompName(net.nodeVec[Int(activeMeas[best].busIdx)].comp), " ", round(r[best]; sigdigits = 3), " pu")
end

## the report row of a released transformer without a changer of the
## released role (regulator 0)
function _tap_no_changer_row(br::Branch, ob::Int, mrid::AbstractString)::NamedTuple
  return (branch = ob, regulator = 0, name = getCompName(br.comp), changer = "", mrid = String(mrid), changer_mrid = "", kind = :none, role = :none, side = :none, mode = br.tap_est_mode, alpha_deg = nothing, du = nothing, neutral_step = 0, s_est = 0.0, du_est = nothing, ratio_est = br.tap_ratio, angle_est_deg = br.phase_shift_deg, s_fixed = 0.0, step_model = 0.0, step_fixed = 0.0, deviation_steps = 0.0, sigma_s = nothing, round_distance_sigma = nothing, step_from_target = nothing, out_of_range = false, fixed = false, reactance_source = :nominal, status = :not_estimated, reason = :no_changer, note = :none, step_assumed = false, notice = _tap_notice_text(:no_changer))
end

# the vertex share of a phase regulator at its estimate: |Re(f'/f)| /
# |f'/f|; below 0.1 only the magnitude effect pins it
function _tap_magnitude_share(ch::SETapChanger, s::Float64)::Float64
  re, ab = _tap_factor_log_slope(ch, s)
  return ab == 0.0 ? 0.0 : abs(re) / ab
end

const _TAP_NOTICE_TEXTS = Dict{Symbol,String}(
  # the import carried no step, the changer started at its middle position
  # (SETapChanger.step_assumed)
  :step_assumed => "step not transferred, middle position assumed",
  :radial_no_voltage_pin => "not estimated (rank test): bridge transformer without a voltage measurement on its cut-off side, the tap would absorb the downstream voltage",
  :no_angle_anchor => "not estimated (rank test): no angle measurement (Va, Ia) or meshed path ties its far side to the reference, and the ratio regulator of the same transformer covers its magnitude effect",
  :symmetric_start => "not estimated (rank test): pure angle regulator without an angle anchor (no Va or Ia row, no meshed path), its step is not determinable",
  :not_separable => "not estimated (rank test): its effect on the branch ratio equals that of another regulator of the same transformer",
  :not_observable => "not estimated (rank test): the other columns of the Jacobian explain its column at the start point, no structural cause found",
  :zero_column => "not estimated: the measurements do not see this regulator (Jacobian column zero within its finite-difference noise)",
  :not_observable_at_estimate => "not estimated (rank test at the estimate): at the converged state with the fixed steps the other columns of the Jacobian explain its column, the estimate is withdrawn and the transferred step stays",
  :weak_at_estimate => "not estimated (rank test at the estimate): determinable only at the converged state, and estimated again there with sigma_s not below 0.5 steps, the estimate is not taken over and the transferred step stays",
  :not_converged_at_estimate => "not estimated (recheck at the estimate): estimated again after a changer of the same transformer was withdrawn to its transferred step, that run did not converge, the transferred step stays",
  :rank_budget => "not estimated: the number of released regulators puts the rank test above its work budget (one minute on a large net), no regulator of the run is estimated, the transferred step stays",
  :mirror_step_fits => "not estimated: no angle anchor ties its far side, and its mirror step fits the measurements as well as the fixed step (J difference below 9), the estimate is withdrawn and the transferred step stays",
  :model_mismatch => "frozen: the branch ratio and angle do not lie on the curves of its tap changer data",
  :disabled => "not estimated: tap estimation is switched off (state_estimation.tap_estimation.enabled = false), the transferred step stays",
  :excluded => "not estimated: named in state_estimation.tap_estimation.exclude, the transferred step stays and is not calculated",
  :machine_transformer => "not estimated: machine (generator step-up) transformer (state_estimation.tap_estimation.include_machine_transformers = false), its step is calculated by a converged run with update_net = true",
  :q_not_measured_in_group => "not calculated: machine transformer whose machine node has no Q measurement on every bus, the transferred step stays",
  :fallback_not_converged => "not estimated: the run with tap estimation did not converge, the result is the run without tap estimation, the transferred step stays",
  :fallback_worse_fit => "not estimated: the fixed steps of the run with tap estimation fitted the measurements worse than the transferred steps (worse_fit), the result is the run without tap estimation, the transferred step stays",
  :no_changer => "not released: the transformer has no tap changer of the released kind",
  :relative_only => "parallel member: estimated only relative to the other member",
  :relative_to_frozen => "not fixed: estimated only relative to the frozen phase regulator of",
  :conditional => "ratio step holds only together with the frozen phase position",
  :sign_from_start => "only the magnitude effect pins this regulator, the sign of its step comes from the start position",
  :vertex => "only the magnitude effect pins this regulator, its step is not reliable",
  :reimport_guard => "a re-import would ignore this changer at this step (its factor leaves 0.5 to 2)",
  # appended by runse! after the topology precheck, no reason and no note
  :mesh_without_p_measurement => "lies in a mesh without P measurement: the step may not be unique",
)

## notice text of a freeze reason or note
_tap_notice_text(key::Symbol)::String = get(_TAP_NOTICE_TEXTS, key, "")

# `text` as the new notice of row `r`, after the step-assumed text when the
# import carried no step for the changer of the row (field step_assumed)
_tap_assumed_notice(r::NamedTuple, text::AbstractString)::String = get(r, :step_assumed, false) ? join(filter(!isempty, [_tap_notice_text(:step_assumed), String(text)]), "; ") : String(text)

## The freeze text of changer `c` of `ts`, true for its kind: a column the noise test froze says that no
## measurement sees it; `:symmetric_start` and `:no_angle_anchor` name what
## the changer moves. A pure angle changer (symmetrical, linear, a table of
## constant ratio) moves no magnitude anywhere on its band; an asymmetrical
## changer and a table of varying ratio do, only not at the start step (a
## turning point of |t|, a flat piece or a vertex of the table ratio).
function _tap_freeze_text(ts::TapBranchState, c::Int)::String
  r = ts.frozen[c]
  r === :none && return ""
  ch = ts.set.changers[c]
  r === :not_observable && ts.cause[c] === :zero_column && return _tap_notice_text(:zero_column)
  pure = !_tap_changes_magnitude(ch)
  if r === :symmetric_start && !pure
    what = ch.kind === :tabular ? "the table ratio is flat or turns at the start step" : string("|t| turns at the start step (alpha ", round(something(ch.alpha_deg, 0.0); digits = 2), " deg)")
    return string("not estimated (rank test): ", what, " and no angle anchor (no Va or Ia row, no meshed path) ties its far side; the changer moves the magnitude away from this step, but from this start its step is not determinable")
  elseif r === :no_angle_anchor && pure
    return "not estimated (rank test): pure angle regulator, no angle measurement (Va, Ia) or meshed path ties its far side to the reference"
  end
  return _tap_notice_text(r)
end

# a step for the notices: whole steps as integers, others with two decimals
function _tap_notice_step(x::Float64)::String
  r = round(x; digits = 2)
  return isinteger(r) ? string(Int(r)) : string(r)
end

## Notice of an estimate `s` (relative to neutral) outside the band of
## changer `ch` (row flag out_of_range). For an
## asymmetrical changer with alpha not 0 the text also names the mirror
## step 2 s* - s across the turning point s* of |t|
## (`_tap_asym_turning_point`): |t| is the same there, and a changer pinned
## only by the magnitude of its ratio can end on either side. Steps are
## model steps (neutral + s). No ';' in the text: the warning line of
## runse! and run.log is built from the joined notice.
## `verdict` (`_tap_band_verdict`, with the distance `dist` in sigma_s)
## adds what the fixation made of it: `:band_end` the end step, `:out_of_range`
## the transferred step `s0`; `:inside` adds nothing (a row that was not
## fixed, or was withdrawn later).
function _tap_out_of_range_text(ch::SETapChanger, s::Float64; verdict::Symbol = :inside, dist::Union{Nothing,Float64} = nothing, s0::Float64 = 0.0)::String
  smin, smax = _tap_band(ch)
  txt = string("estimate at step ", _tap_notice_step(ch.neutral + s), " lies outside the band ", _tap_notice_step(ch.neutral + smin), " to ", _tap_notice_step(ch.neutral + smax))
  if ch.kind === :asymmetrical && ch.du !== nothing && ch.du != 0.0 && abs(sind(ch.alpha_deg)) > 1e-9
    sstar = _tap_asym_turning_point(ch)
    txt = string(txt, ", |t| has the same value at the mirror step ", _tap_notice_step(ch.neutral + 2.0 * sstar - s), " across its turning point at step ", _tap_notice_step(ch.neutral + sstar))
  end
  verdict === :inside && return txt
  endStep = _tap_notice_step(ch.neutral + (s > smax ? smax : smin))
  distText = dist === nothing ? "sigma_s not available" : string(@sprintf("%.3g", dist), " sigma_s beyond the band end step ", endStep)
  verdict === :band_end && return string(txt, ", ", distText, ", fixed to the end step ", endStep)
  return string(txt, ", ", distText, ", not taken over, the transferred step ", _tap_notice_step(ch.neutral + s0), " stays")
end

# distance (in sigma_s of the step) of an estimate outside the band from
# the band end up to which the fixation takes the end step; farther out
# the estimate is not taken over (status :out_of_range)
const _TAP_BAND_SIGMA_MAX = 2.0

## Band verdict of the continuous estimate `s` (relative to neutral) of
## changer `ch` with the step sigma `sg` (`nothing` = none from the
## covariance): `:inside` when it lies at most half a step outside the band
## (its nearest step is in the band, row flag out_of_range false),
## `:band_end` when it lies farther out but at most `_TAP_BAND_SIGMA_MAX`
## sigma_s beyond the band end step (the fixation takes that step), else
## `:out_of_range` (the transferred step stays). Without a sigma_s the
## distance is not known and the verdict is `:out_of_range`. Returns
## `(verdict, distance in sigma_s or nothing)`.
function _tap_band_verdict(ch::SETapChanger, s::Float64, sg)::Tuple{Symbol,Union{Nothing,Float64}}
  smin, smax = _tap_band(ch)
  (smin - 0.5 <= s <= smax + 0.5) && return (:inside, nothing)
  d = s > smax ? s - smax : smin - s
  ds = (sg === nothing || !(sg > 0.0)) ? nothing : d / sg
  return ((ds !== nothing && ds <= _TAP_BAND_SIGMA_MAX) ? :band_end : :out_of_range, ds)
end

## The changers of the estimation map `map` at the continuous estimate `x`
## (sigma per column in `sigma`) among the columns `cols` whose estimate
## lies too far outside the band (`_tap_band_verdict` `:out_of_range`): each
## is frozen with `:out_of_range`, so the fixation (`_tap_fix_columns`) puts
## it at its transferred step instead of the band end. Returns
## `(state, changer, note, tree_phase)` of each, with the note and flag the
## freeze cleared, so a caller that puts its state back can undo it.
function _tap_mark_out_of_range!(map::TapStateMap, x::Vector{Float64}, sigma::AbstractDict, cols::AbstractSet{Int})::Vector{Tuple{TapBranchState,Int,Symbol,Bool}}
  out = Tuple{TapBranchState,Int,Symbol,Bool}[]
  for ts in map.trafos, c in eachindex(ts.cols)
    col = ts.cols[c]
    (col != 0 && col in cols && ts.frozen[c] === :none && _tap_fixable(ts, c)) || continue
    first(_tap_band_verdict(ts.set.changers[c], x[col], get(sigma, col, nothing))) === :out_of_range || continue
    push!(out, (ts, c, ts.notes[c], ts.tree_phase[c]))
    _tap_freeze!(ts, c, :out_of_range)
  end
  return out
end

## The warning line of a tap row whose estimate lies outside the band
## (`t.out_of_range`): regulator name and its notice. runse! logs it once
## per such row, the SE service writes it into run.log.
_tap_out_of_range_line(t::NamedTuple)::String = string("tap regulator ", t.name, isempty(t.changer) ? "" : "/" * t.changer, " estimated outside its band: ", t.notice)

## notice text of note :relative_to_frozen with the names of the frozen
## phase regulators (`TapBranchState.rel_ref`)
_tap_relative_text(names::AbstractString)::String = string(_tap_notice_text(:relative_to_frozen), " ", names, ", its continuous value stays in the fit")

## the notice of one regulator: its freeze reason, its note and the notices
## of the estimate (vertex, re-import guard), joined
function _tap_regulator_notice(keys::Symbol...)::String
  parts = String[_tap_notice_text(k) for k in keys if k !== :none]
  return join(filter(!isempty, parts), "; ")
end

## text of the weak-column hint (no effect on the freeze; factor from
## `_tap_weak_hints!`)
_tap_weak_text(factor::Float64)::String = "weak column, factor $(round(factor; sigdigits = 3)) below the strongest"

## Write the changer positions `s_new` onto branch `br` (the write-back of
## update_taps, the diagnostics stamp, the island mirror's source): ratio and
## shift from the set; with an impedance source of the transformer under
## the mode `impedance` (`_tap_reactance_source` not `:nominal`, as in
## _build_tap_state_map) the impedance by the factor ratio
## F(s_new) / F(s_from) on the live values, r_base_pu/x_base_pu included,
## else r and x stay; a stored set gets the new steps. A typed set
## moves the model steps instead and resolves the winding; the resolver
## always rewrites r_pu and x_pu from the model, so without an applied
## impedance source the values of the run are put back, and under `off`
## r of the run (the resolver's analytic factor is not taken). With
## `continuous = true` a typed set is written like a stored one (ratio,
## shift and impedance from the set at `s_new`, the model steps stay), so a
## position between two steps is not rounded: the diagnostics stamp of a
## regulator that stayed a state (note `:relative_to_frozen`).
function _tap_apply_position!(br::Branch, set::SETapChangerSet, s_new::Vector{Float64}; s_from::Vector{Float64}, impedance::Symbol, continuous::Bool = false)
  applied = _tap_reactance_source(set, impedance) !== :nominal
  if !continuous && set.origin === :typed && br.tap_winding !== nothing
    w = br.tap_winding
    c = 1
    # a model step written off its assumed middle position is a step now
    if w.taps !== nothing && w.taps.tapStepPercent != 0.0
      w.taps.step = round(Int, set.changers[c].neutral + s_new[c])
      s_new[c] == s_from[c] || (br.ratio_step_assumed = false)
      c += 1
    end
    if w.phase_taps !== nothing && c <= length(set.changers)
      w.phase_taps.step = round(Int, set.changers[c].neutral + s_new[c])
      s_new[c] == s_from[c] || (br.phase_step_assumed = false)
    end
    r_run, x_run = br.r_pu, br.x_pu
    resolve_branch_taps!(br, w; tap_changer_model = br.tap_correction)
    applied || (br.r_pu = r_run; br.x_pu = x_run)
    # under off the resolver's analytic factor on r is not taken (the
    # estimation kept r), the x of the curve or table is
    impedance === :off && (br.r_pu = r_run)
    return br
  end
  R, Phi = _tap_set_ratio_shift(set, s_new)
  br.tap_ratio = R
  br.phase_shift_deg = Phi
  # the same for a set of the branch fields: its ratio changer is the
  # ratio changer of the fields, every other one the phase changer
  if set.origin === :fields
    for c in eachindex(set.changers)
      s_new[c] == s_from[c] && continue
      set.changers[c].role === :ratio ? (br.ratio_step_assumed = false) : (br.phase_step_assumed = false)
    end
  end
  if applied
    tsi = br.tap_step_impedance
    fs = _tap_set_under(set, impedance)
    Ff = _tap_impedance_factors(fs, s_from; s0 = s_from, tsi = tsi)
    Fn = _tap_impedance_factors(fs, s_new; s0 = s_from, tsi = tsi)
    br.r_pu *= Fn[1] / Ff[1]
    br.x_pu *= Fn[2] / Ff[2]
    br.r_base_pu = br.r_pu
    br.x_base_pu = br.x_pu
    set_branch_shunt!(br; g_from_pu = br.g_from_pu * (Fn[3] / Ff[3]), b_from_pu = br.b_from_pu * (Fn[4] / Ff[4]), g_to_pu = br.g_to_pu * (Fn[3] / Ff[3]), b_to_pu = br.b_to_pu * (Fn[4] / Ff[4]))
  end
  if set.origin === :cgmes || set.origin === :powsybl
    chs = SETapChanger[_tap_changer_at(set.changers[c], set.changers[c].neutral + s_new[c]) for c in eachindex(set.changers)]
    br.se_tap_changers = SETapChangerSet(set.ratio_n, set.shift_n_deg, chs, set.origin)
  end
  return br
end

# the same changer at another import position; a step written back off the
# assumed middle position is no longer an assumption
_tap_changer_at(ch::SETapChanger, step::Float64) = SETapChanger(ch.name, ch.mrid, ch.kind, ch.role, ch.side, ch.additive, ch.du, ch.alpha_deg, ch.dphi_deg, ch.step_abs, ch.neutral, ch.low, ch.high, step, ch.table, ch.outside_du, ch.x_curve, ch.x_curve_form, ch.impedance, ch.rx_correction, ch.step_assumed && step == ch.step)

## Write the end steps of an SE result onto the branches of `net` (the
## diagnostics of `validate_measurements` and the bad-data rounds of the SE
## service evaluate the residuals there). Per regulator: an estimated one at
## its fixed step, a frozen one at its start position (s_fixed of its row),
## the other changers of the set at their start, one that stayed a state
## (note `:relative_to_frozen`) at its continuous value. A machine
## transformer row (status `:calculated`) stays at its start position: its
## s_fixed is the step calculated after the run, which the run never solved
## with (`_tap_machine_rows`). A result without a fixed tap run (no tap
## rows, not fixed, the fallback) leaves the model.
function _tap_stamp_end_steps!(net::Net, result, impedance::Symbol)
  (result.tapEstimates !== nothing && result.tapFixation !== nothing && result.tapFixation.fixed) || return net
  for b in unique(tr.branch for tr in result.tapEstimates)
    (1 <= b <= length(net.branchVec)) || continue
    br = net.branchVec[b]
    set = _effective_tap_changers(br)
    set === nothing && continue
    sstart = _tap_s_from_branch(set, br)
    sstart === nothing && continue
    s_w = copy(sstart)
    cont = false
    for tr in result.tapEstimates
      (tr.branch == b && 1 <= tr.regulator <= length(s_w) && tr.status !== :calculated) || continue
      s_w[tr.regulator] = tr.s_fixed
      cont |= tr.note === :relative_to_frozen
    end
    s_w == sstart || _tap_apply_position!(br, set, s_w; s_from = sstart, impedance = impedance, continuous = cont)
  end
  return net
end

## Copy the tap position (and with it the impedance a write-back changed)
## from the island copy `src` onto the parent branch `dst`: ratio and
## shift, r/x with their equipment base, the shunt arms, the stored set,
## the step-assumed marks, and for a typed winding the model steps with the resolved neutral point.
function _copy_tap_position!(dst::Branch, src::Branch)
  dst.tap_ratio = src.tap_ratio
  dst.phase_shift_deg = src.phase_shift_deg
  dst.r_pu = src.r_pu
  dst.x_pu = src.x_pu
  dst.r_base_pu = src.r_base_pu
  dst.x_base_pu = src.x_base_pu
  set_branch_shunt!(dst; g_from_pu = src.g_from_pu, b_from_pu = src.b_from_pu, g_to_pu = src.g_to_pu, b_to_pu = src.b_to_pu)
  dst.se_tap_changers = src.se_tap_changers
  # the marks of a step the import did not carry travel with the step
  dst.ratio_step_assumed = src.ratio_step_assumed
  dst.phase_step_assumed = src.phase_step_assumed
  if src.taps_derived && src.tap_winding !== nothing && dst.tap_winding !== nothing
    src.tap_winding.taps === nothing || dst.tap_winding.taps === nothing || (dst.tap_winding.taps.step = src.tap_winding.taps.step)
    src.tap_winding.phase_taps === nothing || dst.tap_winding.phase_taps === nothing || (dst.tap_winding.phase_taps.step = src.tap_winding.phase_taps.step)
    dst.ratio = src.ratio
    dst.angle = src.angle
  end
  return dst
end

## The model data an SE run may write back (update_taps, update_shunts),
## taken before run 1 so the fallback run 2 starts from the same model: an
## island that converged in run 1 writes its steps and shunts into the
## parent net even when another island fails and the whole net falls back.
## Empty when the run writes nothing back.
function _se_writeback_snapshot(net::Net, cfg::StateEstimationConfig)
  taps = cfg.update_taps ? [deepcopy(br) for br in net.branchVec] : Branch[]
  shunts = cfg.update_shunts ? [(sh.y_pu_shunt, sh.B_shunt) for sh in net.shuntVec] : Tuple{ComplexF64,Float64}[]
  return (taps = taps, shunts = shunts)
end

## Put the snapshot of `_se_writeback_snapshot` back onto `net`.
function _se_writeback_restore!(net::Net, snap)
  for (k, br) in enumerate(snap.taps)
    _copy_tap_position!(net.branchVec[k], br)
  end
  for (k, (y, b)) in enumerate(snap.shunts)
    net.shuntVec[k].y_pu_shunt = y
    net.shuntVec[k].B_shunt = b
  end
  return net
end

## 95 % quantile of the chi-square distribution with k degrees of freedom
## (the Delta-J notice): bisection on the regularized lower
## incomplete gamma P(k/2, x/2) by its series; k is a small count of tap
## states
function _chi2_quantile95(k::Int)::Float64
  k >= 1 || return 0.0
  a = k / 2
  # log Gamma(a) for a whole or half-whole a
  lg = isinteger(a) ? sum(log, 1:(Int(a) - 1); init = 0.0) : (0.5 * log(pi) + sum(j -> log(j - 0.5), 1:Int(a - 0.5); init = 0.0))
  function cdf(x)
    y = x / 2
    y <= 0.0 && return 0.0
    term = 1.0 / a
    sum_ = term
    for n in 1:1000
      term *= y / (a + n)
      sum_ += term
      term < 1e-16 * sum_ && break
    end
    return exp(a * log(y) - y - lg) * sum_
  end
  lo, hi = 0.0, max(10.0, 4.0 * k + 20.0)
  for _ in 1:200
    mid = (lo + hi) / 2
    cdf(mid) < 0.95 ? (lo = mid) : (hi = mid)
  end
  return (lo + hi) / 2
end

"""
    setTapEstimation!(net; trafo, mode = :ratio, alpha_deg = nothing, enabled = true)

Release the tap changers of a transformer as state-estimation states
(0.10.0; one step state per changer since 0.60.0). The changers come from
the branch's changer set: what an import stored, the typed models of a
resolved winding, or the nameplate fields. `mode` selects them by role:
`:ratio` the ratio changers, `:pst` the phase changers, `:both` all; `:pst`
and `:both` need a declared phase changer, `:ratio` a ratio changer
(`ArgumentError` otherwise). The direction alpha of an asymmetrical phase
changer comes from the import and is never estimated: `alpha_deg` is
optional and, when given, has to equal it (`ArgumentError` naming both, or
when the selection has no asymmetrical changer). `enabled = false` clears
the release. Neither call writes the stored direction `tap_est_alpha_deg`.

A branch that does not lie on the curves of its changer data is no error
here: the release is recorded and the estimator freezes the regulators with
`:model_mismatch` and a notice. The release does not override the
configuration: a transformer named in `state_estimation.tap_estimation.exclude`
stays out (`:excluded`), and so does a machine transformer unless
`include_machine_transformers` is set (`:machine_transformer`, calculated
after the run).

Returns `(branch, mode, alpha_deg, s0)`: `alpha_deg` of the first selected
asymmetrical changer or `nothing`, `s0` the start positions (steps relative
to neutral) of the selected changers in set order, `nothing` on a model
mismatch, empty for `enabled = false`.
"""
function setTapEstimation!(net::Net; trafo, mode::Symbol = :ratio, alpha_deg::Union{Nothing,Real} = nothing, enabled::Bool = true)
  k = _resolve_trafo_branch(net, trafo)
  br = net.branchVec[k]
  if !enabled
    br.tap_est_mode = :none
    return (branch = k, mode = :none, alpha_deg = nothing, s0 = Float64[])
  end
  mode in (:ratio, :pst, :both) || throw(ArgumentError("setTapEstimation!: mode must be :ratio, :pst, or :both (got :$(mode))"))
  set = _effective_tap_changers(br)
  (set === nothing || isempty(set.changers)) && throw(ArgumentError("setTapEstimation!: branch $(k) carries no tap changer data"))
  hasR = any(ch -> ch.role === :ratio, set.changers)
  hasP = any(ch -> ch.role === :phase, set.changers)
  mode === :ratio && !hasR && throw(ArgumentError("setTapEstimation!: branch $(k) has no ratio tap changer, mode :ratio has nothing to release"))
  mode !== :ratio && !hasP && throw(ArgumentError("setTapEstimation!: branch $(k) has no declared phase tap changer, mode :$(mode) needs one (nameplate phase data, a typed phase model or a stored changer set)"))
  sel = [c for c in eachindex(set.changers) if _tap_role_released(set.changers[c].role, mode)]
  asym = [c for c in sel if set.changers[c].kind === :asymmetrical]
  if alpha_deg !== nothing
    isempty(asym) && throw(ArgumentError("setTapEstimation!: alpha_deg applies to asymmetrical tap changers only; branch $(k) mode :$(mode) selects none"))
    for c in asym
      a = set.changers[c].alpha_deg
      abs(Float64(alpha_deg) - a) <= 1e-9 || throw(ArgumentError("setTapEstimation!: alpha_deg = $(alpha_deg) differs from the direction $(a) degrees the import gives the tap changer $(set.changers[c].name) of branch $(k); alpha is never estimated"))
    end
  end
  br.tap_est_mode = mode
  s0 = _tap_s_from_branch(set, br)
  return (branch = k, mode = mode, alpha_deg = isempty(asym) ? nothing : set.changers[asym[1]].alpha_deg, s0 = s0 === nothing ? nothing : s0[sel])
end

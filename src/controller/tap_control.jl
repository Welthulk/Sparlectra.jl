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

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 1.3.2026
# file: src/controller/tap_control.jl
# purpose: transformer tap controllers for the outer control loop
#          (addPowerTransformerControl!, run_tap_controllers_outer!) with
#          ratio and phase tap probing plus controller report output

@inline _voltage_within_deadband(vm::Float64, target_vm_pu::Float64, deadband_vm_pu::Float64)::Bool = abs(vm - target_vm_pu) <= deadband_vm_pu
@inline _voltage_control_error(vm::Float64, target_vm_pu::Float64, metric::Symbol)::Float64 = metric == :vm2 ? vm^2 - target_vm_pu^2 : vm - target_vm_pu

# Cross-type check of a voltage controller (machine, shunt, HVDC terminal)
# against a tap controller on the same bus. The two fight each other in the
# outer loop unless they stand at different ranks (control_priority) and
# agree on the target within both deadbands: then the lower rank waits
# while the higher one holds the bus in its band and takes over at its
# limit. With different targets the lower rank pulls the voltage out of the
# higher one's band whenever that one is in band, and both walk against
# each other until one sits at its limit, at any ranks. A target on the
# edge of the band agrees: 1.01 against 1.0 +/- 0.01 differs by
# 0.010000000000000009 in floats, so the comparison allows 1e-9.
function _fights_tap_controller(tap::PowerTransformerControl, priority::Int, target_vm_pu::Float64, deadband_vm_pu::Float64)::Bool
  control_priority(tap) == priority && return true
  return abs(tap.target_vm_pu - target_vm_pu) > min(tap.deadband_vm_pu, deadband_vm_pu) + 1e-9
end

# Ratio window the controller moves its tap in: its own control range
# (tap_min/tap_max of the controller) where set, the physical range of the
# branch otherwise.
_tap_control_range(ctrl::PowerTransformerControl, br::Branch)::Tuple{Float64,Float64} = (something(ctrl.tap_min, br.tap_min), something(ctrl.tap_max, br.tap_max))

# Check a control range against the transformer and store it on the
# controller. Both ends lie inside the physical range of the branch
# (tap_min/tap_max, 1e-12 slack for ratios read back from a file) and the
# lower end below the upper one; a range on a controller without ratio
# control is refused. `nothing` keeps the physical limit on that side.
function _set_tap_control_range!(ctrl::PowerTransformerControl, br::Branch, tap_min::Union{Nothing,Real}, tap_max::Union{Nothing,Real})
  (tap_min === nothing && tap_max === nothing) && return ctrl
  who = string("PowerTransformerControl ", br.comp.cName)
  ctrl.control_ratio || throw(ArgumentError(string(who, ": a control range (tap_min/tap_max) needs ratio control (control_ratio = true)")))
  lo = tap_min === nothing ? br.tap_min : Float64(tap_min)
  hi = tap_max === nothing ? br.tap_max : Float64(tap_max)
  for (key, v) in (("tap_min", tap_min), ("tap_max", tap_max))
    v === nothing && continue
    (isfinite(v) && br.tap_min - 1e-12 <= v <= br.tap_max + 1e-12) || throw(ArgumentError(string(who, ": control range ", key, " = ", v, " lies outside the physical tap range ", br.tap_min, " .. ", br.tap_max)))
  end
  lo < hi || throw(ArgumentError(string(who, ": control range tap_min = ", lo, " must be below tap_max = ", hi)))
  ctrl.tap_min = tap_min === nothing ? nothing : lo
  ctrl.tap_max = tap_max === nothing ? nothing : hi
  return ctrl
end

"""
    _tap_controllers(net)

Collect all tap controllers configured in `net` from transformer windings.
"""
function _tap_controllers(net::Net)::Vector{PowerTransformerControl}
  controllers = PowerTransformerControl[]
  seen = IdDict{PowerTransformerControl,Bool}()

  for trafo in net.trafos
    for winding in (trafo.side1, trafo.side2, trafo.side3)
      isnothing(winding) && continue
      for ctrl in winding.controls
        haskey(seen, ctrl) && continue
        seen[ctrl] = true
        push!(controllers, ctrl)
      end
    end
  end
  return controllers
end

"""
    clearTapControllers!(net)

Remove all tap controllers from transformer windings in `net`.
"""
function clearTapControllers!(net::Net)
  for trafo in net.trafos
    empty!(trafo.side1.controls)
    empty!(trafo.side2.controls)
    if !isnothing(trafo.side3)
      empty!(trafo.side3.controls)
    end
  end
  return net
end

"""
    _find_trafo_branch(net, name)

Internal helper to resolve a transformer branch by component name, component id,
or branch index represented as string.
"""
function _find_trafo_branch(net::Net, name::String)::Branch
  for br in net.branchVec
    if br.ratio != 0.0 && (br.comp.cName == name || br.comp.cID == name || string(br.branchIdx) == name)
      return br
    end
  end
  error("PowerTransformerControl: transformer $(name) not found. Use branch component name/id or branch index as String.")
end

"""
    get_bus_vm_pu(net, bus_name)

Return solved voltage magnitude in p.u. for bus `bus_name`.
"""
function get_bus_vm_pu(net::Net, bus_name::String)
  bus = geNetBusIdx(net = net, busName = bus_name)
  return net.nodeVec[bus]._vm_pu
end

"""
    get_branch_p_from_to_mw(net, from_bus, to_bus)

Return active power flow in MW for the oriented branch direction
`from_bus -> to_bus`.
"""
function get_branch_p_from_to_mw(net::Net, from_bus::String, to_bus::String)
  br = getNetBranch(net = net, fromBus = from_bus, toBus = to_bus)
  from = geNetBusIdx(net = net, busName = from_bus)
  return br.fromBus == from ? br.fBranchFlow.pFlow : br.tBranchFlow.pFlow
end

"""
    get_branch_q_from_to_mvar(net, from_bus, to_bus)

Return reactive power flow in MVAr for the oriented branch direction
`from_bus -> to_bus`.
"""
function get_branch_q_from_to_mvar(net::Net, from_bus::String, to_bus::String)
  br = getNetBranch(net = net, fromBus = from_bus, toBus = to_bus)
  from = geNetBusIdx(net = net, busName = from_bus)
  return br.fromBus == from ? br.fBranchFlow.qFlow : br.tBranchFlow.qFlow
end

@inline _controller_mode_label(mode::Symbol) = String(mode)

# the at-limit flag of a tap controller: ratio or shift at an end of its
# control range, a model step or table point at the edge of the window
function _tap_at_range_end(ctrl::PowerTransformerControl, br)::Bool
  lo, hi = _tap_control_range(ctrl, br)
  (isapprox(br.tap_ratio, lo; atol = 1e-12) || isapprox(br.tap_ratio, hi; atol = 1e-12)) && return true
  (_model_step_at_window_edge(ctrl, br) || _table_point_at_window_edge(ctrl, br)) && return true
  return isapprox(br.phase_shift_deg, br.phase_min_deg; atol = 1e-12) || isapprox(br.phase_shift_deg, br.phase_max_deg; atol = 1e-12)
end
function control_name(ctrl::PowerTransformerControl)
  base = strip(ctrl.trafo)
  mode = ctrl.mode == :voltage ? "OLTC" : (ctrl.mode == :branch_active_power ? "PST" : "SCHRAEG")
  return isempty(base) ? mode : string(base, " ", mode)
end
control_enabled(ctrl::PowerTransformerControl) = ctrl.enabled
control_initialize!(::PowerTransformerControl, ::Net, context) = NoControlState()
control_status(ctrl::PowerTransformerControl, ::AbstractControlState)::Symbol = ctrl.status
control_is_converged(ctrl::PowerTransformerControl, ::AbstractControlState)::Bool = ctrl.converged
control_is_blocked(ctrl::PowerTransformerControl, ::AbstractControlState)::Bool = ctrl.at_limit
function control_evaluate!(ctrl::PowerTransformerControl, net::Net, ::AbstractControlState, context)
  br = _find_trafo_branch(net, ctrl.trafo)
  ctrl.at_limit = false
  converged_v = true
  converged_p = true
  if ctrl.mode in (:voltage, :voltage_and_branch_active_power)
    vm = get_bus_vm_pu(net, ctrl.target_bus)
    ctrl.achieved_vm_pu = vm
    converged_v = _voltage_within_deadband(vm, ctrl.target_vm_pu, ctrl.deadband_vm_pu)
  end
  if ctrl.mode in (:branch_active_power, :voltage_and_branch_active_power)
    p = get_branch_p_from_to_mw(net, ctrl.target_branch[1], ctrl.target_branch[2])
    ctrl.achieved_p_mw = p
    converged_p = abs(p - ctrl.p_target_mw) <= ctrl.deadband_p_mw
  end
  ctrl.converged = converged_v && converged_p
  # a ratio outside the control range (the case ratio left out by a
  # narrowed range) holds the controller open until it stands inside
  ctrl.converged = ctrl.converged && !_ratio_outside_control_range(ctrl, br)
  # a locked controller (QSTS dead time) proposes nothing, so apply never
  # sets its flag; the rank gate of run_control! reads it from here
  ctrl.locked && (ctrl.at_limit = _tap_at_range_end(ctrl, br))
  ctrl.status = ctrl.converged ? :converged : :active
  ctrl.outer_iters = context.outer_iteration
  return nothing
end
function control_propose_update!(ctrl::PowerTransformerControl, net::Net, ::AbstractControlState, context)
  br = _find_trafo_branch(net, ctrl.trafo)
  old_ratio = br.tap_ratio
  old_phase = br.phase_shift_deg
  new_ratio = old_ratio
  new_phase = old_phase
  # The probe solves run under the run's own power-flow configuration. Its
  # start strategy matters because imported networks (CGMES/MATPOWER) with
  # off-nominal ratio branches routinely fail from a flat start while
  # converging from their stored profile; its island switches matter because
  # a probe with islands off refused every net with two AC islands that both
  # have branches, although the run itself solves them island by island.
  probe_cfg = context.pf_config
  if _ratio_outside_control_range(ctrl, br)
    # the range wins over the setpoint: the ratio goes onto the nearer end
    # of its range (on a tabular changer the table point nearest it)
    # whatever the voltage error, without a probe, and the phase waits
    # for the next pass
    lo, hi = _tap_control_range(ctrl, br)
    _warn_start_outside_control_range(ctrl, br, lo, hi)
    return (old_ratio = old_ratio, new_ratio = _ratio_toward_control_range(br, lo, hi), old_phase = old_phase, new_phase = old_phase)
  end
  if !ctrl.converged && ctrl.control_ratio && ctrl.mode in (:voltage, :voltage_and_branch_active_power)
    direction = _ratio_probe_direction(net, br, ctrl, probe_cfg)
    direction == 0.0 && (direction = -1.0)
    Δ = ctrl.is_discrete ? br.tap_step : 0.25 * br.tap_step
    e_v = _voltage_control_error(ctrl.achieved_vm_pu, ctrl.target_vm_pu, ctrl.voltage_error_metric)
    up = (e_v < 0.0) ? direction : -direction
    lo, hi = _tap_control_range(ctrl, br)
    # a tabular changer steps onto its next table ratio inside the control
    # range; tap_step is only the mean step there and would stop between
    # two points
    new_ratio = (ctrl.is_discrete && br.tap_ratio_points !== nothing) ? _ratio_points_shift(br, old_ratio, up > 0.0 ? 1 : -1; window = (lo, hi)) : clamp(old_ratio + up * Δ, lo, hi)
  end
  if !ctrl.converged && ctrl.control_phase && ctrl.mode in (:branch_active_power, :voltage_and_branch_active_power)
    direction = _phase_probe_direction(net, br, ctrl, probe_cfg)
    direction == 0.0 && (direction = -1.0)
    Δ = ctrl.is_discrete ? br.phase_step_deg : 0.25 * br.phase_step_deg
    e_p = ctrl.achieved_p_mw - ctrl.p_target_mw
    new_phase = clamp(old_phase + ((e_p < 0.0) ? direction * Δ : -direction * Δ), br.phase_min_deg, br.phase_max_deg)
  end
  return (old_ratio = old_ratio, new_ratio = new_ratio, old_phase = old_phase, new_phase = new_phase)
end
# Typed phase-tap model of the winding that carries this controller:
# `nothing` for controllers on windings without one (MATPOWER general-case
# PSTs, CGMES-flattened PSTs), which keep static reactance.
function _phase_tap_model_of(net::Net, ctrl::PowerTransformerControl)::Union{Nothing,PhaseTapChangerModel}
  for tf in net.trafos
    for w in (tf.side1, tf.side2, tf.side3)
      w === nothing && continue
      any(c -> c === ctrl, w.controls) && return w.phase_taps
    end
  end
  return nothing
end

# X(α) of the controlled PST at a continuous controller angle, in pu.
# Formula models (symmetrical/asymmetrical with x_min/x_max) evaluate
# calcPhaseTapReactance directly at the continuous angle; tabular models snap
# to the nearest table row by angle (no interpolation, the #261 Stage-4
# decision) and use its per-step x_pu. Returns `nothing` when the model
# carries no reactance data, which keeps the coupling strictly opt-in.
function _phase_tap_reactance_at(net::Net, ctrl::PowerTransformerControl, angle_deg::Float64)::Union{Nothing,Float64}
  model = _phase_tap_model_of(net, ctrl)
  model === nothing && return nothing
  if model.kind === :tabular
    model.table === nothing && return nothing
    best = nothing
    bestd = Inf
    for p in model.table
      d = abs(p.angle_deg - angle_deg)
      if d < bestd
        best = p
        bestd = d
      end
    end
    return best === nothing ? nothing : best.x_pu
  end
  (model.x_min === nothing || model.x_max === nothing) && return nothing
  return calcPhaseTapReactance(model, angle_deg)
end

# The typed model of a derived branch for one actuator (`:ratio` or
# `:phase`), `nothing` when the branch works on the legacy grid.
function _derived_tap_model(br::Branch, which::Symbol)
  (br.taps_derived && br.tap_winding !== nothing) || return nothing
  return which === :ratio ? br.tap_winding.taps : br.tap_winding.phase_taps
end

# Value of a model step of a derived branch: the resolved ratio (:ratio)
# or shift (:phase). The tap estimation fixes a typed model on these values.
function _model_step_value(br::Branch, which::Symbol, step::Int)::Float64
  w = br.tap_winding
  res = which === :ratio ? resolve_winding_taps(w; tap_changer_model = br.tap_correction, ratio_step = step) : resolve_winding_taps(w; tap_changer_model = br.tap_correction, phase_step = step)
  return which === :ratio ? res.ratio : res.shift_deg
end

# Controllers move STEPS on a modelled transformer (0.20.0): the model step
# of a derived branch that brings the resolved ratio (:ratio) or shift
# (:phase) closest to `target`, among the current step and its two
# neighbours; moves the step, re-resolves the branch fields through the
# resolver and reports whether it moved. A neighbour whose value lies
# outside `window` (the control range of the controller) is no candidate.
# A current step outside the window jumps to the step nearest the window
# instead (_model_step_toward_window), whatever `target` says: a tap that
# starts outside its range goes back in with its first move.
function _move_model_step!(br::Branch, which::Symbol, target::Float64; window::Tuple{Float64,Float64} = (-Inf, Inf))::Bool
  m = _derived_tap_model(br, which)
  m === nothing && return false
  outside = v -> max(window[1] - v, v - window[2], 0.0)
  best = m.step
  v_now = _model_step_value(br, which, m.step)
  if outside(v_now) > 1e-12
    best = _model_step_toward_window(br, which, window)
  else
    bestd = abs(v_now - target)
    for cand in (m.step - 1, m.step + 1)
      (cand < m.lowStep || cand > m.highStep) && continue
      v = _model_step_value(br, which, cand)
      outside(v) <= 1e-12 || continue
      d = abs(v - target)
      if d < bestd - 1e-12
        best = cand
        bestd = d
      end
    end
  end
  best == m.step && return false
  m.step = best
  # the controller set this step: no longer the assumed middle position
  _clear_tap_step_assumed!(br, which)
  resolve_branch_taps!(br, br.tap_winding; tap_changer_model = br.tap_correction)
  return true
end

# Ratio `nsteps` table points away from `ratio` on a branch with
# `tap_ratio_points` (positive: higher ratio), stopping at the first or last
# point. A ratio between two points (fractional start step) counts the next
# point in the direction as the first step. Only points inside `window`
# (the control range of the controller) count, so a move stops at the last
# point inside it; the default window takes every point.
function _ratio_points_shift(br::Branch, ratio::Float64, nsteps::Int; window::Tuple{Float64,Float64} = (-Inf, Inf))::Float64
  pts = br.tap_ratio_points
  tol = 1e-12
  inside = p -> window[1] - tol <= p <= window[2] + tol
  r = ratio
  for _ in 1:abs(nsteps)
    k = nsteps > 0 ? findfirst(p -> p > r + tol && inside(p), pts) : findlast(p -> p < r - tol && inside(p), pts)
    k === nothing && break
    r = pts[k]
  end
  return r
end

# Ratio `nsteps` steps away from `ratio` on the grid a controller moves on:
# table points on a tabular changer, whole `tap_step`s on the uniform grid,
# kept inside `window`. Without a window a table ends at its last point and
# the uniform grid at tap_min/tap_max. The tap estimation fixes a CGMES
# changer on the same positions.
function _ratio_grid_shift(br::Branch, ratio::Float64, nsteps::Real; window::Union{Nothing,Tuple{Float64,Float64}} = nothing)::Float64
  br.tap_ratio_points !== nothing && return _ratio_points_shift(br, ratio, Int(nsteps); window = something(window, (-Inf, Inf)))
  lo, hi = something(window, (br.tap_min, br.tap_max))
  return clamp(ratio + nsteps * br.tap_step, lo, hi)
end

# Signed number of steps a ratio move from `old` to `new` makes on the grid
# of `br`: table points passed on a tabular changer, whole `tap_step`s on
# the uniform grid, the sign alone without a step.
function _ratio_steps_moved(br::Branch, old::Float64, new::Float64)::Int
  if br.tap_ratio_points !== nothing
    lo, hi = minmax(old, new)
    n = count(p -> p > lo + 1e-12 && p <= hi + 1e-12, br.tap_ratio_points)
    return new >= old ? n : -n
  end
  return br.tap_step > 0.0 ? round(Int, (new - old) / br.tap_step) : Int(sign(new - old))
end

# Move the model step of a derived branch by a signed number of steps,
# clamped to the model range; returns whether the step changed.
function _shift_model_steps!(br::Branch, which::Symbol, nsteps::Int)::Bool
  m = _derived_tap_model(br, which)
  m === nothing && return false
  target = clamp(m.step + nsteps, m.lowStep, m.highStep)
  target == m.step && return false
  m.step = target
  resolve_branch_taps!(br, br.tap_winding; tap_changer_model = br.tap_correction)
  return true
end

# On a modelled transformer the step grid rarely meets a control range
# exactly: a tap inside its range stands at that range's limit when a
# neighbouring model step resolves to a ratio outside it. Without a control
# range the check stays the ratio comparison with tap_min/tap_max alone.
function _model_step_at_window_edge(ctrl::PowerTransformerControl, br::Branch)::Bool
  (ctrl.tap_min === nothing && ctrl.tap_max === nothing) && return false
  m = _derived_tap_model(br, :ratio)
  m === nothing && return false
  lo, hi = _tap_control_range(ctrl, br)
  (lo - 1e-12 <= br.tap_ratio <= hi + 1e-12) || return false
  for cand in (m.step - 1, m.step + 1)
    (cand < m.lowStep || cand > m.highStep) && continue
    r = resolve_winding_taps(br.tap_winding; tap_changer_model = br.tap_correction, ratio_step = cand).ratio
    (lo - 1e-12 <= r <= hi + 1e-12) || return true
  end
  return false
end

# The same on a tabular changer (tap_ratio_points): a ratio inside the
# control range stands at its limit when no table point inside the range
# lies beyond it in one direction. Without a control range the ratio
# comparison with tap_min/tap_max (the first and last point) decides.
function _table_point_at_window_edge(ctrl::PowerTransformerControl, br::Branch)::Bool
  (ctrl.tap_min === nothing && ctrl.tap_max === nothing) && return false
  br.tap_ratio_points === nothing && return false
  lo, hi = _tap_control_range(ctrl, br)
  (lo - 1e-12 <= br.tap_ratio <= hi + 1e-12) || return false
  r = br.tap_ratio
  return _ratio_points_shift(br, r, 1; window = (lo, hi)) == r || _ratio_points_shift(br, r, -1; window = (lo, hi)) == r
end

# The model step nearest to `window`: the step with the least distance
# outside the window over the whole model range, among equals the one whose
# value is closest to the current step's. A window between two grid points
# gives the step next to it; the current step wins a tie, so a step that is
# already nearest stays.
function _model_step_toward_window(br::Branch, which::Symbol, window::Tuple{Float64,Float64})::Int
  m = _derived_tap_model(br, which)
  outside = v -> max(window[1] - v, v - window[2], 0.0)
  v_now = _model_step_value(br, which, m.step)
  best = m.step
  best_out = outside(v_now)
  best_d = 0.0
  for step in m.lowStep:m.highStep
    v = _model_step_value(br, which, step)
    o = outside(v)
    d = abs(v - v_now)
    if o < best_out - 1e-12 || (abs(o - best_out) <= 1e-12 && d < best_d - 1e-12)
      best, best_out, best_d = step, o, d
    end
  end
  return best
end

# Controllers whose start outside the control range a run has reported. A
# QSTS run solves its base case and every step on its own copy of the same
# net, each with its own control loop, so every step would report the same
# start again; the run sets one shared record here for all of them (its
# parallel steps see it through the scope). A run_control! outside such a
# scope has none and reports once anyway, as its first pass moves the
# ratio into the range.
struct _TapStartWarnings
  seen::Set{String}
  lock::ReentrantLock
end
_TapStartWarnings() = _TapStartWarnings(Set{String}(), ReentrantLock())
const _TAP_START_WARNINGS = Base.ScopedValues.ScopedValue{Union{Nothing,_TapStartWarnings}}(nothing)

# true when the warning for controller `key` is the first of the run
function _first_tap_start_warning!(key::String)::Bool
  w = _TAP_START_WARNINGS[]
  w === nothing && return true
  return lock(w.lock) do
    key in w.seen ? false : (push!(w.seen, key); true)
  end
end

# Warning for a ratio that starts outside the control range and goes onto
# its nearer end in this pass: transformer, start ratio and step, the limit
# it lies beyond and the ratio and step it goes to, once per controller and
# run. The plain step grid counts positions from the ratio 1.0, as the
# report rows do; a modelled transformer names its model step, a tabular
# changer the table point.
function _warn_start_outside_control_range(ctrl::PowerTransformerControl, br::Branch, lo::Float64, hi::Float64)
  _first_tap_start_warning!(control_name(ctrl)) || return nothing
  m = _derived_tap_model(br, :ratio)
  position = function (ratio::Float64, step::Union{Nothing,Int})
    step !== nothing && return string("model step ", step)
    if br.tap_ratio_points !== nothing
      k = findfirst(q -> abs(q - ratio) <= 1e-12, br.tap_ratio_points)
      return k === nothing ? "between table points" : string("table point ", k, " of ", length(br.tap_ratio_points))
    end
    p = _tap_position(ratio, 1.0, br.tap_step)
    return ismissing(p) ? "no tap position" : @sprintf("tap position %+d", p)
  end
  start_label = position(br.tap_ratio, m === nothing ? nothing : m.step)
  if m === nothing
    new_ratio = _ratio_toward_control_range(br, lo, hi)
    new_label = position(new_ratio, nothing)
  else
    new_step = _model_step_toward_window(br, :ratio, (lo, hi))
    new_ratio = _model_step_value(br, :ratio, new_step)
    new_label = position(new_ratio, new_step)
  end
  limit = br.tap_ratio < lo ? string("below its control range (tap_min = ", lo, ")") : string("above its control range (tap_max = ", hi, ")")
  @warn @sprintf("PowerTransformerControl: transformer %s starts at ratio %.5f (%s) %s and goes to %.5f (%s) in the first pass.", br.comp.cName, br.tap_ratio, start_label, limit, new_ratio, new_label)
  return nothing
end

# A ratio outside the control range of a voltage controller (the case
# ratio left out by a narrowed range). On a modelled transformer it counts
# only while a model step nearer to the range exists, so a range between
# two grid points does not hold the controller open.
function _ratio_outside_control_range(ctrl::PowerTransformerControl, br::Branch)::Bool
  (ctrl.tap_min === nothing && ctrl.tap_max === nothing) && return false
  (ctrl.control_ratio && ctrl.mode in (:voltage, :voltage_and_branch_active_power)) || return false
  lo, hi = _tap_control_range(ctrl, br)
  lo - 1e-12 <= br.tap_ratio <= hi + 1e-12 && return false
  m = _derived_tap_model(br, :ratio)
  m === nothing || return _model_step_toward_window(br, :ratio, (lo, hi)) != m.step
  return _ratio_toward_control_range(br, lo, hi) != br.tap_ratio
end

# Ratio a tap outside its control range goes to on a branch without a tap
# model: the nearer end of the range on the plain step grid, on a tabular
# changer the table point nearest the range (among equals the one closest
# to the current ratio, the current ratio winning a tie), so the tap stays
# on its table.
function _ratio_toward_control_range(br::Branch, lo::Float64, hi::Float64)::Float64
  pts = br.tap_ratio_points
  pts === nothing && return clamp(br.tap_ratio, lo, hi)
  outside = v -> max(lo - v, v - hi, 0.0)
  best = br.tap_ratio
  best_out = outside(best)
  best_d = 0.0
  for p in pts
    o = outside(p)
    d = abs(p - br.tap_ratio)
    if o < best_out - 1e-12 || (abs(o - best_out) <= 1e-12 && d < best_d - 1e-12)
      best, best_out, best_d = p, o, d
    end
  end
  return best
end

# The absolute write paths of a tap position on a branch without a tap
# model: the ratio goes to both ratio fields (the solver reads tap_ratio, the
# branch flows read ratio), the shift to both angle fields with the X(alpha)
# of the controller's phase tap model where it has one (the next solve
# re-stamps the Y-bus from x_pu). A modelled branch moves its model steps
# instead (_set_model_tap_steps!). control_apply_update! and the carry of a
# sequential QSTS run (qsts_sequential.jl) write through these. A write
# that moves the tap takes the mark of an assumed step off the changer
# (_clear_tap_step_assumed!); the probes write the fields directly and put
# them back, so they keep it.
function _set_ratio_tap!(br::Branch, ratio::Float64)
  ratio == br.tap_ratio || _clear_tap_step_assumed!(br, :ratio)
  br.tap_ratio = ratio
  br.ratio = ratio
  return nothing
end
function _set_phase_shift!(net::Net, ctrl::PowerTransformerControl, br::Branch, shift_deg::Float64)
  shift_deg == br.phase_shift_deg || _clear_tap_step_assumed!(br, :phase)
  br.phase_shift_deg = shift_deg
  br.angle = shift_deg
  x_new = _phase_tap_reactance_at(net, ctrl, shift_deg)
  x_new === nothing || (br.x_pu = x_new)
  return nothing
end
# the model steps of a derived branch (typemin(Int) keeps an actuator as it
# is), re-resolved once through the resolver as _move_model_step! does
function _set_model_tap_steps!(br::Branch, ratio_step::Int, phase_step::Int)
  changed = false
  for (which, step) in ((:ratio, ratio_step), (:phase, phase_step))
    step == typemin(Int) && continue
    m = _derived_tap_model(br, which)
    m === nothing && continue
    m.step == step || _clear_tap_step_assumed!(br, which)
    m.step = step
    changed = true
  end
  changed && resolve_branch_taps!(br, br.tap_winding; tap_changer_model = br.tap_correction)
  return changed
end

function control_apply_update!(ctrl::PowerTransformerControl, net::Net, ::AbstractControlState, update::NamedTuple, context)::Bool
  br = _find_trafo_branch(net, ctrl.trafo)
  moved = false
  if update.new_ratio != update.old_ratio
    if _derived_tap_model(br, :ratio) !== nothing
      moved = _move_model_step!(br, :ratio, update.new_ratio; window = _tap_control_range(ctrl, br)) || moved
    else
      _set_ratio_tap!(br, update.new_ratio)
      moved = true
    end
    # master/slave group (#322): mirror the master's move onto every
    # follower STEP-synchronously (whole steps of the follower's own
    # tap_step, so units with different neutral ratios stay aligned in
    # positions, not in absolute ratios), clamped to the follower's range.
    for fname in ctrl.followers
      fbr = _find_trafo_branch(net, fname)
      fmodel = _derived_tap_model(fbr, :ratio)
      if fmodel !== nothing
        # a modelled follower moves the same number of steps in the
        # direction of the master's ratio change on its own step sign
        nsteps = _ratio_steps_moved(br, update.old_ratio, update.new_ratio)
        # the probes shift and shift back through the same helper, so the
        # mark of an assumed step comes off here, not in the helper
        _shift_model_steps!(fbr, :ratio, nsteps * Int(sign(fmodel.tapStepPercent == 0.0 ? 1.0 : fmodel.tapStepPercent))) && _clear_tap_step_assumed!(fbr, :ratio)
        continue
      end
      if br.tap_ratio_points !== nothing || fbr.tap_ratio_points !== nothing
        # a tabular master or follower: count steps on the master's grid and
        # make as many on the follower's own grid
        nsteps = _ratio_steps_moved(br, update.old_ratio, update.new_ratio)
        fnew = _ratio_grid_shift(fbr, fbr.tap_ratio, nsteps)
      elseif br.tap_step > 0.0 && fbr.tap_step > 0.0
        nsteps = round((update.new_ratio - update.old_ratio) / br.tap_step)
        fnew = _ratio_grid_shift(fbr, fbr.tap_ratio, nsteps)
      else
        fnew = clamp(fbr.tap_ratio + (update.new_ratio - update.old_ratio), fbr.tap_min, fbr.tap_max)
      end
      _set_ratio_tap!(fbr, fnew)
    end
  end
  if update.new_phase != update.old_phase
    if _derived_tap_model(br, :phase) !== nothing
      # the resolver applies the model's X(alpha) with the step
      moved = _move_model_step!(br, :phase, update.new_phase) || moved
    else
      # tap-dependent reactance X(α) with the angle
      _set_phase_shift!(net, ctrl, br, update.new_phase)
      moved = true
    end
  end
  # tap changer tables with r/x/g/b deviations: the impedance follows the
  # step the move reached, on the master and on every follower (a follower
  # mirrors a ratio change even when the master's own step stayed)
  if update.new_ratio != update.old_ratio || update.new_phase != update.old_phase
    _apply_tap_step_impedance!(br)
    for fname in ctrl.followers
      _apply_tap_step_impedance!(_find_trafo_branch(net, fname))
    end
  end
  ctrl.at_limit = _tap_at_range_end(ctrl, br)
  ctrl.status = ctrl.at_limit ? :at_limit : ctrl.status
  return moved
end
# Filter by the per-controller row name, not the transformer id: with two
# disjoint-actuator controllers on one transformer (split Schrägregelung),
# a transformer-id filter would duplicate both rows into each controller.
function control_report_rows(ctrl::PowerTransformerControl, net::Net, ::AbstractControlState, context)
  br = _find_trafo_branch(net, ctrl.trafo)
  own = string(br.comp.cName, " ", _controller_type_label(ctrl, br))
  return filter(row -> row.controller_name == own, buildTapControllerReportRows(net))
end
function control_trace_rows(ctrl::PowerTransformerControl, net::Net, ::AbstractControlState, context)
  br = _find_trafo_branch(net, ctrl.trafo)
  return [(
    outer_iteration = context.outer_iteration,
    controller_name = control_name(ctrl),
    controller_type = "PowerTransformerControl",
    transformer_id = br.comp.cName,
    mode = string(ctrl.mode),
    status = ctrl.status,
    converged = ctrl.converged,
    at_limit = ctrl.at_limit,
    achieved_vm_pu = isnothing(ctrl.achieved_vm_pu) ? missing : ctrl.achieved_vm_pu,
    target_vm_pu = isnothing(ctrl.target_vm_pu) ? missing : ctrl.target_vm_pu,
    achieved_p_mw = isnothing(ctrl.achieved_p_mw) ? missing : ctrl.achieved_p_mw,
    target_p_mw = isnothing(ctrl.p_target_mw) ? missing : ctrl.p_target_mw,
    tap_ratio = br.tap_ratio,
    phase_shift_deg = br.phase_shift_deg,
  )]
end

function control_element_descriptor(ctrl::PowerTransformerControl, net::Net)::Union{Nothing,NamedTuple}
  br = _find_trafo_branch(net, ctrl.trafo)
  ratio = ctrl.control_ratio && br.has_ratio_tap
  phase = ctrl.control_phase && br.has_phase_tap
  actuator = ratio && phase ? :tap_ratio_and_phase_shift : ratio ? :tap_ratio : phase ? :phase_shift_deg : :none
  quantity = ctrl.mode === :voltage ? :bus_voltage : ctrl.mode === :branch_active_power ? :branch_active_power : :bus_voltage_and_branch_active_power
  target = ctrl.target_bus !== nothing ? ctrl.target_bus : ctrl.target_branch === nothing ? "" : string(ctrl.target_branch[1], "->", ctrl.target_branch[2])
  target_value = ctrl.target_vm_pu !== nothing ? ctrl.target_vm_pu : ctrl.p_target_mw
  amin, amax = phase && !ratio ? (br.phase_min_deg, br.phase_max_deg) : _tap_control_range(ctrl, br)
  # a master of a group names its unit count: "T1 (+2 followers)"
  element_label = isempty(ctrl.followers) ? br.comp.cName : string(br.comp.cName, " (+", length(ctrl.followers), " follower", length(ctrl.followers) == 1 ? "" : "s", ")")
  return (
    name = control_name(ctrl),
    element = element_label,
    device = _controller_type_label(ctrl, br),
    actuator = actuator,
    actuator_min = amin,
    actuator_max = amax,
    quantity = quantity,
    target = target,
    target_value = target_value === nothing ? missing : target_value,
    discrete = ctrl.is_discrete,
    enabled = ctrl.enabled,
    status = ctrl.status,
    converged = ctrl.converged,
    at_limit = ctrl.at_limit,
  )
end

function _controller_type_label(ctrl::PowerTransformerControl, br::Branch)::String
  ratio = ctrl.control_ratio && br.has_ratio_tap
  phase = ctrl.control_phase && br.has_phase_tap
  if ratio && phase
    return "OLTC+PST"
  elseif ratio
    return "OLTC"
  elseif phase
    return "PST"
  end
  return "fixed transformer"
end

function _controller_status_label(ctrl::PowerTransformerControl)::String
  !ctrl.enabled && return "inactive"
  ctrl.converged && return "converged"
  ctrl.at_limit && return "at_limit"
  ctrl.status == :max_outer_iters && return "not_converged"
  return "active"
end

@inline function _tap_position(value::Float64, neutral::Float64, step::Float64)
  step <= 0.0 && return missing
  return round(Int, (value - neutral) / step)
end

# The achieved P of a report row comes from the stored branch flows, and only
# calcNetLosses! writes them: a bare runpf! leaves them unset. The row then
# reports missing, like the branch table reports 0.0 for the same branch,
# instead of failing on the unset flow.
function _reported_branch_p_mw(net::Net, from_bus::String, to_bus::String)
  br = getNetBranch(net = net, fromBus = from_bus, toBus = to_bus)
  (isnothing(br.fBranchFlow) || isnothing(br.tBranchFlow)) && return missing
  p = get_branch_p_from_to_mw(net, from_bus, to_bus)
  return isnothing(p) ? missing : p
end

"""
    buildTapControllerReportRows(net::Net)

Build typed, machine-readable tap-controller report rows with engineering fields
for textual and DataFrame-based reporting. `ratio_tap_min`/`ratio_tap_max`
give the range the controller moves the ratio in: its control range where
set, the physical range of the transformer otherwise. `achieved_p_mw` reads
the stored branch flows; it is `missing` until `calcNetLosses!` has computed
them (a bare `runpf!` does not).
"""
function buildTapControllerReportRows(net::Net)::Vector{NamedTuple}
  rows = NamedTuple[]
  for ctrl in _tap_controllers(net)
    br = _find_trafo_branch(net, ctrl.trafo)
    target_bus = ctrl.target_bus
    target_branch = ctrl.target_branch
    achieved_vm = isnothing(target_bus) ? missing : get_bus_vm_pu(net, target_bus)
    achieved_p = isnothing(target_branch) ? missing : _reported_branch_p_mw(net, target_branch[1], target_branch[2])
    push!(
      rows,
      (
        controller_name = string(br.comp.cName, " ", _controller_type_label(ctrl, br)),
        transformer_id = br.comp.cName,
        transformer_branch_index = br.branchIdx,
        from_bus = br.fromBus,
        to_bus = br.toBus,
        control_type = _controller_type_label(ctrl, br),
        mode = _controller_mode_label(ctrl.mode),
        target_bus = isnothing(target_bus) ? missing : target_bus,
        target_vm_pu = isnothing(ctrl.target_vm_pu) ? missing : ctrl.target_vm_pu,
        achieved_vm_pu = achieved_vm,
        target_branch_from = isnothing(target_branch) ? missing : target_branch[1],
        target_branch_to = isnothing(target_branch) ? missing : target_branch[2],
        p_target_mw = isnothing(ctrl.p_target_mw) ? missing : ctrl.p_target_mw,
        achieved_p_mw = achieved_p,
        tap_ratio = br.tap_ratio,
        phase_shift_deg = br.phase_shift_deg,
        ratio_tap_position = ctrl.is_discrete && br.has_ratio_tap ? (_derived_tap_model(br, :ratio) !== nothing ? _derived_tap_model(br, :ratio).step : _tap_position(br.tap_ratio, 1.0, br.tap_step)) : missing,
        phase_tap_position = ctrl.is_discrete && br.has_phase_tap ? (_derived_tap_model(br, :phase) !== nothing ? _derived_tap_model(br, :phase).step : _tap_position(br.phase_shift_deg, 0.0, br.phase_step_deg)) : missing,
        ratio_tap_min = br.has_ratio_tap ? _tap_control_range(ctrl, br)[1] : missing,
        ratio_tap_max = br.has_ratio_tap ? _tap_control_range(ctrl, br)[2] : missing,
        ratio_tap_step = br.has_ratio_tap ? br.tap_step : missing,
        phase_tap_min_deg = br.has_phase_tap ? br.phase_min_deg : missing,
        phase_tap_max_deg = br.has_phase_tap ? br.phase_max_deg : missing,
        phase_tap_step_deg = br.has_phase_tap ? br.phase_step_deg : missing,
        discrete = ctrl.is_discrete,
        converged = ctrl.converged,
        at_limit = ctrl.at_limit,
        status = _controller_status_label(ctrl),
        tap_model_kind = br.taps_derived && br.tap_winding !== nothing ? _tap_model_kind_label(br.tap_winding) : "-",
        power_direction = isnothing(target_branch) ? missing : string(target_branch[1], " -> ", target_branch[2]),
      ),
    )
  end
  return rows
end

"""
    addPowerTransformerControl!(net; ...)

Add and validate a transformer tap controller.

Validation rules:
- Each actuator (ratio tap, phase tap) of a transformer may be driven by at
  most one active controller. Two controllers may coexist on one transformer
  when they drive disjoint actuators, the classic Schrägregelung split of a
  voltage controller (ratio tap) plus an active-power controller (phase tap).
- Required target fields must be set according to `mode`.
- `control_ratio` / `control_phase` must match the chosen mode.

`priority::Int = 0` is the rank in the outer loop (`control_priority`): a
controller with a larger number waits while this one still moves its tap.
A tap counts as at its limit when its ratio stands at `tap_min`/`tap_max`
(or its shift at the end of the phase range); a last partial step onto a
limit that is off the step grid still counts as a move.

`tap_min`/`tap_max` (ratios, default `nothing`) give the controller its own
control range inside the physical range of the transformer: the controller
moves the ratio only within it and stands at its limit at either end. A
value outside the physical range, a lower end not below the upper one, or a
range on a controller without ratio control throws an `ArgumentError`
naming the transformer. A side left at `nothing` keeps the physical limit.
A ratio outside the range (for example the case ratio when the range
leaves it out) goes onto the nearer end of the range in the first pass,
whatever the voltage error, and the controller does not count as converged
until it is inside; on a modelled transformer it goes to the model step
nearest the range, on a tabular changer to the table point nearest it.
That move logs a warning naming the transformer, the start ratio and
step, the limit and the ratio it goes to, once per controller and run; a
QSTS run (`runQsts`) is one run over its base case and all steps. The
probe that finds
the direction of a step still uses the physical range; followers of a
group keep their own physical range.
"""
function addPowerTransformerControl!(
  net::Net;
  trafo::String,
  followers::Vector{String} = String[],
  mode::Symbol,
  target_bus::Union{Nothing,String} = nothing,
  target_branch::Union{Nothing,Tuple{String,String}} = nothing,
  target_vm_pu::Union{Nothing,Float64} = nothing,
  p_target_mw::Union{Nothing,Float64} = nothing,
  q_target_mvar::Union{Nothing,Float64} = nothing,
  control_ratio::Bool = true,
  control_phase::Bool = false,
  is_discrete::Bool = true,
  deadband_vm_pu::Float64 = 1e-3,
  deadband_p_mw::Float64 = 0.5,
  voltage_error_metric::Symbol = :vm,
  max_outer_iters::Int = 20,
  enabled::Bool = true,
  priority::Int = 0,
  tap_min::Union{Nothing,Real} = nothing,
  tap_max::Union{Nothing,Real} = nothing,
)
  br = _find_trafo_branch(net, trafo)
  trafo_obj = findfirst(t -> t.comp.cFrom_bus == br.comp.cFrom_bus && t.comp.cTo_bus == br.comp.cTo_bus, net.trafos)
  # Per-actuator exclusivity: a second active controller on the same
  # transformer is allowed only when it drives a disjoint actuator set
  # (e.g. voltage/ratio + active-power/phase on one Schrägregler unit).
  # A transformer that FOLLOWS in a master/slave group counts as occupied:
  # its ratio tap is driven by the group's master.
  for c in _tap_controllers(net)
    c.enabled || continue
    for fname in c.followers
      fbr = _find_trafo_branch(net, fname)
      if fbr.branchIdx == br.branchIdx && control_ratio
        error("PowerTransformerControl: transformer $(trafo) already follows the group of $(c.trafo); a follower's ratio tap is driven by its master.")
      end
    end
    cbr = _find_trafo_branch(net, c.trafo)
    cbr.branchIdx == br.branchIdx || continue
    if (control_ratio && c.control_ratio) || (control_phase && c.control_phase)
      error("PowerTransformerControl: transformer $(trafo) already has an active controller for this actuator (ratio/phase); use disjoint actuators or disable the existing controller.")
    end
  end
  # two independent tap controllers steering ONE bus voltage fight each
  # other through circulating reactive power (#322); the group (followers)
  # is the supported way to regulate parallel units
  if mode in (:voltage, :voltage_and_branch_active_power) && target_bus !== nothing
    for c in _tap_controllers(net)
      c.enabled || continue
      cbr = _find_trafo_branch(net, c.trafo)
      cbr.branchIdx == br.branchIdx && continue
      if c.target_bus == target_bus && c.target_vm_pu !== nothing
        @warn "PowerTransformerControl: a tap controller ($(c.trafo)) already regulates the voltage of target bus $(target_bus). Two independent controllers on one voltage fight each other via circulating reactive power; regulate parallel units as a GROUP (followers = [...]) instead."
      end
    end
  end
  # master/slave group validation (#322): voltage/ratio groups only; every
  # follower needs its own ratio-tap machinery, must be a different unit,
  # and must not be regulated or followed elsewhere
  if !isempty(followers)
    mode == :voltage || error("PowerTransformerControl: followers are supported for mode = :voltage groups (parallel transformers on one busbar); got mode = $(mode).")
    seen_followers = Set{Int}()
    for fname in followers
      fbr = _find_trafo_branch(net, fname)
      fbr.branchIdx == br.branchIdx && error("PowerTransformerControl: follower $(fname) is the master transformer itself.")
      fbr.has_ratio_tap || error("PowerTransformerControl: follower $(fname) has no ratio-tap machinery (has_ratio_tap = false).")
      fbr.branchIdx in seen_followers && error("PowerTransformerControl: follower $(fname) listed twice.")
      push!(seen_followers, fbr.branchIdx)
      for c in _tap_controllers(net)
        c.enabled || continue
        cbr = _find_trafo_branch(net, c.trafo)
        cbr.branchIdx == fbr.branchIdx && c.control_ratio && error("PowerTransformerControl: follower $(fname) already has its own active ratio controller; a unit cannot regulate and follow at once.")
        for other in c.followers
          obr = _find_trafo_branch(net, other)
          obr.branchIdx == fbr.branchIdx && error("PowerTransformerControl: follower $(fname) already follows the group of $(c.trafo).")
        end
      end
    end
  end
  mode in (:voltage, :branch_active_power, :voltage_and_branch_active_power) || error("PowerTransformerControl: unsupported mode=$(mode)")

  if mode in (:voltage, :voltage_and_branch_active_power)
    isnothing(target_bus) && error("PowerTransformerControl: target_bus is required for mode=$(mode)")
    isnothing(target_vm_pu) && error("PowerTransformerControl: target_vm_pu is required for mode=$(mode)")
    !control_ratio && error("PowerTransformerControl: control_ratio must be true for voltage control")
  end
  if mode in (:branch_active_power, :voltage_and_branch_active_power)
    isnothing(target_branch) && error("PowerTransformerControl: target_branch is required for mode=$(mode)")
    isnothing(p_target_mw) && error("PowerTransformerControl: p_target_mw is required for mode=$(mode)")
    !control_phase && error("PowerTransformerControl: control_phase must be true for branch active power control")
  end

  ctrl = PowerTransformerControl(;
    trafo = trafo,
    followers = followers,
    mode = mode,
    target_bus = target_bus,
    target_branch = target_branch,
    target_vm_pu = target_vm_pu,
    p_target_mw = p_target_mw,
    q_target_mvar = q_target_mvar,
    control_ratio = control_ratio,
    control_phase = control_phase,
    is_discrete = is_discrete,
    deadband_vm_pu = deadband_vm_pu,
    deadband_p_mw = deadband_p_mw,
    voltage_error_metric = voltage_error_metric,
    max_outer_iters = max_outer_iters,
    enabled = enabled,
    priority = priority,
  )
  _set_tap_control_range!(ctrl, br, tap_min, tap_max)
  if isnothing(trafo_obj)
    error("PowerTransformerControl: no transformer object found for branch $(br.branchIdx)")
  end
  push!(net.trafos[trafo_obj].side1.controls, ctrl)
  return net
end

"""
    addTapController!(net; kwargs...)

Attach a tap controller to a transformer; alias of
[`addPowerTransformerControl!`](@ref).
"""
addTapController!(net::Net; kwargs...) = addPowerTransformerControl!(net; kwargs...)

# One probe solve, returns the solver status (0 = converged). On a net with
# several AC islands a non-converged island does not come back as a status:
# the island path throws an ErrorException ("AC island N power-flow solve
# failed"). That is the expected failure here. It counts as a failed probe
# solve, as in the run's own solve (_runpf_with_config!), so the probe
# answers the same way on one island and on several. A non-converged
# one-island probe solve returns its status without a message, so the
# island case gets none either. Any other exception propagates.
function _probe_solve_status(net::Net, config::PowerFlowConfig)::Int
  _, erg = try
    _runpf_config_once!(net, config)
  catch err
    err isa ErrorException || rethrow()
    (0, 1)
  end
  return erg
end

"""
    _phase_probe_direction(net, br, ctrl, config) -> Float64

Internal helper: determines the empirical sign of `ΔP_from_to` for a positive
phase increment (`+phase_step_deg`) on the controlled transformer and branch.
Both probe solves run under the run's power-flow configuration `config`
(start mode, island switches, Q limits, Newton settings), one attempt each,
without the rescue ladder and the DC fallback of the run's own solve.
Returns `-1`, `0`, or `+1`: `-1` when the first solve fails, `0` when the
perturbed solve fails (also a failed island on a net with several AC
islands). The branch values the probe changes are restored in every case.
"""
function _phase_probe_direction(net::Net, br::Branch, ctrl::PowerTransformerControl, config::PowerFlowConfig)
  # power_flow.tol_MW is converted against this net's base, as in the run's own solve
  probe_cfg = _resolve_tolerance_for_net(config, net)
  oldphi = br.phase_shift_deg
  step = br.phase_step_deg
  erg = _probe_solve_status(net, probe_cfg)
  erg != 0 && return -1.0
  # branch flows live in fBranchFlow/tBranchFlow, which only calcNetLosses!
  # rewrites: without the refresh both probe reads return the same stale
  # values and the estimated direction is always 0
  calcNetLosses!(net)
  p0 = get_branch_p_from_to_mw(net, ctrl.target_branch[1], ctrl.target_branch[2])
  oldx = br.x_pu
  derived = _derived_tap_model(br, :phase) !== nothing
  if derived
    # a modelled transformer is probed one model step up; the answer is
    # expressed per degree of shift so the proposal in degrees keeps its sign
    _shift_model_steps!(br, :phase, 1) || return 0.0
  else
    br.phase_shift_deg = clamp(oldphi + step, br.phase_min_deg, br.phase_max_deg)
    br.angle = br.phase_shift_deg
    # perturb the reactance consistently with the apply step (X(α) coupling);
    # without a typed model the probe keeps today's static-x behaviour
    x_probe = _phase_tap_reactance_at(net, ctrl, br.phase_shift_deg)
    x_probe === nothing || (br.x_pu = x_probe)
  end
  dphi = br.phase_shift_deg - oldphi
  p1 = p0
  try
    if _probe_solve_status(net, probe_cfg) == 0
      # the flows are read before the restore: calcNetLosses! uses the
      # branch parameters, which must still be the probe's
      calcNetLosses!(net)
      p1 = get_branch_p_from_to_mw(net, ctrl.target_branch[1], ctrl.target_branch[2])
    end
  finally
    # the probe value never stays on the branch, also when the solve throws
    if derived
      _shift_model_steps!(br, :phase, -1)
    else
      br.phase_shift_deg = oldphi
      br.angle = oldphi
      br.x_pu = oldx
    end
  end
  return sign(p1 - p0) * (dphi == 0.0 ? 1.0 : sign(dphi))
end

"""
    _ratio_probe_direction(net, br, ctrl, config) -> Float64

Internal helper: determines the empirical sign of `ΔVm_target` for a positive
ratio increment (`+tap_step`) on the controlled transformer and target bus.
Both probe solves run under the run's power-flow configuration `config`
(start mode, island switches, Q limits, Newton settings), one attempt each,
without the rescue ladder and the DC fallback of the run's own solve.
Returns `-1`, `0`, or `+1`: `-1` when the first solve fails, `0` when the
perturbed solve fails (also a failed island on a net with several AC
islands). The branch values the probe changes are restored in every case.
"""
function _ratio_probe_direction(net::Net, br::Branch, ctrl::PowerTransformerControl, config::PowerFlowConfig)
  # power_flow.tol_MW is converted against this net's base, as in the run's own solve
  probe_cfg = _resolve_tolerance_for_net(config, net)
  oldratio = br.tap_ratio
  step = br.tap_step
  erg = _probe_solve_status(net, probe_cfg)
  erg != 0 && return -1.0
  vm0 = get_bus_vm_pu(net, ctrl.target_bus)
  derived = _derived_tap_model(br, :ratio) !== nothing
  if derived
    # a modelled transformer is probed one model step up, the answer per
    # unit of ratio change
    _shift_model_steps!(br, :ratio, 1) || return 0.0
  else
    newratio = clamp(oldratio + step, br.tap_min, br.tap_max)
    if isapprox(newratio, oldratio; atol = 1e-12)
      return 0.0
    end
    br.tap_ratio = newratio
    br.ratio = newratio
  end
  dratio = br.tap_ratio - oldratio
  vm1 = vm0
  try
    _probe_solve_status(net, probe_cfg) == 0 && (vm1 = get_bus_vm_pu(net, ctrl.target_bus))
  finally
    # the probe value never stays on the branch, also when the solve throws
    if derived
      _shift_model_steps!(br, :ratio, -1)
    else
      br.tap_ratio = oldratio
      br.ratio = oldratio
    end
  end
  return sign(vm1 - vm0) * (dratio == 0.0 ? 1.0 : sign(dratio))
end

"""
    run_tap_controllers_outer!(net; ...)

Removed from the public control API; always throws. Use
`run_control!(net; pf_config = ..., control_config = ...)` instead.
"""
function run_tap_controllers_outer!(net::Net; kwargs...)
  throw(ArgumentError("run_tap_controllers_outer! has been removed from the public control API; use run_control!(net; pf_config=..., control_config=...)"))
end

"""
    printTapControllerSummary(io, net)

Print a compact controller summary block for all configured tap controllers.
"""
function printTapControllerSummary(io::IO, net::Net)
  if isempty(_tap_controllers(net))
    println(io, "\nControl")
    println(io, "-------")
    println(io, "Transformer controls: none")
    return
  end
  rows = buildTapControllerReportRows(net)
  println(io, "\nTransformer Control Summary")
  println(io, "---------------------------")
  println(io, "Power sign convention: achieved_p_mw is positive in the configured target branch direction (from -> to).")
  for row in rows
    println(io, row.controller_name, " (", row.transformer_id, ", ", row.from_bus, " -> ", row.to_bus, ")")
    println(io, "  controller type    : ", row.control_type)
    println(io, "  mode               : ", row.mode)
    println(io, "  target bus         : ", ismissing(row.target_bus) ? "-" : row.target_bus)
    println(io, "  target Vm          : ", ismissing(row.target_vm_pu) ? "-" : @sprintf("%.4f pu", row.target_vm_pu))
    println(io, "  achieved Vm        : ", ismissing(row.achieved_vm_pu) ? "-" : @sprintf("%.4f pu", row.achieved_vm_pu))
    println(io, "  target branch      : ", ismissing(row.target_branch_from) ? "-" : string(row.target_branch_from, " -> ", row.target_branch_to))
    println(io, "  target P           : ", ismissing(row.p_target_mw) ? "-" : @sprintf("%.3f MW", row.p_target_mw))
    println(io, "  achieved P         : ", ismissing(row.achieved_p_mw) ? "-" : @sprintf("%.3f MW", row.achieved_p_mw))
    println(io, "  tap ratio          : ", @sprintf("%.5f", row.tap_ratio))
    println(io, "  phase shift        : ", @sprintf("%.5f deg", row.phase_shift_deg))
    println(io, "  tap position       : ", ismissing(row.ratio_tap_position) ? "-" : @sprintf("%+d", row.ratio_tap_position))
    println(io, "  phase position     : ", ismissing(row.phase_tap_position) ? "-" : @sprintf("%+d", row.phase_tap_position))
    println(io, "  ratio range        : ", ismissing(row.ratio_tap_min) ? "-" : @sprintf("%.5f .. %.5f", row.ratio_tap_min, row.ratio_tap_max))
    println(io, "  ratio step         : ", ismissing(row.ratio_tap_step) ? "-" : @sprintf("%.5f", row.ratio_tap_step))
    println(io, "  phase range        : ", ismissing(row.phase_tap_min_deg) ? "-" : @sprintf("%.5f .. %.5f deg", row.phase_tap_min_deg, row.phase_tap_max_deg))
    println(io, "  phase step         : ", ismissing(row.phase_tap_step_deg) ? "-" : @sprintf("%.5f deg", row.phase_tap_step_deg))
    println(io, "  discrete           : ", row.discrete)
    println(io, "  converged          : ", row.converged)
    println(io, "  at_limit           : ", row.at_limit)
    println(io, "  status             : ", row.status)
    if !row.converged && row.at_limit
      println(io, "  status detail      : target not fully reached because tap/phase limit was hit")
    end
  end
end

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

# file: src/series_reactance_control.jl
# purpose: TCSC-like series-reactance controller (outer loop, issue #297):
#          a variable series reactance on a line branch steers the branch
#          active power onto a target. The actuator is br.x_pu; every
#          accepted step changes one branch stamp and the outer loop
#          re-stamps the Y-bus before the next solve (same mechanism as the
#          X(alpha) coupling of phase-shifting transformers).

# Minimum admissible series impedance magnitude in p.u. The physical TCSC
# has a resonance region between capacitive and inductive operation that a
# real device must not enter; the model is a continuous clamped reactance,
# and this guard is the documented stand-in: it rejects reactance ranges
# whose impedance magnitude |r + jx| comes closer to zero than this bound
# anywhere in the range (the minimum over x lies at x = 0 when the range
# crosses sign, |z| = |r| there). Splitting the range into two disjoint
# admissible intervals is a possible later refinement, see the issue text.
const _SERIES_CTRL_EPS_Z = 1e-4

"""
    SeriesReactanceControl <: AbstractOuterController

TCSC-like series compensation on a line branch. The actuator is the branch
series reactance `x_pu` within `[x_min_pu, x_max_pu]`; the controlled
quantity is the active power carried by that same branch, measured in the
configured `fromBus` to `toBus` direction (the tap controller's
`achieved_p_mw` convention). In range the controller holds `p_target_mw`
via secant iteration on the scalar map x to P; at a range end the branch
behaves as a fixed compensated line and the controller reports honest
`at_limit`. Negative reactance (net capacitive branch) is admissible, the
impedance-magnitude guard `_SERIES_CTRL_EPS_Z` protects the singular
neighborhood of x = -0 (design choice, see the theory page).
"""
mutable struct SeriesReactanceControl <: AbstractOuterController
  name::String
  fromBus::String
  toBus::String
  branch_idx::Int
  p_target_mw::Float64
  x_min_pu::Float64
  x_max_pu::Float64
  x_pu::Float64
  deadband_p_mw::Float64
  max_outer_iters::Int
  enabled::Bool
  status::Symbol
  converged::Bool
  at_limit::Bool
  achieved_p_mw::Union{Nothing,Float64}
  prev_x_pu::Union{Nothing,Float64}
  prev_p_mw::Union{Nothing,Float64}
  outer_iters::Int
end

"""
    _series_reactance_controllers(net) -> Vector{SeriesReactanceControl}

Collect the series-reactance controllers stored on `net`. They share the
generic outer-controller registry (`net.machineControls`) with the machine
and shunt controllers; guarded with `hasproperty` because
`collect_outer_controllers` accepts `net::Any`.
"""
function _series_reactance_controllers(net)::Vector{SeriesReactanceControl}
  hasproperty(net, :machineControls) || return SeriesReactanceControl[]
  return SeriesReactanceControl[c for c in net.machineControls if c isa SeriesReactanceControl]
end

"""
    clearSeriesReactanceControllers!(net)

Remove all series-reactance controllers from `net` (other controllers stay).
"""
function clearSeriesReactanceControllers!(net::Net)
  filter!(c -> !(c isa SeriesReactanceControl), net.machineControls)
  return net
end

# range validation shared by registration: both range ends and, when the
# range crosses zero, the crossing itself must keep the series impedance
# magnitude above the guard
function _series_ctrl_check_range(r_pu::Float64, x_min_pu::Float64, x_max_pu::Float64)
  x_min_pu < x_max_pu || error("SeriesReactanceControl: x_min_pu ($(x_min_pu)) must be below x_max_pu ($(x_max_pu))")
  for x in (x_min_pu, x_max_pu)
    abs(complex(r_pu, x)) > _SERIES_CTRL_EPS_Z || error("SeriesReactanceControl: series impedance magnitude |$(r_pu) + j$(x)| is within the exclusion guard eps_z = $(_SERIES_CTRL_EPS_Z); the model refuses the resonance neighborhood.")
  end
  if x_min_pu < 0.0 < x_max_pu
    abs(r_pu) > _SERIES_CTRL_EPS_Z || error("SeriesReactanceControl: the reactance range crosses zero and the branch resistance |$(r_pu)| alone is within the exclusion guard eps_z = $(_SERIES_CTRL_EPS_Z); the model refuses the resonance neighborhood.")
  end
  return nothing
end

"""
    addSeriesReactanceControl!(net; fromBus, toBus, p_target_mw, x_min_pu, x_max_pu, ...)

Add a TCSC-like series-reactance controller to the line branch
`fromBus` to `toBus`. The outer control loop moves the branch series
reactance `x_pu` within `[x_min_pu, x_max_pu]` until the branch carries
`p_target_mw` in the `fromBus` to `toBus` direction.

# Arguments
- `net::Net`: the network.
- `fromBus::String`, `toBus::String`: terminals of the controlled line
  branch; the pair also fixes the measurement direction of the flow.
- `p_target_mw::Float64`: active-power target for the branch in MW.
- `x_min_pu::Float64`, `x_max_pu::Float64`: admissible reactance range in
  p.u. Negative values (net capacitive branch) are allowed.
- `deadband_p_mw::Float64 = 0.5`: convergence band around the target.
- `max_outer_iters::Int = 20`, `enabled::Bool = true`: outer-loop
  budget and switch.
- `name`: optional controller name (defaults to `TCSC_<from>_<to>`).

Fails with an error for a missing branch, a transformer branch (taps own
transformer reactance, see the X(alpha) coupling of the tap controller),
an inverted range, a range whose series impedance magnitude enters the
exclusion guard `eps_z`, a starting `x_pu` outside the range, or a second
series controller on the same branch.
"""
function addSeriesReactanceControl!(
  net::Net;
  fromBus::String,
  toBus::String,
  p_target_mw::Float64,
  x_min_pu::Float64,
  x_max_pu::Float64,
  deadband_p_mw::Float64 = 0.5,
  max_outer_iters::Int = 20,
  enabled::Bool = true,
  name::Union{Nothing,String} = nothing,
)
  br = getNetBranch(net = net, fromBus = fromBus, toBus = toBus)
  if br === nothing
    # the flow measurement resolves the branch in the registered orientation,
    # so a reversed registration would break at evaluation time; catch it
    # here with an actionable message instead
    rev = getNetBranch(net = net, fromBus = toBus, toBus = fromBus)
    rev === nothing && error("SeriesReactanceControl: no branch between $(fromBus) and $(toBus).")
    error("SeriesReactanceControl: the branch is stored in the orientation $(toBus) to $(fromBus); register the controller with that bus order (the pair also fixes the measurement direction of achieved_p_mw).")
  end
  # repo convention (see _find_trafo_branch): transformer branches carry a
  # nonzero ratio, line branches carry ratio 0. Taps own transformer
  # reactance, so the series controller only accepts line branches.
  br.ratio == 0.0 || error("SeriesReactanceControl: branch $(fromBus) to $(toBus) is a transformer branch; transformer reactance is owned by the tap machinery (X(alpha)), attach the controller to a line branch instead.")
  _series_ctrl_check_range(br.r_pu, x_min_pu, x_max_pu)
  (x_min_pu <= br.x_pu <= x_max_pu) || error("SeriesReactanceControl: the branch reactance x_pu = $(br.x_pu) starts outside [$(x_min_pu), $(x_max_pu)]; choose a range containing the uncompensated reactance.")
  any(c -> c isa SeriesReactanceControl && c.branch_idx == br.branchIdx, net.machineControls) && error("SeriesReactanceControl: branch $(fromBus) to $(toBus) already has a series-reactance controller.")

  ctrl = SeriesReactanceControl(
    something(name, string("TCSC_", fromBus, "_", toBus)),
    fromBus,
    toBus,
    br.branchIdx,
    p_target_mw,
    x_min_pu,
    x_max_pu,
    br.x_pu,
    deadband_p_mw,
    max_outer_iters,
    enabled,
    :idle,
    false,
    false,
    nothing,
    nothing,
    nothing,
    0,
  )
  push!(net.machineControls, ctrl)
  return ctrl
end

# resolve the controlled branch defensively (same reasoning as the shunt
# controller: indices are push-only stable, but a copied controller may point
# elsewhere; fail loudly instead of actuating the wrong element)
function _find_controlled_series_branch(net::Net, ctrl::SeriesReactanceControl)::Branch
  (1 <= ctrl.branch_idx <= length(net.branchVec)) || error("SeriesReactanceControl $(ctrl.name): branch index $(ctrl.branch_idx) out of range")
  br = net.branchVec[ctrl.branch_idx]
  from = geNetBusIdx(net = net, busName = ctrl.fromBus)
  to = geNetBusIdx(net = net, busName = ctrl.toBus)
  ((br.fromBus == from && br.toBus == to) || (br.fromBus == to && br.toBus == from)) || error("SeriesReactanceControl $(ctrl.name): branch $(ctrl.branch_idx) does not connect $(ctrl.fromBus) and $(ctrl.toBus)")
  return br
end

control_name(ctrl::SeriesReactanceControl) = string(ctrl.name, " TCSC")
control_enabled(ctrl::SeriesReactanceControl) = ctrl.enabled
control_initialize!(::SeriesReactanceControl, ::Net, context) = NoControlState()
control_status(ctrl::SeriesReactanceControl, ::AbstractControlState)::Symbol = ctrl.status
control_is_converged(ctrl::SeriesReactanceControl, ::AbstractControlState)::Bool = ctrl.converged
control_is_blocked(ctrl::SeriesReactanceControl, ::AbstractControlState)::Bool = ctrl.at_limit

function control_evaluate!(ctrl::SeriesReactanceControl, net::Net, ::AbstractControlState, context)
  p = get_branch_p_from_to_mw(net, ctrl.fromBus, ctrl.toBus)
  ctrl.achieved_p_mw = p
  ctrl.converged = p !== nothing && abs(p - ctrl.p_target_mw) <= ctrl.deadband_p_mw
  ctrl.status = ctrl.converged ? :converged : (ctrl.at_limit ? :at_limit : :active)
  ctrl.outer_iters = context.outer_iteration
  return nothing
end

# Secant strategy as in the machine and shunt controllers, with one
# difference: the sign of dP/dx is not hard-coded (in a meshed network it
# depends on which side of the parallel paths the branch sits). The
# bootstrap probes toward the side with more headroom while the sensitivity
# is unmeasured; the secant then uses the measured signed slope, guarded by
# a minimum magnitude.
function control_propose_update!(ctrl::SeriesReactanceControl, net::Net, ::AbstractControlState, context)
  x0 = ctrl.x_pu
  ctrl.converged && return (old_x = x0, new_x = x0)
  p = ctrl.achieved_p_mw
  p === nothing && return (old_x = x0, new_x = x0)
  err = ctrl.p_target_mw - p

  dx = if ctrl.prev_x_pu !== nothing && ctrl.prev_p_mw !== nothing && abs(x0 - ctrl.prev_x_pu) > 1e-9
    slope = (p - ctrl.prev_p_mw) / (x0 - ctrl.prev_x_pu)
    abs(slope) > _MACHINE_CTRL_MIN_SLOPE ? err / slope : _series_ctrl_bootstrap_step(ctrl, x0)
  else
    _series_ctrl_bootstrap_step(ctrl, x0)
  end
  new_x = clamp(x0 + dx, ctrl.x_min_pu, ctrl.x_max_pu)
  return (old_x = x0, new_x = new_x)
end

# direction-agnostic probe: step toward the side with more headroom so a
# measurable (x, P) pair exists for the secant, bounded by the bootstrap
# fraction of that headroom
function _series_ctrl_bootstrap_step(ctrl::SeriesReactanceControl, x0::Float64)::Float64
  up = ctrl.x_max_pu - x0
  down = ctrl.x_min_pu - x0
  headroom = abs(up) >= abs(down) ? up : down
  return _MACHINE_CTRL_BOOTSTRAP_FRACTION * headroom
end

function control_apply_update!(ctrl::SeriesReactanceControl, net::Net, ::AbstractControlState, update::NamedTuple, context)::Bool
  moved = abs(update.new_x - update.old_x) > 1e-9
  if moved
    br = _find_controlled_series_branch(net, ctrl)
    # the actuator: one branch stamp changes, the next outer-loop solve
    # re-stamps the Y-bus (verified mechanism, tap_control X(alpha) path)
    br.x_pu = update.new_x
    ctrl.prev_x_pu = update.old_x
    ctrl.prev_p_mw = ctrl.achieved_p_mw
    ctrl.x_pu = update.new_x
  end
  # clamped reactance = fixed compensated line; honest at_limit only when
  # clamped and not converged, the shunt controller semantics
  ctrl.at_limit = !ctrl.converged && (isapprox(ctrl.x_pu, ctrl.x_min_pu; atol = 1e-9) || isapprox(ctrl.x_pu, ctrl.x_max_pu; atol = 1e-9))
  ctrl.at_limit && (ctrl.status = :at_limit)
  return moved
end

function control_element_descriptor(ctrl::SeriesReactanceControl, net::Net)::Union{Nothing,NamedTuple}
  return (
    name = control_name(ctrl),
    element = string("branch@", ctrl.fromBus, "-", ctrl.toBus),
    device = "TCSC (series compensation)",
    actuator = :series_x_pu,
    actuator_min = ctrl.x_min_pu,
    actuator_max = ctrl.x_max_pu,
    quantity = :branch_active_power,
    target = string(ctrl.fromBus, "->", ctrl.toBus),
    target_value = ctrl.p_target_mw,
    discrete = false,
    enabled = ctrl.enabled,
    status = ctrl.status,
    converged = ctrl.converged,
    at_limit = ctrl.at_limit,
  )
end

function control_report_rows(ctrl::SeriesReactanceControl, net::Net, ::AbstractControlState, context)
  return [(
    controller_name = control_name(ctrl),
    controller_type = "SeriesReactanceControl",
    branch = string(ctrl.fromBus, "->", ctrl.toBus),
    status = ctrl.status,
    converged = ctrl.converged,
    at_limit = ctrl.at_limit,
    achieved_p_mw = ctrl.achieved_p_mw === nothing ? missing : ctrl.achieved_p_mw,
    p_target_mw = ctrl.p_target_mw,
    x_pu = ctrl.x_pu,
    x_min_pu = ctrl.x_min_pu,
    x_max_pu = ctrl.x_max_pu,
  )]
end

"""
    printSeriesReactanceControllerSummary(io::IO, net::Net)
    printSeriesReactanceControllerSummary(net::Net)

Engineering-style summary of the registered series-reactance (TCSC)
controllers, one block per controller: branch and direction, target versus
achieved flow, reactance and range, deadband, and the honest limit flags.
Prints nothing when no series controller is registered; the classical
result print calls it unconditionally, mirroring the machine controller
summary.
"""
function printSeriesReactanceControllerSummary(io::IO, net::Net)
  ctrls = _series_reactance_controllers(net)
  isempty(ctrls) && return
  println(io, "\nSeries Reactance Control Summary (TCSC)")
  println(io, "---------------------------------------")
  for c in ctrls
    println(io, control_name(c), " (branch ", c.fromBus, " -> ", c.toBus, ")")
    println(io, "  target P           : ", @sprintf("%.3f MW", c.p_target_mw))
    println(io, "  achieved P         : ", c.achieved_p_mw === nothing ? "-" : @sprintf("%.3f MW", c.achieved_p_mw))
    println(io, "  series reactance   : ", @sprintf("%.5f pu", c.x_pu))
    println(io, "  reactance range    : ", @sprintf("%.5f .. %.5f pu", c.x_min_pu, c.x_max_pu))
    println(io, "  deadband           : ", @sprintf("%.3f MW", c.deadband_p_mw))
    println(io, "  converged          : ", c.converged)
    println(io, "  at_limit           : ", c.at_limit)
    println(io, "  status             : ", c.status)
    if !c.converged && c.at_limit
      println(io, "  status detail      : target not reached, reactance clamped at the range end (fixed compensated line)")
    end
  end
end
printSeriesReactanceControllerSummary(net::Net) = printSeriesReactanceControllerSummary(stdout, net)

function control_trace_rows(ctrl::SeriesReactanceControl, net::Net, ::AbstractControlState, context)
  return [(
    outer_iteration = context.outer_iteration,
    controller_name = control_name(ctrl),
    controller_type = "SeriesReactanceControl",
    branch = string(ctrl.fromBus, "->", ctrl.toBus),
    status = ctrl.status,
    converged = ctrl.converged,
    at_limit = ctrl.at_limit,
    achieved_p_mw = ctrl.achieved_p_mw === nothing ? missing : ctrl.achieved_p_mw,
    p_target_mw = ctrl.p_target_mw,
    x_pu = ctrl.x_pu,
  )]
end

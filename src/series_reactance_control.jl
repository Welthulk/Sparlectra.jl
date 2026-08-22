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
#          X(alpha) coupling of phase-shifting transformers). Includes the
#          SSSC limit mode (#297 Draft F): the bound is the injectable
#          series voltage, so the admissible reactance deviation scales
#          inversely with the branch current, |x - x_base| <= V_inj,max/|I|,
#          re-evaluated every outer iteration.

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

Two limit modes (`limit_mode`):
- `:reactance_range` (default): fixed `[x_min_pu, x_max_pu]`, the TCSC
  reactance window.
- `:injected_voltage`: SSSC behavior (issue #297 Draft F). The device
  injects a series voltage in quadrature with the line current; stationary
  it acts as a reactance deviation from the natural line reactance
  `x_base_pu`, bounded by the injectable voltage magnitude:
  `|x - x_base_pu| <= v_inj_max_pu / |I|`. The bound is re-evaluated from
  the solved branch current before every outer step, so at high loading the
  usable reactance window SHRINKS (less compensation available exactly when
  the current is large) — the defining contrast to the TCSC's fixed window.
  An at-limit SSSC keeps adjusting while its current-dependent bound still
  moves and only parks `at_limit` once the bound has settled.
"""
mutable struct SeriesReactanceControl <: AbstractOuterController
  name::String
  fromBus::String
  toBus::String
  branch_idx::Int
  p_target_mw::Float64
  x_min_pu::Float64           # live bound in :injected_voltage mode (refreshed per outer iteration)
  x_max_pu::Float64
  limit_mode::Symbol          # :reactance_range (TCSC) or :injected_voltage (SSSC)
  v_inj_max_pu::Union{Nothing,Float64} # injectable series voltage; nothing in :reactance_range mode
  x_base_pu::Float64          # natural (uncompensated) branch reactance at registration
  x_pu::Float64
  deadband_p_mw::Float64
  max_outer_iters::Int
  enabled::Bool
  status::Symbol
  converged::Bool
  at_limit::Bool
  achieved_p_mw::Union{Nothing,Float64}
  i_pu::Union{Nothing,Float64}          # last measured branch current, scales the SSSC bounds
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
  p.u. Negative values (net capacitive branch) are allowed. TCSC mode only.
- `v_inj_max_pu::Union{Nothing,Float64}`: SSSC mode (issue #297 Draft F).
  The maximum injectable series voltage in p.u.; the reactance window
  becomes current-dependent, `|x - x_base| <= v_inj_max_pu / |I|` around
  the natural branch reactance, refreshed every outer iteration. Mutually
  exclusive with `x_min_pu`/`x_max_pu` (the injectable voltage IS the
  limit). Theory in [FACTS Devices](@ref facts_devices).
- `deadband_p_mw::Float64 = 0.5`: convergence band around the target.
- `max_outer_iters::Int = 20`, `enabled::Bool = true`: outer-loop
  budget and switch.
- `name`: optional controller name (defaults to `TCSC_<from>_<to>`, or
  `SSSC_<from>_<to>` in injected-voltage mode).

Fails with an error for a missing branch, a transformer branch (taps own
transformer reactance, see the X(alpha) coupling of the tap controller),
an inverted or missing range (TCSC mode), a doubly specified limit (both a
range and `v_inj_max_pu`), a non-positive injectable voltage, a range whose
series impedance magnitude enters the exclusion guard `eps_z`, a starting
`x_pu` outside the range, or a second series controller on the same branch.
"""
function addSeriesReactanceControl!(
  net::Net;
  fromBus::String,
  toBus::String,
  p_target_mw::Float64,
  x_min_pu::Union{Nothing,Float64} = nothing,
  x_max_pu::Union{Nothing,Float64} = nothing,
  v_inj_max_pu::Union{Nothing,Float64} = nothing,
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
  # limit-mode resolution: an injectable-voltage bound switches the
  # controller into SSSC mode; a fixed window and the voltage bound are
  # mutually exclusive by construction (a device has one or the other).
  local limit_mode::Symbol
  local xlo::Float64
  local xhi::Float64
  if v_inj_max_pu !== nothing
    (x_min_pu !== nothing || x_max_pu !== nothing) && error("SeriesReactanceControl: x_min_pu/x_max_pu and v_inj_max_pu are mutually exclusive — the injectable series voltage is the limit in SSSC mode.")
    v_inj_max_pu > 0.0 || error("SeriesReactanceControl: v_inj_max_pu must be positive, got $(v_inj_max_pu)")
    limit_mode = :injected_voltage
    # no current is known before the first solve: start with the zero-width
    # window at the natural reactance; control_evaluate! widens it from the
    # first solved operating point (propose only runs after evaluate)
    xlo = br.x_pu
    xhi = br.x_pu
  else
    (x_min_pu === nothing || x_max_pu === nothing) && error("SeriesReactanceControl: pass x_min_pu and x_max_pu (TCSC reactance window) or v_inj_max_pu (SSSC injectable-voltage limit).")
    limit_mode = :reactance_range
    _series_ctrl_check_range(br.r_pu, x_min_pu, x_max_pu)
    (x_min_pu <= br.x_pu <= x_max_pu) || error("SeriesReactanceControl: the branch reactance x_pu = $(br.x_pu) starts outside [$(x_min_pu), $(x_max_pu)]; choose a range containing the uncompensated reactance.")
    xlo = x_min_pu
    xhi = x_max_pu
  end
  any(c -> c isa SeriesReactanceControl && c.branch_idx == br.branchIdx, net.machineControls) && error("SeriesReactanceControl: branch $(fromBus) to $(toBus) already has a series-reactance controller.")

  ctrl = SeriesReactanceControl(
    something(name, string(limit_mode === :injected_voltage ? "SSSC_" : "TCSC_", fromBus, "_", toBus)),
    fromBus,
    toBus,
    br.branchIdx,
    p_target_mw,
    xlo,
    xhi,
    limit_mode,
    v_inj_max_pu,
    br.x_pu,
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

control_name(ctrl::SeriesReactanceControl) = string(ctrl.name, ctrl.limit_mode === :injected_voltage ? " SSSC" : " TCSC")
control_enabled(ctrl::SeriesReactanceControl) = ctrl.enabled
control_initialize!(::SeriesReactanceControl, ::Net, context) = NoControlState()
control_status(ctrl::SeriesReactanceControl, ::AbstractControlState)::Symbol = ctrl.status
control_is_converged(ctrl::SeriesReactanceControl, ::AbstractControlState)::Bool = ctrl.converged
control_is_blocked(ctrl::SeriesReactanceControl, ::AbstractControlState)::Bool = ctrl.at_limit

# SSSC live-bounds constants: the current floor keeps the window finite on a
# (numerically) currentless branch — physically an SSSC injects no voltage
# without current, so the model is unconstrained there anyway; the tracking
# tolerance decides when a moving bound has settled (pu reactance).
const _SSSC_MIN_CURRENT_PU = 1e-4
const _SSSC_BOUND_TRACK_TOL_PU = 1e-6

function control_evaluate!(ctrl::SeriesReactanceControl, net::Net, ::AbstractControlState, context)
  p = get_branch_p_from_to_mw(net, ctrl.fromBus, ctrl.toBus)
  ctrl.achieved_p_mw = p
  if ctrl.limit_mode === :injected_voltage && p !== nothing
    # branch current magnitude at the measured (registered from) side:
    # |I_pu| = |S_pu| / |V_pu|
    q = get_branch_q_from_to_mvar(net, ctrl.fromBus, ctrl.toBus)
    vm = get_bus_vm_pu(net, ctrl.fromBus)
    s_pu = sqrt(p^2 + (q === nothing ? 0.0 : q)^2) / net.baseMVA
    i_pu = (isfinite(vm) && vm > 0.0) ? s_pu / vm : s_pu
    ctrl.i_pu = i_pu
    # the SSSC window: |x - x_base| <= v_inj_max / |I|, floored current so a
    # light-load branch gets a wide (not infinite) window
    w = ctrl.v_inj_max_pu / max(i_pu, _SSSC_MIN_CURRENT_PU)
    xlo = ctrl.x_base_pu - w
    xhi = ctrl.x_base_pu + w
    # eps_z guard as a clamp (not an error: a transient current must never
    # abort the run): when the branch resistance alone does not clear the
    # guard, exclude the near-resonance reactance band and keep the side of
    # the window that contains the natural reactance
    if abs(_find_controlled_series_branch(net, ctrl).r_pu) <= _SERIES_CTRL_EPS_Z
      x_guard = _SERIES_CTRL_EPS_Z
      if ctrl.x_base_pu >= 0.0
        xlo = max(xlo, x_guard)
        xhi = max(xhi, xlo)
      else
        xhi = min(xhi, -x_guard)
        xlo = min(xlo, xhi)
      end
    end
    ctrl.x_min_pu = xlo
    ctrl.x_max_pu = xhi
    # release a parked controller whose bound moved: the previous clamp no
    # longer sits on the current limit, so there is tracking work left
    if ctrl.at_limit && min(abs(ctrl.x_pu - xlo), abs(ctrl.x_pu - xhi)) > _SSSC_BOUND_TRACK_TOL_PU
      ctrl.at_limit = false
    end
  end
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
    device = ctrl.limit_mode === :injected_voltage ? "SSSC (VSC)" : "TCSC (series compensation)",
    actuator = :series_x_pu,
    # live bounds in :injected_voltage mode: the window of the LAST evaluated
    # operating point, x_base +- v_inj_max/|I| at the solved branch current
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
    limit_mode = ctrl.limit_mode,
    v_inj_max_pu = ctrl.v_inj_max_pu === nothing ? missing : ctrl.v_inj_max_pu,
    x_base_pu = ctrl.x_base_pu,
    i_pu = ctrl.i_pu === nothing ? missing : ctrl.i_pu,
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
    if c.limit_mode === :injected_voltage
      println(io, "  limit mode         : SSSC injected-voltage limit, V_inj,max = ", @sprintf("%.4f pu", c.v_inj_max_pu))
      println(io, "  live x window      : ", @sprintf("%.5f .. %.5f pu (x_base %.5f, |I| = %s pu)", c.x_min_pu, c.x_max_pu, c.x_base_pu, c.i_pu === nothing ? "-" : @sprintf("%.4f", c.i_pu)))
      if c.i_pu !== nothing
        println(io, "  injected voltage   : ", @sprintf("%.4f pu (of %.4f pu available)", abs(c.x_pu - c.x_base_pu) * c.i_pu, c.v_inj_max_pu))
      end
    else
      println(io, "  reactance range    : ", @sprintf("%.5f .. %.5f pu", c.x_min_pu, c.x_max_pu))
    end
    println(io, "  deadband           : ", @sprintf("%.3f MW", c.deadband_p_mw))
    println(io, "  converged          : ", c.converged)
    println(io, "  at_limit           : ", c.at_limit)
    println(io, "  status             : ", c.status)
    if !c.converged && c.at_limit
      detail = c.limit_mode === :injected_voltage ? "target not reached, injectable series voltage exhausted (window tracks the branch current)" : "target not reached, reactance clamped at the range end (fixed compensated line)"
      println(io, "  status detail      : ", detail)
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

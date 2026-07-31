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

# file: src/shunt_control.jl
# purpose: SVC-style variable-shunt voltage controller (outer loop): a
#          susceptance actuator holds the local bus voltage; at a limit the
#          device keeps its clamped susceptance, so the reactive output then
#          follows V² — the honest SVC limit region, unlike a machine's
#          constant-Q limit.

"""
    ShuntVoltageControl <: AbstractOuterController

SVC-like variable-shunt voltage controller. The actuator is the shunt
susceptance, expressed as MVAr at 1.0 p.u. of the bus nominal voltage
(MATPOWER `Bs` convention, capacitive positive) within
`[bs_min_mvar, bs_max_mvar]`. In range the controller holds
`target_vm_pu` at its own bus via secant iteration; at a limit the
susceptance stays clamped and the injected reactive power follows the bus
voltage squared through the Y-bus stamp — the constant-B region of a real
SVC. Reported honestly via `at_limit`.
"""
mutable struct ShuntVoltageControl <: AbstractOuterController
  name::String
  bus::String
  shunt_idx::Int
  target_vm_pu::Float64
  bs_min_mvar::Float64
  bs_max_mvar::Float64
  bs_mvar::Float64
  deadband_vm_pu::Float64
  max_outer_iters::Int
  enabled::Bool
  status::Symbol
  converged::Bool
  at_limit::Bool
  achieved_vm_pu::Union{Nothing,Float64}
  prev_bs_mvar::Union{Nothing,Float64}
  prev_vm_pu::Union{Nothing,Float64}
  outer_iters::Int
end

"""
    _shunt_controllers(net) -> Vector{ShuntVoltageControl}

Collect the shunt voltage controllers stored on `net`. They share the
generic outer-controller registry (`net.machineControls`) with the machine
controllers; guarded with `hasproperty` because `collect_outer_controllers`
accepts `net::Any`.
"""
function _shunt_controllers(net)::Vector{ShuntVoltageControl}
  hasproperty(net, :machineControls) || return ShuntVoltageControl[]
  return ShuntVoltageControl[c for c in net.machineControls if c isa ShuntVoltageControl]
end

"""
    clearShuntControllers!(net)

Remove all shunt voltage controllers from `net` (machine controllers stay).
"""
function clearShuntControllers!(net::Net)
  filter!(c -> !(c isa ShuntVoltageControl), net.machineControls)
  return net
end

"""
    addShuntVoltageControl!(net; bus, target_vm_pu, bs_min_mvar, bs_max_mvar, ...)

Add an SVC-style variable-shunt voltage controller at `bus`. Creates its own
shunt element (initially at `bs_start_mvar`) whose susceptance the outer
control loop moves within `[bs_min_mvar, bs_max_mvar]` to hold
`target_vm_pu` at the bus.

# Arguments
- `net::Net`: the network.
- `bus::String`: the supported bus. Must be a PQ bus — a PV or slack bus is
  already voltage-held by another unit.
- `target_vm_pu::Float64`: voltage target at `bus` in p.u.
- `bs_min_mvar`, `bs_max_mvar`: susceptance range as MVAr at 1.0 p.u.
  (inductive negative, capacitive positive).
- `bs_start_mvar::Float64 = 0.0`: initial susceptance, clamped to the range.
- `deadband_vm_pu::Float64 = 1e-3`: convergence band around the target.
- `max_outer_iters::Int = 20`, `enabled::Bool = true`: outer-loop budget/switch.
- `name`: optional controller name (defaults to `SVC_<bus>`).

Fails with an error for a missing bus, an inverted susceptance range, a
non-PQ bus, or a second shunt controller on the same bus. Warns when a
transformer tap controller already regulates the same bus voltage (two
controllers steering one voltage fight each other).
"""
function addShuntVoltageControl!(
  net::Net;
  bus::String,
  target_vm_pu::Float64,
  bs_min_mvar::Float64,
  bs_max_mvar::Float64,
  bs_start_mvar::Float64 = 0.0,
  deadband_vm_pu::Float64 = 1e-3,
  max_outer_iters::Int = 20,
  enabled::Bool = true,
  name::Union{Nothing,String} = nothing,
)
  busIdx = geNetBusIdx(net = net, busName = bus)
  target_vm_pu > 0.0 || error("ShuntVoltageControl: target_vm_pu must be positive, got $(target_vm_pu)")
  bs_min_mvar < bs_max_mvar || error("ShuntVoltageControl: bs_min_mvar ($(bs_min_mvar)) must be below bs_max_mvar ($(bs_max_mvar))")
  node = net.nodeVec[busIdx]
  node._nodeType == PQ || error("ShuntVoltageControl: bus $(bus) is $(node._nodeType) — its voltage is already held; a shunt controller needs a PQ bus.")
  any(c -> c isa ShuntVoltageControl && c.bus == bus, net.machineControls) && error("ShuntVoltageControl: bus $(bus) already has a shunt voltage controller.")
  # Cross-type check, same reasoning as for machine controllers: a
  # tap-regulated bus stays PQ, so only this warning can flag the double
  # regulation.
  for tf in net.trafos
    for w in (tf.side1, tf.side2, tf.side3)
      w === nothing && continue
      for ctrl in w.controls
        if ctrl.target_bus == bus && ctrl.target_vm_pu !== nothing
          @warn "ShuntVoltageControl: a transformer tap controller already regulates the voltage of bus $(bus) — two controllers on one voltage fight each other; expect oscillating outer iterations or reconfigure one of them."
        end
      end
    end
  end

  bs0 = clamp(bs_start_mvar, bs_min_mvar, bs_max_mvar)
  addShuntMatpower!(net = net, busName = bus, Gs = 0.0, Bs = bs0)
  ctrl = ShuntVoltageControl(
    something(name, string("SVC_", bus)),
    bus,
    length(net.shuntVec),
    target_vm_pu,
    bs_min_mvar,
    bs_max_mvar,
    bs0,
    deadband_vm_pu,
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

# resolve the controlled shunt defensively (same reasoning as the machine
# controller: indices are push-only stable, but a copied controller may point
# elsewhere — fail loudly instead of actuating the wrong element)
function _find_controlled_shunt(net::Net, ctrl::ShuntVoltageControl)::Shunt
  (1 <= ctrl.shunt_idx <= length(net.shuntVec)) || error("ShuntVoltageControl $(ctrl.name): shunt index $(ctrl.shunt_idx) out of range")
  sh = net.shuntVec[ctrl.shunt_idx]
  sh.busIdx == geNetBusIdx(net = net, busName = ctrl.bus) || error("ShuntVoltageControl $(ctrl.name): shunt $(ctrl.shunt_idx) is not connected to bus $(ctrl.bus)")
  return sh
end

control_name(ctrl::ShuntVoltageControl) = string(ctrl.name, " SVC")
control_enabled(ctrl::ShuntVoltageControl) = ctrl.enabled
control_initialize!(::ShuntVoltageControl, ::Net, context) = NoControlState()
control_status(ctrl::ShuntVoltageControl, ::AbstractControlState)::Symbol = ctrl.status
control_is_converged(ctrl::ShuntVoltageControl, ::AbstractControlState)::Bool = ctrl.converged
control_is_blocked(ctrl::ShuntVoltageControl, ::AbstractControlState)::Bool = ctrl.at_limit

function control_evaluate!(ctrl::ShuntVoltageControl, net::Net, ::AbstractControlState, context)
  vm = get_bus_vm_pu(net, ctrl.bus)
  ctrl.achieved_vm_pu = vm
  ctrl.converged = _voltage_within_deadband(vm, ctrl.target_vm_pu, ctrl.deadband_vm_pu)
  ctrl.status = ctrl.converged ? :converged : (ctrl.at_limit ? :at_limit : :active)
  ctrl.outer_iters = context.outer_iteration
  return nothing
end

# Same secant strategy as the machine controller: bounded bootstrap step
# while the sensitivity is unmeasured, secant afterwards. dV/dBs is positive
# (a capacitive shunt raises the local voltage), so the slope guard reuses
# the machine controller's threshold.
function control_propose_update!(ctrl::ShuntVoltageControl, net::Net, ::AbstractControlState, context)
  b0 = ctrl.bs_mvar
  ctrl.converged && return (old_bs = b0, new_bs = b0)
  vm = ctrl.achieved_vm_pu
  vm === nothing && return (old_bs = b0, new_bs = b0)
  err = ctrl.target_vm_pu - vm

  Δb = if ctrl.prev_bs_mvar !== nothing && ctrl.prev_vm_pu !== nothing && abs(b0 - ctrl.prev_bs_mvar) > 1e-9
    slope = (vm - ctrl.prev_vm_pu) / (b0 - ctrl.prev_bs_mvar)
    slope > _MACHINE_CTRL_MIN_SLOPE ? err / slope : _shunt_ctrl_bootstrap_step(ctrl, b0, err)
  else
    _shunt_ctrl_bootstrap_step(ctrl, b0, err)
  end
  new_bs = clamp(b0 + Δb, ctrl.bs_min_mvar, ctrl.bs_max_mvar)
  return (old_bs = b0, new_bs = new_bs)
end

function _shunt_ctrl_bootstrap_step(ctrl::ShuntVoltageControl, b0::Float64, err::Float64)::Float64
  headroom = err > 0.0 ? (ctrl.bs_max_mvar - b0) : (ctrl.bs_min_mvar - b0)
  return _MACHINE_CTRL_BOOTSTRAP_FRACTION * headroom
end

function control_apply_update!(ctrl::ShuntVoltageControl, net::Net, ::AbstractControlState, update::NamedTuple, context)::Bool
  moved = abs(update.new_bs - update.old_bs) > 1e-9
  if moved
    sh = _find_controlled_shunt(net, ctrl)
    # MATPOWER stamping convention (addShuntMatpower!): y_pu = (Gs + jBs)/baseMVA
    sh.y_pu_shunt = ComplexF64(0.0, update.new_bs) / net.baseMVA
    sh.G_shunt = real(sh.y_pu_shunt)
    sh.B_shunt = imag(sh.y_pu_shunt)
    sh.p_shunt = 0.0
    sh.q_shunt = update.new_bs
    ctrl.prev_bs_mvar = update.old_bs
    ctrl.prev_vm_pu = ctrl.achieved_vm_pu
    ctrl.bs_mvar = update.new_bs
  end
  # clamped susceptance = constant-B limit region; the injected Q keeps
  # following V² through the Y-bus, which is exactly the SVC device behavior
  ctrl.at_limit = !ctrl.converged && (isapprox(ctrl.bs_mvar, ctrl.bs_min_mvar; atol = 1e-9) || isapprox(ctrl.bs_mvar, ctrl.bs_max_mvar; atol = 1e-9))
  ctrl.at_limit && (ctrl.status = :at_limit)
  return moved
end

function control_element_descriptor(ctrl::ShuntVoltageControl, net::Net)::Union{Nothing,NamedTuple}
  return (
    name = control_name(ctrl),
    element = string("shunt@", ctrl.bus),
    device = "SVC (variable shunt)",
    actuator = :shunt_bs_mvar,
    actuator_min = ctrl.bs_min_mvar,
    actuator_max = ctrl.bs_max_mvar,
    quantity = :bus_voltage,
    target = ctrl.bus,
    target_value = ctrl.target_vm_pu,
    discrete = false,
    enabled = ctrl.enabled,
    status = ctrl.status,
    converged = ctrl.converged,
    at_limit = ctrl.at_limit,
  )
end

function control_report_rows(ctrl::ShuntVoltageControl, net::Net, ::AbstractControlState, context)
  return [(
    controller_name = control_name(ctrl),
    controller_type = "ShuntVoltageControl",
    bus = ctrl.bus,
    status = ctrl.status,
    converged = ctrl.converged,
    at_limit = ctrl.at_limit,
    achieved_vm_pu = ctrl.achieved_vm_pu === nothing ? missing : ctrl.achieved_vm_pu,
    target_vm_pu = ctrl.target_vm_pu,
    bs_mvar = ctrl.bs_mvar,
    bs_min_mvar = ctrl.bs_min_mvar,
    bs_max_mvar = ctrl.bs_max_mvar,
  )]
end

function control_trace_rows(ctrl::ShuntVoltageControl, net::Net, ::AbstractControlState, context)
  return [(
    outer_iteration = context.outer_iteration,
    controller_name = control_name(ctrl),
    controller_type = "ShuntVoltageControl",
    bus = ctrl.bus,
    status = ctrl.status,
    converged = ctrl.converged,
    at_limit = ctrl.at_limit,
    achieved_vm_pu = ctrl.achieved_vm_pu === nothing ? missing : ctrl.achieved_vm_pu,
    target_vm_pu = ctrl.target_vm_pu,
    bs_mvar = ctrl.bs_mvar,
  )]
end

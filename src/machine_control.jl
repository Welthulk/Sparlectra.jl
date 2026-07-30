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

# Author: Udo Schmitz (https://github.com/Welthulk)
# file: src/machine_control.jl
# purpose: outer-loop remote voltage control for machines — a generator held
#          as PQ injection whose reactive output is adjusted between power-flow
#          solves until the voltage at a *different* (remote) bus reaches its
#          target, bounded by the machine's reactive limits.

"""
    MachineVoltageControl <: AbstractOuterController

Remote voltage controller for a single machine (generator prosumer).

The machine stays a PQ injection for the inner Newton-Raphson solve; between
solves the controller moves the machine's reactive output `q_mvar` (clamped to
`[qmin_mvar, qmax_mvar]`) until the voltage magnitude at `target_bus` is
within `deadband_vm_pu` of `target_vm_pu`. This is the outer-loop counterpart
of a PV bus whose regulated node is not the machine's own connection point
(CGMES: a `RegulatingControl` whose terminal sits at a foreign bus).

The update rule is a secant iteration on the scalar map Q ↦ Vm(target_bus):
the first move is a bounded fraction of the remaining reactive headroom in the
physically expected direction (more injected Q raises the voltage), every
following move uses the measured sensitivity of the two previous operating
points. Reaching a reactive limit with the target still outside the deadband
parks the controller `at_limit` — the exact analogue of PV→PQ switching under
Q limits, and honest about what the machine can actually deliver.

Runtime fields (`status`, `converged`, `at_limit`, `achieved_vm_pu`, …) are
owned by `run_control!`; construct instances via [`addMachineVoltageControl!`](@ref).
"""
mutable struct MachineVoltageControl <: AbstractOuterController
  name::String                # display name (machine/equipment name)
  bus::String                 # machine's own bus — where the Q actuator lives
  prosumer_idx::Int           # index into net.prosumpsVec (push-only, stable)
  target_bus::String          # regulated (remote) bus
  target_vm_pu::Float64
  deadband_vm_pu::Float64
  qmin_mvar::Float64
  qmax_mvar::Float64
  max_outer_iters::Int
  enabled::Bool
  # runtime state, owned by the control loop
  status::Symbol
  converged::Bool
  at_limit::Bool
  achieved_vm_pu::Union{Nothing,Float64}
  q_mvar::Float64                       # current reactive setpoint (injection, MVAr)
  prev_q_mvar::Union{Nothing,Float64}   # previous secant point
  prev_vm_pu::Union{Nothing,Float64}
  outer_iters::Int
end

"""
    _machine_controllers(net) -> Vector{MachineVoltageControl}

Collect the machine voltage controllers stored on `net`. Guarded with
`hasproperty` because `collect_outer_controllers` accepts `net::Any`.
"""
function _machine_controllers(net)::Vector{MachineVoltageControl}
  hasproperty(net, :machineControls) || return MachineVoltageControl[]
  return MachineVoltageControl[c for c in net.machineControls if c isa MachineVoltageControl]
end

"""
    clearMachineControllers!(net)

Remove all machine voltage controllers from `net`.
"""
function clearMachineControllers!(net::Net)
  empty!(net.machineControls)
  return net
end

# resolve the controlled prosumer defensively: indices are push-only stable,
# but a controller copied onto a reduced net (island extraction) may point
# past the end or at a different machine — fail loudly instead of actuating
# the wrong unit.
function _find_machine_prosumer(net::Net, ctrl::MachineVoltageControl)::ProSumer
  (1 <= ctrl.prosumer_idx <= length(net.prosumpsVec)) || error("MachineVoltageControl $(ctrl.name): prosumer index $(ctrl.prosumer_idx) out of range")
  ps = net.prosumpsVec[ctrl.prosumer_idx]
  isGenerator(ps) || error("MachineVoltageControl $(ctrl.name): prosumer $(ctrl.prosumer_idx) is not a generator")
  busIdx = geNetBusIdx(net = net, busName = ctrl.bus)
  getPosumerBusIndex(ps) == busIdx || error("MachineVoltageControl $(ctrl.name): prosumer $(ctrl.prosumer_idx) is not connected to bus $(ctrl.bus)")
  return ps
end

"""
    addMachineVoltageControl!(net; bus, target_bus, target_vm_pu, ...)

Add a remote voltage controller for the generator at `bus` that regulates the
voltage magnitude at `target_bus`.

# Arguments
- `net::Net`: the network.
- `bus::String`: the machine's own bus. The generator there must be a plain PQ
  injection (not voltage-regulating) — its reactive output is the actuator.
- `target_bus::String`: the regulated bus. Must be a PQ bus: a PV or slack
  target is already voltage-held by another unit and leaves this controller
  nothing to regulate (rejected with an error).
- `target_vm_pu::Float64`: voltage target at `target_bus` in p.u.
- `qmin_mvar`, `qmax_mvar`: reactive actuator range in MVAr. Default to the
  machine's own `minQ`/`maxQ`; required explicitly when the machine carries no
  scalar limits.
- `deadband_vm_pu::Float64 = 1e-3`: convergence band around the target.
- `prosumer_index::Union{Nothing,Int}`: which prosumer to control when several
  generators sit at `bus`; defaults to the single generator there.
- `max_outer_iters::Int = 20`, `enabled::Bool = true`: outer-loop budget/switch.

Fails with an error for a missing bus, a missing or ambiguous machine, a
voltage-regulating (PV) machine, a non-PQ target bus, an inverted Q range, or
a second active controller on the same machine or target bus.
"""
function addMachineVoltageControl!(
  net::Net;
  bus::String,
  target_bus::String,
  target_vm_pu::Float64,
  qmin_mvar::Union{Nothing,Float64} = nothing,
  qmax_mvar::Union{Nothing,Float64} = nothing,
  deadband_vm_pu::Float64 = 1e-3,
  prosumer_index::Union{Nothing,Int} = nothing,
  name::Union{Nothing,String} = nothing,
  max_outer_iters::Int = 20,
  enabled::Bool = true,
)
  busIdx = geNetBusIdx(net = net, busName = bus)
  targetIdx = geNetBusIdx(net = net, busName = target_bus)
  bus == target_bus && error("MachineVoltageControl: target_bus equals the machine bus $(bus) — local voltage control is a PV setpoint (addProsumer! with vm_pu), not a remote controller.")
  target_vm_pu > 0.0 || error("MachineVoltageControl: target_vm_pu must be positive, got $(target_vm_pu)")

  # The regulated bus must be free to move: a PV/slack target is held by
  # another voltage-controlling unit and the two controllers would fight.
  tnode = net.nodeVec[targetIdx]
  tnode._nodeType == PQ || error("MachineVoltageControl: target bus $(target_bus) is $(tnode._nodeType) — its voltage is already held; a remote controller needs a PQ target.")

  # resolve the machine
  ps_idx = if prosumer_index !== nothing
    (1 <= prosumer_index <= length(net.prosumpsVec)) || error("MachineVoltageControl: prosumer_index $(prosumer_index) out of range")
    ps = net.prosumpsVec[prosumer_index]
    isGenerator(ps) || error("MachineVoltageControl: prosumer $(prosumer_index) is not a generator")
    getPosumerBusIndex(ps) == busIdx || error("MachineVoltageControl: prosumer $(prosumer_index) is not connected to bus $(bus)")
    prosumer_index
  else
    gens = [i for (i, p) in enumerate(net.prosumpsVec) if isGenerator(p) && getPosumerBusIndex(p) == busIdx]
    isempty(gens) && error("MachineVoltageControl: no generator at bus $(bus)")
    length(gens) > 1 && error("MachineVoltageControl: $(length(gens)) generators at bus $(bus) — pass prosumer_index to pick one.")
    gens[1]
  end
  ps = net.prosumpsVec[ps_idx]
  # A voltage-regulating machine is (or makes its bus) PV: its reactive output
  # is a solver outcome, not a settable injection.
  ps.isRegulated && error("MachineVoltageControl: the machine at $(bus) is voltage-regulating (PV) — a remote controller needs a PQ machine.")

  qlo = something(qmin_mvar, ps.minQ === nothing ? nothing : ps.minQ)
  qhi = something(qmax_mvar, ps.maxQ === nothing ? nothing : ps.maxQ)
  (qlo === nothing || qhi === nothing) && error("MachineVoltageControl: machine at $(bus) has no scalar minQ/maxQ — pass qmin_mvar/qmax_mvar explicitly.")
  qlo < qhi || error("MachineVoltageControl: empty reactive range [$(qlo), $(qhi)] MVAr")

  # per-actuator and per-target exclusivity, mirroring the tap-controller rule
  for c in _machine_controllers(net)
    c.enabled || continue
    c.prosumer_idx == ps_idx && error("MachineVoltageControl: machine at $(bus) already has an active controller")
    c.target_bus == target_bus && error("MachineVoltageControl: target bus $(target_bus) is already regulated by the controller of $(c.name)")
  end

  q0 = ps.qVal === nothing ? 0.0 : ps.qVal
  ctrl = MachineVoltageControl(
    something(name, string("machine@", bus)),
    bus,
    ps_idx,
    target_bus,
    target_vm_pu,
    deadband_vm_pu,
    qlo,
    qhi,
    max_outer_iters,
    enabled,
    :active,
    false,
    false,
    nothing,
    clamp(q0, qlo, qhi),
    nothing,
    nothing,
    0,
  )
  push!(net.machineControls, ctrl)
  return net
end

# --- AbstractOuterController protocol ---------------------------------------

control_name(ctrl::MachineVoltageControl) = string(ctrl.name, " RVC")
control_enabled(ctrl::MachineVoltageControl) = ctrl.enabled
control_initialize!(::MachineVoltageControl, ::Net, context) = NoControlState()
control_status(ctrl::MachineVoltageControl, ::AbstractControlState)::Symbol = ctrl.status
control_is_converged(ctrl::MachineVoltageControl, ::AbstractControlState)::Bool = ctrl.converged
control_is_blocked(ctrl::MachineVoltageControl, ::AbstractControlState)::Bool = ctrl.at_limit

function control_evaluate!(ctrl::MachineVoltageControl, net::Net, ::AbstractControlState, context)
  vm = get_bus_vm_pu(net, ctrl.target_bus)
  ctrl.achieved_vm_pu = vm
  ctrl.converged = _voltage_within_deadband(vm, ctrl.target_vm_pu, ctrl.deadband_vm_pu)
  ctrl.status = ctrl.converged ? :converged : (ctrl.at_limit ? :at_limit : :active)
  ctrl.outer_iters = context.outer_iteration
  return nothing
end

# Secant fraction of the remaining Q headroom for the first (bootstrap) move.
# Deliberately below 1.0: the sensitivity is unknown at that point, and a
# short first step only costs one extra outer iteration — the secant update
# extrapolates past it immediately.
const _MACHINE_CTRL_BOOTSTRAP_FRACTION = 0.25
# Below this measured sensitivity (p.u. per MVAr) the secant step is
# meaningless — the target bus is electrically insensitive to this machine.
const _MACHINE_CTRL_MIN_SLOPE = 1e-6

function control_propose_update!(ctrl::MachineVoltageControl, net::Net, ::AbstractControlState, context)
  q0 = ctrl.q_mvar
  ctrl.converged && return (old_q = q0, new_q = q0)
  vm = ctrl.achieved_vm_pu
  vm === nothing && return (old_q = q0, new_q = q0)
  err = ctrl.target_vm_pu - vm

  Δq = if ctrl.prev_q_mvar !== nothing && ctrl.prev_vm_pu !== nothing && abs(q0 - ctrl.prev_q_mvar) > 1e-9
    slope = (vm - ctrl.prev_vm_pu) / (q0 - ctrl.prev_q_mvar)
    # A negative or vanishing measured slope contradicts the physics of a
    # working actuator (injecting Q raises the voltage); fall back to the
    # bootstrap step rather than stepping toward the wrong bound.
    slope > _MACHINE_CTRL_MIN_SLOPE ? err / slope : _machine_ctrl_bootstrap_step(ctrl, q0, err)
  else
    _machine_ctrl_bootstrap_step(ctrl, q0, err)
  end
  new_q = clamp(q0 + Δq, ctrl.qmin_mvar, ctrl.qmax_mvar)
  return (old_q = q0, new_q = new_q)
end

# First move without a measured sensitivity: a bounded fraction of the
# remaining headroom toward the physically expected direction.
function _machine_ctrl_bootstrap_step(ctrl::MachineVoltageControl, q0::Float64, err::Float64)::Float64
  headroom = err > 0.0 ? (ctrl.qmax_mvar - q0) : (ctrl.qmin_mvar - q0)
  return _MACHINE_CTRL_BOOTSTRAP_FRACTION * headroom
end

function control_apply_update!(ctrl::MachineVoltageControl, net::Net, ::AbstractControlState, update::NamedTuple, context)::Bool
  moved = abs(update.new_q - update.old_q) > 1e-9
  if moved
    ps = _find_machine_prosumer(net, ctrl)
    # The solver reads bus-level sums (node._qƩGen), the prosumer keeps the
    # per-machine value — update both by the same delta so they stay coherent
    # with any other injections on the bus.
    Δ = update.new_q - (ps.qVal === nothing ? 0.0 : ps.qVal)
    ps.qVal = update.new_q
    node = net.nodeVec[geNetBusIdx(net = net, busName = ctrl.bus)]
    addGenPower!(node = node, p = nothing, q = Δ)
    ctrl.prev_q_mvar = update.old_q
    ctrl.prev_vm_pu = ctrl.achieved_vm_pu
    ctrl.q_mvar = update.new_q
  end
  ctrl.at_limit = !ctrl.converged && (isapprox(ctrl.q_mvar, ctrl.qmin_mvar; atol = 1e-9) || isapprox(ctrl.q_mvar, ctrl.qmax_mvar; atol = 1e-9))
  ctrl.at_limit && (ctrl.status = :at_limit)
  return moved
end

function control_report_rows(ctrl::MachineVoltageControl, net::Net, ::AbstractControlState, context)
  return buildMachineControllerReportRows(net; only = ctrl)
end

function control_trace_rows(ctrl::MachineVoltageControl, net::Net, ::AbstractControlState, context)
  return [(
    outer_iteration = context.outer_iteration,
    controller_name = control_name(ctrl),
    controller_type = "MachineVoltageControl",
    machine_bus = ctrl.bus,
    target_bus = ctrl.target_bus,
    status = ctrl.status,
    converged = ctrl.converged,
    at_limit = ctrl.at_limit,
    achieved_vm_pu = ctrl.achieved_vm_pu === nothing ? missing : ctrl.achieved_vm_pu,
    target_vm_pu = ctrl.target_vm_pu,
    q_mvar = ctrl.q_mvar,
  )]
end

function _machine_controller_status_label(ctrl::MachineVoltageControl)::String
  !ctrl.enabled && return "inactive"
  ctrl.converged && return "converged"
  ctrl.at_limit && return "at_limit"
  return "active"
end

"""
    buildMachineControllerReportRows(net; only=nothing) -> Vector{NamedTuple}

Typed report rows for machine voltage controllers, one per controller;
`only` restricts the output to a single controller instance.
"""
function buildMachineControllerReportRows(net::Net; only::Union{Nothing,MachineVoltageControl} = nothing)::Vector{NamedTuple}
  rows = NamedTuple[]
  for ctrl in _machine_controllers(net)
    only === nothing || ctrl === only || continue
    achieved = get_bus_vm_pu(net, ctrl.target_bus)
    push!(
      rows,
      (
        controller_name = control_name(ctrl),
        controller_type = "MachineVoltageControl",
        machine_bus = ctrl.bus,
        prosumer_index = ctrl.prosumer_idx,
        target_bus = ctrl.target_bus,
        target_vm_pu = ctrl.target_vm_pu,
        achieved_vm_pu = achieved === nothing ? missing : achieved,
        q_mvar = ctrl.q_mvar,
        qmin_mvar = ctrl.qmin_mvar,
        qmax_mvar = ctrl.qmax_mvar,
        deadband_vm_pu = ctrl.deadband_vm_pu,
        converged = ctrl.converged,
        at_limit = ctrl.at_limit,
        status = _machine_controller_status_label(ctrl),
      ),
    )
  end
  return rows
end

"""
    printMachineControllerSummary(io, net)

Print a compact summary block for all configured machine voltage controllers.
"""
function printMachineControllerSummary(io::IO, net::Net)
  ctrls = _machine_controllers(net)
  isempty(ctrls) && return
  println(io, "\nMachine Voltage Control Summary")
  println(io, "-------------------------------")
  for row in buildMachineControllerReportRows(net)
    println(io, row.controller_name, " (machine at ", row.machine_bus, " -> bus ", row.target_bus, ")")
    println(io, "  target Vm          : ", @sprintf("%.4f pu", row.target_vm_pu))
    println(io, "  achieved Vm        : ", ismissing(row.achieved_vm_pu) ? "-" : @sprintf("%.4f pu", row.achieved_vm_pu))
    println(io, "  reactive output    : ", @sprintf("%.3f MVAr", row.q_mvar))
    println(io, "  reactive range     : ", @sprintf("%.3f .. %.3f MVAr", row.qmin_mvar, row.qmax_mvar))
    println(io, "  deadband           : ", @sprintf("%.4f pu", row.deadband_vm_pu))
    println(io, "  converged          : ", row.converged)
    println(io, "  at_limit           : ", row.at_limit)
    println(io, "  status             : ", row.status)
    if !row.converged && row.at_limit
      println(io, "  status detail      : target not reached because a reactive limit was hit")
    end
  end
end

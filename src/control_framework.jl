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
# Date: 26.5.2026
# file: src/control_framework.jl

abstract type AbstractOuterController end
abstract type AbstractControlState end
abstract type AbstractControlUpdate end

Base.@kwdef struct ControlConfig
  enabled::Bool = true
  max_outer_iterations::Int = 20
  trace::Bool = true
  log_iterations::Bool = true
  stop_on_pf_failure::Bool = true
  controllers::Vector{Any} = Any[]
end

Base.@kwdef struct ControlRunResult
  status::Symbol = :no_controllers
  converged::Bool = false
  outer_iterations::Int = 0
  powerflow_solves::Int = 0
  last_pf_iterations::Int = 0
  last_pf_status::Symbol = :not_run
  controllers::Vector{NamedTuple} = NamedTuple[]
  trace::Vector{NamedTuple} = NamedTuple[]
end

struct NoControlState <: AbstractControlState end
struct NoControlUpdate <: AbstractControlUpdate end

control_name(controller::AbstractOuterController) = string(typeof(controller))
control_enabled(controller::AbstractOuterController) = getfield(controller, :enabled)
control_initialize!(::AbstractOuterController, ::Any, context) = NoControlState()
control_evaluate!(::AbstractOuterController, ::Any, ::AbstractControlState, context) = nothing
control_propose_update!(::AbstractOuterController, ::Any, ::AbstractControlState, context) = NoControlUpdate()
control_apply_update!(::AbstractOuterController, ::Any, ::AbstractControlState, ::AbstractControlUpdate, context)::Bool = false
control_is_converged(::AbstractOuterController, ::AbstractControlState)::Bool = false
control_is_blocked(::AbstractOuterController, ::AbstractControlState)::Bool = false
control_status(::AbstractOuterController, ::AbstractControlState)::Symbol = :active
control_report_rows(::AbstractOuterController, ::Any, ::AbstractControlState)::Vector{NamedTuple} = NamedTuple[]
control_trace_rows(::AbstractOuterController, ::Any, ::AbstractControlState)::Vector{NamedTuple} = NamedTuple[]

collect_outer_controllers(net::Any)::Vector{AbstractOuterController} = AbstractOuterController[_tap_controllers(net)...]

function run_control!(net::Any; controllers::Vector{<:AbstractOuterController} = collect_outer_controllers(net), pf_config = nothing, control_config::ControlConfig = ControlConfig(), verbose::Int = 0, performance_profile = nothing)
  if isempty(controllers)
    return ControlRunResult(status = :no_controllers)
  end
  if !control_config.enabled
    return ControlRunResult(status = :disabled)
  end
  if all(c -> !control_enabled(c), controllers)
    return ControlRunResult(status = :disabled)
  end
  pf_runner = () -> runpf!(net; config = pf_config, verbose = verbose, performance_profile = performance_profile)
  states = AbstractControlState[control_initialize!(ctrl, net, (; pf_config, control_config, verbose, performance_profile)) for ctrl in controllers]
  ite, erg = pf_runner()
  solves = 1
  if erg != 0
    return ControlRunResult(status = :pf_failed, converged = false, powerflow_solves = solves, last_pf_iterations = ite, last_pf_status = :failed)
  end
  calcNetLosses!(net)
  calcLinkFlowsKCL!(net)

  active_ids = findall(control_enabled, controllers)
  isempty(active_ids) && return ControlRunResult(status = :no_active_controllers, converged = true, powerflow_solves = solves, last_pf_iterations = ite, last_pf_status = :ok)
  max_outer = min(control_config.max_outer_iterations, isnothing(pf_config) ? control_config.max_outer_iterations : pf_config.max_iter)
  status = :max_outer_iterations
  trace = NamedTuple[]
  outer_iterations = 0
  for it in 1:max_outer
    outer_iterations = it
    for i in active_ids
      control_evaluate!(controllers[i], net, states[i], (; outer_iteration = it))
    end
    all_done = all(i -> control_is_converged(controllers[i], states[i]) || control_is_blocked(controllers[i], states[i]), active_ids)
    if all_done
      status = :converged
      break
    end
    moved = false
    for i in active_ids
      upd = control_propose_update!(controllers[i], net, states[i], (; outer_iteration = it))
      moved = control_apply_update!(controllers[i], net, states[i], upd, (; outer_iteration = it)) || moved
      if control_config.trace
        append!(trace, control_trace_rows(controllers[i], net, states[i]))
      end
    end
    if !moved
      status = :blocked
      break
    end
    ite, erg = pf_runner()
    solves += 1
    erg != 0 && return ControlRunResult(status = :pf_failed, converged = false, outer_iterations = it, powerflow_solves = solves, last_pf_iterations = ite, last_pf_status = :failed, trace = trace)
    calcNetLosses!(net)
    calcLinkFlowsKCL!(net)
  end
  rows = NamedTuple[]
  for (i, ctrl) in enumerate(controllers)
    append!(rows, control_report_rows(ctrl, net, states[i]))
  end
  return ControlRunResult(status = status, converged = status == :converged, outer_iterations = outer_iterations, powerflow_solves = solves, last_pf_iterations = ite, last_pf_status = :ok, controllers = rows, trace = trace)
end

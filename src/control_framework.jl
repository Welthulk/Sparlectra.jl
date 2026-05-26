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
control_report_rows(::AbstractOuterController, ::Any, ::AbstractControlState, context)::Vector{NamedTuple} = NamedTuple[]
control_trace_rows(::AbstractOuterController, ::Any, ::AbstractControlState, context)::Vector{NamedTuple} = NamedTuple[]
collect_outer_controllers(net::Any)::Vector{AbstractOuterController} = AbstractOuterController[_tap_controllers(net)...]
function control_max_outer_iterations(ctrl::AbstractOuterController)::Int
  hasproperty(ctrl, :max_outer_iters) || return typemax(Int)
  return Int(getproperty(ctrl, :max_outer_iters))
end
function latest_control_result(net::Any)
  hasproperty(net, :control_result) || return nothing
  ref = getproperty(net, :control_result)
  return ref isa Base.RefValue ? ref[] : ref
end

"""
    run_control!(net; ...)

Outer-loop orchestration layer around `runpf!` for controller-based network updates.
It does not replace `runpf!`; it coordinates built-in or user-defined controllers and
returns a `ControlRunResult`.
"""
function run_control!(net::Any; controllers::Vector{<:AbstractOuterController} = collect_outer_controllers(net), pf_config = nothing, control_config::ControlConfig = ControlConfig(), verbose::Int = 0, performance_profile = nothing)
  net.control_result[] = nothing
  if isempty(controllers)
    result = ControlRunResult(status = :no_controllers)
    net.control_result[] = result
    return result
  end
  if !control_config.enabled
    result = ControlRunResult(status = :disabled)
    net.control_result[] = result
    return result
  end
  if all(c -> !control_enabled(c), controllers)
    result = ControlRunResult(status = :disabled)
    net.control_result[] = result
    return result
  end
  pf_runner = () -> runpf!(net; config = pf_config, verbose = verbose, performance_profile = performance_profile)
  context = (pf_config = pf_config, control_config = control_config, verbose = verbose, performance_profile = performance_profile, outer_iteration = 0)
  states = AbstractControlState[control_initialize!(ctrl, net, context) for ctrl in controllers]
  ite, erg = pf_runner()
  solves = 1
  if erg != 0
    result = ControlRunResult(status = :pf_failed, converged = false, powerflow_solves = solves, last_pf_iterations = ite, last_pf_status = :failed)
    net.control_result[] = result
    return result
  end
  calcNetLosses!(net)
  calcLinkFlowsKCL!(net)

  active_ids = findall(control_enabled, controllers)
  if isempty(active_ids)
    result = ControlRunResult(status = :no_active_controllers, converged = true, powerflow_solves = solves, last_pf_iterations = ite, last_pf_status = :ok)
    net.control_result[] = result
    return result
  end
  # Outer-loop budget is controlled only by ControlConfig and controller-specific
  # limits. Inner Newton/PF iteration limits stay inside pf_config/runpf!.
  controller_limit = minimum(control_max_outer_iterations(controllers[i]) for i in active_ids)
  max_outer = min(control_config.max_outer_iterations, controller_limit)
  status = :max_outer_iterations
  trace = NamedTuple[]
  outer_iterations = 0
  for it in 1:max_outer
    outer_iterations = it
    for i in active_ids
      ctx = (pf_config = pf_config, control_config = control_config, verbose = verbose, performance_profile = performance_profile, outer_iteration = it)
      control_evaluate!(controllers[i], net, states[i], ctx)
    end
    all_done = all(i -> control_is_converged(controllers[i], states[i]) || control_is_blocked(controllers[i], states[i]), active_ids)
    if all_done
      status = :converged
      break
    end
    moved = false
    for i in active_ids
      ctx = (pf_config = pf_config, control_config = control_config, verbose = verbose, performance_profile = performance_profile, outer_iteration = it)
      upd = control_propose_update!(controllers[i], net, states[i], ctx)
      moved = control_apply_update!(controllers[i], net, states[i], upd, ctx) || moved
      if control_config.trace
        append!(trace, control_trace_rows(controllers[i], net, states[i], ctx))
      end
    end
    if !moved
      status = :blocked
      break
    end
    ite, erg = pf_runner()
    solves += 1
    if erg != 0
      result = ControlRunResult(status = :pf_failed, converged = false, outer_iterations = it, powerflow_solves = solves, last_pf_iterations = ite, last_pf_status = :failed, trace = trace)
      net.control_result[] = result
      return result
    end
    calcNetLosses!(net)
    calcLinkFlowsKCL!(net)
  end
  rows = NamedTuple[]
  for (i, ctrl) in enumerate(controllers)
    append!(rows, control_report_rows(ctrl, net, states[i], context))
  end
  result = ControlRunResult(status = status, converged = status == :converged, outer_iterations = outer_iterations, powerflow_solves = solves, last_pf_iterations = ite, last_pf_status = :ok, controllers = rows, trace = trace)
  net.control_result[] = result
  return result
end

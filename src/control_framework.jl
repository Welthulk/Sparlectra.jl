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
  ite, erg = runpf!(net; config = pf_config, verbose = verbose)
  solves = 1
  if erg != 0
    return ControlRunResult(status = :pf_failed, converged = false, powerflow_solves = solves, last_pf_iterations = ite, last_pf_status = :failed)
  end
  calcNetLosses!(net)
  calcLinkFlowsKCL!(net)

  iterations, erg2 = run_tap_controllers_outer!(net; max_ite = isnothing(pf_config) ? 30 : pf_config.max_ite, tol = isnothing(pf_config) ? 1e-6 : pf_config.tol, verbose = verbose, method = isnothing(pf_config) ? :rectangular : pf_config.method)
  status = erg2 == 0 ? :converged : :pf_failed
  rows = NamedTuple[]
  trace = NamedTuple[]
  for ctrl in controllers
    st = getfield(ctrl, :status)
    push!(rows, (name = control_name(ctrl), status = st, enabled = control_enabled(ctrl), converged = getfield(ctrl, :converged), at_limit = getfield(ctrl, :at_limit), outer_iters = getfield(ctrl, :outer_iters)))
    append!(trace, control_trace_rows(ctrl, net, NoControlState()))
  end
  return ControlRunResult(status = status, converged = status == :converged, outer_iterations = iterations, powerflow_solves = solves + iterations, last_pf_iterations = ite, last_pf_status = status == :converged ? :ok : :failed, controllers = rows, trace = trace)
end

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

"""
    SparlectraRunResult

Typed result returned by [`run_sparlectra`](@ref). The solved network is always
available as `net`, while the remaining fields provide stable framework-level
solver, validation, control-loop, timing, and rectangular diagnostic metadata.
"""
struct SparlectraRunResult
  net::Net
  outcome::Symbol
  numerical_converged::Bool
  solution_available::Bool
  limit_validation_status::Symbol
  final_converged::Bool
  reason::Symbol
  reason_text::String
  iterations::Int
  elapsed_s::Float64
  solver_elapsed_s::Union{Nothing,Float64}
  final_mismatch::Float64
  method::Symbol
  control_status::Symbol
  performance_profile::Any
  diagnostics::NamedTuple
end

_rect_status_get(rect_status, name::Symbol, default) = rect_status !== nothing && hasproperty(rect_status, name) ? getproperty(rect_status, name) : default

function _rectangular_status_diagnostics(rect_status)::NamedTuple
  return (
    nr_converged = _rect_status_get(rect_status, :nr_converged, false),
    active_set_converged = _rect_status_get(rect_status, :active_set_converged, false),
    q_limit_active_set_ok = _rect_status_get(rect_status, :q_limit_active_set_ok, false),
    pv_q_limit_violations = _rect_status_get(rect_status, :pv_q_limit_violations, 0),
    ref_q_limit_violations = _rect_status_get(rect_status, :ref_q_limit_violations, 0),
    final_pv_voltage_residual = _rect_status_get(rect_status, :final_pv_voltage_residual, NaN),
    pv_pq_switching_events = _rect_status_get(rect_status, :pv_pq_switching_events, 0),
    qlimit_active_set_changes = _rect_status_get(rect_status, :qlimit_active_set_changes, 0),
    qlimit_reenable_events = _rect_status_get(rect_status, :qlimit_reenable_events, 0),
    oscillating_buses = _rect_status_get(rect_status, :oscillating_buses, 0),
    guarded_narrow_q_pv_buses = _rect_status_get(rect_status, :guarded_narrow_q_pv_buses, 0),
    branch_quality_status = _rect_status_get(rect_status, :branch_quality_status, :not_checked),
    branch_quality_reason = _rect_status_get(rect_status, :branch_quality_reason, :not_checked),
    branch_quality_metrics = _rect_status_get(rect_status, :branch_quality_metrics, nothing),
    wrong_branch_detection = _rect_status_get(rect_status, :wrong_branch_detection, :off),
    wrong_branch_status = _rect_status_get(rect_status, :wrong_branch_status, :not_checked),
    wrong_branch_reason = _rect_status_get(rect_status, :wrong_branch_reason, :not_checked),
    wrong_branch_low_vm_count = _rect_status_get(rect_status, :wrong_branch_low_vm_count, 0),
    wrong_branch_high_vm_count = _rect_status_get(rect_status, :wrong_branch_high_vm_count, 0),
    wrong_branch_angle_spread_deg = _rect_status_get(rect_status, :wrong_branch_angle_spread_deg, NaN),
    wrong_branch_max_branch_angle_deg = _rect_status_get(rect_status, :wrong_branch_max_branch_angle_deg, NaN),
    wrong_branch_branch_angle_violation_count = _rect_status_get(rect_status, :wrong_branch_branch_angle_violation_count, 0),
    wrong_branch_worst_branch_angle_deg = _rect_status_get(rect_status, :wrong_branch_worst_branch_angle_deg, NaN),
    wrong_branch_rescue_attempted = _rect_status_get(rect_status, :wrong_branch_rescue_attempted, false),
    wrong_branch_rescue_reason = _rect_status_get(rect_status, :wrong_branch_rescue_reason, :disabled),
    wrong_branch_rescue_used = _rect_status_get(rect_status, :wrong_branch_rescue_used, false),
    wrong_branch_rescue_attempts = _rect_status_get(rect_status, :wrong_branch_rescue_attempts, 0),
    wrong_branch_rescue_profile = _rect_status_get(rect_status, :wrong_branch_rescue_profile, :none),
    current_iteration_enabled = _rect_status_get(rect_status, :current_iteration_enabled, false),
    current_iteration_attempted = _rect_status_get(rect_status, :current_iteration_attempted, false),
    current_iteration_accepted = _rect_status_get(rect_status, :current_iteration_accepted, false),
    current_iteration_iterations = _rect_status_get(rect_status, :current_iteration_iterations, 0),
    current_iteration_initial_mismatch = _rect_status_get(rect_status, :current_iteration_initial_mismatch, NaN),
    current_iteration_final_mismatch = _rect_status_get(rect_status, :current_iteration_final_mismatch, NaN),
    current_iteration_reason = _rect_status_get(rect_status, :current_iteration_reason, :disabled),
    current_iteration_artifact = _rect_status_get(rect_status, :current_iteration_artifact, ""),
    island_wise_all_converged = _rect_status_get(rect_status, :island_wise_all_converged, false),
    post_merge_validation_status = _rect_status_get(rect_status, :post_merge_validation_status, :not_applicable),
    post_merge_final_mismatch = _rect_status_get(rect_status, :post_merge_final_mismatch, NaN),
    post_merge_mismatch_status = _rect_status_get(rect_status, :post_merge_mismatch_status, :not_applicable),
  )
end

function _profile_start_projection_diagnostics(performance_profile)::NamedTuple
  performance_profile isa AbstractDict || return NamedTuple()
  haskey(performance_profile, :start_projection_summary) || return NamedTuple()
  summary = performance_profile[:start_projection_summary]
  summary isa NamedTuple || return NamedTuple()
  return summary
end

function _rectangular_solution_available(rect_status)::Bool
  rect_status === nothing && return false
  Bool(rect_status.numerical_converged) || return false

  outcome = Symbol(rect_status.status)
  reason = Symbol(rect_status.reason)
  outcome in (:wrong_branch_detected, :angle_spread_exceeded, :branch_angle_exceeded, :wrong_branch_rescue_not_implemented) && return false
  reason in (:wrong_branch_detected, :angle_spread_exceeded, :branch_angle_exceeded, :wrong_branch_rescue_not_implemented, :rescue_requested_but_not_available) && return false
  return true
end

function _rectangular_run_status(rect_status)
  status = rect_status === nothing ? :solver_error : Symbol(rect_status.status)
  outcome = status
  numerical_converged = rect_status !== nothing && Bool(rect_status.numerical_converged)
  return (
    outcome = outcome,
    numerical_converged = numerical_converged,
    solution_available = _rectangular_solution_available(rect_status),
    limit_validation_status = rect_status === nothing || !numerical_converged ? :skip : Bool(rect_status.q_limit_active_set_ok) ? :ok : :fail,
    final_converged = rect_status !== nothing && Bool(rect_status.final_converged),
    reason = rect_status === nothing ? :solver_error : Symbol(rect_status.reason),
    reason_text = rect_status === nothing ? "solver status unavailable" : String(rect_status.reason_text),
    final_mismatch = rect_status === nothing ? NaN : Float64(rect_status.final_mismatch),
  )
end

_is_accepted_control_status(status::Symbol)::Bool = status in (:none, :converged, :disabled, :no_active_controllers, :no_controllers)

function _compose_framework_status(base_status, control_status::Symbol)
  _is_accepted_control_status(control_status) && return base_status
  if control_status === :pf_failed
    meaningful_rejection = base_status.outcome !== :converged || base_status.reason !== :none || !base_status.final_converged
    meaningful_rejection && base_status.numerical_converged && return merge(base_status, (final_converged = false,))
    return merge(base_status, (outcome = :pf_failed, final_converged = false, reason = :pf_failed, reason_text = "power flow failed during controlled framework run"))
  end
  if !base_status.numerical_converged
    return merge(base_status, (outcome = :pf_failed, final_converged = false, reason = :pf_failed, reason_text = "power flow failed during controlled framework run"))
  end
  outcome = control_status === :blocked ? :control_blocked : control_status === :max_outer_iterations ? :control_max_outer_iterations : :control_not_converged
  return merge(base_status, (outcome = outcome, final_converged = false, reason = outcome, reason_text = "control loop ended with status $(control_status)"))
end

function _build_sparlectra_result(net::Net, cfg::SparlectraConfig, execution, performance_profile)::SparlectraRunResult
  rect_status = cfg.powerflow.method === :rectangular ? rectangular_pf_status(net) : nothing
  status = if cfg.powerflow.method === :rectangular
    _rectangular_run_status(rect_status)
  else
    converged = execution.erg == 0
    (outcome = converged ? :converged : :not_converged, numerical_converged = converged, solution_available = converged, limit_validation_status = :skip, final_converged = converged, reason = converged ? :none : :nr_mismatch_not_converged, reason_text = converged ? "none" : "NR mismatch did not converge", final_mismatch = NaN)
  end
  status = _compose_framework_status(status, execution.control_status)
  status = _append_island_failure_message(status, performance_profile)
  diagnostics = cfg.powerflow.method === :rectangular ? merge(_rectangular_status_diagnostics(rect_status), _profile_start_projection_diagnostics(performance_profile)) : NamedTuple()
  return SparlectraRunResult(net, status.outcome, status.numerical_converged, status.solution_available, status.limit_validation_status, status.final_converged, status.reason, status.reason_text, execution.iterations, execution.elapsed_s, execution.solver_elapsed_s, status.final_mismatch, cfg.powerflow.method, execution.control_status, performance_profile, diagnostics)
end

function _sparlectra_status_row(result::SparlectraRunResult)
  base = (converged = result.final_converged, erg = result.final_converged ? 0 : 1, numerical_converged = result.numerical_converged, iterations = result.iterations, elapsed_s = result.elapsed_s, solver_elapsed_s = result.solver_elapsed_s, method = result.method, outcome = result.outcome, control_status = result.control_status, numerical_solution = result.numerical_converged ? "OK" : "FAIL", solution_available = result.solution_available, limit_validation_status = result.limit_validation_status, final_converged = result.final_converged, reason = result.reason, reason_text = result.reason_text, final_mismatch = result.final_mismatch)
  return merge(result.diagnostics, base)
end

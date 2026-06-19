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
#
# This file is included inside module Sparlectra. Do not add a module wrapper here.
#
# Finalization helpers translate solver internals into stable result metadata.
# They must not influence the Newton iteration or active-set decisions.
# Rectangular power-flow final status and diagnostic helpers.
# Date: 29.5.2026
# file: src/powerflow_rectangular/rectangular_final_status.jl

"""
Finalize Q-limit diagnostics after the numerical rectangular solve.

# Why
Numerical NR convergence alone is not sufficient when PV Q-limits are active.
This helper applies the post-solve acceptance policy for remaining PV violations:
- reject when violations remain and bounded acceptance is disabled/exceeded,
- accept with a dedicated reason when bounded acceptance is enabled.

# Returns
NamedTuple with:
- `qlimit_summary`
- `converged`
- `rejection_reason`
"""
function _finalize_rectangular_qlimit_summary(
  net,
  V,
  Sbus_pu,
  bus_types,
  qmin_pu,
  qmax_pu,
  converged,
  rejection_reason;
  verbose::Int,
  qlimit_trace_enabled::Bool,
  q_hyst_pu::Float64,
  tol::Float64,
  pv_table_rows::Int,
  qlimit_guard_accept_bounded_violations::Bool,
  qlimit_guard_max_remaining_violations::Int,
)
  final_Qload_pu = build_qload_pu(net)
  # Print summary only in verbose/trace mode; stay silent otherwise.
  qlimit_summary_io = (verbose > 0 || qlimit_trace_enabled) ? stdout : devnull

  qlimit_summary = _print_rectangular_qlimit_summary(qlimit_summary_io, net, V, Sbus_pu, bus_types, qmin_pu, qmax_pu, final_Qload_pu; q_hyst_pu = q_hyst_pu, tolerance_pu = tol, max_rows = pv_table_rows, max_console_rows = pv_table_rows)

  remaining_pv_violations = qlimit_summary.pv_violations
  bounded_ok = qlimit_guard_accept_bounded_violations && remaining_pv_violations <= qlimit_guard_max_remaining_violations

  converged_ = converged
  rejection_reason_ = rejection_reason

  if remaining_pv_violations > 0 && !bounded_ok
    # Hard rejection: active PV-limit violations remain after the solve.
    verbose > 0 && @warn "Rectangular NR active-set failed because active PV Q-limit violations remain after the numerical solve." pv_violations = qlimit_summary.pv_violations ref_violations = qlimit_summary.ref_violations
    converged_ = false
    rejection_reason_ = :remaining_pv_q_limit_violations
  elseif remaining_pv_violations > 0 && bounded_ok
    # Soft acceptance: bounded residual violations accepted by policy.
    rejection_reason_ = :bounded_q_limit_violations_accepted
  end

  return (qlimit_summary = qlimit_summary, converged = converged_, rejection_reason = rejection_reason_)
end

"""
Finalize wrong-branch diagnostics and map them to convergence policy.

# Why
Wrong-branch detection is evaluated after a numerically converged state exists.
This helper centralizes policy mapping (`:off`, `:fail`, `:rescue`) so callers
receive a normalized result and consistent rejection reason.

# Returns
NamedTuple with:
- `branch_quality`
- `converged`
- `rejection_reason`
- `wrong_branch_rescue_attempted`
- `wrong_branch_rescue_reason`
"""
function _finalize_rectangular_wrong_branch_diagnostics(
  V,
  bus_types,
  Vset,
  slack_idx,
  converged,
  rejection_reason;
  wrong_branch_detection::Symbol,
  wrong_branch_min_vm_pu::Float64,
  wrong_branch_max_vm_pu::Float64,
  wrong_branch_max_angle_spread_deg::Float64,
  wrong_branch_max_branch_angle_deg::Float64,
  wrong_branch_min_low_vm_count::Int,
  net,
)
  branch_quality = _wrong_branch_not_checked_result()
  wrong_branch_rescue_attempted = false
  wrong_branch_rescue_reason = :disabled

  if wrong_branch_detection != :off
    branch_quality = _check_wrong_branch_solution(
      V,
      bus_types,
      Vset,
      slack_idx;
      min_vm_pu = wrong_branch_min_vm_pu,
      max_vm_pu = wrong_branch_max_vm_pu,
      max_angle_spread_deg = wrong_branch_max_angle_spread_deg,
      max_branch_angle_deg = wrong_branch_max_branch_angle_deg,
      min_low_vm_count = wrong_branch_min_low_vm_count,
      net = net,
    )

    if wrong_branch_detection == :fail && branch_quality.status == :warn
      # Escalate suspicious state to hard failure when policy requires it.
      branch_quality = (; branch_quality..., status = :fail, reason = :wrong_branch_detected)
      converged = false
      rejection_reason = :wrong_branch_detected
    elseif wrong_branch_detection == :rescue && branch_quality.status == :warn
      # Rescue mode is currently not implemented in this path.
      branch_quality = (; branch_quality..., status = :wrong_branch_rescue_not_implemented, reason = :rescue_requested_but_not_available)
      converged = false
      rejection_reason = :wrong_branch_rescue_not_implemented
      wrong_branch_rescue_attempted = false
      wrong_branch_rescue_reason = :rescue_requested_but_not_available
    end
  end

  return (branch_quality = branch_quality, converged = converged, rejection_reason = rejection_reason, wrong_branch_rescue_attempted = wrong_branch_rescue_attempted, wrong_branch_rescue_reason = wrong_branch_rescue_reason)
end

"""
Build the final rectangular solver status payload.

# Why
A stable status schema is required for reporting, tests, and UI consumers.
This helper consolidates convergence flags, Q-limit statistics, and wrong-branch
metrics into one normalized NamedTuple.

# Returns
NamedTuple with:
- `final_reason`
- `final_status`
- `status` (full diagnostics payload)
"""
function _build_rectangular_final_status(
  net,
  numerical_converged::Bool,
  q_limit_active_set_ok::Bool,
  converged::Bool,
  rejection_reason::Symbol,
  qlimit_summary,
  final_pv_voltage_residual::Float64,
  history,
  qlimit_active_set_changes::Int,
  qlimit_reenable_events::Int,
  oscillating_buses::Int,
  guarded_qlimit_buses,
  branch_quality,
  wrong_branch_detection::Symbol,
  wrong_branch_rescue_attempted::Bool,
  wrong_branch_rescue_reason::Symbol,
  mismatch_diagnostics = NamedTuple(),
)
  final_reason = converged ? :none : rejection_reason

  # Reclassify plain mismatch failures when active-set churn is clearly visible.
  if final_reason == :nr_mismatch_not_converged && qlimit_active_set_changes >= 3
    final_reason = :nr_mismatch_not_converged_active_set_unstable
  end

  final_status = _rectangular_solver_status_symbol(numerical_converged, q_limit_active_set_ok, converged, final_reason)

  status = (
    numerical_converged = numerical_converged,
    nr_converged = numerical_converged,
    active_set_converged = q_limit_active_set_ok,
    q_limit_active_set_ok = q_limit_active_set_ok,
    final_converged = converged,
    status = final_status,
    reason = final_reason,
    reason_text = _rectangular_rejection_reason_text(final_reason),
    pv_q_limit_violations = isnothing(qlimit_summary) ? 0 : qlimit_summary.pv_violations,
    ref_q_limit_violations = isnothing(qlimit_summary) ? 0 : qlimit_summary.ref_violations,
    final_pv_voltage_residual = final_pv_voltage_residual,
    final_mismatch = isempty(history) ? Inf : history[end],
    initial_mismatch = isempty(history) ? NaN : history[1],
    nr_initial_mismatch = isempty(history) ? NaN : history[1],
    nr_final_mismatch = isempty(history) ? Inf : history[end],
    mismatch_diagnostics...,
    pv_pq_switching_events = length(net.qLimitLog),
    qlimit_active_set_changes = qlimit_active_set_changes,
    qlimit_reenable_events = qlimit_reenable_events,
    oscillating_buses = oscillating_buses,
    guarded_narrow_q_pv_buses = length(guarded_qlimit_buses),
    branch_quality_status = branch_quality.status,
    branch_quality_reason = branch_quality.reason,
    branch_quality_metrics = branch_quality,
    wrong_branch_detection = wrong_branch_detection,
    wrong_branch_status = branch_quality.status,
    wrong_branch_reason = branch_quality.reason,
    wrong_branch_low_vm_count = branch_quality.low_vm_count,
    wrong_branch_high_vm_count = branch_quality.high_vm_count,
    wrong_branch_angle_spread_deg = branch_quality.angle_spread_deg,
    wrong_branch_max_branch_angle_deg = branch_quality.max_branch_angle_deg,
    wrong_branch_branch_angle_violation_count = branch_quality.branch_angle_violation_count,
    wrong_branch_worst_branch_angle_deg = branch_quality.worst_branch_angle_deg,
    wrong_branch_rescue_attempted = wrong_branch_rescue_attempted,
    wrong_branch_rescue_reason = wrong_branch_rescue_reason,
    # Reserved fields keep downstream status consumers schema-stable.
    wrong_branch_rescue_used = false,
    wrong_branch_rescue_attempts = 0,
    wrong_branch_rescue_profile = :none,
  )

  return (final_reason = final_reason, final_status = final_status, status = status)
end

"""
Store final rectangular solver status and optionally print convergence summary.

# Why
Status persistence and user-visible summary should be emitted from one place to
avoid divergence between stored diagnostics and console output.
"""
function _store_and_print_rectangular_final_status!(net, status, verbose::Int)
  stored_status = _set_rectangular_pf_status!(net, status)
  if verbose > 0
    _print_rectangular_convergence_summary(stdout, stored_status)
  end
  return stored_status
end

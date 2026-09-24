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
# purpose: post-solve status finalization: Q-limit summary acceptance policy,
#          wrong-branch diagnostics, final status NamedTuple assembly, and
#          status storage plus printing

"""
    _final_q_limit_class(dev_pu; tol, q_hyst_pu, final_q_accept_pu) -> Symbol

The class of one overshoot `dev_pu` (`Q - Qmax` or `Qmin - Q`, in pu) of a
PV bus after the solve: `:ok` up to the Newton tolerance, `:within_hysteresis`
up to the switching hysteresis (the switching logic tolerates that much on
purpose, so the final check has to as well), `:bounded` up to
`final_q_accept_pu` (accepted, with a warning), `:violation` beyond it (the
run is not accepted). With both thresholds at zero the check is as strict as
it was up to 0.17.2.
"""
function _final_q_limit_class(dev_pu::Float64; tol::Float64, q_hyst_pu::Float64, final_q_accept_pu::Float64)::Symbol
  dev_pu <= tol && return :ok
  dev_pu <= q_hyst_pu && return :within_hysteresis
  dev_pu <= final_q_accept_pu && return :bounded
  return :violation
end

"""
    classify_final_q_limits(net, Sbus_pu, bus_types, qmin_pu, qmax_pu; q_hyst_pu, final_q_accept_pu, tol) -> NamedTuple

The one final Q-limit check every enforcement mode ends with: which PV buses
lie over Qmax or under Qmin at the solved state, by how much, and what that
means for the run. Judged by the SIZE of the overshoot, not by the number of
buses: `power_flow.qlimits.hysteresis_pu` and
`power_flow.qlimits.final_q_accept_pu` are the two thresholds, see
[`_final_q_limit_class`](@ref).

Up to 0.17.2 the active set counted every PV overshoot beyond the Newton
tolerance as a remaining violation, while the classic modes stopped their
outer loop at the hysteresis and ran no final check at all: the same
operating point (a released machine 0.68 MVAr over Qmax inside a 1 MVAr
hysteresis) was rejected by one mode and accepted by the other.

`Sbus_pu` are the bus injections of the solved state (`V conj(Y V)`), the
machine output is the injection plus the bus load. Returns
`(rows, status, violations, bounded, within_hysteresis, max_dev_pu, buses)`:
one row `(bus, busI, side, q_pu, limit_pu, dev_pu, class)` per PV bus beyond
a limit by more than `tol`, sorted by the overshoot, `status` the worst class
as a run status (`:ok`, `:within_hysteresis`, `:bounded_q_limit_violation`,
`:remaining_pv_q_limit_violations`), `buses` the rows as `busI:class` joined
by `;` for the run metadata.
"""
function classify_final_q_limits(net::Net, Sbus_pu::AbstractVector{ComplexF64}, bus_types::AbstractVector{Symbol}, qmin_pu::AbstractVector, qmax_pu::AbstractVector; q_hyst_pu::Float64, final_q_accept_pu::Float64, tol::Float64)
  qload_pu = build_qload_pu(net)
  rows = NamedTuple{(:bus, :busI, :side, :q_pu, :limit_pu, :dev_pu, :class),Tuple{Int,Int,Symbol,Float64,Float64,Float64,Symbol}}[]
  for bus in eachindex(bus_types)
    bus_types[bus] == :PV || continue
    (bus <= length(qmin_pu) && bus <= length(qmax_pu) && bus <= length(Sbus_pu)) || continue
    q = imag(Sbus_pu[bus]) + (bus <= length(qload_pu) ? qload_pu[bus] : 0.0)
    side, limit, dev = :none, 0.0, 0.0
    if isfinite(qmax_pu[bus]) && q > qmax_pu[bus] + tol
      side, limit, dev = :high, Float64(qmax_pu[bus]), q - qmax_pu[bus]
    elseif isfinite(qmin_pu[bus]) && q < qmin_pu[bus] - tol
      side, limit, dev = :low, Float64(qmin_pu[bus]), qmin_pu[bus] - q
    end
    side === :none && continue
    push!(rows, (bus = bus, busI = _qlimit_original_bus_id(net, bus), side = side, q_pu = q, limit_pu = limit, dev_pu = dev,
      class = _final_q_limit_class(dev; tol = tol, q_hyst_pu = q_hyst_pu, final_q_accept_pu = final_q_accept_pu)))
  end
  sort!(rows; by = r -> -r.dev_pu)
  violations = count(r -> r.class === :violation, rows)
  bounded = count(r -> r.class === :bounded, rows)
  within = count(r -> r.class === :within_hysteresis, rows)
  status = violations > 0 ? :remaining_pv_q_limit_violations : bounded > 0 ? :bounded_q_limit_violation : within > 0 ? :within_hysteresis : :ok
  max_dev = isempty(rows) ? 0.0 : maximum(r.dev_pu for r in rows)
  buses = join((string(r.busI, ":", r.class) for r in rows), ";")
  return (rows = rows, status = status, violations = violations, bounded = bounded, within_hysteresis = within, max_dev_pu = max_dev, buses = buses)
end

"""
    final_q_check_line(status, rows, baseMVA) -> String

The final Q-limit check in one line, the same words in the text report,
`run.log`, the run metadata and on the Web UI result page: the run status and
one clause per row with bus, side, overshoot in pu and MVAr and class.
"""
function final_q_check_line(status::Symbol, rows, baseMVA::Real)::String
  status === :not_evaluated && return "Final Q-limit check: not evaluated (no converged solution)."
  isempty(rows) && return "Final Q-limit check: ok (no PV bus beyond its reactive limits)."
  parts = [@sprintf("bus %d %s by %.4f pu / %.2f MVAr (%s)", r.busI, r.side === :high ? "over Qmax" : "under Qmin", r.dev_pu, r.dev_pu * baseMVA, String(r.class)) for r in rows]
  return string("Final Q-limit check: ", String(status), ", ", join(parts, "; "), ".")
end

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
  final_q_accept_pu::Float64,
  tol::Float64,
  pv_table_rows::Int,
  qlimit_guard_accept_bounded_violations::Bool,
  qlimit_guard_max_remaining_violations::Int,
)
  final_Qload_pu = build_qload_pu(net)
  # Print summary only in verbose/trace mode; stay silent otherwise.
  qlimit_summary_io = (verbose > 0 || qlimit_trace_enabled) ? stdout : devnull

  qlimit_summary = _print_rectangular_qlimit_summary(qlimit_summary_io, net, V, Sbus_pu, bus_types, qmin_pu, qmax_pu, final_Qload_pu; q_hyst_pu = q_hyst_pu, tolerance_pu = tol, final_q_accept_pu = final_q_accept_pu, max_rows = pv_table_rows, max_console_rows = pv_table_rows)

  # the decision comes from the two-threshold classification (size of the
  # overshoot); the count-based guard below stays what it was, a policy to
  # accept a number of real violations
  final_q_check = classify_final_q_limits(net, Sbus_pu, bus_types, qmin_pu, qmax_pu; q_hyst_pu = q_hyst_pu, final_q_accept_pu = final_q_accept_pu, tol = tol)
  remaining_pv_violations = final_q_check.violations
  bounded_ok = qlimit_guard_accept_bounded_violations && remaining_pv_violations <= qlimit_guard_max_remaining_violations

  converged_ = converged
  rejection_reason_ = rejection_reason

  if remaining_pv_violations > 0 && !bounded_ok
    # Hard rejection: active PV-limit violations remain after the solve.
    verbose > 0 && @warn "Rectangular NR active-set failed because active PV Q-limit violations remain after the numerical solve." pv_violations = remaining_pv_violations ref_violations = qlimit_summary.ref_violations
    converged_ = false
    rejection_reason_ = :remaining_pv_q_limit_violations
  elseif remaining_pv_violations > 0 && bounded_ok
    # Soft acceptance: bounded residual violations accepted by policy.
    rejection_reason_ = :bounded_q_limit_violations_accepted
  elseif final_q_check.bounded > 0
    # accepted by size, said out loud: an overshoot beyond the hysteresis
    # but within final_q_accept_pu is a run the user should look at
    for r in final_q_check.rows
      r.class === :bounded || continue
      @warn "Final Q-limit check: bounded violation accepted" bus = r.busI side = r.side dev_pu = r.dev_pu dev_MVAr = r.dev_pu * net.baseMVA final_q_accept_pu = final_q_accept_pu
    end
    rejection_reason_ = :bounded_q_limit_violation
  end

  return (qlimit_summary = qlimit_summary, final_q_check = final_q_check, converged = converged_, rejection_reason = rejection_reason_)
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
  final_q_check,
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
  if !isempty(history) && !isfinite(history[end])
    final_reason = :nr_nonfinite
  end

  final_status = _rectangular_solver_status_symbol(numerical_converged, q_limit_active_set_ok, converged, final_reason)
  finite_history = filter(isfinite, history)

  status = (
    numerical_converged = numerical_converged,
    nr_converged = numerical_converged,
    active_set_converged = q_limit_active_set_ok,
    q_limit_active_set_ok = q_limit_active_set_ok,
    final_converged = converged,
    status = final_status,
    reason = final_reason,
    reason_text = _rectangular_rejection_reason_text(final_reason),
    pv_q_limit_violations = isnothing(final_q_check) ? (isnothing(qlimit_summary) ? 0 : qlimit_summary.pv_violations) : final_q_check.violations,
    ref_q_limit_violations = isnothing(qlimit_summary) ? 0 : qlimit_summary.ref_violations,
    # the final Q-limit check (classify_final_q_limits): the same fields for
    # every enforcement mode, :not_evaluated when no converged state exists
    final_q_check_status = isnothing(final_q_check) ? :not_evaluated : final_q_check.status,
    final_q_check_rows = isnothing(final_q_check) ? NamedTuple[] : final_q_check.rows,
    final_q_check_max_dev_pu = isnothing(final_q_check) ? 0.0 : final_q_check.max_dev_pu,
    final_q_check_buses = isnothing(final_q_check) ? "" : final_q_check.buses,
    final_pv_voltage_residual = final_pv_voltage_residual,
    final_mismatch = isempty(history) ? Inf : history[end],
    initial_mismatch = isempty(history) ? NaN : history[1],
    best_mismatch = isempty(finite_history) ? NaN : minimum(finite_history),
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

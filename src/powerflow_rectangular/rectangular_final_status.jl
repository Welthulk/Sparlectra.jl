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
  # a converged run with power_flow.qlimits.enabled = false has nothing to
  # check; saying "no converged solution" there was wrong (run.log)
  status === :qlimits_disabled && return "Final Q-limit check: not evaluated (Q-limit handling disabled)."
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
  wrong_branch_min_vn_kV::Float64 = 100.0,
  wrong_branch_low_vm_share::Float64 = 0.05,
  wrong_branch_max_bus_angle_deg::Float64 = 120.0,
  wrong_branch_max_plain_steps::Int = 20,
  wrong_branch_collapse_vm_pu::Float64 = 0.5,
  wrong_branch_max_reference_branch_angle_deg::Float64 = 90.0,
  plain_steps::Int = 0,
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
      min_vn_kV = wrong_branch_min_vn_kV,
      low_vm_share = wrong_branch_low_vm_share,
      max_bus_angle_deg = wrong_branch_max_bus_angle_deg,
      max_plain_steps = wrong_branch_max_plain_steps,
      collapse_vm_pu = wrong_branch_collapse_vm_pu,
      max_reference_branch_angle_deg = wrong_branch_max_reference_branch_angle_deg,
      plain_steps = plain_steps,
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

# The message of a run stopped by the switching cap: which setting stopped
# it, at which value, and the buses that reached it, so the reader can raise
# the cap or look at those machines (0.20.2: a stored max_switches of
# 2 was invisible in the Web UI, the old text said "max switching exceeded").
function _max_switching_reason_text(net, switch_counts::AbstractDict, max_switches::Int)::String
  cap = max(max_switches, 1)
  hit = sort!([bus for (bus, n) in switch_counts if n >= cap])
  # the bus names the user gave (busDict keys), not the component names
  bus_name = Dict{Int,String}(idx => String(name) for (name, idx) in net.busDict)
  names = String[string(get(bus_name, bus, "bus $(bus)"), " (", switch_counts[bus], "x)") for bus in hit]
  return string("Q-limit switching stopped: power_flow.qlimits.guard.max_switches = ", max_switches, " was reached by ", length(hit), " bus(es): ", join(names, ", "), ". Raise the value, or set power_flow.qlimits.guard.freeze_after_repeated_switching = false to let these buses keep switching.")
end

# The reasons of a run that used up power_flow.max_iter without converging,
# see _max_iter_stop_reason.
const QLIMIT_MAX_ITER_STOP_REASONS = (:iteration_limit, :type_identification_divergence, :numerical_divergence)

# How many of the last iterations _max_iter_stop_reason looks at.
const MAX_ITER_STOP_WINDOW = 3

"""
    _max_iter_stop_reason(pre_switch_mismatch, type_switched, tol; window = 3) -> Symbol

The reason of a Newton run that used up `power_flow.max_iter` without
converging, from the last `window` iterations. `pre_switch_mismatch[k]` is
the largest mismatch of iteration `k` measured before the Q-limit handling
of that iteration switched any bus, `type_switched[k]` whether it switched
one (PV to PQ or back). Checked in this order:

- `:type_identification_divergence`: a bus type changed in the window and
  the smallest pre-switch mismatch there is at most `10 tol`. The equations
  are almost solved, the PV/PQ assignment does not settle (the
  identification divergence of Zhao, Chiang, Li, Ju, PSCC 2008);
- `:iteration_limit`: no type change in the window, the mismatch fell
  from each iteration to the next over the full window, and the last
  mismatch is below the first mismatch of the run: the run was still
  converging when the limit came;
- `:numerical_divergence`: everything else: the mismatch did not fall,
  types changed at a large mismatch, the run fell over the window but ends
  above the mismatch it started with (it climbed first), or the run was
  shorter than the window.

The window and the factor 10 are fixed, there is no configuration key.
"""
function _max_iter_stop_reason(pre_switch_mismatch::AbstractVector{<:Real}, type_switched::AbstractVector{Bool}, tol::Real; window::Int = MAX_ITER_STOP_WINDOW)::Symbol
  n = min(length(pre_switch_mismatch), length(type_switched))
  n == 0 && return :numerical_divergence
  first_k = max(n - window + 1, 1)
  recent_switch = any(type_switched[k] for k in first_k:n)
  recent_switch && minimum(pre_switch_mismatch[k] for k in first_k:n) <= 10 * tol && return :type_identification_divergence
  falling = n >= window && all(pre_switch_mismatch[k] < pre_switch_mismatch[k-1] for k in (first_k+1):n)
  # a run that climbed away from its start and falls again at the end is
  # no run close to convergence (sp_case60 with a narrow band climbs to
  # 1e+3 and more and then falls for a while)
  below_start = pre_switch_mismatch[n] < pre_switch_mismatch[1]
  !recent_switch && falling && below_start && return :iteration_limit
  return :numerical_divergence
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
  final_q_check_skip::Symbol = :no_converged_solution;
  switch_counts::AbstractDict = Dict{Int,Int}(),
  max_switches::Int = 0,
  pre_switch_history::AbstractVector{<:Real} = Float64[],
  type_switched::AbstractVector{Bool} = Bool[],
)
  final_reason = converged ? :none : rejection_reason

  # a run at max_iter carries one of QLIMIT_MAX_ITER_STOP_REASONS, set by
  # the solver from the last iterations (_max_iter_stop_reason); the former
  # reclassification as _active_set_unstable counted switches over the
  # whole run, without regard to the mismatch
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
    reason_text = final_reason === :max_switching_exceeded ? _max_switching_reason_text(net, switch_counts, max_switches) : _rectangular_rejection_reason_text(final_reason),
    pv_q_limit_violations = isnothing(final_q_check) ? (isnothing(qlimit_summary) ? 0 : qlimit_summary.pv_violations) : final_q_check.violations,
    ref_q_limit_violations = isnothing(qlimit_summary) ? 0 : qlimit_summary.ref_violations,
    # the final Q-limit check (classify_final_q_limits): the same fields for
    # every enforcement mode, :not_evaluated when no converged state exists,
    # :qlimits_disabled when the run converged with Q-limit handling off
    final_q_check_status = isnothing(final_q_check) ? (final_q_check_skip === :qlimits_disabled ? :qlimits_disabled : :not_evaluated) : final_q_check.status,
    final_q_check_rows = isnothing(final_q_check) ? NamedTuple[] : final_q_check.rows,
    final_q_check_max_dev_pu = isnothing(final_q_check) ? 0.0 : final_q_check.max_dev_pu,
    final_q_check_buses = isnothing(final_q_check) ? "" : final_q_check.buses,
    final_pv_voltage_residual = final_pv_voltage_residual,
    final_mismatch = isempty(history) ? Inf : history[end],
    # the largest mismatch of each iteration before its Q-limit handling
    # switched a bus (mismatch_history holds the value after the switch,
    # the mismatch of the equations the next Newton step solves), and the
    # iterations in which a bus type changed
    mismatch_history_pre_switch = collect(Float64, pre_switch_history),
    type_switch_iterations = [k for k in eachindex(type_switched) if type_switched[k]],
    initial_mismatch = isempty(history) ? NaN : history[1],
    best_mismatch = isempty(finite_history) ? NaN : minimum(finite_history),
    nr_initial_mismatch = isempty(history) ? NaN : history[1],
    nr_final_mismatch = isempty(history) ? Inf : history[end],
    mismatch_diagnostics...,
    pv_pq_switching_events = qlimit_pv_to_pq_event_count(net),
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
    wrong_branch_level_kV = branch_quality.level_kV,
    wrong_branch_level_low_vm_count = branch_quality.level_low_vm_count,
    wrong_branch_level_bus_count = branch_quality.level_bus_count,
    wrong_branch_max_bus_angle_deg = branch_quality.max_bus_angle_deg,
    wrong_branch_plain_steps = branch_quality.plain_steps,
    wrong_branch_plain_steps_exceeded = branch_quality.plain_steps_exceeded,
    wrong_branch_angle_spread_deg = branch_quality.angle_spread_deg,
    wrong_branch_max_branch_angle_deg = branch_quality.max_branch_angle_deg,
    wrong_branch_branch_angle_violation_count = branch_quality.branch_angle_violation_count,
    wrong_branch_worst_branch_angle_deg = branch_quality.worst_branch_angle_deg,
    # issue #462: the largest angle across a branch at the reference bus
    # (NaN when not checked) and that branch in case bus numbers
    wrong_branch_reference_branch_angle_deg = branch_quality.reference_branch_angle_deg,
    wrong_branch_reference_branch = branch_quality.reference_branch,
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

A wrong-branch voltage finding after a start without angles behind a phase
shifter (`_wrong_branch_pst_start_hint`) is logged as a warning at every
verbosity level, since the hint names the configuration key that avoids it.

# Why
Status persistence and user-visible summary should be emitted from one place to
avoid divergence between stored diagnostics and console output.
"""
function _store_and_print_rectangular_final_status!(net, status, verbose::Int)
  stored_status = _set_rectangular_pf_status!(net, status)
  if verbose > 0
    _print_rectangular_convergence_summary(stdout, stored_status)
  end
  hint = _wrong_branch_pst_start_hint(stored_status)
  isempty(hint) || @warn string("Wrong-branch check: ", uppercase(String(stored_status.wrong_branch_status)), " (", stored_status.wrong_branch_reason, "). ", hint)
  return stored_status
end

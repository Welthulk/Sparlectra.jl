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

# Diagnostic and Q-limit artifacts are API/Web UI contracts. Keep filenames,
# log markers, and CSV schemas stable unless the artifact views and tests are
# updated together.

function _write_api_timing_summary(io::IO, result::SparlectraRunResult, config::SparlectraConfig, phases::AbstractDict = Dict{Symbol,Float64}())
  benchmark_median = result.performance_profile isa AbstractDict ? get(result.performance_profile, :benchmark_median_s, nothing) : nothing
  println(io, "Timing")
  println(io, "------")
  result.solver_elapsed_s === nothing || println(io, "Solver time : ", round(result.solver_elapsed_s; digits = 6), " s")
  haskey(phases, :artifact_writing) && println(io, "Output time : ", round(Float64(phases[:artifact_writing]); digits = 6), " s")
  haskey(phases, :total) && println(io, "Total time  : ", round(Float64(phases[:total]); digits = 6), " s")
  if config.benchmark.enabled
    println(io, "benchmark_median:     ", benchmark_median === nothing ? "n/a" : "$(round(Float64(benchmark_median); digits = 6)) s")
    println(io, "benchmark_samples:    ", config.benchmark.samples)
  end
  println(io, "iterations:           ", result.iterations)
  println(io, "final_mismatch:       ", isfinite(result.final_mismatch) ? result.final_mismatch : (isnan(result.final_mismatch) ? "NaN" : "unavailable"))
  println(io, "final_mismatch_status:", isfinite(result.final_mismatch) ? "finite" : (isnan(result.final_mismatch) ? "nonfinite" : "not_reported_by_solver"))
  println(io, "final_status:         ", result.outcome)
  println(io, "final_outcome:        ", result.final_converged ? "converged" : result.reason)
  if !result.numerical_converged
    println(io, "q_limit_validity:     last-iteration diagnostic only; NR did not converge")
  end
  println(io)
  return nothing
end

_final_mismatch_status(result::SparlectraRunResult)::String = isfinite(result.final_mismatch) ? "finite" : (isnan(result.final_mismatch) ? "nonfinite" : "not_reported_by_solver")
_final_mismatch_available(result::SparlectraRunResult)::Bool = isfinite(result.final_mismatch)
_final_mismatch_reason(result::SparlectraRunResult)::String = _final_mismatch_available(result) ? "available" : (isnan(result.final_mismatch) ? "solver returned nonfinite final mismatch" : "solver returned unavailable final mismatch")
_diag_value(result::SparlectraRunResult, name::Symbol, default = "unavailable") = hasproperty(result.diagnostics, name) ? getproperty(result.diagnostics, name) : default

function _write_powerflow_mismatch_diagnostics(io::IO, result::SparlectraRunResult; mode::Symbol)
  println(io, "Mismatch summary")
  println(io, "----------------")
  println(io, "initial_mismatch: ", _diag_value(result, :initial_mismatch, _diag_value(result, :nr_initial_mismatch, "unavailable")))
  println(io, "final_mismatch: ", isfinite(result.final_mismatch) ? result.final_mismatch : (isnan(result.final_mismatch) ? "NaN" : "unavailable"))
  println(io, "final_mismatch_status: ", _final_mismatch_status(result))
  println(io, "final_mismatch_available: ", _final_mismatch_available(result))
  _final_mismatch_available(result) || println(io, "final_mismatch_reason_if_unavailable: ", _final_mismatch_reason(result))
  println(io, "iterations: ", result.iterations)
  println(io, "converged: ", result.final_converged)
  println(io, "failure_reason: ", result.final_converged ? "none" : result.reason)
  for key in (:selected_start_candidate, :raw_mismatch, :dc_mismatch, :projected_mismatch,
              :current_iteration_initial_mismatch, :current_iteration_final_mismatch,
              :current_iteration_accepted, :current_iteration_reason)
    println(io, key, ": ", _diag_value(result, key))
  end
  if mode === :full
    for key in (:best_blend_mismatch, :start_projection_mismatch_before,
                :start_projection_mismatch_after, :current_iteration_candidate_mismatch,
                :nr_initial_mismatch, :nr_final_mismatch, :max_active_power_mismatch,
                :max_reactive_power_mismatch,
                :max_voltage_residual_or_setpoint_residual_where_available,
                :requested_angle_mode, :requested_voltage_mode, :selection_reason,
                :fallback_to_raw, :raw_fallback_reason, :dc_angle_start_built,
                :dc_angle_start_valid, :dc_angle_start_applied, :dc_angle_min_deg,
                :dc_angle_max_deg, :dc_angle_spread_deg, :dc_angle_mean_deg,
                :dc_angle_std_deg, :dc_angle_clipped_count, :dc_angle_clip_limit_deg,
                :dc_max_branch_angle_deg, :dc_branch_angle_violation_count,
                :worst_dc_branch_from_bus, :worst_dc_branch_to_bus,
                :worst_dc_branch_angle_deg, :dc_voltage_magnitude_min,
                :dc_voltage_magnitude_max, :dc_mismatch_ratio_vs_raw,
                :requested_dc_worse_than_raw, :dc_mismatch_growth_factor)
      println(io, key, ": ", _diag_value(result, key))
    end
  end
  println(io)
  return nothing
end

_diag_symbol(result::SparlectraRunResult, name::Symbol, default::Symbol = :unavailable)::Symbol = (v = _diag_value(result, name, default); v isa Symbol ? v : default)

# Human-readable classification for a worst-mismatch equation symbol (:P, :Q,
# or :voltage_setpoint), used only in narrative prose, not in the key:value
# lines consumed by tests/tooling.
function _mismatch_equation_label(equation::Symbol)::String
  equation === :P && return "active-power balance"
  equation === :Q && return "reactive-power balance"
  equation === :voltage_setpoint && return "PV/REF voltage setpoint"
  return string(equation)
end

"""
    _write_diagnosis_narrative(io, result)

Write the "Diagnosis" and "Recommendations" sections of `diagnose.log`: a short
prose summary of what the worst-mismatch bus, mismatch trend, and autodamp
behavior indicate, followed by a "Worst-mismatch bus" data section (reusing
the same `top_mismatch_rows` already computed once per solve by
[`_rectangular_mismatch_diagnostics`](@ref)) and a "Branch anomalies" section
from [`_branch_anomaly_diagnostics`](@ref) for the branches incident to that
bus. This turns a flat key/value dump into a report that names a likely cause
and a next step, instead of only reporting numbers.
"""
function _write_diagnosis_narrative(io::IO, result::SparlectraRunResult)
  println(io, "Diagnosis")
  println(io, "---------")
  if result.final_converged
    println(io, "The run converged in ", result.iterations, " iteration(s); final mismatch ", result.final_mismatch, " is within tolerance. No root-cause analysis is needed.")
    println(io)
    return nothing
  end

  bus_id = _diag_value(result, :worst_mismatch_bus_id, "unavailable")
  equation = _diag_symbol(result, :worst_mismatch_equation)
  mismatch_value = _diag_value(result, :worst_mismatch_value, NaN)
  trend = _diag_symbol(result, :mismatch_history_trend)
  floor_hits = _diag_value(result, :autodamp_floor_hits, 0)
  autodamp_failure = _diag_value(result, :autodamp_failure, false)

  if bus_id == "unavailable"
    println(io, "The run did not converge (reason: ", result.reason, "), and no per-bus mismatch breakdown is available for this failure mode.")
  else
    println(io, "The run did not converge (reason: ", result.reason, "). The largest remaining residual is at bus ", bus_id, " (", _mismatch_equation_label(equation), "), with a mismatch of ", mismatch_value, ".")
  end

  trend_text = if trend === :diverging_to_nonfinite
    "The mismatch trend is diverging to a non-finite value — later iterations moved away from a solution rather than toward one."
  elseif trend === :oscillatory
    "The mismatch trend is oscillatory — successive iterations alternate between improving and worsening, which often indicates step control is fighting an unstable region rather than making steady progress."
  elseif trend === :stagnant
    "The mismatch trend is stagnant — recent iterations made little to no progress, which often indicates the step control is stuck rather than diverging."
  elseif trend === :monotonic
    "The mismatch trend is monotonically improving but did not reach tolerance within the configured iteration budget."
  else
    ""
  end
  isempty(trend_text) || println(io, trend_text)

  if autodamp_failure === true
    println(io, "Autodamp diagnostics show repeated floor hits (", floor_hits, ") together with non-improving trial steps, consistent with the step-control damping being unable to make progress from the current start point.")
  end
  println(io)

  bus_index = _diag_value(result, :worst_mismatch_bus_index, "unavailable")
  branch_rows = bus_index isa Integer ? _branch_anomaly_diagnostics(result.net, bus_index) : NamedTuple[]
  flagged_rows = filter(row -> !isempty(row.flags), branch_rows)

  println(io, "Recommendations")
  println(io, "---------------")
  recommendations = String[]
  if !isempty(flagged_rows)
    push!(recommendations, "Branch anomalies were found at the worst-mismatch bus (see \"Branch anomalies at worst-mismatch bus\" below) — review the flagged branch parameters and the MATPOWER import conventions (matpower_import.auto_profile/shift_sign/shift_unit/ratio) before adjusting solver settings.")
  end
  if trend === :diverging_to_nonfinite || autodamp_failure === true
    push!(recommendations, "Try a different start profile (power_flow.start_mode.angle_mode/voltage_mode) before tuning step control further — a poor start is a common cause of floor-hitting autodamp behavior.")
  end
  if trend === :oscillatory
    push!(recommendations, "Consider enabling power_flow.merit (line search) or power_flow.trust_region (step_mode=:dogleg) to stabilize oscillatory steps.")
  end
  if trend === :stagnant
    push!(recommendations, "Check power_flow.qlimits and wrong-branch detection settings — a stagnant trend can indicate the active-set/outer loop, not the Newton step itself, is stuck.")
  end
  push!(recommendations, "Run run_fixed_reference_self_check on this case to check whether the residual is already present at the case's own stored VM/VA (a network-model issue) rather than introduced by the solver's start guess.")
  for (i, rec) in enumerate(recommendations)
    println(io, i, ". ", rec)
  end
  println(io)

  println(io, "Worst-mismatch bus")
  println(io, "-------------------")
  println(io, "worst_mismatch_bus_id: ", bus_id)
  println(io, "worst_mismatch_bus_index: ", bus_index)
  println(io, "worst_mismatch_equation: ", equation)
  println(io, "worst_mismatch_value: ", mismatch_value)
  println(io, "mismatch_history_trend: ", trend)
  println(io, "mismatch_history_first: ", _diag_value(result, :mismatch_history_first, NaN))
  println(io, "mismatch_history_last: ", _diag_value(result, :mismatch_history_last, NaN))
  println(io, "mismatch_history_best: ", _diag_value(result, :mismatch_history_best, NaN))
  println(io, "autodamp_floor_hits: ", floor_hits)
  println(io, "autodamp_failure: ", autodamp_failure)
  println(io, "top_mismatch_rows: ", _format_top_mismatch_rows(_diag_value(result, :top_mismatch_rows, NamedTuple[])))
  println(io)

  println(io, "Branch anomalies at worst-mismatch bus")
  println(io, "---------------------------------------")
  println(io, isempty(branch_rows) ? "unavailable" : _format_branch_anomaly_rows(branch_rows))
  println(io)
  return nothing
end

"""
    _write_execution_failure_diagnostics(path, reason, message)

Write a best-effort `diagnose.log` when a run fails before a
`SparlectraRunResult` is available — e.g. `run_sparlectra` throwing for an
AC-island failure with a non-finite mismatch before the first Newton
iteration (`stage=before_nr`), or an unsupported MATPOWER DC-line
convention. Without this, `run_diagnostics = true` silently produced no
`diagnose.log` at all for these early/exception failures, since the normal
`_write_powerflow_diagnostics` call site is only reached after
`run_sparlectra` returns successfully. Rich per-bus/branch diagnostics need a
solved network, which this failure mode never reaches; per-island detail,
when available, lives in `ac_island_<id>_solver.log`.
"""
function _write_execution_failure_diagnostics(path::AbstractString, reason::AbstractString, message::AbstractString)
  open(path, "w") do io
    println(io, "Sparlectra PowerFlow diagnostics")
    println(io, "================================")
    println(io, "outcome: execution_failure")
    println(io, "reason: ", reason)
    println(io)
    println(io, "Diagnosis")
    println(io, "---------")
    println(io, "The run failed before a solved network/result was available, so no per-bus/branch mismatch breakdown can be reported here.")
    println(io)
    println(io, message)
    println(io)
    println(io, "Recommendations")
    println(io, "---------------")
    println(io, "1. If this names an AC island, check that island's ac_island_<id>_solver.log and ac_island_solver_summary.csv for its start-of-solve mismatch and settings.")
    println(io, "2. Run run_fixed_reference_self_check on this case to check whether a large residual is already present at the case's own stored VM/VA (a network-model issue) rather than introduced by the solver's start guess.")
  end
  return path
end

function _write_powerflow_diagnostics(path::AbstractString, result::SparlectraRunResult; diagnostic_fn = nothing, mode::Symbol = :compact)
  open(path, "w") do io
    println(io, "Sparlectra PowerFlow diagnostics")
    println(io, "================================")
    try
      if diagnostic_fn === nothing
        println(io, "outcome: ", result.outcome)
        println(io, "reason: ", result.reason_text)
        println(io)
        _write_diagnosis_narrative(io, result)
        _write_powerflow_mismatch_diagnostics(io, result; mode = mode)
        if mode === :full
          printQLimitLog(result.net; io)
          println(io)
          printPVQLimitsTable(result.net; io)
          println(io)
          _write_final_q_limit_validation(io, result)
        else
          println(io, "Q-limit validity: ", result.numerical_converged ? "valid final validation available in q_limit.log" : "last-iteration diagnostic only; NR did not converge")
        end
      else
        diagnostic_fn(io, result)
      end
    catch err
      println(io)
      println(io, "Diagnostic generation failed; the PowerFlow result remains valid.")
      println(io, sprint(showerror, err))
    end
  end
  return path
end

function _write_final_q_limit_validation(io::IO, result::SparlectraRunResult)
  pv_violations = hasproperty(result.diagnostics, :pv_q_limit_violations) ? result.diagnostics.pv_q_limit_violations : 0
  ref_violations = hasproperty(result.diagnostics, :ref_q_limit_violations) ? result.diagnostics.ref_q_limit_violations : 0
  if result.numerical_converged && pv_violations == 0 && ref_violations == 0
    println(io, "Final Q-limit validation")
    # Do not report outer-loop switching candidates as final violations.
    # The final validation section must reflect the final bus types and clamped Q values.
    println(io, "  Final PV/REF Q-limit validation: OK")
    println(io, "  PV/REF Q-limit violations: 0")
    println(io, "  Outer-loop selected violations, if any, are reported in the event sections above.")
    return (q_violations = 0, v_violations = 0)
  end
  return printFinalLimitValidation(result.net; io, converged = result.numerical_converged)
end

function _write_q_limit_detail_artifacts(output_path::AbstractString, net::Net; format = "technical")::Vector{String}
  artifacts = String[]
  events = [(iteration = ev.iter, bus = ev.bus, side = String(ev.side)) for ev in net.qLimitLog]
  if !isempty(events)
    _write_namedtuple_csv(joinpath(output_path, "q_limit_events.csv"), events, (:iteration, :bus, :side); format = format)
    push!(artifacts, "q_limit_events.csv")
  end
  rows = NamedTuple[]
  snapshot_rows = isempty(net.qLimitInitialPVRows) ? snapshotPVQLimits!(net) : net.qLimitInitialPVRows
  for row in snapshot_rows
    bus = row.bus
    push!(rows, (
      bus = bus,
      qmin_pu = row.qmin_MVAr / net.baseMVA,
      qmax_pu = row.qmax_MVAr / net.baseMVA,
      qmin_MVAr = row.qmin_MVAr,
      qmax_MVAr = row.qmax_MVAr,
    ))
  end
  if !isempty(rows)
    _write_namedtuple_csv(joinpath(output_path, "q_limit_initial_limits.csv"), rows, (:bus, :qmin_pu, :qmax_pu, :qmin_MVAr, :qmax_MVAr); format = format)
    push!(artifacts, "q_limit_initial_limits.csv")
  end
  rect_status = rectangular_pf_status(net)
  if rect_status !== nothing && hasproperty(rect_status, :matpower_outer_loop)
    outer_rows = collect(rect_status.matpower_outer_loop)
    if !isempty(outer_rows)
      rows = [
        (
          outer_iteration = row.outer_iter,
          mode = String(row.mode),
          selected_count = count(other -> other.outer_iter == row.outer_iter, outer_rows),
          max_violation_MVAr = maximum(other.violation_mvar for other in outer_rows if other.outer_iter == row.outer_iter),
          ref_changed = row.ref_changed,
          bus = row.bus_i,
          generator_index = row.gen_index,
          violation_side = String(row.violation_side),
          action = String(row.action),
        )
        for row in outer_rows
      ]
      columns = (:outer_iteration, :mode, :selected_count, :max_violation_MVAr, :ref_changed, :bus, :generator_index, :violation_side, :action)
      _write_namedtuple_csv(joinpath(output_path, "q_limit_classic_outer_loop.csv"), rows, columns; format = format)
      push!(artifacts, "q_limit_classic_outer_loop.csv")
    end
  end
  return artifacts
end

function _write_q_limit_log_artifact(output_path::AbstractString, result::SparlectraRunResult, metadata::AbstractDict)::String
  path = joinpath(output_path, Q_LIMIT_LOG_ARTIFACT)
  open(path, "w") do io
    println(io, "Sparlectra Q-limit diagnostics")
    println(io, "==============================")
    println(io)
    _write_resolved_q_limit_options(io, metadata)
    println(io, "Q-limit detail artifact  : ", Q_LIMIT_LOG_ARTIFACT)
    println(io)
    if get(metadata, "qlimits_enabled", false) == false
      println(io, "Q-limit guard effective  : false")
      println(io, "Q-limit enforcement      : skipped")
      println(io, "Q-limit active set       : not run")
      println(io, "PV->PQ Q-limit events    : 0")
      println(io, "Guarded PV->PQ locks     : 0")
      return Q_LIMIT_LOG_ARTIFACT
    end

    println(io, "Initial PV Q-limit table")
    println(io, "------------------------")
    printPVQLimitsTable(result.net; io, max_rows = typemax(Int))
    println(io)
    println(io, "Q-limit guard actions")
    println(io, "---------------------")
    guarded = hasproperty(result.diagnostics, :guarded_narrow_q_pv_buses) ? result.diagnostics.guarded_narrow_q_pv_buses : 0
    println(io, "Guarded PV->PQ locks     : ", guarded)
    println(io)
    println(io, "PV->PQ and PQ->PV event details")
    println(io, "-------------------------------")
    printQLimitLog(result.net; io, max_rows = typemax(Int))
    println(io)
    rect_status = rectangular_pf_status(result.net)
    if rect_status !== nothing && hasproperty(rect_status, :qlimit_enforcement_mode) && rect_status.qlimit_enforcement_mode in (:classic_simultaneous, :classic_one_at_a_time)
      println(io, "Classical Q-limit outer loop")
      println(io, "--------------------------------------")
      println(io, "Q-limit handling enabled    : true")
      println(io, "active-set switching used   : false")
      println(io, "classical outer-loop used   : true")
      println(io, "mode                       : ", rect_status.qlimit_enforcement_mode)
      println(io, "base_pf_converged          : ", getproperty(rect_status, :base_pf_converged))
      println(io, "qlimit_enforcement_started : ", getproperty(rect_status, :qlimit_enforcement_started))
      println(io, "final_outcome              : ", getproperty(rect_status, :final_outcome))
      println(io, "outer_iterations           : ", getproperty(rect_status, :matpower_outer_iterations))
      println(io, "PV->PQ Q-limit events      : ", length(result.net.qLimitLog))
      println(io, "PQ->PV re-enable events    : 0")
      println(io)
    else
      println(io, "Q-limit handling mode")
      println(io, "---------------------")
      println(io, "Q-limit handling enabled    : true")
      println(io, "active-set switching used   : true")
      println(io, "classical outer-loop used   : false")
      println(io, "PV->PQ Q-limit events      : ", length(result.net.qLimitLog))
      println(io, "PQ->PV re-enable events    : ", hasproperty(result.diagnostics, :qlimit_reenable_events) ? result.diagnostics.qlimit_reenable_events : 0)
      println(io)
    end
    println(io, "Final PV Q-limit active-set check")
    println(io, "---------------------------------")
    pv_violations = hasproperty(result.diagnostics, :pv_q_limit_violations) ? result.diagnostics.pv_q_limit_violations : 0
    ref_violations = hasproperty(result.diagnostics, :ref_q_limit_violations) ? result.diagnostics.ref_q_limit_violations : 0
    _write_final_q_limit_validation(io, result)
    final_ref_status = result.numerical_converged ? (ref_violations > 0 ? "WARN" : "OK") : "SKIPPED (NR did not converge)"
    final_pvref_status = result.numerical_converged ? (pv_violations > 0 ? "FAIL" : (ref_violations > 0 ? "WARN" : "OK")) : "SKIPPED (NR did not converge)"
    println(io, "Final REF Q-limit diagnostic: ", final_ref_status)
    println(io, "  REF violations: ", ref_violations)
    println(io, "Final PV/REF Q-limit check: ", final_pvref_status)
    println(io, "  remaining violations: ", pv_violations + ref_violations, " (PV ", pv_violations, ", REF ", ref_violations, ")")
    println(io)
    println(io, "Q-Limit Active-Set Summary")
    println(io, "--------------------------")
    summary_status = (
      numerical_converged = result.numerical_converged,
      final_mismatch = result.final_mismatch,
      q_limit_active_set_ok = hasproperty(result.diagnostics, :q_limit_active_set_ok) ? result.diagnostics.q_limit_active_set_ok : false,
      pv_pq_switching_events = hasproperty(result.diagnostics, :pv_pq_switching_events) ? result.diagnostics.pv_pq_switching_events : length(result.net.qLimitLog),
      qlimit_active_set_changes = hasproperty(result.diagnostics, :qlimit_active_set_changes) ? result.diagnostics.qlimit_active_set_changes : 0,
      qlimit_reenable_events = hasproperty(result.diagnostics, :qlimit_reenable_events) ? result.diagnostics.qlimit_reenable_events : 0,
      oscillating_buses = hasproperty(result.diagnostics, :oscillating_buses) ? result.diagnostics.oscillating_buses : 0,
      guarded_narrow_q_pv_buses = guarded,
      status = result.outcome,
    )
    _print_qlimit_active_set_summary(io, summary_status)
  end
  return Q_LIMIT_LOG_ARTIFACT
end

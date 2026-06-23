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

using UUIDs: uuid4

const _SPARLECTRA_API_SCHEMA_VERSION = "1.0"
const WEBUI_PERFORMANCE_TIMING_VALUES = (:off, :compact, :full)
const Q_LIMIT_LOG_ARTIFACT = "q_limit.log"

mutable struct PowerFlowPhaseTimingRecorder
  timings::Vector{Dict{String,Any}}
  active_index::Union{Nothing,Int}
end

PowerFlowPhaseTimingRecorder() = PowerFlowPhaseTimingRecorder(Dict{String,Any}[], nothing)

function _api_elapsed_seconds(start_ns::UInt64)::Float64
  return (time_ns() - start_ns) / 1.0e9
end

function _api_datetime_string(value::Dates.DateTime)::String
  return Dates.format(value, dateformat"yyyy-mm-ddTHH:MM:SS.sss") * "Z"
end

function _phase_elapsed_seconds(started_at::AbstractString, ended_at::AbstractString)::Float64
  start = Dates.DateTime(replace(String(started_at), "Z" => ""), dateformat"yyyy-mm-ddTHH:MM:SS.sss")
  finish = Dates.DateTime(replace(String(ended_at), "Z" => ""), dateformat"yyyy-mm-ddTHH:MM:SS.sss")
  return max(0.0, Dates.value(finish - start) / 1000)
end

function _complete_active_phase!(recorder::PowerFlowPhaseTimingRecorder, status::AbstractString = "completed")
  index = recorder.active_index
  index === nothing && return nothing
  timing = recorder.timings[index]
  haskey(timing, "ended_at") && timing["ended_at"] !== nothing && return nothing
  ended_at = _api_datetime_string(Dates.now(Dates.UTC))
  timing["ended_at"] = ended_at
  timing["elapsed_seconds"] = _phase_elapsed_seconds(timing["started_at"], ended_at)
  timing["status"] = String(status)
  recorder.active_index = nothing
  return timing
end

function _start_service_phase!(recorder::PowerFlowPhaseTimingRecorder, phase::AbstractString; status::AbstractString = "running")
  phase_name = String(phase)
  if recorder.active_index !== nothing && recorder.timings[recorder.active_index]["phase"] == phase_name
    return recorder.timings[recorder.active_index]
  end
  _complete_active_phase!(recorder)
  timing = Dict{String,Any}("phase" => phase_name, "started_at" => _api_datetime_string(Dates.now(Dates.UTC)), "ended_at" => nothing, "elapsed_seconds" => nothing, "status" => String(status))
  push!(recorder.timings, timing)
  recorder.active_index = length(recorder.timings)
  return timing
end

function _phase_elapsed_lookup(phase_timings::AbstractVector, phase::AbstractString)
  for timing in phase_timings
    get(timing, "phase", "") == phase && return get(timing, "elapsed_seconds", nothing)
  end
  return nothing
end

function _api_timing_mode(value)::Symbol
  mode = value isa Symbol ? value : Symbol(lowercase(strip(String(value))))
  mode in WEBUI_PERFORMANCE_TIMING_VALUES || throw(ArgumentError("Performance timing mode must be one of $(join(WEBUI_PERFORMANCE_TIMING_VALUES, ", ")); got $(repr(value))."))
  return mode
end

function _write_api_timing_summary(io::IO, result::SparlectraRunResult, config::SparlectraConfig, phases::AbstractDict = Dict{Symbol,Float64}())
  benchmark_median = result.performance_profile isa AbstractDict ? get(result.performance_profile, :benchmark_median_s, nothing) : nothing
  representative_time = result.performance_profile isa AbstractDict ? get(result.performance_profile, :representative_elapsed_s, result.elapsed_s) : result.elapsed_s
  println(io, "Timing")
  println(io, "------")
  println(io, "Solver time : ", result.solver_elapsed_s === nothing ? "n/a" : "$(round(result.solver_elapsed_s; digits = 6)) s")
  haskey(phases, :artifact_writing) && println(io, "Output time : ", round(Float64(phases[:artifact_writing]); digits = 6), " s")
  println(io, "Wall time   : ", round(Float64(representative_time); digits = 6), " s")
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

function _write_powerflow_diagnostics(path::AbstractString, result::SparlectraRunResult; diagnostic_fn = nothing, mode::Symbol = :compact)
  open(path, "w") do io
    println(io, "Sparlectra PowerFlow diagnostics")
    println(io, "================================")
    try
      if diagnostic_fn === nothing
        println(io, "outcome: ", result.outcome)
        println(io, "reason: ", result.reason_text)
        println(io)
        _write_powerflow_mismatch_diagnostics(io, result; mode = mode)
        if mode === :full
          printQLimitLog(result.net; io)
          println(io)
          printPVQLimitsTable(result.net; io)
          println(io)
          printFinalLimitValidation(result.net; io, converged = result.numerical_converged)
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
    if result.numerical_converged
      printFinalLimitValidation(result.net; io, converged = true)
    else
      printFinalLimitValidation(result.net; io, converged = false)
    end
    pv_violations = hasproperty(result.diagnostics, :pv_q_limit_violations) ? result.diagnostics.pv_q_limit_violations : 0
    ref_violations = hasproperty(result.diagnostics, :ref_q_limit_violations) ? result.diagnostics.ref_q_limit_violations : 0
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

function _write_service_phase_summary(io::IO, phase_timings::AbstractVector)
  println(io, "Phase timings")
  println(io, "-------------")
  for timing in phase_timings
    elapsed = get(timing, "elapsed_seconds", nothing)
    elapsed_text = elapsed === nothing ? "running" : string(round(Float64(elapsed); digits = 6), " s")
    println(io, "  ", rpad(String(get(timing, "phase", "unknown")) * ":", 31), elapsed_text, " (", get(timing, "status", "unknown"), ")")
  end
  println(io)
  return nothing
end

function _write_large_case_timing_summary(io::IO, case_path::AbstractString, phase_timings::AbstractVector, result::Union{Nothing,SparlectraRunResult})
  println(io, "Large case timing summary")
  println(io, "-------------------------")
  println(io, "case_file: ", basename(case_path))
  println(io, "case_extension: ", lowercase(splitext(case_path)[2]))
  println(io, "case_size_bytes: ", isfile(case_path) ? filesize(case_path) : "n/a")
  if result !== nothing && result.net !== nothing
    println(io, "n_bus: ", length(result.net.nodeVec))
    println(io, "n_branch: ", length(result.net.branchVec))
  end
  for (label, phase) in (
    ("reading_matpower_case_seconds", "reading_matpower_case"),
    ("building_sparlectra_net_seconds", "building_sparlectra_net"),
    ("building_ybus_seconds", "building_ybus"),
    ("preparing_start_values_seconds", "preparing_start_values"),
    ("start_projection_seconds", "start_projection"),
    ("solving_powerflow_seconds", "solving_powerflow"),
    ("artifact_seconds", "writing_artifacts"),
    ("total_service_seconds", "total_service"),
  )
    elapsed = _phase_elapsed_lookup(phase_timings, phase)
    elapsed === nothing || println(io, label, ": ", round(Float64(elapsed); digits = 6))
  end
  println(io)
  return nothing
end

function _resolved_matpower_import_runtime_options(config::SparlectraConfig)::Dict{String,Any}
  return Dict{String,Any}(
    "matpower_auto_profile" => String(config.matpower.auto_profile),
    "matpower_ratio" => String(config.matpower.ratio),
    "matpower_shift_sign" => config.matpower.shift_sign,
    "matpower_shift_unit" => String(config.matpower.shift_unit),
    "matpower_bus_shunt_model" => String(config.matpower.bus_shunt_model),
    "matpower_pv_voltage_source" => String(config.matpower.pv_voltage_source),
    "matpower_compare_reference" => String(config.matpower.compare_voltage_reference),
  )
end

function _resolved_q_limit_runtime_options(config::SparlectraConfig)::Dict{String,Any}
  qlimits_enabled = !config.powerflow.qlimits.ignore_q_limits
  preview_mode = config.output.console_q_limit_events === :full ? "full" : "summary"
  return Dict{String,Any}(
    "qlimits_enabled" => qlimits_enabled,
    "qlimit_enforcement_mode" => String(config.powerflow.qlimits.enforcement_mode),
    "qlimit_guard_enabled" => qlimits_enabled && config.powerflow.qlimits.guard,
    "logfile_diagnostics" => String(config.output.logfile_diagnostics),
    "console_q_limit_events" => String(config.output.console_q_limit_events),
    "q_limit_preview_mode" => preview_mode,
    "q_limit_runlog_max_rows" => config.output.console_q_limit_events === :summary || config.output.console_q_limit_events === :off ? 0 : config.output.console_max_rows,
    "q_limit_detail_artifacts" => Q_LIMIT_LOG_ARTIFACT,
  )
end

_metadata_kwargs(metadata::AbstractDict) = (; (Symbol(key) => value for (key, value) in metadata)...)

function _write_resolved_q_limit_options(io::IO, metadata::AbstractDict)
  println(io, "Resolved Q-limit options")
  println(io, "------------------------")
  println(io, "Q-limit handling enabled : ", metadata["qlimits_enabled"])
  println(io, "Q-limit enforcement mode : ", get(metadata, "qlimit_enforcement_mode", "active_set"))
  println(io, "Q-limit guard enabled    : ", metadata["qlimit_guard_enabled"])
  println(io, "Q-limit preview mode     : ", metadata["q_limit_preview_mode"])
  println(io, "Q-limit runlog max rows  : ", metadata["q_limit_runlog_max_rows"])
  println(io, "Q-limit detail artifact  : ", Q_LIMIT_LOG_ARTIFACT)
  println(io)
  return nothing
end

function _profile_elapsed_call_tuple(value)
  if value isa NamedTuple && haskey(value, :elapsed_s)
    return (elapsed = Float64(value.elapsed_s), calls = Int(get(value, :calls, 1)))
  elseif value isa AbstractDict && haskey(value, :elapsed_s)
    return (elapsed = Float64(value[:elapsed_s]), calls = Int(get(value, :calls, 1)))
  elseif value isa Number
    return (elapsed = Float64(value), calls = 1)
  end
  return nothing
end

function _profile_aggregate(profile::AbstractDict, keys)::NamedTuple
  total = 0.0
  calls = 0
  max_elapsed = 0.0
  for key in keys
    haskey(profile, key) || continue
    item = _profile_elapsed_call_tuple(profile[key])
    item === nothing && continue
    total += item.elapsed
    calls += item.calls
    max_elapsed = max(max_elapsed, item.elapsed)
  end
  return (total = total, calls = calls, max = max_elapsed)
end

function _write_compact_internal_timing_aggregates(io::IO, result::SparlectraRunResult)
  result.performance_profile isa AbstractDict || return nothing
  profile = result.performance_profile
  aggregates = (
    building_ybus = _profile_aggregate(profile, (:ybus_assembly, :ybus_expand_isolated)),
    solver_initialization = _profile_aggregate(profile, (:solver_initial_voltage, :solver_initial_injections)),
    start_projection = _profile_aggregate(profile, (:start_projection_dc_linear_solve, :start_projection_apply)),
    newton_iteration = _profile_aggregate(profile, (:iteration_controlled_injections, :iteration_mismatch, :iteration_newton_step)),
    q_limit_processing = _profile_aggregate(profile, (:qlimit_iteration, :qlimit_generation_update)),
    linear_solve = _profile_aggregate(profile, (:newton_step_linear_solve, :start_projection_dc_linear_solve)),
  )
  any(item -> item.calls > 0, values(aggregates)) || return nothing
  println(io)
  println(io, "Compact internal timing aggregates")
  println(io, "----------------------------------")
  for (name, item) in pairs(aggregates)
    item.calls > 0 || continue
    println(io, name, "_count: ", item.calls)
    println(io, name, "_total_seconds: ", round(item.total; digits = 6))
    name === :linear_solve && println(io, "linear_solve_max_seconds: ", round(item.max; digits = 6))
  end
  return nothing
end

function _write_performance_log(path::AbstractString, mode::Symbol, phases::AbstractDict, result::SparlectraRunResult, phase_timings::AbstractVector = Dict{String,Any}[])
  open(path, "w") do io
    println(io, "Sparlectra single-run phase timing")
    println(io, "==================================")
    println(io, "mode: ", mode)
    println(io, "This file describes one API/Web UI run; benchmark mode measures repeated solves.")
    println(io)
    for phase in (:request_parse, :case_resolution, :api_config_build, :case_loading_network_solver, :solver, :postprocessing, :artifact_writing, :total)
      haskey(phases, phase) || continue
      println(io, rpad(String(phase) * ":", 31), round(Float64(phases[phase]); digits = 6), " s")
    end
    if !isempty(phase_timings)
      println(io)
      _write_service_phase_summary(io, phase_timings)
    end
    _write_compact_internal_timing_aggregates(io, result)
    if mode === :full && result.performance_profile isa AbstractDict
      println(io)
      println(io, "Available internal performance profile")
      println(io, "--------------------------------------")
      for key in sort!(collect(keys(result.performance_profile)); by = string)
        println(io, key, ": ", result.performance_profile[key])
      end
    end
  end
  return path
end

function _resolve_detailed_csv_format(value)::NamedTuple
  name = String(value)
  name == "technical" && return (name = name, delimiter = ',', decimal_separator = '.', thousands_separator = "")
  name == "excel_de" && return (name = name, delimiter = ';', decimal_separator = ',', thousands_separator = ".")
  name == "excel_us" && return (name = name, delimiter = ',', decimal_separator = '.', thousands_separator = ",")
  throw(ArgumentError("Unsupported detailed_result_csv_format \"$(name)\". Expected technical, excel_de, or excel_us."))
end

function _group_csv_integer(text::AbstractString, separator::AbstractString)::String
  isempty(separator) && return String(text)

  sign = startswith(text, "-") ? "-" : ""
  digits = isempty(sign) ? String(text) : text[2:end]

  # Thousands grouping is only valid for integer digit strings. Leave textual
  # values such as Bool strings unchanged if they accidentally reach this helper.
  (isempty(digits) || !all(isdigit, digits)) && return String(text)

  first_group = mod(length(digits), 3)
  first_group == 0 && (first_group = 3)
  groups = String[digits[1:first_group]]
  for start = (first_group+1):3:length(digits)
    push!(groups, digits[start:(start+2)])
  end
  return sign * join(groups, separator)
end

function _format_csv_number(value::Integer, format)::String
  return _group_csv_integer(string(value), format.thousands_separator)
end

function _format_csv_number(value::AbstractFloat, format)::String
  isnan(value) && return "NaN"
  isinf(value) && return signbit(value) ? "-Inf" : "Inf"
  technical = @sprintf("%.15g", value)
  if format.name != "technical"
    exponent_marker = findfirst(character -> character in ('e', 'E'), technical)
    if exponent_marker !== nothing
      mantissa = technical[begin:prevind(technical, exponent_marker)]
      exponent = parse(Int, technical[nextind(technical, exponent_marker):end])
      sign = startswith(mantissa, "-") ? "-" : ""
      unsigned = isempty(sign) ? mantissa : mantissa[2:end]
      dot_index = findfirst(==('.'), unsigned)
      fractional_digits = dot_index === nothing ? 0 : ncodeunits(unsigned) - dot_index
      digits = replace(unsigned, "." => "")
      decimal_position = ncodeunits(digits) - fractional_digits + exponent
      if decimal_position <= 0
        technical = sign * "0." * repeat("0", -decimal_position) * digits
      elseif decimal_position >= ncodeunits(digits)
        technical = sign * digits * repeat("0", decimal_position - ncodeunits(digits))
      else
        technical = sign * digits[1:decimal_position] * "." * digits[(decimal_position+1):end]
      end
    end
    dot_index = findfirst(==('.'), technical)
    if dot_index !== nothing
      last_nonzero = lastindex(technical)
      while last_nonzero > dot_index && technical[last_nonzero] == '0'
        last_nonzero = prevind(technical, last_nonzero)
      end
      technical = last_nonzero == dot_index ? technical[begin:prevind(technical, dot_index)] : technical[begin:last_nonzero]
      isempty(technical) && (technical = "0")
      technical == "-0" && (technical = "0")
    end
  end
  exponent_marker = format.name == "technical" ? findfirst(character -> character in ('e', 'E'), technical) : nothing
  mantissa = exponent_marker === nothing ? technical : technical[begin:prevind(technical, exponent_marker)]
  exponent = exponent_marker === nothing ? "" : technical[exponent_marker:end]
  dot_index = findfirst(==('.'), mantissa)
  integer_text = dot_index === nothing ? mantissa : mantissa[begin:prevind(mantissa, dot_index)]
  integer_part = _group_csv_integer(integer_text, format.thousands_separator)
  fractional_part = dot_index === nothing ? "" : string(format.decimal_separator, mantissa[nextind(mantissa, dot_index):end])
  return integer_part * fractional_part * exponent
end

struct CsvFormatRuntime
  name::String
  delimiter::Char
  decimal_separator::Char
  thousands_separator::String
end

CsvFormatRuntime(format) = CsvFormatRuntime(String(format.name), format.delimiter, format.decimal_separator, String(format.thousands_separator))

function _csv_needs_quotes(text::AbstractString, delimiter::Char)::Bool
  for character in text
    character in (delimiter, '"', '\r', '\n') && return true
  end
  return false
end

function write_csv_cell!(io::IO, value, delimiter::Char, fmt::CsvFormatRuntime)
  value === missing && return nothing
  value === nothing && return nothing
  text = value isa Bool ? string(value) : value isa Integer ? _format_csv_number(value, fmt) : value isa AbstractFloat ? _format_csv_number(value, fmt) : string(value)
  if _csv_needs_quotes(text, delimiter)
    print(io, '"')
    for character in text
      character == '"' && print(io, '"')
      print(io, character)
    end
    print(io, '"')
  else
    print(io, text)
  end
  return nothing
end

function write_csv_row_direct!(io::IO, delimiter::Char, fmt::CsvFormatRuntime, values...)
  first = true
  for value in values
    first || print(io, delimiter)
    write_csv_cell!(io, value, delimiter, fmt)
    first = false
  end
  println(io)
  return nothing
end

function _format_csv_value(value, format)::String
  value === missing && return ""
  value === nothing && return ""

  # Bool is an Integer subtype in Julia. Format it before Integer values so
  # Excel-oriented thousands grouping cannot turn true/false into t.rue/fa.lse.
  value isa Bool && return string(value)

  value isa Integer && return _format_csv_number(value, format)
  value isa AbstractFloat && return _format_csv_number(value, format)
  return string(value)
end

function _csv_field(value, delimiter::Char, format = _resolve_detailed_csv_format("technical"))::String
  value === missing && return ""
  value === nothing && return ""
  text = _format_csv_value(value, format)
  if any(character -> character in (delimiter, '"', '\r', '\n'), text)
    return "\"" * replace(text, "\"" => "\"\"") * "\""
  end
  return text
end

const _DETAILED_CSV_BUFFER_INITIAL_BYTES_DEFAULT = 8 * 1024 * 1024
const _DETAILED_CSV_BUFFER_MAX_BYTES_DEFAULT = 64 * 1024 * 1024
const _DETAILED_CSV_STREAMING_THRESHOLD_ROWS_DEFAULT = 100_000
const _DETAILED_CSV_DIRECT_THRESHOLD_BUSES_DEFAULT = 10_000

function _csv_write_options(config = nothing)::NamedTuple
  output = config isa SparlectraConfig ? config.output : config isa OutputConfig ? config : nothing
  mode = output === nothing ? :auto : output.detailed_result_csv_write_mode
  mode in OUTPUT_DETAILED_RESULT_CSV_WRITE_MODE_VALUES || throw(ArgumentError("Unsupported detailed_result_csv_write_mode \"$(mode)\". Expected auto, buffered, or streaming."))
  initial_bytes = output === nothing ? _DETAILED_CSV_BUFFER_INITIAL_BYTES_DEFAULT : output.detailed_result_csv_buffer_initial_bytes
  max_bytes = output === nothing ? _DETAILED_CSV_BUFFER_MAX_BYTES_DEFAULT : output.detailed_result_csv_buffer_max_bytes
  threshold_rows = output === nothing ? _DETAILED_CSV_STREAMING_THRESHOLD_ROWS_DEFAULT : output.detailed_result_csv_streaming_threshold_rows
  initial_bytes = initial_bytes < 0 ? _DETAILED_CSV_BUFFER_INITIAL_BYTES_DEFAULT : initial_bytes
  max_bytes = max_bytes <= 0 ? _DETAILED_CSV_BUFFER_MAX_BYTES_DEFAULT : max_bytes
  threshold_rows = threshold_rows <= 0 ? _DETAILED_CSV_STREAMING_THRESHOLD_ROWS_DEFAULT : threshold_rows
  return (mode = mode, initial_bytes = initial_bytes, max_bytes = max_bytes, threshold_rows = threshold_rows)
end

function _estimated_namedtuple_csv_bytes(rows::AbstractVector, columns)::Int
  return length(join(String.(columns), ',')) + 1 + length(rows) * max(1, length(columns)) * 24
end

function _write_namedtuple_csv_row(io::IO, row, columns, delimiter::Char, resolved_format)
  println(io, join((_csv_field(getproperty(row, column), delimiter, resolved_format) for column in columns), delimiter))
  return nothing
end

function _write_csv_row_values(io::IO, values, delimiter::Char, resolved_format)
  first = true
  for value in values
    first || print(io, delimiter)
    print(io, _csv_field(value, delimiter, resolved_format))
    first = false
  end
  println(io)
  return nothing
end

function _write_namedtuple_csv_buffered(path::AbstractString, rows::AbstractVector, columns, delimiter::Char, resolved_format; initial_bytes::Integer = _DETAILED_CSV_BUFFER_INITIAL_BYTES_DEFAULT)
  buffer = IOBuffer(; sizehint = max(0, Int(initial_bytes)))
  println(buffer, join(String.(columns), delimiter))
  for row in rows
    _write_namedtuple_csv_row(buffer, row, columns, delimiter, resolved_format)
  end
  open(path, "w") do io
    write(io, take!(buffer))
  end
  return path
end

function _write_namedtuple_csv_streaming(path::AbstractString, rows::AbstractVector, columns, delimiter::Char, resolved_format)
  open(path, "w") do io
    println(io, join(String.(columns), delimiter))
    for row in rows
      _write_namedtuple_csv_row(io, row, columns, delimiter, resolved_format)
    end
  end
  return path
end

function _select_namedtuple_csv_write_mode(rows::AbstractVector, columns; config = nothing, estimated_rows::Integer = length(rows))::Symbol
  options = _csv_write_options(config)
  options.mode === :buffered && return :buffered
  options.mode === :streaming && return :streaming
  estimated_rows > options.threshold_rows && return :streaming
  _estimated_namedtuple_csv_bytes(rows, columns) > options.max_bytes && return :streaming
  return :buffered
end

function _write_namedtuple_csv(path::AbstractString, rows::AbstractVector, columns; delimiter::Char = ',', format = nothing, config = nothing, estimated_rows::Integer = length(rows))
  delimiter in (',', ';') || throw(ArgumentError("CSV delimiter must be ',' or ';'."))
  resolved_format = format === nothing ? (name = "custom", delimiter = delimiter, decimal_separator = '.', thousands_separator = "") : _resolve_detailed_csv_format(format)
  resolved_format.delimiter == delimiter || throw(ArgumentError("CSV delimiter does not match detailed CSV format $(resolved_format.name)."))
  options = _csv_write_options(config)
  mode = _select_namedtuple_csv_write_mode(rows, columns; config, estimated_rows)
  mode === :streaming && return _write_namedtuple_csv_streaming(path, rows, columns, delimiter, resolved_format)
  return _write_namedtuple_csv_buffered(path, rows, columns, delimiter, resolved_format; initial_bytes = min(options.initial_bytes, options.max_bytes))
end

function _detailed_csv_export_options(config = nothing)::NamedTuple
  output = config isa SparlectraConfig ? config.output : config isa OutputConfig ? config : nothing
  exporter = output === nothing ? :auto : output.detailed_result_csv_exporter
  threshold_buses = output === nothing ? _DETAILED_CSV_DIRECT_THRESHOLD_BUSES_DEFAULT : output.detailed_result_csv_direct_threshold_buses
  exporter in OUTPUT_DETAILED_RESULT_CSV_EXPORTER_VALUES || throw(ArgumentError("Unsupported detailed_result_csv_exporter \"$(exporter)\". Expected auto, report, or direct."))
  return (exporter = exporter, threshold_buses = threshold_buses <= 0 ? _DETAILED_CSV_DIRECT_THRESHOLD_BUSES_DEFAULT : threshold_buses)
end

function _select_detailed_csv_exporter(net::Net; config = nothing)::Symbol
  options = _detailed_csv_export_options(config)
  options.exporter === :report && return :report
  options.exporter === :direct && return :direct
  return length(net.nodeVec) >= options.threshold_buses ? :direct : :report
end

function _complex_voltage_rows(node_rows::AbstractVector, format)
  return [
    begin
      angle = deg2rad(Float64(row.va_deg))
      v_re = Float64(row.vm_pu) * cos(angle)
      v_im = Float64(row.vm_pu) * sin(angle)
      merge(row, (v_re = round(v_re; sigdigits = 10), v_im = round(v_im; sigdigits = 10), v_complex = string(_format_csv_number(v_re, format), signbit(v_im) ? " - j" : " + j", _format_csv_number(abs(v_im), format))))
    end for row in node_rows
  ]
end

function _branch_csv_values(br)
  f = br.fBranchFlow
  t = br.tBranchFlow
  p_from = isnothing(f) || isnothing(f.pFlow) ? 0.0 : f.pFlow
  q_from = isnothing(f) || isnothing(f.qFlow) ? 0.0 : f.qFlow
  p_to = isnothing(t) || isnothing(t.pFlow) ? 0.0 : t.pFlow
  q_to = isnothing(t) || isnothing(t.qFlow) ? 0.0 : t.qFlow
  rated = isnothing(br.sn_MVA) ? 0.0 : br.sn_MVA
  overloaded = rated > 0.0 && max(abs(p_from), abs(p_to)) > rated
  return (br.comp.cName, br.branchIdx, br.fromBus, br.toBus, br.status, p_from, q_from, p_to, q_to, _default0(br.pLosses), _default0(br.qLosses), rated, overloaded)
end

function _check_csv_abort(abort_checker, row_count::Integer)
  if abort_checker !== nothing && row_count % 1000 == 0
    abort_checker()
  end
  return nothing
end

function _write_detailed_result_csv_artifacts_direct!(artifacts::Vector{String}, output_path::AbstractString, net::Net, result::SparlectraRunResult, config; format = "technical", abort_checker = nothing, timing_metadata = nothing)::Vector{String}
  resolved_format = _resolve_detailed_csv_format(format)
  runtime_format = CsvFormatRuntime(resolved_format)
  bus_columns = (:bus, :bus_name, :type, :vm_pu, :va_deg, :vn_kV, :v_re, :v_im, :v_complex, :v_kV, :p_gen_MW, :q_gen_MVar, :p_load_MW, :q_load_MVar, :q_limit_hit, :control)
  branch_columns = (:branch, :branch_index, :from_bus, :to_bus, :status, :p_from_MW, :q_from_MVar, :p_to_MW, :q_to_MVar, :p_loss_MW, :q_loss_MVar, :rated_MVA, :overloaded)
  busNameByIdx = _bus_name_by_idx(net)
  power_components = _bus_power_component_cache(net)
  cache_start = time_ns()
  control_labels = _bus_control_flag_cache(net)
  control_cache_elapsed = _api_elapsed_seconds(cache_start)
  nodes_sorted = sort(net.nodeVec, by = x -> x.busIdx)
  abort_checker === nothing || abort_checker()
  total_start = time_ns()
  bus_rows = 0
  branch_rows = 0
  current_artifact = artifacts[1]
  current_rows = Ref(0)
  progress_callback = timing_metadata isa AbstractDict ? get(timing_metadata, :progress_callback, nothing) : nothing
  _emit_csv_progress(event::AbstractString; fields...) = progress_callback === nothing ? nothing : progress_callback(event; fields...)
  bus_elapsed = 0.0
  branch_elapsed = 0.0
  _emit_csv_progress(
    "detailed_result_csv_export_started";
    exporter = "direct",
    write_mode = "streaming",
    csv_format = resolved_format.name,
    bus_rows = length(nodes_sorted),
    branch_rows = length(net.branchVec),
  )
  bus_start = time_ns()

  try
    open(joinpath(output_path, artifacts[1]), "w") do io
      println(io, join(String.(bus_columns), resolved_format.delimiter))
      row_count = 0
      for n in nodes_sorted
        row_count += 1
        current_rows[] = row_count
      angle = deg2rad(Float64(n._va_deg))
      v_re = Float64(n._vm_pu) * cos(angle)
      v_im = Float64(n._vm_pu) * sin(angle)
      p_gen, q_gen, p_load, q_load = _bus_power_components(power_components, n.busIdx)
      v_re_rounded = round(v_re; sigdigits = 10)
      v_im_rounded = round(v_im; sigdigits = 10)
      v_complex = string(_format_csv_number(v_re, runtime_format), signbit(v_im) ? " - j" : " + j", _format_csv_number(abs(v_im), runtime_format))
      write_csv_row_direct!(
        io,
        resolved_format.delimiter,
        runtime_format,
        n.busIdx,
        get(busNameByIdx, n.busIdx, n.comp.cName),
        toString(n._nodeType),
        n._vm_pu,
        n._va_deg,
        n.comp.cVN,
        v_re_rounded,
        v_im_rounded,
        v_complex,
        n.comp.cVN * n._vm_pu,
        p_gen,
        q_gen,
        p_load,
        q_load,
        haskey(net.qLimitEvents, n.busIdx),
        _cached_control_label(control_labels, n.busIdx),
      )
      _check_csv_abort(abort_checker, row_count)
      end
      bus_rows = row_count
    end
    bus_elapsed = _api_elapsed_seconds(bus_start)
    _emit_csv_progress("detailed_result_csv_file_written"; artifact = artifacts[1], rows = bus_rows, elapsed_s = bus_elapsed)

    try
      current_artifact = artifacts[2]
      current_rows[] = 0
      branch_start = time_ns()
      open(joinpath(output_path, artifacts[2]), "w") do io
        println(io, join(String.(branch_columns), resolved_format.delimiter))
        row_count = 0
        for br in net.branchVec
          row_count += 1
          current_rows[] = row_count
          f = br.fBranchFlow
          t = br.tBranchFlow
          p_from = isnothing(f) || isnothing(f.pFlow) ? 0.0 : f.pFlow
          q_from = isnothing(f) || isnothing(f.qFlow) ? 0.0 : f.qFlow
          p_to = isnothing(t) || isnothing(t.pFlow) ? 0.0 : t.pFlow
          q_to = isnothing(t) || isnothing(t.qFlow) ? 0.0 : t.qFlow
          rated = isnothing(br.sn_MVA) ? 0.0 : br.sn_MVA
          overloaded = rated > 0.0 && max(abs(p_from), abs(p_to)) > rated
          write_csv_row_direct!(
            io,
            resolved_format.delimiter,
            runtime_format,
            br.comp.cName,
            br.branchIdx,
            br.fromBus,
            br.toBus,
            br.status,
            p_from,
            q_from,
            p_to,
            q_to,
            _default0(br.pLosses),
            _default0(br.qLosses),
            rated,
            overloaded,
          )
          _check_csv_abort(abort_checker, row_count)
        end
        branch_rows = row_count
      end
      branch_elapsed = _api_elapsed_seconds(branch_start)
      _emit_csv_progress("detailed_result_csv_file_written"; artifact = artifacts[2], rows = branch_rows, elapsed_s = branch_elapsed)
    catch err
      timing_metadata !== nothing && (timing_metadata[:partial_error] = "branch_flows.csv failed: $(sprint(showerror, err))")
      _emit_csv_progress("detailed_result_csv_export_partial"; artifact = artifacts[2], rows_written = current_rows[], elapsed_s = _api_elapsed_seconds(total_start), error = timing_metadata === nothing ? sprint(showerror, err) : timing_metadata[:partial_error])
      return [artifacts[1]]
    end
  catch err
    _emit_csv_progress("detailed_result_csv_export_aborted"; artifact = current_artifact, rows_written = current_rows[], elapsed_s = _api_elapsed_seconds(total_start))
    rethrow()
  end
  if timing_metadata !== nothing
    timing_metadata[:exporter] = :direct
    timing_metadata[:write_mode] = :streaming
    timing_metadata[:format] = resolved_format.name
    timing_metadata[:formatting_mode] = :direct_low_allocation
    timing_metadata[:control_label_cache_s] = control_cache_elapsed
    timing_metadata[:bus_rows] = bus_rows
    timing_metadata[:branch_rows] = branch_rows
    timing_metadata[:bus_export_s] = bus_elapsed
    timing_metadata[:branch_export_s] = branch_elapsed
    timing_metadata[:total_export_s] = _api_elapsed_seconds(total_start)
  end
  return artifacts
end

function _write_detailed_result_csv(output_path::AbstractString, result::SparlectraRunResult; format = "technical", config = nothing, abort_checker = nothing, timing_metadata = nothing)::Vector{String}
  resolved_format = _resolve_detailed_csv_format(format)
  result.net === nothing && throw(ArgumentError("PowerFlow result does not contain a network for detailed CSV export."))
  artifacts = ["bus_voltages_complex.csv", "branch_flows.csv"]
  exporter = _select_detailed_csv_exporter(result.net; config)
  exporter === :direct && return _write_detailed_result_csv_artifacts_direct!(artifacts, output_path, result.net, result, config; format = resolved_format.name, abort_checker, timing_metadata)
  report = buildACPFlowReport(result.net; ct = result.elapsed_s, ite = result.iterations, converged = result.final_converged)
  bus_columns = (:bus, :bus_name, :type, :vm_pu, :va_deg, :vn_kV, :v_re, :v_im, :v_complex, :v_kV, :p_gen_MW, :q_gen_MVar, :p_load_MW, :q_load_MVar, :q_limit_hit, :control)
  branch_columns = (:branch, :branch_index, :from_bus, :to_bus, :status, :p_from_MW, :q_from_MVar, :p_to_MW, :q_to_MVar, :p_loss_MW, :q_loss_MVar, :rated_MVA, :overloaded)
  estimated_rows = length(report.nodes) + length(report.branches)
  _write_namedtuple_csv(joinpath(output_path, artifacts[1]), _complex_voltage_rows(report.nodes, resolved_format), bus_columns; delimiter = resolved_format.delimiter, format = resolved_format.name, config, estimated_rows)
  try
    _write_namedtuple_csv(joinpath(output_path, artifacts[2]), report.branches, branch_columns; delimiter = resolved_format.delimiter, format = resolved_format.name, config, estimated_rows)
  catch err
    timing_metadata !== nothing && (timing_metadata[:partial_error] = "branch_flows.csv failed: $(sprint(showerror, err))")
    return [artifacts[1]]
  end
  return artifacts
end

function _csv_solution_quality(raw_result::SparlectraRunResult)::String
  return raw_result.final_converged && raw_result.solution_available ? "converged" : "not_converged_last_iterate"
end

function _write_numerical_outcome_summary(io::IO, raw_result::SparlectraRunResult)
  raw_result.final_converged && raw_result.solution_available && return nothing
  println(io)
  println(io, "Numerical outcome: not_converged")
  println(io, "Primary reason: ", raw_result.reason_text)
  println(io, "Final mismatch: ", isfinite(raw_result.final_mismatch) ? raw_result.final_mismatch : (isnan(raw_result.final_mismatch) ? "NaN" : "unavailable"))
  println(io, "Final mismatch status: ", _final_mismatch_status(raw_result))
  qlimit_ok = hasproperty(raw_result.diagnostics, :q_limit_active_set_ok) ? getproperty(raw_result.diagnostics, :q_limit_active_set_ok) : nothing
  qlimit_ok === nothing || println(io, "Q-limit active-set convergence: ", qlimit_ok ? "yes" : "no")
  for (label, field) in (
    ("PV->PQ switching events", :pv_pq_switching_events),
    ("Q-limit active-set changes", :qlimit_active_set_changes),
    ("Q-limit re-enable events", :qlimit_reenable_events),
    ("Guarded narrow-Q PV buses", :guarded_narrow_q_pv_buses),
  )
    hasproperty(raw_result.diagnostics, field) && println(io, label, ": ", getproperty(raw_result.diagnostics, field))
  end
  return nothing
end

function _update_effective_matpower_raw!(raw::Dict{String,Any}, cfg::SparlectraConfig)
  mat_raw = get!(raw, "matpower_import", Dict{String,Any}())
  mat_raw isa Dict{String,Any} || return raw
  mat_raw["auto_profile"] = String(cfg.matpower.auto_profile)
  mat_raw["ratio"] = String(cfg.matpower.ratio)
  mat_raw["shift_sign"] = cfg.matpower.shift_sign
  mat_raw["shift_unit"] = String(cfg.matpower.shift_unit)
  mat_raw["bus_shunt_model"] = String(cfg.matpower.bus_shunt_model)
  mat_raw["pv_voltage_source"] = String(cfg.matpower.pv_voltage_source)
  mat_raw["compare_voltage_reference"] = String(cfg.matpower.compare_voltage_reference)
  return raw
end

function _write_matpower_auto_profile_artifact(output_path::AbstractString, profile, cfg::SparlectraConfig; casefile::AbstractString)::String
  artifact = joinpath(output_path, "matpower_auto_profile.log")
  open(artifact, "w") do io
    println(io, "MATPOWER import auto-profile artifact")
    println(io, "====================================")
    write_matpower_import_auto_profile(io, profile, cfg; casefile = casefile)
  end
  return artifact
end

function _final_outcome_payload(raw_result::SparlectraRunResult)::Dict{String,Any}
  return Dict{String,Any}(
    "converged" => raw_result.final_converged,
    "numerical_converged" => raw_result.numerical_converged,
    "solution_available" => raw_result.solution_available,
    "iterations" => raw_result.iterations,
    "final_mismatch" => isfinite(raw_result.final_mismatch) ? raw_result.final_mismatch : nothing,
    "reason" => String(raw_result.reason),
  )
end

function _runtime_metadata_payload(; case_path::AbstractString, lifecycle::AbstractDict = Dict{String,Any}())
  payload = Dict{String,Any}(
    "runtime_request" => Dict{String,Any}(
      "casefile" => basename(case_path),
      "casefile_path" => String(case_path),
    ),
  )
  haskey(lifecycle, "configured_default_casefile") && (payload["configured_default_casefile"] = lifecycle["configured_default_casefile"])
  haskey(lifecycle, "runtime_casefile") && (payload["runtime_casefile"] = lifecycle["runtime_casefile"])
  for key in ("solver_status", "artifact_status", "run_status", "last_phase", "last_heartbeat", "final_outcome")
    haskey(lifecycle, key) && (payload[key] = lifecycle[key])
  end
  return payload
end

function _write_run_metadata_artifact(output_path::AbstractString; case_path::AbstractString, lifecycle::AbstractDict = Dict{String,Any}())::String
  artifact = joinpath(output_path, "run_metadata.yaml")
  _write_yaml_file(artifact, _runtime_metadata_payload(case_path = case_path, lifecycle = lifecycle))
  return artifact
end

function _api_result(;
  run_id::String = string(uuid4()),
  schema_version::String = _SPARLECTRA_API_SCHEMA_VERSION,
  status::Symbol,
  success::Bool,
  converged = nothing,
  solution_available::Bool = false,
  iterations = nothing,
  final_mismatch = nothing,
  reason = nothing,
  message = nothing,
  casefile = nothing,
  config_file = nothing,
  output_dir::String,
  logfile = nothing,
  result_file = nothing,
  artifacts = SparlectraApiArtifact[],
  service_phase_timings = Dict{String,Any}[],
  metadata = Dict{String,Any}(),
  raw_result = nothing,
)
  timings = Dict{String,Any}[Dict{String,Any}(String(key) => value for (key, value) in timing) for timing in service_phase_timings]
  return SparlectraApiResult(run_id, schema_version, status, success, converged, solution_available, iterations, final_mismatch, reason, message, casefile, config_file, output_dir, logfile, result_file, artifacts, timings, Dict{String,Any}(String(key) => value for (key, value) in metadata), raw_result)
end

function _write_api_result_file(result::SparlectraApiResult)
  result.result_file === nothing || write(result.result_file, to_json(result), '\n')
  return result
end

function _refresh_api_artifacts(result::SparlectraApiResult)::SparlectraApiResult
  return _api_result(
    run_id = result.run_id,
    schema_version = result.schema_version,
    status = result.status,
    success = result.success,
    converged = result.converged,
    solution_available = result.solution_available,
    iterations = result.iterations,
    final_mismatch = result.final_mismatch,
    reason = result.reason,
    message = result.message,
    casefile = result.casefile,
    config_file = result.config_file,
    output_dir = result.output_dir,
    logfile = result.logfile,
    result_file = result.result_file,
    artifacts = collect_sparlectra_api_artifacts(result.output_dir),
    service_phase_timings = result.service_phase_timings,
    metadata = result.metadata,
    raw_result = result.raw_result,
  )
end

function _finalize_api_result(result::SparlectraApiResult)::SparlectraApiResult
  _write_api_result_file(result)
  refreshed = _refresh_api_artifacts(result)
  _write_api_result_file(refreshed)
  refreshed = _refresh_api_artifacts(refreshed)
  _write_api_result_file(refreshed)
  return refreshed
end

function _api_failure(reason::String, message::String; run_id::String = string(uuid4()), casefile, config_file, output_dir::String, logfile::String, result_file::String, service_phase_timings = Dict{String,Any}[], metadata = Dict{String,Any}())::SparlectraApiResult
  open(logfile, "a") do io
    println(io, "Sparlectra API failure: ", reason)
    println(io, message)
  end
  failure_metadata = merge(Dict{String,Any}("failure_reason" => reason), Dict{String,Any}(String(key) => value for (key, value) in metadata))
  result = _api_result(run_id = run_id, status = :failed, success = false, reason = reason, message = message, casefile = casefile, config_file = config_file, output_dir = output_dir, logfile = logfile, result_file = result_file, service_phase_timings = service_phase_timings, metadata = failure_metadata)
  return _finalize_api_result(result)
end

function _api_execution_failure(reason::String, message::String; run_id::String, casefile, config_file, output_dir::String, logfile::String, result_file::String, phase_recorder::PowerFlowPhaseTimingRecorder, performance_timing = :off, metadata = Dict{String,Any}())::SparlectraApiResult
  _complete_active_phase!(phase_recorder, "failed")
  _start_service_phase!(phase_recorder, "finalizing_failed")
  _complete_active_phase!(phase_recorder, "failed")
  timing_mode = try
    _api_timing_mode(performance_timing)
  catch
    :off
  end
  if timing_mode !== :off
    open(joinpath(output_dir, "performance.log"), "a") do io
      println(io, "Sparlectra API failure: ", reason)
      println(io, message)
      println(io)
      _write_service_phase_summary(io, phase_recorder.timings)
    end
  end
  return _api_failure(reason, message; run_id, casefile, config_file, output_dir, logfile, result_file, service_phase_timings = phase_recorder.timings, metadata = metadata)
end

"""
    run_sparlectra_api(; casefile, config_file, output_dir, config_overrides=Dict(),
                        performance_timing=:off, run_diagnostics=false,
                        detailed_result_csv=false,
                        detailed_result_csv_format="technical",
                        detailed_result_csv_semicolon=false) -> SparlectraApiResult

Run one MATPOWER power-flow case through a stable, non-interactive API contract.
The function validates GUI overrides, writes `effective_config.yaml`, delegates
the numerical work to [`run_sparlectra`](@ref), captures textual output in
`run.log`, writes `result.json`, discovers all generated files, and returns
structured status and artifact metadata. The input configuration template is
never modified. `performance_timing` may be `:off`, `:compact`, or `:full` and
writes a single-run `performance.log`; `run_diagnostics=true` captures existing
PowerFlow diagnostic printers in `diagnose.log`; and `detailed_result_csv=true`
writes Excel-friendly bus-voltage and branch-flow CSV artifacts. Optional
`detailed_result_csv_format` accepts `technical`, `excel_de`, or `excel_us`.
The legacy `detailed_result_csv_semicolon=true` maps to `excel_de` when the
explicit format is omitted. Artifact generation does not change PowerFlow run
success.
"""
function run_sparlectra_api(;
  casefile::AbstractString,
  config_file::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH,
  output_dir::AbstractString,
  config_overrides::AbstractDict = Dict{String,Any}(),
  performance_timing = :off,
  run_diagnostics::Bool = false,
  detailed_result_csv::Bool = false,
  detailed_result_csv_format = nothing,
  detailed_result_csv_semicolon::Bool = false,
)::SparlectraApiResult
  return _run_sparlectra_api(
    casefile = casefile,
    config_file = config_file,
    output_dir = output_dir,
    config_overrides = config_overrides,
    performance_timing = performance_timing,
    run_diagnostics = run_diagnostics,
    detailed_result_csv = detailed_result_csv,
    detailed_result_csv_format = detailed_result_csv_format,
    detailed_result_csv_semicolon = detailed_result_csv_semicolon,
    run_id = string(uuid4()),
  )
end

function _run_sparlectra_api(;
  casefile::AbstractString,
  config_file::AbstractString,
  output_dir::AbstractString,
  config_overrides::AbstractDict,
  performance_timing = :off,
  run_diagnostics::Bool = false,
  detailed_result_csv::Bool = false,
  detailed_result_csv_format = nothing,
  detailed_result_csv_semicolon::Bool = false,
  phase_timings::AbstractDict = Dict{Symbol,Float64}(),
  run_id::String,
  cancellation_token = nothing,
  phase_callback = phase -> nothing,
  operation_callback = (event; fields...) -> nothing,
)::SparlectraApiResult
  total_start = time_ns()
  phases = Dict{Symbol,Float64}(Symbol(key) => Float64(value) for (key, value) in phase_timings)
  phase_recorder = PowerFlowPhaseTimingRecorder()
  emit_phase = phase -> begin
    _start_service_phase!(phase_recorder, String(phase))
    phase_callback(String(phase))
  end
  emit_phase("total_service")
  timing_mode = try
    _api_timing_mode(performance_timing)
  catch err
    output_path = abspath(output_dir)
    mkpath(output_path)
    logfile = joinpath(output_path, "run.log")
    result_file = joinpath(output_path, "result.json")
    touch(logfile)
    return _api_failure("invalid_performance_timing", sprint(showerror, err); run_id, casefile = abspath(casefile), config_file = abspath(config_file), output_dir = output_path, logfile, result_file)
  end
  output_path = abspath(output_dir)
  mkpath(output_path)
  case_path = abspath(casefile)
  config_path = abspath(config_file)
  logfile = joinpath(output_path, "run.log")
  result_file = joinpath(output_path, "result.json")
  touch(logfile)
  emit_phase("preparing_configuration")
  _check_powerflow_cancelled!(cancellation_token)

  isfile(config_path) || return _api_failure("config_file_not_found", "Configuration file not found: $(config_path)"; run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)

  config_start = time_ns()
  nested_overrides = try
    validate_gui_config_overrides(config_overrides)
  catch err
    return _api_failure("invalid_config_override", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)
  end

  config = nothing
  effective_raw = nothing
  try
    config, effective_raw = _load_api_config(config_path, nested_overrides)
  catch err
    return _api_failure("invalid_configuration", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)
  end
  effective_config = joinpath(output_path, "effective_config.yaml")
  _write_yaml_file(effective_config, effective_raw)
  _write_run_metadata_artifact(output_path; case_path = case_path)
  _check_powerflow_cancelled!(cancellation_token)
  phases[:api_config_build] = _api_elapsed_seconds(config_start)
  isfile(case_path) || return _api_failure("casefile_not_found", "MATPOWER case file not found: $(case_path)"; run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)

  raw_result = nothing
  qlimit_metadata = _resolved_q_limit_runtime_options(config)
  qlimit_metadata["runtime_casefile"] = basename(case_path)
  qlimit_metadata["runtime_casefile_path"] = case_path
  qlimit_metadata["configured_default_casefile"] = config.matpower.case
  matpower_metadata = _resolved_matpower_import_runtime_options(config)
  operation_callback("powerflow_effective_options"; run_id = run_id, case = basename(case_path), _metadata_kwargs(qlimit_metadata)..., _metadata_kwargs(matpower_metadata)...)
  emit_phase("reading_matpower_case")
  api_performance_profile = Dict{Symbol,Any}(:cancellation_check => () -> _check_powerflow_cancelled!(cancellation_token), :phase_callback => phase -> emit_phase(String(phase)), :output_dir => output_path)
  execution_start = time_ns()
  try
    open(logfile, "a") do io
      _write_resolved_q_limit_options(io, qlimit_metadata)
      with_logger(ConsoleLogger(io)) do
        redirect_stdout(io) do
          redirect_stderr(io) do
            cd(output_path) do
              raw_result = run_sparlectra(casefile = basename(case_path), path = dirname(case_path), config = config, performance_profile = api_performance_profile)
            end
          end
        end
      end
    end
  catch err
    err isa PowerFlowAborted && rethrow()
    if err isa MatpowerIO.UnsupportedMatpowerDclineError
      details = Dict{String,Any}(String(key) => value for (key, value) in err.details)
      details["solver_status"] = "aborted"
      details["service_status"] = "failed"
      details["run_status"] = "failed"
      details["last_phase"] = "reading_matpower_case"
      operation_callback("matpower_dcline_detected"; run_id = run_id, _metadata_kwargs(details)...)
      operation_callback("matpower_dcline_unsupported"; run_id = run_id, _metadata_kwargs(details)...)
      operation_callback("powerflow_aborted_unsupported_matpower_dcline"; run_id = run_id, _metadata_kwargs(details)...)
      _write_run_metadata_artifact(output_path; case_path = case_path, lifecycle = details)
      return _api_execution_failure("unsupported_matpower_dcline", err.message; run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file, phase_recorder, performance_timing, metadata = details)
    end
    message = sprint(showerror, err, catch_backtrace())
    reason = get(phase_recorder.timings[phase_recorder.active_index === nothing ? length(phase_recorder.timings) : phase_recorder.active_index], "phase", "") == "loading_julia_case" ? "loading_julia_case_failed" : "execution_error"
    return _api_execution_failure(reason, message; run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file, phase_recorder, performance_timing)
  end
  auto_profile_result = raw_result.performance_profile isa AbstractDict ? get(raw_result.performance_profile, :matpower_auto_profile_result, nothing) : nothing
  auto_profile_casefile = raw_result.performance_profile isa AbstractDict ? get(raw_result.performance_profile, :matpower_auto_profile_casefile, case_path) : case_path
  if auto_profile_result !== nothing
    config = auto_profile_result.config
    _update_effective_matpower_raw!(effective_raw, config)
    _write_yaml_file(effective_config, effective_raw)
    _write_matpower_auto_profile_artifact(output_path, auto_profile_result, config; casefile = String(auto_profile_casefile))
    operation_callback("powerflow_effective_options"; run_id = run_id, case = basename(case_path), _metadata_kwargs(qlimit_metadata)..., _metadata_kwargs(_resolved_matpower_import_runtime_options(config))...)
  end
  phases[:case_loading_network_solver] = _api_elapsed_seconds(execution_start)
  raw_result.solver_elapsed_s === nothing || (phases[:solver] = raw_result.solver_elapsed_s)
  if raw_result.performance_profile isa AbstractDict && haskey(raw_result.performance_profile, :postprocess_losses_and_flows)
    phases[:postprocessing] = Float64(raw_result.performance_profile[:postprocess_losses_and_flows])
  end

  artifact_start = time_ns()
  emit_phase("writing_diagnostics")
  _check_powerflow_cancelled!(cancellation_token)
  emit_phase("writing_artifacts")
  operation_callback("powerflow_lifecycle_status"; run_id = run_id, solver_status = "completed", artifact_status = "running", run_status = "finalizing", last_phase = "writing_artifacts")
  run_diagnostics && _write_powerflow_diagnostics(joinpath(output_path, "diagnose.log"), raw_result; mode = config.output.logfile_diagnostics)
  q_limit_artifacts = raw_result.net !== nothing ? [_write_q_limit_log_artifact(output_path, raw_result, qlimit_metadata)] : String[]
  if (run_diagnostics || detailed_result_csv) && raw_result.net !== nothing
    append!(q_limit_artifacts, _write_q_limit_detail_artifacts(output_path, raw_result.net; format = "technical"))
  end
  csv_artifacts = String[]
  csv_export_error = nothing
  csv_export_status = detailed_result_csv ? "pending" : "disabled"
  csv_export_skip_reason = nothing
  csv_format = nothing
  csv_timing_metadata = Dict{Symbol,Any}(:progress_callback => ((event; fields...) -> operation_callback(event; run_id = run_id, fields...)))
  if detailed_result_csv
    csv_format_name = detailed_result_csv_format === nothing ? (detailed_result_csv_semicolon ? "excel_de" : "technical") : String(detailed_result_csv_format)
    csv_format = try
      _resolve_detailed_csv_format(csv_format_name)
    catch err
      return _api_failure("invalid_detailed_result_csv_format", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)
    end
  end
  if !detailed_result_csv
    csv_export_skip_reason = "csv_disabled"
  elseif raw_result.net === nothing
    csv_export_status = "skipped"
    csv_export_skip_reason = "no_network_available"
  else
    emit_phase("writing_csv_artifacts")
    operation_callback("powerflow_lifecycle_status"; run_id = run_id, solver_status = "completed", artifact_status = "running", run_status = "finalizing", last_phase = "writing_csv_artifacts")
    try
      csv_artifacts = _write_detailed_result_csv(output_path, raw_result; format = csv_format.name, config, abort_checker = () -> _check_powerflow_cancelled!(cancellation_token), timing_metadata = csv_timing_metadata)
      csv_export_status = haskey(csv_timing_metadata, :partial_error) ? "partial" : (raw_result.final_converged && raw_result.solution_available ? "exported" : "exported_diagnostic")
    catch err
      csv_export_error = sprint(showerror, err)
      csv_export_status = "failed"
      operation_callback("detailed_result_csv_export_failed"; run_id = run_id, error = csv_export_error, final_converged = raw_result.final_converged, solution_available = raw_result.solution_available, has_network = raw_result.net !== nothing, csv_format = csv_format.name, status = "failed")
      open(logfile, "a") do io
        println(io, "Detailed result CSV export failed: ", csv_export_error)
      end
    end
  end
  if detailed_result_csv && csv_export_status == "skipped"
    operation_callback("detailed_result_csv_export_skipped"; run_id = run_id, reason = csv_export_skip_reason, final_converged = raw_result.final_converged, solution_available = raw_result.solution_available, has_network = raw_result.net !== nothing, csv_format = csv_format.name, status = "skipped")
    open(logfile, "a") do io
      println(io, "Detailed result CSV export skipped: ", csv_export_skip_reason)
      println(io, "  final_converged: ", raw_result.final_converged)
      println(io, "  solution_available: ", raw_result.solution_available)
      println(io, "  has_network: ", raw_result.net !== nothing)
      println(io, "  csv_format: ", csv_format.name)
    end
  end
  _check_powerflow_cancelled!(cancellation_token)
  _complete_active_phase!(phase_recorder)
  emit_phase("finalizing_success")
  _complete_active_phase!(phase_recorder)
  phases[:artifact_writing] = _api_elapsed_seconds(artifact_start)
  phases[:total] = _api_elapsed_seconds(total_start)
  open(logfile, "a") do io
    println(io)
    println(io, "API run summary")
    println(io, "===============")
    _write_api_timing_summary(io, raw_result, config, phases)
    println(io, "Case file format: ", lowercase(splitext(case_path)[2]))
    println(io, "Case file size: ", filesize(case_path), " bytes")
    println(io)
    _write_service_phase_summary(io, phase_recorder.timings)
    _write_large_case_timing_summary(io, case_path, phase_recorder.timings, raw_result)
    if !isempty(q_limit_artifacts)
      println(io, "PV Q-limit details")
      println(io, "------------------")
      println(io, "full details        : ", Q_LIMIT_LOG_ARTIFACT)
      println(io, "full initial limits : ", "q_limit_initial_limits.csv" in q_limit_artifacts ? "q_limit_initial_limits.csv" : "not generated")
      println(io, "full events         : ", "q_limit_events.csv" in q_limit_artifacts ? "q_limit_events.csv" : "not generated")
      println(io, "classic outer loop  : ", "q_limit_classic_outer_loop.csv" in q_limit_artifacts ? "q_limit_classic_outer_loop.csv" : "not generated")
      println(io)
    end
    if config.output.logfile_results === :full
      println(io, "Full run details")
      println(io, "----------------")
      println(io, "casefile: ", case_path)
      println(io, "runtime_casefile: ", basename(case_path))
      println(io, "configured_default_casefile: ", config.matpower.case)
      println(io, "config_file: ", config_path)
      println(io, "output_dir: ", output_path)
      println(io, "diagnostics_artifact: ", run_diagnostics ? "diagnose.log" : "disabled")
      println(io, "performance_artifact: ", timing_mode === :off ? "disabled" : "performance.log")
      println(io, "Q-limit handling enabled : ", !config.powerflow.qlimits.ignore_q_limits)
      println(io, "Q-limit enforcement mode : ", String(config.powerflow.qlimits.enforcement_mode))
      println(io, "q_limit_detail_artifacts: ", isempty(q_limit_artifacts) ? "none" : join(q_limit_artifacts, ", "))
      println(io, "detailed_result_csv_status: ", csv_export_status)
      csv_export_skip_reason === nothing || println(io, "detailed_result_csv_skip_reason: ", csv_export_skip_reason)
      if detailed_result_csv && csv_export_skip_reason === nothing
        println(io, "detailed_result_csv_solution_quality: ", _csv_solution_quality(raw_result))
        raw_result.final_converged && raw_result.solution_available || println(io, "detailed_result_csv_warning: values are from the last Newton iterate and are not a converged power-flow solution")
      end
      println(io, "detailed_result_csv_artifacts: ", detailed_result_csv ? (isempty(csv_artifacts) ? csv_export_status : join(csv_artifacts, ", ")) : "disabled")
      println(io, "detailed_result_csv_format: ", detailed_result_csv ? csv_format.name : "disabled")
      println(io, "detailed_result_csv_exporter: ", detailed_result_csv && raw_result.net !== nothing ? _select_detailed_csv_exporter(raw_result.net; config) : "disabled")
      println(io, "detailed_result_csv_write_mode: ", detailed_result_csv ? config.output.detailed_result_csv_write_mode : "disabled")
      println(io, "detailed_result_csv_actual_write_mode: ", detailed_result_csv && raw_result.net !== nothing && _select_detailed_csv_exporter(raw_result.net; config) === :direct ? "streaming" : (detailed_result_csv ? config.output.detailed_result_csv_write_mode : "disabled"))
      println(io, "detailed_result_csv_bus_rows: ", detailed_result_csv && raw_result.net !== nothing ? length(raw_result.net.nodeVec) : "disabled")
      println(io, "detailed_result_csv_branch_rows: ", detailed_result_csv && raw_result.net !== nothing ? length(raw_result.net.branchVec) : "disabled")
      println(io, "detailed_result_csv_delimiter: ", detailed_result_csv ? (csv_format.delimiter == ';' ? "semicolon" : "comma") : "disabled")
      println(io, "detailed_result_csv_decimal_separator: ", detailed_result_csv ? (csv_format.decimal_separator == ',' ? "comma" : "dot") : "disabled")
      println(io, "detailed_result_csv_thousands_separator: ", detailed_result_csv ? (csv_format.thousands_separator == "," ? "comma" : csv_format.thousands_separator == "." ? "dot" : "none") : "disabled")
      if haskey(csv_timing_metadata, :exporter)
        println(io, "Detailed CSV exporter       : ", csv_timing_metadata[:exporter])
        println(io, "CSV write mode              : ", csv_timing_metadata[:write_mode])
        println(io, "CSV format                  : ", csv_timing_metadata[:format])
        println(io, "Bus CSV rows                : ", csv_timing_metadata[:bus_rows])
        println(io, "Branch CSV rows             : ", csv_timing_metadata[:branch_rows])
        println(io, "Bus control label cache time: ", @sprintf("%.6f s", csv_timing_metadata[:control_label_cache_s]))
        println(io, "Bus CSV export time         : ", @sprintf("%.6f s", csv_timing_metadata[:bus_export_s]))
        println(io, "Branch CSV export time      : ", @sprintf("%.6f s", csv_timing_metadata[:branch_export_s]))
        println(io, "CSV total export time       : ", @sprintf("%.6f s", csv_timing_metadata[:total_export_s]))
        println(io, "CSV formatting mode         : ", csv_timing_metadata[:formatting_mode])
      end
      csv_export_error === nothing || println(io, "detailed_result_csv_error: ", csv_export_error)
      haskey(csv_timing_metadata, :partial_error) && println(io, "detailed_result_csv_error: ", csv_timing_metadata[:partial_error])
      println(io)
      print_effective_config(io, config)
      println(io)
      println(io, "Available status diagnostics")
      println(io, "----------------------------")
      for key in keys(raw_result.diagnostics)
        println(io, key, ": ", getproperty(raw_result.diagnostics, key))
      end
      println(io)
    end
    _write_numerical_outcome_summary(io, raw_result)
  end
  timing_mode === :off || _write_performance_log(joinpath(output_path, "performance.log"), timing_mode, phases, raw_result, phase_recorder.timings)
  _check_powerflow_cancelled!(cancellation_token)

  numerical_success = raw_result.final_converged && raw_result.solution_available
  mismatch = isfinite(raw_result.final_mismatch) ? raw_result.final_mismatch : nothing
  final_outcome = _final_outcome_payload(raw_result)
  rect_status = raw_result.net === nothing ? nothing : rectangular_pf_status(raw_result.net)
  current_iteration_metadata = rect_status === nothing ? Dict{String,Any}(
    "current_iteration_enabled" => false,
    "current_iteration_attempted" => false,
    "current_iteration_accepted" => false,
    "current_iteration_iterations" => 0,
    "current_iteration_initial_mismatch" => NaN,
    "current_iteration_final_mismatch" => NaN,
    "current_iteration_reason" => "unavailable",
    "current_iteration_artifact" => "",
  ) : Dict{String,Any}(
    "current_iteration_enabled" => getproperty(rect_status, :current_iteration_enabled),
    "current_iteration_attempted" => getproperty(rect_status, :current_iteration_attempted),
    "current_iteration_accepted" => getproperty(rect_status, :current_iteration_accepted),
    "current_iteration_iterations" => getproperty(rect_status, :current_iteration_iterations),
    "current_iteration_initial_mismatch" => getproperty(rect_status, :current_iteration_initial_mismatch),
    "current_iteration_final_mismatch" => getproperty(rect_status, :current_iteration_final_mismatch),
    "current_iteration_reason" => String(getproperty(rect_status, :current_iteration_reason)),
    "current_iteration_artifact" => getproperty(rect_status, :current_iteration_artifact),
  )
  classic_outer_loop_passes = rect_status !== nothing && hasproperty(rect_status, :matpower_outer_iterations) ? rect_status.matpower_outer_iterations : 0
  pv_to_pq_events = raw_result.net === nothing ? 0 : length(raw_result.net.qLimitLog)
  active_set_events = config.powerflow.qlimits.enforcement_mode === :active_set ? pv_to_pq_events : 0
  final_metadata = merge(Dict{String,Any}(
    "solver_status" => "completed",
    "service_status" => "completed",
    "numerical_status" => numerical_success ? "converged" : "not_converged",
    "artifact_status" => csv_export_error === nothing ? (csv_export_status == "partial" ? "partial" : "completed") : "failed",
    "run_status" => numerical_success ? "completed" : "completed_nonconverged",
    "last_phase" => "finalizing_success",
    "last_heartbeat" => Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ"),
    "final_outcome" => final_outcome,
    "detailed_result_csv_status" => csv_export_status,
    "detailed_result_csv_skip_reason" => csv_export_skip_reason,
    "detailed_result_csv_error" => csv_export_error,
    "detailed_result_csv_artifacts" => csv_artifacts,
    "detailed_result_csv_solution_quality" => detailed_result_csv && csv_export_skip_reason === nothing ? _csv_solution_quality(raw_result) : nothing,
    "detailed_result_csv_warning" => detailed_result_csv && csv_export_skip_reason === nothing && !numerical_success ? "values are from the last Newton iterate and are not a converged power-flow solution" : nothing,
    "q_limit_active_set_events" => active_set_events,
    "q_limit_pv_to_pq_events" => pv_to_pq_events,
    "q_limit_classic_outer_loop_passes" => classic_outer_loop_passes,
    "webui_request_settings" => merge(Dict{String,Any}(String(key) => value for (key, value) in config_overrides), Dict{String,Any}(
      "casefile" => casefile,
      "config_file" => config_file,
      "performance_timing" => String(performance_timing),
      "run_diagnostics" => run_diagnostics,
      "detailed_result_csv" => detailed_result_csv,
      "detailed_result_csv_format" => detailed_result_csv_format === nothing ? "technical" : String(detailed_result_csv_format),
    )),
  ), qlimit_metadata, current_iteration_metadata)
  haskey(csv_timing_metadata, :partial_error) && (final_metadata["detailed_result_csv_error"] = csv_timing_metadata[:partial_error])
  _write_run_metadata_artifact(output_path; case_path = case_path, lifecycle = final_metadata)
  message = numerical_success ? "PowerFlow run completed." : "PowerFlow run completed, but numerical solver did not converge."
  result = _api_result(
    run_id = run_id,
    status = numerical_success ? :succeeded : :not_converged,
    success = numerical_success,
    converged = raw_result.numerical_converged,
    solution_available = raw_result.solution_available,
    iterations = raw_result.iterations,
    final_mismatch = mismatch,
    reason = String(raw_result.reason),
    message = message,
    casefile = case_path,
    config_file = config_path,
    output_dir = output_path,
    logfile = logfile,
    result_file = result_file,
    service_phase_timings = phase_recorder.timings,
    metadata = final_metadata,
    raw_result = raw_result,
  )
  return _finalize_api_result(result)
end

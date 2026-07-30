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

# Run finalization helpers own run.log/performance.log details so the API
# orchestrator records phase boundaries without embedding every artifact format.

"""
    _capture_run_output(f, io; live=false)

Run `f()` with stdout/stderr/logging captured into `io` (the archived
`run.log` stream). With `live = true` (config `output.console_live`) the same
bytes are additionally mirrored to the process's real stdout while the run
executes — console and archive become two views of one stream instead of the
console being silenced. The archive content is identical in both modes.

Implementation note: `redirect_stdout` needs an OS-level stream, so the live
mode tees through a `Pipe` pumped by an async task; the pump is drained
before returning so `io` is complete when the caller closes the logfile.
"""
function _capture_run_output(f::Function, io::IO; live::Bool = false)
  if !live
    return with_logger(ConsoleLogger(io)) do
      redirect_stdout(io) do
        redirect_stderr(io) do
          f()
        end
      end
    end
  end
  real_out = stdout
  pipe = Pipe()
  Base.link_pipe!(pipe; reader_supports_async = true, writer_supports_async = true)
  pump = @async begin
    while !eof(pipe)
      data = readavailable(pipe)
      isempty(data) && continue
      write(io, data)
      write(real_out, data)
      flush(real_out)
    end
  end
  try
    return with_logger(ConsoleLogger(pipe)) do
      redirect_stdout(pipe) do
        redirect_stderr(pipe) do
          f()
        end
      end
    end
  finally
    close(pipe.in)
    wait(pump)
  end
end

function _write_service_phase_summary(io::IO, phase_timings::AbstractVector)
  println(io, "Phase timings")
  println(io, "-------------")
  # Repeated phases (newton_iteration / q_limit_processing / linear_solve fire
  # once per NR iteration) are aggregated into one line with a count — the
  # readable narrative stays a handful of lines per run. The raw per-entry
  # sequence remains available in result.json (`service_phase_timings`).
  order = String[]
  agg = Dict{String,NamedTuple{(:count, :total, :running, :status),Tuple{Int,Float64,Bool,String}}}()
  for timing in phase_timings
    phase = String(get(timing, "phase", "unknown"))
    elapsed = get(timing, "elapsed_seconds", nothing)
    prev = get(agg, phase, nothing)
    prev === nothing && push!(order, phase)
    agg[phase] = (
      count = (prev === nothing ? 0 : prev.count) + 1,
      total = (prev === nothing ? 0.0 : prev.total) + (elapsed === nothing ? 0.0 : Float64(elapsed)),
      running = (prev === nothing ? false : prev.running) || elapsed === nothing,
      status = String(get(timing, "status", "unknown")),
    )
  end
  for phase in order
    item = agg[phase]
    elapsed_text = item.running && item.count == 1 ? "running" : string(round(item.total; digits = 6), " s")
    suffix = item.count > 1 ? string(" (", item.status, ", ×", item.count, ")") : string(" (", item.status, ")")
    println(io, "  ", rpad(phase * ":", 31), elapsed_text, suffix)
  end
  println(io)
  return nothing
end

function _service_total_timing(total_start_ns::UInt64; status::AbstractString = "completed")
  elapsed = max(0.0, (time_ns() - total_start_ns) / 1.0e9)
  ended_at = _api_datetime_string(Dates.now(Dates.UTC))
  return Dict{String,Any}(
    "phase" => "total_service",
    "started_at" => nothing,
    "ended_at" => ended_at,
    "elapsed_seconds" => isfinite(elapsed) ? elapsed : 0.0,
    "status" => String(status),
    "timing_role" => "envelope",
  )
end

function _finalize_service_timings!(phase_recorder::PowerFlowPhaseTimingRecorder, total_start_ns::UInt64; status::AbstractString = "completed")
  _complete_active_phase!(phase_recorder, status)
  filter!(timing -> get(timing, "phase", "") != "total_service", phase_recorder.timings)
  push!(phase_recorder.timings, _service_total_timing(total_start_ns; status = status))
  phase_recorder.active_index = nothing
  return phase_recorder.timings
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
  if result !== nothing && result.solver_elapsed_s !== nothing
    println(io, "solver_elapsed_s: ", round(Float64(result.solver_elapsed_s); digits = 6))
  end
  for (label, phase) in (
    ("reading_matpower_case_seconds", "reading_matpower_case"),
    ("building_sparlectra_net_seconds", "building_sparlectra_net"),
    ("building_ybus_seconds", "building_ybus"),
    ("preparing_start_values_seconds", "preparing_start_values"),
    ("start_projection_seconds", "start_projection"),
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
    "matpower_apply_bus_names" => config.matpower.apply_bus_names,
    "matpower_apply_branch_names" => config.matpower.apply_branch_names,
    "matpower_apply_branch_kind" => config.matpower.apply_branch_kind,
    "matpower_import_for001_contingencies" => config.matpower.import_for001_contingencies,
    "matpower_dcline_mode" => String(config.matpower.matpower_dcline_mode),
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

# `_perf_profile_time!` accumulates fine-grained sub-timings (per phase:
# calls / summed seconds / summed allocation bytes) under
# `performance_profile[:timings]`. Islands share one profile dict, so these
# sums cover all islands of the run. Without this section they were invisible
# in performance.log (the raw-dump section suppresses the dict as too long),
# which is how 3 s of "solver_finalization" stayed unexplained in issue #288.
function _write_internal_subtimings(io::IO, result::SparlectraRunResult)
  result.performance_profile isa AbstractDict || return nothing
  timings = get(result.performance_profile, :timings, nothing)
  timings isa AbstractDict && !isempty(timings) || return nothing
  println(io)
  println(io, "Internal sub-timings (all islands summed)")
  println(io, "-----------------------------------------")
  entries = sort!(collect(timings); by = item -> -item.second.elapsed_s)
  for (phase, item) in entries
    line = string("  ", rpad(String(phase) * ":", 42), lpad(string(round(item.elapsed_s; digits = 3)), 9), " s  (calls: ", item.calls)
    println(io, line, item.bytes > 0 ? string(", alloc: ", Base.format_bytes(item.bytes), ")") : ")")
  end
  return nothing
end

# Some `performance_profile` entries (per-island diagnostics/statuses) hold
# deeply nested vectors/dicts already duplicated in dedicated per-island
# artifacts (`ac_island_<id>_solver.log`, `ac_island_solver_summary.csv`).
# Dumped raw, their default Julia `string(...)` representation is one
# multi-thousand-character line — unreadable, and heavy enough that some text
# editors struggle to open the file. Internal callback closures carry no
# diagnostic content (their printed form is just an opaque `#N` gensym), so
# they are omitted entirely rather than summarized.
const _PERFORMANCE_LOG_VERBOSE_VALUE_MAX_CHARS = 500

function _compact_performance_profile_value(key::Symbol, value)::Union{Nothing,String}
  value isa Function && return nothing
  if key === :ac_island_diagnostics || key === :ac_island_solver_statuses
    return "$(length(value)) island entry/entries (full detail in ac_island_<id>_solver.log / ac_island_solver_summary.csv)"
  end
  text = string(value)
  length(text) > _PERFORMANCE_LOG_VERBOSE_VALUE_MAX_CHARS && return "$(typeof(value)) ($(length(text)) chars, compact representation omitted; see structured result metadata)"
  return text
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
    _write_internal_subtimings(io, result)
    if mode === :full && result.performance_profile isa AbstractDict
      println(io)
      println(io, "Available internal performance profile")
      println(io, "--------------------------------------")
      for key in sort!(collect(keys(result.performance_profile)); by = string)
        summary = _compact_performance_profile_value(key, result.performance_profile[key])
        summary === nothing || println(io, key, ": ", summary)
      end
    end
  end
  return path
end

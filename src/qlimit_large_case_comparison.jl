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

const QLIMIT_LARGE_CASE_MODES = (:active_set, :classic_simultaneous, :classic_one_at_a_time)
const QLIMIT_LARGE_CASE_DEFAULTS = ("case13659pegase.m", "case_SyntheticUSA.m")
const QLIMIT_LARGE_CASE_SUMMARY_BASENAME = "qlimit_large_case_mode_comparison"
const QLIMIT_ARTIFACT_NAMES = ("q_limit.log", "q_limit_events.csv", "q_limit_initial_limits.csv", "q_limit_classic_outer_loop.csv")

"""
    compare_qlimit_large_case_modes(; cases=..., output_root=..., case_cache_dir=...)

Resolve each requested MATPOWER case through Sparlectra's case download helper,
then run the standard PowerFlow API once for each Q-limit enforcement mode.
Missing or externally unavailable optional large cases are recorded as skipped
rows instead of aborting the whole diagnostic. The function writes compact CSV
and JSON summaries and returns their paths plus all row dictionaries.
"""
function compare_qlimit_large_case_modes(;
  cases::AbstractVector{<:AbstractString} = collect(QLIMIT_LARGE_CASE_DEFAULTS),
  output_root::AbstractString = joinpath(default_webui_output_root(), "diagnostics", "qlimit_large_case_modes"),
  case_cache_dir::AbstractString = default_webui_case_cache_dir(),
  config_file::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH,
  modes::AbstractVector{Symbol} = collect(QLIMIT_LARGE_CASE_MODES),
  resolver = _resolve_qlimit_large_case,
  runner = run_sparlectra_api,
  io::IO = stdout,
)
  outroot = abspath(output_root)
  mkpath(outroot)
  rows = Vector{Dict{String,Any}}()
  for case in cases
    resolved = try
      resolver(case; outdir = case_cache_dir)
    catch err
      Dict{String,Any}("available" => false, "path" => "", "status" => "skipped", "source" => "unavailable", "reason" => sprint(showerror, err))
    end
    if !Bool(get(resolved, "available", false))
      reason = String(get(resolved, "reason", "case could not be resolved or downloaded"))
      println(io, "Skipping ", case, ": ", reason)
      push!(rows, _qlimit_skip_row(case, resolved, modes, reason))
      continue
    end
    for mode in modes
      run_dir = joinpath(outroot, _safe_summary_token(case), String(mode))
      row = _run_qlimit_mode_row(case, String(resolved["path"]), String(get(resolved, "source", "resolved")), mode, run_dir, config_file, runner)
      push!(rows, row)
    end
  end
  csv_path = joinpath(outroot, QLIMIT_LARGE_CASE_SUMMARY_BASENAME * ".csv")
  json_path = joinpath(outroot, QLIMIT_LARGE_CASE_SUMMARY_BASENAME * ".json")
  _write_qlimit_summary_csv(csv_path, rows)
  _write_qlimit_summary_json(json_path, rows)
  return (rows = rows, csv_path = csv_path, json_path = json_path, output_root = outroot)
end

function _resolve_qlimit_large_case(case::AbstractString; outdir::AbstractString)
  local_path = joinpath(outdir, case)
  existed = isfile(local_path)
  path = ensure_casefile(case; outdir = outdir, to_jl = false)
  return Dict{String,Any}("available" => true, "path" => path, "status" => existed ? "present" : "downloaded", "source" => existed ? "cache" : "download")
end

function _run_qlimit_mode_row(case::AbstractString, case_path::AbstractString, source::AbstractString, mode::Symbol, run_dir::AbstractString, config_file::AbstractString, runner)
  overrides = Dict{String,Any}(
    "matpower_import.auto_profile" => "apply",
    "matpower_import.pv_voltage_source" => "gen_vg",
    "matpower_import.compare_voltage_reference" => "imported_setpoint",
    "matpower_import.bus_shunt_model" => "admittance",
    "power_flow.max_iter" => 80,
    "power_flow.tol" => 1.0e-8,
    "power_flow.start_mode.voltage_mode" => "classic",
    "power_flow.start_mode.angle_mode" => "classic",
    "power_flow.qlimits.enabled" => true,
    "power_flow.qlimits.enforcement_mode" => String(mode),
    "output.logfile_results" => "full",
  )
  started = time_ns()
  try
    result = runner(casefile = case_path, config_file = config_file, output_dir = run_dir, config_overrides = overrides, run_diagnostics = true, performance_timing = :compact)
    elapsed = (time_ns() - started) / 1.0e9
    return _qlimit_result_row(case, case_path, source, mode, result, elapsed)
  catch err
    return _qlimit_api_error_row(case, case_path, source, mode, run_dir, sprint(showerror, err, catch_backtrace()), (time_ns() - started) / 1.0e9)
  end
end

function _qlimit_result_row(case, case_path, source, mode, result::SparlectraApiResult, elapsed_s)
  metadata = result.metadata
  artifacts = Dict(artifact.name => artifact.path for artifact in result.artifacts if artifact.name in QLIMIT_ARTIFACT_NAMES && artifact.exists)
  return Dict{String,Any}(
    "case" => case,
    "resolved_case_path" => case_path,
    "case_source" => source,
    "case_status" => "available",
    "mode" => String(mode),
    "api_status" => "completed",
    "run_status" => String(get(metadata, "run_status", String(result.status))),
    "numerical_status" => String(get(metadata, "numerical_status", result.success ? "converged" : "not_converged")),
    "reason" => something(result.reason, ""),
    "message" => something(result.message, ""),
    "iterations" => result.iterations,
    "final_mismatch" => result.final_mismatch,
    "runtime_wall_s" => elapsed_s,
    "pv2pq_event_count" => get(metadata, "q_limit_pv_to_pq_events", ""),
    "q_limit_active_set_changes" => get(metadata, "q_limit_active_set_events", ""),
    "classic_outer_loop_passes" => get(metadata, "q_limit_classic_outer_loop_passes", ""),
    "q_limit_log" => get(artifacts, "q_limit.log", ""),
    "q_limit_events_csv" => get(artifacts, "q_limit_events.csv", ""),
    "q_limit_initial_limits_csv" => get(artifacts, "q_limit_initial_limits.csv", ""),
    "q_limit_classic_outer_loop_csv" => get(artifacts, "q_limit_classic_outer_loop.csv", ""),
    "run_directory" => result.output_dir,
  )
end

function _qlimit_api_error_row(case, case_path, source, mode, run_dir, reason, elapsed_s)
  return Dict{String,Any}(
    "case" => case, "resolved_case_path" => case_path, "case_source" => source, "case_status" => "available",
    "mode" => String(mode), "api_status" => "failed", "run_status" => "failed", "numerical_status" => "not_run",
    "reason" => reason, "message" => "PowerFlow API call failed before producing a diagnostic result.",
    "iterations" => "", "final_mismatch" => "", "runtime_wall_s" => elapsed_s,
    "pv2pq_event_count" => "", "q_limit_active_set_changes" => "", "classic_outer_loop_passes" => "",
    "q_limit_log" => "", "q_limit_events_csv" => "", "q_limit_initial_limits_csv" => "", "q_limit_classic_outer_loop_csv" => "",
    "run_directory" => run_dir,
  )
end

function _qlimit_skip_row(case, resolved, modes, reason)
  return Dict{String,Any}(
    "case" => case, "resolved_case_path" => String(get(resolved, "path", "")), "case_source" => String(get(resolved, "source", "unavailable")),
    "case_status" => "skipped", "mode" => join(String.(modes), "|"), "api_status" => "skipped", "run_status" => "skipped",
    "numerical_status" => "not_run", "reason" => reason, "message" => "Optional large case was not available.",
    "iterations" => "", "final_mismatch" => "", "runtime_wall_s" => "", "pv2pq_event_count" => "",
    "q_limit_active_set_changes" => "", "classic_outer_loop_passes" => "", "q_limit_log" => "", "q_limit_events_csv" => "",
    "q_limit_initial_limits_csv" => "", "q_limit_classic_outer_loop_csv" => "", "run_directory" => "",
  )
end

const _QLIMIT_SUMMARY_COLUMNS = [
  "case", "resolved_case_path", "case_source", "case_status", "mode", "api_status", "run_status", "numerical_status",
  "reason", "message", "iterations", "final_mismatch", "runtime_wall_s", "pv2pq_event_count", "q_limit_active_set_changes",
  "classic_outer_loop_passes", "q_limit_log", "q_limit_events_csv", "q_limit_initial_limits_csv", "q_limit_classic_outer_loop_csv", "run_directory",
]

function _csv_escape(value)
  text = value === nothing ? "" : string(value)
  return occursin(r"[\",\n\r]", text) ? string('"', replace(text, "\"" => "\"\""), '"') : text
end

function _write_qlimit_summary_csv(path::AbstractString, rows)
  open(path, "w") do io
    println(io, join(_QLIMIT_SUMMARY_COLUMNS, ","))
    for row in rows
      println(io, join((_csv_escape(get(row, column, "")) for column in _QLIMIT_SUMMARY_COLUMNS), ","))
    end
  end
  return path
end

function _write_qlimit_summary_json(path::AbstractString, rows)
  open(path, "w") do io
    _write_json(io, rows)
    println(io)
  end
  return path
end

_safe_summary_token(case::AbstractString) = replace(splitext(basename(case))[1], r"[^A-Za-z0-9_.-]" => "_")

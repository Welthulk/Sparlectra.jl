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
const QLIMIT_LARGE_CASE_START_PROFILES = (:classic_start, :robust_dc_start, :autoprofile_dc_start)
const QLIMIT_LARGE_CASE_SUMMARY_BASENAME = "qlimit_large_case_mode_comparison"
const QLIMIT_ARTIFACT_NAMES = ("q_limit.log", "q_limit_events.csv", "q_limit_initial_limits.csv", "q_limit_classic_outer_loop.csv", "effective_config.yaml", "matpower_auto_profile.log")

const _QLIMIT_START_PROFILE_OVERRIDES = Dict{Symbol,Dict{String,Any}}(
  :classic_start => Dict{String,Any}(
    "power_flow.start_mode.voltage_mode" => "classic",
    "power_flow.start_mode.angle_mode" => "classic",
    "power_flow.start_mode.profile_source" => "",
  ),
  :robust_dc_start => Dict{String,Any}(
    "power_flow.start_mode.voltage_mode" => "profile_blend",
    "power_flow.start_mode.angle_mode" => "dc",
    "power_flow.start_mode.profile_source" => "matpower_reference",
  ),
  :autoprofile_dc_start => Dict{String,Any}(
    "matpower_import.auto_profile" => "apply",
    "matpower_import.pv_voltage_source" => "gen_vg",
    "matpower_import.compare_voltage_reference" => "imported_setpoint",
    "matpower_import.bus_shunt_model" => "admittance",
    "power_flow.start_mode.voltage_mode" => "profile_blend",
    "power_flow.start_mode.angle_mode" => "dc",
    "power_flow.start_mode.profile_source" => "matpower_reference",
    "power_flow.autodamp" => true,
    "power_flow.qlimits.enabled" => true,
    "power_flow.qlimits.guard.max_switches" => 10,
  ),
)

"""
    compare_qlimit_large_case_modes(; cases=..., start_profiles=..., modes=..., io=devnull, verbose=false, ...)

Resolve each requested MATPOWER case through Sparlectra's case download helper,
then run the standard PowerFlow API once for each requested `case × start_profile × qlimit_mode`
combination. Missing or externally unavailable optional large cases are recorded
as skipped rows for every requested combination instead of aborting the whole
diagnostic. The function writes row-based CSV and JSON summaries and returns
their paths plus all row dictionaries.
"""
function compare_qlimit_large_case_modes(;
  cases::AbstractVector{<:AbstractString} = collect(QLIMIT_LARGE_CASE_DEFAULTS),
  output_root::AbstractString = joinpath(default_webui_output_root(), "diagnostics", "qlimit_large_case_modes"),
  case_cache_dir::AbstractString = default_webui_case_cache_dir(),
  config_file::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH,
  modes::AbstractVector{Symbol} = collect(QLIMIT_LARGE_CASE_MODES),
  start_profiles::AbstractVector{Symbol} = collect(QLIMIT_LARGE_CASE_START_PROFILES),
  verbose_logs::Bool = false,
  logfile_results::AbstractString = verbose_logs ? "full" : "compact",
  resolver = _resolve_qlimit_large_case,
  runner = run_sparlectra_api,
  io::IO = devnull,
  verbose::Bool = false,
)
  outroot = abspath(output_root)
  mkpath(outroot)
  _validate_qlimit_modes!(modes)
  _validate_qlimit_start_profiles!(start_profiles)
  rows = Vector{Dict{String,Any}}()
  total_runs = length(cases) * length(start_profiles) * length(modes)
  run_index = 0
  _progress(io, "Q-limit large-case comparison")
  _progress(io, "output root: $(outroot)")
  _progress(io, "")
  _progress(io, "Cases:")
  for (i, case) in enumerate(cases)
    _progress(io, "  $(i). $(case)")
  end
  _progress(io, "")
  _progress(io, "Profiles:")
  for profile in start_profiles
    _progress(io, "  " * rpad(String(profile), 24) * _profile_description(profile))
  end
  _progress(io, "")
  _progress(io, "Modes: $(join(String.(modes), ", "))")
  _progress(io, "Log detail: $(logfile_results)")
  _progress(io, "")
  for (case_index, case) in enumerate(cases)
    _progress(io, "[case $(case_index)/$(length(cases))] resolving $(case) ...")
    resolved = try
      resolver(case; outdir = case_cache_dir)
    catch err
      Dict{String,Any}("available" => false, "path" => "", "status" => "error", "source" => "unavailable", "reason" => string(typeof(err), ": ", sprint(showerror, err)))
    end
    if !Bool(get(resolved, "available", false))
      reason = String(get(resolved, "reason", "case could not be resolved or downloaded"))
      _progress(io, "[case $(case_index)/$(length(cases))] skipped: $(reason)")
      for profile in start_profiles, mode in modes
        run_index += 1
        push!(rows, _qlimit_skip_row(case, resolved, profile, mode, reason))
      end
      continue
    end
    case_path = String(resolved["path"])
    _progress(io, "[case $(case_index)/$(length(cases))] $(String(get(resolved, "status", "resolved"))): $(case_path)")
    for profile in start_profiles, mode in modes
      run_index += 1
      run_dir = joinpath(outroot, _safe_summary_token(case), String(profile), String(mode))
      _progress(io, "[$(run_index)/$(total_runs)] $(case) | $(profile) | $(mode) ...")
      row = _run_qlimit_mode_row(case, case_path, String(get(resolved, "status", "resolved")), String(get(resolved, "source", "resolved")), profile, mode, run_dir, config_file, runner, logfile_results)
      push!(rows, row)
      _progress(io, _qlimit_completion_line(row); always = verbose || get(row, "run_status", "") != "completed")
    end
  end
  csv_path = joinpath(outroot, QLIMIT_LARGE_CASE_SUMMARY_BASENAME * ".csv")
  json_path = joinpath(outroot, QLIMIT_LARGE_CASE_SUMMARY_BASENAME * ".json")
  _write_qlimit_summary_csv(csv_path, rows)
  _write_qlimit_summary_json(json_path, rows)
  _print_qlimit_grouped_summary(io, rows)
  _progress(io, "")
  _progress(io, "Artifacts:")
  _progress(io, "  CSV : $(csv_path)")
  _progress(io, "  JSON: $(json_path)")
  return (rows = rows, csv_path = csv_path, json_path = json_path, output_root = outroot)
end

function _progress(io::IO, msg::AbstractString; always::Bool = true)
  always || return nothing
  println(io, msg)
  flush(io)
  return nothing
end

function _validate_qlimit_modes!(modes)
  invalid = setdiff(Symbol.(modes), QLIMIT_LARGE_CASE_MODES)
  isempty(invalid) || throw(ArgumentError("unsupported Q-limit modes: $(join(String.(invalid), ", "))"))
  return nothing
end

function _validate_qlimit_start_profiles!(profiles)
  invalid = [profile for profile in profiles if !haskey(_QLIMIT_START_PROFILE_OVERRIDES, profile)]
  isempty(invalid) || throw(ArgumentError("unsupported start profiles: $(join(String.(invalid), ", "))"))
  return nothing
end

function _resolve_qlimit_large_case(case::AbstractString; outdir::AbstractString)
  local_path = joinpath(outdir, case)
  existed = isfile(local_path)
  path = ensure_casefile(case; outdir = outdir, to_jl = false)
  return Dict{String,Any}("available" => true, "path" => path, "status" => existed ? "cached" : "downloaded", "source" => existed ? "cache" : "download")
end

function _profile_overrides(profile::Symbol)
  return copy(_QLIMIT_START_PROFILE_OVERRIDES[profile])
end

function _profile_description(profile::Symbol)
  overrides = _profile_overrides(profile)
  parts = String[]
  haskey(overrides, "matpower_import.auto_profile") && push!(parts, "auto_profile=$(overrides["matpower_import.auto_profile"])")
  push!(parts, "V=$(get(overrides, "power_flow.start_mode.voltage_mode", ""))")
  push!(parts, "angle=$(get(overrides, "power_flow.start_mode.angle_mode", ""))")
  source = get(overrides, "power_flow.start_mode.profile_source", "")
  isempty(String(source)) || push!(parts, "source=$(source)")
  return join(parts, ", ")
end

function _run_qlimit_mode_row(case::AbstractString, case_path::AbstractString, case_status::AbstractString, source::AbstractString, profile::Symbol, mode::Symbol, run_dir::AbstractString, config_file::AbstractString, runner, logfile_results::AbstractString)
  profile_overrides = _profile_overrides(profile)
  overrides = merge(Dict{String,Any}(
    "matpower_import.auto_profile" => "apply",
    "matpower_import.pv_voltage_source" => "gen_vg",
    "matpower_import.compare_voltage_reference" => "imported_setpoint",
    "matpower_import.bus_shunt_model" => "admittance",
    "power_flow.max_iter" => 80,
    "power_flow.tol" => 1.0e-8,
    "power_flow.qlimits.enabled" => true,
    "power_flow.qlimits.enforcement_mode" => String(mode),
    "output.logfile_results" => logfile_results,
  ), profile_overrides)
  started = time_ns()
  try
    result = runner(casefile = case_path, config_file = config_file, output_dir = run_dir, config_overrides = overrides, run_diagnostics = true, performance_timing = :compact)
    elapsed = (time_ns() - started) / 1.0e9
    return _qlimit_result_row(case, case_path, case_status, source, profile, mode, result, elapsed, profile_overrides)
  catch err
    return _qlimit_api_error_row(case, case_path, case_status, source, profile, mode, run_dir, string(typeof(err), ": ", sprint(showerror, err)), (time_ns() - started) / 1.0e9, profile_overrides)
  end
end

function _qlimit_base_row(case, case_path, case_status, source, profile, mode, profile_overrides)
  return Dict{String,Any}(
    "case" => case,
    "case_path" => case_path,
    "case_resolve_status" => case_status,
    "case_source" => source,
    "start_profile" => String(profile),
    "start_voltage_mode" => get(profile_overrides, "power_flow.start_mode.voltage_mode", ""),
    "start_angle_mode" => get(profile_overrides, "power_flow.start_mode.angle_mode", ""),
    "start_profile_source" => get(profile_overrides, "power_flow.start_mode.profile_source", ""),
    "qlimit_mode" => String(mode),
    "mode" => String(mode),
  )
end

function _qlimit_result_row(case, case_path, case_status, source, profile, mode, result::SparlectraApiResult, elapsed_s, profile_overrides)
  metadata = result.metadata
  artifacts = Dict(artifact.name => artifact.path for artifact in result.artifacts if artifact.name in QLIMIT_ARTIFACT_NAMES && artifact.exists)
  effective_path = get(artifacts, "effective_config.yaml", joinpath(result.output_dir, "effective_config.yaml"))
  auto_log_path = get(artifacts, "matpower_auto_profile.log", joinpath(result.output_dir, "matpower_auto_profile.log"))
  qlimit_log_path = get(artifacts, "q_limit.log", "")
  fallback = _parse_run_log_fallback(result.logfile)
  matpower = _extract_effective_matpower_import(effective_path, auto_log_path)
  row = _qlimit_base_row(case, case_path, case_status, source, profile, mode, profile_overrides)
  run_status, numerical_status, error_class, error_message = _classify_api_result(result, fallback, matpower, profile_overrides)
  merge!(row, Dict{String,Any}(
    "api_completed" => true,
    "run_status" => run_status,
    "numerical_status" => numerical_status,
    "iterations" => something(result.iterations, get(fallback, "iterations", "")),
    "final_mismatch" => something(result.final_mismatch, get(fallback, "final_mismatch", "")),
    "solver_time" => get(metadata, "solver_elapsed_s", ""),
    "wall_time" => elapsed_s,
    "pv2pq_events" => get(metadata, "q_limit_pv_to_pq_events", ""),
    "pv2pq_buses" => get(metadata, "q_limit_pv_to_pq_buses", ""),
    "qlimit_active_set_changes" => get(metadata, "q_limit_active_set_events", ""),
    "qlimit_reenable_events" => get(metadata, "q_limit_reenable_events", ""),
    "qlimit_classic_outer_iterations" => get(metadata, "q_limit_classic_outer_loop_passes", ""),
    "run_directory" => result.output_dir,
    "matpower_auto_profile_requested" => String(get(profile_overrides, "matpower_import.auto_profile", "apply")),
    "matpower_auto_profile_effective" => get(matpower, "auto_profile", ""),
    "matpower_ratio" => get(matpower, "ratio", ""),
    "matpower_shift_unit" => get(matpower, "shift_unit", ""),
    "matpower_shift_sign" => get(matpower, "shift_sign", ""),
    "matpower_bus_shunt_model" => get(matpower, "bus_shunt_model", ""),
    "matpower_pv_voltage_source" => get(matpower, "pv_voltage_source", ""),
    "matpower_compare_voltage_reference" => get(matpower, "compare_voltage_reference", ""),
    "auto_profile_log_path" => isfile(auto_log_path) ? auto_log_path : "",
    "effective_config_path" => isfile(effective_path) ? effective_path : "",
    "q_limit_log_path" => qlimit_log_path,
    "q_limit_log" => qlimit_log_path,
    "q_limit_events_csv" => get(artifacts, "q_limit_events.csv", ""),
    "q_limit_initial_limits_csv" => get(artifacts, "q_limit_initial_limits.csv", ""),
    "q_limit_classic_outer_loop_csv" => get(artifacts, "q_limit_classic_outer_loop.csv", ""),
    "error_class" => error_class,
    "error_message" => error_message,
  ))
  return row
end

function _classify_api_result(result::SparlectraApiResult, fallback, matpower, profile_overrides)
  requested = String(get(profile_overrides, "matpower_import.auto_profile", "apply"))
  if requested == "apply" && !_auto_profile_confirmed(matpower)
    return ("completed", "invalid_diagnostic", "invalid_diagnostic_configuration", "auto_profile_apply_requested_but_effective_profile_not_confirmed")
  end
  if !result.success && !result.converged && result.iterations === nothing
    message = String(result.message)
    isempty(message) && (message = String(get(fallback, "final_status", "")))
    return ("api_failed_before_solve", "error", "api_failed_before_solve", message)
  end
  status = result.success || result.converged ? "api_completed_and_converged" : "api_completed_but_not_converged"
  numerical = result.converged ? "converged" : "not_converged"
  message = result.success || result.converged ? "" : String(result.message)
  isempty(message) && (message = String(get(fallback, "final_status", "")))
  return (status, numerical, result.success || result.converged ? "" : "api_completed_but_not_converged", message)
end

_auto_profile_confirmed(matpower) = get(matpower, "auto_profile", "") == "apply" && !isempty(string(get(matpower, "shift_unit", ""))) && !isempty(string(get(matpower, "shift_sign", ""))) && !isempty(string(get(matpower, "ratio", "")))

function _qlimit_api_error_row(case, case_path, case_status, source, profile, mode, run_dir, reason, elapsed_s, profile_overrides)
  row = _qlimit_base_row(case, case_path, case_status, source, profile, mode, profile_overrides)
  merge!(row, _qlimit_empty_result_fields("api_exception", "error", false, run_dir, reason; wall_time = elapsed_s, error_class = "api_exception"))
  return row
end

function _qlimit_skip_row(case, resolved, profile, mode, reason)
  profile_overrides = _profile_overrides(profile)
  row = _qlimit_base_row(case, String(get(resolved, "path", "")), String(get(resolved, "status", "skipped")), String(get(resolved, "source", "unavailable")), profile, mode, profile_overrides)
  merge!(row, _qlimit_empty_result_fields("case_resolve_failed", "skipped", false, "", reason; error_class = "case_resolve_failed"))
  return row
end

function _qlimit_empty_result_fields(run_status, numerical_status, api_completed, run_dir, error_message; wall_time = "", error_class = "")
  return Dict{String,Any}(
    "api_completed" => api_completed,
    "run_status" => run_status,
    "numerical_status" => numerical_status,
    "iterations" => "",
    "final_mismatch" => "",
    "solver_time" => "",
    "wall_time" => wall_time,
    "pv2pq_events" => "",
    "pv2pq_buses" => "",
    "qlimit_active_set_changes" => "",
    "qlimit_reenable_events" => "",
    "qlimit_classic_outer_iterations" => "",
    "run_directory" => run_dir,
    "matpower_auto_profile_requested" => "apply",
    "matpower_auto_profile_effective" => "",
    "matpower_ratio" => "",
    "matpower_shift_unit" => "",
    "matpower_shift_sign" => "",
    "matpower_bus_shunt_model" => "",
    "matpower_pv_voltage_source" => "",
    "matpower_compare_voltage_reference" => "",
    "auto_profile_log_path" => "",
    "effective_config_path" => "",
    "q_limit_log_path" => "",
    "q_limit_log" => "",
    "q_limit_events_csv" => "",
    "q_limit_initial_limits_csv" => "",
    "q_limit_classic_outer_loop_csv" => "",
    "error_class" => error_class,
    "error_message" => error_message,
  )
end

function _qlimit_completion_line(row)
  reason = isempty(String(get(row, "error_message", ""))) ? "" : ", reason=$(row["error_message"])"
  return "       status: $(row["numerical_status"]), iter=$(row["iterations"]), mismatch=$(row["final_mismatch"]), profile=$(_matpower_profile_label(row)), time=$(round(Float64(get(row, "wall_time", 0.0)); digits = 2))s$(reason)\n       dir: $(row["run_directory"])"
end

const _QLIMIT_SUMMARY_COLUMNS = [
  "case", "case_path", "case_resolve_status", "case_source", "start_profile", "start_voltage_mode", "start_angle_mode", "start_profile_source",
  "qlimit_mode", "mode", "api_completed", "run_status", "numerical_status", "iterations", "final_mismatch", "solver_time", "wall_time", "matpower_auto_profile_requested", "matpower_auto_profile_effective", "matpower_ratio", "matpower_shift_unit", "matpower_shift_sign", "matpower_bus_shunt_model", "matpower_pv_voltage_source", "matpower_compare_voltage_reference", "auto_profile_log_path", "effective_config_path", "pv2pq_events",
  "pv2pq_buses", "qlimit_active_set_changes", "qlimit_reenable_events", "qlimit_classic_outer_iterations", "run_directory", "q_limit_log",
  "q_limit_log_path", "q_limit_events_csv", "q_limit_initial_limits_csv", "q_limit_classic_outer_loop_csv", "error_class", "error_message",
]

function _csv_escape(value)
  text = value === nothing ? "" : string(value)
  return occursin(r"[\",\n\r]", text) ? string('"', replace(text, "\"" => "\"\""), '"') : text
end

function _extract_effective_matpower_import(effective_config_path::AbstractString, auto_profile_log_path::AbstractString)
  values = Dict{String,Any}(
    "auto_profile" => "",
    "ratio" => "",
    "shift_unit" => "",
    "shift_sign" => "",
    "bus_shunt_model" => "",
    "pv_voltage_source" => "",
    "compare_voltage_reference" => "",
  )
  if isfile(effective_config_path)
    try
      cfg = load_sparlectra_config(effective_config_path; reload = true)
      values["auto_profile"] = String(cfg.matpower.auto_profile)
      values["ratio"] = String(cfg.matpower.ratio)
      values["shift_unit"] = String(cfg.matpower.shift_unit)
      values["shift_sign"] = cfg.matpower.shift_sign
      values["bus_shunt_model"] = String(cfg.matpower.bus_shunt_model)
      values["pv_voltage_source"] = String(cfg.matpower.pv_voltage_source)
      values["compare_voltage_reference"] = String(cfg.matpower.compare_voltage_reference)
    catch
      # Keep row generation robust; the missing fields are classified by the caller.
    end
  end
  if isfile(auto_profile_log_path)
    text = read(auto_profile_log_path, String)
    for (field, pattern) in (
      "ratio" => r"matpower_import\.ratio\s*[:=]\s*([A-Za-z0-9_.+-]+)",
      "shift_unit" => r"matpower_import\.shift_unit\s*[:=]\s*([A-Za-z0-9_.+-]+)",
      "shift_sign" => r"matpower_import\.shift_sign\s*[:=]\s*([A-Za-z0-9_.+-]+)",
    )
      if isempty(string(values[field]))
        m = match(pattern, text)
        m === nothing || (values[field] = m.captures[1])
      end
    end
    occursin("Auto-profile recommendation:", text) && isempty(string(values["auto_profile"])) && (values["auto_profile"] = "apply")
  end
  return values
end

function _parse_run_log_fallback(path::AbstractString)
  data = Dict{String,Any}()
  isfile(path) || return data
  text = read(path, String)
  for pattern in (r"final[_ ]status\s*[:=]\s*([^\n\r]+)"i, r"status\s*[:=]\s*([^\n\r]+)"i)
    m = match(pattern, text)
    m === nothing || (data["final_status"] = strip(m.captures[1]); break)
  end
  m = match(r"iterations?\s*[:=]\s*(\d+)"i, text)
  m === nothing || (data["iterations"] = parse(Int, m.captures[1]))
  m = match(r"(?:final[_ ]?)?mismatch\s*[:=]\s*([0-9.eE+-]+)"i, text)
  m === nothing || (data["final_mismatch"] = parse(Float64, m.captures[1]))
  return data
end

function _matpower_profile_label(row)
  unit = string(get(row, "matpower_shift_unit", ""))
  sign = string(get(row, "matpower_shift_sign", ""))
  ratio = string(get(row, "matpower_ratio", ""))
  isempty(unit) && isempty(sign) && isempty(ratio) && return "unconfirmed"
  return "$(unit)/$(sign)" * (isempty(ratio) ? "" : ", ratio=$(ratio)")
end

function _print_qlimit_grouped_summary(io::IO, rows)
  _progress(io, "")
  _progress(io, "Summary by case")
  for case in unique(String(row["case"]) for row in rows)
    _progress(io, case)
    case_rows = filter(row -> row["case"] == case, rows)
    for profile in unique(String(row["start_profile"]) for row in case_rows)
      _progress(io, "  $(profile)")
      profile_rows = filter(row -> row["start_profile"] == profile, case_rows)
      for row in profile_rows
        status = String(row["numerical_status"]) == "converged" ? "CONVERGED" : uppercase(replace(String(row["numerical_status"]), "_" => " "))
        reason = isempty(String(get(row, "error_message", ""))) ? "" : "   reason=$(row["error_message"])"
        _progress(io, "    " * rpad(String(row["mode"]), 24) * rpad(status, 23) * "iter=$(row["iterations"])   mismatch=$(row["final_mismatch"])   profile=$(_matpower_profile_label(row))$(reason)")
      end
    end
    _progress(io, "")
  end
  return nothing
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

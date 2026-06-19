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

function _effective_config_with_runtime_case(effective_raw, case_path::AbstractString, config::SparlectraConfig)
  raw = deepcopy(effective_raw)
  runtime = get!(raw, "runtime", Dict{String,Any}())
  runtime isa AbstractDict || (runtime = raw["runtime"] = Dict{String,Any}())
  # The effective config may contain default form choices, but diagnostics need
  # the resolved runtime casefile so support analysis can identify what was run.
  runtime["casefile"] = String(case_path)
  runtime["case_name"] = splitext(basename(case_path))[1]
  runtime["case_source"] = "webui_mpower_data"
  runtime["configured_default_casefile"] = config.matpower.case
  return raw
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

include("run_failures.jl")

include("run_diagnostic_artifacts.jl")
include("run_finalization.jl")

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

include("run_csv_exports.jl")

function _csv_solution_quality(raw_result::SparlectraRunResult)::String
  return raw_result.final_converged && raw_result.solution_available ? "converged" : "not_converged_last_iterate"
end

include("run_lifecycle_metadata.jl")

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

include("run_matpower_artifacts.jl")

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
  _write_yaml_file(effective_config, _effective_config_with_runtime_case(effective_raw, case_path, config))
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
    _write_yaml_file(effective_config, _effective_config_with_runtime_case(effective_raw, case_path, config))
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
      # CSV export is delegated so API orchestration only chooses when exports
      # happen; the helper preserves filenames, schemas, progress events, and
      # partial-export behavior expected by the Web UI.
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
  # Lifecycle metadata is a service/API contract: Last Errors, run listings, and
  # clients consume these status fields even when the numerical result is
  # non-converged but artifact finalization succeeded.
  final_metadata = _build_success_lifecycle_metadata(
    raw_result,
    config;
    numerical_success,
    final_outcome,
    csv_export_status,
    csv_export_skip_reason,
    csv_export_error,
    csv_artifacts,
    detailed_result_csv,
    config_overrides,
    casefile,
    config_file,
    performance_timing,
    run_diagnostics,
    detailed_result_csv_format,
    qlimit_metadata,
    csv_timing_metadata,
  )
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

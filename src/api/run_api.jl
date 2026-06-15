using UUIDs: uuid4

const _SPARLECTRA_API_SCHEMA_VERSION = "1.0"
const WEBUI_PERFORMANCE_TIMING_VALUES = (:off, :compact, :full)

function _api_elapsed_seconds(start_ns::UInt64)::Float64
  return (time_ns() - start_ns) / 1.0e9
end

function _api_timing_mode(value)::Symbol
  mode = value isa Symbol ? value : Symbol(lowercase(strip(String(value))))
  mode in WEBUI_PERFORMANCE_TIMING_VALUES || throw(ArgumentError("Performance timing mode must be one of $(join(WEBUI_PERFORMANCE_TIMING_VALUES, ", ")); got $(repr(value))."))
  return mode
end

function _write_api_timing_summary(io::IO, result::SparlectraRunResult, config::SparlectraConfig)
  benchmark_median = result.performance_profile isa AbstractDict ? get(result.performance_profile, :benchmark_median_s, nothing) : nothing
  representative_time = result.performance_profile isa AbstractDict ? get(result.performance_profile, :representative_elapsed_s, result.elapsed_s) : result.elapsed_s
  println(io, "Timing")
  println(io, "------")
  println(io, "solver_time:          ", result.solver_elapsed_s === nothing ? "n/a" : "$(round(result.solver_elapsed_s; digits = 6)) s")
  println(io, "representative_time:  ", round(Float64(representative_time); digits = 6), " s")
  if config.benchmark.enabled
    println(io, "benchmark_median:     ", benchmark_median === nothing ? "n/a" : "$(round(Float64(benchmark_median); digits = 6)) s")
    println(io, "benchmark_samples:    ", config.benchmark.samples)
  end
  println(io, "iterations:           ", result.iterations)
  println(io, "final_mismatch:       ", isfinite(result.final_mismatch) ? result.final_mismatch : "n/a")
  println(io, "final_status:         ", result.outcome)
  println(io, "final_outcome:        ", result.final_converged ? "converged" : result.reason)
  println(io)
  return nothing
end

function _write_powerflow_diagnostics(path::AbstractString, result::SparlectraRunResult; diagnostic_fn = nothing)
  open(path, "w") do io
    println(io, "Sparlectra PowerFlow diagnostics")
    println(io, "================================")
    try
      if diagnostic_fn === nothing
        println(io, "outcome: ", result.outcome)
        println(io, "reason: ", result.reason_text)
        println(io)
        printQLimitLog(result.net; io)
        println(io)
        printPVQLimitsTable(result.net; io)
        println(io)
        printFinalLimitValidation(result.net; io)
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

function _write_performance_log(path::AbstractString, mode::Symbol, phases::AbstractDict, result::SparlectraRunResult)
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
  first_group = mod(length(digits), 3)
  first_group == 0 && (first_group = 3)
  groups = String[digits[1:first_group]]
  for start in (first_group + 1):3:length(digits)
    push!(groups, digits[start:(start + 2)])
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
  if format.name != "technical" && occursin(r"[eE]", technical)
    exponent_marker = findfirst(character -> character in ('e', 'E'), technical)
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
      technical = sign * digits[1:decimal_position] * "." * digits[(decimal_position + 1):end]
    end
  end
  if format.name != "technical" && occursin('.', technical)
    technical = replace(technical, r"0+$" => "")
    technical = replace(technical, r"\.$" => "")
    isempty(technical) && (technical = "0")
    technical == "-0" && (technical = "0")
  end
  mantissa_exponent = split(technical, r"(?=[eE])"; limit = 2)
  mantissa = mantissa_exponent[1]
  exponent = format.name == "technical" && length(mantissa_exponent) == 2 ? mantissa_exponent[2] : ""
  parts = split(mantissa, '.'; limit = 2)
  integer_part = _group_csv_integer(parts[1], format.thousands_separator)
  fractional_part = length(parts) == 2 ? format.decimal_separator * parts[2] : ""
  return integer_part * fractional_part * exponent
end

function _format_csv_value(value, format)::String
  value === missing && return ""
  value === nothing && return ""
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

function _write_namedtuple_csv(path::AbstractString, rows::AbstractVector, columns; delimiter::Char = ',', format = nothing)
  delimiter in (',', ';') || throw(ArgumentError("CSV delimiter must be ',' or ';'."))
  resolved_format = format === nothing ?
                    (name = "custom", delimiter = delimiter, decimal_separator = '.', thousands_separator = "") :
                    _resolve_detailed_csv_format(format)
  resolved_format.delimiter == delimiter || throw(ArgumentError("CSV delimiter does not match detailed CSV format $(resolved_format.name)."))
  open(path, "w") do io
    println(io, join(String.(columns), delimiter))
    for row in rows
      println(io, join((_csv_field(getproperty(row, column), delimiter, resolved_format) for column in columns), delimiter))
    end
  end
  return path
end

function _complex_voltage_rows(node_rows::AbstractVector, format)
  return [begin
    angle = deg2rad(Float64(row.va_deg))
    v_re = Float64(row.vm_pu) * cos(angle)
    v_im = Float64(row.vm_pu) * sin(angle)
    merge(
      row,
      (
        v_re = round(v_re; sigdigits = 10),
        v_im = round(v_im; sigdigits = 10),
        v_complex = string(
          _format_csv_number(v_re, format),
          signbit(v_im) ? " - j" : " + j",
          _format_csv_number(abs(v_im), format),
        ),
      ),
    )
  end for row in node_rows]
end

function _write_detailed_result_csv(output_path::AbstractString, result::SparlectraRunResult; format = "technical")::Vector{String}
  resolved_format = _resolve_detailed_csv_format(format)
  result.net === nothing && throw(ArgumentError("PowerFlow result does not contain a network for detailed CSV export."))
  report = buildACPFlowReport(
    result.net;
    ct = result.elapsed_s,
    ite = result.iterations,
    converged = result.final_converged,
  )
  bus_columns = (
    :bus, :bus_name, :type, :vm_pu, :va_deg, :vn_kV, :v_re, :v_im,
    :v_complex, :v_kV, :p_gen_MW, :q_gen_MVar, :p_load_MW, :q_load_MVar,
    :q_limit_hit, :control,
  )
  branch_columns = (
    :branch, :branch_index, :from_bus, :to_bus, :status, :p_from_MW,
    :q_from_MVar, :p_to_MW, :q_to_MVar, :p_loss_MW, :q_loss_MVar,
    :rated_MVA, :overloaded,
  )
  artifacts = ["bus_voltages_complex.csv", "branch_flows.csv"]
  _write_namedtuple_csv(joinpath(output_path, artifacts[1]), _complex_voltage_rows(report.nodes, resolved_format), bus_columns; delimiter = resolved_format.delimiter, format = resolved_format.name)
  _write_namedtuple_csv(joinpath(output_path, artifacts[2]), report.branches, branch_columns; delimiter = resolved_format.delimiter, format = resolved_format.name)
  return artifacts
end

function _api_result(; run_id::String = string(uuid4()), schema_version::String = _SPARLECTRA_API_SCHEMA_VERSION, status::Symbol, success::Bool, converged = nothing, solution_available::Bool = false, iterations = nothing, final_mismatch = nothing, reason = nothing, message = nothing, casefile = nothing, config_file = nothing, output_dir::String, logfile = nothing, result_file = nothing, artifacts = SparlectraApiArtifact[], raw_result = nothing)
  return SparlectraApiResult(run_id, schema_version, status, success, converged, solution_available, iterations, final_mismatch, reason, message, casefile, config_file, output_dir, logfile, result_file, artifacts, raw_result)
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

function _api_failure(reason::String, message::String; run_id::String = string(uuid4()), casefile, config_file, output_dir::String, logfile::String, result_file::String)::SparlectraApiResult
  open(logfile, "a") do io
    println(io, "Sparlectra API failure: ", reason)
    println(io, message)
  end
  result = _api_result(run_id = run_id, status = :failed, success = false, reason = reason, message = message, casefile = casefile, config_file = config_file, output_dir = output_dir, logfile = logfile, result_file = result_file)
  return _finalize_api_result(result)
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
)::SparlectraApiResult
  total_start = time_ns()
  phases = Dict{Symbol,Float64}(Symbol(key) => Float64(value) for (key, value) in phase_timings)
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
  phase_callback("preparing_configuration")
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
  _check_powerflow_cancelled!(cancellation_token)
  phases[:api_config_build] = _api_elapsed_seconds(config_start)
  isfile(case_path) || return _api_failure("casefile_not_found", "MATPOWER case file not found: $(case_path)"; run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)

  raw_result = nothing
  phase_callback("loading_case")
  api_performance_profile = Dict{Symbol,Any}(
    :cancellation_check => () -> _check_powerflow_cancelled!(cancellation_token),
    :phase_callback => phase -> phase_callback(String(phase)),
  )
  execution_start = time_ns()
  try
    open(logfile, "a") do io
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
    message = sprint(showerror, err, catch_backtrace())
    return _api_failure("execution_error", message; run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)
  end
  phases[:case_loading_network_solver] = _api_elapsed_seconds(execution_start)
  raw_result.solver_elapsed_s === nothing || (phases[:solver] = raw_result.solver_elapsed_s)
  if raw_result.performance_profile isa AbstractDict && haskey(raw_result.performance_profile, :postprocess_losses_and_flows)
    phases[:postprocessing] = Float64(raw_result.performance_profile[:postprocess_losses_and_flows])
  end

  artifact_start = time_ns()
  phase_callback("diagnostics")
  _check_powerflow_cancelled!(cancellation_token)
  phase_callback("artifact_writing")
  run_diagnostics && _write_powerflow_diagnostics(joinpath(output_path, "diagnose.log"), raw_result)
  csv_artifacts = String[]
  csv_export_error = nothing
  csv_format = nothing
  if detailed_result_csv
    csv_format_name = detailed_result_csv_format === nothing ?
                      (detailed_result_csv_semicolon ? "excel_de" : "technical") :
                      String(detailed_result_csv_format)
    csv_format = try
      _resolve_detailed_csv_format(csv_format_name)
    catch err
      return _api_failure("invalid_detailed_result_csv_format", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)
    end
  end
  if detailed_result_csv && raw_result.final_converged && raw_result.solution_available
    try
      csv_artifacts = _write_detailed_result_csv(output_path, raw_result; format = csv_format.name)
    catch err
      csv_export_error = sprint(showerror, err)
      open(logfile, "a") do io
        println(io, "Detailed result CSV export failed: ", csv_export_error)
      end
    end
  end
  _check_powerflow_cancelled!(cancellation_token)
  open(logfile, "a") do io
    println(io)
    println(io, "API run summary")
    println(io, "===============")
    _write_api_timing_summary(io, raw_result, config)
    if config.output.logfile_results === :full
      println(io, "Full run details")
      println(io, "----------------")
      println(io, "casefile: ", case_path)
      println(io, "config_file: ", config_path)
      println(io, "output_dir: ", output_path)
      println(io, "diagnostics_artifact: ", run_diagnostics ? "diagnose.log" : "disabled")
      println(io, "performance_artifact: ", timing_mode === :off ? "disabled" : "performance.log")
      println(io, "detailed_result_csv_artifacts: ", detailed_result_csv ? (isempty(csv_artifacts) ? "failed" : join(csv_artifacts, ", ")) : "disabled")
      println(io, "detailed_result_csv_format: ", detailed_result_csv ? csv_format.name : "disabled")
      println(io, "detailed_result_csv_delimiter: ", detailed_result_csv ? (csv_format.delimiter == ';' ? "semicolon" : "comma") : "disabled")
      println(io, "detailed_result_csv_decimal_separator: ", detailed_result_csv ? (csv_format.decimal_separator == ',' ? "comma" : "dot") : "disabled")
      println(io, "detailed_result_csv_thousands_separator: ", detailed_result_csv ? (csv_format.thousands_separator == "," ? "comma" : csv_format.thousands_separator == "." ? "dot" : "none") : "disabled")
      csv_export_error === nothing || println(io, "detailed_result_csv_error: ", csv_export_error)
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
  end
  phases[:artifact_writing] = _api_elapsed_seconds(artifact_start)
  phases[:total] = _api_elapsed_seconds(total_start)
  timing_mode === :off || _write_performance_log(joinpath(output_path, "performance.log"), timing_mode, phases, raw_result)
  _check_powerflow_cancelled!(cancellation_token)
  phase_callback("finalizing_success")

  success = raw_result.final_converged && raw_result.solution_available
  mismatch = isfinite(raw_result.final_mismatch) ? raw_result.final_mismatch : nothing
  result = _api_result(
    run_id = run_id,
    status = success ? :succeeded : :failed,
    success = success,
    converged = raw_result.numerical_converged,
    solution_available = raw_result.solution_available,
    iterations = raw_result.iterations,
    final_mismatch = mismatch,
    reason = success ? nothing : String(raw_result.reason),
    message = success ? nothing : raw_result.reason_text,
    casefile = case_path,
    config_file = config_path,
    output_dir = output_path,
    logfile = logfile,
    result_file = result_file,
    raw_result = raw_result,
  )
  return _finalize_api_result(result)
end

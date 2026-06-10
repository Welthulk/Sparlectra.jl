using UUIDs: uuid4

const _SPARLECTRA_API_SCHEMA_VERSION = "1.0"

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

function _api_failure(reason::String, message::String; casefile, config_file, output_dir::String, logfile::String, result_file::String)::SparlectraApiResult
  open(logfile, "a") do io
    println(io, "Sparlectra API failure: ", reason)
    println(io, message)
  end
  result = _api_result(status = :failed, success = false, reason = reason, message = message, casefile = casefile, config_file = config_file, output_dir = output_dir, logfile = logfile, result_file = result_file)
  return _finalize_api_result(result)
end

"""
    run_sparlectra_api(; casefile, config_file, output_dir, config_overrides=Dict()) -> SparlectraApiResult

Run one MATPOWER power-flow case through a stable, non-interactive API contract.
The function validates GUI overrides, writes `effective_config.yaml`, delegates
the numerical work to [`run_sparlectra`](@ref), captures textual output in
`run.log`, writes `result.json`, discovers all generated files, and returns
structured status and artifact metadata. The input configuration template is
never modified.
"""
function run_sparlectra_api(;
  casefile::AbstractString,
  config_file::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH,
  output_dir::AbstractString,
  config_overrides::AbstractDict = Dict{String,Any}(),
)::SparlectraApiResult
  output_path = abspath(output_dir)
  mkpath(output_path)
  case_path = abspath(casefile)
  config_path = abspath(config_file)
  logfile = joinpath(output_path, "run.log")
  result_file = joinpath(output_path, "result.json")
  touch(logfile)

  isfile(config_path) || return _api_failure("config_file_not_found", "Configuration file not found: $(config_path)"; casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)

  nested_overrides = try
    validate_gui_config_overrides(config_overrides)
  catch err
    return _api_failure("invalid_config_override", sprint(showerror, err); casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)
  end

  config = nothing
  effective_raw = nothing
  try
    config, effective_raw = _load_api_config(config_path, nested_overrides)
  catch err
    return _api_failure("invalid_configuration", sprint(showerror, err); casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)
  end

  effective_config = joinpath(output_path, "effective_config.yaml")
  _write_yaml_file(effective_config, effective_raw)
  isfile(case_path) || return _api_failure("casefile_not_found", "MATPOWER case file not found: $(case_path)"; casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)

  raw_result = nothing
  try
    open(logfile, "a") do io
      with_logger(ConsoleLogger(io)) do
        redirect_stdout(io) do
          redirect_stderr(io) do
            cd(output_path) do
              raw_result = run_sparlectra(casefile = basename(case_path), path = dirname(case_path), config = config)
            end
          end
        end
      end
    end
  catch err
    message = sprint(showerror, err, catch_backtrace())
    return _api_failure("execution_error", message; casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)
  end

  success = raw_result.final_converged && raw_result.solution_available
  mismatch = isfinite(raw_result.final_mismatch) ? raw_result.final_mismatch : nothing
  result = _api_result(
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

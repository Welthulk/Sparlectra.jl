const _POWERFLOW_SERVICE_RUNS = Dict{String,SparlectraApiResult}()
const _POWERFLOW_SERVICE_LOCK = ReentrantLock()


function _registered_powerflow_run(run_id::AbstractString)
  id = String(run_id)
  return lock(_POWERFLOW_SERVICE_LOCK) do
    get(_POWERFLOW_SERVICE_RUNS, id, nothing)
  end
end

function _register_powerflow_run!(result::SparlectraApiResult)
  lock(_POWERFLOW_SERVICE_LOCK) do
    _POWERFLOW_SERVICE_RUNS[result.run_id] = result
  end
  return result
end

include("service_json.jl")
include("artifact_registry.jl")
include("run_index.jl")
include("webui_jobs.jl")

function _prefer_julia_casefile(casefile::AbstractString; emit_julia_case_fn = FetchMatpowerCase.emit_julia_case)::String
  case_path = abspath(casefile)
  extension = lowercase(splitext(case_path)[2])
  extension == ".jl" && return case_path
  extension == ".m" || throw(ArgumentError("Unsupported casefile extension: $(casefile) (expected .m or .jl)"))

  julia_path = first(splitext(case_path)) * ".jl"
  isfile(julia_path) && return julia_path
  try
    emitted_path = abspath(emit_julia_case_fn(case_path, dirname(case_path)))
    if isfile(emitted_path)
      return emitted_path
    end
    @warn "MATPOWER Julia case generation did not produce a file; using the .m case" casefile = case_path emitted_path
  catch err
    @warn "MATPOWER Julia case generation failed; using the .m case" casefile = case_path exception = (err, catch_backtrace())
  end
  return case_path
end

function _resolve_powerflow_casefile(
  casefile::AbstractString,
  case_directory::AbstractString;
  ensure_casefile_fn = FetchMatpowerCase.ensure_casefile,
  emit_julia_case_fn = FetchMatpowerCase.emit_julia_case,
)::String
  requested = strip(casefile)
  isempty(requested) && throw(ArgumentError("PowerFlow casefile must not be empty."))
  occursin(r"^[A-Za-z][A-Za-z0-9+.-]*://", requested) && throw(ArgumentError("MATPOWER case URLs are not accepted."))

  extension = lowercase(splitext(requested)[2])
  extension in (".m", ".jl") || throw(ArgumentError("Unsupported casefile extension: $(requested) (expected .m or .jl)"))
  if isfile(requested)
    return _prefer_julia_casefile(requested; emit_julia_case_fn)
  end
  occursin(r"[\\/]", requested) && throw(ArgumentError("Case file not found: $(requested)"))

  trusted_directory = abspath(case_directory)
  mkpath(trusted_directory)
  resolved = try
    ensure_casefile_fn(requested; outdir = trusted_directory, to_jl = true)
  catch err
    fallback_name = extension == ".m" ? requested : first(splitext(requested)) * ".m"
    fallback_path = joinpath(trusted_directory, fallback_name)
    if isfile(fallback_path)
      @warn "MATPOWER case resolution could not generate the Julia case; using the downloaded .m case" casefile = fallback_path exception = (err, catch_backtrace())
      return abspath(fallback_path)
    end
    rethrow()
  end
  isfile(resolved) || throw(ArgumentError("Resolved MATPOWER case file not found: $(resolved)"))
  return _prefer_julia_casefile(resolved; emit_julia_case_fn)
end

"""
    start_powerflow_run(request::AbstractDict; case_directory=nothing) -> Dict{String,Any}

Start one local PowerFlow service run above [`run_sparlectra_api`](@ref). The
request must provide `casefile`, `config_file`, and `output_root`; optional
`config_overrides` are forwarded to the programmatic API. A unique run ID is
chosen before execution, and all generated files are written beneath
`output_root/run_id`. Completed API runs, including failed runs with a
`result.json`, are registered in memory and in the persistent run index. Public
service failures are returned as structured dictionaries.

When `case_directory` is provided by a trusted caller, bare `.m` or `.jl` case
names are resolved there through [`ensure_casefile`](@ref). Generated Julia
cases are preferred for execution, while existing path behavior is retained.
"""
function start_powerflow_run(request::AbstractDict; case_directory::Union{Nothing,AbstractString} = nothing, case_resolver = _resolve_powerflow_casefile)::Dict{String,Any}
  request_start = time_ns()
  phases = Dict{Symbol,Float64}()
  cancellation_token = _service_request_value(request, "cancellation_token", nothing)
  phase_callback = _service_request_value(request, "phase_callback", phase -> nothing)
  set_phase = phase -> (phase_callback(String(phase)); _check_powerflow_cancelled!(cancellation_token))
  set_phase("resolving_case")
  _check_powerflow_cancelled!(cancellation_token)
  casefile = _service_request_value(request, "casefile")
  casefile isa AbstractString && !isempty(strip(casefile)) || return _service_failure("missing_casefile", "PowerFlow service request requires a nonempty casefile.")
  if case_directory === nothing
    _check_powerflow_cancelled!(cancellation_token)
    isfile(casefile) || return _service_failure("missing_casefile", "MATPOWER case file not found: $(abspath(casefile))")
    casefile = abspath(casefile)
  else
    resolution_start = time_ns()
    _check_powerflow_cancelled!(cancellation_token)
    set_phase("checking_case_cache")
    casefile = try
      set_phase("preparing_case_file")
      case_resolver(casefile, case_directory)
    catch err
      return _service_failure("invalid_casefile", sprint(showerror, err))
    end
    phases[:case_resolution] = _api_elapsed_seconds(resolution_start)
    _check_powerflow_cancelled!(cancellation_token)
  end
  set_phase("preparing_configuration")

  config_file = _service_request_value(request, "config_file")
  config_file isa AbstractString && !isempty(strip(config_file)) || return _service_failure("missing_config_file", "PowerFlow service request requires a nonempty config_file.")
  isfile(config_file) || return _service_failure("missing_config_file", "Configuration file not found: $(abspath(config_file))")

  output_root = _service_request_value(request, "output_root")
  output_root isa AbstractString && !isempty(strip(output_root)) || return _service_failure("invalid_request", "PowerFlow service request requires a nonempty output_root.")

  config_overrides = _service_request_value(request, "config_overrides", Dict{String,Any}())
  config_overrides isa AbstractDict || return _service_failure("invalid_request", "config_overrides must be dictionary-like.")
  performance_timing = _service_request_value(request, "performance_timing", :off)
  run_diagnostics = _service_request_value(request, "run_diagnostics", false)
  run_diagnostics isa Bool || return _service_failure("invalid_request", "run_diagnostics must be boolean.")
  detailed_result_csv = _service_request_value(request, "detailed_result_csv", false)
  detailed_result_csv isa Bool || return _service_failure("invalid_request", "detailed_result_csv must be boolean.")
  detailed_result_csv_semicolon = _service_request_value(request, "detailed_result_csv_semicolon", false)
  detailed_result_csv_semicolon isa Bool || return _service_failure("invalid_request", "detailed_result_csv_semicolon must be boolean.")
  detailed_result_csv_format = _service_request_value(request, "detailed_result_csv_format", nothing)
  if detailed_result_csv && detailed_result_csv_format !== nothing
    detailed_result_csv_format isa AbstractString || return _service_failure("invalid_request", "detailed_result_csv_format must be a string.")
    try
      _resolve_detailed_csv_format(detailed_result_csv_format)
    catch err
      return _service_failure("invalid_request", sprint(showerror, err))
    end
  end
  phases[:request_parse] = _api_elapsed_seconds(request_start)

  requested_run_id = _service_request_value(request, "run_id", nothing)
  run_id = requested_run_id === nothing ? string(uuid4()) : String(requested_run_id)
  _safe_powerflow_run_id(run_id) || return _service_failure("unsafe_run_id", "Unsafe PowerFlow run ID rejected.")
  root = abspath(output_root)
  output_dir = joinpath(root, run_id)
  result = try
    _run_sparlectra_api(
      casefile = casefile,
      config_file = config_file,
      output_dir = output_dir,
      config_overrides = config_overrides,
      performance_timing = performance_timing,
      run_diagnostics = run_diagnostics,
      detailed_result_csv = detailed_result_csv,
      detailed_result_csv_format = detailed_result_csv_format,
      detailed_result_csv_semicolon = detailed_result_csv_semicolon,
      phase_timings = phases,
      run_id = run_id,
      cancellation_token = cancellation_token,
      phase_callback = phase_callback,
    )
  catch err
    err isa PowerFlowAborted && rethrow()
    return _service_failure("execution_error", sprint(showerror, err, catch_backtrace()); run_id = run_id)
  end

  try
    lock(_POWERFLOW_SERVICE_LOCK) do
      _POWERFLOW_SERVICE_RUNS[result.run_id] = result
      _write_powerflow_run_index!(root, result)
    end
  catch err
    lock(_POWERFLOW_SERVICE_LOCK) do
      delete!(_POWERFLOW_SERVICE_RUNS, result.run_id)
    end
    return _service_failure("run_index_error", sprint(showerror, err, catch_backtrace()); run_id = run_id)
  end
  return to_dict(result)
end

function get_powerflow_result(run_id::AbstractString)::Dict{String,Any}
  result = _registered_powerflow_run(run_id)
  result === nothing && return _service_failure("run_not_found", "No PowerFlow run found for run_id $(run_id)."; run_id = run_id)
  return to_dict(result)
end

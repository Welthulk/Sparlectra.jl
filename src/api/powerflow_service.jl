const _POWERFLOW_SERVICE_RUNS = Dict{String,SparlectraApiResult}()
const _POWERFLOW_SERVICE_LOCK = ReentrantLock()

function _service_failure(reason::AbstractString, message::AbstractString; run_id = nothing)::Dict{String,Any}
  failure = Dict{String,Any}(
    "status" => "failed",
    "success" => false,
    "reason" => String(reason),
    "message" => String(message),
  )
  run_id === nothing || (failure["run_id"] = String(run_id))
  return failure
end

function _service_request_value(request::AbstractDict, key::String, default = nothing)
  haskey(request, key) && return request[key]
  symbol_key = Symbol(key)
  haskey(request, symbol_key) && return request[symbol_key]
  return default
end

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

"""
    start_powerflow_run(request::AbstractDict) -> Dict{String,Any}

Start one local PowerFlow service run above [`run_sparlectra_api`](@ref). The
request must provide `casefile`, `config_file`, and `output_root`; optional
`config_overrides` are forwarded to the programmatic API. A unique run ID is
chosen before execution, and all generated files are written beneath
`output_root/run_id`. Public service failures are returned as structured
dictionaries.
"""
function start_powerflow_run(request::AbstractDict)::Dict{String,Any}
  casefile = _service_request_value(request, "casefile")
  casefile isa AbstractString && !isempty(strip(casefile)) || return _service_failure("missing_casefile", "PowerFlow service request requires a nonempty casefile.")
  isfile(casefile) || return _service_failure("missing_casefile", "MATPOWER case file not found: $(abspath(casefile))")

  config_file = _service_request_value(request, "config_file")
  config_file isa AbstractString && !isempty(strip(config_file)) || return _service_failure("missing_config_file", "PowerFlow service request requires a nonempty config_file.")
  isfile(config_file) || return _service_failure("missing_config_file", "Configuration file not found: $(abspath(config_file))")

  output_root = _service_request_value(request, "output_root")
  output_root isa AbstractString && !isempty(strip(output_root)) || return _service_failure("invalid_request", "PowerFlow service request requires a nonempty output_root.")

  config_overrides = _service_request_value(request, "config_overrides", Dict{String,Any}())
  config_overrides isa AbstractDict || return _service_failure("invalid_request", "config_overrides must be dictionary-like.")

  run_id = string(uuid4())
  output_dir = joinpath(abspath(output_root), run_id)
  result = try
    _run_sparlectra_api(
      casefile = casefile,
      config_file = config_file,
      output_dir = output_dir,
      config_overrides = config_overrides,
      run_id = run_id,
    )
  catch err
    return _service_failure("execution_error", sprint(showerror, err, catch_backtrace()); run_id = run_id)
  end

  _register_powerflow_run!(result)
  return to_dict(result)
end

"""
    get_powerflow_result(run_id::AbstractString) -> Dict{String,Any}

Return the transport-safe result metadata for a run started by
[`start_powerflow_run`](@ref), or a structured `run_not_found` failure.
"""
function get_powerflow_result(run_id::AbstractString)::Dict{String,Any}
  result = _registered_powerflow_run(run_id)
  result === nothing && return _service_failure("run_not_found", "No PowerFlow run found for run_id $(run_id)."; run_id = run_id)
  return to_dict(result)
end

"""
    list_powerflow_artifacts(run_id::AbstractString)

List transport-safe artifact metadata for a registered PowerFlow run. Unknown
run IDs return a structured `run_not_found` failure dictionary.
"""
function list_powerflow_artifacts(run_id::AbstractString)
  result = _registered_powerflow_run(run_id)
  result === nothing && return _service_failure("run_not_found", "No PowerFlow run found for run_id $(run_id)."; run_id = run_id)
  return [to_dict(artifact) for artifact in result.artifacts]
end

function _unsafe_artifact_name(name::String)::Bool
  isempty(name) && return true
  isabspath(name) && return true
  occursin('\\', name) && return true
  occursin(r"^[A-Za-z]:", name) && return true
  parts = split(replace(name, '\\' => '/'), '/'; keepempty = true)
  return any(part -> isempty(part) || part == "." || part == "..", parts)
end

function _artifact_belongs_to_run(artifact::SparlectraApiArtifact, output_dir::String)::Bool
  artifact.exists && isfile(artifact.path) || return false
  root = realpath(output_dir)
  path = realpath(artifact.path)
  relative = relpath(path, root)
  return first(splitpath(relative)) != ".."
end

"""
    resolve_powerflow_artifact(run_id::AbstractString, artifact_name::AbstractString)

Resolve an artifact by its metadata `name`, never by an arbitrary filesystem
path. Absolute paths, traversal components, Windows-style paths, missing files,
and artifacts escaping the selected run directory are returned as structured
failures.
"""
function resolve_powerflow_artifact(run_id::AbstractString, artifact_name::AbstractString)
  result = _registered_powerflow_run(run_id)
  result === nothing && return _service_failure("run_not_found", "No PowerFlow run found for run_id $(run_id)."; run_id = run_id)

  name = String(artifact_name)
  _unsafe_artifact_name(name) && return _service_failure("unsafe_artifact_name", "Unsafe artifact name rejected: $(name)"; run_id = run_id)

  index = findfirst(artifact -> artifact.name == name, result.artifacts)
  index === nothing && return _service_failure("artifact_not_found", "No artifact named $(name) belongs to PowerFlow run $(run_id)."; run_id = run_id)
  artifact = result.artifacts[index]
  isfile(artifact.path) || return _service_failure("artifact_not_found", "Artifact $(name) is no longer available for PowerFlow run $(run_id)."; run_id = run_id)
  _artifact_belongs_to_run(artifact, result.output_dir) || return _service_failure("unsafe_artifact_name", "Artifact $(name) does not resolve inside PowerFlow run $(run_id)."; run_id = run_id)
  return artifact
end

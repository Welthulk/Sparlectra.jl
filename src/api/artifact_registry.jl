function _safe_powerflow_run_id(run_id::AbstractString)::Bool
  id = String(run_id)
  return !isempty(id) && !_unsafe_artifact_name(id) && occursin(r"^[A-Za-z0-9][A-Za-z0-9._-]*$", id)
end

function _path_is_within(path::AbstractString, root::AbstractString)::Bool
  normalized_path = normpath(abspath(path))
  normalized_root = normpath(abspath(root))
  relative = relpath(normalized_path, normalized_root)
  relative == "." && return true
  return first(splitpath(relative)) != ".."
end

function _existing_path_is_within(path::AbstractString, root::AbstractString)::Bool
  _path_is_within(path, root) || return false
  ispath(path) || return true
  isdir(root) || return false
  return _path_is_within(realpath(path), realpath(root))
end

function list_powerflow_artifacts(run_id::AbstractString)
  result = _registered_powerflow_run(run_id)
  result === nothing && return _service_failure("run_not_found", "No PowerFlow run found for run_id $(run_id)."; run_id = run_id)
  return [to_dict(artifact) for artifact in collect_sparlectra_api_artifacts(result.output_dir)]
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

  artifacts = collect_sparlectra_api_artifacts(result.output_dir)
  index = findfirst(artifact -> artifact.name == name, artifacts)
  index === nothing && return _service_failure("artifact_not_found", "No artifact named $(name) belongs to PowerFlow run $(run_id)."; run_id = run_id)
  artifact = artifacts[index]
  isfile(artifact.path) || return _service_failure("artifact_not_found", "Artifact $(name) is no longer available for PowerFlow run $(run_id)."; run_id = run_id)
  _artifact_belongs_to_run(artifact, result.output_dir) || return _service_failure("unsafe_artifact_name", "Artifact $(name) does not resolve inside PowerFlow run $(run_id)."; run_id = run_id)
  return artifact
end


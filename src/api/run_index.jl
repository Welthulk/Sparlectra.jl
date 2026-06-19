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

const POWERFLOW_RUN_INDEX_FILENAME = "powerflow_runs_index.json"

function _powerflow_index_path(output_root::AbstractString)::String
  return joinpath(abspath(output_root), POWERFLOW_RUN_INDEX_FILENAME)
end

function _empty_powerflow_run_index()::Dict{String,Any}
  return Dict{String,Any}("schema_version" => _SPARLECTRA_API_SCHEMA_VERSION, "runs" => Any[])
end

function _indexed_run_paths(entry::AbstractDict, output_root::AbstractString)
  run_id = get(entry, "run_id", nothing)
  output_dir = get(entry, "output_dir", nothing)
  result_file = get(entry, "result_file", nothing)
  run_id isa AbstractString && !_unsafe_artifact_name(String(run_id)) || return (valid = false, reason = "unsafe_run_id", run_id = run_id, output_dir = output_dir, result_file = result_file)
  output_dir isa AbstractString || return (valid = false, reason = "invalid_output_dir", run_id = run_id, output_dir = output_dir, result_file = result_file)
  result_file isa AbstractString || return (valid = false, reason = "invalid_result_file", run_id = run_id, output_dir = output_dir, result_file = result_file)

  root = normpath(abspath(output_root))
  directory = normpath(abspath(output_dir))
  result_path = normpath(abspath(result_file))
  expected_directory = joinpath(root, String(run_id))
  directory == expected_directory || return (valid = false, reason = "unsafe_output_dir", run_id = run_id, output_dir = directory, result_file = result_path)
  _existing_path_is_within(directory, root) || return (valid = false, reason = "unsafe_output_dir", run_id = run_id, output_dir = directory, result_file = result_path)
  result_path == joinpath(directory, "result.json") || return (valid = false, reason = "unsafe_result_file", run_id = run_id, output_dir = directory, result_file = result_path)
  _existing_path_is_within(result_path, directory) || return (valid = false, reason = "unsafe_result_file", run_id = run_id, output_dir = directory, result_file = result_path)
  return (valid = true, reason = nothing, run_id = String(run_id), output_dir = directory, result_file = result_path)
end

function _powerflow_run_index_entry(result::SparlectraApiResult)::Dict{String,Any}
  return Dict{String,Any}(
    "run_id" => result.run_id,
    "schema_version" => result.schema_version,
    "status" => String(result.status),
    "success" => result.success,
    "output_dir" => result.output_dir,
    "result_file" => result.result_file,
    "casefile" => result.casefile,
    "config_file" => result.config_file,
    "final_mismatch" => result.final_mismatch,
    "iterations" => result.iterations,
    "timestamp" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
  )
end

function _powerflow_lifecycle_metadata(data::AbstractDict)::Dict{String,Any}
  metadata = get(data, "metadata", Dict{String,Any}())
  lifecycle = metadata isa AbstractDict ? Dict{String,Any}(String(key) => value for (key, value) in metadata) : Dict{String,Any}()
  for key in ("solver_status", "artifact_status", "run_status", "last_phase", "last_heartbeat", "final_outcome")
    haskey(data, key) && (lifecycle[key] = data[key])
  end
  return lifecycle
end

"""
    load_powerflow_run_index(output_root::AbstractString) -> Dict{String,Any}

Load the persistent PowerFlow run index beneath `output_root`. A missing or
unreadable index returns an empty transport-safe index. Invalid top-level
content is reported through `reason` and `message` without throwing.
"""
function load_powerflow_run_index(output_root::AbstractString)::Dict{String,Any}
  index_path = _powerflow_index_path(output_root)
  isfile(index_path) || return _empty_powerflow_run_index()
  try
    parsed = _parse_service_json(read(index_path, String))
    parsed isa AbstractDict || error("PowerFlow run index must contain a JSON object")
    runs = get(parsed, "runs", Any[])
    runs isa AbstractVector || error("PowerFlow run index field 'runs' must be an array")
    return Dict{String,Any}(
      "schema_version" => string(get(parsed, "schema_version", _SPARLECTRA_API_SCHEMA_VERSION)),
      "runs" => Any[run for run in runs],
    )
  catch err
    index = _empty_powerflow_run_index()
    index["reason"] = "invalid_run_index"
    index["message"] = sprint(showerror, err)
    return index
  end
end

function _validated_powerflow_run_entry(entry, output_root::AbstractString)::Dict{String,Any}
  entry isa AbstractDict || return Dict{String,Any}("available" => false, "reason" => "invalid_run_entry")
  listed = Dict{String,Any}(String(key) => _api_transport_value(value) for (key, value) in entry)
  paths = _indexed_run_paths(entry, output_root)
  if !haskey(listed, "timestamp") && paths.valid && isfile(paths.result_file)
    listed["timestamp"] = Dates.format(Dates.unix2datetime(stat(paths.result_file).mtime), "yyyy-mm-dd HH:MM:SS")
  end
  listed["available"] = paths.valid && isdir(paths.output_dir) && isfile(paths.result_file)
  if !paths.valid
    listed["reason"] = paths.reason
  elseif !isdir(paths.output_dir)
    listed["reason"] = "output_dir_not_found"
  elseif !isfile(paths.result_file)
    listed["reason"] = "result_file_not_found"
  end
  return listed
end

"""
    list_powerflow_runs(output_root::AbstractString) -> Vector{Dict{String,Any}}

List runs from the persistent index for a future run-history UI. Each entry has
an `available` flag; missing directories, missing `result.json` files, and
unsafe indexed paths are described by a structured `reason` instead of raising.
"""
function list_powerflow_runs(output_root::AbstractString)::Vector{Dict{String,Any}}
  index = load_powerflow_run_index(output_root)
  runs = [_validated_powerflow_run_entry(entry, output_root) for entry in index["runs"]]
  root = abspath(output_root)
  active_jobs = lock(_POWERFLOW_SERVICE_LOCK) do
    [_webui_job_snapshot(job) for job in values(_POWERFLOW_WEBUI_JOBS) if job["output_root"] == root && get(job, "status", "") in _POWERFLOW_WEBUI_ACTIVE_STATES]
  end
  indexed_ids = Set(string(get(run, "run_id", "")) for run in runs)
  for job in active_jobs
    job["available"] = true
    job["timestamp"] = Dates.format(job["started_at"], "yyyy-mm-dd HH:MM:SS")
    job["run_id"] in indexed_ids || push!(runs, job)
  end
  return runs
end

function _write_powerflow_run_entries!(output_root::AbstractString, runs::AbstractVector)
  root = abspath(output_root)
  mkpath(root)
  contents = Dict{String,Any}("schema_version" => _SPARLECTRA_API_SCHEMA_VERSION, "runs" => runs)
  index_path = _powerflow_index_path(root)
  temporary_path = index_path * ".tmp"
  # Write the run index atomically enough for local Web UI restarts: a complete
  # JSON file replaces the previous index only after serialization succeeds.
  open(temporary_path, "w") do io
    _write_json(io, contents)
    println(io)
  end
  mv(temporary_path, index_path; force = true)
  return index_path
end

function _write_powerflow_run_index!(output_root::AbstractString, result::SparlectraApiResult)
  index = load_powerflow_run_index(output_root)
  runs = Any[entry for entry in index["runs"] if !(entry isa AbstractDict && get(entry, "run_id", nothing) == result.run_id)]
  push!(runs, _powerflow_run_index_entry(result))
  sort!(runs; by = entry -> entry isa AbstractDict ? string(get(entry, "run_id", "")) : "")
  return _write_powerflow_run_entries!(output_root, runs)
end

"""
    refresh_powerflow_run_registry!(output_root::AbstractString) -> Dict{String,Any}

Reload valid PowerFlow service runs from the persistent index under
`output_root` into the in-process registry. Runs that were active when the Web
UI stopped are marked as stale aborted results before registration, while
invalid or missing entries are reported in `unavailable_runs`.
"""
function refresh_powerflow_run_registry!(output_root::AbstractString)::Dict{String,Any}
  root = abspath(output_root)
  index = load_powerflow_run_index(root)
  recovered = SparlectraApiResult[]
  stale_recovered = SparlectraApiResult[]
  unavailable = Dict{String,Any}[]
  if haskey(index, "reason")
    push!(unavailable, Dict{String,Any}(
      "reason" => index["reason"],
      "message" => get(index, "message", nothing),
    ))
  end
  for entry in index["runs"]
    if !(entry isa AbstractDict)
      push!(unavailable, Dict{String,Any}("reason" => "invalid_run_entry"))
      continue
    end
    paths = _indexed_run_paths(entry, root)
    if !paths.valid
      push!(unavailable, Dict{String,Any}("run_id" => _api_transport_value(paths.run_id), "reason" => paths.reason))
      continue
    end
    if !isfile(paths.result_file)
      push!(unavailable, Dict{String,Any}("run_id" => paths.run_id, "reason" => "result_file_not_found"))
      continue
    end
    try
      data = _parse_service_json(read(paths.result_file, String))
      data isa AbstractDict || error("result.json must contain a JSON object")
      if lowercase(string(get(data, "status", ""))) in _POWERFLOW_WEBUI_ACTIVE_STATES
        # A restarted Web UI cannot know whether an active run reached a valid
        # solution, so persist one explicit interrupted status before recovery.
        last_phase = string(get(data, "last_phase", get(data, "current_phase", "unknown")))
        data["status"] = "interrupted_unknown"
        data["success"] = false
        data["solution_available"] = false
        data["reason"] = "webui_restart_lost_live_task"
        data["message"] = "Run state was recovered after Web UI restart. No live solver task is attached anymore; partial artifacts may be available."
        data["solver_status"] = string(get(data, "solver_status", last_phase in ("postprocessing_result", "writing_artifacts", "writing_csv_artifacts", "finalizing_success") ? "completed" : "running"))
        data["artifact_status"] = string(get(data, "artifact_status", startswith(last_phase, "writing") ? "running" : "not_started"))
        data["run_status"] = "interrupted_unknown"
        data["final_outcome"] = "webui_restart_lost_live_task"
        data["last_phase"] = last_phase
        data["last_heartbeat"] = Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ")
        open(paths.result_file, "w") do io
          _write_json(io, data)
          println(io)
        end
        push!(stale_recovered, _reconstruct_powerflow_result(data, paths))
      end
      push!(recovered, _reconstruct_powerflow_result(data, paths))
    catch err
      push!(unavailable, Dict{String,Any}(
        "run_id" => paths.run_id,
        "reason" => "invalid_result_file",
        "message" => sprint(showerror, err),
      ))
    end
  end

  lock(_POWERFLOW_SERVICE_LOCK) do
    for (run_id, result) in collect(_POWERFLOW_SERVICE_RUNS)
      _path_is_within(result.output_dir, root) && delete!(_POWERFLOW_SERVICE_RUNS, run_id)
    end
    for result in recovered
      _POWERFLOW_SERVICE_RUNS[result.run_id] = result
    end
  end
  return Dict{String,Any}(
    "status" => isempty(unavailable) ? "succeeded" : "partial",
    "success" => isempty(unavailable),
    "loaded_runs" => [result.run_id for result in recovered],
    "runs" => recovered,
    "stale_recovered_runs" => stale_recovered,
    "unavailable_runs" => unavailable,
  )
end

"""
    delete_powerflow_run(run_id::AbstractString; output_root::AbstractString) -> Dict{String,Any}

Delete one registered PowerFlow run beneath `output_root`. The run ID must be a
safe index entry whose directory is exactly `<output_root>/<run_id>`; arbitrary
paths and unregistered directories are never removed.
"""
function delete_powerflow_run(run_id::AbstractString; output_root::AbstractString)::Dict{String,Any}
  id = String(run_id)
  _unsafe_artifact_name(id) && return _service_failure("unsafe_run_id", "Unsafe PowerFlow run ID rejected."; run_id = id)
  root = abspath(output_root)
  index = load_powerflow_run_index(root)
  entry_index = findfirst(entry -> entry isa AbstractDict && get(entry, "run_id", nothing) == id, index["runs"])
  entry_index === nothing && return _service_failure("run_not_found", "PowerFlow run is not registered for this output root."; run_id = id)
  entry = index["runs"][entry_index]
  paths = _indexed_run_paths(entry, root)
  paths.valid || return _service_failure(paths.reason, "PowerFlow run has an unsafe indexed path and was not deleted."; run_id = id)

  try
    ispath(paths.output_dir) && rm(paths.output_dir; recursive = true)
  catch err
    return _service_failure("delete_failed", sprint(showerror, err); run_id = id)
  end
  remaining = Any[entry for (entry_position, entry) in enumerate(index["runs"]) if entry_position != entry_index]
  _write_powerflow_run_entries!(root, remaining)
  lock(_POWERFLOW_SERVICE_LOCK) do
    result = get(_POWERFLOW_SERVICE_RUNS, id, nothing)
    result === nothing || (_path_is_within(result.output_dir, root) && delete!(_POWERFLOW_SERVICE_RUNS, id))
  end
  return Dict{String,Any}("status" => "succeeded", "success" => true, "run_id" => id)
end

"""
    delete_all_powerflow_runs(; output_root::AbstractString) -> Dict{String,Any}

Delete all safely registered PowerFlow run directories beneath `output_root`.
Entries that cannot be validated or removed remain indexed and are reported in
`failed_runs`; deletion never follows an indexed path outside the output root.
"""
function delete_all_powerflow_runs(; output_root::AbstractString)::Dict{String,Any}
  root = abspath(output_root)
  index = load_powerflow_run_index(root)
  deleted = String[]
  failed = Dict{String,Any}[]
  remaining = Any[]
  for entry in index["runs"]
    if !(entry isa AbstractDict)
      push!(remaining, entry)
      push!(failed, Dict{String,Any}("reason" => "invalid_run_entry"))
      continue
    end
    paths = _indexed_run_paths(entry, root)
    if !paths.valid
      push!(remaining, entry)
      push!(failed, Dict{String,Any}("run_id" => _api_transport_value(paths.run_id), "reason" => paths.reason))
      continue
    end
    try
      ispath(paths.output_dir) && rm(paths.output_dir; recursive = true)
      push!(deleted, paths.run_id)
    catch err
      push!(remaining, entry)
      push!(failed, Dict{String,Any}("run_id" => paths.run_id, "reason" => "delete_failed", "message" => sprint(showerror, err)))
    end
  end
  _write_powerflow_run_entries!(root, remaining)
  lock(_POWERFLOW_SERVICE_LOCK) do
    for run_id in deleted
      result = get(_POWERFLOW_SERVICE_RUNS, run_id, nothing)
      result === nothing || (_path_is_within(result.output_dir, root) && delete!(_POWERFLOW_SERVICE_RUNS, run_id))
    end
  end
  return Dict{String,Any}(
    "status" => isempty(failed) ? "succeeded" : "partial",
    "success" => isempty(failed),
    "deleted_runs" => deleted,
    "failed_runs" => failed,
  )
end

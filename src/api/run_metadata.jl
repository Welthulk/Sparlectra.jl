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

# file: src/api/run_metadata.jl
#
# This file contains API run-result and metadata assembly helpers. Keeping them
# separate lets run_api.jl focus on orchestration.

function _runtime_metadata_payload(; case_path::AbstractString, lifecycle::AbstractDict = Dict{String,Any}())
  payload = Dict{String,Any}(
    "runtime_request" => Dict{String,Any}(
      "casefile" => basename(case_path),
      "casefile_path" => String(case_path),
    ),
  )
  haskey(lifecycle, "configured_default_casefile") && (payload["configured_default_casefile"] = lifecycle["configured_default_casefile"])
  haskey(lifecycle, "runtime_casefile") && (payload["runtime_casefile"] = lifecycle["runtime_casefile"])
  for key in ("solver_status", "artifact_status", "run_status", "last_phase", "last_heartbeat", "final_outcome")
    haskey(lifecycle, key) && (payload[key] = lifecycle[key])
  end
  return payload
end

function _write_run_metadata_artifact(output_path::AbstractString; case_path::AbstractString, lifecycle::AbstractDict = Dict{String,Any}())::String
  artifact = joinpath(output_path, "run_metadata.yaml")
  _write_yaml_file(artifact, _runtime_metadata_payload(case_path = case_path, lifecycle = lifecycle))
  return artifact
end

function _api_result(;
  run_id::String = string(uuid4()),
  schema_version::String = _SPARLECTRA_API_SCHEMA_VERSION,
  status::Symbol,
  success::Bool,
  converged = nothing,
  solution_available::Bool = false,
  iterations = nothing,
  final_mismatch = nothing,
  reason = nothing,
  message = nothing,
  casefile = nothing,
  config_file = nothing,
  output_dir::String,
  logfile = nothing,
  result_file = nothing,
  artifacts = SparlectraApiArtifact[],
  service_phase_timings = Dict{String,Any}[],
  metadata = Dict{String,Any}(),
  raw_result = nothing,
)
  timings = Dict{String,Any}[Dict{String,Any}(String(key) => value for (key, value) in timing) for timing in service_phase_timings]
  return SparlectraApiResult(run_id, schema_version, status, success, converged, solution_available, iterations, final_mismatch, reason, message, casefile, config_file, output_dir, logfile, result_file, artifacts, timings, Dict{String,Any}(String(key) => value for (key, value) in metadata), raw_result)
end

function _write_api_result_file(result::SparlectraApiResult)
  result.result_file === nothing || write(result.result_file, to_json(result), '\n')
  return result
end

function _api_result_with_artifacts(result::SparlectraApiResult, artifacts::Vector{SparlectraApiArtifact})::SparlectraApiResult
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
    artifacts = artifacts,
    service_phase_timings = result.service_phase_timings,
    metadata = result.metadata,
    raw_result = result.raw_result,
  )
end

function _refresh_api_artifacts(result::SparlectraApiResult)::SparlectraApiResult
  return _api_result_with_artifacts(result, collect_sparlectra_api_artifacts(result.output_dir))
end

function _write_api_result_file_with_discovered_artifacts(result::SparlectraApiResult)::SparlectraApiResult
  result.result_file === nothing && return result
  discovered = collect_sparlectra_api_artifacts(result.output_dir)
  refreshed = _api_result_with_artifacts(result, discovered)
  _write_api_result_file(refreshed)
  return refreshed
end

function _finalize_api_result(result::SparlectraApiResult)::SparlectraApiResult
  _write_api_result_file(result)
  refreshed = _refresh_api_artifacts(result)
  if haskey(refreshed.metadata, "artifact_status") && refreshed.metadata["artifact_status"] == "not_started" && !isempty(refreshed.artifacts)
    refreshed = _api_result(
      run_id = refreshed.run_id,
      status = refreshed.status,
      success = refreshed.success,
      converged = refreshed.converged,
      solution_available = refreshed.solution_available,
      iterations = refreshed.iterations,
      final_mismatch = refreshed.final_mismatch,
      reason = refreshed.reason,
      message = refreshed.message,
      casefile = refreshed.casefile,
      config_file = refreshed.config_file,
      output_dir = refreshed.output_dir,
      logfile = refreshed.logfile,
      result_file = refreshed.result_file,
      artifacts = refreshed.artifacts,
      service_phase_timings = refreshed.service_phase_timings,
      metadata = merge(refreshed.metadata, Dict{String,Any}("artifact_status" => "completed")),
      raw_result = refreshed.raw_result,
    )
  end
  _write_api_result_file(refreshed)
  return _write_api_result_file_with_discovered_artifacts(refreshed)
end

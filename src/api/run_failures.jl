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

# file: src/api/run_failures.jl
#
# This file contains API failure construction helpers. It keeps backend failure
# metadata consistent across CLI/API/Web UI callers.

function _api_failure(reason::String, message::String; run_id::String = string(uuid4()), casefile, config_file, output_dir::String, logfile::String, result_file::String, service_phase_timings = Dict{String,Any}[], metadata = Dict{String,Any}())::SparlectraApiResult
  open(logfile, "a") do io
    println(io, "Sparlectra API failure: ", reason)
    println(io, message)
  end
  generated_artifacts = collect_sparlectra_api_artifacts(output_dir)
  artifact_status = isempty(generated_artifacts) ? "not_started" : "completed"
  failure_metadata = merge(
    Dict{String,Any}(
      "failure_reason" => reason,
      "artifact_status" => artifact_status,
      "solver_status" => "failed",
      "service_status" => "failed",
      "run_status" => "failed",
      "numerical_status" => "failed",
    ),
    Dict{String,Any}(String(key) => value for (key, value) in metadata),
  )
  if artifact_status == "completed"
    failure_metadata["artifact_status"] = "completed"
  end
  result = _api_result(run_id = run_id, status = :failed, success = false, reason = reason, message = message, casefile = casefile, config_file = config_file, output_dir = output_dir, logfile = logfile, result_file = result_file, service_phase_timings = service_phase_timings, metadata = failure_metadata)
  return _finalize_api_result(result)
end

function _api_execution_failure(reason::String, message::String; run_id::String, casefile, config_file, output_dir::String, logfile::String, result_file::String, phase_recorder::PowerFlowPhaseTimingRecorder, performance_timing = :off, metadata = Dict{String,Any}(), total_start_ns::Union{Nothing,UInt64} = nothing)::SparlectraApiResult
  _complete_active_phase!(phase_recorder, "failed")
  _start_service_phase!(phase_recorder, "finalizing_failed")
  if total_start_ns === nothing
    _complete_active_phase!(phase_recorder, "failed")
  else
    _finalize_service_timings!(phase_recorder, total_start_ns; status = "failed")
  end
  timing_mode = try
    _api_timing_mode(performance_timing)
  catch
    :off
  end
  if timing_mode !== :off
    open(joinpath(output_dir, "performance.log"), "a") do io
      println(io, "Sparlectra API failure: ", reason)
      println(io, message)
      println(io)
      _write_service_phase_summary(io, phase_recorder.timings)
    end
  end
  return _api_failure(reason, message; run_id, casefile, config_file, output_dir, logfile, result_file, service_phase_timings = phase_recorder.timings, metadata = metadata)
end

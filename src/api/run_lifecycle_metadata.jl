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

# API run lifecycle metadata helpers.
#
# These helpers centralize Web UI-visible lifecycle fields. Keep key names and
# status strings stable because the service index, Last Errors, and clients use
# this metadata to explain completed, failed, and non-converged runs.

function _current_iteration_lifecycle_metadata(rect_status)::Dict{String,Any}
  rect_status === nothing && return Dict{String,Any}(
    "current_iteration_enabled" => false,
    "current_iteration_attempted" => false,
    "current_iteration_accepted" => false,
    "current_iteration_iterations" => 0,
    "current_iteration_initial_mismatch" => NaN,
    "current_iteration_final_mismatch" => NaN,
    "current_iteration_reason" => "unavailable",
    "current_iteration_artifact" => "",
  )
  return Dict{String,Any}(
    "current_iteration_enabled" => getproperty(rect_status, :current_iteration_enabled),
    "current_iteration_attempted" => getproperty(rect_status, :current_iteration_attempted),
    "current_iteration_accepted" => getproperty(rect_status, :current_iteration_accepted),
    "current_iteration_iterations" => getproperty(rect_status, :current_iteration_iterations),
    "current_iteration_initial_mismatch" => getproperty(rect_status, :current_iteration_initial_mismatch),
    "current_iteration_final_mismatch" => getproperty(rect_status, :current_iteration_final_mismatch),
    "current_iteration_reason" => String(getproperty(rect_status, :current_iteration_reason)),
    "current_iteration_artifact" => getproperty(rect_status, :current_iteration_artifact),
  )
end

function _build_success_lifecycle_metadata(raw_result::SparlectraRunResult, config::SparlectraConfig; numerical_success::Bool, final_outcome::Dict{String,Any}, csv_export_status::AbstractString, csv_export_skip_reason, csv_export_error, csv_artifacts::Vector{String}, detailed_result_csv::Bool, config_overrides, casefile, config_file, performance_timing, run_diagnostics::Bool, detailed_result_csv_format, qlimit_metadata::AbstractDict, csv_timing_metadata::AbstractDict)::Dict{String,Any}
  rect_status = raw_result.net === nothing ? nothing : rectangular_pf_status(raw_result.net)
  current_iteration_metadata = _current_iteration_lifecycle_metadata(rect_status)
  classic_outer_loop_passes = rect_status !== nothing && hasproperty(rect_status, :matpower_outer_iterations) ? rect_status.matpower_outer_iterations : 0
  pv_to_pq_events = raw_result.net === nothing ? 0 : length(raw_result.net.qLimitLog)
  active_set_events = config.powerflow.qlimits.enforcement_mode === :active_set ? pv_to_pq_events : 0
  metadata = merge(Dict{String,Any}(
    "solver_status" => "completed",
    "service_status" => "completed",
    "numerical_status" => numerical_success ? "converged" : "not_converged",
    "artifact_status" => csv_export_error === nothing ? (csv_export_status == "partial" ? "partial" : "completed") : "failed",
    "run_status" => numerical_success ? "completed" : "completed_nonconverged",
    "last_phase" => "finalizing_success",
    "last_heartbeat" => Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ"),
    "final_outcome" => final_outcome,
    "detailed_result_csv_status" => csv_export_status,
    "detailed_result_csv_skip_reason" => csv_export_skip_reason,
    "detailed_result_csv_error" => csv_export_error,
    "detailed_result_csv_artifacts" => csv_artifacts,
    "detailed_result_csv_solution_quality" => detailed_result_csv && csv_export_skip_reason === nothing ? _csv_solution_quality(raw_result) : nothing,
    "detailed_result_csv_warning" => detailed_result_csv && csv_export_skip_reason === nothing && !numerical_success ? "values are from the last Newton iterate and are not a converged power-flow solution" : nothing,
    "q_limit_active_set_events" => active_set_events,
    "q_limit_pv_to_pq_events" => pv_to_pq_events,
    "q_limit_classic_outer_loop_passes" => classic_outer_loop_passes,
    "webui_request_settings" => merge(Dict{String,Any}(String(key) => value for (key, value) in config_overrides), Dict{String,Any}(
      "casefile" => casefile,
      "config_file" => config_file,
      "performance_timing" => String(performance_timing),
      "run_diagnostics" => run_diagnostics,
      "detailed_result_csv" => detailed_result_csv,
      "detailed_result_csv_format" => detailed_result_csv_format === nothing ? "technical" : String(detailed_result_csv_format),
    )),
  ), qlimit_metadata, current_iteration_metadata)
  # Partial CSV exports are still successful API runs, but the Web UI needs the
  # partial file error in the stable lifecycle field used by Last Errors.
  haskey(csv_timing_metadata, :partial_error) && (metadata["detailed_result_csv_error"] = csv_timing_metadata[:partial_error])
  return metadata
end

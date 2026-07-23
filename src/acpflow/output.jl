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

function _sparlectra_result_mode(net::Net, out_cfg::OutputConfig)::Symbol
  mode = out_cfg.logfile_results === :compact ? :summary : out_cfg.logfile_results
  return length(net.nodeVec) >= out_cfg.result_table_large_case_threshold_buses && out_cfg.logfile_results != :full ? out_cfg.result_table_large_case_mode : mode
end

function _postprocess_sparlectra_result!(result::SparlectraRunResult, cfg::SparlectraConfig; emit_output::Bool = true)
  result.solution_available || return result
  phase_callback = result.performance_profile isa AbstractDict ? get(result.performance_profile, :phase_callback, phase -> nothing) : phase -> nothing
  phase_callback("postprocessing_result")
  if cfg.powerflow.solver === :dc
    # DC branch flows are already final (lossless, computed directly from
    # theta) — calcNetLosses!/calcLinkFlowsKCL! are AC-π-model-based and
    # would silently overwrite them, and printACPFlowResults would error
    # against DC's empty AC diagnostics.
    out_cfg = cfg.output
    if emit_output && out_cfg.logfile_results !== :off
      _perf_profile_time!(result.performance_profile, :result_output) do
        printDcPowerFlowResults(result.net, result.elapsed_s)
      end
    end
    return result
  end
  _perf_profile_time!(result.performance_profile, :postprocess_losses_and_flows) do
    calcNetLosses!(result.net)
    calcLinkFlowsKCL!(result.net)
  end
  out_cfg = cfg.output
  if emit_output && out_cfg.logfile_results !== :off
    max_rows = out_cfg.result_table_max_rows > 0 ? out_cfg.result_table_max_rows : nothing
    _perf_profile_time!(result.performance_profile, :result_output) do
      printACPFlowResults(result.net, result.elapsed_s, result.iterations, cfg.powerflow.tol, false, ""; converged = result.final_converged, solver = result.diagnostics.solver, solver_time_s = result.solver_elapsed_s, result_mode = _sparlectra_result_mode(result.net, out_cfg), max_rows = max_rows)
    end
  end
  return result
end

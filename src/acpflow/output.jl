function _sparlectra_result_mode(net::Net, out_cfg::OutputConfig)::Symbol
  mode = out_cfg.logfile_results === :compact ? :summary : out_cfg.logfile_results
  return length(net.nodeVec) >= out_cfg.result_table_large_case_threshold_buses && out_cfg.logfile_results != :full ? out_cfg.result_table_large_case_mode : mode
end

function _postprocess_sparlectra_result!(result::SparlectraRunResult, cfg::SparlectraConfig; emit_output::Bool = true)
  result.solution_available || return result
  phase_callback = result.performance_profile isa AbstractDict ? get(result.performance_profile, :phase_callback, phase -> nothing) : phase -> nothing
  phase_callback("postprocessing_result")
  _perf_profile_time!(result.performance_profile, :postprocess_losses_and_flows) do
    calcNetLosses!(result.net)
    calcLinkFlowsKCL!(result.net)
  end
  out_cfg = cfg.output
  if emit_output && out_cfg.logfile_results !== :off
    max_rows = out_cfg.result_table_max_rows > 0 ? out_cfg.result_table_max_rows : nothing
    _perf_profile_time!(result.performance_profile, :result_output) do
      printACPFlowResults(result.net, result.elapsed_s, result.iterations, cfg.powerflow.tol, false, ""; converged = result.final_converged, solver = result.method, solver_time_s = result.solver_elapsed_s, result_mode = _sparlectra_result_mode(result.net, out_cfg), max_rows = max_rows)
    end
  end
  return result
end

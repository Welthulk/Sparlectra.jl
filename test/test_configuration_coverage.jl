using Sparlectra
using Test

function _leaf_paths(x; prefix = "")
  paths = String[]
  if x isa AbstractDict
    for key in sort!(collect(keys(x)); by = k -> String(k))
      k = String(key)
      child = x[key]
      child_prefix = isempty(prefix) ? k : string(prefix, ".", k)
      append!(paths, _leaf_paths(child; prefix = child_prefix))
    end
  else
    push!(paths, prefix)
  end
  return paths
end

function run_configuration_coverage_tests()
  @testset "Configuration YAML leaf coverage" begin
    yaml = Sparlectra.load_yaml_dict(joinpath(@__DIR__, "..", "src", "configuration.yaml.example"))
    leaves = Set(_leaf_paths(yaml))

    mapped_or_reserved = Set([
      "power_flow.method", "power_flow.sparse", "power_flow.flatstart", "power_flow.tol", "power_flow.max_iter", "power_flow.autodamp", "power_flow.autodamp_min",
      "power_flow.start_mode.angle_mode", "power_flow.start_mode.voltage_mode", "power_flow.start_mode.start_projection", "power_flow.start_mode.try_dc_start", "power_flow.start_mode.try_blend_scan", "power_flow.start_mode.branch_guard", "power_flow.start_mode.measure_candidates", "power_flow.start_mode.accept_unmeasured_dc_start", "power_flow.start_mode.reuse_import_data", "power_flow.start_mode.blend_lambdas", "power_flow.start_mode.dc_angle_limit_deg",
      "power_flow.qlimits.enabled", "power_flow.qlimits.start_iter", "power_flow.qlimits.start_mode", "power_flow.qlimits.auto_q_delta_pu", "power_flow.qlimits.hysteresis_pu", "power_flow.qlimits.cooldown_iters", "power_flow.qlimits.trace_buses", "power_flow.qlimits.lock_pv_to_pq_buses",
      "power_flow.qlimits.guard.enabled", "power_flow.qlimits.guard.min_q_range_pu", "power_flow.qlimits.guard.narrow_range_mode", "power_flow.qlimits.guard.zero_range_mode", "power_flow.qlimits.guard.violation_mode", "power_flow.qlimits.guard.violation_threshold_pu", "power_flow.qlimits.guard.max_switches", "power_flow.qlimits.guard.max_remaining_violations", "power_flow.qlimits.guard.accept_bounded_violations", "power_flow.qlimits.guard.freeze_after_repeated_switching", "power_flow.qlimits.guard.log",
      "state_estimation.enabled", "state_estimation.method", "state_estimation.sparse", "state_estimation.tol", "state_estimation.max_iter", "state_estimation.flatstart", "state_estimation.jac_eps", "state_estimation.update_net", "state_estimation.observability.enabled",
      "matpower_import.case", "matpower_import.auto_profile", "matpower_import.auto_profile_log", "matpower_import.pv_voltage_source", "matpower_import.pv_voltage_mismatch_tol_pu", "matpower_import.compare_voltage_reference", "matpower_import.bus_shunt_model", "matpower_import.shift_unit", "matpower_import.shift_sign", "matpower_import.ratio", "matpower_import.enable_pq_gen_controllers", "matpower_import.preallocate_network", "matpower_import.preallocate_min_buses",
      "performance.enabled", "performance.level", "performance.print_to_console", "performance.write_to_logfile", "performance.show_allocations", "performance.show_iteration_table", "performance.compact_logging", "performance.skip_reference_comparison", "performance.skip_expensive_diagnostics", "performance.skip_branch_neighborhood_report", "performance.max_diagnostic_rows",
      "runtime.print_thread_config", "runtime.julia_threads", "runtime.blas_threads",
      "diagnostics.log_effective_config", "diagnostics.matpower_reference", "diagnostics.branch_shift_conventions", "diagnostics.negative_branch_impedance", "diagnostics.pv_voltage_references", "diagnostics.residual_clusters", "diagnostics.nodal_balance_breakdown", "diagnostics.branch_neighborhood", "diagnostics.compact_console", "diagnostics.detailed_log",
      "output.console_summary", "output.console_auto_profile", "output.console_diagnostics", "output.console_q_limit_events", "output.console_max_rows", "output.logfile_results", "output.result_table_max_rows", "output.result_table_large_case_threshold_buses", "output.result_table_large_case_mode", "output.logfile_diagnostics", "output.logfile_performance", "output.logfile_warnings",
      "benchmark.enabled", "benchmark.methods", "benchmark.seconds", "benchmark.samples", "benchmark.show_once", "benchmark.show_once_output", "benchmark.show_once_max_nodes",
      "extensions.reserved",
    ])

    uncovered = setdiff(leaves, mapped_or_reserved)
    @test isempty(uncovered)
  end
end

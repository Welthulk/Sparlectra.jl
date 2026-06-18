using Sparlectra
using Test
using Logging

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

function _canonical_yaml_leaf_keys()
  yaml = Sparlectra.load_yaml_dict(joinpath(@__DIR__, "..", "src", "configuration.yaml.example"))
  return Set(_leaf_paths(yaml))
end

function _auto_profile_shift_case()
  base = 100.0
  bus = [
    1.0 3.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 110.0 1.0 1.1 0.9
    2.0 1.0 0.0 0.0 0.0 0.0 1.0 1.0 -5.0 110.0 1.0 1.1 0.9
  ]
  branch = [1.0 2.0 0.01 0.10 0.0 999.0 999.0 999.0 1.0 10.0 1.0 -360.0 360.0]
  ybus = Sparlectra.MatpowerIO.build_ybus_matpower(bus, branch, base; matpower_shift_unit = :deg, matpower_shift_sign = 1.0, matpower_ratio = :normal)
  v = bus[:, 8] .* cis.(bus[:, 9] .* (pi / 180.0))
  scalc = v .* conj.(ybus * v)
  bus[2, 3] = -real(scalc[2]) * base
  bus[2, 4] = -imag(scalc[2]) * base
  gen = [1.0 0.0 0.0 999.0 -999.0 1.0 100.0 1.0 999.0 0.0]
  return Sparlectra.MatpowerIO.MatpowerCase("auto_profile_shift", base, bus, gen, branch, nothing, nothing)
end

function _auto_profile_pv_mismatch_case()
  mpc = _auto_profile_shift_case()
  bus = copy(mpc.bus)
  bus[2, 2] = 2.0
  bus[2, 8] = 1.00
  gen = [
    1.0 0.0 0.0 999.0 -999.0 1.0 100.0 1.0 999.0 0.0
    2.0 0.0 0.0 999.0 -999.0 1.04 100.0 1.0 999.0 0.0
  ]
  return Sparlectra.MatpowerIO.MatpowerCase("auto_profile_pv", mpc.baseMVA, bus, gen, mpc.branch, nothing, nothing)
end

function run_configuration_coverage_tests()
  @testset "Configuration YAML key coverage" begin
    leaves = _canonical_yaml_leaf_keys()

    mapped_keys = Set([
      "power_flow.method", "power_flow.flatstart", "power_flow.tol", "power_flow.max_iter", "power_flow.autodamp", "power_flow.autodamp_min", "power_flow.wrong_branch_detection", "power_flow.wrong_branch_rescue", "power_flow.wrong_branch_min_vm_pu", "power_flow.wrong_branch_max_vm_pu", "power_flow.wrong_branch_max_angle_spread_deg", "power_flow.wrong_branch_max_branch_angle_deg", "power_flow.wrong_branch_min_low_vm_count", "power_flow.wrong_branch_rescue_max_attempts", "power_flow.rectangular_workspace_reuse", "power_flow.rectangular_preallocate_workspace", "power_flow.rectangular_workspace_min_buses",
      "power_flow.start_mode.angle_mode", "power_flow.start_mode.voltage_mode", "power_flow.start_mode.profile_source", "power_flow.start_mode.start_projection", "power_flow.start_mode.try_dc_start", "power_flow.start_mode.try_blend_scan", "power_flow.start_mode.branch_guard", "power_flow.start_mode.measure_candidates", "power_flow.start_mode.accept_unmeasured_dc_start", "power_flow.start_mode.reuse_import_data", "power_flow.start_mode.blend_lambdas", "power_flow.start_mode.dc_angle_limit_deg",
      "power_flow.qlimits.enabled", "power_flow.qlimits.start_iter", "power_flow.qlimits.start_mode", "power_flow.qlimits.auto_q_delta_pu", "power_flow.qlimits.hysteresis_pu", "power_flow.qlimits.cooldown_iters", "power_flow.qlimits.trace_buses", "power_flow.qlimits.lock_pv_to_pq_buses",
      "power_flow.qlimits.guard.enabled", "power_flow.qlimits.guard.min_q_range_pu", "power_flow.qlimits.guard.narrow_range_mode", "power_flow.qlimits.guard.zero_range_mode", "power_flow.qlimits.guard.violation_mode", "power_flow.qlimits.guard.violation_threshold_pu", "power_flow.qlimits.guard.max_switches", "power_flow.qlimits.guard.max_remaining_violations", "power_flow.qlimits.guard.accept_bounded_violations", "power_flow.qlimits.guard.freeze_after_repeated_switching", "power_flow.qlimits.guard.log",
      "state_estimation.enabled", "state_estimation.method", "state_estimation.tol", "state_estimation.max_iter", "state_estimation.flatstart", "state_estimation.jac_eps", "state_estimation.update_net", "state_estimation.observability.enabled",
      "matpower_import.case", "matpower_import.cases", "matpower_import.auto_profile", "matpower_import.auto_profile_log", "matpower_import.pv_voltage_source", "matpower_import.pv_voltage_mismatch_tol_pu", "matpower_import.compare_voltage_reference", "matpower_import.bus_shunt_model", "matpower_import.shift_unit", "matpower_import.shift_sign", "matpower_import.ratio", "matpower_import.enable_pq_gen_controllers", "matpower_import.preallocate_network", "matpower_import.preallocate_min_buses",
      "performance.enabled", "performance.level", "performance.print_to_console", "performance.write_to_logfile", "performance.show_allocations", "performance.show_iteration_table", "performance.compact_logging", "performance.representative_warmup_runs", "performance.compare_cold_warm", "performance.skip_reference_comparison", "performance.skip_expensive_diagnostics", "performance.skip_branch_neighborhood_report", "performance.max_diagnostic_rows",
      "runtime.print_thread_config", "runtime.julia_threads", "runtime.blas_threads",
      "diagnostics.log_effective_config", "diagnostics.console_summary", "diagnostics.console_auto_profile", "diagnostics.console_diagnostics", "diagnostics.console_q_limit_events", "diagnostics.console_max_rows", "diagnostics.logfile_diagnostics",
      "output.console_summary", "output.console_auto_profile", "output.console_diagnostics", "output.console_q_limit_events", "output.console_max_rows", "output.logfile_results", "output.result_table_max_rows", "output.result_table_large_case_threshold_buses", "output.result_table_large_case_mode", "output.detailed_result_csv_write_mode", "output.detailed_result_csv_exporter", "output.detailed_result_csv_direct_threshold_buses", "output.detailed_result_csv_buffer_initial_bytes", "output.detailed_result_csv_buffer_max_bytes", "output.detailed_result_csv_streaming_threshold_rows", "output.logfile_diagnostics", "output.logfile_performance", "output.logfile_warnings",
      "benchmark.enabled", "benchmark.methods", "benchmark.seconds", "benchmark.samples", "benchmark.show_once", "benchmark.show_once_output", "benchmark.show_once_max_nodes",
      "control.enabled", "control.max_outer_iterations", "control.trace", "control.log_iterations", "control.stop_on_pf_failure", "control.controllers",
    ])
    reserved_keys = Set(["extensions.reserved"])
    mapped_or_reserved = union(mapped_keys, reserved_keys)

    @test isempty(setdiff(leaves, mapped_or_reserved))
    @test isempty(setdiff(reserved_keys, leaves))

    expected_consumers = Dict(
      "power_flow.tol" => :PowerFlowConfig,
      "power_flow.max_iter" => :PowerFlowConfig,
      "power_flow.autodamp" => :PowerFlowConfig,
      "power_flow.autodamp_min" => :PowerFlowConfig,
      "power_flow.wrong_branch_detection" => :PowerFlowConfig,
      "power_flow.start_mode.angle_mode" => :StartModeConfig,
      "power_flow.start_mode.voltage_mode" => :StartModeConfig,
      "power_flow.start_mode.profile_source" => :StartModeConfig,
      "power_flow.qlimits.start_iter" => :QLimitConfig,
      "matpower_import.shift_sign" => :MatpowerImportConfig,
      "matpower_import.cases" => :MatpowerImportConfig,
      "output.logfile_performance" => :OutputConfig,
      "benchmark.enabled" => :BenchmarkConfig,
      "state_estimation.method" => :StateEstimationConfig,
      "runtime.julia_threads" => :RuntimeConfig,
      "diagnostics.console_diagnostics" => :DiagnosticsConfig,
      "diagnostics.logfile_diagnostics" => :DiagnosticsConfig,
      "extensions.reserved" => :Reserved,
    )
    @test all(haskey(expected_consumers, key) for key in keys(expected_consumers))
    @test expected_consumers["extensions.reserved"] === :Reserved
  end

  @testset "Configuration forwarding with sentinel values" begin
    cfgfile = tempname() * ".yaml"
    write(cfgfile, """
power_flow:
  tol: 1.0e-7
  max_iter: 17
  autodamp: true
  autodamp_min: 0.17
  wrong_branch_detection: fail
  wrong_branch_rescue: true
  wrong_branch_min_vm_pu: 0.65
  wrong_branch_rescue_max_attempts: 1
  start_mode:
    angle_mode: dc
    voltage_mode: profile_blend
    profile_source: matpower_reference
    start_projection: true
    try_dc_start: false
    try_blend_scan: false
    blend_lambdas: [0.11, 0.22]
    dc_angle_limit_deg: 33.0
  qlimits:
    enabled: true
    start_iter: 4
    cooldown_iters: 6
matpower_import:
  shift_unit: rad
  shift_sign: -1.0
  ratio: reciprocal
output:
  logfile_performance: full
benchmark:
  enabled: false
""")
    cfg = Sparlectra.load_sparlectra_config(cfgfile; reload = true)
    @test cfg.powerflow.tol == 1.0e-7
    @test cfg.powerflow.max_iter == 17
    @test cfg.powerflow.autodamp === true
    @test cfg.powerflow.autodamp_min == 0.17
    @test cfg.powerflow.wrong_branch_detection === :fail
    @test cfg.powerflow.wrong_branch_rescue === true
    @test cfg.powerflow.wrong_branch_min_vm_pu == 0.65
    @test cfg.powerflow.wrong_branch_rescue_max_attempts == 1
    @test cfg.powerflow.start_mode.angle_mode === :dc
    @test cfg.powerflow.start_mode.voltage_mode === :profile_blend
    @test cfg.powerflow.start_mode.profile_source === :matpower_reference
    @test cfg.powerflow.start_mode.start_projection === true
    @test cfg.powerflow.start_mode.try_dc_start === false
    @test cfg.powerflow.start_mode.try_blend_scan === false
    @test cfg.powerflow.start_mode.blend_lambdas == [0.11, 0.22]
    @test cfg.powerflow.start_mode.dc_angle_limit_deg == 33.0
    @test cfg.powerflow.qlimits.start_iter == 4
    @test cfg.powerflow.qlimits.cooldown_iters == 6
    @test cfg.matpower.shift_unit === :rad
    @test cfg.matpower.shift_sign == -1.0
    @test cfg.matpower.ratio === :reciprocal
    @test cfg.output.logfile_performance === :full
    @test cfg.benchmark.enabled === false
  end

  @testset "MATPOWER auto-profile decision rules" begin
    mpc = _auto_profile_shift_case()
    cfg = Sparlectra.SparlectraConfig(Dict(
      "matpower_import" => Dict("auto_profile" => "off", "shift_unit" => "rad"),
    ))
    off = Sparlectra.matpower_import_auto_profile(mpc, cfg; mode = :off)
    @test isempty(off.rows)
    @test off.config.matpower.shift_unit === :rad

    rec = Sparlectra.matpower_import_auto_profile(mpc, cfg; mode = :recommend)
    @test rec.config.matpower.shift_unit === :rad
    shift_unit_row = only(row for row in rec.rows if row.option == "matpower_import.shift_unit")
    @test shift_unit_row.recommended == "deg"
    @test shift_unit_row.action === :recommend

    applied = Sparlectra.matpower_import_auto_profile(mpc, cfg; mode = :apply)
    @test applied.config.matpower.shift_unit === :deg
    applied_row = only(row for row in applied.rows if row.option == "matpower_import.shift_unit")
    @test applied_row.action === :applied

    ambiguous_cfg = Sparlectra.SparlectraConfig(Dict("matpower_import" => Dict("auto_profile" => "apply")))
    ambiguous = Sparlectra.matpower_import_auto_profile(_auto_profile_pv_mismatch_case(), ambiguous_cfg; mode = :apply)
    compare_row = only(row for row in ambiguous.rows if row.option == "matpower_import.compare_voltage_reference")
    @test compare_row.recommended == "hybrid"
    @test compare_row.action === :recommend
    @test ambiguous.config.matpower.compare_voltage_reference === :imported_setpoint

    io = IOBuffer()
    Sparlectra.print_matpower_import_auto_profile(io, ambiguous.rows)
    Sparlectra.print_matpower_import_auto_profile_effective_options(io, ambiguous.config)
    text = String(take!(io))
    @test occursin("MATPOWER auto-profile recommendations", text)
    @test occursin("Final effective MATPOWER auto-profile options", text)
    @test occursin("matpower_import.compare_voltage_reference", text)
  end

  @testset "Configuration value-domain validation" begin
    tol_bad = tempname() * ".yaml"
    write(tol_bad, "power_flow:\n  tol: 0\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(tol_bad; reload = true)

    damping_bad = tempname() * ".yaml"
    write(damping_bad, "power_flow:\n  autodamp_min: 0\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(damping_bad; reload = true)

    angle_limit_bad = tempname() * ".yaml"
    write(angle_limit_bad, "power_flow:\n  start_mode:\n    dc_angle_limit_deg: 0\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(angle_limit_bad; reload = true)

    bool_ok_true = tempname() * ".yaml"
    write(bool_ok_true, "power_flow:\n  autodamp: true\n")
    @test Sparlectra.load_sparlectra_config(bool_ok_true; reload = true).powerflow.autodamp === true

    bool_ok_false = tempname() * ".yaml"
    write(bool_ok_false, "power_flow:\n  autodamp: false\n")
    @test Sparlectra.load_sparlectra_config(bool_ok_false; reload = true).powerflow.autodamp === false

    enum_bad = tempname() * ".yaml"
    write(enum_bad, "power_flow:\n  start_mode:\n    angle_mode: invalid_angle\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(enum_bad; reload = true)
    wrong_branch_bad = tempname() * ".yaml"
    write(wrong_branch_bad, "power_flow:\n  wrong_branch_detection: maybe\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(wrong_branch_bad; reload = true)
  end

  @testset "Removed diagnostics keys are rejected" begin
    removed_diag_keys = (
      "matpower_reference",
      "branch_shift_conventions",
      "negative_branch_impedance",
      "pv_voltage_references",
      "residual_clusters",
      "nodal_balance_breakdown",
      "branch_neighborhood",
      "detailed_log",
    )
    for key in removed_diag_keys
      cfg_bad = tempname() * ".yaml"
      write(cfg_bad, "diagnostics:\n  $(key): true\n")
      @test_throws ArgumentError Sparlectra.load_sparlectra_config(cfg_bad; reload = true)
    end
  end

end

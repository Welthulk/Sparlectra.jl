# Copyright 2023-2026 Udo Schmitz
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

# file: test/test_configuration_coverage.jl
# purpose: tests configuration YAML key coverage, sentinel forwarding,
#          value-domain validation, deprecated and removed key handling,
#          MATPOWER auto-profile decision rules, and console_live capture
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
      "power_flow.method", "power_flow.solver", "power_flow.linear_solver", "power_flow.apslf.order", "power_flow.apslf.use_pade", "power_flow.apslf.nr_polish", "power_flow.apslf_start.enabled", "power_flow.apslf_start.order", "power_flow.flatstart", "power_flow.tol", "power_flow.max_iter", "power_flow.autodamp", "power_flow.autodamp_min", "power_flow.auto_slack", "power_flow.rescue", "power_flow.dc.fallback", "power_flow.wrong_branch_detection", "power_flow.wrong_branch_rescue", "power_flow.wrong_branch_min_vm_pu", "power_flow.wrong_branch_max_vm_pu", "power_flow.wrong_branch_max_angle_spread_deg", "power_flow.wrong_branch_max_branch_angle_deg", "power_flow.wrong_branch_min_low_vm_count", "power_flow.wrong_branch_rescue_max_attempts", "power_flow.rectangular_workspace_reuse", "power_flow.rectangular_preallocate_workspace", "power_flow.rectangular_workspace_min_buses",
      "power_flow.islands.enabled", "power_flow.islands.mode", "power_flow.islands.reference_policy", "power_flow.islands.diagnostic_continue_after_failure",
      "power_flow.distributed_slack.enabled", "power_flow.distributed_slack.p_mode", "power_flow.distributed_slack.respect_p_limits", "power_flow.distributed_slack.fallback", "power_flow.distributed_slack.weights",
      "power_flow.external_grid.enabled", "power_flow.external_grid.source", "power_flow.external_grid.sk_MVA", "power_flow.external_grid.rx",
      "power_flow.start_mode.angle_mode", "power_flow.start_mode.voltage_mode", "power_flow.start_mode.profile_source", "power_flow.start_mode.start_projection", "power_flow.start_mode.try_dc_start", "power_flow.start_mode.try_blend_scan", "power_flow.start_mode.branch_guard", "power_flow.start_mode.measure_candidates", "power_flow.start_mode.accept_unmeasured_dc_start", "power_flow.start_mode.dc_seed_unconditional", "power_flow.start_mode.reuse_import_data", "power_flow.start_mode.blend_lambdas", "power_flow.start_mode.dc_angle_limit_deg",
      "power_flow.start_current_iteration.enabled", "power_flow.start_current_iteration.max_iter", "power_flow.start_current_iteration.tol", "power_flow.start_current_iteration.damping", "power_flow.start_current_iteration.accept_only_if_improved", "power_flow.start_current_iteration.min_improvement_factor", "power_flow.start_current_iteration.vm_min_pu", "power_flow.start_current_iteration.vm_max_pu", "power_flow.start_current_iteration.max_angle_step_deg", "power_flow.start_current_iteration.only_for_large_cases",
      "power_flow.merit.enabled", "power_flow.merit.armijo_c1", "power_flow.merit.scale_p", "power_flow.merit.scale_q", "power_flow.merit.scale_v", "power_flow.merit.fallback_max_mismatch",
      "power_flow.trust_region.enabled", "power_flow.trust_region.initial_radius", "power_flow.trust_region.min_radius", "power_flow.trust_region.max_radius", "power_flow.trust_region.eta_accept", "power_flow.trust_region.shrink_factor", "power_flow.trust_region.expand_factor", "power_flow.trust_region.expand_threshold", "power_flow.trust_region.step_mode",
      "power_flow.qlimits.enabled", "power_flow.qlimits.enforcement_mode", "power_flow.qlimits.start_iter", "power_flow.qlimits.start_mode", "power_flow.qlimits.auto_q_delta_pu", "power_flow.qlimits.hysteresis_pu", "power_flow.qlimits.cooldown_iters", "power_flow.qlimits.trace_buses", "power_flow.qlimits.lock_pv_to_pq_buses",
      "power_flow.qlimits.guard.enabled", "power_flow.qlimits.guard.min_q_range_pu", "power_flow.qlimits.guard.narrow_range_mode", "power_flow.qlimits.guard.zero_range_mode", "power_flow.qlimits.guard.violation_mode", "power_flow.qlimits.guard.violation_threshold_pu", "power_flow.qlimits.guard.max_switches", "power_flow.qlimits.guard.max_remaining_violations", "power_flow.qlimits.guard.accept_bounded_violations", "power_flow.qlimits.guard.freeze_after_repeated_switching", "power_flow.qlimits.guard.log",
      "state_estimation.enabled", "state_estimation.method", "state_estimation.tol", "state_estimation.max_iter", "state_estimation.flatstart", "state_estimation.jac_eps", "state_estimation.update_net", "state_estimation.pmu_ref_offset", "state_estimation.observability.enabled",
      "matpower_import.case", "matpower_import.cases", "matpower_import.auto_profile", "matpower_import.auto_profile_log", "matpower_import.pv_voltage_source", "matpower_import.pv_voltage_mismatch_tol_pu", "matpower_import.compare_voltage_reference", "matpower_import.bus_shunt_model", "matpower_import.shift_unit", "matpower_import.shift_sign", "matpower_import.ratio", "matpower_import.enable_pq_gen_controllers", "matpower_import.preallocate_network", "matpower_import.preallocate_min_buses", "matpower_import.apply_bus_names", "matpower_import.apply_branch_names", "matpower_import.apply_branch_kind", "matpower_import.import_for001_contingencies", "matpower_import.matpower_dcline_mode", "matpower_import.net_cache.enabled",
      "cgmes_import.path", "cgmes_import.base_mva", "cgmes_import.require_boundary", "cgmes_import.tap_control", "cgmes_import.machine_control", "cgmes_import.ignore_connected", "cgmes_import.vset_min_pu", "cgmes_import.vset_max_pu", "cgmes_import.multi_slack", "cgmes_import.start_values", "cgmes_import.placeholder_guards", "cgmes_import.infer_base_voltages", "cgmes_import.hvdc_mode",
      "short_circuit.c_factor",
      "transformer.tap_changer_model",
      "matpower_export.write_solution",
      "performance.enabled", "performance.level", "performance.print_to_console", "performance.write_to_logfile", "performance.show_allocations", "performance.show_iteration_table", "performance.compact_logging", "performance.representative_warmup_runs", "performance.compare_cold_warm", "performance.skip_reference_comparison", "performance.skip_expensive_diagnostics", "performance.skip_branch_neighborhood_report", "performance.max_diagnostic_rows",
      "runtime.print_thread_config", "runtime.julia_threads", "runtime.blas_threads", "runtime.casefile", "runtime.case_name", "runtime.case_source", "runtime.configured_default_casefile",
      "runtime.parallel.enabled", "runtime.parallel.max_tasks", "runtime.parallel.min_work_items",
      "diagnostics.log_effective_config",
      "output.console_summary", "output.console_live", "output.console_auto_profile", "output.console_diagnostics", "output.console_q_limit_events", "output.console_max_rows", "output.logfile_results", "output.result_table_max_rows", "output.result_table_large_case_threshold_buses", "output.result_table_large_case_mode", "output.detailed_result_csv_write_mode", "output.detailed_result_csv_exporter", "output.detailed_result_csv_direct_threshold_buses", "output.detailed_result_csv_buffer_initial_bytes", "output.detailed_result_csv_buffer_max_bytes", "output.detailed_result_csv_streaming_threshold_rows", "output.logfile_diagnostics", "output.logfile_performance", "output.logfile_warnings",
      "benchmark.enabled", "benchmark.methods", "benchmark.seconds", "benchmark.samples", "benchmark.show_once", "benchmark.show_once_output", "benchmark.show_once_max_nodes",
      "contingency.rescue_ladder",
      "control.enabled", "control.max_outer_iterations", "control.trace", "control.log_iterations", "control.stop_on_pf_failure", "control.controllers",
      "webui.show_case_settings_notice", "webui.warmup",
    ])
    reserved_keys = Set(["extensions.reserved"])
    mapped_or_reserved = union(mapped_keys, reserved_keys)

    @test isempty(setdiff(leaves, mapped_or_reserved))
    @test isempty(setdiff(reserved_keys, leaves))

    expected_consumers = Dict(
      "power_flow.tol" => :PowerFlowConfig,
      "power_flow.solver" => :PowerFlowConfig,
      "power_flow.apslf.order" => :ApslfConfig,
      "power_flow.apslf_start.enabled" => :ApslfStartConfig,
      "power_flow.max_iter" => :PowerFlowConfig,
      "power_flow.autodamp" => :PowerFlowConfig,
      "power_flow.autodamp_min" => :PowerFlowConfig,
      "power_flow.wrong_branch_detection" => :PowerFlowConfig,
      "power_flow.start_mode.angle_mode" => :StartModeConfig,
      "power_flow.start_mode.voltage_mode" => :StartModeConfig,
      "power_flow.start_mode.profile_source" => :StartModeConfig,
      "power_flow.qlimits.enforcement_mode" => :QLimitConfig,
      "power_flow.qlimits.start_iter" => :QLimitConfig,
      "power_flow.merit.enabled" => :MeritLineSearchConfig,
      "power_flow.merit.armijo_c1" => :MeritLineSearchConfig,
      "power_flow.trust_region.enabled" => :TrustRegionConfig,
      "power_flow.trust_region.initial_radius" => :TrustRegionConfig,
      "matpower_import.shift_sign" => :MatpowerImportConfig,
      "matpower_import.cases" => :MatpowerImportConfig,
      "transformer.tap_changer_model" => :TransformerConfig,
      "matpower_export.write_solution" => :MatpowerExportConfig,
      "output.logfile_performance" => :OutputConfig,
      "benchmark.enabled" => :BenchmarkConfig,
      "state_estimation.method" => :StateEstimationConfig,
      "runtime.julia_threads" => :RuntimeConfig,
      "diagnostics.log_effective_config" => :DiagnosticsConfig,
      "extensions.reserved" => :Reserved,
      "webui.show_case_settings_notice" => :WebUIConfig,
    )
    @test all(haskey(expected_consumers, key) for key in keys(expected_consumers))
    @test expected_consumers["extensions.reserved"] === :Reserved
  end

  @testset "Test runner output mode helpers" begin
    @test selected_test_profile(String[], Dict("SPARLECTRA_TEST_PROFILE" => "extended")) === :extended
    @test selected_test_profile(["fast"], Dict("SPARLECTRA_TEST_PROFILE" => "extended")) === :fast
    @test selected_test_profile(["--verbose", "all"], Dict{String,String}()) === :all
    @test sparlectra_test_verbose(["--verbose"], Dict{String,String}())
    @test sparlectra_test_verbose(String[], Dict("SPARLECTRA_TEST_VERBOSE" => "1"))
    @test !sparlectra_test_verbose(String[], Dict{String,String}())

    quiet_path = tempname()
    open(quiet_path, "w+") do io
      redirect_stdio(stdout = io) do
        @test quiet_test_output(verbose = false) do
          println("Runtime casefile: hidden")
          return :quiet_result
        end === :quiet_result
      end
      seekstart(io)
      @test isempty(read(io, String))
    end

    verbose_path = tempname()
    open(verbose_path, "w+") do io
      redirect_stdio(stdout = io) do
        @test quiet_test_output(verbose = true) do
          println("Runtime casefile: visible")
          return :verbose_result
        end === :verbose_result
      end
      seekstart(io)
      @test occursin("Runtime casefile: visible", read(io, String))
    end
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
matpower_export:
  write_solution: false
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
    @test cfg.matpower_export.write_solution === false
  end

  @testset "Q-limit enforcement mode user YAML keys" begin
    for mode in (:active_set, :classic_simultaneous, :classic_one_at_a_time)
      cfgfile = tempname() * ".yaml"
      write(cfgfile, """
power_flow:
  qlimits:
    enforcement_mode: $(mode)
""")
      cfg = Sparlectra.load_sparlectra_config(cfgfile; reload = true)
      @test cfg.powerflow.qlimits.enforcement_mode === mode
    end
    for (legacy, canonical) in ((:matpower_simultaneous, :classic_simultaneous), (:matpower_one_at_a_time, :classic_one_at_a_time))
      cfgfile = tempname() * ".yaml"
      write(cfgfile, """
power_flow:
  qlimits:
    enforcement_mode: $(legacy)
""")
      cfg = Sparlectra.load_sparlectra_config(cfgfile; reload = true)
      @test cfg.powerflow.qlimits.enforcement_mode === canonical
    end
    err = try
      Sparlectra.QLimitConfig(Dict("enforcement_mode" => "definitely_not_supported"))
      nothing
    catch caught
      caught
    end
    @test err isa ArgumentError
    @test occursin("classic_simultaneous", sprint(showerror, err))
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

    apply_cfg = Sparlectra.SparlectraConfig(Dict(
      "matpower_import" => Dict("auto_profile" => "apply", "shift_unit" => "rad"),
    ))
    applied = Sparlectra.matpower_import_auto_profile(mpc, apply_cfg; mode = :apply)
    @test applied.config.matpower.shift_unit === :deg
    applied_row = only(row for row in applied.rows if row.option == "matpower_import.shift_unit")
    @test applied_row.action === :applied
    @test any(pair -> first(pair) === :shift_unit && last(pair) === :deg, applied.applied)
    convention_rows = Sparlectra._matpower_import_auto_profile_convention_scan(mpc)
    @test length(convention_rows) == 8

    ambiguous_cfg = Sparlectra.SparlectraConfig(Dict("matpower_import" => Dict("auto_profile" => "apply")))
    ambiguous = Sparlectra.matpower_import_auto_profile(_auto_profile_pv_mismatch_case(), ambiguous_cfg; mode = :apply)
    compare_row = only(row for row in ambiguous.rows if row.option == "matpower_import.compare_voltage_reference")
    @test compare_row.recommended == "hybrid"
    @test compare_row.action === :applied
    @test ambiguous.config.matpower.compare_voltage_reference === :hybrid
    @test any(pair -> first(pair) === :compare_voltage_reference && last(pair) === :hybrid, ambiguous.applied)

    io = IOBuffer()
    Sparlectra.print_matpower_import_auto_profile(io, ambiguous.rows)
    Sparlectra.print_matpower_import_auto_profile_effective_options(io, ambiguous.config)
    text = String(take!(io))
    @test occursin("MATPOWER auto-profile recommendations", text)
    @test occursin("Final effective MATPOWER auto-profile options", text)
    @test occursin("matpower_import.compare_voltage_reference", text)

    explicit_start_cfg = Sparlectra.SparlectraConfig(Dict(
      "matpower_import" => Dict("auto_profile" => "apply"),
      "power_flow" => Dict(
        "start_mode" => Dict(
          "angle_mode" => "classic",
          "voltage_mode" => "classic",
          "profile_source" => "flat",
        ),
        "qlimits" => Dict("start_mode" => "iteration"),
      ),
    ))
    fragile_scan = mpc -> [
      (shift_unit = :deg, shift_sign = 1.0, ratio = :normal, stats = (; ok = true), score = 0.01),
      (shift_unit = :deg, shift_sign = -1.0, ratio = :normal, stats = (; ok = true), score = 0.01),
      (shift_unit = :deg, shift_sign = 1.0, ratio = :reciprocal, stats = (; ok = true), score = 0.01),
      (shift_unit = :deg, shift_sign = -1.0, ratio = :reciprocal, stats = (; ok = true), score = 0.01),
      (shift_unit = :rad, shift_sign = 1.0, ratio = :normal, stats = (; ok = true), score = 0.01),
      (shift_unit = :rad, shift_sign = -1.0, ratio = :normal, stats = (; ok = true), score = 0.01),
      (shift_unit = :rad, shift_sign = 1.0, ratio = :reciprocal, stats = (; ok = true), score = 0.01),
      (shift_unit = :rad, shift_sign = -1.0, ratio = :reciprocal, stats = (; ok = true), score = 0.01),
    ]
    conservative_mpc = _auto_profile_pv_mismatch_case()
    conservative_mpc.gen[:, 4] .= 1.0
    conservative_mpc.gen[:, 5] .= 0.0
    conservative = Sparlectra.matpower_import_auto_profile(conservative_mpc, explicit_start_cfg; mode = :apply, convention_scan = fragile_scan)
    @test conservative.config.powerflow.start_mode.angle_mode === :classic
    @test conservative.config.powerflow.start_mode.voltage_mode === :classic
    @test conservative.config.powerflow.start_mode.profile_source === :flat
    @test conservative.config.powerflow.qlimits.start_mode === :iteration
    start_angle_row = only(row for row in conservative.rows if row.option == "power_flow.start_mode.angle_mode")
    qlimit_row = only(row for row in conservative.rows if row.option == "power_flow.qlimits.start_mode")
    @test start_angle_row.recommended == "dc"
    @test start_angle_row.action === :skipped
    @test qlimit_row.action === :skipped
    io = IOBuffer()
    Sparlectra.write_matpower_import_auto_profile(io, conservative, conservative.config; casefile = "synthetic_fragile.m")
    conservative_text = String(take!(io))
    @test occursin("power_flow.start_mode.angle_mode", conservative_text)
    @test occursin("skipped", conservative_text)
    @test occursin("power_flow.start_mode.angle_mode: classic", conservative_text)

    io = IOBuffer()
    Sparlectra.write_matpower_import_auto_profile(io, applied, applied.config; casefile = "case13659pegase.m")
    apply_text = String(take!(io))
    @test occursin("Runtime casefile: case13659pegase.m", apply_text)
    @test occursin("Original MATPOWER import options:", apply_text)
    @test occursin("Auto-profile recommendation:", apply_text)
    @test occursin("Final effective MATPOWER import options:", apply_text)
    @test occursin("auto_profile = apply", apply_text)
    @test occursin("shift_unit   = deg", apply_text)

    oom_cfg = Sparlectra.SparlectraConfig(Dict(
      "matpower_import" => Dict("auto_profile" => "apply", "shift_unit" => "rad"),
    ))
    oom_result = @test_logs (:warn, r"matpower_auto_profile_scan_skipped") Sparlectra.matpower_import_auto_profile(
      mpc,
      oom_cfg;
      mode = :apply,
      convention_scan = mpc -> throw(OutOfMemoryError()),
    )
    @test oom_result.config.matpower.shift_unit === :rad
    @test isempty(oom_result.applied)
    @test any(row -> occursin("matpower_auto_profile_scan_skipped", row.reason), oom_result.rows)
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

    # _copy_sparlectra_with_powerflow must carry EVERY field (regression: the
    # former hand-written keyword list silently reset cgmes and webui to
    # defaults on every copy).
    cgmes_cfg = Sparlectra.SparlectraConfig(cgmes = Sparlectra.CGMESImportConfig(start_values = :sv, base_mva = 50.0))
    cgmes_copied = Sparlectra._copy_sparlectra_with_powerflow(cgmes_cfg, cgmes_cfg.powerflow)
    @test cgmes_copied.cgmes.start_values === :sv
    @test cgmes_copied.cgmes.base_mva == 50.0
    @test all(getfield(cgmes_copied, f) == getfield(cgmes_cfg, f) for f in fieldnames(Sparlectra.SparlectraConfig))

    # cgmes_import.start_values: default flat, sv accepted, invalid rejected
    # with the key name in the message.
    # auto is the shipped default: deliveries with an SvVoltage state are
    # solved from it, others fall back to the flat start (resolved per run).
    @test Sparlectra.CGMESImportConfig().start_values === :auto
    flat_forced = tempname() * ".yaml"
    write(flat_forced, "cgmes_import:\n  start_values: flat\n")
    @test Sparlectra.load_sparlectra_config(flat_forced; reload = true).cgmes.start_values === :flat
    sv_ok = tempname() * ".yaml"
    write(sv_ok, "cgmes_import:\n  start_values: sv\n")
    @test Sparlectra.load_sparlectra_config(sv_ok; reload = true).cgmes.start_values === :sv
    sv_bad = tempname() * ".yaml"
    write(sv_bad, "cgmes_import:\n  start_values: warm\n")
    err = try
      Sparlectra.load_sparlectra_config(sv_bad; reload = true)
      nothing
    catch e
      e
    end
    @test err isa ArgumentError
    @test occursin("cgmes_import.start_values", sprint(showerror, err))

    # matpower_import.matpower_dcline_mode: the stale early-template value
    # reject_active must not brick old user configs — it loads as
    # pf_injections with a deprecation warning; ignore_inactive stays a
    # deliberate choice and loads verbatim.
    dcline_stale = tempname() * ".yaml"
    write(dcline_stale, "matpower_import:\n  matpower_dcline_mode: reject_active\n")
    stale_cfg = @test_logs (:warn, r"reject_active is deprecated in configuration files") Sparlectra.load_sparlectra_config(dcline_stale; reload = true)
    @test stale_cfg.matpower.matpower_dcline_mode === :pf_injections
    dcline_deliberate = tempname() * ".yaml"
    write(dcline_deliberate, "matpower_import:\n  matpower_dcline_mode: ignore_inactive\n")
    @test Sparlectra.load_sparlectra_config(dcline_deliberate; reload = true).matpower.matpower_dcline_mode === :ignore_inactive

    # cgmes_import.placeholder_guards: default warn_skip, strict accepted,
    # invalid rejected with the key name in the message.
    @test Sparlectra.CGMESImportConfig().placeholder_guards === :warn_skip
    pg_ok = tempname() * ".yaml"
    write(pg_ok, "cgmes_import:\n  placeholder_guards: strict\n")
    @test Sparlectra.load_sparlectra_config(pg_ok; reload = true).cgmes.placeholder_guards === :strict
    pg_bad = tempname() * ".yaml"
    write(pg_bad, "cgmes_import:\n  placeholder_guards: ignore\n")
    pg_err = try
      Sparlectra.load_sparlectra_config(pg_bad; reload = true)
      nothing
    catch e
      e
    end
    @test pg_err isa ArgumentError
    @test occursin("cgmes_import.placeholder_guards", sprint(showerror, pg_err))
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

  @testset "Q-limit start mode public values" begin
    for mode in ("iteration", "auto", "iteration_or_auto")
      cfg = Sparlectra.SparlectraConfig(Dict(
        "power_flow" => Dict("qlimits" => Dict("start_mode" => mode)),
      ))
      @test cfg.powerflow.qlimits.start_mode === Symbol(mode)
    end
  end

  @testset "Configuration refresh" begin
    stale = tempname() * ".yaml"
    write(stale, "power_flow:\n  tol: 1.0e-6\n  start_mode:\n    voltage_mode: bus_vm_va_blend\n  qlimits:\n    enabled: true\n")
    dry = Sparlectra.refresh_sparlectra_config_file(stale)
    @test dry.success
    @test dry.changed
    @test !dry.written
    @test "power_flow.qlimits.enforcement_mode" in dry.missing_keys
    @test "power_flow.start_mode.voltage_mode" in dry.normalized_keys
    @test occursin("tol: 1.0e-6", dry.refreshed_text)
    @test occursin("voltage_mode: profile_blend", dry.refreshed_text)
    @test occursin("profile_source: matpower_reference", dry.refreshed_text)

    written = Sparlectra.refresh_sparlectra_config_file(stale; write = true)
    @test written.success
    @test written.written
    @test written.backup_path !== nothing
    @test isfile(written.backup_path)
    cfg = Sparlectra.load_sparlectra_config(stale; reload = true)
    @test cfg.powerflow.tol == 1.0e-6
    @test cfg.powerflow.start_mode.voltage_mode === :profile_blend
    @test cfg.powerflow.start_mode.profile_source === :matpower_reference
    @test cfg.powerflow.qlimits.enforcement_mode === :active_set

    for (legacy, canonical) in (("matpower_simultaneous", "classic_simultaneous"), ("matpower_one_at_a_time", "classic_one_at_a_time"))
      p = tempname() * ".yaml"
      write(p, "power_flow:\n  qlimits:\n    enforcement_mode: $(legacy)\n")
      result = Sparlectra.refresh_sparlectra_config_file(p)
      @test "power_flow.qlimits.enforcement_mode" in result.normalized_keys
      @test occursin("enforcement_mode: $(canonical)", result.refreshed_text)
    end

    dup = tempname() * ".yaml"
    write(dup, "output:\n  detailed_result_csv_exporter: auto\n  detailed_result_csv_exporter: direct\n")
    dup_result = Sparlectra.refresh_sparlectra_config_file(dup; write = true)
    @test !dup_result.success
    @test !dup_result.written
    @test "output.detailed_result_csv_exporter" in dup_result.duplicate_keys
    @test occursin("direct", read(dup, String))
  end

  @testset "console_live capture tees to console and archive identically" begin
    # live=false: output only in the archive stream. The archive must be a
    # real OS stream here (redirect_stdout rejects IOBuffer) — exactly what
    # production passes (the open run.log IOStream).
    quiet_path, quiet_io = mktemp()
    result = redirect_stdout(devnull) do
      Sparlectra._capture_run_output(quiet_io) do
        println("captured line")
        42
      end
    end
    close(quiet_io)
    @test result == 42
    @test occursin("captured line", read(quiet_path, String))

    # live=true: the same bytes land in the archive AND on the (outer) console.
    # The read must happen AFTER the redirect block: redirect_stdout dups the
    # pipe's write end onto fd 1, so closing outer.in inside the block leaves
    # the dup open and read() would never see EOF (deadlock).
    archive = IOBuffer()
    outer = Pipe()
    Base.link_pipe!(outer; reader_supports_async = true, writer_supports_async = true)
    redirect_stdout(outer) do
      Sparlectra._capture_run_output(archive; live = true) do
        println("teed line")
      end
    end
    close(outer.in)
    console_text = read(outer, String)
    @test occursin("teed line", String(take!(archive)))
    @test occursin("teed line", console_text)

    # error path: output written BEFORE a mid-run throw must be flushed into
    # the archive (pump drained in the finally block), and the exception
    # must propagate — exactly what a crash post-mortem needs from run.log
    archive2 = IOBuffer()
    thrown = redirect_stdout(devnull) do
      try
        Sparlectra._capture_run_output(archive2; live = true) do
          println("before crash")
          error("boom")
        end
        false
      catch err
        occursin("boom", sprint(showerror, err))
      end
    end
    @test thrown
    @test occursin("before crash", String(take!(archive2)))

    # config surface: default off, parseable on
    @test !Sparlectra.OutputConfig().console_live
    @test Sparlectra.OutputConfig(Dict("output" => Dict("console_live" => true))).console_live
  end

  @testset "Deprecated diagnostics.* keys load with a warning, not an error" begin
    # Regression (2026-07-30): stored user/webui configs still carry the old
    # diagnostics.console_* duplicates of output.*; after their removal from
    # the default file the unknown-key validation rejected every such config
    # ("Unknown Sparlectra configuration key: diagnostics.console_diagnostics")
    # and bricked the Web UI start. Deprecated keys must warn and be ignored.
    p = tempname() * ".yaml"
    write(p, "diagnostics:\n  console_diagnostics: full\n  console_max_rows: 50\n  log_effective_config: true\n")
    cfg = @test_logs (:warn, r"diagnostics\.console_diagnostics is deprecated") (:warn, r"diagnostics\.console_max_rows is deprecated") match_mode = :any Sparlectra.load_sparlectra_config(p; reload = true)
    @test cfg isa Sparlectra.SparlectraConfig
    @test cfg.diagnostics.log_effective_config
    # genuinely unknown keys still fail loudly
    bad = tempname() * ".yaml"
    write(bad, "diagnostics:\n  no_such_key: 1\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(bad; reload = true)
    # the config-refresh path migrates stored files by dropping the dead keys
    # (scoped to the diagnostics block — output.console_diagnostics is the
    # legitimate owner and stays in the refreshed text)
    refreshed = Sparlectra.refresh_sparlectra_config_file(p)
    @test "diagnostics.console_diagnostics" in refreshed.normalized_keys
    @test "diagnostics.console_max_rows" in refreshed.normalized_keys
    diag_block = match(r"(?m)^diagnostics:\n((?:^  .*\n?)*)", refreshed.refreshed_text)
    @test diag_block !== nothing
    @test !occursin("console_diagnostics", diag_block.captures[1])
    @test occursin("log_effective_config", diag_block.captures[1])
  end
end

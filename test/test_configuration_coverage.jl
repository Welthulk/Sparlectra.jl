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
  yaml = Sparlectra.load_yaml_dict(joinpath(@__DIR__, "..", "src", "config", "configuration.yaml.example"))
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
      "power_flow.method", "power_flow.mode", "power_flow.solver", "power_flow.linear_solver", "power_flow.apslf.order", "power_flow.apslf.use_pade", "power_flow.apslf.nr_polish", "power_flow.apslf_start.enabled", "power_flow.apslf_start.order", "power_flow.flatstart", "power_flow.tol", "power_flow.max_iter", "power_flow.autodamp", "power_flow.autodamp_min", "power_flow.auto_slack", "power_flow.rescue", "power_flow.dc.fallback", "power_flow.wrong_branch_detection", "power_flow.wrong_branch_rescue", "power_flow.wrong_branch_min_vm_pu", "power_flow.wrong_branch_max_vm_pu", "power_flow.wrong_branch_max_angle_spread_deg", "power_flow.wrong_branch_max_branch_angle_deg", "power_flow.wrong_branch_min_low_vm_count", "power_flow.wrong_branch_rescue_max_attempts", "power_flow.rectangular_workspace_reuse", "power_flow.rectangular_preallocate_workspace", "power_flow.rectangular_workspace_min_buses",
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
      "matpower_import.pv_voltage_source", "matpower_import.pv_voltage_mismatch_tol_pu", "matpower_import.compare_voltage_reference", "matpower_import.shift_unit", "matpower_import.shift_sign", "matpower_import.ratio", "matpower_import.enable_pq_gen_controllers", "matpower_import.apply_bus_names", "matpower_import.apply_branch_names", "matpower_import.apply_branch_kind", "matpower_import.import_for001_contingencies", "matpower_import.matpower_dcline_mode",
      "model.bus_shunt_model", "model.tap_changer_model", "model.auto_profile", "model.auto_profile_log", "model.net_cache_enabled", "model.preallocate_network", "model.preallocate_min_buses",
      "cgmes_import.path", "cgmes_import.base_mva", "cgmes_import.require_boundary", "cgmes_import.tap_control", "cgmes_import.machine_control", "cgmes_import.ignore_connected", "cgmes_import.vset_min_pu", "cgmes_import.vset_max_pu", "cgmes_import.multi_slack", "cgmes_import.start_values", "cgmes_import.placeholder_guards", "cgmes_import.infer_base_voltages", "cgmes_import.hvdc_mode",
      "short_circuit.c_factor", "short_circuit.sweep_method", "short_circuit.takahashi_min_buses",
      "matpower_export.write_solution",
      "performance.enabled", "performance.level", "performance.print_to_console", "performance.write_to_logfile", "performance.show_allocations", "performance.show_iteration_table", "performance.compact_logging", "performance.representative_warmup_runs", "performance.compare_cold_warm", "performance.skip_reference_comparison", "performance.skip_expensive_diagnostics", "performance.skip_branch_neighborhood_report", "performance.max_diagnostic_rows",
      "config_version", "scope",
      "runtime.case", "runtime.cases", "runtime.print_thread_config", "runtime.julia_threads", "runtime.blas_threads", "runtime.casefile", "runtime.case_name", "runtime.case_source", "runtime.configured_default_casefile",
      "runtime.parallel.enabled", "runtime.parallel.max_tasks", "runtime.parallel.min_work_items",
      "diagnostics.log_effective_config",
      "output.console_summary", "output.console_live", "output.console_auto_profile", "output.console_diagnostics", "output.console_q_limit_events", "output.console_max_rows", "output.logfile_results", "output.result_table_max_rows", "output.result_table_large_case_threshold_buses", "output.result_table_large_case_mode", "output.detailed_result_csv_write_mode", "output.detailed_result_csv_exporter", "output.detailed_result_csv_direct_threshold_buses", "output.detailed_result_csv_buffer_initial_bytes", "output.detailed_result_csv_buffer_max_bytes", "output.detailed_result_csv_streaming_threshold_rows", "output.logfile_diagnostics", "output.logfile_performance", "output.logfile_warnings", "output.startup_latency_hint",
      "benchmark.enabled", "benchmark.methods", "benchmark.seconds", "benchmark.samples", "benchmark.show_once", "benchmark.show_once_output", "benchmark.show_once_max_nodes",
      "contingency.rescue_ladder", "contingency.screening.mode", "contingency.screening.margin_pct",
      "control.enabled", "control.max_outer_iterations", "control.trace", "control.log_iterations", "control.stop_on_pf_failure", "control.controllers",
      "webui.show_case_settings_notice", "webui.operation_log_retention_days",
      # task_config_arrival_v0100: these were in the typed configuration and
      # documented, but missing from the template, which made them
      # unreachable ("Unknown Sparlectra configuration key") for every user
      # who followed the documentation
      "power_flow.tol_MW", "power_flow.dc.angle_reference_deg", "power_flow.dc.ignore_out_of_service",
      "state_estimation.imag_activation_iteration", "state_estimation.ia_current_floor_A",
      "state_estimation.update_shunts", "state_estimation.update_taps",
      "state_estimation.report_residual_correlation", "state_estimation.k_eliminate",
      "state_estimation.robust_mode", "state_estimation.robust", "state_estimation.robust_start_iteration",
      "state_estimation.robust_k1", "state_estimation.robust_k2", "state_estimation.k_suppress",
      "state_estimation.suppression_sigma", "state_estimation.max_eliminations",
      "state_estimation.rank_tol_factor",
      "state_estimation.takahashi_min_states", "state_estimation.topology_precheck",
      "state_estimation.topology_open_flow_k", "state_estimation.topology_dead_flow_k",
      "state_estimation.topology_voltage_k", "state_estimation.topology_kcl_k",
      "state_estimation.topology_cluster_min",
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
      "runtime.cases" => :RuntimeConfig,
      "model.tap_changer_model" => :ModelConfig,
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

  @testset "configuration version, scope, and case precedence" begin
    dir = mktempdir()
    # a file without config_version reads as version 0, with one warning and
    # the 0 => 1 aliases applied
    v0 = joinpath(dir, "v0.yaml")
    write(v0, "matpower_import:\n  bus_shunt_model: voltage_dependent_injection\n  case: case57.m\ntransformer:\n  tap_changer_model: impedance_correction\n")
    cfg0 = Sparlectra.load_sparlectra_config(v0; reload = true)
    @test cfg0.model.bus_shunt_model === :voltage_dependent_injection
    @test cfg0.model.tap_changer_model === :impedance_correction
    @test cfg0.runtime.case == "case57.m"
    # ONE alias warning per file, not one per key (2026-09-07): a user file
    # with six legacy names produced six boxed warnings at every start, which
    # is what a Windows start reported as "still all those warnings". The
    # information that matters is WHICH names are still in use, and that fits
    # in a single line. TestLogger, not the message count on the console:
    # maxlog suppression is a property of the logger, the emission is not.
    v0_many = joinpath(dir, "v0_many.yaml")
    write(v0_many, "matpower_import:\n  auto_profile: true\n  auto_profile_log: false\n  preallocate_network: true\n  preallocate_min_buses: 500\n  case: case57.m\ntransformer:\n  tap_changer_model: impedance_correction\n")
    logger = Test.TestLogger(min_level = Logging.Warn)
    Logging.with_logger(logger) do
      Sparlectra.load_sparlectra_config(v0_many; reload = true)
    end
    alias_records = [r for r in logger.logs if occursin("legacy key name", r.message)]
    @test length(alias_records) == 1
    # and it names every one of them, so nothing is lost by collapsing
    for key in ("matpower_import.auto_profile", "matpower_import.auto_profile_log",
                "matpower_import.preallocate_network", "matpower_import.preallocate_min_buses",
                "matpower_import.case", "transformer.tap_changer_model")
      @test (key, occursin(key, alias_records[1].message)) == (key, true)
    end
    # after the documented one-time refresh the file loads without any of it
    Sparlectra.refresh_sparlectra_config_file(v0_many; write = true, backup = false)
    quiet = Test.TestLogger(min_level = Logging.Warn)
    migrated = Logging.with_logger(quiet) do
      Sparlectra.load_sparlectra_config(v0_many; reload = true)
    end
    @test isempty([r for r in quiet.logs if occursin("legacy key name", r.message)])
    @test migrated.model.tap_changer_model === :impedance_correction
    @test migrated.runtime.case == "case57.m"
    @test migrated.model.preallocate_min_buses == 500
    # a version newer than the running Sparlectra is an error, not a guess
    v9 = joinpath(dir, "v9.yaml")
    write(v9, "config_version: 99\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(v9; reload = true)
    # a case-scope file cannot pose as the general configuration
    vc = joinpath(dir, "vc.yaml")
    write(vc, "config_version: 1\nscope: case\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(vc; reload = true)

    # case configuration file: header binding and scope policing
    case_m = joinpath(dir, "case57.m")
    write(case_m, "% fixture\n")
    cc = Sparlectra.case_config_path(case_m)
    @test basename(cc) == "case57.config.yaml"
    write(cc, "config_version: 1\nscope: case\ncase: other.m\npower_flow:\n  max_iter: 43\n")
    err_mismatch = try
      Sparlectra.load_case_config(case_m)
      nothing
    catch err
      err
    end
    @test err_mismatch isa ArgumentError
    @test occursin("case_config_mismatch", sprint(showerror, err_mismatch))
    # a session-scope key in a case configuration is a hard error
    write(cc, "config_version: 1\nscope: case\ncase: case57.m\noutput:\n  console_summary: false\n")
    err_scope = try
      Sparlectra.load_case_config(case_m)
      nothing
    catch err
      err
    end
    @test err_scope isa ArgumentError
    @test occursin("not case scope", sprint(showerror, err_scope))
    # the canonical .scf.json double extension binds <stem>.config.yaml
    @test basename(Sparlectra.case_config_path(joinpath(dir, "case57.scf.json"))) == "case57.config.yaml"

    # precedence, one key on each level (D5), highest first: override,
    # case configuration file, deprecated in-file block, general file
    general = joinpath(dir, "general.yaml")
    write(general, "config_version: 1\nscope: general\npower_flow:\n  max_iter: 41\n")
    scf_case = joinpath(dir, "case57.scf.json")
    scf_root = string(
      "{\"version\": \"1.0\", \"type\": \"input\", \"is_batch\": false, \"attributes\": {}, \"data\": {}, ",
      "\"sparlectra\": {\"format_version\": \"", Sparlectra.SCF_FORMAT_VERSION, "\", \"config\": {\"power_flow.max_iter\": 42}}}",
    )
    write(scf_case, scf_root)
    write(cc, "config_version: 1\nscope: case\ncase: case57.scf.json\npower_flow:\n  max_iter: 43\n")
    with_override = Sparlectra.resolve_config(general, scf_case, Dict{String,Any}("power_flow.max_iter" => 44))
    @test with_override.config.powerflow.max_iter == 44
    with_case_file = Sparlectra.resolve_config(general, scf_case)
    @test with_case_file.config.powerflow.max_iter == 43
    rm(cc)
    with_block = Sparlectra.resolve_config(general, scf_case)
    @test with_block.config.powerflow.max_iter == 42
    # issue #1 point 1 (decided: defaults, format independent per stage-2
    # review): the PRESENCE of a case configuration file makes the case
    # self-contained, so a case-scope key the file does not set resolves
    # to the packaged default, NOT to the general file, for every input
    # format alike; machine-scope keys still come from the general file,
    # and a case WITHOUT a config file keeps the full chain
    write(general, "config_version: 1\nscope: general\npower_flow:\n  max_iter: 41\noutput:\n  console_max_rows: 33\n")
    write(scf_case, string("{\"version\": \"1.0\", \"type\": \"input\", \"is_batch\": false, \"attributes\": {}, \"data\": {}, ", "\"sparlectra\": {\"format_version\": \"", Sparlectra.SCF_FORMAT_VERSION, "\"}}"))
    tmpl_max_iter = Sparlectra.SparlectraConfig(Sparlectra.load_yaml_dict(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH)).powerflow.max_iter
    @test tmpl_max_iter != 41
    # no case config file: the general file participates, .scf.json included
    @test Sparlectra.resolve_config(general, scf_case).config.powerflow.max_iter == 41
    # with a case config file that sets ANOTHER key, the missing key falls
    # to the packaged default; the machine-scope key still applies
    write(cc, "config_version: 1\nscope: case\ncase: case57.scf.json\npower_flow:\n  tol: 1.0e-6\n")
    two_stage = Sparlectra.resolve_config(general, scf_case)
    @test two_stage.config.powerflow.max_iter == tmpl_max_iter
    @test two_stage.config.powerflow.tol == 1.0e-6
    @test two_stage.config.output.console_max_rows == 33
    # the same rule on a MATPOWER case: format independent
    plain = joinpath(dir, "plain57.m")
    write(plain, "% fixture\n")
    plain_cc = Sparlectra.case_config_path(plain)
    write(plain_cc, "config_version: 1\nscope: case\ncase: plain57.m\npower_flow:\n  tol: 1.0e-6\n")
    @test Sparlectra.resolve_config(general, plain).config.powerflow.max_iter == tmpl_max_iter
    rm(plain_cc)
    # a case level or an explicit override still beats the default
    write(cc, "config_version: 1\nscope: case\ncase: case57.scf.json\npower_flow:\n  max_iter: 43\n")
    @test Sparlectra.resolve_config(general, scf_case).config.powerflow.max_iter == 43
    @test Sparlectra.resolve_config(general, scf_case, Dict{String,Any}("power_flow.max_iter" => 44)).config.powerflow.max_iter == 44
    rm(cc)
    plain = joinpath(dir, "plain57.m")
    write(plain, "% fixture\n")
    with_general = Sparlectra.resolve_config(general, plain)
    @test with_general.config.powerflow.max_iter == 41
    # write_case_config round trip: header written, case scope enforced,
    # empty settings remove the file
    written = Sparlectra.write_case_config(plain, Dict{String,Any}("power_flow.max_iter" => 77); form = Dict{String,Any}("se_tol" => 1.0e-7))
    content = read(written, String)
    @test occursin("scope: case", content)
    @test occursin("case: plain57.m", content)
    @test Sparlectra.load_case_config(plain) == Dict{String,Any}("power_flow.max_iter" => 77)
    @test Sparlectra.resolve_config(general, plain).config.powerflow.max_iter == 77
    # YAML 1.1 lookalikes: the enum value "off" is quoted on write and
    # survives the round trip as a string (bare `off` came back as the
    # boolean false)
    written_off = Sparlectra.write_case_config(plain, Dict{String,Any}("model.auto_profile" => "off"))
    @test occursin("auto_profile: \"off\"", read(written_off, String))
    @test Sparlectra.load_case_config(plain)["model.auto_profile"] == "off"
    Sparlectra.write_case_config(plain, Dict{String,Any}())
    @test !isfile(written)
  end

  @testset "YAML scalar write/read round trip" begin
    # the writer's self-inverse contract: every string written by
    # _yaml_scalar_text must come back from the repository's own reader
    # with identical type and value (the bare enum value "off" once came
    # back as boolean false). Table per the review: the YAML 1.1 boolean
    # family, null spellings, number forms, dates, colon/hash strings,
    # symbol and list lookalikes, padded values, ordinary identifiers.
    for s in ("off", "on", "yes", "no", "Off", "TRUE", "false", "null", "~", "y", "n", "007", "-5", "+5", ".5", "1e3", "1.0", "0x1f", "2026-09-02", ":symbolish", "[1, 2]", " padded ", "with: colon", "hash # inside", "recommend", "case57.m", "excel_de")
      text = Sparlectra._yaml_scalar_text(s)
      back = Sparlectra.parse_yaml_scalar(text)
      @test back isa AbstractString
      @test String(back) == s
    end
    # non-strings keep their types through the same pair
    for v in (true, false, 42, 1.5)
      @test Sparlectra.parse_yaml_scalar(Sparlectra._yaml_scalar_text(v)) === v
    end
  end

  @testset "Test runner output mode helpers" begin
    @test selected_test_profile(String[], Dict("SPARLECTRA_TEST_PROFILE" => "extended")) === :extended
    @test selected_test_profile(["fast"], Dict("SPARLECTRA_TEST_PROFILE" => "extended")) === :fast
    @test selected_test_profile(["--verbose", "all"], Dict{String,String}()) === :all
    @test sparlectra_test_verbose(["--verbose"], Dict{String,String}())
    @test sparlectra_test_verbose(String[], Dict("SPARLECTRA_TEST_VERBOSE" => "1"))
    @test !sparlectra_test_verbose(String[], Dict{String,String}())

    quiet_path = test_scratch_path()
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

    # immediate failure surfacing: a captured "Test Failed at" block is
    # replayed right away instead of hiding until the suite's end (the
    # third swallowed-information fix after SKIPPED lines and stage 2);
    # driven with a SYNTHETIC capture, because a real failing @test here
    # would fail this suite
    fail_capture = test_scratch_path()
    write(fail_capture, string(
      "some solver output\n",
      "group name: Test Failed at /repo/test/test_x.jl:42\n",
      "  Expression: occursin(\"a\", \"b\")\n",
      "   Evaluated: occursin(\"a\", \"b\")\n",
      "\n",
      "more output after the block\n",
    ))
    fail_out = IOBuffer()
    @test _surface_test_failures(fail_capture, fail_out) == 1
    fail_text = String(take!(fail_out))
    @test occursin("surfaced immediately", fail_text)
    @test occursin("Test Failed at /repo/test/test_x.jl:42", fail_text)
    @test occursin("Evaluated: occursin", fail_text)
    @test !occursin("more output after the block", fail_text)
    clean_out = IOBuffer()
    write(fail_capture, "all green, nothing to surface\n")
    @test _surface_test_failures(fail_capture, clean_out) == 0
    @test isempty(String(take!(clean_out)))
    rm(fail_capture; force = true)

    verbose_path = test_scratch_path()
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
    cfgfile = test_scratch_path(".yaml")
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
      cfgfile = test_scratch_path(".yaml")
      write(cfgfile, """
power_flow:
  qlimits:
    enforcement_mode: $(mode)
""")
      cfg = Sparlectra.load_sparlectra_config(cfgfile; reload = true)
      @test cfg.powerflow.qlimits.enforcement_mode === mode
    end
    for (legacy, canonical) in ((:matpower_simultaneous, :classic_simultaneous), (:matpower_one_at_a_time, :classic_one_at_a_time))
      cfgfile = test_scratch_path(".yaml")
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
      "model" => Dict("auto_profile" => "off"),
      "matpower_import" => Dict("shift_unit" => "rad"),
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
      "model" => Dict("auto_profile" => "apply"),
      "matpower_import" => Dict("shift_unit" => "rad"),
    ))
    applied = Sparlectra.matpower_import_auto_profile(mpc, apply_cfg; mode = :apply)
    @test applied.config.matpower.shift_unit === :deg
    applied_row = only(row for row in applied.rows if row.option == "matpower_import.shift_unit")
    @test applied_row.action === :applied
    @test any(pair -> first(pair) === :shift_unit && last(pair) === :deg, applied.applied)
    convention_rows = Sparlectra._matpower_import_auto_profile_convention_scan(mpc)
    @test length(convention_rows) == 8

    ambiguous_cfg = Sparlectra.SparlectraConfig(Dict("model" => Dict("auto_profile" => "apply")))
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
      "model" => Dict("auto_profile" => "apply"),
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
    # The full block is what `output.console_auto_profile = :full` asks for.
    # It is no longer the default: a run where every check says "keep" used
    # to print about seventy console lines of MATPOWER option names, which
    # made a Sparlectra study read like a MATPOWER front end (maintainer,
    # 2026-09-05).
    verbose_cfg = Sparlectra.SparlectraConfig(Dict("output" => Dict("console_auto_profile" => "full")))
    io = IOBuffer()
    Sparlectra.write_matpower_import_auto_profile(io, conservative, verbose_cfg; casefile = "synthetic_fragile.m")
    conservative_text = String(take!(io))
    @test occursin("power_flow.start_mode.angle_mode", conservative_text)
    @test occursin("skipped", conservative_text)

    io = IOBuffer()
    Sparlectra.write_matpower_import_auto_profile(io, applied, verbose_cfg; casefile = "case13659pegase.m")
    apply_text = String(take!(io))
    @test occursin("Runtime casefile: case13659pegase.m", apply_text)
    @test occursin("Original MATPOWER import options:", apply_text)
    @test occursin("Auto-profile recommendation:", apply_text)
    @test occursin("Final effective MATPOWER import options:", apply_text)
    @test occursin("shift_unit   = deg", apply_text)

    # the DEFAULT is compact: one line when nothing changes, and only the
    # rows that actually recommend something when they do
    io = IOBuffer()
    Sparlectra.write_matpower_import_auto_profile(io, applied, applied.config; casefile = "case13659pegase.m")
    compact_text = String(take!(io))
    @test occursin("Import conventions:", compact_text)
    @test !occursin("Original MATPOWER import options", compact_text)
    @test !occursin("Runtime casefile", compact_text)
    @test occursin("shift_unit", compact_text)
    @test count(==('\n'), compact_text) <= 6
    # and :off says nothing at all
    quiet_cfg = Sparlectra.SparlectraConfig(Dict("output" => Dict("console_auto_profile" => "off")))
    io = IOBuffer()
    Sparlectra.write_matpower_import_auto_profile(io, applied, quiet_cfg; casefile = "case13659pegase.m")
    @test isempty(String(take!(io)))

    oom_cfg = Sparlectra.SparlectraConfig(Dict(
      "model" => Dict("auto_profile" => "apply"),
      "matpower_import" => Dict("shift_unit" => "rad"),
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
    tol_bad = test_scratch_path(".yaml")
    write(tol_bad, "power_flow:\n  tol: 0\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(tol_bad; reload = true)

    damping_bad = test_scratch_path(".yaml")
    write(damping_bad, "power_flow:\n  autodamp_min: 0\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(damping_bad; reload = true)

    angle_limit_bad = test_scratch_path(".yaml")
    write(angle_limit_bad, "power_flow:\n  start_mode:\n    dc_angle_limit_deg: 0\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(angle_limit_bad; reload = true)

    bool_ok_true = test_scratch_path(".yaml")
    write(bool_ok_true, "power_flow:\n  autodamp: true\n")
    @test Sparlectra.load_sparlectra_config(bool_ok_true; reload = true).powerflow.autodamp === true

    bool_ok_false = test_scratch_path(".yaml")
    write(bool_ok_false, "power_flow:\n  autodamp: false\n")
    @test Sparlectra.load_sparlectra_config(bool_ok_false; reload = true).powerflow.autodamp === false

    enum_bad = test_scratch_path(".yaml")
    write(enum_bad, "power_flow:\n  start_mode:\n    angle_mode: invalid_angle\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(enum_bad; reload = true)
    wrong_branch_bad = test_scratch_path(".yaml")
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
    flat_forced = test_scratch_path(".yaml")
    write(flat_forced, "cgmes_import:\n  start_values: flat\n")
    @test Sparlectra.load_sparlectra_config(flat_forced; reload = true).cgmes.start_values === :flat
    sv_ok = test_scratch_path(".yaml")
    write(sv_ok, "cgmes_import:\n  start_values: sv\n")
    @test Sparlectra.load_sparlectra_config(sv_ok; reload = true).cgmes.start_values === :sv
    sv_bad = test_scratch_path(".yaml")
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
    dcline_stale = test_scratch_path(".yaml")
    write(dcline_stale, "config_version: 1\nscope: general\nmatpower_import:\n  matpower_dcline_mode: reject_active\n")
    stale_cfg = @test_logs (:warn, r"reject_active is deprecated in configuration files") Sparlectra.load_sparlectra_config(dcline_stale; reload = true)
    @test stale_cfg.matpower.matpower_dcline_mode === :pf_injections
    dcline_deliberate = test_scratch_path(".yaml")
    write(dcline_deliberate, "matpower_import:\n  matpower_dcline_mode: ignore_inactive\n")
    @test Sparlectra.load_sparlectra_config(dcline_deliberate; reload = true).matpower.matpower_dcline_mode === :ignore_inactive

    # cgmes_import.placeholder_guards: default warn_skip, strict accepted,
    # invalid rejected with the key name in the message.
    @test Sparlectra.CGMESImportConfig().placeholder_guards === :warn_skip
    pg_ok = test_scratch_path(".yaml")
    write(pg_ok, "cgmes_import:\n  placeholder_guards: strict\n")
    @test Sparlectra.load_sparlectra_config(pg_ok; reload = true).cgmes.placeholder_guards === :strict
    pg_bad = test_scratch_path(".yaml")
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

  @testset "Tolerance can be stated in MW" begin
    # power_flow.tol_MW existed since task_tol_watts and was reachable only
    # by editing a YAML file, and even that failed because the template did
    # not carry the key. The form now states ONE value with a unit: two
    # fields side by side read like two competing tolerances (maintainer,
    # 2026-09-06).
    spec = Sparlectra._webui_option_spec("power_flow_tol_unit")
    @test spec.config_key === nothing
    @test spec.allowed_values == ("pu", "MW")
    html = Sparlectra.render_settings_page(output_root = mktempdir())
    @test occursin("name=\"power_flow_tol_unit\"", html)
    @test !occursin("power_flow_tol_mw", html)

    # pu (the default): the value becomes power_flow.tol and nothing else
    pu_req = Sparlectra.powerflow_webui_request(
      Dict{String,Any}("casefile" => "case14.m", "power_flow_tol" => "1e-8", "power_flow_tol_unit" => "pu");
      default_output_root = mktempdir())
    @test pu_req["config_overrides"]["power_flow.tol"] == 1.0e-8
    @test !haskey(pu_req["config_overrides"], "power_flow.tol_MW")

    # MW: the same field becomes power_flow.tol_MW, and the per-unit key is
    # NOT sent, so a run can never carry both bounds
    mw_req = Sparlectra.powerflow_webui_request(
      Dict{String,Any}("casefile" => "case14.m", "power_flow_tol" => "0.001", "power_flow_tol_unit" => "MW");
      default_output_root = mktempdir())
    @test mw_req["config_overrides"]["power_flow.tol_MW"] == 0.001
    @test !haskey(mw_req["config_overrides"], "power_flow.tol")

    # the key arrives in the typed configuration, from YAML and as override
    cfg = Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("tol_MW" => 1.0)))
    @test cfg.powerflow.tol_MW == 1.0
    @test Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.tol_MW" => 0.002))["power_flow"]["tol_MW"] == 0.002
  end

  @testset "Every configuration key arrives (task_config_arrival_v0100)" begin
    # Four settings were found in one day that exist, are documented, are
    # shown, and do not act. This closes the class instead of the cases: the
    # key list comes from the typed configuration itself, so a key added
    # later cannot be forgotten here.
    #
    # Two properties per key:
    #   REACHABLE - a user YAML carrying the key loads at all. Everything is
    #     validated against configuration.yaml.example, so a key missing
    #     THERE is refused as "unknown" no matter how well documented it is.
    #     That is how 22 state-estimation keys plus power_flow.tol_MW were
    #     unusable while their documentation described them in full.
    #   ARRIVES - a value that differs from the default reaches the typed
    #     configuration. Only checked where a differing value can be formed
    #     mechanically (numbers, booleans); symbols and collections have
    #     value domains a generic test must not invent, so for those the
    #     reachability check is the assertion.
    section_yaml = Dict("powerflow" => "power_flow", "matpower" => "matpower_import",
                        "cgmes" => "cgmes_import", "shortcircuit" => "short_circuit")
    is_leaf(v) = v isa Union{Symbol,AbstractString,Number,Bool,Nothing,AbstractVector,AbstractDict,Tuple}
    function collect_keys!(out, prefix, x)
      for f in fieldnames(typeof(x))
        v = getfield(x, f)
        path = isempty(prefix) ? String(f) : string(prefix, ".", f)
        is_leaf(v) ? push!(out, (path, v)) : collect_keys!(out, path, v)
      end
      return out
    end
    base = Sparlectra.SparlectraConfig(Dict{String,Any}())
    entries = Tuple{String,Any}[]
    for f in fieldnames(typeof(base))
      v = getfield(base, f)
      is_leaf(v) && continue
      String(f) == "user_set_keys" && continue
      collect_keys!(entries, get(section_yaml, String(f), String(f)), v)
    end

    # Explicit exclusions, one reason each. A silent exclusion or an empty
    # reason is not acceptable (same rule the reference-page coverage uses).
    # Three entries, and each one names a mechanism rather than an opinion.
    # The two `sparse` fields used to sit here as "internal, always true";
    # a field that can never be anything else is not a setting, so they were
    # deleted from the configuration structures instead (the old KEY is
    # still refused by name, so a stored file says what happened).
    excluded = Dict(
      "power_flow.qlimits.ignore_q_limits" => "internal inverse of power_flow.qlimits.enabled, never written by a user",
      "power_flow.start_mode.flatstart" => "set through power_flow.flatstart; the constructor files it under start_mode",
      "power_flow.distributed_slack.weights" => "free-form mapping (block style), present in the template and validated separately",
    )
    for (key, reason) in excluded
      @test !isempty(strip(reason))
      @test any(e -> first(e) == key, entries)
    end

    function nest_value(path, value)
      d::Any = value
      for p in Iterators.reverse(split(path, "."))
        d = Dict{String,Any}(String(p) => d)
      end
      return d
    end
    function load_with(path, value)
      file = test_scratch_path(".yaml")
      open(file, "w") do io
        println(io, "config_version: 1")
        Sparlectra._write_yaml_dict(io, nest_value(path, value isa Symbol ? String(value) : value))
      end
      return Sparlectra.load_sparlectra_config(file; reload = true)
    end
    read_key(cfg, path) = begin
      cur = cfg
      for p in split(path, ".")
        cur = getfield(cur, Symbol(p))
      end
      cur
    end
    # yaml path of a field path: the constructors file "guard_enabled" under
    # "guard.enabled" and so on, so the underscore variant is tried too
    function yaml_variants(path)
      parts = split(path, ".")
      out = [path]
      last = parts[end]
      i = findfirst('_', last)
      i === nothing || push!(out, join(vcat(parts[1:end-1], [last[1:i-1], last[i+1:end]]), "."))
      return out
    end

    allowed_by_key = Dict{String,Any}()
    for spec in Sparlectra.WEBUI_OPTION_SPECS
      spec.config_key === nothing && continue
      isempty(spec.allowed_values) && continue
      allowed_by_key[String(spec.config_key)] = spec.allowed_values
    end
    @test length(allowed_by_key) >= 15

    unreachable = String[]
    arrived = 0
    reach_only = 0
    for (path, value) in entries
      haskey(excluded, path) && continue
      reached = false
      for ypath in yaml_variants(path)
        try
          load_with(ypath, value)
          reached = true
        catch
        end
        reached && break
      end
      reached || (push!(unreachable, path); continue)
      # A differing value where one can be formed without inventing a
      # domain. Numbers and booleans are mechanical. For symbols the domain
      # is NOT guessed: it is read from the option spec that already
      # declares the allowed values, so a symbol key with a form field is
      # checked for arrival too (output.console_auto_profile, one of the
      # four triggers, is exactly such a key).
      spec_values = get(allowed_by_key, path, ())
      differing = value isa Bool ? !value :
                  value isa Integer ? value + 1 :
                  value isa AbstractFloat ? (isfinite(value) ? value + 1.0 : nothing) :
                  value isa Symbol ? begin
                    other = [v for v in spec_values if Symbol(String(v)) !== value]
                    isempty(other) ? nothing : String(first(other))
                  end : nothing
      if differing === nothing
        reach_only += 1
        continue
      end
      moved = false
      for ypath in yaml_variants(path)
        try
          cfg = load_with(ypath, differing)
          got = read_key(cfg, path)
          got isa Symbol && (got = String(got))
          isequal(got, differing) && (moved = true)
        catch
          # a domain validation may refuse default+1 (ordered knees, ranges);
          # the key is reachable, which is what this loop proves
        end
        moved && break
      end
      moved ? (arrived += 1) : (reach_only += 1)
    end

    @test isempty(unreachable)
    # the loop must actually cover the configuration, not a handful of keys
    @test length(entries) >= 200
    # 70 of the 232 reachable keys carry a number or a boolean, so a
    # differing value can be formed mechanically; the rest are symbols and
    # collections whose value domains a generic test must not invent
    @test arrived >= 60
    @test arrived + reach_only + length(excluded) == length(entries)

    # the four settings that started this task, by name
    for (path, differing, reader) in (
        ("power_flow.tol_MW", 0.002, cfg -> cfg.powerflow.tol_MW),
        ("state_estimation.k_suppress", 5.5, cfg -> cfg.state_estimation.k_suppress),
        ("state_estimation.rank_tol_factor", 25.0, cfg -> cfg.state_estimation.rank_tol_factor),
        ("output.console_auto_profile", "full", cfg -> cfg.output.console_auto_profile))
      cfg = load_with(path, differing)
      got = reader(cfg)
      @test (got isa Symbol ? String(got) : got) == differing
    end
  end

  @testset "Web UI keys are reachable in both directions" begin
    # tol_MW failed the first direction (a key the UI was allowed to set,
    # with no field to set it in) and stale fields fail the second. Both are
    # asserted, because each direction hides a different defect.
    specs = [s for s in Sparlectra.WEBUI_OPTION_SPECS if s.config_key !== nothing]
    have = Set(String(s.config_key) for s in specs)
    @test length(specs) >= 60

    # A: every form field maps to a key that exists AND may be set from the
    # GUI. A typo or a renamed key fails here instead of at run time.
    for spec in specs
      key = String(spec.config_key)
      @test key in Sparlectra.GUI_EDITABLE_CONFIG_KEYS
      probe = try
        Sparlectra.validate_gui_config_overrides(Dict{String,Any}(key => spec.default))
        true
      catch err
        # a value-domain complaint proves the KEY was accepted; only an
        # unknown-key or not-editable complaint is a finding here
        msg = sprint(showerror, err)
        !occursin("Unknown", msg) && !occursin("not allowed for GUI editing", msg)
      end
      @test probe
    end

    # B: every GUI-editable key either has a field or says why it has none.
    # Without this rule a key can be declared editable and never surface,
    # which is exactly what happened to power_flow.tol_MW.
    fieldless_reasons = Dict(
      "matpower_import.apply_branch_kind" => "import detail, set per case in the case configuration file",
      "matpower_import.apply_branch_names" => "import detail, set per case in the case configuration file",
      "matpower_import.import_for001_contingencies" => "DTF import detail, chosen by the FOR002 selection on the Case page",
      "model.net_cache_enabled" => "process-level cache switch, not a per-run choice",
      "power_flow.islands.diagnostic_continue_after_failure" => "island diagnostics detail, YAML and API only",
      "power_flow.islands.enabled" => "island solving follows the network, not a form choice",
      "power_flow.islands.mode" => "island solving follows the network, not a form choice",
      "power_flow.islands.reference_policy" => "island solving follows the network, not a form choice",
      "power_flow.method" => "the form offers the solver choice as power_flow.solver; method has one supported value",
      "short_circuit.sweep_method" => "short-circuit performance switch, YAML and API only",
      "power_flow.tol_MW" => "set through the tolerance VALUE field plus its unit selector (pu or MW), so one number cannot claim two units",
    )
    for key in sort(collect(Sparlectra.GUI_EDITABLE_CONFIG_KEYS))
      key in have && continue
      # the output.* group is console and CSV fine tuning: deliberately not
      # on the run form, otherwise every run would carry twenty display
      # switches. One reason for the whole group.
      startswith(key, "output.") && continue
      @test haskey(fieldless_reasons, key)
      haskey(fieldless_reasons, key) && @test !isempty(strip(fieldless_reasons[key]))
    end
    # and no stale reason: every named key must still be GUI-editable and
    # still have no field, otherwise the list rots
    for key in keys(fieldless_reasons)
      @test key in Sparlectra.GUI_EDITABLE_CONFIG_KEYS
      @test !(key in have)
    end
  end

  @testset "Form defaults do not drift from the code they mirror" begin
    # Step 5 of task_config_arrival_v0100: no numeric literal in a service or
    # API signature may duplicate something the configuration owns. After the
    # SE service moved to configuration-resolved keywords, the remaining
    # literals belong to the measurement GENERATOR and to the server start,
    # neither of which is a configuration key. They are still mirrored by
    # form defaults, so the two are pinned against each other here: a change
    # on one side without the other fails this test instead of silently
    # giving a run different numbers than the form promised.
    @test Sparlectra._webui_option_default("gen_passive_sigma") == 0.05
    @test Sparlectra._webui_option_default("gen_seed") == 42
    # and the ones that DO have a configuration key take it from there.
    # max_eliminations joined them after the review of 2026-09-06: a form
    # field whose value no configuration key can set is exactly the
    # asymmetry this task removed, so the elimination budget became
    # state_estimation.max_eliminations instead of a service literal.
    se = Sparlectra.state_estimation_config()
    @test Sparlectra._webui_option_default("se_max_eliminations") == se.max_eliminations
    @test Sparlectra._webui_option_default("se_k_suppress") == se.k_suppress
    @test Sparlectra._webui_option_default("se_k_eliminate") == se.k_eliminate
    @test Sparlectra._webui_option_default("se_max_iter") == se.max_iter
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
      cfg_bad = test_scratch_path(".yaml")
      write(cfg_bad, "diagnostics:\n  $(key): true\n")
      @test_throws ArgumentError Sparlectra.load_sparlectra_config(cfg_bad; reload = true)
    end
  end

  @testset "Stored configurations survive removed keys" begin
    # The startup warm-up is gone (0.10.0), but every stored user and Web UI
    # configuration written before that still carries `webui.warmup`. Without
    # the silent-removal entry each of them would fail to load with "Unknown
    # Sparlectra configuration key" - the Web UI would not start at all.
    cfg_path = test_scratch_path(".yaml")
    write(cfg_path, "webui:\n  warmup: true\n  operation_log_retention_days: 7\n")
    cfg = Sparlectra.load_sparlectra_config(cfg_path; reload = true)
    @test cfg.webui.operation_log_retention_days == 7
    @test !hasproperty(cfg.webui, :warmup)
    # as an override it follows the normal unknown-key path: it is no longer
    # a config key and not GUI-editable either
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("webui.warmup" => true))
  end

  @testset "Q-limit start mode public values" begin
    for mode in ("iteration", "auto", "iteration_or_auto")
      cfg = Sparlectra.SparlectraConfig(Dict(
        "power_flow" => Dict("qlimits" => Dict("start_mode" => mode)),
      ))
      @test cfg.powerflow.qlimits.start_mode === Symbol(mode)
    end
  end

  @testset "Configuration-derived network parameters reach every construction path" begin
    # The table from the task, executed. Four paths set these differently once,
    # and a case built through the wrong one ran with hysteresis 0 whatever the
    # configuration said. `bus_shunt_model` is consumed while the shunts are
    # built, so it is checked on the constructed network, not stamped.
    cfg = Sparlectra.SparlectraConfig(Dict(
      "power_flow" => Dict("qlimits" => Dict("cooldown_iters" => 3, "hysteresis_pu" => 0.05),
        "start_mode" => Dict("flatstart" => true)),
      "model" => Dict("bus_shunt_model" => "voltage_dependent_injection"),
    ))
    @test cfg.powerflow.qlimits.cooldown_iters == 3
    @test cfg.powerflow.qlimits.hysteresis_pu == 0.05

    # the stamping function itself: one source of truth for all paths
    stamped = Sparlectra.Net(name = "t", baseMVA = 100.0)
    Sparlectra._apply_config_net_parameters!(stamped, cfg)
    @test stamped.cooldown_iters == 3
    @test stamped.q_hyst_pu == 0.05
    @test stamped.flatstart

    # MATPOWER and the case format run through the shared import helper
    warmup = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m"))
    if isfile(warmup)
      mp = Sparlectra._import_sparlectra_net(warmup, nothing, cfg)
      @test mp.cooldown_iters == 3
      @test mp.q_hyst_pu == 0.05
      @test mp.bus_shunt_model === :voltage_dependent_injection
      d = mktempdir()
      scf = exportSCF(mp; file = joinpath(d, "params.scf.json"))
      sc = Sparlectra._import_sparlectra_net(scf, nothing, cfg)
      @test sc.cooldown_iters == 3
      @test sc.q_hyst_pu == 0.05
    end

    # DTF: the importer takes the model into its constructor, the switching
    # parameters are stamped afterwards, exactly as the service does it
    dtf = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "FOR001.DAT"))
    isfile(dtf) || println("      configuration coverage: DTF stamping leg SKIPPED (data/mpower/FOR001.DAT not present)")
    if isfile(dtf)
      dnet = Sparlectra.DTFImporter.build_net(Sparlectra.DTFImporter.read_dtf(dtf); bus_shunt_model = cfg.model.bus_shunt_model)
      @test dnet.bus_shunt_model === :voltage_dependent_injection
      Sparlectra._apply_config_net_parameters!(dnet, cfg)
      @test dnet.cooldown_iters == 3
      @test dnet.q_hyst_pu == 0.05
    end
  end

  @testset "Configuration refresh" begin
    stale = test_scratch_path(".yaml")
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
      p = test_scratch_path(".yaml")
      write(p, "power_flow:\n  qlimits:\n    enforcement_mode: $(legacy)\n")
      result = Sparlectra.refresh_sparlectra_config_file(p)
      @test "power_flow.qlimits.enforcement_mode" in result.normalized_keys
      @test occursin("enforcement_mode: $(canonical)", result.refreshed_text)
    end

    # regression 2026-09-02: refresh on a version-0 file must MOVE the
    # aliased model keys with their values, not stamp config_version 1 and
    # leave them behind (that wrote files the loader then rejected as
    # unknown keys)
    v0 = test_scratch_path(".yaml")
    write(v0, "matpower_import:\n  auto_profile_log: false\n  bus_shunt_model: voltage_dependent_injection\n  net_cache:\n    enabled: true\ntransformer:\n  tap_changer_model: impedance_correction\n")
    v0_result = Sparlectra.refresh_sparlectra_config_file(v0; write = true)
    @test v0_result.success
    @test "matpower_import.auto_profile_log" in v0_result.normalized_keys
    @test "transformer.tap_changer_model" in v0_result.normalized_keys
    v0_raw = Sparlectra.load_yaml_dict(v0)
    @test Sparlectra._dotted_config_get(v0_raw, "matpower_import.auto_profile_log") === nothing
    @test Sparlectra._dotted_config_get(v0_raw, "matpower_import.net_cache") === nothing
    # the emptied transformer section must not linger as an unknown key
    @test Sparlectra._dotted_config_get(v0_raw, "transformer") === nothing
    v0_cfg = Sparlectra.load_sparlectra_config(v0; reload = true)
    @test v0_cfg.model.auto_profile_log == false
    @test v0_cfg.model.bus_shunt_model === :voltage_dependent_injection
    @test v0_cfg.model.net_cache_enabled == true
    @test v0_cfg.model.tap_changer_model === :impedance_correction

    # half-migrated file as the broken refresh wrote it: config_version 1,
    # model block filled with template defaults, user values still under
    # the old names. The repair moves the user value over the default.
    hybrid = test_scratch_path(".yaml")
    write(hybrid, "config_version: 1\nmodel:\n  auto_profile_log: true\n  bus_shunt_model: admittance\nmatpower_import:\n  auto_profile_log: false\n  bus_shunt_model: voltage_dependent_injection\n")
    hybrid_result = Sparlectra.refresh_sparlectra_config_file(hybrid; write = true)
    @test hybrid_result.success
    @test "matpower_import.auto_profile_log" in hybrid_result.normalized_keys
    hybrid_cfg = Sparlectra.load_sparlectra_config(hybrid; reload = true)
    @test hybrid_cfg.model.auto_profile_log == false
    @test hybrid_cfg.model.bus_shunt_model === :voltage_dependent_injection

    # an explicitly set (non-default) target wins over the stale source
    explicit = test_scratch_path(".yaml")
    write(explicit, "config_version: 1\nmodel:\n  bus_shunt_model: voltage_dependent_injection\nmatpower_import:\n  bus_shunt_model: admittance\n")
    explicit_result = Sparlectra.refresh_sparlectra_config_file(explicit; write = true)
    @test explicit_result.success
    explicit_cfg = Sparlectra.load_sparlectra_config(explicit; reload = true)
    @test explicit_cfg.model.bus_shunt_model === :voltage_dependent_injection

    dup = test_scratch_path(".yaml")
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
    p = test_scratch_path(".yaml")
    write(p, "diagnostics:\n  console_diagnostics: full\n  console_max_rows: 50\n  log_effective_config: true\n")
    cfg = @test_logs (:warn, r"diagnostics\.console_diagnostics is deprecated") (:warn, r"diagnostics\.console_max_rows is deprecated") match_mode = :any Sparlectra.load_sparlectra_config(p; reload = true)
    @test cfg isa Sparlectra.SparlectraConfig
    @test cfg.diagnostics.log_effective_config
    # genuinely unknown keys still fail loudly
    bad = test_scratch_path(".yaml")
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

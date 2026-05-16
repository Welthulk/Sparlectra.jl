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

# file: test/test_matpower_example.jl

using Sparlectra
using Test

function _is_gitlab_ci()::Bool
  return lowercase(get(ENV, "GITLAB_CI", "false")) == "true"
end

function _ensure_matpower_case_for_example_tests(case_name::AbstractString)
  local_case = joinpath(Sparlectra.MPOWER_DIR, case_name)
  isfile(local_case) && return abspath(local_case)
  _is_gitlab_ci() && return nothing
  try
    return Sparlectra.FetchMatpowerCase.ensure_casefile(case_name; outdir = Sparlectra.MPOWER_DIR, to_jl = false, overwrite = false)
  catch
    return nothing
  end
end

function run_matpower_example_tests()
  @testset "MATPOWER example keyword forwarding" begin
    source = read(joinpath(@__DIR__, "..", "src", "examples", "matpower_import.jl"), String)
    normalized_source = replace(source, "\r\n" => "\n")
    signature_match = match(r"(?s)function bench_run_acpflow\(;.*?\n\)", normalized_source)
    @test !isnothing(signature_match)
    signature = signature_match.match

    expected_keywords = [
      "autodamp",
      "autodamp_min",
      "start_projection",
      "start_projection_try_dc_start",
      "start_projection_try_blend_scan",
      "start_projection_blend_lambdas",
      "start_projection_dc_angle_limit_deg",
      "qlimit_trace_buses",
      "qlimit_guard",
      "qlimit_guard_min_q_range_pu",
      "qlimit_guard_zero_range_mode",
      "qlimit_guard_narrow_range_mode",
      "qlimit_guard_max_switches",
      "qlimit_guard_freeze_after_repeated_switching",
      "qlimit_guard_accept_bounded_violations",
      "qlimit_guard_max_remaining_violations",
      "qlimit_guard_log",
      "matpower_ratio",
      "reference_override",
      "diagnose_branch_neighborhood",
      "diagnose_branch_neighborhood_buses",
      "diagnose_branch_neighborhood_depth",
      "diagnose_branch_neighborhood_maxlines",
      "diagnose_residual_clusters",
      "diagnose_residual_cluster_threshold_mw",
      "diagnose_residual_cluster_maxlines",
      "diagnose_nodal_balance_breakdown",
      "diagnose_nodal_balance_buses",
      "diagnose_nodal_balance_maxlines",
      "diagnose_negative_branch_impedance",
      "diagnose_negative_branch_impedance_fail_on_negative_r",
    ]

    for keyword in expected_keywords
      @test occursin(keyword, signature)
      @test occursin(keyword * " = " * keyword, normalized_source)
    end

    @test occursin("Base.invokelatest(", normalized_source)
    @test occursin(r"Base\.invokelatest\(\s*getfield\(@__MODULE__, :bench_run_acpflow\);", normalized_source)
    @test occursin("Base.invokelatest(getfield(@__MODULE__, :main))", normalized_source)
    @test occursin("SPARLECTRA_MATPOWER_IMPORT_NO_MAIN", normalized_source)
    @test occursin("function _enable_pq_gen_controllers_for_method(method::Symbol, requested::Bool)::Bool", normalized_source)
    @test occursin("return requested && method === :rectangular", normalized_source)
    @test occursin("!enable_pq_gen_controllers && m === :rectangular", normalized_source)
    @test occursin("function _print_converged_loss_summary(io::IO, method::Symbol, status, net::Sparlectra.Net)", normalized_source)
    @test occursin("losses P=", normalized_source)
    @test occursin("status_ref = status_ref", normalized_source)
    @test occursin("Base.invokelatest(getfield(@__MODULE__, :_print_pv_voltage_reference_diagnostics)", normalized_source)
    @test occursin("Sparlectra.Slack", normalized_source)
    @test occursin("Sparlectra.PV", normalized_source)
    @test occursin("matpower_auto_profile", normalized_source)
    @test occursin("function _matpower_auto_profile", normalized_source)
    @test occursin("console_summary::Bool", signature)
    @test occursin("console_auto_profile::Symbol", signature)
    @test occursin("console_q_limit_events::Symbol", signature)
    @test occursin("function _print_matpower_auto_profile_compact", normalized_source)
    @test occursin("function _print_matpower_run_summary", normalized_source)

    example_path = joinpath(@__DIR__, "..", "src", "examples", "matpower_import.jl")
    old_no_main = get(ENV, "SPARLECTRA_MATPOWER_IMPORT_NO_MAIN", nothing)
    ENV["SPARLECTRA_MATPOWER_IMPORT_NO_MAIN"] = "1"
    mod = Module(:MatpowerImportExampleSmoke)
    try
      Base.include(mod, example_path)

      example_yaml = joinpath(@__DIR__, "..", "src", "examples", "matpower_import.yaml.example")
      @test isfile(example_yaml)
      yaml_cfg = Base.invokelatest(() -> getfield(mod, :load_yaml_config)(example_yaml))
      example_case = String(yaml_cfg["case"])
      @test example_case == "case118.m"
      example_case_path = _ensure_matpower_case_for_example_tests(example_case)
      if isnothing(example_case_path)
        @test_skip "MATPOWER example case missing; dataset download skipped or unavailable"
      else
        @test isfile(example_case_path)
      end
      old_yaml = get(ENV, "SPARLECTRA_MATPOWER_IMPORT_YAML", nothing)
      try
        ENV["SPARLECTRA_MATPOWER_IMPORT_YAML"] = example_yaml
        @test Base.invokelatest(() -> getfield(mod, :_yaml_path_from_inputs)()) == example_yaml
      finally
        if isnothing(old_yaml)
          delete!(ENV, "SPARLECTRA_MATPOWER_IMPORT_YAML")
        else
          ENV["SPARLECTRA_MATPOWER_IMPORT_YAML"] = old_yaml
        end
      end

      cfg = Base.invokelatest(() -> getfield(mod, :bench_config_for_case)("case14.m", Dict{String,Any}(
        "autodamp" => true,
        "autodamp_min" => 0.002,
        "start_projection" => true,
        "start_projection_try_dc_start" => false,
        "start_projection_try_blend_scan" => false,
        "start_projection_blend_lambdas" => [0.1, 0.9],
        "start_projection_dc_angle_limit_deg" => 45.0,
        "qlimit_trace_buses" => [40712],
        "qlimit_guard" => true,
        "qlimit_guard_min_q_range_pu" => 2.5e-4,
        "qlimit_guard_zero_range_mode" => "lock_pq",
        "qlimit_guard_narrow_range_mode" => "prefer_pq",
        "qlimit_guard_max_switches" => 4,
        "qlimit_guard_freeze_after_repeated_switching" => false,
        "qlimit_guard_accept_bounded_violations" => true,
        "qlimit_guard_max_remaining_violations" => 2,
        "qlimit_guard_log" => false,
        "diagnose_matpower_reference" => true,
        "diagnose_branch_shift_conventions" => true,
        "diagnose_branch_neighborhood" => true,
        "diagnose_branch_neighborhood_buses" => [7934, 6621],
        "diagnose_branch_neighborhood_depth" => 2,
        "diagnose_branch_neighborhood_maxlines" => 25,
        "diagnose_residual_clusters" => true,
        "diagnose_residual_cluster_threshold_mw" => 5.0,
        "diagnose_residual_cluster_maxlines" => 4,
        "diagnose_nodal_balance_breakdown" => true,
        "diagnose_nodal_balance_buses" => [511, 2852],
        "diagnose_nodal_balance_maxlines" => 250,
        "diagnose_nodal_balance_include_branches" => false,
        "diagnose_nodal_balance_include_generators" => false,
        "diagnose_nodal_balance_include_shunts" => false,
        "diagnose_negative_branch_impedance" => true,
        "diagnose_negative_branch_impedance_maxlines" => 7,
        "diagnose_negative_branch_impedance_fail_on_negative_r" => false,
        "diagnose_negative_branch_impedance_fail_on_negative_x" => false,
        "diagnose_negative_branch_impedance_warn_threshold_abs_r" => 0.001,
        "diagnose_negative_branch_impedance_warn_threshold_abs_x" => 0.002,
        "diagnose_maxlines" => 3,
        "matpower_shift_unit" => "rad",
        "matpower_shift_sign" => -1.0,
        "matpower_ratio" => "reciprocal",
        "reference_override" => true,
        "reference_vm_pu" => 1.01,
        "reference_va_deg" => -3.0,
        "log_effective_config" => true,
        "matpower_auto_profile" => "recommend",
        "matpower_auto_profile_log" => true,
        "console_summary" => true,
        "console_auto_profile" => "compact",
        "console_diagnostics" => "summary",
        "console_q_limit_events" => "summary",
        "console_max_rows" => 7,
        "logfile_diagnostics" => "full",
      )))
      @test cfg.autodamp === true
      @test cfg.autodamp_min == 0.002
      @test cfg.start_projection === true
      @test cfg.start_projection_try_dc_start === false
      @test cfg.start_projection_try_blend_scan === false
      @test cfg.start_projection_blend_lambdas == [0.1, 0.9]
      @test cfg.start_projection_dc_angle_limit_deg == 45.0
      @test cfg.qlimit_trace_buses == [40712]
      @test cfg.qlimit_guard === true
      @test cfg.qlimit_guard_min_q_range_pu == 2.5e-4
      @test cfg.qlimit_guard_zero_range_mode === :lock_pq
      @test cfg.qlimit_guard_narrow_range_mode === :prefer_pq
      @test cfg.qlimit_guard_max_switches == 4
      @test cfg.qlimit_guard_freeze_after_repeated_switching === false
      @test cfg.qlimit_guard_accept_bounded_violations === true
      @test cfg.qlimit_guard_max_remaining_violations == 2
      @test cfg.qlimit_guard_log === false
      cfg_io = IOBuffer()
      Base.invokelatest(() -> getfield(mod, :_print_effective_config)(cfg_io, cfg))
      cfg_text = String(take!(cfg_io))
      @test occursin("qlimit_guard_min_q_range_pu: 0.00025", cfg_text)
      @test occursin("qlimit_guard_zero_range_mode: lock_pq", cfg_text)
      @test cfg.qlimit_guard_min_q_range_pu == 2.5e-4 # regression: non-default YAML guard value reaches the effective config used by bench_run_acpflow.
      @test_throws ErrorException Sparlectra.run_acpflow(casefile = "__missing_qlimit_guard_smoke__.m", qlimit_guard = cfg.qlimit_guard, qlimit_guard_min_q_range_pu = cfg.qlimit_guard_min_q_range_pu)
      @test cfg.diagnose_matpower_reference === true
      @test cfg.diagnose_branch_shift_conventions === true
      @test cfg.diagnose_branch_neighborhood === true
      @test cfg.diagnose_branch_neighborhood_buses == [7934, 6621]
      @test cfg.diagnose_branch_neighborhood_depth == 2
      @test cfg.diagnose_branch_neighborhood_maxlines == 25
      @test cfg.diagnose_residual_clusters === true
      @test cfg.diagnose_residual_cluster_threshold_mw == 5.0
      @test cfg.diagnose_residual_cluster_maxlines == 4
      @test cfg.diagnose_nodal_balance_breakdown === true
      @test cfg.diagnose_nodal_balance_buses == [511, 2852]
      @test cfg.diagnose_nodal_balance_maxlines == 250
      @test cfg.diagnose_nodal_balance_include_branches === false
      @test cfg.diagnose_nodal_balance_include_generators === false
      @test cfg.diagnose_nodal_balance_include_shunts === false
      @test cfg.diagnose_negative_branch_impedance === true
      @test cfg.diagnose_negative_branch_impedance_maxlines == 7
      @test cfg.diagnose_negative_branch_impedance_fail_on_negative_r === false
      @test cfg.diagnose_negative_branch_impedance_fail_on_negative_x === false
      @test cfg.diagnose_negative_branch_impedance_warn_threshold_abs_r == 0.001
      @test cfg.diagnose_negative_branch_impedance_warn_threshold_abs_x == 0.002
      @test cfg.diagnose_maxlines == 3
      @test cfg.matpower_shift_unit == "rad"
      @test cfg.matpower_shift_sign == -1.0
      @test cfg.matpower_ratio == "reciprocal"
      @test cfg.reference_override === true
      @test cfg.reference_vm_pu == 1.01
      @test cfg.reference_va_deg == -3.0
      @test cfg.log_effective_config === true
      @test cfg.matpower_auto_profile === :recommend
      @test cfg.matpower_auto_profile_log === true
      @test cfg.console_summary === true
      @test cfg.console_auto_profile === :compact
      @test cfg.console_diagnostics === :summary
      @test cfg.console_q_limit_events === :summary
      @test cfg.console_max_rows == 7
      @test cfg.logfile_diagnostics === :full
      @test occursin("matpower_shift_sign", normalized_source)
      @test occursin("matpower_ratio", normalized_source)
      @test occursin("reference_override", normalized_source)
      @test occursin("_print_reference_override_status", normalized_source)
      @test occursin("_print_effective_config", normalized_source)
      @test occursin("function _print_matpower_reference_diagnostics", normalized_source)
      @test occursin("function _print_matpower_branch_neighborhood_diagnostics", normalized_source)
      @test occursin("function _print_residual_cluster_diagnostics", normalized_source)
      @test occursin("function _print_negative_branch_impedance_diagnostics", normalized_source)
      @test occursin("function _print_nodal_balance_breakdown", normalized_source)
      @test occursin("Branch-shift convention scan", normalized_source)
      mktemp() do yaml_path, yaml_io
        write(yaml_io, "diagnose_nodal_balance_buses:\n  - 7934\n  - 6621\ndiagnose_negative_branch_impedance: true\n")
        close(yaml_io)
        block_cfg = Base.invokelatest(() -> getfield(mod, :load_yaml_config)(yaml_path))
        @test block_cfg["diagnose_nodal_balance_buses"] == [7934, 6621]
        @test block_cfg["diagnose_negative_branch_impedance"] === true
      end
      @test Base.invokelatest(() -> getfield(mod, :_enable_pq_gen_controllers_for_method)(:rectangular, true)) === true
      @test Base.invokelatest(() -> getfield(mod, :_enable_pq_gen_controllers_for_method)(:polar, true)) === false
      skipped_compare_summary = Base.invokelatest(() -> getfield(mod, :_show_once_summary_row)(:rectangular, (; converged = true), nothing, false; compare_available = false))
      @test skipped_compare_summary.converged === true
      @test skipped_compare_summary.iterations == -1
      @test isnan(skipped_compare_summary.elapsed_s)
      @test isnan(skipped_compare_summary.max_dvm)
      @test isnan(skipped_compare_summary.max_dva)
      @test isnan(skipped_compare_summary.slack_delta_va)
      @test skipped_compare_summary.angle_alignment === :none
      @test skipped_compare_summary.cmp_ok === false
      @test skipped_compare_summary.compare_status === :skip
      @test skipped_compare_summary.numerical_solution === :ok
      @test skipped_compare_summary.q_limit_active_set === :ok
      @test skipped_compare_summary.final_converged === true
      @test skipped_compare_summary.reason_text == "none"
      compared_summary = Base.invokelatest(() -> getfield(mod, :_show_once_summary_row)(:rectangular, (; converged = true, iterations = 4, elapsed_s = 0.125), (; max_dvm = 0.01, max_dva = 0.2), true; compare_available = true))
      @test compared_summary.converged === true
      @test compared_summary.iterations == 4
      @test compared_summary.elapsed_s == 0.125
      @test compared_summary.max_dvm == 0.01
      @test compared_summary.max_dva == 0.2
      @test isnan(compared_summary.slack_delta_va)
      @test compared_summary.angle_alignment === :none
      compared_summary_with_alignment = Base.invokelatest(() -> getfield(mod, :_show_once_summary_row)(:rectangular, (; converged = true, iterations = 4, elapsed_s = 0.125), (; max_dvm = 0.01, max_dva = 0.2, slack_delta_va = -30.0, angle_alignment = :slack), true; compare_available = true))
      @test compared_summary_with_alignment.slack_delta_va == -30.0
      @test compared_summary_with_alignment.angle_alignment === :slack
      @test compared_summary.cmp_ok === true
      @test compared_summary.compare_status === :ok
      rejected_summary = Base.invokelatest(() -> getfield(mod, :_show_once_summary_row)(:rectangular, (; converged = false, numerical_converged = true, q_limit_active_set_ok = false, final_converged = false, reason_text = "remaining PV Q-limit violations"), (; max_dvm = 0.036, max_dva = 0.44), true; compare_available = true))
      @test rejected_summary.numerical_solution === :ok
      @test rejected_summary.q_limit_active_set === :fail
      @test rejected_summary.final_converged === false
      @test rejected_summary.compare_status === :ok
      @test rejected_summary.reason_text == "remaining PV Q-limit violations"
      compared_warn_summary = Base.invokelatest(() -> getfield(mod, :_show_once_summary_row)(:rectangular, (; converged = true, iterations = 4, elapsed_s = 0.125), (; max_dvm = 0.01, max_dva = 0.2, compare_status = :warn), true; compare_available = true))
      @test compared_warn_summary.compare_status === :warn
      mpc_seeded = (; bus = hcat(collect(1.0:3.0), fill(1.0, 3), zeros(3, 5), [1.02, 1.01, 0.99], [0.0, -1.0, -2.0]))
      @test Base.invokelatest(() -> getfield(mod, :_warn_if_flatstart_uses_only_voltage_setpoints)("case1951rte.m", (; opt_flatstart = true), mpc_seeded)) === nothing
      @test_logs (:info, r"opt_flatstart=false uses stored MATPOWER voltage magnitudes and angles") Base.invokelatest(() -> getfield(mod, :_warn_if_flatstart_uses_only_voltage_setpoints)("case1951rte.m", (; opt_flatstart = false), mpc_seeded))

      large_bus = zeros(1000, 13)
      for r in axes(large_bus, 1)
        large_bus[r, 1] = r
        large_bus[r, 2] = r == 1 ? 3.0 : 1.0
        large_bus[r, 8] = 1.0
      end
      large_mpc = (; baseMVA = 100.0, bus = large_bus, gen = zeros(0, 21), branch = zeros(0, 13))
      auto_apply_cfg = Base.invokelatest(() -> getfield(mod, :bench_config_for_case)("case_large.m", Dict{String,Any}(
        "matpower_auto_profile" => "apply",
        "opt_flatstart" => false,
      )))
      auto_apply = Base.invokelatest(() -> getfield(mod, :_matpower_auto_profile)(large_mpc, auto_apply_cfg, Dict{String,Any}(
        "matpower_auto_profile" => "apply",
        "opt_flatstart" => false,
      )))
      @test auto_apply.mode === :apply
      # Explicit YAML values win; inferred large-case start-projection options are still applied.
      @test auto_apply.cfg.opt_flatstart === false
      @test :opt_flatstart in auto_apply.preserved
      @test auto_apply.cfg.start_projection === true
      @test auto_apply.cfg.flatstart_angle_mode === :dc
      @test auto_apply.cfg.flatstart_voltage_mode === :bus_vm_va_blend
      auto_recommend_cfg = Base.invokelatest(() -> getfield(mod, :bench_config_for_case)("case_large.m", Dict{String,Any}("matpower_auto_profile" => "recommend")))
      auto_recommend = Base.invokelatest(() -> getfield(mod, :_matpower_auto_profile)(large_mpc, auto_recommend_cfg, Dict{String,Any}("matpower_auto_profile" => "recommend")))
      @test auto_recommend.mode === :recommend
      @test auto_recommend.cfg.start_projection === false
      auto_io = IOBuffer()
      @test Base.invokelatest(() -> getfield(mod, :_print_matpower_auto_profile)(auto_io, auto_apply)) === nothing
      auto_text = String(take!(auto_io))
      @test occursin("MATPOWER auto-profile", auto_text)
      @test occursin("preserved explicit value", auto_text)
      auto_compact_io = IOBuffer()
      @test Base.invokelatest(() -> getfield(mod, :_print_matpower_auto_profile_compact)(auto_compact_io, auto_apply)) === nothing
      auto_compact_text = String(take!(auto_compact_io))
      @test occursin("MATPOWER auto-profile: apply", auto_compact_text)
      @test occursin("branch shift", auto_compact_text)
      @test !occursin("reason:", auto_compact_text)

      diagnostic_mpc = (;
        baseMVA = 100.0,
        bus = [1.0 3.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 110.0 1.0 1.1 0.9;
               2.0 1.0 0.0 0.0 0.0 0.0 1.0 1.0 -0.1 110.0 1.0 1.1 0.9],
        gen = [1.0 0.0 0.0 0.0 0.0 1.0 100.0 1.0 0.0 0.0 zeros(1, 11)...],
        branch = [1.0 2.0 0.01 0.10 0.02 0.0 0.0 0.0 1.0 0.01 1.0 0.0 0.0],
      )
      slack_ref = Base.invokelatest(() -> getfield(mod, :_matpower_slack_reference)(diagnostic_mpc))
      @test slack_ref.bus == 1
      @test slack_ref.vm_pu == 1.0
      @test slack_ref.va_deg == 0.0
      ref_io = IOBuffer()
      @test Base.invokelatest(() -> getfield(mod, :_print_reference_override_status)(ref_io, diagnostic_mpc; reference_override = false, reference_vm_pu = 1.1, reference_va_deg = 45.0)) === nothing
      @test occursin("override: disabled", String(take!(ref_io)))
      diag_io = IOBuffer()
      @test Base.invokelatest(() -> getfield(mod, :_print_matpower_reference_diagnostics)(diag_io, diagnostic_mpc; matpower_shift_sign = -1.0, matpower_shift_unit = "rad", diagnose_branch_shift_conventions = true, maxlines = 2)) === nothing
      @test occursin("Branch-shift convention scan", String(take!(diag_io)))
      neighborhood_io = IOBuffer()
      @test Base.invokelatest(() -> getfield(mod, :_print_matpower_reference_diagnostics)(neighborhood_io, diagnostic_mpc; matpower_shift_sign = -1.0, matpower_shift_unit = "rad", diagnose_branch_neighborhood = true, diagnose_branch_neighborhood_buses = [1, 2], diagnose_branch_neighborhood_depth = 1, diagnose_branch_neighborhood_maxlines = 5, maxlines = 2)) === nothing
      neighborhood_text = String(take!(neighborhood_io))
      @test occursin("MATPOWER branch-neighborhood fixed-reference diagnostics", neighborhood_text)
      @test occursin("branch row 1: f_bus=1  t_bus=2", neighborhood_text)
      @test occursin("connects two selected high-residual buses", neighborhood_text)
      @test occursin("fixed-reference flow Sf=", neighborhood_text)

      negative_mpc = (;
        baseMVA = 100.0,
        bus = [1.0 3.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 110.0 1.0 1.1 0.9;
               2.0 1.0 10.0 2.0 0.0 0.0 1.0 0.99 -1.0 110.0 1.0 1.1 0.9],
        gen = [1.0 10.0 2.0 10.0 -10.0 1.0 100.0 1.0 20.0 0.0 zeros(1, 11)...],
        branch = [1.0 2.0 -0.01 0.10 0.02 0.0 0.0 0.0 1.05 5.0 1.0 -360.0 360.0],
      )
      stamp = Base.invokelatest(() -> getfield(mod, :_matpower_branch_stamp)(view(negative_mpc.branch, 1, :)))
      expected_y = inv(-0.01 + 0.10im)
      @test stamp.br_r == -0.01
      @test isapprox(real(expected_y), -0.99009900990099; atol = 1e-14, rtol = 0.0)
      @test isapprox(stamp.Ytt - 0.01im, expected_y; atol = 1e-12, rtol = 0.0)
      ybus_negative = Sparlectra.MatpowerIO.build_ybus_matpower(negative_mpc.bus, negative_mpc.branch, negative_mpc.baseMVA)
      @test isapprox(ybus_negative[2, 2] - 0.01im, expected_y; atol = 1e-12, rtol = 0.0)
      negative_net = Sparlectra.createNetFromMatPowerCase(mpc = negative_mpc)
      @test negative_net.branchVec[1].r_pu == -0.01
      neg_io = IOBuffer()
      @test Base.invokelatest(() -> getfield(mod, :_print_negative_branch_impedance_diagnostics)(neg_io, negative_mpc; diagnose_branch_neighborhood_buses = [1, 2], diagnose_nodal_balance_buses = [2], fail_on_negative_r = false)) === nothing
      neg_text = String(take!(neg_io))
      @test occursin("Negative branch resistance/reactance detected", neg_text)
      @test occursin("<negative BR_R>", neg_text)
      @test occursin("Sparlectra preserves the signed impedance values", neg_text)
      @test occursin("y=1/(R+jX)=", neg_text)
      @test_throws ErrorException Base.invokelatest(() -> getfield(mod, :_print_negative_branch_impedance_diagnostics)(IOBuffer(), negative_mpc; fail_on_negative_r = true))
      nodal_diag = Base.invokelatest(() -> getfield(mod, :_matpower_reference_residuals)(negative_mpc))
      nodal_io = IOBuffer()
      @test Base.invokelatest(() -> getfield(mod, :_print_nodal_balance_breakdown)(nodal_io, negative_mpc, nodal_diag; buses = [1, 2], maxlines = 5)) === nothing
      nodal_text = String(take!(nodal_io))
      @test occursin("S_spec", nodal_text)
      @test occursin("S_branch_sum", nodal_text)
      @test occursin("residual dS=S_branch_sum-S_spec", nodal_text)
      @test occursin("formula: S_spec = S_gen_online - S_load - S_shunt", nodal_text)

      # BUS_I 196 uses a small equivalent Pd adjustment so this focused fixture
      # reproduces the known fixed-reference residual without embedding case300.m.
      case300_neighborhood_mpc = (;
        name = "case300_reference_neighborhood",
        baseMVA = 100.0,
        bus = [196.0 1.0 15.2823184571025 3.0 0.0 0.0 1.0 0.9695 -25.32 115.0 3.0 1.06 0.94;
               197.0 1.0 43.0 14.0 0.0 0.0 1.0 0.9907 -23.72 115.0 3.0 1.06 0.94;
               210.0 1.0 0.0 0.0 0.0 0.0 1.0 0.9788 -24.22 115.0 3.0 1.06 0.94;
               204.0 1.0 72.0 24.0 0.0 0.0 1.0 0.9718 -25.7 66.0 3.0 1.06 0.94;
               2040.0 1.0 0.0 0.0 0.0 0.0 1.0 0.9653 -14.94 115.0 3.0 1.06 0.94],
        gen = zeros(0, 21),
        branch = [196.0 197.0 0.014 0.04 0.004 0.0 0.0 0.0 0.0 0.0 1.0 -360.0 360.0;
                  196.0 210.0 0.03 0.081 0.01 0.0 0.0 0.0 0.0 0.0 1.0 -360.0 360.0;
                  204.0 2040.0 0.02 0.204 -0.012 0.0 0.0 0.0 1.05 0.0 1.0 -360.0 360.0;
                  196.0 2040.0 0.0001 0.02 0.0 0.0 0.0 0.0 1.0 0.0 1.0 -360.0 360.0],
      )
      case300_diag = Base.invokelatest(() -> getfield(mod, :_matpower_reference_residuals)(case300_neighborhood_mpc))
      busrow = Dict(Int(case300_diag.bus[r, 1]) => r for r in axes(case300_diag.bus, 1))
      @test isapprox(real(case300_diag.mis[busrow[2040]]), 9.269150; atol = 1e-6, rtol = 0.0)
      @test isapprox(real(case300_diag.mis[busrow[196]]), -9.260983; atol = 1e-6, rtol = 0.0)
      case300_diag_io = IOBuffer()
      @test Base.invokelatest(() -> getfield(mod, :_print_matpower_reference_diagnostics)(case300_diag_io, case300_neighborhood_mpc; maxlines = 3)) === nothing
      case300_diag_text = String(take!(case300_diag_io))
      @test occursin("raw MATPOWER-style fixed-reference check", case300_diag_text)
      @test occursin("not necessarily a Sparlectra solver or import error", case300_diag_text)
      cluster_io = IOBuffer()
      @test Base.invokelatest(() -> getfield(mod, :_print_matpower_reference_diagnostics)(cluster_io, case300_neighborhood_mpc; diagnose_residual_clusters = true, diagnose_residual_cluster_threshold_mw = 100.0, diagnose_residual_cluster_maxlines = 5, maxlines = 3)) === nothing
      cluster_text = String(take!(cluster_io))
      @test occursin("MATPOWER fixed-reference residual clusters", cluster_text)
      @test occursin("cluster 1:", cluster_text)
      @test occursin("196, 197, 2040", cluster_text)
      @test occursin("row 4 196->2040", cluster_text)

      pv_diag_mpc = Sparlectra.MatpowerIO.MatpowerCase(
        "pv_diagnostic",
        100.0,
        [
          1 3 0.0 0.0 0.0 0.0 1 1.00 0.0 110.0 1 1.1 0.9
          2 2 0.0 0.0 0.0 0.0 1 1.02 0.0 110.0 1 1.1 0.9
          3 1 20.0 5.0 0.0 0.0 1 0.98 -2.0 110.0 1 1.1 0.9
        ],
        [
          1 0.0 0.0 999.0 -999.0 1.0 100.0 1 999.0 0.0 zeros(1, 11)...
          2 50.0 0.0 999.0 -999.0 1.04 100.0 1 999.0 0.0 zeros(1, 11)...
        ],
        [
          1 2 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 1 -60.0 60.0
          2 3 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 1 -60.0 60.0
        ],
        nothing,
        nothing,
      )
      pv_diag_net = Sparlectra.createNetFromMatPowerCase(mpc = pv_diag_mpc, matpower_pv_voltage_source = :gen_vg)
      pv_diag_io = IOBuffer()
      @test Base.invokelatest(() -> getfield(mod, :_print_pv_voltage_reference_diagnostics)(pv_diag_io, pv_diag_mpc, pv_diag_net; matpower_pv_voltage_source = :gen_vg, compare_voltage_reference = :hybrid, maxlines = 3)) === nothing
      pv_diag_text = String(take!(pv_diag_io))
      @test occursin("Post-solve active PV/REF setpoint residual", pv_diag_text)
      @test occursin("BUS.VM fallback", pv_diag_text)

      warn_logfile, warn_io = mktemp()
      close(warn_io)
      warn_mpc = deepcopy(pv_diag_mpc)
      warn_mpc.bus[3, 2] = 2.0
      warn_status = Base.invokelatest(() -> getfield(mod, :_MatpowerImportLogStatus)())
      rows = Base.invokelatest(() -> getfield(mod, :_with_log_table)(warn_logfile, warn_status) do
        Sparlectra.MatpowerIO.pv_voltage_reference_rows(warn_mpc; warn = true)
      end)
      warn_text = read(warn_logfile, String)
      @test length(rows) == 3
      @test warn_status.warnings == 1
      @test warn_status.errors == 0
      @test occursin("Reported warnings/errors", warn_text)
      @test occursin("MATPOWER PV/REF buses without online generators: 1; using BUS.VM fallback", warn_text)
      @test occursin("BUS_I=3", warn_text)
      @test !occursin("has no online generator; falling back", warn_text)
      rm(warn_logfile; force = true)

      logfile, io = mktemp()
      close(io)
      result = Base.invokelatest(() -> getfield(mod, :bench_run_acpflow)(;
        casefile = "case14.m",
        methods = Symbol[],
        mpc = nothing,
        logfile = logfile,
        show_diff = false,
        tol_vm = 0.02,
        tol_va = 2.0,
        autodamp = true,
        autodamp_min = 0.002,
        start_projection = true,
        start_projection_try_dc_start = false,
        start_projection_try_blend_scan = false,
        start_projection_blend_lambdas = [0.1, 0.9],
        start_projection_dc_angle_limit_deg = 45.0,
        qlimit_trace_buses = [40712],
        qlimit_guard = cfg.qlimit_guard,
        qlimit_guard_min_q_range_pu = cfg.qlimit_guard_min_q_range_pu,
        qlimit_guard_zero_range_mode = cfg.qlimit_guard_zero_range_mode,
        qlimit_guard_narrow_range_mode = cfg.qlimit_guard_narrow_range_mode,
        qlimit_guard_max_switches = cfg.qlimit_guard_max_switches,
        qlimit_guard_freeze_after_repeated_switching = cfg.qlimit_guard_freeze_after_repeated_switching,
        qlimit_guard_accept_bounded_violations = cfg.qlimit_guard_accept_bounded_violations,
        qlimit_guard_max_remaining_violations = cfg.qlimit_guard_max_remaining_violations,
        qlimit_guard_log = cfg.qlimit_guard_log,
        matpower_shift_unit = "rad",
        matpower_shift_sign = -1.0,
        matpower_ratio = "reciprocal",
        reference_override = true,
        reference_vm_pu = 1.01,
        reference_va_deg = -3.0,
        diagnose_branch_neighborhood = true,
        diagnose_residual_clusters = true,
        diagnose_residual_cluster_threshold_mw = 5.0,
        diagnose_residual_cluster_maxlines = 4,
        diagnose_branch_neighborhood_buses = [7934, 6621],
        diagnose_branch_neighborhood_depth = 2,
        diagnose_branch_neighborhood_maxlines = 25,
        diagnose_nodal_balance_breakdown = true,
        diagnose_nodal_balance_buses = [511, 2852],
        diagnose_nodal_balance_maxlines = 2,
        diagnose_negative_branch_impedance = true,
        diagnose_negative_branch_impedance_fail_on_negative_r = false,
        log_effective_config = true,
        effective_config = cfg,
        show_once = false,
        benchmark = false,
        console_summary = true,
        console_diagnostics = :compact,
        console_q_limit_events = :summary,
        console_max_rows = 3,
        logfile_diagnostics = :full,
      ))
      @test result == Dict{Symbol,Any}()
      rm(logfile; force = true)
    finally
      if isnothing(old_no_main)
        delete!(ENV, "SPARLECTRA_MATPOWER_IMPORT_NO_MAIN")
      else
        ENV["SPARLECTRA_MATPOWER_IMPORT_NO_MAIN"] = old_no_main
      end
    end
  end

  @testset "MATPOWER SHIFT convention options" begin
    bus = [1.0 3.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 110.0 1.0 1.1 0.9;
           2.0 1.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 110.0 1.0 1.1 0.9]
    branch_rad = [1.0 2.0 0.01 0.10 0.02 0.0 0.0 0.0 1.0 (pi / 6) 1.0 -360.0 360.0]
    branch_deg = copy(branch_rad)
    branch_deg[1, 10] = -30.0
    y_rad = Sparlectra.MatpowerIO.build_ybus_matpower(bus, branch_rad, 100.0; matpower_shift_unit = "rad", matpower_shift_sign = -1.0)
    y_deg = Sparlectra.MatpowerIO.build_ybus_matpower(bus, branch_deg, 100.0)
    @test isapprox(y_rad, y_deg; atol = 1e-12, rtol = 0.0)
  end

  @testset "MATPOWER nominal TAP transformer preservation" begin
    bus = [196.0 3.0 10.0 3.0 0.0 0.0 1.0 0.9695 -25.32 115.0 3.0 1.06 0.94;
           2040.0 1.0 0.0 0.0 0.0 0.0 1.0 0.9653 -14.94 115.0 3.0 1.06 0.94;
           3000.0 1.0 0.0 0.0 0.0 0.0 1.0 0.9700 -20.00 115.0 3.0 1.06 0.94]
    gen = [196.0 0.0 0.0 10.0 -10.0 1.0 100.0 1.0 10.0 0.0 zeros(1, 11)...]
    branch = [196.0 2040.0 0.0001 0.02 0.0 0.0 0.0 0.0 1.0 0.0 1.0 -360.0 360.0;
              2040.0 3000.0 0.0001 0.02 0.0 0.0 0.0 0.0 0.0 0.0 1.0 -360.0 360.0]
    mpc = (; name = "nominal_tap_transformer", baseMVA = 100.0, bus = bus, gen = gen, branch = branch, gencost = nothing, bus_name = nothing)

    net = Sparlectra.createNetFromMatPowerCase(mpc = mpc)
    @test length(net.trafos) == 1
    @test length(net.linesAC) == 1
    @test net.branchVec[1].ratio == 1.0
    @test net.branchVec[1].tap_ratio == 1.0
    @test net.branchVec[2].ratio == 0.0
    @test isapprox(Matrix(Sparlectra.createYBUS(net = net, sparse = false)), Sparlectra.MatpowerIO.build_ybus_matpower(bus, branch, 100.0); atol = 1e-12, rtol = 0.0)
  end

  @testset "MATPOWER TAP ratio convention options" begin
    bus = [1.0 3.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 110.0 1.0 1.1 0.9;
           2.0 1.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 110.0 1.0 1.1 0.9]
    branch = [1.0 2.0 0.01 0.10 0.02 0.0 0.0 0.0 2.0 0.0 1.0 0.0 0.0]
    branch_recip = copy(branch)
    branch_recip[1, 9] = 0.5
    y_recip_option = Sparlectra.MatpowerIO.build_ybus_matpower(bus, branch, 100.0; matpower_ratio = "reciprocal")
    y_recip_data = Sparlectra.MatpowerIO.build_ybus_matpower(bus, branch_recip, 100.0)
    @test isapprox(y_recip_option, y_recip_data; atol = 1e-12, rtol = 0.0)

    mpc = (;
      name = "ratio_convention",
      baseMVA = 100.0,
      bus = bus,
      gen = [1.0 0.0 0.0 10.0 -10.0 1.0 100.0 1.0 10.0 0.0 zeros(1, 11)...],
      branch = branch,
      gencost = nothing,
      bus_name = nothing,
    )
    net_normal = Sparlectra.createNetFromMatPowerCase(mpc = mpc, matpower_ratio = "normal")
    net_recip = Sparlectra.createNetFromMatPowerCase(mpc = mpc, matpower_ratio = "reciprocal")
    @test net_normal.branchVec[1].ratio == 2.0
    @test net_recip.branchVec[1].ratio == 0.5
  end

  @testset "MATPOWER local Julia case resolution" begin
    mktempdir() do dir
      case_path = joinpath(dir, "local_case.jl")
      write(case_path, """
      (
        name = \"local_case\",
        baseMVA = 100.0,
        bus = [1 3 0 0 0 0 1 1.0 0.0 110 1 1.1 0.9],
        gen = [1 0 0 10 -10 1.0 100 1 10 0],
        branch = zeros(0, 13),
        gencost = nothing,
        bus_name = nothing,
      )
      """)

      resolved = Sparlectra.FetchMatpowerCase.ensure_casefile("local_case.jl"; outdir = dir, to_jl = true, overwrite = false)
      @test resolved == abspath(case_path)

      old_pwd = pwd()
      try
        cd(dirname(@__DIR__))
        relative_path = relpath(case_path, pwd())
        mpc = Sparlectra.MatpowerIO.read_case(relative_path; legacy_compat = false)
        @test mpc.name == "local_case"
        @test size(mpc.bus) == (1, 13)
      finally
        cd(old_pwd)
      end
    end
  end
  return true
end

if abspath(PROGRAM_FILE) == @__FILE__
  Base.invokelatest(run_matpower_example_tests)
end

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

# file: test/test_legacy_configuration_coverage.jl
# purpose: legacy profile: the testsets of test/test_configuration_coverage.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_configuration_coverage.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function test_configuration_qlimit_enforcement_keys()
  @testset "Q-limit enforcement mode user YAML keys" begin (function ()
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
    # `off` is the disabled state, not a mode: the Web UI settings pages up
    # to 0.11.1 saved the word, and the loader rejected the file afterwards
    let cfgfile = test_scratch_path(".yaml")
      write(cfgfile, "power_flow:\n  qlimits:\n    enabled: true\n    enforcement_mode: off\n")
      cfg = Sparlectra.load_sparlectra_config(cfgfile; reload = true)
      @test cfg.powerflow.qlimits.ignore_q_limits
      @test cfg.powerflow.qlimits.enforcement_mode === :active_set
      @test Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.qlimits.enforcement_mode" => "off"))["power_flow"]["qlimits"]["enforcement_mode"] == "off"
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
    # #464: power_flow.qlimits.classic_max_passes reaches the solver keyword
    # qlimit_max_outer from a configuration file. The shipped Zeng case needs
    # two Q-limit updates in classic_one_at_a_time. One pass stops the loop at
    # the limit: the second update found there is not applied, one pass is
    # counted, and the run is not converged (status reason
    # max_outer_iterations, the violation left in the final Q check). Two
    # passes converge.
    zeng = joinpath(dirname(@__DIR__), "data", "scf", "case14_zeng_p306_one_at_a_time_C.scf.json")
    quiet = Sparlectra.OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)
    for (passes, outcome, converged, reason) in ((1, :max_outer_iterations, false, :max_outer_iterations), (2, :converged, true, :none))
      cfgfile = test_scratch_path(".yaml")
      write(cfgfile, "power_flow:\n  rescue: false\n  qlimits:\n    enforcement_mode: classic_one_at_a_time\n    classic_max_passes: $(passes)\n")
      cfg = Sparlectra.load_sparlectra_config(cfgfile; reload = true)
      @test cfg.powerflow.qlimits.classic_max_passes == passes
      r = run_sparlectra(net = Sparlectra.importSCF(zeng), config = Sparlectra.SparlectraConfig(powerflow = cfg.powerflow, output = quiet))
      st = Sparlectra.rectangular_pf_status(r.net)
      # passes counted and clamps applied are the same number: one per pass
      @test (st.final_outcome, st.matpower_outer_iterations, length(r.net.qLimitLog)) == (outcome, passes, passes)
      @test (r.final_converged, r.reason) == (converged, reason)
      @test r.numerical_converged
      @test (st.final_q_check_status === :remaining_pv_q_limit_violations) == !converged
    end
    @test Sparlectra.QLimitConfig().classic_max_passes == 30
    @test_throws ArgumentError Sparlectra.QLimitConfig(Dict{String,Any}("classic_max_passes" => 0))
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.qlimits.classic_max_passes" => 0))
    err = try
      Sparlectra.QLimitConfig(Dict("enforcement_mode" => "definitely_not_supported"))
      nothing
    catch caught
      caught
    end
    @test err isa ArgumentError
    @test occursin("classic_simultaneous", sprint(showerror, err))
    # one default per key: the packaged template, the struct (a library run
    # without a configuration file) and the loader fallback (a user file
    # that omits the key) agree field by field; guard.max_switches was 3 in
    # the template and 10 in the struct up to 0.20.2
    tmpl = Sparlectra.SparlectraConfig(Sparlectra.load_yaml_dict(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH)).powerflow.qlimits
    from_struct = Sparlectra.QLimitConfig()
    from_empty = Sparlectra.QLimitConfig(Dict{String,Any}())
    for f in fieldnames(Sparlectra.QLimitConfig)
      @test (f, getfield(from_struct, f)) == (f, getfield(tmpl, f))
      @test (f, getfield(from_empty, f)) == (f, getfield(tmpl, f))
    end
    # 0.30.0 agreement for the ten release keys: template = struct = loader
    # fallback = the keyword defaults the solver entry points resolve to
    # (read back from the performance profile of a config-less library run
    # on the shipped sp_case118) = the effective configuration of a service
    # run. The keyword defaults read DEFAULT_QLIMIT_CONFIG and
    # DEFAULT_NEWTON_UPDATE; the net and the MATPOWER importer read the
    # hysteresis and cooldown from the same constant, so final_q_accept_pu
    # ("auto" = 2 * hysteresis) agrees too.
    tmpl_cfg = Sparlectra.SparlectraConfig(Sparlectra.load_yaml_dict(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH))
    empty_cfg = Sparlectra.SparlectraConfig(Dict{String,Any}("config_version" => 1))
    expected = (newton_update = :polar, start_iter = 3, start_mode = :iteration_or_auto, guard = true, guard_min_q_range_pu = 0.02,
      guard_narrow_range_mode = :lock_pq, guard_violation_mode = :lock_pq, guard_max_switches = 3, final_q_accept_pu = 0.02, auto_profile = :off)
    for cfg in (tmpl_cfg, empty_cfg, Sparlectra.SparlectraConfig())
      q = cfg.powerflow.qlimits
      @test (cfg.powerflow.newton_update, q.start_iter, q.start_mode, q.guard, q.guard_min_q_range_pu, q.guard_narrow_range_mode, q.guard_violation_mode, q.guard_max_switches, q.final_q_accept_pu, cfg.model.auto_profile) == values(expected)
    end
    @test Sparlectra.DEFAULT_NEWTON_UPDATE === expected.newton_update
    # the batch path (N-1 and scenarios through the contingency service)
    # takes its Q-limit keywords from the same mapping as the single run
    # (_qlimit_solver_kwargs); before 0.30.2 it forwarded none of them. The
    # template and an empty configuration map to the release defaults here,
    # and the service test "scenario service sources (step 5)" shows that a
    # configured lock list changes an N-1 result.
    for cfg in (tmpl_cfg, empty_cfg)
      qk = Sparlectra._qlimit_solver_kwargs(cfg.powerflow.qlimits)
      @test (qk.qlimit_start_iter, qk.qlimit_start_mode, qk.qlimit_guard, qk.qlimit_guard_min_q_range_pu, qk.qlimit_guard_narrow_range_mode, qk.qlimit_guard_violation_mode, qk.qlimit_guard_max_switches) ==
            (expected.start_iter, expected.start_mode, expected.guard, expected.guard_min_q_range_pu, expected.guard_narrow_range_mode, expected.guard_violation_mode, expected.guard_max_switches)
    end
    lib_net = Sparlectra.createNetFromMatPowerFile(filename = joinpath(dirname(@__DIR__), "data", "mpower", "sp_case118.m"))
    prof = Dict{Symbol,Any}()
    redirect_stdout(devnull) do
      runpf!(lib_net, 30, 1e-8, 0; performance_profile = prof)
    end
    @test (prof[:newton_update], prof[:qlimit_start_iter], prof[:qlimit_start_mode], prof[:qlimit_guard], prof[:qlimit_guard_min_q_range_pu], prof[:qlimit_guard_narrow_range_mode], prof[:qlimit_guard_violation_mode], prof[:qlimit_guard_max_switches], prof[:final_q_accept_pu]) ==
          (expected.newton_update, expected.start_iter, expected.start_mode, expected.guard, expected.guard_min_q_range_pu, expected.guard_narrow_range_mode, expected.guard_violation_mode, expected.guard_max_switches, expected.final_q_accept_pu)
    # 0.30.2: the three dishonest-Newton keys agree the same way (template,
    # struct, loader fallback, solver keyword defaults)
    jr_expected = (false, 10.0, 3)
    for cfg in (tmpl_cfg, empty_cfg, Sparlectra.SparlectraConfig())
      pf = cfg.powerflow
      @test (pf.jacobian_reuse, pf.jacobian_reuse_min_reduction, pf.jacobian_reuse_max_steps) == jr_expected
    end
    @test (Sparlectra.DEFAULT_JACOBIAN_REUSE, Sparlectra.DEFAULT_JACOBIAN_REUSE_MIN_REDUCTION, Sparlectra.DEFAULT_JACOBIAN_REUSE_MAX_STEPS) == jr_expected
    @test (prof[:jacobian_reuse], prof[:jacobian_reuse_min_reduction], prof[:jacobian_reuse_max_steps]) == jr_expected
    @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict{String,Any}("config_version" => 1, "power_flow" => Dict{String,Any}("jacobian_reuse_min_reduction" => 1.0)))
    @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict{String,Any}("config_version" => 1, "power_flow" => Dict{String,Any}("jacobian_reuse_max_steps" => 0)))
    # Deliberate exception to the one-default rule above, which covers the
    # solver behaviour: output.logfile_results is log verbosity. A library run
    # without a configuration file stays quiet (off), the packaged template
    # that the Web UI and the services start from logs the full result table.
    # Pinned here so the difference cannot drift unnoticed.
    @test Sparlectra.OutputConfig().logfile_results === :off
    @test empty_cfg.output.logfile_results === :off
    @test tmpl_cfg.output.logfile_results === :full
    # power_flow.power_mode_lu (0.30.2) the same way, and an unknown value is
    # refused by the loader and the override validation
    for cfg in (tmpl_cfg, empty_cfg, Sparlectra.SparlectraConfig())
      @test cfg.powerflow.power_mode_lu === :auto
    end
    @test Sparlectra.DEFAULT_POWER_MODE_LU === :auto
    @test prof[:power_mode_lu] === :auto
    @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict{String,Any}("config_version" => 1, "power_flow" => Dict{String,Any}("power_mode_lu" => "pardiso")))
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.power_mode_lu" => "pardiso"))
  end)() end
  return nothing
end

function test_configuration_value_domain_validation()
  @testset "Configuration value-domain validation" begin (function ()
    tol_bad = test_scratch_path(".yaml")
    write(tol_bad, "power_flow:\n  tol: 0\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(tol_bad; reload = true)

    damping_bad = test_scratch_path(".yaml")
    write(damping_bad, "power_flow:\n  autodamp_min: 0\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(damping_bad; reload = true)

    fit_bad = test_scratch_path(".yaml")
    write(fit_bad, "model:\n  auto_profile_max_fit_pu: -0.1\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(fit_bad; reload = true)

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

    # power_flow.plausibility.vmin_warn: 0 switches the check off and is
    # valid, a negative threshold is refused in a file and as an override
    @test Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.plausibility.vmin_warn" => 0.0))["power_flow"]["plausibility"]["vmin_warn"] == 0.0
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.plausibility.vmin_warn" => -0.1))
    vw_bad = test_scratch_path(".yaml")
    write(vw_bad, "power_flow:\n  plausibility:\n    vmin_warn: -0.1\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(vw_bad; reload = true)
  end)() end
  return nothing
end

function test_configuration_webui_keys_both_directions()
  @testset "Web UI keys are reachable in both directions" begin (function ()
    # tol_MW failed the first direction (a key the UI was allowed to set,
    # with no field to set it in) and stale fields fail the second. Both are
    # asserted, because each direction hides a different defect.
    # collect FIRST, filter second. WEBUI_OPTION_SPECS is an NTuple of 97
    # elements: Julia unrolls tuple iteration, and filtering makes the result
    # length unknown, so a comprehension straight off the tuple ends in
    # Base.grow_to! over a 97-way unrolled generator. Measured,
    # that single line cost 5.61 s of the 5.93 s this group spent inferring
    # generator machinery, and it was the largest item in the whole group.
    # collect() on a tuple with a concrete eltype knows both length and type,
    # and filter() on the resulting Vector is one specialization.
    specs = filter(s -> s.config_key !== nothing, collect(SparlectraApp.WEBUI_OPTION_SPECS))
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
      "power_flow.islands.reference_policy" => "island solving follows the network, not a form choice",
      "power_flow.method" => "the form offers the solver choice as power_flow.solver; method has one supported value",
      "short_circuit.sweep_method" => "short-circuit performance switch, YAML and API only",
      "short_circuit.converter_current_factor" => "converter fault current without maxIFault, a modelling choice per case; set via configuration.yaml, the case sidecar or the API",
      # the benchmark mode belongs to run_matpower_case in a script (with
      # `using BenchmarkTools`); the Web UI has no benchmark option since 0.20.2
      "benchmark.enabled" => "benchmark mode of run_matpower_case, scripts only; no Web UI option",
      "benchmark.samples" => "benchmark mode of run_matpower_case, scripts only; no Web UI option",
      "benchmark.seconds" => "benchmark mode of run_matpower_case, scripts only; no Web UI option",
      "power_flow.tol_MW" => "set through the tolerance VALUE field plus its unit selector (pu or MW), so one number cannot claim two units",
      "state_estimation.robust" => "issue #377: case scope like the other estimator options, but no dedicated UI control (the form exposes robust_mode instead); set via configuration.yaml or the case sidecar",
      "state_estimation.topology_precheck" => "issue #377: case scope, no dedicated UI control; set via configuration.yaml or the case sidecar",
      "state_estimation.linear_solver" => "estimator linear-algebra backend, a performance switch; set via configuration.yaml or the case sidecar",
      "state_estimation.symmetric_gain" => "estimator linear-algebra switch (Cholesky route for the gain matrix), a performance choice; set via configuration.yaml or the case sidecar",
      "qsts.mode" => "the Web UI QSTS button sets the mode for its own run only; the key itself is set via configuration.yaml or the case sidecar",
      "qsts.base_step" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.snapshot_controls" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.sequential_start" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.resync_every" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.dead_time_steps" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.hunting_window" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.vm_min_pu" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.vm_max_pu" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.loading_max_pct" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.series.min_structural_overlap" => "series check of CGMES snapshots, a library function without a Web UI control; set via configuration.yaml or the case sidecar",
      "qsts.series.max_injection_churn" => "series check of CGMES snapshots, a library function without a Web UI control; set via configuration.yaml or the case sidecar",
      "qsts.series.max_gap" => "series check of CGMES snapshots, a library function without a Web UI control; set via configuration.yaml or the case sidecar",
      "cgmes_import.map_energy_source" => "CGMES import switch for deliveries that carry generation as EnergySource; set via configuration.yaml or the case sidecar",
      "cgmes_import.map_power_electronics" => "CGMES import switch for CGMES 3.0 deliveries with PowerElectronicsConnection; set via configuration.yaml or the case sidecar",
      "cgmes_import.shunt_control" => "CGMES import switch for the bank regulators of a delivery; set via configuration.yaml or the case sidecar",
      "power_flow.plausibility.vmin_warn" => "threshold of the plausibility warning after a converged run; set via configuration.yaml or the case sidecar",
      "state_estimation.sigma_by_voltage_level" => "a table, not a scalar field: the sigma table of the State Estimation section posts it (se_sigma_lvl_* fields), the run takes it as an override and Save settings writes it into the case configuration",
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
  end)() end
  return nothing
end

function test_configuration_removed_diagnostics_rejected()
  @testset "Removed diagnostics keys warn and are ignored" begin (function ()
    # a stored file that still carries one of the removed diagnostics
    # switches loads; the key is warned about by name and dropped
    # (unknown or obsolete keys never fail a run)
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
      write(cfg_bad, "diagnostics:\n  $(key): true\n  log_effective_config: true\n")
      cfg = @test_logs (:warn, Regex("Unknown Sparlectra configuration key diagnostics\\.$(key) is ignored")) match_mode = :any Sparlectra.load_sparlectra_config(cfg_bad; reload = true)
      @test cfg.diagnostics.log_effective_config
      @test !hasproperty(cfg.diagnostics, Symbol(key))
    end
  end)() end
  return nothing
end

function test_configuration_refresh()
  @testset "Configuration refresh" begin (function ()
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

    # form fields of the Web UI benchmark option (taken out in 0.20.2) in an
    # older stored configuration: removed and named, the refresh does not
    # stop; the benchmark section of the configuration stays. Since 0.30.2
    # the rest of the `form` block goes too: the Web UI reads form defaults
    # from case files only, so in the configuration file it is never read
    retired = test_scratch_path(".yaml")
    write(retired, "config_version: 1\nscope: general\nbenchmark:\n  samples: 7\nform:\n  benchmark_samples: 10\n  benchmark_seconds: 1.0\n  gen_seed: 3\n")
    retired_result = Sparlectra.refresh_sparlectra_config_file(retired)
    @test retired_result.success
    @test retired_result.changed
    @test retired_result.removed_keys == ["form", "form.benchmark_samples", "form.benchmark_seconds"]
    @test count(w -> occursin("the Web UI has no benchmark option", w), retired_result.warnings) == 2
    @test !occursin("benchmark_samples", retired_result.refreshed_text)
    @test !occursin("gen_seed", retired_result.refreshed_text)
    @test occursin("samples: 7", retired_result.refreshed_text)

    # keys no loader reads any more are deleted by a refresh and named:
    # a removed key (islands.mode), a silently removed one (webui.warmup) and
    # a key this version does not know, and the form block; the refreshed
    # text loads without an ignored-key warning
    superfluous = test_scratch_path(".yaml")
    write(superfluous, "config_version: 1\nscope: general\npower_flow:\n  tol: 1.0e-7\n  islands:\n    mode: solve_parallel\n  no_such_key: 3\nwebui:\n  warmup: true\nform:\n  gen_seed: 4\n")
    superfluous_result = Sparlectra.refresh_sparlectra_config_file(superfluous)
    @test superfluous_result.success
    @test sort(superfluous_result.removed_keys) == ["form", "power_flow.islands.mode", "power_flow.no_such_key", "webui.warmup"]
    @test count(w -> startswith(w, "Removed power_flow.islands.mode:"), superfluous_result.warnings) == 1
    @test !occursin("solve_parallel", superfluous_result.refreshed_text) && !occursin("no_such_key", superfluous_result.refreshed_text) && !occursin("warmup: true", superfluous_result.refreshed_text)
    @test occursin("tol: 1.0e-7", superfluous_result.refreshed_text) && !occursin("gen_seed", superfluous_result.refreshed_text)
    reloaded = test_scratch_path(".yaml")
    write(reloaded, superfluous_result.refreshed_text)
    reloaded_cfg = @test_logs min_level = Logging.Warn Sparlectra.load_sparlectra_config(reloaded; reload = true)
    @test reloaded_cfg.powerflow.tol == 1.0e-7

    for (legacy, canonical) in (("matpower_simultaneous", "classic_simultaneous"), ("matpower_one_at_a_time", "classic_one_at_a_time"))
      p = test_scratch_path(".yaml")
      write(p, "power_flow:\n  qlimits:\n    enforcement_mode: $(legacy)\n")
      result = Sparlectra.refresh_sparlectra_config_file(p)
      @test "power_flow.qlimits.enforcement_mode" in result.normalized_keys
      @test occursin("enforcement_mode: $(canonical)", result.refreshed_text)
    end

    # regression: refresh on a version-0 file must MOVE the
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
  end)() end
  return nothing
end

function run_legacy_configuration_coverage_tests()
  # the same wrapper as in run_configuration_coverage_tests: the warnings
  # these functions provoke on purpose are expected, any other one is not
  for testfn in (
    test_configuration_qlimit_enforcement_keys,
    test_configuration_value_domain_validation,
    test_configuration_webui_keys_both_directions,
    test_configuration_removed_diagnostics_rejected,
    test_configuration_refresh,
  )
    run_with_expected_warnings(testfn, CONFIG_EXPECTED_WARNINGS)
  end
  return nothing
end

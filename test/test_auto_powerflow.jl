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

# file: test/test_auto_powerflow.jl
# purpose: automatic power-flow mode: feature extraction on the four case
#          classes, ordered-rule decision engine, step-control exclusion,
#          escalation ladder with cap and skip reasons, Q-limit dimension,
#          hint generator ids, precedence, DC fallback, startup hint

function _autopf_features(f) # synthetic AutoPfFeatures with defaults
  d = Dict{Symbol,Any}(
    :n_bus => 100, :n_branch => 120, :n_ac_islands => 1, :rx_median => 0.3,
    :rx_p90 => 0.6, :share_r_gt_x => 0.1, :n_phase_shifters => 0, :n_pv => 20,
    :pv_share => 0.2, :n_q_zero_range => 0, :n_q_narrow_range => 0,
    :n_gen_with_q_limits => 10, :share_pv_with_limits => 1.0,
    :has_start_profile => false, :start_profile_plausible => false,
    :n_voltage_levels => 2,
  )
  merge!(d, Dict{Symbol,Any}(pairs(f)))
  return Sparlectra.AutoPfFeatures((d[k] for k in fieldnames(Sparlectra.AutoPfFeatures))...)
end

# a produced option set must never activate merit and trust-region together
# (transitive exclusion through autodamp is a validation error)
function _autopf_step_control_ok(options::AbstractDict)
  merit = get(options, "power_flow.merit.enabled", false) === true
  trust = get(options, "power_flow.trust_region.enabled", false) === true
  auto = get(options, "power_flow.autodamp", false) === true
  return !(merit && trust) && (!merit || auto) && !(trust && auto)
end

function run_auto_powerflow_tests()
  @testset "auto powerflow" begin
    @testset "feature extraction on the four case classes" begin
      cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
      # class 1: small transmission case (load_fixture_net: the shipped
      # sp_case14 replaces the downloaded case14; 14 buses, machines with
      # Q bands, one synchronous island, all the features this classifier
      # reads)
      net14 = Sparlectra._se_import_case_net(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")), cfg)
      f14 = Sparlectra.collect_auto_pf_features(net14)
      @test f14.n_bus == 14
      @test f14.n_ac_islands == 1
      @test f14.n_gen_with_q_limits > 0
      @test f14.rx_median > 0.0
      @test 0.0 <= f14.share_r_gt_x <= 1.0
      # class 2: phase-shifter case (tracked PST demo case)
      netpst = Sparlectra._se_import_case_net(abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m")), cfg)
      fpst = Sparlectra.collect_auto_pf_features(netpst)
      @test fpst.n_phase_shifters >= 1
      # class 3: resistive network (R > X on every branch; test-local mutation
      # AFTER the read-only extraction contract was already exercised above)
      netres = Sparlectra._se_import_case_net(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")), cfg)
      for br in netres.branchVec
        br.r_pu = 2.0 * abs(br.x_pu) + 0.01
      end
      fres = Sparlectra.collect_auto_pf_features(netres)
      @test fres.share_r_gt_x == 1.0
      @test fres.rx_median > 1.0
      # class 4: multi-island (open enough branches to isolate bus 14)
      netisl = Sparlectra._se_import_case_net(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")), cfg)
      for br in netisl.branchVec
        (Int(br.fromBus) == 14 || Int(br.toBus) == 14) && (br.status = 0)
      end
      fisl = Sparlectra.collect_auto_pf_features(netisl)
      @test fisl.n_ac_islands >= 2
      # extraction is read-only on the untouched net: repeated calls agree
      @test Sparlectra.collect_auto_pf_features(net14) == f14
    end

    @testset "decision rules (first match wins, exclusion by construction)" begin
      small = Sparlectra.select_auto_pf_strategy(_autopf_features((n_bus = 14,)))
      @test small.profile === :small_default
      assisted = Sparlectra.select_auto_pf_strategy(_autopf_features((n_bus = 300, has_start_profile = true, start_profile_plausible = true)))
      @test assisted.profile === :profile_assisted
      @test assisted.options["power_flow.start_mode.voltage_mode"] == "profile_blend"
      trans = Sparlectra.select_auto_pf_strategy(_autopf_features((n_bus = 300, share_r_gt_x = 0.05)))
      @test trans.profile === :transmission_flat
      @test trans.options["power_flow.start_mode.angle_mode"] == "dc"
      res = Sparlectra.select_auto_pf_strategy(_autopf_features((n_bus = 300, share_r_gt_x = 0.7, rx_median = 1.5)))
      @test res.profile === :resistive
      @test res.options["power_flow.start_current_iteration.enabled"] === true
      hard = Sparlectra.select_auto_pf_strategy(_autopf_features((n_bus = 5000,)))
      @test hard.profile === :hard_case
      # first match wins: a plausible profile beats the hard_case trigger
      both = Sparlectra.select_auto_pf_strategy(_autopf_features((n_bus = 5000, has_start_profile = true, start_profile_plausible = true)))
      @test both.profile === :profile_assisted
      # start-value core rule: plausible profile beats flat start everywhere
      @test !haskey(both.options, "power_flow.start_mode.flatstart")
      for d in (small, assisted, trans, res, hard, both)
        @test _autopf_step_control_ok(d.options)
        # the solver-internal rescue ladder stays active (calibration
        # sweep: case13659pegase converges only through it)
        @test !haskey(d.options, "power_flow.rescue")
        @test !isempty(d.reasons)
      end
      # islands gate independent of profile
      isl = Sparlectra.select_auto_pf_strategy(_autopf_features((n_bus = 14, n_ac_islands = 3)))
      @test isl.options["power_flow.islands.enabled"] === true
      # Q-limit dimension: no limits -> subsystem off; zero ranges -> locks
      noq = Sparlectra.select_auto_pf_strategy(_autopf_features((n_bus = 14, n_gen_with_q_limits = 0)))
      @test noq.options["power_flow.qlimits.enabled"] === false
      zq = Sparlectra.select_auto_pf_strategy(_autopf_features((n_bus = 14, n_q_zero_range = 2)))
      @test zq.options["power_flow.qlimits.guard_zero_range_mode"] == "lock_pq"
      @test zq.options["power_flow.qlimits.guard_narrow_range_mode"] == "lock_pq"
    end

    @testset "escalation ladder, cap, precedence, service artifacts" begin
      # load_fixture_net: the shipped sp_case14 replaces the downloaded
      # case14. Its start_state converges in one iteration, so the forced
      # non-convergence adds the user flatstart override (which the auto
      # mode must respect anyway); warmup_casePST is unusable here because
      # its hard_case profile hits the ladder stage-combination finding
      # recorded in the report (merit left on while autodamp is switched
      # off)
      case = abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"))
      # forced non-convergence (user-set max_iter survives every stage) runs
      # the ladder to its cap with honest skip reasons
      cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true, overrides = Dict{String,Any}("power_flow" => Dict{String,Any}("mode" => "auto", "max_iter" => 1, "flatstart" => true)))
      r = run_sparlectra(casefile = case, config = cfg)
      @test !r.final_converged
      rec = Sparlectra.auto_pf_record(r.net)
      @test rec !== nothing
      stages = [a.stage for a in rec.attempts]
      @test stages[1] === :base
      @test :L1_step_control in stages
      @test count(a -> !isempty(a.solver) && a.stage !== :base, rec.attempts) <= Sparlectra.AUTO_PF_THRESHOLDS.escalation_max_stages
      l1 = rec.attempts[findfirst(a -> a.stage === :L1_step_control, rec.attempts)]
      @test any(c -> occursin("merit", c) || occursin("trust_region", c) || occursin("autodamp", c), l1.changed)
      l5 = rec.attempts[findfirst(a -> a.stage === :L5_qlimit_mode, rec.attempts)]
      @test isempty(l5.solver) && l5.qlimit_evidence.skipped == "skipped_no_qlimit_evidence"
      # tolerance is untouched by every stage
      @test cfg.powerflow.tol == Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true).powerflow.tol
      # precedence: the user's explicit key survives and the conflict is logged
      cfg2 = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true, overrides = Dict{String,Any}("power_flow" => Dict{String,Any}("mode" => "auto", "autodamp" => false)))
      r2 = run_sparlectra(casefile = case, config = cfg2)
      rec2 = Sparlectra.auto_pf_record(r2.net)
      @test any(c -> startswith(c, "power_flow.autodamp"), rec2.conflicts)
      # service path: metadata mirror plus the decision-log artifact
      root = mktempdir()
      rs = start_powerflow_run(Dict{String,Any}("casefile" => case, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "config_overrides" => Dict{String,Any}("power_flow.mode" => "auto")))
      @test rs["status"] == "succeeded"
      @test rs["metadata"]["auto_mode_enabled"] === true
      @test haskey(rs["metadata"], "auto_profile")
      @test haskey(rs["metadata"], "auto_final_solver")
      logp = joinpath(rs["output_dir"], "auto_mode_decision.log")
      @test isfile(logp)
      logtxt = read(logp, String)
      @test occursin("selected profile", logtxt)
      @test occursin("attempts:", logtxt)
      # conflict text reaches the artifact too
      rs2 = start_powerflow_run(Dict{String,Any}("casefile" => case, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "config_overrides" => Dict{String,Any}("power_flow.mode" => "auto", "power_flow.autodamp" => false)))
      @test occursin("kept at the user value", read(joinpath(rs2["output_dir"], "auto_mode_decision.log"), String))
      # manual runs only carry the disabled flag
      rman = start_powerflow_run(Dict{String,Any}("casefile" => case, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root))
      @test rman["metadata"]["auto_mode_enabled"] === false
      @test !haskey(rman["metadata"], "auto_profile")
    end

    @testset "DC fallback: honest labeling, off by default" begin
      # the service path cannot override power_flow.flatstart (not a GUI
      # key), so forcing real flat-start non-convergence at max_iter=1
      # needs a start_state-free copy of the shipped case
      case = joinpath(mktempdir(), "sp_case14_flat.scf.json")
      raw = Sparlectra.scf_json_parse(read(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")), String))
      delete!(raw["sparlectra"], "start_state")
      write(case, Sparlectra.scf_json_string(raw))
      root = mktempdir()
      # enabled: the ladder ends in the labeled DC approximation
      rdc = start_powerflow_run(Dict{String,Any}("casefile" => case, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "config_overrides" => Dict{String,Any}("power_flow.mode" => "auto", "power_flow.max_iter" => 1, "power_flow.dc.fallback" => true)))
      @test rdc["status"] != "succeeded"
      @test rdc["converged"] === false
      @test rdc["metadata"]["dc_fallback_solution"] === true
      @test occursin("DC fallback", rdc["message"])
      # default: clean failure without a DC result
      rnf = start_powerflow_run(Dict{String,Any}("casefile" => case, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "config_overrides" => Dict{String,Any}("power_flow.mode" => "auto", "power_flow.max_iter" => 1)))
      @test rnf["metadata"]["dc_fallback_solution"] === false
      @test occursin("skipped_dc_fallback_disabled", read(joinpath(rnf["output_dir"], "auto_mode_decision.log"), String))
    end

    @testset "stage builders: L4 locks, L5 gate, L6 seed-only" begin
      pf = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true).powerflow
      stages = Sparlectra._auto_pf_escalation_stages(pf)
      byid = Dict(s.id => s for s in stages)
      ev_chatter = (switching_events = 30, active_set_changes = 5, reenable_events = 0, switch_counts = Dict(5 => pf.qlimits.guard_max_switches + 2, 7 => 1))
      pf4, changed4 = byid[:L4_qlimit_relax].transform(pf, ev_chatter)
      @test 5 in pf4.qlimits.lock_pv_to_pq_buses
      @test !(7 in pf4.qlimits.lock_pv_to_pq_buses)
      @test pf4.qlimits.guard_freeze_after_repeated_switching
      @test any(c -> occursin("lock_pv_to_pq_buses", c), changed4)
      # L5 gate: evidence runs it, no evidence skips it with a stable reason
      ev_none = (switching_events = 0, active_set_changes = 0, reenable_events = 0, switch_counts = Dict{Int,Int}())
      @test byid[:L5_qlimit_mode].gate(ev_chatter, pf) == (true, "")
      @test byid[:L5_qlimit_mode].gate(ev_none, pf)[2] == "skipped_no_qlimit_evidence"
      pf5, _ = byid[:L5_qlimit_mode].transform(pf, ev_chatter)
      @test pf5.qlimits.enforcement_mode !== pf.qlimits.enforcement_mode
      # L6: APSLF is a seed only; the transform enables apslf_start and the
      # solver stays rectangular (label test). The stage used to be gated on
      # the AnalyticLoadFlow extension being loaded; the package is a required
      # dependency now, so the gate is always open
      @test byid[:L6_apslf_seed].gate(ev_none, pf) == (true, "")
      pf6, _ = byid[:L6_apslf_seed].transform(pf, ev_none)
      @test pf6.apslf_start.enabled
      @test pf6.solver === :rectangular
      @test occursin("apslf-seeded", Sparlectra._auto_pf_solver_label(pf6, nothing))
    end

    @testset "hint generator: stable ids from synthetic diagnostics" begin
      ev(sc; se = 0, asc = 0) = (switching_events = se, active_set_changes = asc, reenable_events = 0, switch_counts = sc)
      f_no = _autopf_features((n_gen_with_q_limits = 0,))
      ids(h) = Set(x.id for x in h)
      h1 = Sparlectra.auto_pf_hints(f_no, ev(Dict{Int,Int}()), NamedTuple[]; converged = true)
      @test :q_no_limits_defined in ids(h1)
      f_zero = _autopf_features((n_q_zero_range = 2, n_q_narrow_range = 1))
      h2 = Sparlectra.auto_pf_hints(f_zero, ev(Dict{Int,Int}()), NamedTuple[]; converged = true)
      @test :q_zero_range_locked in ids(h2)
      @test :q_narrow_range in ids(h2)
      h3 = Sparlectra.auto_pf_hints(_autopf_features(NamedTuple()), ev(Dict(3 => 12)), NamedTuple[]; converged = true)
      @test :q_switching_frozen in ids(h3)
      attempts = NamedTuple[(stage = :L5_qlimit_mode, converged = true)]
      h4 = Sparlectra.auto_pf_hints(_autopf_features(NamedTuple()), ev(Dict{Int,Int}()), attempts; converged = true)
      @test :q_mode_pinned in ids(h4)
      h5 = Sparlectra.auto_pf_hints(_autopf_features(NamedTuple()), ev(Dict(4 => 2); se = 9), NamedTuple[]; converged = false)
      @test :q_failure_with_switching in ids(h5)
      for h in (h1..., h2..., h3..., h4..., h5...)
        @test h.id isa Symbol
        @test !isempty(h.text)
      end
    end

    @testset "startup latency hint: once, suppressible, flavor-gated" begin
      cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
      if Sparlectra.webui_runtime_flavor().kind === :native
        Sparlectra._STARTUP_LATENCY_HINT_SHOWN[] = false
        @test_logs (:info, r"First run in a fresh Julia process") Sparlectra._maybe_print_startup_latency_hint(cfg)
        # second call in the same process stays silent
        @test_logs Sparlectra._maybe_print_startup_latency_hint(cfg)
        # suppression key silences the first call too
        Sparlectra._STARTUP_LATENCY_HINT_SHOWN[] = false
        cfg_off = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true, overrides = Dict{String,Any}("output" => Dict{String,Any}("startup_latency_hint" => false)))
        @test_logs Sparlectra._maybe_print_startup_latency_hint(cfg_off)
      else
        println("      auto_powerflow: startup-hint print path SKIPPED (session runs from sysimage/app); covering the gate only")
        Sparlectra._STARTUP_LATENCY_HINT_SHOWN[] = false
        @test_logs Sparlectra._maybe_print_startup_latency_hint(cfg)
      end
      Sparlectra._STARTUP_LATENCY_HINT_SHOWN[] = true
    end
  end

  @testset "auto profile precedence level (D11)" begin
    base = Sparlectra.SparlectraConfig(Dict{String,Any}())
    pairs = Pair{Symbol,Any}[:ratio => :reciprocal, :shift_unit => :rad]
    # nothing explicitly set: both recommendations apply
    lev = Sparlectra.apply_auto_profile_level(base, pairs)
    @test lev.config.matpower.ratio === :reciprocal
    @test lev.config.matpower.shift_unit === :rad
    @test isempty(lev.skipped)
    # a key set by the user to a NON-DEFAULT value survives; the rest applies
    user = Sparlectra._copy_sparlectra_with_user_keys(Sparlectra.SparlectraConfig(Dict{String,Any}("matpower" => Dict{String,Any}("shift_sign" => -1.0))), Set(["matpower_import.shift_sign"]))
    pairs2 = Pair{Symbol,Any}[:shift_sign => 1.0, :ratio => :reciprocal]
    lev2 = Sparlectra.apply_auto_profile_level(user, pairs2)
    @test lev2.config.matpower.shift_sign == -1.0
    @test lev2.config.matpower.ratio === :reciprocal
    @test length(lev2.skipped) == 1
    @test first(first(lev2.skipped)) === :shift_sign
    # a key merely PRESENT at its default (refresh-written YAML) is no
    # choice: the recommendation still applies
    present = Sparlectra._copy_sparlectra_with_user_keys(base, Set(["matpower_import.ratio"]))
    lev3 = Sparlectra.apply_auto_profile_level(present, pairs)
    @test lev3.config.matpower.ratio === :reciprocal
    @test isempty(lev3.skipped)
    # unknown field: hard error instead of a silent rewrite
    @test_throws ArgumentError Sparlectra.apply_auto_profile_level(base, Pair{Symbol,Any}[:tol => 1.0e-9])
    # resolve_config carries the level; an explicit override wins over it
    mktempdir() do d
      case = joinpath(d, "caseX.m")
      write(case, "% fixture\n")
      r = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, case; auto_profile_overrides = Dict{String,Any}("matpower_import.ratio" => "reciprocal"))
      @test r.config.matpower.ratio === :reciprocal
      @test haskey(r.auto_profile_config, "matpower_import.ratio")
      r2 = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, case, Dict{String,Any}("matpower_import.ratio" => "normal"); auto_profile_overrides = Dict{String,Any}("matpower_import.ratio" => "reciprocal"))
      @test r2.config.matpower.ratio === :normal
      @test !haskey(r2.auto_profile_config, "matpower_import.ratio")
    end
  end
  return nothing
end

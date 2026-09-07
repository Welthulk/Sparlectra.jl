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

# file: test/test_scenarios.jl
# purpose: scenario task step 1: the patch model, its load-time validation
#          (rejections carry the scenario name and the op index), the
#          sparlectra.scenarios file round trip, the legacy contingencies
#          mapping, and the N-1 expansion equality with the existing
#          generators, on the tracked MATPOWER warmup case
#          (load_fixture_net: the downloaded case14 is gone; the length
#          assertions are all RELATIVE to the case, so only the fixture
#          changed).

function _scenario_case14()
  cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
  mpc = Sparlectra.MatpowerIO.read_case(abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m")); legacy_compat = true)
  case = Sparlectra.convert_case(MatpowerAdapter(), mpc, Sparlectra.matpower_adapter_options(cfg))
  net = Sparlectra.build_net(case; config = cfg)
  return case, net
end

function run_scenario_patch_tests()
  @testset "scenario patch model" begin
    case, net = _scenario_case14()
    index = ScenarioIndex(case)
    branch_id = first(r.id for r in case.data.line)
    trafo_id = first(r.id for r in case.data.generic_branch)
    gen_id = first(r.id for r in case.data.sym_gen)
    load_id = first(r.id for r in case.data.sym_load)

    @testset "validation accepts the D1 surface" begin
      good = ScenarioSet(scenarios = [
        Scenario(name = "line out", ops = [PatchOp(op = :status, target = :branch, id = branch_id, value = 0.0)]),
        Scenario(name = "gen setpoint", ops = [PatchOp(op = :set, target = :generator, id = gen_id, field = :p, value = 25.0)]),
        Scenario(name = "load scale", weight = 2.0, ops = [PatchOp(op = :scale, target = :load, id = load_id, factor = 1.2)]),
        Scenario(name = "tap step", ops = [PatchOp(op = :set, target = :transformer, id = trafo_id, field = :tap_pos, value = 2.0)]),
      ])
      @test validate_scenarios(good, index) === good
    end

    @testset "rejections name the scenario and the op index" begin
      unknown = ScenarioSet(scenarios = [Scenario(name = "ghost", ops = [PatchOp(op = :status, target = :branch, id = 999999, value = 0.0)])])
      err = try
        validate_scenarios(unknown, index)
        ""
      catch e
        sprint(showerror, e)
      end
      @test occursin("\"ghost\"", err)
      @test occursin("op 1", err)
      @test occursin("unknown component id", err)

      wrongfield = ScenarioSet(scenarios = [Scenario(name = "bad field", ops = [PatchOp(op = :set, target = :load, id = load_id, field = :tap_pos, value = 1.0)])])
      err2 = try
        validate_scenarios(wrongfield, index)
        ""
      catch e
        sprint(showerror, e)
      end
      @test occursin("\"bad field\"", err2)
      @test occursin("tap_pos applies to a transformer", err2)

      empty_ops = ScenarioSet(scenarios = [Scenario(name = "hollow", ops = PatchOp[])])
      @test_throws ArgumentError validate_scenarios(empty_ops, index)

      dup = ScenarioSet(scenarios = [
        Scenario(name = "twin", ops = [PatchOp(op = :status, target = :branch, id = branch_id, value = 0.0)]),
        Scenario(name = "twin", ops = [PatchOp(op = :status, target = :branch, id = branch_id, value = 1.0)]),
      ])
      @test_throws ArgumentError validate_scenarios(dup, index)

      mismatch = ScenarioSet(scenarios = [Scenario(name = "class clash", ops = [PatchOp(op = :status, target = :generator, id = branch_id, value = 0.0)])])
      err3 = try
        validate_scenarios(mismatch, index)
        ""
      catch e
        sprint(showerror, e)
      end
      @test occursin("is a branch, not a generator", err3)
    end

    @testset "file round trip of the scenarios block" begin
      set = ScenarioSet(
        scenarios = [
          Scenario(name = "double outage", weight = 0.5, ops = [
            PatchOp(op = :status, target = :branch, id = branch_id, value = 0.0),
            PatchOp(op = :status, target = :transformer, id = trafo_id, value = 0.0),
          ]),
          Scenario(name = "load up", ops = [PatchOp(op = :scale, target = :load, id = load_id, field = :q, factor = 1.1)]),
        ],
        mode = :explicit,
        exclusions = ["kept out"],
      )
      case.sparlectra.scenarios = scenario_set_dict(set)
      path = joinpath(mktempdir(), "scenario_case.scf.json")
      write_scf_json(case, path)
      back = read_scf_json(path)
      got = scf_case_scenarios(back)
      @test got !== nothing
      @test got.mode === :explicit
      @test got.exclusions == ["kept out"]
      @test length(got.scenarios) == 2
      @test got.scenarios[1].name == "double outage"
      @test got.scenarios[1].weight == 0.5
      @test length(got.scenarios[1].ops) == 2
      @test got.scenarios[1].ops[2].target === :transformer
      @test got.scenarios[2].ops[1].factor == 1.1
      @test got.scenarios[2].ops[1].field === :q
      @test validate_scenarios(got, ScenarioIndex(back)) === got
      case.sparlectra.scenarios = Dict{String,Any}()
    end

    @testset "legacy contingencies map onto the scenario model" begin
      raw = Dict{String,Any}(
        "mode" => "explicit",
        "cases" => Any[Dict{String,Any}("name" => "old style", "outages" => Any[Dict{String,Any}("component" => branch_id)])],
      )
      mapped = Sparlectra.scenario_set_from_contingencies(raw, case)
      @test length(mapped.scenarios) == 1
      @test mapped.scenarios[1].name == "old style"
      @test mapped.scenarios[1].ops[1].op === :status
      @test mapped.scenarios[1].ops[1].id == branch_id
      @test mapped.scenarios[1].ops[1].target === :branch
      @test validate_scenarios(mapped, index) === mapped
    end

    @testset "apply and restore: tap patch semantics" begin
      # a tap patch on an UNREGULATED transformer changes the ratio the
      # solver sees and restores it (maintainer decision 1: first-class)
      tapped = findfirst(br -> br.has_ratio_tap, net.branchVec)
      if tapped === nothing
        println("      scenarios: tap-patch semantics SKIPPED (warmup case build carries no ratio tap changer)")
      else
        tap_id = first(id for (id, k) in index.kind_by_id if k === :transformer && get(index.internal_by_id, id, 0) == tapped)
        br = net.branchVec[tapped]
        before_ratio = Sparlectra.calcBranchRatio(br)
        undo = apply!(net, [PatchOp(op = :set, target = :transformer, id = tap_id, field = :tap_pos, value = 2.0)], index)
        @test Sparlectra.calcBranchRatio(net.branchVec[tapped]) != before_ratio
        restore!(net, undo)
        @test Sparlectra.calcBranchRatio(net.branchVec[tapped]) == before_ratio
      end
    end

    @testset "tap patch on a REGULATED transformer is rejected, controller named" begin
      tapped = findfirst(br -> br.has_ratio_tap, net.branchVec)
      if tapped === nothing
        println("      scenarios: regulated-tap rejection SKIPPED (warmup case build carries no ratio tap changer)")
      else
        tap_id = first(id for (id, k) in index.kind_by_id if k === :transformer && get(index.internal_by_id, id, 0) == tapped)
        trafo_k = count(i -> Sparlectra._scf_is_transformer(net.branchVec[i]), 1:tapped)
        ctrl = Sparlectra.PowerTransformerControl(trafo = "T_reg", mode = :voltage, target_bus = Sparlectra._scf_bus_name(net, Int(net.branchVec[tapped].toBus)), target_vm_pu = 1.0)
        push!(net.trafos[trafo_k].side1.controls, ctrl)
        regulated = ScenarioSet(scenarios = [Scenario(name = "fights the controller", ops = [PatchOp(op = :set, target = :transformer, id = tap_id, field = :tap_pos, value = 1.0)])])
        err = try
          validate_scenarios(regulated, index, net)
          ""
        catch e
          sprint(showerror, e)
        end
        @test occursin("regulated by controller", err)
        @test occursin("T_reg", err)
        pop!(net.trafos[trafo_k].side1.controls)
        @test validate_scenarios(regulated, index, net) === regulated
      end
    end

    @testset "bitwise restore over every N-1 scenario" begin
      # THE step-2 acceptance test: for every N-1 branch and generator
      # scenario, apply! then restore! leaves the working copy bitwise
      # equal to the base on every field of every component type, checked
      # with the reflection comparison the SCF round trip uses
      # load_fixture_net: the MATPOWER leg runs on the tracked warmup case,
      # the breadth (couplers, FACTS, controllers, HVDC) comes from the
      # shipped demo cases; nothing downloads on a fresh install
      for casefile in ("warmup_casePST.m",)
        cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
        mpc = Sparlectra.MatpowerIO.read_case(abspath(joinpath(dirname(@__DIR__), "data", "mpower", casefile)); legacy_compat = true)
        scfcase = Sparlectra.convert_case(MatpowerAdapter(), mpc, Sparlectra.matpower_adapter_options(cfg))
        base = Sparlectra.build_net(scfcase; config = cfg)
        work = Sparlectra.build_net(scfcase; config = cfg)
        idx = ScenarioIndex(scfcase)
        scenarios = expand_scenarios(ScenarioSet(mode = :n1_all), work, idx)
        for s in scenarios
          undo = apply!(work, s.ops, idx)
          restore!(work, undo)
        end
        diffs = scf_roundtrip_field_diffs(base, work)
        @test (casefile, diffs) == (casefile, String[])
      end
      for demo in ("sp_case60", "sp_case188")
        scfcase = read_scf_json(abspath(joinpath(dirname(@__DIR__), "data", "scf", demo * ".scf.json")))
        base = Sparlectra.build_net(scfcase)
        work = Sparlectra.build_net(scfcase)
        idx = ScenarioIndex(scfcase)
        for s in expand_scenarios(ScenarioSet(mode = :n1_all), work, idx)
          undo = apply!(work, s.ops, idx)
          restore!(work, undo)
        end
        @test (demo, scf_roundtrip_field_diffs(base, work)) == (demo, String[])
      end
      let fixture = abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_casePST.scf.json"))
        scfcase = read_scf_json(fixture)
        base = Sparlectra.build_net(scfcase)
        work = Sparlectra.build_net(scfcase)
        idx = ScenarioIndex(scfcase)
        scenarios = expand_scenarios(ScenarioSet(mode = :n1_all), work, idx)
        @test !isempty(scenarios)
        for s in scenarios
          undo = apply!(work, s.ops, idx)
          restore!(work, undo)
        end
        @test ("warmup_casePST", scf_roundtrip_field_diffs(base, work)) == ("warmup_casePST", String[])
      end
    end

    @testset "engine: worker resets bitwise between evaluations" begin
      # step-3 acceptance for the REUSED working copy: after every
      # evaluation (branch outage, generator outage, patch scenario, all
      # WITH solves) the worker equals the engine template bitwise on
      # every field the SCF reflection comparison checks; state bleeding
      # from one scenario into the next would surface exactly here
      engine = Sparlectra.ScenarioEngine(net; index = index)
      worker = Sparlectra.ScenarioWorker(deepcopy(engine.template))
      branch_internal = first(v for (id, v) in index.internal_by_id if index.kind_by_id[id] === :branch)
      gen_id = first(id for (id, k) in index.kind_by_id if k === :generator)
      load_id2 = first(id for (id, k) in index.kind_by_id if k === :load)
      items = Any[
        Sparlectra._ScenarioOutageItem("branch out", 1.0, :branch, branch_internal),
        Sparlectra._ScenarioOutageItem("gen out", 1.0, :gen, index.internal_by_id[gen_id]),
        Sparlectra._ScenarioPatchItem(Scenario(name = "load up", ops = [PatchOp(op = :scale, target = :load, id = load_id2, factor = 1.3)])),
      ]
      for it in items
        r = Sparlectra.evaluate!(engine, worker, it)
        @test r isa Sparlectra.ContingencyResult
        @test (it, scf_roundtrip_field_diffs(engine.template, worker.net)) == (it, String[])
      end
    end

    @testset "engine: runScenarios! N-1 equals runContingencies!" begin
      cases14 = vcat(Sparlectra.generateN1Branches(net), Sparlectra.generateN1Generators(net))
      via_cases = Sparlectra.runContingencies!(net, cases14)
      via_scenarios = runScenarios!(net, ScenarioSet(mode = :n1_all); index = index)
      d = mktempdir()
      a = Sparlectra.writeContingencyResultsCSV(joinpath(d, "cases.csv"), via_cases)
      b = Sparlectra.writeContingencyResultsCSV(joinpath(d, "scenarios.csv"), via_scenarios)
      @test read(a, String) == read(b, String)
      # a general patch scenario runs the same ladder and metrics
      gen_id = first(id for (id, k) in index.kind_by_id if k === :generator)
      patched = runScenarios!(net, [Scenario(name = "setpoint", ops = [PatchOp(op = :set, target = :generator, id = gen_id, field = :p, value = 30.0)])]; index = index)
      @test length(patched) == 1
      @test patched[1].converged
      @test patched[1].error === nothing
    end

    @testset "screening surface (step 4)" begin
      cases14 = vcat(Sparlectra.generateN1Branches(net), Sparlectra.generateN1Generators(net))
      # :off keeps the historical result type and rows bit-identical
      off = Sparlectra.runContingencies!(net, cases14)
      @test off isa Vector{Sparlectra.ContingencyResult}
      off2 = Sparlectra.runContingencies!(net, cases14; screening_mode = :off)
      d = mktempdir()
      @test read(Sparlectra.writeContingencyResultsCSV(joinpath(d, "a.csv"), off), String) == read(Sparlectra.writeContingencyResultsCSV(joinpath(d, "b.csv"), off2), String)
      # :flag returns the ScenarioResult surface, forwards the inner fields,
      # and never screens away a case the full run reports as violating
      flagged = Sparlectra.runContingencies!(net, cases14; screening_mode = :flag, screening_margin_pct = 10.0)
      @test flagged isa Vector{ScenarioResult}
      @test [r.name for r in flagged] == [r.name for r in off]
      violating = Set(r.name for r in off if !isempty(r.overloads) || !isempty(r.voltage_violations) || !r.converged)
      @test isempty(intersect(violating, Set(r.name for r in flagged if r.screened)))
      for r in flagged
        if r.screened
          @test r.start_used === :screen
          @test r.screening_estimate !== nothing
          @test r.error === nothing
        else
          # a flagged or unscreenable case carries the FULL run row
          @test r.result isa Sparlectra.ContingencyResult
          @test r.start_used !== :screen
        end
      end
      # the screening CSV appends exactly the two D8 columns
      csv = read(Sparlectra.writeContingencyResultsCSV(joinpath(d, "flag.csv"), flagged), String)
      header = first(split(csv, '\n'))
      @test endswith(header, ";screened;screening_estimate")
      # :only never runs a full solve where an estimate exists; islanding
      # cases still carry real solves (start_used != :screen)
      only_res = Sparlectra.runContingencies!(net, cases14; screening_mode = :only, screening_margin_pct = 10.0)
      @test any(r -> r.screened, only_res)
      @test all(r -> r.screened || r.start_used !== :screen, only_res)
    end

    @testset "N-1 expansion equals the existing generators on case14" begin
      expanded = expand_scenarios(ScenarioSet(mode = :n1_branches), net, index)
      reference = Sparlectra.generateN1Branches(net)
      @test length(expanded) == length(reference)
      @test [s.name for s in expanded] == [c.name for c in reference]
      @test all(length(s.ops) == 1 && s.ops[1].op === :status && s.ops[1].value == 0.0 for s in expanded)
      gens = expand_scenarios(ScenarioSet(mode = :n1_generators), net, index)
      gref = Sparlectra.generateN1Generators(net)
      @test length(gens) == length(gref)
      @test [s.name for s in gens] == [c.name for c in gref]
      alln1 = expand_scenarios(ScenarioSet(mode = :n1_all), net, index)
      @test length(alln1) == length(reference) + length(gref)
      excluded = expand_scenarios(ScenarioSet(mode = :n1_branches, exclusions = [reference[1].element]), net, index)
      @test length(excluded) == length(reference) - 1
    end
  end
  return nothing
end

function run_scenario_engine_extended_tests()
  @testset "scenario engine case118 fixture" begin
    # THE step-3 gate: the engine's N-1 result CSV on case118 (all
    # branches, all generators) is byte identical to the CSV the
    # pre-engine per-case-deepcopy implementation produced at 7438b6e
    # (test/fixtures/contingency_case118_n1_7438b6e.csv, generated once
    # from a worktree of that commit). The net is built with the pinned
    # packaged-default import options the fixture generation used.
    fixture = abspath(joinpath(@__DIR__, "fixtures", "contingency_case118_n1_7438b6e.csv"))
    case_path = large_case_path("case118.m")
    if case_path === nothing
      println("      scenario engine: case118 CSV gate SKIPPED (case118.m not in the large-case directory)")
    else
      @test isfile(fixture)
      net = Sparlectra.createNetFromMatPowerFile(filename = case_path, flatstart = false, enable_pq_gen_controllers = true, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)
      cases = vcat(Sparlectra.generateN1Branches(net), Sparlectra.generateN1Generators(net))
      @test length(cases) == 240
      results = Sparlectra.runContingencies!(net, cases)
      out = Sparlectra.writeContingencyResultsCSV(joinpath(mktempdir(), "case118_n1.csv"), results)
      @test read(out, String) == read(fixture, String)

      # step-4 acceptance: with screening :flag and margin 10, every case the
      # full run reports as a limit violation (or failure) gets the full run,
      # no false negatives; the false-positive count and the full-run share
      # are REPORTED, not gated (case118 is a stressed system)
      flagged = Sparlectra.runContingencies!(net, cases; screening_mode = :flag, screening_margin_pct = 10.0)
      @test flagged isa Vector{ScenarioResult}
      violating = Set(r.name for r in results if !isempty(r.overloads) || !isempty(r.voltage_violations) || !r.converged)
      false_negatives = intersect(violating, Set(r.name for r in flagged if r.screened))
      @test (sort(collect(false_negatives)), length(violating) > 0) == (String[], true)
      full_runs = count(r -> !r.screened, flagged)
      fp = count(r -> !r.screened && !(r.name in violating), flagged)
      println("      scenario engine case118 screening: ", length(flagged) - full_runs, " screened, ", full_runs, " full runs (", round(100 * full_runs / length(flagged); digits = 1), " %), ", fp, " false positives, 0 false negatives")

      # step 4b: the same acceptance with DISTRIBUTED SLACK enabled; the
      # screening runs on the augmented system (lambda state, renormalized
      # participation on a generator outage) and must not screen away any
      # case the distributed-slack full run reports as violating or failed
      ds_kw = (distributed_slack_enabled = true,)
      ds_full = Sparlectra.runContingencies!(net, cases; ds_kw...)
      t_ds = @timed Sparlectra.runContingencies!(net, cases; screening_mode = :flag, screening_margin_pct = 10.0, ds_kw...)
      ds_flagged = t_ds.value
      ds_violating = Set(r.name for r in ds_full if !isempty(r.overloads) || !isempty(r.voltage_violations) || !r.converged)
      ds_fneg = intersect(ds_violating, Set(r.name for r in ds_flagged if r.screened))
      @test (sort(collect(ds_fneg)), length(ds_flagged)) == (String[], length(cases))
      ds_screened = count(r -> r.screened, ds_flagged)
      @test ds_screened > 0
      println("      scenario engine case118 distributed-slack screening: ", ds_screened, " of ", length(cases), " screened (", round(100 * ds_screened / length(cases); digits = 1), " %), :flag ", round(t_ds.time; digits = 2), " s, 0 false negatives")
    end
  end
  @testset "scenario engine sp_case60 screening acceptance" begin
    # load_fixture_net: the shipped operated grid judges screening QUALITY
    # on every install (the case300 anchor below stays cache-gated with its
    # historical trust-gate story): no false negatives at margin 10, with
    # and without distributed slack
    net = Sparlectra.importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case60.scf.json")))
    cases = vcat(Sparlectra.generateN1Branches(net), Sparlectra.generateN1Generators(net))
    full = Sparlectra.runContingencies!(net, cases)
    flagged = Sparlectra.runContingencies!(net, cases; screening_mode = :flag, screening_margin_pct = 10.0)
    violating = Set(r.name for r in full if !isempty(r.overloads) || !isempty(r.voltage_violations) || !r.converged)
    @test sort(collect(intersect(violating, Set(r.name for r in flagged if r.screened)))) == String[]
    @test count(r -> r.screened, flagged) > 0
    ds_full = Sparlectra.runContingencies!(net, cases; distributed_slack_enabled = true)
    ds_flagged = Sparlectra.runContingencies!(net, cases; screening_mode = :flag, screening_margin_pct = 10.0, distributed_slack_enabled = true)
    ds_violating = Set(r.name for r in ds_full if !isempty(r.overloads) || !isempty(r.voltage_violations) || !r.converged)
    @test sort(collect(intersect(ds_violating, Set(r.name for r in ds_flagged if r.screened)))) == String[]
    @test count(r -> r.screened, ds_flagged) > 0
    println("      scenario engine sp_case60: RAN (", length(cases), " cases, ", count(r -> r.screened, flagged), " screened, 0 false negatives, ds ", count(r -> r.screened, ds_flagged), " screened)")
  end

  @testset "scenario engine sp_case60 CSV byte fixture" begin
    # load_fixture_net: the engine's N-1 result CSV on the shipped operated
    # grid is byte identical to the tracked fixture generated at bf0b076
    # (the place-name demo cases). Unlike the case118 gate above (an
    # INDEPENDENT oracle from the pre-engine implementation, which stays
    # cache-gated), this fixture guards regressions of the engine's values
    # and CSV format on every install without a download; sp_case60 keeps
    # the tracked file small (10k, the 188 variant carried 27-bus violation
    # lists per row).
    fixture = abspath(joinpath(@__DIR__, "fixtures", "contingency_sp_case60_n1_bf0b076.csv"))
    @test isfile(fixture)
    net = Sparlectra.importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case60.scf.json")))
    cases = vcat(Sparlectra.generateN1Branches(net), Sparlectra.generateN1Generators(net))
    results = Sparlectra.runContingencies!(net, cases)
    out = Sparlectra.writeContingencyResultsCSV(joinpath(mktempdir(), "sp60_n1.csv"), results)
    @test read(out, String) == read(fixture, String)
  end

  @testset "scenario service sources (step 5)" begin
    # the service request picks the scenario source and screening mode: an
    # SCF case with its own scenarios block runs via file_block, an external
    # scenario JSON via external_file, and the two D8 columns plus the
    # screening metadata arrive exactly when screening is active
    cfgpath = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH
    cfg = Sparlectra.load_sparlectra_config(cfgpath; reload = true)
    mpc = Sparlectra.MatpowerIO.read_case(abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m")); legacy_compat = true)
    case = Sparlectra.convert_case(MatpowerAdapter(), mpc, Sparlectra.matpower_adapter_options(cfg))
    index = ScenarioIndex(case)
    branch_id = first(r.id for r in case.data.line)
    gen_id = first(id for (id, k) in index.kind_by_id if k === :generator)
    set = ScenarioSet(scenarios = [
      Scenario(name = "line out", ops = [PatchOp(op = :status, target = :branch, id = branch_id, value = 0.0)]),
      Scenario(name = "gen setpoint", ops = [PatchOp(op = :set, target = :generator, id = gen_id, field = :p, value = 10.0)]),
    ])
    case.sparlectra.scenarios = scenario_set_dict(set)
    dir = mktempdir()
    case_file = joinpath(dir, "step5_case.scf.json")
    write_scf_json(case, case_file)
    redirect_stdout(devnull) do
      res = Sparlectra._run_contingency_service(case_file, cfgpath, joinpath(dir, "run_fb"), "step5_fb", "branch"; scenario_source = "file_block", screening_mode = "flag")
      d1 = Sparlectra.to_dict(res)
      @test d1["status"] == "succeeded"
      @test d1["metadata"]["contingency_screening_mode"] == "flag"
      @test d1["metadata"]["contingency_cases_source"] == "scenario_file_block"
      @test haskey(d1["metadata"], "contingency_screened")
      csv1 = readlines(joinpath(dir, "run_fb", "contingency_n1.csv"))
      @test endswith(first(csv1), ";screened;screening_estimate")
      @test length(csv1) == 3
      # external scenario JSON, screening explicitly off: classic CSV columns
      scen_file = joinpath(dir, "scenarios.json")
      write(scen_file, Sparlectra.scf_json_string(scenario_set_dict(set)))
      res2 = Sparlectra._run_contingency_service(case_file, cfgpath, joinpath(dir, "run_ext"), "step5_ext", "branch"; scenario_source = "external_file", scenario_file = scen_file, screening_mode = "off")
      d2 = Sparlectra.to_dict(res2)
      @test d2["status"] == "succeeded"
      @test d2["metadata"]["contingency_screening_mode"] == "off"
      @test d2["metadata"]["contingency_cases_source"] == "scenario_external_file"
      csv2 = readlines(joinpath(dir, "run_ext", "contingency_n1.csv"))
      @test !endswith(first(csv2), ";screened;screening_estimate")
      # an n1 source works on any format and records itself
      res3 = Sparlectra._run_contingency_service(case_file, cfgpath, joinpath(dir, "run_n1"), "step5_n1", "branch"; scenario_source = "n1_generators", screening_mode = "off")
      d3 = Sparlectra.to_dict(res3)
      @test d3["status"] == "succeeded"
      @test d3["metadata"]["contingency_cases_source"] == "n1_generators"
      # rejections carry invalid_request
      bad = Sparlectra.to_dict(Sparlectra._run_contingency_service(case_file, cfgpath, joinpath(dir, "run_bad"), "step5_bad", "branch"; scenario_source = "bogus"))
      @test bad["reason"] == "invalid_request"
      # a CGMES delivery with an id-addressed scenario source is rejected
      # WITH the way out named (export as SCF once, then reference its ids)
      cgdir = mktempdir()
      cg = Sparlectra.to_dict(Sparlectra._run_contingency_service(cgdir, cfgpath, joinpath(dir, "run_cg"), "step5_cg", "branch"; scenario_source = "file_block"))
      @test cg["reason"] == "invalid_request"
      @test occursin("Export the case as SCF once", cg["message"])
      @test occursin("reference its component ids", cg["message"])
      @test occursin("n1_all", cg["message"])
      bad2 = Sparlectra.to_dict(Sparlectra._run_contingency_service(case_file, cfgpath, joinpath(dir, "run_bad2"), "step5_bad2", "branch"; scenario_source = "external_file"))
      @test bad2["reason"] == "invalid_request"
    end
  end

  @testset "scenario engine case300 screening acceptance" begin
    # the operated-grid acceptance case per the test-network rule
    # (2026-09-03): case300 judges screening QUALITY, no false negatives at
    # margin 10 with the 0.005 pu trust gate (its outage 57-63 forced the
    # gate: buses 63/64/526 collapse to 0.84 pu behind a 0.0104 pu one-step
    # residual); loud SKIPPED when the untracked case file is absent
    case_path = large_case_path("case300.m")
    if case_path === nothing
      println("      scenario engine case300: SKIPPED (case300.m not in the large-case directory)")
    else
      net = Sparlectra.createNetFromMatPowerFile(filename = case_path, flatstart = false, enable_pq_gen_controllers = true, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)
      cases = vcat(Sparlectra.generateN1Branches(net), Sparlectra.generateN1Generators(net))
      full = Sparlectra.runContingencies!(net, cases)
      flagged = Sparlectra.runContingencies!(net, cases; screening_mode = :flag, screening_margin_pct = 10.0)
      violating = Set(r.name for r in full if !isempty(r.overloads) || !isempty(r.voltage_violations) || !r.converged)
      false_negatives = intersect(violating, Set(r.name for r in flagged if r.screened))
      @test sort(collect(false_negatives)) == String[]
      screened = count(r -> r.screened, flagged)
      fp = count(r -> !r.screened && !(r.name in violating), flagged)
      println("      scenario engine case300: RAN (", length(cases), " cases, ", screened, " screened, ", fp, " false positives, 0 false negatives)")
    end
  end

  @testset "scenarios workshop runs with its assertions" begin
    # scenario task step 7: the Literate workshop is executable Julia and
    # carries an @assert next to every printed number; RUNNING it here is
    # what keeps the notebook from drifting silently. Gated on the local
    # case118 (the workshop's Colab path downloads instead).
    case_path = large_case_path("case118.m")
    if case_path === nothing
      println("      scenarios workshop: SKIPPED (case118.m not in the large-case directory)")
    else
      workshop = abspath(joinpath(dirname(@__DIR__), "docs", "lit", "workshop_scenarios.jl"))
      @test isfile(workshop)
      mod = Module(:WorkshopScenariosRun)
      Base.include(mod, workshop)
      @test true
    end
  end

  @testset "scenario engine case1354 scaling timings" begin
    # the pegase case measures RUNTIME SCALING only (test-network rule
    # 2026-09-03: an OPF test instance whose base carries overloads by
    # construction says nothing about screening quality); seconds for
    # :off, :flag, :only, gated on SPARLECTRA_LARGE_CASES_DIR with a loud
    # SKIPPED line (canonical pegase conventions: rad, shift sign -1.0,
    # ratio normal)
    large_dir = get(ENV, "SPARLECTRA_LARGE_CASES_DIR", "")
    case_path = isempty(large_dir) ? "" : joinpath(large_dir, "case1354pegase.m")
    if isempty(large_dir) || !isfile(case_path)
      println("      scenario engine case1354: SKIPPED (", isempty(large_dir) ? "SPARLECTRA_LARGE_CASES_DIR not set" : string("case1354pegase.m not found under ", large_dir), ")")
    else
      net = Sparlectra.createNetFromMatPowerFile(filename = case_path, matpower_shift_unit = :rad, matpower_shift_sign = -1.0, matpower_ratio = :normal)
      cases = vcat(Sparlectra.generateN1Branches(net), Sparlectra.generateN1Generators(net))
      Sparlectra.runContingencies!(net, cases[1:2])
      t_off = @timed Sparlectra.runContingencies!(net, cases)
      t_flag = @timed Sparlectra.runContingencies!(net, cases; screening_mode = :flag, screening_margin_pct = 10.0)
      t_only = @timed Sparlectra.runContingencies!(net, cases; screening_mode = :only, screening_margin_pct = 10.0)
      @test length(t_off.value) == length(cases)
      println("      scenario engine case1354: RAN (", length(cases), " cases, off/flag/only = ", round(t_off.time; digits = 2), "/", round(t_flag.time; digits = 2), "/", round(t_only.time; digits = 2), " s, threads = ", Threads.nthreads(), ")")
    end
  end
  return nothing
end

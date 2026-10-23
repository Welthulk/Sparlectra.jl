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

# file: test/test_legacy_scenarios.jl
# purpose: legacy profile: the testsets of test/test_scenarios.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here from test_scenarios.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_scenarios_tests()
  # cut from testset: scenario patch model
  (function ()
      case, net = _scenario_case14()
      index = ScenarioIndex(case)

      @testset "bitwise restore over every N-1 scenario" begin (function ()
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
          @test ("sp_casePST", scf_roundtrip_field_diffs(base, work)) == ("sp_casePST", String[])
        end
      end)() end

      @testset "engine: worker resets bitwise between evaluations" begin (function ()
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
          Sparlectra._ScenarioOutageItem(Sparlectra._CaseLabel("branch out", :branch, ["branch out"], 1.0), [(:branch, branch_internal)]),
          Sparlectra._ScenarioOutageItem(Sparlectra._CaseLabel("gen out", :gen, ["gen out"], 1.0), [(:gen, index.internal_by_id[gen_id])]),
          Sparlectra._ScenarioPatchItem(Scenario(name = "load up", ops = [PatchOp(op = :scale, target = :load, id = load_id2, factor = 1.3)]), Sparlectra._CaseLabel("load up", :explicit, ["load up"], 1.0)),
        ]
        for it in items
          r = Sparlectra.evaluate!(engine, worker, it)
          @test r isa Sparlectra.ContingencyResult
          @test (it, scf_roundtrip_field_diffs(engine.template, worker.net)) == (it, String[])
        end
      end)() end

      @testset "events on the outage path (#355: D1, links, kinds)" begin (function ()
        case60 = read_scf_json(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case60.scf.json")))
        net60 = Sparlectra.build_net(case60)
        idx60 = ScenarioIndex(case60)
        out(id) = PatchOp(op = :status, target = idx60.kind_by_id[id], id = id, value = 0.0)
        # D1: two generator outages in one scenario. The ops were applied in
        # op order by position, so the second deletion hit the unit that had
        # moved into the slot ([a, b] gave 0.8407 pu and 49.2 %, [b, a] 0.8664
        # pu and 87.5 % on a35aa395). Both orders are now the same event.
        gens = sort([id for (id, k) in idx60.kind_by_id if k === :generator]; by = id -> idx60.internal_by_id[id])
        a, b = gens[1], gens[2]
        neutral = first(id for (id, k) in idx60.kind_by_id if k === :load)
        scen = [
          Scenario(name = "ab", ops = [out(a), out(b)]),
          Scenario(name = "ba", ops = [out(b), out(a)]),
          # the same event through the patch path (a neutral op routes it there)
          Scenario(name = "ab patched", ops = [out(a), out(b), PatchOp(op = :scale, target = :load, id = neutral, factor = 1.0)]),
          Scenario(name = "ba patched", ops = [PatchOp(op = :scale, target = :load, id = neutral, factor = 1.0), out(b), out(a)]),
        ]
        r = redirect_stdout(() -> runScenarios!(net60, scen; index = idx60, parallel_enabled = false), devnull)
        @test all(x -> x.converged, r)
        @test _scenario_values_equal(r[1], r[2])
        @test _scenario_values_equal(r[1], r[3])
        @test _scenario_values_equal(r[3], r[4])
        # the patch path removes exactly the two units, in either op order
        before = [Sparlectra.getCompName(ps.comp) for ps in net60.prosumpsVec]
        work = deepcopy(net60)
        for ops in (scen[1].ops, scen[2].ops)
          undo = apply!(work, ops, idx60)
          @test length(work.prosumpsVec) == length(before) - 2
          kept = Set(objectid(ps) for ps in work.prosumpsVec)
          restore!(work, undo)
          removed = [k for (k, ps) in enumerate(work.prosumpsVec) if !(objectid(ps) in kept)]
          @test removed == sort([idx60.internal_by_id[a], idx60.internal_by_id[b]])
        end
        @test scf_roundtrip_field_diffs(net60, work) == String[]

        # a link outage (the Neubach coupler) runs on the outage path and
        # splits the station; its row names the link and reports the kind
        link_case = only(generateLinkOutages(net60; names = contingency_names(net60, case60)))
        @test link_case.kind === :link
        rl = only(redirect_stdout(() -> runContingencies!(net60, [link_case]; names = contingency_names(net60, case60), parallel_enabled = false), devnull))
        @test rl.kind === :link && rl.elements == link_case.elements
        @test rl.converged
        # markIsolatedBuses! with links = true marks by the link groups of the
        # power flow: with every branch at one coupler bus open, that bus
        # hangs on the other through the coupler and stays in service, the
        # default (branches only) isolates it; with the branches at both
        # coupler buses open the group has no branch and both are isolated
        cpl = deepcopy(net60)
        cpl_link = only(l for l in cpl.linkVec if l.status == 1)
        cpl_buses = [Int(cpl_link.fromBus), Int(cpl_link.toBus)]
        for br in cpl.branchVec
          (Int(br.fromBus) == cpl_buses[1] || Int(br.toBus) == cpl_buses[1]) && Sparlectra.setBranchStatus!(br, false)
        end
        @test !any(in(Sparlectra.markIsolatedBuses!(net = deepcopy(cpl), links = true)), cpl_buses)
        @test cpl_buses[1] in Sparlectra.markIsolatedBuses!(net = deepcopy(cpl))
        for br in cpl.branchVec
          (Int(br.fromBus) in cpl_buses || Int(br.toBus) in cpl_buses) && Sparlectra.setBranchStatus!(br, false)
        end
        @test all(in(Sparlectra.markIsolatedBuses!(net = cpl, links = true)), cpl_buses)
        # links = true marks by the link groups of the power flow, bus by bus
        # equal to the solver's isolation: C behind B in service, the
        # branch-less group E-F and the group G-H whose only branch runs
        # inside it isolated, J behind an open link isolated, the chain
        # D-M-N in service
        grp = Net(name = "link_groups", baseMVA = 100.0)
        foreach(b -> addBus!(net = grp, busName = b, vn_kV = 110.0), ("A", "B", "C", "D", "E", "F", "G", "H", "J", "M", "N"))
        for (f, t) in (("A", "B"), ("A", "D"), ("B", "D"), ("G", "H"))
          addPIModelACLine!(net = grp, fromBus = f, toBus = t, r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
        end
        for (f, t, s) in (("B", "C", 1), ("E", "F", 1), ("G", "H", 1), ("B", "J", 0), ("D", "M", 1), ("M", "N", 1))
          addLink!(net = grp, fromBus = f, toBus = t, status = s)
        end
        addProsumer!(net = grp, busName = "A", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "A")
        foreach(b -> addProsumer!(net = grp, busName = b, type = "ENERGYCONSUMER", p = 2.0, q = 0.5), ("B", "C", "D", "E", "G", "N"))
        @test first(validate!(net = grp))
        wnet, reps, _ = Sparlectra._merged_pf_net(grp)
        solver_iso = [i for i in eachindex(grp.nodeVec) if reps[i] in wnet.isoNodes]
        @test [grp.busDict[b] for b in ("E", "F", "G", "H", "J")] == solver_iso
        @test Sparlectra.markIsolatedBuses!(net = deepcopy(grp), links = true) == solver_iso
        # the patch path of runScenarios! counts a cut-off link group whole:
        # opening A-B of the radial net cuts off B and C (15 MW), in the ring
        # no branch outage cuts off C
        for (ring, shed) in ((false, [15.0, 20.0]), (true, [0.0, 0.0, 0.0]))
          lcase = read_scf_json(exportSCF(link_group_net(; ring = ring); file = joinpath(mktempdir(), "link_group.scf.json")))
          lidx = ScenarioIndex(lcase)
          lscen = [Scenario(name = string("out ", id), ops = [PatchOp(op = :status, target = :branch, id = id, value = 0.0)]) for id in sort!([id for (id, k) in lidx.kind_by_id if k === :branch])]
          lres = redirect_stdout(() -> runScenarios!(Sparlectra.build_net(lcase), lscen; index = lidx, parallel_enabled = false), devnull)
          @test [r.shed_load_mw for r in lres] == shed
          @test count(r -> r.error !== nothing && occursin("B, C cut off", r.error), lres) == (ring ? 0 : 1)
        end
        # a generated double circuit equals the hand-written explicit scenario
        # of the same two outages, value for value
        names60 = contingency_names(net60, case60)
        dc = first(generateDoubleCircuits(net60; names = names60))
        rdc = only(redirect_stdout(() -> runContingencies!(net60, [dc]; names = names60, parallel_enabled = false), devnull))
        sdc = scenario_from_case(dc, net60, idx60; names = names60)
        rsc = only(redirect_stdout(() -> runScenarios!(net60, [Scenario(name = "by hand", ops = sdc.ops)]; index = idx60, parallel_enabled = false), devnull))
        @test rdc.kind === :double_circuit && rsc.kind === :explicit
        @test rdc.elements == rsc.elements
        @test _scenario_values_equal(rdc, rsc)
      end)() end

      @testset "correlated event generators on the demo cases (#355)" begin (function ()
        demo(name) = (c = read_scf_json(abspath(joinpath(dirname(@__DIR__), "data", "scf", name * ".scf.json"))); n = Sparlectra.build_net(c); (c, n, contingency_names(n, c)))
        # sp_case14: one double circuit (ids 26/27); its name group is the same
        # member set and is removed as a duplicate
        c14, n14, m14 = demo("sp_case14")
        dc14 = generateDoubleCircuits(n14; names = m14)
        @test length(dc14) == 1
        @test sort(dc14[1].elements) == ["Ostheim_Neuwiese_1", "Ostheim_Neuwiese_2"]
        both14 = generateCorrelatedEvents(n14; names = m14, double_circuits = true, name_groups = true)
        @test [c.kind for c in both14] == [:double_circuit]
        # sp_case60: six double circuits, the coupler, Neubach as two busbar
        # sections and one station, the name group across the two sections
        c60, n60, m60 = demo("sp_case60")
        @test length(generateDoubleCircuits(n60; names = m60)) == 6
        @test [c.elements for c in generateLinkOutages(n60; names = m60)] == [["Link_1"]]
        bb = [c for c in generateBusbarEvents(n60; names = m60) if startswith(c.name, "BB Neubach_110")]
        @test sort([c.name for c in bb]) == ["BB Neubach_110", "BB Neubach_110b"]
        @test all(c -> "Link_1" in c.elements, bb)
        sub = only(c for c in generateSubstationEvents(n60; names = m60) if c.name == "SUB Neubach_110")
        @test !("Link_1" in sub.elements) && length(sub.elements) == 4
        ng = only(c for c in generateNameGroups(n60; names = m60) if c.name == "NG Neubach_Waldsee")
        @test length(ng.elements) == 2 && all(e -> startswith(e, "Neubach_Waldsee#"), ng.elements)
        # sp_case188: seven double circuits, the coupler, the Eichried_Falkwik
        # name group and the Falkwik station
        c188, n188, m188 = demo("sp_case188")
        @test length(generateDoubleCircuits(n188; names = m188)) == 7
        @test length(generateLinkOutages(n188; names = m188)) == 1
        @test sort(only(c for c in generateNameGroups(n188; names = m188) if c.name == "NG Eichried_Falkwik").elements) == ["Eichried_Falkwik_1", "Eichried_Falkwik_2"]
        @test any(c -> c.name == "SUB Falkwik_380", generateSubstationEvents(n188; names = m188))
        # no name group joins branches of different nominal voltage, on any
        # demo case; a plain net without file names has no name groups
        for (c, n, m) in (demo("sp_case5"), (c14, n14, m14), (c60, n60, m60), (c188, n188, m188))
          for ev in generateNameGroups(n, names = m)
            ks = [Sparlectra._resolve_named_branch(m, e) for e in ev.elements]
            vns = unique(max(Sparlectra.getNodeVn(n.nodeVec[Int(n.branchVec[k].fromBus)]), Sparlectra.getNodeVn(n.nodeVec[Int(n.branchVec[k].toBus)])) for k in ks)
            @test length(vns) == 1
          end
        end
        @test isempty(redirect_stderr(() -> generateNameGroups(n60), devnull))
        # one member set appears once, as the first class: a single-bus station
        # is :busbar with busbars on, :substation with stations only
        ev = generateCorrelatedEvents(n60; names = m60, busbars = true, substations = true, links = true, name_groups = true)
        @test length(unique(Set(Sparlectra._contingency_case_members(n60, m60, c)) for c in ev)) == length(ev)
        @test count(c -> c.kind === :substation, ev) == 1
        subs_only = generateCorrelatedEvents(n60; names = m60, double_circuits = false, substations = true)
        @test all(c -> c.kind === :substation, subs_only) && length(subs_only) > 1
        # N-2 pairs: every unordered pair of closed branches, the pair that is
        # also a double circuit listed once (as the event), and the cap stops
        # before any case is built
        @test length(generateBranchPairs(n14; names = m14)) == 136
        @test length(generateBranchPairs(n60; names = m60)) == 2701
        with_pairs = generateCorrelatedEvents(n14; names = m14, pairs = true)
        @test count(c -> c.kind === :double_circuit, with_pairs) == 1
        @test count(c -> c.kind === :pair, with_pairs) == 135
        @test all(c -> c.kind !== :pair || length(c.elements) == 2, with_pairs)
        @test_throws ArgumentError generateBranchPairs(n14; names = m14, max_pairs = 135)
      end)() end

  end)()
  # cut from the runner function
  @testset "scenario engine: screening of events and pairs (#355)" begin (function ()
    # Every generated event and every N-2 pair is screened like a single
    # outage (all branch members out of the Y-bus, all units out of the
    # injections in one step). Acceptance: no false negative against the full
    # solve, on the shipped SCF cases and sp_case118, events and pairs.
    violating(rs) = Set(r.name for r in rs if !isempty(r.overloads) || !isempty(r.voltage_violations) || !r.converged)
    function no_false_negative(net, cases, names)
      full = redirect_stdout(() -> runContingencies!(net, cases; names = names), devnull)
      flagged = redirect_stdout(() -> runContingencies!(net, cases; names = names, screening_mode = :flag, screening_margin_pct = 10.0), devnull)
      missed = intersect(violating(full), Set(r.name for r in flagged if r.screened))
      return (missed = sort(collect(missed)), screened = count(r -> r.screened, flagged), total = length(cases))
    end
    # the same judgement for the large lists without solving every case
    # twice: a false negative can only be a case the screening KEEPS as an
    # estimate (the decision _evaluate_item takes), so only those are solved
    # in full (sp_case188 has 23436 pairs)
    function no_false_negative_kept(net, cases, names)
      eng = redirect_stdout(() -> Sparlectra.ScenarioEngine(net; screening_mode = :flag, screening_margin_pct = 10.0, names = names), devnull)
      items = Sparlectra._engine_items_from_cases(eng.template, names, cases)
      kept = ContingencyCase[]
      for (c, it) in zip(cases, items)
        it isa Sparlectra._ScenarioOutageItem || continue
        scr = Sparlectra._screen_outage(eng, it)
        (scr.estimate !== nothing && !scr.flagged) && push!(kept, c)
      end
      full = isempty(kept) ? ContingencyResult[] : redirect_stdout(() -> runContingencies!(net, kept; names = names), devnull)
      return (missed = sort(collect(violating(full))), screened = length(kept), total = length(cases))
    end
    scf(name) = (c = read_scf_json(abspath(joinpath(dirname(@__DIR__), "data", "scf", name * ".scf.json"))); n = Sparlectra.build_net(c); (n, contingency_names(n, c)))
    events(n, m) = generateCorrelatedEvents(n; names = m, busbars = true, substations = true, links = true, name_groups = true)
    counts = String[]
    for name in ("sp_case14", "sp_case60", "sp_case188")
      n, m = scf(name)
      check = name == "sp_case188" ? no_false_negative_kept : no_false_negative
      ev = check(n, events(n, m), m)
      @test ev.missed == String[]
      pr = check(n, generateBranchPairs(n; names = m), m)
      @test pr.missed == String[]
      push!(counts, "$(name) events $(ev.screened)/$(ev.total), pairs $(pr.screened)/$(pr.total)")
    end
    n118 = Sparlectra.createNetFromMatPowerFile(filename = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "sp_case118.m")), flatstart = false, enable_pq_gen_controllers = true, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)
    m118 = contingency_names(n118)
    ev118 = no_false_negative_kept(n118, events(n118, m118), m118)
    @test ev118.missed == String[]
    push!(counts, "sp_case118 events $(ev118.screened)/$(ev118.total)")
    # the reasons: a coupler is never screened, an event that splits the
    # graph neither, a double circuit of the meshed ring is
    n60, m60 = scf("sp_case60")
    eng = redirect_stdout(() -> Sparlectra.ScenarioEngine(n60; screening_mode = :flag, names = m60), devnull)
    reason(case) = Sparlectra._screen_outage(eng, only(Sparlectra._engine_items_from_cases(eng.template, m60, [case]))).reason
    @test reason(only(generateLinkOutages(n60; names = m60))) === :link
    @test reason(first(c for c in generateBusbarEvents(n60; names = m60) if c.name == "BB Altau_220")) === :splits_graph
    @test reason(first(generateDoubleCircuits(n60; names = m60))) in (:ok, :structural_pv_gen)
    # above 32 changed Jacobian rows the dense correction is not tried
    J0 = sparse(1.0LinearAlgebra.I, 40, 40)
    @test Sparlectra._screen_woodbury_solve(lu(J0), J0, J0 + sparse(0.5LinearAlgebra.I, 40, 40), ones(40)) === :rank
    println("      scenario engine event screening: RAN (", join(counts, "; "), ", 0 false negatives)")
  end)() end

  @testset "scenario service sources (step 5)" begin (function ()
    # the service request picks the scenario source and screening mode: an
    # SCF case with its own scenarios block runs via file_block, an external
    # scenario JSON via external_file, and the two screening columns plus the
    # screening metadata arrive exactly when screening is active
    cfgpath = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH
    cfg = Sparlectra.load_sparlectra_config(cfgpath; reload = true)
    mpc = Sparlectra.MatpowerIO.read_case(abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m")); legacy_compat = true)
    case = Sparlectra.convert_case(MatpowerAdapter(), mpc, Sparlectra.matpower_adapter_options(cfg))
    index = ScenarioIndex(case)
    branch_id = first(r.id for r in case.data.line)
    gen_id = first(id for (id, k) in index.kind_by_id if k === :generator)
    set = ScenarioSet(scenarios = [Scenario(name = "line out", ops = [PatchOp(op = :status, target = :branch, id = branch_id, value = 0.0)]), Scenario(name = "gen setpoint", ops = [PatchOp(op = :set, target = :generator, id = gen_id, field = :p, value = 10.0)])])
    case.sparlectra.scenarios = scenario_set_dict(set)
    dir = mktempdir()
    case_file = joinpath(dir, "step5_case.scf.json")
    write_scf_json(case, case_file)
    redirect_stdout(devnull) do
      res = SparlectraApp._run_contingency_service(case_file, cfgpath, joinpath(dir, "run_fb"), "step5_fb", "branch"; scenario_source = "file_block", screening_mode = "flag")
      d1 = SparlectraApp.to_dict(res)
      @test d1["status"] == "succeeded"
      @test d1["metadata"]["contingency_screening_mode"] == "flag"
      @test d1["metadata"]["contingency_cases_source"] == "scenario_file_block"
      @test haskey(d1["metadata"], "contingency_screened")
      csv1 = run_csv(joinpath(dir, "run_fb", "contingency_n1.csv"))
      @test csv1.header[(end-1):end] == ["screened", "screening_estimate"]
      @test length(csv1.rows) == 2
      # a requested screening the network does not allow (the shipped Q(U)
      # case carries a voltage-dependent controller) says why in run.log
      qu = abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14_qu.scf.json"))
      resqu = SparlectraApp._run_contingency_service(qu, cfgpath, joinpath(dir, "run_qu"), "step5_qu", "branch"; screening_mode = "flag")
      @test SparlectraApp.to_dict(resqu)["status"] == "succeeded"
      @test any(l -> occursin("screening: flag requested but not available: the network has voltage-dependent controllers", l), readlines(joinpath(dir, "run_qu", "run.log")))
      # external scenario JSON, screening explicitly off: classic CSV columns
      scen_file = joinpath(dir, "scenarios.json")
      write(scen_file, Sparlectra.scf_json_string(scenario_set_dict(set)))
      res2 = SparlectraApp._run_contingency_service(case_file, cfgpath, joinpath(dir, "run_ext"), "step5_ext", "branch"; scenario_source = "external_file", scenario_file = scen_file, screening_mode = "off")
      d2 = SparlectraApp.to_dict(res2)
      @test d2["status"] == "succeeded"
      @test d2["metadata"]["contingency_screening_mode"] == "off"
      @test d2["metadata"]["contingency_cases_source"] == "scenario_external_file"
      csv2 = readlines(joinpath(dir, "run_ext", "contingency_n1.csv"))
      # screening off: the classic columns, ending with the error column
      @test run_csv(joinpath(dir, "run_ext", "contingency_n1.csv")).header[end] == "error"
      # #456: the scenario source solves with the configured distributed
      # slack (it used to get none while run.log announced it): the same
      # scenarios with power_flow.distributed_slack on give other rows
      res2ds = SparlectraApp._run_contingency_service(case_file, cfgpath, joinpath(dir, "run_ext_ds"), "step5_ext_ds", "branch"; scenario_source = "external_file", scenario_file = scen_file, screening_mode = "off", config_overrides = Dict{String,Any}("power_flow.distributed_slack.enabled" => true))
      @test SparlectraApp.to_dict(res2ds)["status"] == "succeeded"
      csv2ds = readlines(joinpath(dir, "run_ext_ds", "contingency_n1.csv"))
      @test length(csv2ds) == length(csv2)
      @test csv2ds[2:end] != csv2[2:end]
      # the N-1 path honours power_flow.qlimits.* like the single run: the
      # service forwarded none of the Q-limit settings before 0.30.2, so a
      # branch N-1 of sp_case14 with every PV bus locked to PQ gave the same
      # rows as without (measured 0 of 17 rows different; with the
      # forwarding 14 of 17 differ)
      sp14 = joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")
      nb14 = length(Sparlectra.importSCF(sp14).nodeVec)
      for (sub, ov) in (("n1_plain", Dict{String,Any}()), ("n1_lock", Dict{String,Any}("power_flow.qlimits.lock_pv_to_pq_buses" => collect(1:nb14))))
        @test SparlectraApp.to_dict(SparlectraApp._run_contingency_service(sp14, cfgpath, joinpath(dir, sub), sub, "branch"; scenario_source = "n1_branches", screening_mode = "off", config_overrides = ov))["status"] == "succeeded"
      end
      csv_plain = readlines(joinpath(dir, "n1_plain", "contingency_n1.csv"))
      csv_lock = readlines(joinpath(dir, "n1_lock", "contingency_n1.csv"))
      @test length(csv_lock) == length(csv_plain)
      @test csv_lock[2:end] != csv_plain[2:end]
      # an n1 source works on any format and records itself
      res3 = SparlectraApp._run_contingency_service(case_file, cfgpath, joinpath(dir, "run_n1"), "step5_n1", "branch"; scenario_source = "n1_generators", screening_mode = "off")
      d3 = SparlectraApp.to_dict(res3)
      @test d3["status"] == "succeeded"
      @test d3["metadata"]["contingency_cases_source"] == "n1_generators"
      # rejections carry invalid_request
      bad = SparlectraApp.to_dict(SparlectraApp._run_contingency_service(case_file, cfgpath, joinpath(dir, "run_bad"), "step5_bad", "branch"; scenario_source = "bogus"))
      @test bad["reason"] == "invalid_request"
      # a CGMES delivery with an id-addressed scenario source is rejected
      # WITH the way out named (export as SCF once, then reference its ids)
      cgdir = mktempdir()
      cg = SparlectraApp.to_dict(SparlectraApp._run_contingency_service(cgdir, cfgpath, joinpath(dir, "run_cg"), "step5_cg", "branch"; scenario_source = "file_block"))
      @test cg["reason"] == "invalid_request"
      @test occursin("Export the case as SCF once", cg["message"])
      @test occursin("reference its component ids", cg["message"])
      @test occursin("n1_all", cg["message"])
      bad2 = SparlectraApp.to_dict(SparlectraApp._run_contingency_service(case_file, cfgpath, joinpath(dir, "run_bad2"), "step5_bad2", "branch"; scenario_source = "external_file"))
      @test bad2["reason"] == "invalid_request"
    end
  end)() end

  @testset "service: correlated events, pairs, display names (#355)" begin (function ()
    cfgpath = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH
    dir = mktempdir()
    sp14 = joinpath(dir, "sp_case14.scf.json")
    cp(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")), sp14)
    run(sub, ov; kw...) = redirect_stdout(devnull) do
      SparlectraApp.to_dict(SparlectraApp._run_contingency_service(sp14, cfgpath, joinpath(dir, sub), sub, "branch"; scenario_source = "n1_branches", screening_mode = "off", config_overrides = ov, kw...))
    end
    # events on: the double circuit runs next to the 17 single outages, the
    # rows carry the file names, run.log and the metadata count by kind
    d = run("ev", Dict{String,Any}("contingency.events.enabled" => true))
    @test d["status"] == "succeeded"
    @test d["metadata"]["contingency_cases_by_kind"] == Dict{String,Any}("branch" => 17, "double_circuit" => 1)
    log = read(joinpath(dir, "ev", "run.log"), String)
    @test occursin("correlated events: 1 double_circuit", log)
    csv = read(joinpath(dir, "ev", "contingency_n1.csv"), String)
    @test occursin("Ostheim_Neuwiese_1|Ostheim_Neuwiese_2", csv)
    @test occursin("DC Ostheim_110-Neuwiese_110", csv)
    # pairs: 136 pairs of 17 branches, the double-circuit pair listed once
    # (as the event); the cap stops the run before solving
    dp = run("pairs", Dict{String,Any}("contingency.events.enabled" => true, "contingency.pairs.enabled" => true))
    @test dp["metadata"]["contingency_cases_by_kind"] == Dict{String,Any}("branch" => 17, "double_circuit" => 1, "pair" => 135)
    dcap = run("cap", Dict{String,Any}("contingency.pairs.enabled" => true, "contingency.pairs.max_pairs" => 100))
    @test dcap["status"] == "failed" && dcap["reason"] == "contingency_too_many_pairs"
    @test occursin("136 pairs", dcap["message"])
    # a MATPOWER case runs plain N-1 with events switched on, and says why
    mp = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m"))
    dm = redirect_stdout(devnull) do
      SparlectraApp.to_dict(SparlectraApp._run_contingency_service(mp, cfgpath, joinpath(dir, "mp"), "mp", "branch"; scenario_source = "n1_branches", screening_mode = "off", config_overrides = Dict{String,Any}("contingency.events.enabled" => true)))
    end
    @test dm["status"] == "succeeded"
    @test occursin("correlated events and pairs: not run, they need a Sparlectra Case Format case", read(joinpath(dir, "mp", "run.log"), String))
  end)() end

  return true
end

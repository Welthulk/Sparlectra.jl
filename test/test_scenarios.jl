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
# purpose: the patch model, its load-time validation
#          (rejections carry the scenario name and the op index), the
#          sparlectra.scenarios file round trip, the legacy contingencies
#          mapping, and the N-1 expansion equality with the existing
#          generators, on the tracked MATPOWER warmup case
#          (load_fixture_net: the downloaded case14 is gone; the length
#          assertions are all RELATIVE to the case, so only the fixture
#          changed).

# the field-by-field round-trip comparison of the SCF tests, shared here
isdefined(@__MODULE__, :scf_roundtrip_field_diffs) || include("test_scf_support.jl")

# every value of two result rows except the label (name, kind, elements)
_scenario_values_equal(x, y) = all(f -> f in (:name, :kind, :elements) || isequal(getproperty(x, f), getproperty(y, f)), fieldnames(Sparlectra.ContingencyResult))

function _scenario_case14()
  cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
  mpc = Sparlectra.MatpowerIO.read_case(abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m")); legacy_compat = true)
  case = Sparlectra.convert_case(MatpowerAdapter(), mpc, Sparlectra.matpower_adapter_options(cfg))
  net = Sparlectra.build_net(case; config = cfg)
  return case, net
end

# A hand-built QSTS index over a built net, as the CGMES series union will
# carry it: dense ids over buses, branches (ratio 0 is a line), prosumers
# (the reference unit stands in as the external grid), shunts, links and
# tap controllers. `ids` maps (group, internal index) to the id, the group
# being :bus, :branch, :prosumer, :shunt, :link or :tap.
function _qsts_test_index(net::Sparlectra.Net; parked::Dict{Int,Sparlectra.ProSumer} = Dict{Int,Sparlectra.ProSumer}(), origin::Symbol = :cgmes_series)
  kinds = Dict{Int,Symbol}()
  internals = Dict{Int,Int}()
  names = Dict{Int,String}()
  mrids = Dict{String,Int}()
  ids = Dict{Tuple{Symbol,Int},Int}()
  function add!(group::Symbol, kind::Symbol, internal::Int)
    id = length(kinds) + 1
    kinds[id] = kind
    internals[id] = internal
    names[id] = string(kind, " ", internal)
    mrids[string("_", kind, "_", internal)] = id
    ids[(group, internal)] = id
    return nothing
  end
  foreach(b -> add!(:bus, :bus, b), eachindex(net.nodeVec))
  foreach(i -> add!(:branch, net.branchVec[i].ratio != 0.0 ? :transformer : :branch, i), eachindex(net.branchVec))
  foreach(i -> add!(:prosumer, net.prosumpsVec[i].referencePri !== nothing ? :external_grid : Sparlectra.isGenerator(net.prosumpsVec[i]) ? :generator : :load, i), eachindex(net.prosumpsVec))
  foreach(i -> add!(:shunt, :shunt, i), eachindex(net.shuntVec))
  foreach(i -> add!(:link, :link, i), eachindex(net.linkVec))
  foreach(i -> add!(:tap, :tap_controller, i), eachindex(Sparlectra._tap_controllers(net)))
  return ScenarioIndex(kinds, internals, names, mrids, parked, origin), ids
end

function run_scenario_patch_tests()
  @testset "scenario patch model" begin (function ()
    case, net = _scenario_case14()
    index = ScenarioIndex(case)
    branch_id = first(r.id for r in case.data.line)
    trafo_id = first(r.id for r in case.data.generic_branch)
    gen_id = first(r.id for r in case.data.sym_gen)
    load_id = first(r.id for r in case.data.sym_load)

    @testset "rejections name the scenario and the op index" begin (function ()
      # the positive case first: every op kind on its own target passes
      good = ScenarioSet(
        scenarios = [
          Scenario(name = "line out", ops = [PatchOp(op = :status, target = :branch, id = branch_id, value = 0.0)]),
          Scenario(name = "gen setpoint", ops = [PatchOp(op = :set, target = :generator, id = gen_id, field = :p, value = 25.0)]),
          Scenario(name = "load scale", weight = 2.0, ops = [PatchOp(op = :scale, target = :load, id = load_id, factor = 1.2)]),
          Scenario(name = "tap step", ops = [PatchOp(op = :set, target = :transformer, id = trafo_id, field = :tap_pos, value = 2.0)]),
        ],
      )
      @test validate_scenarios(good, index) === good

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

      dup = ScenarioSet(scenarios = [Scenario(name = "twin", ops = [PatchOp(op = :status, target = :branch, id = branch_id, value = 0.0)]), Scenario(name = "twin", ops = [PatchOp(op = :status, target = :branch, id = branch_id, value = 1.0)])])
      @test_throws ArgumentError validate_scenarios(dup, index)

      mismatch = ScenarioSet(scenarios = [Scenario(name = "class clash", ops = [PatchOp(op = :status, target = :generator, id = branch_id, value = 0.0)])])
      err3 = try
        validate_scenarios(mismatch, index)
        ""
      catch e
        sprint(showerror, e)
      end
      @test occursin("is a branch, not a generator", err3)
    end)() end

    @testset "file round trip of the scenarios block" begin (function ()
      set = ScenarioSet(
        scenarios = [
          Scenario(name = "double outage", weight = 0.5, ops = [PatchOp(op = :status, target = :branch, id = branch_id, value = 0.0), PatchOp(op = :status, target = :transformer, id = trafo_id, value = 0.0)]),
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
    end)() end

    @testset "every shipped SCF case resolves its components" begin (function ()
      # ScenarioIndex reads the internal indices from the build-order fields
      # of the extra block. Three hand-written feeder cases carried none, so
      # every op on them failed in apply! ("component id 6 has no internal
      # index"), the scenarios shipped in feeder3_hardpv_vsPGM included. Each
      # index must exist and point at the element on the row's buses.
      local dir, c, idx, n, bus_of, undo
      at(vec, ix, id) = vec[ix.internal_by_id[id]]
      ends(x) = Set(Int.((x.fromBus, x.toBus)))
      dir = abspath(joinpath(dirname(@__DIR__), "data", "scf"))
      for f in sort!(filter(name -> endswith(name, ".scf.json"), readdir(dir)))
        c = read_scf_json(joinpath(dir, f))
        idx = ScenarioIndex(c)
        @test (f, sort!([id for id in keys(idx.kind_by_id) if !haskey(idx.internal_by_id, id)])) == (f, Int[])
        n = Sparlectra.build_net(c)
        bus_of = Dict(id => k for (k, id) in enumerate(Sparlectra.scf_node_ids_in_build_order(c)))
        @test (f, [r.id for r in vcat(c.data.source, c.data.sym_gen, c.data.sym_load) if Sparlectra._scf_prosumer_bus(at(n.prosumpsVec, idx, r.id)) != bus_of[r.node]]) == (f, Int[])
        @test (f, [r.id for r in vcat(c.data.line, c.data.generic_branch) if ends(at(n.branchVec, idx, r.id)) != Set((bus_of[r.from_node], bus_of[r.to_node]))]) == (f, Int[])
        @test (f, [r.id for r in c.data.shunt if at(n.shuntVec, idx, r.id).busIdx != bus_of[r.node]]) == (f, Int[])
        @test (f, [r.id for r in c.data.link if ends(at(n.linkVec, idx, r.id)) != Set((bus_of[r.from_node], bus_of[r.to_node]))]) == (f, Int[])
        if f == "feeder3_hardpv_vsPGM.scf.json"
          # the case's own scenarios block applies
          for s in scf_case_scenarios(c).scenarios
            undo = apply!(n, s.ops, idx)
            @test (s.name, isempty(undo.entries)) == (s.name, false)
            restore!(n, undo)
          end
        end
      end
    end)() end

    @testset "apply! keeps the bus load sums for their readers" begin (function ()
      # after a load outage the island totals (detect_ac_islands), the
      # passive-bus test (findPassiveBuses) and the link-flow allocation
      # (calcLinkFlowsKCL!) read the bus load sums: a patched net must give
      # what a net built without that load gives, and restore! the base
      path = abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case60.scf.json"))
      # a load moved to a link end, so the link-flow allocation sees it;
      # one load scaled and the moved one taken out (the reference case
      # writes the same change into the case data)
      function edited(; outage::Bool)
        c = read_scf_json(path)
        moved = c.data.sym_load[1]
        moved.node = first(c.data.link).from_node
        if outage
          # the reference drops the row: build_net keeps a sym_load whose
          # status is 0, so the status alone is no outage of the case data
          other = c.data.sym_load[2]
          other.p_specified *= 2.0
          other.q_specified *= 2.0
          deleteat!(c.data.sym_load, 1)
        end
        return c
      end
      # `local` is load-bearing: the enclosing test function holds `case`,
      # `net` and `index` of case14, and a plain assignment in this nested
      # function would overwrite them for every later testset
      local case, ops, idx, base, net, ref, undo, solved
      case = edited(outage = false)
      ops = [PatchOp(op = :scale, target = :load, id = case.data.sym_load[2].id, factor = 2.0), PatchOp(op = :status, target = :load, id = case.data.sym_load[1].id, value = 0.0)]
      idx = ScenarioIndex(case)
      base = Sparlectra.build_net(case)
      net = Sparlectra.build_net(case)
      ref = Sparlectra.build_net(edited(outage = true))
      undo = apply!(net, ops, idx)
      sums(n) = [(something(nd._pƩLoad, 0.0), something(nd._qƩLoad, 0.0)) for nd in n.nodeVec]
      @test all(isapprox.(first.(sums(net)), first.(sums(ref)); atol = 1e-9)) && all(isapprox.(last.(sums(net)), last.(sums(ref)); atol = 1e-9))
      island_load(n) = sort!([round(r.total_load_p_mw; digits = 9) for r in Sparlectra.detect_ac_islands(n).rows])
      @test island_load(net) == island_load(ref)
      @test sort!(collect(Sparlectra.findPassiveBuses(net))) == sort!(collect(Sparlectra.findPassiveBuses(ref)))
      restore!(net, undo)
      @test scf_roundtrip_field_diffs(base, net) == String[]
      # the link flow after a solve: patched copy against the reference
      solved = Sparlectra.build_net(case)
      apply!(solved, ops, idx)
      for n in (solved, ref)
        redirect_stdout(() -> runpf!(n, 40, 1e-9, 0), devnull)
        calcLinkFlowsKCL!(n)
      end
      @test isapprox(first(solved.linkVec).pFlow_MW, first(ref.linkVec).pFlow_MW; atol = 1e-6)
      @test isapprox(first(solved.linkVec).qFlow_MVar, first(ref.linkVec).qFlow_MVar; atol = 1e-6)
    end)() end

    @testset "apply and restore: tap patch semantics" begin (function ()
      # a tap patch on an UNREGULATED transformer changes the ratio the
      # solver sees and restores it (first-class)
      tapped = findfirst(br -> br.has_ratio_tap, net.branchVec)
      if tapped === nothing
        println("      scenarios: tap-patch semantics SKIPPED (warmup case build carries no ratio tap changer)")
      else
        tap_id = first(id for (id, k) in index.kind_by_id if k === :transformer && get(index.internal_by_id, id, 0) == tapped)
        br = net.branchVec[tapped]
        before_ratio = Sparlectra.calcBranchRatio(br)
        # a step marked as the assumed middle position: the stated tap_pos
        # takes the mark off, restore! puts it back
        br.ratio_step_assumed = true
        undo = apply!(net, [PatchOp(op = :set, target = :transformer, id = tap_id, field = :tap_pos, value = 2.0)], index)
        @test Sparlectra.calcBranchRatio(net.branchVec[tapped]) != before_ratio
        cleared = !br.ratio_step_assumed
        restore!(net, undo)
        @test Sparlectra.calcBranchRatio(net.branchVec[tapped]) == before_ratio
        @test cleared && br.ratio_step_assumed
        br.ratio_step_assumed = false
        # a CGMES changer (tap_ratio_grid :linear) steps ratio + n * step;
        # the fraction rule moved its ratio the other way (3W CGMES fixture,
        # 8 steps of 1 % from 1.0: 0.926 instead of 1.08)
        _, lnet = _scenario_case14()
        lbr = lnet.branchVec[tapped]
        Sparlectra.applyTapNameplate!(lbr; tap_step = 0.01, tap_min_step = -10.0, tap_max_step = 10.0, tap_ratio_grid = :linear)
        lundo = apply!(lnet, [PatchOp(op = :set, target = :transformer, id = tap_id, field = :tap_pos, value = 8.0)], index)
        @test lbr.tap_ratio == lbr.ratio + 8 * 0.01
        restore!(lnet, lundo)
        @test lbr.tap_ratio == lbr.ratio
      end
    end)() end

    @testset "set angle_deg on the external grid" begin (function ()
      # validate_scenarios took the op and apply! threw "set field angle_deg
      # does not apply to this appliance". The value is the source's
      # u_ref_angle; on the reference unit every bus angle of the island
      # moves by the same amount and the magnitudes stay.
      local fcase, fidx, src, ops, fbase, fwork, k, bus, undo
      fcase = read_scf_json(abspath(joinpath(dirname(@__DIR__), "data", "scf", "feeder3_hardpv_A.scf.json")))
      fidx = ScenarioIndex(fcase)
      src = only(fcase.data.source).id
      ops = [PatchOp(op = :set, target = :external_grid, id = src, field = :angle_deg, value = 10.0)]
      @test validate_scenarios(ScenarioSet(scenarios = [Scenario(name = "angle", ops = ops)]), fidx) isa ScenarioSet
      fbase = Sparlectra.build_net(fcase)
      fwork = Sparlectra.build_net(fcase)
      k = fidx.internal_by_id[src]
      bus = Sparlectra._scf_prosumer_bus(fwork.prosumpsVec[k])
      undo = apply!(fwork, ops, fidx)
      @test fwork.prosumpsVec[k].va_deg == 10.0 && fwork.nodeVec[bus]._va_deg == 10.0
      restore!(fwork, undo)
      @test isequal(fwork.prosumpsVec[k].va_deg, fbase.prosumpsVec[k].va_deg) && isequal(fwork.nodeVec[bus]._va_deg, fbase.nodeVec[bus]._va_deg)
      apply!(fwork, ops, fidx)
      for n in (fbase, fwork)
        redirect_stdout(() -> runpf!(n, 40, 1e-9, 0), devnull)
      end
      # the solved angle comes back through angle(V), so the last bit may differ
      @test isapprox(fwork.nodeVec[bus]._va_deg, 10.0; atol = 1e-9)
      @test all(isapprox(fwork.nodeVec[i]._va_deg - fbase.nodeVec[i]._va_deg, 10.0; atol = 1e-9) for i in eachindex(fbase.nodeVec))
      @test all(isapprox(fwork.nodeVec[i]._vm_pu, fbase.nodeVec[i]._vm_pu; atol = 1e-9) for i in eachindex(fbase.nodeVec))
    end)() end

    @testset "tap patch on a REGULATED transformer is rejected, controller named" begin (function ()
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
    end)() end

    @testset "N-1 expansion equals the existing generators on case14" begin (function ()
      # the expansion names the cases after the file (contingency_names
      # with the index), the generator does the same with that map
      expanded = expand_scenarios(ScenarioSet(mode = :n1_branches), net, index)
      file_names = contingency_names(net, index)
      reference = Sparlectra.generateN1Branches(net; names = file_names)
      @test length(expanded) == length(reference)
      @test [s.name for s in expanded] == [c.name for c in reference]
      @test all(length(s.ops) == 1 && s.ops[1].op === :status && s.ops[1].value == 0.0 for s in expanded)
      gens = expand_scenarios(ScenarioSet(mode = :n1_generators), net, index)
      gref = Sparlectra.generateN1Generators(net; names = file_names)
      @test length(gens) == length(gref)
      @test [s.name for s in gens] == [c.name for c in gref]
      alln1 = expand_scenarios(ScenarioSet(mode = :n1_all), net, index)
      @test length(alln1) == length(reference) + length(gref)
      excluded = expand_scenarios(ScenarioSet(mode = :n1_branches, exclusions = [only(reference[1].elements)]), net, index)
      @test length(excluded) == length(reference) - 1
      # one result-level comparison: runScenarios! on :n1_all writes the same
      # CSV as runContingencies! on the generated cases, both named after
      # the case file
      via_cases = Sparlectra.runContingencies!(net, vcat(reference, gref); names = file_names)
      via_scenarios = runScenarios!(net, ScenarioSet(mode = :n1_all); index = index)
      d = mktempdir()
      a = Sparlectra.writeContingencyResultsCSV(joinpath(d, "cases.csv"), via_cases)
      b = Sparlectra.writeContingencyResultsCSV(joinpath(d, "scenarios.csv"), via_scenarios)
      @test read(a, String) == read(b, String)
      # a general patch scenario runs the same ladder and metrics
      patch_gen = first(id for (id, k) in index.kind_by_id if k === :generator)
      patched = runScenarios!(net, [Scenario(name = "setpoint", ops = [PatchOp(op = :set, target = :generator, id = patch_gen, field = :p, value = 30.0)])]; index = index)
      @test (length(patched), patched[1].converged, patched[1].error) == (1, true, nothing)

      # a three-winding transformer expands to ONE scenario that opens all
      # three legs (through SCF, which carries the star point and the leg
      # ids), and that scenario gives the same result as the N-1 case
      path = tempname() * ".scf.json"
      Sparlectra.exportSCF(three_winding_star_net(); file = path)
      starcase = read_scf_json(path)
      starnet = Sparlectra.build_net(starcase)
      staridx = ScenarioIndex(starcase)
      t3s = only(s for s in expand_scenarios(ScenarioSet(mode = :n1_branches), starnet, staridx) if s.name == "B2_B3_B4")
      @test length(t3s.ops) == 3
      @test all(op -> op.op === :status && op.value == 0.0 && op.target === :transformer, t3s.ops)
      rs = only(redirect_stdout(() -> Sparlectra.runScenarios!(starnet, [t3s]; index = staridx), devnull))
      rc = only(redirect_stdout(() -> Sparlectra.runContingencies!(starnet, [c for c in Sparlectra.generateN1Branches(starnet) if c.kind === :transformer3w]), devnull))
      @test rs.converged && rc.converged
      @test rs.min_vm_pu == rc.min_vm_pu
      @test rs.voltage_violations == rc.voltage_violations
      @test rs.error == rc.error
      # 0.40.0: the legs are resolved when the item is built (member item);
      # the same three ops forced through the patch path (a neutral scale
      # op) give the bitwise same values, and the row lists the legs
      @test rs.kind === :transformer3w && rc.kind === :transformer3w
      @test length(rs.elements) == 3 && rs.elements == rc.elements
      neutral = first(id for (id, k) in staridx.kind_by_id if k in (:load, :generator))
      t3p = Scenario(name = "patched", ops = vcat(t3s.ops, [PatchOp(op = :scale, target = staridx.kind_by_id[neutral], id = neutral, factor = 1.0)]))
      rp = only(redirect_stdout(() -> Sparlectra.runScenarios!(starnet, [t3p]; index = staridx), devnull))
      @test _scenario_values_equal(rs, rp)
    end)() end

    @testset "QSTS vocabulary on a cgmes_series index" begin (function ()
      # QSTS stage 3: the delta vocabulary of the CGMES series union, run
      # on the shipped sp_case60 with a hand-built index of origin
      # :cgmes_series. Every row applies, checks its effect and restores
      # the net BITWISE. `local` is load-bearing: the enclosing test
      # function holds `case`, `net` and `index` of the warmup case.
      local qpath, qcase, base, work, qidx, ids, scfidx, rows, set1, parked, pidx, tmpl, tmpl_ref, undo, mixed, nets, iso_before, applied, edited, eref, patched, plain
      qpath = abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case60.scf.json"))
      qcase = read_scf_json(qpath)
      base = Sparlectra.build_net(qcase)
      work = Sparlectra.build_net(qcase)
      qidx, ids = _qsts_test_index(work)
      scfidx, _ = _qsts_test_index(work; origin = :scf)
      # one tap data record for transformer 70: a table with three ratios
      # and an impedance per step, as a CGMES snapshot carries it
      local tsi = Sparlectra.TapStepImpedance(0.01, 0.1, 0.0, 0.0, 0.0, 0.0, [0.9, 1.0, 1.1], [(1.02, 1.03, 1.0, 1.0), (1.0, 1.0, 1.0, 1.0), (0.98, 0.97, 1.0, 1.0)], Float64[], NTuple{4,Float64}[])
      qidx.tap_data[ids[(:branch, 70)]] = [Sparlectra.BranchTapData(true, 0.9, 1.1, 0.1, [0.9, 1.0, 1.1], false, -30.0, 30.0, 1.25, tsi)]
      set1(group, i, field, value) = PatchOp(op = :set, target = qidx.kind_by_id[ids[(group, i)]], id = ids[(group, i)], field = field, value = value)
      # the case: 1 the reference unit (bus 1), 2 a regulated generator
      # alone at bus 4, 6 the only load at bus 13, 14 a load next to the
      # unregulated generator 28 at bus 36; branch 1 a line, 70 a
      # transformer; shunt 1 at bus 13; two tap controllers
      @test base.prosumpsVec[1].referencePri == 1 && Sparlectra.isGenerator(base.prosumpsVec[2]) && base.nodeVec[4]._nodeType == Sparlectra.PV
      @test !Sparlectra.isGenerator(base.prosumpsVec[6]) && base.nodeVec[13]._pƩLoad == 25.0 && base.branchVec[70].ratio != 0.0
      @test length(Sparlectra._tap_controllers(base)) == 2 && isempty(base.isoNodes)
      parked_state(ps) = ps.proSumptionType == Sparlectra.Consumption && ps.pVal === 0.0 && ps.qVal === 0.0 && !ps.isRegulated &&
                         ps.referencePri === nothing && ps.referencePriority == 0 && ps.participationFactor === nothing && ps.minQ === nothing && ps.maxQ === nothing &&
                         all(f -> getfield(ps, f) === nothing, (:vset_adjust, :vstep_pu, :tap_steps_down, :tap_steps_up, :quController, :puController)) &&
                         !Sparlectra.isGenerator(ps) && !Sparlectra.isRegulating(ps) && !Sparlectra.isSlack(ps)
      ld(n, b) = (something(n.nodeVec[b]._pƩLoad, 0.0), something(n.nodeVec[b]._qƩLoad, 0.0))
      rows = [
        ("bus isolated", [set1(:bus, 13, :isolated, 1.0)],
          n -> n.nodeVec[13]._nodeType == Sparlectra.Isolated && n.nodeVec[13]._vm_pu === 0.0 && n.nodeVec[13]._va_deg === 0.0 && n.isoNodes == [13]),
        ("bus start values", [set1(:bus, 20, :vm_start, 1.031), set1(:bus, 20, :va_start, -4.5)],
          n -> n.nodeVec[20]._vm_pu == 1.031 && n.nodeVec[20]._va_deg == -4.5),
        ("line status only", [set1(:branch, 1, :status, 0.0)],
          n -> (n.branchVec[1].status, n.branchVec[1].from_status, n.branchVec[1].to_status) == (0, 1, 1)),
        ("line terminals", [set1(:branch, 1, :from_status, 0.0), set1(:branch, 1, :to_status, 0.0)],
          n -> (n.branchVec[1].status, n.branchVec[1].from_status, n.branchVec[1].to_status) == (1, 0, 0)),
        ("transformer status and terminals", [set1(:branch, 70, :status, 0.0), set1(:branch, 70, :from_status, 0.0), set1(:branch, 70, :to_status, 0.0)],
          n -> (n.branchVec[70].status, n.branchVec[70].from_status, n.branchVec[70].to_status) == (0, 0, 0)),
        ("transformer ratio", [set1(:branch, 70, :ratio, 1.05)],
          n -> n.branchVec[70].ratio == 1.05 && n.branchVec[70].tap_ratio == 1.05),
        # branch parameters: r and x write the base field too, the terminal
        # shunt values keep the totals their sums, the shift writes angle
        # (read by the flows) and phase_shift_deg (read by the solver)
        ("line impedance", [set1(:branch, 1, :r_pu, 0.07), set1(:branch, 1, :x_pu, -0.02)],
          n -> (br = n.branchVec[1]; (br.r_pu, br.r_base_pu, br.x_pu, br.x_base_pu) == (0.07, 0.07, -0.02, -0.02))),
        ("line charging per terminal", [set1(:branch, 1, :b_from_pu, 0.03), set1(:branch, 1, :b_to_pu, 0.01), set1(:branch, 1, :g_from_pu, 0.001)],
          n -> (br = n.branchVec[1]; (br.b_from_pu, br.b_to_pu, br.b_pu, br.g_from_pu, br.g_to_pu, br.g_pu) == (0.03, 0.01, 0.03 + 0.01, 0.001, 0.0, 0.001))),
        ("transformer impedance, magnetizing, shift", [set1(:branch, 70, :x_pu, 0.1), set1(:branch, 70, :g_to_pu, 0.002), set1(:branch, 70, :b_to_pu, -0.004), set1(:branch, 70, :shift_deg, -3.5)],
          # the from arm keeps its nameplate magnetizing (sp_case60 since 0.50.0), so the
          # sums are that arm plus the new to arm
          n -> (br = n.branchVec[70]; (br.x_pu, br.x_base_pu, br.r_pu, br.g_to_pu, br.b_to_pu, br.g_pu, br.b_pu, br.angle, br.phase_shift_deg) == (0.1, 0.1, base.branchVec[70].r_pu, 0.002, -0.004, base.branchVec[70].g_from_pu + 0.002, base.branchVec[70].b_from_pu + -0.004, -3.5, -3.5))),
        ("transformer tap data", [set1(:branch, 70, :tap_data, 1.0)],
          n -> (br = n.branchVec[70]; (br.tap_min, br.tap_max, br.tap_step, br.tap_ratio_points, br.has_phase_tap, br.tap_step_impedance === tsi) == (0.9, 1.1, 0.1, [0.9, 1.0, 1.1], false, true))),
        ("shunt conductance", [set1(:shunt, 1, :g_pu, 0.02)],
          n -> n.shuntVec[1].y_pu_shunt == ComplexF64(0.02, imag(base.shuntVec[1].y_pu_shunt)) && n.shuntVec[1].G_shunt == 0.02),
        ("generator limits, regulation, reference", [set1(:prosumer, 2, :qmin, -50.0), set1(:prosumer, 2, :qmax, 60.0), set1(:prosumer, 2, :regulated, 0.0),
            set1(:prosumer, 2, :reference, 1.0), set1(:prosumer, 2, :reference_priority, 2.0), set1(:prosumer, 2, :participation, 0.5)],
          n -> (ps = n.prosumpsVec[2]; ps.minQ == -50.0 && ps.maxQ == 60.0 && !ps.isRegulated && ps.referencePri == 4 && ps.referencePriority == 2 && ps.participationFactor == 0.5 && 4 in n.slackVec)),
        ("generator limits to none", [set1(:prosumer, 2, :qmin, nothing), set1(:prosumer, 2, :qmax, nothing), set1(:prosumer, 2, :participation, nothing)],
          n -> n.prosumpsVec[2].minQ === nothing && n.prosumpsVec[2].maxQ === nothing),
        ("external grid gives up the reference", [set1(:prosumer, 1, :reference, 0.0)],
          n -> n.prosumpsVec[1].referencePri === nothing && isempty(n.slackVec)),
        ("tap controller fields", [set1(:tap, 1, :enabled, 0.0), set1(:tap, 1, :target_vm_pu, 1.02), set1(:tap, 1, :deadband_vm_pu, 0.005), set1(:tap, 1, :p_target_mw, 25.0), set1(:tap, 1, :deadband_p_mw, 2.0), set1(:tap, 2, :target_vm_pu, nothing)],
          n -> (c = Sparlectra._tap_controllers(n); !c[1].enabled && c[1].target_vm_pu == 1.02 && c[1].deadband_vm_pu == 0.005 && c[1].p_target_mw == 25.0 && c[1].deadband_p_mw == 2.0 && c[2].target_vm_pu === nothing)),
        ("park a load", [set1(:prosumer, 6, :in_service, 0.0)],
          n -> parked_state(n.prosumpsVec[6]) && ld(n, 13) == (0.0, 0.0) && length(n.prosumpsVec) == length(base.prosumpsVec)),
        ("park a load next to a generator", [set1(:prosumer, 14, :in_service, 0.0)],
          n -> parked_state(n.prosumpsVec[14]) && ld(n, 36) == ld(base, 36) .- (3.0, 0.7) && n.nodeVec[36]._pƩGen == 4.0),
        ("park a generator", [set1(:prosumer, 2, :in_service, 0.0)],
          n -> parked_state(n.prosumpsVec[2]) && n.nodeVec[4]._nodeType == Sparlectra.PQ && n.nodeVec[4]._pƩGen === nothing && ld(n, 4) == ld(base, 4)),
        ("park the external grid", [set1(:prosumer, 1, :in_service, 0.0)],
          n -> parked_state(n.prosumpsVec[1]) && isempty(n.slackVec)),
      ]
      for (label, ops, check) in rows
        @test validate_scenarios(ScenarioSet(scenarios = [Scenario(name = label, ops = ops)]), qidx) isa ScenarioSet
        undo = apply!(work, ops, qidx)
        @test (label, check(work)) == (label, true)
        restore!(work, undo)
        @test (label, _qsts_bitwise_diffs(base, work)) == (label, String[])
        @test (label, scf_roundtrip_field_diffs(base, work)) == (label, String[])
        # an SCF index rejects every op of the vocabulary with the op named,
        # and apply! on it keeps today's rejection
        for op in ops
          err = try
            validate_scenarios(ScenarioSet(scenarios = [Scenario(name = label, ops = [op])]), scfidx)
            ""
          catch e
            e isa ArgumentError || rethrow()
            sprint(showerror, e)
          end
          @test (label, op.field, occursin(string("set ", op.field, " is QSTS vocabulary"), err)) == (label, op.field, true)
          @test_throws ArgumentError apply!(work, [op], scfidx)
        end
      end
      @test _qsts_bitwise_diffs(base, work) == String[]

      # the bus op writes a NEW isoNodes vector: a later markIsolatedBuses!
      # refills that one in place, and restore! puts the old one back
      iso_before = work.isoNodes
      undo = apply!(work, [set1(:bus, 13, :isolated, 1.0)], qidx)
      @test work.isoNodes !== iso_before
      Sparlectra.markIsolatedBuses!(net = work)
      # the op isolates a bus that keeps its branches: the marking leaves
      # its type alone and only rebuilds isoNodes
      local b13 = qidx.internal_by_id[ids[(:bus, 13)]]
      local brs13 = [j for (j, br) in enumerate(work.branchVec) if br.status == 1 && (br.fromBus == b13 || br.toBus == b13)]
      @test !isempty(brs13) && Sparlectra.isIsolated(work.nodeVec[b13]) && !(b13 in work.isoNodes)
      restore!(work, undo)
      @test work.isoNodes === iso_before && isempty(iso_before)
      @test _qsts_bitwise_diffs(base, work) == String[]
      # a bus the marking isolated itself is typed back once a branch reaches
      # it again, unless an op isolated it in between
      local cut = deepcopy(work)
      foreach(j -> setBranchStatus!(cut.branchVec[j], false), brs13)
      @test b13 in Sparlectra.markIsolatedBuses!(net = cut)
      local undo13 = apply!(cut, [set1(:bus, 13, :isolated, 1.0)], qidx)
      foreach(j -> setBranchStatus!(cut.branchVec[j], true), brs13)
      Sparlectra.markIsolatedBuses!(net = cut)
      @test Sparlectra.isIsolated(cut.nodeVec[b13])
      restore!(cut, undo13)
      Sparlectra.markIsolatedBuses!(net = cut)
      @test !Sparlectra.isIsolated(cut.nodeVec[b13]) && !(b13 in cut.isoNodes)

      # validation of the vocabulary on the QSTS index itself
      reject(ops) = try
        validate_scenarios(ScenarioSet(scenarios = [Scenario(name = "bad", ops = ops)]), qidx)
        ""
      catch e
        e isa ArgumentError || rethrow()
        sprint(showerror, e)
      end
      @test occursin("set ratio applies to a transformer, not to a branch", reject([set1(:branch, 1, :ratio, 1.0)]))
      @test occursin("value must be 0 or 1", reject([set1(:bus, 13, :isolated, 0.5)]))
      @test occursin("whole number", reject([set1(:prosumer, 2, :reference_priority, 1.5)]))
      @test occursin("parked values of id", reject([set1(:prosumer, 6, :in_service, 1.0)]))
      @test occursin("second set qmin", reject([set1(:prosumer, 2, :qmin, -1.0), set1(:prosumer, 2, :qmin, -2.0)]))
      @test occursin("set x_pu: needs a finite value", reject([set1(:branch, 1, :x_pu, Inf)]))
      @test occursin("needs tap data record 2 of id", reject([set1(:branch, 70, :tap_data, 2.0)]))
      @test occursin("record number >= 1", reject([set1(:branch, 70, :tap_data, 0.0)]))
      @test_throws ArgumentError ScenarioIndex(qidx.kind_by_id, qidx.internal_by_id, qidx.name_by_id, qidx.id_by_mrid, qidx.parked_values, :matpower)
      # apply! rejects r = x = 0 and an r/x write on a typed tap model
      # before any write, so the net stays as it was
      applied(ops) = try
        apply!(work, ops, qidx)
        ""
      catch e
        e isa ArgumentError || rethrow()
        sprint(showerror, e)
      end
      @test occursin("r_pu and x_pu would both be 0", applied([set1(:branch, 1, :r_pu, 0.0), set1(:branch, 1, :x_pu, 0.0)]))
      # shift_deg passes the validation on a :branch (a phase-stamped
      # EquivalentBranch has ratio 1), on a line apply! rejects it
      @test reject([set1(:branch, 1, :shift_deg, 1.0)]) == ""
      @test occursin("is a line; shift_deg applies to a transformer", applied([set1(:branch, 1, :shift_deg, 1.0)]))
      work.branchVec[70].taps_derived = true
      @test occursin("typed tap model", applied([set1(:branch, 70, :x_pu, 0.1)]))
      work.branchVec[70].taps_derived = false
      # a single r or x op, the common case of a series delta, is judged
      # against the net value of the other field: line 1 (r 0.06, x 0.3)
      # takes x = 0 alone and r = 0 alone, with x already 0 it rejects r = 0
      for fz in (:x_pu, :r_pu)
        undo = apply!(work, [set1(:branch, 1, fz, 0.0)], qidx)
        @test (fz, getfield(work.branchVec[1], fz)) == (fz, 0.0)
        restore!(work, undo)
        @test (fz, _qsts_bitwise_diffs(base, work)) == (fz, String[])
      end
      work.branchVec[1].x_pu = 0.0
      @test occursin("r_pu and x_pu would both be 0", applied([set1(:branch, 1, :r_pu, 0.0)]))
      @test work.branchVec[1].r_pu == base.branchVec[1].r_pu
      work.branchVec[1].x_pu = base.branchVec[1].x_pu
      @test _qsts_bitwise_diffs(base, work) == String[]

      # unpark: the template holds generator 2 and load 6 parked and bus 20
      # isolated, the index their values from before; the delta brings all
      # three back, restore! returns the template bitwise
      tmpl = Sparlectra.build_net(qcase)
      parked = Dict{Int,Sparlectra.ProSumer}(ids[(:prosumer, i)] => deepcopy(tmpl.prosumpsVec[i]) for i in (2, 6))
      pidx, _ = _qsts_test_index(tmpl; parked = parked)
      apply!(tmpl, [set1(:prosumer, 2, :in_service, 0.0), set1(:prosumer, 6, :in_service, 0.0), set1(:bus, 20, :isolated, 1.0)], pidx)
      @test parked_state(tmpl.prosumpsVec[2]) && parked_state(tmpl.prosumpsVec[6]) && tmpl.isoNodes == [20]
      tmpl_ref = deepcopy(tmpl)
      undo = apply!(tmpl, [set1(:prosumer, 2, :in_service, 1.0), set1(:prosumer, 6, :in_service, 1.0), set1(:bus, 20, :isolated, 0.0),
                           set1(:bus, 20, :vm_start, base.nodeVec[20]._vm_pu), set1(:bus, 20, :va_start, base.nodeVec[20]._va_deg)], pidx)
      for i in (2, 6), f in (:proSumptionType, :pVal, :qVal, :isRegulated, :referencePri, :referencePriority, :participationFactor, :minQ, :maxQ,
                             :vset_adjust, :vstep_pu, :tap_steps_down, :tap_steps_up, :quController, :puController)
        @test (i, f, isequal(getfield(tmpl.prosumpsVec[i], f), getfield(base.prosumpsVec[i], f))) == (i, f, true)
      end
      @test all(isapprox.(ld(tmpl, 13), ld(base, 13); atol = 1e-12))
      @test [tmpl.nodeVec[b]._nodeType for b in (4, 13, 20)] == [base.nodeVec[b]._nodeType for b in (4, 13, 20)]
      @test isempty(tmpl.isoNodes) && tmpl.nodeVec[4]._pƩGen == base.nodeVec[4]._pƩGen
      restore!(tmpl, undo)
      @test _qsts_bitwise_diffs(tmpl_ref, tmpl) == String[]
      @test scf_roundtrip_field_diffs(tmpl_ref, tmpl) == String[]

      # the order of the ops in a delta does not change the result: QSTS
      # ops mixed with SCF ops (setpoints, scaling, b_pu, a link, a
      # generator outage that deletes), a scale and a set q on one load,
      # a generator vm_pu and a start value on the same bus
      mixed = [set1(:prosumer, 6, :in_service, 0.0), set1(:prosumer, 14, :in_service, 0.0), set1(:prosumer, 28, :in_service, 0.0),
               set1(:prosumer, 2, :qmin, -40.0), set1(:prosumer, 2, :vm_pu, 1.02), set1(:bus, 4, :vm_start, 1.015),
               set1(:bus, 13, :isolated, 1.0), set1(:bus, 20, :vm_start, 0.99), set1(:branch, 1, :status, 0.0), set1(:branch, 70, :ratio, 0.98),
               set1(:branch, 1, :b_from_pu, 0.02), set1(:branch, 1, :b_to_pu, 0.01), set1(:branch, 70, :x_pu, 0.1), set1(:branch, 70, :shift_deg, 1.5),
               set1(:shunt, 1, :g_pu, 0.01), set1(:shunt, 1, :b_pu, 0.2), set1(:tap, 1, :enabled, 0.0), set1(:tap, 2, :deadband_vm_pu, 0.002),
               PatchOp(op = :scale, target = :load, id = ids[(:prosumer, 7)], factor = 1.1), set1(:prosumer, 7, :q, 2.5), set1(:prosumer, 8, :p, 17.0),
               PatchOp(op = :status, target = :link, id = ids[(:link, 1)], value = 0.0), PatchOp(op = :status, target = :generator, id = ids[(:prosumer, 3)], value = 0.0)]
      @test validate_scenarios(ScenarioSet(scenarios = [Scenario(name = "mixed", ops = mixed)]), qidx) isa ScenarioSet
      nets = [Sparlectra.build_net(qcase) for _ in 1:3]
      apply!(nets[1], mixed, qidx)
      apply!(nets[2], reverse(mixed), qidx)
      apply!(nets[3], vcat(mixed[2:2:end], mixed[1:2:end]), qidx)
      @test length(nets[1].prosumpsVec) == length(base.prosumpsVec) - 1
      @test nets[1].nodeVec[4]._vm_pu == 1.015 && nets[1].nodeVec[13]._vm_pu === 0.0
      @test _qsts_bitwise_diffs(nets[1], nets[2]) == String[]
      @test _qsts_bitwise_diffs(nets[1], nets[3]) == String[]

      # the solver reads what the branch ops write: line 1 with other r, x
      # and charging, patched, equals the case built with those values
      # bitwise, and so does its power flow, which differs from the
      # unpatched one
      edited = read_scf_json(qpath)
      let eidx = ScenarioIndex(edited), row = only(r for r in edited.data.line if get(eidx.internal_by_id, r.id, 0) == 1)
        row.r1 *= 1.2
        row.x1 *= 1.5
        row.c1 = 5.0e-7
      end
      eref = Sparlectra.build_net(edited)
      patched = Sparlectra.build_net(qcase)
      apply!(patched, [set1(:branch, 1, f, getfield(eref.branchVec[1], f)) for f in (:r_pu, :x_pu, :b_from_pu, :b_to_pu)], qidx)
      @test eref.branchVec[1].b_pu > 0.0 && _qsts_bitwise_diffs(eref, patched) == String[]
      plain = Sparlectra.build_net(qcase)
      @test [last(redirect_stdout(() -> runpf!(n, 40, 1e-10, 0), devnull)) for n in (eref, patched, plain)] == [0, 0, 0]
      @test isequal([(nd._vm_pu, nd._va_deg) for nd in patched.nodeVec], [(nd._vm_pu, nd._va_deg) for nd in eref.nodeVec])
      @test maximum(abs(patched.nodeVec[b]._vm_pu - plain.nodeVec[b]._vm_pu) for b in eachindex(plain.nodeVec)) > 1e-6
    end)() end

  end)() end
  return nothing
end

function run_scenario_engine_extended_tests()
  @testset "scenario engine sp_case60 screening acceptance" begin (function ()
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
    # the :flag surface: rows in the order of the full run; a screened row
    # carries the estimate, any other row the full run row
    @test flagged isa Vector{ScenarioResult}
    @test [r.name for r in flagged] == [r.name for r in full]
    for r in flagged
      if r.screened
        @test r.start_used === :screen
        @test r.screening_estimate !== nothing
        @test r.error === nothing
      else
        @test r.result isa Sparlectra.ContingencyResult
        @test r.start_used !== :screen
      end
    end
    # the screening CSV appends exactly the two screening columns
    header = run_csv(Sparlectra.writeContingencyResultsCSV(joinpath(mktempdir(), "flag.csv"), flagged)).header
    @test header[(end-1):end] == ["screened", "screening_estimate"]
    # :only never runs a full solve where an estimate exists; islanding
    # cases still carry real solves
    only_res = Sparlectra.runContingencies!(net, cases; screening_mode = :only, screening_margin_pct = 10.0)
    @test any(r -> r.screened, only_res)
    @test all(r -> r.screened || r.start_used !== :screen, only_res)
    ds_full = Sparlectra.runContingencies!(net, cases; distributed_slack_enabled = true)
    ds_flagged = Sparlectra.runContingencies!(net, cases; screening_mode = :flag, screening_margin_pct = 10.0, distributed_slack_enabled = true)
    ds_violating = Set(r.name for r in ds_full if !isempty(r.overloads) || !isempty(r.voltage_violations) || !r.converged)
    @test sort(collect(intersect(ds_violating, Set(r.name for r in ds_flagged if r.screened)))) == String[]
    @test count(r -> r.screened, ds_flagged) > 0
    # warm active set (0.30.2, default since 0.30.3): outages start from the
    # base case's PV/PQ state; on the shipped grid the violation list is the
    # same as from the file's PV/PQ state (the cold run, named explicitly
    # since the default changed)
    warm = Sparlectra.runContingencies!(net, cases; warm_active_set = true, warm_note = Ref(""))
    cold = Sparlectra.runContingencies!(net, cases; warm_active_set = false)
    viol_list(rs) = [(r.name, r.converged, sort([o.name for o in r.overloads]), sort(r.voltage_violations)) for r in rs]
    @test viol_list(warm) == viol_list(cold)
    # cold check (0.30.3): with a margin that covers every outage, each warm
    # result is solved cold as well and the less favourable one counts, so
    # no outage ends more favourable than in the cold run
    checked = Sparlectra.runContingencies!(net, cases; warm_active_set = true, warm_cold_check = true, warm_cold_check_margin_pu = 1.0, warm_note = Ref(""))
    @test all(!Sparlectra._contingency_less_favourable(c, k) for (c, k) in zip(cold, checked))
    # the order of "less favourable": a failed result never wins, then more
    # violations, then the lower lowest voltage, then the higher loading
    row(conv, nviol, vmin, load) = Sparlectra.ContingencyResult("x", :explicit, String[], 1.0, conv, 3, :warm, 1.05, vmin, load, 0.0, Sparlectra.OverloadRecord[], ["b$(i)" for i in 1:nviol], 1, 0.0, nothing)
    # the switch (contingency.warm_cold_check, default off): with the same
    # margin a tight warm result is checked when on and never by default
    check_on = Sparlectra.ScenarioEngine(net; warm_cold_check = true, warm_cold_check_margin_pu = 1.0, warm_note = Ref(""))
    check_off = Sparlectra.ScenarioEngine(net; warm_cold_check_margin_pu = 1.0, warm_note = Ref(""))
    @test Sparlectra._warm_needs_cold_check(check_on, row(true, 0, 0.95, 50.0))
    @test !Sparlectra._warm_needs_cold_check(check_off, row(true, 0, 0.95, 50.0))
    @test !Sparlectra._warm_needs_cold_check(check_off, row(true, 2, 0.85, 120.0))
    # iteration limits (contingency.max_iter, 0.31.0): maxIte caps every
    # outage solve, base_maxIte the base case; from start angles 0 the base
    # case needs more than one step
    flat_angles = deepcopy(net)
    foreach(n -> (n._va_deg = 0.0), flat_angles.nodeVec)
    quietly(f) = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        redirect_stdout(f, devnull)
    end
    @test !quietly(() -> Sparlectra.ScenarioEngine(flat_angles; maxIte = 1, warm_note = Ref(""))).base_converged
    split_limits = quietly(() -> Sparlectra.ScenarioEngine(flat_angles; maxIte = 1, base_maxIte = 30, warm_note = Ref("")))
    @test split_limits.base_converged && (split_limits.maxIte, split_limits.base_maxIte) == (1, 30)
    capped = quietly(() -> Sparlectra.runContingencies!(net, cases[1:10]; maxIte = 1, base_maxIte = 30, warm_note = Ref("")))
    @test all(r -> !r.converged && r.iterations <= 1, capped)
    @test Sparlectra._contingency_less_favourable(row(true, 1, 0.97, 50.0), row(true, 0, 0.92, 90.0))
    @test Sparlectra._contingency_less_favourable(row(true, 0, 0.94, 50.0), row(true, 0, 0.95, 90.0))
    @test Sparlectra._contingency_less_favourable(row(true, 0, 0.95, 91.0), row(true, 0, 0.95, 90.0))
    @test !Sparlectra._contingency_less_favourable(row(false, 0, NaN, NaN), row(true, 0, 0.95, 90.0))
    # rounding noise of two solves of the same solution decides nothing
    @test !Sparlectra._contingency_less_favourable(row(true, 0, 0.95 - 1e-9, 90.0 + 1e-6), row(true, 0, 0.95, 90.0))
    @test Sparlectra._contingency_less_favourable(row(true, 0, 0.95, 90.0), row(false, 0, NaN, NaN))
    # change-based loading flag (0.30.2): a branch the base case already
    # loads at 95 percent used to flag every outage (the network-wide
    # maximum stayed above 100 - margin); now only outages that change it
    # enough are flagged, still without a false negative
    eng = Sparlectra.ScenarioEngine(net)
    top_k = argmax(k -> isnan(eng.base_loadings[k]) ? -Inf : eng.base_loadings[k], eachindex(net.branchVec))
    top = net.branchVec[top_k]
    top_loading = eng.base_loadings[top_k]
    top.sn_MVA = top.sn_MVA * top_loading / 95.0
    hot_full = Sparlectra.runContingencies!(net, cases)
    hot_flagged = Sparlectra.runContingencies!(net, cases; screening_mode = :flag, screening_margin_pct = 10.0)
    hot_violating = Set(r.name for r in hot_full if !isempty(r.overloads) || !isempty(r.voltage_violations) || !r.converged)
    @test sort(collect(intersect(hot_violating, Set(r.name for r in hot_flagged if r.screened)))) == String[]
    @test count(r -> r.screened, hot_flagged) > 0
    println("      scenario engine sp_case60: RAN (", length(cases), " cases, ", count(r -> r.screened, flagged), " screened, 0 false negatives, ds ", count(r -> r.screened, ds_flagged), " screened, one branch at 95 percent ", count(r -> r.screened, hot_flagged), " screened)")
  end)() end

  @testset "scenario engine sp_case300 screening acceptance" begin (function ()
    # the operated-grid acceptance case per the test-network rule
    # (the test-network rule): case300 judges screening QUALITY, no false negatives at
    # margin 10 with the 0.005 pu trust gate (its outage 57-63 forced the
    # gate: buses 63/64/526 collapse to 0.84 pu behind a 0.0104 pu one-step
    # residual); runs on the shipped data/mpower/sp_case300.m (synthetic
    # 300-bus operated grid, 411 branches, 69 generators)
    case_path = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "sp_case300.m"))
    @test isfile(case_path)
    begin
      net = Sparlectra.createNetFromMatPowerFile(filename = case_path, flatstart = false, enable_pq_gen_controllers = true, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)
      cases = vcat(Sparlectra.generateN1Branches(net), Sparlectra.generateN1Generators(net))
      full = Sparlectra.runContingencies!(net, cases)
      flagged = Sparlectra.runContingencies!(net, cases; screening_mode = :flag, screening_margin_pct = 10.0)
      violating = Set(r.name for r in full if !isempty(r.overloads) || !isempty(r.voltage_violations) || !r.converged)
      false_negatives = intersect(violating, Set(r.name for r in flagged if r.screened))
      @test sort(collect(false_negatives)) == String[]
      screened = count(r -> r.screened, flagged)
      fp = count(r -> !r.screened && !(r.name in violating), flagged)
      println("      scenario engine sp_case300: RAN (", length(cases), " cases, ", screened, " screened, ", fp, " false positives, 0 false negatives)")
    end
  end)() end

  return nothing
end

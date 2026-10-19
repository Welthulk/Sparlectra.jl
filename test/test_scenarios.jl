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

# Fixture comparison for contingency result CSVs. Byte identity holds on the
# platform that generated a fixture but not across BLAS/UMFPACK builds: on
# Windows the last one or two significant digits of min_vm_pu and the
# loadings differ, and a CRLF checkout adds a byte per line. So the shape is
# compared exactly (header, row count, every non-numeric field: names,
# converged, start_used, error texts, violation lists) and numeric fields
# within rtol 1e-9, NaN matching NaN. The writer's format is still guarded:
# a changed delimiter, column order or quoting fails the shape comparison.
function _csv_fixture_fields(line::AbstractString)
  fields = String[]
  buf = IOBuffer()
  quoted = false
  i = firstindex(line)
  while i <= lastindex(line)
    c = line[i]
    if quoted
      if c == '"'
        nxt = nextind(line, i)
        if nxt <= lastindex(line) && line[nxt] == '"'
          write(buf, '"')
          i = nxt
        else
          quoted = false
        end
      else
        write(buf, c)
      end
    elseif c == '"'
      quoted = true
    elseif c == ','
      push!(fields, String(take!(buf)))
    else
      write(buf, c)
    end
    i = nextind(line, i)
  end
  push!(fields, String(take!(buf)))
  return fields
end

function _csv_fixture_matches(out_path::AbstractString, fixture_path::AbstractString; rtol::Float64 = 1e-9)::Bool
  a = split(replace(read(out_path, String), "\r\n" => "\n"), '\n'; keepempty = false)
  b = split(replace(read(fixture_path, String), "\r\n" => "\n"), '\n'; keepempty = false)
  if length(a) != length(b)
    println("      csv fixture: $(length(a)) lines written, $(length(b)) in the fixture")
    return false
  end
  for (lineno, (la, lb)) in enumerate(zip(a, b))
    fa = _csv_fixture_fields(la)
    fb = _csv_fixture_fields(lb)
    if length(fa) != length(fb)
      println("      csv fixture line $(lineno): $(length(fa)) fields written, $(length(fb)) in the fixture")
      return false
    end
    for (col, (va, vb)) in enumerate(zip(fa, fb))
      xa = tryparse(Float64, va)
      xb = tryparse(Float64, vb)
      same = if xa !== nothing && xb !== nothing
        (isnan(xa) && isnan(xb)) || isapprox(xa, xb; rtol = rtol, atol = 1e-12)
      else
        va == vb
      end
      if !same
        println("      csv fixture line $(lineno) column $(col): written \"$(va)\", fixture \"$(vb)\"")
        return false
      end
    end
  end
  return true
end

# every value of two result rows except the label (name, kind, elements)
_scenario_values_equal(x, y) = all(f -> f in (:name, :kind, :elements) || isequal(getproperty(x, f), getproperty(y, f)), fieldnames(Sparlectra.ContingencyResult))

function _scenario_case14()
  cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
  mpc = Sparlectra.MatpowerIO.read_case(abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m")); legacy_compat = true)
  case = Sparlectra.convert_case(MatpowerAdapter(), mpc, Sparlectra.matpower_adapter_options(cfg))
  net = Sparlectra.build_net(case; config = cfg)
  return case, net
end

function run_scenario_patch_tests()
  @testset "scenario patch model" begin (function ()
    case, net = _scenario_case14()
    index = ScenarioIndex(case)
    branch_id = first(r.id for r in case.data.line)
    trafo_id = first(r.id for r in case.data.generic_branch)
    gen_id = first(r.id for r in case.data.sym_gen)
    load_id = first(r.id for r in case.data.sym_load)

    @testset "validation accepts the patch surface" begin (function ()
      good = ScenarioSet(
        scenarios = [
          Scenario(name = "line out", ops = [PatchOp(op = :status, target = :branch, id = branch_id, value = 0.0)]),
          Scenario(name = "gen setpoint", ops = [PatchOp(op = :set, target = :generator, id = gen_id, field = :p, value = 25.0)]),
          Scenario(name = "load scale", weight = 2.0, ops = [PatchOp(op = :scale, target = :load, id = load_id, factor = 1.2)]),
          Scenario(name = "tap step", ops = [PatchOp(op = :set, target = :transformer, id = trafo_id, field = :tap_pos, value = 2.0)]),
        ],
      )
      @test validate_scenarios(good, index) === good
    end)() end

    @testset "rejections name the scenario and the op index" begin (function ()
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
        undo = apply!(net, [PatchOp(op = :set, target = :transformer, id = tap_id, field = :tap_pos, value = 2.0)], index)
        @test Sparlectra.calcBranchRatio(net.branchVec[tapped]) != before_ratio
        restore!(net, undo)
        @test Sparlectra.calcBranchRatio(net.branchVec[tapped]) == before_ratio
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

    @testset "engine: runScenarios! N-1 equals runContingencies!" begin (function ()
      # both entries name the rows after the case file when given its map
      file_names = contingency_names(net, index)
      cases14 = vcat(Sparlectra.generateN1Branches(net; names = file_names), Sparlectra.generateN1Generators(net; names = file_names))
      via_cases = Sparlectra.runContingencies!(net, cases14; names = file_names)
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

  end)() end
  return nothing
end

function run_scenario_engine_extended_tests()
  @testset "scenario engine sp_case118 fixture" begin (function ()
    # THE step-3 gate: the engine's N-1 result CSV on case118 (all
    # branches, all generators) is byte identical to the CSV the
    # pre-engine per-case-deepcopy implementation produced at 7438b6e
    # (regenerated at e93eafc, test/fixtures/contingency_case118_n1_e93eafc.csv,
    # for issue #376: the writer's default format changed from a hardcoded
    # semicolon+dot-decimal mix to the run-wide "technical" format, comma
    # delimiter and dot decimal with no thousands grouping - a values-only
    # regenerate, not a computation change). The net is built with the
    # pinned packaged-default import options the fixture generation used.
    # since sp_case118 the gate is a REGRESSION fixture of the engine on the
    # shipped synthetic 118-bus case (same role as the sp_case60 fixture
    # below, 240 rows instead of 82), no longer the external pre-engine
    # oracle on the downloaded case118. The fixture is created by the first
    # run that finds it missing and printed LOUDLY so it gets committed; from
    # then on every run compares against it. Regenerated at 0.30.0 (polar
    # Newton update, template Q-limit defaults): 3 of 240 contingencies end
    # in another active set, the rest differ in iteration counts only.
    # Regenerated at 0.30.2 (island runs read the reactive limits of their
    # own buses): exactly the 3 outages that split the net (109-105,
    # 114-93, 94-107) changed, in voltages and loading, 94-107 also in
    # iterations (12 to 10); every non-islanding row is unchanged.
    # Regenerated at 0.30.3 (warm active set default): 215 of 240 rows
    # changed in the iteration count only (the warm start saves steps), the
    # values in the last digits; one real change, B_ACL_345_68_81 (min Vm
    # 0.98014 against 0.98070 cold, 240.3 against 239.2 percent). The
    # optional cold check would keep this warm result as well (it is the
    # less favourable one), so the fixture holds with the check on or off.
    # Regenerated for #475 (joint extra release of twice-clamped machines on
    # a converged iterate): one cell, B_ACL_345_68_81 iterations 16 to 22;
    # its group release does not hold and is undone, every value stays.
    # Regenerated at 0.40.0 (#355): the columns kind and elements after the
    # name; every other cell unchanged (compared column by column).
    # Regenerated at 0.50.0: the three islanding outages (109-105, 114-93,
    # 94-107) converge in 8 instead of 10 iterations, their values move in
    # the 12th digit; every other row is unchanged.
    fixture = abspath(joinpath(@__DIR__, "fixtures", "contingency_sp_case118_n1.csv"))
    case_path = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "sp_case118.m"))
    @test isfile(case_path)
    begin
      net = Sparlectra.createNetFromMatPowerFile(filename = case_path, flatstart = false, enable_pq_gen_controllers = true, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)
      cases = vcat(Sparlectra.generateN1Branches(net), Sparlectra.generateN1Generators(net))
      @test length(cases) == 240
      results = Sparlectra.runContingencies!(net, cases)
      out = Sparlectra.writeContingencyResultsCSV(joinpath(mktempdir(), "sp_case118_n1.csv"), results)
      if !isfile(fixture)
        cp(out, fixture)
        println("      scenario engine sp_case118 fixture: CREATED ", fixture, " (commit it; this run compared nothing)")
      end
      @test _csv_fixture_matches(out, fixture)

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
      println("      scenario engine sp_case118 screening: ", length(flagged) - full_runs, " screened, ", full_runs, " full runs (", round(100 * full_runs / length(flagged); digits = 1), " %), ", fp, " false positives, 0 false negatives")

      # the same acceptance with DISTRIBUTED SLACK enabled; the
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
      println("      scenario engine sp_case118 distributed-slack screening: ", ds_screened, " of ", length(cases), " screened (", round(100 * ds_screened / length(cases); digits = 1), " %), :flag ", round(t_ds.time; digits = 2), " s, 0 false negatives")
    end
  end)() end
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

  @testset "scenario engine sp_case60 CSV byte fixture" begin (function ()
    # load_fixture_net: the engine's N-1 result CSV on the shipped operated
    # grid is byte identical to the tracked fixture, generated at bf0b076
    # (the place-name demo cases) and regenerated at e93eafc for issue #376
    # (the writer's default format changed from a hardcoded semicolon+dot mix
    # to the run-wide "technical" format; a values-only regenerate, see the
    # case118 gate's comment above). Unlike the case118 gate above (an
    # INDEPENDENT oracle from the pre-engine implementation, which stays
    # cache-gated), this fixture guards regressions of the engine's values
    # and CSV format on every install without a download; sp_case60 keeps
    # the tracked file small (10k, the 188 variant carried 27-bus violation
    # lists per row).
    # Regenerated at 0.30.0 for the polar Newton update and the template's
    # Q-limit defaults (iteration counts one lower on most rows, end states
    # at the 1e-8 level): like the case118 gate, the first run that finds
    # the fixture missing creates it from the engine's own CSV and says so.
    # Regenerated at 0.40.0 (#355): the columns kind and elements, the case
    # and overload names from the case file (contingency_names with the
    # typed case), every numeric cell unchanged.
    # Changed at 0.41.0 (link groups): the row of Neubach_Waldsee#29 no
    # longer reports Neubach_110b cut off; the bus hangs on Neubach_110a
    # through the coupler and the solve feeds it. Every numeric cell
    # unchanged.
    # Regenerated at 0.50.0: demo case rebuilt, transformer magnetizing from
    # nameplate since 0.32.1. Only min_vm_pu, max_vm_pu and
    # max_branch_loading_pct move (at most 7.7e-4 pu and 0.077 percent;
    # first row min_vm_pu old 0.921871975379177, new 0.921540854811351);
    # iterations, islands, overloads and violations unchanged.
    fixture = abspath(joinpath(@__DIR__, "fixtures", "contingency_sp_case60_n1.csv"))
    scase = read_scf_json(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case60.scf.json")))
    net = Sparlectra.build_net(scase)
    file_names = contingency_names(net, scase)
    cases = vcat(Sparlectra.generateN1Branches(net; names = file_names), Sparlectra.generateN1Generators(net; names = file_names))
    results = Sparlectra.runContingencies!(net, cases; names = file_names)
    out = Sparlectra.writeContingencyResultsCSV(joinpath(mktempdir(), "sp60_n1.csv"), results)
    if !isfile(fixture)
      cp(out, fixture)
      println("      scenario engine sp_case60 fixture: CREATED ", fixture, " (commit it; this run compared nothing)")
    end
    @test isfile(fixture)
    @test _csv_fixture_matches(out, fixture)
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

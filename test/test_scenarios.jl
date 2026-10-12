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

    @testset "the removed contingencies block fails with the way out" begin (function ()
      # 0.40.0: the legacy block is no longer read; a file that still
      # carries it names the scenarios block and the editor instead of
      # running something else
      path = joinpath(mktempdir(), "legacy.scf.json")
      write_scf_json(case, path)
      root = Sparlectra.scf_json_parse(read(path, String))
      root["sparlectra"]["contingencies"] = Dict{String,Any}("mode" => "explicit", "cases" => Any[Dict{String,Any}("name" => "old style", "outages" => Any[Dict{String,Any}("component" => branch_id)])])
      write(path, Sparlectra.scf_json_string(root))
      err = try
        read_scf_json(path)
        nothing
      catch e
        e
      end
      @test err isa ArgumentError
      @test occursin("sparlectra.scenarios", err.msg) && occursin("scenario editor", err.msg)
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
      end
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
        @test ("warmup_casePST", scf_roundtrip_field_diffs(base, work)) == ("warmup_casePST", String[])
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

    @testset "screening surface (step 4)" begin (function ()
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
      # the screening CSV appends exactly the two screening columns
      csv = read(Sparlectra.writeContingencyResultsCSV(joinpath(d, "flag.csv"), flagged), String)
      header = first(split(csv, '\n'))
      # default format is "technical" (comma delimiter) since issue #376;
      # see "printer and CSV writer" in test_contingency.jl for the excel_de
      # (semicolon) coverage of the same writer.
      @test endswith(header, ",screened,screening_estimate")
      # :only never runs a full solve where an estimate exists; islanding
      # cases still carry real solves (start_used != :screen)
      only_res = Sparlectra.runContingencies!(net, cases14; screening_mode = :only, screening_margin_pct = 10.0)
      @test any(r -> r.screened, only_res)
      @test all(r -> r.screened || r.start_used !== :screen, only_res)
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
      csv1 = readlines(joinpath(dir, "run_fb", "contingency_n1.csv"))
      # default format is "technical" (comma delimiter) since issue #376
      @test endswith(first(csv1), ",screened,screening_estimate")
      @test length(csv1) == 3
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
      @test !endswith(first(csv2), ",screened,screening_estimate")
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
    # weights saved under the generated names (before 0.40.0) still apply
    net14 = Sparlectra.importSCF(sp14)
    old_name = only(Sparlectra.generateN1Branches(net14)[1].elements)
    new_name = only(Sparlectra.generateN1Branches(net14; names = contingency_names(net14, read_scf_json(sp14)))[1].elements)
    @test old_name != new_name
    wfile = joinpath(dir, "sp_case14.contingency-weights.csv")
    write(wfile, "name;weight\n$(old_name);2.5\n")
    dw = run("w", Dict{String,Any}(); weights_path = wfile)
    @test dw["metadata"]["contingency_weighted_cases"] == 1
    wlog = read(joinpath(dir, "w", "run.log"), String)
    @test occursin("1 weight(s) stored under the generated names applied", wlog)
    @test !occursin("match no branch", wlog)
    @test occursin("$(new_name),branch,$(new_name),2.5,", read(joinpath(dir, "w", "contingency_n1.csv"), String))
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

  @testset "scenarios workshop runs with its assertions" begin (function ()
    # the Literate workshop is executable Julia and
    # carries an @assert next to every printed number; RUNNING it here is
    # what keeps the notebook from drifting silently. The workshop runs on
    # the shipped SCF cases sp_case14 and sp_case60 (0.40.0: correlated
    # events), so it needs no download.
    workshop = abspath(joinpath(dirname(@__DIR__), "docs", "lit", "workshop_scenarios.jl"))
    @test isfile(workshop)
    @test isfile(joinpath(dirname(@__DIR__), "data", "scf", "sp_case60.scf.json"))
    # a module with `include`, like a notebook's Main: the first cell
    # includes the workshop's warm-up file (a plain Module() has no include)
    sym = gensym(:WorkshopScenariosRun)
    mod = Core.eval(Main, :(baremodule $sym
      using Base
      eval(x) = Core.eval($sym, x)
      include(x) = Base.include($sym, abspath(x))
    end))
    redirect_stdout(devnull) do
      Base.include(mod, workshop)
    end
    @test true
    println("      scenarios workshop: RAN")
  end)() end

  return nothing
end

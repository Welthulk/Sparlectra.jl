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

# file: test/test_contingency.jl
# purpose: N-1 contingency batch API (multi-core Phase 4): case generation
#          (all branches, transformer filter, FOR001 mapping), serial vs
#          parallel result equality on case14, islanding-without-reference
#          and non-convergence reported (not thrown), base net immutability,
#          table printer and CSV writer.

using Sparlectra
using Test

_contingency_results_equal(a::ContingencyResult, b::ContingencyResult) = all(isequal(getfield(a, f), getfield(b, f)) for f in fieldnames(ContingencyResult))

function run_contingency_tests()
  @testset "N-1 contingency batch (Phase 4)" begin
    case14 = joinpath(Sparlectra.MPOWER_DIR, "case14.m")

    @testset "case generation" begin
      net = createNetFromMatPowerFile(filename = case14)
      cases = generateN1Branches(net)
      @test length(cases) == length(net.branchVec)
      @test all(c.kind === :branch for c in cases)
      lines_only = generateN1Branches(net; include_transformers = false)
      @test length(lines_only) == 17   # case14: 17 lines, 3 transformers
      @test_throws ArgumentError ContingencyCase("x", :bus, "x")

      empty!(net.for001Contingencies)
      append!(net.for001Contingencies, [cases[1].element, "NO_SUCH_BRANCH"])
      for001 = generateContingenciesFromFOR001(net)
      @test length(for001) == 2
      @test for001[2].element == "NO_SUCH_BRANCH"

      # parallel circuits share a component name; every circuit must get its
      # own disambiguated case and remove ITS branch, not always the first
      twin = Net(name = "n1_twin", baseMVA = 100.0)
      for b in ("A", "B")
        addBus!(net = twin, busName = b, vn_kV = 110.0)
      end
      addProsumer!(net = twin, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = twin, busName = "B", type = "ENERGYCONSUMER", p = 20.0, q = 5.0)
      addPIModelACLine!(net = twin, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = twin, fromBus = "A", toBus = "B", r_pu = 0.02, x_pu = 0.16, b_pu = 0.0, status = 1)
      ok_twin, _ = validate!(net = twin)
      @test ok_twin
      twin_names = [getCompName(b.comp) for b in twin.branchVec]
      if twin_names[1] == twin_names[2]
        twin_cases = generateN1Branches(twin)
        @test length(twin_cases) == 2
        @test twin_cases[1].element != twin_cases[2].element
        @test occursin("#", twin_cases[1].element)
        twin_results = runContingencies!(twin, twin_cases; parallel_enabled = false)
        @test all(r -> r.converged, twin_results)
        # removing circuit 1 (low impedance) is a different case than
        # removing circuit 2: the remaining line differs, so do the voltages
        @test twin_results[1].min_vm_pu != twin_results[2].min_vm_pu
        println("contingency parallel-circuit disambiguation: RAN (shared component name)")
      else
        println("contingency parallel-circuit disambiguation: SKIPPED (builder gives parallel circuits unique names: $(twin_names))")
      end

      # screening filters (#331 Phase 2): a 380/110 kV net bridged by a
      # transformer, with rated lines, gives deterministic filter counts
      fnet = Net(name = "n1_filter", baseMVA = 100.0)
      addBus!(net = fnet, busName = "A", vn_kV = 380.0)
      addBus!(net = fnet, busName = "B", vn_kV = 380.0)
      addBus!(net = fnet, busName = "C", vn_kV = 110.0)
      addBus!(net = fnet, busName = "D", vn_kV = 110.0)
      addProsumer!(net = fnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = fnet, busName = "D", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
      addPIModelACLine!(net = fnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1, ratedS = 1200.0)
      addPIModelTrafo!(net = fnet, fromBus = "B", toBus = "C", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1, ratio = 1.0, shift_deg = 0.0)
      addPIModelACLine!(net = fnet, fromBus = "C", toBus = "D", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1, ratedS = 300.0)
      ok_f, _ = validate!(net = fnet)
      @test ok_f
      @test length(generateN1Branches(fnet)) == 3
      @test length(generateN1Branches(fnet; include_transformers = false)) == 2
      # voltage screen keeps the transformer (its high side is 380 kV) and the
      # 380 kV line, drops the 110 kV line
      @test length(generateN1Branches(fnet; min_vn_kV = 200.0)) == 2
      # rating screen: only the 1200 MVA line clears 500; the trafo is unrated
      @test length(generateN1Branches(fnet; min_sn_MVA = 500.0)) == 1
      @test length(generateN1Branches(fnet; min_sn_MVA = 1.0)) == 2
      # name pattern, regex and substring
      @test length(generateN1Branches(fnet; name_pattern = r"_ACL_")) == 2
      @test length(generateN1Branches(fnet; name_pattern = "_2WT_")) == 1
      @test length(generateN1Branches(fnet; name_pattern = "110")) == 1

      # per-case weight: default 1.0, explicit via the 4-arg constructor,
      # validated non-negative and finite
      @test ContingencyCase("x").weight == 1.0
      @test ContingencyCase("x", :branch, "x", 2.5).weight == 2.5
      @test_throws ArgumentError ContingencyCase("x", :branch, "x", -1.0)
    end

    @testset "case14 full N-1: serial and parallel identical" begin
      net = createNetFromMatPowerFile(filename = case14)
      before_branches = length(net.branchVec)
      before_vm = [n._vm_pu for n in net.nodeVec]
      cases = generateN1Branches(net)
      serial = runContingencies!(net, cases; parallel_enabled = false)
      @test length(serial) == length(cases)
      @test count(r -> r.converged, serial) >= length(cases) - 1
      for max_tasks in (1, Threads.nthreads())
        par = runContingencies!(net, cases; parallel_enabled = true, parallel_max_tasks = max_tasks, parallel_min_work_items = 2)
        @test all(_contingency_results_equal(serial[i], par[i]) for i in eachindex(serial))
      end
      # the base net is never mutated: same branch count, same voltages
      @test length(net.branchVec) == before_branches
      @test isequal([n._vm_pu for n in net.nodeVec], before_vm)
      if Threads.nthreads() > 1
        println("contingency parallel batch: RAN with ", Threads.nthreads(), " threads")
      else
        println("contingency parallel batch: fallback-only run (single-threaded test process); the threaded path runs in the --threads=4 battery")
      end
    end

    @testset "islanding without reference is reported, not thrown" begin
      net = Net(name = "n1_island", baseMVA = 100.0)
      for b in ("A", "B", "C", "D")
        addBus!(net = net, busName = b, vn_kV = 110.0)
      end
      addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
      # island {C, D} after cutting B-C carries only loads: no promotable
      # generator, so the matpower_like reference selection must fail
      addProsumer!(net = net, busName = "C", type = "ENERGYCONSUMER", p = 5.0, q = 1.0)
      addProsumer!(net = net, busName = "D", type = "ENERGYCONSUMER", p = 4.0, q = 1.0)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "B", toBus = "C", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "C", toBus = "D", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      ok, _ = validate!(net = net)
      @test ok
      cut = getCompName(getNetBranch(net = net, fromBus = "B", toBus = "C").comp)
      results = runContingencies!(net, [ContingencyCase(cut)]; parallel_enabled = false)
      @test length(results) == 1
      @test !results[1].converged
      # the message is specific: the {C, D} island carries only loads (5 + 4 MW),
      # so it names the load-only cause and the disconnected load in MW
      @test occursin("islanded", results[1].error)
      @test occursin("load-only", results[1].error)
      @test occursin("9.0 MW", results[1].error)
      # the disconnected load is also a structured field, not only in the message
      @test results[1].shed_load_mw == 9.0
      @test results[1].island_count >= 2

      # counterpart: an island WITH a PV generator is promoted to its own
      # reference by the matpower_like policy and the case solves
      net_pv = Net(name = "n1_island_pv", baseMVA = 100.0)
      for b in ("A", "B", "C", "D")
        addBus!(net = net_pv, busName = b, vn_kV = 110.0)
      end
      addProsumer!(net = net_pv, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net_pv, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
      addProsumer!(net = net_pv, busName = "C", type = "GENERATOR", p = 5.0, vm_pu = 1.0, qMin = -10.0, qMax = 10.0)
      addProsumer!(net = net_pv, busName = "D", type = "ENERGYCONSUMER", p = 4.0, q = 1.0)
      addPIModelACLine!(net = net_pv, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net_pv, fromBus = "B", toBus = "C", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net_pv, fromBus = "C", toBus = "D", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      ok_pv, _ = validate!(net = net_pv)
      @test ok_pv
      cut_pv = getCompName(getNetBranch(net = net_pv, fromBus = "B", toBus = "C").comp)
      promoted = runContingencies!(net_pv, [ContingencyCase(cut_pv)]; parallel_enabled = false)
      @test promoted[1].converged
      @test promoted[1].island_count == 2
    end

    @testset "non-convergence and unknown elements are reported" begin
      net = Net(name = "n1_diverge", baseMVA = 100.0)
      for b in ("A", "B", "C")
        addBus!(net = net, busName = b, vn_kV = 110.0)
      end
      addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net, busName = "C", type = "ENERGYCONSUMER", p = 350.0, q = 100.0)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "B", toBus = "C", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "C", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      ok, _ = validate!(net = net)
      @test ok
      direct = getCompName(getNetBranch(net = net, fromBus = "A", toBus = "C").comp)
      results = runContingencies!(net, [ContingencyCase(direct), ContingencyCase("NO_SUCH_BRANCH")]; parallel_enabled = false)
      @test length(results) == 2
      @test !results[1].converged
      @test results[1].error !== nothing
      @test !results[2].converged
      @test occursin("unknown branch", results[2].error)
    end

    @testset "printer and CSV writer" begin
      net = createNetFromMatPowerFile(filename = case14)
      cases = generateN1Branches(net)[1:5]
      results = runContingencies!(net, cases; parallel_enabled = false)
      txt = sprint(io -> printContingencyResults(io, results; max_rows = 3))
      @test occursin("N-1 contingency results", txt)
      @test occursin(results[2].name, txt)
      @test occursin("more row(s) not shown", txt)
      csv = joinpath(mktempdir(), "n1.csv")
      @test writeContingencyResultsCSV(csv, results) == csv
      lines = readlines(csv)
      @test length(lines) == length(results) + 1
      @test startswith(lines[1], "name;weight;converged;iterations;start_used")
    end

    @testset "case weights (#331 Phase 2)" begin
      net = createNetFromMatPowerFile(filename = case14)
      cases = generateN1Branches(net)[1:5]

      # applyContingencyWeights: looked up by name, default for unlisted names
      weights = Dict(cases[1].name => 3.0, cases[2].name => 0.25)
      weighted = applyContingencyWeights(cases, weights; default = 0.5)
      @test length(weighted) == length(cases)
      @test weighted[1].weight == 3.0
      @test weighted[2].weight == 0.25
      @test weighted[3].weight == 0.5           # unlisted name -> default

      # readContingencyWeightsCSV: semicolon, header + comment lines skipped
      dir = mktempdir()
      p1 = joinpath(dir, "w.csv")
      open(p1, "w") do io
        println(io, "name;weight")
        println(io, "# outage rates")
        println(io, cases[1].name, ";4.0")
        println(io, cases[2].name, ";1.5")
      end
      w1 = readContingencyWeightsCSV(p1)
      @test length(w1) == 2
      @test w1[cases[1].name] == 4.0
      # comma delimiter, no header line
      p2 = joinpath(dir, "w2.csv")
      open(p2, "w") do io
        println(io, cases[1].name, ",2.0")
      end
      @test readContingencyWeightsCSV(p2)[cases[1].name] == 2.0
      # a negative weight is rejected
      p3 = joinpath(dir, "w3.csv")
      open(p3, "w") do io
        println(io, cases[1].name, ";-2.0")
      end
      @test_throws ArgumentError readContingencyWeightsCSV(p3)

      # the weight is carried through to the result and shown in the table
      results = runContingencies!(net, weighted; parallel_enabled = false)
      @test all(results[i].weight == weighted[i].weight for i in eachindex(results))
      @test occursin("weight", sprint(io -> printContingencyResults(io, results)))

      # weight 0 is a pure ranking weight, NOT a skip switch: solving a case
      # with weight 0 gives the same outcome as weight 1, only the weight field
      # differs, so nothing is skipped. Dropping an outage is the filters' job.
      r1 = runContingencies!(net, [ContingencyCase(cases[1].name, :branch, cases[1].element, 1.0)]; parallel_enabled = false)[1]
      r0 = runContingencies!(net, [ContingencyCase(cases[1].name, :branch, cases[1].element, 0.0)]; parallel_enabled = false)[1]
      @test r0.weight == 0.0
      @test r1.weight == 1.0
      @test r0.converged == r1.converged    # same solve outcome, not skipped
      @test r0.iterations == r1.iterations
      @test r0.error == r1.error
    end

    @testset "generator outages (#331 Phase 3)" begin
      net = createNetFromMatPowerFile(filename = case14)
      gcases = generateN1Generators(net)
      # case14 has 5 generators; every case is kind = :gen
      @test length(gcases) == 5
      @test all(c -> c.kind === :gen, gcases)
      # generator names are not unique, so they are disambiguated by index
      @test all(c -> occursin("#", c.element), gcases)

      # Pg and name filters
      @test length(generateN1Generators(net; min_pg_MW = 1e9)) == 0
      @test length(generateN1Generators(net; name_pattern = r"NoSuchGen")) == 0

      # a generator outage removes only that unit; the slack picks up the loss.
      # Removing case14's single slack leaves the net reference-less (reported,
      # not thrown); auto_slack promotes the strongest survivor so all solve.
      res = runContingencies!(net, gcases; parallel_enabled = false)
      @test count(r -> r.converged, res) == 4
      @test any(r -> !r.converged && occursin("no slack bus", r.error), res)
      res_auto = runContingencies!(net, gcases; parallel_enabled = false, auto_slack = true)
      @test count(r -> r.converged, res_auto) == 5

      # distributed slack flows through as a keyword; the batch still solves
      res_ds = runContingencies!(net, gcases; parallel_enabled = false, auto_slack = true, distributed_slack_enabled = true)
      @test count(r -> r.converged, res_ds) == 5

      # serial vs parallel identity holds for generator cases too
      gs = runContingencies!(net, gcases; parallel_enabled = false)
      gp = runContingencies!(net, gcases; parallel_enabled = true, parallel_max_tasks = 2, parallel_min_work_items = 2)
      @test all(_contingency_results_equal(gs[i], gp[i]) for i in eachindex(gs))

      # unknown generator element is reported, not thrown
      unk = runContingencies!(net, [ContingencyCase("NO_SUCH_GEN", :gen, "NO_SUCH_GEN")]; parallel_enabled = false)
      @test !unk[1].converged
      @test occursin("unknown generator", unk[1].error)

      # stranded-generation fixture (reviewer request, non-pegase): two separate
      # areas, area 2 = bus C (a regulating PV generator AND a fixed-injection PQ
      # generator) feeding load at D. Both areas carry a reference in the base
      # case. A generator outage that removes C's PV generator leaves area 2 with
      # injection (the PQ unit) but no voltage reference -> the parked
      # stranded-generation outcome, reproduced without pegase.
      sg = Net(name = "n1_stranded_gen", baseMVA = 100.0)
      for b in ("A", "A2", "C", "D")
        addBus!(net = sg, busName = b, vn_kV = 110.0)
      end
      addProsumer!(net = sg, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = sg, busName = "A2", type = "ENERGYCONSUMER", p = 5.0, q = 1.0)
      addPIModelACLine!(net = sg, fromBus = "A", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addProsumer!(net = sg, busName = "C", type = "GENERATOR", p = 30.0, vm_pu = 1.02)  # PV (regulating, the area reference)
      addProsumer!(net = sg, busName = "C", type = "GENERATOR", p = 12.0, q = 3.0)         # PQ (fixed injection)
      addProsumer!(net = sg, busName = "D", type = "ENERGYCONSUMER", p = 20.0, q = 6.0)
      addPIModelACLine!(net = sg, fromBus = "C", toBus = "D", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      ok_sg, _ = validate!(net = sg)
      @test ok_sg
      sres = runContingencies!(sg, generateN1Generators(sg); parallel_enabled = false)
      smsgs = [r.error === nothing ? "converged" : r.error for r in sres]
      # PV generator removed, PQ generator survives: injection but no reference
      @test any(m -> occursin("generation stranded", m), smsgs)
      # slack generator removed: area 1 becomes a load-only island
      @test any(m -> occursin("load-only", m), smsgs)
      # PQ generator removed, the PV reference survives: area 2 still solves
      @test any(m -> m == "converged", smsgs)
    end

    @testset "overload reporting (#331 Phase 4)" begin
      # two parallel rated lines feed a load; removing one overloads the other
      net = Net(name = "n1_overload", baseMVA = 100.0)
      addBus!(net = net, busName = "A", vn_kV = 220.0)
      addBus!(net = net, busName = "B", vn_kV = 220.0)
      addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 100.0, q = 20.0)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.005, x_pu = 0.04, b_pu = 0.0, status = 1, ratedS = 60.0)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.005, x_pu = 0.04, b_pu = 0.0, status = 1, ratedS = 60.0)
      ok_o, _ = validate!(net = net)
      @test ok_o
      cases = generateN1Branches(net)
      @test length(cases) == 2
      res = runContingencies!(net, cases; parallel_enabled = false)
      # the surviving line is overloaded, reported as a structured OverloadRecord
      over = [r for r in res if !isempty(r.overloads)]
      @test !isempty(over)
      @test isa(over[1].overloads[1], OverloadRecord)
      @test over[1].overloads[1].loading_pct > 100.0
      @test over[1].max_branch_loading_pct == over[1].overloads[1].loading_pct
      # the record carries the base loading and the delta to it
      o1 = over[1].overloads[1]
      @test isfinite(o1.loading_base_pct) && o1.loading_base_pct < 100.0   # base within limits
      @test o1.delta_pct ≈ o1.loading_pct - o1.loading_base_pct
      @test o1.sn_MVA == 60.0
      @test o1.s_MVA ≈ o1.loading_pct / 100.0 * o1.sn_MVA
      # severity = weight * max(0, max_loading - 100) for a converged overload
      @test over[1].severity ≈ over[1].weight * (over[1].max_branch_loading_pct - 100.0)
      @test over[1].severity > 0.0

      # aggregate report over the batch
      rep = buildContingencyReport(res)
      @test rep.n_cases == 2
      @test rep.n_with_overload >= 1
      @test rep.worst_loading_pct > 100.0
      @test rep.worst_severity > 0.0
      @test rep.worst_severity_case != ""
      @test !isempty(rep.top_overloaded)
      txt = sprint(printContingencyReport, rep)
      @test occursin("N-1 contingency report", txt)
      @test occursin("most-overloaded branches", txt)

      # shed_load_mw is a structured field and feeds the report's shed total.
      # Cutting M-L1 strands the load-only island {L1, L2} (7 + 3 MW); the slack
      # side {S, M} keeps its reference (S stays connected via S-M), so this is
      # the island-validation path, not a whole-net "no slack" failure.
      isl = Net(name = "n1_shed", baseMVA = 100.0)
      for b in ("S", "M", "L1", "L2")
        addBus!(net = isl, busName = b, vn_kV = 110.0)
      end
      addProsumer!(net = isl, busName = "S", type = "EXTERNALNETWORKINJECTION", referencePri = "S", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = isl, busName = "M", type = "ENERGYCONSUMER", p = 2.0, q = 1.0)
      addProsumer!(net = isl, busName = "L1", type = "ENERGYCONSUMER", p = 7.0, q = 2.0)
      addProsumer!(net = isl, busName = "L2", type = "ENERGYCONSUMER", p = 3.0, q = 1.0)
      addPIModelACLine!(net = isl, fromBus = "S", toBus = "M", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = isl, fromBus = "M", toBus = "L1", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = isl, fromBus = "L1", toBus = "L2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      validate!(net = isl)
      cut = getCompName(getNetBranch(net = isl, fromBus = "M", toBus = "L1").comp)
      sres = runContingencies!(isl, [ContingencyCase(cut)]; parallel_enabled = false)
      @test sres[1].shed_load_mw == 10.0
      srep = buildContingencyReport(sres)
      @test srep.total_shed_load_mw == 10.0
      @test srep.n_islanded == 1
      # a failed/islanded case has NaN severity so it sorts to the TOP
      @test isnan(sres[1].severity)
      combined = vcat(sres, res)
      out = sprint(io -> printContingencyResults(io, combined))          # default :severity sort
      @test findfirst(sres[1].name, out)[1] < findfirst(over[1].name, out)[1]
      # :none keeps input order (islanded case first, as given)
      out_input = sprint(io -> printContingencyResults(io, combined; sort_by = :none))
      @test findfirst(sres[1].name, out_input)[1] < findfirst(res[end].name, out_input)[1]
    end

    @testset "Web UI contingency service (#331 Phase 5)" begin
      cfg = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH
      root = mktempdir()
      dicts = Dict{String,Any}()
      for kind in ("branch", "gen")
        od = joinpath(root, "ct_$(kind)")
        res = redirect_stdout(devnull) do
          Sparlectra._run_contingency_service(case14, cfg, od, "ct_$(kind)", kind)
        end
        d = Sparlectra.to_dict(res)
        dicts[kind] = d
        @test d["status"] == "succeeded"
        @test isfile(joinpath(od, "contingency_n1.csv"))       # artifacts written to the run dir
        @test isfile(joinpath(od, "run.log"))
        md = d["metadata"]
        @test md["run_mode"] == "contingency"
        @test md["contingency_kind"] == kind
        @test md["contingency_cases"] > 0
        # explicit run-status keys so the result view renders like a completed run
        @test md["run_status"] == "completed"
        @test Sparlectra._webui_contingency_summary(d) !== nothing
      end
      # a generator outage on case14 removes the only slack: reported (not thrown),
      # counted as no_slack, and named in the summary badge so it does not read
      # as a tool failure
      @test dicts["gen"]["metadata"]["contingency_no_slack"] >= 1
      @test occursin("slack", Sparlectra._webui_contingency_summary(dicts["gen"]))
      # a non-contingency run gets no contingency summary row
      @test Sparlectra._webui_contingency_summary(Dict("metadata" => Dict("run_mode" => "powerflow"))) === nothing
      # invalid kind is rejected, not thrown
      bad = redirect_stdout(devnull) do
        Sparlectra._run_contingency_service(case14, cfg, joinpath(root, "ct_bad"), "ct_bad", "nonsense")
      end
      @test Sparlectra.to_dict(bad)["status"] == "failed"
      @test Sparlectra.to_dict(bad)["reason"] == "invalid_request"

      # per-case weights (#331 Phase 5 follow-up): applied when a weight file is
      # present next to the case, warned (not fatal) on unmatched names, and
      # absent by default. A copy in a temp dir keeps the packaged case clean.
      wdir = mktempdir()
      w14 = joinpath(wdir, "case14.m")
      cp(case14, w14)
      wnames = redirect_stdout(devnull) do
        [c.name for c in generateN1Branches(Sparlectra._import_sparlectra_net(w14, nothing, Sparlectra.load_sparlectra_config(cfg; reload = true)))]
      end
      wf = Sparlectra._webui_case_weights_path(w14)
      open(wf, "w") do io
        println(io, "name;weight")
        println(io, wnames[1], ";4.0")
        println(io, "NO_SUCH_ELEMENT;2.0")
      end
      wmd = redirect_stdout(devnull) do
        Sparlectra.to_dict(Sparlectra._run_contingency_service(w14, cfg, joinpath(root, "wt"), "wt", "branch"; weights_path = wf))["metadata"]
      end
      @test wmd["contingency_weights_applied"] == true
      @test wmd["contingency_weighted_cases"] == 1
      @test any(occursin("match no", l) for l in readlines(joinpath(root, "wt", "run.log")))
      nmd = redirect_stdout(devnull) do
        Sparlectra.to_dict(Sparlectra._run_contingency_service(w14, cfg, joinpath(root, "nw"), "nw", "branch"))["metadata"]
      end
      @test nmd["contingency_weights_applied"] == false
    end

    @testset "start-value ladder (#331 Phase 1)" begin
      net = createNetFromMatPowerFile(filename = case14)
      cases = generateN1Branches(net)

      # default [:warm]: converged cases report :warm, failed cases :none
      base = runContingencies!(net, cases)
      @test all(r -> r.converged ? r.start_used === :warm : r.start_used === :none, base)
      # the table gains a start column, the CSV a start_used field
      @test occursin("start", sprint(io -> printContingencyResults(io, base[1:3])))

      # ladder validation: non-empty, allowed stages only, no duplicates
      @test_throws ArgumentError runContingencies!(net, cases[1:1]; rescue_ladder = Symbol[])
      @test_throws ArgumentError runContingencies!(net, cases[1:1]; rescue_ladder = [:warm, :warm])
      @test_throws ArgumentError runContingencies!(net, cases[1:1]; rescue_ladder = [:warm, :bogus])
      # the config section validates the same way
      @test Sparlectra.SparlectraConfig(Dict("contingency" => Dict("rescue_ladder" => ["warm", "dc"]))).contingency.rescue_ladder == [:warm, :dc]
      @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("contingency" => Dict("rescue_ladder" => ["warm", "warm"])))

      # a longer ladder never loses a case that already converged warm, and its
      # extra stages only ever help (bitwise-equal for the warm-converging cases)
      longer = runContingencies!(net, cases; rescue_ladder = [:warm, :dc, :flat])
      @test count(r -> r.converged, longer) >= count(r -> r.converged, base)
      @test all(base[i].converged ? _contingency_results_equal(base[i], longer[i]) : true for i in eachindex(base))

      # retry_flat_start is a deprecated alias for appending :flat; converged
      # cases are unchanged
      retried = runContingencies!(net, cases; retry_flat_start = true)
      @test all(base[i].converged ? _contingency_results_equal(base[i], retried[i]) : true for i in eachindex(base))

      # alias collision: retry_flat_start = true with a ladder that ALREADY
      # contains :flat must not trip the duplicate validation; the alias is a
      # no-op there (no throw, identical results to the plain ladder)
      with_flat = runContingencies!(net, cases; rescue_ladder = [:warm, :flat])
      with_flat_alias = runContingencies!(net, cases; rescue_ladder = [:warm, :flat], retry_flat_start = true)
      @test all(_contingency_results_equal(with_flat[i], with_flat_alias[i]) for i in eachindex(with_flat))

      # serial vs parallel identity holds with a multi-stage ladder too
      serial = runContingencies!(net, cases; rescue_ladder = [:warm, :dc], parallel_enabled = false)
      par = runContingencies!(net, cases; rescue_ladder = [:warm, :dc], parallel_enabled = true, parallel_max_tasks = 2, parallel_min_work_items = 2)
      @test all(_contingency_results_equal(serial[i], par[i]) for i in eachindex(serial))
    end

    @testset "non-converged base case runs the solver rescue ladder (#331)" begin
      # a mild net seeded with a deliberately bad start: the plain solve from
      # that seed diverges, but the solver rescue ladder (alternate start)
      # recovers it. runContingencies! must rescue the base rather than fall
      # straight to the flat template.
      function _bad_seed_net()
        n = Net(name = "n1_base_rescue", baseMVA = 100.0)
        for b in ("S", "M", "L")
          addBus!(net = n, busName = b, vn_kV = 110.0)
        end
        addProsumer!(net = n, busName = "S", type = "EXTERNALNETWORKINJECTION", referencePri = "S", vm_pu = 1.0, va_deg = 0.0)
        addProsumer!(net = n, busName = "L", type = "ENERGYCONSUMER", p = 40.0, q = 15.0)
        addPIModelACLine!(net = n, fromBus = "S", toBus = "M", r_pu = 0.02, x_pu = 0.10, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = n, fromBus = "M", toBus = "L", r_pu = 0.02, x_pu = 0.10, b_pu = 0.0, status = 1)
        validate!(net = n)
        n.flatstart = false
        for node in n.nodeVec
          getNodeType(node) === Sparlectra.Slack && continue
          node._vm_pu = 0.02
          node._va_deg = -179.0
        end
        return n
      end
      # pre-conditions: plain fails, the rescue ladder converges
      _, erg_plain = runpf!(_bad_seed_net(), 30, 1e-8, 0; islands_enabled = true)
      @test erg_plain != 0
      _, erg_rescue = runpf!(_bad_seed_net(), PowerFlowConfig(rescue = true, islands_enabled = true, max_iter = 30, tol = 1e-8))
      @test erg_rescue == 0

      # the batch must take the rescue branch (no flat-fallback warning) and
      # produce converged, warm-started contingencies
      base_net = _bad_seed_net()
      ccases = generateN1Branches(base_net)
      logs, results = Test.collect_test_logs() do
        runContingencies!(base_net, ccases; parallel_enabled = false)
      end
      # the base was rescued, so the FLAT-fallback warning must NOT fire (the
      # direct proof that the rescue branch was taken before the flat fallback)
      @test !any(occursin("did not converge even with the solver rescue ladder", string(l.message)) for l in logs)
      # the batch ran the outages on the rescued (warm) base: at least one
      # outage keeps the net solvable and converges warm from it (had the base
      # fallen to the flat template, there would be no rescued warm state)
      @test length(results) == length(ccases)
      @test any(r -> r.converged && r.start_used === :warm, results)
    end
  end
end

# Extended profile: a small N-1 identity slice on the real case1354pegase
# (needs the Web UI case cache; states RAN or SKIPPED). Serial vs parallel
# field equality on a branch subset with max_tasks = 1 and auto, plus the
# retry_flat_start smoke (must not change converged cases).
function run_contingency_extended_tests()
  @testset "N-1 identity on case1354pegase (extended)" begin
    case_path = joinpath(homedir(), ".local", "state", "sparlectra", "webui", "data", "mpower", "case1354pegase.m")
    if !isfile(case_path)
      println("contingency extended case1354: SKIPPED (case not cached under the Web UI data directory)")
      @test true
      return
    end
    # canonical pegase import convention (from the case sidecar): angles in rad,
    # shift sign -1.0, tap ratio normal (the default; NOT reciprocal). Made
    # explicit so the fixture documents the convention rather than relying on it.
    net = createNetFromMatPowerFile(filename = case_path, matpower_shift_unit = :rad, matpower_shift_sign = -1.0, matpower_ratio = :normal)
    cases = generateN1Branches(net)[1:60]
    serial = runContingencies!(net, cases; parallel_enabled = false)
    @test length(serial) == 60
    for max_tasks in (1, Threads.nthreads())
      par = runContingencies!(net, cases; parallel_enabled = true, parallel_max_tasks = max_tasks, parallel_min_work_items = 2)
      @test all(_contingency_results_equal(serial[i], par[i]) for i in eachindex(serial))
    end
    retried = runContingencies!(net, cases; parallel_enabled = true, parallel_min_work_items = 2, retry_flat_start = true)
    @test all(!serial[i].converged || _contingency_results_equal(serial[i], retried[i]) for i in eachindex(serial))
    println("contingency extended case1354: RAN (", count(r -> r.converged, serial), "/60 converged, threads = ", Threads.nthreads(), ")")
  end
end

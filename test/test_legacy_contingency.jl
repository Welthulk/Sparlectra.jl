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

# file: test/test_legacy_contingency.jl
# purpose: legacy profile: the testsets of test/test_contingency.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_contingency.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_contingency_tests()
  # cut from testset: N-1 contingency batch (Phase 4)
  (function ()
      # load_fixture_net: the service legs run on the shipped SCF case
      service_case = abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"))

      @testset "case generation" begin (function ()
        net = load_fixture_net("sp_case14")
        cases = generateN1Branches(net)
        @test length(cases) == length(net.branchVec)
        @test all(c.kind === :branch for c in cases)
        lines_only = generateN1Branches(net; include_transformers = false)
        @test length(lines_only) == 15   # sp_case14: 15 lines, 2 transformers
        @test_throws ArgumentError ContingencyCase("x", :bus, "x")

        empty!(net.for001Contingencies)
        append!(net.for001Contingencies, [only(cases[1].elements), "NO_SUCH_BRANCH"])
        for001 = generateContingenciesFromFOR001(net)
        @test length(for001) == 2
        @test only(for001[2].elements) == "NO_SUCH_BRANCH"

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
          @test only(twin_cases[1].elements) != only(twin_cases[2].elements)
          @test occursin("#", only(twin_cases[1].elements))
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
      end)() end

      @testset "printer and CSV writer" begin (function ()
        net = load_fixture_net("sp_case14")
        cases = generateN1Branches(net)[1:5]
        results = runContingencies!(net, cases; parallel_enabled = false)
        txt = sprint(io -> printContingencyResults(io, results; max_rows = 3))
        @test occursin("Contingency results", txt)
        @test occursin(results[2].name, txt)
        @test occursin("more row(s) not shown", txt)
        csv = joinpath(mktempdir(), "n1.csv")
        # the default under the portable locale (`auto` follows the locale)
        @test withenv(() -> writeContingencyResultsCSV(csv, results), "LC_ALL" => "C.UTF-8") == csv
        lines = readlines(csv)
        @test length(lines) == length(results) + 1
        # Since issue #376 the writer follows the same run-wide CSV format as
        # every other artifact; the default is now "technical" (comma
        # delimiter), matching output.csv_format's default. "excel_de"
        # reproduces the historical hardcoded semicolon delimiter.
        @test startswith(lines[1], "name,kind,elements,weight,converged,iterations,start_used")
        csv_de = joinpath(mktempdir(), "n1_de.csv")
        @test writeContingencyResultsCSV(csv_de, results; format = "excel_de") == csv_de
        lines_de = readlines(csv_de)
        @test startswith(lines_de[1], "name;kind;elements;weight;converged;iterations;start_used")
      end)() end

      @testset "case weights (#331 Phase 2)" begin (function ()
        net = load_fixture_net("sp_case14")
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
        r1 = runContingencies!(net, [ContingencyCase(cases[1].name, :branch, only(cases[1].elements), 1.0)]; parallel_enabled = false)[1]
        r0 = runContingencies!(net, [ContingencyCase(cases[1].name, :branch, only(cases[1].elements), 0.0)]; parallel_enabled = false)[1]
        @test r0.weight == 0.0
        @test r1.weight == 1.0
        @test r0.converged == r1.converged    # same solve outcome, not skipped
        @test r0.iterations == r1.iterations
        @test r0.error == r1.error
      end)() end

      @testset "generator outages (#331 Phase 3)" begin (function ()
        net = load_fixture_net("sp_case14")
        gcases = generateN1Generators(net)
        # sp_case14 has 3 generator-type units (external grid plus two
        # machines); every case is kind = :gen
        @test length(gcases) == 3
        @test all(c -> c.kind === :gen, gcases)
        # generator names are not unique, so they are disambiguated by index
        @test all(c -> occursin("#", only(c.elements)), gcases)

        # Pg and name filters
        @test length(generateN1Generators(net; min_pg_MW = 1e9)) == 0
        @test length(generateN1Generators(net; name_pattern = r"NoSuchGen")) == 0

        # a generator outage removes only that unit; the slack picks up the loss.
        # Removing sp_case14's single slack leaves the net reference-less. With
        # auto_slack = false that ends the case (reported, not thrown); by
        # default the strongest survivor takes over and all three solve.
        res = runContingencies!(net, gcases; parallel_enabled = false, auto_slack = false)
        @test count(r -> r.converged, res) == 2
        @test any(r -> !r.converged && occursin("no slack bus", r.error), res)
        res_auto = redirect_stdout(devnull) do
          runContingencies!(net, gcases; parallel_enabled = false)
        end
        @test count(r -> r.converged, res_auto) == 3
        @test any(r -> r.error !== nothing && occursin("reference taken over by bus", r.error), res_auto)

        # distributed slack flows through as a keyword; the batch still solves
        res_ds = redirect_stdout(devnull) do
          runContingencies!(net, gcases; parallel_enabled = false, distributed_slack_enabled = true)
        end
        @test count(r -> r.converged, res_ds) == 3

        # serial vs parallel identity for generator cases: run_contingency_parallel_tests

        # unknown generator element is reported, not thrown
        unk = runContingencies!(net, [ContingencyCase("NO_SUCH_GEN", :gen, "NO_SUCH_GEN")]; parallel_enabled = false)
        @test !unk[1].converged
        @test occursin("unknown generator", unk[1].error)

        # reference hand-over fixture (non-pegase): two separate areas, area
        # 2 = bus C (a regulating PV generator AND a fixed-injection PQ
        # generator) feeding load at D. Both areas carry a reference in the base
        # case. A generator outage that removes C's PV generator leaves area 2
        # with injection (the PQ unit) but no voltage reference: the island
        # takes that unit as its reference, with or without auto_slack (islands
        # always find themselves a reference); only the load-only area stays
        # islanded.
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
        sres = redirect_stdout(devnull) do
          runContingencies!(sg, generateN1Generators(sg); parallel_enabled = false, auto_slack = false)
        end
        smsgs = [r.error === nothing ? "converged" : r.error for r in sres]
        # slack generator removed: area 1 becomes a load-only island
        @test any(m -> occursin("load-only", m), smsgs)
        # PV generator removed, the PQ unit survives and takes the reference;
        # PQ generator removed, the PV reference survives: both solve
        @test count(r -> r.converged, sres) == 2
        @test !any(m -> occursin("generation stranded", m), smsgs)
        # the default (auto_slack) ends on the same three outcomes
        sauto = redirect_stdout(devnull) do
          runContingencies!(sg, generateN1Generators(sg); parallel_enabled = false)
        end
        @test count(r -> r.converged, sauto) == 2
        @test only(r for r in sauto if !r.converged).error == "islanded: load-only, 5.0 MW load disconnected"

        # 0.30.2 reference priority on the shipped two-island case: the outage
        # of the north slack hands its island to the unit with the stated
        # priority (Brunntal_110, 2), not to the largest unit (Kesselau_110),
        # which takes over once the priorities are cleared
        pnet = importSCF(joinpath(dirname(@__DIR__), "data", "scf", "two_islands_prio.scf.json"))
        pslack = findfirst(Sparlectra.isSlack, pnet.prosumpsVec)
        pname = getCompName(pnet.prosumpsVec[pslack].comp)
        pcase = only(c for c in generateN1Generators(pnet) if only(c.elements) in (pname, string(pname, "#", pslack)))
        pplain = deepcopy(pnet)
        foreach(ps -> ps.referencePriority = 0, pplain.prosumpsVec)
        for (n, bus) in ((pnet, "Brunntal_110"), (pplain, "Kesselau_110"))
          r = only(redirect_stdout(devnull) do
            runContingencies!(n, [pcase]; parallel_enabled = false)
          end)
          @test r.converged
          @test occursin("reference taken over by bus $(bus) ", r.error)
        end
      end)() end

      @testset "overload reporting (#331 Phase 4)" begin (function ()
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
        @test occursin("Contingency report", txt)
        @test occursin("most-overloaded branches", txt)
        # parallel circuits share a name: three overload records of one name in
        # ONE contingency count that contingency once (case6495rte reported a
        # branch overloaded by 21045 of 9019 contingencies)
        rec = OverloadRecord("PAR", 112.0, 90.0, 22.0, 112.0, 100.0)
        twins = [ContingencyResult("out$(k)", :branch, ["out$(k)"], 1.0, true, 3, :warm, 1.0, 0.95, 112.0, 12.0, [rec, rec, rec], String[], 1, 0.0, nothing) for k in 1:2]
        @test buildContingencyReport(twins).top_overloaded == ["PAR" => 2]

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
      end)() end

      @testset "Web UI contingency service (#331 Phase 5)" begin (function ()
        cfg = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH
        root = mktempdir()
        dicts = Dict{String,Any}()
        for kind in ("branch", "gen")
          od = joinpath(root, "ct_$(kind)")
          res = redirect_stdout(devnull) do
            SparlectraApp._run_contingency_service(service_case, cfg, od, "ct_$(kind)", kind)
          end
          d = SparlectraApp.to_dict(res)
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
          @test SparlectraApp._webui_contingency_summary(d) !== nothing
        end
        # a generator outage on sp_case14 removes the only slack: the strongest
        # remaining unit takes over, the case solves, the row names the bus and
        # the summary counts it
        @test dicts["gen"]["metadata"]["contingency_no_slack"] == 0
        @test dicts["gen"]["metadata"]["contingency_reference_taken_over"] >= 1
        @test dicts["gen"]["metadata"]["contingency_converged"] == dicts["gen"]["metadata"]["contingency_cases"]
        @test occursin("on a reference another unit took over", SparlectraApp._webui_contingency_summary(dicts["gen"]))
        @test occursin("reference taken over by bus", read(joinpath(dicts["gen"]["output_dir"], "contingency_n1.csv"), String))
        # a non-contingency run gets no contingency summary row
        @test SparlectraApp._webui_contingency_summary(Dict("metadata" => Dict("run_mode" => "powerflow"))) === nothing
        # invalid kind is rejected, not thrown
        bad = redirect_stdout(devnull) do
          SparlectraApp._run_contingency_service(service_case, cfg, joinpath(root, "ct_bad"), "ct_bad", "nonsense")
        end
        @test SparlectraApp.to_dict(bad)["status"] == "failed"
        @test SparlectraApp.to_dict(bad)["reason"] == "invalid_request"

        # per-case weights (#331 Phase 5 follow-up): applied when a weight file is
        # present next to the case, warned (not fatal) on unmatched names, and
        # absent by default. A copy in a temp dir keeps the packaged case clean.
        wdir = mktempdir()
        w14 = joinpath(wdir, "sp_case14.scf.json")
        cp(service_case, w14)
        wnames = redirect_stdout(devnull) do
          [c.name for c in generateN1Branches(Sparlectra._import_sparlectra_net(w14, nothing, Sparlectra.load_sparlectra_config(cfg; reload = true)))]
        end
        wf = SparlectraApp._webui_case_weights_path(w14)
        open(wf, "w") do io
          println(io, "name;weight")
          println(io, wnames[1], ";4.0")
          println(io, "NO_SUCH_ELEMENT;2.0")
        end
        wmd = redirect_stdout(devnull) do
          SparlectraApp.to_dict(SparlectraApp._run_contingency_service(w14, cfg, joinpath(root, "wt"), "wt", "branch"; weights_path = wf))["metadata"]
        end
        @test wmd["contingency_weights_applied"] == true
        @test wmd["contingency_weighted_cases"] == 1
        @test any(occursin("match no", l) for l in readlines(joinpath(root, "wt", "run.log")))
        nmd = redirect_stdout(devnull) do
          SparlectraApp.to_dict(SparlectraApp._run_contingency_service(w14, cfg, joinpath(root, "nw"), "nw", "branch"))["metadata"]
        end
        @test nmd["contingency_weights_applied"] == false
      end)() end

      @testset "start-value ladder (#331 Phase 1)" begin (function ()
        net = load_fixture_net("sp_case14")
        cases = generateN1Branches(net)

        # default [:warm]: converged cases report :warm, failed cases :none
        base = runContingencies!(net, cases)
        @test all(r -> r.converged ? r.start_used === :warm : r.start_used === :none, base)
        # :warm starts from the solved base voltages even when the imported
        # net carries a flat-start flag (power_flow.flatstart = true): the
        # flag only changes how the BASE case starts, so every outage takes
        # the same iterations and lands on the same solution. Before 0.31.1
        # the flag turned each :warm stage into a plain flat start
        # (case6495rte: no outage converged).
        flagged = load_fixture_net("sp_case14")
        flagged.flatstart = true
        warm_flag = runContingencies!(flagged, cases)
        @test [r.iterations for r in warm_flag] == [r.iterations for r in base]
        @test all(isapprox(a.min_vm_pu, b.min_vm_pu; atol = 1e-9) || (isnan(a.min_vm_pu) && isnan(b.min_vm_pu)) for (a, b) in zip(warm_flag, base))
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
        # serial vs parallel identity with a multi-stage ladder:
        # run_contingency_parallel_tests
      end)() end

      @testset "non-converged base case runs the solver rescue ladder (#331)" begin (function ()
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

        # #456: the rescue and the :apslf stage solve on the batch's own
        # keywords. sp_case14 with distributed slack, seeded so badly that the
        # base solve fails: the rescued template must carry the distributed
        # solution (the flat-start reference), not the single-slack one the
        # rescue used to fall back to (template angles were 4.8 degrees off the
        # distributed reference and equal to the single-slack one).
        function _bad_seed_14()
          n = load_fixture_net("sp_case14")
          n.flatstart = false
          for node in n.nodeVec
            getNodeType(node) === Sparlectra.Slack && continue
            node._vm_pu = 0.02
            node._va_deg = -179.0
          end
          return n
        end
        function _flat_ref(; kw...)
          n = load_fixture_net("sp_case14")
          n.flatstart = true
          _, erg = runpf!(n, 30, 1e-8, 0; islands_enabled = true, kw...)
          return (n, erg)
        end
        # the participation warnings of pg-weighted sharing on this case are
        # expected and not the subject here
        _, (ref_ds, ref_ss, erg_seed, engine) = Test.collect_test_logs() do
          redirect_stdout(devnull) do
            ds = _flat_ref(distributed_slack_enabled = true)
            ss = _flat_ref()
            _, e = runpf!(_bad_seed_14(), 30, 1e-8, 0; islands_enabled = true, distributed_slack_enabled = true)
            eng = Sparlectra.ScenarioEngine(_bad_seed_14(); pf_kwargs = (; distributed_slack_enabled = true))
            (ds, ss, e, eng)
          end
        end
        # pre-conditions: the seed defeats the plain solve, and the two slack
        # models give clearly different angles on this case
        @test ref_ds[2] == 0 && ref_ss[2] == 0
        ref_ds, ref_ss = ref_ds[1], ref_ss[1]
        @test erg_seed != 0
        dva(a, b) = maximum(abs(a.nodeVec[i]._va_deg - b.nodeVec[i]._va_deg) for i in eachindex(a.nodeVec))
        @test dva(ref_ds, ref_ss) > 1.0
        @test engine.base_converged
        @test dva(engine.template, ref_ds) < 1e-4
        # the :apslf stage reaches the same distributed-slack solution as the
        # :warm stage on every outage both solve (it used to differ by up to
        # 30 loading points: a single-slack solve)
        net14 = load_fixture_net("sp_case14")
        cases14 = generateN1Branches(net14)
        _, (warm14, apslf14) = Test.collect_test_logs() do
          redirect_stdout(devnull) do
            w = runContingencies!(net14, cases14; parallel_enabled = false, distributed_slack_enabled = true, rescue_ladder = [:warm])
            a = runContingencies!(net14, cases14; parallel_enabled = false, distributed_slack_enabled = true, rescue_ladder = [:apslf])
            (w, a)
          end
        end
        both = [i for i in eachindex(warm14) if warm14[i].converged && apslf14[i].converged]
        @test !isempty(both)
        @test all(i -> apslf14[i].start_used === :apslf, both)
        @test maximum(abs(warm14[i].max_branch_loading_pct - apslf14[i].max_branch_loading_pct) for i in both) < 1e-3
      end)() end

  end)()
  return true
end

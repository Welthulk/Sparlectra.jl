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
#          parallel result equality on the shipped sp_case14 (load_fixture_net), islanding-without-reference
#          and non-convergence reported (not thrown), base net immutability,
#          table printer and CSV writer.

using Sparlectra
using Test

_contingency_results_equal(a::ContingencyResult, b::ContingencyResult) = all(isequal(getfield(a, f), getfield(b, f)) for f in fieldnames(ContingencyResult))

function run_contingency_tests()
  @testset "N-1 contingency batch (Phase 4)" begin (function ()

    run_contingency_parallel_tests()

    @testset "islanding without reference is reported, not thrown" begin (function ()
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
    end)() end

    @testset "non-convergence and unknown elements are reported" begin (function ()
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
    end)() end

    @testset "an outage that cuts off a bus is reported as islanding" begin (function ()
      # A (reference) feeds B, and B feeds the radial buses C (a load) and G
      # (a unit). The outage of B-C leaves C without a connection: the
      # remaining network solves, and the result names the bus and the load
      # it lost instead of reading as a clean case. The outage of B-G takes
      # the unit out of service the same way.
      net = Net(name = "n1_cut_off", baseMVA = 100.0)
      for b in ("A", "B", "C", "G")
        addBus!(net = net, busName = b, vn_kV = 110.0)
      end
      addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 20.0, q = 5.0)
      addProsumer!(net = net, busName = "C", type = "ENERGYCONSUMER", p = 7.0, q = 2.0)
      addProsumer!(net = net, busName = "G", type = "GENERATOR", p = 4.0, q = 1.0)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "B", toBus = "C", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "B", toBus = "G", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      ok, _ = validate!(net = net)
      @test ok
      res = runContingencies!(net, generateN1Branches(net); parallel_enabled = false)
      @test length(res) == 4
      @test all(r -> r.converged, res)
      load_cut = only(r for r in res if r.shed_load_mw > 0.0)
      @test load_cut.shed_load_mw == 7.0
      @test startswith(load_cut.error, "islanded: bus C cut off, 7.0 MW load disconnected")
      # the same cut-off with a Q(U) controller anywhere in the net: the
      # controller rebuilds the injections every iteration, and the cut-off
      # load used to come back as a residual no step could remove (case6495rte
      # with its P(U)/Q(U) controllers: the outage stagnated at that load)
      qnet = deepcopy(net)
      addProsumer!(net = qnet, busName = "B", type = "GENERATOR", p = 2.0, q = 0.0, qu_controller = QUController(VoltageCharacteristic([(0.95, 0.05), (1.05, -0.05)]); qmin_pu = -0.1, qmax_pu = 0.1))
      qres = runContingencies!(qnet, [c for c in generateN1Branches(qnet) if only(c.elements) == only(r.name for r in res if r === load_cut)]; parallel_enabled = false)
      @test only(qres).converged
      @test startswith(only(qres).error, "islanded: bus C cut off, 7.0 MW load disconnected")
      unit_cut = only(r for r in res if r.error !== nothing && occursin("bus G cut off", r.error))
      @test occursin("1 generating unit(s) with 4.0 MW out of service", unit_cut.error)
      @test unit_cut.shed_load_mw == 0.0
      # the two parallel feeders cut off nothing
      @test count(r -> r.error === nothing, res) == 2
      report = buildContingencyReport(res)
      @test report.n_islanded == 2
      @test report.total_shed_load_mw == 7.0
      # the outage of the feeding unit removes the only reference: with
      # the unit at G takes over and the row names its bus; with
      # auto_slack = false the case ends on the missing reference
      units = generateN1Generators(net)
      handed = redirect_stdout(devnull) do
        runContingencies!(net, units; parallel_enabled = false)
      end
      taken = [r for r in handed if r.error !== nothing && occursin("reference taken over by bus G", r.error)]
      @test length(taken) == 1
      @test taken[1].converged
      @test !occursin("island", taken[1].error)
      plain = runContingencies!(net, units; parallel_enabled = false, auto_slack = false)
      @test any(r -> !r.converged && occursin("no slack bus", r.error), plain)

      # a three-winding transformer (star equivalent) trips as a whole: one
      # case named after it, no single-winding case unless asked for; with
      # the transformer out the isolated star point is dropped silently and
      # only the real bus behind it (B4, 2 MW) is reported
      star = three_winding_star_net()
      grp = only(Sparlectra._three_winding_groups(star))
      legs = Set(getCompName(star.branchVec[k].comp) for k in grp.legs)
      cases3 = generateN1Branches(star)
      t3 = only(c for c in cases3 if c.kind === :transformer3w)
      @test t3.name == "B2_B3_B4"
      @test !any(c -> only(c.elements) in legs, cases3)
      @test count(c -> only(c.elements) in legs, generateN1Branches(star; three_winding_legs = true)) == 3
      r3 = only(redirect_stdout(() -> runContingencies!(star, [t3]; parallel_enabled = false), devnull))
      @test r3.converged
      @test r3.shed_load_mw == 2.0
      @test r3.error == "islanded: bus B4 cut off, 2.0 MW load disconnected; the remaining network solved"
      # the star point never counts for the voltage band: an upper limit just
      # below its solved voltage flags the visible buses above it, not it
      runpf!(star, 30, 1e-8, 0)
      vstar = star.nodeVec[grp.star]._vm_pu
      m = Sparlectra._contingency_metrics(star, 0.9, vstar - 1e-6, fill(NaN, length(star.branchVec)), contingency_names(star))
      @test !("Aux3WT_B2_B3_B4" in m.violations)
      @test !isempty(m.violations)
      @test m.vmin == minimum(star.nodeVec[k]._vm_pu for k in eachindex(star.nodeVec) if k != grp.star)

      # a bus behind a closed link (C at B) is solved with its link group:
      # in the ring no branch outage cuts it off, opening the link does;
      # in the radial net the outage of A-B cuts off B and C together. The
      # import types C Isolated (branches only), which hid C's load from
      # every cut-off and counted C as cut off nowhere.
      ring = link_group_net(; ring = true)
      rres = runContingencies!(ring, vcat(generateN1Branches(ring), generateLinkOutages(ring)); parallel_enabled = false)
      @test all(r -> r.converged, rres)
      @test [(r.kind, r.shed_load_mw, r.error) for r in rres] == [(:branch, 0.0, nothing), (:branch, 0.0, nothing), (:branch, 0.0, nothing),
                                                                  (:link, 5.0, "islanded: bus C cut off, 5.0 MW load disconnected; the remaining network solved")]
      radial = link_group_net(; ring = false)
      # branch order: A-B, A-D
      @test [(r.shed_load_mw, r.error) for r in runContingencies!(radial, generateN1Branches(radial); parallel_enabled = false)] ==
            [(15.0, "islanded: buses B, C cut off, 15.0 MW load disconnected; the remaining network solved"),
             (20.0, "islanded: bus D cut off, 20.0 MW load disconnected; the remaining network solved")]
    end)() end

  end)() end
end

# The N-1 batch on worker threads against the serial batch.
# run_contingency_tests runs it in the pf profile; the threads profile
# (tools/run_gates.sh threads) repeats it with four threads.
function run_contingency_parallel_tests()
  @testset "sp_case14 full N-1: serial and parallel identical" begin (function ()
    net = load_fixture_net("sp_case14")
    before_branches = length(net.branchVec)
    before_vm = [n._vm_pu for n in net.nodeVec]
    cases = generateN1Branches(net)
    serial = runContingencies!(net, cases; parallel_enabled = false)
    @test length(serial) == length(cases)
    # sp_case14 islands its load-only feeder on three N-1 outages by
    # design (the demo fixture pins 14 of 17 converged)
    @test count(r -> r.converged, serial) >= length(cases) - 3
    for max_tasks in (1, Threads.nthreads())
      par = runContingencies!(net, cases; parallel_enabled = true, parallel_max_tasks = max_tasks, parallel_min_work_items = 2)
      @test all(_contingency_results_equal(serial[i], par[i]) for i in eachindex(serial))
    end
    # progress callback (Web UI N-1 counter): one call per
    # finished outage, serial and parallel, each value of 1:total exactly
    # once with the right total, and the results unchanged by it.
    # Repeated 5x per path (threads rule: a single parallel run can pass
    # by luck); the callback runs on worker threads, so it locks.
    for parallel in (false, true), rep in 1:5
      seen = Tuple{Int,Int}[]
      seen_lock = ReentrantLock()
      cb = (done, total) -> lock(() -> push!(seen, (done, total)), seen_lock)
      with_progress = runContingencies!(net, cases; parallel_enabled = parallel, parallel_max_tasks = Threads.nthreads(), parallel_min_work_items = 2, progress = cb)
      @test length(seen) == length(cases)
      @test all(t -> t[2] == length(cases), seen)
      @test sort!(first.(seen)) == collect(1:length(cases))
      @test all(_contingency_results_equal(serial[i], with_progress[i]) for i in eachindex(serial))
    end
    # the base net is never mutated: same branch count, same voltages
    @test length(net.branchVec) == before_branches
    @test isequal([n._vm_pu for n in net.nodeVec], before_vm)
    parallel_comparison_report("N-1 batch on sp_case14 (5 repeats with progress)")
  end)() end

  @testset "generator outages: serial and parallel identical" begin (function ()
    net = load_fixture_net("sp_case14")
    gcases = generateN1Generators(net)
    gs = runContingencies!(net, gcases; parallel_enabled = false)
    gp = runContingencies!(net, gcases; parallel_enabled = true, parallel_max_tasks = 2, parallel_min_work_items = 2)
    @test all(_contingency_results_equal(gs[i], gp[i]) for i in eachindex(gs))
    parallel_comparison_report("N-1 generator outages on sp_case14")
  end)() end

  @testset "multi-stage ladder: serial and parallel identical" begin (function ()
    net = load_fixture_net("sp_case14")
    cases = generateN1Branches(net)
    serial = runContingencies!(net, cases; rescue_ladder = [:warm, :dc], parallel_enabled = false)
    par = runContingencies!(net, cases; rescue_ladder = [:warm, :dc], parallel_enabled = true, parallel_max_tasks = 2, parallel_min_work_items = 2)
    @test all(_contingency_results_equal(serial[i], par[i]) for i in eachindex(serial))
    parallel_comparison_report("N-1 with ladder [:warm, :dc] on sp_case14")
  end)() end
end

# Extended profile: a small N-1 identity slice on the real case1354pegase
# on the shipped sp_case1354. Serial vs parallel
# field equality on a branch subset with max_tasks = 1 and auto, plus a
# flat stage in the ladder (must not change converged cases).
function run_contingency_extended_tests()
  @testset "N-1 identity on sp_case1354 (extended)" begin (function ()
    # the shipped data/mpower/sp_case1354.m (synthetic 1354-bus grid, 2040
    # branches), so the identity check runs on every checkout
    case_path = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "sp_case1354.m"))
    @test isfile(case_path)
    net = createNetFromMatPowerFile(filename = case_path, flatstart = false, enable_pq_gen_controllers = true, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)
    cases = generateN1Branches(net)[1:60]
    serial = runContingencies!(net, cases; parallel_enabled = false)
    @test length(serial) == 60
    for max_tasks in (1, Threads.nthreads())
      par = runContingencies!(net, cases; parallel_enabled = true, parallel_max_tasks = max_tasks, parallel_min_work_items = 2)
      @test all(_contingency_results_equal(serial[i], par[i]) for i in eachindex(serial))
    end
    retried = runContingencies!(net, cases; parallel_enabled = true, parallel_min_work_items = 2, rescue_ladder = [:warm, :flat])
    @test all(!serial[i].converged || _contingency_results_equal(serial[i], retried[i]) for i in eachindex(serial))
    println("contingency extended case1354: RAN (", count(r -> r.converged, serial), "/60 converged, threads = ", Threads.nthreads(), ")")
    parallel_comparison_report("N-1 slice of 60 outages on sp_case1354")
  end)() end
end

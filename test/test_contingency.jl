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
      @test results[1].error == "islanded without reference"
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
      @test startswith(lines[1], "name;converged;iterations")
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
    net = createNetFromMatPowerFile(filename = case_path, matpower_shift_unit = :rad, matpower_shift_sign = -1.0)
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

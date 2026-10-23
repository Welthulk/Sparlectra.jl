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

# file: test/test_parallel_foundation.jl
# purpose: thread-safety foundation of the parallel execution paths: the
#          per-worker performance-profile child/merge helpers, the condest
#          thunk on a copied Net, and the island fan-out against the serial
#          island loop.

using Sparlectra
using Test

function run_parallel_foundation_tests()
  @testset "Parallel foundation (Phase 1)" begin (function ()
    @testset "performance-profile child and merge" begin (function ()
      recorder_calls = Ref(0)
      parent = Dict{Symbol,Any}(:enabled => true, :show_allocations => false, :output_dir => "/tmp/x", :phase_callback => _ -> (recorder_calls[] += 1), :cancellation_check => () -> nothing)

      @test Sparlectra._perf_profile_child(nothing) === nothing
      @test Sparlectra._perf_profile_child(Dict{Symbol,Any}(:enabled => false)) === nothing

      child = Sparlectra._perf_profile_child(parent)
      @test child isa Dict{Symbol,Any}
      @test child[:enabled] === true
      @test child[:output_dir] == "/tmp/x"
      @test haskey(child, :cancellation_check)
      # orchestrator-only: the recorder callback must NOT reach workers
      @test !haskey(child, :phase_callback)
      @test !haskey(child, :diagnostic_artifact_prefix)

      # worker writes: timings, iteration rows, per-worker scalars
      Sparlectra._perf_profile_add!(child, :ybus_assembly, 0.25, 100)
      Sparlectra._perf_profile_add!(child, :ybus_assembly, 0.25, 100)
      Sparlectra._perf_profile_push_iteration!(child, (iter = 1, mismatch = 0.1))
      child[:linear_solver_backend] = :umfpack_reuse

      # parent already carries one timing of the same phase (serial part)
      Sparlectra._perf_profile_add!(parent, :ybus_assembly, 1.0, 50)
      Sparlectra._perf_profile_merge!(parent, child, "ac_island_2_")

      row = parent[:timings][:ybus_assembly]
      @test row.calls == 3
      @test isapprox(row.elapsed_s, 1.5; atol = 1e-12)
      @test row.bytes == 250
      @test length(parent[:iterations]) == 1
      @test parent[:ac_island_2_linear_solver_backend] === :umfpack_reuse
      @test !haskey(parent, :linear_solver_backend)
      @test recorder_calls[] == 0

      # the orchestrating fan-out accounts its wall clock separately
      Sparlectra._perf_profile_add!(parent, :parallel_wall_time, 0.75, 0)
      @test parent[:timings][:parallel_wall_time].elapsed_s == 0.75

      # disabled parent: merge is a no-op
      off = Dict{Symbol,Any}(:enabled => false)
      @test Sparlectra._perf_profile_merge!(off, child, "x_") === off
      @test !haskey(off, :timings)

      # Phase 2 contract: call sites never branch on the disabled case, the
      # helpers swallow nothing children (and a nothing parent) themselves
      before = deepcopy(parent)
      @test Sparlectra._perf_profile_merge!(parent, nothing, "ac_island_9_") === parent
      @test length(parent[:timings]) == length(before[:timings])
      @test Sparlectra._perf_profile_merge!(nothing, child, "x_") === nothing
    end)() end

    @testset "deepcopy carries the condest thunk" begin (function ()
      # the solver status stores a lazy condition-estimate closure; a copied
      # Net must keep a WORKING thunk (its own memo Ref), not a dead one
      net = createTest3BusNet()
      _, erg = runpf!(net, 20, 1e-8, 0)
      @test erg == 0
      copied = deepcopy(net)
      # evaluate on the copy FIRST (fresh memo) and then on the original
      kappa_copy = condestJacobian(copied)
      kappa_orig = condestJacobian(net)
      @test isfinite(kappa_copy) && kappa_copy > 0
      @test kappa_copy == kappa_orig
      # second evaluation hits the memo and stays identical
      @test condestJacobian(copied) == kappa_copy
    end)() end

    run_parallel_island_tests()
  end)() end
end

# The island fan-out of the power flow against the serial island loop.
# run_parallel_foundation_tests runs it in the pf profile; the threads profile
# (tools/run_gates.sh threads) repeats it with four threads.
function run_parallel_island_tests()
  @testset "two islands of 300 buses: parallel equals serial (0.30.2)" begin (function ()
    # two copies of the shipped sp_case300 in one net (bus numbers of the
    # second shifted by 1000, names suffixed): two islands of 300 buses,
    # each with its own slack, above the default threshold of 200 buses
    mpc = Sparlectra.MatpowerIO.read_case(joinpath(dirname(@__DIR__), "data", "mpower", "sp_case300.m"); legacy_compat = true)
    off = 1000.0
    bus2 = copy(mpc.bus); bus2[:, 1] .+= off
    gen2 = copy(mpc.gen); gen2[:, 1] .+= off
    br2 = copy(mpc.branch); br2[:, 1] .+= off; br2[:, 2] .+= off
    names2 = mpc.bus_name === nothing ? nothing : vcat(mpc.bus_name, [n * "_B" for n in mpc.bus_name])
    two = Sparlectra.MatpowerIO.MatpowerCase("two_islands_300", mpc.baseMVA, vcat(mpc.bus, bus2), vcat(mpc.gen, gen2), vcat(mpc.branch, br2), nothing, names2, nothing, nothing, nothing, nothing, nothing)
    solve = function (parallel::Bool)
      net = Sparlectra.createNetFromMatPowerCase(mpc = two)
      profile = Dict{Symbol,Any}(:enabled => true)
      it, erg = runpf!(net, 30, 1e-8, 0; islands_enabled = true, islands_parallel_enabled = parallel, performance_profile = profile)
      return (; net, profile, it, erg)
    end
    serial = solve(false)
    @test serial.erg == 0
    @test length(Sparlectra.detect_ac_islands(serial.net).rows) == 2
    if Threads.nthreads() > 1
      for _ = 1:3
        par = solve(true)
        @test haskey(par.profile[:timings], :parallel_wall_time)
        @test par.erg == 0 && par.it == serial.it
        @test all(par.net.nodeVec[i]._vm_pu == serial.net.nodeVec[i]._vm_pu && par.net.nodeVec[i]._va_deg == serial.net.nodeVec[i]._va_deg for i in eachindex(serial.net.nodeVec))
        @test length(par.net.qLimitLog) == length(serial.net.qLimitLog)
        @test Set(keys(par.profile[:ac_island_solver_statuses])) == Set(keys(serial.profile[:ac_island_solver_statuses]))
        @test all(par.profile[:ac_island_solver_statuses][k].status === serial.profile[:ac_island_solver_statuses][k].status for k in keys(serial.profile[:ac_island_solver_statuses]))
      end
    end
    parallel_comparison_report("power flow, two islands of 300 buses (3 repeats)")
  end)() end

  @testset "island parallel identity (Phase 2; by size since 0.30.2)" begin (function ()
    # multi-island fixture: n disconnected 3-bus feeders, distinct loads
    function build_multi_island(n)
      net = Net(name = "mi_$(n)", baseMVA = 100.0)
      for k in 1:n
        a = "A$(k)"; b = "B$(k)"; c = "C$(k)"
        for bus in (a, b, c)
          addBus!(net = net, busName = bus, vn_kV = 110.0)
        end
        addProsumer!(net = net, busName = a, type = "EXTERNALNETWORKINJECTION", referencePri = a, vm_pu = 1.0, va_deg = 0.0)
        addProsumer!(net = net, busName = b, type = "ENERGYCONSUMER", p = 10.0 + k, q = 3.0)
        addProsumer!(net = net, busName = c, type = "ENERGYCONSUMER", p = 5.0 + k, q = 1.0)
        addPIModelACLine!(net = net, fromBus = a, toBus = b, r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = net, fromBus = b, toBus = c, r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      end
      validate!(net = net)
      return net
    end

    # 0.30.2: there is no serial mode to choose; the serial reference is
    # runtime.parallel off (the master switch), the parallel run lowers
    # the island size threshold to 1 (the fixture's islands have 3 buses)
    solve_mode = function (mode; max_tasks = Threads.nthreads())
      net = build_multi_island(5)
      profile = Dict{Symbol,Any}(:enabled => true)
      it, erg = runpf!(net, 25, 1e-8, 0; islands_enabled = true, islands_parallel_enabled = mode === :parallel, islands_parallel_min_buses = 1, islands_parallel_max_tasks = max_tasks, performance_profile = profile)
      return (; net, profile, it, erg)
    end

    serial = solve_mode(:serial)
    # BITWISE identity required (Phase 0 review item 5), for max_tasks=1
    # (falls back to the serial loop, same functions) and max_tasks=auto
    for max_tasks in (1, Threads.nthreads())
      par = solve_mode(:parallel; max_tasks = max_tasks)
      @test par.erg == 0
      @test par.it == serial.it
      for i in eachindex(serial.net.nodeVec)
        @test par.net.nodeVec[i]._vm_pu == serial.net.nodeVec[i]._vm_pu
        @test par.net.nodeVec[i]._va_deg == serial.net.nodeVec[i]._va_deg
      end
      # phase-name set identical to serial (parallel adds only the wall clock)
      par_phases = Set(keys(par.profile[:timings]))
      delete!(par_phases, :parallel_wall_time)
      @test par_phases == Set(keys(serial.profile[:timings]))
      statuses = par.profile[:ac_island_solver_statuses]
      @test length(statuses) == 5
      @test all(st.status === :converged for (_, st) in statuses)
    end

    if Threads.nthreads() > 1
      # the parallel path actually engaged: wall clock recorded, per-island
      # scalars merged under their island prefix, no shared-slot leftovers
      par = solve_mode(:parallel)
      @test haskey(par.profile[:timings], :parallel_wall_time)
      @test haskey(par.profile, :ac_island_1_linear_solver_backend)
      @test !haskey(par.profile, :diagnostic_artifact_prefix)
    end
    parallel_comparison_report("power flow, five islands of 3 buses")

    # failure semantics under parallel: all islands report their REAL
    # status even with diagnostic_continue_after_failure = false
    net_bad = build_multi_island(4)
    # island 2 gets an absurd load so its solve fails
    addProsumer!(net = net_bad, busName = "B2", type = "ENERGYCONSUMER", p = 1.0e6, q = 1.0e6)
    profile_bad = Dict{Symbol,Any}(:enabled => true)
    @test_throws Exception runpf!(net_bad, 15, 1e-8, 0; islands_enabled = true, islands_parallel_enabled = true, islands_parallel_min_buses = 1, islands_parallel_max_tasks = Threads.nthreads(), performance_profile = profile_bad)
    statuses_bad = profile_bad[:ac_island_solver_statuses]
    if Threads.nthreads() > 1
      # parallel: every island already ran, all four report their REAL status
      @test length(statuses_bad) == 4
      @test all(st.status !== :skipped_after_previous_failure for (_, st) in statuses_bad)
      @test any(st.status !== :converged for (_, st) in statuses_bad)
      @test count(st.status === :converged for (_, st) in statuses_bad) == 3
    else
      # single-threaded fallback = serial loop: stops at the failing island,
      # so exactly the islands up to and including the failure are recorded
      @test length(statuses_bad) == 2
    end

    # Q limits per island (0.30.2): each island net got the per-bus limit
    # vectors of the whole net, so an island solve read the limits of
    # buses 1..n of the net (rectangular) or dropped them for the length
    # mismatch (APSLF). Two islands of three buses: island A carries a
    # machine with a wide band, island B the machine B2 that clamps at
    # Qmax when B is solved alone. The two-island run must clamp exactly
    # as the single-island runs, with both solvers; the rectangular
    # solver also on its parallel path (APSLF solves islands serially).
    function qlimit_islands(parts)
      net = Net(name = "qlimit_islands", baseMVA = 100.0)
      if :A in parts
        for b in ("A1", "A2", "A3")
          addBus!(net = net, busName = b, vn_kV = 110.0, vm_pu = 1.02)
        end
        addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "A1")
        addProsumer!(net = net, busName = "A2", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 0.0, vm_pu = 1.01, qMin = -100.0, qMax = 100.0)
        addProsumer!(net = net, busName = "A3", type = "ENERGYCONSUMER", p = 40.0, q = 12.0)
        addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = net, fromBus = "A2", toBus = "A3", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      end
      if :B in parts
        for b in ("B1", "B2", "B3")
          addBus!(net = net, busName = b, vn_kV = 110.0, vm_pu = 1.02)
        end
        addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "B1")
        addProsumer!(net = net, busName = "B2", type = "SYNCHRONOUSMACHINE", p = 10.0, q = 0.0, vm_pu = 1.05, qMin = -10.0, qMax = 10.0)
        addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 80.0, q = 25.0)
        addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.02, x_pu = 0.12, b_pu = 0.0, status = 1)
      end
      return net
    end
    quiet = OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)
    qsolve = function (solver::Symbol, parts; parallel::Bool = false, lock = Int[])
      net = qlimit_islands(parts)
      if solver === :rectangular
        _, erg = runpf!(net, 30, 1e-10, 0; islands_enabled = true, islands_parallel_enabled = parallel, islands_parallel_min_buses = 1, islands_parallel_max_tasks = Threads.nthreads(), lock_pv_to_pq_buses = lock)
        @test erg == 0
      else
        r = run_sparlectra(net = net, config = SparlectraConfig(powerflow = PowerFlowConfig(solver = :apslf, tol = 1e-10, islands = Sparlectra.IslandPowerFlowConfig(enabled = true)), output = quiet))
        @test r.final_converged
      end
      return net
    end
    clamped(net) = Dict(name => net.qLimitEvents[b] for (name, b) in net.busDict if haskey(net.qLimitEvents, b))
    vm_of(net, names) = [net.nodeVec[net.busDict[n]]._vm_pu for n in names]
    for solver in (:rectangular, :apslf)
      only_a = qsolve(solver, (:A,))
      only_b = qsolve(solver, (:B,))
      @test clamped(only_b) == Dict("B2" => :max)
      @test isempty(clamped(only_a))
      runs = solver === :rectangular && Threads.nthreads() > 1 ? vcat([false], fill(true, 5)) : [false]
      for parallel in runs
        both = qsolve(solver, (:A, :B); parallel = parallel)
        @test (solver, parallel, clamped(both)) == (solver, parallel, Dict("B2" => :max))
        @test length(both.qLimitLog) == length(only_a.qLimitLog) + length(only_b.qLimitLog)
        @test isapprox(vm_of(both, ("B1", "B2", "B3")), vm_of(only_b, ("B1", "B2", "B3")); atol = 1e-9)
        @test isapprox(vm_of(both, ("A1", "A2", "A3")), vm_of(only_a, ("A1", "A2", "A3")); atol = 1e-9)
      end
    end
    # APSLF solves islands serially; only the rectangular solver runs parallel
    parallel_comparison_report("Q limits per island, rectangular solver (5 repeats)")
    # bus lists of internal positions go to the islands in their own
    # numbering: locking A2 (position 2 of the net) must not lock B2,
    # position 2 of the second island
    locked = qsolve(:rectangular, (:A, :B); lock = [2])
    @test isapprox(vm_of(locked, ("B1", "B2", "B3")), vm_of(qsolve(:rectangular, (:B,)), ("B1", "B2", "B3")); atol = 1e-9)
    # and A2 itself runs as PQ (#458), converted inside its island at
    # iteration 0: its voltage leaves the 1.01 pu setpoint, and the
    # conversion comes back to the whole net as a Q-limit entry (the
    # island path carries switched buses as entries and voltages, not as
    # node types; the same holds for B2's clamp)
    a2 = locked.busDict["A2"]
    @test abs(locked.nodeVec[a2]._vm_pu - 1.01) > 1e-3
    @test any(e -> e.bus == a2 && e.iter == 0, locked.qLimitLog)
    @test isapprox(vm_of(qsolve(:rectangular, (:A, :B)), ("A2",))[1], 1.01; atol = 1e-9)
  end)() end
end

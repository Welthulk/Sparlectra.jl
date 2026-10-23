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

# file: test/test_legacy_solver_interface.jl
# purpose: legacy profile: the testsets of test/test_solver_interface.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_solver_interface.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_solver_interface_tests()
  # cut from testset: Solver interface
    @testset "Wrong-branch helper classification" begin (function ()
      @test isapprox(Sparlectra._circular_angle_spread_deg([179.0, -179.0]), 2.0; atol = 1e-8)
      @test isapprox(Sparlectra._circular_angle_spread_deg([-170.0, 170.0]), 20.0; atol = 1e-8)
      Vok = ComplexF64[1.0 + 0im, 1.01 + 0.01im, 1.03 - 0.02im]
      bus_types = [:PQ, :PV, :Slack]
      vset = [1.0, 1.01, 1.03]
      net = createTest3BusNet()
      ok = Sparlectra._check_wrong_branch_solution(net, Vok, bus_types, vset, 3; min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 180.0, max_branch_angle_deg = 90.0, min_low_vm_count = 1)
      @test ok.status == :ok

      Vlow = ComplexF64[0.40 + 0im, 1.0 + 0im, 1.03 + 0im]
      low = Sparlectra._check_wrong_branch_solution(net, Vlow, bus_types, vset, 3; min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 180.0, max_branch_angle_deg = 90.0, min_low_vm_count = 1)
      @test low.status == :warn
      @test low.reason == :low_voltage_magnitude

      Vbad = ComplexF64[NaN + 0im, 1.0 + 0im, 1.03 + 0im]
      bad = Sparlectra._check_wrong_branch_solution(net, Vbad, bus_types, vset, 3; min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 180.0, max_branch_angle_deg = 90.0, min_low_vm_count = 1)
      @test bad.status == :fail
      @test bad.reason == :nonfinite_voltage

      branch_net = Net(name = "branch_wrap_test", baseMVA = 100.0)
      addBus!(net = branch_net, busName = "B1", vn_kV = 110.0)
      addBus!(net = branch_net, busName = "B2", vn_kV = 110.0)
      addBus!(net = branch_net, busName = "B3", vn_kV = 110.0)
      addPIModelACLine!(net = branch_net, fromBus = "B1", toBus = "B2", r_pu = 0.0, x_pu = 0.1, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = branch_net, fromBus = "B2", toBus = "B3", r_pu = 0.0, x_pu = 0.1, b_pu = 0.0, status = 1)
      setNetBranchStatus!(net = branch_net, branchNr = 2, status = 0)
      branch_net.branchVec[1].phase_shift_deg = 10.0
      branch_net.branchVec[1].angle = 10.0

      Vwrap = ComplexF64[exp(im * deg2rad(179.0)), exp(im * deg2rad(-179.0)), exp(im * deg2rad(0.0))]
      wrap = Sparlectra._check_wrong_branch_solution(Vwrap, [:PQ, :PQ, :Slack], [1.0, 1.0, 1.0], 3; net = branch_net, min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 360.0, max_branch_angle_deg = 20.0, min_low_vm_count = 1)
      @test wrap.status == :ok
      @test wrap.branch_angle_violation_count == 0
      @test !isnothing(wrap.worst_branch)
      @test wrap.worst_branch.angle_check_basis == :effective_bus_angle_minus_phase_shift
      @test isapprox(wrap.worst_branch.angle_diff_raw_deg, 2.0; atol = 1e-8)
      @test isapprox(wrap.worst_branch.angle_diff_effective_deg, 12.0; atol = 1e-8)

      strict = Sparlectra._check_wrong_branch_solution(Vwrap, [:PQ, :PQ, :Slack], [1.0, 1.0, 1.0], 3; net = branch_net, min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 360.0, max_branch_angle_deg = 11.0, min_low_vm_count = 1)
      @test strict.status == :warn
      @test strict.reason == :branch_angle_exceeded
      @test strict.branch_angle_violation_count == 1

      wrap_spread = Sparlectra._check_wrong_branch_solution(Vwrap, [:PQ, :PQ, :Slack], [1.0, 1.0, 1.0], 3; net = nothing, min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 10.0, max_branch_angle_deg = Inf, min_low_vm_count = 1)
      @test wrap_spread.status == :warn
      @test wrap_spread.reason == :angle_spread_exceeded
      @test isapprox(wrap_spread.angle_spread_deg, 181.0; atol = 1e-8)

      Vwide = ComplexF64[exp(im * deg2rad(120.0)), exp(im * deg2rad(-120.0)), exp(im * deg2rad(0.0))]
      wide_spread = Sparlectra._check_wrong_branch_solution(Vwide, [:PQ, :PQ, :Slack], [1.0, 1.0, 1.0], 3; net = nothing, min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 30.0, max_branch_angle_deg = Inf, min_low_vm_count = 1)
      @test wide_spread.status == :warn
      @test wide_spread.reason == :angle_spread_exceeded
      @test isapprox(wide_spread.angle_spread_deg, 240.0; atol = 1e-8)

      # Highest-voltage-level scope: a voltage dip below the top level is
      # normal operating spread and must not trigger (regression: healthy
      # 45 kV feeders of a real CGMES snapshot were flagged SUSPECT while the
      # 380 kV level was clean); the same dip on the top level still triggers.
      mixed_net = Net(name = "wrong_branch_levels", baseMVA = 100.0)
      addBus!(net = mixed_net, busName = "HV1", vn_kV = 380.0)
      addBus!(net = mixed_net, busName = "HV2", vn_kV = 380.0)
      addBus!(net = mixed_net, busName = "MV1", vn_kV = 45.0)
      mixed_types = [:Slack, :PQ, :PQ]
      mixed_vset = [1.0, 1.0, 1.0]
      lower_dip = Sparlectra._check_wrong_branch_solution(mixed_net, ComplexF64[1.0, 1.01, 0.90], mixed_types, mixed_vset, 1; min_vm_pu = 0.95, max_vm_pu = 1.30, max_angle_spread_deg = 180.0, max_branch_angle_deg = Inf, min_low_vm_count = 1)
      @test lower_dip.status == :ok
      @test lower_dip.low_vm_count == 0
      # reported band refers to the checked (top) level only
      @test isapprox(lower_dip.min_vm_pu, 1.0; atol = 1e-12)
      top_dip = Sparlectra._check_wrong_branch_solution(mixed_net, ComplexF64[1.0, 0.90, 1.0], mixed_types, mixed_vset, 1; min_vm_pu = 0.95, max_vm_pu = 1.30, max_angle_spread_deg = 180.0, max_branch_angle_deg = Inf, min_low_vm_count = 1)
      @test top_dip.status == :warn
      @test top_dip.reason == :low_voltage_magnitude
      @test top_dip.low_vm_count == 1
    end)() end

    @testset "Q-limit guard and active-set status" begin (function ()
      net = createTest3BusNet()
      bus = geNetBusIdx(net = net, busName = "STATION1")
      qmin_pu = [-Inf, 0.0, -Inf]
      qmax_pu = [Inf, 0.0, Inf]
      bus_types = [:PQ, :PV, :Slack]
      S = ComplexF64[0.0 + 0.0im, 0.1 + 0.0im, 0.0 + 0.0im]
      guarded = Sparlectra._apply_qlimit_guard_to_rectangular_active_set!(
        net,
        bus_types,
        S,
        zeros(Float64, 3),
        qmin_pu,
        qmax_pu;
        min_q_range_pu = 1e-4,
        zero_range_mode = :lock_pq,
        narrow_range_mode = :prefer_pq,
        log = false,
        verbose = 0,
      )
      @test guarded == [bus]
      @test bus_types[bus] == :PQ
      @test length(net.qLimitLog) == 1

      # The locked bus starts from the voltage its own fixed injection gives
      # against its neighbours, not from the machine setpoint: its bus
      # equation holds at the start, every other bus keeps its start value.
      Y3 = Sparlectra.createYBUS(net = net, sparse = true, printYBUS = false)
      V3 = ComplexF64[1.02, 1.04, 1.0]
      V3_before = copy(V3)
      sweeps = Sparlectra._start_guarded_buses_from_fixed_injection!(V3, Y3, S, guarded)
      @test sweeps >= 1
      @test abs(V3[bus] * conj((Y3 * V3)[bus]) - S[bus]) < 1e-9
      @test all(V3[k] == V3_before[k] for k in eachindex(V3) if k != bus)

      # Regression (case_SyntheticUSA island 1 did not converge with the
      # guard on): every PV machine of sp_case1354 with a zero Q band at 0
      # and a setpoint of 1.06 pu it cannot hold, its bus voltage at that
      # setpoint (as the importer starts a PV bus). Locked at the setpoint
      # Newton stopped at the iteration limit, damped or not; started from
      # the fixed injection it converges in 7 iterations and every guarded
      # bus ends PQ with zero generator Q. The run uses autodamp = true, the
      # template default, on purpose: with 259 of 260 machines locked the
      # network is almost all PQ and undamped Newton fails here from either
      # start. With the bus voltages left at 1.0 pu below the setpoint (a
      # state the importer does not produce) the old start converged
      # undamped in 10 iterations and the new one does not, because the
      # corrected magnitudes make the first full step larger; damped, both
      # starts converge there (measured).
      zb_net = Sparlectra.createNetFromMatPowerFile(filename = joinpath(dirname(@__DIR__), "data", "mpower", "sp_case1354.m"))
      zb_pv = Set(k for (k, n) in enumerate(zb_net.nodeVec) if Sparlectra.getNodeType(n) == Sparlectra.PV)
      for ps in zb_net.prosumpsVec
        (Sparlectra.isGenerator(ps) && ps.comp.cFrom_bus in zb_pv) || continue
        ps.minQ = 0.0
        ps.maxQ = 0.0
        ps.vm_pu = 1.06
      end
      foreach(k -> Sparlectra.setVmVa!(node = zb_net.nodeVec[k], vm_pu = 1.06), zb_pv)
      Sparlectra.buildQLimits!(zb_net)
      _, zb_erg = runpf!(zb_net, 30, 1e-8, 0; qlimits_enabled = true, qlimit_guard = true, qlimit_guard_min_q_range_pu = 0.02, autodamp = true)
      @test zb_erg == 0
      @test Sparlectra.rectangular_pf_status(zb_net).guarded_narrow_q_pv_buses == length(zb_pv)
      @test all(Sparlectra.getNodeType(zb_net.nodeVec[k]) == Sparlectra.PQ for k in zb_pv)
      @test all(abs(something(zb_net.nodeVec[k]._qƩGen, NaN)) < 1e-9 for k in zb_pv)

      bus_types[bus] = :PV
      changed, reenabled = Sparlectra.active_set_q_limits!(
        net,
        2,
        3;
        get_qreq_pu = _ -> 0.2,
        is_pv = b -> (bus_types[b] == :PV),
        make_pq! = (b, _qclamp, _side) -> (bus_types[b] = :PQ),
        make_pv! = b -> (bus_types[b] = :PV),
        qmin_pu = qmin_pu,
        qmax_pu = qmax_pu,
        pv_orig_mask = trues(3),
        allow_reenable = false,
        q_hyst_pu = 0.0,
        cooldown_iters = 0,
        qlimit_guard_max_switches = 1,
        qlimit_guard_freeze_after_repeated_switching = true,
      )
      @test !changed
      @test !reenabled
      @test bus_types[bus] == :PV
      @test length(net.qLimitLog) == 1

      function zero_range_pv_net()
        guarded_net = createTest3BusNet()
        guarded_bus = geNetBusIdx(net = guarded_net, busName = "STATION1")
        guarded_net.nodeVec[guarded_bus]._nodeType = Sparlectra.PV
        empty!(guarded_net.qmin_pu)
        append!(guarded_net.qmin_pu, [-Inf, 0.0, -Inf])
        empty!(guarded_net.qmax_pu)
        append!(guarded_net.qmax_pu, [Inf, 0.0, Inf])
        return guarded_net
      end

      # switching the guard off leaves the zero-headroom bus unlocked
      opt_out_net = zero_range_pv_net()
      redirect_stdout(devnull) do
        runpf!(opt_out_net; config = PowerFlowConfig(max_iter = 0, qlimits = QLimitConfig(guard = false)))
      end
      @test isempty(opt_out_net.qLimitLog)

      # a run stopped by the switching cap names the key, its value and the
      # buses that reached it (the guard's lock counts as the bus's first
      # switch, so a cap of 1 stops the run)
      capped_net = zero_range_pv_net()
      redirect_stdout(devnull) do
        runpf!(capped_net; config = PowerFlowConfig(qlimits = QLimitConfig(guard = true, guard_max_switches = 1)))
      end
      capped = Sparlectra.rectangular_pf_status(capped_net)
      @test capped.reason === :max_switching_exceeded
      @test occursin("power_flow.qlimits.guard.max_switches = 1", capped.reason_text)
      @test occursin("STATION1 (1x)", capped.reason_text)
    end)() end

    @testset "Independently solved AC islands do not overwrite each other's diagnostic logs" begin (function ()
      # Regression for a real bug: islands share one performance_profile/output_dir,
      # so fixed-name diagnostic artifacts (merit_linesearch.log, trust_region.log,
      # current_iteration_start.log, apslf_start.log) were silently overwritten by
      # whichever island solved last, making per-island diagnosis on a failing
      # multi-island run impossible/misleading.
      net = Net(name = "island_log_prefix_test", baseMVA = 100.0)
      for busName in ("A1", "A2", "B1", "B2")
        addBus!(net = net, busName = busName, vn_kV = 110.0)
      end
      addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
      addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "A1")
      addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
      addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
      addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
      refreshBusTypesFromProsumers!(net)

      mktempdir() do outdir
        profile = Dict{Symbol,Any}(:output_dir => outdir)
        cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 40, tol = 1e-8, islands_enabled = true, autodamp = true, merit = Sparlectra.MeritLineSearchConfig(enabled = true)))
        result = run_sparlectra(net = net, config = cfg, performance_profile = profile)
        @test result.final_converged

        @test !isfile(joinpath(outdir, "merit_linesearch.log"))
        @test isfile(joinpath(outdir, "ac_island_1_merit_linesearch.log"))
        @test isfile(joinpath(outdir, "ac_island_2_merit_linesearch.log"))

        log1 = read(joinpath(outdir, "ac_island_1_merit_linesearch.log"), String)
        log2 = read(joinpath(outdir, "ac_island_2_merit_linesearch.log"), String)
        @test log1 != log2
        @test occursin("merit_initial: 0.005450000000000001", log1)
        @test occursin("merit_initial: 0.03445", log2)

        artifacts = SparlectraApp.collect_sparlectra_api_artifacts(outdir)
        island1_artifact = only(a for a in artifacts if a.name == "ac_island_1_merit_linesearch.log")
        island2_artifact = only(a for a in artifacts if a.name == "ac_island_2_merit_linesearch.log")
        @test island1_artifact.kind === :merit_linesearch
        @test island2_artifact.kind === :merit_linesearch
      end

      # Non-island (single-network) runs keep the unprefixed filename.
      mktempdir() do outdir
        profile = Dict{Symbol,Any}(:output_dir => outdir)
        single_net = createCIGRE()
        _, erg = runpf_rectangular!(single_net, 30, 1e-8, 0; autodamp = true, merit_enabled = true, performance_profile = profile)
        @test erg == 0
        @test isfile(joinpath(outdir, "merit_linesearch.log"))
        @test !any(startswith(f, "ac_island_") for f in readdir(outdir))
      end
    end)() end

    # A solved state can satisfy every limit and still be non-physical: a
    # machine pinned at Qmax whose voltage ends ABOVE its setpoint cannot be
    # an operating point, because a machine at its upper reactive limit has
    # nothing left to raise the voltage with. The counter is what makes the
    # difference between two enforcement modes visible at all.
    @testset "Q-V characteristic check finds non-physical generator states" begin (function ()
      # STATION1 carries the machine with vm_pu = 1.027273 and the Q band
      net = createTest3BusNet(qlim_min = -20.0, qlim_max = 20.0)
      bus = Sparlectra.geNetBusIdx(net = net, busName = "STATION1")

      # the machine sits at Qmax and the voltage ended ABOVE its setpoint
      net.qLimitEvents[bus] = :max
      net.nodeVec[bus]._qƩGen = 20.0
      Sparlectra.setVmVa!(node = net.nodeVec[bus], vm_pu = 1.0300, va_deg = 0.0)
      rows = qvCharacteristicViolations(net)
      @test length(rows) == 1
      @test rows[1].bus == bus
      @test rows[1].side === :max
      @test rows[1].significant
      @test rows[1].dv_pu > 0

      # same clamp, voltage BELOW the setpoint: that is the physical case
      Sparlectra.setVmVa!(node = net.nodeVec[bus], vm_pu = 1.0200, va_deg = 0.0)
      @test isempty(qvCharacteristicViolations(net))

      # voltage on the wrong side again, but the machine is no longer at the
      # limit: a stale switching event is not a violation
      Sparlectra.setVmVa!(node = net.nodeVec[bus], vm_pu = 1.0300, va_deg = 0.0)
      net.nodeVec[bus]._qƩGen = 19.0
      @test isempty(qvCharacteristicViolations(net))

      # within the band around the setpoint the state counts, but not as
      # significant: a converged solve may land on either side by rounding
      net.nodeVec[bus]._qƩGen = 20.0
      Sparlectra.setVmVa!(node = net.nodeVec[bus], vm_pu = 1.027273 + 1.0e-5, va_deg = 0.0)
      band_rows = qvCharacteristicViolations(net)
      @test length(band_rows) == 1
      @test !band_rows[1].significant

      # printing: the section is always there, so a missing line can never be
      # mistaken for a clean result
      Sparlectra.setVmVa!(node = net.nodeVec[bus], vm_pu = 1.0300, va_deg = 0.0)
      io = IOBuffer()
      printQVCharacteristicCheck(net; io = io)
      out = String(take!(io))
      @test occursin("Q-V characteristic:", out)
      @test occursin("non-physical generator state", out)
      io = IOBuffer()
      printQVCharacteristicCheck(net; io = io, converged = false)
      @test occursin("skipped (no converged solution)", String(take!(io)))
      io = IOBuffer()
      empty!(net.qLimitEvents)
      printQVCharacteristicCheck(net; io = io)
      @test occursin("no non-physical generator states", String(take!(io)))

      # The check must read the CLAMPED reactive power, the same value the
      # result print shows as Qg. It read the node sum instead, and the
      # classical outer loop applies its clamp in the solver's specification
      # without writing it back to the node: on the Zeng case the check saw
      # 54.17 MVAr at a 50 MVAr limit and reported nothing while the print
      # showed 50.000. Under every enforcement mode a clamped
      # machine's Q as the check sees it is its limit.
      for mode in (:active_set, :classic_simultaneous, :classic_one_at_a_time)
        fnet = Net(name = "feeder_$(mode)", baseMVA = 100.0)
        addBus!(net = fnet, busName = "B1", vn_kV = 110.0)
        addBus!(net = fnet, busName = "B2", vn_kV = 110.0)
        addBus!(net = fnet, busName = "B3", vn_kV = 110.0)
        addProsumer!(net = fnet, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)
        addPIModelACLine!(net = fnet, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = fnet, fromBus = "B2", toBus = "B3", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
        addProsumer!(net = fnet, busName = "B2", type = "SYNCHRONOUSMACHINE", p = 10.0, q = 0.0, vm_pu = 1.0, qMin = -10.0, qMax = 10.0)
        addProsumer!(net = fnet, busName = "B3", type = "ENERGYCONSUMER", p = 80.0, q = 25.0)
        @test validate!(net = fnet)[1]
        _, ferg = runpf!(fnet, 40, 1e-8, 0; qlimits_enabled = true, qlimit_enforcement_mode = mode)
        @test ferg == 0
        fb2 = Sparlectra.geNetBusIdx(net = fnet, busName = "B2")
        @test get(fnet.qLimitEvents, fb2, nothing) === :max
        @test isapprox(Sparlectra._qv_effective_qgen(fnet, fb2, 10.0), 10.0; atol = 1e-6)

        # issue #374: the RESULT TABLE (not only the Q-V check) must show
        # the same clamped Qg for the PQ* bus, under every enforcement
        # mode; buildACPFlowReport feeds the console table, the API's
        # detailed CSV, and bus_powers.csv, so this one row covers all three.
        report = buildACPFlowReport(fnet)
        b2_row = only(r for r in report.nodes if r.bus == fb2)
        @test isapprox(b2_row.q_gen_MVar, 10.0; atol = 1e-6)
      end

      # issue #374 (part 2): a PV bus that stays PV (never clamped) must
      # also show its SOLVED Qg, not the nameplate q the case declared -
      # the pre-fix code only took the solved _qƩGen when the nameplate
      # summed to exactly zero, so a normal non-zero nameplate silently won.
      pvnet = Net(name = "pv_stays_pv", baseMVA = 100.0)
      addBus!(net = pvnet, busName = "B1", vn_kV = 110.0)
      addBus!(net = pvnet, busName = "B2", vn_kV = 110.0)
      addBus!(net = pvnet, busName = "B3", vn_kV = 110.0)
      addProsumer!(net = pvnet, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)
      addPIModelACLine!(net = pvnet, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = pvnet, fromBus = "B2", toBus = "B3", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      # nameplate q = 5.0 (deliberately non-zero and far from the solved
      # value), a wide-open Q band so the machine never clamps
      addProsumer!(net = pvnet, busName = "B2", type = "SYNCHRONOUSMACHINE", p = 10.0, q = 5.0, vm_pu = 1.0, qMin = -100.0, qMax = 100.0)
      addProsumer!(net = pvnet, busName = "B3", type = "ENERGYCONSUMER", p = 40.0, q = 15.0)
      @test validate!(net = pvnet)[1]
      _, pverg = runpf!(pvnet, 40, 1e-8, 0)
      @test pverg == 0
      pb2 = Sparlectra.geNetBusIdx(net = pvnet, busName = "B2")
      @test get(pvnet.qLimitEvents, pb2, nothing) === nothing
      pv_report = buildACPFlowReport(pvnet)
      pv_b2_row = only(r for r in pv_report.nodes if r.bus == pb2)
      @test !isapprox(pv_b2_row.q_gen_MVar, 5.0; atol = 1e-3)
      @test isapprox(pv_b2_row.q_gen_MVar, pvnet.nodeVec[pb2]._qƩGen; atol = 1e-9)
    end)() end

    # The Q-limit outcome of a solve on a working copy did not reach the net
    # (before 0.60.0): the bus types the solve set and the scheduled Q the
    # classic modes clamp on the unit stayed on the merged copy of a net
    # with a closed link and on the island copies of an island-wise run. The
    # clamped bus stayed PV on the net, and the report took the Q of a PV bus
    # from the node, which the classic modes leave at the base solve.
    # sp_case60 (one coupler, Neubach_110 to Neubach_110b) with a 4.0 .. 5.5
    # MVAr band at Birkkamp_110 and Westtal_110: Birkkamp_110 ends at Qmax
    # below its U_set 1.01 (U 0.99663 pu under the active set, 1.00312 pu
    # under both classic modes) and was PV on the net with 7.258 MVAr in the
    # report under the classic modes. Truth: the same net with the coupler
    # resolved by hand (link open, the branches of Neubach_110b on
    # Neubach_110), solved without a merge, and the generator Q of every bus
    # from its solved voltages. two_islands_prio island-wise under
    # classic_simultaneous with a -5 .. 5 MVAr band at Brunntal_110 and
    # Kesselau_110: both clamped (report 9.645 and -30.325 MVAr before), and
    # Fennwiese_110, the reference island 2 promotes for want of one, clamped
    # at its Qmax of 30 MVAr with the reference moved to Dornberg_110
    # (report 41.80 MVAr before).
    @testset "Q-limit outcome of a merged or island-wise solve reaches the net" begin (function ()
      scf = joinpath(dirname(@__DIR__), "data", "scf")
      # generator Q of every bus from the solved voltages of a net without
      # closed links: injection V conj(Y V) plus what the bus draws itself
      q_from_voltages = function (n)
        active = [b for b in eachindex(n.nodeVec) if !(b in n.isoNodes)]
        Y = Sparlectra.createYBUS(net = n)
        @test size(Y, 1) == length(active)
        V = [n.nodeVec[b]._vm_pu * cis(deg2rad(n.nodeVec[b]._va_deg)) for b in active]
        _, q_draw = Sparlectra._bus_consumption_MVA(n)
        s = V .* conj.(Y * V) .* n.baseMVA
        return Dict(b => imag(s[k]) + q_draw[b] for (k, b) in enumerate(active))
      end
      narrow = function (resolved::Bool)
        n = redirect_stdout(() -> Sparlectra.importSCF(joinpath(scf, "sp_case60.scf.json")), devnull)
        setQLimits!(net = n, qmin_MVar = 4.0, qmax_MVar = 5.5, busName = ["Birkkamp_110", "Westtal_110"])
        l = only(n.linkVec)
        member = l.toBus
        if resolved
          @test !any(ps -> getPosumerBusIndex(ps) == member, n.prosumpsVec) && !any(sh -> sh.busIdx == member, n.shuntVec)
          l.status = 0
          for br in n.branchVec
            br.fromBus == member && (br.fromBus = l.fromBus)
            br.toBus == member && (br.toBus = l.fromBus)
          end
          empty!(n.isoNodes)
          markIsolatedBuses!(net = n, log = false)
        end
        return n, member
      end
      for (mode, u_bk) in ((:active_set, 0.99663), (:classic_simultaneous, 1.00312), (:classic_one_at_a_time, 1.00312))
        closed, member = narrow(false)
        resolved, _ = narrow(true)
        _, erg_c = @test_logs min_level = Logging.Warn runpf!(closed, 30, 1e-8, 0; method = :rectangular, qlimit_enforcement_mode = mode)
        _, erg_r = @test_logs min_level = Logging.Warn runpf!(resolved, 30, 1e-8, 0; method = :rectangular, qlimit_enforcement_mode = mode)
        @test (mode, erg_c, erg_r) == (mode, 0, 0)
        bk = geNetBusIdx(net = closed, busName = "Birkkamp_110")
        @test (mode, getNodeType(closed.nodeVec[bk]), closed.qLimitEvents[bk]) == (mode, Sparlectra.PQ, :max)
        @test isapprox(closed.nodeVec[bk]._vm_pu, u_bk; atol = 1e-5)
        others = [b for b in eachindex(closed.nodeVec) if b != member]
        @test (mode, getNodeType.(closed.nodeVec[others])) == (mode, getNodeType.(resolved.nodeVec[others]))
        @test closed.qLimitEvents == resolved.qLimitEvents
        @test maximum(abs(closed.nodeVec[b]._vm_pu - resolved.nodeVec[b]._vm_pu) for b in others) < 1e-9
        q_true = q_from_voltages(resolved)
        @test isapprox(q_true[bk], 5.5; atol = 1e-6)
        rep_c = Dict(r.bus => r for r in buildACPFlowReport(closed).nodes)
        rep_r = Dict(r.bus => r for r in buildACPFlowReport(resolved).nodes)
        for b in sort(collect(keys(q_true)))
          @test (mode, b, isapprox(rep_c[b].q_gen_MVar, q_true[b]; atol = 1e-6), isapprox(rep_r[b].q_gen_MVar, q_true[b]; atol = 1e-6)) == (mode, b, true, true)
        end
        # the stability check and the Q-V check give the verdict of the
        # resolved net; Westtal_110 fails under the classic modes (PV at
        # 6.136 MVAr above its band, within the Q hysteresis of the scan)
        st_c = Sparlectra.qlimitStabilityReport(closed)
        st_r = Sparlectra.qlimitStabilityReport(resolved)
        @test [(r.bus, r.end_type, r.reached_limit, r.ok) for r in st_c] == [(r.bus, r.end_type, r.reached_limit, r.ok) for r in st_r]
        @test all(isapprox(a.q_MVar, b.q_MVar; atol = 1e-6) && isapprox(a.q_report_MVar, a.q_MVar; atol = 1e-6) for (a, b) in zip(st_c, st_r))
        @test (mode, [r.bus_name for r in st_c if !r.ok]) == (mode, mode === :active_set ? String[] : ["Westtal_110"])
        @test qvCharacteristicViolations(closed) == qvCharacteristicViolations(resolved) == NamedTuple[]
      end

      n = redirect_stdout(() -> Sparlectra.importSCF(joinpath(scf, "two_islands_prio.scf.json")), devnull)
      setQLimits!(net = n, qmin_MVar = -5.0, qmax_MVar = 5.0, busName = ["Brunntal_110", "Kesselau_110"])
      _, erg = @test_logs (:info, r"Classical Q-limit enforcement changed reference bus") min_level = Logging.Info redirect_stdout(() -> runpf!(n, 30, 1e-8, 0; method = :rectangular, islands_enabled = true, qlimit_enforcement_mode = :classic_simultaneous), devnull)
      @test erg == 0
      q_true = q_from_voltages(n)
      rep = Dict(r.bus => r for r in buildACPFlowReport(n).nodes)
      for (name, type, side, q) in (("Brunntal_110", Sparlectra.PQ, :max, 5.0), ("Kesselau_110", Sparlectra.PQ, :min, -5.0), ("Fennwiese_110", Sparlectra.PQ, :max, 30.0), ("Dornberg_110", Sparlectra.Slack, nothing, nothing))
        b = geNetBusIdx(net = n, busName = name)
        @test (name, getNodeType(n.nodeVec[b]), get(n.qLimitEvents, b, nothing)) == (name, type, side)
        q === nothing || @test isapprox(q_true[b], q; atol = 1e-6)
      end
      for b in sort(collect(keys(q_true)))
        @test (b, isapprox(rep[b].q_gen_MVar, q_true[b]; atol = 1e-6)) == (b, true)
      end
    end)() end

    @testset "AC rescue ladder and DC fallback" begin (function ()
      # NOTE: private inner name — an anonymous fixture assigning to a name
      # that also exists in the enclosing testset would rebind that local.
      mkrescuenet = function (; p_mw::Float64 = 300.0)
        rnet = Net(name = "rescue_fixture", baseMVA = 100.0)
        addBus!(net = rnet, busName = "S", vn_kV = 110.0)
        addBus!(net = rnet, busName = "L", vn_kV = 110.0)
        addPIModelACLine!(net = rnet, fromBus = "S", toBus = "L", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
        addProsumer!(net = rnet, busName = "S", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "S")
        addProsumer!(net = rnet, busName = "L", type = "ENERGYCONSUMER", p = p_mw, q = 50.0)
        return rnet
      end
      # A poisoned start state the plain solve cannot recover from within the
      # iteration budget; the ladder restarts from this exact state.
      poison! = rnet -> Sparlectra.setVmVa!(node = rnet.nodeVec[2], vm_pu = 0.01, va_deg = 179.0)

      # Without rescue the poisoned run fails.
      plain = mkrescuenet()
      poison!(plain)
      cfg_plain = Sparlectra.PowerFlowConfig(max_iter = 5, tol = 1e-8, rescue = false, start_mode = Sparlectra.StartModeConfig(flatstart = false))
      @test last(runpf!(plain, cfg_plain)) == 1

      # With rescue a ladder strategy converges and is recorded.
      rescued = mkrescuenet()
      poison!(rescued)
      cfg_rescue = Sparlectra.PowerFlowConfig(max_iter = 5, tol = 1e-8, rescue = true, start_mode = Sparlectra.StartModeConfig(flatstart = false))
      profile = Dict{Symbol,Any}()
      rite, rerg = runpf!(rescued, cfg_rescue; performance_profile = profile)
      @test rerg == 0
      @test rite > 0
      @test get(profile, :ac_rescue_strategy, :missing) in (:alternate_start, :autodamp, :dc_seed, :settled_qlimits)
      # the solver status says the solution came from the rescue, with the
      # failed first solve's start and iterations
      rst = Sparlectra.rectangular_pf_status(rescued)
      @test rst.rescue_used === true
      @test rst.rescue_strategy === profile[:ac_rescue_strategy]
      @test rst.rescue_first_flatstart === false
      @test rst.rescue_first_iterations == 5

      # The ladder must offer the globalization package for Q-limit-noisy
      # systems: merit line search, damping, and Q-limit switching held back
      # until the reactive requests settle. Asserted on
      # the variant list so a large-case regression cannot silently drop it.
      variants = Sparlectra._rescue_config_variants(Sparlectra.PowerFlowConfig())
      settled = only(cfg for (name, cfg) in variants if name === :settled_qlimits)
      @test settled.merit.enabled === true
      @test settled.autodamp === true
      @test settled.qlimits.start_mode === :auto
      # honest: the limits themselves stay enforced, only the timing changes
      @test settled.qlimits.ignore_q_limits === false

      # A genuinely infeasible case (load far beyond the transfer limit):
      # rescue exhausts the ladder, the DC fallback leaves a usable DC state,
      # and the AC status honestly stays non-converged.
      infeasible = mkrescuenet(p_mw = 5000.0)
      cfg_fb = Sparlectra.PowerFlowConfig(max_iter = 10, tol = 1e-8, rescue = true, dc = Sparlectra.DcPowerFlowConfig(fallback = true))
      profile_fb = Dict{Symbol,Any}()
      ferg = last(runpf!(infeasible, cfg_fb; performance_profile = profile_fb))
      @test ferg == 1
      @test get(profile_fb, :ac_rescue_strategy, :missing) == :none
      @test Sparlectra.rectangular_pf_status(infeasible).rescue_strategy === :none
      @test get(profile_fb, :dc_fallback_applied, false) === true
      @test infeasible.nodeVec[2]._vm_pu == 1.0
      @test infeasible.nodeVec[2]._va_deg < 0.0

      # The rescue ladder is ON by default (a user should get a result, not a
      # divergence message); the DC fallback stays opt-in because it swaps the
      # model. Explicitly disabling rescue reproduces the raw failure.
      @test Sparlectra.PowerFlowConfig().rescue === true
      @test Sparlectra.PowerFlowConfig().dc.fallback === false
      off = mkrescuenet()
      poison!(off)
      cfg_off = Sparlectra.PowerFlowConfig(max_iter = 5, tol = 1e-8, rescue = false, start_mode = Sparlectra.StartModeConfig(flatstart = false))
      profile_off = Dict{Symbol,Any}()
      @test last(runpf!(off, cfg_off; performance_profile = profile_off)) == 1
      @test !haskey(profile_off, :ac_rescue_strategy)
      @test !hasproperty(Sparlectra.rectangular_pf_status(off), :rescue_used)
      @test !haskey(profile_off, :dc_fallback_applied)
    end)() end

      @testset "Active-set Q-limits switch on the converged iterate (post-convergence round)" begin (function ()
        # Regression: a converged iterate is always "ready" for the Q-limit
        # check. Previously the start gating could outlast convergence (NR
        # faster than qlimit_start_iter, or :auto never settling), the check
        # was skipped on the converged iterate, and the run exited with its
        # violations unevaluated — rejected as remaining_pv_q_limit_violations
        # instead of being enforced. The gate here sits far beyond the ~5
        # iterations this net needs, so ONLY the post-convergence round can
        # perform the (legitimate) PV→PQ switch.
        net = createTest3BusNet(cooldown = 2, hyst_pu = 0.01, qlim_min = -15.0, qlim_max = 15.0)
        _, erg = runpf!(net, 40, 1e-8, 0; method = :rectangular, qlimit_start_iter = 50)
        @test erg == 0
        st = Sparlectra.rectangular_pf_status(net)
        @test st.final_converged == true
        @test st.reason == :none
        @test length(net.qLimitLog) == 1
        # same net, same limits, default gating — the reference behavior the
        # post-convergence round must reproduce
        net_ref = createTest3BusNet(cooldown = 2, hyst_pu = 0.01, qlim_min = -15.0, qlim_max = 15.0)
        _, erg_ref = runpf!(net_ref, 40, 1e-8, 0; method = :rectangular)
        @test erg_ref == 0
        @test length(net_ref.qLimitLog) == 1
      end)() end

      @testset "Classical Q-limit violation accounting excludes fixed generator-type units" begin (function ()
        # Regression: _matpower_q_limit_violations shared the FULL bus reactive
        # requirement equally among all non-clamped generator-type units. A
        # fixed injection at a PV bus (never-regulating: a Stage-0 HVDC
        # converter, a kept boundary equivalent, a PQ-degraded SVC) then
        # "violated" its own zero-width limits, was clamped, and from the next
        # round on its reactive power was loaded onto the remaining regulator —
        # false violations and wrong clamps on any network with a fixed and a
        # regulating unit on one bus (nothing CGMES-specific).
        #
        # Here the bus needs ≈50 MVAr: 40 come from the fixed unit, ≈10 from
        # the regulator whose ±30 MVAr limits are comfortable. Correct
        # accounting finds NO violation; the old equal share attributed 25/25
        # and clamped both units in turn.
        net = Net(name = "qlimit_fixed_unit", baseMVA = 100.0)
        addBus!(net = net, busName = "B1", vn_kV = 110.0)
        addBus!(net = net, busName = "B2", vn_kV = 110.0)
        addBus!(net = net, busName = "B3", vn_kV = 110.0)
        addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = net, fromBus = "B1", toBus = "B3", r_pu = 0.02, x_pu = 0.12, b_pu = 0.0, status = 1)
        addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
        addProsumer!(net = net, busName = "B2", type = "GENERATOR", p = 20.0, vm_pu = 1.0, qMin = -30.0, qMax = 30.0)
        # the fixed unit: generator-class, never regulating, zero-width Q range
        addProsumer!(net = net, busName = "B2", type = "EXTERNALNETWORKINJECTION", p = 0.0, q = 40.0, qMin = 40.0, qMax = 40.0, isRegulated = false)
        addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 60.0, q = 52.0)
        ok, _ = validate!(net = net)
        @test ok

        _, erg = runpf!(net, 30, 1e-8, 0; method = :rectangular, qlimit_enforcement_mode = :classic_one_at_a_time, qlimit_max_outer = 10)
        @test erg == 0
        st = Sparlectra.rectangular_pf_status(net)
        @test st.final_outcome == :converged
        # neither the fixed unit nor the regulator may have been clamped
        @test isempty(st.matpower_outer_loop)
        @test isempty(net.qLimitLog)
        # the PV bus must still hold its setpoint (nothing degraded it to PQ)
        @test Sparlectra.getNodeType(net.nodeVec[Sparlectra.geNetBusIdx(net = net, busName = "B2")]) == Sparlectra.PV
      end)() end

      @testset "Classical Q-limit outer loop does not overwrite each round's diagnostic artifacts" begin (function ()
        # Regression: the classical outer loop's inner runpf_rectangular! calls all
        # shared performance_profile[:output_dir] with no per-round prefix, so
        # trust_region.log (like merit_linesearch.log/current_iteration_start.log)
        # from an earlier outer round was silently overwritten by the next round's
        # solve -- same bug class as the AC-island diagnostic-overwrite fix.
        mktempdir() do tmpdir
          perf = Dict{Symbol,Any}(:output_dir => tmpdir)
          net = createTest3BusNet(cooldown = 2, hyst_pu = 0.01, qlim_min = -15.0, qlim_max = 15.0)
          _, erg = runpf!(net, 30, 1e-8, 0; method = :rectangular, autodamp = false, trust_region_enabled = true, qlimit_enforcement_mode = :classic_one_at_a_time, qlimit_max_outer = 5, performance_profile = perf)
          @test erg == 0
          @test length(net.qLimitLog) == 1
          files = readdir(tmpdir)
          @test "qlimit_outer_0_trust_region.log" in files
          @test "qlimit_outer_1_trust_region.log" in files
          @test !("trust_region.log" in files)
          @test !haskey(perf, :diagnostic_artifact_prefix)
        end
      end)() end

  return true
end

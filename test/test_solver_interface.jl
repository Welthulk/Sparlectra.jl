# Copyright 2023–2026 Udo Schmitz
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

# file: test/test_solver_interface.jl
# purpose: tests the external solver interface (PFModel/PFSolution) and the
#          rectangular solver's control paths: Q-limit modes, merit line
#          search, trust region, AC islands, and the rescue ladder

using Test
using Sparlectra



"""
    test_external_solver_interface() -> Bool

Minimal unit test for the external solver interface:
- PFModel construction
- canonical mismatch evaluation
- solution write-back into Net

Returns `true` if all checks pass.
"""
function test_external_solver_interface()::Bool

  @testset "External solver interface (PFModel/PFSolution)" begin (function ()

    # Existing 3-bus test network from testgrid.jl
    net = createTest3BusNet()

    # Build canonical PF model
    model = buildPfModel(
      net;
      opt_sparse = true,
      flatstart = true,
      include_limits = false,
      verbose = 0,
    )

    # --- Structural checks ---------------------------------------------------
    @test size(model.Ybus, 1) == length(model.busIdx_net)
    @test size(model.Ybus, 2) == length(model.busIdx_net)

    @test length(model.busIdx_net) == 3
    @test model.busIdx_net == [1, 2, 3]
    @test model.busType == [:PQ, :PV, :Slack]
    @test model.slack_idx == 3

    # --- Initial state -------------------------------------------------------
    @test length(model.V0) == 3
    @test all(isfinite.(abs.(model.V0)))

    # --- Canonical mismatch --------------------------------------------------
    r0 = mismatchInf(model, model.V0)
    @test isfinite(r0)
    @test r0 ≥ 0.0

    # --- Apply solution back into Net ---------------------------------------
    sol = PFSolution(
      V = copy(model.V0),
      converged = true,
      iters = 0,
      residual_inf = r0,
    )

    applyPfSolution!(
      net,
      model,
      sol;
      write_pq_results = false,
      verbose = 0,
    )

    for (k, busIdx) in enumerate(model.busIdx_net)
      node = net.nodeVec[busIdx]
      @test isapprox(node._vm_pu, abs(model.V0[k]); atol = 1e-12)
      @test isapprox(node._va_deg, rad2deg(angle(model.V0[k])); atol = 1e-10)
    end

  end)() end

  return true
end

function run_solver_interface_tests()
  # Covers solver-interface integration and option behavior:
  # external model API, Q-limit reporting/autocorrection, PV->PQ locking, and final-limit reporting.
  @testset "Solver interface" begin (function ()
    @testset "tolerance physical equivalent" begin (function ()
      # the run log, the diagnostics and the docs show tol * baseMVA in a
      # readable unit; the default 1e-8 pu at the 100 MVA base is 1 W
      @test Sparlectra.format_tolerance_physical(1.0e-8, 100.0) == "1.0e-8 pu, equals 1 W at 100.0 MVA base"
      @test Sparlectra.format_tolerance_physical(1.0e-5, 100.0) == "1.0e-5 pu, equals 1 kW at 100.0 MVA base"
      @test Sparlectra.format_tolerance_physical(0.01, 100.0) == "0.01 pu, equals 1 MW at 100.0 MVA base"
      @test occursin("W at 50.0 MVA base", Sparlectra.format_tolerance_physical(1.0e-8, 50.0))
    end)() end

    @testset "power_flow.tol_MW wins over tol at the case base (part B)" begin (function ()
      # absent: nothing changes, the per-unit default stands
      base_cfg = Sparlectra.SparlectraConfig()
      @test base_cfg.powerflow.tol_MW === nothing
      @test base_cfg.powerflow.tol == 1.0e-8
      net = Sparlectra.importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json")))
      @test Sparlectra._resolve_tolerance_for_net(base_cfg.powerflow, net) === base_cfg.powerflow
      # set: converted with the NETWORK base, and it outranks tol
      cfg = Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("tol" => 1.0e-5, "tol_MW" => 0.002)))
      @test cfg.powerflow.tol_MW == 0.002
      eff = Sparlectra._resolve_tolerance_for_net(cfg.powerflow, net)
      @test eff.tol == 0.002 / net.baseMVA
      @test eff.tol != cfg.powerflow.tol
      # validation: zero and negative are errors NAMING the key, an
      # implausibly coarse bound is a warning and stays usable
      for bad in (0.0, -1.0)
        err = try
          Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("tol_MW" => bad)))
          nothing
        catch e
          e
        end
        @test err isa ArgumentError
        @test occursin("power_flow.tol_MW", sprint(showerror, err))
      end
      coarse = Test.@test_logs (:warn, r"very coarse convergence bound") match_mode = :any Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("tol_MW" => 5.0)))
      @test coarse.powerflow.tol_MW == 5.0
      # the key travels through the override surface like any other
      @test "power_flow.tol_MW" in Sparlectra.GUI_EDITABLE_CONFIG_KEYS
    end)() end

    test_external_solver_interface()
    @testset "Rectangular PF status lives on the Net" begin (function ()
      # the status is a Net field, so distinct Nets are independent and
      # deepcopy carries the status with the copy (thread safety)
      net = createTest3BusNet()
      @test Sparlectra.rectangular_pf_status(net) === nothing
      status = (status = :test_status,)
      @test Sparlectra._set_rectangular_pf_status!(net, status) === status
      @test Sparlectra.rectangular_pf_status(net) === status

      other = createTest3BusNet()
      @test Sparlectra.rectangular_pf_status(other) === nothing

      copied = deepcopy(net)
      @test Sparlectra.rectangular_pf_status(copied).status === :test_status
      Sparlectra._set_rectangular_pf_status!(copied, (status = :copy_only,))
      @test Sparlectra.rectangular_pf_status(net) === status
    end)() end
    @testset "Flat-start voltage setpoints" begin (function ()
      net = createTest3BusNet()
      net.nodeVec[1]._vm_pu = 0.94
      net.nodeVec[1]._va_deg = -7.0
      net.nodeVec[2]._vm_pu = 1.04
      net.nodeVec[2]._va_deg = -3.0
      net.nodeVec[3]._vm_pu = 1.06
      net.nodeVec[3]._va_deg = 2.0

      Vflat, slack_idx = Sparlectra.initialVrect(net; flatstart = true)
      @test slack_idx == 3
      # every non-slack bus starts at the slack angle (2 degrees)
      @test abs(Vflat[1]) == 1.0
      @test isapprox(angle(Vflat[1]), deg2rad(2.0); atol = 1e-12)
      @test abs(Vflat[2]) == 1.04
      @test isapprox(angle(Vflat[2]), deg2rad(2.0); atol = 1e-12)
      @test abs(Vflat[3]) == 1.06
      @test isapprox(rad2deg(angle(Vflat[3])), 2.0; atol = 1e-12)

      model = Sparlectra.buildPfModel(net; flatstart = true, include_limits = false, start_projection = false)
      @test model.busIdx_net == [1, 2, 3]
      @test abs(model.V0[1]) == 1.0
      @test angle(model.V0[1]) == 0.0
      @test abs(model.V0[2]) == 1.04
      @test angle(model.V0[2]) == 0.0
      @test abs(model.V0[3]) == 1.06
      @test isapprox(rad2deg(angle(model.V0[3])), 2.0; atol = 1e-12)

      Vseeded, _ = Sparlectra.initialVrect(net; flatstart = false)
      @test isapprox(abs(Vseeded[1]), 0.94; atol = 1e-12)
      @test isapprox(rad2deg(angle(Vseeded[1])), -7.0; atol = 1e-12)
      @test isapprox(abs(Vseeded[2]), 1.04; atol = 1e-12)
      @test isapprox(rad2deg(angle(Vseeded[2])), -3.0; atol = 1e-12)
      @test isapprox(abs(Vseeded[3]), 1.06; atol = 1e-12)
      @test isapprox(rad2deg(angle(Vseeded[3])), 2.0; atol = 1e-12)

      iso_net = Net(name = "iso_flatstart", baseMVA = 100.0)
      for busName in ("B1", "B2", "B3", "B4")
        addBus!(net = iso_net, busName = busName, vn_kV = 110.0)
      end
      add_test_line_pu_c!(net = iso_net, fromBus = "B1", toBus = "B3", length_km = 10.0, r_pu = 0.0, x_pu = 0.4, c_nf_per_km = 9.55, tanδ = 0.0)
      add_test_line_pu_c!(net = iso_net, fromBus = "B3", toBus = "B4", length_km = 10.0, r_pu = 0.0, x_pu = 0.4, c_nf_per_km = 9.55, tanδ = 0.0)
      add_test_line_pu_c!(net = iso_net, fromBus = "B4", toBus = "B1", length_km = 10.0, r_pu = 0.0, x_pu = 0.4, c_nf_per_km = 9.55, tanδ = 0.0)
      addProsumer!(net = iso_net, busName = "B4", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.06, va_deg = 5.0, referencePri = "B4")
      addProsumer!(net = iso_net, busName = "B3", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 0.0, vm_pu = 1.04)
      addProsumer!(net = iso_net, busName = "B1", type = "ENERGYCONSUMER", p = 30.0, q = 10.0)
      setNodeVoltage!(net = iso_net, busName = "B4", vm_pu = 1.06, va_deg = 5.0)
      setPVBusVset!(net = iso_net, busName = "B3", vm_pu = 1.04)
      markIsolatedBuses!(net = iso_net, log = false)

      iso_model = Sparlectra.buildPfModel(iso_net; flatstart = true, include_limits = false, start_projection = false)
      @test iso_net.isoNodes == [2]
      @test iso_model.busIdx_net == [1, 3, 4]
      @test iso_model.busType == [:PQ, :PV, :Slack]
      @test abs.(iso_model.V0) == [1.0, 1.04, 1.06]
      @test all(isapprox.(rad2deg.(angle.(iso_model.V0)), [0.0, 0.0, 5.0]; atol = 1e-12))
      @test iso_model.Vset == [1.0, 1.04, 1.06]
    end)() end
    # Checks the Q-limit switch log plus sign-validation/autocorrect behavior.
    @testset "Q-limit reporting and validation options" begin (function ()
      net = createTest3BusNet()
      net.nodeVec[1]._nodeType = Sparlectra.PV
      net.nodeVec[2]._nodeType = Sparlectra.PV
      net.nodeVec[3]._nodeType = Sparlectra.PV
      empty!(net.qmin_pu)
      append!(net.qmin_pu, [-0.1, 0.05, -0.2])
      empty!(net.qmax_pu)
      append!(net.qmax_pu, [0.2, 0.1, 0.3])

      switch_io = IOBuffer()
      active_set_q_limits!(
        net,
        10,
        length(net.nodeVec);
        get_qreq_pu = bus -> bus == 2 ? 0.2 : 0.0,
        is_pv = bus -> bus == 2,
        make_pq! = (bus, qclamp, side) -> nothing,
        make_pv! = bus -> nothing,
        qmin_pu = fill(-Inf, length(net.nodeVec)),
        qmax_pu = [Inf, 0.1, Inf],
        pv_orig_mask = trues(length(net.nodeVec)),
        allow_reenable = false,
        q_hyst_pu = 0.0,
        cooldown_iters = 0,
        verbose = 1,
        io = switch_io,
      )
      switch_log = String(take!(switch_io))
      @test occursin("Q=0.200000 pu", switch_log)
      @test occursin("Qmax=0.100000 pu", switch_log)
      @test occursin("MVAr", switch_log)

      res = validate_q_limit_signs!(net.qmin_pu, net.qmax_pu; io = IOBuffer(), autocorrect = true, warn = false)
      @test res.flagged >= 1
      @test res.corrected >= 1
      @test net.qmin_pu[2] <= 0.0
      @test net.qmax_pu[2] >= 0.0
    end)() end

    # power_flow.qlimits.lock_pv_to_pq_buses (#458): a listed PV bus runs as
    # PQ from the start, at its scheduled Q (33.2 MVAr at STATION1) clamped
    # into its limits, logged at iteration 0, never released back to PV. The
    # list used to only exclude the bus from PV->PQ switching: with limits of
    # +-1 MVAr the bus stayed PV at 32.9 MVAr and the run failed. Guard off,
    # so nothing but the list reclassifies the bus; both enforcement families.
    @testset "PV->PQ lock option" begin (function ()
      vset = 1.027273
      for mode in (:active_set, :classic_simultaneous)
        solve(lock, limited; qlimits_enabled = true) = begin
          n = createTest3BusNet()
          limited && setQLimits!(net = n, qmin_MVar = -1.0, qmax_MVar = 1.0, busName = "STATION1")
          _, erg = runpf!(n, 20, 1e-6, 0; method = :rectangular, qlimit_guard = false, qlimit_enforcement_mode = mode, qlimits_enabled = qlimits_enabled, lock_pv_to_pq_buses = lock)
          (n, erg)
        end
        # limits of +-1 MVAr: unlocked the active set clamps the machine
        # during the iteration, locked it starts as PQ at the clamped value
        net_unlocked, erg_unlocked = solve(Int[], true)
        @test (mode, erg_unlocked) == (mode, 0)
        @test getNodeType(net_unlocked.nodeVec[2]) == Sparlectra.PQ
        @test first(net_unlocked.qLimitLog).iter > 0
        net_locked, erg_locked = solve([2], true)
        @test (mode, erg_locked) == (mode, 0)
        @test getNodeType(net_locked.nodeVec[2]) == Sparlectra.PQ
        @test [(e.iter, e.bus, e.side) for e in net_locked.qLimitLog] == [(0, 2, :max)]
        @test isapprox(net_locked.nodeVec[2]._vm_pu, net_unlocked.nodeVec[2]._vm_pu; atol = 1e-6)
        # no limits at all: unlocked the machine holds its setpoint as PV,
        # locked it runs as PQ at its scheduled 33.2 MVAr and the voltage
        # floats off the setpoint
        free_unlocked, erg_fu = solve(Int[], false)
        @test erg_fu == 0
        @test getNodeType(free_unlocked.nodeVec[2]) == Sparlectra.PV
        @test isapprox(free_unlocked.nodeVec[2]._vm_pu, vset; atol = 1e-6)
        free_locked, erg_fl = solve([2], false)
        @test erg_fl == 0
        @test getNodeType(free_locked.nodeVec[2]) == Sparlectra.PQ
        @test isapprox(free_locked.nodeVec[2]._qƩGen, 33.2; atol = 1e-9)
        @test abs(free_locked.nodeVec[2]._vm_pu - vset) > 1e-4
        # the list is a Q-limit setting: with Q limits ignored it does nothing
        ignored, erg_ig = solve([2], false; qlimits_enabled = false)
        @test erg_ig == 0
        @test getNodeType(ignored.nodeVec[2]) == Sparlectra.PV
      end
      # the configuration path reaches the same conversion, and the final
      # check accepts the run (it used to end final_converged = false)
      quiet = OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)
      cfg_net = createTest3BusNet()
      setQLimits!(net = cfg_net, qmin_MVar = -1.0, qmax_MVar = 1.0, busName = "STATION1")
      r = redirect_stdout(devnull) do
        run_sparlectra(net = cfg_net, config = SparlectraConfig(powerflow = PowerFlowConfig(qlimits = QLimitConfig(guard = false, lock_pv_to_pq_buses = [2])), output = quiet))
      end
      @test r.final_converged
      @test getNodeType(cfg_net.nodeVec[2]) == Sparlectra.PQ
      @test isapprox(cfg_net.nodeVec[2]._qƩGen, 1.0; atol = 1e-9)
      # Q-limit stability check (Zhao et al. 2008, (6)/(7)): the lock holds
      # STATION1 as PQ at its 50 MVAr Qmax (schedule 60 MVAr) although the
      # machine needs 31.5 MVAr to hold its setpoint, so the voltage rises
      # above it: a machine at Qmax above U_set violates (7), the one
      # violation of the run, named with its lock. The +-1 MVAr run above
      # sits at Qmax below U_set and passes. Unlocked the machine holds U_set
      # as PV inside its band: no violation.
      @test r.diagnostics.qlimit_stability_violations == 0
      stability_run(lock) = begin
        n = createTest3BusNet()
        setQLimits!(net = n, qmin_MVar = -100.0, qmax_MVar = 50.0, busName = "STATION1")
        for ps in n.prosumpsVec
          Sparlectra.getPosumerBusIndex(ps) == 2 && Sparlectra.isGenerator(ps) && (ps.qVal = 60.0)
        end
        redirect_stdout(devnull) do
          run_sparlectra(net = n, config = SparlectraConfig(powerflow = PowerFlowConfig(qlimits = QLimitConfig(guard = false, lock_pv_to_pq_buses = lock)), output = quiet))
        end
      end
      r_lock = @test_logs (:warn, r"Q-limit stability check: 1 voltage-regulated bus\(es\).*STATION1 \(PQ at Qmax.*locked by lock_pv_to_pq_buses") match_mode = :any stability_run([2])
      @test r_lock.final_converged
      @test r_lock.diagnostics.qlimit_stability_violations == 1
      row = only(r_lock.diagnostics.qlimit_stability_report)
      @test (row.bus, row.end_type, row.reached_limit, row.check_6, row.check_7, row.ok, row.locked_by) == (2, :PQ, :max, true, false, false, :lock_pv_to_pq_buses)
      @test row.u_pu > row.uset_pu + 1e-3
      @test isapprox(row.q_MVar, 50.0; atol = 1e-6)
      @test qlimit_stability_line(r_lock.diagnostics.qlimit_stability_status, r_lock.diagnostics.qlimit_stability_report) == "Q-limit stability violations: 1 (bus STATION1)"
      r_free = @test_logs min_level = Logging.Warn stability_run(Int[])
      @test r_free.final_converged
      @test r_free.diagnostics.qlimit_stability_violations == 0
      @test only(r_free.diagnostics.qlimit_stability_report).end_type === :PV
      # auto mode: the check takes the Q-limit configuration of the attempt
      # that produced the result, not the user configuration. Without Q
      # limits at any generator the auto gate switches the handling off, so
      # the run is not checked (the user configuration said on and gave 0
      # violations on a run that never judged a limit)
      n_auto = createTest3BusNet()
      r_auto = @test_logs min_level = Logging.Warn redirect_stdout(() -> run_sparlectra(net = n_auto, config = SparlectraConfig(powerflow = PowerFlowConfig(mode = :auto), output = quiet)), devnull)
      @test r_auto.final_converged
      @test Sparlectra.auto_pf_record(n_auto).final_qlimits.ignore_q_limits
      @test (r_auto.diagnostics.qlimit_stability_status, r_auto.diagnostics.qlimit_stability_violations) == (:qlimits_disabled, 0)
      # a regulated machine (band -1 .. 1 MVAr) on a bus joined to the slack
      # bus by a closed link, added before the units as a CGMES import adds
      # its closed switches. The slack balances the cluster, so the 10.8
      # MVAr the machine carries is a share of the cluster Q that no limit
      # binds: the bus is out of the check like the slack itself. It was
      # reported as a PV machine outside its band on one island and left out
      # with a second island (INSEL_A/INSEL_B, its own reference) in the net.
      linked_slack_run(island) = begin
        n = Net(name = "link3", baseMVA = 100.0)
        for b in ("ASTADT", "STATION1", "VERBUND", "VERBUND_B")
          addBus!(net = n, busName = b, vn_kV = 110.0)
        end
        island && for b in ("INSEL_A", "INSEL_B")
          addBus!(net = n, busName = b, vn_kV = 110.0)
        end
        for (f, t) in (("ASTADT", "STATION1"), ("ASTADT", "VERBUND"), ("VERBUND", "STATION1"))
          add_test_line_pu_c!(net = n, fromBus = f, toBus = t, length_km = 25.0, r_pu = 0.0, x_pu = 0.4, c_nf_per_km = 9.55, tanδ = 0.0)
        end
        island && add_test_line_pu_c!(net = n, fromBus = "INSEL_A", toBus = "INSEL_B", length_km = 25.0, r_pu = 0.0, x_pu = 0.4, c_nf_per_km = 9.55, tanδ = 0.0)
        addLink!(net = n, fromBus = "VERBUND", toBus = "VERBUND_B", status = 1)
        addProsumer!(net = n, busName = "VERBUND", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.018182, va_deg = 0.0, referencePri = "VERBUND")
        addProsumer!(net = n, busName = "VERBUND_B", type = "SYNCHRONOUSMACHINE", p = 10.0, q = 0.0, vm_pu = 1.018182, qMin = -1.0, qMax = 1.0)
        addProsumer!(net = n, busName = "STATION1", type = "SYNCHRONOUSMACHINE", p = 70.0, q = 33.2, vm_pu = 1.027273)
        addProsumer!(net = n, busName = "ASTADT", type = "ENERGYCONSUMER", p = 100.0, q = 30.0)
        if island
          addProsumer!(net = n, busName = "INSEL_A", type = "SYNCHRONOUSMACHINE", p = 10.0, q = 0.0, vm_pu = 1.01, qMin = -30.0, qMax = 30.0)
          addProsumer!(net = n, busName = "INSEL_B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
        end
        redirect_stdout(() -> run_sparlectra(net = n, config = SparlectraConfig(output = quiet)), devnull)
      end
      for island in (false, true)
        r_link = @test_logs min_level = Logging.Warn linked_slack_run(island)
        @test (island, r_link.final_converged, r_link.diagnostics.qlimit_stability_violations) == (island, true, 0)
        @test (island, [row.bus_name for row in r_link.diagnostics.qlimit_stability_report]) == (island, ["STATION1"])
      end
      # a Q(U) generator next to the regulating machine of STATION1 (band
      # -10 .. 10 MVAr, scheduled at 0 MVAr, needs 31.5 to hold U_set). The
      # bus band counts the Q(U) unit at its Q at the held voltage, -1.091
      # MVAr (characteristic 0.02 pu at 0.95, -0.02 pu at 1.05), so it is
      # -11.091 .. 8.909 MVAr. The result report keeps the schedule of the
      # regulating machine at such a bus, so the check reads Q from the
      # voltages. Clamped at Qmax below U_set the bus passes at 8.909 MVAr;
      # the report value, 1.179 MVAr (schedule plus the Q(U) unit at the
      # solved 0.9705 pu), made it a PQ machine inside its band off U_set.
      qu_run(qmax) = begin
        n = createTest3BusNet()
        setQLimits!(net = n, qmin_MVar = -qmax, qmax_MVar = qmax, busName = "STATION1")
        for ps in n.prosumpsVec
          Sparlectra.getPosumerBusIndex(ps) == 2 && Sparlectra.isGenerator(ps) && (ps.qVal = 0.0)
        end
        addProsumer!(net = n, busName = "STATION1", type = "GENERATOR", p = 0.0, q = 0.0, qu_controller = QUController(VoltageCharacteristic([(0.95, 0.02), (1.05, -0.02)]); qmin_pu = -0.05, qmax_pu = 0.05))
        res = redirect_stdout(devnull) do
          run_sparlectra(net = n, config = SparlectraConfig(powerflow = PowerFlowConfig(qlimits = QLimitConfig(guard = false)), output = quiet))
        end
        (n, res)
      end
      n_qu, r_qu = @test_logs min_level = Logging.Warn qu_run(10.0)
      @test (r_qu.final_converged, r_qu.diagnostics.qlimit_stability_violations) == (true, 0)
      row = only(r_qu.diagnostics.qlimit_stability_report)
      @test (row.end_type, row.reached_limit, row.ok, row.q_source) == (:PQ, :max, true, :voltages)
      @test isapprox(row.qmax_MVar, 8.90908; atol = 1e-5)
      @test isapprox(row.q_MVar, row.qmax_MVar; atol = 1e-6)
      @test isapprox(row.q_report_MVar, 1.17889; atol = 1e-5)
      @test row.u_pu < row.uset_pu - 1e-3
      # free as PV in a wide band (+-100 MVAr) the machine holds U_set at the
      # Q the solver writes to the node; with the band narrowed to +-10 MVAr
      # on the solved state that Q is beyond Qmax and the bus fails. The
      # report value (-1.091 MVAr, the Q(U) unit alone) passed it.
      n_pv, r_pv = qu_run(100.0)
      @test (r_pv.final_converged, only(r_pv.diagnostics.qlimit_stability_report).end_type) == (true, :PV)
      setQLimits!(net = n_pv, qmin_MVar = -10.0, qmax_MVar = 10.0, busName = "STATION1")
      row = only(Sparlectra.qlimitStabilityReport(n_pv))
      @test (row.end_type, row.in_band, row.ok, row.q_source) == (:PV, false, false, :voltages)
      @test isapprox(row.q_MVar, n_pv.nodeVec[2]._qƩGen; atol = 1e-6)
      @test isapprox(row.q_MVar, 31.4525; atol = 1e-3)
      @test isapprox(row.q_report_MVar, -1.09092; atol = 1e-5)
    end)() end

    @testset "Q-limit guard: a narrow band locks only at a reached limit" begin (function ()
      # the machine at STATION1 needs q_free to hold its setpoint; a band of
      # 1 MVAr (0.01 pu, below the 0.02 pu guard threshold) around q_free
      # holds it, a band of 1 MVAr below q_free does not
      solve(band, mode) = begin
        n = createTest3BusNet()
        band === nothing || setQLimits!(net = n, qmin_MVar = band[1], qmax_MVar = band[2], busName = "STATION1")
        _, erg = runpf!(n, 20, 1e-6, 0; method = :rectangular, qlimit_enforcement_mode = :active_set, qlimit_guard = true, qlimit_guard_narrow_range_mode = mode)
        (n, erg)
      end
      free, erg_free = solve(nothing, :lock_pq)
      @test erg_free == 0
      q_free = free.nodeVec[2]._qƩGen
      vset = free.nodeVec[2]._vm_pu
      # inside the band: lock_pq leaves the machine PV at its setpoint, no
      # Q-limit event; prefer_pq runs it as PQ at the band middle from the
      # start (iteration 0)
      inside, erg_in = solve((q_free - 0.5, q_free + 0.5), :lock_pq)
      @test erg_in == 0
      @test getNodeType(inside.nodeVec[2]) == Sparlectra.PV
      @test isempty(inside.qLimitLog)
      @test isapprox(inside.nodeVec[2]._vm_pu, vset; atol = 1e-6)
      middle, erg_mid = solve((q_free - 0.5, q_free + 0.5), :prefer_pq)
      @test erg_mid == 0
      @test getNodeType(middle.nodeVec[2]) == Sparlectra.PQ
      @test [(e.iter, e.bus) for e in middle.qLimitLog] == [(0, 2)]
      # above the band: lock_pq clamps the machine at Qmax when the active
      # set reaches it (an iteration after 0) and keeps it there
      above, erg_ab = solve((q_free - 1.5, q_free - 0.5), :lock_pq)
      @test erg_ab == 0
      @test getNodeType(above.nodeVec[2]) == Sparlectra.PQ
      @test isapprox(above.nodeVec[2]._qƩGen, q_free - 0.5; atol = 1e-9)
      @test all(e -> e.iter > 0 && e.side === :max, above.qLimitLog)

      # Q-limit stability check on sp_case60 with a 4.0 .. 5.5 MVAr band
      # (0.015 pu, below the guard threshold) at the regulated machines of
      # Birkkamp_110 and Westtal_110 (generators 159 and 160 of the SCF
      # file). Birkkamp_110 needs 7.5 MVAr and ends at Qmax with U 1.0026
      # below U_set 1.01, (6) and (7) hold; Westtal_110 holds U_set at
      # 5.33 MVAr inside its band. With the guard the lock at the first
      # reached limit is named, without it the active set clamps the bus on
      # its own: 0 violations either way. The net has a closed coupler, so Q
      # comes from the solved voltages, and the result report gives the same
      # 5.500 MVAr under every mode (before 0.60.0 it gave 7.389 MVAr under
      # classic_simultaneous: the merged solve left the bus PV on the net).
      c60 = joinpath(dirname(@__DIR__), "data", "scf", "sp_case60.scf.json")
      narrow_run(mode, guard; shunt_q = 0.0) = begin
        n = redirect_stdout(() -> Sparlectra.importSCF(c60), devnull)
        setQLimits!(net = n, qmin_MVar = 4.0, qmax_MVar = 5.5, busName = ["Birkkamp_110", "Westtal_110"])
        shunt_q == 0.0 || addShunt!(net = n, busName = "Westtal_110", pShunt = 0.0, qShunt = shunt_q)
        cfg = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, c60, Dict{String,Any}("power_flow.qlimits.enforcement_mode" => mode, "power_flow.qlimits.guard.enabled" => guard)).config
        res = @test_logs min_level = Logging.Warn redirect_stdout(() -> run_sparlectra(net = n, config = cfg), devnull)
        (n, res)
      end
      for (mode, guard, lock) in (("active_set", true, :guard_lock_pq), ("active_set", false, :none), ("classic_simultaneous", true, :none))
        n, res = narrow_run(mode, guard)
        @test (mode, guard, res.final_converged, res.diagnostics.qlimit_stability_violations) == (mode, guard, true, 0)
        rows = Dict(r.bus_name => r for r in res.diagnostics.qlimit_stability_report)
        bk = rows["Birkkamp_110"]
        @test (mode, bk.end_type, bk.reached_limit, bk.ok, bk.locked_by, bk.q_source) == (mode, :PQ, :max, true, lock, :voltages)
        @test isapprox(bk.q_MVar, 5.5; atol = 1e-6)
        @test bk.u_pu < bk.uset_pu - 1e-3
        wt = rows["Westtal_110"]
        @test (wt.end_type, wt.reached_limit, wt.ok) == (:PV, :none, true)
        @test isapprox(bk.q_report_MVar, bk.q_MVar; atol = 1e-6)
      end
      # Q from the solved voltages counts what the bus draws itself, so a
      # shunt gives the same Q whether it sits in the Y-bus or enters as a
      # voltage-dependent injection: a 1 MVAr shunt at Westtal_110 (PV at
      # 4.306 MVAr), the solved state read again with the shunt switched to
      # the injection model (the state comes from the Y-bus run).
      # Without the shunt term the second reading was 1.040 MVAr higher.
      n, res = narrow_run("active_set", false; shunt_q = 1.0)
      wt_y = only(r for r in res.diagnostics.qlimit_stability_report if r.bus_name == "Westtal_110")
      n.shuntVec[end].model = :VoltageDependentInjection
      wt_v = only(r for r in Sparlectra.qlimitStabilityReport(n) if r.bus_name == "Westtal_110")
      @test (res.final_converged, wt_y.end_type, wt_y.ok, wt_y.q_source, wt_v.q_source) == (true, :PV, true, :voltages, :voltages)
      @test isapprox(wt_v.q_MVar, wt_y.q_MVar; atol = 1e-9)
      @test isapprox(wt_y.q_MVar, wt_y.q_report_MVar; atol = 1e-6)
    end)() end

    # The reason of a run at max_iter (Zhao, Chiang, Li, Ju, PSCC 2008: the
    # numerical divergence and the identification divergence of the PV/PQ
    # switching), from the last three iterations: a type change there at a
    # pre-switch mismatch of at most 10 tol is :type_identification_divergence,
    # no type change, a falling mismatch and a last mismatch below the first
    # of the run :iteration_limit, the rest :numerical_divergence. The status
    # symbol stays :not_converged.
    @testset "Q-limit stop reason at max_iter" begin (function ()
      # the rule alone, tol 1e-6
      @test Sparlectra._max_iter_stop_reason([1.0, 0.1, 0.01], [false, false, false], 1e-6) === :iteration_limit
      @test Sparlectra._max_iter_stop_reason([1.0, 0.1, 0.2], [false, false, false], 1e-6) === :numerical_divergence
      @test Sparlectra._max_iter_stop_reason([1.0, 1e-3, 1e-12], [false, false, true], 1e-6) === :type_identification_divergence
      @test Sparlectra._max_iter_stop_reason([1.0, 0.1, 0.01], [false, false, true], 1e-6) === :numerical_divergence
      # a switch before the window does not count; a run that falls over the
      # window but ends above its first mismatch is no iteration limit
      @test Sparlectra._max_iter_stop_reason([2.0, 1.0, 0.1, 0.01], [true, false, false, false], 1e-6) === :iteration_limit
      @test Sparlectra._max_iter_stop_reason([1e-12, 1.0, 0.1, 0.01], [true, false, false, false], 1e-6) === :numerical_divergence
      @test Sparlectra._max_iter_stop_reason([1.0, 0.1], [false, false], 1e-6) === :numerical_divergence

      # the settings of the PV/PQ study: every automatic machinery off
      # (rescue, start projection, damping, merit, trust region, controls)
      study_cfg(mode, max_iter; guard = true, start_iter = nothing, dc_fallback = false) = begin
        q = Dict{String,Any}("enabled" => true, "enforcement_mode" => mode, "guard" => Dict{String,Any}("enabled" => guard))
        start_iter === nothing || (q["start_mode"] = "iteration"; q["start_iter"] = start_iter)
        SparlectraConfig(Dict{String,Any}("control" => Dict{String,Any}("enabled" => false),
          "power_flow" => Dict{String,Any}("autodamp" => false, "dc" => Dict{String,Any}("fallback" => dc_fallback), "max_iter" => max_iter,
            "merit" => Dict{String,Any}("enabled" => false), "mode" => "manual", "qlimits" => q, "rescue" => false, "solver" => "rectangular",
            "apslf_start" => Dict{String,Any}("enabled" => false), "start_current_iteration" => Dict{String,Any}("enabled" => false),
            "start_mode" => Dict{String,Any}("start_projection" => false), "tol" => 1.0e-6, "trust_region" => Dict{String,Any}("enabled" => false),
            "wrong_branch_detection" => "off")))
      end
      stop_run(net, cfg) = @test_logs min_level = Logging.Warn redirect_stdout(() -> run_sparlectra(net = net, config = cfg), devnull)
      scf = joinpath(dirname(@__DIR__), "data", "scf")
      zeng() = redirect_stdout(() -> Sparlectra.importSCF(joinpath(scf, "case14_zeng_p306_activeSet_A.scf.json")), devnull)

      # case14_zeng, active set switching only on a converged iterate
      # (start_iter 99): max_iter 3 stops at 7.5e-3 with the mismatch falling
      # (3.06, 0.358, 7.5e-3) and no switch; max_iter 5 converges in
      # iteration 5 (1.2e-11 before the switch), then buses 1, 2, 3 and 6
      # switch and the mismatch of the new types is 0.26
      r3 = stop_run(zeng(), study_cfg("active_set", 3; start_iter = 99))
      st3 = Sparlectra.rectangular_pf_status(r3.net)
      @test (r3.outcome, r3.reason, st3.status, st3.type_switch_iterations) == (:not_converged, :iteration_limit, :not_converged, Int[])
      @test isapprox(st3.mismatch_history_pre_switch[end], 7.4993e-3; rtol = 1e-3)
      r5 = stop_run(zeng(), study_cfg("active_set", 5; start_iter = 99))
      st5 = Sparlectra.rectangular_pf_status(r5.net)
      @test (r5.outcome, r5.reason, st5.status, st5.type_switch_iterations) == (:not_converged, :type_identification_divergence, :not_converged, [5])
      # the pre-switch mismatch is kept beside the one after the switch
      @test st5.mismatch_history_pre_switch[end] <= 1e-10 && st5.mismatch_history[end] > 0.1
      @test st5.mismatch_history_pre_switch[1:4] == st5.mismatch_history[1:4]
      # the stability report of the last iterate: the four switched buses
      @test r5.diagnostics.qlimit_stability_status === :last_iterate
      @test [(row.bus, row.end_type) for row in r5.diagnostics.qlimit_stability_report] == [(1, :PQ), (2, :PQ), (3, :PQ), (6, :PQ)]
      @test r5.diagnostics.qlimit_stability_violations == count(row -> !row.ok, r5.diagnostics.qlimit_stability_report)
      @test endswith(qlimit_stability_line(:last_iterate, r5.diagnostics.qlimit_stability_report), "at the last iterate (no converged solution)")

      # the default switching start: the active set switches buses 1, 2, 3
      # and 6 in iteration 3 at 7.5e-3, far above 10 tol, so max_iter 3 is
      # a numerical divergence; with max_iter 5 it releases 2, 3 and 6 in
      # iteration 4 and clamps 3 and 6 again in iteration 5, the releases
      # are in the Q-limit log
      @test stop_run(zeng(), study_cfg("active_set", 3)).reason === :numerical_divergence
      r5d = stop_run(zeng(), study_cfg("active_set", 5))
      @test r5d.reason === :numerical_divergence
      @test [(e.iter, e.bus, e.kind) for e in r5d.net.qLimitLog] ==
            [(3, 1, :hit), (3, 2, :hit), (3, 3, :hit), (3, 6, :hit), (4, 2, :release), (4, 3, :release), (4, 6, :release), (5, 3, :hit), (5, 6, :hit)]
      @test r5d.diagnostics.pv_pq_switching_events == 6
      # the classic modes switch between passes; the base solve at max_iter 3
      # changes no type and falls (3.06, 0.358, 7.5e-3)
      for mode in ("classic_simultaneous", "classic_one_at_a_time")
        @test (mode, stop_run(zeng(), study_cfg(mode, 3)).reason) == (mode, :iteration_limit)
      end

      # sp_case60 with a -1.0 .. 0.5 MVAr band at Birkkamp_110 and
      # Westtal_110 (generators 159 and 160): no strategy converges. The
      # mismatch climbs to 1e+3 .. 1e+6 between the switches; a controlled
      # run keeps the solver's reason, and the net (it has a closed coupler,
      # the solve runs on a contracted copy) carries the Q-limit log of the
      # failed solve. The active set without the guard ends falling with no
      # switch in the last three iterations (540, 135, 33.7 against 2.3e-10
      # in iteration 1): above the first mismatch, so :numerical_divergence
      # (the rule without the comparison with the first mismatch gave
      # :iteration_limit)
      c60 = joinpath(scf, "sp_case60.scf.json")
      for (mode, guard) in (("active_set", true), ("active_set", false), ("classic_simultaneous", true), ("classic_one_at_a_time", true))
        n = redirect_stdout(() -> Sparlectra.importSCF(c60), devnull)
        setQLimits!(net = n, qmin_MVar = -1.0, qmax_MVar = 0.5, busName = ["Birkkamp_110", "Westtal_110"])
        r = stop_run(n, study_cfg(mode, 80; guard = guard))
        st = Sparlectra.rectangular_pf_status(r.net)
        @test (mode, r.outcome, r.reason, st.status) == (mode, :pf_failed, :numerical_divergence, :not_converged)
        @test endswith(r.reason_text, "(controlled framework run)")
        @test (mode, [(e.bus, e.side, e.kind) for e in r.net.qLimitLog[1:2]]) == (mode, [(11, :max, :hit), (25, :max, :hit)])
        @test (mode, r.diagnostics.qlimit_stability_status, length(r.diagnostics.qlimit_stability_report)) == (mode, :last_iterate, 3)
        if !guard
          h = st.mismatch_history_pre_switch
          @test h[end] < h[end-1] < h[end-2] && h[end] > h[1] && all(k -> k < length(h) - 2, st.type_switch_iterations)
        end
      end
      # the same run with power_flow.dc.fallback: the net then holds the DC
      # state (vm = 1 pu), not the last iterate, so the stability check does
      # not run on it
      n = redirect_stdout(() -> Sparlectra.importSCF(c60), devnull)
      setQLimits!(net = n, qmin_MVar = -1.0, qmax_MVar = 0.5, busName = ["Birkkamp_110", "Westtal_110"])
      r = stop_run(n, study_cfg("active_set", 80; guard = false, dc_fallback = true))
      st = Sparlectra.rectangular_pf_status(r.net)
      @test (r.reason, st.dc_fallback_applied, r.diagnostics.qlimit_stability_status, r.diagnostics.qlimit_stability_report) == (:numerical_divergence, true, :not_checked, NamedTuple[])
      @test st.qlimit_stability_line == "Q-limit stability violations: not checked (no converged solution)"
    end)() end

    @testset "Rectangular Q-limit trace diagnostics" begin (function ()
      net = createTest3BusNet()
      setQLimits!(net = net, qmin_MVar = -1.0, qmax_MVar = 1.0, busName = "STATION1")
      trace_path, trace_io = mktemp()
      close(trace_io)
      open(trace_path, "w") do io
        redirect_stdout(io) do
          runpf!(net, 3, 1e-8, 0; method = :rectangular, qlimit_start_iter = 99, qlimit_trace_buses = [2])
        end
      end
      trace_text = read(trace_path, String)
      rm(trace_path; force = true)
      @test occursin("Q-limit trace enabled", trace_text)
      @test occursin("BUS_I=2", trace_text)
      @test occursin("qlimit_start_iter has not been reached", trace_text)
      pq_io = IOBuffer()
      Sparlectra._print_rectangular_qlimit_trace(
        pq_io,
        net,
        4,
        2,
        :PQ,
        ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im],
        [1.0, 1.0, 1.0],
        NaN,
        [-Inf, -0.01, -Inf],
        [Inf, 0.01, Inf];
        q_hyst_pu = 0.0,
        qlimit_start_iter = 2,
        qlimit_start_mode = :iteration_or_auto,
        qlimit_iter_ready = true,
        qlimit_auto_ready = true,
        qlimit_ready = true,
        qlimit_check_active = true,
        converged_before_switching = false,
        qlimit_mode = :switch_to_pq,
        lock_mask = falses(3),
      )
      pq_trace_text = String(take!(pq_io))
      @test occursin("decision=keep PQ", pq_trace_text)
      @test occursin("reason=bus is no longer an active PV/REF bus", pq_trace_text)
      @test !occursin("decision=keep PV reason=bus is not an active PV/REF bus", pq_trace_text)

      ignore_io = IOBuffer()
      Sparlectra._print_rectangular_qlimit_trace(
        ignore_io, net, 1, 2, :PV, ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im], [1.0, 1.0, 1.0], 0.02, [-Inf, -0.01, -Inf], [Inf, 0.01, Inf];
        q_hyst_pu = 0.0, qlimit_start_iter = 1, qlimit_start_mode = :iteration, qlimit_iter_ready = true, qlimit_auto_ready = false, qlimit_ready = true, qlimit_check_active = true, converged_before_switching = false, qlimit_mode = :switch_to_pq, lock_mask = Bool[false, true, false], qlimit_lock_reason = :ignore_q_limits,
      )
      @test occursin("ignore_q_limits=true disables Q-limit switching", String(take!(ignore_io)))

      final_io = IOBuffer()
      res = Sparlectra._print_rectangular_qlimit_summary(
        final_io,
        net,
        ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im],
        ComplexF64[0.0 + 0.0im, 0.0 + 0.02im, 0.0 + 0.0im],
        [:PQ, :PV, :Slack],
        [-Inf, -0.01, -0.05],
        [Inf, 0.01, 0.05],
        zeros(Float64, 3);
        q_hyst_pu = 0.0,
        tolerance_pu = 1e-6,
        max_rows = 30,
      )
      final_text = String(take!(final_io))
      @test res.violating == 1
      @test res.pv_violations == 1
      @test res.ref_violations == 0
      @test occursin("Final PV Q-limit active-set check: FAIL", final_text)
      @test occursin("active PV buses checked", final_text)
      @test occursin("active REF buses checked", final_text)
      @test occursin("PV violations: 1", final_text)
      @test occursin("REF violations: 0", final_text)
      @test occursin("Final PV/REF Q-limit check: FAIL", final_text)
      @test occursin("Qcalc pu", final_text)
      @test occursin("switched", final_text)
      @test occursin("return band", final_text)

      ref_io = IOBuffer()
      ref_res = Sparlectra._print_rectangular_qlimit_summary(
        ref_io,
        net,
        ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im],
        ComplexF64[0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.08im],
        [:PQ, :PV, :Slack],
        [-Inf, -0.01, -0.05],
        [Inf, 0.01, 0.05],
        zeros(Float64, 3);
        q_hyst_pu = 0.0,
        tolerance_pu = 1e-6,
        max_rows = 30,
      )
      ref_text = String(take!(ref_io))
      @test ref_res.pv_violations == 0
      @test ref_res.ref_violations == 1
      @test ref_res.status == :warn
      @test occursin("Final PV Q-limit active-set check: OK", ref_text)
      @test occursin("Final REF Q-limit diagnostic: WARN", ref_text)
      @test occursin("PV violations: 0", ref_text)
      @test occursin("REF violations: 1", ref_text)
      @test occursin("Final PV/REF Q-limit check: WARN", ref_text)

  limited_io = IOBuffer()
  limited_res = Sparlectra._print_rectangular_qlimit_summary(
    limited_io,
    net,
    ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im],
    ComplexF64[0.0 + 0.0im, 0.0 + 0.02im, 0.0 + 0.08im],
    [:PQ, :PV, :Slack],
    [-Inf, -0.01, -0.05],
    [Inf, 0.01, 0.05],
    zeros(Float64, 3);
    q_hyst_pu = 0.0,
    tolerance_pu = 1e-6,
    max_rows = 30,
    max_console_rows = 1,
  )
  limited_text = String(take!(limited_io))
  @test limited_res.violating == 2
  @test occursin("1 more violation rows omitted", limited_text)
  @test occursin("max_rows/max_console_rows", limited_text)

  summary_io = IOBuffer()
  Sparlectra._print_rectangular_convergence_summary(summary_io, (
    numerical_converged = true,
    q_limit_active_set_ok = false,
    final_converged = false,
    reason_text = "remaining PV Q-limit violations",
  ))
  summary_text = String(take!(summary_io))
  @test occursin("numerical_solution=OK", summary_text)
  @test occursin("q_limit_active_set=FAIL", summary_text)
  @test occursin("final_converged=false", summary_text)
  @test occursin("reason=remaining PV Q-limit violations", summary_text)
  @test Sparlectra._rectangular_solver_status_symbol(true, false, false, :remaining_pv_q_limit_violations) == :converged_limits_failed
  @test Sparlectra._rectangular_solver_status_symbol(false, false, false, :singular_newton_step) == :singular_jacobian
  @test Sparlectra._rectangular_solver_status_symbol(false, false, false, :nr_mismatch_not_converged) == :not_converged
end)() end

    @testset "Wrong-branch detection modes are distinct" begin (function ()
      base_cfg = PowerFlowConfig(max_iter = 40, tol = 1e-8, start_mode = StartModeConfig(flatstart = true), wrong_branch_min_vm_pu = 1.20)

      net_warn = createTest3BusNet()
      _, erg_warn = runpf!(net_warn; config = PowerFlowConfig(; NamedTuple{fieldnames(PowerFlowConfig)}(getfield.(Ref(base_cfg), fieldnames(PowerFlowConfig)))..., wrong_branch_detection = :warn))
      st_warn = Sparlectra.rectangular_pf_status(net_warn)
      @test erg_warn == 0
      @test st_warn.wrong_branch_status == :warn

      net_fail = createTest3BusNet()
      _, erg_fail = runpf!(net_fail; config = PowerFlowConfig(; NamedTuple{fieldnames(PowerFlowConfig)}(getfield.(Ref(base_cfg), fieldnames(PowerFlowConfig)))..., wrong_branch_detection = :fail))
      st_fail = Sparlectra.rectangular_pf_status(net_fail)
      @test erg_fail != 0
      @test st_fail.wrong_branch_status == :fail

      net_rescue = createTest3BusNet()
      _, erg_rescue = runpf!(net_rescue; config = PowerFlowConfig(; NamedTuple{fieldnames(PowerFlowConfig)}(getfield.(Ref(base_cfg), fieldnames(PowerFlowConfig)))..., wrong_branch_detection = :rescue))
      st_rescue = Sparlectra.rectangular_pf_status(net_rescue)
      @test erg_rescue != 0
      @test st_rescue.wrong_branch_status == :wrong_branch_rescue_not_implemented
      @test st_rescue.wrong_branch_reason == :rescue_requested_but_not_available
      @test st_rescue.wrong_branch_rescue_attempted === false
    end)() end

    # #456 class (a setting lost on the way to the solver runs on the
    # callee's default without any message): every solver option of the
    # positional runpf! is set once to a value other than its default and
    # must arrive at the keyword runpf_rectangular! on all three paths (plain
    # net, merged link cluster, island-wise solve). A new runpf! keyword must
    # either get a value in the table or be named as local. A process of its
    # own: the probe replaces the keyword method there and stops the solve
    # on arrival.
    @testset "positional runpf! hands every solver option to the solver" begin (function ()
      script = test_scratch_path(".jl")
      write(script, raw"""
      using Sparlectra
      const S = Sparlectra

      # the keyword lists of the three methods on the path, read before the probe
      # replaces the keyword method
      decl(f, sig) = Set(Base.kwarg_decl(only(methods(f, sig))))
      kw_runpf = decl(S.runpf!, (S.Net, Int, Float64, Int))
      kw_rect_pos = decl(S.runpf_rectangular!, (S.Net, Int, Float64, Int))
      kw_rect = decl(S.runpf_rectangular!, (S.Net,))

      # keywords the positional runpf! acts on itself and does not hand on
      runpf_local = Set([:method, :validate_limits_after_pf, :q_limit_violation_headroom, :auto_slack, :islands_enabled, :islands_mode, :islands_reference_policy,
        :islands_diagnostic_continue_after_failure, :islands_parallel_enabled, :islands_parallel_max_tasks, :islands_parallel_min_work_items, :islands_parallel_min_buses])

      function feeder!(net, k)
        a, b, c = "A$(k)", "B$(k)", "C$(k)"
        foreach(bus -> addBus!(net = net, busName = bus, vn_kV = 110.0), (a, b, c))
        addProsumer!(net = net, busName = a, type = "EXTERNALNETWORKINJECTION", referencePri = a, vm_pu = 1.0, va_deg = 0.0)
        addProsumer!(net = net, busName = b, type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
        addProsumer!(net = net, busName = c, type = "ENERGYCONSUMER", p = 5.0, q = 1.0)
        addPIModelACLine!(net = net, fromBus = a, toBus = b, r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = net, fromBus = b, toBus = c, r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      end
      function build(kind)
        net = Net(name = "probe_$(kind)", baseMVA = 100.0)
        feeder!(net, 1)
        kind === :island && feeder!(net, 2)
        if kind === :merged
          addBus!(net = net, busName = "C1X", vn_kV = 110.0)
          addLink!(net = net, fromBus = "C1", toBus = "C1X", status = 1)
        end
        validate!(net = net)
        return net
      end

      struct ProbeStop <: Exception end
      const RECEIVED = Ref{Any}(nothing)
      # the probe takes the place of the keyword method: it records what arrives
      # and stops the solve
      Core.eval(S, :(function runpf_rectangular!(net::Net; kw...)
        Main.RECEIVED[] = NamedTuple(kw)
        throw(Main.ProbeStop())
      end))
      function received(f)
        RECEIVED[] = nothing
        try
          f()
        catch err
          err isa ProbeStop || rethrow()
        end
        return RECEIVED[]
      end

      defaults = received(() -> runpf!(build(:plain), 10))
      profile = Dict{Symbol,Any}()
      nondefault = (;
        damp = 0.5, autodamp = true, autodamp_min = 0.07, newton_update = :rectangular, power_mode = true, power_mode_lu = :klu,
        jacobian_reuse = true, jacobian_reuse_min_reduction = 12.5, jacobian_reuse_max_steps = 4, qlimit_warm_clamps = Dict("G1" => :max),
        merit_enabled = true, merit_armijo_c1 = 2.0e-4, merit_scale_p = 1.1, merit_scale_q = 1.2, merit_scale_v = 1.3, merit_fallback_max_mismatch = false,
        trust_region_enabled = true, trust_region_initial_radius = 0.9, trust_region_min_radius = 2.0e-4, trust_region_max_radius = 9.0, trust_region_eta_accept = 0.15,
        trust_region_shrink_factor = 0.4, trust_region_expand_factor = 2.5, trust_region_expand_threshold = 0.8, trust_region_step_mode = :dogleg,
        opt_flatstart = !defaults.opt_flatstart, pv_table_rows = 31,
        lock_pv_to_pq_buses = [2], qlimit_mode = :adjust_vset, qlimit_max_outer = 7,
        start_projection = true, start_projection_try_dc_start = false, start_projection_try_blend_scan = false, start_projection_branch_guard = false,
        start_projection_measure_candidates = false, start_projection_accept_unmeasured_dc_start = true, start_projection_blend_lambdas = [0.3],
        start_projection_dc_angle_limit_deg = 45.0, start_projection_requested_angle_mode = :dc, start_projection_requested_voltage_mode = :all_bus_vm,
        start_projection_ratio_profile = true, start_rotate_behind_pst = true,
        start_current_iteration_enabled = true, start_current_iteration_max_iter = 11, start_current_iteration_tol = 2.0e-3, start_current_iteration_damping = 0.6,
        start_current_iteration_accept_only_if_improved = false, start_current_iteration_min_improvement_factor = 0.97, start_current_iteration_vm_min_pu = 0.55,
        start_current_iteration_vm_max_pu = 1.45, start_current_iteration_max_angle_step_deg = 25.0, start_current_iteration_only_for_large_cases = true,
        apslf_start_enabled = true, apslf_start_order = 20,
        qlimit_start_iter = defaults.qlimit_start_iter + 2, qlimit_start_mode = :auto, qlimit_auto_q_delta_pu = 3.0e-4, qlimit_trace_buses = [1], qlimit_lock_reason = :ignore_q_limits,
        qlimit_guard = !defaults.qlimit_guard, qlimit_guard_min_q_range_pu = 0.03, qlimit_guard_zero_range_mode = :prefer_pq, qlimit_guard_narrow_range_mode = :delayed_switch,
        qlimit_guard_log = false, qlimit_guard_max_switches = defaults.qlimit_guard_max_switches + 1, qlimit_guard_accept_bounded_violations = true,
        qlimit_guard_max_remaining_violations = 2, qlimit_guard_freeze_after_repeated_switching = false, qlimit_guard_violation_mode = :ignore,
        qlimit_guard_violation_threshold_pu = 2.0e-4, qlimits_enabled = false, qlimit_enforcement_mode = :classic_simultaneous,
        wrong_branch_detection = :fail, wrong_branch_rescue = true, wrong_branch_min_vm_pu = 0.65, wrong_branch_max_vm_pu = 1.35, wrong_branch_max_angle_spread_deg = 170.0,
        wrong_branch_max_branch_angle_deg = 80.0, wrong_branch_min_low_vm_count = 2, wrong_branch_min_vn_kV = 90.0, wrong_branch_low_vm_share = 0.06,
        wrong_branch_max_bus_angle_deg = 110.0, wrong_branch_max_plain_steps = 15, wrong_branch_collapse_vm_pu = 0.45, wrong_branch_max_reference_branch_angle_deg = 85.0,
        wrong_branch_rescue_max_attempts = 3,
        rectangular_workspace_reuse = false, rectangular_preallocate_workspace = :on, rectangular_workspace_min_buses = 500, linear_solver = :umfpack,
        performance_profile = profile,
        distributed_slack_enabled = true, distributed_slack_p_mode = :pmax_weighted, distributed_slack_respect_p_limits = false, distributed_slack_fallback = :ref_only,
        distributed_slack_weights = Dict("A1" => 1.0),
      )
      forwarded = setdiff(kw_runpf, runpf_local)
      names(s) = isempty(s) ? "none" : join(sort(collect(s)), ", ")
      # every keyword of runpf! is named as local or has a value in the table, and
      # the positional and the keyword runpf_rectangular! declare each forwarded one
      println("unclassified: ", names(symdiff(Set(keys(nondefault)), forwarded)), "; local not declared: ", names(setdiff(runpf_local, kw_runpf)))
      println("not declared on the path: ", names(union(setdiff(forwarded, kw_rect_pos), setdiff(forwarded, kw_rect))))
      println("arrive without keywords: ", names(symdiff(Set(keys(defaults)), union(forwarded, Set([:method, :maxiter, :tol, :verbose])))), " extra or missing")
      println("table values equal to the default: ", names([k for k in keys(nondefault) if isequal(nondefault[k], defaults[k])]))
      # the island call hands its own values for these four
      island_override = (; wrong_branch_rescue = false, wrong_branch_rescue_max_attempts = 0, lock_pv_to_pq_buses = [2], performance_profile = profile)
      for kind in (:plain, :merged, :island)
        got = received(() -> runpf!(build(kind), 11, 1.0e-5, 1; nondefault..., islands_enabled = kind === :island))
        expect = kind === :island ? merge(nondefault, island_override) : nondefault
        lost = [k for k in keys(nondefault) if !haskey(got, k) || !(k === :performance_profile ? got[k] === expect[k] : isequal(got[k], expect[k]))]
        println("path ", kind, ": ", length(keys(nondefault)) - length(lost), " options arrive, lost: ", names(lost), ", positional ", (got.maxiter, got.tol, got.verbose, got.method) == (11, 1.0e-5, 1, :rectangular))
      end
      """)
      library = normpath(joinpath(@__DIR__, ".."))
      buffer = IOBuffer()
      proc = run(pipeline(Cmd(`$(Base.julia_cmd()) --startup-file=no --project=$(library) $(script)`; ignorestatus = true); stdout = buffer, stderr = buffer))
      out = String(take!(buffer))
      line(prefix) = (m = match(Regex("^" * prefix * "(.*)\$", "m"), out); m === nothing ? "missing in the child output:\n" * out : String(m.captures[1]))
      @test success(proc)
      @test line("unclassified: ") == "none; local not declared: none"
      @test line("not declared on the path: ") == "none"
      @test line("arrive without keywords: ") == "none extra or missing"
      @test line("table values equal to the default: ") == "none"
      for kind in ("plain", "merged", "island")
        @test endswith(line("path $(kind): "), " options arrive, lost: none, positional true")
      end
    end)() end

    @testset "Typed power-flow config entry points" begin (function ()
      # no `sparse` keyword any more: sparse matrices are mandatory since
      # 0.9.x, so the field was deleted rather than kept as a constant true
      pf_config = PowerFlowConfig(max_iter = 40, tol = 1e-8, start_mode = StartModeConfig(flatstart = true))

      net_direct = createTest3BusNet()
      _, erg_direct = runpf!(net_direct; config = pf_config)
      @test erg_direct == 0

      net_project = createTest3BusNet()
      _, erg_project = runpf!(net_project; config = SparlectraConfig(powerflow = pf_config))
      @test erg_project == 0

      net_rectangular = createTest3BusNet()
      _, erg_rectangular = runpf_rectangular!(net_rectangular, 20, 1e-8, 0; method = :rectangular)
      @test erg_rectangular == 0

      cfg_output_off = SparlectraConfig(powerflow = pf_config, output = OutputConfig(logfile_results = :off))
      net_runner = createTest3BusNet()
      result = run_sparlectra(net = net_runner, config = cfg_output_off)
      @test result isa SparlectraRunResult
      @test result.net === net_runner
      @test result.numerical_converged
      @test result.solution_available
      @test result.control_status === :none
      @test result.final_converged
      @test result.outcome in (:converged, :converged_with_limit_warnings)
      @test result.limit_validation_status in (:ok, :skip)
      @test result.reason isa Symbol
      @test result.reason_text isa String
      @test result.iterations >= 0
      @test result.final_mismatch isa Float64
      @test result.diagnostics isa NamedTuple

      fallback_result = SparlectraRunResult(
        result.net,
        :converged,
        true,
        true,
        :skip,
        true,
        :none,
        "none",
        0,
        0.0,
        nothing,
        NaN,
        :classic,
        :none,
        nothing,
        (converged = false, erg = 1, numerical_converged = false, solution_available = false, final_converged = false, outcome = :diagnostic_override, reason = :diagnostic_override),
      )
      fallback_row = Sparlectra._sparlectra_status_row(fallback_result)
      @test fallback_row.converged === true
      @test fallback_row.erg == 0
      @test fallback_row.numerical_converged === true
      @test fallback_row.solution_available === true
      @test fallback_row.final_converged === true
      @test fallback_row.outcome === :converged
      @test fallback_row.reason === :none
      @test !hasproperty(fallback_row, :q_limit_active_set_ok)

      net_alias = createTest3BusNet()
      alias_result = run_acpflow(net = net_alias, config = cfg_output_off)
      @test alias_result isa SparlectraRunResult
      @test alias_result.net === net_alias
      @test alias_result.outcome == result.outcome
      @test alias_result.method == result.method

      @test_throws TypeError run_sparlectra(net = createTest3BusNet(), config = pf_config)
      @test_throws TypeError run_acpflow(net = createTest3BusNet(), config = pf_config)

      net_with_explicit_cfg = createTest3BusNet()
      explicit_output = mktemp() do path, io
        redirect_stdout(io) do
          @test run_sparlectra(net = net_with_explicit_cfg, config = cfg_output_off).numerical_converged
        end
        flush(io)
        return read(path, String)
      end
      @test !occursin("AC Power Flow Results", explicit_output)

      cfg_output_full = SparlectraConfig(powerflow = pf_config, output = OutputConfig(logfile_results = :full))
      net_with_output = createTest3BusNet()
      full_output = mktemp() do path, io
        redirect_stdout(io) do
          @test run_sparlectra(net = net_with_output, config = cfg_output_full).numerical_converged
        end
        flush(io)
        return read(path, String)
      end
      @test occursin("AC Power Flow Results", full_output)
      @test Sparlectra._sparlectra_result_mode(net_with_output, OutputConfig(logfile_results = :compact, result_table_large_case_threshold_buses = 1, result_table_large_case_mode = :summary)) === :summary
    end)() end

    @testset "Framework status composition" begin (function ()
      rejected_outcomes = (:wrong_branch_detected, :angle_spread_exceeded, :branch_angle_exceeded, :wrong_branch_rescue_not_implemented)
      for outcome in rejected_outcomes
        rect_status = (status = outcome, numerical_converged = true, reason = :none)
        @test !Sparlectra._rectangular_solution_available(rect_status)
      end
      rejected_reasons = (:wrong_branch_detected, :angle_spread_exceeded, :branch_angle_exceeded, :wrong_branch_rescue_not_implemented, :rescue_requested_but_not_available)
      for reason in rejected_reasons
        rect_status = (status = :converged, numerical_converged = true, reason = reason)
        @test !Sparlectra._rectangular_solution_available(rect_status)
      end
      @test !Sparlectra._rectangular_solution_available(nothing)
      @test !Sparlectra._rectangular_solution_available((status = :not_converged, numerical_converged = false, reason = :nr_mismatch_not_converged))
      @test Sparlectra._rectangular_solution_available((status = :converged_limits_failed, numerical_converged = true, reason = :remaining_q_limit_violations))

      base = (outcome = :converged, numerical_converged = true, solution_available = true, limit_validation_status = :ok, final_converged = true, reason = :none, reason_text = "none", final_mismatch = 0.0)
      for accepted_status in (:none, :converged, :disabled, :no_active_controllers, :no_controllers)
        accepted = Sparlectra._compose_framework_status(base, accepted_status)
        @test accepted.final_converged
        @test accepted.outcome === :converged
      end

      blocked = Sparlectra._compose_framework_status(base, :blocked)
      @test blocked.numerical_converged
      @test blocked.solution_available
      @test !blocked.final_converged
      @test blocked.outcome === :control_blocked
      @test blocked.reason === :control_blocked
      @test occursin("blocked", blocked.reason_text)

      exhausted = Sparlectra._compose_framework_status(base, :max_outer_iterations)
      @test exhausted.outcome === :control_max_outer_iterations
      @test !exhausted.final_converged

      unknown = Sparlectra._compose_framework_status(base, :custom_terminal_status)
      @test unknown.outcome === :control_not_converged
      @test !unknown.final_converged

      rejected_base = merge(base, (outcome = :wrong_branch_detected, solution_available = false, final_converged = false, reason = :wrong_branch_detected, reason_text = "wrong branch detected"))
      rejected = Sparlectra._compose_framework_status(rejected_base, :pf_failed)
      @test rejected.numerical_converged
      @test !rejected.solution_available
      @test !rejected.final_converged
      @test rejected.outcome === :wrong_branch_detected
      @test rejected.reason === :wrong_branch_detected
      @test occursin("wrong branch", lowercase(rejected.reason_text))

      # a failed solve keeps the solver's reason under the control frame,
      # the outcome is :pf_failed; a status without a reason gets :pf_failed
      failed_base = merge(base, (outcome = :not_converged, numerical_converged = false, solution_available = false, final_converged = false, reason = :numerical_divergence, reason_text = "NR mismatch did not converge (numerical divergence)", final_mismatch = Inf))
      failed = Sparlectra._compose_framework_status(failed_base, :pf_failed)
      @test !failed.numerical_converged
      @test !failed.solution_available
      @test !failed.final_converged
      @test failed.outcome === :pf_failed
      @test failed.reason === :numerical_divergence
      @test failed.reason_text == "NR mismatch did not converge (numerical divergence) (controlled framework run)"
      unknown_failed = Sparlectra._compose_framework_status(merge(failed_base, (reason = :none, reason_text = "none")), :pf_failed)
      @test (unknown_failed.outcome, unknown_failed.reason) == (:pf_failed, :pf_failed)
    end)() end

    @testset "Configured MATPOWER batch parsing and deterministic local execution" begin (function ()
      cfg_cases = SparlectraConfig(Dict("runtime" => Dict("cases" => [" case_a.m ", "case_b.m"])))
      @test cfg_cases.runtime.cases == ["case_a.m", "case_b.m"]
      @test configured_matpower_cases(cfg_cases) == ["case_a.m", "case_b.m"]

      cfg_single = SparlectraConfig(Dict("runtime" => Dict("case" => " case_single.m ")))
      @test configured_matpower_cases(cfg_single) == ["case_single.m"]

      cfg_precedence = SparlectraConfig(Dict("runtime" => Dict("case" => "case_single.m", "cases" => ["case_a.m", "case_b.m"])))
      @test configured_matpower_cases(cfg_precedence) == ["case_a.m", "case_b.m"]
      @test_throws ArgumentError SparlectraConfig(Dict("runtime" => Dict("cases" => "case_a.m")))
      @test_throws ArgumentError SparlectraConfig(Dict("runtime" => Dict("cases" => ["case_a.m", " "])))
      @test_throws ArgumentError run_sparlectra_cases(config = SparlectraConfig())
      @test_throws ArgumentError run_sparlectra_cases(config = cfg_single, performance_profile = Dict{Symbol,Any}())

      mktempdir() do tmpdir
        for (case_name, load_mw) in [("case_a.m", 40), ("case_b.m", 60)]
          function_name = splitext(case_name)[1]
          write(joinpath(tmpdir, case_name), """
function mpc = $(function_name)
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
1 3 0 0 0 0 1 1.0 0 110 1 1.1 0.9;
2 1 $(load_mw) 15 0 0 1 1.0 0 110 1 1.1 0.9;
];
mpc.gen = [
1 100 0 300 -300 1.02 100 1 300 0;
];
mpc.branch = [
1 2 0.01 0.05 0.0 999 999 999 0 0 1 -360 360;
];
""")
        end
        active_before = active_sparlectra_config()
        results = run_sparlectra_cases(config = cfg_precedence, path = tmpdir)
        @test length(results) == 2
        @test [result.net.name for result in results] == ["case_a", "case_b"]
        @test all(result -> result isa SparlectraRunResult, results)
        @test all(result -> result.numerical_converged, results)
        @test active_sparlectra_config() === active_before
      end
    end)() end

    @testset "Framework runner input validation" begin (function ()
      net = createTest3BusNet()
      @test_throws ArgumentError run_sparlectra(net = net, casefile = "case3.m")
      @test_throws ArgumentError run_sparlectra()
      @test_throws ArgumentError run_sparlectra(net = net, path = "unused")
      @test_throws ArgumentError run_acpflow(net = net, casefile = "case3.m")
      @test_throws ArgumentError run_acpflow()
      @test_throws ArgumentError run_acpflow(net = net, path = "unused")
    end)() end

    @testset "Rectangular damping validation" begin (function ()
      @test_throws ErrorException Sparlectra._validate_rectangular_damping(1.0, 0.0)
      @test_throws ErrorException Sparlectra._validate_rectangular_damping(0.2, 0.3)
      @test_throws ErrorException Sparlectra._validate_rectangular_damping(0.0, 0.01)
      @test_throws ErrorException Sparlectra._validate_rectangular_damping(1.01, 0.05)
    end)() end

    @testset "Merit value and Armijo acceptance (unit)" begin (function ()
      Y = ComplexF64[0.0 - 10.0im 0.0 + 10.0im; 0.0 + 10.0im 0.0 - 10.0im]
      V = ComplexF64[1.0 + 0.0im, 1.0 + 0.0im]
      S = ComplexF64[0.0 + 0.0im, -1.0 - 0.2im]
      bus_types = [:Slack, :PQ]
      Vset = [1.0, 1.0]
      F0 = Sparlectra.mismatch_rectangular(Y, V, S, bus_types, Vset, 1)
      non_slack = Sparlectra.non_slack_indices(2, 1)

      expected_f0 = 0.5 * (F0[1]^2 + F0[2]^2)
      @test Sparlectra._rectangular_merit_value(F0, non_slack, bus_types, 1.0, 1.0, 1.0) ≈ expected_f0
      scaled = Sparlectra._rectangular_merit_value(F0, non_slack, bus_types, 2.0, 1.0, 1.0)
      @test scaled ≈ 0.5 * ((2.0 * F0[1])^2 + F0[2]^2)

      oversized_delta = [0.0, -1.0]
      merit_log = NamedTuple[]
      alpha, Vtrial, trial_mismatch = Sparlectra.choose_rectangular_autodamp(
        Y, V, S, oversized_delta, F0;
        slack_idx = 1, damp = 1.0, autodamp_min = 1e-3, bus_types = bus_types, Vset = Vset,
        merit_enabled = true, armijo_c1 = 1.0e-4, fallback_max_mismatch = true, merit_log = merit_log,
      )
      @test length(merit_log) == 1
      @test merit_log[1].accept_reason == :armijo
      @test merit_log[1].directional_derivative ≈ -2 * merit_log[1].f_before
      @test merit_log[1].directional_derivative < 0

      # merit-disabled (default) path stays bit-identical to the classic algorithm.
      alpha_classic, Vtrial_classic, trial_mismatch_classic = Sparlectra.choose_rectangular_autodamp(
        Y, V, S, oversized_delta, F0; slack_idx = 1, damp = 1.0, autodamp_min = 1e-3, bus_types = bus_types, Vset = Vset,
      )
      @test alpha == alpha_classic
      @test Vtrial == Vtrial_classic
      @test trial_mismatch == trial_mismatch_classic

      # active_set_changed=true skips the merit comparison for this call.
      skip_log = NamedTuple[]
      Sparlectra.choose_rectangular_autodamp(
        Y, V, S, oversized_delta, F0;
        slack_idx = 1, damp = 1.0, autodamp_min = 1e-3, bus_types = bus_types, Vset = Vset,
        merit_enabled = true, active_set_changed = true, merit_log = skip_log,
      )
      @test skip_log[1].accept_reason == :active_set_skip

      # fallback_max_mismatch=false skips straight to the conservative fallback when Armijo fails.
      strict_log = NamedTuple[]
      Sparlectra.choose_rectangular_autodamp(
        Y, V, S, oversized_delta, F0;
        slack_idx = 1, damp = 1.0, autodamp_min = 1e-3, bus_types = bus_types, Vset = Vset,
        merit_enabled = true, armijo_c1 = 0.4999, fallback_max_mismatch = false, merit_log = strict_log,
      )
      @test strict_log[1].accept_reason in (:armijo, :fallback_conservative)
    end)() end

    @testset "Functional: merit.enabled=true converges and logs an armijo acceptance" begin (function ()
      mktempdir() do tmpdir
        perf = Dict{Symbol,Any}(:output_dir => tmpdir)
        net = createCIGRE()
        _, erg = runpf_rectangular!(net, 30, 1e-8, 0; autodamp = true, merit_enabled = true, performance_profile = perf)
        @test erg == 0
        st = Sparlectra.rectangular_pf_status(net)
        @test st.merit_enabled == true
        @test st.merit_used_iterations >= 1

        logpath = joinpath(tmpdir, "merit_linesearch.log")
        @test isfile(logpath)
        log_text = read(logpath, String)
        @test occursin("accept_reason=armijo", log_text)
        @test occursin("merit_enabled: true", log_text)
      end
    end)() end

    @testset "Active-set change skips merit comparison for that iteration" begin (function ()
      net = createTest3BusNet()
      setQLimits!(net = net, qmin_MVar = -1.0, qmax_MVar = 1.0, busName = "STATION1")
      _, erg = runpf!(net, 20, 1e-6, 0; method = :rectangular, autodamp = true, merit_enabled = true)
      @test erg == 0
      st = Sparlectra.rectangular_pf_status(net)
      @test st.merit_active_set_skip_count > 0
    end)() end

    @testset "Config validation (merit and trust_region)" begin (function ()
      # power_flow blocks the loader must refuse
      refused = [
        "  autodamp: true\n  merit:\n    armijo_c1: 0.5\n",
        "  autodamp: true\n  merit:\n    scale_p: -1.0\n",
        # merit needs autodamp
        "  autodamp: false\n  merit:\n    enabled: true\n",
        "  trust_region:\n    initial_radius: 1.0\n    min_radius: 1.0\n",
        "  trust_region:\n    shrink_factor: 1.5\n",
        "  trust_region:\n    expand_factor: 0.5\n",
        # autodamp and the trust region exclude each other
        "  autodamp: true\n  trust_region:\n    enabled: true\n",
        "  trust_region:\n    step_mode: bogus\n",
      ]
      for block in refused
        path = test_scratch_path(".yaml")
        write(path, "config_version: 1\npower_flow:\n" * block)
        @test_throws ArgumentError Sparlectra.load_sparlectra_config(path; reload = true)
      end
      # power_flow blocks that load, with the values they must give
      accepted = [
        ("  autodamp: true\n  merit:\n    enabled: true\n", pf -> pf.merit.enabled == true),
        ("  autodamp: false\n  trust_region:\n    enabled: true\n    initial_radius: 2.0\n",
         pf -> pf.trust_region.enabled == true && pf.trust_region.initial_radius == 2.0 && pf.trust_region.step_mode === :scaled),
        ("  autodamp: false\n  trust_region:\n    enabled: true\n    step_mode: dogleg\n", pf -> pf.trust_region.step_mode === :dogleg),
      ]
      for (block, check) in accepted
        path = test_scratch_path(".yaml")
        write(path, "config_version: 1\npower_flow:\n" * block)
        @test (block, check(Sparlectra.load_sparlectra_config(path; reload = true).powerflow)) == (block, true)
      end
    end)() end

    @testset "Unit: dogleg gradient, Cauchy point, and path branch selection" begin (function ()
      # Simple diagonal J so g = Jᵀ F, Jg, α*, and p_C are hand-checkable.
      J = sparse([2.0 0.0; 0.0 4.0])
      F0 = [1.0, 2.0]
      g = Sparlectra._rectangular_merit_gradient(J, F0)
      @test g == [2.0, 8.0]  # Jᵀ F with J diagonal: g_i = J_ii * F_i

      p_C = Sparlectra._rectangular_cauchy_point(J, g)
      Jg = J * g
      expected_alpha = dot(g, g) / dot(Jg, Jg)
      @test p_C ≈ (-expected_alpha) .* g
      # Steepest descent: p_C must point exactly opposite to g.
      @test p_C ≈ -norm(p_C) / norm(g) .* g

      p_N = [3.0, 0.0]
      pN_norm = norm(p_N)
      pC_norm = norm(p_C)

      # Branch 1: Newton step already inside a generous radius.
      step1, reason1 = Sparlectra._rectangular_dogleg_step(p_N, pN_norm, p_C, pC_norm, 100.0)
      @test reason1 === :dogleg_newton
      @test step1 == p_N

      # Branch 2: even the Cauchy point exceeds a tiny radius.
      tiny_radius = pC_norm / 10
      step2, reason2 = Sparlectra._rectangular_dogleg_step(p_N, pN_norm, p_C, pC_norm, tiny_radius)
      @test reason2 === :dogleg_cauchy
      @test norm(step2) ≈ tiny_radius
      @test step2 ≈ (tiny_radius / pC_norm) .* p_C

      # Branch 3: radius strictly between ‖p_C‖ and ‖p_N‖ forces interpolation,
      # and the interpolated step must land exactly on the radius boundary.
      mid_radius = (pC_norm + pN_norm) / 2
      @test pC_norm < mid_radius < pN_norm
      step3, reason3 = Sparlectra._rectangular_dogleg_step(p_N, pN_norm, p_C, pC_norm, mid_radius)
      @test reason3 === :dogleg_interp
      @test norm(step3) ≈ mid_radius atol = 1e-9
    end)() end

    @testset "Unit: dogleg active_set_changed skips the comparison without moving the radius" begin (function ()
      Y = ComplexF64[0.0 - 10.0im 0.0 + 10.0im; 0.0 + 10.0im 0.0 - 10.0im]
      V = ComplexF64[1.0 + 0.0im, 1.0 + 0.0im]
      S = ComplexF64[0.0 + 0.0im, -1.0 - 0.2im]
      bus_types = [:Slack, :PQ]
      Vset = [1.0, 1.0]
      F0 = Sparlectra.mismatch_rectangular(Y, V, S, bus_types, Vset, 1)
      J = Sparlectra.build_rectangular_jacobian_pq_pv(sparse(Y), V, bus_types, Vset, 1; dPinj_dVm = zeros(2), dQinj_dVm = zeros(2))
      δx = Sparlectra.solve_linear(J, -F0; allow_pinv = true)

      radius_ref = Ref(1.0)
      tr_log = NamedTuple[]
      Vout = Sparlectra.choose_rectangular_trust_region_step(
        Y, V, S, δx, F0, J;
        slack_idx = 1, bus_types = bus_types, Vset = Vset, radius_ref = radius_ref,
        min_radius = 1e-4, max_radius = 10.0, eta_accept = 0.1, shrink_factor = 0.5,
        expand_factor = 2.0, expand_threshold = 0.75, diagnostics = tr_log,
        step_mode = :dogleg, active_set_changed = true,
      )
      @test Vout isa Vector{ComplexF64}
      @test length(tr_log) == 1
      @test tr_log[1].accept_reason === :active_set_skip
      @test tr_log[1].accepted == true
      @test radius_ref[] == 1.0  # radius left unchanged by the active-set-skip path
    end)() end

    @testset "Unit: choose_rectangular_trust_region_step acceptance and collapse" begin (function ()
      Y = ComplexF64[0.0 - 10.0im 0.0 + 10.0im; 0.0 + 10.0im 0.0 - 10.0im]
      V = ComplexF64[1.0 + 0.0im, 1.0 + 0.0im]
      S = ComplexF64[0.0 + 0.0im, -1.0 - 0.2im]
      bus_types = [:Slack, :PQ]
      Vset = [1.0, 1.0]
      F0 = Sparlectra.mismatch_rectangular(Y, V, S, bus_types, Vset, 1)
      J = Sparlectra.build_rectangular_jacobian_pq_pv(sparse(Y), V, bus_types, Vset, 1; dPinj_dVm = zeros(2), dQinj_dVm = zeros(2))
      δx = Sparlectra.solve_linear(J, -F0; allow_pinv = true)

      radius_ref = Ref(1.0)
      tr_log = NamedTuple[]
      Vout = Sparlectra.choose_rectangular_trust_region_step(
        Y, V, S, δx, F0, J;
        slack_idx = 1, bus_types = bus_types, Vset = Vset, radius_ref = radius_ref,
        min_radius = 1e-4, max_radius = 10.0, eta_accept = 0.1, shrink_factor = 0.5,
        expand_factor = 2.0, expand_threshold = 0.75, diagnostics = tr_log,
      )
      @test Vout isa Vector{ComplexF64}
      @test length(tr_log) == 1
      @test tr_log[1].accepted == true
      @test tr_log[1].collapsed == false

      collapse_radius_ref = Ref(1e-5)
      collapse_log = NamedTuple[]
      @test_throws Sparlectra.RectangularTrustRegionCollapsed Sparlectra.choose_rectangular_trust_region_step(
        Y, V, S, δx, F0, J;
        slack_idx = 1, bus_types = bus_types, Vset = Vset, radius_ref = collapse_radius_ref,
        min_radius = 1e-4, max_radius = 10.0, eta_accept = 0.1, shrink_factor = 0.5,
        expand_factor = 2.0, expand_threshold = 0.75, diagnostics = collapse_log,
      )
      @test length(collapse_log) == 1
      @test collapse_log[1].collapsed == true
    end)() end

    @testset "Functional: trust_region.enabled=true converges and logs accepted steps" begin (function ()
      # per step mode the accept reason its log must show, and for dogleg
      # the scaled reason it must not show
      for (step_mode, accepted, absent) in ((:scaled, "accept_reason=scaled", nothing), (:dogleg, r"accept_reason=dogleg_(newton|interp|cauchy)", "accept_reason=scaled"))
        mktempdir() do tmpdir
          perf = Dict{Symbol,Any}(:output_dir => tmpdir)
          net = createCIGRE()
          _, erg = runpf_rectangular!(net, 30, 1e-8, 0; autodamp = false, trust_region_enabled = true, trust_region_step_mode = step_mode, performance_profile = perf)
          @test erg == 0
          st = Sparlectra.rectangular_pf_status(net)
          @test st.trust_region_enabled == true
          @test st.tr_step_count >= 1
          @test st.tr_collapsed == false

          logpath = joinpath(tmpdir, "trust_region.log")
          @test isfile(logpath)
          log_text = read(logpath, String)
          @test occursin("accepted=true", log_text)
          @test occursin("trust_region_enabled: true", log_text)
          @test occursin(accepted, log_text)
          absent === nothing || @test !occursin(absent, log_text)
        end
      end
    end)() end

    @testset "Radius collapse reports :trust_region_collapsed on an unsolvable case" begin (function ()
      net = Net(name = "tr_collapse_demo", baseMVA = 100.0, flatstart = true)
      addBus!(net = net, busName = "SLACK", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
      addBus!(net = net, busName = "LOAD", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
      addPIModelACLine!(net = net, fromBus = "SLACK", toBus = "LOAD", r_pu = 0.01, x_pu = 0.5, b_pu = 0.0, status = 1)
      addProsumer!(net = net, busName = "SLACK", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "SLACK")
      addProsumer!(net = net, busName = "LOAD", type = "ENERGYCONSUMER", p = 200.0, q = 60.0)
      ok, _ = validate!(net = net)
      @test ok

      _, erg = runpf_rectangular!(net, 30, 1e-8, 0; autodamp = false, trust_region_enabled = true, trust_region_min_radius = 1e-3)
      @test erg != 0
      st = Sparlectra.rectangular_pf_status(net)
      @test st.reason == :trust_region_collapsed
      @test st.tr_collapsed == true
    end)() end


    @testset "Current-iteration rejection log reports candidate guard details" begin (function ()
      mktempdir() do tmpdir
        run_guarded_current_iteration_start = getfield(Sparlectra, :_run_guarded_current_iteration_start)
        Ybus = ComplexF64[1 -1; -1 1]
        original = ComplexF64[1.0 + 0.0im, 1.0 + 0.0im]
        S = ComplexF64[0.0 + 0.0im, 10.0 + 0.0im]
        bus_types = [:Slack, :PQ]
        profile = Dict{Symbol, Any}(:output_dir => tmpdir)

        restored, summary = run_guarded_current_iteration_start(
          Ybus,
          original,
          S,
          bus_types,
          [1.0, 1.0],
          1;
          enabled = true,
          max_iter = 1,
          tol = 1.0e-6,
          damping = 1.0,
          accept_only_if_improved = true,
          min_improvement_factor = 0.98,
          vm_min_pu = 0.5,
          vm_max_pu = 1.5,
          max_angle_step_deg = 30.0,
          only_for_large_cases = false,
          large_case_min_buses = 10,
          performance_profile = profile,
        )

        @test restored == original
        @test summary.current_iteration_attempted === true
        @test summary.current_iteration_accepted === false
        @test summary.current_iteration_reason === :voltage_magnitude_guard
        log_text = read(joinpath(tmpdir, "current_iteration_start.log"), String)
        for needle in (
          "current_iteration_attempted: true",
          "current_iteration_accepted: false",
          "current_iteration_reason: voltage_magnitude_guard",
          "candidate_voltage_magnitude_min:",
          "candidate_voltage_magnitude_max: 11.0",
          "candidate_voltage_low_count:",
          "candidate_voltage_high_count: 1",
          "candidate_voltage_worst_high_bus: 2",
          "candidate_voltage_worst_high_value: 11.0",
          "candidate_max_angle_step_deg:",
          "rejected_at_iteration: 1",
          "rejection_stage: candidate_guard",
          "original_start_values_restored: true",
          "restored_voltage_magnitude_min: 1.0",
          "restored_voltage_magnitude_max: 1.0",
        )
          @test occursin(needle, log_text)
        end
      end
    end)() end

    @testset "AC island detection and independent solving" begin (function ()
      function two_island_net(; second_ref::Symbol = :slack)
        net = Net(name = "two_island", baseMVA = 100.0)
        for name in ("A1", "A2", "B1", "B2")
          addBus!(net = net, busName = name, vn_kV = 110.0)
        end
        addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
        addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "A1")
        addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
        if second_ref === :slack
          addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
        elseif second_ref === :pv
          addProsumer!(net = net, busName = "B1", type = "SYNCHRONOUSMACHINE", p = 12.0, q = 0.0, vm_pu = 1.0, isRegulated = true)
        elseif second_ref === :gen
          # a unit without voltage control: neither slack nor PV
          addProsumer!(net = net, busName = "B1", type = "SYNCHRONOUSMACHINE", p = 12.0, q = 1.0, isRegulated = false)
        end
        addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 8.0, q = 2.0)
        refreshBusTypesFromProsumers!(net)
        return net
      end

      net = two_island_net()
      report = Sparlectra.detect_ac_islands(net)
      @test length(report.rows) == 2
      @test all(row -> row.n_branch == 1, report.rows)
      default_solved = deepcopy(net)
      default_iterations, default_erg = runpf!(default_solved; config = PowerFlowConfig(max_iter = 40))
      @test default_erg == 0
      @test default_iterations > 0

      mktempdir() do tmpdir
        solved = deepcopy(net)
        profile = Dict{Symbol,Any}(:output_dir => tmpdir)
        iterations, erg = runpf!(solved; config = PowerFlowConfig(max_iter = 40, islands_enabled = true), performance_profile = profile)
        @test erg == 0
        @test isfile(joinpath(tmpdir, "ac_islands.csv"))
        for artifact in ("ac_islands.csv", "ac_island_solver_summary.csv", "ac_island_1_solver.log", "q_limit.log", "matpower_dcline.csv")
          @test !isfile(joinpath(pwd(), artifact))
        end
        @test count(!isempty, split(read(joinpath(tmpdir, "ac_islands.csv"), String), '\n')) == 3
        @test all(node -> isfinite(node._vm_pu) && isfinite(node._va_deg), solved.nodeVec)
        rect_status = Sparlectra.rectangular_pf_status(solved)
        @test rect_status.status == :converged
        @test rect_status.reason == :none
        @test rect_status.reason_text == "All AC islands converged independently."
        @test rect_status.island_wise_all_converged === true
        @test rect_status.post_merge_validation_status == :not_applicable
        @test rect_status.post_merge_mismatch_status == :not_applicable
        @test isnan(rect_status.post_merge_final_mismatch)
        @test !occursin("AC island 1 power-flow solve failed", rect_status.reason_text)
        @test iterations == sum(getproperty(status, :iterations) for status in values(profile[:ac_island_solver_statuses]))
        @test rect_status.iterations == iterations
      end

      # without an output directory (library call, N-1 and scenario
      # workers) no island artifact is written; the former fallback was a
      # fixed file in tempdir() that every concurrent solve overwrote
      tmp_artifact = joinpath(tempdir(), "ac_islands.csv")
      tmp_before = isfile(tmp_artifact) ? mtime(tmp_artifact) : nothing
      _, no_dir_erg = runpf!(deepcopy(net); config = PowerFlowConfig(max_iter = 40, islands_enabled = true))
      @test no_dir_erg == 0
      @test (isfile(tmp_artifact) ? mtime(tmp_artifact) : nothing) == tmp_before

      pv_ref_net = two_island_net(second_ref = :pv)
      pv_report = Sparlectra.detect_ac_islands(pv_ref_net)
      @test pv_report.rows[2].status == "promote_pv_ref"
      _, pv_erg = runpf!(pv_ref_net; config = PowerFlowConfig(max_iter = 40, islands_enabled = true))
      @test pv_erg == 0

      # an island always finds itself a reference: a generating unit without
      # voltage control is promoted whatever auto_slack says (only a
      # load-only island has no candidate and is rejected)
      gen_ref_net = two_island_net(second_ref = :gen)
      gen_report = Sparlectra.detect_ac_islands(gen_ref_net; promote_generators = true)
      @test gen_report.rows[2].status == "promote_generator_ref"
      _, gen_erg = redirect_stdout(devnull) do
        runpf!(gen_ref_net; config = PowerFlowConfig(max_iter = 40, islands_enabled = true, auto_slack = false))
      end
      @test gen_erg == 0
      @test all(node -> isfinite(node._vm_pu) && isfinite(node._va_deg), gen_ref_net.nodeVec)

      no_ref_net = two_island_net(second_ref = :none)
      no_ref_report = Sparlectra.detect_ac_islands(no_ref_net; promote_generators = true)
      @test no_ref_report.rows[2].status == "missing_ref"
      @test_throws ErrorException runpf!(no_ref_net; config = PowerFlowConfig(max_iter = 40, islands_enabled = true))

      # 0.30.2 reference priority on the shipped two-island case: the south
      # island has no slack; the unit with the best stated priority
      # (Fennwiese_110, 1) wins over the largest one, and with the priorities
      # cleared reference_candidate_rank decides (Eschwald_110, the largest),
      # no longer the smallest PV bus index (Dornberg_110)
      prio_file = joinpath(dirname(@__DIR__), "data", "scf", "two_islands_prio.scf.json")
      prio_net = importSCF(prio_file)
      south_ref(n) = only(r.chosen_ref_bus for r in Sparlectra.detect_ac_islands(n).rows if n.busDict["Gruenau_110"] in r.buses)
      @test south_ref(prio_net) == prio_net.busDict["Fennwiese_110"]
      plain_net = deepcopy(prio_net)
      foreach(ps -> ps.referencePriority = 0, plain_net.prosumpsVec)
      @test south_ref(plain_net) == plain_net.busDict["Eschwald_110"]
      _, prio_erg = redirect_stdout(devnull) do
        runpf!(prio_net, 40, 1e-8, 0; islands_enabled = true)
      end
      @test prio_erg == 0
      # the chosen reference holds the island's angle
      @test prio_net.nodeVec[prio_net.busDict["Fennwiese_110"]]._va_deg == 0.0
      # the optional SCF field survives the round trip byte for byte, except
      # the "created_by" stamp: it names the package version that wrote the
      # file, so a version bump alone must not fail the comparison
      unstamped(text) = replace(text, r"\"created_by\": \"Sparlectra [^\"]*\"" => "\"created_by\": \"\"")
      @test unstamped(read(exportSCF(importSCF(prio_file); file = joinpath(mktempdir(), "prio.scf.json"), case_name = "two_islands_prio"), String)) == unstamped(read(prio_file, String))
      @test_throws ArgumentError addProsumer!(net = plain_net, busName = "Gruenau_110", type = "GENERATOR", p = 1.0, referencePriority = -1)
    end)() end

    # A network with a single AC island never enters the independent per-island
    # solve branch above, so it must not fall back to reporting the pre-solve
    # placeholder (iterations=0, stage=pre_solve_validation, island_wise_all_converged=false)
    # once it has actually converged.
    @testset "Single-island run reports real convergence, not pre-solve placeholders" begin (function ()
      single_island_net = Net(name = "single_island", baseMVA = 100.0)
      addBus!(net = single_island_net, busName = "A1", vn_kV = 110.0)
      addBus!(net = single_island_net, busName = "A2", vn_kV = 110.0)
      addPIModelACLine!(net = single_island_net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
      addProsumer!(net = single_island_net, busName = "A1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "A1")
      addProsumer!(net = single_island_net, busName = "A2", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
      refreshBusTypesFromProsumers!(single_island_net)

      mktempdir() do tmpdir
        profile = Dict{Symbol,Any}(:output_dir => tmpdir)
        cfg = SparlectraConfig(powerflow = PowerFlowConfig(islands_enabled = true), output = OutputConfig(logfile_results = :off))
        result = run_sparlectra(; net = single_island_net, config = cfg, performance_profile = profile)
        @test result.outcome == :converged
        @test result.iterations > 0
        @test result.diagnostics.island_wise_all_converged === true

        rect_status = Sparlectra.rectangular_pf_status(single_island_net)
        @test rect_status.island_wise_all_converged === true

        island_log = read(joinpath(tmpdir, "ac_island_1_solver.log"), String)
        @test occursin("iterations: $(result.iterations)", island_log)
        @test occursin("stage: post_solve_validation", island_log)
        @test !occursin("iterations: 0", island_log)
        @test !occursin("stage: pre_solve_validation", island_log)
      end
    end)() end

# Ensures final-limit validation remains robust when q-generation data is partially missing.
    @testset "Final limit validation tolerates missing qgen" begin (function ()
      net = createTest3BusNet()
      net.nodeVec[3]._qƩGen = nothing
      io = IOBuffer()
      res = printFinalLimitValidation(net; q_headroom = 0.20, io = io)
      out = String(take!(io))
      @test haskey(res, :q_violations)
      @test haskey(res, :v_violations)
      @test occursin("Final Q-limit validation", out)
      @test occursin("Qgen [MVAr]", out) || occursin("Q-limits: no violations", out)
      empty!(net.qmin_pu)
      append!(net.qmin_pu, [-Inf, -0.01, -Inf])
      empty!(net.qmax_pu)
      append!(net.qmax_pu, [Inf, 0.01, Inf])
      net.nodeVec[2]._qƩGen = 5.0
      # 5 MVAr at bus 2 against limits of +-1 MVAr
      @test printFinalLimitValidation(net; q_headroom = 0.0, io = IOBuffer()).q_violations >= 1
      io = IOBuffer()
      printFinalLimitValidation(net; q_headroom = 0.20, io = io, converged = false)
      @test occursin("Last-iteration Q-limit diagnostic (NR did not converge; values are not a valid final solution)", String(take!(io)))
    end)() end

  end)() end
end

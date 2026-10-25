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

# file: test/test_legacy_voltage_dependent_control.jl
# purpose: legacy profile: the testsets of test/test_voltage_dependent_control.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_voltage_dependent_control.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_voltage_dependent_control_tests()
  # cut from testset: Voltage dependent P(U)/Q(U) control
    # A Q(U) unit makes the solver rebuild S every iteration; the rebuild
    # used to drop the Q-limit clamp, so the regulating unit ran at its
    # scheduled 3 MVAr instead of its 5 MVAr limit. Q(U) unit on the clamped
    # bus itself and on another bus, and the classic outer loop, which
    # clamps by fixing the unit's qVal. The clamped bus reports its physical
    # generator Q (injection plus load), not the band value at 1.04 pu.
    # Last case: the bus is on lock_pv_to_pq_buses with a schedule of
    # -20 MVAr inside the band; the regulating unit keeps its schedule (it
    # ran at -10 when the lock counted the Q(U) unit with its qVal of 2 MVAr
    # instead of its -8 MVAr at 1.04 pu).
    @testset "Q-limit clamp survives the Q(U) rebuild" begin (function ()
      for (qu_bus, mode, lock, q_reg, q_reg_expected) in (("PV", :active_set, false, 3.0, 5.0), ("Q", :active_set, false, 3.0, 5.0), ("PV", :classic_simultaneous, false, 3.0, 5.0), ("PV", :active_set, true, -20.0, -20.0))
        net = Net(name = "qu_clamp_$(qu_bus)_$(mode)_$(lock)", baseMVA = 100.0)
        addBus!(net = net, busName = "Slack", vn_kV = 110.0)
        addBus!(net = net, busName = "PV", vn_kV = 110.0)
        addBus!(net = net, busName = "Q", vn_kV = 110.0)
        addPIModelACLine!(net = net, fromBus = "Slack", toBus = "PV", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = net, fromBus = "Slack", toBus = "Q", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
        addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
        addProsumer!(net = net, busName = "PV", type = "SYNCHRONOUSMACHINE", p = 20.0, q = q_reg, vm_pu = 1.04, isRegulated = true, qMin = -30.0, qMax = 5.0)
        addProsumer!(net = net, busName = "PV", type = "ENERGYCONSUMER", p = 30.0, q = 20.0)
        qu = QUController(make_characteristic([(0.95, 0.10), (1.05, -0.10)]); qmin_MVAr = -50.0, qmax_MVAr = 50.0, sbase_MVA = 100.0)
        addProsumer!(net = net, busName = qu_bus, type = "GENERATOR", p = 5.0, q = 2.0, qMin = -50.0, qMax = 50.0, qu_controller = qu)

        b = geNetBusIdx(net = net, busName = "PV")
        ite, erg = runpf!(net, 30, 1e-10, 0; method = :rectangular, qlimit_enforcement_mode = mode, lock_pv_to_pq_buses = lock ? [b] : Int[])
        @test erg == 0
        @test haskey(net.qLimitEvents, b)
        V = buildVoltageVector(net)
        qgen_phys = imag(Sparlectra.calc_injections(Sparlectra.createYBUS(net = net), V)[b]) * net.baseMVA + 20.0
        q_qu = qu_bus == "PV" ? first(evaluate_controller(qu, abs(V[b]))) * net.baseMVA : 0.0
        # the regulating unit sits at its limit, or at its schedule when locked
        @test isapprox(qgen_phys - q_qu, q_reg_expected; atol = 1e-6)
        @test isapprox(net.nodeVec[b]._qƩGen, qgen_phys; atol = 1e-6)
      end
    end)() end

    # The node kept the scheduled Q of a Q(U) unit on a PQ bus (and the
    # scheduled P of a P(U) unit) after the solve, while the solve and the
    # report used the characteristic: sp_case118 with buses 4, 8 and 10 as
    # PQ, their units scheduled at 15 MVAr and a characteristic of +30 MVAr
    # at 0.95 pu to -30 MVAr at 1.05 pu showed 15.0 on the node and -1.4229,
    # -5.6674, -10.7799 MVAr in the report. Bus 12 stays PV with a P(U)
    # characteristic of 100 MW at 0.95 pu to 60 MW at 1.05 pu: the node and
    # the report showed its scheduled 122 MW where the unit runs at 78 MW at
    # 1.005 pu. Hand form of the characteristics (linear between the two
    # points): Q = 30 - 600 (U - 0.95) MVAr, P = 100 - 400 (U - 0.95) MW.
    # The generator Q from the solved voltages (injection plus bus load)
    # checks that the solve itself ran on the characteristic. The write-back
    # happens only after convergence: a run stopped at the iteration limit
    # keeps the schedule on the node.
    # Third run: branch 114-93 open (island of buses 89, 93, 95, 113, 116),
    # a closed link 4-5 and an island-wise solve, so bus 4 sits in a
    # two-bus PQ cluster of the merged copy; the original node 4 kept its
    # 15 MVAr while the report showed the characteristic (4.27 MVAr). A
    # fixed 10 MVAr unit at bus 4 and a fixed 50 MW unit at bus 12 get
    # their schedule from distributeBusResults!, the controlled units their
    # characteristic: spread over the schedules, the fixed units got 1.02
    # MVAr and 37.2 MW.
    @testset "Q(U) and P(U) generation written back to the node" begin (function ()
      qu_buses = (4, 8, 10)
      build = function (; island::Bool = false)
        mpc = Sparlectra.MatpowerIO.read_case(joinpath(pkgdir(Sparlectra), "data", "mpower", "sp_case118.m"); legacy_compat = true)
        for r in axes(mpc.bus, 1)
          Int(mpc.bus[r, 1]) in qu_buses && (mpc.bus[r, 2] = 1.0)
        end
        for g in axes(mpc.gen, 1)
          Int(mpc.gen[g, 1]) in qu_buses && (mpc.gen[g, 3] = 15.0)
        end
        if island
          for r in axes(mpc.branch, 1)
            Set(Int.(mpc.branch[r, 1:2])) == Set((93, 114)) && (mpc.branch[r, 11] = 0.0)
          end
        end
        n = redirect_stdout(() -> Sparlectra.createNetFromMatPowerCase(mpc = mpc, enable_pq_gen_controllers = true), devnull)
        for ps in n.prosumpsVec
          b = getPosumerBusIndex(ps)
          if b in qu_buses && has_qu_controller(ps)
            ps.quController = QUController(make_characteristic([(0.95, 0.3), (1.05, -0.3)]); qmin_MVAr = -30.0, qmax_MVAr = 30.0, sbase_MVA = n.baseMVA)
          elseif b == 12 && Sparlectra.isGenerator(ps)
            ps.puController = PUController(make_characteristic([(0.95, 1.0), (1.05, 0.6)]); pmin_pu = 0.0, pmax_pu = 2.0)
          end
        end
        if island
          name = Dict(v => k for (k, v) in n.busDict)
          addLink!(net = n, fromBus = name[4], toBus = name[5])
          addProsumer!(net = n, busName = name[4], type = "GENERATOR", p = 0.0, q = 10.0, qMin = -300.0, qMax = 300.0)
          addProsumer!(net = n, busName = name[12], type = "GENERATOR", p = 50.0, q = 0.0, qMin = 0.0, qMax = 0.0)
        end
        return n
      end
      net_stopped = build()
      _, erg_stopped = runpf!(net_stopped, 1, 1e-8, 0; method = :rectangular)
      @test erg_stopped == 1
      @test all(net_stopped.nodeVec[b]._qƩGen == 15.0 for b in qu_buses)
      @test net_stopped.nodeVec[12]._pƩGen == 122.0
      net = build()
      _, erg = runpf!(net, 30, 1e-8, 0; method = :rectangular)
      @test erg == 0
      report = buildACPFlowReport(net)
      V = buildVoltageVector(net)
      s_inj = Sparlectra.calc_injections(Sparlectra.createYBUS(net = net), V) .* net.baseMVA
      p_load, q_load = Sparlectra._bus_consumption_MVA(net)
      for b in qu_buses
        node = net.nodeVec[b]
        @test getNodeType(node) == Sparlectra.PQ
        @test !haskey(net.qLimitEvents, b)
        q_hand = 30.0 - 600.0 * (node._vm_pu - 0.95)
        row = only(r for r in report.nodes if r.bus == b)
        @test isapprox(node._qƩGen, q_hand; atol = 1e-9)
        @test isapprox(row.q_gen_MVar, q_hand; atol = 1e-9)
        @test isapprox(imag(s_inj[b]) + q_load[b], q_hand; atol = 1e-5)
      end
      node12 = net.nodeVec[12]
      @test getNodeType(node12) == Sparlectra.PV
      p_hand = 100.0 - 400.0 * (node12._vm_pu - 0.95)
      @test isapprox(p_hand, 78.0; atol = 1e-6)
      row12 = only(r for r in report.nodes if r.bus == 12)
      @test isapprox(node12._pƩGen, p_hand; atol = 1e-9)
      @test isapprox(row12.p_gen_MW, p_hand; atol = 1e-9)
      @test isapprox(real(s_inj[12]) + p_load[12], p_hand; atol = 1e-5)

      net_isl = build(island = true)
      _, erg_isl = runpf!(net_isl, 30, 1e-8, 0; method = :rectangular, islands_enabled = true)
      @test erg_isl == 0
      report_isl = buildACPFlowReport(net_isl)
      distributeBusResults!(net_isl)
      units(b) = [ps for ps in net_isl.prosumpsVec if Sparlectra.isGenerator(ps) && getPosumerBusIndex(ps) == b]
      for b in qu_buses
        node = net_isl.nodeVec[b]
        q_hand = 30.0 - 600.0 * (node._vm_pu - 0.95)
        q_fixed = b == 4 ? 10.0 : 0.0
        row = only(r for r in report_isl.nodes if r.bus == b)
        @test isapprox(node._qƩGen, q_hand + q_fixed; atol = 1e-9)
        @test isapprox(row.q_gen_MVar, q_hand + q_fixed; atol = 1e-9)
        for ps in units(b)
          @test isapprox(ps.qRes, has_qu_controller(ps) ? q_hand : 10.0; atol = 1e-9)
        end
      end
      @test net_isl.nodeVec[5]._qƩGen === nothing
      node12 = net_isl.nodeVec[12]
      p_hand = 100.0 - 400.0 * (node12._vm_pu - 0.95)
      @test isapprox(node12._pƩGen, p_hand + 50.0; atol = 1e-9)
      for ps in units(12)
        @test isapprox(ps.pRes, has_pu_controller(ps) ? p_hand : 50.0; atol = 1e-9)
      end

      # Fourth run, one island with a closed coupler: runpf! refused P(U),
      # Q(U) and voltage-dependent shunts on a net with a closed link
      # (before 0.60.0), so the third run had to go island-wise. sp_case60
      # (coupler Neubach_110 to Neubach_110b) with the regulated machines of
      # Birkkamp_110 and Westtal_110 on Q(U) (Qmax at U_set - 1 %, Qmin at
      # U_set + 1 %, band -35 .. 35 MVAr) and P(U) (P0 at U_set - 2 %, 0.8 P0
      # at U_set + 2 %), and a 5 MVAr voltage-dependent shunt at
      # Neubach_110b, against the same net with the coupler resolved by hand
      # (link open, the branches of Neubach_110b on Neubach_110, the shunt
      # there too). Both solve to U 1.00723 and 1.01760 pu at the two
      # machines, Q 9.596 and 8.246 MVAr on the characteristic.
      c60 = function (resolved::Bool; solve::Bool = true, regulated::Bool = false)
        n = redirect_stdout(() -> Sparlectra.importSCF(joinpath(pkgdir(Sparlectra), "data", "scf", "sp_case60.scf.json")), devnull)
        if resolved
          l = only(n.linkVec)
          l.status = 0
          for br in n.branchVec
            br.fromBus == l.toBus && (br.fromBus = l.fromBus)
            br.toBus == l.toBus && (br.toBus = l.fromBus)
          end
          empty!(n.isoNodes)
          markIsolatedBuses!(net = n, log = false)
        end
        for name in ("Birkkamp_110", "Westtal_110")
          b = geNetBusIdx(net = n, busName = name)
          ps = only(p for p in n.prosumpsVec if Sparlectra.isGenerator(p) && getPosumerBusIndex(p) == b && Sparlectra.isRegulating(p))
          us, p0 = ps.vm_pu, ps.pVal
          ps.quController = QUController(make_characteristic([(0.99us, 0.35), (1.01us, -0.35)]); qmin_MVAr = -35.0, qmax_MVAr = 35.0, sbase_MVA = n.baseMVA)
          ps.puController = PUController(make_characteristic([(0.98us, p0 / n.baseMVA), (1.02us, 0.8p0 / n.baseMVA)]); pmin_MW = 0.0, pmax_MW = p0, sbase_MVA = n.baseMVA)
          regulated && continue
          ps.isRegulated = false
          ps.vm_pu = nothing
        end
        refreshBusTypesFromProsumers!(n)
        Sparlectra._buildQLimits!(n)
        addShunt!(net = n, busName = resolved ? "Neubach_110" : "Neubach_110b", pShunt = 0.0, qShunt = -5.0, bus_shunt_model = :voltage_dependent_injection)
        solve || return n
        _, erg = runpf!(n, 30, 1e-10, 0; method = :rectangular)
        @test (resolved, erg) == (resolved, 0)
        return n
      end
      closed = c60(false)
      resolved = c60(true)
      @test length(closed.linkVec) == 1 && closed.linkVec[1].status == 1
      member = closed.linkVec[1].toBus
      @test geNetBusIdx(net = closed, busName = "Neubach_110b") == member
      u(n, b) = n.nodeVec[b]._vm_pu * cis(deg2rad(n.nodeVec[b]._va_deg))
      @test maximum(abs(u(closed, b) - u(resolved, b == member ? closed.linkVec[1].fromBus : b)) for b in eachindex(closed.nodeVec)) < 1e-9
      @test isapprox(closed.shuntVec[end].q_shunt, resolved.shuntVec[end].q_shunt; atol = 1e-9)
      report_c = buildACPFlowReport(closed)
      # generator Q from the voltages of the resolved net (its isolated
      # Neubach_110b is out of the Y-bus): injection plus what the bus draws
      active = [b for b in eachindex(resolved.nodeVec) if !(b in resolved.isoNodes)]
      Y = Sparlectra.createYBUS(net = resolved)
      @test size(Y, 1) == length(active)
      v_res = [u(resolved, b) for b in active]
      s_res = Dict(zip(active, v_res .* conj.(Y * v_res) .* resolved.baseMVA))
      _, q_draw = Sparlectra._bus_consumption_MVA(resolved)
      for (name, us, u_exp, q_exp) in (("Birkkamp_110", 1.01, 1.00723, 9.596), ("Westtal_110", 1.02, 1.01760, 8.246))
        b = geNetBusIdx(net = closed, busName = name)
        node = closed.nodeVec[b]
        q_hand = 35.0 - 70.0 * (node._vm_pu - 0.99us) / (0.02us)
        p_hand = 40.0 - 8.0 * (node._vm_pu - 0.98us) / (0.04us)
        row = only(r for r in report_c.nodes if r.bus == b)
        @test (name, getNodeType(node)) == (name, Sparlectra.PQ)
        @test isapprox(node._vm_pu, u_exp; atol = 1e-5)
        @test isapprox(q_hand, q_exp; atol = 1e-3)
        @test isapprox(imag(s_res[b]) + q_draw[b], q_hand; atol = 1e-6)
        @test isapprox(node._qƩGen, q_hand; atol = 1e-9) && isapprox(row.q_gen_MVar, q_hand; atol = 1e-9)
        @test isapprox(node._pƩGen, p_hand; atol = 1e-9) && isapprox(row.p_gen_MW, p_hand; atol = 1e-9)
      end
      # The same two nets under power_flow.solver = apslf, which applies no
      # Q(U)/P(U) characteristic and no voltage-dependent shunt on any net
      # and says so in one warning. With the coupler closed the APSLF solve
      # refused the run after that warning (before 0.60.0: its own refusal
      # of voltage-dependent injections on merged nets outlived the caller's,
      # which became the warning in 0.40.0). The machines keep their voltage
      # regulation here: with regulation off APSLF solves their static Q
      # schedules, and that point does not converge, closed or resolved
      # (mismatch 0.0563 pu both). Both solve the static set points to the
      # same voltages.
      quiet = OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)
      cfg_ap = SparlectraConfig(powerflow = PowerFlowConfig(solver = :apslf, tol = 1e-10), output = quiet)
      ap = map((false, true)) do r
        n = c60(r; solve = false, regulated = true)
        res = @test_logs (:warn, r"power_flow.solver=apslf runs without controllers: .* 2 Q\(U\)/P\(U\) controller\(s\) and 1 voltage-dependent shunt") match_mode = :any redirect_stdout(() -> run_sparlectra(net = n, config = cfg_ap), devnull)
        @test (r, res.final_converged) == (r, true)
        return n
      end
      @test length(ap[1].linkVec) == 1 && ap[1].linkVec[1].status == 1
      @test maximum(abs(u(ap[1], b) - u(ap[2], b == member ? ap[1].linkVec[1].fromBus : b)) for b in eachindex(ap[1].nodeVec)) < 1e-9
      # The generation the closed net reports at each generator bus outside
      # the cluster equals injection plus draw from the voltages of the
      # resolved net. Before 0.60.0 the APSLF solve of a net with a closed
      # link wrote it on the merged copy only: on plain sp_case60
      # Birkkamp_110 reported 0.0 MVAr against 7.258 MVAr with the coupler
      # resolved, and the generation of the slack Altau_220 differed by
      # 19.99 (MW or MVAr, the larger of the two). At Birkkamp_110 and
      # Westtal_110 the node carries the solved Q (8.228 and 6.118 MVAr; 7.258
      # and 5.711 before APSLF took its injections like Newton, which brings
      # the voltage-dependent shunt of sp_case60 in at its 1 pu value);
      # the report evaluates their Q(U) characteristic at the held voltage
      # there (0 MVAr at U_set), which APSLF never applied, on the resolved
      # net as well, so the report rows are compared for the other buses.
      report_ap = buildACPFlowReport(ap[1])
      active_ap = [b for b in eachindex(ap[2].nodeVec) if !(b in ap[2].isoNodes)]
      Y_ap = Sparlectra.createYBUS(net = ap[2])
      v_ap = [u(ap[2], b) for b in active_ap]
      s_ap = Dict(zip(active_ap, v_ap .* conj.(Y_ap * v_ap) .* ap[2].baseMVA))
      p_draw_ap, q_draw_ap = Sparlectra._bus_consumption_MVA(ap[2]; controllers = false)
      gen_buses = unique(getPosumerBusIndex(p) for p in ap[2].prosumpsVec if Sparlectra.isGenerator(p))
      @test geNetBusIdx(net = ap[1], busName = "Altau_220") in gen_buses
      qu_buses = [geNetBusIdx(net = ap[1], busName = name) for name in ("Birkkamp_110", "Westtal_110")]
      for (b, q_exp) in zip(qu_buses, (8.228, 6.118))
        @test isapprox(ap[1].nodeVec[b]._qƩGen, imag(s_ap[b]) + q_draw_ap[b]; atol = 1e-6)
        @test isapprox(ap[1].nodeVec[b]._qƩGen, q_exp; atol = 1e-3)
      end
      for b in setdiff(gen_buses, (member, ap[1].linkVec[1].fromBus), qu_buses)
        row = only(r for r in report_ap.nodes if r.bus == b)
        @test (b, isapprox(row.p_gen_MW, real(s_ap[b]) + p_draw_ap[b]; atol = 1e-6), isapprox(row.q_gen_MVar, imag(s_ap[b]) + q_draw_ap[b]; atol = 1e-6)) == (b, true, true)
      end
    end)() end

  return true
end

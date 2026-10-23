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

# file: test/test_legacy_upfc_control.jl
# purpose: legacy profile: the testsets of test/test_upfc_control.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_upfc_control.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_upfc_control_tests()
  # cut from testset: UPFC composite controller (#325)
  (function ()
    # two parallel corridors between source A and sink B (the series-reactance
    # loop fixture) plus a PQ machine at M2: the shunt converter of the UPFC on
    # corridor A-M2. Series target on A->M2, shunt target at the load bus B.
    function _build_upfc_net(; load_p = 80.0, load_q = 20.0)
      net = Net(name = "upfc_test_loop", baseMVA = 100.0)
      for b in ("A", "M1", "M2", "B")
        addBus!(net = net, busName = b, vn_kV = 110.0)
      end
      addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = load_p, q = load_q)
      addProsumer!(net = net, busName = "M2", type = "GENERATOR", p = 0.0, q = 0.0)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "M1", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "M1", toBus = "B", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "M2", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "M2", toBus = "B", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
      ok, msg = validate!(net = net)
      ok || error("test net invalid: $msg")
      return net
    end
    _solved_state(net) = ([n._vm_pu for n in net.nodeVec], [n._va_deg for n in net.nodeVec])
      # ---- full model (#326): DC-link-coupled UPFC ---------------------------
      # a meshed corridor with the shunt converter at the SENDING bus I and the
      # UPFC on the line I->J; the parallel path S->L lets the UPFC steer the
      # I->J flow. Diagram:
      #   S(slack) --- I ==[UPFC]== J --- L(load)      S ---------------- L
      #                (shunt at I)                     (parallel path)
      function _build_upfc_mesh(; load_p = 90.0, load_q = 30.0)
        net = Net(name = "upfc_full_mesh", baseMVA = 100.0)
        for b in ("S", "I", "J", "L")
          addBus!(net = net, busName = b, vn_kV = 110.0)
        end
        addProsumer!(net = net, busName = "S", type = "EXTERNALNETWORKINJECTION", referencePri = "S", vm_pu = 1.0, va_deg = 0.0)
        addProsumer!(net = net, busName = "I", type = "GENERATOR", p = 0.0, q = 0.0)   # shunt converter
        addProsumer!(net = net, busName = "L", type = "ENERGYCONSUMER", p = load_p, q = load_q)
        addPIModelACLine!(net = net, fromBus = "S", toBus = "I", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = net, fromBus = "I", toBus = "J", r_pu = 0.02, x_pu = 0.18, b_pu = 0.0, status = 1)  # UPFC line
        addPIModelACLine!(net = net, fromBus = "J", toBus = "L", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = net, fromBus = "S", toBus = "L", r_pu = 0.02, x_pu = 0.16, b_pu = 0.0, status = 1)  # parallel path
        ok, msg = validate!(net = net)
        ok || error("mesh net invalid: $msg")
        return net
      end

      @testset "YAML type upfc and double-apply no-op" begin (function ()
        # one row per model: the YAML registration gives the same solved
        # state as the programmatic twin, and a second apply adds nothing
        # (regression: not a FieldError, not a duplicate-controller error;
        # the idempotency check reads fromBus/toBus plus the upfc_group marker)
        rows = [
          (build = _build_upfc_net,
           register = net -> addUpfcControl!(net; fromBus = "A", toBus = "M2", shunt_bus = "M2", target_bus = "B", target_vm_pu = 0.99, p_target_mw = 35.0, v_inj_max_pu = 0.08, s_max_mva = 40.0, name = "upfc_main"),
           entry = Dict{String,Any}("type" => "upfc", "name" => "upfc_main", "from_bus" => "A", "to_bus" => "M2", "shunt_bus" => "M2", "target_bus" => "B", "target_vm_pu" => 0.99, "p_target_mw" => 35.0, "v_inj_max_pu" => 0.08, "s_max_mva" => 40.0),
           solve = net -> run_control!(net),
           # the quadrature composite registers a series and a shunt controller
           registered = net -> (length(net.machineControls),
                                [(c.name, c.upfc_group) for c in net.machineControls if c isa Sparlectra.SeriesReactanceControl],
                                [c.name for c in net.machineControls if c isa MachineVoltageControl]),
           expected = (2, [("upfc_main_series", "upfc_main")], ["upfc_main_shunt"]),
           # same registration, same arithmetic: Vm and Va to 1e-12
           state = net -> vcat(_solved_state(net)...), atol = 1e-12),
          (build = _build_upfc_mesh,
           register = net -> addUpfcControl!(net; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I", p_target_mw = 38.0, q_target_mvar = 15.0, q_shunt_mvar = 10.0, v_inj_max_pu = 0.60, s_max_mva = 120.0, name = "upfc_full", max_outer_iters = 80),
           entry = Dict{String,Any}("type" => "upfc", "name" => "upfc_full", "model" => "full", "from_bus" => "I", "to_bus" => "J", "shunt_bus" => "I", "p_target_mw" => 38.0, "q_target_mvar" => 15.0, "q_shunt_mvar" => 10.0, "v_inj_max_pu" => 0.60, "s_max_mva" => 120.0),
           solve = net -> run_control!(net; control_config = ControlConfig(max_outer_iterations = 80)),
           # the full model is one multi-actuator controller
           registered = net -> (length(Sparlectra._upfc_full_controllers(net)), [(u.name, u.upfc_group) for u in Sparlectra._upfc_full_controllers(net)], String[]),
           expected = (1, [("upfc_full", "upfc_full")], String[]),
           # the YAML entry has no outer iteration budget of its own: sorted Vm to 1e-9
           state = net -> sort([n._vm_pu for n in net.nodeVec]), atol = 1e-9),
        ]
        for row in rows
          net_prog = row.build()
          row.register(net_prog)
          row.solve(net_prog)
          net_yaml = row.build()
          cfg = ControlConfig(controllers = Any[row.entry])
          @test applyConfiguredControllers!(net_yaml, cfg) == 1
          @test row.registered(net_yaml) == row.expected
          row.solve(net_yaml)
          @test maximum(abs.(row.state(net_prog) .- row.state(net_yaml))) <= row.atol
          @test applyConfiguredControllers!(net_yaml, cfg) == 0
          @test row.registered(net_yaml) == row.expected
        end

        # structural validation: unknown key and missing required key
        @test_throws ArgumentError Sparlectra._validate_controller_entries([Dict{String,Any}("type" => "upfc", "name" => "x", "from_bus" => "A", "to_bus" => "M2", "shunt_bus" => "M2", "target_bus" => "B", "target_vm_pu" => 0.99, "p_target_mw" => 35.0, "v_inj_max_pu" => 0.08, "typo_key" => 1)])
        @test_throws ArgumentError Sparlectra._validate_controller_entries([Dict{String,Any}("type" => "upfc", "name" => "x", "from_bus" => "A", "to_bus" => "M2")])
      end)() end

      @testset "full UPFC: independent P and Q on one line (#326)" begin (function ()
        # the headline capability the #325 quadrature composite cannot deliver:
        # distinct P AND Q targets on the same line reached simultaneously.
        # tight deadbands + a generous outer budget so the assertion is on the
        # fixed point, not on the deadband width
        for (pt, qt, qsh) in ((40.0, 10.0, 0.0), (35.0, 20.0, -10.0), (38.0, 15.0, 10.0))
          net = _build_upfc_mesh()
          r = addUpfcControl!(net; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I",
                              p_target_mw = pt, q_target_mvar = qt, q_shunt_mvar = qsh,
                              v_inj_max_pu = 0.60, s_max_mva = 120.0,
                              deadband_p_mw = 1e-3, deadband_q_mvar = 1e-3, max_outer_iters = 80)
          runpf!(net, 30, 1e-8, 0)
          cres = run_control!(net; control_config = ControlConfig(max_outer_iterations = 80))
          u = r.upfc
          @test cres.status == :converged
          @test u.converged
          # both targets hit simultaneously (independent P and Q)
          @test abs(u.achieved_p_mw - pt) <= 0.05
          @test abs(u.achieved_q_mvar - qt) <= 0.05
          # the series injection carries an active component (the phase-shifter
          # DOF): P_se is nonzero here, unlike the quadrature composite
          @test abs(u.p_se_mw) > 1e-3
          # DC-link balance holds on the re-solved state
          @test u.dc_residual_mw <= 0.05
          # the pair is lossless to the network: the series converter's active
          # power is stamped into the branch as |I_s|^2 Re(z_add), and the shunt
          # converter injects exactly that at the bus. Until 0.32.0 the shunt
          # ABSORBED it as well (P_sh had the sign of the absorbed power), so a
          # lossless device drew 2 |I_s|^2 Re(z_add) from the network.
          brij = getNetBranch(net = net, fromBus = "I", toBus = "J")
          vc(b) = (n = net.nodeVec[geNetBusIdx(net = net, busName = b)]; n._vm_pu * cis(deg2rad(n._va_deg)))
          i_s = (vc("I") - vc("J")) / complex(brij.r_pu, brij.x_pu)
          z_add = complex(brij.r_pu - brij.r_base_pu, brij.x_pu - brij.x_base_pu)
          @test isapprox(u.p_sh_mw, abs2(i_s) * real(z_add) * net.baseMVA; atol = 0.05)
          @test isapprox(u.p_se_mw, -u.p_sh_mw; atol = 1e-9)
          # shunt reactive setpoint delivered (within its coupled rating)
          @test isapprox(u.q_sh_mvar, qsh; atol = 1e-6)
          @test u.upfc_group == "UPFC_I_J"
          @test u.series_phase == :free
          # one controller with the series-voltage actuator on the branch,
          # not a pair of controllers as in the quadrature composite
          @test length(Sparlectra._upfc_full_controllers(net)) == 1
          @test [e.quantity for e in controllableElements(net) if e.actuator == :upfc_series_voltage] == [:branch_pq]
        end
      end)() end

      @testset "full UPFC: quadrature series reduces to the SSSC (#326)" begin (function ()
        # series_phase = :quadrature constrains V_se perpendicular to I_s, so
        # P_se = 0 exactly and the series converter is a pure reactance change.
        # With the shunt inactive (q_shunt = 0, P_se = 0) the full UPFC in
        # quadrature equals a standalone SSSC targeting the same line P.
        netq = _build_upfc_mesh()
        rq = addUpfcControl!(netq; model = :full, series_phase = :quadrature, fromBus = "I", toBus = "J", shunt_bus = "I",
                             p_target_mw = 40.0, q_target_mvar = 0.0, q_shunt_mvar = 0.0,
                             v_inj_max_pu = 0.60, s_max_mva = 40.0, deadband_p_mw = 1e-3, max_outer_iters = 80)
        runpf!(netq, 30, 1e-8, 0)
        run_control!(netq; control_config = ControlConfig(max_outer_iterations = 80))
        @test rq.upfc.p_se_mw == 0.0                 # quadrature: no active injection
        @test abs(rq.upfc.achieved_p_mw - 40.0) <= 0.05

        # standalone SSSC with a wide reactance window (so neither device clamps)
        nets = _build_upfc_mesh()
        addSeriesReactanceControl!(nets; fromBus = "I", toBus = "J", p_target_mw = 40.0, x_min_pu = -0.10, x_max_pu = 0.60, deadband_p_mw = 1e-3)
        runpf!(nets, 30, 1e-8, 0)
        run_control!(nets; control_config = ControlConfig(max_outer_iterations = 80))
        vq = sort([n._vm_pu for n in netq.nodeVec])
        vs = sort([n._vm_pu for n in nets.nodeVec])
        # both reach the same line P via a reactance change, so the solved states
        # agree (the small residual is the two secants' deadband slop)
        @test maximum(abs.(vq .- vs)) <= 5e-3
        @test isapprox(get_branch_p_from_to_mw(nets, "I", "J"), rq.upfc.achieved_p_mw; atol = 0.1)
      end)() end

      @testset "full UPFC: classical result reports the controller (#326)" begin (function ()
        # the classical printACPFlowResults must count the UPFC and print its
        # summary block (the FACTS controllers used to be invisible there)
        net = _build_upfc_mesh()
        addUpfcControl!(net; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I",
                        p_target_mw = 40.0, q_target_mvar = 10.0, q_shunt_mvar = 0.0,
                        v_inj_max_pu = 0.30, s_max_mva = 120.0, deadband_p_mw = 1e-2, deadband_q_mvar = 1e-2, max_outer_iters = 80)
        run_control!(net; control_config = ControlConfig(max_outer_iterations = 80))
        calcNetLosses!(net)
        txt = mktempdir() do d
          cd(d) do
            printACPFlowResults(net, 0.0, 1, 1e-8, true)
            read("result_$(net.name).txt", String)
          end
        end
        @test occursin("UPFC: 1", txt)                              # counted on the Controllers line
        @test occursin("UPFC Control Summary", txt)                 # engineering summary block
        @test occursin("DC-link residual", txt)                     # the coupling quantity is reported
        @test occursin("line P target/achieved", txt)
        # the controlled branch row is marked (Ctrl = UPFC, its P target, status)
        @test any(l -> occursin("I -> J", l) && occursin("UPFC", l) && occursin("40.000", l), split(txt, "\n"))
        # the standalone summary printer works too
        sub = sprint(printUpfcFullControllerSummary, net)
        @test occursin("UPFC_I_J", sub)
        # the quadrature composite (SSSC + STATCOM) is now counted honestly too:
        # the machine (STATCOM) side used to be uncounted on the Controllers line
        netq = _build_upfc_net()
        addUpfcControl!(netq; fromBus = "A", toBus = "M2", shunt_bus = "M2", target_bus = "B",
                        target_vm_pu = 0.99, p_target_mw = 35.0, v_inj_max_pu = 0.08, s_max_mva = 40.0)
        run_control!(netq)
        calcNetLosses!(netq)
        qtxt = mktempdir() do d
          cd(d) do
            printACPFlowResults(netq, 0.0, 1, 1e-8, true)
            read("result_$(netq.name).txt", String)
          end
        end
        @test occursin("MachV: 1", qtxt)   # the STATCOM shunt side
        @test occursin("TCSC: 1", qtxt)    # the SSSC series side
      end)() end

      @testset "full UPFC: low-current guard keeps z_add finite (#326)" begin (function ()
        # z_add = V_se / I_s is ill-conditioned as the line current vanishes; the
        # `_UPFC_MIN_CURRENT_PU` floor must keep every result finite (no NaN/Inf)
        # and leave the branch at its base impedance on a (near) dead line.
        function _chain(load_mw)
          m = Net(name = "upfc_chain", baseMVA = 100.0)
          for b in ("S", "I", "J", "L")
            addBus!(net = m, busName = b, vn_kV = 110.0)
          end
          addProsumer!(net = m, busName = "S", type = "EXTERNALNETWORKINJECTION", referencePri = "S", vm_pu = 1.0, va_deg = 0.0)
          addProsumer!(net = m, busName = "I", type = "GENERATOR", p = 0.0, q = 0.0)
          addProsumer!(net = m, busName = "L", type = "ENERGYCONSUMER", p = load_mw, q = load_mw * 0.3)
          addPIModelACLine!(net = m, fromBus = "S", toBus = "I", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
          addPIModelACLine!(net = m, fromBus = "I", toBus = "J", r_pu = 0.02, x_pu = 0.18, b_pu = 0.0, status = 1)
          addPIModelACLine!(net = m, fromBus = "J", toBus = "L", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
          ok, msg = validate!(net = m)
          ok || error("chain net invalid: $msg")
          return m
        end
        for load in (5.0, 0.05, 0.0005)
          net = _chain(load)
          r = addUpfcControl!(net; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I",
                              p_target_mw = 0.0, q_target_mvar = 0.0, q_shunt_mvar = 0.0,
                              v_inj_max_pu = 0.30, s_max_mva = 120.0,
                              deadband_p_mw = 1e-3, deadband_q_mvar = 1e-3, max_outer_iters = 60)
          run_control!(net; control_config = ControlConfig(max_outer_iterations = 60))
          u = r.upfc
          brij = getNetBranch(net = net, fromBus = "I", toBus = "J")
          @test all(isfinite, (real(u.v_se_pu), imag(u.v_se_pu), brij.r_pu, brij.x_pu, u.p_se_mw))
          @test isfinite(u.dc_residual_mw)
        end
        # the smallest load drives |I_s| below the 1e-4 pu floor; the guard then
        # holds the branch at its base impedance (no blow-up)
        net = _chain(0.0005)
        r = addUpfcControl!(net; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I",
                            p_target_mw = 0.0, q_target_mvar = 0.0, q_shunt_mvar = 0.0,
                            v_inj_max_pu = 0.30, s_max_mva = 120.0, deadband_p_mw = 1e-3, deadband_q_mvar = 1e-3, max_outer_iters = 60)
        run_control!(net; control_config = ControlConfig(max_outer_iterations = 60))
        brij = getNetBranch(net = net, fromBus = "I", toBus = "J")
        @test abs(r.upfc.i_s_pu) < 1e-4                 # below the guard floor
        @test isapprox(brij.r_pu, 0.02; atol = 1e-9)    # base impedance kept
        @test isapprox(brij.x_pu, 0.18; atol = 1e-9)
      end)() end

      @testset "full UPFC: SC and export use the physical base impedance (#329)" begin (function ()
        # The full model stamps the series converter's active injection onto the
        # LIVE branch impedance (its resistance part can go negative). Short
        # circuit and the interchange exports must read the PHYSICAL BASE impedance
        # (r_base_pu/x_base_pu) instead, so a UPFC control run no longer
        # contaminates fault or export data. The base net (no UPFC) and the
        # compensated net must then agree.
        _feeder(bus) = (mrid = bus, name = "F_" * bus, bus = bus, maxInitialSymShCCurrent_A = 20_000.0, minInitialSymShCCurrent_A = 16_000.0, maxR1ToX1Ratio = 0.1, minR1ToX1Ratio = 0.1, maxR0ToX0Ratio = nothing, maxZ0ToZ1Ratio = nothing, ikSecond = nothing, governorSCD = nothing)
        scd = Sparlectra.CGMESImporter.CGMESShortCircuitData([_feeder("S")], NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[])
        # base network (no UPFC): the short circuit reference result
        base = _build_upfc_mesh()
        rb = runShortCircuit!(base, scd; case = :max)
        @test length(rb.rows) >= 1
        base_brij = getNetBranch(net = base, fromBus = "I", toBus = "J")

        # after a full-UPFC control run the I->J branch carries a live r that
        # differs from its physical base (the sign of Re(z_add) depends on the
        # operating point: until 0.32.1 it was negative here as a consequence of
        # the inverted DC-link balance), while the base is preserved
        net = _build_upfc_mesh()
        addUpfcControl!(net; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I",
                        p_target_mw = 48.0, q_target_mvar = 8.0, q_shunt_mvar = 0.0,
                        v_inj_max_pu = 0.30, s_max_mva = 120.0, deadband_p_mw = 1e-2, deadband_q_mvar = 1e-2, max_outer_iters = 80)
        run_control!(net; control_config = ControlConfig(max_outer_iterations = 80))
        brij = getNetBranch(net = net, fromBus = "I", toBus = "J")
        @test !isapprox(brij.r_pu, brij.r_base_pu; atol = 1e-9)  # live operating point differs from the base
        @test brij.r_base_pu > 0.0                               # physical base preserved
        @test isapprox(brij.r_base_pu, base_brij.r_pu; atol = 1e-12)
        @test isapprox(brij.x_base_pu, base_brij.x_pu; atol = 1e-12)

        # short circuit now SUCCEEDS on the compensated net (it reads the base) and
        # matches the base-net fault currents bus for bus
        ru = runShortCircuit!(net, scd; case = :max)
        @test length(ru.rows) == length(rb.rows)
        for (a, b) in zip(sort(rb.rows; by = r -> r.bus_idx), sort(ru.rows; by = r -> r.bus_idx))
          @test a.bus_idx == b.bus_idx
          @test (isnan(a.ik_kA) && isnan(b.ik_kA)) || isapprox(a.ik_kA, b.ik_kA; atol = 1e-9)
        end

        # both interchange exports now SUCCEED on the compensated net, writing the
        # base impedance instead of rejecting the negative-r operating point
        mktempdir() do d
          base_ok = _build_upfc_mesh()
          runpf!(base_ok, 30, 1e-8, 0)
          calcNetLosses!(base_ok)
          writeCGMESFiles(base_ok; path = mkpath(joinpath(d, "base_cgmes")))
          writeMatpowerCasefile(base_ok, joinpath(d, "base.m"))
          calcNetLosses!(net)
          @test (writeCGMESFiles(net; path = mkpath(joinpath(d, "upfc_cgmes"))); true)
          @test (writeMatpowerCasefile(net, joinpath(d, "upfc.m")); true)
        end

        # restoreBaseImpedances! and clearUpfcFullControllers! return the live
        # field to the physical base (#329)
        restoreBaseImpedances!(net)
        @test isapprox(brij.r_pu, brij.r_base_pu; atol = 1e-12)
        @test isapprox(brij.x_pu, brij.x_base_pu; atol = 1e-12)
        clearUpfcFullControllers!(net)
        @test isempty([c for c in net.machineControls if c isa Sparlectra.UpfcFullControl])
      end)() end

  end)()
  return true
end

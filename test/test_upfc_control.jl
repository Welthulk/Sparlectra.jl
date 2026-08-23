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

# file: test/test_upfc_control.jl
# purpose: tests for the UPFC stationary quadrature composite (issue #325):
#          exact equivalence with the manually registered SSSC+STATCOM pair,
#          both limit characteristics at their clamps, all-or-nothing
#          registration, result-row pairing, and the YAML type upfc with the
#          double-apply no-op.

function run_upfc_control_tests()
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

  @testset "UPFC composite controller (#325)" begin
    @testset "composite equals the manual SSSC+STATCOM pair" begin
      # net A: the two controllers registered by hand, in the same order the
      # composite uses (series first: run_control! executes in registration
      # order, so the order is part of the identity)
      net_manual = _build_upfc_net()
      s_manual = addSeriesReactanceControl!(net_manual; fromBus = "A", toBus = "M2", p_target_mw = 35.0, v_inj_max_pu = 0.08, name = "pair_series")
      addMachineVoltageControl!(net_manual; bus = "M2", target_bus = "B", target_vm_pu = 0.99, s_max_mva = 40.0, name = "pair_shunt")
      m_manual = last(net_manual.machineControls)
      r_manual = run_sparlectra(net = net_manual)

      # net B: one composite call with the identical parameters
      net_upfc = _build_upfc_net()
      upfc = addUpfcControl!(net_upfc; fromBus = "A", toBus = "M2", shunt_bus = "M2", target_bus = "B", target_vm_pu = 0.99, p_target_mw = 35.0, v_inj_max_pu = 0.08, s_max_mva = 40.0)
      r_upfc = run_sparlectra(net = net_upfc)

      @test r_manual.final_converged
      @test r_upfc.final_converged
      # identical code paths must give identical numbers, an exactness
      # check, not a tolerance check (1e-12 per the acceptance criterion)
      vm_a, va_a = _solved_state(net_manual)
      vm_b, va_b = _solved_state(net_upfc)
      @test maximum(abs.(vm_a .- vm_b)) <= 1e-12
      @test maximum(abs.(va_a .- va_b)) <= 1e-12
      @test abs(upfc.series.x_pu - s_manual.x_pu) <= 1e-12
      @test abs(upfc.series.achieved_p_mw - s_manual.achieved_p_mw) <= 1e-12
      @test abs(upfc.shunt.q_mvar - m_manual.q_mvar) <= 1e-12
      @test abs(upfc.shunt.achieved_vm_pu - m_manual.achieved_vm_pu) <= 1e-12
      @test upfc.series.converged == s_manual.converged
      @test upfc.shunt.converged == m_manual.converged
      @test upfc.series.at_limit == s_manual.at_limit
      @test upfc.shunt.at_limit == m_manual.at_limit

      # composite identity: names, group marker, return value
      @test upfc.name == "UPFC_A_M2"
      @test upfc.series.name == "UPFC_A_M2_series"
      @test upfc.shunt.name == "UPFC_A_M2_shunt"
      @test upfc.series.upfc_group == "UPFC_A_M2"
      @test upfc.shunt.upfc_group == "UPFC_A_M2"
      @test upfc.series isa Sparlectra.SeriesReactanceControl
      @test upfc.shunt isa MachineVoltageControl

      # result rows: one row per actuator, device strings name the pair
      rows = controllableElements(net_upfc)
      @test length(rows) == 2
      series_row = only([r for r in rows if r.actuator == :series_x_pu])
      shunt_row = only([r for r in rows if r.actuator == :machine_q_mvar])
      @test series_row.device == "UPFC series (VSC pair, stationary quadrature model)"
      @test shunt_row.device == "UPFC shunt (VSC pair, stationary quadrature model)"
      # the manual pair keeps the plain device vocabulary
      manual_rows = controllableElements(net_manual)
      @test only([r for r in manual_rows if r.actuator == :series_x_pu]).device == "SSSC (VSC)"
      @test only([r for r in manual_rows if r.actuator == :machine_q_mvar]).device == "STATCOM (VSC)"
    end

    @testset "both converter limits at their clamps" begin
      # shunt side driven into its rating: a small converter cannot lift the
      # sagging load bus to 1.02 pu; delivered Q tracks V * S_max (the
      # STATCOM acceptance, live bound from the solved terminal voltage)
      net = _build_upfc_net()
      upfc = addUpfcControl!(net; fromBus = "A", toBus = "M2", shunt_bus = "M2", target_bus = "B", target_vm_pu = 1.02, p_target_mw = 35.0, v_inj_max_pu = 0.08, s_max_mva = 12.0)
      run_control!(net)
      vt = get_bus_vm_pu(net, "M2")
      @test !upfc.shunt.converged
      @test upfc.shunt.at_limit
      @test upfc.shunt.q_mvar ≈ vt * 12.0 atol = 2e-2
      @test upfc.shunt.qmax_mvar ≈ vt * 12.0 atol = 1e-9

      # series side driven into its rating: a tight injectable voltage pins
      # the reactance window; the EFFECTIVE injected voltage
      # |x - x_base| * |I| sits at v_inj_max (the SSSC acceptance)
      net2 = _build_upfc_net()
      upfc2 = addUpfcControl!(net2; fromBus = "A", toBus = "M2", shunt_bus = "M2", target_bus = "B", target_vm_pu = 0.99, p_target_mw = 35.0, v_inj_max_pu = 0.01, s_max_mva = 40.0)
      run_control!(net2)
      @test !upfc2.series.converged
      @test upfc2.series.at_limit
      v_inj = abs(upfc2.series.x_pu - upfc2.series.x_base_pu) * upfc2.series.i_pu
      @test v_inj ≈ 0.01 rtol = 5e-2
      p_lim = get_branch_p_from_to_mw(net2, "A", "M2")
      @test p_lim < 35.0 - upfc2.series.deadband_p_mw
    end

    @testset "composite-level rejections leave the net untouched" begin
      # every rejected call must leave zero controllers behind
      # (all-or-nothing rule)
      rejections = [
        # shunt_bus not an end of the corridor
        () -> (net = _build_upfc_net(); (net, () -> addUpfcControl!(net; fromBus = "A", toBus = "M2", shunt_bus = "M1", target_bus = "B", target_vm_pu = 0.99, p_target_mw = 35.0, v_inj_max_pu = 0.08, s_max_mva = 40.0))),
        # both ratings given
        () -> (net = _build_upfc_net(); (net, () -> addUpfcControl!(net; fromBus = "A", toBus = "M2", shunt_bus = "M2", target_bus = "B", target_vm_pu = 0.99, p_target_mw = 35.0, v_inj_max_pu = 0.08, s_max_mva = 40.0, i_max_ka = 0.2))),
        # neither rating given
        () -> (net = _build_upfc_net(); (net, () -> addUpfcControl!(net; fromBus = "A", toBus = "M2", shunt_bus = "M2", target_bus = "B", target_vm_pu = 0.99, p_target_mw = 35.0, v_inj_max_pu = 0.08))),
        # non-positive injectable voltage
        () -> (net = _build_upfc_net(); (net, () -> addUpfcControl!(net; fromBus = "A", toBus = "M2", shunt_bus = "M2", target_bus = "B", target_vm_pu = 0.99, p_target_mw = 35.0, v_inj_max_pu = 0.0, s_max_mva = 40.0))),
        # no branch between the named buses
        () -> (net = _build_upfc_net(); (net, () -> addUpfcControl!(net; fromBus = "A", toBus = "B", shunt_bus = "A", target_bus = "M1", target_vm_pu = 0.99, p_target_mw = 35.0, v_inj_max_pu = 0.08, s_max_mva = 40.0))),
      ]
      for build in rejections
        net, call = build()
        @test_throws ErrorException call()
        @test length(net.machineControls) == 0
      end

      # rollback: the series side registers, then the shunt side fails
      # (target equals the machine bus, rejected by addMachineVoltageControl!);
      # the composite must remove the series controller again and rethrow
      net = _build_upfc_net()
      @test_throws ErrorException addUpfcControl!(net; fromBus = "A", toBus = "M2", shunt_bus = "M2", target_bus = "M2", target_vm_pu = 0.99, p_target_mw = 35.0, v_inj_max_pu = 0.08, s_max_mva = 40.0)
      @test length(net.machineControls) == 0
    end

    @testset "YAML type upfc and double-apply no-op" begin
      net_prog = _build_upfc_net()
      addUpfcControl!(net_prog; fromBus = "A", toBus = "M2", shunt_bus = "M2", target_bus = "B", target_vm_pu = 0.99, p_target_mw = 35.0, v_inj_max_pu = 0.08, s_max_mva = 40.0, name = "upfc_main")
      run_control!(net_prog)

      net_yaml = _build_upfc_net()
      cfg = ControlConfig(controllers = Any[Dict{String,Any}("type" => "upfc", "name" => "upfc_main", "from_bus" => "A", "to_bus" => "M2", "shunt_bus" => "M2", "target_bus" => "B", "target_vm_pu" => 0.99, "p_target_mw" => 35.0, "v_inj_max_pu" => 0.08, "s_max_mva" => 40.0)])
      @test applyConfiguredControllers!(net_yaml, cfg) == 1
      @test length(net_yaml.machineControls) == 2
      series = only([c for c in net_yaml.machineControls if c isa Sparlectra.SeriesReactanceControl])
      shunt = only([c for c in net_yaml.machineControls if c isa MachineVoltageControl])
      @test series.name == "upfc_main_series"
      @test shunt.name == "upfc_main_shunt"
      @test series.upfc_group == "upfc_main"
      run_control!(net_yaml)
      # same registration -> same numbers as the programmatic twin
      vm_a, va_a = _solved_state(net_prog)
      vm_b, va_b = _solved_state(net_yaml)
      @test maximum(abs.(vm_a .- vm_b)) <= 1e-12
      @test maximum(abs.(va_a .- va_b)) <= 1e-12

      # double apply is a no-op (regression: not a FieldError, not a
      # duplicate-controller error; the idempotency check reads
      # fromBus/toBus plus the upfc_group marker)
      @test applyConfiguredControllers!(net_yaml, cfg) == 0
      @test length(net_yaml.machineControls) == 2

      # structural validation: unknown key and missing required key
      @test_throws ArgumentError Sparlectra._validate_controller_entries([Dict{String,Any}("type" => "upfc", "name" => "x", "from_bus" => "A", "to_bus" => "M2", "shunt_bus" => "M2", "target_bus" => "B", "target_vm_pu" => 0.99, "p_target_mw" => 35.0, "v_inj_max_pu" => 0.08, "typo_key" => 1)])
      @test_throws ArgumentError Sparlectra._validate_controller_entries([Dict{String,Any}("type" => "upfc", "name" => "x", "from_bus" => "A", "to_bus" => "M2")])
    end

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

    @testset "full UPFC: independent P and Q on one line (#326)" begin
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
        # shunt reactive setpoint delivered (within its coupled rating)
        @test isapprox(u.q_sh_mvar, qsh; atol = 1e-6)
        @test u.upfc_group == "UPFC_I_J"
        @test u.series_phase == :free
      end
    end

    @testset "full UPFC: quadrature series reduces to the SSSC (#326)" begin
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
    end

    @testset "full UPFC: result rows and element view (#326)" begin
      net = _build_upfc_mesh()
      r = addUpfcControl!(net; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I",
                          p_target_mw = 38.0, q_target_mvar = 15.0, q_shunt_mvar = 10.0,
                          v_inj_max_pu = 0.60, s_max_mva = 120.0, deadband_p_mw = 1e-3, deadband_q_mvar = 1e-3, max_outer_iters = 80)
      runpf!(net, 30, 1e-8, 0)
      run_control!(net; control_config = ControlConfig(max_outer_iterations = 80))
      els = [e for e in controllableElements(net) if e.actuator == :upfc_series_voltage]
      @test length(els) == 1
      @test els[1].device == "UPFC (full, DC-link coupled)"
      @test els[1].element == "branch@I-J"
      @test els[1].quantity == :branch_pq
      # one controller (not a pair): full model is a single multi-actuator device
      @test length(Sparlectra._upfc_full_controllers(net)) == 1
    end

    @testset "full UPFC: classical result reports the controller (#326)" begin
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
      # base net without the UPFC keeps the byte-stable "none" anchor
      base = _build_upfc_mesh()
      btxt = mktempdir() do d
        cd(d) do
          runpf!(base, 30, 1e-8, 0)
          calcNetLosses!(base)
          printACPFlowResults(base, 0.0, 1, 1e-8, true)
          read("result_$(base.name).txt", String)
        end
      end
      @test occursin("Transformer controls: none", btxt)
      @test !occursin("UPFC Control Summary", btxt)

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
    end

    @testset "full UPFC: low-current guard keeps z_add finite (#326)" begin
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
    end

    @testset "full UPFC: negative-r branch is rejected by short circuit (#326)" begin
      # the full model maps the series converter's active injection to a
      # NEGATIVE branch resistance, persisted in place. A subsequent IEC 60909
      # short circuit must FAIL LOUDLY on that non-physical branch instead of
      # silently computing wrong fault currents (the converter is bypassed
      # under fault; the SC needs the physical impedance).
      _feeder(bus) = (mrid = bus, name = "F_" * bus, bus = bus, maxInitialSymShCCurrent_A = 20_000.0, minInitialSymShCCurrent_A = 16_000.0, maxR1ToX1Ratio = 0.1, minR1ToX1Ratio = 0.1, maxR0ToX0Ratio = nothing, maxZ0ToZ1Ratio = nothing, ikSecond = nothing, governorSCD = nothing)
      scd = Sparlectra.CGMESImporter.CGMESShortCircuitData([_feeder("S")], NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[])
      # base network (no UPFC): the short circuit runs normally
      base = _build_upfc_mesh()
      rb = runShortCircuit!(base, scd; case = :max)
      @test length(rb.rows) >= 1
      # after a full-UPFC control run the I->J branch carries a negative r
      net = _build_upfc_mesh()
      addUpfcControl!(net; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I",
                      p_target_mw = 48.0, q_target_mvar = 8.0, q_shunt_mvar = 0.0,
                      v_inj_max_pu = 0.30, s_max_mva = 120.0, deadband_p_mw = 1e-2, deadband_q_mvar = 1e-2, max_outer_iters = 80)
      run_control!(net; control_config = ControlConfig(max_outer_iterations = 80))
      @test getNetBranch(net = net, fromBus = "I", toBus = "J").r_pu < 0.0
      err = try
        runShortCircuit!(net, scd; case = :max)
        nothing
      catch e
        e
      end
      @test err !== nothing
      @test occursin("negative series resistance", sprint(showerror, err))

      # symmetric with the SC guard: neither interchange export may write the
      # compensated (negative-r) operating point silently
      mktempdir() do d
        base_ok = _build_upfc_mesh()
        runpf!(base_ok, 30, 1e-8, 0)
        calcNetLosses!(base_ok)
        bdir = mkpath(joinpath(d, "base_cgmes"))
        writeCGMESFiles(base_ok; path = bdir)                 # base net exports fine
        writeMatpowerCasefile(base_ok, joinpath(d, "base.m"))
        udir = mkpath(joinpath(d, "upfc_cgmes"))
        cg = try
          writeCGMESFiles(net; path = udir)
          nothing
        catch e
          e
        end
        @test cg !== nothing && occursin("negative series resistance", sprint(showerror, cg))
        mp = try
          writeMatpowerCasefile(net, joinpath(d, "upfc.m"))
          nothing
        catch e
          e
        end
        @test mp !== nothing && occursin("negative series resistance", sprint(showerror, mp))
      end
    end

    @testset "full UPFC: registration validation (#326)" begin
      # model = :full needs q_target_mvar and rejects the quadrature-only keys
      net = _build_upfc_mesh()
      @test_throws ErrorException addUpfcControl!(net; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I", p_target_mw = 40.0, v_inj_max_pu = 0.6, s_max_mva = 120.0)  # no q_target_mvar
      @test length(net.machineControls) == 0
      @test_throws ErrorException addUpfcControl!(net; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I", p_target_mw = 40.0, q_target_mvar = 10.0, v_inj_max_pu = 0.6, s_max_mva = 120.0, target_bus = "J", target_vm_pu = 1.0)  # :full rejects target_bus/target_vm_pu
      @test length(net.machineControls) == 0
      # both / neither shunt rating
      @test_throws ErrorException addUpfcControl!(net; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I", p_target_mw = 40.0, q_target_mvar = 10.0, v_inj_max_pu = 0.6, s_max_mva = 120.0, i_max_ka = 0.5)
      @test_throws ErrorException addUpfcControl!(net; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I", p_target_mw = 40.0, q_target_mvar = 10.0, v_inj_max_pu = 0.6)
      # non-positive injected-voltage limit; no branch; no generator at shunt bus
      @test_throws ErrorException addUpfcControl!(net; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I", p_target_mw = 40.0, q_target_mvar = 10.0, v_inj_max_pu = 0.0, s_max_mva = 120.0)
      @test_throws ErrorException addUpfcControl!(net; model = :full, fromBus = "I", toBus = "S", shunt_bus = "I", p_target_mw = 40.0, q_target_mvar = 10.0, v_inj_max_pu = 0.6, s_max_mva = 120.0)  # no I-S? (there is S-I; wrong order) -> orientation error
      @test_throws ErrorException addUpfcControl!(net; model = :full, fromBus = "I", toBus = "J", shunt_bus = "J", p_target_mw = 40.0, q_target_mvar = 10.0, v_inj_max_pu = 0.6, s_max_mva = 120.0)  # no generator at J
      @test length(net.machineControls) == 0
      # bad model symbol
      @test_throws ErrorException addUpfcControl!(net; model = :bogus, fromBus = "I", toBus = "J", shunt_bus = "I", p_target_mw = 40.0, v_inj_max_pu = 0.6, s_max_mva = 120.0)
    end

    @testset "full UPFC: YAML type upfc with model=full (#326)" begin
      net_prog = _build_upfc_mesh()
      addUpfcControl!(net_prog; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I",
                      p_target_mw = 38.0, q_target_mvar = 15.0, q_shunt_mvar = 10.0,
                      v_inj_max_pu = 0.60, s_max_mva = 120.0, name = "upfc_full", max_outer_iters = 80)
      run_control!(net_prog; control_config = ControlConfig(max_outer_iterations = 80))

      net_yaml = _build_upfc_mesh()
      cfg = ControlConfig(controllers = Any[Dict{String,Any}("type" => "upfc", "name" => "upfc_full", "model" => "full",
        "from_bus" => "I", "to_bus" => "J", "shunt_bus" => "I", "p_target_mw" => 38.0, "q_target_mvar" => 15.0,
        "q_shunt_mvar" => 10.0, "v_inj_max_pu" => 0.60, "s_max_mva" => 120.0)])
      @test applyConfiguredControllers!(net_yaml, cfg) == 1
      @test length(Sparlectra._upfc_full_controllers(net_yaml)) == 1
      u = only(Sparlectra._upfc_full_controllers(net_yaml))
      @test u.name == "upfc_full"
      @test u.upfc_group == "upfc_full"
      run_control!(net_yaml; control_config = ControlConfig(max_outer_iterations = 80))
      vm_a = sort([n._vm_pu for n in net_prog.nodeVec])
      vm_b = sort([n._vm_pu for n in net_yaml.nodeVec])
      @test maximum(abs.(vm_a .- vm_b)) <= 1e-9
      # double apply is a no-op
      @test applyConfiguredControllers!(net_yaml, cfg) == 0
      @test length(Sparlectra._upfc_full_controllers(net_yaml)) == 1
    end
  end
end

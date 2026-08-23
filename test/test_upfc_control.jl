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
  end
end

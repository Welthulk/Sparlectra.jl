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

# file: test/test_series_reactance_control.jl
# purpose: tests for the TCSC-like series-reactance controller (issue #297):
#          loop-network target tracking, honest at_limit, bit-identical
#          baseline without the controller, element-row vocabulary,
#          registration validation, and deadband behavior.

function run_series_reactance_control_tests()
  # two parallel corridors between source A and sink B; the corridor
  # reactance ratio 1:2 puts one third of the transfer on the controlled
  # corridor at baseline, so a higher target needs a visible reactance move
  function _build_loop_net()
    net = Net(name = "tcsc_test_loop", baseMVA = 100.0)
    for b in ("A", "M1", "M2", "B")
      addBus!(net = net, busName = b, vn_kV = 110.0)
    end
    addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
    addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 80.0, q = 20.0)
    addPIModelACLine!(net = net, fromBus = "A", toBus = "M1", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "M1", toBus = "B", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "A", toBus = "M2", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "M2", toBus = "B", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
    ok, msg = validate!(net = net)
    ok || error("test net invalid: $msg")
    return net
  end

  # net with a transformer branch for the rejection test
  function _build_trafo_net()
    net = Net(name = "tcsc_test_trafo", baseMVA = 100.0)
    for b in ("A", "T", "B")
      addBus!(net = net, busName = b, vn_kV = 110.0)
    end
    addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
    addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 30.0, q = 10.0)
    addPIModelTrafo!(net = net, fromBus = "A", toBus = "T", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "T", toBus = "B", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
    return net
  end

  _solved_state(net) = ([n._vm_pu for n in net.nodeVec], [n._va_deg for n in net.nodeVec])

  @testset "Series reactance controller (#297)" begin
    @testset "loop network: target reached within deadband" begin
      net = _build_loop_net()
      ctrl = addSeriesReactanceControl!(net; fromBus = "A", toBus = "M2", p_target_mw = 35.0, x_min_pu = 0.02, x_max_pu = 0.30)
      result = run_sparlectra(net = net)
      @test result.final_converged
      p = get_branch_p_from_to_mw(net, "A", "M2")
      @test abs(p - 35.0) <= ctrl.deadband_p_mw
      @test ctrl.converged
      @test !ctrl.at_limit
      @test ctrl.status == :converged
      @test ctrl.x_min_pu <= ctrl.x_pu <= ctrl.x_max_pu
      # the actuator moved: baseline reactance was 0.20
      @test ctrl.x_pu < 0.20
      br = getNetBranch(net = net, fromBus = "A", toBus = "M2")
      @test isapprox(br.x_pu, ctrl.x_pu; atol = 1e-12)
      @test ctrl.achieved_p_mw !== nothing
    end

    @testset "out-of-range target: honest at_limit" begin
      net = _build_loop_net()
      ctrl = addSeriesReactanceControl!(net; fromBus = "A", toBus = "M2", p_target_mw = 70.0, x_min_pu = 0.02, x_max_pu = 0.30)
      result = run_sparlectra(net = net)
      # the solve itself succeeds, only the controller reports the miss
      @test result.final_converged
      @test isapprox(ctrl.x_pu, 0.02; atol = 1e-9)
      @test ctrl.at_limit
      @test !ctrl.converged
      p = get_branch_p_from_to_mw(net, "A", "M2")
      @test p < 70.0
    end

    @testset "bit-identical baseline without the controller" begin
      net_plain = _build_loop_net()
      _, erg_a = runpf!(net_plain, 30, 1e-10, 0)
      @test erg_a == 0
      vm_a, va_a = _solved_state(net_plain)

      # a registered but disabled controller must not touch the solve
      net_off = _build_loop_net()
      ctrl = addSeriesReactanceControl!(net_off; fromBus = "A", toBus = "M2", p_target_mw = 35.0, x_min_pu = 0.02, x_max_pu = 0.30, enabled = false)
      _, erg_b = runpf!(net_off, 30, 1e-10, 0)
      @test erg_b == 0
      vm_b, va_b = _solved_state(net_off)
      @test vm_a == vm_b
      @test va_a == va_b
      @test ctrl.x_pu == 0.20
    end

    @testset "element row vocabulary" begin
      net = _build_loop_net()
      addSeriesReactanceControl!(net; fromBus = "A", toBus = "M2", p_target_mw = 35.0, x_min_pu = 0.02, x_max_pu = 0.30)
      rows = [r for r in controllableElements(net) if r.actuator == :series_x_pu]
      @test length(rows) == 1
      row = rows[1]
      @test row.quantity == :branch_active_power
      @test row.element == "branch@A-M2"
      @test row.device == "TCSC (series compensation)"
      @test row.target == "A->M2"
      @test row.target_value == 35.0
      @test row.actuator_min == 0.02
      @test row.actuator_max == 0.30
      @test row.discrete == false
    end

    @testset "registration validation" begin
      net = _build_trafo_net()
      # transformer branch rejected (taps own transformer reactance)
      @test_throws ErrorException addSeriesReactanceControl!(net; fromBus = "A", toBus = "T", p_target_mw = 10.0, x_min_pu = 0.02, x_max_pu = 0.30)

      net2 = _build_loop_net()
      # missing branch
      @test_throws ErrorException addSeriesReactanceControl!(net2; fromBus = "A", toBus = "B", p_target_mw = 10.0, x_min_pu = 0.02, x_max_pu = 0.30)
      # reversed orientation is caught with a hint, not silently accepted
      @test_throws ErrorException addSeriesReactanceControl!(net2; fromBus = "M2", toBus = "A", p_target_mw = 10.0, x_min_pu = 0.02, x_max_pu = 0.30)
      # inverted range
      @test_throws ErrorException addSeriesReactanceControl!(net2; fromBus = "A", toBus = "M2", p_target_mw = 10.0, x_min_pu = 0.30, x_max_pu = 0.02)
      # starting reactance outside the range
      @test_throws ErrorException addSeriesReactanceControl!(net2; fromBus = "A", toBus = "M2", p_target_mw = 10.0, x_min_pu = 0.25, x_max_pu = 0.40)
      # range end inside the impedance-magnitude guard: |r + jx| must stay
      # above eps_z; r = 0.02 keeps the sign crossing itself admissible,
      # so provoke the guard with a range end exactly at -r on a branch
      # with tiny resistance
      net3 = Net(name = "tcsc_test_eps", baseMVA = 100.0)
      for b in ("A", "B2")
        addBus!(net = net3, busName = b, vn_kV = 110.0)
      end
      addProsumer!(net = net3, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net3, busName = "B2", type = "ENERGYCONSUMER", p = 10.0, q = 2.0)
      addPIModelACLine!(net = net3, fromBus = "A", toBus = "B2", r_pu = 1e-6, x_pu = 0.10, b_pu = 0.0, status = 1)
      # zero crossing with near-zero resistance: |z| dips to |r| = 1e-6 < eps_z
      @test_throws ErrorException addSeriesReactanceControl!(net3; fromBus = "A", toBus = "B2", p_target_mw = 5.0, x_min_pu = -0.05, x_max_pu = 0.15)
      # range end at x = 0 with near-zero resistance is rejected as well
      @test_throws ErrorException addSeriesReactanceControl!(net3; fromBus = "A", toBus = "B2", p_target_mw = 5.0, x_min_pu = 0.0, x_max_pu = 0.15)
      # a second controller on the same branch is rejected
      net4 = _build_loop_net()
      addSeriesReactanceControl!(net4; fromBus = "A", toBus = "M2", p_target_mw = 35.0, x_min_pu = 0.02, x_max_pu = 0.30)
      @test_throws ErrorException addSeriesReactanceControl!(net4; fromBus = "A", toBus = "M2", p_target_mw = 20.0, x_min_pu = 0.02, x_max_pu = 0.30)
    end

    @testset "deadband: met target moves nothing" begin
      net = _build_loop_net()
      # baseline corridor-2 flow is about 27.0 MW; a target inside the
      # deadband around the baseline must leave the actuator untouched
      _, erg0 = runpf!(net, 30, 1e-10, 0)
      @test erg0 == 0
      calcNetLosses!(net)
      p0 = get_branch_p_from_to_mw(net, "A", "M2")
      ctrl = addSeriesReactanceControl!(net; fromBus = "A", toBus = "M2", p_target_mw = p0, x_min_pu = 0.02, x_max_pu = 0.30, deadband_p_mw = 0.5)
      result = run_sparlectra(net = net)
      @test result.final_converged
      @test ctrl.converged
      @test ctrl.x_pu == 0.20
      @test ctrl.prev_x_pu === nothing
    end
  end
end

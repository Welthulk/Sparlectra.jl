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

# file: test/test_pv_voltage_residuals.jl

function _create_pv_voltage_regression_net(; vset::Float64 = 1.05)
  net = Net(name = "pv_voltage_residual_regression", baseMVA = 100.0)

  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)
  addBus!(net = net, busName = "PV", vn_kV = 110.0)

  addACLine!(net = net, fromBus = "Slack", toBus = "Load", length = 10.0, r = 0.02, x = 0.20, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "Load", toBus = "PV", length = 10.0, r = 0.02, x = 0.20, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "Slack", toBus = "PV", length = 12.0, r = 0.02, x = 0.25, c_nf_per_km = 0.0, tanδ = 0.0)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = 40.0, q = 12.0)
  addProsumer!(net = net, busName = "PV", type = "SYNCHRONOUSMACHINE", p = 25.0, vm_pu = vset, va_deg = 0.0, qMin = -500.0, qMax = 500.0)

  # Keep the voltage setpoint explicit: some network-construction paths preserve
  # an existing default bus voltage when a regulating prosumer is added.
  setPVBusVset!(net = net, busName = "PV", vm_pu = vset)
  return net
end

function _rectangular_one_step_voltage(initial_vm::Float64, vset::Float64)
  net = _create_pv_voltage_regression_net(vset = vset)
  Y = Sparlectra.createYBUS(net = net, sparse = true, printYBUS = false)
  S = Sparlectra.buildComplexSVec(net)
  V0, slack_idx = Sparlectra.initialVrect(net; flatstart = true)
  bus_types = [:Slack, :PQ, :PV]
  Vset = [1.0, 1.0, vset]
  V0[3] = initial_vm + 0.0im

  V1 = Sparlectra.complex_newton_step_rectangular(Y, V0, S; slack_idx = slack_idx, bus_types = bus_types, Vset = Vset)
  F0 = Sparlectra.mismatch_rectangular(Y, V0, S, bus_types, Vset, slack_idx)
  return abs(V1[3]), F0[4]
end

function _create_phase_shifted_pv_angle_regression_net()
  net = Net(name = "phase_shifted_pv_angle_regression", baseMVA = 100.0)

  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "PV", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)

  addPIModelTrafo!(net = net, fromBus = "Slack", toBus = "PV", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1, ratio = 1.0, shift_deg = 30.0)
  addACLine!(net = net, fromBus = "PV", toBus = "Load", length = 10.0, r = 0.02, x = 0.20, c_nf_per_km = 0.0, tanδ = 0.0)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 50.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "PV", type = "SYNCHRONOUSMACHINE", p = 30.0, vm_pu = 1.03, qMin = -500.0, qMax = 500.0)
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = 35.0, q = 10.0)

  pv_idx = geNetBusIdx(net = net, busName = "PV")
  setVmVa!(node = net.nodeVec[pv_idx], vm_pu = 0.97, va_deg = 20.0)
  setPVBusVset!(net = net, busName = "PV", vm_pu = 1.03)
  return net
end

function run_pv_voltage_residual_tests()
  @testset "PV voltage residual sign conventions" begin
    initial_vm = 1.0
    vset = 1.05

    @testset "voltage magnitude replacement preserves an existing angle" begin
      V_old = 0.97 * cis(deg2rad(-100.0))
      V_new = Sparlectra._apply_voltage_magnitude_preserving_angle(V_old, 1.03)

      @test isapprox(abs(V_new), 1.03; atol = 1e-12, rtol = 0.0)
      @test isapprox(rad2deg(angle(V_new)), -100.0; atol = 1e-12, rtol = 0.0)
      @test Sparlectra._apply_voltage_magnitude_preserving_angle(0.0 + 0.0im, 1.03) == 1.03 + 0.0im

      # Regression guard: assigning only the real magnitude is the old bug pattern.
      V_wrong = 1.03 + 0.0im
      @test !isapprox(angle(V_wrong), angle(V_old); atol = 1e-12, rtol = 0.0)
    end

    @testset "rectangular PV row moves toward setpoint after one Newton step" begin
      vm_after, residual = _rectangular_one_step_voltage(initial_vm, vset)
      @test residual < 0.0
      @test abs(vm_after - vset) < abs(initial_vm - vset)
      @test isapprox(vm_after, vset; atol = 1e-5)
    end

    @testset "solver paths converge with PV voltage setpoint" begin
      net = _create_pv_voltage_regression_net(vset = vset)
      _, erg = runpf!(net, 40, 1e-9, 0; method = :rectangular)
      @test erg == 0
      @test isapprox(net.nodeVec[3]._vm_pu, vset; atol = 1e-7)
      @test getNodeType(net.nodeVec[3]) == Sparlectra.PV

      for kwargs in ((; method = :polar_full), (; method = :classic), (; method = :polar))
        unsupported_net = _create_pv_voltage_regression_net(vset = vset)
        @test_throws ArgumentError runpf!(unsupported_net, 40, 1e-9, 0; kwargs...)
      end
    end

    @testset "rectangular MATPOWER flat-start seed does not replace imported PV setpoint" begin
      mpc = _synthetic_pv_vg_mismatch_case()
      net = Sparlectra.createNetFromMatPowerCase(mpc = mpc, log = false, flatstart = true, matpower_pv_voltage_source = :gen_vg)
      slack_idx = geNetBusIdx(net = net, busName = "1")
      pv_idx = geNetBusIdx(net = net, busName = "2")

      # Mimic flatstart_voltage_mode=profile_blend with profile_source: matpower_reference: node voltages are start
      # guesses, while regulating prosumer vm_pu stores the imported setpoint.
      net.nodeVec[slack_idx]._vm_pu = 0.5 * (net.nodeVec[slack_idx]._vm_pu + 1.00)
      net.nodeVec[pv_idx]._vm_pu = 0.5 * (net.nodeVec[pv_idx]._vm_pu + 1.02)

      _, erg = runpf!(
        net,
        40,
        1e-9,
        0;
        method = :rectangular,
        opt_flatstart = false,
        start_projection = true,
        start_projection_try_dc_start = true,
        start_projection_try_blend_scan = true,
      )

      rows = Sparlectra.MatpowerIO.pv_voltage_reference_rows(mpc; net = net, matpower_pv_voltage_source = :gen_vg)
      pv_row = rows[findfirst(row -> row.busI == 2, rows)]
      @test erg == 0
      @test getNodeType(net.nodeVec[pv_idx]) == Sparlectra.PV
      @test isapprox(pv_row.imported_vset, 1.04; atol = 1e-12)
      @test isapprox(net.nodeVec[pv_idx]._vm_pu, pv_row.imported_vset; atol = 1e-7)
      @test abs(pv_row.dvm_vset) <= 1e-7
      @test abs(pv_row.dvm_bus) > 0.015
    end

    @testset "phase-shifted PV start preserves angle and avoids wrong branch" begin
      net = _create_phase_shifted_pv_angle_regression_net()
      slack_idx = geNetBusIdx(net = net, busName = "Slack")
      pv_idx = geNetBusIdx(net = net, busName = "PV")

      V_start, _ = Sparlectra.initialVrect(net; flatstart = false)
      V_flat, _ = Sparlectra.initialVrect(net; flatstart = true)
      @test isapprox(abs(V_start[pv_idx]), 1.03; atol = 1e-12, rtol = 0.0)
      @test isapprox(rad2deg(angle(V_start[pv_idx])), 20.0; atol = 1e-12, rtol = 0.0)
      @test abs(rad2deg(angle(V_start[pv_idx]))) > 1.0
      @test isapprox(rad2deg(angle(V_start[slack_idx])), 50.0; atol = 1e-12, rtol = 0.0)

      # Intentional flat start: only non-reference bus angles are initialized to zero.
      @test isapprox(angle(V_flat[pv_idx]), 0.0; atol = 1e-12, rtol = 0.0)
      @test isapprox(rad2deg(angle(V_flat[slack_idx])), 50.0; atol = 1e-12, rtol = 0.0)

      _, erg = runpf!(
        net,
        50,
        1e-9,
        0;
        method = :rectangular,
        opt_flatstart = false,
        start_projection = true,
        wrong_branch_detection = :fail,
        wrong_branch_min_vm_pu = 0.5,
        wrong_branch_max_branch_angle_deg = 90.0,
      )
      status = Sparlectra.rectangular_pf_status(net)
      V_result = Sparlectra.buildVoltageVector(net)

      @test erg == 0
      @test status.numerical_converged === true
      @test status.final_converged === true
      @test status.wrong_branch_status === :ok
      @test minimum(abs.(V_result)) > 0.5
      @test isapprox(abs(V_result[pv_idx]), 1.03; atol = 1e-7, rtol = 0.0)
      @test abs(rad2deg(angle(V_result[pv_idx]))) > 1.0
    end
  end
end

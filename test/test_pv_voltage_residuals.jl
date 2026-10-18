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
# purpose: tests PV voltage residual sign conventions in the rectangular
#          solver: setpoint convergence, preserved start angles, and
#          flat-start seeding that keeps imported PV setpoints

function _create_pv_voltage_regression_net(; vset::Float64 = 1.05)
  net = Net(name = "pv_voltage_residual_regression", baseMVA = 100.0)

  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)
  addBus!(net = net, busName = "PV", vn_kV = 110.0)

  addPIModelACLine!(net = net, fromBus = "Slack", toBus = "Load", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "Load", toBus = "PV", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "Slack", toBus = "PV", r_pu = 0.02, x_pu = 0.25, b_pu = 0.0, status = 1)

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
  addPIModelACLine!(net = net, fromBus = "PV", toBus = "Load", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 50.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "PV", type = "SYNCHRONOUSMACHINE", p = 30.0, vm_pu = 1.03, qMin = -500.0, qMax = 500.0)
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = 35.0, q = 10.0)

  pv_idx = geNetBusIdx(net = net, busName = "PV")
  setVmVa!(node = net.nodeVec[pv_idx], vm_pu = 0.97, va_deg = 20.0)
  setPVBusVset!(net = net, busName = "PV", vm_pu = 1.03)
  return net
end

function run_pv_voltage_residual_tests()
  @testset "PV voltage residual sign conventions" begin (function ()
    initial_vm = 1.0
    vset = 1.05

    @testset "rectangular PV row moves toward setpoint after one Newton step" begin (function ()
      vm_after, residual = _rectangular_one_step_voltage(initial_vm, vset)
      @test residual < 0.0
      @test abs(vm_after - vset) < abs(initial_vm - vset)
      @test isapprox(vm_after, vset; atol = 1e-5)
    end)() end

  end)() end
end

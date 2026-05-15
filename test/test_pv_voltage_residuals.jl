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

function _full_polar_identity_one_step_voltage(initial_vm::Float64, vset::Float64)
  net = _create_pv_voltage_regression_net(vset = vset)
  Y = Sparlectra.createYBUS(net = net, sparse = true, printYBUS = false)
  busVec, _ = Sparlectra.getBusData(net.nodeVec, net.baseMVA, true)
  _, slackIdx = Sparlectra.getBusTypeVec(busVec)
  n_pq = count(bus -> bus.type == Sparlectra.PQ, busVec)
  n_pv = count(bus -> bus.type == Sparlectra.PV, busVec)
  pv_idx = findfirst(bus -> bus.type == Sparlectra.PV, busVec)
  @test pv_idx !== nothing

  busVec[pv_idx].vm_pu = initial_vm
  Vset = [(bus.type == Sparlectra.PV) ? vset : 1.0 for bus in busVec]
  Δ = Sparlectra.residuum_full_withPV(Y, busVec, Vset, n_pq, n_pv, false)
  J = Sparlectra.calcJacobian_withPVIdentity(Y, busVec, Sparlectra.adjacentBranches(Y, false), slackIdx, n_pq, n_pv; sparse = true)
  Δx = J \ Δ

  pv_row = 2 * ((pv_idx >= slackIdx) ? (pv_idx - 1) : pv_idx)
  return initial_vm * (1.0 + Δx[pv_row]), Δ[pv_row], Δx[pv_row]
end

function _rectangular_one_step_voltage(initial_vm::Float64, vset::Float64)
  net = _create_pv_voltage_regression_net(vset = vset)
  Y = Sparlectra.createYBUS(net = net, sparse = true, printYBUS = false)
  S = Sparlectra.buildComplexSVec(net)
  V0, slack_idx = Sparlectra.initialVrect(net; flatstart = true)
  bus_types = [:Slack, :PQ, :PV]
  Vset = [1.0, 1.0, vset]
  V0[3] = initial_vm + 0.0im

  V1 = Sparlectra.complex_newton_step_rectangular(Y, V0, S; slack_idx = slack_idx, bus_types = bus_types, Vset = Vset, use_sparse = true)
  F0 = Sparlectra.mismatch_rectangular(Y, V0, S, bus_types, Vset, slack_idx)
  return abs(V1[3]), F0[4]
end

function run_pv_voltage_residual_tests()
  @testset "PV voltage residual sign conventions" begin
    initial_vm = 1.0
    vset = 1.05

    @testset "full polar PV identity row moves toward setpoint after one Newton step" begin
      vm_after, residual, relative_step = _full_polar_identity_one_step_voltage(initial_vm, vset)
      @test residual > 0.0
      @test relative_step > 0.0
      @test abs(vm_after - vset) < abs(initial_vm - vset)
      @test isapprox(vm_after, vset; atol = 1e-12)
    end

    @testset "rectangular PV row moves toward setpoint after one Newton step" begin
      vm_after, residual = _rectangular_one_step_voltage(initial_vm, vset)
      @test residual < 0.0
      @test abs(vm_after - vset) < abs(initial_vm - vset)
      @test isapprox(vm_after, vset; atol = 1e-5)
    end

    @testset "solver paths converge with PV voltage setpoint" begin
      for (method, kwargs) in [
        (:rectangular, (; opt_sparse = true)),
        (:rectangular, (; opt_fd = true, opt_sparse = true)),
        (:polar_full, (; opt_sparse = true)),
        (:classic, (; opt_sparse = true)),
      ]
        net = _create_pv_voltage_regression_net(vset = vset)
        _, erg = runpf!(net, 40, 1e-9, 0; method = method, kwargs...)
        @test erg == 0
        @test isapprox(net.nodeVec[3]._vm_pu, vset; atol = 1e-7)
        @test getNodeType(net.nodeVec[3]) == Sparlectra.PV
      end
    end

    @testset "rectangular MATPOWER flat-start seed does not replace imported PV setpoint" begin
      mpc = _synthetic_pv_vg_mismatch_case()
      net = Sparlectra.createNetFromMatPowerCase(mpc = mpc, log = false, flatstart = true, matpower_pv_voltage_source = :gen_vg)
      slack_idx = geNetBusIdx(net = net, busName = "1")
      pv_idx = geNetBusIdx(net = net, busName = "2")

      # Mimic flatstart_voltage_mode=bus_vm_va_blend: node voltages are start
      # guesses, while regulating prosumer vm_pu stores the imported setpoint.
      net.nodeVec[slack_idx]._vm_pu = 0.5 * (net.nodeVec[slack_idx]._vm_pu + 1.00)
      net.nodeVec[pv_idx]._vm_pu = 0.5 * (net.nodeVec[pv_idx]._vm_pu + 1.02)

      _, erg = runpf!(
        net,
        40,
        1e-9,
        0;
        method = :rectangular,
        opt_sparse = true,
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
  end
end

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

# file: test/test_dc_powerflow.jl
#
# Tests the standalone DC power flow (rundcpf!, power_flow.solver=:dc):
# reference-value comparisons, phase-shifter (Pfinj) coverage, island
# handling, config/dispatch wiring, and the seed_ac_start convenience path.
#
# Reference angles/flows below were computed with a from-scratch,
# Sparlectra-code-independent plain-LinearAlgebra script (textbook DC-PF /
# MATPOWER makeBdc formula: b=1/(x*tap), Bbus[f,f]+=b etc., theta_slack=0),
# cross-checked by the exact power-balance identity (slack injection ==
# total load - other generation), which matches to machine precision in
# every case below. This is a self-consistency check against an
# independently-implemented reference, not an actual MATPOWER/Octave binary
# comparison (neither is available in this environment).

using Sparlectra
using Test

"""
    _dc3_net(; shift_deg=0.0) -> Net

Small 3-bus programmatic fixture (slack B1, PV generator + load at B2, pure
load at B3) used as the hand-verified DC-PF reference case. `shift_deg != 0`
replaces the B2-B3 line with a phase-shifting transformer of the same
series reactance, to exercise `Pfinj`.
"""
function _dc3_net(; shift_deg::Float64 = 0.0)::Net
  net = Net(name = "dc3", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.02, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B3", r_pu = 0.03, x_pu = 0.12, b_pu = 0.0, status = 1)
  if shift_deg == 0.0
    addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.02, x_pu = 0.10, b_pu = 0.0, status = 1)
  else
    addPIModelTrafo!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.02, x_pu = 0.10, b_pu = 0.0, status = 1, ratio = 1.0, shift_deg = shift_deg)
  end
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "B1")
  addProsumer!(net = net, busName = "B2", type = "GENERATOR", p = 25.0, q = 0.0, vm_pu = 1.01)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 20.0, q = 10.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 45.0, q = 15.0)
  refreshBusTypesFromProsumers!(net)
  return net
end

function run_dc_powerflow_tests()
  @testset "DC power flow (rundcpf!)" begin
    @testset "Reference values: 3-bus hand-verified fixture" begin
      net = _dc3_net()
      report = rundcpf!(net)
      @test report.metadata.converged
      @test report.metadata.solver === :dc
      va = Dict(row.bus_name => row.va_deg for row in report.nodes)
      @test isapprox(va["B1"], 0.0; atol = 1e-9)
      @test isapprox(va["B2"], -0.656992; atol = 1e-4)
      @test isapprox(va["B3"], -1.764710; atol = 1e-4)
      pf = Dict((row.from_bus, row.to_bus) => row.p_from_MW for row in report.branches)
      @test isapprox(pf[(1, 2)], 14.3333; atol = 1e-3)
      @test isapprox(pf[(1, 3)], 25.6667; atol = 1e-3)
      @test isapprox(pf[(2, 3)], 19.3333; atol = 1e-3)
      # Slack MW injection is the *solved* value (315-...-style residual),
      # not the (unused, for DC) specified generator setpoint.
      slack_row = only(filter(row -> row.bus_name == "B1", report.nodes))
      @test isapprox(slack_row.p_gen_MW, 40.0; atol = 1e-6)
      # Lossless model: p_to_MW is exactly -p_from_MW for every branch.
      @test all(isapprox(row.p_to_MW, -row.p_from_MW; atol = 1e-12) for row in report.branches)
    end

    @testset "Reference values: MATPOWER case9" begin
      net = createNetFromMatPowerFile(filename = ensure_casefile("case9.m"))
      report = rundcpf!(net)
      @test report.metadata.converged
      va_expected = [0.0, 9.796019, 5.06056, -2.211159, -3.738091, 2.206657, 0.822441, 3.959011, -4.0634]
      va = [row.va_deg for row in sort(report.nodes, by = r -> r.bus)]
      @test all(isapprox(a, b; atol = 1e-3) for (a, b) in zip(va, va_expected))
      pf_expected = Dict((1, 4) => 67.0, (4, 5) => 28.9674, (5, 6) => -61.0326, (3, 6) => 85.0, (6, 7) => 23.9674, (7, 8) => -76.0326, (8, 2) => -163.0, (8, 9) => 86.9674, (9, 4) => -38.0326)
      for row in report.branches
        key = (row.from_bus, row.to_bus)
        haskey(pf_expected, key) && @test isapprox(row.p_from_MW, pf_expected[key]; atol = 1e-3)
      end
      slack_row = only(filter(row -> row.bus == 1, report.nodes))
      @test isapprox(slack_row.p_gen_MW, 67.0; atol = 1e-6)
    end

    @testset "Reference values: MATPOWER case14 (off-nominal-tap transformers)" begin
      net = createNetFromMatPowerFile(filename = ensure_casefile("case14.m"))
      report = rundcpf!(net)
      @test report.metadata.converged
      va_expected = [0.0, -5.012011, -12.953663, -10.583667, -9.093894, -14.852079, -13.907055, -13.907055, -15.694689, -15.974123, -15.61885, -15.967077, -16.139704, -17.188288]
      va = [row.va_deg for row in sort(report.nodes, by = r -> r.bus)]
      @test all(isapprox(a, b; atol = 1e-3) for (a, b) in zip(va, va_expected))
      slack_row = only(filter(row -> row.bus == 1, report.nodes))
      @test isapprox(slack_row.p_gen_MW, 219.0; atol = 1e-6)
    end

    @testset "Phase-shifter coverage (Pfinj)" begin
      net = _dc3_net(shift_deg = 5.0)
      report = rundcpf!(net)
      @test report.metadata.converged
      va = Dict(row.bus_name => row.va_deg for row in report.nodes)
      @test isapprox(va["B2"], 0.676342; atol = 1e-4)
      @test isapprox(va["B3"], -3.764710; atol = 1e-4)
      pf = Dict((row.from_bus, row.to_bus) => row.p_from_MW for row in report.branches)
      @test isapprox(pf[(1, 2)], -14.755488; atol = 1e-3)
      @test isapprox(pf[(1, 3)], 54.755488; atol = 1e-3)
      @test isapprox(pf[(2, 3)], -9.755488; atol = 1e-3)
      # Total real power is conserved (lossless): shifting redistributes flow,
      # it never changes the slack's required injection.
      slack_row = only(filter(row -> row.bus_name == "B1", report.nodes))
      @test isapprox(slack_row.p_gen_MW, 40.0; atol = 1e-6)
    end

    @testset "Property: net injection balances to zero (lossless)" begin
      for (casefile, atol) in (("case9.m", 1e-6), ("case14.m", 1e-6))
        net = createNetFromMatPowerFile(filename = ensure_casefile(casefile))
        report = rundcpf!(net)
        total = sum(row.p_gen_MW - row.p_load_MW for row in report.nodes)
        @test isapprox(total, 0.0; atol = atol)
      end
    end

    @testset "Island coverage: two independent islands" begin
      island_net = Net(name = "dc_islands", baseMVA = 100.0)
      for busName in ("A1", "A2", "B1", "B2")
        addBus!(net = island_net, busName = busName, vn_kV = 110.0)
      end
      addPIModelACLine!(net = island_net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = island_net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
      addProsumer!(net = island_net, busName = "A1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "A1")
      addProsumer!(net = island_net, busName = "A2", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
      addProsumer!(net = island_net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
      addProsumer!(net = island_net, busName = "B2", type = "ENERGYCONSUMER", p = 8.0, q = 2.0)
      refreshBusTypesFromProsumers!(island_net)

      report = rundcpf!(island_net)
      @test report.metadata.converged
      pf = Dict((row.from_bus, row.to_bus) => row.p_from_MW for row in report.branches)
      # Each island's slack supplies exactly its own island's load (no cross-island coupling).
      @test isapprox(only(v for (k, v) in pf if k[1] in (1, 2) || k[2] in (1, 2)), 10.0; atol = 1e-6)
      @test isapprox(only(v for (k, v) in pf if k[1] in (3, 4) || k[2] in (3, 4)), 8.0; atol = 1e-6)
    end

    @testset "Config validation: power_flow.solver=dc / power_flow.dc.*" begin
      default_cfg = Sparlectra.SparlectraConfig()
      @test default_cfg.powerflow.dc.angle_reference_deg == 0.0
      @test default_cfg.powerflow.dc.ignore_out_of_service === true

      cfg_dc = Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("solver" => "dc", "dc" => Dict("angle_reference_deg" => 3.5))))
      @test cfg_dc.powerflow.solver === :dc
      @test cfg_dc.powerflow.dc.angle_reference_deg == 3.5

      @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("solver" => "not_a_solver")))
    end

    @testset "Controller + DC solver rejection" begin
      mpc = Sparlectra.MatpowerIO.read_case(ensure_casefile("case14.m"))
      net = Sparlectra.createNetFromMatPowerCase(mpc = mpc, flatstart = true)
      addPowerTransformerControl!(net; trafo = "B_2WT_1_4_7", mode = :voltage, target_bus = "7", target_vm_pu = 1.0, control_ratio = true)
      @test length(Sparlectra.collect_outer_controllers(net)) == 1

      cfg = Sparlectra.SparlectraConfig(powerflow = Sparlectra.PowerFlowConfig(solver = :dc), output = OutputConfig(logfile_results = :off))
      @test_throws ArgumentError run_sparlectra(net = net, config = cfg)
    end

    @testset "run_sparlectra dispatch: power_flow.solver=:dc" begin
      cfg = Sparlectra.SparlectraConfig(powerflow = Sparlectra.PowerFlowConfig(solver = :dc), output = OutputConfig(logfile_results = :off))
      result = run_sparlectra(casefile = "case9.m", path = joinpath(dirname(ensure_casefile("case9.m"))), config = cfg)
      @test result.method === :dc
      @test result.final_converged
      @test result.diagnostics.solver === :dc
      slack_row = only(filter(n -> n.busIdx == 1, result.net.nodeVec))
      @test isapprox(slack_row._pƩGen, 67.0; atol = 1e-6)
    end

    @testset "angle_reference_deg is an exact uniform shift" begin
      net0 = _dc3_net()
      report0 = rundcpf!(net0)
      netr = _dc3_net()
      reportr = rundcpf!(netr; angle_reference_deg = 10.0)
      for (row0, rowr) in zip(report0.nodes, reportr.nodes)
        @test isapprox(rowr.va_deg, row0.va_deg + 10.0; atol = 1e-9)
      end
      # Flows are unaffected by a uniform reference shift (only angle *differences* matter).
      for (b0, br) in zip(report0.branches, reportr.branches)
        @test isapprox(br.p_from_MW, b0.p_from_MW; atol = 1e-9)
      end
    end

    @testset "seed_ac_start: DC-seeded AC Newton-Raphson solve" begin
      net = createNetFromMatPowerFile(filename = ensure_casefile("case9.m"))
      report = rundcpf!(net; seed_ac_start = true)
      @test report.metadata.seed_ac_start === true
      @test report.metadata.ac_converged === true
      @test report.metadata.ac_iterations isa Int
      # net now holds the AC-converged solution, not the DC one: voltage
      # magnitudes are no longer flattened to 1.0 everywhere.
      @test any(!isapprox(n._vm_pu, 1.0; atol = 1e-9) for n in net.nodeVec)
      @test dc_pf_status(net) !== nothing
      @test dc_pf_status(net).numerical_converged
    end

    @testset "DC power-flow status registry is separate from the AC one" begin
      net = _dc3_net()
      rundcpf!(net)
      @test dc_pf_status(net) !== nothing
      @test dc_pf_status(net).solver === :dc
      @test Sparlectra.rectangular_pf_status(net) === nothing
    end
  end
end

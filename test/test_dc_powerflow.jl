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
# purpose: tests the standalone DC power flow (rundcpf!): reference values
#          for hand-verified and MATPOWER cases, phase shifters, islands,
#          config and dispatch wiring, and DC-seeded AC starts
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
  @testset "Phase-shifter coverage (Pfinj)" begin (function ()
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
  end)() end

  @testset "Island coverage: two independent islands" begin (function ()
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
  end)() end

  @testset "Controller + DC solver rejection" begin (function ()
    # load_fixture_net: the shipped sp_case14 carries a REAL declared tap
    # controller, so the rejection is tested against the file's own
    # controller instead of one hand-attached to a downloaded case
    net = Sparlectra.importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")))
    @test length(Sparlectra.collect_outer_controllers(net)) == 1

    cfg = Sparlectra.SparlectraConfig(powerflow = Sparlectra.PowerFlowConfig(solver = :dc), output = OutputConfig(logfile_results = :off))
    @test_throws ArgumentError run_sparlectra(net = net, config = cfg)
  end)() end

  @testset "run_sparlectra dispatch: power_flow.solver=:dc" begin (function ()
    cfg = Sparlectra.SparlectraConfig(powerflow = Sparlectra.PowerFlowConfig(solver = :dc), output = OutputConfig(logfile_results = :off))
    result = run_sparlectra(casefile = "sp_case5.scf.json", path = joinpath(dirname(@__DIR__), "data", "scf"), config = cfg)
    @test result.method === :dc
    @test result.final_converged
    @test result.diagnostics.solver === :dc
    # hand-checkable: 32 MW load minus the 12 MW PV machine, lossless DC
    slack_row = result.net.nodeVec[result.net.busDict["Sandau_110"]]
    @test isapprox(slack_row._pƩGen, 20.0; atol = 1e-6)
  end)() end

  @testset "angle_reference_deg is an exact uniform shift" begin (function ()
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
    # a DC run fills the DC status of the net and leaves the AC one empty
    @test dc_pf_status(net0) !== nothing
    @test dc_pf_status(net0).solver === :dc
    @test Sparlectra.rectangular_pf_status(net0) === nothing
  end)() end

  @testset "seed_ac_start: DC-seeded AC Newton-Raphson solve" begin (function ()
    # load_fixture_net: property-style, runs on the shipped case
    net = load_fixture_net("sp_case5")
    report = rundcpf!(net; seed_ac_start = true)
    @test report.metadata.seed_ac_start === true
    @test report.metadata.ac_converged === true
    @test report.metadata.ac_iterations isa Int
    # net now holds the AC-converged solution, not the DC one: voltage
    # magnitudes are no longer flattened to 1.0 everywhere.
    @test any(!isapprox(n._vm_pu, 1.0; atol = 1e-9) for n in net.nodeVec)
    @test dc_pf_status(net) !== nothing
    @test dc_pf_status(net).numerical_converged
  end)() end

  @testset "Config validation: solver=dc, dc.angle_reference_deg, dc_seed_unconditional" begin (function ()
    cfg_dc = Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("solver" => "dc", "dc" => Dict("angle_reference_deg" => 3.5))))
    @test cfg_dc.powerflow.solver === :dc
    @test cfg_dc.powerflow.dc.angle_reference_deg == 3.5

    cfg_seed = Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("start_mode" => Dict("dc_seed_unconditional" => true))))
    @test cfg_seed.powerflow.start_mode.dc_seed_unconditional === true
    @test cfg_seed.powerflow.solver === :rectangular

    # Redundant with solver=dc (the solver already is the standalone DC power flow).
    @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("solver" => "dc", "start_mode" => Dict("dc_seed_unconditional" => true))))
    # Unsupported for solver=apslf.
    @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("solver" => "apslf", "start_mode" => Dict("dc_seed_unconditional" => true))))
    # Mutually exclusive with the other rectangular-NR start-value source.
    @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("apslf_start" => Dict("enabled" => true), "start_mode" => Dict("dc_seed_unconditional" => true))))
  end)() end

  @testset "closed busbar couplers are contracted" begin (function ()
    # the two sections of a CLOSED coupler are one DC node: on sp_case60
    # they get identical angles by construction (solved apart, they were
    # 3.1 deg off)
    rep60 = rundcpf!(load_fixture_net("sp_case60"))
    va = Dict(row.bus_name => row.va_deg for row in rep60.nodes)
    @test va["Neubach_110"] == va["Neubach_110b"]
    # warmup_casePST: the second section hangs on the link alone, so a
    # missing contraction splits off an island without reference; it
    # solves, and the lossless injections balance to zero
    rep_pst = rundcpf!(load_fixture_net("warmup_casePST"))
    @test rep_pst.metadata.converged
    @test isapprox(sum(row.p_gen_MW - row.p_load_MW for row in rep_pst.nodes), 0.0; atol = 1e-6)
  end)() end

end

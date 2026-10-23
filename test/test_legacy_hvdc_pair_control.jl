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

# file: test/test_legacy_hvdc_pair_control.jl
# purpose: legacy profile: the testsets of test/test_hvdc_pair_control.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_hvdc_pair_control.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_hvdc_pair_control_tests()
  # cut from testset: HVDC pair controller (#297 Draft B)
    @testset "disabled controller keeps the baseline bit-identical" begin (function ()
      net_a, _, _ = _build_hvdc_two_area_net("hvdc_base_a")
      net_b, _, _ = _build_hvdc_two_area_net("hvdc_base_b")
      addHvdcPairControl!(net_b; from_bus = "B2", to_bus = "B4", p_transfer_mw = 999.0, enabled = false)
      # the fixture is two islands by design; the raw legacy runpf! call
      # needs island solving enabled explicitly (PowerFlowConfig defaults
      # to enabled, the positional API does not)
      ite_a, erg_a = runpf!(net_a, 40, 1e-9, 0; method = :rectangular, islands_enabled = true)
      ite_b, erg_b = runpf!(net_b, 40, 1e-9, 0; method = :rectangular, islands_enabled = true)
      @test erg_a == 0 && erg_b == 0
      @test ite_a == ite_b
      vm_a = [something(n._vm_pu, NaN) for n in net_a.nodeVec]
      vm_b = [something(n._vm_pu, NaN) for n in net_b.nodeVec]
      @test vm_a == vm_b
    end)() end

    @testset "YAML declaration round trip (#305 schema)" begin (function ()
      # declared via control.controllers, run through the real config path,
      # compared against the programmatic twin
      net_prog, f_idx, t_idx = _build_hvdc_two_area_net("hvdc_yaml_prog")
      addHvdcPairControl!(net_prog; from_bus = "B2", to_bus = "B4", p_transfer_mw = 90.0, loss_mw = 5.0)
      _hvdc_run!(net_prog)
      net_yaml, yf_idx, yt_idx = _build_hvdc_two_area_net("hvdc_yaml_decl")
      cfg = ControlConfig(controllers = Any[Dict{String,Any}("type" => "hvdc_pair", "name" => "b2b_main", "from_bus" => "B2", "to_bus" => "B4", "p_transfer_mw" => 90.0, "loss_mw" => 5.0)])
      @test applyConfiguredControllers!(net_yaml, cfg) == 1
      ctrl = only(Sparlectra._hvdc_pair_controllers(net_yaml))
      @test ctrl.name == "b2b_main"
      _hvdc_run!(net_yaml)
      @test net_yaml.prosumpsVec[yf_idx].pVal == net_prog.prosumpsVec[f_idx].pVal == -90.0
      @test net_yaml.prosumpsVec[yt_idx].pVal == net_prog.prosumpsVec[t_idx].pVal == 85.0
      # idempotent re-apply
      @test applyConfiguredControllers!(net_yaml, cfg) == 0
    end)() end

    @testset "MATPOWER dcline paired_control mode" begin (function ()
      # two islands joined by one dcline row (17 columns, so the effective
      # PT is recomputed from LOSS0/LOSS1: 80 - (3 + 0.0125 * 80) = 76)
      case = mktempdir()
      path = joinpath(case, "case_b2b.m")
      write(
        path,
        """
function mpc = case_b2b
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
1 3 0 0 0 0 1 1.0 0 380 1 1.1 0.9;
2 1 40 10 0 0 1 1.0 0 380 1 1.1 0.9;
3 2 0 0 0 0 1 1.0 0 380 1 1.1 0.9;
4 1 50 12 0 0 1 1.0 0 380 1 1.1 0.9;
];
mpc.gen = [
1 60 0 300 -300 1.0 100 1 300 0;
3 40 0 300 -300 1.0 100 1 300 0;
];
mpc.branch = [
1 2 0.01 0.08 0 250 250 250 0 0 1 -360 360;
3 4 0.01 0.08 0 250 250 250 0 0 1 -360 360;
];
mpc.dcline = [
2 4 1 80 76 0 0 1.0 1.0 0 100 -50 50 -50 50 3 0.0125;
];
""",
      )
      net_paired = createNetFromMatPowerFile(filename = path, matpower_dcline_mode = :paired_control)
      ctrl = only(Sparlectra._hvdc_pair_controllers(net_paired))
      @test ctrl.name == "DCLINE_1"
      # r0.9.9: the persistent link record exists and carries the controller
      @test length(net_paired.hvdcLinks) == 1
      @test net_paired.hvdcLinks[1].controller_name == "DCLINE_1"
      @test net_paired.hvdcLinks[1].source == :matpower && net_paired.hvdcLinks[1].kind == :p2p
      @test ctrl.p_transfer_mw == 80.0
      @test ctrl.loss_mw == 3.0
      @test ctrl.loss_fraction == 0.0125
      result = _hvdc_run!(net_paired)
      @test result.converged == true
      # invariant matches the MATPOWER loss model exactly
      m = only(net_paired.matpowerDclineMetadata)
      @test net_paired.prosumpsVec[m.from_prosumer].pVal == -80.0
      @test net_paired.prosumpsVec[m.to_prosumer].pVal == 80.0 - (3.0 + 0.0125 * 80.0)
      # consistency anchor: with the target equal to the row's PF, the
      # controlled solve reproduces the pf_injections solve
      net_fixed = createNetFromMatPowerFile(filename = path, matpower_dcline_mode = :pf_injections)
      # r0.9.9: Stage-0 import records the link without a controller
      @test length(net_fixed.hvdcLinks) == 1
      @test net_fixed.hvdcLinks[1].controller_name === nothing
      _, erg = runpf!(net_fixed, 40, 1e-9, 0; method = :rectangular, islands_enabled = true)
      @test erg == 0
      vm_fixed = [something(n._vm_pu, NaN) for n in net_fixed.nodeVec]
      vm_paired = [something(n._vm_pu, NaN) for n in net_paired.nodeVec]
      @test isapprox(vm_fixed, vm_paired; atol = 1e-9)
      # a changed transfer target actually moves the flow
      ctrl.p_transfer_mw = 40.0
      ctrl.p_applied = false
      result2 = _hvdc_run!(net_paired)
      @test result2.converged == true
      @test net_paired.prosumpsVec[m.from_prosumer].pVal == -40.0
      @test net_paired.prosumpsVec[m.to_prosumer].pVal == 40.0 - (3.0 + 0.0125 * 40.0)
    end)() end

  return true
end

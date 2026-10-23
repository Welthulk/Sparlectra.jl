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

# file: test/test_legacy_external_grid.jl
# purpose: legacy profile: the testsets of test/test_external_grid.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here from test_external_grid.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_external_grid_tests()
  # cut from testset: External grid (issue #299)
    @testset "PF invariance: carried SC data changes no power-flow result" begin (function ()
      ref = _eg_fixture()
      ite_ref, erg_ref = runpf!(ref, 30, 1e-8, 0)
      @test erg_ref == 0

      # Same fixture, but the slack is created through addExternalGrid! with
      # finite short-circuit data. The PF side must be exactly the manual
      # slack-prosumer path.
      net = Net(name = "extgrid_fixture", baseMVA = 100.0)
      for b in ("B1", "B2", "B3", "B4")
        addBus!(net = net, busName = b, vn_kV = 110.0)
      end
      addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.06, b_pu = 0.02, status = 1)
      addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.02, x_pu = 0.09, b_pu = 0.02, status = 1)
      addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.015, x_pu = 0.08, b_pu = 0.02, status = 1)
      addPIModelACLine!(net = net, fromBus = "B4", toBus = "B1", r_pu = 0.012, x_pu = 0.07, b_pu = 0.02, status = 1)
      addExternalGrid!(net = net, busName = "B1", vm_pu = 1.02, sk_max_MVA = 3000.0, rx_max = 0.1)
      addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 30.0, vm_pu = 1.01, isRegulated = true)
      addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 40.0, q = 12.0)
      addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 35.0, q = 8.0)

      @test net.slackVec == ref.slackVec
      ite, erg = runpf!(net, 30, 1e-8, 0)
      @test erg == 0
      @test ite == ite_ref
      for (a, b) in zip(_eg_state(net), _eg_state(ref))
        @test a[1] ≈ b[1] atol = 1e-12
        @test a[2] ≈ b[2] atol = 1e-12
      end

      # Conversion into the feeder record: full ENI contract, at add time.
      @test length(net.sc_sources.external_network_injections) == 1
      eni = net.sc_sources.external_network_injections[1]
      @test eni.bus == "B1"
      @test eni.mrid == "native-eni-B1"
      @test eni.maxInitialSymShCCurrent_A ≈ 3000.0e6 / (sqrt(3.0) * 110.0e3) rtol = 1e-12
      @test eni.minInitialSymShCCurrent_A === nothing
      @test eni.maxR1ToX1Ratio == 0.1
      @test eni.minR1ToX1Ratio === nothing
      @test eni.maxR0ToX0Ratio === nothing && eni.governorSCD === nothing
    end)() end

    @testset "internal impedance: non-ideal source (stage 2)" begin (function ()
      build = (; internal, sk) -> begin
        net = Net(name = "extgrid_int", baseMVA = 100.0)
        addBus!(net = net, busName = "B1", vn_kV = 110.0)
        addBus!(net = net, busName = "B2", vn_kV = 110.0)
        addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
        addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 50.0, q = 10.0)
        addExternalGrid!(net = net, busName = "B1", vm_pu = 1.02, sk_max_MVA = sk, rx_max = 0.1, internal_impedance = internal)
        net
      end

      ideal = build(internal = false, sk = 3000.0)
      _, erg_i = runpf!(ideal, 30, 1e-8, 0)
      @test erg_i == 0
      vm_ideal = getNodeVm(ideal.nodeVec[1])

      # Stiff limit: a huge Sk must reproduce the ideal slack. The series
      # impedance scales as baseMVA/Sk, so sk = 1e10 puts the residual
      # voltage drop safely below the 1e-8 pu assertion (measured ~1.6e-8
      # at the original 1e9; impedance, not tolerance, is the knob).
      stiff = build(internal = true, sk = 1.0e10)
      @test length(stiff.nodeVec) == 3
      @test stiff.slackVec == [3]                       # slack moved to the internal bus
      @test Sparlectra._is_aux_node(stiff.nodeVec[3])           # the internal bus is auxiliary
      _, erg_s = runpf!(stiff, 30, 1e-8, 0)
      @test erg_s == 0
      @test abs(getNodeVm(stiff.nodeVec[1]) - vm_ideal) < 1e-8

      # Realistic Sk: terminal voltage droops and the angle shifts; the
      # terminal bus is now an ordinary solved bus.
      weak = build(internal = true, sk = 500.0)
      _, erg_w = runpf!(weak, 30, 1e-8, 0)
      @test erg_w == 0
      @test getNodeVm(weak.nodeVec[1]) < vm_ideal
      @test abs(weak.nodeVec[1]._va_deg) > 0.1

      # The feeder record stays anchored at the physical connection bus and
      # the auxiliary branch is inert in the SC network: the fault at B1
      # reproduces the declared current exactly.
      res = runShortCircuit!(weak; case = :max)
      row = first(r for r in res.rows if r.bus == "B1")
      @test row.ik_kA ≈ 500.0 / (sqrt(3.0) * 110.0) rtol = 1e-6

      # One internal-impedance external grid per bus.
      @test_throws ArgumentError addExternalGrid!(net = weak, busName = "B1", sk_max_MVA = 100.0, internal_impedance = true)
    end)() end

  return true
end

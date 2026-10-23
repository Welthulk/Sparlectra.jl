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

# file: test/test_legacy_cgmes_export.jl
# purpose: legacy profile: the testsets of test/test_cgmes_export.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_cgmes_export.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_cgmes_export_tests()
  # cut from testset: CGMES export identity
    @testset "regulated tap group exports one shared TapChangerControl" begin (function ()
      # #322 export half: master and follower reference the SAME control
      # and both carry controlEnabled, so a reimport regroups them instead
      # of seeing independent (fighting) controllers
      n = Net(name = "tccpar", baseMVA = 100.0)
      for (b, vn) in (("H1", 110.0), ("H2", 110.0), ("L1", 20.0), ("L2", 20.0))
        addBus!(net = n, busName = b, vn_kV = vn)
      end
      addProsumer!(net = n, busName = "H1", type = "EXTERNALNETWORKINJECTION", referencePri = "H1", vm_pu = 1.02, va_deg = 0.0)
      addProsumer!(net = n, busName = "L1", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
      addProsumer!(net = n, busName = "L2", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
      addPIModelACLine!(net = n, fromBus = "H1", toBus = "H2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = n, fromBus = "L1", toBus = "L2", r_pu = 0.02, x_pu = 0.1, b_pu = 0.0, status = 1)
      add2WTrafo!(net = n, fromBus = "H1", toBus = "L1", sn_mva = 63.0, vk_percent = 12.0, vkr_percent = 0.4, pfe_kw = 0.0, i0_percent = 0.0)
      add2WTrafo!(net = n, fromBus = "H1", toBus = "L1", sn_mva = 63.0, vk_percent = 12.0, vkr_percent = 0.4, pfe_kw = 0.0, i0_percent = 0.0)
      for tf in n.trafos
        tf.side1.taps = PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = -9, highStep = 9, neutralStep = 0, voltageIncrement_kV = 110.0 * 0.00625)
      end
      for br in n.branchVec
        br.ratio == 0.0 && continue
        br.has_ratio_tap = true
        br.tap_step = 0.00625
      end
      m = string(n.branchVec[3].branchIdx)
      fl = string(n.branchVec[4].branchIdx)
      addPowerTransformerControl!(n; trafo = m, followers = [fl], mode = :voltage, target_bus = "L1", target_vm_pu = 1.01, deadband_vm_pu = 0.004)
      dir = mktempdir()
      files = writeCGMESFiles(n; path = dir, created = _EXPORT_STAMP)
      eq = read(files[1], String)
      ssh = read(files[3], String)
      @test count("cim:TapChangerControl rdf:ID", eq) == 1
      @test count("TapChanger.TapChangerControl rdf:resource", eq) == 2
      @test count("TapChanger.controlEnabled>true", ssh) == 2
      @test occursin("RegulatingControl.targetDeadband", ssh)
      res = importCGMES(path = dir, name = "tcc_back", tap_control = true)
      ctrls = Sparlectra._tap_controllers(res.net)
      @test length(ctrls) == 1
      c = only(ctrls)
      @test length(c.followers) == 1
      @test isapprox(something(c.target_vm_pu, NaN), 1.01; atol = 1e-9)
      @test isapprox(c.deadband_vm_pu, 0.004; atol = 1e-9)
    end)() end

  return true
end

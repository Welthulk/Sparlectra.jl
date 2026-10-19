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

# file: test/test_legacy_tap_changer_model.jl
# purpose: legacy profile: the testsets of test/test_tap_changer_model.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_tap_changer_model.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_tap_changer_model_tests()
  # cut from the runner function
  @testset "calcTransformerRatio regression on tapped 2WT" begin (function ()
    taps = PowerTransformerTaps(Vn_kV = 110.0, step = 4, lowStep = -4, highStep = 6, neutralStep = 1, voltageIncrement_kV = 1.1)
    w1 = PowerTransformerWinding(110.0, 0.0, 0.12, 0.0, 0.0, 1.0, 0.0, 110.0, 100.0, taps, true, nothing)
    w2 = PowerTransformerWinding(20.0, 0.0, 0.0, 0.0, 0.0, nothing, 0.0, 20.0, 100.0, nothing, true, nothing)
    trafo_comp = Sparlectra.getBranchComp(110.0, 1, 2, 1, "2WT")
    trafo = PowerTransformer(trafo_comp, true, w1, w2, nothing, Sparlectra.Ratio)
    @test trafo.tapSideNumber == 1

    ratio = calcTransformerRatio(trafo)

    base_ratio = (trafo.side1.Vn * trafo.side2.ratedU) / (trafo.side2.Vn * trafo.side1.ratedU)
    corr = 1.0 + (taps.step - taps.neutralStep) * taps.tapStepPercent / 100.0
    expected = base_ratio / corr

    @test isapprox(ratio, expected; atol = 1e-12)
  end)() end

  @testset "Base.show contains convention and fixed highStep typo" begin (function ()
    taps = PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = -4, highStep = 6, neutralStep = 1, voltageIncrement_kV = 1.1)
    txt = sprint(show, taps)
    @test occursin("convention=", txt)
    @test occursin("highStep=", txt)
    @test !occursin("hightStep=", txt)
  end)() end

  return true
end

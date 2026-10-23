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

# file: test/test_tap_changer_model.jl
# purpose: hand calculations for the ratio tap-changer model: neutral voltage,
#          calcRatioTapCorrection, calcRatioTapRange, the Branch tap limits
#          and calcTransformerRatio on a tapped 2WT

function run_tap_changer_model_tests()
  @testset "PowerTransformerTaps: convention check and neutral voltage" begin (function ()
    @test_throws ArgumentError PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = -4, highStep = 6, neutralStep = 1, voltageIncrement_kV = 1.1, convention = :bogus)

    # neutralU from the ratio: 20 kV * 1.05 = 21 kV
    taps_ratio = PowerTransformerTaps(Vn_kV = 20.0, step = 0, lowStep = -2, highStep = 2, neutralStep = 0, voltageIncrement_kV = 0.4, neutralU_ratio = 1.05)
    @test isapprox(taps_ratio.neutralU, 21.0; atol = 1e-12)
    @test isapprox(taps_ratio.neutralU_ratio, 1.05; atol = 1e-12)
    # ratio from neutralU: 19.5 kV / 20 kV = 0.975
    taps_neutral = PowerTransformerTaps(Vn_kV = 20.0, step = 0, lowStep = -2, highStep = 2, neutralStep = 0, voltageIncrement_kV = 0.4, neutralU = 19.5)
    @test isapprox(taps_neutral.neutralU_ratio, 0.975; atol = 1e-12)
  end)() end

  @testset "ratio tap hand calculations" begin (function ()
    # 1.1 kV per step on 110 kV is 1 % per step, neutral step 1
    taps = PowerTransformerTaps(Vn_kV = 110.0, step = 2, lowStep = -4, highStep = 6, neutralStep = 1, voltageIncrement_kV = 1.1)
    @test isapprox(taps.tapStepPercent, 1.0; atol = 1e-12)

    # correction 1 + (n - n0) * 1 %
    @test isapprox(calcRatioTapCorrection(taps; step = taps.neutralStep), 1.0; atol = 1e-12)
    @test isapprox(calcRatioTapCorrection(taps; step = 4), 1.03; atol = 1e-12)
    @test isapprox(calcRatioTapCorrection(taps; step = -2), 0.97; atol = 1e-12)
    @test isapprox(calcRatioTapCorrection(taps), 1.01; atol = 1e-12)   # default step = taps.step = 2

    # band from step -4 (1 - 5 %) to step 6 (1 + 5 %), 1 % per step
    for (lo, hi) in ((-4, 6), (6, -4))   # the second is an inverted range (tapSign -1)
      t = PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = lo, highStep = hi, neutralStep = 1, voltageIncrement_kV = 1.1)
      @test t.tapSign == (hi > lo ? 1 : -1)
      r = calcRatioTapRange(t)
      @test isapprox(r.tap_min, 0.95; atol = 1e-12)
      @test isapprox(r.tap_max, 1.05; atol = 1e-12)
      @test isapprox(r.tap_step, 0.01; atol = 1e-12)
    end

    # 110/20 kV 2WT with the taps on winding 1
    tapped_2wt = function (step::Int)
      t = PowerTransformerTaps(Vn_kV = 110.0, step = step, lowStep = -4, highStep = 6, neutralStep = 1, voltageIncrement_kV = 1.1)
      w1 = PowerTransformerWinding(110.0, 0.0, 0.12, 0.0, 0.0, 1.0, 0.0, 110.0, 100.0, t, true, nothing)
      w2 = PowerTransformerWinding(20.0, 0.0, 0.0, 0.0, 0.0, nothing, 0.0, 20.0, 100.0, nothing, true, nothing)
      return PowerTransformer(Sparlectra.getBranchComp(110.0, 1, 2, 1, "2WT"), true, w1, w2, nothing, Sparlectra.Ratio)
    end

    # a Branch built on the tapped 2WT takes the same band
    trafo = tapped_2wt(0)
    br = Branch(branchIdx = 1, from = 1, to = 2, baseMVA = 100.0, branch = trafo, id = 1, ratio = 1.0, side = 1, vn_kV = 110.0, values_are_pu = true)
    @test isapprox(br.tap_min, 0.95; atol = 1e-12)
    @test isapprox(br.tap_max, 1.05; atol = 1e-12)
    @test isapprox(br.tap_step, 0.01; atol = 1e-12)

    # calcTransformerRatio at step 4: nominal ratio (110/20)/(110/20) = 1,
    # divided by the correction 1.03
    trafo4 = tapped_2wt(4)
    @test trafo4.tapSideNumber == 1
    @test isapprox(calcTransformerRatio(trafo4), 1.0 / 1.03; atol = 1e-12)
  end)() end
end

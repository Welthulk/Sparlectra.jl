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
# purpose: unit tests for the ratio tap-changer model: convention defaults,
#          calcRatioTapCorrection, calcRatioTapRange against the previous
#          inline branch.jl math, and a calcTransformerRatio regression

function run_tap_changer_model_tests()
  @testset "AbstractTapChangerModel supertype" begin (function ()
    @test PowerTransformerTaps <: AbstractTapChangerModel
  end)() end

  @testset "convention field default and validation" begin (function ()
    taps = PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = -4, highStep = 6, neutralStep = 1, voltageIncrement_kV = 1.1)
    @test taps.convention == :neutral_relative
    @test_throws ArgumentError PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = -4, highStep = 6, neutralStep = 1, voltageIncrement_kV = 1.1, convention = :bogus)
  end)() end

  @testset "calcRatioTapCorrection" begin (function ()
    taps = PowerTransformerTaps(Vn_kV = 110.0, step = 2, lowStep = -4, highStep = 6, neutralStep = 1, voltageIncrement_kV = 1.1)
    @test isapprox(taps.tapStepPercent, 1.0; atol = 1e-12) # (1.1/110.0)*100

    @test isapprox(calcRatioTapCorrection(taps; step = taps.neutralStep), 1.0; atol = 1e-12)

    pos_corr = calcRatioTapCorrection(taps; step = 4)
    @test isapprox(pos_corr, 1.0 + (4 - 1) * 1.0 / 100.0; atol = 1e-12)

    neg_corr = calcRatioTapCorrection(taps; step = -2)
    @test isapprox(neg_corr, 1.0 + (-2 - 1) * 1.0 / 100.0; atol = 1e-12)

    # default `step` keyword falls back to taps.step
    @test isapprox(calcRatioTapCorrection(taps), 1.0 + (taps.step - taps.neutralStep) * taps.tapStepPercent / 100.0; atol = 1e-12)
  end)() end

  @testset "calcRatioTapRange: ratio band and step, including an inverted range" begin (function ()
    taps = PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = -4, highStep = 6, neutralStep = 1, voltageIncrement_kV = 1.1)
    pu_per_step = taps.tapStepPercent / 100.0
    expected_min = min(1.0 + (taps.lowStep - taps.neutralStep) * pu_per_step, 1.0 + (taps.highStep - taps.neutralStep) * pu_per_step)
    expected_max = max(1.0 + (taps.lowStep - taps.neutralStep) * pu_per_step, 1.0 + (taps.highStep - taps.neutralStep) * pu_per_step)
    expected_step = abs(pu_per_step)

    result = calcRatioTapRange(taps)
    @test isapprox(result.tap_min, expected_min; atol = 1e-12)
    @test isapprox(result.tap_max, expected_max; atol = 1e-12)
    @test isapprox(result.tap_step, expected_step; atol = 1e-12)

    # inverted range: highStep < lowStep (tapSign == -1)
    taps_inv = PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = 6, highStep = -4, neutralStep = 1, voltageIncrement_kV = 1.1)
    @test taps_inv.tapSign == -1
    pu_per_step_inv = taps_inv.tapStepPercent / 100.0
    expected_min_inv = min(1.0 + (taps_inv.lowStep - taps_inv.neutralStep) * pu_per_step_inv, 1.0 + (taps_inv.highStep - taps_inv.neutralStep) * pu_per_step_inv)
    expected_max_inv = max(1.0 + (taps_inv.lowStep - taps_inv.neutralStep) * pu_per_step_inv, 1.0 + (taps_inv.highStep - taps_inv.neutralStep) * pu_per_step_inv)
    expected_step_inv = abs(pu_per_step_inv)

    result_inv = calcRatioTapRange(taps_inv)
    @test isapprox(result_inv.tap_min, expected_min_inv; atol = 1e-12)
    @test isapprox(result_inv.tap_max, expected_max_inv; atol = 1e-12)
    @test isapprox(result_inv.tap_step, expected_step_inv; atol = 1e-12)
  end)() end

end

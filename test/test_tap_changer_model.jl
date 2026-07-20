# file: test/test_tap_changer_model.jl

function run_tap_changer_model_tests()
  @testset "AbstractTapChangerModel supertype" begin
    @test PowerTransformerTaps <: AbstractTapChangerModel
  end

  @testset "convention field default and validation" begin
    taps = PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = -4, highStep = 6, neutralStep = 1, voltageIncrement_kV = 1.1)
    @test taps.convention == :neutral_relative
    @test_throws ArgumentError PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = -4, highStep = 6, neutralStep = 1, voltageIncrement_kV = 1.1, convention = :bogus)
  end

  @testset "calcRatioTapCorrection" begin
    taps = PowerTransformerTaps(Vn_kV = 110.0, step = 2, lowStep = -4, highStep = 6, neutralStep = 1, voltageIncrement_kV = 1.1)
    @test isapprox(taps.tapStepPercent, 1.0; atol = 1e-12) # (1.1/110.0)*100

    @test isapprox(calcRatioTapCorrection(taps; step = taps.neutralStep), 1.0; atol = 1e-12)

    pos_corr = calcRatioTapCorrection(taps; step = 4)
    @test isapprox(pos_corr, 1.0 + (4 - 1) * 1.0 / 100.0; atol = 1e-12)

    neg_corr = calcRatioTapCorrection(taps; step = -2)
    @test isapprox(neg_corr, 1.0 + (-2 - 1) * 1.0 / 100.0; atol = 1e-12)

    # default `step` keyword falls back to taps.step
    @test isapprox(calcRatioTapCorrection(taps), 1.0 + (taps.step - taps.neutralStep) * taps.tapStepPercent / 100.0; atol = 1e-12)
  end

  @testset "calcRatioTapRange matches previous inline branch.jl math" begin
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
  end

  @testset "calcTransformerRatio regression on tapped 2WT" begin
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
  end

  @testset "Base.show contains convention and fixed highStep typo" begin
    taps = PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = -4, highStep = 6, neutralStep = 1, voltageIncrement_kV = 1.1)
    txt = sprint(show, taps)
    @test occursin("convention=", txt)
    @test occursin("highStep=", txt)
    @test !occursin("hightStep=", txt)
  end
end

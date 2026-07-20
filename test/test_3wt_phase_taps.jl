# file: test/test_3wt_phase_taps.jl

function _3wt_phase_taps_fixtures()
  tmp1 = TransformerModelParameters(sn_MVA = 100.0, vk_percent = 12.0, vkr_percent = 0.5, pk_kW = 0.0, i0_percent = 0.1, p0_kW = 0.0)
  tmp2 = TransformerModelParameters(sn_MVA = 80.0, vk_percent = 12.0, vkr_percent = 0.5, pk_kW = 0.0, i0_percent = 0.1, p0_kW = 0.0)
  tmp3 = TransformerModelParameters(sn_MVA = 20.0, vk_percent = 12.0, vkr_percent = 0.5, pk_kW = 0.0, i0_percent = 0.1, p0_kW = 0.0)
  tapSettings = PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = -8, highStep = 8, neutralStep = 0, voltageIncrement_kV = 1.5)
  return tmp1, tmp2, tmp3, tapSettings
end

function run_3wt_phase_taps_tests()
  @testset "create3WTWindings! existing docstring example: no phase taps" begin
    tmp1, tmp2, tmp3, tapSettings = _3wt_phase_taps_fixtures()
    w1, w2, w3 = create3WTWindings!(u_kV = [110.0, 20.0, 10.0], sn_MVA = [100.0, 80.0, 20.0], addEx_Side = [tmp1, tmp2, tmp3], sh_deg = [0.0, 0.0, 0.0], tap_side = 1, tap = tapSettings)

    @test isnothing(w1.phase_taps)
    @test isnothing(w2.phase_taps)
    @test isnothing(w3.phase_taps)

    # Hand-derived expectations for the docstring example (tap_side = 1, no new keywords).
    @test w1.Vn == 110.0 && w2.Vn == 20.0 && w3.Vn == 10.0
    @test w1.ratedU == 110.0 && w2.ratedU == 20.0 && w3.ratedU == 10.0
    @test w1.ratedS == 100.0 && w2.ratedS == 80.0 && w3.ratedS == 20.0
    @test w1.shift_degree == 0.0 && w2.shift_degree == 0.0 && w3.shift_degree == 0.0
    @test w1.isPu_RXGB == false && w2.isPu_RXGB == false && w3.isPu_RXGB == false
    @test w1.modelData === tmp1 && w2.modelData === tmp2 && w3.modelData === tmp3

    # Pre-existing tap_side side-selection convention (documented [1,2,3], 0 = no tap) currently
    # never attaches `tap` for tap_side != 0, because the loop reuses/overwrites the `tap` local
    # binding on every non-matching side before the matching side is reached (see docs/dev report
    # for issue #261 Stage: 3WT phase taps). This test snapshots that existing behaviour; it is
    # intentionally not "fixed" here.
    @test isnothing(w1.taps)
    @test isnothing(w2.taps)
    @test isnothing(w3.taps)
  end

  @testset "phase_tap_side attaches PhaseTapChangerModel to exactly one winding" begin
    tmp1, tmp2, tmp3, tapSettings = _3wt_phase_taps_fixtures()
    psc = PhaseTapChangerModel(kind = :asymmetrical, step = 0, lowStep = -8, highStep = 8, neutralStep = 0, winding_connection_angle_deg = 60.0)
    w1, w2, w3 = create3WTWindings!(u_kV = [110.0, 20.0, 10.0], sn_MVA = [100.0, 80.0, 20.0], addEx_Side = [tmp1, tmp2, tmp3], sh_deg = [0.0, 0.0, 0.0], tap_side = 1, tap = tapSettings, phase_tap_side = 2, phase_taps = psc)

    @test isnothing(w1.phase_taps)
    @test isnothing(w3.phase_taps)
    @test !isnothing(w2.phase_taps)

    m = w2.phase_taps
    @test m.kind === :asymmetrical
    @test m.winding_connection_angle_deg == 60.0
    @test m.step == 0
    @test m.lowStep == -8
    @test m.highStep == 8
    @test m.neutralStep == 0
  end

  @testset "ratio tap and phase tap can share the same winding" begin
    tmp1, tmp2, tmp3, tapSettings = _3wt_phase_taps_fixtures()
    psc = PhaseTapChangerModel(kind = :asymmetrical, step = 0, lowStep = -8, highStep = 8, neutralStep = 0, winding_connection_angle_deg = 60.0)
    # tap_side = 0 is the only value for which the pre-existing loop actually preserves `tap`
    # (it matches on the very first iteration, before the destructive overwrite occurs) — see the
    # snapshot test above. Used here only to exercise "ratio tap and phase tap coexist on winding 1".
    w1, w2, w3 = create3WTWindings!(u_kV = [110.0, 20.0, 10.0], sn_MVA = [100.0, 80.0, 20.0], addEx_Side = [tmp1, tmp2, tmp3], sh_deg = [0.0, 0.0, 0.0], tap_side = 0, tap = tapSettings, phase_tap_side = 1, phase_taps = psc)

    @test !isnothing(w1.taps)
    @test !isnothing(w1.phase_taps)
    @test isnothing(w2.taps) && isnothing(w2.phase_taps)
    @test isnothing(w3.taps) && isnothing(w3.phase_taps)
  end

  @testset "create3WTWindings! phase-tap keyword validation" begin
    tmp1, tmp2, tmp3, tapSettings = _3wt_phase_taps_fixtures()
    psc = PhaseTapChangerModel(kind = :asymmetrical, step = 0, lowStep = -8, highStep = 8, neutralStep = 0, winding_connection_angle_deg = 60.0)

    # phase_tap_side set without a model
    @test_throws ArgumentError create3WTWindings!(u_kV = [110.0, 20.0, 10.0], sn_MVA = [100.0, 80.0, 20.0], addEx_Side = [tmp1, tmp2, tmp3], sh_deg = [0.0, 0.0, 0.0], tap_side = 1, tap = tapSettings, phase_tap_side = 2)

    # model set without phase_tap_side
    @test_throws ArgumentError create3WTWindings!(u_kV = [110.0, 20.0, 10.0], sn_MVA = [100.0, 80.0, 20.0], addEx_Side = [tmp1, tmp2, tmp3], sh_deg = [0.0, 0.0, 0.0], tap_side = 1, tap = tapSettings, phase_taps = psc)

    # phase_tap_side out of range
    @test_throws ArgumentError create3WTWindings!(u_kV = [110.0, 20.0, 10.0], sn_MVA = [100.0, 80.0, 20.0], addEx_Side = [tmp1, tmp2, tmp3], sh_deg = [0.0, 0.0, 0.0], tap_side = 1, tap = tapSettings, phase_tap_side = 4, phase_taps = psc)
  end
end

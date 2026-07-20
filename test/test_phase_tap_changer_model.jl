# file: test/test_phase_tap_changer_model.jl

function run_phase_tap_changer_model_tests()
  @testset "PhaseTapChangerModel constructor validation" begin
    @test_throws ArgumentError PhaseTapChangerModel(kind = :bogus, step = 0, lowStep = -5, highStep = 5, neutralStep = 0)

    # Stage 4: :tabular is now implemented (previously rejected with a "Stage 4" hint).
    table = [TapTablePoint(step = s, ratio = 1.0, angle_deg = Float64(s)) for s in -5:5]
    tabular_model = PhaseTapChangerModel(kind = :tabular, step = 0, lowStep = -5, highStep = 5, neutralStep = 0, table = table)
    @test tabular_model.kind === :tabular

    @test_throws ArgumentError PhaseTapChangerModel(kind = :asymmetrical, step = 0, lowStep = -5, highStep = 5, neutralStep = 0)
    # ψ given: no error
    PhaseTapChangerModel(kind = :asymmetrical, step = 0, lowStep = -5, highStep = 5, neutralStep = 0, winding_connection_angle_deg = 90.0)
  end

  @testset "calcPhaseTapAngleRatio :symmetrical" begin
    m = PhaseTapChangerModel(kind = :symmetrical, step = 0, lowStep = -10, highStep = 10, neutralStep = 0, voltage_step_increment = 0.01)

    neutral_result = calcPhaseTapAngleRatio(m; step = 0)
    @test neutral_result.effective_ratio == 1.0
    @test neutral_result.effective_shift_deg == 0.0

    # u = 0.01, n-n0 = 5 -> f = 0.05, alpha = 2*atand(0.025)
    expected_alpha_deg = 2.0 * atand(0.025)
    off_neutral_result = calcPhaseTapAngleRatio(m; step = 5)
    @test isapprox(off_neutral_result.effective_shift_deg, -expected_alpha_deg; atol = 1e-12) # default :reciprocal_from_side negates
    @test off_neutral_result.effective_ratio == 1.0

    for step in (-10, -5, 0, 5, 10)
      @test calcPhaseTapAngleRatio(m; step = step).effective_ratio == 1.0
    end

    m_direct = PhaseTapChangerModel(kind = :symmetrical, step = 0, lowStep = -10, highStep = 10, neutralStep = 0, voltage_step_increment = 0.01, convention = :direct_regulating_vector)
    direct_result = calcPhaseTapAngleRatio(m_direct; step = 5)
    @test isapprox(direct_result.effective_shift_deg, expected_alpha_deg; atol = 1e-12)
  end

  @testset "calcPhaseTapAngleRatio :asymmetrical delegates to calcSkewAngleTap" begin
    u = 0.02
    m_qb = PhaseTapChangerModel(kind = :asymmetrical, step = 4, lowStep = -8, highStep = 8, neutralStep = 0, voltage_step_increment = u, winding_connection_angle_deg = 90.0)
    f_qb = calcPhaseTapFraction(m_qb)
    direct_qb = calcSkewAngleTap(tap_fraction = f_qb, skew_angle_deg = 90.0)
    result_qb = calcPhaseTapAngleRatio(m_qb)
    @test result_qb.effective_ratio == direct_qb.effective_ratio
    @test result_qb.effective_shift_deg == direct_qb.effective_shift_deg
    @test result_qb.regulating_vector == direct_qb.regulating_vector

    m_lon = PhaseTapChangerModel(kind = :asymmetrical, step = 4, lowStep = -8, highStep = 8, neutralStep = 0, voltage_step_increment = u, winding_connection_angle_deg = 0.0)
    f_lon = calcPhaseTapFraction(m_lon)
    direct_lon = calcSkewAngleTap(tap_fraction = f_lon, skew_angle_deg = 0.0)
    result_lon = calcPhaseTapAngleRatio(m_lon)
    @test result_lon.effective_shift_deg == 0.0
    @test result_lon.effective_ratio == direct_lon.effective_ratio
  end

  @testset "calcPhaseTapReactance" begin
    m_sym_no_x = PhaseTapChangerModel(kind = :symmetrical, step = 0, lowStep = -10, highStep = 10, neutralStep = 0, voltage_step_increment = 0.01)
    @test isnothing(calcPhaseTapReactance(m_sym_no_x, 1.0))

    m_sym_x = PhaseTapChangerModel(kind = :symmetrical, step = 0, lowStep = -10, highStep = 10, neutralStep = 0, voltage_step_increment = 0.01, x_min = 0.1, x_max = 0.2)
    alphamax_sym = calcPhaseTapAngleRatio(m_sym_x; step = m_sym_x.highStep).effective_shift_deg
    @test calcPhaseTapReactance(m_sym_x, 0.0) == 0.1
    @test isapprox(calcPhaseTapReactance(m_sym_x, alphamax_sym), 0.2; atol = 1e-12)
    alpha_mid_sym = calcPhaseTapAngleRatio(m_sym_x; step = 5).effective_shift_deg
    expected_mid_sym = 0.1 + (0.2 - 0.1) * (sind(alpha_mid_sym / 2.0) / sind(alphamax_sym / 2.0))^2
    @test isapprox(calcPhaseTapReactance(m_sym_x, alpha_mid_sym), expected_mid_sym; atol = 1e-12)

    m_asym_x = PhaseTapChangerModel(kind = :asymmetrical, step = 0, lowStep = -10, highStep = 10, neutralStep = 0, voltage_step_increment = 0.01, winding_connection_angle_deg = 90.0, x_min = 0.15, x_max = 0.3)
    alphamax_asym = calcPhaseTapAngleRatio(m_asym_x; step = m_asym_x.highStep).effective_shift_deg
    @test calcPhaseTapReactance(m_asym_x, 0.0) == 0.15
    @test isapprox(calcPhaseTapReactance(m_asym_x, alphamax_asym), 0.3; atol = 1e-12)
    alpha_mid_asym = calcPhaseTapAngleRatio(m_asym_x; step = 5).effective_shift_deg
    expected_mid_asym = 0.15 + (0.3 - 0.15) * (tand(alpha_mid_asym) / tand(alphamax_asym))^2
    @test isapprox(calcPhaseTapReactance(m_asym_x, alpha_mid_asym), expected_mid_asym; atol = 1e-12)
  end

  @testset "DTF _dtf_effective_transformer_tap equivalence with pre-migration formula" begin
    nominal_voltages_kv = [110.0]
    from_bus = Sparlectra.DTFImporter.DTFBus("", 1, 1, 1, "PV", 110.0, 0.0, 0.0, 0.0, 10.0, 2.0, -5.0, 5.0)
    to_bus = Sparlectra.DTFImporter.DTFBus("", 2, 2, 1, "SLACK", 110.0, 0.0, 0.0, 0.0, 20.0, 3.0, -10.0, 10.0)
    longitudinal_range_percent = 12.0
    added_voltage_angle_deg = 15.0
    max_tap_step = 8
    actual_tap_step = 3
    control = Sparlectra.DTFImporter.DTFTransformerControl(
      "", 1, "PV", "SLACK", "A", "", "PV", "SLACK",
      110.0, 110.0, longitudinal_range_percent, added_voltage_angle_deg, max_tap_step, actual_tap_step,
      nothing, nothing, nothing,
    )
    case = Sparlectra.DTFImporter.DTFCase(
      "synthetic_phase_tap",
      100.0,
      Sparlectra.DTFImporter.DTFParams("", Float64[]),
      ["synthetic_phase_tap"],
      nominal_voltages_kv,
      Sparlectra.DTFImporter.DTFSize("", 2, 1, 0, 1, "SLACK"),
      [Sparlectra.DTFImporter.DTFBranch("", 1, 'T', 1, "A", "PV", "SLACK", 1.21, 6.05, 0.0, 0.0, nothing)],
      Sparlectra.DTFImporter.DTFCompensation[],
      [control],
      [from_bus, to_bus],
      Sparlectra.DTFImporter.DTFOutage[],
      Sparlectra.DTFImporter.DTFTrailingRecord[],
    )
    branch = only(case.branches)

    result = Sparlectra.DTFImporter._dtf_effective_transformer_tap(case, branch, control, from_bus, to_bus)

    tap_fraction_old = (longitudinal_range_percent / 100.0) * actual_tap_step / max_tap_step
    tap_result_old = calcSkewAngleTap(tap_fraction = tap_fraction_old, skew_angle_deg = added_voltage_angle_deg)

    @test result.model == :skew_angle
    @test result.tap_fraction == tap_fraction_old
    @test result.shift_deg == tap_result_old.effective_shift_deg
    @test result.effective_complex == tap_result_old.regulating_vector
    @test result.ratio == tap_result_old.effective_ratio # base_ratio == 1.0 for the default :neutral_one transformer_ratio_mode
  end

  @testset "PowerTransformerWinding without phase_taps defaults to nothing" begin
    w_kw = PowerTransformerWinding(Vn_kV = 110.0)
    @test isnothing(w_kw.phase_taps)

    w_pos = PowerTransformerWinding(110.0, 0.0, 0.12)
    @test isnothing(w_pos.phase_taps)
  end
end

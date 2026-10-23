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

# file: test/test_phase_tap_table.jl
# purpose: unit tests for the tabular phase-tap model: TapTablePoint and
#          constructor validation, step-range derivation, calcPhaseTapTable
#          lookup, and the tabular reactance and fraction

function run_phase_tap_table_tests()
  @testset "TapTablePoint / PhaseTapChangerModel :tabular constructor validation" begin (function ()
    # empty table
    @test_throws ArgumentError PhaseTapChangerModel(kind = :tabular, step = 0, neutralStep = 0, table = TapTablePoint[])

    # duplicate steps
    dup_table = [TapTablePoint(step = 0, ratio = 1.0, angle_deg = 0.0), TapTablePoint(step = 0, ratio = 1.01, angle_deg = 1.0)]
    @test_throws ArgumentError PhaseTapChangerModel(kind = :tabular, step = 0, neutralStep = 0, table = dup_table)

    # unsorted steps
    unsorted_table = [TapTablePoint(step = 1, ratio = 1.0, angle_deg = 0.0), TapTablePoint(step = 0, ratio = 1.0, angle_deg = 0.0)]
    @test_throws ArgumentError PhaseTapChangerModel(kind = :tabular, step = 0, neutralStep = 0, table = unsorted_table)

    valid_table = [TapTablePoint(step = s, ratio = 1.0 + 0.01 * s, angle_deg = Float64(s)) for s in -3:3]

    # neutralStep not in table
    @test_throws ArgumentError PhaseTapChangerModel(kind = :tabular, step = 0, neutralStep = 99, table = valid_table)

    # table on non-tabular kind
    @test_throws ArgumentError PhaseTapChangerModel(kind = :symmetrical, step = 0, lowStep = -3, highStep = 3, neutralStep = 0, voltage_step_increment = 0.01, table = valid_table)

    # increment/psi/x fields set on :tabular
    @test_throws ArgumentError PhaseTapChangerModel(kind = :tabular, step = 0, neutralStep = 0, table = valid_table, voltage_step_increment = 0.01)
    @test_throws ArgumentError PhaseTapChangerModel(kind = :tabular, step = 0, neutralStep = 0, table = valid_table, winding_connection_angle_deg = 90.0)
    @test_throws ArgumentError PhaseTapChangerModel(kind = :tabular, step = 0, neutralStep = 0, table = valid_table, x_min = 0.1)
    @test_throws ArgumentError PhaseTapChangerModel(kind = :tabular, step = 0, neutralStep = 0, table = valid_table, x_max = 0.2)

    # inconsistent explicit lowStep/highStep
    @test_throws ArgumentError PhaseTapChangerModel(kind = :tabular, step = 0, lowStep = -99, neutralStep = 0, table = valid_table)
    @test_throws ArgumentError PhaseTapChangerModel(kind = :tabular, step = 0, highStep = 99, neutralStep = 0, table = valid_table)
  end)() end

  @testset "lowStep/highStep auto-derivation from table" begin (function ()
    table = [TapTablePoint(step = s, ratio = 1.0, angle_deg = Float64(s)) for s in -4:6]
    m = PhaseTapChangerModel(kind = :tabular, step = 0, neutralStep = 0, table = table)
    @test m.lowStep == -4
    @test m.highStep == 6

    # explicit, consistent lowStep/highStep also works
    m2 = PhaseTapChangerModel(kind = :tabular, step = 0, lowStep = -4, highStep = 6, neutralStep = 0, table = table)
    @test m2.lowStep == -4
    @test m2.highStep == 6
  end)() end

  @testset "calcPhaseTapTable exact lookup" begin (function ()
    table = [
      TapTablePoint(step = -1, ratio = 0.99, angle_deg = -2.0, x_pu = 0.11),
      TapTablePoint(step = 0, ratio = 1.0, angle_deg = 0.0), # x_pu omitted -> nothing
      TapTablePoint(step = 1, ratio = 1.01, angle_deg = 2.0, x_pu = 0.13),
    ]
    m = PhaseTapChangerModel(kind = :tabular, step = 0, neutralStep = 0, table = table)

    r_minus1 = calcPhaseTapTable(m; step = -1)
    @test r_minus1.effective_ratio == 0.99
    @test r_minus1.effective_shift_deg == -2.0
    @test r_minus1.x_pu == 0.11

    r0 = calcPhaseTapTable(m; step = 0)
    @test r0.effective_ratio == 1.0
    @test isnothing(r0.x_pu)

    r_default_step = calcPhaseTapTable(m) # defaults to m.step == 0
    @test r_default_step.effective_ratio == r0.effective_ratio

    @test_throws ArgumentError calcPhaseTapTable(m; step = 99)
    # regulating_vector of a table row is the inverse of the formula path's
    # (ratio, shift) mapping: a formula model read back through its own table
    # gives the formula's vector again
    m_formula = PhaseTapChangerModel(kind = :asymmetrical, step = 3, lowStep = -8, highStep = 8, neutralStep = 0, voltage_step_increment = 0.02, winding_connection_angle_deg = 25.0)
    steps = m_formula.lowStep:m_formula.highStep
    formula = Dict(s => calcPhaseTapAngleRatio(m_formula; step = s) for s in steps)
    m_tab = PhaseTapChangerModel(kind = :tabular, step = 3, neutralStep = 0, table = [TapTablePoint(step = s, ratio = formula[s].effective_ratio, angle_deg = formula[s].effective_shift_deg) for s in steps])
    @test all(isapprox(calcPhaseTapAngleRatio(m_tab; step = s).regulating_vector, formula[s].regulating_vector; atol = 1e-12) for s in steps)
  end)() end

  @testset "calcPhaseTapReactance and calcPhaseTapFraction on :tabular" begin (function ()
    table = [
      TapTablePoint(step = -2, ratio = 1.0, angle_deg = 0.0, x_pu = 0.10),
      TapTablePoint(step = 0, ratio = 1.0, angle_deg = 0.0, x_pu = 0.20),
      TapTablePoint(step = 2, ratio = 1.0, angle_deg = 0.0), # x_pu === nothing
    ]
    m = PhaseTapChangerModel(kind = :tabular, step = -2, neutralStep = 0, table = table)
    @test calcPhaseTapReactance(m, 12345.0) == 0.10 # alpha_deg is ignored; row at m.step is used

    m0 = PhaseTapChangerModel(kind = :tabular, step = 0, neutralStep = 0, table = table)
    @test calcPhaseTapReactance(m0, -9999.0) == 0.20

    m2 = PhaseTapChangerModel(kind = :tabular, step = 2, neutralStep = 0, table = table)
    @test isnothing(calcPhaseTapReactance(m2, 0.0))

    # a table has no step fraction
    @test_throws ArgumentError calcPhaseTapFraction(m0)
  end)() end
end

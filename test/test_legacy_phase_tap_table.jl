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

# file: test/test_legacy_phase_tap_table.jl
# purpose: legacy profile: the testsets of test/test_phase_tap_table.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_phase_tap_table.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_phase_tap_table_tests()
  # cut from the runner function
  @testset "calcPhaseTapAngleRatio :tabular reproduces the formula path bit-identically" begin (function ()
    m_formula = PhaseTapChangerModel(kind = :asymmetrical, step = 3, lowStep = -8, highStep = 8, neutralStep = 0, voltage_step_increment = 0.02, winding_connection_angle_deg = 25.0)
    formula_results = Dict(s => calcPhaseTapAngleRatio(m_formula; step = s) for s in m_formula.lowStep:m_formula.highStep)

    table = [TapTablePoint(step = s, ratio = formula_results[s].effective_ratio, angle_deg = formula_results[s].effective_shift_deg) for s in m_formula.lowStep:m_formula.highStep]
    m_tabular = PhaseTapChangerModel(kind = :tabular, step = 3, neutralStep = 0, table = table)

    for s in m_formula.lowStep:m_formula.highStep
      formula_result = formula_results[s]
      tabular_result = calcPhaseTapAngleRatio(m_tabular; step = s)
      # ratio/shift are stored verbatim in the table row, so these are bit-identical.
      @test tabular_result.effective_ratio == formula_result.effective_ratio
      @test tabular_result.effective_shift_deg == formula_result.effective_shift_deg
      # regulating_vector is reconstructed from (ratio, shift_deg) via deg2rad/cis, a
      # trig round-trip that is not bit-identical at the ULP level; this validates the
      # derivation formula itself, not bit-for-bit equality.
      @test isapprox(tabular_result.regulating_vector, formula_result.regulating_vector; atol = 1e-12)
    end
  end)() end

  return true
end

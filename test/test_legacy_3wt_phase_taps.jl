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

# file: test/test_legacy_3wt_phase_taps.jl
# purpose: legacy profile: the testsets of test/test_3wt_phase_taps.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_3wt_phase_taps.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_3wt_phase_taps_tests()
  # cut from the runner function
  @testset "create3WTWindings! existing docstring example: no phase taps" begin (function ()
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
    # binding on every non-matching side before the matching side is reached (Issue #261,
    # Stage: 3WT phase taps). This test snapshots that existing behaviour; it is
    # intentionally not "fixed" here.
    @test isnothing(w1.taps)
    @test isnothing(w2.taps)
    @test isnothing(w3.taps)
  end)() end

  return true
end

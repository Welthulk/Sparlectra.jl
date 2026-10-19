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

# file: test/extended/test_legacy_dtf_importer.jl
# purpose: legacy profile: the testsets of test/extended/test_dtf_importer.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from extended/test_dtf_importer.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_dtf_importer_tests()
  # cut from the runner function
  @testset "native DTF full local fixture" begin
    if !isfile(DTF_FIXTURE)
      @info "Skipping full FOR001 fixture validation; external reference network is not tracked. Place a local file at test/fixtures/dtf/FOR001.DAT or data/DTF/FOR001.DAT for manual validation."
      @test_skip "FOR001 full reference fixture not available"
      return
    end
    case = Sparlectra.DTFImporter.read_dtf(DTF_FIXTURE)
    @test case.size.NGES == 13
    @test length(case.branches) == 27
    @test count(b -> b.kind == 'T', case.branches) == 5
    net = Sparlectra.createNetFromDTFFile(DTF_FIXTURE)
    @test length(net.branchVec) == 27
    @test Sparlectra.bus_shunt_totals_pu(net).total_g_pu ≈ 0.0
  end

  return true
end

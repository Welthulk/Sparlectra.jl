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

# file: test/test_legacy_measurement_table.jl
# purpose: legacy profile: the testsets of test/test_measurement_table.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_measurement_table.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_measurement_table_tests()
  # cut from the runner function
  @testset "measurement table: bad data on the SmallGrid fixture" begin
    (function ()
      result = _smallgrid_import()
      rows = readMeasurementTable(_SMALLGRID_TABLE)
      s1 = importMeasurementTable(rows, cgmesMeasurementKeys(result)).sets[1]
      marked = Set(_table_row_id(r) for r in rows if r.timestamp == s1.timestamp && occursin("bad_data", r.source))
      @test length(marked) == 4
      applyMeasurementSet!(result.net, s1)
      d = with_state_estimation_config(update_net = false, max_eliminations = 10) do
        runse_diagnostics(result.net, Measurement[m for m in result.net.measurements])
      end
      # the eliminations stop when the chi-square band test passes, and
      # exactly the marked rows are eliminated
      @test d.stop_reason === :consistent
      @test Set(e.id for e in d.eliminations) == marked
    end)()
  end

  return true
end

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

# file: test/test_legacy_matpower_metadata.jl
# purpose: legacy profile: the testsets of test/test_matpower_metadata.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_matpower_metadata.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_matpower_metadata_tests()
  # cut from the runner function
  @testset "Generator costs: import, export, SCF (#471)" begin (function ()
    # the shipped sp_case9 carries a gencost block: model 2 and model 1 P
    # rows plus Q rows, padded to the width of the piecewise linear row
    src = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "sp_case9.m"))
    m0 = Sparlectra.MatpowerIO.read_case(src)
    @test size(m0.gencost) == (6, 10)
    gens(net) = [ps for ps in net.prosumpsVec if Sparlectra.isGenerator(ps)]
    mktempdir() do tmp
      net = Sparlectra.createNetFromMatPowerFile(filename = src)
      @test all(ps -> ps.gencost !== nothing && ps.gencost.q !== nothing, gens(net))
      # MATPOWER round trip: the matrix comes back identical
      f1 = joinpath(tmp, "rt.m")
      writeMatpowerCasefile(net, f1; write_solution = false)
      @test Sparlectra.MatpowerIO.read_case(f1).gencost == m0.gencost
      # SCF keeps the rows (extra.gencost) and writes them back unchanged
      f2 = joinpath(tmp, "rt.scf.json")
      exportSCF(net; file = f2)
      net2 = importSCF(f2)
      @test [(g.gencost.p, g.gencost.q) for g in gens(net2)] == [(g.gencost.p, g.gencost.q) for g in gens(net)]
      f3 = joinpath(tmp, "rt2.m")
      writeMatpowerCasefile(net2, f3; write_solution = false)
      @test Sparlectra.MatpowerIO.read_case(f3).gencost == m0.gencost
      # the converter path (MATPOWER to SCFCase to Net) assigns the same rows
      cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
      net3 = Sparlectra.build_net(Sparlectra.convert_case(MatpowerAdapter(), m0, Sparlectra.matpower_adapter_options(cfg)); config = cfg)
      @test [(g.gencost.p, g.gencost.q) for g in gens(net3)] == [(g.gencost.p, g.gencost.q) for g in gens(net)]
      # a case without gencost writes none
      n118 = Sparlectra.createNetFromMatPowerFile(filename = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "sp_case118.m")))
      f4 = joinpath(tmp, "n118.m")
      writeMatpowerCasefile(n118, f4; write_solution = false)
      m118 = Sparlectra.MatpowerIO.read_case(f4)
      @test m118.gencost === nothing
      # a MATPOWER import without apply_branch_kind stores the kind :unknown
      # on its transformers; the export must still write them as "T" with
      # their ratio, not as lines with TAP 0 (nine transformers, 0.935 to
      # 0.985), and the lines (ratio 0) as "L"
      @test m118.branch[:, 9] == [br.ratio for br in n118.branchVec]
      @test m118.branch_kind == [br.ratio != 0.0 ? "T" : "L" for br in n118.branchVec]
      @test count(==("T"), m118.branch_kind) == length(n118.trafos) == 9
      # a generator added in Sparlectra has no cost row: the block is left
      # out with a warning instead of shifting the rows onto other units
      addProsumer!(net = net, busName = only([k for (k, v) in net.busDict if v == 5]), type = "GENERATOR", p = 5.0, q = 0.0)
      f5 = joinpath(tmp, "added.m")
      @test_logs (:warn, r"mpc.gencost not written") match_mode = :any writeMatpowerCasefile(net, f5; write_solution = false)
      @test Sparlectra.MatpowerIO.read_case(f5).gencost === nothing
      # Q rows on only some units: the P rows are written, the Q rows left out
      net6 = Sparlectra.createNetFromMatPowerFile(filename = src)
      g6 = gens(net6)[2]
      g6.gencost = GenCost(g6.gencost.p)
      f6 = joinpath(tmp, "partial_q.m")
      @test_logs (:warn, r"only some generators carry a Q cost row") match_mode = :any writeMatpowerCasefile(net6, f6; write_solution = false)
      @test Sparlectra.MatpowerIO.read_case(f6).gencost == m0.gencost[1:3, :]
      # a block that cannot be assigned or holds an invalid row is not
      # imported, with a warning (the case itself still imports)
      bad = replace(read(src, String), "\t2\t0\t0\t3\t0.002\t0\t0\t0\t0\t0;\n" => "")
      fb = joinpath(tmp, "bad_rows.m")
      write(fb, bad)
      nb = @test_logs (:warn, r"mpc.gencost has 5 rows for 3 generators") match_mode = :any Sparlectra.createNetFromMatPowerFile(filename = fb)
      @test all(ps -> ps.gencost === nothing, gens(nb))
      # GenCost itself is strict: a row that is not model 1 or 2 is an error
      @test_throws ArgumentError GenCost([3.0, 0.0, 0.0, 2.0, 1.0, 0.0])
      @test_throws ArgumentError GenCost([1.0, 0.0, 0.0, 2.0, 10.0, 400.0])
    end
  end)() end

  return true
end

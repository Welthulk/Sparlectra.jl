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

# file: examples/others/exp_matpower_gencost.jl
# purpose: MATPOWER generator costs (mpc.gencost) on the shipped sp_case9:
#          the rows each generator keeps, the export that writes them back
#          unchanged (also after an SCF round trip), and what happens when a
#          generator added in Sparlectra carries no cost row.

using Sparlectra

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

gens(net::Net) = [ps for ps in net.prosumpsVec if Sparlectra.isGenerator(ps)]

"""
    main() -> NamedTuple

1. import `data/mpower/sp_case9.m` and list the cost rows of each generator
   (P row and Q row, MATPOWER model 1 or 2);
2. export it with `writeMatpowerCasefile` and read the file back: the
   `mpc.gencost` matrix is the same;
3. the same after an SCF round trip (`exportSCF`/`importSCF`);
4. add a generator without costs: the export leaves the block out and says
   why, because a partial block would shift the rows onto other units; give
   the new unit a `GenCost` and the block is written again.
"""
function main()
  print_example_banner("examples/others/exp_matpower_gencost.jl", "MATPOWER generator costs kept through import, export and SCF")
  src = joinpath(pkgdir(Sparlectra), "data", "mpower", "sp_case9.m")
  original = Sparlectra.MatpowerIO.read_case(src).gencost
  net = redirect_stdout(devnull) do
    createNetFromMatPowerFile(filename = src)
  end
  rows = [(name = Sparlectra.getCompName(g.comp), model_p = Int(g.gencost.p[1]), p = g.gencost.p, q = g.gencost.q) for g in gens(net)]

  out = mktempdir()
  roundtrip = joinpath(out, "sp_case9_costs.m")
  redirect_stdout(devnull) do
    writeMatpowerCasefile(net, roundtrip; write_solution = false)
  end
  same_matpower = Sparlectra.MatpowerIO.read_case(roundtrip).gencost == original

  scf = joinpath(out, "sp_case9_costs.scf.json")
  exportSCF(net; file = scf)
  via_scf = importSCF(scf)
  roundtrip2 = joinpath(out, "sp_case9_costs_scf.m")
  redirect_stdout(devnull) do
    writeMatpowerCasefile(via_scf, roundtrip2; write_solution = false)
  end
  same_scf = Sparlectra.MatpowerIO.read_case(roundtrip2).gencost == original

  # a unit added in Sparlectra: no cost row, so the export writes no block
  bus5 = only([k for (k, v) in net.busDict if v == 5])
  addProsumer!(net = net, busName = bus5, type = "GENERATOR", p = 5.0, q = 0.0)
  partial = joinpath(out, "sp_case9_added.m")
  redirect_stdout(devnull) do
    writeMatpowerCasefile(net, partial; write_solution = false)
  end
  block_after_add = Sparlectra.MatpowerIO.read_case(partial).gencost !== nothing

  # give it costs (a linear P cost of 30 per MWh, a zero Q cost) and the block is back
  gens(net)[end].gencost = GenCost([2.0, 0.0, 0.0, 2.0, 30.0, 0.0], [2.0, 0.0, 0.0, 2.0, 0.0, 0.0])
  complete = joinpath(out, "sp_case9_added_costs.m")
  redirect_stdout(devnull) do
    writeMatpowerCasefile(net, complete; write_solution = false)
  end
  rows_after_costs = size(Sparlectra.MatpowerIO.read_case(complete).gencost, 1)
  return (rows = rows, same_matpower = same_matpower, same_scf = same_scf, block_after_add = block_after_add, rows_after_costs = rows_after_costs)
end

result = run_example(main)
println()
println("cost rows per generator (MATPOWER model 1 = piecewise linear, 2 = polynomial):")
for r in result.rows
  println("  ", rpad(r.name, 24), " model ", r.model_p, "  P ", r.p, "  Q ", r.q)
end
println()
println("MATPOWER export, re-read: same gencost matrix = ", result.same_matpower)
println("after an SCF round trip:  same gencost matrix = ", result.same_scf)
println("generator added without costs: block written = ", result.block_after_add, " (a warning names the unit)")
println("after GenCost on the new unit: ", result.rows_after_costs, " cost rows (P and Q for 4 generators)")

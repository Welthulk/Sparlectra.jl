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
# Date: 2026-10-12
# file: examples/others/exp_case_scenarios.jl
# purpose: the scenarios a case file carries, run the way the Web UI button
#          "Scenarios" runs them: a scenarios block written into a copy of the
#          shipped sp_case60 (two hand-written events, a load increase and the
#          generated busbar sections of the block's `events` key), read back
#          from the file and run alone, without the N-1 list. Ends with a
#          ranked table and the report by kind.
#          julia --project=. examples/others/exp_case_scenarios.jl

using Sparlectra
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

const SCF_CASE = joinpath(dirname(dirname(pathof(Sparlectra))), "data", "scf", "sp_case60.scf.json")

# the SCF id of a component named in the file; a name the file uses twice
# (the two legs of Neubach_Waldsee) gives the lower id
id_of(index, name) = minimum(id for (id, n) in index.name_by_id if n == name)
out(index, name) = PatchOp(op = :status, target = index.kind_by_id[id_of(index, name)], id = id_of(index, name), value = 0.0)

function main()
  print_example_banner("examples/others/exp_case_scenarios.jl", "the scenarios block of a case file, run alone (the Web UI Scenarios button)")

  # 1. a copy of the shipped case: the scenarios belong to the file
  workdir = mktempdir()
  path = joinpath(workdir, "sp_case60.scf.json")
  cp(SCF_CASE, path)
  case = read_scf_json(path)
  index = ScenarioIndex(case)

  # 2. the block: events by hand, a load increase, and generated busbar
  #    sections (the block asks for them with its `events` key)
  loads = [id for (id, k) in index.kind_by_id if k === :load]
  block = ScenarioSet(
    scenarios = [
      Scenario(name = "Neubach double line", weight = 2.0, ops = [out(index, "Altau_Neubach_1"), out(index, "Altau_Neubach_2")]),
      Scenario(name = "Neubach coupler and one Waldsee leg", ops = [out(index, "Link_1"), out(index, "Neubach_Waldsee")]),
      Scenario(name = "evening peak +15 %", ops = [PatchOp(op = :scale, target = :load, id = id, factor = 1.15) for id in loads]),
    ],
    events = [:busbars],
  )
  case.sparlectra.scenarios = scenario_set_dict(block)
  write_scf_json(case, path)
  println("written: ", basename(path), " with ", length(block.scenarios), " scenarios and the event class busbars")

  # 3. read the file back and run its block, nothing else
  case = read_scf_json(path)
  net = build_net(case)
  index = ScenarioIndex(case)
  set = scf_case_scenarios(case)
  results = redirect_stdout(() -> runScenarios!(net, set; index = index), devnull)
  println("ran ", length(results), " cases: ", length(set.scenarios), " from the file, ", length(results) - length(set.scenarios), " generated busbar sections")
  println()

  # 4. the worst cases first: failures, then severity, then the lowest voltage
  ranked = sort(results; by = r -> (r.converged ? 1 : 0, -something(r.severity, 0.0), r.min_vm_pu))
  @printf("%-34s %-11s %-9s %9s %11s %9s  %s\n", "case", "kind", "converged", "Vmin[pu]", "loading[%]", "shed[MW]", "out")
  println("-"^110)
  for r in first(ranked, 12)
    @printf("%-34s %-11s %-9s %9s %11s %9.1f  %s\n", first(r.name, 34), String(r.kind), r.converged ? "yes" : "NO",
      isnan(r.min_vm_pu) ? "-" : @sprintf("%.4f", r.min_vm_pu), isnan(r.max_branch_loading_pct) ? "-" : @sprintf("%.1f", r.max_branch_loading_pct),
      r.shed_load_mw, first(join(r.elements, ", "), 40))
  end
  length(ranked) > 12 && println("... ", length(ranked) - 12, " more")
  println()

  # 5. the report: counts by outcome and by kind
  printContingencyReport(buildContingencyReport(results))
  return nothing
end

run_example(main)

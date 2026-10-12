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

# file: docs/lit/warmup/workshop_scenarios.jl
# purpose: compile warm-up of the scenarios workshop
#          (docs/lit/workshop_scenarios.jl): runs every path the chapters
#          use once on the small shipped sp_case5, so the chapters do not
#          stall on first-call compilation. Included by the first code cell
#          of the notebook after it has set up the study cases. Not part of
#          the library.

"""
    warmup()

Run every compile-heavy path of the scenarios workshop once on the 5-bus
case and print how long it took. The output of the paths is discarded.
"""
function warmup()
  println("warm-up: compiles the paths of this workshop once; how long it takes depends on the machine (a Colab session is several times slower than a desktop)")
  ## every path of the chapters once on the 5-bus case, output discarded:
  ## the SCF read, the name map, N-1, events by hand and generated,
  ## the scenarios block with events, pairs, the report and the CSV
  t_paths = @elapsed redirect_stdout(devnull) do
    wdir = mktempdir()
    wpath = joinpath(dirname(dirname(pathof(Sparlectra))), "data", "scf", "sp_case5.scf.json")
    wcase = read_scf_json(wpath)
    wnet = build_net(wcase)
    windex = ScenarioIndex(wcase)
    wnames = contingency_names(wnet, wcase)
    n1 = runContingencies!(wnet, generateN1Branches(wnet; names=wnames); names=wnames)
    line = first(id for (id, k) in windex.kind_by_id if k === :branch)
    gen = first(id for (id, k) in windex.kind_by_id if k === :generator)
    hand = Scenario(name="by hand", ops=[PatchOp(op=:status, target=:branch, id=line, value=0.0), PatchOp(op=:status, target=:generator, id=gen, value=0.0)])
    runScenarios!(wnet, [hand]; index=windex)
    events = generateCorrelatedEvents(wnet; names=wnames, busbars=true, substations=true, links=true, name_groups=true)
    generateDoubleCircuits(wnet; names=wnames)
    generateLinkOutages(wnet; names=wnames)
    generateNameGroups(wnet; names=wnames)
    runContingencies!(wnet, events; names=wnames)
    wcase.sparlectra.scenarios = scenario_set_dict(ScenarioSet(scenarios=[hand], events=[:double_circuits, :busbars]))
    wfile = joinpath(wdir, "sp_case5.scf.json")
    write_scf_json(wcase, wfile)
    runScenarios!(wnet, scf_case_scenarios(read_scf_json(wfile)); index=windex)
    pairs = runContingencies!(wnet, generateBranchPairs(wnet; names=wnames); names=wnames)
    buildContingencyReport(vcat(n1, pairs))
    writeContingencyResultsCSV(joinpath(wdir, "w.csv"), pairs)
  end
  println("warm-up: ", round(t_paths; digits=2), " s (SCF, names, N-1, events, pairs, report); everything warm")
  return nothing
end

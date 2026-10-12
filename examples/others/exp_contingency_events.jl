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
# file: examples/others/exp_contingency_events.jl
# purpose: correlated contingency events and N-2 pairs (#355), route B of the
#          scenarios workshop: the same steps through the API. On the shipped
#          SCF cases sp_case14 and sp_case60: N-1 with the file names, the
#          double circuit by hand and generated, the Neubach coupler, busbar
#          sections and station, the name group across two sections, two
#          generators in one event in both orders, the N-2 pairs with the
#          emergent ones, and the report by kind.
#          julia --project=. examples/others/exp_contingency_events.jl

using Sparlectra

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

const SCF_DIR = joinpath(dirname(dirname(pathof(Sparlectra))), "data", "scf")

# a shipped SCF case with its network, id index and the file names
function study(name)
  case = read_scf_json(joinpath(SCF_DIR, name * ".scf.json"))
  net = build_net(case)
  return (case = case, net = net, index = ScenarioIndex(case), names = contingency_names(net, case))
end

harmful(r) = !r.converged || (r.error !== nothing && occursin("island", r.error)) || !isempty(r.overloads) || !isempty(r.voltage_violations)

function main()
  print_example_banner("examples/others/exp_contingency_events.jl", "correlated events and N-2 pairs on sp_case14 and sp_case60")
  c14 = study("sp_case14")
  c60 = study("sp_case60")

  # 1. N-1 with the names of the file
  n1 = runContingencies!(c14.net, generateN1Branches(c14.net; names = c14.names); names = c14.names)
  println("sp_case14 N-1: ", length(n1), " outages, ", count(r -> r.converged, n1), " converged")

  # 2. the double circuit, by hand (SCF ids) and generated: the same row values
  ids = sort([id for (id, nm) in c14.index.name_by_id if startswith(nm, "Ostheim_Neuwiese_")])
  by_hand = Scenario(name = "both Ostheim-Neuwiese lines", ops = [PatchOp(op = :status, target = :branch, id = id, value = 0.0) for id in ids])
  r_hand = only(runScenarios!(c14.net, [by_hand]; index = c14.index))
  r_gen = only(runContingencies!(c14.net, generateDoubleCircuits(c14.net; names = c14.names); names = c14.names))
  println(r_gen.name, ": vmin ", round(r_gen.min_vm_pu; digits = 4), " pu (by hand ", round(r_hand.min_vm_pu; digits = 4), " pu)")

  # 3. the events the configuration switches on (contingency.events.*)
  cfg = ContingencyEventsConfig(enabled = true, busbars = true)
  events = generateCorrelatedEvents(c14.net; names = c14.names, double_circuits = cfg.double_circuits, busbars = cfg.busbars, max_members = cfg.max_members)
  println("sp_case14 events: ", length(events), " (", count(c -> c.kind === :busbar, events), " busbar sections)")

  # 4. Neubach: the coupler, one section, the whole station
  neubach = vcat(generateLinkOutages(c60.net; names = c60.names),
    [c for c in generateBusbarEvents(c60.net; names = c60.names) if startswith(c.name, "BB Neubach_110")],
    [c for c in generateSubstationEvents(c60.net; names = c60.names) if c.name == "SUB Neubach_110"])
  for r in runContingencies!(c60.net, neubach; names = c60.names)
    println(rpad(r.name, 20), rpad(String(r.kind), 12), "out: ", join(r.elements, ", "), "; ", something(r.error, "nothing cut off"))
  end

  # 5. the name group across two sections
  ng = only(c for c in generateNameGroups(c60.net; names = c60.names) if c.name == "NG Neubach_Waldsee")
  println(ng.name, ": ", join(ng.elements, ", "))

  # 6. two generators in one event: the order of the ops does not matter
  gens = sort([id for (id, k) in c60.index.kind_by_id if k === :generator]; by = id -> c60.index.internal_by_id[id])
  unit(id) = PatchOp(op = :status, target = :generator, id = id, value = 0.0)
  ab, ba = runScenarios!(c60.net, [Scenario(name = "a then b", ops = [unit(gens[1]), unit(gens[2])]), Scenario(name = "b then a", ops = [unit(gens[2]), unit(gens[1])])]; index = c60.index)
  println("two units out: vmin ", round(ab.min_vm_pu; digits = 4), " pu in either order (", ab.min_vm_pu == ba.min_vm_pu, ")")

  # 7. N-2 pairs: harmful although both singles are harmless
  singles = generateN1Branches(c60.net; names = c60.names)
  single_harm = Dict(only(c.elements) => harmful(r) for (c, r) in zip(singles, runContingencies!(c60.net, singles; names = c60.names)))
  pairs = generateBranchPairs(c60.net; names = c60.names)
  results = runContingencies!(c60.net, pairs; names = c60.names)
  emergent = count(((p, r),) -> harmful(r) && !any(e -> single_harm[e], p.elements), zip(pairs, results))
  println("sp_case60 pairs: ", length(pairs), ", harmful ", count(harmful, results), ", harmful with two harmless singles ", emergent)

  # 8. the report by kind
  printContingencyReport(buildContingencyReport(runContingencies!(c60.net, vcat(singles, neubach); names = c60.names)))
  return nothing
end

run_example(main)

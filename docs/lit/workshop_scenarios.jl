# Copyright 2023-2026 Udo Schmitz                                             #src
#                                                                             #src
# Licensed under the Apache License, Version 2.0 (the "License");             #src
# you may not use this file except in compliance with the License.            #src
# You may obtain a copy of the License at                                     #src
#                                                                             #src
#     http://www.apache.org/licenses/LICENSE-2.0                              #src
#                                                                             #src
# Unless required by applicable law or agreed to in writing, software         #src
# distributed under the License is distributed on an "AS IS" BASIS,           #src
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    #src
# See the License for the specific language governing permissions and         #src
# limitations under the License.                                              #src
#                                                                             #src
# file: docs/lit/workshop_scenarios.jl                                        #src
# purpose: Literate.jl source of the scenarios workshop notebook and its      #src
#          Documenter page: outages as data in a Sparlectra Case Format       #src
#          file (route A) and through the API (route B, the example           #src
#          examples/others/exp_contingency_events.jl): N-1, a double circuit  #src
#          by hand and generated, correlated events from the scenarios block  #src
#          and from the configuration, busbar sections and stations, name     #src
#          groups, two generators in one event, N-2 pairs, and reading the    #src
#          results. Every printed number sits next to an @assert, and the     #src
#          test suite RUNS this file, so the notebook cannot drift silently.  #src
#          Regenerate the committed outputs with                              #src
#          `julia --project=docs docs/generate_notebooks.jl`.                 #src

# # Scenarios and correlated events: outages as data
#
# > **Level: Advanced**, builds on the N-1 chapter of the tour.
#
# [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_scenarios.ipynb)
#
# > **Note:** This workshop was created with AI assistance and is reviewed
# > and curated by the maintainer; it is not a fully machine-generated text.
#
# An N-1 study asks which single outage hurts. Real faults often take more
# than one element: both circuits on the same towers, every line at a
# faulted busbar, a whole station. This workshop treats such events as
# data. Route A writes them into the scenarios block of a Sparlectra Case
# Format (SCF) file, so the case carries its study. Route B builds them
# with the generators of the API, the way the example
# `examples/others/exp_contingency_events.jl` does. Both routes run on the
# same engine path as N-1 and give the same rows.
#
# > **Note:** On Google Colab the install cell takes a minute or two on a
# > fresh session (package download and precompilation). Colab's Julia
# > version may change over time; this notebook targets Julia >= 1.12.

#nb # ## Setup (Colab)
#nb # This cell installs Sparlectra from GitHub (branch `main`) into a fresh
#nb # temporary environment. The isolation matters on Colab: the shared
#nb # default environment ships many preinstalled packages, and installing
#nb # anything there triggers precompilation of that whole stack. Run this
#nb # cell first, once per session; it takes a minute or two.
#nb # > **Colab is slow here.** The install cell below fetches Sparlectra and
#nb # > its dependencies and precompiles them on Colab's two cores; that takes a
#nb # > minute or two and prints one line per package as it goes. Sparlectra
#nb # > itself precompiles in well under a minute; the first cell that
#nb # > imports a case and solves it then compiles for about a minute more,
#nb # > every later cell runs at speed. Wait for the install cell to finish
#nb # > before running the next one; a second click on it starts the whole
#nb # > install again.
#nb using Pkg
#nb Pkg.activate(temp = true)
#nb ## no APSLF in this workshop: skip the precompile workload of AnalyticLoadFlow (set before Pkg.add, the precompile processes inherit ENV)
#nb ENV["ANALYTICLOADFLOW_PRECOMPILE_WORKLOAD"] = "off"
#nb Pkg.add(url = "https://github.com/Welthulk/Sparlectra.jl", rev = "main")
#nb ## To test another branch, set rev to its name, e.g. rev = "dev/r0.10.0".
#nb ## For the latest registered release use: Pkg.add("Sparlectra")

# ## Warm-up and the study cases
#
# Julia compiles each function on first use. This cell loads the package
# and sets up what the chapters share, so any chapter can run on its own
# after it; a code cell that uses a name from another cell says so in a
# comment. Its last line, `warmup()`, runs every path the chapters use once
# on the small shipped `sp_case5`, so nothing stalls mid-workshop; its code
# is in `docs/lit/warmup/workshop_scenarios.jl` of the package. How long the
# warm-up takes depends on the machine.
#
# The study cases are two shipped Sparlectra Case Format (SCF) files:
# `sp_case14`, a meshed 110 kV level with one double circuit, and
# `sp_case60`, a 220/110 kV grid with a double ring, a station of two
# busbar sections joined by a coupler, and a line pair that carries the
# same name on both legs. Scenarios are a feature of this format: the
# file carries them, and the result rows use the names written in the file.

using Sparlectra

## the two study cases, each as file path, typed case, network, id index
## and the display names of the file
data = joinpath(dirname(dirname(pathof(Sparlectra))), "data", "scf")
function study(name)
    path = joinpath(data, name * ".scf.json")
    case = read_scf_json(path)
    net = build_net(case)
    return (path=path, case=case, net=net, index=ScenarioIndex(case), names=contingency_names(net, case))
end
c14 = study("sp_case14")
c60 = study("sp_case60")
println("sp_case14: ", length(c14.net.nodeVec), " buses, ", length(c14.net.branchVec), " branches; sp_case60: ", length(c60.net.nodeVec), " buses, ", length(c60.net.branchVec), " branches, ", length(c60.net.linkVec), " coupler")
@assert length(c14.net.branchVec) == 17
@assert length(c60.net.branchVec) == 74
@assert length(c60.net.linkVec) == 1

## compile every path the chapters use once; the warm-up code is in
## docs/lit/warmup/workshop_scenarios.jl of the installed package
include(pkgdir(Sparlectra, "docs", "lit", "warmup", "workshop_scenarios.jl"))
warmup()

# ## 1. The N-1 baseline
#
# One case per in-service branch. `names = c14.names` makes the rows carry
# the line names of the file (`Ostheim_Neuwiese_1`), the names the imported
# network carries as well.

## uses: c14 (warm-up cell)
n1 = runContingencies!(c14.net, generateN1Branches(c14.net; names=c14.names); names=c14.names)
println("N-1 on sp_case14: ", length(n1), " branch outages, ", count(r -> r.converged, n1), " converged, ", count(r -> r.error !== nothing && occursin("island", r.error), n1), " cut a bus off")
@assert length(n1) == 17
@assert count(r -> r.converged, n1) == 14
@assert all(r -> r.kind === :branch && length(r.elements) == 1, n1)

# ## 2. The double circuit: by hand and generated
#
# The two lines between Ostheim and Neuwiese hang on the same towers; one
# fault can take both. Route A writes that event into the file: a scenario
# with two status operations on the SCF component ids. Route B asks
# `generateDoubleCircuits` for it. Both run as one outage of two branches
# and give the same row values.

## uses: c14 (warm-up cell)
ids = sort([id for (id, nm) in c14.index.name_by_id if startswith(nm, "Ostheim_Neuwiese_")])
by_hand = Scenario(name="both Ostheim-Neuwiese lines", ops=[PatchOp(op=:status, target=:branch, id=id, value=0.0) for id in ids])
r_hand = only(runScenarios!(c14.net, [by_hand]; index=c14.index))
dc = only(generateDoubleCircuits(c14.net; names=c14.names))
r_gen = only(runContingencies!(c14.net, [dc]; names=c14.names))
println("sp_case14, both lines out, power flow result:")
println("  case                          kind            lowest bus Vm [pu]   highest branch loading [%]")
for r in (r_hand, r_gen)
    println("  ", rpad(r.name, 30), rpad(String(r.kind), 16), lpad(round(r.min_vm_pu; digits=4), 18), lpad(round(r.max_branch_loading_pct; digits=1), 29))
end
@assert ids == [26, 27]
@assert r_hand.kind === :explicit && r_gen.kind === :double_circuit
@assert r_hand.elements == r_gen.elements == ["Ostheim_Neuwiese_1", "Ostheim_Neuwiese_2"]
@assert r_hand.min_vm_pu == r_gen.min_vm_pu && r_hand.max_branch_loading_pct == r_gen.max_branch_loading_pct

# ## 3. Generated events from the file and from the configuration
#
# A scenarios block can ask for the events instead of listing them: the
# key `events` names the classes. We write it into a copy of the case file,
# read it back and run it. The service does the same for a plain N-1 run
# when `contingency.events.enabled` is set; the class switches are the
# other `contingency.events.*` keys, here as a `ContingencyEventsConfig`.

## uses: c14 (warm-up cell); by_hand (2. The double circuit: by hand and generated)
@isdefined(by_hand) || error("Run the section \"2. The double circuit: by hand and generated\" first: it sets up by_hand.")
workdir = mktempdir()
copy14 = joinpath(workdir, "sp_case14.scf.json")
set = ScenarioSet(scenarios=[by_hand], events=[:double_circuits, :busbars])
c14.case.sparlectra.scenarios = scenario_set_dict(set)
write_scf_json(c14.case, copy14)
back = scf_case_scenarios(read_scf_json(copy14))
from_file = runScenarios!(c14.net, back; index=c14.index)
kinds = Dict(k => count(r -> r.kind === k, from_file) for k in unique(r.kind for r in from_file))
println("scenarios block read back from the file: ", length(from_file), " cases (", join(["$(k) $(n)" for (k, n) in sort(collect(kinds))], ", "), ")")
cfg = ContingencyEventsConfig(enabled=true, double_circuits=true, busbars=true)
from_cfg = generateCorrelatedEvents(c14.net; names=c14.names, double_circuits=cfg.double_circuits, busbars=cfg.busbars, max_members=cfg.max_members)
println("events generated from the configuration (double circuits, busbars): ", length(from_cfg))
@assert back.events == [:double_circuits, :busbars]
@assert kinds[:explicit] == 1 && kinds[:double_circuit] == 1 && kinds[:busbar] == 12
@assert length(from_cfg) == 13

# One of the busbar events, `BB Ostheim_110`, cuts off the slack bus of
# sp_case14. For that event the engine makes another generator the slack
# and says so in the `auto_slack:` line above (`Gen_110` at
# `Dallwitz_110`). The other events keep the original slack.
#
# In the Web UI the button **Scenarios** runs this block (it appears for a
# case file with a scenarios block); the Settings page carries the
# `contingency.events.*` switches in its N-1 box, which add the events to a
# **Contingency (N-1)** run.

# ## 4. Neubach: the coupler, one section, the whole station
#
# Station Neubach on sp_case60 has two 110 kV busbar sections joined by a
# coupler (a bus link). Opening the coupler alone splits the station into
# two buses. A busbar fault on one section opens its lines and the
# coupler; the other section stays. A station fault opens every line at
# both sections. Loads and units stay in the model; what the event cuts
# off is reported (`shed_load_mw`, the note in `error`).

## uses: c60 (warm-up cell)
neubach = vcat(generateLinkOutages(c60.net; names=c60.names),
    [c for c in generateBusbarEvents(c60.net; names=c60.names) if startswith(c.name, "BB Neubach_110")],
    [c for c in generateSubstationEvents(c60.net; names=c60.names) if c.name == "SUB Neubach_110"])
for r in runContingencies!(c60.net, neubach; names=c60.names)
    println(rpad(r.name, 20), rpad(String(r.kind), 12), "out: ", join(r.elements, ", "), "; shed load ", r.shed_load_mw, " MW; ", something(r.error, "nothing cut off"))
end
@assert [c.kind for c in neubach] == [:link, :busbar, :busbar, :substation]
@assert "Link_1" in neubach[2].elements && "Link_1" in neubach[3].elements
@assert !("Link_1" in neubach[4].elements)

# ## 5. A name group across two sections
#
# The line pair `Neubach_Waldsee` carries the same name on both legs and
# starts at the two different sections of Neubach. Its endpoints differ,
# so `generateDoubleCircuits` does not see it; only the names say it is
# one corridor. Two legs with the same name get `#<index>` in the output.

## uses: c60 (warm-up cell)
ng = only(c for c in generateNameGroups(c60.net; names=c60.names) if c.name == "NG Neubach_Waldsee")
println(ng.name, ": ", join(ng.elements, ", "))
@assert length(ng.elements) == 2 && all(e -> startswith(e, "Neubach_Waldsee#"), ng.elements)
@assert !any(c -> Set(c.elements) == Set(ng.elements), generateDoubleCircuits(c60.net; names=c60.names))

# ## 6. Two generators in one event
#
# An event may take several units. The order of the operations does not
# matter: units leave the working copy in descending position order, so
# no position shifts under another one.

## uses: c60 (warm-up cell)
gens = sort([id for (id, k) in c60.index.kind_by_id if k === :generator]; by=id -> c60.index.internal_by_id[id])
unit(id) = PatchOp(op=:status, target=:generator, id=id, value=0.0)
ab, ba = runScenarios!(c60.net, [Scenario(name="a then b", ops=[unit(gens[1]), unit(gens[2])]), Scenario(name="b then a", ops=[unit(gens[2]), unit(gens[1])])]; index=c60.index)
println("sp_case60, units gens[1] and gens[2] out, lowest bus Vm [pu] (power flow result):")
println("  order a then b:   ", round(ab.min_vm_pu; digits=4))
println("  order b then a:   ", round(ba.min_vm_pu; digits=4))
@assert ab.min_vm_pu == ba.min_vm_pu && ab.max_branch_loading_pct == ba.max_branch_loading_pct

# `gens[1]` is the slack unit of sp_case60 (`Gen_220` at `Altau_220`), so
# both scenarios lose the slack. Each run prints one `auto_slack:` line
# naming the unit that takes over, `Gen_110` at `Birkkamp_110` in both
# orders.

# ## 7. N-2 pairs
#
# `generateBranchPairs` lists every pair of closed branches: 2701 on
# sp_case60. A pair is harmful when it fails, cuts a bus off, overloads a
# branch or leaves the voltage band. Most harmful pairs contain a single
# outage that is harmful on its own. The interesting ones are harmful
# although both singles are harmless.

## uses: c60 (warm-up cell)
harmful(r) = !r.converged || (r.error !== nothing && occursin("island", r.error)) || !isempty(r.overloads) || !isempty(r.voltage_violations)
singles = generateN1Branches(c60.net; names=c60.names)
## the island reports of 2701 solves would fill the page; the runs are quiet
single_results = redirect_stdout(() -> runContingencies!(c60.net, singles; names=c60.names), devnull)
single_harm = Dict(only(c.elements) => harmful(r) for (c, r) in zip(singles, single_results))
pairs = generateBranchPairs(c60.net; names=c60.names)
pair_results = redirect_stdout(() -> runContingencies!(c60.net, pairs; names=c60.names), devnull)
emergent = [p for (p, r) in zip(pairs, pair_results) if harmful(r) && !any(e -> single_harm[e], p.elements)]
circuits = Set(Set(c.elements) for c in generateDoubleCircuits(c60.net; names=c60.names))
println("N-2 on sp_case60, number of branch pairs:")
for (label, n) in (("all pairs", length(pairs)), ("harmful", count(harmful, pair_results)),
    ("harmful although both singles are harmless", length(emergent)),
    ("of these, double circuits ($(length(circuits)) in the case)", count(p -> Set(p.elements) in circuits, emergent)))
    println("  ", rpad(label * ":", 46), lpad(n, 5))
end
@assert length(pairs) == 2701
## 0.50.0: harmful old 1106, new 1115; emergent old 124, new 133; demo case rebuilt, transformer magnetizing from nameplate since 0.32.1 #src
@assert count(harmful, pair_results) == 1115
@assert length(emergent) == 133
@assert count(p -> Set(p.elements) in circuits, emergent) == 0

# None of the six double circuits is among them: the ring carries each
# corridor twice, and losing both legs moves the flow to the other half of
# the ring. The 133 are mostly voltage events: two lines whose outages
# each lower the same buses a little, together below 0.9 pu. That is why a
# rule that keeps only pairs with a single near a limit misses some of
# them; the full list is the safe choice on a network of this size.

# ## 8. Reading the results
#
# Every row says what it is (`kind`) and what was out (`elements`, by the
# names of the file). The report counts the cases per kind; the CSV writes
# both columns after the name, the elements joined by `|`.

## uses: c60 (warm-up cell); neubach (4. Neubach: the coupler, one section, the whole station); workdir (3. Generated events from the file and from the configuration)
@isdefined(neubach) || error("Run the section \"4. Neubach: the coupler, one section, the whole station\" first: it sets up neubach.")
@isdefined(workdir) || error("Run the section \"3. Generated events from the file and from the configuration\" first: it sets up workdir.")
report = buildContingencyReport(redirect_stdout(() -> runContingencies!(c60.net, vcat(generateN1Branches(c60.net; names=c60.names), neubach); names=c60.names), devnull))
println("sp_case60 report, cases by kind: ", join(["$(k) $(n)" for (k, n) in report.n_by_kind], ", "))
## format="technical" fixes comma and decimal point; the default follows the regional settings
csv = writeContingencyResultsCSV(joinpath(workdir, "events.csv"), redirect_stdout(() -> runContingencies!(c60.net, neubach; names=c60.names), devnull); format="technical")
println(first(readlines(csv)))
@assert report.n_by_kind == [:branch => 74, :link => 1, :busbar => 2, :substation => 1]
@assert startswith(first(readlines(csv)), "name,kind,elements,")

# The Web UI shows the same: a `kind` column, the elements and the
# overloaded branches in the line below each row, and the count of events
# in the run summary.

# ## 9. SCF only
#
# Scenarios and events belong to the Sparlectra Case Format. A MATPOWER,
# CGMES or PowSyBl case runs a plain N-1 (the service says so in
# `run.log`), and converting an SCF case to one of these formats does not
# carry the scenarios block.
#
# ## Where to go from here
#
# - [N-1 Contingency Analysis](../contingency.md): correlated events, N-2
#   pairs, screening and the display names.
# - [Sparlectra Case Format](../scf.md): the scenarios block with `kind`
#   and `events`.
# - The example `examples/others/exp_contingency_events.jl` runs the same
#   steps as a script.

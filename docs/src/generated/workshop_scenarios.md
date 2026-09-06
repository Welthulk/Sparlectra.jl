```@meta
EditURL = "../../lit/workshop_scenarios.jl"
```

# Scenarios and screening: outages as data

> **Level: Advanced**, builds on the N-1 chapter of the tour.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_scenarios.ipynb)

> **Note:** This workshop was created with AI assistance and is reviewed
> and curated by the maintainer; it is not a fully machine-generated text.

An N-1 study answers "which single outage hurts". This workshop goes one
step further: outages and operating changes become DATA, a scenarios
block inside a Sparlectra Case Format (SCF) file, so a case ships the
study it was meant for. You will expand the classic N-1 list, write
three scenarios by hand (a double outage, a load scaling, a generator
setpoint change), run everything through `runScenarios!`, switch the
contingency screening from `:off` to `:flag` and read what the screened
share does and does not promise, and finally write the block back into
a copy of the case file and read it again.

> **Note:** On Google Colab the install cell takes a few minutes on a
> fresh session (package download and precompilation). Colab's Julia
> version may change over time; this notebook targets Julia >= 1.12.

## Load case118 and turn it into an SCF case

case118 is the IEEE 118-bus system. We fetch the MATPOWER file (the
repository's local copy when present, a download otherwise), build the
network once through the MATPOWER adapter, and export it as an SCF file
into a scratch directory, because THIS workshop edits the case file.

````@example workshop_scenarios
using Sparlectra

workdir = mktempdir()
local_case = joinpath(dirname(dirname(pathof(Sparlectra))), "data", "mpower", "case118.m")
case_m = isfile(local_case) ? local_case : Sparlectra.FetchMatpowerCase.ensure_casefile("case118.m"; outdir = workdir, to_jl = false)
net = createNetFromMatPowerFile(filename = case_m, flatstart = false, enable_pq_gen_controllers = true, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)
scf_path = joinpath(workdir, "case118_workshop.scf.json")
exportSCF(net; file = scf_path, case_name = "case118_workshop", source_reference = "case118.m")
scfcase = Sparlectra.read_scf_json(scf_path)
index = ScenarioIndex(scfcase)
println("buses: ", length(net.nodeVec), ", branches: ", length(net.branchVec))
@assert length(net.nodeVec) == 118
@assert length(net.branchVec) == 186
````

## The N-1 expansion

The N-1 modes are the special case of the scenario model: one status
patch per in-service branch or generator, expanded through the same
generators `runContingencies!` uses.

````@example workshop_scenarios
n1 = expand_scenarios(ScenarioSet(mode = :n1_all), net, index)
n1_branches = expand_scenarios(ScenarioSet(mode = :n1_branches), net, index)
println("N-1 all: ", length(n1), " scenarios (", length(n1_branches), " branch outages, ", length(n1) - length(n1_branches), " generator outages)")
@assert length(n1) == 240
@assert length(n1_branches) == 186
@assert all(length(s.ops) == 1 for s in n1)
````

## Three scenarios by hand

A scenario is a named, weighted, ordered list of patch operations on SCF
component ids. The ids come from the typed case; the `extra` block maps
them to the reference names.

````@example workshop_scenarios
branch_ids = [r.id for r in scfcase.data.line]
gen_id = first(id for (id, k) in index.kind_by_id if k === :generator)
load_ids = [id for (id, k) in index.kind_by_id if k === :load]
handmade = ScenarioSet(scenarios = [
  Scenario(name = "double outage", weight = 2.0, ops = [
    PatchOp(op = :status, target = :branch, id = branch_ids[1], value = 0.0),
    PatchOp(op = :status, target = :branch, id = branch_ids[2], value = 0.0),
  ]),
  Scenario(name = "load up 20 percent", ops = [PatchOp(op = :scale, target = :load, id = load_ids[1], factor = 1.2)]),
  Scenario(name = "gen setpoint 25 MW", ops = [PatchOp(op = :set, target = :generator, id = gen_id, field = :p, value = 25.0)]),
])
validate_scenarios(handmade, index)
results = runScenarios!(net, handmade; index = index)
for r in results
  println(r.name, ": converged = ", r.converged, ", vmin = ", round(r.min_vm_pu; digits = 4), " pu, iterations = ", r.iterations)
end
@assert length(results) == 3
@assert all(r.converged for r in results)
@assert all(isfinite(r.min_vm_pu) for r in results)
````

The multi-op scenario runs the SAME solver ladder and metrics as an N-1
case; the difference is only where the patch list comes from.

## Screening: off against flag

With screening `:off` (the default everywhere) every scenario gets the
full solve, byte for byte the historical result. With `:flag` every
non-islanding single outage is estimated first with one Woodbury-corrected
Newton step on the factorized base Jacobian, and only scenarios whose
estimate comes close to a limit run fully. Screening is a deliberate
opt-in: a linearized step cannot see a bus-type switch, so a class of
generator outages (reactive capability loss) is flagged structurally,
and the honest check before using `:flag` on your own network is exactly
the pair of numbers printed here.

````@example workshop_scenarios
full_set = ScenarioSet(mode = :n1_all)
t_off = @elapsed off_results = runScenarios!(net, full_set; index = index)
t_flag = @elapsed flag_results = runScenarios!(net, full_set; index = index, screening_mode = :flag)
n_screened = count(r -> r.screened, flag_results)
println("off: ", length(off_results), " full solves in ", round(t_off; digits = 2), " s")
println("flag: ", n_screened, " of ", length(flag_results), " screened (", round(100 * n_screened / length(flag_results); digits = 1), " %) in ", round(t_flag; digits = 2), " s")
@assert off_results isa Vector{Sparlectra.ContingencyResult}
@assert flag_results isa Vector{ScenarioResult}
@assert length(off_results) == 240
@assert length(flag_results) == 240
@assert n_screened > 0
````

no violating case may ever be screened away (the acceptance criterion)

````@example workshop_scenarios
violating = Set(r.name for r in off_results if !isempty(r.overloads) || !isempty(r.voltage_violations) || !r.converged)
@assert isempty(intersect(violating, Set(r.name for r in flag_results if r.screened)))
````

## One screened and one flagged case in detail

A screened row carries `start_used = :screen` and the estimate; a
flagged row carries the FULL run plus the estimate that flagged it.

````@example workshop_scenarios
screened_row = first(r for r in flag_results if r.screened)
full_row = only(r for r in off_results if r.name == screened_row.name)
println("screened ", screened_row.name, ": estimated vmin ", round(screened_row.min_vm_pu; digits = 4), " pu against full-run vmin ", round(full_row.min_vm_pu; digits = 4), " pu")
@assert screened_row.start_used === :screen
@assert screened_row.screening_estimate !== nothing
@assert abs(screened_row.min_vm_pu - full_row.min_vm_pu) < 0.02

flagged_row = first(r for r in flag_results if !r.screened && r.screening_estimate !== nothing)
println("flagged ", flagged_row.name, ": full solve ran (start = ", flagged_row.start_used, "), estimate said vmin ", round(flagged_row.screening_estimate.vmin_pu; digits = 4), " pu")
@assert flagged_row.start_used !== :screen
@assert flagged_row.converged || flagged_row.error !== nothing
````

## The scenarios block travels with the case

Writing the block into the SCF file makes the study part of the case:
the Web UI, the service (`scenario_source = file_block`) and the CLI all
read it from there.

````@example workshop_scenarios
scfcase.sparlectra.scenarios = scenario_set_dict(handmade)
Sparlectra.write_scf_json(scfcase, scf_path)
back = scf_case_scenarios(Sparlectra.read_scf_json(scf_path))
println("re-read scenarios: ", join((s.name for s in back.scenarios), ", "))
@assert back !== nothing
@assert [s.name for s in back.scenarios] == [s.name for s in handmade.scenarios]
@assert back.scenarios[1].weight == 2.0
@assert length(back.scenarios[1].ops) == 2
````

## Where to go from here

- [N-1 Contingency Analysis](../contingency.md): the execution model, the
  screening section with the calibration table, and why `:off` is the
  default.
- [Sparlectra Case Format](../scf.md): the scenarios block in the format
  specification.
- The Web UI's scenario editor (`/powerflow/scenarios`) edits exactly the
  block this workshop wrote by hand.


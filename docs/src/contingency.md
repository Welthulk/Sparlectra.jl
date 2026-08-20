# N-1 Contingency Analysis

The contingency batch API evaluates single-branch outages (N-1) on a solved
base case: which line or transformer outage violates voltage limits,
overloads other branches, or splits the network. It is the third parallel
execution surface of the multi-core work (after island solving and
short-circuit sweeps): every contingency is an independent solve on its own
copy of the network, so the batch scales with Julia threads.

## Quick start

```julia
using Sparlectra

net = createNetFromMatPowerFile(filename = "case14.m")
cases = generateN1Branches(net)                    # one case per in-service branch
results = runContingencies!(net, cases)            # parallel when Julia has threads
printContingencyResults(results; max_rows = 20)
writeContingencyResultsCSV("n1.csv", results)
```

Start Julia with `julia --threads=auto` to use all cores; with one thread
the batch runs serially and produces identical results. The runnable
showcase is `examples/others/exp_contingency_n1.jl` (case1354pegase, all
1991 branches, serial vs parallel side by side).

## Execution model

- The base `net` is NEVER mutated. `runContingencies!` first solves a
  template copy of the base case; every contingency then works on a
  `deepcopy` of that template (the copy cost is measured in the analysis
  notes; a structural-sharing copy is a known follow-up), removes its
  branch via `removeBranch!`, re-marks isolated buses, and solves.
- Warm start: because the template carries the solved base voltages, every
  contingency solve starts from the base operating point. When the base
  case itself does not converge, the batch falls back to flat starts and
  prints one warning (on large imported cases a flat start frequently
  diverges, so treat that warning as "fix the base case first").
- Failures are REPORTED, never thrown: a case that does not converge, that
  splits off an island without any reference or promotable generator
  (`error = "islanded without reference"`), or whose element cannot be
  resolved, is returned as a `ContingencyResult` with `converged = false`
  and an actionable `error` line. An island WITH a PV generator is promoted
  to its own reference by the `matpower_like` policy and solves normally;
  `island_count` reports the post-outage island count either way.
- Parallel circuits: two branches between the same buses can share one
  component name. `generateN1Branches` disambiguates them as
  `"<name>#<branchIdx>"` so each circuit gets its own outage case.
- Parallelism: the case list fans out over `runtime.parallel.max_tasks`
  chunks of Julia tasks (gated by `runtime.parallel.enabled` and
  `min_work_items`; the `parallel_*` keywords override the active
  configuration per call). Results are returned in input order and are
  identical to the serial run.

## Evaluation per case

- Voltage band: buses outside `[vm_min_pu, vm_max_pu]` (defaults 0.9/1.1)
  are listed in `voltage_violations`; the envelope is reported as
  `min_vm_pu`/`max_vm_pu` over the non-isolated buses.
- Branch loading: for every branch with a finite `sn_MVA` rating (MATPOWER
  `RATE_A`, CGMES `ratedS`), loading is `100 · max(|S_from|, |S_to|) /
  sn_MVA`; the worst value lands in `max_branch_loading_pct` and branches
  above 100 percent in `overloaded_branches`. Without any rated branch the
  loading is `NaN` — no rating model is invented.
- Case sources: `generateN1Branches(net; include_transformers = true)`
  enumerates the in-service branches (transformers recognizable by
  component type or a nonzero winding ratio), and
  `generateContingenciesFromFOR001(net)` maps imported MATPOWER FOR001
  contingency names; unresolvable names become failed result rows instead
  of disappearing.

## API

```@docs
ContingencyCase
ContingencyResult
runContingencies!
generateN1Branches
generateContingenciesFromFOR001
printContingencyResults
writeContingencyResultsCSV
```

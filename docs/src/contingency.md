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
  `deepcopy` of that template (about 9 ms per case on case1354pegase),
  removes its branch via `removeBranch!`, re-marks isolated buses, and
  solves.
- Warm start: because the template carries the solved base voltages, every
  contingency solve starts from the base operating point. When the base
  case itself does not converge, it is retried through the full solver
  rescue ladder (`runpf!` with `rescue = true`) before the batch falls back
  to flat starts with one warning (on large imported cases a flat start
  frequently diverges, so treat that warning as "fix the base case first").
- Per-case start-value ladder (`contingency.rescue_ladder`, issue #331): an
  ordered, duplicate-free subset of `(:warm, :apslf, :dc, :flat)`, default
  `[:warm]`. Each stage is ONE bounded `runpf!` on the case-local net with a
  distinct start recipe, tried in order until one converges; the stage that
  won is reported as `start_used`:
  - `:warm` starts from the base-case (template) voltages;
  - `:apslf` uses an APSLF start (needs `using AnalyticLoadFlow`; the stage
    is dropped with a warning when the extension is not loaded);
  - `:dc` uses flat magnitudes with DC-projected start angles;
  - `:flat` forces `flatstart = true`.
  The heavier solver-config variants (`:settled_qlimits`, `:autodamp`) are
  used only for the base case above, not per case. `retry_flat_start` is a
  deprecated alias for appending `:flat` to the ladder (kept one minor
  cycle). Measured on case1354pegase (1991 N-1 branch outages, imported with
  the canonical PEGASE convention `shift_unit = :rad`, `shift_sign = -1.0`,
  `ratio = :normal`, and run with `maxIte = 80` and
  `trust_region_enabled = true`, everything else left at the solver defaults):
  the ladder converges 1936 of 1991, the warm start carrying almost all of
  them. Every stage shares that `maxIte`, but note that only `:warm`, `:flat`,
  and `:dc` forward the extra `runpf!` solver keywords (`trust_region_enabled`
  and the like); the `:apslf` stage runs through a config-driven solve without
  them. Of the remaining 55, 53 outages split off a load-only island (no
  generator or source that could serve as an angle reference), reported as
  `islanded without reference`: that is the correct topological outcome, a
  load-shed event of up to about 240 MW, not a solver failure (see the
  load-only-island note below). Only 2 outages, `B_ACL_380_3145_2918` and
  `B_2WT_380_3145_7770`, are genuine non-convergences: each stays connected
  (`island_count = 1`) but hits the iteration limit on every ladder stage, so
  they are ill-conditioned cases, not start-value problems. The ladder only
  ever improves convergence; it cannot give a reference to a load-only island,
  so it is most useful where a warm start is not already dominant.
- Failures are REPORTED, never thrown: a case that does not converge, that
  splits off an island without any reference or promotable generator
  (`error = "islanded without reference"`), or whose element cannot be
  resolved, is returned as a `ContingencyResult` with `converged = false`
  and an actionable `error` line. An island WITH a Slack bus or a PV generator
  is promoted to its own reference by the `matpower_like` policy and solves
  normally; `island_count` reports the post-outage island count either way.
- Load-only islands are a result, not a defect. An island that carries only
  loads (or only fixed-injection PQ generation) has no voltage or frequency
  regulation and therefore no valid angle reference. Declaring an arbitrary PQ
  bus its slack would give a mathematically valid but physically meaningless
  solution, so the honest N-1 statement, matching MATPOWER's convention that an
  island needs a REF or PV bus, names the cause: a load-only outage reports
  `error = "islanded: load-only, X MW load disconnected"`, while an island that
  strands generation without a voltage-controlled source reports
  `"islanded without reference: X MW load, Y MW generation stranded ..."`. On
  case1354pegase all 53 reference-less outages are the load-only case (no
  stranded generation), so no start recipe or slack-search could rescue them.
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
  above 100 percent in `overloads` as [`OverloadRecord`](@ref)s (name, loading,
  the base-case loading and the `delta_pct` to it, `s_MVA`, `sn_MVA`, worst
  first). The delta is the readable part: 105 percent after an outage is a
  different event when the branch sat at 40 percent before than at 98 percent.
  Base loadings come from the solved base case, computed once and reused for
  every contingency. Without any rated branch the loading is `NaN`, and no
  rating model is invented.
- Severity: each result carries `severity = weight · max(0, max loading - 100)`,
  `NaN` for a failed case. `printContingencyResults` ranks by it (failures
  first, then the heaviest weighted overloads), so on a 1991-row list the worst
  contingency is on the first page instead of lost in input order. This is what
  gives the per-case `weight` its consumer.
- Case sources: `generateN1Branches(net; include_transformers = true)`
  enumerates the in-service branches (transformers recognizable by
  component type or a nonzero winding ratio), and
  `generateContingenciesFromFOR001(net)` maps imported MATPOWER FOR001
  contingency names; unresolvable names become failed result rows instead
  of disappearing.
- Screening filters keep the case list to the outages worth simulating on a
  large grid: `generateN1Branches(net; min_vn_kV, min_sn_MVA, name_pattern)`.
  `min_vn_kV` keeps a branch when the higher of its two endpoint voltages
  clears the threshold (so an EHV/HV transformer stays in), `min_sn_MVA`
  keeps only branches carrying a rating at or above the threshold, and
  `name_pattern` keeps names matching a `Regex` or containing a substring.
  All default to "no filter".
- Case weights: each [`ContingencyCase`](@ref) carries a `weight` (default
  `1.0`) for a probability- or severity-weighted ranking; it is carried
  unchanged into the [`ContingencyResult`](@ref) and shown in the table and
  CSV. Attach outage rates in bulk with `applyContingencyWeights(cases,
  weights)`, reading the name-to-weight map from a two-column CSV via
  `readContingencyWeightsCSV`.

## Generator outages

`generateN1Generators(net; min_pg_MW, name_pattern)` builds one
`kind = :gen` case per in-service generator (any injection: generator,
external-grid feed-in, or synchronous machine), filtered optionally by output
`|Pg|` or name. A generator outage removes ONLY that unit's injection, so the
topology is unchanged and the lost active power must be picked up elsewhere:

- By default the slack bus absorbs it. Pass `distributed_slack_enabled = true`
  to share the loss over the surviving participants instead (the keyword flows
  through to `runpf!`; it needs a surviving reference and does not itself
  supply one).
- Removing a bus's last voltage-regulating unit demotes it to PQ. Removing the
  system's ONLY slack leaves it reference-less and is reported as
  `no slack bus registered` (the honest N-1 finding that the slack unit is
  critical). For meaningful generator N-1 on a real grid, pass
  `auto_slack = true` so the solver promotes the strongest surviving generator
  to slack, mirroring real frequency control.
- When a generator outage strands a SEPARATE island that still carries
  injection but no reference (its only regulating unit is the one removed), the
  result names it: `islanded without reference: X MW load, Y MW generation
  stranded`. This is the same reporting the load-only note above describes, now
  reached from the generation side.

## Reporting

Per-case output is the fixed-width `printContingencyResults` table (severity
ranked by default, `sort_by = :none` for input order) and the
`writeContingencyResultsCSV` file (input order, carrying loading, severity, the
overloads, shed load, and weight). For a batch overview,
`buildContingencyReport(results)` folds the results into a
[`ContingencyReport`](@ref): counts by outcome (converged, islanded,
non-converged, with overload, with voltage violation), the total and worst load
shed, the single worst branch loading, the worst weighted severity, and the
branches overloaded by the most contingencies. `printContingencyReport` prints
it as a compact summary block. The structured fields
(`overloads::Vector{OverloadRecord}` with base and delta, `shed_load_mw`,
`severity`) let a consumer rank and total without parsing the message strings.

## Web UI

The "Contingency (N-1)" button on the power-flow run form runs the batch on the
loaded case, with a selector for branch or generator outages (the "generator"
option submits the run parameter `contingency_kind = gen`, the same `:gen` kind
`generateN1Generators` uses). It reuses the same
service, run directory, and case cache as a normal power-flow run (a mode flag
on the single run endpoint, not a separate workflow), so the case is resolved
and built through the shared config-driven import path; `rescue_ladder` is read
from `contingency.rescue_ladder`. The run writes `contingency_n1.csv` and a
`run.log` report and shows an outcome summary. A generator outage that removes
the system's only slack is named as the expected N-1 finding, not a tool
failure; rerun the case with `auto_slack = true` if you want the solver to
promote a surviving generator instead.

### Weights in the Web UI

The "edit N-1 weights" link next to the outage-kind selector opens a per-case
weights editor. Weights live next to the case as
`<case-stem>.contingency-weights.csv`, the exact two-column format
`readContingencyWeightsCSV` parses; the file travels with the case (it is hidden
from the case list and deleted with the case). The editor seeds a table with the
case's real element names (so you never type them blind), offers a raw-CSV text
area, and a file upload. Uploading REPLACES an existing weight file (a weight
list is a working document, unlike a case import), and the upload is validated
before it is stored, so a malformed CSV is rejected with the line number and the
old file is left untouched. Rows left at exactly `1.0` are omitted when the table
is saved, so the file stays a diff of the default. A run picks the weights up
automatically when the file exists, with no extra option: its presence is the
switch. Names in the file that match no element are reported in the `run.log`
(a weight list can outlive a case edit), never fatal.

A weight only reorders the severity ranking; it never skips a case (to leave an
outage out, filter it at generation time, not with a weight of `0`). Grid data
(MATPOWER, CGMES) carries no outage rates, so these weights are always
user-supplied: typical sources are branch failure rates (per 100 km and year)
or an importance rating by voltage level.

## API

```@docs
ContingencyCase
ContingencyResult
OverloadRecord
ContingencyReport
runContingencies!
generateN1Branches
generateN1Generators
generateContingenciesFromFOR001
applyContingencyWeights
readContingencyWeightsCSV
printContingencyResults
writeContingencyResultsCSV
buildContingencyReport
printContingencyReport
```

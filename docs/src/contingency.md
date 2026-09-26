# N-1 Contingency Analysis

The contingency batch API evaluates outages or scenarios on a solved base
case: which line, transformer or generator outage, or which combined patch
scenario, violates voltage limits, overloads other branches, or splits the
network. Case lists come from the N-1 generators (`generateN1Branches`,
`generateN1Generators`), from imported FOR001 metadata, or from a case
file's scenarios block (`runScenarios!`). The batch fans out over Julia
threads, one reused working copy per chunk.

## Quick start

```julia
using Sparlectra

net = createNetFromMatPowerFile(filename = "case14.m")
cases = generateN1Branches(net)                    # one case per in-service branch
results = runContingencies!(net, cases)            # parallel when Julia has threads
printContingencyResults(results; max_rows = 20)
writeContingencyResultsCSV("n1.csv", results)
```

Start Julia with `julia --threads=auto` to use all cores; one thread runs
serially with identical results. Screening is off by default (see below).
The runnable showcase is `examples/others/exp_contingency_n1.jl`
(case1354pegase, all branches, serial and parallel side by side); its base
case already carries overloads, so it demonstrates throughput, not N-1
evaluation on an operated grid.

## Execution model

- The base `net` is never mutated: `runContingencies!` solves a template
  copy of the base case, and every case is evaluated on a working copy
  that is reset to the template state after each case.
- A branch outage takes the branch's own charging arms with it (a
  transformer's magnetizing admittance sits on the branch since 0.20.0,
  see [Branch model](branchmodel.md)); a bus shunt part that a MATPOWER
  reimport recorded for the branch (`mpc.sparlectra.branch_shunts`) is
  removed from the bus shunt for that case as well.
- Warm start from the base solution: every case starts from the solved
  base operating point. A base case that does not converge is retried
  through the solver rescue ladder (`runpf!` with `rescue = true`) before
  the batch falls back to flat starts with one warning; treat that warning
  as "fix the base case first", a flat start on a large imported case
  frequently diverges. Per case, `contingency.rescue_ladder` (an ordered,
  duplicate-free subset of `(:warm, :apslf, :dc, :flat)`, default
  `[:warm]`) runs one bounded `runpf!` per stage until one converges; the
  winning stage is reported as `start_used`.
- Per-chunk working copies: the case list fans out over
  `runtime.parallel.max_tasks` chunks (gated by `runtime.parallel.enabled`
  and `min_work_items`; the `parallel_*` keywords override the
  configuration per call), each chunk on one reused working copy.
- Results are reproducible bitwise: they come back in input order,
  identical to the serial run.

| Ladder stage | Start values | Solver keywords |
|---|---|---|
| `:warm` | the template (base-case) voltages | forwards the extra `runpf!` keywords (`trust_region_enabled` and the like) |
| `:apslf` | an APSLF start (AnalyticLoadFlow.jl is a required dependency, so the stage is always available) | a config-driven solve without the extra keywords |
| `:dc` | flat magnitudes with DC-projected start angles | forwards the extra keywords |
| `:flat` | `flatstart = true` | forwards the extra keywords; `retry_flat_start` is a deprecated alias for appending `:flat` |

All stages share `maxIte`. The solver-config variants (`:settled_qlimits`,
`:autodamp`) apply to the base case only, and the ladder only improves
convergence: it cannot give a reference to a load-only island.

**Notes**

- Failures are reported, never thrown: a case that does not converge, that
  splits off an island without any reference or promotable generator
  (`error = "islanded without reference"`), or whose element cannot be
  resolved comes back as a [`ContingencyResult`](@ref) with
  `converged = false` and an `error` line. An island with a Slack bus or a
  PV generator is promoted to its own reference by the `matpower_like`
  policy and solves normally; `island_count` reports the post-outage
  island count either way.
- Load-only islands are a result, not a defect: an island with only loads
  (or only fixed-injection PQ generation) has no voltage or frequency
  regulation and no valid angle reference, and an arbitrary PQ slack would
  be physically meaningless. As in MATPOWER (an island needs a REF or PV
  bus), such an outage reports
  `error = "islanded: load-only, X MW load disconnected"`, and an island
  that strands generation without a voltage-controlled source reports
  `"islanded without reference: X MW load, Y MW generation stranded ..."`.
- Parallel circuits: two branches between the same buses can share one
  component name; `generateN1Branches` disambiguates them as
  `"<name>#<branchIdx>"`.
- A scenario file or block addresses components by SCF ids; the SCF
  converter supplies that index while the electrics remain the imported
  net.

## Evaluation per case

Every case is checked against the voltage band and the branch ratings of
the base case; a case weight only reorders the severity ranking and never
skips a case (filter at generation time instead of using a weight of `0`).
Grid data carries no outage rates, so weights are user-supplied: branch
failure rates (per 100 km and year) or an importance rating by voltage
level.

| Criterion | Keyword or config key | Result field |
|---|---|---|
| Voltage band | `vm_min_pu`, `vm_max_pu` (defaults 0.9/1.1) | `voltage_violations` (buses outside the band); `min_vm_pu`, `max_vm_pu` (envelope over the non-isolated buses) |
| Branch loading | `sn_MVA` rating of the branch (MATPOWER `RATE_A`, CGMES `ratedS`); loading is `100 · max(\|S_from\|, \|S_to\|) / sn_MVA` | `max_branch_loading_pct` (worst value; `NaN` without any rated branch, no rating model is invented); `overloads`, branches above 100 percent as [`OverloadRecord`](@ref)s (name, loading, base-case loading and the `delta_pct` to it, `s_MVA`, `sn_MVA`, worst first) |
| Severity | `severity = weight · max(0, max loading - 100)` | `severity` (`NaN` for a failed case; `printContingencyResults` ranks by it, failures first) |
| Convergence | `maxIte`, `contingency.rescue_ladder` | `converged`, `iterations`, `start_used` |
| Islands | reference policy `matpower_like` | `island_count`, `shed_load_mw`, `error` |
| Weight | `weight` on [`ContingencyCase`](@ref) (default `1.0`); `applyContingencyWeights(cases, weights)` attaches outage rates in bulk, `readContingencyWeightsCSV` reads the name-to-weight map from a two-column CSV | `weight`, carried unchanged into the result, the table and the CSV |

| Case source | Call | Selection |
|---|---|---|
| in-service branches | `generateN1Branches(net; include_transformers = true)` | transformers recognized by component type or a nonzero winding ratio; filters `min_vn_kV` (kept when the higher endpoint voltage clears the threshold, so an EHV/HV transformer stays in), `min_sn_MVA` (rated at or above), `name_pattern` (a `Regex` or a substring); all default to no filter |
| generators | `generateN1Generators(net; min_pg_MW, name_pattern)` | one `kind = :gen` case per in-service generator, external-grid feed-in or synchronous machine, filtered by `\|Pg\|` or name |
| FOR001 metadata | `generateContingenciesFromFOR001(net)` | imported MATPOWER FOR001 contingency names; unresolvable names become failed result rows |
| scenarios block or JSON | `runScenarios!` | the per-scenario `weight` of the block; a weight file applies to case-list runs only (the outage-kind selector and the `n1_*` scenario sources) |

## Screening

With `screening_mode = :flag` (keyword on `runContingencies!` and
`runScenarios!`, or `contingency.screening.mode` for the service and the
Web UI) every non-islanding single outage is first estimated with one
Woodbury-corrected Newton step on the factorized base Jacobian; only
scenarios whose estimate comes within `screening.margin_pct` (default 10)
of a limit get the full solve. Switch it on only after checking, on your
own network, the screening share and that the flagged margins carry the
violations (the reported `screened` count and reason breakdown).

| Mode | Effect | Result |
|---|---|---|
| `:off` (default) | every case gets the full solve | results and CSV are the full-solve output |
| `:flag` | estimate first, full solve only for flagged cases | screened rows carry `start_used = :screen`, `screened = true` and the estimate; the CSV appends the `screened` and `screening_estimate` columns |
| `:only` | the estimates without full runs | same columns |

**Notes**

- Always solved in full: a generator outage at a bus with regulating
  generation (PV or already Q-clamped), bridge outages, slack units, and
  scenarios whose one-step residual stays above the 0.005 pu trust gate or
  does not shrink.
- Never screened: patch scenarios and distributed-slack last-participant
  outages.

!!! details "Why it is built this way"
    A linearized step cannot see a bus-type switch. On case300 the outage
    of a generator with a zero schedule left a one-step residual of 4e-13
    and a bus minimum at 0.93 pu, while the full solve dropped it to
    0.87 pu: the bus lost its reactive capability and fell to PQ. No
    residual gate sees that class, so the engine flags it structurally
    (regulating generation at the outaged bus always gets the full solve).

    The 0.005 pu trust gate comes from the widest grown network in the
    calibration (case300, five networks and 1478 N-1 cases with margin 10,
    zero false negatives), not from an average; on a base case that
    already sits at its limits the screening share is honestly zero.

## Generator outages

A generator outage removes only the unit's injection; the lost active
power is picked up elsewhere:

- By default the slack bus absorbs it. `distributed_slack_enabled = true`
  shares the loss over the surviving participants (forwarded to `runpf!`;
  it needs a surviving reference and does not supply one).
- Removing a bus's last voltage-regulating unit demotes it to PQ. Removing
  the only slack is reported as `no slack bus registered`: the N-1 answer
  that this unit is critical. For generator N-1 on a real grid, pass
  `auto_slack = true` so the solver promotes the strongest surviving
  generator to slack, mirroring frequency control.
- A separate island that keeps injection but loses its only regulating unit
  reports `islanded without reference: X MW load, Y MW generation stranded`.

## Reporting

| Output | Call | Content |
|---|---|---|
| table | `printContingencyResults` | fixed-width, severity ranked by default, `sort_by = :none` for input order |
| CSV | `writeContingencyResultsCSV` | input order, with loading, severity, overloads, shed load and weight |
| summary | `buildContingencyReport(results)`, `printContingencyReport` | a [`ContingencyReport`](@ref): counts by outcome (converged, islanded, non-converged, with overload, with voltage violation), total and worst load shed, the worst branch loading, the worst weighted severity, and the branches overloaded by the most contingencies |

The structured fields (`overloads::Vector{OverloadRecord}` with base and
delta, `shed_load_mw`, `severity`) let a consumer rank and total without
parsing message strings.

## Web UI

The "Contingency (N-1)" button on the power-flow run form runs the batch
on the loaded case, with a selector for branch or generator outages (the
"generator" option submits `contingency_kind = gen`). It reuses the
service, run directory and case cache of a normal run, takes
`rescue_ladder` from `contingency.rescue_ladder`, writes
`contingency_n1.csv` and a `run.log` report, and shows an outcome summary
that names a generator outage removing the only slack (rerun with
`auto_slack = true`). Per-case weights are edited next to the outage-kind
selector and stored beside the case as
`<case-stem>.contingency-weights.csv`; see [Web UI](webui.md).

## API

The contingency types and functions are documented on the
[Contingency reference page](reference_contingency.md).

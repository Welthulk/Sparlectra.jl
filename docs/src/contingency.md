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
  through the solver rescue ladder (the strategies of `runpf!` with
  `rescue = true`, on the same solver keywords as the base solve, so
  distributed slack and the Q-limit settings hold there too) before
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
| `:apslf` | an APSLF start (AnalyticLoadFlow.jl is a required dependency, so the stage is always available) | forwards the extra keywords |
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
  has no voltage or frequency regulation and no valid angle reference, and
  an arbitrary PQ slack would be physically meaningless. Such an outage
  reports `error = "islanded: load-only, X MW load disconnected"`. An
  island that keeps a generating unit always finds itself a reference: the
  bus of its best voltage-controlled unit first, otherwise of its best
  generating unit (a fixed-injection PQ unit included), independent of
  `auto_slack`; the row names the bus that took over. Only an island whose generation is static compensation
  alone reports `"islanded without reference: ... generation stranded"`.
- Lost reference: `runContingencies!`, `runScenarios!`, the service and
  the Web UI solve with `auto_slack` (the default). When an outage removes
  the reference (the slack unit, or
  the branch that ties it to an island), the best remaining unit takes
  over by the order every reference choice uses: the stated reference
  priority first (`referencePriority`, 1 is the strongest, see
  [Reference priority](slack_vs_source.md#Reference-priority)), then
  [`reference_candidate_rank`](@ref) (an external network injection first,
  then a unit that regulates the voltage of its own bus, then the size);
  an island without a voltage-controlled unit takes its best generating
  unit. Setting priorities gives the replacement of the slack a
  user-controlled order. The row names the bus ("reference taken over by bus
  ..."). Only an island without any generating unit stays without a
  reference. `auto_slack = false` ends such a case on the missing
  reference instead.
- Slack model: the service passes `power_flow.distributed_slack` of the run
  configuration to every post-outage solve. An outage that takes away the
  path of the reference unit's output can be without a solution on a
  single slack and have one when the units share the mismatch (line 1-2 of
  the IEEE 14-bus case).
- Cut-off buses: the outage of a radial branch leaves its far bus without
  a connection. The remaining network solves (`converged = true`), and the
  result carries `error = "islanded: bus B cut off, X MW load disconnected,
  ..."` with `shed_load_mw` set to that load; a unit at such a bus is named
  as out of service. The case counts as islanded in the summary.
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

## [Screening](@id contingency_screening)

With `screening_mode = :flag` (keyword on `runContingencies!` and
`runScenarios!`, or `contingency.screening.mode` for the service and the
Web UI) every non-islanding single outage is first estimated with one
Woodbury-corrected Newton step on the factorized base Jacobian; only
scenarios whose estimate comes close to a limit get the full solve. For a
voltage, close means within `screening.margin_pct` (default 10) percent of
the band width; for a branch, its estimated loading plus the change the
outage causes on it (at least 1 point, at most the margin) reaches 100
percent, because the error of the estimate grows with that change. A
branch the base case already loads near its limit therefore flags only the
outages that move it, not every outage of the network. Switch it on only after checking, on your
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

## [Warm active set](@id contingency_warm_active_set)

The base case of an N-1 batch is solved with Q limits, and on a large
network many machines end clamped at a limit. Without the warm start
every outage starts from the base voltages but with the file's bus
types, so it clamps the same machines again (a jump of the mismatch and
about twice the Newton steps). The warm active set
(`contingency.warm_active_set`, default `true` since 0.30.3; keyword
`warm_active_set` on `runContingencies!` and `runScenarios!`; Web UI: the
checkbox in the "N-1 screening" fieldset of the Settings page) starts every
outage and scenario with the base case's clamped machines as PQ at the
limit they reached. The active set stays free: a clamped machine whose
voltage recovers goes back to PV, new clamps follow the normal rules.

With reactive limits a power flow can have two valid solutions: the same
outage, solved from different starting bus types, can end with a machine
clamped at its limit in one solution and holding its voltage in the other,
both within every Q limit. The warm start can therefore reach a different
solution than the cold one, and on the measured cases it was the more
favourable one. The cold check is an option for that case
(`contingency.warm_cold_check`, default `false`; keyword `warm_cold_check`;
Web UI: "Cold check of tight outages" next to the margin field): when on, a
warm result whose lowest voltage comes within
`contingency.warm_cold_check_margin_pu` (default `0.02` pu) of the lower
voltage limit, or that has an overload or a voltage violation, is
solved a second time from the file's PV/PQ state, and the less favourable of
the two counts (more violations, then the lower lowest voltage, then the
higher loading). Off, every outage reports its warm result.

| Use | Effect |
|---|---|
| default | warm active set `true`, cold check `false` |
| applies when | the base case converged with Q limits on; otherwise one line in the run log says why it was not applied |
| result | the same where the limited solution is unique; where it is not, the warm one (with the cold check on, an outage near a limit or with a violation is checked cold) |
| an outage that does not converge from the warm state | is solved again at once from the file's PV/PQ state; `start_used` is `warm_cold` and the note says so |
| with the cold check on, an outage it judges less favourable cold | the cold result counts; `start_used` is `warm_cold_check` and the note names the warm result |

On case_ACTIVSg2000 (full branch N-1, 3206 outages) the warm start needs
about half the Newton steps (8.0 to 4.2 per outage) with the same
convergence and the same violation list. On five outages of unit
transformers the two starts end on different solutions: seven machines hold
their voltage warm and stay clamped cold. Cold, they sit at their lower Q
limit with the voltage below their setpoint, which is not a consistent
limited state (a machine clamping twice is not released again, #475); the
warm solution is consistent and is reported by default. With the cold check
on, these five report the cold result, the less favourable of the two. On
case300 one outage converged only from the file's state and another only from the warm
state: near the limit of solvability the start decides which limited
solution, if any, Newton reaches. The retry above keeps every
outage the plain start solves.

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

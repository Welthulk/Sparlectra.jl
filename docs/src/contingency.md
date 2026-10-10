# N-1 Contingency Analysis

The contingency batch API evaluates outages or scenarios on a solved base
case: which line, transformer or generator outage, or which patch
scenario, violates voltage limits, overloads other branches or splits the
network. Case lists come from `generateN1Branches` and
`generateN1Generators`, from imported FOR001 metadata, or from a case
file's scenarios block (`runScenarios!`).

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
serially with identical results. Screening is off by default. The
runnable showcase is `examples/others/exp_contingency_n1.jl` (the shipped
`sp_case1354.m`, all branches, serial and parallel); its base case already
carries overloads, so it shows throughput, not N-1 evaluation on an
operated grid.

## Execution model

The base `net` is never mutated: the batch solves a template copy, and
every case runs on a working copy reset to the template after each case.
The case list fans out over `runtime.parallel.max_tasks` chunks (gated by
`runtime.parallel.enabled` and `min_work_items`, overridable per call
with the `parallel_*` keywords), one working copy per chunk; results
come back in input order, bitwise identical to the serial run. A branch
outage takes the branch's own charging arms with it
([Branch model](branchmodel.md)) and a bus shunt part a MATPOWER reimport
recorded for the branch (`mpc.sparlectra.branch_shunts`).

Every case starts from the solved base operating point. A base case that
does not converge is retried through the solver rescue ladder (`runpf!`
with `rescue = true`) before the batch falls back to flat starts with one
warning: fix the base case first. Per case, `contingency.rescue_ladder`
(an ordered subset of `(:warm, :apslf, :dc, :flat)`, default `[:warm]`)
runs one bounded `runpf!` per stage until one converges; the winning
stage is reported as `start_used`. All stages share `maxIte` and forward
the extra `runpf!` keywords; the solver-config variants
(`:settled_qlimits`, `:autodamp`) apply to the base case only, and the
ladder cannot give a reference to a load-only island.

| Ladder stage | Start values |
|---|---|
| `:warm` | the template (base-case) voltages |
| `:apslf` | an APSLF start (AnalyticLoadFlow.jl is a required dependency, so the stage is always available) |
| `:dc` | flat magnitudes with DC-projected start angles |
| `:flat` | `flatstart = true`; `retry_flat_start` is a deprecated alias for appending `:flat` |

Failures are reported, never thrown: a case that does not converge, that
islands without any reference or promotable generator, or whose element
cannot be resolved comes back as a [`ContingencyResult`](@ref) with
`converged = false` and an `error` line; `island_count` reports the
post-outage island count either way.

- A load-only island has no regulation and no angle reference:
  `error = "islanded: load-only, X MW load disconnected"`. An island
  whose generation is static compensation alone reports
  `"islanded without reference: ... generation stranded"`.
- An island that keeps a generating unit takes its own reference, the bus
  of its best voltage-controlled unit, otherwise of its best generating
  unit, in the order every reference choice uses
  ([Reference priority](slack_vs_source.md#Reference-priority), then
  [`reference_candidate_rank`](@ref)); the row names the bus. The batch
  entries, the service and the Web UI solve with `auto_slack` (the
  default); `auto_slack = false` ends a case that lost the whole-network
  reference on the missing reference.
- The outage of a radial branch cuts off its far bus: the rest solves
  (`converged = true`), the result carries
  `error = "islanded: bus B cut off, X MW load disconnected, ..."` with
  `shed_load_mw`, and the case counts as islanded.
- The service passes `power_flow.distributed_slack` to every post-outage
  solve; an outage that cuts the path of the reference unit's output can
  have a solution only with distributed slack (line 1-2 of the IEEE
  14-bus case).
- Parallel circuits sharing one name are disambiguated as
  `"<name>#<branchIdx>"`; a scenario file or block addresses components
  by SCF ids.

## Evaluation per case

Every case is checked against the voltage band and the branch ratings of
the base case; a weight only reorders the severity ranking and never
skips a case (filter at generation time instead). Grid data carries no
outage rates, so weights are user-supplied: failure rates per 100 km and
year, or an importance rating by voltage level.

| Criterion | Keyword or config key | Result field |
|---|---|---|
| Voltage band | `vm_min_pu`, `vm_max_pu` (defaults 0.9/1.1) | `voltage_violations`; `min_vm_pu`, `max_vm_pu` over the non-isolated buses (a three-winding star point never counts) |
| Branch loading | `sn_MVA` rating (MATPOWER `RATE_A`, CGMES `ratedS`); loading is `100 · max(\|S_from\|, \|S_to\|) / sn_MVA` | `max_branch_loading_pct` (`NaN` without any rated branch); `overloads` as [`OverloadRecord`](@ref)s (name, loading, base-case loading, `delta_pct`, `s_MVA`, `sn_MVA`, worst first) |
| Severity | `severity = weight · max(0, max loading - 100)` | `severity` (`NaN` for a failed case; the table ranks by it, failures first) |
| Convergence | `maxIte`, `contingency.rescue_ladder` | `converged`, `iterations`, `start_used` |
| Islands | reference policy `matpower_like` | `island_count`, `shed_load_mw`, `error` |
| Weight | `weight` on [`ContingencyCase`](@ref) (default `1.0`); `applyContingencyWeights(cases, weights)` attaches rates in bulk, `readContingencyWeightsCSV` reads a two-column CSV | `weight`, carried into result, table and CSV |

| Case source | Call | Selection |
|---|---|---|
| in-service branches | `generateN1Branches(net; include_transformers = true)` | transformers by component type or a nonzero winding ratio; filters `min_vn_kV` (the higher endpoint counts), `min_sn_MVA`, `name_pattern` (`Regex` or substring); a three-winding transformer is one `kind = :transformer3w` case ([Three-winding transformers](@ref contingency_three_winding)) |
| generators | `generateN1Generators(net; min_pg_MW, name_pattern)` | one `kind = :gen` case per in-service generator, external-grid feed-in or synchronous machine |
| FOR001 metadata | `generateContingenciesFromFOR001(net)` | imported FOR001 contingency names; unresolvable names become failed rows |
| scenarios block or JSON | `runScenarios!` | the per-scenario `weight` of the block; a weight file applies to case-list runs only |

## [Three-winding transformers](@id contingency_three_winding)

A three-winding transformer is a star equivalent: three legs meeting in
an auxiliary star-point bus. It trips as a whole, so `generateN1Branches`
lists it as one case of kind `:transformer3w`, named after the
transformer, that opens all legs; single-winding outages need
`three_winding_legs = true`. The `n1_branches` scenario mode expands it to
one scenario with a status patch per leg, and a case-file contingency that
lists exactly the legs of one transformer runs as that transformer.

The star point is a computational node, not a bus: it never enters
`min_vm_pu`, `max_vm_pu` or `voltage_violations`, and with the transformer
out it is dropped silently (no island, no lost load). The bus tables show
it with the type `AUX`.

## [Iteration limits](@id contingency_iteration_limits)

The base case of a batch is solved with `power_flow.max_iter`, every
outage and scenario with `contingency.max_iter` (default `30`; keywords
`maxIte` and `base_maxIte` of `runContingencies!` and `runScenarios!`;
Web UI: **Outage max iterations** in the N-1 contingency box of the
Settings page). A lower outage limit keeps a large batch fast: an outage
that needs more steps ends as not converged instead of iterating on.

## [Screening](@id contingency_screening)

With `screening_mode = :flag` (keyword on `runContingencies!` and
`runScenarios!`, or `contingency.screening.mode` for the service and the
Web UI) every non-islanding single outage is first estimated with one
Woodbury-corrected Newton step on the factorized base Jacobian; only
cases whose estimate comes close to a limit get the full solve: a voltage
within `screening.margin_pct` (default 10) percent of the band width, or
a branch whose estimated loading plus the change the outage causes on it
(at least 1 point, at most the margin) reaches 100 percent, because the
error of the estimate grows with that change. Switch screening on only
after checking, on your own network, the screening share and that the
flagged margins carry the violations (the reported `screened` count and
reason breakdown).

| Mode | Effect | Result |
|---|---|---|
| `:off` (default) | every case gets the full solve | results and CSV are the full-solve output |
| `:flag` | estimate first, full solve only for flagged cases | screened rows carry `start_used = :screen`, `screened = true` and the estimate; the CSV appends the `screened` and `screening_estimate` columns |
| `:only` | the estimates without full runs | same columns |

Always solved in full: a generator outage at a bus with regulating
generation (PV or already Q-clamped), bridge outages, slack units, and
cases whose one-step residual stays above the 0.005 pu trust gate or does
not shrink. Never screened: patch scenarios and distributed-slack
last-participant outages.

!!! details "Why it is built this way"
    A linearized step cannot see a bus-type switch: on case300 a generator
    outage left a one-step residual of 4e-13 and a bus minimum at 0.93 pu
    while the full solve dropped it to 0.87 pu, the bus having fallen to
    PQ; hence regulating generation at the outaged bus always gets the
    full solve. The 0.005 pu trust gate is the widest value of a
    calibration over five networks and 1478 N-1 cases (margin 10, zero
    false negatives).

## [Warm active set](@id contingency_warm_active_set)

The base case of an N-1 batch is solved with Q limits, and on a large
network many machines end clamped. Without the warm start every outage
starts from the base voltages but with the file's bus types and clamps
the same machines again (a jump of the mismatch, about twice the Newton
steps). The warm active set (`contingency.warm_active_set`, default
`true`; keyword `warm_active_set` on `runContingencies!` and
`runScenarios!`; Web UI: the checkbox in the N-1 contingency box of the
Settings page) starts every outage and scenario with the base case's
clamped machines as PQ at their limit. The active set stays free: a
machine whose voltage recovers goes back to PV, new clamps follow the
normal rules.

With reactive limits a power flow can have two valid solutions, and the
same outage solved from different starting bus types can end with a
machine clamped in one and holding its voltage in the other; the warm
start can therefore reach a different solution than the cold one, on the
measured cases the more favourable one. The cold check
(`contingency.warm_cold_check`, default `false`; keyword
`warm_cold_check`; Web UI: **Cold check of tight outages** next to the
margin field) solves a warm result whose lowest voltage comes within
`contingency.warm_cold_check_margin_pu` (default `0.02` pu) of the lower
limit, or that has an overload or a voltage violation, a second time from
the file's PV/PQ state, and the less favourable result counts (more
violations, then the lower lowest voltage, then the higher loading).

| Use | Effect |
|---|---|
| default | warm active set `true`, cold check `false` |
| applies when | the base case converged with Q limits on; otherwise one line in the run log says why not |
| result | the same where the limited solution is unique; where it is not, the warm one (with the cold check on, an outage near a limit or with a violation is checked cold) |
| an outage that does not converge from the warm state | is solved again at once from the file's PV/PQ state; `start_used` is `warm_cold` and the note says so |
| with the cold check on, an outage it judges less favourable cold | the cold result counts; `start_used` is `warm_cold_check` and the note names the warm result |

On case_ACTIVSg2000 (3206 outages) the warm start needs fewer Newton
steps per outage with the same convergence and violation list. Machines that clamped twice in early iterations get one
joint release at the converged point, undone for the group if one of them
hits its limit again, so both starts agree where they used to differ. On
case300 one outage converged only from the file's state and another only
from the warm state: near the limit of solvability the start decides
which limited solution Newton reaches, and the retry above keeps every
outage the plain start solves.

## Generator outages

A generator outage removes only the unit's injection. The lost active
power goes to the slack bus, or `distributed_slack_enabled = true` shares
it over the surviving participants (it needs a surviving reference and
does not supply one). Removing a bus's last voltage-regulating unit
demotes it to PQ. Removing the only slack is reported as
`no slack bus registered`, the N-1 answer that this unit is critical;
with `auto_slack = true` the strongest surviving generator is promoted to
slack, mirroring frequency control. An island that keeps injection but
loses its only regulating unit reports
`islanded without reference: X MW load, Y MW generation stranded`.

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

The run form runs the batch on the loaded case (outage kind, scenario
source, screening, weights editor), writes `contingency_n1.csv` and a
`run.log` report, and names a generator outage that removes the only
slack in its summary (rerun with `auto_slack = true`). Controls:
[Scenarios and N-1](webui_reference.md#Scenarios-and-N-1).

## API

The contingency types and functions are documented on the
[Contingency reference page](reference_contingency.md).

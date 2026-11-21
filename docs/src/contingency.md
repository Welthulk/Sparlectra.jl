# N-1 Contingency Analysis

The contingency batch API evaluates outages or scenarios on a solved base
case: which line, transformer or generator outage, which group of elements
that trip together, or which patch scenario violates voltage limits,
overloads other branches or splits the network. Case lists come from
`generateN1Branches` and `generateN1Generators`, from the event generators
([Correlated events](@ref contingency_events)), from N-2 pairs
([N-2 pairs](@ref contingency_pairs)), from imported FOR001 metadata, or
from a case file's scenarios block (`runScenarios!`).

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
operated grid. `examples/others/exp_contingency_events.jl` runs correlated
events and N-2 pairs on the shipped SCF cases.

## Execution model

The base `net` is never mutated: the batch solves a template copy, and
every case runs on a working copy reset to the template after each case.
On networks from 2000 buses the cases run in worker processes of the
batch pool, below that on threads (`runtime.parallel.executor` or the
keyword `executor`, see [Batches in worker processes](@ref parallel_execution_processes)).
On threads the case list fans out over `runtime.parallel.max_tasks`
chunks, one working copy per chunk; in processes every worker holds one
working copy and takes the next case when it is free. Both are gated by
`runtime.parallel.enabled` and `min_work_items`, overridable per call
with the `parallel_*` keywords; results come back in input order, bitwise
identical to the serial run. A branch
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
| `:flat` | `flatstart = true` |

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
  `shed_load_mw`, and the case counts as islanded. Buses joined by closed
  links count as one node, as in the power flow: a bus that hangs on B
  through a closed link alone is cut off with B, and in service as long
  as B is.
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
one scenario with a status patch per leg. The result row lists the legs in
`elements`.

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
steps per outage with the same convergence and violation list. Machines that clamped twice in early iterations get a
joint release at the converged point; one that hits its limit again stays
clamped and the group is released again without it, so both starts agree
where they used to differ. On
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

## [Correlated events](@id contingency_events)

Some elements trip together: two circuits on the same towers, every
branch at a faulted busbar, a whole station. Such an event is one case
with several elements, solved in one run on the same path as N-1. Five
generators build them from the topology; each returns a list of
[`ContingencyCase`](@ref)s:

| Class | Generator | Members | Name |
|---|---|---|---|
| `:double_circuit` | [`generateDoubleCircuits`](@ref) | closed branches with the same two endpoints | `DC <bus>-<bus>` |
| `:busbar` | [`generateBusbarEvents`](@ref) | one bus: its closed branches and its closed couplers; the other sections of the station stay | `BB <bus>` |
| `:substation` | [`generateSubstationEvents`](@ref) | a station (buses joined by closed links): every closed branch at any of its buses | `SUB <bus>` |
| `:link` | [`generateLinkOutages`](@ref) | one closed coupler | the link name |
| `:name_group` | [`generateNameGroups`](@ref) | branches whose file names differ only in a circuit number (`_1`, `_2`), with the same voltage and the same endpoint stations | `NG <stem>` |

[`generateCorrelatedEvents`](@ref) runs the selected classes and keeps
every member set once, as the first class in the order of the table
(double circuit, name group, busbar, station, link): a double circuit that
is also a name group stays a double circuit, a station with one bus is a
busbar event. Events with more than `max_members` elements are skipped; an
event of one branch is no event (that is the N-1 case). The legs of a
three-winding transformer are never members.

A busbar or station event opens branches and couplers only. Loads, units
and shunts at the buses stay in the model; the run reports what the event
cuts off (`shed_load_mw`, stranded units, or an island with its own
reference), as for a radial outage. Opening a coupler splits the station:
the two sections run as separate buses.

Name groups need the file names of an SCF case. A MATPOWER, CGMES or
PowSyBl case keeps its generated names, so `generateNameGroups` returns an
empty list there.

```julia
case = read_scf_json(joinpath(pkgdir(Sparlectra), "data", "scf", "sp_case60.scf.json"))
net = build_net(case)
names = contingency_names(net, case)                 # the file names
cases = vcat(generateN1Branches(net; names = names),
             generateCorrelatedEvents(net; names = names, busbars = true, links = true, name_groups = true))
results = runContingencies!(net, cases; names = names)
```

The service adds the events to a branch N-1 of an SCF case (the
`n1_all` and `n1_branches` sources and the branch outage run) when
`contingency.events.enabled` is set; the class switches and the size limit
are the other `contingency.events.*` keys
([Configuration](configuration.md)). A generator-only run and the other
formats run plain N-1, and `run.log` says why. A scenarios block asks for
events with its `events` key ([Sparlectra Case Format](scf.md)).
Scenarios, and with them the events, are a feature of the SCF case: a
conversion to MATPOWER, CGMES or PowSyBl does not carry them.

## [N-2 pairs](@id contingency_pairs)

[`generateBranchPairs`](@ref) lists every unordered pair of closed branches
as a case of kind `:pair`, named `P <a> + <b>`. The count grows with the
square of the branch count: 74 branches give 2701 pairs, 217 give 23436.
More than `max_pairs` pairs (default 50000) throw before anything is
solved; the run never truncates silently. In the service,
`contingency.pairs.enabled` adds the pairs to a branch N-1 of an SCF case;
a pair with the same member set as an event is listed once, as the event.

A rule that skips pairs which cannot be harmful would make this usable on
larger networks. None ships: the candidates tested so far missed harmful
pairs on the shipped demo cases, so `contingency.pairs.rule` knows `all`
only.

## [Display names](@id contingency_names)

Contingency output names elements by [`contingency_names`](@ref): for an
SCF case the file names (`sparlectra.extra[<id>].name`), for every other
case the generated component names. A name shared by several elements
(parallel circuits, two legs named alike) carries `#<index>`. The case
names, the `elements` of a row, the overloaded branches, the report and the
Web UI table use these names; the network's own component names are not
changed, so the power flow, state estimation and every export keep the
generated names. A weight file saved under the generated names before
0.40.0 still applies: the run maps the old names to the cases and says so
in `run.log`, and the weights editor writes the file names.

## Reporting

| Output | Call | Content |
|---|---|---|
| table | `printContingencyResults` | fixed-width, severity ranked by default, `sort_by = :none` for input order |
| CSV | `writeContingencyResultsCSV` | input order: `name`, `kind`, `elements` (joined by `\|`), weight, convergence, voltages, loading, severity, overloads, voltage violations, islands, shed load, error |
| summary | `buildContingencyReport(results)`, `printContingencyReport` | a [`ContingencyReport`](@ref): counts by outcome (converged, islanded, non-converged, with overload, with voltage violation) and by kind, total and worst load shed, the worst branch loading, the worst weighted severity, and the branches overloaded by the most contingencies |

The structured fields (`overloads::Vector{OverloadRecord}` with base and
delta, `shed_load_mw`, `severity`) let a consumer rank and total without
parsing message strings.

## [Injection batch](@id contingency_injection_batch)

[`runInjectionBatch`](@ref) solves many power flows of one network that
differ only in the load: N scenarios, each with a load factor per bus
([`InjectionScenarios`](@ref)). [`random_load_scenarios`](@ref) draws the
factors uniformly from a band around 1 with a seed; [`load_ramp_scenarios`](@ref)
gives one factor for all buses that rises over N steps, a stand-in for a
day profile of 96 quarter hours. Every scenario starts from the start
iterate of the net, or with `start = :base` from the solved base case (flat
start and start helpers off, as in the QSTS steps; the base case is solved
on the CPU for both backends); the net itself is not changed.

```julia
net = importSCF("data/scf/sp_case60.scf.json")
sc = random_load_scenarios(net, 32; band = 0.2, seed = 1)
r = runInjectionBatch(net, sc)
count(r.converged), r.iterations, r.qlimit_violations
```

The batch runs with Q limits off and without the rescue ladder: a scenario
that does not converge is reported, not retried. Each scenario reports
how many PV buses would leave their reactive limits by more than
`qlimit_tol_pu` (default `1e-4` pu), so a caller can rerun those with the
limits on.

On the CPU (`backend = :cpu`) every scenario is a `runpf!` on a working
copy of the net, reset to the net after each scenario. On networks from
2000 buses the scenarios run in worker processes of the batch pool, one
working copy per worker, each worker taking the next scenario when it is
free; below that the copies are reused within a chunk and the chunks run
on Julia threads (`executor` and `max_processes`, `runtime.parallel`, the
`parallel_*` keywords, see [Batches in worker processes](@ref parallel_execution_processes)).
BLAS runs with one thread inside the batch (`blas_threads`), in the main
process and in every worker. The result of a scenario is the same on
every executor and equals the single `runpf!` with the same loads and
configuration bit for bit if that `runpf!` also runs with BLAS at one
thread. With more BLAS threads UMFPACK can round differently, and on
larger networks the single `runpf!` then differs from the batch in the
last bits. With `power_flow.power_mode` on, every working copy solves the
base case once before its first scenario and keeps its Ybus and LU
analysis across its scenarios (a chunk on threads, all scenarios of a
worker process); every copy starts without the power-mode cache of the
net.

On the GPU (`backend = :gpu`, after `using CUDA, CUDSS`) the Ybus, the
Jacobian pattern and one cuDSS analysis are built once, and all scenarios
iterate side by side in device memory. `gpu_strategy = :sequential`
(default) refactorizes and solves the scenarios one after another with
the one analysis in cuDSS's deterministic mode, so a repeated batch gives
the same bits. `:nonuniform` and `:uniform` use the cuDSS batch modes; the
uniform batch exists only in the non-deterministic mode in cuDSS 1.0.0 and
skips the converged scenarios unless `gpu_uniform_mask = false`. The
precision follows `power_flow.compute_backend.gpu_precision`.
`gpu_two_phase` runs the scenarios in two phases instead: Float32 on the
device down to a mismatch of `gpu_switch_mismatch` (default `1e-5`) or
until the steps stall, then Float64 steps with the last Float32 LU
(`:g5`; `:g5r` refactorizes or goes to Float64 when a step reduces the
mismatch by less than half) or a Float64 Newton on the CPU (`:g6`).
A network the GPU path does not cover (links, more than one island,
Q(U) or P(U) controllers with a voltage-dependent characteristic,
distributed slack, damping or start helpers), or a device error such as
running out of GPU memory, runs the batch on the CPU, and the result's
`note` names the reason. With `power_flow.compute_backend.gpu_fallback_to_cpu
= false` a refusal throws and a device error is passed on.

## Web UI

**Contingency (N-1)** runs the batch on the loaded case (outage kind,
N-1 list, screening, weights editor), writes `contingency_n1.csv` and a
`run.log` report, and names a generator outage that removes the only
slack in its summary (rerun with `auto_slack = true`). For a case file
with a scenarios block, **Scenarios** runs that block alone. Controls:
[Scenarios and N-1](webui_reference.md#Scenarios-and-N-1).

## API

The contingency types and functions are documented on the
[Contingency reference page](reference_contingency.md).

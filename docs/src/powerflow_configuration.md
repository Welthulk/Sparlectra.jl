# Power-Flow Configuration

## Start-source concepts

| Start source | Meaning | Requires previous solution? | Typical use | Risk |
|---|---|---:|---|---|
| Flat start | Synthetic `1.0 pu / 0°` style initialization | no | Small/simple cases | Can converge poorly on stressed grids |
| DC angle start | Synthetic DC-like angle estimate | no | Large transmission-like grids | Ignores full AC effects |
| MATPOWER imported reference | `BUS.VM/VA`, `GEN.VG`, imported setpoint logic | maybe | Validation against imported cases | Can resemble using known answer |
| Historical values | Previous trusted operating point | yes | Repeated operations | Stale operating point |
| SCADA/state-estimation values | Measured/estimated state | yes | Operations workflows | Measurement noise/outliers |
| Profile values | Explicit configured profile values | maybe | Reproducible studies | May not match topology changes |

## [Solver core options](@id pf-solver-core)

The options every AC run reads: formulation and mode, flat start, tolerance and iteration cap, damping, slack promotion, the rescue ladder and the DC fallback. The Web UI form carries the ones a run is usually tuned with; the rest is YAML only.

| YAML path | Type | Default | Allowed values | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `power_flow.method` | Symbol/String | `rectangular` | `rectangular` | AC solver formulation. | Always (current core). | N/A | Fixed single implementation. | Must match `benchmark.methods`. |
| `power_flow.mode` | Symbol | `manual` | `manual`, `auto` | `manual` keeps every option exactly as configured. `auto` inspects the imported network (size, AC islands, R/X profile, phase shifters, start profile, generator Q limits) and fills the start-value, step-control, and Q-limit strategy for keys the user did not set explicitly; explicit keys (overrides, or file values differing from the template default) always win. On non-convergence a bounded escalation ladder retries (step-control switch, guarded pre-solve, full start projection, Q-limit relax/mode switch, APSLF-seeded NR, optional DC fallback via `power_flow.dc.fallback`); the tolerance is never changed. Decisions, conflicts, attempts, and hints land in the `auto_mode_decision.log` artifact and the result metadata (`auto_profile`, `auto_final_stage`, `auto_final_solver`, `auto_hints`). | Integrations that should just work without reading this page (see [Integration Guide](integration.md)). | Strict reproducibility studies where every option must be pinned. | Feature extraction is one cheap pass; escalations only run on failed solves. | The solver-internal `power_flow.rescue` ladder stays active inside every auto attempt; the escalation stages add the strategies it does not cover. |
| `power_flow.flatstart` | Bool | `false` | `true`, `false` | Flat start (1.0 pu, 0 degrees); `false` keeps imported start voltages (MATPOWER `VM`/`VA`, CGMES `SvVoltage`, SCF `start_state`). The Settings page offers it as the one start switch: while it is set, a run switches `apslf_start`, `dc_seed_unconditional` and `start_current_iteration` off and treats both start modes as `classic` (the stored values stay, `run.log` names the overrides). On CGMES runs an explicit `cgmes_import.start_values` (`sv` or `flat`) wins over this key; under `auto` a set flat start is honoured, see [CGMES Import](cgmes_import.md). | Synthetic starts, method studies. | With a projected start (`start_mode.voltage_mode` or `angle_mode` not `classic`) the projection wins and the start is not flat. | Low. | Filed under `start_mode` in the configuration structure. |
| `power_flow.tol` | Float64 | `1.0e-8` | positive real | Convergence bound for the largest single bus mismatch (infinity norm, active and reactive alike; PV rows contribute their voltage residual). Physically `tol * baseMVA`: 1e-8 pu equals 1 W at a 100 MVA base, and the run log and diagnostics print the equivalent. | Accuracy-sensitive studies. | Overly tight on big batches. | Tighter means more iterations. | `max_iter`. |
| `power_flow.tol_MW` | Float64 or unset | unset | positive real | The same bound in physical units. When set it wins over `power_flow.tol` and is converted with the case base at run time (`tol = tol_MW / baseMVA`), because the base belongs to the network, not to the configuration. Unset changes nothing. | Stating the accuracy in MW instead of per unit; the run form offers it next to the per-unit tolerance (leave it empty to use per unit). | A value near the network's own power scale accepts an unsolved state; values at or above 1 MW are warned about. | Zero or negative is refused by name. | `tol`. |
| `power_flow.max_iter` | Int | `80` | positive integer | Iteration cap. | Hard cases. | Very low values. | Upper runtime bound. | `tol`, `qlimits`. |
| `power_flow.autodamp` | Bool | `true` | `true`, `false` | Adaptive damping. | Difficult convergence. | Strict algorithm comparison. | Small overhead, often fewer failures. | `autodamp_min`. |
| `power_flow.autodamp_min` | Float64 | `0.05` | positive real | Minimum damping factor. | Stabilizing hard cases. | Near-zero damping on easy grids. | Lower can increase iterations. | Active only with `autodamp=true`. |
| `power_flow.auto_slack` | Bool | `false` | `true`, `false` | Promote a reference when the case registers no slack: external network injections first, then the largest generator (`ratedS`, then `maxP`, then dispatch). The promotion is logged; without it the run aborts with the no-slack error. | Cases whose data carries no reference (edited or partial models). | Data-quality checks: a silently chosen slack can mask an import problem. | None on cases that already have a slack. | `ensureSlack!` is the underlying API; the CGMES importer applies the same ranking at import time. |
| `power_flow.rescue` | Bool | `true` | `true`, `false` | After a non-converged AC solve, retry from the original start state with a fixed strategy ladder: `alternate_start` (toggle the flat-start flag), `autodamp` (adaptive damping, skipped when already active), `dc_seed` (flat magnitudes with DC-projected start angles), `settled_qlimits` (Armijo merit line search, low damping floor, and Q-limit switching held back until the reactive requests settle via `qlimits.start_mode = :auto`). The first converging strategy wins; its name is logged and recorded as `:ac_rescue_strategy` in the performance profile. `settled_qlimits` targets large systems whose early PV/PQ switching destabilises the iteration. Note the interaction behind that: `qlimits.start_mode = iteration_or_auto` (the shipped default) is an **or**, so a small `qlimits.start_iter` always wins and the `auto` criterion never gets a chance; on large systems, prefer `auto`. | Batch/suite runs over cases with difficult start states. | Solver comparisons and diagnostics: a rescued run hides *which* start failed. | Only failed runs pay for retries; converging runs are untouched. | Config-driven paths only (`runpf!(net, cfg)`, service, Web UI). Distinct from `wrong_branch_rescue`, which reacts to a *converged but implausible* solution. |
| `power_flow.dc.fallback` | Bool | `false` | `true`, `false` | When the AC solve (and the rescue ladder, if enabled) did not converge, run the standalone DC power flow: the net then carries DC angles and branch P flows (vm = 1 pu, no reactive results). The AC status honestly stays non-converged (`erg = 1`); the fallback is logged and recorded as `:dc_fallback_applied` in the performance profile. | Getting *some* flow picture out of a case that AC-diverges. | Any study that needs voltages or reactive power. | One extra linear solve on failed runs. | Uses `power_flow.dc.*` settings; distinct from `power_flow.solver: dc`, which always runs DC only. |
| `power_flow.linear_solver` | Symbol/String | `umfpack_reuse` | `umfpack`, `umfpack_reuse` | Sparse linear-algebra backend for the Newton step of the rectangular solver. `umfpack_reuse` keeps UMFPACK's factorization but reuses the symbolic analysis of the first iteration via `lu!` (analyze once, refactor per iteration). It re-analyzes automatically on active-set pattern changes and falls back to the `umfpack` chain on any factorization error. | `umfpack_reuse` on large cases where the linear solve dominates runtime. | Expecting a third backend: `klu` is not offered and fails validation. | `umfpack_reuse` cuts the repeated symbolic-analysis cost of `umfpack`. | Distinct from `power_flow.solver`, which selects the power-flow *method* (rectangular/apslf/dc); `linear_solver` only affects the rectangular Newton step. With the reuse backend the Jacobian is assembled with a structural (value-independent) sparsity pattern so the analysis stays reusable. Diagnostic counters appear in the solver status (`linear_solver_analyze_count`, `..._refactor_count`, `..._fallback_count`). |

## [Solver selection (rectangular vs. APSLF)](@id pf-solver-selection)

`power_flow.solver` selects the executing solver; `power_flow.method` fixes
the AC formulation family (always `rectangular`). `apslf` routes the run
through the external-solver bridge (`buildPfModel`,
`solvePf(ApslfSolver(...))`, `applyPfSolution!`) to AnalyticLoadFlow.jl's
analytic power-series solver (a required dependency; see `apslf_solver`).

```yaml
power_flow:
  solver: rectangular   # rectangular | apslf
  apslf:
    order: 24
    use_pade: true
    nr_polish: false
    convergence_radius: true
  apslf_start:
    enabled: false
    order: 40
```

| YAML path | Type | Default | Allowed values | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `power_flow.solver` | Symbol/String | `rectangular` | `rectangular`, `apslf` | Selects the executing solver. | `apslf` to use the analytic power-series solver instead of NR. | `apslf` without AnalyticLoadFlow.jl loaded (raises a clear error). | `apslf` skips the NR iteration loop entirely. | Rejects `apslf_start.enabled=true` when set to `apslf`. |
| `power_flow.apslf.order` | Int | `24` | `>= 1` | Highest power-series coefficient computed. | Higher for stressed/large-angle cases. | Unnecessarily high orders on easy cases (cost). | Higher order increases solve cost. | `use_pade`. |
| `power_flow.apslf.use_pade` | Bool | `true` | `true`, `false` | Evaluate the voltage series via Padé `[L/M]` approximants instead of direct Taylor summation. | Default; improves convergence radius. | Direct Taylor comparison studies. | Small evaluation overhead, usually better accuracy per order. | `order`. |
| `power_flow.apslf.nr_polish` | Bool | `false` | `true`, `false` | Run a Newton-Raphson polishing step on the series result. | Debugging the series result against a Newton finish. | Default off: the series alone is a load-flow solution. | Adds a small number of NR iterations. | Only active when `solver=apslf`. |
| `power_flow.apslf.convergence_radius` | Bool | `true` | `true`, `false` | Evaluate the APSLF convergence radius: the distance `dmin` of the nearest Padé pole to the evaluation point `s = 1`, with the bus that owns it and a GRN/YEL/RED level (AnalyticLoadFlow `stability_from_Vcoeff`). Reported in the result header (`APSLF radius`), the run metadata (`apslf_convergence_radius`) and on the runs page next to the Jacobian condition. | Default; judges how far the series solution sits from its continuation limit. | Large networks where the evaluation (about the cost of the solve) is not wanted. | Comparable to the solve itself. | Only active when `solver=apslf`. |
| `power_flow.apslf_start.enabled` | Bool | `false` | `true`, `false` | Use the APSLF solver as a start-value generator ahead of the rectangular NR solve (guarded, like `start_current_iteration`). | Difficult NR starts. | `solver=apslf` (rejected: start generator only makes sense ahead of NR). | Adds one series solve before NR. | `solver`, `apslf_start.order`. |
| `power_flow.apslf_start.order` | Int | `40` | `>= 1` | Series order used by the start-value generator. | Same considerations as `apslf.order`. | Unnecessarily high orders for a start-only pass. | Higher order increases pre-solve cost. | `apslf_start.enabled`. |

`power_flow.apslf_start` has no `use_pade`/`nr_polish` fields (polishing is
left to the NR solve, `nr_polish=false`) and always runs unconstrained
(`Qmin`/`Qmax` are not passed): `power_flow.qlimits.*` governs only the
rectangular NR solve that follows.

### [APSLF start values](@id pf-apslf-start)

`power_flow.apslf_start.enabled` uses the AnalyticLoadFlow.jl-backed APSLF
solver as a guarded start-value generator ahead of the rectangular
Newton-Raphson solve, with the same insertion point and accept/reject guard
style as the current-iteration pre-solve: the candidate is only adopted
when it strictly improves the rectangular mismatch, otherwise the original
start values are restored. Default: disabled. Diagnostic artifact:
`apslf_start.log`.

This mode always runs with no NR polish and no Q-limit enforcement, and
neither is configurable here:

- NR polish is always off internally (`nr_polish=false`) because the
  downstream rectangular Newton-Raphson solve performs that polishing step
  itself.
- Q-limits are always unconstrained during this pre-solve, independent of
  `power_flow.qlimits.enabled` or any other Q-limit setting.
  `power_flow.qlimits.*` only governs the rectangular NR solve that
  follows; the generator's only job is producing a better starting voltage
  profile, not enforcing reactive limits.

Requires AnalyticLoadFlow.jl to be loaded; mutually exclusive with
`power_flow.solver = apslf` (rejected at configuration time: the
start-value generator only makes sense ahead of the NR solve).

`power_flow.apslf_start.order` sets the highest power-series coefficient
used by the start-value generator, with the same considerations as
`power_flow.apslf.order`: higher orders can improve the series
approximation but cost more before the candidate is even evaluated for
acceptance. Default: 40. No effect unless
`power_flow.apslf_start.enabled = true`.

## Solver selection (DC power flow)

`power_flow.solver: dc` selects the standalone DC power flow (MATPOWER
`rundcpf`/`makeBdc`-equivalent): a linear screening model from branch series
reactance only (`B'`, no `r`, no shunt, no line charging), transformer
`phase_shift_deg` as a phase-shift injection vector; `Vm` is `1.0 pu`
everywhere and there are no losses. [`rundcpf!`](@ref) calls it directly,
independent of `run_sparlectra`/`power_flow.solver`.

```yaml
power_flow:
  solver: dc   # rectangular | apslf | dc
  dc:
    angle_reference_deg: 0.0
    ignore_out_of_service: true
```

| YAML path | Type | Default | Allowed values | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `power_flow.solver` | Symbol/String | `rectangular` | `rectangular`, `apslf`, `dc` | Selects the executing solver. | `dc` for a fast linear screening solve or as an AC start-value source. | `dc` when Vm/Q/loss results are required (the model doesn't define them). | `dc` is a single direct linear solve, no iteration. | Rejects active controllers (the outer-loop tap/PST controllers and the Q(U)/P(U) controllers that act inside the Newton step), mirrors `apslf`. |
| `power_flow.dc.angle_reference_deg` | Float64 | `0.0` | any real | Uniform angle offset added to every bus after the slack-referenced solve; the slack bus itself is fixed at this reference. | Matching an external reference-angle convention. | N/A | None (exact post-hoc shift, not a re-solve). | None; mathematically independent of the rest of the DC solve. |
| `power_flow.dc.ignore_out_of_service` | Bool | `true` | `true` | Documents that `status == 0` branches are always excluded from `B'`. | Always (current fixed behavior). | N/A | N/A | Not currently a live toggle. |

**Notes**

- `rundcpf!` accepts a `seed_ac_start::Bool=false` keyword (not a YAML
  option): when `true`, a successful DC solve re-seeds and runs the AC
  rectangular solve from the DC angles (Slack/PV magnitude setpoints
  restored first); `net` then holds the AC solution, the returned
  `DcPowerFlowReport` still reflects the DC step, the AC outcome sits in
  `report.metadata.ac_converged`/`ac_iterations`/`ac_elapsed_s`.
- DC results use a dedicated `DcPowerFlowReport` (angles and lossless
  branch flows) and `dc_pf_status(net)`, a registry separate from
  `rectangular_pf_status`, so a DC result is never mistaken for an AC one.
- With `power_flow.solver = :dc` and multiple AC islands,
  `ac_island_solver_summary.csv` (see below) stays empty for DC-solved
  islands: it reads AC-only `rectangular_pf_status` fields. The DC solve
  and its `dc_pf_status` result are not affected.

## AC island diagnostics

`power_flow.islands` enables structural AC-island diagnostics for networks
with several disconnected AC components. No DC line becomes an AC branch
and no artificial admittance bridge is added; the diagnostics report the
AC topology the rectangular solver sees after import
(`matpower_import.matpower_dcline_mode: pf_injections` remains a fixed
terminal-injection approximation).

```yaml
power_flow:
  islands:
    enabled: true
    mode: solve_independent
    reference_policy: matpower_like
    diagnostic_continue_after_failure: true
```

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `power_flow.islands.enabled` | Bool | `true` | `true`, `false` | Write AC island diagnostics before and after the solve. |
| `power_flow.islands.mode` | Symbol/String | `solve_independent` | `solve_independent`, `solve_parallel` | Island solve mode. `solve_parallel` runs the detected islands concurrently on Julia threads (largest island first), gated by `runtime.parallel.enabled`, `max_tasks`, and `min_work_items`; with one thread, a disabled switch, or too few islands it falls back to the identical serial loop. Results are bitwise identical to `solve_independent`; the fan-out wall clock is recorded as `parallel_wall_time`. One semantic difference: with `diagnostic_continue_after_failure: false` a parallel run cannot skip islands after the first failure (they are already in flight), so every island reports its actual status and the failure is raised after all islands finish; the serial mode keeps today's immediate stop. |
| `power_flow.islands.reference_policy` | Symbol/String | `matpower_like` | `matpower_like` | Select the in-island REF bus when present, otherwise report the first PV/PQ bus that would be promoted for diagnostics. |
| `power_flow.islands.diagnostic_continue_after_failure` | Bool | `true` | `true`, `false` | Keep diagnostics for all detected islands even when the combined run fails. |

**Artifacts** (run output directory):

| Artifact | Content |
|---|---|
| `ac_island_solver_summary.csv` | One row per detected AC island: selected reference bus, PV/PQ/REF counts, propagated solver settings, final status, final mismatch, failure reason. Islands the solver never attempted (de-energized single-bus islands, islands skipped after an earlier failure) appear with `final_status=not_attempted`, `failure_reason=not_attempted`, `stage=not_attempted`, `iterations=0`, zeroed switching statistics and `unavailable` mismatch fields. |
| `ac_island_<id>_solver.log`, `ac_island_<id>_mismatch_history.csv` | Per-island details with the same topology, setting and solve-status fields, written only for islands the solver actually attempted. |
| `q_limit_processing_status` column | Q-limit-specific outcome only: `disabled` when Q-limits are off, `not_attempted` for islands never solved, otherwise the solver-recorded outcome (currently `unavailable`, the solver stores no dedicated Q-limit processing status). It never mirrors `failure_reason`; Q-limit activity is in the switching-statistics columns (`pv_pq_switching_events`, `qlimit_active_set_changes`, ...). |

**Notes**

- When an island-aware run fails, the API/result message names the first
  failing island: id, bus and branch counts, reference bus, PV/PQ/REF
  counts, iteration count, final mismatch, mismatch status (`finite`,
  `nonfinite`, `NaN`, or `Inf`), failure reason, stage, start-projection
  setting, and the artifact path, taken from that island's own solver
  record; the generic `before_nr`/`during_nr` stage heuristic is used only
  when no per-island record exists.
- Island-wise solving is structural support, not a convergence guarantee:
  if any island or the combined solve fails, the run is reported as failed
  and the island artifacts are kept.

For SyntheticUSA island diagnostics, start with Q-limit handling disabled
so the baseline tests topology, reference selection and start projection
before active-set effects:

```yaml
matpower_import:
  matpower_dcline_mode: pf_injections
  auto_profile: apply
  compare_voltage_reference: hybrid

power_flow:
  tol: 1.0e-5
  max_iter: 80
  autodamp: true
  autodamp_min: 0.01

  islands:
    enabled: true
    mode: solve_independent
    reference_policy: matpower_like
    diagnostic_continue_after_failure: true

  qlimits:
    enabled: false

  start_current_iteration:
    enabled: false

  start_mode:
    angle_mode: dc
    voltage_mode: profile_blend
    profile_source: matpower_reference
    start_projection: true
    try_dc_start: true
    try_blend_scan: true
    branch_guard: true
    measure_candidates: true
    reuse_import_data: true
```

Then enable `power_flow.qlimits.enabled`: per-island logs add
`q_limit_processing_status`, switching events, active-set changes, reenable
events, guarded narrow-range PV buses, the final PV voltage residual and
the mismatch metrics.

## [Distributed active-power slack](@id pf-distributed-slack)

`power_flow.distributed_slack` spreads the active-power imbalance of an
island (load + losses − scheduled generation) over the participating
generators instead of loading it onto the reference bus. Theory, augmented
Newton system and applicability rules: [Solver Guide](solver.md).

```yaml
power_flow:
  distributed_slack:
    enabled: false
    p_mode: pg_weighted
    respect_p_limits: true
    fallback: error
    weights: {}
```

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `power_flow.distributed_slack.enabled` | Bool | `false` | `true`, `false` | Master switch. Disabled runs are bit-identical to the classical single-slack solver; imported participation factors alone never activate the feature. |
| `power_flow.distributed_slack.p_mode` | Symbol/String | `pg_weighted` | `pg_weighted`, `pmax_weighted`, `headroom_weighted`, `imported`, `explicit` | Weight source for the participation factors: scheduled `Pg`, `maxP`, remaining headroom `max(maxP − Pg, 0)`, the imported factor (MATPOWER gen column 21 `APF`, CGMES `GeneratingUnit.normalPF`), or the explicit `weights` table. |
| `power_flow.distributed_slack.respect_p_limits` | Bool | `true` | `true`, `false` | Warn per participant whose corrected output `Pg + alpha·lambda_P` leaves `[minP, maxP]` by more than a relative tolerance (0.01 % of the limit, floor 1e-3 MW; epsilon overshoots of participants scheduled exactly at a limit are numerical noise). Warns only; it does not clamp and re-solve. |
| `power_flow.distributed_slack.fallback` | Symbol/String | `error` | `error`, `ref_only` | Behavior when an island has no valid participant for the chosen mode: abort with an error, or warn and solve that island classically (reference bus absorbs everything). |
| `power_flow.distributed_slack.weights` | Mapping | `{}` | bus name/index → weight ≥ 0 | Only read with `p_mode: explicit`; must then be non-empty with at least one positive weight. Keys resolve against bus names first, then against bus indices written as strings. |

**Notes**

- Candidates are the generator-type prosumers at the island's REF and PV
  buses; fixed injections at PQ buses (Stage-0 HVDC converter injections,
  kept boundary equivalent injections) never participate. Invalid
  candidates (missing data, non-finite or non-positive weight) are dropped
  with a debug log and counted in the result metadata; surviving weights
  are normalized to `sum(alpha) = 1` per island.
- In island-wise runs each island solves with its own `lambda_P`, reported
  in the per-island solver statuses.

| Surface | Fields |
|---|---|
| Structured status (next to the wrong-branch metadata) | `distributed_slack_active`, `distributed_slack_mode`, `distributed_slack_lambda_p_pu`, `distributed_slack_lambda_p_mw`, `distributed_slack_participants`, `distributed_slack_alpha_sum`, `distributed_slack_dropped`, `distributed_slack_p_limit_violations`, and the per-participant table `distributed_slack_participation` with bus, alpha share, correction `dP` and scheduled output; at `verbose > 0` a compact summary with the top participants is printed. |
| `printACPFlowResults` | Columns `dSl alpha` and `Pg eff MW` on participating buses (the `Pg` column keeps the schedule) and a one-line summary with mode and `lambda_P` in the result header. |

## [External grid source](@id pf-external-grid)

`power_flow.external_grid` computes the marked slack bus as a non-ideal
external-grid source: before the solve, `convertSlackToExternalGrid!` moves
the reference voltage to a hidden internal bus `<bus>__extgrid_int` behind
the feeder impedance `z = Un²/Sk''` (split by the R/X ratio), and the former
slack bus becomes an ordinary solved bus whose voltage droops under load.
Theory (both formulations, the stiff limit, the short-circuit effect):
[Slack Bus and External Grid Sources](slack_vs_source.md).

```yaml
power_flow:
  external_grid:
    enabled: false
    source: auto
    sk_MVA: 2000.0
    rx: 0.1
```

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `power_flow.external_grid.enabled` | Bool | `false` | `true`, `false` | Master switch. Off keeps the classical ideal slack. |
| `power_flow.external_grid.source` | Symbol/String | `auto` | `auto`, `config` | Where `Sk''`/`R/X` come from. `auto` prefers the values the case declares on the slack bus (a CGMES `ExternalNetworkInjection`, or a Sparlectra Case Format `source` with `sk`/`rx_ratio`; logged as *declared by the case data*) and falls back to the config numbers (MATPOWER/DTF cases carry no such data). `config` always uses the config numbers. |
| `power_flow.external_grid.sk_MVA` | Float | `2000.0` | > 0 | Initial symmetrical short-circuit power of the feeder. The series impedance is `z_pu = baseMVA/sk_MVA` on the per-voltage-level base. |
| `power_flow.external_grid.rx` | Float | `0.1` | ≥ 0 | R/X ratio of the feeder impedance. |

**Notes**

- Only the primary slack bus is converted; with `multi_slack` every other
  island keeps its reference. The conversion is logged in `run.log`; rescue
  retries cannot stack sources.
- The classical result print names the connection in its `Grid connection:`
  header line and reports the internal reference bus with type `SOURCE`.
  Web UI: the **External grid source** fieldset of the advanced run options.
- `power_flow.external_grid.enabled` is mutually exclusive with
  `power_flow.distributed_slack.enabled` (configuration error).

!!! details "Why source and distributed slack exclude each other"
    Both decide who covers the island's power imbalance, the source by
    importing it through the feeder, the distributed slack by spreading it
    over the island's generators. Combined, the source would be forced to
    a participation share of zero and degenerate to a bare angle anchor.

## [Start mode options](@id pf-start-mode)

Where the Newton-Raphson iteration starts: the angle and voltage initialisation modes, the profile a start can be taken from, the start projection and the DC and blend attempts that run before the first step. A flat start switches all of them off.

| YAML path | Type | Default | Allowed values | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `power_flow.start_mode.angle_mode` | Symbol/String | `dc` | `classic`, `dc`, `bus_va_blend`, `matpower_va` | Angle initialization source mode. | Transmission-like starts. | If trusted measured/historical state exists. | Can reduce iterations. | `try_dc_start`, `dc_angle_limit_deg`. |
| `power_flow.start_mode.voltage_mode` | Symbol/String | `profile_blend` | `classic`, `pv_gen_vg`, `pv_bus_vm`, `all_bus_vm`, `profile_blend` | Voltage-magnitude/angle blend strategy. | Imported-reference assisted starts. | Untrusted imported data. | Small startup overhead. | `blend_lambdas`, `reuse_import_data`. |
| `power_flow.start_mode.profile_source` | Symbol/String | `matpower_reference` | `flat`, `dc`, `bus_metadata`, `historical_profile`, `matpower_reference`, `state_estimation`, `se_snapshot`, `scada_snapshot` | Source for external or model-derived start profiles. `matpower_reference` means imported `BUS.VM`/`BUS.VA` values for regression/benchmarking/import reproduction. | Profile-aware starts. | When source data is unavailable or untrusted. | Minimal parsing overhead. | Keep source explicit for diagnostics and reproducibility. |
| `power_flow.start_mode.profile_source = state_estimation` | Symbol/String | see above | see above | SE chain start: the power flow starts from the estimated voltages of a preceding state estimation; the model injections stay authoritative, the measurement/model difference goes into the slack. Programmatic entry `runpf_from_se!(...; mode = :se_state)`. | Model-authoritative PF after an SE. | No preceding SE result (clear error). | Starts near the solution. | `runse!(updateNet = true)`, `readSEStateCSV!`, the Web UI chain action. |
| `power_flow.start_mode.profile_source = se_snapshot` | Symbol/String | see above | see above | SE chain start with balance takeover: additionally the nodal balances come from the estimation (working net only, the persistent model is never mutated). Converges in 0 or 1 iterations with the slack pickup below tolerance. Programmatic entry `runpf_from_se!(...; mode = :se_snapshot)`. | Snapshot-authoritative PF (reproduce the estimated operating point). | The model injections must stay authoritative (use `state_estimation`). | Immediate convergence. | Result metadata records the slack pickup. |
| `power_flow.start_mode.start_projection` | Bool | `true` | `true`, `false` | Enables start-projection workflow. | Robustness on hard cases. | Minimal-path microbench runs. | Extra startup pass. | Gates start-projection sub-options. |
| `power_flow.start_mode.try_dc_start` | Bool | `true` | `true`, `false` | Try DC candidate start. | Large transmission cases. | Highly resistive distribution cases. | Low overhead. | `dc_angle_limit_deg`. |
| `power_flow.start_mode.try_blend_scan` | Bool | `true` | `true`, `false` | Scan blend candidates. | Mixed/tricky start data. | Easy cases needing speed. | Startup cost ∝ lambda count. | `blend_lambdas`. |
| `power_flow.start_mode.branch_guard` | Bool | `true` | `true`, `false` | Branch sanity guard for candidate starts. | Stability-focused runs. | Rarely disabled. | Low. | Candidate measurement options. |
| `power_flow.start_mode.measure_candidates` | Bool | `true` | `true`, `false` | Score/select among candidates. | Multiple start candidates. | Fastest startup path only. | Low/medium startup overhead. | `try_dc_start`, `try_blend_scan`. |
| `power_flow.start_mode.accept_unmeasured_dc_start` | Bool | `false` | `true`, `false` | Allow DC start without measurement checks. | Synthetic studies. | Measurement-driven workflows. | Can avoid fallback retries. | `try_dc_start`. |
| `power_flow.start_mode.dc_seed_unconditional` | Bool | `false` | `true`, `false` | Unconditionally run a full standalone DC power flow first and seed Newton-Raphson's start angles from it, bypassing `angle_mode`/`try_dc_start`/`measure_candidates` entirely (no quality gate). | You explicitly want the `rundcpf!(seed_ac_start=true)` two-step behavior from the config-driven path, with normal diagnostics/artifacts/Q-limits still applying to the AC solve. | You want the existing measured/guarded candidate selection (the default) to keep the option to fall back away from a bad DC candidate. | One extra DC solve (cheap, linear) before every rectangular NR run. | `power_flow.solver` must be `rectangular` (rejected for `apslf`/`dc`); mutually exclusive with `power_flow.apslf_start.enabled`. |
| `power_flow.start_mode.reuse_import_data` | Bool | `true` | `true`, `false` | Reuse imported MATPOWER references. | Trusted imports. | Uncertain conversion data. | Small reduction in recomputation. | `matpower_import.*` voltage reference keys. |
| `power_flow.start_mode.blend_lambdas` | Vector{Float64} | `[0.25,0.5,0.75]` | real vector (typ. 0..1) | Lambda candidates for blend scan. | Need robust candidate search. | Very large lambda sets. | Linear startup growth with vector size. | `try_blend_scan`. |
| `power_flow.start_mode.dc_angle_limit_deg` | Float64 | `60.0` | positive real | DC-start angle magnitude cap (deg). | Conservative angle starts. | Overly restrictive values. | Negligible. | `try_dc_start`. |

### DC-seeded Newton-Raphson start (`dc_seed_unconditional`)

`angle_mode = dc` (the default) is one candidate inside the
measured/guarded start projection (`project_rectangular_start`), a
lightweight DC-angle estimate with a fallback if it looks bad.
`power_flow.start_mode.dc_seed_unconditional = true` instead runs the
standalone DC power flow (the solver of `power_flow.solver = dc`,
per-island handling included) before Newton-Raphson and always uses its bus
angles, like [`rundcpf!(net; seed_ac_start=true)`](@ref) but inside the
config-driven `run_sparlectra`/Web UI pipeline, so Q-limit handling,
wrong-branch detection and diagnostics/artifacts still apply to the AC
solve.

Voltage magnitudes are unaffected (`Vm = 1.0 pu` in the DC model):
Slack/PV setpoints are preserved, other buses get their magnitude from the
normal `voltage_mode`/`start_projection` handling. Only
with `power_flow.solver = rectangular`; rejected for `apslf`/`dc` and
together with `power_flow.apslf_start.enabled` (mutually exclusive
start-value sources). Web UI: **Use DC start values**; while active, the
**Start angle mode** and **Start voltage mode** selects are grayed out.

## Guarded current-iteration start pre-solve

`power_flow.start_current_iteration` enables an optional guarded
current-injection/current-iteration pre-solve for start values: a
start-value preconditioner, not a new solver; the final AC solve remains
the rectangular Newton-Raphson path.

```text
Start Voltage Mode + Start Angle Mode
→ optional start projection / candidate selection
→ optional current-iteration pre-solve
→ Newton-Raphson power flow
→ optional Q-limit handling / outer loop logic
```

It adds no new `power_flow.start_mode.voltage_mode` or
`power_flow.start_mode.angle_mode` value: it consumes the voltage profile
prepared by start mode and start projection, attempts a limited PQ-bus
current update, and keeps the candidate only when the guards accept it;
otherwise the original start values are restored.

### [Current-iteration start options](@id pf-current-iteration-start)

`power_flow.start_current_iteration.enabled` enables a guarded
current-injection/current-iteration pre-solve before the Newton-Raphson
power-flow solver starts. This is not a separate power-flow solver and it
does not replace Newton-Raphson. It is a start-value preconditioner:
Sparlectra first builds the initial voltage profile from Start Voltage Mode
and Start Angle Mode, then optionally tries a few current-iteration steps to
improve that initial profile. The improved voltage profile is accepted only
if it passes the voltage and angle guards and improves the existing
Sparlectra mismatch metric. If it does not improve the start, the original
start values are restored and Newton-Raphson starts normally. Default:
disabled. Enable this only for difficult cases where the normal start
profile or DC/profile-blend start is not robust enough. Diagnostic
artifact: `current_iteration_start.log`.

```yaml
power_flow:
  start_current_iteration:
    enabled: false
    max_iter: 10
    tol: 1.0e-3
    damping: 0.5
    accept_only_if_improved: true
    min_improvement_factor: 0.98
    vm_min_pu: 0.5
    vm_max_pu: 1.5
    max_angle_step_deg: 30.0
    only_for_large_cases: false
```

| YAML path | Type | Default | Meaning | When to use |
|---|---:|---:|---|---|
| `power_flow.start_current_iteration.enabled` | Bool | `false` | Enable the guarded current-iteration pre-solve. | The pre-solve above; keep it off unless the normal start profile or a DC/profile-blend start is not robust enough. |
| `power_flow.start_current_iteration.max_iter` | Int | `10` | Maximum number of current-iteration pre-solve steps before Newton-Raphson starts. | A higher value gives the pre-solve more chances to reduce the initial mismatch, but costs time and may move the start profile too far from the original initialization; the result is still guarded. Keep this small, increase it only when diagnostics show that the mismatch keeps improving but the pre-solve stops too early. |
| `power_flow.start_current_iteration.tol` | Float64 | `1.0e-3` | Stopping tolerance of the pre-solve: below it the pre-solve stops before `max_iter`. | Controls only the start-value pre-solve, not the final Newton-Raphson tolerance. Use a relatively loose value; the purpose is a better starting point, not a solved power flow. |
| `power_flow.start_current_iteration.damping` | Float64 | `0.5` | Damping factor in `(0, 1]` for the current-iteration voltage update; 1.0 applies the full update. | Smaller values blend the update with the previous voltage and make the pre-solve more conservative, which avoids large voltage or angle jumps. Lower it if the pre-solve is rejected by the voltage or angle guards, raise it only if the pre-solve is stable but improves too slowly. |
| `power_flow.start_current_iteration.accept_only_if_improved` | Bool | `true` | Accept the candidate only when it improves the existing Sparlectra mismatch metric; otherwise the original start values are restored. | Should normally stay enabled. Disabling it is only useful for expert experiments, because it can let a worse start profile enter Newton-Raphson. |
| `power_flow.start_current_iteration.min_improvement_factor` | Float64 | `0.98` | Required improvement ratio when `accept_only_if_improved` is on; the best candidate mismatch must be at most this factor times the original mismatch. | 0.98 means the candidate mismatch must be at least about 2 percent lower than the original. Smaller values require a stronger improvement, values closer to 1.0 accept smaller improvements. Keep it close to 1.0 for a conservative pre-solve, lower it only when tiny improvements are not useful and you want to accept only clearly better starts. |
| `power_flow.start_current_iteration.vm_min_pu` | Float64 | `0.5` | Lower voltage-magnitude guard: a candidate with any bus voltage below it is rejected and the original start values are restored. | Prevents the pre-solve from sending Newton-Raphson into an implausible low-voltage start region. Lowering makes the guard more permissive, raising makes the pre-solve more conservative; `current_iteration_start.log` shows the candidate voltage minima. |
| `power_flow.start_current_iteration.vm_max_pu` | Float64 | `1.5` | Upper voltage-magnitude guard: a candidate with any bus voltage above it is rejected and the original start values are restored. | Prevents unrealistic over-voltage start profiles from entering Newton-Raphson. Lowering makes the guard stricter, raising allows larger candidate voltages; `current_iteration_start.log` shows the candidate voltage maxima. |
| `power_flow.start_current_iteration.max_angle_step_deg` | Float64 | `30.0` | Maximum allowed angle change of a single current-iteration update; a larger jump rejects the candidate and restores the original start values. | Guards against unstable or wrong-branch start profiles. Lower it for a more conservative pre-solve, raise it only when diagnostics show that otherwise plausible candidates are rejected solely by this guard. |
| `power_flow.start_current_iteration.only_for_large_cases` | Bool | `false` | Run the pre-solve only for cases Sparlectra classifies as large enough for this extra start-value preparation (the large-case threshold of the rectangular workspace). | Avoids spending time on small cases where normal start values usually work. Enable it to have the pre-solve available for difficult large MATPOWER cases without changing the behaviour for small examples. |

### Recommended usage

For difficult MATPOWER starts, a conservative setup is:

```yaml
power_flow:
  start_mode:
    voltage_mode: profile_blend
    angle_mode: dc
    profile_source: matpower_reference
  start_current_iteration:
    enabled: true
    max_iter: 10
    damping: 0.5
    accept_only_if_improved: true
```

This may reduce the initial mismatch, but it is experimental and does not
guarantee a rescue (a rejected candidate leaves Newton-Raphson on the
original start values) and does not replace MATPOWER auto-profile
selection, DC angle starts, or start projection.

### Diagnostics and interpretation

When a run has an output directory in its performance profile, the
pre-solve writes `current_iteration_start.log`:

| Group | Fields |
|---|---|
| Outcome | `current_iteration_enabled`, `current_iteration_attempted`, `current_iteration_accepted`, `current_iteration_reason` |
| Mismatch | `initial_mismatch`, `final_mismatch`, `iterations` |
| Candidate voltages | `candidate_voltage_magnitude_min`, `candidate_voltage_magnitude_max`, `candidate_voltage_low_count`, `candidate_voltage_high_count`, `candidate_voltage_worst_low_bus`, `candidate_voltage_worst_high_bus`, and their corresponding values |
| Candidate angles | `candidate_max_angle_step_deg`, `maximum_angle_step_deg` |
| Rejection | `guard_violations`, `rejection_stage`, `rejected_at_iteration`, `original_start_values_restored` |
| Restored ranges after rejection | `restored_voltage_magnitude_min`, `restored_voltage_magnitude_max` |

```text
current_iteration_accepted: true
  The pre-solve result was used as Newton-Raphson start values.

current_iteration_accepted: false and original_start_values_restored: true
  The candidate was rejected and Newton-Raphson started from the original start values.

current_iteration_reason: voltage_magnitude_guard
  At least one candidate voltage was outside vm_min_pu/vm_max_pu.

current_iteration_reason: angle_step_guard
  The candidate changed an angle beyond max_angle_step_deg.

current_iteration_reason: not_improved
  The candidate did not improve the mismatch enough.
```

Other reasons: `disabled`, `skipped_small_case`, `max_iter`,
`tolerance_reached`, `invalid_voltage`, `singular_current_update`, and
`invalid_mismatch`.

**Notes**

- Active-set path: the pre-solve runs before the Newton-Raphson solve;
  classical MATPOWER-style Q-limit outer-loop modes apply it only to the
  first inner solve. It does not change Q-limit switching decisions.
- It can be rejected by the voltage, angle-step, finite-value,
  singular-update or mismatch-improvement guards, and does not help cases
  whose main issue is a model/convention mismatch, a wrong branch, or an
  invalid MATPOWER import convention.

## [Merit-function line search options](@id pf-merit)

`power_flow.merit` enables an optional Armijo sufficient-decrease acceptance
test inside the autodamp backtracking loop of the rectangular
Newton-Raphson solver, not a replacement for Newton-Raphson or autodamp.
Disabled by default (`enabled: false`), which leaves the max-mismatch
autodamp behavior unchanged. Background:
[Merit-Function Line Search](@ref).

```yaml
power_flow:
  autodamp: true
  merit:
    enabled: false
    armijo_c1: 1.0e-4
    scale_p: 1.0
    scale_q: 1.0
    scale_v: 1.0
    fallback_max_mismatch: true
```

| YAML path | Type | Default | Allowed | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `power_flow.merit.enabled` | Bool | `false` | `true`, `false` | Master switch for the Armijo merit-function line search. | Diagnosing or improving step acceptance on difficult flat-start cases where the ∞-norm autodamp criterion accepts poor steps. | Strict comparison against historical autodamp-only behavior. | One extra weighted residual-norm evaluation per already-computed trial mismatch; negligible. | Requires `power_flow.autodamp = true`; validation error otherwise. |
| `power_flow.merit.armijo_c1` | Float64 | `1.0e-4` | real in `(0, 0.5)` | Sufficient-decrease constant `c₁` in the Armijo condition. | Tuning how strict the accepted decrease must be. | N/A | Larger values reject more trials, increasing backtracking steps. | Only active when `merit.enabled = true`. |
| `power_flow.merit.scale_p` | Float64 | `1.0` | positive real | Diagonal weight for active-power (`ΔP`) residual entries in the merit function. | P/Q/V residuals differ in magnitude and one wants a balanced merit value. | Default per-unit systems where residuals are already comparable. | None (evaluated alongside the existing mismatch check). | YAML-only; not exposed in the Web UI. |
| `power_flow.merit.scale_q` | Float64 | `1.0` | positive real | Diagonal weight for reactive-power (`ΔQ`) residual entries (PQ buses). | Same as `scale_p`. | Same as `scale_p`. | None. | YAML-only; not exposed in the Web UI. |
| `power_flow.merit.scale_v` | Float64 | `1.0` | positive real | Diagonal weight for voltage-setpoint (`ΔV`) residual entries (PV buses). | Same as `scale_p`. | Same as `scale_p`. | None. | YAML-only; not exposed in the Web UI. |
| `power_flow.merit.fallback_max_mismatch` | Bool | `true` | `true`, `false` | Behavior when no backtracking trial satisfies the Armijo condition. `true` falls back to the existing max-mismatch criterion (first improving trial, else the conservative best-finite trial); `false` skips straight to the conservative best-finite-trial fallback. | Most cases (`true`, safest). Use `false` only to force the most conservative step whenever Armijo fails. | N/A | `false` can pick smaller steps than `true` would, increasing iteration count. | Only active when `merit.enabled = true`. |

### Recommended usage

```yaml
power_flow:
  autodamp: true
  merit:
    enabled: true
    armijo_c1: 1.0e-4
    fallback_max_mismatch: true
```

Enable `power_flow.merit.enabled` with `power_flow.autodamp` where the
∞-norm autodamp criterion is suspected to accept a step that reduces the
worst-bus mismatch while increasing the overall residual energy. Leave
`scale_p`/`scale_q`/`scale_v` at `1.0` unless P/Q/V residuals differ by
orders of magnitude.

### Diagnostics and interpretation

| Surface | Fields |
|---|---|
| `merit_linesearch.log` (one line per Newton iteration, written when the run has an output directory in its performance profile and `merit.enabled = true`) | `f_before`, `directional_derivative`, `tested_alphas`, `accepted_alpha`, `accept_reason` |
| Solver status | `merit_enabled`, `merit_used_iterations`, `merit_fallback_count`, `merit_active_set_skip_count`, `merit_initial`, `merit_final` |

`accept_reason` values:

```text
armijo
  A backtracking trial satisfied the Armijo sufficient-decrease condition; it was accepted directly.

fallback_max_mismatch
  No trial satisfied Armijo; fallback_max_mismatch=true, and the classic max-mismatch criterion
  found an improving trial, which was accepted instead.

fallback_conservative
  No trial satisfied Armijo (and, if fallback_max_mismatch=true, no trial improved the max
  mismatch either); the most conservative finite trial from the backtracking sweep was accepted.

active_set_skip
  A PV/PQ active-set switch (Q-limit handling) happened during this Newton iteration, so the
  residual vector's entries changed meaning; the merit comparison was skipped for this iteration
  and the classic max-mismatch criterion was used instead.
```

**Notes**

- The merit function `f(x) = 1/2 ‖W F(x)‖²` can have local minima with
  `F(x) ≠ 0`; satisfying Armijo does not by itself guarantee convergence.
- It neither selects the solution branch (high vs. low voltage) nor changes
  the candidate start-value ranking (`start_mode.measure_candidates`),
  which stays mismatch-based.
- PV/PQ active-set switches change the meaning of `f`; the merit comparison
  is skipped for that iteration (`accept_reason = active_set_skip`).

## [Trust-region step control options](@id pf-trust-region)

`power_flow.trust_region` enables an optional scaled-Newton trust-region
alternative to `power_flow.autodamp`: it caps the Newton step norm at an
adaptive radius and accepts or rejects trials by merit-function decrease
rather than by the max-mismatch criterion. Disabled by default
(`enabled: false`); mutually exclusive with `power_flow.autodamp = true`
(both control the step length, enabling both is a configuration error).
Background: [Trust-Region Step Control](@ref).

```yaml
power_flow:
  autodamp: false
  trust_region:
    enabled: false
    initial_radius: 1.0
    min_radius: 1.0e-4
    max_radius: 10.0
    eta_accept: 0.1
    shrink_factor: 0.5
    expand_factor: 2.0
    expand_threshold: 0.75
    step_mode: scaled
```

| YAML path | Type | Default | Allowed | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `power_flow.trust_region.enabled` | Bool | `false` | `true`, `false` | Master switch for scaled-Newton trust-region step control. | Difficult flat-start cases where a merit-decrease step-acceptance rule (rather than max-mismatch backtracking) is preferred. | Together with `power_flow.autodamp = true` (validation error). | One extra weighted-residual evaluation and a sparse matrix-vector product per trial, reusing the already-built Jacobian; no extra factorization. | Requires `power_flow.autodamp = false`. |
| `power_flow.trust_region.step_mode` | Symbol/String | `scaled` | `scaled`, `dogleg` | Trial-step construction: `scaled` rescales the full Newton direction to the radius (default, unchanged behavior); `dogleg` blends the Newton direction with a steepest-descent (Cauchy) step along the dogleg path when the radius shrinks below the Newton step norm. See [Trust-Region Step Control](@ref), "Dogleg step mode". | Cases where the Newton direction is suspected to degrade (become a poor descent direction) partway through the solve, and graceful degradation is preferred over repeated rescaling toward collapse. | Cases already converging cleanly under `scaled`; `dogleg` trades some speed (possible Cauchy-direction crawling) for robustness and has no benefit when `scaled` already works. | Two extra sparse matrix-vector products (`Jᵀ(WF)`, `Jg`) and a closed-form scalar root per accepted/rejected trial; no extra factorization. | Only meaningful when `enabled = true`. Does not change the `autodamp`/`merit` mutual-exclusion rules. |
| `power_flow.trust_region.initial_radius` | Float64 | `1.0` | positive, `> min_radius`, `<= max_radius` | Starting trust-region radius, in per-unit state-vector (2-)norm. | Tuning how large the first trial step may be. | N/A | Larger values risk more rejected/shrunk first steps on hard cases. | Bounded by `min_radius`/`max_radius`. |
| `power_flow.trust_region.min_radius` | Float64 | `1.0e-4` | positive, `< initial_radius` | Radius floor. Falling below it declares non-convergence (`reason = :trust_region_collapsed`) instead of looping indefinitely. | Bounding worst-case retry effort. | Setting it near `initial_radius`, which makes collapse detection overly aggressive. | Smaller values allow more shrink retries before giving up. | Read every retry inside one Newton iteration. |
| `power_flow.trust_region.max_radius` | Float64 | `10.0` | positive, `>= initial_radius` | Radius ceiling; the radius never expands past this value. | Preventing runaway expansion on well-behaved cases. | N/A | Higher ceiling allows larger accepted steps once the model is trusted. | Caps the `expand_factor` growth. |
| `power_flow.trust_region.eta_accept` | Float64 | `0.1` | positive real | Minimum actual/predicted reduction ratio `rho` required to accept a trial step. | Tuning acceptance strictness. | N/A | Higher values reject more trials, increasing shrink/retry iterations. | Compared against `rho` computed from the merit function. |
| `power_flow.trust_region.shrink_factor` | Float64 | `0.5` | `(0, 1)` | Radius multiplier applied on a rejected trial. | Tuning how aggressively the radius shrinks after a bad step. | N/A | Smaller values shrink faster, reaching `min_radius`/collapse sooner. | Applied repeatedly within one Newton iteration until accepted or collapsed. |
| `power_flow.trust_region.expand_factor` | Float64 | `2.0` | `> 1` | Radius multiplier applied on a strongly successful boundary-hitting step. | Letting the radius grow once the model proves reliable. | N/A | Larger values reach `max_radius` faster. | Only applied when the step also hit the radius boundary. |
| `power_flow.trust_region.expand_threshold` | Float64 | `0.75` | `(0, 1)` | `rho` threshold above which an accepted, boundary-hitting step triggers expansion. | Tuning how "good" a step must be before trusting a larger radius. | N/A | Higher values expand less often, growing the radius more conservatively. | Must be `>= eta_accept` in practice for a coherent expand/accept ordering (not separately validated). |

### Diagnostics and interpretation

| Surface | Fields |
|---|---|
| `trust_region.log` (one line per Newton iteration, written when the run has an output directory in its performance profile and `trust_region.enabled = true`) | `radius_before`, `rho`, `tested_radii` (the shrink sequence tried this iteration), `rejected_steps`, `accepted`, `radius_after`, `collapsed`, `accept_reason` |
| Solver status | `trust_region_enabled`, `tr_step_count`, `tr_rejected_steps`, `tr_min_radius`, `tr_max_radius` (the observed radius range across the run, not the configured bounds), `tr_final_radius`, `tr_collapsed`, and the dogleg-branch counters `tr_dogleg_newton_count`, `tr_dogleg_interp_count`, `tr_dogleg_cauchy_count`, `tr_active_set_skip_count` (all `0` in `scaled` mode) |

`accept_reason` values: `scaled` (the only reason of the `scaled` step
mode); `dogleg_newton`/`dogleg_interp`/`dogleg_cauchy` (the accepted point
on the dogleg path, `dogleg` mode only); `active_set_skip` (a PV/PQ switch
this iteration, `dogleg` mode only: the comparison was skipped and one
`scaled`-style trial accepted unconditionally); `none` (collapse, no step
accepted).

```text
tr_collapsed: false
  The trust region never dropped below min_radius; every Newton iteration eventually accepted a step.

tr_collapsed: true
  The radius fell below min_radius without an accepted step in some iteration; the run reports
  reason = :trust_region_collapsed. Lower min_radius, raise initial_radius, or fall back to
  autodamp for this case.

tr_rejected_steps > 0
  At least one trial was rejected (rho < eta_accept) and the radius was shrunk before an
  eventual accept (or collapse). A high count relative to tr_step_count suggests the model
  is a poor local predictor for this case; consider a different start profile.
```

## Step-control mode combinations (autodamp / merit / trust-region)

`power_flow.autodamp`, `power_flow.merit.enabled`, and
`power_flow.trust_region.enabled` jointly select the Newton step length.
`PowerFlowConfig` enforces the rules at load time (`ArgumentError`) for
YAML, API `config_overrides` and the Web UI alike:
`power_flow.merit.enabled = true` requires `autodamp = true` (merit is a test
inside the autodamp loop), `power_flow.trust_region.enabled = true` requires
`autodamp = false` (trust-region replaces that loop), so merit and
trust-region can never both hold.

| `autodamp` | `merit.enabled` | `trust_region.enabled` | Valid? | Behavior |
|---|---|---|---|---|
| `false` | `false` | `false` | Valid | Fixed-damping Newton step (no autodamp backtracking, no merit, no trust-region). |
| `false` | `false` | `true`  | Valid | Trust-region step control alone. |
| `false` | `true`  | `false` | **Invalid** | `ArgumentError`: merit requires `autodamp = true`. |
| `false` | `true`  | `true`  | **Invalid** | `ArgumentError`: merit requires `autodamp = true` (and trust-region already forbids `autodamp = true`). |
| `true`  | `false` | `false` | Valid | Classic autodamp backtracking (max-mismatch acceptance), unchanged default behavior. |
| `true`  | `false` | `true`  | **Invalid** | `ArgumentError`: trust-region requires `autodamp = false`. |
| `true`  | `true`  | `false` | Valid | Autodamp backtracking with the Armijo merit-function acceptance test. |
| `true`  | `true`  | `true`  | **Invalid** | `ArgumentError`: both merit and trust-region reject this combination of `autodamp`. |

The Web UI mirrors these rules: the *Autodamping & merit-function line
search* box and the *Trust-region step control* box exclude each other, and
the merit toggle is disabled whenever autodamping is off.
`power_flow.trust_region.step_mode` (`scaled`/`dogleg`) is orthogonal to
this matrix and only takes effect where `trust_region.enabled = true`.

## [Q-limit options and guard](@id pf-qlimits)

How the reactive-power limits of PV machines are enforced: the master switch, the enforcement mode (active set or the classic outer loop), when the enforcement starts, the hysteresis around a limit, the final check that every mode ends with, and the guard that judges narrow ranges and remaining violations.

| YAML path | Type | Default | Allowed values | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `power_flow.qlimits.enabled` | Bool | `true` | `true`, `false` | Master switch for Q-limit enforcement. | Realistic PV/PQ behavior needed. | Pure unconstrained PF tests. | Can increase switching iterations. | Gates all `qlimits.*`. |
| `power_flow.qlimits.enforcement_mode` | Symbol/String | `active_set` | `active_set`, `classic_simultaneous`, `classic_one_at_a_time`; `off` disables the handling like `enabled: false` | Selects active-set or classical Q-limit switching. | Compare dynamic switching with classical outer-loop enforcement. | Treating legacy `matpower_*` aliases as preferred public values. | Classical modes may require repeated solves. | See [Q-limit Switching Strategy](q_limit_switching_strategy.md). |
| `power_flow.qlimits.start_iter` | Int | `3` | integer | First iteration index for Q limits. | Delay switching noise early. | Very late switching on hard cases. | Affects convergence speed. | `start_mode`, `auto_q_delta_pu`. |
| `power_flow.qlimits.start_mode` | Symbol/String | `iteration_or_auto` | `iteration`, `auto`, `iteration_or_auto` | Activation policy. | Mixed robustness/perf runs. | Mismatched with expected policy. | Small control logic overhead. | `auto_q_delta_pu`. |
| `power_flow.qlimits.auto_q_delta_pu` | Float64 | `1e-4` | nonnegative real | Auto activation threshold. | Fine-tuning switch timing. | Extreme values. | Low. | `start_mode=auto` or `iteration_or_auto`. |
| `power_flow.qlimits.hysteresis_pu` | Float64 | `0.01` | nonnegative real | Hysteresis margin near Q limits. | Reduce switch chattering. | Too large if strict tracking needed. | Can reduce oscillatory iterations. | `cooldown_iters`, guard modes. |
| `power_flow.qlimits.cooldown_iters` | Int | `1` | nonnegative integer | Cooldown iterations after switching. | Reduce repeated toggling. | Too long cooldown on tight limits. | Affects convergence pace. | Hysteresis and freeze behavior. |
| `power_flow.qlimits.reenable_v_hyst_pu` | Float64 | `1e-4` | nonnegative real | Voltage margin of the PQ->PV release of a clamped machine: at Qmax released when `Vm > Vset + margin`, at Qmin when `Vm < Vset - margin`. | Default. | Chattering machines (raise the margin). | None. | Release also needs `hysteresis_pu > 0` or `cooldown_iters > 0`, the cooldown and the one-retry guard. |
| `power_flow.qlimits.final_q_accept_pu` | Float64 or `auto` | `auto` (`2 * hysteresis_pu`) | `auto` or real `>= hysteresis_pu` | Size bound of the final Q-limit check every enforcement mode ends with: an overshoot up to `hysteresis_pu` is within the hysteresis, up to this value bounded (accepted with a warning), beyond it a remaining violation (run not accepted). | Default. | Strict tracking (set both thresholds to 0). | None. | `hysteresis_pu`; see [Q-limits](powerlimits.md). |
| `power_flow.qlimits.trace_buses` | Vector{Int} | `[]` | bus-id vector | Trace selected bus events. | Targeted diagnostics. | Large full-network trace. | Logging overhead if populated. | Output and diagnostics verbosity. |
| `power_flow.qlimits.lock_pv_to_pq_buses` | Vector{Int} | `[]` | bus-id vector | Force listed buses into PQ-lock behavior. | Known problematic buses. | Blindly on all buses. | Can simplify switching dynamics. | Guard modes. |
| `power_flow.qlimits.guard.enabled` | Bool | `true` | `true`, `false` | Enable guard subsystem. | Prevent unstable switching. | Pure baseline comparisons. | Small runtime overhead. | Guard fields below. |
| `power_flow.qlimits.guard.min_q_range_pu` | Float64 | `0.02` | nonnegative real | Range threshold for narrow/zero detection. | Robust zero/narrow range handling. | Too high threshold. | Low. | Narrow/zero modes. |
| `power_flow.qlimits.guard.narrow_range_mode` | Symbol/String | `lock_pq` | `prefer_pq`, `lock_pq` | Action for narrow Q range units. | Convergence protection. | If strict PV control required. | Can reduce oscillations. | Hysteresis/cooldown. |
| `power_flow.qlimits.guard.zero_range_mode` | Symbol/String | `lock_pq` | `lock_pq` | Action for zero Q range units. | Deterministic limit handling. | N/A | Low. | Lock lists and violation mode. |
| `power_flow.qlimits.guard.violation_mode` | Symbol/String | `lock_pq` | `delayed_switch`, `lock_pq` | Action on persistent violations. | Robustness under bad limits. | Aggressive switching studies. | Can add control logic. | Threshold and switch caps. |
| `power_flow.qlimits.guard.violation_threshold_pu` | Float64 | `1e-4` | nonnegative real | Violation threshold. | Tune sensitivity. | Extreme values. | Low. | `violation_mode`. |
| `power_flow.qlimits.guard.max_switches` | Int | `3` | nonnegative integer | Max switches before freeze logic. | Stop chattering. | Too low on valid dynamic cases. | Can reduce wasted iterations. | `freeze_after_repeated_switching`. |
| `power_flow.qlimits.guard.max_remaining_violations` | Int | `0` | nonnegative integer | Allowed violations at guarded exit. | Controlled tolerance policies. | Strict zero-violation policies. | Low. | `accept_bounded_violations`. |
| `power_flow.qlimits.guard.accept_bounded_violations` | Bool | `false` | `true`, `false` | Permit bounded residual violations. | Practical operations tradeoff. | Strict compliance studies. | May reduce retries. | `max_remaining_violations`. |
| `power_flow.qlimits.guard.freeze_after_repeated_switching` | Bool | `true` | `true`, `false` | Freeze after repeated switch cycling. | Anti-chatter behavior. | Cases requiring unrestricted switching. | Can stabilize solves. | `max_switches`. |
| `power_flow.qlimits.guard.log` | Bool | `true` | `true`, `false` | Emit guard logs. | Diagnostics/debugging. | Quiet batch runs. | I/O overhead when enabled. | `output.console_q_limit_events`. |

## Safe configuration refresh

`refresh_sparlectra_config_file(path; write=false)` checks a user YAML
against the current template without modifying it: user values are kept,
missing keys come from `src/config/configuration.yaml.example`, duplicate
keys are reported, and `normalize_deprecated=true` rewrites known aliases.
Writing is explicit (`write=true`) after a timestamped
`.bak-YYYYmmdd-HHMMSS` backup; duplicate keys prevent automatic writes.
Sparlectra never silently rewrites user YAML during startup.

## Validated large cases

"Standard" means the default import conventions and start mode work as-is.

| Case | Size | Converges | Recommended settings |
|---|---:|:---:|---|
| `case145.m` | 145 buses | ✅ (6 it) | standard |
| `case300.m` | 300 buses | ✅ (≈5 it) | DC-seeded `profile_blend` start; note: the case's own stored reference angles are not perfectly power-balanced, so treat reference comparisons as a data diagnostic |
| `case1354pegase.m` | 1 354 buses | ✅ (3 it) | PEGASE shift convention: `matpower_import.shift_sign = -1.0`, `shift_unit = rad` |
| `case1951rte.m` | 1 951 buses | ✅ (6 it) | standard import + DC-seeded blend start; expect legitimate PV→PQ switches |
| `case2869pegase.m` | 2 869 buses | ✅ (3 it) | PEGASE shift convention (as above) |
| `case9241pegase.m` | 9 241 buses | ✅ (4 it) | PEGASE shift convention + DC-seeded `profile_blend` start |
| `case_ACTIVSg10k.m` | 10 000 buses | ✅ | DC-seeded `profile_blend` start (`model.auto_profile = recommend` selects it) |
| `case13659pegase.m` | 13 659 buses | ✅ (6 it) | start from the case's own reference (`angle_mode = matpower_va`, `voltage_mode = all_bus_vm`); the usually-robust DC blend start does **not** converge on this case |
| `case_ACTIVSg25k.m` | 25 000 buses | ✅ (8 it) | DC-seeded blend start + Q-limit guard (many zero/narrow-Q generators) |
| `case_SyntheticUSA.m` | 82 000 buses, 3 islands | ✅ | DC-seeded `profile_blend` start with autodamping; islands solve independently, `mpc.dcline` rows import as fixed injections. A plain flat start does not converge |

`model.auto_profile = recommend` picks the right convention and start
profile for all of the above except `case13659pegase.m`, which needs the
stored-reference start.

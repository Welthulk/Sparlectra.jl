# Power-Flow Configuration

The `power_flow.*` keys. The solver and its start projection, damping and
step control are described in the [Solver Guide](solver.md), which start
recipe fits which case in [Start strategies](@ref start_strategies_page),
file handling and precedence in [Central Configuration](configuration.md).

## [Solver core options](@id pf-solver-core)

The options every AC run reads. The Web UI form carries the ones a run is
usually tuned with; the rest is YAML only.

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `power_flow.method` | Symbol/String | `rectangular` | `rectangular` | AC formulation; must match `benchmark.methods`. |
| `power_flow.mode` | Symbol | `manual` | `manual`, `auto` | `manual` keeps every option as configured. `auto` inspects the network and fills the start-value, step-control and Q-limit strategy for keys you did not set (explicit keys always win); on non-convergence a bounded escalation ladder retries, the tolerance never changes. Decisions land in `auto_mode_decision.log` and the result metadata (`auto_profile`, `auto_final_stage`, `auto_final_solver`, `auto_final_strategy`, `auto_hints`). See the [Integration Guide](integration.md). Not for studies where every option must be pinned. |
| `power_flow.flatstart` | Bool | `false` | `true`, `false` | Flat start: 1.0 pu and 0 degrees, PV and slack buses at their setpoints, three-winding star points at the ratio of their lowest-impedance winding. `false` keeps the imported start voltages (MATPOWER `VM`/`VA`, CGMES `SvVoltage`, SCF `start_state`). The key sets only the profile: the start projection, `dc_seed_unconditional`, `start_current_iteration` and `apslf_start` follow their own keys and run from it, and `angle_mode`/`voltage_mode` count as `classic` while it is on (`run.log` names the override). A bare flat start needs those machines off (the template has `start_projection: true`). On CGMES runs an explicit `cgmes_import.start_values` wins, see [CGMES Import](cgmes_import.md). |
| `power_flow.tol` | Float64 | `1.0e-5` in the configuration file; the Web UI run form posts `1.0e-8` | positive real | Convergence bound for the largest single bus mismatch (infinity norm over active and reactive residuals; PV rows contribute their voltage residual). Physically `tol * baseMVA`: 1e-8 pu is 1 W at a 100 MVA base; the run log prints the equivalent. |
| `power_flow.tol_MW` | Float64 or unset | unset | positive real | The same bound in MW. When set it wins over `tol` and is converted with the case base at run time (`tol = tol_MW / baseMVA`). Values at or above 1 MW are warned about, zero or negative refused. The run form offers it next to the per-unit field. |
| `power_flow.max_iter` | Int | `80` | positive integer | Iteration cap; a run that reaches it ends as non-converged. |
| `power_flow.autodamp` | Bool | `true` | `true`, `false` | Adaptive damping of the Newton step. Off only for strict algorithm comparisons. |
| `power_flow.autodamp_min` | Float64 | `0.05` | positive real | Minimum damping factor while autodamp is active; lower values stabilize hard cases and can add iterations. |
| `power_flow.auto_slack` | Bool | `false` | `true`, `false` | Promote a reference when the case registers no slack: the unit with the best reference priority first, then external network injections, then a unit that regulates its own bus, then the largest unit (`ratedS`, `maxP`, dispatch); see [Reference priority](slack_vs_source.md#Reference-priority). The promotion is logged; without the key a network without slack aborts, which keeps a data error visible. Islands that lost their reference pick one by the same order regardless of this key. `ensureSlack!` is the underlying API; the CGMES importer applies the same ranking at import. |
| `power_flow.rescue` | Bool | `true` | `true`, `false` | After a non-converged AC solve, retry from the original start with a fixed ladder: `alternate_start` (toggle the flat-start flag), `autodamp` (skipped when already on), `dc_seed` (flat magnitudes, DC-projected angles), `settled_qlimits` (merit line search, low damping floor, Q-limit switching held back via `qlimits.start_mode = auto`). The first converging strategy wins: the result header reads `Flatstart: ... (failed, rescued)` with a `Rescue` line (the `Iterations` value is the rescue solve's), the metadata carry `rescue_used`, `rescue_strategy`, `rescue_first_flatstart`, `rescue_first_iterations`, the profile `:ac_rescue_strategy`, the Web UI shows a Rescue row. Only failed runs pay. Config-driven paths only (`runpf!(net, cfg)`, service, Web UI). Off for solver studies that must see the raw failure. Distinct from `wrong_branch_rescue`, which concerns a converged but implausible solution. |
| `power_flow.dc.fallback` | Bool | `false` | `true`, `false` | When the AC solve and the rescue ladder did not converge, run the standalone DC power flow: the net then carries DC angles and branch P flows (`vm = 1 pu`, no reactive results). The AC status stays non-converged (`erg = 1`); the profile records `:dc_fallback_applied`. Uses `power_flow.dc.*`; unlike `solver: dc` it runs only after a failed AC solve. |
| `power_flow.newton_update` | Symbol/String | `polar` | `polar`, `rectangular` | How a Newton step is applied: `polar` as magnitude and angle (MATPOWER's `newtonpf` update; keeps the magnitudes in range on the large rotations of a flat start, case9241pegase converges from the flat start only this way), `rectangular` adds the step to the complex voltage. Same Jacobian, same solution, different iterates. See [Newton Update](solver.md#newton_update). |
| `power_flow.linear_solver` | Symbol/String | `umfpack_reuse` | `umfpack`, `umfpack_reuse` | Sparse backend of the Newton step: `umfpack_reuse` analyzes the pattern once and refactors per iteration (re-analyzes on active-set pattern changes, falls back to the `umfpack` chain on a factorization error), `umfpack` analyzes every iteration. Counters `linear_solver_analyze_count`, `..._refactor_count`, `..._fallback_count` in the solver status. Independent of `power_flow.solver`; `klu` is not offered here. |
| `power_flow.power_mode` | Bool | `false` | `true`, `false` | Repeated solves on one network keep the Ybus, the sparse-LU analysis and the work arrays between solves; see [Power mode](@ref power-mode). |
| `power_flow.power_mode_lu` | Symbol/String | `auto` | `auto`, `klu`, `umfpack` | Sparse LU of power mode: `auto` chooses once per network from KLU's flop estimate of the first Jacobian (KLU up to 5e7 flops, UMFPACK above or when the KLU extension is not loaded). Read only with power mode; see [Power-mode LU](@ref power_mode_lu). |
| `power_flow.jacobian_reuse` | Bool | `false` | `true`, `false` | Dishonest Newton: the next step reuses the factorization while the mismatch falls fast; changes the iteration count, not the solution. See [Dishonest Newton](@ref dishonest_newton). |
| `power_flow.jacobian_reuse_min_reduction` | Float | `10.0` | greater than 1 | Factor by which a step must cut the maximum mismatch for the next step to reuse the factorization. |
| `power_flow.jacobian_reuse_max_steps` | Int | `3` | at least 1 | Reused steps in a row before a forced refactorization. |

## [Solver selection (rectangular vs. APSLF)](@id pf-solver-selection)

`power_flow.solver` selects the executing solver; `power_flow.method` fixes
the AC formulation (`rectangular`). `apslf` routes the run through the
external-solver bridge (`buildPfModel`, `solvePf(ApslfSolver(...))`,
`applyPfSolution!`) to AnalyticLoadFlow.jl's analytic power-series solver,
a required dependency; see [External Solvers](external_solvers.md).

```yaml
power_flow:
  solver: rectangular   # rectangular | apslf | dc
  apslf:
    order: 24
    nr_polish: false
    convergence_radius: true
  apslf_start:
    enabled: false
    order: 40
```

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `power_flow.solver` | Symbol/String | `rectangular` | `rectangular`, `apslf`, `dc` | Executing solver: Newton-Raphson, the power series, or the standalone DC power flow (next section). `apslf` skips the Newton loop and rejects `apslf_start.enabled = true`. |
| `power_flow.apslf.order` | Int | `24` | `>= 1` | Highest power-series coefficient; higher orders help stressed cases and cost more. The series is evaluated with Padé approximants. |
| `power_flow.apslf.nr_polish` | Bool | `false` | `true`, `false` | Newton-Raphson polishing step on the series result, a debugging aid; the series alone is a load-flow solution. |
| `power_flow.apslf.convergence_radius` | Bool | `true` | `true`, `false` | Evaluate the Padé margin: the distance `dmin` of the nearest Padé pole to `s = 1`, with its bus and a GRN/YEL/RED level; reported in the result header (`APSLF radius`), the metadata (`apslf_pade_margin`, `apslf_pade_margin_bus`) and on the runs page. Costs about as much as the solve. The series radius (`apslf_series_radius`, `apslf_series_radius_bus`; below 1 the plain series does not reach `s = 1`) is reported on every APSLF run independent of this key. |
| `power_flow.apslf_start.enabled` | Bool | `false` | `true`, `false` | APSLF as a start-value generator ahead of the Newton solve, see below. Rejected together with `solver: apslf`. |
| `power_flow.apslf_start.order` | Int | `40` | `>= 1` | Series order of the start-value generator; no effect while it is off. |

### [APSLF start values](@id pf-apslf-start)

With `power_flow.apslf_start.enabled` the APSLF solver runs once before
the rectangular Newton-Raphson solve, at the insertion point of the
current-iteration pre-solve and with the same guard: the candidate is
adopted only when it strictly improves the rectangular mismatch, otherwise
the original start values are restored. The generator always runs without
NR polish (the Newton solve that follows does the polishing) and without Q
limits (`power_flow.qlimits.*` governs only the Newton solve); neither is
configurable here. Requires AnalyticLoadFlow.jl; mutually exclusive with
`power_flow.solver = apslf` and `start_mode.dc_seed_unconditional`.
Diagnostic artifact: `apslf_start.log`.

## Solver selection (DC power flow)

`power_flow.solver: dc` selects the standalone DC power flow (MATPOWER
`rundcpf`/`makeBdc` equivalent): a linear model from branch series
reactance only (`B'`, no `r`, no shunts, no line charging), transformer
`phase_shift_deg` as a phase-shift injection vector, `Vm = 1.0 pu`
everywhere, no losses. [`rundcpf!`](@ref) calls it directly, independent of
`power_flow.solver`. Active controllers are rejected (the outer-loop tap and
PST controllers and the Q(U)/P(U) controllers inside the Newton step;
`apslf` runs without them and warns instead).

```yaml
power_flow:
  solver: dc
  dc:
    angle_reference_deg: 0.0
    ignore_out_of_service: true
```

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `power_flow.dc.angle_reference_deg` | Float64 | `0.0` | any real | Uniform angle offset added to every bus after the slack-referenced solve; the slack bus sits at this reference. An exact post-hoc shift, no re-solve. |
| `power_flow.dc.ignore_out_of_service` | Bool | `true` | `true` | Documents that `status == 0` branches are always excluded from `B'`; not a live toggle. |

- `rundcpf!(net; seed_ac_start = true)` (a keyword, not a YAML key) re-seeds
  the AC rectangular solve from the DC angles after a successful DC solve
  (slack and PV magnitudes restored first); `net` then holds the AC
  solution, the returned `DcPowerFlowReport` still describes the DC step,
  and the AC outcome sits in `report.metadata.ac_converged`,
  `ac_iterations`, `ac_elapsed_s`.
- DC results live in `DcPowerFlowReport` and `dc_pf_status(net)`, a
  registry separate from `rectangular_pf_status`, so a DC result is never
  mistaken for an AC one. With several AC islands,
  `ac_island_solver_summary.csv` stays empty for DC-solved islands.

## AC island diagnostics

`power_flow.islands` writes structural AC-island diagnostics for networks
with several disconnected AC components; the islands are solved
independently. No DC line becomes an AC branch and no artificial bridge is
added (`matpower_import.matpower_dcline_mode: pf_injections` stays a fixed
terminal-injection approximation). Reference selection per island:
[Reference priority](slack_vs_source.md#Reference-priority); threads:
[Parallel execution](parallel_execution.md).

```yaml
power_flow:
  islands:
    enabled: true
    parallel_min_buses: 200
    reference_policy: matpower_like
    diagnostic_continue_after_failure: true
```

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `power_flow.islands.enabled` | Bool | `true` | `true`, `false` | Write the island diagnostics before and after the solve. |
| `power_flow.islands.parallel_min_buses` | Int | `200` | at least 1 | Islands of at least this many buses are solved on their own Julia threads when two or more exist and Julia runs with more than one thread (`runtime.parallel.enabled` is the master switch); smaller islands follow serially. Results are identical to a serial solve. Configuration file only, no Web UI control. The former key `power_flow.islands.mode` is ignored with a warning. |
| `power_flow.islands.reference_policy` | Symbol/String | `matpower_like` | `matpower_like` | Keep the in-island REF bus when present; otherwise promote the bus of the island's best voltage-controlled unit, and without one its best generating unit (reference priority first, then `reference_candidate_rank`). Islands without any generating unit fail before NR. The choice and its reason are in `ac_islands.csv`. |
| `power_flow.islands.diagnostic_continue_after_failure` | Bool | `true` | `true`, `false` | Keep the diagnostics of all islands when the combined run fails. |

Artifacts in the run output directory:

| Artifact | Content |
|---|---|
| `ac_islands.csv` | Bus, branch, generator/load, DC-line terminal, power-balance, reference and status diagnostics per island (written with several islands). |
| `ac_island_solver_summary.csv` | One row per island: reference bus, PV/PQ/REF counts, solver settings, final status, final mismatch, failure reason. Islands the solver never attempted (de-energized single-bus islands, islands skipped after an earlier failure) appear with `final_status=not_attempted`, `failure_reason=not_attempted`, `stage=not_attempted`, `iterations=0`, zeroed switching statistics and `unavailable` mismatch fields. |
| `ac_island_<id>_solver.log`, `ac_island_<id>_mismatch_history.csv` | Per-island details with the same fields, written only for attempted islands. |
| `q_limit_processing_status` column | `disabled` with Q limits off, `not_attempted` for islands never solved, otherwise the solver-recorded outcome (currently `unavailable`). It never mirrors `failure_reason`; Q-limit activity is in the switching-statistics columns (`pv_pq_switching_events`, `qlimit_active_set_changes`, ...). |

When an island-aware run fails, the result message names the first failing
island: id, bus and branch counts, reference bus, PV/PQ/REF counts,
iterations, final mismatch and its status (`finite`, `nonfinite`, `NaN`,
`Inf`), failure reason, stage, start-projection setting and artifact path,
from that island's own solver record. Island-wise solving is structural
support, not a convergence guarantee: if any island fails, the run is
reported as failed and the artifacts are kept. On a multi-island case,
first run with `qlimits.enabled: false` so that the baseline tests
topology, reference selection and start projection; then enable it, and
the per-island logs add the switching events, active-set changes, guarded
narrow-range buses and the final PV voltage residual.

## [Distributed active-power slack](@id pf-distributed-slack)

`power_flow.distributed_slack` spreads the active-power imbalance of an
island (load plus losses minus scheduled generation) over the participating
generators instead of the reference bus. Theory, the augmented Newton
system and the applicability rules: [Solver Guide](solver.md).

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
| `power_flow.distributed_slack.enabled` | Bool | `false` | `true`, `false` | Master switch. Off is bit-identical to the single-slack solver; imported participation factors alone never activate it. Mutually exclusive with `external_grid.enabled` (configuration error). |
| `power_flow.distributed_slack.p_mode` | Symbol/String | `pg_weighted` | `pg_weighted`, `pmax_weighted`, `headroom_weighted`, `imported`, `explicit` | Weight source: scheduled `Pg`, `maxP`, the headroom `max(maxP - Pg, 0)`, the imported factor (MATPOWER gen column 21 `APF`, CGMES `GeneratingUnit.normalPF`), or the `weights` table. |
| `power_flow.distributed_slack.respect_p_limits` | Bool | `true` | `true`, `false` | Warn per participant whose corrected output `Pg + alpha * lambda_P` leaves `[minP, maxP]` by more than 0.01 percent of the limit (floor 1e-3 MW; smaller overshoots of units scheduled at a limit are numerical noise). Warns only, no clamp and re-solve. |
| `power_flow.distributed_slack.fallback` | Symbol/String | `error` | `error`, `ref_only` | An island without a valid participant for the chosen mode: abort, or warn and solve it classically with the reference bus absorbing everything. |
| `power_flow.distributed_slack.weights` | Mapping | `{}` | bus name or index to weight `>= 0` | Read with `p_mode: explicit` only; then non-empty with at least one positive weight. Keys resolve against bus names first, then bus indices written as strings. |

Candidates are the generator-type prosumers at the island's REF and PV
buses; fixed injections at PQ buses (HVDC converter injections, kept
boundary equivalents) never participate. Invalid candidates (missing data,
non-finite or non-positive weight) are dropped with a debug log and counted
in the metadata; the surviving weights are normalized to `sum(alpha) = 1`
per island, and in island-wise runs each island solves with its own
`lambda_P`.

| Surface | Fields |
|---|---|
| Structured status (next to the wrong-branch metadata) | `distributed_slack_active`, `distributed_slack_mode`, `distributed_slack_lambda_p_pu`, `distributed_slack_lambda_p_mw`, `distributed_slack_participants`, `distributed_slack_alpha_sum`, `distributed_slack_dropped`, `distributed_slack_p_limit_violations`, the per-participant table `distributed_slack_participation` (bus, alpha, correction `dP`, scheduled output); at `verbose > 0` a summary with the top participants is printed. |
| `printACPFlowResults` | Columns `dSl alpha` and `Pg eff MW` on participating buses (the `Pg` column keeps the schedule) and a header line with mode and `lambda_P`. |

## [External grid source](@id pf-external-grid)

`power_flow.external_grid` computes the marked slack bus as a non-ideal
external-grid source: before the solve, `convertSlackToExternalGrid!` moves
the reference voltage to a hidden internal bus `<bus>__extgrid_int` behind
the feeder impedance `z = Un^2/Sk''` (split by the R/X ratio), and the
former slack bus becomes a solved bus whose voltage droops under load.
Theory, the stiff limit and the short-circuit effect:
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
| `power_flow.external_grid.enabled` | Bool | `false` | `true`, `false` | Master switch; off keeps the ideal slack. Mutually exclusive with `distributed_slack.enabled` (configuration error): both decide who covers the island's imbalance, and combined the source would be forced to a participation share of zero and degenerate to a bare angle anchor. |
| `power_flow.external_grid.source` | Symbol/String | `auto` | `auto`, `config` | Where `Sk''` and `R/X` come from: `auto` prefers the values the case declares on the slack bus (CGMES `ExternalNetworkInjection`, SCF `source` with `sk`/`rx_ratio`; logged as declared by the case data) and falls back to the numbers below (MATPOWER and DTF carry none); `config` always uses the numbers below. |
| `power_flow.external_grid.sk_MVA` | Float | `2000.0` | `> 0` | Initial symmetrical short-circuit power of the feeder; `z_pu = baseMVA/sk_MVA` on the voltage-level base. |
| `power_flow.external_grid.rx` | Float | `0.1` | `>= 0` | R/X ratio of the feeder impedance. |

Only the primary slack bus is converted; with `multi_slack` every other
island keeps its reference. The conversion is logged in `run.log`; rescue
retries cannot stack sources. The result print names the connection in its
`Grid connection:` line and lists the internal reference bus with type
`SOURCE`. Web UI: the **External grid source** fieldset of the advanced run
options.

## [Start mode options](@id pf-start-mode)

Where the Newton-Raphson iteration starts: the angle and voltage modes, the
profile a start is taken from, and the start projection with its DC, blend
and ratio-profile candidates. A flat start (`power_flow.flatstart`) keeps
the projection and the other start machines and treats the two modes as
`classic`. Mechanism: [Solver Guide](solver.md); recipes per case:
[Start strategies](@ref start_strategies_page).

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `power_flow.start_mode.angle_mode` | Symbol/String | `dc` | `classic`, `dc`, `bus_va_blend`, `matpower_va` | Start angles: zero, a DC estimate (capped by `dc_angle_limit_deg`), the mean of the current and the imported angle, or the imported angle. Prefer the imported angles when a trusted solved state exists. |
| `power_flow.start_mode.voltage_mode` | Symbol/String | `profile_blend` | `classic`, `pv_gen_vg`, `pv_bus_vm`, `all_bus_vm`, `profile_blend` | Start magnitudes: setpoints on PV and slack buses only (`classic`, `pv_gen_vg`), the imported magnitude on those buses (`pv_bus_vm`) or on every bus (`all_bus_vm`), or the mean of setpoint and imported magnitude on every bus (`profile_blend`, with `profile_source: matpower_reference` or the SE sources). |
| `power_flow.start_mode.profile_source` | Symbol/String | `matpower_reference` | `flat`, `dc`, `bus_metadata`, `historical_profile`, `matpower_reference`, `state_estimation`, `se_snapshot`, `scada_snapshot` | Profile the blend reads; `matpower_reference` is the imported `BUS.VM`/`BUS.VA`. `state_estimation` starts from the estimated voltages of a preceding estimation while the model injections stay authoritative (the measurement/model difference goes into the slack; `runpf_from_se!(...; mode = :se_state)` after `runse!(updateNet = true)` or `readSEStateCSV!`, also the Web UI chain action). `se_snapshot` also takes the nodal balances from the estimation (working net only, the model is never mutated) and converges in 0 or 1 iterations with the slack pickup recorded in the metadata (`runpf_from_se!(...; mode = :se_snapshot)`). Both need a preceding estimation, else a clear error. |
| `power_flow.start_mode.start_projection` | Bool | `true` | `true`, `false` | Measure candidate starts and take the best; gates the sub-options below. Off only for minimal-path micro-benchmarks. |
| `power_flow.start_mode.try_dc_start` | Bool | `true` | `true`, `false` | Try the DC-angle candidate; little use on highly resistive distribution cases. |
| `power_flow.start_mode.try_blend_scan` | Bool | `true` | `true`, `false` | Try the blend candidates of `blend_lambdas`; the startup cost grows with their number. |
| `power_flow.start_mode.branch_guard` | Bool | `true` | `true`, `false` | Branch sanity guard on candidate starts. |
| `power_flow.start_mode.measure_candidates` | Bool | `true` | `true`, `false` | Score the candidates by mismatch before selection. |
| `power_flow.start_mode.accept_unmeasured_dc_start` | Bool | `false` | `true`, `false` | Accept the DC start without the measurement check (synthetic studies). |
| `power_flow.start_mode.dc_seed_unconditional` | Bool | `false` | `true`, `false` | Run the standalone DC power flow (the solver of `solver: dc`, per island) before Newton-Raphson and always take its angles, bypassing `angle_mode`, `try_dc_start` and `measure_candidates` and their quality gate: the `rundcpf!(seed_ac_start = true)` two-step inside the config-driven pipeline, with Q limits, wrong-branch detection and artifacts still applied to the AC solve. Magnitudes are untouched (slack and PV setpoints kept, the rest from `voltage_mode`). Only with `solver: rectangular`; mutually exclusive with `apslf_start.enabled`. Web UI: **Use DC start values**, which grays out the two mode selects. |
| `power_flow.start_mode.ratio_profile` | Bool | `true` | `true`, `false` | Flat starts only: the projection also measures the ratio profile, the flat profile with every PQ bus magnitude scaled by the product of the off-nominal transformer ratios on its shortest path from the reference (PV and slack keep their setpoints, angles stay 0). It replaces the flat profile under the same rule as the DC and blend candidates (smallest residual 2-norm, at least 10 percent below the flat profile); a requested DC angle start keeps precedence. Without it a flat start drives a circulating current through every off-nominal transformer, which on stiff windings (a CGMES three-winding star equivalent) can keep Newton from converging at all. Off only for studies of the bare flat start. |
| `power_flow.start_mode.reuse_import_data` | Bool | `true` | `true`, `false` | Reuse the imported MATPOWER reference values during start generation (`matpower_import.*` voltage reference keys). |
| `power_flow.start_mode.blend_lambdas` | Vector{Float64} | `[0.25, 0.5, 0.75]` | real vector, typically 0 to 1 | Blend factors of the scan. |
| `power_flow.start_mode.dc_angle_limit_deg` | Float64 | `60.0` | positive real | Magnitude cap of the DC start angles. |

## Guarded current-iteration start pre-solve

`power_flow.start_current_iteration` adds an optional guarded
current-injection pre-solve between the start projection and
Newton-Raphson. It is a start-value preconditioner, not a solver: it takes
the voltage profile prepared by start mode and projection, runs a few
damped PQ-bus current updates, and keeps the result only when it passes
the voltage and angle guards and improves the mismatch metric; otherwise
the original start values are restored. It adds no `voltage_mode` or
`angle_mode` value. Default off; enable it for cases where the normal or
the DC/profile-blend start is not enough. Diagnostic artifact:
`current_iteration_start.log`.

```text
Start voltage mode + start angle mode
-> optional start projection / candidate selection
-> optional current-iteration pre-solve
-> Newton-Raphson power flow
-> optional Q-limit handling / outer loop
```

### [Current-iteration start options](@id pf-current-iteration-start)

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

| YAML path | Type | Default | Meaning |
|---|---:|---:|---|
| `power_flow.start_current_iteration.enabled` | Bool | `false` | Enable the pre-solve. |
| `power_flow.start_current_iteration.max_iter` | Int | `10` | Maximum pre-solve steps. Keep it small; raise it only when the log shows the mismatch still improving when the pre-solve stops. |
| `power_flow.start_current_iteration.tol` | Float64 | `1.0e-3` | Stopping tolerance of the pre-solve, not the Newton tolerance; a loose value is enough, the goal is a better start, not a solved power flow. |
| `power_flow.start_current_iteration.damping` | Float64 | `0.5` | Damping of the voltage update in `(0, 1]`; 1.0 applies the full update. Lower it when the guards reject the candidate, raise it only when the pre-solve is stable but slow. |
| `power_flow.start_current_iteration.accept_only_if_improved` | Bool | `true` | Accept the candidate only when it improves the mismatch metric. Off lets a worse start enter Newton-Raphson; expert experiments only. |
| `power_flow.start_current_iteration.min_improvement_factor` | Float64 | `0.98` | Required ratio of candidate to original mismatch when `accept_only_if_improved` is on; 0.98 demands about 2 percent improvement, smaller values demand more. |
| `power_flow.start_current_iteration.vm_min_pu` | Float64 | `0.5` | Lower voltage guard: any bus below it rejects the candidate. The log shows the candidate minima. |
| `power_flow.start_current_iteration.vm_max_pu` | Float64 | `1.5` | Upper voltage guard: any bus above it rejects the candidate. The log shows the candidate maxima. |
| `power_flow.start_current_iteration.max_angle_step_deg` | Float64 | `30.0` | Largest angle change of one update; a larger jump rejects the candidate. Raise it only when the log shows plausible candidates rejected by this guard alone. |
| `power_flow.start_current_iteration.only_for_large_cases` | Bool | `false` | Run the pre-solve only for cases above the large-case threshold of the rectangular workspace, so small examples stay unchanged. |

The pre-solve runs before the Newton solve; the classic Q-limit outer-loop
modes apply it to the first inner solve only, and it never changes a
switching decision. A rejected candidate leaves Newton-Raphson on the
original start values. It does not help cases whose problem is a model or
convention mismatch, a wrong branch, or an invalid MATPOWER import
convention.

### Diagnostics and interpretation

`current_iteration_start.log` (written when the run has an output
directory in its performance profile):

| Group | Fields |
|---|---|
| Outcome | `current_iteration_enabled`, `current_iteration_attempted`, `current_iteration_accepted`, `current_iteration_reason` |
| Mismatch | `initial_mismatch`, `final_mismatch`, `iterations` |
| Candidate voltages | `candidate_voltage_magnitude_min`, `candidate_voltage_magnitude_max`, `candidate_voltage_low_count`, `candidate_voltage_high_count`, `candidate_voltage_worst_low_bus`, `candidate_voltage_worst_high_bus` with their values |
| Candidate angles | `candidate_max_angle_step_deg`, `maximum_angle_step_deg` |
| Rejection | `guard_violations`, `rejection_stage`, `rejected_at_iteration`, `original_start_values_restored` |
| Restored ranges | `restored_voltage_magnitude_min`, `restored_voltage_magnitude_max` |

`current_iteration_accepted: true` means Newton-Raphson started from the
candidate; `false` with `original_start_values_restored: true` means it
started from the original values. `current_iteration_reason` names why:
`voltage_magnitude_guard` (a voltage outside `vm_min_pu`/`vm_max_pu`),
`angle_step_guard` (a step beyond `max_angle_step_deg`), `not_improved`,
or `disabled`, `skipped_small_case`, `max_iter`, `tolerance_reached`,
`invalid_voltage`, `singular_current_update`, `invalid_mismatch`.

## [Merit-function line search options](@id pf-merit)

`power_flow.merit` adds an Armijo sufficient-decrease test to the autodamp
backtracking loop of the rectangular solver; it replaces neither
Newton-Raphson nor autodamp. Off by default, which leaves the max-mismatch
autodamp unchanged. Requires `power_flow.autodamp = true` (validation error
otherwise). Theory: [Merit-Function Line Search](@ref).

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

| YAML path | Type | Default | Allowed | Meaning |
|---|---:|---:|---|---|
| `power_flow.merit.enabled` | Bool | `false` | `true`, `false` | Master switch. For difficult flat-start cases where the infinity-norm autodamp criterion accepts a step that lowers the worst-bus mismatch but raises the overall residual. One extra weighted residual norm per trial mismatch. |
| `power_flow.merit.armijo_c1` | Float64 | `1.0e-4` | real in `(0, 0.5)` | Sufficient-decrease constant of the Armijo condition; larger values reject more trials. |
| `power_flow.merit.scale_p` | Float64 | `1.0` | positive real | Weight of the active-power residuals in the merit function. YAML only. |
| `power_flow.merit.scale_q` | Float64 | `1.0` | positive real | Weight of the reactive-power residuals (PQ buses). YAML only. |
| `power_flow.merit.scale_v` | Float64 | `1.0` | positive real | Weight of the voltage-setpoint residuals (PV buses). YAML only. |
| `power_flow.merit.fallback_max_mismatch` | Bool | `true` | `true`, `false` | When no trial satisfies Armijo: `true` falls back to the max-mismatch criterion (first improving trial, else the best finite trial), `false` goes straight to the best finite trial, which can mean smaller steps and more iterations. |

Leave the scales at `1.0` unless the P, Q and V residuals differ by orders
of magnitude.

### Diagnostics and interpretation

| Surface | Fields |
|---|---|
| `merit_linesearch.log` (one line per Newton iteration, written when the run has an output directory in its performance profile) | `f_before`, `directional_derivative`, `tested_alphas`, `accepted_alpha`, `accept_reason` |
| Solver status | `merit_enabled`, `merit_used_iterations`, `merit_fallback_count`, `merit_active_set_skip_count`, `merit_initial`, `merit_final` |

`accept_reason` values: `armijo` (a trial satisfied the Armijo condition),
`fallback_max_mismatch` (none did, and the max-mismatch criterion found an
improving trial), `fallback_conservative` (no trial satisfied Armijo or
improved the max mismatch; the most conservative finite trial was taken),
`active_set_skip` (a PV/PQ switch in this iteration changed the meaning of
the residual entries, so the merit comparison was skipped and the
max-mismatch criterion used).

## [Trust-region step control options](@id pf-trust-region)

`power_flow.trust_region` is a scaled-Newton trust region, an alternative
to `power_flow.autodamp`: it caps the Newton step norm at an adaptive
radius and accepts or rejects trials by merit decrease. Off by default;
requires `power_flow.autodamp = false` (both control the step length,
enabling both is a configuration error). Theory and the dogleg mode:
[Trust-Region Step Control](@ref).

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

| YAML path | Type | Default | Allowed | Meaning |
|---|---:|---:|---|---|
| `power_flow.trust_region.enabled` | Bool | `false` | `true`, `false` | Master switch, for difficult flat-start cases where a merit-decrease rule beats max-mismatch backtracking. One extra weighted-residual evaluation and a sparse matrix-vector product per trial on the existing Jacobian, no extra factorization. |
| `power_flow.trust_region.step_mode` | Symbol/String | `scaled` | `scaled`, `dogleg` | `scaled` rescales the full Newton direction to the radius; `dogleg` blends it with a steepest-descent (Cauchy) step when the radius shrinks below the Newton step norm, for cases where the Newton direction degrades partway through the solve. `dogleg` costs two more sparse matrix-vector products per trial and brings nothing where `scaled` already converges. |
| `power_flow.trust_region.initial_radius` | Float64 | `1.0` | positive, `> min_radius`, `<= max_radius` | Starting radius, per-unit 2-norm of the state vector; larger values risk more rejected first steps. |
| `power_flow.trust_region.min_radius` | Float64 | `1.0e-4` | positive, `< initial_radius` | Radius floor; below it the run ends with `reason = :trust_region_collapsed` instead of looping. Near `initial_radius` it declares collapse too early. |
| `power_flow.trust_region.max_radius` | Float64 | `10.0` | positive, `>= initial_radius` | Radius ceiling. |
| `power_flow.trust_region.eta_accept` | Float64 | `0.1` | positive real | Minimum actual/predicted reduction ratio `rho` to accept a trial; higher values reject more. |
| `power_flow.trust_region.shrink_factor` | Float64 | `0.5` | `(0, 1)` | Radius multiplier on a rejected trial, applied until a trial is accepted or the radius collapses. |
| `power_flow.trust_region.expand_factor` | Float64 | `2.0` | `> 1` | Radius multiplier on an accepted step that hit the radius with `rho` above `expand_threshold`. |
| `power_flow.trust_region.expand_threshold` | Float64 | `0.75` | `(0, 1)` | `rho` above which such a step expands the radius; should stay `>= eta_accept` (not validated). |

### Diagnostics and interpretation

| Surface | Fields |
|---|---|
| `trust_region.log` (one line per Newton iteration, written when the run has an output directory in its performance profile) | `radius_before`, `rho`, `tested_radii` (the shrink sequence of this iteration), `rejected_steps`, `accepted`, `radius_after`, `collapsed`, `accept_reason` |
| Solver status | `trust_region_enabled`, `tr_step_count`, `tr_rejected_steps`, `tr_min_radius`, `tr_max_radius` (the observed range, not the configured bounds), `tr_final_radius`, `tr_collapsed`, and the dogleg counters `tr_dogleg_newton_count`, `tr_dogleg_interp_count`, `tr_dogleg_cauchy_count`, `tr_active_set_skip_count` (all `0` in `scaled` mode) |

`accept_reason` values: `scaled` (the only value of the `scaled` mode);
`dogleg_newton`, `dogleg_interp`, `dogleg_cauchy` (the accepted point on
the dogleg path); `active_set_skip` (a PV/PQ switch this iteration,
`dogleg` only: the comparison was skipped and one `scaled`-style trial
accepted unconditionally); `none` (collapse, no step accepted).
`tr_collapsed: true` means the radius fell below `min_radius` in some
iteration without an accepted step: lower `min_radius`, raise
`initial_radius`, or use autodamp for this case. A high `tr_rejected_steps`
relative to `tr_step_count` says the model predicts poorly for this case;
try another start profile.

## Step-control mode combinations (autodamp / merit / trust-region)

`power_flow.autodamp`, `power_flow.merit.enabled` and
`power_flow.trust_region.enabled` jointly select the step length.
`PowerFlowConfig` enforces two rules at load time (`ArgumentError`) for
YAML, API `config_overrides` and the Web UI alike: merit needs
`autodamp = true` (it is a test inside the autodamp loop), trust region
needs `autodamp = false` (it replaces that loop). So merit and trust region
never hold together, `autodamp = true` alone is the classic max-mismatch
backtracking, and `autodamp = false` without trust region is a
fixed-damping Newton step. The Web UI mirrors the rules: the autodamp and
merit box and the trust-region box exclude each other, and the merit
toggle is disabled while autodamp is off. `trust_region.step_mode` is
independent of the rules.

## [Q-limit options and guard](@id pf-qlimits)

How the reactive-power limits of PV machines are enforced: master switch,
enforcement mode, when switching starts, hysteresis and release, the final
check every mode ends with, and the guard for narrow ranges and remaining
violations. Algorithm: [Q-limits](powerlimits.md) and
[Q-limit Switching Strategy](q_limit_switching_strategy.md). In the Web UI
these keys form the Q-limit block of the power-flow form
([Form options](webui_reference.md#webui-form-options)).

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `power_flow.qlimits.enabled` | Bool | `true` | `true`, `false` | Master switch; off is an unconstrained power flow. |
| `power_flow.qlimits.enforcement_mode` | Symbol/String | `active_set` | `active_set`, `classic_simultaneous`, `classic_one_at_a_time`; `off` disables the handling like `enabled: false` | `active_set` switches PV to PQ inside the iteration. The classic modes run an outer loop (solve with switching off, clamp violating generators to the limit, convert their buses to PQ, rerun without PQ to PV re-enable), all violations at once or the largest per pass, and may need repeated solves. The legacy aliases `matpower_simultaneous` and `matpower_one_at_a_time` map to the `classic_*` values. |
| `power_flow.qlimits.start_iter` | Int | `3` | integer | First Newton iteration from which switching may run; the classic modes switch between solves and ignore it. |
| `power_flow.qlimits.start_mode` | Symbol/String | `iteration_or_auto` | `iteration`, `auto`, `iteration_or_auto` | `iteration` starts at `start_iter`, `auto` once the Q change per iteration falls below `auto_q_delta_pu`, `iteration_or_auto` at whichever comes first, so a small `start_iter` always wins there. On large systems whose early switching destabilises the iteration prefer `auto`. |
| `power_flow.qlimits.auto_q_delta_pu` | Float64 | `1e-4` | nonnegative real | Q change per iteration below which the `auto` rule lets switching begin. |
| `power_flow.qlimits.hysteresis_pu` | Float64 | `0.01` | nonnegative real | Margin beyond a limit before a machine switches, and the overshoot the final check still counts as within limits. Larger values damp chattering. |
| `power_flow.qlimits.cooldown_iters` | Int | `1` | nonnegative integer | Iterations a bus waits after a switch before it may switch again. |
| `power_flow.qlimits.reenable_v_hyst_pu` | Float64 | `1e-4` | nonnegative real | Voltage margin of the PQ to PV release of a clamped machine: at Qmax when `Vm > Vset + margin`, at Qmin when `Vm < Vset - margin`. A release also needs `hysteresis_pu > 0` or `cooldown_iters > 0`, the cooldown and the one-retry guard. Raise it for chattering machines. |
| `power_flow.qlimits.final_q_accept_pu` | Float64 or `auto` | `auto` (`2 * hysteresis_pu`) | `auto` or real `>= hysteresis_pu` | Bound of the final Q-limit check every mode ends with: an overshoot up to `hysteresis_pu` is within the hysteresis, up to this value bounded (accepted with a warning), beyond it a remaining violation (run not accepted). Both thresholds at 0 for strict tracking. |
| `power_flow.qlimits.classic_max_passes` | Int | `30` | integer `>= 1` | Pass limit of the classic outer loop: solves after the base solve, one per Q-limit update. Violations found after the last pass stop the loop with `max_outer_iterations` (metadata `q_limit_classic_outer_loop_stop`, `run.log`); the run is then not converged, the last update is not applied, the net keeps the state of the last solve, and `q_limit_classic_outer_loop_passes` counts the passes solved. |
| `power_flow.qlimits.trace_buses` | Vector{Int} | `[]` | bus numbers of the case file (an internal position where no bus carries that number) | Write the switching events of these buses in detail to the Q-limit trace. |
| `power_flow.qlimits.lock_pv_to_pq_buses` | Vector{Int} | `[]` | internal bus positions (1 = first bus); the file-based MATPOWER path maps case bus numbers to positions | Run the listed PV buses as PQ from the start, at their scheduled Q clamped into their limits; logged at iteration 0, never released back to PV. Read by every enforcement mode of the Newton solver, ignored with Q limits off and by APSLF. Auto mode adds the buses that switched too often. |
| `power_flow.qlimits.guard.enabled` | Bool | `true` | `true`, `false` | Before the solve, run machines with a narrow or zero Q range as PQ (the range rules below). The switch cap, freezing, the violation rule and the bounded-violation acceptance act without it. |
| `power_flow.qlimits.guard.min_q_range_pu` | Float64 | `0.02` | nonnegative real | Q range below which a PV bus counts as narrow. The range is judged per PV bus after every in-service machine of the bus is aggregated (Qmin and Qmax summed), so a zero-range machine next to a wide-range one on the same bus is not guarded; buses the case types PQ are never considered, and the reported count is PV buses, not machines. The keyword `qlimit_guard_min_q_range_pu` of `runpf_rectangular!` has the same default. |
| `power_flow.qlimits.guard.narrow_range_mode` | Symbol/String | `lock_pq` | `prefer_pq`, `lock_pq` | Rule for a narrow range: both values run the machine as PQ at the middle of its range from the start of the solve. |
| `power_flow.qlimits.guard.zero_range_mode` | Symbol/String | `lock_pq` | `lock_pq` | Rule for `Qmin == Qmax`: run the machine as PQ at that Q. |
| `power_flow.qlimits.guard.violation_mode` | Symbol/String | `lock_pq` | `delayed_switch`, `lock_pq` | When a violating machine switches to PQ: `lock_pq` once its Q exceeds the limit by `violation_threshold_pu`, `delayed_switch` once it exceeds the limit by `hysteresis_pu`. |
| `power_flow.qlimits.guard.violation_threshold_pu` | Float64 | `1e-4` | nonnegative real | Overshoot from which the `lock_pq` violation rule switches; unused under `delayed_switch`. |
| `power_flow.qlimits.guard.max_switches` | Int | `3` | integer `>= 1` (0 acts as 1) | Switches of one bus after which it counts as oscillating. With freezing on, a converged run that reaches it ends as `max_switching_exceeded`; the message names this key, its value and the buses. |
| `power_flow.qlimits.guard.max_remaining_violations` | Int | `0` | nonnegative integer | PV buses that may still violate a limit at the end when bounded violations are accepted. |
| `power_flow.qlimits.guard.accept_bounded_violations` | Bool | `false` | `true`, `false` | Accept a converged run with up to that many remaining violations instead of rejecting it. |
| `power_flow.qlimits.guard.freeze_after_repeated_switching` | Bool | `true` | `true`, `false` | Freeze a bus that reached `max_switches` in its current type. |
| `power_flow.qlimits.guard.log` | Bool | `true` | `true`, `false` | Print how many narrow-range machines the guard locked as PQ (verbose console output, `output.console_q_limit_events`). |

The defaults are those of the packaged template; a library call without a
configuration file (`runpf!(net, ...)`, `run_sparlectra(net = ...)` with
`SparlectraConfig()`) and the keyword defaults of the solver entry points
use the same values. A library caller who needs the guard off sets
`qlimits.guard.enabled: false` (or `qlimit_guard = false`).

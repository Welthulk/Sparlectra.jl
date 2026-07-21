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

## Solver core options

| YAML path | Type | Default | Allowed values | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `power_flow.method` | Symbol/String | `rectangular` | `rectangular` | AC solver formulation. | Always (current core). | N/A | Fixed single implementation. | Must match `benchmark.methods`. |
| `power_flow.sparse` | Bool | `true` | `true` | Sparse linear algebra mode. | Always (required). | N/A | Scales better for large systems. | Validated together with `method`. |
| `power_flow.flatstart` | Bool | `true` | `true`, `false` | Legacy flatstart toggle. | Synthetic starts. | When start projection and imported starts are used. | Low. | Combined into `start_mode`. |
| `power_flow.tol` | Float64 | `1.0e-5` | positive real | PF tolerance. | Accuracy-sensitive studies. | Overly tight on big batches. | Tighter means more iterations. | `max_iter`. |
| `power_flow.max_iter` | Int | `80` | positive integer | Iteration cap. | Hard cases. | Very low values. | Upper runtime bound. | `tol`, `qlimits`. |
| `power_flow.autodamp` | Bool | `true` | `true`, `false` | Adaptive damping. | Difficult convergence. | Strict algorithm comparison. | Small overhead, often fewer failures. | `autodamp_min`. |
| `power_flow.autodamp_min` | Float64 | `0.05` | positive real | Minimum damping factor. | Stabilizing hard cases. | Near-zero damping on easy grids. | Lower can increase iterations. | Active only with `autodamp=true`. |

## Solver selection (rectangular vs. APSLF)

`power_flow.solver` selects which solver actually computes the power-flow
solution; it is independent of `power_flow.method`, which fixes the AC
formulation family (currently always `rectangular`). `apslf` routes the run
through the external-solver bridge (`buildPfModel` → `solvePf(ApslfSolver(...))`
→ `applyPfSolution!`) to AnalyticLoadFlow.jl's analytic power-series solver
instead of the internal Newton-Raphson loop. It requires the optional
AnalyticLoadFlow.jl dependency to be loaded (`using AnalyticLoadFlow`); see
`apslf_solver`.

```yaml
power_flow:
  solver: rectangular   # rectangular | apslf
  apslf:
    order: 40
    use_pade: true
    nr_polish: true
  apslf_start:
    enabled: false
    order: 40
```

| YAML path | Type | Default | Allowed values | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `power_flow.solver` | Symbol/String | `rectangular` | `rectangular`, `apslf` | Selects the executing solver. | `apslf` to use the analytic power-series solver instead of NR. | `apslf` without AnalyticLoadFlow.jl loaded (raises a clear error). | `apslf` skips the NR iteration loop entirely. | Rejects `apslf_start.enabled=true` when set to `apslf`. |
| `power_flow.apslf.order` | Int | `40` | `>= 1` | Highest power-series coefficient computed. | Higher for stressed/large-angle cases. | Unnecessarily high orders on easy cases (cost). | Higher order increases solve cost. | `use_pade`. |
| `power_flow.apslf.use_pade` | Bool | `true` | `true`, `false` | Evaluate the voltage series via Padé `[L/M]` approximants instead of direct Taylor summation. | Default; improves convergence radius. | Direct Taylor comparison studies. | Small evaluation overhead, usually better accuracy per order. | `order`. |
| `power_flow.apslf.nr_polish` | Bool | `true` | `true`, `false` | Run a Newton-Raphson polishing step on the series result. | Default; tightens the final residual. | Pure series-accuracy studies. | Adds a small number of NR iterations. | Only active when `solver=apslf`. |
| `power_flow.apslf_start.enabled` | Bool | `false` | `true`, `false` | Use the APSLF solver as a start-value generator ahead of the rectangular NR solve (guarded, like `start_current_iteration`). | Difficult NR starts. | `solver=apslf` (rejected: start generator only makes sense ahead of NR). | Adds one series solve before NR. | `solver`, `apslf_start.order`. |
| `power_flow.apslf_start.order` | Int | `40` | `>= 1` | Series order used by the start-value generator. | Same considerations as `apslf.order`. | Unnecessarily high orders for a start-only pass. | Higher order increases pre-solve cost. | `apslf_start.enabled`. |

`power_flow.apslf_start` deliberately has no `use_pade`/`nr_polish` fields:
polishing is left to the downstream NR solve, so the generator always runs
with `nr_polish=false` internally. It also always runs unconstrained
(`Qmin`/`Qmax` are not passed to the series solve), independent of
`power_flow.qlimits.enabled` or any other Q-limit setting: `power_flow.qlimits.*`
only governs the rectangular NR solve that follows, not the guarded start-value
candidate itself. This is a fixed internal behavior, not a separate
configurable option — the generator's only job is producing a better starting
voltage profile; Q-limit enforcement is entirely the downstream NR solve's
responsibility.

## AC island diagnostics

`power_flow.islands` enables structural AC-island diagnostics for imported or
constructed networks that contain multiple disconnected AC components. Island
diagnostics do not model DC lines as AC branches and do not add artificial
admittance bridges. They report the AC topology that the existing rectangular
solver sees after import; `matpower_import.matpower_dcline_mode:
pf_injections` remains a fixed terminal-injection approximation only.

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
| `power_flow.islands.mode` | Symbol/String | `solve_independent` | `solve_independent` | Reserved island solve mode; diagnostics preserve the island-wise solving concept. |
| `power_flow.islands.reference_policy` | Symbol/String | `matpower_like` | `matpower_like` | Select the in-island REF bus when present, otherwise report the first PV/PQ bus that would be promoted for diagnostics. |
| `power_flow.islands.diagnostic_continue_after_failure` | Bool | `true` | `true`, `false` | Keep diagnostics for all detected islands even when the combined run fails. |

When an island-aware run fails, the API/result message identifies the first
failing island with its island id, bus and branch counts, selected reference
bus, PV/PQ/REF counts, iteration count, final mismatch, mismatch status
(`finite`, `nonfinite`, `NaN`, or `Inf`), failure reason, stage
(`before_nr`, `during_nr`, or `after_nr_validation`), start-projection setting,
and the diagnostic artifact path.

Each run output directory receives:

- `ac_island_solver_summary.csv`: one row per detected AC island, including the
  selected reference bus, PV/PQ/REF counts, propagated solver settings, final
  status, final mismatch, and failure reason.
- `ac_island_<id>_solver.log`: compact per-island details with the same
  topology, setting, and solve-status fields.

Island-wise solving is structural support, not a convergence guarantee. A large
island can fail independently while smaller islands are structurally valid. If
any island or the combined solve fails, Sparlectra still reports the combined
run as failed and preserves the island artifacts for diagnosis.

For pure SyntheticUSA island diagnostics, start with Q-limit handling disabled
so the baseline tests island topology, reference selection, start projection,
and rectangular NR behavior before adding active-set effects:

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

After this pure baseline is understood, enable `power_flow.qlimits.enabled` to
diagnose active-set behavior. Per-island logs then include Q-limit state such as
`q_limit_processing_status`, switching-event counts, active-set changes,
reenable events, guarded narrow-range PV buses, final PV voltage residual, and
available mismatch metrics.

## Start mode options

| YAML path | Type | Default | Allowed values | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `power_flow.start_mode.angle_mode` | Symbol/String | `dc` | `classic`, `dc`, `bus_va_blend`, `matpower_va` | Angle initialization source mode. | Transmission-like starts. | If trusted measured/historical state exists. | Can reduce iterations. | `try_dc_start`, `dc_angle_limit_deg`. |
| `power_flow.start_mode.voltage_mode` | Symbol/String | `profile_blend` | `classic`, `pv_gen_vg`, `pv_bus_vm`, `all_bus_vm`, `profile_blend` | Voltage-magnitude/angle blend strategy. | Imported-reference assisted starts. | Untrusted imported data. | Small startup overhead. | `blend_lambdas`, `reuse_import_data`. |
| `power_flow.start_mode.profile_source` | Symbol/String | `matpower_reference` | `flat`, `dc`, `bus_metadata`, `historical_profile`, `matpower_reference`, `state_estimation`, `scada_snapshot` | Source for external or model-derived start profiles. `matpower_reference` means imported `BUS.VM`/`BUS.VA` values for regression/benchmarking/import reproduction. | Profile-aware starts. | When source data is unavailable or untrusted. | Minimal parsing overhead. | Keep source explicit for diagnostics and reproducibility. |
| `power_flow.start_mode.start_projection` | Bool | `true` | `true`, `false` | Enables start-projection workflow. | Robustness on hard cases. | Minimal-path microbench runs. | Extra startup pass. | Gates start-projection sub-options. |
| `power_flow.start_mode.try_dc_start` | Bool | `true` | `true`, `false` | Try DC candidate start. | Large transmission cases. | Highly resistive distribution cases. | Low overhead. | `dc_angle_limit_deg`. |
| `power_flow.start_mode.try_blend_scan` | Bool | `true` | `true`, `false` | Scan blend candidates. | Mixed/tricky start data. | Easy cases needing speed. | Startup cost ∝ lambda count. | `blend_lambdas`. |
| `power_flow.start_mode.branch_guard` | Bool | `true` | `true`, `false` | Branch sanity guard for candidate starts. | Stability-focused runs. | Rarely disabled. | Low. | Candidate measurement options. |
| `power_flow.start_mode.measure_candidates` | Bool | `true` | `true`, `false` | Score/select among candidates. | Multiple start candidates. | Fastest startup path only. | Low/medium startup overhead. | `try_dc_start`, `try_blend_scan`. |
| `power_flow.start_mode.accept_unmeasured_dc_start` | Bool | `false` | `true`, `false` | Allow DC start without measurement checks. | Synthetic studies. | Measurement-driven workflows. | Can avoid fallback retries. | `try_dc_start`. |
| `power_flow.start_mode.reuse_import_data` | Bool | `true` | `true`, `false` | Reuse imported MATPOWER references. | Trusted imports. | Uncertain conversion data. | Small reduction in recomputation. | `matpower_import.*` voltage reference keys. |
| `power_flow.start_mode.blend_lambdas` | Vector{Float64} | `[0.25,0.5,0.75]` | real vector (typ. 0..1) | Lambda candidates for blend scan. | Need robust candidate search. | Very large lambda sets. | Linear startup growth with vector size. | `try_blend_scan`. |
| `power_flow.start_mode.dc_angle_limit_deg` | Float64 | `60.0` | positive real | DC-start angle magnitude cap (deg). | Conservative angle starts. | Overly restrictive values. | Negligible. | `try_dc_start`. |


## Guarded current-iteration start pre-solve

`power_flow.start_current_iteration` enables an optional guarded current-injection/current-iteration pre-solve for start values. It is a start-value preconditioner, not a new power-flow solver. The final AC power-flow solve remains the rectangular Newton-Raphson path.

The rectangular power-flow start sequence is:

```text
Start Voltage Mode + Start Angle Mode
→ optional start projection / candidate selection
→ optional current-iteration pre-solve
→ Newton-Raphson power flow
→ optional Q-limit handling / outer loop logic
```

Current iteration does not introduce a new `start_voltage_mode`, `start_angle_mode`, `power_flow.start_mode.voltage_mode`, or `power_flow.start_mode.angle_mode` value. It consumes the voltage profile prepared by the existing start-mode and start-projection settings, then attempts a limited PQ-bus current update before Newton-Raphson starts. The candidate is used only when the guarded checks accept it. If a guard rejects the candidate or the mismatch does not improve enough, Sparlectra restores the original start values before entering Newton-Raphson.

### Configuration block

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
| `power_flow.start_current_iteration.enabled` | Bool | `false` | Enable the optional guarded current-iteration pre-solve. |
| `power_flow.start_current_iteration.max_iter` | Int | `10` | Maximum number of current-iteration update steps. |
| `power_flow.start_current_iteration.tol` | Float64 | `1.0e-3` | Stops the pre-solve early when the rectangular mismatch of the candidate is at or below this tolerance. |
| `power_flow.start_current_iteration.damping` | Float64 | `0.5` | Damping factor in `(0, 1]` applied to each candidate current update. |
| `power_flow.start_current_iteration.accept_only_if_improved` | Bool | `true` | Require the best candidate mismatch to improve relative to the original start mismatch. |
| `power_flow.start_current_iteration.min_improvement_factor` | Float64 | `0.98` | Acceptance threshold multiplier when improvement checking is active; the best mismatch must be at most `initial_mismatch * min_improvement_factor`. |
| `power_flow.start_current_iteration.vm_min_pu` | Float64 | `0.5` | Lower voltage-magnitude guard for candidate voltages. |
| `power_flow.start_current_iteration.vm_max_pu` | Float64 | `1.5` | Upper voltage-magnitude guard for candidate voltages. |
| `power_flow.start_current_iteration.max_angle_step_deg` | Float64 | `30.0` | Maximum allowed candidate angle step in degrees for a single current-iteration update. |
| `power_flow.start_current_iteration.only_for_large_cases` | Bool | `false` | Attempt the pre-solve only when the bus count reaches the implemented large-case threshold used by the rectangular workspace configuration; smaller cases are skipped with reason `skipped_small_case`. |

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

This combination may reduce the initial mismatch for difficult imported cases, but it is experimental and is not guaranteed to rescue a non-converging case. If the candidate violates guards or does not improve the mismatch enough, the pre-solve is rejected and Newton-Raphson starts from the original start values. It does not replace MATPOWER auto-profile selection, DC angle starts, or start projection.

### Diagnostics and interpretation

When a run has an output directory in its performance profile, the pre-solve writes `current_iteration_start.log`. The artifact records:

- `current_iteration_enabled`, `current_iteration_attempted`, `current_iteration_accepted`, and `current_iteration_reason`.
- `initial_mismatch`, `final_mismatch`, and `iterations`.
- Candidate voltage diagnostics: `candidate_voltage_magnitude_min`, `candidate_voltage_magnitude_max`, `candidate_voltage_low_count`, `candidate_voltage_high_count`, `candidate_voltage_worst_low_bus`, `candidate_voltage_worst_high_bus`, and their corresponding values.
- Candidate angle diagnostics: `candidate_max_angle_step_deg` and `maximum_angle_step_deg`.
- Rejection diagnostics: `guard_violations`, `rejection_stage`, `rejected_at_iteration`, and `original_start_values_restored`.
- Restored voltage ranges after rejection: `restored_voltage_magnitude_min` and `restored_voltage_magnitude_max`.

Important interpretations:

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

Other implementation reasons include `disabled`, `skipped_small_case`, `max_iter`, `tolerance_reached`, `invalid_voltage`, `singular_current_update`, and `invalid_mismatch`.

### Q-limit interaction

For the active-set path, the guarded current-iteration pre-solve is attempted before the Newton-Raphson solve when enabled. For classical MATPOWER-style Q-limit outer-loop modes, it is applied only to the first inner solve; later inner solves in the same outer loop start without repeating current iteration. Current iteration does not directly change Q-limit switching decisions.

### Limitations

- Experimental diagnostic start-value helper.
- Can be rejected by voltage, angle-step, finite-value, singular-update, or mismatch-improvement guards.
- Can reduce the initial mismatch without guaranteeing Newton-Raphson convergence.
- May not help cases whose main issue is model/convention mismatch, wrong branch, or an invalid MATPOWER import convention.
- Does not replace MATPOWER auto-profile, DC start, or start projection.

## Merit-function line search options

`power_flow.merit` enables an optional Armijo sufficient-decrease acceptance criterion inside the existing autodamp backtracking loop of the rectangular Newton-Raphson solver. It is an alternative, opt-in acceptance test — not a replacement for Newton-Raphson and not a replacement for autodamp. Disabled by default (`enabled: false`), which leaves today's max-mismatch autodamp behavior byte-for-byte unchanged. See [Merit-Function Line Search](@ref) in the solver documentation for the theoretical background.

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

Enable `power_flow.merit.enabled` together with `power_flow.autodamp` on cases where the historical ∞-norm autodamp criterion is suspected to accept a step that reduces the worst-bus mismatch while increasing overall residual energy. Leave `scale_p`/`scale_q`/`scale_v` at `1.0` unless P/Q/V residuals are known to differ by orders of magnitude in a particular case.

### Diagnostics and interpretation

When a run has an output directory in its performance profile and `merit.enabled = true`, the solver writes `merit_linesearch.log` with one line per Newton iteration: `f_before`, `directional_derivative`, `tested_alphas`, `accepted_alpha`, and `accept_reason`. The solver status is extended with `merit_enabled`, `merit_used_iterations`, `merit_fallback_count`, `merit_active_set_skip_count`, `merit_initial`, and `merit_final`.

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

### Limitations

- Opt-in line-search criterion, not a general-purpose nonlinear optimizer; it does not replace Newton-Raphson.
- The merit function `f(x) = 1/2 ‖W F(x)‖²` can have local minima with `F(x) ≠ 0`; satisfying Armijo does not by itself guarantee convergence.
- Does not influence which of multiple numerical solution branches (e.g. high-voltage vs. low-voltage) the solver converges to.
- PV/PQ active-set switches make `f` discontinuous in meaning across the switch; the merit comparison is skipped for that iteration (`accept_reason = active_set_skip`).
- Does not change candidate start-value ranking (`start_mode.measure_candidates`), which remains mismatch-based.

## Trust-region step control options

`power_flow.trust_region` enables an optional scaled-Newton trust-region alternative to `power_flow.autodamp`: it caps the Newton step norm at an adaptive radius and accepts/rejects trials by merit-function decrease rather than by the max-mismatch criterion. Disabled by default (`enabled: false`). **Mutually exclusive with `power_flow.autodamp = true`** — both mechanisms control the Newton step length; enabling both is a configuration error. See [Trust-Region Step Control](@ref) in the solver documentation for the theoretical background.

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
```

| YAML path | Type | Default | Allowed | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `power_flow.trust_region.enabled` | Bool | `false` | `true`, `false` | Master switch for scaled-Newton trust-region step control. | Difficult flat-start cases where a merit-decrease step-acceptance rule (rather than max-mismatch backtracking) is preferred. | Together with `power_flow.autodamp = true` (validation error). | One extra weighted-residual evaluation and a sparse matrix-vector product per trial, reusing the already-built Jacobian; no extra factorization. | Requires `power_flow.autodamp = false`. |
| `power_flow.trust_region.initial_radius` | Float64 | `1.0` | positive, `> min_radius`, `<= max_radius` | Starting trust-region radius, in per-unit state-vector (2-)norm. | Tuning how large the first trial step may be. | N/A | Larger values risk more rejected/shrunk first steps on hard cases. | Bounded by `min_radius`/`max_radius`. |
| `power_flow.trust_region.min_radius` | Float64 | `1.0e-4` | positive, `< initial_radius` | Radius floor. Falling below it declares non-convergence (`reason = :trust_region_collapsed`) instead of looping indefinitely. | Bounding worst-case retry effort. | Setting it near `initial_radius`, which makes collapse detection overly aggressive. | Smaller values allow more shrink retries before giving up. | Read every retry inside one Newton iteration. |
| `power_flow.trust_region.max_radius` | Float64 | `10.0` | positive, `>= initial_radius` | Radius ceiling; the radius never expands past this value. | Preventing runaway expansion on well-behaved cases. | N/A | Higher ceiling allows larger accepted steps once the model is trusted. | Caps the `expand_factor` growth. |
| `power_flow.trust_region.eta_accept` | Float64 | `0.1` | positive real | Minimum actual/predicted reduction ratio `rho` required to accept a trial step. | Tuning acceptance strictness. | N/A | Higher values reject more trials, increasing shrink/retry iterations. | Compared against `rho` computed from the merit function. |
| `power_flow.trust_region.shrink_factor` | Float64 | `0.5` | `(0, 1)` | Radius multiplier applied on a rejected trial. | Tuning how aggressively the radius shrinks after a bad step. | N/A | Smaller values shrink faster, reaching `min_radius`/collapse sooner. | Applied repeatedly within one Newton iteration until accepted or collapsed. |
| `power_flow.trust_region.expand_factor` | Float64 | `2.0` | `> 1` | Radius multiplier applied on a strongly successful boundary-hitting step. | Letting the radius grow once the model proves reliable. | N/A | Larger values reach `max_radius` faster. | Only applied when the step also hit the radius boundary. |
| `power_flow.trust_region.expand_threshold` | Float64 | `0.75` | `(0, 1)` | `rho` threshold above which an accepted, boundary-hitting step triggers expansion. | Tuning how "good" a step must be before trusting a larger radius. | N/A | Higher values expand less often, growing the radius more conservatively. | Must be `>= eta_accept` in practice for a coherent expand/accept ordering (not separately validated). |

### Diagnostics and interpretation

When a run has an output directory in its performance profile and `trust_region.enabled = true`, the solver writes `trust_region.log` with one line per Newton iteration: `radius_before`, `rho`, `tested_radii` (the shrink sequence tried this iteration), `rejected_steps`, `accepted`, `radius_after`, and `collapsed`. The solver status is extended with `trust_region_enabled`, `tr_step_count`, `tr_rejected_steps`, `tr_min_radius`, `tr_max_radius` (the observed radius range across the run, not the configured bounds), `tr_final_radius`, and `tr_collapsed`.

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

## Q-limit options and guard

| YAML path | Type | Default | Allowed values | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `power_flow.qlimits.enabled` | Bool | `true` | `true`, `false` | Master switch for Q-limit enforcement. | Realistic PV/PQ behavior needed. | Pure unconstrained PF tests. | Can increase switching iterations. | Gates all `qlimits.*`. |
| `power_flow.qlimits.enforcement_mode` | Symbol/String | `active_set` | `active_set`, `classic_simultaneous`, `classic_one_at_a_time` | Selects active-set or classical Q-limit switching. | Compare dynamic switching with classical outer-loop enforcement. | Treating legacy `matpower_*` aliases as preferred public values. | Classical modes may require repeated solves. | See [Q-limit Switching Strategy](q_limit_switching_strategy.md). |
| `power_flow.qlimits.start_iter` | Int | `3` | integer | First iteration index for Q limits. | Delay switching noise early. | Very late switching on hard cases. | Affects convergence speed. | `start_mode`, `auto_q_delta_pu`. |
| `power_flow.qlimits.start_mode` | Symbol/String | `iteration_or_auto` | `iteration`, `auto`, `iteration_or_auto` | Activation policy. | Mixed robustness/perf runs. | Mismatched with expected policy. | Small control logic overhead. | `auto_q_delta_pu`. |
| `power_flow.qlimits.auto_q_delta_pu` | Float64 | `1e-4` | nonnegative real | Auto activation threshold. | Fine-tuning switch timing. | Extreme values. | Low. | `start_mode=auto` or `iteration_or_auto`. |
| `power_flow.qlimits.hysteresis_pu` | Float64 | `0.01` | nonnegative real | Hysteresis margin near Q limits. | Reduce switch chattering. | Too large if strict tracking needed. | Can reduce oscillatory iterations. | `cooldown_iters`, guard modes. |
| `power_flow.qlimits.cooldown_iters` | Int | `1` | nonnegative integer | Cooldown iterations after switching. | Reduce repeated toggling. | Too long cooldown on tight limits. | Affects convergence pace. | Hysteresis and freeze behavior. |
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

Use `refresh_sparlectra_config_file(path; write=false)` to check an existing user YAML against the current template without modifying it. The refresh helper preserves existing user-provided values, adds missing keys from `src/configuration.yaml.example`, reports duplicate keys, and can normalize known deprecated aliases when `normalize_deprecated=true`.

Writing is explicit: pass `write=true` to update the file. By default a timestamped `.bak-YYYYmmdd-HHMMSS` backup is created before writing, and duplicate keys prevent automatic writes. This refresh mechanism is a maintenance tool only; normal configuration loading still accepts supported Q-limit enforcement legacy aliases, and Sparlectra never silently rewrites user YAML during startup.

State Estimation
=================

Weighted least squares (WLS) state estimation reconstructs the network state
from redundant, noisy measurements, with observability analysis and bad-data
diagnostics. The estimator extensions (FACTS view, links, islands, tap
estimation, topology validation) are on
[State Estimation Extensions](state_estimation_extensions.md), the file
format, the measurement generator and the chained power flow on
[State Estimation Measurements](state_estimation_measurements.md).

## Theory (compact)

Classical WLS formulation:

* State vector: `x = [θ(non-slack); Vm(all buses)]`
* Measurement model: `z = h(x) + e`
* Objective: `J(x) = (z - h(x))' * W * (z - h(x))`

`z` is the measurement vector, `h(x)` the nonlinear prediction from the
network model, `W = diag(1/σ²)` the inverse-variance weighting. The
algorithm linearizes `h(x)` and iterates Newton-style
until the update norm or the residual criteria meet the tolerance. With PMU
voltage-angle measurements the state vector gains the PMU reference-angle
offset α (see below): `x = [θ(non-slack); Vm(all buses); α]`.

## Why the FD measurement Jacobian works

The measurement Jacobian `H = ∂h/∂x` is built by finite differences, column
by column:

```math
\frac{\partial h}{\partial x_k}(x)
\approx
\frac{h(x + \varepsilon e_k) - h(x)}{\varepsilon}.
```

WLS only needs the local first-order sensitivity
($h(x + \delta x) \approx h(x) + H(x)\,\delta x$), so the estimation model is
the same as with an analytic Jacobian; only the derivative is evaluated
numerically.

## Measurement model

Supported measurement types:

* `VmMeas` (bus voltage magnitude, p.u.)
* `VaMeas` (bus voltage angle in degrees, PMU synchrophasor)
* `PinjMeas`, `QinjMeas` (bus injections, MW/MVar)
* `PflowMeas`, `QflowMeas` (branch flows with direction, MW/MVar)
* `ImagMeas` (branch current magnitude in ampere, with direction)
* `IaMeas` (PMU current angle), `ShuntQMeas` (shunt bay reactive power)

Synthetic-data workflow for studies and tests: solve a power flow, create
measurements with `generateMeasurementsFromPF` (sigmas from
`measurementStdDevs`, optional Gaussian noise). Real operation uses field
measurements directly.

### PMU voltage-angle measurements (`VaMeas`)

PMUs measure the bus voltage phasor, magnitude and absolute phase angle,
time-synchronized via GPS (IEEE C37.118). A PMU reports angles against its
GPS clock, so all PMU angles share one unknown rotation against the slack
reference. Sparlectra estimates that rotation as an additional state, the
reference-angle offset `α` (the slack-bus angle in the PMU time base), and
reports it with the result.

**Use**

| | |
|---|---|
| Add | `addVaMeasurement!(net; busName, value, sigma)` (degrees); `addPmuPhasorMeasurement!(net; busName, vm_pu, va_deg)` writes the Vm+Va pair |
| Generate | `generateMeasurementsFromPF(includeVa = true)`; default sigma `measurementStdDevs(va = 0.02)` |
| Config key | `state_estimation.pmu_ref_offset` (keyword `pmuRefOffset` on `runse!`, `validate_measurements`, `evaluate_global_observability`, ...): `:auto` (default) adds the offset state as soon as active `VaMeas` rows exist, `:off` takes PMU angles as slack-referenced |
| Result field | `SEResult.vaRefOffsetDeg` |

**Notes**

- With `:auto`, a PMU time base that coincides with the slack reference gives
  `α` of about 0, so the mode is safe to leave on.
- With `:off` and a wrong time base, `J` inflates by orders of magnitude and
  every PMU angle row shows up in the bad-data ranking (see
  `state_estimation_pmu_angles.jl`).
- The state count becomes `n = (2·nbus − 1) + 1 = 2·nbus`; the first PMU
  angle measurement pins `α` and is critical on its own. Redundancy for `α`
  (and offset-robust bad-data detection on PMU angles) needs at least two
  PMU angle measurements; a `VaMeas` at the slack bus reads `α` directly
  (`θ_slack = 0`).
- `VaMeas` values and sigmas are in degrees (the `_va_deg` convention).
  Angle residuals wrap at the plus/minus 180 degree seam, for `VaMeas` and
  `IaMeas` alike, so two angles across the seam give their physically small
  residual, not one near 360 degrees.
- PMU magnitudes are ordinary `VmMeas` rows with a small sigma (e.g. 0.002
  p.u.).

!!! details "Why it is built this way"
    **Accuracy.** Typical angle standard deviations are 0.01 to 0.05°
    (total vector error below 1 %), against percent-level SCADA power
    measurements. With `w = 1/σ²` PMU angle rows get a very large weight;
    they are modeled as ordinary, well-weighted measurements, not as hard
    constraints.

    **Reference.** The classical estimator fixes `θ_slack = 0` and
    estimates relative angles. The offset state absorbs the PMU rotation:

    ```math
    z_{Va,i} = \theta_i + \alpha + e_i,
    ```

    so the Jacobian gains one column (`∂z_{Va,i}/∂α = 1` on PMU angle rows,
    zero elsewhere) and the network angles `θ_i` stay slack-referenced.

### Branch current-magnitude measurements (`ImagMeas`)

`ImagMeas` models SCADA branch current magnitudes in ampere at one branch
end, predicted as `I = 1000 · |S| / (√3 · Vn · |V|)` from the branch flow at
that end. Current magnitudes are auxiliary: they supplement power
measurements and never replace them. What they buy is bad-data
localizability.

**Use**

| | |
|---|---|
| Add | `addImagMeasurement!(net; ...)` with `direction = :from`/`:to`; the bus-referenced variant `addImagMeasurement!(busName = ...)` (no direction) measures the shunt-bay current |
| Generate | `generateMeasurementsFromPF(includeImag = true)`; default sigma `measurementStdDevs(imag = 10.0)` (a class 1 instrument on a range of a few hundred ampere) |
| Config key | `state_estimation.imag_activation_iteration` (default 2, keyword `imagActivationIteration`): the first WLS iteration that takes current rows |

**Notes**

- No observability: `evaluate_global_observability` and
  `evaluate_local_observability` drop `ImagMeas` rows before building the
  Jacobian, so the system must be observable from power and voltage
  measurements alone.
- Iteration gate: `ImagMeas` rows stay out of the WLS update while the
  iteration count is below `imag_activation_iteration`, so the first
  linearization from a flat start runs on power and voltage rows only.
- Value gate: a current with `value < 3 sigma` is excluded for the whole run
  (one info line per exclusion).
- Bus-referenced currents feed the shunt back-calculation (case B below).

!!! details "Why it is built this way"
    Three rules keep currents auxiliary: observability from power and
    voltage alone, the iteration gate, and the value gate. Near zero current
    the derivative of `|I|` is discontinuous, which is what the value gate
    avoids. What currents buy is bad-data localizability: they raise the
    residual sensitivities `w_ii` of the power measurements (mean `w_ii` of
    the active-power rows from 0.58 to 0.81) and shrink the residual
    correlations.

### PMU current-phasor angles (`IaMeas`)

`IaMeas` is the angle of the branch-end (or shunt-bay) current in degrees,
referenced to the PMU time base like `VaMeas`: the offset α applies to both,
and current phasors alone also activate the offset state. The prediction
shares the complex end current with the magnitude (`I = conj(S/V)` scaled to
ampere); residuals wrap at the plus/minus 180 degree seam.

**Use**

| | |
|---|---|
| Add | `addIaMeasurement!` (mirrors `addImagMeasurement!`); `addCurrentPhasorMeasurement!` writes the `ImagMeas` + `IaMeas` pair in one call |
| Config key | `state_estimation.ia_current_floor_A` (default 10 A): magnitude floor when no paired `ImagMeas` exists |
| Diagnostics | gated rows are listed under `gating_notes` (`:ia_below_current_floor`) |

**Notes**

- An `IaMeas` row is active only while the predicted magnitude passes
  `3 * sigma_I_ref`, the sigma of the paired `ImagMeas` at the same end,
  else `state_estimation.ia_current_floor_A`.
- A row whose paired magnitude falls to the `ImagMeas` value gate leaves
  with it.
- Current angles share the `ImagMeas` iteration gate and are excluded from
  observability.

!!! details "Why it is built this way"
    A current angle is meaningless near zero current, so gated rows are out
    of the solve, of `J` and of the redundancy. A linear PMU-only estimation
    carrying observability through current phasors is future work.

### Shunt estimation and back-calculation (`ShuntQMeas`)

The operating point of switchable shunts (reactors, capacitor banks) is
often wrong or stale in the model. With a direct bay measurement the
susceptance becomes an estimator state (case A); without one, a bay current
plus a measured voltage yields a derived reactive-power pseudo-measurement
(case B).

**Use**

| | |
|---|---|
| Case A, a direct bay measurement exists | `addShuntQMeasurement!(net; busName, value, sigma)` (MVar, sign convention of the solved `q_shunt`) or a bus-referenced `ImagMeas`; `setShuntEstimation!(net; busName = ...)` turns the susceptance into an estimator state |
| Case B, no direct Q measurement | `deriveShuntPseudoMeasurements!(net; sigmaFloor)` before `runse!` emits a `ShuntQMeas` pseudo-measurement (id prefix `SHDERIV`) per shunt bay with a bus-referenced `ImagMeas` and a measured voltage |
| Config key | `state_estimation.update_shunts` (keyword `updateShunts`, default true), see [Write-back policy](@ref) |
| Result field | `SEResult.shuntEstimates`: `B_model`, `B_est`, `delta` and `frozen` per released shunt |
| Artifact | `shunt_estimates.csv` |

**Notes**

- A released shunt without an active direct bay measurement (`ShuntQMeas` or
  bus-referenced `ImagMeas`) stays frozen at the model value with a warning,
  and so does a `B` column no measurement pins.
- Shunts under the voltage-dependent injection mode (`bus_shunt_model`) are
  rejected.
- Case B needs the measured voltage (it enters `Q = -B V^2` squared); a
  shunt without one is skipped with a warning.
- Case B takes the unmeasured bay `P` as 0 and doubles the propagated sigma;
  for filter circuits with a substantial active component use a direct
  `ShuntQMeas` (case A).
- Each released shunt consumes one degree of freedom (`ν = m - n` counts the
  `B` columns).

!!! details "Why it is built this way"
    **Case A.** The state vector becomes
    `x = [θ(non-slack); Vm(all); B_1..B_k; α]`, initialized from the model
    value (never flat). The released admittances leave the Ybus diagonal
    once per run and the predictions re-add the injection analytically,
    `S_sh = |V_i|^2 conj(G + jB) baseMVA`, so the Ybus stays constant across
    the FD perturbations.

    **Case B.** The derivation takes `Q_sh` from `S = sqrt(3) U I`, the sign
    from the model susceptance, and the sigma from first-order error
    propagation (`sigmaFloor` sets a minimum). No bay active-power
    measurement type exists, hence the `P = 0` assumption. `SHDERIV` rows
    are protected from elimination like the `ZI` rows, since eliminating
    one would mask its source measurements.

### Passive / transit buses

Buses without load, generation, or shunt get no hard equality-constraint
block in the WLS solver. Model them with zero-injection pseudo-measurements
`Pinj = 0` and `Qinj = 0` of very small variance: `findPassiveBuses(net)`
detects the buses, `addZeroInjectionMeasurements!(meas; net, sigma=...)`
appends the rows. Otherwise a passive node in a sparse set leaves the
estimator merely critical or weakly redundant.

## Observability

[Observability](observability.md) covers the global and local checks, the
quality labels, critical measurements and their zero residual, `wii` and the
0.3 guideline, the correlation bound, zero-injection handling, the
diagnostics table and the H-matrix demo. Run
`evaluate_global_observability(net; ...)` before estimating (structural
islands, then FD-aware numeric rank; verdict `:observable` / `:critical` /
`:not_observable`), and `evaluate_local_observability(net, cols; ...)` for
placement studies on chosen state columns.

## Integration with the Net workflow

SE runs on the same `Net` as the power flow:

1. Build or import the `Net`
2. Build measurements (SCADA/PMU/custom)
3. Optional for synthetic studies: `runpf!` + `generateMeasurementsFromPF`
4. Check observability (global/local)
5. Run the estimator (`runse!`)
6. Optionally write the estimate back (`state_estimation.update_net = true`)

Power flow computes states from setpoints, SE from measured values;
redundancy enables bad-data detection from residual statistics.

**Diagnostics workflow**

| Call | Returns |
|---|---|
| `validate_measurements(net, measurements)` | runs SE once (thresholds from `state_estimation.k_eliminate` and friends) and returns objective statistics (`value`, `dof`, `zscore`, `within_3sigma`), the largest normalized residual, the full ranking by `\|normalized_residual\|` with `wii` and a `localizable` flag per row, the suspicious list, and optional residual-correlation columns |
| `runse_diagnostics(net, measurements)` | adds sequential bad-data elimination: while the 3σ band test fails and a suspicious measurement exists, the top suspect is deactivated and the diagnostics rerun, up to `state_estimation.max_eliminations` (default 3) times; one trace row per elimination (id, normalized residual before, `wii`, objective before/after) and a `stop_reason`: `:consistent`, `:no_suspicious_left`, `:no_localizable_suspect`, `:max_eliminations`, or `:not_converged` |
| `summarize_se_diagnostics(diag)` | compact summary (`global_consistency`, `reason`, suspicious count) |
| `print_se_diagnostics(diag; io=stdout, topN=10, format=:plain\|:markdown)` | statistics, ranking (with `wii`) and elimination trace |

Zero-injection rows (id prefix `ZI`) and near-exact rows (`sigma <= 1e-6`)
are never eliminated; they encode network structure. `global_consistency`
is `true` when the estimator converged and the objective passes the
Wilson-Hilferty band test; otherwise `false` with a `reason`.

### Write-back policy

Estimation never mutates the model unless one of the write-back keys is
set, and then only from a converged run: `state_estimation.update_net`
(keyword `updateNet`) writes the estimated voltages,
`state_estimation.update_taps` (keyword `updateTaps`, off by default)
the fixed tap positions, `state_estimation.update_shunts` (keyword
`updateShunts`, on by default: a released shunt is released to be
corrected) the estimated susceptances. A frozen regulator or shunt keeps
its exact model value inside a write-back; without the flag the model value
is bitwise untouched. The snapshot power flow of
[PF from the estimate](state_estimation_measurements.md#PF-from-the-estimate)
follows the same rule: its temporary prosumers are removed after the run.

### The Wilson-Hilferty band test

The objective $J = \sum_i (z_i - h_i(\hat{x}))^2 / \sigma_i^2$ says whether
measurements, sigmas and model fit together as a whole. When every sigma
honestly describes its measurement error, `J` is a χ² variable with
`ν = m - n` degrees of freedom (`m` measurements, `n` states). Its
expected value is `E[J] = ν` and its standard deviation `√(2ν)`, hence
the rule of thumb `J ≈ ν` (the run message writes `dof` for `ν`). The run
message leads with `J/dof` and the band verdict; on failure the report
names the direction.

**Verdicts**

| Verdict | Meaning |
|---|---|
| passed (`\|z_WH\| ≤ 3`) | measurements, sigmas and model agree |
| `:high` | `J` too large: the residuals exceed what the sigmas allow, bad data or a model error (a wrong parameter, tap or switch state), the classical alarm |
| `:low` | `J` implausibly small: the sigmas overstate the errors; the normal signature of noise-free synthetic data, not an alarm |
| `:no_redundancy` | `ν = 0`, the test is skipped (the residuals are structurally zero) |
| `small_redundancy` flag | `ν < 30`, informational: reduced power, relative spread `√(2/ν)` |

**Notes**

- Read `J/ν`, not `J`: a bare `J` grows with the number of measurements
  (`J = 104` at `dof = 95` is a ratio of 1.1 and normal).
- The classical three-sigma score `z = (J - ν)/√(2ν)` is still reported; the Wilson-Hilferty verdict
  decides.
- `runse_diagnostics` eliminates only while the band test fails `:high`.
- On multi-island nets every island carries its own verdict next to the
  summed one (see [Island-wise estimation](@ref)).
- When an exhausted elimination still leaves `:high` with the suspects
  clustered at one station, the diagnosis shifts to a suspected topology
  error (see [Topology validation](@ref)).

!!! details "Why it is built this way"
    **The familiar rule.** `J` has the expected value `ν` and the standard
    deviation `σ_J = √(2ν)`, so the textbook check is the three-sigma band
    `ν - 3σ_J ≤ J ≤ ν + 3σ_J`, that is `|z| ≤ 3` for the z-score
    `z = (J - ν) / √(2ν)`. Sparlectra still reports this `z`.

    **What is wrong with it.** The three-sigma band assumes `J` is normally
    distributed. `J` is χ²-distributed, and that distribution is skewed to
    the right: the tail above `ν` is longer than the tail below it. A
    symmetric band therefore lets too many large `J` pass and flags too many
    small ones; at small `ν` its lower bound `ν - 3√(2ν)` even drops below
    zero, where `J` cannot go (`ν = 4`: lower bound `-4.5`).

    **The correction.** Wilson and Hilferty showed that the cube root
    `(J/ν)^{1/3}` is very nearly normal, with mean `1 - 2/(9ν)` and variance
    `2/(9ν)`, usable from `ν ≈ 3` upward. The band test applies the same
    three-sigma rule to that cube root instead of to `J`:

    ```math
    z_{WH} = \frac{(J/\nu)^{1/3} - \left(1 - \tfrac{2}{9\nu}\right)}{\sqrt{2/(9\nu)}}, \qquad |z_{WH}| \le 3
    ```

    Translated back to `J` the band is asymmetric, wider above `ν` and
    narrower below it. For `ν = 95` it runs from about `J = 59` to
    `J = 142`, where the symmetric rule gives 54 to 136; the two disagree
    only near the edges, and there the Wilson-Hilferty verdict is the right
    one.

### [Bad-data thresholds and robust modes](@id se-bad-data)

Both limits work on the normalized residual $r_{n,i} = r_i / \sqrt{\Omega_{ii}}$,
computed with the original sigmas; only the staged knees
`robust_k1`/`robust_k2` stay on the raw ratio $|r_i|/\sigma_i$. Down-weighting
and elimination read the same quantity, so their limits are directly
comparable. Rule of thumb: suppression for online smoothing, elimination for
identification.

**Modes**

| Mode | What it does | Config key | When to use |
|---|---|---|---|
| Elimination | a row is suspicious from `rn >= k_eliminate` on; the elimination workflow removes it, bounded by the budget `max_eliminations`, only while the band test reports `:high` | `state_estimation.k_eliminate` (default 3.0), `state_estimation.max_eliminations` (default 3) | identification of gross errors |
| `off` | plain WLS | `state_estimation.robust_mode = off` | clean data |
| `staged` | two-stage R modification, active from `state_estimation.robust_start_iteration` (default 3) on, with knees `robust_k1`/`robust_k2` (defaults 3/6) on the raw ratio `t = \|z - h(x)\|/σ`: `t ≤ k1` leaves the weight unchanged, `k1 < t ≤ k2` widens tangentially (`σ_mod = σ (2t/k1 - 1)`, continuous at `t = k1`), `t > k2` suppresses the gradient contribution (`σ_mod = \|z - h\|/k1`); stages are recomputed every iteration, so a measurement can recover | `state_estimation.robust_mode = staged`; the Boolean `state_estimation.robust` is an alias for `staged` and applies while `robust_mode` is `off` | the classic behavior |
| `replacement` | pins every row whose normalized residual reaches `k_suppress` to the fixed `suppression_sigma` (in the measurement's unit) during the solve; the row stays in the system | `state_estimation.robust_mode = replacement`, `state_estimation.k_suppress` (default 4.0), `state_estimation.suppression_sigma` (default 2000) | online smoothing |

**One scale for both decisions**

| Limit | Default | What happens at it |
|---|---|---|
| `k_suppress` | 4.0 | the row is down-weighted: it is solved with the replacement sigma and loses its influence, but stays in the system and in every statistic |
| `k_eliminate` | 3.0 | the row is removed from the estimate (elimination workflow, bounded by `max_eliminations`) |

**Result fields**

| | |
|---|---|
| Affected rows | `SEResult.robustRows` (stage 1/2 for staged, stage 3 for replacement, with `t` and `σ_mod/σ`) |
| Active objective | `SEResult.activeObjective` (`J_active`; also the Web UI summary, `run.log`, metadata keys `se_objective_active` / `se_dof_active` / `se_suppressed_rows`) |

**Notes**

- Only localizable rows are touched: a row whose residual sensitivity `wii`
  does not exceed 0.3 (the `localizable` flag of the report) is neither
  down-weighted nor eliminated, however far past its limit `rn` sits.
  Non-localizable suspects are reported and passed over, the trace counts
  them per round (`skipped_unlocalizable`), and if only non-localizable
  suspects remain the workflow stops with `:no_localizable_suspect`: the
  data cannot tell which row is wrong, a measurement gap rather than bad
  data and a case for more telemetry, not for the robust solve.
- Virtual measurements (`σ ≤ 1e-6`, including the `ZI` rows) are exempt in
  every mode.
- The replacement threshold is judged on the converged state, never on
  flat-start transients: the estimator solves plain, builds the suppression
  set from the converged residuals, re-solves with those weights frozen, and
  repeats until the set is stable.
- A row whose `\Omega_{ii}` sits at the numerical floor is not localizable,
  gets `rn = 0`, and is never suppressed.
- Elimination and suppression combine freely; the Web UI warns when
  `k_suppress < k_eliminate`, because suppressed rows then rarely reach the
  elimination.

!!! details "Why it is built this way"
    **One scale.** `suppression_sigma` is not a third limit: it is the sigma
    a down-weighted row is solved with. The normalized residual is the right
    scale: at `wii = 0.3` the residual spread is
    `\sqrt{\Omega_{ii}} = 0.55\,\sigma`, and for a nearly critical row
    `\Omega_{ii}` approaches zero while the raw ratio `|r_i|/\sigma_i` stays
    small, so a raw-ratio limit is weakest where a gross error does the most
    damage. In `staged` mode the knees and the elimination limit thus read
    on different scales, and the form says so.

    **Only localizable rows.** The residual of a non-localizable row is
    structurally small, a large `rn` says little about the row, and
    intervening pushes the error onto the neighbour that shared its
    redundancy (measured: the partner then reached `rn` 58 and was
    eliminated although healthy).

    **Solve weights only.** The modification changes only the solve
    weights: `J`, normalized residuals, `wii`, `K` and the suspicion ranking
    run on the original sigmas, otherwise the suppression would hide the
    very measurement it suppresses. Judging the replacement threshold on
    flat-start transients would suppress healthy tight rows, hence the
    converged-state rule.

    **J versus J_active.** With replacement suppression `J_active` is the
    objective over the rows the estimator trusted, with its own dof. The
    honest `J` says whether the data carry a problem, `J_active` whether
    the estimate is healthy despite it; the band verdict stays on the
    honest `J`, so suppression cannot silence the alarm.

### [Localizability: the residual sensitivity `wii`](@id se-localizability)

`wii = Ω_ii · w_i` is the share of measurement `i`'s own error that reaches
its residual (near 0: nearly critical, the error hides in the state). The
report flags `localizable = wii > wiiThreshold` with the literature
threshold 0.3; below it the largest-normalized-residual logic cannot be
trusted. Background in [Observability](observability.md).

### [Sequential elimination budget](@id se-max-eliminations)

Upper bound for the sequential bad-data elimination: while the chi-square
band test fails and a suspicious measurement (normalized residual at or
above 3) is localizable (`wii` > 0.3), the worst row is deactivated and the
estimation reruns, up to this many times. `0` disables elimination
(diagnostics only).

Protected rows are never eliminated: zero-injection pseudo-measurements
(`ZI`), derived shunt rows (`SHDERIV`), and `LINKAGG` cluster aggregates
encode model knowledge, not telemetry.

The elimination trace (which row, normalized residual before, objective
drop) lands in `se_diagnostics.md`. Configuration key
`state_estimation.max_eliminations` (default 3); the Web UI field
**Sequential elimination budget** sets it per run.

### [Residual correlations (optional K-matrix report)](@id se-k-report)

With the correlation report on, each ranking row carries the maximum
|correlation coefficient| of its normalized residual over all partners,
`K = D^{-1/2} Ω D^{-1/2}`. Above `1/√2 ≈ 0.707` two measurements form a
simple-redundant group: a gross error in one is statistically
indistinguishable from an error in the other, and the report marks the row.
Reporting only, no automatic action.

**Use**

| | |
|---|---|
| Config key | `state_estimation.report_residual_correlation` (default true, keyword `reportResidualCorrelation`) |
| Config key | `state_estimation.takahashi_min_states` (default 200): the state count above which `Ω_ii` comes from the sparse Takahashi path |
| Result field | `omega_path` (the path taken), `state_variances = diag(G⁻¹)` (confidence intervals) |

**Notes**

- With active `ImagMeas` rows the sensitivities `W` and correlations `K`
  are load-flow dependent: valid for the estimated operating point, not
  network constants.
- The correlation columns are skipped with a warning above 20000
  measurement rows (see the size limits below).
- A Takahashi guard failure falls back to the dense path with a warning.

!!! details "Why it is built this way"
    `Ω_ii` (and with it `wii` and the normalized residuals) is computed on
    one of two paths. Above `takahashi_min_states` states (measured
    dense/sparse crossover near 130 to 150) the sparse gain matrix is
    factorized with UMFPACK and the Takahashi selected inverse delivers
    exactly the `G⁻¹` entries the diagonal needs. Below the threshold, for
    K-matrix requests (the correlations need the full `Ω`, m×m dense by
    definition), and on any Takahashi guard failure the dense `pinv` path
    runs instead. Both paths report `state_variances` and are bounded (size
    limits below); the measurement Jacobian itself is assembled sparse
    throughout.

### What runs at which size

Warm timings on one development machine, as orders of magnitude, for Vm
plus P/Q injection sets (the `runse!` column includes the topology
pre-check, the diagnostics column times the residual diagnostics on an
assembled `H`, and the first row pays the compilation for the ones after
it, the same effect as in [Tests](tests.md)):

| Network | m × n | `runse!` | residual diagnostics | notes |
|---|---|---|---|---|
| sp_case60 (60 buses) | 326 × 117 | 0.8 s | dense path (below the 200-state Takahashi threshold) | every diagnostic available |
| sp_case188 (188 buses) | 562 × 373 | 0.1 s | Takahashi, 0.2 s | every diagnostic available |
| case1354pegase | 4062 × 2707 | 0.1 s | Takahashi, 0.1 s | every diagnostic available; criticality from `diag(Omega)` (no size budget) |
| case13659pegase | 40977 × 27317 | 2.0 s | Takahashi, 1.9 s | criticality from `diag(Omega)`; K-matrix report skipped above 20000 rows |

Three announced bounds:

* single-row criticality and the dark-state SVD are skipped with a warning
  above an m·n budget of 300000;
* the dense `pinv` diagnostics fallback stops at 2000 states; beyond it the
  estimator errors with the reason, including why the Takahashi path
  refused;
* the `state_estimation.report_residual_correlation` columns are skipped
  with a warning above 20000 measurement rows.

Large assemblies use finite-difference column coloring: states that share
no measurement row are perturbed together (12 colors on the 60- and 188-bus
demo cases, 28 on case1354, 85 on case13659, so 85 prediction sweeps
instead of 27317), bit-identical to the per-column assembly. Estimator and
diagnostics share one state definition: the diagnostics Jacobian carries the
released shunts' `B` states at their estimated values (frozen shunts at the
model value) and the state count `n` in `ν = m - n` includes them.

## Minimal example

```julia
using Sparlectra
using Random

result = run_sparlectra(casefile = "case9.m")
net = result.net

std = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7)
setMeasurementsFromPF!(
    net;
    includeVm = true,
    includePinj = true,
    includeQinj = true,
    includePflow = true,
    includeQflow = true,
    noise = true,
    stddev = std,
    rng = MersenneTwister(42),
)

gobs = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net)
end
println("Global observability quality: ", gobs.quality)

se = with_state_estimation_config(max_iter = 12, tol = 1e-6, flatstart = true, jac_eps = 1e-6, update_net = true) do
  runse!(net)
end
println("Converged: ", se.converged, ", iterations: ", se.iterations)
```

## Example without PF pre-step (measurement-driven)

```julia
using Sparlectra

net = Net(name = "se_measurement_driven", baseMVA = 100.0)
addBus!(net = net, busName = "B1", vn_kV = 110.0)
addBus!(net = net, busName = "B2", vn_kV = 110.0)
addBus!(net = net, busName = "B3", vn_kV = 110.0)
addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B3", toBus = "B1", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)

ok, msg = validate!(net = net)
ok || error("Validation failed: \$msg")

empty!(net.measurements)
append!(net.measurements, Measurement[
    Measurement(typ = VmMeas, value = 1.01, sigma = 0.002, busIdx = 1, id = "VM_B1"),
    Measurement(typ = VmMeas, value = 0.99, sigma = 0.004, busIdx = 2, id = "VM_B2"),
    Measurement(typ = PinjMeas, value = -25.0, sigma = 1.0, busIdx = 2, id = "PINJ_B2"),
    Measurement(typ = QinjMeas, value = -8.0, sigma = 1.0, busIdx = 2, id = "QINJ_B2"),
    Measurement(typ = PflowMeas, value = 24.0, sigma = 0.8, branchIdx = 1, direction = :from, id = "PF_12"),
    Measurement(typ = PflowMeas, value = 23.5, sigma = 0.8, branchIdx = 1, direction = :from, id = "PF_12_REDUNDANT"),
])

obs = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net)
end
println("Observable quality: ", obs.quality)

se = with_state_estimation_config(max_iter = 12, tol = 1e-6, flatstart = true, jac_eps = 1e-6, update_net = true) do
  runse!(net)
end
println("Converged: ", se.converged)
```

## Adding measurements with helper functions

Helper functions resolve bus names and branch references, so no
`Measurement(...)` has to be constructed by hand:

```julia
using Sparlectra

net = Net(name = "se_helpers", baseMVA = 100.0)
addBus!(net = net, busName = "B1", vn_kV = 110.0)
addBus!(net = net, busName = "B2", vn_kV = 110.0)
addBus!(net = net, busName = "B3", vn_kV = 110.0)
addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B3", toBus = "B1", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)

empty!(net.measurements)
addVmMeasurement!(net; busName = "B1", value = 1.01, sigma = 0.002)
addPinjMeasurement!(net; busName = "B2", value = -25.0, sigma = 1.0)
addQinjMeasurement!(net; busName = "B2", value = -8.0, sigma = 1.0)
addPflowMeasurement!(net; fromBus = "B1", toBus = "B2", value = 24.0, sigma = 0.8, direction = :from)
addQflowMeasurement!(net; branchNr = 1, value = 6.5, sigma = 0.8, direction = :to)

obs = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
  evaluate_global_observability(net)
end
println("Observable quality: ", obs.quality)
```

## PMU example

```julia
using Sparlectra

# ... build net, e.g. from a MATPOWER case ...

# PMU phasor at bus "B2": accurate magnitude + GPS-referenced angle (degrees).
addPmuPhasorMeasurement!(net; busName = "B2", vm_pu = 1.012, va_deg = -3.72)

# Equivalent, with the single-component helpers:
#   addVmMeasurement!(net; busName = "B2", value = 1.012, sigma = 0.002)
#   addVaMeasurement!(net; busName = "B2", value = -3.72, sigma = 0.02)

se = runse!(net)                 # pmu_ref_offset = :auto (default)
println("PMU reference offset α: ", se.vaRefOffsetDeg, " deg")
```

Literature on PMU-based state estimation:

* A. Abur, A. Gómez Expósito: *Power System State Estimation: Theory and
  Implementation* (hybrid SCADA/PMU WLS formulation).
* A. G. Phadke, J. S. Thorp: *Synchronized Phasor Measurements and Their
  Applications* (PMU measurement principle, IEEE C37.118 accuracy classes).
* NASPI TR-006: *Phase Angle Calculations: Considerations and Use Cases*
  (reference-angle handling across PMU installations).

A [Sparlectra Case Format](scf.md) case carries its own measurement set, so
an estimation on it needs no separate CSV; see
[Case files carry their own measurements](@ref se-case-file-measurements).

## Further examples and workshop material

* Tutorial with a 7-bus setup: the state-estimation chapter of the [advanced workshop tour](generated/workshop_tour_advanced.md) and the dedicated [state-estimation notebook](generated/workshop_state_estimation.md)
* Bad-data localization, sequential elimination, robust estimation, and shunt-parameter estimation, hands-on: the [SE diagnostics notebook](generated/workshop_se_diagnostics.md)
* Tap estimation (ratio, phase, cascade, the non-identifiable case), current measurements, PMU phasors, the topology advisory, island-wise estimation, $J$ versus $J_{active}$, and a closing decision tree for "$J$ is too large, what now": the [taps, phasors and topology notebook](generated/workshop_se_taps.md)
* Detailed WLS reporting example script: `examples/state_estimation/state_estimation_wls.jl`
* PMU angle measurements and the reference-offset state α: `examples/state_estimation/state_estimation_pmu_angles.jl`
* Observability-focused scenario script: `examples/state_estimation/state_estimation_observability.jl`
* Passive-bus ZIB comparison example: `examples/state_estimation/state_estimation_passive_bus_zib_comparison.jl`
* Matrix-based observability/redundancy demo: `examples/state_estimation/h_matrix_observability_demo.jl`

## Further reading

* [State Estimation Extensions](state_estimation_extensions.md): FACTS in
  the estimator, links, island-wise estimation, transformer tap estimation,
  topology validation.
* [State Estimation Measurements](state_estimation_measurements.md):
  measurement file format, the measurement generator, case files with
  measurements, power flow from the estimate.

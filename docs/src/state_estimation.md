State Estimation
=================

This page summarizes the state-estimation (SE) functionality in Sparlectra and
shows how it connects to regular network studies.

> **Release status:** State Estimation is a regular Sparlectra feature. The
> implementation provides a practical WLS workflow for studies, examples,
> and applications.

## Theory (compact)

Sparlectra currently provides a classical weighted least-squares (WLS)
formulation:

* State vector: `x = [θ(non-slack); Vm(all buses)]`
* Measurement model: `z = h(x) + e`
* Objective: `J(x) = (z - h(x))' * W * (z - h(x))`

Where:

* `z` is the measurement vector,
* `h(x)` is the nonlinear prediction of each measurement from the network model,
* `W = diag(1/σ²)` is the inverse-variance weighting matrix.

The algorithm linearizes `h(x)` and iterates Newton-style until the update norm
or residual criteria satisfy tolerance.

When PMU voltage-angle measurements are present, the state vector is
augmented by one additional state (the PMU reference-angle offset α, see
below):

* State vector: `x = [θ(non-slack); Vm(all buses); α]`

## Why the FD measurement Jacobian works

The internal helper `_measurement_jacobian_fd` approximates the Jacobian of the
measurement model `h(x)` by finite differences. The key idea is the same as in
power-flow FD Newton methods, but now the nonlinear map is the SE prediction
function rather than the PF mismatch function.

If `h(x)` is differentiable, then for a small perturbation `δx`:

```math
h(x + \delta x) \approx h(x) + H(x)\,\delta x,
```

where `H(x) = \partial h / \partial x` is the measurement Jacobian. Each
Jacobian column can therefore be approximated numerically via

```math
\frac{\partial h}{\partial x_k}(x)
\approx
\frac{h(x + \varepsilon e_k) - h(x)}{\varepsilon}.
```

This works because WLS only needs the **local first-order sensitivity** of the
measurements with respect to the state in order to build the linearized normal
equations. The underlying estimation model does not change; only the derivative
evaluation is numerical instead of analytic.

Conceptually:

* PF FD Jacobian approximates derivatives of the residual map `F(x)`.
* SE FD Jacobian approximates derivatives of the measurement map `h(x)`.

In both cases, the finite-difference step is justified by the same first-order
Taylor approximation.

## Measurement model

The current implementation supports these measurement types:

* `VmMeas` (bus voltage magnitude, p.u.)
* `VaMeas` (bus voltage angle in degrees, PMU synchrophasor — see below)
* `PinjMeas`, `QinjMeas` (bus injections, MW/MVar)
* `PflowMeas`, `QflowMeas` (branch flows with direction, MW/MVar)
* `ImagMeas` (branch current magnitude in ampere, with direction, see below)

### PMU voltage-angle measurements (`VaMeas`)

Phasor measurement units (PMUs) measure the bus voltage phasor — magnitude
**and** absolute phase angle — time-synchronized via GPS (IEEE C37.118).
Two properties distinguish PMU angles from classical SCADA measurements:

1. **They are very accurate.** Typical angle standard deviations are
   0.01–0.05° (total vector error < 1 %), versus percent-level accuracy for
   SCADA power measurements. In the WLS weighting `w = 1/σ²` this gives PMU
   angle rows a very large weight; PMU measurements are therefore modeled as
   ordinary — merely well-weighted — measurements, not as hard constraints.
2. **They are referenced to a common time base, not to the slack bus.**
   The classical estimator fixes `θ_slack = 0` and estimates *relative*
   angles. A PMU reports angles relative to its GPS clock, so all PMU angles
   share one unknown rotation against the estimator's slack reference.

Sparlectra resolves the reference mismatch with an **additional state
variable**: the reference-angle offset `α` (the slack-bus angle expressed in
the PMU time base). Each PMU angle measurement is predicted as

```math
z_{Va,i} = \theta_i + \alpha + e_i,
```

so the measurement Jacobian gains one extra column (`∂z_{Va,i}/∂α = 1` for
every PMU angle row, zero elsewhere). The network angles `θ_i` stay
slack-referenced; `α` absorbs the common rotation and is reported in
`SEResult.vaRefOffsetDeg`.

Behavior is controlled by `state_estimation.pmu_ref_offset`
(keyword `pmuRefOffset` on `runse!`, `validate_measurements`,
`evaluate_global_observability`, ...):

* `:auto` (default): the offset state is added as soon as active `VaMeas`
  measurements exist. If the PMU time base happens to coincide with the
  slack reference, `α` is simply estimated as ≈ 0 — the mode is safe to
  leave on.
* `:off`: no offset state; PMU angles are assumed to be slack-referenced
  already. If that assumption is wrong, the common rotation cannot be
  absorbed: the objective `J` inflates by orders of magnitude and all PMU
  angle rows show up as suspicious in the bad-data ranking (see the
  `state_estimation_pmu_angles.jl` example for a demonstration).

Observability accounting with the offset state:

* The state count becomes `n = (2·nbus − 1) + 1 = 2·nbus`.
* The first PMU angle measurement pins `α` and is therefore *critical* on
  its own; redundancy for `α` (and offset-robust bad-data detection on PMU
  angles) requires at least two PMU angle measurements.
* A `VaMeas` at the slack bus directly reads `α` (since `θ_slack = 0`).

Practical notes:

* `VaMeas` values and sigmas are in **degrees**, matching the network-facing
  `_va_deg` convention; `measurementStdDevs(va = 0.02)` provides the default.
* Angle residuals wrap at the plus/minus 180 degree seam, for voltage
  angles (`VaMeas`) and current angles (`IaMeas`) alike, so a pair of
  angles on opposite sides of the seam produces the small residual it
  physically is instead of one near 360 degrees.
* PMU magnitude measurements need no special treatment — model them as
  `VmMeas` with a small sigma (e.g. 0.002 p.u.).

### Branch current-magnitude measurements (`ImagMeas`)

SCADA systems routinely provide branch current magnitudes. `ImagMeas` models
them in ampere at one branch end (`direction = :from`/`:to`), predicted as
`I = 1000 · |S| / (√3 · Vn · |V|)` from the branch flow at that end. Add them
with `addImagMeasurement!`, generate them synthetically with
`generateMeasurementsFromPF(includeImag = true)`; the default sigma comes
from `measurementStdDevs(imag = 10.0)` (a class 1 instrument on a measuring
range of a few hundred ampere, the typical HV line current transformer).

Current magnitudes are **auxiliary** measurements: they
supplement power measurements but must never replace them. Three rules
enforce this:

1. **No observability.** `evaluate_global_observability` and
   `evaluate_local_observability` exclude `ImagMeas` rows before building the
   Jacobian: the system must be observable from power and voltage
   measurements alone, and the verdict is invariant under adding currents.
2. **Iteration gate.** The estimator keeps `ImagMeas` rows out of the WLS
   update step while the iteration count is below
   `state_estimation.imag_activation_iteration` (default 2, keyword
   `imagActivationIteration`), so the first linearization from a flat start
   runs on power/voltage measurements only.
3. **Value gate.** A current measurement with `value < 3 sigma` is excluded
   for the whole run (one info line per exclusion): near zero current the
   derivative of `|I|` is discontinuous and the row would destabilize the
   iteration.

What currents do buy is **bad-data localizability**: additional current
measurements raise the residual sensitivities `w_ii` of the power
measurements (mean `w_ii` of the active-power rows from
0.58 up to 0.81) and shrink the residual correlations, so gross errors
separate more clearly. They also feed the shunt
back-calculation (case B below); the bus-referenced shunt-bay variant of
`ImagMeas` (`addImagMeasurement!(busName = ...)`, no direction) measures the
current drawn by the shunt at that bus.

### PMU current-phasor angles (`IaMeas`)

`IaMeas` completes the PMU current phasor: the ANGLE of the branch-end (or
shunt-bay) current in degrees, referenced to the PMU time base exactly like
`VaMeas` (the common offset alpha applies to both; current phasors alone
also activate the offset state). `addIaMeasurement!` mirrors
`addImagMeasurement!`, and `addCurrentPhasorMeasurement!` writes the
`ImagMeas` + `IaMeas` pair in one call. The prediction shares the complex
end current with the magnitude (`I = conj(S/V)` scaled to ampere), and
angle residuals wrap exactly at the plus/minus 180 degree seam.

A current angle is meaningless near zero current, so an `IaMeas` row is
active only while the predicted current magnitude passes
`3 * sigma_I_ref`, where `sigma_I_ref` is the sigma of the paired
`ImagMeas` at the same end, else the config floor
`state_estimation.ia_current_floor_A` (default 10 A). A row whose paired
magnitude falls to the `ImagMeas` value gate leaves with it. Gated rows
are excluded from the solve, from `J`, and from the redundancy, and the
diagnostics list them under `gating_notes`
(`:ia_below_current_floor`). Current angles share the `ImagMeas`
iteration gate and are excluded from observability in this stage; a
linear PMU-only estimation that would let current phasors carry
observability is future work.

### Shunt estimation and back-calculation (`ShuntQMeas`)

The actual operating point of switchable shunts (compensation reactors,
capacitor banks) is often wrong or stale in the model. Sparlectra covers the
two concept cases:

**Case A, a direct bay measurement exists.** `ShuntQMeas` is the reactive
power of the shunt bay in MVar (sign convention of the solved `q_shunt`
results, helper `addShuntQMeasurement!`). Releasing a shunt with
`setShuntEstimation!(net; busName = ...)` turns its susceptance into an
estimator state: the state vector becomes
`x = [θ(non-slack); Vm(all); B_1..B_k; α]`, initialized from the model value
(never flat). Internally the released shunts' model admittances leave the
Ybus diagonal once per run and the predictions re-add the injection
analytically, `S_sh = |V_i|^2 conj(G + jB) baseMVA`, with the state `B`; the
Ybus stays constant across the FD perturbations. Guards keep the mechanism
honest: a released shunt without an active direct bay measurement
(`ShuntQMeas` or bus-referenced `ImagMeas`) stays frozen at the model value
with a warning, and so does a `B` column no measurement pins (checked on
the built Jacobian). Shunts under the voltage-dependent injection mode
(`bus_shunt_model`) are rejected in this phase. Results arrive in
`SEResult.shuntEstimates` (`B_model`, `B_est`, `delta`, `frozen` per
released shunt); the model itself is only overwritten behind the explicit
`updateShunts = true` keyword (config `state_estimation.update_shunts`).
Each released shunt consumes one degree of freedom (`ν = m - n` counts the
`B` columns as states).

**Case B, no direct Q measurement.** `deriveShuntPseudoMeasurements!(net)`
is a preprocessing step before `runse!`: for every shunt bay with a current
measurement (bus-referenced `ImagMeas`) and a measured voltage it derives
`Q_sh` from `S = sqrt(3) U I`, takes the sign from the model susceptance
(reactor versus capacitor), and emits a `ShuntQMeas` pseudo-measurement (id
prefix `SHDERIV`) whose sigma comes from first-order error propagation
(`sigmaFloor` optionally sets a minimum). The measured voltage is a hard
prerequisite: it enters `Q = -B V^2` squared, so nominal voltage is not an
acceptable substitute and the shunt is skipped with a warning instead.

Deliberate 0.10.0 restriction: no bay active-power measurement type exists,
so the derivation always takes the bay `P` as 0 (the concept foresaw an
optional bay P) and doubles the propagated sigma to cover the neglected
active part. For reactors and capacitor banks `P` is small against `Q` and
the error is secondary; for filter circuits with a substantial active
component the derived value is biased, use a direct `ShuntQMeas` (case A)
there instead.

`SHDERIV` rows are protected from bad-data elimination like the `ZI`
pseudo-measurements: eliminating one would silently mask its source
measurements.

### Passive / transit buses

For buses without load, generation, or shunt contribution, Sparlectra does
**not** currently introduce a separate hard equality-constraint block in the
WLS solver. Instead, the recommended modeling approach is to add
zero-injection pseudo-measurements

* `Pinj = 0`
* `Qinj = 0`

for those buses. In other words, the physical equality constraint is embedded
through very small-variance measurements in the standard WLS formulation.

Helper functions:

* `findPassiveBuses(net)` detects passive / transit buses from the bus power
  aggregates.
* `addZeroInjectionMeasurements!(meas; net, sigma=...)` appends the matching
  zero-injection pseudo-measurements automatically.

This is especially useful in sparse measurement scenarios, where a passive node
may otherwise leave the estimator merely critical or weakly redundant.

At the moment, this is the supported way to model ZIB behavior in Sparlectra.
There is not yet a separate hard-constraint solver block for zero-injection
buses.

Typical synthetic-data workflow (for studies/tests):

1. Solve a power flow to get a physically consistent reference state.
2. Create synthetic measurements using `generateMeasurementsFromPF`.
3. Configure standard deviations via `measurementStdDevs`.
4. Optionally add Gaussian noise.

In real operation, SE uses field measurements directly and does not require a
preceding power-flow run to create data.

## Observability

Observability has its own page: [Observability](observability.md)
covers what the global and local checks answer, the quality labels,
critical measurements and their zero residual, `wii` and the 0.3
guideline, the correlation bound, zero-injection handling, how to read
the diagnostics table, and the H-matrix demo as the didactic entry.
Short form for this page's flow: run
`evaluate_global_observability(net; ...)` before estimating (the
two-stage check, structural islands then FD-aware numeric rank, yields
`:observable` / `:critical` / `:not_observable`), and use
`evaluate_local_observability(net, cols; ...)` for placement studies
on chosen state columns.

## Integration with the Net workflow

SE is designed to run on the same `Net` data model used for power flow:

1. Build/import `Net`
2. Build measurements (SCADA/PMU/custom)
3. Optional for synthetic studies: run `runpf!` + `generateMeasurementsFromPF`
4. Check observability (global/local)
5. Run estimator (`runse!`)
6. Optionally write estimates back into the network (`updateNet = true`)

Conceptually, SE is the measurement-driven counterpart of power flow:

* Power flow computes states from setpoints.
* SE computes states from measured values.
* Measurement redundancy improves robustness and enables bad-data detection
  using residual statistics.

Sparlectra exposes a public diagnostics workflow for bad-data and statistical
consistency checks:

* `validate_measurements(net, measurements; normalizedThreshold=3.0, wiiThreshold=0.3, ...)`
  runs SE once and returns:
  * objective statistics (`value`, `dof`, `zscore`, `within_3sigma`)
  * largest normalized residual
  * full measurement ranking by `|normalized_residual|`, each row carrying
    the residual sensitivity `wii` and a `localizable` flag (see below)
  * suspicious measurement list based on `normalizedThreshold`
  * optional residual-correlation columns (see below)
* `runse_diagnostics(net, measurements; max_eliminations=3, ...)` extends
  this with **sequential bad-data elimination**: while the χ²-like 3σ band
  test fails and a suspicious measurement exists, the top suspect is
  deactivated and the diagnostics rerun, up to `max_eliminations` times.
  The result carries one trace row per elimination (id, normalized residual
  before, `wii`, objective before/after) and a `stop_reason`:
  `:consistent`, `:no_suspicious_left`, `:no_localizable_suspect`,
  `:max_eliminations`, or `:not_converged`.
  Zero-injection pseudo-measurements (id prefix `ZI`) and near-exact
  measurements (`sigma <= 1e-6`) are never eliminated; they encode network
  structure, not telemetry. The old keyword `deactivate_and_rerun = true`
  remains as an alias for `max_eliminations = 1`.
* `summarize_se_diagnostics(diag)` creates a compact interpretation summary
  (`global_consistency`, `reason`, suspicious count).
* `print_se_diagnostics(diag; io=stdout, topN=10, format=:plain|:markdown)`
  pretty-prints the statistics, the ranking (with the `wii` column), and the
  elimination trace for reports.

Meaning of `global_consistency`:

* `true`: estimator converged and the objective passes the Wilson-Hilferty
  band test (globally plausible data/model fit).
* `false`: non-convergence, or the band test failed (see the reason).

### The Wilson-Hilferty band test

**What the test asks.** The WLS objective
$J = \sum_i (z_i - h_i(\hat{x}))^2 / \sigma_i^2$ is the one number that
says whether measurements, sigmas, and model fit together AS A WHOLE.
When every declared sigma honestly describes its measurement error, `J` is
a χ² variable with `ν = m - n` degrees of freedom: its expected value IS
`ν`, with a spread of about `√(2ν)`. That gives the rule of thumb the
result pages print: expected `J ≈ dof`. A `J` far above says the
residuals are larger than the sigmas allow (bad data or a model error); a
`J` far below says the residuals are smaller than the sigmas promise (the
sigmas are overstated, nothing is wrong with the network).

**Why a transformation.** The χ² distribution is right-skewed, strongly
so at small `ν`: the symmetric rule `|J - ν| ≤ 3√(2ν)` misjudges both
tails there (its lower bound can even go negative, so an implausibly small
`J` always "passed"). The Wilson-Hilferty transformation fixes this with
one cube root: `(J/ν)^{1/3}` is very nearly NORMAL with mean `1 - 2/(9ν)`
and variance `2/(9ν)`, usable from `ν ≈ 3` upward. The decision statistic
is its z-score

```math
z_{WH} = \frac{(J/\nu)^{1/3} - \left(1 - \tfrac{2}{9\nu}\right)}{\sqrt{2/(9\nu)}}
```

and the band test passes when `|z_WH| ≤ 3`. In `J` terms the acceptance
interval is asymmetric, as the skewed distribution demands: for `ν = 95`
it runs from about `J = 59` to `J = 142`, shifted upward against the
symmetric rule's 54 to 136, which mislabels both tails. The legacy
symmetric z-score is still reported for continuity but no longer decides.

**Read $J/\nu$, not $J$.** A bare `J` grows with the number of measurements,
so the same healthy set reports `J = 49` on a small case and `J = 211` on a
larger one. `J = 104` at `dof = 95` is a ratio of 1.1 and perfectly normal.
The run message therefore leads with `J/dof` and the band verdict, and keeps
the raw `J` and `dof` behind it.

On failure the report distinguishes the two directions, because they mean
opposite things:

* `:high`: `J` too large. Bad data or a model error (a wrong parameter, a
  wrong tap, a wrong switch state), the classical alarm.
* `:low`: `J` implausibly small. The measurement sigmas overstate the
  errors; this is the normal signature of noise-free synthetic data and
  NOT an alarm. The old symmetric test silently passed this case, hiding
  overstated sigmas.

Special cases: `ν = 0` skips the test with reason `:no_redundancy` (the
residuals are structurally zero, bad data is invisible); `ν < 30` sets the
informational `small_redundancy` flag (the test stays valid, its power is
reduced, relative spread `√(2/ν)`).

**Consequences in the workflow.** The band verdict is the driver of the
bad-data loop: `runse_diagnostics` eliminates suspicious measurements ONLY
WHILE the band test fails `:high` and stops as soon as the objective is
plausible again (a passing `J` means there is nothing left to justify
removing a measurement). On multi-island nets every island carries its own
verdict next to the summed one, so a single bad island cannot hide inside
an otherwise healthy sum (see [Island-wise estimation](@ref)). And when an
exhausted elimination still leaves `:high` with the suspects clustered at
one station, the diagnosis shifts from bad data to a suspected topology
error (see [Topology validation](@ref)).

### Bad-data thresholds and robust modes

Both limits work on the NORMALIZED residual
`rn = r_i / \sqrt{\Omega_{ii}}`, computed with the ORIGINAL sigmas, and
are exposed in the Web UI's estimator options. Only the staged knees
`robust_k1`/`robust_k2` stay on the raw ratio `|r|/\sigma` (see below).

**Elimination limit.** `state_estimation.k_eliminate` (keyword
`normalizedThreshold`, default 3.0) bounds the candidacy of the sequential
elimination: a row is suspicious from `rn >= k_eliminate` on. It works
together with the elimination budget (`max_eliminations`) and only while
the band test reports `:high`.

**Robust mode.** `state_estimation.robust_mode` selects the solve-weight
family: `off` (plain WLS), `staged`, or `replacement`. The legacy Boolean
`state_estimation.robust` remains an alias for `staged` and applies while
`robust_mode` is `off`.

`staged` is the two-stage R modification (active from
`state_estimation.robust_start_iteration`, default 3) with configurable
knees `robust_k1`/`robust_k2` (defaults 3/6 reproduce the classic
behavior bitwise). These knees read the RAW ratio `t = |z - h(x)|/σ`, not
the normalized residual of the two limits above: `t ≤ k1` leaves the weight
unchanged, `k1 < t ≤ k2`
widens tangentially (`σ_mod = σ (2t/k1 - 1)`, continuous at `t = k1`),
`t > k2` suppresses the gradient contribution (`σ_mod = |z - h|/k1`).
Stages are recomputed every iteration, so a measurement can recover.

`replacement` pins every row whose normalized residual
`rn = r_i/\sqrt{\Omega_{ii}}` reaches `state_estimation.k_suppress`
(default 4.0) to the fixed `state_estimation.suppression_sigma` (default
2000, in the measurement's unit) during the solve; the row stays in the
system. The threshold is judged on the CONVERGED state, never on
flat-start transients: the estimator solves plain first, builds the
suppression set from the converged residuals, re-solves with those weights
frozen, and repeats until the set is stable. (A per-iteration threshold
would suppress healthy tight rows, such as passive-node balances, during
the transient and let the estimate drift.)

**One scale for both decisions.** Down-weighting and elimination read the
same quantity, the normalized residual, so their two limits are directly
comparable:

| Limit | Default | What happens at it |
|---|---|---|
| `k_suppress` | 4.0 | the row is DOWN-WEIGHTED: it is solved with the replacement sigma and loses its influence, but stays in the system and in every statistic |
| `k_eliminate` | 3.0 | the row is REMOVED from the estimate (elimination workflow, bounded by `max_eliminations`) |

`suppression_sigma` is not a third limit: it is the sigma a down-weighted
row is solved with, in the unit of the measurement.

Until 0.10.0 the suppression judged the RAW ratio `|r_i|/\sigma_i`
instead. That measure differs from `rn` exactly where it matters: with
`wii = \Omega_{ii} w_i`, a row at the localizability guideline
`wii = 0.3` has `\sqrt{\Omega_{ii}} = 0.55\,\sigma`, and for a nearly
critical row `\Omega_{ii}` approaches zero while the raw ratio stays
small. The suppression was therefore weakest precisely where a gross error
does the most damage. The `wii` guard carries over unchanged: a row whose
`\Omega_{ii}` sits at the numerical floor is not localizable, gets
`rn = 0`, and is never suppressed.

**Only localizable rows are touched.** One precondition for BOTH decisions:
a row whose residual sensitivity `wii` does not exceed 0.3 (the guideline
the report prints as `localizable`) is neither down-weighted nor
eliminated, however far past its limit the `rn` sits. At `wii \approx 0` the row is nearly critical: its
residual is structurally small, so a large `rn` says little about the row
itself, and intervening there does not remove bad data, it removes the
little information the row still carries and pushes the error onto the
neighbour that shared its redundancy. Measured on the warm-up PST case,
suppressing such a row (`wii` 0.05, raw residual 1.4 sigma, `rn` 5.8) drove
its partner to `rn` 58 and had that HEALTHY row eliminated.

The same bound guards the ELIMINATION, and there it matters more: removing
a row is the sharper intervention, so it must not have the weaker
precondition. A critical measurement has a structurally zero residual, the
estimate follows it completely and nothing contradicts it, so a gross error
there is invisible and the row can never earn its removal; removing it
anyway only costs observability. Non-localizable suspects are therefore
reported and passed over, the trace counts them per round
(`skipped_unlocalizable`), and if suspects remain while none of them is
localizable the workflow stops with `:no_localizable_suspect` instead of
removing the next best row and calling the case solved. That stop reason
says the data cannot tell WHICH row is wrong, which is a measurement gap
around those rows, not bad data.

The staged knees `robust_k1`/`robust_k2` deliberately keep their classic
`|r|/\sigma` scale. In `staged` mode the two limits on screen therefore
read on DIFFERENT scales: the knees on `|r|/\sigma`, the elimination
limit on `rn`. They are not directly comparable there, and the form says
so. In `replacement` mode both share the `rn` scale.

In every mode virtual measurements (`σ ≤ 1e-6`, including the `ZI`
pseudo-measurements) are exempt, and the modification changes ONLY the
solve weights. Every statistic (the objective `J`, normalized residuals,
`wii`, `K`, the suspicion ranking) runs on the original sigmas; otherwise
the suppression would hide the very measurement it suppresses.
`SEResult.robustRows` lists the affected rows (stage 1/2 for staged,
stage 3 for replacement, with `t` and `σ_mod/σ`). Elimination and
suppression combine freely; the Web UI warns (without rejecting) when
`k_suppress < k_eliminate`, because suppressed rows then rarely reach the
elimination. Rule of thumb: suppression for online smoothing, elimination
for identification.

When replacement suppression removed rows from the state, the result
additionally reports **J_active** (`SEResult.activeObjective`, the Web UI
summary, `run.log`, and the metadata keys `se_objective_active` /
`se_dof_active` / `se_suppressed_rows`): the objective over the rows the
estimator actually trusted, with its own dof. The pair reads as "the
honest J says the DATA carry a problem, J_active says whether the
ESTIMATE is healthy despite it". The 3σ band verdict deliberately stays
on the honest `J`: judged on J_active, suppression could silence the very
alarm it should raise.

### Localizability: the residual sensitivity `wii`

`wii = Ω_ii · w_i` is the diagonal of the residual sensitivity matrix (the
share of measurement `i`'s own error that reaches its residual). `wii ≈ 1`
means an error there shows up almost fully in the residual (well
localizable); `wii ≈ 0` marks a nearly critical measurement whose error
hides in the state estimate. The report flags `localizable = wii >
wiiThreshold` with the literature threshold 0.3: below
it, the largest-normalized-residual logic cannot be trusted to point at the
truly faulty meter.

### Residual correlations (optional K-matrix report)

With `state_estimation.report_residual_correlation = true` (keyword
`reportResidualCorrelation`), each ranking row additionally carries the
maximum |correlation coefficient| of its normalized residual over all
partners, `K = D^{-1/2} Ω D^{-1/2}`. Above `1/√2 ≈ 0.707` two measurements
form a simple-redundant group: a gross error in one is statistically
indistinguishable from an error in the other, and the
report marks the row. This is reporting only, no automatic action.

Two caveats:

* With active `ImagMeas` rows, the sensitivities `W` and correlations `K`
  become load-flow dependent: they are valid for the
  estimated operating point, not network constants.
* `Ω_ii` (and with it `wii` and the normalized residuals) is computed on one
  of two paths, reported as `omega_path` in the diagnostics: above
  `state_estimation.takahashi_min_states` states (default 200, the measured
  dense/sparse crossover sits near 130 to 150) the sparse gain matrix is
  factorized with UMFPACK and the shared **Takahashi selected inverse**
  delivers exactly the `G⁻¹` entries the diagonal needs (every required
  index pair is a structural nonzero of `G` and therefore inside the factor
  pattern; measured 4.6x over the dense `pinv` at 1800 states). Below the
  threshold, for K-matrix requests (the correlations need the full `Ω`),
  and on any Takahashi guard failure the original dense `pinv` path runs
  instead, with a warning on a fallback, never a wrong number. Both paths
  also report `state_variances = diag(G⁻¹)` (confidence intervals).
* Since the sparse core (task_se_sparse, 0.10.0) the dense conveniences
  are bounded instead of unbounded: the `pinv` fallback runs only up to
  2000 states, and the K-matrix report is refused by name above 20000
  measurement rows (`Ω` is m×m dense by definition). Beyond those bounds
  the estimator states the reason, including why Takahashi refused,
  instead of allocating; the measurement Jacobian itself is assembled
  sparse throughout, which is what makes large networks reachable at
  all.

### What runs at which size

Measured behavior of the 0.10.0 estimator (Vm plus P/Q injection sets),
taken on the maintainer machine with the shipped demo cases and the two
pegase cases; the previous release refused the largest case up front with
a computed 22 GB dense-Jacobian requirement. Read the numbers as orders of
magnitude on comparable hardware, not as a benchmark. The
`runse!` column is the estimation including its topology pre-check; the
diagnostics column times `_residual_diagnostics` on an assembled `H`,
which the pre-check does not enter:

| Network | m × n | `runse!` | residual diagnostics | notes |
|---|---|---|---|---|
| sp_case60 (60 buses) | 326 × 117 | 0.8 s | dense path (below the 200-state Takahashi threshold) | every diagnostic available |
| sp_case188 (188 buses) | 562 × 373 | 0.1 s | Takahashi, 0.2 s | every diagnostic available |
| case1354pegase | 4062 × 2707 | 0.1 s | Takahashi, 0.1 s | criticality classification skipped (`criticality_skipped = true`, warning names the m·n budget) |
| case13659pegase | 40977 × 27317 | 2.0 s | Takahashi, 1.9 s | as above; K-matrix report refused by name above 20000 rows |

Read these times WARM. The first row pays the compilation for the ones
after it, which is why the 60-bus case looks eight times slower than the
188-bus case; the same effect is documented for the test groups in
[Tests](tests.md).

Three bounded behaviors, each announced rather than silent: single-row
criticality and the dark-state SVD are skipped with a warning above an
m·n budget of 300000 (the 188-bus demo case already spends seconds
there), the dense `pinv` diagnostics fallback stops at 2000 states and
beyond it the estimator errors with the reason including why the
Takahashi path refused, and `state_estimation.report_residual_correlation`
is refused by name above 20000 measurement rows. Large assemblies use
finite-difference column coloring derived from the measurement
topology: states that share no measurement row are perturbed together
(12 colors on the 60- and 188-bus demo cases, 28 on case1354, 85 on
case13659), which is what turns the 27317-state Jacobian from 27317
prediction sweeps into 85 while staying bit-identical to the per-column
assembly.

Estimator and diagnostics share the same state definition since 0.10.0:
the diagnostics Jacobian carries the released shunts' `B` states at their
estimated values (frozen shunts stay at the model value, exactly as in the
estimator), and ν counts the `B` states.

### Known limitations (0.10.0)

* **Case B assumes an unmeasured bay `P` of zero** (see the deliberate
  restriction above): biased for filter circuits with a substantial active
  component.
* **Down-weighting only reaches localizable rows.** A row with
  `wii <= 0.3` is never suppressed, even past the limit: its residual is
  structurally small, so its normalized residual says little about the row
  itself, and suppressing it pushes the error onto the neighbour that
  shared its redundancy. Bad data on a nearly critical measurement stays a
  case for more telemetry, not for the robust solve.
* The three bounded behaviors above (dense `pinv` up to 2000 states,
  K-matrix report refused above 20000 rows, criticality and dark-state SVD
  skipped above the m·n budget) are limits as well; they are described
  where they occur rather than repeated here.

This workflow can be used both for automated checks (NamedTuple result
inspection) and for human-readable diagnostics output.

## FACTS in state estimation (`se_view`)

State estimation is a snapshot: `runse!` never invokes the outer control
loop, every controller is frozen at its current operating point. The concept
transformation table applies:

| Device (PF view) | SE view | estimable |
|---|---|---|
| OLTC (in-phase regulator) | transformer with fixed tap | future work |
| PST (phase-shifting regulator) | transformer with fixed tap | future work |
| Schraegregler | transformer with fixed taps | future work |
| SVC / STATCOM | shunt with Q injection | `B` per case A above |
| Compensation reactor / capacitor bank | shunt | `B` per case A/B above |
| Series compensation | fixed branch impedance | no |

`se_view(net)` reports this frozen view without mutating the net: the
registered controllers, the shunt-estimation releases and their
classification, the link clusters the estimator will fuse, and every
measurement the WLS will exclude or aggregate. `print_se_view(view)` renders
it (`format = :markdown` available).

## Links in state estimation

Bus links (impedance-less couplers, `addLink!`) are not part of the Ybus, so
the estimator handles them exactly like the power flow: **SE runs on the
contracted net.** Every closed-link cluster is fused onto its representative
bus (dense renumbering via the island subnet machinery); open links remain
real separations. With `updateNet = true` all cluster members receive the
representative's estimated voltage, mirroring the PF result sync.

Measurement rules on a fused cluster:

* `VmMeas`/`VaMeas` on a member remap exactly (voltage equality).
* `PinjMeas`/`QinjMeas` aggregate only as a whole: when every cluster member
  carries an active injection measurement of the kind, one `LINKAGG` row
  with the summed value and root-sum-square sigma replaces them;
  otherwise all member injections of that kind are excluded with a warning
  (a partial sum would bias the fused balance). `LINKAGG` rows appear in the
  diagnostics with `measurement_index = 0` and are never eliminated.
* `ShuntQMeas` and bus-referenced `ImagMeas` follow their device to the
  representative; branch measurements keep their branch.
* A branch that collapses inside a cluster (shorted by the closed coupler)
  loses its flow measurements with a warning.

**Link flow measurements are allocation inputs, never WLS rows.** A P/Q flow
measurement on a link (`addPflowMeasurement!(net; linkNr = ...)`, positive
from `link.fromBus` to `link.toBus`) says nothing about the fused system
state; it constrains only the flow split. `calcLinkFlowsSE!(net)` extends
the KCL allocation to a weighted least-squares split per link component:
nodal balance rows keep weight 1, each link measurement adds a row with
weight `1/sigma^2`, and without measurements the result equals
`calcLinkFlowsKCL!` exactly. The returned rows carry `source = :kcl` or
`:measured_ls` plus the measurement residual per link. Run it after
`runse!(...; updateNet = true)`, which fills the branch flows and shunt
powers the allocation consumes.

## Island-wise estimation

A net can split into several synchronous AC islands: multi-area MATPOWER
cases, and virtually every CGMES delivery (stub islands from disconnected
stations or boundary remnants are routine). Mirroring the power flow's
island-wise solving, `runse!` partitions the measurement set onto the
islands and estimates **every island that carries measurements with its own
angle reference**; unmeasured islands are skipped and reported. The rules:

* **Closed links fuse islands.** Partitioning respects the link
  contraction: islands connected by a closed busbar coupler are one
  estimation group (the contraction handles the cluster as before); an
  OPEN link remains a real separation.
* **Reference per island.** An island with its own slack keeps it. An
  island without one gets the MATPOWER-like promotion the power flow uses
  (`detect_ac_islands`' chosen reference bus). An unpowered stub island
  without any candidate still gets an ANGLE PIN on its first bus: unlike
  the PF slack, the estimator's reference only fixes the angle, the bus
  magnitude stays an ordinary state. The pin is prosumer-backed so bus-type
  refreshes cannot wipe it.
* **Merged result.** On a multi-island net `SEResult.voltages` is indexed
  by the ORIGINAL bus numbering (unestimated islands keep their start
  values), `objectiveJ`/`dof` are the sums over the estimated islands
  (independent chi-squares add), and `SEResult.islands` carries one row per
  island: estimated/skipped, iterations, per-island `J`/`dof` and its own
  Wilson-Hilferty verdict (`band_reason`, `z_wh`), so a single bad island
  cannot hide inside the summed band test. The printed diagnostics and
  `se_diagnostics.md` show the same per-island table.
* **Observability and diagnostics partition the same way.**
  `evaluate_global_observability` judges each measured island on its own
  subnet (quality = worst measured island, note `:island_partition`,
  per-island results in `islands`); `validate_measurements` and the
  sequential elimination in `runse_diagnostics` work across islands with
  measurement indices in the caller's vector, unchanged.
* **Chain.** `runpf_from_se!` enables island-wise PF solving on
  multi-island nets automatically; the slack pickup sums over the island
  references.

## Transformer tap estimation

A transformer whose actual tap position disagrees with the model (a stale
SCADA position, a local operation) poisons the estimate around it: the band
test reports `:high` and the suspicious measurements cluster at the
transformer although every telemetry row is healthy. Releasing the tap as
an estimation state resolves this.

**Release.** `setTapEstimation!(net; trafo, mode = :ratio, alpha_deg,
enabled = true)` marks a transformer branch (index or component name). The
released tap becomes one or two additional WLS states with the cascade
model

```math
t(r_1, r_2) = \frac{t_{base}}{(1 + r_1)\,(1 + r_2\, e^{j\alpha})}
```

built from the same `calcSkewAngleTap`/`calcAdmittance` convention as the
importers: `mode = :ratio` releases the longitudinal regulator $r_1$,
`:pst` the skew regulator $r_2$ with its fixed nameplate direction
`alpha_deg`, `:both` releases the pair. The regulator states start from the
CURRENT branch position (a contradictory position for the declared mode is
an error, not silently absorbed). The stamped admittance terms leave the
Ybus once and every prediction re-adds the branch at the cascade position
of the current state, so at the initial position the released model
reproduces the stamped power flow to machine precision.

**Mandatory fixation run.** A tap changer sits on a mechanical grid, not on
a continuum. After the estimation converges, each released regulator is
rounded to its nearest mechanical step (regulator 1 on the fraction grid
`tap_step`, regulator 2 on the shift grid `phase_step_deg`; out-of-range
positions are flagged and clamped) and ONE more run is solved in which the
tap is no state variable any more. The reported
`SEResult.voltages`/`objectiveJ`/`dof` are those of the fixed final run.

**Write-back.** Estimation never silently overwrites model data: the FIXED
mechanical positions reach the branch fields only with
`state_estimation.update_taps = true` (keyword `updateTaps`), only from a
converged fixed run, and a frozen regulator keeps its exact model value
inside the written cascade. Without the flag the model tap is bitwise
untouched (same protection as `update_shunts`).

**Machine transformers.** A generator step-up transformer is the mirror
case: its machine terminal is AVR-held and not independently observed, so
releasing its tap would only absorb the voltage control. Instead,
`calcMachineTrafoTapFromSE(net; trafo, v_machine_pu, p_mw, q_mvar)`
back-calculates the position AFTER the estimation from the estimated
network-side voltage, the AVR setpoint, the dispatch P, and the MEASURED
machine reactive power (under AVR control Q is telemetry, never a schedule
value). Release and back-calculation are mutually exclusive per
transformer; the result reports the electrical and nearest mechanical step
plus the reactive-power residual of the fit, and writes nothing back.

**Off-grid note.** A chi-square band failure that only APPEARS through the
fixation is not bad data: when the continuous fit (taps still states)
passes or undershoots the band and the fixed run fails it `:high`,
`tapFixation.offgrid_residual` is set and the result summary says so. The
J jump then comes from rounding to the mechanical step: the true position
sits between steps, or the step table (`tap_step`, neutral position) is
wrong. Chasing measurements with the bad-data workflow is the wrong move
in that situation.

**Reporting.** `SEResult.tapEstimates` carries one row per released
transformer: branch index and name (mRID for CGMES cases; MATPOWER/DTF
cases are addressed by branch and bus numbers), the continuous electrical
step, the fixed mechanical step, and range flags. `SEResult.tapFixation`
reports `j_before`/`dof_before` (taps still states) versus
`j_after`/`dof_after` (taps fixed): a tiny `j_before` with a large
`j_after` means the true position sits between mechanical steps; both tiny
means the fixed step explains the measurements. On multi-island nets the
rows carry the island id and the fixation sums aggregate like the island
chi-squares.

Note that the fixation makes `dof` go UP, not down — that is correct and
worth spelling out because the intuition "fixing = less freedom = fewer
degrees of freedom" points the wrong way. `dof = m - n` counts the
REDUNDANCY of the residuals (and is the expected value of `J` for healthy
noise), not the model's freedom. A released tap is one extra estimated
STATE: it consumes one measurement's worth of information to determine
itself, so `dof_before = m - (voltage states + taps)`. The fixation nails
the tap to its mechanical step and removes it from the state vector, so
that measurement information becomes pure check information again:
`dof_after = m - voltage states = dof_before + released taps`.
Consistently, `J` moves the other way: the continuous run can adjust one
more direction (smaller `J`, smaller expected value), the fixed run cannot
(slightly larger `J`, larger expected value) — e.g. `J 58.5 (dof 41) ->
J 62.5 (dof 42)`, both inside their band. Only a fixation that jumps OUT
of the band is a finding, and that is exactly what the
`:offgrid_tap_residual` note flags.

**Fallback when nothing settles.** The guards below freeze single taps that
cannot be pinned. They cannot see the case where the set carries the
voltages fine but is collectively too thin for the released taps: the run
then does not converge at all. A non-convergence WITH released taps
therefore freezes every released tap back to its model position and repeats
the estimation once. The run log names it (`repeating WITHOUT tap
estimation`) and the result metadata carries
`se_tap_estimation_fallback = true`, because the tap positions in such a
result are model values, not estimates. Measured on case300 with 98
released taps: no result at all before, a converged estimate after.

**Release guards.** Tap states are deliberately outside the observability
count (like the gated current rows), so dedicated guards protect the
release instead. Estimating a tap needs redundancy around the transformer:
a meshed path or measurements on both sides. A BRIDGE transformer (its
removal cuts the net) whose cut-off side carries no active voltage
measurement makes the tap and the downstream voltage indistinguishable;
such a release is frozen at the current position and reported with reason
`radial_no_voltage_pin`. Any remaining tap state the active measurements
cannot pin passes the same per-column numerical test the released shunt
states use and freezes with reason `not_observable`; a `:both` release can
freeze partially (one regulator stays a state). A transformer whose
regulators are all frozen leaves the solve with its model stamp restored,
bitwise identical to no release, and still appears in `tapEstimates` with
its reason. Machine (generator step-up) transformers, recognized as a
generator bus hanging on the transformer alone, are skipped by the Web UI
mass release; releasing one explicitly via `setTapEstimation!` remains
possible and deliberate.

In the Web UI the "estimate taps" option of the estimator run releases
every other in-service transformer with a ratio tap changer (mode
`:ratio`); the result page shows the per-transformer table (frozen taps
carry their reason in the status column) and the fixation J drop, and
`se_tap_estimates.csv` lands in the run artifacts. PST releases (`:pst`,
`:both`) are an API-level choice.

## Topology validation

A wrong service state (a breaker recorded closed that is actually open, or
the reverse) poisons an estimation differently from bad telemetry: the
errors cluster around one station and no elimination can cure them.
Sparlectra validates the topology in three stages, ALL ADVISORY: findings
and recommendations only, never a mutation of a status, a measurement, or
the model. Out-of-service elements are exempt from the plausibility
checks; only the status-contradiction check looks at them, because a
measured flow over an open element IS the topology error.

**Stage 1, pre-checks.** `validate_topology(net, measurements)` runs pure
and linear before any estimation (automatic at the start of `runse!` with
`state_estimation.topology_precheck = true`, the default): a measured flow
over an OPEN element (`:open_element_with_flow`), a CLOSED branch reading
dead at both measured ends while its neighbourhood carries load
(`:closed_element_without_flow`, low severity), voltage measurements
disagreeing across a closed link (`:closed_link_voltage_mismatch`), and
the node balance at COMPLETELY measured nodes (`:kcl_violation`; shunt
buses use the measured voltage for the shunt term, partially measured
nodes are skipped, never guessed). Thresholds are sigma multiples
(`state_estimation.topology_open_flow_k` and friends). The check RESULT is
logged on every service run, the clean case included; findings land in
`SEResult.topologyFindings` and warn.

**Stage 2, classification.** When the sequential elimination EXHAUSTS its
budget while the band test still fails `:high` and the surviving suspects
cluster at one station (at least `state_estimation.topology_cluster_min`,
station = closed-link contraction cluster), `runse_diagnostics` reports
`:topology_error_suspected_at_station` instead of the
eliminations-exhausted interpretation: this fingerprint separates a wrong
service state from bad telemetry. A curable gross error never produces a
topology finding (the elimination succeeds first, test-enforced).
Stage-1 findings at the same station are noted as `precheck_agreement`;
with correlation columns enabled a fully correlated cluster carries
`correlated_group`.

**Stage 3, hypothesis test.** `test_topology_hypotheses(net, measurements)`
toggles each candidate's service state on a WORKING COPY, re-runs the
estimation, and reports J and the band z-score before versus after,
ranked. A hypothesis is supported when the `:high` failure disappears
under the toggle (inside the band, or below it, which a noise-free set
produces). Several supported hypotheses are flagged ambiguous instead of
hidden. Candidates come from the stage-2 stations and the stage-1 status
contradictions (`candidates = :auto`, capped by `max_candidates` with
per-candidate timing), or from an explicit list. The input net is bitwise
untouched, and NOTHING is ever switched automatically: the output is a
recommendation list for the operator. In the Web UI the test sits behind
an explicit button on the SE result page and writes
`topology_hypotheses.md`.

## Measurement files and the SE chain

**Measurement CSV v1.** `writeMeasurementsCSV(net; file)` and
`readMeasurementsCSV!(net; file, replace = true)` persist measurement sets:

    # sparlectra-measurements v1
    type,bus,from_bus,to_bus,branch_nr,link_nr,direction,value,sigma,active,id

Exactly one location group per row (`bus` | `from_bus`+`to_bus`+`branch_nr`
| `link_nr`), buses by name (matched whitespace-insensitively, CGMES names
can carry padding), decimal point, UTF-8. Generated sets carry an explicit
case binding IN the file (a `# case: <name>` comment): the association
survives renames and re-uploads, the Web UI stars cases with a bound set
and labels foreign sets, and the SE service refuses a set bound to a
different case up front instead of failing row by row. Sets without a
binding (older files) still run, with a log note. The Web UI additionally
offers download, re-upload, and an inline editor (small files, saved
atomically, first line must stay the version comment); bad-data findings
land machine-readable in `se_bad_data.csv` (suspicious and eliminated
measurements with their network locations, plus per row whether the
robust solve SUPPRESSED it with the fixed replacement sigma or
down-weighted it in the staged modification; suppressed rows below the
elimination threshold are listed too). `branch_nr` disambiguates
parallel branches between the same bus pair (routine in CGMES nets); the
reader verifies it against the named endpoints and rejects a mismatch.
Files with the first-cut header (no `branch_nr` column) are still read;
their branch rows resolve by bus pair and reject parallels (same rule as
`addPflowMeasurement!`). The import is atomic with line-precise errors
(`file:line: reason`); an unknown version comment or type name rejects the
whole file. Writing and re-reading the file reproduces every row exactly,
for every measurement type including link rows. On CGMES-sourced nets the generator references buses by their
preserved ENTSO-E mRID (UUID; `writeMeasurementsCSV(busReference = :mrid)`,
reader resolves names, component ids, and mRIDs alike); MATPOWER and DTF
nets have no mRIDs and keep bus numbers/names. A generated file carries a
structured `sparlectra-taps v1` comment table with the transformer tap
positions of the generation state (electrical, fixed, and transferred step
per transformer, deviation percentage), which the SE page renders and
offers for download together with the file.

**PF from the estimate.** `runpf_from_se!(net, maxIte, tol, verbose; mode)`
starts a power flow from the last SE result (requires
`runse!(...; updateNet = true)` or a state restored via `readSEStateCSV!`):

* `mode = :se_state` (config `profile_source = state_estimation`): start
  from the estimated voltages; the model injections stay authoritative, the
  measurement/model difference goes into the (possibly distributed) slack.
* `mode = :se_snapshot` (config `profile_source = se_snapshot`): also take
  the nodal balances from the estimation. The takeover works through
  temporary delta load prosumers (the solver assembles injections from the
  prosumers, not from node aggregates; the slack is exempt since its
  injection is the balancing variable), removed after the run; the
  persistent model is never mutated. With consistent PV setpoints the PF
  converges in 0 or 1 iterations and the returned `slack_pickup_mw/mvar`
  stays below tolerance. This guarantee holds for the FROZEN operating
  point: with active outer-loop controllers (taps, machine control,
  Q-limits) the snapshot only guarantees the mismatch-free start; the
  controls then iterate normally and legitimately move the operating point,
  which shows up as a small residual pickup.

`writeSEStateCSV`/`readSEStateCSV!` persist the estimated state
(`se_state.csv`), which is how the Web UI chain hands an SE run's result to
a later power-flow or N-1 run (the N-1 report metadata records the SE run
id). Use `:se_state` when the model must stay authoritative and `:se_snapshot`
when the estimated operating point itself is the study base.

**Web UI.** The state-estimation section of the Runs page
(`/powerflow#state-estimation`; the former `/stateestimation` URL
redirects there) uses the shared case selection, MATPOWER or CGMES, and
selects a measurement set (uploaded `.csv`
files are offered only when the content sniff finds the v1 version comment;
the set generated for the selected case is preselected, sets generated for
another case are labeled), runs observability (traffic light,
structural-island note), the WLS solve, and the diagnostics; artifacts
(`measurements.csv`, `se_diagnostics.md`, `se_view.md`, `se_state.csv`,
`shunt_estimates.csv`) land in the run history with kind `se`, and the
result page offers "run power flow from this estimate"; its summary shows
the chi-square plausibility at a glance: `J`, the expected value
(`E[J] = dof` for healthy noise), and whether `J` lies inside the
Wilson-Hilferty 3-sigma band (with the `:high`/`:low` reason when it does
not). The summary always describes the state AFTER the elimination
workflow: eliminated rows are out of `J` and `dof` (a found and removed
gross error therefore brings `J` back to about `dof`), while suppressed
rows stay in the honest `J` and additionally surface as `J_active`. A
"Reset saved settings" button under the generator deletes the per-case
settings sidecar, so all forms of the case fall back to their defaults. The demo generator solves the case once (island-wise where the
delivery has multiple islands) and writes `<case>.measurements.csv`; the
measurement sigmas are entered per quantity in PERCENT OF THE MEASURED
VALUE (U, P, Q, and I, with a checkbox enabling the current rows):
percent of reading keeps one accuracy setting meaningful across voltage
levels, where an absolute MW sigma cannot fit a 400 kV corridor and a
30 kV feeder at once (`generateMeasurementsFromPF(relativeSigma = true)`
with the `measurementSigmaFloors` floors guarding near-zero readings; the
0.05 MW/MVar power floor models the transducer range term and keeps
zero-injection buses from stiffening the flat-start solve).
Seeded noise is ON by default: without it the measurements match the
model exactly and `J` lands near 0 instead of near `dof`, which reads
like a broken statistic on first sight (the band note then says
`:low`, "not an alarm"). The noise plus a bad-data injection (k times
sigma; 10 is
clearly detectable) provide a reproducible test vector for the
elimination and robust workflows: the "bad data rows" field says HOW MANY
measurements get the gross error, the seed decides WHICH rows (drawn from
the telemetry rows; protected `ZI` constraints are never corrupted), and
the generation message names every corrupted row. A tap-deviation option
generates the measurements from a state whose transformers run a chosen
WHOLE number of mechanical tap steps off the model position (a tap
changer has no half positions; the model file stays untouched): the
"transformers (max)" field caps how many transformers get the deviation,
the seed decides which. Drawn are ONLY transformers whose deviation the
tap estimation can absorb (non-machine, with a declared changer); the cap
lands there and the generation message notes when fewer transformers were
eligible than requested. A machine (generator step-up) transformer is
used only when the case has nothing else, because its deviation is
skipped by the mass release by design and would leave an unexplainable
model error (`J` far above `dof`); the estimator run warns in that case.
Each hit is named in the generation message. The same seed always reproduces the identical set, including
the same bad-data rows and the same transformers. Off-grid positions
still occur in real data, and the estimator reports them honestly through
the off-grid note after the fixation; they are just not a generator
setting.

### Measurement generator v2

The generator's truth state, flow placement, and passive-node handling are
configurable in the Web UI form:

**Truth state.** `fresh solve` (default) solves the case now (tolerance
tightened to at most 1e-8, island-wise); the confirmation names the
iterations and tolerance of the pre-solve. `from run` adopts the solved
voltages of a successful run of the SAME case from the run history: SE
runs through `se_state.csv`, power-flow runs through the detail CSV
`bus_voltages_complex.csv` (any of the CSV export formats). Nothing is
re-solved, the adopted values are bit-exact against the run artifact, and
a missing artifact, a foreign case, or an unknown run id rejects up front
with a readable message. The tap-deviation option is locked with `from
run` (a run state is a finished snapshot); the service enforces the lock
even when the form graying is bypassed.

**Flow measurements per branch.** `both ends` (default) measures P/Q (and
I) at from AND to. `one end (balance-aware)` keeps exactly one flow group
per branch; the end at a bus WITH injection telemetry is preferred (`ZI`
pseudo-rows do not count), and the from end wins when both or neither
qualify. The choice is deterministic and documented per branch in the
set's `# flow_end,...` comments.

**Passive nodes.** Buses without generation, load, and shunt
(`findPassiveBuses`) get explicit zero balance rows `Pinj/Qinj = 0` at the
configured sigma (default 0.05 MW/MVar). Small sigmas at many passive
nodes stiffen the flat start; for hard balances enable the checkbox
instead, which writes protected zero-injection constraints
(`addZeroInjectionMeasurements!`, prefix `ZI`, excluded from elimination and robust modification, at the configured sigma but never tighter than `ZERO_INJECTION_SIGMA` = 0.001 MW) and no duplicate injection rows. The 1 kW floor is measured, not chosen: at
1e-6 such a row weighs a million times a normal power measurement, and the
squaring in `G = H' W H` then pushes the normal equations past double
precision. On the 25000-bus set that showed as an estimate that did not
converge at all, `J/dof = 3e17` and zero-injection residuals of 182 GW, the
exact opposite of what the tight sigma was meant to enforce.

**Delta file.** Generated sets carry the noise-free truth value of every
row as `# truth_value,...` comments. An SE run on such a set writes
`se_deltas.csv`: per measurement the measured, truth, and estimated
values with all three deltas (measured-truth, estimated-truth,
estimated-measured), the sigma, the normalized residual of the estimate,
and the elimination flag, plus one row per released tap comparing the
set's documented deviation with the fixed step. The set's provenance
(truth source, flow-end choices, passive handling) lands in the
`# generator: v2` comment block.

Naming triad of the service layer: the request key is `se_mode`, the SE
run's metadata reports `run_mode = "se"`, and the chained power flow reports
`run_mode = "powerflow_se_start"` (with `se_run_id` and `se_start_mode`);
an SE-started N-1 keeps `run_mode = "contingency"` plus the same two SE
reference keys.

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

gobs = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
println("Global observability quality: ", gobs.quality)

se = runse!(net; maxIte = 12, tol = 1e-6, flatstart = true, jacEps = 1e-6, updateNet = true)
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

obs = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
println("Observable quality: ", obs.quality)

se = runse!(net; maxIte = 12, tol = 1e-6, flatstart = true, jacEps = 1e-6, updateNet = true)
println("Converged: ", se.converged)
```

## Adding measurements with helper functions

Instead of constructing each `Measurement(...)` manually, you can build the
measurement vector with helper functions that resolve bus names and branch
references for you:

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

obs = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
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

* A. Abur, A. Gómez Expósito: *Power System State Estimation — Theory and
  Implementation* (hybrid SCADA/PMU WLS formulation).
* A. G. Phadke, J. S. Thorp: *Synchronized Phasor Measurements and Their
  Applications* (PMU measurement principle, IEEE C37.118 accuracy classes).
* NASPI TR-006: *Phase Angle Calculations: Considerations and Use Cases*
  (reference-angle handling across PMU installations).

## Case files carry their own measurements

A [Sparlectra Case Format](scf.md) case stores its measurement set inside the
file. An estimation on such a case therefore needs no separate CSV: the rows
that came with the case are used, and the run still writes `measurements.csv`
into its output directory so the result stays reproducible from its own
artifacts. Naming a measurement set explicitly still works and takes
precedence. Every other format continues to require one.

A set generated for the case takes precedence over the one the file carries:
generating is a deliberate act and the newer data. The generated set REPLACES
what the case brought along instead of being added to it.

**Ideal versus measured values.** A measurement set records whether it
carries noise, and a case file keeps that statement in
`sparlectra.measurements.provenance`. This matters because a noise-free set
returns $J = 0$ by construction: nothing contradicts anything, and the zero
says nothing about the estimate. The state-estimation page therefore names
the state of the carried set, and **Add noise to this set**
([`addMeasurementNoise!`](@ref)) turns ideal values into realistic ones
without recomputing anything: each value is perturbed with the sigma its own
row declares, the operating point and the sigmas stay as they are. On the
warm-up case that moves $J$ from 0 to about 43 at 42 degrees of freedom.

The Web UI never preselects a measurement set that belongs to a different
case: only a set bound to the selected case is offered as the default, a case
file with its own measurements offers those first and says how many it
carries, and otherwise the page asks for a set instead of arming one that can
only fail the binding check. A set bound to the case an SCF file was exported
FROM stays usable, because the file records its own source.

The provenance carries more than the noise statement: the per-row truth
values and the tap deviations the generator applied travel with the case as
well. They are what lets a run on a case file produce the same
`se_deltas.csv` and the same tap warning a run from a CSV set produces.

**A quantity measured twice.** Redundant transducers are normal in
telemetry, but a set in which EVERY quantity appears twice is a corrupted
file, usually one generation appended onto rows that were already there. Such
a set inflates $J$ without any single residual looking wrong, which reads as
a topology error and sends the diagnosis in the wrong direction. Reading a
set therefore reports how many rows repeat an already measured quantity
(active rows only), and the run states it in the log and in its message.
Regenerating the set is the fix.

## Further examples and workshop material

* Extended tutorial and a simple 7-bus setup: the state-estimation chapter of the [advanced workshop tour](generated/workshop_tour_advanced.md) and the dedicated [state-estimation notebook](generated/workshop_state_estimation.md)
* Bad-data localization, sequential elimination, robust estimation, and shunt-parameter estimation, hands-on: the [SE diagnostics notebook](generated/workshop_se_diagnostics.md)
* Tap estimation (ratio, phase, cascade, and the case that is not identifiable), current measurements, PMU phasors, the topology advisory, island-wise estimation, and $J$ versus $J_{active}$: the [taps, phasors and topology notebook](generated/workshop_se_taps.md), which closes with a decision tree for "$J$ is too large, what now"
* Detailed WLS reporting example script: `examples/state_estimation/state_estimation_wls.jl`
* PMU angle measurements and the reference-offset state α: `examples/state_estimation/state_estimation_pmu_angles.jl`
* Observability-focused scenario script: `examples/state_estimation/state_estimation_observability.jl`
* Passive-bus ZIB comparison example: `examples/state_estimation/state_estimation_passive_bus_zib_comparison.jl`
* Matrix-based observability/redundancy demo: `examples/state_estimation/h_matrix_observability_demo.jl`

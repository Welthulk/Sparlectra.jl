# Observability

Observability answers one question before any estimate exists: does the
measurement set determine the network state, and if not, where does it
fail? The checks run standalone, before and independently of `runse!`.

## What the global check answers

`evaluate_global_observability(net; ...)` decides whether the state is
determined by the active measurements in `net.measurements`: all bus
angles and magnitudes, the circulation states of measured coupler rings,
and the PMU reference offset when it is active. Current rows (`ImagMeas`,
`IaMeas`) are left out, and the voltage range check of `runse!` applies.
Released shunt susceptances and tap steps are not part of this check.
Whether a tap step is observable is decided per changer by `runse!` before
its first iteration, see
[Observability per changer](@ref se-tap-observability).

A net with more than one network island (islands joined by a closed link
count as one) is judged per measured island, on its own subnet with its own
angle reference. The quality is that of the worst measured island. `notes`
holds `:island_partition`, `islands` holds the result per island, and
unmeasured islands are counted and left out.

The verdict comes from the decoupled models
([below](@ref se_observability_gaps)). The AC measurement Jacobian gives
`numerical_rank`, the criticality and the weak class. It is built by
forward differences, so its error floor is of the order of `jacEps` rather
than machine epsilon. With `tol = nothing` the rank tolerance is

```math
\text{tol} = f \cdot \varepsilon_J \cdot \sigma_{\max}
```

with `f = state_estimation.rank_tol_factor` (default 10.0) and `jacEps` as
$\varepsilon_J$; an explicitly passed `tol` wins. The rows and then the
columns of the Jacobian are normalized to unit norm, and $\sigma_{\max}$ is
the largest singular value of that matrix. Column normalization makes the
tolerance mean the same on every network (a voltage column carries entries
around 1, an angle column the admittances of a line). Rows are normalized
instead of divided by their sigma: the sigma would lift zero-injection rows
orders of magnitude above the others, and without the row step a bus
coupler entered as a branch with `x = 1e-5` pu dominates the columns of its
two buses. The local check (`evaluate_local_observability`) takes its
submatrix from the raw Jacobian and its $\sigma_{\max}$ from that submatrix.

The result reports the measurement count `m`, the state count `n`, the
redundancy and the redundancy ratio

```math
r = m - n, \qquad \rho = \frac{m}{n},
```

the structural and numerical flags, and the quality label.

## The quality labels

- `:good`: full rank with usable redundancy; the estimate and its
  diagnostics are trustworthy.
- `:critical`: full rank but with critical measurements in the set (see
  below); the estimate exists, parts of it are unprotected against bad
  data.
- `:weakly_observable`: no null direction, but a branch flow has a
  standard deviation $\sqrt{g_k^{\mathsf T} G^{-1} g_k}$ above
  `state_estimation.observability.weak_flow_sigma_mw` (20 MW) or a
  magnitude one above `weak_vm_sigma_pu` (0.01 pu), from the covariance of
  the estimate. Before the estimation the class is provisional, at the
  stored state. `runse!` determines it at the estimate
  (`weakObservability`): under `abort` after a converged run with
  `update_net` unless `weak_check_after_estimation = false`, under
  `observable_islands` and `scheduled` always.
- `:not_observable`: a null direction of the decoupled P-theta or Q-V
  model; see [When the set is not observable](@ref se_observability_gaps).

The answer to `:not_observable` is a measurement. The decoupled models use
an exact rank tolerance, so `rank_tol_factor` does not move the verdict.
`unobservable_state_columns` names the states the set does not pin down
(per network island in `islands` when the net has several); a flow or
injection measurement at one of those places fixes the cause.

## [When the set is not observable](@id se_observability_gaps)

**Islands.** A network island is a topologically separate part with its
own reference and no gap. An observable island lies inside one: the buses
whose states the measurements fix once one of them is the reference,
labeled `k.j`. Not observable is decided on the decoupled linear P-theta
model (P rows and PMU angles against the angles) and Q-V model (Q and Vm
rows against the magnitudes). `identify_observable_islands` (printed by
`print_observable_islands`) returns the islands after Monticelli and Wu,
the unobservable branches between them, angle and magnitude deficit and
the region of dark angles. An unmeasured branch in one island is observable.

**Estimating the islands one by one.** With
`state_estimation.observability.on_unobservable = observable_islands`
`runse!` estimates each island on its own; an island without the slack
has relative angles (`observableIslands`). Dark-magnitude and weak buses
leave their island. Unestimated buses and unobservable flows are NaN.
When no island can be estimated, `runse!` throws `SENotObservableError`.
Under the default `abort` the library does not check; the SE service and
the Web UI stop and name the islands.

**Closing with scheduled values.** `on_unobservable = scheduled` first
takes replacement values (`restore_observability`): zero injections at
passive buses, scheduled P and Q, then U from the stored state, a set
point, the median of the measured U of the level or the nominal voltage,
each only when it lowers a deficit. P and Q get the sigma
`restore_sigma_rel` (30 %) of the value, at least that of a meter; U the
sigma of a meter, the spread of the level for the median and
`restore_vm_sigma_nominal_pu` at the nominal voltage. The rows
(`restoredRows`, origin `:pseudo_restore`) are never eliminated; the Web UI
shows the run yellow. A branch or bus left above the weak bounds gives
"weakly observable after restoring". A replacement value is mostly
critical, no residual shows its error.

**Suggesting meter places.** `suggest_measurement_placement(net, meas;
kind)` adds one meter per step, within a network island, up to
`placement_max_steps`: `injection` (P, Q or U), `flow` (P and Q at a branch
end) or `pmu` (an angle; the first one may only fix the PMU offset). For
`injection` and `pmu` the steps equal the deficit (plus the offset step),
the minimum; `flow` is greedy and not guaranteed minimal. Examples 13 to 16
of the [state-estimation workshop](generated/workshop_state_estimation.md)
and `examples/state_estimation/state_estimation_observable_islands.jl`.

## What the local check answers

`evaluate_local_observability(net, cols; ...)` asks the same question
for a chosen subset of state columns, for example one bus angle and one
magnitude: the instrument for manual sensor and PMU placement studies.

## Critical measurements, and why their residual is zero

A measurement is critical when removing it makes the system
unobservable: its information exists nowhere else in the set. The
estimate then reproduces it exactly, its residual is zero regardless of
its error, and no residual-based test can flag it: a gross error in a
critical measurement lands fully in the state. The diagnostics report
classifies critical versus redundant measurements for this reason.

### How the critical measurements are found

The residual sensitivity matrix

```math
\Omega = R - H\,G^{-1}H^{\mathsf T}, \qquad G = H^{\mathsf T} W H
```

says how much of a measurement error becomes visible in that
measurement's own residual; for a critical measurement `Omega_ii` is
exactly zero. Only the diagonal of `H G^-1 H'` is needed, never the full
inverse; the Takahashi selected inverse delivers it from the
factorization of `G` (Takahashi above `takahashi_min_states` states, the
dense path below). Weights do not move the zeros, only the scale in
between, so the observability check runs with unit weights and no sigma
has to be known.

The threshold is dimensionless on the share of a row's own error that
reaches its residual,

```math
w_{ii} = \Omega_{ii}\, w_i ,
```

and a row is critical at

```math
w_{ii} \le \max\!\left( \left(\frac{\text{tol}}{\sigma_{\max}}\right)^{2},\; 10^{-8} \right),
```

the FD-aware rank tolerance made relative to `sigma_max` and squared,
floored at `1e-8`: far above the rounding of the selected inverse and far
below any real redundancy. The floor is needed because a structurally
critical row can come back at `1e-12` from fill-in and rounding instead of
at zero.

The per-row rank tests remain as `state_estimation.criticality_method =
rank` for a cross-check: one rank decision per measurement row (a dense
SVD on small sets, a sparse QR on large ones), budgeted at 300000 rows times
states and skipped above it (`criticality_skipped = true`). The selected
inverse answers the same question from one factorization for all rows,
has no size budget, and grades the information: `Omega_ii = 0` means
critical, a normalized `wii` below the 0.3 guideline nearly critical.

### Rank source

The report names the rank source in `rank_method`.
`state_estimation.rank_method = decomposition` (the default) is the SVD of
the Jacobian on small sets; on large ones a sparse QR finds a deficit and
the null-space iteration on the gain matrix sizes it. `pivots` reads the rank from the
LDLt factorization of the gain matrix `H'H`, which the criticality pass
needs anyway, as the number of pivots above the squared rank tolerance;
exact for a full-rank Jacobian, not on a deficit (the pivots of an LDLt are
not the eigenvalues of the gain matrix, and a singular gain matrix cannot
be factorized at all), so every deficit goes to the decomposition and is
reported as `rank_method = decomposition`. Both settings state the same
rank on every case; choose `pivots` where the decomposition is the measured
cost on a very large state vector. The result carries `criticality_wii`
per active row and `criticality_method`; the run log lists the critical
and the nearly critical rows.

## `w_ii`, the localizability indicator

`wii` is the diagonal of the residual sensitivity matrix: the share of
measurement `i`'s own error that reaches its own residual. `wii` near 1
means an error there shows up almost fully in the residual (well
localizable); `wii` near 0 marks a nearly critical measurement whose
error hides in the state estimate. The report flags
`localizable = wii > wiiThreshold` with the threshold 0.3 from the
bad-data literature: below it, the largest-normalized-residual logic
cannot be trusted to point at the faulty meter.

## The correlation bound for singly redundant groups

The optional K-matrix report marks rows whose normalized residual
correlates with a partner above `1/sqrt(2)` (about 0.707): a simply
redundant group, in which a gross error cannot be attributed to one of
the two. Formula, configuration and result fields:
[Residual correlations](@ref se-k-report).

## How zero-injection buses enter

The zero-injection rows of passive buses are rows of H like measurements
and raise observability and redundancy around the bus. That holds for
weighted rows and for equality constraints alike; how the two enter the
estimate is described in
[Zero injections and coupler balances](@ref se_zero_injection).

## Reading the diagnostics table

Each active measurement gets one row: type, location, sigma, the
normalized residual, `wii` with the localizable flag, the
critical/redundant classification and, when the K report is on, the
correlation mark. Read the band test first, then the largest normalized
residual among localizable rows, then `wii` before trusting any single-row
verdict; treat critical rows as unprotected rather than clean. Ranking
semantics, elimination trace and robust interplay:
[State Estimation](state_estimation.md).

## The H matrix as the didactic entry

`measurement_jacobian(net; ...)` returns the labeled Jacobian behind
both checks: `H` with one described row per active measurement (a link
flow row has the location `link <k>`) and named state columns (`Va(bus)`,
`Vm(bus)`, then `circP(j)` and `circQ(j)` for the circulation states of
measured coupler rings, plus `alpha` when PMU angle measurements activate
the reference-offset state). It shows which rows touch which columns
before any SVD runs. The example suite
`examples/run_state_estimation_suite.jl` writes such a measurement-matrix
page with a stability verdict from rank, redundancy and `cond(H)`; the
[state estimation workshop](generated/workshop_state_estimation.md)
walks the same matrix on a small network.

# Observability

Observability answers one question before any estimate exists: does the
measurement set determine the network state, and if not, where does it
fail? The checks run standalone, before and independently of `runse!`.

## What the global check answers

`evaluate_global_observability(net; ...)` decides whether the complete
state (all bus angles and magnitudes, plus released extra states) is
determined by the active measurements in `net.measurements`, in two
stages:

1. **Structural stage**: `detect_ac_islands` runs on the contracted SE
   net. More than one synchronous island containing measured buses
   yields `:not_observable` with `:structural_islands` in `notes`,
   because no measurement ties the islands' angle references together.
2. **FD-aware numeric rank**: the measurement Jacobian is built by
   forward differences, so its error floor is of the order of `jacEps`
   rather than machine epsilon. With `tol = nothing` the rank tolerance
   defaults to

   ```math
   \text{tol} = f \cdot \varepsilon_J \cdot \sigma_{\max}
   ```

   with `f = state_estimation.rank_tol_factor` (default 10.0) and
   `jacEps` as $\varepsilon_J$; an explicitly passed `tol` always wins.
   The global rank runs on a Jacobian whose rows and then columns are
   normalized to unit norm, and $\sigma_{\max}$ is the largest singular
   value of that normalized matrix; rows are not divided by their sigma.
   The local check (`evaluate_local_observability`) takes its submatrix
   from the raw finite-difference Jacobian and its $\sigma_{\max}$ from
   that submatrix.

The result reports the measurement count `m`, the state count `n`, the
redundancy and the redundancy ratio

```math
r = m - n, \qquad \rho = \frac{m}{n},
```

the structural and numerical flags, and the quality label.

Why normalize: against the raw matrix the tolerance means different things
on different networks (a voltage column carries entries around 1, an angle
column the admittances of a line), and on a large network genuine small
singular values drop below a cut measured against the largest column.
Column normalization removes that dependence, and the verdict becomes
insensitive to the factor over several orders of magnitude. Rows are not
divided by their sigma, which would lift zero-injection pseudo-measurements
(sigma `1e-6`) six orders of magnitude above ordinary rows; they are
normalized to unit norm, before the columns, because a bus coupler entered
as a branch with near-zero impedance (`x = 1e-5` pu) otherwise dominates
the columns of its two buses and pushes the feeding line below the cut.

## The quality labels

- `:observable`: full rank with usable redundancy; the estimate and its
  diagnostics are trustworthy.
- `:critical`: full rank but with critical measurements in the set (see
  below); the estimate exists, parts of it are unprotected against bad
  data.
- `:not_observable`: rank deficiency or structural islands; no estimate
  determines the whole state, and the notes say which stage failed.

The answer to `:not_observable` is a measurement, not a tolerance.
`unobservable_state_columns` names the states the set does not pin down;
a flow or injection measurement at one of those places fixes the cause.
Lowering `rank_tol_factor` only moves the line at which a direction
counts as present, and the estimate then silently rests on the start
values in the unmeasured directions.

## What the local check answers

`evaluate_local_observability(net, cols; ...)` asks the same question
for a chosen subset of state columns, for example one bus angle and one
magnitude: the instrument for manual sensor and PMU placement studies.
There is no observable-island decomposition, no automatic
pseudo-measurement restoration and no placement optimizer, see the
[feature matrix](feature_matrix.md).

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
SVD up to 2000 states, a sparse QR above), budgeted at 300000 rows times
states and skipped above it (`criticality_skipped = true`). The selected
inverse answers the same question from one factorization for all rows,
has no size budget, and grades the information: `Omega_ii = 0` means
critical, a normalized `wii` below the 0.3 guideline nearly critical.

### Rank source

The report names the rank source in `rank_method`.
`state_estimation.rank_method = decomposition` (the default) is the SVD of
the Jacobian below 2000 states and the sparse QR factorization above: the
exact numerical rank, deficit included. `pivots` reads the rank from the
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

A passive bus (no generation, no load, no shunt) contributes exact
knowledge: its injections are zero. Sparlectra models that as tightly
weighted zero-injection pseudo-measurements, which raise observability
and redundancy around the bus; they are protected from elimination and
robust down-weighting, because removing exact knowledge is never the
right reaction to a residual. The weight is finite (a hard-constraint
solver block is not implemented, see the feature matrix), so extreme
weights trade against conditioning.

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
both checks: `H` with one described row per active measurement and
named state columns (`Va(bus)`, `Vm(bus)`, plus `alpha` when PMU angle
measurements activate the reference-offset state), which shows which
rows touch which columns before any SVD runs. The example suite
`examples/run_state_estimation_suite.jl` writes such a measurement-matrix
page with a stability verdict from rank, redundancy and `cond(H)`; the
[state estimation workshop](generated/workshop_state_estimation.md)
walks the same matrix on a small network.

# Observability

Observability answers one question before any estimate exists: does the
measurement set determine the network state at all, and if not, where
does it fail? State-estimation quality depends on it more than on any
solver knob, which is why the checks run standalone, before and
independently of `runse!`.

## What the global check answers

`evaluate_global_observability(net; ...)` assesses whether the complete
state (all bus angles and magnitudes, plus released extra states) is
determined by the active measurements in `net.measurements`. The check
is two-staged:

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
   The rank runs on a column-normalized Jacobian; rows are not divided
   by their sigma.

Typical metrics of the result: the measurement count `m`, the state
count `n`, the redundancy and the redundancy ratio

```math
r = m - n, \qquad \rho = \frac{m}{n},
```

the structural and numerical flags, and the quality label.

!!! details "Why the rank uses a column-normalized Jacobian and an FD-aware tolerance"
    An eps-scale SVD tolerance would count FD noise as rank and can miss
    a structurally unobservable direction, hence the tolerance is tied to
    `jacEps`.

    The tolerance is relative to `sigma_max` and would otherwise mean
    different things on different networks: a voltage column carries
    entries around 1, an angle column the admittances of a 765 kV line,
    and on a large network the genuine small singular values drop below
    a cut that is measured against the largest column. Column
    normalization removes that dependence, and with it the verdict is
    insensitive to the factor over several orders of magnitude; only a
    very large factor starts cutting into genuine directions.

    Rows are deliberately not divided by their sigma. Rank is invariant
    under positive row scaling, so sigma answers nothing here, while
    dividing by it lifts zero-injection pseudo-measurements (sigma
    `1e-6`) six orders of magnitude above ordinary rows and lets them
    dominate the tolerance.

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
adding a flow or injection measurement at one of those places fixes the
cause. Lowering `rank_tol_factor` only moves the line at which a
direction counts as present, so a set that genuinely lacks information
passes as observable and the estimate silently rests on the start
values in those directions.

## What the local check answers

`evaluate_local_observability(net, cols; ...)` asks the same question
for a chosen subset of state columns, for example one bus angle and one
magnitude. That is the placement instrument: instead of an optimizer,
Sparlectra lets you probe exactly the states you care about and see
whether the present set pins them, which is what manual sensor and
PMU placement studies need. The scope is stated in the
[feature matrix](feature_matrix.md): there is no observable-island
decomposition, no automatic pseudo-measurement restoration, and no
placement optimizer; local observability on a chosen subset plus
explicit critical-measurement reporting is the offered alternative.

## Critical measurements, and why their residual is zero

A measurement is critical when removing it makes the system
unobservable: its information exists nowhere else in the set. The
estimate then reproduces it exactly, so its residual is exactly zero
regardless of its error, and no residual-based test can ever flag it.
That is the practical reason redundancy matters: a gross error in a
critical measurement is invisible to the diagnostics and lands fully
in the state. The diagnostics report classifies critical versus
redundant measurements for exactly this reason.

### How the critical measurements are found

The residual sensitivity matrix

```math
\Omega = R - H\,G^{-1}H^{\mathsf T}, \qquad G = H^{\mathsf T} W H
```

says how much of a measurement error becomes visible in that
measurement's own residual. For a critical measurement `Omega_ii` is
exactly zero: the estimator adapts to the row completely, its residual
is always zero, an error there stays invisible. Only the diagonal of
`H G^-1 H'` is needed, never the full inverse; the Takahashi selected
inverse delivers that diagonal from the factorization of `G` (Takahashi
above `takahashi_min_states` states, the dense path below). Weights do
not move the zeros, only the scale in between, so the observability
check runs it with unit weights and no sigma has to be known.

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
below any real redundancy.

The per-row rank tests are available as
`state_estimation.criticality_method = rank` (budgeted at 300000 rows
times states, `criticality_skipped = true` above) as the cross-check;
the test suite runs both methods against each other on a boundary case
(a demo set with one flow measurement removed).

!!! details "Why the selected inverse replaces the per-row rank tests"
    The literal reading of the definition strikes every measurement row
    in turn and decides the rank of the remaining Jacobian again (a dense
    SVD up to 2000 states, a sparse QR above). With m measurements that
    is m rank decisions, which is why that method is budgeted by rows
    times states and skipped above the budget, leaving large networks
    without a criticality flag.

    `Omega_ii = 0` is the statement of the rank test, reached from the
    other side, and one factorization answers it for all rows at once:
    the classification costs a fraction of the rank tests; there is no
    size budget, every set the estimator itself can factorize gets its
    flag; the information is graded at no extra cost, because
    `Omega_ii = 0` means critical while the normalized `wii` below the 0.3
    guideline means nearly critical, where the rank test can only answer
    yes or no; and the diagnostics report and the classification compute
    the same quantity in one code path.

    The open point is the tolerance. The rank test is binary, `Omega_ii`
    is a floating-point number, and a structurally critical row can come
    back at `1e-12` from fill-in and rounding instead of at zero. Hence
    the dimensionless threshold on `wii` above.

### Rank source

The rank itself has two sources, and the report names the one used in
`rank_method`.

`state_estimation.rank_method = decomposition` (the default) is the SVD of
the Jacobian below 2000 states and the sparse QR factorization above. It
gives the exact numerical rank, deficit included, and on every shipped case
it costs about the same as the alternative.

`pivots` reads the rank from the LDLt factorization of the gain matrix
`H'H`, the factorization the criticality pass needs anyway, as the number
of pivots above the squared rank tolerance. That count is exact when the
Jacobian has full rank. On a deficit it is not (the pivots of an LDLt are
not the eigenvalues of the gain matrix, and a singular gain matrix cannot
be factorized at all), so the pivot rule hands every deficit to the
decomposition and reports `rank_method = decomposition` for it.

Both settings state the same rank on every case. Choose `pivots` where the
decomposition is the measured cost on a very large state vector;
otherwise the default is the exact answer at no extra cost.

The result carries `criticality_wii` per active row and
`criticality_method`; the state-estimation run log lists the critical
rows and the nearly critical ones. The rank decision for the
observability verdict itself: up to 2000 states the exact dense SVD,
above that a sparse QR factorization with the identical FD-aware
tolerance.

## `w_ii`, the localizability indicator

`wii` is the diagonal of the residual sensitivity matrix: the share of
measurement `i`'s own error that reaches its own residual. `wii` near 1
means an error there shows up almost fully in the residual (well
localizable); `wii` near 0 marks a nearly critical measurement whose
error hides in the state estimate. The report flags
`localizable = wii > wiiThreshold` with the literature threshold 0.3:
below it, the largest-normalized-residual logic cannot be trusted to
point at the truly faulty meter. The threshold is the established
guideline from the bad-data literature, not a tuned constant.

## The correlation bound for singly redundant groups

With the optional K-matrix report, each ranking row carries the maximum
absolute correlation of its normalized residual over all partners,

```math
K = D^{-1/2}\, \Omega\, D^{-1/2} .
```

Above `1/sqrt(2)` (about 0.707) two measurements form a simply redundant
group: a gross error in one is statistically indistinguishable from an
error in the other, and the report marks the row. This is reporting
only, never automatic action.

## How zero-injection buses enter

A passive bus (no generation, no load, no shunt) contributes exact
knowledge: its injections are zero. Sparlectra models that as tightly
weighted zero-injection pseudo-measurements, which raise observability
and redundancy around the bus; they are protected from elimination and
robust down-weighting, because removing exact knowledge is never the
right reaction to a residual. The weight is finite (a hard-constraint
solver block is deliberately not implemented, see the feature matrix),
so extreme weights versus conditioning is a real trade recorded in
[State Estimation](state_estimation.md).

## Reading the diagnostics table

Each active measurement gets one row: its identity (type, location,
sigma), the normalized residual against the residual covariance, the
`wii` column with the localizable flag, the critical/redundant
classification, and, when the K report is on, the correlation mark.
Read it in this order: the band test first (is the set consistent at
all), then the largest normalized residual among localizable rows,
then `wii` before trusting any single-row verdict, and treat critical
rows as unprotected rather than clean. The full ranking semantics,
elimination trace and robust interplay live in
[State Estimation](state_estimation.md).

## The H matrix as the didactic entry

`measurement_jacobian(net; ...)` returns the labeled Jacobian behind
both checks: `H` with one described row per active measurement and
named state columns (`Va(bus)`, `Vm(bus)`, plus `alpha` when PMU angle
measurements activate the reference-offset state). Seeing which rows
touch which columns makes structural observability concrete before any
SVD runs. The state-estimation example suite
(`examples/run_state_estimation_suite.jl`) writes such a
measurement-matrix page, including a stability verdict from rank,
redundancy and `cond(H)`, on every run; the
[state estimation workshop](generated/workshop_state_estimation.md)
walks the same matrix on a small network and shows observability
breaking down live.

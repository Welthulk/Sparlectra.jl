# Observability

Observability answers one question before any estimate exists: does the
measurement set determine the network state at all, and if not, where
does it fail? State-estimation quality depends on it more than on any
solver knob, which is why the checks run standalone, before and
independently of `runse!`.

## What the global check answers

`evaluate_global_observability(net; ...)` assesses whether the COMPLETE
state (all bus angles and magnitudes, plus released extra states) is
determined by the active measurements in `net.measurements`. Since
0.10.0 the check is two-staged:

1. **Structural stage**: `detect_ac_islands` runs on the contracted SE
   net. More than one synchronous island containing measured buses
   yields `:not_observable` with `:structural_islands` in `notes`,
   because no measurement ties the islands' angle references together.
2. **FD-aware numeric rank**: the measurement Jacobian is built by
   forward differences, so its error floor is O(`jacEps`) rather than
   machine epsilon. An eps-scale SVD tolerance would count FD noise as
   rank and can miss a structurally unobservable direction; with
   `tol = nothing` the rank tolerance therefore defaults to
   `state_estimation.rank_tol_factor * jacEps * sigma_max` (factor
   default 10.0), and an explicitly passed `tol` always wins.

   The rank runs on a **column-normalized** Jacobian, because that
   tolerance is relative to `sigma_max` and would otherwise mean
   different things on different networks: a voltage column carries
   entries around 1, an angle column the admittances of a 765 kV line,
   and on a large network the genuine small singular values drop below
   a cut that is measured against the largest column. Rows are
   deliberately **not** divided by their sigma. Rank is invariant under
   positive row scaling, so sigma answers nothing here, while dividing
   by it lifts zero-injection pseudo-measurements (sigma `1e-6`) six
   orders of magnitude above ordinary rows and lets them dominate the
   tolerance. Measured on a 25000-bus set of 203916 rows for 49999
   states: rank deficit 191 unscaled, 12308 with row scaling, 0 with
   column normalization alone. With the columns normalized the verdict
   is insensitive to the factor: across `sp_case5` to that 25000-bus
   set the deficit stays 0 from factor 100 down to 0.01, and only a
   factor of 1000 starts cutting into genuine directions.

Typical metrics of the result: the measurement count `m`, the state
count `n`, redundancy `r = m - n` and the redundancy ratio
`rho = m / n`, the structural and numerical flags, and the quality
label.

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
for a CHOSEN subset of state columns, for example one bus angle and one
magnitude. That is the placement instrument: instead of an optimizer,
Sparlectra lets you probe exactly the states you care about and see
whether the present set pins them, which is what manual sensor- and
PMU-placement studies need. The honest scope is stated in the
[feature matrix](feature_matrix.md): there is no observable-island
decomposition, no automatic pseudo-measurement restoration, and no
placement optimizer; local observability on a chosen subset plus
explicit critical-measurement reporting is the offered alternative.

## Critical measurements, and why their residual is zero

A measurement is CRITICAL when removing it makes the system
unobservable: its information exists nowhere else in the set. The
estimate then reproduces it exactly, so its residual is exactly zero
regardless of its error, and no residual-based test can ever flag it.
That is the practical reason redundancy matters: a gross error in a
critical measurement is invisible to the diagnostics and lands fully
in the state. The diagnostics report classifies critical versus
redundant measurements for exactly this reason.

The classification is one rank test per measurement row, which stops
paying for itself on large systems: above a measured budget (rows
times states beyond 300000; the 188-bus demo case already spends
seconds here) the check is skipped with a warning, the result carries
`criticality_skipped = true`, and the quality label then reflects
observability and redundancy only. The rank decision itself scales:
up to 2000 states it is the exact dense SVD it always was, above that
a sparse QR factorization with the identical FD-aware tolerance.

## `w_ii`, the localizability indicator

`wii = Omega_ii * w_i` is the diagonal of the residual sensitivity
matrix: the share of measurement `i`'s own error that reaches its own
residual. `wii` near 1 means an error there shows up almost fully in
the residual (well localizable); `wii` near 0 marks a nearly critical
measurement whose error hides in the state estimate. The report flags
`localizable = wii > wiiThreshold` with the literature threshold 0.3:
below it, the largest-normalized-residual logic cannot be trusted to
point at the truly faulty meter. The threshold is the established
guideline from the bad-data literature, not a tuned constant.

## The correlation bound for singly redundant groups

With the optional K-matrix report, each ranking row carries the maximum
absolute correlation of its normalized residual over all partners,
`K = D^(-1/2) Omega D^(-1/2)`. Above `1/sqrt(2)` (about 0.707) two
measurements form a simply redundant group: a gross error in one is
statistically indistinguishable from an error in the other, and the
report marks the row. This is reporting only, never automatic action.

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
all), then the largest normalized residual AMONG LOCALIZABLE rows,
then `wii` before trusting any single-row verdict, and treat critical
rows as unprotected rather than clean. The full ranking semantics,
elimination trace and robust interplay live in
[State Estimation](state_estimation.md).

## The H matrix as the didactic entry

`measurement_jacobian(net; ...)` returns the labeled Jacobian behind
both checks: `H` with one described row per active measurement and
named state columns (`Va(bus)`, `Vm(bus)`, plus `alpha` when PMU angle
measurements activate the reference-offset state). It is the didactic
entry point: seeing which rows touch which columns makes structural
observability concrete before any SVD runs. The state-estimation
example suite (`examples/run_state_estimation_suite.jl`) writes such a
measurement-matrix page, including a stability verdict from rank,
redundancy and `cond(H)`, on every run; the
[state estimation workshop](generated/workshop_state_estimation.md)
walks the same matrix on a small network and shows observability
breaking down live.

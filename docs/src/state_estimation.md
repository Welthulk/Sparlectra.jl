State Estimation
=================

Weighted least squares (WLS) state estimation reconstructs the network state
from redundant, noisy measurements. This page covers the estimator, its
measurement types, the bad-data diagnostics and the result report.
Extensions (FACTS view, links, islands, taps, topology):
[State Estimation Extensions](state_estimation_extensions.md). Files,
generator, chained power flow:
[State Estimation Measurements](state_estimation_measurements.md). Rank and
critical measurements: [Observability](observability.md). Keys and
defaults:
[State-Estimation Configuration](state_estimation_configuration.md).

## Theory (compact)

* State vector: `x = [θ(non-slack); Vm(all buses)]`
* Measurement model: `z = h(x) + e`
* Objective: `J(x) = (z - h(x))' * W * (z - h(x))`

`z` is the measurement vector, `h(x)` the prediction from the network
model, `W = diag(1/σ²)` the inverse-variance weighting. The solver
linearizes `h(x)` and iterates Newton-style until the update norm or the
residual criteria meet the tolerance. PMU voltage-angle measurements add
the reference-angle offset α: `x = [θ(non-slack); Vm(all buses); α]`.

### [Start value](@id se_start_value)

The Gauss-Newton iteration starts from the flat start
(`state_estimation.flatstart = true`, the default): every |V| at 1 pu,
the slack magnitude from the net, and every angle 0. Inside the estimate the
slack sits at 0; the result voltages (`SEResult.voltages`), the write-back
(`update_net`), the CSV files and the Web UI show the angles in the
reference of the model, so the slack keeps its stored angle (for example
the SV angle of a CGMES import) and a flat start is every angle at the
slack angle there. With `flatstart = false` the warm start takes the state
stored in the net against the stored slack angle: the last estimate written
back, a power flow result or the imported state. A start from measured
voltages is not offered.

A synthetic study generates its measurements from a power flow whose
voltages are still in the net; a warm start then begins on the solution
and says nothing about convergence. The examples therefore start flat:
with `flatstart = true`, or with a hand reset of the non-slack buses to
1 pu and 0 deg before a warm start. Examples that only evaluate the
Jacobian at a given state (condition number, diagnostics benchmark) use the
warm start on purpose.

Start ladder (`state_estimation.start_ladder`, default `true`): when the
run from the flat start does not converge, `runse!` solves once more from
the voltages stored in the net, the case state of the import (MATPOWER
Vm/Va, CGMES SV voltages, the SCF start state). A net whose stored voltages
equal the flat start gets no second try. `SEResult.start` names the start
of the result (`flat`, `stored`, `stored_after_flat`), the SE service writes
one run.log line and the metadata `se_start` (with a staged start also
`se_tap_staged_start_from`). On case6495rte with a generated set the flat
start did not converge in 50 iterations and the case voltages converged in
3.

With `state_estimation.start_mode.rotate_behind_pst = true` (default
`false`) the flat start turns a bus that the slack reaches only across
phase-shifting transformers by their shift sum, the same angles as the
power flow flat start under `power_flow.start_mode.rotate_behind_pst`
([Start strategies](start_strategies.md)); a bus that a path without
shifter reaches keeps 0.

## Why the FD measurement Jacobian works

The measurement Jacobian `H = ∂h/∂x` is built by finite differences, column
by column:

```math
\frac{\partial h}{\partial x_k}(x)
\approx
\frac{h(x + \varepsilon e_k) - h(x)}{\varepsilon}.
```

WLS needs only the local first-order sensitivity
($h(x + \delta x) \approx h(x) + H(x)\,\delta x$), so the estimation model is
the same as with an analytic Jacobian; only the derivative is evaluated
numerically.

## Measurement model

Measurement types:

* `VmMeas` (bus voltage magnitude, p.u.)
* `VaMeas` (bus voltage angle in degrees, PMU synchrophasor)
* `PinjMeas`, `QinjMeas` (bus injections, MW/MVar)
* `PflowMeas`, `QflowMeas` (branch flows with direction, MW/MVar; on a closed
  bus link see [Link flow measurements](@ref se_link_measurements))
* `ImagMeas` (branch current magnitude in ampere, with direction)
* `IaMeas` (PMU current angle), `ShuntQMeas` (shunt bay reactive power)

Synthetic sets come from a solved power flow through
`generateMeasurementsFromPF`, see
[Measurement generator](state_estimation_measurements.md#Measurement-generator).

### PMU voltage-angle measurements (`VaMeas`)

A PMU reports the bus voltage phasor against its GPS clock (IEEE C37.118),
so all PMU angles share one unknown rotation against the angle reference
of the model.
Sparlectra estimates that rotation as an additional state, the
reference-angle offset `α`, the angle of the PMU time base against the
angle reference of the model.

**Use**

| | |
|---|---|
| Add | `addVaMeasurement!(net; busName, value, sigma)` (degrees); `addPmuPhasorMeasurement!(net; busName, vm_pu, va_deg)` writes the Vm+Va pair |
| Generate | `generateMeasurementsFromPF(includeVa = true)`; default sigma `measurementStdDevs(va = 0.02)` |
| Config key | `state_estimation.pmu_ref_offset` (keyword `pmuRefOffset`): `:auto` (default) adds the offset state as soon as active `VaMeas` rows exist and gives `α` of about 0 when the time bases coincide, `:off` takes PMU angles in the angle reference of the model |
| Result field | `SEResult.vaRefOffsetDeg` |

**Notes**

- With `:off` and a wrong time base, `J` inflates by orders of magnitude
  and every PMU angle row lands in the bad-data ranking.
- The state count becomes `n = (2·nbus − 1) + 1 = 2·nbus`; the first PMU
  angle measurement pins `α` and is critical on its own, redundancy for `α`
  needs at least two. A `VaMeas` at the slack bus reads the stored slack
  angle plus `α`.
- Values and sigmas are in degrees. Angle residuals wrap at the plus/minus
  180 degree seam (`VaMeas` and `IaMeas`).
- PMU magnitudes are ordinary `VmMeas` rows with a small sigma (e.g. 0.002
  p.u.).

The offset state absorbs the PMU rotation,

```math
z_{Va,i} = \theta_i + \alpha + e_i,
```

so the Jacobian gains one column (`∂z_{Va,i}/∂α = 1` on PMU angle rows,
zero elsewhere), with `θ_i` in the angle reference of the model. PMU angle
rows carry a very large weight `1/σ²` but remain ordinary measurements,
not hard constraints.

### Branch current-magnitude measurements (`ImagMeas`)

`ImagMeas` models SCADA branch current magnitudes in ampere at one branch
end, predicted as `I = 1000 · |S| / (√3 · Vn · |V|)` from the branch flow at
that end. Currents are auxiliary: they never replace power measurements;
they raise the residual sensitivities `w_ii` of the power rows and so make
bad data localizable.

**Use**

| | |
|---|---|
| Add | `addImagMeasurement!(net; ...)` with `direction = :from`/`:to`; the bus-referenced variant `addImagMeasurement!(busName = ...)` (no direction) measures the shunt-bay current |
| Generate | `generateMeasurementsFromPF(includeImag = true)`; default sigma `measurementStdDevs(imag = 10.0)` |
| Config key | `state_estimation.imag_activation_iteration` (default 2, keyword `imagActivationIteration`): the first WLS iteration that takes current rows, so the first linearization from a flat start runs on power and voltage rows only |

**Notes**

- No observability: both observability checks drop `ImagMeas` rows before
  building the Jacobian; the system must be observable from power and
  voltage measurements alone.
- Value gate: a current with `value < 3 sigma` is excluded for the whole run
  (one info line per exclusion); the derivative of `|I|` is discontinuous
  near zero current.
- Bus-referenced currents feed the shunt back-calculation (case B below).

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

An `IaMeas` row is active only while the predicted magnitude passes
`3 * sigma_I_ref`, where `sigma_I_ref` is the sigma of the paired
`ImagMeas` at the same end, or `state_estimation.ia_current_floor_A`
without one; a row whose paired magnitude fails the `ImagMeas` value gate
leaves with it. Current angles share the
`ImagMeas` iteration gate and are excluded from observability; gated rows
are out of the solve, of `J` and of the redundancy.

### Shunt estimation and back-calculation (`ShuntQMeas`)

The operating point of switchable shunts (reactors, capacitor banks) is
often stale in the model. With a direct bay measurement the susceptance
becomes an estimator state (case A); without one, a bay current plus a
measured voltage yields a derived reactive-power pseudo-measurement
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

- A `ShuntQMeas` (and a bus-referenced `ImagMeas`) measures the bay of the
  bus: its prediction is the sum over every in-service shunt at that bus,
  as the power flow stamps them. With several shunts at a bus,
  `setShuntEstimation!` releases the one registered last; the others stay
  at their model values.
- A released shunt without an active direct bay measurement (`ShuntQMeas` or
  bus-referenced `ImagMeas`) stays frozen at the model value with a warning,
  and so does a `B` column no measurement pins.
- Shunts under the voltage-dependent injection mode (`bus_shunt_model`) are
  rejected.
- Case B needs the measured voltage (it enters `Q = -B V^2` squared); a
  shunt without one is skipped with a warning. The unmeasured bay `P` is
  taken as 0 and the propagated sigma doubled; filter circuits with a
  substantial active component need a direct `ShuntQMeas` (case A).
- Each released shunt consumes one degree of freedom (`ν = m - n` counts the
  `B` columns).

In case A the state vector becomes
`x = [θ(non-slack); Vm(all); B_1..B_k; α]`, initialized from the model
value; the predictions add the injection analytically,
`S_sh = |V_i|^2 conj(G + jB) baseMVA`, so the Ybus stays constant across the
FD perturbations. Case B takes `Q_sh` from `S = sqrt(3) U I`, the sign from
the model susceptance, the sigma from first-order error propagation
(`sigmaFloor` sets a minimum); `SHDERIV` rows are protected from
elimination like `ZI` rows.

### [Link flow measurements](@id se_link_measurements)

A P or Q measurement on a closed bus link (a coupler, `addPflowMeasurement!`
with `linkNr`, or a `P`/`Q` row on a link terminal of the measurement table)
is a WLS row with weight $1/\sigma^2$, like a branch flow. The link has no
Ybus entry, so the row predicts the flow from the nodal balances of the link
cluster: injection specification of each member minus its shunt draw minus
the branch flows leaving it, at the estimated state. These balances stay
equations (see [Zero injections and coupler balances](@ref se_zero_injection));
only the measurement is weighted. On a single coupler or any
tree of couplers that fixes the flow, on a ring of couplers the measured
links add one circulation state per independent loop they see. A gross
error on a link row gets a normalized residual and `wii` and is eliminated
by the sequential bad-data loop like any other row. Before its elimination
the error pulls the state in proportion to the weights. After the estimate
`calcLinkFlowsSE!` gives the link flows: the balance split at the estimated
state plus the estimated circulations. See
[State Estimation Extensions](state_estimation_extensions.md) for the
method.

### [Zero injections and coupler balances](@id se_zero_injection)

A passive bus (no load, generation or shunt) injects nothing, and so does the fused bus of a closed
coupler cluster whose members are all passive (with load or generation the fused bus carries the sum of
the member injection rows as a weighted `LINKAGG` row). A busbar joined to the network by closed
couplers only (no branch of its own, typed `Isolated` by the import) is a member of its cluster and
gets its zero-injection rows like any passive bus; without them the cluster would have no balance.
With `state_estimation.zero_injection = on`
(default, weighted) each balance is a pseudo-measurement with sigma
`state_estimation.zero_injection_sigma_mva` in MVA (default 0.001 MVA, 1 kW). It counts in $J$ and in the ranking, is never
eliminated and leaves an injection of the order of the sigma. With `constraint` the balances are
equality constraints $c(x) = 0$, solved in an extended system with Lagrange multipliers $\lambda$
(KKT). The balance then holds exactly, and the rows leave $J$ while the degrees of freedom stay $m - n$:

```math
\begin{bmatrix} H_M^T W H_M & C^T \\ C & 0 \end{bmatrix} \begin{bmatrix} \Delta x \\ \lambda \end{bmatrix} = \begin{bmatrix} H_M^T W r_M \\ -c(x) \end{bmatrix}, \qquad C = \partial c / \partial x
```

$\lambda_k$ is the force of the measurements against the balance at bus $k$; a large value means a
measurement error nearby. A small sigma makes the gain matrix ill-conditioned. The KKT system scales
the constraint rows by $a = \max \operatorname{diag}(G) / \max |C|$ with $G = H_M^T W H_M$, so both blocks are of the same
order; the solution does not depend on $a$, the condition number does. `SEResult.conditioning` gives the
condition number of the run. On sp_case14 with a noisy generated set and the passive buses as zero
injections it is about 6.9e9 for the gain matrix under `on` with the default sigma and 2.6e4 for the
KKT matrix under `constraint`. Inside a cluster the link flows $f = A^{+} b(x) + N c$ meet every member balance by
construction, and the fused balance is the `LINKAGG` row, so there is no separate link multiplier.

## Observability

Run `evaluate_global_observability(net; ...)` before estimating (structural
islands, then FD-aware numeric rank; verdict `:good` / `:critical` /
`:not_observable`), and `evaluate_local_observability(net, cols; ...)` for
placement studies on chosen state columns; see
[Observability](observability.md).

## Integration with the Net workflow

SE runs on the same `Net` as the power flow:

1. Build or import the `Net`
2. Build measurements (SCADA/PMU/custom), or `runpf!` plus
   `generateMeasurementsFromPF` for synthetic studies
3. Check observability (global/local)
4. Run the estimator (`runse!`)
5. Optionally write the estimate back (`state_estimation.update_net = true`)

**Diagnostics workflow**

| Call | Returns |
|---|---|
| `validate_measurements(net, measurements)` | one SE run; objective statistics (`value`, `dof`, the classical score `zscore`, the Wilson-Hilferty score `z_wh`, and `within_3sigma`, which is the Wilson-Hilferty band test `\|z_WH\| ≤ 3`), the ranking by `\|normalized_residual\|` with `wii` and a `localizable` flag per row, the suspicious list, optional residual-correlation columns |
| `runse_diagnostics(net, measurements)` | adds the [sequential elimination](@ref se-max-eliminations): one trace row per elimination (id, normalized residual before, `wii`, objective before/after) and a `stop_reason`: `:consistent`, `:band_low`, `:no_suspicious_left`, `:no_localizable_suspect`, `:max_eliminations` or `:not_converged`; the counter `range_check` of the [range check](@ref se-range-check) |
| `summarize_se_diagnostics(diag)` | `global_consistency` (`true` when the estimator converged and the objective passes the band test, else `false` with a `reason`), suspicious count |
| `print_se_diagnostics(diag; io=stdout, topN=10, format=:plain\|:markdown)` | statistics, ranking (with `wii`) and elimination trace |

### Write-back policy

Estimation never mutates the model unless a write-back key is set, and
then only from a converged run: `state_estimation.update_net` (keyword
`updateNet`) writes the estimated voltages, `state_estimation.update_taps`
(keyword `updateTaps`, off by default) the fixed tap positions,
`state_estimation.update_shunts` (keyword `updateShunts`, on by default)
the estimated susceptances. A tap changer with `status = not_estimated`
(its `reason` says why) and a frozen shunt keep their exact model value; without the flag the model value is bitwise untouched. The
snapshot power flow of
[PF from the estimate](state_estimation_measurements.md#PF-from-the-estimate)
removes its temporary prosumers after the run for the same reason.

### [Tap estimation switch and fallback](@id se-tap-fallback)

`state_estimation.tap_estimation.enabled` (default `true`, Web UI "estimate
taps") switches the tap estimation. Off, no changer becomes a state, even
one released with `setTapEstimation!`: every changer keeps its transferred
step (SSH profile or measurement file), and each row of `tapEstimates` and
`se_tap_estimates.csv` says `status = not_estimated`, `reason = disabled`.
The rows carry `status` (`estimated`, `out_of_range` for an estimate too
far outside the band, `not_estimated`, `calculated` for a machine
transformer) and `reason` for every changer.

**Fallback.** When the run with tap estimation does not converge (the
estimation run, the fixation run, or non-finite values), `runse!` solves
once more without tap estimation on the same measurements and start
values, and that second run is the result. run.log of the SE service gets
two lines:

    tap estimation: run 1 not converged (estimation_not_converged, 40 iterations), falling back
    run 2 without tap estimation: converged, J = 812.4, dof = 96

`tapFixation.fallback_used` and the metadata key `se_tap_fallback_used` are
set, every changer row has `reason = fallback_not_converged` (an excluded
or machine transformer keeps its own reason, see
[Tap estimation per transformer](@ref se-tap-exclude)), and the result page
shows one line above the tap table. If the second run does not converge
either, the run fails with both reasons.
The SE service falls back the same way when its staged start rejects a
converged estimate: J after the fixation above the J of the pre-solve with
the transferred steps (relative margin 1e-6, acceptance rule `worse_fit`).
The first line then reads `run 1 rejected (worse_fit, ...)` and the rows
have `reason = fallback_worse_fit`.
When the staged start with the taps frozen does not converge, the SE
service starts no released run and falls back at once (reason
`staged_not_converged`, see
[Staged start and fallback](state_estimation_extensions.md)).
`state_estimation.tap_estimation.fallback = false` (Web UI "tap fallback")
switches off both and returns the first run as it is.

**Bad data on the end steps.** With released taps the SE service runs the
bad-data diagnostics on the end steps of the result run, not on free tap
states. A round is: the result run (taps estimated and fixed, or the
fallback above), then `runse_diagnostics` with every tap frozen at its end
step (the fixed step; a tap not estimated, and every tap after the
fallback, at its transferred step). When it eliminates a row, the next
round repeats the result run with released taps on the cleaned set; the
rounds end with the first round that eliminates nothing. One row is
eliminated per round, so every row is judged on the end steps of a run that
no longer contains the rows removed before it. `max_eliminations` counts
over all rounds. A gross error that drives the released run into the
fallback is judged on the transferred steps. One that the free taps absorb
in part is judged on the steps it pulled; once it is eliminated, the
repeated run estimates the taps without it. run.log names each eliminated row with the result it was found on, the
metadata key `se_bad_data_rounds` gives the number of rounds. Without a
released tap the sequence is one round: diagnostics, rows out, final run.

### [Tap estimation per transformer](@id se-tap-exclude)

Every single transformer can be kept out of the tap estimation.
`state_estimation.tap_estimation.exclude` lists transformers by name or
CGMES mRID. Such a transformer keeps its transferred step, its rows say
`status = not_estimated`, `reason = excluded`, and its step is not
calculated either. An entry that names no transformer of the net gives a
warning naming it (`runse!`, `runse_diagnostics`; the SE service writes
it to run.log) and has no effect.

Machine (generator step-up) transformers are kept out by default: the AVR
sets the machine terminal voltage, nothing observes it independently, and
a tap state would absorb whatever the AVR does. After a converged run with
`update_net = true` their step is calculated from the AVR setpoint, the
dispatch and the measured machine P and Q at the machine bus (the method of
`calcMachineTrafoTapFromSE`): `status = calculated`,
`reason = machine_transformer`, nothing is written back. Without a result
(no setpoint, or a run without `update_net`) the row stays `not_estimated`
with that reason and its notice says why.
`state_estimation.tap_estimation.include_machine_transformers = true`
estimates machine transformers like the others.
A machine transformer has a terminal that connects to nothing but this
transformer and carries a generator. The buses that closed bus links tie
together count as one terminal: a generator behind a closed coupler is a
machine transformer's generator, and a closed link is no further
connection. The calculation then takes the setpoint and schedule of the
generators on those buses and the sum of their measured P and Q
injections. A sum is taken only when every one of these buses has an
injection row, as in the link aggregation of the estimator. Without a Q
row on every bus the step is not calculated: under voltage control the
scheduled Q is not the machine's Q. The row stays `not_estimated` with
`reason = q_not_measured_in_group` and the transferred step, and one info
line (run.log in the SE service) names the transformer and the buses
without Q. Without a P row on every bus the scheduled P is used and the
row carries `note = p_from_schedule` (column `note` of
`se_tap_estimates.csv`).
With `include_machine_transformers = true` such a
transformer is released like any other and the rank test decides. A
machine transformer is a bridge; without a voltage measurement on the
machine node it stays estimated as long as flows at both ends and the Q
of the machine tell the step apart through the losses, otherwise it is
kept out (`reason = radial_no_voltage_pin`). Both rules are decided
once on the net as given, so the bad-data diagnostics, which solve on a
copy with the closed links merged, keep out the same transformers as the
final run.

Both keys act on every released changer, also one released with
`setTapEstimation!`, and also in the run of the fallback. With the switch
off every changer says `disabled` and nothing is calculated. They are case
scope.
In the Web UI the Case page lists the transformers with tap changers with
one "estimate" checkbox each (machine transformers unchecked by default);
**Save case options** writes the unchecked ones into `exclude` and sets
`include_machine_transformers` when a machine transformer is checked. A
machine transformer the case configuration names in `exclude` stays there
on a save (its row says "excluded (not calculated)").
run.log of the SE service has one line per excluded and per machine
transformer.

### [Observability per changer](@id se-tap-observability)

Before the first iteration the estimator decides from the measurement set
which steps it can estimate: a rank test on the measurement Jacobian at
the start point, with the step as a state and every column scaled to norm
1. A changer whose column the other columns explain (for example a phase
changer on a bridge without an angle anchor: its angle and every angle
behind it move together) is not estimated. It keeps its transferred step
(SSH profile or measurement file), its row says `status = not_estimated`
with the reason (`radial_no_voltage_pin`, `not_separable`,
`no_angle_anchor`, `symmetric_start` name the structural cause,
`not_observable` that none applies), and the notice gives the measured
minus estimated voltage at its terminals as a hint. No angle measurement is
needed: a phase changer in a mesh is estimated when the flows of the mesh
fix it. After the fixation the rank test runs once more at the converged
state: a changer it drops there, or a phase changer without an angle
anchor whose mirror step fits as well, is withdrawn
(`not_observable_at_estimate`, `mirror_step_fits`), and a changer that
becomes determinable is estimated again with the other changers of its
transformer, also in a run where the start test kept every changer out;
that re-estimate is taken over only with `sigma_s` below 0.5 steps
(`weak_at_estimate` otherwise). The other changers of the transformer of a
withdrawn changer, or of one kept at its transferred step as
`out_of_range`, are estimated again with it at that step.
For a changer kept out and for one of a run that fell back, the column
`step_from_target` gives the step that meets the control target of its tap
controller (voltage setpoint, or active power over the transformer),
calculated from the estimate after a run with `update_net = true`. It is
not applied; a notice says when it differs from the transferred step by
more than one step. Every estimated step carries `sigma_s`, its standard
deviation from the covariance of the estimate, and
`round_distance_sigma`, the distance to the step it rounds to in units of
`sigma_s`. Details: [Transformer tap estimation](@ref se-tap-estimation).

### The Wilson-Hilferty band test

The objective $J = \sum_i (z_i - h_i(\hat{x}))^2 / \sigma_i^2$ says whether
measurements, sigmas and model fit together. When every sigma describes its
measurement error, `J` is a χ² variable with `ν = m - n` degrees of freedom
(`m` the active rows after the contraction, zero-injection
pseudo-measurements included, current-angle rows gated at the final state
excluded; `n` the state count: angles, magnitudes, the PMU offset and every
released shunt or tap state). `E[J] = ν`, standard deviation `√(2ν)`, hence
the rule of thumb `J ≈ ν` (the run message writes `dof` for `ν`) and the
leading number of the run message, `J/dof`.

**Verdicts**

| Verdict | Meaning |
|---|---|
| passed (`\|z_WH\| ≤ 3`) | measurements, sigmas and model agree |
| `:high` | `J` too large: the residuals exceed what the sigmas allow, bad data or a model error (a wrong parameter, tap or switch state), the classical alarm |
| `:low` | `J` implausibly small: the sigmas overstate the errors; the normal signature of noise-free synthetic data, not an alarm. Not possible at `ν = 1` and `2`: there the lower bound of the band is below zero |
| `:no_redundancy` | `ν = 0`, the test is skipped (the residuals are structurally zero) |
| `small_redundancy` flag | `ν < 30`, informational: reduced power, relative spread `√(2/ν)` |

**Notes**

- Read `J/ν`, not `J`: a bare `J` grows with the number of measurements
  (`J = 104` at `dof = 95` is a ratio of 1.1 and normal).
- The classical score `z = (J - ν)/√(2ν)` is still reported; the
  Wilson-Hilferty verdict decides.
- `runse_diagnostics` eliminates only while the band test fails.
- On multi-island nets every island carries its own verdict next to the
  summed one (see [Island-wise estimation](@ref)).
- An exhausted elimination that still leaves `:high` with the suspects
  clustered at one station points at a topology error (see
  [Topology validation](@ref)).

The textbook band `|z| ≤ 3` on `z = (J - ν) / √(2ν)` assumes a normal `J`.
`J` is χ²-distributed and skewed to the right, so the symmetric band flags
too many large `J` and passes too many small ones, and below `ν = 18` its
lower bound goes negative. Wilson and Hilferty showed that the cube root
`(J/ν)^{1/3}` is very nearly normal, with mean `1 - 2/(9ν)` and variance
`2/(9ν)`, usable from `ν ≈ 3` upward; the band test applies the three-sigma
rule to that cube root:

```math
z_{WH} = \frac{(J/\nu)^{1/3} - \left(1 - \tfrac{2}{9\nu}\right)}{\sqrt{2/(9\nu)}}, \qquad |z_{WH}| \le 3
```

Translated back to `J` the band is wider above `ν` and narrower below it.

Against exact χ² quantiles the Wilson-Hilferty upper bound is within 3 %
even at `ν = 1` (false-alarm rate 0.12 % instead of 0.135 %). The lower
bound is noticeably too low below `ν ≈ 10` and below zero for `ν ≤ 2`,
which affects only the `:low` verdict. The Delta-J notice of the tap
fixation does not use this approximation: it compares the J rise with the
exact 95 % quantile of χ² with `k` degrees of freedom, `k` the step
states the fixation removes.

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
| Elimination | a row is suspicious from `rn >= k_eliminate` on; the elimination workflow removes it, bounded by the budget `max_eliminations`, only while the band test fails | `state_estimation.k_eliminate` (default 3.0), `state_estimation.max_eliminations` (default 5) | identification of gross errors |
| `off` | plain WLS | `state_estimation.robust_mode = off` | clean data |
| `staged` | two-stage R modification, active from `state_estimation.robust_start_iteration` (default 3) on, with knees `robust_k1`/`robust_k2` (defaults 3/6) on the raw ratio `t = \|z - h(x)\|/σ`: `t ≤ k1` leaves the weight unchanged, `k1 < t ≤ k2` widens tangentially (`σ_mod = σ (2t/k1 - 1)`, continuous at `t = k1`), `t > k2` suppresses the gradient contribution (`σ_mod = \|z - h\|/k1`); at `t = k2` the sigma steps down from `σ (2 k2/k1 - 1)` to `σ k2/k1` (3σ to 2σ with 3/6), as in the formulas before the knees were configurable; stages are recomputed every iteration, so a measurement can recover | `state_estimation.robust_mode = staged`; the Boolean `state_estimation.robust` is an alias for `staged` and applies while `robust_mode` is `off` | the classic behavior |
| `replacement` | pins every row whose normalized residual reaches `k_suppress` to the fixed `suppression_sigma` (in the measurement's unit) during the solve; the row stays in the system and in every statistic | `state_estimation.robust_mode = replacement`, `state_estimation.k_suppress` (default 4.0), `state_estimation.suppression_sigma` (default 2000) | online smoothing |

Defaults and ranges of every key:
[State-Estimation Configuration](state_estimation_configuration.md).

**Result fields**

| | |
|---|---|
| Affected rows | `SEResult.robustRows` (stage 1/2 for staged, stage 3 for replacement, with `t` and `σ_mod/σ`) |
| Active objective | `SEResult.activeObjective` (`J_active`; also the Web UI summary, `run.log`, metadata keys `se_objective_active` / `se_dof_active` / `se_suppressed_rows`) |

**Notes**

- Only localizable rows are touched: a row whose `wii` does not exceed 0.3
  (the `localizable` flag) is neither down-weighted nor eliminated, however
  large its `rn`, because intervening would push the error onto the
  neighbour that shares its redundancy. Non-localizable suspects are
  reported and passed over (`skipped_unlocalizable` per round in the
  trace); if only such suspects remain the workflow stops with
  `:no_localizable_suspect`, a measurement gap rather than bad data. A
  band failure `:low` (J below the band, the sigmas overstate the errors)
  is not bad data either: the elimination does not start and reports
  `:band_low`. A row whose $\Omega_{ii}$ sits at the numerical floor gets
  `rn = 0`.
- Virtual measurements (`σ ≤ 1e-6`), the `ZI` rows (sigma
  `state_estimation.zero_injection_sigma_mva`, default 1e-3 MW, also the
  aggregate of a coupler cluster of them) and the derived shunt rows
  (`SHDERIV`) are exempt in every mode: they are never eliminated and never
  down-weighted by `replacement` or `staged`. They carry network structure
  or other measurements, not a meter of their own.
- The replacement threshold is judged on the converged state, never on
  flat-start transients: the estimator solves plain, builds the suppression
  set from the converged residuals, re-solves with those weights frozen and
  repeats until the set is stable.
- The modification changes only the solve weights: `J`, normalized
  residuals, `wii`, `K` and the ranking run on the original sigmas, so a
  suppression cannot hide the measurement it suppresses. `J_active` is the
  objective over the rows the estimator trusted, with its own dof; the band
  verdict stays on the honest `J`.
- `suppression_sigma` is the sigma a down-weighted row is solved with, not a
  third limit. The normalized residual is the right scale because for a
  nearly critical row $\Omega_{ii}$ approaches zero while the raw ratio
  $|r_i|/\sigma_i$ stays small (at `wii = 0.3` the residual spread is
  $\sqrt{\Omega_{ii}} = 0.55\,\sigma$).

### [Range check of the voltage measurements](@id se-range-check)

Before every WLS (`runse!`, `validate_measurements`, `runse_diagnostics`,
hence also the SE service and the Web UI) each active voltage measurement
outside 0.5 to 1.5 pu, 0 kV and non-finite values included, is checked
against the switching state. On a bus that is not energized (an isolated
bus, or an island without a source) 0 kV is right: the bus and its branches
stay out of the estimation, the rows of the bus, its branches and its links
are inactive for the run, reason `deenergized`. On an energized bus the
value is a gross error. It is not eliminated but replaced by a pseudo
measurement at the imported state of the bus (the `SvVoltage` of a CGMES
import, otherwise 1.0 pu) with its sigma times
`state_estimation.bad_data.substitute_sigma_factor` (default 3), reason
`range_check_substituted`. Each finding gets one line in the log and in
`run.log`; the diagnostics carry the counter `range_check` (`substituted`,
`deenergized`) and mark the row in the ranking (field
`range_check = :range_check_substituted`), the SE service writes the reason
into `measurements_verdict.csv` and the counts into the metadata
(`se_range_check_substituted`, `se_range_check_deenergized`). Only then run
the WLS, the residual analysis and the elimination.

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

Protected rows (`ZI`, `SHDERIV`, `LINKAGG` aggregates and any row with
`sigma <= 1e-6`) encode model knowledge, not telemetry, and are never
eliminated. The trace lands in `se_diagnostics.md`. Configuration key
`state_estimation.max_eliminations` (default 5, one row per round, so also
the round limit of the loop); the Web UI field **Sequential elimination
budget** sets it per run. Pseudo measurements of the
[range check](@ref se-range-check) are never eliminated.

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

With active `ImagMeas` rows the sensitivities `W` and correlations `K` are
load-flow dependent: valid for the estimated operating point, not network
constants. Below `takahashi_min_states` states, for K-matrix requests (the
correlations need the full m×m `Ω`) and after a Takahashi guard failure
(with a warning) `Ω_ii` comes from the dense `pinv` path; above 2000 states
that fallback is a blockwise sparse solve.

### What runs at which size

Which path the residual diagnostics take, for Vm plus P/Q injection sets,
with the default `report_residual_correlation = true`. Up to 2000 states
the correlation report needs the full `Ω` and pins the dense path; above
2000 states `Ω_ii` comes from the Takahashi path and the correlations from
a blockwise solve:

| Network | m × n | residual diagnostics | notes |
|---|---|---|---|
| sp_case60 (60 buses) | 326 × 117 | dense path (below the 200-state Takahashi threshold) | every diagnostic available |
| sp_case188 (188 buses) | 562 × 373 | dense path (the correlation report needs the full `Ω`; Takahashi with `report_residual_correlation = false`) | every diagnostic available |
| case1354pegase | 4062 × 2707 | Takahashi | every diagnostic available; criticality from `diag(Omega)` (no size budget) |
| case13659pegase | 40977 × 27317 | Takahashi | criticality from `diag(Omega)`; K-matrix report skipped above 20000 rows |

Three bounds:

* the per-row criticality tests (`criticality_method = rank`) are skipped
  with a warning above an m·n budget of 300000, the dark-state SVD behind
  `unobservable_state_columns` above 2000 states or measurements;
* the dense `pinv` diagnostics path stops at 2000 states; above it the
  diagonal comes from the Takahashi path or, when that refuses, from a
  blockwise sparse solve;
* the `state_estimation.report_residual_correlation` columns are skipped
  with a warning above 20000 measurement rows.

Large assemblies perturb states that share no measurement row together
(column coloring), bit-identical to the per-column assembly. The
diagnostics Jacobian carries the released shunts' `B` states at their
estimated values, and `n` in `ν = m - n` includes them.

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

The helpers resolve bus names and branch references; the `Measurement`
constructor is the hand-built equivalent (last row).

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
ok || error("Validation failed: $msg")

empty!(net.measurements)
addVmMeasurement!(net; busName = "B1", value = 1.01, sigma = 0.002)
addVmMeasurement!(net; busName = "B2", value = 0.99, sigma = 0.004)
addPinjMeasurement!(net; busName = "B2", value = -25.0, sigma = 1.0)
addQinjMeasurement!(net; busName = "B2", value = -8.0, sigma = 1.0)
addPflowMeasurement!(net; fromBus = "B1", toBus = "B2", value = 24.0, sigma = 0.8, direction = :from)
addQflowMeasurement!(net; branchNr = 1, value = 6.5, sigma = 0.8, direction = :to)
push!(net.measurements, Measurement(typ = PflowMeas, value = 23.5, sigma = 0.8, branchIdx = 1, direction = :from, id = "PF_12_REDUNDANT"))

# observability check and runse! as in the minimal example
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

Literature on PMU-based state estimation: A. Abur, A. Gómez Expósito,
*Power System State Estimation: Theory and Implementation*; A. G. Phadke,
J. S. Thorp, *Synchronized Phasor Measurements and Their Applications*;
NASPI TR-006, *Phase Angle Calculations: Considerations and Use Cases*.

## Further examples and workshop material

* 7-bus tutorial: [advanced workshop tour](generated/workshop_tour_advanced.md), [state-estimation notebook](generated/workshop_state_estimation.md)
* Bad-data localization, elimination, robust estimation, shunt estimation: the [SE diagnostics notebook](generated/workshop_se_diagnostics.md)
* Tap estimation, current measurements, PMU phasors, topology advisory, islands, $J$ versus $J_{active}$: the [taps, phasors and topology notebook](generated/workshop_se_taps.md)
* Scripts in `examples/state_estimation/`: `state_estimation_wls.jl`, `state_estimation_pmu_angles.jl`, `state_estimation_observability.jl`, `state_estimation_passive_bus_zib_comparison.jl`, `state_estimation_zero_injection_constraint.jl`, `h_matrix_observability_demo.jl`
* A [Sparlectra Case Format](scf.md) case carries its own measurement set, see [Case files carry their own measurements](@ref se-case-file-measurements)

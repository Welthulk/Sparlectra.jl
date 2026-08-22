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
* Angle residuals are not wrap-corrected; with realistic transmission-grid
  angle spreads (≪ 180°) this is not a limitation.
* PMU magnitude measurements need no special treatment — model them as
  `VmMeas` with a small sigma (e.g. 0.002 p.u.).

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

SE quality depends strongly on observability.

### Global observability

Use `evaluate_global_observability(net; ...)` to assess if the complete
state can be estimated from the active measurements stored in `net.measurements`.

Typical metrics:

* Number of measurements `m`
* Number of states `n`
* Redundancy `r = m - n`
* Redundancy ratio `ρ = m / n`
* Numerical/structural observability flags
* Quality label (e.g. `:observable`, `:critical`, `:not_observable`)

### Local observability

Use `evaluate_local_observability(net, cols; ...)` to assess a selected subset
of state columns (for example one bus angle and one bus magnitude).

This is useful for sensor-placement studies and for identifying vulnerable areas.

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

* `validate_measurements(net, measurements; normalizedThreshold=3.0, ...)`
  runs SE once and returns:
  * objective statistics (`value`, `dof`, `zscore`, `within_3sigma`)
  * largest normalized residual
  * full measurement ranking by `|normalized_residual|`
  * suspicious measurement list based on `normalizedThreshold`
* `runse_diagnostics(net, measurements; deactivate_and_rerun=true, ...)`
  extends this with a deactivate-and-rerun step for the top suspicious
  measurement. Note: a single rerun can improve the objective significantly,
  but `global_consistency` may still stay `false` when multiple bad or
  mismodeled measurements remain active.
* `summarize_se_diagnostics(diag)` creates a compact interpretation summary
  (`global_consistency`, `reason`, suspicious count).
* `print_se_diagnostics(diag; io=stdout, topN=10, format=:plain|:markdown)`
  pretty-prints the statistics and ranking for reports.

Meaning of `global_consistency`:

* `true`: estimator converged and the objective is within the 3σ plausibility
  band (globally plausible data/model fit).
* `false`: either non-convergence or implausibly large objective value
  (possible bad data, wrong uncertainty model, or topology/model mismatch).

This workflow can be used both for automated checks (NamedTuple result
inspection) and for human-readable diagnostics output.

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

## Further examples and workshop material

* Extended tutorial and a simple 7-bus setup: [workshop tour](generated/workshop_tour.md)
* Detailed WLS reporting example script: `examples/state_estimation/state_estimation_wls.jl`
* PMU angle measurements and the reference-offset state α: `examples/state_estimation/state_estimation_pmu_angles.jl`
* Observability-focused scenario script: `examples/state_estimation/state_estimation_observability.jl`
* Passive-bus ZIB comparison example: `examples/state_estimation/state_estimation_passive_bus_zib_comparison.jl`
* Matrix-based observability/redundancy demo: `examples/state_estimation/h_matrix_observability_demo.jl`

## H-matrix observability demo (A..E)

If you want to study observability directly on Jacobian-like matrices `H` without
building a full network first, use:

```bash
julia --project=. examples/state_estimation/h_matrix_observability_demo.jl
```

The script evaluates each matrix with:

* Structural observability via sparsity matching
* Numerical observability via rank test
* Per-row redundancy classification (critical vs. redundant)
* Local observability on selected state-column subsets

Included demo matrices:

* `H_A`: fully observable with duplicate information (clear redundancy).
* `H_B`: minimal square observable case (`m = n`), therefore every row is critical.
* `H_C`: structurally observable but numerically fragile (near linear dependence).
* `H_D`: sparse case highlighting matching behavior and extra measurements.
* `H_E`: incidence-like matrix paired with a toy graph/spanning-tree interpretation.

This is intended as a compact didactic companion to
`evaluate_global_observability` / `evaluate_local_observability`.

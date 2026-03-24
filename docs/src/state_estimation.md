State Estimation
=================

This page summarizes the state-estimation (SE) functionality in Sparlectra and
shows how it connects to regular network studies.

> **Release status:** State Estimation is currently **experimental**. The
> current implementation is intended as a first practical WLS workflow for
> studies, examples, and early application feedback.

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

* `VmMeas` (bus voltage magnitude)
* `PinjMeas`, `QinjMeas` (bus injections)
* `PflowMeas`, `QflowMeas` (branch flows with direction)

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
  measurement.
* `summarize_se_diagnostics(diag)` creates a compact interpretation summary
  (`global_consistency`, `reason`, suspicious count).
* `print_se_diagnostics(diag; io=stdout, topN=10, format=:plain|:markdown)`
  pretty-prints the statistics and ranking for reports.

This workflow can be used both for automated checks (NamedTuple result
inspection) and for human-readable diagnostics output.

## Minimal example

```julia
using Sparlectra
using Random

net = run_acpflow(casefile = "case9.m")

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
addBus!(net = net, busName = "B1", busType = "Slack", vn_kV = 110.0)
addBus!(net = net, busName = "B2", busType = "PQ", vn_kV = 110.0)
addBus!(net = net, busName = "B3", busType = "PQ", vn_kV = 110.0)
addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B3", toBus = "B1", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0)

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
addBus!(net = net, busName = "B1", busType = "Slack", vn_kV = 110.0)
addBus!(net = net, busName = "B2", busType = "PQ", vn_kV = 110.0)
addBus!(net = net, busName = "B3", busType = "PQ", vn_kV = 110.0)
addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B3", toBus = "B1", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0)

empty!(net.measurements)
addVmMeasurement!(net; busName = "B1", value = 1.01, sigma = 0.002)
addPinjMeasurement!(net; busName = "B2", value = -25.0, sigma = 1.0)
addQinjMeasurement!(net; busName = "B2", value = -8.0, sigma = 1.0)
addPflowMeasurement!(net; fromBus = "B1", toBus = "B2", value = 24.0, sigma = 0.8, direction = :from)
addQflowMeasurement!(net; branchNr = 1, value = 6.5, sigma = 0.8, direction = :to)

obs = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
println("Observable quality: ", obs.quality)
```

## Further examples and workshop material

* Extended tutorial and a simple 7-bus setup: [Workshop](workshop.md)
* Detailed WLS reporting example script: `src/examples/state_estimation_wls.jl`
* Observability-focused scenario script: `src/examples/state_estimation_observability.jl`
* Passive-bus ZIB comparison example: `src/examples/state_estimation_passive_bus_zib_comparison.jl`
* Matrix-based observability/redundancy demo: `src/examples/h_matrix_observability_demo.jl`

## H-matrix observability demo (A..E)

If you want to study observability directly on Jacobian-like matrices `H` without
building a full network first, use:

```bash
julia --project=. src/examples/h_matrix_observability_demo.jl
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

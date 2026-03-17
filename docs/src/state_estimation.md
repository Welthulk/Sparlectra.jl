State Estimation
=================

This page summarizes the state-estimation (SE) functionality in Sparlectra and
shows how it connects to regular network studies.

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

## Measurement model

The current implementation supports these measurement types:

* `VmMeas` (bus voltage magnitude)
* `PinjMeas`, `QinjMeas` (bus injections)
* `PflowMeas`, `QflowMeas` (branch flows with direction)

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

Use `evaluate_global_observability(net, meas; ...)` to assess if the complete
state can be estimated from active measurements.

Typical metrics:

* Number of measurements `m`
* Number of states `n`
* Redundancy `r = m - n`
* Redundancy ratio `ρ = m / n`
* Numerical/structural observability flags
* Quality label (e.g. `:observable`, `:critical`, `:not_observable`)

### Local observability

Use `evaluate_local_observability(net, meas, cols; ...)` to assess a selected
subset of state columns (for example one bus angle and one bus magnitude).

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

## Minimal example

```julia
using Sparlectra
using Random

net = run_acpflow(casefile = "case9.m")

std = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7)
meas = generateMeasurementsFromPF(
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

gobs = evaluate_global_observability(net, meas; flatstart = true, jacEps = 1e-6)
println("Global observability quality: ", gobs.quality)

se = runse!(net, meas; maxIte = 12, tol = 1e-6, flatstart = true, jacEps = 1e-6, updateNet = true)
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
ok || error("Validation failed: $msg")

meas = Measurement[
    Measurement(typ = VmMeas, value = 1.01, sigma = 0.002, busIdx = 1, id = "VM_B1"),
    Measurement(typ = VmMeas, value = 0.99, sigma = 0.004, busIdx = 2, id = "VM_B2"),
    Measurement(typ = PinjMeas, value = -25.0, sigma = 1.0, busIdx = 2, id = "PINJ_B2"),
    Measurement(typ = QinjMeas, value = -8.0, sigma = 1.0, busIdx = 2, id = "QINJ_B2"),
    Measurement(typ = PflowMeas, value = 24.0, sigma = 0.8, branchIdx = 1, direction = :from, id = "PF_12"),
    Measurement(typ = PflowMeas, value = 23.5, sigma = 0.8, branchIdx = 1, direction = :from, id = "PF_12_REDUNDANT"),
]

obs = evaluate_global_observability(net, meas; flatstart = true, jacEps = 1e-6)
println("Observable quality: ", obs.quality)

se = runse!(net, meas; maxIte = 12, tol = 1e-6, flatstart = true, jacEps = 1e-6, updateNet = true)
println("Converged: ", se.converged)
```

## Further examples and workshop material

* Extended tutorial and a simple 7-bus setup: [Workshop](workshop.md)
* Detailed WLS reporting example script: `src/examples/state_estimation_wls.jl`
* Observability-focused scenario script: `src/examples/state_estimation_observability.jl`

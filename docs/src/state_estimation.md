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

Typical workflow for generating measurements:

1. Solve a power flow to get a physically consistent reference state.
2. Create synthetic measurements using `generateMeasurementsFromPF`.
3. Configure standard deviations via `measurementStdDevs`.
4. Optionally add Gaussian noise.

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
2. Run AC power flow (`runpf!`)
3. Build measurements (`generateMeasurementsFromPF` or custom)
4. Check observability (global/local)
5. Run estimator (`runse!`)
6. Optionally write estimates back into the network (`updateNet = true`)

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

## Further examples and workshop material

* Extended tutorial and a simple 7-bus setup: [Workshop](workshop.md)
* Detailed WLS reporting example script: `src/examples/state_estimation_wls.jl`
* Observability-focused scenario script: `src/examples/state_estimation_observability.jl`

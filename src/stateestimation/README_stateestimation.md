# State estimation (developer notes)

## Purpose

Classical WLS state estimation: the measurement model, the estimator with
observability and bad-data diagnostics, transformer tap estimation as
extra states, and topology validation.

## Include order

```julia
include("stateestimation/measurements.jl")
include("stateestimation/tap_estimation.jl")     # after numerics/takahashi.jl
include("stateestimation/state_estimation.jl")
include("stateestimation/topology_validation.jl")
```

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `measurements.jl` | `MeasurementType`, `Measurement`, all add helpers, synthetic generation, CSV IO | `addVmMeasurement!`, `addPmuPhasorMeasurement!`, `generateMeasurementsFromPF` |
| `state_estimation.jl` | The WLS estimator, observability, bad-data diagnostics, result views | `runse!`, `runse_diagnostics`, `se_view` |
| `tap_estimation.jl` | Transformer tap ratios as extra estimator states | `setTapEstimation!`, `calcMachineTrafoTapFromSE` |
| `topology_validation.jl` | Station-level topology hypothesis testing | `validate_topology`, `test_topology_hypotheses` |

## Conventions

- Measurement standard deviations and sigma floors come from
  `measurementStdDevs`/`measurementSigmaFloors`; synthetic measurements
  generated from a power flow are marked as ideal.
- PMU angle measurements may carry a time base that does not line up with the
  slack reference; the estimator solves for the offset as an extra state
  (`state_estimation.pmu_ref_offset`).
- The J/dof statistic leads the diagnostics presentation; keep new
  diagnostics consistent with `runse_diagnostics` output.

## Reference

Rendered API documentation: `docs/src/reference_stateestimation.md`.
User-facing documentation: `docs/src/state_estimation.md`.

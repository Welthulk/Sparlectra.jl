# Examples Overview

This page summarizes the most relevant runnable examples in `examples/`.

## Power flow and network operation

- `using_links.jl`  
  Minimal link workflow (open/close links) and resulting PF behavior.
- `exp_configuration.jl`
  Demonstrates central configuration loading and typed runtime configuration usage.
- `exp_synthetic_tiled_grid_pf_perf.jl`
  Synthetic tiled-grid PF performance example.

## Transformer and tap control

- `tap_control_demo_grid.jl`  
  Lightweight three-controller demo (OLTC + PST + Schrägregler) using `run_acpflow(net = ...)`,
  central Sparlectra configuration (`examples/configuration.yaml` or
  `SPARLECTRA_CONFIGURATION_YAML`) plus demo-specific
  `examples/tap_control_demo_grid.yaml`, and `latest_control_result(net)` for
  controller rows and trace rows. Optional classic output:
  `SPARLECTRA_TAP_DEMO_CLASSIC=1`; optional raw control-result rows:
  `SPARLECTRA_TAP_DEMO_RAW=1`.

## Voltage-dependent and Q-limit controls

- `example_voltage_dependent_control_rectangular.jl`  
  Demonstrates P(U)/Q(U) behavior with the rectangular solver workflow.
- `example_q_limit_voltage_adjustment.jl`  
  Shows `qlimit_mode = :adjust_vset` and related run variants.

## Reporting and exports

- `using_netreports.jl`  
  Builds and prints structured `ACPFlowReport` outputs.
- `export_solution.jl`  
  Exports solved network data for downstream usage and writes outputs under
  `examples/_out/export_solution/...`.

## State estimation

- `state_estimation_wls.jl`  
  Baseline weighted least squares state estimation run.
- `state_estimation_manual_measurements.jl`  
  Manual measurement setup and estimation workflow.
- `state_estimation_observability.jl`  
  Observability-focused scenario and diagnostics.
- `state_estimation_passive_bus_zib_comparison.jl`  
  Passive-bus / ZIB handling comparison.
- `usage_state_estimation_diagnostics.jl`  
  Practical diagnostics usage workflow.
- `h_matrix_observability_demo.jl`  
  Matrix-level observability/redundancy exploration.

## Utility / analysis scripts

- `matpower_import.jl`
  Configurable MATPOWER import utility script using the central
  `src/configuration.yaml.example` defaults plus optional
  `examples/configuration.yaml` local overrides.
  Supports logfile-only MATPOWER reference diagnostics, effective configuration logging, and
  configurable branch `SHIFT` sign/unit handling for phase-shifter convention checks.
  The YAML options can also diagnose and choose PV/REF voltage references when
  MATPOWER `BUS.VM` values differ from generator `GEN.VG` setpoints; use
  `compare_voltage_reference: imported_setpoint` to validate against the actual
  imported PV/REF setpoints, or `compare_voltage_reference: bus_vm` to reproduce
  stored bus-matrix voltages strictly.
- `network_analyzer.jl`  
  Developer/exploration script for network diagnostics and inspection helpers.
- `visul_chaos.jl`  
  Developer/exploration stress/visual helper script.

## How to run

From repository root:

```bash
julia --project=. examples/tap_control_demo_grid.jl
```

General pattern:

```bash
julia --project=. examples/<example_name>.jl
```

# Examples Overview

This page summarizes the most relevant runnable examples in `src/examples/`.

## Power flow and network operation

- `using_links.jl`  
  Minimal link workflow (open/close links) and resulting PF behavior.
- `compare_stands_loadflow.jl`  
  Compares Sparlectra load-flow outputs against reference expectations.
- `import_case_xxxyyyy_loadflow.jl`  
  End-to-end import and load-flow execution flow.
- `exp_autodamped_rectangular_pf.jl`: Demonstrates `autodamp = true` for
  residual-based rectangular Newton backtracking.
- `exp_start_projection_rectangular_pf.jl`: Demonstrates `start_projection = true`
  with DC-angle and blend-scan candidates.

## Transformer and tap control

- `tap_control_demo_grid.jl` (+ `tap_control_demo_grid.yaml.example`)  
  Configurable demo grid for multi-controller experiments and logging.

## Voltage-dependent and Q-limit controls

- `example_voltage_dependent_control_rectangular.jl`  
  Demonstrates P(U)/Q(U) behavior with the rectangular solver workflow.
- `example_q_limit_voltage_adjustment.jl`  
  Shows `qlimit_mode = :adjust_vset` and related run variants.

## Reporting and exports

- `using_netreports.jl`  
  Builds and prints structured `ACPFlowReport` outputs.
- `export_solution.jl`  
  Exports solved network data for downstream usage.

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

- `matpower_import.jl` (+ `matpower_import.yaml.example`)  
  Configurable MATPOWER import utility script.
  Supports logfile-only MATPOWER reference diagnostics, effective YAML logging, and
  configurable branch `SHIFT` sign/unit handling for phase-shifter convention checks.
  The YAML options can also diagnose and choose PV/REF voltage references when
  MATPOWER `BUS.VM` values differ from generator `GEN.VG` setpoints; use
  `compare_voltage_reference: imported_setpoint` to validate against the actual
  imported PV/REF setpoints, or `compare_voltage_reference: bus_vm` to reproduce
  stored bus-matrix voltages strictly.
- `network_analyzer.jl`  
  Network diagnostics and inspection helpers.
- `visul_chaos.jl`  
  Stress/visual exploration helper script.

## How to run

From repository root:

```bash
julia --project=. src/examples/example_transformer_tap.jl
```

General pattern:

```bash
julia --project=. src/examples/<example_name>.jl
```

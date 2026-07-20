# Examples Overview

This page summarizes the most relevant runnable examples in `examples/`.

## Power flow and network operation

- `using_links.jl`  
  Minimal link workflow (open/close links) and resulting PF behavior.
- `exp_synthetic_tiled_grid_pf_perf.jl`
  Synthetic tiled-grid PF performance example.
- `exp_configured_matpower_cases.jl`
  Runs ordered `matpower_import.cases` entries sequentially through `run_sparlectra_cases`.
- `exp_transformer_loss_extension.jl`
  Demonstrates Sparlectra's MATPOWER transformer-loss extension by exporting a
  small transformer case with active no-load conductance metadata and reimporting
  it without double-counting the equivalent bus shunts.
- `exp_programmatic_api.jl`
  Runs one MATPOWER case through the GUI-ready `run_sparlectra_api` contract and lists explicit artifacts.
- `exp_current_iteration_start.jl`
  Demonstrates enabling the guarded current-iteration start pre-solve through normal API configuration overrides and prints its metadata/artifact status.
- `exp_powerflow_service.jl`
  Starts a local service run, looks up its serialized result by run ID, and lists its artifacts without an HTTP server.

## Transformer and tap control

- `exp_transformer_tap_changer_model.jl`
  Imports the same off-nominal-tap MATPOWER transformer once with
  `transformer.tap_changer_model = :ideal` and once with
  `:impedance_correction`, and compares the resulting series impedance and
  power-flow solution.
- `exp_3wt_phase_taps.jl`
  Builds the same three-winding transformer with `create3WTWindings!` in
  three separate configurations, demonstrating the `phase_tap_side`/
  `phase_taps` keywords (Issue #261): **Case 1** OLTC only (ratio tap on the
  HV winding); **Case 2** PST/Schrägregler only (asymmetrical
  `PhaseTapChangerModel` on the MV winding, no ratio tap anywhere); **Case 3**
  (bonus) both combined on the HV winding (`tap_side == phase_tap_side`).
  Data-model only — resolving a phase-tap model into the AUX-bus branch and
  outer-loop control of a single 3WT winding are not wired yet.
- `tap_control_demo_grid.jl`  
  Lightweight three-controller demo (OLTC + PST + Schrägregler) using `run_sparlectra(net = ...)`,
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

- `export_solution.jl`  
  Exports solved network data for downstream usage and writes outputs under
  `examples/_out/export_solution/...`.

## DTF validation (Testnetz13 / FOR001-FOR002)

External FOR001/FOR002 validation datasets are not shipped with Sparlectra;
place local files under `data/DTF/` or pass explicit `--dtf-file`/`--for002-file`
paths. See the tests page for the full mode/case reference.

- `validate_dtf_suite.jl`  
  Unified CLI runner for all DTF checks: import audit, native base-case
  validation against FOR002, DTF-listed outage validation, and the
  DTF -> existing MATPOWER export/import roundtrip. Select checks with
  `--mode=all|audit|base|outages|matpower` and datasets with `--case=A,B,...`.
  The former single-purpose scripts `validate_dtf_for002_testnetz13.jl`,
  `validate_dtf_for002_outages_testnetz13.jl`, and
  `validate_dtf_matpower_export_testnetz13.jl` were consolidated into this
  suite; the shared implementations live in `examples/internal/dtf_validation_*.jl`,
  each of which is also directly runnable as a single-purpose CLI entry point.
- `dtf_validation_report.jl`  
  Cross-case report over the local FOR001/FOR002 case set: transformer loss
  decomposition, voltage-transfer diagnostics, and transformer-ratio-mode
  comparisons, written as CSV/Markdown under `examples/_out/dtf_validation/`.
- `for002_matpower_metadata_validation.jl`  
  Standalone diagnostic for FOR002/MATPOWER metadata fixtures (bus/branch
  name normalization and comparison artifacts).

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
- `matpower_import_multi_config.jl`  
  Developer diagnostic that runs one MATPOWER case through several
  configuration YAML files and compares the final rectangular PF status
  (auto-profile mode, wrong-branch detection, solver start settings).
- `qlimit_large_case_mode_comparison.jl`  
  Manual Q-limit comparison workflow over optional large MATPOWER cases
  (`case × start_profile × qlimit_mode`); excluded from normal test profiles.
- `network_analyzer.jl`  
  Developer/exploration script for network diagnostics and inspection helpers.

## How to run

From repository root:

```bash
julia --project=. examples/tap_control_demo_grid.jl
```

General pattern:

```bash
julia --project=. examples/<example_name>.jl
```

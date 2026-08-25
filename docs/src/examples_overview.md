# Examples Overview

All runnable examples live under `examples/`, grouped by topic:
`examples/powerflow/`, `examples/others/`, `examples/state_estimation/`,
`examples/dtf/`, shared infrastructure in `examples/internal/`, experimental
material in `examples/experimental/`. Run any example directly:

```bash
julia --project=. examples/<folder>/<example>.jl
```

or run a whole topic through its suite runner (fresh subprocess per example,
summary at the end): `run_powerflow_suite.jl`, `run_others_suite.jl`,
`run_state_estimation_suite.jl`, `run_val_dtf_suite.jl`,
`run_cgmes_suite.jl`, `run_short_circuit_suite.jl`,
`run_parallel_suite.jl` (parallel-vs-serial demos with identity checks).

## Power flow and network operation

| Example | Folder | Demonstrates | Suite |
|---|---|---|---|
| `matpower_import.jl` | `powerflow` | Minimal MATPOWER import and solve | powerflow |
| `matpower_import_multi_config.jl` | `powerflow` | One case under several configurations | powerflow |
| `exp_configured_matpower_cases.jl` | `powerflow` | Ordered `matpower_import.cases` batches via `run_sparlectra_cases` | powerflow |
| `exp_programmatic_api.jl` | `powerflow` | GUI-ready `run_sparlectra_api` contract, explicit artifacts | powerflow |
| `exp_powerflow_service.jl` | `powerflow` | Local service run, result lookup by run ID, artifact listing (no HTTP server) | powerflow |
| `exp_distributed_slack_modes.jl` | `powerflow` | Classical single slack vs distributed slack (`pg_weighted` vs imported `APF` shares, `lambda_P` metadata) | standalone |
| `exp_external_grid_comparison.jl` | `powerflow` | Grid-connection modeling comparison on an 8-bus ring (issue #299): ideal slack at two alternative buses, non-ideal external-grid source with internal impedance, and distributed slack, tabulated per bus | powerflow |
| `exp_dc_powerflow.jl` | `powerflow` | Standalone DC power flow (`rundcpf!`), DC-seeded AC start | powerflow |
| `exp_condition_number.jl` | `powerflow` | Jacobian condition-number estimate at the solved operating point (`condestJacobian`, `reportCondition`), incl. a stressed near-singular variant | powerflow |
| `exp_hvdc_b2b_pairing.jl` | `others` | Back-to-back HVDC pairing controller: Stage-0 snapshot vs steerable transfer on a two-island net (#297 Draft B) | others |
| `exp_hvdc_meshed_ac_tie.jl` | `others` | HVDC pair in parallel to an AC tie: one reference per synchronous island, setpoint pair as parallel PQ path, `island_feed` rejected ([theory](hvdc_back_to_back.md)) | others |
| `exp_parallel_islands.jl` | `powerflow` | `islands.mode solve_parallel` vs `solve_independent` on an 8-island net, wall clocks side by side, bitwise-identical voltages | parallel |
| `exp_parallel_sc_sweep.jl` | `others` | IEC 60909-0 all-bus sweep serial vs threaded chunks (8000 fault locations), row-identical results | parallel |
| `exp_contingency_n1.jl` | `others` | Full branch N-1 on case1354pegase (`runContingencies!`): serial vs parallel, warm-start evidence, top-10 worst contingencies ([theory](contingency.md)) | parallel |
| `exp_open_terminal_line.jl` | `others` | One-sided open line: full charging draw at the closed bus and the Ferranti rise at the open end ([theory](branchmodel.md)) | others |
| `exp_current_iteration_start.jl` | `powerflow` | Guarded current-iteration start pre-solve via config overrides | powerflow |
| `exp_diagnose_self_check.jl` | `others` | `run_fixed_reference_self_check` and the narrative `diagnose.log` report | others |
| `exp_short_circuit.jl` | `others` | `runShortCircuit!` (IEC 60909-0) on a hand-built feeder+machine net: Ik'' max/min and the safety flag on defaulted data | short_circuit |
| `exp_short_circuit_reference.jl` | `others` | PASS/FAIL check against the analytic IEC 60909-0 reference values | short_circuit |
| `exp_short_circuit_cgmes.jl` | `others` | `runShortCircuit!` on the ENTSO-E MicroGrid BE delivery from the local test-set cache | short_circuit |
| `exp_synthetic_tiled_grid_pf_perf.jl` | `powerflow` | Synthetic tiled-grid PF performance study | powerflow |
| `qlimit_large_case_mode_comparison.jl` | `powerflow` | Q-limit enforcement modes on large cases | powerflow |
| `apslf_demo.jl` | `powerflow` | APSLF analytic solver as standalone/primary/start-value backend | powerflow |
| `mc_probabilistic_powerflow.jl` | `powerflow` | Monte-Carlo load scaling on `case14.m` (N = 1000): per-bus Vm statistics, band violations, convergence rate | powerflow |

## Transformers and controllers

| Example (`examples/others/`) | Demonstrates | Suite |
|---|---|---|
| `exp_transformer_tap_changer_model.jl` | `tap_changer_model = :ideal` vs `:impedance_correction` on an off-nominal tap | others |
| `exp_transformer_loss_extension.jl` | MATPOWER transformer-loss extension export/reimport round trip | others |
| `exp_3wt_phase_taps.jl` | 3WT with `phase_tap_side`/`phase_taps`: OLTC-only, PST-only, combined (data model only) | others |
| `tap_control_demo_grid.jl` | Three controllers at once (OLTC + PST + combined) via `run_sparlectra(net = …)` and `latest_control_result` | others |
| `tap_control_schraeg_two_controllers.jl` | Split combined regulation: independent voltage/ratio and power/phase controllers on one unit | others |
| `exp_pst_reactance_coupling.jl` | PST control with tap-dependent series reactance X(α) vs. a static-reactance run | others |
| `exp_svc_shunt_voltage_control.jl` | SVC variable-shunt voltage control (in-range + honest at_limit) with the controllable-element view | others |
| `exp_tcsc_series_reactance_control.jl` | TCSC series-reactance control steering a loop-network flow split onto a branch target ([theory](series_compensation.md)), incl. honest `at_limit` | others |
| `exp_facts_limit_modes.jl` | FACTS limit characteristics side by side ([theory](facts.md)): constant-Q box vs STATCOM `V*S_max` vs SVC `V^2*B`, plus SSSC injected-voltage window vs TCSC fixed window (#297 Drafts A/E/F) | others |
| `exp_auto_slack_selection.jl` | Automatic slack selection (`power_flow.auto_slack` / `ensureSlack!`) when a case registers no voltage reference | others |
| `exp_ac_rescue_dc_fallback.jl` | Non-convergence handling: `power_flow.rescue` strategy ladder plus `power_flow.dc.fallback` standalone-DC result | others |
| `exp_cgmes_import_analysis.jl` | `analyzeCGMES` report naming the missing declared dependency of an incomplete CGMES delivery | others |
| `exp_cgmes_infer_base_voltages.jl` | `cgmes_import.infer_base_voltages`: reconstruct missing nominal voltages from the SV state and solve | others |
| `exp_cgmes_topology_processor.jl` | Node-breaker import without a TP profile (#314): derived bus partition compared bus for bus against the shipped TP on MiniGrid | others |
| `machine_remote_voltage_control.jl` | Remote voltage control via machine reactive power ([theory](remote_voltage_control.md)), incl. the honest `at_limit` outcome | others |
| `using_links.jl` | Busbar coupler as bus link, open/close behavior | others |
| `network_analyzer.jl` | Topology analysis before/after removing a branch | others |
| `export_solution.jl` | Solver-agnostic `PFModel`/`PFSolution` export | others |

## Voltage-dependent and Q-limit control

| Example (`examples/powerflow/`) | Demonstrates | Suite |
|---|---|---|
| `example_voltage_dependent_control_rectangular.jl` | P(U)/Q(U) droop behavior in the rectangular solver | powerflow |
| `example_q_limit_voltage_adjustment.jl` | `qlimit_mode = :adjust_vset` run variants | powerflow |

## CGMES

| Example | Demonstrates | Suite |
|---|---|---|
| `run_cgmes_suite.jl` | Guided walkthrough (diagnose → import → solve → SV validation) plus the full sweep over every ENTSO-E/ReliCapGrid test set with a result table | cgmes |
| `experimental/cgmes_fetch_testsets.jl` | One-time fetch of the ENTSO-E test-set package into the local cache | standalone |
| `experimental/cgmes_export_demo.jl` | CGMES export (`writeCGMESFiles`) on a small net | others (optional) |
| `experimental/val_realgrid_remeasure.jl` | RealGrid measurement ladder: baseline, Q-limits, distributed slack, machine control | standalone |

## DTF validation (Testnetz13 / FOR001-FOR002)

External FOR001/FOR002 datasets are not shipped; place files under
`data/DTF/` or pass `--dtf-file`/`--for002-file`.

| Example | Demonstrates | Suite |
|---|---|---|
| `run_val_dtf_suite.jl` | Unified CLI: import audit, base-case validation, outage validation, MATPOWER round trip (`--mode`, `--case`) | val_dtf |
| `dtf/dtf_validation_report.jl` | Cross-case CSV/Markdown report: loss decomposition, voltage-transfer diagnostics, ratio-mode comparisons | standalone |
| `dtf/for002_matpower_metadata_validation.jl` | FOR002/MATPOWER metadata fixtures (name normalization, comparison artifacts) | standalone |

## State estimation

| Example (`examples/state_estimation/`) | Demonstrates | Suite |
|---|---|---|
| `state_estimation_wls.jl` | Baseline WLS run | state_estimation |
| `state_estimation_manual_measurements.jl` | Manual measurement setup | state_estimation |
| `state_estimation_observability.jl` | Observability scenario and diagnostics | state_estimation |
| `state_estimation_passive_bus_zib_comparison.jl` | Passive-bus / ZIB handling comparison | state_estimation |
| `state_estimation_pmu_angles.jl` | PMU angle measurements and the reference-offset state α | state_estimation |
| `usage_state_estimation_diagnostics.jl` | Practical diagnostics workflow | state_estimation |
| `h_matrix_observability_demo.jl` | Matrix-level observability/redundancy exploration | state_estimation |
| `mc_state_estimation_study.jl` | Monte-Carlo WLS error statistics on the 7-bus [workshop-tour network](generated/workshop_tour.md) (M = 500) | state_estimation |

The **suite** column names the runner that executes the example
(`standalone` = run directly, not part of a suite; `optional` = skipped when
its inputs are unavailable). The suite runners print a per-example pass/fail
summary and are themselves exercised by the extended test profile.

# Examples Overview

Runnable examples live under `examples/powerflow/`, `examples/others/`
(also the shared suite infrastructure), `examples/state_estimation/`,
`examples/dtf/` and `examples/cgmes/`. Run one directly:

```bash
julia --project=. examples/<folder>/<example>.jl
julia --project=app examples/powerflow/exp_programmatic_api.jl   # examples of the service layer: the app environment
```

or a whole topic through its suite runner (fresh subprocess per example,
summary at the end): `run_powerflow_suite.jl`, `run_others_suite.jl`,
`run_state_estimation_suite.jl`, `run_val_dtf_suite.jl`,
`run_cgmes_suite.jl`, `run_short_circuit_suite.jl` and
`run_parallel_suite.jl` (parallel vs serial with identity checks). The
**Suite** column names the runner (`standalone`: run directly;
`optional`: skipped when its inputs are unavailable); the runners run in
the extended test profile.

## Power flow and network operation

| Example | Folder | Demonstrates | Suite |
|---|---|---|---|
| `matpower_import.jl` | `powerflow` | MATPOWER import and solve | powerflow |
| `matpower_import_multi_config.jl` | `powerflow` | one case under several configurations | powerflow |
| `exp_configured_matpower_cases.jl` | `powerflow` | `runtime.cases` batches via `run_sparlectra_cases` | powerflow |
| `exp_programmatic_api.jl` | `powerflow` | the `run_sparlectra_api` contract | powerflow |
| `exp_powerflow_service.jl` | `powerflow` | local service run, result and artifact lookup by run ID | powerflow |
| `exp_distributed_slack_modes.jl` | `powerflow` | single vs distributed slack (`pg_weighted`, imported `APF` shares) | standalone |
| `exp_external_grid_comparison.jl` | `powerflow` | ideal slack, external-grid source and distributed slack on an 8-bus ring | powerflow |
| `exp_dc_powerflow.jl` | `powerflow` | standalone DC power flow (`rundcpf!`), DC-seeded AC start | powerflow |
| `exp_condition_number.jl` | `powerflow` | Jacobian condition estimate (`condestJacobian`, `reportCondition`) | powerflow |
| `exp_hvdc_b2b_pairing.jl` | `others` | HVDC pairing controller on a two-island net ([theory](hvdc_back_to_back.md)) | others |
| `exp_hvdc_meshed_ac_tie.jl` | `others` | HVDC pair in parallel to an AC tie, `island_feed` rejected ([theory](hvdc_back_to_back.md)) | others |
| `exp_parallel_islands.jl` | `powerflow` | islands on threads vs serial, bitwise-identical voltages | parallel |
| `exp_parallel_sc_sweep.jl` | `others` | all-bus short-circuit sweep serial vs threaded, row-identical | parallel |
| `exp_contingency_n1.jl` | `others` | full branch N-1 on the shipped `sp_case1354.m`, serial vs parallel ([theory](contingency.md)) | parallel |
| `exp_open_terminal_line.jl` | `others` | one-sided open line: charging draw and Ferranti rise ([theory](branchmodel.md)) | others |
| `exp_current_iteration_start.jl` | `powerflow` | guarded current-iteration start pre-solve via config overrides | powerflow |
| `exp_diagnose_self_check.jl` | `others` | `run_fixed_reference_self_check` and the `diagnose.log` report | others |
| `exp_short_circuit.jl` | `others` | `runShortCircuit!` on a hand-built net, safety flag on defaulted data | short_circuit |
| `exp_short_circuit_reference.jl` | `others` | PASS/FAIL check against analytic IEC 60909-0 values | short_circuit |
| `exp_short_circuit_cgmes.jl` | `others` | `runShortCircuit!` on the ENTSO-E MicroGrid BE delivery (local test-set cache) | short_circuit |
| `exp_synthetic_tiled_grid_pf_perf.jl` | `powerflow` | synthetic tiled-grid performance study | powerflow |
| `qlimit_large_case_mode_comparison.jl` | `powerflow` | Q-limit enforcement modes on large cases | powerflow |
| `apslf_demo.jl` | `powerflow` | APSLF as standalone, primary or start-value backend | powerflow |
| `mc_probabilistic_powerflow.jl` | `powerflow` | Monte-Carlo load scaling on the shipped `sp_case14` (N = 1000) | powerflow |

## Transformers and controllers

| Example (`examples/others/`) | Demonstrates | Suite |
|---|---|---|
| `exp_transformer_tap_changer_model.jl` | `tap_changer_model = :ideal` vs `:impedance_correction` on an off-nominal tap | others |
| `exp_transformer_loss_extension.jl` | MATPOWER transformer-loss extension round trip | others |
| `exp_3wt_phase_taps.jl` | 3WT with `phase_tap_side`/`phase_taps`: OLTC-only, PST-only, combined (data model only) | others |
| `tap_control_demo_grid.jl` | three controllers at once (OLTC, PST, combined), `latest_control_result` | others |
| `tap_control_schraeg_two_controllers.jl` | split combined regulation: two independent controllers on one unit | others |
| `exp_pst_reactance_coupling.jl` | PST control with tap-dependent series reactance X(α) vs a static reactance | others |
| `exp_svc_shunt_voltage_control.jl` | SVC variable-shunt voltage control, in range and `at_limit` | others |
| `exp_tcsc_series_reactance_control.jl` | TCSC series-reactance flow control ([theory](series_compensation.md)) | others |
| `exp_facts_limit_modes.jl` | FACTS limit characteristics side by side ([theory](facts.md)) | others |
| `exp_facts_base_impedance.jl` | physical base impedance after a full-UPFC run ([theory](facts.md)) | others |
| `exp_auto_slack_selection.jl` | automatic slack selection (`power_flow.auto_slack` / `ensureSlack!`) | others |
| `exp_matpower_gencost.jl` | MATPOWER generator costs on the shipped `sp_case9`, written back unchanged (also via SCF) | others |
| `exp_island_reference_priority.jl` | reference priority on the shipped `data/scf/two_islands_prio.scf.json`, in the power flow and in N-1 | others |
| `exp_powsybl_iidm_import.jl` | PowSyBl IIDM import (`read_iidm_tables`) checked against the OpenLoadFlow reference ([theory](powsybl_import.md)) | others |
| `exp_jacobian_block.jl` | the 2x2 Jacobian block between two buses (`jacobianBlock`), with a network diagram | others |
| `exp_tap_influence_zone.jl` | one tap moved by a few steps: the buses whose voltage changes beyond a threshold | others |
| `exp_tap_sensitivity.jl` | bus-voltage sensitivity to one tap from the Jacobian (`dV/dtau = -J^-1 dF/dtau`), with a Markdown report | others |
| `exp_ac_rescue_dc_fallback.jl` | the `power_flow.rescue` ladder plus the `power_flow.dc.fallback` result | others |
| `exp_cgmes_import_analysis.jl` | `analyzeCGMES` on an incomplete delivery | others |
| `exp_cgmes_infer_base_voltages.jl` | `cgmes_import.infer_base_voltages` from the SV state | others |
| `exp_cgmes_topology_processor.jl` | node-breaker import without a TP profile (MiniGrid) | others |
| `machine_remote_voltage_control.jl` | remote voltage control via machine reactive power ([theory](remote_voltage_control.md)) | others |
| `using_links.jl` | busbar coupler as bus link, open/close behavior | others |
| `network_analyzer.jl` | topology analysis before and after removing a branch | others |
| `export_solution.jl` | solver-agnostic `PFModel`/`PFSolution` export | others |

## Voltage-dependent and Q-limit control

| Example (`examples/powerflow/`) | Demonstrates | Suite |
|---|---|---|
| `example_voltage_dependent_control_rectangular.jl` | P(U)/Q(U) droop behavior in the rectangular solver | powerflow |
| `example_q_limit_voltage_adjustment.jl` | `qlimit_mode = :adjust_vset` run variants | powerflow |
| `example_qlimit_reenable_voltage_rule.jl` | PQ->PV release of the active-set mode on the Zeng/Chiang 14-bus case | powerflow |

## CGMES

| Example | Demonstrates | Suite |
|---|---|---|
| `run_cgmes_suite.jl` | guided walkthrough plus the sweep over every ENTSO-E/ReliCapGrid test set | cgmes |
| `cgmes/cgmes_fetch_testsets.jl` | one-time fetch of the ENTSO-E test-set package into the local cache | standalone |
| `cgmes/cgmes_export_demo.jl` | CGMES export (`writeCGMESFiles`) on a small net | others (optional) |
| `cgmes/val_realgrid_remeasure.jl` | RealGrid measurement ladder: baseline, Q-limits, distributed slack, machine control | standalone |

## DTF validation (Testnetz13 / FOR001-FOR002)

External FOR001/FOR002 datasets are not shipped; place files under
`data/DTF/` or pass `--dtf-file`/`--for002-file`.

| Example | Demonstrates | Suite |
|---|---|---|
| `run_val_dtf_suite.jl` | unified CLI: audit, base-case and outage validation, MATPOWER round trip (`--mode`, `--case`) | val_dtf |
| `dtf/dtf_validation_report.jl` | cross-case CSV/Markdown validation report | standalone |
| `dtf/for002_matpower_metadata_validation.jl` | FOR002/MATPOWER metadata fixtures | standalone |

## State estimation

| Example (`examples/state_estimation/`) | Demonstrates | Suite |
|---|---|---|
| `state_estimation_wls.jl` | baseline WLS run | state_estimation |
| `state_estimation_manual_measurements.jl` | manual measurement setup | state_estimation |
| `state_estimation_observability.jl` | observability scenario and diagnostics | state_estimation |
| `state_estimation_passive_bus_zib_comparison.jl` | passive-bus / ZIB handling comparison | state_estimation |
| `state_estimation_pmu_angles.jl` | PMU angle measurements and the reference-offset state α | state_estimation |
| `usage_state_estimation_diagnostics.jl` | diagnostics workflow | state_estimation |
| `state_estimation_imag_bad_data.jl` | `ImagMeas` raising localizability (`wii`), elimination trace | state_estimation |
| `state_estimation_shunt_estimation.jl` | shunt estimation as a state (case A) and as a `SHDERIV` pseudo-measurement (case B) | state_estimation |
| `state_estimation_links_facts.jl` | SE on the contracted net (links): `LINKAGG`, `se_view` | state_estimation |
| `state_estimation_robust.jl` | robust R modification suppressing a gross error | state_estimation |
| `state_estimation_chain.jl` | measurement CSV v1 round trip and the SE-started PF (`runpf_from_se!`) | state_estimation |
| `bench_takahashi_diagnostics.jl` | benchmark: dense `pinv` vs Takahashi selected inverse (`takahashi_min_states`) | (manual, not in suite) |
| `h_matrix_observability_demo.jl` | matrix-level observability/redundancy exploration | state_estimation |
| `mc_state_estimation_study.jl` | Monte-Carlo WLS error statistics on the 7-bus [workshop-tour network](generated/workshop_tour.md) (M = 500) | state_estimation |

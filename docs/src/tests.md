# Test Suite

Sparlectra uses profile-aware test loading through `test/runtests.jl`.
Profile selection precedence is:
1. CLI argument (`julia --project=. test/runtests.jl <profile>`)
2. `SPARLECTRA_TEST_PROFILE`
3. default `fast`

## Test profiles

| Profile | Command | Scope | Intended use |
|---|---|---|---|
| `fast` (default) | `julia --project=. test/runtests.jl fast` | Small unit tests and focused integration regressions | Normal development and default CI |
| `extended` | `julia --project=. test/runtests.jl extended` | Extended-only integration, lifecycle, documentation, repository-hygiene, example, and stress tests | Before merge or after broad integration changes |
| `all` | `julia --project=. test/runtests.jl all` | Fast followed by extended | Complete local or CI verification |

## Fast versus extended ownership

| Area | Fast coverage | Extended coverage |
|---|---|---|
| Core/model | Small constructors, transformer checks, bus/prosumer semantics, link behavior, representative rectangular PF/Q-limit regressions, and small sparse fallbacks. Large sparse Ybus and large MATPOWER matrix-body checks are extended-only. | Large sparse Ybus smoke, large MATPOWER matrix-body scanner, synthetic/stress grids, and longer integration examples. |
| API | Serialization and transport helpers, path and validation safety, one successful small API run, one pre-solver failure, one numerical/island failure, Solver-time and Total-time contracts, critical DC-line default and strict-rejection smokes, and one small independent-island regression. | Exhaustive CSV/export matrices, repeated artifact inventories, island artifact-content matrices, persistent history/delete/reload lifecycle coverage, and repeated presentation/performance-log modes. |
| Web UI | Form parsing and backend validation, result rendering, active and terminal timing cards, commit-span omission, tolerance-step hook, path traversal rejection, DTF upload role classification, primary-case and FOR002 selector filtering, and stubbed route checks without a real asynchronous solver run. | Case-profile persistence with asynchronous jobs, real run/result polling, artifact preview/download/ZIP/history/delete lifecycle checks, browser-launcher matrices, socket/server lifecycle, and Markdown/help/documentation cross-products. |
| Documentation/hygiene | No repository-wide documentation/help scan in fast; only focused source-level smoke checks tied to edited paths. | Configuration documentation consistency and normalized tracked-path/content repository hygiene scans. |

`Pkg.test()` uses the same test runner and therefore the default `fast` profile unless `SPARLECTRA_TEST_PROFILE` is set:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Fast profile groups

| Group | Files | Main checks |
|---|---|---|
| `core_model` | `test/testgrid.jl` | Core net construction and validation, small inline MATPOWER import helpers, file-based MATPOWER projected-start normalization, PV/PQ lock-ID resolution, Q-limit enforcement-mode parsing and classical base-failure/no-reenable dispatch checks, rectangular nonfinite mismatch and status-row diagnostic preservation, Jacobian condition-number estimator checks (exact diagonal case, Hager bound, near-singular flagging, input rejection), link handling, shunts, reporting/output checks, and summary-file output regression. Large sparse Ybus checks and large MATPOWER matrix-body scanner coverage are extended-only. |
| `terminal_status` | `test/test_terminal_status.jl` | Per-terminal branch status (r0.9.10, one-sided open branches): the dangling-node reference anchor (a `:open_to`/`:open_from` reduction matches a solve with an explicit zero-injection auxiliary bus to 1e-10, for a line and for a transformer with off-nominal ratio and phase shift), bitwise equivalence of both-flags-open with `status = 0` in the Y-bus, solver toggles between runs (sparsity-pattern invalidation), isolation of the open-side bus in the island report, the MATPOWER export substitution (BR_STATUS 0 plus exact `Y_in` bus shunt, `open_terminal=` marker, voltage-identical reimport), the result surface (`open@to` marker, `Open terminals` header count, `terminal_state`/open-end voltage report columns), and the state-estimation guard rejecting flow measurements on partially open branches |
| `configuration` | `test/test_configuration_coverage.jl`, `test/test_runner_helpers.jl` | Configuration-key coverage, safe refresh checks including current-iteration start defaults, forwarding checks, start-voltage value-domain validation (including `cgmes_import.start_values` and `cgmes_import.placeholder_guards` parsing: defaults `flat`/`warn_skip`, `sv`/`strict` accepted, invalid values rejected with the key name), and test-runner output-mode helper checks |
| `matpower_metadata` | `test/test_matpower_metadata.jl` | MATPOWER parser metadata fields, legacy bus sorting of `bus_name`, opt-in bus-name import, branch-kind override, branch metadata retention, FOR001 contingency mapping, default/opt-in `mpc.dcline` behavior including PF/loss/Q/V/Q-limit mapping and voltage-control safeguards, `writeMatpowerCasefile` `write_solution` column count/VM-VA/marker behavior, the `mpc.sparlectra.tap_changer_model` roundtrip marker preventing double tap-impedance correction, and `net.name` staying the case name (not a bus's original name) when `mpc.bus_name` metadata is present |
| `programmatic_api` | `test/test_api.jl` | Focused GUI-ready power-flow API coverage for serialization and transport helpers, path and validation safety (including run-index acceptance of symlinked output roots, skipped on Windows), one successful small API run, one pre-solver failure, one numerical/island failure, Solver-time and Total-time contracts, critical DC-line default and strict-rejection smokes, and one small independent-island regression. Exhaustive CSV/export matrices, full artifact inventories, persistent restart/history/delete lifecycle matrices, repeated performance-log modes, complete island artifact-content matrices, and the fixed-reference self-check contract (verbatim-start summary in `self_check.log`, per-bus `self_check_residuals.csv`, pinned case14 residual regression) are extended-only in `programmatic_api_extended`. |
| `webui` | `test/test_webui.jl` | Focused Web UI coverage for form parsing and backend validation, result rendering, active and terminal timing cards, commit-span omission, tolerance-step hook, path traversal rejection, DTF upload role classification, primary-case and FOR002 selector filtering, the CGMES-export run option (form checkbox + help topic, `export_cgmes` request flag incl. hidden-false/absent-field defaults, and the result-page summary row for completed/failed/absent export metadata), the last-edit-wins precedence between the configuration file and saved case settings (a newer YAML wins for its own keys and sets the notice flag, an older one keeps the sidecar values), and stubbed route checks without a real asynchronous solver run. Real asynchronous PowerFlow job lifecycles, artifact preview/download/ZIP/history/delete matrices, browser-launcher platform matrices, socket/server lifecycle checks, Markdown help/documentation cross-product validation, and repeated real MATPOWER runs are extended-only in `webui_extended`. |
| `state_estimation` | `test/test_state_estimation.jl` | WLS state-estimation behavior and observability-oriented regressions |
| `dc_powerflow` | `test/test_dc_powerflow.jl` | Standalone DC power flow (`rundcpf!`, `power_flow.solver=:dc`): reference-value comparisons against an independently-computed script (3-bus hand-verified fixture, MATPOWER case9/case14) for bus angles and branch flows, phase-shifter (`Pfinj`) coverage via a synthetic phase-shifting-transformer fixture, lossless power-balance property checks, two-island independent-solve coverage, `power_flow.dc.*`/`power_flow.solver` config validation, active-outer-loop-controller rejection (mirrors the `apslf` rejection test), `run_sparlectra` dispatch (`method=:dc`, `diagnostics.solver=:dc`), `angle_reference_deg` exact-uniform-shift behavior, the `dc_pf_status`/`rectangular_pf_status` registry separation, and `rundcpf!(net; seed_ac_start=true)` DC-seeded AC Newton-Raphson convergence |
| `distributed_slack` | `test/test_distributed_slack.jl` | Distributed active-power slack: disabled-path bit-identity to the classical solver, per-bus residual conservation (`alpha .* lambda_P` pattern, REF/load buses untouched), REF-equivalence (`explicit` weight fully on the reference bus reproduces the classical solution and `lambda_P` equals the classically REF-absorbed power), mode agreement (`pg_weighted` vs matching `explicit`), MATPOWER `APF` column import (`:imported` shares proportional to APF; APF presence never changes disabled results), fallback behavior (`error` throws, `ref_only` warns and matches classical), Q-limit interaction (PV→PQ switch keeps participation), the config-load regression for the `weights: {}` placeholder plus block-style weight tables (the minimal YAML reader has no flow-mapping support), and independent per-island `lambda_P` via the AC-island solve path. The CGMES `normalPF` arrival assertions live in the extended `cgmes_importer` group (MicroGrid zero-factor mapping, MiniGrid positive-factor arrival, cache-gated with explicit skip message). |
| `island_diagnostics` | `test/test_island_diagnostics.jl` | AC island diagnostics reporting regressions (source: a 6209-bus CGMES delivery with 158 single-bus islands): the island failure message reports the failing island's own iteration count and stage (never the `iterations=0`/`stage=before_nr` fallback while a per-island record exists), the failing island is looked up in `:ac_island_solver_statuses` when the combined status is untagged, never-attempted islands render as `not_attempted` with `unavailable` mismatch fields instead of inheriting the failed island's statistics, per-island `solver.log`/`mismatch_history.csv` artifacts are written only for attempted islands, `q_limit_processing_status` no longer aliases `failure_reason`, and the single-island combined-status fallback keeps working; also the no-slack diagnostics (all-isolated vs. isolated-slack vs. never-registered messages) and the `power_flow.auto_slack` promotion (keyword and config path, ENI-over-generator ranking, no-op on a registered slack) |
| `short_circuit` | `test/test_short_circuit.jl` | IEC 60909-0 short circuit (`runShortCircuit!`): hand-derived analytic reference cases with the full derivation in the test comments (feeder reproduces its declared current at the connection point with the c-factor cancelling; series line impedance added on top; synchronous-machine x''d conversion incl. §6.6.3 fictitious resistance; kappa/i_p with the method-b cap; line charging proven irrelevant), IEC Table-1 c-factor selection and the scalar override, the safety-flag contract as a three-way assertion (documented default + log warning + per-row flag), asynchronous-machine skip with island-wide Ik''max lower-bound flag (max case only), `:no_source`/`:isolated` row statuses, `short_circuit.c_factor` config coverage (default, valid override, rejected out-of-band value), and the WebUI data-check gating helper (`_webui_case_has_short_circuit_data`: marker found / not found / cached / unresolvable path stays enabled); the Phase 3 parallel all-bus sweep identity guard (serial vs parallel rows `isequal`-identical for max/min case with `max_tasks=1` and auto on a three-island fixture and on bundled case14 with `buses = :all`; RAN/SKIPPED stated per run) |
| `parallel_foundation` | `test/test_parallel_foundation.jl` | Thread-safety foundation and parallel island path of the multi-core work (Phases 1 and 2): `runtime.parallel.*` configuration defaults and validation (`enabled`, `max_tasks` auto/integer-string, `min_work_items`), the `parallel_max_tasks` resolver, the per-worker performance-profile helpers (`_perf_profile_child` excludes the orchestrator-only `:phase_callback`, `_perf_profile_merge!` sums timings under unchanged phase names, appends iteration rows, prefixes per-worker scalars, and swallows `nothing` children/parents), the solver status as a Net field (registry globals gone, AC/DC fields separate, `deepcopy` carries a working condest thunk), UMFPACK `copy(F)` finalizer safety, the startup summary `parallel:` line, the Web UI `runtime_parallel_enabled` checkbox wiring, and the `solve_parallel` island identity guard: bitwise-equal voltages and iteration counts vs `solve_independent` with `max_tasks=1` and auto on a five-island fixture, identical phase-name sets plus `parallel_wall_time`, per-island prefixed scalars, and the parallel failure semantics (all islands report their real status, no skips). Threaded assertions state RAN/fallback depending on `Threads.nthreads()`; the `--threads=4` battery exercises the parallel branch. |
| `external_grid` | `test/test_external_grid.jl` | External-grid element (issue #299): the external-grid/distributed-slack mutual-exclusion configuration error (both cover the same imbalance; YAML pair rejected, single option loads), the duck-typing contract guard (`NativeShortCircuitData` field-identical to `CGMESShortCircuitData` by name, order, and type — machine-checked, not comment-checked), PF invariance (an `addExternalGrid!` slack with finite `Sk''` reproduces the manual slack path to machine precision — carried SC data changes no power-flow result, incl. the full ENI record contract), SC hand calculation on a single 110 kV bus (`Zq`, `Ik''`, `Sk''`, kappa/i\_p over two R/X values, c cancelling), min-case semantics (`sk_min` used unflagged with the rx\_min→rx\_max default, explicit `rx_min` wins, missing `sk_min` skips the feeder with the engine's `:no_source` flag), parallel-feeder stacking with unique, removal-safe mrids (suffix continues after the highest surviving id, incl. the same-prefix-bus guard), input validation (`ArgumentError` table), the net-copy regression (`sc_sources` survives `deepcopy`, SC on the copy matches), the `convertSlackToExternalGrid!` demotion (self-consistent bus types right after the demotion, conversion note contents, no-slack/non-slack-bus `ArgumentError`s), the connection statement in the classical result print (`Grid connection:` prose line naming slack vs. source incl. `Sk''`/R-X and the internal slack bus; type column `SOURCE` instead of `SLACK` on the internal bus, also in the structured report rows), and the `internal_impedance` variant (slack moves to the tagged auxiliary bus, stiff `Sk''→∞` limit reproduces the ideal solution below 1e-8 pu, realistic `Sk''` droops voltage and shifts angles, the feeder record stays anchored at the physical connection bus, one internal-impedance grid per bus) |
| `controls` | `test/test_voltage_dependent_control.jl`, `test/test_transformer_phase_shift.jl`, `test/test_tap_controller.jl`, `test/test_series_reactance_control.jl`, `test/test_hvdc_pair_control.jl`, `test/test_tap_changer_model.jl`, `test/test_phase_tap_changer_model.jl`, `test/test_phase_tap_table.jl` | Voltage-dependent controls, transformer phase-shift control, tap-controller behavior, YAML controller instantiation (#305: round-trip via the real YAML reader against the programmatic twin, load-time and apply-time validation, idempotent re-apply), machine remote voltage control (`MachineVoltageControl`: API validation, secant convergence onto a remote target, honest `at_limit` on the reactive bounds), successful baseline PF preservation when controls are disabled, the `AbstractTapChangerModel`/`convention` field plus `calcRatioTapCorrection`/`calcRatioTapRange` consolidation, the CGMES `PhaseTapChangerModel`/`calcPhaseTapFraction`/`calcPhaseTapAngleRatio`/`calcPhaseTapReactance` formulas plus the DTF importer migration onto them, and the tabular `TapTablePoint`/`kind = :tabular`/`calcPhaseTapTable` override path including the formula-vs-table round-trip regression; the PST X(α) coupling (winding-resolved model lookup, continuous-angle formula evaluation, nearest-row tabular mapping, accepted moves update `x_pu` to the model value at the final angle, devices without reactance data stay static, probe direction consistent with and without tracking incl. angle/reactance restore) and the cross-type warning when a tap controller already regulates a machine controller's target bus; the SVC shunt voltage controller (API validation, in-range secant regulation onto the setpoint with the actuated shunt carrying the same susceptance, honest constant-B `at_limit` on an unreachable target, `clearShuntControllers!`) and the generic controllable-element view (`controllableElements` rows for tap, machine, and shunt controllers incl. `ControlRunResult.elements`); the TCSC series-reactance controller (`SeriesReactanceControl`: loop-network target tracking within the deadband, honest `at_limit` clamping on an unreachable target, bit-identical baseline with a disabled controller, element-row vocabulary `:series_x_pu`/`:branch_active_power`, registration validation incl. transformer rejection and the `eps_z` impedance guard, and the no-move deadband case); the HVDC pair controller (`HvdcPairControl`, #297 Draft B: registration validation incl. slack and PV guards, exact pairing invariant on a two-island fixture, reversed transfer, rating clamp with honest `at_limit`, per-side voltage-target secant incl. unreachable-target clamping, bit-identical baseline with a disabled controller, MATPOWER `paired_control` mode with the `pf_injections` consistency anchor, YAML `hvdc_pair` round trip, element-row vocabulary `:hvdc_p_transfer_mw`/`:hvdc_transfer`, the grid-forming `mode = :island_feed` incl. mode-dependent validation, the island-draw mirror, the rating clamp, and the YAML `mode: island_feed` declaration; r0.9.9 additions: the persistent `HvdcLink` records (`net.hvdcLinks` fills for MATPOWER Stage-0/paired and `addHvdcLink!`, controller attach/detach round trip), the `HVDC Link Flows` result surface (table, header count, converter-loss line, `ACPFlowReport.hvdc_links`), and the one-reference-per-synchronous-island rule (multi-reference error naming the buses, setpoint pair as parallel PQ path next to an AC tie, `invalid_topology` for a demoted `island_feed` reference)) |

## Extended profile additions

The `extended` profile is extended-only: it does not run the fast profile first. Use `all` when a single invocation must execute both fast and extended suites.

Current extended-only groups are:

- `core_model_extended`
- `powerflow_rectangular`
- `factorized_linear_solver`
- `3wt_phase_taps`
- `programmatic_api_extended`
- `webui_extended`
- `repository_hygiene`
- `apslf`
- `cgmes_importer`
- `cgmes_export`
- `dtf_extended`

| Extended addition | File | Main checks |
|---|---|---|
| `legacy/remove` | `test/testremove.jl` | Remove/delete behavior and consistency after structural edits |
| `powerflow_rectangular` | `test/test_solver_interface.jl` | Rectangular power-flow API behavior, sparse-only solver path, AC-island detection/reference validation/independent solving including aggregate all-island convergence status and iteration accounting, Q-limit and typed configuration entry checks, current-iteration start rejection diagnostics on a tiny synthetic case, including accepted and rejected framework control-status composition, ordered local MATPOWER batch execution, the thin `run_acpflow` alias, legacy-keyword rejection, the optional merit-function Armijo line search (value/scaling, acceptance/fallback/active-set-skip reasons, config validation, merit.enabled=false regression, and a converging flatstart/Q-limit-switching functional case), the optional scaled-Newton trust-region step control (config validation, acceptance/collapse unit coverage, trust_region.enabled=false regression, a converging functional case with log/diagnostics assertions, and radius-collapse non-convergence reporting), and wrong-branch detection output visibility across ACPFlowReport metadata, the AC island diagnostics CSV, and the console summary line, including the highest-voltage-level scope regression (a dip below the top level never triggers; the same dip on the top level does), and the AC rescue ladder / DC fallback (`power_flow.rescue` recovers a poisoned start and records the winning strategy, `power_flow.dc.fallback` leaves a DC state with the AC status honestly non-converged, defaults leave both mechanisms off) |
| `factorized_linear_solver` | `test/test_factorized_linear_solver.jl` | umfpack_reuse linear-solver selection and factorization-reuse behavior for the rectangular power-flow Jacobian path, including rejection of the removed klu value |
| `3wt_phase_taps` | `test/test_3wt_phase_taps.jl` | `create3WTWindings!`'s `phase_tap_side`/`phase_taps` keywords — existing-behaviour snapshot, single-winding attachment, ratio+phase coexisting on one winding, and validation errors |
| `pv_voltage_residuals` | `test/test_pv_voltage_residuals.jl` | PV-voltage residual behavior, angle-preserving voltage-setpoint starts, phase-shifted PV integration coverage, and related solver diagnostics |
| `matpower_examples` | `test/test_matpower_example.jl` | MATPOWER example runner path, output routing, performance/profile rendering, and runtime configuration forwarding, removed start-voltage alias rejection, and canonical profile-blend parsing |
| `synthetic_grids` | `test/test_synthetic_grids.jl` | Synthetic network generation and larger synthetic-grid regression coverage |
| `configuration_docs` | `test/test_configuration_docs.jl` | Configuration documentation and docs/config consistency checks |
| `cgmes_importer` | `test/test_cgmes_importer.jl` | CGMES importer: generic reader semantics on a synthetic in-memory fixture (profiles, overlays, references, DifferenceModel skip), the base-voltage inference (`infer_base_voltages` reconstructs stripped catalogs from the SV state, one summary warning, net still solves; without the option the typed abort stays) and the self-loop line guard (both terminals on one topological node map to a shunt notice, not a branch), the import-failure analysis (`importFailureAnalysis` names supplied models, missing vs. satisfied `md:Model.DependentOn` prerequisites, and the verdict; the `import_analysis_mode` service run succeeds on an importable delivery, fails with `import_analysis_not_importable` on a missing declared dependency with the `import_analysis.txt` artifact written, rejects non-CGMES cases, and excludes the other modes), `ReactiveCapabilityCurve` interpolation, remote-regulating machines on a synthetic fixture (default held-PV fallback, `machine_control = true` controller wiring + control-loop convergence, and the voltage-held-target fallback), test-set alias/border bookkeeping, placeholder guards + `AsynchronousMachine` Stage-0 mapping (FullGrid warns about its X.99 shunt/tap-row fillers and solves from a flat start; MiniGrid's three motors close the former SV gap, asserted SV-tight), and — gated on the local ENTSO-E cache with explicit skip messages — MicroGrid/Assembled import + `runpf!` + `compareWithSV` validation, Stage-2 tap controllers, PSEI PST fixtures (tabular table wiring and the documented end-2 angle-sign deviation), the RealGrid tabular-PST SV-start-mismatch regression, the CGMES fixed-reference self-check (SV voltages reach the solver verbatim even against a `flatstart: true` base config, `self_check.log`/`self_check_residuals.csv` artifacts, `cgmes.log` no-SV bus count), `cgmes_import.start_values` selection (flat/sv decision line in `run.log`/`cgmes.log`, precedence over a hostile `power_flow.flatstart`, mandatory `sv_compare.csv`/`sv_compare_flows.csv` artifacts + summary metadata in both modes, MATPOWER negative case with no effect/artifact/metadata), the `runShortCircuit!` MicroGrid plausibility sweep (finite positive Ik'' on every bus, max ≥ min per bus, no motor/defaulted-data flags — consistent with `shortCircuitCoverage`), the `short_circuit_mode` service run (CSV artifacts with the short-circuit schema, coverage report in `run.log`, MATPOWER rejected with `short_circuit_requires_cgmes`, mutual exclusion with `diagnose_mode`, real-delivery data-check), and `normalPF` → `participationFactor` arrival (MicroGrid zero-factor mapping, MiniGrid positive factor); the FullGrid placeholder-guard block additionally asserts `cgmes_import.placeholder_guards = strict` aborts the import with the offending object named |
| `cgmes_export` | `test/test_cgmes_export.jl` | CGMES export identity and roundtrip contract (`writeCGMESFiles`, EQ+TP+SSH+SV + optional delivery zip): structural keys in `net.cgmes_ids` with lexicographic bus-pair normalization and first-seen parallel-line numbering, deterministic uuid5 minting (id equals uuid5 over the key), byte-identical re-export and independent-build determinism with a pinned `created` stamp, renaming a component never changes its mRID, the duplicate-mRID guard aborting before any file is written (both offending keys named), the power-flow-identical self-roundtrip over every exported class (transformer with off-nominal ratio + 2° phase shift + ratio-tap machinery, PV machine, slack injection, SVC with rating-derived Q limits, load, shunt, bus link → export → `importCGMES` → `runpf` matches the original within 1e-6 and starts from the exported SV in ≤ 2 iterations, empty notice list), tool provenance in every file (`Generated by Sparlectra.jl v<version>` comment + `md:Model.description`, stamp = pinned `created`), SSH content (machine sign convention, `referencePriority` on the slack unit, `RegulatingControl.targetValue` in kV, shunt sections, `bPerSection`) and SV validation (`compareWithSV` of the re-import: dvm/dva < 1e-12, flow rows < 1e-9), zip re-import, and — gated on the local ENTSO-E cache with an explicit skip message — the MicroGrid roundtrip proving imported TN/ACL/PT/PT3/EC/SM/SH mRIDs survive an export (canonical form without the RDF underscore prefix) even after renaming every line, the 3W star reassembled into one three-end transformer (star TN absent from the export), the re-import rebuilding the identical electrical model (every branch parameter equal incl. the PST angle), the short-circuit evaluation of the re-imported delivery reproducing the original Ik'' rows and flags for both cases, `cgmesLineShortCircuitData` feeding harvested zero-sequence line attributes into the EQ profile, and the `export_cgmes` service-run path (exactly one delivery-zip artifact, no loose profile files, the zip re-imports with the original mRIDs, `cgmes_export_*` metadata incl. notices, `run.log` line) |
| `apslf` | `test/test_apslf.jl` | APSLF (AnalyticLoadFlow.jl) integration: `power_flow.solver`/`apslf`/`apslf_start` config validation (including the `solver=apslf` + `apslf_start.enabled` conflict), outer-loop-controller rejection with `solver=apslf` (no AnalyticLoadFlow.jl needed for this check), Web UI form parsing of the new fields into the effective config, and the "AnalyticLoadFlow.jl not installed" error path. The adapter-mapping (`PFModel` -> AnalyticLoadFlow spec, PF ordering), standalone case14 convergence-vs-NR, and start-value-generator guard/`nr_polish=false` tests additionally require AnalyticLoadFlow.jl to be resolvable in the active environment (a weak dependency, not part of the default project); they skip with an informational message otherwise instead of failing. |

### `dtf_extended` group

`dtf_extended` is a normal `extended`-only group. It owns exactly five DTF test files and invokes each of their run functions once:

| File | Run function | Main checks |
|---|---|---|
| `test/extended/test_dtf_importer.jl` | `run_dtf_importer_tests` | DTF format parser and direct Net-builder coverage, including voltage-level-index branch conversion, transformer controls, bus-type semantics, parsed outage metadata, and the persisted typed phase-tap model on controlled windings (FOR001E skew/longitudinal regulators carry `phase_taps` with the winding connection angle) |
| `test/extended/test_dtf_for002_validation_example.jl` | `run_dtf_for002_validation_example_tests` | DTF format -> current Julia compatibility module `DTFImporter` -> `Net` -> power-flow validation example smoke coverage against FOR002 diagnostics; verifies generated CSV/Markdown artifacts, lightweight default result/concise CLI output, explicit detailed diagnostics mode, and does not invoke the fast suite |
| `test/extended/test_dtf_for002_outage_validation_example.jl` | `run_dtf_for002_outage_validation_example_tests` | DTF-listed outage validation against FOR002 outage reference reports for the Testnetz13 cases |
| `test/extended/test_dtf_matpower_export_validation_example.jl` | `run_dtf_matpower_export_validation_example_tests` | DTF format -> `Sparlectra.Net` -> existing `writeMatpowerCasefile` -> MATPOWER import roundtrip validation for the Testnetz13 base case and outages listed by the DTF file |
| `test/extended/test_dtf_api_webui_integration.jl` | `run_dtf_api_webui_integration_tests` | DTF format input through the API and Web UI paths |

`extended` includes `dtf_extended`, and `all` runs `fast` once followed by the complete `extended` profile once, so the DTF tests run exactly once per `extended` or `all` invocation. Optional external FOR001/FOR002 fixtures continue to skip cleanly (as `Broken`, not failed) when absent.

`test/extended/runtests_extended.jl` remains an optional standalone DTF-focused runner. It includes the same five files directly and is useful for iterating on DTF coverage without running the rest of the extended profile:

```bash
julia --project=. test/extended/runtests_extended.jl
```

## Native DTF validation examples with FOR002 reference reports

The native DTF examples with FOR002 reference reports validate Testnetz13 through the direct DTF format path:
the current Julia compatibility module `DTFImporter` reads the DTF format with `DTFImporter.read_dtf` -> `DTFImporter.build_net` -> `runpf!`. They deliberately avoid MATPOWER import/export and the generated FOR001 builder so that native DTF format parsing, Net construction, and solved branch-flow reporting are exercised directly. FOR002 is used as a legacy textual reference report.

External FOR001/FOR002 validation datasets are not shipped with Sparlectra.
Place local validation files under `data/DTF/` or pass explicit paths with
`--dtf-file` and `--for002-file`. If the files are absent, the optional
validation scripts stop with a clear missing-data message instead of
substituting unrelated data or claiming that validation was executed.

`examples/run_val_dtf_suite.jl` is the shared CLI runner for all three checks
(plus the DTF import audit). It bundles cases and modes into one command and
dispatches to library modules under `examples/internal/` (`dtf_validation_base.jl`,
`dtf_validation_outages.jl`, `dtf_validation_matpower.jl`, `dtf_validation_audit.jl`),
so there is one implementation of each check shared by the CLI, the suite, and
the extended test files.

Run the base-case validation with:

```bash
julia --project=. examples/run_val_dtf_suite.jl --mode=base --case=A --data-dir=data/DTF --output-dir=examples/_out/dtf_for002_native_validation --write-csv=true --write-markdown=true
```

Run the outage validation with:

```bash
julia --project=. examples/run_val_dtf_suite.jl --mode=outages --case=A --data-dir=data/DTF --output-dir=examples/_out/dtf_for002_native_outages --write-csv=true --write-markdown=true
```

Each command writes concise console output plus CSV and Markdown files in the requested `--output-dir`. The Markdown files summarize the run, while CSV files keep row-level bus, generator, branch, KCL, state-residual, and metric diagnostics.

The MATPOWER export validation example checks a different question: whether an
already-built native DTF `Net` can be exported by Sparlectra's existing
MATPOWER exporter and re-imported without materially changing the solved
Sparlectra result. Its required path is
`DTFImporter.read_dtf` -> `DTFImporter.build_net` -> `writeMatpowerCasefile` ->
`createNetFromMatPowerFile`; it does not implement a DTF-specific exporter and
does not use FOR002 as the primary roundtrip reference. Run it with:

```bash
julia --project=. examples/run_val_dtf_suite.jl --mode=matpower --case=A --data-dir=data/DTF --output-dir=examples/_out/dtf_matpower_export_testnetz13 --write-csv=true --write-markdown=true --write-matpower=true --run-outages=true
```

Each `examples/internal/dtf_validation_*.jl` module is also directly runnable
as its own single-purpose CLI entry point (for example
`julia --project=. examples/internal/dtf_validation_base.jl --dtf-file=... --for002-file=...`),
which is what the extended example tests drive for exact console-output
coverage; prefer the suite above for everyday case/mode selection.

The command writes `dtf_matpower_export_summary.md`,
`dtf_matpower_export_metrics.csv`, bus/branch/generator comparison CSV files,
and exported MATPOWER `.m` cases in the selected output directory. The exporter
uses Sparlectra's established optional MATPOWER metadata fields
`mpc.bus_name`, `mpc.branch_name`, `mpc.branch_kind`, and
`mpc.for001_contingencies` when that information is available, so DTF bus names,
stable branch labels, DTF branch kind (`L`/`T`), and the parsed outage cards can
be recovered by `MatpowerIO.read_case_m`. Standard MATPOWER bus, generator, and
branch data remain sufficient for solving files that do not contain that
metadata. The roundtrip proves preservation through Sparlectra's MATPOWER
export/import path; it does not independently certify agreement with external
FOR002 outage blocks.

The roundtrip import disables MATPOWER PQ-generator controller reinterpretation
for this diagnostic so DTF PQ generators remain fixed injections. The exporter
also writes MATPOWER TAP as `0.0` for line rows and as the explicit Sparlectra
ratio for transformer rows, matching MATPOWER's line/transformer convention.

Important metrics:

- `converged`: whether `runpf!` reported successful convergence.
- `iterations`: Newton iterations used by the native solve.
- `final mismatch`: infinity-norm mismatch recomputed from the solved state.
- `max voltage deviation`: largest voltage magnitude or angle difference from FOR002 printed bus values.
- `max branch P/Q deviation`: largest directed branch endpoint-flow difference from FOR002.
- `max generator/slack Q deviation`: difference between solved generator/slack reactive output and FOR002 generator-Q reporting.
- `state residual P/Q`: injection residual from forcing FOR002 printed voltage magnitudes/angles into the native Ybus and comparing the calculated injections with the FOR002 bus table. Reported only in the CSV diagnostics (`dtf_state_residual.csv`, `dtf_validation_metrics.csv`, and the outage equivalents), not in console or Markdown summaries.

Current Testnetz13 interpretation:

- Branch-flow deviations are small and are the strongest validation signal for the native DTF path.
- Slack Q is solved by the power flow and should not be compared with the specified input Q as if it were fixed.
- State residuals are gross-error diagnostics only: their rounding-noise floor (FOR002 prints rounded voltages/angles, and transformer-adjacent buses amplify tiny differences) sits far above real model deviations, so they are kept in the row-level CSVs but excluded from console and Markdown summaries.
- Outage validation currently executes only the outages listed in DTF.
- FOR002 may contain more outage blocks than FOR001 lists; unmatched FOR002 blocks are treated as reference text, not executed scenarios.

## Offline and runtime expectations

The default `fast` profile is intended to be offline-safe and should not download MATPOWER cases or run benchmark loops.
The Q-limit large-case comparison workflow is a manual diagnostic tool: it resolves optional large MATPOWER cases through Sparlectra's case registry/download support, compares `case × start_profile × qlimit_mode`, and uses its CSV/JSON summaries as the primary output. Real runs for cases such as `case13659pegase.m` and `case_SyntheticUSA.m` are intentionally excluded from the fast tests; automated coverage uses stubs for case resolution and API execution.
The experimental large-case Q-limit comparison test block is suppressed from the normal fast profile; run `test/test_qlimit_large_case_comparison.jl` manually when maintaining that tool.

Default fast-profile output is intentionally compact: the runner prints the selected profile, one `[n/8]` marker per group, and Julia's final test summary. MATPOWER import diagnostics, auto-profile tables, runtime casefile banners, Q-limit tables, and similar artifact-oriented diagnostic blocks are suppressed in normal test stdout so progress remains scannable.

Fast profile example on Windows / Julia 1.12.6: 934 tests passed in approximately 95 seconds. Runtime is machine-dependent and is not a CI threshold.

Use an explicit verbose opt-in when debugging a noisy test path:

```bash
SPARLECTRA_TEST_VERBOSE=1 julia --project=. test/runtests.jl
julia --project=. test/runtests.jl --verbose
```

Verbose mode does not change the selected profile. Existing profile selection remains available, for example:

```bash
julia --project=. test/runtests.jl fast --verbose
julia --project=. test/runtests.jl extended --verbose
```

The `fast` profile runs normal unit and focused integration tests. The `extended` profile runs only extended repository-hygiene, documentation-coverage, example, and fixture-heavy checks. The `all` profile runs fast followed by extended. Set `SPARLECTRA_TEST_GC_BETWEEN_GROUPS=1` to request a GC cycle after each completed group; per-group output reports elapsed seconds, allocated MiB, and GC seconds.

The `extended` profile may include MATPOWER/example/output-heavy tests and native DTF diagnostic-example checks with FOR002 reference reports. These tests stay isolated from the default profile.

Use `fast` during normal development. Use `extended` before merging changes that affect configuration, MATPOWER import, output formatting, performance reporting, or broader integration paths.

## Fast-profile volume review

The fast profile currently contains a mix of true unit/smoke coverage and several integration-style service/UI paths:

- True unit or focused smoke tests: configuration key/value validation, MATPOWER auto-profile decision rules on tiny synthetic cases, rectangular solver API checks with small fixtures, core model invariants, state-estimation smoke/regression cases, and control-loop unit/regression checks.
- Integration-style tests that remain in fast because they protect recent public behavior: focused API service request/metadata smoke checks, Web UI form rendering and stubbed routing, allowlisted documentation/help routing, operation-log safety, run deletion safety, and one small API service run.
- Heavier or broader tests already isolated in extended: MATPOWER example runner coverage, synthetic-grid regressions, configuration documentation consistency, PV residual integration coverage, structural remove/delete behavior, real asynchronous Web UI/API job lifecycles, and artifact lifecycle matrices.
- Expensive or duplicate candidates to watch: repeated API/Web UI artifact-generation parity checks, broad service-path status/recovery assertions that overlap between `test_api.jl` and `test_webui.jl`, and any future large-case or repeated auto-profile scans. These should move to `extended` if they become slow, require network/cached large cases, or duplicate a smaller fast regression.

No tests were moved in this review. The current fast-profile volume is acceptable as long as default stdout remains quiet and the service/UI cases continue to use small offline fixtures. Future large-case regressions such as `case13659pegase.m` should use a small reproducible proxy in fast and keep the real large-case check in extended/manual verification unless the case is already cached and cheap.

## Pre-merge verification gate (config / MATPOWER / output / performance / docs changes)

For branches that touch central configuration, MATPOWER runner behavior, output routing/formatting, performance reporting, or documentation/config consistency, complete this checklist before merge:

### Bash

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl
SPARLECTRA_TEST_PROFILE=extended julia --project=. test/runtests.jl
julia --project=docs docs/make.jl
```

### PowerShell

```powershell
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl

$env:SPARLECTRA_TEST_PROFILE="extended"
julia --project=. test/runtests.jl
Remove-Item Env:SPARLECTRA_TEST_PROFILE

julia --project=docs docs/make.jl
```

This keeps the default local workflow fast while making the extended profile and docs build an explicit pre-merge gate for integration-heavy changes.

## Output-summary regression note

The fast profile includes a regression for `printACPFlowResults(...; toFile=true, result_mode=:summary)`.
It verifies that the result file is closed/flushed before the function returns and that the summary contains Q-limit counter labels.
Equivalent environment-variable usage remains supported:
```bash
SPARLECTRA_TEST_PROFILE=extended julia --project=. test/runtests.jl
```

## Progress output

The runner prints a lightweight progress view, for example:

```text
Test framework: fast
[1/8] core_model
[2/8] powerflow_rectangular
...
[8/8] controls
```

Julia's final `Test Summary` remains unchanged and visible at the end.

## Rectangular/Q-limit diagnostics in tests

The rectangular convergence and Q-limit active-set diagnostic block is not printed in normal test runs.
Those diagnostics remain available only through explicit diagnostic requests (for example solver `verbose > 0` paths used during focused debugging).

## Experimental/internal DTF Web/API input path

The PowerFlow service and Web UI include an experimental/internal DTF format input path for diagnostics and validation. This path is deliberately cautious and is not yet announced as a public supported file format. Use `case_format = :dtf_for001` for explicit native input, or `case_format = :auto` only when the FOR001 markers are unambiguous; ambiguous `.DAT` files are rejected instead of being silently interpreted.

The native API path uses `DTFImporter.read_dtf` -> `DTFImporter.build_net` -> `Sparlectra.Net` -> `run_sparlectra`/`runpf!` and does not solve through a MATPOWER intermediate conversion. The Web UI places the selector in the advanced/internal **Input format** section with the cautious label "DTF diagnostics (experimental/internal)". FOR002 is treated as a reference/result file rather than a runnable input case: FOR002-like `.DAT` files are hidden from the primary case selector, can be selected or typed as an optional FOR002 reference file when available in the case cache, and must not be auto-paired with FOR001. DTF-listed outages can be requested explicitly as all outages or selected outage labels/indices; the default remains base-case-only.

Generated artifacts use the existing PowerFlow run artifact mechanism. Stable DTF-path artifact filenames include `dtf_import_summary.md`, `dtf_import_summary.csv`, `dtf_for002_base_comparison.md`, `dtf_for002_base_metrics.csv`, `dtf_native_matpower_export.m`, and per-outage files such as `dtf_outage_1_summary.md` and `dtf_outage_1_metrics.csv`.

Limitation: DC lines, HVDC links, and active MATPOWER `mpc.dcline` data are not modeled by this native DTF/MATPOWER power-flow path. They fail clearly with structured unsupported-DC-line diagnostics instead of being approximated or dropped.

# Test Suite

Sparlectra uses profile-aware test loading through `test/runtests.jl`.
Profile selection precedence is:
1. CLI argument (`julia --project=. test/runtests.jl <profile>`)
2. `SPARLECTRA_TEST_PROFILE`
3. default `fast`

## Test profiles

| Profile | Command | Scope | Intended use |
|---|---|---|---|
| `fast` (default) | `julia --project=. test/runtests.jl fast` | Everything that checks numerics, the network model, and the case formats (risk-based split) | Normal development and default CI |
| `extended` | `julia --project=. test/runtests.jl extended` | Surfaces, artifacts, documentation, fixtures: Web UI, service lifecycle, examples, hygiene, CGMES/DTF fixture suites | Before merge or after broad integration changes |
| `all` | `julia --project=. test/runtests.jl all` | Fast followed by extended | Complete local or CI verification |

## Group timings measure compile distribution, not test cost

Per-group seconds in a profile run are dominated by Julia compilation, and
which group pays it depends on run order: the first group that touches a
code path compiles it for everyone after. Measured on the maintainer
machine (same process, second run warm):

| group | cold | warm |
|---|---|---|
| `scf` | 138.0 s | 0.6 s |
| `programmatic_api_extended` | 138.4 s | 0.3 s |
| `state_estimation` | 38.1 s | 0.2 s |
| `auto_powerflow` | 71.6 s | 3.2 s |

Judge a "group got slower" observation against warm times or the sum over
both profiles, never against a single group's cold seconds; deleting tests
barely moves these numbers. The lever for wall-clock time is compilation
reuse (a sysimage-based runner, or a cached Julia compile directory in CI).

Note: CI invokes the fast profile, so since the risk-based resort the
numerics groups (rectangular solver, factorized linear solver, PV voltage
residuals, 3WT phase taps) are covered in CI. The CI workflow file itself
follows in a separate workflow-only commit (workflow files never mix with
content commits).

## Fast versus extended ownership

| Area | Fast coverage | Extended coverage |
|---|---|---|
| Core/model | Small constructors, transformer checks, bus/prosumer semantics, link behavior, representative rectangular PF/Q-limit regressions, and small sparse fallbacks. Large sparse Ybus and large MATPOWER matrix-body checks are extended-only. | Large sparse Ybus smoke, large MATPOWER matrix-body scanner, synthetic/stress grids, and longer integration examples. |
| API | Serialization and transport helpers, path and validation safety, one successful small API run, one pre-solver failure, one numerical/island failure, Solver-time and Total-time contracts, critical DC-line default and strict-rejection smokes, and one small independent-island regression. | Exhaustive CSV/export matrices, repeated artifact inventories, island artifact-content matrices, persistent history/delete/reload lifecycle coverage, and repeated presentation/performance-log modes. |
| Web UI | Form parsing and backend validation, result rendering, active and terminal timing cards, commit-span omission, tolerance-step hook, path traversal rejection, DTF upload role classification, primary-case and FOR002 selector filtering, the buildSysimage one-call dry run, the sysimage launcher decision (all four staleness reasons including a src file newer than the image, and the no-terminal default that builds), the agreement between that launcher verdict and `Sparlectra.webui_sysimage_problem` across every fixture state, the sysimage build-progress reader (running, gone-stale, done, failed) and the Sysimage page in each of those states including its rebuild guard, a parse check of the generated app CLI (run/se/n1 commands, config flags) driven from the checkout tool, and stubbed route checks without a real asynchronous solver run. | Case-profile persistence with asynchronous jobs, real run/result polling, artifact preview/download/ZIP/history/delete lifecycle checks, browser-launcher matrices, socket/server lifecycle, and Markdown/help/documentation cross-products. |
| Documentation/hygiene | No repository-wide documentation/help scan in fast; only focused source-level smoke checks tied to edited paths. | Configuration documentation consistency and normalized tracked-path/content repository hygiene scans, plus the reference-page coverage assertion (every src/ Julia file on exactly one reference page or on the exclusion list), a scan that no docstring is separated from its definition by a blank line (Julia then drops it silently, and only the documentation build notices, by failing on an unresolved `@ref`) and, in the maintainer's working repository, a boundary check that the tracked tree carries no private working files (it reports plainly that it did not run where that tooling is absent). |

`Pkg.test()` uses the same test runner and therefore the default `fast` profile unless `SPARLECTRA_TEST_PROFILE` is set:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Large case files

Large MATPOWER and DTF cases are not part of the repository. The test suite
looks for them in the shared case directory resolved by
`Sparlectra.large_cases_dir()`:

1. `SPARLECTRA_LARGE_CASES_DIR`, if set (override for CI and special setups)
2. otherwise the Web UI user case directory
   (`~/.local/state/sparlectra/webui/data/mpower` on Linux,
   `%LOCALAPPDATA%\Sparlectra\WebUI\data\mpower` on Windows)

Test groups that need a large case check availability and report `SKIPPED`
when the file is absent. The fast profile does not depend on any large case.
Tests only read from this directory; generated files go to `mktempdir`.
`data/mpower/` in the repository holds fixtures only, and the gate check
fails on untracked files there.

**The assertion count therefore depends on that directory, and a measurement
has to say which count it was taken with.** With the large cases present the
fast profile runs 7179 assertions; without them, which is CI and any checkout
whose user directory holds no cases, it runs 7121. Both are correct runs of the
same suite, and the `SKIPPED` lines say which one happened. A compile-time
baseline table does not, so it has to record the number alongside its timings,
or the next measurement compares two different suites.

`SPARLECTRA_TEST_SKIP_GROUPS` (comma-separated group names, default empty)
excludes profile groups from a run. It exists for the app build workload
(`tools/app_workload_runtime.jl`), which traces the fast profile as a
PackageCompiler precompile run; every skipped group is printed at the start
of the run so a skip can never be mistaken for a pass. Do not use it to get
a green development or CI run.

`tools/run_gates.sh` starts Julia with `--startup-file=no`. A gate is not an
interactive session, and a personal `startup.jl` that loads Revise invalidates
472 method instances of the stock Julia system image, which the run then infers
again for nothing. Use the same flag when running `test/runtests.jl` by hand on
a machine whose `startup.jl` loads packages.

The **sysimage** workload does not trace the test suite at all any more: it
runs a curated list of interactive paths once each (see
[Sysimage](sysimage.md)). Tracing the whole fast profile dominated the
build time without covering a Web UI path the curated list misses.
`SPARLECTRA_SYSIMAGE_TRACE_TESTS=1` restores the old behavior for a
maintainer comparison.

## Fast profile groups

| Group | Files | Main checks |
|---|---|---|
| `core_model` | `test/testgrid.jl` | Core net construction and validation, small inline MATPOWER import helpers, file-based MATPOWER projected-start normalization, PV/PQ lock-ID resolution, Q-limit enforcement-mode parsing and classical base-failure/no-reenable dispatch checks, rectangular nonfinite mismatch and status-row diagnostic preservation, Jacobian condition-number estimator checks (exact diagonal case, Hager bound, near-singular flagging, input rejection), link handling, shunts, reporting/output checks, and summary-file output regression. Large sparse Ybus checks and large MATPOWER matrix-body scanner coverage are extended-only. |
| `terminal_status` | `test/test_terminal_status.jl` | Per-terminal branch status (r0.9.10, one-sided open branches): the dangling-node reference anchor (a `:open_to`/`:open_from` reduction matches a solve with an explicit zero-injection auxiliary bus to 1e-10, for a line and for a transformer with off-nominal ratio and phase shift), bitwise equivalence of both-flags-open with `status = 0` in the Y-bus, solver toggles between runs (sparsity-pattern invalidation), isolation of the open-side bus in the island report, the MATPOWER export substitution (BR_STATUS 0 plus exact `Y_in` bus shunt, `open_terminal=` marker, voltage-identical reimport), the result surface (`open@to` marker, `Open terminals` header count, `terminal_state`/open-end voltage report columns), and the state-estimation guard rejecting flow measurements on partially open branches |
| `powerflow_rectangular` | `test/test_solver_interface.jl` | Rectangular power-flow API behavior, sparse-only solver path, AC-island detection/reference validation/independent solving including aggregate all-island convergence status and iteration accounting, Q-limit and typed configuration entry checks, current-iteration start rejection diagnostics on a tiny synthetic case, including accepted and rejected framework control-status composition, ordered local MATPOWER batch execution, the thin `run_acpflow` alias, legacy-keyword rejection, the optional merit-function Armijo line search (value/scaling, acceptance/fallback/active-set-skip reasons, config validation, merit.enabled=false regression, and a converging flatstart/Q-limit-switching functional case), the optional scaled-Newton trust-region step control (config validation, acceptance/collapse unit coverage, trust_region.enabled=false regression, a converging functional case with log/diagnostics assertions, and radius-collapse non-convergence reporting), and wrong-branch detection output visibility across ACPFlowReport metadata, the AC island diagnostics CSV, and the console summary line, including the highest-voltage-level scope regression (a dip below the top level never triggers; the same dip on the top level does), and the AC rescue ladder / DC fallback (`power_flow.rescue` recovers a poisoned start and records the winning strategy, `power_flow.dc.fallback` leaves a DC state with the AC status honestly non-converged, defaults leave both mechanisms off) |
| `factorized_linear_solver` | `test/test_factorized_linear_solver.jl` | umfpack_reuse linear-solver selection and factorization-reuse behavior for the rectangular power-flow Jacobian path, including rejection of the removed klu value |
| `pv_voltage_residuals` | `test/test_pv_voltage_residuals.jl` | PV-voltage residual behavior, angle-preserving voltage-setpoint starts, phase-shifted PV integration coverage, and related solver diagnostics |
| `3wt_phase_taps` | `test/test_3wt_phase_taps.jl` | `create3WTWindings!`'s `phase_tap_side`/`phase_taps` keywords — existing-behaviour snapshot, single-winding attachment, ratio+phase coexisting on one winding, and validation errors |
| `configuration` | `test/test_configuration_coverage.jl`, `test/test_runner_helpers.jl` | Configuration-key coverage, safe refresh checks including current-iteration start defaults, forwarding checks, start-voltage value-domain validation (including `cgmes_import.start_values` and `cgmes_import.placeholder_guards` parsing: defaults `flat`/`warn_skip`, `sv`/`strict` accepted, invalid values rejected with the key name), and test-runner output-mode helper checks |
| `matpower_metadata` | `test/test_matpower_metadata.jl` | MATPOWER parser metadata fields, legacy bus sorting of `bus_name`, opt-in bus-name import, branch-kind override, branch metadata retention, FOR001 contingency mapping, default/opt-in `mpc.dcline` behavior including PF/loss/Q/V/Q-limit mapping and voltage-control safeguards, `writeMatpowerCasefile` `write_solution` column count/VM-VA/marker behavior, the `mpc.sparlectra.tap_changer_model` roundtrip marker preventing double tap-impedance correction, the `mpc.sparlectra.links` busbar-coupler extension (parse, link contraction in the power flow, export roundtrip, malformed-block rejection), the `mpc.sparlectra.tap_changers` nameplate extension (neutral-vs-live mapping onto the branch tap fields, ratio-less phase shifters, the Delta-u additional-voltage stepper with cascade-derived live position and roundtrip, exclusive phase grids, generated-flow consistency, export roundtrip, rejection matrix), and `net.name` staying the case name (not a bus's original name) when `mpc.bus_name` metadata is present |
| `programmatic_api` | `test/test_api.jl` | Focused GUI-ready power-flow API coverage for serialization and transport helpers, path and validation safety (including run-index acceptance of symlinked output roots, skipped on Windows), one successful small API run, one pre-solver failure, one numerical/island failure, Solver-time and Total-time contracts, critical DC-line default and strict-rejection smokes, and one small independent-island regression. Exhaustive CSV/export matrices, full artifact inventories, persistent restart/history/delete lifecycle matrices, repeated performance-log modes, complete island artifact-content matrices, and the fixed-reference self-check contract (verbatim-start summary in `self_check.log`, per-bus `self_check_residuals.csv`, the sp_case14 residual property plus the cache-gated historical case14 anchor) are extended-only in `programmatic_api_extended`. |
| `state_estimation` | `test/test_state_estimation.jl` | WLS state-estimation behavior and observability-oriented regressions, incl. `measurement_jacobian` (labeled H: row/column counts against the global check, flow rows named by oriented branch, bus rows by bus, every row touches a state, rank agrees with `evaluate_global_observability`), the `ImagMeas` current-magnitude rules (prediction roundtrip, 3 sigma value gate, shunt-bay variant, observability invariance), bad-data localization stage 1 (`wii` exposure, sequential elimination with trace/stop reasons, ZIB/SHDERIV protection, optional K-matrix report), shunt estimation (case A recovery with freeze guards, write-back gate, voltage sensitivity; case B back-calculation with sign and sigma propagation), the phase-3 link rules (SE on the contracted net with member sync, LINKAGG aggregation and partial-cluster exclusion, link measurements never in the WLS, bitwise KCL regression of the W2 allocation plus measured ring split, `se_view` classification), the phase-4 hardening (robust R stages with solve/diagnosis separation incl. the generalized knees k1/k2 with bitwise 3/6 equivalence and continuity at t = k1, the replacement suppression judged on the converged state with statistics on original sigmas, the k_eliminate exposure, Wilson-Hilferty band test with `:high`/`:low` reasons, diagnostics state unification, two-stage FD-aware observability), the phase-5 chain (measurement CSV v1 roundtrip with atomic import incl. `branch_nr` parallel-branch disambiguation and the legacy header, `runpf_from_se!` start modes with slack pickup, untouched model loads, and the measurement/model-discrepancy case separating `:se_state` from `:se_snapshot`), island-wise SE (two-island partition with per-island references and band verdicts, unmeasured-island skip, closed-link island fusion, merged chain registration), the transformer tap estimation (PF-equivalence gate at 1e-12 for the unstamped cascade overlay in all three release modes plus the live-tap KCL consistency regression, generated measurements stay self-consistent when tap_ratio leaves neutral; ratio/PST/two-regulator roundtrips recovering injected mechanical steps with the mandatory fixation run incl. the off-grid half-step case and the `:offgrid_tap_residual` band-flip note; the release guards: machine-transformer detection, structural bridge freeze without a far-side voltage bitwise equal to no release, partial `:both` column re-packing with bitwise stamp restore; and the two-island aggregation with per-island rows, summed fixation J/dof, and bitwise cross-island isolation; the completion coverage: `updateTaps` write-back with bitwise protection, the PMU current-phasor synergy shrinking the phase-regulator spread, bad data next to a released trafo (elimination names the row, tap stays within half a step, diagnostics on the tap-unified net), the machine-trafo back-calculation `calcMachineTrafoTapFromSE` incl. mutual exclusion, and the `se_view` release listing), the Takahashi diagnostics (dense/sparse path equivalence on three case classes, FD exact-zero sparsity, guard fallback with warning, K pinning the dense path), and the DTF case through the SE service (the requested format travels with the request; without it the run fails naming `dtf_for001` instead of a bare "unknown format") |
| `topology_validation` | `test/test_topology_validation.jl` | Topology validation (all advisory): the stage-1 precheck matrix (each finding kind fires exactly its synthetic case, a clean set fires none) incl. the false-positive guards (legitimately unloaded branch with quiet neighbours, open link disagreeing freely, partially measured node skipped, shunt buses balanced with the measured voltage), the advisory contract (findings present AND the estimation still runs), the stage-2 fingerprint (model-closed/actually-open transformer: elimination exhausts, band `:high`, the finding names the right station; a curable 10-sigma error never classifies as topology), the stage-3 hypothesis test (`:auto` candidates, the true toggle ranked first and supported, honest ambiguity on the small net, the input bitwise untouched, the island case), and the singular-normal-equations regression (dense LS fallback with a raised `svd_max_n`, parallel released taps converging on a least-squares split instead of crashing) |
| `auto_powerflow` | Automatic power-flow mode: read-only feature extraction on four case classes, ordered-rule strategy selection with step-control exclusion, the escalation ladder (cap, skip reasons, L4 locks, L5 evidence gate, seed-only L6), DC-fallback labeling, precedence, decision-log artifact and metadata mirror, hint-generator ids, the startup latency hint, and the auto-profile precedence level (D11: explicit non-default user keys survive every recommendation, applied and skipped pairs reported, resolve_config carries the level and explicit overrides beat it). |
| `scf` | Sparlectra Case Format writer: deterministic ids with byte-identical repeat and rebuild exports, the PGM component subset, SI unit conversion against hand-computed references (line ohm/farad on the to-side base, live k/theta, appliance W/var), the namespaced block (roles, tap cascade with the Delta-u controller, names, config, opt-in start state), measurement-to-sensor pairing with complete row accounting, the FACTS base-impedance guard, the round-trip identity (byte-identical re-export, identical Y-bus, measurements, power-flow iterations and voltages, identical state-estimation objective and state), the permanent v1 fixture, reader validation (unknown keys, unsupported components, dangling references, duplicate ids, batch datasets, malformed JSON), the controller round trip through the declarative schema, the study definitions with their load-time validation, the study runs themselves (an N-1 case list and exclusions taken from the file, a rejected multi-outage case, a short-circuit run on the file's own sources with its case, c-factor and explicit bus sweep, and the missing-source path), the round trip on real MATPOWER cases (case14 and case57: byte identity, component names, bus order, tap ratios, shunt signs, Y-bus, unlimited ratings; cache-gated with a spoken skip since load_fixture_net, because their quirks have no shipped equivalent), the own-name-vs-external_id channel on the shipped sp_case5 (reference name applied, external_id restored, re-export preserves both), power-grid-model interoperability (a file with no namespaced block whose sources become the reference, its solved voltages checked against an independently derived model), the strict PGM writer mode, the three-winding nameplate form with its named types and its rejections, the shunt-state and PGM fault blocks, the case-file upload including a rejected foreign JSON, the download route with its refusals (directory escape, missing file), and the service case resolver, SCF as a first-class case format (detection, selector, framework and service runs, case-file configuration precedence), the logging surface in the configuration allowlist, and the Web UI export route including its failure paths. |
| `dc_powerflow` | `test/test_dc_powerflow.jl` | Standalone DC power flow (`rundcpf!`, `power_flow.solver=:dc`): reference-value comparisons against an independently-computed script (3-bus hand-verified fixture always; the MATPOWER case9/case14 oracle tables cache-gated with spoken skips since load_fixture_net) for bus angles and branch flows, the hand-checkable sp_case5 DC dispatch (slack exactly 32 minus 12 MW), phase-shifter (`Pfinj`) coverage via a synthetic phase-shifting-transformer fixture, lossless power-balance property checks on the tracked warmup case and the shipped sp_case5, two-island independent-solve coverage, `power_flow.dc.*`/`power_flow.solver` config validation, active-outer-loop-controller rejection (mirrors the `apslf` rejection test), `run_sparlectra` dispatch (`method=:dc`, `diagnostics.solver=:dc`), `angle_reference_deg` exact-uniform-shift behavior, the `dc_pf_status`/`rectangular_pf_status` registry separation, `rundcpf!(net; seed_ac_start=true)` DC-seeded AC Newton-Raphson convergence, and the one-specification-source contract (#323: a node-sum edit moves neither solver, a prosumer edit moves both, DC injections equal the real part of `buildComplexSVec`) |
| `distributed_slack` | `test/test_distributed_slack.jl` | Distributed active-power slack: disabled-path bit-identity to the classical solver, per-bus residual conservation (`alpha .* lambda_P` pattern, REF/load buses untouched), REF-equivalence (`explicit` weight fully on the reference bus reproduces the classical solution and `lambda_P` equals the classically REF-absorbed power), mode agreement (`pg_weighted` vs matching `explicit`), MATPOWER `APF` column import (`:imported` shares proportional to APF; APF presence never changes disabled results), fallback behavior (`error` throws, `ref_only` warns and matches classical), Q-limit interaction (PV→PQ switch keeps participation), the config-load regression for the `weights: {}` placeholder plus block-style weight tables (the minimal YAML reader has no flow-mapping support), and independent per-island `lambda_P` via the AC-island solve path. The CGMES `normalPF` arrival assertions live in the extended `cgmes_importer` group (MicroGrid zero-factor mapping, MiniGrid positive-factor arrival, cache-gated with explicit skip message). |
| `island_diagnostics` | `test/test_island_diagnostics.jl` | AC island diagnostics reporting regressions (source: a 6209-bus CGMES delivery with 158 single-bus islands): the island failure message reports the failing island's own iteration count and stage (never the `iterations=0`/`stage=before_nr` fallback while a per-island record exists), the failing island is looked up in `:ac_island_solver_statuses` when the combined status is untagged, never-attempted islands render as `not_attempted` with `unavailable` mismatch fields instead of inheriting the failed island's statistics, per-island `solver.log`/`mismatch_history.csv` artifacts are written only for attempted islands, `q_limit_processing_status` no longer aliases `failure_reason`, and the single-island combined-status fallback keeps working; also the no-slack diagnostics (all-isolated vs. isolated-slack vs. never-registered messages) and the `power_flow.auto_slack` promotion (keyword and config path, ENI-over-generator ranking, no-op on a registered slack) |
| `short_circuit` | `test/test_short_circuit.jl` | IEC 60909-0 short circuit (`runShortCircuit!`): hand-derived analytic reference cases with the full derivation in the test comments (feeder reproduces its declared current at the connection point with the c-factor cancelling; series line impedance added on top; synchronous-machine x''d conversion incl. §6.6.3 fictitious resistance; kappa/i_p with the method-b cap; line charging proven irrelevant), IEC Table-1 c-factor selection and the scalar override, the safety-flag contract as a three-way assertion (documented default + log warning + per-row flag), asynchronous-machine skip with island-wide Ik''max lower-bound flag (max case only), `:no_source`/`:isolated` row statuses, `short_circuit.c_factor` config coverage (default, valid override, rejected out-of-band value), and the WebUI data-check gating helper (`_webui_case_has_short_circuit_data`: marker found / not found / cached / unresolvable path stays enabled); the Phase 3 parallel all-bus sweep identity guard (serial vs parallel rows `isequal`-identical for max/min case with `max_tasks=1` and auto on a three-island fixture and on the shipped sp_case60 with `buses = :all`; RAN/SKIPPED stated per run; plus the `sweep_method = :takahashi` selected inverse: machine-precision agreement with `:solves` on values, identical statuses/flags, serial and parallel, deterministic repeats, invalid method rejected) |
| `parallel_foundation` | `test/test_parallel_foundation.jl` | Thread-safety foundation and parallel island path of the multi-core work (Phases 1 and 2): `runtime.parallel.*` configuration defaults and validation (`enabled`, `max_tasks` auto/integer-string, `min_work_items`), the `parallel_max_tasks` resolver, the per-worker performance-profile helpers (`_perf_profile_child` excludes the orchestrator-only `:phase_callback`, `_perf_profile_merge!` sums timings under unchanged phase names, appends iteration rows, prefixes per-worker scalars, and swallows `nothing` children/parents), the solver status as a Net field (registry globals gone, AC/DC fields separate, `deepcopy` carries a working condest thunk), UMFPACK `copy(F)` finalizer safety, the startup summary `parallel:` line, the Web UI `runtime_parallel_enabled` checkbox wiring, and the `solve_parallel` island identity guard: bitwise-equal voltages and iteration counts vs `solve_independent` with `max_tasks=1` and auto on a five-island fixture, identical phase-name sets plus `parallel_wall_time`, per-island prefixed scalars, and the parallel failure semantics (all islands report their real status, no skips). Threaded assertions state RAN/fallback depending on `Threads.nthreads()`; the `--threads=4` battery exercises the parallel branch. |
| `scenarios` | `test/test_scenarios.jl` | Scenario patch model (scenario task step 1): D1 validation surface (status, set, scale with their target classes and fields), rejections with scenario name and op index (unknown id, class mismatch, inapplicable field, empty ops, duplicate names), the sparlectra.scenarios file round trip through the typed case, the legacy contingencies-to-scenarios mapping, the N-1 expansion equality with generateN1Branches and generateN1Generators on the tracked warmup case including exclusions, and step 2: apply!/restore! with the undo log (tap patch on an unregulated transformer changes the solver-visible ratio and restores it, a tap patch on a regulated transformer is rejected with the controller named, and the acceptance test: for the tracked warmup case, the shipped sp_case60/sp_case188 and the tracked SCF fixture, every N-1 branch and generator scenario leaves the working copy BITWISE equal to the base after apply! plus restore!, checked with scf_roundtrip_field_diffs), and step 3: the engine worker resets bitwise between evaluations WITH solves (branch outage, generator outage, patch scenario against the engine template), and `runScenarios!` with mode n1_all writes the identical result CSV as `runContingencies!` on the generated case list (the case118 byte gate against the 7438b6e fixture is the extended group `scenario_engine`); step 4 (screening surface): `:off` keeps the historical `ContingencyResult` rows and CSV bytes, `:flag` returns `ScenarioResult` rows with no violating case screened away, screened rows read `start_used = :screen` with an estimate attached, the screening CSV appends exactly the `screened`/`screening_estimate` columns, and `:only` reports estimates without full solves where an estimate exists |
| `demo_cases` | `test/test_demo_cases.jl` | The shipped demo cases (`data/scf/sp_case*`, task_demo_cases_v0100) as regression fixtures: each case is read through `import_case` like any other SCF file and all five run kinds are executed and compared against the tracked fixture in `test/fixtures/demo_cases/` (power flow: iterations, losses, per-bus vm/va to 1e-6 pu keyed by SCF node id; state estimation on the measurements the file itself carries: convergence, dof, objective; IEC 60909 short circuit at the fixture's two buses from the file's own feeder record; N-1 summary counts over `generateN1Branches`; the file's explicit scenarios through `runScenarios!` with per-scenario convergence and max loading), plus the self-built provenance line in the file, the per-case size budget, and the presence of the bad-data measurement variant |
| `contingency` | `test/test_contingency.jl` | N-1 contingency batch API (multi-core Phase 4): case generation (all in-service branches, transformer filter via component type or nonzero ratio, parallel-circuit disambiguation `name#branchIdx` with a shared-name fixture, FOR001 mapping incl. unresolvable names, #331 Phase 2 screening filters `min_vn_kV`/`min_sn_MVA`/`name_pattern` on a 380/110 kV fixture, and per-case weight validation), sp_case14 full N-1 with serial vs parallel field equality (`max_tasks=1` and auto) and base-net immutability, islanding semantics (a load-only island reports the specific `islanded: load-only, X MW` cause, a PV island is promoted and solves), non-convergence and unknown elements reported not thrown, the table printer plus CSV writer, and the #331 Phase 1 start-value ladder (`contingency.rescue_ladder`: `start_used` reported per case, ladder validation of the empty/duplicate/unknown-stage inputs at both the keyword and the config-section level, a longer ladder never losing a warm-converging case and staying serial-vs-parallel identical, the deprecated `retry_flat_start` alias, and a non-converged base case routed through the solver rescue ladder before the flat fallback, verified via a bad-seed fixture whose plain solve diverges and whose no-flat-fallback warning is asserted), and the #331 Phase 2 case weights (`applyContingencyWeights` by name with a default, `readContingencyWeightsCSV` for semicolon/comma/header/comment/negative rows, weight 0 as a pure ranking weight that still solves, and the weight carried into the result and shown in the table/CSV), and the #331 Phase 3 generator outages (`generateN1Generators` count/kind/Pg-name filters, a generator outage solving with the slack absorbing the loss, the slack-unit outage reported as "no slack bus" by default and rescued by `auto_slack`, `distributed_slack_enabled` flowing through, serial-vs-parallel identity, unknown-generator reporting, and a non-pegase stranded-generation fixture where removing a bus's PV unit leaves its island with injection but no reference), and the #331 Phase 4 overload reporting (structured `OverloadRecord` loadings with base loading and delta on a two-parallel-line overload fixture, `severity = weight * max(0, loading - 100)` with the result table severity-ranked failures-first, `shed_load_mw` as a structured field on a load-only island, and `buildContingencyReport`/`printContingencyReport` aggregating outcome counts, worst loading, worst severity, and total load shed), and the #331 Phase 5 Web UI service (`_run_contingency_service` for branch and generator kinds on the shipped sp_case14: succeeded status, `contingency_n1.csv` + `run.log` artifacts, `run_mode`/kind/counts metadata, the slack-unit outage named in the summary badge, and an invalid kind rejected not thrown) |
| `external_grid` | `test/test_external_grid.jl` | External-grid element (issue #299): the external-grid/distributed-slack mutual-exclusion configuration error (both cover the same imbalance; YAML pair rejected, single option loads), the duck-typing contract guard (`NativeShortCircuitData` field-identical to `CGMESShortCircuitData` by name, order, and type — machine-checked, not comment-checked), PF invariance (an `addExternalGrid!` slack with finite `Sk''` reproduces the manual slack path to machine precision — carried SC data changes no power-flow result, incl. the full ENI record contract), SC hand calculation on a single 110 kV bus (`Zq`, `Ik''`, `Sk''`, kappa/i\_p over two R/X values, c cancelling), min-case semantics (`sk_min` used unflagged with the rx\_min→rx\_max default, explicit `rx_min` wins, missing `sk_min` skips the feeder with the engine's `:no_source` flag), parallel-feeder stacking with unique, removal-safe mrids (suffix continues after the highest surviving id, incl. the same-prefix-bus guard), input validation (`ArgumentError` table), the net-copy regression (`sc_sources` survives `deepcopy`, SC on the copy matches), the `convertSlackToExternalGrid!` demotion (self-consistent bus types right after the demotion, conversion note contents, no-slack/non-slack-bus `ArgumentError`s), the connection statement in the classical result print (`Grid connection:` prose line naming slack vs. source incl. `Sk''`/R-X and the internal slack bus; type column `SOURCE` instead of `SLACK` on the internal bus, also in the structured report rows), and the `internal_impedance` variant (slack moves to the tagged auxiliary bus, stiff `Sk''→∞` limit reproduces the ideal solution below 1e-8 pu, realistic `Sk''` droops voltage and shifts angles, the feeder record stays anchored at the physical connection bus, one internal-impedance grid per bus) |
| `controls` | `test/test_voltage_dependent_control.jl`, `test/test_transformer_phase_shift.jl`, `test/test_tap_controller.jl`, `test/test_series_reactance_control.jl`, `test/test_upfc_control.jl`, `test/test_hvdc_pair_control.jl`, `test/test_tap_changer_model.jl`, `test/test_phase_tap_changer_model.jl`, `test/test_phase_tap_table.jl` | Voltage-dependent controls, transformer phase-shift control, tap-controller behavior, YAML controller instantiation (#305: round-trip via the real YAML reader against the programmatic twin, load-time and apply-time validation, idempotent re-apply), machine remote voltage control (`MachineVoltageControl`: API validation, secant convergence onto a remote target, honest `at_limit` on the reactive bounds), the STATCOM current-limit mode (#297 Draft A: registration validation incl. rating exclusivity and the `i_max_ka` conversion, at-limit `Q = V * S_max` tracking with live bounds across operating points, deadband-level in-range equivalence with the constant-Q mode, inert mode fields on constant-Q controllers, `STATCOM (VSC)` element vocabulary) and the SVC-vs-STATCOM limit contrast on one sagging corridor (#297 Draft E acceptance: delivered Q asserted against `V^2 * B` versus `V * S_max` numerically), the FACTS config surface (`s_max_mva` STATCOM and `v_inj_max_pu` SSSC entries instantiate via `applyConfiguredControllers!`, re-apply skips both; regression: the series idempotency check crashed on wrong field names instead of skipping), successful baseline PF preservation when controls are disabled, the `AbstractTapChangerModel`/`convention` field plus `calcRatioTapCorrection`/`calcRatioTapRange` consolidation, the CGMES `PhaseTapChangerModel`/`calcPhaseTapFraction`/`calcPhaseTapAngleRatio`/`calcPhaseTapReactance` formulas plus the DTF importer migration onto them, and the tabular `TapTablePoint`/`kind = :tabular`/`calcPhaseTapTable` override path including the formula-vs-table round-trip regression; the PST X(α) coupling (winding-resolved model lookup, continuous-angle formula evaluation, nearest-row tabular mapping, accepted moves update `x_pu` to the model value at the final angle, devices without reactance data stay static, probe direction consistent with and without tracking incl. angle/reactance restore) and the cross-type warning when a tap controller already regulates a machine controller's target bus; the master/slave transformer group (#322: circulating-Q physics on misaligned parallel units, group stepping with aligned taps and near-equal reactive split under a group-sized deadband, follower exclusivity/self-follow/double-follow/mode validation, the same-target-bus warning, `followers` through the config surface, `(+N followers)` element label; the CGMES side lives in the importer suite: a patched MicroGrid copy with two RatioTapChangers on ONE TapChangerControl imports as one master with one follower); the SVC shunt voltage controller (API validation, in-range secant regulation onto the setpoint with the actuated shunt carrying the same susceptance, honest constant-B `at_limit` on an unreachable target, `clearShuntControllers!`); the MSC/MSR discrete bank mode (#324: registration validation incl. step sign and no-admissible-block range, start-value snapping onto the block grid, whole-block truncated stepping that parks with `status = :parked` before crossing the target without exhausting the outer budget, `at_limit` on the outermost block, `MSC/MSR (switched shunt bank)` element vocabulary with `discrete = true` and step fields on the report row, `step_mvar` through `applyConfiguredControllers!` incl. idempotent re-apply, and continuous-SVC regression) and the generic controllable-element view (`controllableElements` rows for tap, machine, and shunt controllers incl. `ControlRunResult.elements`); the TCSC series-reactance controller (`SeriesReactanceControl`: loop-network target tracking within the deadband, honest `at_limit` clamping on an unreachable target, bit-identical baseline with a disabled controller, element-row vocabulary `:series_x_pu`/`:branch_active_power`, registration validation incl. transformer rejection and the `eps_z` impedance guard, and the no-move deadband case; SSSC injected-voltage mode, #297 Draft F: limit-form exclusivity validation, converged operation inside the live current-dependent window `x_base ± v_inj_max/|I|`, pinned operation with the effective injected voltage at `v_inj_max`, `SSSC (VSC)` element vocabulary with live window bounds, mode fields on the report row, summary printer naming the SSSC limit, and TCSC-mode regression); the UPFC stationary quadrature composite (`addUpfcControl!`, #325: exact equivalence with the manually registered SSSC+STATCOM pair on a two-corridor fixture, both limit characteristics at their clamps, composite-level rejections and the series-side rollback leaving the net untouched, `UPFC series/shunt (VSC pair, ...)` result-row pairing, YAML type `upfc` incl. the double-apply no-op) and the full DC-link-coupled model (`model = :full`, #326: independent line P and Q reached simultaneously with the DC-link residual closed on a meshed sending-bus fixture, the series active injection `P_se` nonzero, `series_phase = :quadrature` reducing to a standalone SSSC with `P_se = 0`, the `UPFC (full, DC-link coupled)` element row, per-model registration validation, the YAML `model: full` round trip with the double-apply no-op, the low-current `z_add` guard keeping the branch finite at its base impedance on a near-dead line, and the negative-resistance branch left by the full model being rejected loudly by a subsequent IEC 60909 short circuit while the base network runs); the HVDC pair controller (`HvdcPairControl`, #297 Draft B: registration validation incl. slack and PV guards, exact pairing invariant on a two-island fixture, reversed transfer, rating clamp with honest `at_limit`, per-side voltage-target secant incl. unreachable-target clamping, bit-identical baseline with a disabled controller, MATPOWER `paired_control` mode with the `pf_injections` consistency anchor, YAML `hvdc_pair` round trip, element-row vocabulary `:hvdc_p_transfer_mw`/`:hvdc_transfer`, the grid-forming `mode = :island_feed` incl. mode-dependent validation, the island-draw mirror, the rating clamp, and the YAML `mode: island_feed` declaration; r0.9.9 additions: the persistent `HvdcLink` records (`net.hvdcLinks` fills for MATPOWER Stage-0/paired and `addHvdcLink!`, controller attach/detach round trip), the `HVDC Link Flows` result surface (table, header count, converter-loss line, `ACPFlowReport.hvdc_links`), and the one-reference-per-synchronous-island rule (multi-reference error naming the buses, setpoint pair as parallel PQ path next to an AC tie, `invalid_topology` for a demoted `island_feed` reference)) |

## Extended profile additions

The `extended` profile is extended-only: it does not run the fast profile first. Use `all` when a single invocation must execute both fast and extended suites.

Current extended-only groups are:

- `core_model_extended`
- `contingency_extended`
- `scenario_engine`
- `example_infra`
- `programmatic_api_extended`
- `webui`
- `webui_extended`
- `repository_hygiene`
- `apslf`
- `cgmes_importer`
- `cgmes_export`
- `dtf_extended`

| Extended addition | File | Main checks |
|---|---|---|
| `webui` | `test/test_webui.jl` | Focused Web UI coverage for form parsing and backend validation, result rendering, active and terminal timing cards, commit-span omission, tolerance-step hook, path traversal rejection, DTF upload role classification, primary-case and FOR002 selector filtering, the CGMES-export run option (form checkbox + help topic, `export_cgmes` request flag incl. hidden-false/absent-field defaults, and the result-page summary row for completed/failed/absent export metadata), the last-edit-wins precedence between the configuration file and saved case settings (a newer YAML wins for its own keys and sets the notice flag, an older one keeps the sidecar values), the SE phase-5 chain (measurement-set upload classification via content sniff, the SE section of the Runs page (`/stateestimation` is a redirect onto its anchor) and demo generator, the service-level SE run with artifacts and history kind `se`, the SE-started PF and N-1 with the SE run id in the metadata and result parity against a manual run), the measurement generator v2 chain (truth state fresh-solve/from-run with bit-exact state adoption and the rejection paths, deterministic balance-aware single flow ends, passive nodes at a configured sigma or as protected ZI constraints, the se_deltas.csv artifact, and the bad-data threshold surface with staged/legacy service parity and invalid-mode rejection; plus the Delta-u PST chain on the tracked warmup case: generator targets the additional-voltage stepper, the mass release estimates it :pst along the nameplate psi, the fixation lands on the injected step, and the machine trafo stays calculated), the buildSysimage one-call dry run, the sysimage launcher decision from tools/sysimage_launcher.jl (image/metadata present, Julia version, Manifest hash, a src file newer than the image, and the no-terminal default) together with a parity check that `Sparlectra.webui_sysimage_problem` returns the identical verdict for every one of those states (the two implementations cannot share code, because the launcher runs before the package is loaded), the sysimage build-progress contract (a build that stopped reporting for more than a minute counts as gone, so a crashed builder cannot lock the refresh button) and the Sysimage page rendering idle, running, and failed including its build-log tail, a parse check of the generated app CLI driven from the checkout tool, and stubbed route checks without a real asynchronous solver run. Real asynchronous PowerFlow job lifecycles, artifact preview/download/ZIP/history/delete matrices, browser-launcher platform matrices, socket/server lifecycle checks, Markdown help/documentation cross-product validation, repeated real MATPOWER runs, and the scenario editor flows (scenario task step 6: a three-op scenario created through the form on the shipped sp_case60 bundle (staged copy switched to explicit-list mode), saved, reloaded and run with `:flag` showing the row and the screening columns; the result page's N-1 table showing a screened row with its estimate detail after `:flag`, no screened column after `:off`, and a failed case's error text in the row; a `tap_pos` op on a regulated transformer rejected with the controller named and NOT written to the file; a MATPOWER case getting the export hint with no editor, no `file_block` option on the main form, and an SCF case getting both) are extended-only in `webui_extended`. |
| `contingency_extended` | `test/test_contingency.jl` | N-1 identity slice on the real case1354pegase (first 60 branches; runs only when `SPARLECTRA_LARGE_CASES_DIR` points at the case, otherwise an explicit SKIPPED line, never a counting pass): serial vs parallel field equality with `max_tasks=1` and auto, plus the `retry_flat_start` invariance on converged cases |
| `scenario_engine` | `test/test_scenarios.jl` | Scenario task step 3, the engine gate: the N-1 result CSV on case118 (all branches, all generators, 240 cases through `runContingencies!` on the reused-worker engine) is byte identical to `test/fixtures/contingency_case118_n1_7438b6e.csv`, the CSV the pre-engine per-case-deepcopy implementation produced at commit 7438b6e (SKIPPED line when the untracked `data/mpower/case118.m` is absent); step 4, the screening acceptance on the same batch: with `screening_mode = :flag` and margin 10 no case the full run reports as violating (or failed) is screened away (false-positive count and full-run share are printed, not gated); step 4b, the same acceptance with distributed slack enabled (augmented lambda system, renormalized participation on a generator outage, screened share and seconds printed); the case300 acceptance per the test-network rule (operated grid, no false negatives at margin 10 with the 0.005 pu trust gate; loud SKIPPED without the local case file); the case1354pegase run measures runtime SCALING only (`:off`/`:flag`/`:only` seconds, gated on `SPARLECTRA_LARGE_CASES_DIR`, no screening-quality claim, see the test-network rule); step 7, the scenarios workshop (`docs/lit/workshop_scenarios.jl`) is EXECUTED here in a fresh module, so the `@assert` next to every number it prints gates the notebook against silent drift (loud SKIPPED without the local case118); step 5, the service sources: an SCF case with its own scenarios block runs via `scenario_source = file_block` (screening metadata plus the two CSV columns), an external scenario JSON via `external_file` with `screening_mode = off` keeping the classic columns, `n1_generators` records itself as the case source, unknown sources or a missing scenario file reject as `invalid_request`, and a CGMES delivery with an id-addressed scenario source is rejected with the way out named (export as SCF once, then reference its component ids) |
| `example_infra` | `test/test_example_suites.jl` | Example suite runner infrastructure: CLI parsing, skip and filter logic, CSV and Markdown escaping, registry integrity, and a --help smoke test of the generated runners |
| `matpower_examples` | `test/test_matpower_example.jl` | MATPOWER example runner path, output routing, performance/profile rendering, and runtime configuration forwarding, removed start-voltage alias rejection, and canonical profile-blend parsing |
| `synthetic_grids` | `test/test_synthetic_grids.jl` | Synthetic network generation and larger synthetic-grid regression coverage |
| `configuration_docs` | `test/test_configuration_docs.jl` | Configuration documentation and docs/config consistency checks |
| `cgmes_importer` | `test/test_cgmes_importer.jl` | CGMES importer: generic reader semantics on a synthetic in-memory fixture (profiles, overlays, references, DifferenceModel skip), the base-voltage inference (`infer_base_voltages` reconstructs stripped catalogs from the SV state, one summary warning, net still solves; without the option the typed abort stays) and the self-loop line guard (both terminals on one topological node map to a shunt notice, not a branch), the import-failure analysis (`importFailureAnalysis` names supplied models, missing vs. satisfied `md:Model.DependentOn` prerequisites, and the verdict; the `import_analysis_mode` service run succeeds on an importable delivery, fails with `import_analysis_not_importable` on a missing declared dependency with the `import_analysis.txt` artifact written, rejects non-CGMES cases, and excludes the other modes), `ReactiveCapabilityCurve` interpolation, remote-regulating machines on a synthetic fixture (default held-PV fallback, `machine_control = true` controller wiring + control-loop convergence, and the voltage-held-target fallback), test-set alias/border bookkeeping, placeholder guards + `AsynchronousMachine` Stage-0 mapping (FullGrid warns about its X.99 shunt/tap-row fillers and solves from a flat start; MiniGrid's three motors close the former SV gap, asserted SV-tight), and, gated on the local ENTSO-E cache with explicit skip messages, MicroGrid/Assembled import + `runpf!` + `compareWithSV` validation, Stage-2 tap controllers, PSEI PST fixtures (tabular table wiring and the documented end-2 angle-sign deviation), the RealGrid tabular-PST SV-start-mismatch regression, the CGMES fixed-reference self-check (SV voltages reach the solver verbatim even against a `flatstart: true` base config, `self_check.log`/`self_check_residuals.csv` artifacts, `cgmes.log` no-SV bus count), `cgmes_import.start_values` selection (flat/sv decision line in `run.log`/`cgmes.log`, precedence over a hostile `power_flow.flatstart`, mandatory `sv_compare.csv`/`sv_compare_flows.csv` artifacts + summary metadata in both modes, MATPOWER negative case with no effect/artifact/metadata), the `runShortCircuit!` MicroGrid plausibility sweep (finite positive Ik'' on every bus, max ≥ min per bus, no motor/defaulted-data flags, consistent with `shortCircuitCoverage`), the `short_circuit_mode` service run (CSV artifacts with the short-circuit schema, coverage report in `run.log`, MATPOWER rejected with `short_circuit_requires_cgmes`, mutual exclusion with `diagnose_mode`, real-delivery data-check), and `normalPF` → `participationFactor` arrival (MicroGrid zero-factor mapping, MiniGrid positive factor); the FullGrid placeholder-guard block additionally asserts `cgmes_import.placeholder_guards = strict` aborts the import with the offending object named; the node-breaker topology processor (#314, cache-gated): MiniGrid/SmallGrid/FullGrid imports without their TP file reproduce the TP partition (equal bus/branch/link counts, both sides flat-started, sorted-Vm multisets identical to 1e-8 for MiniGrid/SmallGrid, FullGrid within 0.05 for the documented one-CN-two-TNs data finding), the TP path stays processor-free, the T1 dead-fragment tolerance (branch counts equal, TP Vm multiset contained in the derived one), and the no-topology abort (bus-branch EQ without TP raises `CGMESImportError` naming the missing topology information, boundary TNs not counting) |
| `cgmes_export` | `test/test_cgmes_export.jl` | CGMES export identity and roundtrip contract (`writeCGMESFiles`, EQ+TP+SSH+SV + optional delivery zip): structural keys in `net.cgmes_ids` with lexicographic bus-pair normalization and first-seen parallel-line numbering, deterministic uuid5 minting (id equals uuid5 over the key), byte-identical re-export and independent-build determinism with a pinned `created` stamp, renaming a component never changes its mRID, the duplicate-mRID guard aborting before any file is written (both offending keys named), the power-flow-identical self-roundtrip over every exported class (transformer with off-nominal ratio + 2° phase shift + ratio-tap machinery, PV machine, slack injection, SVC with rating-derived Q limits, load, shunt, bus link → export → `importCGMES` → `runpf` matches the original within 1e-6 and starts from the exported SV in ≤ 2 iterations, empty notice list), tool provenance in every file (`Generated by Sparlectra.jl v<version>` comment + `md:Model.description`, stamp = pinned `created`), SSH content (machine sign convention, `referencePriority` on the slack unit, `RegulatingControl.targetValue` in kV, shunt sections, `bPerSection`) and SV validation (`compareWithSV` of the re-import: dvm/dva < 1e-12, flow rows < 1e-9), the regulated tap group (master + follower export ONE shared TapChangerControl with controlEnabled on both, reimport regroups them with target and deadband intact), zip re-import, and — gated on the local ENTSO-E cache with an explicit skip message — the MicroGrid roundtrip proving imported TN/ACL/PT/PT3/EC/SM/SH mRIDs survive an export (canonical form without the RDF underscore prefix) even after renaming every line, the 3W star reassembled into one three-end transformer (star TN absent from the export), the re-import rebuilding the identical electrical model (every branch parameter equal incl. the PST angle), the short-circuit evaluation of the re-imported delivery reproducing the original Ik'' rows and flags for both cases, `cgmesLineShortCircuitData` feeding harvested zero-sequence line attributes into the EQ profile, and the `export_cgmes` service-run path (exactly one delivery-zip artifact, no loose profile files, the zip re-imports with the original mRIDs, `cgmes_export_*` metadata incl. notices, `run.log` line) |
| `apslf` | `test/test_apslf.jl` | APSLF (AnalyticLoadFlow.jl) integration: `power_flow.solver`/`apslf`/`apslf_start` config validation (including the `solver=apslf` + `apslf_start.enabled` conflict), outer-loop-controller rejection with `solver=apslf`, Web UI form parsing of the new fields into the effective config, the adapter mapping (`PFModel` -> AnalyticLoadFlow spec, PF ordering), the standalone sp_case5 convergence against NR, and the start-value-generator guard including `nr_polish=false`. Nothing here is conditional: AnalyticLoadFlow.jl is a required dependency since 0.10.0, so the whole surface runs in every profile that includes the file. |

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
dispatches to library modules under `examples/dtf/` (`dtf_validation_base.jl`,
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

Each `examples/dtf/dtf_validation_*.jl` module is also directly runnable
as its own single-purpose CLI entry point (for example
`julia --project=. examples/dtf/dtf_validation_base.jl --dtf-file=... --for002-file=...`),
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
New tests use the SHIPPED demo cases (`data/scf/sp_case5`, `sp_case14`, `sp_case60`, `sp_case188`, see [Shipped Demo Cases](demo_cases.md)) instead of downloaded ones; a download belongs only in tests of the download and cache machinery itself (maintainer rule 2026-09-03). The load_fixture_net step (2026-09-04) completed that migration: `load_fixture_net` serves the shipped sp_ cases (legacy case9/case14/case57 resolve cache-only for externally anchored oracle tables, whose testsets gate on `fixture_net_available` and speak their skips), and a fresh installation downloads nothing. Measured with the MATPOWER cache set aside: the fast profile drops from 6873 to 6829 assertions (4 spoken skips: the case9/case14 DC oracle tables, the real-MATPOWER round-trip legs, the DTF stamping leg) and the extended profile from 5580 to 5567 (4 additional spoken skips: the case118 byte gate, the case300 screening anchor, the scenarios workshop, the case14 self-check anchor), with no file downloaded, so CI without a cache covers all but 57 assertions of the full local run. The state-estimation suite additionally guards the sparse core: the FD column coloring must produce a bit-identical Jacobian (colptr/rowval/nzval) against the per-column assembly on sp_case60, sp_case188 and the cache-gated case1354pegase, the structural coupling superset must contain the numeric pattern at two operating points per case, and the criticality budget and Takahashi LDLt second attempt carry their own regressions.
The Q-limit large-case comparison workflow is a manual diagnostic tool: it resolves optional large MATPOWER cases through Sparlectra's case registry/download support, compares `case × start_profile × qlimit_mode`, and uses its CSV/JSON summaries as the primary output. Real runs for cases such as `case13659pegase.m` and `case_SyntheticUSA.m` are intentionally excluded from the fast tests; automated coverage uses stubs for case resolution and API execution.
The experimental large-case Q-limit comparison test block is suppressed from the normal fast profile; run `test/test_qlimit_large_case_comparison.jl` manually when maintaining that tool.

Large measurement networks are picked by what a test actually judges (rule 2026-09-03). The pegase cases (`case1354pegase`, `case13659pegase`) are convergence and scaling instances from the OPF world: their base states carry overloads and non-converging outages by construction (case1354's base holds a branch above 108 percent, and 56 of its 2251 N-1 cases do not converge). They stay the oracle for convergence, import conventions (rad, shift sign -1.0), start ladders, Q-limit switching, and runtime/memory scaling, and nothing else. Anything that presupposes an operated grid, screening quality, N-1 loading margins, voltage bands, or scenario evaluation, is judged on RealGrid (ENTSO-E CGMES, 6051 buses, local cache, imported through the CGMES adapter) and `case300` as the small local MATPOWER case; ACTIVSg2000 joins when it becomes available locally. A screening run that flags every pegase case says nothing about the screening.

Default fast-profile output is intentionally compact: the runner prints the selected profile, one line for the include phase, one `[n/8]` marker per group, and Julia's final test summary. Each group's `PASS` line carries wall time, **compile and recompile time** (from `@timed`, Julia 1.11 and newer), allocations and GC time, so a slow group can be told apart from a group that merely compiled a lot: the two need opposite remedies. The include line measures what the groups cannot, namely parsing and compiling the testset bodies before any group runs. MATPOWER import diagnostics, auto-profile tables, runtime casefile banners, Q-limit tables, and similar artifact-oriented diagnostic blocks are suppressed in normal test stdout so progress remains scannable. Two exceptions surface immediately instead of being captured away. Availability-gated skips: captured lines containing `SKIPPED` are re-printed under the group's `PASS` line, so a gated testset (startup-hint under a sysimage session, the pegase slice without `SPARLECTRA_LARGE_CASES_DIR`) never skips invisibly and a lower-than-documented assertion count has its explanation right in the log. And test failures: a plain `@test` failure throws nothing until the outer suite aggregates, so the runner scans each group's capture for `Test Failed at` / `Error During Test at` blocks and replays them to stderr right under the group line (capped per block and per group), which spares the `--verbose` rerun that diagnosing a swallowed failure used to need.

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

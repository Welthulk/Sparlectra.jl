# Test Suite

`test/runtests.jl` loads tests by profile. Profile selection precedence:
1. CLI argument (`julia --project=. test/runtests.jl <profile>`)
2. `SPARLECTRA_TEST_PROFILE`
3. default `fast`

`Pkg.test()` uses the same runner and therefore the `fast` profile unless
`SPARLECTRA_TEST_PROFILE` is set:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Test profiles

Every test group belongs to exactly one of eight base profiles; `extended`
is the union of the five that are neither `fast`, `install` nor
`workshops`, and `all` runs everything. The documentation build is a gate of
its own (`sh tools/run_gates.sh docs`), not a profile: a change to Markdown
or a docstring needs the docs build and nothing else.

| Profile | Command | Scope | When to run |
|---|---|---|---|
| `fast` (default) | `julia --project=. test/runtests.jl fast` | Network model, terminal status, the rectangular Newton core and its linear solver, DC power flow, distributed slack, islands, external grid, MATPOWER metadata, the small API smoke | Every pull request and every commit; CI runs this one |
| `pf` | `julia --project=. test/runtests.jl pf` | Auto mode, short circuit, parallel foundation, contingencies, scenarios, all controllers, the extended grid and contingency sets, the scenario engine, APSLF | After a change to the power-flow path, its controllers or the scenario engine |
| `se` | `julia --project=. test/runtests.jl se` | State estimation, observability, topology validation | After a change to the estimator or the measurement model |
| `config` | `julia --project=. test/runtests.jl config` | Configuration surface and validation, configuration documentation coverage, repository hygiene | After a change to a configuration key, a default or the documentation of one |
| `webui` | `julia --project=. test/runtests.jl webui` | The local browser UI, both files | Before every Web UI pull request; `fast` proves nothing about `src/webui` |
| `extd` | `julia --project=. test/runtests.jl extd` | Shipped demo cases, SCF, service lifecycle, MATPOWER examples, example infrastructure, net cache, synthetic grids | After a change to a format, an importer, the service layer or the examples |
| `adapters` | `julia --project=. test/runtests.jl adapters` | The format importers with their own fixtures: CGMES import and export, DTF, PowSyBl | After a change to an importer, an adapter, a fixture or the shared network helpers they call |
| `install` | `julia --project=. test/runtests.jl install` | The installation path: a copy of the checkout without the application manifest, `start_webui.jl --env-only` twice, the first start sets up and compiles, the second does nothing; plus the environment-dependent PowSyBl extension check | Before a release; about a minute of cold compile, so not part of `extended` |
| `workshops` | `julia --project=. test/runtests.jl workshops` | Every Literate workshop under `docs/lit` top to bottom in a fresh module, with the assert next to each printed number | After a change to a workshop or to an API a workshop uses; before the notebooks are regenerated |
| `extended` | `julia --project=. test/runtests.jl extended` | `pf`, `se`, `config`, `webui`, `extd` and `adapters` in that order | Before a merge |
| `all` | `julia --project=. test/runtests.jl all` | `fast`, `extended`, then `install` and `workshops` | Declaring a branch merge-ready, together with the docs build |
The service layer and the Web UI are the `SparlectraApp` package under
`app/`; the runner puts that directory on the load path and loads the
package for every profile, so `webui`, the service part of `extd` and the
application smoke test of `fast` (one power flow and one state estimation
through the service API) run from the library checkout with the commands
above.

Group membership (the table `TEST_PROFILES` in `test/runtests.jl` is the
source; the runner refuses a group in two profiles or in none):

| Profile | Groups |
|---|---|
| `fast` | `core_model`, `terminal_status`, `powerflow_rectangular`, `factorized_linear_solver`, `pv_voltage_residuals`, `3wt_phase_taps`, `dc_powerflow`, `distributed_slack`, `island_diagnostics`, `external_grid`, `matpower_metadata`, `programmatic_api` |
| `pf` | `auto_powerflow`, `short_circuit`, `parallel_foundation`, `contingency`, `scenarios`, `controls`, `core_model_extended`, `contingency_extended`, `scenario_engine`, `apslf` |
| `se` | `state_estimation`, `observability`, `topology_validation` |
| `config` | `configuration`, `configuration_docs`, `repository_hygiene` |
| `webui` | `webui`, `webui_extended` |
| `extd` | `demo_cases`, `scf`, `programmatic_api_extended`, `matpower_examples`, `example_infra`, `net_cache`, `synthetic_grids` |
| `adapters` | `cgmes_importer`, `cgmes_export`, `dtf_extended`, `powsybl_importer` |
| `install` | `install` |
| `workshops` | `workshops` |
## Runner switches and output

- `SPARLECTRA_TEST_SKIP_GROUPS` (comma-separated group names) drops named
  groups from any profile and prints every skipped group at the start of the
  run. It exists for the app build workload
  (`tools/app_workload_runtime.jl`); do not use it to get a green
  development or CI run.
- `SPARLECTRA_TEST_SHOW_WARNINGS=1` lets `@warn` records from Sparlectra
  modules reach the console. By default the runner drops them: the tests
  exercise deprecated configuration values, version-less YAML files,
  non-case keys handed to `write_case_config` and advisory topology
  prechecks on purpose. Errors, log records from other packages, and
  `@test_logs` assertions are not affected.
- `SPARLECTRA_TEST_VERBOSE=1` or `--verbose`
  (`julia --project=. test/runtests.jl fast --verbose`) shows the full test
  output; the profile selection is unchanged.
- `SPARLECTRA_TEST_GC_BETWEEN_GROUPS=1` requests a GC cycle after each
  completed group.
- `SPARLECTRA_TIMING_TESTS=1` enables the runtime-scaling legs (see
  `scenario_engine` below).
- `SPARLECTRA_SYSIMAGE_TRACE_TESTS=1` makes the sysimage workload trace the
  fast profile instead of its curated list of interactive paths (see
  [Sysimage](sysimage.md)), for a comparison only.
- `tools/run_gates.sh` starts Julia with `--startup-file=no`: a personal
  `startup.jl` that loads Revise invalidates hundreds of method instances of
  the system image, which the run then infers again. Use the same flag when
  running `test/runtests.jl` by hand on such a machine.

Default output is compact: the selected profile, one line for the include
phase (parsing and compiling the testset bodies), one `[n/8]` marker per
group, and Julia's final test summary. Each group's `PASS` line carries wall
time, compile and recompile time, allocations and GC time, so a slow group
can be told apart from a group that merely compiled a lot. Import
diagnostics, auto-profile tables, Q-limit tables and the rectangular/Q-limit
diagnostic blocks are suppressed (they remain available through solver
`verbose > 0`). Two things surface immediately: captured lines containing
`SKIPPED` are re-printed under the group's `PASS` line, so an
availability-gated testset never skips invisibly; and `Test Failed at` /
`Error During Test at` blocks are replayed to stderr under the group line.

```text
Test framework: fast
[1/8] core_model
[2/8] powerflow_rectangular
...
[8/8] controls
```

## Reference timings

Measured warm on the development machine (16 cores, Julia 1.13, package
images built with the full precompile workload, each profile in its own
process, the second run with the compile cache of both packages in place,
also for the child processes some tests start; the first run after a source
change does not count): fast 96 s, pf 101 s, se 48 s, config 22 s, webui
120 s, extd 149 s, install about 75 s. A profile that grows by more than ten
percent against these numbers needs a cause before the change is committed;
the numbers are updated here on purpose, never silently. A run right after
a source change can show more, because the sysimage dry-run test starts a
child process with its own compile cache.

The image matters as much as the code: on the default image (precompile
workload off, the install default) every profile compiles the solver paths
itself and takes roughly 1.8 times these numbers (webui test group 185 s
instead of 69 s, measured on one commit). `tools/run_gates.sh` therefore
rebuilds both images with `SPARLECTRA_PRECOMPILE_WORKLOAD=full` when the
loaded image was not built with it (the mode is baked into the image as
`PRECOMPILE_WORKLOAD_MODE`; an image without the constant counts as
unknown and is rebuilt), announces the rebuild with its duration, and a
comparison against the reference is only valid on such a run. Measured
this way on 2026-09-26: fast 77 s, extended (every area profile in one
process) 315 s. CI runs the fast profile on the default image, so the
install path (no workload, every path compiled on first use) keeps a gate
of its own.

Per-group seconds are dominated by Julia compilation, and which group pays
depends on run order: the first group that touches a code path compiles it
for everyone after. Judge a "group got slower" observation against warm
times, never against a single group's cold seconds. Two conventions keep the
compile share down and are expected of new tests and new service code:

- Every `@testset "..." begin ... end` body is wrapped as
  `begin (function () ... end)() end`. Julia infers one function body at a
  time, and the cost grows faster than the length. A testset body that
  contains a `return` of its own scope, a `const`, a type definition,
  `using`/`import` or `@eval` cannot be wrapped and stays as it is.
- A long function that takes loosely typed arguments (a callback, a
  `nothing`-or-value setting, an `AbstractDict`) is compiled once per
  distinct argument type combination. The API entry, the state-estimation
  service and the configuration context take a thin keyword wrapper that
  puts those values behind one concrete carrier struct (or a
  `@nospecialize`); the options struct in `src/api/run_api.jl` is the
  pattern.

## Fast versus extended ownership

"Extended" here means the six profiles that are not `fast` (`pf`, `se`, `config`, `webui`, `extd`, `adapters`).

| Area | Fast coverage | Extended coverage |
|---|---|---|
| Core/model | Small constructors, transformer checks, bus/prosumer semantics, link behavior, representative rectangular PF/Q-limit regressions, and small sparse fallbacks. Large sparse Ybus and large MATPOWER matrix-body checks are extended-only. | Large sparse Ybus smoke, large MATPOWER matrix-body scanner, synthetic/stress grids, and longer integration examples. |
| API | Serialization and transport helpers, path and validation safety, one successful small API run, one pre-solver failure, one numerical/island failure, Solver-time and Total-time contracts, critical DC-line default and strict-rejection smokes, and one small independent-island regression. | Exhaustive CSV/export matrices, repeated artifact inventories, island artifact-content matrices, persistent history/delete/reload lifecycle coverage, and repeated presentation/performance-log modes. |
| Web UI | Form parsing and backend validation, result rendering, active and terminal timing cards, commit-span omission, tolerance-step hook, path traversal rejection, DTF upload role classification, primary-case and FOR002 selector filtering, the buildSysimage one-call dry run, the sysimage launcher decision (all four staleness reasons including a src file newer than the image, and the no-terminal default that does not build), the agreement between that launcher verdict and `Sparlectra.webui_sysimage_problem` across every fixture state, the sysimage build-progress reader (running, gone-stale, done, failed) and the Sysimage page in each of those states including its rebuild guard, a parse check of the generated app CLI (run/se/n1 commands, config flags) driven from the checkout tool, and stubbed route checks without a real asynchronous solver run. | Case-profile persistence with asynchronous jobs, real run/result polling, artifact preview/download/ZIP/history/delete lifecycle checks, browser-launcher matrices, socket/server lifecycle, and Markdown/help/documentation cross-products. |
| Documentation/hygiene | No repository-wide documentation/help scan in fast; only focused source-level smoke checks tied to edited paths. | Configuration documentation consistency and normalized tracked-path/content repository hygiene scans, plus the reference-page coverage assertion (every src/ Julia file on exactly one reference page or on the exclusion list), a scan that no docstring is separated from its definition by a blank line (Julia then drops it silently, and only the documentation build notices, by failing on an unresolved `@ref`), plus additional check files next to the hygiene test that are included when present and named in the group report. |
## Case files and fixtures

Large MATPOWER and DTF cases are not part of the repository. Tests look for
them in the directory resolved by `Sparlectra.large_cases_dir()`:

1. `SPARLECTRA_LARGE_CASES_DIR`, if set (override for CI and special setups)
2. otherwise the Web UI user case directory
   (`~/.local/state/sparlectra/webui/data/mpower` on Linux,
   `%LOCALAPPDATA%\Sparlectra\WebUI\data\mpower` on Windows)

Test groups that need a large case check availability and report `SKIPPED`
when the file is absent; the fast profile does not depend on any large case.
Tests only read from this directory; generated files go to `mktempdir`.
`data/mpower/` in the repository holds fixtures only, and the gate check
fails on untracked files there. The assertion count depends on that
directory, so a measurement has to say which count it was taken with; the
`SKIPPED` lines say which run happened.

New tests use the shipped demo cases (`data/scf/sp_case5`, `sp_case14`,
`sp_case60`, `sp_case188`, see [Shipped Demo Cases](demo_cases.md)) and the
synthetic MATPOWER cases under `data/mpower` (`sp_case9`, `sp_case118`,
`sp_case300`, `sp_case1354`); `load_fixture_net` serves them. A download
belongs only in tests of the download and cache machinery itself; a fresh
installation downloads nothing. What still gates on the case cache says so
with a `SKIPPED` line: the SCF round-trip legs on real MATPOWER files
(numeric bus names, unlimited ratings, a neutral tap outside its band), the
case9/case14 DC oracle tables, the case118 byte gate, the case300 screening
anchor, and the DTF deliveries.

No test downloads a CGMES delivery or reads the cache under `data/CGMES`;
the ENTSO-E conformity packages and the ReliCapGrid models are third-party
data that is not checked in. CGMES coverage runs on deliveries exported by
Sparlectra itself: `tools/gen_cgmes_fixtures.jl` imports sp_case14 (OLTC
voltage controller), sp_case118 (54 machines with asymmetric Q limits and
synchronous condensers) and sp_casePST (phase-shifting transformer, bus
link), solves them with Q-limits off and writes EQ, TP, SSH and SV under
`data/cgmes_demo/<case>/` with a fixed header stamp, so a regeneration is
byte identical; the SV profile is the solution of the delivery itself. Not
covered, because the exporter cannot produce it: node-breaker topology,
boundary sets and X-node equivalents, DifferenceModel files, CGMES 3.0,
HVDC converters, machine short-circuit attributes, and the quirks of real
deliveries.

Large measurement networks are picked by what a test judges. The pegase
cases (`case1354pegase`, `case13659pegase`) are OPF instances whose base
states carry overloads and non-converging outages by construction; they are
the oracle for convergence, import conventions (rad, shift sign -1.0), start
ladders, Q-limit switching and runtime/memory scaling, and nothing else.
Anything that presupposes an operated grid (screening quality, N-1 loading
margins, voltage bands, scenario evaluation) is judged on RealGrid (ENTSO-E
CGMES, local cache, via the CGMES adapter) and `case300` as the small local
MATPOWER case.

The Q-limit large-case comparison (`test/test_qlimit_large_case_comparison.jl`)
is a manual diagnostic tool outside every profile: it resolves optional
large MATPOWER cases through the case registry and compares
`case × start_profile × qlimit_mode` with CSV/JSON summaries. Run it by hand
when maintaining that tool.

## Group reference: model, solver, estimator and format groups

| Group | Files | Main checks |
|---|---|---|
| `core_model` | `test/testgrid.jl` | The active-set voltage-side release (the shipped Zeng/Chiang case reaches the physical solution with at least one PQ->PV release and a clean Q-V check, the classic one-at-a-time mode ending on a physical point as well (no bound on the voltage difference between the modes: switching variants need not coincide), while a margin of 1.0 reproduces the old state with zero releases and a Q-V finding; the feeder machine clamped below its setpoint stays clamped; the shipped case118 keeps its clamped set and iteration count with and without the rule; the unit test on `active_set_q_limits!` covers both sides, the margin and the Q band fallback without voltages; the final Q-limit check with its two thresholds: every class on synthetic injections, the bounded warning, the strict 0/0 regression, the configuration default and validation, the same verdict for the three enforcement modes on sp_case14, and the classic run's within-hysteresis row at bus 6 of the Zeng case); core net construction and validation, small inline MATPOWER import helpers, file-based MATPOWER projected-start normalization, PV/PQ lock-ID resolution, Q-limit enforcement-mode parsing and classical base-failure/no-reenable dispatch checks, rectangular nonfinite mismatch and status-row diagnostic preservation, Jacobian condition-number estimator checks (exact diagonal case, Hager bound, near-singular flagging, input rejection), link handling, shunts, reporting/output checks, and summary-file output regression. Large sparse Ybus checks and large MATPOWER matrix-body scanner coverage are extended-only. Since 0.20.0 `test_piline_g.jl` also covers the per-terminal shunt split: an explicit split on a line and a transformer stamps each arm on its own end (the from arm through the ratio squared), the flow helpers use the arm of the end they leave, the two setters keep the totals as the sums, a partial keyword set is refused, and no code path outside `branch.jl` writes a branch total alone (grep over `src`). |
| `terminal_status` | `test/test_terminal_status.jl` | Per-terminal branch status (one-sided open branches): the dangling-node reference anchor (a `:open_to`/`:open_from` reduction matches a solve with an explicit zero-injection auxiliary bus to 1e-10, for a line and for a transformer with off-nominal ratio and phase shift), bitwise equivalence of both-flags-open with `status = 0` in the Y-bus, solver toggles between runs (sparsity-pattern invalidation), isolation of the open-side bus in the island report, the MATPOWER export substitution (BR_STATUS 0 plus exact `Y_in` bus shunt, `open_terminal=` marker, voltage-identical reimport), the result surface (`open@to` marker, `Open terminals` header count, `terminal_state`/open-end voltage report columns), and the state-estimation guard rejecting flow measurements on partially open branches |
| `powerflow_rectangular` | `test/test_solver_interface.jl` | Rectangular power-flow API behavior, sparse-only solver path, AC-island detection/reference validation/independent solving including aggregate all-island convergence status and iteration accounting, Q-limit and typed configuration entry checks, current-iteration start rejection diagnostics on a tiny synthetic case, including accepted and rejected framework control-status composition, ordered local MATPOWER batch execution, the thin `run_acpflow` alias, legacy-keyword rejection, the optional merit-function Armijo line search (value/scaling, acceptance/fallback/active-set-skip reasons, config validation, merit.enabled=false regression, and a converging flatstart/Q-limit-switching functional case), the optional scaled-Newton trust-region step control (config validation, acceptance/collapse unit coverage, trust_region.enabled=false regression, a converging functional case with log/diagnostics assertions, and radius-collapse non-convergence reporting), and wrong-branch detection output visibility across ACPFlowReport metadata, the AC island diagnostics CSV, and the console summary line, including the highest-voltage-level scope regression (a dip below the top level never triggers; the same dip on the top level does), and the AC rescue ladder / DC fallback (`power_flow.rescue` recovers a poisoned start and records the winning strategy, `power_flow.dc.fallback` leaves a DC state with the AC status honestly non-converged, defaults leave both mechanisms off) |
| `factorized_linear_solver` | `test/test_factorized_linear_solver.jl` | umfpack_reuse linear-solver selection and factorization-reuse behavior for the rectangular power-flow Jacobian path, including rejection of the `klu` value |
| `pv_voltage_residuals` | `test/test_pv_voltage_residuals.jl` | PV-voltage residual behavior, angle-preserving voltage-setpoint starts, phase-shifted PV integration coverage, and related solver diagnostics |
| `3wt_phase_taps` | `test/test_3wt_phase_taps.jl` | `create3WTWindings!`'s `phase_tap_side`/`phase_taps` keywords: existing-behaviour snapshot, single-winding attachment, ratio+phase coexisting on one winding, and validation errors |
| `configuration` | `test/test_configuration_coverage.jl`, `test/test_runner_helpers.jl` | Configuration-key coverage, safe refresh checks including current-iteration start defaults, forwarding checks, start-voltage value-domain validation (including `cgmes_import.start_values` and `cgmes_import.placeholder_guards` parsing: defaults `flat`/`warn_skip`, `sv`/`strict` accepted, invalid values rejected with the key name), and test-runner output-mode helper checks |
| `matpower_metadata` | `test/test_matpower_metadata.jl` | MATPOWER parser metadata fields, legacy bus sorting of `bus_name`, opt-in bus-name import, branch-kind override, branch metadata retention, FOR001 contingency mapping, default/opt-in `mpc.dcline` behavior including PF/loss/Q/V/Q-limit mapping and voltage-control safeguards, `writeMatpowerCasefile` `write_solution` column count/VM-VA/marker behavior, the `mpc.sparlectra.tap_changer_model` roundtrip marker preventing double tap-impedance correction, the `mpc.sparlectra.links` busbar-coupler extension (parse, link contraction in the power flow, export roundtrip, malformed-block rejection), the `mpc.sparlectra.tap_changers` nameplate extension (neutral-vs-live mapping onto the branch tap fields, ratio-less phase shifters, the Delta-u additional-voltage stepper with cascade-derived live position and roundtrip, exclusive phase grids, generated-flow consistency, export roundtrip, rejection matrix), and `net.name` staying the case name (not a bus's original name) when `mpc.bus_name` metadata is present |
| `programmatic_api` | `test/test_api.jl` | Focused GUI-ready power-flow API coverage for serialization and transport helpers, path and validation safety (including run-index acceptance of symlinked output roots, skipped on Windows), one successful small API run, one pre-solver failure, one numerical/island failure, Solver-time and Total-time contracts, critical DC-line default and strict-rejection smokes, one small independent-island regression, and the net-parameter stamp guard (each importer stamps the configured switching parameters exactly once: MATPOWER and SCF on shipped cases, CGMES on the checked-in sp_casePST delivery under `data/cgmes_demo`, DTF when `data/DTF/FOR001.DAT` is present with a spoken skip otherwise). Exhaustive CSV/export matrices, full artifact inventories, persistent restart/history/delete lifecycle matrices, repeated performance-log modes, complete island artifact-content matrices, and the fixed-reference self-check contract (verbatim-start summary in `self_check.log`, per-bus `self_check_residuals.csv`, the sp_case14 residual property plus the cache-gated historical case14 anchor) are extended-only in `programmatic_api_extended`. |
| `programmatic_api_extended` | `test/test_api_extended.jl` | Extended coverage of `run_sparlectra_api` and the local PowerFlow service: case resolution, Q-limit metadata and the `q_limit.log` artifact, the artifact list with kinds and MIME types, the one CSV format for every artifact of a run, the diagnose self-check report, run deletion safety. Since 0.20.0 also the tap-changer model precedence artifact: an SCF case with a typed phase model run through `run_sparlectra_api` writes `tap_models.log` (kind `tap_model_log`, named in `run.log`), the same network without the model writes none. |
| `state_estimation` | `test/test_state_estimation.jl` | WLS state-estimation behavior (baseline convergence and voltage accuracy, measurement add helpers, passive-bus zero-injection handling, PMU angle measurements), the `ImagMeas` current-magnitude rules (prediction roundtrip, 3 sigma value gate, shunt-bay variant, observability invariance), bad-data localization stage 1 (`wii` exposure, sequential elimination with trace/stop reasons, ZIB/SHDERIV protection, optional K-matrix report), shunt estimation (case A recovery with freeze guards, write-back gate, voltage sensitivity; case B back-calculation with sign and sigma propagation), the phase-3 link rules (SE on the contracted net with member sync, LINKAGG aggregation and partial-cluster exclusion, link measurements never in the WLS, bitwise KCL regression of the W2 allocation plus measured ring split, `se_view` classification), the phase-4 hardening (robust R stages with solve/diagnosis separation incl. the generalized knees k1/k2 with bitwise 3/6 equivalence and continuity at t = k1, the replacement suppression judged on the converged state with statistics on original sigmas, the k_eliminate exposure, Wilson-Hilferty band test with `:high`/`:low` reasons, diagnostics state unification), the phase-5 chain (measurement CSV v1 roundtrip with atomic import incl. `branch_nr` parallel-branch disambiguation and the legacy header, `runpf_from_se!` start modes with slack pickup, untouched model loads, and the measurement/model-discrepancy case separating `:se_state` from `:se_snapshot`), island-wise SE (two-island partition with per-island references and band verdicts, unmeasured-island skip, closed-link island fusion, merged chain registration), the transformer tap estimation (PF-equivalence gate at 1e-12 for the unstamped cascade overlay in all three release modes plus the live-tap KCL consistency regression, generated measurements stay self-consistent when tap_ratio leaves neutral; ratio/PST/two-regulator roundtrips recovering injected mechanical steps with the mandatory fixation run incl. the off-grid half-step case and the `:offgrid_tap_residual` band-flip note; the release guards: machine-transformer detection, structural bridge freeze without a far-side voltage bitwise equal to no release, partial `:both` column re-packing with bitwise stamp restore; and the two-island aggregation with per-island rows, summed fixation J/dof, and bitwise cross-island isolation; the completion coverage: `updateTaps` write-back with bitwise protection, the PMU current-phasor synergy shrinking the phase-regulator spread, bad data next to a released trafo (elimination names the row, tap stays within half a step, diagnostics on the tap-unified net), the machine-trafo back-calculation `calcMachineTrafoTapFromSE` incl. mutual exclusion, and the `se_view` release listing), the tap fallback on non-convergence, the run configuration as an argument, and the DTF case through the SE service (the requested format travels with the request; without it the run fails naming `dtf_for001` instead of a bare "unknown format") |
| `observability` | `test/test_observability.jl` | The observability side of state estimation on nets built once per run (3-bus, sp_case14, the link net closed/open/with shunt): global and local observability metrics on full, reduced and sparse sets incl. `measurement_jacobian` (labeled H: row/column counts against the global check, flow rows named by oriented branch, bus rows by bus, every row touches a state, rank agrees with `evaluate_global_observability`), a numerical deficit never overruled by a full structural matching (thinned sp_case14 draws), the rank verdict independent of measurement scale and `rank_tol_factor` forwarded, null-space dark states (the column-restricted local check is necessary but not sufficient, the 7-bus spanning-tree set names exactly the dropped bus), the matrix-only helpers (rank, matching, redundancy, local column selection, error paths), the two-stage FD-aware rank tolerance on the linked net (closed link observable, open link as measured island partition, explicit `tol` wins, island-internal angle pair fails the local check), and the Takahashi diagnostics (dense/sparse path equivalence on three case classes, FD exact-zero sparsity, guard fallback with warning, K pinning the dense path, LDLt second attempt, bit-identical FD column coloring on the shipped cases and cache-gated on case1354pegase) with criticality from diag(Omega) (no budget, a 600-row identity has every row critical, the per-row rank tests stay as `criticality_method = rank` with their budget skip, and both methods agree on sp_case5 with a flow row removed) |
| `topology_validation` | `test/test_topology_validation.jl` | Topology validation (all advisory): the stage-1 precheck matrix (each finding kind fires exactly its synthetic case, a clean set fires none) incl. the false-positive guards (lightly loaded branch expected near zero and flat start state without expectation, sp_case14 shipped set clean, open link disagreeing freely, partially measured node skipped, shunt buses balanced with the measured voltage), the advisory contract (findings present AND the estimation still runs), the stage-2 fingerprint (model-closed/actually-open transformer: elimination exhausts, band `:high`, the finding names the right station; a curable 10-sigma error never classifies as topology), the stage-3 hypothesis test (`:auto` candidates, the true toggle ranked first and supported, honest ambiguity on the small net, the input bitwise untouched, the island case), and the singular-normal-equations regression (dense LS fallback with a raised `svd_max_n`, parallel released taps converging on a least-squares split instead of crashing) |
| `auto_powerflow` | Automatic power-flow mode: read-only feature extraction on four case classes, ordered-rule strategy selection with step-control exclusion, the escalation ladder (cap, skip reasons, L4 locks, L5 evidence gate, seed-only L6), DC-fallback labeling, precedence, decision-log artifact and metadata mirror, hint-generator ids, the startup latency hint, and the auto-profile precedence level (explicit non-default user keys survive every recommendation, applied and skipped pairs reported, resolve_config carries the level and explicit overrides beat it). |
| `scf` | Sparlectra Case Format writer: the one result-CSV writer (`write_result_csv` under technical and excel_de, the SE state CSV round trip in excel_de, the AC island report and every CSV of an excel_de power-flow run carrying the semicolon); deterministic ids with byte-identical repeat and rebuild exports, the PGM component subset, SI unit conversion against hand-computed references (line ohm/farad on the to-side base, live k/theta, appliance W/var), the namespaced block (roles, tap cascade with the Delta-u controller, names, config, opt-in start state), measurement-to-sensor pairing with complete row accounting, the FACTS base-impedance guard, the round-trip identity (byte-identical re-export, identical Y-bus, measurements, power-flow iterations and voltages, identical state-estimation objective and state), the permanent v1 fixture, reader validation (unknown keys, unsupported components, dangling references, duplicate ids, batch datasets, malformed JSON), the controller round trip through the declarative schema, the study definitions with their load-time validation, the study runs themselves (an N-1 case list and exclusions taken from the file, a rejected multi-outage case, a short-circuit run on the file's own sources with its case, c-factor and explicit bus sweep, and the missing-source path), the round trip on real MATPOWER cases (case14 and case57: byte identity, component names, bus order, tap ratios, shunt signs, Y-bus, unlimited ratings; cache-gated with a spoken skip since load_fixture_net, because their quirks have no shipped equivalent), the own-name-vs-external_id channel on the shipped sp_case5 (reference name applied, external_id restored, re-export preserves both), power-grid-model interoperability (a file with no namespaced block whose sources become the reference, its solved voltages checked against an independently derived model), the strict PGM writer mode, the three-winding nameplate form with its named types and its rejections, the shunt-state and PGM fault blocks, the case-file upload including a rejected foreign JSON, the download route with its refusals (directory escape, missing file), and the service case resolver, SCF as a first-class case format (detection, selector, framework and service runs, case-file configuration precedence), the logging surface in the configuration allowlist, and the Web UI export route including its failure paths. |
| `dc_powerflow` | `test/test_dc_powerflow.jl` | Standalone DC power flow (`rundcpf!`, `power_flow.solver=:dc`): reference-value comparisons against an independently-computed script (3-bus hand-verified fixture always; the MATPOWER case9/case14 oracle tables cache-gated with spoken skips since load_fixture_net) for bus angles and branch flows, the hand-checkable sp_case5 DC dispatch (slack exactly 32 minus 12 MW), phase-shifter (`Pfinj`) coverage via a synthetic phase-shifting-transformer fixture, lossless power-balance property checks on the tracked warmup case and the shipped sp_case5, two-island independent-solve coverage, `power_flow.dc.*`/`power_flow.solver` config validation, active-outer-loop-controller rejection (mirrors the `apslf` rejection test), `run_sparlectra` dispatch (`method=:dc`, `diagnostics.solver=:dc`), `angle_reference_deg` exact-uniform-shift behavior, the `dc_pf_status`/`rectangular_pf_status` registry separation, `rundcpf!(net; seed_ac_start=true)` DC-seeded AC Newton-Raphson convergence, and the one-specification-source contract (a node-sum edit moves neither solver, a prosumer edit moves both, DC injections equal the real part of `buildComplexSVec`) |
| `distributed_slack` | `test/test_distributed_slack.jl` | Distributed active-power slack: disabled-path bit-identity to the classical solver, per-bus residual conservation (`alpha .* lambda_P` pattern, REF/load buses untouched), REF-equivalence (`explicit` weight fully on the reference bus reproduces the classical solution and `lambda_P` equals the classically REF-absorbed power), mode agreement (`pg_weighted` vs matching `explicit`), MATPOWER `APF` column import (`:imported` shares proportional to APF; APF presence never changes disabled results), fallback behavior (`error` throws, `ref_only` warns and matches classical), Q-limit interaction (PV→PQ switch keeps participation), the config-load regression for the `weights: {}` placeholder plus block-style weight tables (the minimal YAML reader has no flow-mapping support), and independent per-island `lambda_P` via the AC-island solve path. The CGMES `normalPF` arrival through `importCGMES` is not covered by a test. |
| `island_diagnostics` | `test/test_island_diagnostics.jl` | AC island diagnostics reporting regressions (source: a 6209-bus CGMES delivery with 158 single-bus islands): the island failure message reports the failing island's own iteration count and stage (never the `iterations=0`/`stage=before_nr` fallback while a per-island record exists), the failing island is looked up in `:ac_island_solver_statuses` when the combined status is untagged, never-attempted islands render as `not_attempted` with `unavailable` mismatch fields instead of inheriting the failed island's statistics, per-island `solver.log`/`mismatch_history.csv` artifacts are written only for attempted islands, `q_limit_processing_status` does not alias `failure_reason`, and the single-island combined-status fallback keeps working; also the no-slack diagnostics (all-isolated vs. isolated-slack vs. never-registered messages) and the `power_flow.auto_slack` promotion (keyword and config path, ENI-over-generator ranking, no-op on a registered slack) |
| `short_circuit` | `test/test_short_circuit.jl` | IEC 60909-0 short circuit (`runShortCircuit!`): hand-derived analytic reference cases with the full derivation in the test comments (feeder reproduces its declared current at the connection point with the c-factor cancelling; series line impedance added on top; synchronous-machine x''d conversion incl. §6.6.3 fictitious resistance; kappa/i_p with the method-b cap; line charging proven irrelevant), IEC Table-1 c-factor selection and the scalar override, the safety-flag contract as a three-way assertion (documented default + log warning + per-row flag), asynchronous-machine skip with island-wide Ik''max lower-bound flag (max case only), `:no_source`/`:isolated` row statuses, `short_circuit.c_factor` config coverage (default, valid override, rejected out-of-band value), and the WebUI data-check gating helper (`_webui_case_has_short_circuit_data`: marker found / not found / cached / unresolvable path stays enabled); the Phase 3 parallel all-bus sweep identity guard (serial vs parallel rows `isequal`-identical for max/min case with `max_tasks=1` and auto on a three-island fixture and on the shipped sp_case60 with `buses = :all`; RAN/SKIPPED stated per run; plus the `sweep_method = :takahashi` selected inverse: machine-precision agreement with `:solves` on values, identical statuses/flags, serial and parallel, deterministic repeats, invalid method rejected) |
| `parallel_foundation` | `test/test_parallel_foundation.jl` | Thread-safety foundation and parallel island path of the multi-core work (Phases 1 and 2): `runtime.parallel.*` configuration defaults and validation (`enabled`, `max_tasks` auto/integer-string, `min_work_items`), the `parallel_max_tasks` resolver, the per-worker performance-profile helpers (`_perf_profile_child` excludes the orchestrator-only `:phase_callback`, `_perf_profile_merge!` sums timings under unchanged phase names, appends iteration rows, prefixes per-worker scalars, and swallows `nothing` children/parents), the solver status as a Net field (registry globals gone, AC/DC fields separate, `deepcopy` carries a working condest thunk), UMFPACK `copy(F)` finalizer safety, the startup summary `parallel:` line, the Web UI `runtime_parallel_enabled` checkbox wiring, and the `solve_parallel` island identity guard: bitwise-equal voltages and iteration counts vs `solve_independent` with `max_tasks=1` and auto on a five-island fixture, identical phase-name sets plus `parallel_wall_time`, per-island prefixed scalars, and the parallel failure semantics (all islands report their real status, no skips). Threaded assertions state RAN/fallback depending on `Threads.nthreads()`; the `--threads=4` battery exercises the parallel branch. |
| `scenarios` | `test/test_scenarios.jl` | Scenario patch model: validation surface (status, set, scale with their target classes and fields), rejections with scenario name and op index (unknown id, class mismatch, inapplicable field, empty ops, duplicate names), the sparlectra.scenarios file round trip through the typed case, the legacy contingencies-to-scenarios mapping, the N-1 expansion equality with generateN1Branches and generateN1Generators on the tracked warmup case including exclusions, apply!/restore! with the undo log (tap patch on an unregulated transformer changes the solver-visible ratio and restores it, a tap patch on a regulated transformer is rejected with the controller named, and the acceptance test: for the tracked warmup case, the shipped sp_case60/sp_case188 and the tracked SCF fixture, every N-1 branch and generator scenario leaves the working copy BITWISE equal to the base after apply! plus restore!, checked with scf_roundtrip_field_diffs), the engine worker resets bitwise between evaluations WITH solves (branch outage, generator outage, patch scenario against the engine template), and `runScenarios!` with mode n1_all writes the identical result CSV as `runContingencies!` on the generated case list (the case118 byte gate against the 7438b6e fixture is the extended group `scenario_engine`); the screening surface: `:off` keeps the historical `ContingencyResult` rows and CSV bytes, `:flag` returns `ScenarioResult` rows with no violating case screened away, screened rows read `start_used = :screen` with an estimate attached, the screening CSV appends exactly the `screened`/`screening_estimate` columns, and `:only` reports estimates without full solves where an estimate exists |
| `demo_cases` | `test/test_demo_cases.jl` | The shipped demo cases (`data/scf/sp_case*`) as regression fixtures: each case is read through `import_case` like any other SCF file and all five run kinds are executed and compared against the tracked fixture in `test/fixtures/demo_cases/` (power flow: iterations, losses, per-bus vm/va to 1e-6 pu keyed by SCF node id; state estimation on the measurements the file itself carries: convergence, dof, objective; IEC 60909 short circuit at the fixture's two buses from the file's own feeder record; N-1 summary counts over `generateN1Branches`; the file's explicit scenarios through `runScenarios!` with per-scenario convergence and max loading), plus the self-built provenance line in the file, the per-case size budget, and the presence of the bad-data measurement variant |
| `contingency` | `test/test_contingency.jl` | N-1 contingency batch API (multi-core Phase 4): case generation (all in-service branches, transformer filter via component type or nonzero ratio, parallel-circuit disambiguation `name#branchIdx` with a shared-name fixture, FOR001 mapping incl. unresolvable names, Phase 2 screening filters `min_vn_kV`/`min_sn_MVA`/`name_pattern` on a 380/110 kV fixture, and per-case weight validation), sp_case14 full N-1 with serial vs parallel field equality (`max_tasks=1` and auto) and base-net immutability, islanding semantics (a load-only island reports the specific `islanded: load-only, X MW` cause, a PV island is promoted and solves), non-convergence and unknown elements reported not thrown, the table printer plus CSV writer, and the Phase 1 start-value ladder (`contingency.rescue_ladder`: `start_used` reported per case, ladder validation of the empty/duplicate/unknown-stage inputs at both the keyword and the config-section level, a longer ladder never losing a warm-converging case and staying serial-vs-parallel identical, the deprecated `retry_flat_start` alias, and a non-converged base case routed through the solver rescue ladder before the flat fallback, verified via a bad-seed fixture whose plain solve diverges and whose no-flat-fallback warning is asserted), and the Phase 2 case weights (`applyContingencyWeights` by name with a default, `readContingencyWeightsCSV` for semicolon/comma/header/comment/negative rows, weight 0 as a pure ranking weight that still solves, and the weight carried into the result and shown in the table/CSV), and the Phase 3 generator outages (`generateN1Generators` count/kind/Pg-name filters, a generator outage solving with the slack absorbing the loss, the slack-unit outage reported as "no slack bus" by default and rescued by `auto_slack`, `distributed_slack_enabled` flowing through, serial-vs-parallel identity, unknown-generator reporting, and a non-pegase stranded-generation fixture where removing a bus's PV unit leaves its island with injection but no reference), and the Phase 4 overload reporting (structured `OverloadRecord` loadings with base loading and delta on a two-parallel-line overload fixture, `severity = weight * max(0, loading - 100)` with the result table severity-ranked failures-first, `shed_load_mw` as a structured field on a load-only island, and `buildContingencyReport`/`printContingencyReport` aggregating outcome counts, worst loading, worst severity, and total load shed), and the Phase 5 Web UI service (`_run_contingency_service` for branch and generator kinds on the shipped sp_case14: succeeded status, `contingency_n1.csv` + `run.log` artifacts, `run_mode`/kind/counts metadata, the slack-unit outage named in the summary badge, and an invalid kind rejected not thrown) |
| `external_grid` | `test/test_external_grid.jl` | External-grid element: the external-grid/distributed-slack mutual-exclusion configuration error (both cover the same imbalance; YAML pair rejected, single option loads), the duck-typing contract guard (`NativeShortCircuitData` field-identical to `CGMESShortCircuitData` by name, order, and type; machine-checked, not comment-checked), PF invariance (an `addExternalGrid!` slack with finite `Sk''` reproduces the manual slack path to machine precision; carried SC data changes no power-flow result, incl. the full ENI record contract), SC hand calculation on a single 110 kV bus (`Zq`, `Ik''`, `Sk''`, kappa/i\_p over two R/X values, c cancelling), min-case semantics (`sk_min` used unflagged with the rx\_min→rx\_max default, explicit `rx_min` wins, missing `sk_min` skips the feeder with the engine's `:no_source` flag), parallel-feeder stacking with unique, removal-safe mrids (suffix continues after the highest surviving id, incl. the same-prefix-bus guard), input validation (`ArgumentError` table), the net-copy regression (`sc_sources` survives `deepcopy`, SC on the copy matches), the `convertSlackToExternalGrid!` demotion (self-consistent bus types right after the demotion, conversion note contents, no-slack/non-slack-bus `ArgumentError`s), the connection statement in the classical result print (`Grid connection:` prose line naming slack vs. source incl. `Sk''`/R-X and the internal slack bus; type column `SOURCE` instead of `SLACK` on the internal bus, also in the structured report rows), and the `internal_impedance` variant (slack moves to the tagged auxiliary bus, stiff `Sk''→∞` limit reproduces the ideal solution below 1e-8 pu, realistic `Sk''` droops voltage and shifts angles, the feeder record stays anchored at the physical connection bus, one internal-impedance grid per bus) |
| `controls` | `test/test_voltage_dependent_control.jl`, `test/test_transformer_phase_shift.jl`, `test/test_tap_controller.jl`, `test/test_series_reactance_control.jl`, `test/test_upfc_control.jl`, `test/test_hvdc_pair_control.jl`, `test/test_tap_changer_model.jl`, `test/test_phase_tap_changer_model.jl`, `test/test_phase_tap_table.jl` | the control-loop console contract (inner-solver diagnostic blocks once per run, one summary line per later pass, no wrong-branch block with the detection off, `control.verbose_passes` for the full per-pass output, the result header naming the passes), Voltage-dependent controls, transformer phase-shift control, tap-controller behavior, YAML controller instantiation (round-trip via the real YAML reader against the programmatic twin, load-time and apply-time validation, idempotent re-apply), machine remote voltage control (`MachineVoltageControl`: API validation, secant convergence onto a remote target, honest `at_limit` on the reactive bounds), the STATCOM current-limit mode (registration validation incl. rating exclusivity and the `i_max_ka` conversion, at-limit `Q = V * S_max` tracking with live bounds across operating points, deadband-level in-range equivalence with the constant-Q mode, inert mode fields on constant-Q controllers, `STATCOM (VSC)` element vocabulary) and the SVC-vs-STATCOM limit contrast on one sagging corridor (acceptance: delivered Q asserted against `V^2 * B` versus `V * S_max` numerically), the FACTS config surface (`s_max_mva` STATCOM and `v_inj_max_pu` SSSC entries instantiate via `applyConfiguredControllers!`, re-apply skips both; regression: the series idempotency check crashed on wrong field names instead of skipping), successful baseline PF preservation when controls are disabled, the `AbstractTapChangerModel`/`convention` field plus `calcRatioTapCorrection`/`calcRatioTapRange` consolidation, the CGMES `PhaseTapChangerModel`/`calcPhaseTapFraction`/`calcPhaseTapAngleRatio`/`calcPhaseTapReactance` formulas plus the DTF importer migration onto them, and the tabular `TapTablePoint`/`kind = :tabular`/`calcPhaseTapTable` override path including the formula-vs-table round-trip regression; the PST X(α) coupling (winding-resolved model lookup, continuous-angle formula evaluation, nearest-row tabular mapping, accepted moves update `x_pu` to the model value at the final angle, devices without reactance data stay static, probe direction consistent with and without tracking incl. angle/reactance restore) and the cross-type warning when a tap controller already regulates a machine controller's target bus; the master/slave transformer group (circulating-Q physics on misaligned parallel units, group stepping with aligned taps and near-equal reactive split under a group-sized deadband, follower exclusivity/self-follow/double-follow/mode validation, the same-target-bus warning, `followers` through the config surface, `(+N followers)` element label; the CGMES side lives in the importer suite: a patched MicroGrid copy with two RatioTapChangers on ONE TapChangerControl imports as one master with one follower); the SVC shunt voltage controller (API validation, in-range secant regulation onto the setpoint with the actuated shunt carrying the same susceptance, honest constant-B `at_limit` on an unreachable target, `clearShuntControllers!`); the MSC/MSR discrete bank mode (registration validation incl. step sign and no-admissible-block range, start-value snapping onto the block grid, whole-block truncated stepping that parks with `status = :parked` before crossing the target without exhausting the outer budget, `at_limit` on the outermost block, `MSC/MSR (switched shunt bank)` element vocabulary with `discrete = true` and step fields on the report row, `step_mvar` through `applyConfiguredControllers!` incl. idempotent re-apply, and continuous-SVC regression) and the generic controllable-element view (`controllableElements` rows for tap, machine, and shunt controllers incl. `ControlRunResult.elements`); the TCSC series-reactance controller (`SeriesReactanceControl`: loop-network target tracking within the deadband, honest `at_limit` clamping on an unreachable target, bit-identical baseline with a disabled controller, element-row vocabulary `:series_x_pu`/`:branch_active_power`, registration validation incl. transformer rejection and the `eps_z` impedance guard, and the no-move deadband case; SSSC injected-voltage mode, limit-form exclusivity validation, converged operation inside the live current-dependent window `x_base ± v_inj_max/|I|`, pinned operation with the effective injected voltage at `v_inj_max`, `SSSC (VSC)` element vocabulary with live window bounds, mode fields on the report row, summary printer naming the SSSC limit, and TCSC-mode regression); the UPFC stationary quadrature composite (`addUpfcControl!`: exact equivalence with the manually registered SSSC+STATCOM pair on a two-corridor fixture, both limit characteristics at their clamps, composite-level rejections and the series-side rollback leaving the net untouched, `UPFC series/shunt (VSC pair, ...)` result-row pairing, YAML type `upfc` incl. the double-apply no-op) and the full DC-link-coupled model (`model = :full`: independent line P and Q reached simultaneously with the DC-link residual closed on a meshed sending-bus fixture, the series active injection `P_se` nonzero, `series_phase = :quadrature` reducing to a standalone SSSC with `P_se = 0`, the `UPFC (full, DC-link coupled)` element row, per-model registration validation, the YAML `model: full` round trip with the double-apply no-op, the low-current `z_add` guard keeping the branch finite at its base impedance on a near-dead line, and the negative-resistance branch left by the full model being rejected loudly by a subsequent IEC 60909 short circuit while the base network runs); the HVDC pair controller (`HvdcPairControl`: registration validation incl. slack and PV guards, exact pairing invariant on a two-island fixture, reversed transfer, rating clamp with honest `at_limit`, per-side voltage-target secant incl. unreachable-target clamping, bit-identical baseline with a disabled controller, MATPOWER `paired_control` mode with the `pf_injections` consistency anchor, YAML `hvdc_pair` round trip, element-row vocabulary `:hvdc_p_transfer_mw`/`:hvdc_transfer`, the grid-forming `mode = :island_feed` incl. mode-dependent validation, the island-draw mirror, the rating clamp, and the YAML `mode: island_feed` declaration; the persistent `HvdcLink` records (`net.hvdcLinks` fills for MATPOWER Stage-0/paired and `addHvdcLink!`, controller attach/detach round trip), the `HVDC Link Flows` result surface (table, header count, converter-loss line, `ACPFlowReport.hvdc_links`), and the one-reference-per-synchronous-island rule (multi-reference error naming the buses, setpoint pair as parallel PQ path next to an AC tie, `invalid_topology` for a demoted `island_feed` reference)) Since 0.20.0 a typed `:asymmetrical` phase model on a winding drives the branch (the model shift at its step, the precedence line in `net.tapModelNotices`) and a phase controller moves the model step (3 to 5) to the flow of a branch built at step 5. |
## Group reference: surface, service and fixture groups

| Extended addition | File | Main checks |
|---|---|---|
| `workshops` | `test/test_workshops.jl` | Every `docs/lit/workshop_*.jl` runs top to bottom in its own module with output and logging silenced; a failing `@assert` in a workshop fails the test. The generated notebooks and pages cannot drift from the library unnoticed. About 50 s. |
| `install` | `test/test_install.jl` | The installation path of a fresh checkout: a copy without the application manifest, `start_webui.jl --env-only` run twice as child processes offline against the depot; the first start reports the environment setup and the compile, the second reports both environments up to date and the packages compiled. |
| `webui` | `test/test_webui.jl` | Focused Web UI coverage for form parsing and backend validation, result rendering, active and terminal timing cards, commit-span omission, tolerance-step hook, path traversal rejection, DTF upload role classification, primary-case and FOR002 selector filtering, the CGMES-export run option (form checkbox + help topic, `export_cgmes` request flag incl. hidden-false/absent-field defaults, and the result-page summary row for completed/failed/absent export metadata), the last-edit-wins precedence between the configuration file and saved case settings (a newer YAML wins for its own keys and sets the notice flag, an older one keeps the sidecar values), the SE phase-5 chain (measurement-set upload classification via content sniff, the SE section of the Runs page (`/stateestimation` is a redirect onto its anchor) and demo generator, the service-level SE run with artifacts and history kind `se`, the SE-started PF and N-1 with the SE run id in the metadata and result parity against a manual run), the measurement generator v2 chain (truth state fresh-solve/from-run with bit-exact state adoption and the rejection paths, deterministic balance-aware single flow ends, passive nodes at a configured sigma or as protected ZI constraints, the se_deltas.csv artifact, and the bad-data threshold surface with staged/legacy service parity and invalid-mode rejection; plus the Delta-u PST chain on the tracked warmup case: generator targets the additional-voltage stepper, the mass release estimates it :pst along the nameplate psi, the fixation lands on the injected step, and the machine trafo stays calculated), the buildSysimage one-call dry run, the sysimage launcher decision from tools/sysimage_launcher.jl (image/metadata present, Julia version, Manifest hash, a src file newer than the image, and the no-terminal default) together with a parity check that `Sparlectra.webui_sysimage_problem` returns the identical verdict for every one of those states (the two implementations cannot share code, because the launcher runs before the package is loaded), the sysimage build-progress contract (a build that stopped reporting for more than a minute counts as gone, so a crashed builder cannot lock the refresh button) and the Sysimage page rendering idle, running, and failed including its build-log tail, a parse check of the generated app CLI driven from the checkout tool, and stubbed route checks without a real asynchronous solver run. Real asynchronous PowerFlow job lifecycles, artifact preview/download/ZIP/history/delete matrices, browser-launcher platform matrices, socket/server lifecycle checks, Markdown help/documentation cross-product validation, repeated real MATPOWER runs, and the scenario editor flows (a three-op scenario created through the form on the shipped sp_case60 bundle (staged copy switched to explicit-list mode), saved, reloaded and run with `:flag` showing the row and the screening columns; the result page's N-1 table showing a screened row with its estimate detail after `:flag`, no screened column after `:off`, and a failed case's error text in the row; a `tap_pos` op on a regulated transformer rejected with the controller named and NOT written to the file; a MATPOWER case getting the export hint with no editor, no `file_block` option on the main form, and an SCF case getting both) are extended-only in `webui_extended`. |
| `contingency_extended` | `test/test_contingency.jl` | N-1 identity slice on the shipped sp_case1354 (first 60 branches): serial vs parallel field equality with `max_tasks=1` and auto, plus the `retry_flat_start` invariance on converged cases |
| `scenario_engine` | `test/test_scenarios.jl` | The engine gate: the N-1 result CSV on the shipped sp_case118 (all branches, all generators, 240 cases through `runContingencies!` on the reused-worker engine) is byte identical to `test/fixtures/contingency_case118_n1_7438b6e.csv`; the screening acceptance on the same batch: with `screening_mode = :flag` and margin 10 no case the full run reports as violating (or failed) is screened away (false-positive count and full-run share are printed, not gated); the same acceptance with distributed slack enabled (augmented lambda system, renormalized participation on a generator outage, screened share and seconds printed); the sp_case300 acceptance per the test-network rule (operated grid, no false negatives at margin 10 with the 0.005 pu trust gate; loud SKIPPED without the local case file); the sp_case1354 run measures runtime SCALING only and runs with `SPARLECTRA_TIMING_TESTS=1` (a SKIPPED line otherwise) (`:off`/`:flag`/`:only` seconds, gated on `SPARLECTRA_LARGE_CASES_DIR`, no screening-quality claim, see the test-network rule); the scenarios workshop (`docs/lit/workshop_scenarios.jl`) is EXECUTED here in a fresh module on the shipped `data/mpower/sp_case118.m`, so the `@assert` next to every number it prints gates the notebook against silent drift on every checkout; the service sources: an SCF case with its own scenarios block runs via `scenario_source = file_block` (screening metadata plus the two CSV columns), an external scenario JSON via `external_file` with `screening_mode = off` keeping the classic columns, `n1_generators` records itself as the case source, unknown sources or a missing scenario file reject as `invalid_request`, and a CGMES delivery with an id-addressed scenario source is rejected with the way out named (export as SCF once, then reference its component ids) |
| `example_infra` | `test/test_example_suites.jl` | Example suite runner infrastructure: CLI parsing, skip and filter logic, CSV and Markdown escaping, registry integrity, and a --help smoke test of the generated runners |
| `matpower_examples` | `test/test_matpower_example.jl` | MATPOWER example runner path, output routing, performance/profile rendering, and runtime configuration forwarding, removed start-voltage alias rejection, and canonical profile-blend parsing |
| `synthetic_grids` | `test/test_synthetic_grids.jl` | Synthetic network generation and larger synthetic-grid regression coverage |
| `configuration_docs` | `test/test_configuration_docs.jl` | Configuration documentation and docs/config consistency checks |
| `cgmes_importer` | `test/test_cgmes_importer.jl` | CGMES importer on synthetic in-memory deliveries and on the checked-in deliveries under `data/cgmes_demo` (exported by `tools/gen_cgmes_fixtures.jl` from sp_case14, sp_case118 and sp_casePST). Synthetic: generic reader semantics (profiles, overlays, references, DifferenceModel skip), the base-voltage inference (`infer_base_voltages` reconstructs stripped catalogs from the SV state, one summary warning, net still solves; without the option the typed abort stays), the self-loop line guard, the import-failure analysis (`importFailureAnalysis` names supplied models, missing vs. satisfied `md:Model.DependentOn` prerequisites, and the verdict; the `import_analysis_mode` service run succeeds on an importable delivery, fails with `import_analysis_not_importable` on a missing declared dependency with the `import_analysis.txt` artifact written, rejects non-CGMES cases, and excludes the other modes), `ReactiveCapabilityCurve` interpolation, result-table alignment with long names, the Web UI documentation wiring, test-set alias/border bookkeeping (table consistency only, no fetch), the synthetic 3W star delivery (3W-leg tap controller from `TapChangerControl` with control-loop convergence, `NonlinearShuntCompensator` mapping, the multi-valued reference notice), remote-regulating machines (default held-PV fallback, `machine_control = true` controller wiring + control-loop convergence, the voltage-held-target fallback), the `short_circuit_mode` rejections (MATPOWER rejected with `short_circuit_requires_cgmes`, mutual exclusion with `diagnose_mode`), `cgmes_import.start_values` having no effect on a MATPOWER run, and the inverted `vset_min_pu`/`vset_max_pu` band as an `ArgumentError`. Checked-in deliveries: `summarizeCGMES` (version 2.4.15, no unresolved reference, no boundary hint), `importCGMES` on the folder and on a zip packed at run time (bus, branch and link counts, slack, no bus without SV, no importer warning), the SV profile and the source case's solution reproduced from the SV start and from a flat start (`compareWithSV` max |dvm| and |dva| below 1e-6, flow rows below 1e-6 MW, bus-by-bus agreement with the source net solved the same way, Q-limits off), the OLTC of sp_case14 arriving as a voltage tap controller only with `tap_control = true` (same target bus, target and deadband as the source case, transformer on the same buses with its tap machinery, control loop converged within the deadband), the PST of sp_casePST keeping its negative shift and its bus link, the sp_case118 synchronous condensers arriving as regulated zero-P units and the asymmetric `minQ`/`maxQ` pairs arriving as the Stage-1 hull, the SCF export of a delivery (mRIDs as `external_id`, deterministic file, re-import to the same solution), `run_sparlectra_api` dispatch on a delivery zip (format and version metadata, `cgmes.log` sections, SV comparison artifacts of the auto start, `start_values` flat and sv winning over a hostile `power_flow.flatstart`, `cgmes.log` written when the solve fails), the `import_analysis_mode` and `short_circuit_mode` service runs (the short-circuit run completes with every row flagged as a lower bound and names the missing machine data, because the exporter writes no x''_d, ratedS or ratedU), and the fixed-reference self-check (SV start verbatim under a `flatstart: true` base config, start residual at noise level, residual CSV, transformer adjacency, `no-SV buses: 0`). Not covered, because the exporter cannot produce it: node-breaker topology, boundary sets, DifferenceModel, CGMES 3.0, HVDC, machine short-circuit attributes, the quirks of real deliveries |
| `cgmes_export` | `test/test_cgmes_export.jl` | CGMES export identity and roundtrip contract (`writeCGMESFiles`, EQ+TP+SSH+SV + optional delivery zip): structural keys in `net.cgmes_ids` with lexicographic bus-pair normalization and first-seen parallel-line numbering, deterministic uuid5 minting (id equals uuid5 over the key), byte-identical re-export and independent-build determinism with a pinned `created` stamp, renaming a component never changes its mRID, the duplicate-mRID guard aborting before any file is written (both offending keys named), the power-flow-identical self-roundtrip over every exported class (transformer with off-nominal ratio + 2° phase shift + ratio-tap machinery, PV machine, slack injection, SVC with rating-derived Q limits, load, shunt, bus link → export → `importCGMES` → `runpf` matches the original within 1e-6 and starts from the exported SV in ≤ 2 iterations, empty notice list), tool provenance in every file (`Generated by Sparlectra.jl v<version>` comment + `md:Model.description`, stamp = pinned `created`), SSH content (machine sign convention, `referencePriority` on the slack unit, `RegulatingControl.targetValue` in kV, shunt sections, `bPerSection`) and SV validation (`compareWithSV` of the re-import: dvm/dva < 1e-12, flow rows < 1e-9), the regulated tap group (master + follower export ONE shared TapChangerControl with controlEnabled on both, reimport regroups them with target and deadband intact), and zip re-import, all on nets built in memory; and the export-import-export identity on the checked-in deliveries (each fixture imported, solved from its SV start with Q-limits off, re-exported with the fixture's header stamp, and compared object by object against the fixture: header verbatim, every mRID present, every attribute equal with numeric values to 1e-9 relative and SvPowerFlow rows to 1e-6 MW). Named as the fields the importer does not carry, and excluded: the RatioTapChanger and TapChangerControl objects (the importer folds the ratio-tap step into the branch ratio; the exporter writes tap machinery from the winding record only), `PowerTransformerEnd.ratedU` of an end that carried a ratio tap changer (the fixture absorbed the step correction into it), `SynchronousMachine.minQ`/`maxQ` (read as the hull), and the names of shunts and their terminals (the importer names shunts by its own bus index) |
| `powsybl_importer` | `test/test_powsybl_importer.jl` | PowSyBl adapter on the five shipped bundles under `test/fixtures/powsybl` (ieee14, ieee57, four_substations, micro_grid_be, eurostag_tie_lines): bundle reader failures (missing required column named with its table, wrong format version), the bundle round trip (every column equal with NaN, element types kept), builder structure (bus count against the bus-breaker rows, links against retained switches, branch count against the built elements, island count against the synchronous components, paired dangling lines skipped into their tie line, bus-view id in the external-id channel), the power flow against the OpenLoadFlow bus voltages within the task bands (ieee14 and ieee57 1e-6 pu and 1e-4 deg, the others 1e-5 pu and 1e-3 deg, run with islands, OLF's balance through the imported participation factors, no Q-limit hysteresis, remote regulation on micro_grid_be) plus a flat-start convergence, the full Y-bus identity of ieee14 against MATPOWER case14 (cache-gated on `case14.m`, RAN or SKIPPED stated), the branch flows of ieee14 against the OLF `p1`/`q1` columns (flow sign, the one-sided shunt excess that the builder carries as a bus shunt added back), patched fixtures on four_substations (one-sided open line, open retained switch as an open link, remote regulation held local with the report entry and as PQ), the slack override with its `override` reason, and the `import_case` dispatch on the bundle and on the `.xiidm` file with the service entry on both. Since 0.20.0 the native IIDM reader (`read_iidm_tables`) is checked on every fixture: its tables equal the bundle's in every schema column except the state columns (rows joined by the index columns), the power flow from the `.xiidm` lands in the same bands as from the bundle, micro_grid_be's stale file state (solved at other tap positions) is dropped with the notice in the report while ieee14's state seeds the start, and a compressed file or a non-IIDM XML is refused with the construct named; the shunt count equals the `shunt_compensators` rows (no branch-derived bus shunt), the N-1 outage of the four_substations PST takes exactly the PST two-port entries out of the Y-bus with the shunt list untouched, the MATPOWER export of micro_grid_be carries `mpc.sparlectra.branch_shunts` and reads back with the Y-bus at 1e-12 and the parts removable per branch, the SE self-test on ieee14 and micro_grid_be passes the global test, and the three bundles under `data/powsybl` are offered and staged by the Web UI selector. |
| `apslf` | `test/test_apslf.jl` | APSLF (AnalyticLoadFlow.jl) integration: the PV-bus regression (ring3 and sp_case5 through `run_sparlectra`, NR against APSLF without polish to 1e-6 pu, radius GRN; guards the series result itself), the AnalyticLoadFlow version guard (`check_apslf_version`), the APSLF workshop (`docs/lit/workshop_apslf.jl`) executed with its assertions on the shipped `data/mpower/sp_case118.m` (no download, no gate), the residual judged against the final active set (a machine clamped at Qmax solves, is logged as a Q-limit event and the run converges), the convergence radius in the status/header and its off switch, `power_flow.solver`/`apslf`/`apslf_start` config validation (including the `solver=apslf` + `apslf_start.enabled` conflict), outer-loop-controller rejection with `solver=apslf`, Web UI form parsing of the new fields into the effective config, the adapter mapping (`PFModel` -> AnalyticLoadFlow spec, PF ordering), the standalone sp_case5 convergence against NR, and the start-value-generator guard including `nr_polish=false`. Nothing here is conditional: AnalyticLoadFlow.jl is a required dependency, so the whole surface runs in every profile that includes the file. |
### `dtf_extended` group

`dtf_extended` is an `adapters` group that owns five DTF test files and invokes
each run function once:

| File | Run function | Main checks |
|---|---|---|
| `test/extended/test_dtf_importer.jl` | `run_dtf_importer_tests` | DTF format parser and direct Net-builder coverage, including voltage-level-index branch conversion, transformer controls, bus-type semantics, parsed outage metadata, and the persisted typed phase-tap model on controlled windings (FOR001E skew/longitudinal regulators carry `phase_taps` with the winding connection angle) |
| `test/extended/test_dtf_for002_validation_example.jl` | `run_dtf_for002_validation_example_tests` | DTF format -> current Julia compatibility module `DTFImporter` -> `Net` -> power-flow validation example smoke coverage against FOR002 diagnostics; verifies generated CSV/Markdown artifacts, lightweight default result/concise CLI output, explicit detailed diagnostics mode, and does not invoke the fast suite |
| `test/extended/test_dtf_for002_outage_validation_example.jl` | `run_dtf_for002_outage_validation_example_tests` | DTF-listed outage validation against FOR002 outage reference reports for the Testnetz13 cases |
| `test/extended/test_dtf_matpower_export_validation_example.jl` | `run_dtf_matpower_export_validation_example_tests` | DTF format -> `Sparlectra.Net` -> existing `writeMatpowerCasefile` -> MATPOWER import roundtrip validation for the Testnetz13 base case and outages listed by the DTF file |
| `test/extended/test_dtf_api_webui_integration.jl` | `run_dtf_api_webui_integration_tests` | DTF format input through the API and Web UI paths |
Optional external FOR001/FOR002 fixtures skip cleanly (as `Broken`, not
failed) when absent. `test/extended/runtests_extended.jl` is an optional
standalone runner that includes the same five files:

```bash
julia --project=. test/extended/runtests_extended.jl
```

## DTF validation examples with FOR002 reference reports

The native DTF examples validate Testnetz13 through the direct DTF path
`DTFImporter.read_dtf` -> `DTFImporter.build_net` -> `runpf!`, without a
MATPOWER intermediate; FOR002 is a legacy textual reference report. External
FOR001/FOR002 datasets are not shipped: place them under `data/DTF/` or pass
`--dtf-file` and `--for002-file`; without them the scripts stop with a
missing-data message.

`examples/run_val_dtf_suite.jl` is the CLI runner for the base, outage and
MATPOWER-export checks plus the DTF import audit. It dispatches to
`examples/dtf/` (`dtf_validation_base.jl`, `dtf_validation_outages.jl`,
`dtf_validation_matpower.jl`, `dtf_validation_audit.jl`), one implementation
per check shared by the CLI, the suite and the extended tests; each module
also runs on its own (for example
`julia --project=. examples/dtf/dtf_validation_base.jl --dtf-file=... --for002-file=...`).

Base-case validation:

```bash
julia --project=. examples/run_val_dtf_suite.jl --mode=base --case=A --data-dir=data/DTF --output-dir=examples/_out/dtf_for002_native_validation --write-csv=true --write-markdown=true
```

Outage validation:

```bash
julia --project=. examples/run_val_dtf_suite.jl --mode=outages --case=A --data-dir=data/DTF --output-dir=examples/_out/dtf_for002_native_outages --write-csv=true --write-markdown=true
```

MATPOWER export validation (`DTFImporter.build_net` -> `writeMatpowerCasefile`
-> `createNetFromMatPowerFile`, checking that the roundtrip does not
materially change the solved result; FOR002 is not its reference):

```bash
julia --project=. examples/run_val_dtf_suite.jl --mode=matpower --case=A --data-dir=data/DTF --output-dir=examples/_out/dtf_matpower_export_testnetz13 --write-csv=true --write-markdown=true --write-matpower=true --run-outages=true
```

Each command writes console output plus CSV and Markdown files in
`--output-dir`: the Markdown files summarize the run, the CSV files keep
row-level bus, generator, branch, KCL, state-residual and metric
diagnostics. The export validation writes `dtf_matpower_export_summary.md`,
`dtf_matpower_export_metrics.csv`, comparison CSV files and the exported
`.m` cases. The exporter uses the optional MATPOWER metadata fields
`mpc.bus_name`, `mpc.branch_name`, `mpc.branch_kind` and
`mpc.for001_contingencies`, so DTF bus names, branch kind (`L`/`T`) and the
parsed outage cards survive `MatpowerIO.read_case_m`. The roundtrip import
keeps DTF PQ generators as fixed injections (no controller
reinterpretation), and the exporter writes TAP `0.0` for line rows and the
explicit ratio for transformer rows.

Metrics:

- `converged`: whether `runpf!` reported convergence.
- `iterations`: Newton iterations of the native solve.
- `final mismatch`: infinity-norm mismatch recomputed from the solved state.
- `max voltage deviation`: largest voltage magnitude or angle difference from the FOR002 bus values.
- `max branch P/Q deviation`: largest directed branch endpoint-flow difference from FOR002.
- `max generator/slack Q deviation`: solved generator/slack reactive output against the FOR002 generator-Q reporting.
- `state residual P/Q`: injection residual from forcing the FOR002 voltages into the native Ybus; CSV only (`dtf_state_residual.csv`, `dtf_validation_metrics.csv` and the outage equivalents), because its rounding-noise floor sits far above real model deviations.

Testnetz13 interpretation: branch-flow deviations are the strongest
validation signal; slack Q is solved by the power flow and must not be
compared with the specified input Q; outage validation executes only the
outages listed in the DTF file, and unmatched FOR002 outage blocks are
reference text, not executed scenarios.

## Pre-merge verification gate

For branches that touch central configuration, MATPOWER runner behavior,
output routing/formatting, performance reporting, or documentation/config
consistency:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl
SPARLECTRA_TEST_PROFILE=extended julia --project=. test/runtests.jl
julia --project=docs docs/make.jl
```

PowerShell:

```powershell
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl

$env:SPARLECTRA_TEST_PROFILE="extended"
julia --project=. test/runtests.jl
Remove-Item Env:SPARLECTRA_TEST_PROFILE

julia --project=docs docs/make.jl
```

## DTF input through the API and the Web UI

The PowerFlow service and the Web UI accept the DTF format as an internal
input path for diagnostics and validation. Use `case_format = :dtf_for001`
for explicit native input, or `case_format = :auto` only when the FOR001
markers are unambiguous; ambiguous `.DAT` files are rejected. The path is
`DTFImporter.read_dtf` -> `DTFImporter.build_net` -> `Sparlectra.Net` ->
`run_sparlectra`/`runpf!`, without a MATPOWER intermediate. In the Web UI
the selector sits in the **Input format** section. FOR002 is a reference
file, not a runnable case: FOR002-like `.DAT` files are hidden from the
primary case selector, can be chosen as the optional FOR002 reference file,
and are never auto-paired with FOR001. DTF-listed outages can be requested
as all outages or by label/index; the default is base-case only.

Artifacts use the PowerFlow run artifact mechanism: `dtf_import_summary.md`,
`dtf_import_summary.csv`, `dtf_for002_base_comparison.md`,
`dtf_for002_base_metrics.csv`, `dtf_native_matpower_export.m`, and
per-outage files such as `dtf_outage_1_summary.md` and
`dtf_outage_1_metrics.csv`. DC lines, HVDC links and active MATPOWER
`mpc.dcline` data are not modeled by this path; they fail with structured
unsupported-DC-line diagnostics instead of being approximated or dropped.

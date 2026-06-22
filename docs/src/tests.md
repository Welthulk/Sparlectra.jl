# Test Suite

Sparlectra uses profile-aware test loading through `test/runtests.jl`.
Profile selection precedence is:
1. CLI argument (`julia --project=. test/runtests.jl <profile>`)
2. `SPARLECTRA_TEST_PROFILE`
3. default `fast`

## Test profiles

| Profile | Command | Scope | Intended use |
|---|---|---|---|
| `fast` (default) | `julia --project=. test/runtests.jl fast` | Core offline tests | Normal local development and default CI smoke |
| `extended` | `julia --project=. test/runtests.jl extended` | Fast + integration/heavier tests | Before merge and after configuration, MATPOWER, or integration changes |
| `all` | `julia --project=. test/runtests.jl all` | Currently alias for `extended` | Reserved for future all-only suites and CI matrix clarity |

`Pkg.test()` uses the same test runner and therefore the default `fast` profile unless `SPARLECTRA_TEST_PROFILE` is set:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Fast profile groups

| Group | Files | Main checks |
|---|---|---|
| `core_model` | `test/testgrid.jl` | Core net construction and validation, inline MATPOWER import helpers including large matrix-block scanning, file-based MATPOWER projected-start normalization, PV/PQ lock-ID resolution, Q-limit enforcement-mode parsing and classical base-failure/no-reenable dispatch checks, rectangular status-row diagnostic preservation, link handling, shunts, reporting/output checks, and summary-file output regression |
| `powerflow_rectangular` | `test/test_solver_interface.jl` | Rectangular power-flow API behavior, sparse-only solver path, Q-limit and typed configuration entry checks, including accepted and rejected framework control-status composition, ordered local MATPOWER batch execution, the thin `run_acpflow` alias, and legacy-keyword rejection |
| `configuration` | `test/test_configuration_coverage.jl`, `test/test_runner_helpers.jl` | Configuration-key coverage, safe refresh checks, forwarding checks, value-domain validation, and test-runner output-mode helper checks |
| `programmatic_api` | `test/test_api.jl` | GUI-ready power-flow API stable unique run IDs, schema versioning, success/failure status, override validation, effective configuration with runtime casefile metadata, classic/full log distinction, timing/status summaries, service phase timing persistence, performance, diagnostic, and detailed CSV artifacts (including buffered/streaming write-mode selection, direct-writer equivalence/timing metadata, non-converged solution-available CSV export, bus-control label cache regression coverage, CSV progress events, comma/semicolon quoting, bounded Q-limit logging/detail artifacts, pre-mutation Q-limit snapshots, non-converged Q-limit validation skip reporting, and contained diagnostic failure), artifact discovery including stale-metadata CSV rescans, Dict/NamedTuple/JSON/YAML serialization, and local service persistent indexing, timestamp fallback, restart recovery, safe single/all-run deletion, run lookup, and artifact/path safety behavior |
| `qlimit_large_case_comparison` | `test/test_qlimit_large_case_comparison.jl` | Offline-safe orchestration tests for the large-case Q-limit mode diagnostic: cached case, successful download resolution, skipped unavailable case, all three mode schedules, and summary rows for skipped, converged, and non-converged outcomes |
| `webui` | `test/test_webui.jl` | Public Web UI binding, standalone app-window command selection, application-root discovery and case/configuration dropdowns, case-settings sidecar profile saves next to runtime MATPOWER cases from freshly completed in-memory runs, type-safe reload/override/rejection behavior for saved profile values, sidecar values flowing through GUI override validation into effective configuration, version rendering, performance/diagnostic/detailed-CSV controls and help (including default comma and opt-in Excel semicolon format), hidden temporary warm-up behavior, allowlisted form mapping including checked and unchecked Q-limit checkbox propagation, complete Markdown-backed form-help coverage, history-preserving secondary-page navigation, allowlisted documentation-link rewriting, documentation-page whitelisting and traversal rejection, enlarged readable artifact/help/docs structure and CSS formatting, service-backed runs, shared logo/header rendering, exact PNG asset routing, result/artifact/history HTML with timestamps, truncated large artifact previews, status badges, phase-aware abort visibility without per-iteration operation-log spam, refresh/delete controls, authoritative output-root handling, explicit, Ctrl-C, and request-aware heartbeat shutdown, listener reuse, absence of direct solver calls, and loopback HTTP smoke requests |
| `state_estimation` | `test/test_state_estimation.jl` | WLS state-estimation behavior and observability-oriented regressions |
| `controls` | `test/test_voltage_dependent_control.jl`, `test/test_transformer_phase_shift.jl`, `test/test_tap_controller.jl` | Voltage-dependent controls, transformer phase-shift control, tap-controller behavior, and successful baseline PF preservation when controls are disabled |

## Extended profile additions

| Extended addition | File | Main checks |
|---|---|---|
| `remove` | `test/testremove.jl` | Remove/delete behavior and consistency after structural edits |
| `pv_voltage_residuals` | `test/test_pv_voltage_residuals.jl` | PV-voltage residual behavior, angle-preserving voltage-setpoint starts, phase-shifted PV integration coverage, and related solver diagnostics |
| `matpower_example` | `test/test_matpower_example.jl` | MATPOWER example runner path, output routing, performance/profile rendering, and runtime configuration forwarding |
| `synthetic_grids` | `test/test_synthetic_grids.jl` | Synthetic network generation and larger synthetic-grid regression coverage |
| `configuration_docs` | `test/test_configuration_docs.jl` | Configuration documentation and docs/config consistency checks |

## Offline and runtime expectations

The default `fast` profile is intended to be offline-safe and should not download MATPOWER cases or run benchmark loops.
Default fast-profile output is intentionally compact: the runner prints the selected profile, one `[n/8]` marker per group, and Julia's final test summary. MATPOWER import diagnostics, auto-profile tables, runtime casefile banners, Q-limit tables, and similar artifact-oriented diagnostic blocks are suppressed in normal test stdout so progress remains scannable.

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

The `extended` profile may include MATPOWER/example/output-heavy tests. These tests stay isolated from the default profile.

Use `fast` during normal development. Use `extended` before merging changes that affect configuration, MATPOWER import, output formatting, performance reporting, or broader integration paths.

## Fast-profile volume review

The fast profile currently contains a mix of true unit/smoke coverage and several integration-style service/UI paths:

- True unit or focused smoke tests: configuration key/value validation, MATPOWER auto-profile decision rules on tiny synthetic cases, rectangular solver API checks with small fixtures, core model invariants, state-estimation smoke/regression cases, and control-loop unit/regression checks.
- Integration-style tests that remain in fast because they protect recent public behavior: API service request/metadata/artifact smoke checks, Web UI form rendering and routing, allowlisted documentation/help routing, operation-log safety, run deletion safety, and small service-backed Web UI/API runs.
- Heavier or broader tests already isolated in extended: MATPOWER example runner coverage, synthetic-grid regressions, configuration documentation consistency, PV residual integration coverage, and structural remove/delete behavior.
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
[1/5] core_model
[2/5] powerflow_rectangular
...
```

Julia's final `Test Summary` remains unchanged and visible at the end.

## Rectangular/Q-limit diagnostics in tests

The rectangular convergence and Q-limit active-set diagnostic block is not printed in normal test runs.
Those diagnostics remain available only through explicit diagnostic requests (for example solver `verbose > 0` paths used during focused debugging).

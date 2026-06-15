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
| `core_model` | `test/testgrid.jl` | Core net construction and validation, inline MATPOWER import helpers, file-based MATPOWER projected-start normalization, PV/PQ lock-ID resolution, and rectangular status-row diagnostic preservation, link handling, shunts, reporting/output checks, and summary-file output regression |
| `powerflow_rectangular` | `test/test_solver_interface.jl` | Rectangular power-flow API behavior, sparse-only solver path, Q-limit and typed configuration entry checks, including accepted and rejected framework control-status composition, ordered local MATPOWER batch execution, the thin `run_acpflow` alias, and legacy-keyword rejection |
| `configuration` | `test/test_configuration_coverage.jl` | Configuration-key coverage, forwarding checks, and value-domain validation |
| `programmatic_api` | `test/test_api.jl` | GUI-ready power-flow API stable unique run IDs, schema versioning, success/failure status, override validation, effective configuration, classic/full log distinction, timing/status summaries, service phase timing persistence, performance, diagnostic, and detailed CSV artifacts (including comma/semicolon quoting and contained diagnostic failure), artifact discovery, Dict/NamedTuple/JSON/YAML serialization, developer large-case Web UI cache benchmark helper smoke coverage, and local service persistent indexing, timestamp fallback, restart recovery, safe single/all-run deletion, run lookup, and artifact/path safety behavior |
| `webui` | `test/test_webui.jl` | Public Web UI binding, standalone app-window command selection, application-root discovery and case/configuration dropdowns, version rendering, performance/diagnostic/detailed-CSV controls and help (including default comma and opt-in Excel semicolon format), hidden temporary warm-up behavior, allowlisted form mapping, complete Markdown-backed form-help coverage, history-preserving secondary-page navigation, allowlisted documentation-link rewriting, documentation-page whitelisting and traversal rejection, enlarged readable artifact/help/docs structure and CSS formatting, service-backed runs, shared logo/header rendering, exact PNG asset routing, result/artifact/history HTML with timestamps, status badges, phase-aware abort visibility without per-iteration operation-log spam, refresh/delete controls, authoritative output-root handling, explicit, Ctrl-C, and request-aware heartbeat shutdown, listener reuse, absence of direct solver calls, and loopback HTTP smoke requests |
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

The `extended` profile may include MATPOWER/example/output-heavy tests. These tests stay isolated from the default profile.

Use `fast` during normal development. Use `extended` before merging changes that affect configuration, MATPOWER import, output formatting, performance reporting, or broader integration paths.

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

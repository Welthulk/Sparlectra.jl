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
| `powerflow_rectangular` | `test/test_solver_interface.jl` | Rectangular power-flow API behavior, sparse-only solver path, Q-limit and typed configuration entry checks, current-iteration start rejection diagnostics on a tiny synthetic case, including accepted and rejected framework control-status composition, ordered local MATPOWER batch execution, the thin `run_acpflow` alias, and legacy-keyword rejection |
| `configuration` | `test/test_configuration_coverage.jl`, `test/test_runner_helpers.jl` | Configuration-key coverage, safe refresh checks including current-iteration start defaults, forwarding checks, start-voltage value-domain validation, and test-runner output-mode helper checks |
| `matpower_metadata` | `test/test_matpower_metadata.jl` | MATPOWER parser metadata fields, legacy bus sorting of `bus_name`, opt-in bus-name import, branch-kind override, branch metadata retention, FOR001 contingency mapping, and default/opt-in `mpc.dcline` behavior |
| `programmatic_api` | `test/test_api.jl` | GUI-ready power-flow API stable unique run IDs, schema versioning, success/failure status including unsupported active MATPOWER `mpc.dcline` abort metadata/log/Web UI rendering, override validation including current-iteration start options, effective configuration with runtime casefile metadata, classic/full log distinction, timing/status summaries, service phase timing persistence, performance, diagnostic, and detailed CSV artifacts (including buffered/streaming write-mode selection, direct-writer equivalence/timing metadata, non-converged solution-available CSV export, bus-control label cache regression coverage, CSV progress events, comma/semicolon quoting, bounded Q-limit logging/detail artifacts, pre-mutation Q-limit snapshots, non-converged Q-limit validation skip reporting, and contained diagnostic failure), artifact discovery including stale-metadata CSV rescans, Dict/NamedTuple/JSON/YAML serialization, and local service persistent indexing, timestamp fallback, restart recovery, safe single/all-run deletion, run lookup, and artifact/path safety behavior |
| `webui` | `test/test_webui.jl` | Public Web UI binding, standalone app-window command selection including Linux app-window and generic opener fallbacks, application-root discovery and case/configuration dropdowns, case-settings sidecar profile saves next to runtime MATPOWER cases from freshly completed in-memory runs, type-safe reload/override/rejection behavior for saved profile values, sidecar values flowing through GUI override validation into effective configuration, version rendering, performance/diagnostic/detailed-CSV controls and help (including default comma and opt-in Excel semicolon format), hidden temporary warm-up behavior, allowlisted form mapping including start-voltage dropdown values, collapsible current-iteration start controls with hidden-false checkbox propagation, and checked/unchecked Q-limit checkbox propagation, complete Markdown-backed form-help coverage, history-preserving secondary-page navigation, allowlisted documentation-link rewriting, documentation-page whitelisting and traversal rejection, enlarged readable artifact/help/docs structure and CSS formatting, service-backed runs, shared logo/header rendering, exact PNG asset routing, result/artifact/history HTML with timestamps, truncated large artifact previews, status badges, phase-aware abort visibility without per-iteration operation-log spam, refresh/delete controls, authoritative output-root handling, explicit, Ctrl-C, and request-aware heartbeat shutdown, listener reuse, absence of direct solver calls, and loopback HTTP smoke requests |
| `state_estimation` | `test/test_state_estimation.jl` | WLS state-estimation behavior and observability-oriented regressions |
| `controls` | `test/test_voltage_dependent_control.jl`, `test/test_transformer_phase_shift.jl`, `test/test_tap_controller.jl` | Voltage-dependent controls, transformer phase-shift control, tap-controller behavior, and successful baseline PF preservation when controls are disabled |

## Extended profile additions

| Extended addition | File | Main checks |
|---|---|---|
| `remove` | `test/testremove.jl` | Remove/delete behavior and consistency after structural edits |
| `pv_voltage_residuals` | `test/test_pv_voltage_residuals.jl` | PV-voltage residual behavior, angle-preserving voltage-setpoint starts, phase-shifted PV integration coverage, and related solver diagnostics |
| `matpower_example` | `test/test_matpower_example.jl` | MATPOWER example runner path, output routing, performance/profile rendering, and runtime configuration forwarding, removed start-voltage alias rejection, and canonical profile-blend parsing |
| `synthetic_grids` | `test/test_synthetic_grids.jl` | Synthetic network generation and larger synthetic-grid regression coverage |
| `configuration_docs` | `test/test_configuration_docs.jl` | Configuration documentation and docs/config consistency checks |
| `dtf_importer` | `test/extended/test_dtf_importer.jl` | Native DTF/FOR001 parser and direct Net-builder coverage, including voltage-level-index branch conversion, transformer controls, bus-type semantics, and parsed outage metadata |
| `dtf_for002_validation_example` | `test/extended/test_dtf_for002_validation_example.jl` | Native FOR001/DTF -> `DTFImporter` -> `Net` -> power-flow validation example smoke coverage against FOR002 diagnostics; verifies generated CSV/Markdown artifacts, lightweight default result/concise CLI output, explicit detailed diagnostics mode, and does not invoke the fast suite |
| `dtf_matpower_export_validation_example` | `test/extended/test_dtf_matpower_export_validation_example.jl` | DTF -> `Sparlectra.Net` -> existing `writeMatpowerCasefile` -> MATPOWER import roundtrip validation for the Testnetz13 base case and DTF-listed outages |

## Native DTF/FOR002 validation examples

The native DTF/FOR002 examples validate Testnetz13 through the direct DTF path:
`DTFImporter.read_dtf` -> `DTFImporter.build_net` -> `runpf!`. They deliberately avoid MATPOWER import/export and the generated FOR001 builder so that native DTF parsing, Net construction, and solved branch-flow reporting are exercised directly. FOR002 is used as a legacy textual reference report.

Run the base-case validation with:

```bash
julia --project=. examples/validate_dtf_for002_testnetz13.jl --dtf-file=test/fixtures/dtf/FOR001.DAT --for002-file=examples/FOR002.DAT --output-dir=examples/_out/dtf_for002_native_validation --write-csv=true --write-markdown=true
```

Run the outage validation with:

```bash
julia --project=. examples/validate_dtf_for002_outages_testnetz13.jl --dtf-file=test/fixtures/dtf/FOR001.DAT --for002-file=examples/FOR002.DAT --output-dir=examples/_out/dtf_for002_native_outages --write-csv=true --write-markdown=true
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
julia --project=. examples/validate_dtf_matpower_export_testnetz13.jl --dtf-file=test/fixtures/dtf/FOR001.DAT --for002-file=examples/FOR002.DAT --output-dir=examples/_out/dtf_matpower_export_testnetz13 --write-csv=true --write-markdown=true --write-matpower=true --run-outages=true
```

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
- `state residual P/Q`: injection residual from forcing FOR002 printed voltage magnitudes/angles into the native Ybus and comparing the calculated injections with the FOR002 bus table.

Current Testnetz13 interpretation:

- Branch-flow deviations are small and are the strongest validation signal for the native DTF path.
- Slack Q is solved by the power flow and should not be compared with the specified input Q as if it were fixed.
- State residuals are diagnostic, not hard pass/fail criteria; they are sensitive to FOR002 rounding and transformer-adjacent bus voltage/angle differences.
- Outage validation currently executes only the outages listed in FOR001/DTF.
- FOR002 may contain more outage blocks than FOR001 lists; unmatched FOR002 blocks are treated as reference text, not executed scenarios.

## Offline and runtime expectations

The default `fast` profile is intended to be offline-safe and should not download MATPOWER cases or run benchmark loops.
The Q-limit large-case comparison workflow is a manual diagnostic tool: it resolves optional large MATPOWER cases through Sparlectra's case registry/download support, compares `case × start_profile × qlimit_mode`, and uses its CSV/JSON summaries as the primary output. Real runs for cases such as `case13659pegase.m` and `case_SyntheticUSA.m` are intentionally excluded from the fast tests; automated coverage uses stubs for case resolution and API execution.
The experimental large-case Q-limit comparison test block is suppressed from the normal fast profile; run `test/test_qlimit_large_case_comparison.jl` manually when maintaining that tool.

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

The `extended` profile may include MATPOWER/example/output-heavy tests and native DTF/FOR002 diagnostic-example checks. These tests stay isolated from the default profile.

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

## Experimental/internal DTF/FOR001 Web/API input path

The PowerFlow service and Web UI include an experimental/internal DTF/FOR001 input path for diagnostics and validation. This path is deliberately cautious and is not yet announced as a public supported file format. Use `case_format = :dtf_for001` for explicit native input, or `case_format = :auto` only when the FOR001 markers are unambiguous; ambiguous `.DAT` files are rejected instead of being silently interpreted.

The native API path uses `DTFImporter.read_dtf` -> `DTFImporter.build_net` -> `Sparlectra.Net` -> `run_sparlectra`/`runpf!` and does not solve through a MATPOWER intermediate conversion. The Web UI places the selector in the advanced/internal **Input format** section with the cautious label "DTF/FOR001 diagnostics (experimental/internal)". FOR002 can be supplied as an optional legacy reference report for diagnostics. DTF-listed outages can be requested explicitly as all outages or selected outage labels/indices; the default remains base-case-only.

Generated artifacts use the existing PowerFlow run artifact mechanism. Stable DTF artifact names include `dtf_import_summary.md`, `dtf_import_summary.csv`, `dtf_for002_base_comparison.md`, `dtf_for002_base_metrics.csv`, `dtf_native_matpower_export.m`, and per-outage files such as `dtf_outage_1_summary.md` and `dtf_outage_1_metrics.csv`.

Limitation/TODO: DC lines, HVDC links, and active MATPOWER `mpc.dcline` data are not modeled by this native DTF/MATPOWER power-flow path. They fail clearly with structured unsupported-DC-line diagnostics instead of being approximated or dropped.

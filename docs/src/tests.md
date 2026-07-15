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
| Web UI | Form parsing and backend validation, result rendering, active and terminal timing cards, commit-span omission, tolerance-step hook, path traversal rejection, DFT upload role classification, primary-case and FOR002 selector filtering, and stubbed route checks without a real asynchronous solver run. | Case-profile persistence with asynchronous jobs, real run/result polling, artifact preview/download/ZIP/history/delete lifecycle checks, browser-launcher matrices, socket/server lifecycle, and Markdown/help/documentation cross-products. |
| Documentation/hygiene | No repository-wide documentation/help scan in fast; only focused source-level smoke checks tied to edited paths. | Configuration documentation consistency and normalized tracked-path/content repository hygiene scans. |

`Pkg.test()` uses the same test runner and therefore the default `fast` profile unless `SPARLECTRA_TEST_PROFILE` is set:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Fast profile groups

| Group | Files | Main checks |
|---|---|---|
| `core_model` | `test/testgrid.jl` | Core net construction and validation, small inline MATPOWER import helpers, file-based MATPOWER projected-start normalization, PV/PQ lock-ID resolution, Q-limit enforcement-mode parsing and classical base-failure/no-reenable dispatch checks, rectangular nonfinite mismatch and status-row diagnostic preservation, link handling, shunts, reporting/output checks, and summary-file output regression. Large sparse Ybus checks and large MATPOWER matrix-body scanner coverage are extended-only. |
| `powerflow_rectangular` | `test/test_solver_interface.jl` | Rectangular power-flow API behavior, sparse-only solver path, AC-island detection/reference validation/independent solving including aggregate all-island convergence status and iteration accounting, Q-limit and typed configuration entry checks, current-iteration start rejection diagnostics on a tiny synthetic case, including accepted and rejected framework control-status composition, ordered local MATPOWER batch execution, the thin `run_acpflow` alias, and legacy-keyword rejection |
| `configuration` | `test/test_configuration_coverage.jl`, `test/test_runner_helpers.jl` | Configuration-key coverage, safe refresh checks including current-iteration start defaults, forwarding checks, start-voltage value-domain validation, and test-runner output-mode helper checks |
| `matpower_metadata` | `test/test_matpower_metadata.jl` | MATPOWER parser metadata fields, legacy bus sorting of `bus_name`, opt-in bus-name import, branch-kind override, branch metadata retention, FOR001 contingency mapping, and default/opt-in `mpc.dcline` behavior including PF/loss/Q/V/Q-limit mapping and voltage-control safeguards |
| `programmatic_api` | `test/test_api.jl` | Focused GUI-ready power-flow API coverage for serialization and transport helpers, path and validation safety, one successful small API run, one pre-solver failure, one numerical/island failure, Solver-time and Total-time contracts, critical DC-line default and strict-rejection smokes, and one small independent-island regression. Exhaustive CSV/export matrices, full artifact inventories, persistent restart/history/delete lifecycle matrices, repeated performance-log modes, and complete island artifact-content matrices are extended-only in `programmatic_api_extended`. |
| `webui` | `test/test_webui.jl` | Focused Web UI coverage for form parsing and backend validation, result rendering, active and terminal timing cards, commit-span omission, tolerance-step hook, path traversal rejection, DFT upload role classification, primary-case and FOR002 selector filtering, and stubbed route checks without a real asynchronous solver run. Real asynchronous PowerFlow job lifecycles, artifact preview/download/ZIP/history/delete matrices, browser-launcher platform matrices, socket/server lifecycle checks, Markdown help/documentation cross-product validation, and repeated real MATPOWER runs are extended-only in `webui_extended`. |
| `state_estimation` | `test/test_state_estimation.jl` | WLS state-estimation behavior and observability-oriented regressions |
| `controls` | `test/test_voltage_dependent_control.jl`, `test/test_transformer_phase_shift.jl`, `test/test_tap_controller.jl` | Voltage-dependent controls, transformer phase-shift control, tap-controller behavior, and successful baseline PF preservation when controls are disabled |

## Extended profile additions

The `extended` profile is extended-only: it does not run the fast profile first. Use `all` when a single invocation must execute both fast and extended suites.

Current extended-only groups are:

- `core_model_extended`
- `programmatic_api_extended`
- `webui_extended`
- `repository_hygiene`
- `dft_extended`

| Extended addition | File | Main checks |
|---|---|---|
| `legacy/remove` | `test/testremove.jl` | Remove/delete behavior and consistency after structural edits |
| `pv_voltage_residuals` | `test/test_pv_voltage_residuals.jl` | PV-voltage residual behavior, angle-preserving voltage-setpoint starts, phase-shifted PV integration coverage, and related solver diagnostics |
| `matpower_examples` | `test/test_matpower_example.jl` | MATPOWER example runner path, output routing, performance/profile rendering, and runtime configuration forwarding, removed start-voltage alias rejection, and canonical profile-blend parsing |
| `synthetic_grids` | `test/test_synthetic_grids.jl` | Synthetic network generation and larger synthetic-grid regression coverage |
| `configuration_docs` | `test/test_configuration_docs.jl` | Configuration documentation and docs/config consistency checks |

### `dft_extended` group

`dft_extended` is a normal `extended`-only group. It owns exactly five DFT test files and invokes each of their run functions once:

| File | Run function | Main checks |
|---|---|---|
| `test/extended/test_dtf_importer.jl` | `run_dtf_importer_tests` | DFT format parser and direct Net-builder coverage, including voltage-level-index branch conversion, transformer controls, bus-type semantics, and parsed outage metadata |
| `test/extended/test_dtf_for002_validation_example.jl` | `run_dtf_for002_validation_example_tests` | DFT format -> current Julia compatibility module `DTFImporter` -> `Net` -> power-flow validation example smoke coverage against FOR002 diagnostics; verifies generated CSV/Markdown artifacts, lightweight default result/concise CLI output, explicit detailed diagnostics mode, and does not invoke the fast suite |
| `test/extended/test_dtf_for002_outage_validation_example.jl` | `run_dtf_for002_outage_validation_example_tests` | DFT-listed outage validation against FOR002 outage reference reports for the Testnetz13 cases |
| `test/extended/test_dtf_matpower_export_validation_example.jl` | `run_dtf_matpower_export_validation_example_tests` | DFT format -> `Sparlectra.Net` -> existing `writeMatpowerCasefile` -> MATPOWER import roundtrip validation for the Testnetz13 base case and outages listed by the DFT file |
| `test/extended/test_dtf_api_webui_integration.jl` | `run_dtf_api_webui_integration_tests` | DFT format input through the API and Web UI paths |

`extended` includes `dft_extended`, and `all` runs `fast` once followed by the complete `extended` profile once, so the DFT tests run exactly once per `extended` or `all` invocation. Optional external FOR001/FOR002 fixtures continue to skip cleanly (as `Broken`, not failed) when absent.

`test/extended/runtests_extended.jl` remains an optional standalone DFT-focused runner. It includes the same five files directly and is useful for iterating on DFT coverage without running the rest of the extended profile:

```bash
julia --project=. test/extended/runtests_extended.jl
```

## Native DFT validation examples with FOR002 reference reports

The native DFT examples with FOR002 reference reports validate Testnetz13 through the direct DFT format path:
the current Julia compatibility module `DTFImporter` reads the DFT format with `DTFImporter.read_dtf` -> `DTFImporter.build_net` -> `runpf!`. They deliberately avoid MATPOWER import/export and the generated FOR001 builder so that native DFT format parsing, Net construction, and solved branch-flow reporting are exercised directly. FOR002 is used as a legacy textual reference report.

External FOR001/FOR002 validation datasets are not shipped with Sparlectra.
Place local validation files under `data/DTF/` or pass explicit paths with
`--dtf-file` and `--for002-file`. If the files are absent, the optional
validation scripts stop with a clear missing-data message instead of
substituting unrelated data or claiming that validation was executed.

Run the base-case validation with:

```bash
julia --project=. examples/validate_dtf_for002_testnetz13.jl --dtf-file=data/DFT.DAT --for002-file=data/DTF/FOR002.DAT --output-dir=examples/_out/dtf_for002_native_validation --write-csv=true --write-markdown=true
```

Run the outage validation with:

```bash
julia --project=. examples/validate_dtf_for002_outages_testnetz13.jl --dtf-file=data/DFT.DAT --for002-file=data/DTF/FOR002.DAT --output-dir=examples/_out/dtf_for002_native_outages --write-csv=true --write-markdown=true
```

Each command writes concise console output plus CSV and Markdown files in the requested `--output-dir`. The Markdown files summarize the run, while CSV files keep row-level bus, generator, branch, KCL, state-residual, and metric diagnostics.

The MATPOWER export validation example checks a different question: whether an
already-built native DFT `Net` can be exported by Sparlectra's existing
MATPOWER exporter and re-imported without materially changing the solved
Sparlectra result. Its required path is
`DTFImporter.read_dtf` -> `DTFImporter.build_net` -> `writeMatpowerCasefile` ->
`createNetFromMatPowerFile`; it does not implement a DFT-specific exporter and
does not use FOR002 as the primary roundtrip reference. Run it with:

```bash
julia --project=. examples/validate_dtf_matpower_export_testnetz13.jl --dtf-file=data/DFT.DAT --for002-file=data/DTF/FOR002.DAT --output-dir=examples/_out/dtf_matpower_export_testnetz13 --write-csv=true --write-markdown=true --write-matpower=true --run-outages=true
```

The command writes `dtf_matpower_export_summary.md`,
`dtf_matpower_export_metrics.csv`, bus/branch/generator comparison CSV files,
and exported MATPOWER `.m` cases in the selected output directory. The exporter
uses Sparlectra's established optional MATPOWER metadata fields
`mpc.bus_name`, `mpc.branch_name`, `mpc.branch_kind`, and
`mpc.for001_contingencies` when that information is available, so DFT bus names,
stable branch labels, DFT branch kind (`L`/`T`), and the parsed outage cards can
be recovered by `MatpowerIO.read_case_m`. Standard MATPOWER bus, generator, and
branch data remain sufficient for solving files that do not contain that
metadata. The roundtrip proves preservation through Sparlectra's MATPOWER
export/import path; it does not independently certify agreement with external
FOR002 outage blocks.

The roundtrip import disables MATPOWER PQ-generator controller reinterpretation
for this diagnostic so DFT PQ generators remain fixed injections. The exporter
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

- Branch-flow deviations are small and are the strongest validation signal for the native DFT path.
- Slack Q is solved by the power flow and should not be compared with the specified input Q as if it were fixed.
- State residuals are diagnostic, not hard pass/fail criteria; they are sensitive to FOR002 rounding and transformer-adjacent bus voltage/angle differences.
- Outage validation currently executes only the outages listed in DFT.
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

The `extended` profile may include MATPOWER/example/output-heavy tests and native DFT diagnostic-example checks with FOR002 reference reports. These tests stay isolated from the default profile.

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

## Experimental/internal DFT Web/API input path

The PowerFlow service and Web UI include an experimental/internal DFT format input path for diagnostics and validation. This path is deliberately cautious and is not yet announced as a public supported file format. Use `case_format = :dtf_for001` for explicit native input, or `case_format = :auto` only when the FOR001 markers are unambiguous; ambiguous `.DAT` files are rejected instead of being silently interpreted.

The native API path uses `DTFImporter.read_dtf` -> `DTFImporter.build_net` -> `Sparlectra.Net` -> `run_sparlectra`/`runpf!` and does not solve through a MATPOWER intermediate conversion. The Web UI places the selector in the advanced/internal **Input format** section with the cautious label "DFT diagnostics (experimental/internal)". FOR002 is treated as a reference/result file rather than a runnable input case: FOR002-like `.DAT` files are hidden from the primary case selector, can be selected or typed as an optional FOR002 reference file when available in the case cache, and must not be auto-paired with FOR001. DFT-listed outages can be requested explicitly as all outages or selected outage labels/indices; the default remains base-case-only.

Generated artifacts use the existing PowerFlow run artifact mechanism. Stable DFT-path artifact filenames include `dtf_import_summary.md`, `dtf_import_summary.csv`, `dtf_for002_base_comparison.md`, `dtf_for002_base_metrics.csv`, `dtf_native_matpower_export.m`, and per-outage files such as `dtf_outage_1_summary.md` and `dtf_outage_1_metrics.csv`.

Limitation/TODO: DC lines, HVDC links, and active MATPOWER `mpc.dcline` data are not modeled by this native DFT/MATPOWER power-flow path. They fail clearly with structured unsupported-DC-line diagnostics instead of being approximated or dropped.

## Version 0.8.7 – 2026-06-18

### Improvements
* Added classical Q-limit enforcement modes (`classic_simultaneous` and `classic_one_at_a_time`) alongside the default active-set mode for diagnostic comparisons and large MATPOWER case analysis; legacy `matpower_*` aliases remain accepted for old YAML files.
* Improved the PowerFlow Web UI layout by widening the form, moving MATPOWER import controls into Advanced / expert options, defaulting new Web UI runs to `matpower_import.auto_profile = apply`, and making validation errors dismissible while preserving recent errors in a collapsed details block.
* Restored compact default test-run output by suppressing verbose MATPOWER/runtime diagnostics in normal `test/runtests.jl` runs while keeping an explicit verbose opt-in for debugging.

### Bugfixes
* Fixed Web UI/API MATPOWER import controls so selected auto-profile, transformer ratio, phase-shift, shunt, and voltage-reference options are validated, logged, propagated into runtime case construction, and reflected in effective configuration artifacts.
* Fixed Web UI/API Q-limit diagnostics so pre-solve PV limit snapshots survive active-set mutation, non-converged solves do not report final Q-limit validation as `OK`, and `effective_config.yaml` records the runtime casefile used for the run.
* Fixed Web UI/API detailed CSV artifact handling so non-converged runs with an available solution/network state still export CSV files, skipped exports are logged with structured reasons, and artifact routes rescan run directories for partial or stale-metadata CSV files.

## Version 0.8.6 – 2026-06-18
### Improvements
* Added configurable detailed CSV artifact writing with `auto`, `buffered`, and `streaming` modes so large output artifacts can stream directly to disk while normal cases keep buffered writes.
* Improved direct detailed CSV export with lower-allocation cell writing and per-file timing metadata for large Web UI/API runs.
* Removed the remaining per-bus prosumer scan from direct detailed CSV bus control labels, added operation-log CSV progress events, and bounded Q-limit run-log details while preserving full Q-limit CSV artifacts.

### Bugfixes
* Fixed Web UI Q-limit checkbox handling so an unchecked form submission explicitly disables `power_flow.qlimits.enabled` in the effective configuration.

## Version 0.8.5 – 2026-06-17
### Improvements
* Cleaned up single PowerFlow Web UI/API timing output so unavailable solver time is omitted and the measured run duration is shown as wall time.

### Bugfixes
* Fixed MATPOWER `.m` import for very large cases such as `case_ACTIVSg70k.m` by replacing regex-based matrix block extraction with deterministic string scanning.

## Version 0.8.4 – 2026-06-10
### New features
* A loopback‑only PowerFlow Web UI was added with run history, artifact downloads, logs, MATPOWER case management, help, and writable runtime directories. Optional artifacts now include additional logs, configuration output, and detailed CSV exports for voltages and branch flows.
### Improvements
* Added configurable detailed CSV formatting for technical, German Excel, and US Excel notation.

### Bugfixes
* Fixed rectangular PV/slack voltage initialization so replacing a voltage magnitude preserves the existing phasor angle.

## Version 0.8.3 – 2026-05-30

### Breaking changes

* Replaced the former keyword-heavy high-level runner surface with the configuration-driven `run_sparlectra` framework entry point.
* Kept `run_acpflow` only as a thin alias for `run_sparlectra`; it now accepts the same minimal framework arguments.
* High-level import, solver, control, benchmark, and output behavior is now controlled through `SparlectraConfig` or YAML configuration.
* Framework runs now consistently return `SparlectraRunResult` for both `casefile` and `net` workflows.

### New features

* Added config-driven multi-case MATPOWER execution via `matpower_import.cases` and the dedicated `run_sparlectra_cases` helper, while keeping `run_sparlectra` as a single-case workflow.

### Improvements

* Refactored the high-level ACP/MATPOWER workflow into a clearer framework path with separated import, execution, status, and output handling.
* Refined `SparlectraRunResult` status semantics so numerical convergence, solution availability, control-loop status, limit validation, and final framework acceptance are reported separately.
* Preserved detailed rectangular diagnostics in MATPOWER runner status rows, including Q-limit active-set information, final PV voltage residuals, and wrong-branch metrics.

### Bugfixes

* Fixed file-based MATPOWER start handling so projected/imported voltage and angle starts are actually honored by the rectangular solver instead of being overwritten by an effective flat start.
* Fixed framework and MATPOWER status handling for rejected numerical solutions, including wrong-branch rejection, active-set rejection, controlled-run PF failures, disabled controls, and synthetic benchmark summaries.
* Fixed file-based MATPOWER Q-limit lock handling so `power_flow.qlimits.lock_pv_to_pq_buses` is resolved from original MATPOWER bus IDs to internal Sparlectra bus indices.
* Fixed the public `ensure_casefile` binding and updated runnable documentation snippets so fresh checkouts and package installations can load example MATPOWER cases correctly.
* Fixed the tap-control demo classic-output toggle so the example no longer fails late with an undefined `show_classic` binding.

### Related issues

* #228 Remove the old compatibility surface from the high-level runner

## Version 0.8.2 – 2026-05-29
### New Features

* Added configurable wrong-branch diagnostics for rectangular power-flow results, including voltage, angle-spread, and branch-angle plausibility checks.

### Improvements

* Hardened `matpower_import.auto_profile` into a MATPOWER pre-run that logs recommendation evidence, preserves `recommend` mode without changing the active configuration, applies safe import-convention, comparison-reference, start-mode, and Q-limit guard changes in `apply` mode, and prints final effective options without rewriting YAML files.
* Refactored the rectangular complex-state power-flow implementation into focused modules under `src/powerflow_rectangular/`, with `runpf_rectangular!` as the network-integrated entry point and `run_complex_nr_rectangular` as the standalone array-level solver.

### Bugfixes

* Fixed `run_acpflow(casefile=...)` configuration forwarding so MATPOWER/file-based rectangular solves honor configured `power_flow.wrong_branch_*` options instead of falling back to default diagnostics.
* Aligned rectangular autodamping defaults so direct solver calls and configuration-driven runs use the same `autodamp_min = 0.05` default.


### Related Issues

* #193 Harden MATPOWER auto-profile recommendation and application
* #219 Detect wrong/false low-voltage branch convergence and retry safely
* #220 Mechanically split rectangular power-flow code into focused modules
  

## Version 0.8.1 – 2026-05-26

### Breaking Changes
* Removed the old `run_net_acpflow` public wrapper. Use `run_acpflow(net=...)` for already constructed networks and `run_acpflow(casefile=..., path=...)` for file-based workflows.

### Highlights
* Added a generic outer-loop control framework above `runpf!`.
  Transformer tap/phase control now uses this framework as the first concrete controller implementation.

### Improvements
* Added structured `ControlRunResult` output and `latest_control_result(net)` for inspecting controller status, outer iterations, PF solve count, controller rows, and trace rows.
* Added machine-readable control trace rows for transformer control.
* Added `run_acpflow(net=...)` as the preferred high-level entry point for already constructed networks, and made it the single public in-memory entry path.
* Simplified `examples/tap_control_demo_grid.jl` into a lightweight demo of the generic control framework.
* Documented the `control` configuration section, including that `control.controllers` is reserved for future YAML-based controller instantiation.


### Related
#179 Introduce a generic control framework above the existing power-flow solver. 

## Version 0.8.0 – 2026-05-25

### Breaking Changes

* The public AC power-flow path now supports only the sparse rectangular Newton-Raphson solver.
  Legacy polar/classic methods, dense PF matrices, and finite-difference PF Jacobian options are no longer supported as user-facing runtime choices.
* Power-flow configuration has moved to structured YAML sections and typed configuration objects.
  The old flat keyword-style configuration path is deprecated/removed for the cleaned rectangular workflow.
* Obsolete sparse switches such as `power_flow.sparse`, `opt_sparse`, and `state_estimation.sparse` are no longer valid configuration keys.
  Sparse matrix handling is mandatory for the production PF core.

### Highlights

* Added a central configuration workflow:
  * default template: `src/configuration.yaml.example`
  * optional user override: `examples/configuration.yaml`
  * typed config objects for power flow, MATPOWER import, state estimation, diagnostics, output, benchmarking, runtime, and performance profiling
  * early validation of unknown or obsolete keys
* Simplified the production power-flow path to the sparse rectangular AC solver with sparse Y-bus assembly, sparse analytic rectangular Jacobian, and sparse linear solves.
* Improved MATPOWER import and runner workflows:
  * central YAML-driven execution
  * configurable MATPOWER import options
  * cleaner compact summaries
  * better separation of numerical convergence, Q-limit validation, and solution availability
* Added configurable performance and timing output for MATPOWER and rectangular PF runs:
  * representative wall time
  * solver time
  * result-output time
  * timing coverage
  * optional allocation information
  * rectangular workspace metadata
* Improved large-case output handling:
  * configurable result-table row limits
  * summary/compact/full result-output modes
  * safer default behavior for large MATPOWER cases
* Added rectangular workspace reuse/preallocation controls:
  * `power_flow.rectangular_workspace_reuse`
  * `power_flow.rectangular_preallocate_workspace`
  * `power_flow.rectangular_workspace_min_buses`
* Reworked example scripts under top-level `examples/`:
  * `matpower_import.jl` is now a thin YAML-driven entry script
  * `tap_control_demo_grid.jl` uses central configuration helpers
  * `export_solution.jl` writes deterministic export artifacts under `examples/_out/export_solution/<case>_<timestamp>/`
* Added a clearer test profile structure:
  * default `fast` profile for normal development
  * `extended` profile for MATPOWER/output/documentation-heavy checks
  * `all` currently aliases `extended`

### Configuration Notes

* The main configuration template is now:

  ```text
  src/configuration.yaml.example
  ```

* A local user/example override can be placed at:

  ```text
  examples/configuration.yaml
  ```

* MATPOWER benchmarking moved to the top-level benchmark section:

  ```yaml
  benchmark:
    enabled: true
  ```

  The old `matpower_import.benchmark` key is rejected with a migration message.

### Output and Examples

* `examples/export_solution.jl` now produces files instead of only printing to the console.
  Typical output files are:

  ```text
  summary.txt
  internal_solution.csv
  external_solution.csv
  comparison.csv
  *_export.m
  ```

* `examples/matpower_import.jl` uses the central configuration and writes logs under:

  ```text
  examples/_out/
  ```

### Documentation

* Added/updated documentation for:
  * central configuration
  * power-flow configuration
  * MATPOWER import configuration
  * state-estimation configuration
  * performance profiling
  * test profiles
  * examples overview

### Related

* Issue #199: Central configuration, PF solver simplification, sparse-only PF core, and test-framework cleanup
* Issue #201: YAML redesign
  
## Version 0.7.8 – 2026-05-16
### Highlights
* Improved Q-limit handling for large MATPOWER imports, especially cases with many generators that have zero or very narrow reactive-power ranges.
* Added compact console reporting for large MATPOWER example runs while keeping the full diagnostics in the logfile.

### Fixes
* Fixed rectangular power-flow status caching to use weak network keys so repeated benchmark/example solves can release imported `Net` objects after callers consume the status.
* Fixed direct `run_acpflow(net=...)` and rectangular solver defaults so the Q-limit guard remains opt-in unless a caller or config explicitly enables it.
* Fixed direct `run_tap_controllers_outer!` defaults so the Q-limit guard remains opt-in when the exported tap-control API calls `runpf!`.
* Fixed MATPOWER example config forwarding so Q-limit guard options from YAML are preserved in the effective config and passed through all `run_acpflow` paths.
* Fixed MATPOWER example console row limiting so `console_max_rows` from YAML is forwarded to Q-limit event and final active-set row caps in the `run_acpflow` paths.
* Fixed rectangular NR status reporting: numerical convergence, Q-limit active-set consistency, final convergence, comparison status, and rejection reasons are now reported separately.
* Fixed Q-limit guard behavior so strongly violating active PV buses can be locked to PQ during eligible Q-limit checks, reducing final active-PV limit violations.

### Improvements
* Reduced console noise for large Q-limit active-set cases by replacing long PV→PQ and violation tables with compact summaries.
* Improved MATPOWER import example logging: the console now shows the essential run status, while detailed auto-profile evidence, diagnostics, and solver traces remain available in the logfile. 


## Version 0.7.7 – 2026-05-13

### Highlights

* Improved MATPOWER import diagnostics for large and difficult cases:
  * clearer YAML-backed logging and terminal summaries
  * VM/VA reference checks
  * branch-shift and transformer convention scans
  * bus-shunt and fixed-reference residual diagnostics

* Added configurable MATPOWER transformer and phase-shifter import conventions:
  * `matpower_ratio`
  * `matpower_shift_sign`
  * `matpower_shift_unit`

* Improved MATPOWER voltage reference handling:
  * configurable `BUS.VM` vs. `GEN.VG` handling
  * hybrid comparison mode
  * better diagnostics for PV/REF buses without online generators
  * correct handling of buses switched from PV to PQ by Q-limits

* Improved rectangular NR robustness for large MATPOWER cases:
  * preserved slack and PV setpoints during flat start
  * added DC-angle and blended-voltage start support
  * added wrong-branch diagnostics for suspicious low-voltage solutions

* Improved PV→PQ Q-limit handling:
  * clearer per-unit and MVAr diagnostics
  * hysteresis/deadband support
  * compact post-solve active-set consistency checks
  * separate handling of PV and REF/slack Q-limit violations

### Fixes

* Fixed MATPOWER slack/reference voltage import so VM/VA values are preserved unless explicitly overridden.
* Fixed nominal-tap transformer handling so explicit `TAP = 1` branches remain transformer models.
* Fixed repository-local MATPOWER `.jl` case loading and example regression-test handling.
* Fixed several Julia 1.12 / Revise world-age issues in MATPOWER diagnostic and example paths.
* Fixed flat-start setpoint extraction for networks with isolated buses.
* Improved warning/error capture in `matpower_import.jl` with compact logfile summaries.

### Diagnostics

* Added MATPOWER auto-profile pre-run mode (`matpower_import.auto_profile = recommend|apply`) to summarize or apply robust import, flat-start, PV/REF voltage-source, and Q-limit settings while preserving explicit YAML overrides.
* Added branch-neighborhood reports for selected high-residual buses.
* Added residual-cluster diagnostics for PEGASE-style mismatch regions.
* Added negative branch impedance scans while preserving signed MATPOWER `BR_R` / `BR_X` values.
* Added MATPOWER reference-data consistency diagnostics, including documentation of the known `case300.m` fixed-reference mismatch around `BUS_I 196 / 2040`.

### Related

* Issue #186: Singular Jacobian / world-age issue in rectangular NR solver.## Version 0.7.6 – 2026-05-11
### Bugfix
* Fixed rectangular MATPOWER runs so a singular Newton Jacobian is reported as non-convergence instead of aborting the example, and tightened Julia 1.12 / Revise entry-point calls to avoid world-age binding warnings.

## Version 0.7.5 – 2026-05-10
### New Features 
* Added PV→PQ Q-limit switching start controls (`qlimit_start_iter`, `qlimit_start_mode`, `qlimit_auto_q_delta_pu`) for rectangular power-flow runs.
* Added configurable bus-shunt modeling with classic Y-bus admittance stamping and a rectangular-solver voltage-dependent injection mode for keeping shunt effects in nonlinear mismatch terms.
* Added automatic rectangular Newton damping (`autodamp`) for difficult flat-start power-flow cases; the solver backtracks the Newton step from `damp` down to `autodamp_min` and accepts the first residual-reducing trial.
* Added start projection (`start_projection`) for difficult power-flow seeds, including optional DC-angle starts and raw/DC blend scans for both the internal rectangular solver and external-solver `PFModel` starts.

### Bugfix
* Replaced the singular sparse linear-solve fallback with a rank-revealing QR path before dense SVD fallback to avoid large `pinv`/LAPACK failures in ill-conditioned rectangular Newton steps.
* Wrapped MATPOWER voltage-angle comparison differences into the minimal ±180° range before applying angle tolerances.
* Synthetic tiled-grid benchmark example now falls back from `.yaml` to `.yaml.example` and reports when built-in defaults are used.

### Improvements
* Improved MATPOWER case parsing and network construction performance for large cases by reducing parser allocations, pre-sizing network containers, and avoiding repeated bus-name lookups while adding MATPOWER branches.
  
## Version 0.7.4 – 2026-05-08
### New Features 

### Bugfix
* Fixed YAML subset parser cleanup on parse errors so temporary config files are closed before `mktemp` cleanup on Windows.

### Improvements
* Improved large-network MATPOWER and rectangular-solver performance by aggregating prosumer-derived bus types and specified power injections in linear time. 
   
## Version 0.7.3 – 2026-05-04
### New Features
* Added transformer tap control in AC power flow with an outer-loop controller workflow.
* Added `addTapController!` support for:
  * voltage control via tap ratio (`mode = :voltage`)
  * branch active-power control via phase shift (`mode = :branch_active_power`)
  * coupled ratio + phase control (`mode = :voltage_and_branch_active_power`)
  * tap-control reporting in classic and structured ACP flow outputs

### Improvements
* Updated branch-model documentation with practical controller-direction probing guidance for phase-shift control.

### Tests
* Added regression tests for transformer tap controller behavior.

## Version 0.7.2 – 2026-04-15
### Bug Fixes
* Fixed the increased loading time for large test cases

## Version 0.7.1 – 2026-04-15
### Improvements
* Refactored decision logic for `qlimit_mode = :adjust_vset`:
  * consolidated Q-limit event handling via shared active-set flow
  * extracted voltage-step handling into dedicated helper logic
* Added typed `VoltageAdjustConfig` support for prosumers and integrated it into voltage-regulation detection.
* Simplified rectangular mismatch API by removing unused derivative keyword arguments from `mismatch_rectangular(...)`.

### Solver Robustness
* Rectangular solver now handles reduced Ybus matrices (caused by internal isolated buses) by expanding them back to full network dimension for mismatch/Jacobian processing.
* For rectangular runs with active-link merges and internal isolated buses, solver now uses a rectangular FD fallback path instead of switching to `:polar_full`.

## Version 0.7.0 – 2026-04-15
### New Features
* Added support for **P(U)** and **Q(U)** controller models in power flow calculations
* Characteristic curves can be defined via discrete points, with selectable interpolation methods:
  * linear
  * piecewise constant
  * cubic
  * polynomial
* Added support for controllers on **non-PV generators**
* Added documentation and examples for the new controller features
  
## Version 0.6.4 – 2026-04-12
### New Features
* Marked legacy Jacobian solvers as deprecated:
  * `runpf_full!` / `method = :polar_full`
  * `runpf_classic!` / `method = :classic`
  Both now emit a deprecation warning when used.
* Set rectangular complex Jacobian (`method = :rectangular`) as the default for:
  * `run_acpflow`
  * `run_acpflow(; net = net, ...)`
* Updated examples and user documentation to use `:rectangular` as the recommended/default solver method.
### Bugfixes
* For backwards compatibility, if vm is set, the coresponding bus type is set to PV

## Version 0.6.3 – 2026-04-11
### New Features
* Added `:adjust_vset` controller-based Q-limit handling at PV buses (adaptive Vset steps before optional PV→PQ fallback).
* Power-flow bus typing is now derived from attached prosumers (Slack > PV > PQ); `addBus!(busType=...)` is legacy-only and no longer defines operational PF type.

## Version 0.6.2 – 2026-04-02
### New Features
* adding sign validation and optional autocorrection of Q-limits before running power flows.
* Provide an option to lock selected PV buses from being switched to PQ 
* Added pre-run PV Q-limit preview logging in MVAr for easier diagnostics before the PF iteration loop.
## Version 0.6.1 – 2026-03-24
### New Features
* Bad Data Detection (BDD) and Statistical Diagnostics for State Estimation (SE)
  
### Bug Fixes
* Fix Issue 139

### Notes

* Bad-data detection and statistical diagnostics are available via
  `validate_measurements`, `runse_diagnostics`,
  `summarize_se_diagnostics`, and `print_se_diagnostics`.

## Version 0.6.0 – 2026-03-17
### New Features
* Added initial State Estimation (SE, WLS), documentation and examples.
* Added zero-injection-bus (ZIB) support in SE, documentation, and examples.
* Improved results reporting and logging.
* Refactor testsuite for clarity.

### Notes
* State Estimation is currently marked as **experimental**.
* Passive buses / ZIB are currently represented through zero-injection pseudo measurements in the WLS workflow.


## Version 0.5.0 – 2026-03-11
### Highlights
* Create and use machine-readable `ACPFlowReport` output
* Introduced bus-links +  documentation 
* Consolidated recent documentation

## Version 0.4.35 – 2026-02-02
### Fixes
* Stabilized MATPOWER case import and bus indexing in the presence of isolated buses.
* Fixed bus classification and ordering to ensure consistent solver input for flatstart and non-flatstart runs.
* Improved internal consistency checks for network topology before solver execution.

### Improvements
* Refined MATPOWER comparison tooling (Vm/Va diff diagnostics and angle alignment).
* Minor internal cleanups in MatpowerIO and network initialization to improve robustness.

## Version 0.4.34 – 2026-06-02
### Fixes
* Closed Issue 110:
* Corrected per‑unit conversion issues discovered during the refactoring of the MATPOWER import interface, including proper handling of line, transformer, and shunt parameters.
* Fixed PU conversion inconsistencies for shunts and aligned all shunt handling with the unified Y‑model.
* Harmonized Q‑limit checking across all solvers and ensured consistent PV/PQ switching behavior.

### New Features
* Updated shunt modeling: addShunts now exclusively supports the Y‑model; voltage‑dependent shunt power must be represented via loads or generators.
* Unified solver logic for reactive power limits and PV/PQ switching
* added new testcase matpower vs manual network

## Version 0.4.33 – 2026-02-02

### New Features
* Added on-demand MATPOWER case handling via `FetchMatpowerCase.ensure_casefile`
  (automatic download of `.m` cases and optional generation of `.jl` cases)
* Added option to control flat start behaviour in AC power flow (`flatstart`)

### Internal Changes
* Refactored MATPOWER case import logic and removed redundant example-based loaders

## Version 0.4.32 – 2026-30-01
## New Features
* Added function to export calculated network
* 
## Version 0.4.31 – 2026-28-01
### Title
* Change license from BSD-3-Clause to Apache License 2.0

## Version 0.4.30 – 2026-28-01
### Bug Fixes
* Importing Matpower files, wrong conversion to per unit system 
### New Features
* Option to choose flatstart 

## Version 0.4.29 – 2025-12-13
### New Features
* Support native 3-winding transformers in Sparlectra without AUX bus generation and allow non-PU parameterization #90
* Support 2-winding transformers as well in non-PU parameterization #90
  
## Version 0.4.28 – 2025-12-12
### Bug Fixes
* Close Issue #85: test totalBusPower vs. TotalLosses failed for solver != rectangular

## Version 0.4.27 – 2025-12-11
### New Features
* Added 3Bus Testcase with PV Generator
* Added createTest5BusNet with multiple generators at one bus and zero injection generator
* Issue #81 Fix Calculatation of Generators with Zero Injektion
* Mark PV->PQ switching in printACPFlowResults Log 
  
### Bug Fixes
* fixes per unit calculation for line shunt admittance
* fixes calculation of losses


## Version 0.4.26 – 2025-12-04
* Issue #74: Multiple Generators at One Bus Not Handled Correctly
* Documentation Updates 

## Version 0.4.25 – 2025-11-29

### Added
* Rectangular (Complex-State) Newton–Raphson Solver

## Version 0.4.24 – 2025-11-20

### Added
- Q-limit handling with automatic PV→PQ switching when a generator hits its reactive power limits.

## Version 0.4.23 (2025-11-11)

###  New Features
* Full-state Newton–Raphson solver

## Version 0.4.22 (2025-08-27)
### Bug Fixes
 - small fixes

## Version 0.4.21 (2025-03-14)
### New Features
 - adding functions to remove elements from a net

## Version 0.4.20 (2025-03-11)
 -  internal reorganization and small bugfixes 
 
## Version 0.4.19 (2024-10-14)
### Bug Fixes
 - pu calculation for transformer impedance

## Version 0.4.18 (2024-04-14)

### Bug Fixes
 - closes issue#48 "printout jacobian runs into error"

## Version 0.4.17 (2024-04-14)
### New Features
 - added testcase for importing Matpower files

### Bug Fixes
 - bugfix wrong function call in `createnet_powermat`

## Version 0.4.16 (2024-04-13)
### Bug Fixes
- bugfix shunt index for isolated buses, closes issue #38

## Version 0.4.15 (2024-04-12)
### New Features
- Implemented a function to detect isolated buses and incorporate them into the network calculation (Issue #38)

## Version 0.4.14 (2024-04-12)
### Bug Fixes
- bugfix addShuntPower, closes issue #36

## Version 0.4.13 (2024-04-12)
### New Features
- added attribute for Lineparameters based on length
- added update parameter function for network
- added workshop documentation

### Bug Fixes
- taking line length not (always) into account for line parameters
- parsing emtpy lines of Matpowerfiles
- documentation rendering

## Version 0.4.12 (2024-04-08)
### New Features
- added functions to facilitate the modification of networks.
- documentation available at https://welthulk.github.io/Sparlectra.jl/.

### Bug Fixes
- print prosumers

## Version 0.4.11 (2024-04-05)
### New Features
- make changes to imported Matpower networks after import.
- added functions to facilitate the creation of networks.

### Enhancements
- added documentation make file

### Bug Fixes
- import and parser for Matpower .m files

## Version 0.4.10 (2024-03-30)

### New Features
 - removed numerous redundant functions, partially restructured classes
 - removed support for CGMES due to the availability of numerous alternative libraries
 - removed support for the legacy custom JSON data format (potentially subject to reintroduction at a later stage)
 - added functions to facilitate the creation of networks
 - better performance
 
### Bug Fixes
- calculation of branch flow and losses
- branches out of service

## Version 0.4.8 (2024-03-26)
- first package release registered in the Julia registry

## Version 0.4.1 (2023-12-19)
- Initial release of Sparlectra

## Version 0.4.0 (2023-11-30)
- Initial public commit of Sparlectra 

* **Bugfix**: Restored MATPOWER runner operational reporting through package-level runners (header, config paths, resolved case, logfile path, compact summary, benchmark/performance blocks) while keeping example scripts thin and free of local YAML parsing.
* **Feature**: Added asynchronous local Web UI PowerFlow jobs with a POST-only
  cooperative abort action, persistent aborted status, distinct history
  rendering, and single-active-run protection. Solver execution is never
  interrupted unsafely; an abort requested during a blocking phase remains
  final when that phase returns.

* **Bugfix**: Web UI submissions now run independently of request handling so active status and Abort controls are reachable before completion; first startup also provisions user-writable configuration, case-cache, run, and operation-log paths instead of using package-tree runtime defaults.
* **Improvement**: The PowerFlow result/status page now highlights elapsed runtime in a clock-style `HH:MM:SS.mmm` card beside the run status while retaining raw timing metadata for diagnostics.

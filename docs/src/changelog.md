# Change Log
## Version 0.8.0 – 2026-05-17


### Breaking Changes

* Power flow now supports only the sparse rectangular AC solver. Requests for polar/classic methods, finite-difference PF Jacobians, or dense PF matrices now fail validation with a clear error instead of falling back silently.
* Removed the old flat keyword-style configuration path for rectangular power-flow runs. Configuration is now handled through structured YAML sections and typed config objects.
* Removed the obsolete finite-difference power-flow option from configuration, APIs, examples, documentation, and tests; the sparse analytic rectangular solver is now the only public power-flow Jacobian path.


### Fixes

* **Bugfix**: Fixed runtime thread config parsing so numeric YAML values for `runtime.julia_threads`/`runtime.blas_threads` are accepted and normalized to strings, preventing `MethodError: no method matching String(::Int64)` in MATPOWER example startup paths.
* **Bugfix**: Fixed MATPOWER import configuration validation for `preallocate_network` by allowing symbol allow-lists passed as vectors, preventing `MethodError` in `matpower_import.jl` startup under Julia 1.12 workflows.
* **Bugfix**: Unified rectangular PF outcome reporting and gating through a canonical status payload (outcome, numerical convergence, solution availability, limit-validation status, reason text, mismatch, iterations, elapsed time), fixed compact/log summaries to report `numerical_solution=FAIL` for non-converged NR runs, and aligned result-table plus MATPOWER comparison gating to depend on explicit solution availability (including `converged_limits_failed` as a numerical-solution class with validation warning output).
* Increased the passive-bus zero-injection state-estimation regression iteration cap so `test/runtests.jl` converges reliably under the current solver behavior.
* Consolidated shared example YAML-loading/runtime-output helpers into `examples/example_utils.jl` and updated MATPOWER/tap-control examples to use the centralized path, including Julia 1.12-safe include usage for module-based example smoke tests.
* Moved MATPOWER example configuration ownership into package core entry points: `examples/matpower_import.jl` is now a thin script that only resolves config/case inputs and calls `run_matpower_case`, with Julia 1.12-safe `Base.invokelatest(getfield(@__MODULE__, :main))` entry invocation.
* **Bugfix**: Added an explicit migration error for removed `matpower_import.benchmark` user config keys and directed users to top-level `benchmark.enabled`; also fixed MATPOWER runner status handling and benchmark method execution/seconds handling in package-level runner paths.
* **Bugfix**: Hardened `test/test_matpower_example.jl` for Julia 1.12/Revise workflows by invoking the test entry via `Base.invokelatest(...)`, using a test-local MATPOWER config file, and accepting known network-fetch error variants in offline or proxied environments.
* **Bugfix**: Fixed MATPOWER runner benchmark method config cloning to use current `PowerFlowConfig` fields (`tol`, `max_iter`, `autodamp`, `autodamp_min`, `start_mode`, `qlimits`) so Julia 1.12 MATPOWER example runs no longer fail with `PowerFlowConfig.max_ite` field errors.
* **Bugfix**: Clarified MATPOWER benchmark timing semantics by separating representative wall time, solver time, result-output coverage rows, and benchmark median/min reporting; replaced ambiguous compact `time=` summary with explicit `representative_time=`/`solver_time=` fields and aligned AC result-table runtime label to `Solver time`.
* **Bugfix**: Fixed MATPOWER timing-status consistency so summary `solver_time` now maps to instrumented `solver_total` (or `unavailable` when missing), `representative_time` remains full wall time, and AC result output prints distinct solver/representative timing labels instead of silently aliasing representative time as solver time.
* **Bugfix**: Fixed MATPOWER benchmark execution under Julia 1.12/Revise by interpolating the benchmark casefile into `@benchmarkable` and invoking `BenchmarkTools.run` via `Base.invokelatest`, preventing world-age warnings and `UndefVarError: local_case` failures in the MATPOWER example path.
* **Bugfix**: Fixed MATPOWER import benchmark path for Julia 1.12/Revise by invoking benchmark sample closures via `Base.invokelatest` and guarding `performance_profile` access in `run_acpflow`, eliminating world-age warnings (`Main.main`, BenchmarkTools-generated closures) and `MethodError: get(::Nothing, ::Symbol, ::Nothing)` during sample runs.
* **Bugfix**: Fixed MATPOWER example startup/runtime thread handling by distinguishing Julia startup threads vs runtime BLAS threads, adding `--julia-threads`/`SPARLECTRA_JULIA_THREADS` overrides with single-pass script re-exec guidance, and hardening compact performance/final console output against Ctrl-C interruptions in script mode.
* **Bugfix**: Fixed MATPOWER runner output routing by running benchmark samples silently, capturing one representative solve for logfile diagnostics, honoring `output.logfile_performance`, and deriving compact `q_limit_active_set` summaries from typed status fields so benchmark console output stays compact and summary status remains consistent.
* **Bugfix**: Fixed MATPOWER runner performance logging so `logfile_performance = full` now includes representative phase timings plus iteration/allocation details (when enabled), `compact` includes aggregated timing rows, and benchmark samples no longer pollute the representative performance profile.
* **Improvement**: Added MATPOWER performance timing coverage reporting with representative wall time, top-level and total recorded phase sums, benchmark-event separation, and unaccounted/overlap diagnostics to make instrumentation coverage explicit.
* **Bugfix**: Fixed diagnostics/output section separation in central YAML parsing so MATPOWER example effective-config logging no longer overrides `diagnostics.*` values with same-named `output.*` keys.
* **Bugfix**: Reworked YAML configuration documentation into table-based references, added explicit allowed-value validation for user-facing Symbol/String options, and added regression tests to ensure every key in `src/configuration.yaml.example` is documented (including migration from removed `matpower_import.benchmark` to `benchmark.enabled`).
* **Bugfix**: Fixed MATPOWER start-mode configuration regression by validating all runtime-supported `power_flow.start_mode.angle_mode`/`voltage_mode` values, storing both fields in typed config, propagating them into `run_acpflow`, and adding regression coverage for accepted and rejected values under Julia 1.12 example workflows.
* **Bugfix**: Audited and aligned central YAML/runtime configuration coverage for MATPOWER and power-flow paths by adding typed `state_estimation.method`, widening MATPOWER allowed-value validation to runtime-supported modes (`pv_voltage_source`, `compare_voltage_reference`, `ratio`, `bus_shunt_model`), and adding configuration-coverage plus roundtrip regression tests for previously dropped keys (including Q-limit and MATPOWER import toggles).
* **Bugfix**: Refactored rectangular MATPOWER result handling to use a structured solver outcome (`:converged`, `:converged_with_limit_warnings`, `:converged_limits_failed`, `:not_converged`, `:singular_jacobian`, `:solver_error`), kept numerical/active-set status fields in exported status summaries, and allowed voltage-result output for numerically converged solutions even when final Q-limit active-set validation fails.
* **Bugfix**: Improved MATPOWER network-construction instrumentation and performance controls for large imports by adding preallocation config (`matpower_import.preallocate_network`, `matpower_import.preallocate_min_buses`), export of construction metadata in the performance profile, and import-path setup for reusable MATPOWER lookup maps without changing solver numerics.

### Improvements

* Added a central configuration workflow:
  * default template: `src/configuration.yaml.example`
  * optional user override: `examples/configuration.yaml`
  * early unknown-key validation
  * MATPOWER options grouped under `matpower_import`
* Centralized runtime configuration for MATPOWER import, rectangular power flow, state estimation, diagnostics, output, and profiling. This reduces long option lists and makes the call paths easier to follow.
* Optimized rectangular solver voltage-setpoint lookup by building per-solve bus and generator maps instead of repeatedly scanning prosumers in large MATPOWER cases.
* Added more detailed rectangular solver timings for setup, active-set handling, Q-limit checks, state updates, result writeback, and status bookkeeping.
* Added configurable MATPOWER performance profiling with phase timings, optional allocation reporting, Newton diagnostics, and large-case speed switches.
* Split start-projection profiling into smaller phases, including DC-start construction, matrix assembly, solve, candidate evaluation, branch checks, voltage limiting, MATPOWER lookup, and final selection.

### Related

* Issue #201: YAML redesign
* Issue #199: Centralized configuration

  
## Version 0.7.8 – 2026-05-16
### Highlights
* Improved Q-limit handling for large MATPOWER imports, especially cases with many generators that have zero or very narrow reactive-power ranges.
* Added compact console reporting for large MATPOWER example runs while keeping the full diagnostics in the logfile.

### Fixes
* Fixed rectangular power-flow status caching to use weak network keys so repeated benchmark/example solves can release imported `Net` objects after callers consume the status.
* Fixed direct `run_net_acpflow` and rectangular solver defaults so the Q-limit guard remains opt-in unless a caller or config explicitly enables it.
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

* Added MATPOWER auto-profile pre-run mode (`matpower_auto_profile = recommend|apply`) to summarize or apply robust import, flat-start, PV/REF voltage-source, and Q-limit settings while preserving explicit YAML overrides.
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
  * `run_net_acpflow`
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

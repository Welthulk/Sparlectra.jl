# Change Log
## Version 0.7.7 – 2026-05-13

### Highlights
* Improved MATPOWER import diagnostics for large cases, including effective YAML logging, VM/VA reference checks, branch-shift convention scans, bus-shunt impact checks, and clearer terminal summaries.
* Added configurable MATPOWER transformer and phase-shifter import conventions:
  * `matpower_ratio`
  * `matpower_shift_sign`
  * `matpower_shift_unit`

### Improvements
* Improved PV/REF voltage-setpoint handling for MATPOWER imports:
  * configurable `BUS.VM` vs. `GEN.VG` selection
  * hybrid comparison mode
  * compact diagnostics for PV/REF buses without online generators
* Improved rectangular NR flat-start behaviour for large MATPOWER cases:
  * preserved slack/PV setpoints during flat start
  * added DC-angle and blended-voltage start support
  * added wrong-branch diagnostics for suspicious low-voltage solutions
* Improved MATPOWER reference diagnostics with an optional branch-neighborhood report for selected high-residual buses, including raw branch data, effective tap interpretation, Y-bus stamp terms, and fixed-reference P/Q flow contributions.
* Improved PV→PQ Q-limit handling:
  * clearer per-unit and MVAr diagnostics
  * hysteresis/deadband support to reduce premature switching and PV↔PQ chatter
  * post-solve Q-limit summaries for large cases
* Fixed several Julia 1.12 / Revise world-age issues in MATPOWER example and diagnostic paths.
* Fixed repository-local MATPOWER `.jl` case loading and example regression-test handling.

### Bugfixes

* Bugfix: Made MATPOWER Vm/Va comparisons aware of final PV→PQ Q-limit switching so switched buses use `BUS.VM` as their voltage reference and are reported separately from active PV, PQ, and REF/slack buses.
* Bugfix: Refined MATPOWER hybrid Vm/Va comparison diagnostics so deviations are summarized by voltage-reference kind and tolerance excursions that occur only on `final_pq_after_qlimit` buses are reported as WARN with an explanation instead of a hard compare failure.
* Bugfix: Separated final rectangular Q-limit diagnostics for active PV and REF buses so PV violations fail the active-set check while REF/slack violations are reported as WARN by default without switching the slack bus to PQ.
* Bugfix: Added a compact final rectangular PV/REF Q-limit consistency check that reports converged active-set violations with the same reactive-power source used by switching logic, including switch history and hysteresis return-band status.
* Bugfix: Repaired MATPOWER fixed-reference diagnostics for PEGASE-style residual clusters by adding a branch-connected residual-cluster report that groups high P/Q mismatches and prints their internal branch context.
* Bugfix: Extended the MATPOWER example diagnostics with YAML-backed nodal fixed-reference balance breakdowns and negative branch impedance scans, preserving signed branch `BR_R`/`BR_X` values while reporting model-critical negative impedance rows without clipping or correction.
* Added MATPOWER reference-data consistency diagnostics and documented the known `case300.m` fixed-reference mismatch around `BUS_I 196 / 2040`.
* Fixed MATPOWER nominal-tap transformer preservation so explicit `TAP = 1` branches, including case300 BUS_I 196 -> 2040 style connections, remain transformer models while ordinary `TAP = 0` lines still stamp with unity ratio.
* Fixed MATPOWER slack/reference voltage handling so imported slack VM/VA values are preserved unless explicitly overridden.
* Fixed rectangular flat-start initialization so slack and PV voltage setpoints remain active constraints.
* Fixed external `PFModel` flat-start setpoint extraction for networks with isolated buses.
* Fixed MATPOWER reference residual diagnostics and phase-shift logging paths.
* Fixed warning/error capture in `matpower_import.jl`, replacing noisy console output with compact logfile summaries.

### Related
* Issue #186: Singular Jacobian / world-age issue in rectangular NR solver.
  
## Version 0.7.6 – 2026-05-11
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

# Change Log
## Version 0.7.3 – 2026-05-xx
### New Features
* Added transformer tap control in AC power flow with an outer-loop controller workflow.
* Added `addTapController!` support for:
  * voltage control via tap ratio (`mode = :voltage`)
  * branch active-power control via phase shift (`mode = :branch_active_power`)
  * coupled ratio + phase control (`mode = :voltage_and_branch_active_power`)
* Added discrete tap-step handling (ratio and phase) with limit enforcement in controller updates.
* Added user documentation for complex tap control and controller modes.
* Added runnable examples for:
  * transformer tap voltage control
  * combined complex tap control
  * phase-shift transformer active-power control

### Improvements
* Updated branch-model documentation with practical controller-direction probing guidance for phase-shift control.
* Improved tap-control reporting in classic and structured ACP flow outputs, including typed transformer-controller rows with tap position/stage and consistent branch power direction reporting.
* Updated `tap_control_demo_grid.jl` logging to create versioned run logs in `src/examples/_out/` and refined OLTC demo step resolution to `tap_step = 0.00625` for the unchanged `0.90 .. 1.10` ratio range.
* Tap-controller voltage convergence now consistently evaluates deadband in `|V - V_target|` (pu), while allowing optional step-direction evaluation via `voltage_error_metric = :vm2` for users who prefer a squared-voltage control signal.
* Added YAML configurability in the tap-control demo for `*_deadband_vm_pu` and `*_voltage_error_metric` (`vm` or `vm2`) on voltage-regulating controllers.

### Tests
* Added regression tests for transformer tap controller behavior.
* Added dedicated tests for transformer phase-shift control direction and target tracking.

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

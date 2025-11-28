# Change Log

## **Version 0.4.25 – 2025-11-28**

### **Added**

* **Rectangular (Complex-State) Newton–Raphson Solver**

  * New solver path via `runpf_rectangular!` and `run_complex_nr_rectangular_for_net!`.
  * Supports PQ, PV and Slack buses using *rectangular complex voltage state* (Vr, Vi).
  * Supports both **analytic** Jacobian (Wirtinger-based) and **finite-difference** Jacobian.
  * Optional **sparse analytic Jacobian** via `build_rectangular_jacobian_pq_pv_sparse`.
* **PV→PQ active-set logic** now implemented also for the rectangular solver.

  * Automatically clips Q to limits and switches bus type during iteration.
  * Writes Q values back for former PV buses (same behavior as full/polar solver).
* **Unified solver options**

  * New flags: `opt_fd`, `opt_sparse`, `method=:rectangular`.
  * Integrated in `run_acpflow` and `run_net_acpflow`.

### **Changed**

* Rectangular solver now mirrors the behavior of `calcNewtonRaphson_withPVIdentity!` regarding Q-limits.
* Improved solver interface: shared options for sparse/FD variants across all NR solvers.
* Clean-up of deprecated code and migration of complex-state NR into `jacobian_complex.jl`.

### **Fixed**

* Sign inconsistency in rectangular Q-residuals after PV→PQ switching.
* Correct propagation of updated Q values into network structures after PV→PQ conversion.

### **Notes**

* Sparse rectangular Jacobian is supported but yields limited speedup (dense FD path still heavy).
* Full documentation update planned next.


## Version 0.4.24 – 2025-11-20

### Added
- Q-limit handling with automatic PV→PQ switching when a generator hits its reactive power limits.
- New field `isPVactive::Union{Nothing,Bool}` to track current PV/PQ status.
- Q-limit logging utilities (recording, resetting, inspecting limit hits).
- New Jacobian debugging helpers (`zero_a_row!`, `print_jacobian`).

### Changed
- Q-limits are now collected from prosumer/generator data and aggregated per bus.
- Newton–Raphson flow and Jacobian updated to use the new active-set Q-limit logic.
- API cleanup related to deprecated PV-generator Q-limit functions.


## Version 0.4.23 (2025-11-11)

###  **New Features**

* **Full-state Newton–Raphson solver** (`calcNewtonRaphson_withPVIdentity!`)

  * Solves for all non-slack bus states (θ, V) in a unified system.
  * PV buses use an identity constraint row enforcing `V − Vset = 0`.
* **Active Q-limit handling**

  * Automatic PV → PQ switching when reactive generation hits min/max limits.
  * Optional hysteresis band (`q_hyst_pu`) and cooldown iterations for stable re-enabling.
  * Limit events logged in `net.qLimitEvents` and accessible via `printQLimitLog()`.
* **New convenience wrappers**

  * `runpf_full!()` for running the full-system power flow.
  * Q-limit utilities: `buildQLimits!`, `resetQLimitLog!`, `getQLimits_pu()`.
* **Improved result tracking**

  * `Net` objects now record per-bus and total P/Q residuals after each iteration.
  * PV/Q-limit state changes visible in `net` result fields.


### **Technical Changes**

* PV identity rows replace Q-balance equations in the Jacobian.
* Robust fallback defaults for missing fields (`q_hyst_pu`, `cooldown_iters`).
* Q-limit clipping operates directly in p.u., using aggregated generator data.
* Improved compatibility with `printACPFlowResults()` and other result tools.


### **Known Limitations (to be addressed next)**

1. **Q-limits are not yet parameterizable at network creation.**
   They must still be constructed via `buildQLimits!()` after the network is loaded.
2. **MATPOWER import currently ignores Qmin/Qmax columns.**
3. **Jacobian still uses `sin()`/`cos()` formulations.**
   Future versions will use a rectangular (non-trigonometric) derivative form for speed and stability.
4. **Branch flow and loss post-processing** must still be called manually (e.g. `calcBranchFlow(net)`).

---

### **Developer Notes**
* Compatible Julia version: ≥ 1.11.1


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

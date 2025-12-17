# Change Log
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

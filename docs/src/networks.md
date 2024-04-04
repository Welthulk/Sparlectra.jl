API Network 
=============
#### Usage Guide for `Net` Module

#### Introduction
The `Net` module provides functionality for creating and manipulating power system network models in Julia. It includes features for defining buses, branches, transformers, prosumers, and shunts, as well as methods for running power flow analysis.

#### Overview of `Net` Module

The `Net` module consists of the following components:

- `Net` struct: Represents a power system network.
- Functions for adding components to the network:
  - `addBus!`: Adds a bus to the network.
  - `addShunt!`: Adds a shunt to a bus in the network.
  - `addACLine!`: Adds an AC line segment between two buses in the network.
  - `add2WTrafo!`: Adds a two-winding transformer between two buses in the network.
  - `addPIModellTrafo!`: Adds a transformer modeled as a PI-branch between two buses in the network.
  - `addProsumer!`: Adds a prosumer (generator or load) to a bus in the network.
- Validation function:
  - `validate`: Validates the network configuration.
- Utility functions:
  - `geNetBusIdx`: Gets the index of a bus in the network.
  - `getNetOrigBusIdx`: Gets the original index of a bus in the network.
  - `hasBusInNet`: Checks if a bus is present in the network.
  - `get_bus_vn_kV`: Gets the nominal voltage of a bus in the network.
  - `getBusType`: Gets the type of a bus in the network.
  - `lockNet!`: Locks the network to prevent further modifications.
  - `setTotalLosses!`: Sets the total losses in the network.
  - `getTotalLosses`: Gets the total losses in the network.


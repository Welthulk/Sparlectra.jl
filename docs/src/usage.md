# Getting Started

## Installation

```julia
using Pkg
Pkg.add("Sparlectra")
```

## Importing a Matpower Network Configuration

This example demonstrates how to import a network configuration from a Matpower case file and run a power flow analysis on it. The casefile is located in the `data` directory of the package. The `run_acpflow` function is used to run the power flow analysis.

### Example
```julia
using Sparlectra
using Logging
global_logger(ConsoleLogger(stderr, Logging.Warn))  # Set global logger to log to stderr with WARN level

file = "case3.m"
tol = 1e-8
ite = 10
verbose = 0 # 0: no output, 1: iteration norm, 2: + Y-Bus, 3: + Jacobian, 4: + Power Flow

# Call acpflow function with input parameters and measure execution time
run_acpflow(max_ite= ite,tol = tol, casefile=file)
```


> Note: Network Data
While contributions to the project are appreciated, please note that providing support for individualized network data issues is beyond the scope of this project, as it is not maintained by an organization. Users are encouraged to take initiative in resolving such issues independently and sharing their results with the community.





## Usage Guide for `Net` Module

### Introduction
The `Net` module provides functionality for creating and manipulating power system network models in Julia. It includes features for defining buses, branches, transformers, prosumers, and shunts, as well as methods for running power flow analysis.

### Overview of `Net` Module

The `Net` module consists of the following components:

- `Net` struct: Represents a power system network.
- Functions for adding components to the network:
  - `addBus!`: Adds a bus to the network.
  - `addShunt!`: Adds a shunt to a bus in the network.
  - `addACLine!`: Adds an AC line segment between two buses in the network.
  - `add2WTrafo!`: Adds a two-winding transformer between two buses in the network.
  - `addProsumer!`: Adds a prosumer (generator or load) to a bus in the network.
- Validation function:
  - `validate`: Validates the network configuration.
- Utility functions:
  - `geNetBusIdx`: Gets the index of a bus in the network.


### Creating a Network
```julia
# Import the Net module
  using Sparlectra
  using BenchmarkTools
  # Create a network
  # Slack bus B1 at 220 kV with 1.0 pu voltage and 0.0 pu angle
  # PQ bus B2 at 220 kV with 1.0 pu voltage and 0.0 pu angle
  # PQ bus B3 at 22 kV with 1.0 pu voltage and 0.0 pu angle
  # Shunt at bus B3 with 0.0 kW and 150.0 kVar
  # AC line from bus B1 to B2 at 220 kV with 100 km length, 0.0653 ohm/km resistance, 0.398 ohm/km reactance, 9.08 nF/km capacitance, and 0.0 power factor
  # 2W transformer from bus B2 to B3 with 1000 MVA rating, 13.0% voltage ratio, 0.28% resistance, 20.0 kW power factor, and 0.06% no-load current
  # Energy consumer at bus B3 with 285.0 kW and 200.0 kVar
  # Synchronous machine at bus B1 with 1.02 pu voltage and 0.0 pu angle
  #
  # G->1* ----L---- 2 --T-- 3->L
  #                 x
  #   
  net = Net(name = "testnet", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", busType = "Slack", vn_kV = 220.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "B2", busType = "PQ", vn_kV = 220.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "B3", busType = "PQ", vn_kV = 22.0, vm_pu = 1.0, va_deg = 0.0)
  addShunt!(net = net, busName = "B3", pShunt = 0.0, qShunt = 150.0)
  addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  add2WTrafo!(net = net, fromBus = "B2", toBus = "B3", sn_mva = 1000.0, vk_percent = 13.0, vkr_percent = 0.28, pfe_kw = 20.0, i0_percent = 0.06)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 285.0, q = 200.0)
  addProsumer!(net = net, busName = "B1", type = "SYNCHRONOUSMACHINE", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)

  # Run power flow
  tol = 1e-6
  maxIte = 10  
  etime = @elapsed begin
    ite, erg = runpf!(net, maxIte, tol, 0)
  end
  if erg != 0
    @warn "Power flow did not converge"        
  else
    calcNetLosses!(net)
    printACPFlowResults(net, etime, ite, tol)
  end
  
```
### Exporting a Matpower Network Configuration
Once you have created a network configuration, you can export it to a Matpower case file using the following function:
```julia
  writeMatpowerCasefile(myNet, filename)
```








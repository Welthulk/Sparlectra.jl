# Getting Started

## Installation
### First-time Julia Users
SPARLECTRA requires the Julia programming language. You can download Julia from [here](https://julialang.org/downloads/).

### Installing Program Files
To add and register the Sparlectra package, open the Julia REPL and press `]` to enter the Package Manager. Then, execute the following command:

```
add https://github.com/Welthulk/Sparlectra.jl
```	

The Julia package manager will automatically download the required packages.

> **Hint**: 
For development purposes, you can clone the package from the Git repository and place it in your user directory (e.g., `~/.julia/dev` or `%HOME%/.julia/dev`).

## Running the Power Flow Example using VS Code
Navigate to the `test` folder and open the `testparser.jl` file in VS Code. Press `CTRL+F5` to run the file. Please note that during the first execution, the program will be compiled, which may require additional time. The output will be displayed in the terminal.

# Usage Guide for `Net` Module

## Introduction
The `Net` module provides functionality for creating and manipulating power system network models in Julia. It includes features for defining buses, branches, transformers, prosumers, and shunts, as well as methods for running power flow analysis.

## Overview of `Net` Module

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

## Example Usage

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
## Importing a Matpower Network Configuration
You can import a network configuration from a Matpower case file using the following function:
```julia
  net = createNetFromMatPowerFile(filename)
```
For example, to import a network configuration from a file named `case2.m` located in the `data/mpower` directory, you can use the following code:
```julia
  filename = "case2.m"
  path = joinpath(pwd(), "data", "mpower", filename)  
  net = createNetFromMatPowerFile(path)
  ...
```
> Note: Not all Matpower case file features are supported.

## Exporting a Matpower Network Configuration
Once you have created a network configuration, you can export it to a Matpower case file using the following function:
```julia
  writeMatpowerCasefile(myNet, filename)
```







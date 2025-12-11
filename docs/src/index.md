# Sparlectra

Sparlectra is a Julia package for the simulation of electrical power systems. It primarily features a program for calculating load flow using the Newton-Raphson method. The focus is to provide valuable insights into load flow calculations for both students and ambitious professionals.

## Features

- **Power Flow Analysis**: AC power flow calculation using the Newton-Raphson method
- **Network Modeling**: Comprehensive modeling of electrical network components
  - Buses (PQ, PV, Slack types)
  - Transmission lines
  - Transformers (2-winding and 3-winding)
  - Generators and loads
  - Shunts
- **File Import/Export**: Support for Matpower (.m) file format
- **Network Modification**: Tools for creating, modifying, and analyzing networks
- **Analysis Tools**: Calculate losses, detect isolated buses, analyze branch flows

## Installation

For installation, run the following command in the Julia REPL:
```julia
import Pkg
Pkg.add("Sparlectra")
```

## Quick Start

Here's a simple example to get you started:

```julia
using Sparlectra

# Create a new network
net = Net(name = "example", baseMVA = 100.0)

# Add buses
addBus!(net = net, busName = "B1", busType = "PQ", vn_kV = 110.0)
addBus!(net = net, busName = "B2", busType = "PQ", vn_kV = 110.0)
addBus!(net = net, busName = "B3", busType = "Slack", vn_kV = 110.0)

# Add a transmission line
addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 10.0, r = 0.01, x = 0.1)

# Add a load
addProsumer!(net = net, busName = "B1", type = "ENERGYCONSUMER", p = 1.0, q = 0.2)

# Add a generator
addProsumer!(net = net, busName = "B3", type = "SYNCHRONMASCHINE", referencePri = "B3")

# Run power flow
ite, erg = runpf!(net, 10, 1e-6, 0)
if erg == 0
    calcNetLosses!(net)
    printACPFlowResults(net, 0.0, ite, 1e-6)
end
```

## Documentation Structure

- **[Changelog](changelog.md)**: Version history and updates
- **[Networks](networks.md)**: Creating and manipulating network models
- **[Branch Model](branchmodel.md)**: Details of the network branch model
- **[Import/Export](import.md)**: Importing and exporting network configurations
- **[Component Removal](remove_functions.md)**: Removing components from networks
- **[Workshop](workshop.md)**: Guided exercises and examples
- **[Function Reference](reference.md)**: Complete API documentation
- **[Powerlimit Guide](powerlimits_solvers.md)**: Handling of power limits

## Contributors

- [Udo Schmitz](https://www.linkedin.com/in/udo-schmitz-315536250/)

## License

Sparlectra is open-source software.

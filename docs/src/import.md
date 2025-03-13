# Import and Export Network Configuration

Sparlectra provides functionality to import network models from Matpower case files and export networks to Matpower format. This page documents the functions for reading, importing, and exporting network configurations.

## Importing Matpower Files

### Basic Import

The simplest way to import a Matpower case file is using the `run_acpflow` function, which both imports the network and runs a power flow analysis:

```julia
using Sparlectra
using Logging

# Import a Matpower case file and run power flow
net = run_acpflow(
    casefile = "case3.m",       # Case file name
    max_ite = 10,               # Maximum iterations
    tol = 1e-8,                 # Convergence tolerance
    verbose = 0,                # Verbosity level (0: none, 1: iterations, 2+: more details)
    printResultToFile = false   # Whether to print results to a file
)
```

### Import Without Power Flow

If you want to import a network without immediately running power flow, you can use the `createNetFromMatPowerFile` function:

```julia
using Sparlectra

# Path to the Matpower file
file_path = "path/to/case5.m"

# Import the network configuration
net = createNetFromMatPowerFile(file_path, false)  # Second parameter controls logging
```

### Import Parser

The `casefileparser` function parses Matpower case files and returns the raw data arrays:

```julia
using Sparlectra

file = "case9.m"

# Parse the file and get raw data
case_name, baseMVA, busData, genData, branchData = casefileparser(file)

# Now you can work with the raw data arrays
println("Case name: $case_name")
println("Base MVA: $baseMVA")
println("Number of buses: $(size(busData, 1))")
println("Number of generators: $(size(genData, 1))")
println("Number of branches: $(size(branchData, 1))")
```

## Exporting Networks to Matpower Format

You can export a Sparlectra network to a Matpower case file using the `writeMatpowerCasefile` function:

```julia
using Sparlectra

# First create or import a network
net = Net(name = "export_example", baseMVA = 100.0)
# Add components to the network...

# Export the network to a Matpower case file
filepath = "path/to/output/export_example.m"
writeMatpowerCasefile(net, filepath)
```

## Running Power Flow on Imported Networks

After importing a network, you can run a power flow analysis:

```julia
using Sparlectra

# Import a network from a Matpower file
net = createNetFromMatPowerFile("case5.m", false)

# Run power flow
tol = 1e-6
max_ite = 10
verbose = 0

# Run power flow
ite, erg = runpf!(net, max_ite, tol, verbose)

# Check results and calculate losses
if erg == 0
  calcNetLosses!(net)
  printACPFlowResults(net, 0.0, ite, tol)
else
  @warn "Power flow did not converge after $ite iterations"
end
```

## Working with Custom File Paths

You can specify custom paths for import and export operations:

```julia
using Sparlectra

# Import from a specific path
path = "C:/Users/YourUsername/Documents"
file = "custom_case.m"
net = run_acpflow(
    casefile = file,
    path = path,
    max_ite = 10,
    tol = 1e-6
)

# Export to a specific path
output_path = joinpath(path, "exported_case.m")
writeMatpowerCasefile(net, output_path)
```

## Matpower File Format Notes

- Sparlectra currently reads bus, generator, and branch data from Matpower files
- Additional Matlab functions within the .m file are not supported
- The Matpower format version supported is version 2

!!! Note 
    Please note that providing support for individualized network data issues is beyond the scope of this project, as it is not maintained by an organization. Users are encouraged to take initiative in resolving such issues independently and sharing their results with the community.

## API Reference

```@autodocs 
  Modules = [Sparlectra]   
  Pages = ["import.jl", "createnet_powermat.jl", "exportMatPower.jl", "run_acpflow.jl"]
  Order = [:type, :function]
```
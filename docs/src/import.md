Import and Export Netowrks Module
=============
#### Importing a Matpower Network Configuration

This example demonstrates how to import a network configuration from a Matpower case file and run a power flow analysis on it. The casefile is located in the `data` directory of the package. The `run_acpflow` function is used to run the power flow analysis.

##### Example
```julia
using Sparlectra
using Logging

file = "case3.m"
tol = 1e-8
ite = 10
verbose = 0 # 0: no output, 1: iteration norm, 2: + Y-Bus, 3: + Jacobian, 4: + Power Flow

# Call acpflow function with input parameters and measure execution time
run_acpflow(max_ite= ite,tol = tol, casefile=file)
```

!!! Note 
    Please note that providing support for individualized network data issues is beyond the scope of this project, as it is not maintained by an organization. Users are encouraged to take initiative in resolving such issues independently and sharing their results with the community.



```@autodocs 
  Modules = [Sparlectra]   
  Pages = ["import.jl","createnet_powermat.jl", "exportMatPower.jl","run_acpflow.jl"]
  Order = [:type, :function]
```  


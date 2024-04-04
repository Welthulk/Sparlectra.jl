Sparlectra
=============

Sparlectra is a Julia package for the simulation of electrical power systems. It primarily features a program for calculating load flow using the Newton-Raphson method. The focus is to provide valuable insights into load flow calculations for both students and ambitious professionals. The package supports the import and export of Matpower .m files, although currently it only reads bus, generator, and branch data from these files. Please note that additional Matlab functions within the .m file are not supported. Additionally, you can create your own network using easy-to-use functions provided by the package.

---

#### Installation
For installation, run the following command in the Julia REPL:
```julia
import Pkg
Pkg.add("Sparlectra")
```
---
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

---

#### Creating a Network
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
---

#### Exporting a Matpower Network Configuration
Once you have created a network configuration, you can export it to a Matpower case file using the following function:
```julia
  writeMatpowerCasefile(myNet, filename)
```
---

#### Contributors
- [Udo Schmitz](https://www.linkedin.com/in/udo-schmitz-315536250/) - Electrical Engineer
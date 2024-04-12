Workshop
=============
This workshop will guide you through the basic steps of using Sparlectra.jl to create and manipulate power system networks.
## Loading data from a file

```julia
using Sparlectra
using Logging
global_logger(ConsoleLogger(stderr, Logging.Warn))  # Set global logger to log to stderr with WARN level

file = "caseXYZ.m"
path="C:/Users/YourUsername/Documents"

printResultToFile = false
tol = 1e-6
ite = 10
verbose = 0       # 0: no output, 1: iteration norm, 2: + Y-Bus, 3: + Jacobian, 4: + Power Flow

net = run_acpflow(max_ite= ite,tol = tol, path=path, casefile=file, verbose=verbose, printResultToFile = printResultToFile)
```

## Creating a network from scratch and exporting it to a file

```julia
using Sparlectra
using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))  # Set global logger to log to stderr with WARN level
tol = 1e-8
ite = 10
verbose = 0       # 0: no output, 1: iteration norm, 2: + Y-Bus, 3: + Jacobian, 4: + Power Flow
writeCase = true
print_results = true

net = Net(name = "workshop_case5", baseMVA = 100.0)
addBus!(net = net, busName = "B1", busType = "PQ", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
addBus!(net = net, busName = "B2", busType = "PQ", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
addBus!(net = net, busName = "B3", busType = "PQ", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
addBus!(net = net, busName = "B4", busType = "PQ", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
addBus!(net = net, busName = "B5", busType = "Slack", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)

addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 25.0, r = 0.2, x = 0.39)
addACLine!(net = net, fromBus = "B1", toBus = "B3", length = 25.0, r = 0.2, x = 0.39)
addACLine!(net = net, fromBus = "B2", toBus = "B4", length = 25.0, r = 0.2, x = 0.39)
addACLine!(net = net, fromBus = "B3", toBus = "B4", length = 25.0, r = 0.2, x = 0.39)
addACLine!(net = net, fromBus = "B4", toBus = "B5", length = 25.0, r = 0.2, x = 0.39)

addProsumer!(net = net, busName = "B1", type = "ENERGYCONSUMER", p = 1.0, q = 2.0)
addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 1.0, q = 2.0)
addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 1.0, q = 2.0)

addProsumer!(net = net, busName = "B5", type = "SYNCHRONMASCHINE", referencePri = "B5", vm_pu = 1.0, va_deg = 0.0)
addProsumer!(net = net, busName = "B1", type = "GENERATOR", p = 1.1, q = 2.0)

result, msg = validate(net = net)
if !result
  @warn msg
  return false
end

if writeCase
  filename = "workshop_case5"
  path="C:/Users/YourUsername/Documents"  
  writeMatpowerCasefile(net, path)
end

# Run power flow
tol = 1e-6
maxIte = 10
etime = @elapsed begin
  ite, erg = runpf!(net, maxIte, tol, verbose)
end

if erg != 0
  @warn "Power flow did not converge"  
elseif print_results
  calcNetLosses!(net)
  printACPFlowResults(net, etime, ite, tol)
end
```
## Loading a file and manipulating the network

```julia
using Sparlectra

file = "case5.m" # load a case file from the data folder
net = run_acpflow(casefile=file) # run the power flow calculation

setBranchStatus!(net=net, fromBus="1", toBus="2", status=0) # set the status of the branch between Bus1 and 
run_net_acpflow(net = net) # run the power flow calculation again
addBusShuntPower!(net = net, busName = "1", p = 0.0, q = 1.0) # Update the power of Bus1 to 0.0 MW and 1.0 MVar
run_net_acpflow(net = net) # run the power flow calculation again

filename="_case5a.m" # define the filename, do not overwrite the original file and change the case-name inside the file after writing
jpath = joinpath(pwd(), "data", "mpower", filename)  
writeMatpowerCasefile(net, jpath) # write the net to a matpower case file

```
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

result, msg = validate!(net = net)
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

brVec = getNetBranchNumberVec(net = net, fromBus = "1", toBus = "2")  
setNetBranchStatus!(net = net, branchNr = brVec[1], status = 0) # set the status of the branch between Bus1 and Bus2 to 0
run_net_acpflow(net = net) # run the power flow calculation again
addBusShuntPower!(net = net, busName = "1", p = 0.0, q = 1.0) # Update the power of Bus1 to 0.0 MW and 1.0 MVar
run_net_acpflow(net = net) # run the power flow calculation again

filename="_case5a.m" # define the filename, do not overwrite the original file and change the case-name inside the file after writing
jpath = joinpath(pwd(), "data", "mpower", filename)  
writeMatpowerCasefile(net, jpath) # write the net to a matpower case file

```
## How to enable rectangular NR with Q-limits in your script
### 1. Prepare / load a network
You can build the workshop network from the earlier example and
reuse the resulting `net`.

### 2. Define PV buses and Q-limits

Make sure your generator buses are of type `PV` and set
reactive power limits (in MVar)

```julia
# Example: set bus types if needed
setBusType!(net, "B5", "Slack")   # Slack bus
setBusType!(net, "B1", "PV")      # PV bus with generator
...
addProsumer!(net = net, busName = "B1", type = "SYNCHRONOUSMACHINE", p = 10.0, q = 10.0, vm_pu = 1.03, va_deg = 0.0, qMax = 50, qMin = -50)
```

### 3. Validate the network

```julia
result, msg = validate!(net = net)
if !result
    @error "Network validation failed: $msg"
    return
end
#lockNet!(net)  # Prevent further structural modifications
```

### 4. Run rectangular NR with Q-limits

Use `runpf!` with `method = :rectangular` to select the complex-state solver.
You can also switch between analytic and finite-difference Jacobian and control
sparsity and damping.

```julia
maxIte   = 20
tol      = 1e-8
verbose  = 1      # 0: silent, 1: NR norm, 2+: more detail
damp     = 0.2    # damping for Newton steps
opt_fd   = false  # use finite-difference Jacobian if true
opt_sparse = true # use sparse Ybus / Jacobian where applicable

etime = @elapsed begin
    ite, status = runpf!(
        net,
        maxIte,
        tol,
        verbose;
        method     = :rectangular,
        damp       = damp,
        opt_fd     = opt_fd,
        opt_sparse = opt_sparse,
    )
end

if status != 0
    @warn "Rectangular NR did not converge (status = $status)"
    return
end
```

Internally this calls `runpf_rectangular!` / `run_complex_nr_rectangular_for_net!`,
which:

* uses `mismatch_rectangular` for PQ/PV constraints,
* enforces Q-limits with an active set (PV→PQ if limits are hit, optional PQ→PV),
* writes final voltages and bus injections back into `net`.

### 5. Distribute bus results to prosumers (loads / generators)

After a successful run, node powers have been computed. To see individual
loads and generators, distribute the bus results to all attached `ProSumer`s:

```julia
# Map solved bus results back down to all prosumers
distribute_all_bus_results!(net)
```

If multiple generators are connected to the same bus and one of them is at
its Q-limit, Sparlectra uses a simple “water-filling” style redistribution
so that:

* total bus P/Q stays consistent with the solved power flow,
* individual generator Q stays within its limits,
* remaining reactive power is redistributed among non-limited generators at
  that bus.

### 6. Print results (branches, buses, prosumers, Q-limit log)

You can now print the usual branch results plus the new prosumer summary and
Q-limit events.

```julia
# Branch flows and losses
calcNetLosses!(net)
printACPFlowResults(net, etime, ite, tol)

# Prosumers: list loads and generators per bus (after water-filling)
printProsumerResults(net)

# Optional: list all buses that hit Q-limits during the NR iterations
printQLimitLog(net)
```

This gives you a compact workflow:

1. Load / build the network
2. Set PV buses and Q-limits
3. Run rectangular NR via `runpf!(…; method = :rectangular, …)`
4. Call `distribute_all_bus_results!(net)`
5. Use `printACPFlowResults`, `printProsumerResults`, and `printQLimitLog`
   to inspect the outcome.

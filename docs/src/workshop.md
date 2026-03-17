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

## Working with links (bus couplers / sectionalizers)

Links are impedance-less topological connections between two buses. They are
useful to model busbar couplers or sectionalizers.

Important behavior:

* Links are **not** electrical branches and are therefore not stamped into Y-Bus.
* During `runpf!`, active links (`status = 1`) are treated as ideal couplers.
* Linked buses share the same solved voltage state when the link is closed.
* Link flows are reported after the solve via KCL-based allocation.

```julia
using Sparlectra

net = Net(name = "workshop_links", baseMVA = 100.0)

# Two buses that can be coupled by a link
addBus!(net = net, busName = "Bus1",  busType = "PQ",    vn_kV = 110.0)
addBus!(net = net, busName = "Bus1a", busType = "PQ",    vn_kV = 110.0)
addBus!(net = net, busName = "Bus4",  busType = "PQ",    vn_kV = 110.0)
addBus!(net = net, busName = "Bus5",  busType = "Slack", vn_kV = 110.0)

# Regular branches
addPIModelACLine!(net = net, fromBus = "Bus1",  toBus = "Bus4", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "Bus1a", toBus = "Bus4", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "Bus4",  toBus = "Bus5", r_pu = 0.006, x_pu = 0.050, b_pu = 0.0)

# Add a bus link (closed)
linkNr = addLink!(net = net, fromBus = "Bus1", toBus = "Bus1a", status = 1)

# Injections / loads
addProsumer!(net = net, busName = "Bus1", type = "GENERATOR", p = 45.0, q = 0.0, vm_pu = 1.01)
addProsumer!(net = net, busName = "Bus5", type = "EXTERNALNETWORKINJECTION", referencePri = "Bus5", vm_pu = 1.02, va_deg = 0.0)
addProsumer!(net = net, busName = "Bus1a", type = "LOAD", p = 30.0, q = 10.0)

ite, status, etime = run_net_acpflow(
    net = net,
    max_ite = 25,
    tol = 1e-8,
    method = :polar_full,
    opt_sparse = true,
    show_results = false,
)

# Toggle link state and rerun
setNetLinkStatus!(net = net, linkNr = linkNr, status = 0)  # open link
ite2, status2, etime2 = run_net_acpflow(
    net = net,
    max_ite = 25,
    tol = 1e-8,
    method = :polar_full,
    opt_sparse = true,
    show_results = false,
)

# Build machine-readable report object (new reporting workflow)
report = buildACPFlowReport(
    net;
    ct = etime2,
    ite = ite2,
    tol = 1e-8,
    converged = (status2 == 0),
    solver = :polar_full,
)

println(report)
println("Link rows in report: ", length(report.links))

# Optional classic text report
printACPFlowResults(net, etime2, ite2, 1e-8)
```

Notes:

* Do not connect links to a slack bus.
* Linked buses should use the same bus type (e.g. both PQ).
* Use links for topology switching logic, not for physical impedance modeling.

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
distributeBusResults!(net)
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
4. Call `distributeBusResults!(net)`
5. Use `printACPFlowResults`, `printProsumerResults`, and `printQLimitLog`
   to inspect the outcome.


## State Estimation (SE) in Sparlectra

Short method overview:

* Sparlectra uses a classical **Weighted Least Squares (WLS)** state estimator.
* The estimated state is `x = [θ(non-slack); Vm(all buses)]`, i.e. bus voltage
  angles (except slack) and voltage magnitudes.
* Measurements `z` (voltage, injections, branch flows) are linked to the
  network model through `h(x)`.
* The algorithm iteratively minimizes the weighted residual objective:
  `J(x) = (z - h(x))' * W * (z - h(x))`.

### Supported measurement types

The current SE workflow supports:

* `VmMeas` (bus voltage magnitude)
* `PinjMeas`, `QinjMeas` (bus active/reactive injection)
* `PflowMeas`, `QflowMeas` (branch active/reactive flow, direction-aware)

In most studies, measurements are generated from a converged power-flow result
with `generateMeasurementsFromPF`, optionally with Gaussian noise and
measurement-specific standard deviations.

Important: this power-flow step is only required when you want synthetic
measurements (teaching, validation, benchmarking, Monte-Carlo studies).
In real operation, SE starts from actual SCADA/PMU measurements.

### Observability

Sparlectra provides two levels:

* **Global observability** for the complete state (`evaluate_global_observability`)
* **Local observability** for selected state columns (`evaluate_local_observability`)

Useful indicators include:

* Jacobian rank / reduced Jacobian rank
* Redundancy `r = m - n` and ratio `ρ = m / n`
* Quality labels such as `:observable`, `:critical`, `:not_observable`

### Integration into network workflows

SE is integrated into the same `Net` object used by power flow:

1. Build or import a network
2. Define measurements (real telemetry or synthetic data)
3. Optional for studies only: run power flow to generate synthetic measurements
4. Check observability
5. Run `runse!`
6. Optionally write estimated states back to `net` (`updateNet = true`)

Conceptual mapping:

* Power flow computes system states from setpoints and model assumptions.
* State estimation computes system states from measured values and the model.
* With redundant measurements, SE can statistically suppress bad/noisy data and
  provide a more robust state than any single measurement channel.

## State-estimation example: simple 7-bus network

```julia
using Sparlectra
using Random

# 1) Build a simple 7-bus network
net = Net(name = "workshop_se_7bus", baseMVA = 100.0)

addBus!(net = net, busName = "B1", busType = "Slack", vn_kV = 110.0, vm_pu = 1.02, va_deg = 0.0)
for i in 2:7
    addBus!(net = net, busName = "B$(i)", busType = "PQ", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
end

# Ring + cross-connections
addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B4", toBus = "B5", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B5", toBus = "B6", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B6", toBus = "B7", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B7", toBus = "B1", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B2", toBus = "B5", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B3", toBus = "B6", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0)

# Source / generation / loads
addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)
addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 60.0, q = 10.0)
addProsumer!(net = net, busName = "B2", type = "LOAD", p = 35.0, q = 10.0)
addProsumer!(net = net, busName = "B4", type = "LOAD", p = 45.0, q = 15.0)
addProsumer!(net = net, busName = "B5", type = "LOAD", p = 25.0, q = 8.0)
addProsumer!(net = net, busName = "B6", type = "LOAD", p = 30.0, q = 10.0)
addProsumer!(net = net, busName = "B7", type = "LOAD", p = 20.0, q = 6.0)

ok, msg = validate!(net = net)
ok || error("Validation failed: $msg")

# 2) Solve reference power flow
ite_pf, status_pf = runpf!(net, 40, 1e-10, 0; method = :polar_full, opt_sparse = true)
status_pf == 0 || error("Power flow did not converge")

# 3) Build synthetic measurements (with light noise)
std = measurementStdDevs(vm = 1e-3, pinj = 0.8, qinj = 0.8, pflow = 0.5, qflow = 0.5)
meas = generateMeasurementsFromPF(
    net;
    includeVm = true,
    includePinj = true,
    includeQinj = true,
    includePflow = true,
    includeQflow = true,
    noise = true,
    stddev = std,
    rng = MersenneTwister(7),
)

# 4) Observability checks
gobs = evaluate_global_observability(net, meas; flatstart = true, jacEps = 1e-6)
println("Global quality: ", gobs.quality, ", redundancy: ", gobs.redundancy)

state_map = buildSEStateMap(net)
θ4 = state_map.theta_col_by_bus[4]
v4 = state_map.vm_col_by_bus[4]
lobs = evaluate_local_observability(net, meas, [θ4, v4]; flatstart = true, jacEps = 1e-6)
println("Local quality@bus4: ", lobs.quality, ", redundancy: ", lobs.redundancy)

# 5) Run state estimation
se = runse!(net, meas; maxIte = 12, tol = 1e-6, flatstart = true, jacEps = 1e-6, updateNet = true)
println("SE converged: ", se.converged, ", iterations: ", se.iterations, ", J: ", se.objectiveJ)
```

Workshop tips:

* Start with `noise = false`, then increase realism gradually.
* Deactivate selected measurements (`active = false`) and re-check
  global/local observability.
* Pay attention to `quality = :critical` for fragile configurations.

## State-estimation example: network + measurement set (no PF pre-step)

This variant reflects a real measurement-driven workflow: you have a network
model and a telemetry set, then run SE directly.

```julia
using Sparlectra

# 1) Build the same 7-bus network model
net = Net(name = "workshop_se_7bus_meas_only", baseMVA = 100.0)

addBus!(net = net, busName = "B1", busType = "Slack", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
for i in 2:7
    addBus!(net = net, busName = "B$(i)", busType = "PQ", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
end

addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B4", toBus = "B5", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B5", toBus = "B6", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B6", toBus = "B7", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "B7", toBus = "B1", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0)

ok, msg = validate!(net = net)
ok || error("Validation failed: $msg")

# 2) Example telemetry set (e.g. SCADA snapshot)
meas = Measurement[
    Measurement(typ = VmMeas,   value = 1.019, sigma = 0.002, busIdx = 1, id = "VM_B1"),
    Measurement(typ = VmMeas,   value = 0.997, sigma = 0.004, busIdx = 3, id = "VM_B3"),
    Measurement(typ = VmMeas,   value = 0.989, sigma = 0.004, busIdx = 5, id = "VM_B5"),

    Measurement(typ = PinjMeas, value = -33.0, sigma = 1.0, busIdx = 2, id = "PINJ_B2"),
    Measurement(typ = QinjMeas, value = -9.0,  sigma = 1.0, busIdx = 2, id = "QINJ_B2"),
    Measurement(typ = PinjMeas, value = -44.0, sigma = 1.0, busIdx = 4, id = "PINJ_B4"),
    Measurement(typ = QinjMeas, value = -14.0, sigma = 1.0, busIdx = 4, id = "QINJ_B4"),

    Measurement(typ = PflowMeas, value = 28.0, sigma = 0.8, branchIdx = 1, direction = :from, id = "PF_1_FROM"),
    Measurement(typ = QflowMeas, value = 8.0,  sigma = 0.8, branchIdx = 1, direction = :from, id = "QF_1_FROM"),
    Measurement(typ = PflowMeas, value = 21.0, sigma = 0.8, branchIdx = 3, direction = :from, id = "PF_3_FROM"),
    Measurement(typ = PflowMeas, value = 20.6, sigma = 0.8, branchIdx = 3, direction = :from, id = "PF_3_FROM_REDUNDANT"),
]

# 3) Observability and estimation
gobs = evaluate_global_observability(net, meas; flatstart = true, jacEps = 1e-6)
println("Global quality: ", gobs.quality, ", redundancy: ", gobs.redundancy)

se = runse!(net, meas; maxIte = 12, tol = 1e-6, flatstart = true, jacEps = 1e-6, updateNet = true)
println("SE converged: ", se.converged, ", iterations: ", se.iterations, ", J: ", se.objectiveJ)

# Optional: inspect largest normalized residuals for bad-data screening
for i in eachindex(se.normalizedResiduals)
    println("meas[", i, "] |r|/σ = ", abs(se.normalizedResiduals[i]))
end
```

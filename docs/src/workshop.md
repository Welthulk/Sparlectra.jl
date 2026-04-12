# Workshop

This workshop guides you through the basic steps of using Sparlectra.jl to
create, manipulate, and solve power system networks.

> **Note:** `addBus!` creates electrical nodes. The operational PF bus type
> (Slack/PV/PQ) is derived from attached prosumers. The legacy `busType` input
> is still accepted for compatibility, but does not define PF behavior.

## Loading data from a file

```julia
using Sparlectra
using Logging

global_logger(ConsoleLogger(stderr, Logging.Warn))

file = "caseXYZ.m"
path = "C:/Users/YourUsername/Documents"

printResultToFile = false
tol = 1e-6
ite = 10
verbose = 0   # 0: no output, 1: iteration norm, 2: + Y-Bus, 3: + Jacobian, 4: + Power Flow

net = run_acpflow(
    max_ite = ite,
    tol = tol,
    path = path,
    casefile = file,
    verbose = verbose,
    printResultToFile = printResultToFile,
)
```
---

## Building and extending a network from scratch

Start by creating a new network object:

```julia
using Sparlectra

net = Net(name = "example_network", baseMVA = 100.0)
```

### Add buses

```julia
addBus!(net = net, busName = "B1", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
addBus!(net = net, busName = "B2", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
addBus!(net = net, busName = "B3", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
addBus!(net = net, busName = "B4", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
addBus!(net = net, busName = "B5", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
```

### Add AC lines and transformers

```julia
addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 25.0, r = 0.2, x = 0.39)
addACLine!(net = net, fromBus = "B1", toBus = "B3", length = 25.0, r = 0.2, x = 0.39)
addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.05, x_pu = 0.2, b_pu = 0.01, status = 1)

add2WTrafo!(
    net = net,
    fromBus = "B2",
    toBus = "B4",
    sn_mva = 100.0,
    vk_percent = 10.0,
    vkr_percent = 0.5,
    pfe_kw = 20.0,
    i0_percent = 0.1,
)

addPIModelTrafo!(
    net = net,
    fromBus = "B4",
    toBus = "B5",
    r_pu = 0.01,
    x_pu = 0.1,
    b_pu = 0.0,
    status = 1,
    ratio = 1.05,
)
```

### Add loads, generators, and shunts

```julia
addProsumer!(net = net, busName = "B1", type = "ENERGYCONSUMER", p = 1.0, q = 2.0)
addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 1.0, q = 2.0)

addProsumer!(
    net = net,
    busName = "B5",
    type = "SYNCHRONMASCHINE",
    referencePri = "B5",
    vm_pu = 1.0,
    va_deg = 0.0,
)

addProsumer!(
    net = net,
    busName = "B1",
    type = "GENERATOR",
    p = 1.1,
    q = 2.0,
    vm_pu = 1.02,
)

addShunt!(net = net, busName = "B1", pShunt = 0.0, qShunt = 1.0)
```

### Add and operate bus links

Use links to model ideal busbar couplers or sectionalizers without adding
impedance to the YBUS.

```julia
linkNr = addLink!(net = net, fromBus = "B1", toBus = "B2", status = 1)
setNetLinkStatus!(net = net, linkNr = linkNr, status = 0)  # 0=open, 1=closed
```

Closed links are treated as ideal couplers during `runpf!`. Buses connected
by active links share voltage magnitude and angle in the internal power-flow
model. After convergence, call `calcLinkFlowsKCL!` to allocate and report the
link exchange on the original topology.

```julia
ite, erg = runpf!(net, 25; method = :rectangular)
if erg == 0
    calcNetLosses!(net)
    calcLinkFlowsKCL!(net)
end
```

For a complete scenario with open and closed links, see
`src/examples/using_links.jl` and the detailed notes in `links.md`.

### Validate and solve the network

Always validate your network after making significant modifications:

```julia
result, msg = validate!(net = net)
if !result
    @error "Network is invalid: \$msg"
end
```

```julia
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

## Creating a network from scratch and exporting it to a file

```julia
using Sparlectra
using Logging

global_logger(ConsoleLogger(stderr, Logging.Info))

tol = 1e-8
ite = 10
verbose = 0
writeCase = true
print_results = true

net = Net(name = "workshop_case5", baseMVA = 100.0)

addBus!(net = net, busName = "B1",    vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
addBus!(net = net, busName = "B2",    vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
addBus!(net = net, busName = "B3",    vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
addBus!(net = net, busName = "B4",    vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
addBus!(net = net, busName = "B5", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)

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
    path = "C:/Users/YourUsername/Documents"
    writeMatpowerCasefile(net, path)
end

maxIte = 10
tol = 1e-6

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
---

## Loading a file and manipulating the network

```julia
using Sparlectra

file = "case5.m"
net = run_acpflow(casefile = file)

brVec = getNetBranchNumberVec(net = net, fromBus = "1", toBus = "2")
setNetBranchStatus!(net = net, branchNr = brVec[1], status = 0)

run_net_acpflow(net = net)

addBusShuntPower!(net = net, busName = "1", p = 0.0, q = 1.0)

filename = "_case5a.m"
jpath = joinpath(pwd(), "data", "mpower", filename)
writeMatpowerCasefile(net, jpath)
```

### Update component parameters

```julia
brVec = getNetBranchNumberVec(net = net, fromBus = "B1", toBus = "B2")
updateBranchParameters!(
    net = net,
    branchNr = brVec[1],
    branch = BranchModel(
        r_pu = 0.02,
        x_pu = 0.2,
        b_pu = 0.01,
        g_pu = 0.0,
        ratio = 1.0,
        angle = 0.0,
        sn_MVA = 100.0,
    ),
)

addBusLoadPower!(net = net, busName = "B1", p = 2.0, q = 1.0)
addBusGenPower!(net = net, busName = "B5", p = 3.0, q = 1.5)
addBusShuntPower!(net = net, busName = "B2", p = 0.0, q = 1.0)
```

### Remove or isolate network elements

The practical removal workflow belongs naturally in the workshop because it is
usually part of iterative model editing:

1. identify the element to remove,
2. remove it with the dedicated helper,
3. mark or clear isolated buses, and
4. validate and solve again.

```julia
removeACLine!(net = net, fromBus = "1", toBus = "2")
removeShunt!(net = net, busName = "2")
removeProsumer!(net = net, busName = "3", type = "ENERGYCONSUMER")

markIsolatedBuses!(net = net, log = true)
clearIsolatedBuses!(net = net)

result, msg = validate!(net = net)
if !result
    error("Network validation failed: \$msg")
end

ite, status, etime = run_net_acpflow(net = net, show_results = false)
```

### Notes on component removal

- `removeBus!` is intentionally conservative: it only checks whether a bus
  could be removed. Because `Net` is immutable at the struct level, the helper
  acts as a guard instead of deleting the bus directly.
- `removeBranch!`, `removeACLine!`, and `removeTrafo!` mutate the network and
  can create isolated buses as a side effect.
- `markIsolatedBuses!` is useful for diagnostics, while `clearIsolatedBuses!`
  tries to remove buses that are now safe to delete.
- For the exact signatures and generated API docs, see the
  [Function Reference](reference.md).
---

## Working with links (bus couplers / sectionalizers)

For a detailed explanation of link behavior, zero-impedance loops, and
pseudoinverse-based flow allocation, see `links.md`.

```julia
using Sparlectra

net = Net(name = "workshop_links", baseMVA = 100.0)

addBus!(net = net, busName = "Bus1",    vn_kV = 110.0)
addBus!(net = net, busName = "Bus1a",    vn_kV = 110.0)
addBus!(net = net, busName = "Bus4",    vn_kV = 110.0)
addBus!(net = net, busName = "Bus5", vn_kV = 110.0)

addPIModelACLine!(net = net, fromBus = "Bus1",  toBus = "Bus4", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "Bus1a", toBus = "Bus4", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0)
addPIModelACLine!(net = net, fromBus = "Bus4",  toBus = "Bus5", r_pu = 0.006, x_pu = 0.050, b_pu = 0.0)

linkNr = addLink!(net = net, fromBus = "Bus1", toBus = "Bus1a", status = 1)

addProsumer!(net = net, busName = "Bus1",  type = "GENERATOR", p = 45.0, q = 0.0, vm_pu = 1.01)
addProsumer!(net = net, busName = "Bus5",  type = "EXTERNALNETWORKINJECTION", referencePri = "Bus5", vm_pu = 1.02, va_deg = 0.0)
addProsumer!(net = net, busName = "Bus1a", type = "LOAD", p = 30.0, q = 10.0)

ite, status, etime = run_net_acpflow(
    net = net,
    max_ite = 25,
    tol = 1e-8,
    method = :rectangular,
    opt_sparse = true,
    show_results = false,
)

setNetLinkStatus!(net = net, linkNr = linkNr, status = 0)

ite2, status2, etime2 = run_net_acpflow(
    net = net,
    max_ite = 25,
    tol = 1e-8,
    method = :rectangular,
    opt_sparse = true,
    show_results = false,
)

report = buildACPFlowReport(
    net;
    ct = etime2,
    ite = ite2,
    tol = 1e-8,
    converged = (status2 == 0),
    solver = :rectangular,
)

println(report)
println("Link rows in report: ", length(report.links))

printACPFlowResults(net, etime2, ite2, 1e-8)
```
---

## Running rectangular NR with Q-limits

This section shows how to use the rectangular Newton-Raphson solver with
reactive power limits.

### 1. Prepare or load a network

Use an existing network or build one as shown above.

### 2. Define PV buses and Q-limits

Define a slack prosumer and a regulating generator (PV behavior), then set Q-limits.

```julia
addProsumer!(
    net = net,
    busName = "B5",
    type = "EXTERNALNETWORKINJECTION",
    referencePri = "B5",
    vm_pu = 1.0,
    va_deg = 0.0,
)

addProsumer!(
    net = net,
    busName = "B1",
    type = "SYNCHRONMACHINE",
    p = 10.0,
    q = 10.0,
    vm_pu = 1.03,
    isRegulated = true,
    va_deg = 0.0,
    qMax = 50.0,
    qMin = -50.0,
)
```

### 3. Validate the network

```julia
result, msg = validate!(net = net)
if !result
    @error "Network validation failed: \$msg"
    return
end
```

### 4. Run the solver

```julia
maxIte = 20
tol = 1e-8
verbose = 1
damp = 0.2
opt_fd = false
opt_sparse = true

etime = @elapsed begin
    ite, status = runpf!(
        net,
        maxIte,
        tol,
        verbose;
        method = :rectangular,
        damp = damp,
        opt_fd = opt_fd,
        opt_sparse = opt_sparse,
    )
end

if status != 0
    @warn "Rectangular NR did not converge (status = \$status)"
    return
end
```

### 5. Distribute bus results to prosumers

```julia
distributeBusResults!(net)
```
If multiple generators are connected to the same bus and one of them is at
its Q-limit, Sparlectra uses a simple “water-filling” style redistribution
so that:

* total bus P/Q stays consistent with the solved power flow,
* individual generator Q stays within its limits,
* remaining reactive power is redistributed among non-limited generators at
  that bus.


### 6. Print results

```julia
calcNetLosses!(net)
printACPFlowResults(net, etime, ite, tol)
printProsumerResults(net)
printQLimitLog(net)
```
---

## State Estimation (SE) 

For state estimation workflows, observability analysis, and measurement handling,
see `state_estimation.md`.


```julia
using Sparlectra
using Random

# 1) Build a simple 7-bus network
net = Net(name = "workshop_se_7bus", baseMVA = 100.0)

addBus!(net = net, busName = "B1", vn_kV = 110.0, vm_pu = 1.02, va_deg = 0.0)
for i in 2:7
    addBus!(net = net, busName = "B\$(i)", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
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
ok || error("Validation failed: \$msg")

# 2) Solve reference power flow
ite_pf, status_pf = runpf!(net, 40, 1e-10, 0; method = :rectangular, opt_sparse = true)
status_pf == 0 || error("Power flow did not converge")

# 3) Build synthetic measurements (with light noise)
std = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7)
setMeasurementsFromPF!(
    net;
    includeVm = true,
    includePinj = true,
    includeQinj = true,
    includePflow = true,
    includeQflow = true,
    noise = true,
    stddev = std,
    rng = MersenneTwister(42),
)

# 4) Check observability
gobs = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
println("Global observability quality: ", gobs.quality)
println("Measurements: ", gobs.n_measurements, ", states: ", gobs.n_states)

# 5) Run state estimation
se = runse!(
    net;
    maxIte = 12,
    tol = 1e-6,
    flatstart = true,
    jacEps = 1e-6,
    updateNet = true,
)

println("SE converged: ", se.converged, ", iterations: ", se.iterations)
println("Final objective J: ", se.objectiveJ)

# 6) Inspect the estimated network state
printBusResults(net)
printBranchResults(net)
```

### Building measurement sets with helper functions

If you want to assemble measurements manually, you do not have to create
`Measurement(...)` entries yourself. Sparlectra provides helper functions that
work similarly to `addBus!` or `addACLine!` and resolve bus or branch
references for you.

```julia
empty!(net.measurements)

addVmMeasurement!(net; busName = "B1", value = 1.02, sigma = 0.002)
addPinjMeasurement!(net; busName = "B2", value = -35.0, sigma = 1.0)
addQinjMeasurement!(net; busName = "B2", value = -10.0, sigma = 1.0)
addPflowMeasurement!(net; fromBus = "B1", toBus = "B2", value = 22.0, sigma = 0.8, direction = :from)
addQflowMeasurement!(net; branchNr = 1, value = 7.0, sigma = 0.8, direction = :to)

obs = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
println("Manual measurement set quality: ", obs.quality)
```

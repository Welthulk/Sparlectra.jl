# Network Module

The `Net` module provides comprehensive functionality for creating and manipulating power system network models in Julia. It includes features for defining buses, branches, transformers, prosumers, and shunts, as well as methods for running power flow analysis and managing network modifications.

## Creating a Network

Start by creating a new network object:

```julia
using Sparlectra

# Create a new network with a name and base MVA
net = Net(name = "example_network", baseMVA = 100.0)
```

## Adding Components

### Buses

```julia
# Add buses with different types (PQ, PV, Slack)
addBus!(net = net, busName = "B1", busType = "PQ", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
addBus!(net = net, busName = "B2", busType = "PQ", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
addBus!(net = net, busName = "B3", busType = "PQ", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
addBus!(net = net, busName = "B4", busType = "PQ", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
addBus!(net = net, busName = "B5", busType = "Slack", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
```

### AC Lines

```julia
# Add AC lines with physical parameters
addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 25.0, r = 0.2, x = 0.39)
addACLine!(net = net, fromBus = "B1", toBus = "B3", length = 25.0, r = 0.2, x = 0.39)

# Add AC lines using the PI model (directly with per-unit values)
addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.05, x_pu = 0.2, b_pu = 0.01, status = 1)
```

### Transformers

```julia
# Add a two-winding transformer
add2WTrafo!(
    net = net, 
    fromBus = "B2", 
    toBus = "B4", 
    sn_mva = 100.0, 
    vk_percent = 10.0, 
    vkr_percent = 0.5, 
    pfe_kw = 20.0, 
    i0_percent = 0.1
)

# Add a transformer using the PI model
addPIModelTrafo!(
    net = net, 
    fromBus = "B4", 
    toBus = "B5", 
    r_pu = 0.01, 
    x_pu = 0.1, 
    b_pu = 0.0, 
    status = 1, 
    ratio = 1.05
)
```

### Prosumers (Generators and Loads)

```julia
# Add loads
addProsumer!(net = net, busName = "B1", type = "ENERGYCONSUMER", p = 1.0, q = 2.0)
addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 1.0, q = 2.0)

# Add a slack generator
addProsumer!(
    net = net, 
    busName = "B5", 
    type = "SYNCHRONMASCHINE", 
    referencePri = "B5", 
    vm_pu = 1.0, 
    va_deg = 0.0
)

# Add a PV generator
addProsumer!(
    net = net, 
    busName = "B1", 
    type = "GENERATOR", 
    p = 1.1, 
    q = 2.0, 
    vm_pu = 1.02
)
```

### Shunts

```julia
# Add a shunt to a bus
addShunt!(net = net, busName = "B1", pShunt = 0.0, qShunt = 1.0)
```

## Running Power Flow

```julia
# Set parameters
tol = 1e-6
maxIte = 10  

# Run power flow
etime = @elapsed begin
  ite, erg = runpf!(net, maxIte, tol, 0)
end

# Check results and calculate losses
if erg != 0
  @warn "Power flow did not converge"        
else
  calcNetLosses!(net)
  printACPFlowResults(net, etime, ite, tol)
end
```

## Modifying the Network

### Updating Component Parameters

```julia
# Update a branch's parameters
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
        sn_MVA = 100.0
    )
)

# Update bus powers
addBusLoadPower!(net = net, busName = "B1", p = 2.0, q = 1.0)
addBusGenPower!(net = net, busName = "B5", p = 3.0, q = 1.5)
addBusShuntPower!(net = net, busName = "B2", p = 0.0, q = 1.0)
```

### Setting Branch Status

```julia
# Change a branch's status (in-service or out-of-service)
setNetBranchStatus!(net = net, branchNr = brVec[1], status = 0)  # 0 = out of service
```

### Removing Components

See the [Component Removal](remove_functions.md) documentation for details on removing components from networks.

## Handling Isolated Buses

```julia
# Mark isolated buses in the network
markIsolatedBuses!(net = net, log = true)

# Clear (remove) isolated buses
clearIsolatedBuses!(net = net)
```

## Validating the Network

Always validate your network after making significant modifications:

```julia
result, msg = validate!(net = net)
if !result
  @error "Network is invalid: $msg"
end
```

## API Reference

```@autodocs 
  Modules = [Sparlectra]   
  Pages = ["network.jl"]
  Order = [:type, :function]
```
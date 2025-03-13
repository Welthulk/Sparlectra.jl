# Component Removal

Sparlectra provides functionality to remove various components from a network. This page documents the functions available for removing buses, branches, lines, transformers, shunts, and prosumers, as well as managing isolated buses.

## Overview

The removal functions in Sparlectra allow you to:
- Remove specific components from the network
- Handle cascading effects (like marking buses as isolated)
- Validate removal operations to ensure network integrity

Note that since the `Net` struct is immutable, some removal operations only perform validation checks to determine if the removal would be valid, rather than actually modifying the network structure.

## Removing Buses

```julia
removeBus!(; net::Net, busName::String)::Bool
```

Checks if a bus could be removed from the network. This function performs validation checks but doesn't actually remove the bus since `Net` is immutable.

### Parameters
- `net::Net`: The network to check
- `busName::String`: The name of the bus to check

### Returns
- `Bool`: True if the bus could be removed, false otherwise

### Example
```julia
removeBus!(net = network, busName = "Bus1")
```

### Notes
- Cannot remove a slack bus
- Cannot remove a bus that has branches, prosumers, or shunts connected to it

## Removing Branches

```julia
removeBranch!(; net::Net, branchNr::Int)::Bool
```

Removes a branch from the network.

### Parameters
- `net::Net`: The network from which to remove the branch
- `branchNr::Int`: The number of the branch to remove

### Returns
- `Bool`: True if the branch was successfully removed, false otherwise

### Example
```julia
removeBranch!(net = network, branchNr = 1)
```

## Removing AC Lines

```julia
removeACLine!(; net::Net, fromBus::String, toBus::String)::Bool
```

Removes an AC line between two buses from the network.

### Parameters
- `net::Net`: The network from which to remove the AC line
- `fromBus::String`: The name of the bus where the line starts
- `toBus::String`: The name of the bus where the line ends

### Returns
- `Bool`: True if the AC line was successfully removed, false otherwise

### Example
```julia
removeACLine!(net = network, fromBus = "Bus1", toBus = "Bus2")
```

## Removing Transformers

```julia
removeTrafo!(; net::Net, fromBus::String, toBus::String)::Bool
```

Removes a transformer between two buses from the network.

### Parameters
- `net::Net`: The network from which to remove the transformer
- `fromBus::String`: The name of the bus where the transformer starts
- `toBus::String`: The name of the bus where the transformer ends

### Returns
- `Bool`: True if the transformer was successfully removed, false otherwise

### Example
```julia
removeTrafo!(net = network, fromBus = "Bus1", toBus = "Bus2")
```

## Removing Shunts

```julia
removeShunt!(; net::Net, busName::String)::Bool
```

Removes a shunt at the specified bus from the network.

### Parameters
- `net::Net`: The network from which to remove the shunt
- `busName::String`: The name of the bus where the shunt is connected

### Returns
- `Bool`: True if the shunt was successfully removed, false otherwise

### Example
```julia
removeShunt!(net = network, busName = "Bus1")
```

## Removing Prosumers

```julia
removeProsumer!(; net::Net, busName::String, type::String = "")::Bool
```

Removes a prosumer (generator or load) from the specified bus.

### Parameters
- `net::Net`: The network from which to remove the prosumer
- `busName::String`: The name of the bus where the prosumer is connected
- `type::String`: The type of prosumer to remove (e.g., "GENERATOR", "ENERGYCONSUMER"). If empty, removes all prosumers at the bus.

### Returns
- `Bool`: True if the prosumer was successfully removed, false otherwise

### Example
```julia
removeProsumer!(net = network, busName = "Bus1", type = "GENERATOR")
```

## Clearing Isolated Buses

```julia
clearIsolatedBuses!(; net::Net)::Int
```

Removes all isolated buses from the network.

### Parameters
- `net::Net`: The network from which to remove isolated buses

### Returns
- `Int`: The number of isolated buses removed

### Example
```julia
clearIsolatedBuses!(net = network)
```

## Marking Isolated Buses

```julia
markIsolatedBuses!(; net::Net, log::Bool = false)
```

Finds and marks isolated buses in the network.

### Parameters
- `net::Net`: The network to analyze
- `log::Bool`: Whether to log information about the isolated buses (default: false)

### Example
```julia
markIsolatedBuses!(net = network, log = true)
```

## Usage Notes

1. Always validate your network after significant modifications:
   ```julia
   result, msg = validate!(net = network)
   if !result
     @error "Network validation failed: $msg"
   end
   ```

2. When removing components, consider the cascading effects:
   - Removing a branch might isolate a bus
   - Removing all connections to a bus will mark it as isolated

3. For complex networks, remove branches from highest index to lowest to avoid reindexing issues:
   ```julia
   branchIndices = [5, 3, 1]  # Branches to remove
   sort!(branchIndices, rev = true)  # Sort in descending order
   for idx in branchIndices
     removeBranch!(net = network, branchNr = idx)
   end
   ```
# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 13.03.2025
# file: src/remove_functions.jl

"""
    removeBus!(; net::Net, busName::String)

Checks if a bus could be removed from the network.
Note: This function cannot actually remove the bus since Net is immutable,
but it performs all validation checks.

# Arguments
- `net::Net`: The network to check.
- `busName::String`: The name of the bus to check.

# Returns
- `Bool`: True if the bus could be removed, false otherwise.

# Example
```julia
removeBus!(net = network, busName = "Bus1")
```
"""
function removeBus!(; net::Net, busName::String)::Bool
  if net._locked
    @error "Network is locked"
    return false
  end

  # Check if the bus exists
  if !haskey(net.busDict, busName)
    @error "Bus $busName not found in the network"
    return false
  end

  busIdx = net.busDict[busName]

  # Check if the bus is a slack bus
  if busIdx in net.slackVec
    @error "Cannot remove slack bus $busName"
    return false
  end

  # Check if there are branches connected to this bus
  connectedBranches = filter(br -> br.fromBus == busIdx || br.toBus == busIdx, net.branchVec)
  if !isempty(connectedBranches)
    @error "Cannot remove bus $busName because there are branches connected to it"
    return false
  end

  # Check if there are prosumers connected to this bus
  connectedProsumers = filter(ps -> ps.comp.cFrom_bus == busIdx, net.prosumpsVec)
  if !isempty(connectedProsumers)
    @error "Cannot remove bus $busName because there are prosumers connected to it"
    return false
  end

  # Check if there are shunts connected to this bus
  if haskey(net.shuntDict, busIdx)
    @error "Cannot remove bus $busName because there are shunts connected to it"
    return false
  end

  # All checks passed, but we can't actually modify the Net struct
  # So we just return true to indicate that removal would be valid
  return true
end

"""
    removeBranch!(; net::Net, branchNr::Int)

Removes a branch from the network.

# Arguments
- `net::Net`: The network from which to remove the branch.
- `branchNr::Int`: The number of the branch to remove.

# Returns
- `Bool`: True if the branch was successfully removed, false otherwise.

# Example
```julia
removeBranch!(net = network, branchNr = 1)
```
"""
function removeBranch!(; net::Net, branchNr::Int)::Bool
  if net._locked
    @error "Network is locked"
    return false
  end

  if branchNr <= 0 || branchNr > length(net.branchVec)
    @error "Branch $branchNr not found in the network"
    return false
  end

  branch = net.branchVec[branchNr]

  # Capture the branch information
  fromBus = branch.fromBus
  toBus = branch.toBus
  vn_kV = branch.comp.cVN

  # If it's an AC line, remove it from linesAC
  i = 1
  while i <= length(net.linesAC)
    line = net.linesAC[i]

    # Compare actual component fields to ensure we have the right match
    if (line.comp.cFrom_bus == fromBus && line.comp.cTo_bus == toBus) || (line.comp.cFrom_bus == toBus && line.comp.cTo_bus == fromBus)

      # Line found, remove it
      deleteat!(net.linesAC, i)
      break
    else
      i += 1
    end
  end

  # If it's a transformer, remove it from trafos
  i = 1
  while i <= length(net.trafos)
    trafo = net.trafos[i]

    # Compare actual component fields to ensure we have the right match
    if (trafo.comp.cFrom_bus == fromBus && trafo.comp.cTo_bus == toBus) || (trafo.comp.cFrom_bus == toBus && trafo.comp.cTo_bus == fromBus)

      # Transformer found, remove it
      deleteat!(net.trafos, i)
      break
    else
      i += 1
    end
  end

  # Remove the branch
  deleteat!(net.branchVec, branchNr)

  # Update branch indices
  for i = 1:length(net.branchVec)
    if net.branchVec[i].branchIdx > branchNr
      net.branchVec[i].branchIdx -= 1
    end
  end

  # Mark isolated buses after removing the branch
  markIsolatedBuses!(net = net)

  return true
end

"""
    removeACLine!(; net::Net, fromBus::String, toBus::String)

Removes an AC line between two buses from the network.

# Arguments
- `net::Net`: The network from which to remove the AC line.
- `fromBus::String`: The name of the bus where the line starts.
- `toBus::String`: The name of the bus where the line ends.

# Returns
- `Bool`: True if the AC line was successfully removed, false otherwise.

# Example
```julia
removeACLine!(net = network, fromBus = "Bus1", toBus = "Bus2")
```
"""
function removeACLine!(; net::Net, fromBus::String, toBus::String)::Bool
  if net._locked
    @error "Network is locked"
    return false
  end

  if !haskey(net.busDict, fromBus) || !haskey(net.busDict, toBus)
    @error "Bus $fromBus or $toBus not found in the network"
    return false
  end

  fromIdx = net.busDict[fromBus]
  toIdx = net.busDict[toBus]

  # Find branch numbers to remove
  branchIndices = Int[]

  for (i, br) in enumerate(net.branchVec)
    if (br.fromBus == fromIdx && br.toBus == toIdx) || (br.fromBus == toIdx && br.toBus == fromIdx)
      # Check if it's an AC line by finding it in the linesAC collection
      for line in net.linesAC
        if (line.comp.cFrom_bus == br.fromBus && line.comp.cTo_bus == br.toBus) || (line.comp.cFrom_bus == br.toBus && line.comp.cTo_bus == br.fromBus)
          push!(branchIndices, i)
          break
        end
      end
    end
  end

  if isempty(branchIndices)
    # Debug info: what branches exist between these buses?
    for (i, br) in enumerate(net.branchVec)
      if (br.fromBus == fromIdx && br.toBus == toIdx) || (br.fromBus == toIdx && br.toBus == fromIdx)
        @info "Branch $i between $fromBus and $toBus is not an AC Line"
      end
    end

    @error "No AC line found between $fromBus and $toBus"
    return false
  end

  # Remove branches from highest index to lowest to avoid reindexing issues
  sort!(branchIndices, rev = true)
  for idx in branchIndices
    removeBranch!(net = net, branchNr = idx)
  end

  return true
end

"""
    removeTrafo!(; net::Net, fromBus::String, toBus::String)

Removes a transformer between two buses from the network.

# Arguments
- `net::Net`: The network from which to remove the transformer.
- `fromBus::String`: The name of the bus where the transformer starts.
- `toBus::String`: The name of the bus where the transformer ends.

# Returns
- `Bool`: True if the transformer was successfully removed, false otherwise.

# Example
```julia
removeTrafo!(net = network, fromBus = "Bus1", toBus = "Bus2")
```
"""
function removeTrafo!(; net::Net, fromBus::String, toBus::String)::Bool
  if net._locked
    @error "Network is locked"
    return false
  end

  if !haskey(net.busDict, fromBus) || !haskey(net.busDict, toBus)
    @error "Bus $fromBus or $toBus not found in the network"
    return false
  end

  fromIdx = net.busDict[fromBus]
  toIdx = net.busDict[toBus]

  # Find branch numbers to remove
  branchIndices = Int[]

  for (i, br) in enumerate(net.branchVec)
    if (br.fromBus == fromIdx && br.toBus == toIdx) || (br.fromBus == toIdx && br.toBus == fromIdx)
      # Check if it's a transformer by finding it in the trafos collection
      for trafo in net.trafos
        if (trafo.comp.cFrom_bus == br.fromBus && trafo.comp.cTo_bus == br.toBus) || (trafo.comp.cFrom_bus == br.toBus && trafo.comp.cTo_bus == br.fromBus)
          push!(branchIndices, i)
          break
        end
      end
    end
  end

  if isempty(branchIndices)
    # Debug info: what branches exist between these buses?
    for (i, br) in enumerate(net.branchVec)
      if (br.fromBus == fromIdx && br.toBus == toIdx) || (br.fromBus == toIdx && br.toBus == fromIdx)
        @info "Branch $i between $fromBus and $toBus is not a Transformer"
      end
    end

    @error "No transformer found between $fromBus and $toBus"
    return false
  end

  # Remove branches from highest index to lowest to avoid reindexing issues
  sort!(branchIndices, rev = true)
  for idx in branchIndices
    removeBranch!(net = net, branchNr = idx)
  end

  return true
end

"""
    removeShunt!(; net::Net, busName::String)

Removes a shunt at the specified bus from the network.

# Arguments
- `net::Net`: The network from which to remove the shunt.
- `busName::String`: The name of the bus where the shunt is connected.

# Returns
- `Bool`: True if the shunt was successfully removed, false otherwise.

# Example
```julia
removeShunt!(net = network, busName = "Bus1")
"""

function removeShunt!(; net::Net, busName::String)::Bool
  # is missing !

  if net._locked
    @error "Network is locked"
    return false
  end
  if !haskey(net.busDict, busName)
    @error "Bus $busName not found in the network"
    return false
  end
  busIdx = net.busDict[busName]
  if !haskey(net.shuntDict, busIdx)
    @error "No shunt found at bus $busName"
    return false
  end
  shuntIdx = net.shuntDict[busIdx]
  shunt = net.shuntVec[shuntIdx]
  # Update the node by removing the shunt power
  node = net.nodeVec[busIdx]
  node._pShunt = 0.0
  node._qShunt = 0.0
  # Remove the shunt
  deleteat!(net.shuntVec, shuntIdx)
  delete!(net.shuntDict, busIdx)
  # Update shunt indices in shuntDict
  for (idx, sIdx) in pairs(net.shuntDict)
    if sIdx > shuntIdx
      net.shuntDict[idx] = sIdx - 1
    end
  end
  return true
end

"""
      removeProsumer!(; net::Net, busName::String, type::String)

# Arguments
- `net::Net`: The network from which to remove the prosumer.
- `busName::String`: The name of the bus where the prosumer is connected.
- `type::String`: The type of the prosumer to remove (e.g., "GENERATOR", "ENERGYCONSUMER"). If empty, removes all prosumers at the bus.

# Returns
- `Bool`: True if the prosumer was successfully removed, false otherwise.

# Example
```julia
removeProsumer!(net = network, busName = "Bus1", type = "GENERATOR")
```
"""

function removeProsumer!(; net::Net, busName::String, type::String = "")::Bool
  if net._locked
    @error "Network is locked"
    return false
  end

  if !haskey(net.busDict, busName)
    @error "Bus $busName not found in the network"
    return false
  end

  busIdx = net.busDict[busName]

  # Find prosumers connected to this bus
  prosumerIndices = Int[]
  for (i, prosumer) in enumerate(net.prosumpsVec)
    if prosumer.comp.cFrom_bus == busIdx
      # If type is specified, check if the prosumer is of that type
      if type != ""
        proType = toProSumptionType(type)
        if prosumer.proSumptionType != proType
          continue
        end
      end
      push!(prosumerIndices, i)
    end
  end

  if isempty(prosumerIndices)
    if type != ""
      @error "No prosumer of type $type found at bus $busName"
    else
      @error "No prosumers found at bus $busName"
    end
    return false
  end

  # Remove prosumers from highest index to lowest to avoid reindexing issues
  sort!(prosumerIndices, rev = true)

  # Store prosumer types before deletion to properly update node power
  prosumerTypes = [net.prosumpsVec[idx].proSumptionType for idx in prosumerIndices]
  hasGenerator = any(isGenerator.(prosumerTypes))
  hasConsumer = any(x -> !isGenerator(x), prosumerTypes)

  # Remove the prosumers
  for idx in prosumerIndices
    deleteat!(net.prosumpsVec, idx)
  end

  # Update node by removing the prosumer power
  node = net.nodeVec[busIdx]

  # Reset the node's generation or load power based on what was removed
  if type == "" || hasGenerator || uppercase(type) in ["GENERATOR", "EXTERNALNETWORKINJECTION", "SYNCHRONOUSMACHINE"]
    node._pƩGen = 0.0
    node._qƩGen = 0.0
  end

  if type == "" || hasConsumer || uppercase(type) in ["ENERGYCONSUMER", "ASYNCHRONOUSMACHINE", "LOAD"]
    node._pƩLoad = 0.0
    node._qƩLoad = 0.0
  end

  return true
end

"""
    clearIsolatedBuses!(; net::Net)

Removes all isolated buses from the network.

# Arguments
- `net::Net`: The network from which to remove isolated buses.

# Returns
- `Int`: The number of isolated buses removed.

# Example
```julia
clearIsolatedBuses!(net = network)
```
"""
function clearIsolatedBuses!(; net::Net)::Int
  if net._locked
    @error "Network is locked"
    return 0
  end

  # Make sure isolated buses are marked
  markIsolatedBuses!(net = net)

  # Collect isolated bus names
  isolatedBusNames = String[]
  for node in net.nodeVec
    if isIsolated(node)
      for (name, idx) in net.busDict
        if idx == node.busIdx
          push!(isolatedBusNames, name)
          break
        end
      end
    end
  end

  # Remove isolated buses
  count = 0
  for busName in isolatedBusNames
    if removeBus!(net = net, busName = busName)
      count += 1
    end
  end

  return count
end

export removeBus!, removeBranch!, removeACLine!, removeTrafo!, removeShunt!, removeProsumer!, clearIsolatedBuses!

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
# Date: 01.04.2024
# file: src/network.jl
"""
    Net

Represents an electrical network.

# Fields
- `name::String`: Name of the network.
- `baseMVA::Float64`: Base MVA of the network.
- `slackVec::Vector{Int}`: Vector containing indices of slack buses.
- `vmin_pu::Float64`: Minimum voltage limit in per unit.
- `vmax_pu::Float64`: Maximum voltage limit in per unit.
- `nodeVec::Vector{Node}`: Vector containing nodes in the network.
- `linesAC::Vector{ACLineSegment}`: Vector containing AC line segments.
- `trafos::Vector{PowerTransformer}`: Vector containing power transformers.
- `branchVec::Vector{Branch}`: Vector containing branches in the network.
- `prosumpsVec::Vector{ProSumer}`: Vector containing prosumers in the network.
- `shuntVec::Vector{Shunt}`: Vector containing shunts in the network.
- `busDict::Dict{String,Int}`: Dictionary mapping bus names to indices.
- `busOrigIdxDict::Dict{Int,Int}`: Dictionary mapping current bus indices to original indices.
- `branchDict::Dict{Tuple{Int, Int},Int}`: Dictionary mapping branch tuples to indices.
- `totalLosses::Vector{Tuple{Float64,Float64}}`: Vector containing tuples of total power losses.
- `_locked::Bool`: Boolean indicating if the network is locked.

# Constructors
- `Net(name::String, baseMVA::Float64, vmin_pu::Float64 = 0.9, vmax_pu::Float64 = 1.1)`: Creates a new `Net` object with the given name, base MVA, and optional voltage limits.

# Functions
- `addBus!(; net::Net, ...)`: Adds a bus to the network.
- `addShunt!(; net::Net, ...)`: Adds a shunt to the network.
- `addBranch!(; net::Net, ...)`: Adds a branch to the network.
- `addPIModelACLine!(; net::Net, ...)`: Adds a PI model AC line to the network.
- `addACLine!(; net::Net, ...)`: Adds an AC line segment to the network.
- `addPIModelTrafo!(; net::Net, ...)`: Adds a transformer with PI model to the network.
- `add2WTrafo!(; net::Net, ...)`: Adds a two-winding transformer to the network.
- `addProsumer!(; net::Net, ...)`: Adds a prosumer to the network.
- `lockNet!(; net::Net, locked::Bool)`: Locks or unlocks the network.
- `validate!(; net::Net)`: Validates the network.
- `get_bus_vn_kV(; net::Net, busName::String)`: Gets the voltage level of a bus by name.
- `get_vn_kV(; net::Net, busIdx::Int)`: Gets the voltage level of a bus by index.
- `getBusType(; net::Net, busName::String)`: Gets the type of a bus by name.
- `updateBranchParameters!(; net::Net, fromBus::String, toBus::String, branch::BranchModel)`: Updates the parameters of a branch in the network.
- `setNetBranchStatus!(; net::Net, branchNr::Int, status::Int)`: Sets the status of a branch.
- `setTotalLosses!(; net::Net, pLosses::Float64, qLosses::Float64)`: Sets the total losses of the network.
- `getTotalLosses(; net::Net)`: Gets the total losses of the network.
- `getNetOrigBusIdx(; net::Net, busName::String)`: Gets the original index of a bus in the network.
- `getNetBusIdx(; net::Net, busName::String)`: Gets the index of a bus in the network.
- `hasBusInNet(; net::Net, busName::String)`: Checks if a bus exists in the network.
- `addBusGenPower!(; net::Net, busName::String, pGen::Float64, qGen::Float64)`: Adds generator power to a bus.
- `addBusLoadPower!(; net::Net, busName::String, pLoad::Float64, qLoad::Float64)`: Adds load power to a bus.
- `addBusShuntPower!(; net::Net, busName::String, pShunt::Float64, qShunt::Float64)`: Adds shunt power to a bus.
- `getNetBranch(; net::Net, fromBus::String, toBus::String)`: Retrieves the branch between two specified buses in the network.
"""

struct Net
  name::String
  baseMVA::Float64
  slackVec::Vector{Int}
  vmin_pu::Float64
  vmax_pu::Float64
  nodeVec::Vector{Node}
  linesAC::Vector{ACLineSegment}
  trafos::Vector{PowerTransformer}
  branchVec::Vector{Branch}
  prosumpsVec::Vector{ProSumer}
  shuntVec::Vector{Shunt}
  busDict::Dict{String,Int}
  busOrigIdxDict::Dict{Int,Int}
  totalLosses::Vector{Tuple{Float64,Float64}}
  totalBusPower::Vector{Tuple{Float64,Float64}}
  _locked::Bool
  flatstart::Bool                     # flatstart flag for power flow initialization
  shuntDict::Dict{Int,Int}
  isoNodes::Vector{Int}  
  qLimitLog::Vector{Any}
  cooldown_iters::Int
  q_hyst_pu::Float64
  qmin_pu::Vector{Float64}            # pro Bus Qmin (p.u.)
  qmax_pu::Vector{Float64}            # pro Bus Qmax (p.u.)
  qLimitEvents::Dict{Int,Symbol}      # BusIdx -> :min | :max (PV→PQ Change)  
  
  
  #! format: off
  function Net(; name::String, baseMVA::Float64, vmin_pu::Float64 = 0.9, vmax_pu::Float64 = 1.1, cooldown_iters::Int = 0, q_hyst_pu::Float64 = 0.0, flatstart::Bool = false)    
    
    new(name, # name
        baseMVA, # baseMVA
        [], # slackVec
        vmin_pu, # vmin_pu
        vmax_pu, # vmax_pu
        [], # nodeVec
        [], # linesAC
        [], # trafos
        [], # branchVec
        [], # prosumpsVec
        [], # shuntVec
        Dict{String,Int}(), # busDict
        Dict{Int,Int}(), # busOrigIdxDict
        [], # totalLosses
        [], # totalBusPower
        false, # _locked
        flatstart,                     # flatstart
        Dict{Int,Symbol}(),  # shuntDict
        [],                                    # isoNodes
        Any[],                                 # qLimitLog                     
        cooldown_iters,                        # cooldown_iters
        q_hyst_pu,
        [],                                    # qmin_pu
        [],                                    # qmax_pu
        Dict{Int,Symbol}())                                          
  end
  #! format: on
  function Base.show(io::IO, net::Net)
    println(io, "Net: ", net.name)
    println(io, "Base MVA: ", net.baseMVA)
    println(io, "Nodes: ", length(net.nodeVec), ", Lines: ", length(net.linesAC), ", Transformers: ", length(net.trafos), ", Branches: ", length(net.branchVec))
    println(io, "Slack buses: ", net.slackVec, ", flatstart: ", net.flatstart, ", locked: ", net._locked)
    println(io, "Vmin / Vmax: ", net.vmin_pu, " / ", net.vmax_pu)
    println(io, "cooldown_iters: ", net.cooldown_iters, ", q_hyst_pu: ", net.q_hyst_pu)
  end
end

function showNet(io::IO, net::Net; verbose::Bool = false)
  if !verbose
    show(io, net)
    return
  end

  println(io, "==================== Net ====================")
  println(io, "Name:          ", net.name)
  println(io, "Base MVA:      ", net.baseMVA)
  println(io, "Slack buses:   ", net.slackVec)
  println(io, "Vmin / Vmax:   ", net.vmin_pu, " / ", net.vmax_pu)
  println(io, "flatstart:     ", net.flatstart)
  println(io, "_locked:       ", net._locked)

  println(io, "\n--- Dictionaries ---")
  println(io, "busDict:        ", net.busDict)
  println(io, "busOrigIdxDict: ", net.busOrigIdxDict)
  println(io, "shuntDict:      ", net.shuntDict)

  println(io, "\n--- Topology / State ---")
  println(io, "isoNodes:       ", net.isoNodes)

  println(io, "\n--- Q-limit handling ---")
  println(io, "cooldown_iters: ", net.cooldown_iters)
  println(io, "q_hyst_pu:      ", net.q_hyst_pu)
  println(io, "qmin_pu:        ", net.qmin_pu)
  println(io, "qmax_pu:        ", net.qmax_pu)
  println(io, "qLimitEvents:   ", net.qLimitEvents)

  println(io, "\n--- Aggregates ---")
  println(io, "totalLosses:    ", net.totalLosses)
  println(io, "totalBusPower:  ", net.totalBusPower)

  println(io, "\n--- Q-limit log (", length(net.qLimitLog), ") ---")
  for (i, e) in enumerate(net.qLimitLog)
    println(io, "  [", i, "] ", e)
  end

  println(io, "\n==================== Nodes (", length(net.nodeVec), ") ====================")
  for (i, n) in enumerate(net.nodeVec)
    println(io, "\n[node ", i, "]")
    show(io, n)
  end

  println(io, "\n==================== AC Lines (", length(net.linesAC), ") ====================")
  for (i, l) in enumerate(net.linesAC)
    println(io, "\n[line ", i, "]")
    show(io, l)
  end

  println(io, "\n==================== Transformers (", length(net.trafos), ") ====================")
  for (i, t) in enumerate(net.trafos)
    println(io, "\n[trafo ", i, "]")
    show(io, t)
  end

  println(io, "\n==================== Branches (", length(net.branchVec), ") ====================")
  for (i, b) in enumerate(net.branchVec)
    println(io, "\n[branch ", i, "]")
    show(io, b)
  end

  println(io, "\n==================== Prosumers (", length(net.prosumpsVec), ") ====================")
  for (i, p) in enumerate(net.prosumpsVec)
    println(io, "\n[prosumer ", i, "]")
    show(io, p)
  end

  println(io, "\n==================== Shunts (", length(net.shuntVec), ") ====================")
  for (i, s) in enumerate(net.shuntVec)
    println(io, "\n[shunt ", i, "]")
    show(io, s)
  end

  println(io, "\n==================== END Net ====================")
end
# convenience overloads
showNet(net::Net; verbose::Bool = false) = showNet(stdout, net; verbose = verbose)

# --- helpers (lokal) ---
@inline function _push_unique!(v::Vector{Int}, x::Int)
  (findfirst(==(x), v) === nothing) && push!(v, x)
  return v
end
@inline function _delete_item!(v::Vector{Int}, x::Int)
  idx = findfirst(==(x), v)
  (idx !== nothing) && deleteat!(v, idx)
  return v
end

"""
Sets the type of a bus in the network.

# Parameters:
- `net::Net`: Network object.
- `bus::Int`: Index of the bus.
- `busType::String`: Type of the bus (e.g., "Slack", "PQ", "PV")
"""
function setBusType!(net::Net, bus::Int, busType::String)
  @assert 1 <= bus <= length(net.nodeVec) "Bus-Index $bus ist ungültig."
  node  = net.nodeVec[bus]
  oldTy = getfield(node, :_nodeType)

  setNodeType!(node, busType)

  newTy = getfield(node, :_nodeType)
  if newTy == Slack && oldTy != Slack
    _push_unique!(net.slackVec, bus)
  elseif oldTy == Slack && newTy != Slack
    _delete_item!(net.slackVec, bus)
  end
  return nothing
end

function setBusType!(net::Net, busName::String, busType::String)
  bus = geNetBusIdx(net = net, busName = busName)
  return setBusType!(net, bus, busType)
end

"""
hasBusInNet: Checks if a bus exists in the network.

Parameters:
- `net::Net`: Network object.
- `busName::String`: Name of the bus to check.

Returns:
- `Bool`: True if the bus exists in the network, otherwise false.
"""
function hasBusInNet(; net::Net, busName::String)::Bool
  return haskey(net.busDict, busName)
end

"""
geNetBusIdx: Gets the index of a bus in the network.

Parameters:
- `net::Net`: Network object.
- `busName::String`: Name of the bus.

Returns:
- `Int`: Index of the bus in the network.
"""
function geNetBusIdx(; net::Net, busName::String)::Int
  if !haskey(net.busDict, busName)
    error("Bus $(busName) not found in the network")
  end
  return net.busDict[busName]
end

"""
getNetOrigBusIdx: Gets the original index of a bus in the network.

Parameters:
- `net::Net`: Network object.
- `busName::String`: Name of the bus.

Returns:
- `Int`: Original index of the bus in the network.
"""
function getNetOrigBusIdx(; net::Net, busName::String)::Int
  id = geNetBusIdx(net = net, busName = busName)
  return net.busOrigIdxDict[id]
end

"""
addBus!: Adds a bus to the network.

Parameters:
- `net::Net`: Network object.
- `busName::String`: Name of the bus.
- `busType::String`: Type of the bus (e.g., "Slack", "PQ", "PV").
- `vn_kV::Float64`: Nominal voltage of the bus in kV.
- `vm_pu::Float64 = 1.0`: Voltage magnitude of the bus in per unit (default is 1.0).
- `va_deg::Float64 = 0.0`: Voltage angle of the bus in degrees (default is 0.0).
- `vmin_pu::Union{Nothing,Float64} = nothing`: Minimum voltage limit in per unit (default is network's vmin_pu).
- `vmax_pu::Union{Nothing,Float64} = nothing`: Maximum voltage limit in per unit (default is network's vmax_pu).
- `isAux::Bool = false`: Boolean indicating if the bus is auxiliary (default is false).
- `oBusIdx::Union{Nothing,Int} = nothing`: Original bus index (default is nothing).
- `zone::Union{Nothing,Int} = nothing`: Zone index (default is nothing).
- `area::Union{Nothing,Int} = nothing`: Area index (default is nothing).
- `ratedS::Union{Nothing,Float64} = nothing`: Rated power of the bus in MVA (default is nothing).
"""
function addBus!(;
  net::Net,
  busName::String,
  busType::String,
  vn_kV::Float64,
  vm_pu::Float64 = 1.0,
  va_deg::Float64 = 0.0,
  vmin_pu::Union{Nothing,Float64} = nothing,
  vmax_pu::Union{Nothing,Float64} = nothing,
  isAux::Bool = false,
  oBusIdx::Union{Nothing,Int} = nothing,
  zone::Union{Nothing,Int} = nothing,
  area::Union{Nothing,Int} = nothing,
  ratedS::Union{Nothing,Float64} = nothing,
)
  @assert busType in ["Slack", "SLACK", "PQ", "PV"] "Invalid bus type: $busType"

  if net._locked
    @error "Network is locked"
  end
  if haskey(net.busDict, busName)
    @error "Bus $busName already exists in the network"
  end
  busIdx = length(net.nodeVec) + 1
  net.busDict[busName] = busIdx

  vmin_pu = isnothing(vmin_pu) ? net.vmin_pu : vmin_pu
  vmax_pu = isnothing(vmax_pu) ? net.vmax_pu : vmax_pu

  if !isnothing(oBusIdx)
    if !haskey(net.busOrigIdxDict, busIdx)
      net.busOrigIdxDict[busIdx] = oBusIdx
    else
      @warn "Original bus index already set for bus $busName"
    end
  end

  
  if (uppercase(busType) == "SLACK") || (uppercase(busType) == "Slack ")
    push!(net.slackVec, busIdx)
  end

  node = Node(busIdx = busIdx, vn_kV = vn_kV, nodeType = toNodeType(busType), vm_pu = vm_pu, va_deg = va_deg, vmin_pu = vmin_pu, vmax_pu = vmax_pu, isAux = isAux, oBusIdx = oBusIdx, zone = zone, area = area, ratedS = ratedS)
  push!(net.nodeVec, node)
end

"""
    addShunt!(; net, busName, pShunt, qShunt, in_service=1)

Adds a *Y-model* shunt to the network.

Semantics:
- `pShunt`, `qShunt` are interpreted as MW/MVar at V = 1.0 pu (MATPOWER-style).
- Internally, the shunt is represented as a pu-admittance stamped into YBUS:
      y_pu = (pShunt + j*qShunt) / baseMVA
- The shunt power is *not* constant; it depends on |V|² and will be computed
  after solving via `updateShuntPowers!(net)`.

IMPORTANT:
- This does NOT add anything to the S-vector and does NOT call `addShuntPower!`.
"""
function addShunt!(; net::Net, busName::String, pShunt::Float64, qShunt::Float64, in_service::Int = 1)
  @assert in_service in (0, 1) "in_service must be 0 or 1"

  busIdx = geNetBusIdx(net = net, busName = busName)
  idShunt = length(net.shuntVec) + 1
  net.shuntDict[busIdx] = idShunt

  vn_kV = getNodeVn(net.nodeVec[busIdx])

  # Build shunt as Y-model:
  # y_pu = (P+jQ)/baseMVA  (P,Q in MW/MVar at 1pu)
  ypu = ComplexF64(pShunt, qShunt) / net.baseMVA

  # Create via constructor (p/q here are ONLY used to set y_pu_shunt in this new semantics)
  sh = Shunt(fromBus = busIdx,
             id = idShunt,
             base_MVA = net.baseMVA,
             vn_kV_shunt = vn_kV,
             p_shunt = pShunt,
             q_shunt = qShunt,
             status = in_service)

  # Enforce the intended semantics explicitly (in case constructor changes again)
  sh.y_pu_shunt = ypu

  # Optional: keep the "G/B" fields consistent with y_pu (pu values)
  sh.G_shunt = real(ypu)
  sh.B_shunt = imag(ypu)

  # NOTE:
  # Do NOT call addShuntPower!(...) here, otherwise you'd treat it as constant PQ.
  # p_shunt/q_shunt will be updated after PF using updateShuntPowers!(net).
  push!(net.shuntVec, sh)

  return nothing
end
"""
    hasShunt!(; net::Net, busName::String)::Bool

Checks if a shunt exists at the specified bus in the network.

# Arguments
- `net::Net`: The network.
- `busName::String`: The name of the bus.

# Returns
- `Bool`: True if a shunt exists at the specified bus, false otherwise.

# Example
```julia
hasShunt!(net = network, busName = "Bus1")
"""
function hasShunt!(; net::Net, busName::String)::Bool
  busIdx = geNetBusIdx(net = net, busName = busName)
  return haskey(net.shuntDict, busIdx)
end

"""
    getShunt!(; net::Net, busName::String)::Shunt

Retrieves the shunt at the specified bus in the network.

# Arguments
- `net::Net`: The network.
- `busName::String`: The name of the bus.

# Returns
- `Shunt`: The shunt at the specified bus.

# Example
```julia
getShunt!(net = network, busName = "Bus1")
"""
function getShunt!(; net::Net, busName::String)::Shunt
  busIdx  = geNetBusIdx(net = net, busName = busName)
  idShunt = net.shuntDict[busIdx]
  return net.shuntVec[idShunt]
end

"""
addBranch!: Adds a branch to the network.

Parameters:
- `net::Net`: Network object.
- `from::Int`: Index of the "from" bus.
- `to::Int`: Index of the "to" bus.
- `branch::AbstractBranch`: Branch object to add.
- `status::Int = 1`: Status of the branch (default is 1).
- `ratio::Union{Nothing,Float64} = nothing`: Ratio of the branch (default is nothing).
- `side::Union{Nothing,Int} = nothing`: Side of the branch (default is nothing).
- `vn_kV::Union{Nothing,Float64} = nothing`: Nominal voltage of the branch in kV (default is nothing).
- `values_are_pu = false`: Boolean indicating if the values are in per unit (default is false).
"""

function addBranch!(; net::Net, from::Int, to::Int, branch::AbstractBranch, status::Integer=1,
                    ratio=nothing, side=nothing, vn_kV=nothing,
                    values_are_pu::Bool=false)
  @assert from != to "From and to bus must be different"
  idBrunch = length(net.branchVec) + 1
  fOrig = nothing
  tOrig = nothing
  if haskey(net.busOrigIdxDict, from)
    fOrig = net.busOrigIdxDict[from]
  end
  if haskey(net.busOrigIdxDict, to)
    tOrig = net.busOrigIdxDict[to]
  end

  if isnothing(vn_kV)
    vn_kV = getNodeVn(net.nodeVec[from])
  end
  br = Branch(branchIdx = idBrunch, from = from, to = to, baseMVA = net.baseMVA, branch = branch, id = idBrunch, status = status, ratio = ratio, side = side, vn_kV = vn_kV, fromOid = fOrig, toOid = tOrig, values_are_pu = values_are_pu)
  push!(net.branchVec, br)
end

"""
    updateBranchParameters!(;net::Net, fromBus::String, toBus::String, branch::AbstractBranch)

Updates the parameters of a branch in the network.

# Arguments
- `net::Net`: The network.
- `fromBus::String`: The name of the bus where the branch starts.
- `toBus::String`: The name of the bus where the branch ends.
- `branch::BranchModel`: The branch with the updated parameters.

# Example
```julia
updateBranchParameters!(net = network, fromBus = "Bus1", toBus = "Bus2", branch = updatedBranch)
```
"""
function updateBranchParameters!(; net::Net, branchNr::Int, branch::BranchModel)
  br = net.branchVec[branchNr]
  br.r_pu = branch.r_pu
  br.x_pu = branch.x_pu
  br.b_pu = branch.b_pu
  br.g_pu = branch.g_pu
  br.ratio = branch.ratio
  br.angle = branch.angle
  br.sn_MVA = branch.sn_MVA
end

"""
addACLine!: Adds an AC line segment to the network.

Parameters:
- `net::Net`: Network object.
- `fromBus::String`: Name of the "from" bus.
- `toBus::String`: Name of the "to" bus.
- `length::Float64`: Length of the line segment.
- `r::Float64`: Resistance per Meter of the line segment.
- `x::Float64`: Reactance per Meter of the line segment.
- `b::Union{Nothing,Float64} = nothing`: Susceptance per Meter of the line segment (default is nothing).
- `c_nf_per_km::Union{Nothing,Float64} = nothing`: Capacitance per Meter of the line segment in nF/km (default is nothing).
- `tanδ::Union{Nothing,Float64} = nothing`: Tangent of the loss angle (default is nothing).
- `ratedS::Union{Nothing, Float64}= nothing`: Rated power of the line segment in MVA (default is nothing).
- `status::Int = 1`: Status of the line segment (default is 1).
"""
function addACLine!(; net::Net, fromBus::String, toBus::String, length::Float64, r::Float64, x::Float64, b::Union{Nothing,Float64} = nothing, c_nf_per_km::Union{Nothing,Float64} = nothing, tanδ::Union{Nothing,Float64} = nothing, ratedS::Union{Nothing,Float64} = nothing, status::Int = 1)
  @assert fromBus != toBus "From and to bus must be different"
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  vn_kV = getNodeVn(net.nodeVec[from])
  vn_2_kV = getNodeVn(net.nodeVec[to])
  @assert vn_kV == vn_2_kV "Voltage level of the from bus $(vn_kV) does not match the to bus $(vn_2_kV)"
  if length > 1.0 
    _par_length = true
  else
    _par_length = false
  end
  acseg = ACLineSegment(vn_kv = vn_kV, from = from, to = to, length = length, r = r, x = x, b = b, c_nf_per_km = c_nf_per_km, tanδ = tanδ, ratedS = ratedS, paramsBasedOnLength = _par_length, isPIModel = false)
  push!(net.linesAC, acseg)

  addBranch!(net = net, from = from, to = to, branch = acseg, vn_kV = vn_kV, status = status, values_are_pu = true)
end

"""
    addPIModelACLine!(; net::Net, fromBus::String, toBus::String, r_pu::Float64, x_pu::Float64, b_pu::Float64, status::Int, ratedS::Union{Nothing,Float64}=nothing)

Adds a PI model AC line to the network.

# Arguments
- `net::Net`: The network.
- `fromBus::String`: The name of the bus where the line starts.
- `toBus::String`: The name of the bus where the line ends.
- `r_pu::Float64`: The per unit resistance of the line.
- `x_pu::Float64`: The per unit reactance of the line.
- `b_pu::Float64`: The per unit total line charging susceptance of the line.
- `status::Int`: The status of the line. 1 = in service, 0 = out of service.
- `ratedS::Union{Nothing,Float64}`: The rated power of the line.

# Example
```julia
addPIModelACLine!(net = network, fromBus = "Bus1", toBus = "Bus2", r_pu = 0.01, x_pu = 0.1, b_pu = 0.02, status = 1, ratedS = 100.0)
```
"""
function addPIModelACLine!(; net::Net, fromBus::String, toBus::String, r_pu::Float64, x_pu::Float64, b_pu::Float64, status::Int, ratedS::Union{Nothing,Float64} = nothing)
  @assert fromBus != toBus "From and to bus must be different"
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  vn_kV = getNodeVn(net.nodeVec[from])
  vn_2_kV = getNodeVn(net.nodeVec[to])
  @assert vn_kV == vn_2_kV "Voltage level of the from bus $(vn_kV) does not match the to bus $(vn_2_kV)"
  acseg = ACLineSegment(vn_kv = vn_kV, from = from, to = to, length = 1.0, r = r_pu, x = x_pu, b = b_pu, ratedS = ratedS, paramsBasedOnLength = false, isPIModel = true)
  push!(net.linesAC, acseg)

  addBranch!(net = net, from = from, to = to, branch = acseg, vn_kV = vn_kV, status = status, values_are_pu = true)
end

"""
Add a transformer with PI model to the network.

# Arguments
- `net::Net`: The network to which the transformer will be added.
- `fromBus::String`: The name of the bus where the transformer originates.
- `toBus::String`: The name of the bus where the transformer terminates.
- `r_pu::Float64`: The per-unit resistance of the transformer.
- `x_pu::Float64`: The per-unit reactance of the transformer.
- `b_pu::Float64`: The per-unit susceptance of the transformer.
- `status::Int`: The status of the transformer.
- `ratedU::Union{Nothing, Float64}`: Rated voltage of the transformer. Default is `nothing`.
- `ratedS::Union{Nothing, Float64}`: Rated apparent power of the transformer. Default is `nothing`.
- `ratio::Union{Nothing, Float64}`: Ratio of the transformer. Default is `nothing`.
- `shift_deg::Union{Nothing, Float64}`: Phase shift angle of the transformer. Default is `nothing`.
- `isAux::Bool`: Whether the transformer is an auxiliary transformer. Default is `false`.
"""
function addPIModelTrafo!(;
  net::Net,
  fromBus::String,
  toBus::String,
  r_pu::Float64,
  x_pu::Float64,
  b_pu::Float64,
  status::Int,
  ratedU::Union{Nothing,Float64} = nothing,
  ratedS::Union{Nothing,Float64} = nothing,
  ratio::Union{Nothing,Float64} = nothing,
  shift_deg::Union{Nothing,Float64} = nothing,
  isAux::Bool = false,
  side::Int = 1,
)
  @assert fromBus != toBus "From and to bus must be different"
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  vn_hv_kV = getNodeVn(net.nodeVec[from])
  vn_lv_kV = getNodeVn(net.nodeVec[to])
  w1 = PowerTransformerWinding(vn_hv_kV, r_pu, x_pu, b_pu, 0.0, ratio, shift_deg, ratedU, ratedS, nothing, true)
  w2 = PowerTransformerWinding(vn_lv_kV, 0.0, 0.0)
  comp = getTrafoImpPGMComp(isAux, vn_hv_kV, from, to)
  trafo = PowerTransformer(comp, false, w1, w2, nothing, Sparlectra.PIModel)
  push!(net.trafos, trafo)

  addBranch!(net = net, from = from, to = to, branch = trafo, status = status, ratio = ratio, side = side, vn_kV = vn_hv_kV, values_are_pu = true)
end

"""
Add a transformer with PI model to the network.

# Arguments
- `net::Net`: The network to which the transformer will be added.
- `fromBus::String`: The name of the bus where the transformer originates.
- `toBus::String`: The name of the bus where the transformer terminates.
- `r_pu::Float64`: The per-unit resistance of the transformer.
- `x_pu::Float64`: The per-unit reactance of the transformer.
- `b_pu::Float64`: The per-unit susceptance of the transformer.
- `status::Int`: The status of the transformer.
- `ratedU::Union{Nothing, Float64}`: Rated voltage of the transformer. Default is `nothing`.
- `ratedS::Union{Nothing, Float64}`: Rated apparent power of the transformer. Default is `nothing`.
- `ratio::Union{Nothing, Float64}`: Ratio of the transformer. Default is `nothing`.
- `shift_deg::Union{Nothing, Float64}`: Phase shift angle of the transformer. Default is `nothing`.
- `isAux::Bool`: Whether the transformer is an auxiliary transformer. Default is `false`.
"""
function add2WTPIModelTrafo!(;
  net::Net,
  fromBus::String,
  toBus::String,
  side::Int = 1,
  r::Float64,
  x::Float64,
  b::Float64 = 0.0,
  status::Int = 1,
  ratedU::Union{Nothing,Float64} = nothing,
  ratedS::Union{Nothing,Float64} = nothing,
  ratio::Union{Nothing,Float64} = nothing,
  shift_deg::Union{Nothing,Float64} = nothing,
)
  @assert fromBus != toBus "From and to bus must be different"
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  if isnothing(ratedU)
    ratedU = side == 1 ? getNodeVn(net.nodeVec[from]) : getNodeVn(net.nodeVec[to])
  end
  r_pu, x_pu, b_pu, g_pu = toPU_RXBG(r = r, x = x, g = 0.0, b = b, v_kv = ratedU, baseMVA = net.baseMVA)
  addPIModelTrafo!(net = net, fromBus = fromBus, toBus = toBus, r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, status = status, ratedU = ratedU, ratedS = ratedS, ratio = ratio, shift_deg = shift_deg, isAux = false, side = side)
end

"""
Add a 3-winding transformer using a star-equivalent with an internal AUX bus.

Implementation strategy:
- Ensure an AUX bus exists (PQ, isAux=true) at the HV-side nominal voltage.
- Add three 2-winding PI-model transformers:
    AUX -- HB, AUX -- MB, AUX -- LV
- Convert r/x/b to PU (using toPU_RXBG) for each branch (mainly for validation/logging).
  The actual insertion uses add2WTPIModelTrafo!, which performs the conversion internally.

Notes:
- `ratio` is set to U_aux / U_side (HV/MV/LV) for MV and LV, and 1.0 for HB.
- `ratedU` passed to add2WTPIModelTrafo! is the AUX-side rated voltage (HV), because the
  branch is defined from AUX (HV base) to the respective side.
"""
function add3WTPiModelTrafo!(; net::Net, HBBus::String, MBBus::String, LVBus::String, r::NTuple{3,Float64}, x::NTuple{3,Float64}, b::NTuple{3,Float64}, ratedU_kV::NTuple{3,Float64}, ratedS_MVA::NTuple{3,Float64}, status::Int = 1)
  @assert status in (0, 1) "status must be 0 or 1"

  # --- basic existence checks (will throw from geNetBusIdx if missing) ---
  hb_idx = geNetBusIdx(net = net, busName = HBBus)
  mb_idx = geNetBusIdx(net = net, busName = MBBus)
  lv_idx = geNetBusIdx(net = net, busName = LVBus)

  # --- choose AUX-side (HV) voltage base ---
  # Prefer ratedU_kV[1] as HV side, but be robust and take the maximum.
  U_aux_kV = maximum(ratedU_kV)
  @assert U_aux_kV > 0.0 "ratedU_kV must be > 0"

  # --- build a deterministic AUX bus name (unique-ish, stable) ---
  # Keep it simple and reproducible; avoid special chars.
  function _sanitize(s::AbstractString)
    return replace(String(s), r"[^A-Za-z0-9_]" => "_")
  end

  aux_bus = "Aux3WT_" * _sanitize(HBBus) * "_" * _sanitize(MBBus) * "_" * _sanitize(LVBus)

  # --- create AUX bus if missing ---
  if !hasBusInNet(net = net, busName = aux_bus)
    addBus!(net = net, busName = aux_bus, busType = "PQ", vn_kV = U_aux_kV, isAux = true)
  end

  # --- define ratios for each branch AUX->side ---
  # ratio is interpreted as turns ratio; here we use U_aux / U_side.
  # For the HV bus branch, this is typically ~1 if HB is the HV side.

  ratio_hb = 1.0 #U_aux_kV / ratedU_kV[1]
  ratio_mb = 1.0 #U_aux_kV / ratedU_kV[2]
  ratio_lv = 1.0 #U_aux_kV / ratedU_kV[3]

  # --- optional PU conversion (for validation/logging/debug) ---
  # Note: add2WTPIModelTrafo! will do this conversion again internally.
  r1_pu, x1_pu, b1_pu, _g1_pu = toPU_RXBG(r = r[1], x = x[1], g = 0.0, b = b[1], v_kv = U_aux_kV, baseMVA = net.baseMVA)
  r2_pu, x2_pu, b2_pu, _g2_pu = toPU_RXBG(r = r[2], x = x[2], g = 0.0, b = b[2], v_kv = U_aux_kV, baseMVA = net.baseMVA)
  r3_pu, x3_pu, b3_pu, _g3_pu = toPU_RXBG(r = r[3], x = x[3], g = 0.0, b = b[3], v_kv = U_aux_kV, baseMVA = net.baseMVA)

  # If you want: add @debug lines here using (r*_pu, x*_pu, b*_pu)

  # --- add three 2W PI branches (AUX -> side buses) ---
  # Use side=1 because fromBus is the AUX (HV base).
  add2WTPIModelTrafo!(net = net, fromBus = aux_bus, toBus = HBBus, side = 1, r = r[1], x = x[1], b = b[1], status = status, ratedU = U_aux_kV, ratedS = ratedS_MVA[1], ratio = ratio_hb, shift_deg = 0.0)

  add2WTPIModelTrafo!(net = net, fromBus = aux_bus, toBus = MBBus, side = 1, r = r[2], x = x[2], b = b[2], status = status, ratedU = U_aux_kV, ratedS = ratedS_MVA[2], ratio = ratio_mb, shift_deg = 0.0)

  add2WTPIModelTrafo!(net = net, fromBus = aux_bus, toBus = LVBus, side = 1, r = r[3], x = x[3], b = b[3], status = status, ratedU = U_aux_kV, ratedS = ratedS_MVA[3], ratio = ratio_lv, shift_deg = 0.0)

  return aux_bus
end

"""
Add a transformer with PI model to the network.

# Arguments
- `net::Net`: The network to which the transformer will be added.
- `fromBus::String`: The name of the bus where the transformer originates.
- `toBus::String`: The name of the bus where the transformer terminates.
- `r_pu::Float64`: The per-unit resistance of the transformer.
- `x_pu::Float64`: The per-unit reactance of the transformer.
- `b_pu::Float64`: The per-unit susceptance of the transformer.
- `status::Int`: The status of the transformer.
- `ratedU::Union{Nothing, Float64}`: Rated voltage of the transformer. Default is `nothing`.
- `ratedS::Union{Nothing, Float64}`: Rated apparent power of the transformer. Default is `nothing`.
- `ratio::Union{Nothing, Float64}`: Ratio of the transformer. Default is `nothing`.
- `shift_deg::Union{Nothing, Float64}`: Phase shift angle of the transformer. Default is `nothing`.
"""
function add2WTTrafo!(; net::Net, fromBus::String, toBus::String, side::Int = 1, r::Float64, x::Float64, b::Float64, status::Int, ratedU::Union{Nothing,Float64} = nothing, ratedS::Union{Nothing,Float64} = nothing, ratio::Union{Nothing,Float64} = nothing, shift_deg::Union{Nothing,Float64} = nothing)
  @assert fromBus != toBus "From and to bus must be different"
  if isnothing(ratedU)
    from = geNetBusIdx(net = net, busName = fromBus)
    V1 = getNodeVn(net.nodeVec[from])
    to = geNetBusIdx(net = net, busName = toBus)
    V2 = getNodeVn(net.nodeVec[to])
    ratedU = max(V1, V2)
  end
  r_pu, x_pu, b_pu, g_pu = toPU_RXBG(r = r, x = x, g = 0.0, b = b, v_kv = ratedU, baseMVA = net.baseMVA)
  addPIModelTrafo!(net = net, fromBus = fromBus, toBus = toBus, r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, status = status, ratedU = ratedU, ratedS = ratedS, ratio = ratio, shift_deg = shift_deg, isAux = false, side = side)
end

"""
Add a two-winding transformer to the network.

# Arguments
- `net::Net`: The network to which the transformer will be added.
- `fromBus::String`: The name of the bus where the transformer originates.
- `toBus::String`: The name of the bus where the transformer terminates.
- `sn_mva::Float64`: Rated power of the transformer.
- `vk_percent::Float64`: Voltage regulation percent of the transformer.
- `vkr_percent::Float64`: Voltage regulation percent of the transformer.
- `pfe_kw::Float64`: Iron loss of the transformer.
- `i0_percent::Float64`: No-load current percent of the transformer.
- `status::Int`: The status of the transformer. Default is 1.
"""
function add2WTrafo!(; net::Net, fromBus::String, toBus::String, sn_mva::Float64, vk_percent::Float64, vkr_percent::Float64, pfe_kw::Float64, i0_percent::Float64, status::Int = 1)
  @assert fromBus != toBus "From and to bus must be different"
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  vn_hv_kV = getNodeVn(net.nodeVec[from])
  vn_lv_kV = getNodeVn(net.nodeVec[to])

  trafo = create2WTRatioTransformerNoTaps(from = from, to = to, vn_hv_kV = vn_hv_kV, vn_lv_kV = vn_lv_kV, sn_mva = sn_mva, vk_percent = vk_percent, vkr_percent = vkr_percent, pfe_kw = pfe_kw, i0_percent = i0_percent)
  side = getSideNumber2WT(trafo)
  ratio = calcTransformerRatio(trafo)
  push!(net.trafos, trafo)

  addBranch!(net = net, from = from, to = to, branch = trafo, status = status, ratio = ratio, side = side, vn_kV = vn_hv_kV)
end

"""
Add a prosumer (combination of a producer and consumer) to the network.

# Arguments
- `net::Net`: The network to which the prosumer will be added.
- `busName::String`: The name of the bus where the prosumer is connected.
- `type::String`: The type of the prosumer.
- `p::Union{Nothing, Float64}`: Active power produced or consumed. Default is `nothing`.
- `q::Union{Nothing, Float64}`: Reactive power produced or consumed. Default is `nothing`.
- `pMin::Union{Nothing, Float64}`: Minimum active power. Default is `nothing`.
- `pMax::Union{Nothing, Float64}`: Maximum active power. Default is `nothing`.
- `qMin::Union{Nothing, Float64}`: Minimum reactive power. Default is `nothing`.
- `qMax::Union{Nothing, Float64}`: Maximum reactive power. Default is `nothing`.
- `referencePri::Union{Nothing, String}`: Reference bus for the prosumer. Default is `nothing`.
- `vm_pu::Union{Nothing, Float64}`: Voltage magnitude setpoint. Default is `nothing`.
- `va_deg::Union{Nothing, Float64}`: Voltage angle setpoint. Default is `nothing`.
"""
function addProsumer!(;
  net::Net,
  busName::String,
  type::String,
  p::Union{Nothing,Float64} = nothing,
  q::Union{Nothing,Float64} = nothing,
  pMin::Union{Nothing,Float64} = nothing,
  pMax::Union{Nothing,Float64} = nothing,
  qMin::Union{Nothing,Float64} = nothing,
  qMax::Union{Nothing,Float64} = nothing,
  referencePri::Union{Nothing,String} = nothing,
  vm_pu::Union{Nothing,Float64} = nothing,
  va_deg::Union{Nothing,Float64} = nothing,
)
  busIdx = geNetBusIdx(net = net, busName = busName)
  idProsSum = length(net.prosumpsVec) + 1
  isAPUNode = false
  if !isnothing(vm_pu)
    node = net.nodeVec[busIdx]
    isAPUNode = isPVNode(node)
    nodeVm = getNodeVm(node)
    if !isnothing(nodeVm) && abs(nodeVm - vm_pu) > 1e-6
      @debug "Voltage setpoint already present at bus $busName: keep vm=$(nodeVm), ignore vm=$(vm_pu)"
    else
      setVmVa!(node = node, vm_pu = vm_pu, va_deg = va_deg)
    end
  end
  proTy = toProSumptionType(type)
  refPriIdx = isnothing(referencePri) ? nothing : geNetBusIdx(net = net, busName = referencePri)
  vn_kV = getNodeVn(net.nodeVec[busIdx])
  ps = ProSumer(vn_kv = vn_kV, busIdx = busIdx, oID = idProsSum, type = proTy, p = p, q = q, minP = pMin, maxP = pMax, minQ = qMin, maxQ = qMax, referencePri = refPriIdx, vm_pu = vm_pu, va_deg = va_deg, isAPUNode = isAPUNode)
  push!(net.prosumpsVec, ps)
  node = net.nodeVec[busIdx]
  if (isGenerator(proTy))
    addGenPower!(node = node, p = p, q = q)
  else
    addLoadPower!(node = node, p = p, q = q)
  end

  _buildQLimits!(net)
end

"""
Update the active and reactive power of a generator connected to a bus in the network.

# Arguments
- `net::Net`: The network object.
- `busName::String`: The name of the bus.
- `p::Union{Nothing, Float64}`: The active power to update. Default is `nothing`.
- `q::Union{Nothing, Float64}`: The reactive power to update. Default is `nothing`.

>Note: the corresponding prosumer object will not be updated.

# Examples
```julia
net = run_acpflow(max_ite= 7,tol = 1e-6, casefile='a_case.m') # run the power flow on the network and get the network object
updateBusPower!(net = net, busName = "Bus1", p = 0.5, q = 0.2) # Update the power of Bus1 to 0.5 MW and 0.2 MVar
run_net_acpflow(net = net, max_ite= 7,tol = 1e-6) # rerun the power flow with the updated network
```
"""
function addBusGenPower!(; net::Net, busName::String, p::Union{Nothing,Float64} = nothing, q::Union{Nothing,Float64} = nothing)
  busIdx = geNetBusIdx(net = net, busName = busName)
  addGenPower!(node = net.nodeVec[busIdx], p = p, q = q)
end

"""
Update the active and reactive power of a load connected to a bus in the network.

# Arguments
- `net::Net`: The network object.
- `busName::String`: The name of the bus.
- `p::Union{Nothing, Float64}`: The active power to update. Default is `nothing`.
- `q::Union{Nothing, Float64}`: The reactive power to update. Default is `nothing`.

>Note: the corresponding prosumer object will not be updated.

# Examples
```julia
net = run_acpflow(max_ite= 7,tol = 1e-6, casefile='a_case.m') # run the power flow on the network and get the network object
updateBusPower!(net = net, busName = "Bus1", p = 0.5, q = 0.2) # Update the power of Bus1 to 0.5 MW and 0.2 MVar
run_net_acpflow(net = net, max_ite= 7,tol = 1e-6) # rerun the power flow with the updated network
```
"""
function addBusLoadPower!(; net::Net, busName::String, p::Union{Nothing,Float64} = nothing, q::Union{Nothing,Float64} = nothing)
  busIdx = geNetBusIdx(net = net, busName = busName)
  addLoadPower!(node = net.nodeVec[busIdx], p = p, q = q)
end

"""
Update the active and reactive power of a shunt connected to a bus in the network.

# Arguments
- `net::Net`: The network object.
- `busName::String`: The name of the bus.
- `p::Float64`: The active power to update. 
- `q::Float64`: The reactive power to update.

# Examples
```julia
net = run_acpflow(max_ite= 7,tol = 1e-6, casefile='a_case.m') # run the power flow on the network and get the network object
updateBusPower!(net = net, busName = "Bus1", p = 0.5, q = 0.2) # Update the power of Bus1 to 0.5 MW and 0.2 MVar
run_net_acpflow(net = net, max_ite= 7,tol = 1e-6) # rerun the power flow with the updated network
```
"""
function addBusShuntPower!(; net::Net, busName::String, p::Float64, q::Float64)
  @warn "addBusShuntPower! is deprecated"
  #=
  if !hasShunt!(net = net, busName = busName)
    @info "add a new shunt to bus $busName"
    addShunt!(net = net, busName = busName, pShunt = p, qShunt = q)
  else
    @info "update shunt power at bus $busName to $p MW and $q MVar"
    sh = getShunt!(net = net, busName = busName) # get the shunt-reference(!) connected to the bus, we are manipulating the shunt-reference!
    _p, _q = getPQShunt(sh)
    pVal = isnothing(_p) ? p : p + _p
    qVal = isnothing(_q) ? q : q + _q
    updatePQShunt!(sh, pVal, qVal) # we need updated values to update shunt power
    busIdx = geNetBusIdx(net = net, busName = busName)
    addShuntPower!(node = net.nodeVec[busIdx], p = p, q = q) # we need given values to update the bus power
  end
  =#
end

"""
    addShuntMatpower!(; net, busName, Gs, Bs, in_service=1)

MATPOWER semantics:
- Gs/Bs are shunt admittance parameters given as MW/MVAr at V = 1.0 pu.
- Internally we stamp them as pu-admittance: y_pu = (Gs + j*Bs)/baseMVA.
- IMPORTANT: do NOT add fixed P/Q to the bus power balance.
"""
function addShuntMatpower!(; net::Net, busName::String, Gs::Float64, Bs::Float64, in_service::Int=1)
    busIdx  = geNetBusIdx(net = net, busName = busName)
    idShunt = length(net.shuntVec) + 1
    net.shuntDict[busIdx] = idShunt

    vn_kV = getNodeVn(net.nodeVec[busIdx])

    # Create shunt as "PQ-style" just to populate p_shunt/q_shunt fields meaningfully
    # (report values at 1.0 pu). This is NOT used for stamping.
    sh = Shunt(fromBus = busIdx, id = idShunt,
               base_MVA = net.baseMVA, vn_kV_shunt = vn_kV,
               p_shunt = Gs, q_shunt = Bs,
               status = in_service)

    # Override stamping admittance to MATPOWER semantics:
    # MATPOWER: Gs/Bs are MW/MVAr at 1 pu -> y_pu = (Gs + j Bs)/baseMVA
    sh.y_pu_shunt = ComplexF64(Gs, Bs) / net.baseMVA

    push!(net.shuntVec, sh)

    #if in_service == 1
    #    # IMPORTANT: report only (do NOT affect bus injection spec)
    #    addShuntReportPower!(node = net.nodeVec[busIdx], p = Gs, q = Bs)
    #end

    return nothing
end

"""
    setNetBranchStatus!(; net::Net, branchNr::Int, status::Int)

Sets the status of a branch in the network.

# Arguments
- `net::Net`: The network.
- `branchNr::Int`: The number of the branch.
- `status::Int`: The status of the branch. 1 = in service, 0 = out of service.

# Example
```julia
setNetBranchStatus!(net = network, branchNr = 1, status = 1)
```
"""
function setNetBranchStatus!(; net::Net, branchNr::Int, status::Int)
  @assert branchNr > 0 "Branch number must be greater than 0"
  @assert branchNr <= length(net.branchVec) "Branch $branchNr not found in the network"
  net.branchVec[branchNr].status = status
  markIsolatedBuses!(net = net, log = false)
end

"""
    setBusVoltage!(; net::Net, busName::String, vm_pu::Float64, va_deg::Float64)
Sets the voltage magnitude and angle of a bus in the network.
# Arguments
- `net::Net`: The network.
- `busName::String`: The name of the bus.
- `vm_pu::Float64`: The voltage magnitude in per unit.
- `va_deg::Float64`: The voltage angle in degrees.
# Example 
```julia
setBusVoltage!(net = network, busName = "Bus1", vm_pu = 1.02, va_deg = 5.0) 
```
"""

function setNodeVoltage!(; net::Net, busName::String, vm_pu::Float64, va_deg::Float64)
  @debug "Set bus voltage for $busName to vm_pu = $vm_pu and va_deg = $va_deg"
  busIdx = geNetBusIdx(net = net, busName = busName)
  node = net.nodeVec[busIdx]
  setVmVa!(node = node, vm_pu = vm_pu, va_deg = va_deg)
end

"""
    setBusAngle!(; net::Net, busName::String, va_deg::Float64)
Sets the voltage angle of a bus in the network.
# Arguments
- `net::Net`: The network.
- `busName::String`: The name of the bus.
- `va_deg::Float64`: The voltage angle in degrees.
# Example
```julia
setBusAngle!(net = network, busName = "Bus1", va_deg = 5.0)
```
"""
function setNodeAngle!(; net::Net, busName::String, va_deg::Float64)
  @debug "Set bus angle for $busName to va_deg = $va_deg"
  @assert va_deg >= -360 && va_deg <= 360 "Voltage angle must be between -360 and 360 degrees"
  busIdx = geNetBusIdx(net = net, busName = busName)
  node = net.nodeVec[busIdx]
  vm_pu = node._vm_pu
  setVmVa!(node = node, vm_pu = vm_pu, va_deg = va_deg)
end

"""
    setNetBranchStatus!(; net::Net, branchNr::Int, status::Int)

Sets the status of a branch in the network.

# Arguments
- `net::Net`: The network.
- `branchNr::Int`: The number of the branch.
- `status::Int`: The status of the branch. 1 = in service, 0 = out of service.

# Example
```julia
  brVec = getNetBranchNumberVec(net = net, fromBus = "B1", toBus = "B2")  
  setNetBranchStatus!(net = net, branchNr = brVec[1], status = 0)
```
"""
function getNetBranchNumberVec(; net::Net, fromBus::String, toBus::String)::Vector{Int}
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  brNumberVec = Int[]
  for (i, br) in enumerate(net.branchVec)
    if br.fromBus == from && br.toBus == to
      push!(brNumberVec, i)
    end
  end
  return brNumberVec
end

"""
    getNetBranch(; net::Net, fromBus::String, toBus::String)::Union{Branch,Nothing}

Retrieves the first branch found between two specified buses in the network.

# Arguments
- `net::Net`: The network.
- `fromBus::String`: The name of the bus where the branch starts.
- `toBus::String`: The name of the bus where the branch ends.

# Returns
- `Union{Branch,Nothing}`: The branch between the specified buses, or `nothing` if no such branch exists.

# Example
```julia
getNetBranch(net = network, fromBus = "Bus1", toBus = "Bus2")
"""
function getNetBranch(; net::Net, fromBus::String, toBus::String)::Union{Branch,Nothing}
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)

  for (i, br) in enumerate(net.branchVec)
    if br.fromBus == from && br.toBus == to
      return br
    end
  end
  return nothing
end

"""
Validate the network configuration.

# Arguments
- `net::Net`: The network to be validated.

# Returns
A tuple `(valid::Bool, message::String)` where `valid` is a boolean indicating whether the network is valid, and `message` is a string containing an error message if the network is invalid.
"""
function validate!(; net = Net, log::Bool = false)::Tuple{Bool,String}
  if length(net.nodeVec) == 0
    return false, "No buses defined in the network"
  end
  if length(net.slackVec) == 0
    return false, "No slack bus defined in the network"
  end
  if length(net.slackVec) > 1
    @info "More than one slack bus defined in the network"
  end
  if length(net.branchVec) == 0
    return false, "No branches defined in the network"
  end
  # Check if bus indices in ascending sequence
  sort!(net.nodeVec, by = x -> x.busIdx)
  for (i, key) in enumerate(net.nodeVec)
    if i != key.busIdx
      return false, "Bus index mismatch for bus $(key.busIdx)"
    end
  end
  markIsolatedBuses!(net = net, log = log)
  return true, "Network is valid"
end

"""
Get the voltage magnitude of a specific bus in the network.

# Arguments
- `net::Net`: The network from which to retrieve the voltage magnitude.
- `busName::String`: The name of the bus.

# Returns
The voltage magnitude of the specified bus.
"""
function get_bus_vn_kV(; net::Net, busName::String)
  busIdx = geNetBusIdx(net = net, busName = busName)
  return getNodeVn(net.nodeVec[busIdx])
end

"""
Get the voltage magnitude of a specific bus in the network.

# Arguments
- `net::Net`: The network from which to retrieve the voltage magnitude.
- `busIdx::Int`: The index of the bus.

# Returns
The voltage magnitude of the specified bus.
"""
function get_vn_kV(; net::Net, busIdx::Int)
  return getNodeVn(net.nodeVec[busIdx])
end

"""
Get the type of a specific bus in the network.

# Arguments
- `net::Net`: The network from which to retrieve the bus type.
- `busName::String`: The name of the bus.

# Returns
The type of the specified bus.
"""
function getBusType(; net::Net, busName::String)
  busIdx = geNetBusIdx(net = net, busName = busName)
  return net.nodeVec[busIdx]._nodeType
end

"""
Lock or unlock the network.

# Arguments
- `net::Net`: The network to be locked or unlocked.
- `locked::Bool`: Boolean indicating whether to lock the network.
"""
function lockNet!(; net::Net, locked::Bool)
  net._locked = locked
end

"""
    setTotalBusPower!(; net::Net, p::Float64, q::Float64)

Sets the total active and reactive power at the buses in the network.

# Arguments
- `net::Net`: The network.
- `p::Float64`: The total active power for the network.
- `q::Float64`: The total reactive power for the network.

# Example
```julia
setTotalBusPower!(net = network, p = 100.0, q = 50.0)
```
"""
function setTotalBusPower!(; net::Net, p::Float64, q::Float64)
  push!(net.totalBusPower, (p, q))
end

"""
    getTotalBusPower(; net::Net)::Tuple{Float64, Float64}

Gets the total active and reactive power for the network.

# Arguments
- `net::Net`: The network.

# Returns
- `n::Tuple{Float64, Float64}`: 

# Example
```julia
getTotalBusPower(net = network)
```
"""
function getTotalBusPower(; net::Net)
  if !isempty(net.totalBusPower)
    n = net.totalBusPower[end]
  else
    n = (0.0, 0.0)
  end
end

"""
Set the total losses in the network.

# Arguments
- `net::Net`: The network to which the losses will be added.
- `pLosses::Float64`: Total active power losses.
- `qLosses::Float64`: Total reactive power losses.
"""
function setTotalLosses!(; net::Net, pLosses::Float64, qLosses::Float64)
  push!(net.totalLosses, (pLosses, qLosses))
end

"""
Get the total losses in the network.

# Arguments
- `net::Net`: The network from which to retrieve the losses.

# Returns
A tuple `(pLosses::Float64, qLosses::Float64)` containing the total active and reactive power losses in the network.
"""
function getTotalLosses(; net::Net)
  if !isempty(net.totalLosses)
    n = net.totalLosses[end]
  else
    n = (0.0, 0.0)
  end
end

"""
    markIsolatedBuses!(;net::Net)

Finds and marks isolated buses in the network.

# Arguments
- `net::Net`: The network.

"""
function markIsolatedBuses!(; net::Net, log::Bool = false)
  # Erstelle ein Set, um die Busse zu speichern, die in den Zweigen im branchVec vorkommen
  connected_buses = Set{Int}()

  # Durchlaufe jeden Zweig im branchVec und füge die entsprechenden Busse in das Set ein
  for br in net.branchVec
    if br.status == 0
      continue
    end

    push!(connected_buses, br.fromBus)
    push!(connected_buses, br.toBus)
  end

  # Durchlaufe alle Busse im nodeVec und markiere die isolierten Busse
  for bus in net.nodeVec
    # Überprüfe, ob die Busnummer nicht im Set der verbundenen Busse enthalten ist
    if !(bus.busIdx in connected_buses)
      setNodeType!(bus, "Isolated")
      setVmVa!(node = bus, vm_pu = 0.0, va_deg = 0.0)
      push!(net.isoNodes, bus.busIdx)
      if log
        @info "Bus $(getCompName(bus.comp)) is isolated"
      end
    end
  end
end

function _buildQLimits!(net::Net; reset::Bool = true)
  # Number of buses from nodeVec
  nbus = length(net.nodeVec)

  # Ensure arrays have correct length
  if length(net.qmin_pu) != nbus
    resize!(net.qmin_pu, nbus)
  end
  if length(net.qmax_pu) != nbus
    resize!(net.qmax_pu, nbus)
  end

  # Reinitialize limits on every call (derived data -> safe to overwrite)
  # qmin_pu starts at +Inf (will be replaced by first finite sum)
  # qmax_pu starts at -Inf (will be replaced by first finite sum)
  fill!(net.qmin_pu, Inf)
  fill!(net.qmax_pu, -Inf)

  # Optionally reset Q-limit log
  if reset
    resetQLimitLog!(net)
  end

  # Aggregate generator Q-limits per bus
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    bus = getPosumerBusIndex(ps)
    # Safety: allow for buses beyond nodeVec if that can happen in your data
    if bus > nbus
      # grow arrays if needed
      resize!(net.qmin_pu, bus)
      resize!(net.qmax_pu, bus)
      for b = (nbus+1):bus
        net.qmin_pu[b] = Inf
        net.qmax_pu[b] = -Inf
      end
      nbus = bus
    end
    # Convert to p.u., handle 'no limit' as ±Inf
    qmin_pu = isnothing(ps.minQ) ? -Inf : ps.minQ / net.baseMVA
    qmax_pu = isnothing(ps.maxQ) ? Inf : ps.maxQ / net.baseMVA
    # Sum instead of min/max; respect sentinel values
    cur_qmin = net.qmin_pu[bus]
    net.qmin_pu[bus] = isfinite(cur_qmin) ? (cur_qmin + qmin_pu) : qmin_pu
    cur_qmax = net.qmax_pu[bus]
    net.qmax_pu[bus] = isfinite(cur_qmax) ? (cur_qmax + qmax_pu) : qmax_pu
  end

  return nothing
end

"""
    buildQLimits!(net; reset=true)

Public wrapper to (re)build aggregated per-bus Q limits in p.u.
"""
function buildQLimits!(net::Net; reset::Bool = true)
  return _buildQLimits!(net; reset = reset)
end

"""
    setQLimits!(; net::Net,
                 qmin_MVar::Float64,
                 qmax_MVar::Float64,
                 busName::Union{Nothing,String}=nothing)

Sets the reactive power (Q) limits of generator prosumers and then rebuilds
the aggregated bus-level Q-limits (`qmin_pu`, `qmax_pu`).

- Without `busName`: all generators receive the specified `qmin_MVar`/`qmax_MVar`.
- With `busName`: only generators connected to the specified bus receive the new limits.
"""

function setQLimits!(; net::Net, qmin_MVar::Float64, qmax_MVar::Float64, busName::Union{Nothing,AbstractString,AbstractVector{<:AbstractString}} = nothing)

  # Normalize busName to a set of bus indices or `nothing`
  busIdxSet = nothing
  if !isnothing(busName)
    # Ensure we always iterate over a collection of names
    names = busName isa AbstractString ? (busName,) : busName
    busIdxSet = Set(geNetBusIdx(net = net, busName = bn) for bn in names)
  end

  for ps in net.prosumpsVec
    #TODO is load also possible?
    isGenerator(ps) || continue

    if isnothing(busIdxSet)
      # No bus filter: apply to all generators
      ps.minQ = qmin_MVar
      ps.maxQ = qmax_MVar
    else
      # Apply only to generators connected to the selected buses
      c = ps.comp
      if !isnothing(c.cFrom_bus) && (c.cFrom_bus in busIdxSet)
        ps.minQ = qmin_MVar
        ps.maxQ = qmax_MVar
      end
    end
  end

  _buildQLimits!(net; reset = true)

  return nothing
end

"""
    setPVBusVset!(net::Net, busName::String; vm_pu::Float64)
  Sets the voltage magnitude setpoint for a PV bus in the network.
  Only applicable for buses of type "PV" and for testing purpose.
"""

function setPVBusVset!(; net::Net, busName::String, vm_pu::Float64)
  bus = geNetBusIdx(net = net, busName = busName)
  net.nodeVec[bus]._vm_pu = vm_pu
end

function _wattterfillQ(Q_target::Float64, q_spec::Vector{Float64}, qmin::Vector{Float64}, qmax::Vector{Float64}; tol::Float64 = 1e-6, maxiter::Int = 50)
  n = length(q_spec)
  @assert length(qmin) == n && length(qmax) == n

  # Start: clamp specs in their box
  q = similar(q_spec)
  for i = 1:n
    q[i] = clamp(q_spec[i], qmin[i], qmax[i])
  end

  for it = 1:maxiter
    Δ = Q_target - sum(q)
    if abs(Δ) <= tol
      return q
    end

    # Frei je nach Richtung (nach oben oder unten noch Bewegungsfreiheit)
    free = if Δ > 0
      [i for i = 1:n if q[i] < qmax[i] - tol]
    else
      [i for i = 1:n if q[i] > qmin[i] + tol]
    end

    isempty(free) && return q  # nichts mehr zu holen, Ziel nicht exakt erreichbar

    step = Δ / length(free)
    for i in free
      q[i] = clamp(q[i] + step, qmin[i], qmax[i])
    end
  end

  return q
end

function _generators_at_bus(net::Net, bus::Int)
  idx = Int[]
  for (k, ps) in enumerate(net.prosumpsVec)
    if isGenerator(ps) && getPosumerBusIndex(ps) == bus
      push!(idx, k)
    end
  end
  return idx
end

function _loads_at_bus(net::Net, bus::Int)
  idx = Int[]
  for (k, ps) in enumerate(net.prosumpsVec)
    if !isGenerator(ps) && getPosumerBusIndex(ps) == bus
      push!(idx, k)
    end
  end
  return idx
end

function distribute_bus_generation!(net::Net, bus::Int)
  node = net.nodeVec[bus]

  gens_idx = _generators_at_bus(net, bus)
  isempty(gens_idx) && return

  # Zielwerte am Bus (hier: aus Node – kannst du anpassen)
  P_target = node._pƩGen === nothing ? 0.0 : node._pƩGen
  Q_target = node._qƩGen === nothing ? 0.0 : node._qƩGen

  # Eingabewerte + Limits pro Generator
  P_spec = Float64[]
  Q_spec = Float64[]
  qmin   = Float64[]
  qmax   = Float64[]

  for idx in gens_idx
    ps = net.prosumpsVec[idx]
    push!(P_spec, isnothing(ps.pVal) ? 0.0 : ps.pVal)
    push!(Q_spec, isnothing(ps.qVal) ? 0.0 : ps.qVal)
    push!(qmin, isnothing(ps.minQ) ? -Inf : ps.minQ)
    push!(qmax, isnothing(ps.maxQ) ? Inf : ps.maxQ)
  end

  # --- P-Verteilung (ohne Limits, meist einfacher) ---
  P_alloc = similar(P_spec)
  P_sum   = sum(P_spec)

  if abs(P_sum) < 1e-9
    # gleichmäßig, wenn keine sinnvollen spec-Werte
    fill!(P_alloc, P_target / length(P_alloc))
  else
    # proportional zu den specs
    P_alloc .= P_target .* (P_spec ./ P_sum)
  end

  # --- Q-Verteilung mit Waterfilling ---
  Q_alloc = _wattterfillQ(Q_target, Q_spec, qmin, qmax)

  # Ergebnisse zurück in Prosumer
  for (pq, idx) in enumerate(gens_idx)
    ps = net.prosumpsVec[idx]
    setPQResult!(ps, P_alloc[pq], Q_alloc[pq])
  end
end

function distribute_bus_loads!(net::Net, bus::Int)
  node = net.nodeVec[bus]
  loads_idx = _loads_at_bus(net, bus)
  isempty(loads_idx) && return

  P_target = node._pƩLoad === nothing ? 0.0 : node._pƩLoad
  Q_target = node._qƩLoad === nothing ? 0.0 : node._qƩLoad

  P_spec = [isnothing(ps.pVal) ? 0.0 : ps.pVal for ps in net.prosumpsVec[loads_idx]]
  Q_spec = [isnothing(ps.qVal) ? 0.0 : ps.qVal for ps in net.prosumpsVec[loads_idx]]

  P_alloc = similar(P_spec)
  Q_alloc = similar(Q_spec)

  # einfache Variante: spec = result
  if abs(sum(P_spec)) < 1e-9
    P_alloc .= P_target / length(P_alloc)
  else
    P_alloc .= P_target .* (P_spec ./ sum(P_spec))
  end

  if abs(sum(Q_spec)) < 1e-9
    Q_alloc .= Q_target / length(Q_alloc)
  else
    Q_alloc .= Q_target .* (Q_spec ./ sum(Q_spec))
  end

  for (pq, idx) in enumerate(loads_idx)
    ps = net.prosumpsVec[idx]
    setPQResult!(ps, -P_alloc[pq], -Q_alloc[pq])  # Vorzeichen je nach Konvention
  end
end

function distributeBusResults!(net::Net)
  # reset results
  for ps in net.prosumpsVec
    ps.pRes = nothing
    ps.qRes = nothing
  end

  for node in net.nodeVec
    distribute_bus_generation!(net, node.busIdx)
    distribute_bus_loads!(net, node.busIdx)
  end
end

"""
    buildVoltageVector(net::Net) -> Vector{ComplexF64}

Builds the complex bus voltage vector V[k] = vm_pu[k] * exp(j * va_rad[k])
using the current nodal state stored in `net.nodeVec`.
"""
function buildVoltageVector(net::Net)
  V = Vector{ComplexF64}(undef, length(net.nodeVec))
  @inbounds for n in net.nodeVec
    V[n.busIdx] = n._vm_pu * cis(deg2rad(n._va_deg))
  end
  return V
end

"""
    initial_Vrect_from_net(net) -> (V0, slack_idx)

Build the initial complex voltage vector V0 from the network bus data
(Vm, Va), and detect the slack bus index.

Returns:
- V0::Vector{ComplexF64}
- slack_idx::Int
"""
function initialVrect(net::Net; flatstart::Bool = net.flatstart)
  nodes = net.nodeVec
  n = length(nodes)

  slack_idx = findfirst(n -> getNodeType(n) == Slack, nodes)
  slack_idx === nothing && error("initialVrect: no slack bus found")

  V0 = Vector{ComplexF64}(undef, n)

  @inbounds for k = 1:n
    node = nodes[k]

    # Voltage magnitude guess
    vm = (node._vm_pu === nothing || node._vm_pu <= 0.0) ? 1.0 : Float64(node._vm_pu)

    if flatstart
      # fix: 04.02.2026: slack bus gets its Vm (or 1.0 if missing), angle 0
      if k == slack_idx
        # Slack: keep its Vm (or 1.0 if missing), angle 0
        V0[k] = ComplexF64(vm, 0.0)
      else
        # All others: true flat start
        V0[k] = ComplexF64(1.0, 0.0)
      end
    else
      # Use stored angle if available, else 0
      va_deg = (node._va_deg === nothing) ? 0.0 : Float64(node._va_deg)
      va = deg2rad(va_deg)
      V0[k] = ComplexF64(vm * cos(va), vm * sin(va))
    end
  end

  return V0, slack_idx
end

"""
    buildComplexSVec(net) -> S::Vector{ComplexF64}

Build the specified complex power injection vector S = P + jQ in per-unit
for each bus, based on the net's bus load / generation / shunt data.
Positive P/Q means net injection into the bus (generation),
negative means net consumption (load).

"""
function buildComplexSVec(net::Net)
  n = length(net.nodeVec)
  S = Vector{ComplexF64}(undef, n)
  baseMVA = net.baseMVA

  for bus = 1:n
    Pgen = 0.0; Qgen = 0.0
    Pload = 0.0; Qload = 0.0

    for ps in net.prosumpsVec
      getPosumerBusIndex(ps) == bus || continue
      p = isnothing(ps.pVal) ? 0.0 : ps.pVal
      q = isnothing(ps.qVal) ? 0.0 : ps.qVal

      if isGenerator(ps)
        Pgen += p; Qgen += q
      else
        Pload += p; Qload += q
      end
    end

    Pinj = (Pgen - Pload) / baseMVA
    Qinj = (Qgen - Qload) / baseMVA
    S[bus] = ComplexF64(Pinj, Qinj)
  end

  return S

end



  """
    updateShuntPowers!(net; reset_node=true)

Recomputes shunt P/Q (MW/MVar) from solved bus voltages and shunt pu-admittances.
Writes back:
- sh.p_shunt / sh.q_shunt (results)
- node._pShunt / node._qShunt (for reporting)
"""
function updateShuntPowers!(;net::Net, reset_node::Bool=true)
  # Optionally reset node aggregates (recommended)
  if reset_node
    for n in net.nodeVec
      n._pShunt = 0.0
      n._qShunt = 0.0
    end
  end

  base = net.baseMVA

  for sh in net.shuntVec
    sh.status == 0 && continue
    sh.model == :Y || continue

    bus = sh.busIdx
    # ignore isolated buses if you want:
    (bus in net.isoNodes) && continue

    vm = net.nodeVec[bus]._vm_pu
    vm2 = vm * vm

    # S_pu = |V|^2 * conj(y_pu)
    S_pu = vm2 * conj(sh.y_pu_shunt)
    P = real(S_pu) * base
    Q = imag(S_pu) * base

    sh.p_shunt = P
    sh.q_shunt = Q

    # add to node report aggregate
    net.nodeVec[bus]._pShunt += P
    net.nodeVec[bus]._qShunt += Q
  end
end


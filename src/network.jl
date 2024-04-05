# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 01.04.2024
# include-file network.jl
"""
Net: Represents an electrical network.

Fields:
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
- `branchDict::Dict{Tuple{Int, Int},Int} : Dictionary mapping branch tuples to indices.
- `totalLosses::Vector{Tuple{Float64,Float64}}`: Vector containing tuples of total power losses.
- `_locked::Bool`: Boolean indicating if the network is locked.

Constructors:
- `Net(name::String, baseMVA::Float64, vmin_pu::Float64 = 0.9, vmax_pu::Float64 = 1.1)`: Creates a new `Net` object with the given name, base MVA, and optional voltage limits.

Functions:
- `addBus!(; net::Net, ...)`: Adds a bus to the network.
- `addShunt!(; net::Net, ...)`: Adds a shunt to the network.
- `addBranch!(; net::Net, ...)`: Adds a branch to the network.
- `addACLine!(; net::Net, ...)`: Adds an AC line segment to the network.
- `addPIModellTrafo!(; net::Net, ...)`: Adds a transformer with PI model to the network.
- `add2WTrafo!(; net::Net, ...)`: Adds a two-winding transformer to the network.
- `addProsumer!(; net::Net, ...)`: Adds a prosumer to the network.
- `validate(; net::Net)`: Validates the network.
- `get_bus_vn_kV(; net::Net, busName::String)`: Gets the voltage level of a bus by name.
- `get_vn_kV(; net::Net, busIdx::Int)`: Gets the voltage level of a bus by index.
- `getBusType(; net::Net, busName::String)`: Gets the type of a bus by name.
- `lockNet!(; net::Net, locked::Bool)`: Locks or unlocks the network.
- `setBranchStatus!(; net::Net, branchIdx::Int, status::Int)`: Sets the status of a branch.
- `setTotalLosses!(; net::Net, pLosses::Float64, qLosses::Float64)`: Sets the total losses of the network.
- `getTotalLosses(; net::Net)`: Gets the total losses of the network.
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
  branchDict::Dict{Tuple{Int, Int},Int}
  totalLosses::Vector{Tuple{Float64,Float64}}
  _locked::Bool

  function Net(; name::String, baseMVA::Float64, vmin_pu::Float64 = 0.9, vmax_pu::Float64 = 1.1)
    new(name, baseMVA, [], vmin_pu, vmax_pu, [], [], [], [], [], [], Dict{String,Int}(), Dict{Int,Int}(), Dict{Tuple{Int, Int},Int}(), [], false)
  end

  function Base.show(io::IO, o::Net)
    println(io, "Net: ", o.name)
    println(io, "Base MVA: ", o.baseMVA)
    println(io, "Slack: ", o.slackVec)
    println(io, "Vmin: ", o.vmin_pu)
    println(io, "Vmax: ", o.vmax_pu)
    println(io, "Nodes: ", length(o.nodeVec))
    println(io, "Lines: ", length(o.linesAC))
    println(io, "Transformers: ", length(o.trafos))
    println(io, "Branches: ", length(o.branchVec))
    println(io, "Prosumers: ", length(o.prosumpsVec))
    println(io, "Shunts: ", length(o.shuntVec))
  end
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
function addBus!(; net::Net, busName::String, busType::String, vn_kV::Float64, vm_pu::Float64 = 1.0, va_deg::Float64 = 0.0, vmin_pu::Union{Nothing,Float64} = nothing, vmax_pu::Union{Nothing,Float64} = nothing, isAux::Bool = false, 
                   oBusIdx::Union{Nothing,Int} = nothing, zone::Union{Nothing,Int} = nothing, area::Union{Nothing,Int} = nothing, ratedS::Union{Nothing,Float64} = nothing)
  
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

  if (uppercase(busType) == "SLACK")    
    push!(net.slackVec, busIdx)
  end

  node = Node(busIdx = busIdx, vn_kV = vn_kV, nodeType = toNodeType(busType), vm_pu = vm_pu, va_deg = va_deg, vmin_pu = vmin_pu, vmax_pu = vmax_pu, isAux = isAux, oBusIdx = oBusIdx, zone = zone, area = area, ratedS = ratedS)
  push!(net.nodeVec, node)
  
end

"""
addShunt!: Adds a shunt to the network.

Parameters:
- `net::Net`: Network object.
- `busName::String`: Name of the bus to which the shunt is added.
- `pShunt::Float64`: Active power of the shunt in MW.
- `qShunt::Float64`: Reactive power of the shunt in MVar.
- `in_service::Int = 1`: Indicator for shunt's in-service status (default is 1).
"""
function addShunt!(; net::Net, busName::String, pShunt::Float64, qShunt::Float64, in_service::Int = 1)  
  busIdx  = geNetBusIdx(net = net, busName = busName)  
  idShunt = length(net.shuntVec) + 1
    
  vn_kV = getNodeVn(net.nodeVec[busIdx])
  sh = Shunt(fromBus = busIdx, id = idShunt, base_MVA = net.baseMVA, Vn_kV_shunt = vn_kV, p_shunt = pShunt, q_shunt = qShunt, status = in_service)
  push!(net.shuntVec, sh)
  node = net.nodeVec[busIdx]
  if in_service == 1
    addShuntPower!(node = node, p = pShunt, q = qShunt)
  end  
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
"""
function addBranch!(; net::Net, from::Int, to::Int, branch::AbstractBranch, status::Int = 1, ratio::Union{Nothing,Float64} = nothing, side::Union{Nothing,Int} = nothing, vn_kV::Union{Nothing,Float64} = nothing)
  
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
  
  
  br = Branch(from = from, to = to, baseMVA = net.baseMVA, branch = branch, id = idBrunch, status = status, ratio = ratio, side = side, vn_kV = vn_kV, fromOid = fOrig, toOid = tOrig)
  push!(net.branchVec, br)
  index = length(net.branchVec)
  branchTupple = (from, to)
  net.branchDict[branchTupple] = index
end

"""
addACLine!: Adds an AC line segment to the network.

Parameters:
- `net::Net`: Network object.
- `fromBus::String`: Name of the "from" bus.
- `toBus::String`: Name of the "to" bus.
- `length::Float64`: Length of the line segment.
- `r::Float64`: Resistance of the line segment.
- `x::Float64`: Reactance of the line segment.
- `b::Union{Nothing,Float64} = nothing`: Susceptance of the line segment (default is nothing).
- `c_nf_per_km::Union{Nothing,Float64} = nothing`: Capacitance of the line segment in nF/km (default is nothing).
- `tanδ::Union{Nothing,Float64} = nothing`: Tangent of the loss angle (default is nothing).
- `ratedS::Union{Nothing, Float64}= nothing`: Rated power of the line segment in MVA (default is nothing).
- `status::Int = 1`: Status of the line segment (default is 1).
"""
function addACLine!(; net::Net, fromBus::String, toBus::String, length::Float64, r::Float64, x::Float64, b::Union{Nothing,Float64} = nothing, c_nf_per_km::Union{Nothing,Float64} = nothing, tanδ::Union{Nothing,Float64} = nothing, 
                      ratedS::Union{Nothing, Float64}= nothing, status::Int = 1)
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  vn_kV = getNodeVn(net.nodeVec[from])
  vn_2_kV = getNodeVn(net.nodeVec[to])
  @assert vn_kV == vn_2_kV "Voltage level of the from bus $(vn_kV) does not match the to bus $(vn_2_kV)"
  acseg = ACLineSegment(vn_kv = vn_kV, from = from, to = to, length = length, r = r, x = x, b = b, c_nf_per_km = c_nf_per_km, tanδ = tanδ, ratedS = ratedS)
  push!(net.linesAC, acseg)

  addBranch!(net = net, from = from, to = to, branch = acseg, vn_kV = vn_kV, status = status)
  
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
function addPIModellTrafo!(; net::Net, fromBus::String, toBus::String, r_pu::Float64, x_pu::Float64, b_pu::Float64, status::Int, ratedU::Union{Nothing,Float64}=nothing, ratedS::Union{Nothing,Float64} = nothing, 
                            ratio::Union{Nothing,Float64} = nothing, shift_deg::Union{Nothing,Float64} = nothing, isAux::Bool = false)
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  vn_hv_kV = getNodeVn(net.nodeVec[from])
  vn_lv_kV = getNodeVn(net.nodeVec[to])
  w1 = PowerTransformerWinding(vn_hv_kV, r_pu, x_pu, b_pu, 0.0, ratio, shift_deg, ratedU, ratedS, nothing, true )
  w2 = PowerTransformerWinding(vn_lv_kV, 0.0, 0.0)
  comp = getTrafoImpPGMComp(isAux, vn_hv_kV, from, to)
  trafo = PowerTransformer(comp, false, w1, w2, nothing, Sparlectra.PIModel)
  push!(net.trafos, trafo)
  
  addBranch!(net = net, from = from, to = to, branch = trafo, status = status, ratio = ratio, side = 1, vn_kV = vn_hv_kV)
  
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
function addProsumer!(; net::Net, busName::String, type::String, p::Union{Nothing,Float64} = nothing, q::Union{Nothing,Float64} = nothing, pMin::Union{Nothing,Float64} = nothing, pMax::Union{Nothing,Float64} = nothing, 
                        qMin::Union{Nothing,Float64} = nothing, qMax::Union{Nothing,Float64} = nothing, referencePri::Union{Nothing,String} = nothing, vm_pu::Union{Nothing,Float64} = nothing, va_deg::Union{Nothing,Float64} = nothing)
  busIdx = geNetBusIdx(net = net, busName = busName)  
  idProsSum = length(net.prosumpsVec) + 1
  isAPUNode = false   
  if !isnothing(vm_pu)
    node = net.nodeVec[busIdx]
    isAPUNode = isPVNode(node)
    setVmVa!(node = node, vm_pu = vm_pu, va_deg = va_deg)
  end
  proTy = toProSumptionType(type)
  refPriIdx = isnothing(referencePri) ? nothing : geNetBusIdx(net = net, busName = referencePri)
  vn_kV = getNodeVn(net.nodeVec[busIdx])
  prosumer = ProSumer(vn_kv = vn_kV, busIdx = busIdx, oID = idProsSum, type = proTy, p = p, q = q, minP=pMin, maxP=pMax, minQ=qMin, maxQ=qMax, referencePri = refPriIdx, vm_pu = vm_pu, va_deg = va_deg, isAPUNode = isAPUNode)
  push!(net.prosumpsVec, prosumer)

  node = net.nodeVec[busIdx]
  if (isGenerator(proTy))
    addGenPower!(node = node, p = p, q = q)
  else
    addLoadPower!(node = node, p = p, q = q)
  end
end

"""
Set the status of a branch in the network.

# Arguments
- `net::Net`: The network object.
- `fromBus::String`: The name of the bus where the branch originates.
- `toBus::String`: The name of the bus where the branch terminates.
- `status::Int`: The new status of the branch.

# Examples
```julia
net = run_acpflow(max_ite= 7,tol = 1e-6, casefile='a_case.m') # run the power flow on the network and get the network object
setBranchStatus!(net, "Bus1", "Bus2", 1)  # Set the status of the branch from Bus1 to Bus2 to 1.
run_net_acpflow(net = net, max_ite= 7,tol = 1e-6) # rerun the power flow with the updated network
```
"""
function setBranchStatus!(; net::Net, fromBus::String, toBus::String, status::Int)
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  BranchTupple = (from, to)
  index = net.branchDict[BranchTupple]  
  net.branchVec[index].status = status  
end

"""
Validate the network configuration.

# Arguments
- `net::Net`: The network to be validated.

# Returns
A tuple `(valid::Bool, message::String)` where `valid` is a boolean indicating whether the network is valid, and `message` is a string containing an error message if the network is invalid.
"""
function validate(; net = Net)
  
  if length(net.nodeVec) == 0    
    return false, "No buses defined in the network" 
  end
  if length(net.slackVec) == 0    
    return false, "No slack bus defined in the network"
  end
  if length(net.slackVec) > 1
    @warn "More than one slack bus defined in the network"    
  end
  if length(net.branchVec) == 0
    return false, "No branches defined in the network"
  end
  # Check if bus indices in ascending sequence
  sort!(net.nodeVec, by = x -> x.busIdx)
  for (i,key) in enumerate(net.nodeVec)
    if i != key.busIdx      
      return false, "Bus index mismatch for bus $(key.busIdx)"
    end
  end
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
# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 01.04.2024
# include-file network.jl
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
  qmin_pu::Vector{Float64}            # pro Bus Qmin (p.u.)
  qmax_pu::Vector{Float64}            # pro Bus Qmax (p.u.)
  qLimitEvents::Dict{Int,Symbol}      # BusIdx -> :min | :max (PV→PQ-Umschaltung)  
  _locked::Bool
  shuntDict::Dict{Int,Int}
  isoNodes::Vector{Int}
  qLimitLog::Vector{Any}
  cooldown_iters::Int 
  q_hyst_pu::Float64  
  #! format: off
  function Net(; name::String, baseMVA::Float64, vmin_pu::Float64 = 0.9, vmax_pu::Float64 = 1.1)    


    # letzte Zeile im new(...) ersetzen/erweitern:
    new(name, baseMVA, [], vmin_pu, vmax_pu, [], [], [], [], [], [], 
        Dict{String,Int}(), Dict{Int,Int}(), [], [],                  # totalLosses, totalBusPower
        Float64[], Float64[], Dict{Int,Symbol}(),                     # qmin_pu, qmax_pu, qLimitEvents
        false, Dict{Int,Int}(), [],                                   # _locked, shuntDict, isoNodes
        Any[],                                                        # qLimitLog 
        0,                                                            # cooldown_iters (0=aus)
        0.0)                                                          # q_hyst_pu    


  end
  #! format: on
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
    if !isempty(o.qLimitEvents)
        println(io, "Q-limit events: ", length(o.qLimitEvents))
    end
    if !isempty(o.qLimitLog)
        println(io, "Q-limit log entries: ", length(o.qLimitLog))
    end

  end
end
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

# --- positionsbasierte Kern-APIs ---
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
  busIdx = geNetBusIdx(net = net, busName = busName)
  idShunt = length(net.shuntVec) + 1
  net.shuntDict[busIdx] = idShunt
  vn_kV = getNodeVn(net.nodeVec[busIdx])
  sh = Shunt(fromBus = busIdx, id = idShunt, base_MVA = net.baseMVA, vn_kV_shunt = vn_kV, p_shunt = pShunt, q_shunt = qShunt, status = in_service)
  push!(net.shuntVec, sh)
  if in_service == 1
    addShuntPower!(node = net.nodeVec[busIdx], p = pShunt, q = qShunt)
  end
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
"""
function addBranch!(; net::Net, from::Int, to::Int, branch::AbstractBranch, status::Integer = 1, ratio::Union{Nothing,Float64} = nothing, side::Union{Nothing,Int} = nothing, vn_kV::Union{Nothing,Float64} = nothing)
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
  br = Branch(branchIdx = idBrunch, from = from, to = to, baseMVA = net.baseMVA, branch = branch, id = idBrunch, status = status, ratio = ratio, side = side, vn_kV = vn_kV, fromOid = fOrig, toOid = tOrig)
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
  @debug "Adding AC line segment from $fromBus to $toBus"
  @assert fromBus != toBus "From and to bus must be different"
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  vn_kV = getNodeVn(net.nodeVec[from])
  vn_2_kV = getNodeVn(net.nodeVec[to])
  @assert vn_kV == vn_2_kV "Voltage level of the from bus $(vn_kV) does not match the to bus $(vn_2_kV)"
  acseg = ACLineSegment(vn_kv = vn_kV, from = from, to = to, length = length, r = r, x = x, b = b, c_nf_per_km = c_nf_per_km, tanδ = tanδ, ratedS = ratedS, paramsBasedOnLength = false, isPIModel = false)
  push!(net.linesAC, acseg)

  addBranch!(net = net, from = from, to = to, branch = acseg, vn_kV = vn_kV, status = status)
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
  acseg = ACLineSegment(vn_kv = vn_kV, from = from, to = to, length = 1.0, r = r_pu, x = x_pu, b = b_pu, ratedS = ratedS, paramsBasedOnLength = true, isPIModel = true)
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
    setVmVa!(node = node, vm_pu = vm_pu, va_deg = va_deg)
  end
  proTy = toProSumptionType(type)
  refPriIdx = isnothing(referencePri) ? nothing : geNetBusIdx(net = net, busName = referencePri)
  vn_kV = getNodeVn(net.nodeVec[busIdx])
  prosumer = ProSumer(vn_kv = vn_kV, busIdx = busIdx, oID = idProsSum, type = proTy, p = p, q = q, minP = pMin, maxP = pMax, minQ = qMin, maxQ = qMax, referencePri = refPriIdx, vm_pu = vm_pu, va_deg = va_deg, isAPUNode = isAPUNode)
  push!(net.prosumpsVec, prosumer)

  node = net.nodeVec[busIdx]
  if (isGenerator(proTy))
    addGenPower!(node = node, p = p, q = q)
  else
    addLoadPower!(node = node, p = p, q = q)
  end
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
  @debug "Set branch status to $fromBusSwitch and $toBusSwitch"

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

"""
    setBusGeneratorQLimits!(; net::Net, busName::String, qmin_MVar::Float64, qmax_MVar::Float64)

Setzt für **alle Generator-ProSumer** am Bus `busName` die reaktiven Grenzwerte
in MVar und baut anschließend die Bus-aggregierten Q-Limits neu auf.
"""
function setBusGeneratorQLimits!(; net::Net, busName::String, qmin_MVar::Float64, qmax_MVar::Float64)
    busIdx = geNetBusIdx(net = net, busName = busName)
    for ps in net.prosumpsVec
        if isGenerator(ps)
            c = ps.comp
            if !isnothing(c.cFrom_bus) && c.cFrom_bus == busIdx
                ps.minQ = qmin_MVar
                ps.maxQ = qmax_MVar
            end
        end
    end
    buildQLimits!(net)   # aggregierte per-Bus-Limits (p.u.) neu aufbauen
    return nothing
end

"""
    setPVGeneratorQLimitsAll!(; net::Net, qmin_MVar::Float64, qmax_MVar::Float64)

Setzt auf **allen PV-Bussen** die Generator-Q-Grenzen (für alle Generator-ProSumer,
die an einem PV-Bus hängen) und baut anschließend die Aggregation neu auf.
"""
function setPVGeneratorQLimitsAll!(; net::Net, qmin_MVar::Float64, qmax_MVar::Float64)
    pv_set = Set{Int}()
    for n in net.nodeVec
        if isPVNode(n)
            push!(pv_set, n.busIdx)
        end
    end
    for ps in net.prosumpsVec
        if isGenerator(ps)
            c = ps.comp
            if !isnothing(c.cFrom_bus) && (c.cFrom_bus in pv_set)
                ps.minQ = qmin_MVar
                ps.maxQ = qmax_MVar
            end
        end
    end
    buildQLimits!(net)
    return nothing
end

"""
    rebuildQLimits!(; net::Net)

Bequemer Wrapper auf `buildQLimits!`.
"""
rebuildQLimits!(; net::Net) = (buildQLimits!(net); nothing)


# ------------------------------------------------------------
# Helpers to tweak PV setpoints and Q-limits *without* mutating Net fields
# ------------------------------------------------------------

# Robust bus lookup by name (falls Knotenname mal "busName" oder "name" heißt)
# Fallback: akzeptiere "B<idx>"-Strings.
function _find_bus_index_by_name(net::Net, busName::String)::Int
    idx = findfirst(n ->
        (hasproperty(n, :busName) && n.busName == busName) ||
        (hasproperty(n, :name)    && n.name    == busName), net.nodeVec)
    if isnothing(idx)
        # Allow "B12" → 12
        if startswith(busName, "B")
            parsed = tryparse(Int, busName[2:end])
            if parsed !== nothing && 1 <= parsed <= length(net.nodeVec)
                return parsed
            end
        end
        error("Bus name '$busName' not found in net.nodeVec.")
    end
    return idx
end

# ------------------------------------------------------------
# Setze den Spannungs-Sollwert (Vset) eines PV-Busses über das Node-Feld _vm_pu.
# Diese Information greift der Full-NR über Vset = node._vm_pu (falls vorhanden) ab.
# ------------------------------------------------------------
function setPVBusVset!(net::Net; bus::Int, vm_pu::Float64)
    node = net.nodeVec[bus]
    if !isPVNode(node)   # richtiger PV/Slack-Check
        @warn "setPVBusVset!: Bus $bus is not PV; setting _vm_pu anyway (the solver will ignore it for non-PV)."
    end
    setfield!(node, :_vm_pu, vm_pu)
    return nothing
end

function setPVBusVset!(net::Net, busName::String; vm_pu::Float64)
    bus = geNetBusIdx(net = net, busName = busName)
    #bus = _find_bus_index_by_name(net, busName)
    return setPVBusVset!(net; bus=bus, vm_pu=vm_pu)
end

# ------------------------------------------------------------
# Setze Q-Limits (MVar) für alle Generator-ProSumer an einem Bus.
# Achtung: Werte sind in MVar (nicht p.u.); die Aggregation nach p.u. macht limits.jl.
# ------------------------------------------------------------
function setPVBusQLimits!(net::Net; bus::Int, qmin_MVar::Float64, qmax_MVar::Float64)
    if qmax_MVar < qmin_MVar
        error("setPVBusQLimits!: qmax_MVar < qmin_MVar (got $qmax_MVar < $qmin_MVar)")
    end
    for ps in net.prosumpsVec
        try
            if isGenerator(ps) && Sparlectra._prosumer_bus_index(ps) == bus
                setfield!(ps, :minQ, qmin_MVar)
                setfield!(ps, :maxQ, qmax_MVar)
            else
              @warn "setPVBusQLimits!: skipping a non-generator or prosumer at a different bus."
            end
        catch err
            @warn "setPVBusQLimits!: skipping a prosumer due to $err"
        end
    end
    return nothing
end

function setPVBusQLimits!(net::Net, busName::String; qmin_MVar::Float64, qmax_MVar::Float64)
    bus = _find_bus_index_by_name(net, busName)
    return setPVBusQLimits!(net; bus=bus, qmin_MVar=qmin_MVar, qmax_MVar=qmax_MVar)
end

# ------------------------------------------------------------
# Optionales Convenience: Setze denselben Vset für *alle* PV-Busse.
# ------------------------------------------------------------
function setAllPVVset!(net::Net; vm_pu::Float64)
    for n in net.nodeVec
        if hasproperty(n, :type) && n.type == PV
            setfield!(n, :_vm_pu, vm_pu)
        end
    end
    return nothing
end

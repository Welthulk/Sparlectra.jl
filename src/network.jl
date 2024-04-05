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
  _locked::Bool

  function Net(; name::String, baseMVA::Float64, vmin_pu::Float64 = 0.9, vmax_pu::Float64 = 1.1)
    new(name, baseMVA, [], vmin_pu, vmax_pu, [], [], [], [], [], [], Dict{String,Int}(), Dict{Int,Int}(), [], false)
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

function hasBusInNet(; net::Net, busName::String)::Bool
  return haskey(net.busDict, busName)
end

function geNetBusIdx(; net::Net, busName::String)::Int
  if !haskey(net.busDict, busName)
    error("Bus $(busName) not found in the network")
  end   
  return net.busDict[busName]
end

function getNetOrigBusIdx(; net::Net, busName::String)::Int
  id = geNetBusIdx(net = net, busName = busName)
  return net.busOrigIdxDict[id]
end

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
end

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

# Transformer from Matpower Data
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

function get_bus_vn_kV(; net::Net, busName::String)
  busIdx = geNetBusIdx(net = net, busName = busName)
  return getNodeVn(net.nodeVec[busIdx])
end

function get_vn_kV(; net::Net, busIdx::Int)
  return getNodeVn(net.nodeVec[busIdx])
end

function getBusType(; net::Net, busName::String)
  busIdx = geNetBusIdx(net = net, busName = busName)
  return net.nodeVec[busIdx]._nodeType
end

function lockNet!(; net::Net, locked::Bool)  
  net._locked = locked
end

function setTotalLosses!(; net::Net, pLosses::Float64, qLosses::Float64)
  push!(net.totalLosses, (pLosses, qLosses))
end

function getTotalLosses(; net::Net)
  if !isempty(net.totalLosses)
    n = net.totalLosses[end]
  else
    n = (0.0, 0.0)
  end
end
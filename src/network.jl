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
  totalLosses::Vector{Tuple{Float64, Float64}}

  function Net(
    name,
    baseMVA::Float64,
    slack::Integer,
    nodeVec::Vector{Node},
    linesAC::Vector{ACLineSegment},
    trafos::Vector{PowerTransformer},
    branchVec::Vector{Branch},
    prosumpsVec::Vector{ProSumer},
    shuntVec::Vector{Shunt},
  )
    slackVec = Vector{Int}()
    push!(slackVec, slack)
    new(name, baseMVA, slackVec, 0.9, 1.1, nodeVec, linesAC, trafos, branchVec, prosumpsVec, shuntVec, Dict{String,Int}(),[])
  end

  function Net(; name::String, baseMVA::Float64, vmin_pu::Float64 = 0.9, vmax_pu::Float64 = 1.1)
    new(name, baseMVA, [], vmin_pu, vmax_pu, [], [], [], [], [], [], Dict{String,Int}(), [])
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

function geNetBusIdx(; net::Net, busName::String)::Int
  try
    return net.busDict[busName]
  catch
    error("Bus $busName not found in the network")
  end
end

function addBus!(; net::Net, busName::String, busType::String, vn_kV::Float64, vm_pu::Float64 = 1.0, va_deg::Float64 = 0.0, isAux::Bool = false)
  @assert busType in ["Slack", "PQ", "PV"] "Invalid bus type: $busType"
  busIdx = length(net.nodeVec) + 1
  node = Node(busIdx = busIdx, vn_kV = vn_kV, nodeType = toNodeType(busType), vm_pu = vm_pu, va_deg = va_deg, vmin_pu = net.vmin_pu, vmax_pu = net.vmax_pu, isAux = isAux)
  push!(net.nodeVec, node)
  net.busDict[busName] = busIdx
  if (busType == "Slack")
    push!(net.slackVec, busIdx)
  end
end

function addShunt!(; net::Net, busName::String, pShunt::Float64, qShunt::Float64, in_service::Int = 1)
  @assert length(net.nodeVec) > 0 "No buses defined in the network"
  idShunt = length(net.shuntVec) + 1
  busIdx = geNetBusIdx(net = net, busName = busName)
  vn_kV = getNodeVn(net.nodeVec[busIdx])
  sh = Shunt(fromBus = busIdx, id = idShunt, base_MVA = net.baseMVA, Vn_kV_shunt = vn_kV, p_shunt = pShunt, q_shunt = qShunt, status = in_service)
  push!(net.shuntVec, sh)
  node = net.nodeVec[busIdx]
  addShuntPower!(node = node, p = pShunt, q = qShunt)
end

function addBranch!(; net::Net, from::Int, to::Int, branch::AbstractBranch, status::Integer = 1, ratio::Union{Nothing,Float64} = nothing, side::Union{Nothing,Int} = nothing, vn_kV::Union{Nothing,Float64} = nothing)
  idBrunch = length(net.branchVec) + 1
  br = Branch(from = from, to = to, baseMVA = net.baseMVA, branch = branch, id = idBrunch, status = status, ratio = ratio, side = side, vn_kV = vn_kV)
  push!(net.branchVec, br)
end

function addACLine!(; net::Net, fromBus::String, toBus::String, length::Float64, r::Float64, x::Float64, c_nf_per_km::Union{Nothing,Float64} = nothing, tanδ::Union{Nothing,Float64} = nothing)
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  vn_kV = getNodeVn(net.nodeVec[from])
  vn_2_kV = getNodeVn(net.nodeVec[to])
  @assert vn_kV == vn_2_kV "Voltage level of the from bus $(vn_kV) does not match the to bus $(vn_2_kV)"
  acseg = ACLineSegment(vn_kv = vn_kV, from = from, to = to, length = length, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  addBranch!(net = net, from = from, to = to, branch = acseg, status = 1, vn_kV = vn_kV)
  push!(net.linesAC, acseg)
end

function add2WTrafo!(; net::Net, fromBus::String, toBus::String, sn_mva::Float64, vk_percent::Float64, vkr_percent::Float64, pfe_kw::Float64, i0_percent::Float64)
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  vn_hv_kV = getNodeVn(net.nodeVec[from])
  vn_lv_kV = getNodeVn(net.nodeVec[to])
  #@assert n_hv_kV == vn_hv_kV "Voltage level of the from bus $(n_hv_kV) does not match the transformer's high voltage side $(vn_hv_kV)"
  #@assert n_lv_kV == vn_lv_kV "Voltage level of the to bus $(n_lv_kV) does not match the transformer's low voltage side $(vn_lv_kV)"

  trafo = create2WTRatioTransformerNoTaps(from = from, to = to, vn_hv_kv = vn_hv_kV, vn_lv_kv = vn_lv_kV, sn_mva = sn_mva, vk_percent = vk_percent, vkr_percent = vkr_percent, pfe_kw = pfe_kw, i0_percent = i0_percent)
  side = getSideNumber2WT(trafo)
  ratio = calcTransformerRatio(trafo)
  addBranch!(net = net, from = from, to = to, branch = trafo, status = 1, ratio = ratio, side = side, vn_kV = vn_hv_kV)
  push!(net.trafos, trafo)
end

function addProsumer!(; net::Net, busName::String, type::String, p::Union{Nothing,Float64} = nothing, q::Union{Nothing,Float64} = nothing, 
                        referencePri::Union{Nothing,String} = nothing, vm_pu::Union{Nothing,Float64} = nothing, va_deg::Union{Nothing,Float64} = nothing)

  busIdx = geNetBusIdx(net = net, busName = busName)
  idProsSum = length(net.prosumpsVec) + 1
  isAPUNode = false
  if !isnothing(vm_pu) && !isnothing(va_deg)
    isAPUNode = true
    node = net.nodeVec[busIdx]
    setVmVa!(node = node, vm_pu = vm_pu, va_deg = va_deg)
  end
  proTy = toProSumptionType(type)
  refPriIdx = isnothing(referencePri) ? nothing : geNetBusIdx(net = net, busName = referencePri)
  vn_kV = getNodeVn(net.nodeVec[busIdx])
  prosumer = ProSumer(vn_kv = vn_kV, busIdx = busIdx, oID = idProsSum, type = proTy, p = p, q = q, referencePri = refPriIdx, vm_pu = vm_pu, va_deg = va_deg, isAPUNode = isAPUNode)
  push!(net.prosumpsVec, prosumer)

  node = net.nodeVec[busIdx]
  if (isGenerator(proTy))
    addGenPower!(node = node, p = p, q = q)
  else
    addLoadPower!(node = node, p = p, q = q)
  end
end

function validate(; net = Net)
  val = true
  if length(net.nodeVec) == 0
    @warn "No buses defined in the network"
    val = false
  end
  if length(net.slackVec) == 0
    @warn "No slack bus defined in the network"
    val = false
  end
  if length(net.slackVec) > 1
    @warn "More than one slack bus defined in the network"    
  end
  if length(net.branchVec) == 0
    @warn "No branches defined in the network"
    val = false
  end  
  return val
end

function setTotalLosses!(; net::Net, pLosses::Float64, qLosses::Float64)  
  push!(net.totalLosses, (pLosses, qLosses))
end

function getTotalLosses(;net::Net)
  if !isempty(net.totalLosses)
    n = net.totalLosses[end]
  else
    n = (0.0, 0.0)
  end  
end
# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 20.6.2023
# CIGRE HV network
# Source of this network can be found here (Task Force C6.04.02 ): https://www.researchgate.net/publication/271963972_TF_C60402_TB_575_--_Benchmark_Systems_for_Network_Integration_of_Renewable_and_Distributed_Energy_Resources

using Sparlectra
using Sparlectra.ResDataTypes
using Sparlectra.SparlectraNet
using BenchmarkTools
using Logging

function test_acpflow(verbose::Int = 0)
  myNet = testCreateNetworkFromScratch(false)
  result = true
  
  sparse = false
  Y = createYBUS(myNet.branchVec, myNet.shuntVec, (verbose > 2), sparse)
  iterations = 8
  tol = 1e-6
  etime = @elapsed begin
    ite, erg = calcNewtonRaphson!(Y, myNet.nodeVec, myNet.baseMVA, iterations, tol, verbose, sparse)
  end

  if erg == 0
    result = result & true
    if verbose > 0 
      calcNetLosses!(myNet.nodeVec, myNet.branchVec, myNet.baseMVA, (verbose > 1))
      printACPFlowResults(myNet, etime, ite)
    end
  else
    result = result & false
  end

  return result
end

function test_NBI_MDO()
  myNet = testCreateNetworkFromScratch(false)
  result = true

  nodeNumberVec = Vector{Int}()
  branchTupleSet = Set{Tuple}()

  for n in myNet.nodeVec
    if n.busIdx in nodeNumberVec
      @warn "node", n, "is already in set!"
    else
      push!(nodeNumberVec, n.busIdx)
    end
  end

  parallel_branches_found = 0
  for b in myNet.branchVec
    tupple = (b.fromBus, b.toBus)
    if b.isParallel == true
      parallel_branches_found += 1
    elseif tupple in branchTupleSet
      b.isParallel = true
      parallel_branches_found += 1
      @warn "branch", b, "is parallel, but not marked as parallel!"
    else
      push!(branchTupleSet, (b.fromBus, b.toBus))
    end
  end
  result = result & (parallel_branches_found = 1)

  A = getNBI(nodeNumberVec, branchTupleSet)
  order = mdoRCM(length(nodeNumberVec), branchTupleSet)
  if order == [9, 10, 11, 12, 13, 2, 5, 1, 6, 4, 3, 7, 8]
    result = result & true
  else
    @warn "expected MCO-Order: [9, 10, 11, 12, 13, 2, 5, 1, 6, 4, 3, 7, 8]"
    result = false
  end

  return result
end

function testNetwork()::Bool
  myNet = testCreateNetworkFromScratch(false)
  result = true
  result = result & checkNodeIdsUnique(myNet.nodeVec)
  result = result & checkNodeConnections(myNet.nodeVec)
  result = result & checkTransformerIdIsUnique(myNet.trafos)
  result = result & checkLineIdIsUnique(myNet.linesAC)
  result = result & checkProSumerIdIsUnique(myNet.prosumpsVec)
  return result
end

function testCreateNetworkFromScratch(verbose::Bool = false)

  nodeVec = Vector{Node}()
  shuntVec = Vector{Shunt}()
  ACLines = Vector{ACLineSegment}()
  Transformers = Vector{PowerTransformer}()
  proSumersVec = Vector{ProSumer}()
  branchVec = Vector{Branch}()
  proSumersVec = Vector{ProSumer}()
  
  busDict = Dict{String, Int}()
  idBus = 0
  idShunt = 0
  idBrunch = 0
  idProsSum = 0
  Sbase_MVA = 1000.0
  netName = "cigre"

  # length of coupling bus 6b-6a
  lLine_6a6b = 1.0

  function reports()
    @info "Report:"
    @info "Nodes: "
    anzNodes = 0
    for n in nodeVec
      anzNodes += 1
      reportCrusialData(n)
    end
    @info "Number of Nodes: ", anzNodes

    @info "Lines: "
    anzL = 0
    for acseg in ACLines
      anzL += 1
      reportCrusialData(acseg)
    end
    @info "Number of Lines: ", anzL

    @info "Trafos: "
    anzT = 0
    for t in Transformers
      anzT += 1
      reportCrusialData(t)
    end
    @info "Number of Trafos: ", anzT

    @info "Prosumers: "
    anzP = 0
    for p in proSumersVec
      anzP += 1
      reportCrusialData(p)
    end
    @info "Number of Prosumers: ", anzP
  end

  function addBusDict!(busIdx::Int, busName::String)
    busDict[busName] = busIdx
  end

  function getBusIdx(busName::String)::Int
    return busDict[busName]    
  end
  
  

  
  function addBus!(;busName::String, busType::String, vn_kV::Float64, 
                   pLoad::Union{Nothing, Float64}= nothing, qLoad::Union{Nothing,Float64}= nothing, 
                   pGen::Union{Nothing, Float64}= nothing, qGen::Union{Nothing, Float64}= nothing, 
                   pShunt::Union{Nothing, Float64}= nothing, qShunt::Union{Nothing, Float64}= nothing, 
                   vm::Union{Nothing, Float64}= nothing, va::Union{Nothing, Float64}= nothing, isAux::Bool = false)
    vmin=0.9
    vmax=1.1
    idBus += 1
    node = Node(busIdx = idBus, vn_kV = vn_kV, nodeType = toNodeType(busType), vm_pu = vm, va_deg = va, vmin_pu = vmin, vmax_pu = vmax, pƩLoad = pLoad, qƩLoad = qLoad, pShunt = pShunt, qShunt = qShunt, pƩGen = pGen, qƩGen = qGen, isAux = isAux)
    push!(nodeVec, node)   
    addBusDict!(idBus, busName)
  end

  function addShunt!(;busIdx::Int, pShunt::Float64, qShunt::Float64, vn::Float64, in_service::Int = 1)
    idShunt += 1
    sh = Shunt(fromBus = busIdx, id = idShunt, base_MVA = Sbase_MVA, Vn_kV_shunt = vn, p_shunt = pShunt, q_shunt = qShunt, status = in_service)
    push!(shuntVec, sh)
  end

  function addBranch!(;from::Int, to::Int, baseMVA::Float64, branch::AbstractBranch, status::Integer=1,  ratio::Union{Nothing,Float64}=nothing,  side::Union{Nothing,Int}=nothing, vn_kV::Union{Nothing,Float64}=nothing )
    idBrunch += 1
    br = Branch(from = from, to = to, baseMVA = baseMVA, branch = branch, id = idBrunch, status = status, ratio = ratio, side = side, vn_kV = vn_kV)
    push!(branchVec, br)
  end

  function addACLine!(;vn_kv::Float64, from::Int, to::Int, length::Float64, r::Float64, x::Float64, c_nf_per_km::Union{Nothing,Float64} = nothing, tanδ::Union{Nothing,Float64} = nothing)
    acseg = ACLineSegment(vn_kv = vn_kv, from = from, to = to, length = length, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
    addBranch!(from = from, to = to, baseMVA = Sbase_MVA, branch = acseg, status = 1, vn_kV = vn_kv)
    push!(ACLines, acseg)
    
  end

  function add2WTrafo!(;from::Int, to::Int, sn_mva::Float64, vn_hv_kv::Float64, vn_lv_kv::Float64, vk_percent::Float64, vkr_percent::Float64, pfe_kw::Float64, i0_percent::Float64)    
    trafo = create2WTRatioTransformerNoTaps(from = from, to = to, vn_hv_kv = vn_hv_kv, vn_lv_kv = vn_lv_kv, sn_mva = sn_mva, vk_percent = vk_percent, vkr_percent = vkr_percent, pfe_kw = pfe_kw, i0_percent = i0_percent)
    side = getSideNumber2WT(trafo)
    ratio = calcTransformerRatio(trafo)    
    addBranch!(from = from, to = to, baseMVA = Sbase_MVA, branch = trafo, status = 1, ratio = ratio, side= side, vn_kV = vn_hv_kv)
    push!(Transformers, trafo)
  end

  function addProsumer!(;vn_kv::Float64, busIdx::Int, type::String, p::Union{Nothing,Float64} = nothing, q::Union{Nothing,Float64} = nothing, referencePri::Union{Nothing,Integer} = nothing, vm_pu::Union{Nothing,Float64} = nothing, vm_degree::Union{Nothing,Float64} = nothing)
    idProsSum += 1
    isAPUNode = false
    if !isnothing(vm_pu)
      isAPUNode = true
    end  
    proTy = toProSumptionType(type)
    prosumer = ProSumer(vn_kv = vn_kv, busIdx = busIdx, oID = idProsSum, type = proTy, p = p, q = q, referencePri = referencePri, vm_pu = vm_pu, vm_degree = vm_degree, isAPUNode = isAPUNode)
    push!(proSumersVec, prosumer)
  end

  # A R E A 1
  # Bus 1, 220kV  
  addBus!(busName = "B1", busType = "PQ", vn_kV = 220.0, vm = 1.0, va = 0.0)
  # Bus 9, 22kV  
  slackBusIdx = 9
  addBus!(busName = "B9", busType = "Slack", vn_kV =22.0, vm = 1.03, va = 0.0)
  # Bus 7, 22kV  
  addBus!(busName = "B7", busType = "PQ", vn_kV = 22.0,  vm = 1.0, va = 0.0)
  # Bus 2, 220kV  
  addBus!(busName = "B2", busType = "PQ", vn_kV = 220.0,  pLoad = 250.0, qLoad = 200.0, vm = 1.0, va = 0.0)
  # Bus 10, 22kV  
  addBus!(busName = "B10", busType = "PV", vn_kV = 22.0,  pGen = 500.0, vm = 1.03, va = 0.0)
  
  # A R E A 2
  # B6a, 220kV  
  addBus!(busName = "B6a", busType = "PQ", vn_kV = 220.0, pLoad = 435.0, qLoad = 296.0, pShunt = 0.0, qShunt = 180.0, vm = 1.0, va = 0.0)
  # B6b, 220kV   
  addBus!(busName = "B6b", busType = "PQ", vn_kV = 220.0, vm = 1.0, va = 0.0)  
  # Bus 12, 22kV  
  addBus!(busName = "B12", busType = "PV", vn_kV = 22.0, pGen = 300.0, vm = 1.03, va = 0.0)
  
  # A R E A 3
  # Bus 5, 220kV  
  addBus!(busName = "B5", busType = "PQ", vn_kV = 220.0, vm = 1.0, va = 0.0, pLoad = 103.0, qLoad = 62.0, pShunt = 0.0, qShunt = 80.0)
  # Bus 4, 220kV  
  addBus!(busName = "B4", busType = "PQ", vn_kV = 220.0, vm = 1.0, va = 0.0, qShunt = 160.0)
  # Bus 3, 220kV
  addBus!(busName = "B3", busType = "PQ", vn_kV = 220.0, vm = 1.0, va = 0.0, pLoad = 103.0, qLoad = 62.0, isAux = false)
  # Bus 8, 380kV  
  addBus!(busName = "B8", busType = "PQ", vn_kV = 380.0, vm = 1.0, va = 0.0)
  # Bus 11, 22kV  
  addBus!(busName = "B11", busType = "PV", vn_kV = 22.0, pGen = 200.0, vm = 1.03, va = 0.0)
  
  ##### Shunt  
  addShunt!(busIdx = 6, pShunt = 0.0, qShunt = 180.0, vn = 220.0)
  addShunt!(busIdx = 4, pShunt = 0.0, qShunt = 160.0, vn = 220.0)
  addShunt!(busIdx = 5, pShunt = 0.0, qShunt = 80.0, vn = 220.0)

  ###### ACLines   
  # 220kV-Lines
  # (bus1, bus2, length_km=100)
  addACLine!(vn_kv = 220.0, from = getBusIdx("B1"), to = getBusIdx("B2"), length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)  
  # (bus1, bus6a, length_km=300)
  addACLine!(vn_kv = 220.0, from = getBusIdx("B1"), to = getBusIdx("B6a"), length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus2, bus5, length_km=300)
  addACLine!(vn_kv = 220.0, from = getBusIdx("B2"), to = getBusIdx("B5"), length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus3, bus4, length_km=100)
  addACLine!(vn_kv = 220.0, from = getBusIdx("B3"), to = getBusIdx("B4"), length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus3, bus4, length_km=100)
  addACLine!(vn_kv = 220.0, from = getBusIdx("B4"), to = getBusIdx("B4"), length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus4, bus5, length_km=300)
  addACLine!(vn_kv = 220.0, from = getBusIdx("B4"), to = getBusIdx("B5"), length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus4, bus6a, length_km=300)
  addACLine!(vn_kv = 220.0, from = getBusIdx("B4"), to = getBusIdx("B6a"), length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  #(bus6a, bus6b, length_km=length_km_6a_6b)
  addACLine!(vn_kv = 220.0, from = getBusIdx("B6a"), to = getBusIdx("B6b"), length = lLine_6a6b, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  
  # 380kV-Lines
  # (bus7, bus8, length_km=600)
  addACLine!(vn_kv = 380.0, from = getBusIdx("B7"), to = getBusIdx("B8"), length = 600.0, r = 0.0328, x = 0.312, c_nf_per_km = 11.5, tanδ = 0.0)
  
  #### Trafos
  # 'Trafo 1-7' (bus7, bus1, sn_mva=1000,  vn_hv_kv=380, vn_lv_kv=220, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0,   shift_degree=0.0)
  add2WTrafo!(from = getBusIdx("B7"), to=getBusIdx("B1"), sn_mva = 1000.0, vn_hv_kv = 380.0, vn_lv_kv = 220.0, vk_percent = 13.0, vkr_percent = 0.28, pfe_kw = 20.0, i0_percent = 0.06)
  #Trafo 3-8, (bus8, bus3, sn_mva=1000,  vn_hv_kv=380, vn_lv_kv=220, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=0.0)
  add2WTrafo!(from = getBusIdx("B8"), to=getBusIdx("B3"), sn_mva = 1000.0, vn_hv_kv = 380.0, vn_lv_kv = 220.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  # Trafo 9-1 (bus1, bus9, sn_mva=1000,vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0,vk_percent=13.0, pfe_kw=0, i0_percent=0,shift_degree=330.0)
  add2WTrafo!(from = getBusIdx("B1"), to=getBusIdx("B9"), sn_mva = 1000.0, vn_hv_kv = 220.0, vn_lv_kv = 22.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  # Trafo 10-2 (bus2, bus10, sn_mva=1000, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0)
  add2WTrafo!(from = getBusIdx("B2"), to=getBusIdx("B10"), sn_mva = 1000.0, vn_hv_kv = 220.0, vn_lv_kv = 22.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  # Trafo 11-3 (bus3, bus11, sn_mva=1000, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0)
  add2WTrafo!(from = getBusIdx("B3"), to=getBusIdx("B11"), sn_mva = 1000.0, vn_hv_kv = 220.0, vn_lv_kv = 22.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  # 'Trafo 12-6b' (bus6b, bus12, sn_mva=500, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0)
  add2WTrafo!(from = getBusIdx("B6b"), to=getBusIdx("B12"), sn_mva = 500.0, vn_hv_kv = 220.0, vn_lv_kv = 22.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)

  # Loads
  # 'Load 2' (bus2, p_mw=285, q_mvar=200)  
  addProsumer!(vn_kv = 220.0, busIdx = getBusIdx("B2"), type = "ENERGYCONSUMER", p = 285.0, q = 200.0)
  # 'Load 3' (bus3, p_mw=103.0, q_mvar=62.0)
  addProsumer!(vn_kv = 220.0, busIdx = getBusIdx("B3"), type = "ENERGYCONSUMER", p = 103.0, q = 62.0)  
  # 'Load 4' (net_cigre_hv, bus4, p_mw=326.0, q_mvar=244.0, name='Load 4')
  addProsumer!(vn_kv = 220.0, busIdx = getBusIdx("B4"), type = "ENERGYCONSUMER", p = 326.0, q = 244.0)
  # 'Load 5' (net_cigre_hv, bus5, p_mw=103, q_mvar=62, name='Load 5')	
  addProsumer!(vn_kv = 220.0, busIdx = getBusIdx("B5"), type = "ENERGYCONSUMER", p = 103.0, q = 62.0)
  # 'Load 6a' (net_cigre_hv, bus6a, p_mw=435, q_mvar=296, name='Load 6a')
  addProsumer!(vn_kv = 220.0, busIdx = getBusIdx("B6a"), type = "ENERGYCONSUMER", p = 435.0, q = 296.0)

  # Generators
  # 'Generator 9' (bus9, vm_pu=1.03, va_degree=0)
  addProsumer!(vn_kv = 22.0, busIdx = getBusIdx("B9"), type = "EXTERNALNETWORKINJECTION", vm_pu = 1.03, referencePri = 1)
  # 'Generator 10' (bus10, vm_pu=1.03, p_mw=500)
  addProsumer!(vn_kv = 22.0, busIdx = getBusIdx("B10"), type = "SYNCHRONOUSMACHINE", p = 500.0)
  # 'Generator 11' (bus11, vm_pu=1.03, p_mw=200)
  addProsumer!(vn_kv = 22.0, busIdx = getBusIdx("B11"), type = "SYNCHRONOUSMACHINE", p = 200.0)
  # 'Generator 12' (bus12, vm_pu=1.03, p_mw=300)
  addProsumer!(vn_kv = 22.0, busIdx = getBusIdx("B12"), type = "SYNCHRONOUSMACHINE", p = 300.0)

  sort!(nodeVec, lt = busComparison)
  #### reporting
  if verbose
    reports()
  end

  myNet = Net(netName, Sbase_MVA, slackBusIdx, nodeVec, ACLines, Transformers, branchVec, proSumersVec, shuntVec)
  return myNet
end

test_acpflow(0)
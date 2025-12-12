# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 20.6.2023
# CIGRE HV network
# Source of this network can be found here (Task Force C6.04.02 ): https://www.researchgate.net/publication/271963972_TF_C60402_TB_575_--_Benchmark_Systems_for_Network_Integration_of_Renewable_and_Distributed_Energy_Resources

function test_acpflow(verbose::Int = 0; lLine_6a6b::Float64 = 0.01, damp::Float64 = 1.0, method::Symbol = :rectangular, opt_sparse = true)::Bool
  net = createCIGRE(lLine_6a6b)
  tol = 1e-6
  maxIte = 25
  print_results = (verbose > 0)
  result = true
  etime = @elapsed begin
    ite, erg = runpf!(net, maxIte, tol, verbose; method = method, damp = damp, opt_sparse = opt_sparse)
  end
  if erg != 0
    @warn "Power flow did not converge"
    result = false
  end

  if print_results
    V = buildVoltageVector(net)
    calcNetLosses!(net, V)
    distributeBusResults!(net)
    printACPFlowResults(net, etime, ite, tol; converged = (erg == 0), solver = method)
    printProsumerResults(net)
    printQLimitLog(net; sort_by = :bus)
  end

  return result
end

function test_NBI_MDO()
  myNet = createCIGRE()
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

  for b in myNet.branchVec
    tupple = (b.fromBus, b.toBus)
    if tupple ∉ branchTupleSet
      push!(branchTupleSet, tupple)
    end
  end

  A = getNBI(nodeNumberVec, branchTupleSet)
  order = mdoRCM(length(nodeNumberVec), branchTupleSet)
  if order == [2, 5, 8, 7, 13, 3, 12, 11, 1, 4, 6, 9, 10]
    result = result & true
  else
    @warn "expected MDO-Order: [2, 5, 8, 7, 13, 3, 12, 11, 1, 4, 6, 9, 10], result order: $(order)"
    result = false
  end

  return result
end

function testNetwork()::Bool
  myNet = createCIGRE()
  result, msg = validate!(net = myNet)
  return result
end

function getTestFilePathName()
  filename = "cigre.m"
  jpath = joinpath(pwd(), "data", "mpower", filename)
  return jpath
end

function testExportMatpower()
  myNet = createCIGRE()
  case = myNet.name
  writeMatpowerCasefile(myNet, getTestFilePathName())
  return true
end

function testImportMatpower()
  filename = "case5.m"
  path = joinpath(pwd(), "data", "mpower", filename)

  net = createNetFromMatPowerFile(path, false)
  if strip(net.name) != "case5"
    @warn "Failed to import network case5.m from file: $path"
    return false
  end

  if length(net.nodeVec) != 5
    @warn "Expected 5 nodes, found: $(length(net.nodeVec))"
    return false
  end
  if length(net.branchVec) != 5
    @warn "Expected 5 branches, found: $(length(net.branchVec))"
    return false
  end

  return true
end

function rmTestfiles()
  file = getTestFilePathName()
  if isfile(file)
    rm(file)
  end
  return true
end

function testISOBusses()
  Sbase_MVA = 1000.0
  netName = "isobus"
  net = Net(name = netName, baseMVA = Sbase_MVA)
  @debug "Creating $netName test network"
  addBus!(net = net, busName = "B1", busType = "PQ", vn_kV = 220.0)
  addBus!(net = net, busName = "B2", busType = "PV", vn_kV = 220.0)
  addBus!(net = net, busName = "B3", busType = "PQ", vn_kV = 220.0)
  addBus!(net = net, busName = "B4", busType = "PQ", vn_kV = 220.0)
  addBus!(net = net, busName = "B5", busType = "SLACK", vn_kV = 220.0)

  addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)
  addACLine!(net = net, fromBus = "B1", toBus = "B3", length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)
  addACLine!(net = net, fromBus = "B2", toBus = "B4", length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)
  addACLine!(net = net, fromBus = "B3", toBus = "B4", length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)
  addACLine!(net = net, fromBus = "B5", toBus = "B4", length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)
  addACLine!(net = net, fromBus = "B5", toBus = "B4", length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0, ratedS = 250.0, status = 1)

  addShunt!(net = net, busName = "B3", pShunt = 0.0, qShunt = 180.0)
  addShunt!(net = net, busName = "B2", pShunt = 0.0, qShunt = 180.0)
  addShunt!(net = net, busName = "B4", pShunt = 0.0, qShunt = 180.0)

  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 285.0, q = 200.0)
  addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 103.0, q = 62.0)

  addProsumer!(net = net, busName = "B5", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.03, va_deg = 0.0, referencePri = "B5")
  addProsumer!(net = net, busName = "B2", type = "SYNCHRONOUSMACHINE", p = 600.0, vm_pu = 1.03, va_deg = 0.0)

  result, msg = validate!(net = net, log = true)
  if !result
    @warn msg
    return false
  end
  tol = 1e-6
  maxIte = 10
  verbose = 0       # 0: no output, 1: iteration norm, 2: + Y-Bus, 3: + Jacobian, 4: + Power Flow

  etime = @elapsed begin
    ite, erg = runpf!(net, maxIte, tol, verbose)
  end
  if erg != 0
    @warn "Power flow did not converge"
    printACPFlowResults(net, etime, ite, tol; converged = false)
    return false
  end

  brVec = getNetBranchNumberVec(net = net, fromBus = "B1", toBus = "B2")
  setNetBranchStatus!(net = net, branchNr = brVec[1], status = 0)
  brVec = getNetBranchNumberVec(net = net, fromBus = "B1", toBus = "B3")
  setNetBranchStatus!(net = net, branchNr = brVec[1], status = 0)

  if length(net.isoNodes) != 1
    @warn "Expected 1 isolated node, found: $(length(net.isoNodes))"
    return false
  end
  if net.isoNodes[1] != 1
    @warn "Expected isolated node 1, found: $(net.isoNodes[1])"
    return false
  end

  return true
end

function createCIGRE(lLine_6a6b = 0.01)::Net
  Sbase_MVA = 1000.0
  netName = "cigre"
  net = Net(name = netName, baseMVA = Sbase_MVA)
  @debug "Creating $netName test network"
  # A R E A 1
  # Bus 1, 220kV  
  addBus!(net = net, busName = "B1", busType = "PQ", vn_kV = 220.0)
  # Bus 9, 22kV  
  addBus!(net = net, busName = "B9", busType = "Slack", vn_kV = 22.0)
  # Bus 7, 22kV  
  addBus!(net = net, busName = "B7", busType = "PQ", vn_kV = 380.0)
  # Bus 2, 220kV  
  addBus!(net = net, busName = "B2", busType = "PQ", vn_kV = 220.0)
  # Bus 10, 22kV  
  addBus!(net = net, busName = "B10", busType = "PV", vn_kV = 22.0)
  # A R E A 2
  # B6a, 220kV  
  addBus!(net = net, busName = "B6a", busType = "PQ", vn_kV = 220.0)
  # B6b, 220kV   
  addBus!(net = net, busName = "B6b", busType = "PQ", vn_kV = 220.0)
  # Bus 12, 22kV  
  addBus!(net = net, busName = "B12", busType = "PV", vn_kV = 22.0)
  # A R E A 3
  # Bus 5, 220kV  
  addBus!(net = net, busName = "B5", busType = "PQ", vn_kV = 220.0)
  # Bus 4, 220kV  
  addBus!(net = net, busName = "B4", busType = "PQ", vn_kV = 220.0)
  # Bus 3, 220kV
  addBus!(net = net, busName = "B3", busType = "PQ", vn_kV = 220.0)
  # Bus 8, 380kV  
  addBus!(net = net, busName = "B8", busType = "PQ", vn_kV = 380.0)
  # Bus 11, 22kV  
  addBus!(net = net, busName = "B11", busType = "PV", vn_kV = 22.0)

  ##### Shunt    
  addShunt!(net = net, busName = "B6a", pShunt = 0.0, qShunt = 180.0)
  addShunt!(net = net, busName = "B4", pShunt = 0.0, qShunt = 160.0)
  addShunt!(net = net, busName = "B5", pShunt = 0.0, qShunt = 80.0)

  ###### ACLines   
  # 220kV-Lines
  # (bus1, bus2, length_km=100)
  addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus1, bus6a, length_km=300)
  addACLine!(net = net, fromBus = "B1", toBus = "B6a", length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus2, bus5, length_km=300)
  addACLine!(net = net, fromBus = "B2", toBus = "B5", length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus3, bus4, length_km=100)
  addACLine!(net = net, fromBus = "B3", toBus = "B4", length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus3, bus4, length_km=100)
  addACLine!(net = net, fromBus = "B3", toBus = "B4", length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus4, bus5, length_km=300)
  addACLine!(net = net, fromBus = "B4", toBus = "B5", length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus4, bus6a, length_km=300)
  addACLine!(net = net, fromBus = "B4", toBus = "B6a", length = 300.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  #(bus6a, bus6b, length_km=length_km_6a_6b)
  addACLine!(net = net, fromBus = "B6a", toBus = "B6b", length = lLine_6a6b, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)

  # 380kV-Lines
  # (bus7, bus8, length_km=600)
  addACLine!(net = net, fromBus = "B7", toBus = "B8", length = 600.0, r = 0.0328, x = 0.312, c_nf_per_km = 11.5, tanδ = 0.0)

  #### Trafos  
  # 'Trafo 1-7' (bus7, bus1, sn_mva=1000,  vn_hv_kv=380, vn_lv_kv=220, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0,   shift_degree=0.0)
  add2WTrafo!(net = net, fromBus = "B7", toBus = "B1", sn_mva = 1000.0, vk_percent = 13.0, vkr_percent = 0.28, pfe_kw = 20.0, i0_percent = 0.06)
  #Trafo 3-8, (bus8, bus3, sn_mva=1000,  vn_hv_kv=380, vn_lv_kv=220, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=0.0)
  add2WTrafo!(net = net, fromBus = "B8", toBus = "B3", sn_mva = 1000.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  # Trafo 9-1 (bus1, bus9, sn_mva=1000,vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0,vk_percent=13.0, pfe_kw=0, i0_percent=0,shift_degree=330.0)
  add2WTrafo!(net = net, fromBus = "B1", toBus = "B9", sn_mva = 1000.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  # Trafo 10-2 (bus2, bus10, sn_mva=1000, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0)
  add2WTrafo!(net = net, fromBus = "B2", toBus = "B10", sn_mva = 1000.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  # Trafo 11-3 (bus3, bus11, sn_mva=1000, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0)
  add2WTrafo!(net = net, fromBus = "B3", toBus = "B11", sn_mva = 1000.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  # 'Trafo 12-6b' (bus6b, bus12, sn_mva=500, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0)
  add2WTrafo!(net = net, fromBus = "B6b", toBus = "B12", sn_mva = 500.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)

  # Loads
  # 'Load 2' (bus2, p_mw=285, q_mvar=200)  
  # addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 285.0, q = 200.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 285.0, q = 200.0)
  # 'Load 3' (bus3, p_mw=103.0, q_mvar=62.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 80.0, q = 30.0)
  # 'Load 4' (net_cigre_hv, bus4, p_mw=326.0, q_mvar=244.0, name='Load 4')
  addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 320.0, q = 200.0)
  # 'Load 5' (net_cigre_hv, bus5, p_mw=103, q_mvar=62, name='Load 5')	
  addProsumer!(net = net, busName = "B5", type = "ENERGYCONSUMER", p = 80.0, q = 30.0)
  # 'Load 6a' (net_cigre_hv, bus6a, p_mw=435, q_mvar=296, name='Load 6a')
  addProsumer!(net = net, busName = "B6a", type = "ENERGYCONSUMER", p = 435.0, q = 296.0)

  # Generators
  # 'Generator 9' (bus9, vm_pu=1.03, va_degree=0)
  addProsumer!(net = net, busName = "B9", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.03, va_deg = 0.0, referencePri = "B9")
  # 'Generator 10' (bus10, vm_pu=1.03, p_mw=500)
  addProsumer!(net = net, busName = "B10", type = "SYNCHRONOUSMACHINE", p = 500.0, vm_pu = 1.03, va_deg = 0.0)
  # 'Generator 11' (bus11, vm_pu=1.03, p_mw=200)
  addProsumer!(net = net, busName = "B11", type = "SYNCHRONOUSMACHINE", p = 200.0, vm_pu = 1.03, va_deg = 0.0)
  # 'Generator 12' (bus12, vm_pu=1.03, p_mw=300)
  addProsumer!(net = net, busName = "B12", type = "SYNCHRONOUSMACHINE", p = 300.0, vm_pu = 1.03, va_deg = 0.0)

  return net
end

function createTest5BusNet(; cooldown = 0, hyst_pu = 0.0, qlim_min = nothing, qlim_max = nothing, pq_only = false, mul_gens::Bool = true)::Net
  Sbase_MVA = 100.0
  netName = "test5bus"
  r = 0.05
  x = 0.5
  c_nf_per_km = 10.0
  tanδ = 0.0

  if pq_only
    busTypeBus3 = "PQ"
  else
    busTypeBus3 = "PV"
  end

  @debug "Creating $netName test network with qlim_min=$qlim_min, qlim_max=$qlim_max, pq_only=$pq_only"

  Bus5Net = Net(name = netName, baseMVA = Sbase_MVA, cooldown_iters = cooldown, q_hyst_pu = hyst_pu)
  addBus!(net = Bus5Net, busName = "B1", busType = "Slack", vn_kV = 110.0)
  addBus!(net = Bus5Net, busName = "B2", busType = "PQ", vn_kV = 110.0)
  addBus!(net = Bus5Net, busName = "B3", busType = busTypeBus3, vn_kV = 110.0)
  addBus!(net = Bus5Net, busName = "B4", busType = "PQ", vn_kV = 110.0)
  addBus!(net = Bus5Net, busName = "B5", busType = busTypeBus3, vn_kV = 110.0)

  addACLine!(net = Bus5Net, fromBus = "B1", toBus = "B3", length = 20.0, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  addACLine!(net = Bus5Net, fromBus = "B3", toBus = "B5", length = 50.0, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  addACLine!(net = Bus5Net, fromBus = "B5", toBus = "B4", length = 100.0, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  addACLine!(net = Bus5Net, fromBus = "B4", toBus = "B1", length = 100.0, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  addACLine!(net = Bus5Net, fromBus = "B4", toBus = "B2", length = 100.0, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  addACLine!(net = Bus5Net, fromBus = "B2", toBus = "B1", length = 20.0, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)

  addProsumer!(net = Bus5Net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  if mul_gens
    addProsumer!(net = Bus5Net, busName = "B2", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
    addProsumer!(net = Bus5Net, busName = "B2", type = "ENERGYCONSUMER", p = 25.0, q = 7.0)
    addProsumer!(net = Bus5Net, busName = "B3", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 15.0, vm_pu = 1.03, va_deg = 0.0, qMax = 0.6*qlim_max, qMin = 0.6*qlim_min)
    addProsumer!(net = Bus5Net, busName = "B3", type = "SYNCHRONOUSMACHINE", p = 10.0, q = 10.0, vm_pu = 1.03, va_deg = 0.0, qMax = 0.4*qlim_max, qMin = 0.4*qlim_min)
    addProsumer!(net = Bus5Net, busName = "B4", type = "ENERGYCONSUMER", p = 25.0, q = 7.0)
    addProsumer!(net = Bus5Net, busName = "B4", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
    addProsumer!(net = Bus5Net, busName = "B5", type = "ENERGYCONSUMER", p = 10.0, q = 5.0)
    addProsumer!(net = Bus5Net, busName = "B5", type = "ENERGYCONSUMER", p = 15.0, q = 1.0)
    if !pq_only
      addProsumer!(net = Bus5Net, busName = "B5", type = "SYNCHRONOUSMACHINE", p = 0.0, vm_pu = 1.0, qMax = qlim_max, qMin = qlim_min)
    end
  else
    addProsumer!(net = Bus5Net, busName = "B2", type = "ENERGYCONSUMER", p = 50.0, q = 15.0)
    addProsumer!(net = Bus5Net, busName = "B3", type = "SYNCHRONOUSMACHINE", p = 30.0, q = 25.0, vm_pu = 1.0, va_deg = 0.0, qMax = qlim_max, qMin = qlim_min)
    addProsumer!(net = Bus5Net, busName = "B4", type = "ENERGYCONSUMER", p = 50.0, q = 15.0)
    addProsumer!(net = Bus5Net, busName = "B5", type = "ENERGYCONSUMER", p = 25.0, q = 10.0)
  end

  #addShunt!(net = Bus5Net, busName = "B5", pShunt = 0.0, qShunt = -5.7)

  return Bus5Net
end

function createTest3BusNet(; cooldown = 0, hyst_pu = 0.0, qlim_min = nothing, qlim_max = nothing)::Net
  # Simple 3-bus network
  #
  #  ASTADT        STATION1
  # <--|---------------|<--- Generator 
  #    |-------       |
  #            |      |
  #            --------|<---- EXTERNALNETWORKINJECTION
  #                 VERBUND  
  Sbase_MVA = 100.0
  netName = "test3bus"

  r = 0.0
  x = 0.4
  s = 25.0
  c_nf_per_km = 9.55
  tanδ = 0.0

  vm_pu_STATION1 = 1.027273
  vm_pu_VERBUND = 1.018182

  @debug "Creating $netName test network with qlim_min=$qlim_min, qlim_max=$qlim_max"

  Bus3Net = Net(name = netName, baseMVA = Sbase_MVA, cooldown_iters = cooldown, q_hyst_pu = hyst_pu)

  addBus!(net = Bus3Net, busName = "ASTADT", busType = "PQ", vn_kV = 110.0)
  addBus!(net = Bus3Net, busName = "STATION1", busType = "PV", vn_kV = 110.0)
  addBus!(net = Bus3Net, busName = "VERBUND", busType = "SLACK", vn_kV = 110.0)

  addACLine!(net = Bus3Net, fromBus = "ASTADT", toBus = "STATION1", length = s, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  addACLine!(net = Bus3Net, fromBus = "ASTADT", toBus = "VERBUND", length = s, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  addACLine!(net = Bus3Net, fromBus = "VERBUND", toBus = "STATION1", length = s, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)

  addProsumer!(net = Bus3Net, busName = "VERBUND", type = "EXTERNALNETWORKINJECTION", vm_pu = vm_pu_VERBUND, va_deg = 0.0, referencePri = "VERBUND")
  addProsumer!(net = Bus3Net, busName = "STATION1", type = "SYNCHRONOUSMACHINE", p = 70.0, q = 33.2, vm_pu = vm_pu_STATION1, qMax = qlim_max, qMin = qlim_min)
  addProsumer!(net = Bus3Net, busName = "ASTADT", type = "ENERGYCONSUMER", p = 100.0, q = 30.0)

  return Bus3Net
end

function createTest2BusNet(; cooldown = 0, hyst_pu = 0.0, qlim_min = nothing, qlim_max = nothing)::Net
  # Simple 2-bus network
  #   Bus1            Bus2
  #  ->|---------------|<--- Generator
  #   Slack           PV
  Sbase_MVA = 100.0
  netName = "test2bus"
  r = 0.05
  x = 0.5
  c_nf_per_km = 10.0
  tanδ = 0.0

  @debug "Creating $netName test network with qlim_min=$qlim_min, qlim_max=$qlim_max"

  Bus2Net = Net(name = netName, baseMVA = Sbase_MVA, cooldown_iters = cooldown, q_hyst_pu = hyst_pu)
  addBus!(net = Bus2Net, busName = "B1", busType = "Slack", vn_kV = 110.0)
  addBus!(net = Bus2Net, busName = "B2", busType = "PV", vn_kV = 110.0)

  addACLine!(net = Bus2Net, fromBus = "B1", toBus = "B2", length = 20.0, r = r, x = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)

  addProsumer!(net = Bus2Net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  addProsumer!(net = Bus2Net, busName = "B2", type = "SYNCHRONOUSMACHINE", p = 20.0, vm_pu = 1.05, va_deg = 0.0, qMax = qlim_max, qMin = qlim_min)

  return Bus2Net
end

function test_3BusNet(verbose::Int = 0, qlim::Float64 = 15.0, method::Symbol = :rectangular, opt_fd::Bool = false, opt_sparse::Bool = false)
  net = createTest3BusNet(cooldown = 2, hyst_pu = 0.01, qlim_min = -qlim, qlim_max = qlim)
  tol = 1e-12
  maxIte = 50
  print_results = (verbose > 0)
  result = true

  pv_names = ["STATION1"]
  etim = 0.0
  etim = @elapsed begin
    ite, erg = runpf!(net, maxIte, tol, verbose, method = method, opt_fd = opt_fd, opt_sparse = opt_sparse)
    if erg != 0
      @info "Full-system power flow did not converge"
      result = false
    end
  end

  #

  hit = pv_hit_q_limit(net, pv_names)
  calcNetLosses!(net)
  distributeBusResults!(net)
  if print_results
    printACPFlowResults(net, etim, ite, tol; converged = result, solver = method)
    printProsumerResults(net)
    printQLimitLog(net; sort_by = :bus)
  end

  tp = getTotalBusPower(net = net)
  tl = getTotalLosses(net = net)
  if verbose > 1
    @info "Total Power: P=$(tp[1]), Q=$(tp[2])"
    @info "Total Losses: P=$(tl[1]), Q=$(tl[2])"
  end

  @test all(isapprox.(tp, tl; atol = 1e-6))

  if qlim < 33.2
    return hit==true
  else
    return hit==false
  end
end

function test_5BusNet(verbose::Int = 0, qlim::Float64 = 20.0, method::Symbol = :rectangular, opt_fd::Bool = false, opt_sparse::Bool = false)
  net = createTest5BusNet(cooldown = 2, hyst_pu = 0.01, qlim_min = -qlim, qlim_max = qlim)
  tol = 1e-9
  maxIte = 50
  print_results = (verbose > 0)
  result = true

  pv_names = ["B3"]
  etim = 0.0
  etim = @elapsed begin
    ite, erg = runpf!(net, maxIte, tol, verbose, method = method, opt_fd = opt_fd, opt_sparse = opt_sparse)
    if erg != 0
      @info "Full-system power flow did not converge"
      result = false
    end
  end

  #

  hit = pv_hit_q_limit(net, pv_names)

  if print_results
    V = buildVoltageVector(net)
    calcNetLosses!(net, V)
    distributeBusResults!(net)
    printACPFlowResults(net, etim, ite, tol; converged = result, solver = method)
    printProsumerResults(net)
    printQLimitLog(net; sort_by = :bus)
  end

  return hit==true
end
export test_3BusNet, test_5BusNet

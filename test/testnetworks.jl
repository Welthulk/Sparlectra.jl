# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 20.6.2023
# CIGRE HV network (net_cigre_hv)
# Source of this network can be found here (Task Force C6.04.02 ): https://www.researchgate.net/publication/271963972_TF_C60402_TB_575_--_Benchmark_Systems_for_Network_Integration_of_Renewable_and_Distributed_Energy_Resources
using Sparlectra.ResDataTypes
using Sparlectra.SparqlQueryCGMES
using Sparlectra.SparlectraTools
using Sparlectra.SparlectraExport
using Sparlectra.SparlectraNet
using Sparlectra.SparlectraResult
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

  Sbase_MVA = 1000.0
  netName = "testNetWork"

  # Seite 1 = from, Seite 2 = to
  Seite1 = ResDataTypes.toSeitenTyp("S1")
  Seite2 = ResDataTypes.toSeitenTyp("S2")

  # length of coupling bus 6b-6a
  lLine_6a6b = 1.0

  # Netzwerk erstellen
  nodeVec = Vector{ResDataTypes.Node}()
  shuntVec = Vector{ResDataTypes.Shunt}()

  # A R E A 1 start
  # Busses: B9, B10, B1, B2, B7

  # Bus 1, 4 Terminals, 220kV
  # Seite 1, Line 1-2 (net_cigre_hv, bus1, bus2, length_km=100, std_type='Line220kV', name='Line 1-2')
  # Seite 1, Line 1-6a (net_cigre_hv, bus1, bus6a, length_km=300, std_type='Line220kV', name='Line 1-6a')
  # Seite 2, 'Trafo 1-7' (net_cigre_hv, bus7, bus1, sn_mva=1000,  vn_hv_kv=380, vn_lv_kv=220, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0,   shift_degree=0.0, name='Trafo 1-7')
  # Seite 2, 'Trafo 9-1' (net_cigre_hv, bus1, bus9, sn_mva=1000,vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0,vk_percent=13.0, pfe_kw=0, i0_percent=0,shift_degree=330.0, name='Trafo 9-1')

  tVec = ResDataTypes.Terminal[]
  Vn = 220.0

  # ACLINESEGMENT
  cName = "Line 1-2"
  cID = "L1-2_20230620"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)

  tL1_2_B1 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tL1_2_B1)
  # ACLINESEGMENT
  cName = "Line 1-6a"
  cID = "L1-6a_20230620"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)

  tL1_6a_B1 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tL1_6a_B1)
  # POWERTRANSFORMER
  cName = "Trafo 1-7"
  cID = "T1-7_20230620"
  cTrafo = ResDataTypes.toComponentTyp("POWERTRANSFORMER")
  c = ResDataTypes.Component(cID, cName, cTrafo, Vn)
  tT1_7_B1 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tT1_7_B1)
  # POWERTRANSFORMER
  cName = "Trafo 9-1"
  cID = "T9-1_20230620"
  cTrafo = ResDataTypes.toComponentTyp("POWERTRANSFORMER")
  c = ResDataTypes.Component(cID, cName, cTrafo, Vn)
  tT1_9_B1 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tT1_9_B1)

  # Create Node B1
  nName = "B1"
  nID = "B001_20230620"
  busIdx = 1
  nArea = 1
  nZone = 1

  vm = 1.0
  va = 0.0
  p_gen = 0.0
  q_gen = 0.0
  p_load = 0.0
  q_load = 0.0
  p_shunt = 0.0
  q_shunt = 0.0
  cNodeType = ResDataTypes.Busbarsection
  c = Component(nID, nName, cNodeType, Vn)
  nodeType = toNodeType("PQ")
  node = Node(c, tVec, busIdx, nodeType, nothing, Sbase_MVA, nZone, nArea, vm, va, p_load, q_load, p_shunt, q_shunt, p_gen, q_gen)
  push!(nodeVec, node)

  # Bus 9, 2 Terminals, 22kV, Area 1
  # Seite 1, 'Trafo 9-1' (net_cigre_hv, bus1, bus9, sn_mva=1000,vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0,vk_percent=13.0, pfe_kw=0, i0_percent=0,shift_degree=330.0, name='Trafo 9-1')
  # Seite 2, 'Generator 9' (net_cigre_hv, bus9, vm_pu=1.03, va_degree=0, name='Generator 9')
  Vn = 22.0
  tVec = ResDataTypes.Terminal[]
  # ExternalNetworkInjection
  cName = "Generator 9"
  cID = "E9_20230620"
  cExInj = ResDataTypes.toComponentTyp("ExternalNetworkInjection")
  c = ResDataTypes.Component(cID, cName, cExInj, Vn)
  tE1_B9 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tE1_B9)
  # POWERTRANSFORMER
  cName = "Trafo 9-1"
  cID = "T9-1_20230620"
  cTrafo = ResDataTypes.toComponentTyp("POWERTRANSFORMER")
  c = ResDataTypes.Component(cID, cName, cTrafo, Vn)
  tT9_1_B9 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tT9_1_B9)

  nName = "B9"
  nID = "B009_20230620"

  busIdx = 9
  nArea = 1
  nZone = 1
  slackBusIdx = busIdx

  vm = 1.03
  va = 0.0
  p_gen = 0.0
  q_gen = 0.0
  p_load = 0.0
  q_load = 0.0
  p_shunt = 0.0
  q_shunt = 0.0
  cNodeType = ResDataTypes.Busbarsection
  c = Component(nID, nName, cNodeType, Vn)
  nodeType = toNodeType("Slack")
  node = Node(c, tVec, busIdx, nodeType, nothing, 1000.0, nZone, nArea, vm, va, p_load, q_load, p_shunt, q_shunt, p_gen, q_gen)
  push!(nodeVec, node)

  # Bus 7, 2 Terminals, 22kV, Area 1
  # Seite 2, 'Trafo 1-7' (net_cigre_hv, bus1, bus7, sn_mva=1000,vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0,vk_percent=13.0, pfe_kw=0, i0_percent=0,shift_degree=330.0, name='Trafo 1-7')
  # Seite 1, 'Line 7-8' (net_cigre_hv, bus7, bus8, length_km=600, std_type='Line380kV', name='Line 7-8')                     
  tVec = ResDataTypes.Terminal[]
  # POWERTRANSFORMER
  Vn = 380.0
  cName = "Trafo 1-7"
  cID = "T1-7_20230620"
  cTrafo = ResDataTypes.toComponentTyp("POWERTRANSFORMER")
  c = ResDataTypes.Component(cID, cName, cTrafo, Vn)
  tT1_7_B7 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tT1_7_B7)

  # ACLINESEGMENT
  cName = "Line 7-8"
  cID = "L7-8_20230620"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)
  tL7_8_B7 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tL7_8_B7)

  nName = "B7"
  nID = "B007_20230620"
  busIdx = 7
  nArea = 1
  nZone = 1
  vm = 1.0
  va = 0.0
  p_gen = 0.0
  q_gen = 0.0
  p_load = 0.0
  q_load = 0.0
  p_shunt = 0.0
  q_shunt = 0.0
  cNodeType = ResDataTypes.Busbarsection
  c = Component(nID, nName, cNodeType, Vn)
  nodeType = toNodeType("PQ")
  node = Node(c, tVec, busIdx, nodeType, nothing, 1000.0, nZone, nArea, vm, va, p_load, q_load, p_shunt, q_shunt, p_gen, q_gen)
  push!(nodeVec, node)

  # Bus 2, 4 Terminals, 220kV, Area 1
  # Seite 2, 'Line 1-2' (net_cigre_hv, bus1, bus2, length_km=100, std_type='Line220kV', name='Line 1-2')
  # Seite 1, 'Line 2-5' (net_cigre_hv, bus2, bus5, length_km=300, std_type='Line220kV', name='Line 2-5')
  # Seite 1, 'Trafo 10-2' (net_cigre_hv, bus2, bus10, sn_mva=1000, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0, name='Trafo 10-2')
  # Seite 2, 'Load 2' (net_cigre_hv, bus2, p_mw=285, q_mvar=200, name='Load 2')
  Vn = 220.0
  tVec = ResDataTypes.Terminal[]
  # ACLINESEGMENT
  cName = "Line 1-2"
  cID = "L1-2_20230620"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)
  tL1_2_B2 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tL1_2_B2)

  # ACLINESEGMENT
  cName = "Line 2-5"
  cID = "L2-5_20230620"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)
  tL2_5_B2 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tL2_5_B2)

  # POWERTRANSFORMER
  cName = "Trafo 10-2"
  cID = "T10-2_20230620"
  cTrafo = ResDataTypes.toComponentTyp("POWERTRANSFORMER")
  c = ResDataTypes.Component(cID, cName, cTrafo, Vn)
  tT10_2_B2 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tT10_2_B2)

  # EnergyConsumer
  cName = "Load 2"
  cID = "Load2_20230620"
  cLoad = ResDataTypes.toComponentTyp("EnergyConsumer")
  c = ResDataTypes.Component(cID, cName, cLoad, Vn)
  tLoad2_B2 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tLoad2_B2)

  # Bus 2
  nName = "B2"
  nID = "B002_20230620"
  busIdx = 2
  nArea = 1
  nZone = 1

  vm = 1.0
  va = 0.0
  p_gen = 0.0
  q_gen = 0.0
  p_load = 285.0
  q_load = 200.0
  p_shunt = 0.0
  q_shunt = 0.0
  cNodeType = ResDataTypes.Busbarsection
  c = Component(nID, nName, cNodeType, Vn)
  nodeType = toNodeType("PQ")
  node = Node(c, tVec, busIdx, nodeType, nothing, 1000.0, nZone, nArea, vm, va, p_load, q_load, p_shunt, q_shunt, p_gen, q_gen)
  push!(nodeVec, node)

  # Bus 10, 2 Terminals, 22kV, Area 1
  # Seite 2, 'Trafo 10-2' (net_cigre_hv, bus2, bus10, sn_mva=1000, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0, name='Trafo 10-2')
  # Seite 2, 'Generator 10' (net_cigre_hv, bus10, vm_pu=1.03, p_mw=500, name='Generator 10')
  Vn = 22.0
  tVec = ResDataTypes.Terminal[]
  # POWERTRANSFORMER
  cName = "Trafo 10-2"
  cID = "T10-2_20230620"
  cTrafo = ResDataTypes.toComponentTyp("POWERTRANSFORMER")
  c = ResDataTypes.Component(cID, cName, cTrafo, Vn)
  tT10_2_B10 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tT10_2_B10)
  # Generator
  cName = "Generator 10"
  cID = "Gen10_20230620"
  cGen = ResDataTypes.toComponentTyp("Generator")
  c = ResDataTypes.Component(cID, cName, cGen, Vn)
  tGen10_B10 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tGen10_B10)

  # Bus 10
  nName = "B10"
  nID = "B010_20230620"
  busIdx = 10
  nArea = 1
  nZone = 1

  vm = 1.03
  va = 0.0
  p_gen = 500.0
  q_gen = 0.0
  p_load = 0.0
  q_load = 0.0
  p_shunt = 0.0
  q_shunt = 0.0
  cNodeType = ResDataTypes.Busbarsection
  c = Component(nID, nName, cNodeType, Vn)
  nodeType = toNodeType("PV")
  node = Node(c, tVec, busIdx, nodeType, nothing, 1000.0, nZone, nArea, vm, va, p_load, q_load, p_shunt, q_shunt, p_gen, q_gen)
  push!(nodeVec, node)
  # A R E A 1 finished

  # A R E A 2 start

  # B6a 220kV, 4 Terminals
  # Seite 2, 'Line 1-6a' (net_cigre_hv, bus1, bus6a, length_km=300, std_type='Line220kV', name='Line 1-6a')
  # Seite 1, 'Line 4-6a' (net_cigre_hv, bus4, bus6a, length_km=300, std_type='Line220kV', name='Line 4-6a')
  # Seite 1, 'Line 6a-6b' (net_cigre_hv, bus6a, bus6b, length_km=length_km_6a_6b, std_type='Line220kV', name='Line 6a-6b')
  # Seite 1, 'Shunt 6a' (net_cigre_hv, bus6a, p_mw=0.0, q_mvar=-180, name='Shunt 6a')
  # Seite 2, 'Load 6a' (net_cigre_hv, bus6a, p_mw=435, q_mvar=296, name='Load 6a')
  # Area 3 
  Vn = 220.0
  tVec = ResDataTypes.Terminal[]

  # ACLINESEGMENT
  cName = "Line 1-6a"
  cID = "L1-6a_20230620"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)
  tL1_6a_B6a = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tL1_6a_B6a)

  # ACLINESEGMENT
  cName = "Line 4-6a"
  cID = "L4-6a_20230620"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)
  tL4_6a_B6a = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tL4_6a_B6a)

  # ACLINESEGMENT
  cName = "Line 6a-6b"
  cID = "L6a-6b_20230621"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)
  tL6a_6b_B6a = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tL6a_6b_B6a)

  # EnergyConsumer
  cName = "Load 6a"
  cID = "Load6a_20230621"
  cLoad = ResDataTypes.toComponentTyp("EnergyConsumer")
  c = ResDataTypes.Component(cID, cName, cLoad, Vn)
  tLoad6a_B6a = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tLoad6a_B6a)

  # LinearShuntCompensator
  cName = "Shunt 6a"
  cID = "Shunt6a_20230621"
  cShunt = ResDataTypes.toComponentTyp("LinearShuntCompensator")
  c = ResDataTypes.Component(cID, cName, cShunt, Vn)
  tShunt6a_B6a = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tShunt6a_B6a)

  # Bus 6a
  nName = "B6a"
  nID = "B6a_20230621"
  busIdx = 6
  nArea = 1
  nZone = 3

  vm = 1.0
  va = 0.0
  p_gen = 500.0
  q_gen = 0.0
  p_load = 435.0
  q_load = 296.0
  p_shunt = 0.0
  q_shunt = 180.0
  cNodeType = ResDataTypes.Busbarsection
  c = Component(nID, nName, cNodeType, Vn)
  nodeType = toNodeType("PQ")
  node = Node(c, tVec, busIdx, nodeType, nothing, Sbase_MVA, nZone, nArea, vm, va, p_load, q_load, p_shunt, q_shunt, p_gen, q_gen)
  push!(nodeVec, node)

  # B6b 220kV, 2 Terminals, Area 3
  # Seite 1, 'Trafo 12-6b' (net_cigre_hv, bus6b, bus12, sn_mva=500, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0, name='Trafo 12-6b')
  # Seite 2, 'Line 6a-6b' (net_cigre_hv, bus6a, bus6b, length_km=length_km_6a_6b, std_type='Line220kV', name='Line 6a-6b')
  Vn = 220.0
  tVec = ResDataTypes.Terminal[]

  # ACLINESEGMENT
  cName = "Line 6a-6b"
  cID = "L6a-6b_20230621"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)
  tL6a_6b_B6b = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tL6a_6b_B6b)

  # POWERTRANSFORMER
  cName = "Trafo 12-6b"
  cID = "T12-6b_20230621"
  cTrafo = ResDataTypes.toComponentTyp("POWERTRANSFORMER")
  c = ResDataTypes.Component(cID, cName, cTrafo, Vn)
  tT12_6b_B6b = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tT12_6b_B6b)

  # Bus 6b
  nName = "B6b"
  nID = "B6b_20230621"
  busIdx = 13
  nArea = 1
  nZone = 3

  vm = 1.0
  va = 0.0
  p_gen = 0.0
  q_gen = 0.0
  p_load = 0.0
  q_load = 0.0
  p_shunt = 0.0
  q_shunt = 0.0
  cNodeType = ResDataTypes.Busbarsection
  c = Component(nID, nName, cNodeType, Vn)
  nodeType = toNodeType("PQ")
  node = Node(c, tVec, busIdx, nodeType, nothing, Sbase_MVA, nZone, nArea, vm, va, p_load, q_load, p_shunt, q_shunt, p_gen, q_gen)
  push!(nodeVec, node)

  # Bus 12, 22kV, 2 Terminals, Area 3
  # Seite2 'Trafo 12-6b' (net_cigre_hv, bus6b, bus12, sn_mva=500, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0, name='Trafo 12-6b')
  # Seite1 'Generator 12' (net_cigre_hv, bus12, vm_pu=1.03, p_mw=300, name='Generator 12')
  Vn = 22.0
  tVec = ResDataTypes.Terminal[]
  # POWERTRANSFORMER
  cName = "Trafo 12-6b"
  cID = "T12-6b_20230621"
  cTrafo = ResDataTypes.toComponentTyp("POWERTRANSFORMER")
  c = ResDataTypes.Component(cID, cName, cTrafo, Vn)
  tT12_6b_B12 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tT12_6b_B12)

  # Generator
  cName = "Generator 12"
  cID = "Gen12_20230621"
  cGen = ResDataTypes.toComponentTyp("Generator")
  c = ResDataTypes.Component(cID, cName, cGen, Vn)
  tGen12_B12 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tGen12_B12)

  # Bus 12
  nName = "B12"
  nID = "B12_20230621"
  busIdx = 12
  nArea = 1
  nZone = 3

  vm = 1.03
  va = 0.0
  p_gen = 300.0
  q_gen = 0.0
  p_load = 0.0
  q_load = 0.0
  p_shunt = 0.0
  q_shunt = 0.0
  cNodeType = ResDataTypes.Busbarsection
  c = Component(nID, nName, cNodeType, Vn)
  nodeType = toNodeType("PV")
  node = Node(c, tVec, busIdx, nodeType, nothing, Sbase_MVA, nZone, nArea, vm, va, p_load, q_load, p_shunt, q_shunt, p_gen, q_gen)
  push!(nodeVec, node)

  # A R E A 2 finished

  # A R E A 3 start

  # Bus 5, 220kV, 4 Terminals
  # Seite2 'Line 2-5' (net_cigre_hv, bus2, bus5, length_km=300, std_type='Line220kV', name='Line 2-5')
  # Seite1 'Line 4-5' (net_cigre_hv, bus4, bus5, length_km=300,  std_type='Line220kV', name='Line 4-5')
  # Seite2 'Load 5' (net_cigre_hv, bus5, p_mw=103, q_mvar=62, name='Load 5')	
  # Seite1 'Shunt 5' (net_cigre_hv, bus5, p_mw=0.0, q_mvar=-80, name='Shunt 5')
  Vn = 220.0
  tVec = ResDataTypes.Terminal[]
  # ACLINESEGMENT
  cName = "Line 2-5"
  cID = "L2-5_20230620"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)
  tL2_5_B5 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tL2_5_B5)
  # ACLINESEGMENT
  cName = "Line 4-5"
  cID = "L4-5_20230621"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)
  tL4_5_B5 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tL4_5_B5)
  # Load
  cName = "Load 5"
  cID = "Load5_20230621"
  cLoad = ResDataTypes.toComponentTyp("Load")
  c = ResDataTypes.Component(cID, cName, cLoad, Vn)
  tLoad5_B5 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tLoad5_B5)
  # Shunt
  cName = "Shunt 5"
  cID = "Shunt5_20230621"
  cShunt = ResDataTypes.toComponentTyp("Shunt")
  c = ResDataTypes.Component(cID, cName, cShunt, Vn)
  tShunt5_B5 = ResDataTypes.Terminal(c, Seite1)

  push!(tVec, tShunt5_B5)

  # Bus 5
  nName = "B5"
  nID = "B5_20230621"
  busIdx = 5
  nArea = 1
  nZone = 2

  vm = 1.0
  va = 0.0
  p_gen = 0.0
  q_gen = 0.0
  p_load = 103.0
  q_load = 62.0
  p_shunt = 0.0
  q_shunt = 80.0
  cNodeType = ResDataTypes.Busbarsection
  c = Component(nID, nName, cNodeType, Vn)
  nodeType = toNodeType("PQ")
  node = Node(c, tVec, busIdx, nodeType, nothing, Sbase_MVA, nZone, nArea, vm, va, p_load, q_load, p_shunt, q_shunt, p_gen, q_gen)
  push!(nodeVec, node)

  # Bus 4, 220kV, 6 Terminals
  # Seite 2, 'Line 4-5' (net_cigre_hv, bus4, bus5, length_km=300,  std_type='Line220kV', name='Line 4-5')
  # Seite 2, 'Line 4-6a' (net_cigre_hv, bus4, bus6a, length_km=300, std_type='Line220kV', name='Line 4-6a')
  # Seite 1, 'Line 3-4' (net_cigre_hv, bus3, bus4, length_km=100,  std_type='Line220kV', name='Line 3-4')
  # Seite 1, 'Line 3-4_2' (net_cigre_hv, bus3, bus4, length_km=100,  std_type='Line220kV', name='Line 3-4_2') parallele Leitung zu 'Line 3-4'
  # Seite 2, 'Load 4' (net_cigre_hv, bus4, p_mw=0.0, q_mvar=0.0, name='Load 4')
  # Seite 1, 'Shunt 4' (net_cigre_hv, bus4, p_mw=0.0, q_mvar=-160, name='Shunt 4')
  Vn = 220.0
  tVec = ResDataTypes.Terminal[]
  # ACLINESEGMENT
  cName = "Line 4-5"
  cID = "L4-5_20230621"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)
  tL4_5_B4 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tL4_5_B4)
  # ACLINESEGMENT
  cName = "Line 4-6a"
  cID = "L4-6a_20230620"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)
  tL4_6a_B4 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tL4_6a_B4)

  # ACLINESEGMENT
  cName = "Line 3-4"
  cID = "L3-4_20230621"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)
  tL3_4_B4 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tL3_4_B4)
  # ACLINESEGMENT
  cName = "Line 3-4_2"
  cID = "L3-4_2_20230621"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)
  tL3_4_2_B4 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tL3_4_2_B4)
  # Load
  cName = "Load 4"
  cID = "Load4_20230621"
  cLoad = ResDataTypes.toComponentTyp("Load")
  c = ResDataTypes.Component(cID, cName, cLoad, Vn)
  tLoad4_B4 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tLoad4_B4)
  # Shunt
  cName = "Shunt 4"
  cID = "Shunt4_20230621"
  cShunt = ResDataTypes.toComponentTyp("Shunt")
  c = ResDataTypes.Component(cID, cName, cShunt, Vn)
  tShunt4_B4 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tShunt4_B4)
  # Bus 4
  nName = "B4"
  nID = "B4_20230621"
  busIdx = 4
  nArea = 1
  nZone = 2

  vm = 1.0
  va = 0.0
  p_gen = 0.0
  q_gen = 0.0
  p_load = 0.0
  q_load = 0.0
  p_shunt = 0.0
  q_shunt = 160.0
  cNodeType = ResDataTypes.Busbarsection
  c = Component(nID, nName, cNodeType, Vn)
  nodeType = toNodeType("PQ")
  node = Node(c, tVec, busIdx, nodeType, nothing, Sbase_MVA, nZone, nArea, vm, va, p_load, q_load, p_shunt, q_shunt, p_gen, q_gen)
  push!(nodeVec, node)

  # Bus 3, 220kV, 5 Terminals
  # Seite 2, 'Line 3-4' (net_cigre_hv, bus3, bus4, length_km=100,std_type='Line220kV', name='Line 3-4')
  # Seite 2, 'Line 3-4_2' (net_cigre_hv, bus3, bus4, length_km=100,std_type='Line220kV', name='Line 3-4_2') parallele Leitung zu 'Line 3-4'
  # Seite 2, 'Trafo 3-8', (net_cigre_hv, bus8, bus3, sn_mva=1000,  vn_hv_kv=380, vn_lv_kv=220, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=0.0, name='Trafo 3-8')
  # Seite 1, 'Trafo 11-3', (net_cigre_hv, bus3, bus11, sn_mva=1000, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0, name='Trafo 11-3')
  # Seite 1, 'Load 3' (net_cigre_hv, bus3, p_mw=103.0, q_mvar=62.0, name='Load 3')
  Vn = 220.0
  tVec = ResDataTypes.Terminal[]
  # ACLINESEGMENT
  cName = "Line 3-4"
  cID = "L3-4_20230621"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)
  tL3_4_B3 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tL3_4_B3)
  # ACLINESEGMENT
  cName = "Line 3-4_2"
  cID = "L3-4_2_20230621"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)
  tL3_4_2_B3 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tL3_4_2_B3)
  # TRANSFORMER
  cName = "Trafo 3-8"
  cID = "T3-8_20230620"
  cTrafo = ResDataTypes.toComponentTyp("TRANSFORMER")
  c = ResDataTypes.Component(cID, cName, cTrafo, Vn)
  tT3_8_B3 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tT3_8_B3)
  # TRANSFORMER
  cName = "Trafo 11-3"
  cID = "T11-3_20230621"
  cTrafo = ResDataTypes.toComponentTyp("TRANSFORMER")
  c = ResDataTypes.Component(cID, cName, cTrafo, Vn)
  tT3_11_B3 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tT3_11_B3)
  # Load
  cName = "Load 3"
  cID = "Load3_20230621"
  cLoad = ResDataTypes.toComponentTyp("Load")
  c = ResDataTypes.Component(cID, cName, cLoad, Vn)
  tLoad3_B3 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tLoad3_B3)
  # Bus 3
  nName = "B3"
  nID = "B3_20230621"
  busIdx = 3
  nArea = 1
  nZone = 2

  vm = 1.0
  va = 0.0
  p_gen = 0.0
  q_gen = 0.0
  p_load = 103.0
  q_load = 62.0
  p_shunt = 0.0
  q_shunt = 0.0
  cNodeType = ResDataTypes.Busbarsection
  c = Component(nID, nName, cNodeType, Vn)
  nodeType = toNodeType("PQ")
  node = Node(c, tVec, busIdx, nodeType, nothing, Sbase_MVA, nZone, nArea, vm, va, p_load, q_load, p_shunt, q_shunt, p_gen, q_gen)
  push!(nodeVec, node)

  # Bus 8, 380kV, 2 Terminals, Area 2
  # Seite 2, 'Line 7-8' (net_cigre_hv, bus7, bus8, length_km=600, std_type='Line380kV', name='Line 7-8')
  # Seite 1, 'Trafo 3-8' (net_cigre_hv, bus8, bus3, sn_mva=1000,  vn_hv_kv=380, vn_lv_kv=220, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=0.0, name='Trafo 3-8')

  Vn = 380.0
  tVec = ResDataTypes.Terminal[]
  # ACLINESEGMENT
  cName = "Line 7-8"
  cID = "L7-8_20230620"
  cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
  c = ResDataTypes.Component(cID, cName, cLine, Vn)
  tL7_8_B8 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tL7_8_B8)
  # TRANSFORMER
  cName = "Trafo 3-8"
  cID = "T3-8_20230620"
  cTrafo = ResDataTypes.toComponentTyp("TRANSFORMER")
  c = ResDataTypes.Component(cID, cName, cTrafo, Vn)
  tT3_8_B8 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tT3_8_B8)
  # Bus 8
  nName = "B8"
  nID = "B8_20230620"
  busIdx = 8
  nArea = 1
  nZone = 2

  vm = 1.0
  va = 0.0
  p_gen = 0.0
  q_gen = 0.0
  p_load = 0.0
  q_load = 0.0
  p_shunt = 0.0
  q_shunt = 0.0
  cNodeType = ResDataTypes.Busbarsection
  c = Component(nID, nName, cNodeType, Vn)
  nodeType = toNodeType("PQ")
  node = Node(c, tVec, busIdx, nodeType, nothing, Sbase_MVA, nZone, nArea, vm, va, p_load, q_load, p_shunt, q_shunt, p_gen, q_gen)
  push!(nodeVec, node)

  # Bus 11, 22kV, 2 Terminals, Area 2
  # Seite 2, 'Trafo 11-3' (net_cigre_hv, bus3, bus11, sn_mva=1000, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0, name='Trafo 11-3')
  # Seite 1, 'Generator 11' (net_cigre_hv, bus11, vm_pu=1.03, p_mw=200, name='Generator 11')

  Vn = 22.0
  tVec = ResDataTypes.Terminal[]
  # TRANSFORMER
  cName = "Trafo 11-3"
  cID = "T11-3_20230621"
  cTrafo = ResDataTypes.toComponentTyp("TRANSFORMER")
  c = ResDataTypes.Component(cID, cName, cTrafo, Vn)
  tT11_3_B11 = ResDataTypes.Terminal(c, Seite2)
  push!(tVec, tT11_3_B11)
  # GENERATOR
  cName = "Generator 11"
  cID = "G11_20230621"
  cGen = ResDataTypes.toComponentTyp("GENERATOR")
  c = ResDataTypes.Component(cID, cName, cGen, Vn)
  tG11_B11 = ResDataTypes.Terminal(c, Seite1)
  push!(tVec, tG11_B11)
  # Bus 11
  nName = "B11"
  nID = "B11_20230621"
  busIdx = 11
  nArea = 1
  nZone = 2

  vm = 1.03
  va = 0.0
  p_gen = 200.0
  q_gen = 0.0
  p_load = 0.0
  q_load = 0.0
  p_shunt = 0.0
  q_shunt = 0.0
  cNodeType = ResDataTypes.Busbarsection
  c = Component(nID, nName, cNodeType, Vn)
  nodeType = toNodeType("PV")
  node = Node(c, tVec, busIdx, nodeType, nothing, Sbase_MVA, nZone, nArea, vm, va, p_load, q_load, p_shunt, q_shunt, p_gen, q_gen)
  push!(nodeVec, node)

  ##### Shunt  
  sh_p = 0.0
  sh_q = 180.0
  sh_nID = "B6a_20230621"
  sh_BusIdx = 6
  in_service = 1
  vn = 220.0
  cCompT = ResDataTypes.LinearShuntCompensator
  cName = "Sh_B6a"
  cSh = Component(nID, cName, cCompT, vn)
  ratio = 1.0  # (Vn / vn_kv)^2
  y_pu = calcYShunt(p_shunt, q_shunt, ratio, Sbase_MVA)
  sh = ResDataTypes.Shunt(cSh, sh_nID, sh_BusIdx, sh_p, sh_q, y_pu, in_service)
  push!(shuntVec, sh)

  sh_p = 0.0
  sh_q = 160.0
  sh_nID = "B4_20230621"
  sh_BusIdx = 4
  in_service = 1
  sh_vn = 220.0
  cCompT = ResDataTypes.LinearShuntCompensator
  cName = "Sh_B4"
  cSh = Component(nID, cName, cCompT, sh_vn)
  ratio = 1.0  # (Vn / vn_kv)^2
  y_pu = calcYShunt(p_shunt, q_shunt, ratio, Sbase_MVA)
  sh = ResDataTypes.Shunt(cSh, sh_nID, sh_BusIdx, sh_p, sh_q, y_pu, in_service)
  push!(shuntVec, sh)

  sh_p = 0.0
  sh_q = 80.0
  sh_nID = "B5_20230621"
  sh_BusIdx = 5
  in_service = 1
  sh_vn = 220.0
  cCompT = ResDataTypes.LinearShuntCompensator
  cName = "Sh_B5"
  cSh = Component(nID, cName, cCompT, sh_vn)
  ratio = 1.0  # (Vn / vn_kv)^2
  y_pu = calcYShunt(p_shunt, q_shunt, ratio, Sbase_MVA)
  sh = ResDataTypes.Shunt(cSh, sh_nID, sh_BusIdx, sh_p, sh_q, y_pu, in_service)
  push!(shuntVec, sh)

  ###### ACLines 
  # Create Line Vector
  ACLines = SparqlQueryCGMES.ACLineVector()

  b = nothing
  g = nothing
  # 220kV-Lines
  Vn = 220.0
  r = 0.0653
  x = 0.398  
  c_nf_per_km = 9.08

  # (net_cigre_hv, bus1, bus2, length_km=100, std_type='Line220kV', name='Line 1-2')
  cName = "Line 1-2"
  cID = "L1-2_20230620"
  length = 100.0
  asec = ACLineSegment(cID, cName, Vn, length, r, x, b, g, c_nf_per_km)
  push!(ACLines, asec)

  # (net_cigre_hv, bus1, bus6a, length_km=300, std_type='Line220kV', name='Line 1-6a')
  cName = "Line 1-6a"
  cID = "L1-6a_20230620"
  length = 300.0
  asec = ACLineSegment(cID, cName, Vn, length, r, x, b, g, c_nf_per_km)
  push!(ACLines, asec)

  # (net_cigre_hv, bus2, bus5, length_km=300, std_type='Line220kV', name='Line 2-5')
  cName = "Line 2-5"
  cID = "L2-5_20230620"
  length = 300.0
  asec = ACLineSegment(cID, cName, Vn, length, r, x, b, g, c_nf_per_km)
  push!(ACLines, asec)

  # (net_cigre_hv, bus3, bus4, length_km=100, std_type='Line220kV', name='Line 3-4')
  cName = "Line 3-4"
  cID = "L3-4_20230621"
  length = 100.0
  asec = ACLineSegment(cID, cName, Vn, length, r, x, b, g, c_nf_per_km)
  push!(ACLines, asec)

  # (net_cigre_hv, bus3, bus4, length_km=100, std_type='Line220kV', name='Line 3-4_2')
  cName = "Line 3-4_2"
  cID = "L3-4_2_20230621"
  length = 100.0
  asec = ACLineSegment(cID, cName, Vn, length, r, x, b, g, c_nf_per_km)
  push!(ACLines, asec)

  # (net_cigre_hv, bus4, bus5, length_km=300, std_type='Line220kV', name='Line 4-5')
  cName = "Line 4-5"
  cID = "L4-5_20230621"
  length = 300.0
  asec = ACLineSegment(cID, cName, Vn, length, r, x, b, g, c_nf_per_km)
  push!(ACLines, asec)

  # (net_cigre_hv, bus4, bus6a, length_km=300, std_type='Line220kV', name='Line 4-6a')
  cName = "Line 4-6a"
  cID = "L4-6a_20230620"
  length = 300.0
  asec = ACLineSegment(cID, cName, Vn, length, r, x, b, g, c_nf_per_km)
  push!(ACLines, asec)

  #(net_cigre_hv, bus6a, bus6b, length_km=length_km_6a_6b, std_type='Line220kV', name='Line 6a-6b')
  cName = "Line 6a-6b"
  cID = "L6a-6b_20230621"
  length = lLine_6a6b #0.1 # Kupplung?
  asec = ACLineSegment(cID, cName, Vn, length, r, x, b, g, c_nf_per_km)
  push!(ACLines, asec)

  # 380kV-Lines
  Vn = 380.0
  r = 0.0328
  x = 0.312  
  c_nf_per_km = 11.5

  # (net_cigre_hv, bus7, bus8, length_km=600, std_type='Line380kV', name='Line 7-8')
  cName = "Line 7-8"
  cID = "L7-8_20230620"
  length = 600.0
  asec = ACLineSegment(cID, cName, Vn, length, r, x, b, g, c_nf_per_km)
  push!(ACLines, asec)

  #### Trafos
  Transformers = Vector{ResDataTypes.PowerTransformer}()
  sn = 1000.0
  vn_hv = 380.0
  vn_lv = 220.0
  vk_percent = 13.0
  vkr_percent = 0.28
  pfe_kw = 20.0
  i0 = 0.06
  shift_degree = 0.0

  regelungEin = false
  cName = "Trafo 1-7"
  cID = "T1-7_20230620"

  r_hv, x_hv, g_hv, b_hv = calcTrafoParams(sn_mva=sn, vn_hv_kv=vn_hv, vk_percent=vk_percent, vkr_percent=vkr_percent, pfe_kw=pfe_kw, i0_percent=i0)
  
  s1 = PowerTransformerWinding(Vn=vn_hv, r=r_hv, x=x_hv, isPu_RXGB=false, shift_degree=shift_degree)
  s2 = PowerTransformerWinding(Vn=vn_lv, r=0.0, x=0.0, isPu_RXGB=false, shift_degree=shift_degree)
  s3 = nothing
  cmp = ResDataTypes.Component(cID, cName, "POWERTRANSFORMER", vn_hv)
  Trafo = PowerTransformer(cmp, regelungEin, s1, s2, s3)
  push!(Transformers, Trafo)

  cName = "Trafo 3-8"
  cID = "T3-8_20230620"

  regelungEin = false
  cmp = ResDataTypes.Component(cID, cName, "POWERTRANSFORMER", vn_hv)
  Trafo = PowerTransformer(cmp, regelungEin, s1, s2, s3)
  push!(Transformers, Trafo)

  cName = "Trafo 9-1"
  cID = "T9-1_20230620"

  sn = 1000.0
  vn_hv = 220.0
  vn_lv = 22.0
  vk_percent = 18.0
  vkr_percent = 0.32
  shift_degree = 330.0
  r_hv, x_hv, g_hv, b_hv = calcTrafoParams(sn_mva=sn, vn_hv_kv=vn_hv, vk_percent=vk_percent, vkr_percent=vkr_percent, pfe_kw=pfe_kw, i0_percent=i0)
  s1 = PowerTransformerWinding(Vn=vn_hv, r=r_hv, x=x_hv, shift_degree=shift_degree, isPu_RXGB=false)
  s2 = PowerTransformerWinding(Vn=vn_lv, r=0.0, x=0.0,shift_degree=0.0, isPu_RXGB=false)
  s3 = nothing
  cmp = ResDataTypes.Component(cID, cName, "POWERTRANSFORMER", vn_hv)
  Trafo = PowerTransformer(cmp, regelungEin, s1, s2, s3)
  push!(Transformers, Trafo)

  cName = "Trafo 10-2"
  cID = "T10-2_20230620"

  cmp = ResDataTypes.Component(cID, cName, "POWERTRANSFORMER", vn_hv)
  Trafo = PowerTransformer(cmp, regelungEin, s1, s2, s3)
  push!(Transformers, Trafo)

  cName = "Trafo 11-3"
  cID = "T11-3_20230621"

  cmp = ResDataTypes.Component(cID, cName, "POWERTRANSFORMER", vn_hv)
  Trafo = PowerTransformer(cmp, regelungEin, s1, s2, s3)
  push!(Transformers, Trafo)

  cName = "Trafo 12-6b"
  cID = "T12-6b_20230621"

  sn = 500.0
  vn_hv = 220.0
  vn_lv = 22.0
  vk_percent = 16.2
  vkr_percent = 0.28
  shift_degree = 330.0
  r_hv, x_hv, g_hv, b_hv = calcTrafoParams(sn_mva=sn, vn_hv_kv=vn_hv, vk_percent=vk_percent, pk_kw=nothing, vkr_percent=vkr_percent, pfe_kw=pfe_kw, i0_percent=i0)
  s1 = PowerTransformerWinding(Vn=vn_hv, r=r_hv, x=x_hv, shift_degree=shift_degree, isPu_RXGB=false)
  s2 = PowerTransformerWinding(Vn=vn_lv, r=0.0, x=0.0,shift_degree=0.0, isPu_RXGB=false)
  s3 = nothing
  cmp = ResDataTypes.Component(cID, cName, "POWERTRANSFORMER", vn_hv)
  Trafo = PowerTransformer(cmp, regelungEin, s1, s2, s3)
  push!(Transformers, Trafo)

  #### PROSUMER
  #
  proSumersVec = Vector{ResDataTypes.ProSumer}()
  # Loads
  # EnergyConsumer
  Vn = 220.0
  cName = "Load 2"
  cID = "Load2_20230620"
  nID = "B002_20230620"
  ratedS = nothing
  ratedU = nothing
  qPercent = nothing
  p = 285.0
  q = 200.0
  maxP = nothing
  minP = nothing
  maxQ = nothing
  minQ = nothing
  ratedPowerFactor = nothing
  referencePri = nothing
  equipment = ResDataTypes.Component(cID, cName, "ENERGYCONSUMER", Vn)
  pRS = ProSumer(equipment, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri)
  push!(proSumersVec, pRS)

  Vn = 220.0
  cName = "Load 3"
  cID = "Load3_20230621"
  nID = "B3_20230621"
  ratedS = nothing
  ratedU = nothing
  qPercent = nothing
  p = 325.0
  q = 244.0
  maxP = nothing
  minP = nothing
  maxQ = nothing
  minQ = nothing
  ratedPowerFactor = nothing
  referencePri = nothing
  equipment = ResDataTypes.Component(cID, cName, "ENERGYCONSUMER", Vn)
  pRS = ProSumer(equipment, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri)
  push!(proSumersVec, pRS)

  Vn = 220.0
  cName = "Load 4"
  cID = "Load4_20230621"
  nID = "B4_20230621"
  ratedS = nothing
  ratedU = nothing
  qPercent = nothing
  p = 326.0
  q = 244.0
  maxP = nothing
  minP = nothing
  maxQ = nothing
  minQ = nothing
  ratedPowerFactor = nothing
  referencePri = nothing
  equipment = ResDataTypes.Component(cID, cName, "ENERGYCONSUMER", Vn)
  pRS = ProSumer(equipment, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri)
  push!(proSumersVec, pRS)

  Vn = 220.0
  cName = "Load 5"
  cID = "Load5_20230621"
  nID = "B5_20230621"
  ratedS = nothing
  ratedU = nothing
  qPercent = nothing
  p = 103.0
  q = 62.0
  maxP = nothing
  minP = nothing
  maxQ = nothing
  minQ = nothing
  ratedPowerFactor = nothing
  referencePri = nothing
  equipment = ResDataTypes.Component(cID, cName, "ENERGYCONSUMER", Vn)
  pRS = ProSumer(equipment, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri)
  push!(proSumersVec, pRS)

  Vn = 220.0
  cName = "Load 6a"
  cID = "Load6a_20230621"
  nID = "B6a_20230621"
  ratedS = nothing
  ratedU = nothing
  qPercent = nothing
  p = 435.0
  q = 296.0
  maxP = nothing
  minP = nothing
  maxQ = nothing
  minQ = nothing
  ratedPowerFactor = nothing
  referencePri = nothing
  equipment = ResDataTypes.Component(cID, cName, "ENERGYCONSUMER", Vn)
  pRS = ProSumer(equipment, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri)
  push!(proSumersVec, pRS)

  Vn = 22.0
  cName = "Generator 9"
  cID = "E9_20230620"
  nID = "B009_20230620"
  ratedS = nothing
  ratedU = nothing
  qPercent = nothing
  p = nothing
  q = nothing
  maxP = nothing
  minP = nothing
  maxQ = nothing
  minQ = nothing
  ratedPowerFactor = nothing
  referencePri = nothing
  equipment = ResDataTypes.Component(cID, cName, "EXTERNALNETWORKINJECTION", Vn)
  pRS = ProSumer(equipment, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri)
  push!(proSumersVec, pRS)

  Vn = 22.0
  cName = "Generator 10"
  cID = "Gen10_20230620"
  nID = "B010_20230620"

  ratedU = nothing
  qPercent = nothing
  p = 500.0
  q = nothing
  # invented values!
  ratedS = 500.0
  maxP = 800.0
  minP = 10.0
  maxQ = 200.0
  minQ = 0.0
  ratedPowerFactor = nothing
  referencePri = nothing
  equipment = ResDataTypes.Component(cID, cName, "SYNCHRONOUSMACHINE", Vn)
  pRS = ProSumer(equipment, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri)
  push!(proSumersVec, pRS)

  Vn = 22.0
  cName = "Generator 11"
  cID = "G11_20230621"
  nID = "B11_20230621"

  ratedU = nothing
  qPercent = nothing
  p = 200.0
  q = nothing
  # invented values!
  ratedS = 500.0
  maxP = 800.0
  minP = 10.0
  maxQ = 200.0
  minQ = 0.0
  ratedPowerFactor = nothing
  referencePri = nothing
  equipment = ResDataTypes.Component(cID, cName, "SYNCHRONOUSMACHINE", Vn)
  pRS = ProSumer(equipment, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri)
  push!(proSumersVec, pRS)

  Vn = 22.0
  cName = "Generator 12"
  cID = "Gen12_20230621"
  nID = "B12_20230621"

  ratedU = nothing
  qPercent = nothing
  p = 300.0
  q = nothing
  # invented values!
  ratedS = 500.0
  maxP = 800.0
  minP = 10.0
  maxQ = 200.0
  minQ = 0.0
  ratedPowerFactor = nothing
  referencePri = nothing
  equipment = ResDataTypes.Component(cID, cName, "SYNCHRONOUSMACHINE", Vn)
  pRS = ProSumer(equipment, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri)
  push!(proSumersVec, pRS)

  sort!(nodeVec, lt = busComparison)
  #### reporting
  if verbose
    reports()
  end

  branchVec = createBranchVectorFromNodeVector!(nodes=nodeVec, lines=ACLines, trafos=Transformers, Sbase_MVA=Sbase_MVA, prosumps=proSumersVec, shunts=shuntVec, log=true)
  

  myNet = Net(netName, Sbase_MVA, slackBusIdx, nodeVec, ACLines, Transformers, branchVec, proSumersVec, shuntVec)
  return myNet
end

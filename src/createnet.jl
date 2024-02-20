# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 28.6.2023
# include-file createnet.jl

roundUpToNearest100(number) = ceil(number / 100) * 100

"""
Purpose: correction of sequence number in node vector e.g. wrong numbering of "Seite"
"""
function fixSequenceNumberInNodeVec!(nodes::Vector{ResDataTypes.Node}, trafos::Vector{ResDataTypes.PowerTransformer})
  trafoDict = Dict{String,ResDataTypes.PowerTransformer}()
  trafoNameDict = Dict{String,String}()
  for trafo in trafos
    trafoDict[trafo.comp.cID] = trafo
    trafoNameDict[trafo.comp.cName] = trafo.comp.cID
  end

  for node in nodes
    terminals = node.terminals
    for term in terminals
      if term.comp.cTyp == ResDataTypes.Trafo
        # key = term.equipment.ID : TransformatorEnd (ID) <> TransformeID !
        nodename = node.comp.cName
        name = term.comp.cName
        id = nothing
        try
          id = trafoNameDict[name]
        catch
          @warn "key $name not found"
          continue
        end
        val = nothing
        try
          val = trafoDict[id] # trafo
        catch
          @warn "key $name not found"
          continue
        end

        if isnothing(val)
          @warn "node is nothing"
          continue
        end

        side = term.seite
        VnNode = node.comp.cVN

        VnCheck = nothing
        expectedSide = nothing
        VnSide1 = nothing
        VnSide2 = nothing
        VnSide3 = nothing

        if val.isBiWinder
          VnSide1 = val.side1.Vn
          VnSide2 = val.side2.Vn
          if side == ResDataTypes.Seite1
            VnCheck = val.side1.Vn
          elseif side == ResDataTypes.Seite2
            VnCheck = val.side2.Vn
          else
            @warn "VnNode: $VnNode: Seite <$side> for 2WT not expected!"
            continue
          end
        else
          VnSide1 = val.side1.Vn
          VnSide2 = val.side2.Vn
          VnSide3 = val.side3.Vn
          if side == ResDataTypes.Seite1
            VnCheck = val.side1.Vn
          elseif side == ResDataTypes.Seite2
            VnCheck = val.side2.Vn
          elseif side == ResDataTypes.Seite3
            VnCheck = val.side3.Vn
          else
            @warn "VnNode: $VnNode: Seite <$side> for 3WT not expected!"
            continue
          end
        end

        if VnCheck != VnNode
          # expectedSide -> Look for VnNode in the transformer-sides

          if VnNode == VnSide1
            expectedSide = ResDataTypes.Seite1
          elseif VnNode == VnSide2
            expectedSide = ResDataTypes.Seite2
          elseif VnNode == VnSide3
            expectedSide = ResDataTypes.Seite3
          end

          @info "fix wrong sequence number for node: $nodename trafo: $name from: $side to: $expectedSide"
          ResDataTypes.ChageSequnceNumber!(term, expectedSide)

        else
          continue
        end
      end # if term.equipment.Typ == ResDataTypes.Trafo 
    end # for term in terminals 
  end # for node in nodes
end # FixSequenceNumberInNodeVec!

function getBranchForLine(comp::ResDataTypes.AbstractComponent, fromBus::Int, toBus::Int, fromNodeID::String, toNodeID::String, line::ResDataTypes.ACLineSegment, ZBase::Float64)::ResDataTypes.Branch
  r = (line.r === nothing) ? 0.0 : line.r
  x = (line.x === nothing) ? 0.0 : line.x
  b = (line.b === nothing) ? 0.0 : line.b
  rpu = r / ZBase
  xpu = x / ZBase
  bpu = b * ZBase
  rateA = 0.0 # = unlimited
  rateB = 0.0 # = unlimited
  rateC = 0.0 # = unlimited
  ratio = 0.0 # line
  angle = 0.0 # line
  status = 1 # line
  angmin = -360.0 # line
  angmax = 360.0 # line  
  
  b = ResDataTypes.Branch(comp, fromBus, toBus, fromNodeID, toNodeID, rpu, xpu, bpu, 0.0, ratio, angle, status)
  return b
end

function getBranchFor2WT(comp::ResDataTypes.AbstractComponent, fromBus::Int, toBus::Int, fromNodeID::String, toNodeID::String, trafo::ResDataTypes.PowerTransformer, Sbase_MVA::Float64)::ResDataTypes.Branch
  HVSide = trafo.HVSideNumber
  tapSide = trafo.tapSideNumber
  shiftDegree = 0.0
  if HVSide == 1
    HV = trafo.side1
    LV = trafo.side2
    ratedUHV = (trafo.side1.ratedU === nothing) ? trafo.side1.Vn : trafo.side1.ratedU
    ratedULV = (trafo.side2.ratedU === nothing) ? trafo.side2.Vn : trafo.side2.ratedU
    ratedSHV = trafo.side1.ratedS
    ratedSLV = trafo.side2.ratedS

    if !isnothing(trafo.side1.shift_degree)
      shiftDegree = trafo.side1.shift_degree
    end
  else
    HV = trafo.side2
    LV = trafo.side1
    ratedUHV = (trafo.side2.ratedU === nothing) ? trafo.side2.Vn : trafo.side2.ratedU
    ratedULV = (trafo.side1.ratedU === nothing) ? trafo.side1.Vn : trafo.side1.ratedU
    ratedSHV = trafo.side2.ratedS
    ratedSLV = trafo.side1.ratedS
    if !isnothing(trafo.side2.shift_degree)
      shiftDegree = trafo.side2.shift_degree
    end
  end

  HVUbase = HV.Vn
  LVUbase = LV.Vn

  if tapSide == HVSide || tapSide == 0
    r = (HV.r === nothing) ? 0.0 : HV.r
    x = (HV.x === nothing) ? 0.0 : HV.x
    b = (HV.b === nothing) ? 0.0 : HV.b
    ZBase = HVUbase^2 / Sbase_MVA
  else
    r = (LV.r === nothing) ? 0.0 : LV.r
    x = (LV.x === nothing) ? 0.0 : LV.x
    b = (LV.b === nothing) ? 0.0 : LV.b
    ZBase = LVUbase^2 / Sbase_MVA
  end

  rpu = r / ZBase
  xpu = x / ZBase
  bpu = b * ZBase

  rateA = (ratedSHV === nothing) ? 0.0 : ratedSHV
  rateB = 0.0
  rateC = 0.0

  if trafo.isControlled
    if tapSide == 1
      taps = trafo.side1.taps
    elseif tapSide == 2
      taps = trafo.side2.taps
    else
      @assert "$trafo has more than two tap sides!"
    end

    step = taps.step
    nStep = taps.neutralStep
    stepVol = taps.voltageIncrement
    # ratio calculated in percent
    # inverse ratio if tap side is LV
    corr = (step - nStep) * stepVol
    corr = corr / 100.0
    corr = 1.0 + corr
    if trafo.tapSideNumber > 1
      corr = 1.0 / corr
    end
  else
    corr = 1.0
  end

  ratio = (HVUbase * ratedULV) / (LVUbase * ratedUHV)

  ratio = ratio * corr

  angle = shiftDegree
  status = 1

  angmin = -360.0
  angmax = 360.0

  b = ResDataTypes.Branch(comp, fromBus, toBus, fromNodeID, toNodeID, rpu, xpu, bpu, 0.0, ratio, angle, status)

  return b
end

function setBranchFor3WT!(hv_bus::Int, mv_bus::Int, lv_bus::Int, trafo::ResDataTypes.PowerTransformer, auxBus::ResDataTypes.Node, Sbase_MVA::Float64, bVec::Vector{ResDataTypes.Branch}, log::Bool = false)
  cID = trafo.comp.cID
  auxBuxCmp = auxBus.comp
  auxBusIdx = auxBus.busIdx
  ausBusID = auxBus.comp.cID

  vn_aux_kv = auxBus.comp.cVN
  vn_hv_kv = trafo.side1.Vn
  vn_mv_kv = trafo.side2.Vn
  vn_lv_kv = trafo.side3.Vn

  HV_kv = vn_hv_kv
  MV_kv = vn_mv_kv
  LV_kv = vn_lv_kv

  shift_mv_degree = trafo.side2.shift_degree
  shift_lv_degree = trafo.side3.shift_degree

  tapSide = 0
  tap_min = 0
  tap_max = 0
  tap_neutral = 0
  tap_step_percent = 0.0
  tap_pos = nothing
  neutralU = nothing
  regelungEin = false
  if trafo.isControlled
    regelungEin = true
    tapSide = trafo.tapSideNumber
    taps = nothing
    if tapSide == 1
      taps = trafo.side1.taps
    elseif tapSide == 2
      taps = trafo.side2.taps
    elseif tapSide == 3
      taps = trafo.side3.taps
    else
      @assert "$trafo has no taps!"
    end

    tap_pos = taps.step
    tap_min = taps.lowStep
    tap_max = taps.highStep
    tap_neutral = taps.neutralStep
    tap_step_percent = taps.voltageIncrement
    ntap = calcTapCorr(tap_pos, tap_neutral, tap_step_percent)

    neutralU = taps.neutralU

    if tapSide == 1
      neutralU = isnothing(neutralU) ? vn_hv_kv : neutralU
      vn_hv_kv = neutralU * ntap
    elseif tapSide == 2
      neutralU = isnothing(neutralU) ? vn_hv_kv : neutralU
      vn_mv_kv = neutralU * ntap
    elseif tapSide == 3
      neutralU = isnothing(neutralU) ? vn_hv_kv : neutralU
      vn_lv_kv = neutralU * ntap
    end
  end

  @debug "Controlled: $regelungEin, HV: $vn_hv_kv, Tap: $tap_pos, MV: $vn_mv_kv, LV: $vn_lv_kv"

  rk_T1 = isnothing(trafo.side1.r) ? 0.0 : trafo.side1.r
  xk_T1 = isnothing(trafo.side1.x) ? 0.0 : trafo.side1.x
  bm_T1 = isnothing(trafo.side1.b) ? 0.0 : trafo.side1.b
  gm_T1 = isnothing(trafo.side1.g) ? 0.0 : trafo.side1.g

  rk_T2 = isnothing(trafo.side2.r) ? 0.0 : trafo.side2.r
  xk_T2 = isnothing(trafo.side2.x) ? 0.0 : trafo.side2.x
  bm_T2 = isnothing(trafo.side2.b) ? 0.0 : trafo.side2.b
  gm_T2 = isnothing(trafo.side2.g) ? 0.0 : trafo.side2.g

  rk_T3 = isnothing(trafo.side3.r) ? 0.0 : trafo.side3.r
  xk_T3 = isnothing(trafo.side3.x) ? 0.0 : trafo.side3.x
  bm_T3 = isnothing(trafo.side3.b) ? 0.0 : trafo.side3.b
  gm_T3 = isnothing(trafo.side3.g) ? 0.0 : trafo.side3.g

  r_pu, x_pu, b_pu, g_pu = SparlectraNet.calcTwoPortPU(HV_kv, Sbase_MVA, rk_T1, xk_T1, bm_T1, gm_T1)

  if tapSide != 1
    tapPos = tap_neutral
  else
    tapPos = tap_pos
  end
  T1_ratio = calcRatio(HV_kv, vn_hv_kv, HV_kv, vn_hv_kv, tapPos, tap_neutral, tap_step_percent, 1)
  branch_1 = ResDataTypes.Branch(auxBuxCmp, hv_bus, auxBusIdx, cID, ausBusID, r_pu, x_pu, b_pu, g_pu, T1_ratio, 0.0, 1)

  @debug "branch_1", branch_1

  push!(bVec, branch_1)

  r_pu, x_pu, b_pu, g_pu = SparlectraNet.calcTwoPortPU(MV_kv, Sbase_MVA, rk_T2, xk_T2, bm_T2, gm_T2)
  if tapSide != 2
    tapPos = tap_neutral
  else
    tapPos = tap_pos
  end

  T2_ratio = calcRatio(HV_kv, vn_aux_kv, MV_kv, vn_mv_kv, tapPos, tap_neutral, tap_step_percent, 1)
  branch_2 = ResDataTypes.Branch(auxBuxCmp, auxBusIdx, mv_bus, ausBusID, cID, r_pu, x_pu, b_pu, g_pu, T2_ratio, shift_mv_degree, 1)

  @debug "branch_2", branch_2

  push!(bVec, branch_2)

  r_pu, x_pu, b_pu, g_pu = SparlectraNet.calcTwoPortPU(LV_kv, Sbase_MVA, rk_T3, xk_T3, bm_T3, gm_T3)
  if tapSide != 3
    tapPos = tap_neutral
  else
    tapPos = tap_pos
  end

  T3_ratio = calcRatio(HV_kv, vn_aux_kv, LV_kv, vn_lv_kv, tapPos, tap_neutral, tap_step_percent, 1)
  branch_3 = ResDataTypes.Branch(auxBuxCmp, auxBusIdx, lv_bus, ausBusID, cID, r_pu, x_pu, b_pu, g_pu, T3_ratio, shift_lv_degree, 1)

  @debug "branch_3", branch_3

  push!(bVec, branch_3)
end

#FIXME: PhaseShifter  
function createBranchVectorFromNodeVector!(;nodes::Vector{ResDataTypes.Node}, lines::Vector{ResDataTypes.ACLineSegment}, trafos::Vector{ResDataTypes.PowerTransformer}, Sbase_MVA::Float64, 
                                           shunts::Union{Nothing, Vector{ResDataTypes.Shunt}}=nothing, prosumps::Union{Nothing, Vector{ResDataTypes.ProSumer}}=nothing, log::Bool= false)::Vector{ResDataTypes.Branch}
  branchVector = Vector{ResDataTypes.Branch}()
  LineDict = Dict{String,ResDataTypes.ACLineSegment}()
  auxBusDict = Dict{String,ResDataTypes.Node}()
  TrafoDict = Dict{String,ResDataTypes.PowerTransformer}()
  NodeDict = Dict{String,ResDataTypes.Node}()
  
  for n in nodes
    #@show busIdx = n.busIdx, n.comp.cID
    pgm_comp = ImpPGMComp(n.comp, n.busIdx, n.busIdx)
    n.comp = pgm_comp
    if n.comp.cTyp == ResDataTypes.AuxBus
      auxBusDict[n._auxNodeID] = n
    end    
    NodeDict[n.comp.cID] = n
  end
  
  for line in lines
    LineDict[line.comp.cID] = line    
  end

  for s in shunts
    @show "Shunt:", s
    pgm_comp = ImpPGMComp(s.comp, s.busIdx, s.busIdx)
    s.comp = pgm_comp
  end

  for p in prosumps
    if haskey(NodeDict, p.nodeID)
      n = NodeDict[p.nodeID]
      pgm_comp = ImpPGMComp(p.comp, n.busIdx, n.busIdx)
      p.comp = pgm_comp
      if !isnothing(p.referencePri) && p.referencePri > 0
        if isnothing(n._vm_pu) && !isnothing(p.vm_pu) && p.vm_pu > 0.0
          @show n._vm_pu = p.vm_pu
        end
      end

    else
      @warn "ProSumerDict not found: ", p.nodeID
    end  
  end  
  
  for trafo in trafos
    TrafoDict[trafo.comp.cID] = trafo
  end
  
  eTypesArr::Array{ResDataTypes.ComponentTyp} = [ResDataTypes.LineC, ResDataTypes.Trafo]

  for eType in eTypesArr # automatic sorted      
    for n in nodes
      for t in n.terminals
        if t.comp.cTyp == eType
          seite = t.seite
          fNodeID = n.comp.cID
          Ubase = t.comp.cVN
          ZBase = Ubase^2 / Sbase_MVA
          fromBus = n.busIdx
          comp = t.comp

          if seite == ResDataTypes.Seite1
            id = t.comp.cID
            node2, terminal2 = searchNodeByTerminalId(id, nodes, ResDataTypes.Seite2)
            toBus = node2.busIdx
            tNodeID = node2.comp.cID
            tUbase = terminal2.comp.cVN
            if t.comp.cTyp == ResDataTypes.LineC
              @assert(Ubase == tUbase, "Ubase $(Ubase)terminal1  not equal to terminal2 $(tUbase)!")
            end

            if t.comp.cTyp == ResDataTypes.LineC # handle lines to create branches
              if haskey(LineDict, id)
                line = LineDict[id]                
                newID = "#branch_"*comp.cID
                pgm_comp = ImpPGMComp(comp, fromBus, toBus, newID)
                b = getBranchForLine(pgm_comp, fromBus, toBus, fNodeID, tNodeID, line, ZBase)
                push!(branchVector, b)                                
                line.comp = ImpPGMComp(comp, fromBus, toBus)                
              else
                @warn "could not handle Line $id"
              end
            elseif t.comp.cTyp == ResDataTypes.Trafo # handle trafos to create branches
              if haskey(TrafoDict, id)
                trafo = TrafoDict[id]
                if trafo.isBiWinder # 2WT
                  newID = "#branch_"*comp.cID
                  pgm_comp = ImpPGMComp(comp, fromBus, toBus, newID)                  
                  b = getBranchFor2WT(pgm_comp, fromBus, toBus, fNodeID, tNodeID, trafo, Sbase_MVA)                  
                  trafo.comp = ImpPGMComp(comp, fromBus, toBus)
                  #@show "2WT", fromBus, toBus
                  push!(branchVector, b)
                else # 3WT                      
                  if haskey(auxBusDict, trafo.comp.cID)
                    auxBus = auxBusDict[trafo.comp.cID]
                    #@info "3WT - AuxBus", auxBus
                    node3, terminal3 = searchNodeByTerminalId(id, nodes, ResDataTypes.Seite3)
                    setBranchFor3WT!(fromBus, toBus, node3.busIdx, trafo, auxBus, Sbase_MVA, branchVector, true)
                  else
                    @warn "could not find auxBus for 3WT $(id)"
                  end
                end
              else
                @warn "could not handle trafo $id"
              end
            end # if Line & Trafo
          end # if Seite
        end # if Line or Trafo
      end # for Terminal
    end # for Node
  end # eType  

  setParallelBranches!(branchVector)
  return branchVector
end # function

function nodeComparison(node1::ResDataTypes.Node, node2::ResDataTypes.Node)
  node1.comp.cVN >= node2.comp.cVN
end

function createNetFromTripleStore(endpoint::String, sbase_mva::Union{Nothing,Float64}, netName::String, slackBusName::String, log::Bool)::ResDataTypes.Net
  # Create a Sparlectra network from a triple store
  nodes = Sparlectra.SparqlQueryCGMES.getNodes(endpoint, false)
  lines = Sparlectra.SparqlQueryCGMES.getLines(endpoint, false)
  prosumers = Sparlectra.SparqlQueryCGMES.getProSumption(endpoint, false)
  trafos = Sparlectra.SparqlQueryCGMES.getTrafos(endpoint, false)

  debug = false
  if debug
    for n in nodes
      println(n)
    end
    for l in lines
      println(l)
    end
    for p in prosumers
      println(p)
    end
    for t in trafos
      println(t)
    end
  end

  shunts = Vector{ResDataTypes.Shunt}()
  puNodesIDDict = Dict{String,String}()
  nodesDict = Dict{String,ResDataTypes.Node}()

  sort!(nodes, lt = nodeComparison)

  posibleSlackBusID = ""
  for s in prosumers
    if s.referencePri == 1
      posibleSlackBusID = s.nodeID

      #if s.proSumptionType == ResDataTypes.ProSumptionType.PQ
      #  shunt = ResDataTypes.Shunt(s.nodeID, s.pVal, s.qVal)
      #  push!(shunts, shunt)
      #end
    elseif s.isAPUNode
      if log
        @info "createNetFromTripleStore: PU Node - $(s.nodeID)"
      end
      puNodesIDDict[s.nodeID] = s.nodeID
    end
  end

  idx = 0
  slackBusIdx = 0
  slackBusID = ""
  slackName = ""
  slackNameFound = false
  # set bus numbers
  for n in nodes
    idx += 1
    n.busIdx = idx
    n._kidx = idx
    if !slackNameFound && occursin(slackBusName, n.comp.cName)
      slackNameFound = true
      slackBusIdx = idx
      @info "createNetFromTripleStore: slackBusName - $(slackBusName), slackBusIdx: $(slackBusIdx)"
      slackName = n.comp.cName
      slackBusID = n.comp.cID
    end
    nodesDict[n.comp.cID] = n
  end
  auxBusIdx = idx

  if isnothing(sbase_mva)
    sbase_mva = 0.0
    for s in prosumers
      if isnothing(sbase_mva)
        if !isnothing(s.ratedS)
          if sbase_mva < s.ratedS
            sbase_mva = s.ratedS
          end
        elseif !isnothing(s.pVal) && isnothing(s.qVal)
          s_mva = sqrt(s.pVal^2 + s.qVal^2)
          if sbase_mva < s_mva
            sbase_mva = s_mva
          end
        elseif !isnothing(s.maxP) && !isnothing(s.maxQ)
          s_mva = sqrt(s.maxP^2 + s.maxQ^2)
          if sbase_mva < s_mva
            sbase_mva = s_mva
          end
        end
      end
    end
  end
  sbase_mva = roundUpToNearest100(sbase_mva) # round up to nearest 100's
  sign = -1.0
  for s in prosumers
    e = s.comp
    if isMotor(e) || isExternalNetworkInjection(e) || isLoad(e)
      if haskey(nodesDict, s.nodeID)
        n = nodesDict[s.nodeID]
        addAktivePower!(n, sign * s.pVal)
        addReaktivePower!(n, sign * s.qVal)
      else
        @warn "could not find node for motor $(s.nodeID), $(s.comp.cName)"
      end
    end

    if isGenerator(e)
      if haskey(nodesDict, s.nodeID)
        n = nodesDict[s.nodeID]
        addGenAktivePower!(n, sign * s.pVal)
        addGenReaktivePower!(n, sign * s.qVal)
      else
        @warn "could not find node for generator $(s.nodeID), $(s.comp.cName)"
      end
    end

    if isShunt(e)
      setShuntPower!(n, s.pVal, s.qVal)
    end
  end

  # set up aux-busses for 3WT
  for t in trafos
    if !t.isBiWinder
      auxBusIdx += 1
      cName = "aux_#" * string(auxBusIdx)
      nodeID = string(UUIDs.uuid4())
      vn_aux_kv = t.side1.Vn
      ratedS = t.side1.ratedS
      tnodeId = t.comp.cID

      tVec = Vector{Terminal}()
      auxBuxCmp = Component(nodeID, cName, "AUXBUS", vn_aux_kv)
      push!(tVec, Terminal(auxBuxCmp, ResDataTypes.Seite1))
      push!(tVec, Terminal(auxBuxCmp, ResDataTypes.Seite2))
      push!(tVec, Terminal(auxBuxCmp, ResDataTypes.Seite3))

      auxNode = Node(auxBuxCmp, tVec, auxBusIdx, ResDataTypes.PQ, tnodeId, ratedS, 1, 1, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
      if log
        @info "auxNode: $(auxNode)"
      end
      push!(nodes, auxNode)
    end
  end

  for node in nodes
    if !isnothing(node._ratedS)
      if node._ratedS > sbase_mva
        sbase_mva = node._ratedS
      end
    end

    if node.busIdx == slackBusIdx
      node._nodeType = ResDataTypes.Slack
      #@info "found slackBusIdx: $(slackBusIdx), Name: $(slackName), origin name: $(slackBusID)"
    else
      if haskey(puNodesIDDict, node.comp.cID)
        node._nodeType = ResDataTypes.PV
      else
        node._nodeType = ResDataTypes.PQ
      end
    end
  end

  if log
    @info "SBase = $(sbase_mva) MVA"
    @info "choosen SlackBus-Number: $(slackBusIdx), Name: $(slackName), origin name: $(slackBusID)"
    if posibleSlackBusID != ""
      p = nodesDict[posibleSlackBusID]
      busIdx = p.busIdx
      name = p.comp.cName

      @info "posible SlackBus-Number: $(busIdx), Name: $(name), SlackBusID: $(posibleSlackBusID)"
    end
  end

  fixSequenceNumberInNodeVec!(nodes, trafos)

  branchVec = createBranchVectorFromNodeVector!(nodes=nodes, lines=lines, trafos=trafos, Sbase_MVA=sbase_mva, shunts=shunts, prosumps=prosumers)

  renumberNodeIdx!(branchVec, nodes, log)

  net = Net(netName, sbase_mva, slackBusIdx, nodes, lines, trafos, branchVec, prosumers, shunts)

  return net
end

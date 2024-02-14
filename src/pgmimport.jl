# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 08.02.2024
# include-file pgmimport.jl
# power grid model json importer

# base_MVA = 0.0 for default value in case file, otherwise set to desired value
function createNetFromPGM(filename, base_MVA::Float64 = 0.0, log = false, check = false)::ResDataTypes.Net
  @info "create network from PGM-File: $(filename)"

  nodes, lines, wt2, wt3, sym_gens, sym_loads, shunts, source = SparlectraImport.pgmparser(filename)

  @assert wt3 === nothing "3WT transformers not supported yet"

  base_name = basename(filename)
  netName, ext = splitext(base_name)
  @info "Netname: $netName"
  umrech_MVA = 1e-6

  if base_MVA == 0.0
    # if baseMVA is not given, so we search for the maximum sn of the transformers
    sn_max = 0.0
    for trafo in wt2
      sn = float(trafo["sn"])
      if sn > sn_max
        sn_max = sn
      end
    end
    @assert sn_max > 0.0 "no transformer found!"
    baseMVA = sn_max * umrech_MVA
  else
    baseMVA = base_MVA
  end
  @info "baseMVA: $baseMVA MVA"

  VoltageDict = Dict{Integer,Float64}()
  NodeIDDict = Dict{Integer,String}()

  ACLines = Vector{ResDataTypes.ACLineSegment}()
  trafos = Vector{ResDataTypes.PowerTransformer}()
  prosum = Vector{ResDataTypes.ProSumer}()
  shuntVec = Vector{ResDataTypes.Shunt}()
  branchVec = Vector{ResDataTypes.Branch}()
  nodeVec = Vector{ResDataTypes.Node}()

  NodeTerminalsDict = Dict{Integer,Vector{ResDataTypes.Terminal}}()
  NodeParametersDict = Dict{Integer,ResDataTypes.NodeParameters}()
  busVec = Vector{Bus}()
  # search type of the busses PV, PQ, Slack
  bus_types = Dict{Int,Int}()
  for sym_gen in sym_gens
    bus_id = sym_gen["node"]
    bus_type = sym_gen["type"]
    bus_type += 1
    bus_types[bus_id] = bus_type
  end

  for sym_load in sym_loads
    bus_id = sym_load["node"]
    bus_type = sym_load["type"]
    bus_type += 1
    bus_types[bus_id] = bus_type
  end

  # slack bus  
  vm_pu_slack = 1.0
  slackIdx = 0 # Counter for slack index  
  nodeIdSlack = ""
  slackNodeID_num= 0
  for s in source
    slackNodeID_num = Int64(s["id"])
    bus_id = Int64(s["node"])
    vm_pu_slack = float(s["u_ref"])
    bus_type = 3
    bus_types[bus_id] = bus_type
    slackIdx = bus_id
    break # only one slack bus
  end
  @assert slackIdx != 0 "no slack bus found!"

  for bus in nodes
    busIdx = Int64(bus["id"])
    vn_kv = float(bus["u_rated"]) / 1000.0
    name = "Bus_" * string(busIdx)
    nodeID = "#ID_" * name
    
    vm_pu = nothing
    va_deg = nothing
    btype = 1
    if haskey(bus_types, busIdx)
      btype = bus_types[busIdx]
    end   

    if btype == 3
      nodeIdSlack = nodeID
      vm_pu = vm_pu_slack
      va_deg = 0.0        
    end
    
    a_bus = Bus(busIdx, name, nodeID, "", vn_kv, btype, vm_pu, va_deg)
    push!(busVec, a_bus)

    VoltageDict[busIdx] = vn_kv

    # Set up node IDs via bus index
    NodeIDDict[busIdx] = nodeID

    # initialize NodeTerminalsDict
    NodeTerminalsDict[busIdx] = Vector{ResDataTypes.Terminal}()
    # initialize NodeParametersDict, set bus index
    NodeParameters = ResDataTypes.NodeParameters(busIdx)
    NodeParametersDict[busIdx] = NodeParameters
  end
  
  #Check if the bus numbers are consecutive and unique          
  for (i, bus) in enumerate(busVec)
    @assert bus.busIdx == i "bus numbers are not consecutive and unique" 
  end

  @info "$(length(busVec)) busses created..."
  # Set up auxillary buses

  b = nothing
  g = nothing

  for line in lines
    from = Int64(line["from_node"])
    to = Int64(line["to_node"])

    num = Int64(line["id"])
    cName = "Line_" * string(from) * "_" * string(to)
    cID = "#ID_Line" * string(num)
    r = float(line["r1"])
    x = float(line["x1"])
    c_f = float(line["c1"])
    c_nf = c_f * 1e9
    tan_delta = float(line["tan1"])
    if tan_delta != 0.0
      @warn "tan_delta not zero: $tan_delta"
    end
    from_status = line["from_status"]
    to_status = line["to_status"]
    if from_status == 1 && to_status == 1
      inService = 1
    else
      inService = 0
    end
    Vn1 = VoltageDict[from]
    Vn2 = VoltageDict[to]
    Vn = Vn1

    @assert Vn1 == Vn2 "Voltage levels are different: $Vn1 != $Vn2"
    #length not given, so set to 1.0km
    cPGM = ImpPGMComp(cID, cName, toComponentTyp("ACLINESEGMENT"), Vn, from, to)
    asec = ACLineSegment(cPGM, 1.0, r, x, 0.0, 0.0, c_nf, 0.0)
    push!(ACLines, asec)

    t1 = ResDataTypes.Terminal(cPGM, ResDataTypes.Seite1)
    t2 = ResDataTypes.Terminal(cPGM, ResDataTypes.Seite2)

    checkBusNumber(from, busVec)
    checkBusNumber(to, busVec)

    t2Terminal = NodeTerminalsDict[from]
    t1Terminal = NodeTerminalsDict[to]
    push!(t1Terminal, t1)
    push!(t2Terminal, t2)

    bch = c_f * 2.0 * 50.0 * pi
    g_us = 0.0 # tan delta = 0 

    r_pu, x_pu, b_pu, g_pu = calcTwoPortPU(Vn, baseMVA, r, x, bch, g_us)
    fromNodeID = NodeIDDict[from]
    toNodeID = NodeIDDict[to]
    from_org = from
    to_org = to
    branch = Branch(cPGM, from, to, from_org, to_org, fromNodeID, toNodeID, r_pu, x_pu, b_pu, g_pu, 0.0, 0.0, inService)
    push!(branchVec, branch)
  end
  @info "$(length(ACLines)) aclines created..."

  for t in wt2
    sn = float(t["sn"])
    u1 = float(t["u1"])
    u2 = float(t["u2"])

    from_kV = u1 / 1000.0 # rated voltage at the from side
    to_kV = u2 / 1000.0 # rated voltage at the to side
    if from_kV > to_kV
      vn_hv = from_kV
      vn_lv = to_kV
    else
      vn_hv = to_kV
      vn_lv = from_kV
    end
    oID = Int64(t["id"])
    uk = float(t["uk"])
    i0 = float(t["i0"])
    p0_W = float(t["p0"])
    pk_W = float(t["pk"])
    from_node = Int64(t["from_node"])
    to_node = Int64(t["to_node"])

    from_status = t["from_status"]
    to_status = t["to_status"]
    if from_status == 1 && to_status == 1
      inService = 1
    else
      inService = 0
    end

    cName = "Trafo_" * string(from_node) * "_" * string(to_node)
    cID = "#ID_Trafo" * string(oID)

    if check
      checkBusNumber(from_node, busVec)
      checkBusNumber(to_node, busVec)
    end

    tap_side = Int64(t["tap_side"])
    tap_pos = Int64(t["tap_pos"])
    tap_min = Int64(t["tap_min"])
    tap_max = Int64(t["tap_max"])
    tap_nom = Int64(t["tap_nom"])
    tap_size_V = float(t["tap_size"])
    #side_1 = 0, side_2 = 1, side_3 = 2
    tapSeite = (tap_side == 0) ? ((from_kV > to_kV) ? 1 : 2) : (tap_side == 1) ? ((from_kV > to_kV) ? 2 : 1) : @error "tap_side: $tap_side, from_kV: $from_kV, to_kV: $to_kV"

    tap_neutral = div(tap_max + tap_min, 2)
    HV = VoltageDict[from_node]
    LV = VoltageDict[to_node]

    if tap_side == 0
      tapStepPercent = calcTapStepPercent(tap_size_V, u1)
    else
      tapStepPercent = calcTapStepPercent(tap_size_V, u2)
    end
    shift_degree = 0.0
    sn_MVA = sn / 1e6
    ratio = calcRatio(HV, vn_hv, LV, vn_lv, tap_pos, tap_neutral, tapStepPercent, tapSeite)
    r_pu = nothing
    x_pu = nothing
    b_pu = nothing
    g_pu = nothing
    if tapSeite == 1
      Vtab = calcNeutralU(1.0, vn_hv, tap_min, tap_max, tapStepPercent)
      tap = PowerTransformerTaps(tap_pos, tap_min, tap_max, tap_neutral, tapStepPercent, Vtab)
      r, x, b, g = calcTrafoParamsSI(sn_max, u2, uk, sn, pk_W, i0, p0_W)
      r_pu, x_pu, b_pu, g_pu = calcTwoPortPU(vn_hv, baseMVA, r, x, b, g)

      s1 = PowerTransformerWinding(vn_hv, r, x, b, g, shift_degree, from_kV, sn_MVA, tap)
      s2 = PowerTransformerWinding(vn_lv, 0.0, 0.0)
    else
      Vtab = calcNeutralU(neutralU_ratio, vn_lv, tap_min, tap_max, tapStepPercent)
      tap = PowerTransformerTaps(tap_pos, tap_min, tap_max, tap_neutral, tapStepPercent, Vtab)
      r, x, b, g = calcTrafoParamsSI(sn_max, u1, uk, sn, pk_W, i0, p0_W)
      r_pu, x_pu, b_pu, g_pu = calcTwoPortPU(vn_lv, baseMVA, r, x, b, g)

      s1 = PowerTransformerWinding(vn_hv, 0.0, 0.0)
      s2 = PowerTransformerWinding(vn_lv, r, x, b, g, shift_degree, to_kV, sn_MVA, tap)
    end
    #@debug "Trafo: ", cName, ", r: ", r, ", x: ", x, ", b: ", b, ", g: ", g
    #@debug "Trafo: ", cName, ", r_pu: ", r_pu, ", x_pu: ", x_pu, ", b_pu: ", b_pu, ", g_pu: ", g_pu
    s3 = nothing

    
    addEx = TransformesAdditionalParameters(sn, uk, pk_W, i0, p0_W)
    cImpPGMComp = ImpPGMComp(cID, cName, toComponentTyp("POWERTRANSFORMER"), vn_hv, from_node, to_node)
    trafo = PowerTransformer(cImpPGMComp, true, s1, s2, s3, ResDataTypes.Ratio, addEx)
    push!(trafos, trafo)
    #@show trafo
    
    cmp1 = ResDataTypes.Component(cID, cName, "POWERTRANSFORMER", vn_hv)
    t1 = ResDataTypes.Terminal(cmp1, ResDataTypes.Seite1)
    cmp2 = ResDataTypes.Component(cID, cName, "POWERTRANSFORMER", vn_lv)
    t2 = ResDataTypes.Terminal(cmp2, ResDataTypes.Seite2)

    t1Terminal = NodeTerminalsDict[from_node]
    t2Terminal = NodeTerminalsDict[to_node]

    fromNodeID = NodeIDDict[from_node]
    toNodeID = NodeIDDict[to_node]

    push!(t1Terminal, t1)
    push!(t2Terminal, t2)

    branch = Branch(cImpPGMComp, from_node, to_node, fromNodeID, toNodeID, r_pu, x_pu, b_pu, g_pu, ratio, shift_degree, inService)

    push!(branchVec, branch)
  end
  @info "$(length(trafos)) power transformers created..."

  ratedS = nothing
  ratedU = nothing
  qPercent = nothing
  maxP = nothing
  minP = nothing
  maxQ = nothing
  minQ = nothing
  ratedPowerFactor = nothing
  referencePri = nothing
  lanz = 0

  for sym_load in sym_loads
    status = sym_load["status"]
    if status == 0
      continue
    end
    type = sym_load["type"]
    @assert type == 0 "only constant power loads are supported"
    lanz += 1
    id = sym_load["id"]
    bus = sym_load["node"]
    vn = VoltageDict[bus]
    cName = "Load_" * string(id)
    cID = "#ID_Load_" * string(id)
    nID = NodeIDDict[bus] # String
    p = sym_load["p_specified"] * umrech_MVA
    q = sym_load["q_specified"] * umrech_MVA

    comp = ImpPGMComp(cID, cName, toComponentTyp("LOAD"), vn, bus, bus)
    pRS = ProSumer(comp, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri, nothing, nothing)
    push!(prosum, pRS)
    t1 = Terminal(comp, ResDataTypes.Seite2)
    t1Terminal = NodeTerminalsDict[bus]
    push!(t1Terminal, t1)

    nParms = NodeParametersDict[bus]
    nParms.pƩLoad = isnothing(nParms.pƩLoad) ? p : nParms.pƩLoad + p
    nParms.qƩLoad = isnothing(nParms.qƩLoad) ? q : nParms.qƩLoad + q
    NodeParametersDict[bus] = nParms
  end
  @info "$(lanz) loads created..."

  lanz = 0
  for sym_gen in sym_gens
    status = sym_gen["status"]
    if status == 0
      continue
    end
    type = sym_gen["type"]
    @assert type == 0 "only constant power generators are supported"
    lanz += 1
    id = sym_gen["id"]
    bus = sym_gen["node"]
    vn = VoltageDict[bus]
    cName = "Gen_" * string(id)
    cID = "#ID_Gen_" * string(id)
    nID = NodeIDDict[bus] # String
    p = sym_gen["p_specified"] * umrech_MVA
    q = sym_gen["q_specified"] * umrech_MVA
    
    comp = ImpPGMComp(cID, cName, toComponentTyp("GENERATOR"), vn, bus, bus)
    pRS = ProSumer(comp, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri, nothing, nothing)
    push!(prosum, pRS)
    t2 = Terminal(comp, ResDataTypes.Seite1)
    t2Terminal = NodeTerminalsDict[bus]
    push!(t2Terminal, t2)

    nParms = NodeParametersDict[bus]
    nParms.pƩGen = isnothing(nParms.pƩGen) ? p : nParms.pƩGen + p
    nParms.qƩGen = isnothing(nParms.qƩGen) ? q : nParms.qƩGen + q
    NodeParametersDict[bus] = nParms
  end
  @info "$(lanz) generators created..."

  lanz = 0
  for shunt in shunts
    status = shunt["status"]
    if status == 0
      continue
    end
    lanz += 1
    bus = shunt["node"]
    id = shunt["id"]
    vn_kv = VoltageDict[bus]
    cName = "Shunt_" * string(id)
    cID = "#ID_Shunt_" * string(id)
    
    p_shunt, q_shunt = calcPQ_Shunt(Float64(shunt["g1"]), Float64(shunt["b1"] ), vn_kv)
    Y = Complex(p_shunt, q_shunt)
    y_pu = calc_y_pu(Y, baseMVA, vn_kv)
    comp = ImpPGMComp(cID, cName, toComponentTyp("LINEARSHUNTCOMPENSATOR"), vn_kv, bus, bus)
    sh = Shunt(comp, NodeIDDict[bus], bus, p_shunt, q_shunt, y_pu, status)
    push!(shuntVec, sh)
    t1 = Terminal(comp, ResDataTypes.Seite1)
    t1Terminal = NodeTerminalsDict[bus]
    push!(t1Terminal, t1)

    nParms = NodeParametersDict[bus]
    nParms.pShunt = isnothing(nParms.pShunt) ? p_shunt : nParms.pShunt + p_shunt
    nParms.qShunt = isnothing(nParms.qShunt) ? q_shunt : nParms.qShunt + q_shunt
    NodeParametersDict[bus] = nParms
  end
  @info "$(lanz) shunts created..."

  for b in busVec
    #@show b
    busIdx = b.busIdx
    busName = b.name
    nodeType = b.type
    cID = b.id
    Vn = b.vn_kv
    terminals = NodeTerminalsDict[busIdx]

    if log
      println("Bus: ", busName, " busIdx: ", busIdx, " ID: ", cID, " Voltage: ", Vn, "\nTerminals: ", terminals)
    end

    c = ImpPGMComp(cID, busName, ResDataTypes.Busbarsection, Vn, busIdx, busIdx)

    node = Node(c, terminals, busIdx, busIdx, toNodeType(nodeType))

    nParms = NodeParametersDict[busIdx]
    
    if nodeType == 3      
      nParms.vm_pu = b.vm_pu
      nParms.va_deg = b.va_deg
    end
    
    setNodeParameters!(node, nParms)
    if log
      @show node
    end
    push!(nodeVec, node)
  end

  @info "3WT, not implemented yet"
  if check
    checkNodeConnections(nodeVec)
  end

  setParallelBranches!(branchVec)

  net = ResDataTypes.Net(netName, baseMVA, slackIdx, nodeVec, ACLines, trafos, branchVec, prosum, shuntVec)

  return net
end

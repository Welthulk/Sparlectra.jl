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
    c_nf = c_f  * 1e9
    tan_delta = float(line["tan1"])
 
    bch = c_f * 2.0 * pi * 50.0 
    g_us = 0.0 # tan delta = 0 

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
    if inService == 0
      continue
    end

    Vn1 = VoltageDict[from]
    Vn2 = VoltageDict[to]
    Vn = Vn1

    @assert Vn1 == Vn2 "Voltage levels are different: $Vn1 != $Vn2"
    #length not given, so set to 1.0km
    cPGM = ImpPGMComp(cID, cName, toComponentTyp("ACLINESEGMENT"), Vn, from, to)
    asec = ACLineSegment(cPGM, 1.0, r, x, bch, g_us, c_nf, 0.0)
    push!(ACLines, asec)
    #@show p=get_line_parameters(asec)

    t1 = ResDataTypes.Terminal(cPGM, ResDataTypes.Seite1)
    t2 = ResDataTypes.Terminal(cPGM, ResDataTypes.Seite2)

    checkBusNumber(from, busVec)
    checkBusNumber(to, busVec)

    t2Terminal = NodeTerminalsDict[from]
    t1Terminal = NodeTerminalsDict[to]
    push!(t1Terminal, t1)
    push!(t2Terminal, t2)

    r,x,b,g = getRXBG(asec)    
    _g = isnothing(g) ? 0.0 : g
    _b = isnothing(b) ? 0.0 : b
    r_pu, x_pu, b_pu, g_pu = calcTwoPortPU(Vn, baseMVA, r, x, _b, _g)
    
    fromNodeID = NodeIDDict[from]
    toNodeID = NodeIDDict[to]
    from_org = from
    to_org = to
    branch = Branch(cPGM, from, to, from_org, to_org, fromNodeID, toNodeID, r_pu, x_pu, b_pu, g_pu, 0.0, 0.0, inService)
    push!(branchVec, branch)
  end
  @info "$(length(ACLines)) aclines created..."

  for t in wt2
    sn_MVA = float(t["sn"])*umrech_MVA
    vn_hv_kV = float(t["u1"])*1e-3 # rated voltage at the from side
    vn_lv_kV = float(t["u2"])*1e-3 # rated voltage at the to side

    oID = Int64(t["id"])
    vk_percent = float(t["uk"])*100
    i0_percent = float(t["i0"])*100
    p0_kW = float(t["p0"])*1e-3
    pk_kW = float(t["pk"])*1e-3
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
    tap_size_kV = float(t["tap_size"])*1e-3
    tap_neutral = tap_nom
    HV = VoltageDict[from_node]
    LV = VoltageDict[to_node]
    # 0 = from_side, 1 = to_side
    shift_degree = 0.0
    
    addEx = TransformerModelParameters(sn_MVA=sn_MVA, vk_percent=vk_percent, vkr_percent=nothing, pk_kW=pk_kW, i0_percent=i0_percent, p0_kW=p0_kW)

    r_pu = nothing
    x_pu = nothing
    b_pu = nothing
    g_pu = nothing

    s1 = nothing
    s2 = nothing
    s3 = nothing
    if tap_side == 0 
      
      s2 = PowerTransformerWinding(Vn_kV=vn_lv_kV)
      tap = PowerTransformerTaps(Vn_kV=vn_hv_kV, step=tap_pos, lowStep=tap_min, highStep=tap_max, neutralStep=tap_neutral, voltageIncrement_kV=tap_size_kV)
      s1 = PowerTransformerWinding(Vn_kV=vn_hv_kV, modelData=addEx, shift_degree = shift_degree, ratedU = vn_hv_kV, ratedS=sn_MVA, taps=tap) 
      
      ratio = calcRatio(HV, vn_hv_kV, LV, vn_lv_kV, tap.step, tap.neutralStep, tap.tapStepPercent)
    else

      s1 = PowerTransformerWinding(Vn_kV=vn_hv_kV)
      
      tap = PowerTransformerTaps(Vn_kV=vn_lv_kV, step=tap_pos, lowStep=tap_min, highStep=tap_max, neutralStep=tap_neutral, voltageIncrement_kV=tap_size_kV)      
      s2 = PowerTransformerWinding(Vn_kV=vn_lv_kV, modelData=addEx, shift_degree = shift_degree, ratedU = vn_hv_kV, ratedS=sn_MVA, taps=tap)       
      
      ratio = calcRatio(HV, vn_hv_kV, LV, vn_lv_kV, tap.step, tap.neutralStep, tap.tapStepPercent)
    end

      
    cImpPGMComp = ImpPGMComp(cID, cName, toComponentTyp("POWERTRANSFORMER"), vn_hv_kV, from_node, to_node)
    trafo = PowerTransformer(cImpPGMComp, true, s1, s2, s3)

    push!(trafos, trafo)
        
    cmp1 = ResDataTypes.Component(cID, cName, "POWERTRANSFORMER", vn_hv_kV)
    t1 = ResDataTypes.Terminal(cmp1, ResDataTypes.Seite1)
    cmp2 = ResDataTypes.Component(cID, cName, "POWERTRANSFORMER", vn_lv_kV)
    t2 = ResDataTypes.Terminal(cmp2, ResDataTypes.Seite2)

    t1Terminal = NodeTerminalsDict[from_node]
    t2Terminal = NodeTerminalsDict[to_node]

    fromNodeID = NodeIDDict[from_node]
    toNodeID = NodeIDDict[to_node]

    push!(t1Terminal, t1)
    push!(t2Terminal, t2)
    branch = Branch(branchC=cImpPGMComp, baseMVA=baseMVA, fromNodeID=fromNodeID, toNodeID=toNodeID, trafo=trafo, ratio=ratio, status=inService,isParallel=false) 
        
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
    comp = ImpPGMComp(cID, cName, toComponentTyp("LINEARSHUNTCOMPENSATOR"), vn_kv, bus, bus)
    
    sh = Shunt(comp=comp, nodeID=NodeIDDict[bus],  base_MVA = baseMVA, Vn_kV_shunt = vn_kv, g_shunt = Float64(shunt["g1"]), b_shunt = Float64(shunt["b1"] ) )
    p_shunt = sh.p_shunt
    q_shunt = sh.q_shunt

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

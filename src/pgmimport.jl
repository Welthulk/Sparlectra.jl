# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 08.02.2024
# include-file pgmimport.jl
# power grid model json importer

function identifyIsolatedBuses(nodes, lines, wt2, wt3)
  isolated_buses = Set{Int}()
  connected_buses = Set{Int}()

  for o in lines
    if o["from_status"] == 0 || o["to_status"] == 0
      continue
    end
    push!(connected_buses, o["from_node"])
    push!(connected_buses, o["to_node"])
  end

  for o in wt2
    if o["from_status"] == 0 || o["to_status"] == 0
      continue
    end
    push!(connected_buses, o["from_node"])
    push!(connected_buses, o["to_node"])
  end

  for o in wt3
    if o["status_1"] == 0 || o["status_2"] == 0 || o["status_3"] == 0
      continue
    end
    push!(connected_buses, o["node_1"])
    push!(connected_buses, o["node_2"])
    push!(connected_buses, o["node_3"])
  end

  for node in nodes
    node_id = node["id"]
    if !(node_id in connected_buses)
      push!(isolated_buses, node_id)
    end
  end
  @info "isolated busses: $(isolated_buses)"
  return isolated_buses
end

function reassignBusNumbers!(nodes, lines, wt2, wt3, s_gen, s_load, shunt, source)
  isolated_buses = identifyIsolatedBuses(nodes, lines, wt2, wt3)
  vm_pu = 1.0
  slack_bus = 0
  if length(isolated_buses) > 0
    @info "reassignBusNumbers!: Isolated busses found: $(isolated_buses)"
  end

  bus_mapping = Dict{Int,Int}()
  # type of the busses PV, PQ, Slack
  bus_types = Dict{Int,Int}()

  new_bus_number = 1

  for node in nodes
    old_bus_number = node["id"]

    # Überspringe isolierte Busse
    if old_bus_number in isolated_buses
      node["o_id"] = old_bus_number
      @debug "isolated bus found: $(old_bus_number)"      
      node["type"] = 4
      continue
    end

    bus_mapping[old_bus_number] = new_bus_number
    node["id"] = new_bus_number
    node["o_id"] = old_bus_number
    new_bus_number += 1
  end

  for line in lines
    line["o_from"] = line["from_node"]
    line["o_to"] = line["to_node"]
    line["from_node"] = get(bus_mapping, line["from_node"], line["from_node"])
    line["to_node"] = get(bus_mapping, line["to_node"], line["to_node"])
    if !haskey(line, "length")
      line["length"] = 1.0
    end
  end

  for transformer in wt2
    transformer["o_from"] = transformer["from_node"]
    transformer["o_to"] = transformer["to_node"]

    transformer["from_node"] = get(bus_mapping, transformer["from_node"], transformer["from_node"])
    transformer["to_node"] = get(bus_mapping, transformer["to_node"], transformer["to_node"])
  end

  for transformer in wt3
    for i = 1:3
      node_key = "node_$(i)"
      if haskey(transformer, node_key)
        transformer["o_$(node_key)"] = transformer[node_key]
        transformer[node_key] = get(bus_mapping, transformer[node_key], transformer[node_key])
      end
    end
  end

  for gen in s_gen
    gen["o_node"] = gen["node"]
    gen["node"] = get(bus_mapping, gen["node"], gen["node"])

    bus_id = Int(gen["node"])
    bus_type = gen["type"]
    bus_type += 1
    bus_types[bus_id] = bus_type
  end

  for load in s_load
    load["o_node"] = load["node"]
    load["node"] = get(bus_mapping, load["node"], load["node"])

    bus_id = Int(load["node"])
    bus_type = load["type"]
    bus_type += 1
    bus_types[bus_id] = bus_type
  end

  for s in shunt
    s["o_node"] = s["node"]
    s["node"] = get(bus_mapping, s["node"], s["node"])

    bus_id = Int(s["node"])
    if !haskey(bus_types, bus_id)
      bus_types[bus_id] = 1
    end
  end

  for s in source
    s["o_node"] = s["node"]
    s["node"] = get(bus_mapping, s["node"], s["node"])
    vm_pu = float(s["u_ref"])
    slack_bus = Int(s["node"])
    bus_id = Int(s["node"])
    bus_types[bus_id] = 3
  end

  for node in nodes
    bus_id = node["id"]
    if haskey(node, "type")
      type = node["type"]
      @debug "bus type already set for bus $(bus_id), type: $(type)"
      continue
    end

    if haskey(bus_types, bus_id)
      node["type"] = bus_types[bus_id]
    else
      @debug "bus type not found for bus $(bus_id), set to PQ"
      node["type"] = 1
    end
  end
  return vm_pu, slack_bus
end

function searchMaxS(wt2, wt3, gens)
  sn_max = 0.0
  if !isnothing(wt2)
    for trafo in wt2
      sn = float(trafo["sn"])
      if sn > sn_max
        sn_max = sn
      end
    end
  end
  if !isnothing(wt3)
    for t in wt3
      sn_1 = float(t["sn_1"])
      sn_2 = float(t["sn_2"])
      sn_3 = float(t["sn_3"])
      sn = max(sn_1, sn_2, sn_3)
      if sn > sn_max
        sn_max = sn
      end
    end
  end
  if !isnothing(gens)
    for gen in gens
      p = float(gen["p_specified"])
      q = float(gen["q_specified"])
      sn = sqrt(p^2 + q^2)
      if sn > sn_max
        sn_max = sn
      end
    end
  end
  @assert sn_max != 0.0 "no sn found"
  return sn_max
end

# base_MVA = 0.0 for default value in case file, otherwise set to desired value
function createNetFromPGM(filename, base_MVA::Float64 = 0.0, log = false, check = false)::ResDataTypes.Net
  function printRefVals(bus)
    nParms = NodeParametersDict[bus]

    if !(isnothing(nParms.pƩGen) && isnothing(nParms.qƩGen))
      println("bus:", bus)
      println("GEN p:", nParms.pƩGen)
      println("GEN q:", nParms.qƩGen)
    end

    if !(isnothing(nParms.pƩLoad) && isnothing(nParms.qƩLoad))
      println("bus:", bus)
      println("LOAD p:", nParms.pƩLoad)
      println("LOAD q:", nParms.qƩLoad)
    end
  end

  @info "create network from PGM-File: $(filename)"

  nodes, lines, wt2, wt3, sym_gens, sym_loads, shunts, source = pgmparser(filename)

  @info "search for isolated busses (reassign bus numbers...)"
  vm_pu_slack, slackIdx = reassignBusNumbers!(nodes, lines, wt2, wt3, sym_gens, sym_loads, shunts, source)
  @assert slackIdx != 0 "no slack bus found"
  @info "slack: $slackIdx, vm_pu: $vm_pu_slack"
  if check
    @info "check for isolated busses and slack bus"
    for n in nodes
      @debug n
      id   = Int64(n["id"])
      o_id = Int64(n["o_id"])
      type = Int64(n["type"])
      if id != o_id
        @info "bus number re-assigned: $o_id -> $id"
      end

      if type == 4
        @info "isolated bus found: $id"
      elseif type == 3
        @info "slack bus found: $id"
      end
    end
  end

  base_name = basename(filename)
  netName, ext = splitext(base_name)
  @info "Netname: $netName"
  umrech_MVA = 1e-6

  if base_MVA == 0.0
    # if baseMVA is not given, so we search for the maximum sn of the transformers
    baseMVA = searchMaxS(wt2, wt3, sym_gens) * umrech_MVA
  else
    baseMVA = base_MVA
  end
  @info "baseMVA: $baseMVA MVA"

  VoltageDict = Dict{Integer,Float64}()

  AuxBusDict = Dict{String,Integer}()

  ACLines = Vector{ResDataTypes.ACLineSegment}()
  trafos = Vector{ResDataTypes.PowerTransformer}()
  prosum = Vector{ResDataTypes.ProSumer}()
  shuntVec = Vector{ResDataTypes.Shunt}()
  branchVec = Vector{ResDataTypes.Branch}()
  nodeVec = Vector{ResDataTypes.Node}()

  busIDStringDict = Dict{Int,String}()

  NodeTerminalsDict = Dict{Integer,Vector{ResDataTypes.Terminal}}()
  NodeParametersDict = Dict{Integer,ResDataTypes.NodeParameters}()
  busVec = Vector{Bus}()

  busIdx = 0 # busses are numbered from 1 to n ?
  for bus in nodes
    #@assert Int64(bus["id"]) == busIdx + 1 "bus numbers are not consecutive and unique"
    idx = Int64(bus["id"])
    vn_kv = float(bus["u_rated"]) * 1e-3
    btype = Int64(bus["type"])
    if btype == 4
      @info "isolated bus found, bus number: $idx, vn_kv: $vn_kv, type: $btype, skip..."
      continue
    end
    busIdx = idx
    vm_pu = (btype == 3) ? vm_pu_slack : 1.0
    va_deg = 0.0

    a_bus = Bus(busIdx, vn_kv, btype, vm_pu, va_deg)
    busIDStringDict[busIdx] = a_bus.id

    push!(busVec, a_bus)

    VoltageDict[busIdx] = vn_kv
    # initialize NodeTerminalsDict
    NodeTerminalsDict[busIdx] = Vector{ResDataTypes.Terminal}()
    # initialize NodeParametersDict, set bus index
    NodeParameters = ResDataTypes.NodeParameters(busIdx)
    NodeParametersDict[busIdx] = NodeParameters
  end
  @info "$(length(busVec)) busses created..."

  # Set up auxillary buses
  auxAnz = 0
  if !isnothing(wt3)
    for t in wt3
      @debug "wt3 found..."
      busIdx += 1
      auxAnz += 1
      u1 = float(t["u1"]) * 1e-3
      from_node = Int64(t["node_1"])
      to_node = Int64(t["node_2"])
      to_node_3 = Int64(t["node_3"])

      cName, auxID = getWT3AuxBusID(u1, from_node, to_node, to_node_3)
      bus = Bus(busIdx, cName, auxID, "", u1, 1) # PQ-Bus
      push!(busVec, bus)

      VoltageDict[busIdx] = u1
      AuxBusDict[auxID] = busIdx

      NodeTerminalsDict[busIdx] = Vector{ResDataTypes.Terminal}()

      NodeParameters = ResDataTypes.NodeParameters(busIdx)
      NodeParametersDict[busIdx] = NodeParameters
    end
  end
  auxAnz > 0 ? (@info "$auxAnz auxillary busses created...") : (@info "no aux busses created...")

  b = nothing
  g = nothing
  for line in lines
    from = Int64(line["from_node"])
    to = Int64(line["to_node"])
    oID = Int64(line["id"])

    r = float(line["r1"])
    x = float(line["x1"])
    c_nf = float(line["c1"]) * 1e9
    tan_delta = float(line["tan1"])

    if tan_delta != 0.0
      @warn "tan_delta not zero: $tan_delta, not handelt yet!"
    end

    from_status = line["from_status"]
    to_status = line["to_status"]
    inService = (from_status == 1 && to_status == 1) ? 1 : 0
    if inService == 0
      @info "line $oID not in service"
      continue
    end

    Vn1 = VoltageDict[from]
    Vn2 = VoltageDict[to]
    Vn = Vn1

    @assert Vn1 == Vn2 "Voltage levels are different: $Vn1 != $Vn2"
    cImpPGMComp = getLineImpPGMComp(Vn, from, to)

    #length not given, so set to 1.0km
    asec = ACLineSegment(cImpPGMComp, 1.0, r, x, c_nf, tan_delta)
    push!(ACLines, asec)
    #@show p=get_line_parameters(asec)

    t1 = ResDataTypes.Terminal(cImpPGMComp)
    t2 = ResDataTypes.Terminal(cImpPGMComp)

    checkBusNumber(from, busVec)
    checkBusNumber(to, busVec)

    t2Terminal = NodeTerminalsDict[from]
    t1Terminal = NodeTerminalsDict[to]
    push!(t1Terminal, t1)
    push!(t2Terminal, t2)
    
    branch = Branch(vn_kV= Vn, baseMVA=baseMVA, from=from, to=to, branch=asec, id=Int(oID), status=inService)
    push!(branchVec, branch)
  end
  length(ACLines) > 0 ? (@info "$(length(ACLines)) aclines created...") : (@info "no aclines found...")

  for t in wt2
    sn_MVA = float(t["sn"]) * umrech_MVA
    vn_hv_kV = float(t["u1"]) * 1e-3 # rated voltage at the from side
    vn_lv_kV = float(t["u2"]) * 1e-3 # rated voltage at the to side

    oID = Int64(t["id"])
    vk_percent = float(t["uk"]) * 100
    i0_percent = float(t["i0"]) * 100
    p0_kW = float(t["p0"]) * 1e-3
    pk_kW = float(t["pk"]) * 1e-3
    from_node = Int64(t["from_node"])
    to_node = Int64(t["to_node"])

    from_status = t["from_status"]
    to_status = t["to_status"]
    if from_status == 1 && to_status == 1
      inService = 1
    else
      inService = 0
    end
    if inService == 0
      @info "trafo $oID not in service"
      continue
    end

    if check
      checkBusNumber(from_node, busVec)
      checkBusNumber(to_node, busVec)
    end

    tap_side = Int64(t["tap_side"])
    tap_pos = Int64(t["tap_pos"])
    tap_min = Int64(t["tap_min"])
    tap_max = Int64(t["tap_max"])
    tap_nom = Int64(t["tap_nom"])
    tap_size_kV = float(t["tap_size"]) * 1e-3
    tap_neutral = tap_nom
    HV = VoltageDict[from_node]
    LV = VoltageDict[to_node]
    # 0 = from_side, 1 = to_side
    shift_degree = 0.0

    addEx = TransformerModelParameters(sn_MVA = sn_MVA, vk_percent = vk_percent, vkr_percent = nothing, pk_kW = pk_kW, i0_percent = i0_percent, p0_kW = p0_kW)

    s1 = nothing
    s2 = nothing
    s3 = nothing
    
    if tap_side == 0      
      tap = (tap_min == tap_max) ? nothing : PowerTransformerTaps(Vn_kV = vn_hv_kV, step = tap_pos, lowStep = tap_min, highStep = tap_max, neutralStep = tap_neutral, voltageIncrement_kV = tap_size_kV)
      s1 = PowerTransformerWinding(Vn_kV = HV, ratedU = vn_hv_kV, ratedS = sn_MVA, modelData = addEx, shift_degree = shift_degree,  taps = tap)      
      s2 = PowerTransformerWinding(Vn_kV = LV, ratedU = vn_lv_kV, ratedS = sn_MVA)
    else            
      tap = (tap_min == tap_max) ? nothing : PowerTransformerTaps(Vn_kV = vn_lv_kV, step = tap_pos, lowStep = tap_min, highStep = tap_max, neutralStep = tap_neutral, voltageIncrement_kV = tap_size_kV)
      s1 = PowerTransformerWinding(Vn_kV = HV, ratedU = vn_hv_kV, ratedS = sn_MVA)      
      s2 = PowerTransformerWinding(Vn_kV = LV, ratedU = vn_lv_kV, ratedS = sn_MVA, modelData = addEx, shift_degree = shift_degree,  taps = tap)
    end
    cmp = getTrafoImpPGMComp(false, vn_hv_kV, from_node, to_node)
    trafo = PowerTransformer(cmp, true, s1, s2, s3)
    ratio = calcTransformerRatio(trafo)
    push!(trafos, trafo)

    cmp1 = ResDataTypes.Component(cmp.cID, cmp.cName, "POWERTRANSFORMER", vn_hv_kV)
    t1 = ResDataTypes.Terminal(cmp1)
    cmp2 = ResDataTypes.Component(cmp.cID, cmp.cName, "POWERTRANSFORMER", vn_lv_kV)
    t2 = ResDataTypes.Terminal(cmp2)

    t1Terminal = NodeTerminalsDict[from_node]
    t2Terminal = NodeTerminalsDict[to_node]

    push!(t1Terminal, t1)
    push!(t2Terminal, t2)
    branch = Branch(baseMVA=baseMVA, from=from_node, to=to_node, branch=trafo, id=Int(oID), status=Int(inService))
    push!(branchVec, branch)
  end
  (length(trafos) > 0) ? (@info "$(length(trafos)) 2WTs created......") : (@info "no 2WTs found...")

  anz_3wt = 0
  for t in wt3
    anz_3wt += 1
    oID = Int64(t["id"])
    u1 = float(t["u1"]) * 1e-3 # HV
    u2 = float(t["u2"]) * 1e-3 # MV
    u3 = float(t["u3"]) * 1e-3 # LV
    sn1 = float(t["sn_1"]) * umrech_MVA
    sn2 = float(t["sn_2"]) * umrech_MVA
    sn3 = float(t["sn_3"]) * umrech_MVA
    uk12 = float(t["uk_12"]) * 100
    uk13 = float(t["uk_13"]) * 100
    uk23 = float(t["uk_23"]) * 100
    pk12 = float(t["pk_12"]) * 1e-3
    pk13 = float(t["pk_13"]) * 1e-3
    pk23 = float(t["pk_23"]) * 1e-3
    i0 = float(t["i0"]) * 100
    p0 = float(t["p0"]) * 1e-3
    from_node = Int64(t["node_1"])
    to_node = Int64(t["node_2"])
    to_node_3 = Int64(t["node_3"])

    cName = "3WTTrafo_" * string(oID) * "_" * string(from_node) * "_" * string(to_node) * "_" * string(to_node_3)
    cID = "#ID_3WTTrafo_" * string(oID) * "_" * string(from_node) * "_" * string(to_node) * "_" * string(to_node_3)

    tap_side = Int64(t["tap_side"])
    tap_pos = Int64(t["tap_pos"])
    tap_min = Int64(t["tap_min"])
    tap_max = Int64(t["tap_max"])
    tap_nom = Int64(t["tap_nom"])
    tap_size_kV = float(t["tap_size"]) * 1e-3
    tap_neutral = tap_nom

    HV = VoltageDict[from_node]
    MV = VoltageDict[to_node]
    LV = VoltageDict[to_node_3]


    # create a 3WT - Transformer
    addEx_s1 = TransformerModelParameters(sn_MVA = sn1, vk_percent = uk12, vkr_percent = nothing, pk_kW = pk12, i0_percent = i0, p0_kW = p0)
    addEx_s2 = TransformerModelParameters(sn_MVA = sn2, vk_percent = uk13, vkr_percent = nothing, pk_kW = pk13, i0_percent = i0, p0_kW = p0)
    addEx_s3 = TransformerModelParameters(sn_MVA = sn3, vk_percent = uk23, vkr_percent = nothing, pk_kW = pk23, i0_percent = i0, p0_kW = p0)

    neutral_U = tap_side == 0 ? u1 : tap_side == 1 ? u2 : u3

    tap = (tap_min == tap_max) ? nothing : PowerTransformerTaps(Vn_kV = neutral_U, step = tap_pos, lowStep = tap_min, highStep = tap_max, neutralStep = tap_neutral, voltageIncrement_kV = tap_size_kV)

    u = [u1, u2, u3]
    sn = [sn1, sn2, sn3]

    addEx = [addEx_s1, addEx_s2, addEx_s3]
    sh_deg = [0.0, 0.0, 0.0]
    s1, s2, s3 = create3WTWindings!(u_kV = u, sn_MVA = sn, addEx_Side = addEx, sh_deg = sh_deg, tap_side = (tap_side + 1), tap = tap)

    cImpPGMComp = ImpPGMComp3WT(cID, cName, toComponentTyp("POWERTRANSFORMER"), u[1], from_node, to_node, to_node_3)
    wt3Trafo = PowerTransformer(cImpPGMComp, (!isnothing(tap)), s1, s2, s3)
    push!(trafos, wt3Trafo)

    # aux bus
    auxBusName, auxID = getWT3AuxBusID(u1, from_node, to_node, to_node_3)
    # search aux bus    
    AuxBusIdx = AuxBusDict[auxID]
    vn_aux_kv = VoltageDict[AuxBusIdx]

    auxBuxCmp = ImpPGMComp(auxID, auxBusName, toComponentTyp("POWERTRANSFORMER"), vn_aux_kv, from_node, AuxBusIdx)
    # create terminals for HV_Bus -> Aux_Bus
    t_aux = Terminal(auxBuxCmp)
    t_hv  = Terminal(cImpPGMComp)
    # search terminals 
    t1Terminal = NodeTerminalsDict[from_node]
    auxTerminal = NodeTerminalsDict[AuxBusIdx]
    # save
    push!(t1Terminal, t_hv)
    push!(auxTerminal, t_aux)
    inService = t["status_1"] == 1
    
    side = 1
    ratio = 1.0  
    
    branch = Branch(baseMVA=baseMVA, from=from_node, to=to_node, branch=wt3Trafo, id=Int(oID), status=Int(inService), side=side, ratio=ratio)
    push!(branchVec, branch)

    # create branch for the T2 Aux -> MV (u2)
    t2Terminal = NodeTerminalsDict[to_node]
    t_T2_aux = ResDataTypes.Terminal(auxBuxCmp)
    auxBusName2, auxID2 = getWT3AuxBusID(u2, from_node, to_node, to_node_3)
    auxBuxCmp2 = ImpPGMComp(auxID2, auxBusName2, toComponentTyp("POWERTRANSFORMER"), u2, to_node, AuxBusIdx)
    t_T2_mvbus = ResDataTypes.Terminal(auxBuxCmp2)
    push!(auxTerminal, t_T2_aux)
    push!(t2Terminal, t_T2_mvbus)

    inService = t["status_2"] == 1
    side = 2
    ratio = (vn_aux_kv/MV)*(u2/u1)
    @debug "3WT side2 ratio: ", ratio
    branch = Branch(baseMVA=baseMVA, from=to_node, to=AuxBusIdx, branch=wt3Trafo, id=Int(oID), status=Int(inService), side=side, ratio=ratio)    
    push!(branchVec, branch)

    # create branch for the T3 AUX-> LV (u3)
    t3Terminal = NodeTerminalsDict[to_node_3]
    t_T3_aux = ResDataTypes.Terminal(auxBuxCmp)

    auxBusName3, auxID3 = getWT3AuxBusID(u3, from_node, to_node, to_node_3)
    auxBuxCmp3 = ImpPGMComp(auxID3, auxBusName3, toComponentTyp("POWERTRANSFORMER"), u3, AuxBusIdx, to_node_3)
    t_T3_lvbus = ResDataTypes.Terminal(auxBuxCmp3)
    push!(auxTerminal, t_T3_aux)
    push!(t3Terminal, t_T3_lvbus)

    inService = t["status_3"] == 1
    side = 3
    ratio = (vn_aux_kv/LV)*(u3/u1)
    @debug "3WT side3 ratio: ", ratio    
    branch = Branch(baseMVA=baseMVA, from=AuxBusIdx, to=to_node_3, branch=wt3Trafo, id=Int(oID), status=Int(inService), side=side, ratio=ratio)    
    push!(branchVec, branch)
  end
  (anz_3wt > 0) ? (@info "$(anz_3wt) 3WTs created......") : (@info "no 3WTs found...")

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
    oID = sym_load["id"]
    status = sym_load["status"]
    if status == 0
      @info "load $(oID) not in service"
      continue
    end
    type = sym_load["type"]
    @assert type == 0 "only constant power loads are supported"
    lanz += 1

    bus = sym_load["node"]
    vn = VoltageDict[bus]

    p = sym_load["p_specified"] * umrech_MVA
    q = sym_load["q_specified"] * umrech_MVA
    if abs(p) < 1e-8 && abs(q) < 1e-8
      @info "load $(oID) has no power"
      continue
    end
    comp = getProSumPGMComp(vn, bus, false, Int(oID))
    pRS = ProSumer(comp, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri, nothing, nothing)    
    push!(prosum, pRS)
    t1 = Terminal(comp)
    t1Terminal = NodeTerminalsDict[bus]
    push!(t1Terminal, t1)

    nParms = NodeParametersDict[bus]
    nParms.pƩLoad = isnothing(nParms.pƩLoad) ? p : nParms.pƩLoad + p
    nParms.qƩLoad = isnothing(nParms.qƩLoad) ? q : nParms.qƩLoad + q

    NodeParametersDict[bus] = nParms
  end
  (lanz > 0) ? (@info "$(lanz) loads created...") : (@info "no loads found...")

  lanz = 0
  for sym_gen in sym_gens
    type = sym_gen["type"]
    @assert type == 0 "only constant power generators are supported"

    oID = sym_gen["id"]
    bus = sym_gen["node"]

    vn = VoltageDict[bus]
    status = sym_gen["status"]
    if status == 0
      @info "generator $(oID) not in service"
      continue
    end
    lanz += 1
    p = sym_gen["p_specified"] * umrech_MVA
    q = sym_gen["q_specified"] * umrech_MVA
    if abs(p) < 1e-8 && abs(q) < 1e-8
      @info "generator $(oID) has no power"
      continue
    end
    comp = getProSumPGMComp(vn, bus, true, Int(oID))    
    pRS = ProSumer(comp, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri, nothing, nothing)    
    push!(prosum, pRS)
    t2 = Terminal(comp)
    t2Terminal = NodeTerminalsDict[bus]
    push!(t2Terminal, t2)

    nParms = NodeParametersDict[bus]
    nParms.pƩGen = isnothing(nParms.pƩGen) ? p : nParms.pƩGen + p
    nParms.qƩGen = isnothing(nParms.qƩGen) ? q : nParms.qƩGen + q
    NodeParametersDict[bus] = nParms
  end
  (lanz > 0) ? (@info "$(lanz) generators created...") : (@assert false, "no generators found")

  lanz = 0
  for shunt in shunts
    oID = shunt["id"]
    status = shunt["status"]
    if status == 0
      @info "shunt $(oID) not in service"
      continue
    end
    lanz += 1
    bus = shunt["node"]

    vn_kv = VoltageDict[bus]
    # for PGM the voltage is set to 1.0 kV    
    sh = Shunt(fromBus=bus, id=oID, base_MVA=baseMVA, Vn_kV_shunt=vn_kv, g_shunt= Float64(shunt["g1"]), b_shunt= Float64(shunt["b1"]), ratio=1.0, status = status)
    
    p_shunt, q_shunt = getPQShunt(sh)

    push!(shuntVec, sh)
    t1 = Terminal(sh.comp)
    t1Terminal = NodeTerminalsDict[bus]
    push!(t1Terminal, t1)

    nParms = NodeParametersDict[bus]
    nParms.pShunt = isnothing(nParms.pShunt) ? p_shunt : nParms.pShunt + p_shunt
    nParms.qShunt = isnothing(nParms.qShunt) ? q_shunt : nParms.qShunt + q_shunt
    NodeParametersDict[bus] = nParms
  end
  (lanz > 0) ? (@info "$(lanz) shunts created...") : (@info "no shunts found...")

  for b in busVec
    busIdx = b.busIdx
    #busName = b.name
    nodeType = b.type
    #cID = b.id
    Vn = b.vn_kv
    terminals = NodeTerminalsDict[busIdx]
    
    node = Node(busIdx= busIdx, Vn_kV = Vn, t=terminals, nodeType = toNodeType(nodeType))
    nParms = NodeParametersDict[busIdx]

    if nodeType == 3
      nParms.vm_pu = b.vm_pu
      nParms.va_deg = b.va_deg
    end
    
    if check
      printRefVals(busIdx)
    end
    
    setNodeParameters!(node, nParms)
    push!(nodeVec, node)
  end

  if check
    checkNodeConnections(nodeVec)
  end

  length(branchVec) > 0 ? (@info "$(length(branchVec)) branches created...") : (@assert "no branches found!")
  

  net = ResDataTypes.Net(netName, baseMVA, slackIdx, nodeVec, ACLines, trafos, branchVec, prosum, shuntVec)

  return net
end

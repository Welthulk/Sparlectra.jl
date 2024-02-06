# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 01.8.2023
# include-file readnetfromfile.jl
# Bus
struct Bus
  busIdx::Int64
  name::String
  id::String
  wt3id::String
  vn_kv::Float64
  type::Int64

  function Bus(busIdx::Int64, name::String, id::String, wt3id::String, vn_kv::Float64, type::Int64)
    new(busIdx, name, id, wt3id, vn_kv, type)
  end
  function Base.show(io::IO, bus::Bus)
    print(io, "busIdx:", bus.busIdx, ", bus: ", bus.name, ", id: ", bus.id, ", wt3id: ", bus.wt3id, ", vn: ", bus.vn_kv, ", type: ", bus.type)
  end
end

function readTapMod(tap::Dict{String,Any})
  side = tap["tap_side"]
  tap_min = Int(tap["tap_min"])
  tap_max = Int(tap["tap_max"])
  tap_neutral = Int(tap["tap_neutral"])
  tap_step_percent = float(tap["tap_step_percent"])
  shift_degree = float(tap["shift_degree"])
  neutralU_ratio = float(tap["neutralU_ratio"])

  return side, tap_min, tap_max, tap_neutral, tap_step_percent, shift_degree, neutralU_ratio
end

function checkBusNumber(bus::Int64, busVec::Vector{Bus})::Bool
  for b in busVec
    if b.busIdx == bus
      return true
    end
  end
  @warn "bus $bus not found!"
  return false
end

# base_MVA = 0.0 for default value in case file, otherwise set to desired value
function createNetFromFile(filename, base_MVA::Float64 = 0.0, log::Bool = false)::ResDataTypes.Net
  debug = false
  log_println(message) = log && println(message)
  println("create network from file: $(filename)")
  netName, baseMVA, parsed_data = SparlectraImport.jsonparser(filename, debug)
  if base_MVA > 0.0
    baseMVA = base_MVA
  end
  # parsing the data
  if log
    println("Parsing the data...")
    println("Netname: ", netName)
    println("BaseMVA: ", baseMVA, "\n")
  end

  buses = parsed_data["buses"]
  linemod = parsed_data["linemod"]
  lines = parsed_data["lines"]
  vkDepChr = parsed_data["vkDepChr"]
  tapMod = parsed_data["tapmod"]
  trafo2WT = parsed_data["trafowt2"]
  trafo3WT = parsed_data["trafowt3"]
  shunts = parsed_data["shunts"]
  loads = parsed_data["loads"]
  sgens = parsed_data["sgens"]
  vgens = parsed_data["vgens"]
  grids = parsed_data["ext_grids"]

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

  log_println("Buses:")
  busIdx = 0 # Counter for bus index
  slackIdx = 0 # Counter for slack index
  for bus in buses
    busIdx = Int64(bus["bus"])
    nodeID = string(bus["id"])
    name = string(bus["name"])
    vn_kv = float(bus["vn_kv"])
    btype = Int64(bus["type"])

    if bus["type"] == 3
      slackIdx = busIdx
    end

    bus = Bus(busIdx, name, nodeID, "", vn_kv, btype)
    if log
      @show bus
    end
    push!(busVec, bus)

    # Set up voltage levels via bus index
    Voltage = vn_kv
    VoltageDict[busIdx] = Voltage

    # Set up node IDs via bus index
    NodeIDDict[busIdx] = nodeID

    # initialize NodeTerminalsDict
    NodeTerminalsDict[busIdx] = Vector{ResDataTypes.Terminal}()
    # initialize NodeParametersDict, set bus index
    NodeParameters = ResDataTypes.NodeParameters(busIdx)
    NodeParametersDict[busIdx] = NodeParameters
  end
  @assert slackIdx != 0 "no slack bus found!"

  #Check if the bus numbers are consecutive and unique        
  for (i, bus) in enumerate(busVec)
    if bus.busIdx != i
      @error "bus numbers are not consecutive and unique"
    end
  end

  # Set up auxillary buses
  if !isnothing(trafo3WT)
    log_println("aux buses:")
    for t3WT in trafo3WT
      busIdx += 1
      cName = "aux_#" * string(busIdx)
      nodeID = string(UUIDs.uuid4())
      cID = string(t3WT["id"])

      vn_hv_kv = float(t3WT["vn_hv_kv"])
      type = 1 # PQ
      bus = Bus(busIdx, cName, nodeID, cID, vn_hv_kv, type)
      if log
        @show bus
      end
      push!(busVec, bus)

      VoltageDict[busIdx] = vn_hv_kv
      NodeIDDict[busIdx] = nodeID
      # initialize NodeTerminalsDict
      NodeTerminalsDict[busIdx] = Vector{ResDataTypes.Terminal}()
      # initialize NodeParametersDict, set bus index
      NodeParameters = ResDataTypes.NodeParameters(busIdx)
      NodeParametersDict[busIdx] = NodeParameters
    end
  end
  @assert slackIdx != 0 "no slack bus found!"

  log_println("\nLines:")
  b = nothing
  g = nothing
  angMin = -360.0
  angMax = 360.0
  ratio = 0.0 # ratio = 0.0 for ACLineSegment
  angle = 0.0 # angle = 0.0 for ACLineSegment
  mod = "" # 
  for line in lines
    cName = string(line["name"])
    cID = string(line["id"])
    length = float(line["length_km"])
    r = 0.0
    x = 0.0
    c_nf_per_km = 0.0
    g_us = 0.0
    max_i_ka = 0.0
    try
      mod = line["line_model"]
    catch
      mod = ""
    end
    if mod != "" && !isnothing(mod)
      found = false
      for m in linemod
        if m["model"] == mod
          r = float(m["r_ohm_per_km"])
          x = float(m["x_ohm_per_km"])
          c_nf_per_km = float(m["c_nf_per_km"])
          g_us = float(m["g_us_per_km"])
          max_i_ka = float(m["max_i_ka"])
          found = true
          if log
            println("found model: ", mod)
          end
          break
        end
      end
      @assert found == true "line model $mod not found"
    else
      r = float(line["r_ohm_per_km"])
      x = float(line["x_ohm_per_km"])
      c_nf_per_km = float(line["c_nf_per_km"])
      g_us = float(line["g_us_per_km"])
      max_i_ka = float(line["max_i_ka"])
    end

    Bus = line["from_bus"]

    Vn = VoltageDict[Bus]

    status = line["in_service"]
    asec = ACLineSegment(cID, cName, Vn, length, r, x, b, g, c_nf_per_km)
    if log
      println(asec)
    end
    push!(ACLines, asec)

    cLine = ResDataTypes.toComponentTyp("ACLINESEGMENT")
    c = ResDataTypes.Component(cID, cName, cLine, Vn)

    t1 = ResDataTypes.Terminal(c, ResDataTypes.Seite1)
    t2 = ResDataTypes.Terminal(c, ResDataTypes.Seite2)
    from_bus = line["from_bus"]
    to_bus = line["to_bus"]
    checkBusNumber(from_bus, busVec)
    checkBusNumber(to_bus, busVec)

    t1Terminal = NodeTerminalsDict[from_bus]
    t2Terminal = NodeTerminalsDict[to_bus]
    push!(t1Terminal, t1)
    push!(t2Terminal, t2)
    bch = c_nf_per_km * 1e-9 * length * 2.0 * 50.0 * pi

    r_pu, x_pu, b_pu, g_pu = calcTwoPortPU(Vn, baseMVA, r * length, x * length, bch, g_us * length)

    if log
      println("r_pu: ", r_pu, ", x_pu: ", x_pu, ", b_pu: ", b_pu, ", g_pu: ", g_pu)
    end

    fromNodeID = NodeIDDict[from_bus]
    toNodeID = NodeIDDict[to_bus]
    if log
      println("fromNodeID: ", fromNodeID, ", toNodeID: ", toNodeID, ", from_bus: ", from_bus, ", to_bus: ", to_bus, ", r_pu: ", r_pu, ", x_pu: ", x_pu, ", b_pu: ", b_pu, ", g_pu: ", g_pu, ", ratio: ", ratio, ", angle: ", angle, ", status: ", status, ", angMin: ", angMin, ", angMax: ", angMax)
    end
    branch = ResDataTypes.Branch(c, from_bus, to_bus, from_bus, to_bus, fromNodeID, toNodeID, r_pu, x_pu, b_pu, g_pu, ratio, angle, status)
    if log
      @show branch
    end
    push!(branchVec, branch)
  end

  log_println("\nTransformers (2WT):")
  for t in trafo2WT
    sn = float(t["sn_mva"])
    vn_hv = float(t["vn_hv_kv"])
    vn_lv = float(t["vn_lv_kv"])
    vk_percent = float(t["vk_percent"])
    vkr_percent = float(t["vkr_percent"])
    pfe_kw = float(t["pfe_kw"])
    io_percent = float(t["i0_percent"])
    shift_degree = 0.0

    cName = string(t["name"])
    cID = string(t["id"])
    busIdx = t["hv_bus"]
    USIdx = t["lv_bus"]
    checkBusNumber(busIdx, busVec)
    checkBusNumber(USIdx, busVec)
    inService = t["in_service"]
    nID = NodeIDDict[busIdx]

    tMod = 0
    try
      tMod = t["tap_model"]
    catch
    end

    tap_pos = 0
    try
      tap_pos = t["tap_pos"]
    catch
    end

    vk_dependence = ""
    try
      vk_dependence = t["vk_dependence"]
    catch
    end

    tapSeite = 0
    tap_min = 0
    tap_max = 0
    tap_neutral = 0
    tap_step_percent = 0.0
    neutralU_ratio = 1.0
    regelungEin = false
    if tMod > 0
      regelungEin = true
      for tap in tapMod
        m = tap["model"]
        if m == tMod
          side, tap_min, tap_max, tap_neutral, tap_step_percent, shift_degree, neutralU_ratio = readTapMod(tap)

          if side == "HV" || side == "hv" || side == "H" || side == "h" || side == "1" || side == 1
            tapSeite = 1
          else
            tapSeite = 2
          end

          if log
            println("tap_side: ", side, ", tap_min: ", tap_min, ", tap_max: ", tap_max, ", tap_neutral: ", tap_neutral, ", tap_step_percent: ", tap_step_percent, ", neutralU_ratio:", neutralU_ratio)
          end

          break
        end
      end
    end

    # handle vk dependence
    tapVec = []
    vkVec = []
    xDependence = false
    vkcorr = 0.0
    if vk_dependence != raw""
      for dep in vkDepChr
        if vk_dependence in keys(vkDepChr)
          dep = vkDepChr[vk_dependence]

          tapVec = dep["tap"]
          vkVec = dep["vk"]

          xTaps = [Int64(x) for x in tapVec]
          yVKs = [Float64(y) for y in vkVec]

          vkcorr = calcVKDependence(xTaps, yVKs, Float64(tap_pos))

          xDependence = true
          break
        end
      end
    end

    if xDependence
      vk_percent = vkcorr
    end
    r_pu = 0.0
    x_pu = 0.0
    b_pu = 0.0
    g_pu = 0.0
    HV_Voltage = VoltageDict[busIdx]
    LV_Voltage = VoltageDict[USIdx]
    ratio = calcRatio(HV_Voltage, vn_hv, LV_Voltage, vn_lv, tap_pos, tap_neutral, tap_step_percent, tapSeite)
    if regelungEin && tapSeite == 1
      Vtab = calcNeutralU(neutralU_ratio, vn_hv, tap_min, tap_max, tap_step_percent)
      tap = ResDataTypes.PowerTransformerTaps(tap_pos, tap_min, tap_max, tap_neutral, tap_step_percent, Vtab)
      r_hv, x_hv, b_hv, g_hv, = calcTrafoParams(sn, vn_hv, vk_percent, vkr_percent, pfe_kw, io_percent)
      r_pu, x_pu, b_pu, g_pu = calcTwoPortPU(vn_hv, sn, r_hv, x_hv, b_hv, g_hv)
      s1 = PowerTransformerWinding(vn_hv, r_hv, x_hv, b_hv, g_hv, shift_degree, vn_hv, sn, tap)
      s2 = PowerTransformerWinding(vn_lv, 0.0, 0.0)

    elseif regelungEin && tapSeite == 2
      Vtab = calcNeutralU(neutralU_ratio, vn_lv, tap_min, tap_max, tap_step_percent)
      tap = ResDataTypes.PowerTransformerTaps(tap_pos, tap_min, tap_max, tap_neutral, tap_step_percent, Vtab)
      # maibe an error here: vn_hv?
      r_lv, x_lv, b_lv, g_lv = calcTrafoParams(sn, vn_hv, vk_percent, vkr_percent, pfe_kw, io_percent)
      r_pu, x_pu, b_pu, g_pu = calcTwoPortPU(vn_hv, sn, r_lv, x_lv, b_lv, g_lv)
      s1 = PowerTransformerWinding(vn_hv, 0.0, 0.0)
      s2 = PowerTransformerWinding(vn_lv, r_lv, x_lv, b_lv, g_lv, shift_degree, vn_lv, sn, tap)

    else
      r_hv, x_hv, b_hv, g_hv = calcTrafoParams(sn, vn_hv, vk_percent, vkr_percent, pfe_kw, io_percent)
      r_pu, x_pu, b_pu, g_pu = calcTwoPortPU(vn_hv, sn, r_hv, x_hv, b_hv, g_hv)
      s1 = PowerTransformerWinding(vn_hv, r_hv, x_hv, g_hv, b_hv, vn_hv, sn, nothing)
      s2 = PowerTransformerWinding(vn_lv, 0.0, 0.0)
    end
    s3 = nothing
    cmp = ResDataTypes.Component(cID, cName, "POWERTRANSFORMER", vn_hv)
    trafo = PowerTransformer(cmp, regelungEin, s1, s2, s3)
    if log
      @show trafo
    end
    push!(trafos, trafo)

    t1 = ResDataTypes.Terminal(cmp, ResDataTypes.Seite1)
    cmp2 = ResDataTypes.Component(cID, cName, "POWERTRANSFORMER", vn_lv)
    t2 = ResDataTypes.Terminal(cmp2, ResDataTypes.Seite2)

    t1Terminal = NodeTerminalsDict[busIdx]
    t2Terminal = NodeTerminalsDict[USIdx]

    toNodeID = NodeIDDict[USIdx]

    push!(t1Terminal, t1)
    push!(t2Terminal, t2)

    branch = ResDataTypes.Branch(cmp, busIdx, USIdx, nID, toNodeID, r_pu, x_pu, b_pu, g_pu, ratio, shift_degree, inService)
    if log
      @show branch
    end
    push!(branchVec, branch)
  end

  if !isnothing(trafo3WT)
    log_println("\nTransformers (3WT):")
    for t3WT in trafo3WT
      """
       comment: This code section is programmed "straight forward" to better understand the model. 
       But because of this, the section contains many code duplications.         
      """
      hv_bus = t3WT["hv_bus"]
      checkBusNumber(hv_bus, busVec)

      mv_bus = t3WT["mv_bus"]
      checkBusNumber(mv_bus, busVec)

      lv_bus = t3WT["lv_bus"]
      checkBusNumber(lv_bus, busVec)

      cName = t3WT["name"]
      cID = string(t3WT["id"])

      vn_hv_kv = t3WT["vn_hv_kv"]
      vn_mv_kv = t3WT["vn_mv_kv"]
      vn_lv_kv = t3WT["vn_lv_kv"]

      HV_kv = vn_hv_kv
      MV_kv = vn_mv_kv
      LV_kv = vn_lv_kv

      sn_hv_mva = t3WT["sn_hv_mva"]
      sn_mv_mva = t3WT["sn_mv_mva"]
      sn_lv_mva = t3WT["sn_lv_mva"]

      vk_hv_percent = t3WT["vk_hv_percent"]
      vk_mv_percent = t3WT["vk_mv_percent"]
      vk_lv_percent = t3WT["vk_lv_percent"]

      vkr_hv_percent = t3WT["vkr_hv_percent"]
      vkr_mv_percent = t3WT["vkr_mv_percent"]
      vkr_lv_percent = t3WT["vkr_lv_percent"]

      pfe_kw = t3WT["pfe_kw"]
      i0_percent = t3WT["i0_percent"]

      shift_mv_degree = t3WT["shift_mv_degree"]
      shift_lv_degree = t3WT["shift_lv_degree"]
      tap_model = 0
      try
        tap_model = t3WT["tap_model"]
      catch
      end
      tap_pos = nothing
      try
        tap_pos = t3WT["tap_pos"]
      catch
      end
      vk_dependence = nothing
      try
        vk_dependence = t3WT["vk_dependence"]
      catch
      end
      in_service = 1
      try
        in_service = t3WT["in_service"]
      catch
      end
      if log
        println(
          "hv_bus: ",
          hv_bus,
          ", mv_bus: ",
          mv_bus,
          ", lv_bus: ",
          lv_bus,
          ", name: ",
          cName,
          ", id: ",
          cID,
          ", vn_hv_kv: ",
          vn_hv_kv,
          ", vn_mv_kv: ",
          vn_mv_kv,
          ", vn_lv_kv: ",
          vn_lv_kv,
          ", sn_hv_mva: ",
          sn_hv_mva,
          ", sn_mv_mva: ",
          sn_mv_mva,
          ", sn_lv_mva: ",
          sn_lv_mva,
          ", vk_hv_percent: ",
          vk_hv_percent,
          ", vk_mv_percent: ",
          vk_mv_percent,
          ", vk_lv_percent: ",
          vk_lv_percent,
          ", vkr_hv_percent: ",
          vkr_hv_percent,
          ", vkr_mv_percent: ",
          vkr_mv_percent,
          ", vkr_lv_percent: ",
          vkr_lv_percent,
          ", pfe_kw: ",
          pfe_kw,
          ", i0_percent: ",
          i0_percent,
          ", shift_mv_degree: ",
          shift_mv_degree,
          ", shift_lv_degree: ",
          shift_lv_degree,
          ", tap_model: ",
          tap_model,
          ", tap_pos: ",
          tap_pos,
          ", vk_dependence: ",
          vk_dependence,
          ", in_service: ",
          in_service,
        )
      end
      tapSide = 0
      tap_min = 0
      tap_max = 0
      tap_neutral = 0
      tap_step_percent = 0.0
      shift_degree = 0.0
      neutralU_ratio = 1.0
      if isnothing(tap_pos)
        regelungEin = false
      else
        regelungEin = true
      end
      if tap_model > 0
        for tap in tapMod
          m = tap["model"]
          if m == tap_model
            side, tap_min, tap_max, tap_neutral, tap_step_percent, shift_degree, neutralU_ratio = readTapMod(tap)

            if side == "HV" || side == "hv" || side == "H" || side == "h" || side == "1" || side == 1
              tapSide = 1
            elseif side == "MV" || side == "mv" || side == "M" || side == "m" || side == "2" || side == 2
              tapSide = 2
            elseif side == "LV" || side == "lv" || side == "L" || side == "l" || side == "3" || side == 3
              tapSide = 3
            else
              @error "tap side $side not found"
            end

            if log
              println("tap_side: ", side, ", tap_min: ", tap_min, ", tap_max: ", tap_max, ", tap_neutral: ", tap_neutral, ", tap_step_percent: ", tap_step_percent, ", neutralU_ratio:", neutralU_ratio)
            end

            break
          end
        end
      end
      # handle vk dependence
      tapVec = []
      vkVec = []
      xDependence = false
      vkcorr = 0.0
      if vk_dependence != raw"" && !isnothing(vk_dependence)
        for dep in vkDepChr
          if vk_dependence in keys(vkDepChr)
            dep = vkDepChr[vk_dependence]

            tapVec = dep["tap"]
            vkVec = dep["vk"]

            xTaps = [Int64(x) for x in tapVec]
            yVKs = [Float64(y) for y in vkVec]

            vkcorr = calcVKDependence(xTaps, yVKs, Float64(tap_pos))

            xDependence = true
            break
          end
        end
      end

      tap = nothing
      neutralU = 0.0
      ntap = calcTapCorr(tap_pos, tap_neutral, tap_step_percent)

      if tapSide == 1
        neutralU = vn_hv_kv
        vn_hv_kv = neutralU * ntap
        if xDependence
          vk_hv_percent = vkcorr
        end
      elseif tapSide == 2
        neutralU = vn_mv_kv
        vn_mv_kv = neutralU * ntap
        if xDependence
          vk_mv_percent = vkcorr
        end
      elseif tapSide == 3
        neutralU = vn_lv_kv
        vn_lv_kv = neutralU * ntap
        if xDependence
          vk_lv_percent = vkcorr
        end
      end
      rk_T1, xk_T1, bm_T1, gm_T1 = calc3WTParams(1, sn_hv_mva, sn_mv_mva, sn_lv_mva, HV_kv, MV_kv, LV_kv, vk_hv_percent, vk_mv_percent, vk_lv_percent, vkr_hv_percent, vkr_mv_percent, vkr_lv_percent, pfe_kw, i0_percent)

      rk_T2, xk_T2, bm_T2, gm_T2 = calc3WTParams(2, sn_hv_mva, sn_mv_mva, sn_lv_mva, HV_kv, MV_kv, LV_kv, vk_hv_percent, vk_mv_percent, vk_lv_percent, vkr_hv_percent, vkr_mv_percent, vkr_lv_percent, 0.0, 0.0)

      rk_T3, xk_T3, bm_T3, gm_T3 = calc3WTParams(3, sn_hv_mva, sn_mv_mva, sn_lv_mva, HV_kv, MV_kv, LV_kv, vk_hv_percent, vk_mv_percent, vk_lv_percent, vkr_hv_percent, vkr_mv_percent, vkr_lv_percent, 0.0, 0.0)

      if log
        println("rk_T1: ", rk_T1, ", xk_T1: ", xk_T1, ", bm_T1: ", bm_T1, ", gm_T1: ", gm_T1)
        println("rk_T2: ", rk_T2, ", xk_T2: ", xk_T2, ", bm_T2: ", bm_T2, ", gm_T2: ", gm_T2)
        println("rk_T3: ", rk_T3, ", xk_T3: ", xk_T3, ", bm_T3: ", bm_T3, ", gm_T3: ", gm_T3)
      end

      if regelungEin
        tap = ResDataTypes.PowerTransformerTaps(tap_pos, tap_min, tap_max, tap_neutral, tap_step_percent, neutralU)
      end

      if tapSide == 1
        s1 = PowerTransformerWinding(vn_hv_kv, rk_T1, xk_T1, bm_T1, gm_T1, 0.0, neutralU, sn_hv_mva, tap)
        s2 = PowerTransformerWinding(vn_mv_kv, rk_T2, xk_T2, bm_T2, gm_T2, shift_mv_degree, nothing, sn_mv_mva, nothing)
        s3 = PowerTransformerWinding(vn_lv_kv, rk_T3, xk_T3, bm_T3, gm_T3, shift_lv_degree, nothing, sn_lv_mva, nothing)
      elseif tapSide == 2
        s1 = PowerTransformerWinding(vn_mv_kv, rk_T1, xk_T1, bm_T1, gm_T1, 0.0, nothing, sn_hv_mva, nothing)
        s2 = PowerTransformerWinding(vn_hv_kv, rk_T2, xk_T2, bm_T2, gm_T2, shift_mv_degree, neutralU, sn_mv_mva, tap)
        s3 = PowerTransformerWinding(vn_lv_kv, rk_T3, xk_T3, bm_T3, gm_T3, shift_lv_degree, nothing, sn_lv_mva, nothing)
      elseif tapSide == 3
        s1 = PowerTransformerWinding(vn_lv_kv, rk_T1, xk_T1, bm_T1, gm_T1, 0.0, nothing, sn_hv_mva, nothing)
        s2 = PowerTransformerWinding(vn_mv_kv, rk_T2, xk_T2, bm_T2, gm_T2, shift_mv_degree, nothing, sn_mv_mva, nothing)
        s3 = PowerTransformerWinding(vn_hv_kv, rk_T3, xk_T3, bm_T3, gm_T3, shift_lv_degree, neutralU, sn_lv_mva, tap)
      else
        s1 = PowerTransformerWinding(vn_hv_kv, rk_T1, xk_T1, bm_T1, gm_T1, 0.0, vn_hv_kv, sn_hv_mva, nothing)
        s2 = PowerTransformerWinding(vn_mv_kv, rk_T2, xk_T2, bm_T2, gm_T2, shift_mv_degree, vn_mv_kv, sn_mv_mva, nothing)
        s3 = PowerTransformerWinding(vn_lv_kv, rk_T3, xk_T3, bm_T3, gm_T3, shift_lv_degree, vn_lv_kv, sn_lv_mva, nothing)
      end

      t1_cmp = ResDataTypes.Component(cID, cName, "POWERTRANSFORMER", HV_kv)
      trafo = PowerTransformer(t1_cmp, regelungEin, s1, s2, s3)
      if log
        @show trafo
      end
      push!(trafos, trafo)
      AuxBusName = ""
      AuxBusIdx = 0
      AuxBusID = ""
      vn_aux_kv = 0.0
      for b in busVec
        if b.wt3id == cID
          AuxBusIdx = b.busIdx
          AuxBusName = b.name
          AuxBusID = b.id
          vn_aux_kv = b.vn_kv
          break
        end
      end
      auxBuxCmp = ResDataTypes.Component(AuxBusID, AuxBusName, "POWERTRANSFORMER", vn_aux_kv)
      if log
        println("AuxBusIdx: ", AuxBusIdx, ", AuxBusID: ", AuxBusID, ", vn_aux_kv: ", vn_aux_kv, ", AuxBusName: ", AuxBusName)
        @show auxBuxCmp
      end

      # create terminals for T1, HV_BUS -> AuxBus        
      terminal_T1_aux = ResDataTypes.Terminal(auxBuxCmp, ResDataTypes.Seite2)
      terminal_T1_HV = ResDataTypes.Terminal(t1_cmp, ResDataTypes.Seite1)
      # search for terminals of HV_BUS              
      t1Terminal = NodeTerminalsDict[hv_bus]
      AuxTerminal = NodeTerminalsDict[AuxBusIdx]
      # push terminals to terminal vector
      push!(t1Terminal, terminal_T1_HV)
      push!(AuxTerminal, terminal_T1_aux)
      r_pu, x_pu, b_pu, g_pu = calcTwoPortPU(HV_kv, baseMVA, rk_T1, xk_T1, bm_T1, gm_T1)
      if tapSide != 1
        tapPos = tap_neutral
      else
        tapPos = tap_pos
      end
      T1_ratio = calcRatio(HV_kv, vn_hv_kv, HV_kv, vn_hv_kv, tapPos, tap_neutral, tap_step_percent, 1)
      branch = ResDataTypes.Branch(t1_cmp, hv_bus, AuxBusIdx, cID, AuxBusID, r_pu, x_pu, b_pu, g_pu, T1_ratio, shift_degree, in_service)
      if log
        @show branch
      end
      push!(branchVec, branch)

      # create terminals for T2, AuxBux -> MV_BUS

      t2Terminal = NodeTerminalsDict[mv_bus]
      t2_cmp = ResDataTypes.Component(cID, cName, "POWERTRANSFORMER", MV_kv)
      terminal_T2_aux = ResDataTypes.Terminal(auxBuxCmp, ResDataTypes.Seite1)
      terminal_T2_mvbus = ResDataTypes.Terminal(t2_cmp, ResDataTypes.Seite2)
      push!(AuxTerminal, terminal_T2_aux)
      push!(t2Terminal, terminal_T2_mvbus)
      r_pu, x_pu, b_pu, g_pu = calcTwoPortPU(MV_kv, baseMVA, rk_T2, xk_T2, bm_T2, gm_T2)
      if tapSide != 2
        tapPos = tap_neutral
      else
        tapPos = tap_pos
      end

      T2_ratio = calcRatio(HV_kv, vn_aux_kv, MV_kv, vn_mv_kv, tapPos, tap_neutral, tap_step_percent, 1)
      mvID = NodeIDDict[mv_bus]
      branch = ResDataTypes.Branch(t2_cmp, AuxBusIdx, mv_bus, AuxBusID, mvID, r_pu, x_pu, b_pu, g_pu, T2_ratio, shift_degree, in_service)
      if log
        @show branch
      end
      push!(branchVec, branch)

      # create terminals for T3, AuxBux -> LV_BUS

      t3Terminal = NodeTerminalsDict[lv_bus]
      t3_cmp = ResDataTypes.Component(cID, cName, "POWERTRANSFORMER", LV_kv)
      terminal_T3_aux = ResDataTypes.Terminal(auxBuxCmp, ResDataTypes.Seite1)
      terminal_T3_lvbus = ResDataTypes.Terminal(t3_cmp, ResDataTypes.Seite2)
      push!(AuxTerminal, terminal_T3_aux)
      push!(t3Terminal, terminal_T3_lvbus)
      r_pu, x_pu, b_pu, g_pu = calcTwoPortPU(LV_kv, baseMVA, rk_T3, xk_T3, bm_T3, gm_T3)
      if tapSide != 3
        tapPos = tap_neutral
      else
        tapPos = tap_pos
      end

      T3_ratio = calcRatio(HV_kv, vn_aux_kv, LV_kv, vn_lv_kv, tapPos, tap_neutral, tap_step_percent, 1)
      lvID = NodeIDDict[lv_bus]
      branch = ResDataTypes.Branch(t3_cmp, AuxBusIdx, lv_bus, AuxBusID, lvID, r_pu, x_pu, b_pu, g_pu, T3_ratio, shift_degree, in_service)
      if log
        @show branch
      end
      push!(branchVec, branch)
    end
  end

  if !isnothing(loads)
    log_println("\nLoads:")

    for l in loads
      Bus = l["bus"]
      checkBusNumber(Bus, busVec)
      Vn = VoltageDict[Bus]

      cName = string(l["name"])
      cID = string(l["id"])
      nID = NodeIDDict[Bus]

      p = l["p_mw"]
      q = l["q_mvar"]
      ratedS = nothing
      ratedU = nothing
      qPercent = nothing
      maxP = nothing
      minP = nothing
      maxQ = nothing
      minQ = nothing
      ratedPowerFactor = nothing
      referencePri = nothing

      comp = ResDataTypes.Component(cID, cName, "LOAD", Vn)
      pRS = ProSumer(comp, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri, nothing, nothing)
      if log
        @show pRS
      end
      push!(prosum, pRS)

      t1 = ResDataTypes.Terminal(comp, ResDataTypes.Seite1)
      t1Terminal = NodeTerminalsDict[Bus]
      push!(t1Terminal, t1)

      nParms = NodeParametersDict[Bus]
      nParms.pƩLoad = isnothing(nParms.pƩLoad) ? p : nParms.pƩLoad + p
      nParms.qƩLoad = isnothing(nParms.qƩLoad) ? q : nParms.qƩLoad + q
      NodeParametersDict[Bus] = nParms
      if log
        println("nParms: ", nParms)
      end
    end
  end

  # static Generators     
  if !isnothing(sgens)
    log_println("\n\nStatic Generators:")
    for g in sgens
      Bus = g["bus"]
      checkBusNumber(Bus, busVec)
      Vn = VoltageDict[Bus]

      cName = string(g["name"])
      cID = string(g["id"])
      nID = NodeIDDict[Bus] # String

      p = float(g["p_mw"])
      q = float(g["q_mvar"])

      ratedS = nothing
      ratedU = nothing
      qPercent = nothing
      maxP = nothing
      minP = nothing
      maxQ = nothing
      minQ = nothing
      ratedPowerFactor = nothing
      referencePri = nothing

      comp = ResDataTypes.Component(cID, cName, "GENERATOR", Vn)

      pRS = ProSumer(comp, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri, nothing, nothing)
      if log
        @show pRS
      end
      push!(prosum, pRS)

      t2 = ResDataTypes.Terminal(comp, ResDataTypes.Seite2)
      t2Terminal = NodeTerminalsDict[Bus]
      push!(t2Terminal, t2)

      nParms = NodeParametersDict[Bus]
      nParms.pƩGen = p
      nParms.qƩGen = q
      NodeParametersDict[Bus] = nParms
      if log
        println("nParms: ", nParms)
      end
    end
  end

  # Voltage controlled Generators    
  if !isnothing(vgens)
    log_println("\nVoltage controlled Generators:")
    for g in vgens
      vm_pu = float(g["vm_pu"])
      Bus = g["bus"]
      checkBusNumber(Bus, busVec)
      Vn = VoltageDict[Bus]
      cName = string(g["name"])
      cID = string(g["id"])
      nID = NodeIDDict[Bus]
      p = float(g["p_mw"])
      q = nothing
      ratedS = nothing
      ratedU = nothing
      qPercent = nothing
      maxP = nothing
      minP = nothing
      maxQ = float(g["max_q_mvar"])
      minQ = float(g["min_q_mvar"])
      ratedPowerFactor = nothing
      referencePri = nothing
      comp = ResDataTypes.Component(cID, cName, "SYNCHRONOUSMACHINE", Vn)
      pRS = ProSumer(comp, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri, vm_pu, nothing)
      if log
        @show pRS
      end
      push!(prosum, pRS)

      t2 = ResDataTypes.Terminal(comp, ResDataTypes.Seite2)
      t2Terminal = NodeTerminalsDict[Bus]
      push!(t2Terminal, t2)

      nParms = NodeParametersDict[Bus]
      nParms.vm_pu = vm_pu
      nParms.pƩGen = p
      NodeParametersDict[Bus] = nParms
      if log
        println("nParms: ", nParms)
      end
    end
  end
  if !isnothing(grids)
    log_println("\nGrids:")
    for g in grids
      #EXG: bus=1, name=ex1, id=, vm_pu=1.02, va_degree=50.0
      Bus = g["bus"]
      checkBusNumber(Bus, busVec)
      Vn = VoltageDict[Bus]
      cName = string(g["name"])
      cID = string(g["id"])
      nID = NodeIDDict[Bus]
      vm_pu = float(g["vm_pu"])
      vm_degree = float(g["va_degree"])
      p = nothing
      q = nothing
      ratedS = nothing
      ratedU = nothing
      qPercent = nothing
      maxP = nothing
      minP = nothing
      maxQ = nothing
      minQ = nothing
      ratedPowerFactor = nothing
      referencePri = nothing
      comp = ResDataTypes.Component(cID, cName, "EXTERNALNETWORKINJECTION", Vn)
      pRS = ProSumer(comp, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri, vm_pu, vm_degree)
      if log
        @show pRS
      end
      push!(prosum, pRS)

      t2 = ResDataTypes.Terminal(comp, ResDataTypes.Seite2)
      t2Terminal = NodeTerminalsDict[Bus]
      push!(t2Terminal, t2)

      nParms = NodeParametersDict[Bus]
      nParms.vm_pu = vm_pu
      nParms.va_deg = vm_degree
      NodeParametersDict[Bus] = nParms
      if log
        println("nParms: ", nParms)
      end
    end
  end
  if !isnothing(shunts)
    log_println("\nShunts:")
    for s in shunts
      # Attention: 
      # Bus is a used as a short cut to find the busnumber for calculating Y-Bus. 
      # When the bus numbers are rearranged, this number must also be adjusted.
      Bus = s["bus"]
      checkBusNumber(Bus, busVec)
      Vn = VoltageDict[Bus]
      vn_kv = (s["vn_kv"] === "") ? Vn : float(s["vn_kv"])

      cName = string(s["name"])
      cID = string(s["id"])
      nID = NodeIDDict[Bus]
      in_service = Int(s["in_service"])

      if in_service < 1
        state = 0.0
      else
        state = 1.0
      end

      p_shunt = float(s["p_mw"])
      q_shunt = float(s["q_mvar"])

      ratio = (Vn / vn_kv)^2
      y_pu = calcYShunt(p_shunt, q_shunt, ratio, baseMVA)
      comp = ResDataTypes.Component(cID, cName, "LINEARSHUNTCOMPENSATOR", vn_kv)
      sh = ResDataTypes.Shunt(comp, nID, Bus, p_shunt, q_shunt, y_pu, in_service)
      if log
        @show sh
      end
      push!(shuntVec, sh)
      t1 = ResDataTypes.Terminal(comp, ResDataTypes.Seite1)
      t1Terminal = NodeTerminalsDict[Bus]
      push!(t1Terminal, t1)

      nParms = NodeParametersDict[Bus]
      nParms.pShunt = isnothing(nParms.pShunt) ? p_shunt : nParms.pShunt + p_shunt
      nParms.qShunt = isnothing(nParms.qShunt) ? q_shunt : nParms.qShunt + q_shunt
      NodeParametersDict[Bus] = nParms
      if log
        println("nParms: ", nParms)
      end
    end
  end
  log_println("\ncreate network:")

  for b in busVec
    busIdx = b.busIdx
    busName = b.name
    nodeType = b.type
    cID = b.id
    Vn = b.vn_kv
    terminals = NodeTerminalsDict[busIdx]

    if log
      println("Bus: ", busName, " busIdx: ", busIdx, " ID: ", cID, " Voltage: ", Vn, "\nTerminals: ", terminals)
    end
    node = ResDataTypes.Node(cID, busName, Vn, terminals)

    nodeType = ResDataTypes.toNodeType(nodeType)
    setNodeType!(node, nodeType)
    setBusIdx!(node, busIdx)

    if nodeType == nodeType == ResDataTypes.Isolated
      setNodeIdx!(node, 4)
    else
      setNodeIdx!(node, busIdx)
    end

    nParms = NodeParametersDict[busIdx]
    setNodeParameters!(node, nParms)
    if log
      @show node
    end
    push!(nodeVec, node)
  end

  checkNodeConnections(nodeVec)

  net = ResDataTypes.Net(netName, baseMVA, slackIdx, nodeVec, ACLines, trafos, branchVec, prosum, shuntVec)

  return net
end


function createNetFromPGM(filename)
  @info "create network from PGM-File: $(filename)"
  
  nodes, lines, wt2, wt3, sym_gens, sym_loads, shunts, source = SparlectraImport.pgmparser(filename)
  #println("Nodes: ", nodes)
  #println("Source: ", source)
  base_name = basename(filename)
  netName, ext = splitext(base_name)
  @info "Netname: $netName"
  
  # baseMVA is not given, so we search for the maximum sn of the transformers
  sn_max = 0.0
  for trafo in wt2
    sn = float(trafo["sn"])
    if sn > sn_max
      sn_max = sn
    end
  end    
  baseMVA = sn_max/1e6
  @info "baseMVA: $baseMVA MVA"
    
  busVec = Vector{Bus}()
  VoltageDict = Dict{Integer,Float64}()
  NodeIDDict = Dict{Integer,String}()

  ACLines = Vector{ResDataTypes.ACLineSegment}()
  trafos = Vector{ResDataTypes.PowerTransformer}()
  prosum = Vector{ResDataTypes.ProSumer}()
  shuntVec = Vector{ResDataTypes.Shunt}()
  branchVec = Vector{ResDataTypes.Branch}()
  nodeVec = Vector{ResDataTypes.Node}()

  NodeTerminalsDict = Dict{Integer,Vector{Terminal}}()
  NodeParametersDict = Dict{Integer,NodeParameters}()
  bus_types = Dict{Int, Int}()
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
  for s in source
    bus_id = s["node"]
    bus_type = 3
    bus_types[bus_id] = bus_type
    break # only one slack bus
  end
  
  slackIdx = 0 # Counter for slack index
  for bus in nodes
    busIdx = Int64(bus["id"])    
    vn_kv = float(bus["u_rated"])/1000.0
    name = "Bus_"*string(busIdx)
    nodeID = "#ID_"*name
    if haskey(bus_types, busIdx)
      btype = bus_types[busIdx]
      if btype == 3
        slackIdx = busIdx
      end
    else      
      btype = 1
    end
    a_bus = Bus(busIdx, name, nodeID, "", vn_kv, btype)
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
  @assert slackIdx != 0 "no slack bus found!"
  for (i, bus) in enumerate(busVec)
    @assert bus.busIdx == i "bus numbers are not consecutive and unique"        
  end
  
  #for b in busVec
  #  println(b)
  #end

  @info "$(length(busVec)) busses created..."
  # Set up auxillary buses
  @assert wt3 === nothing "3WT transformers not supported yet"
  
  b = nothing
  g = nothing

  for line in lines
    from = Int64(line["from_node"])
    to = Int64(line["to_node"])
    
    num = line["id"]
    cName = "Line_"*string(from)*"_"*string(to)
    cID = "#ID_Line"*string(num)
    r = float(line["r1"])
    x = float(line["x1"])
    c_f = float(line["c1"])
    c_nf = c_f*1e9
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
    
    asec = ACLineSegment(cID, cName, Vn, 1.0, r, x, b, g, c_nf)
    push!(ACLines, asec)

    cLine = toComponentTyp("ACLINESEGMENT")
    c = Component(cID, cName, cLine, Vn)
    t1 = ResDataTypes.Terminal(c, ResDataTypes.Seite1)
    t2 = ResDataTypes.Terminal(c, ResDataTypes.Seite2)
    
    checkBusNumber(from, busVec)
    checkBusNumber(to, busVec)

    t2Terminal = NodeTerminalsDict[from]
    t1Terminal = NodeTerminalsDict[to]
    push!(t1Terminal, t1)
    push!(t2Terminal, t2)
    
    bch = c_f * 2.0 * 50.0 * pi
    g_us = 0.0 # should be calculated later ....
    
    r_pu, x_pu, b_pu, g_pu = calcTwoPortPU(Vn, baseMVA, r , x , bch, g_us)
    fromNodeID = NodeIDDict[from]
    toNodeID = NodeIDDict[to]
    from_org = from
    to_org = to
    branch = Branch(c, from, to, from_org, to_org, fromNodeID, toNodeID, r_pu, x_pu, b_pu, g_pu, 0.0, 0.0, inService)
    push!(branchVec, branch)
  end
  @info "$(length(ACLines)) aclines created..."

  for t in wt2
    sn = float(t["sn"])
    u1 = float(t["u1"])
    u2 = float(t["u2"])
    
    from_kV = u1/1000.0 # rated voltage at the from side
    to_kV = u2/1000.0 # rated voltage at the to side
    if from_kV > to_kV
      vn_hv = from_kV
      vn_lv = to_kV
    else
      vn_hv = to_kV
      vn_lv = from_kV
    end

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
    cName = "Trafo_"*string(from_node)*"_"*string(to_node)
    cID = "#ID_Trafo"*string(t["id"])
    c = Component(cID, cName, "POWERTRANSFORMER", vn_hv)
    t1 = Terminal(c, ResDataTypes.Seite1)
    t2 = Terminal(c, ResDataTypes.Seite2)
    checkBusNumber(from_node, busVec)
    checkBusNumber(to_node, busVec)

    
    tap_side = Int64(t["tap_side"])
    tap_pos = Int64(t["tap_pos"])
    tap_min = Int64(t["tap_min"])
    tap_max = Int64(t["tap_max"])
    tap_nom = Int64(t["tap_nom"])
    tap_size_V = float(t["tap_size"])
    #side_1 = 0, side_2 = 1, side_3 = 2
    tapSeite = (tap_side == 0) ? ((from_kV > to_kV) ? 1 : 2) :
           (tap_side == 1) ? ((from_kV > to_kV) ? 2 : 1) : @error "tap_side: $tap_side, from_kV: $from_kV, to_kV: $to_kV"
    
    
    tap_neutral = div(tap_max + tap_min,2)
    HV = VoltageDict[from_node]
    LV = VoltageDict[to_node]
    
    if tap_side == 0
      tapStepPercent=calcTapStepPercent(tap_size_V, u1)
    else
      tapStepPercent=calcTapStepPercent(tap_size_V, u2)
    end    
    shift_degree = 0.0
    sn_MVA = sn/1e6
    ratio = calcRatio(HV,vn_hv,LV,vn_lv,tap_pos,tap_neutral,tapStepPercent, tapSeite)
    if tapSeite == 1
      Vtab = calcNeutralU(1.0, vn_hv, tap_min, tap_max, tapStepPercent)
      tap = PowerTransformerTaps(tap_pos, tap_min, tap_max, tap_neutral, tapStepPercent, Vtab)
      r,x,b,g = calcTrafoParamsSI(sn_max, u2, uk, sn, pk_W, i0, p0_W )       
      r_pu, x_pu, b_pu, g_pu = calcTwoPortPU(vn_hv, baseMVA, r, x, b, g)
          
      s1 = PowerTransformerWinding(vn_hv, r, x, b, g, shift_degree, from_kV, sn_MVA, tap)
      s2 = PowerTransformerWinding(vn_lv, 0.0, 0.0)
    else
      Vtab = calcNeutralU(neutralU_ratio, vn_lv, tap_min, tap_max, tapStepPercent)
      tap = PowerTransformerTaps(tap_pos, tap_min, tap_max, tap_neutral, tapStepPercent, Vtab)
      r,x,b,g = calcTrafoParamsSI(sn_max, u1, uk, sn, pk_W, i0, p0_W ) 

      s1 = PowerTransformerWinding(vn_hv, 0.0, 0.0)
      s2 = PowerTransformerWinding(vn_lv, r, x, b, g, shift_degree, to_kV, sn_MVA, tap)
    end
    #@debug "Trafo: ", cName, ", r: ", r, ", x: ", x, ", b: ", b, ", g: ", g
    #@debug "Trafo: ", cName, ", r_pu: ", r_pu, ", x_pu: ", x_pu, ", b_pu: ", b_pu, ", g_pu: ", g_pu
    s3 = nothing
    cmp = Component(cID, cName, "POWERTRANSFORMER", vn_hv)
    trafo = PowerTransformer(cmp, true, s1, s2, s3)
    push!(trafos, trafo)
    #@show trafo

    t1 = ResDataTypes.Terminal(cmp, ResDataTypes.Seite1)
    cmp2 = ResDataTypes.Component(cID, cName, "POWERTRANSFORMER", vn_lv)
    t2 = ResDataTypes.Terminal(cmp2, ResDataTypes.Seite2)

    t1Terminal = NodeTerminalsDict[from_node]
    t2Terminal = NodeTerminalsDict[to_node]

    fromNodeID = NodeIDDict[from_node]
    toNodeID = NodeIDDict[to_node]

    push!(t1Terminal, t1)
    push!(t2Terminal, t2)

    branch = Branch(cmp, from_node, to_node, fromNodeID, toNodeID, r_pu, x_pu, b_pu, g_pu, ratio, shift_degree, inService)

    push!(branchVec, branch)


  end
  @info "$(length(trafos)) power transformers created..."
  
  
  for sym_gen in sym_gens
    bus_id = sym_gen["node"]
    vn = VoltageDict[bus_id]
  end
  
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
    status=sym_load["status"]
    if status == 0
      continue
    end
    lanz+=1
    id = sym_load["id"]
    bus = sym_load["node"]
    vn = VoltageDict[bus]
    cName = "Load_"*string(id)
    cID = "#ID_Load_"*string(id)
    nID = NodeIDDict[bus] # String
    p = sym_load["p_specified"]
    q = sym_load["q_specified"]

    comp = Component(cID, cName, "LOAD", vn)
    pRS = ProSumer(comp, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri, nothing, nothing)
    push!(prosum, pRS)
    t1 = Terminal(comp, ResDataTypes.Seite1)
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
    status=sym_gen["status"]
    if status == 0
      continue
    end
    lanz+=1
    id = sym_gen["id"]
    bus = sym_gen["node"]
    vn = VoltageDict[bus]
    cName = "Gen_"*string(id)
    cID = "#ID_Gen_"*string(id)
    nID = NodeIDDict[bus] # String
    p = sym_gen["p_specified"]
    q = sym_gen["q_specified"]
    comp = Component(cID, cName, "GENERATOR", vn)
    pRS = ProSumer(comp, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri, nothing, nothing)
    push!(prosum, pRS)
    t2 = Terminal(comp, ResDataTypes.Seite1)
    t2Terminal = NodeTerminalsDict[bus]
    push!(t2Terminal, t2)

    nParms = NodeParametersDict[bus]
    nParms.pƩGen = p
    nParms.qƩGen = q
    NodeParametersDict[bus] = nParms
  end
  @info "$(lanz) generators created..."
  
  lanz = 0
  for shunt in shunts
    status = shunt["status"]
    if status == 0
      continue
    end
    lanz+=1
    bus = shunt["node"]
    id = shunt["id"]
    vn = VoltageDict[bus]    
    cName = "Shunt_"*string(id)
    cID = "#ID_Shunt_"*string(id)
    g1 = float(shunt["g1"])
    b1 = float(shunt["b1"])
    comp = Component(cID, cName, "LINEARSHUNTCOMPENSATOR", vn)
    Y = Complex(g1, b1)
    y_pu = calc_y_pu(Y, vn, baseMVA)
    comp = Component(cID, cName, "LINEARSHUNTCOMPENSATOR", vn)
    sh = Shunt(comp, NodeIDDict[bus], bus, g1, b1, y_pu, status)
    push!(shuntVec, sh)
    t1 = Terminal(comp, ResDataTypes.Seite1)
    t1Terminal = NodeTerminalsDict[bus]
    push!(t1Terminal, t1)
    # p_shunt and q_shunt are not used in the shunt model
  end
  @info "$(lanz) shunts created..."  

  for b in busVec
    busIdx = b.busIdx
    busName = b.name
    nodeType = b.type
    cID = b.id
    Vn = b.vn_kv
    terminals = NodeTerminalsDict[busIdx]

    
    println("Bus: ", busName, " busIdx: ", busIdx, " ID: ", cID, " Voltage: ", Vn, "\nTerminals: ", terminals)
    
    node = Node(cID, busName, Vn, terminals)

    nodeType = toNodeType(nodeType)
    setNodeType!(node, nodeType)
    setBusIdx!(node, busIdx)

    if nodeType == ResDataTypes.Isolated
      setNodeIdx!(node, 4)
    else
      setNodeIdx!(node, busIdx)
    end

    nParms = NodeParametersDict[busIdx]
    setNodeParameters!(node, nParms)
    
    #@show node
    
    push!(nodeVec, node)
  end
  @warn "TODO: check Seite1/Seite2, check branchVec, check Power Injektion (to big cause of SI), parallel branches, ..."
  @info "3WT, not implemented yet"
  checkNodeConnections(nodeVec)

  net = ResDataTypes.Net(netName, baseMVA, slackIdx, nodeVec, ACLines, trafos, branchVec, prosum, shuntVec)

  return net





end

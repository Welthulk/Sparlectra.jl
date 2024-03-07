# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 05.03.2024
# experimental arae

# outdated functions:
function setParallelBranches!(branches::Vector{Branch})
  branchTupleSet = Set{Tuple}()
  branchDict = Dict{Tuple{Integer,Integer},Vector{Branch}}()

  for b in branches
    tupple = (b.fromBus, b.toBus)

    if tupple in branchTupleSet
      existing_branches = branchDict[tupple]
      push!(existing_branches, b)
    else
      branchDict[tupple] = [b]
      push!(branchTupleSet, tupple)
    end
  end

  for (k, b_vec) in branchDict
    if length(b_vec) > 1
      sum_b_pu = 0.0
      sum_g_pu = 0.0
      sum_z = 0.0

      for b in b_vec
        b.isParallel = true
        b.skipYBus = true
        if b.status == 1          
          sum_b_pu += b.b_pu
          sum_g_pu += b.g_pu
          sum_z += (b.r_pu - b.x_pu * im) / (b.r_pu^2 + b.x_pu^2)
        end
      end

      z_total = 1.0 / sum_z
      r_pu = real(z_total)
      x_pu = imag(z_total)

      last_b = b_vec[end]
      for (i,b) in enumerate(b_vec)
        if b.status == 1
          adjP = AdjElecParams(r_pu = r_pu, x_pu = x_pu, b_pu = sum_b_pu, g_pu = sum_g_pu)
          setAdjElecParam!(adjP, b)
          @debug "branch (2): $(b)"
          if i==1
            @debug "last parallel element (3): $(b)"
            b.skipYBus = false
          end
        end
      end
    end
  end
end


function createNetFromFile(filename, base_MVA::Float64 = 0.0, log::Bool = false, check = false)::ResDataTypes.Net
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

  slack_vm_pu = 1.0
  slack_va_deg = 0.0
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
    @assert bus.busIdx == i "bus numbers are not consecutive and unique"
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
    @assert length > 0 "line length must be greater than 0"
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

    from_bus = line["from_bus"]
    to_bus = line["to_bus"]
    checkBusNumber(from_bus, busVec)
    checkBusNumber(to_bus, busVec)

    cPGM = ImpPGMComp(cID, cName, toComponentTyp("ACLINESEGMENT"), Vn, from_bus, to_bus)
    asec = ACLineSegment(cPGM, length, r, x, b, g, c_nf_per_km, 0.0)
    if log
      println(asec)
    end
    push!(ACLines, asec)

    t1 = ResDataTypes.Terminal(cPGM, ResDataTypes.Seite1)
    t2 = ResDataTypes.Terminal(cPGM, ResDataTypes.Seite2)

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
    branch = ResDataTypes.Branch(cPGM, from_bus, to_bus, from_bus, to_bus, fromNodeID, toNodeID, r_pu, x_pu, b_pu, g_pu, ratio, angle, status)
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
    org_vk_percent = vk_percent
    vkr_percent = nothing
    pk_kW = nothing
    try
      vkr_percent = float(t["vkr_percent"])
      pk_kW = float(t["pk_kw"])
    catch
    end
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

    vk_dependence = nothing
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
    if !isnothing(vk_dependence)
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
    
    mParm=nothing    
    if !isnothing(pk_kW)             
      vkp=org_vk_percent/100.0
      i0p=io_percent
      p=pk_kW*1e3
      s=sn*1e6
      pfe=pfe_kw*1e3
      mParm=TransformerModelParameters(s, vkp, p, i0p, pfe)
    end

    if regelungEin && tapSeite == 1
      Vtab = calcNeutralU(neutralU_ratio, vn_hv, tap_min, tap_max, tap_step_percent)
      tap = ResDataTypes.PowerTransformerTaps(tap_pos, tap_min, tap_max, tap_neutral, tap_step_percent, Vtab)
      
      r_hv, x_hv, b_hv, g_hv, = calcTrafoParams(sn_mva=sn, vn_hv_kv=vn_hv, vk_percent=vk_percent,  pk_kw=pk_kW, vkr_percent=vkr_percent,  pfe_kw=pfe_kw, i0_percent=io_percent)
      
      r_pu, x_pu, b_pu, g_pu = calcTwoPortPU(vn_hv, sn, r_hv, x_hv, b_hv, g_hv)
      #function PowerTransformerWinding(; Vn::Float64, r::Float64, x::Float64, b::Union{Nothing, Float64} = nothing, g::Union{Nothing, Float64} = nothing, shift_degree::Union{Nothing, Float64} = nothing, ratedU::Union{Nothing, Float64} = nothing, ratedS::Union{Nothing, Float64} = nothing, taps::Union{Nothing, PowerTransformerTaps} = nothing, isPu_RXGB::Union{Nothing, Bool} = nothing, modelData::Union{Nothing, TransformerModelParameters} = nothing)
      #  new(Vn, r, x, b, g, shift_degree, ratedU, ratedS, taps, isPu_RXGB, modelData)
      s1 = PowerTransformerWinding(Vn = vn_hv, r = r_hv, x = x_hv, b = b_hv, g = g_hv, shift_degree = shift_degree, ratedU = vn_hv, ratedS = sn, taps = tap, isPu_RXGB = false, modelData = mParm)
      s2 = PowerTransformerWinding(Vn = vn_lv, r = 0.0, x = 0.0)

    elseif regelungEin && tapSeite == 2
      Vtab = calcNeutralU(neutralU_ratio, vn_lv, tap_min, tap_max, tap_step_percent)
      tap = ResDataTypes.PowerTransformerTaps(tap_pos, tap_min, tap_max, tap_neutral, tap_step_percent, Vtab)
      # maibe an error here: vn_hv?
      r_lv, x_lv, b_lv, g_lv = calcTrafoParams(sn_mva=sn, vn_hv_kv=vn_hv, vk_percent=vk_percent, pk_kw=pk_kw, vkr_percent=vkr_percent, pfe_kw=pfe_kw, i0_percent=io_percent)
      r_pu, x_pu, b_pu, g_pu = calcTwoPortPU(vn_hv, sn, r_lv, x_lv, b_lv, g_lv)
      s1 = PowerTransformerWinding(Vn = vn_hv, r = 0.0, x = 0.0)
      s2 = PowerTransformerWinding(Vn = vn_lv, r = r_lv, x = x_lv, b = b_lv, g = g_lv, shift_degree = shift_degree, ratedU = vn_lv, ratedS = sn, taps = tap, isPu_RXGB = false, modelData = mParm)

    else
      r_hv, x_hv, b_hv, g_hv = calcTrafoParams(sn_mva=sn, vn_hv_kv=vn_hv, vk_percent=vk_percent, pk_kw=pk_kw, vkr_percent=vkr_percent,  pfe_kw=pfe_kw, i0_percent=io_percent)
      r_pu, x_pu, b_pu, g_pu = calcTwoPortPU(vn_hv, sn, r_hv, x_hv, b_hv, g_hv)
      s1 = PowerTransformerWinding(Vn = vn_hv, r = r_hv, x = x_hv, g = g_hv, b = b_hv, ratedU = vn_hv, ratedS = sn, isPu_RXGB = false, modelData = mParm)
      s2 = PowerTransformerWinding(Vn = vn_lv, r = 0.0, x = 0.0)
    end
    s3 = nothing

    cImpPGMComp = ImpPGMComp(cID, cName, toComponentTyp("POWERTRANSFORMER"), vn_hv, busIdx, USIdx)
    trafo = PowerTransformer(cImpPGMComp, regelungEin, s1, s2, s3, ResDataTypes.Ratio)
    if log
      @show trafo
    end
    push!(trafos, trafo)

    t1 = ResDataTypes.Terminal(cImpPGMComp, ResDataTypes.Seite1)
    cmp2 = ImpPGMComp(cID, cName, toComponentTyp("POWERTRANSFORMER"), vn_lv, busIdx, USIdx)
    t2 = ResDataTypes.Terminal(cmp2, ResDataTypes.Seite2)

    t1Terminal = NodeTerminalsDict[busIdx]
    t2Terminal = NodeTerminalsDict[USIdx]

    toNodeID = NodeIDDict[USIdx]

    push!(t1Terminal, t1)
    push!(t2Terminal, t2)

    branch = ResDataTypes.Branch(cImpPGMComp, busIdx, USIdx, nID, toNodeID, r_pu, x_pu, b_pu, g_pu, ratio, shift_degree, inService)
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
      
      vkcorr = 1.0
      if vk_dependence != raw"" && !isnothing(vk_dependence)
        for dep in vkDepChr
          if vk_dependence in keys(vkDepChr)
            dep = vkDepChr[vk_dependence]

            tapVec = dep["tap"]
            vkVec = dep["vk"]

            xTaps = [Int64(x) for x in tapVec]
            yVKs = [Float64(y) for y in vkVec]

            vkcorr = calcVKDependence(xTaps, yVKs, Float64(tap_pos))

            
            
          end
        end
      end

      tap = nothing
      neutralU = 0.0
      ntap = calcTapCorr(tap_pos, tap_neutral, tap_step_percent)

      if tapSide == 1
        neutralU = vn_hv_kv
        vn_hv_kv = neutralU * ntap        
        vk_hv_percent = vkcorr        
      elseif tapSide == 2
        neutralU = vn_mv_kv
        vn_mv_kv = neutralU * ntap        
        vk_mv_percent = vkcorr        
      elseif tapSide == 3
        neutralU = vn_lv_kv
        vn_lv_kv = neutralU * ntap        
        vk_lv_percent = vkcorr        
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
      comp = ImpPGMComp(cID, cName, toComponentTyp("LOAD"), Vn, Bus, Bus)
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

      vm_pu = nothing
      vm_degree = nothing

      isPUNode = false

      if Bus == slackIdx
        referencePri = slackIdx
        vm_pu = 1.0
        vm_degree = 0.0
        slack_vm_pu = vm_pu
        slack_va_deg = vm_degree
      end
      comp = ImpPGMComp(cID, cName, toComponentTyp("GENERATOR"), Vn, Bus, Bus)
      pRS = ProSumer(comp, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri, vm_pu, vm_degree, ResDataTypes.Injection, isPUNode)
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
      vm_degree = nothing
      isPUNode = false

      if Bus == slackIdx
        referencePri = slackIdx
        vm_degree = 0.0
        slack_vm_pu = vm_pu
        slack_va_deg = vm_degree
        isPUNode = true
      end
      comp = ImpPGMComp(cID, cName, toComponentTyp("SYNCHRONOUSMACHINE"), Vn, Bus, Bus)
      pRS = ProSumer(comp, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri, vm_pu, vm_degree, ResDataTypes.Injection, isPUNode)
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
      isPUNode = true
      if Bus == slackIdx
        referencePri = slackIdx
        slack_vm_pu = vm_pu
        slack_va_deg = vm_degree
      end

      comp = ImpPGMComp(cID, cName, toComponentTyp("EXTERNALNETWORKINJECTION"), Vn, Bus, Bus)
      pRS = ProSumer(comp, nID, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri, vm_pu, vm_degree, ResDataTypes.Injection, isPUNode)
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

      comp = ImpPGMComp(cID, cName, toComponentTyp("LINEARSHUNTCOMPENSATOR"), vn_kv, Bus, Bus)
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

    c = ImpPGMComp(cID, busName, ResDataTypes.Busbarsection, Vn, busIdx, busIdx)
    node = Node(c, terminals, busIdx, busIdx, toNodeType(nodeType))
    nParms = NodeParametersDict[busIdx]

    if nodeType == 3
      nParms.vm_pu = slack_vm_pu
      nParms.va_deg = slack_va_deg
    end
    setNodeParameters!(node, nParms)

    if log
      @show node
    end
    push!(nodeVec, node)
  end

  if check
    checkNodeConnections(nodeVec)
  end

  setParallelBranches!(branchVec)

  net = ResDataTypes.Net(netName, baseMVA, slackIdx, nodeVec, ACLines, trafos, branchVec, prosum, shuntVec)

  return net
end

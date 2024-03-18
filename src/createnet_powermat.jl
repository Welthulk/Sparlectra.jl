# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 01.10.2023
# include-file createnet_powermat.jl

# get minimum degree ordering

function makeMDO!(busMapDict::Dict{Int,Int}, busData::Matrix{Float64}, branchData::Matrix{Float64}, busDict::Dict{String,Int64}, branchDict::Dict{String,Int64}, checkOnly::Bool)
  nodeIsolateSet = Set{Int}()
  nodeNumberSet = Set{Int}()
  branchTupleSet = Set{Tuple}()
  transDict = Dict{Int,Int}()
  rtransDict = Dict{Int,Int}()

  idx = 0
  col = size(busData, 2)
  for row in eachrow(busData[:, 1:col])
    btype = Int64(row[busDict["type"]])
    if btype == 4 # isolated      
      @debug "isolated bus found, skipped"
      continue
    end
    idx += 1

    transDict[Int64(row[busDict["bus"]])] = idx
    rtransDict[idx] = Int64(row[busDict["bus"]])

    push!(nodeNumberSet, idx)
    push!(nodeIsolateSet, idx)
  end

  col = size(branchData, 2)
  for row in eachrow(branchData[:, 1:col])
    status = Int64(row[branchDict["status"]])
    if status == 1
      _fbus = Int64(row[branchDict["fbus"]])
      _tbus = Int64(row[branchDict["tbus"]])

      fbus = transDict[_fbus]
      tbus = transDict[_tbus]

      push!(branchTupleSet, (fbus, tbus))

      delete!(nodeIsolateSet, fbus)
      delete!(nodeIsolateSet, tbus)
    end
  end

  if !isempty(nodeIsolateSet)
    @error "isolated nodes that are not marked as isolated found: $(nodeIsolateSet)"
  end

  if checkOnly
    for i in nodeNumberSet
      _bus = rtransDict[i]
      busMapDict[_bus] = i
    end
  else
    order = mdoRCM(length(nodeNumberSet), branchTupleSet)
    for i in order
      _bus = rtransDict[i]
      busMapDict[_bus] = i
    end

    @debug begin
      println("result minimum degree ordering: ")
      println(order, "\n")
      println("nodes : ", branchTupleSet)
      println("branch: ", nodeNumberSet)
      println("mdo: ", order)
    end
  end
end

# helper
function _createDict()
  busKeys = ["bus", "type", "Pd", "Qd", "Gs", "Bs", "area", "Vm", "Va", "baseKV", "zone", "Vmax", "Vmin"]
  busDict = Dict{String,Int}()

  for (idx, key) in enumerate(busKeys)
    busDict[key] = idx
  end

  genKeys = ["bus", "Pg", "Qg", "Qmax", "Qmin", "Vg", "mBase", "status", "Pmax", "Pmin", "Pc1", "Pc2", "Qc1min", "Qc1max", "Qc2min", "Qc2max", "ramp_agc", "ramp_10", "ramp_30", "ramp_q", "apf"]
  genDict = Dict{String,Int}()
  for (idx, key) in enumerate(genKeys)
    genDict[key] = idx
  end

  branchKey = ["fbus", "tbus", "r", "x", "b", "rateA", "rateB", "rateC", "ratio", "angle", "status", "angmin", "angmax"]
  branchDict = Dict{String,Int}()
  for (idx, key) in enumerate(branchKey)
    branchDict[key] = idx
  end

  return busDict, genDict, branchDict
end

# base_MVA = 0.0 for default value in case file, otherwise set to desired value
function createNetFromMatPowerFile(filename, base_MVA::Float64 = 0.0, log::Bool = false, mdo::Bool = true)::ResDataTypes.Net
  debug = false
  @debug debug = true
  mFak = 10.0 # approximation factor for load and shunt for qMax, qMin, pMax, pMin

  ACLines = Vector{ACLineSegment}()
  trafos = Vector{PowerTransformer}()
  prosum = Vector{ProSumer}()
  shuntVec = Vector{Shunt}()
  branchVec = Vector{Branch}()
  nodeVec = Vector{Node}()

  NodeDict = Dict{Integer,Node}()
  vnDict = Dict{Integer,Float64}()

  proSumDict = Dict{Integer,ProSumer}()
  busMapDict = Dict{Int,Int}()

  @info "create network from case-file: $(filename)"
  netName, baseMVA, busData, genData, branchData = casefileparser(filename)

  slackIdx = 0
  if base_MVA > 0.0
    baseMVA = base_MVA
  end
  # parsing the data
  if debug
    println("Parsing the data...")
    println("Netname: ", netName)
    println("BaseMVA: ", baseMVA, "\n")
    println("Buses: len=", size(busData))
    col = size(busData, 2)

    for row in eachrow(busData[:, 1:col])
      println(row)
    end
    println("\nGenerators: len=", size(genData))
    col = size(genData, 2)

    for row in eachrow(genData[:, 1:col])
      println(row)
    end
    println("\nBranches: len=", size(branchData))
    col = size(branchData, 2)

    for row in eachrow(branchData[:, 1:col])
      println(row)
    end
  end

  busDict, genDict, branchDict = _createDict()
  makeMDO!(busMapDict, busData, branchData, busDict, branchDict, mdo == false)

  col = size(busData, 2)
  numb = 0
  lastBus = 0
  for row in eachrow(busData[:, 1:col])
    btype = Int64(row[busDict["type"]])
    if btype == 4 # isolated      
      @info "isolated bus $(busIdx) found, skipped"
      continue
    end

    numb += 1
    kIdx = Int64(row[busDict["bus"]]) # original bus number    
    if kIdx <= lastBus
      msg = "bus numbers are not consecutive and unique, lastBus: $(lastBus), current bus: $(kIdx)"
      throw(msg)
    end
    lastBus = kIdx

    busIdx = busMapDict[kIdx]

    vn_kv = float(row[busDict["baseKV"]]) <= 0.0 ? (@warn("Warnung: vn_kv at bus $(kIdx) <= 0"); 1.0) : float(row[busDict["baseKV"]])
    va_deg = float(row[busDict["Va"]])
    vm_pu = float(row[busDict["Vm"]]) <= 0.0 ? (@warn("vm_pu at bus $(kIdx) <= 0"); 1.0) : float(row[busDict["Vm"]])
    pƩLoad = float(row[busDict["Pd"]]) < 0.0 ? (@info("pLoad at bus $(kIdx) p < 0"); float(row[busDict["Pd"]])) : float(row[busDict["Pd"]])
    qƩLoad = float(row[busDict["Qd"]]) < 0.0 ? (@info("qLoad at bus $(kIdx) q < 0"); float(row[busDict["Qd"]])) : float(row[busDict["Qd"]])

    pShunt = float(row[busDict["Gs"]]) < 0.0 ? (@info("pShunt at bus $(kIdx) p < 0"); float(row[busDict["Gs"]])) : float(row[busDict["Gs"]])
    qShunt = float(row[busDict["Bs"]]) < 0.0 ? (@info("qShunt at bus $(kIdx) q < 0"); float(row[busDict["Bs"]])) : float(row[busDict["Bs"]])

    ratedS = baseMVA
    zone = Int64(row[busDict["zone"]])
    area = Int64(row[busDict["area"]])

    node = Node(busIdx = busIdx, Vn_kV = vn_kv, nodeType = toNodeType(btype), ratedS = ratedS, zone = zone, area = area, vm_pu = vm_pu, va_deg = va_deg, pƩLoad = pƩLoad, qƩLoad = qƩLoad, pShunt = pShunt, qShunt = qShunt)

    if btype == 3 && slackIdx == 0
      slackIdx = busIdx
    elseif btype == 3
      msg = "more than one slack bus found!"
      throw(msg)
    end

    vnDict[busIdx] = vn_kv
    NodeDict[busIdx] = node
    push!(nodeVec, node)
    if pShunt != 0.0 || qShunt != 0.0
      shunt = Shunt(fromBus = busIdx, id = kIdx, base_MVA = baseMVA, Vn_kV_shunt = vn_kv, p_shunt = pShunt, q_shunt = qShunt)
      push!(shuntVec, shunt)
    end

    if pƩLoad != 0.0 || qƩLoad != 0.0
      qMax = abs(mFak * qƩLoad <= baseMVA) ? abs(mFak * qƩLoad) : baseMVA  # approximation!
      qMin = -qMax  # approximation!
      pMax = abs(mFak * pƩLoad) <= baseMVA ? abs(mFak * pƩLoad) : base_MVA # approximation!
      pMin = -pMax # approximation!
      referencePri = slackIdx == busIdx ? busIdx : nothing

      vm_degree = 0.0
      p = ProSumer(vn_kv = vn_kv, oID = kIdx, busIdx = busIdx, type = toProSumptionType("LOAD"), p = pƩLoad, q = qƩLoad, maxP = pMax, minP = pMin, maxQ = qMax, minQ = qMin, referencePri = referencePri, vm_pu = vm_pu, vm_degree = vm_degree)
      push!(prosum, p)
    end
  end# read bus

  # branches: 
  col = size(branchData, 2)
  for row in eachrow(branchData[:, 1:col])
    _fbus = Int64(row[branchDict["fbus"]])
    _tbus = Int64(row[branchDict["tbus"]])
    if haskey(busMapDict, _fbus) == false || haskey(busMapDict, _tbus) == false
      @info "branch $(cName) not in service"
      continue
    else
      fbus = busMapDict[_fbus]
      tbus = busMapDict[_tbus]
    end

    r_pu = float(row[branchDict["r"]])
    x_pu = float(row[branchDict["x"]])
    b_pu = float(row[branchDict["b"]])
    ratedS = float(row[branchDict["rateA"]])
    ratio = float(row[branchDict["ratio"]])
    angle = float(row[branchDict["angle"]])
    status = Int64(row[branchDict["status"]])
    vn_kv = vnDict[fbus]

    piModel = BranchModel(r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, g_pu = 0.0, ratio = ratio, angle = angle, sn_MVA = ratedS)
    b = Branch(vn_kV = vn_kv, baseMVA = baseMVA, from = fbus, to = tbus, branch = piModel, id = _fbus, status = status)
    push!(branchVec, b)

    l_km = 1.0
    if ratio == 0.0 # line
      z_base = (vn_kv^2) / baseMVA
      r = r_pu * z_base
      x = x_pu * z_base
      l_km = round((r / 0.1), digits = 1) <= 0.0 ? 1.0 : l_km    # approximation
      line = ACLineSegment(vn_kv = vn_kv, from = fbus, to = tbus, length = l_km, r = r, x = x)
      push!(ACLines, line)
    else # transformer
      vnh_kv = vnDict[fbus]
      vnl_kv = vnDict[tbus]
      c = getTrafoImpPGMComp(false, vnh_kv, fbus, tbus)
      w1 = w2 = nothing
      if ratedS > 0.0
        vnh_kv > vnl_kv ? (ratedU = vnh_kv)  : (ratedU = vnh_kv)
        uk, P_kW, i0, Pfe_kW = recalc_trafo_model_data(baseMVA = baseMVA, Sn_MVA = ratedS, ratedU_kV = ratedU, r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, isPerUnit = true)
        modelData = TransformerModelParameters(sn_MVA = ratedS, vk_percent = uk, pk_kW = P_kW, i0_percent = i0, p0_kW = Pfe_kW)
        vnh_kv > vnl_kv ? (w1 = PowerTransformerWinding(Vn_kV = vnh_kv, modelData = modelData); w2 = PowerTransformerWinding(Vn_kV = vnl_kv)) : (w1 = PowerTransformerWinding(Vn_kV = vnh_kv); w2 = PowerTransformerWinding(Vn_kV = vnl_kv, modelData = modelData))
      else
        vnh_kv > vnl_kv ? (w1 = PowerTransformerWinding(Vn_kV = vnh_kv); w2 = PowerTransformerWinding(Vn_kV = vnl_kv)) : (w1 = PowerTransformerWinding(Vn_kV = vnh_kv); w2 = PowerTransformerWinding(Vn_kV = vnl_kv))
      end
      tap = false
      tType = ResDataTypes.Ratio
      if angle != 0.0
        @info "found phase shifter in casefile $(netName) (fbus: $(fbus), tbus: $(tbus))"
        tType = ResDataTypes.PhaseShifter
      end
      trafo = PowerTransformer(c, tap, w1, w2, nothing, tType)
      push!(trafos, trafo)
    end
  end

  # Generators:   
  col = size(genData, 2)
  for row in eachrow(genData[:, 1:col])
    status = Int64(row[genDict["status"]])
    _bus = Int64(row[genDict["bus"]])
    if status < 1
      @info "generator $(_bus) not in service"
      continue
    end
    pGen = float(row[genDict["Pg"]]) < 0.0 ? (@info("pGen at bus $(_bus) p < 0"); float(row[genDict["Pg"]])) : float(row[genDict["Pg"]])
    qGen = float(row[genDict["Qg"]]) < 0.0 ? (@info("qGen at bus $(_bus) q < 0"); float(row[genDict["Qg"]])) : float(row[genDict["Qg"]])

    if abs(pGen) < 1e-6 && abs(qGen) < 1e-6
      node = NodeDict[_bus]
      if node._nodeType == ResDataTypes.PV || node._nodeType == ResDataTypes.PQ
        @info "generator $(_bus) has no power output, ignored"
        continue
      end
    end

    bus = busMapDict[_bus]
    n = NodeDict[bus]
    isPUNode = false
    if n._nodeType == ResDataTypes.PV
      isPUNode = true
    end

    qMax = float(row[genDict["Qmax"]])
    qMin = float(row[genDict["Qmin"]])
    pMax = float(row[genDict["Pmax"]])
    pMin = float(row[genDict["Pmin"]])

    vm_pu = float(row[genDict["Vg"]])
    mBase = float(row[genDict["mBase"]])

    referencePri = slackIdx == bus ? bus : nothing
    vm_degree = 0.0
    vn_kv = vnDict[bus]
    p = ProSumer(vn_kv = vn_kv, busIdx = bus, oID = _bus, type = toProSumptionType("GENERATOR"), p = pGen, q = qGen, maxP = pMax, minP = pMin, maxQ = qMax, minQ = qMin, referencePri = referencePri, vm_pu = vm_pu, isAPUNode = isPUNode)

    proSumDict[_bus] = p
    push!(prosum, p)

    # set generation power for node
    node = NodeDict[bus]
    addGenPower!(node = node, p = pGen, q = qGen)
    if node._vm_pu != vm_pu
      if log
        @info "node voltage mismatch: Bus $(bus): $(node._vm_pu), Gen $(bus): $(vm_pu)"
      end
      node._vm_pu = vm_pu
    end
  end# read Generators

  sort!(nodeVec, by = x -> x.busIdx)

  net = ResDataTypes.Net(netName, baseMVA, slackIdx, nodeVec, ACLines, trafos, branchVec, prosum, shuntVec)
  return net
end

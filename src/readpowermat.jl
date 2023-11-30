# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 01.10.2023
# include-file readpowermate.jl

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
      @info "isolated bus found, skipped"
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
  mFak = 4.0 # approximation factor for load and shunt for qMax, qMin, pMax, pMin
  elemDict = Dict{Tuple,Int}()
  function addElement!(fromBus, toBus, type, element, vector)
    push!(vector, element)
    index = length(vector)
    key = (fromBus, toBus, type)
    elemDict[key] = index
  end

  function getIndex(fromBus, toBus, type)
    key = (fromBus, toBus, type)
    index = nothing
    if haskey(elemDict, key)
      index = elemDict[key]
    end

    return index
  end

  ACLines = Vector{ResDataTypes.ACLineSegment}()
  trafos = Vector{ResDataTypes.PowerTransformer}()
  prosum = Vector{ResDataTypes.ProSumer}()
  shuntVec = Vector{ResDataTypes.Shunt}()
  branchVec = Vector{ResDataTypes.Branch}()
  nodeVec = Vector{ResDataTypes.Node}()

  NodeIDDict = Dict{Integer,String}()
  NodeDict = Dict{Integer,ResDataTypes.Node}()
  vnDict = Dict{Integer,Float64}()
  NodeTerminalsDict = Dict{Integer,Vector{ResDataTypes.Terminal}}()
  proSumDict = Dict{Integer,ResDataTypes.ProSumer}()
  busMapDict = Dict{Int,Int}()

  println("create network from case-file: $(filename)")
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
      msg= "bus numbers are not consecutive and unique, lastBus: $(lastBus), current bus: $(kIdx)"
      throw(msg) 
    end
    lastBus = kIdx

    busIdx = busMapDict[kIdx]
    cName = "Bus_" * string(kIdx)
    vn_kv = float(row[busDict["baseKV"]]) <= 0.0 ? (@warn("Warnung: vn_kv at bus $(kIdx) <= 0"); 1.0) : float(row[busDict["baseKV"]])
    va_deg = float(row[busDict["Va"]])
    vm_pu = float(row[busDict["Vm"]]) <= 0.0 ? (@warn("vm_pu at bus $(kIdx) <= 0"); 1.0) : float(row[busDict["Vm"]])
    pƩLoad = float(row[busDict["Pd"]]) < 0.0 ? (@info("pLoad at bus $(kIdx) p < 0"); float(row[busDict["Pd"]]) ) : float(row[busDict["Pd"]])
    qƩLoad = float(row[busDict["Qd"]]) < 0.0 ? (@info("qLoad at bus $(kIdx) q < 0"); float(row[busDict["Qd"]]) ) : float(row[busDict["Qd"]])

    pShunt = float(row[busDict["Gs"]]) < 0.0 ? (@info("pShunt at bus $(kIdx) p < 0"); float(row[busDict["Gs"]]) ) : float(row[busDict["Gs"]])
    qShunt = float(row[busDict["Bs"]]) < 0.0 ? (@info("qShunt at bus $(kIdx) q < 0"); float(row[busDict["Bs"]]) ) : float(row[busDict["Bs"]])

    pƩGen = 0.0
    qƩGen = 0.0
    ratedS = baseMVA
    zone = Int64(row[busDict["zone"]])
    area = Int64(row[busDict["area"]])

    cID = string(UUIDs.uuid4())

    cType = ResDataTypes.toComponentTyp("BUS")
    c = ResDataTypes.Component(cID, cName, cType, vn_kv)

    node = Node(c, Vector{ResDataTypes.Terminal}(), busIdx, kIdx, ResDataTypes.toNodeType(btype), nothing, ratedS, zone, area, vm_pu, va_deg, pƩLoad, qƩLoad, pShunt, qShunt, pƩGen, qƩGen)
    if btype == 3 && slackIdx == 0
      slackIdx = busIdx
    elseif btype == 3
      msg = "more than one slack bus found!"
      throw(msg)
    end

    vnDict[busIdx] = vn_kv
    NodeIDDict[busIdx] = cID
    NodeDict[busIdx] = node
    NodeTerminalsDict[busIdx] = Vector{ResDataTypes.Terminal}()
    push!(nodeVec, node)
    if pShunt != 0.0 || qShunt != 0.0
      cShunt = "Sh_" * string(kIdx)#*"_"*string(busIdx)
      c = ResDataTypes.Component(string(UUIDs.uuid4()), cShunt, ResDataTypes.toComponentTyp("SHUNT"), vn_kv)
      state = 1
      y_pu = calcYShunt(pShunt, qShunt, 1.0, baseMVA)
      shunt = ResDataTypes.Shunt(c, cID, busIdx, pShunt, qShunt, y_pu, state)
      push!(shuntVec, shunt)
    end

    if pƩLoad != 0.0 || qƩLoad != 0.0
      shName = "Ld_" * string(kIdx) * "_" * string(busIdx)
      shID = string(UUIDs.uuid4())
      c = ResDataTypes.Component(shID, shName, ResDataTypes.toComponentTyp("LOAD"), vn_kv)

      qMax = mFak * qƩLoad  # approximation!
      qMin = -mFak * qƩLoad # approximation!
      pMax = mFak * pƩLoad  # approximation!
      pMin = -mFak * pƩLoad # approximation!
      nodeId = cID
      ratedPowerFactor = nothing
      referencePri = slackIdx == busIdx ? busIdx : nothing

      vm_degree = 0.0
      qPercent = nothing

      Sber = sqrt(pƩLoad^2 + qƩLoad^2)
      ratedS = Sber
      ratedU = vn_kv
      p = ProSumer(c, nodeId, ratedS, ratedU, qPercent, pƩLoad, qƩLoad, pMin, pMax, qMin, qMax, ratedPowerFactor, referencePri, vm_pu, vm_degree)
      push!(prosum, p)
      t1 = ResDataTypes.Terminal(c, ResDataTypes.Seite1)
      t1Terminal = NodeTerminalsDict[busIdx]
      push!(t1Terminal, t1)
    end

  end# read bus

  # branches: 
  col = size(branchData, 2)
  for row in eachrow(branchData[:, 1:col])
    _fbus = Int64(row[branchDict["fbus"]])
    _tbus = Int64(row[branchDict["tbus"]])
    cName = "Branch_" * string(_fbus) * "_" * string(_tbus)


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
    ratio = float(row[branchDict["ratio"]])
    angle = float(row[branchDict["angle"]])
    status = Int64(row[branchDict["status"]])
    cID = string(UUIDs.uuid4())

    l_km = 1.0
    if status == 1
      cType = ResDataTypes.toComponentTyp("BRANCH")
      vn_kv = vnDict[fbus]
      c = ResDataTypes.Component(cID, cName, cType, vn_kv)
      isParallel = false
      fromNodeId = NodeIDDict[fbus]
      toNodeId = NodeIDDict[tbus]

      index = getIndex(fbus, tbus, cType)
      if index !== nothing
        b2 = branchVec[index]
        rges = 1.0 / (1.0 / b2.r_pu + 1.0 / r_pu)
        xges = 1.0 / (1.0 / b2.x_pu + 1.0 / x_pu)
        bges = b2.b_pu + b_pu

        b2.r_pu = rges
        b2.x_pu = xges
        b2.b_pu = bges

        b2.isParallel = true
        @info "branch $(cName) already in branchVec, update values: r_pu: $(rges), x_pu: $(xges), b_pu: $(bges)"
        if ratio == 0.0
          cType = ResDataTypes.toComponentTyp("ACLINESEGMENT")
          index = getIndex(fbus, tbus, cType)
          if index !== nothing
            l2 = ACLines[index]
            l2.r = rges
            l2.x = xges
            l2.b = bges
          end
        else
          cType = ResDataTypes.toComponentTyp("POWERTRANSFORMER")
          index = getIndex(fbus, tbus, cType)
          if index !== nothing
            t2 = trafos[index]
            t2.side1.r = rges
            t2.side1.x = xges
            t2.side1.b = bges
          end
        end
      else
        b = Branch(c, fbus, tbus, _fbus, _tbus, fromNodeId, toNodeId, r_pu, x_pu, b_pu, 0.0, ratio, angle, status, nothing, nothing, isParallel)
        addElement!(fbus, tbus, cType, b, branchVec)

        t1 = ResDataTypes.Terminal(c, ResDataTypes.Seite1)
        t2 = ResDataTypes.Terminal(c, ResDataTypes.Seite2)
        t1Terminal = NodeTerminalsDict[fbus]
        t2Terminal = NodeTerminalsDict[tbus]
        push!(t1Terminal, t1)
        push!(t2Terminal, t2)
        if ratio == 0.0
          cType = ResDataTypes.toComponentTyp("ACLINESEGMENT")
          lName = "Line_" * string(fbus) * "_" * string(tbus)
          lID = string(UUIDs.uuid4())
          cLine = ResDataTypes.Component(lID, lName, cType, vn_kv)
          z_base = (vn_kv^2) / baseMVA
          r = r_pu * z_base
          l_km = round((r / 0.1), digits = 1) <= 0.0 ? 1.0 : l_km    # approximation
          line = ACLineSegment(cLine, l_km, r_pu, x_pu, b_pu)
          addElement!(fbus, tbus, cType, line, ACLines)
        else
          cType = ResDataTypes.toComponentTyp("POWERTRANSFORMER")
          vnh_kv = vnDict[fbus]
          w1 = PowerTransformerWinding(vnh_kv, r_pu, x_pu, b_pu, nothing, angle, nothing, nothing, nothing)
          vnl_kv = vnDict[tbus]
          w2 = PowerTransformerWinding(vnl_kv, 0.0, 0.0, 0.0, nothing, 0.0, nothing, nothing, nothing)
          cTName = "Trafo_" * string(fbus) * "_" * string(tbus)
          cTID = string(UUIDs.uuid4())
          cTrafo = ResDataTypes.Component(cTID, cTName, cType, vn_kv)
          tap = false
          tType = ResDataTypes.Ratio
          if angle != 0.0
            @info "found phase shifter in casefile $(netName) (fbus: $(fbus), tbus: $(tbus))"
            tType = ResDataTypes.PhaseShifter
          end
          trafo = PowerTransformer(cTrafo, tap, w1, w2, nothing, tType)
          addElement!(fbus, tbus, cType, trafo, trafos)
        end
      end
    else
      @info "branch $(cName) not in service"
    end
  end

  # Generators:   
  col = size(genData, 2)
  for row in eachrow(genData[:, 1:col])
    _bus = Int64(row[genDict["bus"]])
    pGen = float(row[genDict["Pg"]]) < 0.0 ? (@info("pGen at bus $(_bus) p < 0"); float(row[genDict["Pg"]]) ) : float(row[genDict["Pg"]])
    qGen = float(row[genDict["Qg"]]) < 0.0 ? (@info("qGen at bus $(_bus) q < 0"); float(row[genDict["Qg"]]) ) : float(row[genDict["Qg"]])

    
    bus = busMapDict[_bus]
    cName = "Gen_" * string(_bus)

    qMax = float(row[genDict["Qmax"]])
    qMin = float(row[genDict["Qmin"]])
    vm_pu = float(row[genDict["Vg"]])
    mBase = float(row[genDict["mBase"]])
    status = Int64(row[genDict["status"]])
    pMax = float(row[genDict["Pmax"]])
    pMin = float(row[genDict["Pmin"]])

    cID = string(UUIDs.uuid4())
    nodeId = NodeIDDict[bus]
    ratedPowerFactor = nothing

    if status == 1
      if haskey(proSumDict, _bus)
        proSum = proSumDict[_bus]
        _maxP = proSum.maxP
        _minP = proSum.minP
        _maxQ = proSum.maxQ
        _minQ = proSum.minQ
        pMax = max(_maxP, pMax)
        pMin = min(_minP, pMin)
        qMax = max(_maxQ, qMax)
        qMin = min(_minQ, qMin)
        @info "generator $(cName) already in prosumption list, update values: pGen: $(pGen), qGen: $(qGen), pMax: $(pMax), pMin: $(pMin), qMax: $(qMax), qMin: $(qMin)"
      end

      referencePri = slackIdx == bus ? bus : nothing
      vm_degree = 0.0
      qPercent = nothing
      ratedS = mBase
      ratedU = vnDict[bus]
      vn_kv = vnDict[bus]
      c = ResDataTypes.Component(cID, cName, ResDataTypes.toComponentTyp("GENERATOR"), vn_kv)
      p = ProSumer(c, nodeId, ratedS, ratedU, qPercent, pGen, qGen, pMin, pMax, qMin, qMax, ratedPowerFactor, referencePri, vm_pu, vm_degree)
      proSumDict[_bus] = p
      push!(prosum, p)

      # set generation power for node
      node = NodeDict[bus]
      node._pƩGen = pGen
      node._qƩGen = qGen
      if node._vm_pu != vm_pu
        if log
          @info "node voltage mismatch: Bus $(bus): $(node._vm_pu), Gen $(bus): $(vm_pu)"
        end
        node._vm_pu = vm_pu
      end
      t1 = ResDataTypes.Terminal(c, ResDataTypes.Seite1)
      t1Terminal = NodeTerminalsDict[bus]
      push!(t1Terminal, t1)
    else
      @info "generator $(cName) not in service"
    end
  end# read Generators

  for t in NodeTerminalsDict
    busIdx = t[1]
    terminals = t[2]
    node = NodeDict[busIdx]
    node.terminals = terminals
  end

  sort!(nodeVec, by = x -> x.busIdx)
  @debug begin
    for n in nodeVec
      println("original bus number: $(n._kidx), new bus $(n.busIdx) number")
    end
  end

  net = ResDataTypes.Net(netName, baseMVA, slackIdx, nodeVec, ACLines, trafos, branchVec, prosum, shuntVec)
  return net
end




# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 01.10.2023
# include-file createnet_powermat.jl

# helper
#! format: off
function _creaDict()

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
#! format: on
"""
Creates a network from a MatPower case file.

# Arguments
- `filename`: Path to the MatPower case file.
- `log::Bool = false`: Whether to log information (default: false).

# Returns
A Net object representing the network.

"""

createNetFromMatPowerFile(filename::String, log::Bool = false, flatstart::Bool = false) = createNetFromMatPowerFile(; filename = filename, log = log, flatstart = flatstart)

function createNetFromMatPowerFile(; filename::String, log::Bool = false, flatstart::Bool = false)::Net
  function pInfo(msg::String)
    if log
      @info msg
    end
  end
  getIntBusStr(Float64) = string(Int(Float64))
  getFloatBusStr(Int64) = string(Int64)

  debug = false
  @debug debug = true
  mFak = 10.0 # approximation factor for load and shunt for qMax, qMin, pMax, pMin

  mpc        = MatpowerIO.read_case(filename; legacy_compat = true)
  netName    = mpc.name
  baseMVA    = mpc.baseMVA
  busData    = mpc.bus
  genData    = mpc.gen
  branchData = mpc.branch

  slackIdx = 0
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

  busDict, genDict, branchDict = _creaDict()
  flat = flatstart
  myNet = Net(name = string(netName), baseMVA = baseMVA, flatstart = flat)

  col = size(busData, 2)

  for row in eachrow(busData[:, 1:col])
    btype = Int64(row[busDict["type"]])
    kIdx = Int64(row[busDict["bus"]]) # original bus number    
    raw_vn = float(row[busDict["baseKV"]])
    if raw_vn <= 0.0
      log && @warn("Warnung: vn_kv at bus $(kIdx) <= 0")
      vn_kv = 1.0
    else
      vn_kv = raw_vn
    end

    va_deg = float(row[busDict["Va"]])

    vm_pu = float(row[busDict["Vm"]]) <= 0.0 ? (@warn("vm_pu at bus $(kIdx) <= 0"); 1.0) : float(row[busDict["Vm"]])
    pƩLoad = float(row[busDict["Pd"]]) < 0.0 ? (pInfo("pLoad at bus $(kIdx) p < 0"); float(row[busDict["Pd"]])) : float(row[busDict["Pd"]])
    qƩLoad = float(row[busDict["Qd"]]) < 0.0 ? (@debug("qLoad at bus $(kIdx) q < 0"); float(row[busDict["Qd"]])) : float(row[busDict["Qd"]])

    pShunt = float(row[busDict["Gs"]]) < 0.0 ? (pInfo("pShunt at bus $(kIdx) p < 0"); float(row[busDict["Gs"]])) : float(row[busDict["Gs"]])
    qShunt = float(row[busDict["Bs"]]) < 0.0 ? (@debug("qShunt at bus $(kIdx) q < 0"); float(row[busDict["Bs"]])) : float(row[busDict["Bs"]])

    zone = Int64(row[busDict["zone"]])
    area = Int64(row[busDict["area"]])
    vmin = float(row[busDict["Vmin"]])
    vmax = float(row[busDict["Vmax"]])

    if btype == 3 && slackIdx == 0
      slackIdx = kIdx
    end
    busName = string(kIdx)

    addBus!(net = myNet, busName = busName, busType = toString(toNodeType(btype)), vn_kV = vn_kv, vm_pu = vm_pu, va_deg = va_deg, vmin_pu = vmin, vmax_pu = vmax, isAux = false, oBusIdx = kIdx, zone = zone, area = area)

    if pShunt != 0.0 || qShunt != 0.0
      addShunt!(net = myNet, busName = busName, pShunt = pShunt, qShunt = qShunt)
    end

    if pƩLoad != 0.0 || qƩLoad != 0.0
      qMax = abs(mFak * qƩLoad <= baseMVA) ? abs(mFak * qƩLoad) : baseMVA  # approximation!
      qMin = -qMax  # approximation!
      pMax = abs(mFak * pƩLoad) <= baseMVA ? abs(mFak * pƩLoad) : baseMVA # approximation!
      pMin = -pMax # approximation!
      addProsumer!(net = myNet, busName = busName, type = "ENERGYCONSUMER", p = pƩLoad, q = qƩLoad, pMin = pMin, pMax = pMax, qMin = qMin, qMax = qMax)
    end
  end# read bus

  # branches: 
  col = size(branchData, 2)
  for row in eachrow(branchData[:, 1:col])
    _fbus = Int(row[branchDict["fbus"]])
    _tbus = Int(row[branchDict["tbus"]])

    fbus = string(_fbus)
    tbus = string(_tbus)
    vn_kv_fromBus = get_bus_vn_kV(net = myNet, busName = fbus)
    vn_kv_toBus = get_bus_vn_kV(net = myNet, busName = tbus)

    hasBusInNet(net = myNet, busName = fbus) || (@warn("bus $(fbus) not found, branch ignored."); continue)
    hasBusInNet(net = myNet, busName = tbus) || (@warn("bus $(tbus) not found, branch ignored."); continue)

    r_pu = float(row[branchDict["r"]])
    x_pu = float(row[branchDict["x"]])
    b_pu = float(row[branchDict["b"]])
    ratedS_raw = float(row[branchDict["rateA"]])
    ratedS = ratedS_raw > 0.0 ? ratedS_raw : Inf

    #ratedS = float(row[branchDict["rateA"]])
    ratio = float(row[branchDict["ratio"]])
    ratio == 0.0 && (ratio = 1.0) # avoid division by zero
    angle = float(row[branchDict["angle"]])
    status = Int64(row[branchDict["status"]])
    isLine = false
    if (ratio == 1.0 && angle == 0.0 || ratio == 0.0) && vn_kv_fromBus == vn_kv_toBus
      isLine = true
    end

    if isLine
      addPIModelACLine!(net = myNet, fromBus = fbus, toBus = tbus, r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, status = status, ratedS = ratedS)
    else # transformer
      addPIModelTrafo!(net = myNet, fromBus = fbus, toBus = tbus, r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, status = status, ratedS = ratedS, ratio = ratio, shift_deg = angle)
    end
  end

  # Generators:   
  col = size(genData, 2)
  for row in eachrow(genData[:, 1:col])
    status = Int64(row[genDict["status"]])
    _bus = Int(row[genDict["bus"]])
    bus = string(_bus)
    if status < 1
      pInfo("generator $(_bus) not in service, ignored")
      continue
    end

    pGen = float(row[genDict["Pg"]]) < 0.0 ? (pInfo("pGen at bus $(_bus) p < 0"); float(row[genDict["Pg"]])) : float(row[genDict["Pg"]])
    qGen = float(row[genDict["Qg"]]) < 0.0 ? (@debug("qGen at bus $(_bus) q < 0"); float(row[genDict["Qg"]])) : float(row[genDict["Qg"]])

    if abs(pGen) < 1e-6 && abs(qGen) < 1e-6
      busType = getBusType(net = myNet, busName = bus)
      if busType == PQ
        pInfo("generator $(_bus) has no power output, ignored")
        continue
      end
    end

    qMax = float(row[genDict["Qmax"]])
    qMin = float(row[genDict["Qmin"]])
    pMax = float(row[genDict["Pmax"]])
    pMin = float(row[genDict["Pmin"]])
    vm_pu = float(row[genDict["Vg"]])
    mBase = float(row[genDict["mBase"]])

    referencePri = slackIdx == _bus ? bus : nothing

    if mBase != baseMVA
      @warn "generator $(_bus) has different baseMVA than network baseMVA"
    end

    addProsumer!(net = myNet, busName = bus, type = "GENERATOR", vm_pu = vm_pu, referencePri = referencePri, p = pGen, q = qGen, pMax = pMax, pMin = pMin, qMax = qMax, qMin = qMin)
  end# read Generators

  erg, msg = validate!(net = myNet)
  if erg == false
    @error "network is invalid: $msg"
  end
  return myNet
end

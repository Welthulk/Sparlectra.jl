# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 01.10.2023
# file: src/createnet_powermat.jl
# helper
#! format: off


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


#! format: on

"""
    createNetFromMatPowerCase(; mpc, log=false, flatstart=false) -> Net

Builds a Sparlectra `Net` from a MATPOWER-like container `mpc`.

`mpc` can be either:
- a `NamedTuple` with fields `name, baseMVA, bus, gen, branch` (optionally `gencost, bus_name`)
- or a struct with the same field names (e.g. `MatpowerCase`)

All matrices are expected in MATPOWER v2 column conventions.
"""
function createNetFromMatPowerCase(; mpc, log::Bool=false, flatstart::Bool=false)::Net
  # Small logger helper
  pInfo(msg::String) = (log ? (@info msg) : nothing)

  # --- Extract fields from NamedTuple / struct uniformly ---
  name    = hasproperty(mpc, :name)    ? getproperty(mpc, :name)    : "mpc"
  baseMVA = hasproperty(mpc, :baseMVA) ? Float64(getproperty(mpc, :baseMVA)) : error("mpc.baseMVA missing")
  busData = hasproperty(mpc, :bus)     ? Matrix{Float64}(getproperty(mpc, :bus)) : error("mpc.bus missing")
  genData = hasproperty(mpc, :gen)     ? Matrix{Float64}(getproperty(mpc, :gen)) : error("mpc.gen missing")
  brData  = hasproperty(mpc, :branch)  ? Matrix{Float64}(getproperty(mpc, :branch)) : error("mpc.branch missing")

  # --- Legacy-compatible dicts (same as your old importer) ---
  busDict, genDict, branchDict = _createDict()

  if log
    @info "Creating new Net: $(name) with baseMVA=$(baseMVA), flatstart=$(flatstart)"
  end
  
  myNet = Net(name = String(name), baseMVA = baseMVA, flatstart = flatstart)

  # --- Find slack bus index from BUS_TYPE==3 (MATPOWER) ---
  slackIdx = 0
  for row in eachrow(busData)
    btype = Int(row[busDict["type"]])
    if btype == 3
      slackIdx = Int(row[busDict["bus"]])
      break
    end
  end

  # --- Buses + Loads + Shunts (same semantics as your existing importer) ---
  mFak = 10.0

  for row in eachrow(busData)
    btype = Int(row[busDict["type"]])
    kIdx  = Int(row[busDict["bus"]])     # original bus number
    busName = string(kIdx)

    raw_vn = Float64(row[busDict["baseKV"]])
    vn_kv  = raw_vn <= 0.0 ? 1.0 : raw_vn

    va_deg = Float64(row[busDict["Va"]])
    vm_pu  = Float64(row[busDict["Vm"]])
    vm_pu  = vm_pu <= 0.0 ? 1.0 : vm_pu

    pLoad  = Float64(row[busDict["Pd"]])
    qLoad  = Float64(row[busDict["Qd"]])

    pShunt = Float64(row[busDict["Gs"]])
    qShunt = Float64(row[busDict["Bs"]])

    zone = Int(row[busDict["zone"]])
    area = Int(row[busDict["area"]])
    vmin = Float64(row[busDict["Vmin"]])
    vmax = Float64(row[busDict["Vmax"]])


    addBus!(
      net = myNet,
      busName = busName,
      busType = toString(toNodeType(btype)),
      vn_kV = vn_kv,
      vm_pu = vm_pu,
      va_deg = va_deg,
      vmin_pu = vmin,
      vmax_pu = vmax,
      isAux = false,
      oBusIdx = kIdx,
      zone = zone,
      area = area,
    )
    # fix 04.02.2026: MATPOWER -> p.u. admittance components:
    pShunt_pu = pShunt / baseMVA
    qShunt_pu = qShunt / baseMVA
    if pShunt != 0.0 || qShunt != 0.0
      addShunt!(net=myNet, busName=busName, pShunt=pShunt_pu, qShunt=qShunt_pu)
    end

    if pLoad != 0.0 || qLoad != 0.0
      qMax = min(abs(mFak * qLoad), baseMVA)
      qMin = -qMax
      pMax = min(abs(mFak * pLoad), baseMVA)
      pMin = -pMax
      addProsumer!(
        net = myNet,
        busName = busName,
        type = "ENERGYCONSUMER",
        p = pLoad,
        q = qLoad,
        pMin = pMin,
        pMax = pMax,
        qMin = qMin,
        qMax = qMax,
      )
    end
  end

  # --- Branches (line vs transformer decision same as old importer) ---
  for row in eachrow(brData)
    fbus = string(Int(row[branchDict["fbus"]]))
    tbus = string(Int(row[branchDict["tbus"]]))

    hasBusInNet(net = myNet, busName = fbus) || (@warn "bus $(fbus) not found, branch ignored."; continue)
    hasBusInNet(net = myNet, busName = tbus) || (@warn "bus $(tbus) not found, branch ignored."; continue)

    r_pu = Float64(row[branchDict["r"]])
    x_pu = Float64(row[branchDict["x"]])
    b_pu = Float64(row[branchDict["b"]])

    ratedS_raw = Float64(row[branchDict["rateA"]])
    ratedS = ratedS_raw > 0.0 ? ratedS_raw : Inf

    ratio = Float64(row[branchDict["ratio"]])
    ratio = (ratio == 0.0) ? 1.0 : ratio

    angle  = Float64(row[branchDict["angle"]])
    status = Int(row[branchDict["status"]])

    vn_from = get_bus_vn_kV(net = myNet, busName = fbus)
    vn_to   = get_bus_vn_kV(net = myNet, busName = tbus)

    isLine = ((ratio == 1.0 && angle == 0.0) && (vn_from == vn_to))

    if isLine
      addPIModelACLine!(
        net = myNet,
        fromBus = fbus,
        toBus = tbus,
        r_pu = r_pu,
        x_pu = x_pu,
        b_pu = b_pu,
        status = status,
        ratedS = ratedS,
      )
    else
      addPIModelTrafo!(
        net = myNet,
        fromBus = fbus,
        toBus = tbus,
        r_pu = r_pu,
        x_pu = x_pu,
        b_pu = b_pu,
        status = status,
        ratedS = ratedS,
        ratio = ratio,
        shift_deg = angle,
      )
    end
  end

  # --- Generators ---
  for row in eachrow(genData)
    status = Int(row[genDict["status"]])
    status < 1 && continue

    bus = string(Int(row[genDict["bus"]]))

    pGen = Float64(row[genDict["Pg"]])
    qGen = Float64(row[genDict["Qg"]])

    qMax = Float64(row[genDict["Qmax"]])
    qMin = Float64(row[genDict["Qmin"]])
    pMax = Float64(row[genDict["Pmax"]])
    pMin = Float64(row[genDict["Pmin"]])
    vm_pu = Float64(row[genDict["Vg"]])
    mBase = Float64(row[genDict["mBase"]])

    referencePri = (slackIdx == Int(row[genDict["bus"]])) ? bus : nothing
    (mBase != baseMVA) && @warn "generator $(bus) has different baseMVA than network baseMVA"

    addProsumer!(
      net = myNet,
      busName = bus,
      type = "GENERATOR",
      vm_pu = vm_pu,
      referencePri = referencePri,
      p = pGen,
      q = qGen,
      pMax = pMax,
      pMin = pMin,
      qMax = qMax,
      qMin = qMin,
    )
  end

  ok, msg = validate!(net = myNet)
  ok || @error "network is invalid: $msg"

  return myNet
end

function createNetFromMatPowerFile(; filename::String, log::Bool=false, flatstart::Bool=false)::Net
  mpc = MatpowerIO.read_case(filename; legacy_compat=true)
  return createNetFromMatPowerCase(; mpc=mpc, log=log, flatstart=flatstart)
end

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
function createNetFromMatPowerCase(; mpc, log::Bool=false, flatstart::Bool=false, cooldown::Int = 0, q_hyst_pu::Float64 = 0.0, enable_pq_gen_controllers::Bool = true)::Net
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
  BUS_I = busDict["bus"]; BUS_TYPE = busDict["type"]; PD = busDict["Pd"]; QD = busDict["Qd"]; GS = busDict["Gs"]; BS = busDict["Bs"]
  BUS_AREA = busDict["area"]; VM = busDict["Vm"]; VA = busDict["Va"]; BASE_KV = busDict["baseKV"]; BUS_ZONE = busDict["zone"]; VMAX = busDict["Vmax"]; VMIN = busDict["Vmin"]
  GEN_BUS = genDict["bus"]; PG = genDict["Pg"]; QG = genDict["Qg"]; QMAX = genDict["Qmax"]; QMIN = genDict["Qmin"]; VG = genDict["Vg"]
  MBASE = genDict["mBase"]; GEN_STATUS = genDict["status"]; PMAX = genDict["Pmax"]; PMIN = genDict["Pmin"]
  F_BUS = branchDict["fbus"]; T_BUS = branchDict["tbus"]; BR_R = branchDict["r"]; BR_X = branchDict["x"]; BR_B = branchDict["b"]
  RATE_A = branchDict["rateA"]; TAP = branchDict["ratio"]; SHIFT = branchDict["angle"]; BR_STATUS = branchDict["status"]

  mp_bus_type = Dict{Int,Int}()
  sizehint!(mp_bus_type, size(busData, 1))
  for row in eachrow(busData)
    mp_bus_type[Int(row[BUS_I])] = Int(row[BUS_TYPE])
  end

  if log
    @info "Creating new Net: $(name) with baseMVA=$(baseMVA), flatstart=$(flatstart)"
  end
  
  myNet = Net(name = String(name), baseMVA = baseMVA, flatstart = flatstart, cooldown_iters = cooldown, q_hyst_pu = q_hyst_pu)
  nbus = size(busData, 1)
  nbranch = size(brData, 1)
  ngen = size(genData, 1)
  sizehint!(myNet.nodeVec, nbus)
  sizehint!(myNet.busDict, nbus)
  sizehint!(myNet.busOrigIdxDict, nbus)
  sizehint!(myNet.branchVec, nbranch)
  sizehint!(myNet.linesAC, nbranch)
  sizehint!(myNet.trafos, nbranch)
  sizehint!(myNet.prosumpsVec, nbus + ngen)
  sizehint!(myNet.shuntVec, nbus)
  sizehint!(myNet.shuntDict, nbus)
  bus_idx_by_orig = Dict{Int,Int}()
  sizehint!(bus_idx_by_orig, nbus)

  # --- Find slack bus index from BUS_TYPE==3 (MATPOWER) ---
  slackIdx = 0
  for row in eachrow(busData)
    btype = Int(row[BUS_TYPE])
    if btype == 3
      slackIdx = Int(row[BUS_I])
      break
    end
  end

  # --- Buses + Loads + Shunts (same semantics as your existing importer) ---
  mFak = 10.0

  for row in eachrow(busData)
    btype = Int(row[BUS_TYPE])
    kIdx  = Int(row[BUS_I])     # original bus number
    busName = string(kIdx)

    raw_vn = Float64(row[BASE_KV])
    vn_kv  = raw_vn <= 0.0 ? 1.0 : raw_vn

    va_deg = Float64(row[VA])
    vm_pu  = Float64(row[VM])
    vm_pu  = vm_pu <= 0.0 ? 1.0 : vm_pu

    pLoad  = Float64(row[PD])
    qLoad  = Float64(row[QD])

    pShunt = Float64(row[GS])
    qShunt = Float64(row[BS])

    zone = Int(row[BUS_ZONE])
    area = Int(row[BUS_AREA])
    vmin = Float64(row[VMIN])
    vmax = Float64(row[VMAX])


    addBus!(
      net = myNet,
      busName = busName,
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
    #=
    pShunt_pu = pShunt / baseMVA
    qShunt_pu = qShunt / baseMVA
    if pShunt != 0.0 || qShunt != 0.0
      addShunt!(net=myNet, busName=busName, pShunt=pShunt_pu, qShunt=qShunt_pu)
    end
    =#
    # new 04.02.2026: directly add shunt as bus shunt
    pShunt = Float64(row[GS])
    qShunt = Float64(row[BS])

    bus_idx_by_orig[kIdx] = length(myNet.nodeVec)

    if pShunt != 0.0 || qShunt != 0.0
      addShuntMatpower!(net=myNet, busName=busName, Gs=pShunt, Bs=qShunt)
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
        defer_bus_type_refresh = true,
      )
    end
  end

  # --- Branches (line vs transformer decision same as old importer) ---
  for row in eachrow(brData)
    fbus_orig = Int(row[F_BUS])
    tbus_orig = Int(row[T_BUS])
    from_idx = get(bus_idx_by_orig, fbus_orig, 0)
    to_idx = get(bus_idx_by_orig, tbus_orig, 0)

    from_idx != 0 || (@warn "bus $(fbus_orig) not found, branch ignored."; continue)
    to_idx != 0 || (@warn "bus $(tbus_orig) not found, branch ignored."; continue)

    r_pu = Float64(row[BR_R])
    x_pu = Float64(row[BR_X])
    b_pu = Float64(row[BR_B])

    ratedS_raw = Float64(row[RATE_A])
    ratedS = ratedS_raw > 0.0 ? ratedS_raw : Inf

    ratio = Float64(row[TAP])
    ratio = (ratio == 0.0) ? 1.0 : ratio

    angle  = Float64(row[SHIFT])
    status = Int(row[BR_STATUS])

    vn_from = getNodeVn(myNet.nodeVec[from_idx])
    vn_to   = getNodeVn(myNet.nodeVec[to_idx])

    isLine = ((ratio == 1.0 && angle == 0.0) && (vn_from == vn_to))

    if isLine
      _addPIModelACLine_by_idx!(
        net = myNet,
        from = from_idx,
        to = to_idx,
        r_pu = r_pu,
        x_pu = x_pu,
        b_pu = b_pu,
        status = status,
        ratedS = ratedS,
      )
    else
      _addPIModelTrafo_by_idx!(
        net = myNet,
        from = from_idx,
        to = to_idx,
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
  pq_gen_controller_count = 0
  for row in eachrow(genData)
    status = Int(row[GEN_STATUS])
    status < 1 && continue

    bus = string(Int(row[GEN_BUS]))

    pGen = Float64(row[PG])
    qGen = Float64(row[QG])

    qMax = Float64(row[QMAX])
    qMin = Float64(row[QMIN])
    pMax = Float64(row[PMAX])
    pMin = Float64(row[PMIN])
    vm_pu = Float64(row[VG])
    mBase = Float64(row[MBASE])
    btype = get(mp_bus_type, Int(row[GEN_BUS]), 1)

    referencePri = (slackIdx == Int(row[GEN_BUS])) ? bus : nothing
    (mBase != baseMVA) && @debug "generator $(bus) has different mBase than network baseMVA (allowed in MATPOWER)" bus mBase baseMVA

    pu_controller = nothing
    qu_controller = nothing
    if enable_pq_gen_controllers && btype == 1
      p_pu = pGen / baseMVA
      q_pu = qGen / baseMVA
      pu_controller = PUController(make_characteristic([(0.0, p_pu), (2.0, p_pu)]); pmin_MW = pMin, pmax_MW = pMax, sbase_MVA = baseMVA)
      qu_controller = QUController(make_characteristic([(0.0, q_pu), (2.0, q_pu)]); qmin_MVAr = qMin, qmax_MVAr = qMax, sbase_MVA = baseMVA)
      pq_gen_controller_count += 1
    end

    addProsumer!(
      net = myNet,
      busName = bus,
      type = "GENERATOR",
      vm_pu = (btype == 2 || btype == 3) ? vm_pu : nothing,
      referencePri = referencePri,
      p = pGen,
      q = qGen,
      pMax = pMax,
      pMin = pMin,
      qMax = qMax,
      qMin = qMin,
      qu_controller = qu_controller,
      pu_controller = pu_controller,
      isRegulated = (btype == 2),
      defer_bus_type_refresh = true,
    )
  end

  refreshBusTypesFromProsumers!(myNet)
  _buildQLimits!(myNet)

  if enable_pq_gen_controllers && pq_gen_controller_count > 0
    @info "MATPOWER import: PQ generator limits mapped to constant P(U)/Q(U) controllers" count = pq_gen_controller_count
  end

  ok, msg = validate!(net = myNet)
  ok || @error "network is invalid: $msg"

  return myNet
end

function createNetFromMatPowerFile(; filename::String,
    log::Bool=false,
    flatstart::Bool=false,
    cooldown::Int = 0,
    q_hyst_pu::Float64 = 0.0,
    enable_pq_gen_controllers::Bool = true,
    verbose::Int = 0)::Net

  mpc = MatpowerIO.read_case(filename; legacy_compat=true)

  # Build the network first
  net = createNetFromMatPowerCase(; mpc=mpc, log=log, flatstart=flatstart,
                                  cooldown=cooldown, q_hyst_pu=q_hyst_pu, enable_pq_gen_controllers=enable_pq_gen_controllers)

  # Always apply MATPOWER isolated flags (BUS_TYPE==4) onto net
  MatpowerIO.apply_mp_isolated_buses!(net, mpc; verbose=verbose)

  # Only if not flatstart: take VM/VA from mpc.bus as initial guess
  MatpowerIO.apply_mp_bus_vmva_init!(net, mpc; flatstart=flatstart, verbose=verbose)

  return net
end

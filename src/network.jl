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
# Date: 01.04.2024
# file: src/network.jl
# purpose: the central Net container and its construction API: addBus!,
#          addACLine!, add2WTrafo!, addProsumer!, addShunt!, addLink! and
#          friends, plus validation, locking, and bus/branch lookup helpers
"""
    Net

Represents an electrical network.

# Fields
- `name::String`: Name of the network.
- `baseMVA::Float64`: Base MVA of the network.
- `slackVec::Vector{Int}`: Vector containing indices of slack buses.
- `vmin_pu::Float64`: Minimum voltage limit in per unit.
- `vmax_pu::Float64`: Maximum voltage limit in per unit.
- `nodeVec::Vector{Node}`: Vector containing nodes in the network.
- `linesAC::Vector{ACLineSegment}`: Vector containing AC line segments.
- `trafos::Vector{PowerTransformer}`: Vector containing power transformers.
- `branchVec::Vector{Branch}`: Vector containing branches in the network.
- `prosumpsVec::Vector{ProSumer}`: Vector containing prosumers in the network.
- `shuntVec::Vector{Shunt}`: Vector containing shunts in the network.
- `busDict::Dict{String,Int}`: Dictionary mapping bus names to indices.
- `busOrigIdxDict::Dict{Int,Int}`: Dictionary mapping current bus indices to original indices.
- `busOriginalNameDict::Dict{Int,String}`: Optional original source-file bus names keyed by current bus index.
- `branchDict::Dict{Tuple{Int, Int},Int}`: Dictionary mapping branch tuples to indices.
- `totalLosses::Vector{Tuple{Float64,Float64}}`: Vector containing tuples of total power losses.
- `_locked::Bool`: Boolean indicating if the network is locked.
- `measurements::Vector`: State-estimation measurements stored on the network.
- `matpower_branch_metadata::Dict{Int,NamedTuple}`: Optional MATPOWER branch metadata keyed by imported branch index.
- `for001Contingencies::Vector{String}`: Optional FOR001 contingency branch names imported from MATPOWER metadata.
- `matpowerDclineMetadata::Vector{NamedTuple}`: Optional imported MATPOWER DC-line terminal-injection metadata.
- `cgmes_ids::Dict{String,String}`: CGMES mRIDs keyed by structural identity (e.g. `"TN|<bus>"`, `"ACL|<busA>|<busB>|<k>"`); filled by the CGMES importer and reused by the CGMES exporter so exported ids stay stable.

# Constructors
- `Net(name::String, baseMVA::Float64, vmin_pu::Float64 = 0.9, vmax_pu::Float64 = 1.1)`: Creates a new `Net` object with the given name, base MVA, and optional voltage limits.

# Functions
- `addBus!(; net::Net, ...)`: Adds a bus to the network.
- `addShunt!(; net::Net, ...)`: Adds a shunt to the network.
- `addBranch!(; net::Net, ...)`: Adds a branch to the network.
- `addPIModelACLine!(; net::Net, ...)`: Adds a PI model AC line to the network.
- `addACLine!(; net::Net, ...)`: Adds an AC line segment to the network.
- `addPIModelTrafo!(; net::Net, ...)`: Adds a transformer with PI model to the network.
- `add2WTrafo!(; net::Net, ...)`: Adds a two-winding transformer to the network.
- `addProsumer!(; net::Net, ...)`: Adds a prosumer to the network.
- `lockNet!(; net::Net, locked::Bool)`: Locks or unlocks the network.
- `validate!(; net::Net)`: Validates the network.
- `get_bus_vn_kV(; net::Net, busName::String)`: Gets the voltage level of a bus by name.
- `get_vn_kV(; net::Net, busIdx::Int)`: Gets the voltage level of a bus by index.
- `getBusType(; net::Net, busName::String)`: Gets the type of a bus by name.
- `updateBranchParameters!(; net::Net, fromBus::String, toBus::String, branch::BranchModel)`: Updates the parameters of a branch in the network.
- `setNetBranchStatus!(; net::Net, branchNr::Int, status::Int)`: Sets the status of a branch.
- `setTotalLosses!(; net::Net, pLosses::Float64, qLosses::Float64)`: Sets the total losses of the network.
- `getTotalLosses(; net::Net)`: Gets the total losses of the network.
- `getNetOrigBusIdx(; net::Net, busName::String)`: Gets the original index of a bus in the network.
- `getNetBusIdx(; net::Net, busName::String)`: Gets the index of a bus in the network.
- `hasBusInNet(; net::Net, busName::String)`: Checks if a bus exists in the network.
- `addBusGenPower!(; net::Net, busName::String, pGen::Float64, qGen::Float64)`: Adds generator power to a bus.
- `addBusLoadPower!(; net::Net, busName::String, pLoad::Float64, qLoad::Float64)`: Adds load power to a bus.
- `getNetBranch(; net::Net, fromBus::String, toBus::String)`: Retrieves the branch between two specified buses in the network.
"""

mutable struct Net
  name::String
  baseMVA::Float64
  slackVec::Vector{Int}
  vmin_pu::Float64
  vmax_pu::Float64
  nodeVec::Vector{Node}
  linesAC::Vector{ACLineSegment}
  trafos::Vector{PowerTransformer}
  branchVec::Vector{Branch}
  linkVec::Vector{BusLink}
  prosumpsVec::Vector{ProSumer}
  shuntVec::Vector{Shunt}
  busDict::Dict{String,Int}
  busOrigIdxDict::Dict{Int,Int}
  busOriginalNameDict::Dict{Int,String}
  totalLosses::Vector{Tuple{Float64,Float64}}
  totalBusPower::Vector{Tuple{Float64,Float64}}
  _locked::Bool
  flatstart::Bool                     # flatstart flag for power flow initialization
  shuntDict::Dict{Int,Int}
  isoNodes::Vector{Int}
  qLimitLog::Vector{Any}
  cooldown_iters::Int
  q_hyst_pu::Float64
  qmin_pu::Vector{Float64}            # pro Bus Qmin (p.u.)
  qmax_pu::Vector{Float64}            # pro Bus Qmax (p.u.)
  qLimitInitialPVRows::Vector{Any}    # pre-solve PV Q-limit snapshot
  qLimitEvents::Dict{Int,Symbol}      # BusIdx -> :min | :max (PV→PQ Change)  
  measurements::Vector
  bus_shunt_model::Symbol
  control_result::Union{Nothing,ControlRunResult}
  matpower_branch_metadata::Dict{Int,NamedTuple}
  for001Contingencies::Vector{String}
  matpowerDclineMetadata::Vector{NamedTuple}
  # Outer-loop machine voltage controllers (remote Q-regulation). Typed with
  # the abstract controller type because the concrete MachineVoltageControl is
  # defined after this file; filled via addMachineVoltageControl!.
  machineControls::Vector{AbstractOuterController}
  # CGMES mRIDs keyed by structural identity ("TN|<bus>", "ACL|<busA>|<busB>|<k>", ...).
  # Filled by the CGMES importer, read by the CGMES exporter; empty for nets
  # from other sources (the exporter then mints deterministic uuid5 ids).
  cgmes_ids::Dict{String,String}
  # Native short-circuit source records (issue #299): filled by
  # addExternalGrid!, read by runShortCircuit!(net). Carrying data here must
  # never change power-flow results by itself (participationFactor precedent).
  sc_sources::NativeShortCircuitData
  # Persistent HVDC link records (r0.9.9): one entry per link regardless of
  # source (MATPOWER dcline, CGMES DC-topology pair, API) or mode (Stage-0
  # fixed injections or attached pair controller). Read by the result layer;
  # live values always come from the prosumers/controller.
  hvdcLinks::Vector{HvdcLink}
  # Solver status of the most recent rectangular/DC solve on this Net
  # (thread-safety Phase 1: replaces the former global weak-ref registries,
  # which raced under concurrent solves). Contract: one Net is solved by at
  # most one task at a time; deepcopy(net) deliberately carries the status
  # with the copy (the island path relies on per-copy status), and any
  # future saveNet-style serializer must SKIP both fields: they hold solver
  # internals including closures (condition-estimate thunk), not model data.
  _rectangular_pf_status::Any
  _dc_pf_status::Any

  #! format: off
  function Net(; name::String, baseMVA::Float64, vmin_pu::Float64 = 0.9, vmax_pu::Float64 = 1.1, cooldown_iters::Int = 0, q_hyst_pu::Float64 = 0.0, flatstart::Bool = false, bus_shunt_model = :admittance)    
    shunt_model = normalize_bus_shunt_model(bus_shunt_model)
    
    new(name, # name
        baseMVA, # baseMVA
        [], # slackVec
        vmin_pu, # vmin_pu
        vmax_pu, # vmax_pu
        [], # nodeVec
        [], # linesAC
        [], # trafos
        [], # branchVec
        [], # linkVec
        [], # prosumpsVec
        [], # shuntVec
        Dict{String,Int}(), # busDict
        Dict{Int,Int}(), # busOrigIdxDict
        Dict{Int,String}(), # busOriginalNameDict
        [], # totalLosses
        [], # totalBusPower
        false, # _locked
        flatstart,                     # flatstart
        Dict{Int,Int}(),  # shuntDict
        [],                                    # isoNodes
        Any[],                                 # qLimitLog                     
        cooldown_iters,                        # cooldown_iters
        q_hyst_pu,
        [],                                    # qmin_pu
        [],                                    # qmax_pu
        Any[],                                 # qLimitInitialPVRows
        Dict{Int,Symbol}(),
        [],
        shunt_model,
        nothing,
        Dict{Int,NamedTuple}(),
        String[],
        NamedTuple[],
        AbstractOuterController[],                 # machineControls
        Dict{String,String}(),                     # cgmes_ids
        NativeShortCircuitData(),                  # sc_sources
        HvdcLink[],                                # hvdcLinks
        nothing,                                   # _rectangular_pf_status
        nothing)                                   # _dc_pf_status
  end
  #! format: on
  function Base.show(io::IO, net::Net)
    println(io, "Net: ", net.name)
    println(io, "Base MVA: ", net.baseMVA)
    println(io, "Nodes: ", length(net.nodeVec), ", Lines: ", length(net.linesAC), ", Transformers: ", length(net.trafos), ", Branches: ", length(net.branchVec), ", Links: ", length(net.linkVec), ", Prosumers: ", length(net.prosumpsVec), ", Shunts: ", length(net.shuntVec))
    println(io, "Slack buses: ", net.slackVec, ", flatstart: ", net.flatstart, ", locked: ", net._locked)
    println(io, "Vmin / Vmax: ", net.vmin_pu, " / ", net.vmax_pu)
    println(io, "cooldown_iters: ", net.cooldown_iters, ", q_hyst_pu: ", net.q_hyst_pu)
    println(io, "Measurements: ", length(net.measurements))
    println(io, "Tap controllers: ", sum(length, (t.side1.controls for t in net.trafos); init = 0) + sum(length, (t.side2.controls for t in net.trafos); init = 0) + sum((isnothing(t.side3) ? 0 : length(t.side3.controls) for t in net.trafos); init = 0))
    isempty(net.machineControls) || println(io, "Machine controllers: ", length(net.machineControls))
  end
end

const BUS_SHUNT_ADMITTANCE = :admittance
const BUS_SHUNT_VOLTAGE_DEPENDENT_INJECTION = :voltage_dependent_injection
const BUS_SHUNT_MODEL_VALUES = (BUS_SHUNT_ADMITTANCE, BUS_SHUNT_VOLTAGE_DEPENDENT_INJECTION)

"""
    normalize_bus_shunt_model(value) -> Symbol

Validate and normalize the bus-shunt modeling option. Supported values are
`"admittance"` and `"voltage_dependent_injection"`.
"""
function normalize_bus_shunt_model(value)::Symbol
  s = lowercase(String(value))
  s = replace(s, "-" => "_")
  model = Symbol(s)
  model in BUS_SHUNT_MODEL_VALUES || error("Invalid bus_shunt_model=$(value). Supported values are \"admittance\" and \"voltage_dependent_injection\".")
  return model
end

_shunt_component_model(model::Symbol)::Symbol = model == BUS_SHUNT_ADMITTANCE ? :Y : :VoltageDependentInjection

"""
    bus_shunt_totals_pu(net) -> NamedTuple

Return the count and total per-unit conductance/susceptance of in-service bus
shunts. The totals are based on the stored shunt admittances and are useful for
compact import/configuration logging.
"""
function bus_shunt_totals_pu(net::Net)
  g = 0.0
  b = 0.0
  count = 0
  for sh in net.shuntVec
    sh.status == 0 && continue
    count += 1
    g += real(sh.y_pu_shunt)
    b += imag(sh.y_pu_shunt)
  end
  return (count = count, total_g_pu = g, total_b_pu = b)
end

"""
    log_bus_shunt_model(net) -> Nothing

Emit a compact log message with the selected bus-shunt modeling mode, in-service
bus-shunt count, and aggregate per-unit conductance/susceptance.
"""
function log_bus_shunt_model(net::Net)
  totals = bus_shunt_totals_pu(net)
  @info "Bus shunt model" bus_shunt_model = String(net.bus_shunt_model) bus_shunt_count = totals.count bus_shunt_total_g_pu = totals.total_g_pu bus_shunt_total_b_pu = totals.total_b_pu
end

function showNet(io::IO, net::Net; verbose::Bool = false)
  if !verbose
    show(io, net)
    return
  end

  println(io, "==================== Net ====================")
  println(io, "Name:          ", net.name)
  println(io, "Base MVA:      ", net.baseMVA)
  println(io, "Slack buses:   ", net.slackVec)
  println(io, "Vmin / Vmax:   ", net.vmin_pu, " / ", net.vmax_pu)
  println(io, "flatstart:     ", net.flatstart)
  println(io, "_locked:       ", net._locked)

  println(io, "\n--- Dictionaries ---")
  println(io, "busDict:        ", net.busDict)
  println(io, "busOrigIdxDict: ", net.busOrigIdxDict)
  println(io, "shuntDict:      ", net.shuntDict)

  println(io, "\n--- Topology / State ---")
  println(io, "isoNodes:       ", net.isoNodes)

  println(io, "\n--- Q-limit handling ---")
  println(io, "cooldown_iters: ", net.cooldown_iters)
  println(io, "q_hyst_pu:      ", net.q_hyst_pu)
  println(io, "qmin_pu:        ", net.qmin_pu)
  println(io, "qmax_pu:        ", net.qmax_pu)
  println(io, "qLimitEvents:   ", net.qLimitEvents)
  println(io, "measurements:   ", length(net.measurements))

  println(io, "\n--- Aggregates ---")
  println(io, "totalLosses:    ", net.totalLosses)
  println(io, "totalBusPower:  ", net.totalBusPower)

  println(io, "\n--- Q-limit log (", length(net.qLimitLog), ") ---")
  for (i, e) in enumerate(net.qLimitLog)
    println(io, "  [", i, "] ", e)
  end

  println(io, "\n==================== Nodes (", length(net.nodeVec), ") ====================")
  for (i, n) in enumerate(net.nodeVec)
    println(io, "\n[node ", i, "]")
    show(io, n)
  end

  println(io, "\n==================== AC Lines (", length(net.linesAC), ") ====================")
  for (i, l) in enumerate(net.linesAC)
    println(io, "\n[line ", i, "]")
    show(io, l)
  end

  println(io, "\n==================== Transformers (", length(net.trafos), ") ====================")
  for (i, t) in enumerate(net.trafos)
    println(io, "\n[trafo ", i, "]")
    show(io, t)
  end

  println(io, "\n==================== Branches (", length(net.branchVec), ") ====================")
  for (i, b) in enumerate(net.branchVec)
    println(io, "\n[branch ", i, "]")
    show(io, b)
  end

  println(io, "\n==================== Prosumers (", length(net.prosumpsVec), ") ====================")
  for (i, p) in enumerate(net.prosumpsVec)
    println(io, "\n[prosumer ", i, "]")
    show(io, p)
  end

  println(io, "\n==================== Shunts (", length(net.shuntVec), ") ====================")
  for (i, s) in enumerate(net.shuntVec)
    println(io, "\n[shunt ", i, "]")
    show(io, s)
  end

  println(io, "\n==================== Links (", length(net.linkVec), ") ====================")
  for (i, l) in enumerate(net.linkVec)
    println(io, "\n[link ", i, "]")
    show(io, l)
  end
  if length(net.measurements) == 0
    println(io, "\n==================== No measurements ====================")
    return
  else
    Println(io, "\n==================== Mesurements ====================")
    for (i, m) in enumerate(net.measurements)
      println(io, "\n[measurement ", i, "]")
      show(io, m)
    end
  end
  println(io, "\n==================== END Net ====================")
end
# convenience overloads
showNet(net::Net; verbose::Bool = false) = showNet(stdout, net; verbose = verbose)

# --- helpers (lokal) ---
@inline function _push_unique!(v::Vector{Int}, x::Int)
  (findfirst(==(x), v) === nothing) && push!(v, x)
  return v
end
@inline function _delete_item!(v::Vector{Int}, x::Int)
  idx = findfirst(==(x), v)
  (idx !== nothing) && deleteat!(v, idx)
  return v
end

const _setbus_bus_type_warned = Ref(false)
const _legacy_bus_type_trace_keys = ("1", "true", "yes", "on")

@inline function _trace_legacy_bus_type_warnings()::Bool
  return lowercase(get(ENV, "SPARLECTRA_TRACE_LEGACY_BUSTYPE", "0")) in _legacy_bus_type_trace_keys
end

function _legacy_bus_type_caller_hint()::String
  bt = stacktrace()
  for fr in bt
    f = String(fr.file)
    if !occursin("src/network.jl", f)
      return "$(basename(f)):$(fr.line)"
    end
  end
  return "unknown"
end

"""
Sets the type of a bus in the network.

# Parameters:
- `net::Net`: Network object.
- `bus::Int`: Index of the bus.
- `busType::String`: Type of the bus (e.g., "Slack", "PQ", "PV")
"""
function setBusType!(net::Net, bus::Int, busType::String)
  if !_setbus_bus_type_warned[]
    msg = "setBusType! only updates the static node label. Power-flow bus typing is derived from attached prosumers."
    if _trace_legacy_bus_type_warnings()
      msg *= " caller=" * _legacy_bus_type_caller_hint()
    end
    @warn msg
    _setbus_bus_type_warned[] = true
  end
  @assert 1 <= bus <= length(net.nodeVec) "The $bus index is invalid. It must be between 1 and $(length(net.nodeVec))."
  node  = net.nodeVec[bus]
  oldTy = getfield(node, :_nodeType)

  setNodeType!(node, busType)

  newTy = getfield(node, :_nodeType)
  if newTy == Slack && oldTy != Slack
    _push_unique!(net.slackVec, bus)
  elseif oldTy == Slack && newTy != Slack
    _delete_item!(net.slackVec, bus)
  end
  return nothing
end

function setBusType!(net::Net, busName::String, busType::String)
  bus = geNetBusIdx(net = net, busName = busName)
  return setBusType!(net, bus, busType)
end

"""
hasBusInNet: Checks if a bus exists in the network.

Parameters:
- `net::Net`: Network object.
- `busName::String`: Name of the bus to check.

Returns:
- `Bool`: True if the bus exists in the network, otherwise false.
"""
function hasBusInNet(; net::Net, busName::String)::Bool
  return haskey(net.busDict, busName)
end

"""
geNetBusIdx: Gets the index of a bus in the network.

Parameters:
- `net::Net`: Network object.
- `busName::String`: Name of the bus.

Returns:
- `Int`: Index of the bus in the network.
"""
function geNetBusIdx(; net::Net, busName::String)::Int
  if !haskey(net.busDict, busName)
    error("Bus $(busName) not found in the network")
  end
  return net.busDict[busName]
end

"""
getNetOrigBusIdx: Gets the original index of a bus in the network.

Parameters:
- `net::Net`: Network object.
- `busName::String`: Name of the bus.

Returns:
- `Int`: Original index of the bus in the network.
"""
function getNetOrigBusIdx(; net::Net, busName::String)::Int
  id = geNetBusIdx(net = net, busName = busName)
  return net.busOrigIdxDict[id]
end

function getBusProsumers(net::Net, busIdx::Int)::Vector{ProSumer}
  bus_ps = ProSumer[]
  for ps in net.prosumpsVec
    if getPosumerBusIndex(ps) == busIdx
      push!(bus_ps, ps)
    end
  end
  return bus_ps
end

function getEffectiveBusType(net::Net, busIdx::Int)::NodeType
  has_slack = false
  has_regulating = false

  for ps in net.prosumpsVec
    getPosumerBusIndex(ps) == busIdx || continue
    if isSlack(ps)
      has_slack = true
    end
    if isGenerator(ps) && isRegulating(ps)
      has_regulating = true
    end
  end

  if has_slack
    return Slack
  elseif has_regulating
    return PV
  else
    return PQ
  end
end

function getEffectiveBusType(; net::Net, busName::String)::NodeType
  busIdx = geNetBusIdx(net = net, busName = busName)
  return getEffectiveBusType(net, busIdx)
end

function refreshBusTypesFromProsumers!(net::Net)
  nbus = length(net.nodeVec)
  has_slack = falses(nbus)
  has_regulating = falses(nbus)

  @inbounds for ps in net.prosumpsVec
    bus = getPosumerBusIndex(ps)
    (1 <= bus <= nbus) || continue
    if isSlack(ps)
      has_slack[bus] = true
    end
    if isGenerator(ps) && isRegulating(ps)
      has_regulating[bus] = true
    end
  end

  empty!(net.slackVec)
  @inbounds for busIdx in eachindex(net.nodeVec)
    if net.nodeVec[busIdx]._nodeType == Isolated
      continue
    end
    node_type = has_slack[busIdx] ? Slack : (has_regulating[busIdx] ? PV : PQ)
    net.nodeVec[busIdx]._nodeType = node_type
    if node_type == Slack
      push!(net.slackVec, busIdx)
    end
  end
  return nothing
end

"""
addBus!: Adds a bus to the network.

Parameters:
- `net::Net`: Network object.
- `busName::String`: Name of the bus.
- `vn_kV::Float64`: Nominal voltage of the bus in kV.
- `busType::Union{Nothing,String} = nothing`: Legacy static bus label (deprecated, ignored for PF typing).
- `vm_pu::Float64 = 1.0`: Voltage magnitude of the bus in per unit (default is 1.0).
- `va_deg::Float64 = 0.0`: Voltage angle of the bus in degrees (default is 0.0).
- `vmin_pu::Union{Nothing,Float64} = nothing`: Minimum voltage limit in per unit (default is network's vmin_pu).
- `vmax_pu::Union{Nothing,Float64} = nothing`: Maximum voltage limit in per unit (default is network's vmax_pu).
- `isAux::Bool = false`: Boolean indicating if the bus is auxiliary (default is false).
- `oBusIdx::Union{Nothing,Int} = nothing`: Original bus index (default is nothing).
- `zone::Union{Nothing,Int} = nothing`: Zone index (default is nothing).
- `area::Union{Nothing,Int} = nothing`: Area index (default is nothing).
- `ratedS::Union{Nothing,Float64} = nothing`: Rated power of the bus in MVA (default is nothing).
"""
const _addbus_bus_type_warned = Ref(false)

function addBus!(;
  net::Net,
  busName::String,
  vn_kV::Float64,
  busType::Union{Nothing,String} = nothing,
  vm_pu::Float64 = 1.0,
  va_deg::Float64 = 0.0,
  vmin_pu::Union{Nothing,Float64} = nothing,
  vmax_pu::Union{Nothing,Float64} = nothing,
  isAux::Bool = false,
  oBusIdx::Union{Nothing,Int} = nothing,
  zone::Union{Nothing,Int} = nothing,
  area::Union{Nothing,Int} = nothing,
  ratedS::Union{Nothing,Float64} = nothing,
)
  if !isnothing(busType)
    @assert busType in ["Slack", "SLACK", "PQ", "PV"] "Invalid bus type: $busType"
    if !_addbus_bus_type_warned[]
      msg = "addBus!: `busType` is a legacy input. Power-flow bus typing is derived from attached prosumers."
      if _trace_legacy_bus_type_warnings()
        msg *= " caller=" * _legacy_bus_type_caller_hint()
      end
      @warn msg
      _addbus_bus_type_warned[] = true
    end
  end

  if net._locked
    @error "Network is locked"
  end
  if haskey(net.busDict, busName)
    @error "Bus $busName already exists in the network"
  end
  busIdx = length(net.nodeVec) + 1
  net.busDict[busName] = busIdx

  vmin_pu = isnothing(vmin_pu) ? net.vmin_pu : vmin_pu
  vmax_pu = isnothing(vmax_pu) ? net.vmax_pu : vmax_pu

  if !isnothing(oBusIdx)
    if !haskey(net.busOrigIdxDict, busIdx)
      net.busOrigIdxDict[busIdx] = oBusIdx
    else
      @warn "Original bus index already set for bus $busName"
    end
  end

  static_bus_type = isnothing(busType) ? "PQ" : busType
  if uppercase(static_bus_type) == "SLACK"
    push!(net.slackVec, busIdx)
  end

  node = Node(busIdx = busIdx, vn_kV = vn_kV, nodeType = toNodeType(static_bus_type), vm_pu = vm_pu, va_deg = va_deg, vmin_pu = vmin_pu, vmax_pu = vmax_pu, isAux = isAux, oBusIdx = oBusIdx, zone = zone, area = area, ratedS = ratedS)
  push!(net.nodeVec, node)
end

"""
    addShunt!(; net, busName, pShunt, qShunt, in_service=1)

Adds a *Y-model* shunt to the network.

Semantics:
- `pShunt`, `qShunt` are interpreted as MW/MVar at V = 1.0 pu (MATPOWER-style).
- Internally, the shunt is represented as a pu-admittance stamped into YBUS:
      y_pu = (pShunt + j*qShunt) / baseMVA
  when `bus_shunt_model = "admittance"`.
- When `bus_shunt_model = "voltage_dependent_injection"`, the same
  pu-admittance is excluded from YBUS and contributes
  `-|V|^2 * conj(y_pu)` to the specified net injection.
- The shunt power is *not* constant; it depends on |V|² and will be computed
  after solving via `updateShuntPowers!(net)`.

IMPORTANT:
- This does NOT call `addShuntPower!` or add a constant-power shunt load.
- In voltage-dependent injection mode, the solver adds only the voltage-dependent
  equivalent term to the injection/mismatch path.
"""
function addShunt!(; net::Net, busName::String, pShunt::Float64, qShunt::Float64, in_service::Int = 1, bus_shunt_model = net.bus_shunt_model)
  @assert in_service in (0, 1) "in_service must be 0 or 1"
  shunt_model = normalize_bus_shunt_model(bus_shunt_model)

  busIdx = geNetBusIdx(net = net, busName = busName)
  idShunt = length(net.shuntVec) + 1
  net.shuntDict[busIdx] = idShunt

  vn_kV = getNodeVn(net.nodeVec[busIdx])

  # Build shunt as Y-model:
  # y_pu = (P+jQ)/baseMVA  (P,Q in MW/MVar at 1pu)
  ypu = ComplexF64(pShunt, qShunt) / net.baseMVA

  # Create via constructor (p/q here are ONLY used to set y_pu_shunt in this new semantics)
  sh = Shunt(fromBus = busIdx, id = idShunt, base_MVA = net.baseMVA, vn_kV_shunt = vn_kV, p_shunt = pShunt, q_shunt = qShunt, model = _shunt_component_model(shunt_model), status = in_service)

  # Enforce the intended semantics explicitly (in case constructor changes again)
  sh.y_pu_shunt = ypu
  sh.model = _shunt_component_model(shunt_model)

  # Optional: keep the "G/B" fields consistent with y_pu (pu values)
  sh.G_shunt = real(ypu)
  sh.B_shunt = imag(ypu)

  # NOTE:
  # Do NOT call addShuntPower!(...) here, otherwise you'd treat it as constant PQ.
  # p_shunt/q_shunt will be updated after PF using updateShuntPowers!(net).
  push!(net.shuntVec, sh)

  return nothing
end
"""
    hasShunt!(; net::Net, busName::String)::Bool

Checks if a shunt exists at the specified bus in the network.

# Arguments
- `net::Net`: The network.
- `busName::String`: The name of the bus.

# Returns
- `Bool`: True if a shunt exists at the specified bus, false otherwise.

# Example
```julia
hasShunt!(net = network, busName = "Bus1")
"""
function hasShunt!(; net::Net, busName::String)::Bool
  busIdx = geNetBusIdx(net = net, busName = busName)
  return haskey(net.shuntDict, busIdx)
end

"""
    getShunt!(; net::Net, busName::String)::Shunt

Retrieves the shunt at the specified bus in the network.

# Arguments
- `net::Net`: The network.
- `busName::String`: The name of the bus.

# Returns
- `Shunt`: The shunt at the specified bus.

# Example
```julia
getShunt!(net = network, busName = "Bus1")
"""
function getShunt!(; net::Net, busName::String)::Shunt
  busIdx  = geNetBusIdx(net = net, busName = busName)
  idShunt = net.shuntDict[busIdx]
  return net.shuntVec[idShunt]
end

"""
addBranch!: Adds a branch to the network.

Parameters:
- `net::Net`: Network object.
- `from::Int`: Index of the "from" bus.
- `to::Int`: Index of the "to" bus.
- `branch::AbstractBranch`: Branch object to add.
- `status::Int = 1`: Status of the branch (default is 1).
- `ratio::Union{Nothing,Float64} = nothing`: Ratio of the branch (default is nothing).
- `side::Union{Nothing,Int} = nothing`: Side of the branch (default is nothing).
- `vn_kV::Union{Nothing,Float64} = nothing`: Nominal voltage of the branch in kV (default is nothing).
- `values_are_pu = false`: Boolean indicating if the values are in per unit (default is false).
"""

function addBranch!(; net::Net, from::Int, to::Int, branch::AbstractBranch, status::Integer = 1, ratio = nothing, side = nothing, vn_kV = nothing, values_are_pu::Bool = false, from_status::Union{Nothing,Integer} = nothing, to_status::Union{Nothing,Integer} = nothing)
  @assert from != to "From and to bus must be different"
  idBrunch = length(net.branchVec) + 1
  fOrig = nothing
  tOrig = nothing
  if haskey(net.busOrigIdxDict, from)
    fOrig = net.busOrigIdxDict[from]
  end
  if haskey(net.busOrigIdxDict, to)
    tOrig = net.busOrigIdxDict[to]
  end

  if isnothing(vn_kV)
    vn_kV = getNodeVn(net.nodeVec[from])
  end
  br = Branch(branchIdx = idBrunch, from = from, to = to, baseMVA = net.baseMVA, branch = branch, id = idBrunch, status = status, ratio = ratio, side = side, vn_kV = vn_kV, fromOid = fOrig, toOid = tOrig, values_are_pu = values_are_pu, from_status = from_status, to_status = to_status)
  push!(net.branchVec, br)
end

"""
    addLink!(; net::Net, fromBus::String, toBus::String, status::Int = 1)::Int

Adds an impedance-less topological bus link (e.g. busbar coupler) used for
post-power-flow KCL allocation.
"""
function addLink!(; net::Net, fromBus::String, toBus::String, status::Int = 1)::Int
  @assert !net._locked "Network is locked"
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  @assert from != to "fromBus and toBus must be different"
  @assert !(from in net.slackVec || to in net.slackVec) "Links connected to a slack bus are not allowed"
  @assert getNodeType(net.nodeVec[from]) == getNodeType(net.nodeVec[to]) "Linked buses must have the same bus type"

  idLink = length(net.linkVec) + 1
  push!(net.linkVec, BusLink(linkIdx = idLink, fromBus = from, toBus = to, status = status))
  return idLink
end

"""
    setNetLinkStatus!(; net::Net, linkNr::Int, status::Int)

Sets an existing link status (0/1).
"""
function setNetLinkStatus!(; net::Net, linkNr::Int, status::Int)
  @assert 1 <= linkNr <= length(net.linkVec) "linkNr out of range"
  setLinkStatus!(net.linkVec[linkNr], status)
end

"""
    getNetLinks(; net::Net, fromBus::String, toBus::String)::Vector{BusLink}

Returns all links between two buses (both directions).
"""
function getNetLinks(; net::Net, fromBus::String, toBus::String)::Vector{BusLink}
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  return [l for l in net.linkVec if (l.fromBus == from && l.toBus == to) || (l.fromBus == to && l.toBus == from)]
end

"""
    updateBranchParameters!(;net::Net, fromBus::String, toBus::String, branch::AbstractBranch)

Updates the parameters of a branch in the network.

# Arguments
- `net::Net`: The network.
- `fromBus::String`: The name of the bus where the branch starts.
- `toBus::String`: The name of the bus where the branch ends.
- `branch::BranchModel`: The branch with the updated parameters.

# Example
```julia
updateBranchParameters!(net = network, fromBus = "Bus1", toBus = "Bus2", branch = updatedBranch)
```
"""
function updateBranchParameters!(; net::Net, branchNr::Int, branch::BranchModel)
  br = net.branchVec[branchNr]
  br.r_pu = branch.r_pu
  br.x_pu = branch.x_pu
  # a real equipment edit moves the physical base too (#329), so short circuit
  # and export follow the new model; a later series-FACTS stamp overwrites only
  # the live r_pu/x_pu, leaving this base intact
  br.r_base_pu = branch.r_pu
  br.x_base_pu = branch.x_pu
  br.b_pu = branch.b_pu
  br.g_pu = branch.g_pu
  br.ratio = branch.ratio
  br.angle = branch.angle
  br.tap_ratio = branch.ratio == 0.0 ? 1.0 : branch.ratio
  br.phase_shift_deg = branch.ratio == 0.0 ? 0.0 : branch.angle
  br.sn_MVA = branch.sn_MVA
end

"""
addACLine!: Adds an AC line segment to the network.

Parameters:
- `net::Net`: Network object.
- `fromBus::String`: Name of the "from" bus.
- `toBus::String`: Name of the "to" bus.
- `length::Float64`: Length of the line segment.
- `r::Float64`: Resistance per Meter of the line segment.
- `x::Float64`: Reactance per Meter of the line segment.
- `b::Union{Nothing,Float64} = nothing`: Susceptance per Meter of the line segment (default is nothing).
- `c_nf_per_km::Union{Nothing,Float64} = nothing`: Capacitance per Meter of the line segment in nF/km (default is nothing).
- `tanδ::Union{Nothing,Float64} = nothing`: Tangent of the loss angle (default is nothing).
- `ratedS::Union{Nothing, Float64}= nothing`: Rated power of the line segment in MVA (default is nothing).
- `status::Int = 1`: Status of the line segment (default is 1).
"""
function addACLine!(; net::Net, fromBus::String, toBus::String, length::Float64, r::Float64, x::Float64, b::Union{Nothing,Float64} = nothing, c_nf_per_km::Union{Nothing,Float64} = nothing, tanδ::Union{Nothing,Float64} = nothing, ratedS::Union{Nothing,Float64} = nothing, status::Int = 1, from_status::Union{Nothing,Integer} = nothing, to_status::Union{Nothing,Integer} = nothing)
  @assert fromBus != toBus "From and to bus must be different"
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  vn_kV = getNodeVn(net.nodeVec[from])
  vn_2_kV = getNodeVn(net.nodeVec[to])
  @assert vn_kV == vn_2_kV "Voltage level of the from bus $(vn_kV) does not match the to bus $(vn_2_kV)"
  if length > 1.0
    _par_length = true
  else
    _par_length = false
  end
  acseg = ACLineSegment(vn_kv = vn_kV, from = from, to = to, length = length, r = r, x = x, b = b, c_nf_per_km = c_nf_per_km, tanδ = tanδ, ratedS = ratedS, paramsBasedOnLength = _par_length, isPIModel = false)
  push!(net.linesAC, acseg)

  addBranch!(net = net, from = from, to = to, branch = acseg, vn_kV = vn_kV, status = status, values_are_pu = true, from_status = from_status, to_status = to_status)
end

"""
    addPIModelACLine!(; net::Net, fromBus::String, toBus::String, r_pu::Float64, x_pu::Float64, b_pu::Float64, g_pu::Union{Nothing,Float64}=nothing, status::Int, ratedS::Union{Nothing,Float64}=nothing)

Adds a PI model AC line to the network.

# Arguments
- `net::Net`: The network.
- `fromBus::String`: The name of the bus where the line starts.
- `toBus::String`: The name of the bus where the line ends.
- `r_pu::Float64`: The per unit resistance of the line.
- `x_pu::Float64`: The per unit reactance of the line.
- `b_pu::Float64`: The per unit total line charging susceptance of the line.
- `g_pu::Union{Nothing,Float64}`: The per unit total shunt conductance of the line (e.g. CGMES `gch` after conversion). Default is `nothing` (treated as 0.0).
- `status::Int`: The status of the line. 1 = in service, 0 = out of service.
- `ratedS::Union{Nothing,Float64}`: The rated power of the line.

# Example
```julia
addPIModelACLine!(net = network, fromBus = "Bus1", toBus = "Bus2", r_pu = 0.01, x_pu = 0.1, b_pu = 0.02, status = 1, ratedS = 100.0)
```
"""
function _addPIModelACLine_by_idx!(; net::Net, from::Int, to::Int, r_pu::Float64, x_pu::Float64, b_pu::Float64, g_pu::Union{Nothing,Float64} = nothing, status::Int, ratedS::Union{Nothing,Float64} = nothing, from_status::Union{Nothing,Integer} = nothing, to_status::Union{Nothing,Integer} = nothing)
  @assert from != to "From and to bus must be different"
  vn_kV = getNodeVn(net.nodeVec[from])
  vn_2_kV = getNodeVn(net.nodeVec[to])
  @assert vn_kV == vn_2_kV "Voltage level of the from bus $(vn_kV) does not match the to bus $(vn_2_kV)"
  acseg = ACLineSegment(vn_kv = vn_kV, from = from, to = to, length = 1.0, r = r_pu, x = x_pu, b = b_pu, g = g_pu, ratedS = ratedS, paramsBasedOnLength = false, isPIModel = true)
  push!(net.linesAC, acseg)

  addBranch!(net = net, from = from, to = to, branch = acseg, vn_kV = vn_kV, status = status, values_are_pu = true, from_status = from_status, to_status = to_status)
end

function addPIModelACLine!(; net::Net, fromBus::String, toBus::String, r_pu::Float64, x_pu::Float64, b_pu::Float64, g_pu::Union{Nothing,Float64} = nothing, status::Int, ratedS::Union{Nothing,Float64} = nothing, from_status::Union{Nothing,Integer} = nothing, to_status::Union{Nothing,Integer} = nothing)
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  return _addPIModelACLine_by_idx!(net = net, from = from, to = to, r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, g_pu = g_pu, status = status, ratedS = ratedS, from_status = from_status, to_status = to_status)
end

"""
Add a transformer with PI model to the network.

# Arguments
- `net::Net`: The network to which the transformer will be added.
- `fromBus::String`: The name of the bus where the transformer originates.
- `toBus::String`: The name of the bus where the transformer terminates.
- `r_pu::Float64`: The per-unit resistance of the transformer.
- `x_pu::Float64`: The per-unit reactance of the transformer.
- `b_pu::Float64`: The per-unit susceptance of the transformer.
- `status::Int`: The status of the transformer.
- `ratedU::Union{Nothing, Float64}`: Rated voltage of the transformer. Default is `nothing`.
- `ratedS::Union{Nothing, Float64}`: Rated apparent power of the transformer. Default is `nothing`.
- `ratio::Union{Nothing, Float64}`: Ratio of the transformer. Default is `nothing`.
- `shift_deg::Union{Nothing, Float64}`: Phase shift angle of the transformer. Default is `nothing`.
- `isAux::Bool`: Whether the transformer is an auxiliary transformer. Default is `false`.
"""
function _addPIModelTrafo_by_idx!(;
  net::Net,
  from::Int,
  to::Int,
  r_pu::Float64,
  x_pu::Float64,
  b_pu::Float64,
  g_pu::Float64 = 0.0,
  status::Int,
  ratedU::Union{Nothing,Float64} = nothing,
  ratedS::Union{Nothing,Float64} = nothing,
  ratio::Union{Nothing,Float64} = nothing,
  shift_deg::Union{Nothing,Float64} = nothing,
  isAux::Bool = false,
  side::Int = 1,
  controls::Union{Nothing,Vector{PowerTransformerControl}} = nothing,
  from_status::Union{Nothing,Integer} = nothing,
  to_status::Union{Nothing,Integer} = nothing,
)
  @assert from != to "From and to bus must be different"
  vn_hv_kV = getNodeVn(net.nodeVec[from])
  vn_lv_kV = getNodeVn(net.nodeVec[to])
  w1 = PowerTransformerWinding(vn_hv_kV, r_pu, x_pu, b_pu, g_pu, ratio, shift_deg, ratedU, ratedS, nothing, true)
  w2 = PowerTransformerWinding(vn_lv_kV, 0.0, 0.0)
  if !isnothing(controls)
    if side == 1
      w1.controls = copy(controls)
    elseif side == 2
      w2.controls = copy(controls)
    end
  end
  comp = getTrafoImpPGMComp(isAux, vn_hv_kV, from, to)
  trafo = PowerTransformer(comp, false, w1, w2, nothing, Sparlectra.PIModel)
  push!(net.trafos, trafo)

  addBranch!(net = net, from = from, to = to, branch = trafo, status = status, ratio = ratio, side = side, vn_kV = vn_hv_kV, values_are_pu = true, from_status = from_status, to_status = to_status)
  if !isnothing(controls)
    br = net.branchVec[end]
    target_controls = side == 1 ? net.trafos[end].side1.controls : net.trafos[end].side2.controls
    for ctrl in target_controls
      ctrl.trafo = string(br.branchIdx)
    end
  end
end

function addPIModelTrafo!(;
  net::Net,
  fromBus::String,
  toBus::String,
  r_pu::Float64,
  x_pu::Float64,
  b_pu::Float64,
  status::Int,
  ratedU::Union{Nothing,Float64} = nothing,
  ratedS::Union{Nothing,Float64} = nothing,
  ratio::Union{Nothing,Float64} = nothing,
  shift_deg::Union{Nothing,Float64} = nothing,
  isAux::Bool = false,
  side::Int = 1,
  controls::Union{Nothing,Vector{PowerTransformerControl}} = nothing,
  from_status::Union{Nothing,Integer} = nothing,
  to_status::Union{Nothing,Integer} = nothing,
)
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  return _addPIModelTrafo_by_idx!(net = net, from = from, to = to, r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, status = status, ratedU = ratedU, ratedS = ratedS, ratio = ratio, shift_deg = shift_deg, isAux = isAux, side = side, controls = controls, from_status = from_status, to_status = to_status)
end

"""
Add a transformer with PI model to the network.

# Arguments
- `net::Net`: The network to which the transformer will be added.
- `fromBus::String`: The name of the bus where the transformer originates.
- `toBus::String`: The name of the bus where the transformer terminates.
- `r_pu::Float64`: The per-unit resistance of the transformer.
- `x_pu::Float64`: The per-unit reactance of the transformer.
- `b_pu::Float64`: The per-unit susceptance of the transformer.
- `status::Int`: The status of the transformer.
- `ratedU::Union{Nothing, Float64}`: Rated voltage of the transformer. Default is `nothing`.
- `ratedS::Union{Nothing, Float64}`: Rated apparent power of the transformer. Default is `nothing`.
- `ratio::Union{Nothing, Float64}`: Ratio of the transformer. Default is `nothing`.
- `shift_deg::Union{Nothing, Float64}`: Phase shift angle of the transformer. Default is `nothing`.
- `isAux::Bool`: Whether the transformer is an auxiliary transformer. Default is `false`.
"""
function add2WTPIModelTrafo!(;
  net::Net,
  fromBus::String,
  toBus::String,
  side::Int = 1,
  r::Float64,
  x::Float64,
  b::Float64 = 0.0,
  status::Int = 1,
  ratedU::Union{Nothing,Float64} = nothing,
  ratedS::Union{Nothing,Float64} = nothing,
  ratio::Union{Nothing,Float64} = nothing,
  shift_deg::Union{Nothing,Float64} = nothing,
)
  @assert fromBus != toBus "From and to bus must be different"
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  if isnothing(ratedU)
    ratedU = side == 1 ? getNodeVn(net.nodeVec[from]) : getNodeVn(net.nodeVec[to])
  end
  r_pu, x_pu, b_pu, g_pu = toPU_RXBG(r = r, x = x, g = 0.0, b = b, v_kv = ratedU, baseMVA = net.baseMVA)
  addPIModelTrafo!(net = net, fromBus = fromBus, toBus = toBus, r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, status = status, ratedU = ratedU, ratedS = ratedS, ratio = ratio, shift_deg = shift_deg, isAux = false, side = side)
end

"""
Add a 3-winding transformer using a star-equivalent with an internal AUX bus.

Implementation strategy:
- Ensure an AUX bus exists (PQ, isAux=true) at the HV-side nominal voltage.
- Add three 2-winding PI-model transformers:
    AUX -- HB, AUX -- MB, AUX -- LV
- Convert r/x/b to PU (using toPU_RXBG) for each branch (mainly for validation/logging).
  The actual insertion uses add2WTPIModelTrafo!, which performs the conversion internally.

Notes:
- `ratio` is set to U_aux / U_side (HV/MV/LV) for MV and LV, and 1.0 for HB.
- `ratedU` passed to add2WTPIModelTrafo! is the AUX-side rated voltage (HV), because the
  branch is defined from AUX (HV base) to the respective side.
"""
function add3WTPiModelTrafo!(; net::Net, HBBus::String, MBBus::String, LVBus::String, r::NTuple{3,Float64}, x::NTuple{3,Float64}, b::NTuple{3,Float64}, ratedU_kV::NTuple{3,Float64}, ratedS_MVA::NTuple{3,Float64}, status::Int = 1)
  @assert status in (0, 1) "status must be 0 or 1"

  # --- basic existence checks (will throw from geNetBusIdx if missing) ---
  hb_idx = geNetBusIdx(net = net, busName = HBBus)
  mb_idx = geNetBusIdx(net = net, busName = MBBus)
  lv_idx = geNetBusIdx(net = net, busName = LVBus)

  # --- choose AUX-side (HV) voltage base ---
  # Prefer ratedU_kV[1] as HV side, but be robust and take the maximum.
  U_aux_kV = maximum(ratedU_kV)
  @assert U_aux_kV > 0.0 "ratedU_kV must be > 0"

  # --- build a deterministic AUX bus name (unique-ish, stable) ---
  # Keep it simple and reproducible; avoid special chars.
  function _sanitize(s::AbstractString)
    return replace(String(s), r"[^A-Za-z0-9_]" => "_")
  end

  aux_bus = "Aux3WT_" * _sanitize(HBBus) * "_" * _sanitize(MBBus) * "_" * _sanitize(LVBus)

  # --- create AUX bus if missing ---
  if !hasBusInNet(net = net, busName = aux_bus)
    addBus!(net = net, busName = aux_bus, vn_kV = U_aux_kV, isAux = true)
  end

  # --- define ratios for each branch AUX->side ---
  # ratio is interpreted as turns ratio; here we use U_aux / U_side.
  # For the HV bus branch, this is typically ~1 if HB is the HV side.

  ratio_hb = 1.0 #U_aux_kV / ratedU_kV[1]
  ratio_mb = 1.0 #U_aux_kV / ratedU_kV[2]
  ratio_lv = 1.0 #U_aux_kV / ratedU_kV[3]

  # --- optional PU conversion (for validation/logging/debug) ---
  # Note: add2WTPIModelTrafo! will do this conversion again internally.
  r1_pu, x1_pu, b1_pu, _g1_pu = toPU_RXBG(r = r[1], x = x[1], g = 0.0, b = b[1], v_kv = U_aux_kV, baseMVA = net.baseMVA)
  r2_pu, x2_pu, b2_pu, _g2_pu = toPU_RXBG(r = r[2], x = x[2], g = 0.0, b = b[2], v_kv = U_aux_kV, baseMVA = net.baseMVA)
  r3_pu, x3_pu, b3_pu, _g3_pu = toPU_RXBG(r = r[3], x = x[3], g = 0.0, b = b[3], v_kv = U_aux_kV, baseMVA = net.baseMVA)

  # If you want: add @debug lines here using (r*_pu, x*_pu, b*_pu)

  # --- add three 2W PI branches (AUX -> side buses) ---
  # Use side=1 because fromBus is the AUX (HV base).
  add2WTPIModelTrafo!(net = net, fromBus = aux_bus, toBus = HBBus, side = 1, r = r[1], x = x[1], b = b[1], status = status, ratedU = U_aux_kV, ratedS = ratedS_MVA[1], ratio = ratio_hb, shift_deg = 0.0)

  add2WTPIModelTrafo!(net = net, fromBus = aux_bus, toBus = MBBus, side = 1, r = r[2], x = x[2], b = b[2], status = status, ratedU = U_aux_kV, ratedS = ratedS_MVA[2], ratio = ratio_mb, shift_deg = 0.0)

  add2WTPIModelTrafo!(net = net, fromBus = aux_bus, toBus = LVBus, side = 1, r = r[3], x = x[3], b = b[3], status = status, ratedU = U_aux_kV, ratedS = ratedS_MVA[3], ratio = ratio_lv, shift_deg = 0.0)

  return aux_bus
end


"""
Add a two-winding transformer to the network.

# Arguments
- `net::Net`: The network to which the transformer will be added.
- `fromBus::String`: The name of the bus where the transformer originates.
- `toBus::String`: The name of the bus where the transformer terminates.
- `sn_mva::Float64`: Rated power of the transformer.
- `vk_percent::Float64`: Voltage regulation percent of the transformer.
- `vkr_percent::Float64`: Voltage regulation percent of the transformer.
- `pfe_kw::Float64`: Iron loss of the transformer.
- `i0_percent::Float64`: No-load current percent of the transformer.
- `status::Int`: The status of the transformer. Default is 1.
"""
function add2WTrafo!(; net::Net, fromBus::String, toBus::String, sn_mva::Float64, vk_percent::Float64, vkr_percent::Float64, pfe_kw::Float64, i0_percent::Float64, status::Int = 1, controls::Union{Nothing,Vector{PowerTransformerControl}} = nothing, from_status::Union{Nothing,Integer} = nothing, to_status::Union{Nothing,Integer} = nothing)
  @assert fromBus != toBus "From and to bus must be different"
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  vn_hv_kV = getNodeVn(net.nodeVec[from])
  vn_lv_kV = getNodeVn(net.nodeVec[to])

  trafo = create2WTRatioTransformerNoTaps(from = from, to = to, vn_hv_kV = vn_hv_kV, vn_lv_kV = vn_lv_kV, sn_mva = sn_mva, vk_percent = vk_percent, vkr_percent = vkr_percent, pfe_kw = pfe_kw, i0_percent = i0_percent)
  side = getSideNumber2WT(trafo)
  if !isnothing(controls)
    if side == 1
      trafo.side1.controls = copy(controls)
    elseif side == 2
      trafo.side2.controls = copy(controls)
    end
  end
  ratio = calcTransformerRatio(trafo)
  push!(net.trafos, trafo)

  addBranch!(net = net, from = from, to = to, branch = trafo, status = status, ratio = ratio, side = side, vn_kV = vn_hv_kV, from_status = from_status, to_status = to_status)
  if !isnothing(controls)
    br = net.branchVec[end]
    target_controls = side == 1 ? net.trafos[end].side1.controls : net.trafos[end].side2.controls
    for ctrl in target_controls
      ctrl.trafo = string(br.branchIdx)
    end
  end
end

"""
Add a prosumer (combination of a producer and consumer) to the network.

# Arguments
- `net::Net`: The network to which the prosumer will be added.
- `busName::String`: The name of the bus where the prosumer is connected.
- `type::String`: The type of the prosumer.
- `p::Union{Nothing, Float64}`: Active power produced or consumed. Default is `nothing`.
- `q::Union{Nothing, Float64}`: Reactive power produced or consumed. Default is `nothing`.
- `pMin::Union{Nothing, Float64}`: Minimum active power. Default is `nothing`.
- `pMax::Union{Nothing, Float64}`: Maximum active power. Default is `nothing`.
- `qMin::Union{Nothing, Float64}`: Minimum reactive power. Default is `nothing`.
- `qMax::Union{Nothing, Float64}`: Maximum reactive power. Default is `nothing`.
- `referencePri::Union{Nothing, String}`: Reference bus for the prosumer. Default is `nothing`.
- `vm_pu::Union{Nothing, Float64}`: Voltage magnitude setpoint. Default is `nothing`.
- `va_deg::Union{Nothing, Float64}`: Voltage angle setpoint. Default is `nothing`.
- `isRegulated::Bool`: Marks a prosumer as voltage-regulating for PV bus resolution. Default is `false`.
"""
function addProsumer!(;
  net::Net,
  busName::String,
  type::String,
  p::Union{Nothing,Float64} = nothing,
  q::Union{Nothing,Float64} = nothing,
  pMin::Union{Nothing,Float64} = nothing,
  pMax::Union{Nothing,Float64} = nothing,
  qMin::Union{Nothing,Float64} = nothing,
  qMax::Union{Nothing,Float64} = nothing,
  referencePri::Union{Nothing,String} = nothing,
  vm_pu::Union{Nothing,Float64} = nothing,
  va_deg::Union{Nothing,Float64} = nothing,
  vstep_pu::Union{Nothing,Float64} = nothing,
  tap_steps_down::Union{Nothing,Int} = nothing,
  tap_steps_up::Union{Nothing,Int} = nothing,
  qu_controller::Union{Nothing,QUController} = nothing,
  pu_controller::Union{Nothing,PUController} = nothing,
  isRegulated::Bool = false,
  participationFactor::Union{Nothing,Float64} = nothing,
  defer_bus_type_refresh::Bool = false,
)
  busIdx = geNetBusIdx(net = net, busName = busName)
  idProsSum = length(net.prosumpsVec) + 1
  isAPUNode = false
  if !isnothing(vm_pu)
    node = net.nodeVec[busIdx]
    isAPUNode = isPVNode(node)
    nodeVm = getNodeVm(node)
    if !isnothing(nodeVm) && abs(nodeVm - vm_pu) > 1e-6
      @debug "Voltage setpoint already present at bus $busName: keep vm=$(nodeVm), ignore vm=$(vm_pu)"
    else
      setVmVa!(node = node, vm_pu = vm_pu, va_deg = va_deg)
    end
  end
  proTy = toProSumptionType(type)
  refPriIdx = isnothing(referencePri) ? nothing : geNetBusIdx(net = net, busName = referencePri)
  vn_kV = getNodeVn(net.nodeVec[busIdx])
  has_vm_setpoint = !isnothing(vm_pu)
  has_vstep_with_taps = !isnothing(vstep_pu) && (!isnothing(tap_steps_down) || !isnothing(tap_steps_up))
  auto_regulated = isRegulated || has_vm_setpoint || has_vstep_with_taps
  vset_adjust_cfg = if !isnothing(vstep_pu) && !isnothing(tap_steps_down) && !isnothing(tap_steps_up)
    VoltageAdjustConfig(vstep_pu, tap_steps_down, tap_steps_up)
  else
    nothing
  end
  ps = ProSumer(
    vn_kv = vn_kV,
    busIdx = busIdx,
    oID = idProsSum,
    type = proTy,
    p = p,
    q = q,
    minP = pMin,
    maxP = pMax,
    minQ = qMin,
    maxQ = qMax,
    referencePri = refPriIdx,
    vm_pu = vm_pu,
    va_deg = va_deg,
    vstep_pu = vstep_pu,
    tap_steps_down = tap_steps_down,
    tap_steps_up = tap_steps_up,
    vset_adjust = vset_adjust_cfg,
    isRegulated = auto_regulated,
    isAPUNode = isAPUNode,
    quController = qu_controller,
    puController = pu_controller,
    participationFactor = participationFactor,
  )
  push!(net.prosumpsVec, ps)
  node = net.nodeVec[busIdx]
  if (isGenerator(proTy))
    addGenPower!(node = node, p = p, q = q)
  else
    addLoadPower!(node = node, p = p, q = q)
  end

  if !defer_bus_type_refresh
    refreshBusTypesFromProsumers!(net)
    _buildQLimits!(net)
  end
end

"""
    addExternalGrid!(; net, busName, sk_max_MVA, kwargs...)

Add an external grid (IEC 60909-0 network feeder, issue #299): the ideal
voltage source of the superordinate network plus its declared short-circuit
data.

Load-flow side: by default the connection bus becomes the reference (REF)
bus, exactly as a manually added slack-type prosumer does — the external
grid is ideal in the power flow, and the short-circuit attributes change
**no** power-flow result (`participationFactor` precedent). With
`internal_impedance = true` the source instead becomes non-ideal: a hidden
auxiliary bus `<busName>__extgrid_int` carries the reference voltage, and a
series branch with `z_pu = baseMVA / sk_max_MVA` (voltage factor `c = 1` —
the c-factor is a short-circuit concept), split by `rx_max`, connects it to
the terminal bus. The terminal bus then stays an ordinary solved bus. The
auxiliary bus is tagged `isAux`; expect it in reports under its generated
name. Very large `sk_max_MVA` values degrade Jacobian conditioning — use
the default ideal representation when ideal behavior is wanted.

Short-circuit side: `sk_max_MVA` (and optionally `sk_min_MVA`) with the
`rx_max`/`rx_min` ratios are converted **at add time** into the feeder
record consumed by [`runShortCircuit!`](@ref) — anchored at the physical
connection bus also in the `internal_impedance` variant (the auxiliary
branch is a dead end in the short-circuit network and carries no fault
current).

# Arguments
- `net::Net`: the network.
- `busName::String`: connection bus (must exist).
- `vm_pu::Float64 = 1.0`, `va_deg::Float64 = 0.0`: reference voltage.
- `sk_max_MVA::Float64`: declared initial symmetrical short-circuit power,
  maximum case. Mandatory — without it the feeder cannot contribute to any
  short circuit.
- `sk_min_MVA::Union{Nothing,Float64} = nothing`: minimum case; without it
  the `:min` case skips the feeder with the engine's safety flag.
- `rx_max::Float64 = 0.1`: R/X ratio, maximum case (IEC 60909-0 default).
- `rx_min::Union{Nothing,Float64} = nothing`: R/X ratio, minimum case.
  Defaults to `rx_max` when `sk_min_MVA` is given, so a deliberately
  declared minimum feeder does not carry a spurious defaulted-data flag.
- `name::Union{Nothing,String} = nothing`: feeder name (default `busName`).
- `internal_impedance::Bool = false`: non-ideal load-flow variant (above).

Failure behavior: throws `ArgumentError` for a non-positive/non-finite
`sk_max_MVA`, `sk_min_MVA > sk_max_MVA`, or negative R/X ratios. Multiple
external grids on one bus are allowed and stack as parallel feeders.
"""
function addExternalGrid!(;
  net::Net,
  busName::String,
  vm_pu::Float64 = 1.0,
  va_deg::Float64 = 0.0,
  sk_max_MVA::Float64,
  sk_min_MVA::Union{Nothing,Float64} = nothing,
  rx_max::Float64 = 0.1,
  rx_min::Union{Nothing,Float64} = nothing,
  name::Union{Nothing,String} = nothing,
  internal_impedance::Bool = false,
)
  busIdx = geNetBusIdx(net = net, busName = busName)
  (isfinite(sk_max_MVA) && sk_max_MVA > 0.0) || throw(ArgumentError("addExternalGrid!: sk_max_MVA must be finite and > 0; got $(sk_max_MVA)."))
  if sk_min_MVA !== nothing
    (isfinite(sk_min_MVA) && 0.0 < sk_min_MVA <= sk_max_MVA) || throw(ArgumentError("addExternalGrid!: sk_min_MVA must be finite and within (0, sk_max_MVA]; got $(sk_min_MVA) with sk_max_MVA = $(sk_max_MVA)."))
  end
  (isfinite(rx_max) && rx_max >= 0.0) || throw(ArgumentError("addExternalGrid!: rx_max must be finite and >= 0; got $(rx_max)."))
  if rx_min !== nothing
    (isfinite(rx_min) && rx_min >= 0.0) || throw(ArgumentError("addExternalGrid!: rx_min must be finite and >= 0; got $(rx_min)."))
  end

  # Task decision (issue #299): a declared minimum feeder without its own
  # ratio inherits rx_max instead of `nothing`, so the :min case does not
  # flag "no usable R/X ratio" for data the user deliberately provided.
  eff_rx_min = rx_min !== nothing ? rx_min : (sk_min_MVA !== nothing ? rx_max : nothing)

  if internal_impedance
    # Non-ideal variant (Stage 2): the reference voltage moves to a hidden
    # internal bus — the actual slack — behind z = c·Un²/Sk with c = 1. On
    # the per-voltage-level impedance base (z_base = Un²/baseMVA) that is
    # simply z_pu = baseMVA/Sk, independent of Un.
    internalName = busName * "__extgrid_int"
    haskey(net.busDict, internalName) && throw(ArgumentError("addExternalGrid!: internal bus $(internalName) already exists — one internal-impedance external grid per bus."))
    vn_kV = getNodeVn(net.nodeVec[busIdx])
    addBus!(net = net, busName = internalName, vn_kV = vn_kV, vm_pu = vm_pu, va_deg = va_deg, isAux = true)
    z_pu = net.baseMVA / sk_max_MVA
    x_pu = z_pu / sqrt(1.0 + rx_max^2)
    r_pu = rx_max * x_pu
    addPIModelACLine!(net = net, fromBus = internalName, toBus = busName, r_pu = r_pu, x_pu = x_pu, b_pu = 0.0, status = 1)
    addProsumer!(net = net, busName = internalName, type = "EXTERNALNETWORKINJECTION", referencePri = internalName, vm_pu = vm_pu, va_deg = va_deg)
  else
    # Ideal variant (Stage 1): exactly the manual slack path — the prosumer
    # with referencePri marks the connection bus REF via
    # refreshBusTypesFromProsumers!.
    addProsumer!(net = net, busName = busName, type = "EXTERNALNETWORKINJECTION", referencePri = busName, vm_pu = vm_pu, va_deg = va_deg)
  end

  # Feeder record for the short-circuit engine: full CGMES ENI tuple contract
  # (all eleven keys, unknowns as `nothing`). `bus` must be the busDict key —
  # the engine resolves it through net.busDict. The mrid is deterministic and
  # unique per bus; the stamping falls back to it via something(name, mrid),
  # so it must never be nothing. Parallel feeders continue after the highest
  # existing suffix (not the record count) so the id stays unique even after
  # a record is removed and another feeder added.
  ik_max_A = sk_max_MVA * 1.0e6 / (sqrt(3.0) * getNodeVn(net.nodeVec[busIdx]) * 1.0e3)
  ik_min_A = sk_min_MVA === nothing ? nothing : sk_min_MVA * 1.0e6 / (sqrt(3.0) * getNodeVn(net.nodeVec[busIdx]) * 1.0e3)
  mrid_base = "native-eni-$(busName)"
  max_suffix = 0
  for e in net.sc_sources.external_network_injections
    e.bus == busName || continue
    m = String(e.mrid)
    startswith(m, mrid_base) || continue
    rest = m[(lastindex(mrid_base) + 1):end]
    n = isempty(rest) ? 1 : something(tryparse(Int, lstrip(rest, '-')), 0)
    max_suffix = max(max_suffix, n)
  end
  mrid = max_suffix == 0 ? mrid_base : "$(mrid_base)-$(max_suffix + 1)"
  push!(
    net.sc_sources.external_network_injections,
    (
      mrid = mrid,
      name = name === nothing ? busName : name,
      bus = busName,
      maxInitialSymShCCurrent_A = ik_max_A,
      minInitialSymShCCurrent_A = ik_min_A,
      maxR1ToX1Ratio = rx_max,
      minR1ToX1Ratio = eff_rx_min,
      maxR0ToX0Ratio = nothing,
      maxZ0ToZ1Ratio = nothing,
      ikSecond = nothing,
      governorSCD = nothing,
    ),
  )
  return nothing
end

"""
    convertSlackToExternalGrid!(; net, sk_max_MVA, kwargs...) -> String

Replace the marked slack bus by a non-ideal external-grid source (issue
#299): the slack marker moves to a hidden internal bus behind the feeder
impedance, so the former slack bus becomes an ordinary solved bus whose
voltage reacts to loading.

Steps: the slack prosumers at the bus are demoted (reference marker and
voltage regulation removed — their scheduled `p`/`q` injections stay), and
`addExternalGrid!` with `internal_impedance = true` is added at the bus,
carrying the former slack's voltage setpoint as the source's reference
voltage. With multiple slack buses (island references) only the addressed
bus is converted — every other island keeps its reference.

# Arguments
- `net::Net`: the network.
- `busName::Union{Nothing,String} = nothing`: slack bus to convert; default
  is the primary (first registered) slack bus.
- `sk_max_MVA`, `sk_min_MVA`, `rx_max`, `rx_min`, `name`: forwarded to
  [`addExternalGrid!`](@ref).

Returns a short human-readable note describing the conversion (for run
logs). Failure behavior: throws `ArgumentError` when the net has no slack
bus or `busName` is not a slack bus.
"""
function convertSlackToExternalGrid!(;
  net::Net,
  busName::Union{Nothing,String} = nothing,
  sk_max_MVA::Float64,
  sk_min_MVA::Union{Nothing,Float64} = nothing,
  rx_max::Float64 = 0.1,
  rx_min::Union{Nothing,Float64} = nothing,
  name::Union{Nothing,String} = nothing,
)
  isempty(net.slackVec) && throw(ArgumentError("convertSlackToExternalGrid!: the network has no slack bus to convert."))
  busIdx = busName === nothing ? net.slackVec[1] : geNetBusIdx(net = net, busName = busName)
  busIdx in net.slackVec || throw(ArgumentError("convertSlackToExternalGrid!: bus $(busName) is not a slack bus."))
  bus_names = Dict{Int,String}(idx => n for (n, idx) in net.busDict)
  target = bus_names[busIdx]

  # The former slack's voltage setpoint becomes the source's reference
  # voltage; fall back to the node's start voltage when no prosumer carries
  # a setpoint (the ProSumer constructor defaults vm to 1.0 anyway).
  vref = 1.0
  varef = 0.0
  demoted = 0
  for ps in net.prosumpsVec
    (isSlack(ps) && getPosumerBusIndex(ps) == busIdx) || continue
    ps.vm_pu === nothing || (vref = ps.vm_pu)
    ps.va_deg === nothing || (varef = ps.va_deg)
    # Demote: no reference marker, no voltage holding — otherwise the bus
    # would stay PV-held and the source's finite stiffness could never show.
    ps.referencePri = nothing
    ps.isRegulated = false
    demoted += 1
  end
  # Re-derive bus types and Q-limits right after the demotion instead of
  # relying on the addExternalGrid! → addProsumer! call below to do it as a
  # side effect — the demotion must leave a consistent net on its own.
  refreshBusTypesFromProsumers!(net)
  _buildQLimits!(net)

  addExternalGrid!(net = net, busName = target, vm_pu = vref, va_deg = varef, sk_max_MVA = sk_max_MVA, sk_min_MVA = sk_min_MVA, rx_max = rx_max, rx_min = rx_min, name = name, internal_impedance = true)
  return "external grid: slack bus $(target) converted to a non-ideal source (Sk'' = $(sk_max_MVA) MVA, R/X = $(rx_max), $(demoted) slack prosumer(s) demoted; reference moved to $(target)__extgrid_int)"
end

"""
    _apply_external_grid_config!(net, pf::PowerFlowConfig; declared = nothing) -> Union{Nothing,String}

Apply `power_flow.external_grid` to a freshly imported net: when enabled,
convert the marked slack bus into a non-ideal source via
[`convertSlackToExternalGrid!`](@ref). `declared` optionally carries
`(sk_MVA, rx)` read from the case data (CGMES `ExternalNetworkInjection`);
it wins over the config numbers when `source = :auto`. Idempotent: a net
that already contains an external-grid internal bus is left untouched, so
re-runs and rescue retries cannot stack sources. Returns the conversion
note (or `nothing` when nothing was done).
"""
function _apply_external_grid_config!(net::Net, pf; declared::Union{Nothing,NamedTuple} = nothing)::Union{Nothing,String}
  eg = pf.external_grid
  eg.enabled || return nothing
  any(endswith(n, "__extgrid_int") for n in keys(net.busDict)) && return nothing
  if isempty(net.slackVec)
    @warn "power_flow.external_grid.enabled: the net registers no slack bus — nothing to convert (enable power_flow.auto_slack or mark a slack first)."
    return nothing
  end
  use_declared = eg.source === :auto && declared !== nothing
  sk = use_declared ? declared.sk_MVA : eg.sk_MVA
  rx = use_declared ? declared.rx : eg.rx
  note = convertSlackToExternalGrid!(net = net, sk_max_MVA = sk, rx_max = rx)
  note *= use_declared ? " [Sk''/RX declared by the case data]" : " [Sk''/RX from power_flow.external_grid]"
  @info note
  return note
end

"""
Add active/reactive power to the NODE-level generation sum of a bus
(`node._pƩGen`/`node._qƩGen`).

!!! warning "Report layer only — no solver reads this"
    Every solver builds its injections from the PROSUMER objects
    (`buildComplexSVec` reads `net.prosumpsVec`; since 0.9.12 the DC solver
    uses the same source, issue #323), so this edit does NOT change any
    solve. The node sums feed reporting, MATPOWER export, and measurement
    generation. To change what the solvers compute, edit the prosumer
    (e.g. `ps.pVal`/`ps.qVal`) or add one (`addProsumer!`).

# Arguments
- `net::Net`: The network object.
- `busName::String`: The name of the bus.
- `p::Union{Nothing, Float64}`: active power to add in MW. Default `nothing`.
- `q::Union{Nothing, Float64}`: reactive power to add in MVAr. Default `nothing`.
"""
function addBusGenPower!(; net::Net, busName::String, p::Union{Nothing,Float64} = nothing, q::Union{Nothing,Float64} = nothing)
  busIdx = geNetBusIdx(net = net, busName = busName)
  addGenPower!(node = net.nodeVec[busIdx], p = p, q = q)
end

"""
Add active/reactive power to the NODE-level load sum of a bus
(`node._pƩLoad`/`node._qƩLoad`).

!!! warning "Report layer only — no solver reads this"
    Every solver builds its injections from the PROSUMER objects
    (`buildComplexSVec` reads `net.prosumpsVec`; since 0.9.12 the DC solver
    uses the same source, issue #323), so this edit does NOT change any
    solve. The node sums feed reporting, MATPOWER export, and measurement
    generation. To add load the solvers see, add a prosumer instead:
    `addProsumer!(net = net, busName = ..., type = "LOAD", p = ..., q = ...)`.

# Arguments
- `net::Net`: The network object.
- `busName::String`: The name of the bus.
- `p::Union{Nothing, Float64}`: active power to add in MW. Default `nothing`.
- `q::Union{Nothing, Float64}`: reactive power to add in MVAr. Default `nothing`.
"""
function addBusLoadPower!(; net::Net, busName::String, p::Union{Nothing,Float64} = nothing, q::Union{Nothing,Float64} = nothing)
  busIdx = geNetBusIdx(net = net, busName = busName)
  addLoadPower!(node = net.nodeVec[busIdx], p = p, q = q)
end

"""
    addShuntMatpower!(; net, busName, Gs, Bs, in_service=1)

MATPOWER semantics:
- Gs/Bs are shunt admittance parameters given as MW/MVAr at V = 1.0 pu.
- Internally we stamp them as pu-admittance: y_pu = (Gs + j*Bs)/baseMVA.
- IMPORTANT: do NOT add fixed P/Q to the bus power balance.
"""
function addShuntMatpower!(; net::Net, busName::String, Gs::Float64, Bs::Float64, in_service::Int = 1, bus_shunt_model = net.bus_shunt_model)
  shunt_model = normalize_bus_shunt_model(bus_shunt_model)
  busIdx = geNetBusIdx(net = net, busName = busName)
  idShunt = length(net.shuntVec) + 1
  net.shuntDict[busIdx] = idShunt

  vn_kV = getNodeVn(net.nodeVec[busIdx])

  # Create shunt as "PQ-style" just to populate p_shunt/q_shunt fields meaningfully
  # (report values at 1.0 pu). This is NOT used for stamping.
  sh = Shunt(fromBus = busIdx, id = idShunt, base_MVA = net.baseMVA, vn_kV_shunt = vn_kV, p_shunt = Gs, q_shunt = Bs, model = _shunt_component_model(shunt_model), status = in_service)

  # Override stamping admittance to MATPOWER semantics:
  # MATPOWER: Gs/Bs are MW/MVAr at 1 pu -> y_pu = (Gs + j Bs)/baseMVA
  sh.y_pu_shunt = ComplexF64(Gs, Bs) / net.baseMVA
  sh.model = _shunt_component_model(shunt_model)
  sh.G_shunt = real(sh.y_pu_shunt)
  sh.B_shunt = imag(sh.y_pu_shunt)

  push!(net.shuntVec, sh)

  #if in_service == 1
  #    # IMPORTANT: report only (do NOT affect bus injection spec)
  #    addShuntReportPower!(node = net.nodeVec[busIdx], p = Gs, q = Bs)
  #end

  return nothing
end

"""
    setNetBranchStatus!(; net::Net, branchNr::Int, status::Int)

Sets the status of a branch in the network.

# Arguments
- `net::Net`: The network.
- `branchNr::Int`: The number of the branch.
- `status::Int`: The status of the branch. 1 = in service, 0 = out of service.

# Example
```julia
setNetBranchStatus!(net = network, branchNr = 1, status = 1)
```
"""
function setNetBranchStatus!(; net::Net, branchNr::Int, status::Int)
  @assert branchNr > 0 "Branch number must be greater than 0"
  @assert branchNr <= length(net.branchVec) "Branch $branchNr not found in the network"
  # aggregate switch: keeps the per-terminal flags consistent (r0.9.10)
  setBranchStatus!(net.branchVec[branchNr], status == 1)
  markIsolatedBuses!(net = net, log = false)
end

"""
    setBusVoltage!(; net::Net, busName::String, vm_pu::Float64, va_deg::Float64)
Sets the voltage magnitude and angle of a bus in the network.
# Arguments
- `net::Net`: The network.
- `busName::String`: The name of the bus.
- `vm_pu::Float64`: The voltage magnitude in per unit.
- `va_deg::Float64`: The voltage angle in degrees.
# Example 
```julia
setBusVoltage!(net = network, busName = "Bus1", vm_pu = 1.02, va_deg = 5.0) 
```
"""

function setNodeVoltage!(; net::Net, busName::String, vm_pu::Float64, va_deg::Float64)
  @debug "Set bus voltage for $busName to vm_pu = $vm_pu and va_deg = $va_deg"
  busIdx = geNetBusIdx(net = net, busName = busName)
  node = net.nodeVec[busIdx]
  setVmVa!(node = node, vm_pu = vm_pu, va_deg = va_deg)
end

"""
    setBusAngle!(; net::Net, busName::String, va_deg::Float64)
Sets the voltage angle of a bus in the network.
# Arguments
- `net::Net`: The network.
- `busName::String`: The name of the bus.
- `va_deg::Float64`: The voltage angle in degrees.
# Example
```julia
setBusAngle!(net = network, busName = "Bus1", va_deg = 5.0)
```
"""
function setNodeAngle!(; net::Net, busName::String, va_deg::Float64)
  @debug "Set bus angle for $busName to va_deg = $va_deg"
  @assert va_deg >= -360 && va_deg <= 360 "Voltage angle must be between -360 and 360 degrees"
  busIdx = geNetBusIdx(net = net, busName = busName)
  node = net.nodeVec[busIdx]
  vm_pu = node._vm_pu
  setVmVa!(node = node, vm_pu = vm_pu, va_deg = va_deg)
end

"""
    setNetBranchStatus!(; net::Net, branchNr::Int, status::Int)

Sets the status of a branch in the network.

# Arguments
- `net::Net`: The network.
- `branchNr::Int`: The number of the branch.
- `status::Int`: The status of the branch. 1 = in service, 0 = out of service.

# Example
```julia
  brVec = getNetBranchNumberVec(net = net, fromBus = "B1", toBus = "B2")  
  setNetBranchStatus!(net = net, branchNr = brVec[1], status = 0)
```
"""
function getNetBranchNumberVec(; net::Net, fromBus::String, toBus::String)::Vector{Int}
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)
  brNumberVec = Int[]
  for (i, br) in enumerate(net.branchVec)
    if br.fromBus == from && br.toBus == to
      push!(brNumberVec, i)
    end
  end
  return brNumberVec
end

"""
    getNetBranch(; net::Net, fromBus::String, toBus::String)::Union{Branch,Nothing}

Retrieves the first branch found between two specified buses in the network.

# Arguments
- `net::Net`: The network.
- `fromBus::String`: The name of the bus where the branch starts.
- `toBus::String`: The name of the bus where the branch ends.

# Returns
- `Union{Branch,Nothing}`: The branch between the specified buses, or `nothing` if no such branch exists.

# Example
```julia
getNetBranch(net = network, fromBus = "Bus1", toBus = "Bus2")
"""
function getNetBranch(; net::Net, fromBus::String, toBus::String)::Union{Branch,Nothing}
  from = geNetBusIdx(net = net, busName = fromBus)
  to = geNetBusIdx(net = net, busName = toBus)

  for (i, br) in enumerate(net.branchVec)
    if br.fromBus == from && br.toBus == to
      return br
    end
  end
  return nothing
end

"""
Validate the network configuration.

# Arguments
- `net::Net`: The network to be validated.

# Returns
A tuple `(valid::Bool, message::String)` where `valid` is a boolean indicating whether the network is valid, and `message` is a string containing an error message if the network is invalid.
"""
function validate!(; net = Net, log::Bool = false)::Tuple{Bool,String}
  if length(net.nodeVec) == 0
    return false, "No buses defined in the network"
  end
  if length(net.slackVec) == 0
    return false, "No slack bus defined in the network"
  end
  if length(net.slackVec) > 1
    @info "More than one slack bus defined in the network"
  end
  if length(net.branchVec) == 0
    return false, "No branches defined in the network"
  end
  # Check if bus indices in ascending sequence
  sort!(net.nodeVec, by = x -> x.busIdx)
  for (i, key) in enumerate(net.nodeVec)
    if i != key.busIdx
      return false, "Bus index mismatch for bus $(key.busIdx)"
    end
  end
  markIsolatedBuses!(net = net, log = log)
  return true, "Network is valid"
end

"""
Get the voltage magnitude of a specific bus in the network.

# Arguments
- `net::Net`: The network from which to retrieve the voltage magnitude.
- `busName::String`: The name of the bus.

# Returns
The voltage magnitude of the specified bus.
"""
function get_bus_vn_kV(; net::Net, busName::String)
  busIdx = geNetBusIdx(net = net, busName = busName)
  return getNodeVn(net.nodeVec[busIdx])
end

"""
Get the voltage magnitude of a specific bus in the network.

# Arguments
- `net::Net`: The network from which to retrieve the voltage magnitude.
- `busIdx::Int`: The index of the bus.

# Returns
The voltage magnitude of the specified bus.
"""
function get_vn_kV(; net::Net, busIdx::Int)
  return getNodeVn(net.nodeVec[busIdx])
end

"""
Get the type of a specific bus in the network.

# Arguments
- `net::Net`: The network from which to retrieve the bus type.
- `busName::String`: The name of the bus.

# Returns
The type of the specified bus.
"""
function getBusType(; net::Net, busName::String)
  return getEffectiveBusType(net = net, busName = busName)
end

"""
Lock or unlock the network.

# Arguments
- `net::Net`: The network to be locked or unlocked.
- `locked::Bool`: Boolean indicating whether to lock the network.
"""
function lockNet!(; net::Net, locked::Bool)
  net._locked = locked
end

"""
    setTotalBusPower!(; net::Net, p::Float64, q::Float64)

Sets the total active and reactive power at the buses in the network.

# Arguments
- `net::Net`: The network.
- `p::Float64`: The total active power for the network.
- `q::Float64`: The total reactive power for the network.

# Example
```julia
setTotalBusPower!(net = network, p = 100.0, q = 50.0)
```
"""
function setTotalBusPower!(; net::Net, p::Float64, q::Float64)
  push!(net.totalBusPower, (p, q))
end

"""
    getTotalBusPower(; net::Net)::Tuple{Float64, Float64}

Gets the total active and reactive power for the network.

# Arguments
- `net::Net`: The network.

# Returns
- `n::Tuple{Float64, Float64}`: 

# Example
```julia
getTotalBusPower(net = network)
```
"""
function getTotalBusPower(; net::Net)
  if !isempty(net.totalBusPower)
    n = net.totalBusPower[end]
  else
    n = (0.0, 0.0)
  end
end

"""
Set the total losses in the network.

# Arguments
- `net::Net`: The network to which the losses will be added.
- `pLosses::Float64`: Total active power losses.
- `qLosses::Float64`: Total reactive power losses.
"""
function setTotalLosses!(; net::Net, pLosses::Float64, qLosses::Float64)
  push!(net.totalLosses, (pLosses, qLosses))
end

"""
Get the total losses in the network.

# Arguments
- `net::Net`: The network from which to retrieve the losses.

# Returns
A tuple `(pLosses::Float64, qLosses::Float64)` containing the total active and reactive power losses in the network.
"""
function getTotalLosses(; net::Net)
  if !isempty(net.totalLosses)
    n = net.totalLosses[end]
  else
    n = (0.0, 0.0)
  end
end

"""
    markIsolatedBuses!(;net::Net)

Finds and marks isolated buses in the network.

# Arguments
- `net::Net`: The network.

"""
function markIsolatedBuses!(; net::Net, log::Bool = false)
  empty!(net.isoNodes)

  # Create a set to store the buses that appear in the branches in branchVec
  connected_buses = Set{Int}()

  # Iterate through each branch in branchVec and add the corresponding buses
  # to the set. A one-sided open branch (r0.9.10) connects only its CLOSED
  # bus, and only when something can energize that bus: a prosumer there or
  # another closed branch. A dead closed end (no injection, no other
  # connection) stays isolated exactly like before the feature, its
  # charging has no source to draw from.
  partial_closed = Int[]
  for br in net.branchVec
    state = _branch_terminal_state(br)
    if state == :open
      continue
    elseif state == :open_to
      push!(partial_closed, Int(br.fromBus))
    elseif state == :open_from
      push!(partial_closed, Int(br.toBus))
    else
      push!(connected_buses, br.fromBus)
      push!(connected_buses, br.toBus)
    end
  end
  if !isempty(partial_closed)
    # only a reference can balance an island on its own; a machine or load
    # behind a dead partial branch stays de-energized (pre-feature behavior)
    reference_buses = Set{Int}(getPosumerBusIndex(ps) for ps in net.prosumpsVec if ps.referencePri !== nothing)
    for bus in partial_closed
      (bus in connected_buses || bus in reference_buses) && push!(connected_buses, bus)
    end
  end

  # Iterate through all buses in nodeVec and mark the isolated buses
  for bus in net.nodeVec
    # Check whether the bus number is not included in the set of connected buses
    if !(bus.busIdx in connected_buses)
      setNodeType!(bus, "Isolated")
      setVmVa!(node = bus, vm_pu = 0.0, va_deg = 0.0)
      push!(net.isoNodes, bus.busIdx)
      if log
        @info "Bus $(getCompName(bus.comp)) is isolated"
      end
    end
  end

  sort!(net.isoNodes)
end

function _buildQLimits!(net::Net; reset::Bool = true)
  # Number of buses from nodeVec
  nbus = length(net.nodeVec)

  # Ensure arrays have correct length
  if length(net.qmin_pu) != nbus
    resize!(net.qmin_pu, nbus)
  end
  if length(net.qmax_pu) != nbus
    resize!(net.qmax_pu, nbus)
  end

  # Reinitialize limits on every call (derived data -> safe to overwrite)
  # qmin_pu starts at +Inf (will be replaced by first finite sum)
  # qmax_pu starts at -Inf (will be replaced by first finite sum)
  fill!(net.qmin_pu, Inf)
  fill!(net.qmax_pu, -Inf)

  # Optionally reset Q-limit log
  if reset
    resetQLimitLog!(net)
  end

  # Aggregate generator Q-limits per bus
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    bus = getPosumerBusIndex(ps)
    # Safety: allow for buses beyond nodeVec if that can happen in your data
    if bus > nbus
      # grow arrays if needed
      resize!(net.qmin_pu, bus)
      resize!(net.qmax_pu, bus)
      for b = (nbus+1):bus
        net.qmin_pu[b] = Inf
        net.qmax_pu[b] = -Inf
      end
      nbus = bus
    end
    # Convert to p.u., handle 'no limit' as ±Inf
    qmin_pu = isnothing(ps.minQ) ? -Inf : ps.minQ / net.baseMVA
    qmax_pu = isnothing(ps.maxQ) ? Inf : ps.maxQ / net.baseMVA
    # Sum instead of min/max; respect sentinel values
    cur_qmin = net.qmin_pu[bus]
    net.qmin_pu[bus] = isfinite(cur_qmin) ? (cur_qmin + qmin_pu) : qmin_pu
    cur_qmax = net.qmax_pu[bus]
    net.qmax_pu[bus] = isfinite(cur_qmax) ? (cur_qmax + qmax_pu) : qmax_pu
  end

  return nothing
end

"""
    buildQLimits!(net; reset=true)

Public wrapper to (re)build aggregated per-bus Q limits in p.u.
"""
function buildQLimits!(net::Net; reset::Bool = true)
  return _buildQLimits!(net; reset = reset)
end

"""
    setQLimits!(; net::Net,
                 qmin_MVar::Float64,
                 qmax_MVar::Float64,
                 busName::Union{Nothing,String}=nothing)

Sets the reactive power (Q) limits of generator prosumers and then rebuilds
the aggregated bus-level Q-limits (`qmin_pu`, `qmax_pu`).

- Without `busName`: all generators receive the specified `qmin_MVar`/`qmax_MVar`.
- With `busName`: only generators connected to the specified bus receive the new limits.
"""

function setQLimits!(; net::Net, qmin_MVar::Float64, qmax_MVar::Float64, busName::Union{Nothing,AbstractString,AbstractVector{<:AbstractString}} = nothing)

  # Normalize busName to a set of bus indices or `nothing`
  busIdxSet = nothing
  if !isnothing(busName)
    # Ensure we always iterate over a collection of names
    names = busName isa AbstractString ? (busName,) : busName
    busIdxSet = Set(geNetBusIdx(net = net, busName = bn) for bn in names)
  end

  for ps in net.prosumpsVec
    #TODO is load also possible?
    isGenerator(ps) || continue

    if isnothing(busIdxSet)
      # No bus filter: apply to all generators
      ps.minQ = qmin_MVar
      ps.maxQ = qmax_MVar
    else
      # Apply only to generators connected to the selected buses
      c = ps.comp
      if !isnothing(c.cFrom_bus) && (c.cFrom_bus in busIdxSet)
        ps.minQ = qmin_MVar
        ps.maxQ = qmax_MVar
      end
    end
  end

  _buildQLimits!(net; reset = true)

  return nothing
end

"""
    setPVBusVset!(net::Net, busName::String; vm_pu::Float64)
  Sets the voltage magnitude setpoint for a PV bus in the network.
  Only applicable for buses of type "PV" and for testing purpose.
"""

function setPVBusVset!(; net::Net, busName::String, vm_pu::Float64)
  bus = geNetBusIdx(net = net, busName = busName)
  net.nodeVec[bus]._vm_pu = vm_pu
  for ps in net.prosumpsVec
    getPosumerBusIndex(ps) == bus || continue
    (isSlack(ps) || (isGenerator(ps) && isRegulating(ps))) || continue
    ps.vm_pu = vm_pu
  end
end

function _wattterfillQ(Q_target::Float64, q_spec::Vector{Float64}, qmin::Vector{Float64}, qmax::Vector{Float64}; tol::Float64 = 1e-6, maxiter::Int = 50)
  n = length(q_spec)
  @assert length(qmin) == n && length(qmax) == n

  # Start: clamp specs in their box
  q = similar(q_spec)
  for i = 1:n
    q[i] = clamp(q_spec[i], qmin[i], qmax[i])
  end

  for it = 1:maxiter
    Δ = Q_target - sum(q)
    if abs(Δ) <= tol
      return q
    end

    # Frei je nach Richtung (nach oben oder unten noch Bewegungsfreiheit)
    free = if Δ > 0
      [i for i = 1:n if q[i] < qmax[i] - tol]
    else
      [i for i = 1:n if q[i] > qmin[i] + tol]
    end

    isempty(free) && return q  # nichts mehr zu holen, Ziel nicht exakt erreichbar

    step = Δ / length(free)
    for i in free
      q[i] = clamp(q[i] + step, qmin[i], qmax[i])
    end
  end

  return q
end

function _generators_at_bus(net::Net, bus::Int)
  idx = Int[]
  for (k, ps) in enumerate(net.prosumpsVec)
    if isGenerator(ps) && getPosumerBusIndex(ps) == bus
      push!(idx, k)
    end
  end
  return idx
end

function _loads_at_bus(net::Net, bus::Int)
  idx = Int[]
  for (k, ps) in enumerate(net.prosumpsVec)
    if !isGenerator(ps) && getPosumerBusIndex(ps) == bus
      push!(idx, k)
    end
  end
  return idx
end

function distribute_bus_generation!(net::Net, bus::Int)
  node = net.nodeVec[bus]

  gens_idx = _generators_at_bus(net, bus)
  isempty(gens_idx) && return

  # Zielwerte am Bus (hier: aus Node – kannst du anpassen)
  P_target = node._pƩGen === nothing ? 0.0 : node._pƩGen
  Q_target = node._qƩGen === nothing ? 0.0 : node._qƩGen

  # Eingabewerte + Limits pro Generator
  P_spec = Float64[]
  Q_spec = Float64[]
  qmin   = Float64[]
  qmax   = Float64[]

  for idx in gens_idx
    ps = net.prosumpsVec[idx]
    push!(P_spec, isnothing(ps.pVal) ? 0.0 : ps.pVal)
    push!(Q_spec, isnothing(ps.qVal) ? 0.0 : ps.qVal)
    push!(qmin, isnothing(ps.minQ) ? -Inf : ps.minQ)
    push!(qmax, isnothing(ps.maxQ) ? Inf : ps.maxQ)
  end

  # --- P-Verteilung (ohne Limits, meist einfacher) ---
  P_alloc = similar(P_spec)
  P_sum   = sum(P_spec)

  if abs(P_sum) < 1e-9
    # gleichmäßig, wenn keine sinnvollen spec-Werte
    fill!(P_alloc, P_target / length(P_alloc))
  else
    # proportional zu den specs
    P_alloc .= P_target .* (P_spec ./ P_sum)
  end

  # --- Q-Verteilung mit Waterfilling ---
  Q_alloc = _wattterfillQ(Q_target, Q_spec, qmin, qmax)

  # Ergebnisse zurück in Prosumer
  for (pq, idx) in enumerate(gens_idx)
    ps = net.prosumpsVec[idx]
    setPQResult!(ps, P_alloc[pq], Q_alloc[pq])
  end
end

function distribute_bus_loads!(net::Net, bus::Int)
  node = net.nodeVec[bus]
  loads_idx = _loads_at_bus(net, bus)
  isempty(loads_idx) && return

  P_target = node._pƩLoad === nothing ? 0.0 : node._pƩLoad
  Q_target = node._qƩLoad === nothing ? 0.0 : node._qƩLoad

  P_spec = [isnothing(ps.pVal) ? 0.0 : ps.pVal for ps in net.prosumpsVec[loads_idx]]
  Q_spec = [isnothing(ps.qVal) ? 0.0 : ps.qVal for ps in net.prosumpsVec[loads_idx]]

  P_alloc = similar(P_spec)
  Q_alloc = similar(Q_spec)

  # einfache Variante: spec = result
  if abs(sum(P_spec)) < 1e-9
    P_alloc .= P_target / length(P_alloc)
  else
    P_alloc .= P_target .* (P_spec ./ sum(P_spec))
  end

  if abs(sum(Q_spec)) < 1e-9
    Q_alloc .= Q_target / length(Q_alloc)
  else
    Q_alloc .= Q_target .* (Q_spec ./ sum(Q_spec))
  end

  for (pq, idx) in enumerate(loads_idx)
    ps = net.prosumpsVec[idx]
    setPQResult!(ps, -P_alloc[pq], -Q_alloc[pq])  # Vorzeichen je nach Konvention
  end
end

function distributeBusResults!(net::Net)
  # reset results
  for ps in net.prosumpsVec
    ps.pRes = nothing
    ps.qRes = nothing
  end

  for node in net.nodeVec
    distribute_bus_generation!(net, node.busIdx)
    distribute_bus_loads!(net, node.busIdx)
  end
end

"""
    buildVoltageVector(net::Net) -> Vector{ComplexF64}

Builds the complex bus voltage vector V[k] = vm_pu[k] * exp(j * va_rad[k])
using the current nodal state stored in `net.nodeVec`.
"""
function buildVoltageVector(net::Net)
  V = Vector{ComplexF64}(undef, length(net.nodeVec))
  @inbounds for n in net.nodeVec
    V[n.busIdx] = n._vm_pu * cis(deg2rad(n._va_deg))
  end
  return V
end

# Case-specific diagnosis for "no solvable slack": the generic message hid
# the two real situations — a selected slack sitting on an isolated bus
# (single-bus/branch-less deliveries such as one side of a DC border), and a
# net where no slack was ever registered.
function _no_slack_message(where_::String, net::Net)::String
  n = length(net.nodeVec)
  iso = length(net.isoNodes)
  if n > 0 && iso == n
    return "$(where_): all $(n) bus(es) of this network are isolated (no in-service branch) — there is no energized AC system to solve. A branch-less delivery (e.g. a single-bus area or one side of a DC border crossing) is only solvable in its assembled model."
  end
  if !isempty(net.slackVec)
    slack_names = join((getCompName(net.nodeVec[i].comp) for i in net.slackVec if 1 <= i <= n), ", ")
    return "$(where_): the selected slack ($(slack_names)) lies on an isolated bus ($(iso) of $(n) buses isolated) — no energized bus carries the reference. Check the isolation report (ac_islands.csv) or choose a slack in the energized part."
  end
  return "$(where_): no slack bus registered — the network needs one voltage reference (referencePriority, an ExternalNetworkInjection, or a generator marked as slack). Enable power_flow.auto_slack to let the solver pick the strongest candidate itself."
end

# Rank a slack candidate by its dispatchable size. Rated apparent power is the
# most honest capacity figure; active-power capability and the current dispatch
# are the fallbacks when a data source does not carry a rating.
function _slack_candidate_strength(ps::ProSumer)::Float64
  ps.ratedS !== nothing && ps.ratedS > 0.0 && return Float64(ps.ratedS)
  ps.maxP !== nothing && ps.maxP > 0.0 && return Float64(ps.maxP)
  ps.pVal !== nothing && return abs(Float64(ps.pVal))
  return 0.0
end

"""
    ensureSlack!(net::Net; log::Bool = true) -> Union{Nothing,Int}

Make sure the network has a usable voltage reference. When at least one slack
is already registered on a non-isolated bus, nothing changes and `nothing` is
returned. Otherwise the strongest injection candidate on a non-isolated bus is
promoted to slack (`referencePri = 1`) and its bus index is returned.

Candidate ranking: `ExternalNetworkInjection` units win over generators and
synchronous machines (an external grid equivalent is the natural reference);
within each group the largest unit wins, sized by `ratedS`, then `maxP`, then
the current dispatch `|p|`. Static var compensators are never promoted — they
carry no active power. When no candidate exists the network is left untouched
and the regular no-slack error will name the situation.

The promotion mutates the network: the chosen prosumer becomes the reference
and the bus voltage gets a `1.0 pu / 0.0°` setpoint if it has none. Enabled at
solve time via the `power_flow.auto_slack` configuration key (default `false`)
or the `auto_slack` keyword of [`runpf!`](@ref).
"""
function ensureSlack!(net::Net; log::Bool = true)::Union{Nothing,Int}
  nbus = length(net.nodeVec)
  usable(bus) = (1 <= bus <= nbus) && !(bus in net.isoNodes)
  for ps in net.prosumpsVec
    isSlack(ps) && usable(getPosumerBusIndex(ps)) && return nothing
  end
  best = nothing
  best_key = (-1, -Inf)
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    # SVCs are injections in the prosumption typing but cannot carry the
    # island's active-power balance.
    ps.comp.cTyp == StaticVarCompensator && continue
    bus = getPosumerBusIndex(ps)
    usable(bus) || continue
    key = (ps.comp.cTyp == ExternalNetworkInjection ? 1 : 0, _slack_candidate_strength(ps))
    if key > best_key
      best_key = key
      best = ps
    end
  end
  best === nothing && return nothing
  best.referencePri = 1
  bus = getPosumerBusIndex(best)
  node = net.nodeVec[bus]
  if node._vm_pu === nothing || node._vm_pu <= 0.0
    vm = (best.vm_pu !== nothing && best.vm_pu > 0.0) ? Float64(best.vm_pu) : 1.0
    setVmVa!(node = node, vm_pu = vm, va_deg = 0.0)
  end
  refreshBusTypesFromProsumers!(net)
  if log
    kind = best.comp.cTyp == ExternalNetworkInjection ? "external network injection" : "generator"
    strength = _slack_candidate_strength(best)
    println("auto_slack: no usable slack registered — promoted $(kind) '$(getCompName(best.comp))' at bus '$(getCompName(node.comp))' ($(round(strength; digits = 1)) MVA/MW) to slack.")
  end
  return bus
end

"""
    initial_Vrect_from_net(net) -> (V0, slack_idx)

Build the initial complex voltage vector V0 from the network bus data
(Vm, Va), and detect the slack bus index.

Returns:
- V0::Vector{ComplexF64}
- slack_idx::Int
"""
function initialVrect(net::Net; flatstart::Bool = net.flatstart)
  nodes = net.nodeVec
  n = length(nodes)

  slack_idx = findfirst(n -> getNodeType(n) == Slack, nodes)
  slack_idx === nothing && error(_no_slack_message("initialVrect", net))

  V0 = Vector{ComplexF64}(undef, n)

  @inbounds for k = 1:n
    node = nodes[k]

    # Voltage magnitude guess
    vm = (node._vm_pu === nothing || node._vm_pu <= 0.0) ? 1.0 : Float64(node._vm_pu)

    if flatstart
      if k == slack_idx
        # The slack/reference bus remains fixed during the solve. Even for a
        # flat start, keep its imported or explicitly configured reference
        # angle together with its voltage magnitude.
        va_deg = (node._va_deg === nothing) ? 0.0 : Float64(node._va_deg)
        va = deg2rad(va_deg)
        V0[k] = ComplexF64(vm * cos(va), vm * sin(va))
      elseif getNodeType(node) == PV
        # Flat start resets non-reference angles, but keeps voltage magnitude
        # setpoints for PV buses. PQ-bus magnitudes remain the flat 1.0 pu guess.
        V0[k] = ComplexF64(vm, 0.0)
      else
        V0[k] = ComplexF64(1.0, 0.0)
      end
    else
      # Use stored angle if available, else 0
      va_deg = (node._va_deg === nothing) ? 0.0 : Float64(node._va_deg)
      va = deg2rad(va_deg)
      V0[k] = ComplexF64(vm * cos(va), vm * sin(va))
    end
  end

  return V0, slack_idx
end

"""
    buildComplexSVec(net) -> S::Vector{ComplexF64}

Build the specified complex power injection vector S = P + jQ in per-unit
for each bus, based on the net's bus load / generation / shunt data.
Positive P/Q means net injection into the bus (generation),
negative means net consumption (load).

"""
function buildComplexSVec(net::Net)
  n = length(net.nodeVec)
  Pgen = zeros(Float64, n)
  Qgen = zeros(Float64, n)
  Pload = zeros(Float64, n)
  Qload = zeros(Float64, n)
  baseMVA = net.baseMVA

  @inbounds for ps in net.prosumpsVec
    bus = getPosumerBusIndex(ps)
    (1 <= bus <= n) || continue
    p = isnothing(ps.pVal) ? 0.0 : ps.pVal
    q = isnothing(ps.qVal) ? 0.0 : ps.qVal

    if isGenerator(ps)
      Pgen[bus] += p
      Qgen[bus] += q
    else
      Pload[bus] += p
      Qload[bus] += q
    end
  end

  S = Vector{ComplexF64}(undef, n)
  @inbounds for bus in eachindex(S)
    Pinj = (Pgen[bus] - Pload[bus]) / baseMVA
    Qinj = (Qgen[bus] - Qload[bus]) / baseMVA
    S[bus] = ComplexF64(Pinj, Qinj)
  end

  # Voltage-dependent bus shunts are represented as specified injections so
  # they are not double-counted in Ybus. Their equivalent complex-power term is
  # S_sh = |V|^2 * conj(y_sh), so the net specified injection gets the negative
  # of that contribution. Without a voltage argument this base vector uses
  # |V| = 1 pu; the rectangular solver updates it per iteration.
  for sh in net.shuntVec
    sh.status == 0 && continue
    sh.model == :VoltageDependentInjection || continue
    (sh.busIdx in net.isoNodes) && continue
    S[sh.busIdx] -= conj(sh.y_pu_shunt)
  end

  return S
end

function has_voltage_dependent_control(net::Net)::Bool
  for sh in net.shuntVec
    if sh.status != 0 && sh.model == :VoltageDependentInjection
      return true
    end
  end
  for ps in net.prosumpsVec
    if has_qu_controller(ps) || has_pu_controller(ps)
      return true
    end
  end
  return false
end

@inline _safe_vm(vm::Float64; eps::Float64 = 1e-9) = vm > eps ? vm : eps

"""
    buildControlledSVec(net, V) -> (Sspec, dPinj_dVm, dQinj_dVm)

Evaluate bus injections for the current voltage state. If a prosumer carries a
`PUController` and/or `QUController`, the corresponding active/reactive setpoint
is evaluated as a function of local bus voltage magnitude. The derivative
vectors contain `dP_spec/d|V|` and `dQ_spec/d|V|` (both in p.u. per p.u.).
"""
function buildControlledSVec(net::Net, V::Vector{ComplexF64})
  n = length(net.nodeVec)
  length(V) == n || error("buildControlledSVec: voltage vector length mismatch.")

  Sspec = zeros(ComplexF64, n)
  dPinj_dVm = zeros(Float64, n)
  dQinj_dVm = zeros(Float64, n)
  base = net.baseMVA

  for ps in net.prosumpsVec
    bus = getPosumerBusIndex(ps)
    (1 <= bus <= n) || continue

    vm = abs(V[bus])
    vm_safe = _safe_vm(vm)

    p_mw = isnothing(ps.pVal) ? 0.0 : ps.pVal
    q_mvar = isnothing(ps.qVal) ? 0.0 : ps.qVal
    dp_dvm_mw = 0.0
    dq_dvm_mvar = 0.0

    if has_pu_controller(ps)
      p_pu, dp_dvm_pu = evaluate_controller(ps.puController, vm_safe)
      p_mw = p_pu * base
      dp_dvm_mw = dp_dvm_pu * base
    end

    if has_qu_controller(ps)
      q_pu, dq_dvm_pu = evaluate_controller(ps.quController, vm_safe)
      q_mvar = q_pu * base
      dq_dvm_mvar = dq_dvm_pu * base
    end

    sign = isGenerator(ps) ? 1.0 : -1.0
    Sspec[bus] += ComplexF64(sign * p_mw / base, sign * q_mvar / base)
    dPinj_dVm[bus] += sign * dp_dvm_mw / base
    dQinj_dVm[bus] += sign * dq_dvm_mvar / base
  end

  for sh in net.shuntVec
    sh.status == 0 && continue
    sh.model == :VoltageDependentInjection || continue
    bus = sh.busIdx
    (1 <= bus <= n) || continue
    (bus in net.isoNodes) && continue
    vm = abs(V[bus])
    vm_safe = _safe_vm(vm)
    # The shunt contribution follows the same sign as admittance stamping:
    # S_sh,p.u. = |V|^2 * conj(y_sh). The mismatch uses
    # F = S_calc - S_spec, so this contribution is subtracted from S_spec.
    yconj = conj(sh.y_pu_shunt)
    Sspec[bus] -= vm^2 * yconj
    dS_dVm = -2.0 * vm_safe * yconj
    dPinj_dVm[bus] += real(dS_dVm)
    dQinj_dVm[bus] += imag(dS_dVm)
  end

  return Sspec, dPinj_dVm, dQinj_dVm
end

"""
    updateShuntPowers!(net; reset_node=true)

Recomputes shunt P/Q (MW/MVar) from solved bus voltages and shunt pu-admittances.
Writes back:
- sh.p_shunt / sh.q_shunt (results)
- node._pShunt / node._qShunt (for reporting)
"""
function updateShuntPowers!(; net::Net, reset_node::Bool = true)
  # Optionally reset node aggregates (recommended)
  if reset_node
    for n in net.nodeVec
      n._pShunt = 0.0
      n._qShunt = 0.0
    end
  end

  base = net.baseMVA

  for sh in net.shuntVec
    sh.status == 0 && continue
    sh.model in (:Y, :VoltageDependentInjection) || continue

    bus = sh.busIdx
    # ignore isolated buses if you want:
    (bus in net.isoNodes) && continue

    vm = net.nodeVec[bus]._vm_pu
    vm2 = vm * vm

    # S_pu = |V|^2 * conj(y_pu)
    S_pu = vm2 * conj(sh.y_pu_shunt)
    P = real(S_pu) * base
    Q = imag(S_pu) * base

    sh.p_shunt = P
    sh.q_shunt = Q

    # add to node report aggregate
    net.nodeVec[bus]._pShunt += P
    net.nodeVec[bus]._qShunt += Q
  end
end

# Copyright 2023-2026 Udo Schmitz
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
# file: src/adapters/cgmes/cgmes_measurements.jl
# purpose: the CGMES side of the measurement table: the key index that maps
#          terminal, node and Analog mRIDs of an imported delivery to network
#          locations, and the synthetic measurement generator that turns the
#          SV profile of a delivery into a measurement table, and the CGMES
#          method of the pseudo-row generator.

# equipment the mapping turns into prosumers (one terminal, bus injection);
# an EnergySource only with cgmes_import.map_energy_source and a
# PowerElectronicsConnection only with cgmes_import.map_power_electronics,
# otherwise it is in skipped_equipment and resolves to :not_modelled before
# the class test
const _MEASUREMENT_INJECTION_CLASSES = (:EnergyConsumer, :ConformLoad, :NonConformLoad, :StationSupply, :SynchronousMachine, :AsynchronousMachine, :ExternalNetworkInjection, :EquivalentInjection, :StaticVarCompensator, :EnergySource, :PowerElectronicsConnection)
# equipment the mapping turns into bus shunts
const _MEASUREMENT_SHUNT_CLASSES = (:LinearShuntCompensator, :NonlinearShuntCompensator)

# Analog.measurementType values the table carries, and their table type
const _MEASUREMENT_ANALOG_TYPES = Dict(
  "ThreePhaseActivePower" => "P",
  "ThreePhaseReactivePower" => "Q",
  "LineCurrent" => "I",
  "LineToLineVoltage" => "U",
)

"""
    cgmesMeasurementKeys(result::CGMESImportResult) -> MeasurementKeyIndex

The key index of a CGMES import for `importMeasurementTable`: every
terminal, TopologicalNode, ConnectivityNode, BusbarSection and Analog of
the delivery, resolved to the bus, branch end or bus group of `result.net`.

- A terminal of a line or transformer end resolves to that branch end. A
  branch out of service or open at one terminal resolves to
  `:not_modelled` at both ends (the estimator does not take flow
  measurements on a partially open branch yet).
- A terminal of a load, machine, external grid, equivalent injection,
  SVC, power electronics connection or mapped energy source resolves to
  `:injection`, a terminal of a shunt compensator to
  `:shunt`, as long as the equipment is mapped, its terminal connected and
  its bus energized; each such terminal is a member of the group sum at
  its bus.
- Every other terminal (switches, equipment the importer skips or does not
  map, de-energized nodes) resolves to `:not_modelled` with the reason.
- A node resolves to its bus, a ConnectivityNode through its
  TopologicalNode, a BusbarSection through its terminal. A node that is
  not energized resolves to bus 0. A bus that only closed switches join to
  the network (no branch of its own, typed `Isolated` by the import) is
  energized: the solver fuses it with its switch group.
- A bus or branch end without any mRID (the star point of a three-winding
  transformer, the legs ending there, auxiliary buses of the import) takes
  its name key from `nameMeasurementKeys`, in canonical form; a name never
  replaces an mRID key.
- An Analog resolves to its terminal (or, without one, to its
  PowerSystemResource when that is a node or busbar) and its
  `measurementType`: active power `P`, reactive power `Q`, line current
  `I`, line-to-line voltage `U`; other types keep an empty type and are
  reported when a row uses them.
"""
function cgmesMeasurementKeys(result::CGMESImportResult)::Sparlectra.MeasurementKeyIndex
  net = result.net
  store = result.store
  topo = result.topo
  key = Sparlectra._measurement_key
  iso = Sparlectra._unfed_isolated_buses(net)
  pos = Dict(br.branchIdx => k for (k, br) in enumerate(net.branchVec))
  bus_of_tn(tn) = (tn === nothing || !haskey(topo.bus_name, tn)) ? 0 : get(net.busDict, topo.bus_name[tn], 0)
  energized(b) = b > 0 && !(b in iso)
  # node-breaker deliveries often give a terminal only its ConnectivityNode
  # (switch terminals in the ENTSO-E SmallGrid); the TopologicalNode then
  # comes from ConnectivityNode.TopologicalNode
  function tn_of(t::CIMObject)
    tn = tnOfTerminal(t)
    tn === nothing || return tn
    cn = get(t.refs, :ConnectivityNode, nothing)
    cn === nothing && return nothing
    cno = get(store.objects, cn, nothing)
    return cno === nothing ? nothing : get(cno.refs, :TopologicalNode, nothing)
  end

  nodes = Dict{String,Int}()
  for tn in objectsOf(store, :TopologicalNode)
    b = bus_of_tn(tn.mrid)
    nodes[key(tn.mrid)] = energized(b) ? b : 0
  end
  for cn in objectsOf(store, :ConnectivityNode)
    tn = get(cn.refs, :TopologicalNode, nothing)
    tn === nothing && continue
    b = bus_of_tn(tn)
    nodes[key(cn.mrid)] = energized(b) ? b : 0
  end
  for bbs in objectsOf(store, :BusbarSection)
    ts = get(topo.terminals, bbs.mrid, nothing)
    (ts === nothing || isempty(ts)) && continue
    b = bus_of_tn(tn_of(ts[1]))
    nodes[key(bbs.mrid)] = energized(b) ? b : 0
  end

  terminals = Dict{String,Sparlectra.MeasurementTarget}()
  groups = Dict{Tuple{Symbol,Int},Vector{String}}()
  for t in objectsOf(store, :Terminal)
    tk = key(t.mrid)
    b = bus_of_tn(tn_of(t))
    eq_mrid = get(t.refs, :ConductingEquipment, nothing)
    eq = eq_mrid === nothing ? nothing : get(store.objects, eq_mrid, nothing)
    cls = eq === nothing ? :unknown : eq.class
    side = get(result.branch_side_of_terminal, t.mrid, nothing)
    target = if side !== nothing
      bidx, dir = side
      k = get(pos, bidx, 0)
      # the terminal state, not the aggregate status alone: a branch open at
      # one terminal is left out at both ends, since runse! refuses flow
      # measurements on it and would fail the whole estimate
      state = k == 0 ? :open : Sparlectra._branch_terminal_state(net.branchVec[k])
      if state === :open
        Sparlectra.MeasurementTarget(:not_modelled, b, 0, :none, "$(cls) out of service")
      elseif state !== :closed
        Sparlectra.MeasurementTarget(:not_modelled, b, 0, :none, "$(cls) open at one terminal; flow on partially open branches is not supported by the estimator yet")
      else
        br = net.branchVec[k]
        Sparlectra.MeasurementTarget(:branch, Int(dir === :from ? br.fromBus : br.toBus), k, dir, string(cls))
      end
    elseif !energized(b)
      Sparlectra.MeasurementTarget(:not_modelled, b, 0, :none, "$(cls) at a node that is not energized")
    elseif eq_mrid === nothing || eq_mrid in result.skipped_equipment
      Sparlectra.MeasurementTarget(:not_modelled, b, 0, :none, "$(cls) not mapped by the import")
    elseif !terminalConnected(t)
      Sparlectra.MeasurementTarget(:not_modelled, b, 0, :none, "$(cls) terminal disconnected")
    elseif cls in _MEASUREMENT_INJECTION_CLASSES
      push!(get!(Vector{String}, groups, (:injection, b)), tk)
      Sparlectra.MeasurementTarget(:injection, b, 0, :none, string(cls))
    elseif cls in _MEASUREMENT_SHUNT_CLASSES
      push!(get!(Vector{String}, groups, (:shunt, b)), tk)
      Sparlectra.MeasurementTarget(:shunt, b, 0, :none, string(cls))
    elseif (lk = get(result.element_by_mrid, eq_mrid, nothing)) !== nothing && lk[1] === :link
      # a closed switch between two TopologicalNodes is a bus link: P and Q
      # at its terminals are link flow measurements (weighted WLS rows of
      # the estimator). A switch inside one TopologicalNode is no link and
      # stays not modelled below
      l = net.linkVec[lk[2]]
      dir = b == l.fromBus ? :from : (b == l.toBus ? :to : :none)
      dir === :none ? Sparlectra.MeasurementTarget(:not_modelled, b, 0, :none, "$(cls) terminal is at neither bus of its link") : Sparlectra.MeasurementTarget(:link, b, lk[2], dir, "$(cls) bus link")
    else
      Sparlectra.MeasurementTarget(:not_modelled, b, 0, :none, "$(cls) is not a modelled equipment")
    end
    terminals[tk] = target
  end

  analogs = Dict{String,Tuple{String,String,String}}()
  for a in objectsOf(store, :Analog)
    mtype = something(str(a, :measurementType), "")
    typ = get(_MEASUREMENT_ANALOG_TYPES, mtype, "")
    detail = isempty(typ) ? "measurementType '$(mtype)' is not carried by the table" : mtype
    term = get(a.refs, :Terminal, nothing)
    psr = get(a.refs, :PowerSystemResource, nothing)
    loc = term !== nothing ? key(term) : (psr !== nothing && haskey(nodes, key(psr)) ? key(psr) : "")
    if isempty(loc)
      typ = ""
      detail = "Analog has neither a terminal nor a node as location"
    end
    analogs[key(a.mrid)] = (loc, typ, detail)
  end

  # Locations the delivery has no mRID for: the star point of a
  # three-winding transformer and the legs that end there, auxiliary buses
  # the import adds. They take the name key of nameMeasurementKeys (in
  # canonical form), so a measurement there has an address; a name never
  # replaces an mRID key.
  named = Sparlectra.nameMeasurementKeys(net)
  covered_bus = Set{Int}(b for b in values(nodes) if b > 0)
  for (k, b) in named.nodes
    (b > 0 && !(b in covered_bus)) || continue
    ck = key(k)
    haskey(nodes, ck) || haskey(terminals, ck) || (nodes[ck] = b)
  end
  covered_end = Set{Tuple{Int,Symbol}}((t.branchIdx, t.direction) for t in values(terminals) if t.kind === :branch)
  for (k, t) in named.terminals
    (t.kind === :branch && !((t.branchIdx, t.direction) in covered_end)) || continue
    ck = key(k)
    haskey(nodes, ck) || haskey(terminals, ck) || (terminals[ck] = t)
  end

  bus_kV = Float64[nd.comp.cVN for nd in net.nodeVec]
  return Sparlectra.MeasurementKeyIndex(terminals, nodes, analogs, groups, bus_kV)
end

"""
    generateMeasurementTableFromSV(result::CGMESImportResult; kwargs...) -> Vector{NamedTuple}

Synthetic measurement table from the SV profile of a CGMES delivery: U at
every energized TopologicalNode with an `SvVoltage`, and P and Q at both
terminals of every line and transformer end closed at both terminals, computed with the
network model at the SV voltage state (a delivery often carries
`SvPowerFlow` only at injections). Rows are keyed by `terminal_mrid` and
`node_mrid` and sorted by key, so the output does not depend on the order
of the delivery.

Keywords:
- `timestamp = ""`: the timestamp column of every row.
- `include_injections = false`: also write the delivered `SvPowerFlow` of
  injection terminals (P and Q) and shunt terminals (Q).
- `sigma_pq_rel = 0.01`, `sigma_pq_floor = 0.1`: sigma of P and Q as
  `max(sigma_pq_rel * |value|, sigma_pq_floor)` in MW/MVAr.
- `sigma_u_rel = 0.002`: sigma of U relative to the reading (kV).
  `sigma_pq_floor` and `sigma_u_rel` must be positive, `sigma_pq_rel` and
  `bad_data_sigmas` finite and at least 0 (an `ArgumentError` otherwise),
  so every row passes the reader.
- `sigma_by_voltage_level` (default: the active
  `state_estimation.sigma_by_voltage_level`): where the level of the
  measuring bus (the node, the branch end, the injection or shunt) has a
  row with a column for the type, that sigma
  (`max(absolute, rel_pct/100 * |value|)`, MW, MVAr, kV) replaces the
  three keywords above, for the written sigma and the noise alike.
- `noise = true`: add Gaussian noise with that sigma; `false` writes the SV
  values unchanged.
- `bad_data_fraction = 0.0`: share of rows that get a gross error of
  `bad_data_sigmas` sigma (random sign) instead of the noise. They carry
  `source = "synthetic;bad_data"`, all others `source = "synthetic"`.
- `bad_data_sigmas = 20.0`, `seed = 1`: size of the gross error and the
  seed of the random stream; equal arguments give equal rows.

The rows have the fields of the table header and can be written with
`writeMeasurementTable`.
"""
function generateMeasurementTableFromSV(result::CGMESImportResult; timestamp::AbstractString = "", include_injections::Bool = false, sigma_pq_rel::Real = 0.01, sigma_pq_floor::Real = 0.1, sigma_u_rel::Real = 0.002, sigma_by_voltage_level = Sparlectra.state_estimation_config().sigma_by_voltage_level, noise::Bool = true, bad_data_fraction::Real = 0.0, bad_data_sigmas::Real = 20.0, seed::Integer = 1)
  0.0 <= bad_data_fraction <= 1.0 || throw(ArgumentError("generateMeasurementTableFromSV: bad_data_fraction must be in [0, 1], got $(bad_data_fraction)"))
  Sparlectra._check_table_sigmas("generateMeasurementTableFromSV", sigma_pq_rel, sigma_pq_floor, sigma_u_rel, bad_data_sigmas)
  levels = Sparlectra._sigma_levels(sigma_by_voltage_level)
  net = result.net
  store = result.store
  topo = result.topo
  key = Sparlectra._measurement_key
  iso = Sparlectra._unfed_isolated_buses(net)
  pos = Dict(br.branchIdx => k for (k, br) in enumerate(net.branchVec))

  # the SV voltage state per bus; a bus without SvVoltage keeps the start
  # value of the import
  V = ComplexF64[something(nd._vm_pu, 1.0) * cis(deg2rad(something(nd._va_deg, 0.0))) for nd in net.nodeVec]
  # truth values: (key kind, key, type, value, measuring bus)
  truth = Tuple{Symbol,String,String,Float64,Int}[]
  for sv in objectsOf(store, :SvVoltage)
    tn = get(sv.refs, :TopologicalNode, nothing)
    (tn === nothing || !haskey(topo.bus_name, tn)) && continue
    b = get(net.busDict, topo.bus_name[tn], 0)
    (b == 0 || b in iso) && continue
    v = num(sv, :v)
    (v === nothing || v <= 0.0) && continue
    V[b] = (v / topo.vn_kV[tn]) * cis(deg2rad(num(sv, :angle, 0.0)))
    push!(truth, (:node, key(tn), "U", v, b))
  end
  # branch ends: the flow the network model carries at the SV state. The
  # SV profile itself often holds SvPowerFlow only at injections (the
  # ENTSO-E conformity sets do), and the model flow is what an estimator on
  # this model can reproduce exactly.
  for (t, (bidx, dir)) in result.branch_side_of_terminal
    k = get(pos, bidx, 0)
    k == 0 && continue
    br = net.branchVec[k]
    # only branches closed at both ends (runse! refuses flow measurements
    # on a partially open branch)
    (Sparlectra._branch_terminal_state(br) !== :closed || Int(br.fromBus) in iso || Int(br.toBus) in iso) && continue
    s = (dir === :from ? Sparlectra.branchFlow_pu(br, Int(br.fromBus), Int(br.toBus), 1, V) : Sparlectra.branchFlow_pu(br, Int(br.toBus), Int(br.fromBus), 2, V)) * net.baseMVA
    end_bus = Int(dir === :from ? br.fromBus : br.toBus)
    push!(truth, (:terminal, key(t), "P", real(s), end_bus))
    push!(truth, (:terminal, key(t), "Q", imag(s), end_bus))
  end
  # injection terminals: the delivered SvPowerFlow (power into the equipment)
  if include_injections
    keys_index = cgmesMeasurementKeys(result)
    for sv in objectsOf(store, :SvPowerFlow)
      t = get(sv.refs, :Terminal, nothing)
      t === nothing && continue
      target = get(keys_index.terminals, key(t), nothing)
      (target === nothing || !(target.kind in (:injection, :shunt))) && continue
      p = num(sv, :p)
      q = num(sv, :q)
      target.kind === :injection && p !== nothing && push!(truth, (:terminal, key(t), "P", p, target.busIdx))
      q === nothing || push!(truth, (:terminal, key(t), "Q", q, target.busIdx))
    end
  end
  sort!(truth; by = x -> (x[2], x[3]))

  rng = Random.Xoshiro(seed)
  n = length(truth)
  nbad = round(Int, bad_data_fraction * n)
  bad = Set(Random.randperm(rng, n)[1:nbad])
  Row = NamedTuple{(:timestamp, :measurement_mrid, :terminal_mrid, :node_mrid, :type, :value, :sigma, :source, :valid),Tuple{String,String,String,String,String,Float64,Float64,String,Bool}}
  rows = Row[]
  for (i, (kind, k, typ, v, bus)) in enumerate(truth)
    # the level of the measuring bus first, the keywords where it has none
    lv = Sparlectra._sigma_from_level(levels, net.nodeVec[bus].comp.cVN, typ, v)
    sigma = lv !== nothing ? lv : typ == "U" ? Float64(sigma_u_rel) * abs(v) : max(Float64(sigma_pq_rel) * abs(v), Float64(sigma_pq_floor))
    if i in bad
      value = v + (rand(rng, Bool) ? 1.0 : -1.0) * Float64(bad_data_sigmas) * sigma
      source = "synthetic;bad_data"
    else
      value = noise ? v + sigma * randn(rng) : v
      source = "synthetic"
    end
    push!(rows, Row((String(timestamp), "", kind === :terminal ? k : "", kind === :node ? k : "", typ, value, sigma, source, true)))
  end
  return rows
end

# The CGMES method of Sparlectra.generatePseudoMeasurementTable (docstring
# there): the injection members are terminal mRIDs, each resolved through its
# ConductingEquipment to the prosumer the import built for it.
function Sparlectra.generatePseudoMeasurementTable(result::CGMESImportResult, rows::AbstractVector; sigma::Real, sigma_q::Real = sigma, timestamp::Union{Nothing,AbstractString} = nothing, sigma_by_voltage_level = Sparlectra.state_estimation_config().sigma_by_voltage_level)
  key = Sparlectra._measurement_key
  net = result.net
  schedule = Dict{String,NTuple{2,Union{Nothing,Float64}}}()
  for t in objectsOf(result.store, :Terminal)
    eq = get(t.refs, :ConductingEquipment, nothing)
    eq === nothing && continue
    el = get(result.element_by_mrid, eq, nothing)
    # loads (also asynchronous machines), machines, SVCs, equivalent
    # injections and energy sources are :load or :generator, external grids
    # :external_grid; every other kind is no prosumer
    (el === nothing || !(el[1] in (:load, :generator, :external_grid))) && continue
    schedule[key(t.mrid)] = Sparlectra._scheduled_into_equipment(net.prosumpsVec[el[2]])
  end
  return Sparlectra._pseudo_measurement_rows("generatePseudoMeasurementTable", rows, cgmesMeasurementKeys(result), schedule, sigma, sigma_q, timestamp, sigma_by_voltage_level)
end

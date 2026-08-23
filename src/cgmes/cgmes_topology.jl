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

# file: src/cgmes/cgmes_topology.jl
# purpose: layer 4 of the CGMES importer — bus-branch topology from the TP
# profile. One bus per TopologicalNode; equipment lands on buses via
# Terminal.TopologicalNode. Node-breaker (ConnectivityNode) is a later stage.

"""
Bus-branch view of a CGMES delivery (TP profile present).

- `bus_name`: TopologicalNode mRID → unique Sparlectra bus name (concept
  decision D-3: `IdentifiedObject.name`, mRID-suffixed only on collision).
- `vn_kV`: TopologicalNode mRID → nominal voltage.
- `terminals`: ConductingEquipment mRID → its terminals sorted by
  `sequenceNumber` (1-based positions).
"""
struct CGMESTopology
  bus_name::Dict{String,String}
  vn_kV::Dict{String,Float64}
  terminals::Dict{String,Vector{CIMObject}}
end

function _nominalVoltageOfTN(store::CGMESStore, tn::CIMObject)::Union{Nothing,Float64}
  bv = ref(store, tn, :BaseVoltage)
  bv !== nothing && return num(bv, :nominalVoltage)
  # fallback: TN → ConnectivityNodeContainer (VoltageLevel) → BaseVoltage
  container = ref(store, tn, :ConnectivityNodeContainer)
  if container !== nothing
    bv = ref(store, container, :BaseVoltage)
    bv !== nothing && return num(bv, :nominalVoltage)
  end
  return nothing
end

"""
    buildTopology(store) -> CGMESTopology

Derive the bus set and the equipment→bus lookup from the TP profile.
Deterministic bus naming per decision D-3.
"""
function buildTopology(store::CGMESStore; infer_base_voltages::Bool = false, messages::Union{Nothing,Vector{String}} = nothing)::CGMESTopology
  bus_name = Dict{String,String}()
  vn_kV = Dict{String,Float64}()

  inferred_map, sv_seeded, trafo_seeded, propagated = infer_base_voltages ? _inferNominalVoltages(store) : (Dict{String,Float64}(), 0, 0, 0)
  used_inferred = 0
  inferred_levels = Dict{Float64,Int}()

  # name collisions: count names first, suffix only where needed
  namecount = Dict{String,Int}()
  for tn in objectsOf(store, :TopologicalNode)
    n = something(str(tn, :name), tn.mrid)
    namecount[n] = get(namecount, n, 0) + 1
  end
  for tn in objectsOf(store, :TopologicalNode)
    n = something(str(tn, :name), tn.mrid)
    if namecount[n] > 1
      n = string(n, "_", first(lstrip(tn.mrid, '_'), 8))
    end
    bus_name[tn.mrid] = n
    vn = _nominalVoltageOfTN(store, tn)
    if vn === nothing && infer_base_voltages
      vn = get(inferred_map, tn.mrid, nothing)
      if vn !== nothing
        used_inferred += 1
        inferred_levels[vn] = get(inferred_levels, vn, 0) + 1
      end
    end
    if vn === nothing
      # Same analysis as the mapping-level abort: a missing BaseVoltage
      # catalog is almost always a missing boundary file, and the analysis
      # names the declared dependency instead of leaving a bare node error.
      analysis = importFailureAnalysis(store)
      print(analysis)
      throw(CGMESImportError("CGMES topology: TopologicalNode $(n) ($(tn.mrid)) has no resolvable BaseVoltage — see the import analysis (CLI: printed above; Web UI runs: cgmes.log).", analysis))
    end
    vn_kV[tn.mrid] = vn
  end
  if used_inferred > 0 && messages !== nothing
    # One aggregated warning instead of one line per node: reconstructed
    # nominal voltages are substituted data the user must see, but a
    # thousand-line flood would hide everything else.
    levels_text = join(("$(k) kV×$(v)" for (k, v) in sort!(collect(inferred_levels); by = first, rev = true)), ", ")
    push!(messages, "warning: inferred base voltages for $(used_inferred) topological node(s) (SV-seeded $(sv_seeded), transformer-rated $(trafo_seeded), propagated $(propagated); levels: $(levels_text)) — the BaseVoltage catalog is missing (boundary set); nominal voltages are reconstructed, not delivered")
  end

  terminals = Dict{String,Vector{CIMObject}}()
  for t in objectsOf(store, :Terminal)
    eq = get(t.refs, :ConductingEquipment, nothing)
    eq === nothing && continue
    push!(get!(Vector{CIMObject}, terminals, eq), t)
  end
  for v in values(terminals)
    sort!(v; by = t -> something(num(t, :sequenceNumber), 1.0))
  end
  return CGMESTopology(bus_name, vn_kV, terminals)
end

"""TopologicalNode mRID of terminal `t`, or `nothing`."""
tnOfTerminal(t::CIMObject) = get(t.refs, :TopologicalNode, nothing)

"""SSH `connected` flag of terminal `t` (default `true` when absent)."""
terminalConnected(t::CIMObject)::Bool = something(boolval(t, :connected), true)

"""
Bus name and vn for the `seq`-th terminal of equipment `eq`;
returns `nothing` when the terminal or its TN is missing.
"""
function busOfEquipment(topo::CGMESTopology, eq::CIMObject, seq::Int)
  ts = get(topo.terminals, eq.mrid, nothing)
  (ts === nothing || length(ts) < seq) && return nothing
  tn = tnOfTerminal(ts[seq])
  tn === nothing && return nothing
  haskey(topo.bus_name, tn) || return nothing
  return (bus = topo.bus_name[tn], vn_kV = topo.vn_kV[tn], connected = terminalConnected(ts[seq]), tn = tn)
end

# --- node-breaker topology processor (#314) ---------------------------------

"""
    applyNodeBreakerTopologyProcessor!(store, ctx) -> Bool

Derive the bus partition of a node-breaker delivery WITHOUT a TP profile:
aggregate `ConnectivityNode`s across closed, non-retained switches into
derived topological nodes (union-find over the connectivity graph) and
register them as synthetic `TopologicalNode` objects in the store, wiring
every terminal's `Terminal.TopologicalNode` reference. Downstream stages
(`buildTopology`, the whole element mapping) then run untouched.

Semantics, reused unchanged from the mapping layer:
- switch openness via `_switchIsOpen` (SSH `open` overrides EQ
  `normalOpen`; out of service counts as open), classes `_SWITCH_CLASSES`;
- `Switch.retained`: retained closed switches are NOT merged away, their
  connectivity nodes stay separate buses and the existing switch-link
  mapping turns them into bus links;
- nominal voltages resolve through the existing fallback chain: the
  synthetic node carries its group's `ConnectivityNodeContainer`
  (VoltageLevel) and, when resolvable, the `BaseVoltage` directly;
- a group containing a connectivity node that ALREADY maps to a
  topological node (TP_BD boundary nodes) adopts that node instead of a
  synthetic one.

Runs only when connectivity nodes exist and no non-boundary topological
node does; returns whether it ran. TP-carrying deliveries are untouched.
"""
function applyNodeBreakerTopologyProcessor!(store::CGMESStore, ctx)::Bool
  countOf(store, :ConnectivityNode) == 0 && return false
  any(tn -> !(tn.mrid in store.boundary), objectsOf(store, :TopologicalNode)) && return false

  # union-find over connectivity-node mRIDs (path compression)
  parent = Dict{String,String}()
  for cn in objectsOf(store, :ConnectivityNode)
    parent[cn.mrid] = cn.mrid
  end
  function findroot(x::String)::String
    r = x
    while parent[r] != r
      r = parent[r]
    end
    while parent[x] != r
      parent[x], x = r, parent[x]
    end
    return r
  end

  # terminals per equipment (same association buildTopology uses later)
  eq_terminals = Dict{String,Vector{CIMObject}}()
  for t in objectsOf(store, :Terminal)
    eq = get(t.refs, :ConductingEquipment, nothing)
    eq === nothing && continue
    push!(get!(Vector{CIMObject}, eq_terminals, eq), t)
  end

  n_closed = 0
  n_retained = 0
  n_open = 0
  for cls in _SWITCH_CLASSES
    for sw in objectsOf(store, cls)
      ts = get(eq_terminals, sw.mrid, nothing)
      (ts === nothing || length(ts) < 2) && continue
      cn1 = get(ts[1].refs, :ConnectivityNode, nothing)
      cn2 = get(ts[2].refs, :ConnectivityNode, nothing)
      (cn1 === nothing || cn2 === nothing) && continue
      (haskey(parent, cn1) && haskey(parent, cn2)) || continue
      if _switchIsOpen(ctx, sw)
        n_open += 1
        continue
      end
      if something(boolval(sw, :retained), false)
        # a retained closed switch stays a bus coupler: its nodes remain
        # separate derived buses, the link mapping picks the switch up
        n_retained += 1
        continue
      end
      r1, r2 = findroot(cn1), findroot(cn2)
      if r1 != r2
        # smaller mRID wins as representative: deterministic partition
        parent[max(r1, r2)] = min(r1, r2)
        n_closed += 1
      end
    end
  end

  # groups and their existing topological nodes (TP_BD boundary nodes)
  members = Dict{String,Vector{String}}()
  for cn_mrid in keys(parent)
    push!(get!(Vector{String}, members, findroot(cn_mrid)), cn_mrid)
  end

  # busbar names per group make the derived buses read like the station
  # one-line diagram; fall back to the representative node's name
  busbar_name_of_cn = Dict{String,String}()
  for bb in objectsOf(store, :BusbarSection)
    ts = get(eq_terminals, bb.mrid, nothing)
    ts === nothing && continue
    for t in ts
      cn = get(t.refs, :ConnectivityNode, nothing)
      cn === nothing && continue
      nm = str(bb, :name)
      nm === nothing || (busbar_name_of_cn[cn] = nm)
    end
  end

  group_tn = Dict{String,String}()
  n_synth = 0
  n_adopted = 0
  for (root, cns) in members
    sort!(cns)
    existing = String[]
    for cn_mrid in cns
      cn = store.objects[cn_mrid]
      tn = get(cn.refs, :TopologicalNode, nothing)
      tn === nothing || tn in existing || push!(existing, tn)
    end
    if !isempty(existing)
      sort!(existing)
      length(existing) > 1 && push!(ctx.messages, "warning: topology processor: one switch group maps to $(length(existing)) existing topological nodes; keeping $(first(existing))")
      group_tn[root] = first(existing)
      n_adopted += 1
      continue
    end
    # synthesize a topological node for the group
    tn_mrid = string("_DTN", first(cns))
    while haskey(store.objects, tn_mrid)
      tn_mrid = string(tn_mrid, "x")
    end
    name = nothing
    for cn_mrid in cns
      nm = get(busbar_name_of_cn, cn_mrid, nothing)
      nm === nothing || (name = nm; break)
    end
    name === nothing && (name = something(str(store.objects[first(cns)], :name), first(cns)))
    attrs = Dict{Symbol,String}(:name => name)
    refs = Dict{Symbol,String}()
    # voltage chain: CN -> ConnectivityNodeContainer, where the container
    # may be the VoltageLevel directly or a Bay one hop below it
    # (Bay.VoltageLevel -> VoltageLevel.BaseVoltage). Take the first group
    # member that resolves all the way to a BaseVoltage; keep the best
    # container found even when no member resolves fully, so the error
    # message downstream names a real container.
    for cn_mrid in cns
      cn = store.objects[cn_mrid]
      cont = ref(store, cn, :ConnectivityNodeContainer)
      cont === nothing && continue
      if cont.class === :Bay
        vl = ref(store, cont, :VoltageLevel)
        vl === nothing || (cont = vl)
      end
      haskey(refs, :ConnectivityNodeContainer) || (refs[:ConnectivityNodeContainer] = cont.mrid)
      bv = get(cont.refs, :BaseVoltage, nothing)
      if bv !== nothing
        refs[:ConnectivityNodeContainer] = cont.mrid
        refs[:BaseVoltage] = bv
        break
      end
    end
    obj = CIMObject(tn_mrid, :TopologicalNode, attrs, refs, "topology-processor")
    store.objects[tn_mrid] = obj
    push!(get!(Vector{String}, store.byclass, :TopologicalNode), tn_mrid)
    group_tn[root] = tn_mrid
    n_synth += 1
  end

  # wire the terminals (and the connectivity nodes, for completeness) onto
  # their derived topological nodes; existing TP references are never
  # overwritten
  for t in objectsOf(store, :Terminal)
    haskey(t.refs, :TopologicalNode) && continue
    cn = get(t.refs, :ConnectivityNode, nothing)
    cn === nothing && continue
    haskey(parent, cn) || continue
    t.refs[:TopologicalNode] = group_tn[findroot(cn)]
  end
  for cn in objectsOf(store, :ConnectivityNode)
    haskey(cn.refs, :TopologicalNode) && continue
    cn.refs[:TopologicalNode] = group_tn[findroot(cn.mrid)]
  end

  push!(ctx.messages, "topology processor: derived $(n_synth) topological node(s) from $(length(parent)) connectivity node(s) across $(n_closed) closed switch merge(s) (TP profile absent); $(n_retained) retained closed switch(es) kept as bus couplers, $(n_open) open switch(es) split the graph, $(n_adopted) boundary group(s) adopted their TP_BD node")
  return true
end

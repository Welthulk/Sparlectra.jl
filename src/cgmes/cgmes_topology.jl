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

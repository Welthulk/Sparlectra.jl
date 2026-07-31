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

# file: src/cgmes/cgmes_voltage_inference.jl
# purpose: reconstruct missing nominal voltages (cgmes_import.
#          infer_base_voltages) when a delivery ships without its BaseVoltage
#          catalog (the catalog lives in the boundary EQ of real ENTSO-E
#          deliveries): seed from the SV solved voltages (kV) snapped to the
#          standard level series and from transformer rated voltages, then
#          propagate across level-preserving equipment.

# European standard nominal levels in kV, matched against the level library
# of the real ENTSO-E CGM boundary (0.4 … 750 kV). 380 and 400 (and
# 220/275/300/330, 132/135/150) coexist there — snapping picks the
# relatively nearest, so a 393 kV SV value lands on 400 while 385 lands on
# 380, and a 135 kV system is not folded onto 132.
const _STANDARD_NOMINAL_KV = (1150.0, 750.0, 500.0, 400.0, 380.0, 330.0, 300.0, 275.0, 220.0, 150.0, 135.0, 132.0, 110.0, 66.0, 60.0, 45.0, 33.0, 30.0, 20.0, 15.0, 10.0, 6.3, 3.0, 0.69, 0.4)

# Snap a measured/rated kV value to the standard series. Outside the 15 %
# acceptance band the raw value is kept (rounded) — an exotic level is still
# better than refusing the bus, and the summary warning names the substitution
# either way.
function _snapNominalKV(v::Float64)::Float64
  v > 0.0 || return 0.0
  best = first(_STANDARD_NOMINAL_KV)
  best_rel = abs(v - best) / best
  for level in _STANDARD_NOMINAL_KV
    rel = abs(v - level) / level
    rel < best_rel && (best = level; best_rel = rel)
  end
  return best_rel <= 0.15 ? best : round(v; digits = 1)
end

"""
    _inferNominalVoltages(store::CGMESStore) -> (vn::Dict{String,Float64}, sv_n, trafo_n, prop_n)

Reconstruct nominal voltages for topological nodes whose `BaseVoltage` does
not resolve (missing boundary catalog):

1. **SV seed** — `SvVoltage.v` is stored in kV; snapped to the standard
   level series it identifies the node's level directly (de-energized nodes
   with `v ≈ 0` are skipped).
2. **Transformer seed** — `PowerTransformerEnd.ratedU` at the end's terminal
   fills nodes without an SV value.
3. **Propagation** — nodes joined by level-preserving equipment (anything
   but a `PowerTransformer`) share their level; a BFS spreads the seeded
   levels across switches, lines, and equivalent branches.

Returns the per-TN map plus the per-source counts for the summary message.
Nodes that stay unresolved are simply absent from the map — the caller keeps
its regular missing-voltage handling for them.
"""
function _inferNominalVoltages(store::CGMESStore)
  vn = Dict{String,Float64}()

  # 1. SV seeds — the strongest evidence: the solved operating voltage.
  sv_n = 0
  for sv in objectsOf(store, :SvVoltage)
    tn = get(sv.refs, :TopologicalNode, nothing)
    tn === nothing && continue
    v = num(sv, :v, 0.0)
    v > 0.1 || continue
    haskey(vn, tn) && continue
    snapped = _snapNominalKV(v)
    snapped > 0.0 || continue
    vn[tn] = snapped
    sv_n += 1
  end

  # 2. Transformer rated voltages fill nodes the SV state does not cover.
  trafo_n = 0
  for pte in objectsOf(store, :PowerTransformerEnd)
    ratedU = num(pte, :ratedU, 0.0)
    ratedU > 0.0 || continue
    term_id = get(pte.refs, :Terminal, nothing)
    term_id === nothing && continue
    term = get(store.objects, term_id, nothing)
    term === nothing && continue
    tn = get(term.refs, :TopologicalNode, nothing)
    tn === nothing && continue
    haskey(vn, tn) && continue
    snapped = _snapNominalKV(ratedU)
    snapped > 0.0 || continue
    vn[tn] = snapped
    trafo_n += 1
  end

  # 3. Propagate across level-preserving equipment. Terminals grouped per
  #    conducting equipment; transformers are the only level changers here
  #    (DC converters use DCTerminal objects and never join this adjacency).
  adjacency = Dict{String,Vector{String}}()
  eq_tns = Dict{String,Vector{String}}()
  for term in objectsOf(store, :Terminal)
    eq_id = get(term.refs, :ConductingEquipment, nothing)
    tn = get(term.refs, :TopologicalNode, nothing)
    (eq_id === nothing || tn === nothing) && continue
    eq = get(store.objects, eq_id, nothing)
    eq === nothing && continue
    eq.class == :PowerTransformer && continue
    push!(get!(Vector{String}, eq_tns, eq_id), tn)
  end
  for tns in values(eq_tns)
    length(tns) < 2 && continue
    for a in tns, b in tns
      a == b && continue
      push!(get!(Vector{String}, adjacency, a), b)
    end
  end
  queue = collect(keys(vn))
  prop_n = 0
  while !isempty(queue)
    tn = pop!(queue)
    level = vn[tn]
    for neighbour in get(adjacency, tn, String[])
      haskey(vn, neighbour) && continue
      vn[neighbour] = level
      prop_n += 1
      push!(queue, neighbour)
    end
  end

  return vn, sv_n, trafo_n, prop_n
end

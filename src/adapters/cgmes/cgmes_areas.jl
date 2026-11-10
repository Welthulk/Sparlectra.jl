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

# file: src/adapters/cgmes/cgmes_areas.jl
# purpose: the area of each bus of a CGMES import (cgmes_import.area_source,
#          1.2.0): the operator model of a CGM or the region of the
#          substation, the boundary nodes by their neighbours (independent
#          of the file order), and the
#          check of the result against ControlArea/TieFlow

# The container of a TN as the EQ defines it; a Bay stands for its
# VoltageLevel. The topology processor gives its synthetic TNs the container
# of their connectivity nodes, so node-breaker deliveries resolve the same way.
function _cgmesAreaContainer(store::CGMESStore, tn::CIMObject)::Union{Nothing,CIMObject}
  c = ref(store, tn, :ConnectivityNodeContainer)
  (c !== nothing && c.class === :Bay) && (c = something(ref(store, c, :VoltageLevel), c))
  return c
end

# The area of one TN as (key, name), or `nothing` when the TN has none of
# its own (a boundary node, a container without region). `key` groups the
# buses, `name` goes to net.area_names. Nothing here depends on the order of
# the files: the areas are numbered by name (`_assignCGMESAreas!`).
#
# model_authority: the file that defines the container of the TN. In a CGM
# every operator ships its own EQ, so that file names the operator
# (Model.modelingAuthoritySet). The TP file of the same operator is the
# fallback when the container is missing. A file without the header entry
# stands for itself (its name without extension).
function _cgmesAreaKey(store::CGMESStore, tn::CIMObject, source::Symbol, authority::Dict{String,String})
  c = _cgmesAreaContainer(store, tn)
  if source === :model_authority
    src = if c !== nothing && !(c.mrid in store.boundary)
      c.source
    elseif !(tn.mrid in store.boundary)
      tn.source
    else
      nothing
    end
    (src === nothing || !haskey(authority, src)) && return nothing
    a = authority[src]
    name = isempty(a) ? basename_short(src) : a
    return (key = name, name = name)
  end
  (c === nothing || c.mrid in store.boundary) && return nothing
  sub = c.class === :VoltageLevel ? ref(store, c, :Substation) : (c.class === :Substation ? c : nothing)
  sub === nothing && return nothing
  # Substation.Region is the SubGeographicalRegion, its Region the
  # GeographicalRegion
  reg = ref(store, sub, :Region)
  (source === :geographical_region && reg !== nothing) && (reg = ref(store, reg, :Region))
  reg === nothing && return nothing
  return (key = reg.mrid, name = something(str(reg, :name), reg.mrid))
end

# Neighbours of every bus: the other ends of its branches and links.
function _cgmesAreaNeighbours(net)::Vector{Vector{Int}}
  adj = [Int[] for _ in net.nodeVec]
  for br in net.branchVec
    f, t = Int(br.fromBus), Int(br.toBus)
    push!(adj[f], t)
    push!(adj[t], f)
  end
  for l in net.linkVec
    f, t = Int(l.fromBus), Int(l.toBus)
    push!(adj[f], t)
    push!(adj[t], f)
  end
  return adj
end

"""
Set the area of every bus of a CGMES import (`cgmes_import.area_source`).

Each TN takes the area of its container (`_cgmesAreaKey`). The areas are
numbered from 1 in the order of their names (the mRID breaks a tie between
two regions of one name), the names go to `net.area_names`. A bus without
an area of its own (a boundary X node, a per-side split bus, a transformer
star point, a container without region) takes, among the buses at the
other ends of its branches and links that have an area, the area with the
smallest number, i.e. the smallest name. A chain of such buses resolves
step by step, each step from the areas of the step before. So neither
the numbers nor the X nodes depend on the order of the files or the buses.
The import lists every bus assigned this way in one notice. With `:none`
nothing is set. Returns nothing; the result is described in `ctx.messages`.
"""
function _assignCGMESAreas!(net, store::CGMESStore, topo::CGMESTopology, ctx::_MapCtx, source::Symbol)
  source === :none && return nothing
  authority = Dict{String,String}(f.name => f.modeling_authority for f in store.files)
  key_of_bus = Dict{Int,String}()
  name_of = Dict{String,String}()   # key => name
  for (tnm, bus) in topo.bus_name
    k = get(net.busDict, bus, nothing)
    tn = get(store.objects, tnm, nothing)
    (k === nothing || tn === nothing) && continue
    r = _cgmesAreaKey(store, tn, source, authority)
    r === nothing && continue
    key_of_bus[k] = r.key
    name_of[r.key] = r.name
  end
  if isempty(name_of)
    push!(ctx.messages, "notice: areas (cgmes_import.area_source = $(source)): the delivery names none, no bus area set")
    return nothing
  end
  order = sort!(collect(keys(name_of)); by = k -> (name_of[k], k))
  number = Dict{String,Int}(k => i for (i, k) in enumerate(order))
  area = zeros(Int, length(net.nodeVec))
  for (k, key) in key_of_bus
    area[k] = number[key]
  end
  # buses without an area of their own: the smallest area among the
  # neighbours that had one before this step
  adj = _cgmesAreaNeighbours(net)
  taken = Tuple{String,Int}[]
  while true
    step = Tuple{Int,Int}[]
    for k in eachindex(area)
      area[k] == 0 || continue
      cand = [area[b] for b in adj[k] if area[b] != 0]
      isempty(cand) || push!(step, (k, minimum(cand)))
    end
    isempty(step) && break
    for (k, a) in step
      area[k] = a
      push!(taken, (Sparlectra.getCompName(net.nodeVec[k].comp), a))
    end
  end
  for k in eachindex(area)
    area[k] == 0 || (net.nodeVec[k]._area = area[k])
  end
  empty!(net.area_names)
  for (i, key) in enumerate(order)
    net.area_names[i] = name_of[key]
  end
  counts = [count(==(i), area) for i in eachindex(order)]
  push!(ctx.messages, "notice: areas (cgmes_import.area_source = $(source)): $(length(order)) area(s), " * join(["$(i) $(net.area_names[i]) ($(counts[i]) buses)" for i in eachindex(order)], ", "))
  if !isempty(taken)
    sort!(taken)
    push!(ctx.messages, "notice: $(length(taken)) bus(es) without an area of their own (boundary node, split side, star point) took the smallest area among their neighbours: " * join(["$(b) -> $(a)" for (b, a) in taken], ", "))
  end
  n_none = count(==(0), area)
  n_none > 0 && push!(ctx.messages, "notice: $(n_none) bus(es) without an area: no container with $(source) and no neighbour with an area")
  (length(order) == 1 && source === :model_authority) && push!(ctx.messages, "notice: one model authority, all buses in one area; cgmes_import.area_source = sub_geographical_region splits by region")
  _cgmesTieFlowCheck!(ctx, net, store, topo, area)
  return nothing
end

# ControlArea/TieFlow, where the delivery has them, only as a check: every
# TieFlow terminal should sit on a border bus of the areas (an end of a closed
# branch between two areas). A difference is reported, nothing is changed.
function _cgmesTieFlowCheck!(ctx::_MapCtx, net, store::CGMESStore, topo::CGMESTopology, area::Vector{Int})
  countOf(store, :TieFlow) == 0 && return nothing
  border = Set{Int}()
  for br in net.branchVec
    Sparlectra._branch_terminal_state(br) === :closed || continue
    f, t = Int(br.fromBus), Int(br.toBus)
    if area[f] != 0 && area[t] != 0 && area[f] != area[t]
      push!(border, f)
      push!(border, t)
    end
  end
  n = 0
  n_off = 0
  off = String[]
  for tf in objectsOf(store, :TieFlow)
    t = ref(store, tf, :Terminal)
    t === nothing && continue
    tn = tnOfTerminal(t)
    (tn === nothing || !haskey(topo.bus_name, tn)) && continue
    k = get(net.busDict, topo.bus_name[tn], nothing)
    k === nothing && continue
    n += 1
    k in border && continue
    n_off += 1
    push!(off, topo.bus_name[tn])
  end
  n == 0 && return nothing
  unique!(off)
  msg = "notice: areas against ControlArea/TieFlow: $(n - n_off) of $(n) TieFlow terminal(s) at a border bus of the areas"
  isempty(off) || (msg *= "; not at a border: " * join(first(off, 5), ", ") * (length(off) > 5 ? ", ..." : ""))
  push!(ctx.messages, msg)
  return nothing
end

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
# file: src/adapters/cgmes/cgmes_union.jl
# purpose: the union network of a CGMES snapshot series, store side: the
#          node partition over all snapshots (terminal atoms), the switch
#          links with their state per step, the extra `_TPL` couplings, and
#          the union store the unchanged importer maps once.

# Node identity works on TERMINALS, not on connectivity nodes: a TP may carry
# the TN only at the terminal (no ConnectivityNode.TopologicalNode), and
# terminals of one CN may then sit on different TNs (measured on real
# deliveries: equivalent networks re-attached per export). The importer reads
# Terminal.TopologicalNode as well (busOfEquipment), so terminal atoms are
# what the direct import of each snapshot sees.

"""Raised when a snapshot series cannot be turned into a union network."""
struct CGMESSeriesError <: Exception
  msg::String
end
Base.showerror(io::IO, e::CGMESSeriesError) = print(io, "CGMESSeriesError: ", e.msg)

# Topology view of one snapshot, taken from its store and then kept alone,
# so the store itself can be dropped.
struct _SeriesTopo
  label::String
  term_tn::Dict{String,String}        # terminal mRID -> TN mRID of this snapshot's TP
  term_cn::Dict{String,String}        # terminal mRID -> ConnectivityNode mRID
  eq_terms::Dict{String,Vector{String}}  # equipment mRID -> terminal mRIDs by sequence number
  links::Set{String}                  # switches the direct import turns into a bus link
  tn_objects::Dict{String,CIMObject}  # TN mRID -> TN object (name, BaseVoltage, container)
  sv_voltage::Dict{String,CIMObject}  # TN mRID -> its SvVoltage
end

_tpTNs(store::CGMESStore) = [tn for tn in objectsOf(store, :TopologicalNode) if !(tn.mrid in store.boundary)]

function _seriesTopo(store::CGMESStore, label::String)::_SeriesTopo
  # a snapshot without a TP gets its partition from the node-breaker
  # processor, applied to a copy so the caller's store stays untouched
  s = store
  if isempty(_tpTNs(store))
    s = deepcopy(store)
    applyNodeBreakerTopologyProcessor!(s, _MapCtx())
    isempty(_tpTNs(s)) && throw(CGMESSeriesError("snapshot $(label): no TP and no node-breaker connectivity, no node partition"))
  end
  cn_tn = Dict{String,String}()
  for cn in objectsOf(s, :ConnectivityNode)
    tn = get(cn.refs, :TopologicalNode, nothing)
    tn === nothing || (cn_tn[cn.mrid] = tn)
  end
  term_tn = Dict{String,String}(); term_cn = Dict{String,String}()
  seqd = Dict{String,Vector{Tuple{Float64,String}}}()
  for t in objectsOf(s, :Terminal)
    cn = get(t.refs, :ConnectivityNode, nothing)
    cn === nothing || (term_cn[t.mrid] = cn)
    # the terminal's own TN first, as the importer reads it
    tn = something(get(t.refs, :TopologicalNode, nothing), cn === nothing ? nothing : get(cn_tn, cn, nothing), Some(nothing))
    tn === nothing || (term_tn[t.mrid] = tn)
    eq = get(t.refs, :ConductingEquipment, nothing)
    eq === nothing || push!(get!(Vector{Tuple{Float64,String}}, seqd, eq), (something(num(t, :sequenceNumber), 1.0), t.mrid))
  end
  eq_terms = Dict(k => [m for (_, m) in sort!(v)] for (k, v) in seqd)
  # the importer's link rule (_mapSwitches!): closed by _switchIsOpen, both
  # terminals connected, the two TNs differ
  ctx = _MapCtx()
  links = Set{String}()
  for cls in _SWITCH_CLASSES, sw in objectsOf(s, cls)
    ts = get(eq_terms, sw.mrid, String[])
    length(ts) >= 2 || continue
    tn1 = get(term_tn, ts[1], nothing); tn2 = get(term_tn, ts[2], nothing)
    (tn1 === nothing || tn2 === nothing || tn1 == tn2) && continue
    _switchIsOpen(ctx, sw) && continue
    all(m -> terminalConnected(s.objects[m]), ts[1:2]) || continue
    push!(links, sw.mrid)
  end
  tn_objects = Dict(tn.mrid => tn for tn in objectsOf(s, :TopologicalNode))
  sv_voltage = Dict{String,CIMObject}()
  for sv in objectsOf(s, :SvVoltage)
    tn = get(sv.refs, :TopologicalNode, nothing)
    tn === nothing || (sv_voltage[tn] = sv)
  end
  return _SeriesTopo(label, term_tn, term_cn, eq_terms, links, tn_objects, sv_voltage)
end

# The node partition of a series. A union node is the set of terminals that
# share their TN in every snapshot. A terminal absent from a snapshot takes,
# there, the TN its CN has when that CN is not split; otherwise it is free
# in that snapshot ("" in the signature).
struct _UnionPartition
  node_of_term::Dict{String,Int}
  node_sig::Vector{Vector{String}}    # node -> TN per step ("" = free)
  node_terms::Vector{Vector{String}}  # node -> its terminals, sorted
  node_mrid::Vector{String}           # mRID of the union TN
end

function _unionPartition(topos::Vector{_SeriesTopo}, base::Int)::_UnionPartition
  n = length(topos)
  # per step: CN -> its TN when all its terminals agree
  cn_unique = Vector{Dict{String,String}}(undef, n)
  for k in 1:n
    seen = Dict{String,Union{Nothing,String}}()
    for (t, tn) in topos[k].term_tn
      cn = get(topos[k].term_cn, t, nothing)
      cn === nothing && continue
      prev = get(seen, cn, "")
      seen[cn] = prev == "" ? tn : (prev == tn ? tn : nothing)
    end
    cn_unique[k] = Dict{String,String}(cn => tn for (cn, tn) in seen if tn !== nothing)
  end
  terms = sort!(collect(union((keys(t.term_tn) for t in topos)...)))
  cn_of(t) = (for k in 1:n; c = get(topos[k].term_cn, t, nothing); c === nothing || return c; end; nothing)
  node_by_sig = Dict{Vector{String},Int}()
  node_of_term = Dict{String,Int}()
  sigs = Vector{Vector{String}}()
  members = Vector{Vector{String}}()
  for t in terms
    cn = cn_of(t)
    sig = [haskey(topos[k].term_tn, t) ? topos[k].term_tn[t] : (cn === nothing ? "" : get(cn_unique[k], cn, "")) for k in 1:n]
    id = get(node_by_sig, sig, 0)
    if id == 0
      push!(sigs, sig); push!(members, String[])
      id = length(sigs); node_by_sig[sig] = id
    end
    push!(members[id], t)
    node_of_term[t] = id
  end
  # mRID: the base TN when the node is exactly that TN's terminal set (keeps
  # bus names and TN keys of an unswitched series), else _UTN + the smallest
  # terminal mRID, after the _DTN pattern of the node-breaker processor
  base_terms = Dict{String,Int}()
  for tn in values(topos[base].term_tn)
    base_terms[tn] = get(base_terms, tn, 0) + 1
  end
  mrids = Vector{String}(undef, length(sigs))
  for id in eachindex(sigs)
    tn = sigs[id][base]
    whole = tn != "" && get(base_terms, tn, 0) == length(members[id]) && all(t -> get(topos[base].term_tn, t, "") == tn, members[id])
    mrids[id] = whole ? tn : string("_UTN", lstrip(first(members[id]), '_'))
  end
  return _UnionPartition(node_of_term, sigs, members, mrids)
end

# first terminal list with two entries of equipment `m` over the steps
function _seriesEqTerms(topos::Vector{_SeriesTopo}, m::String)::Vector{String}
  for t in topos
    ts = get(t.eq_terms, m, String[])
    length(ts) >= 2 && return ts
  end
  return String[]
end

# Split union nodes at the switches in `extra` (keyword extra_link_switches of
# buildCGMESSeriesUnion). A switch with both terminals in one union node (one
# TN in every step) never becomes a link, so no event could open it. Inside
# such a node the terminals are joined by their CN and by every other switch
# with both terminals in the node, open or closed; the parts left without the
# named switches become union nodes of their own, each named _UTN + its
# smallest terminal like any node that is not a whole base TN, with the
# signature of the node. A named switch between two union nodes needs no
# split. Throws CGMESSeriesError for a name that is not a switch of the
# series, a switch without two terminals on union nodes, a boundary node, a
# terminal without CN in the node, terminals the CNs and switches do not
# join, and a switch whose two sides stay joined by another switch.
function _splitPartition(part::_UnionPartition, topos::Vector{_SeriesTopo}, switch_class::Dict{String,Symbol}, extra::Set{String}, boundary::Set{String})::_UnionPartition
  isempty(extra) && return part
  # node -> the named switches inside it
  inside = Dict{Int,Vector{String}}()
  for m in sort!(collect(extra))
    haskey(switch_class, m) || throw(CGMESSeriesError("extra_link_switches: $(m) is not a switch of the series"))
    ts = _seriesEqTerms(topos, m)
    u1 = isempty(ts) ? 0 : get(part.node_of_term, ts[1], 0)
    u2 = isempty(ts) ? 0 : get(part.node_of_term, ts[2], 0)
    (u1 == 0 || u2 == 0) && throw(CGMESSeriesError("extra_link_switches: switch $(m) has no two terminals on topological nodes"))
    u1 == u2 && push!(get!(Vector{String}, inside, u1), m)
  end
  isempty(inside) && return part
  # the other switches with both terminals in a node to split
  others = Dict{Int,Vector{Vector{String}}}()
  for m in sort!(collect(keys(switch_class)))
    m in extra && continue
    ts = _seriesEqTerms(topos, m)
    isempty(ts) && continue
    u = get(part.node_of_term, ts[1], 0)
    (haskey(inside, u) && get(part.node_of_term, ts[2], 0) == u) && push!(get!(Vector{Vector{String}}, others, u), ts[1:2])
  end
  cn_of(t) = (for tp in topos; c = get(tp.term_cn, t, nothing); c === nothing || return c; end; nothing)
  node_of_term = copy(part.node_of_term)
  sigs = copy(part.node_sig)
  members = copy(part.node_terms)
  mrids = copy(part.node_mrid)
  for u in sort!(collect(keys(inside)))
    named = inside[u]
    part.node_mrid[u] in boundary && throw(CGMESSeriesError("extra_link_switches: switch $(first(named)) sits in boundary node $(part.node_mrid[u])"))
    ts_u = part.node_terms[u]
    pos = Dict(t => i for (i, t) in enumerate(ts_u))
    # union-find over the terminals of the node: `full` with the named
    # switches, `cut` without them
    full = collect(eachindex(ts_u))
    cut = collect(eachindex(ts_u))
    join!(p, a, b) = (ra = _uf_root!(p, a); rb = _uf_root!(p, b); ra == rb || (p[ra] = rb); nothing)
    first_of_cn = Dict{String,Int}()
    for (i, t) in enumerate(ts_u)
      cn = cn_of(t)
      cn === nothing && throw(CGMESSeriesError("extra_link_switches: switch $(first(named)): terminal $(t) of its node has no ConnectivityNode, the split needs node-breaker connectivity"))
      j = get!(first_of_cn, cn, i)
      join!(full, i, j); join!(cut, i, j)
    end
    for ts in get(others, u, Vector{String}[])
      join!(full, pos[ts[1]], pos[ts[2]]); join!(cut, pos[ts[1]], pos[ts[2]])
    end
    for m in named
      ts = _seriesEqTerms(topos, m)
      join!(full, pos[ts[1]], pos[ts[2]])
    end
    length(Set(_uf_root!(full, i) for i in eachindex(ts_u))) == 1 ||
      throw(CGMESSeriesError("extra_link_switches: switch $(first(named)): the terminals of its node are not joined by ConnectivityNodes and switches alone"))
    for m in named
      ts = _seriesEqTerms(topos, m)
      _uf_root!(cut, pos[ts[1]]) == _uf_root!(cut, pos[ts[2]]) &&
        throw(CGMESSeriesError("extra_link_switches: switch $(m): another switch joins its two sides, name that switch as well"))
    end
    # the parts in the order of their smallest terminal (ts_u is sorted); the
    # first keeps the node id, the others are appended
    parts = Dict{Int,Vector{String}}()
    order = Int[]
    for (i, t) in enumerate(ts_u)
      r = _uf_root!(cut, i)
      haskey(parts, r) || push!(order, r)
      push!(get!(Vector{String}, parts, r), t)
    end
    for (j, r) in enumerate(order)
      terms = parts[r]
      if j == 1
        id = u
        members[u] = terms
        mrids[u] = string("_UTN", lstrip(first(terms), '_'))
      else
        push!(sigs, copy(part.node_sig[u]))
        push!(members, terms)
        push!(mrids, string("_UTN", lstrip(first(terms), '_')))
        id = length(sigs)
      end
      for t in terms
        node_of_term[t] = id
      end
    end
  end
  return _UnionPartition(node_of_term, sigs, members, mrids)
end

# Switch links of the union with their state per step, and the extra
# couplings (_TPL) a step needs where two union nodes share its TN without a
# closed switch between them. A switch in `keep` (extra_link_switches) is a
# link even when it conducts in no step.
struct _UnionLinks
  switches::Vector{String}            # switch mRIDs that become links
  switch_nodes::Vector{Tuple{Int,Int}}
  switch_status::Vector{Vector{Int}}  # per switch: status per step
  tpl_nodes::Vector{Tuple{Int,Int}}
  tpl_status::Vector{Vector{Int}}
end

function _uf_root!(p::Vector{Int}, x::Int)
  while p[x] != x
    p[x] = p[p[x]]
    x = p[x]
  end
  return x
end

function _unionLinks(topos::Vector{_SeriesTopo}, part::_UnionPartition, switch_class::Dict{String,Symbol}, keep::Set{String} = Set{String}())::_UnionLinks
  n = length(topos)
  sws = String[]; swn = Tuple{Int,Int}[]; sst = Vector{Vector{Int}}()
  for m in sort!(collect(keys(switch_class)))
    ts = String[]
    for k in 1:n
      ts = get(topos[k].eq_terms, m, String[])
      length(ts) >= 2 && break
    end
    length(ts) >= 2 || continue
    u1 = get(part.node_of_term, ts[1], 0); u2 = get(part.node_of_term, ts[2], 0)
    (u1 == 0 || u2 == 0 || u1 == u2) && continue
    # TP is the authority: a link of the direct import (held, closed), or
    # both terminals on one TN of that step; a switch absent from a step
    # does not conduct there
    st = [m in topos[k].links || (haskey(topos[k].term_tn, ts[1]) && get(topos[k].term_tn, ts[2], "") == topos[k].term_tn[ts[1]]) ? 1 : 0 for k in 1:n]
    (any(==(1), st) || m in keep) || continue
    push!(sws, m); push!(swn, (u1, u2)); push!(sst, st)
  end
  nn = length(part.node_mrid)
  tpl = Dict{Tuple{Int,Int},Vector{Int}}()
  for k in 1:n
    p = collect(1:nn)
    for (i, (a, b)) in enumerate(swn)
      sst[i][k] == 1 || continue
      ra = _uf_root!(p, a); rb = _uf_root!(p, b)
      ra == rb || (p[ra] = rb)
    end
    # union nodes per TN of this step, by component
    comps = Dict{String,Dict{Int,Int}}()   # TN -> component root -> smallest node id
    for id in 1:nn
      tn = part.node_sig[id][k]
      tn == "" && continue
      r = _uf_root!(p, id)
      d = get!(Dict{Int,Int}, comps, tn)
      d[r] = min(get(d, r, id), id)
    end
    for d in values(comps)
      length(d) > 1 || continue
      reps = sort!(collect(values(d)))
      # chain the components by their smallest node; a pair needed in
      # several steps is one coupling with several closed steps
      for j in 2:length(reps)
        key = (reps[1], reps[j])
        v = get!(() -> zeros(Int, n), tpl, key)
        v[k] = 1
      end
    end
  end
  keys_sorted = sort!(collect(keys(tpl)))
  return _UnionLinks(sws, swn, sst, keys_sorted, [tpl[k] for k in keys_sorted])
end

# TN-bound objects: the union store rebuilds them from the partition, so no
# snapshot contributes its own
const _UNION_REBUILT = (:TopologicalNode, :SvVoltage)

# profile set -> file name in the base store. Objects from other snapshots
# keep the base file name of their profile as source: the importer splits
# boundary nodes by source file, and a time stamp in a file name must not
# read as a second area.
_unionBaseFiles(b::CGMESStore) = Dict{Set{Symbol},String}(f.profiles => f.name for f in b.files if !f.skipped)

# Objects the other snapshots add to the base store, collected one snapshot
# at a time so each snapshot store can be dropped after its turn. `added`
# keeps the order of collection with the byclass entry, `multirefs` the
# order of the reference lists; _unionStore replays both in that order, so
# the union store does not depend on whether the stores were held together.
struct _UnionExtras
  objects::Dict{String,CIMObject}                               # key (mRID or `<mRID>#<class>`) -> own copy
  added::Vector{Tuple{Symbol,String}}                           # (byclass entry, key) in order of collection
  multirefs::Vector{Pair{Tuple{String,Symbol},Vector{String}}}  # in order of collection
  alias::Dict{String,String}                                    # `<mRID>#<class>` -> mRID
end
_UnionExtras() = _UnionExtras(Dict{String,CIMObject}(), Tuple{Symbol,String}[], Pair{Tuple{String,Symbol},Vector{String}}[], Dict{String,String}())

# Collect from snapshot store `s` every object that neither the base store
# `b` nor an earlier step has. Steps come in time order, so the first
# appearance wins. An exporter may reuse a mRID for an object of another
# class (measured: a GeneratingUnit that came back as an
# EquivalentInjection); that object enters under `<mRID>#<class>`, and the
# references of its snapshot follow it. Only copies go into `ex`; `s` and
# `b` are not modified.
function _collectExtras!(ex::_UnionExtras, b::CGMESStore, base_file::Dict{Set{Symbol},String}, s::CGMESStore)::_UnionExtras
  file_map = Dict{String,String}(f.name => get(base_file, f.profiles, f.name) for f in s.files)
  renamed = Dict{String,String}()
  for (m, o) in s.objects
    u = get(b.objects, m, nothing)
    u === nothing && (u = get(ex.objects, m, nothing))
    (u === nothing || u.class == o.class) && continue
    renamed[m] = string(m, "#", o.class)
  end
  added = Set{String}()
  for cls in sort!(collect(keys(s.byclass)))
    cls in _UNION_REBUILT && continue
    for m in s.byclass[cls]
      key = get(renamed, m, m)
      (haskey(b.objects, key) || haskey(ex.objects, key)) && continue
      o = s.objects[m]
      refs = isempty(renamed) ? copy(o.refs) : Dict{Symbol,String}(r => get(renamed, v, v) for (r, v) in o.refs)
      ex.objects[key] = CIMObject(key, o.class, copy(o.attrs), refs, get(file_map, o.source, o.source))
      push!(ex.added, (cls, key))
      push!(added, m)
      key == m || (ex.alias[key] = m)
    end
  end
  for ((mm, r), v) in s.multirefs
    mm in added && push!(ex.multirefs, (get(renamed, mm, mm), r) => [get(renamed, x, x) for x in v])
  end
  return ex
end

# The union store: the base snapshot's store `b`, plus the objects the other
# snapshots add (`ex`, see _collectExtras!), with the TNs replaced by union
# TNs, terminals rewired, link switches closed, the never conducting ones
# opened, and the _TPL couplings added as closed breakers. `b` is not
# modified; the objects of `ex` move into the union store, so `ex` is used
# up. Returns the store and the alias map of `ex`.
function _unionStore(b::CGMESStore, ex::_UnionExtras, topos::Vector{_SeriesTopo}, part::_UnionPartition, links::_UnionLinks, base::Int)::Tuple{CGMESStore,Dict{String,String}}
  objects = Dict{String,CIMObject}(m => CIMObject(o.mrid, o.class, copy(o.attrs), copy(o.refs), o.source) for (m, o) in b.objects)
  byclass = Dict{Symbol,Vector{String}}(c => copy(v) for (c, v) in b.byclass)
  multirefs = deepcopy(b.multirefs)
  base_file = _unionBaseFiles(b)
  # the other snapshots' objects in the order of collection, into each
  # container in the same sequence as when all stores were held at once
  for (cls, key) in ex.added
    objects[key] = ex.objects[key]
    push!(get!(Vector{String}, byclass, cls), key)
  end
  for (key, v) in ex.multirefs
    multirefs[key] = v
  end
  # drop the snapshot TNs and their SvVoltages (boundary TNs stay)
  for cls in _UNION_REBUILT
    keep = String[]
    for m in get(byclass, cls, String[])
      o = objects[m]
      tn = cls === :TopologicalNode ? m : get(o.refs, :TopologicalNode, "")
      if tn in b.boundary
        push!(keep, m)
      else
        delete!(objects, m)
      end
    end
    byclass[cls] = keep
  end
  # union TNs with name, base voltage and container of the TN they came from
  # (base step first), and a copy of that TN's SvVoltage as start value
  for (id, mrid) in enumerate(part.node_mrid)
    haskey(objects, mrid) && continue      # a boundary TN kept as it is
    k0 = part.node_sig[id][base] != "" ? base : findfirst(k -> part.node_sig[id][k] != "", eachindex(topos))
    k0 === nothing && continue
    rep = get(topos[k0].tn_objects, part.node_sig[id][k0], nothing)
    rep === nothing && throw(CGMESSeriesError("snapshot $(topos[k0].label): terminal names TN $(part.node_sig[id][k0]), which the TP does not define"))
    objects[mrid] = CIMObject(mrid, :TopologicalNode, copy(rep.attrs), copy(rep.refs), get(base_file, Set([:TP]), rep.source))
    push!(byclass[:TopologicalNode], mrid)
    sv = get(topos[k0].sv_voltage, part.node_sig[id][k0], nothing)
    if sv !== nothing
      svm = string(mrid, "_SV")
      refs = copy(sv.refs); refs[:TopologicalNode] = mrid
      objects[svm] = CIMObject(svm, :SvVoltage, copy(sv.attrs), refs, sv.source)
      push!(get!(Vector{String}, byclass, :SvVoltage), svm)
    end
  end
  # rewire terminals; CN -> TN references no longer describe one node
  for m in get(byclass, :Terminal, String[])
    t = objects[m]
    id = get(part.node_of_term, m, 0)
    if id == 0
      delete!(t.refs, :TopologicalNode)
    else
      t.refs[:TopologicalNode] = part.node_mrid[id]
    end
  end
  for m in get(byclass, :ConnectivityNode, String[])
    m in b.boundary || delete!(objects[m].refs, :TopologicalNode)
  end
  # the angle reference of each base island moves to the union node that
  # carries most terminals of the base TN
  largest = Dict{String,Tuple{Int,String}}()
  for (id, sig) in enumerate(part.node_sig)
    tn = sig[base]; tn == "" && continue
    nterm = length(part.node_terms[id])
    (haskey(largest, tn) && largest[tn][1] >= nterm) || (largest[tn] = (nterm, part.node_mrid[id]))
  end
  for m in get(byclass, :TopologicalIsland, String[])
    isl = objects[m]
    tn = get(isl.refs, :AngleRefTopologicalNode, nothing)
    (tn === nothing || tn in b.boundary) && continue
    haskey(largest, tn) ? (isl.refs[:AngleRefTopologicalNode] = largest[tn][2]) : delete!(isl.refs, :AngleRefTopologicalNode)
    # the member list names snapshot TNs; nothing reads it after the rewiring
    for key in collect(keys(multirefs))
      key[1] == m && delete!(multirefs, key)
    end
    delete!(isl.refs, :TopologicalNodes)
  end
  # switches: links closed, everything else between two union nodes open,
  # so the importer builds exactly the union's links
  linkset = Set(links.switches)
  for cls in _SWITCH_CLASSES, m in get(byclass, cls, String[])
    sw = objects[m]
    if m in linkset
      sw.attrs[:open] = "false"
    else
      ts = String[]
      for k in eachindex(topos)
        ts = get(topos[k].eq_terms, m, String[])
        length(ts) >= 2 && break
      end
      length(ts) >= 2 || continue
      get(part.node_of_term, ts[1], 0) != get(part.node_of_term, ts[2], 0) && (sw.attrs[:open] = "true")
    end
  end
  # _TPL couplings: closed breakers with two terminals on their union nodes
  for (a, c) in links.tpl_nodes
    m = string("_TPL", lstrip(part.node_mrid[a], '_'), "_", lstrip(part.node_mrid[c], '_'))
    objects[m] = CIMObject(m, :Breaker, Dict{Symbol,String}(:name => "_TPL", :open => "false"), Dict{Symbol,String}(), "series-union")
    push!(get!(Vector{String}, byclass, :Breaker), m)
    for (seq, node) in ((1, a), (2, c))
      tm = string(m, "_T", seq)
      objects[tm] = CIMObject(tm, :Terminal, Dict{Symbol,String}(:sequenceNumber => string(seq), :connected => "true"),
                              Dict{Symbol,String}(:ConductingEquipment => m, :TopologicalNode => part.node_mrid[node]), "series-union")
      push!(get!(Vector{String}, byclass, :Terminal), tm)
    end
  end
  # every regulator the series may switch on exists in the union; the delta
  # disables it per step
  for cls in (:RatioTapChanger, :PhaseTapChangerLinear, :PhaseTapChangerSymmetrical, :PhaseTapChangerAsymmetrical, :PhaseTapChangerTabular)
    for m in get(byclass, cls, String[])
      objects[m].attrs[:controlEnabled] = "true"
    end
  end
  for m in get(byclass, :TapChangerControl, String[])
    objects[m].attrs[:enabled] = "true"
  end
  return CGMESStore(objects, byclass, copy(b.files), b.version, copy(b.boundary), multirefs), ex.alias
end

# --- delta of one snapshot against the template ---------------------------

"""
    CGMESSnapshotDelta

One step of a snapshot series expressed against the union template: `ops`
(patch operations of the `:cgmes_series` vocabulary), `control_ops` (ratio
and angle of transformers whose regulator is active in this step; the
sequential mode leaves them to the controllers), `start_ops` (`vm_start`
and `va_start` of every bus the snapshot keeps in service: the start value
of its direct import, the SV voltage or, without one, the neighbour start of
the importer; QSTS mode A starts each step from them), `ssh_steps`
(TapChanger and ShuntCompensator positions of the SSH, reference only) and
`residual` (differences the vocabulary cannot express; empty means the step
is reproduced).
"""
struct CGMESSnapshotDelta
  step::Int
  label::String
  scenario_time::Union{Nothing,DateTime}
  ops::Vector{Sparlectra.PatchOp}
  control_ops::Vector{Sparlectra.PatchOp}
  start_ops::Vector{Sparlectra.PatchOp}
  ssh_steps::Dict{String,Float64}
  residual::Vector{String}
end

const _DELTA_RTOL = 1e-12

# relative comparison with an absolute floor of the same size; `nothing`
# only equals `nothing`
function _delta_differs(a, b)::Bool
  (a === nothing || b === nothing) && return a !== b
  (a isa Bool || b isa Bool) && return a != b
  fa = Float64(a); fb = Float64(b)
  (isnan(fa) || isnan(fb)) && return isnan(fa) != isnan(fb)
  return abs(fa - fb) > _DELTA_RTOL * max(abs(fa), abs(fb), 1.0)
end

# value of an optional quantity in an op: no value travels as `nothing`
_delta_value(x) = x === nothing ? nothing : Float64(x)

# Whether two tap data values differ, through vectors, tuples and the
# immutable records (BranchTapData, TapStepImpedance); `differs` compares
# two numbers: _delta_differs for the delta, !isequal for reusing a record
function _tap_values_differ(a, b, differs)::Bool
  (a === nothing || b === nothing) && return a !== b
  if a isa Number && b isa Number
    return differs(a, b)
  elseif a isa Union{Symbol,AbstractString} || b isa Union{Symbol,AbstractString}
    return a != b
  elseif (a isa AbstractVector && b isa AbstractVector) || (a isa Tuple && b isa Tuple)
    return length(a) != length(b) || any(i -> _tap_values_differ(a[i], b[i], differs), eachindex(a))
  elseif typeof(a) === typeof(b) && isstructtype(typeof(a))
    return any(f -> _tap_values_differ(getfield(a, f), getfield(b, f), differs), fieldnames(typeof(a)))
  end
  return a != b
end

# number of the tap data record `rec` of transformer `id` in the index; an
# equal record is reused, so a series whose steps alternate between two
# tables holds two records
function _tap_record!(index::Sparlectra.ScenarioIndex, id::Int, rec::Sparlectra.BranchTapData)::Int
  recs = get!(Vector{Sparlectra.BranchTapData}, index.tap_data, id)
  k = findfirst(r -> !_tap_values_differ(r, rec, (x, y) -> !isequal(x, y)), recs)
  k === nothing || return k
  push!(recs, rec)
  return length(recs)
end

# Lookup tables for the transformers of the tap controllers of T and R_k:
# transformer name -> branch position (the names _find_trafo_branch accepts,
# the first branch wins as there) and branch position -> mRID (the template
# address with the alias rule of the delta, the R_k mRID)
function _tapControllerMaps(T::Sparlectra.Net, Rk::CGMESImportResult, branch_addr::Dict{Int,String}, alias::Dict{String,String})
  function positions(net)
    pos = Dict{String,Int}()
    for (j, br) in enumerate(net.branchVec)
      br.ratio != 0.0 || continue
      for key in (br.comp.cName, br.comp.cID, string(br.branchIdx))
        get!(pos, string(key), j)
      end
    end
    return pos
  end
  t_addr = Dict{Int,String}(j => get(alias, a, a) for (j, a) in branch_addr)
  r_mrid = Dict{Int,String}()
  for (m, (kind, j)) in Rk.element_by_mrid
    (kind === :branch || kind === :transformer) && (r_mrid[j] = m)
  end
  return positions(T), positions(Rk.net), t_addr, r_mrid
end

# tap model fields the CGMES import does not set; the delta has no op for
# them, a difference is residual
const _DELTA_TAP_MODEL_FIELDS = (:taps_derived, :tap_correction, :phase_du_step, :phase_du_min_step, :phase_du_max_step)

# union bus -> union node id, and R_k bus of a TN
function _bus_node_map(res::CGMESImportResult, node_id::Dict{String,Int})::Dict{Int,Int}
  out = Dict{Int,Int}()
  for (tn, name) in res.topo.bus_name
    id = get(node_id, tn, 0)
    id == 0 && continue
    bus = get(res.net.busDict, name, 0)
    bus == 0 || (out[bus] = id)
  end
  return out
end

function _tn_bus_map(res::CGMESImportResult)::Dict{String,Int}
  out = Dict{String,Int}()
  for (tn, name) in res.topo.bus_name
    bus = get(res.net.busDict, name, 0)
    bus == 0 || (out[tn] = bus)
  end
  return out
end

# link state of step k by link mRID: switch links and _TPL couplings carry
# it from the partition, not from R_k (R_k has links only for held switches)
function _link_status_of_step(part::_UnionPartition, links::_UnionLinks, k::Int)::Dict{String,Int}
  out = Dict{String,Int}()
  for (j, m) in enumerate(links.switches)
    out[m] = links.switch_status[j][k]
  end
  for (j, m) in enumerate(_tpl_mrids(part, links))
    out[m] = links.tpl_status[j][k]
  end
  return out
end

# Delta of step k: the net `T` (the template, or a patched copy of it, with
# the units in `parked_now` parked) against the direct import `Rk` of that
# snapshot, addressed by mRID through `index`. `node_of_bus` maps a union
# bus to its union node, `node_sig[node][k]` is that node's TN in step k.
function _snapshotDelta(T::Sparlectra.Net, index::Sparlectra.ScenarioIndex, parked_now::Set{Int}, node_of_bus::Dict{Int,Int}, node_sig::Vector{Vector{String}},
                        link_status::Dict{String,Int}, alias::Dict{String,String}, Rk::CGMESImportResult, k::Int, label::String, scenario_time;
                        ssh_steps::Dict{String,Float64} = Dict{String,Float64}())::CGMESSnapshotDelta
  PatchOp = Sparlectra.PatchOp
  ops = PatchOp[]; control_ops = PatchOp[]; start_ops = PatchOp[]; residual = String[]
  rk_bus_of_tn = _tn_bus_map(Rk)
  # the R_k bus a union bus lands on in this step (0 = not built there)
  rk_bus(ubus) = (nid = get(node_of_bus, ubus, 0); nid == 0 ? 0 : get(rk_bus_of_tn, node_sig[nid][k], 0))
  # a bus without a TN (the star point of a three-winding transformer) has
  # no union node; it is paired with the R_k bus of the first element seen
  # on it, and every later element must agree with that pairing
  aux = Dict{Int,Int}()
  placed(ubus, rbus) = haskey(node_of_bus, ubus) ? rk_bus(ubus) == rbus : get!(aux, ubus, rbus) == rbus
  set!(target, oid, field, value) = push!(ops, PatchOp(op = :set, target = target, id = oid, field = field, value = value === nothing ? nothing : Float64(value)))
  # ratio and shift of a transformer whose regulator is active in R_k (a
  # group follower included) go to control_ops: a snapshot-wise run applies
  # them, a sequential run leaves the position to the controller
  setc!(oid, field, value) = push!(control_ops, PatchOp(op = :set, target = :transformer, id = oid, field = field, value = Float64(value)))
  rk_controlled = Set{Int}()
  for c in Sparlectra._tap_controllers(Rk.net)
    c.enabled || continue
    for name in vcat([c.trafo], c.followers)
      br = Sparlectra._find_trafo_branch(Rk.net, name)
      j = findfirst(x -> x === br, Rk.net.branchVec)
      j === nothing || push!(rk_controlled, j)
    end
  end
  # R_k buses whose start voltage is part of the state (see the bus ops
  # below): every bus of an electrical island without a slack, and every
  # Slack bus, whose start angle is the angle reference of its island (an
  # island the importer gave a slack keeps the SV angle of that bus, not 0)
  start_fixed = Set{Int}()
  # islands of the power flow: a part behind a switch-only bus shares the
  # slack of the rest and is not an island without slack of its own
  islands, _, refs = Sparlectra._power_flow_islands(Rk.net)
  for (c, slacks) in zip(islands, refs)
    union!(start_fixed, isempty(slacks) ? c : slacks)
  end
  # a bus name of the template and one of R_k on the same bus of this step
  same_bus(a, b) = (a === nothing) == (b === nothing) &&
                   (a === nothing || (get(T.busDict, a, 0) > 0 && get(Rk.net.busDict, b, 0) > 0 && placed(T.busDict[a], Rk.net.busDict[b])))
  # the transformers a tap controller moves (its own and the followers) as
  # sorted mRIDs; `pos` maps a transformer name to its branch position,
  # `mrid_of` a position to its mRID (built on the first tap controller)
  moved(pos, mrid_of, c) = sort!([get(mrid_of, get(pos, n, 0), string("?", n)) for n in vcat([c.trafo], c.followers)])
  tc_maps = nothing
  # branch position -> address, for the legs of a star point
  branch_addr = Dict{Int,String}()
  for (a, bid) in index.id_by_mrid
    index.kind_by_id[bid] in (:branch, :transformer) && (branch_addr[index.internal_by_id[bid]] = a)
  end
  for (addr, id) in sort!(collect(index.id_by_mrid); by = last)
    kind = index.kind_by_id[id]; i = index.internal_by_id[id]
    if kind === :bus
      # isolation of the union bus is the isolation of its R_k bus; a bus
      # R_k did not build is isolated. A star point has no TN: its R_k bus is
      # the star end of a leg R_k keeps in service, and without one R_k
      # builds no star point at all
      rb = 0
      if endswith(addr, "#star")
        for (j, br) in enumerate(T.branchVec)
          (br.fromBus == i || br.toBus == i) || continue
          r = get(Rk.element_by_mrid, get(branch_addr, j, ""), nothing)
          r === nothing && continue
          rbr = Rk.net.branchVec[r[2]]
          rbr.status == 1 || continue
          rb = br.fromBus == i ? rbr.fromBus : rbr.toBus
          break
        end
      else
        rb = rk_bus(i)
      end
      iso_r = rb == 0 || Sparlectra.isIsolated(Rk.net.nodeVec[rb])
      iso_t = Sparlectra.isIsolated(T.nodeVec[i])
      iso_r != iso_t && set!(:bus, id, :isolated, iso_r ? 1 : 0)
      # the start of the R_k bus goes along where the bus comes back, and
      # where it fixes the solution: in an island without a slack the solver
      # promotes a unit to reference and keeps its start voltage (measured),
      # and a slack bus holds its start angle
      iso_r && continue
      nd = Rk.net.nodeVec[rb]
      vm = something(nd._vm_pu, 1.0); va = something(nd._va_deg, 0.0)
      if iso_r != iso_t || rb in start_fixed
        tn = T.nodeVec[i]
        (iso_r != iso_t || _delta_differs(tn._vm_pu, vm)) && set!(:bus, id, :vm_start, vm)
        (iso_r != iso_t || _delta_differs(tn._va_deg, va)) && set!(:bus, id, :va_start, va)
      end
      # mode A starts every step from the state of its own snapshot: the start
      # of every live bus, written whether or not it differs from the
      # template, since an op of the step (a generator vm_pu, an angle) may
      # have moved the bus voltage of the template
      push!(start_ops, PatchOp(op = :set, target = :bus, id = id, field = :vm_start, value = vm))
      push!(start_ops, PatchOp(op = :set, target = :bus, id = id, field = :va_start, value = va))
      continue
    end
    # an element of R_k counts only with the same kind: a reused mRID of
    # another class is a different element (see `alias`)
    r = get(Rk.element_by_mrid, get(alias, addr, addr), nothing)
    (r === nothing || r[1] === kind) || (r = nothing)
    if kind === :branch || kind === :transformer
      tb = T.branchVec[i]
      if r === nothing
        for (f, v) in ((:status, tb.status), (:from_status, tb.from_status), (:to_status, tb.to_status))
          v == 0 || set!(kind, id, f, 0)
        end
        continue
      end
      rb = Rk.net.branchVec[r[2]]
      for (f, a, b) in ((:status, tb.status, rb.status), (:from_status, tb.from_status, rb.from_status), (:to_status, tb.to_status, rb.to_status))
        a == b || set!(kind, id, f, b)
      end
      if kind === :transformer
        ctrl = r[2] in rk_controlled
        _delta_differs(tb.ratio, rb.ratio) && (ctrl ? setc!(id, :ratio, rb.ratio) : set!(:transformer, id, :ratio, rb.ratio))
        # shift_deg writes angle and phase_shift_deg together; R_k with two
        # different values cannot be expressed
        if _delta_differs(tb.phase_shift_deg, rb.phase_shift_deg) || _delta_differs(tb.angle, rb.angle)
          ctrl ? setc!(id, :shift_deg, rb.phase_shift_deg) : set!(:transformer, id, :shift_deg, rb.phase_shift_deg)
          _delta_differs(rb.angle, rb.phase_shift_deg) && push!(residual, string(kind, " ", addr, ": angle ", rb.angle, " differs from phase_shift_deg ", rb.phase_shift_deg))
        end
      elseif tb.ratio != 0.0 && (_delta_differs(tb.phase_shift_deg, rb.phase_shift_deg) || _delta_differs(tb.angle, rb.angle))
        # a phase-stamped EquivalentBranch is a :branch with ratio 1 whose
        # shift follows r21, x21 of the snapshot (_ebPhaseStamp); where R_k
        # builds it as a line (ratio 0, no shift) the shift goes to 0, the
        # same two-port as that line
        set!(:branch, id, :shift_deg, rb.phase_shift_deg)
      elseif tb.ratio == 0.0 && rb.ratio != 0.0
        # a line in the net, a transformer branch in R_k: apply! has no op
        # that turns a line into a transformer branch (_stampableBranches!
        # prepares the template for a stamped EquivalentBranch)
        push!(residual, string(kind, " ", addr, ": a line in the union, a transformer branch with ratio ", rb.ratio, " and shift ", rb.phase_shift_deg, " in the snapshot"))
      end
      # tap changer data: range, step, table ratios and impedance per step.
      # The importer sets them for every changer, the range of a regulated
      # one from its controller, so a step whose EQ changes the table or
      # whose regulated changer is another one of the branch has other data
      # than the template.
      # One record per distinct data set; a line in R_k is the two-port of
      # a stamped template branch (see above) and has no tap data to compare
      if tb.ratio != 0.0 && rb.ratio != 0.0
        rec = Sparlectra.BranchTapData(rb)
        if _tap_values_differ(Sparlectra.BranchTapData(tb), rec, _delta_differs)
          kind === :transformer ? set!(:transformer, id, :tap_data, _tap_record!(index, id, rec)) :
            push!(residual, string(kind, " ", addr, ": tap changer data of a branch that is no transformer"))
        end
      end
      for f in _DELTA_TAP_MODEL_FIELDS
        _tap_values_differ(getfield(tb, f), getfield(rb, f), _delta_differs) && push!(residual, string(kind, " ", addr, ": ", f, " ", getfield(tb, f), " -> ", getfield(rb, f)))
      end
      (tb.tap_winding === nothing) == (rb.tap_winding === nothing) || push!(residual, string(kind, " ", addr, ": typed tap model only on one side"))
      # parameters of the equipment (an export may change them between
      # snapshots, measured on a real series); r and x ops write the live and
      # the physical base value together, a base that differs from the live
      # value in R_k (a FACTS operating point) cannot be expressed
      for f in (:r_pu, :x_pu, :b_from_pu, :g_from_pu, :b_to_pu, :g_to_pu)
        _delta_differs(getfield(tb, f), getfield(rb, f)) && set!(kind, id, f, getfield(rb, f))
      end
      for (f, live) in ((:r_base_pu, :r_pu), (:x_base_pu, :x_pu))
        _delta_differs(getfield(rb, f), getfield(rb, live)) && push!(residual, string(kind, " ", addr, ": ", f, " ", getfield(rb, f), " differs from ", live, " ", getfield(rb, live)))
      end
      _delta_differs(tb.sn_MVA, rb.sn_MVA) && push!(residual, string(kind, " ", addr, ": sn_MVA ", tb.sn_MVA, " -> ", rb.sn_MVA))
      # placement: both ends of the R_k branch on the R_k buses of its union buses
      (placed(tb.fromBus, rb.fromBus) && placed(tb.toBus, rb.toBus)) ||
        push!(residual, string(kind, " ", addr, ": placement (union buses ", tb.fromBus, "/", tb.toBus, ", R_k buses ", rb.fromBus, "/", rb.toBus, ")"))
    elseif kind in (:generator, :load, :external_grid)
      tp = T.prosumpsVec[i]
      parked = id in parked_now
      if r === nothing
        parked || set!(kind, id, :in_service, 0)
        continue
      end
      rp = Rk.net.prosumpsVec[r[2]]
      parked && set!(kind, id, :in_service, 1)
      ref = parked ? index.parked_values[id] : tp
      _delta_differs(ref.pVal, rp.pVal) && set!(kind, id, :p, something(rp.pVal, 0.0))
      _delta_differs(ref.qVal, rp.qVal) && set!(kind, id, :q, something(rp.qVal, 0.0))
      if kind !== :load
        # vm_pu is SCF vocabulary and needs a number; a setpoint that goes
        # away is the :regulated op, a remaining difference is residual
        if _delta_differs(ref.vm_pu, rp.vm_pu)
          rp.vm_pu === nothing ? push!(residual, string(kind, " ", addr, ": voltage setpoint removed")) : set!(kind, id, :vm_pu, rp.vm_pu)
        end
        _delta_differs(ref.minQ, rp.minQ) && set!(kind, id, :qmin, _delta_value(rp.minQ))
        _delta_differs(ref.maxQ, rp.maxQ) && set!(kind, id, :qmax, _delta_value(rp.maxQ))
        ref.isRegulated == rp.isRegulated || set!(kind, id, :regulated, rp.isRegulated ? 1 : 0)
        (ref.referencePri === nothing) == (rp.referencePri === nothing) || set!(kind, id, :reference, rp.referencePri === nothing ? 0 : 1)
        ref.referencePriority == rp.referencePriority || set!(kind, id, :reference_priority, rp.referencePriority)
        _delta_differs(ref.participationFactor, rp.participationFactor) && set!(kind, id, :participation, _delta_value(rp.participationFactor))
      end
      ref.proSumptionType == rp.proSumptionType || push!(residual, string(kind, " ", addr, ": type ", ref.proSumptionType, " -> ", rp.proSumptionType))
      placed(Sparlectra._scf_prosumer_bus(tp), Sparlectra._scf_prosumer_bus(rp)) || push!(residual, string(kind, " ", addr, ": placement"))
    elseif kind === :shunt
      ts = T.shuntVec[i]
      if r === nothing
        ts.status == 0 || push!(ops, PatchOp(op = :status, target = :shunt, id = id, value = 0.0))
        continue
      end
      rs = Rk.net.shuntVec[r[2]]
      ts.status == rs.status || push!(ops, PatchOp(op = :status, target = :shunt, id = id, value = Float64(rs.status)))
      _delta_differs(imag(ts.y_pu_shunt), imag(rs.y_pu_shunt)) && set!(:shunt, id, :b_pu, imag(rs.y_pu_shunt))
      _delta_differs(real(ts.y_pu_shunt), real(rs.y_pu_shunt)) && set!(:shunt, id, :g_pu, real(rs.y_pu_shunt))
    elseif kind === :link
      # switch links and _TPL couplings carry their state per step from the
      # partition, not from R_k: R_k has links only for held switches
      st = get(link_status, addr, nothing)
      if st === nothing
        push!(residual, string("link ", addr, ": not a link of the union partition"))
      else
        T.linkVec[i].status == st || push!(ops, PatchOp(op = :status, target = :link, id = id, value = Float64(st)))
      end
    elseif kind === :tap_controller
      tc = Sparlectra._tap_controllers(T)[i]
      if r === nothing
        tc.enabled && set!(:tap_controller, id, :enabled, 0)
        continue
      end
      rc = Sparlectra._tap_controllers(Rk.net)[r[2]]
      tc.enabled == rc.enabled || set!(:tap_controller, id, :enabled, rc.enabled ? 1 : 0)
      _delta_differs(tc.target_vm_pu, rc.target_vm_pu) && set!(:tap_controller, id, :target_vm_pu, _delta_value(rc.target_vm_pu))
      _delta_differs(tc.deadband_vm_pu, rc.deadband_vm_pu) && set!(:tap_controller, id, :deadband_vm_pu, rc.deadband_vm_pu)
      _delta_differs(tc.p_target_mw, rc.p_target_mw) && set!(:tap_controller, id, :p_target_mw, _delta_value(rc.p_target_mw))
      _delta_differs(tc.deadband_p_mw, rc.deadband_p_mw) && set!(:tap_controller, id, :deadband_p_mw, rc.deadband_p_mw)
      # what the vocabulary has no op for: the kind of control, its ranks
      # and range, the regulated bus or branch, and the transformers it moves
      for f in (:mode, :control_ratio, :control_phase, :is_discrete, :voltage_error_metric, :max_outer_iters, :priority, :tap_min, :tap_max, :q_target_mvar)
        _tap_values_differ(getfield(tc, f), getfield(rc, f), _delta_differs) && push!(residual, string(kind, " ", addr, ": ", f, " ", getfield(tc, f), " -> ", getfield(rc, f)))
      end
      same_bus(tc.target_bus, rc.target_bus) || push!(residual, string(kind, " ", addr, ": regulated bus"))
      ((tc.target_branch === nothing) == (rc.target_branch === nothing) &&
       (tc.target_branch === nothing || (same_bus(tc.target_branch[1], rc.target_branch[1]) && same_bus(tc.target_branch[2], rc.target_branch[2])))) ||
        push!(residual, string(kind, " ", addr, ": regulated branch"))
      tc_maps === nothing && (tc_maps = _tapControllerMaps(T, Rk, branch_addr, alias))
      t_pos, r_pos, t_addr, r_mrid = tc_maps
      moved(t_pos, t_addr, tc) == moved(r_pos, r_mrid, rc) || push!(residual, string(kind, " ", addr, ": transformer or group members"))
    end
  end
  # R_k elements the union does not know
  known = Set{Tuple{String,Symbol}}((get(alias, a, a), index.kind_by_id[i]) for (a, i) in index.id_by_mrid)
  for (m, (kind, _)) in Rk.element_by_mrid
    (m, kind) in known || push!(residual, string(kind, " ", m, ": only in the snapshot"))
  end
  return CGMESSnapshotDelta(k, label, scenario_time, ops, control_ops, start_ops, ssh_steps, residual)
end

_tpl_mrids(part::_UnionPartition, links::_UnionLinks) = [string("_TPL", lstrip(part.node_mrid[a], '_'), "_", lstrip(part.node_mrid[c], '_')) for (a, c) in links.tpl_nodes]

# --- the union of a series ---------------------------------------------------

"""
    CGMESSeriesUnion

The union network of a CGMES snapshot series and one delta per step.
`template` is the union net in the state of the base step; applying
`deltas[k].ops` (plus `control_ops` for a snapshot-wise run) with `index`
gives step `k`. Node identity: `node_of_terminal` maps every terminal mRID
to its union node, `node_mrid[id]` is the union TN behind it, `tn_buses[k]`
maps a TN of step `k` to the union buses it covers. `link_switches` are the
switches that became links, `n_tp_links` the extra `_TPL` couplings a step
needs where two union nodes share its TN without a switch between them.
`n_switch_link_changes[k]` and `n_tpl_link_changes[k]` count the link
status ops of step `k` against the template by kind: only the first are
switch operations; a `_TPL` coupling changes where an export puts terminals
on other topological nodes (equivalents on other connectivity nodes, for
example). Their sum is the number of link ops in `deltas[k].ops`.
`alias` maps the address `<mRID>#<class>` of an object whose mRID an
exporter reused for another class back to that mRID. `sv[k]` is the SV
state of snapshot `k` ([`QstsSVReference`](@ref Sparlectra.QstsSVReference)), the reference of the SV
comparison of [`runQsts`](@ref Sparlectra.runQsts); `import_config` is the `cgmes_import`
configuration the union was built with, which a run must use as well.
"""
struct CGMESSeriesUnion
  template::Sparlectra.Net
  index::Sparlectra.ScenarioIndex
  base_step::Int
  labels::Vector{String}
  scenario_times::Vector{Union{Nothing,DateTime}}
  deltas::Vector{CGMESSnapshotDelta}
  parked::Set{Int}
  node_of_terminal::Dict{String,Int}
  node_mrid::Vector{String}
  node_sig::Vector{Vector{String}}
  node_of_bus::Dict{Int,Int}
  link_status::Vector{Dict{String,Int}}
  alias::Dict{String,String}
  tn_buses::Vector{Dict{String,Vector{Int}}}
  link_switches::Vector{String}
  n_tp_links::Int
  messages::Vector{String}
  sv::Vector{Sparlectra.QstsSVReference}
  import_config::Sparlectra.CGMESImportConfig
  n_switch_link_changes::Vector{Int}
  n_tpl_link_changes::Vector{Int}
end

# id map of the union: every addressed element plus one :bus per union bus;
# ids are dense 1..N over the sorted addresses
function _unionIndex(U::CGMESImportResult)::Sparlectra.ScenarioIndex
  addr = Dict{String,Tuple{Symbol,Int}}(U.element_by_mrid)
  tn_buses = Set{Int}()
  for (tn, name) in U.topo.bus_name
    bus = get(U.net.busDict, name, 0)
    (bus == 0 || haskey(addr, tn)) && continue
    addr[tn] = (:bus, bus)
    push!(tn_buses, bus)
  end
  # the star point of a three-winding transformer has no TN; it is addressed
  # as `<PowerTransformer mRID>#star`, through the end of one of its legs
  for (m, (kind, i)) in U.element_by_mrid
    kind === :transformer || continue
    br = U.net.branchVec[i]
    for b in (br.fromBus, br.toBus)
      b in tn_buses && continue
      o = get(U.store.objects, m, nothing)
      pt = o === nothing ? m : get(o.refs, :PowerTransformer, m)
      addr[string(pt, "#star")] = (:bus, b)
    end
  end
  kinds = Dict{Int,Symbol}(); internals = Dict{Int,Int}(); names = Dict{Int,String}(); ids = Dict{String,Int}()
  for (id, m) in enumerate(sort!(collect(keys(addr))))
    kind, internal = addr[m]
    kinds[id] = kind; internals[id] = internal; ids[m] = id
    o = get(U.store.objects, m, nothing)
    nm = o === nothing ? nothing : str(o, :name)
    nm === nothing || (names[id] = nm)
  end
  return Sparlectra.ScenarioIndex(kinds, internals, names, ids, Dict{Int,Sparlectra.ProSumer}(), :cgmes_series)
end

_series_store(src::CGMESStore) = src
_series_store(src) = redirect_stdout(() -> loadCGMES(src), devnull)

# scenario time of a step: the input's, else the SSH header's, else the EQ
# header's (the rule of checkCGMESSeries)
function _series_time(x::CGMESSeriesInput, s::CGMESStore)::Union{Nothing,DateTime}
  x.scenario_time === nothing || return x.scenario_time
  for prof in (:SSH, :EQ), f in s.files
    prof in f.profiles || continue
    t = _cgmesParseTime(f.scenario_time)
    t === nothing || return t
  end
  return nothing
end

# Pass 1, one snapshot store at a time: each step leaves its scenario time,
# its topology view, its switches and the objects the union store does not
# have yet, then its store goes. Only the base step's store is held until
# the union store is built; it is loaded first because every other step is
# compared against it. Returns the union store, the alias map, the times,
# the partition and the links; no snapshot store outlives this call.
# TapChanger steps and ShuntCompensator sections of one SSH, by mRID: the
# reference positions a sequential run compares its controllers with
function _sshSteps(s::CGMESStore)::Dict{String,Float64}
  out = Dict{String,Float64}()
  for cls in (:RatioTapChanger, :PhaseTapChangerLinear, :PhaseTapChangerSymmetrical, :PhaseTapChangerAsymmetrical, :PhaseTapChangerTabular), o in objectsOf(s, cls)
    v = num(o, :step)
    v === nothing || (out[o.mrid] = v)
  end
  for cls in (:LinearShuntCompensator, :NonlinearShuntCompensator), o in objectsOf(s, cls)
    v = num(o, :sections)
    v === nothing || (out[o.mrid] = v)
  end
  return out
end

function _seriesUnionStore(inputs::AbstractVector{CGMESSeriesInput}, labels::Vector{String}, base::Int, extra::Set{String})
  n = length(inputs)
  b = _series_store(inputs[base].source)
  base_file = _unionBaseFiles(b)
  ex = _UnionExtras()
  times = Vector{Union{Nothing,DateTime}}(undef, n)
  topos = Vector{_SeriesTopo}(undef, n)
  ssh_steps = Vector{Dict{String,Float64}}(undef, n)
  # switch mRID -> class; a later step overwrites an earlier one
  switch_class = Dict{String,Symbol}()
  for k in 1:n
    s = k == base ? b : _series_store(inputs[k].source)
    times[k] = _series_time(inputs[k], s)
    topos[k] = _seriesTopo(s, labels[k])
    ssh_steps[k] = _sshSteps(s)
    for cls in _SWITCH_CLASSES, o in objectsOf(s, cls)
      switch_class[o.mrid] = cls
    end
    k == base || _collectExtras!(ex, b, base_file, s)
  end
  part = _splitPartition(_unionPartition(topos, base), topos, switch_class, extra, b.boundary)
  links = _unionLinks(topos, part, switch_class, extra)
  ustore, alias = _unionStore(b, ex, topos, part, links, base)
  return ustore, alias, times, part, links, ssh_steps
end

# The union import after pass 1. Returns the union net (it becomes the
# template), its index and the union bus -> node map with the outputs of
# pass 1; the union store and the import result end with this call.
function _seriesUnionNet(inputs::AbstractVector{CGMESSeriesInput}, labels::Vector{String}, config, base::Int, extra::Set{String})
  ustore, alias, times, part, links, ssh_steps = _seriesUnionStore(inputs, labels, base, extra)
  U = importCGMES(config; path = ustore, name = string("union of ", labels[base]), ignore_connected = true)
  index = _unionIndex(U)
  node_of_bus = _bus_node_map(U, Dict(m => i for (i, m) in enumerate(part.node_mrid)))
  return U.net, index, node_of_bus, alias, times, part, links, ssh_steps
end

# The SV state of the direct import `Rk` of one snapshot as the reference of
# the snapshot run, with the exclusions of compareWithSV: voltages of every
# TN with an SvVoltage, a bus, an R_k bus that is not isolated or that a
# closed switch joins to a bus with a branch (`_linkFedIsolatedBuses`, the
# step gives it the values of its switch group), `v` and `vn` above 0,
# together with the union buses of that TN (`tn_buses`, a TN without union
# bus is left out); flows of every SvPowerFlow at a line or transformer
# terminal whose equipment the import mapped, whose R_k branch is in service
# and whose terminal bus is not isolated. The R_k branch goes to its template
# position over its mRID with the kind and alias rule of the delta (`id_of`:
# (mRID, kind) -> union id); the side is the one R_k gives the terminal.
function _qstsSVReference(Rk::CGMESImportResult, tn_buses::Dict{String,Vector{Int}}, id_of::Dict{Tuple{String,Symbol},Int}, index::Sparlectra.ScenarioIndex)::Sparlectra.QstsSVReference
  net = Rk.net
  topo = Rk.topo
  iso = Set(net.isoNodes)
  fed = Sparlectra._linkFedIsolatedBuses(net)
  tnb = Vector{Int}[]
  vm = Float64[]
  va = Float64[]
  for sv in objectsOf(Rk.store, :SvVoltage)
    tn = get(sv.refs, :TopologicalNode, nothing)
    tn === nothing && continue
    name = get(topo.bus_name, tn, nothing)
    name === nothing && continue
    bus = get(net.busDict, name, 0)
    (bus == 0 || (bus in iso && !haskey(fed, bus))) && continue
    vn = topo.vn_kV[tn]
    v = num(sv, :v)
    (v === nothing || vn <= 0.0 || v <= 0.0) && continue
    ub = get(tn_buses, tn, nothing)
    ub === nothing && continue
    push!(tnb, sort(ub))
    push!(vm, v / vn)
    push!(va, num(sv, :angle, 0.0))
  end
  position_of = Dict{Int,Int}(br.branchIdx => j for (j, br) in enumerate(net.branchVec))
  element_at = Dict{Int,Tuple{String,Symbol}}()
  for (m, (kind, j)) in Rk.element_by_mrid
    (kind === :branch || kind === :transformer) && (element_at[j] = (m, kind))
  end
  branch = Int[]
  side = Symbol[]
  p_mw = Float64[]
  q_mvar = Float64[]
  for pf in objectsOf(Rk.store, :SvPowerFlow)
    tm = get(pf.refs, :Terminal, nothing)
    tm === nothing && continue
    bs = get(Rk.branch_side_of_terminal, tm, nothing)
    bs === nothing && continue
    t = get(Rk.store.objects, tm, nothing)
    t === nothing && continue
    eqm = get(t.refs, :ConductingEquipment, nothing)
    (eqm === nothing || eqm in Rk.skipped_equipment) && continue
    p = num(pf, :p)
    q = num(pf, :q)
    (p === nothing || q === nothing) && continue
    j = get(position_of, bs[1], 0)
    j == 0 && continue
    rbr = net.branchVec[j]
    rbr.status == 0 && continue
    (bs[2] === :from ? rbr.fromBus : rbr.toBus) in iso && continue
    el = get(element_at, j, nothing)
    el === nothing && continue
    id = get(id_of, el, 0)
    id == 0 && continue
    push!(branch, index.internal_by_id[id])
    push!(side, bs[2])
    push!(p_mw, p)
    push!(q_mvar, q)
  end
  return Sparlectra.QstsSVReference(tnb, vm, va, branch, side, p_mw, q_mvar)
end

# A same-level EquivalentBranch is a line (ratio 0) where r21, x21 equal r,
# x and a transformer branch with ratio 1 and the phase stamp where they
# differ (_mapEquivalentBranches!). An export may change r21, x21 between
# snapshots (measured on a real series: 992 branches between two days), so a
# branch that is a line in the template can be stamped in step k. apply! has
# no op that turns a line into a transformer branch, and the shift of a
# line cannot be set. Such a template branch becomes a transformer branch
# with ratio 1 and shift 0 here, the two-port of the line, built as the
# importer builds a transformer branch; the delta then writes the shift of
# every step, 0 where the step has the line. Called before each delta
# (base step included): a delta computed before the change stays valid,
# since ratio 1 with shift 0 is the line and the delta of a :branch does
# not compare the ratio.
function _stampableBranches!(T::Sparlectra.Net, index::Sparlectra.ScenarioIndex, alias::Dict{String,String}, Rk::CGMESImportResult)
  for (addr, id) in index.id_by_mrid
    index.kind_by_id[id] === :branch || continue
    tb = T.branchVec[index.internal_by_id[id]]
    tb.ratio == 0.0 || continue
    r = get(Rk.element_by_mrid, get(alias, addr, addr), nothing)
    (r === nothing || r[1] !== :branch) && continue
    # only the ratio-1 branch of an EquivalentBranch; anything else is
    # residual of the delta
    Rk.net.branchVec[r[2]].ratio == 1.0 || continue
    tb.ratio = 1.0
    tb.tap_ratio = 1.0
    tb.has_ratio_tap = true
    tb.has_phase_tap = true
  end
  return T
end

# Turn the union net `T` in place into the template: the union in the state
# of the base step `k`. Units the base step lacks are parked, their union
# values go to `index.parked_values`. Returns the parked set, the residual
# of the base step against the bare union, the delta of the base step
# against the finished template and the SV reference of the base step
# (`tn_buses`, `id_of` as for _qstsSVReference); the direct import of the
# base step ends with this call.
function _seriesTemplate!(T::Sparlectra.Net, index::Sparlectra.ScenarioIndex, node_of_bus::Dict{Int,Int}, node_sig::Vector{Vector{String}}, link_status::Dict{String,Int},
                          alias::Dict{String,String}, x::CGMESSeriesInput, label::String, scenario_time, config, k::Int; ssh_steps::Dict{String,Float64} = Dict{String,Float64}(),
                          tn_buses::Dict{String,Vector{Int}} = Dict{String,Vector{Int}}(), id_of::Dict{Tuple{String,Symbol},Int} = Dict{Tuple{String,Symbol},Int}())
  Rb = importCGMES(config; path = x.source, name = label)
  _stampableBranches!(T, index, alias, Rb)
  d0 = _snapshotDelta(T, index, Set{Int}(), node_of_bus, node_sig, link_status, alias, Rb, k, label, scenario_time)
  parked = Set{Int}()
  for op in d0.ops
    if op.op === :set && op.field === :in_service && op.value == 0.0
      index.parked_values[op.id] = deepcopy(T.prosumpsVec[index.internal_by_id[op.id]])
      push!(parked, op.id)
    end
  end
  Sparlectra.apply!(T, vcat(d0.ops, d0.control_ops), index)
  delta = _snapshotDelta(T, index, parked, node_of_bus, node_sig, link_status, alias, Rb, k, label, scenario_time; ssh_steps = ssh_steps)
  return parked, d0.residual, delta, _qstsSVReference(Rb, tn_buses, id_of, index)
end

"""
    buildCGMESSeriesUnion(inputs; config, base_step = config.qsts.base_step, extra_link_switches = String[]) -> CGMESSeriesUnion

Build the union network of one snapshot series (`inputs` as for
[`checkCGMESSeries`](@ref), in time order) and the delta of every step
against it. The template is built on the snapshot of step `base_step`
(default `qsts.base_step` of `config`, 1 = the first snapshot); a step
outside the series throws `CGMESSeriesError`. The unchanged importer maps the union store once with
`ignore_connected = true`, so every element exists; each snapshot is also
imported directly with `config`, and the delta is the comparison of both
nets by mRID; the same import gives the SV state of the step (`sv`, the
reference of the SV comparison of a snapshot run), and the union keeps
`config.cgmes` as `import_config`. Units absent from the base step are
parked in the template.
Node identity uses terminals (see the design notes in the source).

The snapshots are read one at a time in two passes: pass 1 takes from each
store what the union store needs, pass 2 imports each snapshot directly
for its delta. A path source is loaded in each pass, a `CGMESStore` source
is used as given and not modified. Next to the union and one topology view
per step, about two snapshot stores are held at a time.

`extra_link_switches` names switches (mRIDs) that become links although
the TP keeps both their sides on one topological node in every step, so an
event can open them. The union node is split at each such switch along its
ConnectivityNodes and the other switches of the node; the switch is then a
link with state 1 wherever the TP joins its sides. A named switch that
already lies between two union nodes is a link even where it conducts in no
step. Default: none, the union is unchanged.

Throws `CGMESSeriesError` for `cgmes_import.machine_control = true` and
`cgmes_import.hvdc_mode = :paired_control`: their controllers bind to prosumer
positions that a series would have to rebind per step. It also throws for an
entry of `extra_link_switches` that is not a switch of the series, sits in a
boundary node or a node without node-breaker connectivity, or whose two
sides another switch joins (name that switch as well).
"""
function buildCGMESSeriesUnion(inputs::AbstractVector{CGMESSeriesInput}; config, base_step::Int = config.qsts.base_step,
                               extra_link_switches::AbstractVector{<:AbstractString} = String[])::CGMESSeriesUnion
  n = length(inputs)
  n >= 1 || throw(CGMESSeriesError("empty series"))
  1 <= base_step <= n || throw(CGMESSeriesError("base_step $(base_step) outside 1:$(n)"))
  config.cgmes.machine_control && throw(CGMESSeriesError("cgmes_import.machine_control = true is not supported for a series: machine controllers bind to prosumer positions"))
  config.cgmes.hvdc_mode === :paired_control && throw(CGMESSeriesError("cgmes_import.hvdc_mode = paired_control is not supported for a series: the HVDC pairs bind to prosumer positions"))
  labels = [x.label for x in inputs]
  T, index, node_of_bus, alias, times, part, links, ssh_steps = _seriesUnionNet(inputs, labels, config, base_step, Set{String}(String(m) for m in extra_link_switches))
  link_status = [_link_status_of_step(part, links, k) for k in 1:n]
  # per step: TN of that step -> the union buses it covers
  tn_buses = [Dict{String,Vector{Int}}() for _ in 1:n]
  bus_of_node = Dict{Int,Int}((nid, b) for (b, nid) in node_of_bus)
  for (nid, sig) in enumerate(part.node_sig), k in 1:n
    b = get(bus_of_node, nid, 0)
    (b == 0 || sig[k] == "") && continue
    push!(get!(Vector{Int}, tn_buses[k], sig[k]), b)
  end
  # an R_k element (mRID, kind) -> its union id, with the alias rule of the delta
  id_of = Dict{Tuple{String,Symbol},Int}((get(alias, a, a), index.kind_by_id[id]) => id for (a, id) in index.id_by_mrid)
  # template: the union in the state of the base step
  deltas = Vector{CGMESSnapshotDelta}(undef, n)
  sv = Vector{Sparlectra.QstsSVReference}(undef, n)
  parked, base_residual, deltas[base_step], sv[base_step] = _seriesTemplate!(T, index, node_of_bus, part.node_sig, link_status[base_step], alias, inputs[base_step], labels[base_step],
                                                                             times[base_step], config, base_step; ssh_steps = ssh_steps[base_step], tn_buses = tn_buses[base_step], id_of = id_of)
  # pass 2, one direct import at a time; the template is not changed here,
  # so the steps are independent of each other
  for k in 1:n
    k == base_step && continue
    Rk = importCGMES(config; path = inputs[k].source, name = labels[k])
    _stampableBranches!(T, index, alias, Rk)
    deltas[k] = _snapshotDelta(T, index, parked, node_of_bus, part.node_sig, link_status[k], alias, Rk, k, labels[k], times[k]; ssh_steps = ssh_steps[k])
    sv[k] = _qstsSVReference(Rk, tn_buses[k], id_of, index)
  end
  messages = String[]
  isempty(base_residual) || push!(messages, string("base step: ", length(base_residual), " differences the delta cannot express, first: ", first(base_residual)))
  for k in 1:n
    isempty(deltas[k].residual) || push!(messages, string(labels[k], ": ", length(deltas[k].residual), " differences the delta cannot express"))
  end
  # the base delta against the template is empty by construction; anything
  # else means the template is not the base step
  (isempty(deltas[base_step].ops) && isempty(deltas[base_step].control_ops)) ||
    push!(messages, string("base step: the template differs from the base snapshot in ", length(deltas[base_step].ops), " op(s)"))
  # link ops per step by kind: a switch link is switching, a _TPL coupling
  # only follows the terminals an export puts on another TN; a link the
  # index does not address gets no op, so it is not looked up
  switch_ids = Set{Int}(index.id_by_mrid[m] for m in links.switches if haskey(index.id_by_mrid, m))
  tpl_ids = Set{Int}(index.id_by_mrid[m] for m in _tpl_mrids(part, links) if haskey(index.id_by_mrid, m))
  n_switch = [count(op -> op.target === :link && op.id in switch_ids, d.ops) for d in deltas]
  n_tpl = [count(op -> op.target === :link && op.id in tpl_ids, d.ops) for d in deltas]
  return CGMESSeriesUnion(T, index, base_step, labels, times, deltas, parked, part.node_of_term, part.node_mrid, part.node_sig, node_of_bus, link_status, alias, tn_buses,
                          links.switches, length(links.tpl_nodes), messages, sv, config.cgmes, n_switch, n_tpl)
end

# ISO 8601 UTC of a scenario time, "" without one
const _QSTS_TIME_FORMAT = Dates.DateFormat("yyyy-mm-ddTHH:MM:SS")
_qsts_timestamp(t::Union{Nothing,DateTime})::String = t === nothing ? "" : string(Dates.format(t, _QSTS_TIME_FORMAT), "Z")

# The input of a QSTS run on the union: step k applies `ops` and
# `control_ops` of its delta (the SSH positions of the snapshot), a snapshot
# run also its `start_ops` (each step starts from its own snapshot state; a
# sequential run carries the voltages from step to step and leaves them
# out), the delta carries the isolation, the step length is set when every
# step has a time and all gaps are equal.
function _qsts_input(u::CGMESSeriesUnion)::Sparlectra._QstsInput
  n = length(u.deltas)
  minutes = nothing
  if n >= 2 && all(t -> t !== nothing, u.scenario_times)
    gaps = [Dates.value(u.scenario_times[k + 1] - u.scenario_times[k]) / 60_000.0 for k in 1:(n - 1)]
    all(==(first(gaps)), gaps) && first(gaps) > 0.0 && (minutes = first(gaps))
  end
  return Sparlectra._QstsInput(string("union of ", u.labels[u.base_step]), :cgmes_series, u.template, u.index, [vcat(d.ops, d.control_ops) for d in u.deltas],
                               Vector{Sparlectra.PatchOp}[d.start_ops for d in u.deltas], copy(u.labels),
                               String[_qsts_timestamp(t) for t in u.scenario_times], u.base_step, false, Int[length(d.residual) for d in u.deltas],
                               Union{Nothing,Sparlectra.QstsSVReference}[r for r in u.sv], minutes, Dict{String,Any}())
end

"""
    compareSnapshotReproduction(union, k, Rk) -> Vector{String}

Apply step `k` of `union` (`ops` and `control_ops`) to a copy of the
template and compare the result with `Rk`, the direct import of snapshot
`k` (same configuration as the union build), by mRID: every field the delta
vocabulary covers, the parking state, the bus isolation, the link states and
the placement of every element on the buses of its TN. Returns the
differences as text; empty means the step is reproduced.
"""
function compareSnapshotReproduction(u::CGMESSeriesUnion, k::Int, Rk::CGMESImportResult)::Vector{String}
  d = u.deltas[k]
  net = deepcopy(u.template)
  Sparlectra.apply!(net, vcat(d.ops, d.control_ops), u.index)
  parked = copy(u.parked)
  for op in d.ops
    (op.op === :set && op.field === :in_service) || continue
    op.value == 0.0 ? push!(parked, op.id) : delete!(parked, op.id)
  end
  again = _snapshotDelta(net, u.index, parked, u.node_of_bus, u.node_sig, u.link_status[k], u.alias, Rk, k, d.label, d.scenario_time; ssh_steps = d.ssh_steps)
  out = String[]
  for op in again.ops
    push!(out, string("left after apply: ", op.op, " ", op.target, " ", op.id, " ", something(op.field, ""), " ", something(op.value, "")))
  end
  append!(out, again.residual)
  return out
end

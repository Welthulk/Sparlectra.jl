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

# file: src/adapters/powsybl/iidm_reader.jl
# purpose: read an IIDM XML file (.xiidm, .xml) into the same PowsyblTables
#          that pypowsybl's getters deliver, in plain Julia on EzXML: no
#          Python, no bundle. The reader resolves what pypowsybl resolves
#          before it hands out tables: the bus-breaker and bus views of a
#          node-breaker substation (union of nodes through internal
#          connections and closed switches; retained switches stay in the
#          bus-breaker view), the connected and synchronous components, the
#          current tap step of every transformer (rho, alpha and the step
#          corrections of r, x, g, b), the reactive limits at the active
#          setpoint from a capability curve, and the operational limits per
#          side and group. Column names, types and defaults follow
#          POWSYBL_SCHEMA, so build_net_from_powsybl and every bundle test
#          apply unchanged. The state columns (p, q, i, v_mag, v_angle) carry
#          the state written in the file; a bundle written by
#          tools/powsybl_dump.py carries OpenLoadFlow's solution instead.

"""
    read_iidm_tables(path; case = "") -> PowsyblTables

Read an IIDM XML file into PowSyBl tables without Python. `case` names the
network in the manifest (default: the file name without extension). Fails
with an `ArgumentError` naming the construct on a compressed file
(`.bz2`, `.gz`; unpack it first), on an XML file that is not an IIDM
network, and on the tie-line form of IIDM versions before 1.10 (the two
half lines inline); every other IIDM 1.x construct the builder reads is
mapped, unknown elements and extensions are ignored.
"""
function read_iidm_tables(path::AbstractString; case::AbstractString = "")::PowsyblTables
  lowered = lowercase(String(path))
  (endswith(lowered, ".bz2") || endswith(lowered, ".gz")) && throw(ArgumentError("IIDM file $(basename(String(path))) is compressed: unpack it (bzip2 -d / gzip -d) before the import"))
  isfile(path) || throw(ArgumentError("IIDM file $(String(path)) does not exist"))
  doc = EzXML.readxml(String(path))
  root = EzXML.root(doc)
  _iidm_name(root) == "network" || throw(ArgumentError("IIDM file $(basename(String(path))): root element is <$(_iidm_name(root))>, not <network>"))
  name = isempty(case) ? _iidm_case_name(String(path)) : String(case)
  net = _IidmNetwork(String(path), name, _iidm_version(root), _iidm_attrs(root))
  _iidm_collect!(net, root)
  _iidm_resolve_topology!(net)
  _iidm_resolve_components!(net)
  return _iidm_build_tables(net)
end

_iidm_case_name(path::String) = (b = basename(path); replace(replace(b, r"(?i)\.xiidm$" => ""), r"(?i)\.xml$" => ""))

# The local element name: EzXML returns the name without the namespace
# prefix, an explicit prefix is stripped for safety.
_iidm_name(el)::String = (n = EzXML.nodename(el); i = findlast(':', n); i === nothing ? n : n[i+1:end])

function _iidm_attrs(el)::Dict{String,String}
  d = Dict{String,String}()
  for a in EzXML.eachattribute(el)
    d[_iidm_name(a)] = EzXML.nodecontent(a)
  end
  return d
end

function _iidm_version(root)::String
  for (_, uri) in EzXML.namespaces(root)
    m = match(r"schema/iidm/([0-9]+_[0-9]+)$", uri)
    m === nothing || return String(m.captures[1])
  end
  return ""
end

_iidm_str(a::Dict{String,String}, k::String, d::String = "") = get(a, k, d)
function _iidm_float(a::Dict{String,String}, k::String, d::Float64 = NaN)::Float64
  s = get(a, k, nothing)
  s === nothing && return d
  v = tryparse(Float64, strip(s))
  return v === nothing ? d : v
end
function _iidm_int(a::Dict{String,String}, k::String, d::Int = -1)::Int
  s = get(a, k, nothing)
  s === nothing && return d
  v = tryparse(Int, strip(s))
  return v === nothing ? d : v
end
_iidm_bool(a::Dict{String,String}, k::String, d::Bool = false)::Bool = (s = get(a, k, nothing); s === nothing ? d : lowercase(strip(s)) == "true")
# powsybl stores a few factors as Java float; the table carries the float
# widened to double, as pypowsybl does (1.1 arrives as 1.100000023841858)
_iidm_float32(a::Dict{String,String}, k::String, d::Float64 = NaN)::Float64 = (v = _iidm_float(a, k, d); isfinite(v) ? Float64(Float32(v)) : v)

# --- the collected network ------------------------------------------------------

# One equipment terminal: where it sits (voltage level, node or configured
# bus) and what it is for the bus validity rule (:busbar, :branch, :feeder).
mutable struct _IidmTerminal
  equipment::String
  side::String                 # "" for an injection, ONE/TWO/THREE for a branch
  kind::Symbol
  vl::String
  node::Int                    # node-breaker node, -1 in a bus-breaker level
  bus::String                  # configured bus when connected, "" otherwise
  connectable_bus::String
  bb_bus::String               # resolved: bus-breaker view bus ("" when none)
  bus_id::String               # resolved: bus view bus ("" when none)
  connected::Bool              # resolved
end

mutable struct _IidmVoltageLevel
  id::String
  attrs::Dict{String,String}
  substation::String
  node_breaker::Bool
  config_buses::Vector{Dict{String,String}}        # bus-breaker topology, declaration order
  switches::Vector{Dict{String,String}}
  busbars::Vector{Dict{String,String}}
  internal::Vector{Tuple{Int,Int}}
  states::Vector{Dict{String,String}}              # <bus v angle nodes> of a node-breaker level
  terminals::Vector{_IidmTerminal}
  # resolved views
  bb_ids::Vector{String}                           # bus-breaker view buses, in order
  bb_members::Dict{String,Vector{Int}}             # node-breaker: nodes per bus-breaker bus
  bb_of_node::Dict{Int,String}
  bb_of_config::Dict{String,String}                # bus-breaker level: configured bus to itself when it exists
  bb_state::Dict{String,Tuple{Float64,Float64}}
  bb_name::Dict{String,String}
  bus_ids::Vector{String}                          # bus view buses, in order
  bus_of_bb::Dict{String,String}
  bus_state::Dict{String,Tuple{Float64,Float64}}
  bus_name::Dict{String,String}
end

_IidmVoltageLevel(id, attrs, substation, node_breaker) = _IidmVoltageLevel(id, attrs, substation, node_breaker, Dict{String,String}[], Dict{String,String}[], Dict{String,String}[], Tuple{Int,Int}[], Dict{String,String}[], _IidmTerminal[], String[], Dict{String,Vector{Int}}(), Dict{Int,String}(), Dict{String,String}(), Dict{String,Tuple{Float64,Float64}}(), Dict{String,String}(), String[], Dict{String,String}(), Dict{String,Tuple{Float64,Float64}}(), Dict{String,String}())

# An equipment element kept for the table pass: its attributes, its child
# elements of interest (tap changers, limits, curves, steps) and the
# terminals it registered.
struct _IidmEquipment
  kind::Symbol
  attrs::Dict{String,String}
  children::Vector{Tuple{String,Dict{String,String},Vector{Tuple{String,Dict{String,String},Vector{Tuple{String,Dict{String,String}}}}}}}
  terminals::Vector{_IidmTerminal}
  vl::String                   # the voltage level of an injection ("" for a branch)
end

mutable struct _IidmNetwork
  path::String
  case::String
  version::String
  attrs::Dict{String,String}
  substations::Vector{Dict{String,String}}
  vls::Vector{_IidmVoltageLevel}
  vl_index::Dict{String,Int}
  equipment::Vector{_IidmEquipment}
  terminal_index::Dict{Tuple{String,String},_IidmTerminal}
  # resolved components, keyed by bus view bus id
  connected_component::Dict{String,Int}
  synchronous_component::Dict{String,Int}
end

_IidmNetwork(path, case, version, attrs) = _IidmNetwork(path, case, version, attrs, Dict{String,String}[], _IidmVoltageLevel[], Dict{String,Int}(), _IidmEquipment[], Dict{Tuple{String,String},_IidmTerminal}(), Dict{String,Int}(), Dict{String,Int}())

# --- pass 1: collect ------------------------------------------------------------

# Child elements two levels deep, as (name, attrs, grandchildren): enough
# for tap changers with steps and a terminalRef, limit groups with limits
# and temporary limits, curves with points, non-linear shunt sections.
function _iidm_children(el)
  out = Tuple{String,Dict{String,String},Vector{Tuple{String,Dict{String,String},Vector{Tuple{String,Dict{String,String}}}}}}[]
  for c in EzXML.eachelement(el)
    grand = Tuple{String,Dict{String,String},Vector{Tuple{String,Dict{String,String}}}}[]
    for g in EzXML.eachelement(c)
      great = Tuple{String,Dict{String,String}}[]
      for h in EzXML.eachelement(g)
        push!(great, (_iidm_name(h), _iidm_attrs(h)))
      end
      push!(grand, (_iidm_name(g), _iidm_attrs(g), great))
    end
    push!(out, (_iidm_name(c), _iidm_attrs(c), grand))
  end
  return out
end

const _IIDM_INJECTIONS = Dict{String,Tuple{Symbol,Symbol}}(
  "generator" => (:generator, :feeder),
  "battery" => (:battery, :feeder),
  "load" => (:load, :feeder),
  "shuntCompensator" => (:shunt, :feeder),
  "staticVarCompensator" => (:svc, :feeder),
  "vscConverterStation" => (:vsc, :branch),
  "lccConverterStation" => (:lcc, :branch),
  "danglingLine" => (:dangling_line, :branch),
  "boundaryLine" => (:dangling_line, :branch),
  "ground" => (:ground, :feeder),
)

function _iidm_register!(net::_IidmNetwork, t::_IidmTerminal)
  net.terminal_index[(t.equipment, t.side)] = t
  vl = net.vls[net.vl_index[t.vl]]
  push!(vl.terminals, t)
  return t
end

# The terminal of an injection from its attributes (bus/connectableBus or node).
function _iidm_injection_terminal(net::_IidmNetwork, id::String, kind::Symbol, vlid::String, a::Dict{String,String})
  haskey(net.vl_index, vlid) || throw(ArgumentError("IIDM: $(id) sits in voltage level $(vlid), which the file does not declare"))
  t = _IidmTerminal(id, "", kind, vlid, _iidm_int(a, "node", -1), _iidm_str(a, "bus"), _iidm_str(a, "connectableBus"), "", "", false)
  return _iidm_register!(net, t)
end

function _iidm_branch_terminal(net::_IidmNetwork, id::String, side::String, k::Int, a::Dict{String,String})
  vlid = _iidm_str(a, "voltageLevelId$(k)")
  haskey(net.vl_index, vlid) || throw(ArgumentError("IIDM: $(id) side $(k) names voltage level $(repr(vlid)), which the file does not declare"))
  t = _IidmTerminal(id, side, :branch, vlid, _iidm_int(a, "node$(k)", -1), _iidm_str(a, "bus$(k)"), _iidm_str(a, "connectableBus$(k)"), "", "", false)
  return _iidm_register!(net, t)
end

const _IIDM_SIDES = ("ONE", "TWO", "THREE")

function _iidm_collect_voltage_level!(net::_IidmNetwork, el, substation::String)
  a = _iidm_attrs(el)
  id = _iidm_str(a, "id")
  vl = _IidmVoltageLevel(id, a, substation, _iidm_str(a, "topologyKind") == "NODE_BREAKER")
  push!(net.vls, vl)
  net.vl_index[id] = length(net.vls)
  for c in EzXML.eachelement(el)
    n = _iidm_name(c)
    if n == "busBreakerTopology" || n == "nodeBreakerTopology"
      for t in EzXML.eachelement(c)
        tn = _iidm_name(t)
        ta = _iidm_attrs(t)
        if tn == "bus"
          haskey(ta, "nodes") ? push!(vl.states, ta) : push!(vl.config_buses, ta)
        elseif tn == "switch"
          push!(vl.switches, ta)
        elseif tn == "busbarSection"
          push!(vl.busbars, ta)
          _iidm_register!(net, _IidmTerminal(_iidm_str(ta, "id"), "", :busbar, id, _iidm_int(ta, "node", -1), "", "", "", "", false))
        elseif tn == "internalConnection"
          push!(vl.internal, (_iidm_int(ta, "node1"), _iidm_int(ta, "node2")))
        end
      end
    elseif haskey(_IIDM_INJECTIONS, n)
      kind, tkind = _IIDM_INJECTIONS[n]
      ca = _iidm_attrs(c)
      eid = _iidm_str(ca, "id")
      t = _iidm_injection_terminal(net, eid, tkind, id, ca)
      kind in (:battery, :ground) && continue   # counted for the bus validity, no table
      push!(net.equipment, _IidmEquipment(kind, ca, _iidm_children(c), [t], id))
    end
  end
  return nothing
end

function _iidm_collect_branch!(net::_IidmNetwork, el, kind::Symbol, nsides::Int)
  a = _iidm_attrs(el)
  id = _iidm_str(a, "id")
  ts = [_iidm_branch_terminal(net, id, _IIDM_SIDES[k], k, a) for k in 1:nsides]
  push!(net.equipment, _IidmEquipment(kind, a, _iidm_children(el), ts, ""))
  return nothing
end

function _iidm_collect!(net::_IidmNetwork, root)
  # voltage levels first (they may sit under a substation or directly under
  # the network since IIDM 1.6), branches after all levels are known
  branches = Any[]
  for el in EzXML.eachelement(root)
    n = _iidm_name(el)
    if n == "substation"
      sa = _iidm_attrs(el)
      push!(net.substations, sa)
      for c in EzXML.eachelement(el)
        cn = _iidm_name(c)
        if cn == "voltageLevel"
          _iidm_collect_voltage_level!(net, c, _iidm_str(sa, "id"))
        elseif cn in ("twoWindingsTransformer", "threeWindingsTransformer")
          push!(branches, c)
        end
      end
    elseif n == "voltageLevel"
      _iidm_collect_voltage_level!(net, el, "")
    elseif n in ("line", "tieLine", "hvdcLine")
      push!(branches, el)
    end
  end
  for el in branches
    n = _iidm_name(el)
    if n == "line"
      _iidm_collect_branch!(net, el, :line, 2)
    elseif n == "twoWindingsTransformer"
      _iidm_collect_branch!(net, el, :two_wt, 2)
    elseif n == "threeWindingsTransformer"
      _iidm_collect_branch!(net, el, :three_wt, 3)
    elseif n == "tieLine"
      a = _iidm_attrs(el)
      haskey(a, "id_1") && throw(ArgumentError("IIDM file $(basename(net.path)): tie line $(_iidm_str(a, "id")) uses the inline half-line form of IIDM versions before 1.10, which this reader does not map; export the network with a current PowSyBl version"))
      push!(net.equipment, _IidmEquipment(:tie_line, a, _iidm_children(el), _IidmTerminal[], ""))
    elseif n == "hvdcLine"
      push!(net.equipment, _IidmEquipment(:hvdc, _iidm_attrs(el), _iidm_children(el), _IidmTerminal[], ""))
    end
  end
  return nothing
end

# --- pass 2: topology ------------------------------------------------------------

# A small union-find over integer keys.
function _iidm_uf_find(parent::Dict{Int,Int}, x::Int)::Int
  r = x
  while parent[r] != r
    r = parent[r]
  end
  while parent[x] != r
    nx = parent[x]
    parent[x] = r
    x = nx
  end
  return r
end
function _iidm_uf_union!(parent::Dict{Int,Int}, x::Int, y::Int)
  rx = _iidm_uf_find(parent, x)
  ry = _iidm_uf_find(parent, y)
  rx == ry && return nothing
  parent[max(rx, ry)] = min(rx, ry)
  return nothing
end

# Groups of nodes as sorted vectors, ordered by their smallest node.
function _iidm_groups(parent::Dict{Int,Int})::Vector{Vector{Int}}
  by_root = Dict{Int,Vector{Int}}()
  for k in keys(parent)
    push!(get!(by_root, _iidm_uf_find(parent, k), Int[]), k)
  end
  groups = [sort!(g) for g in values(by_root)]
  sort!(groups; by = first)
  return groups
end

# powsybl's bus validity: a bus of the bus view needs a busbar section with
# a feeder, or a branch with two feeders; a bus of the bus-breaker view of
# a node-breaker level needs one terminal; a bus of the bus view of a
# bus-breaker level needs a branch.
function _iidm_bus_view_valid(terminals)::Bool
  feeders = 0
  branches = 0
  busbars = 0
  for t in terminals
    if t.kind == :branch
      branches += 1
      feeders += 1
    elseif t.kind == :feeder
      feeders += 1
    else
      busbars += 1
    end
  end
  return (busbars >= 1 && feeders >= 1) || (branches >= 1 && feeders >= 2)
end

function _iidm_resolve_node_breaker!(vl::_IidmVoltageLevel)
  by_node = Dict{Int,Vector{_IidmTerminal}}()
  for t in vl.terminals
    t.node >= 0 || throw(ArgumentError("IIDM: $(t.equipment) in node-breaker level $(vl.id) has no node"))
    push!(get!(by_node, t.node, _IidmTerminal[]), t)
  end
  nodes = Set{Int}(keys(by_node))
  for s in vl.switches
    push!(nodes, _iidm_int(s, "node1"))
    push!(nodes, _iidm_int(s, "node2"))
  end
  for (n1, n2) in vl.internal
    push!(nodes, n1)
    push!(nodes, n2)
  end
  for st in vl.states
    for tok in split(_iidm_str(st, "nodes"), ',')
      v = tryparse(Int, strip(tok))
      v === nothing || push!(nodes, v)
    end
  end
  # bus-breaker view: closed non-retained switches and internal connections
  # merge nodes; bus view: every closed switch merges
  bb = Dict{Int,Int}(n => n for n in nodes)
  bv = Dict{Int,Int}(n => n for n in nodes)
  for (n1, n2) in vl.internal
    _iidm_uf_union!(bb, n1, n2)
    _iidm_uf_union!(bv, n1, n2)
  end
  for s in vl.switches
    _iidm_bool(s, "open", false) && continue
    n1 = _iidm_int(s, "node1")
    n2 = _iidm_int(s, "node2")
    _iidm_uf_union!(bv, n1, n2)
    _iidm_bool(s, "retained", false) || _iidm_uf_union!(bb, n1, n2)
  end
  for g in _iidm_groups(bb)
    any(n -> haskey(by_node, n), g) || continue
    id = vl.id * "_" * string(first(g))
    push!(vl.bb_ids, id)
    vl.bb_members[id] = g
    vl.bb_name[id] = ""
    for n in g
      vl.bb_of_node[n] = id
    end
  end
  for g in _iidm_groups(bv)
    ts = _IidmTerminal[]
    for n in g
      append!(ts, get(by_node, n, _IidmTerminal[]))
    end
    _iidm_bus_view_valid(ts) || continue
    id = vl.id * "_" * string(first(g))
    push!(vl.bus_ids, id)
    for n in g
      bbid = get(vl.bb_of_node, n, "")
      isempty(bbid) || (vl.bus_of_bb[bbid] = id)
    end
    vl.bus_state[id] = (NaN, NaN)
    for st in vl.states
      toks = split(_iidm_str(st, "nodes"), ',')
      first_node = isempty(toks) ? nothing : tryparse(Int, strip(first(toks)))
      (first_node !== nothing && first_node in g) || continue
      vl.bus_state[id] = (_iidm_float(st, "v"), _iidm_float(st, "angle"))
      break
    end
  end
  for id in vl.bb_ids
    vl.bb_state[id] = get(vl.bus_state, get(vl.bus_of_bb, id, ""), (NaN, NaN))
  end
  for t in vl.terminals
    t.bb_bus = get(vl.bb_of_node, t.node, "")
    t.bus_id = get(vl.bus_of_bb, t.bb_bus, "")
    t.connected = !isempty(t.bus_id)
  end
  return nothing
end

function _iidm_resolve_bus_breaker!(vl::_IidmVoltageLevel)
  index = Dict{String,Int}()
  for (k, b) in enumerate(vl.config_buses)
    id = _iidm_str(b, "id")
    index[id] = k
    push!(vl.bb_ids, id)
    vl.bb_of_config[id] = id
    vl.bb_name[id] = _iidm_str(b, "name")
    vl.bb_state[id] = (_iidm_float(b, "v"), _iidm_float(b, "angle"))
  end
  by_bus = Dict{String,Vector{_IidmTerminal}}()
  for t in vl.terminals
    isempty(t.bus) || push!(get!(by_bus, t.bus, _IidmTerminal[]), t)
  end
  parent = Dict{Int,Int}(k => k for k in 1:length(vl.config_buses))
  for s in vl.switches
    _iidm_bool(s, "open", false) && continue
    k1 = get(index, _iidm_str(s, "bus1"), 0)
    k2 = get(index, _iidm_str(s, "bus2"), 0)
    (k1 > 0 && k2 > 0) && _iidm_uf_union!(parent, k1, k2)
  end
  for g in _iidm_groups(parent)
    branches = 0
    for k in g
      branches += count(t -> t.kind == :branch, get(by_bus, vl.bb_ids[k], _IidmTerminal[]))
    end
    branches >= 1 || continue
    id = vl.id * "_" * string(first(g) - 1)
    push!(vl.bus_ids, id)
    # powsybl names the merged bus after the voltage level's name when it
    # has one ("110.0_0" in the CGMES-derived MicroGrid), else it stays empty
    vl.bus_name[id] = isempty(_iidm_str(vl.attrs, "name")) ? "" : _iidm_str(vl.attrs, "name") * "_" * string(first(g) - 1)
    state = (NaN, NaN)
    for k in g
      vl.bus_of_bb[vl.bb_ids[k]] = id
      st = vl.bb_state[vl.bb_ids[k]]
      (isnan(state[1]) && isfinite(st[1])) && (state = st)
    end
    vl.bus_state[id] = state
  end
  for t in vl.terminals
    t.connected = !isempty(t.bus)
    t.bb_bus = t.connected ? t.bus : t.connectable_bus
    t.bus_id = t.connected ? get(vl.bus_of_bb, t.bus, "") : ""
  end
  return nothing
end

function _iidm_resolve_topology!(net::_IidmNetwork)
  for vl in net.vls
    vl.node_breaker ? _iidm_resolve_node_breaker!(vl) : _iidm_resolve_bus_breaker!(vl)
  end
  return nothing
end

# --- pass 3: components -----------------------------------------------------------

# Connected components over the bus view: AC branches for the synchronous
# components, AC branches plus HVDC links for the connected components.
# Numbered by decreasing size, ties in order of first appearance, as
# powsybl's ComponentsManager numbers them.
function _iidm_resolve_components!(net::_IidmNetwork)
  order = String[]
  for vl in net.vls
    append!(order, vl.bus_ids)
  end
  index = Dict{String,Int}(id => k for (k, id) in enumerate(order))
  ac = Dict{Int,Int}(k => k for k in 1:length(order))
  all = Dict{Int,Int}(k => k for k in 1:length(order))
  join!(uf, a::String, b::String) = ((isempty(a) || isempty(b)) ? nothing : _iidm_uf_union!(uf, index[a], index[b]))
  dangling_bus = Dict{String,String}()
  stations = Dict{String,String}()
  for e in net.equipment
    if e.kind in (:line, :two_wt)
      t1, t2 = e.terminals
      join!(ac, t1.bus_id, t2.bus_id)
      join!(all, t1.bus_id, t2.bus_id)
    elseif e.kind == :three_wt
      buses = [t.bus_id for t in e.terminals if !isempty(t.bus_id)]
      for k in 2:length(buses)
        join!(ac, buses[1], buses[k])
        join!(all, buses[1], buses[k])
      end
    elseif e.kind == :dangling_line
      dangling_bus[_iidm_str(e.attrs, "id")] = e.terminals[1].bus_id
    elseif e.kind in (:vsc, :lcc)
      stations[_iidm_str(e.attrs, "id")] = e.terminals[1].bus_id
    end
  end
  for e in net.equipment
    if e.kind == :tie_line
      b1 = get(dangling_bus, _iidm_tie_half(e.attrs, 1), "")
      b2 = get(dangling_bus, _iidm_tie_half(e.attrs, 2), "")
      join!(ac, b1, b2)
      join!(all, b1, b2)
    elseif e.kind == :hvdc
      join!(all, get(stations, _iidm_str(e.attrs, "converterStation1"), ""), get(stations, _iidm_str(e.attrs, "converterStation2"), ""))
    end
  end
  net.synchronous_component = _iidm_number_components(ac, order)
  net.connected_component = _iidm_number_components(all, order)
  return nothing
end

function _iidm_number_components(uf::Dict{Int,Int}, order::Vector{String})::Dict{String,Int}
  groups = _iidm_groups(uf)                      # ordered by first appearance
  perm = sortperm(groups; by = g -> -length(g))  # stable: ties keep the appearance order
  numbers = Dict{String,Int}()
  for (c, gi) in enumerate(perm)
    for k in groups[gi]
      numbers[order[k]] = c - 1
    end
  end
  return numbers
end

_iidm_tie_half(a::Dict{String,String}, k::Int)::String = (v = _iidm_str(a, "boundaryLineId$(k)"); isempty(v) ? _iidm_str(a, "danglingLineId$(k)") : v)

# --- pass 4: tables ----------------------------------------------------------------

# A row is a Dict of column values; the table takes the schema columns in
# order, filling what a row does not carry with the type's empty value.
const _IidmRow = Dict{Symbol,Any}

function _iidm_table(name::String, rows::Vector{_IidmRow})::NamedTuple
  columns = POWSYBL_SCHEMA[name]
  vecs = Any[]
  for c in columns
    T = c.eltype
    empty = T === String ? "" : T === Float64 ? NaN : T === Int ? -1 : false
    v = _powsybl_vector_type(T)(undef, length(rows))
    for (i, r) in enumerate(rows)
      x = get(r, c.name, empty)
      v[i] = T === String ? String(x) : T === Float64 ? Float64(x) : T === Int ? Int(x) : Bool(x)
    end
    push!(vecs, v)
  end
  return NamedTuple{Tuple(c.name for c in columns)}(Tuple(vecs))
end

# The bus attachment columns of an injection.
function _iidm_attachment!(r::_IidmRow, t::_IidmTerminal)
  r[:voltage_level_id] = t.vl
  r[:bus_id] = t.bus_id
  r[:bus_breaker_bus_id] = t.bb_bus
  r[:node] = t.node
  r[:connected] = t.connected
  return r
end

function _iidm_side!(r::_IidmRow, t::_IidmTerminal, k::Int)
  r[Symbol("voltage_level$(k)_id")] = t.vl
  r[Symbol("bus$(k)_id")] = t.bus_id
  r[Symbol("bus_breaker_bus$(k)_id")] = t.bb_bus
  r[Symbol("node$(k)")] = t.node
  r[Symbol("connected$(k)")] = t.connected
  return r
end

# The terminal a regulatingTerminal / terminalRef child names, or the
# element's own terminal.
function _iidm_ref_terminal(net::_IidmNetwork, children, own::_IidmTerminal, names)
  for (n, a, _) in children
    n in names || continue
    id = _iidm_str(a, "id")
    side = _iidm_str(a, "side")
    t = get(net.terminal_index, (id, side), nothing)
    t === nothing && (t = get(net.terminal_index, (id, ""), nothing))
    t === nothing && (t = get(net.terminal_index, (id, "ONE"), nothing))
    return (t === nothing ? own : t, side)
  end
  return (own, "")
end

# Reactive limits of a generator or VSC station: kind, plain bounds (NaN for
# a curve), the bounds at a given active power from the curve, the points.
function _iidm_reactive_limits(children)
  for (n, a, grand) in children
    if n == "minMaxReactiveLimits"
      return (kind = "MIN_MAX", min_q = _iidm_float(a, "minQ", -floatmax(Float64)), max_q = _iidm_float(a, "maxQ", floatmax(Float64)), points = NTuple{3,Float64}[])
    elseif n == "reactiveCapabilityCurve"
      pts = NTuple{3,Float64}[(_iidm_float(pa, "p"), _iidm_float(pa, "minQ"), _iidm_float(pa, "maxQ")) for (pn, pa, _) in grand if pn == "point"]
      sort!(pts; by = first)
      return (kind = "CURVE", min_q = NaN, max_q = NaN, points = pts)
    end
  end
  return (kind = "MIN_MAX", min_q = -floatmax(Float64), max_q = floatmax(Float64), points = NTuple{3,Float64}[])
end

# powsybl's curve lookup: exact point, clamped outside the range, linear
# between the neighbouring points inside.
function _iidm_curve_at(pts::Vector{NTuple{3,Float64}}, p::Float64, which::Int)::Float64
  isempty(pts) && return NaN
  isfinite(p) || return NaN
  p <= pts[1][1] && return pts[1][which]
  p >= pts[end][1] && return pts[end][which]
  for k in 1:length(pts)-1
    a, b = pts[k], pts[k+1]
    a[1] == p && return a[which]
    if a[1] < p < b[1]
      return a[which] + (b[which] - a[which]) * (p - a[1]) / (b[1] - a[1])
    end
  end
  return pts[end][which]
end

function _iidm_limits_at(lim, p::Float64)
  lim.kind == "CURVE" || return (lim.min_q, lim.max_q)
  return (_iidm_curve_at(lim.points, p, 2), _iidm_curve_at(lim.points, p, 3))
end

# The current step of a tap changer child (ratioTapChanger, phaseTapChanger,
# with a leg suffix on a three-winding transformer): the step attributes at
# the tap position, the position and the step list.
function _iidm_tap_changer(children, name::String)
  for (n, a, grand) in children
    n == name || continue
    low = _iidm_int(a, "lowTapPosition", 0)
    tap = _iidm_int(a, "tapPosition", low)
    steps = Dict{String,String}[sa for (sn, sa, _) in grand if sn == "step"]
    k = tap - low + 1
    (1 <= k <= length(steps)) || throw(ArgumentError("IIDM: tap changer $(name) tap position $(tap) is outside its steps $(low)..$(low + length(steps) - 1)"))
    return (attrs = a, grand = grand, low = low, tap = tap, steps = steps, step = steps[k])
  end
  return nothing
end

_iidm_step_factor(tc, key::String)::Float64 = tc === nothing ? 1.0 : 1.0 + _iidm_float(tc.step, key, 0.0) / 100.0
_iidm_step_rho(tc)::Float64 = tc === nothing ? 1.0 : _iidm_float(tc.step, "rho", 1.0)

function _iidm_tap_changer_rows!(net::_IidmNetwork, ratio_rows, ratio_steps, phase_rows, phase_steps, id::String, side::String, own::_IidmTerminal, rtc, ptc)
  if rtc !== nothing
    a = rtc.attrs
    reg, reg_side = _iidm_ref_terminal(net, rtc.grand, own, ("terminalRef",))
    mode = _iidm_str(a, "regulationMode", "VOLTAGE")
    target_v = haskey(a, "targetV") ? _iidm_float(a, "targetV") : (mode == "VOLTAGE" ? _iidm_float(a, "regulationValue") : NaN)
    push!(ratio_rows, _IidmRow(:id => id, :side => side, :tap => rtc.tap, :solved_tap_position => Float64(_iidm_int(a, "solvedTapPosition", rtc.tap)), :low_tap => rtc.low, :high_tap => rtc.low + length(rtc.steps) - 1, :step_count => length(rtc.steps), :oltc => _iidm_bool(a, "loadTapChangingCapabilities", false), :regulating => _iidm_bool(a, "regulating", false), :target_v => target_v, :target_deadband => _iidm_float(a, "targetDeadband"), :regulating_bus_id => reg.bus_id, :regulated_side => reg_side))
    for (k, s) in enumerate(rtc.steps)
      push!(ratio_steps, _IidmRow(:id => id, :position => rtc.low + k - 1, :side => side, :rho => _iidm_float(s, "rho", 1.0), :r => _iidm_float(s, "r", 0.0), :x => _iidm_float(s, "x", 0.0), :g => _iidm_float(s, "g", 0.0), :b => _iidm_float(s, "b", 0.0)))
    end
  end
  if ptc !== nothing
    a = ptc.attrs
    reg, reg_side = _iidm_ref_terminal(net, ptc.grand, own, ("terminalRef",))
    push!(phase_rows, _IidmRow(:id => id, :side => side, :tap => ptc.tap, :solved_tap_position => Float64(_iidm_int(a, "solvedTapPosition", ptc.tap)), :low_tap => ptc.low, :high_tap => ptc.low + length(ptc.steps) - 1, :step_count => length(ptc.steps), :oltc => _iidm_bool(a, "loadTapChangingCapabilities", false), :regulating => _iidm_bool(a, "regulating", false), :regulation_mode => _iidm_str(a, "regulationMode", "FIXED_TAP"), :regulation_value => _iidm_float(a, "regulationValue"), :target_deadband => _iidm_float(a, "targetDeadband"), :regulating_bus_id => reg.bus_id, :regulated_side => reg_side))
    for (k, s) in enumerate(ptc.steps)
      push!(phase_steps, _IidmRow(:id => id, :position => ptc.low + k - 1, :side => side, :rho => _iidm_float(s, "rho", 1.0), :alpha => _iidm_float(s, "alpha", 0.0), :r => _iidm_float(s, "r", 0.0), :x => _iidm_float(s, "x", 0.0), :g => _iidm_float(s, "g", 0.0), :b => _iidm_float(s, "b", 0.0)))
    end
  end
  return nothing
end

# Operational limits of one element side: the limit groups (IIDM 1.12 and
# later) or the bare currentLimits / apparentPowerLimits / activePowerLimits
# children of older versions (group "DEFAULT").
const _IIDM_LIMIT_TYPES = Dict{String,String}("currentLimits" => "CURRENT", "apparentPowerLimits" => "APPARENT_POWER", "activePowerLimits" => "ACTIVE_POWER")

function _iidm_limit_rows!(rows, children, id::String, side::String, suffix::String, element_type::String, selected::String)
  groups = Tuple{String,Vector{Tuple{String,Dict{String,String},Vector{Tuple{String,Dict{String,String}}}}}}[]
  for (n, a, grand) in children
    if n == "operationalLimitsGroup" * suffix
      push!(groups, (_iidm_str(a, "id"), grand))
    elseif haskey(_IIDM_LIMIT_TYPES, n) && suffix == "" || (endswith(n, suffix) && !isempty(suffix) && haskey(_IIDM_LIMIT_TYPES, n[1:end-length(suffix)]))
      base = isempty(suffix) ? n : n[1:end-length(suffix)]
      temps = Tuple{String,Dict{String,String}}[(gn, ga) for (gn, ga, _) in grand]
      push!(groups, ("DEFAULT", [(base, a, temps)]))
    end
  end
  for (gname, limits) in groups
    for (ln, la, temps) in limits
      type = get(_IIDM_LIMIT_TYPES, ln, nothing)
      type === nothing && continue
      is_selected = isempty(selected) ? gname == "DEFAULT" : gname == selected
      push!(rows, _IidmRow(:element_id => id, :side => side, :type => type, :acceptable_duration => -1, :group_name => gname, :element_type => element_type, :name => "permanent_limit", :value => _iidm_float(la, "permanentLimit"), :fictitious => false, :selected => is_selected))
      for (tn, ta) in temps
        tn == "temporaryLimit" || continue
        push!(rows, _IidmRow(:element_id => id, :side => side, :type => type, :acceptable_duration => _iidm_int(ta, "acceptableDuration", Int(typemax(Int32))), :group_name => gname, :element_type => element_type, :name => _iidm_str(ta, "name"), :value => _iidm_float(ta, "value", floatmax(Float64)), :fictitious => _iidm_bool(ta, "fictitious", false), :selected => is_selected))
      end
    end
  end
  return nothing
end

_iidm_selected_group(a::Dict{String,String}, k::String)::String = (v = _iidm_str(a, "selectedOperationalLimitsGroupIds" * k); isempty(v) ? "DEFAULT" : v)

function _iidm_build_tables(net::_IidmNetwork)::PowsyblTables
  T = Dict{String,Vector{_IidmRow}}(name => _IidmRow[] for name in POWSYBL_TABLE_NAMES)
  for s in net.substations
    push!(T["substations"], _IidmRow(:id => _iidm_str(s, "id"), :name => _iidm_str(s, "name"), :TSO => _iidm_str(s, "tso"), :geo_tags => _iidm_str(s, "geographicalTags"), :country => _iidm_str(s, "country"), :fictitious => _iidm_bool(s, "fictitious")))
  end
  for vl in net.vls
    a = vl.attrs
    push!(T["voltage_levels"], _IidmRow(:id => vl.id, :name => _iidm_str(a, "name"), :substation_id => vl.substation, :nominal_v => _iidm_float(a, "nominalV"), :high_voltage_limit => _iidm_float(a, "highVoltageLimit"), :low_voltage_limit => _iidm_float(a, "lowVoltageLimit"), :fictitious => _iidm_bool(a, "fictitious"), :topology_kind => _iidm_str(a, "topologyKind", "BUS_BREAKER")))
    for id in vl.bus_ids
      v, ang = vl.bus_state[id]
      push!(T["buses"], _IidmRow(:id => id, :name => get(vl.bus_name, id, ""), :v_mag => v, :v_angle => ang, :connected_component => get(net.connected_component, id, -1), :synchronous_component => get(net.synchronous_component, id, -1), :voltage_level_id => vl.id, :fictitious => false, :fictitious_p0 => 0.0, :fictitious_q0 => 0.0))
    end
    for id in vl.bb_ids
      v, ang = vl.bb_state[id]
      bus = get(vl.bus_of_bb, id, "")
      push!(T["bus_breaker_view_buses"], _IidmRow(:id => id, :name => vl.bb_name[id], :v_mag => v, :v_angle => ang, :connected_component => get(net.connected_component, bus, -1), :synchronous_component => get(net.synchronous_component, bus, -1), :voltage_level_id => vl.id, :bus_id => bus, :fictitious => false, :fictitious_p0 => 0.0, :fictitious_q0 => 0.0))
    end
    for s in vl.switches
      if vl.node_breaker
        retained = _iidm_bool(s, "retained", false)
        n1 = _iidm_int(s, "node1")
        n2 = _iidm_int(s, "node2")
        push!(T["switches"], _IidmRow(:id => _iidm_str(s, "id"), :name => _iidm_str(s, "name"), :kind => _iidm_str(s, "kind"), :open => _iidm_bool(s, "open"), :retained => retained, :voltage_level_id => vl.id, :bus_breaker_bus1_id => retained ? get(vl.bb_of_node, n1, "") : "", :bus_breaker_bus2_id => retained ? get(vl.bb_of_node, n2, "") : "", :node1 => n1, :node2 => n2, :fictitious => _iidm_bool(s, "fictitious")))
      else
        push!(T["switches"], _IidmRow(:id => _iidm_str(s, "id"), :name => _iidm_str(s, "name"), :kind => _iidm_str(s, "kind"), :open => _iidm_bool(s, "open"), :retained => true, :voltage_level_id => vl.id, :bus_breaker_bus1_id => _iidm_str(s, "bus1"), :bus_breaker_bus2_id => _iidm_str(s, "bus2"), :node1 => -1, :node2 => -1, :fictitious => _iidm_bool(s, "fictitious")))
      end
    end
    for b in vl.busbars
      id = _iidm_str(b, "id")
      t = net.terminal_index[(id, "")]
      v, ang = get(vl.bb_state, t.bb_bus, (NaN, NaN))
      push!(T["busbar_sections"], _IidmRow(:id => id, :name => _iidm_str(b, "name"), :v => v, :angle => ang, :voltage_level_id => vl.id, :bus_id => t.bus_id, :bus_breaker_bus_id => t.bb_bus, :node => t.node, :connected => t.connected, :fictitious => _iidm_bool(b, "fictitious")))
    end
  end
  dangling = Dict{String,_IidmEquipment}()
  station_line = Dict{String,String}()
  for e in net.equipment
    e.kind == :dangling_line && (dangling[_iidm_str(e.attrs, "id")] = e)
    if e.kind == :hvdc
      station_line[_iidm_str(e.attrs, "converterStation1")] = _iidm_str(e.attrs, "id")
      station_line[_iidm_str(e.attrs, "converterStation2")] = _iidm_str(e.attrs, "id")
    end
  end
  tie_of = Dict{String,String}()
  for e in net.equipment
    e.kind == :tie_line || continue
    tie_of[_iidm_tie_half(e.attrs, 1)] = _iidm_str(e.attrs, "id")
    tie_of[_iidm_tie_half(e.attrs, 2)] = _iidm_str(e.attrs, "id")
  end
  for e in net.equipment
    a = e.attrs
    id = _iidm_str(a, "id")
    name = _iidm_str(a, "name")
    fict = _iidm_bool(a, "fictitious")
    if e.kind == :line
      r = _IidmRow(:id => id, :name => name, :r => _iidm_float(a, "r"), :x => _iidm_float(a, "x"), :g1 => _iidm_float(a, "g1", 0.0), :b1 => _iidm_float(a, "b1", 0.0), :g2 => _iidm_float(a, "g2", 0.0), :b2 => _iidm_float(a, "b2", 0.0), :p1 => _iidm_float(a, "p1"), :q1 => _iidm_float(a, "q1"), :p2 => _iidm_float(a, "p2"), :q2 => _iidm_float(a, "q2"), :fictitious => fict, :selected_limits_group_1 => _iidm_selected_group(a, "1"), :selected_limits_group_2 => _iidm_selected_group(a, "2"))
      _iidm_side!(r, e.terminals[1], 1)
      _iidm_side!(r, e.terminals[2], 2)
      push!(T["lines"], r)
      for k in 1:2
        _iidm_limit_rows!(T["operational_limits"], e.children, id, _IIDM_SIDES[k], string(k), "LINE", _iidm_str(a, "selectedOperationalLimitsGroupIds$(k)"))
      end
    elseif e.kind == :two_wt
      rtc = _iidm_tap_changer(e.children, "ratioTapChanger")
      ptc = _iidm_tap_changer(e.children, "phaseTapChanger")
      u1 = _iidm_float(a, "ratedU1")
      u2 = _iidm_float(a, "ratedU2")
      r = _IidmRow(:id => id, :name => name, :r => _iidm_float(a, "r"), :x => _iidm_float(a, "x"), :g => _iidm_float(a, "g", 0.0), :b => _iidm_float(a, "b", 0.0), :rated_u1 => u1, :rated_u2 => u2, :rated_s => _iidm_float(a, "ratedS"), :p1 => _iidm_float(a, "p1"), :q1 => _iidm_float(a, "q1"), :p2 => _iidm_float(a, "p2"), :q2 => _iidm_float(a, "q2"), :fictitious => fict, :selected_limits_group_1 => _iidm_selected_group(a, "1"), :selected_limits_group_2 => _iidm_selected_group(a, "2"))
      r[:rho] = u2 / u1 * _iidm_step_rho(rtc) * _iidm_step_rho(ptc)
      r[:alpha] = ptc === nothing ? 0.0 : _iidm_float(ptc.step, "alpha", 0.0)
      for key in ("r", "x", "g", "b")
        r[Symbol(key * "_at_current_tap")] = r[Symbol(key)] * _iidm_step_factor(rtc, key) * _iidm_step_factor(ptc, key)
      end
      _iidm_side!(r, e.terminals[1], 1)
      _iidm_side!(r, e.terminals[2], 2)
      push!(T["2_windings_transformers"], r)
      _iidm_tap_changer_rows!(net, T["ratio_tap_changers"], T["ratio_tap_changer_steps"], T["phase_tap_changers"], T["phase_tap_changer_steps"], id, "", e.terminals[1], rtc, ptc)
      for k in 1:2
        _iidm_limit_rows!(T["operational_limits"], e.children, id, _IIDM_SIDES[k], string(k), "TWO_WINDINGS_TRANSFORMER", _iidm_str(a, "selectedOperationalLimitsGroupIds$(k)"))
      end
    elseif e.kind == :three_wt
      u0 = _iidm_float(a, "ratedU0", _iidm_float(a, "ratedU1"))
      r = _IidmRow(:id => id, :name => name, :rated_u0 => u0, :fictitious => fict)
      for k in 1:3
        rtc = _iidm_tap_changer(e.children, "ratioTapChanger$(k)")
        ptc = _iidm_tap_changer(e.children, "phaseTapChanger$(k)")
        uk = _iidm_float(a, "ratedU$(k)")
        for key in ("r", "x", "g", "b")
          v = _iidm_float(a, key * string(k), key in ("g", "b") ? 0.0 : NaN)
          r[Symbol(key * string(k))] = v
          r[Symbol("$(key)$(k)_at_current_tap")] = v * _iidm_step_factor(rtc, key) * _iidm_step_factor(ptc, key)
        end
        r[Symbol("rated_u$(k)")] = uk
        r[Symbol("rated_s$(k)")] = _iidm_float(a, "ratedS$(k)")
        r[Symbol("ratio_tap_position$(k)")] = rtc === nothing ? -99999 : rtc.tap
        r[Symbol("phase_tap_position$(k)")] = ptc === nothing ? -99999 : ptc.tap
        r[Symbol("p$(k)")] = _iidm_float(a, "p$(k)")
        r[Symbol("q$(k)")] = _iidm_float(a, "q$(k)")
        r[Symbol("rho$(k)")] = u0 / uk * _iidm_step_rho(rtc) * _iidm_step_rho(ptc)
        r[Symbol("alpha$(k)")] = ptc === nothing ? 0.0 : _iidm_float(ptc.step, "alpha", 0.0)
        r[Symbol("selected_limits_group_$(k)")] = _iidm_selected_group(a, string(k))
        _iidm_side!(r, e.terminals[k], k)
        _iidm_tap_changer_rows!(net, T["ratio_tap_changers"], T["ratio_tap_changer_steps"], T["phase_tap_changers"], T["phase_tap_changer_steps"], id, _IIDM_SIDES[k], e.terminals[k], rtc, ptc)
        _iidm_limit_rows!(T["operational_limits"], e.children, id, _IIDM_SIDES[k], string(k), "THREE_WINDINGS_TRANSFORMER", _iidm_str(a, "selectedOperationalLimitsGroupIds$(k)"))
      end
      push!(T["3_windings_transformers"], r)
    elseif e.kind == :generator
      t = e.terminals[1]
      lim = _iidm_reactive_limits(e.children)
      target_p = _iidm_float(a, "targetP")
      p = _iidm_float(a, "p")
      lo_t, hi_t = _iidm_limits_at(lim, target_p)
      lo_p, hi_p = _iidm_limits_at(lim, -p)
      reg, _ = _iidm_ref_terminal(net, e.children, t, ("regulatingTerminal",))
      r = _IidmRow(:id => id, :name => name, :energy_source => _iidm_str(a, "energySource", "OTHER"), :target_p => target_p, :min_p => _iidm_float(a, "minP"), :max_p => _iidm_float(a, "maxP"), :min_q => lim.min_q, :max_q => lim.max_q, :min_q_at_target_p => lo_t, :max_q_at_target_p => hi_t, :min_q_at_p => lo_p, :max_q_at_p => hi_p, :rated_s => _iidm_float(a, "ratedS"), :reactive_limits_kind => lim.kind, :target_v => _iidm_float(a, "targetV"), :equivalent_local_target_v => NaN, :target_q => _iidm_float(a, "targetQ"), :voltage_regulator_on => _iidm_bool(a, "voltageRegulatorOn"), :regulated_element_id => reg.equipment, :regulated_bus_id => reg.bus_id, :regulated_bus_breaker_bus_id => reg.bb_bus, :p => p, :q => _iidm_float(a, "q"), :i => NaN, :fictitious => fict, :condenser => _iidm_bool(a, "condenser"))
      _iidm_attachment!(r, t)
      push!(T["generators"], r)
      for (k, pt) in enumerate(lim.points)
        push!(T["reactive_capability_curve_points"], _IidmRow(:id => id, :num => k - 1, :p => pt[1], :min_q => pt[2], :max_q => pt[3]))
      end
    elseif e.kind == :load
      r = _IidmRow(:id => id, :name => name, :type => _iidm_str(a, "loadType", "UNDEFINED"), :p0 => _iidm_float(a, "p0"), :q0 => _iidm_float(a, "q0"), :p => _iidm_float(a, "p"), :q => _iidm_float(a, "q"), :i => NaN, :fictitious => fict)
      _iidm_attachment!(r, e.terminals[1])
      push!(T["loads"], r)
    elseif e.kind == :shunt
      t = e.terminals[1]
      sections = _iidm_int(a, "sectionCount", _iidm_int(a, "currentSectionCount", 0))
      model = "LINEAR"
      g = 0.0
      b = 0.0
      max_sections = _iidm_int(a, "maximumSectionCount", 0)
      if haskey(a, "bPerSection")
        # IIDM before 1.3: the linear model on the element itself
        g = _iidm_float(a, "gPerSection", 0.0) * sections
        b = _iidm_float(a, "bPerSection", 0.0) * sections
        push!(T["linear_shunt_compensator_sections"], _IidmRow(:id => id, :g_per_section => _iidm_float(a, "gPerSection"), :b_per_section => _iidm_float(a, "bPerSection"), :max_section_count => max_sections))
      end
      for (n, ma, grand) in e.children
        if n == "shuntLinearModel"
          max_sections = _iidm_int(ma, "maximumSectionCount", 0)
          g = _iidm_float(ma, "gPerSection", 0.0) * sections
          b = _iidm_float(ma, "bPerSection", 0.0) * sections
          # a missing gPerSection is NaN in the sections table (as pypowsybl
          # reports it) and no conductance in the shunt's g
          push!(T["linear_shunt_compensator_sections"], _IidmRow(:id => id, :g_per_section => _iidm_float(ma, "gPerSection"), :b_per_section => _iidm_float(ma, "bPerSection"), :max_section_count => max_sections))
        elseif n == "shuntNonLinearModel"
          model = "NON_LINEAR"
          secs = [sa for (sn, sa, _) in grand if sn == "section"]
          max_sections = length(secs)
          for k in 1:min(sections, length(secs))
            g += _iidm_float(secs[k], "g", 0.0)
            b += _iidm_float(secs[k], "b", 0.0)
          end
        end
      end
      reg, _ = _iidm_ref_terminal(net, e.children, t, ("regulatingTerminal",))
      r = _IidmRow(:id => id, :name => name, :g => g, :b => b, :model_type => model, :max_section_count => max_sections, :section_count => sections, :solved_section_count => Float64(_iidm_int(a, "solvedSectionCount", sections)), :voltage_regulation_on => _iidm_bool(a, "voltageRegulatorOn"), :target_v => _iidm_float(a, "targetV"), :target_deadband => _iidm_float(a, "targetDeadband"), :regulating_bus_id => reg.bus_id, :p => _iidm_float(a, "p"), :q => _iidm_float(a, "q"), :i => NaN, :fictitious => fict)
      _iidm_attachment!(r, t)
      push!(T["shunt_compensators"], r)
    elseif e.kind == :svc
      t = e.terminals[1]
      mode = _iidm_str(a, "regulationMode", "OFF")
      reg, _ = _iidm_ref_terminal(net, e.children, t, ("regulatingTerminal",))
      r = _IidmRow(:id => id, :name => name, :b_min => _iidm_float(a, "bMin"), :b_max => _iidm_float(a, "bMax"), :target_v => _iidm_float(a, "voltageSetpoint"), :target_q => _iidm_float(a, "reactivePowerSetpoint"), :regulation_mode => mode, :regulating => _iidm_bool(a, "regulating", mode != "OFF"), :regulated_element_id => reg.equipment, :regulated_bus_id => reg.bus_id, :regulated_bus_breaker_bus_id => reg.bb_bus, :p => _iidm_float(a, "p"), :q => _iidm_float(a, "q"), :i => NaN, :fictitious => fict)
      _iidm_attachment!(r, t)
      push!(T["static_var_compensators"], r)
    elseif e.kind == :vsc
      t = e.terminals[1]
      lim = _iidm_reactive_limits(e.children)
      p = _iidm_float(a, "p")
      # the station's active setpoint is the HVDC line's; the bounds at the
      # setpoint are taken at the station's own p (NaN gives the plain bounds)
      lo_p, hi_p = _iidm_limits_at(lim, -p)
      reg, _ = _iidm_ref_terminal(net, e.children, t, ("regulatingTerminal",))
      r = _IidmRow(:id => id, :name => name, :loss_factor => _iidm_float32(a, "lossFactor"), :min_q => lim.min_q, :max_q => lim.max_q, :min_q_at_target_p => lo_p, :max_q_at_target_p => hi_p, :min_q_at_p => lo_p, :max_q_at_p => hi_p, :reactive_limits_kind => lim.kind, :target_v => _iidm_float(a, "voltageSetpoint"), :target_q => _iidm_float(a, "reactivePowerSetpoint"), :voltage_regulator_on => _iidm_bool(a, "voltageRegulatorOn"), :regulated_element_id => reg.equipment, :regulated_bus_id => reg.bus_id, :regulated_bus_breaker_bus_id => reg.bb_bus, :p => p, :q => _iidm_float(a, "q"), :i => NaN, :fictitious => fict, :hvdc_line_id => get(station_line, id, ""))
      _iidm_attachment!(r, t)
      push!(T["vsc_converter_stations"], r)
      for (k, pt) in enumerate(lim.points)
        push!(T["reactive_capability_curve_points"], _IidmRow(:id => id, :num => k - 1, :p => pt[1], :min_q => pt[2], :max_q => pt[3]))
      end
    elseif e.kind == :lcc
      r = _IidmRow(:id => id, :name => name, :power_factor => _iidm_float32(a, "powerFactor"), :loss_factor => _iidm_float32(a, "lossFactor"), :p => _iidm_float(a, "p"), :q => _iidm_float(a, "q"), :i => NaN, :fictitious => fict, :hvdc_line_id => get(station_line, id, ""))
      _iidm_attachment!(r, e.terminals[1])
      push!(T["lcc_converter_stations"], r)
    elseif e.kind == :dangling_line
      t = e.terminals[1]
      key = _iidm_str(a, "pairingKey", _iidm_str(a, "ucteXnodeCode"))
      tie = get(tie_of, id, "")
      r = _IidmRow(:id => id, :name => name, :r => _iidm_float(a, "r"), :x => _iidm_float(a, "x"), :g => _iidm_float(a, "g", 0.0), :b => _iidm_float(a, "b", 0.0), :p0 => _iidm_float(a, "p0"), :q0 => _iidm_float(a, "q0"), :p => _iidm_float(a, "p"), :q => _iidm_float(a, "q"), :i => NaN, :pairing_key => key, :ucte_xnode_code => key, :paired => !isempty(tie), :fictitious => fict, :tie_line_id => tie, :selected_limits_group => _iidm_selected_group(a, ""))
      _iidm_attachment!(r, t)
      push!(T["dangling_lines"], r)
      _iidm_limit_rows!(T["operational_limits"], e.children, id, "NONE", "", "BOUNDARY_LINE", _iidm_str(a, "selectedOperationalLimitsGroupIds"))
    elseif e.kind == :tie_line
      h1 = _iidm_tie_half(a, 1)
      h2 = _iidm_tie_half(a, 2)
      d1 = get(dangling, h1, nothing)
      d2 = get(dangling, h2, nothing)
      (d1 === nothing || d2 === nothing) && throw(ArgumentError("IIDM: tie line $(id) names dangling line $(d1 === nothing ? h1 : h2), which the file does not declare"))
      key = _iidm_str(d1.attrs, "pairingKey", _iidm_str(d1.attrs, "ucteXnodeCode"))
      push!(T["tie_lines"], _IidmRow(:id => id, :name => name, :boundary_line1_id => h1, :dangling_line1_id => h1, :boundary_line2_id => h2, :dangling_line2_id => h2, :pairing_key => key, :ucte_xnode_code => key, :connected1 => d1.terminals[1].connected, :connected2 => d2.terminals[1].connected, :fictitious => fict))
    elseif e.kind == :hvdc
      s1 = _iidm_str(a, "converterStation1")
      s2 = _iidm_str(a, "converterStation2")
      t1 = get(net.terminal_index, (s1, ""), nothing)
      t2 = get(net.terminal_index, (s2, ""), nothing)
      push!(T["hvdc_lines"], _IidmRow(:id => id, :name => name, :converters_mode => _iidm_str(a, "convertersMode"), :target_p => _iidm_float(a, "activePowerSetpoint"), :max_p => _iidm_float(a, "maxP"), :nominal_v => _iidm_float(a, "nominalV"), :r => _iidm_float(a, "r"), :converter_station1_id => s1, :converter_station2_id => s2, :connected1 => t1 !== nothing && t1.connected, :connected2 => t2 !== nothing && t2.connected, :fictitious => fict))
    end
  end
  # The bus state of the file belongs to the positions it was solved at. A
  # tap changer or shunt whose solved position differs from the position
  # the tables carry (micro_grid_be: solvedTapPosition 18 against
  # tapPosition 14) makes that state a start point for another network:
  # from it the Newton iteration diverged and only the rescue ladder
  # reached the solution. The state is dropped then, the import starts
  # flat, and the note travels in the manifest for the import report.
  stale = 0
  for name in ("ratio_tap_changers", "phase_tap_changers")
    for r in T[name]
      (isfinite(r[:solved_tap_position]) && Int(r[:solved_tap_position]) != r[:tap]) && (stale += 1)
    end
  end
  for r in T["shunt_compensators"]
    (isfinite(r[:solved_section_count]) && Int(r[:solved_section_count]) != r[:section_count]) && (stale += 1)
  end
  state_note = ""
  if stale > 0
    state_note = "bus state of the file not used as the start: $(stale) tap changer(s) or shunt(s) were solved at another position than the one imported (solvedTapPosition, solvedSectionCount); the import starts flat"
    for name in ("buses", "bus_breaker_view_buses")
      for r in T[name]
        r[:v_mag] = NaN
        r[:v_angle] = NaN
      end
    end
    for r in T["busbar_sections"]
      r[:v] = NaN
      r[:angle] = NaN
    end
  end
  manifest = Dict{String,Any}(
    "format" => POWSYBL_BUNDLE_FORMAT,
    "format_version" => POWSYBL_BUNDLE_FORMAT_VERSION,
    "case" => net.case,
    "source" => "Sparlectra IIDM reader",
    "source_file" => basename(net.path),
    "network_id" => _iidm_str(net.attrs, "id"),
    "case_date" => _iidm_str(net.attrs, "caseDate"),
    "iidm_version" => net.version,
    "pypowsybl_version" => "",
    "all_attributes" => true,
    "state_note" => state_note,
    "tables" => Dict{String,Any}(name => Dict{String,Any}("rows" => length(T[name])) for name in POWSYBL_TABLE_NAMES),
  )
  kwargs = Dict{Symbol,Any}(POWSYBL_TABLE_FIELDS[name] => _iidm_table(name, T[name]) for name in POWSYBL_TABLE_NAMES)
  return PowsyblTables(; manifest = manifest, kwargs...)
end

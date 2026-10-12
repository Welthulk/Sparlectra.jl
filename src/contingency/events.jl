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
# file: src/contingency/events.jl
# purpose: correlated contingency events (#355): groups of elements that trip
#          together, derived from the topology (double circuits, busbar
#          sections, stations, couplers, file name groups), and N-2 branch
#          pairs. Every generator returns ContingencyCase lists that run on
#          the same engine path as N-1 (one outage item with several
#          members).

# The branches an event may contain: closed (in service, both terminals),
# transformers only with include_transformers, never a leg of a
# three-winding transformer (the transformer is one case of its own).
function _event_branches(net::Net; include_transformers::Bool)::Vector{Int}
  legs = Set{Int}()
  for grp in _three_winding_groups(net)
    union!(legs, grp.legs)
  end
  out = Int[]
  for (k, br) in enumerate(net.branchVec)
    (br.status == 1 && _branch_terminal_state(br) === :closed) || continue
    k in legs && continue
    include_transformers || !_is_transformer_branch(br) || continue
    push!(out, k)
  end
  return out
end

# the higher endpoint voltage of a branch, the level it belongs to
_event_branch_vn(net::Net, br::Branch)::Float64 = max(getNodeVn(net.nodeVec[Int(br.fromBus)]), getNodeVn(net.nodeVec[Int(br.toBus)]))

# one generated event before the filters: kind, name, members as
# (class, index), and the voltage level for min_vn_kV
const _EventDraft = NamedTuple{(:kind, :name, :members, :vn_kV),Tuple{Symbol,String,Vector{Tuple{Symbol,Int}},Float64}}

function _event_cases(drafts::Vector{_EventDraft}, names::ContingencyNames; min_vn_kV::Real, max_members::Int, name_pattern)::Vector{ContingencyCase}
  out = ContingencyCase[]
  for d in drafts
    length(d.members) <= max_members || continue
    d.vn_kV >= min_vn_kV || continue
    (name_pattern === nothing || occursin(name_pattern, d.name)) || continue
    push!(out, ContingencyCase(d.name, d.kind, _contingency_member_elements(names, d.members)))
  end
  return out
end

const _EVENT_DOCS_FILTERS = """
Common keyword arguments:
- `names`: the display names ([`contingency_names`](@ref)); events name
  their members by these element strings and are named after buses.
- `include_transformers` (default `true`): with `false`, transformer
  branches are no members of any event.
- `min_vn_kV` (default `0.0`): keep an event whose highest member voltage
  meets this value.
- `name_pattern` (default `nothing`): keep an event whose name matches this
  `Regex` or contains this substring.
- `max_members` (default `6`): skip an event with more members.

Only closed branches (in service, both terminals closed) are members; the
legs of a three-winding transformer never are (the transformer is one
case of its own). An event with a single branch is no event (that is the
N-1 case) and is skipped.
"""

"""
    generateDoubleCircuits(net; names = contingency_names(net), include_transformers = true,
                           min_vn_kV = 0.0, name_pattern = nothing, max_members = 6)
        -> Vector{ContingencyCase}

One event of kind `:double_circuit` per pair of buses joined by two or more
closed branches: the circuits on the same towers, which a single fault can
take out together. Named `DC <from bus>-<to bus>`.

$(_EVENT_DOCS_FILTERS)
"""
function generateDoubleCircuits(net::Net; names::ContingencyNames = contingency_names(net), include_transformers::Bool = true,
  min_vn_kV::Real = 0.0, name_pattern::Union{Nothing,AbstractString,Regex} = nothing, max_members::Int = 6)::Vector{ContingencyCase}
  groups = Dict{Tuple{Int,Int},Vector{Int}}()
  order = Tuple{Int,Int}[]
  for k in _event_branches(net; include_transformers = include_transformers)
    br = net.branchVec[k]
    key = minmax(Int(br.fromBus), Int(br.toBus))
    haskey(groups, key) || push!(order, key)
    push!(get!(groups, key, Int[]), k)
  end
  drafts = _EventDraft[]
  for key in order
    ks = groups[key]
    length(ks) >= 2 || continue
    br = net.branchVec[ks[1]]
    f, t = Int(br.fromBus), Int(br.toBus)
    push!(drafts, (kind = :double_circuit, name = string("DC ", names.bus[f], "-", names.bus[t]), members = Tuple{Symbol,Int}[(:branch, k) for k in ks], vn_kV = _event_branch_vn(net, br)))
  end
  return _event_cases(drafts, names; min_vn_kV = min_vn_kV, max_members = max_members, name_pattern = name_pattern)
end

# the closed links per bus (index into linkVec)
function _closed_links_by_bus(net::Net)::Dict{Int,Vector{Int}}
  out = Dict{Int,Vector{Int}}()
  for (l, lk) in enumerate(net.linkVec)
    lk.status == 1 || continue
    push!(get!(out, Int(lk.fromBus), Int[]), l)
    push!(get!(out, Int(lk.toBus), Int[]), l)
  end
  return out
end

# the eligible branches per bus
function _event_branches_by_bus(net::Net, eligible::Vector{Int})::Dict{Int,Vector{Int}}
  out = Dict{Int,Vector{Int}}()
  for k in eligible
    br = net.branchVec[k]
    push!(get!(out, Int(br.fromBus), Int[]), k)
    Int(br.toBus) == Int(br.fromBus) || push!(get!(out, Int(br.toBus), Int[]), k)
  end
  return out
end

"""
    generateBusbarEvents(net; names = contingency_names(net), include_transformers = true,
                         min_vn_kV = 0.0, name_pattern = nothing, max_members = 6)
        -> Vector{ContingencyCase}

One event of kind `:busbar` per bus: a fault on one busbar section opens
its closed branches and its closed couplers (links), and the other
sections of the station stay in service. Loads, units and shunts at the bus
stay in the model; the run reports them as cut off (shed load, stranded
units, or an island with its own reference), as for a radial outage.
Named `BB <bus>`. Auxiliary buses (three-winding star points, internal
buses of external grids) are no busbars.

$(_EVENT_DOCS_FILTERS)
"""
function generateBusbarEvents(net::Net; names::ContingencyNames = contingency_names(net), include_transformers::Bool = true,
  min_vn_kV::Real = 0.0, name_pattern::Union{Nothing,AbstractString,Regex} = nothing, max_members::Int = 6)::Vector{ContingencyCase}
  by_bus = _event_branches_by_bus(net, _event_branches(net; include_transformers = include_transformers))
  links = _closed_links_by_bus(net)
  drafts = _EventDraft[]
  for b in eachindex(net.nodeVec)
    _is_aux_node(net.nodeVec[b]) && continue
    ks = get(by_bus, b, Int[])
    ls = get(links, b, Int[])
    (length(ks) + length(ls) >= 2 && !isempty(ks)) || continue
    members = vcat(Tuple{Symbol,Int}[(:branch, k) for k in ks], Tuple{Symbol,Int}[(:link, l) for l in ls])
    push!(drafts, (kind = :busbar, name = string("BB ", names.bus[b]), members = members, vn_kV = getNodeVn(net.nodeVec[b])))
  end
  return _event_cases(drafts, names; min_vn_kV = min_vn_kV, max_members = max_members, name_pattern = name_pattern)
end

"""
    generateSubstationEvents(net; names = contingency_names(net), include_transformers = true,
                             min_vn_kV = 0.0, name_pattern = nothing, max_members = 6)
        -> Vector{ContingencyCase}

One event of kind `:substation` per station: every closed branch incident
to any of its buses opens, so the whole station is out. A station is the
set of buses joined by closed links (busbar sections and their couplers);
a bus without links is a station of its own. Bus elements stay in the model
and are reported as cut off. Named `SUB <station bus>` after the lowest
bus index of the station.

$(_EVENT_DOCS_FILTERS)
"""
function generateSubstationEvents(net::Net; names::ContingencyNames = contingency_names(net), include_transformers::Bool = true,
  min_vn_kV::Real = 0.0, name_pattern::Union{Nothing,AbstractString,Regex} = nothing, max_members::Int = 6)::Vector{ContingencyCase}
  by_bus = _event_branches_by_bus(net, _event_branches(net; include_transformers = include_transformers))
  _, stations = _station_map(net)
  drafts = _EventDraft[]
  for rep in sort!(collect(keys(stations)))
    buses = sort(stations[rep])
    any(b -> _is_aux_node(net.nodeVec[b]), buses) && continue
    ks = sort!(unique(vcat((get(by_bus, b, Int[]) for b in buses)...)))
    length(ks) >= 2 || continue
    push!(drafts, (kind = :substation, name = string("SUB ", names.bus[rep]), members = Tuple{Symbol,Int}[(:branch, k) for k in ks], vn_kV = maximum(b -> getNodeVn(net.nodeVec[b]), buses)))
  end
  return _event_cases(drafts, names; min_vn_kV = min_vn_kV, max_members = max_members, name_pattern = name_pattern)
end

"""
    generateLinkOutages(net; names = contingency_names(net)) -> Vector{ContingencyCase}

One case of kind `:link` per closed bus link (a busbar coupler): the
coupler opens alone, and the two sections it joined run as separate buses.
Named after the link.
"""
function generateLinkOutages(net::Net; names::ContingencyNames = contingency_names(net))::Vector{ContingencyCase}
  return ContingencyCase[ContingencyCase(names.link[l], :link, names.link[l]) for (l, lk) in enumerate(net.linkVec) if lk.status == 1]
end

const DEFAULT_NAME_STEM_PATTERN = ContingencyEventsConfig().name_stem_pattern

"""
    generateNameGroups(net; names = contingency_names(net), stem = Regex(DEFAULT_NAME_STEM_PATTERN),
                       include_transformers = true, min_vn_kV = 0.0, name_pattern = nothing, max_members = 6)
        -> Vector{ContingencyCase}

One event of kind `:name_group` per set of branches whose display names
share a stem: the name with the circuit suffix `stem` matches removed
(`Ostheim_Neuwiese_1` and `_2`, or two legs carrying the same name). Only
branches with the same nominal voltage and the same pair of endpoint
stations group, so a voltage suffix (`Leitung_220`, `Leitung_380`) or a
name reused elsewhere never joins unrelated branches. Named `NG <stem>`.
The grouping needs the file names of an SCF case: with the generated names
of a plain net (`names.generated`) the list is empty and one line says so.

$(_EVENT_DOCS_FILTERS)
"""
function generateNameGroups(net::Net; names::ContingencyNames = contingency_names(net), stem::Regex = Regex(DEFAULT_NAME_STEM_PATTERN),
  include_transformers::Bool = true, min_vn_kV::Real = 0.0, name_pattern::Union{Nothing,AbstractString,Regex} = nothing, max_members::Int = 6)::Vector{ContingencyCase}
  if names.generated
    @info "generateNameGroups: no file names on this network (not an SCF case map), so no name groups." maxlog = 1
    return ContingencyCase[]
  end
  reps, _ = _station_map(net)
  groups = Dict{Tuple{String,Float64,Int,Int},Vector{Int}}()
  order = Tuple{String,Float64,Int,Int}[]
  for k in _event_branches(net; include_transformers = include_transformers)
    br = net.branchVec[k]
    s = replace(names.branch_raw[k], stem => "")
    isempty(s) && continue
    a, b = minmax(reps[Int(br.fromBus)], reps[Int(br.toBus)])
    key = (s, _event_branch_vn(net, br), a, b)
    haskey(groups, key) || push!(order, key)
    push!(get!(groups, key, Int[]), k)
  end
  drafts = _EventDraft[]
  for key in order
    ks = groups[key]
    length(ks) >= 2 || continue
    push!(drafts, (kind = :name_group, name = string("NG ", key[1]), members = Tuple{Symbol,Int}[(:branch, k) for k in ks], vn_kV = key[2]))
  end
  return _event_cases(drafts, names; min_vn_kV = min_vn_kV, max_members = max_members, name_pattern = name_pattern)
end

"""
    generateBranchPairs(net; names = contingency_names(net), rule = :all,
                        include_transformers = true, max_pairs = 50000)
        -> Vector{ContingencyCase}

N-2 pairs of kind `:pair`: with `rule = :all` every unordered pair of
closed branches (the full enumeration; three-winding legs excluded, the
transformer is one case), named `P <a> + <b>` with the element strings of
`names`. More than `max_pairs` pairs throw an `ArgumentError` before any
case is built, so a run never truncates silently. The pair count grows with
the square of the branch count (217 branches give 23436 pairs).
"""
function generateBranchPairs(net::Net; names::ContingencyNames = contingency_names(net), rule::Symbol = :all,
  include_transformers::Bool = true, max_pairs::Int = 50000)::Vector{ContingencyCase}
  rule in CONTINGENCY_PAIR_RULE_VALUES || throw(ArgumentError("generateBranchPairs: rule must be one of $(join(CONTINGENCY_PAIR_RULE_VALUES, ", ")), got :$(rule)."))
  ks = _event_branches(net; include_transformers = include_transformers)
  n = length(ks)
  n_pairs = div(n * (n - 1), 2)
  n_pairs <= max_pairs || throw(ArgumentError("generateBranchPairs: $(n) branches give $(n_pairs) pairs, above contingency.pairs.max_pairs = $(max_pairs). Raise the limit or narrow the branch list; nothing was solved."))
  out = ContingencyCase[]
  sizehint!(out, n_pairs)
  for i in 1:n, j in (i+1):n
    a, b = names.branch[ks[i]], names.branch[ks[j]]
    push!(out, ContingencyCase(string("P ", a, " + ", b), :pair, String[a, b]))
  end
  return out
end

# the member set of a case, for the duplicate removal across classes
function _event_member_key(net::Net, names::ContingencyNames, case::ContingencyCase)
  members = _contingency_case_members(net, names, case)
  members isa AbstractString && return nothing
  return Set(members)
end

"""
    generateCorrelatedEvents(net; names = contingency_names(net),
                             double_circuits = true, busbars = false, substations = false,
                             links = false, name_groups = false,
                             name_stem_pattern = DEFAULT_NAME_STEM_PATTERN, max_members = 6,
                             pairs = false, pair_rule = :all, max_pairs = 50000,
                             include_transformers = true, min_vn_kV = 0.0, name_pattern = nothing)
        -> Vector{ContingencyCase}

The union of the selected event classes in the order double circuits
([`generateDoubleCircuits`](@ref)), name groups
([`generateNameGroups`](@ref)), busbar sections
([`generateBusbarEvents`](@ref)), stations
([`generateSubstationEvents`](@ref)), couplers
([`generateLinkOutages`](@ref)) and, with `pairs`, N-2 pairs
([`generateBranchPairs`](@ref)). Events with the same member set appear
once, as the first class that produced them: a double circuit that is also
a name group stays `:double_circuit`, a single-bus station is `:busbar`
when busbars are on and `:substation` when only stations are, and a pair
that equals an event is dropped.
"""
function generateCorrelatedEvents(net::Net; names::ContingencyNames = contingency_names(net),
  double_circuits::Bool = true, busbars::Bool = false, substations::Bool = false, links::Bool = false, name_groups::Bool = false,
  name_stem_pattern::AbstractString = DEFAULT_NAME_STEM_PATTERN, max_members::Int = 6,
  pairs::Bool = false, pair_rule::Symbol = :all, max_pairs::Int = 50000,
  include_transformers::Bool = true, min_vn_kV::Real = 0.0, name_pattern::Union{Nothing,AbstractString,Regex} = nothing)::Vector{ContingencyCase}
  common = (; names = names, include_transformers = include_transformers, min_vn_kV = min_vn_kV, name_pattern = name_pattern, max_members = max_members)
  lists = Vector{ContingencyCase}[]
  double_circuits && push!(lists, generateDoubleCircuits(net; common...))
  name_groups && push!(lists, generateNameGroups(net; common..., stem = Regex(name_stem_pattern)))
  busbars && push!(lists, generateBusbarEvents(net; common...))
  substations && push!(lists, generateSubstationEvents(net; common...))
  links && push!(lists, generateLinkOutages(net; names = names))
  pairs && push!(lists, generateBranchPairs(net; names = names, rule = pair_rule, include_transformers = include_transformers, max_pairs = max_pairs))
  seen = Set{Set{Tuple{Symbol,Int}}}()
  out = ContingencyCase[]
  for list in lists, case in list
    key = _event_member_key(net, names, case)
    key === nothing && continue
    key in seen && continue
    push!(seen, key)
    push!(out, case)
  end
  return out
end

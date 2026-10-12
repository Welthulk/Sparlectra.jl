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
# file: src/contingency/names.jl
# purpose: display names of the outage elements of a network, the one source
#          of case names, result elements and overload names of contingency
#          output. A plain net uses its generated component names; an SCF
#          case supplies the file names of its `extra` block. The net's own
#          component names are never changed, so power flow, state
#          estimation and every export keep the generated names.

# the case kinds of contingency cases, scenarios and results: single
# elements, the correlated event classes, a user scenario, and N-2 pairs
const CONTINGENCY_SINGLE_KINDS = (:branch, :gen, :transformer3w, :link)
const CONTINGENCY_EVENT_KINDS = (:double_circuit, :busbar, :substation, :name_group)
const CONTINGENCY_KINDS = (CONTINGENCY_SINGLE_KINDS..., CONTINGENCY_EVENT_KINDS..., :explicit, :pair)

"""
    ContingencyNames

Display names of the outage elements of one network, by internal index.
`branch_raw`, `gen_raw` and `link_raw` hold the display name of every
branch, prosumer and bus link (the file name of an SCF case where the file
names the component, the generated component name otherwise); a prosumer
that is no generator carries an empty name. `branch`, `gen` and `link` hold
the element strings contingency output uses: the display name where it is
unique within its class, `"<name>#<index>"` where it is not (parallel
circuits share a generated name, and file names need not be unique either).
The `*_by_element` dictionaries resolve an element string back to its
index. `bus` is the bus name per bus index.

Build it with [`contingency_names`](@ref).
"""
struct ContingencyNames
  branch::Vector{String}
  branch_raw::Vector{String}
  gen::Vector{String}
  gen_raw::Vector{String}
  link::Vector{String}
  link_raw::Vector{String}
  bus::Vector{String}
  branch_by_element::Dict{String,Int}
  gen_by_element::Dict{String,Int}
  link_by_element::Dict{String,Int}
  # true when every name is the generated one (no file names applied)
  generated::Bool
end

# "<name>" where the name is unique among the non-empty entries, otherwise
# "<name>#<index>"; empty raw names stay empty
function _unique_element_names(raw::Vector{String})::Vector{String}
  count = Dict{String,Int}()
  for name in raw
    isempty(name) || (count[name] = get(count, name, 0) + 1)
  end
  return String[isempty(name) ? "" : (count[name] > 1 ? string(name, "#", k) : name) for (k, name) in enumerate(raw)]
end

function _contingency_names(net::Net, branch_raw::Vector{String}, gen_raw::Vector{String}, link_raw::Vector{String}, generated::Bool)::ContingencyNames
  branch = _unique_element_names(branch_raw)
  gen = _unique_element_names(gen_raw)
  link = _unique_element_names(link_raw)
  bus = fill("", length(net.nodeVec))
  for (name, idx) in net.busDict
    1 <= idx <= length(bus) && (bus[idx] = name)
  end
  for k in eachindex(bus)
    isempty(bus[k]) && (bus[k] = getCompName(net.nodeVec[k].comp))
  end
  by(v) = Dict{String,Int}(name => k for (k, name) in enumerate(v) if !isempty(name))
  return ContingencyNames(branch, branch_raw, gen, gen_raw, link, link_raw, bus, by(branch), by(gen), by(link), generated)
end

"""
    contingency_names(net) -> ContingencyNames
    contingency_names(net, index::ScenarioIndex) -> ContingencyNames
    contingency_names(net, case::SCFCase) -> ContingencyNames

The display names of the outage elements of `net`. Without a case the
generated component names are used (MATPOWER, CGMES and PowSyBl cases,
programmatic nets). With the [`ScenarioIndex`](@ref) or the `SCFCase` the
net was built from, every component the file names gets its file name
(`sparlectra.extra[<id>].name`); the rest keep the generated name. Element
strings that are not unique carry `#<index>`. Every case generator, both
batch entries and the service use this one map.
"""
function contingency_names(net::Net)::ContingencyNames
  branch_raw = String[getCompName(br.comp) for br in net.branchVec]
  gen_raw = String[isGenerator(ps) ? getCompName(ps.comp) : "" for ps in net.prosumpsVec]
  link_raw = String[lk.cName for lk in net.linkVec]
  return _contingency_names(net, branch_raw, gen_raw, link_raw, true)
end

# resolve an element string of one class: the exact element string, the
# "<name>#<index>" form whose index carries that display name, or a bare
# display name (the first element carrying it; the historical contract of a
# hand-written case on a shared name)
function _resolve_named_element(by_element::Dict{String,Int}, raw::Vector{String}, element::AbstractString)::Union{Nothing,Int}
  hit = get(by_element, element, nothing)
  hit === nothing || return hit
  hash_pos = findlast('#', element)
  if hash_pos !== nothing
    idx = tryparse(Int, element[(hash_pos+1):end])
    name = element[1:(hash_pos-1)]
    (idx !== nothing && 1 <= idx <= length(raw) && raw[idx] == name) && return idx
    return nothing
  end
  return findfirst(==(element), raw)
end

_resolve_named_branch(names::ContingencyNames, element::AbstractString) = _resolve_named_element(names.branch_by_element, names.branch_raw, element)
_resolve_named_generator(names::ContingencyNames, element::AbstractString) = _resolve_named_element(names.gen_by_element, names.gen_raw, element)
_resolve_named_link(names::ContingencyNames, element::AbstractString) = _resolve_named_element(names.link_by_element, names.link_raw, element)

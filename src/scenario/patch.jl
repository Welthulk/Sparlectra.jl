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

# file: src/scenario/patch.jl
# purpose: the scenario patch model: a scenario is an ordered list of patch operations on SCF
#          component ids, validated against a ScenarioIndex built from the
#          typed case. N-1 is the special case "one status patch per branch
#          or generator" and expands through the existing generators, so
#          the engine sees one kind of input.

# the closed operation list: op => (allowed targets, allowed fields)
const _PATCH_STATUS_TARGETS = (:branch, :transformer, :generator, :load, :shunt, :link)
const _PATCH_SET_TARGETS = (:generator, :load, :shunt, :transformer, :external_grid)
const _PATCH_SET_FIELDS = (:p, :q, :vm_pu, :tap_pos, :b_pu, :angle_deg)
const _PATCH_SCALE_TARGETS = (:load, :generator)

"""
    PatchOp

One patch operation of a scenario: `op` is `:status`,
`:set` or `:scale`, `target` names the component class, `id` the SCF
component id. `:status` carries `value` (0 or 1); `:set` carries `field`
(one of `p`, `q`, `vm_pu`, `tap_pos`, `b_pu`, `angle_deg`) and `value`;
`:scale` carries `factor` and applies it to `p` and `q` together unless
`field` names one of them. `angle_deg` is in degrees: the phase shift of
a transformer, or the voltage angle of an external grid (the SCF
`u_ref_angle`, which the file stores in radians); on the reference unit
it sets the angle reference of its island. No structural operations
exist by design.
"""
Base.@kwdef struct PatchOp
  op::Symbol
  target::Symbol
  id::Int
  field::Union{Nothing,Symbol} = nothing
  value::Union{Nothing,Float64} = nothing
  factor::Union{Nothing,Float64} = nothing
end

"""
    Scenario

A named, weighted, ordered list of patch operations with a `kind` from the
vocabulary of [`ContingencyCase`](@ref): `:explicit` for a user scenario
(the default), a single-element kind for an expanded N-1 outage, an event
class for a generated correlated event, `:pair` for an N-2 pair. The
positional constructor `Scenario(name, weight, ops)` sets `:explicit`.
"""
Base.@kwdef struct Scenario
  name::String
  weight::Float64 = 1.0
  ops::Vector{PatchOp} = PatchOp[]
  kind::Symbol = :explicit
end

Scenario(name::String, weight::Real, ops::Vector{PatchOp}) = Scenario(name, Float64(weight), ops, :explicit)

"""
    ScenarioSet

An ordered vector of scenarios plus the selection `mode` (`:explicit`,
`:n1_branches`, `:n1_generators`, `:n1_all`), the `exclusions` the N-1
expansion honors, and the correlated event classes `events` (any of
`:double_circuits`, `:busbars`, `:substations`, `:links`, `:name_groups`,
`:pairs`) generated next to the listed scenarios. The N-1 modes and the
event classes expand to explicit scenarios at load time
([`expand_scenarios`](@ref)), so the engine sees one kind of input.
"""
Base.@kwdef struct ScenarioSet
  scenarios::Vector{Scenario} = Scenario[]
  mode::Symbol = :explicit
  exclusions::Vector{String} = String[]
  events::Vector{Symbol} = Symbol[]
end

# the event classes a scenario set (document key `events`) and the
# configuration (contingency.events.*, contingency.pairs.*) can ask for
const SCENARIO_EVENT_CLASSES = (:double_circuits, :busbars, :substations, :links, :name_groups, :pairs)

"""
    ScenarioIndex

The id map of one typed case: SCF component id to component class and the
INTERNAL index the built network uses (branch vector position for
branches, prosumer position for appliances, shunt position for shunts,
link position for links), so a patch resolves without name lookups, plus
the file name of every component the case names (`extra[<id>].name`), the
source of the display names of contingency output
([`contingency_names`](@ref)).
"""
struct ScenarioIndex
  kind_by_id::Dict{Int,Symbol}
  internal_by_id::Dict{Int,Int}
  name_by_id::Dict{Int,String}
end

function _scenario_internal(extra::Dict{String,Any}, id::Int, key::String)
  e = get(extra, string(id), nothing)
  e isa AbstractDict || return nothing
  haskey(e, key) || return nothing
  return _scf_int(e[key], string("extra.", key))
end

"""
    ScenarioIndex(case::SCFCase) -> ScenarioIndex

Build the id map of `case`: lines and links map to `:branch` and `:link`,
generic branches to `:transformer`, sources to `:external_grid`, sym_gen
to `:generator`, sym_load to `:load`, shunts to `:shunt`. Internal indices
come from the recorded build-order fields of the extra block.
"""
function ScenarioIndex(case::SCFCase)
  extra = case.sparlectra === nothing ? Dict{String,Any}() : case.sparlectra.extra
  kinds = Dict{Int,Symbol}()
  internals = Dict{Int,Int}()
  d = case.data
  for r in d.line
    kinds[r.id] = :branch
    v = _scenario_internal(extra, r.id, "branch_index")
    v === nothing || (internals[r.id] = v)
  end
  for r in d.generic_branch
    kinds[r.id] = :transformer
    v = _scenario_internal(extra, r.id, "branch_index")
    v === nothing || (internals[r.id] = v)
  end
  # the importer adds the links sorted by id (scf_import.jl), so the
  # internal index is the rank in that order, not the file position
  for (k, r) in enumerate(sort(d.link; by = x -> x.id))
    kinds[r.id] = :link
    internals[r.id] = k
  end
  for (rows, kind) in ((d.source, :external_grid), (d.sym_gen, :generator), (d.sym_load, :load))
    for r in rows
      kinds[r.id] = kind
      v = _scenario_internal(extra, r.id, "element_index")
      v === nothing || (internals[r.id] = v)
    end
  end
  for r in d.shunt
    kinds[r.id] = :shunt
    v = _scenario_internal(extra, r.id, "shunt_index")
    v === nothing || (internals[r.id] = v)
  end
  names = Dict{Int,String}()
  for (id_text, entry) in extra
    entry isa AbstractDict || continue
    name = get(entry, "name", nothing)
    name isa AbstractString && !isempty(name) || continue
    id = tryparse(Int, String(id_text))
    id === nothing || (names[id] = String(name))
  end
  return ScenarioIndex(kinds, internals, names)
end

function contingency_names(net::Net, index::ScenarioIndex)::ContingencyNames
  generated = contingency_names(net)
  branch_raw = copy(generated.branch_raw)
  gen_raw = copy(generated.gen_raw)
  link_raw = copy(generated.link_raw)
  applied = false
  for (id, kind) in index.kind_by_id
    name = get(index.name_by_id, id, nothing)
    internal = get(index.internal_by_id, id, nothing)
    (name === nothing || internal === nothing) && continue
    if kind === :branch || kind === :transformer
      1 <= internal <= length(branch_raw) || continue
      branch_raw[internal] = name
    elseif kind === :generator || kind === :external_grid
      (1 <= internal <= length(gen_raw) && !isempty(gen_raw[internal])) || continue
      gen_raw[internal] = name
    elseif kind === :link
      1 <= internal <= length(link_raw) || continue
      link_raw[internal] = name
    else
      continue
    end
    applied = true
  end
  return _contingency_names(net, branch_raw, gen_raw, link_raw, !applied)
end

contingency_names(net::Net, case::SCFCase)::ContingencyNames = contingency_names(net, ScenarioIndex(case))

_scenario_error(scenario::AbstractString, opidx::Int, msg::AbstractString) = throw(ArgumentError(string("scenario ", repr(String(scenario)), " op ", opidx, ": ", msg)))

"""
    validate_scenarios(set, index) -> ScenarioSet

Load-time validation of the patch model: every operation must name a
known id whose class admits the operation and field; a scenario needs a
non-empty op list and a unique name. Errors carry the scenario name and
the op index. Returns `set` unchanged on success.
"""
function validate_scenarios(set::ScenarioSet, index::ScenarioIndex)::ScenarioSet
  for cls in set.events
    cls in SCENARIO_EVENT_CLASSES || throw(ArgumentError(string("scenario events: unknown class ", cls, "; known: ", join(SCENARIO_EVENT_CLASSES, ", "))))
  end
  seen = Set{String}()
  for s in set.scenarios
    isempty(strip(s.name)) && throw(ArgumentError("scenario names must be non-empty"))
    s.name in seen && throw(ArgumentError(string("duplicate scenario name ", repr(s.name))))
    push!(seen, s.name)
    isempty(s.ops) && throw(ArgumentError(string("scenario ", repr(s.name), ": empty op list")))
    s.kind in CONTINGENCY_KINDS || throw(ArgumentError(string("scenario ", repr(s.name), ": kind must be one of ", join(CONTINGENCY_KINDS, ", "), ", got ", s.kind)))
    (isfinite(s.weight) && s.weight >= 0.0) || throw(ArgumentError(string("scenario ", repr(s.name), ": weight must be a finite value >= 0")))
    for (i, op) in enumerate(s.ops)
      kind = get(index.kind_by_id, op.id, nothing)
      kind === nothing && _scenario_error(s.name, i, string("unknown component id ", op.id))
      kind == op.target || _scenario_error(s.name, i, string("id ", op.id, " is a ", kind, ", not a ", op.target))
      if op.op === :status
        op.target in _PATCH_STATUS_TARGETS || _scenario_error(s.name, i, string("status does not apply to a ", op.target))
        op.value in (0.0, 1.0) || _scenario_error(s.name, i, "status value must be 0 or 1")
      elseif op.op === :set
        op.target in _PATCH_SET_TARGETS || _scenario_error(s.name, i, string("set does not apply to a ", op.target))
        op.field in _PATCH_SET_FIELDS || _scenario_error(s.name, i, string("set field must be one of ", join(_PATCH_SET_FIELDS, ", ")))
        op.value isa Float64 || _scenario_error(s.name, i, "set needs a numeric value")
        if op.field === :tap_pos
          op.target === :transformer || _scenario_error(s.name, i, "tap_pos applies to a transformer")
        elseif op.field === :b_pu
          op.target === :shunt || _scenario_error(s.name, i, "b_pu applies to a shunt")
        elseif op.field === :vm_pu
          op.target in (:generator, :external_grid) || _scenario_error(s.name, i, "vm_pu applies to a generator or the external grid")
        elseif op.field === :angle_deg
          op.target in (:transformer, :external_grid) || _scenario_error(s.name, i, "angle_deg applies to a transformer or the external grid")
        end
      elseif op.op === :scale
        op.target in _PATCH_SCALE_TARGETS || _scenario_error(s.name, i, string("scale does not apply to a ", op.target))
        op.factor isa Float64 && isfinite(op.factor) || _scenario_error(s.name, i, "scale needs a finite factor")
        op.field === nothing || op.field in (:p, :q) || _scenario_error(s.name, i, "scale field must be p or q when given")
      else
        _scenario_error(s.name, i, string("unknown op ", op.op))
      end
    end
  end
  return set
end

# SCF id per internal index of the classes an outage can address
function _scenario_ids_by_internal(index::ScenarioIndex)
  branch = Dict{Int,Int}()
  appliance = Dict{Int,Int}()
  link = Dict{Int,Int}()
  for (id, kind) in index.kind_by_id
    internal = get(index.internal_by_id, id, nothing)
    internal === nothing && continue
    if kind === :branch || kind === :transformer
      branch[internal] = id
    elseif kind === :generator || kind === :external_grid
      appliance[internal] = id
    elseif kind === :link
      link[internal] = id
    end
  end
  return (; branch, appliance, link)
end

"""
    scenario_from_case(case, net, index; names = contingency_names(net, index)) -> Union{Nothing,Scenario}

The scenario form of a [`ContingencyCase`](@ref): one status-0 operation
per outaged element (every leg of a three-winding transformer), with the
case's name, weight and kind. Returns `nothing` when an element does not
resolve to a component of the case (a case list built for another
network).
"""
function scenario_from_case(case::ContingencyCase, net::Net, index::ScenarioIndex; names::ContingencyNames = contingency_names(net, index), ids = _scenario_ids_by_internal(index))::Union{Nothing,Scenario}
  members = _contingency_case_members(net, names, case)
  members isa AbstractString && return nothing
  ops = PatchOp[]
  for (mkind, internal) in members
    id = mkind === :branch ? get(ids.branch, internal, nothing) : mkind === :gen ? get(ids.appliance, internal, nothing) : get(ids.link, internal, nothing)
    id === nothing && return nothing
    push!(ops, PatchOp(op = :status, target = index.kind_by_id[id], id = id, value = 0.0))
  end
  return Scenario(name = case.name, weight = case.weight, ops = ops, kind = case.kind)
end

"""
    expand_scenarios(set, net, index; events = ContingencyEventsConfig(), pairs = ContingencyPairsConfig()) -> Vector{Scenario}

Expand a scenario set to the scenarios the engine runs. The N-1 modes
expand through the EXISTING generators (`generateN1Branches`,
`generateN1Generators`), honoring `exclusions` by element string or case
name; every expanded scenario carries one status patch (a three-winding
transformer one per leg) and the kind of its case. The event classes of
`set.events` add the generated correlated events and N-2 pairs
([`generateCorrelatedEvents`](@ref), [`generateBranchPairs`](@ref)); the
generator options (`max_members`, `name_stem_pattern`, the pair rule and
cap) come from `events` and `pairs`. Explicit scenarios come first, then
the N-1 list, then the events. Case names and elements use the file names
of the case ([`contingency_names`](@ref)).
"""
function expand_scenarios(set::ScenarioSet, net::Net, index::ScenarioIndex; events = ContingencyEventsConfig(), pairs = ContingencyPairsConfig())::Vector{Scenario}
  set.mode in (:explicit, :n1_branches, :n1_generators, :n1_all) || throw(ArgumentError(string("unknown scenario mode ", set.mode)))
  (set.mode === :explicit && isempty(set.events)) && return set.scenarios
  names = contingency_names(net, index)
  ids = _scenario_ids_by_internal(index)
  excluded = Set(set.exclusions)
  out = copy(set.scenarios)
  cases = ContingencyCase[]
  set.mode in (:n1_branches, :n1_all) && append!(cases, generateN1Branches(net; names = names))
  set.mode in (:n1_generators, :n1_all) && append!(cases, generateN1Generators(net; names = names))
  if !isempty(set.events)
    classes = Set(set.events)
    append!(cases, generateCorrelatedEvents(net; names = names,
      double_circuits = :double_circuits in classes, busbars = :busbars in classes, substations = :substations in classes,
      links = :links in classes, name_groups = :name_groups in classes,
      max_members = events.max_members, name_stem_pattern = events.name_stem_pattern,
      pairs = :pairs in classes, pair_rule = pairs.rule, max_pairs = pairs.max_pairs))
  end
  for case in cases
    (case.name in excluded || any(in(excluded), case.elements)) && continue
    sc = scenario_from_case(case, net, index; names = names, ids = ids)
    sc === nothing || push!(out, sc)
  end
  return out
end

## --- the sparlectra.scenarios block ------------------------------------------

function _scenario_op_dict(op::PatchOp)::Dict{String,Any}
  d = Dict{String,Any}("op" => String(op.op), "target" => String(op.target), "id" => op.id)
  op.field === nothing || (d["field"] = String(op.field))
  op.value === nothing || (d["value"] = op.value)
  op.factor === nothing || (d["factor"] = op.factor)
  return d
end

# one scenario in document form; `kind` is written when it is not the
# default :explicit, so a file of user scenarios reads and writes back
# unchanged
function _scenario_dict(s::Scenario)::Dict{String,Any}
  d = Dict{String,Any}("name" => s.name, "weight" => s.weight, "ops" => Any[_scenario_op_dict(op) for op in s.ops])
  s.kind === :explicit || (d["kind"] = String(s.kind))
  return d
end

"""
    scenario_set_dict(set) -> Dict{String,Any}

The document form of a scenario set (`sparlectra.scenarios`).
"""
function scenario_set_dict(set::ScenarioSet)::Dict{String,Any}
  d = Dict{String,Any}("mode" => String(set.mode))
  isempty(set.exclusions) || (d["exclusions"] = copy(set.exclusions))
  isempty(set.events) || (d["events"] = Any[String(e) for e in set.events])
  isempty(set.scenarios) || (d["scenarios"] = Any[_scenario_dict(s) for s in set.scenarios])
  return d
end

function _scenario_from_dict(raw::AbstractDict)::Scenario
  ops = PatchOp[]
  for o in _scf_get(raw, "ops", [])
    o isa AbstractDict || throw(ArgumentError("every scenario op must be an object"))
    push!(ops, PatchOp(
      op = Symbol(String(_scf_require(o, "op", "a scenario op"))),
      target = Symbol(String(_scf_require(o, "target", "a scenario op"))),
      id = _scf_int(_scf_require(o, "id", "a scenario op"), "scenario op id"),
      field = haskey(o, "field") ? Symbol(String(o["field"])) : nothing,
      value = haskey(o, "value") ? _scf_num(o["value"], "scenario op value") : nothing,
      factor = haskey(o, "factor") ? _scf_num(o["factor"], "scenario op factor") : nothing,
    ))
  end
  return Scenario(name = String(_scf_require(raw, "name", "a scenario")), weight = haskey(raw, "weight") ? _scf_num(raw["weight"], "scenario weight") : 1.0, ops = ops, kind = Symbol(String(_scf_get(raw, "kind", "explicit"))))
end

"""
    scenario_set_from_dict(raw) -> ScenarioSet

Read a `sparlectra.scenarios` document block.
"""
function scenario_set_from_dict(raw::AbstractDict)::ScenarioSet
  return ScenarioSet(
    scenarios = Scenario[_scenario_from_dict(s) for s in _scf_get(raw, "scenarios", [])],
    mode = Symbol(String(_scf_get(raw, "mode", "explicit"))),
    exclusions = String[String(x) for x in _scf_get(raw, "exclusions", [])],
    events = Symbol[Symbol(String(x)) for x in _scf_get(raw, "events", [])],
  )
end

"""
    scf_case_scenarios(case) -> Union{Nothing,ScenarioSet}

The scenario set a typed case carries: the `scenarios` block when
present, otherwise `nothing`. Scenarios are a Sparlectra Case Format
feature; converting a case to another format does not carry them.
"""
function scf_case_scenarios(case::SCFCase)::Union{Nothing,ScenarioSet}
  spar = case.sparlectra
  spar === nothing && return nothing
  isempty(spar.scenarios) || return scenario_set_from_dict(spar.scenarios)
  return nothing
end

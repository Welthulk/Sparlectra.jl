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
# purpose: the scenario patch model (scenario task step 1, design decisions
#          D1/D2): a scenario is an ordered list of patch operations on SCF
#          component ids, validated against a ScenarioIndex built from the
#          typed case. N-1 is the special case "one status patch per branch
#          or generator" and expands through the existing generators, so
#          the engine sees one kind of input.

# the closed operation list of D1: op => (allowed targets, allowed fields)
const _PATCH_STATUS_TARGETS = (:branch, :transformer, :generator, :load, :shunt, :link)
const _PATCH_SET_TARGETS = (:generator, :load, :shunt, :transformer, :external_grid)
const _PATCH_SET_FIELDS = (:p, :q, :vm_pu, :tap_pos, :b_pu, :angle_deg)
const _PATCH_SCALE_TARGETS = (:load, :generator)

"""
    PatchOp

One patch operation of a scenario (design decision D1): `op` is `:status`,
`:set` or `:scale`, `target` names the component class, `id` the SCF
component id. `:status` carries `value` (0 or 1); `:set` carries `field`
(one of `p`, `q`, `vm_pu`, `tap_pos`, `b_pu`, `angle_deg`) and `value`;
`:scale` carries `factor` and applies it to `p` and `q` together unless
`field` names one of them. No structural operations exist by design.
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

A named, weighted, ordered list of patch operations (design decision D2).
"""
Base.@kwdef struct Scenario
  name::String
  weight::Float64 = 1.0
  ops::Vector{PatchOp} = PatchOp[]
end

"""
    ScenarioSet

An ordered vector of scenarios plus the selection `mode` (`:explicit`,
`:n1_branches`, `:n1_generators`, `:n1_all`) and the `exclusions` the N-1
expansion honors. The N-1 modes expand to explicit scenarios at load time
([`expand_scenarios`](@ref)), so the engine sees one kind of input.
"""
Base.@kwdef struct ScenarioSet
  scenarios::Vector{Scenario} = Scenario[]
  mode::Symbol = :explicit
  exclusions::Vector{String} = String[]
end

"""
    ScenarioIndex

The id map of one typed case: SCF component id to component class and the
INTERNAL index the built network uses (branch vector position for
branches, prosumer position for appliances, shunt position for shunts),
so a patch resolves without name lookups.
"""
struct ScenarioIndex
  kind_by_id::Dict{Int,Symbol}
  internal_by_id::Dict{Int,Int}
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
  for r in d.link
    kinds[r.id] = :link
    internals[r.id] = findfirst(x -> x.id == r.id, d.link)
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
  return ScenarioIndex(kinds, internals)
end

_scenario_error(scenario::AbstractString, opidx::Int, msg::AbstractString) = throw(ArgumentError(string("scenario ", repr(String(scenario)), " op ", opidx, ": ", msg)))

"""
    validate_scenarios(set, index) -> ScenarioSet

Load-time validation of design decision D1: every operation must name a
known id whose class admits the operation and field; a scenario needs a
non-empty op list and a unique name. Errors carry the scenario name and
the op index. Returns `set` unchanged on success.
"""
function validate_scenarios(set::ScenarioSet, index::ScenarioIndex)::ScenarioSet
  seen = Set{String}()
  for s in set.scenarios
    isempty(strip(s.name)) && throw(ArgumentError("scenario names must be non-empty"))
    s.name in seen && throw(ArgumentError(string("duplicate scenario name ", repr(s.name))))
    push!(seen, s.name)
    isempty(s.ops) && throw(ArgumentError(string("scenario ", repr(s.name), ": empty op list")))
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

"""
    expand_scenarios(set, net, index) -> Vector{Scenario}

Expand the N-1 modes to explicit scenarios through the EXISTING
generators (`generateN1Branches`, `generateN1Generators`), honoring
`exclusions` by case element name; `:explicit` returns the stored
scenarios unchanged. Every expanded scenario carries exactly one status
patch, which is what makes N-1 the special case of the model.
"""
function expand_scenarios(set::ScenarioSet, net::Net, index::ScenarioIndex)::Vector{Scenario}
  set.mode === :explicit && return set.scenarios
  set.mode in (:n1_branches, :n1_generators, :n1_all) || throw(ArgumentError(string("unknown scenario mode ", set.mode)))
  excluded = Set(set.exclusions)
  id_by_internal_branch = Dict{Int,Int}()
  id_by_internal_appliance = Dict{Int,Int}()
  for (id, kind) in index.kind_by_id
    internal = get(index.internal_by_id, id, nothing)
    internal === nothing && continue
    if kind === :branch || kind === :transformer
      id_by_internal_branch[internal] = id
    elseif kind === :generator || kind === :external_grid
      id_by_internal_appliance[internal] = id
    end
  end
  out = Scenario[]
  if set.mode in (:n1_branches, :n1_all)
    for case in generateN1Branches(net)
      case.element in excluded && continue
      internal = _resolve_contingency_branch(net, case.element)
      internal === nothing && continue
      id = get(id_by_internal_branch, internal, nothing)
      id === nothing && continue
      target = index.kind_by_id[id]
      push!(out, Scenario(name = case.name, weight = case.weight, ops = [PatchOp(op = :status, target = target, id = id, value = 0.0)]))
    end
  end
  if set.mode in (:n1_generators, :n1_all)
    for case in generateN1Generators(net)
      case.element in excluded && continue
      internal = _resolve_contingency_generator(net, case.element)
      internal === nothing && continue
      id = get(id_by_internal_appliance, internal, nothing)
      id === nothing && continue
      target = index.kind_by_id[id]
      push!(out, Scenario(name = case.name, weight = case.weight, ops = [PatchOp(op = :status, target = target, id = id, value = 0.0)]))
    end
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

"""
    scenario_set_dict(set) -> Dict{String,Any}

The document form of a scenario set (`sparlectra.scenarios`).
"""
function scenario_set_dict(set::ScenarioSet)::Dict{String,Any}
  d = Dict{String,Any}("mode" => String(set.mode))
  isempty(set.exclusions) || (d["exclusions"] = copy(set.exclusions))
  isempty(set.scenarios) || (d["scenarios"] = Any[Dict{String,Any}("name" => s.name, "weight" => s.weight, "ops" => Any[_scenario_op_dict(op) for op in s.ops]) for s in set.scenarios])
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
  return Scenario(name = String(_scf_require(raw, "name", "a scenario")), weight = haskey(raw, "weight") ? _scf_num(raw["weight"], "scenario weight") : 1.0, ops = ops)
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
  )
end

"""
    scenario_set_from_contingencies(raw, case) -> ScenarioSet

Map a legacy `sparlectra.contingencies` block onto the scenario model
(design decision D3): `mode` and `exclusions` carry over, and every
explicit case becomes one scenario with a status patch per outage
component id.
"""
function scenario_set_from_contingencies(raw::AbstractDict, case::SCFCase)::ScenarioSet
  index = ScenarioIndex(case)
  mode_raw = String(_scf_get(raw, "mode", "explicit"))
  mode = mode_raw == "all_branches" ? :n1_branches : mode_raw == "all_branches_plus" ? :n1_all : :explicit
  scenarios = Scenario[]
  for (k, c) in enumerate(_scf_get(raw, "cases", []))
    c isa AbstractDict || continue
    ops = PatchOp[]
    for o in _scf_get(c, "outages", [])
      id = _scf_int(_scf_require(o, "component", "a contingency outage"), "contingency outage component")
      target = get(index.kind_by_id, id, :branch)
      push!(ops, PatchOp(op = :status, target = target, id = id, value = 0.0))
    end
    push!(scenarios, Scenario(name = String(_scf_get(c, "name", string("contingency_", k))), weight = _scf_num(_scf_get(c, "weight", 1.0), "contingency weight"), ops = ops))
  end
  return ScenarioSet(scenarios = scenarios, mode = mode, exclusions = String[String(x) for x in _scf_get(raw, "exclusions", [])])
end

"""
    scf_case_scenarios(case) -> Union{Nothing,ScenarioSet}

The scenario set a typed case carries: the `scenarios` block when
present, otherwise the legacy `contingencies` block mapped onto the
scenario model, otherwise `nothing`.
"""
function scf_case_scenarios(case::SCFCase)::Union{Nothing,ScenarioSet}
  spar = case.sparlectra
  spar === nothing && return nothing
  isempty(spar.scenarios) || return scenario_set_from_dict(spar.scenarios)
  isempty(spar.contingencies) || return scenario_set_from_contingencies(spar.contingencies, case)
  return nothing
end

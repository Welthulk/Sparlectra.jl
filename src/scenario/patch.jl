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
#          the engine sees one kind of input. An index of a series (the
#          CGMES series union, origin :cgmes_series, or an SCF time series,
#          origin :scf_series) also accepts the QSTS vocabulary.

# the closed operation list: op => (allowed targets, allowed fields)
const _PATCH_STATUS_TARGETS = (:branch, :transformer, :generator, :load, :shunt, :link)
const _PATCH_SET_TARGETS = (:generator, :load, :shunt, :transformer, :external_grid)
const _PATCH_SET_FIELDS = (:p, :q, :vm_pu, :tap_pos, :b_pu, :angle_deg)
const _PATCH_SCALE_TARGETS = (:load, :generator)

# the QSTS vocabulary (stage 3 of the CGMES series): the `:set` fields a
# snapshot delta writes against the union template, field => (targets,
# value rule). Only a series index (origin in _QSTS_ORIGINS) accepts them,
# so the SCF scenario format stays unchanged. The names do not overlap the
# SCF set fields, so the field alone tells the vocabularies apart. r_pu to
# b_to_pu follow an EQ change between two exports; r_pu and x_pu may be
# negative (series capacitors, star legs of a three-winding transformer),
# apply! rejects r_pu = x_pu = 0.
const _QSTS_SET_RULES = Dict{Symbol,Tuple{Vector{Symbol},Symbol}}(
  :isolated => ([:bus], :flag),
  :vm_start => ([:bus], :finite),
  :va_start => ([:bus], :finite),
  :status => ([:branch, :transformer], :flag),
  :from_status => ([:branch, :transformer], :flag),
  :to_status => ([:branch, :transformer], :flag),
  :r_pu => ([:branch, :transformer], :finite),
  :x_pu => ([:branch, :transformer], :finite),
  :g_from_pu => ([:branch, :transformer], :finite),
  :b_from_pu => ([:branch, :transformer], :finite),
  :g_to_pu => ([:branch, :transformer], :finite),
  :b_to_pu => ([:branch, :transformer], :finite),
  :ratio => ([:transformer], :nonzero),
  :shift_deg => ([:branch, :transformer], :finite),
  :tap_data => ([:transformer], :record),
  :g_pu => ([:shunt], :finite),
  :in_service => ([:generator, :load, :external_grid], :flag),
  :qmin => ([:generator, :external_grid], :optional),
  :qmax => ([:generator, :external_grid], :optional),
  :regulated => ([:generator, :external_grid], :flag),
  :reference => ([:generator, :external_grid], :flag),
  :reference_priority => ([:generator, :external_grid], :count),
  :participation => ([:generator, :external_grid], :optional),
  :enabled => ([:tap_controller], :flag),
  :target_vm_pu => ([:tap_controller], :optional),
  :deadband_vm_pu => ([:tap_controller], :nonnegative),
  :p_target_mw => ([:tap_controller], :optional),
  :deadband_p_mw => ([:tap_controller], :nonnegative),
)

# the index origins that take the QSTS vocabulary and the three-pass apply!:
# the CGMES series union and the SCF time series
const _QSTS_ORIGINS = (:cgmes_series, :scf_series)

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

A [`ScenarioIndex`](@ref) with origin `:cgmes_series` or `:scf_series`
also accepts the QSTS `:set` fields: `isolated`, `vm_start`, `va_start` on
a bus; `status`, `from_status`, `to_status`, `r_pu`, `x_pu`, `g_from_pu`,
`b_from_pu`, `g_to_pu`, `b_to_pu` on a branch or transformer and `ratio`,
`shift_deg`, `tap_data` on a transformer; `g_pu` on a shunt; `in_service`
on a generator, load or external grid; `qmin`, `qmax`, `regulated`,
`reference`, `reference_priority`, `participation` on a generator or
external grid; `enabled`, `target_vm_pu`, `deadband_vm_pu`, `p_target_mw`,
`deadband_p_mw` on a tap controller. `qmin`, `qmax`, `participation`, `target_vm_pu` and
`p_target_mw` take `value = nothing` for "no value".

The branch parameters are per unit on the branch base, as the branch
stores them. `r_pu` and `x_pu` write the live and the equipment value
together, the four terminal values keep `g_pu` and `b_pu` their sums,
`ratio` writes `ratio` and `tap_ratio`, `shift_deg` writes `angle` and
`phase_shift_deg` (the SCF `angle_deg` writes `phase_shift_deg` only).
`tap_data` with value `k` writes record `k` of that transformer in
`index.tap_data` (a `BranchTapData`: tap range, step, table ratios,
impedance per step) onto the branch.
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
    BranchTapData

The tap changer data of one transformer branch as the importer sets it: the
ratio tap (`has_ratio_tap`, `tap_min`, `tap_max`, `tap_step`, the table
ratios `tap_ratio_points`), the phase tap (`has_phase_tap`, `phase_min_deg`,
`phase_max_deg`, `phase_step_deg`) and the impedance per tap step
(`tap_step_impedance`). The field names are those of `Branch`.
`BranchTapData(br)` takes the data of a branch. A series index keeps these
records per transformer id (`ScenarioIndex.tap_data`); the QSTS op
`:set :tap_data` with value `k` writes record `k` onto the branch.
"""
struct BranchTapData
  has_ratio_tap::Bool
  tap_min::Float64
  tap_max::Float64
  tap_step::Float64
  tap_ratio_points::Union{Nothing,Vector{Float64}}
  has_phase_tap::Bool
  phase_min_deg::Float64
  phase_max_deg::Float64
  phase_step_deg::Float64
  tap_step_impedance::Union{Nothing,TapStepImpedance}
end

# the ratio points are copied, the record must not follow later writes to
# the branch; TapStepImpedance is immutable and its vectors are never
# written after the import
BranchTapData(br::Branch) = BranchTapData(br.has_ratio_tap, br.tap_min, br.tap_max, br.tap_step, br.tap_ratio_points === nothing ? nothing : copy(br.tap_ratio_points),
                                          br.has_phase_tap, br.phase_min_deg, br.phase_max_deg, br.phase_step_deg, br.tap_step_impedance)

"""
    ScenarioIndex

The id map of one network: component id to component class (`kind_by_id`)
and the INTERNAL index the built network uses (`internal_by_id`: branch
vector position for branches and transformers, prosumer position for
appliances, shunt position for shunts, link position for links, bus number
for `:bus`, position in `_tap_controllers(net)` for `:tap_controller`), so
a patch resolves without name lookups, plus the display name of every
component (`name_by_id`, the SCF `extra[<id>].name`), the source of the
display names of contingency output ([`contingency_names`](@ref)).

Four fields serve the series (QSTS):
- `id_by_mrid`: mRID to component id (CGMES series union);
- `parked_values`: per component id, the values of a unit that is parked
  in the template, copied back by `:set :in_service 1`;
- `tap_data`: per transformer id, the tap changer data records
  (`BranchTapData`) that `:set :tap_data` writes, numbered from 1;
- `origin`: `:scf` for an index of a typed case, `:cgmes_series` for the
  CGMES series union, `:scf_series` for an SCF time series
  ([`SCFSeries`](@ref)). Only a series index accepts the QSTS vocabulary
  ([`validate_scenarios`](@ref)) and switches [`apply!`](@ref) to its
  three passes.

`ScenarioIndex(kind_by_id, internal_by_id, name_by_id)` builds an `:scf`
index with empty mRID, parked-value and tap data maps. The full
constructor takes `tap_data` as a keyword (default empty). Throws
`ArgumentError` for another origin.
"""
struct ScenarioIndex
  kind_by_id::Dict{Int,Symbol}
  internal_by_id::Dict{Int,Int}
  name_by_id::Dict{Int,String}
  id_by_mrid::Dict{String,Int}
  parked_values::Dict{Int,ProSumer}
  tap_data::Dict{Int,Vector{BranchTapData}}
  origin::Symbol
  function ScenarioIndex(kind_by_id::Dict{Int,Symbol}, internal_by_id::Dict{Int,Int}, name_by_id::Dict{Int,String}, id_by_mrid::Dict{String,Int}, parked_values::Dict{Int,ProSumer}, origin::Symbol;
                         tap_data::Dict{Int,Vector{BranchTapData}} = Dict{Int,Vector{BranchTapData}}())
    origin in (:scf, _QSTS_ORIGINS...) || throw(ArgumentError(string("ScenarioIndex: origin must be :scf, :cgmes_series or :scf_series, got :", origin)))
    return new(kind_by_id, internal_by_id, name_by_id, id_by_mrid, parked_values, tap_data, origin)
  end
end

# every construction from a typed case (and every caller from before the
# QSTS fields) is an SCF index without mRIDs or parked units
ScenarioIndex(kind_by_id::Dict{Int,Symbol}, internal_by_id::Dict{Int,Int}, name_by_id::Dict{Int,String}) = ScenarioIndex(kind_by_id, internal_by_id, name_by_id, Dict{String,Int}(), Dict{Int,ProSumer}(), :scf)

function _scenario_internal(extra::Dict{String,Any}, id::Int, key::String)
  e = get(extra, string(id), nothing)
  e isa AbstractDict || return nothing
  haskey(e, key) || return nothing
  return _scf_int(e[key], string("extra.", key))
end

"""
    ScenarioIndex(case::SCFCase; origin = :scf) -> ScenarioIndex

Build the id map of `case`: lines and links map to `:branch` and `:link`,
generic branches to `:transformer`, sources to `:external_grid`, sym_gen
to `:generator`, sym_load to `:load`, shunts to `:shunt`. Internal indices
come from the recorded build-order fields of the extra block.

`origin = :scf_series` builds the index of an SCF time series: the node
ids map to `:bus` as well (internal index from `extra.bus_index`), so the
QSTS bus fields can address them. The parked values stay empty; the
series fills them ([`read_scf_series`](@ref)).
"""
function ScenarioIndex(case::SCFCase; origin::Symbol = :scf)
  origin in (:scf, :scf_series) || throw(ArgumentError(string("ScenarioIndex(case): origin must be :scf or :scf_series, got :", origin)))
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
  if origin === :scf_series
    for r in d.node
      kinds[r.id] = :bus
      v = _scenario_internal(extra, r.id, "bus_index")
      v === nothing || (internals[r.id] = v)
    end
  end
  names = Dict{Int,String}()
  for (id_text, entry) in extra
    entry isa AbstractDict || continue
    name = get(entry, "name", nothing)
    name isa AbstractString && !isempty(name) || continue
    id = tryparse(Int, String(id_text))
    id === nothing || (names[id] = String(name))
  end
  origin === :scf && return ScenarioIndex(kinds, internals, names)
  return ScenarioIndex(kinds, internals, names, Dict{String,Int}(), Dict{Int,ProSumer}(), origin)
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

# the QSTS rule of an op, nothing for the SCF vocabulary
_scenario_qsts_rule(op::PatchOp)::Union{Nothing,Tuple{Vector{Symbol},Symbol}} = op.op === :set ? get(_QSTS_SET_RULES, op.field, nothing) : nothing

# what is wrong with a QSTS value under its rule, nothing when it is fine;
# :optional fields map to Union{Nothing,Float64} fields, so `nothing` is a value
function _scenario_qsts_value_error(rule::Symbol, v::Union{Nothing,Float64})::Union{Nothing,String}
  if rule === :optional
    return (v === nothing || isfinite(v)) ? nothing : "needs a finite value or none"
  end
  v === nothing && return "needs a numeric value"
  rule === :flag && return v in (0.0, 1.0) ? nothing : "value must be 0 or 1"
  isfinite(v) || return "needs a finite value"
  rule === :nonzero && v == 0.0 && return "value must not be 0"
  rule === :nonnegative && v < 0.0 && return "value must be >= 0"
  rule === :count && !(isinteger(v) && v >= 0.0) && return "value must be a whole number >= 0"
  rule === :record && !(isinteger(v) && v >= 1.0) && return "value must be a record number >= 1"
  return nothing
end

# The per-op part of validate_scenarios, in two halves so the caller can
# check something between them (validate_scenarios: the second write of a
# field). Each returns what is wrong with the op, nothing when it passes. The
# SCF series reader runs both and names the CSV line instead of the scenario.

# the class half: vocabulary, id and component class
function _patch_op_class_error(op::PatchOp, index::ScenarioIndex)::Union{Nothing,String}
  # checked before the id: an SCF index knows no bus or tap controller
  # ids, and the error must name the op, not an unknown id
  (_scenario_qsts_rule(op) === nothing || index.origin in _QSTS_ORIGINS) || return string("set ", op.field, " is QSTS vocabulary and needs a ScenarioIndex with origin cgmes_series or scf_series; this index has origin ", index.origin)
  kind = get(index.kind_by_id, op.id, nothing)
  kind === nothing && return string("unknown component id ", op.id)
  kind == op.target || return string("id ", op.id, " is a ", kind, ", not a ", op.target)
  return nothing
end

# the rule half: op, target, field and value
function _patch_op_rule_error(op::PatchOp, index::ScenarioIndex)::Union{Nothing,String}
  qsts_index = index.origin in _QSTS_ORIGINS
  qsts = _scenario_qsts_rule(op)
  if op.op === :status
    op.target in _PATCH_STATUS_TARGETS || return string("status does not apply to a ", op.target)
    op.value in (0.0, 1.0) || return "status value must be 0 or 1"
  elseif qsts !== nothing
    targets, rule = qsts
    op.target in targets || return string("set ", op.field, " applies to a ", join(targets, " or a "), ", not to a ", op.target)
    msg = _scenario_qsts_value_error(rule, op.value)
    msg === nothing || return string("set ", op.field, ": ", msg)
    if op.field === :in_service && op.value == 1.0
      haskey(index.parked_values, op.id) || return string("set in_service 1 needs the parked values of id ", op.id, " in the index")
    end
    if op.field === :tap_data
      n = length(get(index.tap_data, op.id, BranchTapData[]))
      op.value <= n || return string("set tap_data ", Int(op.value), " needs tap data record ", Int(op.value), " of id ", op.id, " in the index (it holds ", n, ")")
    end
  elseif op.op === :set
    op.target in _PATCH_SET_TARGETS || return string("set does not apply to a ", op.target)
    op.field in _PATCH_SET_FIELDS || return string("set field must be one of ", join(qsts_index ? vcat(collect(_PATCH_SET_FIELDS), sort!(collect(keys(_QSTS_SET_RULES)))) : collect(_PATCH_SET_FIELDS), ", "))
    op.value isa Float64 || return "set needs a numeric value"
    if op.field === :tap_pos
      op.target === :transformer || return "tap_pos applies to a transformer"
    elseif op.field === :b_pu
      op.target === :shunt || return "b_pu applies to a shunt"
    elseif op.field === :vm_pu
      op.target in (:generator, :external_grid) || return "vm_pu applies to a generator or the external grid"
    elseif op.field === :angle_deg
      op.target in (:transformer, :external_grid) || return "angle_deg applies to a transformer or the external grid"
    end
  elseif op.op === :scale
    op.target in _PATCH_SCALE_TARGETS || return string("scale does not apply to a ", op.target)
    op.factor isa Float64 && isfinite(op.factor) || return "scale needs a finite factor"
    op.field === nothing || op.field in (:p, :q) || return "scale field must be p or q when given"
  else
    return string("unknown op ", op.op)
  end
  return nothing
end

"""
    validate_scenarios(set, index) -> ScenarioSet

Load-time validation of the patch model: every operation must name a
known id whose class admits the operation and field; a scenario needs a
non-empty op list and a unique name. Errors carry the scenario name and
the op index. Returns `set` unchanged on success.

The QSTS `:set` fields (see [`PatchOp`](@ref)) are accepted only by a
series index (`index.origin` `:cgmes_series` or `:scf_series`); an `:scf`
index rejects them with the op named. With a series index a scenario may
write each (id, op, field) once, `in_service 1` needs the parked
values of its id in `index.parked_values`, and `tap_data k` record `k` of
its id in `index.tap_data`. The branch parameters need a
finite value, `ratio` a nonzero one. Whether `r_pu` and `x_pu` of a branch
are both 0 after the delta depends on the net, so [`apply!`](@ref) checks
that.
"""
function validate_scenarios(set::ScenarioSet, index::ScenarioIndex)::ScenarioSet
  for cls in set.events
    cls in SCENARIO_EVENT_CLASSES || throw(ArgumentError(string("scenario events: unknown class ", cls, "; known: ", join(SCENARIO_EVENT_CLASSES, ", "))))
  end
  qsts_index = index.origin in _QSTS_ORIGINS
  seen = Set{String}()
  for s in set.scenarios
    isempty(strip(s.name)) && throw(ArgumentError("scenario names must be non-empty"))
    s.name in seen && throw(ArgumentError(string("duplicate scenario name ", repr(s.name))))
    push!(seen, s.name)
    isempty(s.ops) && throw(ArgumentError(string("scenario ", repr(s.name), ": empty op list")))
    s.kind in CONTINGENCY_KINDS || throw(ArgumentError(string("scenario ", repr(s.name), ": kind must be one of ", join(CONTINGENCY_KINDS, ", "), ", got ", s.kind)))
    (isfinite(s.weight) && s.weight >= 0.0) || throw(ArgumentError(string("scenario ", repr(s.name), ": weight must be a finite value >= 0")))
    writes = Set{Tuple{Int,Symbol,Symbol}}()
    for (i, op) in enumerate(s.ops)
      msg = _patch_op_class_error(op, index)
      msg === nothing || _scenario_error(s.name, i, msg)
      if qsts_index
        # apply! sorts a QSTS delta into a fixed order, so only a second
        # write of the same field could make the op order matter
        key = (op.id, op.op, something(op.field, :none))
        key in writes && _scenario_error(s.name, i, string("id ", op.id, " gets a second ", op.op, op.field === nothing ? "" : string(" ", op.field), "; a QSTS delta writes each field once"))
        push!(writes, key)
      end
      msg = _patch_op_rule_error(op, index)
      msg === nothing || _scenario_error(s.name, i, msg)
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

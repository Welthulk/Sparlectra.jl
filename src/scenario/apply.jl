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

# file: src/scenario/apply.jl
# purpose: scenario task step 2 (design decision D4): apply! executes a
#          scenario's patch operations on a working copy with an undo log
#          of every change, restore! replays the log in reverse and
#          recomputes the derived state, so the copy returns to the base
#          BITWISE on every field a patch can touch. Status patches use
#          the same mutation functions the contingency runner uses, so the
#          engine's full run reproduces today's results.

# one recorded change: a plain field write, or a structural deletion from
# one of the net's vectors (the position makes the reversed replay exact)
struct UndoEntry
  kind::Symbol            # :field or :deleted
  obj::Any                # the object (field) or the vector (deleted)
  key::Any                # field Symbol, or the index the element lived at
  old::Any                # the previous value, or the deleted element
end

"""
    UndoLog

The reversible record of one [`apply!`](@ref): the field writes and
structural deletions in application order, plus the flag that derived
state (bus types, Q-limit tables, node injection sums, isolated-node
registry) must be recomputed after the replay.
"""
struct UndoLog
  entries::Vector{UndoEntry}
  refresh::Base.RefValue{Bool}
end

UndoLog() = UndoLog(UndoEntry[], Ref(false))

_undo_field!(log::UndoLog, obj, field::Symbol) = push!(log.entries, UndoEntry(:field, obj, field, getfield(obj, field)))

function _undo_set!(log::UndoLog, obj, field::Symbol, value)
  _undo_field!(log, obj, field)
  setfield!(obj, field, value)
  return nothing
end

# wrap a vector deletion so the reversed replay reinserts the exact object
# at the exact position
function _undo_deleteat!(log::UndoLog, vec::AbstractVector, idx::Int)
  push!(log.entries, UndoEntry(:deleted, vec, idx, vec[idx]))
  deleteat!(vec, idx)
  return nothing
end

function _scenario_branch!(log::UndoLog, net::Net, internal::Int, op::PatchOp)
  1 <= internal <= length(net.branchVec) || throw(ArgumentError(string("branch index ", internal, " out of range")))
  br = net.branchVec[internal]
  if op.op === :status
    _undo_set!(log, br, :status, Int(op.value))
    _undo_set!(log, br, :from_status, Int(op.value))
    _undo_set!(log, br, :to_status, Int(op.value))
    log.refresh[] = true
  elseif op.op === :set && op.field === :tap_pos
    # the live ratio on the nameplate grid: tap = neutral / (1 + n * step),
    # the same rule applyTapNameplate! documents
    br.has_ratio_tap || throw(ArgumentError(string("branch ", internal, " has no ratio tap changer")))
    _undo_set!(log, br, :tap_ratio, br.ratio / (1.0 + op.value * br.tap_step))
  elseif op.op === :set && op.field === :angle_deg
    _undo_set!(log, br, :phase_shift_deg, Float64(op.value))
  else
    throw(ArgumentError(string("op ", op.op, " field ", op.field, " does not apply to a branch")))
  end
  return nothing
end

function _scenario_prosumer!(log::UndoLog, net::Net, internal::Int, op::PatchOp)
  1 <= internal <= length(net.prosumpsVec) || throw(ArgumentError(string("prosumer index ", internal, " out of range")))
  ps = net.prosumpsVec[internal]
  if op.op === :status
    if op.value == 0.0
      # the contingency semantics: an outaged unit is REMOVED from the
      # copy, exactly like _remove_contingency_generator!, but with the
      # deletion recorded so restore! can reinsert it
      _undo_deleteat!(log, net.prosumpsVec, internal)
      log.refresh[] = true
    end
  elseif op.op === :set
    if op.field === :p
      _undo_set!(log, ps, :pVal, Float64(op.value))
    elseif op.field === :q
      _undo_set!(log, ps, :qVal, Float64(op.value))
    elseif op.field === :vm_pu
      _undo_set!(log, ps, :vm_pu, Float64(op.value))
      bus = _scf_prosumer_bus(ps)
      if 1 <= bus <= length(net.nodeVec)
        _undo_set!(log, net.nodeVec[bus], :_vm_pu, Float64(op.value))
      end
    else
      throw(ArgumentError(string("set field ", op.field, " does not apply to this appliance")))
    end
    log.refresh[] = true
  elseif op.op === :scale
    if op.field === nothing || op.field === :p
      _undo_set!(log, ps, :pVal, ps.pVal === nothing ? nothing : ps.pVal * op.factor)
    end
    if op.field === nothing || op.field === :q
      _undo_set!(log, ps, :qVal, ps.qVal === nothing ? nothing : ps.qVal * op.factor)
    end
    log.refresh[] = true
  end
  return nothing
end

function _scenario_shunt!(log::UndoLog, net::Net, internal::Int, op::PatchOp)
  1 <= internal <= length(net.shuntVec) || throw(ArgumentError(string("shunt index ", internal, " out of range")))
  sh = net.shuntVec[internal]
  if op.op === :status
    _undo_set!(log, sh, :status, Int(op.value))
  elseif op.op === :set && op.field === :b_pu
    y = sh.y_pu_shunt
    _undo_set!(log, sh, :y_pu_shunt, ComplexF64(real(y), Float64(op.value)))
    _undo_set!(log, sh, :B_shunt, Float64(op.value))
  else
    throw(ArgumentError(string("op ", op.op, " field ", op.field, " does not apply to a shunt")))
  end
  return nothing
end

function _scenario_link!(log::UndoLog, net::Net, internal::Int, op::PatchOp)
  1 <= internal <= length(net.linkVec) || throw(ArgumentError(string("link index ", internal, " out of range")))
  lk = net.linkVec[internal]
  op.op === :status || throw(ArgumentError("only status applies to a link"))
  _undo_set!(log, lk, :status, Int(op.value))
  log.refresh[] = true
  return nothing
end

"""
    apply!(net, ops, index) -> UndoLog

Apply the patch operations to the working copy `net`, addressing
components through the [`ScenarioIndex`](@ref), and record every change.
Derived state (node injection sums, bus types, Q-limit tables, the
isolated-node registry) is refreshed once after the patches, exactly as
the importers do after construction.
"""
function apply!(net::Net, ops::Vector{PatchOp}, index::ScenarioIndex)::UndoLog
  log = UndoLog()
  for op in ops
    internal = get(index.internal_by_id, op.id, nothing)
    internal === nothing && throw(ArgumentError(string("component id ", op.id, " has no internal index")))
    kind = index.kind_by_id[op.id]
    if kind === :branch || kind === :transformer
      _scenario_branch!(log, net, internal, op)
    elseif kind === :generator || kind === :load || kind === :external_grid
      _scenario_prosumer!(log, net, internal, op)
    elseif kind === :shunt
      _scenario_shunt!(log, net, internal, op)
    elseif kind === :link
      _scenario_link!(log, net, internal, op)
    else
      throw(ArgumentError(string("no apply rule for component class ", kind)))
    end
  end
  log.refresh[] && _scenario_refresh_derived!(net)
  return log
end

# the derived state every importer recomputes after construction; restore!
# and apply! run the same sequence, so base and restored copy derive
# identically (bitwise, because the primary fields are bitwise equal)
function _scenario_refresh_derived!(net::Net)
  # a node without generators keeps `nothing` (the constructors never
  # touch its sums); only generator-carrying nodes hold Float sums
  for node in net.nodeVec
    node._pƩGen = nothing
    node._qƩGen = nothing
  end
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    bus = _scf_prosumer_bus(ps)
    (1 <= bus <= length(net.nodeVec)) || continue
    node = net.nodeVec[bus]
    node._pƩGen = something(node._pƩGen, 0.0) + something(ps.pVal, 0.0)
    node._qƩGen = something(node._qƩGen, 0.0) + something(ps.qVal, 0.0)
  end
  refreshBusTypesFromProsumers!(net)
  _buildQLimits!(net)
  return nothing
end

"""
    restore!(net, undo) -> Nothing

Replay the undo log in reverse: reinsert every deleted element at its
recorded position and write back every recorded field value, then refresh
the derived state once. After this the working copy equals the base
bitwise on every field a patch can touch (the step-2 acceptance test
checks exactly that with the SCF reflection comparison).
"""
function restore!(net::Net, undo::UndoLog)
  for entry in Iterators.reverse(undo.entries)
    if entry.kind === :field
      setfield!(entry.obj, entry.key, entry.old)
    elseif entry.kind === :deleted
      insert!(entry.obj, entry.key, entry.old)
    end
  end
  undo.refresh[] && _scenario_refresh_derived!(net)
  return nothing
end

"""
    validate_scenarios(set, index, net) -> ScenarioSet

The net-aware validation of step 2 on top of the structural one: a
`tap_pos` patch is rejected at load time when the transformer carries an
ACTIVE tap controller (the message names it), and accepted only where a
ratio tap changer exists. Tap-step scenarios on unregulated transformers
are a first-class case (maintainer decision 1 of the scenario task).
"""
function validate_scenarios(set::ScenarioSet, index::ScenarioIndex, net::Net)::ScenarioSet
  validate_scenarios(set, index)
  for s in set.scenarios
    for (i, op) in enumerate(s.ops)
      (op.op === :set && op.field === :tap_pos) || continue
      internal = get(index.internal_by_id, op.id, nothing)
      internal === nothing && continue
      (1 <= internal <= length(net.branchVec)) || continue
      br = net.branchVec[internal]
      br.has_ratio_tap || _scenario_error(s.name, i, string("branch ", internal, " has no ratio tap changer"))
      ctrl = _scenario_branch_tap_controller(net, internal)
      ctrl === nothing || _scenario_error(s.name, i, string("transformer ", internal, " is regulated by controller ", _scenario_controller_name(ctrl), "; a tap_pos patch would fight it"))
    end
  end
  return set
end

# tap controllers live on the transformer WINDINGS; net.trafos runs
# parallel to the transformer-classified branches in branch order
function _scenario_branch_tap_controller(net::Net, internal::Int)
  trafo_k = 0
  for i in 1:internal
    _scf_is_transformer(net.branchVec[i]) && (trafo_k += 1)
  end
  (_scf_is_transformer(net.branchVec[internal]) && 1 <= trafo_k <= length(net.trafos)) || return nothing
  trafo = net.trafos[trafo_k]
  for winding in (trafo.side1, trafo.side2, trafo.side3)
    winding === nothing && continue
    winding.controls === nothing && continue
    isempty(winding.controls) || return first(winding.controls)
  end
  return nothing
end

_scenario_controller_name(ctrl) = hasproperty(ctrl, :trafo) && !isempty(String(getproperty(ctrl, :trafo))) ? String(getproperty(ctrl, :trafo)) : string(nameof(typeof(ctrl)), " targeting ", hasproperty(ctrl, :target_bus) ? String(getproperty(ctrl, :target_bus)) : "?")

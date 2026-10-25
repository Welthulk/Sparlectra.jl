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
# purpose: apply! executes a
#          scenario's patch operations on a working copy with an undo log
#          of every change, restore! replays the log in reverse and
#          recomputes the derived state, so the copy returns to the base
#          BITWISE on every field a patch can touch. Status patches use
#          the same mutation functions the contingency runner uses, so the
#          engine's full run reproduces today's results. The QSTS
#          vocabulary of a series index (CGMES series union, SCF time
#          series: park and unpark, bus isolation, pure field writes) runs
#          in three passes.

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

# A scenario that sets the tap position or the angle states the step, so
# the mark of an assumed middle position goes, as when a tap controller
# moves the changer (_clear_tap_step_assumed!). The two fields it writes
# go through the undo log first, so restore! brings the mark back.
function _scenario_clear_step_assumed!(log::UndoLog, br::Branch, which::Symbol)
  _undo_field!(log, br, which === :ratio ? :ratio_step_assumed : :phase_step_assumed)
  _undo_field!(log, br, :se_tap_changers)
  _clear_tap_step_assumed!(br, which)
  return nothing
end

function _scenario_branch!(log::UndoLog, net::Net, internal::Int, op::PatchOp)
  1 <= internal <= length(net.branchVec) || throw(ArgumentError(string("branch index ", internal, " out of range")))
  br = net.branchVec[internal]
  if op.op === :status
    _undo_set!(log, br, :status, Int(op.value))
    _undo_set!(log, br, :from_status, Int(op.value))
    _undo_set!(log, br, :to_status, Int(op.value))
    if op.value == 0.0
      # the bus shunt part a MATPOWER reimport records for this branch leaves
      # with it, as on the outage path (_remove_branch_shunt_parts!); the
      # undo log restores the admittance (D2 of #355)
      for sh in net.shuntVec
        part = get(sh.branch_parts, internal, nothing)
        part === nothing || _undo_set!(log, sh, :y_pu_shunt, sh.y_pu_shunt - part)
      end
    end
    log.refresh[] = true
  elseif op.op === :set && op.field === :tap_pos
    # the live ratio n positions from `ratio` on the grid of the branch, the
    # rule applyTapNameplate! documents: a CGMES changer (tap_ratio_grid
    # :linear) steps ratio + n * step or along its table points, nameplate
    # data steps neutral / (1 + n * step). The fraction rule on a CGMES
    # changer moved the ratio the other way (3W fixture, n = 8 from 1.0:
    # 0.926 instead of 1.08).
    br.has_ratio_tap || throw(ArgumentError(string("branch ", internal, " has no ratio tap changer")))
    if br.tap_ratio_grid === :linear && _tap_r1_model(br) === nothing
      (br.tap_ratio_points === nothing || isinteger(op.value)) || throw(ArgumentError(string("branch ", internal, " steps along a tap table; tap_pos needs a whole position, got ", op.value)))
      _undo_set!(log, br, :tap_ratio, _ratio_grid_shift(br, br.ratio, op.value; window = (-Inf, Inf)))
    else
      _undo_set!(log, br, :tap_ratio, br.ratio / (1.0 + op.value * br.tap_step))
    end
    _scenario_clear_step_assumed!(log, br, :ratio)
  elseif op.op === :set && op.field === :angle_deg
    _undo_set!(log, br, :phase_shift_deg, Float64(op.value))
    _scenario_clear_step_assumed!(log, br, :phase)
  else
    throw(ArgumentError(string("op ", op.op, " field ", op.field, " does not apply to a branch")))
  end
  return nothing
end

# keep the bus load sums (`_pƩLoad`, `_qƩLoad`) in step with a changed or
# removed load: _scenario_refresh_derived! rebuilds only the generator sums,
# and the island totals, the passive-bus test and the link-flow allocation
# read the load sums. The change is a difference on the stored sum (a bus
# can carry load added without a prosumer), written through the undo log so
# restore! returns the sum bitwise.
function _shift_bus_load!(log::UndoLog, net::Net, ps::ProSumer, dp::Float64, dq::Float64)
  isGenerator(ps) && return nothing
  bus = _scf_prosumer_bus(ps)
  (1 <= bus <= length(net.nodeVec)) || return nothing
  node = net.nodeVec[bus]
  dp == 0.0 || _undo_set!(log, node, :_pƩLoad, something(node._pƩLoad, 0.0) + dp)
  dq == 0.0 || _undo_set!(log, node, :_qƩLoad, something(node._qƩLoad, 0.0) + dq)
  return nothing
end

function _scenario_prosumer!(log::UndoLog, net::Net, internal::Int, op::PatchOp)
  1 <= internal <= length(net.prosumpsVec) || throw(ArgumentError(string("prosumer index ", internal, " out of range")))
  ps = net.prosumpsVec[internal]
  if op.op === :status
    # status 0 (the removal) is done by apply! after every other op; status 1
    # on an in-service unit changes nothing
    return nothing
  elseif op.op === :set
    if op.field === :p
      _shift_bus_load!(log, net, ps, Float64(op.value) - something(ps.pVal, 0.0), 0.0)
      _undo_set!(log, ps, :pVal, Float64(op.value))
    elseif op.field === :q
      _shift_bus_load!(log, net, ps, 0.0, Float64(op.value) - something(ps.qVal, 0.0))
      _undo_set!(log, ps, :qVal, Float64(op.value))
    elseif op.field === :vm_pu
      _undo_set!(log, ps, :vm_pu, Float64(op.value))
      bus = _scf_prosumer_bus(ps)
      if 1 <= bus <= length(net.nodeVec)
        _undo_set!(log, net.nodeVec[bus], :_vm_pu, Float64(op.value))
      end
    elseif op.field === :angle_deg
      # the angle of an external grid (SCF source.u_ref_angle), written where
      # the reader puts it: on the unit and, through addProsumer!, on its bus.
      # initialVrect holds a slack bus at that angle, flat start or not, so on
      # the reference unit every angle of its island moves with the value; on
      # a source that is not the reference it is the start angle of its bus.
      _undo_set!(log, ps, :va_deg, Float64(op.value))
      bus = _scf_prosumer_bus(ps)
      if 1 <= bus <= length(net.nodeVec)
        _undo_set!(log, net.nodeVec[bus], :_va_deg, Float64(op.value))
      end
    else
      throw(ArgumentError(string("set field ", op.field, " does not apply to this appliance")))
    end
    log.refresh[] = true
  elseif op.op === :scale
    if op.field === nothing || op.field === :p
      ps.pVal === nothing || _shift_bus_load!(log, net, ps, ps.pVal * op.factor - ps.pVal, 0.0)
      _undo_set!(log, ps, :pVal, ps.pVal === nothing ? nothing : ps.pVal * op.factor)
    end
    if op.field === nothing || op.field === :q
      ps.qVal === nothing || _shift_bus_load!(log, net, ps, 0.0, ps.qVal * op.factor - ps.qVal)
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

## --- the QSTS vocabulary (series index, origin :cgmes_series or :scf_series) --

# the parked state of a generator, load or external grid, written by
# `:set :in_service 0`; `:set :in_service 1` copies the same fields back
# from index.parked_values. The type decides isGenerator, the controller
# fields isRegulating and referencePri isSlack, so a parked unit is a
# 0 MW load in every reader. The unit keeps its position in prosumpsVec,
# which MachineVoltageControl.prosumer_idx and the HVDC records rely on.
const _PARKED_STATE = (
  :proSumptionType => Consumption,
  :pVal => 0.0,
  :qVal => 0.0,
  :isRegulated => false,
  :referencePri => nothing,
  :referencePriority => 0,
  :participationFactor => nothing,
  :minQ => nothing,
  :maxQ => nothing,
  :vset_adjust => nothing,
  :vstep_pu => nothing,
  :tap_steps_down => nothing,
  :tap_steps_up => nothing,
  :quController => nothing,
  :puController => nothing,
)

# the share a prosumer adds to its bus load sum (`_pƩLoad`/`_qƩLoad`); a
# generator adds nothing there, its sums are rebuilt by the refresh
_scenario_load_part(ps::ProSumer)::Tuple{Float64,Float64} = isGenerator(ps) ? (0.0, 0.0) : (something(ps.pVal, 0.0), something(ps.qVal, 0.0))

# move the bus load sum by what a park or unpark changed; called after the
# writes, so whenever the difference is not zero the unit is a load and
# _shift_bus_load! takes it
function _scenario_follow_load_part!(log::UndoLog, net::Net, ps::ProSumer, before::Tuple{Float64,Float64})
  after = _scenario_load_part(ps)
  _shift_bus_load!(log, net, ps, after[1] - before[1], after[2] - before[2])
  return nothing
end

function _park_prosumer!(log::UndoLog, net::Net, ps::ProSumer)
  before = _scenario_load_part(ps)
  for (field, value) in _PARKED_STATE
    _undo_set!(log, ps, field, value)
  end
  _scenario_follow_load_part!(log, net, ps, before)
  return nothing
end

function _unpark_prosumer!(log::UndoLog, net::Net, ps::ProSumer, parked::ProSumer)
  before = _scenario_load_part(ps)
  # the controller and voltage-adjust objects are immutable, so sharing
  # them with the index between worker copies is safe
  for (field, _) in _PARKED_STATE
    _undo_set!(log, ps, field, getfield(parked, field))
  end
  _scenario_follow_load_part!(log, net, ps, before)
  return nothing
end

function _qsts_prosumer!(log::UndoLog, net::Net, internal::Int, op::PatchOp, index::ScenarioIndex)
  1 <= internal <= length(net.prosumpsVec) || throw(ArgumentError(string("prosumer index ", internal, " out of range")))
  ps = net.prosumpsVec[internal]
  f = op.field
  if f === :in_service
    if op.value == 0.0
      _park_prosumer!(log, net, ps)
    else
      parked = get(index.parked_values, op.id, nothing)
      parked === nothing && throw(ArgumentError(string("set in_service 1 on id ", op.id, ": the index holds no parked values for it")))
      _unpark_prosumer!(log, net, ps, parked)
    end
  elseif f === :qmin
    _undo_set!(log, ps, :minQ, op.value)
  elseif f === :qmax
    _undo_set!(log, ps, :maxQ, op.value)
  elseif f === :regulated
    _undo_set!(log, ps, :isRegulated, op.value == 1.0)
  elseif f === :reference
    # referencePri holds the bus of the unit that IS the reference, as the
    # importers write it (cgmes_mapping, createnet_powermat)
    _undo_set!(log, ps, :referencePri, op.value == 1.0 ? getPosumerBusIndex(ps) : nothing)
  elseif f === :reference_priority
    _undo_set!(log, ps, :referencePriority, Int(op.value))
  elseif f === :participation
    _undo_set!(log, ps, :participationFactor, op.value)
  else
    throw(ArgumentError(string("set ", f, " does not apply to an appliance")))
  end
  log.refresh[] = true
  return nothing
end

function _qsts_branch!(log::UndoLog, net::Net, internal::Int, op::PatchOp, index::ScenarioIndex)
  1 <= internal <= length(net.branchVec) || throw(ArgumentError(string("branch index ", internal, " out of range")))
  br = net.branchVec[internal]
  f = op.field
  if f === :status || f === :from_status || f === :to_status
    # a pure field write: the delta carries all three flags of the snapshot,
    # which keep the status invariant there already. The SCF :status op
    # writes all three and also removes MATPOWER bus-shunt parts, which a
    # union delta must not do.
    _undo_set!(log, br, f, Int(op.value))
  elseif f === :ratio
    # the importer passes one ratio to _addPIModelTrafo_by_idx!, and the
    # Branch constructor stores it as both the ratio and the live tap_ratio
    br.ratio != 0.0 || throw(ArgumentError(string("branch ", internal, " is a line; ratio applies to a transformer")))
    _undo_set!(log, br, :ratio, Float64(op.value))
    _undo_set!(log, br, :tap_ratio, Float64(op.value))
  elseif f === :shift_deg
    # the same for the shift: the importer passes one shift_deg, stored as
    # angle and phase_shift_deg. Both move, because the solver reads
    # phase_shift_deg and the branch flows (calcNetLosses!) read angle.
    br.ratio != 0.0 || throw(ArgumentError(string("branch ", internal, " is a line; shift_deg applies to a transformer")))
    _undo_set!(log, br, :angle, Float64(op.value))
    _undo_set!(log, br, :phase_shift_deg, Float64(op.value))
  elseif f === :tap_data
    # the tap changer data of a snapshot as one record (the ratio points and
    # the step impedance table are no single values); the record fields
    # carry the names of the branch fields they replace
    br.ratio != 0.0 || throw(ArgumentError(string("branch ", internal, " is a line; tap_data applies to a transformer")))
    recs = get(index.tap_data, op.id, BranchTapData[])
    k = Int(op.value)
    1 <= k <= length(recs) || throw(ArgumentError(string("set tap_data ", k, " on id ", op.id, ": the index holds ", length(recs), " tap data record(s) for it")))
    rec = recs[k]
    for g in fieldnames(BranchTapData)
      _undo_set!(log, br, g, getfield(rec, g))
    end
  elseif f === :r_pu || f === :x_pu
    # A snapshot changes the equipment, so the base field (r_base_pu,
    # x_base_pu: short circuit and the exports) moves with the live one (the
    # power flow). The live value differs from the base only while a
    # series-FACTS control run holds its operating point there or a typed
    # tap model derives it from the base. The template is an import state
    # without a control run, and _qsts_check_impedance has rejected a
    # taps_derived branch before any write.
    _undo_set!(log, br, f, Float64(op.value))
    _undo_set!(log, br, f === :r_pu ? :r_base_pu : :x_base_pu, Float64(op.value))
  elseif f === :g_from_pu || f === :b_from_pu || f === :g_to_pu || f === :b_to_pu
    # set_branch_shunt! is the one writer of the terminal values and keeps
    # the totals g_pu, b_pu their sums; the undo log takes all six fields
    for g in (:g_from_pu, :b_from_pu, :g_to_pu, :b_to_pu, :g_pu, :b_pu)
      _undo_field!(log, br, g)
    end
    v = Float64(op.value)
    set_branch_shunt!(br; g_from_pu = f === :g_from_pu ? v : br.g_from_pu, b_from_pu = f === :b_from_pu ? v : br.b_from_pu,
                      g_to_pu = f === :g_to_pu ? v : br.g_to_pu, b_to_pu = f === :b_to_pu ? v : br.b_to_pu)
  else
    throw(ArgumentError(string("set ", f, " does not apply to a branch")))
  end
  log.refresh[] = true
  return nothing
end

function _qsts_shunt!(log::UndoLog, net::Net, internal::Int, op::PatchOp)
  1 <= internal <= length(net.shuntVec) || throw(ArgumentError(string("shunt index ", internal, " out of range")))
  sh = net.shuntVec[internal]
  op.field === :g_pu || throw(ArgumentError(string("set ", op.field, " does not apply to a shunt")))
  # the real part only, G_shunt mirrors it as B_shunt mirrors the b_pu op
  _undo_set!(log, sh, :y_pu_shunt, ComplexF64(Float64(op.value), imag(sh.y_pu_shunt)))
  _undo_set!(log, sh, :G_shunt, Float64(op.value))
  return nothing
end

function _qsts_tap_controller!(log::UndoLog, ctrls::Vector{PowerTransformerControl}, internal::Int, op::PatchOp)
  1 <= internal <= length(ctrls) || throw(ArgumentError(string("tap controller index ", internal, " out of range")))
  c = ctrls[internal]
  f = op.field
  if f === :enabled
    _undo_set!(log, c, :enabled, op.value == 1.0)
  elseif f === :target_vm_pu
    _undo_set!(log, c, :target_vm_pu, op.value)
  elseif f === :deadband_vm_pu
    _undo_set!(log, c, :deadband_vm_pu, Float64(op.value))
  elseif f === :p_target_mw
    _undo_set!(log, c, :p_target_mw, op.value)
  elseif f === :deadband_p_mw
    _undo_set!(log, c, :deadband_p_mw, Float64(op.value))
  else
    throw(ArgumentError(string("set ", f, " does not apply to a tap controller")))
  end
  return nothing
end

# A bus op states the island status of the snapshot explicitly, because the
# derived state cannot: refreshBusTypesFromProsumers! skips an Isolated bus,
# and markIsolatedBuses! types back only the buses it isolated itself
# (net._marked_iso), so `Isolated` from an op sticks once set; an op takes
# the bus out of that list for the same reason. `iso` is the working set of
# isolated buses; apply! writes it as a NEW isoNodes vector through the undo
# log once, so restore! puts the old vector object back even if a later
# markIsolatedBuses! has refilled the new one in place.
function _qsts_bus!(log::UndoLog, net::Net, bus::Int, op::PatchOp, iso::Union{Nothing,Set{Int}})::Union{Nothing,Set{Int}}
  1 <= bus <= length(net.nodeVec) || throw(ArgumentError(string("bus ", bus, " out of range")))
  node = net.nodeVec[bus]
  f = op.field
  if f === :isolated
    iso === nothing && (iso = Set{Int}(net.isoNodes))
    # the op owns the type from here: a later marking does not type it back
    bus in net._marked_iso && _undo_set!(log, net, :_marked_iso, filter(!=(bus), net._marked_iso))
    if op.value == 1.0
      # what markIsolatedBuses! writes for an unconnected bus
      _undo_set!(log, node, :_nodeType, Isolated)
      _undo_set!(log, node, :_vm_pu, 0.0)
      _undo_set!(log, node, :_va_deg, 0.0)
      push!(iso, bus)
    else
      # PQ only takes the bus out of the skip in the refresh, which then
      # types it from its prosumers (Slack, PV or PQ)
      _undo_set!(log, node, :_nodeType, PQ)
      delete!(iso, bus)
    end
  elseif f === :vm_start
    _undo_set!(log, node, :_vm_pu, Float64(op.value))
  elseif f === :va_start
    _undo_set!(log, node, :_va_deg, Float64(op.value))
  else
    throw(ArgumentError(string("set ", f, " does not apply to a bus")))
  end
  log.refresh[] = true
  return iso
end

# position of an op in the canonical order of a QSTS delta: pass 1 parks
# and unparks, pass 2 runs every other op, bus start values after the
# appliance ops (a generator vm_pu also writes the bus voltage, the
# snapshot's start value wins) and the isolation last (an isolated bus keeps
# vm and va 0). Inside a pass the ops run by id, op and field, so two writes
# to one bus load sum add up in the same order whatever the delta order.
function _qsts_op_key(op::PatchOp, kind::Symbol)::Tuple{Int,Int,Symbol,Symbol}
  pass = if op.op === :set && op.field === :in_service
    1
  elseif kind === :bus && op.field === :isolated
    4
  elseif kind === :bus
    3
  else
    2
  end
  return (pass, op.id, op.op, something(op.field, :none))
end

# The series impedance checks of a QSTS delta, on the unpatched net before
# any write, so a rejected delta leaves the net as it was: r_pu and x_pu of a
# branch must not both end at 0 (inv(0 + 0im) is NaN in the admittance
# stamp; x_pu = 0 alone is fine for the AC solver, the DC model checks it
# itself), and a branch whose live impedance a typed tap model derives from
# the base (taps_derived) takes no plain r/x write.
function _qsts_check_impedance(net::Net, resolved::Vector{Tuple{PatchOp,Int,Symbol}})
  final = Dict{Int,Tuple{Float64,Float64}}()
  for (op, internal, kind) in resolved
    (op.op === :set && (op.field === :r_pu || op.field === :x_pu) && (kind === :branch || kind === :transformer)) || continue
    1 <= internal <= length(net.branchVec) || throw(ArgumentError(string("branch index ", internal, " out of range")))
    br = net.branchVec[internal]
    br.taps_derived && throw(ArgumentError(string("set ", op.field, " on branch ", internal, ": a typed tap model derives its live impedance from the base, a plain write would break that")))
    r, x = get(final, internal, (br.r_pu, br.x_pu))
    final[internal] = op.field === :r_pu ? (Float64(op.value), x) : (r, Float64(op.value))
  end
  for internal in sort!(collect(keys(final)))
    r, x = final[internal]
    (r == 0.0 && x == 0.0) && throw(ArgumentError(string("set r_pu/x_pu on branch ", internal, ": r_pu and x_pu would both be 0, the series admittance 1/(r + jx) is undefined")))
  end
  return nothing
end

function _apply_qsts!(net::Net, ops::Vector{PatchOp}, index::ScenarioIndex)::UndoLog
  log = UndoLog()
  # every internal index is resolved on the unpatched net, before any op
  resolved = Vector{Tuple{PatchOp,Int,Symbol}}(undef, length(ops))
  for (k, op) in enumerate(ops)
    internal = get(index.internal_by_id, op.id, nothing)
    internal === nothing && throw(ArgumentError(string("component id ", op.id, " has no internal index")))
    resolved[k] = (op, internal, index.kind_by_id[op.id])
  end
  _qsts_check_impedance(net, resolved)
  # stable sort: only a repeated (id, op, field), which validate_scenarios
  # rejects, keeps its delta order
  order = sortperm(resolved; by = r -> _qsts_op_key(r[1], r[3]))
  deletions = Int[]
  iso = nothing
  ctrls = nothing
  for k in order
    op, internal, kind = resolved[k]
    if _scenario_qsts_rule(op) === nothing
      _scenario_apply_op!(log, net, kind, internal, op, deletions)
    elseif kind === :bus
      iso = _qsts_bus!(log, net, internal, op, iso)
    elseif kind === :tap_controller
      ctrls === nothing && (ctrls = _tap_controllers(net))
      _qsts_tap_controller!(log, ctrls, internal, op)
    elseif kind === :branch || kind === :transformer
      _qsts_branch!(log, net, internal, op, index)
    elseif kind === :shunt
      _qsts_shunt!(log, net, internal, op)
    elseif kind === :generator || kind === :load || kind === :external_grid
      _qsts_prosumer!(log, net, internal, op, index)
    else
      throw(ArgumentError(string("no apply rule for set ", op.field, " on component class ", kind)))
    end
  end
  iso === nothing || _undo_set!(log, net, :isoNodes, sort!(collect(iso)))
  _scenario_delete_prosumers!(log, net, deletions)
  # pass 3: one refresh of the derived state, recorded so restore! runs it too
  isempty(ops) || (log.refresh[] = true)
  log.refresh[] && _scenario_refresh_derived!(net)
  return log
end

## --- apply! -------------------------------------------------------------------

# one op of the SCF vocabulary; a prosumer outage is only collected in
# `deletions`, _scenario_delete_prosumers! runs it after every other op
function _scenario_apply_op!(log::UndoLog, net::Net, kind::Symbol, internal::Int, op::PatchOp, deletions::Vector{Int})
  if kind === :branch || kind === :transformer
    _scenario_branch!(log, net, internal, op)
  elseif kind === :generator || kind === :load || kind === :external_grid
    if op.op === :status && op.value == 0.0
      1 <= internal <= length(net.prosumpsVec) || throw(ArgumentError(string("prosumer index ", internal, " out of range")))
      push!(deletions, internal)
    else
      _scenario_prosumer!(log, net, internal, op)
    end
  elseif kind === :shunt
    _scenario_shunt!(log, net, internal, op)
  elseif kind === :link
    _scenario_link!(log, net, internal, op)
  else
    throw(ArgumentError(string("no apply rule for component class ", kind)))
  end
  return nothing
end

# prosumer outages delete from prosumpsVec and shift every later position,
# so they run after all other ops, highest position first (two generator
# outages in ascending order removed the wrong second unit, D1 of #355)
function _scenario_delete_prosumers!(log::UndoLog, net::Net, deletions::Vector{Int})
  for internal in sort!(unique(deletions); rev = true)
    # the contingency semantics: an outaged unit is REMOVED from the copy,
    # like _remove_contingency_generator!, with the deletion recorded so
    # restore! reinserts it; a removed load leaves the bus load sum
    ps = net.prosumpsVec[internal]
    _shift_bus_load!(log, net, ps, -something(ps.pVal, 0.0), -something(ps.qVal, 0.0))
    _undo_deleteat!(log, net.prosumpsVec, internal)
    log.refresh[] = true
  end
  return nothing
end

"""
    apply!(net, ops, index) -> UndoLog

Apply the patch operations to the working copy `net`, addressing
components through the [`ScenarioIndex`](@ref), and record every change.
Derived state (node injection sums, bus types, Q-limit tables, the
isolated-node registry) is refreshed once after the patches, exactly as
the importers do after construction. Prosumer outages (status 0 on a
generator, load or external grid) are applied after every other
operation, highest position first, so the order of the ops does not
change the result; a branch outage also removes the bus shunt part a
MATPOWER reimport records for the branch.

With a series index (`index.origin` `:cgmes_series` or `:scf_series`) the
ops run in three passes: every `set in_service` (park and unpark; a parked
unit keeps its position in `prosumpsVec`), then every other op, then one
refresh. Inside a pass the
order is fixed by id and field, bus start values follow the appliance ops
and bus isolation comes last, so the order of the ops in a delta does not
change the result. Bus isolation ops replace `net.isoNodes` by a new
vector. Every write goes through the undo log, so [`restore!`](@ref)
returns the template bitwise. Before any write, a delta is rejected with
an `ArgumentError` when it leaves `r_pu` and `x_pu` of a branch both 0 or
writes `r_pu`/`x_pu` on a branch whose impedance a typed tap model derives
(`taps_derived`); the net is then unchanged.
"""
function apply!(net::Net, ops::Vector{PatchOp}, index::ScenarioIndex)::UndoLog
  index.origin in _QSTS_ORIGINS && return _apply_qsts!(net, ops, index)
  log = UndoLog()
  # every internal index is resolved on the unpatched net; prosumer outages
  # run after all other ops (_scenario_delete_prosumers!)
  deletions = Int[]
  for op in ops
    internal = get(index.internal_by_id, op.id, nothing)
    internal === nothing && throw(ArgumentError(string("component id ", op.id, " has no internal index")))
    _scenario_apply_op!(log, net, index.kind_by_id[op.id], internal, op, deletions)
  end
  _scenario_delete_prosumers!(log, net, deletions)
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
are a first-class case.
"""
validate_scenarios(set::ScenarioSet, index::ScenarioIndex, net::Net)::ScenarioSet = _validate_scenarios_net(set, index, net, true)

# The net-aware checks behind validate_scenarios(set, index, net). With
# `regulated = false` a tap_pos op on a transformer with a tap controller
# passes: a run whose controllers do not run (QSTS mode A without
# qsts.snapshot_controls) takes it as the recorded position of the step,
# and nothing moves the tap away from it.
function _validate_scenarios_net(set::ScenarioSet, index::ScenarioIndex, net::Net, regulated::Bool)::ScenarioSet
  validate_scenarios(set, index)
  for s in set.scenarios
    for (i, op) in enumerate(s.ops)
      (op.op === :set && op.field === :tap_pos) || continue
      internal = get(index.internal_by_id, op.id, nothing)
      internal === nothing && continue
      (1 <= internal <= length(net.branchVec)) || continue
      br = net.branchVec[internal]
      br.has_ratio_tap || _scenario_error(s.name, i, string("branch ", internal, " has no ratio tap changer"))
      regulated || continue
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

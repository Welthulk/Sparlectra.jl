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

# file: src/qsts/qsts_sequential.jl
# purpose: QSTS mode B (qsts.mode = sequential): the steps of a series in
#          order, each on a fresh copy of the template with the full delta
#          of the step, the outer-loop controllers running, and the state of
#          the tap and bank controllers carried from the last solved step
#          (absolute write of taps and banks, warm start, dead time,
#          hunting memory). Resync steps put taps and banks back on the
#          values of the series. Several series run in parallel, the steps
#          of one series in order. The result is the QstsResult of mode A
#          with controller elements, the switching list and the metrics.

# the loosely typed keywords of a sequential run behind one concrete struct
# (compile rule 2)
struct _QstsSeqOptions
  parallel_enabled::Union{Nothing,Bool}             # the parallel map over the series
  parallel_max_tasks::Union{Nothing,Int}
  parallel_min_work_items::Union{Nothing,Int}
  progress::Any                                     # nothing, or the _BatchProgress of the run, one tick per step
  island_parallel::Bool                             # one series: the islands of a step after runtime.parallel; several: off
  print_island_table::Bool                          # the base case prints the island table (one series only)
end

# the position of one transformer branch of a tap controller (master or
# follower), as the carry writes it back
struct _QstsTapPos
  branch::Int                                       # position in T0.branchVec
  tap_ratio::Float64
  shift_deg::Float64
  ratio_step::Int                                   # model step; typemin(Int) without a ratio model
  phase_step::Int                                   # likewise for the phase model
end

# What one step hands to the next. The list of what carries over is this
# struct; the restart test (a run from step k with the carry of step k - 1
# equals the whole run) proves it complete. Everything else of a step comes
# from T0 again: controller flags and secant memory, machine, series
# reactance, UPFC and HVDC controllers.
struct _QstsCarry
  step::Int                                         # last step run, 0 before step 1
  taps::Vector{_QstsTapPos}                         # per tap controller its master and followers, controller order
  banks::Vector{Float64}                            # bs_mvar per controller (NaN for a tap controller)
  positions::Vector{Int}                            # position per controller (typemin(Int): none, a continuous SVC)
  vm::Vector{Union{Nothing,Float64}}                # warm start: bus voltages of the last solved step
  va::Vector{Union{Nothing,Float64}}
  prev_vm::Vector{Float64}                          # vm_pu values of step `step` when it converged, else NaN
  prev_p::Vector{Float64}                           # p_from_mw likewise
  switch_count::Vector{Int}                         # per controller
  n_hunting::Vector{Int}
  lock_steps::Vector{Int}                           # steps the dead time still locks
  last_moves::Vector{NTuple{2,Tuple{Int,Int}}}      # (step, sign) of the last two counted moves, (0, 0) empty
  last_converged::Int                               # 0 before the first solved step
end

# everything a step reads, built once per series and shared read-only
struct _QstsSeqContext
  input::_QstsInput
  q::QstsConfig
  cold::_QstsContext                                # step configuration before the first warm start
  warm::_QstsContext                                # step configuration with a warm start (flatstart off)
  base_solved::Bool
  base_warm::Tuple{Vector{Union{Nothing,Float64}},Vector{Union{Nothing,Float64}}}   # base solution, else the start values of T0
  snapshot_starts::Vector{Tuple{Vector{Union{Nothing,Float64}},Vector{Union{Nothing,Float64}}}}   # sequential_start = snapshot, else empty
  ctrl_branches::Vector{Vector{Int}}                # per controller: master and followers (empty for a bank)
  ctrl_ids::Vector{String}
  bus_reported::BitVector                           # false for star points, which the band counts leave out
  controls::Bool                                    # control.enabled
end

# one step of a sequential run: the capture of mode A plus the controller
# values and the switching rows of the step
struct _QstsSeqOut
  out::_QstsStepOut
  ctrl::Matrix{Float64}
  switching::Vector{_QstsSwitchRow}
  resync::Bool
end

# one series after the loop, before the assembly into a QstsResult
struct _QstsSeriesRun
  input::_QstsInput
  config::SparlectraConfig
  elements::QstsElements
  base::QstsStepInfo
  outs::Vector{_QstsSeqOut}
  carry::_QstsCarry
  notes::Vector{String}
  warnings::Vector{String}                          # plausibility warnings of the base case, logged by the orchestrator
end

# The controllers whose state a sequential run carries: tap controllers
# (master and followers) and shunt bank controllers (MSC and SVC), in the
# order of collect_outer_controllers. Every copy of T0 lists them in the
# same order.
_qsts_carried_controllers(net::Net)::Vector{AbstractOuterController} = AbstractOuterController[c for c in collect_outer_controllers(net) if c isa PowerTransformerControl || c isa ShuntVoltageControl]

# the actuators a tap controller moves (control_propose_update!), the only
# ones the carry writes: an OLTC leaves the shift of the series alone, a
# PST the ratio
_qsts_moves_ratio(c::PowerTransformerControl)::Bool = c.control_ratio && c.mode in (:voltage, :voltage_and_branch_active_power)
_qsts_moves_phase(c::PowerTransformerControl)::Bool = c.control_phase && c.mode in (:branch_active_power, :voltage_and_branch_active_power)

# the model step whose resolved value lies nearest to `value` (the current
# step wins a tie); typemin(Int) without a model. After a :ratio op the
# model step is stale (the op writes the ratio only), so a position read
# after apply! looks the step up from the value.
function _qsts_nearest_model_step(br::Branch, which::Symbol, value::Float64)::Int
  m = _derived_tap_model(br, which)
  m === nothing && return typemin(Int)
  best = m.step
  bestd = abs(_model_step_value(br, which, m.step) - value)
  for s in m.lowStep:m.highStep
    d = abs(_model_step_value(br, which, s) - value)
    if d < bestd - 1e-12
      best = s
      bestd = d
    end
  end
  return best
end

# the bank susceptance on the net: the controller value while the shunt
# carries its stamp, else the value the series wrote to the shunt (a bank
# whose shunt the series addresses)
function _qsts_bank_value(net::Net, ctrl::ShuntVoltageControl)::Float64
  sh = _find_controlled_shunt(net, ctrl)
  sh.y_pu_shunt == ComplexF64(0.0, ctrl.bs_mvar) / net.baseMVA && return ctrl.bs_mvar
  return imag(sh.y_pu_shunt) * net.baseMVA
end

# taps per branch of every tap controller and the bank value per controller
function _qsts_read_state(net::Net, ctrls::Vector{AbstractOuterController}, ctrl_branches::Vector{Vector{Int}})::Tuple{Vector{_QstsTapPos},Vector{Float64}}
  taps = _QstsTapPos[]
  banks = fill(NaN, length(ctrls))
  for (c, ctrl) in enumerate(ctrls)
    if ctrl isa ShuntVoltageControl
      banks[c] = _qsts_bank_value(net, ctrl)
      continue
    end
    for j in ctrl_branches[c]
      br = net.branchVec[j]
      push!(taps, _QstsTapPos(j, br.tap_ratio, br.phase_shift_deg, _qsts_nearest_model_step(br, :ratio, br.tap_ratio), _qsts_nearest_model_step(br, :phase, br.phase_shift_deg)))
    end
  end
  return taps, banks
end

# the position of a tap controller from its master's tap: the model step,
# else the tap position counted from the ratio 1.0 (shift 0 for a PST), as
# buildTapControllerReportRows counts it
function _qsts_tap_position(ctrl::PowerTransformerControl, br::Branch, p::_QstsTapPos)::Int
  if _qsts_moves_ratio(ctrl)
    p.ratio_step == typemin(Int) || return p.ratio_step
    v = _tap_position(p.tap_ratio, 1.0, br.tap_step)
  else
    p.phase_step == typemin(Int) || return p.phase_step
    v = _tap_position(p.shift_deg, 0.0, br.phase_step_deg)
  end
  ismissing(v) && throw(ArgumentError(string("runQsts: transformer ", br.comp.cName, " of tap controller ", control_name(ctrl), " has no tap step; a sequential run counts tap positions in steps")))
  return v
end

# position per controller (typemin(Int) for a continuous SVC, which has no steps)
function _qsts_positions(net::Net, ctrls::Vector{AbstractOuterController}, ctrl_branches::Vector{Vector{Int}}, taps::Vector{_QstsTapPos}, banks::Vector{Float64})::Vector{Int}
  out = fill(typemin(Int), length(ctrls))
  t = 1
  for (c, ctrl) in enumerate(ctrls)
    if ctrl isa ShuntVoltageControl
      ctrl.step_mvar === nothing || (out[c] = round(Int, banks[c] / ctrl.step_mvar))
      continue
    end
    out[c] = _qsts_tap_position(ctrl, net.branchVec[ctrl_branches[c][1]], taps[t])
    t += length(ctrl_branches[c])
  end
  return out
end

# Steps between the reference (the tap the series gives the step) and the
# simulated tap of the master: the model step difference, else the steps of
# _ratio_steps_moved on the ratio grid (table points or tap_step), for a
# PST whole phase steps.
function _qsts_dpos_ref(ctrl::PowerTransformerControl, br::Branch, ref::_QstsTapPos, sim::_QstsTapPos)::Float64
  if _qsts_moves_ratio(ctrl)
    ref.ratio_step == typemin(Int) || return Float64(sim.ratio_step - ref.ratio_step)
    return Float64(_ratio_steps_moved(br, ref.tap_ratio, sim.tap_ratio))
  end
  ref.phase_step == typemin(Int) || return Float64(sim.phase_step - ref.phase_step)
  return br.phase_step_deg > 0.0 ? Float64(round(Int, (sim.shift_deg - ref.shift_deg) / br.phase_step_deg)) : Float64(sign(sim.shift_deg - ref.shift_deg))
end

# Whether a controller stands at an end of its control range: the ratio at
# tap_min/tap_max of the range (on a modelled or tabular changer also the
# last step inside it), the shift at its limits, a bank on its outermost
# block (a continuous SVC at bs_min/bs_max). The checks of
# control_apply_update!, restricted to the actuators the controller moves.
function _qsts_at_range_end(net::Net, ctrl::PowerTransformerControl, js::Vector{Int})::Bool
  br = net.branchVec[js[1]]
  if _qsts_moves_ratio(ctrl)
    lo, hi = _tap_control_range(ctrl, br)
    (isapprox(br.tap_ratio, lo; atol = 1e-12) || isapprox(br.tap_ratio, hi; atol = 1e-12) || _model_step_at_window_edge(ctrl, br) || _table_point_at_window_edge(ctrl, br)) && return true
  end
  _qsts_moves_phase(ctrl) && (isapprox(br.phase_shift_deg, br.phase_min_deg; atol = 1e-12) || isapprox(br.phase_shift_deg, br.phase_max_deg; atol = 1e-12)) && return true
  return false
end
function _qsts_at_range_end(::Net, ctrl::ShuntVoltageControl, ::Vector{Int})::Bool
  lo, hi = ctrl.step_mvar === nothing ? (ctrl.bs_min_mvar, ctrl.bs_max_mvar) : (_shunt_step_ceil(ctrl.bs_min_mvar, ctrl.step_mvar), _shunt_step_floor(ctrl.bs_max_mvar, ctrl.step_mvar))
  return isapprox(ctrl.bs_mvar, lo; atol = 1e-9) || isapprox(ctrl.bs_mvar, hi; atol = 1e-9)
end
_qsts_at_range_end(::Net, ::AbstractOuterController, ::Vector{Int})::Bool = false

# The carry onto a fresh copy, absolute: model steps through the resolver,
# else ratio and shift (with X(alpha)), then the impedance of a tap table
# for the step; followers take their own carried position; a bank its
# carried susceptance. control_apply_update! is not the path: its model
# step only moves to a neighbour, followers move relative to their own tap,
# and a bank call sets at_limit on the outermost block.
function _qsts_write_carry!(net::Net, ctrls::Vector{AbstractOuterController}, ctrl_branches::Vector{Vector{Int}}, carry::_QstsCarry)
  t = 0
  for (c, ctrl) in enumerate(ctrls)
    if ctrl isa ShuntVoltageControl
      _set_shunt_susceptance!(net, ctrl, carry.banks[c])
      continue
    end
    ctrl isa PowerTransformerControl || continue
    ratio = _qsts_moves_ratio(ctrl)
    for (i, j) in enumerate(ctrl_branches[c])
      t += 1
      p = carry.taps[t]
      br = net.branchVec[j]
      # a follower mirrors the ratio of its master only
      phase = _qsts_moves_phase(ctrl) && i == 1
      (ratio || phase) || continue
      _set_model_tap_steps!(br, ratio ? p.ratio_step : typemin(Int), phase ? p.phase_step : typemin(Int))
      ratio && _derived_tap_model(br, :ratio) === nothing && _set_ratio_tap!(br, p.tap_ratio)
      phase && _derived_tap_model(br, :phase) === nothing && _set_phase_shift!(net, ctrl, br, p.shift_deg)
      _apply_tap_step_impedance!(br)
    end
  end
  return nothing
end

# Resync: taps and banks stay at the values of the series. A modelled
# transformer gets the model step that matches its ratio and shift, since a
# :ratio op leaves the model step stale and the first move of the
# controller would start from it. A bank takes the value its shunt carries.
function _qsts_sync_model_steps!(net::Net, ctrls::Vector{AbstractOuterController}, ctrl_branches::Vector{Vector{Int}})
  for (c, ctrl) in enumerate(ctrls)
    if ctrl isa ShuntVoltageControl
      v = _qsts_bank_value(net, ctrl)
      v == ctrl.bs_mvar || _set_shunt_susceptance!(net, ctrl, v)
      continue
    end
    for j in ctrl_branches[c]
      br = net.branchVec[j]
      for which in (:ratio, :phase)
        m = _derived_tap_model(br, which)
        m === nothing && continue
        m.step = _qsts_nearest_model_step(br, which, which === :ratio ? br.tap_ratio : br.phase_shift_deg)
      end
    end
  end
  return nothing
end

# Parameter drift of a controlled transformer with a tap table: the delta
# writes r, x, g, b of the step at the reference tap, while the table
# rescales its own base values when the tap moves. The base values of the
# copy are put on the step: the values after apply! divided by the factors
# of the reference tap. A value equal to base times factor (the delta did
# not touch it) keeps its base bit for bit.
function _qsts_rebase_tap_table!(net::Net, ctrl_branches::Vector{Vector{Int}})
  for js in ctrl_branches, j in js
    br = net.branchVec[j]
    t = br.tap_step_impedance
    t === nothing && continue
    ra = _tap_step_factors(t.ratio_values, t.ratio_factors, br.tap_ratio)
    ph = _tap_step_factors(t.phase_values, t.phase_factors, br.phase_shift_deg)
    f = (ra[1] * ph[1], ra[2] * ph[2], ra[3] * ph[3], ra[4] * ph[4])
    rebase(value, b0, fk) = (fk == 0.0 || value == b0 * fk) ? b0 : value / fk
    br.tap_step_impedance = TapStepImpedance(rebase(br.r_pu, t.r_pu, f[1]), rebase(br.x_pu, t.x_pu, f[2]), rebase(br.g_from_pu, t.g_from_pu, f[3]),
                                             rebase(br.b_from_pu, t.b_from_pu, f[4]), rebase(br.g_to_pu, t.g_to_pu, f[3]), rebase(br.b_to_pu, t.b_to_pu, f[4]),
                                             t.ratio_values, t.ratio_factors, t.phase_values, t.phase_factors)
  end
  return nothing
end

# Dead time (qsts.dead_time_steps): a controller with lock steps left gets
# locked = true on the copy. run_control! keeps it in its rank group,
# evaluates it and takes no proposal from it; outside its band and not at
# its limit it holds the controllers with a larger priority number, in its
# band or at its limit it lets them act. Returns the controllers of `ctrls`
# the step runs locked.
function _qsts_lock!(ctrls::Vector{AbstractOuterController}, lock_steps::Vector{Int})::BitVector
  locked = BitVector([s > 0 for s in lock_steps])
  for c in findall(locked)
    ctrls[c].locked = true
  end
  return locked
end

# The start voltages of step k: the solution of the last solved step (or
# with sequential_start = snapshot the snapshot solution of step k); a bus
# without a value or at 0 pu there (isolated) takes the base solution.
function _qsts_seq_start(sctx::_QstsSeqContext, carry::_QstsCarry, k::Int)
  bvm, bva = sctx.base_warm
  if !isempty(sctx.snapshot_starts)
    vm, va = sctx.snapshot_starts[k]
  elseif carry.last_converged > 0
    vm, va = carry.vm, carry.va
  else
    return sctx.base_warm
  end
  out_vm = copy(bvm)
  out_va = copy(bva)
  for i in eachindex(vm)
    v = vm[i]
    (v === nothing || v == 0.0) && continue
    out_vm[i] = v
    out_va[i] = va[i]
  end
  return (out_vm, out_va)
end

# the largest absolute difference over the entries where both have a value
function _qsts_max_change(now::AbstractVector{Float64}, before::AbstractVector{Float64})::Float64
  out = NaN
  for i in eachindex(now, before)
    (isnan(now[i]) || isnan(before[i])) && continue
    d = abs(now[i] - before[i])
    out = isnan(out) ? d : max(out, d)
  end
  return out
end

# One step: a fresh copy of T0, the start voltages, the full delta, the
# reference positions read, the tap tables rebased, then the carry written
# (or, in a resync step, the positions of the series kept), dead time
# applied, one solve with the controllers. A converged step gives the
# next carry; a step that does not converge holds the old one, and only its
# dead time runs on.
function _qsts_seq_step(sctx::_QstsSeqContext, T0::Net, k::Int, carry::_QstsCarry)::Tuple{_QstsSeqOut,_QstsCarry}
  q = sctx.q
  input = sctx.input
  resync = k == 1 || (q.resync_every > 0 && (k - 1) % q.resync_every == 0)
  work = deepcopy(T0)
  _restore_start_voltages!(work, _qsts_seq_start(sctx, carry, k))
  apply!(work, input.ops[k], input.index)
  ctrls = sctx.controls ? _qsts_carried_controllers(work) : AbstractOuterController[]
  nc = length(ctrls)
  ref_taps, ref_banks = _qsts_read_state(work, ctrls, sctx.ctrl_branches)
  ref_pos = _qsts_positions(work, ctrls, sctx.ctrl_branches, ref_taps, ref_banks)
  _qsts_rebase_tap_table!(work, sctx.ctrl_branches)
  carried = carry.last_converged > 0
  # a resync step clears the dead time and the hunting memory
  lock_steps = resync ? zeros(Int, nc) : copy(carry.lock_steps)
  last_moves = resync ? fill(((0, 0), (0, 0)), nc) : copy(carry.last_moves)
  if resync
    _qsts_sync_model_steps!(work, ctrls, sctx.ctrl_branches)
  elseif carried
    _qsts_write_carry!(work, ctrls, sctx.ctrl_branches, carry)
  end
  # the position a move starts from: the carried one, before the first
  # solved step the one of the series
  old_pos = carried ? carry.positions : ref_pos
  outside = BitVector([ctrl isa PowerTransformerControl && _ratio_outside_control_range(ctrl, work.branchVec[sctx.ctrl_branches[c][1]]) for (c, ctrl) in enumerate(ctrls)])
  locked = _qsts_lock!(ctrls, lock_steps)
  # isolation and the start of isolated link members as in a snapshot step
  # (_qsts_step): an SCF or event step is marked again, a CGMES delta
  # carries its isolation
  if input.mark_isolated
    markIsolatedBuses!(net = work, log = false, links = true)
  else
    _start_isolated_link_members!(work)
  end
  # the islands of the step, those without a source included; then such
  # islands are isolated (F2), which a CGMES delta carries itself
  n_islands = length(first(_power_flow_islands(work)))
  shed = input.mark_isolated ? _qsts_isolate_sourceless!(work) : Dict{Int,Float64}()
  frozen = _qsts_frozen_buses(work)
  ctx = (sctx.base_solved || carried || !isempty(sctx.snapshot_starts)) ? sctx.warm : sctx.cold
  ex, msg = _qsts_solve!(work, ctx.step_config; print_island_table = false)
  solved = ex !== nothing && ex.erg == 0
  solved && calcNetLosses!(work)
  out = _qsts_capture(ctx, work, k, ex, msg, frozen)
  cvals = fill(NaN, nc, length(QSTS_CONTROLLER_QUANTITIES))
  rows = _QstsSwitchRow[]
  nb = length(work.nodeVec)
  nbr = length(work.branchVec)
  if !solved
    held = _QstsCarry(k, carry.taps, carry.banks, carry.positions, carry.vm, carry.va, fill(NaN, nb), fill(NaN, nbr), carry.switch_count, carry.n_hunting,
                      Int[max(s - 1, 0) for s in lock_steps], last_moves, carry.last_converged)
    return _QstsSeqOut(out, cvals, rows, resync), held
  end
  new_taps, new_banks = _qsts_read_state(work, ctrls, sctx.ctrl_branches)
  new_pos = _qsts_positions(work, ctrls, sctx.ctrl_branches, new_taps, new_banks)
  switch_count = copy(carry.switch_count)
  n_hunting = copy(carry.n_hunting)
  n_moves = 0
  t = 1
  for (c, ctrl) in enumerate(ctrls)
    has_pos = new_pos[c] != typemin(Int)
    steps = has_pos ? new_pos[c] - old_pos[c] : 0
    trigger = resync ? :resync : (outside[c] ? :range : :control)
    steps == 0 || push!(rows, (step = k, element_id = sctx.ctrl_ids[c], old = old_pos[c], new = new_pos[c], steps = steps, trigger = trigger))
    if trigger === :control && steps != 0
      switch_count[c] += abs(steps)
      n_moves += abs(steps)
      hit, last_moves[c] = _qsts_hunting_step(last_moves[c], k, Int(sign(steps)), q.hunting_window)
      hit && (n_hunting[c] += 1)
      lock_steps[c] = q.dead_time_steps
    else
      lock_steps[c] = max(lock_steps[c] - 1, 0)
    end
    cvals[c, _QC_POSITION] = has_pos ? new_pos[c] : NaN
    cvals[c, _QC_REF] = ref_pos[c] == typemin(Int) ? NaN : ref_pos[c]
    cvals[c, _QC_MOVES] = resync ? 0.0 : steps
    cvals[c, _QC_AT_LIMIT] = _qsts_at_range_end(work, ctrl, sctx.ctrl_branches[c]) ? 1.0 : 0.0
    cvals[c, _QC_LOCKED] = locked[c] ? 1.0 : 0.0
    if ctrl isa ShuntVoltageControl
      cvals[c, _QC_BS] = new_banks[c]
      has_pos && (cvals[c, _QC_DPOS_REF] = new_pos[c] - ref_pos[c])
    else
      cvals[c, _QC_DPOS_REF] = _qsts_dpos_ref(ctrl, work.branchVec[sctx.ctrl_branches[c][1]], ref_taps[t], new_taps[t])
      t += length(sctx.ctrl_branches[c])
    end
  end
  # the step quantities of mode B
  stepv = out.step
  q_vm = out.bus[:, _QB_VM]
  q_p = out.branch[:, _QR_P_FROM]
  stepv[_QS_RESYNC] = resync ? 1.0 : 0.0
  stepv[_QS_N_MOVES] = n_moves
  stepv[_QS_N_ISLANDS] = n_islands
  frozen_set = Set(frozen)
  isolated = Set(i for i in eachindex(q_vm) if isnan(q_vm[i]) && !(i in frozen_set))
  stepv[_QS_N_ISOLATED] = length(isolated)
  # the load of the buses without a value (load sign, consumption positive):
  # what the isolation of a sourceless island zeroed, else what the loads of
  # an isolated bus still carry
  bus_load = copy(shed)
  for ps in work.prosumpsVec
    b = Int(ps.comp.cFrom_bus)
    (b in isolated && !haskey(shed, b) && !isGenerator(ps)) || continue
    bus_load[b] = get(bus_load, b, 0.0) + something(ps.pVal, 0.0)
  end
  stepv[_QS_SHED_LOAD] = sum((v for (b, v) in bus_load if b in isolated); init = 0.0)
  if carry.last_converged == k - 1 && k > 1
    stepv[_QS_MAX_DVM_STEP] = _qsts_max_change(q_vm, carry.prev_vm)
    stepv[_QS_MAX_DP_STEP] = _qsts_max_change(q_p, carry.prev_p)
  end
  rep_bus = sctx.bus_reported
  stepv[_QS_N_VM_LOW] = count(i -> rep_bus[i] && !isnan(q_vm[i]) && q_vm[i] < q.vm_min_pu, eachindex(q_vm))
  stepv[_QS_N_VM_HIGH] = count(i -> rep_bus[i] && !isnan(q_vm[i]) && q_vm[i] > q.vm_max_pu, eachindex(q_vm))
  loading = out.branch[:, _QR_LOADING]
  stepv[_QS_N_OVERLOAD] = count(v -> !isnan(v) && v > q.loading_max_pct, loading)
  vm, va = _snapshot_start_voltages(work)
  next = _QstsCarry(k, new_taps, new_banks, new_pos, vm, va, q_vm, q_p, switch_count, n_hunting, lock_steps, last_moves, k)
  return _QstsSeqOut(out, cvals, rows, resync), next
end

# F2 in a step of an SCF or event series: a part of the net that holds no
# source (a unit with P > 0 or voltage regulation; a load never, a unit
# that draws power neither) cannot be solved and must not freeze at its
# start values. As in the import pass (_deEnergizeSourcelessComponents!) its
# injections are zeroed, its branches go out of service and its buses are
# isolated, so the step has no values there; without the zeroing the solver
# still sees the loads and refuses each bus as an island without reference.
# Connectivity runs over in-service branches and closed links. Returns the
# P the de-energized part consumed (loads, and units that drew power) per
# bus, for shed_load_mw.
function _qsts_isolate_sourceless!(net::Net)::Dict{Int,Float64}
  n = length(net.nodeVec)
  parent = collect(1:n)
  function find(i)
    while parent[i] != i
      parent[i] = parent[parent[i]]
      i = parent[i]
    end
    return i
  end
  unite(a, b) = (ra = find(a); rb = find(b); ra == rb || (parent[max(ra, rb)] = min(ra, rb)); nothing)
  for br in net.branchVec
    br.status == 1 && unite(Int(br.fromBus), Int(br.toBus))
  end
  for l in net.linkVec
    l.status == 1 && unite(Int(l.fromBus), Int(l.toBus))
  end
  sourced = falses(n)
  for ps in net.prosumpsVec
    (isGenerator(ps) && (something(ps.pVal, 0.0) > 0.0 || isRegulating(ps))) && (sourced[find(Int(ps.comp.cFrom_bus))] = true)
  end
  dead = Set(i for i in 1:n if !sourced[find(i)] && !isIsolated(net.nodeVec[i]))
  shed = Dict{Int,Float64}()
  isempty(dead) && return shed
  for ps in net.prosumpsVec
    b = Int(ps.comp.cFrom_bus)
    b in dead || continue
    p = something(ps.pVal, 0.0)
    shed[b] = get(shed, b, 0.0) + (isGenerator(ps) ? -p : p)
    ps.pVal = 0.0
    ps.qVal = 0.0
  end
  # all three flags, as the status op writes them: a reader that counts
  # branches by their terminal flags would still see one with status 0 alone
  for br in net.branchVec
    (br.fromBus in dead || br.toBus in dead) || continue
    br.status = 0
    br.from_status = 0
    br.to_status = 0
  end
  for i in sort!(collect(dead))
    nd = net.nodeVec[i]
    setNodeType!(nd, "Isolated")
    setVmVa!(node = nd, vm_pu = 0.0, va_deg = 0.0)
    push!(net.isoNodes, i)
  end
  sort!(unique!(net.isoNodes))
  return shed
end

# the carry before step 1: nothing solved, nothing counted
function _qsts_empty_carry(nc::Int)::_QstsCarry
  return _QstsCarry(0, _QstsTapPos[], fill(NaN, nc), fill(typemin(Int), nc), Union{Nothing,Float64}[], Union{Nothing,Float64}[], Float64[], Float64[],
                    zeros(Int, nc), zeros(Int, nc), zeros(Int, nc), fill(((0, 0), (0, 0)), nc), 0)
end

# The controller elements of a sequential run: address <actuator>#<kind>,
# the actuator a transformer (its branch address), the shunt of a bank when
# the series addresses it (a CGMES bank controller), else the bank's bus.
# Two controllers on one actuator with the same kind would share an
# address; that is refused with both names.
function _qsts_controller_elements(T0::Net, index::ScenarioIndex, e::QstsElements, ctrls::Vector{AbstractOuterController}, ctrl_branches::Vector{Vector{Int}})::QstsElements
  addr_of = Dict{Int,String}(id => a for (a, id) in index.id_by_mrid)
  shunt_addr = Dict{Int,String}()
  for id in sort!(collect(keys(index.kind_by_id)))
    index.kind_by_id[id] === :shunt || continue
    i = get(index.internal_by_id, id, 0)
    haskey(shunt_addr, i) || (shunt_addr[i] = get(addr_of, id, string(id)))
  end
  ids = String[]
  names = String[]
  kinds = Symbol[]
  actuators = String[]
  targets = String[]
  for (c, ctrl) in enumerate(ctrls)
    if ctrl isa PowerTransformerControl
      kind = ctrl.mode === :voltage ? :OLTC : (ctrl.mode === :branch_active_power ? :PST : :SCHRAEG)
      actuator = e.branch_ids[ctrl_branches[c][1]]
      target = ctrl.target_bus === nothing ? "" : e.bus_ids[geNetBusIdx(net = T0, busName = ctrl.target_bus)]
    else
      bank = ctrl::ShuntVoltageControl
      kind = bank.step_mvar === nothing ? :SVC : :MSC
      target = e.bus_ids[geNetBusIdx(net = T0, busName = bank.bus)]
      actuator = get(shunt_addr, bank.shunt_idx, target)
    end
    id = string(actuator, "#", kind)
    k = findfirst(==(id), ids)
    k === nothing || throw(ArgumentError(string("runQsts: controllers ", names[k], " and ", control_name(ctrl), " share the address ", id, "; a sequential run needs one controller per actuator and kind")))
    push!(ids, id)
    push!(names, control_name(ctrl))
    push!(kinds, kind)
    push!(actuators, actuator)
    push!(targets, target)
  end
  return QstsElements(e.bus_ids, e.bus_names, e.bus_vn_kv, e.bus_reported, e.branch_ids, e.branch_names, e.branch_kind, e.branch_from, e.branch_to, e.branch_rating_mva,
                      ids, names, kinds, actuators, targets)
end

# The step configuration of a sequential run: the one of a snapshot run with
# its controllers (qsts.snapshot_controls has no say in mode B), with the
# island parallelism of runtime.parallel when the run has one series.
function _qsts_seq_step_config(config::SparlectraConfig, warm::Bool, island_parallel::Bool)::Tuple{SparlectraConfig,Vector{String}}
  with_controls = _sparlectra_config_with(config; qsts = _struct_with(config.qsts; snapshot_controls = true))
  step_cfg, notes = _qsts_step_config(with_controls, warm)
  if island_parallel
    par = config.runtime.parallel
    step_cfg = _sparlectra_config_with(step_cfg; runtime = _struct_with(step_cfg.runtime; parallel = _struct_with(step_cfg.runtime.parallel; enabled = par.enabled)))
    filter!(t -> !startswith(t, "step configuration: runtime.parallel.enabled"), notes)
  end
  return step_cfg, notes
end

# the steps of one series in order, from the step after the carry's to `stop`
function _qsts_seq_loop(sctx::_QstsSeqContext, T0::Net, carry::_QstsCarry, stop::Int, @nospecialize(progress))::Tuple{Vector{_QstsSeqOut},_QstsCarry}
  outs = _QstsSeqOut[]
  for k in (carry.step + 1):stop
    # a QSTS run is long; the abort token of a Web UI job is checked per step
    sparlectra_check_abort()
    out, carry = _qsts_seq_step(sctx, T0, k, carry)
    push!(outs, out)
    progress === nothing || _progress_tick!(progress)
  end
  return outs, carry
end

# One series: T0 with the configured controllers, the validation, the base
# solve, the optional snapshot run for the starts, then the steps from the
# step after `initial` (nothing: from step 1) to `stop`. The restart test
# passes the carry of an earlier run; every other caller runs all steps.
_run_qsts_series(input::_QstsInput, config::SparlectraConfig, opts::_QstsSeqOptions)::_QstsSeriesRun = _run_qsts_series(input, config, opts, nothing, length(input.ops))
function _run_qsts_series(input::_QstsInput, config::SparlectraConfig, opts::_QstsSeqOptions, initial::Union{Nothing,_QstsCarry}, stop::Int)::_QstsSeriesRun
  n = length(input.ops)
  n >= 1 || throw(ArgumentError("runQsts: the series has no step"))
  1 <= stop <= n || throw(ArgumentError(string("runQsts: stop step ", stop, " outside 1 to ", n)))
  # a flat start is a flat start on the base case and on every step
  config, flat_overrides = _flatstart_profile_config(config)
  q = config.qsts
  (q.sequential_start === :snapshot && !isempty(input.event_list)) &&
    throw(ArgumentError(string("runQsts: qsts.sequential_start = snapshot needs the snapshot of each step; ", input.name, " is an event series without snapshots, use previous")))
  controls = config.control.enabled
  # the configured controllers go onto T0 once (a bank controller appends a
  # shunt), so every step has the same set in the same order
  T0 = deepcopy(input.template)
  controls && !isempty(config.control.controllers) && applyConfiguredControllers!(T0, config.control)
  # a tap_pos op on a regulated transformer is the reference of its step
  # (validation with regulated = false); outside resync steps the carried
  # tap replaces it
  scenarios = Scenario[Scenario(name = string("step ", k), ops = input.ops[k]) for k in 1:n if !isempty(input.ops[k])]
  _validate_scenarios_net(ScenarioSet(scenarios = scenarios), input.index, T0, false)
  notes = String[]
  _qsts_flatstart_note!(notes, flat_overrides)
  n_ref_taps = count(op -> op.op === :set && op.field === :tap_pos && (i = get(input.index.internal_by_id, op.id, 0); 1 <= i <= length(T0.branchVec) && _scenario_branch_tap_controller(T0, i) !== nothing),
                     Iterators.flatten(input.ops))
  n_ref_taps > 0 && push!(notes, string(n_ref_taps, " tap_pos op(s) on transformers with a tap controller are reference positions; outside resync steps the carried tap replaces them"))
  ctrls = controls ? _qsts_carried_controllers(T0) : AbstractOuterController[]
  controls || push!(notes, "control.enabled = false: a sequential run without controllers, taps and banks stay at the values of the series")
  n_other = controls ? length(collect_outer_controllers(T0)) - length(ctrls) : 0
  n_other > 0 && push!(notes, string(n_other, " outer-loop controller(s) other than tap and bank controllers start every step from the template; only taps and banks are carried"))
  ctrl_branches = Vector{Int}[ctrl isa PowerTransformerControl ? Int[findfirst(x -> x === _find_trafo_branch(T0, name), T0.branchVec) for name in vcat([ctrl.trafo], ctrl.followers)] : Int[]
                              for ctrl in ctrls]
  elements = _qsts_controller_elements(T0, input.index, _qsts_elements(T0, input.index), ctrls, ctrl_branches)
  warnings = String[]
  B = deepcopy(T0)
  base_frozen = _qsts_frozen_buses(B)
  start_warnings = _TapStartWarnings()
  ex, msg = Base.ScopedValues.with(() -> _qsts_solve!(B, config; print_island_table = opts.print_island_table), _TAP_START_WARNINGS => start_warnings)
  base_status = _qsts_status(ex)
  base_solved = base_status === :solved
  if base_solved
    calcNetLosses!(B)
    base_warm = _snapshot_start_voltages(B)
    w = _plausibility_diagnostics(B, config.powerflow.plausibility, true).plausibility_warning
    isempty(w) || (push!(notes, string("QSTS base case: ", w)); push!(warnings, string("QSTS base case: ", w)))
  else
    base_warm = _snapshot_start_voltages(T0)
    push!(notes, string("base case ", base_status === :solver_error ? string("failed: ", msg) : string("did not converge", isempty(msg) ? "" : string(": ", msg)), "; the steps start from the start values of the template"))
  end
  base_iterations, base_outer, base_control = _qsts_iterations(B, ex, controls)
  base_label = input.base_step >= 1 ? input.labels[input.base_step] : "base"
  base_time = input.base_step >= 1 ? input.timestamps[input.base_step] : ""
  base = QstsStepInfo(0, base_label, base_time, true, base_status, base_control, base_iterations, base_outer, 0, length(base_frozen),
                      base_solved ? _qsts_slack_buses(B) : Int[], msg)
  cold_cfg, cold_notes = _qsts_seq_step_config(config, false, opts.island_parallel)
  warm_cfg, warm_notes = _qsts_seq_step_config(config, true, opts.island_parallel)
  append!(notes, base_solved ? warm_notes : cold_notes)
  base_solved || cold_cfg.powerflow.start_mode.flatstart == warm_cfg.powerflow.start_mode.flatstart ||
    push!(notes, "from the first converged step on: power_flow.start_mode.flatstart = false (warm start from the step before)")
  starts = Tuple{Vector{Union{Nothing,Float64}},Vector{Union{Nothing,Float64}}}[]
  if q.sequential_start === :snapshot
    # the snapshot run of the same input without controllers; only its bus
    # voltages are read, so its notes and warnings stay out of this run
    a_cfg = _sparlectra_config_with(config; qsts = _struct_with(q; mode = :snapshots, snapshot_controls = false))
    a_opts = _QstsRunOptions(opts.island_parallel ? opts.parallel_enabled : false, opts.parallel_max_tasks, opts.parallel_min_work_items, nothing)
    ra = Logging.with_logger(() -> _run_qsts(input, a_cfg, a_opts), Logging.NullLogger())
    value(x) = isnan(x) ? nothing : x
    for k in 1:n
      push!(starts, (Union{Nothing,Float64}[value(x) for x in ra.bus_values[:, _QB_VM, k]], Union{Nothing,Float64}[value(x) for x in ra.bus_values[:, _QB_VA, k]]))
    end
    push!(notes, "start voltages: the snapshot run of each step without controllers (qsts.sequential_start = snapshot)")
  end
  tap_branch = _qsts_tap_branches(T0)
  sctx = _QstsSeqContext(input, q, _QstsContext(input, cold_cfg, nothing, controls, tap_branch), _QstsContext(input, warm_cfg, nothing, controls, tap_branch),
                         base_solved, base_warm, starts, ctrl_branches, elements.controller_ids, elements.bus_reported, controls)
  carry = initial === nothing ? _qsts_empty_carry(length(ctrls)) : initial
  outs, carry = Base.ScopedValues.with(() -> _qsts_seq_loop(sctx, T0, carry, stop, opts.progress), _TAP_START_WARNINGS => start_warnings)
  return _QstsSeriesRun(input, config, elements, base, outs, carry, notes, warnings)
end

# The QstsResult of a whole series run: values in step order, the
# plausibility notes in step order (their warnings go to `warnings` for the
# orchestrator), the switching rows, the metrics.
function _qsts_seq_result(run::_QstsSeriesRun)::Tuple{QstsResult,Vector{String}}
  input = run.input
  config = run.config
  n = length(input.ops)
  length(run.outs) == n || error("runQsts: the series run holds ", length(run.outs), " of ", n, " steps")
  e = run.elements
  nb = length(e.bus_ids)
  nbr = length(e.branch_ids)
  nc = length(e.controller_ids)
  bus_values = Array{Float64,3}(undef, nb, length(QSTS_BUS_QUANTITIES), n)
  branch_values = Array{Float64,3}(undef, nbr, length(QSTS_BRANCH_QUANTITIES), n)
  step_values = Matrix{Float64}(undef, length(QSTS_STEP_QUANTITIES), n)
  controller_values = Array{Float64,3}(undef, nc, length(QSTS_CONTROLLER_QUANTITIES), n)
  steps = Vector{QstsStepInfo}(undef, n)
  switching = _QstsSwitchRow[]
  notes = copy(run.notes)
  warnings = copy(run.warnings)
  for k in 1:n
    o = run.outs[k]
    (size(o.out.bus, 1) == nb && size(o.out.branch, 1) == nbr) || error("runQsts: step $(k) changed the bus or branch count")
    bus_values[:, :, k] = o.out.bus
    branch_values[:, :, k] = o.out.branch
    step_values[:, k] = o.out.step
    controller_values[:, :, k] = o.ctrl
    steps[k] = o.out.info
    append!(switching, o.switching)
    w = o.out.plausibility_warning
    if !isempty(w)
      label = o.out.info.label
      text = string("QSTS step ", k, isempty(label) ? "" : string(" (", label, ")"), ": ", w)
      push!(notes, text)
      push!(warnings, text)
    end
  end
  n_failed = count(s -> s.status !== :solved, steps)
  n_failed > 0 && push!(notes, string(n_failed, " step(s) did not converge; each held taps, banks and start voltages of the step before"))
  q = config.qsts
  resync = BitVector([o.resync for o in run.outs])
  metrics = _qsts_metric_rows(e, steps, resync, bus_values, branch_values, step_values, controller_values, run.carry.switch_count, run.carry.n_hunting, q,
                              _qsts_step_durations(input.step_minutes, input.timestamps))
  pf = config.powerflow
  r = QstsResult(input.name, input.origin, :sequential, input.base_step, input.step_minutes, run.base, steps, e, bus_values, branch_values, step_values, controller_values,
                 switching, metrics, pf.max_iter, pf.tol, config.qsts.snapshot_controls, config.cgmes.tap_control, config.cgmes.map_energy_source, q, input.event_list, notes)
  return r, warnings
end

# Several series: one task per series over the parallel map of the engine
# (the serial path is the same function), the steps of each series in
# order; the warnings are logged afterwards in series and step order. Each
# series builds its own T0 inside its task and no factorization is shared
# (power mode is off in the steps), so no UMFPACK copy is made under threads.
function _run_qsts_sequential(inputs::Vector{_QstsInput}, config::SparlectraConfig, opts::_QstsSeqOptions)::Vector{QstsResult}
  _qsts_check_config(config, :sequential)
  ns = length(inputs)
  par = config.runtime.parallel
  runs = _engine_parallel_map((_, s) -> _run_qsts_series(inputs[s], config, opts), () -> nothing, ns;
    parallel_enabled = something(opts.parallel_enabled, par.enabled), parallel_max_tasks = something(opts.parallel_max_tasks, parallel_max_tasks(par)),
    parallel_min_work_items = something(opts.parallel_min_work_items, par.min_work_items))
  results = QstsResult[]
  for s in 1:ns
    r, warnings = _qsts_seq_result(runs[s]::_QstsSeriesRun)
    for w in warnings
      @warn w
    end
    push!(results, r)
  end
  return results
end

"""
    runQsts(series::AbstractVector{SCFSeries}; config, parallel_enabled = nothing, parallel_max_tasks = nothing, parallel_min_work_items = nothing, progress = nothing) -> Vector{QstsResult}

QSTS on several SCF series, one [`QstsResult`](@ref) per series. With
`qsts.mode = snapshots` each series is a snapshot run (see
[`runQsts`](@ref)), one after the other. With `qsts.mode = sequential` it is
QSTS mode B; `runQsts(series::SCFSeries; config)` with that mode is the run
on one series.

Mode B solves the steps of a series in order, each on a fresh copy of the
template with the full delta of the step and the outer-loop controllers of
the `control` section running (`control.enabled = false` runs without them).
The state of the tap controllers (masters and followers) and the shunt bank
controllers goes from the last solved step to the next:

1. The base case is solved with the controllers; its solution is the start
   of step 1.
2. Every step: the start voltages (`qsts.sequential_start = previous`: the
   solution of the last solved step, `snapshot`: the solution of the same
   step in a snapshot run without controllers; a bus without a value there
   takes the base solution), the full delta (for a CGMES union without the
   snapshot start values of a snapshot run), the taps and banks the series
   gives the step read as reference, then the carried taps and banks written
   absolutely, then one solve with the controllers.
3. Resync steps (step 1 always, and with `qsts.resync_every = n` every step
   with `(k - 1) % n == 0`) keep the taps and banks of the series and clear
   the dead time and the hunting memory.
4. `qsts.dead_time_steps = d` locks a controller for the d steps after a
   counted move (`enabled = false` in those steps). A locked controller that
   does not stand at an end of its range also holds every controller with a
   larger priority number in that step.
5. A step that does not converge holds taps, banks and start voltages of
   the last solved step; its dead time runs on.

A `tap_pos` op on a transformer with a tap controller is the reference
position of its step (`notes` says so). The step configuration is the one
of a snapshot run with controllers, the island parallelism follows
`runtime.parallel` for one series and is off for several. The checks of a
snapshot run apply, with `qsts.mode = sequential`.

The result adds controller elements (`QstsElements.controller_ids`), the
values of [`QSTS_CONTROLLER_QUANTITIES`](@ref), the step quantities of mode
B in [`QSTS_STEP_QUANTITIES`](@ref), the switching list
([`qstsSwitching`](@ref)) and the metrics ([`qstsMetrics`](@ref)).

Several series run in parallel over Julia threads under `runtime.parallel`
(the keywords override it); the steps of one series run in order. Parallel
and serial runs give identical results. `progress(done, total)` is called
once per finished step over all series in mode B, once per step and series
in mode A.
"""
function runQsts(series::AbstractVector{SCFSeries}; config::SparlectraConfig, parallel_enabled::Union{Nothing,Bool} = nothing, parallel_max_tasks::Union{Nothing,Integer} = nothing,
                 parallel_min_work_items::Union{Nothing,Integer} = nothing, @nospecialize(progress = nothing))::Vector{QstsResult}
  if config.qsts.mode !== :sequential
    return QstsResult[runQsts(s; config = config, parallel_enabled = parallel_enabled, parallel_max_tasks = parallel_max_tasks, parallel_min_work_items = parallel_min_work_items,
                              progress = progress) for s in series]
  end
  _qsts_check_config(config, :sequential)
  return _run_qsts_sequential(_QstsInput[_qsts_scf_input(s, config) for s in series], config, parallel_enabled, parallel_max_tasks, parallel_min_work_items, progress)
end

# the keywords of a vector run behind the options struct: one progress
# counter over the steps of all series, island parallelism and the island
# table only for a single series
function _run_qsts_sequential(inputs::Vector{_QstsInput}, config::SparlectraConfig, parallel_enabled::Union{Nothing,Bool}, parallel_max_tasks::Union{Nothing,Integer},
                              parallel_min_work_items::Union{Nothing,Integer}, @nospecialize(progress))::Vector{QstsResult}
  one = length(inputs) == 1
  tracker = _BatchProgress(progress, sum((length(i.ops) for i in inputs); init = 0))
  opts = _QstsSeqOptions(parallel_enabled, parallel_max_tasks === nothing ? nothing : Int(parallel_max_tasks), parallel_min_work_items === nothing ? nothing : Int(parallel_min_work_items),
                         tracker, one, one)
  return _run_qsts_sequential(inputs, config, opts)
end

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
#
# This file is included inside module Sparlectra. Do not add a module wrapper here.
#
# Rectangular power-flow per-iteration Q-limit active-set helper.

# Date: 29.5.2026
# file: src/powerflow_rectangular/rectangular_qlimit_iteration.jl
# purpose: per-iteration Q-limit active-set handling inside the rectangular
#          Newton loop: PV/PQ switching with hysteresis, cooldown, and guard
#          limits on repeated switching

function _handle_rectangular_qlimit_iteration!(
  net,
  it,
  nb,
  Ybus,
  V,
  S,
  bus_types,
  Vset,
  slack_idx,
  workspace,
  history,
  qmin_pu,
  qmax_pu,
  pv_orig_mask,
  qlimit_mode,
  allow_reenable,
  q_hyst_pu,
  cooldown_iters,
  reenable_v_hyst_pu,
  lock_pv_to_pq_buses,
  qlimit_start_mode,
  qlimit_start_iter,
  qlimit_auto_q_delta_pu,
  qlimit_trace_enabled,
  qlimit_trace_internal,
  qlimit_lock_reason,
  qlimit_max_outer,
  controllers,
  base_vset,
  adjust_counter,
  qlimit_guard_max_switches,
  qlimit_guard_freeze_after_repeated_switching,
  qlimit_guard_violation_mode,
  qlimit_guard_violation_threshold_pu,
  tol,
  verbose,
  max_console_rows,
  performance_profile,
  dslack::Union{Nothing,DistributedSlackState} = nothing,
)
  changed = false
  reenabled = false

  # qreq vector is tracked across iterations for :auto_q_delta start logic.
  _perf_profile_time!(performance_profile, :iteration_qload) do
    copyto!(workspace.qload_pu, build_qload_pu(net))
  end
  Scalc_pu = _perf_profile_time!(performance_profile, :iteration_calc_injections) do
    calc_injections(Ybus, V)
  end
  current_pv_qreq_pu = _perf_profile_time!(performance_profile, :iteration_qreq_vector) do
    fill!(workspace.current_pv_qreq_pu, NaN)
    @inbounds for bus in eachindex(workspace.current_pv_qreq_pu)
      if bus_types[bus] == :PV
        workspace.current_pv_qreq_pu[bus] = imag(Scalc_pu[bus]) + workspace.qload_pu[bus]
      end
    end
    workspace.current_pv_qreq_pu
  end

  qlimit_iter_ready, qlimit_auto_ready, qlimit_ready, converged_this_iter, qlimit_check_active = _perf_profile_time!(performance_profile, :iteration_control_bookkeeping) do
    qlimit_iter_ready_ = it >= qlimit_start_iter
    qlimit_auto_ready_ = false
    if qlimit_start_mode in (:auto, :iteration_or_auto)
      max_q_delta = 0.0
      compared = false
      @inbounds for bus in eachindex(current_pv_qreq_pu)
        if isfinite(current_pv_qreq_pu[bus]) && isfinite(workspace.prev_pv_qreq_pu[bus])
          max_q_delta = max(max_q_delta, abs(current_pv_qreq_pu[bus] - workspace.prev_pv_qreq_pu[bus]))
          compared = true
        end
      end
      qlimit_auto_ready_ = compared && (max_q_delta <= qlimit_auto_q_delta_pu)
    end
    converged_this_iter_ = history[end] <= tol
    # A converged iterate is always "ready": the start gating exists to keep
    # switching off unsettled early iterates, and a converged state is the most
    # settled state there is. Without this, a solve that converges FASTER than
    # the gate allows (measured: the ReliCapGrid CGM island converges in 7
    # iterations, :auto never fires) exits with its violations unevaluated and
    # is then rejected as `remaining_pv_q_limit_violations`: the active-set
    # path had no post-convergence switching round. Switching on the converged
    # iterate is safe: the solver accepts convergence only when the active set
    # did NOT change in that iteration, so the loop simply continues under the
    # existing hysteresis/cooldown/guard machinery until it converges with a
    # stable active set.
    qlimit_ready_ = (qlimit_start_mode == :iteration ? qlimit_iter_ready_ : qlimit_start_mode == :auto ? qlimit_auto_ready_ : (qlimit_iter_ready_ || qlimit_auto_ready_)) || converged_this_iter_
    qlimit_check_active_ = qlimit_ready_
    (qlimit_iter_ready_, qlimit_auto_ready_, qlimit_ready_, converged_this_iter_, qlimit_check_active_)
  end

  if qlimit_trace_enabled
    _perf_profile_time!(performance_profile, :iteration_qlimit_trace) do
      lock_mask_trace = falses(nb)
      for bus in lock_pv_to_pq_buses
        if 1 <= bus <= nb
          lock_mask_trace[bus] = true
        end
      end
      for bus in qlimit_trace_internal
        bus <= nb || continue
        qreq = bus_types[bus] in (:PV, :Slack) ? imag(Scalc_pu[bus]) + workspace.qload_pu[bus] : NaN
        _print_rectangular_qlimit_trace(
          stdout,
          net,
          it,
          bus,
          bus_types[bus],
          V,
          Vset,
          qreq,
          qmin_pu,
          qmax_pu;
          q_hyst_pu = q_hyst_pu,
          qlimit_start_iter = qlimit_start_iter,
          qlimit_start_mode = qlimit_start_mode,
          qlimit_iter_ready = qlimit_iter_ready,
          qlimit_auto_ready = qlimit_auto_ready,
          qlimit_ready = qlimit_ready,
          qlimit_check_active = qlimit_check_active,
          converged_before_switching = converged_this_iter,
          qlimit_mode = qlimit_mode,
          lock_mask = lock_mask_trace,
          qlimit_lock_reason = qlimit_lock_reason,
        )
      end
    end
  end

  F = workspace.rhs_vector
  max_mis = history[end]
  if qlimit_check_active
    # Active-set operation may mutate bus_types and S; mismatch must be recomputed.
    changed, reenabled = _perf_profile_time!(performance_profile, :iteration_qlimit_active_set) do
      active_set_q_limits!(
        net,
        it,
        nb;
        qmin_pu = qmin_pu,
        qmax_pu = qmax_pu,
        pv_orig_mask = pv_orig_mask,
        allow_reenable = qlimit_mode == :switch_to_pq ? allow_reenable : false,
        q_hyst_pu = q_hyst_pu,
        cooldown_iters = cooldown_iters,
        # voltage-side release (#375): the current iterate and the setpoint
        get_vm_pu = bus -> abs(V[bus]),
        get_vset_pu = bus -> Vset[bus],
        v_hyst_pu = reenable_v_hyst_pu,
        # the one-retry guard gives a twice-clamped bus one more release on
        # a converged iterate (#475)
        converged_iterate = converged_this_iter,
        lock_pv_to_pq_buses = lock_pv_to_pq_buses,
        on_violation! = qlimit_mode == :adjust_vset ? ((bus, qreq, side, qclamp) -> _try_adjust_vset_on_q_limit!(net, bus, side, it, controllers, base_vset, Vset, adjust_counter, qlimit_max_outer, verbose)) : nothing,
        verbose = verbose,
        max_console_rows = max_console_rows,
        get_qreq_pu = bus -> begin
          (bus_types[bus] == :Slack) && return 0.0
          return imag(Scalc_pu[bus]) + workspace.qload_pu[bus]
        end,
        is_pv = bus -> (bus_types[bus] == :PV),
        make_pq! = (bus, qclamp_gen_pu, side) -> begin
          bus_types[bus] = :PQ
          qinj_pu = qclamp_gen_pu - workspace.qload_pu[bus]
          S[bus] = ComplexF64(real(S[bus]), qinj_pu)
          net.nodeVec[bus]._qƩGen = qclamp_gen_pu * net.baseMVA
        end,
        make_pv! = (bus) -> begin
          bus_types[bus] = :PV
        end,
        qlimit_guard_max_switches = qlimit_guard_max_switches,
        qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching,
        qlimit_guard_violation_mode = qlimit_guard_violation_mode,
        qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu,
        lock_mask_buffer = workspace.lock_mask,
      )
    end

    if changed || reenabled
      # Recompute at the ACCEPTED lambda: switching changed bus_types/S, not
      # the distributed-slack state (participation is frozen per solve).
      _dslack_set_trial!(dslack, 0.0)
      F = _perf_profile_time!(performance_profile, :iteration_mismatch_after_qlimit) do
        mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx; dslack = dslack)
      end
      max_mis = maximum(abs.(F))
      history[end] = max_mis
    end
  end

  copyto!(workspace.prev_pv_qreq_pu, current_pv_qreq_pu)

  return (; F, max_mis, changed, reenabled, converged_this_iter)
end

# Generator units of a bus whose Q the bus Q band takes at a fixed value
# that follows the unit's own Q(U) characteristic: a non-regulating,
# non-reference unit with a QUController on a bus held by a regulating or
# reference unit. Same rule as _buildQLimits!, which puts such a unit into
# the band with its Q at the held voltage.
@inline _qband_follows_characteristic(ps, held::Bool) = held && !isRegulating(ps) && !isSlack(ps) && has_qu_controller(ps)

"""
    _reapply_qlimit_clamps!(S, dQinj_dVm, net, bus_types, V) -> Int

Apply the Q-limit clamps again after `buildControlledSVec` rebuilt `S` for
the current voltage `V`. With a voltage-dependent controller in the net the
rectangular Newton loop rebuilds `S` every iteration from the units'
schedules, and the rebuild dropped every clamp written before (active set
`make_pq!`, the guard, `lock_pv_to_pq_buses`, the warm clamps): the
regulating unit of a clamped bus was back at its scheduled Q.

A clamped bus is a bus that started the solve as PV, is PQ in `bus_types`
and has an entry in `net.qLimitEvents`. Its clamp value is the generator Q
every clamp site writes to `node._qƩGen` (pu times `baseMVA`), in the terms
of the bus Q band: a Q(U) unit that follows its characteristic (see
`_buildQLimits!`) is in it with its Q at the held voltage. The clamp stands
for the other generator units of the bus; the followers keep their Q at
`V` from the rebuild. Per clamped bus:

`S_q = S_q(rebuild) - Q_other(V) + Q_clamp - Q_follow(Vheld)`

with `Q_other(V)` the Q the rebuild gave the other generator units. Loads,
voltage-dependent shunts and the P part stay as rebuilt. `dQinj_dVm` loses
the derivative of the other units that carry a Q(U) controller, since
their Q is fixed by the clamp.

Returns the number of clamped buses written.
"""
function _reapply_qlimit_clamps!(S::Vector{ComplexF64}, dQinj_dVm::Vector{Float64}, net::Net, bus_types::Vector{Symbol}, V::Vector{ComplexF64})
  isempty(net.qLimitEvents) && return 0
  nb = length(bus_types)
  nodes = net.nodeVec
  clamped = falses(nb)
  nclamped = 0
  for bus in keys(net.qLimitEvents)
    (1 <= bus <= nb && bus_types[bus] == :PQ) || continue
    # the node type is still the type the solve started with (the bus types
    # are synced to the net only after the loop); a bus that started PQ has
    # no clamp, its Q comes from its units
    getNodeType(nodes[bus]) == PV || continue
    nodes[bus]._qƩGen === nothing && continue
    clamped[bus] = true
    nclamped += 1
  end
  nclamped == 0 && return 0

  held = falses(nb)
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    bus = getPosumerBusIndex(ps)
    (1 <= bus <= nb && clamped[bus]) || continue
    (isRegulating(ps) || isSlack(ps)) && (held[bus] = true)
  end

  base = net.baseMVA
  q_other = zeros(Float64, nb)
  dq_other = zeros(Float64, nb)
  q_follow_held = zeros(Float64, nb)
  # held voltages as _buildQLimits! uses them, built only when a follower
  # sits on a clamped bus
  vheld = nothing
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    bus = getPosumerBusIndex(ps)
    (1 <= bus <= nb && clamped[bus]) || continue
    if _qband_follows_characteristic(ps, held[bus])
      vheld === nothing && (vheld = _bus_voltage_setpoints_from_prosumers(net))
      q_follow_held[bus] += first(evaluate_controller(ps.quController, vheld[bus]))
    elseif has_qu_controller(ps)
      # the same evaluation buildControlledSVec put into S
      q_pu, dq_pu = evaluate_controller(ps.quController, _safe_vm(abs(V[bus])))
      q_other[bus] += q_pu
      dq_other[bus] += dq_pu
    else
      q_other[bus] += something(ps.qVal, 0.0) / base
    end
  end

  @inbounds for bus in eachindex(clamped)
    clamped[bus] || continue
    q_clamp = nodes[bus]._qƩGen / base
    S[bus] = ComplexF64(real(S[bus]), imag(S[bus]) - q_other[bus] + q_clamp - q_follow_held[bus])
    dQinj_dVm[bus] -= dq_other[bus]
  end
  return nclamped
end

"""
    _clamped_bus_physical_qgen!(net, Sbus_MVA, q_consumed) -> Int

Write the physical generator Q (MVAr) to `node._qƩGen` of every clamped bus
(PQ with an entry in `net.qLimitEvents`) that carries a Q(U) unit following
its characteristic (see `_reapply_qlimit_clamps!`). The clamp record is the
band value, with the follower at its Q at the held voltage; the clamped bus
does not hold that voltage, so the follower runs at another Q. The physical
value is the solved injection plus what the bus consumes at its solved
voltage magnitude (`q_consumed`, MVAr per bus from `_bus_consumption_MVA`):
loads (a Q(U) load at its characteristic) and voltage-dependent shunts. A
clamped bus without such a unit keeps its clamp record, which is its
physical value already.

Called after the solved voltages are written to the nodes. Returns the
number of buses written.
"""
function _clamped_bus_physical_qgen!(net::Net, Sbus_MVA::AbstractVector, q_consumed::AbstractVector)
  isempty(net.qLimitEvents) && return 0
  nodes = net.nodeVec
  nb = length(nodes)
  clamped = falses(nb)
  for bus in keys(net.qLimitEvents)
    (1 <= bus <= nb && bus <= length(Sbus_MVA)) || continue
    bus in net.isoNodes && continue
    getNodeType(nodes[bus]) == PQ || continue
    clamped[bus] = true
  end
  any(clamped) || return 0

  held = falses(nb)
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    bus = getPosumerBusIndex(ps)
    (1 <= bus <= nb && clamped[bus]) || continue
    (isRegulating(ps) || isSlack(ps)) && (held[bus] = true)
  end
  has_follower = falses(nb)
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    bus = getPosumerBusIndex(ps)
    (1 <= bus <= nb && clamped[bus]) || continue
    _qband_follows_characteristic(ps, held[bus]) && (has_follower[bus] = true)
  end
  any(has_follower) || return 0

  written = 0
  for bus in eachindex(has_follower)
    has_follower[bus] || continue
    nodes[bus]._qƩGen = imag(Sbus_MVA[bus]) + q_consumed[bus]
    written += 1
  end
  return written
end

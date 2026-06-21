# Copyright 2023–2026 Udo Schmitz
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
    qlimit_ready_ = qlimit_start_mode == :iteration ? qlimit_iter_ready_ : qlimit_start_mode == :auto ? qlimit_auto_ready_ : (qlimit_iter_ready_ || qlimit_auto_ready_)
    converged_this_iter_ = history[end] <= tol
    violation_guard_active = qlimit_guard_violation_mode == :lock_pq
    qlimit_check_active_ = qlimit_ready_ && (!converged_this_iter_ || violation_guard_active)
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
      F = _perf_profile_time!(performance_profile, :iteration_mismatch_after_qlimit) do
        mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
      end
      max_mis = maximum(abs.(F))
      history[end] = max_mis
    end
  end

  copyto!(workspace.prev_pv_qreq_pu, current_pv_qreq_pu)

  return (; F, max_mis, changed, reenabled, converged_this_iter)
end

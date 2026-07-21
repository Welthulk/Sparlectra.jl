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
# file: src/powerflow_rectangular/rectangular_qlimit_outer_loop.jl
#

# Rectangular Q-limit outer-loop helpers.
#
# This file owns the MATPOWER-style classical Q-limit outer-loop orchestration.
# It must preserve PV→PQ/ref-bus behavior and diagnostics while keeping the
# Newton iteration driver focused on one rectangular solve attempt.

function _rectangular_qgen_requirements_pu(net::Net, Ybus, V::Vector{ComplexF64})
  Sbus_pu, _ = _compute_rectangular_final_injections(Ybus, V, net.baseMVA)
  qload_pu = build_qload_pu(net)
  return imag.(Sbus_pu) .+ qload_pu
end

function _matpower_q_limit_violations(net::Net, qreq_pu::AbstractVector, bus_types::AbstractVector{Symbol}, fixed_gens::BitVector; tolerance_pu::Float64)
  rows = NamedTuple[]
  for bus in eachindex(bus_types)
    bus_types[bus] in (:PV, :Slack) || continue
    gen_indices = _generators_at_bus(net, bus)
    active_gen_indices = [idx for idx in gen_indices if !fixed_gens[idx]]
    isempty(active_gen_indices) && continue
    # distribute_bus_generation! stores solved per-generator Q when possible.
    # If a network has not populated qRes for a generator, fall back to an equal
    # share of the bus-level regulating Q request so mode dispatch remains
    # deterministic for small synthetic tests.
    fallback_q_pu = qreq_pu[bus] / length(active_gen_indices)
    for gen_idx in active_gen_indices
      ps = net.prosumpsVec[gen_idx]
      qreq = isnothing(ps.qRes) ? fallback_q_pu : ps.qRes / net.baseMVA
      qmax = isnothing(ps.maxQ) ? Inf : ps.maxQ / net.baseMVA
      qmin = isnothing(ps.minQ) ? -Inf : ps.minQ / net.baseMVA
      if isfinite(qmax) && qreq > qmax + tolerance_pu
        push!(rows, (bus = bus, gen_index = gen_idx, side = :max, q_before_pu = qreq, q_limit_pu = qmax, violation_pu = qreq - qmax))
      elseif isfinite(qmin) && qreq < qmin - tolerance_pu
        push!(rows, (bus = bus, gen_index = gen_idx, side = :min, q_before_pu = qreq, q_limit_pu = qmin, violation_pu = qmin - qreq))
      end
    end
  end
  return rows
end

function _matpower_bus_has_regulating_generator(net::Net, bus::Int, fixed_gens::BitVector)
  for gen_idx in _generators_at_bus(net, bus)
    fixed_gens[gen_idx] && continue
    return true
  end
  return false
end

function _select_new_reference_bus!(net::Net, old_ref::Int)
  for (bus, node) in enumerate(net.nodeVec)
    getNodeType(node) == PV || continue
    setNodeType!(node, "Slack")
    @info "Classical Q-limit enforcement changed reference bus" old_ref = old_ref new_ref = bus
    return bus
  end
  return nothing
end

function _run_q_limits_matpower_outer_loop!(
  net::Net;
  mode::Symbol,
  method::Symbol,
  maxiter::Int,
  tol::Float64,
  damp::Float64,
  verbose::Int,
  autodamp::Bool,
  autodamp_min::Float64,
  opt_flatstart::Bool,
  pv_table_rows::Int,
  qlimit_max_outer::Int,
  start_projection::Bool,
  start_projection_try_dc_start::Bool,
  start_projection_try_blend_scan::Bool,
  start_projection_branch_guard::Bool,
  start_projection_measure_candidates::Bool,
  start_projection_accept_unmeasured_dc_start::Bool,
  start_projection_blend_lambdas,
  start_projection_dc_angle_limit_deg::Float64,
  start_projection_requested_angle_mode::Symbol,
  start_projection_requested_voltage_mode::Symbol,
  start_current_iteration_enabled::Bool,
  start_current_iteration_max_iter::Int,
  start_current_iteration_tol::Float64,
  start_current_iteration_damping::Float64,
  start_current_iteration_accept_only_if_improved::Bool,
  start_current_iteration_min_improvement_factor::Float64,
  start_current_iteration_vm_min_pu::Float64,
  start_current_iteration_vm_max_pu::Float64,
  start_current_iteration_max_angle_step_deg::Float64,
  start_current_iteration_only_for_large_cases::Bool,
  apslf_start_enabled::Bool,
  apslf_start_order::Int,
  wrong_branch_detection::Symbol,
  wrong_branch_rescue::Bool,
  wrong_branch_min_vm_pu::Float64,
  wrong_branch_max_vm_pu::Float64,
  wrong_branch_max_angle_spread_deg::Float64,
  wrong_branch_max_branch_angle_deg::Float64,
  wrong_branch_min_low_vm_count::Int,
  wrong_branch_rescue_max_attempts::Int,
  performance_profile,
  rectangular_workspace_reuse::Bool,
  rectangular_preallocate_workspace::Symbol,
  rectangular_workspace_min_buses::Int,
)
  mode in (:classic_simultaneous, :classic_one_at_a_time) || error("Unsupported classical Q-limit mode $(mode).")
  resetQLimitLog!(net)
  fixed_gens = falses(length(net.prosumpsVec))
  outer_rows = NamedTuple[]
  total_iters = 0
  base_pf_converged = false
  qlimit_enforcement_started = false
  final_outcome = :base_pf_not_converged

  for outer in 0:qlimit_max_outer
    iters, erg = runpf_rectangular!(
      net;
      method = method,
      maxiter = maxiter,
      tol = tol,
      damp = damp,
      verbose = verbose,
      autodamp = autodamp,
      autodamp_min = autodamp_min,
      opt_flatstart = outer == 0 ? opt_flatstart : false,
      pv_table_rows = pv_table_rows,
      qlimit_max_outer = qlimit_max_outer,
      start_projection = outer == 0 ? start_projection : false,
      start_projection_try_dc_start = start_projection_try_dc_start,
      start_projection_try_blend_scan = start_projection_try_blend_scan,
      start_projection_branch_guard = start_projection_branch_guard,
      start_projection_measure_candidates = start_projection_measure_candidates,
      start_projection_accept_unmeasured_dc_start = start_projection_accept_unmeasured_dc_start,
      start_projection_blend_lambdas = start_projection_blend_lambdas,
      start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
      start_projection_requested_angle_mode = start_projection_requested_angle_mode,
      start_projection_requested_voltage_mode = start_projection_requested_voltage_mode,
      start_current_iteration_enabled = outer == 0 ? start_current_iteration_enabled : false,
      start_current_iteration_max_iter = start_current_iteration_max_iter,
      start_current_iteration_tol = start_current_iteration_tol,
      start_current_iteration_damping = start_current_iteration_damping,
      start_current_iteration_accept_only_if_improved = start_current_iteration_accept_only_if_improved,
      start_current_iteration_min_improvement_factor = start_current_iteration_min_improvement_factor,
      start_current_iteration_vm_min_pu = start_current_iteration_vm_min_pu,
      start_current_iteration_vm_max_pu = start_current_iteration_vm_max_pu,
      start_current_iteration_max_angle_step_deg = start_current_iteration_max_angle_step_deg,
      start_current_iteration_only_for_large_cases = start_current_iteration_only_for_large_cases,
      apslf_start_enabled = outer == 0 ? apslf_start_enabled : false,
      apslf_start_order = apslf_start_order,
      qlimits_enabled = false,
      qlimit_enforcement_mode = :active_set,
      qlimit_lock_reason = :classical_outer_loop,
      wrong_branch_detection = wrong_branch_detection,
      wrong_branch_rescue = wrong_branch_rescue,
      wrong_branch_min_vm_pu = wrong_branch_min_vm_pu,
      wrong_branch_max_vm_pu = wrong_branch_max_vm_pu,
      wrong_branch_max_angle_spread_deg = wrong_branch_max_angle_spread_deg,
      wrong_branch_max_branch_angle_deg = wrong_branch_max_branch_angle_deg,
      wrong_branch_min_low_vm_count = wrong_branch_min_low_vm_count,
      wrong_branch_rescue_max_attempts = wrong_branch_rescue_max_attempts,
      performance_profile = performance_profile,
      rectangular_workspace_reuse = rectangular_workspace_reuse,
      rectangular_preallocate_workspace = rectangular_preallocate_workspace,
      rectangular_workspace_min_buses = rectangular_workspace_min_buses,
    )
    total_iters += iters
    if erg != 0
      final_outcome = outer == 0 ? :base_pf_not_converged : :pf_not_converged_after_qlimit_update
      break
    end
    outer == 0 && (base_pf_converged = true)

    Yred = createYBUS(net = net, sparse = true, printYBUS = false)
    Ybus = (size(Yred, 1) == length(net.nodeVec)) ? Yred : _expand_ybus_for_isolated_nodes(Yred, length(net.nodeVec), net.isoNodes)
    V, slack_idx = initialVrect(net; flatstart = false)
    qmin_pu, qmax_pu = getQLimits_pu(net)
    bus_types = [getNodeType(node) == Slack ? :Slack : getNodeType(node) == PV ? :PV : :PQ for node in net.nodeVec]
    qreq_pu = _rectangular_qgen_requirements_pu(net, Ybus, V)
    violations = _matpower_q_limit_violations(net, qreq_pu, bus_types, fixed_gens; tolerance_pu = max(tol, net.q_hyst_pu))
    if isempty(violations)
      final_outcome = :converged
      break
    end
    qlimit_enforcement_started = true
    selected = mode == :classic_simultaneous ? violations : [violations[argmax(getfield.(violations, :violation_pu))]]
    ref_changed = false
    for v in selected
      bus = v.bus
      old_type = getNodeType(net.nodeVec[bus])
      net.prosumpsVec[v.gen_index].qVal = v.q_limit_pu * net.baseMVA
      fixed_gens[v.gen_index] = true
      if (old_type == Slack || old_type == PV) && !_matpower_bus_has_regulating_generator(net, bus, fixed_gens)
        setNodeType!(net.nodeVec[bus], "PQ")
        logQLimitHit!(net, outer + 1, bus, v.side)
        if old_type == Slack
          new_ref = _select_new_reference_bus!(net, bus)
          if isnothing(new_ref)
            final_outcome = :no_reference_bus_remaining
            push!(outer_rows, (outer_iter = outer + 1, mode = mode, gen_index = v.gen_index, bus_i = bus, violation_side = v.side, qg_before = v.q_before_pu * net.baseMVA, q_limit = v.q_limit_pu * net.baseMVA, violation_mvar = v.violation_pu * net.baseMVA, action = :clamp_and_convert_no_reference, bus_type_before = old_type, bus_type_after = :PQ, ref_changed = true))
            _set_rectangular_pf_status!(net, (; rectangular_pf_status(net)..., qlimit_enforcement_mode = mode, base_pf_converged = base_pf_converged, qlimit_enforcement_started = qlimit_enforcement_started, final_outcome = final_outcome, matpower_outer_iterations = outer + 1, matpower_outer_loop = outer_rows))
            return total_iters, 1
          end
          ref_changed = true
        end
      end
      push!(outer_rows, (outer_iter = outer + 1, mode = mode, gen_index = v.gen_index, bus_i = bus, violation_side = v.side, qg_before = v.q_before_pu * net.baseMVA, q_limit = v.q_limit_pu * net.baseMVA, violation_mvar = v.violation_pu * net.baseMVA, action = :clamp_and_convert, bus_type_before = old_type, bus_type_after = getNodeType(net.nodeVec[bus]), ref_changed = ref_changed))
    end
    verbose > 0 && @printf(stdout, "Classical Q-limit outer iter %d mode=%s selected=%d max_violation=%.6g MVAr ref_changed=%s\n", outer + 1, String(mode), length(selected), maximum(getfield.(violations, :violation_pu)) * net.baseMVA, string(ref_changed))
    outer == qlimit_max_outer && (final_outcome = :max_outer_iterations)
  end
  st = rectangular_pf_status(net)
  if st !== nothing
    _set_rectangular_pf_status!(net, (; st..., qlimit_enforcement_mode = mode, base_pf_converged = base_pf_converged, qlimit_enforcement_started = qlimit_enforcement_started, final_outcome = final_outcome, matpower_outer_iterations = maximum([0; getfield.(outer_rows, :outer_iter)]), matpower_outer_loop = outer_rows, qlimit_reenable_events = 0))
  end
  return total_iters, final_outcome == :converged ? 0 : 1
end


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
# file: src/powerflow_rectangular/current_iteration_start.jl
#

# Rectangular current-iteration start helpers.
#
# This file contains guarded start-value preconditioning logic only.
# Current-iteration is not a solver; rejected candidates must restore the
# original start values and diagnostics must not influence Newton updates.

function _run_guarded_current_iteration_start(Ybus, Vraw::Vector{ComplexF64}, S::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int; enabled::Bool, max_iter::Int, tol::Float64, damping::Float64, accept_only_if_improved::Bool, min_improvement_factor::Float64, vm_min_pu::Float64, vm_max_pu::Float64, max_angle_step_deg::Float64, only_for_large_cases::Bool, large_case_min_buses::Int, verbose::Int = 0, performance_profile = nothing)
  if !enabled
    summary = _current_iteration_summary(enabled = false, reason = :disabled)
    summary = _record_current_iteration_summary!(performance_profile, summary)
    return Vraw, summary
  end
  if only_for_large_cases && length(Vraw) < large_case_min_buses
    summary = _current_iteration_summary(enabled = true, reason = :skipped_small_case)
    summary = _record_current_iteration_summary!(performance_profile, summary)
    return Vraw, summary
  end
  original = copy(Vraw)
  initial_mismatch = _max_rectangular_mismatch(Ybus, original, S, bus_types, Vset, slack_idx)
  V = copy(original)
  best = copy(original)
  best_mismatch = initial_mismatch
  guard_violations = String[]
  max_angle_step = 0.0
  candidate = copy(original)
  candidate_max_angle_step = 0.0
  rejected_at_iteration = 0
  rejection_stage = :other
  iterations = 0
  reason = :max_iter
  for it in 1:max_iter
    iterations = it
    Vnext = copy(V)
    I = Ybus * V
    @inbounds for k in eachindex(V)
      k == slack_idx && continue
      bus_types[k] == :PQ || continue
      if abs(V[k]) <= eps(Float64) || !isfinite(real(V[k])) || !isfinite(imag(V[k]))
        push!(guard_violations, "non_finite_or_zero_voltage_bus_$(k)")
        reason = :invalid_voltage
        rejected_at_iteration = it
        rejection_stage = :nonfinite_check
        break
      end
      yii = Ybus[k, k]
      if abs(yii) <= eps(Float64)
        push!(guard_violations, "zero_diagonal_admittance_bus_$(k)")
        reason = :singular_current_update
        rejected_at_iteration = it
        rejection_stage = :candidate_guard
        break
      end
      target_i = conj(S[k] / V[k])
      candidate = (target_i - (I[k] - yii * V[k])) / yii
      Vnext[k] = (1.0 - damping) * V[k] + damping * candidate
    end
    isempty(guard_violations) || break
    Vnext[slack_idx] = original[slack_idx]
    @inbounds for k in eachindex(Vnext)
      bus_types[k] == :PV || continue
      Vnext[k] = _apply_voltage_magnitude_preserving_angle(Vnext[k], Vset[k])
    end
    if any(v -> !isfinite(real(v)) || !isfinite(imag(v)), Vnext)
      candidate = copy(Vnext)
      push!(guard_violations, "non_finite_voltage")
      reason = :invalid_voltage
      rejected_at_iteration = it
      rejection_stage = :nonfinite_check
      break
    end
    step_deg = maximum(abs.(angle.(Vnext ./ V))) * 180.0 / π
    candidate_max_angle_step = step_deg
    vm = abs.(Vnext)
    if any(v -> v < vm_min_pu || v > vm_max_pu, vm)
      candidate = copy(Vnext)
      push!(guard_violations, "voltage_magnitude_out_of_range")
      reason = :voltage_magnitude_guard
      rejected_at_iteration = it
      rejection_stage = :candidate_guard
      break
    end
    max_angle_step = max(max_angle_step, step_deg)
    if step_deg > max_angle_step_deg
      candidate = copy(Vnext)
      push!(guard_violations, "angle_step_guard")
      reason = :angle_step_guard
      rejected_at_iteration = it
      rejection_stage = :candidate_guard
      break
    end
    mis = _max_rectangular_mismatch(Ybus, Vnext, S, bus_types, Vset, slack_idx)
    if !isfinite(mis)
      candidate = copy(Vnext)
      push!(guard_violations, "non_finite_mismatch")
      reason = :invalid_mismatch
      rejected_at_iteration = it
      rejection_stage = :mismatch_check
      break
    end
    if isfinite(mis) && mis < best_mismatch
      best = copy(Vnext)
      best_mismatch = mis
    end
    candidate = copy(Vnext)
    V = Vnext
    if mis <= tol
      reason = :tolerance_reached
      break
    end
  end
  threshold = accept_only_if_improved ? initial_mismatch * min_improvement_factor : Inf
  accepted = isfinite(best_mismatch) && (accept_only_if_improved ? best_mismatch <= threshold : true)
  accepted || (reason = isempty(guard_violations) ? :not_improved : reason)
  summary0 = _current_iteration_summary(enabled = true, attempted = true, accepted = accepted, iterations = iterations, initial_mismatch = initial_mismatch, final_mismatch = best_mismatch, reason = reason)
  artifact = _write_current_iteration_start_log(performance_profile, summary0; damping, tol, vm_before = abs.(original), vm_after = abs.(accepted ? best : original), vm_candidate = abs.(accepted ? best : candidate), vm_min_pu, vm_max_pu, candidate_max_angle_step, max_angle_step, restored = !accepted, rejected_at_iteration, rejection_stage, guard_violations)
  summary = merge(summary0, (current_iteration_artifact = artifact,))
  summary = _record_current_iteration_summary!(performance_profile, summary)
  verbose > 0 && @info "Current-iteration start pre-solve" accepted = accepted reason = reason initial_mismatch = initial_mismatch final_mismatch = best_mismatch iterations = iterations
  return accepted ? best : original, summary
end

"""
    _run_guarded_apslf_start(Ybus, Vraw, S, bus_types, Vset, slack_idx; enabled, order, baseMVA, verbose=0, performance_profile=nothing) -> (V, summary)

Guarded start-value preconditioner ahead of the rectangular Newton-Raphson
solve (`power_flow.apslf_start`), analogous to
[`_run_guarded_current_iteration_start`](@ref) but sourcing the candidate
profile from the AnalyticLoadFlow.jl-backed analytic power-series solver
instead of a manual current-injection update.

APSLF is used purely as a start-value improver here, not a solver:
`nr_polish` is always disabled (the downstream Newton-Raphson step performs
the polish) and Q-limits are not considered (matching the scope of the
existing current-iteration pre-solve). The candidate is only adopted when it
strictly improves the rectangular mismatch relative to the incoming start
profile; if AnalyticLoadFlow.jl is not loaded, `apslf_solver` raises its
standard "not installed" error immediately (this is a hard configuration
error, not a numerical guard rejection). Numerical failures of the APSLF
solve itself (non-finite result, internal exception) are caught and treated
as "not improved", restoring the incoming start values.
"""
function _run_guarded_apslf_start(Ybus, Vraw::Vector{ComplexF64}, S::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int; enabled::Bool, order::Int, baseMVA::Float64, verbose::Int = 0, performance_profile = nothing)
  if !enabled
    summary = _apslf_start_summary(enabled = false, reason = :disabled)
    summary = _record_apslf_start_summary!(performance_profile, summary)
    return Vraw, summary
  end
  n = length(Vraw)
  initial_mismatch = _max_rectangular_mismatch(Ybus, Vraw, S, bus_types, Vset, slack_idx)
  solver = apslf_solver(order = order, use_pade = true, nr_polish = false)
  accepted = false
  final_mismatch = initial_mismatch
  reason = :not_improved
  candidate = Vraw
  try
    model = PFModel(
      Ybus = Ybus,
      baseMVA = baseMVA,
      busIdx_net = collect(1:n),
      busType = bus_types,
      slack_idx = slack_idx,
      Vset = Vset,
      Sspec = S,
      V0 = Vraw,
    )
    sol = solvePf(solver, model)
    if length(sol.V) == n && all(v -> isfinite(real(v)) && isfinite(imag(v)), sol.V)
      candidate = sol.V
      trial_mismatch = _max_rectangular_mismatch(Ybus, candidate, S, bus_types, Vset, slack_idx)
      if isfinite(trial_mismatch) && trial_mismatch < initial_mismatch
        accepted = true
        final_mismatch = trial_mismatch
        reason = :improved
      end
    else
      reason = :invalid_voltage
    end
  catch err
    reason = :apslf_solve_error
    verbose > 0 && @warn "APSLF start pre-solve failed; keeping incoming start values" exception = (err, catch_backtrace())
  end
  summary = _apslf_start_summary(enabled = true, attempted = true, accepted = accepted, order = order, initial_mismatch = initial_mismatch, final_mismatch = final_mismatch, reason = reason)
  summary = _record_apslf_start_summary!(performance_profile, summary)
  _write_apslf_start_log(performance_profile, summary)
  verbose > 0 && @info "APSLF start pre-solve" accepted = accepted reason = reason initial_mismatch = initial_mismatch final_mismatch = final_mismatch order = order
  return accepted ? candidate : Vraw, summary
end

@inline function _has_vset_adjust_config(ps::ProSumer)::Bool
  return !isnothing(ps.vset_adjust) || !(isnothing(ps.vstep_pu) && isnothing(ps.tap_steps_down) && isnothing(ps.tap_steps_up))
end

@inline function _bus_label(net::Net, bus::Int)::String
  return getCompName(net.nodeVec[bus].comp)
end

@inline function _resolve_vset_adjust_config(ps::ProSumer)::Union{Nothing,VoltageAdjustConfig}
  if !isnothing(ps.vset_adjust)
    return ps.vset_adjust
  end
  # Legacy per-field voltage-adjust settings are only accepted if complete.
  # Partial definitions are intentionally ignored to avoid ambiguous behavior.
  has_any = _has_vset_adjust_config(ps)
  has_any || return nothing
  all_defined = !isnothing(ps.vstep_pu) && !isnothing(ps.tap_steps_down) && !isnothing(ps.tap_steps_up)
  all_defined || return nothing
  return VoltageAdjustConfig(Float64(ps.vstep_pu), Int(ps.tap_steps_down), Int(ps.tap_steps_up))
end

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

# file: src/powerflow_rectangular/rectangular_diagnostics.jl
#
# This file contains diagnostic assembly helpers for rectangular power-flow runs.
# It must not change numerical iteration behavior.

function _current_iteration_summary(; enabled::Bool, attempted::Bool = false, accepted::Bool = false, iterations::Int = 0, initial_mismatch::Float64 = NaN, final_mismatch::Float64 = NaN, reason::Symbol = enabled ? :disabled : :disabled, artifact::String = "")
  return (
    current_iteration_enabled = enabled,
    current_iteration_attempted = attempted,
    current_iteration_accepted = accepted,
    current_iteration_iterations = iterations,
    current_iteration_initial_mismatch = initial_mismatch,
    current_iteration_final_mismatch = final_mismatch,
    current_iteration_reason = reason,
    current_iteration_artifact = artifact,
  )
end

function _record_current_iteration_summary!(performance_profile, summary)
  performance_profile isa AbstractDict || return summary
  previous = get(performance_profile, :current_iteration_start, nothing)
  if previous !== nothing && getproperty(previous, :current_iteration_attempted) === true && summary.current_iteration_attempted !== true
    # Current-iteration diagnostics describe the start-value pre-solve, not the
    # final Newton attempt. Preserve attempted/accepted/rejected metadata across
    # Q-limit outer-loop retries so final logs do not report a real pre-solve as
    # disabled.
    return previous
  end
  performance_profile[:current_iteration_start] = summary
  return summary
end

function _current_iteration_voltage_guard_diagnostics(vm, vm_min_pu::Float64, vm_max_pu::Float64)
  low = findall(<(vm_min_pu), vm)
  high = findall(>(vm_max_pu), vm)
  low_bus = isempty(low) ? "none" : string(low[argmin(vm[low])])
  low_value = isempty(low) ? "none" : string(minimum(vm[low]))
  high_bus = isempty(high) ? "none" : string(high[argmax(vm[high])])
  high_value = isempty(high) ? "none" : string(maximum(vm[high]))
  return (
    low_count = length(low),
    high_count = length(high),
    worst_low_bus = low_bus,
    worst_low_value = low_value,
    worst_high_bus = high_bus,
    worst_high_value = high_value,
  )
end

function _write_current_iteration_start_log(performance_profile, summary; damping::Float64, tol::Float64, vm_before, vm_after, vm_candidate, vm_min_pu::Float64, vm_max_pu::Float64, candidate_max_angle_step::Float64, max_angle_step::Float64, restored::Bool, rejected_at_iteration::Int = 0, rejection_stage::Symbol = :other, guard_violations = String[])
  performance_profile isa AbstractDict || return ""
  output_dir = get(performance_profile, :output_dir, "")
  output_dir isa AbstractString && !isempty(output_dir) || return ""
  mkpath(output_dir)
  path = joinpath(output_dir, "current_iteration_start.log")
  candidate_guard = _current_iteration_voltage_guard_diagnostics(vm_candidate, vm_min_pu, vm_max_pu)
  open(path, "w") do io
    println(io, "current_iteration_enabled: ", summary.current_iteration_enabled)
    println(io, "current_iteration_attempted: ", summary.current_iteration_attempted)
    println(io, "current_iteration_accepted: ", summary.current_iteration_accepted)
    println(io, "current_iteration_reason: ", summary.current_iteration_reason)
    println(io, "initial_mismatch: ", summary.current_iteration_initial_mismatch)
    println(io, "final_mismatch: ", summary.current_iteration_final_mismatch)
    println(io, "iterations: ", summary.current_iteration_iterations)
    println(io, "damping: ", damping)
    println(io, "tol: ", tol)
    println(io, "voltage_magnitude_min_before: ", isempty(vm_before) ? NaN : minimum(vm_before))
    println(io, "voltage_magnitude_max_before: ", isempty(vm_before) ? NaN : maximum(vm_before))
    println(io, "voltage_magnitude_min_after: ", isempty(vm_after) ? NaN : minimum(vm_after))
    println(io, "voltage_magnitude_max_after: ", isempty(vm_after) ? NaN : maximum(vm_after))
    println(io, "candidate_voltage_magnitude_min: ", isempty(vm_candidate) ? NaN : minimum(vm_candidate))
    println(io, "candidate_voltage_magnitude_max: ", isempty(vm_candidate) ? NaN : maximum(vm_candidate))
    println(io, "candidate_voltage_low_count: ", candidate_guard.low_count)
    println(io, "candidate_voltage_high_count: ", candidate_guard.high_count)
    println(io, "candidate_voltage_worst_low_bus: ", candidate_guard.worst_low_bus)
    println(io, "candidate_voltage_worst_low_value: ", candidate_guard.worst_low_value)
    println(io, "candidate_voltage_worst_high_bus: ", candidate_guard.worst_high_bus)
    println(io, "candidate_voltage_worst_high_value: ", candidate_guard.worst_high_value)
    println(io, "candidate_max_angle_step_deg: ", candidate_max_angle_step)
    println(io, "rejected_at_iteration: ", rejected_at_iteration)
    println(io, "rejection_stage: ", rejection_stage)
    println(io, "maximum_angle_step_deg: ", max_angle_step)
    println(io, "guard_violations: ", isempty(guard_violations) ? "none" : join(guard_violations, ", "))
    println(io, "original_start_values_restored: ", restored)
    println(io, "restored_voltage_magnitude_min: ", restored && !isempty(vm_after) ? minimum(vm_after) : "none")
    println(io, "restored_voltage_magnitude_max: ", restored && !isempty(vm_after) ? maximum(vm_after) : "none")
  end
  return path
end

function _apslf_start_summary(; enabled::Bool, attempted::Bool = false, accepted::Bool = false, order::Int = 0, initial_mismatch::Float64 = NaN, final_mismatch::Float64 = NaN, reason::Symbol = :disabled)
  return (
    apslf_start_enabled = enabled,
    apslf_start_attempted = attempted,
    apslf_start_accepted = accepted,
    apslf_start_order = order,
    apslf_start_initial_mismatch = initial_mismatch,
    apslf_start_final_mismatch = final_mismatch,
    apslf_start_reason = reason,
  )
end

function _record_apslf_start_summary!(performance_profile, summary)
  performance_profile isa AbstractDict || return summary
  previous = get(performance_profile, :apslf_start, nothing)
  if previous !== nothing && getproperty(previous, :apslf_start_attempted) === true && summary.apslf_start_attempted !== true
    # Preserve attempted/accepted metadata across Q-limit outer-loop retries
    # (only the first outer pass runs the APSLF start pre-solve).
    return previous
  end
  performance_profile[:apslf_start] = summary
  return summary
end

function _write_apslf_start_log(performance_profile, summary)
  performance_profile isa AbstractDict || return ""
  output_dir = get(performance_profile, :output_dir, "")
  output_dir isa AbstractString && !isempty(output_dir) || return ""
  mkpath(output_dir)
  path = joinpath(output_dir, "apslf_start.log")
  open(path, "w") do io
    println(io, "apslf_start_enabled: ", summary.apslf_start_enabled)
    println(io, "apslf_start_attempted: ", summary.apslf_start_attempted)
    println(io, "apslf_start_accepted: ", summary.apslf_start_accepted)
    println(io, "apslf_start_order: ", summary.apslf_start_order)
    println(io, "apslf_start_reason: ", summary.apslf_start_reason)
    println(io, "initial_mismatch: ", summary.apslf_start_initial_mismatch)
    println(io, "final_mismatch: ", summary.apslf_start_final_mismatch)
  end
  return path
end

function _rectangular_mismatch_rows(Ybus, V::Vector{ComplexF64}, S::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int; net = nothing, top_n::Int = 10)
  final_F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
  rows = NamedTuple[]
  row = 1
  @inbounds for bus in eachindex(V)
    bus == slack_idx && continue
    bus_id = net === nothing ? bus : _qlimit_original_bus_id(net, bus)
    push!(rows, (row = row, bus_index = bus, bus_id = bus_id, equation = :P, mismatch = final_F[row], abs_mismatch = abs(final_F[row])))
    equation = bus_types[bus] == :PQ ? :Q : :voltage_setpoint
    push!(rows, (row = row + 1, bus_index = bus, bus_id = bus_id, equation = equation, mismatch = final_F[row + 1], abs_mismatch = abs(final_F[row + 1])))
    row += 2
  end
  sort!(rows; by = r -> -r.abs_mismatch)
  top = top_n < 0 ? rows : collect(Iterators.take(rows, min(top_n, length(rows))))
  worst = isempty(rows) ? (row = 0, bus_index = 0, bus_id = 0, equation = :none, mismatch = NaN, abs_mismatch = NaN) : first(rows)
  return (rows = rows, top = top, worst = worst)
end

function _rectangular_mismatch_trend(history::AbstractVector{<:Real}; window::Int = 10)
  isempty(history) && return :unavailable
  finite_history = filter(isfinite, history)
  if length(finite_history) < length(history)
    isempty(finite_history) && return :nonfinite
    return last(finite_history) > first(finite_history) ? :diverging_to_nonfinite : :nonfinite_after_improvement
  end
  tail = collect(Iterators.drop(history, max(length(history) - window, 0)))
  length(tail) < 3 && return :insufficient_history
  diffs = diff(tail)
  tol = max(1e-12, 1e-6 * maximum(abs.(tail)))
  all(d -> d <= tol, diffs) && return :monotonic
  all(abs(d) <= tol for d in diffs) && return :stagnant
  signs = [d > tol ? 1 : d < -tol ? -1 : 0 for d in diffs]
  nonzero = filter(!=(0), signs)
  length(nonzero) >= 2 && any(nonzero[i] != nonzero[i - 1] for i in 2:length(nonzero)) && return :oscillatory
  abs(tail[end] - minimum(tail)) <= tol && return :stagnant
  return :mixed
end

function _rectangular_mismatch_snapshot(label::Symbol, iteration::Int, Ybus, V::Union{Nothing,Vector{ComplexF64}}, S::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int; net = nothing, top_n::Int = 10)
  V === nothing && return (label = label, iteration = iteration, rows = NamedTuple[])
  rows = _rectangular_mismatch_rows(Ybus, V, S, bus_types, Vset, slack_idx; net, top_n).top
  return (label = label, iteration = iteration, rows = rows)
end

function _rectangular_step_statistics(step_diagnostics)
  step_diagnostics isa AbstractVector || return (
    autodamp_step_count = 0,
    autodamp_min_alpha = NaN,
    autodamp_max_alpha = NaN,
    autodamp_mean_alpha = NaN,
    autodamp_floor_hits = 0,
    autodamp_nonimproving_steps = 0,
    autodamp_failure = false,
  )
  isempty(step_diagnostics) && return (
    autodamp_step_count = 0,
    autodamp_min_alpha = NaN,
    autodamp_max_alpha = NaN,
    autodamp_mean_alpha = NaN,
    autodamp_floor_hits = 0,
    autodamp_nonimproving_steps = 0,
    autodamp_failure = false,
  )
  alphas = Float64[getproperty(row, :alpha) for row in step_diagnostics]
  min_alpha = minimum(alphas)
  return (
    autodamp_step_count = length(alphas),
    autodamp_min_alpha = min_alpha,
    autodamp_max_alpha = maximum(alphas),
    autodamp_mean_alpha = sum(alphas) / length(alphas),
    autodamp_floor_hits = count(a -> isapprox(a, min_alpha; atol = 1e-12, rtol = 1e-12), alphas),
    autodamp_nonimproving_steps = count(row -> !getproperty(row, :accepted_improvement), step_diagnostics),
    autodamp_failure = count(row -> !getproperty(row, :accepted_improvement), step_diagnostics) >= max(5, length(step_diagnostics) ÷ 2) && count(a -> isapprox(a, min_alpha; atol = 1e-12, rtol = 1e-12), alphas) >= max(5, length(step_diagnostics) ÷ 2),
  )
end

function _rectangular_mismatch_diagnostics(Ybus, V::Vector{ComplexF64}, S::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int, final_pv_voltage_residual::Float64; net = nothing, history = Float64[], step_diagnostics = nothing, top_n::Int = 10, best_finite_iteration::Int = 0, best_finite_voltage = nothing, last_finite_iteration::Int = 0, last_finite_voltage = nothing, first_nonfinite_iteration::Int = 0, first_nonfinite_voltage = nothing)
  final_F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
  max_active_power_mismatch = 0.0
  max_reactive_power_mismatch = 0.0
  max_voltage_residual_or_setpoint_residual = final_pv_voltage_residual
  row = 1
  @inbounds for bus in eachindex(V)
    bus == slack_idx && continue
    max_active_power_mismatch = max(max_active_power_mismatch, abs(final_F[row]))
    if bus_types[bus] == :PQ
      max_reactive_power_mismatch = max(max_reactive_power_mismatch, abs(final_F[row + 1]))
    else
      max_voltage_residual_or_setpoint_residual = max(max_voltage_residual_or_setpoint_residual, abs(final_F[row + 1]))
    end
    row += 2
  end
  mismatch_rows = _rectangular_mismatch_rows(Ybus, V, S, bus_types, Vset, slack_idx; net, top_n)
  step_stats = _rectangular_step_statistics(step_diagnostics)
  finite_history = filter(isfinite, history)
  primary_mismatch_rows = !isempty(finite_history) && !all(isfinite, final_F) ?
    _rectangular_mismatch_snapshot(:last_finite_iteration, last_finite_iteration, Ybus, last_finite_voltage, S, bus_types, Vset, slack_idx; net, top_n).rows :
    mismatch_rows.top
  top_mismatch_snapshots = [
    _rectangular_mismatch_snapshot(:best_finite_iteration, best_finite_iteration, Ybus, best_finite_voltage, S, bus_types, Vset, slack_idx; net, top_n),
    _rectangular_mismatch_snapshot(:last_finite_iteration, last_finite_iteration, Ybus, last_finite_voltage, S, bus_types, Vset, slack_idx; net, top_n),
  ]
  first_nonfinite_iteration > 0 && push!(top_mismatch_snapshots, _rectangular_mismatch_snapshot(:first_nonfinite_iteration, first_nonfinite_iteration, Ybus, first_nonfinite_voltage, S, bus_types, Vset, slack_idx; net, top_n))
  return (
    max_active_power_mismatch = max_active_power_mismatch,
    max_reactive_power_mismatch = max_reactive_power_mismatch,
    max_voltage_residual_or_setpoint_residual_where_available = max_voltage_residual_or_setpoint_residual,
    worst_mismatch_bus_id = mismatch_rows.worst.bus_id,
    worst_mismatch_bus_index = mismatch_rows.worst.bus_index,
    worst_mismatch_equation = mismatch_rows.worst.equation,
    worst_mismatch_value = mismatch_rows.worst.mismatch,
    top_mismatch_rows = primary_mismatch_rows,
    top_mismatch_rows_by_iteration = top_mismatch_snapshots,
    mismatch_history = collect(history),
    mismatch_history_first = isempty(history) ? NaN : first(history),
    mismatch_history_last = isempty(history) ? NaN : last(history),
    mismatch_history_best = isempty(finite_history) ? NaN : minimum(finite_history),
    best_mismatch = isempty(finite_history) ? NaN : minimum(finite_history),
    mismatch_history_trend = _rectangular_mismatch_trend(history),
    first_nonfinite_iteration = first_nonfinite_iteration > 0 ? first_nonfinite_iteration : nothing,
    last_finite_mismatch = isempty(finite_history) ? NaN : finite_history[end],
    step_stats...,
  )
end

function _merge_current_iteration_diagnostics(status_build, performance_profile)
  if performance_profile isa AbstractDict && haskey(performance_profile, :current_iteration_start)
    ci = performance_profile[:current_iteration_start]
    # Current-iteration diagnostics describe the start-value pre-solve, not the
    # final Newton attempt. Preserve them across Q-limit outer-loop retries so
    # final logs do not report the pre-solve as disabled after it was attempted.
    return merge(status_build, (status = (; status_build.status..., ci...),))
  end
  return status_build
end

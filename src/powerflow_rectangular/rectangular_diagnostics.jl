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

# Independently solved AC islands share one performance_profile/output_dir; the
# solver driver sets :diagnostic_artifact_prefix (e.g. "ac_island_2_") for the
# duration of each island's solve so fixed-name diagnostic artifacts written
# here (current_iteration_start.log, apslf_start.log, merit_linesearch.log,
# trust_region.log) do not silently overwrite each other. Empty/absent outside
# an island split, so filenames are unchanged for the common single-network case.
_diagnostic_artifact_prefix(performance_profile)::String = performance_profile isa AbstractDict ? String(get(performance_profile, :diagnostic_artifact_prefix, "")) : ""

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
  path = joinpath(output_dir, string(_diagnostic_artifact_prefix(performance_profile), "current_iteration_start.log"))
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
  path = joinpath(output_dir, string(_diagnostic_artifact_prefix(performance_profile), "apslf_start.log"))
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

"""
    _branch_anomaly_diagnostics(net, bus_index; top_n=5) -> Vector{NamedTuple}

Scan the in-service branches incident to `bus_index` (an internal bus index,
matching `net.nodeVec`/`V` indexing) for modeling traits that commonly explain
a persistently poor mismatch at that bus: zero impedance, an off-nominal
transformer tap ratio, a large phase shift, or an unusually resistive
(low X/R) branch. Each incident branch is scored by how many traits it
triggers; rows are returned worst-first, capped at `top_n` (`top_n < 0`
returns all incident branches).

This is a passive read of already-imported network data — it does not modify
`net` or influence solver behavior. It exists to turn "bus 42 has the worst
mismatch" into "bus 42 has the worst mismatch, and its transformer branch #17
has a suspicious 0.62 tap ratio" for `diagnose.log`.
"""
_median(xs::AbstractVector{<:Real})::Float64 = (s = sort(xs); n = length(s); isodd(n) ? Float64(s[(n+1)÷2]) : Float64(s[n÷2] + s[n÷2+1]) / 2)

function _branch_anomaly_diagnostics(net::Net, bus_index::Integer; top_n::Int = 5)
  rows = NamedTuple[]
  bus_index < 1 && return rows
  incident = filter(br -> br.status == 1 && (br.fromBus == bus_index || br.toBus == bus_index), net.branchVec)
  isempty(incident) && return rows
  @inbounds for br in incident
    is_transformer = br.ratio != 0.0 || br.has_ratio_tap || br.has_phase_tap
    flags = String[]
    x_over_r = br.r_pu == 0.0 ? Inf : abs(br.x_pu / br.r_pu)
    if br.r_pu == 0.0 && br.x_pu == 0.0
      push!(flags, "zero_impedance")
    end
    if is_transformer && br.ratio != 0.0 && !(0.8 <= br.ratio <= 1.2)
      push!(flags, "off_nominal_tap_ratio")
    end
    if abs(br.phase_shift_deg) > 30.0
      push!(flags, "large_phase_shift")
    end
    if br.r_pu != 0.0 && x_over_r < 0.1
      push!(flags, "unusually_resistive_low_xr")
    end
    # Absolute thresholds above catch textbook-unusual impedances, but a
    # branch with plausible-looking R/X in isolation can still be a modeling
    # error (e.g. a shift/scale mistake) relative to its own neighbors at the
    # same bus. Comparing against sibling branches of the same kind (line vs
    # transformer) catches that without hardcoding a "normal" impedance range.
    # Real transmission lines at one bus routinely span a 3-5x reactance
    # range from length alone (confirmed empirically on case14's own,
    # correctly-imported data), so the outlier ratio below is set well above
    # that natural spread — it is meant to catch import/unit-convention
    # mistakes (typically one to two orders of magnitude), not ordinary
    # engineering variance. Requires >=3 same-kind siblings for a meaningful
    # comparison.
    siblings = filter(other -> other.branchIdx != br.branchIdx && (other.ratio != 0.0 || other.has_ratio_tap || other.has_phase_tap) == is_transformer && other.x_pu > 0.0, incident)
    if length(siblings) >= 3 && br.x_pu > 0.0
      sibling_median_x = _median([other.x_pu for other in siblings])
      if sibling_median_x > 0.0
        ratio_to_siblings = br.x_pu / sibling_median_x
        (ratio_to_siblings > 8.0 || ratio_to_siblings < 1 / 8) && push!(flags, "reactance_outlier_vs_sibling_branches")
      end
    end
    branch_name = hasproperty(br.comp, :cName) ? String(br.comp.cName) : ""
    push!(rows, (
      branch_index = br.branchIdx,
      branch_name = branch_name,
      from_bus = br.fromBus,
      to_bus = br.toBus,
      is_transformer = is_transformer,
      r_pu = br.r_pu,
      x_pu = br.x_pu,
      x_over_r = x_over_r,
      ratio = br.ratio,
      phase_shift_deg = br.phase_shift_deg,
      flags = flags,
    ))
  end
  sort!(rows; by = r -> -length(r.flags))
  top_n < 0 && return rows
  return collect(Iterators.take(rows, min(top_n, length(rows))))
end

function _format_branch_anomaly_rows(rows)::String
  isempty(rows) && return "none"
  lines = String[]
  for row in rows
    label = isempty(row.branch_name) ? "branch #$(row.branch_index)" : "$(row.branch_name) (branch #$(row.branch_index))"
    descriptor = row.is_transformer ? "transformer" : "line"
    detail = "$(label) [$(descriptor), $(row.from_bus)->$(row.to_bus)]: r_pu=$(row.r_pu), x_pu=$(row.x_pu), x/r=$(round(row.x_over_r; digits=3)), ratio=$(row.ratio), phase_shift_deg=$(row.phase_shift_deg)"
    isempty(row.flags) || (detail *= " flags=[$(join(row.flags, ", "))]")
    push!(lines, detail)
  end
  return join(lines, " | ")
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

# Aggregate summary of the merit-function line search across all Newton iterations.
# `merit_step_diagnostics` is one entry per outer NR iteration (or `nothing`/empty
# when `merit_enabled = false`), pushed by `choose_rectangular_autodamp`.
function _merit_linesearch_summary(merit_step_diagnostics, merit_enabled::Bool)
  if !merit_enabled || merit_step_diagnostics === nothing || isempty(merit_step_diagnostics)
    return (
      merit_enabled = merit_enabled,
      merit_used_iterations = 0,
      merit_fallback_count = 0,
      merit_active_set_skip_count = 0,
      merit_initial = NaN,
      merit_final = NaN,
    )
  end
  used = count(row -> getproperty(row, :accept_reason) === :armijo, merit_step_diagnostics)
  fallback = count(row -> getproperty(row, :accept_reason) in (:fallback_max_mismatch, :fallback_conservative), merit_step_diagnostics)
  skipped = count(row -> getproperty(row, :accept_reason) === :active_set_skip, merit_step_diagnostics)
  first_entry = first(merit_step_diagnostics)
  last_entry = last(merit_step_diagnostics)
  merit_final = isfinite(last_entry.f_after) ? last_entry.f_after : last_entry.f_before
  return (
    merit_enabled = merit_enabled,
    merit_used_iterations = used,
    merit_fallback_count = fallback,
    merit_active_set_skip_count = skipped,
    merit_initial = first_entry.f_before,
    merit_final = merit_final,
  )
end

function _write_merit_linesearch_log(performance_profile, merit_step_diagnostics, summary)
  performance_profile isa AbstractDict || return ""
  merit_step_diagnostics === nothing && return ""
  output_dir = get(performance_profile, :output_dir, "")
  output_dir isa AbstractString && !isempty(output_dir) || return ""
  mkpath(output_dir)
  path = joinpath(output_dir, string(_diagnostic_artifact_prefix(performance_profile), "merit_linesearch.log"))
  open(path, "w") do io
    println(io, "merit_enabled: ", summary.merit_enabled)
    println(io, "merit_used_iterations: ", summary.merit_used_iterations)
    println(io, "merit_fallback_count: ", summary.merit_fallback_count)
    println(io, "merit_active_set_skip_count: ", summary.merit_active_set_skip_count)
    println(io, "merit_initial: ", summary.merit_initial)
    println(io, "merit_final: ", summary.merit_final)
    println(io, "")
    println(io, "# iter f_before directional_derivative f_after tested_alphas accepted_alpha accept_reason")
    for (iter, row) in enumerate(merit_step_diagnostics)
      tested = isempty(row.tested_alphas) ? "none" : join(row.tested_alphas, ",")
      println(io, "iter=", iter, " f_before=", row.f_before, " directional_derivative=", row.directional_derivative, " f_after=", row.f_after, " tested_alphas=[", tested, "] accepted_alpha=", row.accepted_alpha, " accept_reason=", row.accept_reason)
    end
  end
  return path
end

function _merge_merit_linesearch_diagnostics(status_build, performance_profile, merit_step_diagnostics, merit_enabled::Bool)
  summary = _merit_linesearch_summary(merit_step_diagnostics, merit_enabled)
  artifact = _write_merit_linesearch_log(performance_profile, merit_step_diagnostics, summary)
  summary = merge(summary, (merit_linesearch_artifact = artifact,))
  return merge(status_build, (status = (; status_build.status..., summary...),))
end

# Aggregate summary of the trust-region step control across all Newton iterations.
# `tr_step_diagnostics` is one entry per outer NR iteration (or `nothing`/empty when
# `trust_region_enabled = false`), pushed by `choose_rectangular_trust_region_step`.
# tr_min_radius/tr_max_radius report the *observed* radius trajectory (the value in
# effect at the start of each iteration), analogous to autodamp_min_alpha/max_alpha.
function _trust_region_summary(tr_step_diagnostics, trust_region_enabled::Bool)
  if !trust_region_enabled || tr_step_diagnostics === nothing || isempty(tr_step_diagnostics)
    return (
      trust_region_enabled = trust_region_enabled,
      tr_step_count = 0,
      tr_rejected_steps = 0,
      tr_min_radius = NaN,
      tr_max_radius = NaN,
      tr_final_radius = NaN,
      tr_collapsed = false,
      tr_dogleg_newton_count = 0,
      tr_dogleg_interp_count = 0,
      tr_dogleg_cauchy_count = 0,
      tr_active_set_skip_count = 0,
    )
  end
  radii = Float64[getproperty(row, :radius_before) for row in tr_step_diagnostics]
  rejected = sum(getproperty(row, :rejected_steps) for row in tr_step_diagnostics)
  collapsed = any(getproperty(row, :collapsed) for row in tr_step_diagnostics)
  final_entry = last(tr_step_diagnostics)
  dogleg_newton_count = count(row -> getproperty(row, :accept_reason) === :dogleg_newton, tr_step_diagnostics)
  dogleg_interp_count = count(row -> getproperty(row, :accept_reason) === :dogleg_interp, tr_step_diagnostics)
  dogleg_cauchy_count = count(row -> getproperty(row, :accept_reason) === :dogleg_cauchy, tr_step_diagnostics)
  active_set_skip_count = count(row -> getproperty(row, :accept_reason) === :active_set_skip, tr_step_diagnostics)
  return (
    trust_region_enabled = trust_region_enabled,
    tr_step_count = length(tr_step_diagnostics),
    tr_rejected_steps = rejected,
    tr_min_radius = minimum(radii),
    tr_max_radius = maximum(radii),
    tr_final_radius = final_entry.radius_after,
    tr_collapsed = collapsed,
    tr_dogleg_newton_count = dogleg_newton_count,
    tr_dogleg_interp_count = dogleg_interp_count,
    tr_dogleg_cauchy_count = dogleg_cauchy_count,
    tr_active_set_skip_count = active_set_skip_count,
  )
end

function _write_trust_region_log(performance_profile, tr_step_diagnostics, summary)
  performance_profile isa AbstractDict || return ""
  tr_step_diagnostics === nothing && return ""
  output_dir = get(performance_profile, :output_dir, "")
  output_dir isa AbstractString && !isempty(output_dir) || return ""
  mkpath(output_dir)
  path = joinpath(output_dir, string(_diagnostic_artifact_prefix(performance_profile), "trust_region.log"))
  open(path, "w") do io
    println(io, "trust_region_enabled: ", summary.trust_region_enabled)
    println(io, "tr_step_count: ", summary.tr_step_count)
    println(io, "tr_rejected_steps: ", summary.tr_rejected_steps)
    println(io, "tr_min_radius: ", summary.tr_min_radius)
    println(io, "tr_max_radius: ", summary.tr_max_radius)
    println(io, "tr_final_radius: ", summary.tr_final_radius)
    println(io, "tr_collapsed: ", summary.tr_collapsed)
    println(io, "tr_dogleg_newton_count: ", summary.tr_dogleg_newton_count)
    println(io, "tr_dogleg_interp_count: ", summary.tr_dogleg_interp_count)
    println(io, "tr_dogleg_cauchy_count: ", summary.tr_dogleg_cauchy_count)
    println(io, "tr_active_set_skip_count: ", summary.tr_active_set_skip_count)
    println(io, "")
    println(io, "# iter radius_before rho tested_radii rejected_steps accepted radius_after collapsed accept_reason")
    for (iter, row) in enumerate(tr_step_diagnostics)
      tested = isempty(row.tested_radii) ? "none" : join(row.tested_radii, ",")
      println(io, "iter=", iter, " radius_before=", row.radius_before, " rho=", row.rho, " tested_radii=[", tested, "] rejected_steps=", row.rejected_steps, " accepted=", row.accepted, " radius_after=", row.radius_after, " collapsed=", row.collapsed, " accept_reason=", row.accept_reason)
    end
  end
  return path
end

function _merge_trust_region_diagnostics(status_build, performance_profile, tr_step_diagnostics, trust_region_enabled::Bool)
  summary = _trust_region_summary(tr_step_diagnostics, trust_region_enabled)
  artifact = _write_trust_region_log(performance_profile, tr_step_diagnostics, summary)
  summary = merge(summary, (trust_region_artifact = artifact,))
  return merge(status_build, (status = (; status_build.status..., summary...),))
end

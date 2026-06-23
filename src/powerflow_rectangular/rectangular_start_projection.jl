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
# Rectangular power-flow start projection and DC-start helpers.
#
# This file is included inside module Sparlectra. Do not add a module wrapper here.

# Date: 29.5.2026
# file: src/powerflow_rectangular/rectangular_start_projection.jl

function _sanitize_rectangular_start(V::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int)
  # Sanitize NaN/Inf/zero seeds before any projection candidate is evaluated.
  # Reason: downstream mismatch/Jacobian code assumes finite voltages with vm > 0.
  Vs = copy(V)
  @inbounds for k in eachindex(Vs)
    Vk = Vs[k]
    vm = abs(Vk)
    if !isfinite(real(Vk)) || !isfinite(imag(Vk)) || vm <= 0.0
      # For regulating buses keep the configured magnitude if available.
      # Fallback to 1.0 pu keeps the candidate physically plausible.
      vm = (bus_types[k] in (:Slack, :PV) && isfinite(Vset[k]) && Vset[k] > 0.0) ? Vset[k] : 1.0
      # Reset angle to 0 for repaired entries; later stages may set better angles.
      Vs[k] = ComplexF64(vm, 0.0)
    end
  end
  # Preserve the original slack complex value exactly.
  # This avoids accidentally changing the global phase reference in sanitization.
  Vs[slack_idx] = V[slack_idx]
  return Vs
end

function _voltage_magnitude_for_projection(Vraw::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, k::Int)
  vm_raw = abs(Vraw[k])
  if bus_types[k] in (:Slack, :PV) && isfinite(Vset[k]) && Vset[k] > 0.0
    # Slack/PV are magnitude-regulated buses in this initialization context.
    return Vset[k]
  elseif isfinite(vm_raw) && vm_raw > 0.0
    # PQ buses keep their seed magnitude when it is numerically valid.
    return vm_raw
  else
    # Last-resort neutral magnitude.
    return 1.0
  end
end

function _dc_angle_start_rectangular(Ybus, Vraw::Vector{ComplexF64}, S::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int; dc_angle_limit_deg::Float64 = 60.0)
  # DC-angle start builds phase guesses from active-power balance on reduced B.
  # Only non-slack buses are solved; slack is the angular reference.
  n = length(Vraw)
  non_slack = non_slack_indices(n, slack_idx)
  nred = length(non_slack)
  nred == 0 && return copy(Vraw)
  pos = build_pos_map(non_slack, n)

  B = if Ybus isa SparseMatrixCSC
    I = Int[]
    J = Int[]
    V = Float64[]
    sizehint!(I, 2 * nnz(Ybus))
    sizehint!(J, 2 * nnz(Ybus))
    sizehint!(V, 2 * nnz(Ybus))
    rv = rowvals(Ybus)
    nz = nonzeros(Ybus)
    @inbounds for j in axes(Ybus, 2)
      for ptr in nzrange(Ybus, j)
        i = rv[ptr]
        i == j && continue
        ri = pos[i]
        ri == 0 && continue
        bij = imag(nz[ptr])
        bij == 0.0 && continue
        if j != slack_idx
          cj = pos[j]
          if cj != 0
            # Off-diagonal reduced B term.
            push!(I, ri)
            push!(J, cj)
            push!(V, -bij)
          end
        end
        # Diagonal accumulation for row i in reduced coordinates.
        push!(I, ri)
        push!(J, ri)
        push!(V, bij)
      end
    end
    sparse(I, J, V, nred, nred)
  else
    B = zeros(Float64, nred, nred)
    @inbounds for i in non_slack
      ri = pos[i]
      for j in axes(Ybus, 2)
        j == i && continue
        bij = imag(Ybus[i, j])
        bij == 0.0 && continue
        if j != slack_idx
          cj = pos[j]
          if cj != 0
            B[ri, cj] -= bij
          end
        end
        B[ri, ri] += bij
      end
    end
    B
  end

  P = real.(S[non_slack])
  θred = _solve_dc_angle_system(B, P, nothing, nred)
  # Clamp relative angles to avoid extreme starts that can destabilize first NR steps.
  limit = deg2rad(dc_angle_limit_deg)
  θslack = angle(Vraw[slack_idx])
  Vdc = similar(Vraw)

  @inbounds for k = 1:n
    vm = _voltage_magnitude_for_projection(Vraw, bus_types, Vset, k)
    if k == slack_idx
      Vdc[k] = Vraw[k]
    else
      θ = θslack + clamp(θred[pos[k]], -limit, limit)
      Vdc[k] = ComplexF64(vm * cos(θ), vm * sin(θ))
    end
  end
  return Vdc
end

function _record_dc_solve_diagnostics!(performance_profile, B, backend::Symbol, condition_warning, nred::Int)
  performance_profile isa AbstractDict || return nothing
  Bool(get(performance_profile, :enabled, false)) || return nothing
  # Store lightweight matrix/solver diagnostics once per DC solve.
  rows, cols = size(B)
  matrix_nnz = B isa SparseMatrixCSC ? nnz(B) : count(!iszero, B)
  density = rows == 0 || cols == 0 ? 0.0 : matrix_nnz / (rows * cols)
  performance_profile[:dc_matrix_size] = (rows, cols)
  performance_profile[:dc_matrix_nnz] = matrix_nnz
  performance_profile[:dc_matrix_density] = density
  performance_profile[:dc_matrix_is_sparse] = (B isa SparseMatrixCSC)
  performance_profile[:dc_solve_backend] = backend
  performance_profile[:dc_solve_reduced_dimension] = nred
  if condition_warning !== nothing
    performance_profile[:dc_solve_condition_warning] = condition_warning
  else
    # Clear stale warning state from previous runs/attempts.
    delete!(performance_profile, :dc_solve_condition_warning)
  end
  return nothing
end

function _sparse_dc_solve(B::SparseMatrixCSC{Float64}, P::Vector{Float64}, performance_profile, nred::Int)
  backend = :sparse_lu_umfpack
  condition_warning = nothing
  try
    # Primary path: sparse LU is usually fastest for well-conditioned B.
    F = lu(B)
    θred = F \ P
    _record_dc_solve_diagnostics!(performance_profile, B, backend, condition_warning, nred)
    return θred
  catch e
    if _is_rectangular_linear_step_failure(e)
      backend = :sparse_qr_fallback
      condition_warning = :singular_or_ill_conditioned_lu
      try
        # First fallback: sparse QR is more robust near singularity.
        θred = qr(B) \ P
        _record_dc_solve_diagnostics!(performance_profile, B, backend, condition_warning, nred)
        return θred
      catch qr_error
        # Second fallback is intentionally restricted to small systems.
        _dense_svd_fallback_allowed(B) || rethrow(qr_error)
        backend = :dense_svd_fallback_small_system
        θred = _svd_pinv_solve(Matrix(B), P)
        _record_dc_solve_diagnostics!(performance_profile, B, backend, condition_warning, nred)
        return θred
      end
    end
    rethrow(e)
  end
end

function _solve_dc_angle_system(B, P::Vector{Float64}, performance_profile, nred::Int)
  # Sparse and dense paths are split to keep backend diagnostics explicit.
  if B isa SparseMatrixCSC{Float64}
    return _sparse_dc_solve(B, P, performance_profile, nred)
  end
  backend = :dense_backslash
  condition_warning = nothing
  try
    θred = B \ P
    # Condition estimate is gated by matrix size to avoid excessive overhead.
    if max(size(B)...) <= 2000
      c = cond(B)
      if !isfinite(c) || c > 1 / sqrt(eps(Float64))
        condition_warning = :large_condition_estimate
      end
    end
    _record_dc_solve_diagnostics!(performance_profile, B, backend, condition_warning, nred)
    return θred
  catch e
    if _is_rectangular_linear_step_failure(e)
      backend = :dense_qr_or_svd_fallback
      condition_warning = :singular_or_ill_conditioned_backslash
      θred = solve_linear(B, P; allow_pinv = true)
      _record_dc_solve_diagnostics!(performance_profile, B, backend, condition_warning, nred)
      return θred
    end
    rethrow(e)
  end
end

function _dc_angle_start_rectangular_profiled(Ybus, Vraw::Vector{ComplexF64}, S::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int; dc_angle_limit_deg::Float64 = 60.0, performance_profile = nothing)
  n = length(Vraw)
  non_slack = _perf_profile_time!(performance_profile, :start_projection_bus_map) do
    non_slack_indices(n, slack_idx)
  end
  nred = length(non_slack)
  nred == 0 && return copy(Vraw)
  pos = _perf_profile_time!(performance_profile, :start_projection_bus_map) do
    build_pos_map(non_slack, n)
  end

  # Profiled variant keeps assembly/solve phases separated for timing analysis.
  B = _perf_profile_time!(performance_profile, :start_projection_dc_matrix_assembly) do
    if Ybus isa SparseMatrixCSC
      I = Int[]
      J = Int[]
      V = Float64[]
      sizehint!(I, 2 * nnz(Ybus))
      sizehint!(J, 2 * nnz(Ybus))
      sizehint!(V, 2 * nnz(Ybus))
      rv = rowvals(Ybus)
      nz = nonzeros(Ybus)
      @inbounds for j in axes(Ybus, 2)
        for ptr in nzrange(Ybus, j)
          i = rv[ptr]
          i == j && continue
          ri = pos[i]
          ri == 0 && continue
          bij = imag(nz[ptr])
          bij == 0.0 && continue
          if j != slack_idx
            cj = pos[j]
            if cj != 0
              push!(I, ri)
              push!(J, cj)
              push!(V, -bij)
            end
          end
          push!(I, ri)
          push!(J, ri)
          push!(V, bij)
        end
      end
      sparse(I, J, V, nred, nred)
    else
      B = zeros(Float64, nred, nred)
      @inbounds for i in non_slack
        ri = pos[i]
        for j in axes(Ybus, 2)
          j == i && continue
          bij = imag(Ybus[i, j])
          bij == 0.0 && continue
          if j != slack_idx
            cj = pos[j]
            if cj != 0
              B[ri, cj] -= bij
            end
          end
          B[ri, ri] += bij
        end
      end
      B
    end
  end

  P = real.(S[non_slack])
  θred = _perf_profile_time!(performance_profile, :start_projection_dc_linear_solve) do
    _solve_dc_angle_system(B, P, performance_profile, nred)
  end
  limit = deg2rad(dc_angle_limit_deg)
  θslack = angle(Vraw[slack_idx])
  Vdc = similar(Vraw)

  _perf_profile_time!(performance_profile, :start_projection_voltage_clipping) do
    @inbounds for k in eachindex(Vraw)
      vm = _voltage_magnitude_for_projection(Vraw, bus_types, Vset, k)
      if k == slack_idx
        Vdc[k] = Vraw[k]
      else
        θ = θslack + clamp(θred[pos[k]], -limit, limit)
        Vdc[k] = ComplexF64(vm * cos(θ), vm * sin(θ))
      end
    end
  end
  return Vdc
end

function _dc_start_quality_diagnostics(Ybus, Vdc::Union{Nothing,Vector{ComplexF64}}, raw_mis::Float64, dc_mis::Float64, slack_idx::Int, dc_angle_limit_deg::Float64)
  if Vdc === nothing || isempty(Vdc)
    return (
      dc_angle_min_deg = missing, dc_angle_max_deg = missing, dc_angle_spread_deg = missing,
      dc_angle_mean_deg = missing, dc_angle_std_deg = missing, dc_angle_clipped_count = 0,
      dc_angle_clip_limit_deg = dc_angle_limit_deg, dc_max_branch_angle_deg = missing,
      dc_branch_angle_violation_count = 0, worst_dc_branch_from_bus = missing,
      worst_dc_branch_to_bus = missing, worst_dc_branch_angle_deg = missing,
      dc_voltage_magnitude_min = missing, dc_voltage_magnitude_max = missing,
      dc_mismatch_ratio_vs_raw = missing, dc_mismatch_growth_factor = missing,
      requested_dc_worse_than_raw = false,
    )
  end
  angles = angle.(Vdc) .* (180.0 / π)
  slack_angle = angles[slack_idx]
  rel_angles = angles .- slack_angle
  vm = abs.(Vdc)
  finite_rel = filter(isfinite, rel_angles)
  mean_angle = isempty(finite_rel) ? NaN : sum(finite_rel) / length(finite_rel)
  std_angle = isempty(finite_rel) ? NaN : sqrt(sum((a - mean_angle)^2 for a in finite_rel) / length(finite_rel))
  clipped_count = count(a -> isfinite(a) && abs(a) >= dc_angle_limit_deg - 1.0e-8, rel_angles)
  worst_from = missing
  worst_to = missing
  worst_angle = 0.0
  violation_count = 0
  limit = dc_angle_limit_deg
  @inbounds for j in axes(Ybus, 2)
    for ptr in (Ybus isa SparseMatrixCSC ? nzrange(Ybus, j) : axes(Ybus, 1))
      i = Ybus isa SparseMatrixCSC ? rowvals(Ybus)[ptr] : ptr
      i < j || continue
      yij = Ybus isa SparseMatrixCSC ? nonzeros(Ybus)[ptr] : Ybus[i, j]
      iszero(yij) && continue
      Δ = abs((angle(Vdc[i]) - angle(Vdc[j])) * 180.0 / π)
      if isfinite(Δ)
        Δ > limit && (violation_count += 1)
        if Δ > worst_angle
          worst_angle = Δ
          worst_from = i
          worst_to = j
        end
      end
    end
  end
  ratio = isfinite(raw_mis) && raw_mis > 0.0 && isfinite(dc_mis) ? dc_mis / raw_mis : missing
  return (
    dc_angle_min_deg = isempty(finite_rel) ? missing : minimum(finite_rel),
    dc_angle_max_deg = isempty(finite_rel) ? missing : maximum(finite_rel),
    dc_angle_spread_deg = isempty(finite_rel) ? missing : maximum(finite_rel) - minimum(finite_rel),
    dc_angle_mean_deg = isempty(finite_rel) ? missing : mean_angle,
    dc_angle_std_deg = isempty(finite_rel) ? missing : std_angle,
    dc_angle_clipped_count = clipped_count,
    dc_angle_clip_limit_deg = dc_angle_limit_deg,
    dc_max_branch_angle_deg = ismissing(worst_from) ? missing : worst_angle,
    dc_branch_angle_violation_count = violation_count,
    worst_dc_branch_from_bus = worst_from,
    worst_dc_branch_to_bus = worst_to,
    worst_dc_branch_angle_deg = ismissing(worst_from) ? missing : worst_angle,
    dc_voltage_magnitude_min = isempty(vm) ? missing : minimum(vm),
    dc_voltage_magnitude_max = isempty(vm) ? missing : maximum(vm),
    dc_mismatch_ratio_vs_raw = ratio,
    dc_mismatch_growth_factor = ratio,
    requested_dc_worse_than_raw = isfinite(raw_mis) && isfinite(dc_mis) && dc_mis > raw_mis,
  )
end

function _blend_voltage_starts(Vraw::Vector{ComplexF64}, Vdc::Vector{ComplexF64}, λ::Float64, slack_idx::Int)
  0.0 <= λ <= 1.0 || error("blend lambda must satisfy 0 ≤ λ ≤ 1 (got $(λ)).")
  V = similar(Vraw)
  @inbounds for k in eachindex(Vraw)
    if k == slack_idx
      # Never blend away the reference bus complex voltage.
      V[k] = Vraw[k]
      continue
    end
    # Convex interpolation in polar components.
    # Simple and robust for startup candidates; exact geodesic interpolation is unnecessary here.
    vm = (1.0 - λ) * abs(Vraw[k]) + λ * abs(Vdc[k])
    θ = (1.0 - λ) * angle(Vraw[k]) + λ * angle(Vdc[k])
    V[k] = ComplexF64(vm * cos(θ), vm * sin(θ))
  end
  return V
end

"""
Selects a numerically robust rectangular voltage start vector for Newton-Raphson.

The function optionally evaluates multiple start candidates derived from the raw seed:
- sanitized raw start,
- DC-angle start,
- optional raw/DC blended starts.

Default behavior chooses the best finite measured start candidate. Requested DC
behavior is different: when `requested_angle_mode == :dc`, a finite and guarded
DC-angle start is used as the requested baseline, and raw is used only if the
DC candidate is invalid/non-finite or rejected by a start guard. The slack-bus
complex voltage is preserved as reference in all generated candidates.

# Arguments
- `Ybus`: Bus admittance matrix (dense or sparse).
- `Vraw::Vector{ComplexF64}`: Raw complex voltage seed.
- `S::Vector{ComplexF64}`: Specified net complex injections in p.u.
- `bus_types::Vector{Symbol}`: Bus types (`:Slack`, `:PV`, `:PQ`, ...).
- `Vset::Vector{Float64}`: Voltage magnitude setpoints used for regulated buses.
- `slack_idx::Int`: Slack-bus index.

# Keywords
- `enabled::Bool=false`: If `false`, returns `Vraw` unchanged.
- `try_dc_start::Bool=true`: Build/evaluate a DC-angle candidate.
- `try_blend_scan::Bool=true`: Evaluate blend candidates between raw and DC start.
- `branch_guard::Bool=true`: Enforce finite-voltage safety fallback to sanitized raw.
- `measure_candidates::Bool=true`: Compare candidates via mismatch metric.
- `accept_unmeasured_dc_start::Bool=false`: Allow DC candidate selection without mismatch evaluation.
- `blend_lambdas=[0.25, 0.5, 0.75]`: Blend weights for candidate generation.
- `dc_angle_limit_deg::Float64=60.0`: Absolute cap for non-slack relative DC angles.
- `verbose::Int=0`: Logging verbosity.
- `performance_profile=nothing`: Optional profiling dictionary for timing/diagnostics.

# Returns
- `Vector{ComplexF64}`: Selected projected start vector for the rectangular NR solve.
"""
function project_rectangular_start(
  Ybus,
  Vraw::Vector{ComplexF64},
  S::Vector{ComplexF64},
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  slack_idx::Int;
  enabled::Bool = false,
  try_dc_start::Bool = true,
  try_blend_scan::Bool = true,
  branch_guard::Bool = true,
  measure_candidates::Bool = true,
  accept_unmeasured_dc_start::Bool = false,
  blend_lambdas::AbstractVector{<:Real} = [0.25, 0.5, 0.75],
  dc_angle_limit_deg::Float64 = 60.0,
  requested_angle_mode::Symbol = :classic,
  requested_voltage_mode::Symbol = :classic,
  verbose::Int = 0,
  performance_profile = nothing,
)
  enabled || return Vraw
  # Projection is opt-in and must be numerically bounded by a positive angle cap.
  dc_angle_limit_deg > 0.0 || error("dc_angle_limit_deg must be > 0 (got $(dc_angle_limit_deg)).")

  t0 = time_ns()
  # Candidate scan always starts from the sanitized raw seed as baseline.
  candidate_count = 1
  raw = _perf_profile_time!(performance_profile, :start_projection_voltage_clipping) do
    _sanitize_rectangular_start(Vraw, bus_types, Vset, slack_idx)
  end
  best = raw
  best_name = :raw
  raw_mis = measure_candidates ? _perf_profile_time!(performance_profile, :start_projection_mismatch_evaluation) do
    _max_rectangular_mismatch(Ybus, raw, S, bus_types, Vset, slack_idx)
  end : NaN
  best_mis = raw_mis
  selection_reason = measure_candidates ? :raw_baseline : :candidate_mismatch_not_measured
  Vdc = nothing
  dc_mis = NaN
  dc_angle_required = requested_angle_mode === :dc
  dc_angle_start_built = false
  dc_angle_start_valid = false
  dc_angle_start_applied = false
  fallback_to_raw = false
  fallback_reason = missing
  best_blend_mis = NaN

  if try_dc_start || dc_angle_required
    Vdc = _perf_profile_time!(performance_profile, :start_projection_dc_start_construction) do
      _dc_angle_start_rectangular_profiled(Ybus, raw, S, bus_types, Vset, slack_idx; dc_angle_limit_deg = dc_angle_limit_deg, performance_profile = performance_profile)
    end
    dc_angle_start_built = true
    dc_angle_start_valid = all(isfinite, real.(Vdc)) && all(isfinite, imag.(Vdc))
    candidate_count += 1
    dc_mis = measure_candidates ? _perf_profile_time!(performance_profile, :start_projection_mismatch_evaluation) do
      _max_rectangular_mismatch(Ybus, Vdc, S, bus_types, Vset, slack_idx)
    end : NaN
    if dc_angle_required && dc_angle_start_valid
      best = Vdc
      best_name = :dc_start
      best_mis = dc_mis
      dc_angle_start_applied = true
      selection_reason = :requested_dc_angle_start
    elseif dc_angle_required
      best_name = :explicit_fallback_raw
      best_mis = raw_mis
      fallback_to_raw = true
      fallback_reason = :invalid_dc_angle_start
      selection_reason = fallback_reason
    elseif measure_candidates && isfinite(dc_mis) && (!isfinite(best_mis) || dc_mis < best_mis)
      best = Vdc
      best_name = :dc_start
      best_mis = dc_mis
      selection_reason = isfinite(raw_mis) ? :finite_improvement : :finite_candidate_replaces_nonfinite_raw
    elseif !measure_candidates && accept_unmeasured_dc_start
      # Optional policy: prefer structured DC start even without mismatch evaluation.
      best = Vdc
      best_name = :dc_start
      selection_reason = :unmeasured_dc_start_forced
    end
  end

  if try_blend_scan && Vdc !== nothing
    _perf_profile_time!(performance_profile, :start_projection_blend_candidate_generation) do
      for λ_raw in blend_lambdas
        λ = Float64(λ_raw)
        Vblend = _blend_voltage_starts(raw, Vdc, λ, slack_idx)
        candidate_count += 1
        blend_mis = measure_candidates ? _perf_profile_time!(performance_profile, :start_projection_mismatch_evaluation) do
          _max_rectangular_mismatch(Ybus, Vblend, S, bus_types, Vset, slack_idx)
        end : NaN
        if isfinite(blend_mis) && (!isfinite(best_blend_mis) || blend_mis < best_blend_mis)
          best_blend_mis = blend_mis
        end
        if !dc_angle_required && measure_candidates && isfinite(blend_mis) && (!isfinite(best_mis) || blend_mis < best_mis)
          best = Vblend
          best_name = Symbol("blend_", λ)
          best_mis = blend_mis
          selection_reason = isfinite(raw_mis) ? :finite_improvement : :finite_candidate_replaces_nonfinite_raw
        end
      end
    end
  end

  if branch_guard
    _perf_profile_time!(performance_profile, :start_projection_branch_guard_checks) do
      # Final safety belt: never return a non-finite candidate into NR.
      # If a candidate became invalid, deterministically fall back to sanitized raw.
      if !all(isfinite, real.(best)) || !all(isfinite, imag.(best))
        best = raw
        best_name = dc_angle_required ? :explicit_fallback_raw : :raw
        best_mis = raw_mis
        selection_reason = dc_angle_required ? :invalid_dc_angle_start : :nonfinite_selected_voltage
        fallback_to_raw = dc_angle_required
        fallback_reason = dc_angle_required ? :invalid_dc_angle_start : fallback_reason
        dc_angle_start_applied = false
      end
    end
  end

  if best_name === :raw && measure_candidates && selection_reason === :raw_baseline
    # Distinguish "raw chosen as initial baseline" from "raw won final comparison".
    selection_reason = :no_finite_improvement
  end
  reported_best_mis = isfinite(best_mis) ? best_mis : missing

  _perf_profile_time!(performance_profile, :start_projection_final_selection) do
    nothing
  end
  dc_quality = _dc_start_quality_diagnostics(Ybus, Vdc, raw_mis, dc_mis, slack_idx, dc_angle_limit_deg)
  if performance_profile isa AbstractDict && Bool(get(performance_profile, :enabled, false))
    # Compact selection summary for diagnostics/UI without storing all candidate vectors.
    performance_profile[:start_projection_summary] = (
      selected = best_name,
      reason = selection_reason,
      candidates = candidate_count,
      best_mismatch = reported_best_mis,
      raw_mismatch = isfinite(raw_mis) ? raw_mis : missing,
      elapsed_s = (time_ns() - t0) / 1e9,
      requested_angle_mode = requested_angle_mode,
      requested_voltage_mode = requested_voltage_mode,
      dc_angle_start_built = dc_angle_start_built,
      dc_angle_start_valid = dc_angle_start_valid,
      dc_angle_start_applied = dc_angle_start_applied,
      selected_start_candidate = best_name,
      selection_reason = selection_reason,
      dc_mismatch = isfinite(dc_mis) ? dc_mis : missing,
      best_blend_mismatch = isfinite(best_blend_mis) ? best_blend_mis : missing,
      projected_mismatch = reported_best_mis,
      fallback_to_raw = fallback_to_raw,
      fallback_reason = fallback_reason,
      raw_fallback_reason = fallback_reason,
      start_projection_mismatch_before = isfinite(raw_mis) ? raw_mis : missing,
      start_projection_mismatch_after = reported_best_mis,
      dc_quality...,
    )
  end

  if verbose > 0
    @info "start projection selected $(best_name)" requested_angle_mode = requested_angle_mode requested_voltage_mode = requested_voltage_mode reason = selection_reason raw_mismatch = (isfinite(raw_mis) ? raw_mis : missing) dc_mismatch = (isfinite(dc_mis) ? dc_mis : missing) projected_mismatch = reported_best_mis dc_angle_start_built = dc_angle_start_built dc_angle_start_valid = dc_angle_start_valid dc_angle_start_applied = dc_angle_start_applied fallback_to_raw = fallback_to_raw fallback_reason = fallback_reason dc_angle_spread_deg = dc_quality.dc_angle_spread_deg dc_max_branch_angle_deg = dc_quality.dc_max_branch_angle_deg dc_mismatch_ratio_vs_raw = dc_quality.dc_mismatch_ratio_vs_raw
  end
  return best
end

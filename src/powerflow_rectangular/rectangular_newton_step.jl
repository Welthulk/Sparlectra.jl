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
# Rectangular power-flow Newton-step and autodamping helpers.
#
# This file is included inside module Sparlectra. Do not add a module wrapper here.

# Date: 29.5.2026
# file: src/powerflow_rectangular/rectangular_newton_step.jl

function _validate_rectangular_damping(damp::Float64, autodamp_min::Float64)
  # Guard rails keep damping in the physically meaningful [0,1] step-length range.
  isfinite(damp) || error("damp must be finite (got $(damp)).")
  isfinite(autodamp_min) || error("autodamp_min must be finite (got $(autodamp_min)).")
  0.0 < damp <= 1.0 || error("damp must satisfy 0 < damp ≤ 1 (got $(damp)).")
  # autodamp_min is the backtracking floor and must not exceed the initial step.
  0.0 < autodamp_min <= damp || error("autodamp_min must satisfy 0 < autodamp_min ≤ damp (got autodamp_min=$(autodamp_min), damp=$(damp)).")
  return nothing
end

# Allocation-light mismatch norm helper used in hot Newton/autodamp paths.
@inline _max_abs_mismatch(F::AbstractVector{<:Real}) = mapreduce(abs, max, F; init = 0.0)

@inline function _apply_rectangular_delta!(Vout::Vector{ComplexF64}, V::Vector{ComplexF64}, δx::Vector{Float64}, slack_idx::Int, non_slack::Vector{Int}, alpha::Float64)
  n = length(V)
  # Start from the current iterate so untouched entries (including potential metadata in V)
  # remain bitwise-consistent unless explicitly updated below.
  copyto!(Vout, V)

  # State layout convention:
  #   δx[1:(n-1)]     -> ΔVr for non-slack buses
  #   δx[n:(2n-2)]    -> ΔVi for non-slack buses
  # where "idx" enumerates non_slack in solver order.
  @inbounds for (idx, bus) in enumerate(non_slack)
    Vout[bus] += ComplexF64(alpha * δx[idx], alpha * δx[(n-1)+idx])
  end

  # Preserve slack reference exactly; non-slack updates only.
  Vout[slack_idx] = V[slack_idx]
  return Vout
end

function _apply_rectangular_delta(V::Vector{ComplexF64}, δx::Vector{Float64}, slack_idx::Int, non_slack::Vector{Int}, alpha::Float64)
  # Convenience wrapper for call sites that need an owned output vector.
  Vnew = similar(V)
  return _apply_rectangular_delta!(Vnew, V, δx, slack_idx, non_slack, alpha)
end

function _max_rectangular_mismatch(Ybus, V::Vector{ComplexF64}, S::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int)
  # Keep mismatch computation centralized so acceptance logic stays consistent
  # with the Newton residual definition.
  F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
  return _max_abs_mismatch(F)
end

# Weighted squared residual norm ‖W F‖² for the merit function f(x) = 1/2 ‖W F(x)‖².
# F is stacked as [ΔP_i, second-equation_i] per non-slack bus i (see mismatch_rectangular);
# the second equation is ΔQ_i for PQ buses and ΔV_i for PV buses, so scale_q/scale_v are
# applied positionally via bus_types, not via the raw F index.
@inline function _rectangular_merit_weighted_sqnorm(F::AbstractVector{Float64}, non_slack::Vector{Int}, bus_types::Vector{Symbol}, scale_p::Float64, scale_q::Float64, scale_v::Float64)
  total = 0.0
  @inbounds for (k, bus) in enumerate(non_slack)
    row = 2k - 1
    wp = scale_p * F[row]
    total += wp * wp
    w2 = bus_types[bus] == :PV ? scale_v : scale_q
    wq = w2 * F[row+1]
    total += wq * wq
  end
  return total
end

# Merit function value f(x) = 1/2 ‖W F(x)‖². See docs/src/solver.md, "Merit-Function Line Search".
@inline _rectangular_merit_value(F::AbstractVector{Float64}, non_slack::Vector{Int}, bus_types::Vector{Symbol}, scale_p::Float64, scale_q::Float64, scale_v::Float64) =
  0.5 * _rectangular_merit_weighted_sqnorm(F, non_slack, bus_types, scale_p, scale_q, scale_v)

# Classic max-mismatch backtracking search: the exact original `choose_rectangular_autodamp`
# algorithm, factored out so both the default (merit-disabled) path and the merit path's
# active-set-skip fallback share one implementation.
function _classic_rectangular_autodamp_search(Ybus, V::Vector{ComplexF64}, S::Vector{ComplexF64}, δx::Vector{Float64}, non_slack::Vector{Int}, current_mismatch::Float64; slack_idx::Int, damp::Float64, autodamp_min::Float64, bus_types::Vector{Symbol}, Vset::Vector{Float64}, diagnostics)
  Vtrial = similar(V)
  best_V = similar(V)

  # Initialize fallback with the minimum admissible step so we always return
  # a conservative finite candidate even when no improving step is found.
  best_alpha = autodamp_min
  _apply_rectangular_delta!(Vtrial, V, δx, slack_idx, non_slack, autodamp_min)
  best_mismatch = _max_rectangular_mismatch(Ybus, Vtrial, S, bus_types, Vset, slack_idx)
  copyto!(best_V, Vtrial)

  alpha = damp
  # Monotone backtracking (α, α/2, α/4, ...) until we hit the configured floor.
  while alpha >= autodamp_min
    _apply_rectangular_delta!(Vtrial, V, δx, slack_idx, non_slack, alpha)
    trial_mismatch = _max_rectangular_mismatch(Ybus, Vtrial, S, bus_types, Vset, slack_idx)

    # Accept first strict improvement over current mismatch for predictable behavior.
    if isfinite(trial_mismatch) && trial_mismatch < current_mismatch
      # Return an owned vector; caller should not alias the internal scratch buffer.
      diagnostics isa AbstractVector && push!(diagnostics, (alpha = alpha, trial_mismatch = trial_mismatch, accepted_improvement = true))
      return alpha, copy(Vtrial), trial_mismatch
    end
    # Track best finite non-improving candidate as a safe fallback.
    if isfinite(trial_mismatch) && trial_mismatch < best_mismatch
      best_alpha = alpha
      best_mismatch = trial_mismatch
      copyto!(best_V, Vtrial)
    end
    alpha *= 0.5
  end

  # Conservative fallback: best finite trial seen during backtracking.
  diagnostics isa AbstractVector && push!(diagnostics, (alpha = best_alpha, trial_mismatch = best_mismatch, accepted_improvement = false))
  return best_alpha, best_V, best_mismatch
end

"""
    choose_rectangular_autodamp(Ybus, V, S, δx, F0; slack_idx, damp, autodamp_min, bus_types, Vset)

Select a Newton step length for the rectangular power-flow solver by backtracking
from `damp` toward `autodamp_min`. The first trial step that reduces the maximum
absolute mismatch is accepted. If no trial reduces the mismatch, the smallest
finite trial is returned so the solver can continue safely with a conservative
step.

When `merit_enabled = true` (opt-in, requires `autodamp = true`), an Armijo
sufficient-decrease criterion on the merit function `f(x) = 1/2 ‖W F(x)‖²`
is tried first for each trial step, in the same backtracking order; see
[Merit-Function Line Search](@ref) for the theoretical background. If
`active_set_changed = true` (a PV/PQ switch happened this Newton iteration),
the merit comparison is skipped for this call and the classic max-mismatch
criterion is used instead, since the residual vector's entries change meaning
across an active-set switch. If no trial satisfies Armijo, `fallback_max_mismatch`
selects whether to fall back to the classic max-mismatch criterion (`true`) or
directly to the conservative best-finite-trial fallback (`false`).

Returns `(alpha, Vtrial, trial_mismatch)`. When `merit_enabled = true`, per-trial
diagnostics (`f_before`, `directional_derivative`, tested/accepted λ, `accept_reason`)
are pushed to `merit_log` if it is an `AbstractVector`.
"""
function choose_rectangular_autodamp(
  Ybus,
  V::Vector{ComplexF64},
  S::Vector{ComplexF64},
  δx::Vector{Float64},
  F0::Vector{Float64};
  slack_idx::Int,
  damp::Float64 = 1.0,
  autodamp_min::Float64 = 0.05,
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  diagnostics = nothing,
  merit_enabled::Bool = false,
  armijo_c1::Float64 = 1.0e-4,
  scale_p::Float64 = 1.0,
  scale_q::Float64 = 1.0,
  scale_v::Float64 = 1.0,
  fallback_max_mismatch::Bool = true,
  active_set_changed::Bool = false,
  merit_log = nothing,
)
  _validate_rectangular_damping(damp, autodamp_min)
  non_slack = non_slack_indices(length(V), slack_idx)
  current_mismatch = _max_abs_mismatch(F0)

  if !merit_enabled
    # Default path: bit-identical to the pre-merit implementation, no extra
    # allocations or computations.
    return _classic_rectangular_autodamp_search(Ybus, V, S, δx, non_slack, current_mismatch; slack_idx = slack_idx, damp = damp, autodamp_min = autodamp_min, bus_types = bus_types, Vset = Vset, diagnostics = diagnostics)
  end

  if active_set_changed
    # PV/PQ switching changed what the residual entries mean at this bus; a
    # merit comparison across the switch is not well-defined for this iteration.
    alpha, Vout, mismatch = _classic_rectangular_autodamp_search(Ybus, V, S, δx, non_slack, current_mismatch; slack_idx = slack_idx, damp = damp, autodamp_min = autodamp_min, bus_types = bus_types, Vset = Vset, diagnostics = diagnostics)
    merit_log isa AbstractVector && push!(merit_log, (f_before = NaN, directional_derivative = NaN, f_after = NaN, tested_alphas = Float64[], accepted_alpha = alpha, accept_reason = :active_set_skip))
    return alpha, Vout, mismatch
  end

  # --- Merit-function Armijo line search --------------------------------
  f0 = _rectangular_merit_value(F0, non_slack, bus_types, scale_p, scale_q, scale_v)
  # ∇f = JᵀWᵀWF and Δx = −J⁻¹F, so ∇fᵀΔx = −‖WF‖² without an extra Jacobian matvec.
  directional_derivative = -2.0 * f0
  tested_alphas = Float64[]

  Vtrial = similar(V)
  best_V = similar(V)
  best_alpha = autodamp_min
  _apply_rectangular_delta!(Vtrial, V, δx, slack_idx, non_slack, autodamp_min)
  F_base = mismatch_rectangular(Ybus, Vtrial, S, bus_types, Vset, slack_idx)
  best_mismatch = _max_abs_mismatch(F_base)
  best_f = _rectangular_merit_value(F_base, non_slack, bus_types, scale_p, scale_q, scale_v)
  copyto!(best_V, Vtrial)

  first_improving_alpha = nothing
  first_improving_V = nothing
  first_improving_mismatch = NaN
  first_improving_f = NaN

  alpha = damp
  while alpha >= autodamp_min
    _apply_rectangular_delta!(Vtrial, V, δx, slack_idx, non_slack, alpha)
    F_trial = mismatch_rectangular(Ybus, Vtrial, S, bus_types, Vset, slack_idx)
    trial_mismatch = _max_abs_mismatch(F_trial)
    f_trial = _rectangular_merit_value(F_trial, non_slack, bus_types, scale_p, scale_q, scale_v)
    push!(tested_alphas, alpha)

    if isfinite(f_trial) && f_trial <= f0 + armijo_c1 * alpha * directional_derivative
      diagnostics isa AbstractVector && push!(diagnostics, (alpha = alpha, trial_mismatch = trial_mismatch, accepted_improvement = isfinite(trial_mismatch) && trial_mismatch < current_mismatch))
      merit_log isa AbstractVector && push!(merit_log, (f_before = f0, directional_derivative = directional_derivative, f_after = f_trial, tested_alphas = copy(tested_alphas), accepted_alpha = alpha, accept_reason = :armijo))
      return alpha, copy(Vtrial), trial_mismatch
    end

    if first_improving_alpha === nothing && isfinite(trial_mismatch) && trial_mismatch < current_mismatch
      first_improving_alpha = alpha
      first_improving_V = copy(Vtrial)
      first_improving_mismatch = trial_mismatch
      first_improving_f = f_trial
    end
    if isfinite(trial_mismatch) && trial_mismatch < best_mismatch
      best_alpha = alpha
      best_mismatch = trial_mismatch
      best_f = f_trial
      copyto!(best_V, Vtrial)
    end
    alpha *= 0.5
  end

  # No trial satisfied the Armijo condition.
  if fallback_max_mismatch && first_improving_alpha !== nothing
    diagnostics isa AbstractVector && push!(diagnostics, (alpha = first_improving_alpha, trial_mismatch = first_improving_mismatch, accepted_improvement = true))
    merit_log isa AbstractVector && push!(merit_log, (f_before = f0, directional_derivative = directional_derivative, f_after = first_improving_f, tested_alphas = copy(tested_alphas), accepted_alpha = first_improving_alpha, accept_reason = :fallback_max_mismatch))
    return first_improving_alpha, first_improving_V, first_improving_mismatch
  end

  diagnostics isa AbstractVector && push!(diagnostics, (alpha = best_alpha, trial_mismatch = best_mismatch, accepted_improvement = false))
  merit_log isa AbstractVector && push!(merit_log, (f_before = f0, directional_derivative = directional_derivative, f_after = best_f, tested_alphas = copy(tested_alphas), accepted_alpha = best_alpha, accept_reason = :fallback_conservative))
  return best_alpha, best_V, best_mismatch
end

"""
    RectangularTrustRegionCollapsed

Thrown by `choose_rectangular_trust_region_step` when the trust-region radius
falls below `min_radius` without ever accepting a step. Caught by the outer
Newton-Raphson loop (`rectangular_network_solver.jl`) alongside the existing
singular-Jacobian failure handling, and reported as non-convergence with
reason `:trust_region_collapsed`.
"""
struct RectangularTrustRegionCollapsed <: Exception end

"""
    choose_rectangular_trust_region_step(Ybus, V, S, δx, F0, J; slack_idx, bus_types, Vset, radius_ref, min_radius, max_radius, eta_accept, shrink_factor, expand_factor, expand_threshold)

Scaled-Newton trust-region step control (dogleg/Steihaug is out of scope).
The full Newton step `δx` is scaled down to the current radius when it
exceeds it; the trial is accepted when the ratio `rho` of actual to
predicted merit-function reduction is at least `eta_accept`. Rejected trials
shrink the radius and retry within the same Newton iteration (bounded by
`min_radius`); the radius persists across iterations via `radius_ref`, which
is mutated in place.

Reuses the exact merit function from the merit-function line search
(`_rectangular_merit_value`, unweighted `W = I`) as `m(x)`, matching
`m(x) = ‖F(x)‖²` up to the constant `1/2` factor common to numerator and
denominator of `rho` — no second merit computation is implemented.

Throws [`RectangularTrustRegionCollapsed`](@ref) if the radius falls below
`min_radius` without an accepted step. Returns the accepted trial voltage
vector.
"""
function choose_rectangular_trust_region_step(
  Ybus,
  V::Vector{ComplexF64},
  S::Vector{ComplexF64},
  δx::Vector{Float64},
  F0::Vector{Float64},
  J;
  slack_idx::Int,
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  radius_ref::Base.RefValue{Float64},
  min_radius::Float64,
  max_radius::Float64,
  eta_accept::Float64,
  shrink_factor::Float64,
  expand_factor::Float64,
  expand_threshold::Float64,
  diagnostics = nothing,
)
  non_slack = non_slack_indices(length(V), slack_idx)
  # W = I: m(x) = 1/2‖F(x)‖², the same function the merit-function line search uses.
  m_x = _rectangular_merit_value(F0, non_slack, bus_types, 1.0, 1.0, 1.0)
  dx_norm = norm(δx)

  Vtrial = similar(V)
  tested_radii = Float64[]
  rejected_steps = 0

  while true
    radius = radius_ref[]
    if radius < min_radius
      diagnostics isa AbstractVector && push!(diagnostics, (radius_before = radius, rho = NaN, tested_radii = copy(tested_radii), rejected_steps = rejected_steps, accepted = false, radius_after = radius, collapsed = true))
      throw(RectangularTrustRegionCollapsed())
    end
    push!(tested_radii, radius)

    # Scaled-Newton: cap the step norm at the radius, direction unchanged (no dogleg).
    scale = dx_norm > radius ? radius / dx_norm : 1.0
    _apply_rectangular_delta!(Vtrial, V, δx, slack_idx, non_slack, scale)
    F_trial = mismatch_rectangular(Ybus, Vtrial, S, bus_types, Vset, slack_idx)
    m_trial = _rectangular_merit_value(F_trial, non_slack, bus_types, 1.0, 1.0, 1.0)

    if isfinite(m_trial)
      # Predicted reduction from the linear model F(x) + J*dx_scaled, using the
      # same merit-value convention as m_x/m_trial for a consistent ratio.
      predicted_F = F0 .+ J * (scale .* δx)
      predicted_reduction = m_x - _rectangular_merit_value(predicted_F, non_slack, bus_types, 1.0, 1.0, 1.0)
      actual_reduction = m_x - m_trial
      rho = actual_reduction / max(predicted_reduction, eps())

      if rho >= eta_accept
        hit_boundary = dx_norm > radius
        if rho >= expand_threshold && hit_boundary
          radius_ref[] = min(radius * expand_factor, max_radius)
        end
        diagnostics isa AbstractVector && push!(diagnostics, (radius_before = radius, rho = rho, tested_radii = copy(tested_radii), rejected_steps = rejected_steps, accepted = true, radius_after = radius_ref[], collapsed = false))
        return copy(Vtrial)
      end
    end

    # Non-finite trial merit or rho < eta_accept: reject, shrink, retry.
    rejected_steps += 1
    radius_ref[] = radius * shrink_factor
  end
end

"""
    complex_newton_step_rectangular(
        Ybus,
        V,
        S;
        slack_idx,
        damp,
        bus_types,
        Vset,
    )

Performs one Newton–Raphson step in rectangular coordinates using the analytic
Jacobian that matches `mismatch_rectangular`.

- State: x = [Vr(non-slack); Vi(non-slack)]
- Residual: F(x) = mismatch_rectangular(...)
"""
function complex_newton_step_rectangular(
  Ybus,
  V::Vector{ComplexF64},
  S::Vector{ComplexF64};
  slack_idx::Int,
  damp::Float64 = 1.0,
  autodamp::Bool = false,
  autodamp_min::Float64 = 0.05,
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  dPinj_dVm::Vector{Float64} = zeros(Float64, length(V)),
  dQinj_dVm::Vector{Float64} = zeros(Float64, length(V)),
  performance_profile = nothing,
  step_diagnostics = nothing,
  merit_enabled::Bool = false,
  armijo_c1::Float64 = 1.0e-4,
  scale_p::Float64 = 1.0,
  scale_q::Float64 = 1.0,
  scale_v::Float64 = 1.0,
  fallback_max_mismatch::Bool = true,
  active_set_changed::Bool = false,
  merit_log = nothing,
  trust_region_enabled::Bool = false,
  tr_radius_ref::Union{Nothing,Base.RefValue{Float64}} = nothing,
  tr_min_radius::Float64 = 1e-4,
  tr_max_radius::Float64 = 10.0,
  tr_eta_accept::Float64 = 0.1,
  tr_shrink_factor::Float64 = 0.5,
  tr_expand_factor::Float64 = 2.0,
  tr_expand_threshold::Float64 = 0.75,
  tr_log = nothing,
)
  n = length(V)
  # Solver assumes state ordering [Vr(non-slack); Vi(non-slack)] consistently
  # across mismatch, Jacobian build, and delta-application helpers.
  @assert length(S) == n
  @assert length(bus_types) == n
  @assert length(Vset) == n
  # Callers can pass preallocated derivative buffers to avoid default allocations.
  @assert length(dPinj_dVm) == n
  @assert length(dQinj_dVm) == n

  # Non-slack bus order defines both state-vector and Jacobian block layout.
  non_slack = non_slack_indices(n, slack_idx)

  # Residual matching the FD variant.
  F0 = _perf_profile_time!(performance_profile, :newton_step_mismatch) do
    mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
  end
  m = length(F0)
  nvar = 2 * (n - 1)
  @assert m == nvar "complex_newton_step_rectangular: mismatch and state dimension differ"

  # Analytic sparse Jacobian aligned with rectangular state ordering.
  J = _perf_profile_time!(performance_profile, :newton_step_jacobian) do
    build_rectangular_jacobian_pq_pv(Ybus, V, bus_types, Vset, slack_idx; dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm)
  end

  # Solve J * δx = -F.
  δx = _perf_profile_time!(performance_profile, :newton_step_linear_solve) do
    solve_linear(J, -F0; allow_pinv = true)
  end
  if autodamp
    # Autodamp path evaluates multiple trial voltages and returns accepted/fallback trial.
    _, Vtrial, _ = _perf_profile_time!(performance_profile, :newton_step_autodamp) do
      choose_rectangular_autodamp(
        Ybus, V, S, δx, F0;
        slack_idx = slack_idx, damp = damp, autodamp_min = autodamp_min, bus_types = bus_types, Vset = Vset, diagnostics = step_diagnostics,
        merit_enabled = merit_enabled, armijo_c1 = armijo_c1, scale_p = scale_p, scale_q = scale_q, scale_v = scale_v,
        fallback_max_mismatch = fallback_max_mismatch, active_set_changed = active_set_changed, merit_log = merit_log,
      )
    end
    return Vtrial
  end

  if trust_region_enabled
    # Scaled-Newton trust-region path: caps ||δx|| at an adaptive radius and
    # accepts/rejects by merit-function decrease. Mutually exclusive with
    # autodamp at the config layer. May throw RectangularTrustRegionCollapsed,
    # caught by the outer NR loop as a non-convergence reason.
    return _perf_profile_time!(performance_profile, :newton_step_trust_region) do
      choose_rectangular_trust_region_step(
        Ybus, V, S, δx, F0, J;
        slack_idx = slack_idx, bus_types = bus_types, Vset = Vset,
        radius_ref = tr_radius_ref, min_radius = tr_min_radius, max_radius = tr_max_radius,
        eta_accept = tr_eta_accept, shrink_factor = tr_shrink_factor, expand_factor = tr_expand_factor,
        expand_threshold = tr_expand_threshold, diagnostics = tr_log,
      )
    end
  end

  # Fixed damping path: single update with validated step length.
  _validate_rectangular_damping(damp, min(autodamp_min, damp))
  Vnext = similar(V)
  step_diagnostics isa AbstractVector && push!(step_diagnostics, (alpha = damp, trial_mismatch = NaN, accepted_improvement = true))
  return _apply_rectangular_delta!(Vnext, V, δx, slack_idx, non_slack, damp)
end

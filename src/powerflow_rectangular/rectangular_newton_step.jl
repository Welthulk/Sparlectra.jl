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

"""
    choose_rectangular_autodamp(Ybus, V, S, δx, F0; slack_idx, damp, autodamp_min, bus_types, Vset)

Select a Newton step length for the rectangular power-flow solver by backtracking
from `damp` toward `autodamp_min`. The first trial step that reduces the maximum
absolute mismatch is accepted. If no trial reduces the mismatch, the smallest
finite trial is returned so the solver can continue safely with a conservative
step.

Returns `(alpha, Vtrial, trial_mismatch)`.
"""
function choose_rectangular_autodamp(Ybus, V::Vector{ComplexF64}, S::Vector{ComplexF64}, δx::Vector{Float64}, F0::Vector{Float64}; slack_idx::Int, damp::Float64 = 1.0, autodamp_min::Float64 = 0.05, bus_types::Vector{Symbol}, Vset::Vector{Float64})
  _validate_rectangular_damping(damp, autodamp_min)
  non_slack = non_slack_indices(length(V), slack_idx)
  current_mismatch = _max_abs_mismatch(F0)

  # Reused trial buffers avoid per-backtracking-step allocations.
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
  return best_alpha, best_V, best_mismatch
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
      choose_rectangular_autodamp(Ybus, V, S, δx, F0; slack_idx = slack_idx, damp = damp, autodamp_min = autodamp_min, bus_types = bus_types, Vset = Vset)
    end
    return Vtrial
  end

  # Fixed damping path: single update with validated step length.
  _validate_rectangular_damping(damp, min(autodamp_min, damp))
  Vnext = similar(V)
  return _apply_rectangular_delta!(Vnext, V, δx, slack_idx, non_slack, damp)
end

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

function _validate_rectangular_damping(damp::Float64, autodamp_min::Float64)
  # Guard rails keep damping in the physically meaningful [0,1] step-length range.
  isfinite(damp) || error("damp must be finite (got $(damp)).")
  isfinite(autodamp_min) || error("autodamp_min must be finite (got $(autodamp_min)).")
  0.0 < damp <= 1.0 || error("damp must satisfy 0 < damp ≤ 1 (got $(damp)).")
  0.0 < autodamp_min <= damp || error("autodamp_min must satisfy 0 < autodamp_min ≤ damp (got autodamp_min=$(autodamp_min), damp=$(damp)).")
  return nothing
end

function _apply_rectangular_delta(V::Vector{ComplexF64}, δx::Vector{Float64}, slack_idx::Int, non_slack::Vector{Int}, alpha::Float64)
  n = length(V)
  Vr = real.(V)
  Vi = imag.(V)
  Vr_new = copy(Vr)
  Vi_new = copy(Vi)

  @inbounds for (idx, bus) in enumerate(non_slack)
    Vr_new[bus] += alpha * δx[idx]
    Vi_new[bus] += alpha * δx[(n-1)+idx]
  end

  Vr_new[slack_idx] = Vr[slack_idx]
  Vi_new[slack_idx] = Vi[slack_idx]
  # Preserve slack reference exactly; non-slack updates only.
  return ComplexF64.(Vr_new, Vi_new)
end

function _max_rectangular_mismatch(Ybus, V::Vector{ComplexF64}, S::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int)
  F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
  return maximum(abs.(F))
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
function choose_rectangular_autodamp(Ybus, V::Vector{ComplexF64}, S::Vector{ComplexF64}, δx::Vector{Float64}, F0::Vector{Float64}; slack_idx::Int, damp::Float64 = 1.0, autodamp_min::Float64 = 1e-3, bus_types::Vector{Symbol}, Vset::Vector{Float64})
  _validate_rectangular_damping(damp, autodamp_min)
  non_slack = non_slack_indices(length(V), slack_idx)
  current_mismatch = maximum(abs.(F0))
  best_alpha = autodamp_min
  best_V = _apply_rectangular_delta(V, δx, slack_idx, non_slack, autodamp_min)
  best_mismatch = _max_rectangular_mismatch(Ybus, best_V, S, bus_types, Vset, slack_idx)

  alpha = damp
  # Backtrack until mismatch improves or we hit the minimum configured step.
  while alpha >= autodamp_min
    Vtrial = _apply_rectangular_delta(V, δx, slack_idx, non_slack, alpha)
    trial_mismatch = _max_rectangular_mismatch(Ybus, Vtrial, S, bus_types, Vset, slack_idx)
    if isfinite(trial_mismatch) && trial_mismatch < current_mismatch
      return alpha, Vtrial, trial_mismatch
    end
    if isfinite(trial_mismatch) && trial_mismatch < best_mismatch
      best_alpha = alpha
      best_V = Vtrial
      best_mismatch = trial_mismatch
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
  autodamp_min::Float64 = 1e-3,
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  dPinj_dVm::Vector{Float64} = zeros(Float64, length(V)),
  dQinj_dVm::Vector{Float64} = zeros(Float64, length(V)),
  performance_profile = nothing,
)
  n = length(V)
  @assert length(S) == n
  @assert length(bus_types) == n
  @assert length(Vset) == n
  @assert length(dPinj_dVm) == n
  @assert length(dQinj_dVm) == n

  # Non-slack indices
  non_slack = non_slack_indices(n, slack_idx)
  # Residual matching the FD variant
  F0 = _perf_profile_time!(performance_profile, :newton_step_mismatch) do
    mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
  end
  m = length(F0)
  nvar = 2 * (n - 1)
  @assert m == nvar "complex_newton_step_rectangular: mismatch and state dimension differ"

  # Analytic rectangular sparse Jacobian
  J = _perf_profile_time!(performance_profile, :newton_step_jacobian) do
    build_rectangular_jacobian_pq_pv(Ybus, V, bus_types, Vset, slack_idx; dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm)
  end

  # Solve J * δx = -F
  δx = _perf_profile_time!(performance_profile, :newton_step_linear_solve) do
    solve_linear(J, -F0; allow_pinv = true)
  end
  if autodamp
    _, Vtrial, _ = _perf_profile_time!(performance_profile, :newton_step_autodamp) do
      choose_rectangular_autodamp(Ybus, V, S, δx, F0; slack_idx = slack_idx, damp = damp, autodamp_min = autodamp_min, bus_types = bus_types, Vset = Vset)
    end
    return Vtrial
  end

  _validate_rectangular_damping(damp, min(autodamp_min, damp))
  return _apply_rectangular_delta(V, δx, slack_idx, non_slack, damp)
end

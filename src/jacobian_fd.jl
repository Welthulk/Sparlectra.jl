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

# file: src/jacobian_fd.jl
# jacobian_fd.jl — Finite-difference Jacobian step for rectangular complex-state NR
#
# This file provides the finite-difference Newton step
#   complex_newton_step_rectangular_fd(...)
#
# It assumes that the following are already available in the module:
#   - mismatch_rectangular(...)
#   - type Net, etc. (only indirectly via callers)
#
# Dependencies:
#   - LinearAlgebra (for pinv and SingularException)

using LinearAlgebra

"""
    complex_newton_step_rectangular_fd(
        Ybus,
        V,
        S;
        slack_idx::Int = 1,
        damp::Float64  = 1.0,
        h::Float64     = 1e-6,
        bus_types::Vector{Symbol},
        Vset::Vector{Float64},
    ) -> Vector{ComplexF64}

Performs one Newton–Raphson step in rectangular coordinates using a
finite-difference Jacobian on the mismatch

    F(V) = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)

which corresponds to a stacked vector of
- ΔP, ΔQ for PQ buses
- ΔP, ΔV for PV buses
over all non-slack buses.

Arguments:
- `Ybus`: bus admittance matrix (n×n, Complex)
- `V`: current complex bus voltage vector (length n)
- `S`: specified complex power injections P + jQ (length n)
- `slack_idx`: index of the slack bus (no equations / no variables)
- `damp`: scalar damping factor for the Newton step (0 < damp ≤ 1)
- `autodamp`: when true, backtrack from `damp` to reduce the current mismatch
- `autodamp_min`: minimum automatic damping factor
- `h`: perturbation step for finite differences
- `bus_types`: vector of bus types (:PQ, :PV, :Slack)
- `Vset`: voltage magnitude setpoints for PV buses

Returns:
- Updated complex voltage vector `V_new` (length n)
"""
function complex_newton_step_rectangular_fd(Ybus, V, S; slack_idx::Int = 1, damp::Float64 = 1.0, autodamp::Bool = false, autodamp_min::Float64 = 1e-3, h::Float64 = 1e-6, bus_types::Vector{Symbol}, Vset::Vector{Float64}, dPinj_dVm::Vector{Float64} = zeros(Float64, length(V)), dQinj_dVm::Vector{Float64} = zeros(Float64, length(V)), performance_profile = nothing)
  n = length(V)

  # Base mismatch F(V)
  F0 = _perf_profile_time!(performance_profile, :newton_step_fd_base_mismatch) do
    mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
  end
  m  = length(F0)  # expected = 2 * (n-1)

  # Non-slack buses
  non_slack = non_slack_indices(n, slack_idx)

  # Variables: Vr(non-slack) and Vi(non-slack)
  nvar = 2 * (n - 1)
  @assert nvar == m "Rectangular FD-Newton: nvar and m should both equal 2*(n-1)"

  J = _perf_profile_time!(performance_profile, :newton_step_fd_jacobian) do
    zeros(Float64, m, nvar)
  end

  Vr = real.(V)
  Vi = imag.(V)

  # Columns 1..(n-1): perturb Vr(non_slack[k])
  for (col_idx, bus) in enumerate(non_slack)
    V_pert = copy(V)
    V_pert[bus] = ComplexF64(Vr[bus] + h, Vi[bus])

    Fp = mismatch_rectangular(Ybus, V_pert, S, bus_types, Vset, slack_idx)
    J[:, col_idx] .= (Fp .- F0) ./ h
  end

  # Columns n..2(n-1): perturb Vi(non_slack[k])
  for (offset, bus) in enumerate(non_slack)
    col_idx = (n - 1) + offset
    V_pert = copy(V)
    V_pert[bus] = ComplexF64(Vr[bus], Vi[bus] + h)

    Fp = mismatch_rectangular(Ybus, V_pert, S, bus_types, Vset, slack_idx)
    J[:, col_idx] .= (Fp .- F0) ./ h
  end

  # Solve J * δx = -F0
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

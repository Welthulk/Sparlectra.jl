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
- `h`: perturbation step for finite differences
- `bus_types`: vector of bus types (:PQ, :PV, :Slack)
- `Vset`: voltage magnitude setpoints for PV buses

Returns:
- Updated complex voltage vector `V_new` (length n)
"""
function complex_newton_step_rectangular_fd(Ybus, V, S; slack_idx::Int = 1, damp::Float64 = 1.0, h::Float64 = 1e-6, bus_types::Vector{Symbol}, Vset::Vector{Float64})
  n = length(V)

  # Base mismatch F(V)
  F0 = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
  m  = length(F0)  # expected = 2 * (n-1)

  # Non-slack buses
  non_slack = collect(1:n)
  deleteat!(non_slack, slack_idx)

  # Variables: Vr(non-slack) and Vi(non-slack)
  nvar = 2 * (n - 1)
  @assert nvar == m "Rectangular FD-Newton: nvar and m should both equal 2*(n-1)"

  J = zeros(Float64, m, nvar)

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
  δx = nothing
  try
    δx = J \ (-F0)
  catch e
    if e isa LinearAlgebra.SingularException
      # Fallback to pseudo-inverse if J is singular/ill-conditioned
      δx = pinv(J) * (-F0)
    else
      rethrow(e)
    end
  end

  # Damping
  δx .*= damp

  # Map update back to Vr/Vi
  Vr_new = copy(Vr)
  Vi_new = copy(Vi)

  @inbounds for (idx, bus) in enumerate(non_slack)
    Vr_new[bus] += δx[idx]
    Vi_new[bus] += δx[(n-1)+idx]
  end

  # Keep slack voltage fixed
  Vr_new[slack_idx] = Vr[slack_idx]
  Vi_new[slack_idx] = Vi[slack_idx]

  return ComplexF64.(Vr_new, Vi_new)
end

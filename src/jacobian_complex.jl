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

# file: src/jacobian_complex.jl
# jacobian_complex.jl — Complex-State Newton-Raphson Power Flow Formulation
#
# This module implements a Newton-Raphson power flow solver using complex voltages
# in rectangular coordinates (Vr + jVi) as state variables, as an alternative to
# the conventional polar formulation (Vm, θ).
#
# Features:
# - Rectangular complex-state Newton-Raphson with PQ and PV bus handling
# - Wirtinger calculus-based Jacobian construction for complex power equations
# - Active-set Q-limit management for PV→PQ switching with optional re-enable
# - Both analytic and finite-difference Jacobian options
# - Hysteresis and cooldown mechanisms for robust PV bus management
# - Direct integration with Sparlectra.jl network data structures
#
# Mathematical Foundation:
# - State vector: x = [Vr(non-slack); Vi(non-slack)] ∈ ℝ^(2(n-1))
# - Complex power: S = V .* conj(Y * V) where Y is the bus admittance matrix  
# - PQ buses: ΔP = Re(S_calc - S_spec), ΔQ = Im(S_calc - S_spec)
# - PV buses: ΔP = Re(S_calc - S_spec), ΔV = |V| - V_set
# - Slack bus voltage is held constant throughout the iteration
#
# Key Functions:
# - run_complex_nr_rectangular_for_net!(): Main solver interface
# - runpf_rectangular!(): Convenience wrapper matching runpf!() signature
# - build_complex_jacobian(): Wirtinger-based Jacobian block construction
# - mismatch_rectangular(): Residual function for PQ/PV bus constraints
#
# Note:
# - The FD Jacobian is mathematically dense; sparse storage does not bring much benefit.
# - The analytic Jacobian currently uses a dense rectangular build, even though the
#   underlying structure is sparse (Ybus-like). A true sparse implementation would
#   require a dedicated builder similar to `calcJacobian(...; sparse=true)`.
#
# References:
# - Wirtinger calculus for complex derivatives

using LinearAlgebra
using SparseArrays
using Printf

"""
    build_complex_jacobian(Ybus, V)

Builds the 2n × 2n Wirtinger-type Jacobian blocks for the complex-state
Newton–Raphson formulation.

Given:
    I = Ybus * V
    S = V .* conj.(I)

We construct the blocks:
    J11 = ∂S/∂V
    J12 = ∂S/∂V*
    J21 = ∂conj(S)/∂V
    J22 = ∂conj(S)/∂V*

Returns:
    J11, J12, J21, J22  (all full matrices, not Diagonal)
"""
function build_complex_jacobian(Ybus, V)
  I = Ybus * V
  n = length(V)

  # J11 = diag(conj(I))
  J11 = Matrix(Diagonal(conj.(I)))

  # J12 = diag(V) * conj(Ybus)
  J12 = Matrix(Diagonal(V)) * conj.(Ybus)

  # J21 = conj(J12)
  J21 = conj.(J12)

  # J22 = conj(J11)
  J22 = conj.(J11)

  return J11, J12, J21, J22
end

"""
    complex_newton_step_rectangular(Ybus, V, S; slack_idx::Int, damp::Float64=1.0)

Perform one Newton-Raphson step in rectangular coordinates using analytic Jacobian.

# Arguments
- `Ybus::AbstractMatrix{ComplexF64}`: Bus admittance matrix (n×n)
- `V::Vector{ComplexF64}`: Current complex bus voltage vector (length n)
- `S::Vector{ComplexF64}`: Specified complex power injections P + jQ (length n)
- `slack_idx::Int`: Index of the slack bus (excluded from equations/variables)
- `damp::Float64=1.0`: Damping factor for Newton step (0 < damp ≤ 1)
- `bus_types::Vector{Symbol}`: Bus type classification (:PQ, :PV, :Slack)
- `Vset::Vector{Float64}`: Voltage magnitude setpoints for PV buses

# Returns
- `Vector{ComplexF64}`: Updated complex voltage vector V_new

# Details
This function performs one Newton-Raphson iteration using:
- State variables: [Vr(non-slack); Vi(non-slack)] - real and imaginary parts
- Mismatch function matching `mismatch_rectangular`:
  * PQ buses: ΔP_i, ΔQ_i (active and reactive power mismatches)
  * PV buses: ΔP_i, ΔV_i where ΔV_i = |V_i| - Vset[i]
- Analytic Jacobian derived from Wirtinger calculus
- Slack bus voltage remains fixed throughout the step

The function builds the analytic Jacobian via `build_rectangular_jacobian_pq_pv`,
solves the linear system J·δx = -F, applies damping, and updates voltages.

# Mathematical Foundation
Uses Wirtinger derivatives for complex power S = V .* conj(Y * V):
- ∂S/∂Vr and ∂S/∂Vi computed analytically
- PV constraint ΔV = |V| - Vset handled via chain rule
- Robust linear solver with pseudoinverse fallback for singular cases
"""

function complex_newton_step_rectangular(Ybus, V, S; slack_idx::Int, damp::Float64 = 1.0)
  n = length(V)

  # Currents and complex power
  I = Ybus * V
  S_calc = V .* conj.(I)
  ΔS = S_calc .- S

  # Real state representation
  Vr = real.(V)
  Vi = imag.(V)

  # Complex Jacobian Jc: n × 2n, mapping [dVr; dVi] → dS (complex)
  Jc = Matrix{ComplexF64}(undef, n, 2n)

  # Derivative wrt Vr: perturb dVr = e_k, dVi = 0
  for k = 1:n
    dVr = zeros(Float64, n)
    dVi = zeros(Float64, n)
    dVr[k] = 1.0
    dV = ComplexF64.(dVr, dVi)

    dI = Ybus * dV
    dS = dV .* conj.(I) .+ V .* conj.(dI)

    Jc[:, k] = dS
  end

  # Derivative wrt Vi: perturb dVi = e_k, dVr = 0
  for k = 1:n
    dVr = zeros(Float64, n)
    dVi = zeros(Float64, n)
    dVi[k] = 1.0
    dV = ComplexF64.(dVr, dVi)

    dI = Ybus * dV
    dS = dV .* conj.(I) .+ V .* conj.(dI)

    Jc[:, n+k] = dS
  end

  # Build real 2n × 2n Jacobian for F = [real(ΔS); imag(ΔS)]
  J = zeros(Float64, 2n, 2n)
  for k = 1:2n
    col            = Jc[:, k]
    J[1:n, k]      .= real.(col)
    J[(n+1):2n, k] .= imag.(col)
  end

  # Real mismatch vector F = [Re(ΔS); Im(ΔS)]
  F = vcat(real.(ΔS), imag.(ΔS))

  # --- Slack-Elimination: nur Nicht-Slack-Busse als Unbekannte/ Gleichungen ---

  # Bus-Indizes ohne Slack
  non_slack = collect(1:n)
  deleteat!(non_slack, slack_idx)

  # Zeilen/Gleichungen: P/Q für alle Nicht-Slack-Busse
  row_idx = vcat(non_slack, n .+ non_slack)  # erst P, dann Q

  # Spalten/Variablen: Vr/Vi für alle Nicht-Slack-Busse
  col_idx = vcat(non_slack, n .+ non_slack)

  Jred = J[row_idx, col_idx]
  Fred = F[row_idx]

  # Solve Jred * δx_red = -Fred (robust gegen Singularität)
  δx_red = nothing
  try
    δx_red = Jred \ (-Fred)
  catch e
    if e isa LinearAlgebra.SingularException
      δx_red = pinv(Jred) * (-Fred)
    else
      rethrow(e)
    end
  end

  # Map back to full state vector (slack has Δ=0)
  δx = zeros(Float64, 2n)
  δx[col_idx] .= δx_red

  # Damping
  δx .*= damp

  δVr = δx[1:n]
  δVi = δx[(n+1):2n]

  V_new = ComplexF64.(Vr .+ δVr, Vi .+ δVi)

  return V_new
end

function run_complex_nr_rectangular(Ybus, V0, S; slack_idx::Int = 1, maxiter::Int = 20, tol::Float64 = 1e-8, verbose::Bool = false, damp::Float64 = 1.0, bus_types::Vector{Symbol}, Vset::Vector{Float64}, use_fd::Bool = false, use_sparse::Bool = false)
  V = copy(V0)
  history = Float64[]

  for iter = 1:maxiter
    F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
    max_mis = maximum(abs.(F))
    push!(history, max_mis)

    if max_mis <= tol
      return V, true, iter, history
    end

    if use_fd
      V = complex_newton_step_rectangular_fd(Ybus, V, S; slack_idx = slack_idx, damp = damp, h = 1e-6, bus_types = bus_types, Vset = Vset)
    else
      V = complex_newton_step_rectangular(Ybus, V, S; slack_idx = slack_idx, damp = damp, bus_types = bus_types, Vset = Vset, use_sparse = use_sparse)
    end
  end

  return V, false, maxiter, history
end

"""
    update_net_voltages_from_complex!(net, V)

Update the bus voltage magnitudes and angles in the network from the
final complex voltages V (in per-unit).
"""
function update_net_voltages_from_complex!(net::Net, V::Vector{ComplexF64})
  nodes = net.nodeVec
  n = length(nodes)
  @assert length(V) == n

  for (k, node) in enumerate(nodes)
    Vk = V[k]
    vm = abs(Vk)
    va_rad = angle(Vk)
    va_deg = rad2deg(va_rad)
    node._vm_pu = vm
    node._va_deg = va_deg
  end
end

"""
    mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx) -> F::Vector{Float64}

Compute the real-valued mismatch vector F(V) for the rectangular
complex-state formulation with PQ and PV buses.

For each non-slack bus i:
- if bus_types[i] == :PQ:
      ΔP_i = Re(S_calc[i]) - Re(S_spec[i])
      ΔQ_i = Im(S_calc[i]) - Im(S_spec[i])

- if bus_types[i] == :PV:
      ΔP_i = Re(S_calc[i]) - Re(S_spec[i])
      ΔV_i = |V[i]| - Vset[i]

F is stacked as [ΔP_2, ΔQ/ΔV_2, ..., ΔP_n, ΔQ/ΔV_n] over all non-slack buses.
"""
function mismatch_rectangular(Ybus, V::Vector{ComplexF64}, S::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int)
  n = length(V)
  @assert length(S) == n
  @assert length(bus_types) == n
  @assert length(Vset) == n

  # Network-based injections for the current state
  I      = Ybus * V
  S_calc = V .* conj.(I)

  # F has 2*(n-1) entries: for each non-slack bus two residuals
  # PQ:  ΔP_i, ΔQ_i
  # PV:  ΔP_i, ΔV_i
  F = zeros(Float64, 2*(n-1))

  row = 1
  @inbounds for i = 1:n
    if i == slack_idx
      continue
    end

    S_ci = S_calc[i]
    S_si = S[i]

    if bus_types[i] == :PQ
      ΔP = real(S_ci) - real(S_si)
      ΔQ = imag(S_ci) - imag(S_si)

      F[row]   = ΔP
      F[row+1] = ΔQ

    elseif bus_types[i] == :PV
      ΔP = real(S_ci) - real(S_si)
      ΔV = abs(V[i]) - Vset[i]

      F[row]   = ΔP
      F[row+1] = ΔV

    else
      error("mismatch_rectangular: unsupported bus type $(bus_types[i]) at bus $i")
    end

    row += 2
  end

  return F
end

"""
    build_rectangular_jacobian_pq_pv_sparse(
        Ybus,
        V,
        bus_types,
        Vset,
        slack_idx,
    ) -> SparseMatrixCSC{Float64}

Builds the analytic rectangular Jacobian corresponding to `mismatch_rectangular`
using the sparsity pattern of `Ybus`.

State vector:
    x = [Vr(non-slack); Vi(non-slack)] ∈ ℝ^(2(n-1))

Residual F(V):
    - PQ buses: ΔP_i, ΔQ_i
    - PV buses: ΔP_i, ΔV_i
    - Slack bus: no equations

Jacobian entries are derived from
    S_i(V) = V_i * conj( (Ybus * V)_i )

Wirtinger-based identities:
    ∂S/∂V   = diag(conj(I)) + diag(V) * conj(Ybus)
    ∂S/∂V*  = diag(V) * conj(Ybus)

Chain rule to rectangular:
    ∂S/∂Vr = ∂S/∂V + ∂S/∂V*
    ∂S/∂Vi = j(∂S/∂V - ∂S/∂V*)

With ΔP_i = Re(ΔS_i), ΔQ_i = Im(ΔS_i), ΔV_i = |V_i| - Vset[i].

Returns:
    J :: SparseMatrixCSC{Float64} with size (2(n-1)) × (2(n-1)).
"""
function build_rectangular_jacobian_pq_pv_sparse(Ybus::SparseMatrixCSC{ComplexF64}, V::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int)
  n = length(V)
  @assert length(bus_types) == n
  @assert length(Vset) == n

  # Network currents
  I = Ybus * V

  # Non-slack bus indices and mapping bus -> position in state vector
  non_slack = collect(1:n)
  deleteat!(non_slack, slack_idx)

  nvar = 2 * (n - 1)   # [Vr(non-slack); Vi(non-slack)]
  m    = nvar          # F has 2 equations per non-slack bus

  # bus index -> position in non_slack (1..n-1), 0 if slack
  pos_non_slack = zeros(Int, n)
  for (k, b) in enumerate(non_slack)
    pos_non_slack[b] = k
  end

  # Row blocks: for each non-slack bus i
  #   row_block[i]   = index of ΔP row for bus i
  #   row_block[i]+1 = index of ΔQ / ΔV row for bus i
  row_block = zeros(Int, n)
  row = 0
  for i = 1:n
    if i == slack_idx
      continue
    end
    row += 2
    row_block[i] = row - 1  # ΔP row
  end

  # Triplet storage
  Iidx = Int[]
  Jidx = Int[]
  Vals = Float64[]
  # Rough capacity hint (purely heuristic)
  sizehint!(Iidx, 16 * nnz(Ybus))
  sizehint!(Jidx, 16 * nnz(Ybus))
  sizehint!(Vals, 16 * nnz(Ybus))

  # Helpers for sparse access
  rv    = rowvals(Ybus)
  nzval = nonzeros(Ybus)

  # --- 1) Contributions from complex power equations (ΔP, ΔQ) ---------------
  #
  # For bus i, bus j:
  #   ∂S_i/∂Vr_j = conj(I_i) * δ_ij + V_i * conj(Y_ij)
  #   ∂S_i/∂Vi_j = j*(conj(I_i) * δ_ij - V_i * conj(Y_ij))
  #
  # Then:
  #   ∂P_i/∂Vr_j = Re(∂S_i/∂Vr_j)
  #   ∂P_i/∂Vi_j = Re(∂S_i/∂Vi_j)
  #   ∂Q_i/∂Vr_j = Im(∂S_i/∂Vr_j)
  #   ∂Q_i/∂Vi_j = Im(∂S_i/∂Vi_j)
  #
  # For PV buses, only the ΔP row uses these derivatives; the second row is ΔV.

  for j = 1:n
    col_pos = pos_non_slack[j]
    if col_pos == 0
      # Slack bus column -> no state variable
      continue
    end

    colVr = col_pos
    colVi = (n - 1) + col_pos

    for ptr in nzrange(Ybus, j)
      i = rv[ptr]

      # Slack bus has no equations
      if i == slack_idx
        continue
      end

      rb = row_block[i]
      if rb == 0
        continue
      end

      rowP = rb          # ΔP row for bus i
      rowQ = rb + 1      # ΔQ/ΔV row for bus i

      Yij = nzval[ptr]

      # Base contributions from V_i * conj(Y_ij)
      S_dVr = V[i] * conj(Yij)
      S_dVi = im * (-V[i] * conj(Yij))

      # Diagonal term from conj(I_i) * δ_ij
      if i == j
        Ii = I[i]
        S_dVr += conj(Ii)
        S_dVi += im * conj(Ii)
      end

      # Real / imaginary parts
      ∂P_Vr = real(S_dVr)
      ∂P_Vi = real(S_dVi)
      ∂Q_Vr = imag(S_dVr)
      ∂Q_Vi = imag(S_dVi)

      # First equation: always ΔP_i for PQ and PV
      if abs(∂P_Vr) > 0.0
        push!(Iidx, rowP);
        push!(Jidx, colVr);
        push!(Vals, ∂P_Vr)
      end
      if abs(∂P_Vi) > 0.0
        push!(Iidx, rowP);
        push!(Jidx, colVi);
        push!(Vals, ∂P_Vi)
      end

      # Second equation:
      #   - PQ: ΔQ_i -> uses Q derivatives
      #   - PV: ΔV_i -> no contribution from S, handled separately
      bt = bus_types[i]
      if bt == :PQ
        if abs(∂Q_Vr) > 0.0
          push!(Iidx, rowQ);
          push!(Jidx, colVr);
          push!(Vals, ∂Q_Vr)
        end
        if abs(∂Q_Vi) > 0.0
          push!(Iidx, rowQ);
          push!(Jidx, colVi);
          push!(Vals, ∂Q_Vi)
        end
      elseif bt == :PV
        # nothing here, ΔV row added below
      else
        error("build_rectangular_jacobian_pq_pv_sparse: unsupported bus type $(bt) at bus $i")
      end
    end
  end

  # --- 2) Contributions for ΔV_i = |V_i| - Vset[i] on PV buses --------------
  #
  # Only depends on local Vr_i, Vi_i:
  #   |V_i| = sqrt(Vr_i^2 + Vi_i^2)
  #   ∂|V_i|/∂Vr_i = Vr_i / |V_i|
  #   ∂|V_i|/∂Vi_i = Vi_i / |V_i|

  for i = 1:n
    if i == slack_idx || bus_types[i] != :PV
      continue
    end

    rb   = row_block[i]
    rowV = rb + 1  # second row for that bus

    pos = pos_non_slack[i]
    if pos == 0
      continue
    end

    vm = abs(V[i])
    if vm == 0.0
      continue
    end

    dVr = real(V[i]) / vm
    dVi = imag(V[i]) / vm

    colVr = pos
    colVi = (n - 1) + pos

    if abs(dVr) > 0.0
      push!(Iidx, rowV);
      push!(Jidx, colVr);
      push!(Vals, dVr)
    end
    if abs(dVi) > 0.0
      push!(Iidx, rowV);
      push!(Jidx, colVi);
      push!(Vals, dVi)
    end
  end

  return sparse(Iidx, Jidx, Vals, m, nvar)
end

"""
    build_rectangular_jacobian_pq_pv_dense(
        Ybus, V, bus_types, Vset, slack_idx
    ) -> J::Matrix{Float64}

Build the analytic rectangular Jacobian for the mismatch vector `F(V)`
defined in `mismatch_rectangular`.

- State vector: x = [Vr(non-slack); Vi(non-slack)]
- Rows: for each non-slack bus i
    * PQ: [ΔP_i; ΔQ_i]
    * PV: [ΔP_i; ΔV_i]  with ΔV_i = |V_i| - Vset[i]

`bus_types` and `Vset` must be consistent with `mismatch_rectangular`.
"""
function build_rectangular_jacobian_pq_pv_dense(Ybus, V::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int)
  n = length(V)
  @assert length(bus_types) == n
  @assert length(Vset) == n

  # --- 1) Wirtinger blocks for S(V) = V .* conj(Ybus * V)
  J11, J12, J21, J22 = build_complex_jacobian(Ybus, V)

  # --- 2) Full 2n×2n rectangular J for ΔP/ΔQ wrt Vr/Vi (all buses)
  # Rows: [ΔP_1..ΔP_n; ΔQ_1..ΔQ_n]
  # Cols: [Vr_1..Vr_n; Vi_1..Vi_n]
  Jrect_full = zeros(Float64, 2n, 2n)

  @inbounds for j = 1:n
    col_sum  = J11[:, j] .+ J12[:, j]  # corresponds to dS/dVr_j
    col_diff = J11[:, j] .- J12[:, j]  # used for dS/dVi_j

    # dP/dVr_j, dQ/dVr_j
    @views Jrect_full[1:n, j]      .= real.(col_sum)
    @views Jrect_full[(n+1):2n, j] .= imag.(col_sum)

    # dP/dVi_j, dQ/dVi_j
    @views Jrect_full[1:n, n+j]      .= -imag.(col_diff)
    @views Jrect_full[(n+1):2n, n+j] .= real.(col_diff)
  end

  # --- 3) Reduce to non-slack variables and rows matching mismatch_rectangular

  non_slack = collect(1:n)
  deleteat!(non_slack, slack_idx)

  nvar = 2 * (n - 1)
  m    = 2 * (n - 1)
  @assert nvar == m

  # Column indices in the full rectangular J that correspond to
  # [Vr(non_slack); Vi(non_slack)]
  col_idx_full = vcat(non_slack, n .+ non_slack)

  J = zeros(Float64, m, nvar)

  row = 1
  @inbounds for i = 1:n
    if i == slack_idx
      continue
    end

    # First row for this bus: ΔP_i
    rowP_full = i                  # P row index in full J
    @views J[row, :] .= Jrect_full[rowP_full, col_idx_full]

    # Second row: ΔQ_i (PQ) or ΔV_i (PV)
    if bus_types[i] == :PQ
      rowQ_full = n + i          # Q row index in full J
      @views J[row+1, :] .= Jrect_full[rowQ_full, col_idx_full]

    elseif bus_types[i] == :PV
      # ΔV_i = |V_i| - Vset[i]
      J[row+1, :] .= 0.0

      pos = findfirst(==(i), non_slack)
      if pos !== nothing
        vm = abs(V[i])
        if vm > 0.0
          dVr = real(V[i]) / vm
          dVi = imag(V[i]) / vm

          # Columns in reduced J:
          #   Vr_i -> index pos
          #   Vi_i -> index (n-1) + pos
          J[row+1, pos]       = dVr
          J[row+1, (n-1)+pos] = dVi
        end
      end
    else
      error("build_rectangular_jacobian_pq_pv: unsupported bus type $(bus_types[i]) at bus $i")
    end

    row += 2
  end

  return J
end

"""
    build_rectangular_jacobian_pq_pv(
        Ybus,
        V,
        bus_types,
        Vset,
        slack_idx;
        use_sparse::Bool = false,
    )

Dispatches to either the dense or sparse rectangular Jacobian builder matching
`mismatch_rectangular`.

- If `use_sparse == true` and `Ybus` is a `SparseMatrixCSC{ComplexF64}`, the
  sparse builder is used.
- Otherwise, the dense builder is used.
"""
function build_rectangular_jacobian_pq_pv(Ybus, V::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int; use_sparse::Bool = false)
  if use_sparse && Ybus isa SparseMatrixCSC{ComplexF64}
    return build_rectangular_jacobian_pq_pv_sparse(Ybus, V, bus_types, Vset, slack_idx)
  else
    return build_rectangular_jacobian_pq_pv_dense(Ybus, V, bus_types, Vset, slack_idx)
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
        use_sparse=false,
    )

Performs one Newton–Raphson step in rectangular coordinates using the analytic
Jacobian that matches `mismatch_rectangular`.

- State: x = [Vr(non-slack); Vi(non-slack)]
- Residual: F(x) = mismatch_rectangular(...)
"""
function complex_newton_step_rectangular(Ybus, V::Vector{ComplexF64}, S::Vector{ComplexF64}; slack_idx::Int, damp::Float64 = 1.0, bus_types::Vector{Symbol}, Vset::Vector{Float64}, use_sparse::Bool = false)
  n = length(V)
  @assert length(S) == n
  @assert length(bus_types) == n
  @assert length(Vset) == n

  # Non-slack indices
  non_slack = collect(1:n)
  deleteat!(non_slack, slack_idx)

  # Residual matching the FD variant
  F0 = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
  m = length(F0)
  nvar = 2 * (n - 1)
  @assert m == nvar "complex_newton_step_rectangular: mismatch and state dimension differ"

  # Analytic Jacobian (dense or sparse)
  J = build_rectangular_jacobian_pq_pv(Ybus, V, bus_types, Vset, slack_idx; use_sparse = use_sparse)

  # Solve J * δx = -F
  δx = nothing
  try
    δx = J \ (-F0)
  catch e
    if e isa LinearAlgebra.SingularException
      δx = pinv(Matrix(J)) * (-F0)
    else
      rethrow(e)
    end
  end

  # Damping
  δx .*= damp

  # Map update back to Vr/Vi
  Vr = real.(V)
  Vi = imag.(V)

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

"""
    run_complex_nr_rectangular_for_net!(net; maxiter=20, tol=1e-8, damp=0.2, verbose=0, use_fd=false)

Run a complex-state Newton-Raphson power flow in rectangular coordinates on a Sparlectra network.

# Arguments
- `net::Net`: Network object containing bus, branch, and generation data
- `maxiter::Int=20`: Maximum number of Newton-Raphson iterations
- `tol::Float64=1e-8`: Convergence tolerance for maximum mismatch
- `damp::Float64=0.2`: Damping factor for Newton step (0 < damp ≤ 1)
- `verbose::Int=0`: Verbosity level (0=quiet, 1=basic info, 2=detailed)
- `use_fd::Bool=false`: Use finite-difference Jacobian instead of analytic

# Returns
- `Tuple{Int, Int}`: (iterations_used, error_code)
  - error_code: 0=converged, 1=max_iterations_reached

# Details
This function implements a complete power flow solver using complex voltages in rectangular 
coordinates (Vr + jVi) as state variables, providing an alternative to conventional 
polar formulations.

## Mathematical Foundation
- **State vector**: x = [Vr(non-slack); Vi(non-slack)] ∈ ℝ^(2(n-1))
- **Complex power**: S = V .* conj(Y * V) where Y is the bus admittance matrix
- **PQ buses**: ΔP = Re(S_calc - S_spec), ΔQ = Im(S_calc - S_spec)
- **PV buses**: ΔP = Re(S_calc - S_spec), ΔV = |V| - V_set
- **Slack bus**: Voltage held constant throughout iterations

## Active Set Q-Limit Management
- **PV→PQ switching**: When reactive power demand violates generator Q-limits
- **Optional PQ→PV re-enable**: With hysteresis band and cooldown mechanisms
- **Robust handling**: Guards against inappropriate switching of non-generator buses

## Algorithm Steps
1. **Initialization**: Extract voltages, build Y-bus, classify bus types
2. **Power specification**: Build S = P + jQ from network loads/generation
3. **Iterative solution**: Newton-Raphson with mismatch function for PQ/PV constraints
4. **Q-limit enforcement**: Active-set management during iterations
5. **Result update**: Write final voltages and computed powers back to network

## Network Integration
- **Input**: Uses `net.nodeVec` for bus data, `net.baseMVA` for per-unit conversion
- **Output**: Updates `node._vm_pu`, `node._va_deg`, `node._pƩGen`, `node._qƩGen`
- **Q-limits**: Integrates with `net.qLimitEvents` and `net.qLimitLog` for tracking
- **Compatibility**: Maintains same interface as `runpf!()` for easy substitution

# See Also
- `runpf_rectangular!()`: Convenience wrapper matching `runpf!()` signature
- `mismatch_rectangular()`: Core mismatch function for PQ/PV constraints
- `build_rectangular_jacobian_pq_pv()`: Analytic Jacobian construction
"""

function run_complex_nr_rectangular_for_net!(net::Net; maxiter::Int = 20, tol::Float64 = 1e-8, damp::Float64 = 0.2, verbose::Int = 0, use_fd::Bool = false, opt_sparse::Bool = true, opt_flatstart::Bool = net.flatstart)
  if verbose > 1
    @info "Running complex rectangular NR power flow... use_fd=$use_fd, opt_sparse=$opt_sparse"
  end

  nodes = net.nodeVec
  n     = length(nodes)
  Sbase = net.baseMVA
  if !opt_sparse
    sparse = n > 1000
  else
    sparse = opt_sparse
  end
  Ybus = createYBUS(net = net, sparse = sparse, printYBUS = (verbose > 1))

  # 1) Initial complex voltages V0 and slack index
  V0, slack_idx = initialVrect(net; flatstart = opt_flatstart)

  # 2) Specified complex power injections S (p.u.), sign convention wie im Polar-NR
  S = buildComplexSVec(net)

  # 3) Bus types and PV setpoints from Node data
  bus_types = Vector{Symbol}(undef, n)
  Vset      = Vector{Float64}(undef, n)

  @inbounds for (k, node) in enumerate(nodes)
    BusType = getNodeType(node)
    if BusType == Slack
      bus_types[k] = :Slack
    elseif BusType == PV
      bus_types[k] = :PV
    elseif BusType == PQ
      bus_types[k] = :PQ
    else
      error("run_complex_nr_rectangular_for_net!: unsupported bus type at bus $k")
    end

    Vset[k] = isnothing(node._vm_pu) ? 1.0 : node._vm_pu
  end

  # 4) Q-limit data 
  qmin_pu, qmax_pu = getQLimits_pu(net)
  resetQLimitLog!(net)

  cooldown_iters = hasfield(typeof(net), :cooldown_iters) ? net.cooldown_iters : 0
  q_hyst_pu      = hasfield(typeof(net), :q_hyst_pu) ? net.q_hyst_pu : 0.0

  # Original PV buses at start (Guard for PQ->PV re-enable)
  pv_orig = [k for k = 1:n if bus_types[k] == :PV]

  # 5) NR-Loop
  V         = copy(V0)
  history   = Float64[]
  converged = false
  iters     = 0

  if verbose > 1
    @info "Starting rectangular complex NR power flow..."
    @info "Initial complex voltages V0:" V0
    @info "Slack bus index:" slack_idx
    @info "maxiter = $maxiter, tol = $tol, damp = $damp"
  end

  for it = 1:maxiter
    iters = it

    # Mismatch with current bus_types & S
    F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
    max_mis = maximum(abs.(F))
    push!(history, max_mis)

    (verbose > 1) && @debug "Rectangular NR iteration" iter=it max_mismatch=max_mis

    if max_mis <= tol
      converged = true
      break
    end

    # --- Q-Limit Active Set: PV -> PQ, optional PQ -> PV -------------------
    changed   = false
    reenabled = false

    if it > 1
      I      = Ybus * V
      S_calc = V .* conj.(I)

      # 4a) PV -> PQ
      @inbounds for k = 1:n
        if k == slack_idx
          continue
        end
        if bus_types[k] != :PV
          continue
        end

        qreq   = imag(S_calc[k])    # current Q requirement (p.u.)
        busIdx = k                  # here busIdx = position in nodeVec

        if qreq > qmax_pu[k]
          # Bus becomes PQ, Q clamped to Qmax
          bus_types[k] = :PQ
          S[k] = complex(real(S[k]), net.qmax_pu[k])
          changed = true

          net.qLimitEvents[busIdx] = :max
          logQLimitHit!(net, it, busIdx, :max)
          (verbose>2) && @printf "PV->PQ Bus %d: Q=%.4f > Qmax=%.4f (it=%d)\n" busIdx qreq qmax_pu[k] it

        elseif qreq < qmin_pu[k]
          bus_types[k] = :PQ
          S[k] = complex(real(S[k]), net.qmin_pu[k])
          changed = true

          net.qLimitEvents[busIdx] = :min
          logQLimitHit!(net, it, busIdx, :min)
          (verbose>2) && @printf "PV->PQ Bus %d: Q=%.4f < Qmin=%.4f (it=%d)\n" busIdx qreq qmin_pu[k] it
        end
      end

      if changed
        F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
        max_mis = maximum(abs.(F))
      end

      # 4b) Optional PQ -> PV Re-enable (Hysterese + Cooldown)            
      if q_hyst_pu > 0.0 || cooldown_iters > 0
        @inbounds for k = 1:n
          # Guards:
          # (1) currently PQ
          # (2) was originally PV
          # (3) has a Q-limit event in this calculation
          if bus_types[k] != :PQ
            continue
          end
          if !(k in pv_orig)
            continue
          end
          if !haskey(net.qLimitEvents, k)
            continue
          end

          qreq = imag(S_calc[k])
          lo = net.qmin_pu[k] + q_hyst_pu
          hi = net.qmax_pu[k] - q_hyst_pu
          ready = (qreq > lo) && (qreq < hi)

          if ready && (cooldown_iters > 0)
            last_it = lastQLimitIter(net, k)
            if !isnothing(last_it) && (it - last_it) < cooldown_iters
              ready = false
            end
          end

          if ready
            bus_types[k] = :PV
            delete!(net.qLimitEvents, k)
            reenabled = true
            (verbose>0) && @printf "PQ->PV Bus %d: Q=%.4f in [%.4f, %.4f] (cooldown=%d)\n" k qreq lo hi cooldown_iters
          end
        end

        if reenabled
          F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
          max_mis = maximum(abs.(F))
        end
      end
    end

    # --- Newton-Stepp (FD oder analytisch) -----------------------------
    if use_fd
      V = complex_newton_step_rectangular_fd(Ybus, V, S; slack_idx = slack_idx, damp = damp, h = 1e-6, bus_types = bus_types, Vset = Vset)
    else
      V = complex_newton_step_rectangular(Ybus, V, S; slack_idx = slack_idx, damp = damp, bus_types = bus_types, Vset = Vset, use_sparse = sparse)
    end

    # Keep slack voltage fixed (safety belt)
    V[slack_idx] = V0[slack_idx]
  end

  # 6) Update voltages back to network
  update_net_voltages_from_complex!(net, V)

  # 7) Compute bus injections from final voltages
  Ibus     = Ybus * V
  Sbus_pu  = V .* conj.(Ibus)
  Sbus_MVA = Sbus_pu .* Sbase

  @debug "Final Voltages Mag = " [abs.(V)...]
  @debug "Final Voltages Ang = " [angle.(V) .* (180.0 / π)...]

  isoNodes = net.isoNodes

  @inbounds for (k, node) in enumerate(nodes)
    # Skip isolated nodes (if any)
    if k in isoNodes
      continue
    end

    Sbus      = Sbus_MVA[k]
    Pbus_MW   = real(Sbus)
    Qbus_MVar = imag(Sbus)

    # Slack bus: always write P and Q generation from the solved state
    if node._nodeType == Sparlectra.Slack
      node._pƩGen = Pbus_MW
      node._qƩGen = Qbus_MVar

    elseif node._nodeType == Sparlectra.PV
      node._qƩGen = Qbus_MVar
      if haskey(net.qLimitEvents, k)
        @debug "Bus $(k) hit Q limit; set as PQ bus (rectangular solver)."
      end
    end
    # PQ buses / pure loads: do not touch _pƩLoad / _pƩGen here.
    # The original load/generation specification remains intact.
  end

  # 8) Update total bus power (sum of complex injections in p.u.)
  p = (sum(real.(Sbus_pu))) * Sbase
  q = (sum(imag.(Sbus_pu))) * Sbase

  if verbose > 1
    @info "Set total bus power to p = $p MW and q = $q MVar"
  end

  setTotalBusPower!(net = net, p = p, q = q)

  return iters, converged ? 0 : 1
end

"""
    runpf_rectangular!(net, maxIte, tolerance=1e-6, verbose=0)

Runs a rectangular complex-state Newton–Raphson power flow on `net::Net`.

Returns:
    (iterations::Int, status::Int)
where `status == 0` indicates convergence.
"""
function runpf_rectangular!(net::Net, maxIte::Int, tolerance::Float64 = 1e-6, verbose::Int = 0; opt_fd::Bool = false, opt_sparse::Bool = true, damp = 1.0, opt_flatstart::Bool = net.flatstart)
  iters, erg = run_complex_nr_rectangular_for_net!(net; maxiter = maxIte, tol = tolerance, damp = damp, verbose = verbose, use_fd = opt_fd, opt_sparse = opt_sparse, opt_flatstart = opt_flatstart)
  return iters, erg
end

"""
    runpf!(net, maxIte, tolerance=1e-6, verbose=0; method=:polar_full)

Unified AC power flow interface.

Arguments:
- `net::Net`: network
- `maxIte::Int`: maximum iterations
- `tolerance::Float64`: mismatch tolerance
- `verbose::Int`: verbosity level
- `method::Symbol`: `:polar_full` (default) or `:rectangular`

Returns:
    (iterations::Int, status::Int)

where `status == 0` indicates convergence.
"""
function runpf!(net::Net, maxIte::Int, tolerance::Float64 = 1e-6, verbose::Int = 0; method::Symbol = :rectangular, opt_fd::Bool = false, opt_sparse::Bool = true, opt_flatstart::Bool = net.flatstart, damp = 1.0)
  #@info "Running AC Power Flow using method: $(method)"
  if method === :polar_full
    return runpf_full!(net, maxIte, tolerance, verbose; opt_sparse    = opt_sparse, opt_flatstart = opt_flatstart)
  elseif method === :rectangular
    return runpf_rectangular!(net, maxIte, tolerance, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, damp = damp, opt_flatstart = opt_flatstart)
  elseif method === :classic
    return runpf_classic!(net, maxIte, tolerance, verbose, opt_sparse, opt_flatstart)
  else
    error("runpf!: unknown method $(method). Use :polar_full or :rectangular.")
  end
end

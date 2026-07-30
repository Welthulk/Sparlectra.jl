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
# Rectangular power-flow analytic Jacobian builders.
#
# This file is included inside module Sparlectra. Do not add a module wrapper here.
# Date: 29.5.2026
# file: src/powerflow_rectangular/rectangular_jacobian_builders.jl  

"""
    build_rectangular_jacobian_pq_pv_sparse(
        Ybus::SparseMatrixCSC{ComplexF64},
        V::Vector{ComplexF64},
        bus_types::Vector{Symbol},
        Vset::Vector{Float64},
        slack_idx::Int;
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
# Emitter receiving each Jacobian entry during triplet assembly (build path).
struct _JacobianTripletEmitter
  Iidx::Vector{Int}
  Jidx::Vector{Int}
  Vals::Vector{Float64}
end

@inline function (e::_JacobianTripletEmitter)(row::Int, col::Int, val::Float64)
  push!(e.Iidx, row)
  push!(e.Jidx, col)
  push!(e.Vals, val)
  return nothing
end

# Emitter replaying the recorded emit order into an existing CSC's nzval
# (in-place refresh path). `overflow` flags an emit count beyond the recorded
# position map; the caller then falls back to a structural rebuild.
mutable struct _JacobianPosMapEmitter
  nz::Vector{Float64}
  pos_map::Vector{Int}
  k::Int
  overflow::Bool
end

@inline function (e::_JacobianPosMapEmitter)(row::Int, col::Int, val::Float64)
  k = e.k + 1
  e.k = k
  if k <= length(e.pos_map)
    @inbounds e.nz[e.pos_map[k]] += val
  else
    e.overflow = true
  end
  return nothing
end

function build_rectangular_jacobian_pq_pv_sparse(
  Ybus::SparseMatrixCSC{ComplexF64},
  V::Vector{ComplexF64},
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  slack_idx::Int;
  dPinj_dVm::Vector{Float64} = zeros(Float64, length(V)),
  dQinj_dVm::Vector{Float64} = zeros(Float64, length(V)),
  vm_eps::Float64 = 1e-9,
  structural_pattern::Bool = false,
  assembly::Union{Nothing,RectangularJacobianAssembly} = nothing,
  dslack::Union{Nothing,DistributedSlackState} = nothing,
)
  # vm_eps avoids unstable derivatives when |V| is very close to zero.
  # structural_pattern=true stores entries even when their current value is
  # exactly zero, so the sparsity pattern depends only on the Ybus structure
  # and the bus types — not on the iterate. Factorization-reuse backends
  # (KLU) require this invariance; the default path keeps dropping numeric
  # zeros to stay bit-for-bit identical to the historical behavior.
  # `assembly` (requires structural_pattern) reuses triplet buffers and, once
  # a pattern is recorded, refreshes the retained CSC's nzval in place.
  n = length(V)
  @assert length(bus_types) == n
  @assert length(Vset) == n
  @assert length(dPinj_dVm) == n
  @assert length(dQinj_dVm) == n
  @assert assembly === nothing || structural_pattern "assembly reuse requires the structural sparsity pattern"

  non_slack = non_slack_indices(n, slack_idx)
  pos_non_slack = build_pos_map(non_slack, n)
  # Distributed slack appends one state (lambda_P) and one equation (the REF
  # bus's P residual) — see mismatch_rectangular. The participant set is fixed
  # per solve, so the extended pattern is as iteration-stable as the base one.
  dsl = dslack !== nothing
  nvar = 2 * (n - 1) + (dsl ? 1 : 0)   # [Vr(non-slack); Vi(non-slack); lambda?]
  m = nvar          # F has 2 equations per non-slack bus (+ the REF-P row)

  # Row blocks follow mismatch_rectangular layout for non-slack buses.
  #   row_block[i]   = index of ΔP row for bus i
  #   row_block[i]+1 = index of ΔQ / ΔV row for bus i
  row_block = zeros(Int, n)
  row = 0
  for i = 1:n
    if i == slack_idx
      continue
    end
    row += 2
    row_block[i] = row - 1  # ΔP row for bus i
  end

  if assembly !== nothing && assembly.valid && assembly.J !== nothing && size(assembly.J, 1) == m
    # In-place refresh: same structural pattern as the recorded build, so the
    # emit order is identical and each entry lands at its recorded nzval slot.
    Jexist = assembly.J
    resize!(assembly.Ibuf, n)
    mul!(assembly.Ibuf, Ybus, V)
    nz = nonzeros(Jexist)
    fill!(nz, 0.0)
    emitter = _JacobianPosMapEmitter(nz, assembly.pos_map, 0, false)
    _assemble_rectangular_jacobian!(emitter, Ybus, V, bus_types, slack_idx, assembly.Ibuf, row_block, pos_non_slack, non_slack, dPinj_dVm, dQinj_dVm, vm_eps, true, n, dslack)
    if !emitter.overflow && emitter.k == length(assembly.pos_map)
      return Jexist
    end
    # Emit count drifted (should be unreachable): rebuild structurally.
    assembly.valid = false
  end

  I = Ybus * V

  # Triplet storage: fresh vectors, or the reused assembly buffers.
  if assembly === nothing
    Iidx = Int[]
    Jidx = Int[]
    Vals = Float64[]
  else
    Iidx = assembly.Iidx
    Jidx = assembly.Jidx
    Vals = assembly.Vals
    empty!(Iidx)
    empty!(Jidx)
    empty!(Vals)
  end
  # Structural upper bound: 4 entries per Ybus nonzero (P/Q × Vr/Vi) plus the
  # PV ΔV rows and injection chain terms (≤ 6 per bus).
  cap = 4 * nnz(Ybus) + 6 * n
  sizehint!(Iidx, cap)
  sizehint!(Jidx, cap)
  sizehint!(Vals, cap)

  emit = _JacobianTripletEmitter(Iidx, Jidx, Vals)
  _assemble_rectangular_jacobian!(emit, Ybus, V, bus_types, slack_idx, I, row_block, pos_non_slack, non_slack, dPinj_dVm, dQinj_dVm, vm_eps, structural_pattern, n, dslack)

  J = sparse(Iidx, Jidx, Vals, m, nvar)

  if assembly !== nothing
    # Record the emit-order → nzval position map for in-place refreshes.
    pos_map = assembly.pos_map
    resize!(pos_map, length(Iidx))
    colptr = J.colptr
    rowval = J.rowval
    @inbounds for t in eachindex(Iidx)
      c = Jidx[t]
      lo = colptr[c]
      hi = colptr[c+1] - 1
      pos_map[t] = searchsortedfirst(view(rowval, lo:hi), Iidx[t]) + lo - 1
    end
    assembly.J = J
    assembly.valid = true
  end

  return J
end

# Shared emit-order-deterministic assembly of all Jacobian entries. Both the
# triplet build and the in-place nzval refresh replay exactly this sequence;
# any change to loop order or emit conditions must keep the two paths in sync
# by construction (they run the same code).
function _assemble_rectangular_jacobian!(
  emit::E,
  Ybus::SparseMatrixCSC{ComplexF64},
  V::Vector{ComplexF64},
  bus_types::Vector{Symbol},
  slack_idx::Int,
  I::AbstractVector{ComplexF64},
  row_block::Vector{Int},
  pos_non_slack::Vector{Int},
  non_slack::Vector{Int},
  dPinj_dVm::Vector{Float64},
  dQinj_dVm::Vector{Float64},
  vm_eps::Float64,
  structural_pattern::Bool,
  n::Int,
  dslack::Union{Nothing,DistributedSlackState} = nothing,
) where {E}
  # Distributed slack: the appended REF-P row and the lambda column live at
  # index 2(n-1)+1. Their emits are woven into the SAME deterministic sequence
  # as everything else (REF-row network entries inside the main column loop,
  # constant alpha entries in a fixed block at the end), so the in-place
  # pos_map replay stays valid — dslack is constant for one solve, and any
  # toggle rebuilds the assembly like an active-set switch does.
  dsl = dslack !== nothing
  row_ref = 2 * (n - 1) + 1
  col_lambda = 2 * (n - 1) + 1
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

      # Slack bus has no equations — except with distributed slack, where its
      # P residual is the appended last row.
      if i == slack_idx
        if dsl
          Yij = nzval[ptr]
          S_dVr = V[i] * conj(Yij)
          S_dVi = im * (-V[i] * conj(Yij))
          # j == slack cannot occur here: the slack column carries no state
          # variable and was skipped above (col_pos == 0).
          emit(row_ref, colVr, real(S_dVr))
          emit(row_ref, colVi, real(S_dVi))
        end
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
      if structural_pattern || abs(∂P_Vr) > 0.0
        emit(rowP, colVr, ∂P_Vr)
      end
      if structural_pattern || abs(∂P_Vi) > 0.0
        emit(rowP, colVi, ∂P_Vi)
      end

      # Second equation:
      #   - PQ: ΔQ_i -> uses Q derivatives
      #   - PV: ΔV_i -> no contribution from S, handled separately
      bt = bus_types[i]
      if bt == :PQ
        if structural_pattern || abs(∂Q_Vr) > 0.0
          emit(rowQ, colVr, ∂Q_Vr)
        end
        if structural_pattern || abs(∂Q_Vi) > 0.0
          emit(rowQ, colVi, ∂Q_Vi)
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
    if vm == 0.0 && !structural_pattern
      continue
    end
    vm_div = vm > 0.0 ? vm : vm_eps

    dVr = real(V[i]) / vm_div
    dVi = imag(V[i]) / vm_div

    colVr = pos
    colVi = (n - 1) + pos

    if structural_pattern || abs(dVr) > 0.0
      emit(rowV, colVr, dVr)
    end
    if structural_pattern || abs(dVi) > 0.0
      emit(rowV, colVi, dVi)
    end
  end

  # Local chain-rule terms for voltage-dependent specified injections.
  # For ΔP = Pcalc - Pspec(|V|): subtract dPspec/d|V| * d|V|/dVr and d|V|/dVi.
  # For PQ second row ΔQ = Qcalc - Qspec(|V|): analogous subtraction.
  for i in non_slack
    pos = pos_non_slack[i]
    rb = row_block[i]
    pos == 0 && continue
    rb == 0 && continue

    vm = abs(V[i])
    vm_safe = vm > vm_eps ? vm : vm_eps
    dvm_dvr = real(V[i]) / vm_safe
    dvm_dvi = imag(V[i]) / vm_safe
    colVr = pos
    colVi = (n - 1) + pos

    # These duplicate triplets merge additively into the diagonal-column
    # entries already stored by the S-block, so pushing them unconditionally
    # in structural mode never introduces iterate-dependent pattern entries.
    dP = dPinj_dVm[i]
    if structural_pattern || dP != 0.0
      emit(rb, colVr, -dP * dvm_dvr)
      emit(rb, colVi, -dP * dvm_dvi)
    end

    if bus_types[i] == :PQ
      dQ = dQinj_dVm[i]
      if structural_pattern || dQ != 0.0
        emit(rb + 1, colVr, -dQ * dvm_dvr)
        emit(rb + 1, colVi, -dQ * dvm_dvi)
      end
    end
  end

  # --- 4) Distributed slack: lambda column (constant -alpha at participant P
  # rows and the REF row). The participant set is fixed per solve, so these
  # entries are pattern-stable regardless of structural_pattern.
  if dsl
    alpha = dslack.alpha
    for i = 1:n
      i == slack_idx && continue
      a = alpha[i]
      a == 0.0 && continue
      rb = row_block[i]
      rb == 0 && continue
      emit(rb, col_lambda, -a)
    end
    emit(row_ref, col_lambda, -alpha[slack_idx])
  end

  return nothing
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
function build_rectangular_jacobian_pq_pv_dense(Ybus, V::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int; dPinj_dVm::Vector{Float64} = zeros(Float64, length(V)), dQinj_dVm::Vector{Float64} = zeros(Float64, length(V)), vm_eps::Float64 = 1e-9, dslack::Union{Nothing,DistributedSlackState} = nothing)
  dslack === nothing || error("build_rectangular_jacobian_pq_pv_dense: distributed slack requires the sparse Jacobian path (power_flow.sparse = true).")
  # Dense variant is primarily for diagnostics/parity checks; sparse is default path.
  n = length(V)
  @assert length(bus_types) == n
  @assert length(Vset) == n
  @assert length(dPinj_dVm) == n
  @assert length(dQinj_dVm) == n

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

  pos_map = build_pos_map(non_slack, n)

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

      pos = pos_map[i]
      if pos != 0
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

  for i in non_slack
    rowP = 2 * pos_map[i] - 1
    pos = pos_map[i]
    vm = abs(V[i])
    vm_safe = vm > vm_eps ? vm : vm_eps
    dvm_dvr = real(V[i]) / vm_safe
    dvm_dvi = imag(V[i]) / vm_safe

    J[rowP, pos] -= dPinj_dVm[i] * dvm_dvr
    J[rowP, (n-1)+pos] -= dPinj_dVm[i] * dvm_dvi

    if bus_types[i] == :PQ
      rowQ = rowP + 1
      J[rowQ, pos] -= dQinj_dVm[i] * dvm_dvr
      J[rowQ, (n-1)+pos] -= dQinj_dVm[i] * dvm_dvi
    end
  end

  return J
end

"""
    build_rectangular_jacobian_pq_pv(
        Ybus::SparseMatrixCSC{ComplexF64},
        V,
        bus_types,
        Vset,
        slack_idx;
    )
"""
function build_rectangular_jacobian_pq_pv(
  Ybus::SparseMatrixCSC{ComplexF64},
  V::Vector{ComplexF64},
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  slack_idx::Int;
  dPinj_dVm::Vector{Float64} = zeros(Float64, length(V)),
  dQinj_dVm::Vector{Float64} = zeros(Float64, length(V)),
  vm_eps::Float64 = 1e-9,
  structural_pattern::Bool = false,
  assembly::Union{Nothing,RectangularJacobianAssembly} = nothing,
  dslack::Union{Nothing,DistributedSlackState} = nothing,
)
  # Default rectangular solver path uses sparse Jacobians for scale and parity.
  return build_rectangular_jacobian_pq_pv_sparse(Ybus, V, bus_types, Vset, slack_idx; dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm, vm_eps = vm_eps, structural_pattern = structural_pattern, assembly = assembly, dslack = dslack)
end

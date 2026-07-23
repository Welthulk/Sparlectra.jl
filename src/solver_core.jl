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

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 2024-06-01
# Purpose: core functions for power flow solver
# file: src/solver_core.jl

# ------------------------------------------------------------
# Helper: bus currents and complex power injections
# S = V .* conj(I), I = Y * V
# ------------------------------------------------------------

using LinearAlgebra
using SparseArrays

"""
    calc_currents(Y, V) -> I

Compute bus current injections I = Y * V.
"""
function calc_currents(Y::AbstractMatrix{ComplexF64}, V::AbstractVector{ComplexF64})
  return Y * V
end

"""
    calc_injections(Y, V) -> S

Compute complex bus power injections S = V .* conj(Y*V).
"""
function calc_injections(Y::AbstractMatrix{ComplexF64}, V::AbstractVector{ComplexF64})
  I = calc_currents(Y, V)
  return V .* conj.(I)
end

"""
    solve_sparse_system(A::SparseMatrixCSC, b; context=:powerflow)

Solve a sparse linear system without densifying `A`. The primary path uses the
standard sparse direct solve (UMFPACK for sparse floating-point matrices where
available), with sparse QR as a singular-step fallback. If both sparse paths
fail, the caller receives a clear error instead of an accidental dense SVD or
`Matrix(A)` fallback in the power-flow core.
"""
const DEFAULT_DENSE_FALLBACK_MAX_N = 64

function solve_sparse_system(A::SparseMatrixCSC, b; context::Symbol = :powerflow, diagnostics = nothing)
  backend = :sparse_backslash
  try
    x = A \ b
    _record_sparse_solve_backend!(diagnostics, backend, context, size(A))
    return x
  catch e
    _is_linear_singularity_error(e) || rethrow(e)
    try
      backend = :sparse_qr_fallback
      x = qr(A) \ b
      _record_sparse_solve_backend!(diagnostics, backend, context, size(A))
      return x
    catch qr_error
      _is_linear_fallback_error(qr_error) || rethrow(qr_error)
      _dense_svd_fallback_allowed(A) || throw(LinearAlgebra.SingularException(0))
      backend = :dense_svd_fallback_small_system
      x = _svd_pinv_solve(Matrix(A), b)
      _record_sparse_solve_backend!(diagnostics, backend, context, size(A))
      return x
    end
  end
end

"""
    solve_linear(A, b; allow_pinv=true)

Solve `A*x = b`. Sparse matrices are handled by `solve_sparse_system` and are
never converted to dense storage. Dense callers may still use the small-system
SVD fallback when explicitly allowed.
"""
function solve_linear(A, b; allow_pinv::Bool = true)
  if A isa SparseMatrixCSC
    return solve_sparse_system(A, b; context = :powerflow)
  end
  try
    return A \ b
  catch e
    if allow_pinv && _is_linear_singularity_error(e)
      try
        return qr(A) \ b
      catch qr_error
        if _is_linear_fallback_error(qr_error)
          _dense_svd_fallback_allowed(A) || rethrow(qr_error)
          return _svd_pinv_solve(A, b)
        end
        rethrow(qr_error)
      end
    end
    rethrow(e)
  end
end

_is_linear_singularity_error(e) = e isa LinearAlgebra.SingularException
_is_linear_fallback_error(e) = e isa LinearAlgebra.LAPACKException || e isa ArgumentError || e isa LinearAlgebra.SingularException

function _allfinite_matrix(A)::Bool
  if A isa SparseArrays.AbstractSparseMatrixCSC
    return all(isfinite, nonzeros(A))
  end
  return all(isfinite, A)
end

function _dense_svd_fallback_allowed(A)::Bool
  return max(size(A)...) <= DEFAULT_DENSE_FALLBACK_MAX_N
end

function _record_sparse_solve_backend!(diagnostics, backend::Symbol, context::Symbol, dims::Tuple{Int, Int})
  diagnostics isa AbstractDict || return nothing
  diagnostics[:backend] = backend
  diagnostics[:context] = context
  diagnostics[:dims] = dims
  return nothing
end

function _regularized_normal_solve(A::AbstractMatrix, b::AbstractVector)
  At = A'
  normal_matrix = At * A
  normal_rhs = At * b
  scale = _matrix_abs_scale(A)
  λ0 = max(sqrt(eps(Float64)) * max(scale * scale, 1.0), eps(Float64))
  for factor in (1.0, 1e2, 1e4, 1e6)
    λ = λ0 * factor
    x = (normal_matrix + λ * I) \ normal_rhs
    all(isfinite, x) && return x
  end
  error("regularized linear solve produced non-finite values")
end

function _matrix_abs_scale(A)::Float64
  if A isa SparseArrays.AbstractSparseMatrixCSC
    vals = nonzeros(A)
    isempty(vals) && return 0.0
    return Float64(maximum(abs, vals))
  end
  isempty(A) && return 0.0
  return Float64(maximum(abs, A))
end

function _svd_pinv_solve(A::AbstractMatrix, b)
  F = svd(A; full = false, alg = LinearAlgebra.QRIteration())
  isempty(F.S) && return zeros(eltype(b), size(A, 2))
  cutoff = max(size(A)...) * eps(real(eltype(F.S))) * maximum(F.S)
  Ub = F.U' * b
  Sinv_b = similar(Ub)
  @inbounds for i in eachindex(Ub)
    Sinv_b[i] = F.S[i] > cutoff ? Ub[i] / F.S[i] : zero(Ub[i])
  end
  return F.Vt' * Sinv_b
end
"""
    non_slack_indices(n, slack_idx) -> Vector{Int}
"""
non_slack_indices(n::Int, slack_idx::Int) = [i for i = 1:n if i != slack_idx]

"""
    build_pos_map(non_slack, n) -> pos::Vector{Int}
pos[i]=k, wenn i==non_slack[k], sonst 0.
"""
function build_pos_map(non_slack::Vector{Int}, n::Int)
  pos = zeros(Int, n)
  @inbounds for (k, b) in enumerate(non_slack)
    pos[b] = k
  end
  return pos
end

"""
    reduce_susceptance_to_nonslack(Bsrc, non_slack, pos, slack_idx, nred; value_of=imag) -> B_reduced

Eliminate the slack row/column from a full `n×n` susceptance-like source
matrix and return the `(n-1)×(n-1)` reduced system used by any DC-angle
linear solve (`_solve_dc_angle_system`). `value_of` extracts the real
susceptance entry from `Bsrc[i,j]`: `imag` for a complex Ybus (the
rectangular-solver DC-angle start-value heuristic, which reduces the full AC
admittance matrix), `identity` for an already-real B′ (a true DC power-flow
matrix, series reactance only). Sparse `Bsrc` in gives a sparse result;
dense in gives dense out. Shared between [`rundcpf!`](@ref)'s B′ reduction
and the rectangular solver's internal DC-angle start projection so the two
call sites never diverge into two slack-elimination implementations.
"""
function reduce_susceptance_to_nonslack(Bsrc, non_slack::Vector{Int}, pos::Vector{Int}, slack_idx::Int, nred::Int; value_of::Function = imag)
  if Bsrc isa SparseMatrixCSC
    I = Int[]
    J = Int[]
    V = Float64[]
    sizehint!(I, 2 * nnz(Bsrc))
    sizehint!(J, 2 * nnz(Bsrc))
    sizehint!(V, 2 * nnz(Bsrc))
    rv = rowvals(Bsrc)
    nz = nonzeros(Bsrc)
    @inbounds for j in axes(Bsrc, 2)
      for ptr in nzrange(Bsrc, j)
        i = rv[ptr]
        i == j && continue
        ri = pos[i]
        ri == 0 && continue
        bij = value_of(nz[ptr])
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
    return sparse(I, J, V, nred, nred)
  else
    B = zeros(Float64, nred, nred)
    @inbounds for i in non_slack
      ri = pos[i]
      for j in axes(Bsrc, 2)
        j == i && continue
        bij = value_of(Bsrc[i, j])
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
    return B
  end
end

"""
    slack_elimination_indices(n, slack_idx)

Return (non_slack, row_idx, col_idx) for [P;Q] and [Vr;Vi] layout of size 2n.
"""
function slack_elimination_indices(n::Int, slack_idx::Int)
  non_slack = non_slack_indices(n, slack_idx)
  idx = vcat(non_slack, n .+ non_slack)
  return non_slack, idx, idx
end

"""
    extract_bus_types_and_vset(net) -> (bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int)

bus_types uses :Slack/:PV/:PQ, Vset is vm setpoint (PV) or current vm default.
"""
function extract_bus_types_and_vset(net::Net)
  refreshBusTypesFromProsumers!(net)
  n         = length(net.nodeVec)
  bus_types = Vector{Symbol}(undef, n)
  Vset      = Vector{Float64}(undef, n)

  slack_idx = 0
  @inbounds for (k, node) in enumerate(net.nodeVec)
    nt = getNodeType(node)
    if nt == Slack
      bus_types[k] = :Slack
      slack_idx = k
    elseif nt == PV
      bus_types[k] = :PV
    elseif nt == PQ
      bus_types[k] = :PQ
    else
      error("extract_bus_types_and_vset: unsupported node type at bus $k")
    end
    Vset[k] = isnothing(node._vm_pu) ? 1.0 : node._vm_pu
  end
  slack_idx == 0 && error("extract_bus_types_and_vset: no slack bus found")
  return bus_types, Vset, slack_idx
end

"""
    build_qload_pu(net) -> Qload_pu::Vector{Float64}

Aggregated reactive load per bus in p.u. (sum of non-generator prosmps qVal / baseMVA).
"""
function build_qload_pu(net::Net)
  nb = length(net.nodeVec)
  Qload = zeros(Float64, nb)
  @inbounds for ps in net.prosumpsVec
    isGenerator(ps) && continue
    bus = getPosumerBusIndex(ps)
    (1 <= bus <= nb) || continue
    q = isnothing(ps.qVal) ? 0.0 : ps.qVal
    Qload[bus] += q / net.baseMVA
  end
  return Qload
end

"""
    build_voltage_vector(busVec) -> V::Vector{ComplexF64}
"""
build_voltage_vector(busVec::Vector{BusData}) = ComplexF64[bus.vm_pu * exp(im * bus.va_rad) for bus in busVec]

"""
    compute_sbus_and_totals(Y, V, baseMVA) -> (Sbus_pu, p_MW, q_MVar)
"""
function compute_sbus_and_totals(Y, V, baseMVA::Float64)
  Sbus_pu = calc_injections(Y, V)
  p = sum(real.(Sbus_pu)) * baseMVA
  q = sum(imag.(Sbus_pu)) * baseMVA
  return Sbus_pu, p, q
end

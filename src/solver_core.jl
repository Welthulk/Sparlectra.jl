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
    solve_linear(A, b; allow_pinv=true)

Solve `A*x = b`. If the direct solve reports a singular matrix and
`allow_pinv=true`, first try a rank-revealing QR solve. For smaller systems,
fall back to dense SVD if QR also fails. Large sparse systems deliberately do
not densify here; callers that can continue safely should catch the remaining
singular linear-step error and report non-convergence instead.
"""
function solve_linear(A, b; allow_pinv::Bool = true)
  try
    return A \ b
  catch e
    if allow_pinv && _is_linear_singularity_error(e)
      try
        # Prefer a rank-revealing QR fallback.  It avoids forming a dense SVD for
        # large sparse Jacobians and is more robust than `pinv` on singular steps.
        return qr(A) \ b
      catch qr_error
        if qr_error isa LinearAlgebra.LAPACKException || qr_error isa ArgumentError || qr_error isa LinearAlgebra.SingularException
          return _svd_pinv_solve(Matrix(A), b)
        end
        rethrow(qr_error)
      end
    end
    rethrow(e)
  end
end

function _svd_pinv_solve(A::AbstractMatrix, b)
  F = svd(A; full = false, alg = LinearAlgebra.QRIteration())
  isempty(F.S) && return zeros(eltype(b), size(A, 2))
  cutoff = max(size(A)...) * eps(real(eltype(F.S))) * maximum(F.S)
  Sinv_b = similar(F.U' * b)
  Ub = F.U' * b
  @inbounds for i in eachindex(Ub)
    Sinv_b[i] = F.S[i] > cutoff ? Ub[i] / F.S[i] : zero(Ub[i])
  end
  return F.Vt' * Sinv_b
end

_is_linear_singularity_error(e) = e isa LinearAlgebra.SingularException
_is_linear_fallback_error(e) = e isa LinearAlgebra.LAPACKException || e isa ArgumentError || e isa LinearAlgebra.SingularException

function _allfinite_matrix(A)::Bool
  if A isa SparseArrays.AbstractSparseMatrixCSC
    return all(isfinite, nonzeros(A))
  end
  return all(isfinite, A)
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
  return max(size(A)...) <= 2000
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

# Copyright 2023-2026 Udo Schmitz
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

# file: src/numerics/takahashi.jl
# purpose: shared Takahashi/Erisman-Tinney selected inverse on a UMFPACK
#          factorization: one backward pass over the filled factor pattern
#          instead of one triangular solve per column. Element-type generic
#          (ComplexF64 for the short-circuit sweep, Float64 for the SE
#          diagnostics). Entry points: takahashi_diag (the diagonal, exactly
#          the semantics the short-circuit module shipped as _takahashi_diag)
#          and takahashi_selected_inverse (every inv(A)[i, j] inside the
#          factor pattern, nothing outside).

"""
    TakahashiSelectedInverse{T}

Queryable selected inverse from [`takahashi_selected_inverse`](@ref):
`S[i, j]` returns `inv(A)[i, j]` when the (permuted) position lies inside
the filled factor pattern, and `nothing` outside it. Internally the entries
live on the pattern of `(L + U)^T` in permuted coordinates; the stored
permutation and row-scaling context map original indices back
(`inv(A)[i, j] = Z[pinv[i], pinv[j]] * Rs[j]`, since `L*U == (Rs .* A)[p, p]`
scales the ROWS of `A`).
"""
struct TakahashiSelectedInverse{T}
  Z::Dict{Tuple{Int,Int},T}
  pinv::Vector{Int}
  Rs::Vector{Float64}
end

@inline function Base.getindex(S::TakahashiSelectedInverse{T}, i::Int, j::Int)::Union{Nothing,T} where {T}
  z = get(S.Z, (S.pinv[i], S.pinv[j]), nothing)
  z === nothing && return nothing
  return z * S.Rs[j]
end

Base.size(S::TakahashiSelectedInverse) = (length(S.pinv), length(S.pinv))

## Shared backward pass: builds the Z entries on the filled pattern of
## (L + U)^T. Returns (Z, ok, info); every applicability guard of the
## original short-circuit implementation is kept verbatim so both entry
## points stay in sync:
## - p == q (symmetric pivot ordering; otherwise original-diagonal positions
##   leave the factor pattern),
## - nonzero pivots in U,
## - unit-lower L,
## - pattern-closure miss counter (defensive; a miss means the Erisman-
##   Tinney self-containment assumption was violated).
function _takahashi_pattern_pass(F)
  p = F.p
  q = F.q
  L = F.L
  T = eltype(L)
  p == q || return (Dict{Tuple{Int,Int},T}(), false, "unsymmetric pivot ordering (p != q)")
  U = F.U
  n = size(L, 1)
  dU = Vector{T}(undef, n)
  for j = 1:n
    dU[j] = U[j, j]
  end
  any(iszero, dU) && return (Dict{Tuple{Int,Int},T}(), false, "zero pivot in U")
  all(isone, diag(L)) || return (Dict{Tuple{Int,Int},T}(), false, "L is not unit-lower")

  # row access of U and L via their CSC transposes (column i = row i)
  Ut = SparseMatrixCSC(transpose(U))
  Lt = SparseMatrixCSC(transpose(L))
  Utc = Ut.colptr
  Utr = Ut.rowval
  Utv = Ut.nzval
  Lp = L.colptr
  Lr = L.rowval
  Lv = L.nzval
  Ltc = Lt.colptr
  Ltr = Lt.rowval

  # Z entries live on the pattern of (L + U)^T; per Erisman-Tinney the
  # recurrences restricted to that (filled) pattern are self-contained
  Z = Dict{Tuple{Int,Int},T}()
  sizehint!(Z, nnz(L) + nnz(U))
  misses = 0
  uppers = Int[]
  for j = n:-1:1
    # lower entries of column j (positions i > j with U[j, i] != 0):
    # Z[i, j] = -sum_{k > j} Z[i, k] * L[k, j]
    for t in Utc[j]:(Utc[j+1]-1)
      i = Utr[t]
      i > j || continue
      acc = zero(T)
      for s2 in Lp[j]:(Lp[j+1]-1)
        k = Lr[s2]
        k > j || continue
        z = get(Z, (i, k), nothing)
        z === nothing ? (misses += 1) : (acc += z * Lv[s2])
      end
      Z[(i, j)] = -acc
    end
    # diagonal: Z[j, j] = 1/d_j - sum_{k > j} (U[j, k]/d_j) * Z[k, j]
    accd = zero(T)
    for t in Utc[j]:(Utc[j+1]-1)
      k = Utr[t]
      k > j || continue
      z = get(Z, (k, j), nothing)
      z === nothing ? (misses += 1) : (accd += (Utv[t] / dU[j]) * z)
    end
    Z[(j, j)] = inv(dU[j]) - accd
    # upper entries of column j (positions i < j with L[j, i] != 0),
    # descending i: Z[i, j] = -sum_{k > i} (U[i, k]/d_i) * Z[k, j]
    empty!(uppers)
    for t in Ltc[j]:(Ltc[j+1]-1)
      Ltr[t] < j && push!(uppers, Ltr[t])
    end
    sort!(uppers; rev = true)
    for i in uppers
      acc = zero(T)
      for t in Utc[i]:(Utc[i+1]-1)
        k = Utr[t]
        k > i || continue
        z = get(Z, (k, j), nothing)
        z === nothing ? (misses += 1) : (acc += (Utv[t] / dU[i]) * z)
      end
      Z[(i, j)] = -acc
    end
  end
  misses == 0 || return (Dict{Tuple{Int,Int},T}(), false, "pattern closure violated at $(misses) lookups")
  return (Z, true, "ok")
end

"""
    takahashi_diag(F) -> (diagZ::Vector, ok::Bool, info::String)

Compute the full diagonal of `inv(A)` from the UMFPACK factorization `F`
of `A` via the Takahashi/Erisman-Tinney selected-inverse recurrences on
the filled factor pattern: one backward pass over `nnz(L) + nnz(U)`
entries instead of one triangular solve per column (measured 34x to 264x
over the serial short-circuit sweep between n = 2000 and n = 16000,
agreeing with direct solves to about 1e-15 relative). The element type
follows the factorization (ComplexF64 for the short-circuit Ybus, Float64
for the SE gain matrix).

Applicability guard: requires the symmetric pivot ordering (`F.p == F.q`)
that UMFPACK's symmetric strategy picks on structurally symmetric
matrices; with an unsymmetric ordering the original-diagonal positions
leave the factor pattern. Returns `ok = false` with the reason in `info`
in that case (or on a zero pivot, a non-unit-lower L, or a pattern-closure
violation, counted defensively); callers must fall back to direct solves.
"""
function takahashi_diag(F)
  T = eltype(F.L)
  Z, ok, info = _takahashi_pattern_pass(F)
  ok || return (T[], false, info)
  # map back: with L*U == (Rs .* A)[p, q] and p == q,
  # inv(A)[i, i] = Rs[i] * Z[r, r] at r = pinv[i]
  pinv_ = invperm(F.p)
  n = length(pinv_)
  diagZ = Vector{T}(undef, n)
  for i = 1:n
    z = get(Z, (pinv_[i], pinv_[i]), nothing)
    z === nothing && return (T[], false, "diagonal position missing from pattern")
    diagZ[i] = F.Rs[i] * z
  end
  return (diagZ, true, "ok")
end

"""
    takahashi_selected_inverse(F) -> (S::Union{Nothing,TakahashiSelectedInverse}, ok::Bool, info::String)

Same backward pass as [`takahashi_diag`](@ref) (identical guards and code
path), but returns the FULL selected inverse on the factor pattern as a
queryable [`TakahashiSelectedInverse`](@ref): `S[i, j]` yields
`inv(A)[i, j]` for any position inside the pattern and `nothing` outside.
This is what the SE diagnostics need for `Omega_ii`: every `G^-1` entry
they look up sits at a structural nonzero of `G` itself, which is contained
in the filled pattern of its factors, so no in-pattern lookup can miss.
"""
function takahashi_selected_inverse(F)
  Z, ok, info = _takahashi_pattern_pass(F)
  ok || return (nothing, false, info)
  return (TakahashiSelectedInverse(Z, invperm(F.p), Vector{Float64}(F.Rs)), true, "ok")
end

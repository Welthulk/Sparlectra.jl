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
# Purpose: factorization-reuse linear-solver backends (KLU, UMFPACK lu!) for
#          the rectangular Newton step
# file: src/powerflow_rectangular/rectangular_klu_solver.jl

using KLU: KLU, klu, klu!

"""
    AbstractNewtonSolverContext

Common supertype of the factorization-reuse contexts for the rectangular
Newton step (`power_flow.linear_solver`). All concrete contexts hold the
reusable factorization, a fingerprint of the analyzed sparsity pattern
(`colptr`/`rowval` copies), and diagnostic counters. Lifetime is one
`runpf_rectangular!` invocation (= one island solve on the island path); a
context must never be shared across islands or threads.
"""
abstract type AbstractNewtonSolverContext end

"""
    KLUNewtonContext

Reuse context for the SuiteSparse KLU backend (`linear_solver = :klu`):
symbolic analysis once per sparsity pattern, `klu!` numeric refactorization
per iteration.
"""
mutable struct KLUNewtonContext <: AbstractNewtonSolverContext
  fact::Union{Nothing,KLU.KLUFactorization{Float64,Int64}}
  nvar::Int
  colptr::Vector{Int64}
  rowval::Vector{Int64}
  analyze_count::Int
  refactor_count::Int
  fallback_count::Int
end

KLUNewtonContext() = KLUNewtonContext(nothing, 0, Int64[], Int64[], 0, 0, 0)

"""
    UmfpackReuseNewtonContext

Reuse context for the UMFPACK backend (`linear_solver = :umfpack_reuse`):
the symbolic analysis of the first `lu` factorization is reused via
`lu!(F, J)` in subsequent iterations. Same multifrontal factorization as the
default `umfpack` path, minus the repeated symbolic-analysis cost.
"""
mutable struct UmfpackReuseNewtonContext <: AbstractNewtonSolverContext
  fact::Union{Nothing,SparseArrays.UMFPACK.UmfpackLU{Float64,Int64}}
  nvar::Int
  colptr::Vector{Int64}
  rowval::Vector{Int64}
  analyze_count::Int
  refactor_count::Int
  fallback_count::Int
end

UmfpackReuseNewtonContext() = UmfpackReuseNewtonContext(nothing, 0, Int64[], Int64[], 0, 0, 0)

_newton_full_factorization(::KLUNewtonContext, J::SparseMatrixCSC{Float64,Int64}) = klu(J)
_newton_full_factorization(::UmfpackReuseNewtonContext, J::SparseMatrixCSC{Float64,Int64}) = lu(J)
_newton_refactorize!(ctx::KLUNewtonContext, J::SparseMatrixCSC{Float64,Int64}) = klu!(ctx.fact, J)
_newton_refactorize!(ctx::UmfpackReuseNewtonContext, J::SparseMatrixCSC{Float64,Int64}) = lu!(ctx.fact, J)

"""
    solve_newton_factorized!(ctx, J, rhs; pattern_changed=false) -> x

Solve `J * x = rhs` through the backend of `ctx`, reusing the symbolic
analysis stored in the context across Newton iterations. A full analyze +
factor runs when `pattern_changed` is set (PV↔PQ active-set switch), when no
factorization exists yet, or when the structural fingerprint of `J` differs
from the analyzed pattern; otherwise only the numeric refactorization runs.
The structural guard is mandatory: silent pattern drift fed into a
refactorization could produce wrong numbers instead of an exception.

Any backend failure (structure mismatch, bad pivots, singularity) resets the
context and falls back to `solve_sparse_system(J, rhs; context = :powerflow)`
so the existing UMFPACK/QR/SVD chain and the outer loop's
`:singular_newton_step` handling keep working unchanged.
"""
function solve_newton_factorized!(ctx::AbstractNewtonSolverContext, J::SparseMatrixCSC{Float64,Int64}, rhs::AbstractVector{Float64}; pattern_changed::Bool = false)
  need_analyze =
    pattern_changed ||
    ctx.fact === nothing ||
    size(J, 1) != ctx.nvar ||
    ctx.colptr != J.colptr ||
    ctx.rowval != J.rowval
  try
    if need_analyze
      ctx.fact = _newton_full_factorization(ctx, J)
      ctx.nvar = size(J, 1)
      ctx.colptr = copy(J.colptr)
      ctx.rowval = copy(J.rowval)
      ctx.analyze_count += 1
    else
      _newton_refactorize!(ctx, J)
      ctx.refactor_count += 1
    end
    return ctx.fact \ rhs
  catch e
    e isa InterruptException && rethrow(e)
    # Bad pivots, singularity, or a structure mismatch that slipped past the
    # guard: drop the factorization and delegate to the shared fallback chain.
    ctx.fact = nothing
    ctx.fallback_count += 1
    return solve_sparse_system(J, rhs; context = :powerflow)
  end
end

"""
    solve_newton_klu!(ctx, J, rhs; pattern_changed=false) -> x

KLU-specific name kept for readability at KLU call sites; identical to
[`solve_newton_factorized!`](@ref).
"""
solve_newton_klu!(ctx::KLUNewtonContext, J::SparseMatrixCSC{Float64,Int64}, rhs::AbstractVector{Float64}; pattern_changed::Bool = false) =
  solve_newton_factorized!(ctx, J, rhs; pattern_changed = pattern_changed)

"""
    newton_linear_solver_context(linear_solver::Symbol)

Map a validated `power_flow.linear_solver` value to a fresh per-solve
context: `:klu` → [`KLUNewtonContext`](@ref), `:umfpack_reuse` →
[`UmfpackReuseNewtonContext`](@ref), `:umfpack` → `nothing` (default direct
`solve_linear` path).
"""
function newton_linear_solver_context(linear_solver::Symbol)
  linear_solver === :klu && return KLUNewtonContext()
  linear_solver === :umfpack_reuse && return UmfpackReuseNewtonContext()
  return nothing
end

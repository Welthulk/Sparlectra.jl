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

# Author: Udo Schmitz (https://github.com/Welthulk)
# file: src/powerflow_rectangular/rectangular_factorized_solver.jl
# purpose: factorization-reuse linear-solver backend (UMFPACK lu!) for the
#          rectangular Newton step. A KLU backend lived here until 0.9.10;
#          it was removed after measurements showed no benefit on power-flow
#          Jacobians (slower factorization on real cases, unsafe under
#          shared-factor threading).

"""
    AbstractNewtonSolverContext

Common supertype of the factorization-reuse contexts for the rectangular
Newton step (`power_flow.linear_solver`). A concrete context holds the
reusable factorization, a fingerprint of the analyzed sparsity pattern
(`colptr`/`rowval` copies), and diagnostic counters. Lifetime is one
`runpf_rectangular!` invocation (= one island solve on the island path); a
context must never be shared across islands or threads.
"""
abstract type AbstractNewtonSolverContext end

"""
    RectangularJacobianAssembly

Reusable buffers for the rectangular Jacobian build (issue #292): the COO
triplet vectors, the retained `SparseMatrixCSC`, and a position map from
emit order to `nzval` index. Under the structural (value-independent)
sparsity pattern of the reuse backend, the emit order is identical in every
iteration with an unchanged active set, so subsequent builds can overwrite
`nzval` in place instead of allocating fresh triplets plus a new CSC.
`valid` is cleared on PV↔PQ active-set changes (pattern changes) and by the
runtime guard when a replay does not reproduce the recorded emit count.

`fixed_pattern` (`power_flow.jacobian_fixed_pattern`, 1.2.3) selects the
fixed Jacobian pattern: PV and PQ rows carry the same entries, so `valid` is
not cleared on active-set changes and the linear context keeps its analysis.
`factor_bus_types` holds, in fixed mode only, the bus types of the Jacobian
of the last full factorization; a backend whose refactorization keeps the
old pivots (KLU) factorizes afresh when the types differ from them.
`dslack_buses` holds the non-slack buses with a nonzero participation
factor of the recorded build (the rows of the lambda column); a replay with
another participant set rebuilds, because the same emit count would write
-alpha into the old slots.
`layout_n`, `layout_slack`, `non_slack`, `pos_non_slack` and `row_block`
hold the state layout of the bus count and slack index last seen
(`_assembly_layout!`), so a Newton iteration does not build them again.
"""
mutable struct RectangularJacobianAssembly
  Iidx::Vector{Int}
  Jidx::Vector{Int}
  Vals::Vector{Float64}
  J::Union{Nothing,SparseMatrixCSC{Float64,Int64}}
  pos_map::Vector{Int}
  Ibuf::Vector{ComplexF64}
  valid::Bool
  fixed_pattern::Bool
  factor_bus_types::Vector{Symbol}
  dslack_buses::Vector{Int}
  layout_n::Int
  layout_slack::Int
  non_slack::Vector{Int}
  pos_non_slack::Vector{Int}
  row_block::Vector{Int}
end

RectangularJacobianAssembly() = RectangularJacobianAssembly(Int[], Int[], Float64[], nothing, Int[], ComplexF64[], false, false, Symbol[], Int[], -1, -1, Int[], Int[], Int[])

# The non-slack bus list, its position map and the row block of each bus
# (the index of its P row) for n buses and this slack index, the same
# vectors non_slack_indices, build_pos_map and the row loop of the Jacobian
# builder give. Built again only when n or the slack index differ from the
# last call. The vectors are read only by their users.
function _assembly_layout!(asm::RectangularJacobianAssembly, n::Int, slack_idx::Int)
  if asm.layout_n != n || asm.layout_slack != slack_idx
    asm.non_slack = non_slack_indices(n, slack_idx)
    asm.pos_non_slack = build_pos_map(asm.non_slack, n)
    row_block = zeros(Int, n)
    row = 0
    for i = 1:n
      i == slack_idx && continue
      row += 2
      row_block[i] = row - 1
    end
    asm.row_block = row_block
    asm.layout_n = n
    asm.layout_slack = slack_idx
  end
  return asm.non_slack, asm.pos_non_slack, asm.row_block
end

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
  rhs::Vector{Float64}
  sol::Vector{Float64}
  assembly::RectangularJacobianAssembly
  # why this CPU context serves a solve for which power_flow.compute_backend
  # asked for the GPU (empty when the CPU was the configured backend)
  backend_note::String
end

UmfpackReuseNewtonContext() = UmfpackReuseNewtonContext(nothing, 0, Int64[], Int64[], 0, 0, 0, Float64[], Float64[], RectangularJacobianAssembly(), "")

_newton_full_factorization(::UmfpackReuseNewtonContext, J::SparseMatrixCSC{Float64,Int64}) = lu(J)
_newton_refactorize!(ctx::UmfpackReuseNewtonContext, J::SparseMatrixCSC{Float64,Int64}) = lu!(ctx.fact, J)
# the backend name a context reports (the KLU extension adds :klu)
_newton_context_backend(::UmfpackReuseNewtonContext) = :umfpack_reuse
# Whether the numeric refactorization of a context chooses its pivots anew.
# UMFPACK's lu! runs threshold partial pivoting on the kept symbolic
# analysis, so a PV/PQ switch under the fixed Jacobian pattern needs no new
# analysis. KLU's klu! (klu_l_refactor) reuses the old pivot sequence, which
# may be wrong for switched rows; it and any context not known to pivot get
# false and factorize afresh after a switch.
_newton_refactor_pivots(::AbstractNewtonSolverContext) = false
_newton_refactor_pivots(::UmfpackReuseNewtonContext) = true

"""
    solve_newton_factorized!(ctx, J, rhs; pattern_changed=false) -> x

Solve `J * x = rhs` through the backend of `ctx`, reusing the symbolic
analysis stored in the context across Newton iterations. A full analyze +
factor runs when `pattern_changed` is set (PV↔PQ active-set switch; under the
fixed Jacobian pattern only for a backend whose refactorization keeps the old
pivots, see `_newton_refactor_pivots`), when no
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
    # The context-owned solution buffer is overwritten on the next call; the
    # Newton step consumes δx within one iteration and never retains it.
    resize!(ctx.sol, length(rhs))
    ldiv!(ctx.sol, ctx.fact, rhs)
    return ctx.sol
  catch e
    e isa InterruptException && rethrow(e)
    # Bad pivots, singularity, or a structure mismatch that slipped past the
    # guard: drop the factorization and delegate to the shared fallback chain.
    ctx.fact = nothing
    ctx.fallback_count += 1
    return solve_sparse_system(J, rhs; context = :powerflow)
  end
end

# the CPU context of a validated power_flow.linear_solver value:
# :umfpack_reuse gives an UmfpackReuseNewtonContext, :umfpack nothing (the
# direct solve_linear path); newton_linear_solver_context adds the GPU
_cpu_newton_linear_solver_context(linear_solver::Symbol) = linear_solver === :umfpack_reuse ? UmfpackReuseNewtonContext() : nothing

# Decides and sets the fixed Jacobian pattern (power_flow.jacobian_fixed_pattern)
# on the context of one solve of `nb` buses (one island on the island path)
# and returns `(used, reason)`. `mode` is :auto, :on or :off; :auto uses the
# pattern from `min_buses` buses on. Without a context (linear solver umfpack)
# there is no analysis to keep: the pattern is not used, and only a
# configured :on says so in the console. A change of the decision
# invalidates the recorded assembly: its pos_map belongs to the other emit
# sequence (a power-mode net solved first with and then without it).
# Reasons: :off, :no_context, :on, :auto_threshold, :auto_below_threshold.
function _set_jacobian_fixed_pattern!(ctx, mode::Symbol, min_buses::Int, nb::Int, linear_solver::Symbol, verbose::Int)
  mode in JACOBIAN_FIXED_PATTERN_VALUES || throw(ArgumentError("jacobian_fixed_pattern must be one of $(JACOBIAN_FIXED_PATTERN_VALUES), got :$(mode)"))
  used, reason = if mode === :off
    false, :off
  elseif ctx === nothing
    if mode === :on && verbose > 0
      println("Fixed Jacobian pattern ignored: linear_solver $(linear_solver) keeps no factorisation; use umfpack_reuse.")
    end
    false, :no_context
  elseif mode === :on
    true, :on
  else
    nb >= min_buses ? (true, :auto_threshold) : (false, :auto_below_threshold)
  end
  ctx === nothing && return used, reason
  asm = ctx.assembly
  if asm.fixed_pattern != used
    asm.fixed_pattern = used
    asm.valid = false
    empty!(asm.factor_bus_types)
  end
  return used, reason
end

# The status and run-log text of the fixed-pattern decision: whether the
# pattern was used, the configured value as written in the configuration
# (auto, true, false) and the reason.
function _jacobian_fixed_pattern_line(mode::Symbol, reason::Symbol, used::Bool, nb::Int, min_buses::Int, linear_solver::Symbol)::String
  configured = mode === :auto ? "auto" : mode === :on ? "true" : "false"
  why = if reason === :auto_threshold
    ", $(nb) buses, threshold $(min_buses)"
  elseif reason === :auto_below_threshold
    ", $(nb) buses below threshold $(min_buses)"
  elseif reason === :no_context
    ", linear_solver $(linear_solver) keeps no factorisation"
  else
    ""
  end
  return string(used ? "on" : "off", " (", configured, why, ")")
end


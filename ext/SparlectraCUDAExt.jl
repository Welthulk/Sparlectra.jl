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
#
# file: ext/SparlectraCUDAExt.jl
# purpose: the Newton linear solve on an NVIDIA GPU with cuDSS
#          (power_flow.compute_backend). Loaded by `using CUDA, CUDSS` next
#          to Sparlectra; the base package does not depend on CUDA. The
#          Jacobian is assembled on the CPU, its values go to the GPU in
#          every step; the symbolic analysis runs once per sparsity pattern,
#          every further step only refactorizes. Precision float64, or
#          float32_mixed: factorization and solve in Float32, mismatch and
#          voltage update stay Float64 on the CPU side. cuDSS runs in its
#          deterministic mode unless gpu_fast_mode is on. Inside
#          Sparlectra._with_device_newton (GPU study) the whole Newton loop
#          runs on the GPU instead (SparlectraCUDAExt/newton_loop.jl).
#          runInjectionBatch with backend :gpu runs N load scenarios of one
#          network side by side on the device (SparlectraCUDAExt/batch_loop.jl),
#          or with gpu_two_phase in a Float32 phase and a Float64 phase
#          (SparlectraCUDAExt/batch_two_phase.jl).

module SparlectraCUDAExt

using Sparlectra
using CUDA
using CUDA.CUSPARSE: CuSparseMatrixCSR
using CUDSS
using SparseArrays
using LinearAlgebra

# The cuDSS factorization of one Jacobian pattern. cuDSS takes CSR, the
# Jacobian is CSC: `perm[k]` is the CSC value index that lands at CSR value
# position k (computed once per pattern), so a refactorization only gathers
# and uploads the values. `owner` is the context whose timers it feeds.
mutable struct CudssFactorization{T}
  solver::CudssSolver{T}
  A::CuSparseMatrixCSR{T,Int32}
  perm::Vector{Int}
  nzbuf::Vector{T}
  x::CuVector{T}
  b::CuVector{T}
  hostbuf::Vector{T}
  owner::Any
end

mutable struct CudssNewtonContext{T} <: Sparlectra.AbstractGpuNewtonContext
  fact::Union{Nothing,CudssFactorization{T}}
  nvar::Int
  colptr::Vector{Int64}
  rowval::Vector{Int64}
  analyze_count::Int
  refactor_count::Int
  fallback_count::Int
  rhs::Vector{Float64}
  sol::Vector{Float64}
  assembly::Sparlectra.RectangularJacobianAssembly
  precision::Symbol
  fallback_to_cpu::Bool
  fast_mode::Bool
  cpu::Union{Nothing,Sparlectra.UmfpackReuseNewtonContext}
  note::String
  fallback_reason::String
  # host <-> device copies and GPU work (analysis, factorization, solve),
  # each measured to completion
  transfer_s::Float64
  gpu_s::Float64
  # the device Newton loop ran for this solve, and its timings
  device_newton::Bool
  device_timings::Union{Nothing,Sparlectra.DeviceNewtonTimings}
end

function CudssNewtonContext{T}(precision::Symbol, fallback_to_cpu::Bool, fast_mode::Bool, note::String) where {T}
  return CudssNewtonContext{T}(nothing, 0, Int64[], Int64[], 0, 0, 0, Float64[], Float64[], Sparlectra.RectangularJacobianAssembly(), precision, fallback_to_cpu, fast_mode, nothing, note, "", 0.0, 0.0, false, nothing)
end

# seconds of `f()` with the device work finished
function _timed(f)
  t0 = time_ns()
  result = f()
  CUDA.synchronize()
  return result, (time_ns() - t0) / 1.0e9
end

# CSR structure of the CSC matrix J and the value permutation CSC -> CSR
function _csr_structure(J::SparseMatrixCSC{Float64,Int64})
  n = size(J, 1)
  idx = SparseMatrixCSC(n, size(J, 2), J.colptr, J.rowval, collect(1:nnz(J)))
  Jt = copy(transpose(idx))
  return Int32.(Jt.colptr), Int32.(Jt.rowval), Jt.nzval
end

# cuDSS runs non-deterministic by default: two factorizations and solves of
# one matrix differ in the last bits (benchmark/gpu/check_device_newton_gpu.jl
# part 1), and on case_ACTIVSg70k that decides which buses end at their Q
# limit. Every solver of the extension is created with deterministic_mode = 1
# unless power_flow.compute_backend.gpu_fast_mode is on.
function _cudss_configure!(solver, fast_mode::Bool)
  fast_mode || cudss_set(solver, "deterministic_mode", 1)
  return solver
end

# cuDSS replaces a zero pivot by pivot_epsilon and the solve returns a
# finite but wrong step; "npivots" counts the replaced pivots, "info" stays
# 0 (CUDSS.jl 1.0.0, benchmark/gpu/check_device_newton_gpu.jl part 1; a row
# of stored zeros gives a non-finite solution instead, caught by the finite
# check). A perturbed pivot is a failed GPU solve: the error sends the step
# to UMFPACK, which throws for a singular Jacobian as in the CPU path. The
# device Newton loop sets pivot_epsilon to 0 instead (newton_loop.jl).
function _check_pivots(solver, precision::Symbol)
  np = cudss_get(solver, "npivots")
  np > 0 && error("cuDSS perturbed $(np) pivot(s) (singular or ill-conditioned Jacobian in $(precision))")
  return nothing
end

function _cudss_factorize(ctx::CudssNewtonContext{T}, J::SparseMatrixCSC{Float64,Int64}) where {T}
  n = size(J, 1)
  rowptr, colval, perm = _csr_structure(J)
  nzbuf = T.(J.nzval[perm])
  (A, x, b), t_up = _timed() do
    (CuSparseMatrixCSR{T,Int32}(CuVector{Int32}(rowptr), CuVector{Int32}(colval), CuVector{T}(nzbuf), (n, n)), CUDA.zeros(T, n), CUDA.zeros(T, n))
  end
  ctx.transfer_s += t_up
  solver, t_fact = _timed() do
    s = _cudss_configure!(CudssSolver(A, "G", 'F'), ctx.fast_mode)
    cudss("analysis", s, x, b)
    cudss("factorization", s, x, b)
    s
  end
  ctx.gpu_s += t_fact
  _check_pivots(solver, ctx.precision)
  return CudssFactorization{T}(solver, A, perm, nzbuf, x, b, Vector{T}(undef, n), ctx)
end

function _cudss_refactorize!(ctx::CudssNewtonContext{T}, J::SparseMatrixCSC{Float64,Int64}) where {T}
  F = ctx.fact
  @inbounds for k in eachindex(F.perm)
    F.nzbuf[k] = T(J.nzval[F.perm[k]])
  end
  _, t_up = _timed(() -> copyto!(F.A.nzVal, F.nzbuf))
  ctx.transfer_s += t_up
  _, t_fact = _timed() do
    cudss_update(F.solver, F.A)
    cudss("refactorization", F.solver, F.x, F.b)
  end
  ctx.gpu_s += t_fact
  _check_pivots(F.solver, ctx.precision)
  return F
end

Sparlectra._newton_full_factorization(ctx::CudssNewtonContext, J::SparseMatrixCSC{Float64,Int64}) = _cudss_factorize(ctx, J)
Sparlectra._newton_refactorize!(ctx::CudssNewtonContext, J::SparseMatrixCSC{Float64,Int64}) = _cudss_refactorize!(ctx, J)

# sol = J \ rhs on the GPU; rhs and sol stay Float64 on the host
function LinearAlgebra.ldiv!(sol::AbstractVector{Float64}, F::CudssFactorization{T}, rhs::AbstractVector{Float64}) where {T}
  ctx = F.owner
  F.hostbuf .= T.(rhs)
  _, t_in = _timed(() -> copyto!(F.b, F.hostbuf))
  _, t_solve = _timed(() -> cudss("solve", F.solver, F.x, F.b))
  _, t_out = _timed(() -> copyto!(F.hostbuf, F.x))
  sol .= F.hostbuf
  ctx.transfer_s += t_in + t_out
  ctx.gpu_s += t_solve
  return sol
end

function _gpu_status()
  CUDA.functional() || return (status = :not_functional, reason = "CUDA.functional() is false (no usable driver or device; see CUDA.versioninfo())", device = "")
  return (status = :available, reason = "", device = CUDA.name(CUDA.device()))
end

function _gpu_context(precision::Symbol, fallback_to_cpu::Bool, fast_mode::Bool, note::String)
  precision === :float32_mixed && return CudssNewtonContext{Float32}(precision, fallback_to_cpu, fast_mode, note)
  return CudssNewtonContext{Float64}(precision, fallback_to_cpu, fast_mode, note)
end

include("SparlectraCUDAExt/device_arith.jl")
include("SparlectraCUDAExt/newton_map.jl")
include("SparlectraCUDAExt/newton_kernels.jl")
include("SparlectraCUDAExt/newton_loop.jl")
include("SparlectraCUDAExt/batch_loop.jl")
include("SparlectraCUDAExt/batch_two_phase.jl")

# the status of a solve that ran the device Newton loop: the stage 1 extras
# plus the device timings and the host copies of the loop
function Sparlectra._newton_context_extras(ctx::CudssNewtonContext)
  base = invoke(Sparlectra._newton_context_extras, Tuple{Sparlectra.AbstractGpuNewtonContext}, ctx)
  tm = ctx.device_timings
  (ctx.device_newton && tm !== nothing) || return base
  return merge(base, (
    gpu_newton = true,
    gpu_newton_build_s = tm.build_s,
    gpu_newton_mismatch_s = tm.mismatch_s,
    gpu_newton_refactor_s = tm.refactor_s,
    gpu_newton_solve_s = tm.solve_s,
    gpu_newton_update_s = tm.update_s,
    gpu_newton_analysis_s = tm.analysis_s,
    gpu_newton_factor_first_s = tm.factor_first_s,
    gpu_newton_transfer_s = tm.upload_s + tm.d2h_s + tm.h2d_s,
    gpu_newton_host_qlimit_s = tm.host_qlimit_s,
    gpu_newton_core_s = tm.core_s,
    gpu_newton_first_call = tm.first_call,
    gpu_newton_d2h_bytes = tm.d2h_bytes,
    gpu_newton_d2h_count = tm.d2h_count,
    gpu_newton_h2d_bytes = tm.h2d_bytes,
    gpu_newton_h2d_count = tm.h2d_count,
  ))
end

function __init__()
  Sparlectra._GPU_BACKEND[] = (status = _gpu_status, context = _gpu_context, newton = _gpu_newton, batch = _gpu_injection_batch)
  return nothing
end

end

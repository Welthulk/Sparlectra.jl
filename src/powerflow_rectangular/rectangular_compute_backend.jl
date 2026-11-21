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

# file: src/powerflow_rectangular/rectangular_compute_backend.jl
# purpose: power_flow.compute_backend: the GPU status, the GPU advisor, and
#          the choice of the Newton linear-solver context (UMFPACK on the CPU
#          or the cuDSS context of ext/SparlectraCUDAExt.jl). The core never
#          loads CUDA; without the extension every GPU request runs on the
#          CPU with the reason recorded.

# The compute backend of the solve in progress. _runpf_config_once! sets it
# from power_flow.compute_backend for the duration of one run; the island
# tasks of that run inherit it. Outside a configured run (a direct
# runpf_rectangular! call) it is the CPU default.
const _COMPUTE_BACKEND = Base.ScopedValues.ScopedValue(ComputeBackendConfig())

# run `f` with `cb` as the compute backend of every Newton solve inside it
_with_compute_backend(@nospecialize(f), cb::ComputeBackendConfig) = Base.ScopedValues.with(f, _COMPUTE_BACKEND => cb)

# The CUDA extension registers itself here in its __init__ with a
# NamedTuple of three functions:
#   status()                            -> (status, reason, device)
#   context(precision, fallback, fast_mode, note)
#                                       -> an AbstractGpuNewtonContext
#   newton(input::DeviceNewtonInput)    -> the device Newton loop
#                                          (rectangular_device_newton.jl)
const _GPU_BACKEND = Ref{Any}(nothing)

"""
    gpu_status() -> NamedTuple

Whether Newton steps can run on an NVIDIA GPU, as
`(status, reason, device)`: `status` is `:unavailable` (the CUDA extension is
not loaded: load `CUDA` and `CUDSS` next to Sparlectra), `:not_functional`
(the extension is loaded but `CUDA.functional()` is false; `reason` says
why) or `:available` (`device` names the GPU). Never throws.
"""
function gpu_status()::NamedTuple{(:status, :reason, :device),Tuple{Symbol,String,String}}
  backend = _GPU_BACKEND[]
  backend === nothing && return (status = :unavailable, reason = "CUDA extension not loaded (load CUDA and CUDSS next to Sparlectra)", device = "")
  return backend.status()
end

"""
    GpuAdvisorRow

One measured range of the GPU advisor: for networks of `min_buses` to
`max_buses` buses, `batch` scenarios and `precision`, the GPU took `factor`
times less time than the CPU in steady state (transfer included). `run`
names the study run (date, commit, machine) the factor comes from.
"""
struct GpuAdvisorRow
  precision::Symbol
  batch::Int
  min_buses::Int
  max_buses::Int
  factor::Float64
  run::String
end

# Measured ranges, filled from the GPU study; empty
# until the study has run, so the advisor reports :no_measurement.
const GPU_ADVISOR_ROWS = GpuAdvisorRow[]

# a GPU factor at or above this counts as faster (the study's criterion)
const GPU_ADVISOR_MIN_FACTOR = 1.2

"""
    gpu_advice(nbus; batch = 1, precision = :float64) -> NamedTuple

The advisor's verdict for a network of `nbus` buses solved in batches of
`batch` scenarios: `(decision, factor, run)` with `decision` one of
`:gpu_faster` (measured factor at least 1.2), `:gpu_slower` or
`:no_measurement` (no study row covers this size and batch). `factor` is
the measured CPU time over GPU time (`NaN` without a measurement) and `run`
the study run it comes from. Based only on measurements, never on the
hardware datasheet.
"""
function gpu_advice(nbus::Integer; batch::Integer = 1, precision::Symbol = :float64)::NamedTuple{(:decision, :factor, :run),Tuple{Symbol,Float64,String}}
  for row in GPU_ADVISOR_ROWS
    (row.precision === precision && row.batch == batch && row.min_buses <= nbus <= row.max_buses) || continue
    return (decision = row.factor >= GPU_ADVISOR_MIN_FACTOR ? :gpu_faster : :gpu_slower, factor = row.factor, run = row.run)
  end
  return (decision = :no_measurement, factor = NaN, run = "")
end

"""
    AbstractGpuNewtonContext <: AbstractNewtonSolverContext

A Newton linear-solver context that factorizes on a GPU (implemented by the
CUDA extension). Besides the fields every context has (`fact`, `nvar`,
`colptr`, `rowval`, the counters, `rhs`, `sol`, `assembly`) it carries
`precision`, `fallback_to_cpu`, `fast_mode` (cuDSS without its
deterministic mode), `cpu` (the UMFPACK context it switched to after a GPU
error, `nothing` before), `note` (why the GPU was chosen),
`fallback_reason` and the timings `transfer_s` and `gpu_s`.
"""
abstract type AbstractGpuNewtonContext <: AbstractNewtonSolverContext end

# The GPU solve with the CPU as the way out: after a GPU error (an
# exception of CUDA or cuDSS, or a non-finite solution) the context
# switches to UMFPACK for the rest of this solve when fallback is allowed,
# and records why; otherwise the error goes up.
function solve_newton_factorized!(ctx::AbstractGpuNewtonContext, J::SparseMatrixCSC{Float64,Int64}, rhs::AbstractVector{Float64}; pattern_changed::Bool = false)
  ctx.cpu === nothing || return solve_newton_factorized!(ctx.cpu, J, rhs; pattern_changed = pattern_changed)
  need_analyze = pattern_changed || ctx.fact === nothing || size(J, 1) != ctx.nvar || ctx.colptr != J.colptr || ctx.rowval != J.rowval
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
    resize!(ctx.sol, length(rhs))
    ldiv!(ctx.sol, ctx.fact, rhs)
    all(isfinite, ctx.sol) || error("the GPU solve returned a non-finite step (singular or ill-conditioned Jacobian in $(ctx.precision))")
    return ctx.sol
  catch e
    e isa InterruptException && rethrow(e)
    ctx.fallback_to_cpu || rethrow(e)
    msg = sprint(showerror, e)
    ctx.fallback_reason = string("GPU error, the rest of the solve ran on the CPU: ", first(msg, 300))
    @warn ctx.fallback_reason
    ctx.fact = nothing
    ctx.fallback_count += 1
    ctx.cpu = UmfpackReuseNewtonContext()
    return solve_newton_factorized!(ctx.cpu, J, rhs; pattern_changed = true)
  end
end

_newton_context_backend(ctx::AbstractGpuNewtonContext) = ctx.cpu === nothing ? Symbol(:cudss_, ctx.precision) : :umfpack_reuse

# What a context adds to the solver status and the performance profile
# about the compute backend: nothing for a plain CPU context, the decision
# for a CPU context that stands in for a requested GPU, and backend,
# precision, cuDSS mode, reason, transfer time and GPU time for a GPU
# context.
_newton_context_extras(::Any) = NamedTuple()
_newton_context_extras(ctx::UmfpackReuseNewtonContext) = isempty(ctx.backend_note) ? NamedTuple() : (compute_backend = :cpu, compute_backend_reason = ctx.backend_note)
function _newton_context_extras(ctx::AbstractGpuNewtonContext)
  return (
    compute_backend = ctx.cpu === nothing ? :gpu : :gpu_then_cpu,
    gpu_precision = ctx.precision,
    gpu_cudss_mode = ctx.fast_mode ? :fast : :deterministic,
    compute_backend_reason = isempty(ctx.fallback_reason) ? ctx.note : ctx.fallback_reason,
    gpu_transfer_s = ctx.transfer_s,
    gpu_solve_s = ctx.gpu_s,
  )
end

# GPU or CPU for a solve of `nbus` buses under `cb`, with the reason
function _compute_backend_choice(cb::ComputeBackendConfig, nbus::Integer)::NamedTuple{(:use_gpu, :reason),Tuple{Bool,String}}
  st = gpu_status()
  st.status === :available || return (use_gpu = false, reason = string("GPU not used: ", st.reason))
  cb.mode === :gpu && return (use_gpu = true, reason = string("GPU requested: ", st.device))
  advice = gpu_advice(nbus; precision = cb.gpu_precision)
  advice.decision === :gpu_faster && return (use_gpu = true, reason = string("auto: GPU expected faster by ", round(advice.factor; sigdigits = 2), " (", advice.run, ")"))
  advice.decision === :gpu_slower && return (use_gpu = false, reason = string("auto: GPU expected slower (factor ", round(advice.factor; sigdigits = 2), ", ", advice.run, ")"))
  return (use_gpu = false, reason = string("auto: no GPU measurement for ", nbus, " buses, CPU"))
end

"""
    newton_linear_solver_context(linear_solver::Symbol; nbus = 0)

A fresh per-solve context for the Newton linear solve. With the compute
backend `cpu` (the default, see [`ComputeBackendConfig`](@ref)) it maps
`power_flow.linear_solver`: `:umfpack_reuse` gives an
[`UmfpackReuseNewtonContext`](@ref), `:umfpack` gives `nothing` (the direct
`solve_linear` path). With `gpu`, or `auto` when the advisor expects the GPU
to be faster for `nbus` buses, it is the cuDSS context of the CUDA
extension. When the GPU is not available the solve runs on an
`UmfpackReuseNewtonContext` that records the reason, or, with
`gpu_fallback_to_cpu = false` and mode `gpu`, an `ArgumentError` names it.
"""
function newton_linear_solver_context(linear_solver::Symbol; nbus::Integer = 0)
  cb = _COMPUTE_BACKEND[]
  cb.mode === :cpu && return _cpu_newton_linear_solver_context(linear_solver)
  choice = _compute_backend_choice(cb, nbus)
  choice.use_gpu && return _GPU_BACKEND[].context(cb.gpu_precision, cb.gpu_fallback_to_cpu, cb.gpu_fast_mode, choice.reason)
  (cb.mode === :gpu && !cb.gpu_fallback_to_cpu) && throw(ArgumentError(string("power_flow.compute_backend.mode = gpu with gpu_fallback_to_cpu = false: ", choice.reason)))
  # the CPU stands in for the requested GPU (or auto chose it): the reuse
  # context carries the reason into the solver status
  # once per process and reason (a run starts several solves, the Web UI
  # many runs); every run names it in its result header ("Compute backend")
  cb.mode === :gpu && @warn choice.reason maxlog = 1 _id = Symbol(choice.reason)
  ctx = UmfpackReuseNewtonContext()
  ctx.backend_note = choice.reason
  return ctx
end

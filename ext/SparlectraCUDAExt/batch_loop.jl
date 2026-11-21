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
# file: ext/SparlectraCUDAExt/batch_loop.jl
# purpose: the injection batch on the GPU (runInjectionBatch with backend
#          :gpu): N scenarios of one network, one Ybus map, one cuDSS
#          analysis, the Newton iterations of all scenarios side by side in
#          device memory. Mismatch, Jacobian and voltage update are the
#          batch kernels of newton_kernels.jl, one launch per phase over all
#          active scenarios. Per iteration the host reads one download, the
#          mismatch norm of every scenario, and retires the converged and
#          the failed ones. Q limits are off, so the bus types never change.
#          Three linear-solve strategies:
#            :sequential  one cuDSS solver, the analysis of the start
#                         Jacobian reused for every scenario: per scenario a
#                         refactorization and a solve on its slice of the
#                         value vector (pointer update, no copy).
#                         Deterministic mode unless gpu_fast_mode.
#            :nonuniform  the cuDSS non-uniform batch (one matrix per
#                         scenario) on the Jacobians with their columns
#                         permuted by the matching of the start Jacobian;
#                         refactorizes every scenario every iteration (no
#                         mask in cuDSS 1.0.0). Deterministic mode unless
#                         gpu_fast_mode.
#            :uniform     the cuDSS uniform batch (ubatch_size N) on the
#                         same permuted Jacobians, converged scenarios
#                         masked out (or, with uniform_mask off, solved
#                         along every iteration). Fast mode only: cuDSS
#                         1.0.0 refuses the deterministic mode with
#                         ubatch_size > 1.
#          A device error (out of memory, cuDSS) with fallback_to_cpu on is
#          returned as (error = message,), so the core runs the batch on the
#          CPU.
#          The batch runs without matching (cuDSS refuses it there), and the
#          as-built Jacobian has zero diagonal entries, hence the column
#          permutation. A permuted system A[:, q] y = b keeps b and gives
#          x[q] = y.

# --- column permutation of the batch strategies ---------------------------------

# The CSR structure of J[:, q] and, per position of the permuted value
# vector, the position in the as-built one. q comes from the cuDSS matching
# of the start Jacobian (single solver, default matching).
struct _ColumnPermutation
  q::Vector{Int32}
  rowptr::Vector{Int32}
  colval::Vector{Int32}
  src::Vector{Int32}
end

function _matching_permutation(dmap::DeviceNewtonMap, nz1::CuVector{T}, nvar::Int) where {T}
  s = CudssSolver(dmap.jrowptr, dmap.jcol, nz1, "G", 'F')
  x = CUDA.zeros(T, nvar)
  b = CUDA.ones(T, nvar)
  cudss("analysis", s, x, b)
  pm = CUDA.zeros(Int32, nvar)
  cudss_set(s, "perm_matching", pm)
  cudss_get(s, "perm_matching")
  q = Array(pm)
  isperm(q) || error("cuDSS perm_matching is not a permutation of 1:$(nvar)")
  h = dmap.host
  qinv = invperm(q)
  rowptr = copy(h.jrowptr)
  colval = similar(h.jcol)
  src = Vector{Int32}(undef, length(h.jcol))
  order = Int[]
  for r = 1:nvar
    lo = h.jrowptr[r]
    hi = h.jrowptr[r+1] - 1
    resize!(order, hi - lo + 1)
    order .= lo:hi
    sort!(order; by = p -> qinv[h.jcol[p]])
    for (k, p) in enumerate(order)
      colval[lo+k-1] = qinv[h.jcol[p]]
      src[lo+k-1] = p
    end
  end
  return _ColumnPermutation(Int32.(q), rowptr, colval, src)
end

# nzp[k, s] = nz[(s-1) nnzJ + src[k]] for the active scenarios
function _k_gather_permuted!(nzp, nz, src, active, nnzJ)
  s = blockIdx().y
  @inbounds active[s] == 0x00 && return nothing
  k = (blockIdx().x - 1) * blockDim().x + threadIdx().x
  k > length(src) && return nothing
  @inbounds nzp[k, s] = nz[(Int64(s) - 1) * Int64(nnzJ) + Int64(src[k])]
  return nothing
end

# sol[q[k], s] = solp[k, s] for the active scenarios
function _k_scatter_solution!(sol, solp, q, active)
  s = blockIdx().y
  @inbounds active[s] == 0x00 && return nothing
  k = (blockIdx().x - 1) * blockDim().x + threadIdx().x
  k > length(q) && return nothing
  @inbounds sol[q[k], s] = solp[k, s]
  return nothing
end

# Voltage-dependent bus shunts, one thread per bus and scenario: the
# specified injection S = Spro - abs(V)^2 * ysh and the chain terms
# dP, dQ = -2 vm_safe * ysh, as buildControlledSVec (vm_safe guards abs(V)
# below vm_eps). A bus without such a shunt keeps S = Spro and zero terms.
# Generic in the real type of V (Float32 in phase 1 of the two-phase batch).
function _k_vdep_shunt!(S, dP, dQ, Spro, ysh, V, active, vm_eps)
  s = blockIdx().y
  @inbounds active[s] == 0x00 && return nothing
  i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
  i > size(V, 1) && return nothing
  R = real(eltype(V))
  y = @inbounds ysh[i]
  sp = @inbounds Spro[i, s]
  if y == 0
    @inbounds S[i, s] = sp
    @inbounds dP[i, s] = zero(R)
    @inbounds dQ[i, s] = zero(R)
    return nothing
  end
  v = @inbounds V[i, s]
  vm = _hypot_cr(real(v), imag(v))
  vm_safe = vm > vm_eps ? vm : vm_eps
  vm2 = _m(vm, vm)
  @inbounds S[i, s] = Complex{R}(_sb(real(sp), _m(vm2, real(y))), _sb(imag(sp), _m(vm2, imag(y))))
  g = _m(-R(2), vm_safe)
  @inbounds dP[i, s] = _m(g, real(y))
  @inbounds dQ[i, s] = _m(g, imag(y))
  return nothing
end

# the prosumer part of S and the shunt admittances on the device
struct _VdepShunts
  Spro::CuMatrix{ComplexF64}
  ysh::CuVector{ComplexF64}
end

function _launch_vdep!(st::DeviceNewtonState, vd::_VdepShunts)
  @cuda threads = _THREADS blocks = (_blocks(st.map.host.n), st.N) _k_vdep_shunt!(st.S, st.dP, st.dQ, vd.Spro, vd.ysh, st.V, st.active, 1e-9)
  return nothing
end

# --- linear-solve strategies --------------------------------------------------

abstract type _BatchLinear end

# :sequential: one solver whose matrix, right-hand side and solution point
# at the slice of the scenario being solved
mutable struct _SequentialLinear{T} <: _BatchLinear
  solver::CudssSolver{T}
  X::CudssMatrix{T}
  B::CudssMatrix{T}
  analyses::Int
  refactorizations::Int
end

# :nonuniform and :uniform: the permuted values, the permuted solution and
# the batch solver
mutable struct _PermutedBatchLinear{T,S} <: _BatchLinear
  perm::_ColumnPermutation
  rowptr::CuVector{Int32}
  colval::CuVector{Int32}
  src::CuVector{Int32}
  q::CuVector{Int32}
  nzp::CuMatrix{T}
  solp::CuMatrix{T}
  solver::S
  X::Any
  B::Any
  uniform::Bool
  mask::Bool
  analyses::Int
  refactorizations::Int
end

_slice(nz::CuVector, s::Int, nnzJ::Int) = view(nz, ((s-1)*nnzJ+1):(s*nnzJ))

function _sequential_linear(st::DeviceNewtonState{T}, fast_mode::Bool) where {T}
  h = st.map.host
  nnzJ = h.nnzJ
  solver = _cudss_configure!(CudssSolver(st.map.jrowptr, st.map.jcol, _slice(st.nz, 1, nnzJ), "G", 'F'), fast_mode)
  # no pivot perturbation: a zero pivot gives a non-finite step, which the
  # next mismatch norm retires as a failed scenario (as the device loop)
  cudss_set(solver, "pivot_epsilon", zero(T))
  X = CudssMatrix(view(st.sol, :, 1))
  B = CudssMatrix(view(st.rhs, :, 1))
  cudss("analysis", solver, X, B)
  cudss("factorization", solver, X, B)
  return _SequentialLinear{T}(solver, X, B, 1, 0)
end

function _permuted_linear(st::DeviceNewtonState{T}, strategy::Symbol, fast_mode::Bool, use_mask::Bool) where {T}
  h = st.map.host
  N = st.N
  nvar = h.nvar
  perm = _matching_permutation(st.map, _slice(st.nz, 1, h.nnzJ), nvar)
  rowptr = CuVector(perm.rowptr)
  colval = CuVector(perm.colval)
  src = CuVector(perm.src)
  q = CuVector(perm.q)
  nzp = CUDA.zeros(T, h.nnzJ, N)
  solp = CUDA.zeros(T, nvar, N)
  all_active = CUDA.ones(UInt8, N)
  @cuda threads = _THREADS blocks = (_blocks(h.nnzJ), N) _k_gather_permuted!(nzp, st.nz, src, all_active, Int64(h.nnzJ))
  CUDA.unsafe_free!(all_active)
  if strategy === :uniform
    solver = CudssSolver(rowptr, colval, vec(nzp), "G", 'F')
    cudss_set(solver, "ubatch_size", N)
    cudss_set(solver, "matching_alg", "none")
    # the deterministic mode is refused with ubatch_size > 1 (cuDSS 1.0.0)
    cudss_set(solver, "pivot_epsilon", zero(T))
    # the mask is part of the solver data from the analysis on; later
    # iterations only change its entries. Without it every iteration
    # refactorizes and solves all N.
    use_mask && cudss_set(solver, "ubatch_mask", ones(Cint, N))
    X = CudssMatrix(T, nvar; nbatch = N)
    cudss_update(X, vec(solp))
    B = CudssMatrix(T, nvar; nbatch = N)
    cudss_update(B, vec(st.rhs))
  else
    A = [CuSparseMatrixCSR{T,Int32}(rowptr, colval, view(nzp, :, s), (nvar, nvar)) for s = 1:N]
    solver = _cudss_configure!(CudssBatchedSolver(A, "G", 'F'), fast_mode)
    cudss_set(solver, "matching_alg", "none")
    cudss_set(solver, "pivot_epsilon", zero(T))
    X = [view(solp, :, s) for s = 1:N]
    B = [view(st.rhs, :, s) for s = 1:N]
  end
  cudss("analysis", solver, X, B)
  cudss("factorization", solver, X, B)
  return _PermutedBatchLinear{T,typeof(solver)}(perm, rowptr, colval, src, q, nzp, solp, solver, X, B, strategy === :uniform, strategy === :uniform && use_mask, 1, 0)
end

# The Newton solves of iteration `it` for the active scenarios (host list
# `act`, device flags st.active). In iteration 1 every scenario has the
# start Jacobian of the analysis, so no refactorization runs.
function _batch_solve!(L::_SequentialLinear{T}, st::DeviceNewtonState{T}, act::Vector{Int}, it::Int, et::_EventTimer) where {T}
  nnzJ = st.map.host.nnzJ
  for s in act
    _timed_phase(et, :refactor) do
      cudss_update(L.X, view(st.sol, :, s))
      cudss_update(L.B, view(st.rhs, :, s))
      if it > 1
        cudss_update(L.solver, st.map.jrowptr, st.map.jcol, _slice(st.nz, s, nnzJ))
        cudss("refactorization", L.solver, L.X, L.B)
        L.refactorizations += 1
      end
    end
    _timed_phase(() -> cudss("solve", L.solver, L.X, L.B), et, :solve)
  end
  return nothing
end

function _batch_solve!(L::_PermutedBatchLinear{T}, st::DeviceNewtonState{T}, act::Vector{Int}, it::Int, et::_EventTimer) where {T}
  h = st.map.host
  N = st.N
  if L.mask
    mask = zeros(Cint, N)
    mask[act] .= 1
    cudss_set(L.solver, "ubatch_mask", mask)
  end
  if it > 1
    _timed_phase(et, :refactor) do
      @cuda threads = _THREADS blocks = (_blocks(h.nnzJ), N) _k_gather_permuted!(L.nzp, st.nz, L.src, st.active, Int64(h.nnzJ))
      cudss("refactorization", L.solver, L.X, L.B)
    end
    L.refactorizations += L.mask ? length(act) : N
  end
  _timed_phase(et, :solve) do
    cudss("solve", L.solver, L.X, L.B)
    @cuda threads = _THREADS blocks = (_blocks(h.nvar), N) _k_scatter_solution!(st.sol, L.solp, L.q, st.active)
  end
  return nothing
end

_linear_free!(::_SequentialLinear) = nothing
function _linear_free!(L::_PermutedBatchLinear)
  for a in (L.rowptr, L.colval, L.src, L.q, L.nzp, L.solp)
    CUDA.unsafe_free!(a)
  end
  return nothing
end

# --- the batch loop ---------------------------------------------------------------

# the CUDSS.jl spelling of a cuDSS algorithm value: 0 is "default" for
# reordering_alg and "none" for matching_alg (CUDSS_MATCHING_ALG_NONE; the
# matching default of cuDSS 1.0 is 6, "algo6"), so the batch modes, which
# set matching_alg to "none", report it as set
_alg_name(a::Integer, zero_name::String = "default") = a == 0 ? zero_name : string("algo", a)

function _gpu_injection_batch(input::Sparlectra.InjectionBatchGpuInput)
  T = input.precision === :float32_mixed ? Float32 : Float64
  try
    input.two_phase === :off || return _gpu_two_phase_batch(input)
    return _gpu_injection_batch(input, T)
  catch e
    # the device errors the CPU fallback is meant for; the device memory of
    # the batch is freed in the finally block of the batch loop
    (_is_gpu_error(e) && input.fallback_to_cpu) || rethrow(e)
    return (error = first(first(split(sprint(showerror, e), '\n')), 300),)
  end
end

function _gpu_injection_batch(input::Sparlectra.InjectionBatchGpuInput, ::Type{T}) where {T}
  t_total = time_ns()
  n = length(input.V0)
  N = size(input.S, 2)
  strategy = input.strategy
  fast_mode = input.fast_mode || strategy === :uniform
  t0 = time_ns()
  hmap = HostNewtonMap(input.Ybus, input.slack_idx)
  map_s = (time_ns() - t0) / 1.0e9
  t0 = time_ns()
  st = DeviceNewtonState{T}(DeviceNewtonMap(hmap), N)
  L = nothing
  vd = any(!iszero, input.ysh) ? _VdepShunts(CuMatrix(input.S), CuVector(input.ysh)) : nothing
  try
    copyto!(st.V, repeat(input.V0, 1, N))
    copyto!(st.S, input.S)
    copyto!(st.Vset, repeat(input.Vset, 1, N))
    copyto!(st.btype, repeat(_bus_type_codes!(Vector{UInt8}(undef, n), input.bus_types), 1, N))
    # one Newton step of full length, as runpf! with damp = 1.0
    fill!(st.alpha, 1.0)
    # J at the start iterate, equal for all scenarios (one V0; the shunt
    # chain terms depend on V only): the analysis runs on scenario 1
    vd === nothing || _launch_vdep!(st, vd)
    _launch_mismatch!(st)
    _launch_jacobian!(st, 1e-9)
    CUDA.synchronize()
    upload_s = (time_ns() - t0) / 1.0e9
    t0 = time_ns()
    L = strategy === :sequential ? _sequential_linear(st, fast_mode) : _permuted_linear(st, strategy, fast_mode, input.uniform_mask)
    CUDA.synchronize()
    analysis_s = (time_ns() - t0) / 1.0e9
    reordering = _alg_name(cudss_get(L.solver, "reordering_alg"))
    matching = _alg_name(cudss_get(L.solver, "matching_alg"), "none")
    # device memory in use after the analysis (all processes on the device)
    used_bytes = CUDA.total_memory() - CUDA.free_memory()
    used_peak = used_bytes

    t_loop = time_ns()
    et = _EventTimer()
    polar = input.newton_update === :polar
    stat_h = zeros(Float64, N)
    active_h = ones(UInt8, N)
    converged = falses(N)
    iterations = zeros(Int, N)
    final_mismatch = fill(NaN, N)
    act = collect(1:N)
    d2h = 0
    h2d = 0
    for it = 1:input.maxiter
      vd === nothing || _timed_phase(() -> _launch_vdep!(st, vd), et, :injections)
      _timed_phase(() -> _launch_mismatch!(st), et, :mismatch)
      # the one download of the iteration (synchronizes)
      copyto!(stat_h, st.stat)
      d2h += 1
      retired = false
      for s in act
        m = stat_h[s]
        iterations[s] = it
        final_mismatch[s] = m
        if !isfinite(m)
          active_h[s] = 0x00
          retired = true
        elseif m <= input.tol
          converged[s] = true
          active_h[s] = 0x00
          retired = true
        end
      end
      if retired
        filter!(s -> active_h[s] == 0x01, act)
        copyto!(st.active, active_h)
        h2d += 1
      end
      (isempty(act) || it == input.maxiter) && break
      it > 1 && _timed_phase(() -> _launch_jacobian!(st, 1e-9), et, :build)
      _batch_solve!(L, st, act, it, et)
      # device memory in use after the solves of the iteration (a driver
      # query, no synchronization); its maximum is the peak of the loop
      used_peak = max(used_peak, CUDA.total_memory() - CUDA.free_memory())
      _timed_phase(() -> _launch_update!(st, polar), et, :update)
    end
    V = Array(st.V)
    CUDA.synchronize()
    loop_s = (time_ns() - t_loop) / 1.0e9
    stats = Dict{Symbol,Any}(
      :gpu_strategy => strategy,
      :gpu_precision => input.precision,
      :gpu_cudss_mode => fast_mode ? :fast : :deterministic,
      :gpu_map_s => map_s,
      :gpu_upload_s => upload_s,
      :gpu_analysis_s => analysis_s,
      :gpu_loop_s => loop_s,
      :gpu_mismatch_s => _phase_seconds(et, :mismatch),
      :gpu_injections_s => _phase_seconds(et, :injections),
      :gpu_vdep_shunts => vd !== nothing,
      :gpu_build_s => _phase_seconds(et, :build),
      :gpu_refactor_s => _phase_seconds(et, :refactor),
      :gpu_solve_s => _phase_seconds(et, :solve),
      :gpu_update_s => _phase_seconds(et, :update),
      :gpu_total_s => (time_ns() - t_total) / 1.0e9,
      :gpu_analyses => L.analyses,
      :gpu_refactorizations => L.refactorizations,
      :gpu_d2h_count => d2h,
      :gpu_h2d_count => h2d,
      :gpu_used_bytes_after_analysis => used_bytes,
      :gpu_used_bytes_peak => used_peak,
      :gpu_device => CUDA.name(CUDA.device()),
      :gpu_cudss_version => string(CUDSS.version()),
      :gpu_reordering_alg => reordering,
      :gpu_matching_alg => matching,
      :gpu_uniform_mask => L isa _PermutedBatchLinear && L.mask,
    )
    return (V = V, converged = Vector{Bool}(converged), iterations = iterations, final_mismatch = final_mismatch, stats = stats)
  finally
    L === nothing || _linear_free!(L)
    vd === nothing || (CUDA.unsafe_free!(vd.Spro); CUDA.unsafe_free!(vd.ysh))
    _free!(st)
  end
end

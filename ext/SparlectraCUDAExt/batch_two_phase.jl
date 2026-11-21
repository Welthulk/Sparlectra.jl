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
# file: ext/SparlectraCUDAExt/batch_two_phase.jl
# purpose: the two-phase variants of the injection batch on the GPU
#          (runInjectionBatch with gpu_two_phase; with N = 1 the single
#          run).
#          Phase 1: V, mismatch, Jacobian and LU in Float32 on the device,
#          until the Float32 mismatch norm is at most switch_tol, or on
#          stagnation (two steps in a row with rho = norm_k / norm_(k-1)
#          above 0.5).
#          Phase 2:
#            :g5   V and mismatch in Float64 on the device, steps with the
#                  last Float32 LU of phase 1 (no refactorization) down to
#                  the batch tolerance. After 5 steps without convergence
#                  one Float32 refactorization at the current iterate,
#                  counted; then on with that LU.
#            :g5r  as :g5, but by the rate: rho above 0.5 gives one counted
#                  Float32 refactorization; rho above 0.5 again after it
#                  gives Float64 (Jacobian and LU) to the end, counted
#                  separately.
#            :g6   none on the device: the phase-1 iterate goes back to the
#                  core, which finishes with a Float64 Newton on the CPU.
#          rho is logged per step in both phases. The scenarios run one
#          after another, as in the sequential strategy: one cuDSS solver
#          holds the factors of one scenario, and phase 2 needs the last
#          phase-1 LU of its own scenario. One analysis, on the Float32
#          Jacobian at the start iterate (equal for every scenario), serves
#          all of them; the work arrays hold one scenario and are reused.

const _RHO_LIMIT = 0.5
const _PHASE2_CHORD_STEPS = 5

# the Float32 arrays of phase 1, one scenario
struct _Phase1Arrays
  yval::CuVector{ComplexF32}
  V0::CuMatrix{ComplexF32}
  V::CuMatrix{ComplexF32}
  I::CuMatrix{ComplexF32}
  S::CuMatrix{ComplexF32}
  Spro::CuMatrix{ComplexF32}
  ysh::CuVector{ComplexF32}
  Vset::CuMatrix{Float32}
  dP::CuMatrix{Float32}
  dQ::CuMatrix{Float32}
  F::CuMatrix{Float32}
  partial::CuMatrix{Float32}
  stat::CuVector{Float32}
  alpha::CuVector{Float32}
end

function _Phase1Arrays(hmap::HostNewtonMap, V0::Vector{ComplexF64}, Vset::Vector{Float64}, ysh::Vector{ComplexF64})
  n = hmap.n
  nb = _blocks(n)
  v0 = CuMatrix(reshape(ComplexF32.(V0), n, 1))
  return _Phase1Arrays(CuVector(ComplexF32.(hmap.yval)), v0, copy(v0), CUDA.zeros(ComplexF32, n, 1), CUDA.zeros(ComplexF32, n, 1),
    CUDA.zeros(ComplexF32, n, 1), CuVector(ComplexF32.(ysh)), CuMatrix(reshape(Float32.(Vset), n, 1)), CUDA.zeros(Float32, n, 1),
    CUDA.zeros(Float32, n, 1), CUDA.zeros(Float32, hmap.nvar, 1), CUDA.zeros(Float32, nb, 1), CUDA.zeros(Float32, 1), CUDA.ones(Float32, 1))
end

_free!(a::_Phase1Arrays) = foreach(f -> CUDA.unsafe_free!(getfield(a, f)), fieldnames(_Phase1Arrays))

# the Float64 Newton of :g5r after its second slow step, one scenario
mutable struct _Float64Linear
  nz::CuVector{Float64}
  rhs::CuMatrix{Float64}
  sol::CuMatrix{Float64}
  solver::Any
  X::Any
  B::Any
end

# mismatch, Jacobian and update launches on explicit arrays (one scenario),
# with the kernels of newton_kernels.jl
function _tp_mismatch!(dmap::DeviceNewtonMap, yval, V, I, S, Vset, btype, F, rhs, partial, stat, active)
  nb = _blocks(dmap.host.n)
  @cuda threads = _THREADS blocks = (nb, 1) _k_current_mismatch!(I, F, rhs, partial, dmap.yrowptr, dmap.ycol, yval, V, S, Vset, btype, dmap.pos, active)
  @cuda threads = _THREADS blocks = 1 _k_reduce_max!(stat, partial, Int32(nb), active)
  return nothing
end

function _tp_jacobian!(dmap::DeviceNewtonMap, nz, yval, V, I, dP, dQ, btype, active, vm_eps)
  h = dmap.host
  nnzJ = Int64(h.nnzJ)
  @cuda threads = _THREADS blocks = (_blocks(h.nnzY), 1) _k_jacobian_offdiag!(nz, dmap.yrow, dmap.jpos, dmap.mcnt, yval, V, btype, active, nnzJ)
  @cuda threads = _THREADS blocks = (_blocks(h.n - 1), 1) _k_jacobian_diag!(nz, I, V, yval, dmap.ydiag, dmap.dpos, dmap.mcnt, dP, dQ, btype, dmap.nonslack, active, vm_eps, nnzJ)
  return nothing
end

function _tp_update!(dmap::DeviceNewtonMap, V, sol, alpha, active, polar::Bool)
  @cuda threads = _THREADS blocks = (_blocks(dmap.host.n - 1), 1) _k_update!(V, sol, alpha, dmap.nonslack, active, polar)
  return nothing
end

function _tp_vdep!(S, dP, dQ, Spro, ysh, V, active, vm_eps)
  @cuda threads = _THREADS blocks = (_blocks(size(V, 1)), 1) _k_vdep_shunt!(S, dP, dQ, Spro, ysh, V, active, vm_eps)
  return nothing
end

_rate(norms::Vector{Float64}) = length(norms) < 2 ? NaN : norms[end] / norms[end-1]

function _gpu_two_phase_batch(input::Sparlectra.InjectionBatchGpuInput)
  t_total = time_ns()
  variant = input.two_phase
  n = length(input.V0)
  N = size(input.S, 2)
  fast_mode = input.fast_mode
  polar = input.newton_update === :polar
  vm_eps = 1e-9
  t0 = time_ns()
  hmap = HostNewtonMap(input.Ybus, input.slack_idx)
  map_s = (time_ns() - t0) / 1.0e9
  t0 = time_ns()
  # phase-2 arrays (V and mismatch Float64, Jacobian, rhs and solution
  # Float32), one scenario; its map serves phase 1 too
  st = DeviceNewtonState{Float32}(DeviceNewtonMap(hmap), 1)
  dmap = st.map
  p1 = nothing
  L = nothing
  L64 = nothing
  has_vd = any(!iszero, input.ysh)
  Smat = CuMatrix(input.S)
  Spro64 = CUDA.zeros(ComplexF64, n, 1)
  V0slack = CuVector([input.V0[input.slack_idx]])
  ysh64 = CuVector(input.ysh)
  try
    p1 = _Phase1Arrays(hmap, input.V0, input.Vset, input.ysh)
    copyto!(st.Vset, reshape(input.Vset, n, 1))
    copyto!(st.btype, reshape(_bus_type_codes!(Vector{UInt8}(undef, n), input.bus_types), n, 1))
    fill!(st.alpha, 1.0)
    # the Float32 Jacobian at the start iterate, equal for every scenario
    # (S enters only the mismatch): analysis and first factorization
    copyto!(Spro64, 1, Smat, 1, n)
    p1.Spro .= ComplexF32.(Spro64)
    has_vd ? _tp_vdep!(p1.S, p1.dP, p1.dQ, p1.Spro, p1.ysh, p1.V, st.active, Float32(vm_eps)) : copyto!(p1.S, p1.Spro)
    _tp_mismatch!(dmap, p1.yval, p1.V, p1.I, p1.S, p1.Vset, st.btype, p1.F, st.rhs, p1.partial, p1.stat, st.active)
    _tp_jacobian!(dmap, st.nz, p1.yval, p1.V, p1.I, p1.dP, p1.dQ, st.btype, st.active, Float32(vm_eps))
    CUDA.synchronize()
    upload_s = (time_ns() - t0) / 1.0e9
    t0 = time_ns()
    L = _sequential_linear(st, fast_mode)
    CUDA.synchronize()
    analysis_s = (time_ns() - t0) / 1.0e9
    analyses = 1
    reordering = _alg_name(cudss_get(L.solver, "reordering_alg"))
    matching = _alg_name(cudss_get(L.solver, "matching_alg"), "none")
    used_after = CUDA.total_memory() - CUDA.free_memory()
    used_peak = used_after

    V = Matrix{ComplexF64}(undef, n, N)
    converged = falses(N)
    iterations = zeros(Int, N)
    final_mismatch = fill(NaN, N)
    p1_steps = zeros(Int, N)
    p2_steps = zeros(Int, N)
    reason = fill(:none, N)
    p2_refactor = zeros(Int, N)
    to_f64 = falses(N)
    rho_log = [Float64[] for _ = 1:N]
    p1_s = 0.0
    p2_s = 0.0
    refactor_count = 0
    f64_refactor_count = 0
    norm32 = zeros(Float32, 1)
    norm64 = zeros(Float64, 1)
    vcol = Vector{ComplexF64}(undef, n)
    t_loop = time_ns()
    for s = 1:N
      ts = time_ns()
      # --- phase 1, all Float32 ---
      copyto!(p1.V, p1.V0)
      copyto!(Spro64, 1, Smat, (s - 1) * n + 1, n)
      p1.Spro .= ComplexF32.(Spro64)
      has_vd || copyto!(p1.S, p1.Spro)
      norms = Float64[]
      while true
        has_vd && _tp_vdep!(p1.S, p1.dP, p1.dQ, p1.Spro, p1.ysh, p1.V, st.active, Float32(vm_eps))
        _tp_mismatch!(dmap, p1.yval, p1.V, p1.I, p1.S, p1.Vset, st.btype, p1.F, st.rhs, p1.partial, p1.stat, st.active)
        copyto!(norm32, p1.stat)
        m = Float64(norm32[1])
        push!(norms, m)
        push!(rho_log[s], _rate(norms))
        if !isfinite(m)
          reason[s] = :failed
        elseif m <= input.switch_tol
          reason[s] = :tol
        elseif length(norms) >= 3 && norms[end] > _RHO_LIMIT * norms[end-1] && norms[end-1] > _RHO_LIMIT * norms[end-2]
          reason[s] = :stagnation
        elseif length(norms) >= input.maxiter
          reason[s] = :maxiter
        end
        reason[s] === :none || break
        # the factors of the analysis are the start Jacobian, which only the
        # first step of scenario 1 finds in the solver
        if !(s == 1 && length(norms) == 1)
          _tp_jacobian!(dmap, st.nz, p1.yval, p1.V, p1.I, p1.dP, p1.dQ, st.btype, st.active, Float32(vm_eps))
          cudss("refactorization", L.solver, L.X, L.B)
          L.refactorizations += 1
        end
        cudss("solve", L.solver, L.X, L.B)
        _tp_update!(dmap, p1.V, st.sol, p1.alpha, st.active, polar)
        p1_steps[s] += 1
      end
      evals = length(norms)
      t1 = time_ns()
      p1_s += (t1 - ts) / 1.0e9
      if variant === :g6 || reason[s] in (:failed, :maxiter)
        # g6 finishes on the CPU; a failed or exhausted phase 1 ends here
        copyto!(vcol, ComplexF64.(Array(p1.V)))
        # the slack keeps its Float64 value (phase 1 holds it rounded)
        vcol[input.slack_idx] = input.V0[input.slack_idx]
        V[:, s] .= vcol
        iterations[s] = evals
        final_mismatch[s] = norms[end]
        continue
      end
      # --- phase 2, V and mismatch Float64 ---
      st.V .= ComplexF64.(p1.V)
      # the slack is never updated: its Float64 value, not the rounded one
      # of phase 1 (else phase 2 solves for a slack voltage off by up to
      # one Float32 rounding)
      copyto!(st.V, input.slack_idx, V0slack, 1, 1)
      has_vd || copyto!(st.S, Spro64)
      norms2 = Float64[]
      since_refactor = 0
      while true
        has_vd && _tp_vdep!(st.S, st.dP, st.dQ, Spro64, ysh64, st.V, st.active, vm_eps)
        rhs = to_f64[s] ? L64.rhs : st.rhs
        _tp_mismatch!(dmap, dmap.yval, st.V, st.I, st.S, st.Vset, st.btype, st.F, rhs, st.partial, st.stat, st.active)
        copyto!(norm64, st.stat)
        m = norm64[1]
        push!(norms2, m)
        rho = _rate(norms2)
        push!(rho_log[s], rho)
        evals += 1
        final_mismatch[s] = m
        isfinite(m) || break
        if m <= input.tol
          converged[s] = true
          break
        end
        evals >= input.maxiter && break
        if variant === :g5
          if since_refactor == _PHASE2_CHORD_STEPS && p2_refactor[s] == 0
            _tp_jacobian!(dmap, st.nz, dmap.yval, st.V, st.I, st.dP, st.dQ, st.btype, st.active, vm_eps)
            cudss("refactorization", L.solver, L.X, L.B)
            L.refactorizations += 1
            p2_refactor[s] += 1
            refactor_count += 1
            since_refactor = 0
          end
        elseif !to_f64[s] && rho > _RHO_LIMIT
          # :g5r: the first slow step refactorizes in Float32, the next one
          # goes to Float64
          if p2_refactor[s] == 0
            _tp_jacobian!(dmap, st.nz, dmap.yval, st.V, st.I, st.dP, st.dQ, st.btype, st.active, vm_eps)
            cudss("refactorization", L.solver, L.X, L.B)
            L.refactorizations += 1
            p2_refactor[s] += 1
            refactor_count += 1
          else
            to_f64[s] = true
          end
        end
        if to_f64[s]
          if L64 === nothing
            L64 = _Float64Linear(CUDA.zeros(Float64, hmap.nnzJ), CUDA.zeros(Float64, hmap.nvar, 1), CUDA.zeros(Float64, hmap.nvar, 1), nothing, nothing, nothing)
          end
          # the mismatch of this iterate again, now into the Float64 rhs
          rhs === L64.rhs || _tp_mismatch!(dmap, dmap.yval, st.V, st.I, st.S, st.Vset, st.btype, st.F, L64.rhs, st.partial, st.stat, st.active)
          _tp_jacobian!(dmap, L64.nz, dmap.yval, st.V, st.I, st.dP, st.dQ, st.btype, st.active, vm_eps)
          if L64.solver === nothing
            L64.solver = _cudss_configure!(CudssSolver(dmap.jrowptr, dmap.jcol, L64.nz, "G", 'F'), fast_mode)
            cudss_set(L64.solver, "pivot_epsilon", 0.0)
            L64.X = CudssMatrix(view(L64.sol, :, 1))
            L64.B = CudssMatrix(view(L64.rhs, :, 1))
            cudss("analysis", L64.solver, L64.X, L64.B)
            cudss("factorization", L64.solver, L64.X, L64.B)
            analyses += 1
          else
            cudss("refactorization", L64.solver, L64.X, L64.B)
          end
          f64_refactor_count += 1
          cudss("solve", L64.solver, L64.X, L64.B)
          _tp_update!(dmap, st.V, L64.sol, st.alpha, st.active, polar)
        else
          cudss("solve", L.solver, L.X, L.B)
          _tp_update!(dmap, st.V, st.sol, st.alpha, st.active, polar)
        end
        p2_steps[s] += 1
        since_refactor += 1
      end
      iterations[s] = evals
      copyto!(vcol, view(Array(st.V), :, 1))
      V[:, s] .= vcol
      p2_s += (time_ns() - t1) / 1.0e9
      used_peak = max(used_peak, CUDA.total_memory() - CUDA.free_memory())
    end
    CUDA.synchronize()
    loop_s = (time_ns() - t_loop) / 1.0e9
    stats = Dict{Symbol,Any}(
      :gpu_strategy => :sequential,
      :gpu_precision => variant,
      :gpu_cudss_mode => fast_mode ? :fast : :deterministic,
      :gpu_two_phase => variant,
      :gpu_switch_tol => input.switch_tol,
      :gpu_map_s => map_s,
      :gpu_upload_s => upload_s,
      :gpu_analysis_s => analysis_s,
      :gpu_loop_s => loop_s,
      :gpu_phase1_s => p1_s,
      :gpu_phase2_s => p2_s,
      :gpu_total_s => (time_ns() - t_total) / 1.0e9,
      :gpu_analyses => analyses,
      :gpu_refactorizations => L.refactorizations + f64_refactor_count,
      :gpu_phase2_refactorizations => refactor_count,
      :gpu_f64_refactorizations => f64_refactor_count,
      :gpu_used_bytes_after_analysis => used_after,
      :gpu_used_bytes_peak => used_peak,
      :gpu_device => CUDA.name(CUDA.device()),
      :gpu_cudss_version => string(CUDSS.version()),
      :gpu_reordering_alg => reordering,
      :gpu_matching_alg => matching,
      :gpu_uniform_mask => false,
      :gpu_max_abs_y => maximum(abs, hmap.yval),
      :tp_phase1_steps => p1_steps,
      :tp_phase2_steps => p2_steps,
      :tp_switch_reason => reason,
      :tp_phase2_refactorizations => p2_refactor,
      :tp_float64 => Vector{Bool}(to_f64),
      :tp_rho => rho_log,
    )
    return (V = V, converged = Vector{Bool}(converged), iterations = iterations, final_mismatch = final_mismatch, stats = stats)
  finally
    p1 === nothing || _free!(p1)
    if L64 !== nothing
      foreach(CUDA.unsafe_free!, (L64.nz, L64.rhs, L64.sol))
    end
    foreach(CUDA.unsafe_free!, (Smat, Spro64, V0slack, ysh64))
    _free!(st)
  end
end

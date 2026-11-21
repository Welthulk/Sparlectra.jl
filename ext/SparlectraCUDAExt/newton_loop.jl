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
# file: ext/SparlectraCUDAExt/newton_loop.jl
# purpose: the device Newton loop (GPU study, stage 1b): mismatch, Jacobian
#          build, cuDSS refactorization and solve and the voltage update in
#          device memory, the Q-limit decision unchanged on the CPU. V stays
#          on the device; per iteration the host reads, in one download, the
#          mismatch norm, with Q limits on the values the Q-limit decision
#          reads (Q of the PV buses, Q and abs(V) of the release candidates)
#          and after a Q-limit switch the norm after the switch. Entry for
#          the solver: _gpu_newton; entry for the proof script:
#          _device_jacobian_mismatch.

# the GPU errors a fallback to the CPU is meant for; anything else (a bug in
# a kernel, InvalidIRError, KernelException) is rethrown
_is_gpu_error(e) = e isa CUDA.CuError || e isa CUDA.OutOfGPUMemoryError || e isa CUDSS.CUDSSError

# Device arrays of N scenarios, column s is scenario s. T is the type of the
# Jacobian values, the right-hand side and the solution (Float64 for G3,
# Float32 for G4); everything else is Float64.
struct DeviceNewtonState{T}
  map::DeviceNewtonMap
  N::Int
  V::CuMatrix{ComplexF64}
  I::CuMatrix{ComplexF64}
  S::CuMatrix{ComplexF64}
  Vset::CuMatrix{Float64}
  dP::CuMatrix{Float64}
  dQ::CuMatrix{Float64}
  btype::CuMatrix{UInt8}
  F::CuMatrix{Float64}
  partial::CuMatrix{Float64}
  # the mismatch norm per scenario
  stat::CuVector{Float64}
  active::CuVector{UInt8}
  recompute::CuVector{UInt8}
  alpha::CuVector{Float64}
  nz::CuVector{T}
  rhs::CuMatrix{T}
  sol::CuMatrix{T}
  # device copies of the iterate of the last finite mismatch (also the
  # iterate before a step not checked yet) and of the best one; the host
  # reads them only at the end of the solve or for a CPU redo
  Vlast::CuMatrix{ComplexF64}
  Vbest::CuMatrix{ComplexF64}
  # Q-limit exchange of a single scenario: requests (b: Q of bus b, -b:
  # abs(V_b)); the sample: the norm, the values of the requests, and after a
  # switch the norm after the switch behind them; the switch buffer
  req::CuVector{Int32}
  sample::CuVector{Float64}
  swbuf::CuVector{Float64}
end

# qcap: buses that can be PV during the solve (Q-limit exchange size)
function DeviceNewtonState{T}(dmap::DeviceNewtonMap, N::Int; qcap::Int = 0) where {T}
  h = dmap.host
  n = h.n
  nb = _blocks(n)
  cz() = CUDA.zeros(ComplexF64, n, N)
  return DeviceNewtonState{T}(
    dmap, N, cz(), cz(), cz(),
    CUDA.zeros(Float64, n, N), CUDA.zeros(Float64, n, N), CUDA.zeros(Float64, n, N),
    CUDA.zeros(UInt8, n, N), CUDA.zeros(Float64, h.nvar, N), CUDA.zeros(Float64, nb, N),
    CUDA.zeros(Float64, N), CUDA.ones(UInt8, N), CUDA.zeros(UInt8, N), CUDA.ones(Float64, N),
    CUDA.zeros(T, h.nnzJ * N), CUDA.zeros(T, h.nvar, N), CUDA.zeros(T, h.nvar, N),
    cz(), cz(),
    CUDA.zeros(Int32, max(2qcap, 1)), CUDA.zeros(Float64, 2 + 2qcap), CUDA.zeros(Float64, max(4qcap, 4)),
  )
end

function _free!(st::DeviceNewtonState)
  for a in (st.V, st.I, st.S, st.Vset, st.dP, st.dQ, st.btype, st.F, st.partial, st.stat, st.active, st.recompute, st.alpha, st.nz, st.rhs, st.sol, st.Vlast, st.Vbest, st.req, st.sample, st.swbuf)
    CUDA.unsafe_free!(a)
  end
  m = st.map
  for a in (m.yrowptr, m.ycol, m.yrow, m.yval, m.jpos, m.mcnt, m.dpos, m.ydiag, m.pos, m.nonslack, m.jrowptr, m.jcol)
    CUDA.unsafe_free!(a)
  end
  return nothing
end

# --- launches ---------------------------------------------------------------

function _launch_mismatch!(st::DeviceNewtonState)
  m = st.map
  n = m.host.n
  nb = _blocks(n)
  @cuda threads = _THREADS blocks = (nb, st.N) _k_current_mismatch!(st.I, st.F, st.rhs, st.partial, m.yrowptr, m.ycol, m.yval, st.V, st.S, st.Vset, st.btype, m.pos, st.active)
  @cuda threads = _THREADS blocks = st.N _k_reduce_max!(st.stat, st.partial, Int32(nb), st.active)
  return nothing
end

function _launch_recompute!(st::DeviceNewtonState)
  m = st.map
  nb = _blocks(m.host.n)
  @cuda threads = _THREADS blocks = (nb, st.N) _k_mismatch_from_current!(st.F, st.rhs, st.partial, st.I, st.V, st.S, st.Vset, st.btype, m.pos, st.recompute)
  @cuda threads = _THREADS blocks = st.N _k_reduce_max!(st.stat, st.partial, Int32(nb), st.recompute)
  return nothing
end

function _launch_jacobian!(st::DeviceNewtonState, vm_eps::Float64)
  m = st.map
  h = m.host
  nnzJ = Int64(h.nnzJ)
  @cuda threads = _THREADS blocks = (_blocks(h.nnzY), st.N) _k_jacobian_offdiag!(st.nz, m.yrow, m.jpos, m.mcnt, m.yval, st.V, st.btype, st.active, nnzJ)
  @cuda threads = _THREADS blocks = (_blocks(h.n - 1), st.N) _k_jacobian_diag!(st.nz, st.I, st.V, m.yval, m.ydiag, m.dpos, m.mcnt, st.dP, st.dQ, st.btype, m.nonslack, st.active, vm_eps, nnzJ)
  return nothing
end

function _launch_update!(st::DeviceNewtonState, polar::Bool)
  m = st.map
  @cuda threads = _THREADS blocks = (_blocks(m.host.n - 1), st.N) _k_update!(st.V, st.sol, st.alpha, m.nonslack, st.active, polar)
  return nothing
end

# norm and Q-limit values into st.sample (after _launch_mismatch!)
function _launch_sample!(st::DeviceNewtonState, nreq::Int)
  @cuda threads = _THREADS blocks = _blocks(nreq) _k_qlimit_sample!(st.sample, st.req, Int32(nreq), st.I, st.V, st.stat)
  return nothing
end

function _launch_scatter!(st::DeviceNewtonState, m::Int)
  @cuda threads = _THREADS blocks = _blocks(m) _k_switch_scatter!(st.S, st.btype, st.swbuf, Int32(m))
  return nothing
end

# --- event timing -------------------------------------------------------------

# CUDA event pairs per phase, recorded on the task stream (the stream the
# cuDSS handle was created on) and read once after the loop, so the timing
# adds no host synchronization inside an iteration.
struct _EventTimer
  pairs::Dict{Symbol,Vector{Tuple{CuEvent,CuEvent}}}
end
_EventTimer() = _EventTimer(Dict{Symbol,Vector{Tuple{CuEvent,CuEvent}}}())

function _timed_phase(f, et::_EventTimer, phase::Symbol)
  e1 = CuEvent()
  CUDA.record(e1)
  f()
  e2 = CuEvent()
  CUDA.record(e2)
  push!(get!(Vector{Tuple{CuEvent,CuEvent}}, et.pairs, phase), (e1, e2))
  return nothing
end

function _phase_seconds(et::_EventTimer, phase::Symbol)
  haskey(et.pairs, phase) || return 0.0
  return sum(Float64(CUDA.elapsed(e1, e2)) for (e1, e2) in et.pairs[phase]; init = 0.0)
end

# --- proof entry --------------------------------------------------------------

"""
    _device_jacobian_mismatch(Ybus, slack, Vmat, btypes, Smat, Vsetmat, dPmat, dQmat; T = Float64, vm_eps = 1e-9)

For the proof script: the mismatch and the Jacobian of N scenarios (column s
of each matrix, `btypes` a vector of N bus type vectors) built by the same
map and kernels as the device Newton loop, as host copies:
`(map, jrowptr, jcol, J, F, I, vm, norm)` with `J` of size nnzJ x N
(element type T), `F` of size nvar x N, `I` and `vm` (abs(V)) n x N and
`norm` the device mismatch norm per scenario.
"""
function _device_jacobian_mismatch(Ybus::SparseMatrixCSC{ComplexF64}, slack::Int, Vmat::Matrix{ComplexF64}, btypes::Vector{Vector{Symbol}}, Smat::Matrix{ComplexF64}, Vsetmat::Matrix{Float64}, dPmat::Matrix{Float64}, dQmat::Matrix{Float64}; T::Type = Float64, vm_eps::Float64 = 1e-9)
  hmap = HostNewtonMap(Ybus, slack)
  dmap = DeviceNewtonMap(hmap)
  N = size(Vmat, 2)
  st = DeviceNewtonState{T}(dmap, N)
  try
    copyto!(st.V, Vmat)
    copyto!(st.S, Smat)
    copyto!(st.Vset, Vsetmat)
    copyto!(st.dP, dPmat)
    copyto!(st.dQ, dQmat)
    codes = Matrix{UInt8}(undef, hmap.n, N)
    for s = 1:N
      _bus_type_codes!(view(codes, :, s), btypes[s])
    end
    copyto!(st.btype, codes)
    _launch_mismatch!(st)
    _launch_jacobian!(st, vm_eps)
    vm = CUDA.zeros(Float64, hmap.n, N)
    @cuda threads = _THREADS blocks = (_blocks(hmap.n), N) _k_vm!(vm, st.V)
    CUDA.synchronize()
    J = reshape(Array(st.nz), hmap.nnzJ, N)
    return (map = hmap, jrowptr = hmap.jrowptr, jcol = hmap.jcol, J = J, F = Array(st.F), I = Array(st.I), vm = Array(vm), norm = Array(st.stat))
  finally
    _free!(st)
  end
end

# For the proof script: the Q-limit sample kernel on V (N = 1) for the
# requests `req` (b: Q of bus b, -b: abs(V_b)), as the loop runs it after
# the mismatch kernel. Returns (norm, values).
function _device_qlimit_sample(Ybus::SparseMatrixCSC{ComplexF64}, slack::Int, V::Vector{ComplexF64}, btypes::Vector{Symbol}, S::Vector{ComplexF64}, Vset::Vector{Float64}, req::Vector{Int32})
  hmap = HostNewtonMap(Ybus, slack)
  n = hmap.n
  st = DeviceNewtonState{Float64}(DeviceNewtonMap(hmap), 1; qcap = cld(length(req), 2))
  try
    copyto!(st.V, reshape(V, n, 1))
    copyto!(st.S, reshape(S, n, 1))
    copyto!(st.Vset, reshape(Vset, n, 1))
    copyto!(st.btype, reshape(_bus_type_codes!(Vector{UInt8}(undef, n), btypes), n, 1))
    copyto!(st.req, 1, req, 1, length(req))
    _launch_mismatch!(st)
    _launch_sample!(st, length(req))
    out = Array(st.sample)
    return (norm = out[1], values = out[2:(1+length(req))])
  finally
    _free!(st)
  end
end

# The kernels of the device loop with the argument types of an actual launch
# on `Ybus`, for the static FMA check of the proof script (CUDA.code_ptx
# needs the device-side types). Returns (name, kernel, argument tuple type).
function _device_kernel_signatures(Ybus::SparseMatrixCSC{ComplexF64}, slack::Int; T::Type = Float64)
  hmap = HostNewtonMap(Ybus, slack)
  st = DeviceNewtonState{T}(DeviceNewtonMap(hmap), 1)
  m = st.map
  tt(args...) = Tuple{map(a -> typeof(CUDA.cudaconvert(a)), args)...}
  try
    return [
      (:_k_current_mismatch!, _k_current_mismatch!, tt(st.I, st.F, st.rhs, st.partial, m.yrowptr, m.ycol, m.yval, st.V, st.S, st.Vset, st.btype, m.pos, st.active)),
      (:_k_mismatch_from_current!, _k_mismatch_from_current!, tt(st.F, st.rhs, st.partial, st.I, st.V, st.S, st.Vset, st.btype, m.pos, st.recompute)),
      (:_k_reduce_max!, _k_reduce_max!, tt(st.stat, st.partial, Int32(1), st.active)),
      (:_k_jacobian_offdiag!, _k_jacobian_offdiag!, tt(st.nz, m.yrow, m.jpos, m.mcnt, m.yval, st.V, st.btype, st.active, Int64(1))),
      (:_k_jacobian_diag!, _k_jacobian_diag!, tt(st.nz, st.I, st.V, m.yval, m.ydiag, m.dpos, m.mcnt, st.dP, st.dQ, st.btype, m.nonslack, st.active, 1e-9, Int64(1))),
      (:_k_vm!, _k_vm!, tt(st.F, st.V)),
      (:_k_qlimit_sample!, _k_qlimit_sample!, tt(st.sample, st.req, Int32(1), st.I, st.V, st.stat)),
      (:_k_switch_scatter!, _k_switch_scatter!, tt(st.S, st.btype, st.swbuf, Int32(1))),
      (:_k_update!, _k_update!, tt(st.V, st.sol, st.alpha, m.nonslack, st.active, true)),
    ]
  finally
    _free!(st)
  end
end

# The update kernel on given V and solution columns (N = 1), for the proof:
# the iterate after one step, to compare with _apply_rectangular_delta!.
function _device_update(Ybus::SparseMatrixCSC{ComplexF64}, slack::Int, V::Vector{ComplexF64}, dx::Vector{Float64}, alpha::Float64, polar::Bool)
  hmap = HostNewtonMap(Ybus, slack)
  dmap = DeviceNewtonMap(hmap)
  st = DeviceNewtonState{Float64}(dmap, 1)
  try
    copyto!(st.V, reshape(V, :, 1))
    copyto!(st.sol, reshape(dx, :, 1))
    fill!(st.alpha, alpha)
    _launch_update!(st, polar)
    return vec(Array(st.V))
  finally
    _free!(st)
  end
end

# --- the solver loop ----------------------------------------------------------

# set once the first device solve of the process has run (kernel compilation)
const _DEVICE_NEWTON_RAN = Ref(false)

# Test hook of the recovery check (benchmark/gpu/check_device_newton_gpu.jl):
# a function (point, iteration, state) -> Bool, called at every point of the
# loop where GPU work can fail. It may throw a GPU error, damage the device
# state (a non-finite solution, a kernel fault that kills the context) or
# return true to take that point's failure branch (:check a failed step,
# :cpu_step a singular CPU step). `nothing` outside the test.
const _FAULT_HOOK = Ref{Any}(nothing)

function _fault(point::Symbol, it::Int, st)
  h = _FAULT_HOOK[]
  h === nothing && return false
  return h(point, it, st)::Bool
end

# What one device solve holds on the GPU, filled step by step by the
# pre-flight, so a GPU error at any point releases exactly what exists.
mutable struct _DeviceSolve
  state::Any
  solver::Any
  X::Any
  B::Any
end
_DeviceSolve() = _DeviceSolve(nothing, nothing, nothing, nothing)

# Host side of the Q-limit exchange of one solve. `qset` are the buses that
# can be PV during the solve (PV at the start or PV in the network data: the
# handler only turns PV buses to PQ and only originally PV buses back).
# `codes` and `Sdev` mirror the bus type codes and S on the device.
mutable struct _QLimitLink
  qset::Vector{Int}
  req::Vector{Int32}
  vals::Vector{Float64}
  sample::Sparlectra.QLimitSample
  codes::Vector{UInt8}
  Sdev::Vector{ComplexF64}
  swbuf::Vector{Float64}
end

function _QLimitLink(input::Sparlectra.DeviceNewtonInput, codes::Vector{UInt8})
  n = length(input.V)
  mask = input.qlimit_args.pv_orig_mask
  qset = [b for b = 1:n if input.bus_types[b] === :PV || (b <= length(mask) && mask[b])]
  return _QLimitLink(qset, Int32[], zeros(2 + 2length(qset)), Sparlectra.QLimitSample(fill(NaN, n), fill(NaN, n)), codes, copy(input.S), Float64[])
end

# The requests of the next iterations: Q of every PV bus (the PV->PQ check
# and the :auto start rule read it), and for every release candidate (PQ
# now, PV in the network data, a clamp in qLimitEvents) its abs(V) for the
# voltage rule and its Q for the band rule and the release record. Without
# a cooldown (cooldown_iters <= 0) the handler can release a bus in the
# same call that clamped it, so every PV bus that is PV in the network data
# also sends abs(V). Device loop runs have qlimit_mode :switch_to_pq, so
# allow_reenable applies as given.
function _qlimit_requests!(link::_QLimitLink, input::Sparlectra.DeviceNewtonInput)
  empty!(link.req)
  qa = input.qlimit_args
  same_call_release = qa.allow_reenable && qa.cooldown_iters <= 0
  for b in link.qset
    orig_pv = b <= length(qa.pv_orig_mask) && qa.pv_orig_mask[b]
    if input.bus_types[b] === :PV
      push!(link.req, Int32(b))
      same_call_release && orig_pv && push!(link.req, Int32(-b))
    elseif qa.allow_reenable && orig_pv && haskey(input.net.qLimitEvents, b)
      push!(link.req, Int32(b))
      push!(link.req, Int32(-b))
    end
  end
  return length(link.req)
end

# the downloaded values into the sample the handler reads; every other
# entry is NaN, so a value from an earlier iteration is never read
function _qlimit_unpack!(link::_QLimitLink)
  q = link.sample.q
  vm = link.sample.vm
  fill!(q, NaN)
  fill!(vm, NaN)
  @inbounds for (k, r) in enumerate(link.req)
    r > 0 ? (q[r] = link.vals[1+k]) : (vm[-r] = link.vals[1+k])
  end
  return nothing
end

# The buses of `qset` whose type code or S on the host differ from the
# device after the Q-limit handler, packed for _k_switch_scatter!; the
# mirrors are updated. Returns the number of buses.
function _switch_pack!(link::_QLimitLink, input::Sparlectra.DeviceNewtonInput)
  empty!(link.swbuf)
  m = 0
  for b in link.qset
    c = _bus_type_code(input.bus_types[b])
    s = input.S[b]
    (c == link.codes[b] && s === link.Sdev[b]) && continue
    link.codes[b] = c
    link.Sdev[b] = s
    push!(link.swbuf, Float64(b), Float64(c), real(s), imag(s))
    m += 1
  end
  return m
end

# The solver state the device loop changes in place (the Q-limit handler and
# the injection rebuild), taken before the first iteration. After a context
# loss nothing can be read from the device; the state is put back and the
# CPU loop runs the whole solve, as after a pre-flight failure.
struct _StartState
  S::Vector{ComplexF64}
  bus_types::Vector{Symbol}
  dP::Vector{Float64}
  dQ::Vector{Float64}
  events::Dict{Int,Symbol}
  nlog::Int
  qset::Vector{Int}
  qgen::Vector{Union{Nothing,Float64}}
  prev_q::Vector{Float64}
end

function _StartState(input::Sparlectra.DeviceNewtonInput, qset::Vector{Int})
  net = input.net
  return _StartState(copy(input.S), copy(input.bus_types), copy(input.dPinj_dVm), copy(input.dQinj_dVm), copy(net.qLimitEvents), length(net.qLimitLog), qset, [net.nodeVec[b]._qƩGen for b in qset], copy(input.workspace.prev_pv_qreq_pu))
end

function _restore!(input::Sparlectra.DeviceNewtonInput, s::_StartState)
  net = input.net
  copyto!(input.S, s.S)
  copyto!(input.bus_types, s.bus_types)
  copyto!(input.dPinj_dVm, s.dP)
  copyto!(input.dQinj_dVm, s.dQ)
  empty!(net.qLimitEvents)
  merge!(net.qLimitEvents, s.events)
  resize!(net.qLimitLog, s.nlog)
  for (k, b) in enumerate(s.qset)
    net.nodeVec[b]._qƩGen = s.qgen[k]
  end
  copyto!(input.workspace.prev_pv_qreq_pu, s.prev_q)
  return nothing
end

function _device_preflight!(ds::_DeviceSolve, input::Sparlectra.DeviceNewtonInput, ctx::CudssNewtonContext{T}, tm::Sparlectra.DeviceNewtonTimings, h2d!::F) where {T,F}
  n = length(input.V)
  _fault(:preflight, 0, nothing)
  t0 = time_ns()
  hmap = HostNewtonMap(input.Ybus, input.slack_idx)
  tm.map_s = (time_ns() - t0) / 1.0e9
  t0 = time_ns()
  codes = _bus_type_codes!(Vector{UInt8}(undef, n), input.bus_types)
  link = _QLimitLink(input, codes)
  dmap = DeviceNewtonMap(hmap)
  st = DeviceNewtonState{T}(dmap, 1; qcap = input.qlimits_enabled ? length(link.qset) : 0)
  ds.state = st
  h2d!(st.V, input.V, n)
  h2d!(st.S, input.S, n)
  h2d!(st.Vset, input.Vset, n)
  h2d!(st.dP, input.dPinj_dVm, n)
  h2d!(st.dQ, input.dQinj_dVm, n)
  h2d!(st.btype, codes, n)
  fill!(st.alpha, input.damp)
  nreq = input.qlimits_enabled ? _qlimit_requests!(link, input) : 0
  nreq > 0 && h2d!(st.req, link.req, nreq)
  # J(V) at the start iterate: cuDSS analyses with values (matching and
  # scaling), and the factorization serves iteration 1 when nothing changes
  # before its step
  _launch_mismatch!(st)
  _launch_jacobian!(st, 1e-9)
  _fault(:preflight_kernels, 0, st)
  # a fault in the uploads or kernels surfaces here, before a cuDSS object
  # exists whose finalizer would crash on a dead context
  CUDA.synchronize()
  tm.upload_s = (time_ns() - t0) / 1.0e9
  solver = _cudss_configure!(CudssSolver(dmap.jrowptr, dmap.jcol, st.nz, "G", 'F'), ctx.fast_mode)
  ds.solver = solver
  # no pivot perturbation: a zero pivot gives a non-finite solution, which
  # the loop sees as a non-finite mismatch norm (UMFPACK fails on the same
  # matrices), instead of a finite wrong step only "npivots" would show
  cudss_set(solver, "pivot_epsilon", zero(T))
  ds.X = CudssMatrix(vec(st.sol))
  ds.B = CudssMatrix(vec(st.rhs))
  t0 = time_ns()
  cudss("analysis", solver, ds.X, ds.B)
  CUDA.synchronize()
  tm.analysis_s = (time_ns() - t0) / 1.0e9
  t0 = time_ns()
  cudss("factorization", solver, ds.X, ds.B)
  CUDA.synchronize()
  tm.factor_first_s = (time_ns() - t0) / 1.0e9
  ctx.analyze_count += 1
  return (link = link, nreq = nreq)
end

# The device loop owns its history vectors; the caller appends them.
mutable struct _LoopLocals
  iters::Int
  converged::Bool
  rejection_reason::Symbol
  history::Vector{Float64}
  pre_switch_history::Vector{Float64}
  type_switched::Vector{Bool}
  step_diagnostics::Vector{NamedTuple}
  best_mis::Float64
  best_it::Int
  last_it::Int
  first_nf_it::Int
  first_nf_v::Union{Nothing,Vector{ComplexF64}}
  # host copies of the best and the last finite iterate, filled at the end
  vbest::Union{Nothing,Vector{ComplexF64}}
  vlast::Union{Nothing,Vector{ComplexF64}}
  changes::Int
  reenables::Int
  switch_its::Int
  steps::Int
  # the iteration whose history entry still waits for the norm after its
  # Q-limit switch (0 for none); that norm comes with the next download
  switch_norm_due::Int
end

_LoopLocals() = _LoopLocals(0, false, :nr_mismatch_not_converged, Float64[], Float64[], Bool[], NamedTuple[], Inf, 0, 0, 0, nothing, nothing, nothing, 0, 0, 0, 0, 0)

# the iterate measured in iteration `it` is the last finite one, and the
# best one if its norm is below the best so far
function _record_measured!(L::_LoopLocals, it::Int, max_mis::Float64)
  L.last_it = it
  if max_mis < L.best_mis
    L.best_mis = max_mis
    L.best_it = it
  end
  return nothing
end

# The norm after the switch of iteration L.switch_norm_due arrived with a
# download: into its history entry and, if that iteration's step is the
# pending one, into the step's profile row. Returns the row.
function _switch_norm_arrived!(L::_LoopLocals, v::Float64, pending::Int, row::NamedTuple)
  due = L.switch_norm_due
  L.history[due] = v
  L.switch_norm_due = 0
  return due == pending ? merge(row, (max_mismatch = v,)) : row
end

# A recovery path for the step of iteration `it` whose norm after the switch
# has not arrived: computed on the CPU from the iterate measured in `it`
# (`Vm`, read from Vlast) with the types and S after the switch. The device
# mismatch is bitwise the CPU one (check_jacobian_gpu.jl), so the entry is
# the value the download would have brought. Returns the row.
function _settle_switch_norm!(L::_LoopLocals, input::Sparlectra.DeviceNewtonInput, Vm::Vector{ComplexF64}, it::Int, row::NamedTuple)
  L.switch_norm_due == it || return row
  F = Sparlectra.mismatch_rectangular(input.Ybus, Vm, input.S, input.bus_types, input.Vset, input.slack_idx)
  return _switch_norm_arrived!(L, maximum(abs.(F)), it, row)
end

function _gpu_newton(input::Sparlectra.DeviceNewtonInput)
  return _gpu_newton(input, input.linear_ctx)
end

function _gpu_newton(input::Sparlectra.DeviceNewtonInput, ctx::CudssNewtonContext{T}) where {T}
  t_core = time_ns()
  tm = Sparlectra.DeviceNewtonTimings()
  tm.first_call = !_DEVICE_NEWTON_RAN[]
  ctx.device_timings = tm
  ds = _DeviceSolve()
  d2h!(dst, src, len) = (t0 = time_ns(); copyto!(dst, 1, src, 1, len); tm.d2h_s += (time_ns() - t0) / 1.0e9; tm.d2h_bytes += len * sizeof(eltype(dst)); tm.d2h_count += 1; dst)
  h2d!(dst, src, len) = (t0 = time_ns(); copyto!(dst, 1, src, 1, len); tm.h2d_s += (time_ns() - t0) / 1.0e9; tm.h2d_bytes += len * sizeof(eltype(src)); tm.h2d_count += 1; dst)
  pre = try
    _device_preflight!(ds, input, ctx, tm, h2d!)
  catch e
    _release!(ds)
    (_is_gpu_error(e) && ctx.fallback_to_cpu) || rethrow(e)
    _switch_to_cpu!(ctx, string("GPU error before the device Newton loop, the solve ran on the CPU: ", first(sprint(showerror, e), 300)))
    return nothing
  end
  # the pre-flight uploads are counted in upload_s, not per iteration
  tm.h2d_s = 0.0
  tm.h2d_bytes = 0
  tm.h2d_count = 0
  _DEVICE_NEWTON_RAN[] = true
  ctx.device_newton = true
  try
    return _device_loop!(input, ctx, ds, pre.link, pre.nreq, tm, t_core, d2h!, h2d!)
  finally
    _release!(ds)
  end
end

# A kernel fault (illegal address) leaves the context dead: every later call
# returns the same error. CUDA.jl then keeps the arrays alive with a logged
# error, but the finalizer of the cuDSS data segfaults in cudssDataDestroy
# and takes the process down at the next GC (measured with CUDSS.jl 1.0.0,
# benchmark/gpu/check_device_newton_gpu.jl part 4). On a dead context the
# data handle is cleared, so the finalizer gets an invalid-value error
# instead, and nothing is freed.
function _release!(ds::_DeviceSolve)
  if _context_dead()
    ds.solver === nothing || (ds.solver.data.data = C_NULL)
    return nothing
  end
  ds.state === nothing || _free!(ds.state)
  return nothing
end

# true when the context carries a sticky error; synchronize is the probe
function _context_dead()
  try
    CUDA.synchronize()
  catch e
    e isa CUDA.CuError || rethrow(e)
    return true
  end
  return false
end

function _switch_to_cpu!(ctx::CudssNewtonContext, reason::String)
  ctx.fallback_reason = reason
  @warn reason
  ctx.fallback_count += 1
  ctx.cpu = Sparlectra.UmfpackReuseNewtonContext()
  return nothing
end

# The step of iteration `step_it` on the CPU from `Vstart`, after the GPU
# step failed or a GPU error hit it; the rest of the solve stays on the CPU.
# `row` is the performance-profile row of that iteration. A singular CPU
# step ends the solve as in the CPU loop (:singular_newton_step, iterations
# counted up to the failed step, no profile row). Returns the iterate the
# solve continues from and the first iteration left to the CPU loop.
function _redo_step_on_cpu!(input::Sparlectra.DeviceNewtonInput, ctx::CudssNewtonContext, L::_LoopLocals, st::DeviceNewtonState, step_it::Int, Vstart::Vector{ComplexF64}, changed::Bool, row::NamedTuple, why::String)
  _switch_to_cpu!(ctx, string(why, " in iteration ", step_it, ", the rest of the solve ran on the CPU"))
  Vnew = _fault(:cpu_step, step_it, st) ? nothing : Sparlectra._device_newton_cpu_step(input, Vstart, L.step_diagnostics, changed)
  if Vnew === nothing
    L.converged = false
    L.rejection_reason = :singular_newton_step
    L.iters = step_it
    return Vstart, input.maxiter + 1
  end
  L.steps += 1
  Sparlectra._perf_profile_push_iteration!(input.performance_profile, row)
  return Vnew, step_it + 1
end

# Host copies of the best and the last finite iterate for the outcome, read
# from the device at the end of the solve. `Vh` is a host copy of the device
# iterate and `vh_it` the iteration whose mismatch measured it (0 if none):
# a record of that iteration is copied from it instead of downloaded.
function _host_records!(L::_LoopLocals, st::DeviceNewtonState, Vh::Vector{ComplexF64}, vh_it::Int, d2h!::F) where {F}
  _fault(:host_records, L.iters, st)
  n = length(Vh)
  if L.last_it > 0
    L.vlast = L.last_it == vh_it ? copy(Vh) : d2h!(Vector{ComplexF64}(undef, n), st.Vlast, n)
  end
  if L.best_it > 0
    L.vbest = L.best_it == L.last_it ? copy(L.vlast) : d2h!(Vector{ComplexF64}(undef, n), st.Vbest, n)
  end
  return nothing
end

# One device Newton solve. V stays on the device; the device keeps copies of
# the iterate of the last finite mismatch (Vlast, also the iterate before a
# step not checked yet) and of the best one (Vbest). Per iteration the host
# reads one download: the mismatch norm and, with Q limits on, the Q-limit
# values and the norm after the switch of the previous iteration, if it had
# one. V goes to the host each iteration only when voltage-dependent
# injections are rebuilt from it on the CPU or the Q-limit trace prints it
# (`host_v`); then the norm is a second download, because the mismatch
# needs the injections the CPU builds from that V. The norm could also come
# from the CPU mismatch on the host copy of V (bitwise the same value), but
# measured on case6495rte that was not faster: the CPU mismatch takes longer
# than the download, which only waits for the injection upload and the
# mismatch kernel.
#
# A GPU step is checked when the next iteration reads its mismatch norm: a
# non-finite norm means a non-finite step (cuDSS runs without pivot
# perturbation, so a singular J gives one). Until then it is `pending`, and
# a failed check or a GPU error in that window redoes it on the CPU from
# Vlast. A GPU error with no pending step lands in one of three places: before
# the history entry of iteration `it` is pushed (stage :measure, and
# :records while the device copies of the measured iterate run) the CPU loop
# takes iteration `it` from V; after the Q-limit handler (:switch, :step)
# only the step of `it` is left and is done on the CPU from Vlast; after the
# loop (:final) only V and the records are read again. If the error killed the
# context, nothing can be read from the device: the state of the solve goes
# back to its start (_StartState) and the CPU loop runs the whole solve.
function _device_loop!(input::Sparlectra.DeviceNewtonInput, ctx::CudssNewtonContext{T}, ds::_DeviceSolve, link::_QLimitLink, nreq::Int, tm::Sparlectra.DeviceNewtonTimings, t_core::UInt64, d2h!::F1, h2d!::F2) where {T,F1,F2}
  st = ds.state
  solver = ds.solver
  X = ds.X
  B = ds.B
  net = input.net
  n = length(input.V)
  et = _EventTimer()
  L = _LoopLocals()
  qlim = input.qlimits_enabled
  host_v = input.has_vdep_control || (qlim && input.qlimit_args.qlimit_trace_enabled)
  normbuf = zeros(Float64, 1)
  # host copy of V, kept current only in host_v mode
  Vcur = copy(input.V)
  # the start of the solve, for a context loss
  Vstart0 = copy(input.V)
  start = _StartState(input, link.qset)
  polar = input.newton_update === :polar
  maxiter = input.maxiter
  # the factorization of the pre-flight is J of iteration 1 as long as
  # nothing changes S, the types or the injection derivatives before the step
  fact_current = !input.has_vdep_control
  # the GPU step not checked yet: its iteration (0 for none), whether the
  # active set changed in it, and its performance-profile row
  pending = 0
  pending_changed = false
  pending_row = (iteration = 0, max_mismatch = NaN, qlimit_changed = false, qlimit_reenabled = false)
  changed = false
  reenabled = false
  it = 0
  # the iteration whose mismatch measured the current device V (0 after a step)
  measured = 0
  max_mis = NaN
  stage = :measure
  step_failed = string("the GPU step returned a non-finite iterate (singular or ill-conditioned Jacobian in ", ctx.precision, ")")

  try
    for it_ = 1:maxiter
      it = it_
      stage = :measure
      changed = false
      reenabled = false
      input.set_phase("newton_iteration")
      L.iters = it
      if host_v
        if pending > 0
          # V goes to the host for the injections or the trace, which
          # also checks the step
          d2h!(Vcur, st.V, n)
          if _fault(:check, it, st) || !all(isfinite, Vcur)
            return _failed_step!(input, ctx, L, st, pending, pending_changed, pending_row, step_failed, tm, et, t_core, d2h!)
          end
          _confirm_step!(L, input, pending_row)
          pending = 0
        end
        if input.has_vdep_control
          t0 = time_ns()
          Snew, dPn, dQn = Sparlectra.buildControlledSVec(net, Vcur)
          copyto!(input.S, Snew)
          copyto!(input.dPinj_dVm, dPn)
          copyto!(input.dQinj_dVm, dQn)
          Sparlectra._reapply_qlimit_clamps!(input.S, input.dQinj_dVm, net, input.bus_types, Vcur)
          tm.host_injections_s += (time_ns() - t0) / 1.0e9
          _fault(:vdep, it, st)
          h2d!(st.S, input.S, n)
          h2d!(st.dP, input.dPinj_dVm, n)
          h2d!(st.dQ, input.dQinj_dVm, n)
          copyto!(link.Sdev, input.S)
        end
      end
      _fault(:mismatch, it, st)
      _timed_phase(() -> _launch_mismatch!(st), et, :mismatch)
      if qlim && !host_v
        _launch_sample!(st, nreq)
        # the norm, the Q-limit values and, after a switch in the previous
        # iteration, the norm after it: one download, one synchronization
        due = L.switch_norm_due
        d2h!(link.vals, st.sample, 1 + nreq + (due > 0 ? 1 : 0))
        max_mis = link.vals[1]
        due > 0 && (pending_row = _switch_norm_arrived!(L, link.vals[2+nreq], pending, pending_row))
      else
        d2h!(normbuf, st.stat, 1)
        max_mis = normbuf[1]
      end
      measured = it
      if pending > 0
        # device mode: the norm of the stepped iterate checks the step
        if _fault(:check, it, st) || !isfinite(max_mis)
          return _failed_step!(input, ctx, L, st, pending, pending_changed, pending_row, step_failed, tm, et, t_core, d2h!)
        end
        _confirm_step!(L, input, pending_row)
        pending = 0
      end
      Sparlectra.sparlectra_check_abort()
      if isfinite(max_mis)
        # the device copies of the measured iterate come before the history
        # entry: a GPU error in them leaves the whole iteration to the CPU
        # loop, which pushes the entry itself
        stage = :records
        _fault(:records, it, st)
        copyto!(st.Vlast, st.V)
        max_mis < L.best_mis && copyto!(st.Vbest, st.V)
        _record_measured!(L, it, max_mis)
        stage = :measure
      end
      push!(L.history, max_mis)
      push!(L.pre_switch_history, max_mis)
      push!(L.type_switched, false)
      if !isfinite(max_mis)
        L.first_nf_it == 0 && (L.first_nf_it = it)
        L.converged = false
        L.rejection_reason = :nr_nonfinite
        break
      end

      # --- Q limits on the CPU, the unchanged handler -------------------
      # no GPU work between the download above and the end of the handler:
      # a GPU error is either before (stage :measure) or after it
      input.set_phase("q_limit_processing")
      converged_this_iter = max_mis <= input.tol
      if qlim
        host_v || _qlimit_unpack!(link)
        t0 = time_ns()
        qr = Sparlectra._handle_rectangular_qlimit_iteration!(net, it, n, input.Ybus, host_v ? Vcur : link.sample, input.S, input.bus_types, input.Vset, input.slack_idx, input.workspace, L.history, input.qlimit_args...)
        tm.host_qlimit_s += (time_ns() - t0) / 1.0e9
        changed = qr.changed
        reenabled = qr.reenabled
        converged_this_iter = qr.converged_this_iter
      end
      stage = :switch
      # counted before the device follows the switch, so a GPU error in the
      # uploads below leaves the counters and type_switched as the CPU
      # loop has them
      changed && (L.changes += 1)
      reenabled && (L.reenables += 1)
      (changed || reenabled) && (L.switch_its += 1; L.type_switched[end] = true)
      if changed || reenabled
        # with V on the device the history entry of this iteration waits for
        # the norm after the switch; with V on the host the handler set it
        host_v || (L.switch_norm_due = it)
        _fault(:switch, it, st)
        m = _switch_pack!(link, input)
        if m > 0
          h2d!(st.swbuf, link.swbuf, 4m)
          _launch_scatter!(st, m)
        end
        if !host_v
          nreq = _qlimit_requests!(link, input)
          nreq > 0 && h2d!(st.req, link.req, nreq)
        end
        fill!(st.recompute, 0x01)
        _timed_phase(() -> _launch_recompute!(st), et, :mismatch)
        # the norm after the switch, which the CPU loop records in history,
        # goes behind the Q-limit values of the next download (no host
        # synchronization here); with V on the host the handler computed it
        # on the CPU (bitwise the same value, check_jacobian_gpu.jl)
        host_v || copyto!(st.sample, 2 + nreq, st.stat, 1, 1)
      end
      input.check_cancel()
      if converged_this_iter && !(changed || reenabled)
        L.converged = true
        L.rejection_reason = :none
        break
      end

      # --- Newton step on the device -------------------------------------
      stage = :step
      input.set_phase("linear_solve")
      Sparlectra._validate_rectangular_damping(input.damp, min(input.autodamp_min, input.damp))
      _fault(:step, it, st)
      if !(fact_current && it == 1 && !(changed || reenabled))
        _timed_phase(() -> _launch_jacobian!(st, 1e-9), et, :build)
        _timed_phase(() -> cudss("refactorization", solver, X, B), et, :refactor)
        ctx.refactor_count += 1
      end
      fact_current = false
      _timed_phase(() -> cudss("solve", solver, X, B), et, :solve)
      _fault(:update, it, st)
      _timed_phase(() -> _launch_update!(st, polar), et, :update)
      measured = 0
      # max_mismatch is the value after the Q-limit handling, as in the CPU loop
      pending = it
      pending_changed = changed || reenabled
      pending_row = (iteration = it, max_mismatch = L.history[end], qlimit_changed = changed, qlimit_reenabled = reenabled)
      input.check_cancel()
    end
    # the loop is over (converged, a non-finite norm, or max_iter used up
    # with the step of maxiter pending) and every history entry is pushed
    stage = :final
    _fault(:final, it, st)
    # max_iter used up right after a switch: its norm was not read yet
    if L.switch_norm_due > 0
      d2h!(link.vals, st.sample, 2 + nreq)
      pending_row = _switch_norm_arrived!(L, link.vals[2+nreq], pending, pending_row)
    end
    # the end of the solve: V to the host; a step the loop ended with (max
    # iterations used up) is checked on it
    Vh = d2h!(Vector{ComplexF64}(undef, n), st.V, n)
    if pending > 0
      if _fault(:check, maxiter + 1, st) || !all(isfinite, Vh)
        return _failed_step!(input, ctx, L, st, pending, pending_changed, pending_row, step_failed, tm, et, t_core, d2h!)
      end
      _confirm_step!(L, input, pending_row)
      pending = 0
    end
    L.first_nf_it > 0 && L.first_nf_it == measured && (L.first_nf_v = copy(Vh))
    _host_records!(L, st, Vh, measured, d2h!)
    return _outcome(input, ctx, L, Vh, maxiter + 1, tm, et, t_core, true)
  catch e
    (_is_gpu_error(e) && ctx.fallback_to_cpu) || rethrow(e)
    why = string("GPU error (", first(sprint(showerror, e), 300), ")")
    if _context_dead()
      # nothing on the device can be read: the CPU loop runs the solve from
      # its start. Performance-profile rows of the device iterations stay.
      _switch_to_cpu!(ctx, string(why, " in iteration ", it, ", the CUDA context is lost, the solve ran again on the CPU from its start"))
      _restore!(input, start)
      return _outcome(input, ctx, _LoopLocals(), Vstart0, 1, tm, et, t_core, false)
    end
    if pending > 0
      # the step of iteration `pending` is not checked: it is redone
      return _failed_step!(input, ctx, L, st, pending, pending_changed, pending_row, why, tm, et, t_core, d2h!)
    end
    if stage === :final
      # the loop ended and no step is pending (converged, a non-finite norm
      # in iteration 1 or in host_v mode, or the last step already checked):
      # only V and the records are left to read. They are read again; the
      # CPU loop gets no iteration, so no history entry and no step is added.
      _switch_to_cpu!(ctx, string(why, " in iteration ", it, " after the device loop ended, V read again"))
      Vh = d2h!(Vector{ComplexF64}(undef, n), st.V, n)
      L.first_nf_it > 0 && L.first_nf_it == measured && (L.first_nf_v = copy(Vh))
      _host_records!(L, st, Vh, measured, d2h!)
      return _outcome(input, ctx, L, Vh, maxiter + 1, tm, et, t_core, false)
    end
    if stage === :measure || stage === :records
      # no history entry of this iteration yet: the CPU loop takes the whole
      # iteration from V (unchanged since its step). In :records the device
      # copies may be half done; the host records then name V itself as the
      # last (and maybe best) iterate, which is what the copies would say.
      stage === :records && _record_measured!(L, it, max_mis)
      _switch_to_cpu!(ctx, string(why, " in iteration ", it, ", the rest of the solve ran on the CPU"))
      Vh = d2h!(Vector{ComplexF64}(undef, n), st.V, n)
      _host_records!(L, st, Vh, measured, d2h!)
      return _outcome(input, ctx, L, Vh, it, tm, et, t_core, false)
    end
    # the handler ran: only the step of this iteration is left, from the
    # measured iterate in Vlast (the update kernel may have written V)
    row = (iteration = it, max_mismatch = L.history[end], qlimit_changed = changed, qlimit_reenabled = reenabled)
    Vm = d2h!(Vector{ComplexF64}(undef, n), st.Vlast, n)
    row = _settle_switch_norm!(L, input, Vm, it, row)
    _host_records!(L, st, Vm, it, d2h!)
    Vout, next_it = _redo_step_on_cpu!(input, ctx, L, st, it, Vm, changed || reenabled, row, why)
    return _outcome(input, ctx, L, Vout, next_it, tm, et, t_core, false)
  end
end

# The pending GPU step of iteration `step_it` failed its check or a GPU error
# hit it: redone on the CPU from the iterate before it (Vlast, measured in
# `step_it`).
function _failed_step!(input::Sparlectra.DeviceNewtonInput, ctx::CudssNewtonContext, L::_LoopLocals, st::DeviceNewtonState, step_it::Int, step_changed::Bool, row::NamedTuple, why::String, tm::Sparlectra.DeviceNewtonTimings, et::_EventTimer, t_core::UInt64, d2h!::F) where {F}
  ctx.fallback_to_cpu || error(why)
  n = length(input.V)
  Vm = d2h!(Vector{ComplexF64}(undef, n), st.Vlast, n)
  row = _settle_switch_norm!(L, input, Vm, step_it, row)
  _host_records!(L, st, Vm, step_it, d2h!)
  Vout, next_it = _redo_step_on_cpu!(input, ctx, L, st, step_it, Vm, step_changed, row, why)
  return _outcome(input, ctx, L, Vout, next_it, tm, et, t_core, false)
end

# a GPU step passed its check: its rows as the CPU loop pushes them after a
# successful step (a fixed step adds no autodamp row, as on the CPU)
function _confirm_step!(L::_LoopLocals, input::Sparlectra.DeviceNewtonInput, row::NamedTuple)
  L.steps += 1
  Sparlectra._perf_profile_push_iteration!(input.performance_profile, row)
  return nothing
end

# `device_ok` is false after a GPU error or a failed step: the event times
# are then not read (the context may be dead), the device phase times stay 0.
function _outcome(input::Sparlectra.DeviceNewtonInput, ctx::CudssNewtonContext, L::_LoopLocals, V::Vector{ComplexF64}, next_it::Int, tm::Sparlectra.DeviceNewtonTimings, et::_EventTimer, t_core::UInt64, device_ok::Bool)
  if device_ok
    _fault(:outcome, L.iters, nothing)
    CUDA.synchronize()
    tm.mismatch_s = _phase_seconds(et, :mismatch)
    tm.build_s = _phase_seconds(et, :build)
    tm.refactor_s = _phase_seconds(et, :refactor)
    tm.solve_s = _phase_seconds(et, :solve)
    tm.update_s = _phase_seconds(et, :update)
  end
  # own kernels per iteration: mismatch, reduction, two Jacobian kernels,
  # update (plus cuDSS refactorization and solve; plus the Q-limit sample
  # with Q limits on)
  tm.launches_per_iteration = 5 + (input.qlimits_enabled ? 1 : 0)
  tm.core_s = (time_ns() - t_core) / 1.0e9
  pp = input.performance_profile
  for (key, val) in ((:gpu_newton_map, tm.map_s), (:gpu_newton_upload, tm.upload_s), (:gpu_newton_analysis, tm.analysis_s), (:gpu_newton_factor_first, tm.factor_first_s), (:gpu_newton_mismatch, tm.mismatch_s), (:gpu_newton_build, tm.build_s), (:gpu_newton_refactor, tm.refactor_s), (:gpu_newton_solve, tm.solve_s), (:gpu_newton_update, tm.update_s), (:gpu_newton_d2h, tm.d2h_s), (:gpu_newton_h2d, tm.h2d_s), (:gpu_newton_host_qlimit, tm.host_qlimit_s), (:gpu_newton_host_injections, tm.host_injections_s), (:gpu_newton_core, tm.core_s))
    Sparlectra._perf_profile_add!(pp, key, val)
  end
  # V, vbest and vlast are separate host buffers of this solve
  return Sparlectra.DeviceNewtonOutcome(
    V, next_it, L.iters, L.converged, L.rejection_reason,
    L.history, L.pre_switch_history, L.type_switched, L.step_diagnostics,
    L.best_mis, L.best_it, L.best_it > 0 ? L.vbest : nothing, L.last_it, L.last_it > 0 ? L.vlast : nothing, L.first_nf_it, L.first_nf_v,
    L.changes, L.reenables, L.switch_its, L.steps,
  )
end

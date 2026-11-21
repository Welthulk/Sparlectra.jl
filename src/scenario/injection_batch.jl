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

# file: src/scenario/injection_batch.jl
# purpose: the injection batch: N power flows of one network that differ
#          only in the load P and Q per bus (random scaling, a ramp as a
#          stand-in for a day profile). Q limits are off in batch mode; each
#          scenario reports how many PV buses would violate a limit. The CPU
#          path runs every scenario through runpf! on the per-chunk workers
#          of the scenario engine (_engine_parallel_map, worker reset by
#          _reset_scenario_worker!); the GPU path (CUDA extension) builds the
#          Ybus and the Jacobian pattern once and runs all scenarios on the
#          device with one cuDSS analysis.

"""
    InjectionScenarios

N scenarios of one network that differ in the load injections:
`load_scale[b, s]` multiplies P and Q of every load at bus `b` (internal
bus index) in scenario `s`. Generators keep their setpoints, the slack
takes the difference. `kind` names how the factors were made
(`:random`, `:ramp` or `:custom`) and `seed` the random seed (`nothing`
for a ramp or custom factors), so a study can record it.

See [`random_load_scenarios`](@ref), [`load_ramp_scenarios`](@ref) and
[`runInjectionBatch`](@ref).
"""
struct InjectionScenarios
  load_scale::Matrix{Float64}
  kind::Symbol
  seed::Union{Nothing,Int}
end

InjectionScenarios(load_scale::AbstractMatrix{<:Real}; kind::Symbol = :custom, seed::Union{Nothing,Int} = nothing) = InjectionScenarios(Matrix{Float64}(load_scale), kind, seed)

Base.length(s::InjectionScenarios) = size(s.load_scale, 2)

"""
    random_load_scenarios(net, N; band = 0.2, seed = 1) -> InjectionScenarios

`N` scenarios with an independent load factor per bus and scenario, drawn
uniformly from `[1 - band, 1 + band]` with `Random.Xoshiro(seed)`. The same
`seed`, `N` and bus count give the same factors.
"""
function random_load_scenarios(net::Net, N::Integer; band::Real = 0.2, seed::Integer = 1)
  N >= 1 || throw(ArgumentError("random_load_scenarios: N must be at least 1, got $(N)"))
  0 <= band < 1 || throw(ArgumentError("random_load_scenarios: band must lie in [0, 1), got $(band)"))
  rng = Random.Xoshiro(seed)
  n = length(net.nodeVec)
  f = 1.0 .+ Float64(band) .* (2.0 .* rand(rng, n, Int(N)) .- 1.0)
  return InjectionScenarios(f, :random, Int(seed))
end

"""
    load_ramp_scenarios(net, N = 96; low = 0.8, high = 1.2) -> InjectionScenarios

`N` scenarios with one load factor for all buses, rising in equal steps
from `low` (scenario 1) to `high` (scenario `N`): a stand-in for a QSTS
day profile with 96 quarter hours.
"""
function load_ramp_scenarios(net::Net, N::Integer = 96; low::Real = 0.8, high::Real = 1.2)
  N >= 1 || throw(ArgumentError("load_ramp_scenarios: N must be at least 1, got $(N)"))
  n = length(net.nodeVec)
  steps = N == 1 ? [Float64(low)] : [Float64(low) + (Float64(high) - Float64(low)) * (s - 1) / (N - 1) for s = 1:N]
  return InjectionScenarios(repeat(reshape(steps, 1, :), n, 1), :ramp, nothing)
end

"""
    InjectionBatchResult

The result of [`runInjectionBatch`](@ref), one column or entry per
scenario: `vm` (pu) and `va_deg` (n x N, internal bus order), `converged`,
`iterations` (mismatch evaluations, the start included, as the length of
the solver's mismatch history), `final_mismatch` (pu), and
`qlimit_violations`, the number of PV buses whose reactive output at the
solution lies outside its limits by more than `qlimit_tol_pu` (Q limits are
off in batch mode). `backend` is `:cpu` or `:gpu`, `strategy` the GPU
linear-solve strategy (`:none` on the CPU), `note` the reason when a GPU
request ran on the CPU, `stats` counts and times of the run (seconds;
`setup_s` and `solve_s`; on the CPU the `executor` and its `workers`,
`worker_copy_s`, the creation of the working copies, and `worker_warm_s`,
the base solve of the workers in power mode, both inside `solve_s` and
summed over the workers (on threads and serially; in worker processes both
are part of `transfer_s`), and per scenario the `runpf!` time
`cpu_runpf_s`; on the GPU the device phases, the number of cuDSS analyses
and refactorizations, the cuDSS library version and reordering
algorithm).
"""
struct InjectionBatchResult
  vm::Matrix{Float64}
  va_deg::Matrix{Float64}
  converged::Vector{Bool}
  iterations::Vector{Int}
  final_mismatch::Vector{Float64}
  qlimit_violations::Vector{Int}
  backend::Symbol
  strategy::Symbol
  note::String
  stats::Dict{Symbol,Any}
end

"""
    InjectionBatchGpuInput

What the GPU batch of the CUDA extension needs, built once per batch on
the CPU: the Ybus, the start iterate `V0` shared by all scenarios, the
specified injections `S` (n x N), the bus types (no switching, Q limits
off), the voltage setpoints, the slack, the Newton settings and the
linear-solve `strategy` (`:sequential`, `:nonuniform`, `:uniform`), the
precision (`:float64`, `:float32_mixed`), `fast_mode` (cuDSS
deterministic mode off), `uniform_mask` (the uniform batch skips the
retired scenarios through `ubatch_mask`; off, it refactorizes and solves
all N every iteration), `fallback_to_cpu` (a GPU error returns its
message instead of throwing, so the batch can run on the CPU), and
`two_phase` (`:off`, `:g5`, `:g5r`, `:g6`) with `switch_tol`, the Float32
mismatch norm at which the two-phase variants leave phase 1 (see
[`runInjectionBatch`](@ref)). `ysh` is the conjugate admittance of the
voltage-dependent bus shunts per bus (zero without one): with such shunts
`S` holds the prosumer part only, and the device rebuilds
`S - abs(V)^2 * ysh` and the chain terms in every iteration, as
`buildControlledSVec` does on the CPU.
"""
struct InjectionBatchGpuInput
  Ybus::SparseMatrixCSC{ComplexF64,Int}
  V0::Vector{ComplexF64}
  S::Matrix{ComplexF64}
  bus_types::Vector{Symbol}
  Vset::Vector{Float64}
  slack_idx::Int
  maxiter::Int
  tol::Float64
  newton_update::Symbol
  strategy::Symbol
  precision::Symbol
  fast_mode::Bool
  ysh::Vector{ComplexF64}
  uniform_mask::Bool
  fallback_to_cpu::Bool
  two_phase::Symbol
  switch_tol::Float64
end

# the batch configuration: the effective power-flow configuration of the
# net with Q limits off and the rescue ladder and the DC fallback off (the
# GPU path has neither; a failed scenario is reported, not retried)
function _injection_batch_config(net::Net, config)
  effective = _runpf_effective_config(net, config)
  pf = effective === nothing ? powerflow_config() : (effective isa SparlectraConfig ? effective.powerflow : effective)
  qlim = _copy_qlimits_with(pf.qlimits; ignore_q_limits = true)
  dc = pf.dc
  dc_off = DcPowerFlowConfig(; NamedTuple{fieldnames(DcPowerFlowConfig)}(getfield.(Ref(dc), fieldnames(DcPowerFlowConfig)))..., fallback = false)
  return _copy_powerflow_with(pf; qlimits = qlim, rescue = false, dc = dc_off)
end

# Scale the loads of `work` by the factors of one scenario: P and Q of every
# load at bus b times f[b], and the bus load sums by the same amounts (as
# the :scale patch does, apply.jl). Only scalar fields change, so
# _reset_scenario_worker! puts them back.
function _apply_load_scale!(work::Net, f::AbstractVector{Float64})
  nodes = work.nodeVec
  n = length(nodes)
  length(f) == n || throw(DimensionMismatch("load_scale has $(length(f)) rows, the net $(n) buses"))
  for ps in work.prosumpsVec
    isGenerator(ps) && continue
    bus = getPosumerBusIndex(ps)
    (1 <= bus <= n) || continue
    fb = f[bus]
    if ps.pVal !== nothing
      p = ps.pVal
      ps.pVal = p * fb
      nodes[bus]._pƩLoad = something(nodes[bus]._pƩLoad, 0.0) + (p * fb - p)
    end
    if ps.qVal !== nothing
      q = ps.qVal
      ps.qVal = q * fb
      nodes[bus]._qƩLoad = something(nodes[bus]._qƩLoad, 0.0) + (q * fb - q)
    end
  end
  return work
end

# Undo _apply_load_scale! on `work`: the load P and Q of the prosumers and
# the bus load sums back to the template values (the only fields it
# changes), cheaper than the full worker reset when nothing was solved.
function _restore_load_scale!(work::Net, template::Net)
  for (w, t) in zip(work.prosumpsVec, template.prosumpsVec)
    isGenerator(t) && continue
    w.pVal = t.pVal
    w.qVal = t.qVal
  end
  for (w, t) in zip(work.nodeVec, template.nodeVec)
    w._pƩLoad = t._pƩLoad
    w._qƩLoad = t._qƩLoad
  end
  return work
end

# What the Q-limit report of every scenario needs, from a copy of the net
# (getQLimits_pu fills the limit tables of the net it reads): the Ybus as
# the solver builds it, the PV buses at the start and their limits in pu.
struct _BatchQLimitData
  Ybus::SparseMatrixCSC{ComplexF64,Int}
  pv::Vector{Int}
  qmin::Vector{Float64}
  qmax::Vector{Float64}
  tol::Float64
end

function _batch_ybus(net::Net)
  n = length(net.nodeVec)
  Yred = createYBUS(net = net, sparse = true, printYBUS = false)
  return size(Yred, 1) == n ? Yred : _expand_ybus_for_isolated_nodes(Yred, n, net.isoNodes)
end

function _BatchQLimitData(probe::Net, tol::Float64)
  qmin, qmax = getQLimits_pu(probe)
  pv = [k for (k, nd) in enumerate(probe.nodeVec) if getNodeType(nd) == PV]
  return _BatchQLimitData(_batch_ybus(probe), pv, copy(qmin), copy(qmax), tol)
end

# PV buses of one scenario whose reactive output (calculated injection plus
# load, as the Q-limit handler forms it) lies outside its limits by more
# than the tolerance
function _batch_qlimit_violations(q::_BatchQLimitData, V::AbstractVector{ComplexF64}, qload::Vector{Float64})
  Scalc = calc_injections(q.Ybus, V)
  count = 0
  for b in q.pv
    qreq = imag(Scalc[b]) + qload[b]
    (qreq > q.qmax[b] + q.tol || qreq < q.qmin[b] - q.tol) && (count += 1)
  end
  return count
end

# one scenario on a CPU worker: scale, runpf! with the batch configuration,
# read the state, reset the worker to the template. `start_v`, if given,
# replaces the start voltages of the worker (phase 2 of the two-phase
# variant :g6 starts from the phase-1 iterate of the GPU).
function _injection_cpu_one(work::Net, template::Net, scen::InjectionScenarios, cfg::PowerFlowConfig, qdata::_BatchQLimitData, s::Int; start_v::Union{Nothing,AbstractVector{ComplexF64}} = nothing)
  t0 = time_ns()
  _apply_load_scale!(work, view(scen.load_scale, :, s))
  start_v === nothing || _restore_start_voltages!(work, (abs.(start_v), rad2deg.(angle.(start_v))))
  qload = build_qload_pu(work)
  converged = false
  t_solve = 0.0
  try
    t1 = time_ns()
    runpf!(work; config = cfg)
    t_solve = (time_ns() - t1) / 1.0e9
    st = rectangular_pf_status(work)
    converged = st !== nothing && st.final_converged === true
  catch err
    # the one named failure: on a net with more than one AC island the
    # island-wise solve throws an ErrorException for an island that does
    # not converge or has no reference (a single island returns erg = 1
    # instead); the island status stays on `work` and the scenario is
    # reported as not converged. Anything else is rethrown.
    err isa ErrorException || rethrow()
  end
  st = rectangular_pf_status(work)
  n = length(work.nodeVec)
  vm = Vector{Float64}(undef, n)
  va = Vector{Float64}(undef, n)
  for (k, nd) in enumerate(work.nodeVec)
    vm[k] = something(nd._vm_pu, NaN)
    va[k] = something(nd._va_deg, NaN)
  end
  hist = st === nothing ? Float64[] : st.mismatch_history
  fm = st === nothing ? NaN : Float64(st.final_mismatch)
  V = vm .* cis.(deg2rad.(va))
  nviol = converged ? _batch_qlimit_violations(qdata, V, qload) : -1
  _reset_scenario_worker!(work, template)
  return (vm = vm, va = va, converged = converged, iterations = length(hist), final_mismatch = fm, qviol = nviol, solve_s = t_solve, item_s = (time_ns() - t0) / 1.0e9)
end

# Run f() with the session fields of `net` (solver status, DC status, power
# cache) set aside, put back afterwards: the working copies of a batch then
# carry no factorization (so they can be made concurrently, and serialized
# for a worker process) and no status of an earlier solve. Every scenario
# solve writes its own status.
function _with_session_fields_aside(f, net::Net)
  saved = (net._rectangular_pf_status, net._dc_pf_status, net._power_cache)
  net._rectangular_pf_status = nothing
  net._dc_pf_status = nothing
  net._power_cache = nothing
  try
    return f()
  finally
    net._rectangular_pf_status, net._dc_pf_status, net._power_cache = saved
  end
end

# Why the GPU path cannot take this batch (empty string: it can). The GPU
# batch covers the plain Newton of one island with one slack: no links,
# no isolated buses, no Q(U) or P(U) controller with a voltage-dependent
# characteristic (constant ones and voltage-dependent bus shunts are
# covered), no distributed slack,
# no damping, start projection or other start helpers, no Jacobian reuse.
function _injection_gpu_refusal(probe::Net, cfg::PowerFlowConfig, Ybus::SparseMatrixCSC{ComplexF64,Int})::String
  backend = _GPU_BACKEND[]
  (backend !== nothing && hasproperty(backend, :batch)) || return "the CUDA extension is not loaded (using CUDA, CUDSS)"
  status = backend.status()
  status.status === :available || return string("no usable GPU: ", status.reason)
  any(lk -> lk.status != 0, probe.linkVec) && return "the network has links in service"
  isempty(probe.isoNodes) || return "the network has isolated buses"
  size(Ybus, 1) == length(probe.nodeVec) || return "the Ybus does not cover every bus"
  count(nd -> getNodeType(nd) == Slack, probe.nodeVec) == 1 || return "the network has more or less than one slack bus"
  any(nd -> getNodeType(nd) == Isolated, probe.nodeVec) && return "the network has isolated buses"
  _ybus_connected(Ybus) || return "the network has more than one island"
  any(ps -> (has_qu_controller(ps) && !_constant_controller(ps.quController)) || (has_pu_controller(ps) && !_constant_controller(ps.puController)), probe.prosumpsVec) && return "Q(U) and P(U) controllers with a voltage-dependent characteristic are not in the GPU batch"
  cfg.distributed_slack.enabled && return "distributed slack is not in the GPU batch"
  (cfg.autodamp || cfg.merit.enabled || cfg.trust_region.enabled) && return "damping (autodamp, merit, trust region) runs on the CPU"
  cfg.start_mode.start_projection && return "the start projection runs on the CPU"
  (cfg.start_current_iteration.enabled || cfg.apslf_start.enabled) && return "the start helpers (current iteration, APSLF start) run on the CPU"
  cfg.jacobian_reuse && return "Jacobian reuse is not in the GPU batch"
  return ""
end

# A controller whose characteristic has one value at every point gives that
# value at any voltage with slope 0 (linear and spline interpolation; the
# MATPOWER import maps PQ generator limits to such controllers), so its
# injection is a constant.
_constant_controller(ctrl) = ctrl.characteristic.interpolation !== :polynomial && all(p -> p[2] == ctrl.characteristic.points[1][2], ctrl.characteristic.points)

# every bus reachable from bus 1 through the Ybus pattern
function _ybus_connected(Y::SparseMatrixCSC)
  n = size(Y, 1)
  n == 0 && return true
  seen = falses(n)
  stack = [1]
  seen[1] = true
  while !isempty(stack)
    j = pop!(stack)
    for k in nzrange(Y, j)
      i = rowvals(Y)[k]
      seen[i] && continue
      seen[i] = true
      push!(stack, i)
    end
  end
  return all(seen)
end

# The GPU input: start iterate, bus types and setpoints as runpf_rectangular!
# forms them (initialVrect, _bus_voltage_setpoints_from_prosumers, the slack
# magnitude fix, the auxiliary-bus flat start), and S per scenario from
# buildComplexSVec on the scaled probe, so each column is the vector the CPU
# solve of that scenario builds.
function _injection_gpu_input(probe::Net, template::Net, scen::InjectionScenarios, cfg::PowerFlowConfig, Ybus::SparseMatrixCSC{ComplexF64,Int}, tol::Float64, strategy::Symbol, uniform_mask::Bool, two_phase::Symbol, switch_tol::Float64, qloads::Matrix{Float64})
  n = length(probe.nodeVec)
  start = cfg.start_mode
  V0, slack_idx = initialVrect(probe; flatstart = start.flatstart, rotate_behind_pst = start.rotate_behind_pst)
  Vset = _bus_voltage_setpoints_from_prosumers(probe)
  bus_types = Symbol[getNodeType(nd) == Slack ? :Slack : getNodeType(nd) == PV ? :PV : :PQ for nd in probe.nodeVec]
  V0[slack_idx] = _apply_voltage_magnitude_preserving_angle(V0[slack_idx], Vset[slack_idx])
  start.flatstart && _aux_bus_flat_start!(V0, probe, slack_idx; rotate_behind_pst = start.rotate_behind_pst)
  N = length(scen)
  # voltage-dependent bus shunts (with constant controllers the only
  # voltage-dependent injection the GPU batch covers): their conjugate
  # admittance summed per bus
  ysh = zeros(ComplexF64, n)
  for sh in probe.shuntVec
    (sh.status == 0 || sh.model != :VoltageDependentInjection || !(1 <= sh.busIdx <= n) || sh.busIdx in probe.isoNodes) && continue
    ysh[sh.busIdx] += conj(sh.y_pu_shunt)
  end
  vdep_ctrl = any(ps -> has_qu_controller(ps) || has_pu_controller(ps), probe.prosumpsVec)
  S = Matrix{ComplexF64}(undef, n, N)
  for s = 1:N
    _apply_load_scale!(probe, view(scen.load_scale, :, s))
    # the prosumer part: with controllers (constant, see the refusal) their
    # value replaces P or Q as in buildControlledSVec, which also subtracts
    # the shunts at V0; buildComplexSVec subtracts them at abs(V) = 1.
    # Added back, the column is the prosumer part (to rounding).
    S[:, s] = vdep_ctrl ? buildControlledSVec(probe, V0)[1] .+ abs.(V0) .^ 2 .* ysh : buildComplexSVec(probe) .+ ysh
    qloads[:, s] = build_qload_pu(probe)
    _restore_load_scale!(probe, template)
  end
  cb = cfg.compute_backend
  return InjectionBatchGpuInput(Ybus, V0, S, bus_types, Vset, slack_idx, cfg.max_iter, tol, cfg.newton_update, strategy, cb.gpu_precision, cb.gpu_fast_mode, ysh, uniform_mask, cb.gpu_fallback_to_cpu, two_phase, switch_tol)
end

# `cfg` with compute_backend.mode = :cpu: the solves of the batch that run
# through runpf! (the CPU scenarios, the base solve of start = :base, phase 2
# of :g6) stay on the CPU whatever mode the caller set, since mode = :gpu
# would hand each of them to the single-run cuDSS context
_injection_cpu_config(cfg::PowerFlowConfig) = cfg.compute_backend.mode === :cpu ? cfg : _copy_powerflow_with(cfg; compute_backend = _struct_with(cfg.compute_backend; mode = :cpu))

# the configuration of a solve that starts from a converged state: flat
# start, start projection, DC seed and current iteration off, as the QSTS
# steps after a converged base (_qsts_step_config)
function _injection_warm_config(cfg::PowerFlowConfig)
  sm = cfg.start_mode
  return _copy_powerflow_with(cfg;
    start_mode = _copy_start_mode_with(sm; flatstart = false, start_projection = false, dc_seed_unconditional = false),
    start_current_iteration = _copy_start_current_iteration_with(cfg.start_current_iteration; enabled = false))
end

# start = :base: solve the base case (the loads of `net`) once with the
# batch configuration and start every scenario from its solution, as the
# QSTS steps after a converged base do (_qsts_step_config: flat start,
# start projection, DC seed and current iteration off). The base solve runs
# on the CPU for every backend: with compute_backend.mode = :gpu it would
# run in the single-run cuDSS context, and the GPU variants would start
# from another base iterate than the CPU batch. Returns the template with
# the base voltages as start values, the warm configuration (compute
# backend as given) and the base solve time and iteration count.
function _injection_warm_template(net::Net, cfg::PowerFlowConfig)
  base = deepcopy(net)
  converged = false
  failure = ""
  t0 = time_ns()
  try
    runpf!(base; config = _injection_cpu_config(cfg))
    st = rectangular_pf_status(base)
    converged = st !== nothing && st.final_converged === true
  catch err
    # the named failure of the island-wise solve (see _injection_cpu_one);
    # its text goes into the error below
    err isa ErrorException || rethrow()
    failure = string(": ", first(split(sprint(showerror, err), '\n')))
  end
  base_s = (time_ns() - t0) / 1.0e9
  converged || error("runInjectionBatch: start = :base needs a converged base case, and the base case did not converge with the batch configuration", failure)
  st = rectangular_pf_status(base)
  template = deepcopy(net)
  _restore_start_voltages!(template, _snapshot_start_voltages(base))
  return template, _injection_warm_config(cfg), base_s, length(st.mismatch_history)
end

# Power mode on the CPU: every worker solves the base case once before its
# first scenario, so its cache (Ybus, LU analysis, workspace, the auto LU
# decision) is warm from the first scenario on, as in a long QSTS run. On
# threads the workers warm in their chunk tasks, side by side (distinct
# nets; a copied net holds no factorization, see the PowerModeCache
# deepcopy), in a worker process once when it loads the batch. A base case
# that fails leaves the worker cold; the scenarios then pay the analysis.
function _warm_power_mode_worker!(work::Net, template::Net, cfg::PowerFlowConfig)
  try
    runpf!(work; config = cfg)
  catch err
    # the named failure of the island-wise solve (see _injection_cpu_one);
    # the worker stays cold (see above)
    err isa ErrorException || rethrow()
  end
  _reset_scenario_worker!(work, template)
  return nothing
end

# the process batch: the worker state is the template, its working copy
# (warmed in power mode), the batch configuration and the Q-limit data; a
# scenario gets its column of load factors and, for phase 2 of :g6, its
# start voltages (else nothing)
function _injection_process_setup(data)
  template = data.net
  work = deepcopy(template)
  data.cfg.power_mode && _warm_power_mode_worker!(work, template, data.cfg)
  return (template = template, work = work, cfg = data.cfg, qdata = data.qdata)
end

function _injection_process_one(state, s::Int, item::Tuple{Vector{Float64},Union{Nothing,Vector{ComplexF64}}})
  f, v = item
  scen = InjectionScenarios(reshape(f, :, 1), :custom, nothing)
  return _injection_cpu_one(state.work, state.template, scen, state.cfg, state.qdata, 1; start_v = v)
end

# The executor of `nitems` CPU scenarios on a net of `nbus` buses
# (_resolve_batch_executor) as (mode, workers, process ids), the process
# pool started when it is needed; the pool start counts into the time of
# the caller.
function _injection_executor!(stats::Dict{Symbol,Any}, nitems::Int, nbus::Int, executor, max_processes, parallel_enabled, parallel_max_tasks, parallel_min_work_items)
  mode, k = _resolve_batch_executor(executor, nitems, nbus; parallel_enabled = parallel_enabled, parallel_max_tasks = parallel_max_tasks, parallel_min_work_items = parallel_min_work_items, max_processes = max_processes)
  stats[:executor] = mode
  stats[:workers] = k
  pids = Int[]
  if mode === :processes
    t_pool = time_ns()
    pids = start_batch_processes!(k)[1:k]
    stats[:pool_start_s] = (time_ns() - t_pool) / 1.0e9
  end
  return (mode = mode, k = k, pids = pids)
end

# The scenarios `idx` on the CPU with the executor `ex`
# (_injection_executor!): in worker processes (the template sent once per
# worker), or on threads or serially with one working copy per chunk, made
# in its task; power mode workers are warmed first. The session fields of
# the template are set aside meanwhile, so no copy holds a factorization.
# BLAS runs at `blas_threads` for the duration (in the workers, or restored
# afterwards). `start_v` (n x N) gives each scenario its own start
# voltages. Process stats go into `stats`. Returns the rows in the order of
# `idx`, the summed times of the worker copies and of the warm-ups (threads
# and serial; in processes both are part of `transfer_s`) and the BLAS
# thread count before the batch.
function _injection_cpu_map(template::Net, scen::InjectionScenarios, cfg::PowerFlowConfig, qdata::_BatchQLimitData, idx::Vector{Int}, ex, stats::Dict{Symbol,Any}, blas_threads; start_v::Union{Nothing,Matrix{ComplexF64}} = nothing)
  blas_before = BLAS.get_num_threads()
  copy_s = Threads.Atomic{Float64}(0.0)
  warm_s = Threads.Atomic{Float64}(0.0)
  rows = _with_session_fields_aside(template) do
    if ex.mode === :processes
      t_ser = time_ns()
      payload = _serialize_payload((net = template, cfg = cfg, qdata = qdata))
      stats[:serialize_s] = (time_ns() - t_ser) / 1.0e9
      stats[:payload_mb] = length(payload) / 2^20
      item = j -> (scen.load_scale[:, idx[j]], start_v === nothing ? nothing : start_v[:, idx[j]])
      prows, pstats = _process_parallel_map(_injection_process_setup, payload, _injection_process_one, item, length(idx), ex.pids; blas_threads = blas_threads)
      merge!(stats, pstats)
      return prows
    end
    prepare = cfg.power_mode ? (w -> (local tw = time_ns(); _warm_power_mode_worker!(w.net, template, cfg); Threads.atomic_add!(warm_s, (time_ns() - tw) / 1.0e9))) : nothing
    try
      blas_threads === nothing || BLAS.set_num_threads(blas_threads)
      # each chunk task makes its own working copy and warms it: the session
      # fields are set aside, so no copy touches a factorization and the
      # copies may run at the same time. copy_s and warm_s sum their times.
      return _engine_parallel_map((worker, j) -> _injection_cpu_one(worker.net, template, scen, cfg, qdata, idx[j]; start_v = start_v === nothing ? nothing : view(start_v, :, idx[j])),
        () -> (local tw = time_ns(); local w = ScenarioWorker(deepcopy(template)); Threads.atomic_add!(copy_s, (time_ns() - tw) / 1.0e9); w), length(idx);
        parallel_enabled = ex.mode === :threads, parallel_max_tasks = ex.k, parallel_min_work_items = 1, workers_in_tasks = true, prepare = prepare)
    finally
      blas_threads === nothing || BLAS.set_num_threads(blas_before)
    end
  end
  return rows, copy_s[], warm_s[], blas_before
end

# Phase 2 of :g6: the scenarios whose phase 1 on the GPU reached the switch
# (tolerance or stagnation) finish with runpf! on the CPU from the phase-1
# iterate, with the warm configuration and on the CPU in any case; the
# others keep the phase-1 result as not converged. Iterations add up the
# mismatch evaluations of both phases.
function _injection_g6_finish(out, template::Net, scen::InjectionScenarios, cfg::PowerFlowConfig, qdata::_BatchQLimitData, stats::Dict{Symbol,Any}, t_solve::UInt64, blas_threads, executor, max_processes, parallel_enabled, parallel_max_tasks, parallel_min_work_items)
  N = length(scen)
  reason = out.stats[:tp_switch_reason]
  idx = [s for s = 1:N if reason[s] in (:tol, :stagnation)]
  warm = _injection_cpu_config(_injection_warm_config(cfg))
  t_cpu = time_ns()
  ex = _injection_executor!(stats, length(idx), length(template.nodeVec), executor, max_processes, parallel_enabled, parallel_max_tasks, parallel_min_work_items)
  rows, copy_s, warm_s, blas_before = _injection_cpu_map(template, scen, warm, qdata, idx, ex, stats, blas_threads; start_v = out.V)
  stats[:cpu_phase2_s] = (time_ns() - t_cpu) / 1.0e9
  stats[:worker_copy_s] = copy_s
  stats[:worker_warm_s] = warm_s
  stats[:blas_threads] = blas_threads === nothing ? blas_before : blas_threads
  vm = abs.(out.V)
  va = rad2deg.(angle.(out.V))
  converged = falses(N)
  iterations = copy(out.iterations)
  final_mismatch = copy(out.final_mismatch)
  qviol = fill(-1, N)
  p2 = zeros(Int, N)
  for (k, s) in enumerate(idx)
    r = rows[k]
    vm[:, s] = r.vm
    va[:, s] = r.va
    converged[s] = r.converged
    iterations[s] += r.iterations
    p2[s] = max(r.iterations - 1, 0)
    final_mismatch[s] = r.final_mismatch
    qviol[s] = r.qviol
  end
  stats[:tp_phase2_steps] = p2
  stats[:solve_s] = (time_ns() - t_solve) / 1.0e9
  return InjectionBatchResult(vm, va, Vector{Bool}(converged), iterations, final_mismatch, qviol, :gpu, :sequential, "", stats)
end

const _INJECTION_GPU_STRATEGIES = (:sequential, :nonuniform, :uniform)
const _INJECTION_STARTS = (:net, :base)
const _INJECTION_TWO_PHASE = (:off, :g5, :g5r, :g6)

"""
    runInjectionBatch(net, scenarios::InjectionScenarios; config = nothing, backend = :cpu, start = :net, gpu_strategy = :sequential, gpu_uniform_mask = true, gpu_two_phase = :off, gpu_switch_mismatch = 1e-5, qlimit_tol_pu = 1e-4, blas_threads = 1, executor = nothing, max_processes = nothing, parallel_enabled = nothing, parallel_max_tasks = nothing, parallel_min_work_items = nothing) -> InjectionBatchResult

Solve every scenario of `scenarios` on `net` and return the states, the
convergence flags and the Q-limit report ([`InjectionBatchResult`](@ref)).
`net` is not changed.

The configuration is the one `runpf!(net)` would use (`config`, else the
configuration the net was imported with, else the active one), with three
changes for the batch: Q limits off (each scenario reports instead how many
PV buses would leave their limits by more than `qlimit_tol_pu`), and the
rescue ladder and the DC fallback off (a failed scenario is reported, not
retried).

`start = :net` starts every scenario from the start iterate of `net` (its
imported voltages, or a flat start if the configuration asks for one).
`start = :base` first solves the base case (the loads of `net`) with the
batch configuration on the CPU, for either backend, so the CPU and the GPU
batch start from the same bits, and starts every scenario from that solution, with
flat start, start projection and the start helpers off, as the QSTS steps
after a converged base; the base solve counts into `stats[:setup_s]`.

`backend = :cpu` runs each scenario through `runpf!` on reused working
copies of `net` (reset after every scenario), in worker processes, on Julia
threads or serially: `executor` (`:auto`, `:threads`, `:processes`, else
`runtime.parallel.executor`) picks, `:auto` takes processes from 2000 buses
on and threads below ([Parallel Execution](@ref parallel_execution_processes)).
The process batch uses the worker pool of [`start_batch_processes!`](@ref)
(started on first use, at most `max_processes` or
`runtime.parallel.max_processes` workers) and sends the net once per
worker; the thread batch makes one working copy per chunk inside its task.
`parallel_enabled = false` or fewer scenarios than `min_work_items` run
serially. The results are the same bit for bit on every executor. With
`power_flow.power_mode` on, every worker solves the base case once before
its first scenario (`stats[:worker_warm_s]` sums these solves on threads),
so its power-mode caches (Ybus, LU analysis) are warm from its first
scenario on and kept across its scenarios; the copies start without the
cache and the solver status of `net`. `blas_threads` sets the BLAS thread
count of the batch (in the workers, or restored afterwards; `nothing`
leaves it). The GC is paused while a scenario solves (`runtime.gc_pause`).

`backend = :gpu` needs the CUDA extension (`using CUDA, CUDSS`): the Ybus,
the Jacobian pattern and one cuDSS analysis are built once, and all
scenarios run their Newton iterations on the device, the Jacobian of all
active scenarios in one kernel launch. `gpu_strategy` picks the linear
solve: `:sequential` (default; one cuDSS solver, refactorization and solve
per scenario, deterministic mode), `:nonuniform` (cuDSS non-uniform batch
on column-permuted Jacobians, deterministic mode) or `:uniform` (cuDSS
uniform batch on column-permuted Jacobians, fast mode only, since cuDSS
1.0.0 refuses the deterministic mode there; `gpu_uniform_mask = false`
refactorizes and solves the retired scenarios as well). The precision
follows `power_flow.compute_backend.gpu_precision`. If the GPU cannot take
the batch (no extension or device, more than one island, links, Q(U) or
P(U) controllers with a voltage-dependent characteristic, distributed
slack, damping or start helpers) or fails on the device (out of device
memory, a cuDSS error), the batch runs on the CPU and `note` says why,
as long as `power_flow.compute_backend.gpu_fallback_to_cpu` is true (the
default); with it false a refusal throws an `ArgumentError` and a device
error is rethrown.

`gpu_two_phase` (GPU and `:sequential` only, study variants) replaces the
precision setting by two phases. Phase 1 runs V, mismatch, Jacobian and LU
in Float32 on the device until the Float32 mismatch norm is at most
`gpu_switch_mismatch`, or until two steps in a row reduce it by less than
half. Phase 2: `:g5` keeps V and mismatch in Float64 on the device and
steps with the last Float32 LU of phase 1, with one Float32
refactorization after 5 steps without convergence; `:g5r` refactorizes in
Float32 at the first step that reduces the mismatch by less than half and
goes over to a Float64 Newton on the device at the next such step; `:g6`
finishes with a Float64 Newton on the CPU (`runpf!` from the phase-1
iterate, flat start and start helpers off, on the CPU workers as with
`backend = :cpu`). `stats` holds per scenario the steps of both phases,
the reason of the switch, the refactorizations of phase 2 and the
reduction rate per step (`:tp_*` keys). `iterations` counts the mismatch
evaluations of both phases (for `:g6` the CPU history, its start
included).

Throws an `ArgumentError` for an unknown backend, start, strategy or
two-phase variant, for `gpu_two_phase` with the CPU backend or another
strategy than `:sequential`, or for a scenario matrix that does not match
the bus count, and an error when `start = :base` and the base case does
not converge.
"""
function runInjectionBatch(
  net::Net,
  scenarios::InjectionScenarios;
  config::Union{Nothing,PowerFlowConfig,SparlectraConfig} = nothing,
  backend::Symbol = :cpu,
  start::Symbol = :net,
  gpu_strategy::Symbol = :sequential,
  gpu_uniform_mask::Bool = true,
  gpu_two_phase::Symbol = :off,
  gpu_switch_mismatch::Float64 = 1e-5,
  qlimit_tol_pu::Float64 = 1e-4,
  blas_threads::Union{Nothing,Int} = 1,
  executor::Union{Nothing,Symbol} = nothing,
  max_processes::Union{Nothing,Int} = nothing,
  parallel_enabled::Union{Nothing,Bool} = nothing,
  parallel_max_tasks::Union{Nothing,Int} = nothing,
  parallel_min_work_items::Union{Nothing,Int} = nothing,
)
  backend in (:cpu, :gpu) || throw(ArgumentError("runInjectionBatch: backend must be :cpu or :gpu, got $(backend)"))
  start in _INJECTION_STARTS || throw(ArgumentError("runInjectionBatch: start must be one of $(_INJECTION_STARTS), got $(start)"))
  gpu_strategy in _INJECTION_GPU_STRATEGIES || throw(ArgumentError("runInjectionBatch: gpu_strategy must be one of $(_INJECTION_GPU_STRATEGIES), got $(gpu_strategy)"))
  gpu_two_phase in _INJECTION_TWO_PHASE || throw(ArgumentError("runInjectionBatch: gpu_two_phase must be one of $(_INJECTION_TWO_PHASE), got $(gpu_two_phase)"))
  gpu_two_phase === :off || (backend === :gpu && gpu_strategy === :sequential) || throw(ArgumentError("runInjectionBatch: gpu_two_phase needs backend = :gpu and gpu_strategy = :sequential"))
  n = length(net.nodeVec)
  size(scenarios.load_scale, 1) == n || throw(ArgumentError("runInjectionBatch: load_scale has $(size(scenarios.load_scale, 1)) rows, the net $(n) buses"))
  N = length(scenarios)
  t_setup = time_ns()
  stats = Dict{Symbol,Any}(:scenarios => N, :start => start)
  # one assignment each: the worker closures below capture them
  template, cfg, base_s, base_iterations = start === :base ? _injection_warm_template(net, _injection_batch_config(net, config)) : (net, _injection_batch_config(net, config), 0.0, 0)
  stats[:base_solve_s] = base_s
  stats[:base_iterations] = base_iterations
  tol = _resolve_tolerance_for_net(cfg, template).tol
  stats[:tol] = tol
  # one copy for the read-only data of the batch (Ybus, Q limits, GPU input)
  probe = deepcopy(template)
  qdata = _BatchQLimitData(probe, qlimit_tol_pu)
  note = ""
  if backend === :gpu
    fallback = cfg.compute_backend.gpu_fallback_to_cpu
    why = _injection_gpu_refusal(probe, cfg, qdata.Ybus)
    isempty(why) || fallback || throw(ArgumentError(string("runInjectionBatch: backend = :gpu with power_flow.compute_backend.gpu_fallback_to_cpu = false: ", why)))
    if isempty(why)
      qloads = Matrix{Float64}(undef, n, N)
      input = _injection_gpu_input(probe, template, scenarios, cfg, qdata.Ybus, tol, gpu_strategy, gpu_uniform_mask, gpu_two_phase, gpu_switch_mismatch, qloads)
      stats[:setup_s] = (time_ns() - t_setup) / 1.0e9
      t_solve = time_ns()
      out = _GPU_BACKEND[].batch(input)
      if hasproperty(out, :error)
        # a device error with gpu_fallback_to_cpu on (the extension rethrows
        # it otherwise); the batch runs on the CPU below
        why = string("the GPU batch failed: ", out.error)
        stats[:gpu_failed_s] = (time_ns() - t_solve) / 1.0e9
      else
        merge!(stats, out.stats)
        gpu_two_phase === :g6 && return _injection_g6_finish(out, template, scenarios, cfg, qdata, stats, t_solve, blas_threads, executor, max_processes, parallel_enabled, parallel_max_tasks, parallel_min_work_items)
        vm = abs.(out.V)
        va = rad2deg.(angle.(out.V))
        qviol = [out.converged[s] ? _batch_qlimit_violations(qdata, view(out.V, :, s), qloads[:, s]) : -1 for s = 1:N]
        stats[:solve_s] = (time_ns() - t_solve) / 1.0e9
        return InjectionBatchResult(vm, va, out.converged, out.iterations, out.final_mismatch, qviol, :gpu, gpu_strategy, "", stats)
      end
    end
    note = string("GPU batch not used, the batch ran on the CPU: ", why)
    @warn note
  end
  # the CPU path solves on the CPU whatever compute_backend.mode says (a GPU
  # request that fell back would otherwise run the single-run GPU path)
  cpu_cfg = _injection_cpu_config(cfg)
  ex = _injection_executor!(stats, N, n, executor, max_processes, parallel_enabled, parallel_max_tasks, parallel_min_work_items)
  stats[:setup_s] = (time_ns() - t_setup) / 1.0e9
  t_solve = time_ns()
  rows, copy_s, warm_s, blas_before = _injection_cpu_map(template, scenarios, cpu_cfg, qdata, collect(1:N), ex, stats, blas_threads)
  stats[:solve_s] = (time_ns() - t_solve) / 1.0e9
  stats[:worker_copy_s] = copy_s
  stats[:worker_warm_s] = warm_s
  stats[:cpu_runpf_s] = Float64[r.solve_s for r in rows]
  stats[:cpu_item_s] = Float64[r.item_s for r in rows]
  stats[:julia_threads] = Threads.nthreads()
  stats[:blas_threads] = blas_threads === nothing ? blas_before : blas_threads
  vm = Matrix{Float64}(undef, n, N)
  va = Matrix{Float64}(undef, n, N)
  for s = 1:N
    vm[:, s] = rows[s].vm
    va[:, s] = rows[s].va
  end
  return InjectionBatchResult(vm, va, Bool[r.converged for r in rows], Int[r.iterations for r in rows], Float64[r.final_mismatch for r in rows], Int[r.qviol for r in rows], :cpu, :none, note, stats)
end

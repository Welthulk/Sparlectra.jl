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

# file: src/scenario/process_pool.jl
# purpose: the process batch: a pool of Julia worker processes (Distributed
#          stdlib, one thread each) that runs the batches of the scenario
#          engine and the injection batch. Threads of one process serialize
#          on two locks inside the sparse LU (the OpenBLAS buffer lock and
#          the GC safepoint wait of the other threads); a process has its
#          own of both. The pool lives until it is stopped or the process
#          ends; a batch sends its data once per worker.

# Distributed is loaded on the first use of the pool, not with the package
# (it adds load time every session would pay). Its functions are newer than
# the calling code, hence invokelatest.
const _DISTRIBUTED_PKGID = Base.PkgId(Base.UUID("8ba89e20-285c-5b6f-9357-94700520ee1b"), "Distributed")

function _distributed()::Module
  m = get(Base.loaded_modules, _DISTRIBUTED_PKGID, nothing)
  m === nothing || return m
  return Base.require(_DISTRIBUTED_PKGID)
end

_dcall(f::Symbol, args...; kwargs...) = Base.invokelatest(getfield(_distributed(), f), args...; kwargs...)
_dcall_type(name::Symbol) = getfield(_distributed(), name)

# the worker processes of the pool (process ids), only those it started
const _BATCH_POOL = Int[]
const _BATCH_POOL_LOCK = ReentrantLock()
# a batch key per batch, so the state of two batches on one worker never mixes
const _BATCH_KEY = Threads.Atomic{Int}(0)

# The command line of a worker: the image, the project, the load and depot
# path and the code options of this process, so a worker loads the same
# Sparlectra and compiles the same code (bit-identical results); one thread.
function _batch_worker_exeflags()::Cmd
  opts = Base.JLOptions()
  flags = String["--startup-file=no", "-t", "1", "--color=no", "-O$(opts.opt_level)"]
  push!(flags, "-J", unsafe_string(opts.image_file))
  opts.check_bounds == 1 && push!(flags, "--check-bounds=yes")
  opts.check_bounds == 2 && push!(flags, "--check-bounds=no")
  project = Base.active_project()
  project === nothing || push!(flags, "--project=$(project)")
  return Cmd(flags)
end

function _batch_worker_env()::Vector{Pair{String,String}}
  sep = Sys.iswindows() ? ";" : ":"
  return ["JULIA_LOAD_PATH" => join(LOAD_PATH, sep), "JULIA_DEPOT_PATH" => join(DEPOT_PATH, sep)]
end

# The extensions of this process that change the numbers of a solve, by
# extension name and trigger package: with KLU loaded, power mode
# factorizes with KLU (power_mode_linear_solver_backend), and a worker
# without it would solve with UMFPACK and round differently from the serial
# and the thread batch. Their trigger packages are loaded on every worker of
# the pool when this process has the extension (a worker that has them
# already returns at once); the CUDA extension is left out on purpose, the
# batch workers solve on the CPU.
const _BATCH_WORKER_EXTENSIONS = [(:SparlectraKLUExt, Base.PkgId(Base.UUID("ef3ab10e-7fda-4108-b977-705223b18434"), "KLU"))]

function _load_batch_worker_extensions!(pids::Vector{Int})
  isempty(pids) && return nothing
  for (ext, pkg) in _BATCH_WORKER_EXTENSIONS
    Base.get_extension(@__MODULE__, ext) === nothing && continue
    _dcall(:remotecall_eval, Main, pids, :(Base.require($(pkg))))
  end
  return nothing
end

"""
    start_batch_processes!(n::Integer) -> Vector{Int}

Start worker processes for the process batch until the pool has at least `n`
(one Julia thread each, the project, load path and system image of this
process), load Sparlectra on them and return the process ids of the pool.
The pool lives until [`stop_batch_processes!`](@ref) or the end of the
process. A batch that runs in processes starts the pool by itself; calling
this first keeps the start (a few seconds) out of the first batch. When
this process has the KLU extension loaded, every worker of the pool loads
KLU as well (power mode then factorizes with KLU on both sides).
Throws an `ErrorException` when a worker loads another Sparlectra than this
process; the workers started by a call that throws are stopped again.
"""
function start_batch_processes!(n::Integer)
  n >= 1 || throw(ArgumentError("start_batch_processes!: n must be at least 1, got $(n)"))
  lock(_BATCH_POOL_LOCK) do
    alive = Set{Int}(_dcall(:procs))
    filter!(in(alive), _BATCH_POOL)
    missing_n = n - length(_BATCH_POOL)
    if missing_n > 0
      pids = _dcall(:addprocs, missing_n; exeflags = _batch_worker_exeflags(), env = _batch_worker_env())
      # the new workers join the pool only when they loaded this
      # Sparlectra; on any failure before that they are removed again, so
      # no worker outlives the call unseen by stop_batch_processes!
      joined = false
      try
        pkg = Base.PkgId(@__MODULE__)
        _dcall(:remotecall_eval, Main, pids, :(Base.require($(pkg))))
        here = pathof(@__MODULE__)
        for p in pids
          there = _dcall(:remotecall_fetch, () -> pathof(Base.root_module(pkg)), p)
          there == here || error("start_batch_processes!: worker $(p) loaded Sparlectra from $(there), this process from $(here)")
        end
        append!(_BATCH_POOL, pids)
        joined = true
      finally
        joined || _dcall(:rmprocs, pids; waitfor = 30.0)
      end
    end
    _load_batch_worker_extensions!(_BATCH_POOL)
    return copy(_BATCH_POOL)
  end
end

"""
    stop_batch_processes!() -> Nothing

Stop the worker processes of the process batch (see
[`start_batch_processes!`](@ref)). Worker processes started otherwise are
not touched.
"""
function stop_batch_processes!()
  lock(_BATCH_POOL_LOCK) do
    if !isempty(_BATCH_POOL)
      alive = Set{Int}(_dcall(:procs))
      pids = filter(in(alive), _BATCH_POOL)
      isempty(pids) || _dcall(:rmprocs, pids; waitfor = 30.0)
      empty!(_BATCH_POOL)
    end
  end
  return nothing
end

"""
    batch_process_ids() -> Vector{Int}

The process ids of the running workers of the process batch (empty when
the pool is not started). A worker that has exited is not listed (and
leaves the pool).
"""
function batch_process_ids()
  lock(_BATCH_POOL_LOCK) do
    # Distributed is loaded with the first worker; without workers there
    # is nothing to check
    isempty(_BATCH_POOL) || filter!(in(Set{Int}(_dcall(:procs))), _BATCH_POOL)
    return copy(_BATCH_POOL)
  end
end

# --- worker side ---------------------------------------------------------------

# the state of the running batches on a worker process, by batch key
const _PROCESS_WORKER_STATE = Dict{Int,Any}()

# The worker state of the last batch, kept for the next batch with the same
# payload (same digest, same setup function): a batch on the same net then
# sends and loads nothing. Every item resets its working copy to the
# template, so a kept state equals a new one. A state is used by one batch
# at a time (`busy`); a second batch with the same payload at the same time
# builds its own.
mutable struct _ProcessWorkerCache
  digest::Vector{UInt8}
  setup::Any
  state::Any
  busy::Bool
end

const _PROCESS_WORKER_CACHE = Ref{Union{Nothing,_ProcessWorkerCache}}(nothing)

function _process_worker_settings!(blas_threads::Union{Nothing,Int}, gc_settings::Tuple{Bool,Int})
  blas_threads === nothing || BLAS.set_num_threads(blas_threads)
  _GC_PAUSE_OVERRIDE[] = gc_settings
  return nothing
end

# Take the kept state for this batch when it was built from the same
# payload and no other batch uses it. Runs on a worker process.
function _process_worker_reuse!(key::Int, @nospecialize(setup), digest::Vector{UInt8}, blas_threads::Union{Nothing,Int}, gc_settings::Tuple{Bool,Int})::Bool
  c = _PROCESS_WORKER_CACHE[]
  (c !== nothing && !c.busy && c.setup === setup && c.digest == digest) || return false
  c.busy = true
  _process_worker_settings!(blas_threads, gc_settings)
  _PROCESS_WORKER_STATE[key] = c.state
  return true
end

# Load the data of one batch on a worker: deserialize the payload, let
# `setup` build the worker state from it, set the BLAS threads and the GC
# pause settings of the main process, and keep the state for the next batch
# with this payload (it replaces the kept one unless that one is in use).
# Runs on a worker process.
function _process_worker_setup!(key::Int, @nospecialize(setup), payload::Vector{UInt8}, digest::Vector{UInt8}, blas_threads::Union{Nothing,Int}, gc_settings::Tuple{Bool,Int})
  t0 = time_ns()
  c = _PROCESS_WORKER_CACHE[]
  # drop the kept state before the new one is built, so two nets are not
  # in memory at once
  c === nothing || c.busy || (_PROCESS_WORKER_CACHE[] = nothing; c = nothing; GC.gc())
  data = Serialization.deserialize(IOBuffer(payload))
  state = setup(data)
  _process_worker_settings!(blas_threads, gc_settings)
  _PROCESS_WORKER_STATE[key] = state
  c === nothing && (_PROCESS_WORKER_CACHE[] = _ProcessWorkerCache(digest, setup, state, true))
  return (time_ns() - t0) / 1.0e9
end

# The abort token of a batch in worker processes. A worker solves its item
# on its only thread and reads no message of the main process until the
# item ends, so the abort travels as a file: the main process creates it
# when its own token is set, and the solver loop on the worker checks for it
# once per Newton iteration through the armed token (sparlectra_check_abort),
# as a thread does with the token in memory. The workers run on this machine
# (addprocs without a host), so they see the file. One stat call per
# iteration costs nothing next to the factorization.
struct _FileAbortFlag
  path::String
end

Base.getindex(f::_FileAbortFlag) = isfile(f.path)

# One item on a worker: the pause brackets the item, so the GC of this
# process can run between two items. A non-empty `abort_path` arms the
# file token for the item.
function _process_worker_eval(key::Int, @nospecialize(eval_one), idx::Int, arg, abort_path::String)
  state = _PROCESS_WORKER_STATE[key]
  isempty(abort_path) && return with_gc_paused(() -> eval_one(state, idx, arg))
  sparlectra_arm_abort_token!(_FileAbortFlag(abort_path))
  try
    return with_gc_paused(() -> eval_one(state, idx, arg))
  finally
    sparlectra_disarm_abort_token!()
  end
end

# End a batch on a worker: the kept state is free for the next batch, a
# state that is not kept is dropped and collected.
function _process_worker_release!(key::Int)
  state = pop!(_PROCESS_WORKER_STATE, key, nothing)
  c = _PROCESS_WORKER_CACHE[]
  kept = c !== nothing && state !== nothing && c.state === state
  kept && (c.busy = false)
  isempty(_PROCESS_WORKER_STATE) && (_GC_PAUSE_OVERRIDE[] = nothing)
  kept || GC.gc()
  return nothing
end

# --- main process side ---------------------------------------------------------

_serialize_payload(x)::Vector{UInt8} = (io = IOBuffer(); Serialization.serialize(io, x); take!(io))

# The error a task of a @sync block ended with (the first one), so an abort
# or a worker error reaches the caller as itself and not wrapped.
function _first_task_error(err)
  err isa CompositeException || return err
  e = first(err.exceptions)
  return e isa TaskFailedException ? e.task.exception : e
end

"""
    _process_parallel_map(setup, payload, eval_one, item_arg, n, pids; blas_threads, progress) -> (results, stats)

The process counterpart of `_engine_parallel_map`: item `k` of `1:n` is
`eval_one(state, k, item_arg(k))` on a worker process, where `state =
setup(deserialize(payload))` was built once per worker; a worker keeps
the state of its last batch, and a later batch with the same payload (same
SHA-256) and setup reuses it without sending the payload. `setup`
and `eval_one` must be functions defined at the top level of Sparlectra
(they are sent by name). Items go to the next free worker (results do not
depend on the worker: every item starts from the same state). The abort
token is checked here before every item and watched while the items run:
when it is set, the workers stop their running items within one Newton
iteration (`_FileAbortFlag`) and the call throws `PowerFlowAborted`.
`progress` is counted here in the main process. The workers drop
the batch state in a `finally`. `stats` holds `transfer_s` (sending and
loading the payload on all workers, wall time), `setup_max_s` (the
slowest worker's load) and `workers_reused` (workers that kept the state).
"""
function _process_parallel_map(@nospecialize(setup), payload::Vector{UInt8}, @nospecialize(eval_one), @nospecialize(item_arg), n::Int, pids::Vector{Int}; blas_threads::Union{Nothing,Int} = 1, @nospecialize(progress = nothing))
  results = Vector{Any}(undef, n)
  key = Threads.atomic_add!(_BATCH_KEY, 1) + 1
  tracker = _BatchProgress(progress, n)
  gc_settings = _gc_pause_settings()
  # a worker's share of the guard limit: the workers run side by side
  gc_settings = (gc_settings[1], gc_settings[2] ÷ max(1, length(pids)))
  t0 = time_ns()
  # the payload digest names the data a worker kept from an earlier batch
  digest = sha256(payload)
  setup_s = zeros(length(pids))
  reused = falses(length(pids))
  stats = Dict{Symbol,Any}()
  # the abort token of the caller (a Web UI job), passed on to the workers
  # as a file that exists once the token is set; a watcher task of this
  # thread creates it
  token = _armed_abort_token()
  abort_path = token === nothing ? "" : tempname()
  watching = Ref(token !== nothing)
  watcher = token === nothing ? nothing : @async begin
    while watching[]
      if token[]
        touch(abort_path)
        break
      end
      sleep(0.05)
    end
  end
  try
    try
      @sync for (i, p) in enumerate(pids)
        @async begin
          reused[i] = _dcall(:remotecall_fetch, _process_worker_reuse!, p, key, setup, digest, blas_threads, gc_settings)
          reused[i] || (setup_s[i] = _dcall(:remotecall_fetch, _process_worker_setup!, p, key, setup, payload, digest, blas_threads, gc_settings))
        end
      end
      stats[:transfer_s] = (time_ns() - t0) / 1.0e9
      stats[:setup_max_s] = maximum(setup_s)
      stats[:workers_reused] = count(reused)
      # the tasks of this block all run on this thread (@async), so the
      # shared counter needs no lock
      next = Ref(0)
      failed = Ref(false)
      @sync for p in pids
        @async while true
          (failed[] || next[] >= n) && break
          next[] += 1
          local idx = next[]
          try
            sparlectra_check_abort()
            results[idx] = _dcall(:remotecall_fetch, _process_worker_eval, p, key, eval_one, idx, item_arg(idx), abort_path)
          catch
            # any failure (abort, worker error) stops the other workers
            # after their current item; the error itself goes on
            failed[] = true
            rethrow()
          end
          _progress_tick!(tracker)
        end
      end
    catch err
      # a failed worker task, an abort or a worker error ends the batch;
      # an abort reaches the caller as PowerFlowAborted (on a worker it ends
      # as a RemoteException), any other error as the first one itself
      token !== nothing && token[] && throw(PowerFlowAborted())
      throw(_first_task_error(err))
    end
  finally
    watching[] = false
    watcher === nothing || wait(watcher)
    @sync for p in pids
      @async try
        _dcall(:remotecall_wait, _process_worker_release!, p, key)
      catch err
        # a worker that died took its batch state with it
        err isa _dcall_type(:ProcessExitedException) || rethrow()
      end
    end
    isempty(abort_path) || rm(abort_path; force = true)
  end
  return results, stats
end

# runtime.parallel.max_processes resolved to a number; `auto` is the
# number of CPU threads, capped by the free memory for a net of `buses`
# buses (_process_memory_cap)
function _parallel_max_processes(cfg::ParallelRuntimeConfig, buses::Int)::Int
  v = lowercase(strip(cfg.max_processes))
  v == "auto" && return min(Sys.CPU_THREADS, _process_memory_cap(buses))
  return parse(Int, v)
end

# The workers 80 percent of the free memory holds. Peak resident memory per
# worker in the scaling run with the GC pause on: 3.1 GB on 9241 buses,
# 3.2 GB on 25k, 3.7 GB on 70k (the loaded package, the code compiled on the
# first items, the net twice and one solve), so 3 GB plus 10 kB per bus.
_process_memory_cap(buses::Int)::Int = max(1, floor(Int, 0.8 * Sys.free_memory() / (3.0e9 + 1.0e4 * buses)))

# The executor of a batch of `n` items on a network of `buses` buses:
# (:serial | :threads | :processes, number of tasks or processes). An
# explicit `executor` (keyword, else runtime.parallel.executor) wins.
# `auto` takes processes for networks of at least _PROCESS_BATCH_MIN_BUSES
# buses when the pool already holds the workers the batch uses, or when the
# batch is large enough to pay the pool start (items x buses at least
# _PROCESS_BATCH_COLD_WORK); threads otherwise. The rule is from the
# executor crossover on 1.3.0 with 24 workers (264 items, pool running):
# processes are 3.0 to 4.2 times faster from 1354 buses on, below 1000
# buses the order changes with the batch size; a batch that starts the pool
# (about 6 s) gains only from about 450 items on 2869 buses and 1100 to 3200
# items on 1354 to 2000 buses, i.e. from about 3e6 item-buses. When the free
# memory holds fewer than two workers, `auto` takes threads.
# `parallel_enabled = false` and batches below min_work_items run serially
# whatever the executor.
const _PROCESS_BATCH_MIN_BUSES = 1000
const _PROCESS_BATCH_COLD_WORK = 3.0e6

function _resolve_batch_executor(executor::Union{Nothing,Symbol}, n::Int, buses::Int; parallel_enabled::Union{Nothing,Bool} = nothing, parallel_max_tasks::Union{Nothing,Int} = nothing, parallel_min_work_items::Union{Nothing,Int} = nothing, max_processes::Union{Nothing,Int} = nothing)
  on, cap, min_items = _resolve_parallel_runtime(parallel_enabled, parallel_max_tasks, parallel_min_work_items)
  pcfg = runtime_config().parallel
  ex = executor === nothing ? Symbol(pcfg.executor) : executor
  ex in (:auto, :threads, :processes) || throw(ArgumentError("executor must be :auto, :threads or :processes, got :$(ex)"))
  (on && n >= min_items && n >= 2) || return (:serial, 1)
  auto = ex === :auto
  k = min(n, something(max_processes, _parallel_max_processes(pcfg, buses)))
  if auto
    # a running pool costs nothing more; a cold one must be paid by the batch
    warm = length(_BATCH_POOL) >= k
    ex = (buses >= _PROCESS_BATCH_MIN_BUSES && (warm || n * buses >= _PROCESS_BATCH_COLD_WORK)) ? :processes : :threads
  end
  if ex === :processes
    k >= 2 && return (:processes, k)
    auto || return (:serial, 1)
  end
  # the task cap as before: one chunk per task, also above nthreads()
  return (Threads.nthreads() > 1 && cap > 1) ? (:threads, cap) : (:serial, 1)
end

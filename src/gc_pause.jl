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

# file: src/gc_pause.jl
# purpose: the garbage collector pause of a run (runtime.gc_pause): the GC is
#          off while a power flow or a batch runs and on again when the last
#          run of the process has ended, also after an error. A resident
#          memory guard collects once when a run grows past its limit.

# How Julia switches the GC, and why the bookkeeping looks like this:
# GC.enable(false) sets a flag of the CALLING THREAD and counts the threads
# with the flag set; the GC runs only while that count is zero, and only the
# same thread can clear its flag again. A run therefore pauses the GC on its
# own thread and is pinned to that thread (sticky task) until it ends, so
# the switch back happens where it was made. Runs on several threads (a
# threaded batch, Web UI jobs) each hold the flag of their thread; the GC
# comes back when the last of them ends. `total` counts the active runs of
# the process, `depth` the runs per thread, `prev` the threads whose flag a
# run set, with the GC state before.
mutable struct _GcPauseState
  const lock::ReentrantLock
  const depth::Dict{Int,Int}
  const prev::Dict{Int,Bool}
  total::Int
  # the guard: limit of the resident memory in bytes (0: no guard), the
  # resident size at which the next collection is due, and whether the guard
  # has asked the threads to give their flag back
  limit_bytes::Int
  @atomic next_collect_bytes::Int
  @atomic tripped::Bool
  @atomic active::Bool
  @atomic peak_resident_bytes::Int
  guard_collections::Int
end

const _GC_PAUSE = _GcPauseState(ReentrantLock(), Dict{Int,Int}(), Dict{Int,Bool}(), 0, 0, 0, false, false, 0, 0)

# Process-wide settings that replace runtime.gc_pause, set by a worker
# process of the process batch for the duration of a batch (a worker does
# not have the configuration of the main process). `nothing`: use the
# active configuration.
const _GC_PAUSE_OVERRIDE = Ref{Union{Nothing,Tuple{Bool,Int}}}(nothing)

# The resident set size of this process in bytes (libuv, on every platform
# Julia runs on). 0 when libuv cannot read it, which switches the guard off.
function _resident_bytes()::Int
  r = Ref{Csize_t}(0)
  rc = ccall(:uv_resident_set_memory, Cint, (Ref{Csize_t},), r)
  return rc == 0 ? Int(r[]) : 0
end

# runtime.gc_pause resolved to (enabled, limit in bytes): `auto` is half of
# the physical memory, 0 GB switches the guard off.
function _gc_pause_settings()::Tuple{Bool,Int}
  ov = _GC_PAUSE_OVERRIDE[]
  ov === nothing || return ov
  cfg = runtime_config().gc_pause
  return (cfg.enabled, _gc_pause_limit_bytes(cfg.max_resident_gb))
end

function _gc_pause_limit_bytes(max_resident_gb::AbstractString)::Int
  v = lowercase(strip(max_resident_gb))
  v == "auto" && return Int(Sys.total_memory() ÷ 2)
  return round(Int, parse(Float64, v) * 2^30)
end

"""
    _gc_pause_enter!(; enabled, limit_bytes) -> token

Start a run with the GC paused (see `src/gc_pause.jl`): pins the calling
task to its thread, counts the run and switches the GC off on this thread
if no run of the thread has done so. The first run of the process sets the
memory limit of the guard. `enabled = false` returns `nothing` and changes
nothing. Every token must be passed to `_gc_pause_exit!` in a
`finally` block of the same task.
"""
function _gc_pause_enter!(; enabled::Union{Nothing,Bool} = nothing, limit_bytes::Union{Nothing,Int} = nothing)
  on, limit = _gc_pause_settings()
  something(enabled, on) || return nothing
  ct = current_task()
  sticky = ct.sticky
  ct.sticky = true
  tid = Threads.threadid()
  st = _GC_PAUSE
  lock(st.lock) do
    st.depth[tid] = get(st.depth, tid, 0) + 1
    st.total += 1
    if st.total == 1
      st.limit_bytes = something(limit_bytes, limit)
      @atomic st.next_collect_bytes = st.limit_bytes
      @atomic st.tripped = false
      @atomic st.active = true
    end
    if !(@atomic st.tripped) && !haskey(st.prev, tid)
      st.prev[tid] = GC.enable(false)
    end
  end
  return (tid = tid, sticky = sticky)
end

"""
    _gc_pause_exit!(token) -> Nothing

End a run started by `_gc_pause_enter!`: the thread gives its GC
flag back when its last run ends, the GC is on again when no run of the
process is left, and the task is no longer pinned (unless it was before).
"""
function _gc_pause_exit!(token)
  token === nothing && return nothing
  st = _GC_PAUSE
  tid = token.tid
  lock(st.lock) do
    # the task is pinned, so this is the thread that set the flag; only
    # code that unpins the task inside a run gets here on another thread
    same_thread = Threads.threadid() == tid
    same_thread || @warn "GC pause: the run ended on thread $(Threads.threadid()), it started on thread $(tid); the GC flag of thread $(tid) stays set"
    d = st.depth[tid] - 1
    d == 0 ? delete!(st.depth, tid) : (st.depth[tid] = d)
    st.total -= 1
    same_thread && d == 0 && haskey(st.prev, tid) && GC.enable(pop!(st.prev, tid))
    if st.total == 0
      @atomic st.active = false
      @atomic st.tripped = false
    elseif (@atomic st.tripped) && isempty(st.prev)
      _gc_collect_once!(st)
    end
  end
  current_task().sticky = token.sticky
  return nothing
end

# The guard has fired and no thread holds a flag any more: one full
# collection, the next one due a quarter of the limit above what is left,
# and the runs of this thread pause again. Called under st.lock.
function _gc_collect_once!(st::_GcPauseState)
  GC.gc()
  st.guard_collections += 1
  rss = _resident_bytes()
  @atomic st.next_collect_bytes = max(st.limit_bytes, rss + st.limit_bytes ÷ 4)
  @atomic st.tripped = false
  tid = Threads.threadid()
  get(st.depth, tid, 0) > 0 && !haskey(st.prev, tid) && (st.prev[tid] = GC.enable(false))
  return nothing
end

"""
    _gc_pause_guard!() -> Nothing

The memory guard, called once per Newton iteration and per batch item:
nothing while no run pauses the GC or the resident memory is below the
limit. Above it the thread gives its GC flag back; the thread that gives
back the last flag collects once, and the runs pause again from the next
item on. With one run in the process (a worker process, a single power
flow) that is one collection right here.
"""
function _gc_pause_guard!()
  st = _GC_PAUSE
  (@atomic st.active) || return nothing
  due = @atomic st.next_collect_bytes
  due > 0 || return nothing
  rss = _resident_bytes()
  peak = @atomic st.peak_resident_bytes
  rss > peak && (@atomic st.peak_resident_bytes = rss)
  rss > due || return nothing
  lock(st.lock) do
    (@atomic st.active) || return
    rss > (@atomic st.next_collect_bytes) || return
    @atomic st.tripped = true
    tid = Threads.threadid()
    haskey(st.prev, tid) && GC.enable(pop!(st.prev, tid))
    isempty(st.prev) && _gc_collect_once!(st)
  end
  return nothing
end

"""
    gc_pause_status() -> NamedTuple

The state of the GC pause of this process: `active_runs` (runs that pause
the GC now), `paused` (whether a run holds the GC off now), `limit_bytes` (the
guard limit of the current or last pause, 0 without guard),
`guard_collections` (collections of the guard since the process started)
and `peak_resident_bytes` (the highest resident size a guard check saw).
"""
function gc_pause_status()
  st = _GC_PAUSE
  return lock(st.lock) do
    (active_runs = st.total, paused = !isempty(st.prev), limit_bytes = st.limit_bytes, guard_collections = st.guard_collections, peak_resident_bytes = @atomic(st.peak_resident_bytes))
  end
end

"""
    with_gc_paused(f; enabled = nothing, limit_bytes = nothing)

Run `f()` with the GC paused as a run of the process (see
`runtime.gc_pause` in [Parallel Execution](@ref parallel_execution_gc_pause)):
off when the first run starts, on again when the last one ends, also when
`f` throws. `enabled` and `limit_bytes` replace the configured values.
Returns what `f` returns.
"""
function with_gc_paused(f; enabled::Union{Nothing,Bool} = nothing, limit_bytes::Union{Nothing,Int} = nothing)
  token = _gc_pause_enter!(; enabled = enabled, limit_bytes = limit_bytes)
  try
    return f()
  finally
    _gc_pause_exit!(token)
  end
end

# Wrap the body of a function definition in a GC pause:
#   @_gc_paused_run function runpf!(...) ... end
# The docstring above the macro call documents the expanded definition.
macro _gc_paused_run(def)
  (def isa Expr && def.head in (:function, :(=)) && length(def.args) == 2) || error("@_gc_paused_run expects a function definition")
  body = def.args[2]
  def.args[2] = quote
    local _gc_pause_token = _gc_pause_enter!()
    try
      $(body)
    finally
      _gc_pause_exit!(_gc_pause_token)
    end
  end
  return esc(def)
end

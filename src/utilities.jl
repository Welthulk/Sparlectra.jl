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

# file: src/utilities.jl
# purpose: performance-profile timing helpers
#          (_perf_profile_time! and friends)
function _perf_profile_enabled(profile)::Bool
  profile isa AbstractDict || return false
  return Bool(get(profile, :enabled, false))
end

function _perf_profile_wants_allocations(profile)::Bool
  _perf_profile_enabled(profile) || return false
  return Bool(get(profile, :show_allocations, false))
end

function _perf_profile_add!(profile, phase::Symbol, elapsed_s::Real, bytes::Integer = 0)
  _perf_profile_enabled(profile) || return nothing
  timings = get!(profile, :timings) do
    Dict{Symbol,NamedTuple{(:calls,:elapsed_s,:bytes),Tuple{Int,Float64,Int}}}()
  end
  prev = get(timings, phase, (calls = 0, elapsed_s = 0.0, bytes = 0))
  timings[phase] = (calls = prev.calls + 1, elapsed_s = prev.elapsed_s + Float64(elapsed_s), bytes = prev.bytes + Int(bytes))
  return nothing
end

function _perf_profile_push_iteration!(profile, row)
  _perf_profile_enabled(profile) || return nothing
  push!(get!(profile, :iterations, NamedTuple[]), row)
  return nothing
end

function _perf_profile_time!(f::F, profile, phase::Symbol) where {F}
  if !_perf_profile_enabled(profile)
    return f()
  elseif _perf_profile_wants_allocations(profile)
    t0 = time_ns()
    completed = false
    timed = nothing
    try
      timed = @timed f()
      completed = true
      _perf_profile_add!(profile, phase, timed.time, timed.bytes)
      return timed.value
    finally
      completed || _perf_profile_add!(profile, phase, (time_ns() - t0) / 1e9, 0)
    end
  else
    t0 = time_ns()
    completed = false
    try
      value = f()
      completed = true
      _perf_profile_add!(profile, phase, (time_ns() - t0) / 1e9, 0)
      return value
    finally
      completed || _perf_profile_add!(profile, phase, (time_ns() - t0) / 1e9, 0)
    end
  end
end

# Keys a parallel worker's child profile inherits from the parent by value.
# Deliberately WITHOUT :phase_callback (mutates a shared recorder; the
# orchestrator alone invokes it, per the Phase 0 review) and WITHOUT
# :diagnostic_artifact_prefix (per-worker value, passed explicitly).
# :cancellation_check stays: it closes over a Threads.Atomic{Bool}, safe to
# poll from any task.
const _PERF_PROFILE_CHILD_KEYS = (:enabled, :show_allocations, :capture_initial_residual_rows, :output_dir, :cancellation_check)

"""
    _perf_profile_child(profile) -> Dict{Symbol,Any} | nothing

Create a fresh per-worker profile for a parallel work item (thread-safety
Phase 1). Returns `nothing` when `profile` is not an enabled profile Dict;
otherwise a new Dict carrying only the read-only seed keys
(`_PERF_PROFILE_CHILD_KEYS`). Workers write timings, iteration rows, and
scalar diagnostics into their child only; the orchestrator folds children
back with [`_perf_profile_merge!`](@ref) after `fetch`, so the parent Dict
is never touched concurrently.
"""
function _perf_profile_child(profile)
  _perf_profile_enabled(profile) || return nothing
  child = Dict{Symbol,Any}()
  for key in _PERF_PROFILE_CHILD_KEYS
    haskey(profile, key) && (child[key] = profile[key])
  end
  return child
end

"""
    _perf_profile_merge!(parent, child, prefix::AbstractString)

Fold a worker's child profile into the parent (serial, orchestrator-only).
`:timings` rows are SUMMED into the parent's `:timings` under their original
phase names, so the phase-name set of a parallel run stays identical to the
serial run (calls/elapsed/bytes accumulate exactly like the serial loop
did). `:iterations` rows are appended in call order. Every other child key
is copied under `Symbol(prefix, key)` so per-worker scalars (backend names,
matrix stats, start-projection summaries) stop overwriting each other.
The wall-clock of the orchestrating fan-out is accounted separately by the
caller via `_perf_profile_add!(parent, :parallel_wall_time, elapsed, 0)`,
which keeps serial-sum versus wall-clock visible in performance.log.
"""
function _perf_profile_merge!(parent, child, prefix::AbstractString)
  _perf_profile_enabled(parent) || return parent
  child isa AbstractDict || return parent
  for (key, value) in child
    key in _PERF_PROFILE_CHILD_KEYS && continue
    if key === :timings && value isa AbstractDict
      timings = get!(parent, :timings) do
        Dict{Symbol,Any}()
      end
      for (phase, row) in value
        prev = get(timings, phase, (calls = 0, elapsed_s = 0.0, bytes = 0))
        timings[phase] = (calls = prev.calls + row.calls, elapsed_s = prev.elapsed_s + Float64(row.elapsed_s), bytes = prev.bytes + Int(row.bytes))
      end
    elseif key === :iterations && value isa AbstractVector
      append!(get!(parent, :iterations, NamedTuple[]), value)
    else
      parent[Symbol(prefix, key)] = value
    end
  end
  return parent
end

function _solver_elapsed_from_profile(profile)
  timings = profile isa AbstractDict ? get(profile, :timings, nothing) : nothing
  row = timings isa AbstractDict ? get(timings, :solver_total, nothing) : nothing
  if row isa NamedTuple && hasproperty(row, :elapsed_s)
    elapsed = Float64(row.elapsed_s)
    return isfinite(elapsed) ? max(0.0, elapsed) : nothing
  end
  legacy_row = profile isa AbstractDict ? get(profile, :solver_total, nothing) : nothing
  elapsed =
    legacy_row isa NamedTuple && hasproperty(legacy_row, :elapsed_s) ? Float64(legacy_row.elapsed_s) :
    legacy_row isa Number ? Float64(legacy_row) : nothing
  return elapsed !== nothing && isfinite(elapsed) ? max(0.0, elapsed) : nothing
end

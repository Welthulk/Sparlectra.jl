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

# file: src/qsts/qsts_metrics.jl
# purpose: metrics of a sequential QSTS run (qsts.mode = sequential): the
#          hunting pattern of a controller as a pure function, violation
#          episodes of the voltage band and the loading limit, loss energy,
#          the metric rows per run, controller, bus and branch, and the
#          column tables qstsSwitching and qstsMetrics.

# Hunting patterns in the net moves of one controller: `moves` holds
# (step, sign) per counted move in step order, sign +1 or -1; an entry with
# sign 0 is a reset (a resync step) that no pattern spans. Three consecutive
# moves with alternating sign whose first and last step lie at most k - 1
# apart are one pattern ("up, down, up within k steps"); overlapping triples
# count each, so + - + - gives two. Throws for k < 3, which can hold no
# three moves. The sequential run counts incrementally with the same rule
# (_qsts_hunting_step); the two must agree.
function _qsts_hunting(moves::AbstractVector{Tuple{Int,Int}}, k::Int)::Int
  k >= 3 || throw(ArgumentError(string("qsts.hunting_window must be >= 3 (a hunting pattern is three moves); got ", k, ".")))
  n = 0
  for i in 3:length(moves)
    (s1, a), (_, b), (s3, c) = moves[i - 2], moves[i - 1], moves[i]
    (a != 0 && b != 0 && c != 0) || continue
    (b == -a && c == -b && s3 - s1 <= k - 1) && (n += 1)
  end
  return n
end

# The incremental form of _qsts_hunting for one new counted move (step,
# sign) after the last two moves `last` ((0, 0) is no move): whether the
# three form a pattern, and the new pair of last moves.
function _qsts_hunting_step(last::NTuple{2,Tuple{Int,Int}}, step::Int, sign::Int, k::Int)::Tuple{Bool,NTuple{2,Tuple{Int,Int}}}
  (s1, a), (_, b) = last
  hit = a != 0 && b != 0 && b == -a && sign == -b && step - s1 <= k - 1
  return hit, (last[2], (step, sign))
end

# Minutes per step: the fixed step length, else the distance to the next
# timestamp (the last step takes the one before it); nothing without a time
# reference (no timestamps, an unreadable one, or a single step without a
# step length).
function _qsts_step_durations(step_minutes::Union{Nothing,Float64}, timestamps::Vector{String})::Union{Nothing,Vector{Float64}}
  n = length(timestamps)
  step_minutes === nothing || return fill(step_minutes, n)
  n >= 2 || return nothing
  times = DateTime[]
  for ts in timestamps
    t = _qsts_parse_time(ts)
    t === nothing && return nothing
    push!(times, t)
  end
  d = Float64[Dates.value(times[k + 1] - times[k]) / 60000.0 for k in 1:(n - 1)]
  push!(d, d[end])
  return d
end

# an ISO 8601 UTC time (YYYY-MM-DDThh:mm[:ss[.s]]Z), nothing for anything else
function _qsts_parse_time(ts::AbstractString)::Union{Nothing,DateTime}
  endswith(ts, "Z") || return nothing
  body = ts[1:(end - 1)]
  for f in (dateformat"yyyy-mm-ddTHH:MM:SS.s", dateformat"yyyy-mm-ddTHH:MM:SS", dateformat"yyyy-mm-ddTHH:MM")
    t = tryparse(DateTime, body, f)
    t === nothing || return t
  end
  return nothing
end

# Episodes of one element over the steps: `excess[k]` is the amount beyond
# the limit (> 0 in violation), NaN where the step has no value, which ends
# an episode like a step inside the limit. Returns (episodes, steps,
# minutes or NaN without durations, largest excess).
function _qsts_episodes(excess::AbstractVector{Float64}, durations::Union{Nothing,Vector{Float64}})::Tuple{Int,Int,Float64,Float64}
  episodes = 0
  steps = 0
  minutes = 0.0
  worst = 0.0
  inside = false
  for k in eachindex(excess)
    x = excess[k]
    if !isnan(x) && x > 0.0
      inside || (episodes += 1)
      inside = true
      steps += 1
      durations === nothing || (minutes += durations[k])
      worst = max(worst, x)
    else
      inside = false
    end
  end
  return episodes, steps, durations === nothing ? NaN : minutes, worst
end

# The metric rows of a sequential run, in this order: the run (element_id
# empty), each controller, each reported bus with a violation of the band,
# each branch with an overload. A row with NaN (minutes and energy without
# a time reference) is left out.
function _qsts_metric_rows(e::QstsElements, steps::Vector{QstsStepInfo}, resync::BitVector, bus_values::Array{Float64,3}, branch_values::Array{Float64,3}, step_values::Matrix{Float64},
                           controller_values::Array{Float64,3}, switch_count::Vector{Int}, n_hunting::Vector{Int}, q::QstsConfig,
                           durations::Union{Nothing,Vector{Float64}})::Vector{_QstsMetricRow}
  rows = _QstsMetricRow[]
  add(id, metric, value) = isnan(value) || push!(rows, (element_id = id, metric = metric, value = Float64(value)))
  n = length(steps)
  solved = [s.status === :solved for s in steps]
  loss = 0.0
  for k in 1:n
    (solved[k] && durations !== nothing) || continue
    loss += step_values[_QS_P_LOSS, k] * durations[k] / 60.0
  end
  add("", "n_steps_solved", count(solved))
  add("", "n_failed", n - count(solved))
  add("", "switch_count_total", sum(switch_count; init = 0))
  add("", "n_hunting_total", sum(n_hunting; init = 0))
  add("", "n_resync", count(resync))
  add("", "loss_energy_mwh", durations === nothing ? NaN : loss)
  for c in eachindex(e.controller_ids)
    id = e.controller_ids[c]
    d = filter(!isnan, controller_values[c, _QC_DPOS_REF, :])
    add(id, "switch_count", switch_count[c])
    add(id, "n_hunting", n_hunting[c])
    add(id, "steps_at_limit", count(==(1.0), controller_values[c, _QC_AT_LIMIT, :]))
    add(id, "steps_locked", count(==(1.0), controller_values[c, _QC_LOCKED, :]))
    add(id, "max_abs_dpos_ref", isempty(d) ? NaN : maximum(abs, d))
    add(id, "rms_dpos_ref", isempty(d) ? NaN : sqrt(sum(abs2, d) / length(d)))
  end
  episode_rows(id, prefix, unit, excess) = begin
    episodes, nsteps, minutes, worst = _qsts_episodes(excess, durations)
    nsteps == 0 && return nothing
    add(id, string(prefix, "_episodes"), episodes)
    add(id, string(prefix, "_steps"), nsteps)
    add(id, string(prefix, "_minutes"), minutes)
    add(id, string(prefix, "_max_", unit), worst)
  end
  for i in eachindex(e.bus_ids)
    e.bus_reported[i] || continue
    vm = bus_values[i, _QB_VM, :]
    episode_rows(e.bus_ids[i], "vm_low", "pu", Float64[isnan(v) ? NaN : q.vm_min_pu - v for v in vm])
    episode_rows(e.bus_ids[i], "vm_high", "pu", Float64[isnan(v) ? NaN : v - q.vm_max_pu for v in vm])
  end
  for j in eachindex(e.branch_ids)
    episode_rows(e.branch_ids[j], "overload", "pct", Float64[isnan(v) ? NaN : v - q.loading_max_pct for v in branch_values[j, _QR_LOADING, :]])
  end
  return rows
end

"""
    qstsSwitching(r::QstsResult) -> (step, element_id, old, new, steps, trigger)

The tap and bank moves of a sequential run as a named tuple of equally long
column vectors, one row per step and controller whose position changed
(the net move of the step: a search back and forth inside one controller
run is not listed, it shows in `outer_iterations`). `old` is the position
at the end of the step before (the carried one), `new` the position at the
end of the step, `steps = new - old`, positions as in
[`QSTS_CONTROLLER_QUANTITIES`](@ref). `trigger` is `"control"` (counted in
the switch count), `"range"` (the position stood outside the control range
before the solve) or `"resync"` (a resync step, which starts from the
position of the series). Empty for a snapshot run; written as
`qsts_switching.csv` by [`writeQstsResults`](@ref).
"""
function qstsSwitching(r::QstsResult)
  w = r.switching
  return (step = Int[x.step for x in w], element_id = String[x.element_id for x in w], old = Int[x.old for x in w], new = Int[x.new for x in w],
          steps = Int[x.steps for x in w], trigger = String[String(x.trigger) for x in w])
end

"""
    qstsMetrics(r::QstsResult) -> (element_id, metric, value)

The metrics of a sequential run as a named tuple of three column vectors:

- the run (`element_id` empty): `n_steps_solved`, `n_failed`,
  `switch_count_total`, `n_hunting_total`, `n_resync`, `loss_energy_mwh`
  (sum of `total_p_loss_mw` times the step length over the solved steps);
- each controller: `switch_count` (sum of the counted steps),
  `n_hunting` (patterns of three counted moves with alternating sign within
  `qsts.hunting_window` steps; a resync step breaks a pattern),
  `steps_at_limit`, `steps_locked`, `max_abs_dpos_ref`, `rms_dpos_ref`;
- each bus outside `qsts.vm_min_pu` or `qsts.vm_max_pu` and each branch above
  `qsts.loading_max_pct` in at least one step: `<kind>_episodes`
  (consecutive steps in violation; a step without a value ends one),
  `<kind>_steps`, `<kind>_minutes` and `<kind>_max_pu` or `_max_pct` (the
  largest excess beyond the limit), with `<kind>` `vm_low`, `vm_high` or
  `overload`.

The step length is `step_minutes` of the series, else the distance to the
next timestamp; without either the `_minutes` rows and `loss_energy_mwh`
are left out. Empty for a snapshot run; written as `qsts_metrics.csv` by
[`writeQstsResults`](@ref).
"""
function qstsMetrics(r::QstsResult)
  m = r.metrics
  return (element_id = String[x.element_id for x in m], metric = String[x.metric for x in m], value = Float64[x.value for x in m])
end

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
#
# Example: global and local observability / redundancy for state estimation
#
# This script progressively deactivates branch-flow measurements and logs
# global/local redundancy metrics and warnings when local observability is lost.

# file: src/examples/state_estimation_observability.jl

using Sparlectra
using Printf
using Dates
using Random

include("state_estimation_wls.jl")

const OUTDIR_OBS = joinpath(@__DIR__, "_out")
const CASEFILE_OBS = "case9.m"

"""
Return the state-vector column indices for voltage angles and magnitudes,
and identify the slack bus index used to exclude the slack angle state.
"""
function _state_column_groups(net::Sparlectra.Net)
  nbus = length(net.nodeVec)
  # Find the slack bus, because its angle is not part of the estimated state vector.
  slack_idx = 1
  for (i, node) in enumerate(net.nodeVec)
    if getNodeType(node) == :Slack
      slack_idx = i
      break
    end
  end

  angle_cols = Int[]
  p = 0
  # Angle states contain all buses except slack, in the internal estimator order.
  for i in axes(net.nodeVec, 1)
    if i != slack_idx
      p += 1
      push!(angle_cols, p)
    end
  end

  # Voltage-magnitude states are stored in the second block of the state vector.
  vm_cols = collect((nbus):(2*nbus-1))
  return angle_cols, vm_cols, slack_idx
end

"""
Compute the branch degree (number of incident branches) for each bus.
"""
function _bus_degrees(net::Sparlectra.Net)
  degrees = fill(0, length(net.nodeVec))
  for br in net.branchVec
    degrees[br.fromBus] += 1
    degrees[br.toBus] += 1
  end
  return degrees
end

"""
Select one or more random non-slack target buses, preferring non-boundary
buses (degree > 1) to build a meaningful local observability stress test.
"""
function _random_target_buses(net::Sparlectra.Net)
  degrees = _bus_degrees(net)
  candidates = Int[]
  # Prefer non-slack buses with at least two incident branches (non-boundary buses).
  for i in eachindex(net.nodeVec)
    if getNodeType(net.nodeVec[i]) != :Slack && degrees[i] > 1
      push!(candidates, i)
    end
  end

  if isempty(candidates)
    for i in eachindex(net.nodeVec)
      getNodeType(net.nodeVec[i]) != :Slack && push!(candidates, i)
    end
  end

  # Fallback for very small/special cases where only boundary buses exist.
  isempty(candidates) && return [1]
  shuffle!(candidates)
  # Randomly choose one to three target buses to stress a local area, not only one node.
  max_targets = min(3, length(candidates))
  n_targets = rand(1:max_targets)
  return sort!(candidates[1:n_targets])
end

"""
Create a readable label for a measurement used in step-by-step logs.
"""
function _measurement_label(m::Sparlectra.Measurement)
  if !isnothing(m.branchIdx)
    return "branch=$(m.branchIdx),dir=$(m.direction),typ=$(m.typ)"
  elseif !isnothing(m.busIdx)
    return "bus=$(m.busIdx),typ=$(m.typ)"
  end
  return "typ=$(m.typ)"
end

"""
Parse a measurement label string into table-ready columns.
"""
function _measurement_table_row(label::String, state::String)
  parts = split(label, ",")
  data = Dict{String,String}()
  for part in parts
    kv = split(part, "="; limit = 2)
    length(kv) == 2 || continue
    data[kv[1]] = kv[2]
  end

  if haskey(data, "branch")
    return [state, "branch", data["branch"], get(data, "dir", "-"), get(data, "typ", "-")]
  elseif haskey(data, "bus")
    return [state, "bus", data["bus"], "-", get(data, "typ", "-")]
  end
  return [state, "unknown", "-", "-", get(data, "typ", label)]
end

"""
Print deactivated measurements in a compact ASCII table.
"""
function _print_deactivated_table(io::IO, title::String, labels::Vector{String}, state::String)
  @printf(io, "  %s:\n", title)
  if isempty(labels)
    @printf(io, "    (none)\n")
    return nothing
  end

  rows = [_measurement_table_row(label, state) for label in labels]
  headers = ["State", "Scope", "ID", "Dir", "Type"]
  widths = [length(h) for h in headers]
  for row in rows
    for i in eachindex(row)
      widths[i] = max(widths[i], length(row[i]))
    end
  end

  function print_row(vals::Vector{String})
    @printf(io, "    |")
    for i in eachindex(vals)
      @printf(io, " %-*s |", widths[i], vals[i])
    end
    @printf(io, "\n")
  end

  sep = "    +" * join([repeat("-", w + 2) for w in widths], "+") * "+"
  @printf(io, "%s\n", sep)
  print_row(headers)
  @printf(io, "%s\n", sep)
  for row in rows
    print_row(row)
  end
  @printf(io, "%s\n", sep)
  return nothing
end

"""
Map a vector of measurement indices to their corresponding log labels.
"""
function _critical_labels(meas::AbstractVector, idx::Vector{Int})
  labels = String[]
  for i in idx
    if 1 <= i <= length(meas)
      push!(labels, _measurement_label(meas[i]))
    end
  end
  return labels
end

"""
Set the active flag of one measurement by reconstructing the immutable struct.
"""
function _set_measurement_active!(meas::AbstractVector, idx::Int, active::Bool)
  m = meas[idx]
  meas[idx] = Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = active, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id)
  return nothing
end

"""
Collect active flow measurements for one branch, grouped by direction
(:from and :to), with deterministic ordering for reproducible logs.
"""
function _branch_direction_measurement_groups(meas::AbstractVector, branch_idx::Int)
  by_direction = Dict{Symbol,Vector{Int}}()
  for i in eachindex(meas)
    m = meas[i]
    is_flow = (m.typ == Sparlectra.PflowMeas || m.typ == Sparlectra.QflowMeas)
    if m.active && is_flow && !isnothing(m.branchIdx) && m.branchIdx == branch_idx && m.direction != :none
      dir = m.direction
      if !haskey(by_direction, dir)
        by_direction[dir] = Int[]
      end
      push!(by_direction[dir], i)
    end
  end

  groups = Vector{Vector{Int}}()
  for dir in sort!(collect(keys(by_direction)))
    idx = by_direction[dir]
    sort!(idx; by = i -> Int(meas[i].typ))
    push!(groups, idx)
  end
  return groups
end

"""
Collect active injection measurements (P/Q) for a specific bus.
"""
function _bus_injection_measurement_group(meas::AbstractVector, bus_idx::Int)
  idx = Int[]
  for i in eachindex(meas)
    m = meas[i]
    is_inj = (m.typ == Sparlectra.PinjMeas || m.typ == Sparlectra.QinjMeas)
    if m.active && is_inj && !isnothing(m.busIdx) && m.busIdx == bus_idx
      push!(idx, i)
    end
  end
  sort!(idx; by = i -> Int(meas[i].typ))
  return idx
end

"""
Return branch IDs adjacent to a target bus and the neighboring bus IDs.
"""
function _target_neighborhood(net::Sparlectra.Net, target_bus::Int)
  adjacent_branches = Int[]
  neighboring_buses = Int[]
  seen_buses = Set{Int}()

  for br in net.branchVec
    if br.fromBus == target_bus
      push!(adjacent_branches, br.branchIdx)
      if !(br.toBus in seen_buses)
        push!(neighboring_buses, br.toBus)
        push!(seen_buses, br.toBus)
      end
    elseif br.toBus == target_bus
      push!(adjacent_branches, br.branchIdx)
      if !(br.fromBus in seen_buses)
        push!(neighboring_buses, br.fromBus)
        push!(seen_buses, br.fromBus)
      end
    end
  end

  sort!(adjacent_branches)
  sort!(neighboring_buses)
  return adjacent_branches, neighboring_buses
end

"""
Append deactivation groups for the provided branch IDs to the global group list.
"""
function _push_groups_for_branches!(groups::Vector{Tuple{String,Union{Nothing,Int},Vector{Int}}}, meas::AbstractVector, branch_ids::Vector{Int})
  for br_id in branch_ids
    for g in _branch_direction_measurement_groups(meas, br_id)
      push!(groups, ("branch", br_id, g))
    end
  end
  return groups
end

"""
Append deactivation groups for injection measurements at the provided buses.
"""
function _push_groups_for_buses!(groups::Vector{Tuple{String,Union{Nothing,Int},Vector{Int}}}, meas::AbstractVector, bus_ids::Vector{Int})
  for bus_id in bus_ids
    g = _bus_injection_measurement_group(meas, bus_id)
    isempty(g) || push!(groups, ("injection", nothing, g))
  end
  return groups
end

"""
Build the ordered measurement-deactivation plan:
1) focused branch flows near targets,
2) focused injections near targets,
3) remaining branch flows,
4) remaining injections.
"""
function _deactivation_groups(net::Sparlectra.Net, meas::AbstractVector, target_buses::Vector{Int})
  groups = Tuple{String,Union{Nothing,Int},Vector{Int}}[]
  focused_adjacent = Set{Int}()
  focused_neighbors = Set{Int}()

  for target_bus in target_buses
    adjacent_branches, neighboring_buses = _target_neighborhood(net, target_bus)
    # Deactivation order: first branch-flow measurements around target area...
    _push_groups_for_branches!(groups, meas, adjacent_branches)
    # ...then injection measurements on target and neighboring buses.
    _push_groups_for_buses!(groups, meas, [target_bus])
    _push_groups_for_buses!(groups, meas, neighboring_buses)
    union!(focused_adjacent, adjacent_branches)
    union!(focused_neighbors, neighboring_buses)
  end

  active_branch_ids = _active_branch_ids(meas)
  # After focused groups, continue with all remaining branch-flow measurements.
  remaining_branches = [id for id in active_branch_ids if !(id in focused_adjacent)]
  _push_groups_for_branches!(groups, meas, remaining_branches)

  active_injection_buses = Int[]
  seen_inj = Set{Int}()
  for m in meas
    is_inj = (m.typ == Sparlectra.PinjMeas || m.typ == Sparlectra.QinjMeas)
    if m.active && is_inj && !isnothing(m.busIdx)
      bus = something(m.busIdx, 0)
      if bus > 0 && !(bus in seen_inj)
        push!(active_injection_buses, bus)
        push!(seen_inj, bus)
      end
    end
  end
  # Finally deactivate remaining injections outside the focused neighborhood.
  remaining_injection_buses = [b for b in sort!(active_injection_buses) if !(b in target_buses) && !(b in focused_neighbors)]
  _push_groups_for_buses!(groups, meas, remaining_injection_buses)

  return groups, sort!(collect(focused_adjacent)), sort!(collect(focused_neighbors))
end

"""
Return sorted branch IDs that still have at least one active flow measurement.
"""
function _active_branch_ids(meas::AbstractVector)
  ids = Set{Int}()
  for m in meas
    if m.active && (m.typ == Sparlectra.PflowMeas || m.typ == Sparlectra.QflowMeas) && !isnothing(m.branchIdx)
      push!(ids, something(m.branchIdx, 0))
    end
  end
  return sort!(collect(ids))
end

"""
Reset non-slack bus states to a neutral point before each WLS execution.
"""
function _reset_non_slack!(net::Sparlectra.Net)
  for n in net.nodeVec
    if getNodeType(n) != :Slack
      n._vm_pu = 1.0
      n._va_deg = 0.0
    end
  end
end

"""
Run the observability reduction experiment and stream a detailed step log
to the provided IO handle.
"""
function run_state_estimation_observability_example(io::IO)
  mkpath(OUTDIR_OBS)

  local_case = joinpath(Sparlectra.MPOWER_DIR, CASEFILE_OBS)
  case_path = if isfile(local_case)
    local_case
  else
    Sparlectra.FetchMatpowerCase.ensure_casefile(CASEFILE_OBS; outdir = Sparlectra.MPOWER_DIR, to_jl = false, overwrite = false)
  end

  net = createNetFromMatPowerFile(filename = case_path)
  _, erg_pf = runpf!(net, 40, 1e-10, 0; method = :rectangular, opt_sparse = true)
  erg_pf == 0 || error("Power flow did not converge")

  std = measurementStdDevs(vm = 0.1, pinj = 1.5, qinj = 1.5, pflow = 0.7, qflow = 0.9)
  setMeasurementsFromPF!(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = true, stddev = std)
  meas = net.measurements

  angle_cols, vm_cols, slack_idx = _state_column_groups(net)
  # Fixed seed keeps the random bus choice reproducible between runs/log comparisons.
  Random.seed!(42)
  target_buses = _random_target_buses(net)
  local_bus = target_buses[1]
  # Local observability is tracked on one representative target bus.
  local_cols = [angle_cols[1], vm_cols[local_bus]]

  all_branch_ids = _active_branch_ids(meas)
  @printf(io, "Observability reduction example for %s\n", CASEFILE_OBS)
  @printf(io, "Slack bus index: %d\n", slack_idx)
  @printf(io, "Selected target buses (random, non-boundary preferred): %s\n", string(target_buses))
  @printf(io, "Local bus under test: %d\n", local_bus)
  @printf(io, "Local columns under test: %s\n", string(local_cols))
  @printf(io, "Initial active branch IDs: %s\n\n", string(all_branch_ids))

  function log_step(step::Int, removed_branch::Union{Nothing,Int}, removed_measurements::Vector{String}, all_deactivated_measurements::Vector{String})
    # Recompute both global and local observability after every deactivation step.
    global_obs = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
    local_obs = evaluate_local_observability(net, local_cols; flatstart = true, jacEps = 1e-6)

    active_total = count(m -> m.active, meas)
    active_branch = length(_active_branch_ids(meas))
    removed_label = isnothing(removed_branch) ? "none" : string(removed_branch)
    @printf(io, "Step %d | removed branch: %s\n", step, removed_label)
    _print_deactivated_table(io, "Newly deactivated measurements", removed_measurements, "new")
    _print_deactivated_table(io, "All deactivated measurements", all_deactivated_measurements, "off")
    @printf(io, "  Active measurements: %d (active branch IDs: %d)\n", active_total, active_branch)
    @printf(
      io,
      "  Global  -> observable(num/str)=%s/%s, rank=%d, match=%d, r=m-n=%d, rho=m/n=%.3f, nu=%d, quality=%s\n",
      string(global_obs.numerical_observable),
      string(global_obs.structural_observable),
      global_obs.numerical_rank,
      global_obs.structural_matching,
      global_obs.redundancy,
      global_obs.redundancy_ratio,
      global_obs.dof,
      string(global_obs.quality)
    )

    @printf(
      io,
      "  Local   -> observable(num/str)=%s/%s, rank=%d, match=%d, r=m-n=%d, rho=m/n=%.3f, nu=%d, quality=%s\n",
      string(local_obs.numerical_observable),
      string(local_obs.structural_observable),
      local_obs.numerical_rank,
      local_obs.structural_matching,
      local_obs.redundancy,
      local_obs.redundancy_ratio,
      local_obs.dof,
      string(local_obs.quality)
    )

    crit_num = _critical_labels(meas, local_obs.numerical_critical_measurement_indices)
    crit_str = _critical_labels(meas, local_obs.structural_critical_measurement_indices)
    @printf(io, "  Local critical (numerical): %s\n", isempty(crit_num) ? "none" : join(crit_num, " | "))
    @printf(io, "  Local critical (structural): %s\n", isempty(crit_str) ? "none" : join(crit_str, " | "))

    if !(local_obs.numerical_observable && local_obs.structural_observable)
      @printf(io, "  WARNING: Local redundancy not sufficient anymore at bus %d (after removing branch %s).\n", local_bus, removed_label)
    elseif local_obs.quality == :critical
      @printf(io, "  WARNING: Local area at bus %d is CRITICAL (single measurement outage can break solvability).\n", local_bus)
    end

    # Reset state before WLS so each step is evaluated from a comparable start.
    _reset_non_slack!(net)
    se = runse!(net; maxIte = 12, tol = 1e-6, flatstart = false, jacEps = 1e-6, updateNet = false)
    @printf(io, "  WLS J=r'Wr: %.6e, dof=%d, J within 3σ-band: %s\n\n", se.objectiveJ, se.dof, string(se.jWithin3Sigma))

    return global_obs, local_obs
  end

  groups, adjacent_branches, neighboring_buses = _deactivation_groups(net, meas, target_buses)
  @printf(io, "Focused adjacent branch IDs: %s\n", string(adjacent_branches))
  @printf(io, "Neighbor buses around target bus: %s\n", string(neighboring_buses))
  @printf(io, "Strategy: randomly choose one or more connected, non-boundary target buses; prioritize branch-flow and injection measurements around them; keep deactivated measurements off; continue until local quality becomes critical and, after removing critical measurements, not_observable.\n\n")

  deactivated_set = Set{Int}()
  # Step 0 is baseline with all measurements active.
  _, local_obs = log_step(0, nothing, String[], String[])
  reached_critical = local_obs.quality == :critical
  reached_not_observable = local_obs.quality == :not_observable
  step = 0

  for (_, removed_branch, group_idx) in groups
    removed_labels = String[]
    for meas_idx in group_idx
      # Guard: once a measurement is switched off it must stay off.
      if meas[meas_idx].active
        _set_measurement_active!(meas, meas_idx, false)
        push!(removed_labels, _measurement_label(meas[meas_idx]))
        push!(deactivated_set, meas_idx)
      end
    end
    # Skip no-op groups that only contained already inactive measurements.
    isempty(removed_labels) && continue
    step += 1
    all_deactivated_labels = [_measurement_label(meas[i]) for i in sort!(collect(deactivated_set))]
    _, local_obs = log_step(step, removed_branch, removed_labels, all_deactivated_labels)

    reached_critical = reached_critical || (local_obs.quality == :critical)
    reached_not_observable = reached_not_observable || (local_obs.quality == :not_observable)
    reached_not_observable && break
  end

  @printf(io, "Summary: reached critical=%s, reached not_observable=%s\n", string(reached_critical), string(reached_not_observable))

  return nothing
end

"""
Run the observability example and write the output to a timestamped log file.
"""
function run_and_log_state_estimation_observability_example()
  mkpath(OUTDIR_OBS)
  timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
  logfile = joinpath(OUTDIR_OBS, "run_observability_$(CASEFILE_OBS)_$(timestamp).log")
  open(logfile, "w") do io
    run_state_estimation_observability_example(io)
  end
  println("Observability example finished. Log written to: $logfile")
  return logfile
end

run_and_log_state_estimation_observability_example()

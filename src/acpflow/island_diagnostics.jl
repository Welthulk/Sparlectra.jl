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

_island_bus_type_symbol(node) = getNodeType(node) == Slack ? :REF : getNodeType(node) == PV ? :PV : :PQ

function _ac_island_components(net::Net)
  n = length(net.nodeVec)
  adjacency = [Int[] for _ in 1:n]
  branch_counts = Dict{Int,Int}()
  for br in net.branchVec
    br.status == 0 && continue
    f = Int(br.fromBus)
    t = Int(br.toBus)
    (1 <= f <= n && 1 <= t <= n) || continue
    push!(adjacency[f], t)
    push!(adjacency[t], f)
  end
  seen = falses(n)
  islands = Vector{Vector{Int}}()
  for start in 1:n
    seen[start] && continue
    queue = [start]
    seen[start] = true
    buses = Int[]
    while !isempty(queue)
      bus = popfirst!(queue)
      push!(buses, bus)
      for nbr in adjacency[bus]
        seen[nbr] && continue
        seen[nbr] = true
        push!(queue, nbr)
      end
    end
    push!(islands, sort!(buses))
  end
  bus_to_island = Dict{Int,Int}()
  for (island_id, buses) in pairs(islands)
    for bus in buses
      bus_to_island[bus] = island_id
    end
    branch_counts[island_id] = 0
  end
  for br in net.branchVec
    br.status == 0 && continue
    island_id = get(bus_to_island, Int(br.fromBus), 0)
    island_id > 0 && get(bus_to_island, Int(br.toBus), 0) == island_id && (branch_counts[island_id] += 1)
  end
  return [(island_id = island_id, buses = buses, n_branch = get(branch_counts, island_id, 0)) for (island_id, buses) in pairs(islands)]
end

function _island_ref_bus(net::Net, buses::AbstractVector{<:Integer})
  refs = [Int(bus) for bus in buses if getNodeType(net.nodeVec[Int(bus)]) == Slack]
  !isempty(refs) && return (bus = refs[1], promoted = false)
  pvs = [Int(bus) for bus in buses if getNodeType(net.nodeVec[Int(bus)]) == PV]
  !isempty(pvs) && return (bus = pvs[1], promoted = true)
  return (bus = Int(first(buses)), promoted = true)
end

function _mismatch_status(x)::String
  isnan(x) && return "NaN"
  isinf(x) && return "Inf"
  isfinite(x) && return "finite"
  return "nonfinite"
end

_status_property(status, name::Symbol, default) = status !== nothing && hasproperty(status, name) ? getproperty(status, name) : default

function _status_for_island(performance_profile, island_id::Integer, fallback)
  if performance_profile isa AbstractDict
    statuses = get(performance_profile, :ac_island_solver_statuses, nothing)
    statuses isa AbstractDict && haskey(statuses, Int(island_id)) && return statuses[Int(island_id)]
  end
  return fallback
end

function _write_island_mismatch_history_artifact(output_dir::AbstractString, island_id::Integer, row_status)::String
  artifact = joinpath(output_dir, "ac_island_$(island_id)_mismatch_history.csv")
  history = _status_property(row_status, :mismatch_history, Float64[])
  open(artifact, "w") do io
    println(io, "iteration,max_mismatch")
    if history isa AbstractVector
      for (idx, value) in pairs(history)
        println(io, idx, ",", value)
      end
    end
  end
  return artifact
end

function _format_top_mismatch_rows(rows)::String
  rows isa AbstractVector || return "unavailable"
  isempty(rows) && return "none"
  return join(("$(row.bus_id)/$(row.bus_index)/$(row.equation)/$(row.mismatch)" for row in rows), "; ")
end

function _format_top_mismatch_snapshots(snapshots)::String
  snapshots isa AbstractVector || return "unavailable"
  isempty(snapshots) && return "none"
  return join(("$(snapshot.label)@$(snapshot.iteration): $(_format_top_mismatch_rows(snapshot.rows))" for snapshot in snapshots), " | ")
end

function _island_solver_settings(cfg::PowerFlowConfig)
  return (
    max_iter = cfg.max_iter,
    tol = cfg.tol,
    autodamp = cfg.autodamp,
    autodamp_min = cfg.autodamp_min,
    angle_mode = cfg.start_mode.angle_mode,
    voltage_mode = cfg.start_mode.voltage_mode,
    start_projection = cfg.start_mode.start_projection,
    qlimits_enabled = !cfg.qlimits.ignore_q_limits,
    qlimit_enforcement_mode = cfg.qlimits.enforcement_mode,
    wrong_branch_detection = cfg.wrong_branch_detection,
    wrong_branch_rescue = cfg.wrong_branch_rescue,
    start_current_iteration_enabled = cfg.start_current_iteration.enabled,
    rectangular_workspace_reuse = cfg.rectangular_workspace_reuse,
    rectangular_preallocate_workspace = cfg.rectangular_preallocate_workspace,
    rectangular_workspace_min_buses = cfg.rectangular_workspace_min_buses,
  )
end

function _collect_ac_island_diagnostics(net::Net, cfg::PowerFlowConfig)
  settings = _island_solver_settings(cfg)
  return [begin
    counts = Dict(:PQ => 0, :PV => 0, :REF => 0)
    for bus in component.buses
      key = _island_bus_type_symbol(net.nodeVec[Int(bus)])
      counts[key] = get(counts, key, 0) + 1
    end
    ref = _island_ref_bus(net, component.buses)
    (
      island_id = component.island_id,
      n_bus = length(component.buses),
      n_branch = component.n_branch,
      ref_bus = ref.bus,
      ref_promoted = ref.promoted,
      n_pq = get(counts, :PQ, 0),
      n_pv = get(counts, :PV, 0),
      n_ref = get(counts, :REF, 0),
      settings = settings,
    )
  end for component in _ac_island_components(net)]
end

function _write_ac_island_diagnostics!(net::Net, cfg::PowerFlowConfig, performance_profile; status = nothing, iterations::Union{Nothing,Integer} = nothing)
  cfg.islands.enabled || return NamedTuple()
  output_dir = performance_profile isa AbstractDict ? get(performance_profile, :output_dir, tempdir()) : tempdir()
  mkpath(output_dir)
  pre = _collect_ac_island_diagnostics(net, cfg)
  artifacts = String[]
  summary_path = joinpath(output_dir, "ac_island_solver_summary.csv")
  open(summary_path, "w") do io
    println(io, "island_id,n_bus,n_branch,ref_bus,chosen_ref_bus,n_pq,n_pv,n_ref,ref_promoted,initial_mismatch,first_mismatch,last_mismatch,best_mismatch,min_voltage_magnitude,max_voltage_magnitude,max_angle_step,first_nonfinite_iteration,last_finite_mismatch,worst_bus,worst_equation,iterations,final_mismatch,mismatch_status,final_status,failure_reason,stage,exception_type,exception_message,stacktrace_top,start_projection,autodamp,autodamp_min,max_iter,tol,angle_mode,voltage_mode,qlimits_enabled,qlimit_enforcement_mode,q_limit_processing_status,pv_pq_switching_events,qlimit_active_set_changes,qlimit_reenable_events,guarded_narrow_q_pv_buses,final_pv_voltage_residual,wrong_branch_detection,start_current_iteration_enabled,workspace_reuse,workspace_preallocate")
    for row in pre
      row_status = _status_for_island(performance_profile, row.island_id, status)
      final_mismatch = Float64(_status_property(row_status, :final_mismatch, NaN))
      initial_mismatch = _status_property(row_status, :initial_mismatch, "unavailable")
      best_mismatch = _status_property(row_status, :best_mismatch, _status_property(row_status, :final_mismatch, "unavailable"))
      row_iterations = Int(_status_property(row_status, :iterations, something(iterations, 0)))
      reason = _status_property(row_status, :reason, status === nothing ? :not_attempted : :before_nr_not_run)
      final_status = _status_property(row_status, :status, status === nothing ? :not_attempted : :not_converged)
      stage = String(_status_property(row_status, :stage, status === nothing ? :not_attempted : (row_iterations == 0 ? :pre_solve_validation : (final_status == :converged ? :post_solve_validation : :newton_iteration))))
      exception_type = String(_status_property(row_status, :exception_type, ""))
      exception_message = String(_status_property(row_status, :exception_message, ""))
      stacktrace_top = String(_status_property(row_status, :stacktrace_top, ""))
      first_nonfinite_iteration = _status_property(row_status, :first_nonfinite_iteration, "none")
      last_finite_mismatch = _status_property(row_status, :last_finite_mismatch, isfinite(final_mismatch) ? final_mismatch : "unavailable")
      q_limit_processing_status = row.settings.qlimits_enabled ? String(reason) : "disabled"
      pv_pq_switching_events = _status_property(row_status, :pv_pq_switching_events, 0)
      qlimit_active_set_changes = _status_property(row_status, :qlimit_active_set_changes, 0)
      qlimit_reenable_events = _status_property(row_status, :qlimit_reenable_events, 0)
      guarded_narrow_q_pv_buses = _status_property(row_status, :guarded_narrow_q_pv_buses, 0)
      final_pv_voltage_residual = _status_property(row_status, :final_pv_voltage_residual, "unavailable")
      artifact = joinpath(output_dir, "ac_island_$(row.island_id)_solver.log")
      push!(artifacts, artifact)
      history_artifact = _write_island_mismatch_history_artifact(output_dir, row.island_id, row_status)
      push!(artifacts, history_artifact)
      mismatch_status = _mismatch_status(final_mismatch)
      top_mismatch_rows = _status_property(row_status, :top_mismatch_rows, NamedTuple[])
      top_mismatch_snapshots = _status_property(row_status, :top_mismatch_rows_by_iteration, NamedTuple[])
      open(artifact, "w") do log
        println(log, "island_id: ", row.island_id)
        println(log, "n_bus: ", row.n_bus)
        println(log, "n_branch: ", row.n_branch)
        println(log, "ref_bus: ", row.ref_bus)
        println(log, "chosen_ref_bus: ", row.ref_bus)
        println(log, "ref_promoted: ", row.ref_promoted)
        println(log, "n_pq: ", row.n_pq)
        println(log, "n_pv: ", row.n_pv)
        println(log, "n_ref: ", row.n_ref)
        println(log, "iterations: ", row_iterations)
        println(log, "initial_mismatch: ", initial_mismatch)
        println(log, "first_mismatch: ", _status_property(row_status, :mismatch_history_first, initial_mismatch))
        println(log, "last_mismatch: ", final_mismatch)
        println(log, "best_mismatch: ", best_mismatch)
        println(log, "min_voltage_magnitude: unavailable")
        println(log, "max_voltage_magnitude: unavailable")
        println(log, "max_angle_step: unavailable")
        println(log, "first_nonfinite_iteration: ", first_nonfinite_iteration === nothing ? "none" : first_nonfinite_iteration)
        println(log, "last_finite_mismatch: ", last_finite_mismatch)
        println(log, "worst_bus: ", _status_property(row_status, :worst_mismatch_bus_id, "unavailable"))
        println(log, "worst_bus_internal_index: ", _status_property(row_status, :worst_mismatch_bus_index, "unavailable"))
        println(log, "worst_equation: ", _status_property(row_status, :worst_mismatch_equation, "unavailable"))
        println(log, "worst_mismatch_value: ", _status_property(row_status, :worst_mismatch_value, "unavailable"))
        println(log, "top_mismatch_rows: ", _format_top_mismatch_rows(top_mismatch_rows))
        println(log, "top_mismatch_rows_by_iteration: ", _format_top_mismatch_snapshots(top_mismatch_snapshots))
        println(log, "mismatch_history_csv: ", basename(history_artifact))
        println(log, "mismatch_history_trend: ", _status_property(row_status, :mismatch_history_trend, "unavailable"))
        println(log, "autodamp_step_count: ", _status_property(row_status, :autodamp_step_count, 0))
        println(log, "autodamp_min_alpha: ", _status_property(row_status, :autodamp_min_alpha, "unavailable"))
        println(log, "autodamp_max_alpha: ", _status_property(row_status, :autodamp_max_alpha, "unavailable"))
        println(log, "autodamp_mean_alpha: ", _status_property(row_status, :autodamp_mean_alpha, "unavailable"))
        println(log, "autodamp_floor_hits: ", _status_property(row_status, :autodamp_floor_hits, 0))
        println(log, "autodamp_nonimproving_steps: ", _status_property(row_status, :autodamp_nonimproving_steps, 0))
        println(log, "autodamp_failure: ", _status_property(row_status, :autodamp_failure, false))
        if _status_property(row_status, :autodamp_failure, false) === true
          println(log, "autodamp_failure_reason: line-search/autodamp floor was hit repeatedly while trial mismatches were non-improving")
        end
        println(log, "final_mismatch: ", final_mismatch)
        println(log, "mismatch_status: ", mismatch_status)
        println(log, "final_status: ", final_status)
        println(log, "failure_reason: ", reason)
        println(log, "stage: ", stage)
        println(log, "exception_type: ", exception_type)
        println(log, "exception_message: ", exception_message)
        println(log, "stacktrace_top: ", stacktrace_top)
        println(log, "qlimits_enabled: ", row.settings.qlimits_enabled)
        println(log, "qlimit_enforcement_mode: ", row.settings.qlimit_enforcement_mode)
        println(log, "q_limit_processing_status: ", q_limit_processing_status)
        println(log, "pv_pq_switching_events: ", pv_pq_switching_events)
        println(log, "qlimit_active_set_changes: ", qlimit_active_set_changes)
        println(log, "qlimit_reenable_events: ", qlimit_reenable_events)
        println(log, "guarded_narrow_q_pv_buses: ", guarded_narrow_q_pv_buses)
        println(log, "final_pv_voltage_residual: ", final_pv_voltage_residual)
        row_status === nothing && println(log, "mismatch_data_unavailable_reason: island was not attempted before diagnostics were written")
        println(log, "solver_settings: ", row.settings)
      end
      println(io, join((row.island_id, row.n_bus, row.n_branch, row.ref_bus, row.ref_bus, row.n_pq, row.n_pv, row.n_ref, row.ref_promoted, initial_mismatch, "unavailable", final_mismatch, best_mismatch, "unavailable", "unavailable", "unavailable", first_nonfinite_iteration, last_finite_mismatch, "unavailable", "unavailable", row_iterations, final_mismatch, mismatch_status, final_status, reason, stage, _csv_field(exception_type, ','), _csv_field(exception_message, ','), _csv_field(stacktrace_top, ','), row.settings.start_projection, row.settings.autodamp, row.settings.autodamp_min, row.settings.max_iter, row.settings.tol, row.settings.angle_mode, row.settings.voltage_mode, row.settings.qlimits_enabled, row.settings.qlimit_enforcement_mode, q_limit_processing_status, pv_pq_switching_events, qlimit_active_set_changes, qlimit_reenable_events, guarded_narrow_q_pv_buses, final_pv_voltage_residual, row.settings.wrong_branch_detection, row.settings.start_current_iteration_enabled, row.settings.rectangular_workspace_reuse, row.settings.rectangular_preallocate_workspace), ','))
    end
  end
  artifacts_summary = (summary_path, artifacts...)
  performance_profile isa AbstractDict && (performance_profile[:ac_island_diagnostics] = pre; performance_profile[:ac_island_artifacts] = artifacts_summary)
  return (island_count = length(pre), diagnostics = pre, artifacts = artifacts_summary)
end

function _append_island_failure_message(status, performance_profile)
  performance_profile isa AbstractDict || return status
  artifacts = get(performance_profile, :ac_island_artifacts, ())
  diagnostics = get(performance_profile, :ac_island_diagnostics, ())
  isempty(artifacts) && return status
  status.final_converged && return status
  failed_island_id = hasproperty(status, :island_id) ? status.island_id : nothing
  first_island = failed_island_id === nothing ? (isempty(diagnostics) ? nothing : first(diagnostics)) : only(row for row in diagnostics if row.island_id == failed_island_id)
  first_island === nothing && return status
  artifact = basename(first(artifacts))
  iterations = hasproperty(status, :iterations) ? status.iterations : 0
  text = string(
    "AC island ", first_island.island_id, " power-flow solve failed:\n",
    "  buses=", first_island.n_bus, " branches=", first_island.n_branch, " ref=", first_island.ref_bus, "\n",
    "  bus_types: PV=", first_island.n_pv, " PQ=", first_island.n_pq, " REF=", first_island.n_ref, "\n",
    "  iterations=", iterations, " final_mismatch=", status.final_mismatch, "\n",
    "  mismatch_status=", _mismatch_status(status.final_mismatch), "\n",
    "  reason=", status.reason, "\n",
    "  start_projection=", first_island.settings.start_projection, "\n",
    "  stage=", (hasproperty(status, :stage) ? String(status.stage) : (iterations == 0 ? "before_nr" : "during_nr")), "\n",
    "  artifact=", artifact,
  )
  return merge(status, (reason_text = text,))
end

function _island_status_value(row_status, key::Symbol, default)
  row_status === nothing && return default
  return hasproperty(row_status, key) ? getproperty(row_status, key) : default
end

function _islandwise_failure_message(performance_profile)::Union{Nothing,String}
  performance_profile isa AbstractDict || return nothing
  artifacts = get(performance_profile, :ac_island_artifacts, ())
  diagnostics = get(performance_profile, :ac_island_diagnostics, ())
  statuses = get(performance_profile, :ac_island_solver_statuses, nothing)
  isempty(diagnostics) && return nothing
  statuses === nothing && return nothing

  lines = String["Island-wise power-flow failed:"]
  qlimits_enabled = false
  qlimit_modes = Set{String}()
  start_current_iteration_enabled = false
  failed_descriptions = String[]
  for row in diagnostics
    row_status = _status_for_island(performance_profile, row.island_id, nothing)
    final_status = _island_status_value(row_status, :status, :not_attempted)
    reason = _island_status_value(row_status, :reason, final_status == :converged ? :none : :unavailable)
    iterations = _island_status_value(row_status, :iterations, 0)
    final_mismatch = _island_status_value(row_status, :final_mismatch, NaN)
    if final_status == :converged
      push!(lines, "  island $(row.island_id): converged, iterations=$(iterations)")
    else
      push!(lines, "  island $(row.island_id): $(final_status), iterations=$(iterations), final_mismatch=$(final_mismatch), reason=$(reason)")
      push!(failed_descriptions, "island $(row.island_id) $(final_status)")
    end
    qlimits_enabled |= row.settings.qlimits_enabled
    start_current_iteration_enabled |= row.settings.start_current_iteration_enabled
    row.settings.qlimits_enabled && push!(qlimit_modes, String(row.settings.qlimit_enforcement_mode))
  end
  if qlimits_enabled
    modes = isempty(qlimit_modes) ? "unknown" : join(sort!(collect(qlimit_modes)), ", ")
    push!(lines, "Q-limit handling was enabled for this run (mode=$(modes)); $(join(failed_descriptions, " and ")) failed under Q-limit enforcement.")
  elseif start_current_iteration_enabled
    push!(lines, "Island run with current-iteration start: qlimits.enabled=false; start_current_iteration.enabled=true.")
  else
    push!(lines, "Pure island NR diagnostic run: qlimits.enabled=false and start_current_iteration.enabled=false.")
  end
  artifact_names = isempty(artifacts) ? ["ac_island_solver_summary.csv", "ac_island_<id>_solver.log"] : basename.(collect(artifacts))
  push!(lines, "See $(join(artifact_names, " and ")).")
  return join(lines, '\n')
end

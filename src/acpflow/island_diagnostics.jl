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

function _write_ac_island_diagnostics!(net::Net, cfg::PowerFlowConfig, performance_profile; status = nothing)
  cfg.islands.enabled || return NamedTuple()
  output_dir = performance_profile isa AbstractDict ? get(performance_profile, :output_dir, tempdir()) : tempdir()
  mkpath(output_dir)
  pre = _collect_ac_island_diagnostics(net, cfg)
  final_mismatch = status === nothing || !hasproperty(status, :final_mismatch) ? NaN : Float64(status.final_mismatch)
  iterations = status === nothing || !hasproperty(status, :iterations) ? 0 : Int(status.iterations)
  reason = status === nothing || !hasproperty(status, :reason) ? :before_nr_not_run : status.reason
  final_status = status === nothing || !hasproperty(status, :status) ? :not_run : status.status
  stage = status !== nothing && hasproperty(status, :stage) ? String(status.stage) : status === nothing ? "pre_solve_validation" : (iterations == 0 ? "pre_solve_validation" : (final_status == :converged ? "post_solve_validation" : "newton_iteration"))
  exception_type = status !== nothing && hasproperty(status, :exception_type) ? String(status.exception_type) : ""
  exception_message = status !== nothing && hasproperty(status, :exception_message) ? String(status.exception_message) : ""
  stacktrace_top = status !== nothing && hasproperty(status, :stacktrace_top) ? String(status.stacktrace_top) : ""
  artifacts = String[]
  summary_path = joinpath(output_dir, "ac_island_solver_summary.csv")
  open(summary_path, "w") do io
    println(io, "island_id,n_bus,n_branch,ref_bus,chosen_ref_bus,n_pq,n_pv,n_ref,ref_promoted,initial_mismatch,first_mismatch,last_mismatch,best_mismatch,min_voltage_magnitude,max_voltage_magnitude,max_angle_step,first_nonfinite_iteration,worst_bus,worst_equation,iterations,final_mismatch,mismatch_status,final_status,failure_reason,stage,exception_type,exception_message,stacktrace_top,start_projection,autodamp,autodamp_min,max_iter,tol,angle_mode,voltage_mode,qlimits_enabled,qlimit_enforcement_mode,wrong_branch_detection,start_current_iteration_enabled,workspace_reuse,workspace_preallocate")
    for row in pre
      artifact = joinpath(output_dir, "ac_island_$(row.island_id)_solver.log")
      push!(artifacts, artifact)
      mismatch_status = _mismatch_status(final_mismatch)
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
        println(log, "iterations: ", iterations)
        println(log, "initial_mismatch: unavailable")
        println(log, "first_mismatch: unavailable")
        println(log, "last_mismatch: ", final_mismatch)
        println(log, "best_mismatch: unavailable")
        println(log, "min_voltage_magnitude: unavailable")
        println(log, "max_voltage_magnitude: unavailable")
        println(log, "max_angle_step: unavailable")
        println(log, "first_nonfinite_iteration: none")
        println(log, "worst_bus: unavailable")
        println(log, "worst_equation: unavailable")
        println(log, "final_mismatch: ", final_mismatch)
        println(log, "mismatch_status: ", mismatch_status)
        println(log, "final_status: ", final_status)
        println(log, "failure_reason: ", reason)
        println(log, "stage: ", stage)
        println(log, "exception_type: ", exception_type)
        println(log, "exception_message: ", exception_message)
        println(log, "stacktrace_top: ", stacktrace_top)
        println(log, "solver_settings: ", row.settings)
      end
      println(io, join((row.island_id, row.n_bus, row.n_branch, row.ref_bus, row.ref_bus, row.n_pq, row.n_pv, row.n_ref, row.ref_promoted, "unavailable", "unavailable", final_mismatch, "unavailable", "unavailable", "unavailable", "unavailable", "none", "unavailable", "unavailable", iterations, final_mismatch, mismatch_status, final_status, reason, stage, _csv_field(exception_type, ','), _csv_field(exception_message, ','), _csv_field(stacktrace_top, ','), row.settings.start_projection, row.settings.autodamp, row.settings.autodamp_min, row.settings.max_iter, row.settings.tol, row.settings.angle_mode, row.settings.voltage_mode, row.settings.qlimits_enabled, row.settings.qlimit_enforcement_mode, row.settings.wrong_branch_detection, row.settings.start_current_iteration_enabled, row.settings.rectangular_workspace_reuse, row.settings.rectangular_preallocate_workspace), ','))
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

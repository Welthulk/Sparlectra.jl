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

# file: src/acpflow/apslf_execution.jl
# purpose: solver wiring for power_flow.solver = :apslf: runs framework solves
#          through the external analytic power-series solver per AC island
#          while keeping the rectangular status and diagnostics contracts
# Solver-selection wiring for power_flow.solver = :apslf: routes the central
# framework run through the external-solver bridge (buildPfModel ->
# solvePf(ApslfSolver) -> applyPfSolution!, via runpf_external!) instead of the
# internal rectangular Newton-Raphson solver. Mirrors the internal path's
# parallel-link merge (_merged_pf_net) and per-island handling
# (detect_ac_islands/_prepare_island_net/_sync_island_solution!) so status
# classification, timing, and the generic AC island diagnostics writer stay on
# the same contracts as the rectangular route (rectangular_pf_status,
# performance_profile[:ac_island_solver_statuses]).

function _apslf_solver_from_config(pf_cfg::PowerFlowConfig)
  acfg = pf_cfg.apslf
  return apslf_solver(order = acfg.order, nr_polish = acfg.nr_polish, convergence_radius = acfg.convergence_radius)
end

# Status fields carrying the two APSLF radii so the result header, the run
# metadata and the runs page can show them like the Jacobian condition:
# the Padé margin dmin (kept in apslf_convergence_radius, the field name the
# workshops read) and the coefficient-growth radius of the series. `sol` is
# the PFSolution of solvePf on `snet`, the net the model was built from (the
# run's net, or an island net with its own bus numbering); the bus fields are
# case bus numbers (busOrigIdxDict, the numbers of the case file), never the
# internal node index, which pointed at the wrong bus on every case whose bus
# numbers are not 1..n (case2848rte: 1966 for the case bus 1987).
function _apslf_radius_status(sol, snet::Net)::NamedTuple
  st = sol.meta.stability
  sr = sol.meta.series_radius
  case_bus(b::Int) = b >= 1 ? _qlimit_original_bus_id(snet, b) : 0
  return _apslf_radius_fields(st.enabled, Float64(st.dmin), case_bus(st.bus), String(st.level), Float64(sr.radius), case_bus(sr.bus))
end

function _apslf_radius_fields(evaluated::Bool, dmin::Float64, bus::Int, level::String, series_radius::Float64, series_bus::Int)::NamedTuple
  return (
    apslf_convergence_radius = dmin,
    apslf_convergence_bus = bus,
    apslf_convergence_level = level,
    apslf_pade_margin_evaluated = evaluated,
    apslf_series_radius = series_radius,
    apslf_series_radius_bus = series_bus,
    apslf_convergence_line = _apslf_radius_line(evaluated, dmin, bus, level, series_radius, series_bus),
  )
end

function _apslf_radius_line(evaluated::Bool, dmin::Float64, bus::Int, level::String, series_radius::Float64, series_bus::Int)::String
  pade = if !evaluated
    "Pade margin not evaluated (power_flow.apslf.convergence_radius: false)"
  elseif !isfinite(dmin)
    "Pade margin not available"
  else
    string("Pade margin dmin = ", round(dmin; sigdigits = 3), " (nearest Pade pole to s = 1 at bus ", bus, ", level ", level, ")")
  end
  series = if isnan(series_radius)
    "series radius not available"
  elseif isinf(series_radius)
    "series radius unbounded (the voltage coefficients vanish)"
  else
    string("series radius R = ", round(series_radius; sigdigits = 3), " (root test on the voltage coefficients, bus ", series_bus, ")")
  end
  return string(pade, "; ", series)
end

# The run's radii over several islands: each of the two is the smallest one
# of any island, with the bus of that island (they may come from different
# islands). `a` is nothing before the first island.
function _apslf_worst_radius(a, b::NamedTuple)::NamedTuple
  a === nothing && return b
  smaller(x, y) = isfinite(y) && (!isfinite(x) || y < x)
  pade_b = smaller(a.apslf_convergence_radius, b.apslf_convergence_radius)
  series_b = isnan(a.apslf_series_radius) ? !isnan(b.apslf_series_radius) : (!isnan(b.apslf_series_radius) && b.apslf_series_radius < a.apslf_series_radius)
  p = pade_b ? b : a
  q = series_b ? b : a
  return _apslf_radius_fields(a.apslf_pade_margin_evaluated, p.apslf_convergence_radius, p.apslf_convergence_bus, p.apslf_convergence_level, q.apslf_series_radius, q.apslf_series_radius_bus)
end

# Register the machines AnalyticLoadFlow clamped at a reactive limit in the
# net's Q-limit log, so the result table shows them as PQ* with the binding
# side and the Q-V check judges them (the rectangular solver logs the same
# events through active_set_q_limits!). The clamps carry the bus index of the
# net the model was built from; `buses` maps an island net's local index to
# the run's index (the island's row.buses), `nothing` when the model was built
# from a net with the run's numbering. Returns the number of clamps, which is
# the PV->PQ switching count of the status (pv_pq_switching_events).
function _apslf_record_clamps!(net::Net, sol, iters::Int; buses::Union{Nothing,AbstractVector{Int}} = nothing)::Int
  for (bus, side) in sol.meta.qlimit_clamps
    logQLimitHit!(net, iters, buses === nothing ? bus : buses[bus], side)
  end
  return length(sol.meta.qlimit_clamps)
end

# Builds a rectangular_pf_status-compatible NamedTuple. Only the fields actually
# read by _rectangular_run_status/_rectangular_status_diagnostics are set;
# rectangular-NR-specific diagnostics (autodamp, active-set changes and
# re-enable events, start-current-iteration, ...) do not apply to APSLF and
# are left to their generic defaults there. The PV->PQ count
# (pv_pq_switching_events) does apply: the callers pass the clamps of the
# solve through `extra`.
function _apslf_pf_status(converged::Bool, final_mismatch::Float64, iterations::Int; extra::NamedTuple = NamedTuple())
  base = (
    numerical_converged = converged,
    nr_converged = converged,
    active_set_converged = converged,
    q_limit_active_set_ok = converged,
    final_converged = converged,
    status = converged ? :converged : :not_converged,
    reason = converged ? :none : :apslf_not_converged,
    reason_text = converged ? "APSLF analytic power-series solve converged." : "APSLF analytic power-series solve did not converge within power_flow.tol.",
    final_mismatch = final_mismatch,
    nr_final_mismatch = final_mismatch,
    iterations = iterations,
    solver = :apslf,
    stage = :apslf_solve,
    exception_type = "",
    exception_message = "",
    stacktrace_top = "",
  )
  return merge(base, extra)
end

"""
    _run_apslf_powerflow!(net, pf_cfg; verbose=0, performance_profile=nothing) -> (iterations::Int, erg::Int)

Solve `net` with the AnalyticLoadFlow.jl-backed APSLF solver
(`power_flow.solver = :apslf`) instead of the internal rectangular
Newton-Raphson path. Mirrors the internal path's parallel-link merge
(`_merged_pf_net`) and per-island handling (`detect_ac_islands`,
`_prepare_island_net`, `_sync_island_solution!`): one `PFModel` is built and
solved per active AC island via `runpf_external!`. Writes a
`rectangular_pf_status`-compatible status onto `net` so result building,
reporting, logs, and the generic AC island diagnostics writer treat this run
on the same contracts as the NR route.

`power_flow.qlimits.ignore_q_limits` is honored (Q-limits are passed to the
solver unless explicitly disabled); all other rectangular-only start/damping/
Q-limit-switching options do not apply and are ignored, since APSLF always
starts from its own canonical analytic germ and performs its own internal
PV/Q-limit handling.

Callers must reject active controllers before calling this function; it does
not check for them itself. That covers the outer-loop tap/PST controllers and
the voltage-dependent Q(U)/P(U) controllers, which are not outer-loop: they
act inside the rectangular Newton step and have no counterpart here.
"""
function _run_apslf_powerflow!(net::Net, pf_cfg::PowerFlowConfig; verbose::Int = 0, performance_profile = nothing)::Tuple{Int,Int}
  solver = _apslf_solver_from_config(pf_cfg)
  include_limits = !pf_cfg.qlimits.ignore_q_limits

  wnet, reps, has_merges = _merged_pf_net(net)
  refreshBusTypesFromProsumers!(wnet)
  if has_merges && has_voltage_dependent_control(wnet)
    error("power_flow.solver=apslf: voltage-dependent injections, including P(U)/Q(U) controllers and bus_shunt_model=voltage_dependent_injection, are not supported with active-link merge handling. Disable merges or use a topology without internal isolated buses.")
  end
  island_report = detect_ac_islands(wnet)
  # same rule as the rectangular entry: one angle reference per synchronous
  # island, rejected early with an actionable message
  _validate_single_reference_per_island!(wnet, island_report)
  multi_island = length(island_report.rows) > 1 && any(row -> row.n_branch > 0, island_report.rows)

  sync_merges_back! = () -> begin
    for i in eachindex(net.nodeVec)
      src = wnet.nodeVec[reps[i]]
      net.nodeVec[i]._vm_pu = src._vm_pu
      net.nodeVec[i]._va_deg = src._va_deg
    end
    updateShuntPowers!(net = net)
  end

  # a fresh Q-limit log per solve, like the rectangular entry
  resetQLimitLog!(net)
  snapshotPVQLimits!(net)
  if !multi_island
    iters, status, sol = runpf_external!(wnet, solver; tol = pf_cfg.tol, flatstart = wnet.flatstart, include_limits = include_limits, verbose = verbose)
    switches = _apslf_record_clamps!(net, sol, iters)
    _set_rectangular_pf_status!(net, _apslf_pf_status(status == 0, sol.residual_inf, iters; extra = merge(_apslf_radius_status(sol, wnet), (pv_pq_switching_events = switches,))))
    has_merges && sync_merges_back!()
    return iters, status
  end

  # the run's output directory only (see the rectangular solver)
  artifact_dir = performance_profile isa AbstractDict ? get(performance_profile, :output_dir, nothing) : nothing
  if artifact_dir !== nothing
    try
      mkpath(String(artifact_dir))
      island_artifact = joinpath(String(artifact_dir), "ac_islands.csv")
      write_ac_island_report(island_artifact, island_report)
      println("AC island diagnostic artifact: ", island_artifact)
    catch err
      @warn "Unable to write AC island diagnostic artifact" exception = (err, catch_backtrace())
    end
  end
  pf_cfg.islands_enabled || error(AC_ISLAND_DISABLED_MESSAGE)
  pf_cfg.islands_mode === :solve_independent || error("Unsupported power_flow.islands.mode=$(pf_cfg.islands_mode).")
  pf_cfg.islands_reference_policy === :matpower_like || error("Unsupported power_flow.islands.reference_policy=$(pf_cfg.islands_reference_policy).")
  _validate_island_references!(island_report)

  total_iters = 0
  total_switches = 0
  first_failure = nothing
  # the smallest radii over the islands are the run's radii
  worst_radius = nothing
  island_statuses = Dict{Int,Any}()
  performance_profile isa AbstractDict && (performance_profile[:ac_island_solver_statuses] = island_statuses)

  for row in island_report.rows
    if first_failure !== nothing && !pf_cfg.islands.diagnostic_continue_after_failure
      skipped_status = _apslf_pf_status(false, NaN, 0; extra = (island_id = row.island_id, reason = :skipped_after_previous_failure, status = :skipped_after_previous_failure, stage = :skipped_after_previous_failure, reason_text = "previous island failed and diagnostic continuation is disabled"))
      island_statuses[Int(row.island_id)] = skipped_status
      _set_rectangular_pf_status!(net, skipped_status)
      continue
    end
    local inet = nothing
    local it = 0
    try
      inet = _prepare_island_net(wnet, row)
      local status
      local sol
      it, status, sol = runpf_external!(inet, solver; tol = pf_cfg.tol, flatstart = inet.flatstart, include_limits = include_limits, verbose = verbose)
      total_iters += it
      # an island net numbers its buses 1..n_island (busIdx_net of its model
      # is island-local): the clamps go to the run's net through row.buses,
      # the radius buses to case numbers through the island's busOrigIdxDict
      island_switches = _apslf_record_clamps!(net, sol, it; buses = row.buses)
      total_switches += island_switches
      worst_radius = _apslf_worst_radius(worst_radius, _apslf_radius_status(sol, inet))
      if status != 0
        failure_status = _apslf_pf_status(false, sol.residual_inf, it; extra = (island_id = row.island_id, pv_pq_switching_events = island_switches))
        island_statuses[Int(row.island_id)] = failure_status
        _set_rectangular_pf_status!(net, failure_status)
        err = ErrorException("AC island $(row.island_id) power-flow solve failed (solver=apslf):\n  buses=$(row.n_bus) branches=$(row.n_branch) ref=$(row.chosen_ref_bus)\n  bus_types: PV=$(row.n_pv) PQ=$(row.n_pq) REF=$(row.n_ref)\n  iterations=$(it) residual_inf=$(sol.residual_inf)")
        first_failure === nothing && (first_failure = err)
        pf_cfg.islands.diagnostic_continue_after_failure || throw(err)
        continue
      end
      all(isfinite(something(wnode._vm_pu, NaN)) && isfinite(something(wnode._va_deg, NaN)) for wnode in inet.nodeVec) || error("AC island $(row.island_id) produced nonfinite voltage results.")
      success_status = _apslf_pf_status(true, sol.residual_inf, it; extra = (island_id = row.island_id, pv_pq_switching_events = island_switches))
      island_statuses[Int(row.island_id)] = success_status
      _sync_island_solution!(wnet, inet, row)
    catch err
      frames = stacktrace(catch_backtrace())
      top = isempty(frames) ? "" : sprint(show, first(frames))
      failure_status = _apslf_pf_status(false, NaN, it; extra = (island_id = row.island_id, reason = :solver_exception, status = :failed, exception_type = nameof(typeof(err)), exception_message = sprint(showerror, err), stacktrace_top = top))
      island_statuses[Int(row.island_id)] = failure_status
      _set_rectangular_pf_status!(net, failure_status)
      first_failure === nothing && (first_failure = err)
      pf_cfg.islands.diagnostic_continue_after_failure || rethrow()
    end
  end

  if first_failure !== nothing
    island_message = _islandwise_failure_message(performance_profile)
    throw(ErrorException(island_message === nothing ? sprint(showerror, first_failure) : island_message))
  end

  updateShuntPowers!(net = wnet)
  island_final_mismatches = Float64[
    Float64(getproperty(status, :final_mismatch)) for status in values(island_statuses)
    if hasproperty(status, :final_mismatch) && isfinite(Float64(getproperty(status, :final_mismatch)))
  ]
  aggregate_final_mismatch = isempty(island_final_mismatches) ? NaN : maximum(island_final_mismatches)
  aggregate_status = _apslf_pf_status(true, aggregate_final_mismatch, total_iters; extra = merge((
    reason_text = "All AC islands converged independently.",
    island_wise_all_converged = true,
    stage = :island_wise_complete,
    pv_pq_switching_events = total_switches,
  ), worst_radius === nothing ? NamedTuple() : worst_radius))
  _set_rectangular_pf_status!(net, aggregate_status)
  performance_profile isa AbstractDict && (performance_profile[:island_wise_all_converged] = true)

  if has_merges
    sync_merges_back!()
  elseif wnet !== net
    for i in eachindex(net.nodeVec)
      net.nodeVec[i]._vm_pu = wnet.nodeVec[i]._vm_pu
      net.nodeVec[i]._va_deg = wnet.nodeVec[i]._va_deg
    end
  end
  return total_iters, 0
end

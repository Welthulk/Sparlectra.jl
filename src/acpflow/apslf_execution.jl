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

# Status fields for AnalyticLoadFlow's diagnosis of a non-converged solve
# (solvePf meta.diagnosis): the reason, the nearest Padé pole (modulus, case
# bus, conjugate pair or not), the change of the Padé value at s = 1 between
# order N - 2 and N, and one line that carries all of it. The line also
# replaces the generic reason_text, so run.log ("Primary reason") and the
# message of the Web UI result show it. Empty for a converged solve.
function _apslf_diagnosis_status(sol, snet::Net)::NamedTuple
  d = sol.meta.diagnosis
  d === nothing && return NamedTuple()
  bus = d.pole_bus >= 1 ? _qlimit_original_bus_id(snet, d.pole_bus) : 0
  line = _apslf_diagnosis_line(d.reason, d.message, d.pole_min_abs, d.pole_pair, bus, d.pade_change, d.order)
  return (
    apslf_diagnosis_reason = d.reason,
    apslf_pole_min_abs = d.pole_min_abs,
    apslf_pole_min_bus = bus,
    apslf_pole_pair = d.pole_pair,
    apslf_pade_change = d.pade_change,
    apslf_diagnosis_line = line,
    reason_text = string("APSLF analytic power-series solve did not converge within power_flow.tol. ", line),
  )
end

function _apslf_diagnosis_line(reason::Symbol, message::String, pole_abs::Float64, pair::Bool, bus::Int, pade_change::Float64, order::Int)::String
  pole = isfinite(pole_abs) ? string("nearest Pade pole |s| = ", round(pole_abs; sigdigits = 3), pair ? " (conjugate pair)" : " (single)", " at bus ", bus) : "no Pade pole found"
  change = isnan(pade_change) ? "Pade change not evaluated (order below 4)" : string("Pade change order ", order - 2, " to ", order, " = ", round(pade_change; sigdigits = 3), " pu")
  # the message of :pade_not_converged already names the change
  return reason === :pade_not_converged ? string("APSLF diagnosis: ", reason, ", ", message, "; ", pole) : string("APSLF diagnosis: ", reason, ", ", message, "; ", pole, "; ", change)
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
# solve through `extra`. q_limit_active_set_ok starts as `converged`; the
# callers pass the band check of the solved net (_apslf_qband_status)
# through `extra`, which sets it false when a machine ended outside its band.
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

# The limit check of a converged APSLF run, on the solved net after the
# write-back: every voltage-regulated machine inside its reactive band, read
# from the rows of the Q-limit stability check (qlimitStabilityReport, the
# generator Q of the result against the bus band of getQLimits_pu). It
# judges the result, not the band AnalyticLoadFlow was given: with the
# generator band handed over against the bus injection, sp_case118 ended
# with the PV machines at buses 4 and 85 above Qmax and four clamped
# machines beyond their band, and the run reported the limits as ok. The
# same failure now shows as a failed limit status, as in the rectangular
# solver (remaining PV Q-limit violations: converged_limits_failed, erg 1).
# `buses`, when given, restricts the check to these buses of `net` (one
# island). Returns `extra` fields for _apslf_pf_status: empty when the run
# did not converge, ran without Q-limits, or every machine is in its band.
function _apslf_qband_status(net::Net, converged::Bool, include_limits::Bool; rows = nothing, buses = nothing)::NamedTuple
  (converged && include_limits) || return NamedTuple()
  rows === nothing && (rows = qlimitStabilityReport(net))
  bad = [r for r in rows if !r.in_band && (buses === nothing || r.bus in buses)]
  isempty(bad) && return NamedTuple()
  text = string("APSLF solution leaves ", length(bad), " voltage-regulated machine(s) outside their reactive band: ", join((_qlimit_stability_row_text(r) for r in bad), "; "))
  return (q_limit_active_set_ok = false, active_set_converged = false, final_converged = false, status = :converged_limits_failed,
    reason = :remaining_pv_q_limit_violations, reason_text = text, pv_q_limit_violations = length(bad))
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

No controller is evaluated here: the outer-loop tap/PST controllers are not
run, and the Q(U)/P(U) controllers and voltage-dependent shunts, which act
inside the rectangular Newton step, have no counterpart; every unit enters
with its static set point, with and without closed links. The caller
(`_execute_sparlectra_powerflow!`) names what is left out in one warning.
"""
function _run_apslf_powerflow!(net::Net, pf_cfg::PowerFlowConfig; verbose::Int = 0, performance_profile = nothing)::Tuple{Int,Int}
  solver = _apslf_solver_from_config(pf_cfg)
  include_limits = !pf_cfg.qlimits.ignore_q_limits

  # The merged copy carries the static set points of every unit, as an
  # unmerged net does: buildPfModel reads the prosumer schedules (the
  # prosumers follow their bus onto the representative) and the Y-bus and
  # never a Q(U)/P(U) characteristic, on any net. The refusal of
  # voltage-dependent injections on a merged net that stood here (a391e34f)
  # outlived the caller's own refusal, which became a warning in 0.40.0, so
  # a closed coupler turned that warning into an error.
  wnet, reps, has_merges = _merged_pf_net(net)
  refreshBusTypesFromProsumers!(wnet)
  island_report = detect_ac_islands(wnet)
  # same rule as the rectangular entry: one angle reference per synchronous
  # island, rejected early with an actionable message
  _validate_single_reference_per_island!(wnet, island_report)
  multi_island = length(island_report.rows) > 1 && any(row -> row.n_branch > 0, island_report.rows)

  # The solve writes slack P/Q and PV Q onto the merged copy, and the
  # clamps are recorded there too; all of it goes back split over the
  # members as in the rectangular entry (_split_merged_bus_powers!, which
  # rebuilds the Q-limit log of the net from the copy). Before 0.60.0 only
  # the voltages went back and the solved generation stayed on the copy:
  # sp_case60 (coupler Neubach_110 to Neubach_110b) reported 0.0 MVAr at
  # Birkkamp_110, where the same net with the coupler resolved by hand gives
  # 7.258 MVAr, and the generation of the slack Altau_220 differed by 19.99
  # (MW or MVAr, the larger of the two). converged = false: APSLF applies no
  # Q(U)/P(U) characteristic, so a member keeps its schedule instead of the
  # characteristic at the voltage.
  sync_merges_back! = () -> begin
    for i in eachindex(net.nodeVec)
      src = wnet.nodeVec[reps[i]]
      net.nodeVec[i]._vm_pu = src._vm_pu
      net.nodeVec[i]._va_deg = src._va_deg
    end
    _split_merged_bus_powers!(net, wnet, reps; converged = false)
    updateShuntPowers!(net = net)
  end

  # a fresh Q-limit log per solve, like the rectangular entry; the clamps
  # are recorded on the solved net (the merged copy, which is a copy of the
  # net with the log of the last solve)
  resetQLimitLog!(net)
  has_merges && resetQLimitLog!(wnet)
  snapshotPVQLimits!(net)
  if !multi_island
    iters, status, sol = runpf_external!(wnet, solver; tol = pf_cfg.tol, flatstart = wnet.flatstart, include_limits = include_limits, verbose = verbose)
    switches = _apslf_record_clamps!(wnet, sol, iters)
    diag = _apslf_diagnosis_status(sol, wnet)
    run_status = _apslf_pf_status(status == 0, sol.residual_inf, iters; extra = merge(_apslf_radius_status(sol, wnet), diag, (pv_pq_switching_events = switches,)))
    # set before the merged copy goes back: the split reads the status of
    # the net (distributed slack shares)
    _set_rectangular_pf_status!(net, run_status)
    has_merges && sync_merges_back!()
    # the band check reads the generator Q of the net after the split
    band = _apslf_qband_status(net, status == 0, include_limits)
    isempty(band) && return iters, status
    _set_rectangular_pf_status!(net, merge(run_status, band))
    return iters, 1
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
      island_switches = _apslf_record_clamps!(wnet, sol, it; buses = row.buses)
      total_switches += island_switches
      worst_radius = _apslf_worst_radius(worst_radius, _apslf_radius_status(sol, inet))
      if status != 0
        island_diag = _apslf_diagnosis_status(sol, inet)
        failure_status = _apslf_pf_status(false, sol.residual_inf, it; extra = merge(island_diag, (island_id = row.island_id, pv_pq_switching_events = island_switches)))
        island_statuses[Int(row.island_id)] = failure_status
        _set_rectangular_pf_status!(net, failure_status)
        diag_text = hasproperty(island_diag, :apslf_diagnosis_line) ? string("\n  ", island_diag.apslf_diagnosis_line) : ""
        err = ErrorException("AC island $(row.island_id) power-flow solve failed (solver=apslf):\n  buses=$(row.n_bus) branches=$(row.n_branch) ref=$(row.chosen_ref_bus)\n  bus_types: PV=$(row.n_pv) PQ=$(row.n_pq) REF=$(row.n_ref)\n  iterations=$(it) residual_inf=$(sol.residual_inf)$(diag_text)")
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

  # the band check on the net with the results of every island, once; each
  # island status gets the rows of its own buses (row.buses are buses of
  # the merged copy, reps maps a bus of the net onto them)
  band_rows = include_limits ? qlimitStabilityReport(net) : NamedTuple[]
  for row in island_report.rows
    id = Int(row.island_id)
    haskey(island_statuses, id) || continue
    row_buses = Set{Int}(row.buses)
    island_buses = Set{Int}(b for b in eachindex(net.nodeVec) if reps[b] in row_buses)
    island_band = _apslf_qband_status(net, Bool(island_statuses[id].final_converged), include_limits; rows = band_rows, buses = island_buses)
    isempty(island_band) || (island_statuses[id] = merge(island_statuses[id], island_band))
  end
  band = _apslf_qband_status(net, true, include_limits; rows = band_rows)
  isempty(band) && return total_iters, 0
  _set_rectangular_pf_status!(net, merge(aggregate_status, band))
  return total_iters, 1
end

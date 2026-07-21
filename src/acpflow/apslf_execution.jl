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
  return apslf_solver(order = acfg.order, use_pade = acfg.use_pade, nr_polish = acfg.nr_polish)
end

# Builds a rectangular_pf_status-compatible NamedTuple. Only the fields actually
# read by _rectangular_run_status/_rectangular_status_diagnostics are set;
# rectangular-NR-specific diagnostics (autodamp, qlimit active-set switching,
# start-current-iteration, ...) do not apply to APSLF and are left to their
# generic defaults there.
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

Callers must reject active outer-loop controllers (tap/PST/Q(U)/P(U)) before
calling this function; it does not check for them itself.
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
  multi_island = length(island_report.rows) > 1 && any(row -> row.n_branch > 0, island_report.rows)

  sync_merges_back! = () -> begin
    for i in eachindex(net.nodeVec)
      src = wnet.nodeVec[reps[i]]
      net.nodeVec[i]._vm_pu = src._vm_pu
      net.nodeVec[i]._va_deg = src._va_deg
    end
    updateShuntPowers!(net = net)
  end

  if !multi_island
    iters, status, sol = runpf_external!(wnet, solver; tol = pf_cfg.tol, flatstart = wnet.flatstart, include_limits = include_limits, verbose = verbose)
    _set_rectangular_pf_status!(net, _apslf_pf_status(status == 0, sol.residual_inf, iters))
    has_merges && sync_merges_back!()
    return iters, status
  end

  try
    artifact_dir = performance_profile isa AbstractDict ? String(get(performance_profile, :output_dir, tempdir())) : tempdir()
    mkpath(artifact_dir)
    island_artifact = joinpath(artifact_dir, "ac_islands.csv")
    write_ac_island_report(island_artifact, island_report)
    println("AC island diagnostic artifact: ", island_artifact)
  catch err
    @warn "Unable to write AC island diagnostic artifact" exception = (err, catch_backtrace())
  end
  pf_cfg.islands_enabled || error(AC_ISLAND_DISABLED_MESSAGE)
  pf_cfg.islands_mode === :solve_independent || error("Unsupported power_flow.islands.mode=$(pf_cfg.islands_mode).")
  pf_cfg.islands_reference_policy === :matpower_like || error("Unsupported power_flow.islands.reference_policy=$(pf_cfg.islands_reference_policy).")
  _validate_island_references!(island_report)

  total_iters = 0
  first_failure = nothing
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
      if status != 0
        failure_status = _apslf_pf_status(false, sol.residual_inf, it; extra = (island_id = row.island_id,))
        island_statuses[Int(row.island_id)] = failure_status
        _set_rectangular_pf_status!(net, failure_status)
        err = ErrorException("AC island $(row.island_id) power-flow solve failed (solver=apslf):\n  buses=$(row.n_bus) branches=$(row.n_branch) ref=$(row.chosen_ref_bus)\n  bus_types: PV=$(row.n_pv) PQ=$(row.n_pq) REF=$(row.n_ref)\n  iterations=$(it) residual_inf=$(sol.residual_inf)")
        first_failure === nothing && (first_failure = err)
        pf_cfg.islands.diagnostic_continue_after_failure || throw(err)
        continue
      end
      all(isfinite(something(wnode._vm_pu, NaN)) && isfinite(something(wnode._va_deg, NaN)) for wnode in inet.nodeVec) || error("AC island $(row.island_id) produced nonfinite voltage results.")
      success_status = _apslf_pf_status(true, sol.residual_inf, it; extra = (island_id = row.island_id,))
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
  aggregate_status = _apslf_pf_status(true, aggregate_final_mismatch, total_iters; extra = (
    reason_text = "All AC islands converged independently.",
    island_wise_all_converged = true,
    stage = :island_wise_complete,
  ))
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

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
# This file is included inside module Sparlectra. Do not add a module wrapper here.
#
# DC power-flow per-island orchestration and Net write-back.
# Mirrors _run_apslf_powerflow! (src/acpflow/apslf_execution.jl)'s shape:
# detect_ac_islands -> per island _prepare_island_net/solve/_sync_island_solution!
# -> failure aggregation. Parallel-link merge handling (_merged_pf_net, used by
# the AC/APSLF paths) is out of scope for this first iteration.

# file: src/powerflow_dc/dc_network_solver.jl

# Minimal status NamedTuple for a DC solve. Deliberately does not reuse
# _apslf_pf_status/rectangular_pf_status's field shape: DC has no
# nr_converged/q_limit_active_set_ok/wrong_branch_* concepts.
function _dc_pf_status(converged::Bool, iterations::Int; extra::NamedTuple = NamedTuple())
  base = (
    numerical_converged = converged,
    final_converged = converged,
    status = converged ? :converged : :not_converged,
    reason = converged ? :none : :dc_solve_failed,
    reason_text = converged ? "DC power-flow linear solve completed." : "DC power-flow linear solve failed.",
    iterations = iterations,
    solver = :dc,
    stage = :dc_solve,
    exception_type = "",
    exception_message = "",
    stacktrace_top = "",
  )
  return merge(base, extra)
end

# Writes theta (deg) + vm_pu=1.0 into every bus, Pf/Pt (qFlow=0, lossless)
# into every in-service branch, and the computed slack injection, all on
# `inet` (an island-local net, or `net` itself in the single-island case).
# Populating the same four branch-flow fields the AC path uses is what lets
# _sync_island_solution! be reused completely unmodified for the multi-island
# case.
function _write_dc_solution!(inet::Net, slack_idx::Int, theta::Vector{Float64}, Pf::Vector{Float64}, Pt::Vector{Float64}, slack_p_mw::Float64, terms::Vector{_DcBranchTerm})
  @inbounds for (k, node) in enumerate(inet.nodeVec)
    node._vm_pu = 1.0
    node._va_deg = rad2deg(theta[k])
  end
  @inbounds for (l, term) in enumerate(terms)
    br = inet.branchVec[term.branch_idx]
    br.fBranchFlow = BranchFlow(1.0, rad2deg(theta[term.from]), Pf[l], 0.0)
    br.tBranchFlow = BranchFlow(1.0, rad2deg(theta[term.to]), Pt[l], 0.0)
    br.pLosses = 0.0
    br.qLosses = 0.0
  end
  inet.nodeVec[slack_idx]._pƩGen = slack_p_mw
  return nothing
end

_dc_island_local_slack(inet::Net)::Int = something(findfirst(n -> getNodeType(n) == Slack, inet.nodeVec), 0)

"""
    _run_dc_powerflow!(net::Net, pf_cfg::PowerFlowConfig; verbose=0, performance_profile=nothing) -> (iterations::Int, erg::Int)

Solve `net`'s DC power flow (`power_flow.solver = :dc`), respecting
`power_flow.islands`: each AC island (`detect_ac_islands`, reused unmodified)
is solved independently via [`solve_dc_powerflow`](@ref) with its own
`matpower_like` reference bus, and the result is written back via
[`_write_dc_solution!`](@ref) + `_sync_island_solution!` (reused unmodified).
Sets [`dc_pf_status`](@ref) on `net`. `iterations` is always `1` (a direct
linear solve, not an iterative process); `erg == 0` means every island
solved (a linear solve either succeeds or throws — there is no partial
convergence state to report).
"""
function _run_dc_powerflow!(net::Net, pf_cfg::PowerFlowConfig; verbose::Int = 0, performance_profile = nothing)::Tuple{Int,Int}
  refreshBusTypesFromProsumers!(net)
  angle_reference_rad = deg2rad(pf_cfg.dc.angle_reference_deg)
  island_report = detect_ac_islands(net)
  multi_island = length(island_report.rows) > 1 && any(row -> row.n_branch > 0, island_report.rows)

  if !multi_island
    slack_idx = _dc_island_local_slack(net)
    slack_idx == 0 && error("DC power flow requires a Slack/REF bus; none found.")
    theta, Pf, Pt, slack_p_mw, terms = solve_dc_powerflow(net, slack_idx; angle_reference_rad, performance_profile)
    _write_dc_solution!(net, slack_idx, theta, Pf, Pt, slack_p_mw, terms)
    _set_dc_pf_status!(net, _dc_pf_status(true, 1))
    verbose > 0 && println("DC power flow solved: 1 island, ", length(terms), " branches, slack bus ", slack_idx, ", slack injection = ", round(slack_p_mw; digits = 4), " MW")
    return 1, 0
  end

  pf_cfg.islands.enabled || error(AC_ISLAND_DISABLED_MESSAGE)
  pf_cfg.islands.mode === :solve_independent || error("Unsupported power_flow.islands.mode=$(pf_cfg.islands.mode).")
  pf_cfg.islands.reference_policy === :matpower_like || error("Unsupported power_flow.islands.reference_policy=$(pf_cfg.islands.reference_policy).")
  _validate_island_references!(island_report)

  first_failure = nothing
  island_statuses = Dict{Int,Any}()
  performance_profile isa AbstractDict && (performance_profile[:ac_island_solver_statuses] = island_statuses)

  for row in island_report.rows
    if first_failure !== nothing && !pf_cfg.islands.diagnostic_continue_after_failure
      skipped = _dc_pf_status(false, 0; extra = (island_id = row.island_id, reason = :skipped_after_previous_failure, status = :skipped_after_previous_failure, stage = :skipped_after_previous_failure, reason_text = "previous island failed and diagnostic continuation is disabled"))
      island_statuses[Int(row.island_id)] = skipped
      _set_dc_pf_status!(net, skipped)
      continue
    end
    try
      inet = _prepare_island_net(net, row)
      slack_idx = _dc_island_local_slack(inet)
      slack_idx == 0 && error("AC island $(row.island_id) has no Slack/REF bus after preparation.")
      theta, Pf, Pt, slack_p_mw, terms = solve_dc_powerflow(inet, slack_idx; angle_reference_rad, performance_profile)
      _write_dc_solution!(inet, slack_idx, theta, Pf, Pt, slack_p_mw, terms)
      success = _dc_pf_status(true, 1; extra = (island_id = row.island_id,))
      island_statuses[Int(row.island_id)] = success
      _sync_island_solution!(net, inet, row)
    catch err
      frames = stacktrace(catch_backtrace())
      top = isempty(frames) ? "" : sprint(show, first(frames))
      failure = _dc_pf_status(false, 0; extra = (island_id = row.island_id, reason = :solver_exception, status = :failed, exception_type = nameof(typeof(err)), exception_message = sprint(showerror, err), stacktrace_top = top))
      island_statuses[Int(row.island_id)] = failure
      _set_dc_pf_status!(net, failure)
      first_failure === nothing && (first_failure = ErrorException("AC island $(row.island_id) DC power-flow solve failed:\n  buses=$(row.n_bus) branches=$(row.n_branch) ref=$(row.chosen_ref_bus)\n  $(sprint(showerror, err))"))
      pf_cfg.islands.diagnostic_continue_after_failure || rethrow()
    end
  end

  if first_failure !== nothing
    island_message = _islandwise_failure_message(performance_profile)
    throw(ErrorException(island_message === nothing ? sprint(showerror, first_failure) : island_message))
  end

  aggregate = _dc_pf_status(true, 1; extra = (reason_text = "All AC islands solved independently.", island_wise_all_converged = true, stage = :island_wise_complete))
  _set_dc_pf_status!(net, aggregate)
  performance_profile isa AbstractDict && (performance_profile[:island_wise_all_converged] = true)
  verbose > 0 && println("DC power flow solved: ", length(island_report.rows), " islands, all solved independently")
  return 1, 0
end

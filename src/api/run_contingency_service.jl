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

# file: src/api/run_contingency_service.jl
# purpose: Web UI/service N-1 contingency run (issue #331 Phase 5): build the
#          net through the shared config-driven import path (MATPOWER) or
#          importCGMES (CGMES), run runContingencies! for the requested outage
#          kind (branch or generator), write contingency_n1.csv and a run.log
#          narrative from the ContingencyReport, and return a SparlectraApiResult
#          through the normal result conventions. Mirrors run_short_circuit_service.jl
#          (a mode flag on POST /powerflow/run, artifacts collected from the run
#          dir, no cache workflow of its own). rescue_ladder is read from the
#          config; the outage kind is a run parameter, not a config key.

"""
    _run_contingency_service(case_path, config_file, output_dir, run_id, kind; weights_path = nothing) -> SparlectraApiResult

Service backend of the Web UI "Contingency (N-1)" button. Builds the net
(MATPOWER through the shared config-driven import path, CGMES through
`importCGMES`), enumerates single-element outages for `kind` (`"branch"` or
`"gen"`), evaluates them with [`runContingencies!`](@ref) using
`config.contingency.rescue_ladder`, writes `contingency_n1.csv` plus a `run.log`
narrative built from [`buildContingencyReport`](@ref), and returns a
`SparlectraApiResult` through the normal result conventions.

When `weights_path` names an existing per-case weight file it is read with
[`readContingencyWeightsCSV`](@ref) and applied with
[`applyContingencyWeights`](@ref) after case generation; the presence of the
file IS the switch (there is no request key). Weights only reorder the
severity ranking, never skip a case. Names that match no case are warned about
in the `run.log` (a weight list can outlive a case edit), never fatal; a file
that cannot be read leaves the run unweighted rather than failing it. The
metadata reports `contingency_weights_applied` and `contingency_weighted_cases`.

A generator outage that removes the system's only slack is reported as a
non-converged case ("no slack bus registered"); that is the expected N-1 finding
that the unit is critical, named explicitly in the run.log and the run message
so it does not read as a tool failure. Rerun with `auto_slack = true` (a solver
keyword) to have the solver promote a surviving generator instead.

Failure behavior: `contingency_unsupported_format` (not MATPOWER or CGMES),
`contingency_no_cases` (no in-service element of the requested kind),
`invalid_request` (bad `kind`), plus the shared import/config failures.
"""
function _run_contingency_service(case_path::AbstractString, config_file::AbstractString, output_dir::AbstractString, run_id::String, kind::AbstractString; weights_path::Union{Nothing,AbstractString} = nothing)::SparlectraApiResult
  mkpath(output_dir)
  logfile = joinpath(output_dir, "run.log")
  result_file = joinpath(output_dir, "result.json")
  base_metadata = Dict{String,Any}("run_mode" => "contingency", "contingency_kind" => kind)

  if !(kind in ("branch", "gen"))
    return _api_failure("invalid_request", "contingency_kind must be \"branch\" or \"gen\", got \"$(kind)\".", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  config = try
    load_sparlectra_config(config_file; reload = true)
  catch err
    return _api_failure("invalid_configuration", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  # Build the net through the SAME paths the power-flow service uses: the shared
  # config-driven MATPOWER import (so pegase shift/ratio options are honored) or
  # importCGMES. runContingencies! solves the base case itself, so an unsolved
  # net is all we hand it.
  format = _detect_case_format(case_path)
  net = try
    if format === :matpower
      _import_sparlectra_net(case_path, nothing, config)
    elseif format === :cgmes
      cgmes_cfg = config.cgmes
      paths, _ = _cgmes_delivery_paths(case_path, cgmes_cfg)
      importCGMES(
        path = length(paths) == 1 ? paths[1] : paths,
        baseMVA = cgmes_cfg.base_mva,
        require_boundary = cgmes_cfg.require_boundary,
        tap_control = cgmes_cfg.tap_control,
        machine_control = cgmes_cfg.machine_control,
        ignore_connected = cgmes_cfg.ignore_connected,
        vset_min_pu = cgmes_cfg.vset_min_pu,
        vset_max_pu = cgmes_cfg.vset_max_pu,
        multi_slack = cgmes_cfg.multi_slack,
        strict_placeholder_guards = cgmes_cfg.placeholder_guards === :strict,
        infer_base_voltages = cgmes_cfg.infer_base_voltages,
        name = basename(case_path),
      ).net
    else
      return _api_failure("contingency_unsupported_format", "N-1 contingency needs a MATPOWER or CGMES case; got format $(format).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
  catch err
    err isa PowerFlowAborted && rethrow()
    return _api_failure("import_error", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  cases = kind == "gen" ? generateN1Generators(net) : generateN1Branches(net)
  if isempty(cases)
    return _api_failure("contingency_no_cases", "No in-service $(kind == "gen" ? "generator" : "branch") to take out for N-1.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  # Per-case weights (issue #331 Phase 5 follow-up): the file's PRESENCE next to
  # the case is the switch, there is no request key. A weight list may outlive a
  # case edit, so names matching no case are warned about, never fatal; a file
  # that cannot be read leaves the run UNWEIGHTED rather than failing it.
  weights_applied = false
  weighted_cases = 0
  unmatched = String[]
  weights_error = nothing
  if weights_path !== nothing && isfile(weights_path)
    try
      weights = readContingencyWeightsCSV(weights_path)
      case_names = Set(c.name for c in cases)
      unmatched = sort!([String(n) for n in keys(weights) if !(n in case_names)])
      cases = applyContingencyWeights(cases, weights)
      weighted_cases = count(c -> c.weight != 1.0, cases)
      weights_applied = true
    catch err
      weights_error = sprint(showerror, err)
    end
  end

  results = runContingencies!(net, cases; rescue_ladder = config.contingency.rescue_ladder)
  report = buildContingencyReport(results)
  csv = writeContingencyResultsCSV(joinpath(output_dir, "contingency_n1.csv"), results)

  # a slack-unit outage surfaces as "no slack bus registered"; name it so the
  # result page does not read it as a tool failure (see the docstring)
  n_no_slack = count(r -> r.error !== nothing && occursin("no slack bus", r.error), results)

  open(logfile, "a") do io
    println(io, "N-1 contingency (", kind == "gen" ? "generator" : "branch", " outages) on ", basename(case_path))
    println(io, "rescue ladder: ", config.contingency.rescue_ladder)
    if weights_applied
      println(io, "weights: applied ", weighted_cases, " non-default weight(s) from ", basename(String(weights_path)))
      if !isempty(unmatched)
        shown = first(unmatched, 10)
        println(io, "weights: ", length(unmatched), " weight name(s) match no ", kind == "gen" ? "generator" : "branch",
          " (a weight list can outlive a case edit): ", join(shown, ", "), length(unmatched) > 10 ? ", ..." : "")
      end
    elseif weights_error !== nothing
      println(io, "weights: a weight file was present but could not be read (", weights_error, "); the run is UNWEIGHTED")
    end
    println(io)
    printContingencyReport(io, report)
    if n_no_slack > 0
      println(io)
      println(io, n_no_slack, " outage(s) removed the system's only voltage reference and are reported as \"no slack bus registered\". That is the expected N-1 finding that the unit is critical, not a tool error; rerun with auto_slack = true to let the solver promote a surviving generator.")
    end
    println(io, "Artifacts: ", basename(csv))
  end

  metadata = merge(
    base_metadata,
    Dict{String,Any}(
      "input_format_detected" => String(format),
      "contingency_cases" => report.n_cases,
      "contingency_converged" => report.n_converged,
      "contingency_islanded" => report.n_islanded,
      "contingency_nonconverged" => report.n_nonconverged,
      "contingency_no_slack" => n_no_slack,
      "contingency_total_shed_mw" => report.total_shed_load_mw,
      "contingency_worst_loading_pct" => report.worst_loading_pct,
      "contingency_worst_severity" => report.worst_severity,
      # weights: surfaced so a weighted ranking is never silently mistaken for
      # an unweighted one
      "contingency_weights_applied" => weights_applied,
      "contingency_weighted_cases" => weighted_cases,
      # explicit run-status keys so history/result views render it like any
      # completed run
      "artifact_status" => "completed",
      "solver_status" => "completed",
      "service_status" => "completed",
      "run_status" => "completed",
    ),
  )

  message = string(
    "Contingency N-1 (", kind == "gen" ? "generator" : "branch", ") completed - ",
    report.n_converged, " of ", report.n_cases, " converged, ",
    report.n_islanded, " islanded (", round(report.total_shed_load_mw; digits = 1), " MW shed)",
    weights_applied ? ", $(weighted_cases) weighted" : "",
    n_no_slack > 0 ? ", $(n_no_slack) removed the only slack (see run.log)" : "",
    ".",
  )
  result = _api_result(
    run_id = run_id,
    status = :succeeded,
    success = true,
    solution_available = false,
    reason = nothing,
    message = message,
    casefile = String(case_path),
    config_file = String(config_file),
    output_dir = String(output_dir),
    logfile = logfile,
    result_file = result_file,
    metadata = metadata,
  )
  return _finalize_api_result(result)
end

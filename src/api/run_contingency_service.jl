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

Scenario task step 5: the request may carry `scenario_source` (`file_block`
runs the case file's own `scenarios`/`contingencies` block through the
scenario model, `external_file` a scenario JSON named by `scenario_file`,
`n1_all`/`n1_branches`/`n1_generators` the generated case lists) and
`screening_mode` (`off`/`flag`/`only`; absent means the
`contingency.screening` configuration, default `off`). With screening
active the CSV gains the `screened`/`screening_estimate` columns and the
metadata reports `contingency_screening_mode` and `contingency_screened`.

Failure behavior: `contingency_unsupported_format` (not MATPOWER or CGMES),
`contingency_no_cases` (no in-service element of the requested kind),
`invalid_request` (bad `kind`, unknown `scenario_source` or
`screening_mode`, missing `scenario_file`), plus the shared import/config
failures.
"""

# Resolve the case file's contingency definition into runner cases. Component
# ids are resolved through `extra[<id>].name`, which is why that name is
# mandatory wherever anything refers to an object. The runner takes ONE
# element per case, so a multi-outage (common-mode or N-2) entry is rejected
# by name instead of being silently reduced to its first element.
function _scf_contingency_cases_from_file(case_path::AbstractString, study::AbstractDict, net::Net)::Vector{ContingencyCase}
  names = scf_extra_names(case_path)
  name_of(id) = get(names, Int(id)) do
    throw(ArgumentError("SCF contingencies: component $(id) has no name in `extra`; the case list cannot be resolved."))
  end
  # The element kind follows from what the name resolves to in the network,
  # not from a second field in the file: a name that matches neither a branch
  # nor a generator is a broken case list and says so here.
  branch_names = Set(getCompName(br.comp) for br in net.branchVec)
  gen_cases = Dict(c.element => c for c in generateN1Generators(net))
  out = ContingencyCase[]
  for case in get(study, "cases", [])
    outages = get(case, "outages", [])
    length(outages) == 1 || throw(ArgumentError("SCF contingencies: case $(get(case, "name", "?")) lists $(length(outages)) outages; the N-1 runner takes one element per case (common-mode and N-2 studies are not supported yet)."))
    element = name_of(_scf_int(only(outages)["component"], "a contingency outage component"))
    kind = element in branch_names ? :branch : haskey(gen_cases, element) ? :gen : throw(ArgumentError("SCF contingencies: $(repr(element)) is neither an in-service branch nor a generator of this network."))
    name = String(get(case, "name", element))
    weight = haskey(case, "weight") ? Float64(case["weight"]) : 1.0
    push!(out, ContingencyCase(name, kind, element, weight))
  end
  return out
end

function _run_contingency_service(case_path::AbstractString, config_file::AbstractString, output_dir::AbstractString, run_id::String, kind::AbstractString; weights_path::Union{Nothing,AbstractString} = nothing, se_state_file::Union{Nothing,AbstractString} = nothing, se_run_id::Union{Nothing,AbstractString} = nothing, se_start_mode::AbstractString = "se_state", scenario_source::Union{Nothing,AbstractString} = nothing, scenario_file::Union{Nothing,AbstractString} = nothing, screening_mode::Union{Nothing,AbstractString} = nothing, screening_margin_pct::Union{Nothing,Real} = nothing)::SparlectraApiResult
  mkpath(output_dir)
  logfile = joinpath(output_dir, "run.log")
  result_file = joinpath(output_dir, "result.json")
  base_metadata = Dict{String,Any}("run_mode" => "contingency", "contingency_kind" => kind)
  if se_run_id !== nothing
    # SE-started base case (SE phase 5 chain): noted up front so even a
    # failure result names the referenced SE run
    base_metadata["se_run_id"] = String(se_run_id)
    base_metadata["se_start_mode"] = String(se_start_mode)
  end

  if !(kind in ("branch", "gen"))
    return _api_failure("invalid_request", "contingency_kind must be \"branch\" or \"gen\", got \"$(kind)\".", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  # scenario task step 5: the request may pick the scenario source; nothing
  # keeps the historical kind/contingencies behavior
  if scenario_source !== nothing && !(scenario_source in ("file_block", "external_file", "n1_all", "n1_branches", "n1_generators"))
    return _api_failure("invalid_request", "scenario_source must be one of file_block, external_file, n1_all, n1_branches, n1_generators; got \"$(scenario_source)\".", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  if scenario_source == "external_file" && !(scenario_file isa AbstractString && isfile(scenario_file))
    return _api_failure("invalid_request", "scenario_source external_file needs an existing scenario_file (a JSON carrying the scenarios block).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  # The same precedence the power-flow path uses (resolve_config, D5), minus
  # request overrides (a study run has no config form of its own): case
  # configuration file, the case file's deprecated block, general file,
  # defaults.
  config = try
    resolve_config(config_file, case_path).config
  catch err
    return _api_failure(_config_resolve_reason(err), sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  # Build the net through the SAME paths the power-flow service uses: the shared
  # config-driven MATPOWER import (so pegase shift/ratio options are honored) or
  # importCGMES. runContingencies! solves the base case itself, so an unsolved
  # net is all we hand it.
  format = _detect_case_format(case_path)
  format in (:matpower, :scf, :cgmes) || return _api_failure("contingency_unsupported_format", "N-1 contingency needs a MATPOWER, CGMES, or Sparlectra Case Format case; got format $(format).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  # id-addressed scenario sources need the typed case; a CGMES delivery
  # carries no referencable component ids, and the rejection NAMES the way
  # out (checked by a test on the message text). Fails fast, before the
  # import.
  if scenario_source in ("file_block", "external_file") && !(format in (:scf, :matpower))
    return _api_failure("invalid_request", "scenario_source $(scenario_source) resolves component ids against the typed case, which a CGMES delivery does not carry. Export the case as SCF once (the Web UI's \"Export as SCF case file\" button or exportSCF), then run the scenario source against that .scf.json and reference its component ids. The n1_all/n1_branches/n1_generators scenario sources work on the CGMES case directly.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  imported = try
    import_case(case_path, config; run_kind = :contingency)
  catch err
    err isa PowerFlowAborted && rethrow()
    return _api_failure("import_error", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  net = imported.net
  # step 3a: the run continues on the effective config of the import (CGMES
  # start-value decision, auto-profile rewrites), like the power-flow service
  config = imported.config
  start_decision = get(imported.provenance, "cgmes_start_decision", nothing)
  start_decision === nothing || open(logfile, "a") do io
    println(io, start_decision)
  end

  # SE-started base case (SE phase 5 chain): restore the estimated state
  # before the base solve. The service net is a throwaway import, so the
  # snapshot's balance takeover mutates nothing persistent. runContingencies!
  # then solves the base from the estimated voltages (flat start off).
  if se_state_file !== nothing
    try
      readSEStateCSV!(net; file = se_state_file)
      if se_start_mode == "se_snapshot"
        st = _se_start_state(net)
        for (i, n) in enumerate(net.nodeVec)
          n._pƩLoad = something(n._pƩGen, 0.0) - st.pinj[i]
          n._qƩLoad = something(n._qƩGen, 0.0) - st.qinj[i]
        end
      end
      net.flatstart = false
    catch err
      return _api_failure("se_start_failed", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
  end

  # scenario task step 5: screening mode from the request, falling back to
  # the contingency.screening configuration (whose default is :flag); the
  # margin always comes from the configuration
  screen_mode = screening_mode === nothing ? config.contingency.screening_mode : Symbol(lowercase(String(screening_mode)))
  if !(screen_mode in CONTINGENCY_SCREENING_MODE_VALUES)
    return _api_failure("invalid_request", "screening_mode must be off, flag, or only; got \"$(screening_mode)\".", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  screen_margin = screening_margin_pct === nothing ? config.contingency.screening_margin_pct : Float64(screening_margin_pct)
  if !(isfinite(screen_margin) && screen_margin >= 0.0)
    return _api_failure("invalid_request", "screening_margin_pct must be a finite value >= 0; got $(screen_margin).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  base_metadata["contingency_screening_mode"] = String(screen_mode)

  # Scenario source (step 5). file_block and external_file run through the
  # scenario model addressed by SCF component ids, so they need the typed
  # case as their ID INDEX (SCF directly, MATPOWER through the explicit
  # converter). task_import_direct: the converted case serves as the
  # addressing index ONLY; the electrics that run are the imported net,
  # so no conversion artifact reaches the results. The n1_* sources and
  # the historical kind/contingencies path stay on the generated case
  # lists, which every format supports.
  results = nothing
  scf_cases_source = "generated"
  cases = ContingencyCase[]
  if scenario_source in ("file_block", "external_file")
    scenario_result = try
      scfcase = format === :scf ? read_scf_json(case_path) : convert_case(MatpowerAdapter(), MatpowerIO.read_case(case_path; legacy_compat = true), matpower_adapter_options(config))
      idx = ScenarioIndex(scfcase)
      set = if scenario_source == "file_block"
        s = scf_case_scenarios(scfcase)
        s === nothing && throw(ArgumentError("scenario_source file_block: the case file carries neither a scenarios nor a contingencies block."))
        s
      else
        scenario_set_from_dict(scf_json_parse(read(String(scenario_file), String)))
      end
      runScenarios!(net, set; index = idx, rescue_ladder = config.contingency.rescue_ladder, screening_mode = screen_mode, screening_margin_pct = screen_margin)
    catch err
      err isa PowerFlowAborted && rethrow()
      return _api_failure("invalid_case_file", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
    results = scenario_result
    scf_cases_source = string("scenario_", scenario_source)
    base_metadata["contingency_cases_source"] = scf_cases_source
    if isempty(results)
      return _api_failure("contingency_no_cases", "The scenario source $(scenario_source) produced no scenario to run.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
  elseif scenario_source in ("n1_all", "n1_branches", "n1_generators")
    cases = scenario_source == "n1_branches" ? generateN1Branches(net) : scenario_source == "n1_generators" ? generateN1Generators(net) : vcat(generateN1Branches(net), generateN1Generators(net))
    scf_cases_source = String(scenario_source)
    base_metadata["contingency_cases_source"] = scf_cases_source
  else
    # A Sparlectra Case Format case may define the study itself (#342): the
    # file's `contingencies` block then decides which cases run, so shipping a
    # case ships the study it was meant for. Explicit request keys still pick
    # the element kind; the file only supplies the case list.
    scf_study = imported.studies.contingencies
    cases = kind == "gen" ? generateN1Generators(net) : generateN1Branches(net)
    if !isempty(scf_study)
      mode = String(get(scf_study, "mode", "explicit"))
      listed = try
        _scf_contingency_cases_from_file(case_path, scf_study, net)
      catch err
        return _api_failure("invalid_case_file", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
      end
      excluded = Set(String(x) for x in get(scf_study, "exclude", []))
      if mode == "explicit"
        cases = listed
        scf_cases_source = "case_file_explicit"
      elseif mode == "all_branches"
        cases = [c for c in cases if !(c.name in excluded) && !(c.element in excluded)]
        scf_cases_source = "case_file_all_branches"
      else
        generated = [c for c in cases if !(c.name in excluded) && !(c.element in excluded)]
        names = Set(c.name for c in generated)
        cases = vcat(generated, [c for c in listed if !(c.name in names)])
        scf_cases_source = "case_file_all_branches_plus"
      end
      base_metadata["contingency_cases_source"] = scf_cases_source
    end
  end
  if results === nothing && isempty(cases)
    return _api_failure("contingency_no_cases", scf_cases_source == "generated" ? "No in-service $(kind == "gen" ? "generator" : "branch") to take out for N-1." : "The case list source $(scf_cases_source) left no case to run.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  # Per-case weights (issue #331 Phase 5 follow-up): the file's PRESENCE next to
  # the case is the switch, there is no request key. A weight list may outlive a
  # case edit, so names matching no case are warned about, never fatal; a file
  # that cannot be read leaves the run UNWEIGHTED rather than failing it.
  weights_applied = false
  weighted_cases = 0
  unmatched = String[]
  weights_error = nothing
  if results === nothing && weights_path !== nothing && isfile(weights_path)
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

  if results === nothing
    results = runContingencies!(net, cases; rescue_ladder = config.contingency.rescue_ladder, screening_mode = screen_mode, screening_margin_pct = screen_margin)
  end
  n_screened = eltype(results) === ScenarioResult ? count(r -> r.screened, results) : 0
  report = buildContingencyReport(results)
  csv = writeContingencyResultsCSV(joinpath(output_dir, "contingency_n1.csv"), results)

  # a slack-unit outage surfaces as "no slack bus registered"; name it so the
  # result page does not read it as a tool failure (see the docstring)
  n_no_slack = count(r -> r.error !== nothing && occursin("no slack bus", r.error), results)

  open(logfile, "a") do io
    println(io, "N-1 contingency (", kind == "gen" ? "generator" : "branch", " outages) on ", basename(case_path))
    if startswith(scf_cases_source, "scenario_")
      println(io, "scenario source: ", replace(scf_cases_source, "scenario_" => ""), " (", length(results), " scenario(s))")
    elseif startswith(scf_cases_source, "case_file")
      println(io, "case list: from the case file's contingencies block (mode ", replace(scf_cases_source, "case_file_" => ""), ", ", length(cases), " case(s))")
    elseif scf_cases_source != "generated"
      println(io, "case list: ", scf_cases_source, " (", length(cases), " case(s))")
    end
    se_run_id !== nothing && println(io, "base case: SE-started (", se_start_mode, ") from SE run ", se_run_id)
    println(io, "rescue ladder: ", config.contingency.rescue_ladder)
    if screen_mode === :off
      println(io, "screening: off (every scenario fully solved)")
    else
      println(io, "screening: ", screen_mode, " (margin ", screen_margin, " %), ", n_screened, " of ", length(results), " scenario(s) screened; screened rows carry first-order estimates, not full solves")
    end
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
      # screening (step 5): the mode that ran and how many rows are estimates
      "contingency_screened" => n_screened,
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
    n_screened > 0 ? ", $(n_screened) screened" : "",
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

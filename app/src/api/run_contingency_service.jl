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

# power mode (0.30.1) with its sparse LU choice (0.30.2) for the worker
# solves of the scenario engine and N-1, only when the configuration
# switches power mode on
_power_mode_kwargs(pf) = pf.power_mode ? (; power_mode = true, power_mode_lu = pf.power_mode_lu) : (;)

# dishonest Newton (0.30.2) for the worker solves of the scenario engine and
# N-1, only when the configuration switches it on
_jacobian_reuse_kwargs(pf) = pf.jacobian_reuse ? (; jacobian_reuse = true, jacobian_reuse_min_reduction = pf.jacobian_reuse_min_reduction, jacobian_reuse_max_steps = pf.jacobian_reuse_max_steps) : (;)

# the distributed-slack model of the run configuration for the scenario
# engine and N-1: every key of power_flow.distributed_slack, so the base
# case, its rescue and every scenario share the units as configured (the
# scenario-source path passed none of them and the N-1 path left out the
# weights and respect_p_limits, #456)
_distributed_slack_kwargs(pf) = pf.distributed_slack.enabled ? (; distributed_slack_enabled = true, distributed_slack_p_mode = pf.distributed_slack.p_mode, distributed_slack_respect_p_limits = pf.distributed_slack.respect_p_limits, distributed_slack_fallback = pf.distributed_slack.fallback, distributed_slack_weights = pf.distributed_slack.weights) : (;)

# the solver keywords both engine entries (scenario source, generated N-1
# cases) receive from the run configuration
# the Q-limit keywords come from the library's one mapping (the same the
# single run uses), so N-1 and scenarios honour power_flow.qlimits.*
_engine_solver_kwargs(pf) = (; _distributed_slack_kwargs(pf)..., _power_mode_kwargs(pf)..., _jacobian_reuse_kwargs(pf)..., _qlimit_solver_kwargs(pf.qlimits)...)

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
  # a case that lists exactly the legs of one three-winding transformer is
  # that transformer as a whole (it trips as one element)
  star_names = Dict{Int,String}(idx => name for (name, idx) in net.busDict)
  leg_sets = Dict(Set(getCompName(net.branchVec[k].comp) for k in g.legs) => g for g in _three_winding_groups(net))
  out = ContingencyCase[]
  for case in get(study, "cases", [])
    outages = get(case, "outages", [])
    if length(outages) > 1
      legs = Set(name_of(_scf_int(o["component"], "a contingency outage component")) for o in outages)
      grp = get(leg_sets, legs, nothing)
      if grp !== nothing
        push!(out, ContingencyCase(String(get(case, "name", grp.name)), :transformer3w, star_names[grp.star], haskey(case, "weight") ? Float64(case["weight"]) : 1.0))
        continue
      end
    end
    length(outages) == 1 || throw(ArgumentError("SCF contingencies: case $(get(case, "name", "?")) lists $(length(outages)) outages; the N-1 runner takes one element per case (common-mode and N-2 studies are not supported yet)."))
    element = name_of(_scf_int(only(outages)["component"], "a contingency outage component"))
    kind = element in branch_names ? :branch : haskey(gen_cases, element) ? :gen : throw(ArgumentError("SCF contingencies: $(repr(element)) is neither an in-service branch nor a generator of this network."))
    name = String(get(case, "name", element))
    weight = haskey(case, "weight") ? Float64(case["weight"]) : 1.0
    push!(out, ContingencyCase(name, kind, element, weight))
  end
  return out
end

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

An outage that removes the reference (the slack unit itself, or the branch
that ties it to an island) does not end the case: the run solves with
`auto_slack`, the best remaining unit takes over (stated reference priority
first, then the strongest unit; an island without a voltage-controlled unit
takes its best generating unit), and the row
names the bus ("reference taken over by bus ..."). Only an island without
any generating unit stays without a reference and is reported as islanded
with the load it loses.

Scenario task step 5: the request may carry `scenario_source` (`file_block`
runs the case file's own `scenarios`/`contingencies` block through the
scenario model, `external_file` a scenario JSON named by `scenario_file`,
`n1_all`/`n1_branches`/`n1_generators` the generated case lists) and
`screening_mode` (`off`/`flag`/`only`; absent means the
`contingency.screening` configuration, default `off`). With screening
active the CSV gains the `screened`/`screening_estimate` columns and the
metadata reports `contingency_screening_mode` and `contingency_screened`.

`progress` (default `nothing`) is handed to [`runContingencies!`](@ref) /
`runScenarios!` unchanged: a thread-safe callable `progress(done, total)`
called once per finished outage or scenario. The Web UI job passes one
that stores the counter in the job snapshot (see `webui_jobs.jl`).

Failure behavior: `contingency_unsupported_format` (not MATPOWER, CGMES,
PowSyBl or Sparlectra Case Format),
`contingency_no_cases` (no in-service element of the requested kind),
`invalid_request` (bad `kind`, unknown `scenario_source` or
`screening_mode`, missing `scenario_file`), plus the shared import/config
failures.
"""
function _run_contingency_service(case_path::AbstractString, config_file::AbstractString, output_dir::AbstractString, run_id::String, kind::AbstractString; config_overrides::AbstractDict = Dict{String,Any}(), weights_path::Union{Nothing,AbstractString} = nothing, se_state_file::Union{Nothing,AbstractString} = nothing, se_run_id::Union{Nothing,AbstractString} = nothing, se_start_mode::AbstractString = "se_state", scenario_source::Union{Nothing,AbstractString} = nothing, scenario_file::Union{Nothing,AbstractString} = nothing, screening_mode::Union{Nothing,AbstractString} = nothing, screening_margin_pct::Union{Nothing,Real} = nothing, @nospecialize(progress = nothing))::SparlectraApiResult
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
  # the request may pick the scenario source; nothing
  # keeps the historical kind/contingencies behavior
  if scenario_source !== nothing && !(scenario_source in ("file_block", "external_file", "n1_all", "n1_branches", "n1_generators"))
    return _api_failure("invalid_request", "scenario_source must be one of file_block, external_file, n1_all, n1_branches, n1_generators; got \"$(scenario_source)\".", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  if scenario_source == "external_file" && !(scenario_file isa AbstractString && isfile(scenario_file))
    return _api_failure("invalid_request", "scenario_source external_file needs an existing scenario_file (a JSON carrying the scenarios block).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  # The same precedence the power-flow path uses (resolve_config), minus
  # request overrides (a study run has no config form of its own): case
  # configuration file, the case file's deprecated block, general file,
  # defaults.
  config = try
    resolve_config(config_file, case_path, config_overrides).config
  catch err
    return _api_failure(_config_resolve_reason(err), sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  # Build the net through the SAME paths the power-flow service uses: the shared
  # config-driven MATPOWER import (so pegase shift/ratio options are honored) or
  # importCGMES. runContingencies! solves the base case itself, so an unsolved
  # net is all we hand it.
  format = _detect_case_format(case_path)
  format in (:matpower, :scf, :cgmes, :powsybl) || return _api_failure("contingency_unsupported_format", "N-1 contingency needs a MATPOWER, CGMES, PowSyBl or Sparlectra Case Format case; got format $(format).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  # id-addressed scenario sources need the typed case; a CGMES delivery
  # carries no referencable component ids, and the rejection NAMES the way
  # out (checked by a test on the message text). Fails fast, before the
  # import.
  if scenario_source in ("file_block", "external_file") && !(format in (:scf, :matpower))
    return _api_failure("invalid_request", "scenario_source $(scenario_source) resolves component ids against the typed case, which a CGMES delivery or a PowSyBl case does not carry. Export the case as SCF once (the Web UI's \"Export as SCF case file\" button or exportSCF), then run the scenario source against that .scf.json and reference its component ids. The n1_all/n1_branches/n1_generators scenario sources work on the CGMES case directly.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
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

  # screening mode from the request, falling back to
  # the contingency.screening configuration (whose default is :flag); the
  # margin always comes from the configuration
  screen_mode = screening_mode === nothing ? config.contingency.screening_mode : Symbol(lowercase(String(screening_mode)))
  if !(screen_mode in CONTINGENCY_SCREENING_MODE_VALUES)
    return _api_failure("invalid_request", "screening_mode must be off, flag, or only; got \"$(screening_mode)\".", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  screen_margin = screening_margin_pct === nothing ? config.contingency.screening_margin_pct : Float64(screening_margin_pct)
  # the warm-active-set line of the scenario engine, written into run.log
  warm_note = Ref("")
  if !(isfinite(screen_margin) && screen_margin >= 0.0)
    return _api_failure("invalid_request", "screening_margin_pct must be a finite value >= 0; got $(screen_margin).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  base_metadata["contingency_screening_mode"] = String(screen_mode)

  # Scenario source (step 5). file_block and external_file run through the
  # scenario model addressed by SCF component ids, so they need the typed
  # case as their ID INDEX (SCF directly, MATPOWER through the explicit
  # converter). The converted case serves as the
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
      # distributed slack, power mode and Jacobian reuse reach the scenario
      # engine exactly as on the N-1 path below (distributed slack used to be
      # missing here while run.log announced it, #456)
      runScenarios!(net, set; index = idx, rescue_ladder = config.contingency.rescue_ladder, maxIte = config.contingency.max_iter, base_maxIte = config.powerflow.max_iter, screening_mode = screen_mode, screening_margin_pct = screen_margin,
        _engine_solver_kwargs(config.powerflow)...,
        warm_active_set = config.contingency.warm_active_set, warm_cold_check = config.contingency.warm_cold_check, warm_cold_check_margin_pu = config.contingency.warm_cold_check_margin_pu, warm_note = warm_note, progress = progress)
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

  # the slack model of the run configuration applies to every post-outage
  # solve: with power_flow.distributed_slack on, the units share the
  # mismatch of an outage as they share it in the base case (the outage of
  # the line that carries the reference unit's output has no solution on a
  # single slack in the IEEE 14-bus case, and one with the shared slack)
  # (power mode: the worker nets of the engine keep their Ybus,
  # factorization and work arrays across the outages when the caller set it)
  dslack = config.powerflow.distributed_slack
  dslack_kwargs = _engine_solver_kwargs(config.powerflow)
  if results === nothing
    # an outage that removes the reference, or splits off an island without
    # one, does not end the case: the strongest remaining unit takes over
    # and the result names it (auto_slack, the default of runContingencies!)
    results = runContingencies!(net, cases; rescue_ladder = config.contingency.rescue_ladder, maxIte = config.contingency.max_iter, base_maxIte = config.powerflow.max_iter, screening_mode = screen_mode, screening_margin_pct = screen_margin, warm_active_set = config.contingency.warm_active_set, warm_cold_check = config.contingency.warm_cold_check, warm_cold_check_margin_pu = config.contingency.warm_cold_check_margin_pu, warm_note = warm_note, progress = progress, dslack_kwargs...)
  end
  n_screened = eltype(results) === ScenarioResult ? count(r -> r.screened, results) : 0
  # a requested screening the network does not allow is named in run.log
  screen_reason = screen_mode === :off ? nothing : _screening_unavailable_reason(net)
  report = buildContingencyReport(results)
  csv = writeContingencyResultsCSV(joinpath(output_dir, "contingency_n1.csv"), results; format = String(config.output.csv_format))

  # a slack-unit outage surfaces as "no slack bus registered"; name it so the
  # result page does not read it as a tool failure (see the docstring)
  n_no_slack = count(r -> r.error !== nothing && occursin("no slack bus", r.error), results)
  # cases that solved on a reference another unit took over
  n_reference_taken = count(r -> r.error !== nothing && occursin("reference taken over", r.error), results)

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
    # warm active set (0.30.2): one line when on, nothing when off
    isempty(warm_note[]) || println(io, warm_note[])
    println(io, "reference: an outage that removes it hands it to the strongest remaining unit (auto_slack)")
    println(io, "slack: ", dslack.enabled ? "distributed (power_flow.distributed_slack, p_mode $(dslack.p_mode))" : "single reference bus per island (power_flow.distributed_slack is off)")
    if screen_mode === :off
      println(io, "screening: off (every scenario fully solved)")
    elseif screen_reason !== nothing
      println(io, "screening: ", screen_mode, " requested but not available: ", screen_reason, "; every scenario fully solved")
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
      println(io, n_no_slack, " outage(s) removed the system's only voltage reference and no generating unit was left to take over (\"no slack bus registered\").")
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
      "contingency_reference_taken_over" => n_reference_taken,
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
    n_reference_taken > 0 ? ", $(n_reference_taken) solved on a reference another unit took over" : "",
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

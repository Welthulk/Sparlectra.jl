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
runs the case file's own `scenarios` block through the scenario model, `external_file` a scenario JSON named by `scenario_file`,
`n1_all`/`n1_branches`/`n1_generators` the generated case lists) and
`screening_mode` (`off`/`flag`/`only`; absent means the
`contingency.screening` configuration, default `off`). With screening
active the CSV gains the `screened`/`screening_estimate` columns and the
metadata reports `contingency_screening_mode` and `contingency_screened`.

`progress` (default `nothing`) is handed to [`runContingencies!`](@ref) /
`runScenarios!` unchanged: a thread-safe callable `progress(done, total)`
called once per finished outage or scenario. The Web UI job passes one
that stores the counter in the job snapshot (see `webui_jobs.jl`).

`phase_callback` (default no-op) receives the run phase by name:
`importing_case`, `solving_base_case`, `running_outages` (or
`running_scenarios` for a file block or external scenario file) and
`writing_artifacts`. The base case is solved inside the batch call, so
the switch to the batch phase comes with the first finished outage.

Failure behavior: `contingency_unsupported_format` (not MATPOWER, CGMES,
PowSyBl or Sparlectra Case Format),
`contingency_no_cases` (no in-service element of the requested kind),
`invalid_request` (bad `kind`, unknown `scenario_source` or
`screening_mode`, missing `scenario_file`), plus the shared import/config
failures.
"""
function _run_contingency_service(case_path::AbstractString, config_file::AbstractString, output_dir::AbstractString, run_id::String, kind::AbstractString; config_overrides::AbstractDict = Dict{String,Any}(), weights_path::Union{Nothing,AbstractString} = nothing, se_state_file::Union{Nothing,AbstractString} = nothing, se_run_id::Union{Nothing,AbstractString} = nothing, se_start_mode::AbstractString = "se_state", scenario_source::Union{Nothing,AbstractString} = nothing, scenario_file::Union{Nothing,AbstractString} = nothing, screening_mode::Union{Nothing,AbstractString} = nothing, screening_margin_pct::Union{Nothing,Real} = nothing, @nospecialize(progress = nothing), @nospecialize(phase_callback = phase -> nothing))::SparlectraApiResult
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
  phase_callback("importing_case")
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
  _record_network_frequency!(base_metadata, net, logfile)

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

  # screening mode and margin from the request, each falling back to the
  # contingency.screening configuration (mode default :off)
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

  # Display names of the output (0.40.0): an SCF case names its components
  # in the file, every other format keeps the generated names. The net's own
  # component names stay untouched. `generated_names` serves the weight
  # aliases below (weights saved under the generated names before 0.40.0).
  typed_case = nothing
  if format === :scf
    typed_case = try
      read_scf_json(case_path)
    catch err
      err isa PowerFlowAborted && rethrow()
      return _api_failure("invalid_case_file", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
  end
  generated_names = contingency_names(net)
  names = typed_case === nothing ? generated_names : contingency_names(net, typed_case)

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
  # The base case is solved inside runContingencies!/runScenarios! before
  # the batch starts, and the service cannot see that boundary. The first
  # finished outage marks it: the wrapper reports the batch phase once,
  # then forwards every tick to the caller's counter. It runs on the batch
  # worker threads; the lock keeps the phase report (the Web UI job writes
  # its dictionary there) from overlapping a counter update of another
  # worker.
  batch_phase = scenario_source in ("file_block", "external_file") ? "running_scenarios" : "running_outages"
  batch_started = Ref(false)
  batch_lock = ReentrantLock()
  batch_progress = (done, total) -> lock(batch_lock) do
    batch_started[] || (batch_started[] = true; phase_callback(batch_phase))
    progress === nothing || progress(done, total)
    return nothing
  end
  if scenario_source in ("file_block", "external_file")
    phase_callback("solving_base_case")
    scenario_result = try
      scfcase = format === :scf ? typed_case : convert_case(MatpowerAdapter(), MatpowerIO.read_case(case_path; legacy_compat = true), matpower_adapter_options(config))
      idx = ScenarioIndex(scfcase)
      set = if scenario_source == "file_block"
        s = scf_case_scenarios(scfcase)
        s === nothing && throw(ArgumentError("scenario_source file_block: the case file carries no scenarios block."))
        s
      else
        scenario_set_from_dict(scf_json_parse(read(String(scenario_file), String)))
      end
      # distributed slack, power mode and Jacobian reuse reach the scenario
      # engine exactly as on the N-1 path below (distributed slack used to be
      # missing here while run.log announced it, #456)
      runScenarios!(net, set; index = idx, rescue_ladder = config.contingency.rescue_ladder, maxIte = config.contingency.max_iter, base_maxIte = config.powerflow.max_iter, screening_mode = screen_mode, screening_margin_pct = screen_margin,
        _engine_solver_kwargs(config.powerflow)...,
        warm_active_set = config.contingency.warm_active_set, warm_cold_check = config.contingency.warm_cold_check, warm_cold_check_margin_pu = config.contingency.warm_cold_check_margin_pu, warm_note = warm_note,
        events = config.contingency.events, pairs = config.contingency.pairs, progress = batch_progress)
    catch err
      err isa PowerFlowAborted && rethrow()
      # the pair cap of a file block that asks for pairs: named as such
      reason = err isa ArgumentError && occursin("max_pairs", err.msg) ? "contingency_too_many_pairs" : "invalid_case_file"
      return _api_failure(reason, sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
    results = scenario_result
    scf_cases_source = string("scenario_", scenario_source)
    base_metadata["contingency_cases_source"] = scf_cases_source
    if isempty(results)
      return _api_failure("contingency_no_cases", "The scenario source $(scenario_source) produced no scenario to run.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
  elseif scenario_source in ("n1_all", "n1_branches", "n1_generators")
    cases = scenario_source == "n1_branches" ? generateN1Branches(net; names = names) : scenario_source == "n1_generators" ? generateN1Generators(net; names = names) : vcat(generateN1Branches(net; names = names), generateN1Generators(net; names = names))
    scf_cases_source = String(scenario_source)
    base_metadata["contingency_cases_source"] = scf_cases_source
  else
    cases = kind == "gen" ? generateN1Generators(net; names = names) : generateN1Branches(net; names = names)
  end

  # Correlated events and N-2 pairs (#355, contingency.events.* and
  # contingency.pairs.*): added to a generated BRANCH list of an SCF case.
  # A generator-only list and the other formats run plain N-1, and run.log
  # says why. A file block asks for its own events (its `events` key).
  ev_cfg = config.contingency.events
  pair_cfg = config.contingency.pairs
  event_note = nothing
  if results === nothing && (ev_cfg.enabled || pair_cfg.enabled)
    branch_list = scenario_source in ("n1_all", "n1_branches") || (scenario_source === nothing && kind == "branch")
    if typed_case === nothing
      event_note = "correlated events and pairs: not run, they need a Sparlectra Case Format case (this case is $(format))"
    elseif !branch_list
      event_note = "correlated events and pairs: not run, the case list holds generators only"
    else
      extra = try
        generateCorrelatedEvents(net; names = names,
          double_circuits = ev_cfg.enabled && ev_cfg.double_circuits, busbars = ev_cfg.enabled && ev_cfg.busbars,
          substations = ev_cfg.enabled && ev_cfg.substations, links = ev_cfg.enabled && ev_cfg.links, name_groups = ev_cfg.enabled && ev_cfg.name_groups,
          name_stem_pattern = ev_cfg.name_stem_pattern, max_members = ev_cfg.max_members,
          pairs = pair_cfg.enabled, pair_rule = pair_cfg.rule, max_pairs = pair_cfg.max_pairs)
      catch err
        (err isa ArgumentError && occursin("max_pairs", err.msg)) || rethrow()
        return _api_failure("contingency_too_many_pairs", err.msg; run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
      end
      # an event equal to an N-1 case (a one-branch station) is not repeated
      single = Set(c.name for c in cases)
      append!(cases, [c for c in extra if !(c.name in single)])
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
  alias_hits = 0
  weighted_cases = 0
  unmatched = String[]
  weights_error = nothing
  if results === nothing && weights_path !== nothing && isfile(weights_path)
    try
      weights = readContingencyWeightsCSV(weights_path)
      # weights saved under the generated names (before the file names of
      # an SCF case were used) still apply, as aliases with one log line
      aliases = _contingency_weight_aliases(net, names, generated_names, cases)
      case_names = Set(c.name for c in cases)
      alias_names = Set(values(aliases))
      unmatched = sort!([String(n) for n in keys(weights) if !(n in case_names) && !(n in alias_names)])
      alias_hits = count(c -> !haskey(weights, c.name) && haskey(weights, get(aliases, c.name, "")), cases)
      cases = applyContingencyWeights(cases, weights; aliases = aliases)
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
    phase_callback("solving_base_case")
    # an outage that removes the reference, or splits off an island without
    # one, does not end the case: the strongest remaining unit takes over
    # and the result names it (auto_slack, the default of runContingencies!)
    results = runContingencies!(net, cases; names = names, rescue_ladder = config.contingency.rescue_ladder, maxIte = config.contingency.max_iter, base_maxIte = config.powerflow.max_iter, screening_mode = screen_mode, screening_margin_pct = screen_margin, warm_active_set = config.contingency.warm_active_set, warm_cold_check = config.contingency.warm_cold_check, warm_cold_check_margin_pu = config.contingency.warm_cold_check_margin_pu, warm_note = warm_note, progress = batch_progress, dslack_kwargs...)
  end
  n_screened = eltype(results) === ScenarioResult ? count(r -> r.screened, results) : 0
  # a requested screening the network does not allow is named in run.log
  screen_reason = screen_mode === :off ? nothing : _screening_unavailable_reason(net)
  phase_callback("writing_artifacts")
  report = buildContingencyReport(results)
  csv = writeContingencyResultsCSV(joinpath(output_dir, "contingency_n1.csv"), results; format = String(config.output.csv_format))

  # a slack-unit outage surfaces as "no slack bus registered"; name it so the
  # result page does not read it as a tool failure (see the docstring)
  n_no_slack = count(r -> r.error !== nothing && occursin("no slack bus", r.error), results)
  # cases that solved on a reference another unit took over
  n_reference_taken = count(r -> r.error !== nothing && occursin("reference taken over", r.error), results)

  open(logfile, "a") do io
    if scenario_source == "file_block"
      println(io, "Scenarios of the case file on ", basename(case_path))
    else
      println(io, "N-1 contingency (", kind == "gen" ? "generator" : "branch", " outages) on ", basename(case_path))
    end
    if startswith(scf_cases_source, "scenario_")
      println(io, "scenario source: ", replace(scf_cases_source, "scenario_" => ""), " (", length(results), " scenario(s))")
    elseif scf_cases_source != "generated"
      println(io, "case list: ", scf_cases_source, " (", length(cases), " case(s))")
    end
    se_run_id !== nothing && println(io, "base case: SE-started (", se_start_mode, ") from SE run ", se_run_id)
    event_note === nothing || println(io, event_note)
    n_events = count(r -> r.kind in CONTINGENCY_EVENT_KINDS || r.kind === :pair, results)
    if n_events > 0
      by_kind = [(k, count(r -> r.kind === k, results)) for k in (CONTINGENCY_EVENT_KINDS..., :pair)]
      println(io, "correlated events: ", join(["$(n) $(k)" for (k, n) in by_kind if n > 0], ", "), " (", n_events, " of ", length(results), " case(s))")
    end
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
      alias_hits > 0 && println(io, "weights: ", alias_hits, " weight(s) stored under the generated names applied to the cases now named after the case file; saving the weights editor writes the new names")
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
      # cases per kind (#355): single outages, correlated events, pairs
      "contingency_cases_by_kind" => Dict{String,Any}(String(k) => n for (k, n) in report.n_by_kind),
      # explicit run-status keys so history/result views render it like any
      # completed run
      "artifact_status" => "completed",
      "solver_status" => "completed",
      "service_status" => "completed",
      "run_status" => "completed",
    ),
  )

  message = string(
    scenario_source == "file_block" ? "Scenarios of the case file completed - " : string("Contingency N-1 (", kind == "gen" ? "generator" : "branch", ") completed - "),
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

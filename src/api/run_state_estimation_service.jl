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

# file: src/api/run_state_estimation_service.jl
# purpose: Web UI/service state-estimation run (SE phase 5): build the net
#          through the shared import paths, read a measurement CSV v1 set,
#          check observability, run runse_diagnostics plus the final
#          runse!(updateNet = true), and persist the SE artifacts including
#          the se_state.csv chain anchor. Also the SE-started power-flow
#          chain run (_run_pf_from_se_service). Mirrors
#          run_contingency_service.jl (mode flag on POST /powerflow/run,
#          artifacts collected from the run dir).

"""
    _se_import_case(case_path, config) -> ImportedCase

Import a case for the state-estimation workflow through the one import
entry point. Throws on unsupported formats and on import errors. Shared by
the SE service run and the Web UI measurement generator so both accept
exactly the same cases, and both run with `ImportedCase.config` (the
CGMES start-value decision and the auto-profile rewrites reach the SE
solves exactly like the power-flow service, step 3a of the adapter task).
"""
function _se_import_case(case_path::AbstractString, config; requested_format::Symbol = :auto)::ImportedCase
  format = _detect_case_format(String(case_path); requested = requested_format)
  # format policy stays here (the wording is this service's contract)
  format in (:scf, :matpower, :cgmes, :dtf_for001) || throw(ArgumentError("State estimation needs a MATPOWER, CGMES, Sparlectra Case Format or DTF case; got format $(format)."))
  return import_case(String(case_path), config; requested_format = requested_format, run_kind = :state_estimation)
end

_se_import_case_net(case_path::AbstractString, config; requested_format::Symbol = :auto)::Net = _se_import_case(case_path, config; requested_format).net

## Format for decisions that only ask "which kind of case is this" (bus
## reference in the measurement CSV, CGMES specifics). A bare `.DAT` is
## ambiguous without an explicit request (FOR001 network vs FOR002
## reference), and that ambiguity must not abort a run that does not even
## depend on the answer, so it resolves to :unknown here.
function _se_case_format(case_path::AbstractString; requested::Symbol = :auto)::Symbol
  return try
    _detect_case_format(String(case_path); requested = requested)
  catch
    :unknown
  end
end

## A phase tap changer counts as DECLARED nameplate data (not a constructor
## default) when the import stored its regulating-vector direction psi in
## tap_est_alpha_deg (mpc.sparlectra.tap_changers, optional psi column,
## default 90). Only declared phase changers join the mass release and the
## generator's deviation targeting.
_declared_phase_tap(br)::Bool = br.has_phase_tap && (br.phase_step_deg > 0.0 || br.phase_du_step > 0.0) && br.tap_est_alpha_deg != 0.0

## from-run truth resolution (measurement generator v2): resolve `run_id` in
## the persistent run index under `run_root`, enforce the case binding, and
## load the solved voltages into `net` WITHOUT re-solving. SE runs read
## se_state.csv (bus by name); PF runs read the detail CSV
## bus_voltages_complex.csv, which exists only when the run wrote detail
## CSVs. A missing artifact rejects up front with the hint to repeat the run
## with the detail CSV enabled; there is deliberately no silent re-solve.
function _se_truth_from_run!(net::Net, run_root::AbstractString, run_id::AbstractString, case_path::AbstractString)
  index = load_powerflow_run_index(run_root)
  entry = nothing
  for e in get(index, "runs", Any[])
    e isa AbstractDict || continue
    if String(get(e, "run_id", "")) == run_id
      entry = e
      break
    end
  end
  entry === nothing && throw(ArgumentError("run $(run_id) not found in the run history"))
  get(entry, "success", false) == true || throw(ArgumentError("run $(run_id) was not successful; the generator only accepts converged runs"))
  bound = basename(String(get(entry, "casefile", "")))
  bound == basename(case_path) || throw(ArgumentError("run $(run_id) belongs to case $(isempty(bound) ? "?" : bound), not $(basename(case_path))"))
  run_mode = String(get(entry, "run_mode", ""))
  run_mode in ("", "se") || throw(ArgumentError("run $(run_id) is a $(run_mode) run; only power-flow and state-estimation runs provide a truth state"))
  outdir = String(get(entry, "output_dir", joinpath(run_root, run_id)))
  ts = String(get(entry, "timestamp", ""))
  if run_mode == "se"
    f = joinpath(outdir, "se_state.csv")
    isfile(f) || throw(ArgumentError("run $(run_id) has no se_state.csv artifact; repeat the state estimation"))
    readSEStateCSV!(net; file = f)
    return (kind = "se", timestamp = ts, source_file = "se_state.csv")
  end
  f = joinpath(outdir, "bus_voltages_complex.csv")
  isfile(f) || throw(ArgumentError("run $(run_id) has no bus_voltages_complex.csv artifact; repeat the run with the detailed result CSV enabled"))
  _se_apply_pf_voltage_csv!(net, f)
  return (kind = "pf", timestamp = ts, source_file = "bus_voltages_complex.csv")
end

## tolerant reader for the PF detail CSV bus_voltages_complex.csv: the file
## is written in a user-selected format (technical/excel_us: ',' delimiter,
## '.' decimals; excel_de: ';' delimiter, ',' decimals). Buses match by
## bus_name first, then by the original bus number; any shape or matching
## failure aborts with a line-precise error instead of guessing.
function _se_apply_pf_voltage_csv!(net::Net, file::AbstractString)
  lines = readlines(file)
  header_ln = findfirst(l -> !startswith(strip(l), "#") && !isempty(strip(l)), lines)
  header_ln === nothing && throw(ArgumentError("$(file): empty file"))
  delim = occursin(";", lines[header_ln]) ? ";" : ","
  cols = [strip(String(c)) for c in split(lines[header_ln], delim)]
  bi = findfirst(==("bus"), cols)
  ni = findfirst(==("bus_name"), cols)
  vi = findfirst(==("vm_pu"), cols)
  ai = findfirst(==("va_deg"), cols)
  (bi === nothing || vi === nothing || ai === nothing) && throw(ArgumentError("$(file): missing bus/vm_pu/va_deg columns"))
  num(s) = delim == ";" ? tryparse(Float64, replace(strip(s), "," => ".")) : tryparse(Float64, strip(s))
  orig_to_idx = Dict{Int,Int}(o => i for (i, o) in net.busOrigIdxDict)
  n = length(net.nodeVec)
  seen = falses(n)
  for ln = (header_ln+1):length(lines)
    line = strip(lines[ln])
    isempty(line) && continue
    fields = split(line, delim)
    length(fields) == length(cols) || throw(ArgumentError("$(file):$(ln): expected $(length(cols)) fields, got $(length(fields)) (a name containing the delimiter?)"))
    idx = nothing
    if ni !== nothing
      name = String(strip(fields[ni]))
      haskey(net.busDict, name) && (idx = net.busDict[name])
    end
    if idx === nothing
      bno = tryparse(Int, strip(fields[bi]))
      bno !== nothing && (idx = get(orig_to_idx, bno, (1 <= bno <= n && isempty(orig_to_idx)) ? bno : nothing))
    end
    idx === nothing && throw(ArgumentError("$(file):$(ln): bus '$(strip(fields[bi]))' not found in the case"))
    vm = num(String(fields[vi]))
    va = num(String(fields[ai]))
    (vm === nothing || va === nothing) && throw(ArgumentError("$(file):$(ln): non-numeric vm_pu/va_deg"))
    setVmVa!(node = net.nodeVec[idx], vm_pu = vm, va_deg = va)
    seen[idx] = true
  end
  for i = 1:n
    seen[i] || getNodeType(net.nodeVec[i]) == Isolated || throw(ArgumentError("$(file): no voltage row for bus $(i); the file does not match this case"))
  end
  return nothing
end

## True when at least one transformer of `net` has its tap released as an
## estimation state.
## The one statement the tap fallback makes, shared by all three surfaces
## that must carry it: the Web UI result summary, the tap table (which the
## fallback REPLACES, because estimated-looking positions that are model
## values are worse than no table), and se_diagnostics.md. One constant so
## the three cannot drift apart; a run whose tap positions are model values
## is otherwise indistinguishable from a successful tap estimation.
const _SE_TAP_FALLBACK_NOTE = "tap estimation did NOT converge: every tap was frozen at its MODEL position, so the tap positions in this result are model values, not estimates. The J of this run therefore measures the model tap positions, and a large J is expected when those positions are wrong."

_any_tap_released(net::Net)::Bool = any(b -> b.tap_est_mode !== :none, net.branchVec)

## Freeze every released tap back to its model position and report how many
## were frozen. Used by the service fallback: a set too thin to pin its taps
## should still deliver an estimate of the voltages.
function _freeze_all_tap_estimation!(net::Net)::Int
  frozen = 0
  for (i, b) in enumerate(net.branchVec)
    b.tap_est_mode === :none && continue
    setTapEstimation!(net; trafo = i, enabled = false)
    frozen += 1
  end
  return frozen
end

## Iterations a run actually needed, not just its LAST solve. With released
## taps the run solves twice: the estimation with the tap states, then the
## fixation run with the taps nailed to their mechanical step. The second is
## short (it starts from a converged state), so reporting only that one hid
## the expensive half: a CGMES run needed 36 to 40 iterations first and
## reported "3", which is how an iteration cap of 30 could look sufficient
## while it broke that run (maintainer, 2026-09-06).
function _se_reported_iterations(res)::Int
  base = res.iterations
  tf = res.tapFixation
  tf === nothing && return base
  return max(base, Int(get(tf, :iterations_estimation, 0)), Int(get(tf, :iterations_fixation, 0)))
end

## shared import for the SE services (same paths as the other services);
## returns the full ImportedCase so the run continues on its effective
## config, and logs the CGMES start decision like the power-flow service
function _se_service_import(case_path, config, run_id, config_file, output_dir, logfile, result_file, base_metadata; requested_format::Symbol = :auto)
  format = _se_case_format(case_path; requested = requested_format)
  # named separately: "unknown" alone would leave the caller guessing, and
  # the way out (naming the format) is not obvious from the file name
  format === :unknown && return nothing, format, _api_failure("se_unsupported_format", "The case format could not be determined. A bare .DAT is ambiguous (FOR001 network case vs FOR002 reference file); pass case_format = :dtf_for001 for a DTF network case.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  format in (:matpower, :cgmes, :scf, :dtf_for001) || return nothing, format, _api_failure("se_unsupported_format", "State estimation needs a MATPOWER, CGMES, Sparlectra Case Format or DTF case; got format $(format).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  imported = try
    _se_import_case(case_path, config; requested_format = requested_format)
  catch err
    err isa PowerFlowAborted && rethrow()
    return nothing, format, _api_failure("import_error", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  decision = get(imported.provenance, "cgmes_start_decision", nothing)
  decision === nothing || open(logfile, "a") do io
    println(io, decision)
  end
  return imported, format, nothing
end

"""
    _run_state_estimation_service(case_path, config_file, output_dir, run_id, measurement_file; kwargs...) -> SparlectraApiResult

Service backend of the Web UI "Run state estimation" action (SE phase 5).
Builds the net through the shared import paths, reads `measurement_file`
(measurement CSV v1, atomic), evaluates global observability (structural
islands plus FD-aware rank, phase 4), runs the bad-data diagnostics
(`runse_diagnostics`, sequential elimination on the configured budget) and
the final `runse!(updateNet = true)`, and writes the SE artifacts:
`measurements.csv` (copy), `se_diagnostics.md`, `se_view.md`,
`shunt_estimates.csv` (when shunts were released), and `se_state.csv` (the
chain anchor a later SE-started power flow consumes via `readSEStateCSV!`).

Options: `max_iter`, `tol`, `flatstart`, `robust`, `max_eliminations`,
`update_shunts`, `report_correlation` mirror the estimator keywords.
`tap_estimation = true` releases the tap of every in-service transformer
that carries a ratio tap changer (`setTapEstimation!` mode `:ratio`) before
the solve; the estimator then fixes each tap to its nearest mechanical step
and reports J before versus after the fixation (`se_tap_estimates.csv`).

Failure reasons: `se_unsupported_format`, `import_error`,
`invalid_measurements` (missing/unreadable/rejected file),
`se_not_observable` (observability quality `:not_observable`),
`se_not_converged`, plus the shared config failure.
"""
## headline mapping for the state-estimation timing file
const _SE_PERF_HEADLINE = (:case_loading_network_solver => "importing_case", :solver => "state_estimation", :postprocessing => "postprocessing_result", :artifact_writing => "writing_artifacts")

# One block per run about the critical measurements (issue #394): which
# method decided, how many rows are critical (their residual is
# structurally zero, so no residual test can ever flag an error there),
# the rows by id, and under :omega the nearly critical rows (wii below the
# 0.3 localizability guideline) as the graded information the rank test
# never had. Long lists are cut; the counts stay complete.
"""
    _se_criticality_lines(net, obs) -> Vector{String}

The criticality report of a state-estimation run as plain lines: the count
per method, the critical rows BY NAME (or "none", so a reader never has to
infer criticality from a wii column), and the nearly critical rows with
their wii. Written to the run log and to `se_diagnostics.md`.
"""
function _se_criticality_lines(net::Net, obs)::Vector{String}
  method = hasproperty(obs, :criticality_method) ? String(obs.criticality_method) : "rank"
  lines = String[]
  if obs.criticality_skipped
    push!(lines, string("critical measurements: classification skipped (method ", method, ", per-row budget exceeded)"))
    return lines
  end
  ids(idx) = [1 <= i <= length(net.measurements) ? String(net.measurements[i].id) : string("#", i) for i in idx]
  num = collect(obs.numerical_critical_measurement_indices)
  str = collect(obs.structural_critical_measurement_indices)
  push!(lines, string("critical measurements (", method, "): ", length(num), " numerical, ", length(str), " structural", hasproperty(obs, :structural_criticality_skipped) && obs.structural_criticality_skipped ? " (structural test skipped above the per-row budget)" : ""))
  # the names are the statement; a count alone sends the reader to the wii
  # column, which cannot tell critical from merely weak
  push!(lines, isempty(num) ? "  critical rows: none (a gross error on any row is detectable)" : string("  critical rows: ", join(first(ids(num), 20), ", "), length(num) > 20 ? string(", ... ", length(num) - 20, " more") : "", " (a gross error there stays invisible)"))
  if hasproperty(obs, :criticality_wii) && !isempty(obs.criticality_wii) && hasproperty(obs, :active_measurement_indices)
    wii = obs.criticality_wii
    act = collect(obs.active_measurement_indices)
    near = [(act[k], wii[k]) for k in eachindex(wii) if 1 <= k <= length(act) && wii[k] > 1e-8 && wii[k] < 0.3]
    if !isempty(near)
      sort!(near; by = x -> x[2])
      push!(lines, string("  nearly critical rows (wii below 0.3, not critical): ", length(near), ": ", join((string(first(ids([i])), " (", round(w; digits = 3), ")") for (i, w) in first(near, 10)), ", "), length(near) > 10 ? ", ..." : ""))
    end
  end
  return lines
end

function _se_log_criticality(io::IO, net::Net, obs)
  for line in _se_criticality_lines(net, obs)
    println(io, line)
  end
  return nothing
end

function _run_state_estimation_service(
  case_path::AbstractString,
  config_file::AbstractString,
  output_dir::AbstractString,
  run_id::String,
  measurement_file::AbstractString;
  # Solver settings default to `nothing`, NOT to a literal: a number here
  # would be a second source of truth next to the configuration, and it was
  # exactly that (service 6.0 against configuration 4.0 for k_suppress,
  # service 30 against configuration 20 for max_iter) which made the same
  # run behave differently depending on the entry point. `nothing` means
  # "not stated by the caller", and the run resolves it below against the
  # effective configuration, so the documented precedence holds: request,
  # then case configuration, then general configuration.
  max_iter::Union{Nothing,Int} = nothing,
  tol::Union{Nothing,Float64} = nothing,
  flatstart::Union{Nothing,Bool} = nothing,
  robust::Union{Nothing,Bool} = nothing,
  max_eliminations::Union{Nothing,Int} = nothing,
  update_shunts::Bool = false,
  report_correlation::Union{Nothing,Bool} = nothing,
  tap_estimation::Bool = false,
  k_eliminate::Union{Nothing,Float64} = nothing,
  robust_mode::Union{Nothing,Symbol} = nothing,
  robust_k1::Union{Nothing,Float64} = nothing,
  robust_k2::Union{Nothing,Float64} = nothing,
  k_suppress::Union{Nothing,Float64} = nothing,
  suppression_sigma::Union{Nothing,Float64} = nothing,
  case_format::Symbol = :auto,
  phase_callback = phase -> nothing,
  # the request's configuration overrides (Web UI run form, API caller):
  # the same top precedence level the power-flow path gives them; without
  # this an SE run read the configuration file only and every per-run
  # setting of the form was lost on the way (the CSV format among them)
  config_overrides::AbstractDict = Dict{String,Any}(),
)::SparlectraApiResult
  mkpath(output_dir)
  logfile = joinpath(output_dir, "run.log")
  result_file = joinpath(output_dir, "result.json")
  base_metadata = Dict{String,Any}("run_mode" => "se")

  # Phase timing, same machinery and same file name the power-flow path uses.
  # Without it a slow estimation could only be guessed at: the topology
  # precheck once cost 45 s of a 55 s run on a 25000-bus network and nothing
  # in the run directory showed it. Phases are coarse on purpose (one line
  # per step, not per iteration), matching the Web UI instrumentation rule.
  phase_recorder = PowerFlowPhaseTimingRecorder()
  # The recorder alone only fills the result metadata. The Web UI status page
  # reads the JOB's phase, so every phase has to be announced as well - the
  # power-flow path does that through the same callback. Without it the page
  # froze at the last phase the service layer set (preparing_configuration).
  se_phase = function (name::AbstractString)
    _start_service_phase!(phase_recorder, name)
    try
      phase_callback(String(name))
    catch err
      # a reporting problem must never take the run down
      err isa PowerFlowAborted && rethrow()
      @debug "state estimation: phase callback failed" phase = name exception = err
    end
    return nothing
  end
  se_total_start = time_ns()

  # the same precedence the power-flow path uses (resolve_config, D5):
  # case configuration file, the case file's deprecated block, general
  # file, defaults
  se_phase("preparing_configuration")
  config = try
    resolve_config(config_file, case_path, config_overrides).config
  catch err
    return _api_failure(_config_resolve_reason(err), sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  # what the caller did not state comes from the EFFECTIVE configuration
  # (case file first, then the general one), never from a literal in this
  # signature
  se_cfg = config.state_estimation
  max_iter = something(max_iter, se_cfg.max_iter)
  tol = something(tol, se_cfg.tol)
  flatstart = something(flatstart, se_cfg.flatstart)
  robust = something(robust, se_cfg.robust)
  k_eliminate = something(k_eliminate, se_cfg.k_eliminate)
  robust_mode = something(robust_mode, se_cfg.robust_mode)
  robust_k1 = something(robust_k1, se_cfg.robust_k1)
  robust_k2 = something(robust_k2, se_cfg.robust_k2)
  k_suppress = something(k_suppress, se_cfg.k_suppress)
  suppression_sigma = something(suppression_sigma, se_cfg.suppression_sigma)
  max_eliminations = something(max_eliminations, se_cfg.max_eliminations)
  report_correlation = something(report_correlation, se_cfg.report_residual_correlation)

  # bad-data threshold surface (GUI-exposed): validated HERE, after the
  # values are resolved, with a user-readable failure instead of a deep
  # solver error. Before, this ran on the keywords alone, so a bad value in
  # a configuration file reached the solver unchecked.
  robust_mode in (:off, :staged, :replacement) || return _api_failure("invalid_request", "se_robust_mode must be off, staged, or replacement (got $(robust_mode)).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  k_eliminate > 0.0 || return _api_failure("invalid_request", "se_k_eliminate must be positive.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  (robust_k1 > 0.0 && robust_k2 >= robust_k1) || return _api_failure("invalid_request", "se_robust_k1 must be positive and se_robust_k2 >= se_robust_k1 (got k1=$(robust_k1), k2=$(robust_k2)).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  max_eliminations >= 0 || return _api_failure("invalid_request", "se_max_eliminations must be >= 0 (got $(max_eliminations)).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  (k_suppress > 0.0 && suppression_sigma > 0.0) || return _api_failure("invalid_request", "se_k_suppress and se_suppression_sigma must be positive.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)


  se_phase("importing_case")
  imported, format, failure = _se_service_import(case_path, config, run_id, config_file, output_dir, logfile, result_file, base_metadata; requested_format = case_format)
  failure !== nothing && return failure
  net = imported.net
  config = imported.config

  # A Sparlectra Case Format case carries its measurements WITH the model
  # (#342), so a run on one needs no separate CSV: the set that came with the
  # case is used, and the run still writes measurements.csv as its artifact so
  # every SE run stays reproducible from its own output.
  from_case_file = format === :scf && !isfile(measurement_file) && !isempty(net.measurements)
  summary = if from_case_file
    open(logfile, "a") do io
      println(io, "measurements: ", length(net.measurements), " row(s) taken from the case file itself (no separate measurement set given)")
    end
    writeMeasurementsCSV(net; file = joinpath(output_dir, "measurements.csv"), busReference = :name)
    (total = length(net.measurements),)
  else
    if !isfile(measurement_file)
      hint = format === :scf ? " The case file carries no measurements either: generate a set, or export the case with its measurements." : ""
      return _api_failure("invalid_measurements", "Measurement file not found: $(measurement_file)." * hint, run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
    # case-binding gate: a set recorded for a different case would fail rows
    # deep in the reader ("bus X not found"); refuse it up front with the
    # actual cause. Sets without a binding (older files) pass with a log note.
    bound_case = _webui_measurement_set_case(String(measurement_file))
    # A case exported to SCF keeps the model of its source case, and it records
    # WHICH source it came from. A set bound to that source therefore fits this
    # run; without this the export silently invalidates every set the user
    # already generated.
    source_of_case = format === :scf ? _scf_source_reference(case_path) : ""
    if !isempty(bound_case) && bound_case != basename(String(case_path)) && bound_case != source_of_case
      return _api_failure("invalid_measurements", "Measurement set $(basename(String(measurement_file))) is bound to case $(bound_case) (recorded in the file), but this run uses $(basename(String(case_path))). Pick the set generated for this case, or regenerate one.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
    isempty(bound_case) && open(logfile, "a") do io
      println(io, "note: the measurement set carries no case binding (# case: comment); older set, association checked by name only")
    end
    read_summary = try
      readMeasurementsCSV!(net; file = measurement_file, replace = true)
    catch err
      return _api_failure("invalid_measurements", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
    cp(measurement_file, joinpath(output_dir, "measurements.csv"); force = true)
    read_summary
  end

  # Carry a duplicate finding into the run: it is the one property of a set that
  # explains a large J without any single measurement looking wrong.
  # Checked on the loaded set, not on the reader result: measurements carried
  # inside a case file take the other branch above and can be doubled just the
  # same.
  duplicate_rows = sum(values(_duplicate_measured_quantities(net.measurements)); init = 0)
  if duplicate_rows > 0
    base_metadata["se_duplicate_rows"] = duplicate_rows
    open(logfile, "a") do io
      println(io, "warning: $(duplicate_rows) of $(length(net.measurements)) measurements repeat a quantity that is already measured; expect an inflated J")
    end
  end

  # observability first (phase 4 two-stage verdict); a not-observable set is
  # a clean rejection, not a solver crash
  obs = try
    evaluate_global_observability(net)
  catch err
    return _api_failure("invalid_measurements", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  base_metadata["se_observability_quality"] = String(obs.quality)
  base_metadata["se_structural_islands"] = :structural_islands in obs.notes
  # critical measurements (issue #394): the classification method, the
  # counts and the rows themselves reach the metadata and the run log
  crit_num = collect(obs.numerical_critical_measurement_indices)
  crit_str = collect(obs.structural_critical_measurement_indices)
  crit_method = hasproperty(obs, :criticality_method) ? String(obs.criticality_method) : "rank"
  base_metadata["se_criticality_method"] = crit_method
  base_metadata["se_critical_measurements"] = length(crit_num)
  base_metadata["se_structural_critical_measurements"] = length(crit_str)
  base_metadata["se_criticality_skipped"] = obs.criticality_skipped
  if hasproperty(obs, :islands)
    base_metadata["se_islands_total"] = length(obs.islands)
    base_metadata["se_islands_measured"] = obs.n_measured_islands
    base_metadata["se_islands_unmeasured"] = obs.n_unmeasured_islands
  end
  if obs.quality == :not_observable
    msg = if hasproperty(obs, :islands)
      bad = [string("island ", i.island, " (", i.n_bus, " buses): ", String(i.result.quality)) for i in obs.islands if i.measured && i.result.quality == :not_observable]
      "The measurement set does not observe every measured island ($(join(bad, "; "))). Each island is estimated with its own reference; add measurements to the failing island or deactivate its rows."
    else
      "The measurement set does not observe the system (quality :not_observable)."
    end
    return _api_failure("se_not_observable", msg, run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  # released transformer taps: every in-service transformer with a
  # ratio tap changer becomes an estimation state (mode :ratio; PST release
  # stays an API-level choice via setTapEstimation!). Released AFTER the
  # observability verdict: tap states are deliberately outside the
  # observability count, like the gated current rows.
  # the generator documents an injected tap deviation inside the set; a run
  # WITHOUT tap estimation then reports a large J that is a model
  # discrepancy, not bad data. Say so before anyone chases measurements.
  # both reads below mine the generator's COMMENTS in the measurement CSV
  # (documented tap deviation, truth values); a set that came with the case
  # file has no such file, so they are skipped rather than guessed
  # Measurements carried inside a case file must yield the same two facts a
  # CSV set yields, or an SCF run silently loses them: which taps the generator
  # deviated (the estimator has to absorb exactly those) and the per-row truth
  # values (the delta file).
  case_provenance = from_case_file ? _scf_measurement_provenance(case_path) : Dict{String,Any}()
  set_tap_devs = if from_case_file
    NamedTuple[(branch = Int(get(d, "branch", 0)), steps = Float64(get(d, "steps", 0.0))) for d in get(case_provenance, "tap_deviations", Any[]) if d isa AbstractDict]
  else
    _measurement_set_tap_deviations(String(measurement_file))
  end
  # A set that DOCUMENTS its tap deviations answers the question the
  # estimator would otherwise have to be told, so those transformers are
  # released automatically - exactly the ones named, nothing else. As a mere
  # hint this produced a J the user could not explain, and the hint appeared
  # only after a converged run. Bad-data elimination is the wrong tool here:
  # the measurements are correct and the model is not, so dropping rows
  # would hide a model error instead of fixing it.
  auto_tap_branches = Int[]
  if !tap_estimation && !isempty(set_tap_devs)
    for d in set_tap_devs
      b = d.branch
      (1 <= b <= length(net.branchVec)) || continue
      try
        setTapEstimation!(net; trafo = b, enabled = true)
        push!(auto_tap_branches, b)
      catch err
        open(logfile, "a") do io
          println(io, "note: tap estimation could not be released automatically on branch ", b, ": ", sprint(showerror, err))
        end
      end
    end
    open(logfile, "a") do io
      if isempty(auto_tap_branches)
        println(io, "note: the measurement set documents a tap deviation (", join(("$(d.steps) step(s) on branch $(d.branch)" for d in set_tap_devs), ", "), ") but no transformer could be released; the estimator cannot absorb it (expect band :high clustered at the transformer)")
      else
        println(io, "tap estimation released automatically on ", length(auto_tap_branches), " transformer(s) named by the measurement set: ", join(auto_tap_branches, ", "))
      end
    end
    base_metadata["se_set_tap_deviation"] = true
    base_metadata["se_auto_released_taps"] = auto_tap_branches
  end

  # topology precheck: the check RESULT is logged on every run,
  # the clean case included (a silent skip would read as a pass); findings
  # are advisory and never block the estimation
  # the configuration of THIS run (#381): the resolved configuration of the
  # case with the caller's overrides folded in, installed for every
  # estimator call below; the registry is untouched afterwards
  run_cfg = _sparlectra_config_with(config; state_estimation = _se_config_with(config.state_estimation; max_iter = max_iter, tol = tol, flatstart = flatstart, robust = robust, k_eliminate = k_eliminate, robust_mode = robust_mode, robust_k1 = robust_k1, robust_k2 = robust_k2, k_suppress = k_suppress, suppression_sigma = suppression_sigma, max_eliminations = max_eliminations, report_residual_correlation = report_correlation, update_shunts = update_shunts, update_net = true))
  se_phase("topology_precheck")
  topo_pre = with_sparlectra_config(() -> validate_topology(net), run_cfg)
  # The findings are advisory and a large synthetic network produces many of
  # them (2170 on a 25000-bus case, 2155 of them "closed branch carries no
  # flow"). Printing every one buries the run log and nobody reads it, so the
  # log keeps a count per kind plus the first few of each; the complete list
  # stays in the metadata and in se_topology.csv.
  open(logfile, "a") do io
    println(io, topo_pre.summary)
    by_kind = Dict{Symbol,Vector{Any}}()
    for f in topo_pre.findings
      push!(get!(Vector{Any}, by_kind, f.kind), f)
    end
    for kind in sort!(collect(keys(by_kind)); by = string)
      fs = by_kind[kind]
      println(io, "  ", kind, ": ", length(fs), " finding(s)")
      for f in first(fs, 5)
        println(io, "    ", f.severity, " at ", f.location, ": ", f.evidence)
      end
      length(fs) > 5 && println(io, "    ... ", length(fs) - 5, " more (full list in se_topology.csv)")
    end
  end
  base_metadata["se_topology_findings"] = [Dict{String,Any}("stage" => String(f.stage), "kind" => String(f.kind), "location" => f.location, "evidence" => f.evidence, "severity" => String(f.severity)) for f in topo_pre.findings]
  base_metadata["se_topology_precheck_summary"] = topo_pre.summary

  tap_released = 0
  tap_skipped_machine = 0
  tap_machine_idxs = Int[]
  if tap_estimation
    for (k, br) in enumerate(net.branchVec)
      hasR = br.has_ratio_tap
      # DECLARED phase changers only (nameplate psi carried in
      # tap_est_alpha_deg by the mpc.sparlectra.tap_changers import): the
      # constructor gives every MATPOWER trafo default phase fields, and
      # mass-releasing those would add meaningless states everywhere
      hasP = _declared_phase_tap(br)
      (br.ratio != 0.0 && (hasR || hasP) && _branch_terminal_state(br) == :closed) || continue
      # machine (generator step-up) transformers are never mass-released:
      # the generator terminal voltage is set by the AVR, not observed, so
      # the tap state would absorb it. An explicit setTapEstimation! on a
      # specific transformer remains the deliberate way around this guard.
      if _is_machine_transformer(net, k)
        tap_skipped_machine += 1
        push!(tap_machine_idxs, k)
        continue
      end
      if hasR && hasP
        setTapEstimation!(net; trafo = k, mode = :both, alpha_deg = br.tap_est_alpha_deg)
      elseif hasP
        setTapEstimation!(net; trafo = k, mode = :pst, alpha_deg = br.tap_est_alpha_deg)
      else
        setTapEstimation!(net; trafo = k, mode = :ratio)
      end
      tap_released += 1
    end
    open(logfile, "a") do io
      tap_released == 0 && println(io, "tap estimation requested, but the case has no in-service transformer with a ratio or declared phase tap changer", tap_skipped_machine > 0 ? " (besides $(tap_skipped_machine) machine transformer(s), which are never mass-released)" : "")
      tap_skipped_machine > 0 && tap_released > 0 && println(io, "tap estimation: $(tap_skipped_machine) machine transformer(s) skipped (generator step-up; release explicitly via setTapEstimation! if intended)")
      # loud warning when the set's documented deviation sits on a SKIPPED
      # machine transformer: estimate taps can then never absorb it and J
      # stays large for a reason the user cannot see otherwise
      for d in set_tap_devs
        (1 <= d.branch <= length(net.branchVec) && _is_machine_transformer(net, d.branch)) || continue
        println(io, "WARNING: the set documents a tap deviation of $(d.steps) step(s) on branch $(d.branch), which is a machine transformer and SKIPPED by the mass release; the estimation cannot absorb it (release it explicitly via setTapEstimation!, or regenerate the set)")
      end
    end
  end

  # diagnostics (report) plus the final state-writing run (chain anchor)
  se_phase("estimation_diagnostics")
  diag = with_sparlectra_config(() -> runse_diagnostics(net), run_cfg)
  # eliminated rows leave the FINAL run: the diagnostics identified them on
  # its own copy, so deactivate them here or the headline J would keep
  # carrying rows the workflow already removed (seen: an injected 10-sigma
  # error, eliminated with stop reason consistent, still pushed the
  # headline to J = 149 at dof 42 instead of J ~ dof)
  for e in diag.eliminations
    (1 <= e.measurement_index <= length(net.measurements)) || continue
    net.measurements[e.measurement_index] = _set_measurement_active(net.measurements[e.measurement_index], false)
  end
  se_phase("state_estimation")
  # kept for the tap fallback below: runse! writes its state into the net
  # whether or not it converged, so the retry needs the state as it was
  # BEFORE the first attempt, not the diverged iterate it gave up on
  v_before_taps = _any_tap_released(net) ? [(nd._vm_pu, nd._va_deg) for nd in net.nodeVec] : nothing
  res = with_sparlectra_config(() -> runse!(net), run_cfg)

  # Released taps are extra states, and a measurement set that carries the
  # voltages fine can still be too thin to pin them: the estimate then does
  # not settle at all and the user gets nothing, although the SAME set
  # estimates cleanly without the taps (maintainer, 2026-09-06, case300 and
  # a CGMES delivery). So a non-convergence WITH released taps is not the
  # final answer: the taps are frozen back to their model position and the
  # estimation is repeated once. The log says it happened, because a silent
  # retry would hide that the reported taps are model values, not estimates.
  tap_fallback_used = false
  if !res.converged && _any_tap_released(net)
    frozen = _freeze_all_tap_estimation!(net)
    # The failed run wrote its state into the net (updateNet = true is not
    # gated on convergence), so without this the retry starts from the
    # DIVERGED iterate. Measured on the svedala set 2026-09-06: the failed
    # tap run left bus angles spread over -178 to +167 degrees where the
    # power flow had -26 to +40. With state_estimation.flatstart = true the
    # retry ignores the start state and the damage stays invisible, which is
    # why this never showed up; with a warm start it decides whether the
    # fallback converges at all.
    if v_before_taps !== nothing
      for (i, nd) in enumerate(net.nodeVec)
        nd._vm_pu, nd._va_deg = v_before_taps[i]
      end
    end
    open(logfile, "a") do io
      println(io, "state estimation did not converge with ", frozen, " released transformer tap(s); repeating WITHOUT tap estimation (the taps keep their model position). A set that cannot pin its taps needs more measurements around those transformers, not more iterations.")
    end
    @info "state estimation: tap estimation switched off after a non-converged run" released = frozen
    res = with_sparlectra_config(() -> runse!(net), run_cfg)
    tap_fallback_used = res.converged
    base_metadata["se_tap_estimation_fallback"] = tap_fallback_used
  end
  if !res.converged
    _write_service_performance_log!(output_dir, phase_recorder, se_total_start; headline = _SE_PERF_HEADLINE, status = "failed", label = "state-estimation")
    return _api_failure("se_not_converged", "State estimation did not converge within $(max_iter) iterations.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  se_phase("postprocessing_result")

  # artifacts
  open(joinpath(output_dir, "se_diagnostics.md"), "w") do io
    # same statement the result page carries: the tap positions below are
    # model values, so the J of this run measures THEM
    tap_fallback_used && println(io, "\n> **Tap estimation fallback.** ", _SE_TAP_FALLBACK_NOTE, "\n")
    print_se_diagnostics(io, diag; topN = 15, format = :markdown)
    # the critical rows by name (maintainer rule 2026-09-21: a detected
    # critical measurement is stated, never left to be read off the wii)
    println(io, "\n## Critical measurements\n")
    for line in _se_criticality_lines(net, obs)
      println(io, "- ", strip(line))
    end
    # topology findings of both stages append to the diagnostics report
    if !isempty(topo_pre.findings) || diag.topology_findings !== nothing
      println(io, "\n## Topology findings (advisory)\n")
      for f in topo_pre.findings
        println(io, "- precheck ", f.severity, ": `", f.kind, "` at ", f.location, " (", f.evidence, ")")
      end
      for tf in something(diag.topology_findings, NamedTuple[])
        println(io, "- classification strong: topology error suspected at station ", tf.location, " (suspects: ", tf.evidence, ")")
      end
      println(io, "\nFindings are advisory: nothing was switched or excluded. Use the hypothesis test on the result page for ranked recommendations.")
    end
  end
  view = se_view(net)
  open(joinpath(output_dir, "se_view.md"), "w") do io
    print_se_view(io, view; format = :markdown)
  end
  # One format for every diagnostic CSV this SE run writes (issue #376),
  # same config key the PowerFlow API uses for its own artifacts.
  se_csv_format = _resolve_detailed_csv_format(String(config.output.csv_format))
  se_csv_delim = se_csv_format.delimiter
  if res.shuntEstimates !== nothing
    open(joinpath(output_dir, "shunt_estimates.csv"), "w") do io
      println(io, join(("bus", "bus_name", "B_model_pu", "B_est_pu", "delta_pu", "frozen"), se_csv_delim))
      for r in res.shuntEstimates
        println(io, join((r.busIdx, _csv_field(String(r.busName), se_csv_delim, se_csv_format), _format_csv_number(Float64(r.B_model), se_csv_format), _format_csv_number(Float64(r.B_est), se_csv_format), _format_csv_number(Float64(r.delta), se_csv_format), r.frozen), se_csv_delim))
      end
    end
  end
  # machine transformers are CALCULATED, never estimated (maintainer
  # directive): their tap is no state variable, so after the estimation the
  # position is back-calculated from the AVR setpoint, the dispatch P, and
  # the MEASURED machine Q (the Qinj telemetry at the machine bus when the
  # set carries it). The evaluation marks these rows as "calculated".
  tap_calc_rows = Dict{String,Any}[]
  if tap_estimation
    trm = _transformer_mrids(net)
    for k in tap_machine_idxs
      mside = _machine_side(net, k)
      mside == 0 && continue
      qrow = findfirst(mm -> mm.active && mm.typ == QinjMeas && mm.busIdx == mside, net.measurements)
      prow = findfirst(mm -> mm.active && mm.typ == PinjMeas && mm.busIdx == mside, net.measurements)
      try
        bt = calcMachineTrafoTapFromSE(net; trafo = k, p_mw = prow === nothing ? nothing : net.measurements[prow].value, q_mvar = qrow === nothing ? nothing : net.measurements[qrow].value)
        push!(tap_calc_rows, Dict{String,Any}("branch" => k, "name" => getCompName(net.branchVec[k].comp), "mrid" => get(trm, k, ""), "mode" => "ratio", "electrical_step" => round(bt.electrical_step; digits = 3), "fixed_step" => bt.fixed_step, "electrical_shift_step" => 0.0, "fixed_shift_step" => 0, "out_of_range" => false, "fixed" => false, "frozen_reason" => "none", "source" => "calculated"))
        open(logfile, "a") do io
          println(io, "machine transformer ", getCompName(net.branchVec[k].comp), ": tap CALCULATED (not estimated) from AVR setpoint and machine telemetry: electrical step ", round(bt.electrical_step; digits = 2), " -> step ", bt.fixed_step, " (Q residual ", round(bt.q_residual_mvar; digits = 2), " MVar)")
        end
      catch err
        open(logfile, "a") do io
          println(io, "machine transformer ", getCompName(net.branchVec[k].comp), ": tap back-calculation not possible (", sprint(showerror, err), ")")
        end
      end
    end
  end

  # bad data at a glance (maintainer request 2026-08-27): every suspicious
  # or eliminated measurement with its network location in one CSV; the
  # diagnostics markdown keeps the full ranking, this is the extract to
  # open first. se_state.csv stays untouched: it is the machine-read chain
  # anchor (bus voltages/balances), not a measurement report.
  suspicious = diag.diagnostics.suspicious_measurements
  # rows the robust solve acted on (final state-writing run): replacement
  # stage 3 = SUPPRESSED with the fixed sigma, staged stages 1/2 =
  # down-weighted. The bad-data extract must say so per row, not just
  # eliminated true/false.
  suppressedIds = Set{String}(String(r.id) for r in something(res.robustRows, NamedTuple[]) if r.stage == 3)
  downweightedIds = Set{String}(String(r.id) for r in something(res.robustRows, NamedTuple[]) if r.stage in (1, 2))
  if !isempty(suspicious) || !isempty(diag.eliminations) || !isempty(suppressedIds)
    nby = _bus_name_by_idx(net)
    open(joinpath(output_dir, "se_bad_data.csv"), "w") do io
      println(io, join(("measurement_index", "id", "type", "bus", "from_bus", "to_bus", "normalized_residual", "wii", "localizable", "eliminated", "suppressed", "downweighted"), se_csv_delim))
      elim = Set(t.id for t in diag.eliminations)
      listed = Set{String}()
      cf(s) = _csv_field(string(s), se_csv_delim, se_csv_format)
      for row in suspicious
        bus = ""
        fb = ""
        tb = ""
        if 1 <= row.measurement_index <= length(net.measurements)
          m = net.measurements[row.measurement_index]
          m.busIdx !== nothing && (bus = get(nby, m.busIdx, string(m.busIdx)))
          if m.branchIdx !== nothing && 1 <= m.branchIdx <= length(net.branchVec)
            br = net.branchVec[m.branchIdx]
            fb = get(nby, Int(br.fromBus), string(Int(br.fromBus)))
            tb = get(nby, Int(br.toBus), string(Int(br.toBus)))
          end
        end
        push!(listed, String(row.id))
        println(
          io,
          join(
            (row.measurement_index, cf(row.id), cf(row.typ), cf(bus), cf(fb), cf(tb), _format_csv_number(Float64(row.normalized_residual), se_csv_format), isnan(row.wii) ? "" : _format_csv_number(Float64(row.wii), se_csv_format), row.localizable, row.id in elim, String(row.id) in suppressedIds, String(row.id) in downweightedIds),
            se_csv_delim,
          ),
        )
      end
      for t in diag.eliminations
        t.id in listed && continue
        push!(listed, String(t.id))
        println(
          io,
          join(
            (t.measurement_index, cf(t.id), cf(t.typ), "", "", "", _format_csv_number(Float64(t.normalized_residual_before), se_csv_format), "", true, true, String(t.id) in suppressedIds, String(t.id) in downweightedIds),
            se_csv_delim,
          ),
        )
      end
      # suppressed rows below the elimination threshold (possible when
      # k_suppress < k_eliminate) still belong in the extract
      for r in something(res.robustRows, NamedTuple[])
        r.stage == 3 || continue
        rid = String(r.id)
        rid in listed && continue
        push!(listed, rid)
        mi = r.measurement_index
        bus = ""
        fb = ""
        tb = ""
        ty = ""
        if 1 <= mi <= length(net.measurements)
          m = net.measurements[mi]
          ty = string(m.typ)
          m.busIdx !== nothing && (bus = get(nby, m.busIdx, string(m.busIdx)))
          if m.branchIdx !== nothing && 1 <= m.branchIdx <= length(net.branchVec)
            br = net.branchVec[m.branchIdx]
            fb = get(nby, Int(br.fromBus), string(Int(br.fromBus)))
            tb = get(nby, Int(br.toBus), string(Int(br.toBus)))
          end
        end
        println(io, join((mi, cf(rid), cf(ty), cf(bus), cf(fb), cf(tb), _format_csv_number(Float64(r.t), se_csv_format), "", "", false, true, false), se_csv_delim))
      end
    end
  end
  # measured/truth/estimated delta file (generator v2): written when the
  # set carries truth_value comments (generated sets). measured and sigma
  # come from the loaded rows, truth is the generator's noise-free value,
  # estimated is measured - residual at the final diagnostics state;
  # eliminated rows keep measured/truth but no estimate contribution flag.
  # Tap rows compare the set's documented generation deviation with the
  # estimation's fixed step (authoritative per-tap table: se_tap_estimates.csv).
  deltas_written = false
  truthById = Dict{String,Float64}()
  for (k, v) in get(case_provenance, "truth_values", Dict{String,Any}())
    v isa Real && (truthById[String(k)] = Float64(v))
  end
  for line in (from_case_file ? String[] : eachline(measurement_file))
    startswith(line, "# truth_value,") || continue
    parts = split(chopprefix(line, "# truth_value,"), ",")
    length(parts) >= 2 || continue
    tv = tryparse(Float64, String(parts[end]))
    tv === nothing && continue
    truthById[String(join(parts[1:(end-1)], ","))] = tv
  end
  if !isempty(truthById)
    nbyD = _bus_name_by_idx(net)
    elimD = Set(t.id for t in diag.eliminations)
    residById = Dict{String,Float64}()
    for row in diag.final_diagnostics.measurement_ranking
      residById[String(row.id)] = row.residual
    end
    open(joinpath(output_dir, "se_deltas.csv"), "w") do io
      println(io, "# sparlectra-se-deltas v1")
      println(io, join(("kind", "id", "type", "bus", "from_bus", "to_bus", "measured", "truth", "estimated", "delta_meas_truth", "delta_est_truth", "delta_est_meas", "sigma", "t_est", "eliminated"), se_csv_delim))
      fmt(v) = _format_csv_number(round(v; sigdigits = 8), se_csv_format)
      cfD(s) = _csv_field(string(s), se_csv_delim, se_csv_format)
      for m in net.measurements
        haskey(truthById, m.id) || continue
        tv = truthById[m.id]
        busD = m.busIdx !== nothing ? get(nbyD, m.busIdx, string(m.busIdx)) : ""
        fbD = ""
        tbD = ""
        if m.branchIdx !== nothing && 1 <= m.branchIdx <= length(net.branchVec)
          brD = net.branchVec[m.branchIdx]
          fbD = get(nbyD, Int(brD.fromBus), string(Int(brD.fromBus)))
          tbD = get(nbyD, Int(brD.toBus), string(Int(brD.toBus)))
        end
        rres = get(residById, m.id, nothing)
        estv = rres === nothing ? nothing : m.value - rres
        println(
          io,
          join(
            ("meas", cfD(m.id), cfD(m.typ), cfD(busD), cfD(fbD), cfD(tbD), fmt(m.value), fmt(tv), estv === nothing ? "" : fmt(estv), fmt(m.value - tv), estv === nothing ? "" : fmt(estv - tv), estv === nothing ? "" : fmt(estv - m.value), fmt(m.sigma), rres === nothing ? "" : fmt(abs(rres) / m.sigma), m.id in elimD),
            se_csv_delim,
          ),
        )
      end
      devByBranch = Dict{Int,Float64}(d.branch => d.steps for d in set_tap_devs)
      for t in something(res.tapEstimates, NamedTuple[])
        dev = get(devByBranch, t.branch, 0.0)
        # a tap is never measured: `measured` and the deltas against it stay
        # empty rather than repeating the truth and pretending a zero residual
        println(io, join(("tap", cfD(t.name), cfD(t.mode), "", "", "", "", fmt(dev), t.fixed_step_1, "", fmt(Float64(t.fixed_step_1) - dev), "", "", "", false), se_csv_delim))
      end
      deltas_written = true
    end
  end

  if res.tapEstimates !== nothing || !isempty(tap_calc_rows)
    open(joinpath(output_dir, "se_tap_estimates.csv"), "w") do io
      println(io, join(("branch", "name", "mrid", "mode", "alpha_deg", "electrical_step", "fixed_step", "electrical_shift_step", "fixed_shift_step", "r1_est", "r2_est", "out_of_range", "fixed", "frozen_reason", "source"), se_csv_delim))
      cfT(s) = _csv_field(string(s), se_csv_delim, se_csv_format)
      for t in something(res.tapEstimates, NamedTuple[])
        println(
          io,
          join(
            (t.branch, cfT(t.name), cfT(t.mrid), cfT(t.mode), _format_csv_number(Float64(t.alpha_deg), se_csv_format), _format_csv_number(round(t.electrical_step_1; digits = 4), se_csv_format), t.fixed_step_1, _format_csv_number(round(t.electrical_step_2; digits = 4), se_csv_format), t.fixed_step_2, t.r1_est, t.r2_est, t.out_of_range, t.fixed, t.frozen_reason == :none ? "" : t.frozen_reason, "estimated"),
            se_csv_delim,
          ),
        )
      end
      for c in tap_calc_rows
        println(io, join((c["branch"], cfT(c["name"]), cfT(c["mrid"]), "ratio", _format_csv_number(0.0, se_csv_format), c["electrical_step"], c["fixed_step"], _format_csv_number(0.0, se_csv_format), 0, "", "", false, false, "", "calculated"), se_csv_delim))
      end
    end
  end
  se_phase("writing_artifacts")
  writeSEStateCSV(net; file = joinpath(output_dir, "se_state.csv"), format = String(config.output.csv_format))

  # band verdict from the POST-elimination report: it must describe the
  # same state as the headline J (the pre-elimination report still tells
  # its story in se_diagnostics.md)
  verdict = diag.final_diagnostics.objective
  open(logfile, "a") do io
    println(io, "State estimation on ", basename(case_path))
    println(io, "measurements: ", summary.total, " rows from ", from_case_file ? "the case file" : basename(String(measurement_file)))
    println(io, "observability: ", obs.quality, :structural_islands in obs.notes ? " (structural islands)" : "")
    _se_log_criticality(io, net, obs)
    println(io, "converged in ", res.iterations, " iteration(s); J = ", round(res.objectiveJ; digits = 6), ", dof = ", res.dof, ", band reason = ", verdict.reason, res.activeObjective === nothing ? "" : string("; J_active = ", round(res.activeObjective.j; digits = 6), " (dof ", res.activeObjective.dof, ", without ", res.activeObjective.suppressed, " suppressed row(s); band verdict stays on J)"))
    if diag.topology_findings !== nothing
      # the stage-2 classification REPLACES the eliminations-exhausted
      # interpretation: this is a suspected topology error, not bad data
      for tf in diag.topology_findings
        println(io, "topology error suspected at station ", tf.location, " (eliminations exhausted, band high; suspects: ", tf.evidence, isempty(tf.notes) ? "" : string("; notes: ", join(tf.notes, ", ")), ")")
      end
      println(io, "suspicious measurements: ", length(diag.diagnostics.suspicious_measurements), "; eliminations: ", length(diag.eliminations), " (stop: ", diag.stop_reason, ", classified as topology, use the hypothesis test)")
    else
      println(io, "suspicious measurements: ", length(diag.diagnostics.suspicious_measurements), "; eliminations: ", length(diag.eliminations), " (stop: ", diag.stop_reason, ")")
    end
    effRobustMode = robust_mode === :off && robust ? :staged : robust_mode
    if effRobustMode === :staged
      println(io, "robust R modification active (staged, k1=", robust_k1, ", k2=", robust_k2, "); ", res.robustRows === nothing ? 0 : length(res.robustRows), " row(s) weighted down in the final iteration")
    elseif effRobustMode === :replacement
      println(io, "bad-data suppression active (replacement, k_suppress=", k_suppress, ", suppression sigma=", suppression_sigma, "); ", res.robustRows === nothing ? 0 : length(res.robustRows), " row(s) suppressed on the converged state; statistics keep the original sigmas")
    end
    if res.tapEstimates !== nothing
      if res.tapFixation !== nothing
        tf = res.tapFixation
        println(io, "tap estimation: ", length(res.tapEstimates), " transformer(s); J before fixation = ", round(tf.j_before; sigdigits = 4), " (dof ", tf.dof_before, "), J after fixation = ", round(tf.j_after; sigdigits = 4), " (dof ", tf.dof_after, ")", tf.fixed ? "" : "; NOT fixed (estimation run did not converge)")
        tf.offgrid_residual && println(io, "note :offgrid_tap_residual: the band failure appears only through the fixation; this is an off-grid true position or a wrong step table, NOT bad data")
      else
        println(io, "tap estimation: every released tap was frozen by the guards; no tap state entered the solve")
      end
      for t in res.tapEstimates
        println(io, "  ", t.name, isempty(t.mrid) ? "" : " (mRID " * t.mrid * ")", ": electrical step ", round(t.electrical_step_1; digits = 2), " -> fixed step ", t.fixed_step_1, t.mode == :ratio ? "" : string(", shift step ", round(t.electrical_step_2; digits = 2), " -> ", t.fixed_step_2), t.out_of_range ? " (OUT OF RANGE, clamped)" : "", t.frozen_reason == :none ? "" : " (FROZEN: $(t.frozen_reason))")
      end
    end
    (!isempty(suspicious) || !isempty(diag.eliminations) || !isempty(suppressedIds)) && println(io, "bad data: se_bad_data.csv (", length(suspicious), " suspicious, ", length(diag.eliminations), " eliminated, ", length(suppressedIds), " suppressed, ", length(downweightedIds), " down-weighted, with network locations)")
    deltas_written && println(io, "deltas: se_deltas.csv (measured vs generator truth vs estimated, per row and per released tap)")
    println(io, "Artifacts: measurements.csv, se_diagnostics.md, se_view.md, se_state.csv", res.shuntEstimates !== nothing ? ", shunt_estimates.csv" : "", res.tapEstimates !== nothing ? ", se_tap_estimates.csv" : "", (!isempty(suspicious) || !isempty(diag.eliminations) || !isempty(suppressedIds)) ? ", se_bad_data.csv" : "", deltas_written ? ", se_deltas.csv" : "")
  end

  metadata = merge(
    base_metadata,
    Dict{String,Any}(
      "input_format_detected" => String(format),
      "se_measurement_rows" => summary.total,
      "se_converged" => res.converged,
      # the maximum over the solves of the run, see _se_reported_iterations
      "se_iterations" => _se_reported_iterations(res),
      "se_iterations_last_solve" => res.iterations,
      "se_objective" => res.objectiveJ,
      "se_dof" => res.dof,
      "se_objective_active" => res.activeObjective === nothing ? nothing : res.activeObjective.j,
      "se_dof_active" => res.activeObjective === nothing ? nothing : res.activeObjective.dof,
      "se_suppressed_rows" => res.activeObjective === nothing ? 0 : res.activeObjective.suppressed,
      "se_band_reason" => String(verdict.reason),
      "se_j_within_3sigma" => res.jWithin3Sigma,
      "se_suspicious" => length(diag.diagnostics.suspicious_measurements),
      "se_eliminations" => length(diag.eliminations),
      "se_robust" => robust,
      "se_deltas" => deltas_written,
      "se_robust_mode" => String(robust_mode === :off && robust ? :staged : robust_mode),
      # the numbers that DECIDED this run, so a result can be read without
      # guessing which of configuration, service or form won (the reason
      # this is here at all: a run estimated with k_suppress 6.0 while the
      # configuration said 4.0, and nothing in the result showed it)
      "se_max_iter" => max_iter,
      "se_tol" => tol,
      "se_k_eliminate" => k_eliminate,
      "se_robust_k1" => robust_k1,
      "se_robust_k2" => robust_k2,
      "se_k_suppress" => k_suppress,
      "se_suppression_sigma" => suppression_sigma,
      "se_update_shunts" => update_shunts,
      "se_islands_estimated" => res.islands === nothing ? 1 : count(i -> i.estimated, res.islands),
      "se_islands_skipped" => res.islands === nothing ? 0 : count(i -> !i.estimated, res.islands),
      "se_tap_estimation" => tap_estimation,
      "se_tap_count" => res.tapEstimates === nothing ? 0 : length(res.tapEstimates),
      "se_tap_skipped_machine" => tap_skipped_machine,
      "se_tap_frozen" => res.tapEstimates === nothing ? 0 : count(t -> t.frozen_reason != :none, res.tapEstimates),
      "se_tap_fixed" => res.tapFixation === nothing ? false : res.tapFixation.fixed,
      "se_tap_j_before" => res.tapFixation === nothing ? nothing : res.tapFixation.j_before,
      "se_tap_dof_before" => res.tapFixation === nothing ? nothing : res.tapFixation.dof_before,
      "se_tap_j_after" => res.tapFixation === nothing ? nothing : res.tapFixation.j_after,
      "se_tap_dof_after" => res.tapFixation === nothing ? nothing : res.tapFixation.dof_after,
      "se_tap_offgrid_residual" => res.tapFixation === nothing ? false : res.tapFixation.offgrid_residual,
      "se_topology_station_findings" => diag.topology_findings === nothing ? nothing : [Dict{String,Any}("location" => tf.location, "evidence" => tf.evidence, "notes" => [String(n) for n in tf.notes]) for tf in diag.topology_findings],
      # per-trafo rows for the result page table (steps, not raw r values)
      "se_tap_estimates" => res.tapEstimates === nothing && isempty(tap_calc_rows) ? nothing : vcat([Dict{String,Any}("branch" => t.branch, "name" => t.name, "mrid" => t.mrid, "mode" => String(t.mode), "electrical_step" => round(t.electrical_step_1; digits = 3), "fixed_step" => t.fixed_step_1, "electrical_shift_step" => round(t.electrical_step_2; digits = 3), "fixed_shift_step" => t.fixed_step_2, "out_of_range" => t.out_of_range, "fixed" => t.fixed, "frozen_reason" => String(t.frozen_reason), "source" => "estimated") for t in something(res.tapEstimates, NamedTuple[])], tap_calc_rows),
      "artifact_status" => "completed",
      "solver_status" => "completed",
      "service_status" => "completed",
      "run_status" => "completed",
    ),
  )

  # J/dof FIRST: J alone says nothing (it grows with the number of
  # measurements), and a reader who sees "J = 104" on a 14-bus case reads an
  # alarm where J/dof = 1.1 says the set is healthy.
  message = string(
    "State estimation completed - ", _se_reported_iterations(res), " iteration(s), J/dof = ",
    res.dof > 0 ? string(round(res.objectiveJ / res.dof; digits = 2)) : "n/a",
    " (J = ", round(res.objectiveJ; digits = 3), ", dof ", res.dof, ", band ", verdict.reason, "), ",
    length(diag.diagnostics.suspicious_measurements), " suspicious, ",
    length(diag.eliminations), " eliminated.",
    diag.topology_findings === nothing ? "" : string(" Topology error suspected at ", join((tf.location for tf in diag.topology_findings), ", "), "; run the hypothesis test."),
    res.tapFixation === nothing ? "" : string(" Tap fixation (", length(res.tapEstimates), " transformer(s)): J ", round(res.tapFixation.j_before; sigdigits = 3), " -> ", round(res.tapFixation.j_after; sigdigits = 3), res.tapFixation.offgrid_residual ? "; off-grid tap residual, not bad data" : "", "."),
    get(base_metadata, "se_set_tap_deviation", false) == true ? (isempty(get(base_metadata, "se_auto_released_taps", Int[])) ? " Hint: the set documents a tap deviation that could not be released automatically; enable 'estimate taps'." : string(" Tap estimation was released automatically on ", length(get(base_metadata, "se_auto_released_taps", Int[])), " transformer(s) named by the set.")) : "",
    duplicate_rows > 0 ? string(" Warning: ", duplicate_rows, " measurement(s) repeat an already measured quantity - that alone inflates J; regenerate the set.") : "",
  )
  # timing file plus the raw phase sequence in result.json, same as the
  # power-flow path, so the Web UI run page can show where the time went
  metadata["service_phase_timings"] = _write_service_performance_log!(output_dir, phase_recorder, se_total_start; headline = _SE_PERF_HEADLINE, label = "state-estimation")
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

"""
    _run_pf_from_se_service(case_path, config_file, output_dir, run_id, se_state_file, se_run_id, se_mode) -> SparlectraApiResult

SE-started power-flow chain run: builds the net,
restores the estimated state from `se_state_file` (the `se_state.csv` of a
preceding SE run), and runs `runpf_from_se!` with `se_mode` (`"se_state"` or
`"se_snapshot"`). The persistent model is never mutated (the snapshot's
balance takeover lives on the working import only). Metadata records the
start mode, the referenced SE run id, and the slack pickup.
"""
function _run_pf_from_se_service(case_path::AbstractString, config_file::AbstractString, output_dir::AbstractString, run_id::String, se_state_file::AbstractString, se_run_id::AbstractString, se_mode::AbstractString)::SparlectraApiResult
  mkpath(output_dir)
  logfile = joinpath(output_dir, "run.log")
  result_file = joinpath(output_dir, "result.json")
  base_metadata = Dict{String,Any}("run_mode" => "powerflow_se_start", "se_run_id" => String(se_run_id), "se_start_mode" => String(se_mode))

  se_mode in ("se_state", "se_snapshot") || return _api_failure("invalid_request", "se_start_mode must be \"se_state\" or \"se_snapshot\", got \"$(se_mode)\".", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)

  # the same precedence the power-flow path uses (resolve_config, D5):
  # case configuration file, the case file's deprecated block, general
  # file, defaults
  config = try
    resolve_config(config_file, case_path).config
  catch err
    return _api_failure(_config_resolve_reason(err), sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  imported, format, failure = _se_service_import(case_path, config, run_id, config_file, output_dir, logfile, result_file, base_metadata)
  failure !== nothing && return failure
  net = imported.net
  config = imported.config

  isfile(se_state_file) || return _api_failure("se_state_missing", "SE state artifact not found for run $(se_run_id); run the state estimation first.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  r = try
    readSEStateCSV!(net; file = se_state_file)
    runpf_from_se!(net, 40, 1e-8, 0; mode = Symbol(se_mode), method = :rectangular)
  catch err
    err isa PowerFlowAborted && rethrow()
    return _api_failure("se_start_failed", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  r.converged || return _api_failure("powerflow_not_converged", "SE-started power flow did not converge ($(r.iterations) iterations).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  calcNetLosses!(net)

  # deviation of the PF solution from the estimated (chosen) SE state: the
  # per-bus answer to "how far does the model pull away from the measured
  # operating point". Full table as artifact, extremes in metadata/log.
  st = _se_start_state(net)
  devrows = NamedTuple[]
  if st !== nothing
    iso = Set(net.isoNodes)
    name_by_idx = _bus_name_by_idx(net)
    for (i, nd) in enumerate(net.nodeVec)
      i in iso && continue
      vm_pf = something(nd._vm_pu, NaN)
      va_pf = something(nd._va_deg, NaN)
      dva = va_pf - st.va[i]
      # angle wrap to (-180, 180]
      dva = mod(dva + 180.0, 360.0) - 180.0
      push!(devrows, (bus = i, name = get(name_by_idx, i, string(i)), mrid = _bus_mrid(net, i; name_by_idx = name_by_idx), vm_se = st.vm[i], vm_pf = vm_pf, dvm = vm_pf - st.vm[i], va_se = st.va[i], va_pf = va_pf, dva = dva))
    end
    devi_format = _resolve_detailed_csv_format(String(config.output.csv_format))
    devi_delim = devi_format.delimiter
    open(joinpath(output_dir, "se_pf_deviation.csv"), "w") do io
      println(io, join(("bus", "name", "mrid", "vm_se_pu", "vm_pf_pu", "dvm_pu", "va_se_deg", "va_pf_deg", "dva_deg"), devi_delim))
      for rw in devrows
        println(
          io,
          join(
            (rw.bus, _csv_field(String(rw.name), devi_delim, devi_format), _csv_field(String(rw.mrid), devi_delim, devi_format), _format_csv_number(Float64(rw.vm_se), devi_format), _format_csv_number(Float64(rw.vm_pf), devi_format), _format_csv_number(Float64(rw.dvm), devi_format), _format_csv_number(Float64(rw.va_se), devi_format), _format_csv_number(Float64(rw.va_pf), devi_format), _format_csv_number(Float64(rw.dva), devi_format)),
            devi_delim,
          ),
        )
      end
    end
  end
  maxdvm = isempty(devrows) ? NaN : maximum(abs(rw.dvm) for rw in devrows)
  maxdva = isempty(devrows) ? NaN : maximum(abs(rw.dva) for rw in devrows)
  busdvm = isempty(devrows) ? "" : devrows[argmax([abs(rw.dvm) for rw in devrows])].name
  busdva = isempty(devrows) ? "" : devrows[argmax([abs(rw.dva) for rw in devrows])].name
  meandvm = isempty(devrows) ? NaN : sum(abs(rw.dvm) for rw in devrows) / length(devrows)

  open(logfile, "a") do io
    println(io, "SE-started power flow (", se_mode, ") on ", basename(case_path))
    println(io, "start state: se_state.csv of SE run ", se_run_id)
    println(io, "converged in ", r.iterations, " iteration(s)")
    println(io, "slack pickup: ", round(r.slack_pickup_mw; digits = 6), " MW / ", round(r.slack_pickup_mvar; digits = 6), " MVar")
    if !isempty(devrows)
      println(io, "deviation PF vs SE state: max |dVm| = ", round(maxdvm; sigdigits = 4), " pu at ", busdvm, ", mean |dVm| = ", round(meandvm; sigdigits = 4), " pu, max |dVa| = ", round(maxdva; sigdigits = 4), " deg at ", busdva, " (full table: se_pf_deviation.csv)")
    end
  end

  metadata = merge(
    base_metadata,
    Dict{String,Any}(
      "iterations" => r.iterations,
      "slack_pickup_mw" => r.slack_pickup_mw,
      "slack_pickup_mvar" => r.slack_pickup_mvar,
      "se_pf_max_dvm_pu" => maxdvm,
      "se_pf_max_dvm_bus" => busdvm,
      "se_pf_mean_dvm_pu" => meandvm,
      "se_pf_max_dva_deg" => maxdva,
      "se_pf_max_dva_bus" => busdva,
      "artifact_status" => "completed",
      "solver_status" => "completed",
      "service_status" => "completed",
      "run_status" => "completed",
    ),
  )

  result = _api_result(
    run_id = run_id,
    status = :succeeded,
    success = true,
    solution_available = true,
    iterations = r.iterations,
    reason = nothing,
    message = string("SE-started power flow (", se_mode, ") converged in ", r.iterations, " iteration(s); slack pickup ", round(r.slack_pickup_mw; digits = 4), " MW."),
    casefile = String(case_path),
    config_file = String(config_file),
    output_dir = String(output_dir),
    logfile = logfile,
    result_file = result_file,
    metadata = metadata,
  )
  return _finalize_api_result(result)
end

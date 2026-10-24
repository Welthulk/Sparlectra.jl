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
#          through the shared import paths, read a measurement table set,
#          check observability, run runse_diagnostics plus the final
#          runse!(updateNet = true) (with released taps in rounds: result
#          run, diagnostics on its end steps), and persist the SE artifacts including
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
solves exactly like the power-flow service).
"""
function _se_import_case(case_path::AbstractString, config; requested_format::Symbol = :auto)::ImportedCase
  format = _detect_case_format(String(case_path); requested = requested_format)
  # format policy stays here (the wording is this service's contract)
  format in (:scf, :matpower, :cgmes, :dtf_for001, :powsybl) || throw(ArgumentError("State estimation needs a MATPOWER, CGMES, Sparlectra Case Format, DTF or PowSyBl case; got format $(format)."))
  return import_case(String(case_path), config; requested_format = requested_format, run_kind = :state_estimation)
end

_se_import_case_net(case_path::AbstractString, config; requested_format::Symbol = :auto)::Net = _se_import_case(case_path, config; requested_format).net

## Format for decisions that only ask "which kind of case is this" (bus
## reference in the measurement CSV, CGMES specifics). A bare `.DAT` is
## ambiguous without an explicit request (FOR001 network vs FOR002
## reference), and that ambiguity must not abort a run that does not even
## depend on the answer, so it resolves to :unknown here.
function _se_case_format(case_path::AbstractString; requested::Symbol = :auto)::Symbol
  return first(_se_case_format_with_reason(case_path; requested = requested))
end

## The detector refuses two things with an ArgumentError: an ambiguous .DAT
## and an explicit format that does not fit the file content. Both messages
## name the way out, so they travel with the `:unknown` verdict instead of
## being replaced by a generic text (until 0.17.2 the state estimation on a
## MATPOWER case with `case_format = scf` reported the .DAT sentence).
function _se_case_format_with_reason(case_path::AbstractString; requested::Symbol = :auto)::Tuple{Symbol,String}
  return try
    (_detect_case_format(String(case_path); requested = requested), "")
  catch err
    err isa ArgumentError || rethrow()
    (:unknown, sprint(showerror, err))
  end
end

## from-run truth resolution (measurement generator v2): resolve `run_id` in
## the persistent run index under `run_root`, enforce the case binding, and
## load the solved voltages into `net` WITHOUT re-solving. SE runs read
## se_state.csv (bus by name); PF runs read the detail CSV
## bus_voltages_complex.csv, which exists only when the run wrote detail
## CSVs. A missing artifact rejects up front with the hint to repeat the run
## with the detail CSV enabled; there is deliberately no silent re-solve.
"""
    _se_truth_run_warning(metadata) -> String

The warning a measurement set inherits from the run it takes its truth
state from: empty for a clean run, otherwise one sentence naming why the
state is no solution of the model (state estimation outside the band,
the fallback without tap estimation, a non-converged run). Flows derived from
such a state violate the node balances, and the set then contradicts the
model before any estimation (a set built from a fallback
run far outside the band produced KCL findings at every passive bus).
"""
function _se_truth_run_warning(metadata::AbstractDict)::String
  reasons = String[]
  band = String(get(metadata, "se_band_reason", ""))
  (isempty(band) || band == "ok") || push!(reasons, "its J was outside the band ($(band))")
  get(metadata, "se_tap_fallback_used", false) == true && push!(reasons, "its taps were not estimated (fallback without tap estimation)")
  status = String(get(metadata, "numerical_status", ""))
  (isempty(status) || status == "converged") || push!(reasons, "it did not converge ($(status))")
  isempty(reasons) && return ""
  return string("the source run is no solution of the model: ", join(reasons, ", "), "; flows derived from its state violate the node balances, prefer truth state 'fresh solve'")
end

function _se_run_metadata(outdir::AbstractString)::Dict{String,Any}
  f = joinpath(outdir, "result.json")
  isfile(f) || return Dict{String,Any}()
  parsed = try
    scf_json_parse(read(f, String))
  catch err
    # an unreadable result file is no reason to refuse the state; the
    # warning is simply not available and the log says why
    @warn "measurement generator: result.json of the truth run could not be read, no source-run warning" file = f exception = err
    return Dict{String,Any}()
  end
  md = parsed isa AbstractDict ? get(parsed, "metadata", Dict{String,Any}()) : Dict{String,Any}()
  return md isa AbstractDict ? Dict{String,Any}(String(k) => v for (k, v) in md) : Dict{String,Any}()
end

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
    return (kind = "se", timestamp = ts, source_file = "se_state.csv", warning = _se_truth_run_warning(_se_run_metadata(outdir)))
  end
  f = joinpath(outdir, "bus_voltages_complex.csv")
  isfile(f) || throw(ArgumentError("run $(run_id) has no bus_voltages_complex.csv artifact; repeat the run with the detailed result CSV enabled"))
  _se_apply_pf_voltage_csv!(net, f)
  return (kind = "pf", timestamp = ts, source_file = "bus_voltages_complex.csv", warning = _se_truth_run_warning(_se_run_metadata(outdir)))
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


## The one statement the tap fallback makes, shared by the surfaces that
## carry it: the Web UI result summary and tap section, and
## se_diagnostics.md. One constant so they cannot drift apart.
const _SE_TAP_FALLBACK_NOTE = "tap estimation did NOT converge: the result is the run without tap estimation, every tap keeps its transferred step (reason fallback_not_converged). The J of this run measures those steps, and a large J is expected when they are wrong."
const _SE_TAP_REJECTED_NOTE = "tap estimation was REJECTED: its estimate, fixed to the mechanical steps, fitted the measurements worse than the transferred steps (worse_fit), so the result is the run without tap estimation and every tap keeps its transferred step (reason fallback_worse_fit). The set does not pin these taps reliably."

## the note of a fallback by its reason (`se_tap_fallback_reason`)
_se_tap_fallback_note(metadata) = get(metadata, "se_tap_fallback_reason", "") == "worse_fit" ? _SE_TAP_REJECTED_NOTE : _SE_TAP_FALLBACK_NOTE

## run.log line 1 of a fallback: the library line for a failed run 1, the
## acceptance rule of the staged start for a converged run 1 it rejected
## (`pre` the staged pre-solve with the taps frozen)
function _se_tap_fallback_line1_service(reason::Symbol, ite::Int, res, pre)::String
  reason === :worse_fit || return _se_tap_fallback_line1(reason, ite)
  return string("tap estimation: run 1 rejected (worse_fit, J = ", round(res.objectiveJ; sigdigits = 6), " at dof ", res.dof, " after the fixation against J = ", round(pre.objectiveJ; sigdigits = 6), " at dof ", pre.dof, " without tap estimation), falling back")
end

## Freeze every released tap back to its model position and report how many
## were frozen (the staged start solves once with the taps frozen).
function _freeze_all_tap_estimation!(net::Net)::Int
  frozen = 0
  for (i, b) in enumerate(net.branchVec)
    b.tap_est_mode === :none && continue
    setTapEstimation!(net; trafo = i, enabled = false)
    frozen += 1
  end
  return frozen
end

## Which transformers are released and in which mode, so a staged start can
## freeze them for one solve and release exactly the same set again (the
## direction of a phase changer comes from its changer data, 0.60.0).
_tap_release_snapshot(net::Net) = [(k, b.tap_est_mode) for (k, b) in enumerate(net.branchVec) if b.tap_est_mode !== :none]

function _restore_tap_release!(net::Net, snap)
  for (k, mode) in snap
    setTapEstimation!(net; trafo = k, mode = mode)
  end
  return length(snap)
end

## Iterations a run actually needed, not just its LAST solve. With released
## taps the run solves twice: the estimation with the tap states, then the
## fixation run with the taps nailed to their mechanical step. The second is
## short (it starts from a converged state), so reporting only that one hid
## the expensive half: a CGMES run needed 36 to 40 iterations first and
## reported "3", which is how an iteration cap of 30 could look sufficient
## while it broke that run.
function _se_reported_iterations(res)::Int
  base = res.iterations
  tf = res.tapFixation
  # after the fallback the reported run is run 2 alone; the iterations of
  # run 1 stand in its run.log line
  (tf === nothing || tf.fallback_used) && return base
  return max(base, Int(get(tf, :iterations_estimation, 0)), Int(get(tf, :iterations_fixation, 0)))
end

## Columns of se_tap_estimates.csv, one row per tap regulator: the
## identification (branch, regulator index, transformer name, changer, both
## mRIDs), then the steps on the changer's own step numbers: the
## transferred step of the input (`step_transferred`), the step the
## estimate takes after rounding (`step_estimated`, the transferred step
## where the row was not estimated or was out of range) and the continuous
## estimate before rounding (`step_estimated_raw`, neutral step plus s).
## Then the changer data (kind, side, alpha, du), Delta U = s du, branch
## ratio and angle at the estimate, the deviation estimated minus
## transferred in steps, the standard deviation of the step from the
## covariance (`sigma_s`) with the distance of the estimate from its
## nearest step in units of it (`round_distance_sigma`, a hint), the step
## from the control target of a regulator kept out automatically
## (`step_from_target`, not applied), whether the step was fixed, the
## source of the transformer impedance (`table`, `analytic`, `nominal`),
## the band flag `out_of_range` (estimate more than half a step outside the
## band; with status `estimated` it was fixed to the band end step), the
## status (`estimated`, `out_of_range`, `not_estimated`, `calculated`),
## reason, note (for example `p_from_schedule` on a calculated machine transformer)
## and the notice text. A state is never written into a step column. The
## metadata rows of the result page carry the same keys plus `role`.
const _SE_TAP_CSV_COLUMNS = ["branch", "regulator", "name", "changer", "mrid", "changer_mrid", "step_transferred", "step_estimated", "step_estimated_raw", "kind", "side", "alpha_deg", "du", "du_est", "ratio_est", "angle_est_deg", "deviation_steps", "sigma_s", "round_distance_sigma", "step_from_target", "fixed", "reactance_source", "out_of_range", "status", "reason", "note", "notice"]

## a step value as an Int when it is a whole step (start positions carry
## round-off from the branch fields), else as it is
_se_tap_step(x::Real) = abs(x - round(x)) < 1e-6 ? round(Int, x) : Float64(x)

## a step for the run log: whole steps as integers, others with two decimals
function _se_tap_step_text(x::Real)::String
  v = _se_tap_step(x)
  return v isa Int ? string(v) : string(round(v; digits = 2))
end

## the signed deviation of a regulator row for the run log
_se_tap_deviation_text(d::Real)::String = d > 0 ? string("+", _se_tap_step_text(d)) : _se_tap_step_text(d)

## the run-log mark of a regulator row whose estimate lies outside the
## band: status out_of_range keeps the transferred step, a fixed row took
## the band end step, any other row (withdrawn later, not fixed) only the mark
function _se_tap_out_of_range_text(t)::String
  t.out_of_range || return ""
  t.status === :out_of_range && return " (OUT OF RANGE, transferred step kept)"
  return t.fixed ? " (OUT OF RANGE, fixed to the band end step)" : " (OUT OF RANGE)"
end

## one tap regulator row of SEResult.tapEstimates as the record of the
## outputs: the CSV columns plus role. The library fields step_model,
## step_fixed and s_est (relative to neutral) become step_transferred,
## step_estimated and step_estimated_raw (absolute, neutral_step + s_est).
## Symbols become strings, the values that do not exist for the changer
## (alpha of a ratio changer, du of a linear changer or a table) stay
## `nothing`.
function _se_tap_row_record(t)::Dict{String,Any}
  return Dict{String,Any}(
    "branch" => t.branch, "regulator" => t.regulator, "name" => t.name, "changer" => t.changer,
    "mrid" => t.mrid, "changer_mrid" => t.changer_mrid, "kind" => String(t.kind), "side" => String(t.side),
    "step_transferred" => _se_tap_step(t.step_model), "step_estimated" => _se_tap_step(t.step_fixed), "step_estimated_raw" => Float64(t.neutral_step + t.s_est),
    "alpha_deg" => t.alpha_deg, "du" => t.du, "du_est" => t.du_est,
    "ratio_est" => t.ratio_est, "angle_est_deg" => t.angle_est_deg, "deviation_steps" => _se_tap_step(t.deviation_steps),
    "sigma_s" => t.sigma_s, "round_distance_sigma" => t.round_distance_sigma, "step_from_target" => t.step_from_target === nothing ? nothing : _se_tap_step(t.step_from_target),
    "out_of_range" => t.out_of_range, "fixed" => t.fixed, "reactance_source" => String(t.reactance_source), "status" => String(t.status), "reason" => String(t.reason), "note" => String(t.note),
    "role" => String(t.role), "notice" => t.notice,
  )
end

## one record as a CSV line in the column order; `nothing` and the symbol
## value "none" (kind, side, reason, note) are empty cells, numbers
## follow the regional CSV format, indices and whole steps are integers
function _se_tap_csv_line(rec::AbstractDict, delim, fmt)::String
  function cell(c)
    v = get(rec, c, nothing)
    v === nothing && return ""
    (v isa Bool || v isa Integer) && return string(v)
    v isa AbstractFloat && return _format_csv_number(v, fmt)
    (c in ("kind", "side", "reason", "note") && v == "none") && return ""
    return _csv_field(String(v), delim, fmt)
  end
  return join((cell(c) for c in _SE_TAP_CSV_COLUMNS), delim)
end

## the regulator rows of SEResult.tapEstimates the estimation saw: rows of
## transformers without a changer of the released kind (regulator 0) and of
## transformers the configuration keeps out (reason excluded,
## machine_transformer or q_not_measured_in_group, counted as
## se_tap_excluded and se_tap_machine) are left out
_se_tap_regulator_rows(rows) = rows === nothing ? NamedTuple[] : [t for t in rows if t.regulator > 0 && !(t.reason in (:excluded, :machine_transformer, :q_not_measured_in_group))]

## number of distinct transformers among the tap regulator rows
_se_tap_transformer_count(rows)::Int = length(unique(t.branch for t in _se_tap_regulator_rows(rows)))

## the comment lines of a measurement CSV without their `# ` prefix, split
## into those before the header (case binding, provenance, taps table) and
## those after the data rows (the generator's truth blocks), so a rewrite of
## the file keeps what the generator recorded
function _measurement_csv_comments(path::AbstractString)
  head = String[]
  tail = String[]
  seen_header = false
  for (ln, raw) in enumerate(eachline(path))
    ln == 1 && continue   # the version comment, written by the writer itself
    s = strip(raw)
    if startswith(s, "#")
      push!(seen_header ? tail : head, String(strip(lstrip(s, '#'))))
    elseif startswith(s, "type")
      seen_header = true
    end
  end
  return head, tail
end

## shared import for the SE services (same paths as the other services);
## returns the full ImportedCase so the run continues on its effective
## config, and logs the CGMES start decision like the power-flow service
function _se_service_import(case_path, config, run_id, config_file, output_dir, logfile, result_file, base_metadata; requested_format::Symbol = :auto)
  format, reason = _se_case_format_with_reason(case_path; requested = requested_format)
  # named separately: "unknown" alone would leave the caller guessing; the
  # detector's own message names the way out (the format to pass, or auto)
  format === :unknown && return nothing, format, _api_failure("se_unsupported_format", string("The case format could not be determined. ", reason), run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  format in (:matpower, :cgmes, :scf, :dtf_for001, :powsybl) || return nothing, format, _api_failure("se_unsupported_format", "State estimation needs a MATPOWER, CGMES, Sparlectra Case Format, DTF or PowSyBl case; got format $(format).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
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
  # the line susceptances from c1 carry the frequency into the estimator
  _record_network_frequency!(base_metadata, imported.net, logfile)
  return imported, format, nothing
end

# A measurement table (one set per timestamp) given as the measurement
# file. Gated by state_estimation.measurement_table. The format of the case
# picks the key index: a CGMES delivery resolves mRIDs through its model
# (cgmesMeasurementKeys), every other format resolves element names
# (nameMeasurementKeys). The run estimates the FIRST set of the file; the
# import report (every set, every row that did not become a measurement) is
# written to measurement_import.txt and summarized in the log. Returns
# `(summary, nothing)` or `(nothing, reason)`.
# The key index a measurement table of this case is read with, and the one
# the converter and the Web UI generator write with: mRIDs for a CGMES case,
# element names for every other format. `nothing` when a CGMES import kept
# no CGMES model.
function _se_measurement_index(imported, format::Symbol)
  format === :cgmes || return nameMeasurementKeys(imported.net)
  cg = get(imported.provenance, "cgmes_result", nothing)
  return cg === nothing ? nothing : cgmesMeasurementKeys(cg)
end

"""
    _se_convert_measurements_v1(case_path, v1_path, out; config_file = DEFAULT_SPARLECTRA_CONFIG_PATH, requested_format = :auto) -> NamedTuple

`convertMeasurementsV1` for a case file: imports the case the way the SE
service does (effective configuration of `config_file` and the case;
`requested_format` as the service's `case_format`, needed for a bare DTF
`.DAT`) and
writes the table with the key index the service reads it with (mRIDs for a
CGMES case, element names otherwise). Returns the result of
`convertMeasurementsV1`. Throws `ArgumentError` when a CGMES import keeps no
CGMES model, and the reader's error for a v1 file it refuses.
"""
function _se_convert_measurements_v1(case_path::AbstractString, v1_path::AbstractString, out::AbstractString; config_file::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH, requested_format::Symbol = :auto)
  config = resolve_config(config_file, case_path).config
  imported = _se_import_case(case_path, config; requested_format = requested_format)
  index = _se_measurement_index(imported, _se_case_format(case_path; requested = requested_format))
  index === nothing && throw(ArgumentError("$(basename(case_path)): the CGMES import kept no CGMES model to key the table by mRID"))
  return convertMeasurementsV1(v1_path, imported.net; out = out, index = index)
end

function _se_load_measurement_table!(net, imported, format::Symbol, measurement_file::String, config, logfile, output_dir::String)
  name = basename(measurement_file)
  # read first: a table the reader refuses (no version line, another
  # version) gets the reader's message whether the switch is on or off.
  # ArgumentError is the reader's refusal of the file as a whole.
  rows = try
    readMeasurementTable(measurement_file)
  catch err
    err isa ArgumentError || rethrow()
    return nothing, sprint(showerror, err)
  end
  config.state_estimation.measurement_table || return nothing, "$(name) is a measurement table. Reading measurement tables is switched off; set state_estimation.measurement_table: true."
  index = _se_measurement_index(imported, format)
  index === nothing && return nothing, "$(name) is a measurement table for a CGMES case, but the import kept no CGMES model to resolve its mRIDs."
  # rows without a sigma cell take the level sigma of the effective
  # configuration (case configuration first), not of the process default;
  # the same holds for the zero-injection sigma (the sigma floor of the
  # import and the sigma of the rows applyMeasurementSet! adds)
  # the sigma every row without a sigma cell gets, for the run artifact
  filled = Dict{Int,Tuple{String,Symbol}}()
  imp = _import_measurement_table(rows, index; sigma_by_voltage_level = config.state_estimation.sigma_by_voltage_level, zero_injection_sigma = config.state_estimation.zero_injection_sigma_mva, filled_sigmas = filled)
  if isempty(imp.sets)
    write(joinpath(output_dir, "measurement_import.txt"), measurementImportReport(imp))
    return nothing, "$(name) holds no measurement row that can be assigned to a set; see measurement_import.txt."
  end
  set = imp.sets[1]
  # applied before the report is written, so the report carries the
  # zero-injection rows of state_estimation.zero_injection for this set
  applyMeasurementSet!(net, set; zero_injection = config.state_estimation.zero_injection !== :off, zero_injection_sigma = config.state_estimation.zero_injection_sigma_mva)
  write(joinpath(output_dir, "measurement_import.txt"), measurementImportReport(imp))
  zi = set.counts[:zero_injection]
  open(logfile, "a") do io
    println(io, "measurements: measurement table ", name, ", ", length(imp.sets), " set(s); estimating set ", isempty(set.timestamp) ? "(no timestamp)" : set.timestamp, ": ", set.counts[:rows], " rows, ", length(set.measurements), " measurements (active ", set.counts[:active], ", inactive ", set.counts[:inactive], "), ", length(set.problems), " row(s) not used, see measurement_import.txt")
    _se_log_zero_injection(io, config, zi)
    length(imp.sets) > 1 && println(io, "note: the other sets of the table are not estimated in this run: ", join((s.timestamp for s in imp.sets[2:end]), ", "))
  end
  return (total = length(set.measurements), zero_injection = zi, sigma_origin = set.sigma_origin, table = (rows = rows, set = set, filled = filled), index = index), nothing
end

# The measurements.csv artifact of an SE run is the measurement set the run
# estimates, written as a measurement table (same writer, same keys as the
# shipped tables), and the run estimates exactly the set that table reads
# back to. A table set keeps its own rows: the rows that became a
# measurement, in file order, with their keys, value text, source and valid
# flag. A row without a sigma cell gets the sigma the import gave it plus
# the source token sigma_level or sigma_default, so the artifact reads back
# the same under any configuration. Zero-injection rows the run added
# (state_estimation.zero_injection) and a set carried by a case file are
# written through the converter's row builder (source se_run or case_file,
# zero-injection rows with the token zero_injection). Measurements without
# a place in the table (a row the table cannot key, such as a flow on an
# unnamed link) are not
# written: they are named as comment lines and in run.log and stay in the
# estimate behind the table set. Returns the artifact rows for the verdict
# file: per data row its file line, key, type and the index of its
# measurement in net.measurements, the key index (for the kV base of U),
# and the id map of a case-file set (old id to table id).
const _SE_ARTIFACT_ROW = NamedTuple{(:timestamp, :measurement_mrid, :terminal_mrid, :node_mrid, :type, :value, :sigma, :source, :valid),Tuple{String,String,String,String,String,String,String,String,Bool}}

function _se_artifact_row(r::MeasurementRow, filled)::_SE_ARTIFACT_ROW
  sigma = _measurement_table_number(r, :sigma)
  source = r.source
  if filled !== nothing
    sigma = filled[1]
    tok = filled[2] === :level ? "sigma_level" : "sigma_default"
    source = isempty(source) ? tok : string(source, ";", tok)
  end
  return _SE_ARTIFACT_ROW((r.timestamp, r.measurement_mrid, r.terminal_mrid, r.node_mrid, r.type, _measurement_table_number(r, :value), sigma, source, r.valid))
end

function _se_write_measurement_artifact!(net::Net, index::MeasurementKeyIndex, file::AbstractString, comments::Vector{String}, logfile, zero_injection_sigma::Real; table = nothing)
  rows = _SE_ARTIFACT_ROW[]
  # the measurement index (in net.measurements as loaded) of every row
  owner = Int[]
  unplaced = String[]
  id_map = Dict{String,String}()
  n_meas = length(net.measurements)
  if table !== nothing
    set = table.set
    id_pos = Dict{String,Int}(m.id => k for (k, m) in enumerate(set.measurements))
    line_owner = Dict{Int,Int}()
    for (id, ls) in set.sources
      for l in ls
        line_owner[l] = id_pos[id]
      end
    end
    for r in table.rows
      haskey(line_owner, r.line) || continue
      push!(rows, _se_artifact_row(r, get(table.filled, r.line, nothing)))
      push!(owner, line_owner[r.line])
    end
    first_added = length(set.measurements) + 1
  else
    first_added = 1
  end
  # rows the run added (zero injection) or the set of a case file
  added = first_added:n_meas
  conv = _measurement_table_from_measurements(net, net.measurements[added], index; source = table === nothing ? "case_file" : "se_run")
  without = Set{Int}(p.line for p in conv.problems)
  # the added rows belong to the set of the table: same timestamp
  ts = table === nothing ? "" : table.set.timestamp
  append!(rows, (_SE_ARTIFACT_ROW((ts, Base.tail(Tuple(r))...)) for r in conv.rows))
  append!(owner, [added[k] for k in eachindex(added) if !(k in without)])
  table === nothing && merge!(id_map, conv.ids)
  for p in conv.problems
    push!(unplaced, string(p.id, ": ", p.problem, ": ", p.detail))
  end
  head = copy(comments)
  write_table(rs, extra) = writeMeasurementTable(file, rs; comments = vcat(head, [string("not in the table: ", u) for u in extra]))
  write_table(rows, unplaced)
  # read back with the run's index; a row that does not read back (a case
  # file row at a bus that is not energized) leaves the table and joins the
  # measurements without a place, and the table is written once more
  # an empty table (no measurement of a case file has a place) has no set
  read_back() = isempty(rows) ? MeasurementSet("", Measurement[], Dict{String,Vector{Int}}(), Dict{Symbol,Int}(), NamedTuple{(:line, :problem, :detail),Tuple{Int,Symbol,String}}[], Dict{String,Symbol}()) : only(_import_measurement_table(readMeasurementTable(file), index; zero_injection_sigma = zero_injection_sigma).sets)
  back = read_back()
  data0 = 2 + length(head) + length(unplaced)
  if !isempty(back.problems)
    bad = Set{Int}(p.line - data0 for p in back.problems)
    for p in back.problems
      push!(unplaced, string("line ", p.line, ": ", p.problem, ": ", p.detail))
    end
    keep = [k for k in eachindex(rows) if !(k in bad)]
    rows = rows[keep]
    owner = owner[keep]
    write_table(rows, unplaced)
    back = read_back()
    data0 = 2 + length(head) + length(unplaced)
    isempty(back.problems) || error("measurements.csv: the run artifact does not read back: $(back.problems[1])")
  end
  carried = Set{Int}(owner)
  rest = Measurement[net.measurements[k] for k in 1:n_meas if !(k in carried)]
  isempty(rest) || open(logfile, "a") do io
    println(io, "note: ", length(rest), " measurement(s) have no place in the measurement table and are not in measurements.csv (named there as comments); they stay in the estimate")
  end
  # the run estimates the set read back from the artifact, then the rest
  empty!(net.measurements)
  append!(net.measurements, back.measurements)
  append!(net.measurements, rest)
  pos = Dict{String,Int}(m.id => k for (k, m) in enumerate(back.measurements))
  line_meas = Dict{Int,Int}()
  for (id, ls) in back.sources
    for l in ls
      line_meas[l] = pos[id]
    end
  end
  lines = collect(data0 .+ (1:length(rows)))
  keys_ = String[isempty(r.measurement_mrid) ? (isempty(r.terminal_mrid) ? r.node_mrid : r.terminal_mrid) : r.measurement_mrid for r in rows]
  return (index = index, lines = lines, keys = keys_, types = String[r.type for r in rows], measurement = Int[line_meas[l] for l in lines], id_map = id_map, rows = length(rows), not_in_table = length(rest))
end

# measurements_verdict.csv: the verdict per data row of measurements.csv,
# joined by the line number of the row in that file (unique; a key with
# its type is not, a bus sum spans the rows of all its members). Residuals
# are in the unit of the table row (MW, MVAr, A, degrees, U in kV), the
# residual of the measurement the row belongs to; a bus sum gives each
# member row the residual of the sum. An eliminated row carries the
# residual and normalized residual of the solve before its removal
# (`eliminations` holds the trace rows of every diagnostics pass of the
# run, with released taps one pass per bad-data round), every other active row those of
# the final diagnostics. A zero-injection
# constraint (state_estimation.zero_injection = constraint) has normalized
# residual 0 and constraint = true. `reason` lists the verdict tokens:
# inactive, eliminated, suppressed, downweighted, suspicious, critical, and
# from the range check of the voltage measurements (`range_rows`)
# range_check_substituted (the row entered the run as a pseudo measurement,
# its residual is that of the substitute) and deenergized (the row sits at
# a bus that is not energized, or at its branches, and left the run).
const _SE_VERDICT_HEADER = ["line", "key", "type", "measurement_id", "residual", "normalized_residual", "wii", "localizable", "eliminated", "constraint", "reason"]

function _se_write_measurement_verdict(file::AbstractString, net::Net, art, diag, eliminations, res, critical::Set{Int}, zi_constraint::Bool, fmt; range_rows = NamedTuple[])
  delim = string(fmt.delimiter)
  cf(s) = _csv_field(string(s), fmt.delimiter, fmt)
  num(x) = x === nothing || isnan(x) ? "" : _format_csv_number(Float64(x), fmt)
  ranking = Dict{Int,Any}(row.measurement_index => row for row in diag.final_diagnostics.measurement_ranking if row.measurement_index >= 1)
  elim = Dict{Int,Any}(t.measurement_index => t for t in eliminations)
  robust = something(res.robustRows, NamedTuple[])
  suppressed = Set{String}(String(r.id) for r in robust if r.stage == 3)
  downweighted = Set{String}(String(r.id) for r in robust if r.stage in (1, 2))
  range_sub = Set{Int}(r.measurement_index for r in range_rows if r.reason === :range_check_substituted)
  range_dead = Set{Int}(j for r in range_rows if r.reason === :deenergized for j in r.deactivated)
  flagged = 0
  open(file, "w") do io
    println(io, join(_SE_VERDICT_HEADER, delim))
    for k in eachindex(art.lines)
      mi = art.measurement[k]
      m = net.measurements[mi]
      # U rows are in kV in the table, VmMeas residuals in pu of the bus
      scale = m.typ == VmMeas ? art.index.bus_kV[m.busIdx] : 1.0
      residual = nothing
      rn = nothing
      wii = nothing
      localizable = false
      suspicious = false
      eliminated = haskey(elim, mi)
      if eliminated
        t = elim[mi]
        residual = t.residual_before * scale
        rn = t.normalized_residual_before
        wii = t.wii
        # elimination requires a localizable row
        localizable = true
      elseif m.active && haskey(ranking, mi)
        row = ranking[mi]
        residual = row.residual * scale
        rn = row.normalized_residual
        wii = row.wii
        localizable = row.localizable
        suspicious = row.suspicious
      end
      constraint = zi_constraint && _zero_injection_constraint(m)
      constraint && (rn = 0.0)
      tokens = String[]
      # the final solve deactivates the eliminated rows in net.measurements;
      # inactive means valid=false in the table
      (m.active || eliminated) || push!(tokens, "inactive")
      mi in range_dead && push!(tokens, "deenergized")
      mi in range_sub && push!(tokens, "range_check_substituted")
      eliminated && push!(tokens, "eliminated")
      m.id in suppressed && push!(tokens, "suppressed")
      m.id in downweighted && push!(tokens, "downweighted")
      suspicious && push!(tokens, "suspicious")
      mi in critical && push!(tokens, "critical")
      any(t -> t in ("eliminated", "suppressed", "downweighted", "suspicious", "range_check_substituted"), tokens) && (flagged += 1)
      println(io, join((art.lines[k], cf(art.keys[k]), cf(art.types[k]), cf(m.id), num(residual), num(rn), num(wii), localizable, eliminated, constraint, cf(join(tokens, ";"))), delim))
    end
  end
  return flagged
end

# state_estimation.sigma_by_voltage_level of one run: per nominal voltage the
# number of measurements and where their sigma came from (row, level,
# default), plus the configured absolute sigmas in pu. `origins` is the
# sigma origin of a measurement table set; every other set (measurement CSV
# v1, case file) carries a sigma in each row. The same rows go to run.log,
# se_diagnostics.md and the result metadata (se_sigma_by_voltage_level).
function _se_sigma_level_report(net::Net, config, origins::AbstractDict)
  summary = _sigma_level_summary(net, net.measurements; origins = origins, levels = config.state_estimation.sigma_by_voltage_level)
  rows = Dict{String,Any}[]
  for l in summary.levels
    push!(rows, Dict{String,Any}(
      "vn_kv" => l.vn_kv, "configured" => l.configured, "measurements" => l.measurements,
      "sigma_row" => l.row, "sigma_level" => l.level, "sigma_default" => l.default,
      "rel_pct" => l.rel_pct,
      "p_pu" => l.sigma_pu.p, "q_pu" => l.sigma_pu.q, "i_pu" => l.sigma_pu.i, "u_pu" => l.sigma_pu.u,
    ))
  end
  return (rows = rows, unused = summary.unused, unlocated = summary.unlocated, configured = !isempty(config.state_estimation.sigma_by_voltage_level))
end

_se_sigma_pu_text(x) = x === nothing ? "default" : string(round(Float64(x); sigdigits = 4))

# the per-level lines of run.log and se_diagnostics.md
function _se_sigma_level_lines(report)::Vector{String}
  lines = String[]
  push!(lines, string("measurement sigma by voltage level", report.configured ? "" : " (`state_estimation.sigma_by_voltage_level` not set: every row without its own sigma takes the default)", ":"))
  for r in report.rows
    push!(lines, string("  ", r["vn_kv"], " kV: ", r["measurements"], " measurement(s), sigma from the row ", r["sigma_row"], ", from the level ", r["sigma_level"], ", default ", r["sigma_default"],
      r["configured"] ? string("; level sigma in pu: P ", _se_sigma_pu_text(r["p_pu"]), ", Q ", _se_sigma_pu_text(r["q_pu"]), ", I ", _se_sigma_pu_text(r["i_pu"]), ", U ", _se_sigma_pu_text(r["u_pu"]), r["rel_pct"] > 0 ? string(", at least ", r["rel_pct"], " % of the value") : "") : "; no row configured"))
  end
  for v in report.unused
    push!(lines, string("  note: ", v, " kV is configured in `state_estimation.sigma_by_voltage_level`, but no bus of this case has that nominal voltage"))
  end
  report.unlocated > 0 && push!(lines, string("  ", report.unlocated, " measurement(s) without a bus or branch are not counted"))
  return lines
end

# state_estimation.zero_injection on a set read from a measurement CSV or a
# case file: the same rows applyMeasurementSet! adds for a table, added
# before measurements.csv is written so the artifact is the estimated set.
# A generated set already carries ZI rows or balance rows at every passive
# bus and gets nothing added. Returns the number of rows added.
function _se_apply_zero_injection!(net::Net, config, logfile)::Int
  added = config.state_estimation.zero_injection !== :off ? _add_missing_zero_injections!(net.measurements, net; sigma = config.state_estimation.zero_injection_sigma_mva) : 0
  open(logfile, "a") do io
    _se_log_zero_injection(io, config, added)
  end
  return added
end

# the log line of the zero-injection switch, the same for every SE path
function _se_log_zero_injection(io::IO, config, added::Int)
  if config.state_estimation.zero_injection !== :off
    mode = config.state_estimation.zero_injection
    # constraint: every ZI row of the set is solved as an equality
    # constraint, the added ones and those the set brought along
    how = mode === :constraint ? ", every ZI row solved as an equality constraint" : string(", sigma ", config.state_estimation.zero_injection_sigma_mva, " MW/Mvar per row")
    println(io, "zero injection: ", added, " row(s) added at ", added ÷ 2, " passive bus(es) without an injection measurement", how, " (state_estimation.zero_injection: ", mode, ")")
  else
    println(io, "zero injection: off, no rows added (state_estimation.zero_injection: off)")
  end
  return nothing
end

# the run.log summary line of the numerics of the reported run: the matrix
# of the last Gauss-Newton step with its condition number (1-norm
# estimate), J and dof, and under zero_injection = constraint where the
# multipliers are
function _se_conditioning_line(res)::String
  c = res.conditioning
  c === nothing && return string("numerics: J = ", round(res.objectiveJ; sigdigits = 6), ", dof = ", res.dof)
  what = c.matrix === :kkt ? "KKT matrix" : "gain matrix"
  kappa = isfinite(c.condition_number) ? string(round(c.condition_number; sigdigits = 3)) : "Inf (singular)"
  zi = res.zeroInjectionMultipliers === nothing ? "" : string(", multipliers of ", length(res.zeroInjectionMultipliers), " bus(es) in se_zero_injection_multipliers.csv")
  return string("numerics: zero_injection ", c.zero_injection, ", condition number of the ", what, " ", kappa, " (1-norm estimate), J = ", round(res.objectiveJ; sigdigits = 6), ", dof = ", res.dof, zi)
end

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

## The bad-data sequence of the SE service on `net`: returns the result
## run, the diagnostics of the whole sequence, every eliminated row (with
## its `round`), the staged-start and fallback state of the last round,
## the number of rounds and the configuration of the last run.
##
## Bad data with released taps runs on the FIXED END STEPS of the result
## run. One round: the result run (taps released,
## estimated and fixed; or the fallback without tap estimation), then the
## diagnostics with every tap frozen at its end step, and the row they
## eliminate leaves the set. When a round eliminated a row, the next round
## repeats the result run with released taps on the cleaned set; the rounds
## end with the first round that eliminates nothing. One row per round: a
## second row of the same pass would be judged on end steps the eliminated
## row had pulled (sp_case14 with -2 kV on the
## 20 kV bus of the deviated transformer: a full pass removed that row and
## two good 20 kV voltages, one row per round only that row). The budget
## state_estimation.max_eliminations counts over all rounds, a round with
## no budget left only diagnoses. Measured on a study case with a
## 0 kV voltage row, before the range check: the released run of round
## 1 does not converge in 50 iterations and falls back, its diagnostics on
## the transferred steps eliminate the 0 kV row, round 2 estimates three
## transformers and gives J 212.27 at dof 924. The order before
## (diagnostics with free taps, then the fallback for good) ended at J
## 414.30 with every tap at its transferred step.
## Without a released tap the result run and the diagnostics see the same
## model, so the plain order (diagnostics, rows out, final run) is the
## same sequence in one round.
function _se_bad_data_rounds!(net::Net, run_cfg, tap_cfg::TapEstimationConfig, tap_released_idxs::Vector{Int}, tap_estimable::Int, tap_staged_start::Bool, logfile::AbstractString, @nospecialize(se_phase))
  tap_rounds = tap_cfg.enabled && tap_estimable > 0
  # rounds of the sequence (1 without a released tap)
  rounds = 0
  staged = false
  # the converged pre-solve of the round, the reference of the acceptance
  # rule below
  staged_pre = nothing
  # the start the staged pre-solve of the last round took (SEResult.start,
  # :none without a staged start)
  staged_start_from = :none
  tap_fallback_used = false
  tap_run1 = (:none, 0)
  # every row the run removed, with the trace row of the pass that removed it
  eliminations = NamedTuple[]
  diag = nothing
  res = nothing
  if !tap_rounds
    diag = with_sparlectra_config(() -> runse_diagnostics(net), run_cfg)
    # eliminated rows leave the FINAL run: the diagnostics identified them
    # on its own copy, so deactivate them here or the headline J would keep
    # carrying rows the workflow already removed (seen: an injected 10-sigma
    # error, eliminated with stop reason consistent, still pushed the
    # headline to J = 149 at dof 42 instead of J ~ dof)
    for e in diag.eliminations
      (1 <= e.measurement_index <= length(net.measurements)) || continue
      net.measurements[e.measurement_index] = _set_measurement_active(net.measurements[e.measurement_index], false)
      push!(eliminations, merge(e, (round = 1,)))
    end
    se_phase("state_estimation")
    res = with_sparlectra_config(() -> runse!(net), run_cfg)
    rounds = 1
  else
    base_cfg = run_cfg
    max_elim_run = base_cfg.state_estimation.max_eliminations
    # every round starts from the same state and model: the start voltages,
    # the model data a converged run writes back (update_shunts,
    # update_taps) and the transferred step of every released branch
    v_start0 = [(nd._vm_pu, nd._va_deg) for nd in net.nodeVec]
    model_start0 = _se_writeback_snapshot(net, base_cfg.state_estimation)
    steps_start0 = [(k, deepcopy(net.branchVec[k])) for k in tap_released_idxs]
    tap_snap = _tap_release_snapshot(net)
    diag_first = nothing
    reverted = nothing
    tap_round = 0
    while true
      tap_round += 1
      if tap_round > 1
        for (i, nd) in enumerate(net.nodeVec)
          nd._vm_pu, nd._va_deg = v_start0[i]
        end
        _se_writeback_restore!(net, model_start0)
        for (k, br) in steps_start0
          _copy_tap_position!(net.branchVec[k], br)
        end
      end
      run_cfg = base_cfg
      rtag = tap_round == 1 ? "" : string(" (round ", tap_round, ")")
      # Staged start for the mass release. From flat voltages the released
      # tap states take large first steps, the per-iteration tap step limit
      # then scales the WHOLE step down, and the estimate crawls:
      # case1354pegase with 240 released transformers needed 63 iterations
      # on a generated measurement set with four deviated taps, so a limit
      # of 50 never reached it and the run reported J/dof 16 instead of
      # 1.02. One solve with the taps frozen first puts the voltages next to
      # the solution; released from there, the same estimate needed 4
      # iterations and reached the same J.
      staged = false
      staged_pre = nothing
      # the staged pre-solve that did not converge (also not from the case
      # voltages of the start ladder): no released run is started on it
      staged_failed = nothing
      se_phase("state_estimation")
      if tap_staged_start
        _freeze_all_tap_estimation!(net)
        v_start = [(nd._vm_pu, nd._va_deg) for nd in net.nodeVec]
        pre = with_sparlectra_config(() -> runse!(net), run_cfg)
        _restore_tap_release!(net, tap_snap)
        # the start ladder of runse! names the start the pre-solve took
        staged_start_from = pre.start
        if pre.converged
          staged = true
          staged_pre = pre
          run_cfg = _sparlectra_config_with(run_cfg; state_estimation = _se_config_with(run_cfg.state_estimation; flatstart = false))
        else
          # the pre-solve wrote its diverged iterate into the net; the
          # fallback below starts from the start state of the round
          for (i, nd) in enumerate(net.nodeVec)
            nd._vm_pu, nd._va_deg = v_start[i]
          end
          staged_failed = pre
        end
        open(logfile, "a") do io
          pre.start === :stored_after_flat && println(io, "tap estimation", rtag, ": staged start: ", _se_start_ladder_line(nothing))
          if pre.converged
            println(io, "tap estimation", rtag, ": staged start, one solve with the ", tap_estimable, " released tap(s) frozen (", pre.iterations, " iteration(s)", pre.start === :stored_after_flat ? " from the case voltages" : "", ", J = ", round(pre.objectiveJ; sigdigits = 6), ", dof ", pre.dof, "), then the taps are estimated from that state")
          else
            println(io, "tap estimation", rtag, ": the staged start did not converge with the taps frozen (", pre.iterations, " iteration(s)", pre.start === :stored_after_flat ? ", also not from the case voltages" : "", "); the released run is not started")
          end
        end
      end
      res = if staged_failed === nothing
        with_sparlectra_config(() -> runse!(net), run_cfg)
      elseif !tap_cfg.fallback
        # tap_estimation.fallback = false: no repeat, the pre-solve that
        # did not converge is the result and the run fails below
        staged_failed
      else
        nothing
      end
      # The fallback of the run with tap estimation (one logic with runse!):
      # when the run fell back, the result of this round is the run without
      # tap estimation, so the diagnostics, the metadata and run.log describe
      # the run that is the result. Its tap rows carry reason
      # fallback_not_converged (fallback_worse_fit for the second trigger
      # below).
      tap_fallback_used = false
      tap_run1 = (:none, 0)
      if res === nothing
        # The staged pre-solve with the taps frozen did not converge, also
        # not from the case voltages: a released run from that state would
        # run its whole iteration budget with the tap states on top (seen on
        # case6495rte, where it then sat in the rank test for hours). The
        # result is the run without tap estimation, from the start state
        # and model of the round, reason staged_not_converged.
        tap_fallback_used = true
        tap_run1 = (:staged_not_converged, staged_failed.iterations)
        run_cfg = _sparlectra_config_with(base_cfg; state_estimation = _se_config_with(base_cfg.state_estimation; tap_estimation = _tap_config_off(tap_cfg)))
        for (i, nd) in enumerate(net.nodeVec)
          nd._vm_pu, nd._va_deg = v_start0[i]
        end
        _se_writeback_restore!(net, model_start0)
        res = _se_tap_fallback_result(with_sparlectra_config(() -> runse!(net), run_cfg), tap_run1...; net = net, measurements = Measurement[m for m in net.measurements], excl = _tap_exclusion_map(net, tap_cfg), update_net = run_cfg.state_estimation.update_net)
        open(logfile, "a") do io
          println(io, _se_tap_fallback_line2(res))
        end
      end
      res_fell_back = !tap_fallback_used && res.tapFixation !== nothing && res.tapFixation.fallback_used
      # Second trigger, the acceptance rule of the staged start
      # (worse_fit): the staged
      # estimate converges where the cold one crawled, but on a set that
      # cannot pin its taps it converges to positions that round badly. The
      # FIXED estimate is kept only if it fits at least as well as the
      # transferred steps in the pre-solve of the same round (same rows).
      # "Larger" carries a relative margin of 1e-6: a tap fixed back onto
      # its transferred step solves the same problem as the pre-solve, and
      # the two J then differ in the last digits only (seen: 49.4686 against
      # 49.4686 rejected without the margin). Its fallback is the same as
      # for a failed run 1, reason worse_fit, rows fallback_worse_fit;
      # fallback = false switches it off.
      tap_worse_fit = !tap_fallback_used && !res_fell_back && staged_pre !== nothing && res.converged && res.tapFixation !== nothing && res.objectiveJ > staged_pre.objectiveJ * (1.0 + 1.0e-6)
      if tap_cfg.fallback && (res_fell_back || tap_worse_fit)
        tap_fallback_used = true
        tap_run1 = res_fell_back ? (res.tapFixation.fallback_reason, res.tapFixation.iterations_estimation) : (:worse_fit, _se_reported_iterations(res))
        run1_res = res
        # the repeat is the run without tap estimation of this round: the
        # configuration of the round without the warm start of the staged
        # pre-solve, from the start voltages and model data of the round.
        # runse! writes its state into the net whether or not it converged
        # (the svedala set left bus angles over -178 to +167 degrees where the
        # power flow had -26 to +40), and the pre-solve and a converged run 1
        # write back (update_shunts, update_taps). Starting from the
        # pre-solve state instead gave a result that differs from the run
        # with tap estimation switched off in the last digits.
        run_cfg = _sparlectra_config_with(base_cfg; state_estimation = _se_config_with(base_cfg.state_estimation; tap_estimation = _tap_config_off(tap_cfg)))
        open(logfile, "a") do io
          println(io, _se_tap_fallback_line1_service(tap_run1..., run1_res, staged_pre))
        end
        for (i, nd) in enumerate(net.nodeVec)
          nd._vm_pu, nd._va_deg = v_start0[i]
        end
        _se_writeback_restore!(net, model_start0)
        # tap_cfg: the settings of run 1, so the excluded and the machine
        # transformers keep their reason in the rows of the fallback run
        res = _se_tap_fallback_result(with_sparlectra_config(() -> runse!(net), run_cfg), tap_run1...; net = net, measurements = Measurement[m for m in net.measurements], excl = _tap_exclusion_map(net, tap_cfg), update_net = run_cfg.state_estimation.update_net)
        open(logfile, "a") do io
          println(io, _se_tap_fallback_line2(res))
        end
      end
      # no end steps without a converged result: the run fails below
      res.converged || break
      # the diagnostics on the end steps of this result: every released
      # branch at its end step (fixed estimate, or the transferred step of
      # the fallback) and frozen, warm from the result state; the steps the
      # run left in the net (update_taps) and the release are put back
      # afterwards
      se_phase("estimation_diagnostics")
      steps_run = [(k, deepcopy(net.branchVec[k])) for k in tap_released_idxs]
      _tap_stamp_end_steps!(net, res, run_cfg.state_estimation.tap_estimation.impedance_correction)
      _freeze_all_tap_estimation!(net)
      budget = clamp(max_elim_run - length(eliminations), 0, 1)
      diag_cfg = _sparlectra_config_with(run_cfg; state_estimation = _se_config_with(run_cfg.state_estimation; max_eliminations = budget, flatstart = false))
      d = with_sparlectra_config(() -> runse_diagnostics(net), diag_cfg)
      for (k, br) in steps_run
        _copy_tap_position!(net.branchVec[k], br)
      end
      _restore_tap_release!(net, tap_snap)
      diag_first === nothing && (diag_first = d)
      d.reverted === nothing || (reverted = d.reverted)
      diag = d
      for e in d.eliminations
        (1 <= e.measurement_index <= length(net.measurements)) || continue
        net.measurements[e.measurement_index] = _set_measurement_active(net.measurements[e.measurement_index], false)
        push!(eliminations, merge(e, (elimination = length(eliminations) + 1, round = tap_round)))
      end
      isempty(d.eliminations) && break
      open(logfile, "a") do io
        println(io, "bad data", rtag, ": ", join((String(e.id) for e in d.eliminations), ", "), " eliminated on the end steps of the result run (", tap_fallback_used ? "fallback without tap estimation" : "taps estimated and fixed",
                ", J = ", round(res.objectiveJ; sigdigits = 6), ", dof ", res.dof, "); the result run is repeated with released taps on the cleaned set")
      end
    end
    # the report of the whole sequence: the first pass shows where the bad
    # rows were (every row active), the last one describes the result
    if diag !== nothing
      diag = (diagnostics = diag_first.diagnostics, rerun = diag_first.rerun, eliminations = eliminations, reverted = reverted, stop_reason = diag.stop_reason,
              final_diagnostics = diag.final_diagnostics, topology_findings = diag.topology_findings, range_check = diag_first.range_check)
    end
    rounds = tap_round
  end
  return (res = res, diag = diag, eliminations = eliminations, staged = staged, staged_start_from = staged_start_from, tap_fallback_used = tap_fallback_used, tap_run1 = tap_run1, rounds = rounds, tap_rounds = tap_rounds, run_cfg = run_cfg)
end

"""
    _run_state_estimation_service(case_path, config_file, output_dir, run_id, measurement_file; kwargs...) -> SparlectraApiResult

Service backend of the Web UI "Run state estimation" action (SE phase 5).
Builds the net through the shared import paths, reads `measurement_file`
(a measurement table; without one, the set an SCF case carries),
evaluates global observability (structural
islands plus FD-aware rank, phase 4), runs the bad-data diagnostics
(`runse_diagnostics`, sequential elimination on the configured budget) and
the final `runse!(updateNet = true)`, and writes the SE artifacts:
`measurements.csv` (the estimated set as a measurement table, comments
carried over; the run estimates the set this file reads back to),
`measurements_verdict.csv` (the verdict per table row), `se_diagnostics.md`, `se_view.md`,
`shunt_estimates.csv` (when shunts were released), and `se_state.csv` (the
chain anchor a later SE-started power flow consumes via `readSEStateCSV!`).

Options: `max_iter`, `tol`, `flatstart`, `robust`, `max_eliminations`,
`update_shunts`, `report_correlation` mirror the estimator keywords.
With `state_estimation.tap_estimation.enabled` (default `true`, the Web UI
"estimate taps" box) every declared tap changer of every in-service
transformer is released (`setTapEstimation!` with the mode its changers
give, `:ratio`, `:pst` or `:both`) before the solve; a transformer that only
carries a ratio has no tap changer, no state and no row; the estimator then fixes each
regulator to its nearest physical step and reports J before versus after
the fixation (`se_tap_estimates.csv`, one row per regulator with `status`
and `reason`). With `false` every changer keeps its transferred step and is
listed with reason `disabled`. Transformers named in
`state_estimation.tap_estimation.exclude` keep their step (reason
`excluded`); machine transformers are calculated after the estimation
(status `calculated`, reason `machine_transformer`) unless
`include_machine_transformers` is set. run.log has one line per such
transformer and one warning per exclude entry without a transformer.
That release is staged: one solve with the taps frozen first, the estimate
from that state. When the run with tap estimation does not converge and
`state_estimation.tap_estimation.fallback` is on, the run is repeated
without tap estimation and that run is the result of the round
(`se_tap_fallback_used`, reason `fallback_not_converged`). The same
fallback runs when the fixed estimate fits worse than the staged pre-solve
(acceptance rule `worse_fit`, reason `fallback_worse_fit`). With released
taps the bad-data diagnostics run on the end steps of that result (every
tap frozen at its fixed step, or at its transferred step after the
fallback); when they eliminate a row (one per round), the next round
repeats the result run with released taps on the cleaned set, until a round
eliminates nothing. `max_eliminations` counts over all rounds
(`se_bad_data_rounds`). The
internal keyword `tap_staged_start = false` keeps the release cold
(regression fixtures only).

Failure reasons: `se_unsupported_format`, `import_error`,
`invalid_measurements` (missing/unreadable/rejected file),
`se_not_observable` (observability quality `:not_observable`),
`se_not_converged`, plus the shared config failure.
"""
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
  # internal: false keeps the mass release COLD (no frozen pre-solve, no
  # acceptance rule). No request field sets it; it exists for the
  # regression fixtures that pin what the diagnostics do on a cold tap run
  # (a reverted elimination is reachable only there) and for comparisons.
  tap_staged_start::Bool = true,
)::SparlectraApiResult
  # the callback and the overrides arrive in a different type from every call
  # site; behind `_SeRunOptions` the 700-line body below has one
  # specialization (the same arrangement as `_run_sparlectra_api`)
  options = _SeRunOptions(max_iter, tol, flatstart, robust, max_eliminations, report_correlation, k_eliminate, robust_mode, robust_k1,
    robust_k2, k_suppress, suppression_sigma, phase_callback, config_overrides, tap_staged_start)
  return _run_state_estimation_service_body(String(case_path), String(config_file), String(output_dir), run_id, String(measurement_file),
    update_shunts, case_format, options)
end

# the optional solver settings keep their Union types as fields: a positional
# argument of Union type is specialized per concrete value type (`nothing` or
# a number, eight variants in one test group), a struct field is not
struct _SeRunOptions
  max_iter::Union{Nothing,Int}
  tol::Union{Nothing,Float64}
  flatstart::Union{Nothing,Bool}
  robust::Union{Nothing,Bool}
  max_eliminations::Union{Nothing,Int}
  report_correlation::Union{Nothing,Bool}
  k_eliminate::Union{Nothing,Float64}
  robust_mode::Union{Nothing,Symbol}
  robust_k1::Union{Nothing,Float64}
  robust_k2::Union{Nothing,Float64}
  k_suppress::Union{Nothing,Float64}
  suppression_sigma::Union{Nothing,Float64}
  phase_callback::Any
  config_overrides::Any
  tap_staged_start::Bool
end

function _run_state_estimation_service_body(
  case_path::String,
  config_file::String,
  output_dir::String,
  run_id::String,
  measurement_file::String,
  update_shunts::Bool,
  case_format::Symbol,
  options::_SeRunOptions,
)::SparlectraApiResult
  max_iter = options.max_iter
  tol = options.tol
  flatstart = options.flatstart
  robust = options.robust
  max_eliminations = options.max_eliminations
  report_correlation = options.report_correlation
  k_eliminate = options.k_eliminate
  robust_mode = options.robust_mode
  robust_k1 = options.robust_k1
  robust_k2 = options.robust_k2
  k_suppress = options.k_suppress
  suppression_sigma = options.suppression_sigma
  phase_callback = options.phase_callback
  tap_staged_start = options.tap_staged_start
  config_overrides = options.config_overrides
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

  # the same precedence the power-flow path uses (resolve_config):
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
  # The two knees belong to the staged mode only: down-weighting starts at
  # k1 and is complete at k2, so k1 above k2 describes no curve. In every
  # other mode the pair is not read; a stored pair that is out of order
  # there must not stop the run, it is replaced by the configured pair and
  # the log says so.
  knees_ok = robust_k1 > 0.0 && robust_k2 >= robust_k1
  if robust_mode === :staged
    knees_ok || return _api_failure("invalid_request", "Staged down-weighting needs its start (\"down-weight from |r|/sigma\", se_robust_k1) positive and not above its end (\"full down-weight at |r|/sigma\", se_robust_k2): got k1=$(robust_k1), k2=$(robust_k2). Lower the first value or raise the second.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  elseif !knees_ok
    open(logfile, "a") do io
      println(io, "note: se_robust_k1=", robust_k1, " and se_robust_k2=", robust_k2, " are out of order; they are not used with down-weighting ", robust_mode, " and were replaced by the configured pair (", se_cfg.robust_k1, ", ", se_cfg.robust_k2, ")")
    end
    robust_k1 = se_cfg.robust_k1
    robust_k2 = se_cfg.robust_k2
  end
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
  from_table = false
  # rows state_estimation.zero_injection added to the set (all paths)
  zi_rows = 0
  summary = if from_case_file
    n_case_rows = length(net.measurements)
    open(logfile, "a") do io
      println(io, "measurements: ", n_case_rows, " row(s) taken from the case file itself (no separate measurement set given)")
    end
    zi_rows = _se_apply_zero_injection!(net, config, logfile)
    case_index = _se_measurement_index(imported, format)
    artifact = _se_write_measurement_artifact!(net, case_index, joinpath(output_dir, "measurements.csv"), [string("case: ", basename(String(case_path))), string("measurement set of state-estimation run ", run_id, ", taken from the case file; verdict per row in measurements_verdict.csv")], logfile, config.state_estimation.zero_injection_sigma_mva)
    (total = n_case_rows,)
  else
    if !isfile(measurement_file)
      hint = format === :scf ? " The case file carries no measurements either: generate a set, or export the case with its measurements." : ""
      return _api_failure("invalid_measurements", "Measurement file not found: $(measurement_file)." * hint, run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
    # case-binding gate (measurement table and CSV v1 alike): a set recorded
    # for a different case would fail rows
    # deep in the reader ("bus X not found"); refuse it up front with the
    # actual cause. Sets without a binding (older files) pass with a log note.
    set_binding = _webui_measurement_set_binding(String(measurement_file))
    bound_case = set_binding.case
    # A case exported to SCF keeps the model of its source case, and it records
    # WHICH source it came from. A set bound to that source therefore fits this
    # run; without this the export silently invalidates every set the user
    # already generated.
    source_of_case = format === :scf ? _scf_source_reference(case_path) : ""
    if !isempty(bound_case) && bound_case != basename(String(case_path)) && bound_case != source_of_case
      return _api_failure("invalid_measurements", "Measurement set $(basename(String(measurement_file))) is bound to case $(bound_case) (recorded in the file), but this run uses $(basename(String(case_path))). Pick the set generated for this case, or regenerate one.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
    # A table bound by the CGMES model id only (the converter writes it) is
    # not refused on a mismatch: the import report counts the rows the case
    # cannot resolve, and that count is the actual test of the binding.
    if isempty(bound_case) && !isempty(set_binding.models)
      case_models = format === :cgmes ? _webui_case_cgmes_model_ids(String(case_path)) : String[]
      matched = any(in(case_models), set_binding.models)
      open(logfile, "a") do io
        println(io, "note: the measurement set is bound by CGMES model ", join(set_binding.models, ", "), matched ? "; this is the EQ model of the case" : string("; the EQ model of the case is ", isempty(case_models) ? "unknown" : join(case_models, ", "), ", rows the case cannot resolve are counted in the import report"))
      end
    elseif isempty(bound_case)
      open(logfile, "a") do io
        println(io, "note: the measurement set carries no case binding (# case: comment); older set, association checked by name only")
      end
    end
    # a table that lost its version line goes the table way too, so the
    # reader names the missing line (the CSV v1 reader would complain about
    # an unknown CSV version instead)
    if isMeasurementTable(String(measurement_file)) || _measurement_table_without_version(String(measurement_file))
      from_table = true
      table_summary, table_error = _se_load_measurement_table!(net, imported, format, String(measurement_file), config, logfile, String(output_dir))
      table_error === nothing || return _api_failure("invalid_measurements", table_error; run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
      zi_rows = table_summary.zero_injection
      # the comment lines of the table (case binding, generator provenance,
      # taps table, truth values) go into the artifact
      table_comments, _ = _measurement_csv_comments(measurement_file)
      # the run lines of an artifact given as the set again are replaced,
      # not stacked
      filter!(c -> !startswith(c, "measurement set of state-estimation run ") && !startswith(c, "not in the table: "), table_comments)
      tset = table_summary.table.set
      push!(table_comments, string("measurement set of state-estimation run ", run_id, ": set ", isempty(tset.timestamp) ? "(no timestamp)" : tset.timestamp, " of ", basename(String(measurement_file)), ", ", sum(length, values(tset.sources); init = 0), " of ", tset.counts[:rows], " row(s) used (the others: measurement_import.txt); verdict per row in measurements_verdict.csv"))
      artifact = _se_write_measurement_artifact!(net, table_summary.index, joinpath(output_dir, "measurements.csv"), table_comments, logfile, config.state_estimation.zero_injection_sigma_mva; table = table_summary.table)
      table_summary
    else
      # the measurement CSV v1 reader was removed in 0.50.0: a v1 file is
      # refused with the converter named, anything else as no table
      name = basename(String(measurement_file))
      msg = _is_measurement_csv_v1(String(measurement_file)) ? "$(name): $(MEASUREMENT_CSV_V1_REMOVED)" : "$(name) is no measurement table (the first line must be $(MEASUREMENT_TABLE_VERSION_LINE))."
      return _api_failure("invalid_measurements", msg; run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
  end

  # where the sigma of every measurement came from, per voltage level; only
  # a table row can be without its own sigma
  sigma_level_report = _se_sigma_level_report(net, config, from_table ? summary.sigma_origin : Dict{String,Symbol}())
  base_metadata["se_sigma_by_voltage_level"] = sigma_level_report.rows
  isempty(sigma_level_report.unused) || (base_metadata["se_sigma_levels_not_in_case"] = sigma_level_report.unused)
  open(logfile, "a") do io
    for line in _se_sigma_level_lines(sigma_level_report)
      println(io, line)
    end
  end

  # Carry a duplicate finding into the run: it is the one property of a set that
  # explains a large J without any single measurement looking wrong.
  # Checked on the loaded set, not on the reader result: measurements carried
  # inside a case file take the other branch above and can be doubled just the
  # same.
  # A measurement table keeps several keys at one location as redundant
  # meters on purpose (the import reports the same key twice as duplicate),
  # so for a table the repeat is stated as a fact, not warned about.
  duplicate_rows = sum(values(_duplicate_measured_quantities(net.measurements)); init = 0)
  if duplicate_rows > 0 && from_table
    open(logfile, "a") do io
      println(io, "note: $(duplicate_rows) of $(length(net.measurements)) measurements are further readings of an already measured quantity (redundant meters in the table)")
    end
  elseif duplicate_rows > 0
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
    # an abort is no rejected measurement set
    err isa PowerFlowAborted && rethrow()
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

  # released transformer taps: every declared tap changer of an in-service
  # transformer becomes an estimation state (the mode from its changer set,
  # _tap_release_mode); a transformer that only carries a ratio has an
  # empty set, gets no state and no tap row. Released AFTER the
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
  # A set that DOCUMENTS its tap deviations (the generator moved taps) can
  # only be estimated with tap estimation on; with the switch off the
  # estimator cannot absorb them and J stays large for a reason the log
  # names here.
  tap_enabled_cfg = config.state_estimation.tap_estimation.enabled
  if !isempty(set_tap_devs)
    base_metadata["se_set_tap_deviation"] = true
    tap_enabled_cfg || open(logfile, "a") do io
      println(io, "note: the measurement set documents a tap deviation (", join(("$(d.steps) step(s) on branch $(d.branch)" for d in set_tap_devs), ", "), ") but tap estimation is switched off (state_estimation.tap_estimation.enabled = false); the estimator cannot absorb it (expect band :high clustered at the transformer)")
    end
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

  # the tap estimation switch of this run: the configuration of the case
  # with the request overrides folded in (the Web UI "estimate taps" box is
  # state_estimation.tap_estimation.enabled)
  tap_cfg = run_cfg.state_estimation.tap_estimation
  tap_enabled = tap_cfg.enabled
  tap_released = 0
  tap_released_idxs = Int[]
  for (k, br) in enumerate(net.branchVec)
    # the release mode from the changer set of the branch (DECLARED changers
    # only; the constructor gives every transformer a default ratio grid and
    # phase fields, and mass-releasing those added meaningless states
    # everywhere, 8617 on case6495rte); no direction: it comes from the
    # changer data.
    # Every changer is released; runse! keeps the excluded and the machine
    # transformers out by the configuration and reports them with that reason
    mode = _tap_release_mode(br)
    (mode !== nothing && _branch_terminal_state(br) == :closed) || continue
    setTapEstimation!(net; trafo = k, mode = mode)
    tap_released += 1
    push!(tap_released_idxs, k)
  end
  # which released transformers stay out (tap_estimation.exclude, machine
  # transformers): one run.log line each below; the exclude entries without
  # a transformer are named once here and left out of the runs, so the
  # many estimator calls of this workflow do not warn again
  # (with the switch off every changer says disabled: no reasons)
  tap_excl = _tap_exclusions(net, tap_cfg; ks = tap_enabled ? tap_released_idxs : Int[])
  if !isempty(tap_excl.unknown)
    tap_cfg = TapEstimationConfig(enabled = tap_cfg.enabled, fallback = tap_cfg.fallback, impedance_correction = tap_cfg.impedance_correction, exclude = String[e for e in tap_cfg.exclude if !(e in tap_excl.unknown)], include_machine_transformers = tap_cfg.include_machine_transformers)
    run_cfg = _sparlectra_config_with(run_cfg; state_estimation = _se_config_with(run_cfg.state_estimation; tap_estimation = tap_cfg))
  end
  tap_n_excluded = count(==(:excluded), values(tap_excl.reasons))
  tap_n_machine = count(==(:machine_transformer), values(tap_excl.reasons))
  # the released transformers the estimation can actually see
  tap_estimable = count(k -> !haskey(tap_excl.reasons, k), tap_released_idxs)
  open(logfile, "a") do io
    for e in tap_excl.unknown
      println(io, "WARNING: state_estimation.tap_estimation.exclude names '", e, "', which is no transformer of this case; the entry has no effect")
    end
    for k in sort!(collect(keys(tap_excl.reasons)))
      nm = getCompName(net.branchVec[k].comp)
      if tap_excl.reasons[k] === :excluded
        println(io, "tap estimation: transformer ", nm, " (branch ", k, ") excluded (state_estimation.tap_estimation.exclude), keeps its transferred step, not calculated")
      else
        println(io, "tap estimation: transformer ", nm, " (branch ", k, ") is a machine transformer, not estimated (include_machine_transformers = false); its step is calculated after the estimation")
      end
    end
    if !tap_enabled
      tap_released > 0 && println(io, "tap estimation: switched off (state_estimation.tap_estimation.enabled = false); ", tap_released, " transformer(s) keep their transferred step")
    else
      # loud warning when the set's documented deviation sits on a
      # transformer the run keeps out: the estimation can never absorb it
      # and J stays large for a reason the user cannot see otherwise
      for d in set_tap_devs
        r = get(tap_excl.reasons, d.branch, :none)
        r === :none && continue
        println(io, "WARNING: the set documents a tap deviation of $(d.steps) step(s) on branch $(d.branch), which is kept out of the tap estimation (", r, "); the estimation cannot absorb it (", r === :excluded ? "remove it from state_estimation.tap_estimation.exclude" : "set state_estimation.tap_estimation.include_machine_transformers", ", or regenerate the set)")
      end
      # a documented deviation on a branch the mass release cannot release
      # (a line, a transformer without changer data, an open branch) is
      # just as unabsorbable; named, or the set's provenance is lost here
      for d in set_tap_devs
        d.branch in tap_released_idxs && continue
        println(io, "note: the measurement set documents a tap deviation of $(d.steps) step(s) on branch $(d.branch), which has no tap changer the run can release; the estimator cannot absorb it (expect band :high clustered at the branch)")
      end
    end
  end

  # range check of the voltage measurements: every
  # estimator call below runs it on net.measurements by itself; the rows
  # are taken here once, with the settings of this run, for run.log, the
  # metadata and the verdict file (same function, same rows)
  range_rows = with_sparlectra_config(() -> last(_se_range_check(net, Measurement[m for m in net.measurements])), run_cfg)
  isempty(range_rows) || open(logfile, "a") do io
    for line in _se_range_check_lines(net, range_rows)
      println(io, line)
    end
  end
  base_metadata["se_range_check_substituted"] = count(r -> r.reason === :range_check_substituted, range_rows)
  base_metadata["se_range_check_deenergized"] = count(r -> r.reason === :deenergized, range_rows)
  # a non-finite measured value goes into the JSON as text
  base_metadata["se_range_check"] = [Dict{String,Any}("id" => r.id, "bus" => r.bus, "value" => isfinite(r.value) ? r.value : string(r.value), "reason" => String(r.reason),
                                                      "substitute" => isfinite(r.substitute) ? r.substitute : nothing, "source" => String(r.source), "sigma_substitute" => isfinite(r.sigma_substitute) ? r.sigma_substitute : nothing,
                                                      "deactivated" => length(r.deactivated)) for r in range_rows]

  # diagnostics (report) plus the final state-writing run (chain anchor),
  # with released taps in rounds on the end steps (_se_bad_data_rounds!)
  se_phase("estimation_diagnostics")
  bd = _se_bad_data_rounds!(net, run_cfg, tap_cfg, tap_released_idxs, tap_estimable, tap_staged_start, logfile, se_phase)
  res, diag, eliminations, staged, tap_fallback_used, tap_run1, tap_rounds, run_cfg = bd.res, bd.diag, bd.eliminations, bd.staged, bd.tap_fallback_used, bd.tap_run1, bd.tap_rounds, bd.run_cfg
  base_metadata["se_bad_data_rounds"] = bd.rounds
  base_metadata["se_tap_staged_start"] = staged
  base_metadata["se_tap_fallback_used"] = tap_fallback_used
  base_metadata["se_tap_fallback_reason"] = String(tap_run1[1])
  # where the estimation started (flat, stored, stored_after_flat: the
  # start ladder): after a converged staged start the start of that
  # pre-solve, the released run then starts from its state; else the start
  # of the result run. The staged pre-solve on its own (none without one).
  base_metadata["se_start"] = String(staged ? bd.staged_start_from : res.start)
  base_metadata["se_tap_staged_start_from"] = String(bd.staged_start_from)
  res.start === :stored_after_flat && open(logfile, "a") do io
    println(io, "start: ", _se_start_ladder_line(nothing), "; the result run starts from the case voltages")
  end
  if !res.converged
    _write_service_performance_log!(output_dir, phase_recorder, se_total_start; headline = _SE_PERF_HEADLINE, status = "failed", label = "state-estimation")
    run1_text = tap_run1[1] === :staged_not_converged ? string("the staged start with the taps frozen (", tap_run1[2], " iterations)") : string("run 1 with tap estimation (", tap_run1[1], ", ", tap_run1[2], " iterations)")
    failure_text = tap_fallback_used ? string("State estimation did not converge: ", run1_text, ", run 2 without tap estimation (", _se_tap_run2_reason(res), ", ", res.iterations, " iterations).") : "State estimation did not converge within $(max_iter) iterations."
    return _api_failure("se_not_converged", failure_text, run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  se_phase("postprocessing_result")

  # artifacts
  open(joinpath(output_dir, "se_diagnostics.md"), "w") do io
    # same statement the result page carries: the tap positions below are
    # model values, so the J of this run measures THEM
    tap_fallback_used && println(io, "\n> **Tap estimation fallback.** ", _se_tap_fallback_note(base_metadata), " The diagnostics below describe that run without tap estimation, the same one the status line reports.\n")
    # The numbers the status line and run.log report come from the FINAL
    # run (res); the diagnostics pass below ran on every active row before
    # any elimination, so its own J counts rows the reported run no longer
    # carries (seen on the MicroGrid: J 3555 at dof 65 in the
    # pass, J 74 at dof 64 reported, one row eliminated in between). The
    # report states the reported numbers first, so it never contradicts the
    # page it belongs to.
    println(io, "\n## Reported result\n")
    # the run, the case and the set this report belongs to, so the file
    # reads on its own (run id, case and set in
    # the report; J/dof next to J, the page's leading number)
    println(io, "- **Run:** ", run_id)
    println(io, "- **Case:** ", basename(case_path))
    println(io, "- **Measurement set:** ", from_case_file ? "carried by the case file" : basename(measurement_file))
    reported_it = _se_reported_iterations(res)
    println(io, "- **Converged:** ", res.converged, " in ", reported_it, " iteration(s)", reported_it != res.iterations ? " (max over the passes; final solve $(res.iterations))" : "")
    println(io, "- **J/dof:** ", round(res.objectiveJ / max(res.dof, 1); digits = 2))
    println(io, "- **Objective J:** ", round(res.objectiveJ; digits = 6), " (dof=", res.dof, ")", tap_fallback_used ? ", without tap estimation (fallback)" : "")
    n_elim = length(eliminations)
    n_rounds = maximum((e.round for e in eliminations); init = 1)
    if n_rounds > 1
      # rows of later rounds were found on the end steps of the repeated
      # result run, not in the pass printed below
      println(io, "- **Eliminated rows:** ", n_elim, " (inactive in the reported run), found on the end steps of the result run in ", n_rounds, " rounds: ",
              join((string("round ", k, ": ", join(("`" * String(e.id) * "`" for e in eliminations if e.round == k), ", ")) for k in 1:n_rounds if any(e -> e.round == k, eliminations)), "; "))
    else
      n_elim > 0 && println(io, "- **Eliminated rows:** ", n_elim, " (identified by the diagnostics pass below; inactive in the reported run, still counted in the pass's own J)")
    end
    (hasproperty(diag, :reverted) && diag.reverted !== nothing) && println(io, "- **Not eliminated:** `", diag.reverted.id, "` (the control solve without it did not converge; the row stays active, see the pass below)")
    println(io, "\nThese are the numbers of the status line and of `run.log`.\n")
    println(io, "## Measurement sigma by voltage level\n")
    sigma_lines = _se_sigma_level_lines(sigma_level_report)
    println(io, sigma_lines[1], "\n")
    for line in sigma_lines[2:end]
      println(io, "- ", strip(line))
    end
    println(io)
    println(io, "## Diagnostics pass (every row active, before any elimination)\n")
    println(io, "The pass below ran on all rows and identified the eliminated ones; its J is therefore not the reported J. It is where a bad row shows.", tap_rounds ? " With released taps it ran on the end steps of the first result run, every tap frozen there." : "", "\n")
    print_se_diagnostics(io, diag; topN = 15, format = :markdown)
    # the critical rows by name (a detected critical measurement is stated, never left to be read off the wii)
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
  # bus links have no Ybus entry: their flows come from the nodal balances
  # of the link clusters (injection specs of the members, shunt draw and
  # branch flows at the estimated state) plus the ring circulations the
  # final run estimated from the link flow measurements, which are WLS rows
  # (a coupler ring has no unique split without them). nothing when no link
  # is closed.
  link_rows = calcLinkFlowsSE!(net)
  view = se_view(net)
  open(joinpath(output_dir, "se_view.md"), "w") do io
    print_se_view(io, view; format = :markdown)
  end
  # One format for every diagnostic CSV this SE run writes (issue #376),
  # same config key the PowerFlow API uses for its own artifacts.
  se_csv_format = _resolve_detailed_csv_format(String(config.output.csv_format))
  se_csv_delim = se_csv_format.delimiter
  n_link_measured = 0
  if link_rows !== nothing
    nbyL = _bus_name_by_idx(net)
    open(joinpath(output_dir, "se_link_flows.csv"), "w") do io
      println(io, join(("link", "name", "from_bus", "to_bus", "p_mw", "q_mvar", "source", "p_residual_mw", "q_residual_mvar"), se_csv_delim))
      cfL(s) = _csv_field(string(s), se_csv_delim, se_csv_format)
      fmtL(v) = isnan(v) ? "" : _format_csv_number(round(v; sigdigits = 8), se_csv_format)
      for r in link_rows
        l = net.linkVec[r.linkIdx]
        r.source === :kcl || (n_link_measured += 1)
        println(io, join((r.linkIdx, cfL(r.name), cfL(get(nbyL, l.fromBus, string(l.fromBus))), cfL(get(nbyL, l.toBus, string(l.toBus))), fmtL(r.pFlow_MW), fmtL(r.qFlow_MVar), r.source, fmtL(r.p_meas_residual), fmtL(r.q_meas_residual)), se_csv_delim))
      end
    end
  end
  # zero_injection = constraint: the Lagrange multiplier of every
  # constraint row, one line per bus (passive bus or fused coupler bus)
  if res.zeroInjectionMultipliers !== nothing
    open(joinpath(output_dir, "se_zero_injection_multipliers.csv"), "w") do io
      println(io, join(("bus", "kind", "lambda_p_per_mw", "lambda_q_per_mvar"), se_csv_delim))
      fmtZ(v) = isnan(v) ? "" : _format_csv_number(round(v; sigdigits = 8), se_csv_format)
      for z in res.zeroInjectionMultipliers
        println(io, join((_csv_field(String(z.bus), se_csv_delim, se_csv_format), z.kind, fmtZ(z.lambda_p), fmtZ(z.lambda_q)), se_csv_delim))
      end
    end
  end
  if res.shuntEstimates !== nothing
    open(joinpath(output_dir, "shunt_estimates.csv"), "w") do io
      println(io, join(("bus", "bus_name", "B_model_pu", "B_est_pu", "delta_pu", "frozen"), se_csv_delim))
      for r in res.shuntEstimates
        println(io, join((r.busIdx, _csv_field(String(r.busName), se_csv_delim, se_csv_format), _format_csv_number(Float64(r.B_model), se_csv_format), _format_csv_number(Float64(r.B_est), se_csv_format), _format_csv_number(Float64(r.delta), se_csv_format), r.frozen), se_csv_delim))
      end
    end
  end
  # machine transformers are CALCULATED, never estimated: runse! keeps them
  # out of the estimation and back-calculates the step from the AVR
  # setpoint, the dispatch P and the MEASURED machine Q (the Qinj telemetry
  # at the machine bus when the set carries it), status calculated. One
  # run.log line per machine transformer with the result or the reason
  # (q_not_measured_in_group: no Q row on every bus of the machine node)
  open(logfile, "a") do io
    for t in something(res.tapEstimates, NamedTuple[])
      (t.reason in (:machine_transformer, :q_not_measured_in_group) && t.role === :ratio) || continue
      if t.status === :calculated
        println(io, "machine transformer ", t.name, ": tap CALCULATED (not estimated) from AVR setpoint and machine telemetry: s ", round(t.s_est; digits = 2), " -> step ", _se_tap_step_text(t.step_fixed), " (", t.notice, ")")
      else
        println(io, "machine transformer ", t.name, ": tap back-calculation not possible (", t.notice, ")")
      end
    end
  end

  # rows the robust solve acted on (final state-writing run): replacement
  # stage 3 = SUPPRESSED with the fixed sigma, staged stages 1/2 =
  # down-weighted; measurements_verdict.csv names them per row
  # suspicious: the rows the FINAL diagnostics flag (after every
  # elimination), the count of the verdict file, run.log, the result page
  # and se_suspicious alike; the pass before the elimination lists its own
  # suspects in se_diagnostics.md
  suspicious = diag.final_diagnostics.suspicious_measurements
  suppressedIds = Set{String}(String(r.id) for r in something(res.robustRows, NamedTuple[]) if r.stage == 3)
  downweightedIds = Set{String}(String(r.id) for r in something(res.robustRows, NamedTuple[]) if r.stage in (1, 2))
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
    # a case-file set got the ids of the table in measurements.csv
    v isa Real && (truthById[get(artifact.id_map, String(k), String(k))] = Float64(v))
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
    elimD = Set(t.id for t in eliminations)
    residById = Dict{String,Float64}()
    for row in diag.final_diagnostics.measurement_ranking
      residById[String(row.id)] = row.residual
    end
    open(joinpath(output_dir, "se_deltas.csv"), "w") do io
      println(io, "# sparlectra-se-deltas v1")
      println(io, join(("kind", "id", "type", "bus", "from_bus", "to_bus", "measured", "truth", "estimated", "delta_meas_truth", "delta_est_truth", "delta_est_meas", "sigma", "t_est", "eliminated"), se_csv_delim))
      fmt(v) = _format_csv_number(round(v; sigdigits = 8), se_csv_format)
      cfD(s) = _csv_field(string(s), se_csv_delim, se_csv_format)
      # ZI rows solved as equality constraints are listed as constraints
      zi_constraint = config.state_estimation.zero_injection === :constraint
      # a row the range check replaced was solved as its substitute with the
      # widened sigma: its residual belongs to the substitute, so the
      # estimate is substitute - residual and t_est uses the widened sigma;
      # `measured` keeps the raw reading of the set
      rangeSubD = Dict{Int,Any}(r.measurement_index => r for r in range_rows if r.reason === :range_check_substituted)
      for (mi, m) in enumerate(net.measurements)
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
        rsub = get(rangeSubD, mi, nothing)
        zD = rsub === nothing ? m.value : rsub.substitute
        sigD = rsub === nothing ? m.sigma : rsub.sigma_substitute
        estv = rres === nothing ? nothing : zD - rres
        println(
          io,
          join(
            (zi_constraint && _zero_injection_constraint(m) ? "constraint" : "meas", cfD(m.id), cfD(m.typ), cfD(busD), cfD(fbD), cfD(tbD), fmt(m.value), fmt(tv), estv === nothing ? "" : fmt(estv), fmt(m.value - tv), estv === nothing ? "" : fmt(estv - tv), estv === nothing ? "" : fmt(estv - m.value), fmt(m.sigma), rres === nothing ? "" : fmt(abs(rres) / sigD), m.id in elimD),
            se_csv_delim,
          ),
        )
      end
      devByBranch = Dict{Int,Float64}(d.branch => d.steps for d in set_tap_devs)
      for t in something(res.tapEstimates, NamedTuple[])
        # the documented deviation belongs to the deviation changer of the
        # branch (the generator moves that one), 0 for the others
        set = 1 <= t.branch <= length(net.branchVec) ? _effective_tap_changers(net.branchVec[t.branch]) : nothing
        dev = (set !== nothing && t.regulator == _tap_deviation_changer(set)) ? get(devByBranch, t.branch, 0.0) : 0.0
        # a tap is never measured: `measured` and the deltas against it stay
        # empty rather than repeating the truth and pretending a zero residual
        println(io, join(("tap", cfD(t.name), cfD(t.mode), "", "", "", "", fmt(dev), fmt(t.deviation_steps), "", fmt(t.deviation_steps - dev), "", "", "", false), se_csv_delim))
      end
      deltas_written = true
    end
  end

  if res.tapEstimates !== nothing
    open(joinpath(output_dir, "se_tap_estimates.csv"), "w") do io
      println(io, join(_SE_TAP_CSV_COLUMNS, se_csv_delim))
      for t in something(res.tapEstimates, NamedTuple[])
        println(io, _se_tap_csv_line(_se_tap_row_record(t), se_csv_delim, se_csv_format))
      end
    end
  end
  # the verdict per row of the measurements.csv artifact, in its own file
  # next to it (the table stays loadable as it is): residual, normalized
  # residual, eliminated, constraint and the reason tokens per table row
  crit_idx = Set{Int}(crit_num)
  union!(crit_idx, crit_str)
  verdict_flagged = _se_write_measurement_verdict(joinpath(output_dir, "measurements_verdict.csv"), net, artifact, diag, eliminations, res, crit_idx, config.state_estimation.zero_injection === :constraint, se_csv_format; range_rows = range_rows)
  se_phase("writing_artifacts")
  writeSEStateCSV(net; file = joinpath(output_dir, "se_state.csv"), format = String(config.output.csv_format))

  # band verdict from the POST-elimination report: it must describe the
  # same state as the headline J (the pre-elimination report still tells
  # its story in se_diagnostics.md)
  verdict = diag.final_diagnostics.objective
  open(logfile, "a") do io
    println(io, "State estimation on ", basename(case_path))
    println(io, "measurements: ", summary.total, " rows from ", from_case_file ? "the case file" : basename(String(measurement_file)), zi_rows > 0 ? string(", plus ", zi_rows, " zero-injection row(s)") : "")
    println(io, "observability: ", obs.quality, :structural_islands in obs.notes ? " (structural islands)" : "")
    _se_log_criticality(io, net, obs)
    println(io, "converged in ", res.iterations, " iteration(s); J = ", round(res.objectiveJ; digits = 6), ", dof = ", res.dof, ", band reason = ", verdict.reason, res.activeObjective === nothing ? "" : string("; J_active = ", round(res.activeObjective.j; digits = 6), " (dof ", res.activeObjective.dof, ", without ", res.activeObjective.suppressed, " suppressed row(s); band verdict stays on J)"))
    println(io, _se_conditioning_line(res))
    if diag.topology_findings !== nothing
      # the stage-2 classification REPLACES the eliminations-exhausted
      # interpretation: this is a suspected topology error, not bad data
      for tf in diag.topology_findings
        println(io, "topology error suspected at station ", tf.location, " (eliminations exhausted, band high; suspects: ", tf.evidence, isempty(tf.notes) ? "" : string("; notes: ", join(tf.notes, ", ")), ")")
      end
      println(io, "suspicious measurements: ", length(suspicious), "; eliminations: ", length(eliminations), " (stop: ", diag.stop_reason, ", classified as topology, use the hypothesis test)")
    else
      println(io, "suspicious measurements: ", length(suspicious), "; eliminations: ", length(eliminations), " (stop: ", diag.stop_reason, ")")
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
        if tf.fallback_used
          println(io, "tap estimation: ", _se_tap_transformer_count(res.tapEstimates), " transformer(s), ", length(_se_tap_regulator_rows(res.tapEstimates)), " regulator(s), not estimated (", tap_run1[1] === :worse_fit ? "fallback_worse_fit" : "fallback_not_converged", ", see the two lines above)")
        else
          println(io, "tap estimation: ", _se_tap_transformer_count(res.tapEstimates), " transformer(s), ", length(_se_tap_regulator_rows(res.tapEstimates)), " regulator(s); J before fixation = ", round(tf.j_before; sigdigits = 4), " (dof ", tf.dof_before, "), J after fixation = ", round(tf.j_after; sigdigits = 4), " (dof ", tf.dof_after, ")", tf.fixed ? "" : "; NOT fixed (the released run did not converge)")
        end
        tf.offgrid_residual && println(io, "note :offgrid_tap_residual: the band failure appears only through the fixation; this is an off-grid true position or a wrong step table, NOT bad data")
        isempty(get(tf, :delta_j_notice, "")) || println(io, "notice: ", tf.delta_j_notice)
      else
        println(io, tap_enabled ? "tap estimation: every released tap was frozen by the guards; no tap state entered the solve" : "tap estimation: switched off, no tap state entered the solve")
      end
      for t in res.tapEstimates
        # one line per regulator: estimate, branch ratio and angle at it,
        # model step -> fixed step with the deviation
        println(io, "  ", t.name, isempty(t.changer) ? "" : "/" * t.changer, isempty(t.mrid) ? "" : " (mRID " * t.mrid * ")", " (", t.kind, "): s ", round(t.s_est; digits = 3), t.du_est === nothing ? "" : string(" (dU ", round(t.du_est; sigdigits = 4), ")"), ", ratio ", round(t.ratio_est; digits = 5), ", angle ", round(t.angle_est_deg; digits = 3), " deg, step ", _se_tap_step_text(t.step_model), " -> ", _se_tap_step_text(t.step_fixed), " (deviation ", _se_tap_deviation_text(t.deviation_steps), ")", t.sigma_s === nothing ? "" : string(", sigma_s ", round(t.sigma_s; sigdigits = 3), t.round_distance_sigma === nothing ? "" : string(", ", round(t.round_distance_sigma; sigdigits = 3), " sigma from the nearest step")), t.step_from_target === nothing ? "" : string(", step from the control target ", _se_tap_step_text(t.step_from_target)), _se_tap_out_of_range_text(t), t.status in (:estimated, :out_of_range) ? "" : t.status === :calculated ? " (CALCULATED, not estimated: $(t.reason))" : " (NOT ESTIMATED: $(t.reason))")
      end
      # rank_budget holds for every regulator of the run with one text: one
      # line with the count instead of one per regulator (8617 on case6495rte)
      nbudget = count(t -> t.reason === :rank_budget, res.tapEstimates)
      nbudget > 0 && println(io, "notice: ", nbudget, " tap regulator(s) ", _tap_notice_text(:rank_budget))
      for t in res.tapEstimates
        (isempty(t.notice) || t.reason === :rank_budget) && continue
        println(io, "notice: tap regulator ", t.name, isempty(t.changer) ? "" : "/" * t.changer, ": ", t.notice)
      end
      # the warning runse! logs for an estimate outside the band, same text
      for t in res.tapEstimates
        t.out_of_range && println(io, "warning: ", _tap_out_of_range_line(t))
      end
    end
    (!isempty(suspicious) || !isempty(eliminations) || !isempty(suppressedIds)) && println(io, "bad data: measurements_verdict.csv (", length(suspicious), " suspicious, ", length(eliminations), " eliminated, ", length(suppressedIds), " suppressed, ", length(downweightedIds), " down-weighted; ", verdict_flagged, " table row(s) flagged)")
    deltas_written && println(io, "deltas: se_deltas.csv (measured vs generator truth vs estimated, per row and per released tap)")
    link_rows !== nothing && println(io, "link flows: se_link_flows.csv (", length(link_rows), " closed link(s), ", n_link_measured, " with link flow measurements or an estimated circulation)")
    println(io, "Artifacts: measurements.csv, measurements_verdict.csv, se_diagnostics.md, se_view.md, se_state.csv", res.shuntEstimates !== nothing ? ", shunt_estimates.csv" : "", res.zeroInjectionMultipliers !== nothing ? ", se_zero_injection_multipliers.csv" : "", res.tapEstimates !== nothing ? ", se_tap_estimates.csv" : "", deltas_written ? ", se_deltas.csv" : "", link_rows !== nothing ? ", se_link_flows.csv" : "")
  end

  metadata = merge(
    base_metadata,
    Dict{String,Any}(
      "input_format_detected" => String(format),
      "se_measurement_rows" => summary.total,
      "se_zero_injection_rows" => zi_rows,
      # the numerics of the reported run: zero-injection mode, the matrix
      # of its last step and its condition number (nothing when singular)
      "se_zero_injection" => String(run_cfg.state_estimation.zero_injection),
      "se_condition_matrix" => res.conditioning === nothing ? nothing : String(res.conditioning.matrix),
      "se_condition_number" => (res.conditioning === nothing || !isfinite(res.conditioning.condition_number)) ? nothing : res.conditioning.condition_number,
      "se_condition_method" => res.conditioning === nothing ? nothing : String(res.conditioning.method),
      "se_zero_injection_multipliers" => res.zeroInjectionMultipliers === nothing ? 0 : length(res.zeroInjectionMultipliers),
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
      "se_suspicious" => length(suspicious),
      "se_eliminations" => length(eliminations),
      "se_robust" => robust,
      "se_deltas" => deltas_written,
      "se_link_flows" => link_rows === nothing ? 0 : length(link_rows),
      "se_link_flows_measured" => n_link_measured,
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
      "se_tap_estimation" => tap_enabled,
      # the impedance rule the run used (configuration, case file or override)
      "se_tap_impedance_correction" => String(run_cfg.state_estimation.tap_estimation.impedance_correction),
      "se_tap_count" => _se_tap_transformer_count(res.tapEstimates),
      "se_tap_regulators" => length(_se_tap_regulator_rows(res.tapEstimates)),
      # transformers kept out by the configuration: named in
      # tap_estimation.exclude, and machine transformers (calculated)
      "se_tap_excluded" => tap_n_excluded,
      "se_tap_machine" => tap_n_machine,
      "se_tap_not_estimated" => count(t -> t.status !== :estimated, _se_tap_regulator_rows(res.tapEstimates)),
      # of those, the estimates too far outside the band (transferred step kept)
      "se_tap_out_of_range" => count(t -> t.status === :out_of_range, _se_tap_regulator_rows(res.tapEstimates)),
      "se_tap_fixed" => res.tapFixation === nothing ? false : res.tapFixation.fixed,
      "se_tap_j_before" => res.tapFixation === nothing ? nothing : res.tapFixation.j_before,
      "se_tap_dof_before" => res.tapFixation === nothing ? nothing : res.tapFixation.dof_before,
      "se_tap_j_after" => res.tapFixation === nothing ? nothing : res.tapFixation.j_after,
      "se_tap_dof_after" => res.tapFixation === nothing ? nothing : res.tapFixation.dof_after,
      "se_tap_offgrid_residual" => res.tapFixation === nothing ? false : res.tapFixation.offgrid_residual,
      "se_tap_j_rise_exceeds_quantile" => res.tapFixation === nothing ? false : res.tapFixation.j_rise_exceeds_quantile,
      "se_tap_sequential_fixation" => res.tapFixation === nothing ? false : res.tapFixation.sequential,
      "se_topology_station_findings" => diag.topology_findings === nothing ? nothing : [Dict{String,Any}("location" => tf.location, "evidence" => tf.evidence, "notes" => [String(n) for n in tf.notes]) for tf in diag.topology_findings],
      # rows for the result page table, one per regulator: the columns of
      # se_tap_estimates.csv plus role
      "se_tap_estimates" => res.tapEstimates === nothing ? nothing : [_se_tap_row_record(t) for t in res.tapEstimates],
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
    length(suspicious), " suspicious, ",
    length(eliminations), " eliminated.",
    diag.topology_findings === nothing ? "" : string(" Topology error suspected at ", join((tf.location for tf in diag.topology_findings), ", "), "; run the hypothesis test."),
    tap_fallback_used ? string(" Tap estimation fallback: run 1 ", tap_run1[1] === :worse_fit ? "rejected" : "not converged", " (", tap_run1[1], ", ", tap_run1[2], " iterations), the result is the run without tap estimation.") :
    res.tapFixation === nothing ? "" : string(" Tap fixation (", _se_tap_transformer_count(res.tapEstimates), " transformer(s)): J ", round(res.tapFixation.j_before; sigdigits = 3), " -> ", round(res.tapFixation.j_after; sigdigits = 3), res.tapFixation.offgrid_residual ? "; off-grid tap residual, not bad data" : "", "."),
    (get(base_metadata, "se_set_tap_deviation", false) == true && !tap_enabled) ? " Hint: the set documents a tap deviation, but tap estimation is switched off; enable 'estimate taps'." : "",
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

  # the same precedence the power-flow path uses (resolve_config):
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

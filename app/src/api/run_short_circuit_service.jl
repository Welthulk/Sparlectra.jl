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
# file: src/api/run_short_circuit_service.jl
# purpose: Web UI/service short-circuit run (issue #277): CGMES
#          import + runShortCircuit! max and min, CSV artifacts, run.log
#          narrative with the coverage report, result.json through the normal
#          API result conventions. No power-flow solve is involved.

# One CSV per case; schema mirrors the ShortCircuitResult rows (the
# short-circuit artifact contract). `reasons` is a "; "-joined list inside
# one CSV field. `format` (technical/excel_de/excel_us, see
# `_resolve_detailed_csv_format`) sets delimiter and number formatting like
# every other CSV artifact of the run (issue #376 follow-up: with excel_de
# the power-flow files used ; and a decimal comma while this one still
# used , and a dot - one run, two formats).
function _write_short_circuit_csv(path::AbstractString, result::ShortCircuitResult; format = result_csv_format())::String
  header = ("bus", "vn_kV", "island", "status", "c", "zk_ohm", "rx_ratio", "ik_kA", "sk_MVA", "kappa", "ip_kA", "flagged", "reasons")
  rows = ((String(row.bus), Float64(row.vn_kV), row.island, row.status, Float64(row.c), Float64(row.zk_ohm), Float64(row.rx_ratio), Float64(row.ik_kA), Float64(row.sk_MVA), Float64(row.kappa), Float64(row.ip_kA), row.contains_defaulted_data, join(row.reasons, "; ")) for row in result.rows)
  write_result_csv(path, header, rows; format = format)
  return path
end

_sc_worst_row(result::ShortCircuitResult) = begin
  ok = [row for row in result.rows if row.status === :ok && isfinite(row.ik_kA)]
  isempty(ok) ? nothing : ok[argmax([row.ik_kA for row in ok])]
end

"""
    _run_short_circuit_scf(case_path, config, ...) -> SparlectraApiResult

Short-circuit run on a Sparlectra Case Format case (#342). The file carries
its own IEC 60909 source data in `sparlectra.components.sc_source` and may
define the study in `sparlectra.short_circuit`: which case is the headline
(`max`/`min`), the `c_factor`, and whether the sweep covers all buses or a
listed selection.

Both cases are always evaluated, so `short_circuit_max.csv` and
`short_circuit_min.csv` exist for every format; the file's `case` decides
which one the headline numbers come from. A `c_factor` in the file loses
against an explicitly configured `short_circuit.c_factor` (a case file sits
below API and CLI overrides), and applies when the configuration is at its
default.
"""
function _run_short_circuit_scf(case_path::AbstractString, config::SparlectraConfig, config_file::AbstractString, output_dir::AbstractString, run_id::String, logfile::String, result_file::String, base_metadata::Dict{String,Any})::SparlectraApiResult
  study, net, config = try
    imported = import_case(case_path, config; run_kind = :short_circuit)
    # step 3a: the run continues on the effective config of the import
    (imported.studies.short_circuit, imported.net, imported.config)
  catch err
    return _api_failure("invalid_case_file", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  sc = net.sc_sources
  n_sources = length(sc.external_network_injections) + length(sc.synchronous_machines) + length(sc.asynchronous_machines) + length(sc.equivalent_injections)
  if n_sources == 0
    return _api_failure("short_circuit_data_missing", "The case file carries no `sparlectra.components.sc_source` entry, so every bus would report :no_source. Export the case from a delivery that has source data, or add the sources to the file.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  sc_state = _sc_source_data_state(sc)
  sc_state.refused && return _api_failure("short_circuit_data_missing", _SC_NO_DATA_MESSAGE, run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)

  headline_case = Symbol(get(study, "case", "max"))
  # An explicit sweep names NODE ids; they resolve through the same reference
  # names the network uses, so an unknown id fails here and not with a row of
  # zeros in the table.
  buses = :all
  # the study block names its buses; a file written in PGM vocabulary states
  # the same thing with `fault` rows, and both are honoured
  selection = String(get(study, "sweep", "all_buses")) == "explicit" ? [Int(id) for id in study["buses"]] : scf_fault_nodes(case_path)
  bus_source = String(get(study, "sweep", "all_buses")) == "explicit" ? "short_circuit.buses" : "data.fault"
  if !isempty(selection)
    names = scf_extra_names(case_path)
    buses = try
      [get(() -> throw(ArgumentError("SCF short_circuit: node $(id) has no name in `extra`; add the name or use a bus name.")), names, id) for id in selection]
    catch err
      return _api_failure("invalid_case_file", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
    end
  end

  c_factor = config.shortcircuit.c_factor
  c_from_file = false
  if !("short_circuit.c_factor" in config.user_set_keys) && haskey(study, "c_factor")
    c_factor = Float64(study["c_factor"])
    c_from_file = true
  end

  # the engine reports every substituted default as a warning; they belong
  # to the result and go to run.log, not to the console
  sc_logger = _ScWarningCollector(String[], Logging.current_logger())
  sc_max, sc_min = try
    Logging.with_logger(sc_logger) do
      (
        runShortCircuit!(net; buses = buses, case = :max, c_factor = c_factor, sweep_method = config.shortcircuit.sweep_method, takahashi_min_buses = config.shortcircuit.takahashi_min_buses),
        runShortCircuit!(net; buses = buses, case = :min, c_factor = c_factor, sweep_method = config.shortcircuit.sweep_method, takahashi_min_buses = config.shortcircuit.takahashi_min_buses),
      )
    end
  catch err
    return _api_failure("short_circuit_error", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  headline = headline_case === :min ? sc_min : sc_max

  if !any(row.status === :ok for row in headline.rows)
    open(logfile, "a") do io
      println(io, "Short-circuit run: the case file's sources produced no usable bus (every row is :no_source or :isolated).")
      printShortCircuitResult(io, headline)
    end
    return _api_failure("short_circuit_data_missing", "The case file carries $(n_sources) source(s), but no bus could be evaluated - see the table in run.log.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  max_csv = _write_short_circuit_csv(joinpath(output_dir, "short_circuit_max.csv"), sc_max; format = String(config.output.csv_format))
  min_csv = _write_short_circuit_csv(joinpath(output_dir, "short_circuit_min.csv"), sc_min; format = String(config.output.csv_format))
  flagged = count(row.contains_defaulted_data for row in headline.rows)
  worst = _sc_worst_row(headline)

  open(logfile, "a") do io
    sc_state.defaults_only && println(io, "WARNING: none of the ", sc_state.total, " source(s) carries short-circuit data; every current below rests on default reactances, not on data of the case.")
    sc_state.partial && println(io, "WARNING: ", sc_state.total - sc_state.with_data, " of ", sc_state.total, " source(s) carry no short-circuit data and enter with default reactances; the rows they feed are flagged.")
    println(io, "Short-circuit run (IEC 60909-0) on ", basename(case_path), " (Sparlectra Case Format)")
    println(io, "study: ", isempty(study) ? "no short_circuit block in the file, defaults used" : "from the file's short_circuit block")
    println(io, "headline case: ", headline_case, buses === :all ? ", all buses" : ", $(length(buses)) selected bus(es) from $(bus_source)")
    println(io, "c-factor: ", c_factor > 0.0 ? string(c_factor, c_from_file ? " (from the case file)" : " (short_circuit.c_factor override)") : "IEC 60909-0 Table 1 by voltage level")
    println(io, "sources: ", n_sources, " (external network injections: ", length(sc.external_network_injections), ", synchronous: ", length(sc.synchronous_machines), ", asynchronous: ", length(sc.asynchronous_machines), ", equivalent injections: ", length(sc.equivalent_injections), ")")
    println(io)
    printShortCircuitResult(io, sc_max)
    println(io)
    printShortCircuitResult(io, sc_min)
    _sc_write_substituted_defaults(io, sc_logger)
    println(io, "Artifacts: ", basename(max_csv), ", ", basename(min_csv))
  end

  metadata = merge(
    base_metadata,
    Dict{String,Any}(
      "input_format_detected" => "scf",
      "sc_study_from_case_file" => !isempty(study),
      "sc_case_selected" => String(Symbol(headline_case)),
      "sc_sources" => n_sources,
      "sc_sources_with_data" => sc_state.with_data,
      "sc_defaults_only" => sc_state.defaults_only,
      "sc_case_rows" => length(headline.rows),
      "sc_flagged_rows" => flagged,
      "sc_c_factor" => c_factor,
      "sc_worst_bus" => worst === nothing ? "" : String(worst.bus),
      "sc_max_ik_kA" => worst === nothing ? NaN : worst.ik_kA,
      "sc_max_ip_kA" => worst === nothing ? NaN : worst.ip_kA,
      "artifact_status" => "completed",
      "solver_status" => "completed",
      "service_status" => "completed",
      "run_status" => "completed",
    ),
  )
  outcome = _sc_outcome(sc_state, flagged, length(headline.rows), "the CSV")
  return _finalize_api_result(_api_result(
    run_id = run_id,
    status = outcome.status,
    success = true,
    solution_available = false,
    reason = outcome.reason,
    message = outcome.message,
    casefile = String(case_path),
    config_file = String(config_file),
    output_dir = String(output_dir),
    logfile = logfile,
    result_file = result_file,
    metadata = metadata,
  ))
end

# The state of the source data, the same rule for every case format. A
# source carries data when the quantity its impedance comes from is in the
# case: x''d of a synchronous machine, the short-circuit current of a
# feeder, the reactance of an equivalent injection, the locked-rotor ratio
# of a motor. A source carries NOTHING when the base of a default is
# missing as well (a machine without rated power; a feeder, an equivalent
# or a motor without its quantity cannot be defaulted at all).
#   every source without anything -> the run refuses (data missing)
#   no source with data           -> warning, the table rests on defaults
#   some sources without data     -> warning, the rows they feed are flagged
#   every source with data        -> success
_sc_positive(x) = x !== nothing && isfinite(Float64(x)) && Float64(x) > 0.0
_sc_field(record, name::Symbol) = hasproperty(record, name) ? getproperty(record, name) : nothing
function _sc_source_data_state(sc)
  total = 0
  with_data = 0
  without_any = 0
  for m in sc.synchronous_machines
    total += 1
    has = _sc_positive(_sc_field(m, :satDirectSubtransX_pu))
    has && (with_data += 1)
    (!has && !_sc_positive(_sc_field(m, :ratedS_MVA))) && (without_any += 1)
  end
  for f in sc.external_network_injections
    total += 1
    has = _sc_positive(_sc_field(f, :maxInitialSymShCCurrent_A)) || _sc_positive(_sc_field(f, :minInitialSymShCCurrent_A))
    has ? (with_data += 1) : (without_any += 1)
  end
  for e in sc.equivalent_injections
    total += 1
    _sc_positive(_sc_field(e, :x_ohm)) ? (with_data += 1) : (without_any += 1)
  end
  for m in sc.asynchronous_machines
    total += 1
    # a motor without its locked-rotor ratio is skipped by the engine,
    # there is no default for it
    _sc_positive(_sc_field(m, :iaIrRatio)) ? (with_data += 1) : (without_any += 1)
  end
  refused = total > 0 && without_any == total
  return (total = total, with_data = with_data, without_any = without_any, refused = refused, defaults_only = total > 0 && with_data == 0 && !refused, partial = with_data > 0 && with_data < total, warning = !refused && with_data < total)
end

# status, reason and message of a finished run, by the state of its data
function _sc_outcome(state, flagged::Int, rows::Int, table::AbstractString)
  state.defaults_only && return (status = :warning, reason = "short_circuit_defaults_only", message = _sc_defaults_message(state, rows))
  state.partial && return (status = :warning, reason = "short_circuit_partial_defaults", message = "Short circuit completed with defaults: $(state.with_data) of $(state.total) source(s) carry short-circuit data, the others enter with default reactances. $(flagged) of $(rows) rows are flagged (reasons in $(table)).")
  flagged > 0 && return (status = :succeeded, reason = "short_circuit_flagged_lower_bound", message = "Short circuit completed - $(flagged) of $(rows) rows carry defaulted/skipped data; Ik'' is a lower bound on those buses (see reasons in $(table)).")
  return (status = :succeeded, reason = nothing, message = "Short circuit completed.")
end

const _SC_NO_DATA_MESSAGE = "None of the short-circuit sources of the case carries data an impedance can be derived from (a machine needs its subtransient reactance or at least its rated power, a feeder its short-circuit current)."

_sc_defaults_message(state, rows::Int) = "Short circuit completed on defaults: none of the $(state.total) source(s) carries short-circuit data, every current rests on default reactances$(state.without_any > 0 ? "; $(state.without_any) source(s) carry no rated power either" : ""). All $(rows) rows are flagged."

# Collects the warnings of a short-circuit evaluation (the engine reports
# every substituted default as one) so the service can write them into
# run.log; every other record goes on to the logger that was active.
struct _ScWarningCollector <: Logging.AbstractLogger
  messages::Vector{String}
  parent::Logging.AbstractLogger
end
Logging.min_enabled_level(logger::_ScWarningCollector) = min(Logging.Warn, Logging.min_enabled_level(logger.parent))
Logging.shouldlog(logger::_ScWarningCollector, level, _module, group, id) = level == Logging.Warn || Logging.shouldlog(logger.parent, level, _module, group, id)
Logging.catch_exceptions(logger::_ScWarningCollector) = Logging.catch_exceptions(logger.parent)
function Logging.handle_message(logger::_ScWarningCollector, level, message, _module, group, id, file, line; kwargs...)
  if level == Logging.Warn
    push!(logger.messages, string(message))
    return nothing
  end
  return Logging.handle_message(logger.parent, level, message, _module, group, id, file, line; kwargs...)
end

# the collected warnings as a block of run.log, each once (the maximum and
# the minimum case report the same substitution)
function _sc_write_substituted_defaults(io::IO, logger::_ScWarningCollector)
  isempty(logger.messages) && return nothing
  println(io)
  println(io, "Substituted defaults:")
  for line in unique(logger.messages)
    println(io, "  ", line)
  end
  return nothing
end

"""
    _run_short_circuit_powsybl(case_path, config, ...) -> SparlectraApiResult

Short-circuit run on a PowSyBl case (an IIDM file or a table bundle). The
sources are the connected generators of the file. A generator with the
IIDM extension `generatorShortCircuit` contributes its subtransient
reactance; one without it is evaluated with the engine's documented
default and every row that depends on it is flagged, so the result of a
file without any short-circuit data is complete but says on every row that
it rests on defaults. Both cases are evaluated on all buses
(`short_circuit_max.csv`, `short_circuit_min.csv`).

A run in which no generator carries short-circuit data ends with status
`warning` (reason `short_circuit_defaults_only`), and the statement leads
the message and `run.log`.

Failure behavior: `import_error` when the case does not import,
`short_circuit_data_missing` when the file carries no connected generator
or when no generator carries short-circuit data or a rated power,
`short_circuit_error` when the evaluation fails.
"""
function _run_short_circuit_powsybl(case_path::AbstractString, config::SparlectraConfig, config_file::AbstractString, output_dir::AbstractString, run_id::String, logfile::String, result_file::String, base_metadata::Dict{String,Any})::SparlectraApiResult
  failure = (reason, message) -> _api_failure(reason, message; run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  imported = try
    import_case(case_path, config; requested_format = :powsybl, run_kind = :short_circuit)
  catch err
    err isa PowerFlowAborted && rethrow()
    return failure("import_error", sprint(showerror, err))
  end
  net = imported.net
  config = imported.config
  machines = net.sc_sources.synchronous_machines
  isempty(machines) && return failure("short_circuit_data_missing", "The PowSyBl case carries no connected generator, so no bus has a short-circuit source.")
  sc_state = _sc_source_data_state(net.sc_sources)
  with_data = sc_state.with_data
  without_any = sc_state.without_any
  sc_state.refused && return failure("short_circuit_data_missing", _SC_NO_DATA_MESSAGE * " In an IIDM file that is the extension generatorShortCircuit and the attribute ratedS of the generators.")

  c_factor = config.shortcircuit.c_factor
  # the engine reports every substituted default as a warning; here they
  # are part of the result (flags in the tables, reasons in the CSV, the
  # count in the message), so they go to run.log and not to the console
  sc_logger = _ScWarningCollector(String[], Logging.current_logger())
  sc_max, sc_min = try
    Logging.with_logger(sc_logger) do
      (
        runShortCircuit!(net; case = :max, c_factor = c_factor, sweep_method = config.shortcircuit.sweep_method, takahashi_min_buses = config.shortcircuit.takahashi_min_buses),
        runShortCircuit!(net; case = :min, c_factor = c_factor, sweep_method = config.shortcircuit.sweep_method, takahashi_min_buses = config.shortcircuit.takahashi_min_buses),
      )
    end
  catch err
    err isa PowerFlowAborted && rethrow()
    return failure("short_circuit_error", sprint(showerror, err))
  end
  if !any(row.status === :ok for row in sc_max.rows)
    open(logfile, "a") do io
      println(io, "Short-circuit run: the generators of the case produced no usable bus (every row is :no_source or :isolated).")
      printShortCircuitResult(io, sc_max)
    end
    return failure("short_circuit_data_missing", "The case carries $(length(machines)) generator(s), but no bus could be evaluated - see the table in run.log.")
  end

  max_csv = _write_short_circuit_csv(joinpath(output_dir, "short_circuit_max.csv"), sc_max; format = String(config.output.csv_format))
  min_csv = _write_short_circuit_csv(joinpath(output_dir, "short_circuit_min.csv"), sc_min; format = String(config.output.csv_format))
  flagged = count(row.contains_defaulted_data for row in sc_max.rows)
  worst = _sc_worst_row(sc_max)
  data_note = with_data == length(machines) ? "every generator carries short-circuit data" : with_data == 0 ? "no generator carries short-circuit data (IIDM extension generatorShortCircuit): every source is evaluated with the default x''d" : "$(with_data) of $(length(machines)) generators carry short-circuit data (IIDM extension generatorShortCircuit), the others are evaluated with the default x''d"
  without_any > 0 && (data_note *= "; $(without_any) of them carry no rated power either")
  # without any short-circuit data the table is complete but rests on the
  # default reactance throughout: the run ends with a warning, and the
  # statement leads the message and the log
  defaults_only = with_data == 0
  open(logfile, "a") do io
    defaults_only && println(io, "WARNING: none of the ", length(machines), " generator(s) carries short-circuit data; every current below rests on the default x''d, not on data of the file.")
    sc_state.partial && println(io, "WARNING: ", sc_state.total - sc_state.with_data, " of ", sc_state.total, " generator(s) carry no short-circuit data and enter with the default x''d; the rows they feed are flagged.")
    println(io, "Short-circuit run (IEC 60909-0) on ", basename(case_path), " (PowSyBl)")
    println(io, "c-factor: ", c_factor > 0.0 ? string(c_factor, " (short_circuit.c_factor override)") : "IEC 60909-0 Table 1 by voltage level")
    println(io, "sources: ", length(machines), " generator(s); ", data_note)
    println(io)
    printShortCircuitResult(io, sc_max)
    println(io)
    printShortCircuitResult(io, sc_min)
    _sc_write_substituted_defaults(io, sc_logger)
    println(io, "Artifacts: ", basename(max_csv), ", ", basename(min_csv))
  end

  metadata = merge(
    base_metadata,
    Dict{String,Any}(
      "input_format_detected" => "powsybl",
      "sc_sources" => length(machines),
      "sc_sources_with_data" => with_data,
      "sc_sources_without_data_and_rating" => without_any,
      "sc_defaults_only" => defaults_only,
      "sc_case_rows" => length(sc_max.rows),
      "sc_flagged_rows" => flagged,
      "sc_c_factor" => c_factor,
      "sc_worst_bus" => worst === nothing ? "" : String(worst.bus),
      "sc_max_ik_kA" => worst === nothing ? NaN : worst.ik_kA,
      "sc_max_ip_kA" => worst === nothing ? NaN : worst.ip_kA,
      "artifact_status" => "completed",
      "solver_status" => "completed",
      "service_status" => "completed",
      "run_status" => "completed",
    ),
  )
  outcome = _sc_outcome(sc_state, flagged, length(sc_max.rows), "short_circuit_max.csv")
  message = defaults_only ? outcome.message * " In an IIDM file the data is the extension generatorShortCircuit of the generators." : outcome.message
  return _finalize_api_result(_api_result(
    run_id = run_id,
    status = outcome.status,
    success = true,
    solution_available = false,
    reason = outcome.reason,
    message = message,
    casefile = String(case_path),
    config_file = String(config_file),
    output_dir = String(output_dir),
    logfile = logfile,
    result_file = result_file,
    metadata = metadata,
  ))
end

"""
    _run_short_circuit_service(case_path, config_file, output_dir, run_id) -> SparlectraApiResult

Service backend of the Web UI "Short circuit" button: imports the CGMES
delivery (boundary resolution shared with the power-flow path), evaluates
`runShortCircuit!` for both cases, writes `short_circuit_max.csv` /
`short_circuit_min.csv` plus a `run.log` narrative including the
short-circuit coverage report, and returns a `SparlectraApiResult` through
the normal result conventions (`result.json`, artifact scan, run registry).

Failure behavior — explicit reasons instead of empty tables:
- `short_circuit_requires_cgmes`: the case is neither a CGMES delivery nor a
  Sparlectra Case Format case (MATPOWER/DTF carry no harvested short-circuit
  data yet, issue #277 follow-up).
- `cgmes_import_error` / `cgmes_boundary_missing`: import failed.
- `short_circuit_data_missing`: the delivery imported but carries no usable
  short-circuit source (every bus would report `:no_source`).
"""
function _run_short_circuit_service(case_path::AbstractString, config_file::AbstractString, output_dir::AbstractString, run_id::String; config_overrides::AbstractDict = Dict{String,Any}())::SparlectraApiResult
  mkpath(output_dir)
  logfile = joinpath(output_dir, "run.log")
  result_file = joinpath(output_dir, "result.json")
  base_metadata = Dict{String,Any}("run_mode" => "short_circuit")

  # The same precedence the power-flow path uses (resolve_config), minus
  # request overrides (a study run has no config form of its own): case
  # configuration file, the case file's deprecated block, general file,
  # defaults.
  config = try
    resolve_config(config_file, case_path, config_overrides).config
  catch err
    return _api_failure(_config_resolve_reason(err), sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  detected = _detect_case_format(case_path)
  # A Sparlectra Case Format case carries its own IEC 60909 source data and
  # may define the study itself (#342), so it runs here like a delivery.
  if detected === :scf
    return _run_short_circuit_scf(case_path, config, config_file, output_dir, run_id, logfile, result_file, base_metadata)
  end
  if detected === :powsybl
    return _run_short_circuit_powsybl(case_path, config, config_file, output_dir, run_id, logfile, result_file, base_metadata)
  end
  if detected !== :cgmes
    return _api_failure("short_circuit_requires_cgmes", "Short-circuit evaluation needs a CGMES delivery, a PowSyBl case or a Sparlectra Case Format case carrying source data — MATPOWER/DTF cases carry no harvested short-circuit source data yet (issue #277 follow-up).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  imported_cgmes = try
    import_case(case_path, config; requested_format = :cgmes, run_kind = :short_circuit)
  catch err
    message = sprint(showerror, err)
    # Same convention as the power-flow path: the typed import error carries
    # the full analysis — persist it where the run's artifacts live.
    if err isa CGMESImporter.CGMESImportError
      try
        open(joinpath(output_dir, "cgmes.log"), "w") do io
          println(io, "# CGMES import report — IMPORT FAILED")
          println(io, "error:  ", message)
          println(io)
          println(io, "## Import analysis")
          print(io, err.analysis)
        end
      catch logerr
        @warn "could not write diagnostic cgmes.log" exception = logerr
      end
    end
    reason = occursin("boundary set missing", message) || occursin("unresolved topology references", message) || occursin("no resolvable BaseVoltage", message) ? "cgmes_boundary_missing" : "cgmes_import_error"
    return _api_failure(reason, message; run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  cgmes_result = imported_cgmes.provenance["cgmes_result"]
  boundary_autodetected = imported_cgmes.provenance["cgmes_boundary_autodetected"]::Bool
  # step 3a: the run continues on the effective config of the import, and
  # the start decision is on record like in the power-flow run log
  config = imported_cgmes.config
  sc_start_decision = get(imported_cgmes.provenance, "cgmes_start_decision", nothing)
  sc_start_decision === nothing || open(logfile, "a") do io
    println(io, sc_start_decision)
  end

  sc_state = _sc_source_data_state(cgmes_result.shortcircuit)
  if sc_state.refused
    open(logfile, "a") do io
      println(io, "Short-circuit run: no source of the delivery carries data an impedance can be derived from.")
      printShortCircuitCoverage(io, cgmes_result.shortcircuit)
    end
    return _api_failure("short_circuit_data_missing", _SC_NO_DATA_MESSAGE * " See the coverage report in run.log.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  c_factor = config.shortcircuit.c_factor
  # the substituted defaults go to run.log, as on the other two paths
  sc_logger = _ScWarningCollector(String[], Logging.current_logger())
  sc_max, sc_min = Logging.with_logger(sc_logger) do
    (
      runShortCircuit!(cgmes_result; case = :max, c_factor = c_factor, sweep_method = config.shortcircuit.sweep_method, takahashi_min_buses = config.shortcircuit.takahashi_min_buses),
      runShortCircuit!(cgmes_result; case = :min, c_factor = c_factor, sweep_method = config.shortcircuit.sweep_method, takahashi_min_buses = config.shortcircuit.takahashi_min_buses),
    )
  end

  # A delivery without any usable source produces only :no_source/:isolated
  # rows — that is missing data, not a zero-current network.
  if !any(row.status === :ok for row in sc_max.rows)
    open(logfile, "a") do io
      println(io, "Short-circuit run: no usable short-circuit source in the delivery.")
      printShortCircuitCoverage(io, cgmes_result.shortcircuit)
    end
    return _api_failure("short_circuit_data_missing", "The delivery imported, but no usable short-circuit source (machine x''_d or feeder Ik) was found — see the coverage report in run.log.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  max_csv = _write_short_circuit_csv(joinpath(output_dir, "short_circuit_max.csv"), sc_max; format = String(config.output.csv_format))
  min_csv = _write_short_circuit_csv(joinpath(output_dir, "short_circuit_min.csv"), sc_min; format = String(config.output.csv_format))

  flagged = count(row.contains_defaulted_data for row in sc_max.rows)
  worst = _sc_worst_row(sc_max)
  open(logfile, "a") do io
    sc_state.defaults_only && println(io, "WARNING: none of the ", sc_state.total, " source(s) carries short-circuit data; every current below rests on default reactances, not on data of the delivery.")
    sc_state.partial && println(io, "WARNING: ", sc_state.total - sc_state.with_data, " of ", sc_state.total, " source(s) carry no short-circuit data and enter with default reactances; the rows they feed are flagged.")
    println(io, "Short-circuit run (IEC 60909-0) on ", basename(case_path), boundary_autodetected ? " (boundary set autodetected)" : "")
    println(io, "c-factor: ", c_factor > 0.0 ? string(c_factor, " (short_circuit.c_factor override)") : "IEC 60909-0 Table 1 by voltage level")
    println(io)
    printShortCircuitResult(io, sc_max)
    println(io)
    printShortCircuitResult(io, sc_min)
    println(io)
    println(io, "Harvested short-circuit source data:")
    printShortCircuitCoverage(io, cgmes_result.shortcircuit)
    _sc_write_substituted_defaults(io, sc_logger)
    println(io, "Artifacts: ", basename(max_csv), ", ", basename(min_csv))
  end

  metadata = merge(
    base_metadata,
    Dict{String,Any}(
      "input_format_detected" => "cgmes",
      "cgmes_buses" => length(cgmes_result.net.nodeVec),
      "sc_sources" => sc_state.total,
      "sc_sources_with_data" => sc_state.with_data,
      "sc_defaults_only" => sc_state.defaults_only,
      "sc_case_rows" => length(sc_max.rows),
      "sc_flagged_rows" => flagged,
      "sc_c_factor" => c_factor,
      "sc_worst_bus" => worst === nothing ? "" : String(worst.bus),
      "sc_max_ik_kA" => worst === nothing ? NaN : worst.ik_kA,
      "sc_max_ip_kA" => worst === nothing ? NaN : worst.ip_kA,
      # explicit run-status keys so the history/result views render like any
      # other completed run
      "artifact_status" => "completed",
      "solver_status" => "completed",
      "service_status" => "completed",
      "run_status" => "completed",
    ),
  )
  # A flagged Ik''max is a lower bound (skipped/defaulted contributions) —
  # surface that as the run message so nobody reads the table unwarned.
  outcome = _sc_outcome(sc_state, flagged, length(sc_max.rows), "short_circuit_max.csv")
  result = _api_result(
    run_id = run_id,
    status = outcome.status,
    success = true,
    solution_available = false,
    reason = outcome.reason,
    message = outcome.message,
    casefile = String(case_path),
    config_file = String(config_file),
    output_dir = String(output_dir),
    logfile = logfile,
    result_file = result_file,
    metadata = metadata,
  )
  return _finalize_api_result(result)
end

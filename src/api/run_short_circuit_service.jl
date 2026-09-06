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
# short-circuit artifact contract). `reasons` is a "; "-joined list inside one CSV field.
function _write_short_circuit_csv(path::AbstractString, result::ShortCircuitResult)::String
  open(path, "w") do io
    println(io, "bus,vn_kV,island,status,c,zk_ohm,rx_ratio,ik_kA,sk_MVA,kappa,ip_kA,flagged,reasons")
    for row in result.rows
      println(io, join((_csv_field(String(row.bus), ','), row.vn_kV, row.island, row.status, row.c, row.zk_ohm, row.rx_ratio, row.ik_kA, row.sk_MVA, row.kappa, row.ip_kA, row.contains_defaulted_data, _csv_field(join(row.reasons, "; "), ',')), ','))
    end
  end
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

  sc_max, sc_min = try
    (
      runShortCircuit!(net; buses = buses, case = :max, c_factor = c_factor, sweep_method = config.shortcircuit.sweep_method, takahashi_min_buses = config.shortcircuit.takahashi_min_buses),
      runShortCircuit!(net; buses = buses, case = :min, c_factor = c_factor, sweep_method = config.shortcircuit.sweep_method, takahashi_min_buses = config.shortcircuit.takahashi_min_buses),
    )
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

  max_csv = _write_short_circuit_csv(joinpath(output_dir, "short_circuit_max.csv"), sc_max)
  min_csv = _write_short_circuit_csv(joinpath(output_dir, "short_circuit_min.csv"), sc_min)
  flagged = count(row.contains_defaulted_data for row in headline.rows)
  worst = _sc_worst_row(headline)

  open(logfile, "a") do io
    println(io, "Short-circuit run (IEC 60909-0) on ", basename(case_path), " (Sparlectra Case Format)")
    println(io, "study: ", isempty(study) ? "no short_circuit block in the file, defaults used" : "from the file's short_circuit block")
    println(io, "headline case: ", headline_case, buses === :all ? ", all buses" : ", $(length(buses)) selected bus(es) from $(bus_source)")
    println(io, "c-factor: ", c_factor > 0.0 ? string(c_factor, c_from_file ? " (from the case file)" : " (short_circuit.c_factor override)") : "IEC 60909-0 Table 1 by voltage level")
    println(io, "sources: ", n_sources, " (external network injections: ", length(sc.external_network_injections), ", synchronous: ", length(sc.synchronous_machines), ", asynchronous: ", length(sc.asynchronous_machines), ", equivalent injections: ", length(sc.equivalent_injections), ")")
    println(io)
    printShortCircuitResult(io, sc_max)
    println(io)
    printShortCircuitResult(io, sc_min)
    println(io, "Artifacts: ", basename(max_csv), ", ", basename(min_csv))
  end

  metadata = merge(
    base_metadata,
    Dict{String,Any}(
      "input_format_detected" => "scf",
      "sc_study_from_case_file" => !isempty(study),
      "sc_case_selected" => String(Symbol(headline_case)),
      "sc_sources" => n_sources,
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
  message = flagged > 0 ? "Short circuit completed - $(flagged) of $(length(headline.rows)) rows carry defaulted/skipped data; Ik'' is a lower bound on those buses (see reasons in the CSV)." : "Short circuit completed."
  return _finalize_api_result(_api_result(
    run_id = run_id,
    status = :succeeded,
    success = true,
    solution_available = false,
    reason = flagged > 0 ? "short_circuit_flagged_lower_bound" : nothing,
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
function _run_short_circuit_service(case_path::AbstractString, config_file::AbstractString, output_dir::AbstractString, run_id::String)::SparlectraApiResult
  mkpath(output_dir)
  logfile = joinpath(output_dir, "run.log")
  result_file = joinpath(output_dir, "result.json")
  base_metadata = Dict{String,Any}("run_mode" => "short_circuit")

  # The same precedence the power-flow path uses (resolve_config, D5), minus
  # request overrides (a study run has no config form of its own): case
  # configuration file, the case file's deprecated block, general file,
  # defaults.
  config = try
    resolve_config(config_file, case_path).config
  catch err
    return _api_failure(_config_resolve_reason(err), sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  detected = _detect_case_format(case_path)
  # A Sparlectra Case Format case carries its own IEC 60909 source data and
  # may define the study itself (#342), so it runs here like a delivery.
  if detected === :scf
    return _run_short_circuit_scf(case_path, config, config_file, output_dir, run_id, logfile, result_file, base_metadata)
  end
  if detected !== :cgmes
    return _api_failure("short_circuit_requires_cgmes", "Short-circuit evaluation needs a CGMES delivery or a Sparlectra Case Format case carrying source data — MATPOWER/DTF cases carry no harvested short-circuit source data yet (issue #277 follow-up).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
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

  c_factor = config.shortcircuit.c_factor
  sc_max = runShortCircuit!(cgmes_result; case = :max, c_factor = c_factor, sweep_method = config.shortcircuit.sweep_method, takahashi_min_buses = config.shortcircuit.takahashi_min_buses)
  sc_min = runShortCircuit!(cgmes_result; case = :min, c_factor = c_factor, sweep_method = config.shortcircuit.sweep_method, takahashi_min_buses = config.shortcircuit.takahashi_min_buses)

  # A delivery without any usable source produces only :no_source/:isolated
  # rows — that is missing data, not a zero-current network.
  if !any(row.status === :ok for row in sc_max.rows)
    open(logfile, "a") do io
      println(io, "Short-circuit run: no usable short-circuit source in the delivery.")
      printShortCircuitCoverage(io, cgmes_result.shortcircuit)
    end
    return _api_failure("short_circuit_data_missing", "The delivery imported, but no usable short-circuit source (machine x''_d or feeder Ik) was found — see the coverage report in run.log.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  max_csv = _write_short_circuit_csv(joinpath(output_dir, "short_circuit_max.csv"), sc_max)
  min_csv = _write_short_circuit_csv(joinpath(output_dir, "short_circuit_min.csv"), sc_min)

  flagged = count(row.contains_defaulted_data for row in sc_max.rows)
  worst = _sc_worst_row(sc_max)
  open(logfile, "a") do io
    println(io, "Short-circuit run (IEC 60909-0) on ", basename(case_path), boundary_autodetected ? " (boundary set autodetected)" : "")
    println(io, "c-factor: ", c_factor > 0.0 ? string(c_factor, " (short_circuit.c_factor override)") : "IEC 60909-0 Table 1 by voltage level")
    println(io)
    printShortCircuitResult(io, sc_max)
    println(io)
    printShortCircuitResult(io, sc_min)
    println(io)
    println(io, "Harvested short-circuit source data:")
    printShortCircuitCoverage(io, cgmes_result.shortcircuit)
    println(io, "Artifacts: ", basename(max_csv), ", ", basename(min_csv))
  end

  metadata = merge(
    base_metadata,
    Dict{String,Any}(
      "input_format_detected" => "cgmes",
      "cgmes_buses" => length(cgmes_result.net.nodeVec),
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
  message = flagged > 0 ? "Short circuit completed — $(flagged) of $(length(sc_max.rows)) rows carry defaulted/skipped data; Ik''max is a lower bound on those buses (see reasons in short_circuit_max.csv)." : "Short circuit completed."
  result = _api_result(
    run_id = run_id,
    status = :succeeded,
    success = true,
    solution_available = false,
    reason = flagged > 0 ? "short_circuit_flagged_lower_bound" : nothing,
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

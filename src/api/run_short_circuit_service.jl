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
    _run_short_circuit_service(case_path, config_file, output_dir, run_id) -> SparlectraApiResult

Service backend of the Web UI "Short circuit" button: imports the CGMES
delivery (boundary resolution shared with the power-flow path), evaluates
`runShortCircuit!` for both cases, writes `short_circuit_max.csv` /
`short_circuit_min.csv` plus a `run.log` narrative including the
short-circuit coverage report, and returns a `SparlectraApiResult` through
the normal result conventions (`result.json`, artifact scan, run registry).

Failure behavior — explicit reasons instead of empty tables:
- `short_circuit_requires_cgmes`: the case is not a CGMES delivery
  (MATPOWER/DTF carry no harvested short-circuit data yet, issue #277
  follow-up).
- `cgmes_import_error` / `cgmes_boundary_missing`: import failed.
- `short_circuit_data_missing`: the delivery imported but carries no usable
  short-circuit source (every bus would report `:no_source`).
"""
function _run_short_circuit_service(case_path::AbstractString, config_file::AbstractString, output_dir::AbstractString, run_id::String)::SparlectraApiResult
  mkpath(output_dir)
  logfile = joinpath(output_dir, "run.log")
  result_file = joinpath(output_dir, "result.json")
  base_metadata = Dict{String,Any}("run_mode" => "short_circuit")

  config = try
    load_sparlectra_config(config_file; reload = true)
  catch err
    return _api_failure("invalid_configuration", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  if _detect_case_format(case_path) !== :cgmes
    return _api_failure("short_circuit_requires_cgmes", "Short-circuit evaluation needs a CGMES delivery — MATPOWER/DTF cases carry no harvested short-circuit source data yet (issue #277 follow-up).", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  cgmes_cfg = config.cgmes
  paths, boundary_autodetected = _cgmes_delivery_paths(case_path, cgmes_cfg)
  cgmes_result = try
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
      name = basename(case_path),
    )
  catch err
    message = sprint(showerror, err)
    reason = occursin("boundary set missing", message) || occursin("unresolved topology references", message) ? "cgmes_boundary_missing" : "cgmes_import_error"
    return _api_failure(reason, message; run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  c_factor = config.shortcircuit.c_factor
  sc_max = runShortCircuit!(cgmes_result; case = :max, c_factor = c_factor)
  sc_min = runShortCircuit!(cgmes_result; case = :min, c_factor = c_factor)

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

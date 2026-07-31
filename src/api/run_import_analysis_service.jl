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
# file: src/api/run_import_analysis_service.jl
# purpose: Web UI/service "Analyze import" run: parse the CGMES delivery
#          without mapping or solving, write the importFailureAnalysis report
#          as import_analysis.txt, and answer through the normal result
#          conventions. A fast pre-check before a (much slower) full import.

"""
    _run_import_analysis_service(case_path, config_file, output_dir, run_id) -> SparlectraApiResult

Service backend of the Web UI "Analyze import" button: loads the CGMES
delivery files (boundary resolution shared with the power-flow path — same
`_cgmes_delivery_paths`), builds the [`importFailureAnalysis`](@ref
Sparlectra.CGMESImporter.importFailureAnalysis) report, and writes it as
`import_analysis.txt` next to `run.log`. No mapping and no power-flow solve
is involved, so the check answers well before a full import would.

Failure reasons: `import_analysis_requires_cgmes` for non-CGMES cases,
`cgmes_read_error` when the delivery files cannot even be parsed.
"""
function _run_import_analysis_service(case_path::AbstractString, config_file::AbstractString, output_dir::AbstractString, run_id::String)::SparlectraApiResult
  mkpath(output_dir)
  logfile = joinpath(output_dir, "run.log")
  result_file = joinpath(output_dir, "result.json")
  base_metadata = Dict{String,Any}("run_mode" => "import_analysis")

  config = try
    load_sparlectra_config(config_file; reload = true)
  catch err
    return _api_failure("invalid_configuration", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  if _detect_case_format(case_path) !== :cgmes
    return _api_failure("import_analysis_requires_cgmes", "Import analysis reads CGMES model headers — MATPOWER/DTF cases carry none.", run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end

  paths, boundary_autodetected = _cgmes_delivery_paths(case_path, config.cgmes)
  store = try
    CGMESImporter.loadCGMES(length(paths) == 1 ? paths[1] : paths)
  catch err
    return _api_failure("cgmes_read_error", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  end
  report = CGMESImporter.importFailureAnalysis(store)

  analysis_file = joinpath(output_dir, "import_analysis.txt")
  open(analysis_file, "w") do io
    print(io, report)
  end

  # Compact stats for the result metadata/summary row; the full detail lives
  # in the report itself.
  stats_nt = CGMESImporter.importabilityStats(store)
  missing_deps = stats_nt.missing_dependencies
  unresolved_n = stats_nt.unresolved_count
  verdict_lines = [l for l in split(report, '\n') if startswith(l, "Verdict:")]
  verdict = isempty(verdict_lines) ? "" : String(replace(first(verdict_lines), "Verdict: " => ""))

  open(logfile, "a") do io
    println(io, "Import analysis on ", basename(case_path), boundary_autodetected ? " (boundary set autodetected)" : "")
    println(io)
    print(io, report)
    println(io)
    println(io, "Artifact: ", basename(analysis_file))
  end

  importable = stats_nt.importable

  stats = Dict{String,Any}(
    "input_format_detected" => "cgmes",
    "import_analysis_files" => length(store.files),
    "import_analysis_missing_dependencies" => length(missing_deps),
    "import_analysis_unresolved_refs" => unresolved_n,
    "import_analysis_verdict" => verdict,
  )

  if !importable
    # The user asked "will this import?" and the answer is no — that is a
    # FAILED run, with the full explanation in import_analysis.txt.
    message = "Import analysis: the delivery will NOT import — $(length(missing_deps)) missing declared dependency(ies), $(unresolved_n) unresolved reference(s). $(verdict) See import_analysis.txt."
    return _api_failure("import_analysis_not_importable", message; run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = merge(base_metadata, stats))
  end

  metadata = merge(
    base_metadata,
    stats,
    Dict{String,Any}(
      "artifact_status" => "completed",
      "solver_status" => "completed",
      "service_status" => "completed",
      "run_status" => "completed",
    ),
  )
  message = unresolved_n == 0 ? "Import analysis completed — the delivery declares no missing dependencies and resolves all references." :
    "Import analysis completed — importable, but $(unresolved_n) non-fatal unresolved reference(s). See import_analysis.txt."
  result = _api_result(
    run_id = run_id,
    status = :succeeded,
    success = true,
    solution_available = false,
    reason = unresolved_n == 0 ? nothing : "import_analysis_nonfatal_gaps",
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

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

# file: src/api/powerflow_service.jl
# purpose: local PowerFlow service layer: in-memory run registry, Web UI case
#          resolution, and the start_powerflow_run / get_powerflow_result
#          service entry points
const _POWERFLOW_SERVICE_RUNS = Dict{String,SparlectraApiResult}()
const _POWERFLOW_SERVICE_LOCK = ReentrantLock()


function _registered_powerflow_run(run_id::AbstractString)
  id = String(run_id)
  return lock(_POWERFLOW_SERVICE_LOCK) do
    get(_POWERFLOW_SERVICE_RUNS, id, nothing)
  end
end

function _register_powerflow_run!(result::SparlectraApiResult)
  lock(_POWERFLOW_SERVICE_LOCK) do
    _POWERFLOW_SERVICE_RUNS[result.run_id] = result
  end
  return result
end

include("service_json.jl")
include("artifact_registry.jl")
include("run_index.jl")
include("webui_jobs.jl")

const GENERATED_MATPOWER_JL_CACHE_MESSAGE = "Generated MATPOWER .jl cache files are not user-selectable. Please use or fetch the corresponding .m MATPOWER case source."

_matpower_cache_jl_bypass_reason() = "generated_jl_cache_hidden_from_webui"

_is_for002_reference_dat(path::AbstractString)::Bool = occursin(r"^for002.*\.dat$"i, basename(strip(String(path))))

_for002_primary_case_message() = "FOR002.DAT is a reference/result file and cannot be used as the primary DTF network input case. Use a runnable DTF network case as the case and enter FOR002.DAT as optional FOR002 reference file."

function _canonical_matpower_source_for_webui(path::AbstractString, case_directory::AbstractString)::String
  case_path = abspath(path)
  extension = lowercase(splitext(case_path)[2])
  extension == ".m" && return case_path
  if extension == ".jl" && dirname(case_path) == abspath(case_directory)
    # Keep Web UI MATPOWER runs on the canonical .m source. Generated .jl files
    # in the cache are internal artifacts and should not become the run input.
    m_path = first(splitext(case_path)) * ".m"
    isfile(m_path) && return abspath(m_path)
    throw(ArgumentError(GENERATED_MATPOWER_JL_CACHE_MESSAGE))
  end
  return case_path
end

function _resolve_powerflow_casefile(
  casefile::AbstractString,
  case_directory::AbstractString;
  ensure_casefile_fn = FetchMatpowerCase.ensure_casefile
)::String
  requested = strip(casefile)
  isempty(requested) && throw(ArgumentError("PowerFlow casefile must not be empty."))
  occursin(r"^[A-Za-z][A-Za-z0-9+.-]*://", requested) && throw(ArgumentError("MATPOWER case URLs are not accepted."))

  # A CGMES delivery is a .zip or an unpacked directory; both are handed to
  # the importer unchanged (its container layer resolves them).
  isdir(requested) && return abspath(requested)
  extension = lowercase(splitext(requested)[2])
  extension in (".m", ".jl", ".dat", ".zip") || throw(ArgumentError("Unsupported casefile extension: $(requested) (expected .m, .jl, .DAT, or .zip)"))
  if isfile(requested)
    if extension == ".zip"
      return abspath(requested)
    end
    if extension == ".dat"
      role = _webui_classify_dat_content(requested)
      _webui_is_runnable_dat_role(role) || throw(ArgumentError("$(basename(requested)) is a $(_webui_dat_role_label(role)) file and cannot be used as the primary PowerFlow case. Choose a runnable DTF network case."))
      return abspath(requested)
    end
    return _canonical_matpower_source_for_webui(requested, case_directory)
  end
  occursin(r"[\\/]", requested) && throw(ArgumentError("Case file not found: $(requested)"))

  trusted_directory = abspath(case_directory)
  mkpath(trusted_directory)
  if extension == ".jl"
    # A cache-local .jl request may only bypass to its matching .m source; a
    # standalone generated cache file is rejected so users see the source case.
    requested_jl = joinpath(trusted_directory, requested)
    requested_m = first(splitext(requested_jl)) * ".m"
    if isfile(requested_m)
      return abspath(requested_m)
    end
    isfile(requested_jl) && throw(ArgumentError(GENERATED_MATPOWER_JL_CACHE_MESSAGE))
  end
  if extension == ".zip"
    local_zip = joinpath(trusted_directory, requested)
    isfile(local_zip) && return abspath(local_zip)
    throw(ArgumentError("Case file not found: $(requested)"))
  end
  if extension == ".dat"
    local_dat = joinpath(trusted_directory, requested)
    if isfile(local_dat)
      role = _webui_classify_dat_content(local_dat)
      _webui_is_runnable_dat_role(role) || throw(ArgumentError("$(requested) is a $(_webui_dat_role_label(role)) file and cannot be used as the primary PowerFlow case. Choose a runnable DTF network case."))
      return abspath(local_dat)
    end
    throw(ArgumentError("Case file not found: $(requested)"))
  end
  local_m = extension == ".m" ? joinpath(trusted_directory, requested) : joinpath(trusted_directory, first(splitext(requested)) * ".m")
  isfile(local_m) && return abspath(local_m)
  resolved = try
    ensure_casefile_fn(extension == ".m" ? requested : first(splitext(requested)) * ".m"; outdir = trusted_directory, to_jl = false)
  catch err
    fallback_name = extension == ".m" ? requested : first(splitext(requested)) * ".m"
    fallback_path = joinpath(trusted_directory, fallback_name)
    if isfile(fallback_path)
      @warn "MATPOWER case resolution could not generate the Julia case; using the downloaded .m case" casefile = fallback_path exception = (err, catch_backtrace())
      return abspath(fallback_path)
    end
    rethrow()
  end
  isfile(resolved) || throw(ArgumentError("Resolved MATPOWER case file not found: $(resolved)"))
  return _canonical_matpower_source_for_webui(resolved, trusted_directory)
end

"""
    start_powerflow_run(request::AbstractDict; case_directory=nothing) -> Dict{String,Any}

Start one local PowerFlow service run above [`run_sparlectra_api`](@ref). The
request must provide `casefile`, `config_file`, and `output_root`; optional
`config_overrides` are forwarded to the programmatic API. A unique run ID is
chosen before execution, and all generated files are written beneath
`output_root/run_id`. Completed API runs, including failed runs with a
`result.json`, are registered in memory and in the persistent run index. Public
service failures are returned as structured dictionaries.

When `case_directory` is provided by a trusted caller, bare `.m` case names are
resolved there through [`ensure_casefile`](@ref) and remain the executed source.
Generated Julia cache files in that directory are hidden from the Web UI service
path: explicit `.jl` requests resolve to a matching `.m` source when present and
are rejected otherwise.
"""
function start_powerflow_run(request::AbstractDict; case_directory::Union{Nothing,AbstractString} = nothing, case_resolver = _resolve_powerflow_casefile)::Dict{String,Any}
  request_start = time_ns()
  phases = Dict{Symbol,Float64}()
  cancellation_token = _service_request_value(request, "cancellation_token", nothing)
  phase_callback = _service_request_value(request, "phase_callback", phase -> nothing)
  operation_callback = _service_request_value(request, "operation_callback", (event; fields...) -> nothing)
  set_phase = phase -> (phase_callback(String(phase)); _check_powerflow_cancelled!(cancellation_token))
  set_phase("resolving_case")
  _check_powerflow_cancelled!(cancellation_token)
  casefile = _service_request_value(request, "casefile")
  casefile isa AbstractString && !isempty(strip(casefile)) || return _service_failure("missing_casefile", "PowerFlow service request requires a nonempty casefile.")
  _is_for002_reference_dat(casefile) && return _service_failure("invalid_casefile", _for002_primary_case_message())
  if case_directory === nothing
    _check_powerflow_cancelled!(cancellation_token)
    isfile(casefile) || return _service_failure("missing_casefile", "MATPOWER case file not found: $(abspath(casefile))")
    casefile = abspath(casefile)
  else
    resolution_start = time_ns()
    _check_powerflow_cancelled!(cancellation_token)
    set_phase("checking_case_cache")
    casefile = try
      set_phase("preparing_case_file")
      case_resolver(casefile, case_directory)
    catch err
      return _service_failure("invalid_casefile", sprint(showerror, err))
    end
    phases[:case_resolution] = _api_elapsed_seconds(resolution_start)
    _check_powerflow_cancelled!(cancellation_token)
  end
  set_phase("preparing_configuration")

  config_file = _service_request_value(request, "config_file")
  config_file isa AbstractString && !isempty(strip(config_file)) || return _service_failure("missing_config_file", "PowerFlow service request requires a nonempty config_file.")
  isfile(config_file) || return _service_failure("missing_config_file", "Configuration file not found: $(abspath(config_file))")

  output_root = _service_request_value(request, "output_root")
  output_root isa AbstractString && !isempty(strip(output_root)) || return _service_failure("invalid_request", "PowerFlow service request requires a nonempty output_root.")

  config_overrides = _service_request_value(request, "config_overrides", Dict{String,Any}())
  config_overrides isa AbstractDict || return _service_failure("invalid_request", "config_overrides must be dictionary-like.")
  config_override_source = String(_service_request_value(request, "config_override_source", "explicit_api_request"))
  config_override_source in ("explicit_api_request", "webui_form_runtime", "case_sidecar", "user_yaml") || return _service_failure("invalid_request", "config_override_source has unsupported value: $(config_override_source)")
  performance_timing = _service_request_value(request, "performance_timing", :off)
  case_format = _service_request_value(request, "case_format", :auto)
  for002_reference_file = _service_request_value(request, "for002_reference_file", nothing)
  if for002_reference_file isa AbstractString && !isempty(strip(for002_reference_file)) && case_directory !== nothing && !isfile(for002_reference_file) && !occursin(r"[\\/]", for002_reference_file)
    cache_for002 = joinpath(abspath(case_directory), strip(for002_reference_file))
    isfile(cache_for002) && (for002_reference_file = abspath(cache_for002))
  end
  run_dtf_outages = _service_request_value(request, "run_dtf_outages", false)
  run_dtf_outages isa Bool || return _service_failure("invalid_request", "run_dtf_outages must be boolean.")
  dtf_outage_selection = _service_request_value(request, "dtf_outage_selection", String[])
  dtf_outage_selection_mode = _service_request_value(request, "dtf_outage_selection_mode", :none)
  compare_for002_outages = _service_request_value(request, "compare_for002_outages", false)
  compare_for002_outages isa Bool || return _service_failure("invalid_request", "compare_for002_outages must be boolean.")
  write_outage_artifacts = _service_request_value(request, "write_outage_artifacts", true)
  write_outage_artifacts isa Bool || return _service_failure("invalid_request", "write_outage_artifacts must be boolean.")
  write_outage_matpower_exports = _service_request_value(request, "write_outage_matpower_exports", false)
  write_outage_matpower_exports isa Bool || return _service_failure("invalid_request", "write_outage_matpower_exports must be boolean.")
  matpower_export_requested = _service_request_value(request, "matpower_export_requested", false)
  matpower_export_requested isa Bool || return _service_failure("invalid_request", "matpower_export_requested must be boolean.")
  run_diagnostics = _service_request_value(request, "run_diagnostics", false)
  run_diagnostics isa Bool || return _service_failure("invalid_request", "run_diagnostics must be boolean.")
  diagnose_mode = _service_request_value(request, "diagnose_mode", false)
  diagnose_mode isa Bool || return _service_failure("invalid_request", "diagnose_mode must be boolean.")
  diagnose_mode && (run_diagnostics = true)
  short_circuit_mode = _service_request_value(request, "short_circuit_mode", false)
  short_circuit_mode isa Bool || return _service_failure("invalid_request", "short_circuit_mode must be boolean.")
  (short_circuit_mode && diagnose_mode) && return _service_failure("invalid_request", "short_circuit_mode and diagnose_mode are mutually exclusive.")
  import_analysis_mode = _service_request_value(request, "import_analysis_mode", false)
  import_analysis_mode isa Bool || return _service_failure("invalid_request", "import_analysis_mode must be boolean.")
  (import_analysis_mode && (diagnose_mode || short_circuit_mode)) && return _service_failure("invalid_request", "import_analysis_mode excludes diagnose_mode and short_circuit_mode.")
  contingency_mode = _service_request_value(request, "contingency_mode", false)
  contingency_mode isa Bool || return _service_failure("invalid_request", "contingency_mode must be boolean.")
  (contingency_mode && (diagnose_mode || short_circuit_mode || import_analysis_mode)) && return _service_failure("invalid_request", "contingency_mode excludes diagnose_mode, short_circuit_mode and import_analysis_mode.")
  # outage kind is a RUN parameter (branch vs generator N-1), not a config key
  contingency_kind = _service_request_value(request, "contingency_kind", "branch")
  (contingency_kind isa AbstractString && contingency_kind in ("branch", "gen")) || return _service_failure("invalid_request", "contingency_kind must be \"branch\" or \"gen\".")
  detailed_result_csv = _service_request_value(request, "detailed_result_csv", false)
  detailed_result_csv isa Bool || return _service_failure("invalid_request", "detailed_result_csv must be boolean.")
  export_cgmes = _service_request_value(request, "export_cgmes", false)
  export_cgmes isa Bool || return _service_failure("invalid_request", "export_cgmes must be boolean.")
  detailed_result_csv_semicolon = _service_request_value(request, "detailed_result_csv_semicolon", false)
  detailed_result_csv_semicolon isa Bool || return _service_failure("invalid_request", "detailed_result_csv_semicolon must be boolean.")
  detailed_result_csv_format = _service_request_value(request, "detailed_result_csv_format", nothing)
  if detailed_result_csv && detailed_result_csv_format !== nothing
    detailed_result_csv_format isa AbstractString || return _service_failure("invalid_request", "detailed_result_csv_format must be a string.")
    try
      _resolve_detailed_csv_format(detailed_result_csv_format)
    catch err
      return _service_failure("invalid_request", sprint(showerror, err))
    end
  end
  phases[:request_parse] = _api_elapsed_seconds(request_start)

  requested_run_id = _service_request_value(request, "run_id", nothing)
  run_id = requested_run_id === nothing ? string(uuid4()) : String(requested_run_id)
  _safe_powerflow_run_id(run_id) || return _service_failure("unsafe_run_id", "Unsafe PowerFlow run ID rejected.")
  root = abspath(output_root)
  output_dir = joinpath(root, run_id)
  if diagnose_mode
    # Diagnose mode forces the fixed-reference self-check settings (see
    # `_self_check_forced_overrides`) on top of the requested config_file. Some
    # of those settings (e.g. start_projection) are not GUI-editable, so they
    # cannot go through `config_overrides`; the merged result is written into
    # the run's own output directory both as the config Sparlectra actually
    # runs with and as a visible, self-documenting run artifact.
    mkpath(output_dir)
    config_file = try
      self_check_path = joinpath(output_dir, "diagnose_self_check_config.yaml")
      _write_yaml_file(self_check_path, _merge_config_overrides(load_yaml_dict(config_file), _self_check_forced_overrides()))
      self_check_path
    catch err
      return _service_failure("invalid_configuration", sprint(showerror, err, catch_backtrace()); run_id = run_id)
    end
    # The Web UI form submits cgmes_start_values as a GUI-editable override;
    # applied on top of the merged self-check config it would re-flatten the
    # start (default "flat") and break the fixed-reference contract. The
    # self-check's forced start_values=sv must win, so drop the form value.
    if config_overrides isa AbstractDict && haskey(config_overrides, "cgmes_import.start_values")
      config_overrides = Dict{String,Any}(k => v for (k, v) in config_overrides if k != "cgmes_import.start_values")
    end
  end
  # Import-analysis runs stop even earlier than short-circuit runs: only the
  # CGMES delivery files are parsed (no mapping, no solve) and the
  # importFailureAnalysis report becomes the run's product — a fast pre-check
  # before a slow full import (see _run_import_analysis_service).
  if import_analysis_mode
    ia_result = try
      _run_import_analysis_service(casefile, config_file, output_dir, run_id)
    catch err
      err isa PowerFlowAborted && rethrow()
      return _service_failure("execution_error", sprint(showerror, err, catch_backtrace()); run_id = run_id)
    end
    try
      lock(_POWERFLOW_SERVICE_LOCK) do
        _POWERFLOW_SERVICE_RUNS[ia_result.run_id] = ia_result
        _write_powerflow_run_index!(root, ia_result)
      end
    catch err
      lock(_POWERFLOW_SERVICE_LOCK) do
        delete!(_POWERFLOW_SERVICE_RUNS, ia_result.run_id)
      end
      return _service_failure("run_index_error", sprint(showerror, err, catch_backtrace()); run_id = run_id)
    end
    return to_dict(ia_result)
  end
  # Short-circuit runs bypass the power-flow pipeline entirely: CGMES import
  # + runShortCircuit! max/min + CSV artifacts, same result/registry
  # conventions (see _run_short_circuit_service). No PF solve is involved.
  if short_circuit_mode
    sc_result = try
      _run_short_circuit_service(casefile, config_file, output_dir, run_id)
    catch err
      err isa PowerFlowAborted && rethrow()
      return _service_failure("execution_error", sprint(showerror, err, catch_backtrace()); run_id = run_id)
    end
    try
      lock(_POWERFLOW_SERVICE_LOCK) do
        _POWERFLOW_SERVICE_RUNS[sc_result.run_id] = sc_result
        _write_powerflow_run_index!(root, sc_result)
      end
    catch err
      lock(_POWERFLOW_SERVICE_LOCK) do
        delete!(_POWERFLOW_SERVICE_RUNS, sc_result.run_id)
      end
      return _service_failure("run_index_error", sprint(showerror, err, catch_backtrace()); run_id = run_id)
    end
    return to_dict(sc_result)
  end
  # N-1 contingency runs (#331 Phase 5) also bypass the single power-flow
  # pipeline: build the net, run the batch for the requested outage kind, write
  # the CSV + report artifacts (see _run_contingency_service). Same result and
  # registry conventions as the short-circuit path.
  if contingency_mode
    # a per-case weight list lives next to the resolved case as
    # <stem>.contingency-weights.csv (issue #331 Phase 5 follow-up). Its PRESENCE
    # is the switch: pass the path only when the file exists, no request key.
    ct_weights_path = try
      wp = _webui_case_weights_path(casefile)
      isfile(wp) ? wp : nothing
    catch
      nothing
    end
    ct_result = try
      _run_contingency_service(casefile, config_file, output_dir, run_id, contingency_kind; weights_path = ct_weights_path)
    catch err
      err isa PowerFlowAborted && rethrow()
      return _service_failure("execution_error", sprint(showerror, err, catch_backtrace()); run_id = run_id)
    end
    try
      lock(_POWERFLOW_SERVICE_LOCK) do
        _POWERFLOW_SERVICE_RUNS[ct_result.run_id] = ct_result
        _write_powerflow_run_index!(root, ct_result)
      end
    catch err
      lock(_POWERFLOW_SERVICE_LOCK) do
        delete!(_POWERFLOW_SERVICE_RUNS, ct_result.run_id)
      end
      return _service_failure("run_index_error", sprint(showerror, err, catch_backtrace()); run_id = run_id)
    end
    return to_dict(ct_result)
  end
  # Phase timings collected before the API handoff become service metadata, not
  # operation-log events for every internal solver step.
  result = try
    _run_sparlectra_api(
      casefile = casefile,
      config_file = config_file,
      output_dir = output_dir,
      case_format = case_format,
      for002_reference_file = for002_reference_file,
      run_dtf_outages = run_dtf_outages,
      dtf_outage_selection = dtf_outage_selection,
      dtf_outage_selection_mode = dtf_outage_selection_mode,
      compare_for002_outages = compare_for002_outages,
      write_outage_artifacts = write_outage_artifacts,
      write_outage_matpower_exports = write_outage_matpower_exports,
      matpower_export_requested = matpower_export_requested,
      config_overrides = config_overrides,
      config_override_source = config_override_source,
      performance_timing = performance_timing,
      run_diagnostics = run_diagnostics,
      detailed_result_csv = detailed_result_csv,
      detailed_result_csv_format = detailed_result_csv_format,
      detailed_result_csv_semicolon = detailed_result_csv_semicolon,
      export_cgmes = export_cgmes,
      phase_timings = phases,
      run_id = run_id,
      cancellation_token = cancellation_token,
      phase_callback = phase_callback,
      operation_callback = operation_callback,
    )
  catch err
    err isa PowerFlowAborted && rethrow()
    return _service_failure("execution_error", sprint(showerror, err, catch_backtrace()); run_id = run_id)
  end
  if diagnose_mode
    # Same self-check summary the programmatic run_fixed_reference_self_check
    # writes: which start machinery was forced off, the start-state residual,
    # and the CGMES SV-coverage caveat — next to diagnose.log.
    try
      _write_self_check_summary(result, result.output_dir)
    catch err
      @warn "could not write self_check.log" exception = err
    end
  end

  try
    lock(_POWERFLOW_SERVICE_LOCK) do
      _POWERFLOW_SERVICE_RUNS[result.run_id] = result
      _write_powerflow_run_index!(root, result)
    end
  catch err
    lock(_POWERFLOW_SERVICE_LOCK) do
      delete!(_POWERFLOW_SERVICE_RUNS, result.run_id)
    end
    return _service_failure("run_index_error", sprint(showerror, err, catch_backtrace()); run_id = run_id)
  end
  return to_dict(result)
end

"""
    get_powerflow_result(run_id::AbstractString) -> Dict{String,Any}

Return serialized metadata for a registered local PowerFlow service run.

If the run ID is unknown, return a structured service failure dictionary instead
of throwing, so Web UI callers can render a stable error response.
"""
function get_powerflow_result(run_id::AbstractString)::Dict{String,Any}
  result = _registered_powerflow_run(run_id)
  result === nothing && return _service_failure("run_not_found", "No PowerFlow run found for run_id $(run_id)."; run_id = run_id)
  return to_dict(result)
end

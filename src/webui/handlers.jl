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

struct SparlectraWebUIResponse
  status::Int
  headers::Vector{Pair{String,String}}
  body::Vector{UInt8}
end

function SparlectraWebUIResponse(status::Integer, body::AbstractString; content_type::AbstractString = "text/html; charset=utf-8", headers = Pair{String,String}[])
  response_headers = Pair{String,String}["Content-Type" => String(content_type)]
  append!(response_headers, headers)
  return SparlectraWebUIResponse(Int(status), response_headers, Vector{UInt8}(codeunits(body)))
end

function _webui_html(body::AbstractString; status::Integer = 200)
  return SparlectraWebUIResponse(status, body)
end

function _webui_redirect(location::AbstractString)
  return SparlectraWebUIResponse(303, ""; headers = ["Location" => String(location)])
end

const _WEBUI_LOGO_PATH = normpath(joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "logo.png"))

function handle_webui_logo()::SparlectraWebUIResponse
  isfile(_WEBUI_LOGO_PATH) || return _webui_html(render_webui_error(404, "Sparlectra.jl logo asset is unavailable."); status = 404)
  return SparlectraWebUIResponse(200, ["Content-Type" => "image/png"], read(_WEBUI_LOGO_PATH))
end

"""Run a PowerFlow request through the Web UI form-to-service boundary."""
function handle_powerflow_run(form::AbstractDict; default_output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, runner = start_powerflow_run, operation_log::AbstractString = default_output_root)::Dict{String,Any}
  request = powerflow_webui_request(form; default_output_root = default_output_root)
  package_case_directory = joinpath(application_root, "data", "mpower")
  requested_case = String(request["casefile"])
  if !isabspath(requested_case) && !occursin('/', requested_case) && !occursin('\\', requested_case)
    cached_case = case_directory === nothing ? "" : joinpath(case_directory, requested_case)
    package_case = joinpath(package_case_directory, requested_case)
    if isfile(cached_case)
      request["casefile"] = cached_case
    elseif isfile(package_case) && case_directory !== nothing
      mkpath(case_directory)
      cp(package_case, cached_case; force = true)
      request["casefile"] = cached_case
    elseif isfile(package_case)
      request["casefile"] = package_case
    end
  end
  effective_case_directory = case_directory === nothing ? package_case_directory : String(case_directory)
  event_callback = (event; fields...) -> record_webui_operation!(operation_log, event; route = "/powerflow/run", method = "POST", user_action = false, fields...)
  return start_webui_powerflow_run(request; case_directory = effective_case_directory, runner, event_callback)
end


function _config_refresh_result_dict(result; config_file::AbstractString = "", downloadable::Bool = false)::Dict{String,Any}
  return Dict{String,Any}(
    "success" => result.success,
    "changed" => result.changed,
    "written" => result.written,
    "backup_path" => result.backup_path === nothing ? "" : String(result.backup_path),
    "missing_keys" => result.missing_keys,
    "normalized_keys" => result.normalized_keys,
    "duplicate_keys" => result.duplicate_keys,
    "warnings" => result.warnings,
    "refreshed_text" => result.refreshed_text,
    "config_file" => String(config_file),
    "downloadable" => downloadable,
  )
end

function handle_powerflow_config_refresh(form::AbstractDict; write::Bool = false, operation_log::AbstractString = "results/powerflow_service")::SparlectraWebUIResponse
  event_prefix = write ? "config_refresh_write" : "config_refresh_check"
  config_text = strip(String(something(_webui_form_value(form, "config_text", ""), "")))
  config_file = strip(String(something(_webui_form_value(form, "config_file", ""), "")))
  record_webui_operation!(operation_log, string(event_prefix, "_started"); route = "/powerflow/config", method = "POST", user_action = true, config_file)
  try
    if !isempty(config_text)
      result = refresh_sparlectra_config_text(config_text)
      data = _config_refresh_result_dict(result; config_file = "browser upload", downloadable = true)
      record_webui_operation!(operation_log, string(event_prefix, "_completed"); route = "/powerflow/config", method = "POST", user_action = true, config_file = "browser upload", changed = result.changed, written = false, backup_path = nothing, missing_key_count = length(result.missing_keys), normalized_key_count = length(result.normalized_keys), duplicate_key_count = length(result.duplicate_keys), message = "refreshed YAML available for download")
      return _webui_html(render_config_refresh_result(data))
    end
    isempty(config_file) && throw(ArgumentError("No configuration file was provided."))
    if write && !isfile(config_file)
      record_webui_operation!(operation_log, "config_refresh_write_rejected"; route = "/powerflow/config", method = "POST", user_action = true, config_file, message = "configuration file is not writable in place")
      return _webui_html(render_webui_error(400, "Configuration refresh writes require a server-local file. Use a server-local selected configuration file for in-place refresh."); status = 400)
    end
    result = refresh_sparlectra_config_file(config_file; write)
    data = _config_refresh_result_dict(result; config_file)
    event = result.success ? string(event_prefix, "_completed") : "config_refresh_write_rejected"
    record_webui_operation!(operation_log, event; route = "/powerflow/config", method = "POST", user_action = true, config_file, changed = result.changed, written = result.written, backup_path = result.backup_path, missing_key_count = length(result.missing_keys), normalized_key_count = length(result.normalized_keys), duplicate_key_count = length(result.duplicate_keys), message = result.success ? "configuration refresh completed" : "duplicate keys require manual review")
    return _webui_html(render_config_refresh_result(data); status = result.success ? 200 : 400)
  catch err
    record_webui_operation!(operation_log, "config_refresh_failed"; route = "/powerflow/config", method = "POST", user_action = true, config_file, message = sprint(showerror, err))
    return _webui_html(render_webui_error(400, sprint(showerror, err)); status = 400)
  end
end
function handle_powerflow_result(run_id::AbstractString)::SparlectraWebUIResponse
  result = get_webui_powerflow_job(run_id)
  status = get(result, "success", false) || get(result, "reason", "") != "run_not_found" ? 200 : 404
  return _webui_html(render_powerflow_result(result); status = status)
end

function _webui_case_profile_scalar(field::AbstractString, value)
  value === nothing && return nothing
  value === missing && return nothing
  value isa Bool && return value
  value isa Integer && return Int(value)
  value isa AbstractFloat && return isfinite(value) ? Float64(value) : throw(ArgumentError("Case-settings field $(field) must be finite."))
  value isa Symbol && return String(value)
  value isa AbstractString && return String(value)
  throw(ArgumentError("Case-settings field $(field) has unsupported value type $(typeof(value))."))
end

function _webui_case_profile_value(field::AbstractString, value)
  value isa AbstractVector || return _webui_case_profile_scalar(field, value)
  return [_webui_case_profile_scalar(field, item) for item in value if !(item === nothing || item === missing)]
end

function _webui_case_profile_setting!(settings::Dict{String,Any}, field::AbstractString, value)
  field in _WEBUI_CASE_PROFILE_FIELDS || throw(ArgumentError("Case-settings field $(field) is not allowed."))
  serialized = _webui_case_profile_value(field, value)
  serialized === nothing && return settings
  settings[field] = serialized
  return settings
end

const _WEBUI_CASE_PROFILE_CONFIG_FIELD_BY_KEY = Dict{String,String}(config_key => field for (config_key, field, _) in _WEBUI_FORM_CONFIG_FIELDS)

function _webui_case_profile_settings(settings_raw::AbstractDict)::Dict{String,Any}
  settings = Dict{String,Any}()
  for (raw_key, raw_value) in settings_raw
    key = String(raw_key)
    field = get(_WEBUI_CASE_PROFILE_CONFIG_FIELD_BY_KEY, key, key)
    field in _WEBUI_CASE_PROFILE_FIELDS || continue
    _webui_case_profile_setting!(settings, field, raw_value)
  end
  return settings
end

function _webui_case_settings_saved_html(path::AbstractString, count::Integer, successful::Bool, override::Bool)::String
  status = successful ? "successful/converged" : (override ? "non-successful, saved via explicit override" : "non-successful")
  body = """
<section class=\"panel case-settings-saved\"><h2>Case settings saved</h2>
<p class=\"alert info\"><strong>Saved sidecar:</strong> <code>$(_webui_escape(path))</code></p>
<ul>
<li><strong>Saved settings:</strong> $(count)</li>
<li><strong>Run status:</strong> $(_webui_escape(status))</li>
</ul>
<p>These settings will be applied when the same case is selected again. Manual edits in the form still override the saved profile for each run.</p>
<p><a href=\"/powerflow?casefile=$(_webui_urlencode(path))\">Open the run form with this case</a></p>
</section>"""
  return _webui_layout("Case settings saved", body; show_back = true)
end

function handle_powerflow_case_settings_save(run_id::AbstractString, form::AbstractDict; output_root::AbstractString, operation_log::AbstractString)::SparlectraWebUIResponse
  result = get_webui_powerflow_job(run_id)
  if get(result, "reason", "") == "run_not_found"
    record_webui_operation!(operation_log, "case_settings_save_failed"; route = "/powerflow/result/$(run_id)/case-settings/save", method = "POST", user_action = true, run_id, status = "rejected", message = "run not found")
    return _webui_html(render_webui_error(404, "PowerFlow run not found."); status = 404)
  end
  successful = _webui_result_successful(result)
  override = _webui_parse_bool(_webui_form_value(form, "override_non_success", false))
  if !successful && !override
    record_webui_operation!(operation_log, "case_settings_save_failed"; route = "/powerflow/result/$(run_id)/case-settings/save", method = "POST", user_action = true, run_id, status = "rejected", message = "non-successful run requires explicit override")
    return _webui_html(render_webui_error(400, "Saving settings from a non-successful run requires the explicit override action."); status = 400)
  end
  metadata = get(result, "metadata", Dict{String,Any}())
  runtime_casefile = ""
  for source in (get(result, "runtime_casefile", nothing),
                 metadata isa AbstractDict ? get(metadata, "runtime_casefile", nothing) : nothing,
                 get(result, "resolved_casefile", nothing),
                 get(result, "casefile", nothing),
                 metadata isa AbstractDict ? get(metadata, "runtime_casefile_path", nothing) : nothing)
    source === nothing && continue
    candidate = String(source)
    if !isempty(candidate)
      runtime_casefile = candidate
      break
    end
  end
  runtime_casefile_path = String(something(
    get(result, "resolved_casefile", nothing),
    get(result, "casefile", nothing),
    metadata isa AbstractDict ? get(metadata, "runtime_casefile_path", nothing) : nothing,
    runtime_casefile,
  ))
  settings_raw = metadata isa AbstractDict ? get(metadata, "webui_request_settings", nothing) : nothing
  if isempty(runtime_casefile) || !(settings_raw isa AbstractDict)
    record_webui_operation!(operation_log, "case_settings_save_failed"; route = "/powerflow/result/$(run_id)/case-settings/save", method = "POST", user_action = true, run_id, status = "rejected", message = "run metadata incomplete")
    return _webui_html(render_webui_error(400, "Run metadata is incomplete; case settings were not saved."); status = 400)
  end
  settings = try
    _webui_case_profile_settings(settings_raw)
  catch err
    message = err isa ArgumentError ? sprint(showerror, err) : "Case-settings profile contains an unsupported value."
    record_webui_operation!(operation_log, "case_settings_save_failed"; route = "/powerflow/result/$(run_id)/case-settings/save", method = "POST", user_action = true, run_id, status = "rejected", message)
    return _webui_html(render_webui_error(400, message); status = 400)
  end
  isempty(settings) && return _webui_html(render_webui_error(400, "No Web UI settings were recorded for this run."); status = 400)
  case_settings_source = !isempty(runtime_casefile_path) ? runtime_casefile_path : runtime_casefile
  key = _webui_normalized_case_key(case_settings_source)
  path = try
    _webui_case_settings_path(output_root, case_settings_source)
  catch err
    message = err isa ArgumentError ? sprint(showerror, err) : "Unsafe MATPOWER case path for case-settings profile."
    record_webui_operation!(operation_log, "case_settings_save_failed"; route = "/powerflow/result/$(run_id)/case-settings/save", method = "POST", user_action = true, run_id, status = "rejected", message)
    return _webui_html(render_webui_error(400, message); status = 400)
  end
  isfile(case_settings_source) || return _webui_html(render_webui_error(400, "The runtime MATPOWER case file is not available; case settings were not saved."); status = 400)
  profile = Dict{String,Any}(
    "schema_version" => 1,
    "profile_kind" => "webui_case_settings",
    "case" => Dict{String,Any}(
      "display_name" => basename(runtime_casefile),
      "profile_path" => path,
      "normalized_case_key" => key,
      "source" => "webui_mpower_data",
      "original_path" => runtime_casefile_path,
    ),
    "saved_from_run" => Dict{String,Any}(
      "run_id" => run_id,
      "status" => successful ? "converged" : string(get(result, "status", "not_converged")),
      "saved_at" => Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ"),
      "override_non_success" => !successful && override,
    ),
    "settings" => settings,
  )
  try
    _write_yaml_file(path, profile)
    record_webui_operation!(operation_log, "case_settings_saved"; route = "/powerflow/result/$(run_id)/case-settings/save", method = "POST", user_action = true, run_id, status = "succeeded", profile_path = path, normalized_case_key = key)
    return _webui_html(_webui_case_settings_saved_html(path, length(settings), successful, !successful && override))
  catch err
    record_webui_operation!(operation_log, "case_settings_save_failed"; route = "/powerflow/result/$(run_id)/case-settings/save", method = "POST", user_action = true, run_id, status = "failed", message = sprint(showerror, err))
    return _webui_html(render_webui_error(500, sprint(showerror, err)); status = 500)
  end
end

function handle_powerflow_abort(run_id::AbstractString)::SparlectraWebUIResponse
  result = abort_webui_powerflow_run(run_id)
  status = get(result, "reason", "") in ("unsafe_run_id", "run_not_found") ? 404 : 303
  status == 303 && return _webui_redirect("/powerflow/result/$(_webui_urlencode(run_id))")
  return _webui_html(render_webui_error(status, get(result, "message", "Abort request failed.")); status)
end

function handle_powerflow_hard_reset(run_id::AbstractString)::SparlectraWebUIResponse
  result = hard_reset_webui_powerflow_run(run_id)
  haskey(result, "reason") && return _webui_html(render_webui_error(409, get(result, "message", "Hard reset failed.")); status = 409)
  return _webui_html(render_webui_hard_reset())
end

function handle_powerflow_artifacts(run_id::AbstractString)::SparlectraWebUIResponse
  artifacts = list_powerflow_artifacts(run_id)
  status = artifacts isa AbstractDict ? 404 : 200
  return _webui_html(render_powerflow_artifacts(run_id, artifacts); status = status)
end

const _WEBUI_TEXT_MIME_TYPES = Set(("text/plain", "application/json", "application/x-yaml", "text/csv", "text/html", "text/markdown"))
const _WEBUI_ARTIFACT_PREVIEW_BYTES = 64 * 1024

function _read_webui_artifact_preview(path::AbstractString; max_bytes::Integer = _WEBUI_ARTIFACT_PREVIEW_BYTES)::NamedTuple
  open(path, "r") do io
    bytes = read(io, max_bytes + 1)
    truncated = length(bytes) > max_bytes
    truncated && resize!(bytes, max_bytes)
    return (content = String(bytes), truncated = truncated)
  end
end

function handle_powerflow_artifact(run_id::AbstractString, artifact_name::AbstractString)::SparlectraWebUIResponse
  artifact = resolve_powerflow_artifact(run_id, artifact_name)
  if artifact isa AbstractDict
    reason = get(artifact, "reason", "artifact_error")
    status = reason in ("unsafe_artifact_name", "artifact_not_found") ? 400 : 404
    return _webui_html(render_webui_error(status, get(artifact, "message", reason)); status = status)
  end
  if artifact.mime_type in _WEBUI_TEXT_MIME_TYPES
    preview = _read_webui_artifact_preview(artifact.path)
    notice = preview.truncated ? "<p class=\"alert warning\">Preview truncated to $(_WEBUI_ARTIFACT_PREVIEW_BYTES) bytes. Use Download for the complete artifact.</p>" : ""
    page = _webui_layout("Artifact: $(artifact.name)", "<section class=\"artifact-text-page\"><p><a class=\"button\" href=\"?download=1\">Download</a></p>$(notice)<pre class=\"artifact-text\">$(_webui_escape(preview.content))</pre></section>"; show_back = true, main_class = "page artifact-page")
    return _webui_html(page)
  end
  bytes = read(artifact.path)
  headers = ["Content-Disposition" => "attachment; filename=\"$(replace(basename(artifact.name), '"' => '_'))\""]
  return SparlectraWebUIResponse(200, Pair{String,String}["Content-Type" => artifact.mime_type; headers], bytes)
end

function handle_powerflow_artifact_download(run_id::AbstractString, artifact_name::AbstractString)::SparlectraWebUIResponse
  artifact = resolve_powerflow_artifact(run_id, artifact_name)
  if artifact isa AbstractDict
    return _webui_html(render_webui_error(400, get(artifact, "message", "Artifact unavailable.")); status = 400)
  end
  headers = Pair{String,String}[
    "Content-Type" => artifact.mime_type,
    "Content-Disposition" => "attachment; filename=\"$(replace(basename(artifact.name), '"' => '_'))\"",
  ]
  return SparlectraWebUIResponse(200, headers, read(artifact.path))
end

function handle_powerflow_history(output_root::AbstractString)::SparlectraWebUIResponse
  return _webui_html(render_powerflow_history(list_powerflow_runs(output_root), output_root; active_run = get_active_webui_powerflow_job()))
end

function handle_webui_operation_log(output_root::AbstractString; download::Bool = false)::SparlectraWebUIResponse
  path = webui_operation_log_path(output_root)
  content = isfile(path) ? read(path, String) : ""
  if download
    headers = ["Content-Disposition" => "attachment; filename=\"$(WEBUI_OPERATION_LOG_FILENAME)\""]
    return SparlectraWebUIResponse(200, content; content_type = "application/x-ndjson; charset=utf-8", headers)
  end
  return _webui_html(render_webui_operation_log(content))
end

function handle_powerflow_refresh(output_root::AbstractString)::Dict{String,Any}
  return refresh_powerflow_run_registry!(output_root)
end

function handle_powerflow_delete(run_id::AbstractString, output_root::AbstractString)::SparlectraWebUIResponse
  job = get_webui_powerflow_job(run_id)
  if get(job, "status", "") in _POWERFLOW_WEBUI_ACTIVE_STATES
    return _webui_html(render_webui_error(409, "This run is still active. Abort it first and wait until it reaches aborted status."); status = 409)
  end
  result = delete_powerflow_run(run_id; output_root)
  get(result, "success", false) && return _webui_redirect("/powerflow/history")
  status = get(result, "reason", "") in ("unsafe_run_id", "unsafe_output_dir", "unsafe_result_file") ? 400 : 404
  return _webui_html(render_webui_error(status, get(result, "message", "Run deletion failed.")); status)
end

function handle_powerflow_delete_all(output_root::AbstractString)::SparlectraWebUIResponse
  result = delete_all_powerflow_runs(; output_root)
  get(result, "success", false) && return _webui_redirect("/powerflow/history")
  message = "Some runs could not be deleted: " * join((string(get(item, "run_id", "unknown"), " (", get(item, "reason", "delete_failed"), ")") for item in result["failed_runs"]), ", ")
  return _webui_html(render_webui_error(500, message); status = 500)
end

function handle_webui_help(topic::AbstractString)::SparlectraWebUIResponse
  metadata = resolve_webui_help_topic(topic)
  metadata === nothing && return _webui_html(render_webui_error(404, "Unknown help topic."); status = 404)
  excerpt = load_webui_help_excerpt(topic)
  excerpt === nothing && return _webui_html(render_webui_error(404, "No help section found for this option."); status = 404)
  return _webui_html(render_webui_help(metadata, excerpt))
end

function handle_webui_docs_index()::SparlectraWebUIResponse
  return _webui_html(render_webui_docs_index(WEBUI_DOC_PAGES))
end

function handle_webui_doc_page(page::AbstractString)::SparlectraWebUIResponse
  metadata = resolve_webui_doc_page(page)
  metadata === nothing && return _webui_html(render_webui_error(404, "Documentation page not found."); status = 404)
  markdown_text = load_webui_markdown_document(page)
  markdown_text === nothing && return _webui_html(render_webui_error(404, "Documentation page is unavailable."); status = 404)
  return _webui_html(render_webui_doc_page(page, metadata, markdown_text))
end

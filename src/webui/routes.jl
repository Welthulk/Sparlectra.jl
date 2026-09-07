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

# file: src/webui/routes.jl
# purpose: Web UI request routing (route_sparlectra_webui) with URL decoding,
#          query and form pair parsing, and configuration notices
function _webui_urldecode(value::AbstractString)::String
  bytes = UInt8[]
  text = String(value)
  index = firstindex(text)
  while index <= lastindex(text)
    char = text[index]
    if char == '%'
      next1 = nextind(text, index)
      next2 = nextind(text, next1)
      next2 <= lastindex(text) || throw(ArgumentError("Invalid URL encoding."))
      push!(bytes, parse(UInt8, text[next1:next2]; base = 16))
      index = nextind(text, next2)
    elseif char == '+'
      push!(bytes, UInt8(' '))
      index = nextind(text, index)
    else
      append!(bytes, codeunits(string(char)))
      index = nextind(text, index)
    end
  end
  return String(bytes)
end

function _webui_parse_pairs(text::AbstractString)::Dict{String,String}
  values = Dict{String,String}()
  isempty(text) && return values
  for pair in split(text, '&')
    parts = split(pair, '='; limit = 2)
    key = _webui_urldecode(parts[1])
    value = length(parts) == 2 ? _webui_urldecode(parts[2]) : ""
    values[key] = value
  end
  return values
end

function _webui_split_target(target::AbstractString)
  parts = split(String(target), '?'; limit = 2)
  return parts[1], length(parts) == 2 ? _webui_parse_pairs(parts[2]) : Dict{String,String}()
end

function _powerflow_config_notice(config_file::AbstractString)
  isempty(strip(config_file)) && return nothing
  isfile(config_file) || return nothing
  try
    result = refresh_sparlectra_config_file(config_file; write = false)
    return (result.changed || !isempty(result.duplicate_keys)) ? _config_refresh_result_dict(result; config_file) : nothing
  catch
    return nothing
  end
end

# Fails open (notice stays visible) on any load error, matching
# _powerflow_config_notice's fail-safe pattern — a broken/missing config file
# should not silently suppress an otherwise-informative notice.
function _powerflow_show_case_settings_notice(config_file::AbstractString)::Bool
  isempty(strip(config_file)) && return true
  isfile(config_file) || return true
  try
    return load_sparlectra_config(config_file; reload = true).webui.show_case_settings_notice
  catch
    return true
  end
end

# One-click dismiss for the "Case-specific settings loaded from ..." notice:
# persists webui.show_case_settings_notice = false into the given config file
# (a nested-dict merge, so other keys/comments in the file are left as-is by
# the merge itself; only the final YAML dump loses comments, same tradeoff as
# the existing "Configuration refresh" writer). Requires a server-local,
# already-existing config file — the same guard the config editor uses.
function _powerflow_dismiss_case_settings_notice!(config_file::AbstractString)
  isempty(strip(config_file)) && throw(ArgumentError("No configuration file was provided."))
  isfile(config_file) || throw(ArgumentError("Configuration refresh writes require a server-local file. Use a server-local selected configuration file."))
  merged = _merge_config_overrides(load_yaml_dict(config_file), Dict{String,Any}("webui" => Dict{String,Any}("show_case_settings_notice" => false)))
  _write_yaml_file(config_file, merged)
  return config_file
end

function route_sparlectra_webui(method::AbstractString, target::AbstractString, form::AbstractDict = Dict{String,String}(); output_root::AbstractString = "results/powerflow_service", runtime = nothing)::SparlectraWebUIResponse
  path, query = _webui_split_target(target)
  verb = uppercase(String(method))
  log_root = runtime === nothing ? output_root : runtime.operation_log
  if verb == "GET" && path == "/assets/logo.png"
    return handle_webui_logo()
  elseif verb == "GET" && startswith(path, "/help/")
    return handle_webui_help(_webui_urldecode(path[(lastindex("/help/") + 1):end]))
  elseif verb == "GET" && path == "/docs"
    return handle_webui_docs_index()
  elseif verb == "GET" && startswith(path, "/docs/")
    return handle_webui_doc_page(_webui_urldecode(path[(lastindex("/docs/") + 1):end]))
  elseif verb == "GET" && path in ("/", "/powerflow")
    _webui_log_route!(log_root, "powerflow_form_opened", verb, path; status = "opened")
    # stage 4A harmonization: without an explicit ?casefile the run page uses
    # the case remembered by the Case page (plain nav links carry no query)
    selected_casefile = get(query, "casefile", "")
    selected_casefile = isempty(selected_casefile) ? _webui_recall_selected_case(output_root) : _webui_remember_selected_case!(output_root, selected_casefile)
    case_profile = isempty(selected_casefile) ? nothing : _webui_load_case_settings(output_root, selected_casefile; case_directory = runtime === nothing ? nothing : runtime.case_directory)
    selected_config_file = get(query, "config_file", runtime === nothing ? "" : runtime.config_file)
    return _webui_html(render_powerflow_form(;
      output_root,
      case_directory = runtime === nothing ? nothing : runtime.case_directory,
      operation_log = runtime === nothing ? webui_operation_log_path(output_root) : runtime.operation_log,
      selected_casefile,
      selected_config_file,
      import_message = get(query, "import_message", ""),
      download_file = get(query, "download", ""),
      error_message = runtime === nothing ? nothing : runtime.startup_config_error,
      config_notice = _powerflow_config_notice(runtime === nothing ? "" : runtime.config_file),
      case_profile,
      show_case_settings_notice = _powerflow_show_case_settings_notice(selected_config_file),
      # stage 4A block 4: the SE section on this page reads its query keys
      # (message, sticky g_* generator values) from the page query
      se_query = query,
    ))
  elseif verb == "GET" && path == "/powerflow/case"
    # stage 4A: case management page (chooser, upload, export, import options)
    _webui_log_route!(log_root, "case_page_opened", verb, path; status = "opened")
    # an explicit ?casefile updates the shared selected-case memory; without
    # one the page reopens on the remembered case
    case_page_casefile = get(query, "casefile", "")
    case_page_casefile = isempty(case_page_casefile) ? _webui_recall_selected_case(output_root) : _webui_remember_selected_case!(output_root, case_page_casefile)
    case_page_profile = isempty(case_page_casefile) ? nothing : _webui_load_case_settings(output_root, case_page_casefile; case_directory = runtime === nothing ? nothing : runtime.case_directory)
    return _webui_html(render_case_page(;
      output_root,
      case_directory = runtime === nothing ? nothing : runtime.case_directory,
      operation_log = runtime === nothing ? webui_operation_log_path(output_root) : runtime.operation_log,
      selected_casefile = case_page_casefile,
      selected_config_file = get(query, "config_file", runtime === nothing ? "" : runtime.config_file),
      import_message = get(query, "import_message", ""),
      download_file = get(query, "download", ""),
      case_profile = case_page_profile,
    ))
  elseif verb == "GET" && path == "/powerflow/settings"
    # stage 4A block 3: solver/output/expert options page
    _webui_log_route!(log_root, "settings_page_opened", verb, path; status = "opened")
    settings_casefile = get(query, "casefile", "")
    settings_casefile = isempty(settings_casefile) ? _webui_recall_selected_case(output_root) : _webui_remember_selected_case!(output_root, settings_casefile)
    settings_profile = isempty(settings_casefile) ? nothing : _webui_load_case_settings(output_root, settings_casefile; case_directory = runtime === nothing ? nothing : runtime.case_directory)
    return _webui_html(render_settings_page(;
      output_root,
      case_directory = runtime === nothing ? nothing : runtime.case_directory,
      operation_log = runtime === nothing ? webui_operation_log_path(output_root) : runtime.operation_log,
      selected_casefile = settings_casefile,
      selected_config_file = get(query, "config_file", runtime === nothing ? "" : runtime.config_file),
      save_message = get(query, "save_message", ""),
      case_profile = settings_profile,
      show_case_settings_notice = _powerflow_show_case_settings_notice(get(query, "config_file", runtime === nothing ? "" : runtime.config_file)),
    ))
  elseif verb == "POST" && path == "/powerflow/settings/save"
    return handle_settings_save(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "POST" && path == "/powerflow/case/options/save"
    return handle_case_options_save(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "POST" && path == "/powerflow/run"
    try
      result = handle_powerflow_run(form; default_output_root = output_root, case_directory = runtime === nothing ? nothing : runtime.case_directory, runner = runtime === nothing ? start_powerflow_run : runtime.runner, operation_log = log_root)
      manual_case = strip(String(something(_webui_form_value(form, "casefile_manual", ""), "")))
      requested_case = isempty(manual_case) ? String(something(_webui_form_value(form, "casefile", ""), "")) : manual_case
      if haskey(result, "run_id") && !haskey(result, "reason")
        get(form, "diagnose_mode", nothing) != "true" || _webui_log_route!(log_root, "diagnose_mode_enabled", verb, path; status = "enabled", run_id = result["run_id"])
        get(form, "performance_timing", "off") == "off" || _webui_log_route!(log_root, "performance_timing_enabled", verb, path; status = String(form["performance_timing"]), run_id = result["run_id"])
        return _webui_redirect("/powerflow/result/$(_webui_urlencode(result["run_id"]))")
      end
      _webui_log_route!(log_root, "powerflow_submit_rejected", verb, path; status = "rejected", run_id = get(result, "run_id", nothing), requested_case, message = get(result, "message", nothing))
      return _webui_html(render_powerflow_result(result); status = 400)
    catch err
      selected_casefile = String(something(_webui_form_value(form, "casefile", ""), ""))
      selected_config_file = String(something(_webui_form_value(form, "config_file", ""), ""))
      _webui_log_route!(log_root, "validation_error", verb, path; status = "rejected", requested_case = selected_casefile, message = sprint(showerror, err))
      return _webui_html(render_powerflow_form(;
        output_root,
        operation_log = log_root,
        error_message = sprint(showerror, err),
        selected_casefile,
        selected_config_file,
        submitted_form = form,
      ); status = 400)
    end
  elseif verb == "POST" && path == "/powerflow/import-cases"
    return handle_powerflow_case_import(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "POST" && path == "/powerflow/resolve-case"
    return handle_powerflow_case_resolve(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "POST" && path == "/powerflow/case-settings/reset"
    return handle_powerflow_case_settings_reset(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "POST" && path == "/powerflow/export-scf"
    return handle_powerflow_export_scf(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "POST" && path == "/powerflow/delete-case"
    return handle_powerflow_case_delete(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "GET" && path == "/powerflow/contingency-weights"
    config_file = get(query, "config_file", runtime === nothing ? DEFAULT_SPARLECTRA_CONFIG_PATH : runtime.config_file)
    return handle_contingency_weights_page(query; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, config_file = config_file, operation_log = log_root)
  elseif verb == "GET" && path == "/powerflow/case/download"
    return handle_powerflow_case_download(query; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory)
  elseif verb == "GET" && path == "/powerflow/contingency-weights/download"
    return handle_contingency_weights_download(query; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory)
  elseif verb == "POST" && path == "/powerflow/contingency-weights/upload"
    return handle_contingency_weights_upload(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "POST" && path == "/powerflow/contingency-weights/save"
    return handle_contingency_weights_save(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "POST" && path == "/powerflow/contingency-weights/reset"
    return handle_contingency_weights_reset(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "GET" && path == "/powerflow/scenarios"
    scen_config = get(query, "config_file", runtime === nothing ? DEFAULT_SPARLECTRA_CONFIG_PATH : runtime.config_file)
    return handle_scenarios_page(query; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, config_file = scen_config, operation_log = log_root)
  elseif verb == "POST" && path == "/powerflow/scenarios/save"
    scen_config = runtime === nothing ? DEFAULT_SPARLECTRA_CONFIG_PATH : runtime.config_file
    return handle_scenarios_save(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, config_file = scen_config, operation_log = log_root)
  elseif verb == "POST" && path == "/powerflow/scenarios/delete"
    return handle_scenarios_delete(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "GET" && startswith(path, "/stateestimation/measurements/download")
    return handle_se_measurement_download(query; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory)
  elseif verb == "GET" && path == "/stateestimation"
    # stage 4A block 4: the SE surface lives on the Runs page. This is a
    # REAL redirect, not a rendering alias: two routes rendering the same
    # surface would recreate exactly the divergence stage 4 removes. The
    # query travels along (case becomes the shared casefile key; message
    # and the sticky g_* generator values pass through), and the anchor
    # lands the browser on the SE section.
    se_pairs = String[]
    for (k, v) in query
      key = String(k) == "case" ? "casefile" : String(k)
      push!(se_pairs, string(_webui_urlencode(key), "=", _webui_urlencode(String(v))))
    end
    se_target = isempty(se_pairs) ? "/powerflow" : string("/powerflow?", join(se_pairs, "&"))
    return _webui_redirect(string(se_target, "#state-estimation"))
  elseif verb == "POST" && path == "/stateestimation/generate-measurements"
    return handle_se_generate_measurements(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "POST" && path == "/stateestimation/add-noise"
    return handle_se_add_noise(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "POST" && path == "/stateestimation/reset-settings"
    return handle_se_reset_settings(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "POST" && path == "/stateestimation/measurements/save"
    return handle_se_measurement_save(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "POST" && path == "/stateestimation/measurements/update-values"
    return handle_se_measurement_update_values(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "POST" && path == "/stateestimation/topology-hypotheses"
    return handle_se_topology_hypotheses(form; output_root, application_root = _webui_application_root(), case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = log_root)
  elseif verb == "POST" && path == "/powerflow/config/check"
    return handle_powerflow_config_refresh(form; write = false, operation_log = log_root)
  elseif verb == "POST" && path == "/powerflow/config/refresh"
    return handle_powerflow_config_refresh(form; write = true, operation_log = log_root)
  elseif verb == "POST" && path == "/powerflow/config/download"
    text = String(something(_webui_form_value(form, "refreshed_text", ""), ""))
    return SparlectraWebUIResponse(200, Pair{String,String}["Content-Type" => "application/x-yaml; charset=utf-8", "Content-Disposition" => "attachment; filename=\"configuration-refreshed.yaml\""], Vector{UInt8}(codeunits(text)))
  elseif verb == "GET" && path == "/powerflow/config/edit"
    config_file = get(query, "config_file", runtime === nothing ? DEFAULT_SPARLECTRA_CONFIG_PATH : runtime.config_file)
    return handle_powerflow_config_editor(config_file)
  elseif verb == "POST" && path == "/powerflow/config/edit"
    return handle_powerflow_config_editor_save(form; operation_log = log_root)
  elseif verb == "POST" && path == "/powerflow/config/dismiss-case-settings-notice"
    config_file = String(something(_webui_form_value(form, "config_file", ""), ""))
    casefile = String(something(_webui_form_value(form, "casefile", ""), ""))
    try
      _powerflow_dismiss_case_settings_notice!(config_file)
      _webui_log_route!(log_root, "case_settings_notice_dismissed", verb, path; status = "dismissed", config_file)
    catch err
      _webui_log_route!(log_root, "case_settings_notice_dismiss_failed", verb, path; status = "rejected", config_file, message = sprint(showerror, err))
      return _webui_html(render_webui_error(400, sprint(showerror, err)); status = 400)
    end
    # back to the Settings page: since stage 4A block 3 the notice (and its
    # dismiss button) renders there
    redirect_target = isempty(casefile) ? "/powerflow/settings" : "/powerflow/settings?casefile=$(_webui_urlencode(casefile))"
    return _webui_redirect(redirect_target)
  elseif verb == "GET" && startswith(path, "/powerflow/result/")
    run_id = _webui_urldecode(path[(lastindex("/powerflow/result/") + 1):end])
    get(query, "autorefresh", "") == "1" || _webui_log_route!(log_root, "powerflow_status_opened", verb, path; status = "opened", run_id)
    return handle_powerflow_result(run_id)
  elseif verb == "POST" && startswith(path, "/powerflow/result/") && endswith(path, "/case-settings/save")
    prefix = "/powerflow/result/"
    suffix = "/case-settings/save"
    encoded_run_id = path[(lastindex(prefix) + 1):(lastindex(path) - lastindex(suffix))]
    run_id = _webui_urldecode(encoded_run_id)
    return handle_powerflow_case_settings_save(run_id, form; output_root, operation_log = log_root)
  elseif verb == "POST" && startswith(path, "/powerflow/abort/")
    run_id = _webui_urldecode(path[(lastindex("/powerflow/abort/") + 1):end])
    prior_status = get(get_webui_powerflow_job(run_id), "status", "")
    response = handle_powerflow_abort(run_id)
    event = prior_status == "aborting" ? "powerflow_abort_already_requested" : prior_status == "aborted" ? "powerflow_abort_ignored" : "powerflow_abort_requested"
    event_status = prior_status == "aborting" ? "already_aborting" : prior_status == "aborted" ? "already_aborted" : response.status == 303 ? "accepted" : "rejected"
    current_phase = get(get_webui_powerflow_job(run_id), "current_phase", nothing)
    _webui_log_route!(log_root, event, verb, path; status = event_status, run_id, current_phase)
    return response
  elseif verb == "POST" && startswith(path, "/powerflow/hard-reset/")
    run_id = _webui_urldecode(path[(lastindex("/powerflow/hard-reset/") + 1):end])
    job = get_webui_powerflow_job(run_id)
    current_phase = get(job, "current_phase", nothing)
    response = handle_powerflow_hard_reset(run_id)
    _webui_log_route!(log_root, "webui_hard_reset_requested", verb, path; status = response.status == 200 ? "accepted" : "rejected", run_id, current_phase)
    if response.status == 200 && runtime !== nothing
      _webui_log_route!(log_root, "webui_shutdown_requested", verb, path; status = "accepted", run_id, current_phase)
      @async begin
        sleep(0.05)
        _webui_request_shutdown!(runtime; reason = :hard_reset)
      end
    end
    return response
  elseif verb == "GET" && startswith(path, "/powerflow/artifacts/")
    return handle_powerflow_artifacts(_webui_urldecode(path[(lastindex("/powerflow/artifacts/") + 1):end]))
  elseif verb == "GET" && startswith(path, "/powerflow/artifact-zip/")
    run_id = _webui_urldecode(path[(lastindex("/powerflow/artifact-zip/") + 1):end])
    response = handle_powerflow_artifacts_zip(run_id)
    _webui_log_route!(log_root, "artifacts_zip_downloaded", verb, path; status = response.status, run_id)
    return response
  elseif verb == "GET" && startswith(path, "/powerflow/artifact/")
    remainder = path[(lastindex("/powerflow/artifact/") + 1):end]
    segments = split(remainder, '/'; limit = 2)
    length(segments) == 2 || return _webui_html(render_webui_error(400, "Artifact route requires a run ID and artifact name."); status = 400)
    run_id, artifact_name = _webui_urldecode.(segments)
    download = get(query, "download", "") == "1"
    response = download ? handle_powerflow_artifact_download(run_id, artifact_name) : handle_powerflow_artifact(run_id, artifact_name)
    _webui_log_route!(log_root, download ? "artifact_downloaded" : "artifact_opened", verb, path; status = response.status, run_id, artifact = artifact_name)
    return response
  elseif verb == "GET" && path == "/powerflow/history"
    _webui_log_route!(log_root, "history_opened", verb, path; status = "opened")
    return handle_powerflow_history(output_root)
  elseif verb == "POST" && path == "/powerflow/refresh"
    handle_powerflow_refresh(output_root)
    _webui_log_route!(log_root, "history_refreshed", verb, path; status = "succeeded")
    return _webui_redirect("/powerflow/history")
  elseif verb == "POST" && startswith(path, "/powerflow/delete/")
    run_id = _webui_urldecode(path[(lastindex("/powerflow/delete/") + 1):end])
    response = handle_powerflow_delete(run_id, output_root)
    event = response.status == 409 ? "run_delete_rejected" : "run_deleted"
    _webui_log_route!(log_root, event, verb, path; status = response.status == 303 ? "succeeded" : "rejected", run_id)
    return response
  elseif verb == "POST" && path == "/powerflow/delete_all"
    response = handle_powerflow_delete_all(output_root)
    _webui_log_route!(log_root, "all_runs_deleted", verb, path; status = response.status == 303 ? "succeeded" : "failed")
    return response
  elseif verb == "GET" && path == "/webui/operation-log"
    _webui_log_route!(log_root, "page_opened", verb, path; status = "opened")
    return handle_webui_operation_log(log_root)
  elseif verb == "GET" && path == "/webui/last-errors"
    _webui_log_route!(log_root, "last_errors_opened", verb, path; status = "opened")
    return _webui_html(render_webui_last_errors(log_root))
  elseif verb == "POST" && path == "/webui/operation-log/clear"
    return handle_webui_operation_log_clear(log_root)
  elseif verb == "GET" && path == "/webui/operation-log/download"
    _webui_log_route!(log_root, "artifact_downloaded", verb, path; status = "succeeded", artifact = WEBUI_OPERATION_LOG_FILENAME)
    return handle_webui_operation_log(log_root; download = true)
  elseif verb == "GET" && path == "/webui/sysimage"
    # the autorefresh poll while a build runs must not fill the operation log
    get(query, "autorefresh", "") == "1" || _webui_log_route!(log_root, "sysimage_page_opened", verb, path; status = "opened")
    return handle_webui_sysimage(; output_root, message = get(query, "message", ""))
  elseif verb == "POST" && path == "/webui/sysimage/rebuild"
    return handle_webui_sysimage_rebuild(; output_root, operation_log = log_root)
  elseif verb == "POST" && path == "/webui/heartbeat"
    runtime === nothing || _webui_record_heartbeat!(runtime)
    return SparlectraWebUIResponse(204, ""; content_type = "text/plain; charset=utf-8")
  elseif verb == "POST" && path == "/webui/shutdown"
    _webui_log_route!(log_root, "webui_shutdown_requested", verb, path; status = runtime === nothing ? "unavailable" : "accepted")
    runtime === nothing && return _webui_html(render_webui_error(503, "Web UI shutdown is unavailable outside a running server."); status = 503)
    @async begin
      sleep(0.05)
      _webui_request_shutdown!(runtime; reason = :explicit_shutdown)
    end
    return _webui_html(render_webui_shutdown())
  elseif verb == "GET" && path == "/static/sparlectra.css"
    return SparlectraWebUIResponse(200, read(joinpath(@__DIR__, "static", "sparlectra.css"), String); content_type = "text/css; charset=utf-8")
  end
  return _webui_html(render_webui_error(404, "Route not found."); status = 404)
end

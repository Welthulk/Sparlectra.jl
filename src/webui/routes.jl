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
    return _webui_html(render_powerflow_form(; output_root, case_directory = runtime === nothing ? nothing : runtime.case_directory, operation_log = runtime === nothing ? webui_operation_log_path(output_root) : runtime.operation_log, selected_config_file = runtime === nothing ? "" : runtime.config_file))
  elseif verb == "POST" && path == "/powerflow/run"
    try
      result = handle_powerflow_run(form; default_output_root = output_root, case_directory = runtime === nothing ? nothing : runtime.case_directory, runner = runtime === nothing ? start_powerflow_run : runtime.runner, operation_log = log_root)
      manual_case = strip(String(something(_webui_form_value(form, "casefile_manual", ""), "")))
      requested_case = isempty(manual_case) ? String(something(_webui_form_value(form, "casefile", ""), "")) : manual_case
      if haskey(result, "run_id") && !haskey(result, "reason")
        get(form, "run_diagnostics", nothing) === nothing || _webui_log_route!(log_root, "diagnostics_enabled", verb, path; status = "enabled", run_id = result["run_id"])
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
        error_message = sprint(showerror, err),
        selected_casefile,
        selected_config_file,
      ); status = 400)
    end
  elseif verb == "GET" && startswith(path, "/powerflow/result/")
    run_id = _webui_urldecode(path[(lastindex("/powerflow/result/") + 1):end])
    _webui_log_route!(log_root, "powerflow_status_opened", verb, path; status = "opened", run_id)
    return handle_powerflow_result(run_id)
  elseif verb == "POST" && startswith(path, "/powerflow/abort/")
    run_id = _webui_urldecode(path[(lastindex("/powerflow/abort/") + 1):end])
    prior_status = get(get_webui_powerflow_job(run_id), "status", "")
    response = handle_powerflow_abort(run_id)
    event = prior_status == "aborting" ? "powerflow_abort_already_requested" : prior_status == "aborted" ? "powerflow_abort_ignored" : "powerflow_abort_requested"
    event_status = prior_status == "aborting" ? "already_aborting" : prior_status == "aborted" ? "already_aborted" : response.status == 303 ? "accepted" : "rejected"
    _webui_log_route!(log_root, event, verb, path; status = event_status, run_id)
    return response
  elseif verb == "GET" && startswith(path, "/powerflow/artifacts/")
    return handle_powerflow_artifacts(_webui_urldecode(path[(lastindex("/powerflow/artifacts/") + 1):end]))
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
  elseif verb == "GET" && path == "/webui/operation-log/download"
    _webui_log_route!(log_root, "artifact_downloaded", verb, path; status = "succeeded", artifact = WEBUI_OPERATION_LOG_FILENAME)
    return handle_webui_operation_log(log_root; download = true)
  elseif verb == "POST" && path == "/webui/heartbeat"
    runtime === nothing || _webui_record_heartbeat!(runtime)
    return SparlectraWebUIResponse(204, ""; content_type = "text/plain; charset=utf-8")
  elseif verb == "POST" && path == "/webui/shutdown"
    _webui_log_route!(log_root, "webui_shutdown_requested", verb, path; status = runtime === nothing ? "unavailable" : "accepted")
    runtime === nothing && return _webui_html(render_webui_error(503, "Web UI shutdown is unavailable outside a running server."); status = 503)
    @async begin
      sleep(0.05)
      _webui_request_shutdown!(runtime)
    end
    return _webui_html(render_webui_shutdown())
  elseif verb == "GET" && path == "/static/sparlectra.css"
    return SparlectraWebUIResponse(200, read(joinpath(@__DIR__, "static", "sparlectra.css"), String); content_type = "text/css; charset=utf-8")
  end
  return _webui_html(render_webui_error(404, "Route not found."); status = 404)
end

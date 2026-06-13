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

function handle_powerflow_result(run_id::AbstractString)::SparlectraWebUIResponse
  result = get_webui_powerflow_job(run_id)
  status = get(result, "success", false) || get(result, "reason", "") != "run_not_found" ? 200 : 404
  return _webui_html(render_powerflow_result(result); status = status)
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

function handle_powerflow_artifact(run_id::AbstractString, artifact_name::AbstractString)::SparlectraWebUIResponse
  artifact = resolve_powerflow_artifact(run_id, artifact_name)
  if artifact isa AbstractDict
    reason = get(artifact, "reason", "artifact_error")
    status = reason in ("unsafe_artifact_name", "artifact_not_found") ? 400 : 404
    return _webui_html(render_webui_error(status, get(artifact, "message", reason)); status = status)
  end
  if artifact.mime_type in _WEBUI_TEXT_MIME_TYPES
    content = read(artifact.path, String)
    page = _webui_layout("Artifact: $(artifact.name)", "<section class=\"artifact-text-page\"><p><a class=\"button\" href=\"?download=1\">Download</a></p><pre class=\"artifact-text\">$(_webui_escape(content))</pre></section>"; show_back = true, main_class = "page artifact-page")
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

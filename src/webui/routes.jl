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

function route_sparlectra_webui(method::AbstractString, target::AbstractString, form::AbstractDict = Dict{String,String}(); output_root::AbstractString = "results/powerflow_service")::SparlectraWebUIResponse
  path, query = _webui_split_target(target)
  verb = uppercase(String(method))
  if verb == "GET" && path == "/assets/logo.png"
    return handle_webui_logo()
  elseif verb == "GET" && startswith(path, "/help/")
    return handle_webui_help(_webui_urldecode(path[(lastindex("/help/") + 1):end]))
  elseif verb == "GET" && path == "/docs"
    return handle_webui_docs_index()
  elseif verb == "GET" && startswith(path, "/docs/")
    return handle_webui_doc_page(_webui_urldecode(path[(lastindex("/docs/") + 1):end]))
  elseif verb == "GET" && path in ("/", "/powerflow")
    return _webui_html(render_powerflow_form(; output_root = output_root))
  elseif verb == "POST" && path == "/powerflow/run"
    try
      result = handle_powerflow_run(form; default_output_root = output_root)
      haskey(result, "run_id") && get(result, "reason", "") != "execution_error" && return _webui_redirect("/powerflow/result/$(_webui_urlencode(result["run_id"]))")
      return _webui_html(render_powerflow_result(result); status = 400)
    catch err
      return _webui_html(render_powerflow_form(; output_root = output_root, error_message = sprint(showerror, err)); status = 400)
    end
  elseif verb == "GET" && startswith(path, "/powerflow/result/")
    return handle_powerflow_result(_webui_urldecode(path[(lastindex("/powerflow/result/") + 1):end]))
  elseif verb == "GET" && startswith(path, "/powerflow/artifacts/")
    return handle_powerflow_artifacts(_webui_urldecode(path[(lastindex("/powerflow/artifacts/") + 1):end]))
  elseif verb == "GET" && startswith(path, "/powerflow/artifact/")
    remainder = path[(lastindex("/powerflow/artifact/") + 1):end]
    segments = split(remainder, '/'; limit = 2)
    length(segments) == 2 || return _webui_html(render_webui_error(400, "Artifact route requires a run ID and artifact name."); status = 400)
    run_id, artifact_name = _webui_urldecode.(segments)
    return get(query, "download", "") == "1" ? handle_powerflow_artifact_download(run_id, artifact_name) : handle_powerflow_artifact(run_id, artifact_name)
  elseif verb == "GET" && path == "/powerflow/history"
    root = get(query, "output_root", String(output_root))
    return handle_powerflow_history(root)
  elseif verb == "POST" && path == "/powerflow/refresh"
    root = get(form, "output_root", String(output_root))
    handle_powerflow_refresh(root)
    return _webui_redirect("/powerflow/history?output_root=$(_webui_urlencode(root))")
  elseif verb == "GET" && path == "/static/sparlectra.css"
    return SparlectraWebUIResponse(200, read(joinpath(@__DIR__, "static", "sparlectra.css"), String); content_type = "text/css; charset=utf-8")
  end
  return _webui_html(render_webui_error(404, "Route not found."); status = 404)
end

function _webui_escape(value)::String
  text = string(something(value, ""))
  return replace(text, '&' => "&amp;", '<' => "&lt;", '>' => "&gt;", '"' => "&quot;", '\'' => "&#39;")
end

function _webui_urlencode(value)::String
  io = IOBuffer()
  for byte in codeunits(String(value))
    char = Char(byte)
    if isascii(char) && (isletter(char) || isdigit(char) || char in ('-', '_', '.', '~'))
      write(io, byte)
    else
      print(io, '%', uppercase(string(byte, base = 16, pad = 2)))
    end
  end
  return String(take!(io))
end

function _webui_option(value, label, selected)
  marker = String(value) == String(selected) ? " selected" : ""
  return "<option value=\"$(_webui_escape(value))\"$(marker)>$(_webui_escape(label))</option>"
end

function _webui_select(name, values, selected)
  options = join((_webui_option(value, replace(String(value), '_' => ' '), selected) for value in values), "")
  return "<select id=\"$(name)\" name=\"$(name)\">$(options)</select>"
end

function _webui_layout(title::AbstractString, content::AbstractString)::String
  return """<!doctype html>
<html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">
<title>$(_webui_escape(title)) · Sparlectra</title><link rel=\"stylesheet\" href=\"/static/sparlectra.css\"></head>
<body><header class=\"site-header\"><a class=\"brand\" href=\"/powerflow\"><img class=\"brand-logo\" src=\"/assets/logo.png\" alt=\"Sparlectra.jl logo\"><span>Sparlectra.jl</span></a><nav><a href=\"/powerflow\">New run</a><a href=\"/powerflow/history\">Run history</a></nav></header>
<main><h1>$(_webui_escape(title))</h1>$(content)</main><footer>Local PowerFlow Web UI · loopback access only</footer></body></html>"""
end

function render_powerflow_form(; output_root::AbstractString = "results/powerflow_service", error_message = nothing)::String
  error_html = error_message === nothing ? "" : "<div class=\"alert error\">$(_webui_escape(error_message))</div>"
  form = """
$(error_html)<p class=\"lede\">Run a local MATPOWER case through the Sparlectra PowerFlow service.</p>
<form method=\"post\" action=\"/powerflow/run\" class=\"panel form-grid\">
<label>MATPOWER case file<input name=\"casefile\" required placeholder=\"/path/to/case9.m\"></label>
<label>Configuration template file<input name=\"config_file\" required value=\"$(_webui_escape(DEFAULT_SPARLECTRA_CONFIG_PATH))\"></label>
<label class=\"span-2\">Output root directory<input name=\"output_root\" required value=\"$(_webui_escape(output_root))\"></label>
<label>PowerFlow tolerance<input name=\"power_flow_tol\" type=\"number\" step=\"any\" min=\"0\" value=\"1e-8\"></label>
<label>Maximum iterations<input name=\"power_flow_max_iter\" type=\"number\" min=\"1\" value=\"80\"></label>
<label class=\"check\"><input name=\"power_flow_autodamp\" type=\"checkbox\" checked> Autodamping enabled</label>
<label>Autodamping minimum<input name=\"power_flow_autodamp_min\" type=\"number\" step=\"any\" min=\"0\" max=\"1\" value=\"0.05\"></label>
<label class=\"check\"><input name=\"power_flow_qlimits_enabled\" type=\"checkbox\" checked> Q-limit handling enabled</label>
<label>Wrong-branch detection$(_webui_select("power_flow_wrong_branch_detection", WRONG_BRANCH_DETECTION_VALUES, :warn))</label>
<label>Start angle mode$(_webui_select("power_flow_start_angle_mode", POWERFLOW_START_ANGLE_MODE_VALUES, :dc))</label>
<label>Start voltage mode$(_webui_select("power_flow_start_voltage_mode", POWERFLOW_START_VOLTAGE_MODE_VALUES, :profile_blend))</label>
<label>Logfile output mode$(_webui_select("output_logfile_results", OUTPUT_LOGFILE_RESULTS_VALUES, :compact))</label>
<label class=\"check\"><input name=\"benchmark_enabled\" type=\"checkbox\"> Benchmark enabled</label>
<label>Benchmark samples<input name=\"benchmark_samples\" type=\"number\" min=\"1\" value=\"10\"></label>
<label>Benchmark seconds<input name=\"benchmark_seconds\" type=\"number\" step=\"any\" min=\"0\" value=\"1.0\"></label>
<div class=\"span-2 actions\"><button type=\"submit\">Start PowerFlow run</button></div></form>"""
  return _webui_layout("PowerFlow run", form)
end

const _WEBUI_RESULT_FIELDS = (
  "run_id", "schema_version", "status", "success", "converged", "solution_available",
  "iterations", "final_mismatch", "reason", "message", "casefile", "config_file", "output_dir",
)

function render_powerflow_result(result::AbstractDict)::String
  run_id = get(result, "run_id", "")
  rows = join(("<tr><th>$(_webui_escape(field))</th><td>$(_webui_escape(get(result, field, nothing)))</td></tr>" for field in _WEBUI_RESULT_FIELDS), "")
  links = isempty(String(run_id)) ? "" : "<div class=\"actions\"><a class=\"button\" href=\"/powerflow/artifacts/$(_webui_urlencode(run_id))\">View artifacts</a></div>"
  return _webui_layout("PowerFlow result", "<section class=\"panel\"><table class=\"details\">$(rows)</table>$(links)</section>")
end

function render_powerflow_artifacts(run_id::AbstractString, artifacts)::String
  if artifacts isa AbstractDict
    return render_powerflow_result(artifacts)
  end
  rows = join((begin
    name = get(artifact, "name", "")
    href = "/powerflow/artifact/$(_webui_urlencode(run_id))/$(_webui_urlencode(name))"
    "<tr><td><a href=\"$(href)\">$(_webui_escape(name))</a></td><td>$(_webui_escape(get(artifact, "kind", "")))</td><td>$(_webui_escape(get(artifact, "mime_type", "")))</td><td>$(_webui_escape(get(artifact, "size_bytes", "")))</td><td>$(_webui_escape(get(artifact, "description", "")))</td></tr>"
  end for artifact in artifacts), "")
  table = "<section class=\"panel\"><p>Run <code>$(_webui_escape(run_id))</code></p><table><thead><tr><th>Name</th><th>Kind</th><th>MIME type</th><th>Bytes</th><th>Description</th></tr></thead><tbody>$(rows)</tbody></table></section>"
  return _webui_layout("Artifacts", table)
end

function render_powerflow_history(runs, output_root::AbstractString)::String
  rows = join((begin
    run_id = get(run, "run_id", "")
    available = get(run, "available", false)
    link = available ? "<a href=\"/powerflow/result/$(_webui_urlencode(run_id))\">$(_webui_escape(run_id))</a>" : _webui_escape(run_id)
    fields = (link, get(run, "status", ""), get(run, "success", ""), available, get(run, "iterations", ""), get(run, "final_mismatch", ""), get(run, "casefile", ""), get(run, "config_file", ""))
    cells = "<td>$(fields[1])</td>" * join(("<td>$(_webui_escape(field))</td>" for field in fields[2:end]), "")
    "<tr>$(cells)</tr>"
  end for run in runs), "")
  content = """<form method=\"post\" action=\"/powerflow/refresh\" class=\"panel inline-form\"><label>Output root<input name=\"output_root\" value=\"$(_webui_escape(output_root))\"></label><button type=\"submit\">Refresh registry</button></form>
<section class=\"panel\"><table><thead><tr><th>Run ID</th><th>Status</th><th>Success</th><th>Available</th><th>Iterations</th><th>Final mismatch</th><th>Case file</th><th>Config file</th></tr></thead><tbody>$(rows)</tbody></table></section>"""
  return _webui_layout("Run history", content)
end

function render_webui_error(status::Integer, message::AbstractString)::String
  return _webui_layout("Request error", "<div class=\"alert error\"><strong>$(_webui_escape(status))</strong> $(_webui_escape(message))</div>")
end

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

"""Render a file-path dropdown while showing only each file's basename."""
function _webui_path_select(name::AbstractString, paths, selected::AbstractString)::String
  selected_path = isempty(selected) ? (isempty(paths) ? "" : first(paths)) : selected
  option_paths = isempty(selected_path) || selected_path in paths ? paths : [selected_path; paths]
  options = join((_webui_option(path, basename(path), selected_path) for path in option_paths), "")
  return "<select id=\"$(name)\" name=\"$(name)\" required>$(options)</select>"
end

function _webui_layout(title::AbstractString, content::AbstractString; show_back::Bool = false)::String
  back_button = show_back ? "<div class=\"page-toolbar\"><a class=\"button back-button\" href=\"/powerflow\" onclick=\"if (document.referrer.startsWith(location.origin)) { history.back(); return false; }\" aria-label=\"Go back to the previous page\">← Back</a></div>" : ""
  return """<!doctype html>
<html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">
<title>$(_webui_escape(title)) · Sparlectra</title><link rel=\"stylesheet\" href=\"/static/sparlectra.css\"></head>
<body><header class=\"site-header\"><a class=\"brand\" href=\"/powerflow\"><img class=\"brand-logo\" src=\"/assets/logo.png\" alt=\"Sparlectra.jl logo\"><span>Sparlectra.jl</span></a><nav><a href=\"/powerflow\">New run</a><a href=\"/powerflow/history\">Run history</a><a href=\"/docs\">Documentation</a><a class=\"project-docs-link\" href=\"https://welthulk.github.io/Sparlectra.jl/\" target=\"_blank\" rel=\"noopener noreferrer\"><svg class=\"github-icon\" viewBox=\"0 0 16 16\" aria-hidden=\"true\"><path fill=\"currentColor\" d=\"M8 0C3.58 0 0 3.64 0 8.13c0 3.59 2.29 6.64 5.47 7.72.4.08.55-.18.55-.39 0-.19-.01-.83-.01-1.51-2.01.38-2.53-.5-2.69-.96-.09-.23-.48-.96-.82-1.15-.28-.15-.68-.53-.01-.54.63-.01 1.08.59 1.23.83.72 1.23 1.87.88 2.33.67.07-.53.28-.88.51-1.08-1.78-.21-3.64-.91-3.64-4.02 0-.89.31-1.62.82-2.19-.08-.2-.36-1.04.08-2.16 0 0 .67-.22 2.2.84A7.4 7.4 0 0 1 8 3.93c.68 0 1.36.09 2 .27 1.53-1.06 2.2-.84 2.2-.84.44 1.12.16 1.96.08 2.16.51.57.82 1.3.82 2.19 0 3.12-1.87 3.81-3.65 4.02.29.25.54.74.54 1.5 0 1.08-.01 1.95-.01 2.22 0 .22.15.47.55.39A8.15 8.15 0 0 0 16 8.13C16 3.64 12.42 0 8 0Z\"/></svg><span>Project Docs</span></a></nav></header>
<main>$(back_button)<h1>$(_webui_escape(title))</h1>$(content)</main><footer>Local PowerFlow Web UI · loopback access only</footer></body></html>"""
end

function _webui_help_link(topic::AbstractString, label::AbstractString)::String
  return "<a class=\"help-link\" href=\"/help/$(_webui_urlencode(topic))\" aria-label=\"Help for $(_webui_escape(label))\" title=\"Help for $(_webui_escape(label))\">?</a>"
end

function _webui_field_label(field::AbstractString, label::AbstractString)::String
  topic = WEBUI_FORM_HELP_TOPICS[String(field)]
  return "<span class=\"field-label\">$(_webui_escape(label)) $(_webui_help_link(topic, label))</span>"
end

function render_powerflow_form(; output_root::AbstractString = "results/powerflow_service", error_message = nothing, application_root::AbstractString = _webui_application_root(), selected_casefile::AbstractString = "", selected_config_file::AbstractString = "")::String
  error_html = error_message === nothing ? "" : "<div class=\"alert error\">$(_webui_escape(error_message))</div>"
  casefiles = _webui_casefile_options(application_root)
  config_files = _webui_config_file_options(application_root)
  case_directory = joinpath(application_root, "data", "mpower")
  config_directory = joinpath(application_root, "examples")
  case_placeholder = joinpath(case_directory, "case.m")
  case_control = isempty(casefiles) ? "<input name=\"casefile\" required value=\"$(_webui_escape(selected_casefile))\" placeholder=\"$(_webui_escape(case_placeholder))\">" : _webui_path_select("casefile", casefiles, selected_casefile)
  config_default = isempty(selected_config_file) ? DEFAULT_SPARLECTRA_CONFIG_PATH : selected_config_file
  config_control = isempty(config_files) ? "<input name=\"config_file\" required value=\"$(_webui_escape(config_default))\">" : _webui_path_select("config_file", config_files, selected_config_file)
  form = """
$(error_html)<p class=\"lede\">Run a local MATPOWER case through the Sparlectra PowerFlow service.</p>
<form id=\"powerflow-run-form\" method=\"post\" action=\"/powerflow/run\" class=\"panel form-grid\" onsubmit=\"this.classList.add('is-submitting'); this.setAttribute('aria-busy', 'true'); this.querySelector('button[type=submit]').disabled = true;\">
<label>$(_webui_field_label("casefile", "MATPOWER case file"))$(case_control)<small class="field-hint">Cases from <code>$(_webui_escape(case_directory))</code></small></label>
<label>$(_webui_field_label("config_file", "Configuration template file"))$(config_control)<small class="field-hint">Configurations from <code>$(_webui_escape(config_directory))</code></small></label>
<label class=\"span-2\">$(_webui_field_label("output_root", "Output root directory"))<input name=\"output_root\" required value=\"$(_webui_escape(output_root))\"></label>
<label>$(_webui_field_label("power_flow_tol", "PowerFlow tolerance"))<input name=\"power_flow_tol\" type=\"number\" step=\"any\" min=\"0\" value=\"1e-8\"></label>
<label>$(_webui_field_label("power_flow_max_iter", "Maximum iterations"))<input name=\"power_flow_max_iter\" type=\"number\" min=\"1\" value=\"80\"></label>
<label class=\"check\"><input name=\"power_flow_autodamp\" type=\"checkbox\" checked>$(_webui_field_label("power_flow_autodamp", "Autodamping enabled"))</label>
<label>$(_webui_field_label("power_flow_autodamp_min", "Autodamping minimum"))<input name=\"power_flow_autodamp_min\" type=\"number\" step=\"any\" min=\"0\" max=\"1\" value=\"0.05\"></label>
<label class=\"check\"><input name=\"power_flow_qlimits_enabled\" type=\"checkbox\" checked>$(_webui_field_label("power_flow_qlimits_enabled", "Q-limit handling enabled"))</label>
<label>$(_webui_field_label("power_flow_wrong_branch_detection", "Wrong-branch detection"))$(_webui_select("power_flow_wrong_branch_detection", WRONG_BRANCH_DETECTION_VALUES, :warn))</label>
<label>$(_webui_field_label("power_flow_start_angle_mode", "Start angle mode"))$(_webui_select("power_flow_start_angle_mode", POWERFLOW_START_ANGLE_MODE_VALUES, :dc))</label>
<label>$(_webui_field_label("power_flow_start_voltage_mode", "Start voltage mode"))$(_webui_select("power_flow_start_voltage_mode", POWERFLOW_START_VOLTAGE_MODE_VALUES, :profile_blend))</label>
<label>$(_webui_field_label("output_logfile_results", "Logfile output mode"))$(_webui_select("output_logfile_results", OUTPUT_LOGFILE_RESULTS_VALUES, :compact))</label>
<label class=\"check\"><input name=\"benchmark_enabled\" type=\"checkbox\">$(_webui_field_label("benchmark_enabled", "Benchmark enabled"))</label>
<label>$(_webui_field_label("benchmark_samples", "Benchmark samples"))<input name=\"benchmark_samples\" type=\"number\" min=\"1\" value=\"10\"></label>
<label>$(_webui_field_label("benchmark_seconds", "Benchmark seconds"))<input name=\"benchmark_seconds\" type=\"number\" step=\"any\" min=\"0\" value=\"1.0\"></label>
<div class=\"span-2 actions\"><button class=\"powerflow-submit\" type=\"submit\"><span class=\"submit-spinner\" aria-hidden=\"true\"></span><span class=\"submit-label\">Start PowerFlow run</span><span class=\"submit-progress-label\" role=\"status\" aria-live=\"polite\">Running PowerFlow…</span></button></div></form>
<script>
window.addEventListener('pageshow', function () {
  const form = document.getElementById('powerflow-run-form');
  if (form === null) return;
  form.classList.remove('is-submitting');
  form.removeAttribute('aria-busy');
  const submitButton = form.querySelector('button[type=submit]');
  if (submitButton !== null) submitButton.disabled = false;
});
</script>"""
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
  return _webui_layout("PowerFlow result", "<section class=\"panel\"><table class=\"details\">$(rows)</table>$(links)</section>"; show_back = true)
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
  return _webui_layout("Artifacts", table; show_back = true)
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
  return _webui_layout("Run history", content; show_back = true)
end

function render_webui_error(status::Integer, message::AbstractString)::String
  return _webui_layout("Request error", "<div class=\"alert error\"><strong>$(_webui_escape(status))</strong> $(_webui_escape(message))</div>"; show_back = true)
end

function render_webui_help(metadata, excerpt::AbstractString)::String
  markdown_html = render_webui_markdown(excerpt; current_page = metadata.page)
  page_url = _webui_urlencode(metadata.page)
  source_file = WEBUI_DOC_PAGES[metadata.page].file
  content = "<section class=\"panel help-page help-panel\">$(markdown_html)<p><a href=\"/docs/$(page_url)\">View the full documentation context</a></p><p class=\"source-note\">Source: <code>$(_webui_escape(source_file))</code></p></section>"
  return _webui_layout(metadata.label, content; show_back = true)
end

function render_webui_docs_index(pages::AbstractDict)::String
  links = join(("<li><a href=\"/docs/$(_webui_urlencode(page))\">$(_webui_escape(metadata.title))</a></li>" for (page, metadata) in sort!(collect(pages); by = first)), "")
  content = "<section class=\"panel docs-page docs-content\"><p>Selected repository documentation pages are rendered directly from their Markdown sources.</p><ul class=\"docs-index\">$(links)</ul></section>"
  return _webui_layout("Documentation", content; show_back = true)
end

function render_webui_doc_page(page::AbstractString, metadata, markdown_text::AbstractString)::String
  content = "<section class=\"panel docs-page docs-content\">$(render_webui_markdown(markdown_text; current_page = page))</section>"
  return _webui_layout(metadata.title, content; show_back = true)
end

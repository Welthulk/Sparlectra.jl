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

const WEBUI_STATUS_AUTO_REFRESH_SECONDS = 2
const _WEBUI_ACTIVE_RUN_STATUSES = Set(("queued", "running", "aborting"))

function _format_elapsed_duration(seconds)::String
  seconds === nothing && return "—"
  elapsed = try
    Float64(seconds)
  catch
    tryparse(Float64, string(seconds))
  end
  (elapsed === nothing || !isfinite(elapsed)) && return "—"
  total_seconds = max(0, floor(Int, elapsed))
  hours, remainder = divrem(total_seconds, 3600)
  minutes, secs = divrem(remainder, 60)
  return lpad(hours, 2, '0') * ":" * lpad(minutes, 2, '0') * ":" * lpad(secs, 2, '0')
end

function _webui_elapsed_seconds(result::AbstractDict, active::Bool)
  elapsed = get(result, "elapsed_seconds", nothing)
  elapsed === nothing || return elapsed
  active || return nothing
  started_at = get(result, "started_at", get(result, "submitted_at", nothing))
  started_at === nothing && return nothing
  started = if started_at isa Dates.DateTime
    started_at
  else
    try
      Dates.DateTime(replace(string(started_at), r"Z$" => ""))
    catch
      return nothing
    end
  end
  return max(0.0, Dates.value(Dates.now(Dates.UTC) - started) / 1000)
end

function _webui_layout(title::AbstractString, content::AbstractString; show_back::Bool = false, main_class::AbstractString = "page", refresh_url = nothing, refresh_seconds::Integer = WEBUI_STATUS_AUTO_REFRESH_SECONDS)::String
  back_button = show_back ? "<div class=\"page-toolbar\"><a class=\"button back-button\" href=\"/powerflow\" onclick=\"if (document.referrer.startsWith(location.origin)) { history.back(); return false; }\" aria-label=\"Go back to the previous page\">← Back</a></div>" : ""
  version_text = "Sparlectra.jl v$(version())"
  refresh_meta = refresh_url === nothing ? "" : "<meta http-equiv=\"refresh\" content=\"$(refresh_seconds); url=$(_webui_escape(refresh_url))\">"
  return """<!doctype html>
<html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">
$(refresh_meta)<title>$(_webui_escape(title)) · Sparlectra</title><link rel=\"stylesheet\" href=\"/static/sparlectra.css\"></head>
<body><header class=\"site-header\"><a class=\"brand\" href=\"/powerflow\"><img class=\"brand-logo\" src=\"/assets/logo.png\" alt=\"Sparlectra.jl logo\"><span>$(_webui_escape(version_text))</span></a><nav><a href=\"/powerflow\">New run</a><a href=\"/powerflow/history\">Run history</a><a href=\"/webui/operation-log\">Operation Log</a><a href=\"/docs\">Documentation</a><a class=\"project-docs-link\" href=\"https://welthulk.github.io/Sparlectra.jl/\" target=\"_blank\" rel=\"noopener noreferrer\"><svg class=\"github-icon\" viewBox=\"0 0 16 16\" aria-hidden=\"true\"><path fill=\"currentColor\" d=\"M8 0C3.58 0 0 3.64 0 8.13c0 3.59 2.29 6.64 5.47 7.72.4.08.55-.18.55-.39 0-.19-.01-.83-.01-1.51-2.01.38-2.53-.5-2.69-.96-.09-.23-.48-.96-.82-1.15-.28-.15-.68-.53-.01-.54.63-.01 1.08.59 1.23.83.72 1.23 1.87.88 2.33.67.07-.53.28-.88.51-1.08-1.78-.21-3.64-.91-3.64-4.02 0-.89.31-1.62.82-2.19-.08-.2-.36-1.04.08-2.16 0 0 .67-.22 2.2.84A7.4 7.4 0 0 1 8 3.93c.68 0 1.36.09 2 .27 1.53-1.06 2.2-.84 2.2-.84.44 1.12.16 1.96.08 2.16.51.57.82 1.3.82 2.19 0 3.12-1.87 3.81-3.65 4.02.29.25.54.74.54 1.5 0 1.08-.01 1.95-.01 2.22 0 .22.15.47.55.39A8.15 8.15 0 0 0 16 8.13C16 3.64 12.42 0 8 0Z\"/></svg><span>Project Docs</span></a><form method="post" action="/webui/shutdown" class="exit-form"><button type="submit" class="exit-button">Stop Web UI</button></form></nav></header>
<main class="$(main_class)">$(back_button)<h1>$(_webui_escape(title))</h1>$(content)</main><footer>$(_webui_escape(version_text)) · Local PowerFlow Web UI · loopback access only</footer>
<script>
(function () {
  const sendHeartbeat = function () {
    fetch('/webui/heartbeat', {method: 'POST', keepalive: true}).catch(function () {});
  };
  sendHeartbeat();
  window.setInterval(sendHeartbeat, 5000);
})();
</script></body></html>"""
end

function _webui_help_link(topic::AbstractString, label::AbstractString)::String
  return "<a class=\"help-link\" href=\"/help/$(_webui_urlencode(topic))\" aria-label=\"Help for $(_webui_escape(label))\" title=\"Help for $(_webui_escape(label))\">?</a>"
end

function _webui_field_label(field::AbstractString, label::AbstractString)::String
  topic = WEBUI_FORM_HELP_TOPICS[String(field)]
  return "<span class=\"field-label\">$(_webui_escape(label)) $(_webui_help_link(topic, label))</span>"
end

function _webui_active_run_banner(active_run)::String
  active_run === nothing && return ""
  run_id = string(get(active_run, "run_id", ""))
  status = lowercase(string(get(active_run, "status", "running")))
  status in ("queued", "running", "aborting") || return ""
  explanation = status == "aborting" ? "<p>Abort requested. A new run may start once this run reaches aborted status.</p>" : ""
  abort_form = status == "aborting" ? "" : "<form method=\"post\" action=\"/powerflow/abort/$(_webui_urlencode(run_id))\"><button type=\"submit\" class=\"danger-button\">Abort</button></form>"
  return "<section class=\"panel active-run-banner\"><strong>PowerFlow run is $(status):</strong> <code>$(_webui_escape(run_id))</code>$(explanation)<div class=\"actions\"><a class=\"button\" href=\"/powerflow/result/$(_webui_urlencode(run_id))\">Open status</a>$(abort_form)</div></section>"
end

function render_powerflow_form(; output_root::AbstractString = "results/powerflow_service", case_directory::Union{Nothing,AbstractString} = nothing, operation_log::AbstractString = webui_operation_log_path(output_root), error_message = nothing, application_root::AbstractString = _webui_application_root(), selected_casefile::AbstractString = "", selected_config_file::AbstractString = "", active_run = get_active_webui_powerflow_job())::String
  error_html = error_message === nothing ? "" : "<div class=\"alert error\">$(_webui_escape(error_message))</div>"
  casefiles = case_directory === nothing ? _webui_casefile_options(application_root) : _webui_casefile_options_in_directory(case_directory)
  bundled_case_directory = joinpath(application_root, "data", "mpower")
  effective_case_directory = case_directory === nothing ? bundled_case_directory : String(case_directory)
  selected_value = strip(selected_casefile)
  existing_value = selected_value in casefiles ? selected_value : ""
  manual_value = isempty(existing_value) ? selected_value : ""
  case_options = join(("<option value=\"$(_webui_escape(casefile))\"$(casefile == existing_value ? " selected" : "")>$(_webui_escape(casefile))</option>" for casefile in casefiles), "")
  case_select = "<select id=\"casefile\" name=\"casefile\"><option value=\"\">-- choose existing case --</option>$(case_options)</select>"
  case_manual = "<input id=\"casefile_manual\" name=\"casefile_manual\" value=\"$(_webui_escape(manual_value))\" placeholder=\"case14.m\">"
  config_default = isempty(selected_config_file) ? DEFAULT_SPARLECTRA_CONFIG_PATH : selected_config_file
  config_control = "<input type=\"hidden\" name=\"config_file\" value=\"$(_webui_escape(config_default))\"><code>$(_webui_escape(config_default))</code>"
  form = """
$(error_html)$(_webui_active_run_banner(active_run))<p class=\"lede\">Run a local MATPOWER case through the Sparlectra PowerFlow service.</p>
<form id=\"powerflow-run-form\" method=\"post\" action=\"/powerflow/run\" class=\"panel form-grid\" onsubmit=\"this.classList.add('is-submitting'); this.setAttribute('aria-busy', 'true'); this.querySelector('button[type=submit]').disabled = true;\">
<label>$(_webui_field_label("casefile", "Existing MATPOWER case"))$(case_select)<small class="field-hint">Cases from <code>$(_webui_escape(effective_case_directory))</code></small></label>
<label><span class="field-label">Or type/download MATPOWER case</span>$(case_manual)<small class="field-hint">Manual input overrides the existing-case selection.</small></label>
<details class=\"citation-box span-2\">
<summary>MATPOWER citation</summary>
<p>Sparlectra supports the <a href=\"https://matpower.org\">MATPOWER</a> case format and MATPOWER test-case data. If you use MATPOWER cases or data in publications, please follow the <a href=\"https://matpower.org/citing/\">MATPOWER citation guidance</a> and cite the software and:</p>
<p>R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas, “MATPOWER: Steady-State Operations, Planning and Analysis Tools for Power Systems Research and Education,” <em>IEEE Transactions on Power Systems</em>, 26(1), 12–19, 2011. <a href=\"https://doi.org/10.1109/TPWRS.2010.2051168\">DOI</a></p>
<p>Some MATPOWER case files, such as ACTIVSg, PEGASE, and RTE cases, may request additional case-specific citations in their file headers.</p>
</details>
<label>$(_webui_field_label("config_file", "Configuration file"))$(config_control)<small class="field-hint">Set at Web UI startup.</small></label>
<p class=\"span-2 output-root-display\"><strong>Output root:</strong> <code>$(_webui_escape(output_root))</code></p>
<p class=\"span-2 output-root-display\"><strong>Config file:</strong> <code>$(_webui_escape(config_default))</code></p>
<p class=\"span-2 output-root-display\"><strong>MATPOWER case cache:</strong> <code>$(_webui_escape(effective_case_directory))</code></p>
<p class=\"span-2 output-root-display\"><strong>Operation log:</strong> <code>$(_webui_escape(operation_log))</code></p>
<label>$(_webui_field_label("power_flow_tol", "PowerFlow tolerance"))<input name=\"power_flow_tol\" type=\"number\" step=\"any\" min=\"0\" value=\"1e-8\"></label>
<label>$(_webui_field_label("power_flow_max_iter", "Maximum iterations"))<input name=\"power_flow_max_iter\" type=\"number\" min=\"1\" value=\"80\"></label>
<label class=\"check\"><input name=\"power_flow_autodamp\" type=\"checkbox\" checked>$(_webui_field_label("power_flow_autodamp", "Autodamping enabled"))</label>
<label>$(_webui_field_label("power_flow_autodamp_min", "Autodamping minimum"))<input name=\"power_flow_autodamp_min\" type=\"number\" step=\"any\" min=\"0\" max=\"1\" value=\"0.05\"></label>
<label class=\"check\"><input name=\"power_flow_qlimits_enabled\" type=\"checkbox\" checked>$(_webui_field_label("power_flow_qlimits_enabled", "Q-limit handling enabled"))</label>
<label>$(_webui_field_label("power_flow_wrong_branch_detection", "Wrong-branch detection"))$(_webui_select("power_flow_wrong_branch_detection", WRONG_BRANCH_DETECTION_VALUES, :warn))</label>
<label>$(_webui_field_label("power_flow_start_angle_mode", "Start angle mode"))$(_webui_select("power_flow_start_angle_mode", POWERFLOW_START_ANGLE_MODE_VALUES, :dc))</label>
<label>$(_webui_field_label("power_flow_start_voltage_mode", "Start voltage mode"))$(_webui_select("power_flow_start_voltage_mode", POWERFLOW_START_VOLTAGE_MODE_VALUES, :profile_blend))</label>
<label>$(_webui_field_label("output_logfile_results", "Logfile output mode"))$(_webui_select("output_logfile_results", OUTPUT_LOGFILE_RESULTS_VALUES, :compact))</label>
<label>$(_webui_field_label("performance_timing", "Performance timing"))$(_webui_select("performance_timing", _WEBUI_PERFORMANCE_TIMING_VALUES, :compact))</label>
<label class=\"check\"><input name=\"run_diagnostics\" type=\"checkbox\">$(_webui_field_label("run_diagnostics", "Run diagnostics"))</label>
<fieldset class=\"span-2 detailed-csv-options\"><legend><label class=\"check\"><input name=\"detailed_result_csv\" type=\"checkbox\">$(_webui_field_label("detailed_result_csv", "Export detailed result CSV files"))</label></legend>
<label class=\"detailed-csv-format\">$(_webui_field_label("detailed_result_csv_format", "CSV format"))$(_webui_select("detailed_result_csv_format", ("technical", "excel_de", "excel_us"), "technical"))<small class=\"field-hint\">technical: comma/decimal point/no grouping; excel_de: semicolon/decimal comma/thousands dot; excel_us: comma/decimal point/thousands comma.</small></label>
</fieldset>
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
  "iterations", "final_mismatch", "reason", "message", "casefile", "resolved_casefile",
  "config_file", "started_at", "elapsed_seconds", "output_dir",
  "current_phase", "phase_started_at", "last_progress_at", "abort_requested_at",
)

function render_powerflow_result(result::AbstractDict)::String
  run_id = get(result, "run_id", "")
  rows = join(("<tr><th>$(_webui_escape(field))</th><td>$(_webui_escape(get(result, field, nothing)))</td></tr>" for field in _WEBUI_RESULT_FIELDS), "")
  status = lowercase(string(get(result, "status", "unknown")))
  active = status in _WEBUI_ACTIVE_RUN_STATUSES
  status_badge = "<span class=\"status-badge $(webui_status_class(result))\">$(_webui_escape(status))</span>"
  elapsed_duration = _format_elapsed_duration(_webui_elapsed_seconds(result, active))
  result_summary = "<div class=\"result-summary\"><div><span class=\"summary-label\">Status</span>$(status_badge)</div><div class=\"runtime-card\"><span class=\"runtime-icon\" aria-hidden=\"true\">⏱</span><span><span class=\"summary-label\">Elapsed time</span><strong>$(_webui_escape(elapsed_duration))</strong></span></div></div>"
  abort_form = status in ("queued", "running") ? "<form method=\"post\" action=\"/powerflow/abort/$(_webui_urlencode(run_id))\"><button type=\"submit\" class=\"danger-button\">Abort run</button></form>" : ""
  active_hint = active ? "<p class=\"status-refresh-hint\">This page refreshes automatically while the run is active.</p>" : ""
  abort_hint = status == "aborting" ? "<p>Waiting for the next safe cancellation point.</p>" * (get(result, "current_phase", "") == "linear_solve" ? "<p>Waiting for current sparse solve to finish.</p>" : "") : ""
  hard_reset = status == "aborting" && get(result, "hard_reset_available", false) ? "<div class=\"alert warning\"><p>Abort is still pending. The calculation may be in a non-interruptible numerical call.</p><form method=\"post\" action=\"/powerflow/hard-reset/$(_webui_urlencode(run_id))\"><button type=\"submit\" class=\"danger-button\">Hard reset Web UI</button></form></div>" : ""
  interrupted = status == "aborted_unknown" ? "<p>The Web UI was restarted or reset while this run was active. This result is not a valid solved PowerFlow result.</p>" : ""
  links = isempty(String(run_id)) ? "" : "$(abort_hint)$(hard_reset)$(interrupted)<div class=\"actions\">$(abort_form)<a class=\"button\" href=\"/powerflow/artifacts/$(_webui_urlencode(run_id))\">View artifacts</a><a class=\"button\" href=\"/powerflow/result/$(_webui_urlencode(run_id))\">Refresh status</a></div>"
  refresh_url = active && !isempty(String(run_id)) ? "/powerflow/result/$(_webui_urlencode(run_id))?autorefresh=1" : nothing
  return _webui_layout("PowerFlow result", "<section class=\"panel\">$(result_summary)<table class=\"details\">$(rows)</table>$(active_hint)$(links)</section>"; show_back = true, refresh_url)
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

function webui_status_class(run::AbstractDict)::String
  status = lowercase(string(get(run, "status", "unknown")))
  success = get(run, "success", nothing)
  status in ("running", "pending", "queued", "aborting") && return "status-running"
  status in ("aborted", "aborted_unknown", "interrupted", "cancelled", "canceled") && return "status-aborted"
  status in ("warning", "partial", "questionable") && return "status-warning"
  status in ("failed", "failure", "error", "not_converged") && return "status-error"
  status in ("succeeded", "success", "converged", "ok") && return success === false ? "status-error" : "status-success"
  success === true && return "status-success"
  success === false && return "status-error"
  return "status-unknown"
end

function _webui_run_timestamp(run::AbstractDict)::String
  timestamp = get(run, "timestamp", "Unknown")
  return isempty(strip(string(timestamp))) ? "Unknown" : string(timestamp)
end

function render_powerflow_history(runs, output_root::AbstractString; active_run = nothing)::String
  ordered_runs = sort!(collect(runs); by = run -> _webui_run_timestamp(run), rev = true)
  rows = join((begin
    run_id = string(get(run, "run_id", ""))
    available = get(run, "available", false)
    link = available ? "<a href=\"/powerflow/result/$(_webui_urlencode(run_id))\">$(_webui_escape(run_id))</a>" : _webui_escape(run_id)
    status = string(get(run, "status", "unknown"))
    status_badge = "<span class=\"status-badge $(webui_status_class(run))\">$(_webui_escape(status))</span>"
    delete_form = "<form method=\"post\" action=\"/powerflow/delete/$(_webui_urlencode(run_id))\" class=\"delete-run-form\"><button type=\"submit\" class=\"danger-button\">Delete</button></form>"
    abort_form = lowercase(status) in ("queued", "running") ? "<form method=\"post\" action=\"/powerflow/abort/$(_webui_urlencode(run_id))\"><button type=\"submit\" class=\"danger-button\">Abort</button></form>" : ""
    fields = (_webui_run_timestamp(run), link, status_badge, available, get(run, "iterations", ""), get(run, "final_mismatch", ""), get(run, "casefile", ""), get(run, "config_file", ""))
    cells = "<td>$(_webui_escape(fields[1]))</td><td>$(fields[2])</td><td>$(fields[3])</td>" * join(("<td>$(_webui_escape(field))</td>" for field in fields[4:end]), "")
    "<tr>$(cells)<td>$(abort_form)$(delete_form)</td></tr>"
  end for run in ordered_runs), "")
  content = "$(_webui_active_run_banner(active_run))<section class=\"panel history-actions\"><p><strong>Output root:</strong> <code>$(_webui_escape(output_root))</code></p><div class=\"actions\"><form method=\"post\" action=\"/powerflow/refresh\"><button type=\"submit\">Refresh registry</button></form><form method=\"post\" action=\"/powerflow/delete_all\"><button type=\"submit\" class=\"danger-button\">Delete all runs</button></form></div></section>\n<section class=\"panel\"><table><thead><tr><th>Date/Time</th><th>Run ID</th><th>Status</th><th>Available</th><th>Iterations</th><th>Final mismatch</th><th>Case file</th><th>Config file</th><th>Delete</th></tr></thead><tbody>$(rows)</tbody></table></section>"
  return _webui_layout("Run history", content; show_back = true)
end

function render_webui_operation_log(content::AbstractString)::String
  controls = "<p><a class=\"button\" href=\"/webui/operation-log/download\">Download operation log</a></p>"
  panel = "<section class=\"artifact-text-page\">$(controls)<pre class=\"artifact-text\">$(_webui_escape(content))</pre></section>"
  return _webui_layout("Operation Log", panel; show_back = true, main_class = "page artifact-page")
end

function render_webui_shutdown()::String
  return "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>Web UI stopped · Sparlectra</title><link rel=\"stylesheet\" href=\"/static/sparlectra.css\"></head><body><main><section class=\"panel shutdown-message\"><h1>Web UI stopped</h1><p>The local Sparlectra Web UI server is shutting down. You may close this window.</p></section></main></body></html>"
end

function render_webui_hard_reset()::String
  return "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>Hard reset requested · Sparlectra</title></head><body><main><section><h1>Hard reset requested</h1><p>The Web UI is shutting down to stop the stuck calculation.</p><p>Restart with <code>julia --project=. start_webui.jl</code>.</p></section></main></body></html>"
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

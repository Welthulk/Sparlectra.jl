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

# file: src/webui/views.jl
# purpose: HTML rendering for all Web UI pages: run form, results, artifacts,
#          history, docs and help, plus escaping and form-control helpers
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
  marker = _webui_form_string(value) == _webui_form_string(selected) ? " selected" : ""
  return "<option value=\"$(_webui_escape(value))\"$(marker)>$(_webui_escape(label))</option>"
end

function _webui_case_import_message(imported::Vector{String}, rejected)::String
  # One line per file with a leading status marker, so a CGMES pre-analysis
  # ("ready to compute" vs "boundary set missing") is readable at a glance
  # instead of hidden inside one long comma-separated sentence.
  lines = String[]
  for entry in imported
    marker = occursin("✅", entry) || occursin("⚠️", entry) || occursin("❌", entry) || occursin("🔗", entry) ? "" : "✔ "
    push!(lines, marker * entry)
  end
  for item in rejected
    push!(lines, "❌ $(first(item)) — $(last(item))")
  end
  isempty(lines) && push!(lines, "No files were selected for import.")
  header = isempty(imported) ? "" : "Imported $(length(imported)) file$(length(imported) == 1 ? "" : "s")"
  if !isempty(rejected)
    header = isempty(header) ? "Rejected $(length(rejected)) file$(length(rejected) == 1 ? "" : "s")" : header * ", rejected $(length(rejected))"
  end
  return isempty(header) ? join(lines, " · ") : header * ": " * join(lines, " · ")
end

"""
    _webui_qlimit_mode_selection(profile_values) -> String

Which entry the Q-limit enforcement-mode control shows. The pair of controls
(a checkbox for `power_flow.qlimits.enabled`, a select for the mode) reads as
one setting to a user, and a mode picked while the handling is off does
nothing: the select therefore shows "off" in that state, and picking "off"
turns the handling off (see the form override mapping).
"""
function _webui_qlimit_mode_selection(profile_values)
  enabled = _webui_form_bool(get(profile_values, "power_flow_qlimits_enabled", _webui_option_default("power_flow_qlimits_enabled")))
  enabled || return "off"
  return _webui_selected(profile_values, "power_flow_qlimits_enforcement_mode", _webui_option_default("power_flow_qlimits_enforcement_mode"))
end

function _webui_select(name, values, selected, extra_attrs::AbstractString = "")
  options = join((_webui_option(value, replace(_webui_form_string(value), '_' => ' '), selected) for value in values), "")
  attrs = isempty(extra_attrs) ? "" : " $(extra_attrs)"
  return "<select id=\"$(name)\" name=\"$(name)\"$(attrs)>$(options)</select>"
end

# A normal form POST gives the browser nothing to show while the server
# works, so a slow action (generating a measurement set solves the case; the
# topology test runs several estimations) is indistinguishable from a page
# that is simply free. Every form carrying `data-busy` gets what the run form
# has: a spinner in the button, the busy text in place of the label, and its
# submit buttons disabled against a second click. Wired into the layout, so a
# new slow form only needs the attribute.
const _WEBUI_BUSY_FORM_SCRIPT = """<script>
(function () {
  var decorate = function (form) {
    var busy = form.getAttribute('data-busy') || 'Working…';
    form.querySelectorAll('button[type=submit]').forEach(function (b) {
      if (b.querySelector('.submit-label') !== null) return;
      var text = b.textContent;
      b.textContent = '';
      var spinner = document.createElement('span');
      spinner.className = 'submit-spinner';
      spinner.setAttribute('aria-hidden', 'true');
      var label = document.createElement('span');
      label.className = 'submit-label';
      label.textContent = text;
      var progress = document.createElement('span');
      progress.className = 'submit-progress-label';
      progress.setAttribute('role', 'status');
      progress.setAttribute('aria-live', 'polite');
      progress.textContent = busy;
      b.appendChild(spinner);
      b.appendChild(label);
      b.appendChild(progress);
    });
    form.addEventListener('submit', function (event) {
      // The submitter's name/value into a hidden input FIRST: a disabled
      // button is dropped from the serialized form data per the HTML spec,
      // which once turned a "Diagnose" click into a normal run.
      var submitter = event.submitter;
      if (submitter && submitter.name) {
        var hidden = form.querySelector('input[type=hidden][data-submitter-value]');
        if (hidden === null) {
          hidden = document.createElement('input');
          hidden.type = 'hidden';
          hidden.setAttribute('data-submitter-value', '');
          form.appendChild(hidden);
        }
        hidden.name = submitter.name;
        hidden.value = submitter.value;
      }
      form.classList.add('is-submitting');
      form.setAttribute('aria-busy', 'true');
      form.querySelectorAll('button[type=submit]').forEach(function (b) { b.disabled = true; });
    });
  };
  var init = function () { document.querySelectorAll('form[data-busy]').forEach(decorate); };
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
  // Back button and bfcache restore a page with the class still set; without
  // this the buttons would stay dead.
  window.addEventListener('pageshow', function () {
    document.querySelectorAll('form[data-busy]').forEach(function (form) {
      form.classList.remove('is-submitting');
      form.removeAttribute('aria-busy');
      form.querySelectorAll('button[type=submit]').forEach(function (b) { b.disabled = false; });
      var stale = form.querySelector('input[type=hidden][data-submitter-value]');
      if (stale !== null) stale.remove();
    });
  });
})();
</script>"""

# A native <input type="file"> renders its button and its "no file selected"
# text in the BROWSER's language, so on a German browser an otherwise English
# page said "Durchsuchen" and "Keine Datei ausgewählt" (maintainer, reported
# 2026-09-09). The native input stays in the form (it is what carries the
# bytes and the `required` check) but is moved off screen; a label styled as
# a button opens it, and a span next to it names what was picked. The span is
# kept current by the script below, wired into the layout, so any page can
# use `_webui_file_input` without further work.
const _WEBUI_FILE_FIELD_SCRIPT = """<script>
(function () {
  var describe = function (input) {
    var field = input.closest('[data-file-field]');
    if (field === null) return;
    var name = field.querySelector('[data-file-field-name]');
    if (name === null) return;
    var files = input.files;
    if (!files || files.length === 0) {
      name.textContent = 'No file selected';
    } else if (files.length === 1) {
      name.textContent = files[0].name;
    } else {
      name.textContent = files.length + ' files selected';
    }
  };
  document.addEventListener('change', function (event) {
    var target = event.target;
    if (target && target.matches && target.matches('input[type=file][data-file-field-input]')) describe(target);
  });
  // bfcache restores the input's selection but not our text
  window.addEventListener('pageshow', function () {
    document.querySelectorAll('input[type=file][data-file-field-input]').forEach(describe);
  });
})();
</script>"""

"""
    _webui_file_input(name; accept, multiple, required, id) -> String

A file picker whose button and status text are in the page's language rather
than the browser's. `id` defaults to `name`; pass one when a page carries two
pickers with the same field name.
"""
function _webui_file_input(name::AbstractString; accept::AbstractString = "", multiple::Bool = false, required::Bool = false, id::AbstractString = name)::String
  input_id = "file-field-$(id)"
  attrs = string(
    isempty(accept) ? "" : " accept=\"$(_webui_escape(accept))\"",
    multiple ? " multiple" : "",
    required ? " required" : "",
  )
  return string(
    "<span class=\"file-field\" data-file-field>",
    "<input id=\"$(input_id)\" class=\"file-field-input\" data-file-field-input type=\"file\" name=\"$(_webui_escape(name))\"$(attrs)>",
    "<label for=\"$(input_id)\" class=\"button secondary-button file-field-button\">$(multiple ? "Choose files" : "Choose file")</label>",
    "<span class=\"file-field-name\" data-file-field-name>No file selected</span>",
    "</span>",
  )
end

"""Render a file-path dropdown while showing only each file's basename."""

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
  milliseconds_total = max(0, round(Int, elapsed * 1000))
  total_seconds, milliseconds = divrem(milliseconds_total, 1000)
  hours, remainder = divrem(total_seconds, 3600)
  minutes, secs = divrem(remainder, 60)
  return lpad(hours, 2, '0') * ":" * lpad(minutes, 2, '0') * ":" * lpad(secs, 2, '0') * "." * lpad(milliseconds, 3, '0')
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

function _webui_phase_elapsed_seconds(result::AbstractDict, phase::AbstractString)
  timings = get(result, "service_phase_timings", Any[])
  timings isa AbstractVector || return nothing
  for timing in timings
    timing isa AbstractDict || continue
    get(timing, "phase", "") == phase || continue
    elapsed = get(timing, "elapsed_seconds", nothing)
    elapsed === nothing && return nothing
    parsed = elapsed isa Number ? Float64(elapsed) : tryparse(Float64, string(elapsed))
    return parsed !== nothing && isfinite(parsed) ? max(0.0, parsed) : nothing
  end
  return nothing
end

function _webui_solver_elapsed_seconds(result::AbstractDict)
  elapsed = get(result, "solver_elapsed_s", get(get(result, "metadata", Dict{String,Any}()), "solver_elapsed_s", nothing))
  elapsed === nothing && return nothing
  parsed = elapsed isa Number ? Float64(elapsed) : tryparse(Float64, string(elapsed))
  return parsed !== nothing && isfinite(parsed) ? max(0.0, parsed) : nothing
end

function _webui_total_elapsed_seconds(result::AbstractDict)
  total = _webui_phase_elapsed_seconds(result, "total_service")
  total !== nothing && return total
  elapsed = get(result, "elapsed_seconds", nothing)
  elapsed === nothing && return nothing
  parsed = elapsed isa Number ? Float64(elapsed) : tryparse(Float64, string(elapsed))
  return parsed !== nothing && isfinite(parsed) ? max(0.0, parsed) : nothing
end

# Everything the fallback Info menu shows is fixed for the life of the
# process, and building it is not free: resolving the case directory probes
# it with a temporary file. Without this memo the result page would do that
# probe on every one of its two-second auto-refreshes.
const _WEBUI_DEFAULT_INFO_MENU = Ref("")

function _webui_default_info_menu()::String
  if isempty(_WEBUI_DEFAULT_INFO_MENU[])
    root = default_webui_output_root()
    _WEBUI_DEFAULT_INFO_MENU[] = _webui_powerflow_info_menu(;
      output_root = root,
      config_file = DEFAULT_SPARLECTRA_CONFIG_PATH,
      case_directory = _webui_case_directory(),
      operation_log = webui_operation_log_path(root),
    )
  end
  return _WEBUI_DEFAULT_INFO_MENU[]
end

function _webui_layout(title::AbstractString, content::AbstractString; show_back::Bool = false, main_class::AbstractString = "page", refresh_url = nothing, refresh_seconds::Integer = WEBUI_STATUS_AUTO_REFRESH_SECONDS, header_info::AbstractString = "")::String
  back_button = show_back ? "<div class=\"page-toolbar\"><a class=\"button back-button\" href=\"/powerflow\" onclick=\"if (document.referrer.startsWith(location.origin)) { history.back(); return false; }\" aria-label=\"Go back to the previous page\">← Back</a></div>" : ""
  version_text = "Sparlectra.jl v$(version())"
  package_path = _sparlectra_package_path()
  commit_sha = _sparlectra_git_commit_sha()
  commit_text = commit_sha === nothing || isempty(strip(String(commit_sha))) ? "" : first(commit_sha, min(7, length(commit_sha)))
  commit_html = isempty(commit_text) ? "" : "<span class=\"runtime-commit\">commit $(_webui_escape(commit_text))</span>"
  # start flavor: standalone app / sysimage (each with its build
  # time) or a plain native session; detection never throws
  flavor = try
    webui_runtime_flavor()
  catch
    (kind = :native, built = nothing)
  end
  flavor_label = flavor.kind === :app ? "standalone app" : flavor.kind === :sysimage ? "sysimage" : "native session"
  flavor_built = flavor.built === nothing ? "" : ", built $(replace(first(String(flavor.built), 16), "T" => " "))"
  flavor_html = "<span class=\"runtime-flavor\">$(_webui_escape(string(flavor_label, flavor_built)))</span>"
  # startup latency hint (native sessions only, dismissible, suppressible
  # via output.startup_latency_hint = false); sysimage/app sessions never
  # see it, and a dismiss sticks per browser via localStorage
  latency_banner = ""
  if flavor.kind === :native
    hint_on = try
      active_sparlectra_config().output.startup_latency_hint
    catch
      true
    end
    hint_on && (latency_banner = "<div class=\"latency-banner\" id=\"latency-banner\"><span>First run in a fresh Julia process includes compilation and is slow; later runs are fast. Keep the process running, or start it again and let it build the sysimage.</span><button type=\"button\" onclick=\"try{localStorage.setItem('sparlectra-latency-dismissed','1')}catch(e){};document.getElementById('latency-banner').remove()\">Dismiss</button></div><script>try{localStorage.getItem('sparlectra-latency-dismissed')==='1'&&document.getElementById('latency-banner').remove()}catch(e){}</script>")
  end
  runtime_info = "<span class=\"runtime-info\" title=\"Package path: $(_webui_escape(package_path))\"><span class=\"runtime-title\">$(_webui_escape(version_text))</span>$(commit_html)$(flavor_html)</span>"
  # The Info control belongs to the HEADER, not to three particular pages.
  # It used to be passed in by the Case, Settings and Runs renderers only, so
  # it vanished the moment the user clicked Operation Log, Run history, Last
  # errors, Docs or a result page (reported 2026-09-07). Everything it shows
  # is a property of the running server, not of the page, so a page that
  # knows better still overrides `header_info` and every other page now gets
  # the same menu built from the defaults those pages would resolve to.
  isempty(header_info) && (header_info = _webui_default_info_menu())
  refresh_attrs = refresh_url === nothing ? "" : " data-refresh-url=\"$(_webui_escape(refresh_url))\" data-refresh-seconds=\"$(refresh_seconds)\""
  return """<!doctype html>
<html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">
<title>$(_webui_escape(title)) · Sparlectra</title><link rel=\"stylesheet\" href=\"/static/sparlectra.css\"></head>
<body><header class=\"site-header\"><a class=\"brand\" href=\"/powerflow\"><img class=\"brand-logo\" src=\"/assets/logo.png\" alt=\"Sparlectra.jl logo\">$(runtime_info)</a><nav><a href=\"/powerflow/case\">Case</a><a href=\"/powerflow/settings\">Settings</a><a href=\"/powerflow\">Runs</a><a href=\"/powerflow/history\">Run history</a><a href=\"/webui/operation-log\">Operation Log</a><a href=\"/docs\">Docs</a><a class=\"project-docs-link\" href=\"https://welthulk.github.io/Sparlectra.jl/\" target=\"_blank\" rel=\"noopener noreferrer\"><svg class=\"github-icon\" viewBox=\"0 0 16 16\" aria-hidden=\"true\"><path fill=\"currentColor\" d=\"M8 0C3.58 0 0 3.64 0 8.13c0 3.59 2.29 6.64 5.47 7.72.4.08.55-.18.55-.39 0-.19-.01-.83-.01-1.51-2.01.38-2.53-.5-2.69-.96-.09-.23-.48-.96-.82-1.15-.28-.15-.68-.53-.01-.54.63-.01 1.08.59 1.23.83.72 1.23 1.87.88 2.33.67.07-.53.28-.88.51-1.08-1.78-.21-3.64-.91-3.64-4.02 0-.89.31-1.62.82-2.19-.08-.2-.36-1.04.08-2.16 0 0 .67-.22 2.2.84A7.4 7.4 0 0 1 8 3.93c.68 0 1.36.09 2 .27 1.53-1.06 2.2-.84 2.2-.84.44 1.12.16 1.96.08 2.16.51.57.82 1.3.82 2.19 0 3.12-1.87 3.81-3.65 4.02.29.25.54.74.54 1.5 0 1.08-.01 1.95-.01 2.22 0 .22.15.47.55.39A8.15 8.15 0 0 0 16 8.13C16 3.64 12.42 0 8 0Z\"/></svg><span>Project Docs</span></a>$(header_info)<a href=\"/webui/last-errors\">Last errors</a><form method="post" action="/webui/shutdown" class="exit-form"><button type="submit" class="exit-button">Stop Web UI</button></form></nav></header>
$(latency_banner)<main class="$(main_class)"$(refresh_attrs)>$(back_button)<h1>$(_webui_escape(title))</h1>$(content)</main><footer>$(_webui_escape(version_text)) · Local PowerFlow Web UI · loopback access only</footer>
<script>
(function () {
  const sendHeartbeat = function () {
    fetch('/webui/heartbeat', {method: 'POST', keepalive: true}).catch(function () {});
  };
  sendHeartbeat();
  window.setInterval(sendHeartbeat, 5000);
  const scheduleAutoRefresh = function () {
    const mainEl = document.querySelector('main[data-refresh-url]');
    if (mainEl === null) return;
    const seconds = parseFloat(mainEl.getAttribute('data-refresh-seconds')) || 2;
    window.setTimeout(function () {
      const target = document.querySelector('main[data-refresh-url]');
      if (target === null) return;
      fetch(target.getAttribute('data-refresh-url')).then(function (response) {
        return response.text();
      }).then(function (html) {
        const newMain = new DOMParser().parseFromString(html, 'text/html').querySelector('main');
        if (newMain === null) return;
        if (newMain.hasAttribute('data-refresh-url')) {
          // Still an interim page (warming up / run active): swap in place to
          // keep polling without stacking history entries.
          target.replaceWith(document.importNode(newMain, true));
          scheduleAutoRefresh();
        } else {
          // Terminal page reached. Do a REAL reload: scripts inside a
          // DOMParser-imported <main> are inert per the HTML spec (never
          // executed), so swapping the final page in place would render it
          // with dead JS — the case-combobox arrow, solver toggles and every
          // other in-page handler would not react until the user forced a
          // navigation (the reported "works only after typing + Enter").
          window.location.reload();
        }
      }).catch(function () {
        scheduleAutoRefresh();
      });
    }, seconds * 1000);
  };
  scheduleAutoRefresh();
})();
</script>$(_WEBUI_BUSY_FORM_SCRIPT)$(_WEBUI_FILE_FIELD_SCRIPT)</body></html>"""
end

function _webui_powerflow_info_menu(; output_root::AbstractString, config_file::AbstractString, case_directory::AbstractString, operation_log::AbstractString)::String
  # Which build is actually serving this page. The header carries the same
  # three facts, but a user reporting a defect reads the info panel, and
  # "am I running the code I just built" has to be answerable there: a
  # sysimage that predates the fix looks exactly like one that contains it.
  commit_sha = _sparlectra_git_commit_sha()
  commit_text = commit_sha === nothing || isempty(strip(String(commit_sha))) ? "unknown" : String(commit_sha)
  flavor = try
    webui_runtime_flavor()
  catch
    (kind = :native, built = nothing)
  end
  flavor_text = if flavor.kind === :native
    "native session (no sysimage)"
  else
    label = flavor.kind === :app ? "standalone app" : "sysimage"
    flavor.built === nothing ? label : string(label, ", built ", replace(first(String(flavor.built), 19), "T" => " "))
  end
  return """
<details class=\"topbar-info-menu\">
<summary>Info</summary>
<div class=\"topbar-info-panel\">
<h2>Run information</h2>
<dl>
<dt>Version</dt><dd><code>Sparlectra.jl v$(_webui_escape(string(version())))</code></dd>
<dt>Commit</dt><dd><code>$(_webui_escape(commit_text))</code></dd>
<dt>Started from</dt><dd><code>$(_webui_escape(flavor_text))</code> <a href=\"/webui/sysimage\">Sysimage</a></dd>
<dt>Output root</dt><dd><code>$(_webui_escape(output_root))</code></dd>
<dt>Config file</dt><dd><code>$(_webui_escape(config_file))</code></dd>
<dt>Case cache</dt><dd><code>$(_webui_escape(case_directory))</code></dd>
<dt>Operation log</dt><dd><code>$(_webui_escape(operation_log))</code></dd>
</dl>
<p class=\"matpower-citation-note\">Using MATPOWER cases or data in publications? Please follow the <a href=\"https://matpower.org/citing/\" target=\"_blank\" rel=\"noopener noreferrer\">MATPOWER citation guidance</a>; the full reference is in the <a href=\"/docs\">documentation</a> (MATPOWER import, Citation).</p>
</div>
</details>"""
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
  # the banner names the calculation, not "PowerFlow" for everything
  kind = string(get(active_run, "kind_label", "PowerFlow run"))
  explanation = status == "aborting" ? "<p>Abort requested; the calculation stops at its next iteration. You can start a new run right away.</p>" : ""
  abort_form = status == "aborting" ? "" : "<form method=\"post\" action=\"/powerflow/abort/$(_webui_urlencode(run_id))\"><button type=\"submit\" class=\"danger-button\">Abort</button></form>"
  return "<section class=\"panel active-run-banner\"><strong>$(_webui_escape(kind)) is $(status):</strong> <code>$(_webui_escape(run_id))</code>$(explanation)<div class=\"actions\"><a class=\"button\" href=\"/powerflow/result/$(_webui_urlencode(run_id))\">Open status</a>$(abort_form)</div></section>"
end

function _webui_config_notice_html(notice)::String
  notice === nothing && return ""
  notice isa AbstractDict || return ""
  changed = get(notice, "changed", false)
  duplicates = get(notice, "duplicate_keys", String[])
  missing = get(notice, "missing_keys", String[])
  normalized = get(notice, "normalized_keys", String[])
  (changed || !isempty(duplicates) || !isempty(missing) || !isempty(normalized)) || return ""
  severity = !isempty(duplicates) ? "warning strong-warning" : !isempty(normalized) ? "warning" : "info"
  detail = !isempty(duplicates) ? " Duplicate YAML keys require manual review before refresh can write." : ""
  return "<div class=\"alert config-notice $(severity)\" role=\"status\"><strong>Configuration notice:</strong> your selected configuration is missing newer options or contains deprecated values.$(detail) Open <a href=\"/powerflow/settings#configuration-maintenance\">the advanced configuration tools</a> to review or refresh it.</div>"
end

function _webui_recent_error_entries(operation_log::AbstractString; limit::Integer = 5)::Vector{Dict{String,Any}}
  operation_log_path = webui_operation_log_path(operation_log)
  isfile(operation_log_path) || return Dict{String,Any}[]
  entries = Dict{String,Any}[]
  for line in eachline(operation_log_path)
    isempty(strip(line)) && continue
    event = try
      _parse_service_json(line)
    catch
      nothing
    end
    event isa AbstractDict || continue
    string(get(event, "event", "")) in ("validation_error", "powerflow_submit_rejected", "powerflow_failed", "internal_error") || continue
    push!(entries, Dict{String,Any}(String(key) => value for (key, value) in event))
  end
  length(entries) <= limit && return entries
  return entries[(end-limit+1):end]
end

function _webui_last_errors_list_html(operation_log::AbstractString)::String
  entries = _webui_recent_error_entries(operation_log)
  isempty(entries) && return "<p>No recent errors.</p>"
  items = join((begin
    timestamp = get(entry, "timestamp", "unknown time")
    route = get(entry, "route", "unknown route")
    casefile = get(entry, "requested_case", get(entry, "casefile", ""))
    run_id = get(entry, "run_id", "")
    reason = get(entry, "reason", "")
    message = get(entry, "message", get(entry, "status", "error"))
    meta = String[]
    isempty(string(casefile)) || push!(meta, "case: $(casefile)")
    isempty(string(run_id)) || push!(meta, "run: $(run_id)")
    isempty(string(reason)) || push!(meta, "reason: $(reason)")
    suffix = isempty(meta) ? "" : " <small>($(_webui_escape(join(meta, ", "))))</small>"
    result_link = isempty(string(run_id)) ? "" : " <a href=\"/powerflow/result/$(_webui_urlencode(string(run_id)))\">Open result</a>"
    "<li><time>$(_webui_escape(timestamp))</time> <code>$(_webui_escape(route))</code>: $(_webui_escape(message))$(suffix)$(result_link)</li>"
  end for entry in reverse(entries)), "")
  return "<ul class=\"last-errors-list\">$(items)</ul>"
end

function render_webui_last_errors(operation_log::AbstractString)::String
  panel = "<section class=\"panel last-errors-panel\">$(_webui_last_errors_list_html(operation_log))</section>"
  return _webui_layout("Last errors", panel; show_back = true)
end

function _webui_error_alert_html(error_message)::String
  error_message === nothing && return ""
  return "<div class=\"alert alert-error error\" role=\"alert\"><span>$(_webui_escape(error_message))</span></div>"
end

"""
    _webui_feedback_modal_html(parts) -> String

Combine one-off feedback/error messages (submission errors, case-import
results, case-settings-loaded notices) into a single dismissible `<dialog>`
popup instead of pushing them inline into the page flow. Returns `""` when
every part is empty, so no empty modal is rendered.
"""
function _webui_feedback_modal_html(parts::AbstractVector{<:AbstractString})::String
  content = join(filter(!isempty, parts), "")
  isempty(content) && return ""
  return "<dialog id=\"feedback-modal\" class=\"feedback-modal\" aria-label=\"PowerFlow notices\"><form method=\"dialog\" class=\"feedback-modal-close-form\"><button type=\"submit\" class=\"feedback-modal-close\" aria-label=\"Close\">×</button></form><div class=\"feedback-modal-body\">$(content)</div></dialog>"
end

"""
    _webui_case_context(; kwargs...) -> NamedTuple

Resolve everything the Web UI pages share about the selected case in one
place: saved-settings precedence (`webui_form_state`), the case list and
effective case directory, the combobox value split (list value vs free-typed
value), the format hint (delegating to `_detect_case_format`, see
`_webui_case_format_hint`), short-circuit button gating, and the SCF-only
scenario gate. Stage 4A extracts this from `render_powerflow_form` so the
Case, Settings, and Runs pages consume one shared context instead of
re-deriving it (and drifting apart) per page.
"""
function _webui_case_context(;
  application_root::AbstractString = _webui_application_root(),
  case_directory::Union{Nothing,AbstractString} = nothing,
  selected_casefile::AbstractString = "",
  selected_config_file::AbstractString = "",
  case_profile = nothing,
  submitted_form = nothing,
  show_case_settings_notice::Bool = true,
)
  profile_values = webui_form_state(; selected_casefile, selected_config_file, sidecar_profile = case_profile, submitted_form, case_directory)
  profile_path = String(get(profile_values, "_profile_path", ""))
  profile_location = isempty(profile_path) ? "the case configuration file" : "<code>$(_webui_escape(profile_path))</code>"
  profile_notice = if isempty(profile_path) || !show_case_settings_notice
    ""
  else
    "<div class=\"alert info case-settings-notice\" role=\"status\"><strong>Case-specific settings loaded from $(profile_location).</strong> Saved Web UI settings prefilled the form. Manual edits on this page override the profile for this run. <form method=\"post\" action=\"/powerflow/case-settings/reset\" class=\"case-settings-notice-dismiss\"><input type=\"hidden\" name=\"casefile\" value=\"$(_webui_escape(selected_casefile))\"><button type=\"submit\" class=\"link-button\" title=\"Delete the saved settings for this case so the form falls back to the configuration defaults. The case file itself is kept.\">Reset saved settings</button></form> <form method=\"post\" action=\"/powerflow/config/dismiss-case-settings-notice\" class=\"case-settings-notice-dismiss\"><input type=\"hidden\" name=\"config_file\" value=\"$(_webui_escape(selected_config_file))\"><input type=\"hidden\" name=\"casefile\" value=\"$(_webui_escape(selected_casefile))\"><button type=\"submit\" class=\"link-button\">Don't show this again</button></form></div>"
  end
  # The case file's own settings are visible, not silent: they moved the
  # controls, and the run will therefore use them. Without the note a user
  # sees a solver or tolerance they never chose and cannot tell why.
  case_file_notice = let fields = get(profile_values, "_case_file_fields", nothing)
    if fields isa AbstractVector && !isempty(fields)
      labels = String[]
      for f in fields
        spec = get(_WEBUI_OPTION_BY_FIELD, String(f), nothing)
        spec === nothing && continue
        push!(labels, string(spec.config_key === nothing ? String(f) : String(spec.config_key), " = ", _webui_form_string(get(profile_values, String(f), ""))))
      end
      items = join("<li><code>$(_webui_escape(l))</code></li>" for l in labels)
      "<div class=\"alert info case-file-settings-notice\" role=\"status\"><strong>This case file brings its own settings.</strong> The form below was prefilled from it, so the run uses what the case ships with:<ul class=\"config-override-list\">$(items)</ul>Edit any control to override it for this run.</div>"
    else
      ""
    end
  end
  casefiles = case_directory === nothing ? _webui_casefile_options(application_root) : _webui_casefile_options_in_directory(case_directory)
  # the shipped SCF demo cases (data/scf) are always offered; on first use
  # the run path stages them into the cache with their sidecars
  casefiles = sort!(unique!(vcat(casefiles, _webui_bundled_scf_options(application_root))); by = lowercase)
  bundled_case_directory = joinpath(application_root, "data", "mpower")
  effective_case_directory = case_directory === nothing ? bundled_case_directory : String(case_directory)
  for002_candidates = _webui_for002_reference_options_in_directory(effective_case_directory)
  selected_value = strip(selected_casefile)
  existing_value = selected_value in casefiles ? selected_value : ""
  manual_value = isempty(existing_value) ? selected_value : ""
  if submitted_form isa AbstractDict
    submitted_existing = strip(_webui_form_string(get(profile_values, "casefile", existing_value)))
    submitted_manual = strip(_webui_form_string(get(profile_values, "casefile_manual", manual_value)))
    existing_value = submitted_existing in casefiles ? submitted_existing : ""
    # The combined case input posts free-typed values under "casefile"; keep
    # them visible on re-render even though they are not in the cache list.
    manual_value = isempty(submitted_manual) && !(submitted_existing in casefiles) ? submitted_existing : submitted_manual
  end
  explicit_case_format = submitted_form isa AbstractDict && _webui_form_value(submitted_form, "case_format", nothing) !== nothing
  effective_case_value = isempty(strip(manual_value)) ? existing_value : manual_value
  format_hint = _webui_case_format_hint(effective_case_value; case_directory = effective_case_directory)
  case_format_value = if explicit_case_format
    strip(_webui_form_string(_webui_form_value(submitted_form, "case_format", "auto")))
  elseif format_hint == :dtf_for001
    "dtf_for001"
  elseif format_hint == :cgmes
    "cgmes"
  else
    # richer preselection (an explicit scf or matpower badge) is 4B
    # material; in 4A the form keeps today's auto default for the rest
    "auto"
  end
  dat_case_assistance = format_hint == :dtf_for001
  # "Short circuit" button gating (selectable only when
  # the data is there): CGMES case + a delivery that actually carries
  # short-circuit source data. An unresolvable path (e.g. a cgmes: alias not
  # fetched yet) stays enabled — the service run reports missing data with an
  # explicit reason either way.
  sc_is_cgmes = case_format_value == "cgmes" || format_hint == :cgmes
  # a case file (#342) carries its own IEC 60909 sources, so it qualifies too
  sc_is_cgmes = sc_is_cgmes || endswith(lowercase(strip(effective_case_value)), ".json")
  sc_has_data = if !sc_is_cgmes
    false
  else
    sc_resolved = try
      _resolve_powerflow_casefile(effective_case_value, effective_case_directory)
    catch
      nothing
    end
    sc_resolved === nothing ? true : _webui_case_has_short_circuit_data(sc_resolved)
  end
  sc_state = !sc_is_cgmes ? "not-cgmes" : (sc_has_data ? "ready" : "missing-data")
  # the scenario gate is a NAMING convention, not format detection: only the
  # canonical <case>.scf.json name carries an editable scenarios block
  scen_is_scf = endswith(lowercase(String(selected_casefile)), ".scf.json")
  return (; profile_values, profile_path, profile_notice, case_file_notice, casefiles, effective_case_directory, for002_candidates, existing_value, manual_value, effective_case_value, format_hint, case_format_value, dat_case_assistance, sc_state, scen_is_scf)
end

# Shared page scripts (stage 4A): the monolith's inline script is split into
# page-scoped building blocks so the Case, Settings, and Runs pages carry
# only what their controls need. Each block installs its own
# DOMContentLoaded listener; they are independent.
const _WEBUI_FEEDBACK_MODAL_SCRIPT = """
<script>
document.addEventListener('DOMContentLoaded', function () {
  const feedbackModal = document.getElementById('feedback-modal');
  if (feedbackModal !== null) {
    feedbackModal.showModal();
    feedbackModal.addEventListener('click', function (event) {
      if (event.target === feedbackModal) {
        feedbackModal.close();
      }
    });
  }
});
</script>"""

const _WEBUI_INFO_MENU_SCRIPT = """
<script>
document.addEventListener('DOMContentLoaded', function () {
  document.addEventListener('click', function (event) {
    document.querySelectorAll('.topbar-info-menu[open]').forEach(function (menu) {
      if (!menu.contains(event.target)) {
        menu.removeAttribute('open');
      }
    });
  });
});
</script>"""

# The case chooser script: editable combobox (open/filter/keyboard), the
# reload-on-selection with loading banner, right-click delete, Enter-resolve
# of unknown names, and the client-side format assistance (DTF hint, CGMES
# vs MATPOWER applicability, format auto-set). Moved verbatim from the run
# monolith in stage 4A with two deliberate changes: the reload target is
# parameterized (the chooser lives on the Case page now) and the resolve
# submit uses the input's OWN form instead of the run form.
function _webui_case_chooser_script(; reload_path::AbstractString = "/powerflow/case")::String
  return """
<script>
document.addEventListener('DOMContentLoaded', function () {
  const caseInput = document.querySelector('input[name="casefile"][data-case-settings-reload="true"]');
  const caseInputInitial = caseInput === null ? '' : caseInput.value.trim();
  const caseComboboxList = document.getElementById('case-combobox-list');
  const caseComboboxToggle = document.getElementById('case-combobox-toggle');
  const caseComboboxOptions = caseComboboxList === null ? [] : Array.prototype.slice.call(caseComboboxList.querySelectorAll('li[data-case-option]'));
  const availableCases = caseComboboxOptions.map(function (item) { return item.getAttribute('data-case-option'); });
  const caseFormat = document.querySelector('select[name="case_format"]');
  const dtfInternalSection = document.querySelector('.dtf-internal-section');
  const datFormatHint = document.getElementById('dtf-dat-format-hint');
  // MATPOWER import conventions do not steer the CGMES importer: gray them
  // out (disable, keep in place) while a CGMES case is selected. Disabled
  // controls drop out of the submitted form, so the save simply keeps the
  // stored values instead of stale Web UI values.
  const updateImportConventionApplicability = function () {
    const isCgmes = caseFormat !== null && caseFormat.value === 'cgmes';
    document.querySelectorAll('[data-matpower-import-field]').forEach(function (el) {
      el.classList.toggle('disabled', isCgmes);
      el.querySelectorAll('input, select').forEach(function (control) { control.disabled = isCgmes; });
    });
    const importHint = document.querySelector('[data-import-conventions-hint]');
    if (importHint !== null) importHint.hidden = !isCgmes;
    // Inverse gating for the CGMES start-values select: it only steers CGMES
    // runs, so it is enabled exactly when the selected case is CGMES and
    // grayed out (kept in place) otherwise.
    document.querySelectorAll('[data-cgmes-start-values-field]').forEach(function (el) {
      el.classList.toggle('disabled', !isCgmes);
      // Checkboxes in this group (require_boundary) gate like the select;
      // their hidden false-fallback inputs must stay enabled so the form
      // still submits an explicit value.
      el.querySelectorAll('select, input[type="checkbox"]').forEach(function (control) { control.disabled = !isCgmes; });
    });
    // "Short circuit" button gating only exists on pages that render the
    // button; on the Case page this is a no-op.
    const scButton = document.querySelector('[data-short-circuit-button]');
    if (scButton !== null) {
      const scDataMissing = scButton.getAttribute('data-sc-state') === 'missing-data';
      scButton.disabled = !isCgmes || scDataMissing;
    }
  };
  const updateDatCaseAssistance = function () {
    const effectiveValue = caseInput === null ? '' : caseInput.value.trim();
    const isDatCase = new RegExp('\\\\.dat\$', 'i').test(effectiveValue);
    // client-side CGMES markers mirror the server heuristic where the browser
    // can: explicit cgmes: alias input or a .zip delivery (directory paths
    // resolve server-side only)
    const isCgmesCase = new RegExp('^cgmes:', 'i').test(effectiveValue) || new RegExp('\\\\.zip\$', 'i').test(effectiveValue);
    if (caseFormat !== null) {
      // Auto-set formats must fall BACK to auto when the typed case stops
      // matching — otherwise a CGMES selection sticks after switching to a
      // MATPOWER case. Only formats this automation set itself are
      // reverted; a manual choice stays untouched.
      if (isDatCase) {
        caseFormat.value = 'dtf_for001';
        caseFormat.dataset.autoFormat = 'dtf_for001';
      } else if (isCgmesCase) {
        caseFormat.value = 'cgmes';
        caseFormat.dataset.autoFormat = 'cgmes';
      } else if (caseFormat.dataset.autoFormat && caseFormat.value === caseFormat.dataset.autoFormat) {
        caseFormat.value = 'auto';
        delete caseFormat.dataset.autoFormat;
      }
    }
    if (dtfInternalSection !== null) {
      dtfInternalSection.classList.toggle('is-dat-selected', isDatCase);
      if (isDatCase) dtfInternalSection.open = true;
    }
    if (datFormatHint !== null) {
      datFormatHint.hidden = !isDatCase;
      datFormatHint.textContent = isDatCase ? '.DAT selected: using internal DTF diagnostics.' : '';
    }
    updateImportConventionApplicability();
  };
  if (caseFormat !== null) {
    caseFormat.addEventListener('change', function () {
      // a manual format choice overrides and clears the automation marker
      delete caseFormat.dataset.autoFormat;
      updateImportConventionApplicability();
    });
  }
  updateDatCaseAssistance();
  // bfcache restores (back-navigation from a result page) skip
  // DOMContentLoaded — re-evaluate the case-dependent gating there too.
  window.addEventListener('pageshow', updateDatCaseAssistance);
  const hideCaseLoadingBanner = function () {
    const banner = document.getElementById('case-loading-banner');
    if (banner !== null) banner.hidden = true;
  };
  // A page restored from the back/forward cache keeps the DOM as it was when
  // the user navigated away — including the loading state. Clear it, or the
  // form looks (and on a stricter style would be) stuck after Back.
  window.addEventListener('pageshow', function () {
    document.querySelectorAll('.case-loading').forEach(function (el) { el.classList.remove('case-loading'); });
    hideCaseLoadingBanner();
  });
  const reloadWithCase = function (value) {
    // Selecting a case reloads the whole page from the server (its saved
    // settings have to be applied). On a large case directory that takes a
    // moment during which the visible fields still show the OLD case's
    // values and then jump — so say plainly that a reload is running instead
    // of letting the user watch settings change under their hands.
    const form = (caseInput !== null && caseInput.form !== null ? caseInput.form : null) || document.querySelector('main');
    if (form !== null) {
      form.classList.add('case-loading');
      // Safety net: if the navigation never happens (blocked, cancelled, or
      // the same page is restored from the bfcache), drop the dimming again
      // so the form never looks stuck.
      window.setTimeout(function () { form.classList.remove('case-loading'); hideCaseLoadingBanner(); }, 15000);
      let banner = document.getElementById('case-loading-banner');
      if (banner === null) {
        banner = document.createElement('div');
        banner.id = 'case-loading-banner';
        banner.className = 'alert info case-loading-banner';
        banner.setAttribute('role', 'status');
        banner.innerHTML = '<span class="submit-spinner" aria-hidden="true"></span> Loading case settings — please wait…';
        form.parentNode.insertBefore(banner, form);
      }
      banner.hidden = false;
    }
    const target = new URL('$(reload_path)', window.location.origin);
    target.searchParams.set('casefile', value);
    const configInput = document.querySelector('input[name="config_file"]');
    if (configInput !== null && configInput.value !== '') {
      target.searchParams.set('config_file', configInput.value);
    }
    window.location.href = target.pathname + target.search;
  };
  if (caseInput !== null) {
    // Editable combobox: the arrow (or ArrowDown, or clicking the field)
    // opens the full cached-case list like the former dropdown, typing
    // filters it. Picking a cached case reloads the page so its saved
    // settings prefill the controls; committing an unknown name/path with
    // Enter resolves it (download/copy into the case cache) without
    // starting a PowerFlow run.
    const setCaseListOpen = function (open) {
      if (caseComboboxList === null) return;
      caseComboboxList.hidden = !open;
      caseInput.setAttribute('aria-expanded', open ? 'true' : 'false');
    };
    const filterCaseList = function (filterText) {
      const needle = String(filterText).trim().toLowerCase();
      caseComboboxOptions.forEach(function (item) {
        item.hidden = needle !== '' && item.getAttribute('data-case-option').toLowerCase().indexOf(needle) === -1;
        item.classList.remove('active');
      });
    };
    const visibleCaseOptions = function () {
      return caseComboboxOptions.filter(function (item) { return !item.hidden; });
    };
    const moveCaseActive = function (step) {
      const visible = visibleCaseOptions();
      if (visible.length === 0) return;
      let index = visible.findIndex(function (item) { return item.classList.contains('active'); });
      index = index === -1 ? (step > 0 ? 0 : visible.length - 1) : (index + step + visible.length) % visible.length;
      caseComboboxOptions.forEach(function (item) { item.classList.remove('active'); });
      visible[index].classList.add('active');
      visible[index].scrollIntoView({ block: 'nearest' });
    };
    const chooseCaseOption = function (item) {
      const value = item.getAttribute('data-case-option');
      caseInput.value = value;
      setCaseListOpen(false);
      updateDatCaseAssistance();
      if (value !== caseInputInitial) reloadWithCase(value);
    };
    if (caseComboboxToggle !== null) {
      caseComboboxToggle.addEventListener('mousedown', function (event) {
        event.preventDefault();
        if (caseComboboxList !== null && caseComboboxList.hidden) {
          filterCaseList('');
          setCaseListOpen(true);
          caseInput.focus();
        } else {
          setCaseListOpen(false);
        }
      });
    }
    caseInput.addEventListener('click', function () {
      filterCaseList('');
      setCaseListOpen(true);
    });
    caseComboboxOptions.forEach(function (item) {
      item.addEventListener('mousedown', function (event) {
        if (event.button !== 0) return;
        event.preventDefault();
        chooseCaseOption(item);
      });
      // Right-click deletes the case file (with confirmation) from the case
      // cache directory without starting a PowerFlow run.
      item.addEventListener('contextmenu', function (event) {
        event.preventDefault();
        const value = item.getAttribute('data-case-option');
        if (!window.confirm('Delete case file "' + value + '" from the case directory?')) return;
        const deleteForm = document.createElement('form');
        deleteForm.method = 'post';
        deleteForm.action = '/powerflow/delete-case';
        const deleteField = document.createElement('input');
        deleteField.type = 'hidden';
        deleteField.name = 'casefile';
        deleteField.value = value;
        deleteForm.appendChild(deleteField);
        document.body.appendChild(deleteForm);
        deleteForm.submit();
      });
    });
    document.addEventListener('mousedown', function (event) {
      const combobox = document.querySelector('[data-case-combobox]');
      if (combobox !== null && !combobox.contains(event.target)) setCaseListOpen(false);
    });
    caseInput.addEventListener('input', function () {
      filterCaseList(caseInput.value);
      setCaseListOpen(true);
      updateDatCaseAssistance();
    });
    caseInput.addEventListener('change', updateDatCaseAssistance);
    caseInput.addEventListener('keydown', function (event) {
      if (event.key === 'ArrowDown' || event.key === 'ArrowUp') {
        event.preventDefault();
        if (caseComboboxList !== null && caseComboboxList.hidden) {
          // Opening the closed list always shows every cached case (like a
          // select); filtering only applies while actively typing.
          filterCaseList('');
          setCaseListOpen(true);
        }
        moveCaseActive(event.key === 'ArrowDown' ? 1 : -1);
        return;
      }
      if (event.key === 'Escape') {
        setCaseListOpen(false);
        return;
      }
      if (event.key !== 'Enter') return;
      event.preventDefault();
      const active = visibleCaseOptions().find(function (item) { return item.classList.contains('active'); });
      if (active !== undefined && caseComboboxList !== null && !caseComboboxList.hidden) {
        chooseCaseOption(active);
        return;
      }
      const value = caseInput.value.trim();
      if (value === '') return;
      setCaseListOpen(false);
      if (availableCases.indexOf(value) >= 0) {
        if (value !== caseInputInitial) reloadWithCase(value);
        return;
      }
      const resolveButton = document.getElementById('resolve-case-button');
      if (resolveButton === null || caseInput.form === null) return;
      if (typeof caseInput.form.requestSubmit === 'function') {
        caseInput.form.requestSubmit(resolveButton);
      } else {
        resolveButton.click();
      }
    });
  }
});
</script>"""
end

# --- stage 4B: adapter option fieldsets generated from options_type -------
#
# The SET of rendered fields per adapter comes from the adapter's option
# struct (options_type) intersected with the :adapter-scope specs, matched
# by the config-key tail; shared model.* fields are claimed by the FIRST
# section that lists them (dedupe), so one page never renders a field
# twice. The order and the wording live in the presentation table below,
# and load-time asserts pin presentation and derivation to each other: a
# new struct field with a spec cannot ship without a presentation entry,
# and a presentation entry cannot outlive its field.
function _webui_adapter_struct_fields(adapter)::Vector{String}
  fields = String[]
  for fname in fieldnames(options_type(adapter))
    tail = String(fname)
    for spec in WEBUI_OPTION_SPECS
      spec.scope == :adapter || continue
      spec.config_key === nothing && continue
      if String(last(split(spec.config_key, "."))) == tail
        push!(fields, spec.field)
        break
      end
    end
  end
  return fields
end

# presentation: label, tooltip, marker attributes, input extras, and (for
# the two CGMES selects) curated option labels; the strings are the ones
# the hand-written fieldsets carried before 4B
const _WEBUI_ADAPTER_FIELD_PRESENTATION = Dict{String,NamedTuple}(
  "cgmes_start_values" => (label = "CGMES start values", title = "auto uses the delivery's own SvVoltage state when it carries one (real deliveries are built around their operating point) and falls back to the flat start otherwise.", attrs = " data-cgmes-start-values-field", input_attrs = "", option_labels = Dict("auto" => "Automatic (SV state when available, default)", "sv" => "Imported SV state", "flat" => "Flat start")),
  "cgmes_require_boundary" => (label = "Require boundary set", title = "Fail the CGMES import when topology references stay unresolved (boundary set missing). Uncheck to import an incomplete delivery anyway (buses without a resolvable BaseVoltage still abort).", attrs = " data-cgmes-start-values-field", input_attrs = "", option_labels = nothing),
  "cgmes_infer_base_voltages" => (label = "Infer missing base voltages", title = "Reconstruct missing nominal voltages when the delivery ships without its BaseVoltage catalog: seeded from the SV voltages and transformer rated voltages, propagated across level-preserving equipment. Substitutions are summarized as a warning. Pair with an unchecked Require boundary set.", attrs = " data-cgmes-start-values-field", input_attrs = "", option_labels = nothing),
  "cgmes_hvdc_mode" => (label = "HVDC converters", title = "How HVDC converters are modeled: fixed Stage-0 injections reproduce the delivery snapshot; a paired controller keeps the two converters of one link coupled (transfer, loss, terminal Q) and makes the transfer steerable.", attrs = " data-cgmes-start-values-field", input_attrs = "", option_labels = Dict("injections" => "Fixed injections (Stage 0, default)", "paired_control" => "Paired controller (steerable link)")),
  "matpower_import_auto_profile" => (label = "MATPOWER auto-profile", title = "", attrs = " data-matpower-import-field", input_attrs = "", option_labels = nothing),
  "matpower_import_ratio" => (label = "Transformer ratio convention", title = "", attrs = " data-matpower-import-field", input_attrs = "", option_labels = nothing),
  "matpower_import_shift_sign" => (label = "Phase-shift sign", title = "", attrs = " data-matpower-import-field", input_attrs = " step=\"2\" min=\"-1\" max=\"1\"", option_labels = nothing),
  "matpower_import_shift_unit" => (label = "Phase-shift unit", title = "", attrs = " data-matpower-import-field", input_attrs = "", option_labels = nothing),
  "matpower_import_bus_shunt_model" => (label = "Bus-shunt model", title = "", attrs = " data-matpower-import-field", input_attrs = "", option_labels = nothing),
  "matpower_import_dcline_mode" => (label = "DC-line mode", title = "How active mpc.dcline rows are modeled: pf_injections adds two fixed terminal injections per row (MATPOWER toggle_dcline equivalent); paired_control additionally couples each pair as a steerable HVDC link controller.", attrs = " data-matpower-import-field", input_attrs = "", option_labels = nothing),
  "matpower_import_pv_voltage_source" => (label = "PV voltage source", title = "", attrs = " data-matpower-import-field", input_attrs = "", option_labels = nothing),
  "matpower_import_compare_voltage_reference" => (label = "Voltage reference comparison", title = "", attrs = " data-matpower-import-field", input_attrs = "", option_labels = nothing),
  "transformer_tap_changer_model" => (label = "Tap-changer model", title = "", attrs = " data-ac-only-field data-matpower-import-field", input_attrs = "", option_labels = nothing),
  "matpower_import_apply_bus_names" => (label = "Apply bus names", title = "Use the case file's mpc.bus_name entries as bus names in results and logs instead of numeric BUS_I identifiers. Requires a bus_name block matching the bus count.", attrs = " data-matpower-import-field", input_attrs = "", option_labels = nothing),
  "matpower_export_write_solution" => (label = "Export Solution", title = "", attrs = "", input_attrs = "", option_labels = nothing),
)

# per-section field order (presentation) plus the extras that widen the
# form surface beyond the conversion struct: the auto-profile
# recommendation, the run-side comparison convention, the model default
# that both the MATPOWER and DTF structs consume, and the export
# convention; each is an :adapter-scope spec without a same-named struct
# field
const _WEBUI_ADAPTER_SECTIONS = (
  (adapter = CGMESAdapter(), key = :cgmes, order = ("cgmes_start_values", "cgmes_require_boundary", "cgmes_infer_base_voltages", "cgmes_hvdc_mode"), extras = ()),
  (adapter = MatpowerAdapter(), key = :matpower, order = ("matpower_import_auto_profile", "matpower_import_ratio", "matpower_import_shift_sign", "matpower_import_shift_unit", "matpower_import_bus_shunt_model", "matpower_import_dcline_mode", "matpower_import_pv_voltage_source", "matpower_import_compare_voltage_reference", "transformer_tap_changer_model", "matpower_import_apply_bus_names", "matpower_export_write_solution"), extras = ("matpower_import_auto_profile", "matpower_import_bus_shunt_model", "matpower_import_compare_voltage_reference", "matpower_export_write_solution")),
  (adapter = DTFAdapter(), key = :dtf, order = (), extras = ()),
  (adapter = PGMAdapter(), key = :pgm, order = (), extras = ()),
)

# derivation-vs-presentation asserts: the ordered presentation of a section
# must equal (struct-derived minus fields claimed earlier) plus extras
let claimed = Set{String}()
  for sec in _WEBUI_ADAPTER_SECTIONS
    derived = [f for f in _webui_adapter_struct_fields(sec.adapter) if !(f in claimed)]
    expected = union(Set(derived), Set(String.(collect(sec.extras))))
    got = Set(String.(collect(sec.order)))
    @assert got == expected "adapter section $(sec.key): presentation order $(sort(collect(got))) does not match derived+extras $(sort(collect(expected)))"
    union!(claimed, got)
    for f in sec.order
      @assert haskey(_WEBUI_ADAPTER_FIELD_PRESENTATION, f) "adapter field $(f) has no presentation entry"
    end
  end
end

# one generated field: reproduces the markup the hand-written fieldsets
# used (hidden-false checkbox pattern, _webui_select markup, help links)
function _webui_adapter_option_html(field::AbstractString, profile_values)::String
  spec = _webui_option_spec(field)
  pres = _WEBUI_ADAPTER_FIELD_PRESENTATION[String(field)]
  title_attr = isempty(pres.title) ? "" : " title=\"$(_webui_escape(pres.title))\""
  label = _webui_field_label(field, pres.label)
  if spec.control == :checkbox
    checked = _webui_checked(profile_values, field, _webui_option_default(field))
    return "<label class=\"check\"$(pres.attrs)$(title_attr)><input name=\"$(field)\" type=\"hidden\" value=\"false\"><input name=\"$(field)\" type=\"checkbox\" value=\"true\"$(checked)>$(label)</label>"
  elseif spec.control == :select
    selected = _webui_selected(profile_values, field, _webui_option_default(field))
    if pres.option_labels !== nothing
      opts = join(("<option value=\"$(v)\"$(String(selected) == String(v) ? " selected" : "")>$(_webui_escape(pres.option_labels[String(v)]))</option>" for v in spec.allowed_values), "")
      return "<label$(pres.attrs)$(title_attr)>$(label)<select name=\"$(field)\">$(opts)</select></label>"
    end
    return "<label$(pres.attrs)$(title_attr)>$(label)$(_webui_select(field, spec.allowed_values, selected))</label>"
  end
  value = _webui_input_value(profile_values, field, _webui_option_default(field))
  return "<label$(pres.attrs)$(title_attr)>$(label)<input name=\"$(field)\" type=\"number\"$(pres.input_attrs) value=\"$(value)\"></label>"
end

"""
    _webui_adapter_options_html(key, profile_values) -> String

The generated option group of one adapter section (stage 4B): the field
set is derived from `options_type(adapter)` and the :adapter-scope specs,
the order and wording come from the presentation table, both pinned by
load-time asserts. Returns "" for an adapter without form options (DTF,
PGM).
"""
function _webui_adapter_options_html(key::Symbol, profile_values)::String
  for sec in _WEBUI_ADAPTER_SECTIONS
    sec.key == key || continue
    isempty(sec.order) && return ""
    # visibility follows the spec section (stage 4B Basic list): :basic
    # fields render directly, :expert fields fold into a nested details
    basics = [f for f in sec.order if _webui_option_spec(f).section == :basic]
    experts = [f for f in sec.order if _webui_option_spec(f).section == :expert]
    parts = String[_webui_adapter_option_html(f, profile_values) for f in basics]
    if !isempty(experts)
      push!(parts, string("<details class=\"span-2 expert-section adapter-expert\"><summary>Advanced</summary>", join((_webui_adapter_option_html(f, profile_values) for f in experts), "\n"), "</details>"))
    end
    return join(parts, "\n")
  end
  return ""
end

"""
    render_case_page(; kwargs...) -> String

The Case page (stage 4A): choose or resolve a case (editable combobox with
right-click delete), upload case files, export the selected case (SCF or
plain PGM) and download it, and edit the per-case IMPORT options (input
format, CGMES import options, MATPOWER import conventions). Saving the
options writes them into the case configuration file next to the case; runs
pick them up through the configuration precedence (`resolve_config` for the
config keys, the request builder's form-block fallback for `case_format`),
not through the run POST.
"""
function render_case_page(;
  output_root::AbstractString = "results/powerflow_service",
  case_directory::Union{Nothing,AbstractString} = nothing,
  operation_log::AbstractString = webui_operation_log_path(output_root),
  error_message = nothing,
  application_root::AbstractString = _webui_application_root(),
  selected_casefile::AbstractString = "",
  selected_config_file::AbstractString = "",
  import_message::AbstractString = "",
  download_file::AbstractString = "",
  case_profile = nothing,
  submitted_form = nothing,
)::String
  ctx = _webui_case_context(; application_root, case_directory, selected_casefile, selected_config_file, case_profile, submitted_form, show_case_settings_notice = false)
  profile_values = ctx.profile_values
  effective_case_directory = ctx.effective_case_directory
  effective_case_value = ctx.effective_case_value
  case_format_value = ctx.case_format_value
  dat_case_assistance = ctx.dat_case_assistance
  error_html = _webui_error_alert_html(error_message)
  # a freshly exported file is offered right here: it was written next to the
  # case (where a run finds it), and this is the way back out to the browser
  download_html = isempty(strip(download_file)) ? "" :
                  " <a class=\"case-download-link\" href=\"/powerflow/case/download?case=$(_webui_urlencode(String(download_file)))\">Download $(_webui_escape(String(download_file)))</a>"
  import_html = isempty(strip(import_message)) && isempty(download_html) ? "" : "<div class=\"alert info case-import-result\" role=\"status\">$(_webui_escape(import_message))$(download_html)</div>"
  dtf_details_attrs = dat_case_assistance ? " class=\"span-2 dtf-internal-section is-dat-selected\" open" : " class=\"span-2 dtf-internal-section\""
  dat_hint_html = dat_case_assistance ? "<p id=\"dtf-dat-format-hint\" class=\"field-hint dat-format-hint span-2\" role=\"status\"><strong>.DAT selected:</strong> using internal DTF diagnostics.</p>" : "<p id=\"dtf-dat-format-hint\" class=\"field-hint dat-format-hint span-2\" role=\"status\" hidden></p>"
  case_options = join((begin
    has_settings = isfile(_webui_case_settings_path(output_root, casefile; case_directory = effective_case_directory))
    label = has_settings ? "$(casefile) ★" : casefile
    "<li role=\"option\" data-case-option=\"$(_webui_escape(casefile))\" title=\"Right-click to delete this case from the case directory\">$(_webui_escape(label))</li>"
  end for casefile in ctx.casefiles), "")
  case_input = "<span class=\"case-combobox\" data-case-combobox><input id=\"casefile\" name=\"casefile\" autocomplete=\"off\" spellcheck=\"false\" role=\"combobox\" aria-expanded=\"false\" aria-controls=\"case-combobox-list\" data-case-settings-reload=\"true\" value=\"$(_webui_escape(effective_case_value))\" placeholder=\"case14.m or /path/to/FOR001.DAT\"><button type=\"button\" id=\"case-combobox-toggle\" class=\"case-combobox-toggle\" aria-label=\"Show available cases\" tabindex=\"-1\">&#9662;</button><ul id=\"case-combobox-list\" class=\"case-combobox-list\" role=\"listbox\" hidden>$(case_options)</ul></span>"
  config_default = isempty(selected_config_file) ? DEFAULT_SPARLECTRA_CONFIG_PATH : selected_config_file
  info_menu = _webui_powerflow_info_menu(; output_root, config_file = config_default, case_directory = effective_case_directory, operation_log)
  import_form = """
<form id=\"case-import-form\" method=\"post\" action=\"/powerflow/import-cases\" enctype=\"multipart/form-data\" class=\"panel form-grid case-import-form\">
<label class=\"span-2\">$(_webui_field_label("casefiles", "Import case files"))$(_webui_file_input("casefiles"; accept = ".m,.M,.dat,.DAT,.zip,.ZIP,.json", multiple = true))</label>
<div class=\"actions span-2\"><button class=\"secondary-button\" type=\"submit\">Import case files</button></div>
</form>
"""
  # chooser + export live in ONE form: Enter in the combobox resolves via the
  # hidden resolve button, the export buttons post the same casefile value
  # through their formaction targets (unchanged endpoints)
  case_form = """
<form id=\"case-select-form\" method=\"post\" action=\"/powerflow/resolve-case\" class=\"panel form-grid case-select-form\">
<input type=\"hidden\" name=\"config_file\" value=\"$(_webui_escape(config_default))\">
<label class="span-2">$(_webui_field_label("casefile", "Case file"))$(case_input)<button type="submit" id="resolve-case-button" formaction="/powerflow/resolve-case" formmethod="post" formnovalidate hidden>Resolve case</button><small class="field-hint">Cases from <code>$(_webui_escape(effective_case_directory))</code> — pick one from the list, or type a case name/path and press Enter to download it into the list.<br>MATPOWER: a case name such as <code>case118.m</code>. CGMES: <code>cgmes:</code> plus one of $(join(("<code>" * a * "</code>" for a in sort(collect(keys(CGMESImporter.CGMES_TESTSET_ALIASES)))), ", ")) — fetches the ENTSO-E test configuration including its boundary set.</small></label>
$(dat_hint_html)
<div class=\"actions case-export-actions\">$(_webui_help_link("webui.scf_export", "Case export")) <button type=\"submit\" class=\"secondary-button\" formaction=\"/powerflow/export-scf\" formmethod=\"post\" formnovalidate title=\"Write the selected case as a Sparlectra Case Format file (.scf.json) into the case directory. One self-describing file: its data section is a valid power-grid-model input dataset, the namespaced sparlectra block carries slack roles, the tap-changer cascade, names and source ids, the measurements found next to the case, and the saved case options.\">Export as SCF case file</button><button type=\"submit\" class=\"secondary-button\" formaction=\"/powerflow/export-scf\" formmethod=\"post\" formnovalidate name=\"scf_strict_pgm\" value=\"true\" title=\"Write the PLAIN power-grid-model dataset (.pgm.json): only the data section, no namespaced sparlectra block. For handing the case to a PGM-only consumer - names, slack roles, tap nameplates, measurements and configuration are not in that file, and a slack generator is rewritten as a PGM source.\">Export as plain PGM</button><a class=\"button secondary-button case-download-link\" href=\"/powerflow/case/download?case=$(_webui_urlencode(selected_casefile))\" title=\"Download the selected case file from the case directory. After an export, pick the .scf.json or .pgm.json in the selector to download it.\">Download selected case</a></div>
</form>
"""
  options_form = """
<form id=\"case-options-form\" method=\"post\" action=\"/powerflow/case/options/save\" class=\"panel form-grid case-options-form\">
<input type=\"hidden\" name=\"casefile\" value=\"$(_webui_escape(effective_case_value))\">
<p class=\"lede span-2\">Import options for this case. Saving writes them into the case configuration file next to the case; every run of this case uses them through the configuration precedence.</p>
<details$(dtf_details_attrs)>
<summary>Input format</summary>
<fieldset>
<label>$(_webui_field_label("case_format", "Case input format"))<select name="case_format"><option value="auto"$(_webui_form_string(case_format_value) == "auto" ? " selected" : "")>Auto</option><option value="matpower"$(_webui_form_string(case_format_value) == "matpower" ? " selected" : "")>MATPOWER</option><option value="dtf_for001"$(_webui_form_string(case_format_value) == "dtf_for001" ? " selected" : "")>DTF diagnostics (experimental/internal)</option><option value="cgmes"$(_webui_form_string(case_format_value) == "cgmes" ? " selected" : "")>CGMES (ENTSO-E, folder or ZIP)</option><option value="scf"$(_webui_form_string(case_format_value) == "scf" ? " selected" : "")>Sparlectra Case Format (.scf.json)</option><option value="pgm"$(_webui_form_string(case_format_value) == "pgm" ? " selected" : "")>power-grid-model JSON (input.json)</option></select></label>
<p class="field-help">SCF and power-grid-model JSON are read by the same importer; the <code>sparlectra</code> block is optional, so a plain power-grid-model dataset loads as well. <em>Auto</em> already resolves every <code>.json</code> to that reader, so these two entries only matter when the extension does not say it.</p>
$(_webui_adapter_options_html(:cgmes, profile_values))
<p class="field-help" data-cgmes-start-values-field>CGMES only: <em>Flat start</em> lets the solver earn the solution itself; <em>Imported SV state</em> starts Newton-Raphson from the delivery's own SvVoltage solution (competing start-value machines are forced off). The SV comparison check (<code>sv_compare.csv</code>) runs either way.</p>
</fieldset>
</details>
<fieldset class=\"import-section\" data-import-conventions-section>
<legend>MATPOWER import conventions</legend>
<p class=\"field-hint span-2\" data-import-conventions-hint hidden>Not applicable to the selected CGMES case: these options steer MATPOWER (and DTF) parsing only. The CGMES importer reads the delivery's own conventions; Export Solution stays available.</p>
$(_webui_adapter_options_html(:matpower, profile_values))
</fieldset>
<div class=\"span-2 actions\"><button class=\"secondary-button\" type=\"submit\">Save case options</button></div>
</form>
"""
  content = """
$(_webui_feedback_modal_html([error_html, import_html]))<p class=\"lede\">Choose, import, export, and configure grid cases. Runs start on the <a href=\"/powerflow$(isempty(strip(effective_case_value)) ? "" : "?casefile=" * _webui_urlencode(effective_case_value))\">Runs page</a>.</p>
$(import_form)
$(case_form)
$(options_form)
$(_WEBUI_FEEDBACK_MODAL_SCRIPT)
$(_webui_case_chooser_script())
$(_WEBUI_INFO_MENU_SCRIPT)"""
  return _webui_layout("Case", content; header_info = info_menu)
end

"""
    _webui_settings_sections_html(; kwargs...) -> String

The solver, output, and expert option sections of the Settings page
(stage 4A block 3): moved VERBATIM from the run monolith. Values prefill
from the shared resolution chain (`profile_values`); the page's save
handler decides where they go (case configuration file or the general
YAML), the run POST no longer carries them.
"""
function _webui_settings_sections_html(; profile_values, config_default, profile_path::AbstractString, selected_casefile::AbstractString, selected_config_file::AbstractString)
  config_maintenance = """
  
<fieldset id="configuration-maintenance" class="config-maintenance">
<legend>$(_webui_field_label("config_maintenance", "Configuration maintenance"))</legend>
<div class="actions span-2"><a class="secondary-button" href="/powerflow/config/edit?config_file=$(_webui_urlencode(config_default))">Configuration Editor</a><button class="secondary-button" type="submit" formaction="/powerflow/config/check" formmethod="post">Check configuration</button><button class="secondary-button" type="submit" formaction="/powerflow/config/refresh" formmethod="post">Refresh configuration</button></div>
</fieldset>
"""
  return """
<label data-ac-only-field title=\"Convergence bound for the largest single bus mismatch (active and reactive alike). Readable in physical units as tol times the case base: 1e-8 pu equals 1 W at a 100 MVA base.\">$(_webui_field_label("power_flow_tol", "Tolerance"))<span class=\"tolerance-field\"><input name=\"power_flow_tol\" type=\"text\" autocomplete=\"off\" spellcheck=\"false\" data-tolerance-input value=\"$(_webui_input_value(profile_values, "power_flow_tol", _webui_option_default("power_flow_tol")))\"><span class=\"tolerance-spin\"><button type=\"button\" class=\"tolerance-spin-up\" data-tolerance-direction=\"up\" aria-label=\"Increase tolerance exponent\">&#9650;</button><button type=\"button\" class=\"tolerance-spin-down\" data-tolerance-direction=\"down\" aria-label=\"Decrease tolerance exponent\">&#9660;</button></span></span></label><label data-ac-only-field class=\"tolerance-unit\" title=\"The unit of the tolerance value. MW states the convergence bound physically: the run converts it with the case's own base (1 MW at a 100 MVA base is 1e-2 pu, 0.001 MW is 1 kW). pu states it per unit, the classical form. The value field to the left holds the number either way.\">$(_webui_field_label("power_flow_tol_unit", "Unit"))<select name=\"power_flow_tol_unit\"><option value=\"pu\"$(_webui_input_value(profile_values, "power_flow_tol_unit", _webui_option_default("power_flow_tol_unit")) == "MW" ? "" : " selected")>pu</option><option value=\"MW\"$(_webui_input_value(profile_values, "power_flow_tol_unit", _webui_option_default("power_flow_tol_unit")) == "MW" ? " selected" : "")>MW</option></select></label>
<label data-nr-only-field>$(_webui_field_label("power_flow_max_iter", "Maximum iterations"))<input name=\"power_flow_max_iter\" type=\"number\" min=\"1\" value=\"$(_webui_input_value(profile_values, "power_flow_max_iter", _webui_option_default("power_flow_max_iter")))\"></label>
<fieldset class=\"distributed-slack-options\" data-nr-only-field>
<legend>Distributed slack</legend>
<label class=\"check span-2\"><input name=\"power_flow_distributed_slack_enabled\" type=\"hidden\" value=\"false\"><input name=\"power_flow_distributed_slack_enabled\" type=\"checkbox\" value=\"true\"$(_webui_checked(profile_values, "power_flow_distributed_slack_enabled", _webui_option_default("power_flow_distributed_slack_enabled")))>$(_webui_field_label("power_flow_distributed_slack_enabled", "Distribute active-power slack over participating generators"))</label>
<label>$(_webui_field_label("power_flow_distributed_slack_p_mode", "Participation mode"))$(_webui_select("power_flow_distributed_slack_p_mode", _webui_option_allowed_values("power_flow_distributed_slack_p_mode"), _webui_selected(profile_values, "power_flow_distributed_slack_p_mode", _webui_option_default("power_flow_distributed_slack_p_mode"))))</label>
<p class=\"field-help\">The REF bus keeps the angle reference; the island's P imbalance is absorbed by participating generators via one λ<sub>P</sub> per island. <code>imported</code> reads participation factors from the case data (MATPOWER <code>APF</code>, CGMES <code>normalPF</code>). Explicit weights are YAML-only (<code>power_flow.distributed_slack.weights</code>).</p>
</fieldset>
<label>$(_webui_field_label("performance_timing", "Performance timing"))$(_webui_select("performance_timing", _webui_option_allowed_values("performance_timing"), _webui_selected(profile_values, "performance_timing", _webui_option_default("performance_timing"))))</label>
<details class=\"span-2 detailed-csv-options\">
<summary>Detailed result CSV export</summary>
<label class=\"check\"><input name=\"detailed_result_csv\" type=\"hidden\" value=\"false\"><input name=\"detailed_result_csv\" type=\"checkbox\" value=\"true\"$(_webui_checked(profile_values, "detailed_result_csv", _webui_option_default("detailed_result_csv")))>$(_webui_field_label("detailed_result_csv", "Export detailed result CSV files"))</label>
<label class=\"detailed-csv-format\">$(_webui_field_label("detailed_result_csv_format", "CSV format"))$(_webui_select("detailed_result_csv_format", _webui_option_allowed_values("detailed_result_csv_format"), _webui_selected(profile_values, "detailed_result_csv_format", _webui_option_default("detailed_result_csv_format"))))</label>
</details>
<label class=\"check span-2\"><input name=\"export_cgmes\" type=\"hidden\" value=\"false\"><input name=\"export_cgmes\" type=\"checkbox\" value=\"true\"$(_webui_checked(profile_values, "export_cgmes", _webui_option_default("export_cgmes")))>$(_webui_field_label("export_cgmes", "Export case as CGMES delivery (EQ+TP+SSH+SV, ZIP)"))</label>
<details class=\"span-2 expert-section\">
<summary>Advanced options</summary>
$(config_maintenance)
<fieldset class=\"step-control-expert\" data-ac-only-field>
<legend>Step control &amp; solver</legend>
<label class=\"check\" title=\"Inspect the imported network and pick the start-value, step-control, and Q-limit strategy automatically; on non-convergence a bounded escalation ladder retries with stronger strategies. Options you set explicitly below always win over the automatic choices. Decisions land in the auto_mode_decision.log artifact of the run.\"><input name=\"power_flow_mode\" type=\"hidden\" value=\"manual\"><input name=\"power_flow_mode\" type=\"checkbox\" value=\"auto\"$(_webui_selected(profile_values, "power_flow_mode", "manual") == "auto" ? " checked" : "")>$(_webui_field_label("power_flow_mode", "Auto mode (network-driven strategy)"))</label>
<details class=\"span-2 step-control-options\" data-step-control-group=\"autodamp\" data-ac-only-field>
<summary>Autodamping &amp; merit-function line search</summary>
<label class=\"check\"><input name=\"power_flow_autodamp\" type=\"hidden\" value=\"false\"><input name=\"power_flow_autodamp\" type=\"checkbox\" value=\"true\" data-autodamp-toggle$(_webui_checked(profile_values, "power_flow_autodamp", _webui_option_default("power_flow_autodamp")))>$(_webui_field_label("power_flow_autodamp", "Autodamping enabled"))</label>
<label>$(_webui_field_label("power_flow_autodamp_min", "Autodamping minimum"))<input name=\"power_flow_autodamp_min\" type=\"number\" step=\"any\" min=\"0\" max=\"1\" data-autodamp-field value=\"$(_webui_input_value(profile_values, "power_flow_autodamp_min", _webui_option_default("power_flow_autodamp_min")))\"></label>
<details class=\"span-2 merit-linesearch-options\">
<summary>Merit-function line search</summary>
<fieldset>
<label class=\"check span-2\"><input name=\"power_flow_merit_enabled\" type=\"hidden\" value=\"false\"><input name=\"power_flow_merit_enabled\" type=\"checkbox\" value=\"true\" data-merit-toggle$(_webui_checked(profile_values, "power_flow_merit_enabled", _webui_option_default("power_flow_merit_enabled")))>$(_webui_field_label("power_flow_merit_enabled", "Enable Armijo merit-function line search (requires Automatic damping)"))</label>
<label>$(_webui_field_label("power_flow_merit_armijo_c1", "Armijo sufficient-decrease constant"))<input name=\"power_flow_merit_armijo_c1\" type=\"number\" step=\"any\" min=\"0\" max=\"0.5\" data-merit-field value=\"$(_webui_input_value(profile_values, "power_flow_merit_armijo_c1", _webui_option_default("power_flow_merit_armijo_c1")))\"></label>
<label class=\"check\"><input name=\"power_flow_merit_fallback_max_mismatch\" type=\"hidden\" value=\"false\"><input name=\"power_flow_merit_fallback_max_mismatch\" type=\"checkbox\" value=\"true\" data-merit-field$(_webui_checked(profile_values, "power_flow_merit_fallback_max_mismatch", _webui_option_default("power_flow_merit_fallback_max_mismatch")))>$(_webui_field_label("power_flow_merit_fallback_max_mismatch", "Fall back to max-mismatch criterion when Armijo is not satisfied"))</label>
<p class=\"field-help span-2\">Residual scaling (<code>scale_p</code>/<code>scale_q</code>/<code>scale_v</code>) is YAML-only and not exposed here.</p>
</fieldset>
</details>
</details>
<details class=\"span-2 step-control-options\" data-step-control-group=\"trust_region\" data-ac-only-field>
<summary>Trust-region step control</summary>
<p class=\"field-help span-2\">Alternative to autodamping (scaled-Newton step control with merit-based step acceptance). Mutually exclusive with autodamping -- enabling one disables the other.</p>
<label class=\"check span-2\"><input name=\"power_flow_trust_region_enabled\" type=\"hidden\" value=\"false\"><input name=\"power_flow_trust_region_enabled\" type=\"checkbox\" value=\"true\" data-trust-region-toggle$(_webui_checked(profile_values, "power_flow_trust_region_enabled", _webui_option_default("power_flow_trust_region_enabled")))>$(_webui_field_label("power_flow_trust_region_enabled", "Enable trust-region step control"))</label>
<label>$(_webui_field_label("power_flow_trust_region_initial_radius", "Initial trust-region radius"))<input name=\"power_flow_trust_region_initial_radius\" type=\"number\" step=\"any\" min=\"0\" max=\"10\" data-trust-region-field value=\"$(_webui_input_value(profile_values, "power_flow_trust_region_initial_radius", _webui_option_default("power_flow_trust_region_initial_radius")))\"></label>
<label>$(_webui_field_label("power_flow_trust_region_eta_accept", "Acceptance ratio (eta)"))<input name=\"power_flow_trust_region_eta_accept\" type=\"number\" step=\"any\" min=\"0\" max=\"1\" data-trust-region-field value=\"$(_webui_input_value(profile_values, "power_flow_trust_region_eta_accept", _webui_option_default("power_flow_trust_region_eta_accept")))\"></label>
<label>$(_webui_field_label("power_flow_trust_region_step_mode", "Step mode"))$(_webui_select("power_flow_trust_region_step_mode", _webui_option_allowed_values("power_flow_trust_region_step_mode"), _webui_selected(profile_values, "power_flow_trust_region_step_mode", _webui_option_default("power_flow_trust_region_step_mode")), "data-trust-region-field"))</label>
<p class=\"field-help span-2\">Radius bounds and shrink/expand factors (<code>min_radius</code>/<code>max_radius</code>/<code>shrink_factor</code>/<code>expand_factor</code>/<code>expand_threshold</code>) are YAML-only and not exposed here.</p>
</details>
<details class=\"span-2 step-control-options\" data-ac-only-field>
<summary>Q-limit handling</summary>
<label class=\"check\"><input name=\"power_flow_qlimits_enabled\" type=\"hidden\" value=\"false\"><input name=\"power_flow_qlimits_enabled\" type=\"checkbox\" value=\"true\"$(_webui_checked(profile_values, "power_flow_qlimits_enabled", _webui_option_default("power_flow_qlimits_enabled")))>$(_webui_field_label("power_flow_qlimits_enabled", "Q-limit handling enabled"))</label>
<label data-nr-only-field>$(_webui_field_label("power_flow_qlimits_enforcement_mode", "Q-limit enforcement mode"))$(_webui_select("power_flow_qlimits_enforcement_mode", vcat(["off"], collect(_webui_option_allowed_values("power_flow_qlimits_enforcement_mode"))), _webui_qlimit_mode_selection(profile_values)))</label>
<p class=\"field-help span-2\">\"off\" switches Q-limit handling off; it is the same thing as clearing the box above, offered here because a mode picked while the handling is off does nothing and looks like it does.</p>
</details>
<details class="span-2 solver-mode-options" open>
<summary>$(_webui_field_label("power_flow_solver", "Solver"))</summary>
<label class="check"><input type="radio" name="power_flow_solver" value="rectangular" data-solver-radio$(_webui_form_string(_webui_selected(profile_values, "power_flow_solver", _webui_option_default("power_flow_solver"))) == "rectangular" ? " checked" : "")>AC (Newton-Raphson, rectangular)</label>
<p class="field-help field-indent">Starts by default from angles produced by a fast DC pre-pass (see <strong>Start angle mode</strong> below, default <code>dc</code>) and then iterates with Newton-Raphson -- this is not a standalone DC solution.</p>
<label class="check"><input type="radio" name="power_flow_solver" value="apslf" data-solver-radio$(_webui_form_string(_webui_selected(profile_values, "power_flow_solver", _webui_option_default("power_flow_solver"))) == "apslf" ? " checked" : "")>APSLF (AnalyticLoadFlow)</label>
<label class="check"><input type="radio" name="power_flow_solver" value="dc" data-solver-radio$(_webui_form_string(_webui_selected(profile_values, "power_flow_solver", _webui_option_default("power_flow_solver"))) == "dc" ? " checked" : "")>DC (linear screening model, replaces Newton-Raphson entirely)</label>
</details>
<fieldset id="apslf-solver-options" class="span-2 apslf-solver-options" data-apslf-solver-options hidden>
<legend>APSLF solver options</legend>
<label class=\"field-indent\">$(_webui_field_label("power_flow_apslf_order", "Highest coefficient (order)"))<input name=\"power_flow_apslf_order\" type=\"number\" min=\"1\" value=\"$(_webui_input_value(profile_values, "power_flow_apslf_order", _webui_option_default("power_flow_apslf_order")))\"></label>
<label class=\"check field-indent\"><input name=\"power_flow_apslf_use_pade\" type=\"hidden\" value=\"false\"><input name=\"power_flow_apslf_use_pade\" type=\"checkbox\" value=\"true\"$(_webui_checked(profile_values, "power_flow_apslf_use_pade", _webui_option_default("power_flow_apslf_use_pade")))>$(_webui_field_label("power_flow_apslf_use_pade", "Padé evaluation"))</label>
<label class=\"check field-indent\"><input name=\"power_flow_apslf_nr_polish\" type=\"hidden\" value=\"false\"><input name=\"power_flow_apslf_nr_polish\" type=\"checkbox\" value=\"true\"$(_webui_checked(profile_values, "power_flow_apslf_nr_polish", _webui_option_default("power_flow_apslf_nr_polish")))>$(_webui_field_label("power_flow_apslf_nr_polish", "NR polish"))</label>
</fieldset>
<details id="apslf-start-options" class="span-2 apslf-start-options" data-apslf-start-options data-ac-only-field>
<summary>Newton-Raphson start values</summary>
<label class=\"check\"><input name=\"power_flow_apslf_start_enabled\" type=\"hidden\" value=\"false\"><input name=\"power_flow_apslf_start_enabled\" type=\"checkbox\" value=\"true\" data-apslf-start-toggle$(_webui_checked(profile_values, "power_flow_apslf_start_enabled", _webui_option_default("power_flow_apslf_start_enabled")))>$(_webui_field_label("power_flow_apslf_start_enabled", "Use APSLF start values"))</label>
<label class=\"field-indent\">$(_webui_field_label("power_flow_apslf_start_order", "Highest coefficient (order)"))<input name=\"power_flow_apslf_start_order\" type=\"number\" min=\"1\" data-apslf-start-order value=\"$(_webui_input_value(profile_values, "power_flow_apslf_start_order", _webui_option_default("power_flow_apslf_start_order")))\"></label>
<label class=\"check\"><input name=\"power_flow_dc_seed_unconditional\" type=\"hidden\" value=\"false\"><input name=\"power_flow_dc_seed_unconditional\" type=\"checkbox\" value=\"true\" data-dc-seed-toggle$(_webui_checked(profile_values, "power_flow_dc_seed_unconditional", _webui_option_default("power_flow_dc_seed_unconditional")))>$(_webui_field_label("power_flow_dc_seed_unconditional", "Use DC start values"))</label>
</details>
<label data-nr-only-field>$(_webui_field_label("power_flow_wrong_branch_detection", "Wrong-branch detection"))$(_webui_select("power_flow_wrong_branch_detection", _webui_option_allowed_values("power_flow_wrong_branch_detection"), _webui_selected(profile_values, "power_flow_wrong_branch_detection", _webui_option_default("power_flow_wrong_branch_detection"))))</label>
<label data-nr-only-field data-dc-seed-inactive-field>$(_webui_field_label("power_flow_start_angle_mode", "Start angle mode"))$(_webui_select("power_flow_start_angle_mode", _webui_option_allowed_values("power_flow_start_angle_mode"), _webui_selected(profile_values, "power_flow_start_angle_mode", _webui_option_default("power_flow_start_angle_mode"))))</label>
<label data-nr-only-field data-dc-seed-inactive-field>$(_webui_field_label("power_flow_start_voltage_mode", "Start voltage mode"))$(_webui_select("power_flow_start_voltage_mode", _webui_option_allowed_values("power_flow_start_voltage_mode"), _webui_selected(profile_values, "power_flow_start_voltage_mode", _webui_option_default("power_flow_start_voltage_mode"))))</label>
<label>$(_webui_field_label("output_logfile_results", "Logfile output mode"))$(_webui_select("output_logfile_results", _webui_option_allowed_values("output_logfile_results"), _webui_selected(profile_values, "output_logfile_results", _webui_option_default("output_logfile_results"))))</label>
</fieldset>
<fieldset class=\"solver-backend-options\" data-nr-only-field>
<legend>Solver backend</legend>
<label>$(_webui_field_label("power_flow_linear_solver", "Linear solver backend"))$(_webui_select("power_flow_linear_solver", _webui_option_allowed_values("power_flow_linear_solver"), _webui_selected(profile_values, "power_flow_linear_solver", _webui_option_default("power_flow_linear_solver"))))</label>
<p class=\"field-help\">Sparse linear-algebra backend for the rectangular Newton step only (independent of the <strong>Solver</strong> choice above). <code>umfpack_reuse</code> reuses the symbolic analysis across iterations via <code>lu!</code>; <code>umfpack</code> is the default behavior.</p>
</fieldset>
<fieldset class=\"external-grid-options\">
<legend>External grid source</legend>
<label class=\"check span-2\"><input name=\"power_flow_external_grid_enabled\" type=\"hidden\" value=\"false\"><input name=\"power_flow_external_grid_enabled\" type=\"checkbox\" value=\"true\"$(_webui_checked(profile_values, "power_flow_external_grid_enabled", _webui_option_default("power_flow_external_grid_enabled")))>$(_webui_field_label("power_flow_external_grid_enabled", "Compute the marked slack as a non-ideal external-grid source"))</label>
<label>$(_webui_field_label("power_flow_external_grid_source", "Sk''/R-X source"))$(_webui_select("power_flow_external_grid_source", _webui_option_allowed_values("power_flow_external_grid_source"), _webui_selected(profile_values, "power_flow_external_grid_source", _webui_option_default("power_flow_external_grid_source"))))</label>
<label>$(_webui_field_label("power_flow_external_grid_sk_mva", "Short-circuit power Sk'' [MVA]"))<input name=\"power_flow_external_grid_sk_mva\" type=\"number\" step=\"any\" min=\"0\" value=\"$(_webui_input_value(profile_values, "power_flow_external_grid_sk_mva", _webui_option_default("power_flow_external_grid_sk_mva")))\"></label>
<label>$(_webui_field_label("power_flow_external_grid_rx", "R/X ratio"))<input name=\"power_flow_external_grid_rx\" type=\"number\" step=\"any\" min=\"0\" value=\"$(_webui_input_value(profile_values, "power_flow_external_grid_rx", _webui_option_default("power_flow_external_grid_rx")))\"></label>
<p class=\"field-help\">The reference voltage moves to a hidden internal bus behind the feeder impedance <code>z = Un²/Sk''</code> — the former slack bus becomes an ordinary bus whose voltage droops under load. <code>auto</code> prefers the Sk''/R-X values a CGMES delivery declares on the slack bus's <code>ExternalNetworkInjection</code> and falls back to the numbers above (MATPOWER/DTF cases carry no such data); <code>config</code> always uses the numbers above. Mutually exclusive with the distributed slack — both decide who covers the imbalance, and combined the source's import would be forced to zero.</p>
</fieldset>
$(isempty(profile_path) ? "" : "<fieldset class=\"saved-case-settings\">
<legend>Saved settings for this case</legend>
<p class=\"field-help\">This case has stored Web UI settings (<code>$(_webui_escape(basename(profile_path)))</code>). They prefill the form and outrank the configuration file for the keys they contain — including ones you may not expect, such as a stored solver choice. Resetting deletes the stored settings only; the case file itself is kept.</p>
<div class=\"actions\"><button type=\"submit\" class=\"secondary-button\" formaction=\"/powerflow/case-settings/reset\" formmethod=\"post\" formnovalidate>Reset saved settings for this case</button></div>
</fieldset>")
<fieldset class=\"startup-options\">
<legend>Startup</legend>
<p class=\"field-help\">Hidden compile runs at startup so the first real run — and the first <strong>Short circuit</strong> click — do not pay the compilation. Costs a few seconds once. Saved in the configuration file; takes effect at the <em>next</em> Web UI start.</p>
</fieldset>
<fieldset class=\"non-convergence-options\" data-nr-only-field>
<legend>Non-convergence handling</legend>
<label class=\"check span-2\"><input name=\"power_flow_rescue\" type=\"hidden\" value=\"false\"><input name=\"power_flow_rescue\" type=\"checkbox\" value=\"true\"$(_webui_checked(profile_values, "power_flow_rescue", _webui_option_default("power_flow_rescue")))>$(_webui_field_label("power_flow_rescue", "Rescue: retry a failed AC solve (alternate start, autodamp, DC-seeded start)"))</label>
<label class=\"check span-2\"><input name=\"runtime_parallel_enabled\" type=\"hidden\" value=\"false\"><input name=\"runtime_parallel_enabled\" type=\"checkbox\" value=\"true\"$(_webui_checked(profile_values, "runtime_parallel_enabled", _webui_option_default("runtime_parallel_enabled")))>$(_webui_field_label("runtime_parallel_enabled", "Parallel execution: use Julia threads for independent work items (islands, sweeps, batches)"))</label>
<label class=\"check span-2\"><input name=\"power_flow_dc_fallback\" type=\"hidden\" value=\"false\"><input name=\"power_flow_dc_fallback\" type=\"checkbox\" value=\"true\"$(_webui_checked(profile_values, "power_flow_dc_fallback", _webui_option_default("power_flow_dc_fallback")))>$(_webui_field_label("power_flow_dc_fallback", "DC fallback: keep a standalone DC result when AC (and rescue) fail"))</label>
<p class=\"field-help\">The rescue ladder restarts from the original start state and logs the winning strategy. The DC fallback leaves angles and branch P flows (vm = 1 pu); the AC status honestly stays non-converged.</p>
</fieldset>
<fieldset class=\"start-current-iteration-options advanced-start-values\" data-nr-only-field>
<legend>Advanced start values</legend>
<label class=\"check span-2\"><input name=\"power_flow_start_current_iteration_enabled\" type=\"hidden\" value=\"false\"><input name=\"power_flow_start_current_iteration_enabled\" type=\"checkbox\" value=\"true\"$(_webui_checked(profile_values, "power_flow_start_current_iteration_enabled", _webui_option_default("power_flow_start_current_iteration_enabled")))>$(_webui_field_label("power_flow_start_current_iteration_enabled", "Enable current-iteration pre-solve"))</label>
<label>$(_webui_field_label("power_flow_start_current_iteration_max_iter", "Current-iteration max iterations"))<input name=\"power_flow_start_current_iteration_max_iter\" type=\"number\" min=\"1\" value=\"$(_webui_input_value(profile_values, "power_flow_start_current_iteration_max_iter", _webui_option_default("power_flow_start_current_iteration_max_iter")))\"></label>
<label>$(_webui_field_label("power_flow_start_current_iteration_tol", "Current-iteration tolerance"))<input name=\"power_flow_start_current_iteration_tol\" type=\"number\" step=\"any\" min=\"0\" value=\"$(_webui_input_value(profile_values, "power_flow_start_current_iteration_tol", _webui_option_default("power_flow_start_current_iteration_tol")))\"></label>
<label>$(_webui_field_label("power_flow_start_current_iteration_damping", "Current-iteration damping"))<input name=\"power_flow_start_current_iteration_damping\" type=\"number\" step=\"any\" min=\"0\" max=\"1\" value=\"$(_webui_input_value(profile_values, "power_flow_start_current_iteration_damping", _webui_option_default("power_flow_start_current_iteration_damping")))\"></label>
<label class=\"check\"><input name=\"power_flow_start_current_iteration_accept_only_if_improved\" type=\"hidden\" value=\"false\"><input name=\"power_flow_start_current_iteration_accept_only_if_improved\" type=\"checkbox\" value=\"true\"$(_webui_checked(profile_values, "power_flow_start_current_iteration_accept_only_if_improved", _webui_option_default("power_flow_start_current_iteration_accept_only_if_improved")))>$(_webui_field_label("power_flow_start_current_iteration_accept_only_if_improved", "Accept only if improved"))</label>
<label>$(_webui_field_label("power_flow_start_current_iteration_min_improvement_factor", "Minimum improvement factor"))<input name=\"power_flow_start_current_iteration_min_improvement_factor\" type=\"number\" step=\"any\" min=\"0\" value=\"$(_webui_input_value(profile_values, "power_flow_start_current_iteration_min_improvement_factor", _webui_option_default("power_flow_start_current_iteration_min_improvement_factor")))\"></label>
<label>$(_webui_field_label("power_flow_start_current_iteration_vm_min_pu", "Minimum voltage guard [pu]"))<input name=\"power_flow_start_current_iteration_vm_min_pu\" type=\"number\" step=\"any\" min=\"0\" value=\"$(_webui_input_value(profile_values, "power_flow_start_current_iteration_vm_min_pu", _webui_option_default("power_flow_start_current_iteration_vm_min_pu")))\"></label>
<label>$(_webui_field_label("power_flow_start_current_iteration_vm_max_pu", "Maximum voltage guard [pu]"))<input name=\"power_flow_start_current_iteration_vm_max_pu\" type=\"number\" step=\"any\" min=\"0\" value=\"$(_webui_input_value(profile_values, "power_flow_start_current_iteration_vm_max_pu", _webui_option_default("power_flow_start_current_iteration_vm_max_pu")))\"></label>
<label>$(_webui_field_label("power_flow_start_current_iteration_max_angle_step_deg", "Maximum angle-step guard [deg]"))<input name=\"power_flow_start_current_iteration_max_angle_step_deg\" type=\"number\" step=\"any\" min=\"0\" value=\"$(_webui_input_value(profile_values, "power_flow_start_current_iteration_max_angle_step_deg", _webui_option_default("power_flow_start_current_iteration_max_angle_step_deg")))\"></label>
<label class=\"check\"><input name=\"power_flow_start_current_iteration_only_for_large_cases\" type=\"hidden\" value=\"false\"><input name=\"power_flow_start_current_iteration_only_for_large_cases\" type=\"checkbox\" value=\"true\"$(_webui_checked(profile_values, "power_flow_start_current_iteration_only_for_large_cases", _webui_option_default("power_flow_start_current_iteration_only_for_large_cases")))>$(_webui_field_label("power_flow_start_current_iteration_only_for_large_cases", "Only for large cases"))</label>
</fieldset>
<fieldset class=\"benchmark-section\">
<legend>Benchmark / repeated timing</legend>
<label>$(_webui_field_label("benchmark_samples", "Benchmark samples (max. repeated measurements)"))<input name=\"benchmark_samples\" type=\"number\" min=\"1\" value=\"$(_webui_input_value(profile_values, "benchmark_samples", _webui_option_default("benchmark_samples")))\"></label>
<label>$(_webui_field_label("benchmark_seconds", "Benchmark max. time budget [s]"))<input name=\"benchmark_seconds\" type=\"number\" step=\"any\" min=\"0\" value=\"$(_webui_input_value(profile_values, "benchmark_seconds", _webui_option_default("benchmark_seconds")))\"></label>
</fieldset>
<label class=\"check span-2\"><input name=\"ignore_webui_settings\" type=\"hidden\" value=\"false\"><input name=\"ignore_webui_settings\" type=\"checkbox\" value=\"true\">$(_webui_field_label("ignore_webui_settings", "Ignore Web UI settings and use configuration defaults"))</label>
</details>
"""
end

# the settings-page script: the CSV language default, the solver and
# step-control toggle machinery, and the tolerance spinner, moved verbatim
# from the run monolith (their controls all render here now)
function _webui_settings_script()::String
  return """
<script>
document.addEventListener('DOMContentLoaded', function () {
  const csvFormat = document.querySelector('select[name="detailed_result_csv_format"]');
  if (csvFormat !== null && csvFormat.value === 'excel_us') {
    const languages = navigator.languages && navigator.languages.length > 0 ? navigator.languages : [navigator.language || ''];
    const defaultFormat = languages.some(function (language) {
      return String(language).toLowerCase().startsWith('de');
    }) ? 'excel_de' : 'excel_us';
    csvFormat.value = defaultFormat;
  }
  const solverRadios = document.querySelectorAll('input[data-solver-radio]');
  const apslfSolverOptions = document.querySelector('[data-apslf-solver-options]');
  const apslfStartOptions = document.querySelector('[data-apslf-start-options]');
  const getSolverMode = function () {
    let value = 'rectangular';
    solverRadios.forEach(function (radio) { if (radio.checked) value = radio.value; });
    return value;
  };
  const isDcMode = function () { return getSolverMode() === 'dc'; };
  const isApslfMode = function () { return getSolverMode() === 'apslf'; };
  // Gray out (disable, but keep visible/in place) a field group that does not apply
  // to the currently selected solver, instead of hiding it: mutually exclusive
  // solver options stay where the user last saw them rather than jumping around.
  const setSolverGroupInactive = function (container, inactive) {
    if (container === null) return;
    container.classList.toggle('disabled', inactive);
    const controls = container.matches('input, select') ? [container] : container.querySelectorAll('input, select');
    controls.forEach(function (control) { control.disabled = inactive; });
  };
  const updateSolverOptions = function () {
    const dc = isDcMode();
    const apslf = isApslfMode();
    setSolverGroupInactive(apslfSolverOptions, dc || !apslf);
    setSolverGroupInactive(apslfStartOptions, dc || apslf);
  };
  const apslfStartToggle = document.querySelector('input[data-apslf-start-toggle]');
  const apslfStartOrderInput = document.querySelector('input[data-apslf-start-order]');
  const dcSeedToggle = document.querySelector('input[data-dc-seed-toggle]');
  const updateApslfStartOrder = function () {
    if (apslfStartOrderInput !== null) apslfStartOrderInput.disabled = apslfStartToggle !== null && !apslfStartToggle.checked;
  };
  // "Use APSLF start values" and "Use DC start values" are two mutually exclusive
  // start-value sources for the same rectangular NR solve (the underlying
  // configuration rejects setting both at once) -- checking one unchecks the other,
  // mirroring the existing autodamp/trust-region exclusivity pattern below.
  const updateStartValueSource = function (changedToggle) {
    if (changedToggle === 'dc_seed' && dcSeedToggle !== null && dcSeedToggle.checked && apslfStartToggle !== null && apslfStartToggle.checked) {
      apslfStartToggle.checked = false;
    } else if (changedToggle === 'apslf_start' && apslfStartToggle !== null && apslfStartToggle.checked && dcSeedToggle !== null && dcSeedToggle.checked) {
      dcSeedToggle.checked = false;
    }
    updateApslfStartOrder();
    updateStepControlOptions();
  };
  if (apslfStartToggle !== null) {
    updateApslfStartOrder();
    apslfStartToggle.addEventListener('change', function () { updateStartValueSource('apslf_start'); });
  }
  if (dcSeedToggle !== null) {
    dcSeedToggle.addEventListener('change', function () { updateStartValueSource('dc_seed'); });
  }
  const autodampToggle = document.querySelector('input[data-autodamp-toggle]');
  const trustRegionToggle = document.querySelector('input[data-trust-region-toggle]');
  const meritToggle = document.querySelector('input[data-merit-toggle]');
  const autodampFields = document.querySelectorAll('[data-autodamp-field]');
  const meritFields = document.querySelectorAll('[data-merit-field]');
  const trustRegionFields = document.querySelectorAll('[data-trust-region-field]');
  const autodampGroup = document.querySelector('[data-step-control-group="autodamp"]');
  const trustRegionGroup = document.querySelector('[data-step-control-group="trust_region"]');
  const nrOnlyFields = document.querySelectorAll('[data-nr-only-field]');
  let updatingStepControl = false;
  const updateStepControlOptions = function (changedToggle) {
    if (updatingStepControl) return;
    updatingStepControl = true;
    if (changedToggle === 'trust_region' && trustRegionToggle !== null && trustRegionToggle.checked && autodampToggle !== null && autodampToggle.checked) {
      autodampToggle.checked = false;
    } else if (changedToggle === 'autodamp' && autodampToggle !== null && autodampToggle.checked && trustRegionToggle !== null && trustRegionToggle.checked) {
      trustRegionToggle.checked = false;
    }
    const dc = isDcMode();
    const apslf = isApslfMode();
    const hideNrOnly = apslf || dc;
    const autodampOn = !hideNrOnly && autodampToggle !== null && autodampToggle.checked;
    const trustRegionOn = !hideNrOnly && trustRegionToggle !== null && trustRegionToggle.checked;
    autodampFields.forEach(function (field) { field.disabled = !autodampOn; });
    if (meritToggle !== null) {
      meritToggle.disabled = !autodampOn;
      if (!autodampOn && meritToggle.checked) meritToggle.checked = false;
    }
    const meritOn = autodampOn && meritToggle !== null && meritToggle.checked;
    meritFields.forEach(function (field) { field.disabled = !meritOn; });
    trustRegionFields.forEach(function (field) { field.disabled = !trustRegionOn; });
    if (autodampGroup !== null) {
      autodampGroup.classList.toggle('disabled', !autodampOn);
    }
    if (trustRegionGroup !== null) {
      trustRegionGroup.classList.toggle('disabled', !trustRegionOn);
    }
    const dcSeedActive = dcSeedToggle !== null && dcSeedToggle.checked;
    nrOnlyFields.forEach(function (container) {
      const dcSeedMakesInactive = dcSeedActive && container.hasAttribute('data-dc-seed-inactive-field');
      setSolverGroupInactive(container, hideNrOnly || dcSeedMakesInactive);
    });
    updatingStepControl = false;
  };
  if (autodampToggle !== null) {
    updateStepControlOptions();
    autodampToggle.addEventListener('change', function () { updateStepControlOptions('autodamp'); });
  }
  if (trustRegionToggle !== null) {
    trustRegionToggle.addEventListener('change', function () { updateStepControlOptions('trust_region'); });
  }
  if (meritToggle !== null) {
    meritToggle.addEventListener('change', function () { updateStepControlOptions('merit'); });
  }
  const acOnlyFields = document.querySelectorAll('[data-ac-only-field]');
  const updateSolverMode = function () {
    const dc = isDcMode();
    acOnlyFields.forEach(function (container) { setSolverGroupInactive(container, dc); });
    updateSolverOptions();
    updateStepControlOptions();
  };
  if (solverRadios.length > 0) {
    updateSolverMode();
    solverRadios.forEach(function (radio) { radio.addEventListener('change', updateSolverMode); });
  }
  const toleranceInput = document.querySelector('input[name="power_flow_tol"][data-tolerance-input]');
  const parseToleranceParts = function (valueText) {
    const trimmed = String(valueText).trim();
    if (trimmed === '') return null;
    const lower = trimmed.toLowerCase();
    const eIndex = lower.indexOf('e');
    const mantissaText = eIndex === -1 ? trimmed : trimmed.slice(0, eIndex);
    const exponentText = eIndex === -1 ? '0' : trimmed.slice(eIndex + 1);
    const rawMantissa = Number(mantissaText);
    const rawExponent = Number(exponentText);
    if (!Number.isFinite(rawMantissa) || rawMantissa <= 0 || !Number.isFinite(rawExponent)) return null;
    if (eIndex === -1) {
      const exponent = Math.floor(Math.log10(rawMantissa));
      return {mantissa: rawMantissa / Math.pow(10, exponent), exponent: exponent};
    }
    return {mantissa: rawMantissa, exponent: rawExponent};
  };
  const formatToleranceParts = function (mantissa, exponent) {
    const rounded = Number(mantissa.toPrecision(12));
    return String(rounded) + 'e' + String(exponent);
  };
  const stepTolerance = function (direction) {
    if (toleranceInput === null) return;
    const parts = parseToleranceParts(toleranceInput.value) || {mantissa: 1, exponent: -8};
    toleranceInput.value = formatToleranceParts(parts.mantissa, parts.exponent + direction);
    toleranceInput.dispatchEvent(new Event('change', {bubbles: true}));
  };
  if (toleranceInput !== null) {
    toleranceInput.addEventListener('keydown', function (event) {
      if (event.key === 'ArrowUp') {
        event.preventDefault();
        stepTolerance(1);
      } else if (event.key === 'ArrowDown') {
        event.preventDefault();
        stepTolerance(-1);
      }
    });
    document.querySelectorAll('.tolerance-spin button[data-tolerance-direction]').forEach(function (button) {
      button.addEventListener('click', function () {
        stepTolerance(button.dataset.toleranceDirection === 'up' ? 1 : -1);
        toleranceInput.focus();
      });
    });
  }
});
</script>"""
end


"""
    render_settings_page(; kwargs...) -> String

The Settings page (stage 4A block 3): solver, start, output, and expert
options plus the configuration block (config file display, maintenance
actions, saved case settings, the ignore switch) on ONE page. Values
prefill from the shared resolution chain; saving writes them to the
chosen target (the selected case's configuration file, or the general
YAML), and runs pick them up through `resolve_config` instead of the run
POST.
"""
function render_settings_page(;
  output_root::AbstractString = "results/powerflow_service",
  case_directory::Union{Nothing,AbstractString} = nothing,
  operation_log::AbstractString = webui_operation_log_path(output_root),
  error_message = nothing,
  application_root::AbstractString = _webui_application_root(),
  selected_casefile::AbstractString = "",
  selected_config_file::AbstractString = "",
  save_message::AbstractString = "",
  case_profile = nothing,
  submitted_form = nothing,
  show_case_settings_notice::Bool = true,
)::String
  ctx = _webui_case_context(; application_root, case_directory, selected_casefile, selected_config_file, case_profile, submitted_form, show_case_settings_notice)
  profile_values = ctx.profile_values
  profile_path = ctx.profile_path
  error_html = _webui_error_alert_html(error_message)
  save_html = isempty(strip(save_message)) ? "" : "<div class=\"alert info settings-save-result\" role=\"status\">$(_webui_escape(save_message))</div>"
  config_default = isempty(selected_config_file) ? DEFAULT_SPARLECTRA_CONFIG_PATH : selected_config_file
  info_menu = _webui_powerflow_info_menu(; output_root, config_file = config_default, case_directory = ctx.effective_case_directory, operation_log)
  case_display = isempty(strip(ctx.effective_case_value)) ? "no case selected" : _webui_escape(ctx.effective_case_value)
  target_case_attrs = isempty(strip(ctx.effective_case_value)) ? " disabled" : " checked"
  target_general_attrs = isempty(strip(ctx.effective_case_value)) ? " checked" : ""
  sections = _webui_settings_sections_html(; profile_values, config_default, profile_path, selected_casefile = String(ctx.effective_case_value), selected_config_file)
  content = """
$(_webui_feedback_modal_html([error_html, save_html, ctx.profile_notice, ctx.case_file_notice]))<p class=\"lede\">Solver, output, and expert options. Values prefill from the effective configuration for <code>$(case_display)</code> (<a href=\"/powerflow/case$(isempty(strip(ctx.effective_case_value)) ? "" : "?casefile=" * _webui_urlencode(ctx.effective_case_value))\">change on the Case page</a>); runs read them through the configuration precedence.</p>
<form id=\"settings-form\" method=\"post\" action=\"/powerflow/settings/save\" class=\"panel form-grid settings-form-card\">
<input type=\"hidden\" name=\"casefile\" value=\"$(_webui_escape(ctx.effective_case_value))\">
<input type=\"hidden\" name=\"config_file\" value=\"$(_webui_escape(config_default))\">
<fieldset class=\"settings-target\">
<legend>Apply to</legend>
<label class=\"check\"><input type=\"radio\" name=\"settings_target\" value=\"this_case\"$(target_case_attrs)>this case (<code>$(case_display)</code>)</label>
<small class=\"settings-target-hint\">case-scope keys go into the case configuration file; machine-scope keys (output, benchmark, runtime, webui) stay out</small>
<label class=\"check\"><input type=\"radio\" name=\"settings_target\" value=\"general\"$(target_general_attrs)>configuration file</label>
<small class=\"settings-target-hint\">all keys, for every case on this machine: <code>$(_webui_escape(config_default))</code></small>
</fieldset>
$(sections)
<div class=\"span-2 actions\"><button class=\"secondary-button\" type=\"submit\">Save settings</button></div>
</form>
$(_webui_settings_script())
$(_WEBUI_FEEDBACK_MODAL_SCRIPT)
$(_WEBUI_INFO_MENU_SCRIPT)"""
  return _webui_layout("Settings", content; header_info = info_menu)
end

function render_powerflow_form(;
  output_root::AbstractString = "results/powerflow_service",
  case_directory::Union{Nothing,AbstractString} = nothing,
  operation_log::AbstractString = webui_operation_log_path(output_root),
  error_message = nothing,
  application_root::AbstractString = _webui_application_root(),
  selected_casefile::AbstractString = "",
  selected_config_file::AbstractString = "",
  active_run = get_active_webui_powerflow_job(),
  config_notice = nothing,
  case_profile = nothing,
  submitted_form = nothing,
  import_message::AbstractString = "",
  download_file::AbstractString = "",
  show_case_settings_notice::Bool = true,
  se_query::AbstractDict = Dict{String,String}(),
)::String
  ctx = _webui_case_context(; application_root, case_directory, selected_casefile, selected_config_file, case_profile, submitted_form, show_case_settings_notice)
  profile_values = ctx.profile_values
  error_html = _webui_error_alert_html(error_message)
  # a freshly exported file is offered right here: it was written next to the
  # case (where a run finds it), and this is the way back out to the browser
  download_html = isempty(strip(download_file)) ? "" :
                  " <a class=\"case-download-link\" href=\"/powerflow/case/download?case=$(_webui_urlencode(String(download_file)))\">Download $(_webui_escape(String(download_file)))</a>"
  import_html = isempty(strip(import_message)) && isempty(download_html) ? "" : "<div class=\"alert info case-import-result\" role=\"status\">$(_webui_escape(import_message))$(download_html)</div>"
  effective_case_directory = ctx.effective_case_directory
  for002_candidates = ctx.for002_candidates
  effective_case_value = ctx.effective_case_value
  dat_case_assistance = ctx.dat_case_assistance
  sc_state = ctx.sc_state
  sc_disabled_attr = sc_state == "ready" ? "" : " disabled"
  sc_title = sc_state == "ready" ? "Balanced short-circuit currents (IEC 60909-0): Ik'' max/min per bus from the delivery's harvested short-circuit data — no power-flow solve involved." :
    sc_state == "missing-data" ? "This case carries no usable short-circuit source data (no machines, feeder short-circuit currents, or equivalent impedances)." :
    "Short-circuit evaluation needs a CGMES delivery with harvested short-circuit data, or a case file carrying sc_source entries."
  # scenario task step 6 (maintainer revision 2026-09-03): the scenario
  # editor and the file_block source are SCF-only in this version; every
  # other format keeps the n1_* sources and gets the export hint instead
  scen_is_scf = ctx.scen_is_scf
  scen_fileblock_option = scen_is_scf ? "<option value=\"file_block\">case file scenarios</option>" : ""
  scen_editor_html = scen_is_scf ?
    "<a class=\"weights-link scenario-editor-link\" href=\"/powerflow/scenarios?case=$(_webui_urlencode(selected_casefile))\" title=\"Edit the case file's scenarios block (list plus form; saved into the SCF file)\">edit scenarios</a>" :
    "<span class=\"field-hint scenario-scf-hint\">Scenarios need an SCF case; export this case as SCF first (Case export above).</span>"
  # the DTF outage details open automatically for a .DAT case, mirroring the
  # Case page's format assistance (the fields inside are RUN parameters and
  # therefore stayed on this page in stage 4A)
  dtf_details_attrs = dat_case_assistance ? " class=\"span-2 dtf-internal-section is-dat-selected\" open" : " class=\"span-2 dtf-internal-section\""
  for002_reference_value = submitted_form isa AbstractDict ? strip(_webui_form_string(_webui_form_value(submitted_form, "for002_reference_file", ""))) : ""
  for002_list_options = join(("<option value=\"$(_webui_escape(candidate))\">$(_webui_escape(candidate))</option>" for candidate in for002_candidates), "")
  for002_list_html = isempty(for002_candidates) ? "" : "<datalist id=\"for002-reference-candidates\">$(for002_list_options)</datalist>"
  for002_list_attr = isempty(for002_candidates) ? "" : " list=\"for002-reference-candidates\""
  config_default = isempty(selected_config_file) ? DEFAULT_SPARLECTRA_CONFIG_PATH : selected_config_file
  config_control = "<input type=\"hidden\" name=\"config_file\" value=\"$(_webui_escape(config_default))\">"
  info_menu = _webui_powerflow_info_menu(; output_root, config_file = config_default, case_directory = effective_case_directory, operation_log)
  notice_html = _webui_config_notice_html(config_notice)
  form = """
$(_webui_feedback_modal_html([error_html, import_html]))$(_webui_active_run_banner(active_run))$(notice_html)<p class=\"lede\">Run a local grid case: power flow, N-1, short circuit, or state estimation.</p>
<form id=\"powerflow-run-form\" data-powerflow-form method=\"post\" action=\"/powerflow/run\" class=\"panel form-grid powerflow-form-card\">
$(config_control)
<label class="span-2">$(_webui_field_label("casefile", "Case file"))<input type="hidden" name="casefile" value="$(_webui_escape(effective_case_value))"><span class="selected-case-display"><code>$(isempty(strip(effective_case_value)) ? "no case selected" : _webui_escape(effective_case_value))</code> <a class="case-page-link" href="/powerflow/case$(isempty(strip(effective_case_value)) ? "" : "?casefile=" * _webui_urlencode(effective_case_value))">change on the Case page</a></span><small class="field-hint">Choosing, importing, exporting, and configuring cases moved to the Case page; this run uses the case named here with its saved case options.</small></label>
<details class=\"span-2 expert-section runs-expert-section\"$(dat_case_assistance ? " open" : "")>
<summary>Advanced run options</summary>
<details$(dtf_details_attrs)>
<summary>DTF outage run (internal diagnostics)</summary>
<fieldset>
<label>$(_webui_field_label("for002_reference_file", "Optional FOR002 reference file"))<input name="for002_reference_file" value=\"$(_webui_escape(for002_reference_value))\" placeholder="examples/FOR002.DAT"$for002_list_attr>$(for002_list_html)</label>
<label><span class="field-label">DTF outage run mode</span><select name="dtf_outage_selection_mode"><option value="none">Run base case only</option><option value="all">Run all DTF outage records</option><option value="selected">Run selected DTF outage records</option></select></label>
<label>$(_webui_field_label("dtf_outage_selection", "Selected DTF outage labels/indices"))<input name="dtf_outage_selection" placeholder="1 or L1 ALPHA S1 -> BETA1 S1"></label>
<label class="check"><input name="write_outage_artifacts" type="hidden" value="false"><input name="write_outage_artifacts" type="checkbox" value="true" checked>Write DTF outage artifacts</label>
<label class="check"><input name="matpower_export_requested" type="hidden" value="false"><input name="matpower_export_requested" type="checkbox" value="true">Write MATPOWER export artifact</label>
<label class="check"><input name="write_outage_matpower_exports" type="hidden" value="false"><input name="write_outage_matpower_exports" type="checkbox" value="true">Write MATPOWER outage exports</label>
</fieldset>
</details>
</details>
<label class=\"check span-2 benchmark-trigger\"><input name=\"benchmark_enabled\" type=\"hidden\" value=\"false\"><input name=\"benchmark_enabled\" type=\"checkbox\" value=\"true\"$(_webui_checked(profile_values, "benchmark_enabled", _webui_option_default("benchmark_enabled")))>$(_webui_field_label("benchmark_enabled", "Enable benchmark measurements"))</label>
<div class=\"span-2 actions\"><button class=\"powerflow-submit\" type=\"submit\"><span class=\"submit-spinner\" aria-hidden=\"true\"></span><span class=\"submit-label\">Start PowerFlow run</span><span class=\"submit-progress-label\" role=\"status\" aria-live=\"polite\">Running PowerFlow…</span></button><button class=\"powerflow-submit diagnose-submit\" type=\"submit\" name=\"diagnose_mode\" value=\"true\" title=\"Run this case in diagnostic mode: evaluates the mismatch at the case's own stored VM/VA (no corrective Newton step) and writes a diagnostic report to diagnose.log.\"><span class=\"submit-spinner\" aria-hidden=\"true\"></span><span class=\"submit-label\">Diagnose</span><span class=\"submit-progress-label\" role=\"status\" aria-live=\"polite\">Running diagnosis…</span></button><button class=\"powerflow-submit short-circuit-submit\" type=\"submit\" name=\"short_circuit_mode\" value=\"true\" data-short-circuit-button data-sc-state=\"$(sc_state)\"$(sc_disabled_attr) title=\"$(_webui_escape(sc_title))\"><span class=\"submit-spinner\" aria-hidden=\"true\"></span><span class=\"submit-label\">Short circuit</span><span class=\"submit-progress-label\" role=\"status\" aria-live=\"polite\">Computing short circuit…</span></button><label class=\"contingency-kind\">N-1 kind <select name=\"contingency_kind\" data-contingency-kind><option value=\"branch\">branch</option><option value=\"gen\">generator</option></select></label><label class=\"scenario-source\" title=\"Where the outage or scenario list comes from: generated N-1 lists, or the SCF case file's own scenarios block.\">Scenarios <select name=\"scenario_source\" data-scenario-source><option value=\"\">generated (N-1 kind)</option><option value=\"n1_all\">N-1 all</option><option value=\"n1_branches\">N-1 branches</option><option value=\"n1_generators\">N-1 generators</option>$(scen_fileblock_option)</select></label><label class=\"screening-mode\" title=\"Contingency screening on the base factorization: off runs every scenario fully (default), flag estimates first and fully solves only flagged scenarios, only reports estimates. Opt in after checking share and margins on your network (see the docs).\">Screening <select name=\"screening_mode\" data-screening-mode><option value=\"\">configured</option><option value=\"off\">off</option><option value=\"flag\">flag</option><option value=\"only\">only</option></select></label><label class=\"screening-margin\" title=\"Flagging margin in percent (empty = configured, default 10): flag when an estimated loading reaches 100 - margin, or an estimated voltage comes within margin percent of the band width to a limit.\">Margin % <input name=\"screening_margin_pct\" type=\"number\" step=\"any\" min=\"0\" class=\"screening-margin-input\" placeholder=\"cfg\"></label>$(scen_editor_html)<a class=\"weights-link\" href=\"/powerflow/contingency-weights?case=$(_webui_urlencode(selected_casefile))\" title=\"Upload or edit the per-case N-1 weight list\">edit N-1 weights</a><button class=\"powerflow-submit contingency-submit\" type=\"submit\" name=\"contingency_mode\" value=\"true\" data-contingency-button title=\"Run an N-1 contingency analysis: take each in-service element (branch or generator, per the selector) out one at a time, check the base case, and report convergence, overloads, voltage violations, load shed, and severity. Writes contingency_n1.csv and a report to run.log.\"><span class=\"submit-spinner\" aria-hidden=\"true\"></span><span class=\"submit-label\">Contingency (N-1)</span><span class=\"submit-progress-label\" role=\"status\" aria-live=\"polite\">Running N-1...</span></button></div></form>
<script>
document.addEventListener('DOMContentLoaded', function () {
  const powerflowForm = document.getElementById('powerflow-run-form');
  if (powerflowForm !== null) {
    const copySubmitterValue = function (submitter) {
      const existing = powerflowForm.querySelector('input[type="hidden"][data-submitter-value]');
      const name = submitter ? submitter.getAttribute('name') : null;
      if (!name) {
        // Plain "Start PowerFlow run" has no name/value of its own to preserve; also
        // drop any stale submitter value left over from a bfcache-restored page so a
        // normal run can never silently inherit a previous "Diagnose" submission.
        if (existing !== null) existing.remove();
        return;
      }
      const hidden = existing !== null ? existing : document.createElement('input');
      if (existing === null) {
        hidden.type = 'hidden';
        hidden.setAttribute('data-submitter-value', 'true');
        powerflowForm.appendChild(hidden);
      }
      hidden.name = name;
      hidden.value = submitter.value;
    };
    // Fallback for browsers without SubmitEvent.submitter (older Safari): capture
    // the clicked submit button before the submit event fires.
    powerflowForm.querySelectorAll('button[type=submit][name]').forEach(function (button) {
      button.addEventListener('click', function () { copySubmitterValue(button); });
    });
    powerflowForm.addEventListener('submit', function (event) {
      // Copy the submitter's name/value into a plain hidden input FIRST: disabling
      // the submit buttons below (needed for double-submit protection) would
      // otherwise drop the submitter itself from the serialized form data per the
      // HTML form-submission spec, silently turning a "Diagnose" click into a
      // normal run.
      copySubmitterValue(event.submitter);
      powerflowForm.classList.add('is-submitting');
      powerflowForm.setAttribute('aria-busy', 'true');
      powerflowForm.querySelectorAll('button[type=submit]').forEach(function (b) { b.disabled = true; });
    });
  }
});

window.addEventListener('pageshow', function () {
  const form = document.getElementById('powerflow-run-form');
  if (form === null) return;
  form.classList.remove('is-submitting');
  form.removeAttribute('aria-busy');
  form.querySelectorAll('button[type=submit]').forEach(function (b) { b.disabled = false; });
  const staleSubmitterValue = form.querySelector('input[type="hidden"][data-submitter-value]');
  if (staleSubmitterValue !== null) staleSubmitterValue.remove();
});
</script>
$(_WEBUI_FEEDBACK_MODAL_SCRIPT)
$(_WEBUI_INFO_MENU_SCRIPT)"""
  # stage 4A block 4: the Runs page carries the run editors and the SE
  # section below the run form. The scenario editor embeds server-side
  # (its row-wiring script must arrive with the page; fetched scripts stay
  # inert), while the weights editor lazy-loads its shared fragment on
  # first open, because seeding the element names builds the net, which is
  # too expensive for every page view; the fragment is script-free, so
  # injecting it is safe.
  runs_case = basename(String(effective_case_value))
  scen_tab = ""
  if scen_is_scf && !isempty(runs_case)
    scen_fragment = _webui_scenarios_editor_tab(runs_case, effective_case_directory)
    isempty(scen_fragment) || (scen_tab = "<details class=\"editor-tab scenario-editor-tab\" id=\"scenario-editor\"><summary>Edit scenarios (case file)</summary>$(scen_fragment)</details>")
  end
  weights_tab = isempty(runs_case) ? "" : string(
    "<details class=\"editor-tab weights-editor-tab\" id=\"n1-weights\" data-weights-fragment=\"/powerflow/contingency-weights?case=$(_webui_urlencode(runs_case))&fragment=1\">",
    "<summary>Edit N-1 weights</summary><div class=\"weights-fragment-slot\"><p>Loading the weights editor…</p></div></details>",
    "<script>document.addEventListener('DOMContentLoaded',function(){var d=document.getElementById('n1-weights');if(d===null)return;var loaded=false;d.addEventListener('toggle',function(){if(!d.open||loaded)return;loaded=true;fetch(d.getAttribute('data-weights-fragment')).then(function(r){return r.text();}).then(function(h){d.querySelector('.weights-fragment-slot').innerHTML=h;}).catch(function(){d.querySelector('.weights-fragment-slot').innerHTML='<p class=\"notice\">Could not load the weights editor; use the edit N-1 weights link above.</p>';});});});</script>",
  )
  se_state = _webui_se_form_state(se_query; output_root, application_root, case_directory, config_file = selected_config_file, selected_fallback = effective_case_value)
  se_html = render_se_form(; se_state...)
  # The page carries the power flow form, scenarios, N-1 weights AND the
  # state estimation, so its heading must not name one of them: someone who
  # clicked "State estimation" (that route redirects here) read "PowerFlow
  # run" above their own section and reported the estimation as mislabelled
  # three times. The sections carry their own headings.
  return _webui_layout("Runs", string(form, scen_tab, weights_tab, se_html); header_info = info_menu)
end

const _WEBUI_RESULT_FIELDS = (
  "run_id",
  "status",
  "success",
  "converged",
  "numerical_converged",
  "solution_available",
  "iterations",
  "final_mismatch",
  "Jacobian condition",
  "reason",
  "message",
  "input_format",
  "input_format_detected",
  "native_dtf_import_used",
  "dtf_bus_count",
  "dtf_branch_count",
  "dtf_outage_count",
  "dtf_slack_bus",
  "for002_reference_used",
  "outage_validation_requested",
  "matpower_export_requested",
  "matpower_export_file",
  "dcline_status",
  "unsupported_dcline_status",
  "dtf_outage_results",
  "Q-limit enforcement mode",
  "Q-limit active-set events",
  "Classical Q-limit outer-loop passes",
  "Runtime casefile",
  "config_file",
  "started_at",
  "elapsed_seconds",
  "output_dir",
  "phase_started_at",
  "last_progress_at",
  "abort_requested_at",
  "service_status",
  "numerical_status",
  "solver_status",
  "artifact_status",
  "run_status",
  "last_phase",
  "last_heartbeat",
)

const _WEBUI_IMPORTANT_RESULT_FIELDS = Set(("converged", "numerical_converged", "solution_available", "iterations", "final_mismatch", "reason"))

function _webui_result_value(result::AbstractDict, field::AbstractString)
  if field == "Q-limit enforcement mode"
    return get(result, "qlimit_enforcement_mode", get(get(result, "metadata", Dict{String,Any}()), "qlimit_enforcement_mode", "n/a"))
  elseif field == "Q-limit active-set events"
    metadata = get(result, "metadata", Dict{String,Any}())
    return get(metadata, "q_limit_active_set_events", get(metadata, "pv_pq_switching_events", "n/a"))
  elseif field == "Classical Q-limit outer-loop passes"
    metadata = get(result, "metadata", Dict{String,Any}())
    return get(metadata, "q_limit_classic_outer_loop_passes", "n/a")
  elseif field == "Runtime casefile"
    return get(result, "runtime_casefile", get(get(result, "metadata", Dict{String,Any}()), "runtime_casefile", "n/a"))
  elseif field == "Jacobian condition"
    # single source of formatting: the metadata carries the identical line
    # the classic result log prints; estimate plus verdict is the fallback
    # for older run records, n/a for runs without a rectangular NR thunk
    metadata = get(result, "metadata", Dict{String,Any}())
    line = get(metadata, "jacobian_condition_line", nothing)
    line isa AbstractString && !isempty(line) && return line
    kappa = get(metadata, "jacobian_condition_estimate", nothing)
    verdict = get(metadata, "jacobian_condition_verdict", nothing)
    kappa isa Real && return string("kappa1(J) = ", round(Float64(kappa), sigdigits = 3), verdict isa AbstractString && !isempty(verdict) ? ", " * verdict : "")
    return "n/a"
  end
  metadata = get(result, "metadata", Dict{String,Any}())
  value = get(result, field, get(metadata, field, nothing))
  value === nothing && return field in _WEBUI_IMPORTANT_RESULT_FIELDS ? "n/a" : ""
  value isa AbstractString && isempty(value) && return field in _WEBUI_IMPORTANT_RESULT_FIELDS ? "n/a" : ""
  return value
end

function _webui_solver_status(result::AbstractDict)::String
  converged = get(result, "converged", get(result, "numerical_converged", nothing))
  converged === true && return "converged"
  converged === false && return "not_converged"
  final_outcome = get(result, "final_outcome", nothing)
  final_outcome isa AbstractDict && haskey(final_outcome, "converged") && return final_outcome["converged"] === true ? "converged" : "not_converged"
  return "n/a"
end

function _webui_result_successful(result::AbstractDict)::Bool
  get(result, "success", false) === true && return true
  lowercase(string(get(result, "status", ""))) in ("succeeded", "success", "converged", "ok") && return true
  return _webui_solver_status(result) == "converged"
end

function _webui_case_settings_save_section(result::AbstractDict)::String
  run_id = String(get(result, "run_id", ""))
  isempty(run_id) && return ""
  metadata = get(result, "metadata", Dict{String,Any}())
  casefile = get(result, "runtime_casefile", get(metadata, "runtime_casefile", get(result, "casefile", "")))
  isempty(String(casefile)) && return ""
  successful = _webui_result_successful(result)
  message = successful ? "This run converged. You can save the current Web UI settings as the default profile for this case." : "This run did not converge. Saving these settings is possible, but they may not be a good default for this case."
  button = successful ? "Save settings for this case" : "Save these settings anyway"
  override = successful ? "" : "<input type=\"hidden\" name=\"override_non_success\" value=\"true\">"
  warning_class = successful ? "info" : "warning"
  return """
<section class=\"panel case-settings-save\"><h2>Case settings</h2>
<p class=\"alert $(warning_class)\">$(_webui_escape(message))</p>
<p><strong>Case:</strong> <code>$(_webui_escape(casefile))</code></p>
<form method=\"post\" action=\"/powerflow/result/$(_webui_urlencode(run_id))/case-settings/save\">
$(override)<button type=\"submit\">$(_webui_escape(button))</button>
</form></section>"""
end

function render_config_refresh_result(result::AbstractDict)::String
  list_html(items) = isempty(items) ? "<li>none</li>" : join(("<li><code>$(_webui_escape(String(item)))</code></li>" for item in items), "")
  status = get(result, "written", false) ? "Configuration refreshed and written." : get(result, "changed", false) ? "Configuration refresh changes are available." : "Configuration is already current."
  config_path = String(get(result, "config_file", ""))
  config_line = isempty(config_path) ? "" : "<p><strong>Configuration file:</strong> <code>$(_webui_escape(config_path))</code></p>"
  backup = isempty(String(get(result, "backup_path", ""))) ? "" : "<p><strong>Backup:</strong> <code>$(_webui_escape(String(result["backup_path"])))</code></p>"
  restart = get(result, "written", false) ? "<p class=\"alert warning\"><strong>Restart or reload the Web UI</strong> before relying on the updated configuration.</p>" : ""
  download = get(result, "downloadable", false) ? "<form method=\"post\" action=\"/powerflow/config/download\"><textarea name=\"refreshed_text\" hidden>$(_webui_escape(String(result["refreshed_text"])))</textarea><button type=\"submit\">Download refreshed YAML</button></form>" : ""
  content = """
<section class=\"panel\"><h1>Configuration refresh</h1><p>$(status)</p>$(config_line)$(backup)$(restart)
<h2>Missing keys added from the template</h2><ul>$(list_html(get(result, "missing_keys", String[])))</ul>
<h2>Deprecated aliases normalized</h2><ul>$(list_html(get(result, "normalized_keys", String[])))</ul>
<h2>Duplicate keys detected</h2><ul>$(list_html(get(result, "duplicate_keys", String[])))</ul>
<h2>Warnings</h2><ul>$(list_html(get(result, "warnings", String[])))</ul>
$(download)<details><summary>Refreshed YAML preview</summary><pre class=\"artifact-text\">$(_webui_escape(String(result["refreshed_text"])))</pre></details>
<p><a class=\"button\" href=\"/powerflow\">Back to PowerFlow</a></p></section>
"""
  return _webui_layout("Configuration refresh", content; show_back = true)
end
"""
    _webui_wrong_branch_badge(result) -> Union{Nothing,String}

Builds a `status-badge` HTML fragment for the wrong-branch detection outcome
carried in `result["metadata"]["wrong_branch_status"]`. Returns `nothing` for
`not_checked` (detection off or never reached) so the summary row is omitted
entirely rather than showing an empty/uninformative badge.
"""
function _webui_wrong_branch_badge(result::AbstractDict)::Union{Nothing,String}
  metadata = get(result, "metadata", Dict{String,Any}())
  status = lowercase(String(get(metadata, "wrong_branch_status", "not_checked")))
  status == "not_checked" && return nothing
  reason = String(get(metadata, "wrong_branch_reason", "unknown"))
  css_class = status == "ok" ? "status-success" : status == "fail" ? "status-error" : "status-warning"
  label = status == "ok" ? "ok" : reason
  return "<span class=\"status-badge $(css_class)\">$(_webui_escape(label))</span>"
end

# Short-circuit runs: compact summary from the run metadata. A flagged
# Ik''max is a lower bound (skipped/defaulted contributions), so flags render
# as a warning badge rather than plain text. Returns nothing for other runs.
function _webui_short_circuit_summary(result::AbstractDict)::Union{Nothing,String}
  metadata = get(result, "metadata", Dict{String,Any}())
  metadata isa AbstractDict || return nothing
  get(metadata, "run_mode", "") == "short_circuit" || return nothing
  worst_bus = string(get(metadata, "sc_worst_bus", ""))
  ik = get(metadata, "sc_max_ik_kA", nothing)
  fmt = x -> x isa Real && isfinite(x) ? string(round(Float64(x); sigdigits = 5)) : "n/a"
  flagged = get(metadata, "sc_flagged_rows", 0)
  rows = get(metadata, "sc_case_rows", 0)
  text = string("worst Ik''max ", fmt(ik), " kA @ ", worst_bus, " (", rows, " buses)")
  badge = flagged isa Real && flagged > 0 ? " <span class=\"status-badge status-warning\">$(flagged) flagged — lower bound</span>" : ""
  return "<code>" * _webui_escape(text) * "</code>" * badge
end

# N-1 contingency runs: one summary row with the outcome counts and the worst
# loading. Names the slack-unit outage (a generator N-1 that removes the only
# reference) so it does not read as a tool failure. Returns nothing otherwise.
function _webui_contingency_summary(result::AbstractDict)::Union{Nothing,String}
  metadata = get(result, "metadata", Dict{String,Any}())
  metadata isa AbstractDict || return nothing
  get(metadata, "run_mode", "") == "contingency" || return nothing
  kind = string(get(metadata, "contingency_kind", "branch")) == "gen" ? "generator" : "branch"
  n = get(metadata, "contingency_cases", 0)
  conv = get(metadata, "contingency_converged", 0)
  isl = get(metadata, "contingency_islanded", 0)
  nonconv = get(metadata, "contingency_nonconverged", 0)
  no_slack = get(metadata, "contingency_no_slack", 0)
  shed = get(metadata, "contingency_total_shed_mw", 0.0)
  worst = get(metadata, "contingency_worst_loading_pct", NaN)
  weighted = get(metadata, "contingency_weights_applied", false) === true
  wcount = get(metadata, "contingency_weighted_cases", 0)
  fmt = x -> x isa Real && isfinite(x) ? string(round(Float64(x); digits = 1)) : "n/a"
  # state the weighting explicitly, so a severity ranking that used a weight file
  # is never read as if it were unweighted
  weight_note = weighted ? ", weighted ($(wcount) case$(wcount == 1 ? "" : "s"))" : ", unweighted"
  # screening (step 6): screened rows carry first-order estimates, not full
  # solves; the share must be visible wherever the counts are
  n_screened = get(metadata, "contingency_screened", 0)
  screen_note = n_screened isa Real && n_screened > 0 ? string(", ", n_screened, " screened (estimates)") : ""
  source = string(get(metadata, "contingency_cases_source", ""))
  source_note = isempty(source) || source == "generated" ? "" : string(", source ", source)
  text = string(kind, " N-1: ", conv, "/", n, " converged, ", isl, " islanded (", fmt(shed), " MW shed), worst loading ", fmt(worst), "%", weight_note, screen_note, source_note)
  badge = ""
  if nonconv isa Real && nonconv > 0
    label = no_slack isa Real && no_slack > 0 ? "$(nonconv) non-converged incl. $(no_slack) that removed the only slack (auto_slack resolves it)" : "$(nonconv) non-converged"
    badge = " <span class=\"status-badge status-warning\">$(label)</span>"
  end
  return "<code>" * _webui_escape(text) * "</code>" * badge
end

# N-1 contingency weights editor (issue #331 Phase 5 follow-up). `elements` is
# the (already capped) list of element names seeded from the case; `stored` maps
# element name to its stored weight; `raw_text` is the current file content for
# the free-text editor. `net_error` set means the seeded table could not be
# built (e.g. a non-N-1 format) so only the raw editor / upload are shown.
function render_contingency_weights_editor(; case::AbstractString, cases::AbstractVector{<:AbstractString}, elements::AbstractVector{<:AbstractString}, stored::AbstractDict, raw_text::AbstractString, message::AbstractString = "", filter::AbstractString = "", total_count::Int = 0, net_error::AbstractString = "")::String
  esc = _webui_escape
  msg_html = isempty(message) ? "" : "<p class=\"notice\">$(esc(message))</p>"
  options = join(("<option value=\"$(esc(c))\"$(c == case ? " selected" : "")>$(esc(c))</option>" for c in cases), "")
  selector = "<form method=\"get\" action=\"/powerflow/contingency-weights\" class=\"panel\"><label>Case <select name=\"case\" onchange=\"this.form.submit()\">$(options)</select></label> <noscript><button type=\"submit\">Load</button></noscript></form>"
  if isempty(case)
    return _webui_layout("N-1 weights", string(selector, msg_html, "<p>Select a case to edit its N-1 contingency weights. Weights only reorder the severity ranking; they never skip a case.</p>"); show_back = true)
  end
  return _webui_layout("N-1 weights: $(case)", string(selector, _webui_weights_editor_fragment(; case, elements, stored, raw_text, message, filter, total_count, net_error)); show_back = true)
end

"""
    _webui_weights_editor_fragment(; case, elements, stored, raw_text, message, filter, total_count, net_error) -> String

The N-1 weights editor itself (upload/download, seeded table, raw CSV),
without page selector and layout: shared by the standalone
/powerflow/contingency-weights page and the Runs page's weights tab
(stage 4A block 4). The tab lazy-loads this fragment through the
standalone route with fragment=1, because seeding the element names
builds the net, which is too expensive for every Runs-page render; the
fragment stays script-free so injecting it is safe.
"""
function _webui_weights_editor_fragment(; case::AbstractString, elements::AbstractVector{<:AbstractString} = String[], stored::AbstractDict = Dict{String,Float64}(), raw_text::AbstractString = "", message::AbstractString = "", filter::AbstractString = "", total_count::Int = 0, net_error::AbstractString = "")::String
  esc = _webui_escape
  msg_html = isempty(message) ? "" : "<p class=\"notice\">$(esc(message))</p>"
  cesc = esc(case)
  cenc = _webui_urlencode(case)
  upload = "<section class=\"panel\"><h2>Upload / download</h2><form method=\"post\" action=\"/powerflow/contingency-weights/upload\" enctype=\"multipart/form-data\"><input type=\"hidden\" name=\"casefile\" value=\"$(cesc)\">$(_webui_file_input("weights_file"; accept = ".csv", required = true)) <button type=\"submit\">Upload (replaces existing)</button></form><p><a href=\"/powerflow/contingency-weights/download?case=$(cenc)\">download current file</a></p><form method=\"post\" action=\"/powerflow/contingency-weights/reset\"><input type=\"hidden\" name=\"casefile\" value=\"$(cesc)\"><button type=\"submit\">reset (delete the weight file)</button></form></section>"
  table_html = if !isempty(net_error)
    "<section class=\"panel\"><h2>Weights table</h2><p class=\"notice\">Element names could not be listed for this case ($(esc(net_error))). Edit the raw CSV below or upload a file.</p></section>"
  else
    rows = join(("<tr><td>$(esc(e))</td><td><input type=\"hidden\" name=\"element\" value=\"$(esc(e))\"><input name=\"weight\" type=\"number\" step=\"any\" min=\"0\" value=\"$(get(stored, e, 1.0))\"></td></tr>" for e in elements), "")
    cap_note = length(elements) < total_count ? "<p class=\"notice\">Showing $(length(elements)) of $(total_count) elements; use the filter or edit the raw CSV.</p>" : ""
    filter_form = "<form method=\"get\" action=\"/powerflow/contingency-weights\"><input type=\"hidden\" name=\"case\" value=\"$(cesc)\"><label>Filter by name <input name=\"filter\" value=\"$(esc(filter))\"></label> <button type=\"submit\">Filter</button></form>"
    "<section class=\"panel\"><h2>Weights table ($(total_count) elements)</h2>$(filter_form)$(cap_note)<form method=\"post\" action=\"/powerflow/contingency-weights/save\"><input type=\"hidden\" name=\"casefile\" value=\"$(cesc)\"><table><thead><tr><th>element</th><th>weight</th></tr></thead><tbody>$(rows)</tbody></table><button type=\"submit\">Save table (rows at 1.0 are omitted)</button></form></section>"
  end
  textarea = "<section class=\"panel\"><h2>Raw CSV</h2><form method=\"post\" action=\"/powerflow/contingency-weights/save\"><input type=\"hidden\" name=\"casefile\" value=\"$(cesc)\"><textarea name=\"weights_text\" rows=\"12\" cols=\"64\">$(esc(raw_text))</textarea><br><button type=\"submit\">Save raw CSV</button></form></section>"
  return string(msg_html, upload, table_html, textarea)
end

# Scenario editor page (scenario task step 6, maintainer revision
# 2026-09-03): SCF cases only; a list of the case file's scenarios (name,
# weight, op count; edit / duplicate / delete / new) and a scenario form
# with op rows (op, target, component filtered by class, then the fields
# the op needs). Validation is server side; errors arrive as `error_text`
# with scenario name and op index and render next to the form. Saving
# writes the scenarios block into the SCF through the existing writer;
# nothing is kept in the browser.
function render_scenarios_editor(;
  case::AbstractString,
  cases::AbstractVector{<:AbstractString},
  scenarios::AbstractVector = [],
  component_options::AbstractVector = [],
  form_scenario = nothing,
  original_name::AbstractString = "",
  message::AbstractString = "",
  error_text::AbstractString = "",
  file_digest::AbstractString = "",
)::String
  esc = _webui_escape
  msg_html = isempty(message) ? "" : "<p class=\"notice\">$(esc(message))</p>"
  err_html = isempty(error_text) ? "" : "<p class=\"notice error scenario-error\">$(esc(error_text))</p>"
  options = join(("<option value=\"$(esc(c))\"$(c == case ? " selected" : "")>$(esc(c))</option>" for c in cases), "")
  selector = "<form method=\"get\" action=\"/powerflow/scenarios\" class=\"panel\"><label>SCF case <select name=\"case\" onchange=\"this.form.submit()\">$(options)</select></label> <noscript><button type=\"submit\">Load</button></noscript></form>"
  if isempty(case)
    return _webui_layout("Scenarios", string(selector, msg_html, "<p>Select an SCF case (.scf.json) to edit its scenarios block. Scenarios need an SCF case; export other formats as SCF first.</p>"); show_back = true)
  end
  return _webui_layout("Scenarios: $(case)", string(selector, _webui_scenarios_editor_fragment(; case, scenarios, component_options, form_scenario, original_name, message, error_text, file_digest)); show_back = true)
end

"""
    _webui_scenarios_editor_fragment(; case, scenarios, component_options, form_scenario, original_name, message, error_text, file_digest) -> String

The scenario editor itself (list plus op form and its row-wiring script),
without page selector and layout: shared by the standalone
/powerflow/scenarios page and the Runs page's scenario tab (stage 4A
block 4). Embedded server-side because the script must arrive with the
initial page (fetched scripts stay inert).
"""
function _webui_scenarios_editor_fragment(;
  case::AbstractString,
  scenarios::AbstractVector = [],
  component_options::AbstractVector = [],
  form_scenario = nothing,
  original_name::AbstractString = "",
  message::AbstractString = "",
  error_text::AbstractString = "",
  file_digest::AbstractString = "",
)::String
  esc = _webui_escape
  msg_html = isempty(message) ? "" : "<p class=\"notice\">$(esc(message))</p>"
  err_html = isempty(error_text) ? "" : "<p class=\"notice error scenario-error\">$(esc(error_text))</p>"
  cesc = esc(case)
  cenc = _webui_urlencode(case)
  list_rows = join((begin
    n = esc(String(s.name))
    ne = _webui_urlencode(String(s.name))
    string(
      "<tr><td>$(n)</td><td>$(s.weight)</td><td>$(length(s.ops))</td><td>",
      "<a href=\"/powerflow/scenarios?case=$(cenc)&scenario=$(ne)\">edit</a> ",
      "<a href=\"/powerflow/scenarios?case=$(cenc)&scenario=$(ne)&mode=duplicate\">duplicate</a> ",
      "<form class=\"inline-form\" method=\"post\" action=\"/powerflow/scenarios/delete\"><input type=\"hidden\" name=\"casefile\" value=\"$(cesc)\"><input type=\"hidden\" name=\"scenario\" value=\"$(n)\"><input type=\"hidden\" name=\"file_digest\" value=\"$(esc(file_digest))\"><button type=\"submit\">delete</button></form>",
      "</td></tr>",
    )
  end for s in scenarios), "")
  list_html = string(
    "<section class=\"panel\"><h2>Scenarios in $(cesc) ($(length(scenarios)))</h2>",
    isempty(scenarios) ? "<p>The case carries no scenarios block yet.</p>" : "<table><thead><tr><th>name</th><th>weight</th><th>ops</th><th>actions</th></tr></thead><tbody>$(list_rows)</tbody></table>",
    "<p><a href=\"/powerflow/scenarios?case=$(cenc)\">new scenario</a> | run from the main page with scenario source \"case file scenarios\"</p></section>",
  )
  comp_opts_html = join(("<option value=\"$(id)\" data-class=\"$(cls)\">$(esc(labeltext))</option>" for (cls, id, labeltext) in component_options), "")
  fs_name = form_scenario === nothing ? "" : String(form_scenario.name)
  fs_weight = form_scenario === nothing ? 1.0 : form_scenario.weight
  fs_ops = form_scenario === nothing ? [nothing] : (isempty(form_scenario.ops) ? [nothing] : form_scenario.ops)
  op_row(op) = begin
    sel(v, w) = v == w ? " selected" : ""
    o = op === nothing ? :status : op.op
    t = op === nothing ? :branch : op.target
    id = op === nothing ? nothing : op.id
    f = op === nothing ? nothing : op.field
    v = op === nothing ? nothing : op.value
    fa = op === nothing ? nothing : op.factor
    string(
      "<tr class=\"scenario-op-row\">",
      "<td><select name=\"op_op\"><option value=\"status\"$(sel(o, :status))>status</option><option value=\"set\"$(sel(o, :set))>set</option><option value=\"scale\"$(sel(o, :scale))>scale</option></select></td>",
      "<td><select name=\"op_target\" class=\"op-target\">",
      join(("<option value=\"$(cls)\"$(sel(t, cls))>$(cls)</option>" for cls in (:branch, :transformer, :generator, :load, :shunt, :link, :external_grid)), ""),
      "</select></td>",
      "<td><select name=\"op_component\" class=\"op-component\"><option value=\"\"></option>",
      join(("<option value=\"$(cid)\" data-class=\"$(cls)\"$(id == cid ? " selected" : "")>$(esc(labeltext))</option>" for (cls, cid, labeltext) in component_options), ""),
      "</select></td>",
      "<td><select name=\"op_field\"><option value=\"\"></option>",
      join(("<option value=\"$(fld)\"$(sel(f, fld))>$(fld)</option>" for fld in (:p, :q, :vm_pu, :tap_pos, :b_pu, :angle_deg)), ""),
      "</select></td>",
      "<td><input name=\"op_value\" type=\"number\" step=\"any\" value=\"$(v === nothing ? "" : v)\" placeholder=\"value\"></td>",
      "<td><input name=\"op_factor\" type=\"number\" step=\"any\" value=\"$(fa === nothing ? "" : fa)\" placeholder=\"factor\"></td>",
      "<td><button type=\"button\" class=\"op-remove\" onclick=\"scenarioRemoveRow(this)\">remove</button></td>",
      "</tr>",
    )
  end
  op_rows = join((op_row(op) for op in fs_ops), "")
  form_html = string(
    "<section class=\"panel\"><h2>",
    isempty(original_name) ? "New scenario" : "Edit scenario: $(esc(original_name))",
    "</h2>$(err_html)<form method=\"post\" action=\"/powerflow/scenarios/save\">",
    "<input type=\"hidden\" name=\"casefile\" value=\"$(cesc)\">",
    "<input type=\"hidden\" name=\"original_name\" value=\"$(esc(original_name))\">",
    "<input type=\"hidden\" name=\"file_digest\" value=\"$(esc(file_digest))\">",
    "<label>Name <input name=\"scenario_name\" value=\"$(esc(fs_name))\" required></label> ",
    "<label>Weight <input name=\"scenario_weight\" type=\"number\" step=\"any\" min=\"0\" value=\"$(fs_weight)\"></label>",
    "<table class=\"scenario-ops\"><thead><tr><th>op</th><th>target</th><th>component</th><th>field (set)</th><th>value</th><th>factor (scale)</th><th></th></tr></thead><tbody id=\"scenario-op-rows\">$(op_rows)</tbody></table>",
    "<button type=\"button\" onclick=\"scenarioAddRow()\">add row</button> ",
    "<button type=\"submit\">Save into case file</button>",
    "</form></section>",
    # the component select filters by the row's target class; without
    # JavaScript the full list stays selectable and the server validates
    "<script>function scenarioFilterRow(row){var t=row.querySelector('.op-target').value;row.querySelectorAll('.op-component option').forEach(function(o){o.hidden=o.value!==''&&o.getAttribute('data-class')!==t;});}function scenarioWireRow(row){row.querySelector('.op-target').addEventListener('change',function(){scenarioFilterRow(row);});scenarioFilterRow(row);}function scenarioAddRow(){var body=document.getElementById('scenario-op-rows');var row=body.rows[0].cloneNode(true);row.querySelectorAll('input').forEach(function(i){i.value='';});scenarioWireRow(row);body.appendChild(row);}function scenarioRemoveRow(btn){var body=document.getElementById('scenario-op-rows');if(body.rows.length>1){btn.closest('tr').remove();}}document.querySelectorAll('#scenario-op-rows tr').forEach(scenarioWireRow);</script>",
  )
  return string(msg_html, list_html, form_html)
end

# Import-analysis runs: one summary row with the verdict and the gap counts.
# Returns nothing for other runs.
function _webui_import_analysis_summary(result::AbstractDict)::Union{Nothing,String}
  metadata = get(result, "metadata", Dict{String,Any}())
  metadata isa AbstractDict || return nothing
  get(metadata, "run_mode", "") == "import_analysis" || return nothing
  missing_deps = get(metadata, "import_analysis_missing_dependencies", 0)
  unresolved = get(metadata, "import_analysis_unresolved_refs", 0)
  verdict = string(get(metadata, "import_analysis_verdict", ""))
  not_importable = string(get(result, "reason", "")) == "import_analysis_not_importable"
  badge = if not_importable
    "<span class=\"status-badge status-error\">not importable — $(missing_deps) missing dependency(ies) · $(unresolved) unresolved reference(s)</span>"
  elseif unresolved isa Real && unresolved > 0
    "<span class=\"status-badge status-warning\">importable · $(unresolved) non-fatal unresolved reference(s)</span>"
  else
    "<span class=\"status-badge status-success\">importable</span>"
  end
  text = isempty(verdict) ? "" : " <code>" * _webui_escape(verdict) * "</code>"
  return badge * text
end

# Optional CGMES export artifact: one summary row from the run metadata.
# Export notices (per-unit content the profiles cannot carry, e.g. a dropped
# transformer phase shift) render as a warning badge so a partial export is
# never mistaken for a full one.
# Returns nothing when the run did not request the export.
function _webui_cgmes_export_summary(result::AbstractDict)::Union{Nothing,String}
  metadata = get(result, "metadata", Dict{String,Any}())
  metadata isa AbstractDict || return nothing
  status = get(metadata, "cgmes_export_status", nothing)
  status === nothing && return nothing
  status = string(status)
  if status == "completed"
    text = string(get(metadata, "cgmes_export_files", ""))
    sc_lines = get(metadata, "cgmes_export_sc_lines", 0)
    sc_lines isa Real && sc_lines > 0 && (text = string(text, " · zero-sequence data on ", Int(sc_lines), " line(s)"))
    notices = string(get(metadata, "cgmes_export_notices", ""))
    badge = isempty(notices) ? "" : " <span class=\"status-badge status-warning\">$(_webui_escape(notices))</span>"
    return "<code>" * _webui_escape(text) * "</code>" * badge
  elseif status == "skipped"
    return "<span class=\"status-badge status-warning\">skipped — $(_webui_escape(string(get(metadata, "cgmes_export_skip_reason", ""))))</span>"
  else
    return "<span class=\"status-badge status-error\">failed — $(_webui_escape(string(get(metadata, "cgmes_export_error", ""))))</span>"
  end
end

# CGMES runs: compact SV-comparison summary from the run metadata (written by
# the mandatory compareWithSV check). Returns nothing for non-CGMES runs —
# the metadata keys only exist on the CGMES path.
function _webui_sv_compare_summary(result::AbstractDict)::Union{Nothing,String}
  metadata = get(result, "metadata", Dict{String,Any}())
  metadata isa AbstractDict || return nothing
  haskey(metadata, "cgmes_sv_compare_status") || return nothing
  status = string(get(metadata, "cgmes_sv_compare_status", "unavailable"))
  fmt = x -> x isa Real && isfinite(x) ? string(round(Float64(x); sigdigits = 3)) : "n/a"
  status == "unavailable" && return "<code>unavailable</code>"
  start_values = string(get(metadata, "cgmes_start_values", "?"))
  ref_offset = get(metadata, "cgmes_sv_compare_va_ref_offset_deg", nothing)
  offset_note = ref_offset isa Real && isfinite(ref_offset) && abs(ref_offset) > 0.01 ? string(" · ref offset ", fmt(ref_offset), "° removed") : ""
  text = string(
    "max|dvm| ", fmt(get(metadata, "cgmes_sv_compare_max_dvm", nothing)), " pu · rms ", fmt(get(metadata, "cgmes_sv_compare_rms_dvm", nothing)),
    " · max|dva| ", fmt(get(metadata, "cgmes_sv_compare_max_dva", nothing)), "° · rms ", fmt(get(metadata, "cgmes_sv_compare_rms_dva", nothing)), "°",
    offset_note,
    " (n=", get(metadata, "cgmes_sv_compare_n", 0), ", ", start_values, " start",
    status == "converged" ? "" : ", " * status,
    ")",
  )
  return "<code>" * _webui_escape(text) * "</code>"
end

## Page title of a finished run: the run KIND, not always "PowerFlow"
## (maintainer 2026-09-04: a state estimation whose result page says
## "PowerFlow result" reads as if the wrong calculation had run, and that
## is exactly how it was reported). The run modes come from the service
## metadata written by the run itself.
function _webui_result_page_title(result::AbstractDict)::String
  metadata = get(result, "metadata", Dict{String,Any}())
  mode = metadata isa AbstractDict ? lowercase(string(get(metadata, "run_mode", ""))) : ""
  mode == "se" && return "State estimation result"
  mode == "powerflow_se_start" && return "PowerFlow result (started from a state estimate)"
  mode == "short_circuit" && return "Short-circuit result"
  mode == "contingency" && return "N-1 contingency result"
  mode == "import_analysis" && return "Import analysis result"
  mode == "diagnose" && return "Diagnostic result"
  return "PowerFlow result"
end

## Method name shown in the run history and on the detail page. Derived
## from the RUN KIND for every kind that does not run a power flow, because
## the history used to fabricate "rectangular" for all of them: a state
## estimation solves weighted least squares, a short circuit is a direct
## IEC 60909 solve, an import analysis solves nothing. Deriving at display
## time also repairs index entries written before the fix (maintainer,
## 2026-09-06).
function _webui_run_method(entry::AbstractDict, stored::AbstractString)
  mode = String(get(entry, "run_mode", ""))
  mode == "se" && return "wls"
  mode == "short_circuit" && return "iec60909"
  mode == "import_analysis" && return ""
  isempty(stored) ? "rectangular" : String(stored)
end

function render_powerflow_result(result::AbstractDict)::String
  run_id = get(result, "run_id", "")
  rows = join(("<tr><th>$(_webui_escape(field))</th><td>$(_webui_escape(_webui_result_value(result, field)))</td></tr>" for field in _WEBUI_RESULT_FIELDS), "")
  status = lowercase(string(get(result, "status", "unknown")))
  active = status in _WEBUI_ACTIVE_RUN_STATUSES
  diagnose_probe = _webui_is_completed_diagnose(result)
  status_text = diagnose_probe ? "diagnosed" : status
  status_badge = _webui_status_badge(webui_status_class(result), status_text, status)
  final_outcome_value = get(result, "final_outcome", nothing)
  stored_solver = final_outcome_value isa AbstractDict ? String(get(final_outcome_value, "solver", "")) : ""
  meta_for_mode = get(result, "metadata", Dict{String,Any}())
  solver_name = _webui_run_method(meta_for_mode isa AbstractDict ? meta_for_mode : Dict{String,Any}(), stored_solver)
  # The case by name and the phase, up top where a reader looks first. The
  # `casefile`/`resolved_casefile` rows below carried the same path twice and
  # broke the layout; `final_outcome` said nothing a reader could use
  # (maintainer 2026-09-09). The path survives in the tooltip.
  # A live job snapshot carries `nothing` for what the run has not produced
  # yet, and `string(nothing)` is the word "nothing": the card said so during
  # every run and named the case only at the end (maintainer 2026-09-10).
  # First non-empty of the resolved path and the requested one; the requested
  # one is known from the first second.
  text = value -> (value === nothing || value === missing) ? "" : String(strip(string(value)))
  case_path = text(get(result, "resolved_casefile", nothing))
  isempty(case_path) && (case_path = text(get(result, "casefile", nothing)))
  case_card = ("Case", "<code title=\"$(_webui_escape(case_path))\">$(_webui_escape(isempty(case_path) ? "n/a" : basename(case_path)))</code>")
  phase = text(get(result, "current_phase", nothing))
  isempty(phase) && (phase = text(get(result, "last_phase", nothing)))
  isempty(phase) && (phase = "n/a")
  phase_card = ("Phase", "<code>$(_webui_escape(phase))</code>")
  summary_rows = if active
    (("Run status", status_badge), case_card, phase_card, ("Elapsed time", "<strong>$(_webui_escape(_format_elapsed_duration(_webui_elapsed_seconds(result, active))))</strong>"))
  else
    base = [("Run status", status_badge), case_card, phase_card, ("Solver", "<code>$(_webui_escape(solver_name))</code>")]
    solver_name == "dc" && push!(base, ("Model", "<span class=\"status-badge status-info\">DC solution</span>"))
    solver_elapsed = _webui_solver_elapsed_seconds(result)
    solver_elapsed === nothing || push!(base, ("Solver time", "<strong>$(_webui_escape(_format_elapsed_duration(solver_elapsed)))</strong>"))
    push!(base, ("Total time", "<strong>$(_webui_escape(_format_elapsed_duration(_webui_total_elapsed_seconds(result))))</strong>"))
    wrong_branch_badge = _webui_wrong_branch_badge(result)
    wrong_branch_badge === nothing || push!(base, ("Wrong-branch check", wrong_branch_badge))
    control_summary = _webui_control_summary(result)
    control_summary === nothing || push!(base, ("Controllers", "<code>$(_webui_escape(control_summary))</code>"))
    sv_summary = _webui_sv_compare_summary(result)
    sv_summary === nothing || push!(base, ("SV comparison", sv_summary))
    sc_summary = _webui_short_circuit_summary(result)
    sc_summary === nothing || push!(base, ("Short circuit", sc_summary))
    ia_summary = _webui_import_analysis_summary(result)
    ia_summary === nothing || push!(base, ("Import analysis", ia_summary))
    contingency_summary = _webui_contingency_summary(result)
    contingency_summary === nothing || push!(base, ("Contingency (N-1)", contingency_summary))
    cgmes_export_summary = _webui_cgmes_export_summary(result)
    cgmes_export_summary === nothing || push!(base, ("CGMES export", cgmes_export_summary))
    auto_summary = _webui_auto_mode_summary(result)
    auto_summary === nothing || push!(base, ("Auto mode", auto_summary))
    se_summary = _webui_se_summary(result)
    se_summary === nothing || push!(base, ("State estimation", se_summary))
    se_start_summary = _webui_se_start_summary(result)
    se_start_summary === nothing || push!(base, ("SE-started PF", se_start_summary))
    Tuple(base)
  end
  result_summary = "<div class=\"result-summary\">" * join(("<div$(label in ("Elapsed time", "Solver time", "Total time") ? " class=\"runtime-card\"" : "")><span class=\"summary-label\">$(label)</span>$(value)</div>" for (label, value) in summary_rows), "") * "</div>"
  abort_form = status in ("queued", "running") ? "<form method=\"post\" action=\"/powerflow/abort/$(_webui_urlencode(run_id))\"><button type=\"submit\" class=\"danger-button\">Abort run</button></form>" : ""
  active_hint = active ? "<p class=\"status-refresh-hint\">This page refreshes automatically while the run is active.</p>" : ""
  abort_hint = status == "aborting" ? "<p>Aborting requested. Current phase: <code>$(_webui_escape(get(result, "current_phase", "unknown")))</code>.</p><p>This phase may need to finish before cancellation is observed.</p>" : ""
  hard_reset =
    status == "aborting" && get(result, "hard_reset_available", false) ?
    "<div class=\"alert warning\"><p>Abort is still pending. The calculation may be in a non-interruptible numerical call.</p><form method=\"post\" action=\"/powerflow/hard-reset/$(_webui_urlencode(run_id))\"><button type=\"submit\" class=\"danger-button\">Hard reset Web UI</button></form></div>" :
    ""
  interrupted =
    status in ("aborted_unknown", "interrupted_unknown") ?
    "<div class=\"alert warning\"><p>Run state was recovered after Web UI restart.</p><p>Last known phase: <code>$(_webui_escape(get(result, "last_phase", get(result, "current_phase", "unknown"))))</code>.</p><p>No live solver task is attached anymore. Partial artifacts may be available.</p></div>" :
    ""
  # bad data must be findable in one click: a direct download when the run
  # produced the se_bad_data.csv extract (suspicious/eliminated rows with
  # their network locations)
  bad_data_link = begin
    od = String(get(result, "output_dir", ""))
    !isempty(od) && isfile(joinpath(od, "se_bad_data.csv")) ? "<a class=\"button\" href=\"/powerflow/artifact/$(_webui_urlencode(run_id))/se_bad_data.csv\">Bad data (se_bad_data.csv)</a>" : ""
  end
  links =
    isempty(String(run_id)) ? "" :
    "$(abort_hint)$(hard_reset)$(interrupted)<div class=\"actions\">$(abort_form)<a class=\"button\" href=\"/powerflow/artifacts/$(_webui_urlencode(run_id))\">View artifacts</a>$(bad_data_link)<a class=\"button\" href=\"/powerflow/artifact-zip/$(_webui_urlencode(run_id))\">Download all artifacts as ZIP</a><a class=\"button\" href=\"/powerflow/result/$(_webui_urlencode(run_id))\">Refresh status</a></div>"
  refresh_url = active && !isempty(String(run_id)) ? "/powerflow/result/$(_webui_urlencode(run_id))?autorefresh=1" : nothing
  save_section = active ? "" : _webui_case_settings_save_section(result)
  se_tap_section = active ? "" : _webui_se_tap_section(result)
  se_topology_section = active ? "" : _webui_se_topology_section(result)
  se_chain_section = active ? "" : _webui_se_chain_section(result)
  n1_table_section = active ? "" : _webui_contingency_table_section(result)
  metadata = get(result, "metadata", Dict{String,Any}())
  override_source = String(get(result, "config_override_source", get(metadata, "config_override_source", "")))
  runtime_notice = override_source == "user_yaml" ? "<div class=\"alert warning\"><strong>Web UI settings ignored.</strong> This run used YAML/default configuration values because the run was submitted with Web UI settings ignored.</div>" : ""
  # A diagnose run is a one-step residual probe. Without this line the page
  # reads like a crashed power flow, and the residual it measured (which IS
  # the answer) looks like the symptom of one.
  diagnose_probe && (runtime_notice *= "<div class=\"alert info\"><strong>This is a diagnosis, not a solve.</strong> The self-check takes exactly one step from the case's own stored voltages (<code>max_iter = 1</code>, no rescue, Q-limit handling off), so a remaining residual is the RESULT and not a failure. The numbers below therefore show one iteration and a non-zero mismatch by design; <code>diagnose.log</code> names the worst bus and what to do about it.</div>")
  return _webui_layout(_webui_result_page_title(result), "<section class=\"panel\">$(result_summary)$(runtime_notice)<table class=\"details\">$(rows)</table>$(active_hint)$(links)$(se_chain_section)</section>$(n1_table_section)$(se_tap_section)$(se_topology_section)$(save_section)"; show_back = true, refresh_url)
end

# N-1 / scenario result table on the run result page (scenario task step 6,
# maintainer correction 2026-09-03): the browser must show WHICH case was
# screened, not only a count. Rows come from the run's contingency_n1.csv
# (14 classic columns, or 16 with the screening pair), ranked like
# printContingencyResults (failures first, then severity descending); the
# CSV artifact keeps input order and is linked as the complete list. A row
# with a screening estimate gets an always-visible detail line (no
# JavaScript): for a screened row the values ARE the one-step estimate,
# for a flagged row the estimate stands next to the full-run values.
# Renders for every case format; only the scenario editor is SCF-only.
const _WEBUI_CONTINGENCY_TABLE_MAX_ROWS = 100

function _webui_contingency_table_section(result::AbstractDict)::String
  metadata = get(result, "metadata", Dict{String,Any}())
  metadata isa AbstractDict || return ""
  get(metadata, "run_mode", "") == "contingency" || return ""
  od = String(get(result, "output_dir", ""))
  csv_path = isempty(od) ? "" : joinpath(od, "contingency_n1.csv")
  isfile(csv_path) || return ""
  lines = readlines(csv_path)
  length(lines) >= 2 || return ""
  has_screening = endswith(first(lines), ";screened;screening_estimate")
  esc = _webui_escape
  rows = NamedTuple[]
  for l in lines[2:end]
    f = split(l, ';')
    length(f) >= 14 || continue
    sev = tryparse(Float64, f[9])
    push!(rows, (
      name = String(f[1]), weight = String(f[2]), converged = String(f[3]) == "true",
      iterations = String(f[4]), start_used = String(f[5]),
      min_vm = String(f[6]), max_vm = String(f[7]), max_loading = String(f[8]),
      severity = sev === nothing ? NaN : sev,
      island_count = String(f[12]), shed = String(f[13]), error = String(f[14]),
      screened = has_screening && length(f) >= 16 && String(f[15]) == "true",
      estimate = has_screening && length(f) >= 16 ? String(f[16]) : "",
    ))
  end
  isempty(rows) && return ""
  # failures first, then severity descending, the printContingencyResults rank
  sort!(rows; by = r -> (isnan(r.severity) ? -Inf : -r.severity))
  shown = first(rows, _WEBUI_CONTINGENCY_TABLE_MAX_ROWS)
  fmtnum(s) = begin
    v = tryparse(Float64, s)
    v === nothing ? esc(s) : (isnan(v) ? "-" : string(round(v; digits = 4)))
  end
  ncols = 11 + (has_screening ? 1 : 0)
  body = join((begin
    detail = if !isempty(r.estimate)
      est = split(r.estimate, '|')
      est_text = length(est) == 3 ? "max loading $(fmtnum(String(est[1]))) %, vmin $(fmtnum(String(est[2]))) pu, vmax $(fmtnum(String(est[3]))) pu" : esc(r.estimate)
      note = r.screened ? "screened: the row values ARE this one-step estimate (no full solve)" : "flagged: estimate before the full run was $(est_text)"
      "<tr class=\"n1-detail\"><td></td><td colspan=\"$(ncols - 1)\">$(note)</td></tr>"
    else
      ""
    end
    string(
      "<tr class=\"", r.converged ? "n1-row" : "n1-row n1-failed", "\">",
      "<td>$(esc(r.name))</td><td>$(fmtnum(r.weight))</td><td>$(r.converged ? "yes" : "NO")</td><td>$(esc(r.iterations))</td><td>$(esc(r.start_used))</td>",
      "<td>$(fmtnum(r.min_vm))</td><td>$(fmtnum(r.max_vm))</td><td>$(fmtnum(r.max_loading))</td><td>$(fmtnum(string(r.severity)))</td>",
      "<td>$(esc(r.island_count))</td><td>$(fmtnum(r.shed))</td>",
      has_screening ? "<td class=\"n1-screened\">$(r.screened ? "yes" : "")</td>" : "",
      isempty(r.error) ? "" : "</tr><tr class=\"n1-detail n1-error\"><td></td><td colspan=\"$(ncols - 1)\">$(esc(r.error))</td>",
      "</tr>",
      detail,
    )
  end for r in shown), "")
  run_id = String(get(result, "run_id", ""))
  csv_link = isempty(run_id) ? "contingency_n1.csv" : "<a href=\"/powerflow/artifact/$(_webui_urlencode(run_id))/contingency_n1.csv\">contingency_n1.csv</a>"
  cap_note = length(rows) > length(shown) ? " Showing the $(length(shown)) worst of $(length(rows)) rows;" : ""
  header_cells = string(
    "<th>case</th><th>weight</th><th>converged</th><th>iter</th><th>start</th><th>Vmin [pu]</th><th>Vmax [pu]</th><th>loading [%]</th><th>severity</th><th>islands</th><th>shed [MW]</th>",
    has_screening ? "<th>screened</th>" : "",
  )
  return string(
    "<section class=\"panel n1-result-table\"><h2>N-1 / scenario results</h2>",
    "<p>Rows ranked by severity, failures first.$(cap_note) the complete list in input order is the CSV artifact $(csv_link).</p>",
    "<div class=\"table-scroll\"><table><thead><tr>", header_cells, "</tr></thead><tbody>", body, "</tbody></table></div></section>",
  )
end

function render_powerflow_artifacts(run_id::AbstractString, artifacts)::String
  if artifacts isa AbstractDict
    return render_powerflow_result(artifacts)
  end
  rows = join(
    (
      begin
        name = get(artifact, "name", "")
        href = "/powerflow/artifact/$(_webui_urlencode(run_id))/$(_webui_urlencode(name))"
        "<tr><td><a href=\"$(href)\">$(_webui_escape(name))</a></td><td>$(_webui_escape(get(artifact, "kind", "")))</td><td>$(_webui_escape(get(artifact, "mime_type", "")))</td><td>$(_webui_escape(get(artifact, "size_bytes", "")))</td><td>$(_webui_escape(get(artifact, "description", "")))</td></tr>"
      end for artifact in artifacts
    ),
    "",
  )
  table = "<section class=\"panel\"><p>Run <code>$(_webui_escape(run_id))</code></p><p><a class=\"button\" href=\"/powerflow/artifact-zip/$(_webui_urlencode(run_id))\">Download all artifacts as ZIP</a></p><table><thead><tr><th>Name</th><th>Kind</th><th>MIME type</th><th>Bytes</th><th>Description</th></tr></thead><tbody>$(rows)</tbody></table></section>"
  return _webui_layout("Artifacts", table; show_back = true)
end

"""
    _webui_is_completed_diagnose(run) -> Bool

A Diagnose run that did its job. The self-check takes exactly ONE step from the
case's own stored voltages (`max_iter = 1`, no rescue), so a remaining residual
is its RESULT, not its failure: that residual is the diagnosis. Reported as a
failed power flow it is unusable, because every diagnosis then looks like a
crash (reported from a live session 2026-09-08).

A diagnose run that could not run at all keeps the failure vocabulary; only the
"completed, did not converge" combination is a finished probe.
"""
function _webui_is_completed_diagnose(run::AbstractDict)::Bool
  # The history rows carry `run_mode` from the index; a result page does not,
  # so the self-check configuration the diagnose flow writes is the marker
  # there. Both paths have to agree, or the badge and the page contradict.
  is_diagnose = lowercase(string(get(run, "run_mode", ""))) == "diagnose"
  if !is_diagnose
    artifacts = get(run, "artifacts", nothing)
    is_diagnose = artifacts isa AbstractVector &&
                  any(a -> a isa AbstractDict && String(get(a, "name", "")) == "diagnose_self_check_config.yaml", artifacts)
  end
  is_diagnose || return false
  status = lowercase(string(get(run, "status", "")))
  run_status = lowercase(string(get(run, "run_status", "")))
  return status == "not_converged" || run_status == "completed_nonconverged"
end

"""
    _webui_status_badge(css_class, label, raw_status) -> String

A status badge. The tooltip carries the raw status ONLY where the label says
something else ("diagnosed" for a `not_converged` self-check); a badge whose
label already is the status gets no tooltip, because repeating the word on
hover tells the reader nothing.
"""
function _webui_status_badge(css_class::AbstractString, label::AbstractString, raw_status::AbstractString)::String
  title = label == raw_status ? "" : " title=\"$(_webui_escape(raw_status))\""
  return "<span class=\"status-badge $(css_class)\"$(title)>$(_webui_escape(label))</span>"
end

function webui_status_class(run::AbstractDict)::String
  status = lowercase(string(get(run, "status", "unknown")))
  success = get(run, "success", nothing)
  _webui_is_completed_diagnose(run) && return "status-info"
  status in ("running", "pending", "queued", "aborting") && return "status-running"
  status in ("aborted", "aborted_unknown", "interrupted", "cancelled", "canceled") && return "status-aborted"
  status in ("warning", "partial", "questionable") && return "status-warning"
  status in ("failed", "failure", "error", "not_converged") && return "status-error"
  get(run, "numerical_status", nothing) == "not_converged" && return "status-error"
  get(run, "converged", nothing) === false && return "status-error"
  status in ("succeeded", "success", "converged", "ok") && return success === false ? "status-error" : "status-success"
  success === true && return "status-success"
  success === false && return "status-error"
  return "status-unknown"
end

function _webui_run_timestamp(run::AbstractDict)::String
  timestamp = get(run, "timestamp", "Unknown")
  return isempty(strip(string(timestamp))) ? "Unknown" : string(timestamp)
end

"""
    _webui_comparable_kind(kind) -> Bool

Whether a run of this kind can go into the side-by-side comparison. The page
reads power-flow artifacts (effective configuration, Q-limit events, branch
flows, bus voltages), so a plain power flow, a diagnose probe and a power flow
started from a state estimate qualify; a state estimation, short circuit, N-1
or import analysis has nothing the page would read. Deliberately NOT decided
by network size: the same case modelled with an external-grid source and with
a slack is exactly the pair one wants side by side (maintainer 2026-09-09).
"""
_webui_comparable_kind(kind::AbstractString)::Bool = lowercase(strip(kind)) in ("", "powerflow", "diagnose", "powerflow_se_start")

"A history row that may be ticked for comparison: comparable kind, finished, and its result still on disk."
function _webui_run_comparable(run::AbstractDict)::Bool
  get(run, "available", false) == true || return false
  _webui_comparable_kind(string(get(run, "run_mode", ""))) || return false
  return !(lowercase(string(get(run, "status", ""))) in _WEBUI_ACTIVE_RUN_STATUSES)
end

function render_powerflow_history(runs, output_root::AbstractString; active_run = nothing)::String
  ordered_runs = sort!(collect(runs); by = run -> _webui_run_timestamp(run), rev = true)
  rows = join((begin
    run_id = string(get(run, "run_id", ""))
    available = get(run, "available", false)
    # Two runs of the same case under different settings are the normal way to
    # look at a Q-limit mode or a solver choice; the history could list them
    # but not put them side by side. The checkbox feeds /powerflow/compare and
    # is offered only where the comparison has something to read.
    pick = _webui_run_comparable(run) ? "<td><input type=\"checkbox\" name=\"run\" value=\"$(_webui_escape(run_id))\" aria-label=\"Select run $(_webui_escape(run_id)) for comparison\"></td>" : "<td></td>"
    link = available ? "<a href=\"/powerflow/result/$(_webui_urlencode(run_id))\">$(_webui_escape(run_id))</a>" : _webui_escape(run_id)
    status = _webui_is_completed_diagnose(run) ? "diagnosed" : string(get(run, "status", "unknown"))
    status_badge = _webui_status_badge(webui_status_class(run), status, string(get(run, "status", "")))
    delete_form = "<form method=\"post\" action=\"/powerflow/delete/$(_webui_urlencode(run_id))\" class=\"delete-run-form\"><button type=\"submit\" class=\"danger-button\">Delete</button></form>"
    abort_form = lowercase(status) in ("queued", "running") ? "<form method=\"post\" action=\"/powerflow/abort/$(_webui_urlencode(run_id))\"><button type=\"submit\" class=\"danger-button\">Abort</button></form>" : ""
    # run kind (SE phase 5): "" = plain power flow, also what pre-0.10 index
    # entries without the field fall back to
    kind = string(get(run, "run_mode", ""))
    kind_label = isempty(kind) ? "powerflow" : kind
    fields = (_webui_run_timestamp(run), link, status_badge, kind_label, available, _webui_run_method(run, String(get(run, "solver", ""))), get(run, "iterations", ""), get(run, "final_mismatch", ""))
    # The two paths no longer fit the row. The name says which case and
    # which configuration; the tooltip keeps the path for whoever needs it
    # (maintainer 2026-09-09).
    path_cell = path -> "<td title=\"$(_webui_escape(string(path)))\">$(_webui_escape(basename(string(path))))</td>"
    cells = "<td>$(_webui_escape(fields[1]))</td><td>$(fields[2])</td><td>$(fields[3])</td>" * join(("<td>$(_webui_escape(field))</td>" for field in fields[4:end]), "") * path_cell(get(run, "casefile", "")) * path_cell(get(run, "config_file", ""))
    "<tr>$(pick)$(cells)<td>$(abort_form)$(delete_form)</td></tr>"
  end for run in ordered_runs), "")
  # The compare button submits the surrounding form; the handler is the place
  # that insists on exactly two runs, so the page stays usable without JS.
  compare_hint = "<p class=\"field-hint\">Tick two runs and press Compare to see them side by side: status, iterations, residual, the configuration keys that differ, the Q-limit events and, where both runs wrote bus voltages, the largest voltage deviation.</p>"
  content = "$(_webui_active_run_banner(active_run))<section class=\"panel history-actions\"><p><strong>Output root:</strong> <code>$(_webui_escape(output_root))</code></p><div class=\"actions\"><form method=\"post\" action=\"/powerflow/refresh\"><button type=\"submit\">Refresh registry</button></form><form method=\"post\" action=\"/powerflow/delete_all\"><button type=\"submit\" class=\"danger-button\">Delete all runs</button></form></div></section>\n<section class=\"panel\"><form method=\"get\" action=\"/powerflow/compare\">$(compare_hint)<div class=\"actions\"><button type=\"submit\" class=\"secondary-button\">Compare selected runs</button></div><table><thead><tr><th>Compare</th><th>Date/Time</th><th>Run ID</th><th>Status</th><th>Kind</th><th>Available</th><th>Solver</th><th>Iterations</th><th>Final mismatch</th><th>Case file</th><th>Config file</th><th>Delete</th></tr></thead><tbody>$(rows)</tbody></table></form></section>"
  return _webui_layout("Run history", content; show_back = true)
end

"""
    _webui_compare_config_diff(a_dir, b_dir) -> Vector{Tuple{String,String,String}}

The configuration keys two runs disagree on, read from the `effective_config.yaml`
each run writes. Line based on purpose: the file is the record of what a run
actually used, and a line diff needs no schema knowledge and cannot go stale
when the configuration grows a section.
"""
function _webui_compare_config_diff(a_dir::AbstractString, b_dir::AbstractString)
  read_keys(dir) = begin
    path = joinpath(dir, "effective_config.yaml")
    out = Dict{String,String}()
    isfile(path) || return out
    # indent stack, so a nested key gets its real dotted path: a plain
    # "remember the last section" would concatenate every section it ever saw
    stack = Tuple{Int,String}[]
    for line in eachline(path)
      isempty(strip(line)) && continue
      startswith(strip(line), "#") && continue
      indent = length(line) - length(lstrip(line))
      m = match(r"^\s*([A-Za-z0-9_.]+):\s*(.*)$", line)
      m === nothing && continue
      key, value = String(m.captures[1]), strip(String(m.captures[2]))
      while !isempty(stack) && stack[end][1] >= indent
        pop!(stack)
      end
      prefix = isempty(stack) ? "" : string(join((s[2] for s in stack), "."), ".")
      if isempty(value)
        push!(stack, (indent, key))
        continue
      end
      out[string(prefix, key)] = value
    end
    return out
  end
  a = read_keys(a_dir)
  b = read_keys(b_dir)
  rows = Tuple{String,String,String}[]
  for key in sort(collect(union(keys(a), keys(b))))
    # `_config_sources` records where each value came from. It mirrors every
    # difference a second time, which reads like two findings instead of one.
    startswith(key, "_config_sources.") && continue
    va = get(a, key, "-")
    vb = get(b, key, "-")
    va == vb || push!(rows, (key, va, vb))
  end
  return rows
end

"""
    _webui_compare_qlimit_buses(dir) -> Vector{Int} or nothing

Clamped buses of a run, from the `q_limit_events.csv` it wrote. The file
follows the run's CSV format like every other artifact (its writer takes the
same `format` argument), so it is read by column name rather than by position.
"""
function _webui_compare_qlimit_buses(dir::AbstractString)
  table = _webui_compare_csv_table(joinpath(dir, "q_limit_events.csv"))
  table === nothing && return nothing
  buses = Int[]
  for row in table.rows
    bus = tryparse(Int, strip(get(row, "bus", "")))
    bus === nothing || push!(buses, bus)
  end
  return sort!(unique!(buses))
end

"""
    _webui_compare_split_csv(line, delimiter) -> Vector{String}

One CSV line into its fields, honoring quoted fields. The detailed export
quotes any cell that contains the delimiter, which the `excel_us` format
produces for every grouped number (`"1,234.5"` in a comma-delimited file).
"""
function _webui_compare_split_csv(line::AbstractString, delimiter::Char)::Vector{String}
  fields = String[]
  buffer = IOBuffer()
  in_quotes = false
  i = firstindex(line)
  while i <= lastindex(line)
    c = line[i]
    if in_quotes
      if c == '"'
        # a doubled quote inside a quoted field is one literal quote
        if i < lastindex(line) && line[nextind(line, i)] == '"'
          print(buffer, '"')
          i = nextind(line, i)
        else
          in_quotes = false
        end
      else
        print(buffer, c)
      end
    elseif c == '"'
      in_quotes = true
    elseif c == delimiter
      push!(fields, String(take!(buffer)))
    else
      print(buffer, c)
    end
    i = nextind(line, i)
  end
  push!(fields, String(take!(buffer)))
  return fields
end

"""
    _webui_compare_csv_table(path) -> NamedTuple or nothing

A detailed-export CSV as named rows, in whichever of the three formats the run
was written with. The delimiter comes from the header line and settles the
number format with it (`run_api.jl`: `technical` = `,` and `.`; `excel_de` =
`;`, `,` decimals, `.` grouping; `excel_us` = `,`, `.` decimals, `,` grouping),
so one rule reads all three. Getting this wrong is not loud: a comma-splitting
reader on a semicolon file finds no columns at all and the page then claims the
two runs have nothing in common (reported 2026-09-08).
"""
function _webui_compare_csv_table(path::AbstractString)
  isfile(path) || return nothing
  lines = readlines(path)
  isempty(lines) && return nothing
  head = strip(lines[1])
  delimiter = count(==(';'), head) > count(==(','), head) ? ';' : ','
  decimal = delimiter == ';' ? ',' : '.'
  thousands = delimiter == ';' ? '.' : ','
  header = _webui_compare_split_csv(head, delimiter)
  rows = Vector{Dict{String,String}}()
  for line in Iterators.drop(lines, 1)
    isempty(strip(line)) && continue
    parts = _webui_compare_split_csv(strip(line), delimiter)
    length(parts) == length(header) || continue
    push!(rows, Dict(header[k] => parts[k] for k in eachindex(header)))
  end
  return (rows = rows, decimal = decimal, thousands = thousands)
end

"A number as the detailed export wrote it, in the format that table uses."
function _webui_compare_number(text::AbstractString, decimal::Char, thousands::Char)
  cleaned = replace(strip(text), string(thousands) => "")
  decimal == '.' || (cleaned = replace(cleaned, decimal => '.'))
  return tryparse(Float64, cleaned)
end

"""
    _webui_compare_voltages(dir) -> Dict or nothing

Bus voltages of a run from `bus_voltages_complex.csv` (detailed CSV export),
keyed by the bus index, which is what two runs of the same case share. The bus
name travels along for the table, because an index alone says nothing to
someone reading the page.
"""
function _webui_compare_voltages(dir::AbstractString)
  table = _webui_compare_csv_table(joinpath(dir, "bus_voltages_complex.csv"))
  table === nothing && return nothing
  out = Dict{String,NamedTuple{(:vm, :va, :name),Tuple{Float64,Float64,String}}}()
  for row in table.rows
    bus = strip(get(row, "bus", ""))
    vm = _webui_compare_number(get(row, "vm_pu", ""), table.decimal, table.thousands)
    va = _webui_compare_number(get(row, "va_deg", ""), table.decimal, table.thousands)
    (isempty(bus) || vm === nothing || va === nothing) && continue
    # the case's own name if it has one, otherwise the index it is keyed by
    name = strip(get(row, "original_bus_name", ""))
    isempty(name) && (name = strip(get(row, "bus_name", "")))
    out[String(bus)] = (vm = vm, va = va, name = isempty(name) ? String(bus) : String(name))
  end
  return isempty(out) ? nothing : out
end

"""
    _webui_compare_losses(dir) -> NamedTuple or nothing

Total branch losses of a run, summed over the `p_loss_MW`/`q_loss_MVar` columns
the run wrote into `branch_flows.csv`. The sum is the only arithmetic this page
does; every number in it is otherwise read as written.
"""
function _webui_compare_losses(dir::AbstractString)
  table = _webui_compare_csv_table(joinpath(dir, "branch_flows.csv"))
  table === nothing && return nothing
  p = 0.0
  q = 0.0
  n = 0
  for row in table.rows
    dp = _webui_compare_number(get(row, "p_loss_MW", ""), table.decimal, table.thousands)
    dq = _webui_compare_number(get(row, "q_loss_MVar", ""), table.decimal, table.thousands)
    dp === nothing && continue
    p += dp
    q += dq === nothing ? 0.0 : dq
    n += 1
  end
  return n == 0 ? nothing : (p_MW = p, q_MVAr = q, branches = n)
end

"""
    render_powerflow_compare(a, b) -> String

Two finished runs side by side. Everything shown comes from what the runs
themselves wrote; nothing is recomputed, so the page also works for runs from
an earlier session.
"""
function render_powerflow_compare(a::AbstractDict, b::AbstractDict)::String
  id(r) = string(get(r, "run_id", ""))
  cell(r, key, default = "-") = _webui_escape(string(get(r, key, default)))
  head = string(
    "<section class=\"panel\"><h2>Two runs side by side</h2><table><thead><tr><th>Property</th>",
    "<th>A: $(cell(a, "run_id"))</th><th>B: $(cell(b, "run_id"))</th></tr></thead><tbody>",
  )
  for (label, key) in (("Case file", "casefile"), ("Status", "status"), ("Converged", "converged"),
    ("Iterations", "iterations"), ("Final mismatch", "final_mismatch"), ("Config file", "config_file"))
    head *= "<tr><td>$(_webui_escape(label))</td><td>$(cell(a, key))</td><td>$(cell(b, key))</td></tr>"
  end
  head *= "</tbody></table><p class=\"actions\"><a href=\"/powerflow/result/$(_webui_urlencode(id(a)))\">Open A</a> &middot; <a href=\"/powerflow/result/$(_webui_urlencode(id(b)))\">Open B</a></p></section>"

  dir_a = String(get(a, "output_dir", ""))
  dir_b = String(get(b, "output_dir", ""))

  diff = _webui_compare_config_diff(dir_a, dir_b)
  cfg_section = if isempty(diff)
    "<section class=\"panel\"><h2>Configuration</h2><p>Both runs used the same effective configuration.</p></section>"
  else
    rows = join(("<tr><td><code>$(_webui_escape(k))</code></td><td>$(_webui_escape(va))</td><td>$(_webui_escape(vb))</td></tr>" for (k, va, vb) in diff), "")
    "<section class=\"panel\"><h2>Configuration differences ($(length(diff)))</h2><table><thead><tr><th>Key</th><th>A</th><th>B</th></tr></thead><tbody>$(rows)</tbody></table></section>"
  end

  qa = _webui_compare_qlimit_buses(dir_a)
  qb = _webui_compare_qlimit_buses(dir_b)
  q_section = if qa === nothing && qb === nothing
    "<section class=\"panel\"><h2>Q-limit events</h2><p>Neither run recorded Q-limit events.</p></section>"
  else
    only_a = setdiff(something(qa, Int[]), something(qb, Int[]))
    only_b = setdiff(something(qb, Int[]), something(qa, Int[]))
    same = isempty(only_a) && isempty(only_b)
    note = same ? "<p>Both runs clamped the same buses.</p>" :
           "<p><strong>The runs clamped different buses.</strong> Only in A: $(_webui_escape(string(only_a))). Only in B: $(_webui_escape(string(only_b))).</p>"
    "<section class=\"panel\"><h2>Q-limit events</h2>$(note)<table><thead><tr><th></th><th>Clamped buses</th></tr></thead><tbody><tr><td>A</td><td>$(_webui_escape(string(something(qa, "no file"))))</td></tr><tr><td>B</td><td>$(_webui_escape(string(something(qb, "no file"))))</td></tr></tbody></table></section>"
  end

  la = _webui_compare_losses(dir_a)
  lb = _webui_compare_losses(dir_b)
  loss_section = if la === nothing || lb === nothing
    "<section class=\"panel\"><h2>Losses</h2><p>At least one run has no <code>branch_flows.csv</code>; it is written when the detailed result CSV export is on.</p></section>"
  else
    dp = la.p_MW - lb.p_MW
    dq = la.q_MVAr - lb.q_MVAr
    # a difference far below the convergence tolerance is the same operating
    # point written twice, not a finding
    verdict = abs(dp) < 1e-6 ? "<p>Both runs end on the same losses.</p>" :
              "<p><strong>The runs differ by $(round(dp; sigdigits = 4)) MW</strong> of active losses ($(round(100 * dp / max(abs(lb.p_MW), eps()); sigdigits = 3)) % of B).</p>"
    body = string(
      "<tr><td>Active losses [MW]</td><td>$(round(la.p_MW; digits = 4))</td><td>$(round(lb.p_MW; digits = 4))</td><td>$(round(dp; sigdigits = 4))</td></tr>",
      "<tr><td>Reactive losses [MVAr]</td><td>$(round(la.q_MVAr; digits = 4))</td><td>$(round(lb.q_MVAr; digits = 4))</td><td>$(round(dq; sigdigits = 4))</td></tr>",
      "<tr><td>Branches summed</td><td>$(la.branches)</td><td>$(lb.branches)</td><td>$(la.branches - lb.branches)</td></tr>",
    )
    "<section class=\"panel\"><h2>Losses</h2>$(verdict)<table><thead><tr><th>Quantity</th><th>A</th><th>B</th><th>A - B</th></tr></thead><tbody>$(body)</tbody></table></section>"
  end

  va = _webui_compare_voltages(dir_a)
  vb = _webui_compare_voltages(dir_b)
  v_section = if va === nothing || vb === nothing
    "<section class=\"panel\"><h2>Bus voltages</h2><p>At least one run has no readable <code>bus_voltages_complex.csv</code>; it is written when the detailed result CSV export is on.</p></section>"
  else
    shared = sort(collect(intersect(keys(va), keys(vb))); by = k -> something(tryparse(Int, k), typemax(Int)))
    if isempty(shared)
      "<section class=\"panel\"><h2>Bus voltages</h2><p>The two runs have no bus in common, so they are not two runs of the same network.</p></section>"
    else
      deltas = [(bus, va[bus].vm - vb[bus].vm, va[bus].va - vb[bus].va) for bus in shared]
      sort!(deltas; by = t -> -abs(t[2]))
      maxdvm = maximum(abs(t[2]) for t in deltas)
      maxdva = maximum(abs(t[3]) for t in deltas)
      if maxdvm < 1e-9 && maxdva < 1e-9
        # a table of zeros over every bus hides the one thing it says
        "<section class=\"panel\"><h2>Bus voltages</h2><p>Both runs end on the same voltages at all $(length(shared)) shared buses.</p></section>"
      else
        rows = join(("<tr><td>$(_webui_escape(bus))</td><td>$(_webui_escape(va[bus].name))</td><td>$(round(va[bus].vm; digits = 6))</td><td>$(round(vb[bus].vm; digits = 6))</td><td>$(round(dvm; sigdigits = 4))</td><td>$(round(dva; sigdigits = 4))</td></tr>" for (bus, dvm, dva) in first(deltas, 10)), "")
        summary = "<p>max |dVm| = <strong>$(round(maxdvm; sigdigits = 4)) pu</strong>, max |dVa| = $(round(maxdva; sigdigits = 4)) deg, over $(length(shared)) shared buses. The ten largest deviations:</p>"
        "<section class=\"panel\"><h2>Bus voltages</h2>$(summary)<table><thead><tr><th>Bus</th><th>Name</th><th>A Vm [pu]</th><th>B Vm [pu]</th><th>dVm</th><th>dVa [deg]</th></tr></thead><tbody>$(rows)</tbody></table></section>"
      end
    end
  end

  return _webui_layout("Compare runs", string(head, cfg_section, q_section, loss_section, v_section); show_back = true)
end

function render_webui_operation_log(content::AbstractString; entries::Integer = 0, bytes::Integer = 0)::String
  # the size is what makes this page unwieldy, so it is stated next to the
  # controls: entries first, because the retention works on entries, not bytes
  size_note = "<p class=\"field-hint\">$(entries) entries, $(round(bytes / 1024; digits = 1)) kB. Entries older than <code>webui.operation_log_retention_days</code> are dropped at every start.</p>"
  controls = string(
    "<p class=\"actions\"><a class=\"button\" href=\"/webui/operation-log/download\">Download operation log</a>",
    "<form method=\"post\" action=\"/webui/operation-log/clear\" onsubmit=\"return confirm('Delete the whole operation log? The file is emptied and a single entry records the deletion.');\">",
    "<button type=\"submit\" class=\"secondary-button\" title=\"Empty the operation log. The file stays in place and records that it was cleared, so the next entries have a starting point.\">Clear operation log</button></form></p>",
    size_note,
  )
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

## Section index of a documentation page (maintainer 2026-09-04: the
## reference pages run to a thousand lines and nobody scrolls that). Built
## from the level-2 headings, so a reader jumps instead of scrolls; the
## anchors are the ones Julia's Markdown renderer emits for headings.
function _webui_doc_page_toc(markdown_text::AbstractString)::String
  entries = String[]
  in_fence = false
  for line in split(String(markdown_text), '\n')
    startswith(line, "```") && (in_fence = !in_fence)
    in_fence && continue
    startswith(line, "## ") || continue
    title = strip(line[4:end])
    isempty(title) && continue
    anchor = lowercase(replace(title, r"[^\w\s-]" => "", r"\s+" => "-"))
    push!(entries, "<li><a href=\"#$(_webui_escape(anchor))\">$(_webui_escape(title))</a></li>")
  end
  length(entries) < 4 && return ""
  return "<details class=\"docs-toc\" open><summary>On this page ($(length(entries)) sections)</summary><ul>$(join(entries, ""))</ul></details>"
end

function render_webui_doc_page(page::AbstractString, metadata, markdown_text::AbstractString)::String
  content = "<section class=\"panel docs-page docs-content\">$(_webui_doc_page_toc(markdown_text))$(render_webui_markdown(markdown_text; current_page = page))</section>"
  return _webui_layout(metadata.title, content; show_back = true)
end

# ---------------------------------------------------------------------------
# State estimation page and result pieces (SE phase 5)
# ---------------------------------------------------------------------------

## result-summary card for an SE run (metadata run_mode == "se")
# auto power-flow mode card: profile, escalation depth, final solver, the
# advisory hint list, and a loud DC-fallback label (never mistakable for AC)
function _webui_auto_mode_summary(result::AbstractDict)::Union{Nothing,String}
  metadata = get(result, "metadata", nothing)
  metadata isa AbstractDict || return nothing
  get(metadata, "auto_mode_enabled", false) === true || return nothing
  head = "profile $(get(metadata, "auto_profile", "?")), $(get(metadata, "auto_escalation_stages", 0)) escalation stage(s), final stage $(get(metadata, "auto_final_stage", "?")), solver $(get(metadata, "auto_final_solver", "?"))"
  s = _webui_escape(head)
  get(metadata, "dc_fallback_solution", false) === true && (s *= "<div class=\"alert warning\">DC fallback approximation; the AC problem did not converge.</div>")
  hints = get(metadata, "auto_hints", String[])
  if hints isa AbstractVector && !isempty(hints)
    s *= "<ul class=\"auto-hints\">" * join(("<li>$(_webui_escape(String(h)))</li>" for h in hints), "") * "</ul>"
  end
  return s
end

"""
    _webui_control_summary(result) -> String or nothing

The controllers the solved network carried, as one line for the summary cards:
tap changers, Q(U) and P(U) machines. Nothing when the run had none, so an
ordinary case keeps its summary short. The counts come from the run's own
metadata (`controllers`), written by the service.
"""
function _webui_control_summary(result::AbstractDict)::Union{Nothing,String}
  metadata = get(result, "metadata", Dict{String,Any}())
  metadata isa AbstractDict || return nothing
  counts = get(metadata, "controllers", nothing)
  counts isa AbstractDict || return nothing
  n(key) = something(tryparse(Int, string(get(counts, key, 0))), 0)
  tap, qu, pu = n("tap"), n("qu"), n("pu")
  tap + qu + pu > 0 || return nothing
  parts = String[]
  qu > 0 && push!(parts, "Q(U) $(qu)")
  pu > 0 && push!(parts, "P(U) $(pu)")
  tap > 0 && push!(parts, "tap $(tap)")
  return join(parts, " · ")
end

function _webui_se_summary(result::AbstractDict)::Union{Nothing,String}
  metadata = get(result, "metadata", Dict{String,Any}())
  metadata isa AbstractDict || return nothing
  get(metadata, "run_mode", "") == "se" || return nothing
  quality = String(get(metadata, "se_observability_quality", "?"))
  light = quality == "good" ? "🟢" : quality == "critical" ? "🟡" : "🔴"
  reason = String(get(metadata, "se_band_reason", "?"))
  parts = String[]
  push!(parts, "$(light) observability $(quality)$(get(metadata, "se_structural_islands", false) ? " (structural islands)" : "")")
  # chi-square plausibility at a glance: E[J] = dof for healthy noise, and
  # the Wilson-Hilferty 3-sigma band verdict with its reason. The two
  # failure directions read very differently and say so explicitly: :high
  # is the alarm (bad data or model mismatch), :low means J is far BELOW
  # the expectation (sigmas overstate the errors, the normal signature of
  # a noise-free synthetic set) and is NOT an alarm.
  within = get(metadata, "se_j_within_3sigma", nothing)
  bandtxt = if within === nothing
    "band $(reason)"
  elseif Bool(within)
    "J within 3σ band ✓"
  elseif reason == "low"
    "J far BELOW the expected value (sigmas overstate the errors; normal for a noise-free set, not an alarm)"
  else
    "J OUTSIDE 3σ band ($(reason))"
  end
  # J/dof leads: J alone grows with the number of measurements, so the same
  # healthy set reads as an alarm on a larger case. The ratio is the number a
  # reader can judge at a glance.
  se_j = Float64(get(metadata, "se_objective", NaN))
  se_dof = get(metadata, "se_dof", "?")
  jratio = se_dof isa Number && se_dof > 0 ? string(round(se_j / se_dof; digits = 2)) : "n/a"
  push!(parts, "$(get(metadata, "se_iterations", "?")) iteration(s), J/dof = $(jratio) (J = $(round(se_j; digits = 3)), dof = $(se_dof)), $(bandtxt)")
  # a set that measures the same quantity twice inflates J on its own; without
  # this note the page shows an alarming J and no cause
  se_dup = get(metadata, "se_duplicate_rows", 0)
  se_dup isa Number && se_dup > 0 && push!(parts, "$(Int(se_dup)) measurement(s) repeat an already measured quantity (regenerate the set: this alone inflates J)")
  # J_active: present only when replacement suppression removed rows from
  # the state; the band verdict above stays on the honest J
  se_ja = get(metadata, "se_objective_active", nothing)
  se_ja isa Number && push!(parts, "J_active = $(round(Float64(se_ja); digits = 3)) (dof = $(get(metadata, "se_dof_active", "?")), without $(get(metadata, "se_suppressed_rows", "?")) suppressed row(s))")
  push!(parts, "$(get(metadata, "se_suspicious", 0)) suspicious, $(get(metadata, "se_eliminations", 0)) eliminated")
  se_rm = String(get(metadata, "se_robust_mode", get(metadata, "se_robust", false) ? "staged" : "off"))
  se_rm == "staged" && push!(parts, "robust R modification on (staged, k1=$(get(metadata, "se_robust_k1", 3.0)), k2=$(get(metadata, "se_robust_k2", 6.0)))")
  # the 6.0 fallback is HISTORY, not a default: runs from before the
  # metadata key existed were solved with the service literal 6.0, so that
  # is the honest value to show for them. Do not "correct" it to the
  # current 4.0, that would misreport old runs.
  se_rm == "replacement" && push!(parts, "bad-data suppression on (replacement, k_suppress=$(get(metadata, "se_k_suppress", 6.0)), sigma=$(get(metadata, "se_suppression_sigma", 2000.0)))")
  # released taps: J of the continuous estimate versus J after the
  # mandatory fixation to the mechanical step (the tap left the state there)
  # The tap fallback must be impossible to miss: a run whose tap positions
  # are MODEL values looks exactly like a successful tap estimation
  # otherwise, and its J measures the model positions. Marked line, not a
  # footnote (maintainer, task_se_tap_bounds_v0100).
  if get(metadata, "se_tap_estimation_fallback", false) == true
    push!(parts, string("WARNING: ", _SE_TAP_FALLBACK_NOTE))
  end
  tapc = Int(get(metadata, "se_tap_count", 0))
  if tapc > 0
    jb = get(metadata, "se_tap_j_before", nothing)
    ja = get(metadata, "se_tap_j_after", nothing)
    fixedtxt = get(metadata, "se_tap_fixed", false) == true ? "" : " (NOT fixed)"
    push!(parts, "tap estimation: $(tapc) transformer(s), J $(jb === nothing ? "?" : round(Float64(jb); sigdigits = 3)) before fixation -> $(ja === nothing ? "?" : round(Float64(ja); sigdigits = 3)) after$(fixedtxt)")
    # off-grid interpretation: a band failure that only
    # appears through the fixation must not be read as bad data
    if get(metadata, "se_tap_offgrid_residual", false) == true
      push!(parts, "off-grid tap residual: the J jump comes from rounding to the mechanical step (true position between steps, or wrong step table), NOT from bad data")
    end
  end
  return _webui_escape(join(parts, " | "))
end

## topology panel on the SE result page (all advisory): the stage-1
## precheck findings, the stage-2 suspected stations, the explicit
## hypothesis-test button (stage 3 never runs automatically), and the
## ranked recommendation table once topology_hypotheses.md exists
function _webui_se_topology_section(result::AbstractDict)::String
  metadata = get(result, "metadata", Dict{String,Any}())
  metadata isa AbstractDict || return ""
  get(metadata, "run_mode", "") == "se" || return ""
  esc = _webui_escape
  pre = get(metadata, "se_topology_findings", nothing)
  stations = get(metadata, "se_topology_station_findings", nothing)
  run_id = String(get(result, "run_id", ""))
  outdir = String(get(result, "output_dir", ""))
  hypfile = isempty(outdir) ? "" : joinpath(outdir, "topology_hypotheses.md")
  havePre = pre isa AbstractVector && !isempty(pre)
  haveSt = stations isa AbstractVector && !isempty(stations)
  haveHyp = !isempty(hypfile) && isfile(hypfile)
  showButton = get(result, "success", false) == true && !isempty(run_id)
  (havePre || haveSt || haveHyp || showButton) || return ""
  parts = String["<section class=\"panel\"><h2>Topology validation (advisory)</h2>"]
  if havePre
    rows = join(("<tr><td>$(esc(String(get(f, "severity", ""))))</td><td><code>$(esc(String(get(f, "kind", ""))))</code></td><td>$(esc(String(get(f, "location", ""))))</td><td>$(esc(String(get(f, "evidence", ""))))</td></tr>" for f in pre), "")
    push!(parts, "<h3>Pre-check findings</h3><div style=\"overflow-x:auto\"><table class=\"details\"><thead><tr><th>severity</th><th>kind</th><th>location</th><th>evidence</th></tr></thead><tbody>$(rows)</tbody></table></div>")
  elseif get(metadata, "se_topology_precheck_summary", "") != ""
    push!(parts, "<p>$(esc(String(get(metadata, "se_topology_precheck_summary", ""))))</p>")
  end
  if haveSt
    rows = join(("<tr><td>$(esc(String(get(f, "location", ""))))</td><td>$(esc(String(get(f, "evidence", ""))))</td><td>$(esc(join(get(f, "notes", String[]), ", ")))</td></tr>" for f in stations), "")
    push!(parts, "<h3>Suspected topology error</h3><p>The bad-data elimination exhausted while the band stayed high and the suspects cluster at one station: this is the fingerprint of a wrong service state, not of bad telemetry.</p><div style=\"overflow-x:auto\"><table class=\"details\"><thead><tr><th>station</th><th>clustered suspects</th><th>notes</th></tr></thead><tbody>$(rows)</tbody></table></div>")
  end
  if haveHyp
    # render the markdown table of the artifact inline (simple pipe parse)
    lines = readlines(hypfile)
    trows = String[]
    thead = ""
    for l in lines
      startswith(l, "|") || continue
      cells = [strip(c) for c in split(l, "|")[2:(end - 1)]]
      all(c -> all(x -> x in ('-', ':', ' '), c), cells) && continue
      if isempty(thead)
        thead = join(("<th>$(esc(String(c)))</th>" for c in cells), "")
      else
        push!(trows, string("<tr>", join(("<td>$(esc(String(c)))</td>" for c in cells), ""), "</tr>"))
      end
    end
    isempty(thead) || push!(parts, "<h3>Hypothesis test (ranked recommendations)</h3><div style=\"overflow-x:auto\"><table class=\"details\"><thead><tr>$(thead)</tr></thead><tbody>$(join(trows, ""))</tbody></table></div><p class=\"field-help\">Recommendations only: nothing has been switched. Full report: topology_hypotheses.md in the artifacts.</p>")
  end
  if showButton
    push!(parts, "<form method=\"post\" action=\"/stateestimation/topology-hypotheses\" data-busy=\"Testing hypotheses, this runs several estimations…\"><input type=\"hidden\" name=\"run_id\" value=\"$(esc(run_id))\"><button type=\"submit\">Test topology hypotheses</button> <span class=\"field-help\">Re-runs the estimation per candidate status toggle on working copies (up to 5 candidates); switches NOTHING.</span></form>")
  end
  push!(parts, "</section>")
  return join(parts, "")
end

## released-tap estimate table on the SE result page: steps on the
## mechanical grid (never raw r values), the fixation J drop, mRIDs for
## CGMES cases (MATPOWER/DTF rows are addressed by branch and name)
function _webui_se_tap_section(result::AbstractDict)::String
  metadata = get(result, "metadata", Dict{String,Any}())
  metadata isa AbstractDict || return ""
  get(metadata, "run_mode", "") == "se" || return ""
  rows = get(metadata, "se_tap_estimates", nothing)
  esc = _webui_escape
  # the fallback replaces the table: showing estimated-looking positions
  # that are model values is worse than showing none
  if get(metadata, "se_tap_estimation_fallback", false) == true
    return string("<section class=\"result-block\"><h3>Transformer tap estimation</h3>",
      "<p class=\"warning\">", esc(_SE_TAP_FALLBACK_NOTE), "</p></section>")
  end
  (rows isa AbstractVector && !isempty(rows)) || return ""
  anymrid = any(!isempty(String(get(t, "mrid", ""))) for t in rows)
  anypst = any(String(get(t, "mode", "ratio")) != "ratio" for t in rows)
  head = string("<tr><th>Branch</th><th>Transformer</th>", anymrid ? "<th>mRID</th>" : "", "<th>Mode</th><th>Electrical step</th><th>Fixed step</th>", anypst ? "<th>Shift step (el.)</th><th>Shift step (fixed)</th>" : "", "<th>Status</th></tr>")
  body = String[]
  for t in rows
    mode = String(get(t, "mode", "ratio"))
    fr = String(get(t, "frozen_reason", "none"))
    # machine transformers are never estimated: their tap is no state
    # variable, the position is back-calculated after the run and marked so
    status = if String(get(t, "source", "estimated")) == "calculated"
      "calculated (machine transformer, not estimated)"
    elseif !(fr in ("none", ""))
      "frozen ($(replace(fr, "_" => " ")))"
    elseif get(t, "out_of_range", false) == true
      "out of range (clamped)"
    elseif get(t, "fixed", false) == true
      "fixed"
    else
      "not fixed"
    end
    push!(body, string(
      "<tr><td>", get(t, "branch", "?"), "</td><td>", esc(String(get(t, "name", ""))), "</td>",
      anymrid ? string("<td>", esc(String(get(t, "mrid", ""))), "</td>") : "",
      "<td>", esc(mode), "</td><td>", get(t, "electrical_step", "?"), "</td><td>", get(t, "fixed_step", "?"), "</td>",
      anypst ? string("<td>", mode == "ratio" ? "-" : get(t, "electrical_shift_step", "?"), "</td><td>", mode == "ratio" ? "-" : get(t, "fixed_shift_step", "?"), "</td>") : "",
      "<td>", esc(status), "</td></tr>",
    ))
  end
  jb = get(metadata, "se_tap_j_before", nothing)
  ja = get(metadata, "se_tap_j_after", nothing)
  db = get(metadata, "se_tap_dof_before", "?")
  da = get(metadata, "se_tap_dof_after", "?")
  jline = jb === nothing || ja === nothing ? "" : "<p>J before fixation: <strong>$(round(Float64(jb); sigdigits = 4))</strong> (dof $(db)) &rarr; J after fixation: <strong>$(round(Float64(ja); sigdigits = 4))</strong> (dof $(da)). After the fixation the tap is no state variable any more; a large J after with a small J before means the estimated position sits between mechanical steps.</p>"
  return string(
    "<section class=\"panel\"><h2>Transformer tap estimates</h2>", jline,
    "<div style=\"overflow-x:auto\"><table class=\"details\"><thead>", head, "</thead><tbody>", join(body, ""), "</tbody></table></div>",
    "<p class=\"field-help\">Full table: se_tap_estimates.csv in the run artifacts.</p></section>",
  )
end

## result-summary card for an SE-started power flow (run_mode == "powerflow_se_start")
function _webui_se_start_summary(result::AbstractDict)::Union{Nothing,String}
  metadata = get(result, "metadata", Dict{String,Any}())
  metadata isa AbstractDict || return nothing
  get(metadata, "run_mode", "") == "powerflow_se_start" || return nothing
  pickup = Float64(get(metadata, "slack_pickup_mw", NaN))
  parts = ["$(get(metadata, "se_start_mode", "?")) from SE run $(get(metadata, "se_run_id", "?")); slack pickup $(round(pickup; digits = 4)) MW"]
  # deviation of the PF solution from the estimated state (how far the
  # model pulls away from the measured operating point)
  maxdvm = Float64(get(metadata, "se_pf_max_dvm_pu", NaN))
  if isfinite(maxdvm)
    push!(parts, "deviation to SE state: max |ΔVm| $(round(maxdvm; sigdigits = 3)) pu at $(get(metadata, "se_pf_max_dvm_bus", "?")), max |Δθ| $(round(Float64(get(metadata, "se_pf_max_dva_deg", NaN)); sigdigits = 3))° at $(get(metadata, "se_pf_max_dva_bus", "?")) (se_pf_deviation.csv)")
  end
  return _webui_escape(join(parts, " | "))
end

## chain action on a successful SE run's result page: dispatch a power flow
## that starts from this estimate (mode select, :se_state default). Plain
## POST to /powerflow/run; deliberately NO onsubmit button-disabling (the
## serialization trap documented in the project learnings).
function _webui_se_chain_section(result::AbstractDict)::String
  metadata = get(result, "metadata", Dict{String,Any}())
  metadata isa AbstractDict || return ""
  get(metadata, "run_mode", "") == "se" || return ""
  get(result, "success", false) || return ""
  run_id = String(get(result, "run_id", ""))
  casefile = String(get(result, "casefile", ""))
  config_file = String(get(result, "config_file", ""))
  isempty(run_id) && return ""
  esc = _webui_escape
  return string(
    "<section class=\"panel\"><h2>Chain: power flow from this estimate</h2>",
    "<p>Starts a power flow from the estimated state. <code>se_state</code> keeps the model injections authoritative (differences go into the slack); <code>se_snapshot</code> also takes the nodal balances from the estimation and converges immediately.</p>",
    "<form method=\"post\" action=\"/powerflow/run\">",
    "<input type=\"hidden\" name=\"casefile\" value=\"$(esc(basename(casefile)))\">",
    "<input type=\"hidden\" name=\"config_file\" value=\"$(esc(config_file))\">",
    "<input type=\"hidden\" name=\"se_start_run_id\" value=\"$(esc(run_id))\">",
    "<label>Start mode <select name=\"se_start_mode\"><option value=\"se_state\" selected>se_state (model authoritative)</option><option value=\"se_snapshot\">se_snapshot (balance takeover)</option></select></label> ",
    "<button type=\"submit\">Run power flow from this estimate</button>",
    "</form>",
    "<p class=\"field-help\">The started run is a normal PF run: the N-1 controls (contingency button, weight editor) work on it unchanged and record this SE run's id in their metadata.</p>",
    "</section>",
  )
end

"""
    render_se_form(; cases, measurements, selected_case, selected_measurement, config_file, message, active_run) -> String

State-estimation page (SE phase 5): case selector (same case sources as the
PowerFlow page), measurement-set selector (content-sniffed CSV v1 files),
the estimator option subset, a demo action generating synthetic measurements
from a solved PF, and the run action posting to `/powerflow/run` with
`se_mode`. Composed from small string pieces like the contingency-weights
editor; no onsubmit button-disabling.
"""
function render_se_form(; cases::Vector{String} = String[], measurements::Vector{String} = String[], case_measurement_rows::Integer = 0, case_measurement_noise::Union{Nothing,Bool} = nothing, selected_case::AbstractString = "", selected_measurement::AbstractString = "", config_file::AbstractString = "", message::AbstractString = "", set_info::Vector{String} = String[], set_provenance::Vector{String} = String[], set_counts::Vector{Tuple{String,Int}} = Tuple{String,Int}[], meas_case::AbstractDict = Dict{String,String}(), starred::AbstractSet = Set{String}(), set_case::AbstractString = "", set_text::AbstractString = "", set_rows::Vector{NamedTuple} = NamedTuple[], gen_values::AbstractDict = Dict{String,String}(), truth_runs::Vector{NamedTuple} = NamedTuple[], active_run = nothing)
  esc = _webui_escape
  # sticky generator inputs: after a generate the redirect carries the
  # submitted values back; the form re-renders with them, not the defaults
  gv(key, default) = esc(String(get(gen_values, key, default)))
  gchk(key) = get(gen_values, key, "") == "true" ? " checked" : ""
  # default-aware variant for fields whose DEFAULT is checked
  gchk2(key, default) = _webui_parse_bool(get(gen_values, key, default)) ? " checked" : ""
  # label sets by their in-file case binding; a set bound to a DIFFERENT
  # case fails the SE run with per-row resolution errors (and the service
  # refuses it up front), so the mismatch must be visible here already
  meas_label = function (name)
    bound = String(get(meas_case, name, ""))
    isempty(bound) && return string(name, " (no case binding)")
    bound == selected_case && return name
    return string(name, " (for ", bound, ")")
  end
  meas_options = join(("<option value=\"$(esc(name))\"$(name == selected_measurement ? " selected" : "")>$(esc(meas_label(name)))</option>" for name in measurements), "")
  # The case file's own set is the FIRST choice when it has one; an empty
  # value means "take what came with the case", which is what the service
  # does. Without one, an empty first entry makes the user pick: preselecting
  # a foreign set armed a run that could only fail the binding check.
  meas_options = if case_measurement_rows > 0
    string("<option value=\"\"$(isempty(selected_measurement) ? " selected" : "")>(from the case file: $(case_measurement_rows) rows)</option>", meas_options)
  elseif isempty(selected_measurement)
    string("<option value=\"\" selected>- pick the set that belongs to this case -</option>", meas_options)
  else
    meas_options
  end
  banner = _webui_active_run_banner(active_run)
  msg_html = isempty(message) ? "" : "<div class=\"alert info\">$(esc(message))</div>"
  # metadata comments of the selected measurement set: the structured tap
  # table (sparlectra-taps v1) renders as a real table, anything else as a
  # plain list; the per-type row counts say WHAT the set contains
  counts_html = isempty(set_counts) ? "" : string("<p><strong>Measured values in this set:</strong> ", join(("$(esc(replace(t, "Meas" => ""))) × $(n)" for (t, n) in set_counts), ", "), " ($(sum(last.(set_counts))) rows; the tap table below is metadata, the measurement rows follow it in the file)</p>")
  info_html = ""
  # binding line + inline editor shared by both info-panel variants
  bind_html = isempty(selected_measurement) ? "" : isempty(set_case) ? "<p><strong>Case binding:</strong> none (older set; the file carries no <code># case:</code> comment). Regenerate or add the line in the editor below.</p>" : "<p><strong>Case binding:</strong> this set belongs to <code>$(esc(set_case))</code> (recorded in the file).</p>"
  # structured value editor: every measurement kind (voltages, flows,
  # balances, currents) editable in place; value/sigma/active only, the
  # identity columns stay read-only. Large sets fall back to text/download.
  table_editor_html = ""
  if !isempty(set_rows)
    trs = String[]
    for r in set_rows
      loc = isempty(r.bus) ? string(r.from_bus, "-", r.to_bus, " #", r.branch_nr, " ", r.direction) : string("bus ", r.bus)
      push!(trs, string(
        "<tr><td><code>", esc(r.id), "</code></td><td>", esc(replace(r.typ, "Meas" => "")), "</td><td>", esc(loc), "</td>",
        "<td><input name=\"v_", r.line, "\" value=\"", esc(r.value), "\" size=\"14\"></td>",
        "<td><input name=\"s_", r.line, "\" value=\"", esc(r.sigma), "\" size=\"10\"></td>",
        "<td><select name=\"a_", r.line, "\"><option value=\"true\"", r.active == "true" ? " selected" : "", ">true</option><option value=\"false\"", r.active == "false" ? " selected" : "", ">false</option></select></td></tr>",
      ))
    end
    capnote = length(set_rows) >= 400 ? "<p class=\"field-help\">Only the first 400 rows are editable here; use the text editor or download for the rest.</p>" : ""
    table_editor_html = string(
      "<details class=\"measurement-editor\"><summary>Edit measurement values (table)</summary>",
      "<form method=\"post\" action=\"/stateestimation/measurements/update-values\" data-busy=\"Updating values from the case…\">",
      "<input type=\"hidden\" name=\"file\" value=\"$(esc(selected_measurement))\">",
      "<input type=\"hidden\" name=\"case\" value=\"$(esc(selected_case))\">",
      "<div style=\"overflow-x:auto;max-height:24rem;overflow-y:auto\"><table class=\"details\"><thead><tr><th>id</th><th>type</th><th>location</th><th>value</th><th>sigma</th><th>active</th></tr></thead><tbody>", join(trs, ""), "</tbody></table></div>",
      capnote,
      "<p><button type=\"submit\">Save values</button> <span class=\"field-help\">Rewrites only value/sigma/active of the changed rows, atomically; a single invalid entry rejects the whole save.</span></p>",
      "</form></details>",
    )
  end
  editor_html = ""
  if !isempty(set_text)
    editor_html = string(
      "<details class=\"measurement-editor\"><summary>Edit measurement file (inline)</summary>",
      "<form method=\"post\" action=\"/stateestimation/measurements/save\">",
      "<input type=\"hidden\" name=\"file\" value=\"$(esc(selected_measurement))\">",
      "<input type=\"hidden\" name=\"case\" value=\"$(esc(selected_case))\">",
      "<textarea name=\"content\" rows=\"18\" spellcheck=\"false\" style=\"width:100%;font-family:monospace;white-space:pre\">$(esc(set_text))</textarea>",
      "<p><button type=\"submit\">Save measurement file</button> <span class=\"field-help\">Saved atomically; the first line must stay <code># sparlectra-measurements v1</code>. For larger files use download and re-upload.</span></p>",
      "</form></details>",
    )
  end
  if !isempty(set_info)
    dl = isempty(selected_measurement) ? "" : "<p><a class=\"button\" href=\"/stateestimation/measurements/download?file=$(_webui_urlencode(selected_measurement))\">Download measurement file</a></p>"
    # generator v2 provenance summary (truth source, flow ends, passive
    # handling), rendered compactly above the taps table
    prov_html = isempty(set_provenance) ? "" : string("<ul class=\"set-provenance\">", join(("<li><code>$(esc(l))</code></li>" for l in set_provenance), ""), "</ul>")
    bind_html = string(bind_html, prov_html)
    if !isempty(set_info) && set_info[1] == "sparlectra-taps v1" && length(set_info) >= 2
      cols = split(set_info[2], ",")
      head = join(("<th>$(esc(String(c)))</th>" for c in cols), "")
      rows = String[]
      for l in set_info[3:end]
        cells = split(l, ","; limit = length(cols))
        push!(rows, string("<tr>", join(("<td>$(esc(String(c)))</td>" for c in cells), ""), "</tr>"))
      end
      info_html = "<section class=\"panel measurement-set-info\"><h2>Measurement set info</h2>$(bind_html)$(counts_html)<h3>Transformer taps at generation</h3><div style=\"overflow-x:auto\"><table class=\"details\"><thead><tr>$(head)</tr></thead><tbody>$(join(rows, ""))</tbody></table></div>$(dl)$(table_editor_html)$(editor_html)</section>"
    else
      info_html = string("<section class=\"panel measurement-set-info\"><h2>Measurement set info</h2>", bind_html, counts_html, "<ul>", join(("<li><code>$(esc(l))</code></li>" for l in set_info), ""), "</ul>$(dl)$(table_editor_html)$(editor_html)</section>")
    end
  elseif !isempty(selected_measurement)
    info_html = "<section class=\"panel measurement-set-info\">$(bind_html)$(counts_html)<p><a class=\"button\" href=\"/stateestimation/measurements/download?file=$(_webui_urlencode(selected_measurement))\">Download measurement file</a></p>$(table_editor_html)$(editor_html)</section>"
  end
  # stage 4A block 4: the SE section shares the page's case selection, so
  # the old per-page case selector is gone; what remains of it is the
  # measurement-set re-upload (same import path as the Case page, returning
  # to this section afterwards)
  upload_form = string(
    "<form method=\"post\" action=\"/powerflow/import-cases\" enctype=\"multipart/form-data\">",
    "<input type=\"hidden\" name=\"return_to\" value=\"stateestimation\">",
    "<input type=\"hidden\" name=\"return_case\" value=\"$(esc(selected_case))\">",
    "<label title=\"Upload a measurement CSV v1 (content-sniffed; lands in the case cache and appears in the set selector)\">Upload measurement file $(_webui_file_input("casefiles"; accept = ".csv", id = "measurements"))</label> ",
    "<button type=\"submit\">Upload</button>",
    "</form>",
  )
  demo_form = isempty(selected_case) ? "" : string(
    "<form method=\"post\" action=\"/stateestimation/generate-measurements\" data-busy=\"Generating, the case is being solved…\">",
    "<input type=\"hidden\" name=\"casefile\" value=\"$(esc(selected_case))\">",
    "<label class=\"check\" title=\"Add seeded Gaussian noise at the sigmas below. Default ON: without noise the measurements match the model exactly and J lands near 0 instead of near dof, which reads like a broken statistic on first sight.\"><input type=\"checkbox\" name=\"noise\" value=\"true\"$(gchk2("noise", "true"))>$(_webui_field_label("se_generator_noise", "noise"))</label> ",
    "<label title=\"Corrupt seed-randomly drawn measurements by k times their sigma. 0 = off; 10 is a good value: clearly detectable bad measurements for the elimination/robust workflow. The count field next to it says how many rows, the seed decides which (documented in the generation message); protected zero-injection rows are never corrupted.\">$(_webui_field_label("se_generator_gross_error", "bad data (k·σ)"))<input type=\"number\" name=\"gross_error_k\" value=\"$(gv("gross_error_k", "0"))\" min=\"0\" max=\"100\" step=\"1\"></label> ",
    "<label title=\"How many measurement rows get the gross error (the seed decides which rows). Only used when bad data k is above 0.\">$(_webui_field_label("se_generator_gross_count", "bad data rows"))<input type=\"number\" name=\"gross_error_count\" value=\"$(gv("gross_error_count", "1"))\" min=\"1\" max=\"1000\" step=\"1\"></label> ",
    "<label title=\"Generate the measurements from a state whose seed-randomly selected in-service transformers run this many MECHANICAL tap steps off the model position (0 = off; whole steps only, a tap changer has no half positions). The deviation lands exactly on each changer's own mechanical grid and is fully recoverable by tap estimation. Creates measurement/model discrepancies around the transformers without touching the model file. Resolve them with the 'estimate taps' option of the estimator run below. Needs truth state 'fresh solve'.\">$(_webui_field_label("se_generator_tap_error", "tap deviation (steps)"))<input type=\"number\" name=\"tap_error_steps\" id=\"gen-tap-steps\" value=\"$(gv("tap_error_steps", "0"))\" min=\"-16\" max=\"16\" step=\"1\"></label> ",
    "<label title=\"At most how many transformers get the tap deviation; the seed decides which. Drawn only from transformers whose deviation the tap estimation can absorb (non-machine, declared changer), so the count is capped there and the message notes when fewer were eligible; machine transformers are used only when the case has nothing else, and the estimator run warns then. Only used when the tap deviation is not 0.\">$(_webui_field_label("se_generator_tap_count", "tap deviation transformers (max)"))<input type=\"number\" name=\"tap_error_count\" id=\"gen-tap-count\" value=\"$(gv("tap_error_count", "1"))\" min=\"1\" max=\"100\" step=\"1\"></label> ",
    "<label title=\"Where the truth state comes from. fresh solve: solves the case now (tol <= 1e-8, island-wise) and generates from that state. from run: adopts the solved voltages of a successful run of THIS case from the run history (PF runs need the detailed result CSV artifact, SE runs use se_state.csv); nothing is re-solved, and the tap deviation is locked.\">$(_webui_field_label("se_generator_truth", "truth state"))<select name=\"gen_truth_source\" id=\"gen-truth-source\"><option value=\"fresh_solve\"$(gv("gen_truth_source", "fresh_solve") == "from_run" ? "" : " selected")>fresh solve</option><option value=\"from_run\"$(gv("gen_truth_source", "fresh_solve") == "from_run" ? " selected" : "")>from run</option></select></label> ",
    "<label title=\"Source run for truth state 'from run': successful power-flow and state-estimation runs of the selected case, newest first\">$(_webui_field_label("se_generator_truth_run", "run"))<select name=\"gen_truth_run_id\" id=\"gen-truth-run\">$(isempty(truth_runs) ? "<option value=\"\">(no successful run of this case)</option>" : join(["<option value=\"$(esc(r.run_id))\"$(r.run_id == String(get(gen_values, "gen_truth_run_id", "")) ? " selected" : "")>$(esc(string(first(r.run_id, 8), "… (", r.kind, ", ", r.timestamp, ")")))</option>" for r in truth_runs], ""))</select></label> ",
    "<label title=\"Flow measurements per branch. both ends: P/Q (and I) at from AND to. one end (balance-aware): exactly one flow group per branch; the end at a bus WITH injection telemetry is preferred, from wins when both or neither have one. The choice is documented per branch in the set comments.\">$(_webui_field_label("se_generator_flow_ends", "flows per branch"))<select name=\"gen_flow_ends\"><option value=\"both\"$(gv("gen_flow_ends", "both") == "one_balance_aware" ? "" : " selected")>both ends</option><option value=\"one_balance_aware\"$(gv("gen_flow_ends", "both") == "one_balance_aware" ? " selected" : "")>one end (balance-aware)</option></select></label> ",
    "<label title=\"Balance sigma (MW/MVar) for passive nodes (no generation, no load, no shunt): their Pinj/Qinj rows become explicit zero balances at this sigma. Small sigmas at many passive nodes stiffen the flat start (a small grid took 32 iterations at 0.01 versus 11 at 0.05); The same value binds the zero-injection constraints when that checkbox is on, with a floor of 0.001 MW (1 kW): a constraint tighter than that does not bind the bus, it makes the normal equations unsolvable (a 25000-bus set with 27268 such rows reached an objective of 3.8e15 per degree of freedom at 1e-6).\">$(_webui_field_label("se_generator_passive_sigma", "passive σ (MW/MVar)"))<input type=\"number\" name=\"gen_passive_sigma\" value=\"$(gv("gen_passive_sigma", "0.05"))\" min=\"0.000001\" step=\"any\"></label> ",
    "<label class=\"check\" title=\"How a passive node (no generation, no load, no shunt) enters the set. OFF: an ordinary balance row Pinj = Qinj = 0 at the sigma on the left, which the estimator weighs against everything else. ON: a protected zero-injection constraint instead (prefix ZI, never eliminated and never down-weighted), written at the sigma on the left but never tighter than 0.001 MW (1 kW). That floor is measured, not cosmetic: at 1e-6 such a row weighs a million times a normal power measurement, and on a 25000-bus set that pushed the normal equations past double precision, so the passive nodes stopped balancing at all (residuals of 182 GW). No duplicate injection rows are written either way.\"><input type=\"checkbox\" name=\"gen_passive_as_zi\" value=\"true\"$(gchk("gen_passive_as_zi"))>$(_webui_field_label("se_generator_passive_zi", "passive as zero-injection constraint"))</label> ",
    "<label title=\"Noise seed: the same seed regenerates the identical set (documented in the set comments), a different seed draws a fresh noise realization\">seed <input type=\"number\" name=\"gen_seed\" value=\"$(gv("gen_seed", "42"))\" min=\"0\" step=\"1\"></label> ",
    "<fieldset class=\"se-options\"><legend>Measurement sigmas $(_webui_help_link("webui.se_generator_sigmas", "Measurement sigmas"))</legend>",
    "<label title=\"Voltage-magnitude accuracy in percent of the measured value (all Vm rows). 0.5 corresponds to a class 0.5 transducer.\">&sigma; U (%) <input type=\"text\" name=\"sigma_u_pct\" value=\"$(gv("sigma_u_pct", "0.5"))\" size=\"6\"></label> ",
    "<label title=\"Also generate branch current-magnitude rows at both ends (currents are auxiliary: gated in the estimator, excluded from observability)\">currents (I) <input type=\"checkbox\" name=\"include_currents\" value=\"true\"$(gchk("include_currents"))></label> ",
    "<label title=\"Current-magnitude accuracy in percent of the measured value (used when currents are enabled)\">&sigma; I (%) <input type=\"text\" name=\"sigma_i_pct\" value=\"$(gv("sigma_i_pct", "1.0"))\" size=\"6\"></label> ",
    "<label title=\"PMU current-phasor ANGLE accuracy in degrees (absolute; an angle passes through zero). 0 = no current-angle rows. Angle rows share the PMU reference offset with voltage angles and are gated near zero current.\">&sigma; Ia (°) <input type=\"text\" name=\"sigma_ia_deg\" value=\"$(gv("sigma_ia_deg", "0"))\" size=\"6\"></label> ",
    "<label title=\"Active-power accuracy in percent of the measured value (injections and branch flows). Percent of reading stays meaningful across voltage levels, unlike an absolute MW sigma.\">&sigma; P (%) <input type=\"text\" name=\"sigma_p_pct\" value=\"$(gv("sigma_p_pct", "1.0"))\" size=\"6\"></label> ",
    "<label title=\"Reactive-power accuracy in percent of the measured value (injections and branch flows)\">&sigma; Q (%) <input type=\"text\" name=\"sigma_q_pct\" value=\"$(gv("sigma_q_pct", "1.0"))\" size=\"6\"></label>",
    "</fieldset>",
    "<button type=\"submit\">Generate synthetic measurements from a solved PF</button>",
    "<p class=\"field-help\">Truth state 'fresh solve' solves the case once (MATPOWER or CGMES); 'from run' adopts a run state unchanged. The CSV v1 lands next to the case and documents truth source, flow-end choices and passive handling in its comments.</p><p class=\"se-recommend\">Recommended: noise on, sigma U 0.5%, P/Q 1%, currents on with sigma I 1%; bad data 10 to exercise the elimination; tap deviation 1 to 2 steps to create a model discrepancy.</p>",
    "<script>(function(){var ts=document.getElementById('gen-truth-source'),tap=document.getElementById('gen-tap-steps'),tapc=document.getElementById('gen-tap-count'),run=document.getElementById('gen-truth-run');function upd(){var fr=ts&&ts.value==='from_run';if(tap){tap.disabled=fr;if(fr){tap.value='0';}}if(tapc){tapc.disabled=fr;}if(run){run.disabled=!fr;}}ts&&ts.addEventListener('change',upd);upd();})();</script>",
    "</form>",
    # per-case settings reset: deletes the saved sidecar profile so both
    # the generator and run forms snap back to the defaults on reload
    "<form method=\"post\" action=\"/stateestimation/reset-settings\">",
    "<input type=\"hidden\" name=\"casefile\" value=\"$(esc(selected_case))\">",
    "<button type=\"submit\" title=\"Delete the saved settings of this case (generator AND run options; the PowerFlow page shares the same per-case profile). All forms fall back to their defaults on the next load.\">Reset saved settings for this case</button>",
    "</form>",
  )
  # a Sparlectra Case Format case carries its measurements with the model, so
  # the run is offered even when the case cache holds no CSV at all
  # a case file that carries its measurements says so where it matters: right
  # above the run action, so nobody generates a set they already have
  # the noise state decides how the resulting J may be read at all, so it is
  # stated here and not left to be inferred from a suspiciously small J
  noise_note = case_measurement_noise === nothing ? " The file does not record whether they carry noise." :
               case_measurement_noise === false ? " They are <strong>ideal (noise-free)</strong> values: an estimation on them returns J = 0 by construction, which is not a quality statement. Use <em>Add noise to this set</em> below to get a realistic run without recomputing a power flow." :
               " They carry measurement noise."
  carries_note = case_measurement_rows > 0 ? "<p class=\"alert info\">This case file carries <strong>$(case_measurement_rows) measurements</strong> with the model. The run uses them as they are; you do not have to generate a set.$(noise_note)</p>" : ""
  # perturbing an existing set needs no power flow, so the action sits with
  # the set it works on and not in the generator fieldset
  noise_form = (case_measurement_rows > 0 || !isempty(measurements)) && !isempty(selected_case) ? string(
    "<form method=\"post\" action=\"/stateestimation/add-noise\" class=\"se-noise-form\" data-busy=\"Adding noise, the case is being read…\">",
    "<input type=\"hidden\" name=\"casefile\" value=\"$(esc(selected_case))\">",
    "<input type=\"hidden\" name=\"measurement_file\" value=\"$(esc(selected_measurement))\">",
    "<label title=\"Seed of the noise draw; the same seed reproduces the same set\">$(_webui_field_label("noise_seed", "noise seed"))<input type=\"number\" name=\"noise_seed\" value=\"42\" min=\"0\"></label> ",
    "<button type=\"submit\" class=\"secondary-button\" title=\"Perturb the values of the set above with each row's own sigma and store the result as a set bound to this case. No power flow is computed and the operating point stays as it is.\">Add noise to this set</button>",
    "</form>",
    "<p class=\"field-help\">Turns an ideal set into a realistic one without recomputing anything: the values are perturbed, the sigmas stay. The result lands next to the case as <code>&lt;case&gt;.noisy.measurements.csv</code> and is preselected afterwards.</p>",
  ) : ""
  run_form = isempty(selected_case) ? "<p>Select a case to run a state estimation.</p>" : (isempty(measurements) && case_measurement_rows == 0) ? "<p>No measurement set found for this case cache. Upload a measurement CSV v1 (Import on the PowerFlow page) or use the demo generator above.</p>" : string(
    carries_note,
    noise_form,
    "<form method=\"post\" action=\"/powerflow/run\">",
    "<input type=\"hidden\" name=\"se_mode\" value=\"true\">",
    "<input type=\"hidden\" name=\"casefile\" value=\"$(esc(selected_case))\">",
    "<input type=\"hidden\" name=\"config_file\" value=\"$(esc(config_file))\">",
    "<p class=\"se-recommend\">Pick the measurement set that belongs to the selected case; the matching set (<case>.measurements.csv) is preselected, sets from other cases are labeled. No set yet? Generate one below.</p><label title=\"Measurement CSV v1 file (content-sniffed)\">$(_webui_field_label("measurement_file", "Measurement set"))<select name=\"measurement_file\">$(meas_options)</select></label>",
    "<fieldset class=\"se-options\"><legend>Estimator options</legend>",
    "<label title=\"Start from a flat voltage profile\">$(_webui_field_label("se_flatstart", "flatstart"))<input type=\"checkbox\" name=\"se_flatstart\" value=\"true\"$(gchk2("se_flatstart", true))></label> ",
    "<label title=\"WLS convergence tolerance on the state step. The default 1e-6 sits AT the noise floor of the finite-difference Jacobian (jac_eps), not below it: a tighter value cannot be reached and the run raises it with a log line. Loosen it for noisy sets, do not tighten it below jac_eps.\">$(_webui_field_label("se_tol", "tol"))<input type=\"text\" name=\"se_tol\" value=\"$(gv("se_tol", "1e-6"))\" size=\"8\"></label> ",
    "<label title=\"Iteration cap, the same limit the service and the configuration use. A converged run usually needs only a few iterations, but a solve with released transformer taps can need close to 40 before it settles, and the count reported afterwards is the one of the LAST solve, not of the most expensive.\">$(_webui_field_label("se_max_iter", "max_iter"))<input type=\"number\" name=\"se_max_iter\" value=\"$(gv("se_max_iter", "50"))\" min=\"1\" max=\"200\"></label> ",
    "<label title=\"What happens to a measurement with a large residual. off: nothing, plain WLS (bad data is only ELIMINATED, see the second limit). staged: down-weighted gradually during the solve, between the two knees k1 and k2. replacement: rows past the down-weight limit get a fixed large sigma and thereby lose their influence. Note on staged: its knees keep the classic |r|/sigma scale, so in this mode the form shows TWO scales side by side, the knees on |r|/sigma and the elimination limit on the normalized residual rn; in replacement mode both limits share the rn scale. All statistics stay on the original sigmas. Rule of thumb: down-weighting for online smoothing, elimination for identification.\">$(_webui_field_label("se_robust_mode", "down-weighting"))<select name=\"se_robust_mode\" id=\"se-robust-mode\"><option value=\"off\"$(gv("se_robust_mode", "off") == "off" ? " selected" : "")>off</option><option value=\"staged\"$(gv("se_robust_mode", "off") == "staged" ? " selected" : "")>staged</option><option value=\"replacement\"$(gv("se_robust_mode", "off") == "replacement" ? " selected" : "")>replacement</option></select></label> ",
    "<label title=\"Second limit: from this NORMALIZED residual rn = r/sqrt(Omega_ii) a measurement is REMOVED from the estimate, not just down-weighted (API keyword normalizedThreshold). Works together with the elimination budget and only while the band test is :high. Same scale as the replacement down-weight limit, so the two are directly comparable.\">$(_webui_field_label("se_k_eliminate", "eliminate from rn"))<input type=\"number\" name=\"se_k_eliminate\" id=\"se-k-eliminate\" value=\"$(gv("se_k_eliminate", "3.0"))\" min=\"0.000001\" step=\"any\"></label> ",
    "<span class=\"se-limit-group\" data-se-limit-group=\"staged\">",
    "<label title=\"First limit, staged mode: below this the weight stays untouched, above it the down-weighting starts and grows. Measured as |r|/sigma, NOT as the normalized residual of the other two limits: the staged knees deliberately keep their classic definition (3/6 reproduces the textbook weighting).\">$(_webui_field_label("se_robust_k1", "down-weight from |r|/sigma"))<input type=\"number\" name=\"se_robust_k1\" value=\"$(gv("se_robust_k1", "3.0"))\" min=\"0.000001\" step=\"any\"></label> ",
    "<label title=\"Staged mode: above this knee the row contributes almost nothing to the gradient any more. Same scale as k1 (|r|/sigma).\">$(_webui_field_label("se_robust_k2", "full down-weight at |r|/sigma"))<input type=\"number\" name=\"se_robust_k2\" value=\"$(gv("se_robust_k2", "6.0"))\" min=\"0.000001\" step=\"any\"></label> ",
    "</span>",
    "<span class=\"se-limit-group\" data-se-limit-group=\"replacement\">",
    "<label title=\"First limit, replacement mode: from this NORMALIZED residual rn = r/sqrt(Omega_ii) a measurement gets the replacement sigma next to it and thereby loses its influence on the estimate; it is NOT removed (that is the elimination limit). Same scale as the elimination limit, so the two are directly comparable.\">$(_webui_field_label("se_k_suppress", "down-weight from rn"))<input type=\"number\" name=\"se_k_suppress\" id=\"se-k-suppress\" value=\"$(gv("se_k_suppress", "4.0"))\" min=\"0.000001\" step=\"any\"></label> ",
    "<label title=\"NOT a limit: the replacement sigma itself, in the unit of the measurement (MW, Mvar, pu). A down-weighted row is solved with this sigma instead of its own, which is why a large value silences it. Statistics and reports keep the original sigma.\">$(_webui_field_label("se_suppression_sigma", "replacement sigma (MW/Mvar)"))<input type=\"number\" name=\"se_suppression_sigma\" value=\"$(gv("se_suppression_sigma", "2000"))\" min=\"0.000001\" step=\"any\"></label> ",
    "</span>",
    "<label title=\"Sequential bad-data elimination budget\">$(_webui_field_label("se_max_eliminations", "max_eliminations"))<input type=\"number\" name=\"se_max_eliminations\" value=\"$(gv("se_max_eliminations", "3"))\" min=\"0\" max=\"20\"></label> ",
    "<span id=\"se-threshold-warning\" class=\"field-help\" style=\"display:none\">Warning: k_suppress is below k_eliminate, so suppressed rows rarely reach the elimination. The run is allowed, but this is usually unintended.</span>",
    # Which limit the selected mode actually uses. Grayed out and disabled,
    # never hidden: hiding makes the form jump and conceals that the setting
    # exists. Same machinery as the solver groups above (classList.toggle
    # plus control.disabled). Two details that would silently lose settings:
    # a disabled control is dropped from the form data per the HTML spec, so
    # the submit handler re-enables everything first (the value has to
    # survive save and reload even while its mode is off), and without
    # JavaScript every field stays usable, as everywhere else in this UI.
    "<script>(function(){var ks=document.getElementById('se-k-suppress'),ke=document.getElementById('se-k-eliminate'),w=document.getElementById('se-threshold-warning');function chk(){if(!ks||!ke||!w)return;var a=parseFloat(ks.value),b=parseFloat(ke.value);w.style.display=(isFinite(a)&&isFinite(b)&&a<b)?'':'none';}ks&&ks.addEventListener('input',chk);ke&&ke.addEventListener('input',chk);chk();" *
    "var mode=document.getElementById('se-robust-mode');var groups=document.querySelectorAll('[data-se-limit-group]');function setInactive(g,inactive){g.classList.toggle('disabled',inactive);g.querySelectorAll('input, select').forEach(function(c){c.disabled=inactive;});}" *
    "function upd(){if(!mode)return;var v=mode.value;groups.forEach(function(g){setInactive(g,g.getAttribute('data-se-limit-group')!==v);});}" *
    "mode&&mode.addEventListener('change',upd);upd();" *
    "var f=mode&&mode.closest('form');f&&f.addEventListener('submit',function(){groups.forEach(function(g){g.querySelectorAll('input, select').forEach(function(c){c.disabled=false;});});});" *
    "})();</script>",
    "<label title=\"Write estimated shunt susceptances back into the model\">$(_webui_field_label("se_update_shunts", "update_shunts"))<input type=\"checkbox\" name=\"se_update_shunts\" value=\"true\"$(gchk("se_update_shunts"))></label> ",
    "<label title=\"Estimate transformer tap positions: releases the tap of every in-service transformer with a ratio tap changer as an extra state, then fixes it to the nearest mechanical step and reruns without the tap state. The result page shows the electrical and fixed steps plus J before/after the fixation. Guarded: machine (generator step-up) transformers are skipped, and a tap the measurement set cannot observe (e.g. a radial transformer without a far-side voltage) is frozen at its current position and reported as frozen instead of absorbing errors.\">$(_webui_field_label("se_tap_estimation", "estimate taps"))<input type=\"checkbox\" name=\"se_tap_estimation\" value=\"true\"$(gchk("se_tap_estimation"))></label> ",
    "<label title=\"Residual-correlation (K matrix) report columns\">$(_webui_field_label("se_report_correlation", "report_correlation"))<input type=\"checkbox\" name=\"se_report_correlation\" value=\"true\"$(gchk("se_report_correlation"))></label>",
    "</fieldset>",
    "<button type=\"submit\">Run state estimation</button>",
    "<p class=\"field-help\">Flow: observability first (traffic light), then the WLS solve, then the bad-data diagnostics and the se_view summary; artifacts land in the run history (kind se).</p>",
    "</form>",
  )
  # stage 4A block 4: this renders the SE SECTION of the Runs page, not a
  # page of its own anymore (GET /stateestimation is a real redirect). The
  # run form stays visible; generator, set info and upload fold into tabs.
  generator_tab = isempty(demo_form) ? "" : "<details class=\"se-tab se-generator\"><summary>Measurement generator</summary>$(demo_form)</details>"
  info_tab = isempty(info_html) ? "" : "<details class=\"se-tab se-set-info\" open><summary>Measurement set details</summary>$(info_html)</details>"
  upload_tab = "<details class=\"se-tab se-upload\"><summary>Upload measurement file</summary>$(upload_form)</details>"
  return "<section class=\"panel se-section\" id=\"state-estimation\"><h2>State estimation</h2>$(banner)$(msg_html)$(run_form)$(generator_tab)$(info_tab)$(upload_tab)</section>"
end

# --- sysimage page -----------------------------------------------------------

"Human-readable size of an existing file, empty string when it is missing."
function _webui_sysimage_size(path::AbstractString)::String
  isfile(path) || return ""
  return string(round(filesize(path) / 1024^2; digits = 1), " MB")
end

"Last `n` non-empty lines of the build log, for the failure case."
function _webui_sysimage_log_tail(path::AbstractString, n::Int = 20)::Vector{String}
  isfile(path) || return String[]
  lines = try
    readlines(path)
  catch err
    # expected failure: the builder holds the file open and may be rotating it
    return ["(build log not readable: $(sprint(showerror, err)))"]
  end
  filter!(l -> !isempty(strip(l)), lines)
  return lines[max(1, end - n + 1):end]
end

"""
    render_webui_sysimage_page(; output_root, message) -> String

The sysimage page: what the image on disk is, whether it is still valid, what
a running build is doing right now, and the button that starts a refresh.

The page auto-refreshes ONLY while a build runs, through the same
`data-refresh-url` mechanism the run status page uses; once the build ends the
rendered page carries no refresh attribute, which makes the client do a real
reload and land on the final state.
"""
function render_webui_sysimage_page(; output_root::AbstractString, message::AbstractString = "")::String
  image = webui_sysimage_path(output_root)
  progress = read_sysimage_build_progress(; output_root)
  active = sysimage_build_active(; output_root)
  build_log = webui_sysimage_build_log_path(output_root)
  problem = try
    webui_sysimage_problem(; image_path = image)
  catch err
    "the validity check failed: $(sprint(showerror, err))"
  end
  flavor = try
    webui_runtime_flavor()
  catch
    (kind = :native, built = nothing)
  end
  flavor_text = flavor.kind === :native ? "native session (no sysimage)" :
                flavor.kind === :app ? "standalone app" :
                flavor.built === nothing ? "sysimage" : "sysimage, built $(replace(first(String(flavor.built), 19), "T" => " "))"

  status_html = if problem === nothing
    "<span class=\"status-badge status-success\">up to date</span>"
  else
    "<span class=\"status-badge status-warning\">needs a rebuild</span> <span class=\"sysimage-reason\">$(_webui_escape(problem))</span>"
  end
  size_text = _webui_sysimage_size(image)

  facts = string(
    "<dl class=\"sysimage-facts\">",
    "<dt>Status</dt><dd>", status_html, "</dd>",
    "<dt>Image file</dt><dd><code>", _webui_escape(image), "</code>", isempty(size_text) ? "" : " (" * _webui_escape(size_text) * ")", "</dd>",
    "<dt>This session</dt><dd><code>", _webui_escape(flavor_text), "</code></dd>",
    "<dt>Build log</dt><dd><code>", _webui_escape(build_log), "</code></dd>",
    "</dl>",
  )

  # A finished rebuild does NOT reach the running process: it kept the image
  # it booted from (see sysimage_restart_pending). Saying that here is the
  # difference between "the refresh worked" and "why is it still slow".
  restart_notice = sysimage_restart_pending(; output_root) ?
                   "<div class=\"alert alert-info\" role=\"status\">A newer image is on disk. This Web UI still runs on the one it started with, so <strong>stop and start it again</strong> to use the new image.</div>" : ""

  # A progress file that exists but does not parse looks EXACTLY like "no
  # build has ever run here", and that is the kind of plausible wrong answer
  # a page must never give: the reader returns nothing for a partial read
  # during a write, which is right for one poll and wrong forever.
  unreadable_progress = progress === nothing && isfile(webui_sysimage_progress_path(output_root)) ?
                        "<div class=\"alert alert-error\" role=\"alert\"><span>The build-progress file exists but cannot be read, so the state of the last build is unknown. Look at the build log, or start a fresh build.</span></div>" : ""

  message_html = isempty(message) ? "" : "<div class=\"alert alert-info\" role=\"status\">$(_webui_escape(message))</div>"

  action_html = if active
    step = get(progress, "step", 0)
    steps = get(progress, "steps", 4)
    phase = String(get(progress, "phase", "working"))
    detail = String(get(progress, "detail", ""))
    elapsed = get(progress, "elapsed_seconds", 0.0)
    minutes = Int(fld(elapsed, 60))
    seconds = Int(floor(elapsed)) % 60
    detail_html = isempty(detail) ? "" : " <span class=\"sysimage-detail\">- $(_webui_escape(detail))</span>"
    string(
      "<section class=\"panel sysimage-progress\">",
      "<h2>Build running</h2>",
      "<p class=\"sysimage-phase\"><strong>[", step, "/", steps, "]</strong> ", _webui_escape(phase), detail_html,
      " <span class=\"sysimage-clock\">", lpad(minutes, 2, '0'), ":", lpad(seconds, 2, '0'), "</span></p>",
      "<progress max=\"", steps, "\" value=\"", step, "\"></progress>",
      "<p class=\"field-help\">The build runs in its own process and survives a browser reload. It keeps working even if you close this page; the previous image stays in use until the new one is finished.</p>",
      "</section>",
    )
  else
    last_html = ""
    if progress !== nothing
      state = String(get(progress, "state", ""))
      note = String(get(progress, "message", ""))
      if state == "failed"
        tail = _webui_sysimage_log_tail(build_log)
        last_html = string(
          "<div class=\"alert alert-error\" role=\"alert\"><span>The last build FAILED: ", _webui_escape(note), "</span></div>",
          isempty(tail) ? "" : "<pre class=\"sysimage-log-tail\">" * _webui_escape(join(tail, "\n")) * "</pre>",
        )
      elseif state == "done"
        last_html = "<p class=\"field-help\">Last build: finished, $(_webui_escape(note)).</p>"
      end
    end
    string(
      "<section class=\"panel sysimage-actions\">",
      "<h2>Refresh the image</h2>",
      last_html,
      "<form method=\"post\" action=\"/webui/sysimage/rebuild\"><button type=\"submit\">Refresh sysimage</button></form>",
      "<p class=\"field-help\">Rebuild the image with the code that is on disk right now. Use this when the image is marked as needing a rebuild, or when a page still paused to compile while you were working: whatever had to be compiled at run time is part of the trace afterwards.</p>",
      "<p class=\"field-help\">It takes a few minutes and runs in the background. The image in use is only replaced at the very end, so this session keeps working while the build runs.</p>",
      "</section>",
    )
  end

  content = string(
    restart_notice,
    unreadable_progress,
    message_html,
    "<section class=\"panel sysimage-status\"><h2>Current image</h2>", facts, "</section>",
    action_html,
  )
  return _webui_layout("Sysimage", content; show_back = true, refresh_url = active ? "/webui/sysimage?autorefresh=1" : nothing)
end

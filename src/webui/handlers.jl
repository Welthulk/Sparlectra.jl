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

# file: src/webui/handlers.jl
# purpose: Web UI response type and action handlers: case import/resolve/
#          delete, run start, config refresh and editor, results, artifacts,
#          and abort/reset endpoints
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
const WEBUI_CASE_IMPORT_MAX_FILE_BYTES = 100 * 1024 * 1024
const WEBUI_CASE_IMPORT_MAX_REQUEST_BYTES = 250 * 1024 * 1024
# a contingency weight list is text (one line per element); the 100 MB case cap
# is wrong here (issue #331 Phase 5 follow-up)
const WEBUI_CONTINGENCY_WEIGHTS_MAX_BYTES = 5 * 1024 * 1024
# cap on RENDERED editor rows: a large case has thousands of elements; showing
# them all as number inputs is unusable, so cap and offer a name filter instead
const WEBUI_CONTINGENCY_WEIGHTS_MAX_ROWS = 200

struct WebUICaseUpload
  filename::String
  data::Vector{UInt8}
end

function handle_webui_logo()::SparlectraWebUIResponse
  isfile(_WEBUI_LOGO_PATH) || return _webui_html(render_webui_error(404, "Sparlectra.jl logo asset is unavailable."); status = 404)
  return SparlectraWebUIResponse(200, ["Content-Type" => "image/png"], read(_WEBUI_LOGO_PATH))
end

function _webui_sanitize_upload_filename(filename::AbstractString)
  raw = String(filename)
  isempty(raw) && return "", "empty filename"
  any(c -> c == '\0' || (iscntrl(c) && c != '\t'), raw) && return basename(raw), "invalid filename"
  base = basename(raw)
  isempty(base) && return "", "empty filename"
  base in (".", "..") && return base, "invalid filename"
  base != raw && return base, "invalid filename"
  occursin('/', base) || occursin('\\', base) ? (base, "invalid filename") : (base, "")
end

function _webui_case_import_uploads(form::AbstractDict)::Vector{WebUICaseUpload}
  value = _webui_form_value(form, "casefiles", WebUICaseUpload[])
  if value isa WebUICaseUpload
    return [value]
  elseif value isa AbstractVector
    return WebUICaseUpload[value...]
  end
  return WebUICaseUpload[]
end

function _webui_write_import_file_atomic(destination::AbstractString, data::Vector{UInt8})
  temp = tempname(dirname(destination))
  try
    open(temp, "w") do io
      write(io, data)
      flush(io)
    end
    mv(temp, destination; force = false)
  catch
    rm(temp; force = true)
    rethrow()
  end
  return nothing
end

"""
    handle_powerflow_case_download(query; output_root, application_root, case_directory)

Send a file from the case directory to the browser. The export buttons write
their result NEXT TO the case (that is where a run can pick it up again), and
this is the way back out: pick the case in the selector, press download, and
the browser saves the file.

Only a plain file inside the case directory is served, addressed by its bare
name: a name with a path separator, a name that escapes the directory, and a
CGMES delivery directory are all refused. Everything else is handed over
verbatim, so a case file, its measurement CSV, or a MATPOWER source all work.
"""
function handle_powerflow_case_download(query::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory = case_directory, application_root = application_root, output_root = output_root)
  requested = strip(String(get(query, "case", "")))
  back(msg) = _webui_redirect("/powerflow/case?import_message=$(_webui_urlencode(msg))")
  isempty(requested) && return back("Select a case before downloading it.")
  # A bare name resolves against the case directory; an absolute path is
  # accepted when it points INTO that directory, because the run form carries
  # the full path back after saving case settings.
  #
  # Both sides are resolved through realpath first. The Web UI state directory
  # is reachable through a symlink (a Flatpak app data path pointing at
  # ~/.local/state), so the form can carry the linked spelling while the
  # runtime holds the resolved one; a purely textual comparison then rejects
  # the user's own case directory.
  resolve(x) = try
    realpath(x)
  catch
    normpath(x)
  end
  path = isabspath(requested) ? resolve(requested) : resolve(joinpath(directory, requested))
  root = string(resolve(directory), Base.Filesystem.path_separator)
  startswith(string(path, Base.Filesystem.path_separator), root) || return back("Invalid case name for the download.")
  isdir(path) && return back("'$(requested)' is a delivery directory, not a single file.")
  isfile(path) || return back("File not found in the case directory: $(requested)")
  # the content type only steers the browser's preview; the disposition is
  # what makes it a download either way
  ext = lowercase(splitext(requested)[2])
  ctype = ext == ".json" ? "application/json; charset=utf-8" :
          ext == ".csv" ? "text/csv; charset=utf-8" :
          ext == ".zip" ? "application/zip" : "text/plain; charset=utf-8"
  return SparlectraWebUIResponse(200, Pair{String,String}["Content-Type" => ctype, "Content-Disposition" => "attachment; filename=\"$(basename(path))\""], read(path))
end

# stage 4A block 4: every SE action lands back on the Runs page's SE
# section (GET /stateestimation is only a redirect, so pointing there would
# just bounce). extra_query carries PRE-ENCODED pairs such as the sticky
# g_* tail of the measurement generator.
function _webui_se_redirect(casefile::AbstractString, message::AbstractString; extra_query::AbstractString = "")
  return _webui_redirect(string("/powerflow?casefile=", _webui_urlencode(casefile), "&message=", _webui_urlencode(message), extra_query, "#state-estimation"))
end

function handle_powerflow_case_import(form::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, operation_log::AbstractString = output_root, max_file_bytes::Integer = WEBUI_CASE_IMPORT_MAX_FILE_BYTES, max_request_bytes::Integer = WEBUI_CASE_IMPORT_MAX_REQUEST_BYTES)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory, application_root, output_root)
  mkpath(directory)
  uploads = _webui_case_import_uploads(form)
  imported = String[]
  imported_roles = Dict{String,String}()
  rejected = Pair{String,String}[]
  if _webui_form_value(form, "case_import_request_oversized", "false") == "true"
    push!(rejected, "request" => "oversized")
  end
  total_bytes = sum((length(upload.data) for upload in uploads); init = 0)
  total_bytes > max_request_bytes && (rejected = [basename(upload.filename) => "oversized" for upload in uploads])
  uploads_to_process = total_bytes > max_request_bytes ? WebUICaseUpload[] : uploads
  for upload in uploads_to_process
    name, reason = _webui_sanitize_upload_filename(upload.filename)
    if !isempty(reason)
      push!(rejected, (isempty(name) ? "(empty)" : name) => reason)
      continue
    elseif !_webui_supported_upload_case_extension(name)
      push!(rejected, name => "unsupported extension")
      continue
    elseif length(upload.data) > max_file_bytes
      push!(rejected, name => "oversized")
      continue
    end
    # a case file is validated BEFORE it is stored: an arbitrary .json must
    # not end up in the case directory just because of its extension
    if lowercase(splitext(name)[2]) == ".json"
      reason = _webui_scf_upload_reason(upload.data)
      if reason !== nothing
        push!(rejected, name => "not a Sparlectra case file ($(reason))")
        continue
      end
    end
    destination = normpath(joinpath(directory, name))
    root = string(normpath(directory), Base.Filesystem.path_separator)
    if destination != normpath(joinpath(directory, basename(destination))) || !startswith(string(destination, Base.Filesystem.path_separator), root)
      push!(rejected, name => "invalid filename")
      continue
    elseif ispath(destination)
      push!(rejected, name => "already exists")
      continue
    end
    try
      _webui_write_import_file_atomic(destination, upload.data)
      push!(imported, name)
      ext = lowercase(splitext(name)[2])
      role = if ext == ".dat"
        _webui_dat_role_label(_webui_classify_dat_content(destination))
      elseif ext == ".zip"
        _webui_cgmes_upload_role(destination, directory)
      elseif ext == ".json"
        "sparlectra_case"
      elseif ext == ".csv"
        # measurement sets (SE phase 5) are detected by the version comment,
        # not the extension; a .csv without it is retained but never offered
        _webui_is_measurement_csv(destination) ? "measurement_set" : "unknown"
      else
        "matpower_case"
      end
      imported_roles[name] = role
    catch
      push!(rejected, name => "write failure")
    end
  end
  selected = ""
  for name in imported
    if _webui_is_user_selectable_case(joinpath(directory, name))
      selected = name
      break
    end
  end
  record_webui_operation!(operation_log, "case_import_completed"; route = "/powerflow/import-cases", method = "POST", user_action = true, selected_count = length(uploads), imported_count = length(imported), rejected_count = length(rejected), imported, rejected = ["$(first(item)): $(last(item))" for item in rejected])
  display_imported = ["$(name) ($(get(imported_roles, name, "unknown")))" for name in imported]
  message = _webui_urlencode(_webui_case_import_message(display_imported, rejected))
  # the SE section hosts its own upload form for measurement sets; return
  # to the Runs page's SE anchor (with its case selection preserved)
  if String(_webui_form_value(form, "return_to", "")) == "stateestimation"
    return _webui_redirect("/powerflow?casefile=$(_webui_urlencode(String(_webui_form_value(form, "return_case", ""))))&message=$(message)#state-estimation")
  end
  query = isempty(selected) ? "" : "?casefile=$(_webui_urlencode(selected))"
  separator = isempty(query) ? "?" : "&"
  return _webui_redirect("/powerflow/case$(query)$(separator)import_message=$(message)")
end

"""
Resolve a manually typed case reference (a bare MATPOWER case name to download,
or a full local path to copy) into the case cache directory, without starting
a PowerFlow run. The resolved file then appears in the "choose existing case"
selector, mirroring `handle_powerflow_case_import` for typed rather than
uploaded cases.
"""
function handle_powerflow_case_resolve(form::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, operation_log::AbstractString = output_root, max_file_bytes::Integer = WEBUI_CASE_IMPORT_MAX_FILE_BYTES)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory, application_root, output_root)
  mkpath(directory)
  # The combined case control posts free-typed values under "casefile";
  # "casefile_manual" is accepted for backward compatibility.
  manual_value = strip(String(something(_webui_form_value(form, "casefile", ""), "")))
  isempty(manual_value) && (manual_value = strip(String(something(_webui_form_value(form, "casefile_manual", ""), ""))))
  if isempty(manual_value)
    message = _webui_urlencode("Enter a case name or path before resolving.")
    return _webui_redirect("/powerflow/case?import_message=$(message)")
  end
  resolved_name = ""
  error_text = ""
  try
    if occursin(r"[\\/]", manual_value)
      # Full/relative path: must already exist locally; copy it into the case
      # cache directory using the same validation as file-upload import.
      isfile(manual_value) || error("case file not found: $(manual_value)")
      name, reason = _webui_sanitize_upload_filename(basename(manual_value))
      isempty(reason) || error("invalid filename ($(reason))")
      _webui_supported_upload_case_extension(name) || error("unsupported extension")
      destination = normpath(joinpath(directory, name))
      root = string(normpath(directory), Base.Filesystem.path_separator)
      (destination == normpath(joinpath(directory, basename(destination))) && startswith(string(destination, Base.Filesystem.path_separator), root)) || error("invalid filename")
      if normpath(manual_value) != destination
        ispath(destination) && error("a case named $(name) already exists in the case directory")
        filesize(manual_value) <= max_file_bytes || error("case file too large")
        _webui_write_import_file_atomic(destination, read(manual_value))
      end
      resolved_name = name
    elseif startswith(lowercase(manual_value), "cgmes:")
      # "cgmes:<alias>": fetch an ENTSO-E CGMES test configuration and pack
      # it (base case plus boundary) as a single ZIP into the case cache.
      alias = strip(manual_value[7:end])
      resolved_path = CGMESImporter.fetchCGMESTestSet(alias; outdir = directory)
      resolved_name = basename(resolved_path)
    else
      # Bare MATPOWER case name: download into the case cache directory.
      resolved_path = ensure_casefile(manual_value; outdir = directory)
      resolved_name = basename(resolved_path)
    end
  catch err
    error_text = sprint(showerror, err)
  end
  if isempty(error_text)
    record_webui_operation!(operation_log, "case_resolve_completed"; route = "/powerflow/resolve-case", method = "POST", user_action = true, requested = manual_value, resolved = resolved_name)
    message = _webui_urlencode("Resolved case: $(resolved_name)")
    return _webui_redirect("/powerflow/case?casefile=$(_webui_urlencode(resolved_name))&import_message=$(message)")
  end
  record_webui_operation!(operation_log, "case_resolve_failed"; route = "/powerflow/resolve-case", method = "POST", user_action = true, requested = manual_value, message = error_text)
  message = _webui_urlencode("Could not resolve case '$(manual_value)': $(error_text)")
  return _webui_redirect("/powerflow/case?import_message=$(message)")
end

"""
Delete a cached case file (and its Web UI case-settings sidecar) from the case
cache directory. Requested from the case combobox via right-click; only bare
filenames inside the case directory are accepted.
"""
function handle_powerflow_case_delete(form::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, operation_log::AbstractString = output_root)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory, application_root, output_root)
  requested = strip(String(something(_webui_form_value(form, "casefile", ""), "")))
  if isempty(requested)
    message = _webui_urlencode("Select a case to delete.")
    return _webui_redirect("/powerflow/case?import_message=$(message)")
  end
  error_text = ""
  if basename(requested) != requested || occursin(r"[\\/]", requested)
    error_text = "invalid case name"
  else
    target = normpath(joinpath(directory, requested))
    if !isfile(target)
      error_text = "case file not found in the case directory"
    else
      try
        # the delete cascade enumerates the companions through the single
        # definition (stage 4B): config, legacy sidecar, weights, and the
        # measurement CSVs, so a deleted case leaves no orphan behind (the
        # measurement sets are new in the cascade; before 4B they stayed)
        companions = case_companion_files(target)
        rm(target)
        for f in companions
          isfile(f) && rm(f; force = true)
        end
      catch err
        error_text = sprint(showerror, err)
      end
    end
  end
  if isempty(error_text)
    record_webui_operation!(operation_log, "case_delete_completed"; route = "/powerflow/delete-case", method = "POST", user_action = true, casefile = requested)
    message = _webui_urlencode("Deleted case: $(requested)")
    return _webui_redirect("/powerflow/case?import_message=$(message)")
  end
  record_webui_operation!(operation_log, "case_delete_failed"; route = "/powerflow/delete-case", method = "POST", user_action = true, casefile = requested, message = error_text)
  message = _webui_urlencode("Could not delete case '$(requested)': $(error_text)")
  return _webui_redirect("/powerflow/case?import_message=$(message)")
end

"""
POST /powerflow/export-scf: write the selected case as a Sparlectra Case
Format file (`<case>.scf.json`, issue #342) into the case directory. The
case is imported through the shared framework import path, so the exported
model is exactly what a run would compute with; the form's own
configuration overrides travel into the file's `sparlectra.config` block,
so the case ships with its working settings. Measurements found next to the
case (`<case>.measurements.csv`) are exported as PGM sensors. Nothing is
solved and nothing on disk is overwritten except the target file.
"""
function handle_powerflow_export_scf(form::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, operation_log::AbstractString = output_root)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory, application_root, output_root)
  requested = strip(String(something(_webui_form_value(form, "casefile", ""), "")))
  back(msg) = _webui_redirect("/powerflow/case?import_message=$(_webui_urlencode(msg))")
  ## on success the redirect also names the written file, so the page can
  ## offer it for download without the user hunting for it in the selector
  done(msg, file) = _webui_redirect("/powerflow/case?import_message=$(_webui_urlencode(msg))&download=$(_webui_urlencode(file))")
  isempty(requested) && return back("Select a case before exporting it as an SCF case file.")
  (basename(requested) == requested && !occursin(r"[\\/]", requested)) || return back("Invalid case name for the SCF export.")
  case_path = normpath(joinpath(directory, requested))
  isfile(case_path) || return back("Case file not found in the case directory: $(requested)")
  # The strict variant writes the plain power-grid-model dataset (no
  # namespaced block), for handing the case to a PGM-only consumer. It gets
  # its own file name so the two never overwrite each other.
  strict = String(_webui_form_value(form, "scf_strict_pgm", "false")) in ("true", "on", "1")
  # the stem drops a format suffix the case already carries, or exporting
  # case14.scf.json as plain PGM would produce case14.scf.pgm.json and every
  # further export would grow another suffix
  stem = first(splitext(basename(requested)))
  (endswith(stem, ".scf") || endswith(stem, ".pgm")) && (stem = first(splitext(stem)))
  out_name = string(stem, strict ? ".pgm.json" : ".scf.json")
  out_path = joinpath(directory, out_name)
  try
    # the request builder validates the same override set the run path uses,
    # so the exported config block can never carry a key a run would reject
    request = powerflow_webui_request(form; default_output_root = output_root)
    overrides = Dict{String,Any}(String(k) => v for (k, v) in get(request, "config_overrides", Dict{String,Any}()))
    config_file = String(get(request, "config_file", DEFAULT_SPARLECTRA_CONFIG_PATH))
    config, _ = _load_api_config(config_file, validate_gui_config_overrides(overrides))
    net = _import_sparlectra_net(case_path, nothing, config)
    sidecar = string(first(splitext(case_path)), ".measurements.csv")
    # the set's own provenance comments travel into the case file: whether the
    # values carry noise decides how a J from them may be read at all
    provenance = Dict{String,Any}()
    if isfile(sidecar)
      readMeasurementsCSV!(net; file = sidecar)
      merge!(provenance, _webui_measurement_set_provenance(sidecar))
    end
    intended = String["power_flow"]
    isempty(net.measurements) || push!(intended, "state_estimation")
    exportSCF(
      net;
      file = out_path,
      case_name = first(splitext(basename(requested))),
      source_format = String(_detect_case_format(case_path)),
      source_reference = basename(requested),
      intended_calculations = intended,
      measurement_provenance = provenance,
      strict_pgm = strict,
    )
    # case-scope settings travel in the case configuration file next to the
    # export, never inside the case file (the writer emits no config block)
    strict || write_case_config(out_path, overrides)
  catch err
    record_webui_operation!(operation_log, "scf_export_failed"; route = "/powerflow/export-scf", method = "POST", user_action = true, casefile = requested, message = sprint(showerror, err))
    return back("SCF export failed for '$(requested)': $(sprint(showerror, err))")
  end
  record_webui_operation!(operation_log, "scf_export_completed"; route = "/powerflow/export-scf", method = "POST", user_action = true, casefile = requested, scf_file = out_name)
  return done("Exported $(out_name) ($(round(filesize(out_path) / 1024; digits = 1)) kB) into the case directory." * (strict ? " Strict PGM: names, roles, tap nameplates, measurements and configuration are NOT in this file." : ""), out_name)
end

"""
Delete a case's saved settings (its case configuration file, plus any
leftover legacy sidecar) without touching the case file. Saved settings
outrank the configuration for their keys, so a stale file can pin a case to
settings the user cannot see or override in the form. This is the way back
to the plain configuration defaults. A deprecated `sparlectra.config` block
inside an SCF case is NOT rewritten here; the form names it when it still
applies.
"""
function handle_powerflow_case_settings_reset(form::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, operation_log::AbstractString = output_root)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory, application_root, output_root)
  requested = strip(String(something(_webui_form_value(form, "casefile", ""), "")))
  if isempty(requested)
    return _webui_redirect("/powerflow/settings?save_message=$(_webui_urlencode("Select a case before resetting its saved settings."))")
  end
  if basename(requested) != requested || occursin(r"[\\/]", requested)
    record_webui_operation!(operation_log, "case_settings_reset_failed"; route = "/powerflow/case-settings/reset", method = "POST", user_action = true, casefile = requested, message = "invalid case name")
    return _webui_redirect("/powerflow/settings?save_message=$(_webui_urlencode("Could not reset settings for '$(requested)': invalid case name"))")
  end
  case_config = _webui_case_settings_path(output_root, requested; case_directory = directory)
  legacy = _webui_legacy_case_settings_path(output_root, requested; case_directory = directory)
  removed = String[]
  try
    if isfile(case_config)
      rm(case_config; force = true)
      push!(removed, "case configuration file")
    end
    if isfile(legacy)
      rm(legacy; force = true)
      push!(removed, "legacy sidecar")
    end
  catch err
    record_webui_operation!(operation_log, "case_settings_reset_failed"; route = "/powerflow/case-settings/reset", method = "POST", user_action = true, casefile = requested, message = sprint(showerror, err))
    return _webui_redirect("/powerflow/settings?casefile=$(_webui_urlencode(requested))&save_message=$(_webui_urlencode("Could not reset settings for '$(requested)': $(sprint(showerror, err))"))")
  end
  if isempty(removed)
    record_webui_operation!(operation_log, "case_settings_reset_noop"; route = "/powerflow/case-settings/reset", method = "POST", user_action = true, casefile = requested, status = "no_saved_settings")
    return _webui_redirect("/powerflow/settings?casefile=$(_webui_urlencode(requested))&save_message=$(_webui_urlencode("No saved settings for '$(requested)': the form already uses the configuration defaults."))")
  end
  record_webui_operation!(operation_log, "case_settings_reset_completed"; route = "/powerflow/case-settings/reset", method = "POST", user_action = true, casefile = requested, cleared = join(removed, ","))
  return _webui_redirect("/powerflow/settings?casefile=$(_webui_urlencode(requested))&save_message=$(_webui_urlencode("Saved settings for '$(requested)' deleted ($(join(removed, ", "))): the form now uses the configuration defaults."))")
end

# --- Contingency weights editor (issue #331 Phase 5 follow-up) ---

# redirect back to the weights editor for a case, carrying a status message
function _webui_weights_redirect(casefile::AbstractString, message::AbstractString)::SparlectraWebUIResponse
  prefix = isempty(casefile) ? "?" : "?case=$(_webui_urlencode(casefile))&"
  return _webui_redirect(string("/powerflow/contingency-weights", prefix, "weights_message=", _webui_urlencode(message)))
end

"""
    handle_contingency_weights_upload(form; output_root, application_root, case_directory, operation_log)

Store (or replace) the per-case N-1 contingency weight file
`<stem>.contingency-weights.csv` next to the case. The target case comes from
the form's `casefile` key; the uploaded file arrives under `casefiles` (the
multipart parser collects every file part there regardless of the input's name,
see `_webui_parse_multipart_form`). The upload is validated with
`readContingencyWeightsCSV` from a temp copy BEFORE anything is stored, so a
malformed CSV is rejected with the parser's line-numbered message and any
existing file is left untouched. Unlike case import, uploading REPLACES an
existing weight file (a weight list is a working document); the replacement is
stated in the redirect message.
"""
function handle_contingency_weights_upload(form::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, operation_log::AbstractString = output_root)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory = case_directory, application_root = application_root, output_root = output_root)
  requested = strip(String(something(_webui_form_value(form, "casefile", ""), "")))
  route = "/powerflow/contingency-weights/upload"
  reject = function (msg::AbstractString)
    record_webui_operation!(operation_log, "contingency_weights_upload_failed"; route = route, method = "POST", user_action = true, casefile = requested, message = msg)
    return _webui_weights_redirect(requested, "Weights upload rejected: $(msg)")
  end
  isempty(requested) && return reject("no case selected")
  (basename(requested) != requested || occursin(r"[\\/]", requested)) && return reject("invalid case name")
  isfile(joinpath(directory, requested)) || return reject("case \"$(requested)\" not found")
  # every uploaded file part arrives under "casefiles" regardless of the input's
  # name attribute (multipart parser); the weights file is here too
  uploads = _webui_case_import_uploads(form)
  isempty(uploads) && return reject("no file uploaded")
  upload = first(uploads)
  lowercase(splitext(basename(upload.filename))[2]) == ".csv" || return reject("weights file must be a .csv")
  length(upload.data) > WEBUI_CONTINGENCY_WEIGHTS_MAX_BYTES && return reject("file exceeds $(div(WEBUI_CONTINGENCY_WEIGHTS_MAX_BYTES, 1024 * 1024)) MB")
  # validate BEFORE storing: parse a temp copy; readContingencyWeightsCSV names
  # the offending line, and nothing is written on failure
  try
    tmp = tempname()
    write(tmp, upload.data)
    try
      readContingencyWeightsCSV(tmp)
    finally
      rm(tmp; force = true)
    end
  catch err
    return reject(first(split(sprint(showerror, err), '\n')))
  end
  dest = try
    _webui_case_weights_path(requested; case_directory = directory)
  catch err
    return reject(first(split(sprint(showerror, err), '\n')))
  end
  startswith(normpath(dest), normpath(directory)) || return reject("invalid destination")
  replaced = isfile(dest)
  # atomic replace (unlike case import, which rejects an existing file)
  temp = tempname(dirname(dest))
  try
    write(temp, upload.data)
    mv(temp, dest; force = true)
  catch err
    rm(temp; force = true)
    return reject(first(split(sprint(showerror, err), '\n')))
  end
  record_webui_operation!(operation_log, "contingency_weights_upload_completed"; route = route, method = "POST", user_action = true, casefile = requested, replaced = replaced)
  return _webui_weights_redirect(requested, replaced ? "Weights replaced for '$(requested)'." : "Weights uploaded for '$(requested)'.")
end

"""
    handle_contingency_weights_page(query; ...) -> SparlectraWebUIResponse

Render the N-1 contingency weights editor for `query["case"]`. Seeds the table
with the case's real element names (`generateN1Branches` + `generateN1Generators`
on the net built through the shared config-driven import), so the user never
types names blind; a large case is capped to `WEBUI_CONTINGENCY_WEIGHTS_MAX_ROWS`
rendered rows with a name filter. A case whose format cannot be listed here (e.g.
CGMES/DTF through this MATPOWER import) falls back to the raw-CSV editor.
"""
function handle_contingency_weights_page(query::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, config_file::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log::AbstractString = output_root)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory = case_directory, application_root = application_root, output_root = output_root)
  cases = _webui_casefile_options_in_directory(directory)
  case = strip(String(get(query, "case", "")))
  message = String(get(query, "weights_message", ""))
  filter_text = strip(String(get(query, "filter", "")))
  _webui_log_route!(operation_log, "contingency_weights_opened", "GET", "/powerflow/contingency-weights"; case = case)
  if isempty(case) || basename(case) != case || occursin(r"[\\/]", case)
    note = isempty(case) ? message : "Invalid case name."
    return _webui_html(render_contingency_weights_editor(; case = "", cases = cases, elements = String[], stored = Dict{String,Float64}(), raw_text = "", message = note))
  end
  wf = _webui_case_weights_path(case; case_directory = directory)
  raw_text = isfile(wf) ? read(wf, String) : ""
  stored = Dict{String,Float64}()
  isfile(wf) && try
    stored = readContingencyWeightsCSV(wf)
  catch
  end
  # seed element names from the case (expensive: build the net once)
  net_error = ""
  all_elements = String[]
  case_path = joinpath(directory, case)
  if isfile(case_path)
    try
      config = load_sparlectra_config(config_file; reload = true)
      net = _import_sparlectra_net(case_path, nothing, config)
      all_elements = vcat([c.name for c in generateN1Branches(net)], [c.name for c in generateN1Generators(net)])
    catch err
      net_error = first(split(sprint(showerror, err), '\n'))
    end
  else
    net_error = "case file not found"
  end
  total = length(all_elements)
  shown = isempty(filter_text) ? all_elements : [e for e in all_elements if occursin(lowercase(filter_text), lowercase(e))]
  capped = first(shown, WEBUI_CONTINGENCY_WEIGHTS_MAX_ROWS)
  # stage 4A block 4: fragment=1 serves the editor body alone for the Runs
  # page's lazy-loading weights tab (same renderer as the full page)
  if get(query, "fragment", "") in ("1", "true")
    return _webui_html(_webui_weights_editor_fragment(; case = case, elements = capped, stored = stored, raw_text = raw_text, message = message, filter = filter_text, total_count = total, net_error = net_error))
  end
  return _webui_html(render_contingency_weights_editor(; case = case, cases = cases, elements = capped, stored = stored, raw_text = raw_text, message = message, filter = filter_text, total_count = total, net_error = net_error))
end

"""
    handle_contingency_weights_save(form; ...) -> SparlectraWebUIResponse

Write the per-case weight file from the editor. The raw-CSV textarea wins when
submitted (empty means clear); otherwise the file is built from the seeded
table's `element`/`weight` pairs, OMITTING rows left at exactly 1.0 so the file
stays a diff of the default. The result is validated with
`readContingencyWeightsCSV` before it replaces the stored file; a table that
produces no non-default row deletes the file.
"""
function handle_contingency_weights_save(form::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, operation_log::AbstractString = output_root)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory = case_directory, application_root = application_root, output_root = output_root)
  requested = strip(String(something(_webui_form_value(form, "casefile", ""), "")))
  route = "/powerflow/contingency-weights/save"
  reject = function (msg::AbstractString)
    record_webui_operation!(operation_log, "contingency_weights_save_failed"; route = route, method = "POST", user_action = true, casefile = requested, message = msg)
    return _webui_weights_redirect(requested, "Save rejected: $(msg)")
  end
  isempty(requested) && return reject("no case selected")
  (basename(requested) != requested || occursin(r"[\\/]", requested)) && return reject("invalid case name")
  dest = try
    _webui_case_weights_path(requested; case_directory = directory)
  catch err
    return reject(first(split(sprint(showerror, err), '\n')))
  end
  startswith(normpath(dest), normpath(directory)) || return reject("invalid destination")
  raw = _webui_form_value(form, "weights_text", nothing)
  text = if raw !== nothing
    String(raw)                               # textarea form; empty means clear
  else
    els = _webui_form_value(form, "element", String[])
    wts = _webui_form_value(form, "weight", String[])
    els = els isa AbstractVector ? String.(els) : (isempty(strip(String(els))) ? String[] : [String(els)])
    wts = wts isa AbstractVector ? String.(wts) : (isempty(strip(String(wts))) ? String[] : [String(wts)])
    lines = String["name;weight"]
    for (e, w) in zip(els, wts)
      v = tryparse(Float64, strip(w))
      v === nothing && return reject("weight for '$(e)' is not a number")
      v == 1.0 && continue                     # omit default rows: a diff, not a dump
      push!(lines, string(e, ";", v))
    end
    length(lines) == 1 ? "" : string(join(lines, "\n"), "\n")
  end
  if isempty(strip(text))
    isfile(dest) && rm(dest; force = true)
    record_webui_operation!(operation_log, "contingency_weights_cleared"; route = route, method = "POST", user_action = true, casefile = requested)
    return _webui_weights_redirect(requested, "Weights cleared for '$(requested)' (no non-default rows).")
  end
  try
    tmp = tempname()
    write(tmp, text)
    try
      readContingencyWeightsCSV(tmp)
    finally
      rm(tmp; force = true)
    end
  catch err
    return reject(first(split(sprint(showerror, err), '\n')))
  end
  temp = tempname(dirname(dest))
  try
    write(temp, text)
    mv(temp, dest; force = true)
  catch err
    rm(temp; force = true)
    return reject(first(split(sprint(showerror, err), '\n')))
  end
  record_webui_operation!(operation_log, "contingency_weights_saved"; route = route, method = "POST", user_action = true, casefile = requested)
  return _webui_weights_redirect(requested, "Weights saved for '$(requested)'.")
end

"""Serve the stored per-case weight file as a CSV download."""
function handle_contingency_weights_download(query::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory = case_directory, application_root = application_root, output_root = output_root)
  case = strip(String(get(query, "case", "")))
  (isempty(case) || basename(case) != case || occursin(r"[\\/]", case)) && return _webui_weights_redirect(case, "invalid case name")
  wf = try
    _webui_case_weights_path(case; case_directory = directory)
  catch
    return _webui_weights_redirect(case, "invalid case name")
  end
  isfile(wf) || return _webui_weights_redirect(case, "no weight file for '$(case)'")
  return SparlectraWebUIResponse(200, Pair{String,String}["Content-Type" => "text/csv; charset=utf-8", "Content-Disposition" => "attachment; filename=\"$(basename(wf))\""], read(wf))
end

"""Delete the per-case weight file (the reset action of the weights editor)."""
function handle_contingency_weights_reset(form::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, operation_log::AbstractString = output_root)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory = case_directory, application_root = application_root, output_root = output_root)
  requested = strip(String(something(_webui_form_value(form, "casefile", ""), "")))
  route = "/powerflow/contingency-weights/reset"
  if isempty(requested) || basename(requested) != requested || occursin(r"[\\/]", requested)
    record_webui_operation!(operation_log, "contingency_weights_reset_failed"; route = route, method = "POST", user_action = true, casefile = requested, message = "invalid case name")
    return _webui_weights_redirect(requested, "invalid case name")
  end
  wf = _webui_case_weights_path(requested; case_directory = directory)
  if !isfile(wf)
    record_webui_operation!(operation_log, "contingency_weights_reset_noop"; route = route, method = "POST", user_action = true, casefile = requested, status = "no_file")
    return _webui_weights_redirect(requested, "No weight file to delete for '$(requested)'.")
  end
  rm(wf; force = true)
  record_webui_operation!(operation_log, "contingency_weights_reset_completed"; route = route, method = "POST", user_action = true, casefile = requested)
  return _webui_weights_redirect(requested, "Weight file deleted for '$(requested)'.")
end

## --- scenario editor (scenario task step 6) ----------------------------------

_webui_scenarios_redirect(case::AbstractString, message::AbstractString) = _webui_redirect(string("/powerflow/scenarios?case=", _webui_urlencode(case), "&scenario_message=", _webui_urlencode(message)))

# valid editor target: a bare .scf.json file name in the case directory
_webui_scenarios_case_ok(case::AbstractString)::Bool = !isempty(case) && basename(case) == case && !occursin(r"[\\/]", case) && endswith(lowercase(case), ".scf.json")

# component options of a typed case for the editor: (class, id, label) with
# "id: name (from-to)" labels resolved through the extra names
function _webui_scenario_component_options(case::SCFCase)
  names = Dict{Int,String}()
  spar = case.sparlectra
  if spar !== nothing
    for (k, v) in spar.extra
      v isa AbstractDict || continue
      haskey(v, "name") || continue
      id = tryparse(Int, String(k))
      id === nothing && continue
      names[id] = String(v["name"])
    end
  end
  node_name(id) = get(names, id, string(id))
  lab(id) = string(id, ": ", get(names, id, "?"))
  lab2(id, a, b) = string(id, ": ", get(names, id, "?"), " (", node_name(a), "-", node_name(b), ")")
  opts = Vector{Tuple{Symbol,Int,String}}()
  d = case.data
  for r in d.line
    push!(opts, (:branch, r.id, lab2(r.id, r.from_node, r.to_node)))
  end
  for r in d.generic_branch
    push!(opts, (:transformer, r.id, lab2(r.id, r.from_node, r.to_node)))
  end
  for r in d.link
    push!(opts, (:link, r.id, lab2(r.id, r.from_node, r.to_node)))
  end
  for r in d.sym_gen
    push!(opts, (:generator, r.id, lab(r.id)))
  end
  for r in d.sym_load
    push!(opts, (:load, r.id, lab(r.id)))
  end
  for r in d.shunt
    push!(opts, (:shunt, r.id, lab(r.id)))
  end
  for r in d.source
    push!(opts, (:external_grid, r.id, lab(r.id)))
  end
  return opts
end

# one scenario from the editor form's parallel op-row fields; rows without a
# component are skipped (a cloned empty row must not fail the save)
function _webui_scenario_from_form(form::AbstractDict)
  aslist(key) = begin
    v = _webui_form_value(form, key, String[])
    v isa AbstractVector ? String.(v) : String[String(v)]
  end
  name = strip(String(something(_webui_form_value(form, "scenario_name", ""), "")))
  wtext = strip(String(something(_webui_form_value(form, "scenario_weight", "1.0"), "1.0")))
  weight = something(tryparse(Float64, wtext), 1.0)
  ops_op = aslist("op_op")
  ops_target = aslist("op_target")
  ops_component = aslist("op_component")
  ops_field = aslist("op_field")
  ops_value = aslist("op_value")
  ops_factor = aslist("op_factor")
  n = maximum(length.((ops_op, ops_target, ops_component, ops_field, ops_value, ops_factor)))
  at(v, i) = i <= length(v) ? strip(v[i]) : ""
  ops = PatchOp[]
  for i in 1:n
    comp = at(ops_component, i)
    isempty(comp) && continue
    id = tryparse(Int, comp)
    id === nothing && throw(ArgumentError(string("op ", i, ": component id ", repr(comp), " is not an integer")))
    fld = at(ops_field, i)
    val = at(ops_value, i)
    fac = at(ops_factor, i)
    push!(ops, PatchOp(
      op = Symbol(at(ops_op, i)),
      target = Symbol(at(ops_target, i)),
      id = id,
      field = isempty(fld) ? nothing : Symbol(fld),
      value = isempty(val) ? nothing : tryparse(Float64, val),
      factor = isempty(fac) ? nothing : tryparse(Float64, fac),
    ))
  end
  return Scenario(name = name, weight = weight, ops = ops)
end

"""
    _webui_scenarios_editor_tab(case, directory) -> String

The scenario editor as an embeddable Runs-page fragment (stage 4A block
4): scenario list plus an empty new-scenario form for the selected SCF
case. Editing or duplicating an existing scenario navigates to the
standalone editor page, which shares the same fragment renderer. Returns
"" when the case is not a readable SCF file; the Runs page then simply
shows no scenario tab (the standalone route still reports the reason).
"""
function _webui_scenarios_editor_tab(case::AbstractString, directory::AbstractString)::String
  _webui_scenarios_case_ok(case) || return ""
  case_path = joinpath(directory, case)
  isfile(case_path) || return ""
  set, opts = try
    c = read_scf_json(case_path)
    sset = scf_case_scenarios(c)
    (sset === nothing ? ScenarioSet() : sset, _webui_scenario_component_options(c))
  catch
    return ""
  end
  return _webui_scenarios_editor_fragment(; case = case, scenarios = set.scenarios, component_options = opts, file_digest = _config_file_hash(case_path))
end

"""
    handle_scenarios_page(query; ...) -> SparlectraWebUIResponse

Render the scenario editor for `query["case"]` (SCF cases only, maintainer
revision 2026-09-03): the case's scenario list plus the scenario form.
`scenario` selects a scenario into the form (`mode=duplicate` copies it
under a new name); without it the form starts a new scenario.
"""
function handle_scenarios_page(query::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, config_file::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log::AbstractString = output_root)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory = case_directory, application_root = application_root, output_root = output_root)
  scf_cases = [c for c in _webui_casefile_options_in_directory(directory) if endswith(lowercase(c), ".scf.json")]
  case = strip(String(get(query, "case", "")))
  message = String(get(query, "scenario_message", ""))
  _webui_log_route!(operation_log, "scenarios_opened", "GET", "/powerflow/scenarios"; case = case)
  if !_webui_scenarios_case_ok(case)
    note = isempty(case) ? message : "Scenarios need an SCF case (.scf.json); export this case as SCF first."
    return _webui_html(render_scenarios_editor(; case = "", cases = scf_cases, message = note))
  end
  case_path = joinpath(directory, case)
  isfile(case_path) || return _webui_html(render_scenarios_editor(; case = "", cases = scf_cases, message = "case file not found: $(case)"))
  scfcase, set, opts = try
    c = read_scf_json(case_path)
    s = scf_case_scenarios(c)
    (c, s === nothing ? ScenarioSet() : s, _webui_scenario_component_options(c))
  catch err
    return _webui_html(render_scenarios_editor(; case = "", cases = scf_cases, message = string("cannot read ", case, ": ", first(split(sprint(showerror, err), '\n')))))
  end
  sel = strip(String(get(query, "scenario", "")))
  mode = strip(String(get(query, "mode", "")))
  form_scenario = nothing
  original_name = ""
  if !isempty(sel)
    i = findfirst(s -> s.name == sel, set.scenarios)
    if i !== nothing
      s = set.scenarios[i]
      if mode == "duplicate"
        form_scenario = Scenario(name = string(s.name, " (copy)"), weight = s.weight, ops = s.ops)
      else
        form_scenario = s
        original_name = s.name
      end
    end
  end
  return _webui_html(render_scenarios_editor(; case = case, cases = scf_cases, scenarios = set.scenarios, component_options = opts, form_scenario = form_scenario, original_name = original_name, message = message, file_digest = _config_file_hash(case_path)))
end

"""
    handle_scenarios_save(form; ...) -> SparlectraWebUIResponse

Validate the submitted scenario against the case (the step-1 structural
rules plus the net-aware step-2 rules, so a tap patch on a regulated
transformer is rejected with the controller named) and write the scenarios
block into the SCF file through the existing writer. On a validation error
the editor re-renders WITH the submitted form content and the error next
to the form; on success it redirects with a message.
"""
function handle_scenarios_save(form::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, config_file::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log::AbstractString = output_root)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory = case_directory, application_root = application_root, output_root = output_root)
  scf_cases = [c for c in _webui_casefile_options_in_directory(directory) if endswith(lowercase(c), ".scf.json")]
  case = strip(String(something(_webui_form_value(form, "casefile", ""), "")))
  original_name = strip(String(something(_webui_form_value(form, "original_name", ""), "")))
  route = "/powerflow/scenarios/save"
  if !_webui_scenarios_case_ok(case)
    record_webui_operation!(operation_log, "scenarios_save_failed"; route = route, method = "POST", user_action = true, casefile = case, message = "invalid case name")
    return _webui_scenarios_redirect(case, "invalid case name")
  end
  case_path = joinpath(directory, case)
  isfile(case_path) || return _webui_scenarios_redirect(case, "case file not found")
  scenario, form_error = try
    (_webui_scenario_from_form(form), "")
  catch err
    (nothing, sprint(showerror, err))
  end
  scfcase = read_scf_json(case_path)
  opts = _webui_scenario_component_options(scfcase)
  set = something(scf_case_scenarios(scfcase), ScenarioSet())
  current_digest = _config_file_hash(case_path)
  if scenario === nothing
    return _webui_html(render_scenarios_editor(; case = case, cases = scf_cases, scenarios = set.scenarios, component_options = opts, original_name = original_name, error_text = form_error, file_digest = current_digest))
  end
  # lost-update guard: two tabs saving into the same SCF file must not
  # silently overwrite each other; the form carries the digest of the file
  # it was rendered from, and a changed file rejects the save with the way
  # out (the re-rendered form keeps the input and the CURRENT digest, so a
  # reviewed resubmit goes through)
  submitted_digest = strip(String(something(_webui_form_value(form, "file_digest", ""), "")))
  if !isempty(submitted_digest) && submitted_digest != current_digest
    stale_msg = "the case file changed since this form was loaded (another tab or process saved); review the list below and save again"
    record_webui_operation!(operation_log, "scenarios_save_stale"; route = route, method = "POST", user_action = true, casefile = case)
    return _webui_html(render_scenarios_editor(; case = case, cases = scf_cases, scenarios = set.scenarios, component_options = opts, form_scenario = scenario, original_name = original_name, error_text = stale_msg, file_digest = current_digest))
  end
  if scfcase.sparlectra === nothing
    return _webui_scenarios_redirect(case, "this file carries no sparlectra block (a plain PGM dataset); export the case as SCF first")
  end
  scenarios = Scenario[s for s in set.scenarios if s.name != original_name]
  push!(scenarios, scenario)
  newset = ScenarioSet(scenarios = scenarios, mode = set.mode, exclusions = set.exclusions)
  # server-side validation: the step-1 structural rules, then the
  # net-aware step-2 rules (regulated-tap rejection names the controller)
  err_text = try
    index = ScenarioIndex(scfcase)
    validate_scenarios(newset, index)
    config = load_sparlectra_config(config_file; reload = true)
    net = build_net(scfcase; config = config)
    validate_scenarios(newset, index, net)
    ""
  catch err
    first(split(sprint(showerror, err), '\n'))
  end
  if !isempty(err_text)
    record_webui_operation!(operation_log, "scenarios_save_rejected"; route = route, method = "POST", user_action = true, casefile = case, message = err_text)
    return _webui_html(render_scenarios_editor(; case = case, cases = scf_cases, scenarios = set.scenarios, component_options = opts, form_scenario = scenario, original_name = original_name, error_text = err_text, file_digest = current_digest))
  end
  scfcase.sparlectra.scenarios = scenario_set_dict(newset)
  write_scf_json(scfcase, case_path)
  record_webui_operation!(operation_log, "scenarios_saved"; route = route, method = "POST", user_action = true, casefile = case, scenario = scenario.name, ops = length(scenario.ops))
  return _webui_scenarios_redirect(case, string("Scenario '", scenario.name, "' saved into ", case, " (", length(scenario.ops), " op(s))."))
end

"""Delete one scenario from the case file's scenarios block."""
function handle_scenarios_delete(form::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, operation_log::AbstractString = output_root)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory = case_directory, application_root = application_root, output_root = output_root)
  case = strip(String(something(_webui_form_value(form, "casefile", ""), "")))
  name = strip(String(something(_webui_form_value(form, "scenario", ""), "")))
  route = "/powerflow/scenarios/delete"
  _webui_scenarios_case_ok(case) || return _webui_scenarios_redirect(case, "invalid case name")
  case_path = joinpath(directory, case)
  isfile(case_path) || return _webui_scenarios_redirect(case, "case file not found")
  # lost-update guard, like the save handler
  submitted_digest = strip(String(something(_webui_form_value(form, "file_digest", ""), "")))
  if !isempty(submitted_digest) && submitted_digest != _config_file_hash(case_path)
    record_webui_operation!(operation_log, "scenarios_delete_stale"; route = route, method = "POST", user_action = true, casefile = case)
    return _webui_scenarios_redirect(case, "the case file changed since the page was loaded (another tab or process saved); review the list and delete again")
  end
  scfcase = read_scf_json(case_path)
  scfcase.sparlectra === nothing && return _webui_scenarios_redirect(case, "this file carries no sparlectra block")
  set = something(scf_case_scenarios(scfcase), ScenarioSet())
  remaining = Scenario[s for s in set.scenarios if s.name != name]
  if length(remaining) == length(set.scenarios)
    return _webui_scenarios_redirect(case, string("No scenario named '", name, "' in ", case, "."))
  end
  scfcase.sparlectra.scenarios = isempty(remaining) ? Dict{String,Any}() : scenario_set_dict(ScenarioSet(scenarios = remaining, mode = set.mode, exclusions = set.exclusions))
  write_scf_json(scfcase, case_path)
  record_webui_operation!(operation_log, "scenarios_deleted"; route = route, method = "POST", user_action = true, casefile = case, scenario = name)
  return _webui_scenarios_redirect(case, string("Scenario '", name, "' deleted from ", case, "."))
end

"""
    _webui_merge_case_config_write(source, keep_updates, form_updates; on_unreadable) -> NamedTuple

Merge updates into the case configuration file next to `source` and write
it: the existing file is the base (`keep_updates` win on conflict), existing
form-block fields survive a save from a page that does not carry them (the
PF result save must not wipe the SE generator options; `case_format` is the
one non-spec form-block field, stage 4A, and survives the same way), and
`form_updates` win over stored form values. An unreadable existing file is
replaced deliberately (the save is the user's explicit request to persist
the current state); `on_unreadable` is called with the error so the caller
can log the loss. Returns `(path, config, form)` of what was written.
"""
function _webui_merge_case_config_write(source::AbstractString, keep_updates::AbstractDict, form_updates::AbstractDict; on_unreadable::Function = err -> nothing)
  existing_config = try
    load_case_config(source)
  catch err
    on_unreadable(err)
    Dict{String,Any}()
  end
  config = merge(existing_config, keep_updates)
  form_fields = Dict{String,Any}()
  for (field, value) in _webui_case_form_defaults(source, nothing)
    (haskey(_WEBUI_OPTION_BY_FIELD, String(field)) || String(field) == "case_format") && (form_fields[String(field)] = value)
  end
  merge!(form_fields, form_updates)
  path = write_case_config(source, config; form = form_fields)
  return (; path, config, form = form_fields)
end

"""
Save the Case page's import options (input format, CGMES import options,
MATPOWER import conventions) into the case configuration file next to the
case (stage 4A). Config-key fields become case-scope configuration entries,
`case_format` goes into the form block; runs pick both up through the
configuration precedence instead of the run POST.
"""
function handle_case_options_save(form::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, operation_log::AbstractString = output_root)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory, application_root, output_root)
  case = strip(String(something(_webui_form_value(form, "casefile", ""), "")))
  route = "/powerflow/case/options/save"
  back = casefile -> "/powerflow/case" * (isempty(casefile) ? "" : "?casefile=" * _webui_urlencode(casefile))
  redirect_message = (casefile, message) -> _webui_redirect(string(back(casefile), isempty(casefile) ? "?" : "&", "import_message=", _webui_urlencode(message)))
  isempty(case) && return redirect_message("", "Select a case before saving case options.")
  source = isabspath(case) ? normpath(case) : _webui_resolve_case_profile_source(case; case_directory = directory)
  if isempty(source) || !(isfile(source) || isdir(source))
    record_webui_operation!(operation_log, "case_options_save_failed"; route, method = "POST", user_action = true, casefile = case, status = "rejected", message = "case not found")
    return redirect_message(case, "Case '$(case)' was not found in the case directory; resolve it first.")
  end
  keep = Dict{String,Any}()
  dropped = String[]
  try
    for (config_key, field, type) in _WEBUI_FORM_CONFIG_FIELDS
      raw = _webui_form_value(form, field, nothing)
      raw === nothing && continue
      # same case-scope split as the result-page save: machine-scope keys
      # (output, benchmark, runtime, webui, export) never enter the case file
      if scf_is_case_config_key(String(config_key))
        keep[String(config_key)] = _webui_parse_form_value(raw, type, field)
      else
        push!(dropped, String(config_key))
      end
    end
  catch err
    record_webui_operation!(operation_log, "case_options_save_failed"; route, method = "POST", user_action = true, casefile = case, status = "rejected", message = sprint(showerror, err))
    return redirect_message(case, "Could not save case options: $(sprint(showerror, err))")
  end
  form_updates = Dict{String,Any}()
  fmt = lowercase(strip(String(something(_webui_form_value(form, "case_format", ""), ""))))
  # scf and pgm name the same reader; both are accepted so the choice
  # made in the form survives the save (see _normalize_case_format)
  fmt in ("auto", "matpower", "dtf_for001", "cgmes", "scf", "pgm") && (form_updates["case_format"] = fmt)
  written = try
    _webui_merge_case_config_write(source, keep, form_updates; on_unreadable = err -> record_webui_operation!(operation_log, "case_options_save_replaced_unreadable"; route, method = "POST", user_action = true, casefile = case, status = "replaced", message = sprint(showerror, err)))
  catch err
    record_webui_operation!(operation_log, "case_options_save_failed"; route, method = "POST", user_action = true, casefile = case, status = "rejected", message = sprint(showerror, err))
    return redirect_message(case, "Could not write the case configuration file: $(sprint(showerror, err))")
  end
  record_webui_operation!(operation_log, "case_options_saved"; route, method = "POST", user_action = true, casefile = case, status = "succeeded", profile_path = written.path, saved_keys = length(written.config), form_fields = length(written.form), machine_scope_dropped = join(sort!(dropped), ","))
  return redirect_message(case, "Saved case options to $(basename(written.path)).")
end

"""
Save the Settings page's options (stage 4A block 3). Target `this_case`
writes case-scope configuration keys into the selected case's
configuration file (machine-scope keys are named and kept out, exactly
like the result-page save) and the request-option fields
(performance_timing, detailed CSV, export_cgmes) into its form block.
Target `general` merges ALL configuration keys into the general YAML
(one `.settings-save.bak` backup beside it); the per-case form-block
fields have no general home and are named as kept out there. Runs pick
everything up through `resolve_config`; the run POST carries none of
these fields any more.
"""
function handle_settings_save(form::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, operation_log::AbstractString = output_root)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory, application_root, output_root)
  case = strip(String(something(_webui_form_value(form, "casefile", ""), "")))
  config_file = strip(String(something(_webui_form_value(form, "config_file", ""), "")))
  isempty(config_file) && (config_file = DEFAULT_SPARLECTRA_CONFIG_PATH)
  target = strip(String(something(_webui_form_value(form, "settings_target", "this_case"), "this_case")))
  route = "/powerflow/settings/save"
  back = msg -> _webui_redirect(string("/powerflow/settings", isempty(case) ? "?" : string("?casefile=", _webui_urlencode(case), "&"), "save_message=", _webui_urlencode(msg)))
  # collect the POSTed option fields once: config keys typed via their
  # specs, request options separately
  config_updates = Dict{String,Any}()
  try
    for (config_key, field, type) in _WEBUI_FORM_CONFIG_FIELDS
      raw = _webui_form_value(form, field, nothing)
      raw === nothing && continue
      config_updates[String(config_key)] = _webui_parse_form_value(raw, type, field)
    end
  catch err
    record_webui_operation!(operation_log, "settings_save_failed"; route, method = "POST", user_action = true, casefile = case, target, status = "rejected", message = sprint(showerror, err))
    return back("Could not save settings: $(sprint(showerror, err))")
  end
  form_updates = Dict{String,Any}()
  for field in _WEBUI_CASE_PROFILE_EXTRA_FIELDS
    spec = _webui_option_spec(field)
    spec.config_key === nothing || continue
    raw = _webui_form_value(form, field, nothing)
    raw === nothing && continue
    form_updates[field] = _webui_parse_form_value(raw, spec.value_type, field)
  end
  if target == "this_case"
    isempty(case) && return back("Select a case first (Case page), or save to the configuration file.")
    source = isabspath(case) ? normpath(case) : _webui_resolve_case_profile_source(case; case_directory = directory)
    (isempty(source) || !(isfile(source) || isdir(source))) && return back("Case '$(case)' was not found in the case directory; resolve it first.")
    keep = Dict{String,Any}()
    dropped = String[]
    for (key, value) in config_updates
      scf_is_case_config_key(key) ? (keep[key] = value) : push!(dropped, key)
    end
    written = try
      _webui_merge_case_config_write(source, keep, form_updates; on_unreadable = err -> record_webui_operation!(operation_log, "settings_save_replaced_unreadable"; route, method = "POST", user_action = true, casefile = case, status = "replaced", message = sprint(showerror, err)))
    catch err
      record_webui_operation!(operation_log, "settings_save_failed"; route, method = "POST", user_action = true, casefile = case, target, status = "rejected", message = sprint(showerror, err))
      return back("Could not write the case configuration file: $(sprint(showerror, err))")
    end
    record_webui_operation!(operation_log, "settings_saved"; route, method = "POST", user_action = true, casefile = case, target, status = "succeeded", profile_path = written.path, saved_keys = length(written.config), form_fields = length(written.form), machine_scope_dropped = join(sort!(dropped), ","))
    note = isempty(dropped) ? "" : " Machine-scope keys kept out (save to the configuration file instead): $(join(sort!(dropped), ", "))."
    return back("Saved settings for this case to $(basename(written.path)).$(note)")
  end
  target == "general" || return back("Unknown settings target '$(target)'.")
  isfile(config_file) || return back("Configuration file not found: $(config_file).")
  nested = try
    validate_gui_config_overrides(config_updates)
  catch err
    record_webui_operation!(operation_log, "settings_save_failed"; route, method = "POST", user_action = true, target, status = "rejected", message = sprint(showerror, err))
    return back("Could not save settings: $(sprint(showerror, err))")
  end
  backup = string(config_file, ".settings-save.bak")
  try
    cp(config_file, backup; force = true)
    merged = _merge_config_overrides(load_yaml_dict(config_file), nested)
    _write_yaml_file(config_file, merged)
  catch err
    record_webui_operation!(operation_log, "settings_save_failed"; route, method = "POST", user_action = true, target, status = "rejected", message = sprint(showerror, err))
    return back("Could not write the configuration file: $(sprint(showerror, err))")
  end
  record_webui_operation!(operation_log, "settings_saved"; route, method = "POST", user_action = true, target, status = "succeeded", config_file, saved_keys = length(config_updates), backup_path = backup)
  note = isempty(form_updates) ? "" : " Per-case options kept out (save with target this case): $(join(sort!(collect(keys(form_updates))), ", "))."
  return back("Saved $(length(config_updates)) key(s) to $(basename(config_file)) (backup: $(basename(backup))).$(note)")
end

"""Run a PowerFlow request through the Web UI form-to-service boundary."""
function handle_powerflow_run(form::AbstractDict; default_output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, runner = start_powerflow_run, operation_log::AbstractString = default_output_root)::Dict{String,Any}
  request = powerflow_webui_request(form; default_output_root = default_output_root, case_directory = case_directory)
  package_case_directory = joinpath(application_root, "data", "mpower")
  requested_case = String(request["casefile"])
  if !isabspath(requested_case) && !occursin('/', requested_case) && !occursin('\\', requested_case)
    cached_case = case_directory === nothing ? "" : joinpath(case_directory, requested_case)
    if isfile(cached_case)
      request["casefile"] = cached_case
    else
      # bundled sources (data/mpower, and data/scf with the shipped demo
      # cases) are staged into the cache with their sidecars on first use
      staged = _webui_stage_bundled_case!(application_root, case_directory, requested_case)
      staged === nothing || (request["casefile"] = staged)
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

function handle_powerflow_config_editor(config_file::AbstractString; message::AbstractString = "", status::Integer = 200, config_text::Union{Nothing,AbstractString} = nothing)::SparlectraWebUIResponse
  path = isempty(strip(config_file)) ? DEFAULT_SPARLECTRA_CONFIG_PATH : String(config_file)
  text = config_text === nothing ? (isfile(path) ? read(path, String) : "") : String(config_text)
  notice = isempty(message) ? "" : "<div class=\"alert info\">$(_webui_escape(message))</div>"
  sidecar_note = "<p class=\"alert warning\">Case-sidecar settings, when present for a selected case, have higher precedence than this global YAML. Delete or refresh the case settings from a run result if they should no longer override the edited file.</p>"
  body = """
<section class="panel">
<h2>Configuration Editor</h2>
$(notice)
<p>Active configuration file: <code>$(_webui_escape(path))</code></p>
$(sidecar_note)
<form method="post" action="/powerflow/config/edit">
<input type="hidden" name="config_file" value="$(_webui_escape(path))">
<textarea name="config_text" rows="30" style="width: 100%; font-family: monospace;">$(_webui_escape(text))</textarea>
<p><button type="submit">Save YAML configuration</button></p>
</form>
</section>
"""
  return _webui_html(_webui_layout("Configuration Editor", body; show_back = true); status)
end

function _validate_powerflow_config_editor_text!(config_text::AbstractString)::Dict{String,Any}
  mktemp() do path, io
    write(io, config_text)
    close(io)
    duplicates = _detect_yaml_duplicate_keys(path)
    if !isempty(duplicates)
      throw(ArgumentError("Duplicate YAML key(s) detected: $(join(duplicates, ", ")). Save refused."))
    end
    load_yaml_dict(path)
    _load_and_validate_config(DEFAULT_SPARLECTRA_CONFIG_PATH, path; cli_overrides = Dict{String,Any}(), overrides = Dict{String,Any}())
    load_sparlectra_config(path; reload = true)
    return load_yaml_dict(path)
  end
end

function _config_editor_error_message(err)::String
  return string("Configuration could not be saved.\n\n", sprint(showerror, err))
end

function _config_editor_backup_path(config_file::AbstractString)::String
  stamp = Dates.format(now(), "yyyymmdd-HHMMSS")
  candidate = string(config_file, ".bak-", stamp)
  !isfile(candidate) && return candidate
  for i in 1:10_000
    numbered = string(candidate, "-", i)
    !isfile(numbered) && return numbered
  end
  throw(ArgumentError("Could not choose a unique backup path for $(config_file)."))
end

function handle_powerflow_config_editor_save(form::AbstractDict; operation_log::AbstractString = "results/powerflow_service")::SparlectraWebUIResponse
  config_file = strip(String(something(_webui_form_value(form, "config_file", ""), "")))
  config_text = String(something(_webui_form_value(form, "config_text", ""), ""))
  tmp_path = nothing
  try
    isempty(config_file) && return _webui_html(render_webui_error(400, "No configuration file was provided."); status = 400)
    parsed = _validate_powerflow_config_editor_text!(config_text)
    parent = dirname(abspath(config_file))
    mkpath(parent)
    tmp_path = tempname(parent; cleanup = false)
    write(tmp_path, _yaml_dict_text(parsed))
    _load_and_validate_config(DEFAULT_SPARLECTRA_CONFIG_PATH, tmp_path; cli_overrides = Dict{String,Any}(), overrides = Dict{String,Any}())
    load_sparlectra_config(tmp_path; reload = true)
    backup_path = _config_editor_backup_path(config_file)
    isfile(config_file) && cp(config_file, backup_path; force = false)
    mv(tmp_path, config_file; force = true)
    tmp_path = nothing
    load_sparlectra_config(config_file; reload = true)
    record_webui_operation!(operation_log, "config_editor_saved"; route = "/powerflow/config/edit", method = "POST", user_action = true, config_file, backup_path)
    return handle_powerflow_config_editor(config_file; message = "Configuration saved. Backup: $(backup_path)")
  catch err
    record_webui_operation!(operation_log, "config_editor_save_failed"; route = "/powerflow/config/edit", method = "POST", user_action = true, config_file, message = sprint(showerror, err))
    return handle_powerflow_config_editor(config_file; message = _config_editor_error_message(err), status = 400, config_text)
  finally
    tmp_path !== nothing && isfile(tmp_path) && rm(tmp_path; force = true)
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
  value isa Union{Bool,Integer,AbstractFloat,Symbol,AbstractString} || throw(ArgumentError("Case-settings field $(field) has unsupported value type $(typeof(value))."))
  return _webui_normalize_case_profile_form_value(field, value)
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
  settings = Dict{String,Any}(spec.field => spec.default for spec in WEBUI_OPTION_SPECS if spec.save_in_case_sidecar)
  for (raw_key, raw_value) in settings_raw
    key = String(raw_key)
    field = get(_WEBUI_CASE_PROFILE_CONFIG_FIELD_BY_KEY, key, key)
    field in _WEBUI_CASE_PROFILE_FIELDS || continue
    _webui_case_profile_setting!(settings, field, raw_value)
  end
  return settings
end

function _webui_case_settings_saved_html(path::AbstractString, casefile::AbstractString, count::Integer, successful::Bool, override::Bool; label::AbstractString = "Saved case configuration:", extra::AbstractString = "")::String
  status = successful ? "successful/converged" : (override ? "non-successful, saved via explicit override" : "non-successful")
  body = """
<section class=\"panel case-settings-saved\"><h2>Case settings saved</h2>
<p class=\"alert info\"><strong>$(_webui_escape(label))</strong> <code>$(_webui_escape(path))</code></p>
<ul>$(extra)
<li><strong>Saved settings:</strong> $(count)</li>
<li><strong>Run status:</strong> $(_webui_escape(status))</li>
</ul>
<p>These settings will be applied when the same case is selected again. Manual edits in the form still override the saved profile for each run.</p>
<p><a href=\"/powerflow?casefile=$(_webui_urlencode(casefile))\">Open the run form with this case</a></p>
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
  # Saved settings live in the case configuration file next to the case
  # (`<stem>.config.yaml`, D8): the case-scope configuration keys as
  # sections, the request-only fields (SE and generator options) in the
  # `form` block. Machine-scope keys (output, benchmark, runtime, webui,
  # export) stay with the machine and are named in the log. No sidecar and
  # no in-case-file block are written any more.
  keep, dropped_keys = _webui_case_config_from_settings(settings_raw)
  # The run's request settings carry only the keys its form submitted, so
  # writing them alone REPLACED the case configuration and dropped every
  # case-scope key without a form control (live case: the converted
  # sidecar's matpower_import.apply_bus_names vanished on save and the
  # measurement set stopped resolving bus names). The existing file is the
  # base, the run's keys win on conflict. An unreadable existing file is
  # replaced deliberately: the save is the user's explicit request to
  # persist the current state, and the operation log names the loss.
  form_updates = Dict{String,Any}()
  for (field, value) in settings
    spec = get(_WEBUI_OPTION_BY_FIELD, String(field), nothing)
    (spec === nothing || spec.config_key !== nothing || !spec.save_in_case_sidecar) && continue
    # `settings` is normalized WITH spec defaults; only fields the request
    # actually carried may overwrite a stored value, or every save resets
    # the other pages' options (the generator seed) to their defaults
    haskey(settings_raw, String(field)) || continue
    form_updates[String(field)] = value
  end
  written = try
    _webui_merge_case_config_write(case_settings_source, keep, form_updates; on_unreadable = err -> record_webui_operation!(operation_log, "case_settings_save_replaced_unreadable"; route = "/powerflow/result/$(run_id)/case-settings/save", method = "POST", user_action = true, run_id, status = "replaced", message = sprint(showerror, err)))
  catch err
    record_webui_operation!(operation_log, "case_settings_save_failed"; route = "/powerflow/result/$(run_id)/case-settings/save", method = "POST", user_action = true, run_id, status = "rejected", message = sprint(showerror, err))
    return _webui_html(render_webui_error(400, "Could not write the case configuration file: $(sprint(showerror, err))"); status = 400)
  end
  written_path = written.path
  keep = written.config
  form_fields = written.form
  # a leftover legacy sidecar would be converted over the fresh save on the
  # next load; remove it now that the case configuration file is the record
  stale_sidecar = false
  try
    legacy = _webui_legacy_case_settings_path(output_root, case_settings_source)
    stale_sidecar = isfile(legacy)
    stale_sidecar && rm(legacy; force = true)
  catch
  end
  record_webui_operation!(operation_log, "case_settings_saved"; route = "/powerflow/result/$(run_id)/case-settings/save", method = "POST", user_action = true, run_id, status = "succeeded", profile_path = written_path, saved_keys = length(keep), form_fields = length(form_fields), machine_scope_dropped = join(dropped_keys, ","), stale_sidecar_removed = stale_sidecar)
  case_dir = _webui_case_directory(; application_root = _webui_application_root(), output_root = output_root)
  back_name = dirname(abspath(case_settings_source)) == abspath(case_dir) ? basename(case_settings_source) : case_settings_source
  return _webui_html(_webui_case_settings_saved_html(written_path, back_name, length(keep) + length(form_fields), successful, !successful && override; label = "Saved case configuration:", extra = string(
    isempty(dropped_keys) ? "" : "<li><strong>Kept on this machine (not case scope):</strong> $(length(dropped_keys)) key(s)</li>",
    stale_sidecar ? "<li><strong>Old settings sidecar removed:</strong> the case configuration file is the record now</li>" : "",
  )))
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

function _zip_crc32(bytes::Vector{UInt8})::UInt32
  crc = 0xffffffff % UInt32
  for byte in bytes
    crc ⊻= UInt32(byte)
    for _ in 1:8
      crc = (crc & 0x00000001) == 1 ? (crc >> 1) ⊻ 0xedb88320 % UInt32 : crc >> 1
    end
  end
  return ~crc
end

function _zip_write_u16(io::IO, value::Integer)
  write(io, UInt8(value & 0xff), UInt8((value >> 8) & 0xff))
end

function _zip_write_u32(io::IO, value::Integer)
  for shift in (0, 8, 16, 24)
    write(io, UInt8((value >> shift) & 0xff))
  end
end

function _build_uncompressed_zip(entries::Vector{Pair{String,Vector{UInt8}}})::Vector{UInt8}
  io = IOBuffer()
  central = IOBuffer()
  for (name, bytes) in entries
    name_bytes = Vector{UInt8}(codeunits(name))
    crc = _zip_crc32(bytes)
    offset = position(io)
    _zip_write_u32(io, 0x04034b50)
    _zip_write_u16(io, 20); _zip_write_u16(io, 0); _zip_write_u16(io, 0)
    _zip_write_u16(io, 0); _zip_write_u16(io, 0)
    _zip_write_u32(io, crc); _zip_write_u32(io, length(bytes)); _zip_write_u32(io, length(bytes))
    _zip_write_u16(io, length(name_bytes)); _zip_write_u16(io, 0)
    write(io, name_bytes); write(io, bytes)
    _zip_write_u32(central, 0x02014b50)
    _zip_write_u16(central, 20); _zip_write_u16(central, 20); _zip_write_u16(central, 0); _zip_write_u16(central, 0)
    _zip_write_u16(central, 0); _zip_write_u16(central, 0)
    _zip_write_u32(central, crc); _zip_write_u32(central, length(bytes)); _zip_write_u32(central, length(bytes))
    _zip_write_u16(central, length(name_bytes)); _zip_write_u16(central, 0); _zip_write_u16(central, 0)
    _zip_write_u16(central, 0); _zip_write_u16(central, 0); _zip_write_u32(central, 0); _zip_write_u32(central, offset)
    write(central, name_bytes)
  end
  central_bytes = take!(central)
  central_offset = position(io)
  write(io, central_bytes)
  _zip_write_u32(io, 0x06054b50)
  _zip_write_u16(io, 0); _zip_write_u16(io, 0); _zip_write_u16(io, length(entries)); _zip_write_u16(io, length(entries))
  _zip_write_u32(io, length(central_bytes)); _zip_write_u32(io, central_offset); _zip_write_u16(io, 0)
  return take!(io)
end

function handle_powerflow_artifacts_zip(run_id::AbstractString)::SparlectraWebUIResponse
  artifacts = list_powerflow_artifacts(run_id)
  artifacts isa AbstractDict && return _webui_html(render_webui_error(404, get(artifacts, "message", "Run not found.")); status = 404)
  entries = Pair{String,Vector{UInt8}}[]
  for artifact in artifacts
    name = String(artifact["name"])
    occursin("..", name) && continue
    (occursin('/', name) || occursin('\\', name)) && continue
    path = String(artifact["path"])
    isfile(path) || continue
    push!(entries, name => read(path))
  end
  isempty(entries) && return _webui_html(render_webui_error(404, "No artifacts are available for this run."); status = 404)
  filename = "sparlectra_run_$(replace(String(run_id), r"[^A-Za-z0-9_.-]" => "_"))_artifacts.zip"
  headers = Pair{String,String}["Content-Type" => "application/zip", "Content-Disposition" => "attachment; filename=\"$(filename)\""]
  return SparlectraWebUIResponse(200, headers, _build_uncompressed_zip(entries))
end

function handle_powerflow_history(output_root::AbstractString)::SparlectraWebUIResponse
  return _webui_html(render_powerflow_history(list_powerflow_runs(output_root), output_root; active_run = get_active_webui_powerflow_job()))
end

"""
    handle_powerflow_compare(run_ids) -> SparlectraWebUIResponse

Two finished runs side by side. Exactly two ids are required, and both must be
loadable; anything else is answered with a message instead of a half-filled
page, because a comparison against a missing run is worse than none.
"""
function handle_powerflow_compare(run_ids::Vector{String})::SparlectraWebUIResponse
  ids = [String(strip(id)) for id in run_ids if !isempty(strip(id))]
  if length(ids) != 2
    return _webui_html(render_webui_error(400,
      "Comparing needs exactly two runs; $(length(ids)) were selected. Go back to the run history, tick two rows and press Compare."); status = 400)
  end
  a = get_webui_powerflow_job(ids[1])
  b = get_webui_powerflow_job(ids[2])
  for (id, entry) in zip(ids, (a, b))
    get(entry, "reason", "") == "run_not_found" && return _webui_html(render_webui_error(404, "Run $(id) was not found."); status = 404)
    # the page reads power-flow artifacts; another kind has nothing it would
    # show, and a half-empty comparison reads like a finding
    metadata = get(entry, "metadata", Dict{String,Any}())
    kind = metadata isa AbstractDict ? string(get(metadata, "run_mode", "")) : ""
    _webui_comparable_kind(kind) || return _webui_html(render_webui_error(400,
      "Run $(id) is a $(kind) run. The comparison reads power-flow results (configuration, Q-limit events, losses, bus voltages), so both runs have to be power-flow runs."); status = 400)
  end
  return _webui_html(render_powerflow_compare(a, b))
end

function handle_webui_operation_log(output_root::AbstractString; download::Bool = false)::SparlectraWebUIResponse
  path = webui_operation_log_path(output_root)
  content = isfile(path) ? read(path, String) : ""
  if download
    headers = ["Content-Disposition" => "attachment; filename=\"$(WEBUI_OPERATION_LOG_FILENAME)\""]
    return SparlectraWebUIResponse(200, content; content_type = "application/x-ndjson; charset=utf-8", headers)
  end
  return _webui_html(render_webui_operation_log(content; entries = count(==(0x0a), codeunits(content)), bytes = ncodeunits(content)))
end

"""
    handle_webui_operation_log_clear(output_root) -> SparlectraWebUIResponse

Empty the operation log from the Web UI. The file stays in place and gets
ONE entry recording that it was cleared, so the log never becomes silently
empty and the following entries have a visible starting point.

Retention drops old entries at every start; this is the manual way to a
fresh log without restarting, for a session that produced more than anyone
wants to read.
"""
function handle_webui_operation_log_clear(output_root::AbstractString)::SparlectraWebUIResponse
  path = webui_operation_log_path(output_root)
  removed = isfile(path) ? count(==(0x0a), codeunits(read(path, String))) : 0
  try
    isfile(path) && open(path, "w") do io
    end
  catch err
    return _webui_html(render_webui_error(500, "Could not clear the operation log: $(sprint(showerror, err))"); status = 500)
  end
  record_webui_operation!(output_root, "operation_log_cleared"; route = "/webui/operation-log/clear", method = "POST", user_action = true, status = "cleared", removed_entries = removed)
  return _webui_redirect("/webui/operation-log")
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

# ---------------------------------------------------------------------------
# State estimation page handlers (SE phase 5)
# ---------------------------------------------------------------------------

"""
    _webui_se_form_state(query; output_root, application_root, case_directory, config_file, selected_fallback) -> NamedTuple

Assembles the state of the state-estimation section embedded on the Runs
page (stage 4A block 4): case and measurement-set discovery, set info and
editors, sticky generator values, and truth-run candidates. `query` carries
the SE-specific keys (`case`, `message`, `g_*`); when `case` is absent, the
page's shared case selection arrives as `selected_fallback`. Returns the
keyword set of `render_se_form`. GET /stateestimation itself is a real
redirect to the Runs page (no second rendering path).
"""
function _webui_se_form_state(query::AbstractDict; output_root::AbstractString, application_root::AbstractString = _webui_application_root(), case_directory = nothing, config_file::AbstractString = "", selected_fallback::AbstractString = "")
  directory = _webui_case_directory(; case_directory = case_directory, application_root = application_root, output_root = output_root)
  cases = _webui_casefile_options_in_directory(directory)
  selected = String(get(query, "case", ""))
  isempty(selected) && (selected = basename(String(selected_fallback)))
  isempty(selected) || selected in cases || (selected = "")
  measurements = _webui_measurement_options_in_directory(directory)
  message = String(get(query, "message", ""))
  # per-set case binding: the in-file `# case:` comment decides; the stem
  # convention <case>.measurements.csv only fills in for older files
  meas_case = Dict{String,String}()
  for name in measurements
    bound = _webui_measurement_set_case(joinpath(directory, name))
    if isempty(bound) && endswith(name, ".measurements.csv")
      stem = String(name[1:(end - length(".measurements.csv"))])
      hits = [c for c in cases if first(splitext(c)) == stem]
      length(hits) == 1 && (bound = hits[1])
    end
    isempty(bound) || (meas_case[name] = bound)
  end
  # a case with at least one bound set gets the star in the selector
  starred = Set{String}(values(meas_case))
  # preselect a set BOUND to the selected case; a foreign set fails the SE
  # run with per-row resolution errors, so it must never be a silent
  # default. WITHOUT a selected case nothing is preselected and no set
  # panel renders (a fresh page must not show set info out of thin air).
  bound_sets = [name for name in measurements if get(meas_case, name, "") == selected]
  # What the case file says about its OWN measurement set: how many rows it
  # carries, and whether those are ideal (noise-free) or measured. Parsed
  # once, because both answers come from the same block.
  case_measurement_rows, case_measurement_noise = try
    if !isempty(selected) && lowercase(splitext(selected)[2]) == ".json"
      block = get(get(scf_json_parse(read(joinpath(directory, selected), String)), "sparlectra", Dict{String,Any}()), "measurements", nothing)
      if block isa AbstractDict
        # one entry of `rows` is one PGM SENSOR, and a sensor carries one or
        # two measurement rows (a P/Q pair). The page must show measurements,
        # not sensors, or it understates the set.
        n = sum(length(get(r, "sparlectra_ids", [])) for r in get(block, "rows", []); init = 0)
        prov = get(block, "provenance", nothing)
        noise = prov isa AbstractDict ? get(prov, "noise", nothing) : nothing
        (n, noise isa Bool ? noise : nothing)
      else
        (0, nothing)
      end
    else
      (0, nothing)
    end
  catch
    (0, nothing)
  end
  # Precedence: a set BOUND to this case wins, because generating one is a
  # deliberate act and the newer data. Among several bound sets the
  # CONVENTION file <case>.measurements.csv wins, exactly as the page
  # promises (before load_fixture_net this was only alphabetical luck:
  # case14... sorted before case9..., warmup_casePST... does not). The
  # measurements carried inside a case file are the default only while no
  # bound set exists; without either, the page asks instead of arming a
  # foreign set.
  convention_set = isempty(selected) ? "" : string(first(splitext(selected)), ".measurements.csv")
  selected_meas = isempty(selected) ? "" : (convention_set in bound_sets ? convention_set : (!isempty(bound_sets) ? first(bound_sets) : ""))
  set_info = isempty(selected_meas) ? String[] : _webui_measurement_set_comments(joinpath(directory, selected_meas))
  # the case-binding comment renders as its own line, not inside the table
  set_case = isempty(selected_meas) ? "" : get(meas_case, selected_meas, "")
  set_info = [l for l in set_info if !startswith(l, "case:")]
  # generator v2 provenance: show the summary lines compactly, drop the
  # bulky per-row truth_value block, and keep the taps table at the front
  # so the structured table renderer still recognizes it
  # `noise:` belongs to the provenance list, and it must be filtered out of
  # set_info like every other provenance line: the taps table has to stay the
  # FIRST entry there, or the structured renderer does not recognize it.
  set_provenance = [l for l in set_info if startswith(l, "generator:") || startswith(l, "seed:") || startswith(l, "noise:") || startswith(l, "source:") || startswith(l, "truth:") || startswith(l, "flow_ends:") || startswith(l, "passive:")]
  n_flow_choices = count(l -> startswith(l, "flow_end,"), set_info)
  n_flow_choices > 0 && push!(set_provenance, "flow-end choices documented for $(n_flow_choices) branch(es)")
  set_info = [l for l in set_info if !(startswith(l, "generator:") || startswith(l, "seed:") || startswith(l, "noise:") || startswith(l, "source:") || startswith(l, "truth:") || startswith(l, "flow_ends:") || startswith(l, "flow_end,") || startswith(l, "passive:") || startswith(l, "truth_value,"))]
  set_counts = isempty(selected_meas) ? Tuple{String,Int}[] : _webui_measurement_set_counts(joinpath(directory, selected_meas))
  # inline editor for small sets: the file exactly as stored, editable and
  # saved back verbatim (atomic); large files stay download-only
  set_text = ""
  if !isempty(selected_meas)
    p = joinpath(directory, selected_meas)
    isfile(p) && filesize(p) <= 262_144 && (set_text = read(p, String))
  end
  # structured value editor rows (voltages, flows, balances, everything)
  set_rows = isempty(selected_meas) ? NamedTuple[] : _webui_measurement_set_rows(joinpath(directory, selected_meas))
  # sticky generator inputs: the generate redirect carries the submitted
  # values back as g_* query keys so the form does not snap to defaults.
  # Below the stickies, the case sidecar restores what was last saved for
  # this case (generator options persist on generate, SE run options on
  # "save settings"); precedence: sticky > sidecar > default.
  gen_values = Dict{String,String}()
  se_profile = Dict{String,Any}()
  if !isempty(selected)
    loaded = _webui_load_case_settings(output_root, joinpath(directory, selected); case_directory = directory)
    loaded isa AbstractDict && (se_profile = loaded)
    for (k, v) in se_profile
      startswith(k, "_") && continue
      gen_values[k] = _webui_form_string(v)
    end
  end
  for (k, v) in query
    startswith(String(k), "g_") && (gen_values[String(k)[3:end]] = String(v))
  end
  # measurement generator v2: successful PF/SE runs of the selected case as
  # candidates for the 'from run' truth state (newest first)
  truth_runs = NamedTuple[]
  if !isempty(selected)
    runidx = load_powerflow_run_index(output_root)
    for e in get(runidx, "runs", Any[])
      e isa AbstractDict || continue
      get(e, "success", false) == true || continue
      String(get(e, "run_mode", "")) in ("", "se") || continue
      basename(String(get(e, "casefile", ""))) == basename(selected) || continue
      push!(truth_runs, (run_id = String(get(e, "run_id", "")), kind = String(get(e, "run_mode", "")) == "se" ? "se" : "pf", timestamp = String(get(e, "timestamp", ""))))
    end
    sort!(truth_runs; by = r -> r.timestamp, rev = true)
  end
  # a case file brings its measurements along; the section then offers the
  # run without any CSV in the cache
  return (cases = cases, measurements = measurements, case_measurement_rows = case_measurement_rows, case_measurement_noise = case_measurement_noise, selected_case = selected, selected_measurement = selected_meas, config_file = config_file, message = message, set_info = set_info, set_provenance = set_provenance, set_counts = set_counts, meas_case = meas_case, starred = starred, set_case = set_case, set_text = set_text, set_rows = set_rows, gen_values = gen_values, truth_runs = truth_runs)
end

"""
    handle_se_measurement_save(form; output_root, application_root, case_directory, operation_log) -> SparlectraWebUIResponse

POST /stateestimation/measurements/save: the inline editor of the SE page.
Overwrites an EXISTING measurement CSV in the case cache with the posted
content, atomically (temp file + rename). Containment mirrors the download:
basename only, `.csv`, target must already exist in the case directory and
be a measurement CSV; the posted content must itself start with the
`# sparlectra-measurements v1` marker and stays below 2 MB. Line endings
are normalized to `\\n` (browsers post textarea content with CRLF).
"""
function handle_se_measurement_save(form::AbstractDict; output_root::AbstractString, application_root::AbstractString = _webui_application_root(), case_directory = nothing, operation_log::AbstractString = output_root)
  directory = _webui_case_directory(; case_directory = case_directory, application_root = application_root, output_root = output_root)
  requested = String(_webui_form_value(form, "file", ""))
  casesel = String(_webui_form_value(form, "case", ""))
  redirectq(msg) = _webui_se_redirect(casesel, msg)
  name = basename(requested)
  (name == requested && !isempty(name) && lowercase(splitext(name)[2]) == ".csv") || return redirectq("invalid measurement file name")
  path = joinpath(directory, name)
  (isfile(path) && _webui_is_measurement_csv(path)) || return redirectq("measurement file not found: $(name)")
  content = replace(String(_webui_form_value(form, "content", "")), "\r\n" => "\n")
  sizeof(content) <= 2_097_152 || return redirectq("content too large for the inline editor; upload the file instead")
  startswith(strip(content), "# sparlectra-measurements v1") || return redirectq("content rejected: the first line must be '# sparlectra-measurements v1'")
  endswith(content, "\n") || (content *= "\n")
  tmp = path * ".tmp"
  write(tmp, content)
  mv(tmp, path; force = true)
  record_webui_operation!(operation_log, "se_measurements_saved"; route = "/stateestimation/measurements/save", method = "POST", user_action = true, measurement_file = name)
  return redirectq("saved $(name)")
end

"""
    handle_se_measurement_download(query; output_root, application_root, case_directory) -> SparlectraWebUIResponse

GET /stateestimation/measurements/download?file=<name>: hand the named
measurement CSV from the case cache to the browser as an attachment. Path
containment: basename only, `.csv`, must live in the case directory and
pass the measurement content sniff.
"""
function handle_se_measurement_download(query::AbstractDict; output_root::AbstractString, application_root::AbstractString = _webui_application_root(), case_directory = nothing)
  directory = _webui_case_directory(; case_directory = case_directory, application_root = application_root, output_root = output_root)
  requested = String(get(query, "file", ""))
  isempty(requested) && return _webui_html(render_webui_error(400, "Missing file parameter."); status = 400)
  name = basename(requested)
  (name == requested && lowercase(splitext(name)[2]) == ".csv") || return _webui_html(render_webui_error(400, "Invalid measurement file name."); status = 400)
  path = joinpath(directory, name)
  (isfile(path) && _webui_is_measurement_csv(path)) || return _webui_html(render_webui_error(404, "Measurement file not found."); status = 404)
  headers = ["Content-Type" => "text/csv; charset=utf-8", "Content-Disposition" => "attachment; filename=\"$(name)\""]
  return SparlectraWebUIResponse(200, headers, read(path))
end

"""
    handle_se_measurement_update_values(form; output_root, application_root, case_directory, operation_log) -> SparlectraWebUIResponse

POST /stateestimation/measurements/update-values: the structured table
editor of the SE page. Rewrites ONLY the value/sigma/active columns of the
addressed data rows (keyed by their file line number), leaving comments,
header, ordering, and all other columns untouched; the write is atomic.
Containment mirrors the text editor: basename, `.csv`, existing measurement
CSV in the case directory. Rejected rows (unparsable value, non-positive
sigma) reject the whole save so a typo never half-applies.
"""
function handle_se_measurement_update_values(form::AbstractDict; output_root::AbstractString, application_root::AbstractString = _webui_application_root(), case_directory = nothing, operation_log::AbstractString = output_root)
  directory = _webui_case_directory(; case_directory = case_directory, application_root = application_root, output_root = output_root)
  requested = String(_webui_form_value(form, "file", ""))
  casesel = String(_webui_form_value(form, "case", ""))
  redirectq(msg) = _webui_se_redirect(casesel, msg)
  name = basename(requested)
  (name == requested && !isempty(name) && lowercase(splitext(name)[2]) == ".csv") || return redirectq("invalid measurement file name")
  path = joinpath(directory, name)
  (isfile(path) && _webui_is_measurement_csv(path)) || return redirectq("measurement file not found: $(name)")
  lines = readlines(path)
  changed = 0
  for (ln, line) in enumerate(lines)
    v = _webui_form_value(form, "v_$(ln)", nothing)
    v === nothing && continue
    s = String(_webui_form_value(form, "s_$(ln)", ""))
    a = String(_webui_form_value(form, "a_$(ln)", "true"))
    parts = split(String(line), ","; limit = 11)
    length(parts) >= 11 || return redirectq("row at line $(ln) is not editable")
    vf = tryparse(Float64, String(v))
    sf = tryparse(Float64, s)
    vf === nothing && return redirectq("line $(ln): value '$(v)' is not a number; nothing was saved")
    (sf === nothing || sf <= 0.0) && return redirectq("line $(ln): sigma '$(s)' must be a positive number; nothing was saved")
    a in ("true", "false") || return redirectq("line $(ln): active must be true or false; nothing was saved")
    if String(parts[8]) != String(v) || String(parts[9]) != s || String(parts[10]) != a
      parts[8] = v
      parts[9] = s
      parts[10] = a
      lines[ln] = join(parts, ",")
      changed += 1
    end
  end
  changed == 0 && return redirectq("no value changed")
  tmp = path * ".tmp"
  open(tmp, "w") do io
    for l in lines
      println(io, l)
    end
  end
  mv(tmp, path; force = true)
  record_webui_operation!(operation_log, "se_measurements_values_updated"; route = "/stateestimation/measurements/update-values", method = "POST", user_action = true, measurement_file = name, changed = changed)
  return redirectq("updated $(changed) measurement value(s) in $(name)")
end

"""
    handle_se_topology_hypotheses(form; output_root, application_root, case_directory, operation_log) -> SparlectraWebUIResponse

POST /stateestimation/topology-hypotheses: stage-3 topology hypothesis test
for a finished SE run, behind an explicit button (never automatic). Re-runs
the estimation per candidate status toggle on WORKING COPIES, writes the
ranked recommendations to `topology_hypotheses.md` in the run directory,
and returns to the result page (which renders the table inline). Nothing
is ever switched: the output is a recommendation list.
"""
function handle_se_topology_hypotheses(form::AbstractDict; output_root::AbstractString, application_root::AbstractString = _webui_application_root(), case_directory = nothing, operation_log::AbstractString = output_root)
  run_id = basename(String(_webui_form_value(form, "run_id", "")))
  back(msg) = _webui_redirect("/powerflow/result/$(_webui_urlencode(run_id))?message=$(_webui_urlencode(msg))")
  isempty(run_id) && return _webui_html(render_webui_error(400, "Missing run id."); status = 400)
  result = get_powerflow_result(run_id)
  metadata = get(result, "metadata", Dict{String,Any}())
  (get(result, "success", false) && get(metadata, "run_mode", "") == "se") || return back("topology hypothesis test needs a successful state-estimation run")
  outdir = String(get(result, "output_dir", ""))
  mf = joinpath(outdir, "measurements.csv")
  isfile(mf) || return back("measurements.csv of this run is not available any more")
  t0 = time_ns()
  rep = try
    config = load_sparlectra_config(String(get(result, "config_file", DEFAULT_SPARLECTRA_CONFIG_PATH)); reload = true)
    # same format hint the run form uses: a bare .DAT is ambiguous and the
    # detector refuses to guess, so the DTF network case would fail here
    se_case = String(get(result, "casefile", ""))
    net = _se_import_case_net(se_case, config; requested_format = _webui_case_format_hint(se_case))
    readMeasurementsCSV!(net; file = mf, replace = true)
    test_topology_hypotheses(net; max_candidates = 5, maxIte = 100, tol = 1e-8)
  catch err
    return back("topology hypothesis test failed: $(sprint(showerror, err))")
  end
  open(joinpath(outdir, "topology_hypotheses.md"), "w") do io
    println(io, "# Topology hypothesis test (advisory)\n")
    println(io, "Base run: J = ", round(rep.base_j; sigdigits = 4), " (dof ", rep.base_dof, ", band ", rep.base_verdict, "). ", rep.n_candidates, " candidate(s) tested.")
    rep.ambiguous && println(io, "\nSEVERAL hypotheses land in the band: the case is ambiguous, ranked by J below.")
    println(io, "\n| element | current status | hypothesis | J before | J after | z before | z after | verdict | time s |")
    println(io, "|---|---|---|---|---|---|---|---|---|")
    for r in rep.recommendations
      println(io, "| ", r.element, " | ", r.current_status, " | ", r.hypothesis, " | ", round(r.j_before; sigdigits = 4), " | ", isnan(r.j_after) ? "-" : round(r.j_after; sigdigits = 4), " | ", round(r.z_before; digits = 2), " | ", isnan(r.z_after) ? "-" : round(r.z_after; digits = 2), " | ", r.verdict, " | ", round(r.elapsed_s; digits = 2), " |")
    end
    println(io, "\nRecommendations only: NOTHING has been switched. A supported hypothesis means the band failure disappears when the element's status is assumed flipped; verify against the real switch state before touching the model.")
  end
  open(joinpath(outdir, "run.log"), "a") do io
    println(io, "topology hypothesis test: ", rep.n_candidates, " candidate(s) in ", round((time_ns() - t0) / 1e9; digits = 2), " s, ", count(r -> r.verdict == :hypothesis_supported, rep.recommendations), " supported", rep.ambiguous ? " (ambiguous)" : "", "; topology_hypotheses.md written")
  end
  record_webui_operation!(operation_log, "se_topology_hypotheses"; route = "/stateestimation/topology-hypotheses", method = "POST", user_action = true, run_id = run_id)
  sup = count(r -> r.verdict == :hypothesis_supported, rep.recommendations)
  return back(string("topology hypothesis test finished: ", rep.n_candidates, " candidate(s), ", sup, " supported", rep.ambiguous ? " (ambiguous, ranked by J)" : "", "; table below"))
end

"""
    handle_se_generate_measurements(form; output_root, application_root, case_directory, operation_log) -> SparlectraWebUIResponse

POST /stateestimation/generate-measurements: demo action. Solves the selected
case (MATPOWER or CGMES, same import paths as the SE run) once and writes a
measurement CSV v1 (`<case>.measurements.csv`) into the case cache, then
redirects back to the SE page with the file offered. Options: per-quantity
sigmas in percent of the measured value (voltage-level independent,
`relativeSigma` generation with the `measurementSigmaFloors` floors);
`noise = true` adds seeded Gaussian noise at those sigmas;
`gross_error_k > 0` additionally corrupts `gross_error_count` seed-randomly
drawn telemetry rows by k times their sigma (a reproducible bad-data test
vector for the elimination/robust workflow); `tap_error_steps` shifts up to
`tap_error_count` seed-randomly drawn estimable transformers by that many
whole mechanical steps for the generation state.
"""
function handle_se_generate_measurements(form::AbstractDict; output_root::AbstractString, application_root::AbstractString = _webui_application_root(), case_directory = nothing, operation_log::AbstractString = output_root)
  directory = _webui_case_directory(; case_directory = case_directory, application_root = application_root, output_root = output_root)
  casefile = String(_webui_form_value(form, "casefile", ""))
  # sticky inputs: every generator value travels back in the redirect as a
  # g_* query key so the re-rendered form keeps what the user typed
  gq = join(("&g_$(k)=$(_webui_urlencode(String(_webui_form_value(form, k, ""))))" for k in ("noise", "gross_error_k", "gross_error_count", "tap_error_steps", "tap_error_count", "sigma_u_pct", "include_currents", "sigma_i_pct", "sigma_ia_deg", "sigma_p_pct", "sigma_q_pct", "gen_truth_source", "gen_truth_run_id", "gen_flow_ends", "gen_passive_sigma", "gen_passive_as_zi", "gen_seed")), "")
  redirectq(msg) = _webui_se_redirect(casefile, msg; extra_query = gq)
  isempty(casefile) && return redirectq("no case selected")
  case_path = joinpath(directory, basename(casefile))
  isfile(case_path) || return redirectq("case not found: $(casefile)")
  noise = String(_webui_form_value(form, "noise", "")) == "true"
  gross_k = something(tryparse(Float64, String(_webui_form_value(form, "gross_error_k", "0"))), 0.0)
  (isfinite(gross_k) && gross_k >= 0.0) || return redirectq("gross error k must be a non-negative number")
  # how many rows get the gross error (the seed decides which rows)
  gross_count = something(tryparse(Int, String(_webui_form_value(form, "gross_error_count", "1"))), 0)
  1 <= gross_count <= 1000 || return redirectq("bad data count must be an integer between 1 and 1000")
  # per-quantity measurement sigmas in PERCENT OF THE MEASURED VALUE (class
  # accuracy), so one setting stays meaningful across voltage levels: an
  # absolute 2 MW would be tight on a 400 kV flow and absurd on a 30 kV
  # feeder, 1 percent of reading fits both. Small per-type floors (see
  # measurementSigmaFloors) guard near-zero readings.
  include_i = String(_webui_form_value(form, "include_currents", "")) == "true"
  sigma_u_pct = something(tryparse(Float64, String(_webui_form_value(form, "sigma_u_pct", "0.5"))), NaN)
  sigma_i_pct = something(tryparse(Float64, String(_webui_form_value(form, "sigma_i_pct", "1.0"))), NaN)
  sigma_p_pct = something(tryparse(Float64, String(_webui_form_value(form, "sigma_p_pct", "1.0"))), NaN)
  sigma_q_pct = something(tryparse(Float64, String(_webui_form_value(form, "sigma_q_pct", "1.0"))), NaN)
  (isfinite(sigma_u_pct) && sigma_u_pct > 0.0) || return redirectq("sigma U must be a positive percentage")
  (!include_i || (isfinite(sigma_i_pct) && sigma_i_pct > 0.0)) || return redirectq("sigma I must be a positive percentage")
  (isfinite(sigma_p_pct) && sigma_p_pct > 0.0) || return redirectq("sigma P must be a positive percentage")
  (isfinite(sigma_q_pct) && sigma_q_pct > 0.0) || return redirectq("sigma Q must be a positive percentage")
  # PMU current ANGLES: absolute degrees (an angle passes through zero, a
  # percentage of it is meaningless); 0 = no angle rows
  sigma_ia_deg = something(tryparse(Float64, String(_webui_form_value(form, "sigma_ia_deg", "0"))), NaN)
  (isfinite(sigma_ia_deg) && sigma_ia_deg >= 0.0) || return redirectq("sigma Ia must be zero (no current-angle rows) or positive degrees")
  # tap deviation: the measurements are generated from a state whose first
  # in-service transformer runs tap_error_steps MECHANICAL steps off the
  # model position, on the same tap-fraction grid the estimator fixes to
  # (r = n * tap_step). The model file is never touched (throwaway import).
  tap_steps = something(tryparse(Float64, String(_webui_form_value(form, "tap_error_steps", "0"))), NaN)
  (isfinite(tap_steps) && abs(tap_steps) <= 16.0) || return redirectq("tap deviation must be a number of steps between -16 and 16")
  # whole steps only: a tap changer has no half
  # positions, and a fractional deviation would leave an off-grid residual
  # that no fixation can remove
  abs(tap_steps - round(tap_steps)) < 1e-9 || return redirectq("tap deviation must be a whole number of steps (a tap changer has no half positions)")
  # how many transformers get the deviation AT MOST (fewer when the case
  # has fewer eligible transformers; the seed decides which)
  tap_count = something(tryparse(Int, String(_webui_form_value(form, "tap_error_count", "1"))), 0)
  1 <= tap_count <= 100 || return redirectq("tap deviation transformer count must be an integer between 1 and 100")
  # measurement generator v2: truth state, flow ends, passive nodes
  truth_source_raw = lowercase(strip(String(_webui_form_value(form, "gen_truth_source", "fresh_solve"))))
  truth_source_raw in ("fresh_solve", "from_run") || return redirectq("truth state must be fresh_solve or from_run")
  truth_source = Symbol(truth_source_raw)
  truth_run_id = String(strip(String(_webui_form_value(form, "gen_truth_run_id", ""))))
  truth_source === :from_run && isempty(truth_run_id) && return redirectq("truth state 'from run' needs a run id")
  # server-side lock behind the GUI graying: a run state is a finished
  # snapshot and is never re-solved with a shifted tap
  truth_source === :from_run && tap_steps != 0.0 && return redirectq("tap deviation requires truth state 'fresh solve'")
  flow_ends_raw = lowercase(strip(String(_webui_form_value(form, "gen_flow_ends", "both"))))
  flow_ends_raw in ("both", "one_balance_aware") || return redirectq("flow measurements per branch must be both or one_balance_aware")
  flow_ends = Symbol(flow_ends_raw)
  passive_sigma = something(tryparse(Float64, String(_webui_form_value(form, "gen_passive_sigma", "0.05"))), NaN)
  (isfinite(passive_sigma) && passive_sigma > 0.0) || return redirectq("passive-node balance sigma must be positive (MW/MVar)")
  passive_as_zi = _webui_parse_bool(String(_webui_form_value(form, "gen_passive_as_zi", "false")))
  # reproducibility seed: the same seed regenerates the identical set, a
  # different one draws a fresh noise realization
  gen_seed = something(tryparse(Int, String(_webui_form_value(form, "gen_seed", "42"))), -1)
  gen_seed >= 0 || return redirectq("seed must be a non-negative integer")
  out_name = string(splitext(basename(casefile))[1], ".measurements.csv")
  out_path = joinpath(directory, out_name)
  # the actual import/solve/generation lives in the SE service layer
  # (_se_generate_measurement_set): the Web UI layer never calls a solver
  # directly, the same architecture rule the run handlers follow
  gen = try
    _se_generate_measurement_set(case_path, out_path; noise = noise, gross_k = gross_k, gross_count = gross_count, tap_steps = tap_steps, tap_count = tap_count, include_i = include_i, sigma_u_pct = sigma_u_pct, sigma_i_pct = sigma_i_pct, sigma_p_pct = sigma_p_pct, sigma_q_pct = sigma_q_pct, sigma_ia_deg = sigma_ia_deg, truth_source = truth_source, run_id = truth_run_id, run_root = output_root, flow_ends = flow_ends, passive_sigma = passive_sigma, passive_as_zi = passive_as_zi, seed = gen_seed)
  catch err
    err isa ArgumentError && return redirectq(err.msg)
    return redirectq("generation failed: $(sprint(showerror, err))")
  end
  record_webui_operation!(operation_log, "se_measurements_generated"; route = "/stateestimation/generate-measurements", method = "POST", user_action = true, casefile = casefile, measurement_file = out_name)
  # persist the generator options in the case sidecar so a case reload
  # restores them (the run-based save path never sees this form)
  try
    _webui_merge_case_settings!(output_root, case_path, Dict{String,Any}("noise" => noise, "gross_error_k" => gross_k, "gross_error_count" => gross_count, "tap_error_steps" => tap_steps, "tap_error_count" => tap_count, "sigma_u_pct" => sigma_u_pct, "include_currents" => include_i, "sigma_i_pct" => sigma_i_pct, "sigma_ia_deg" => sigma_ia_deg, "sigma_p_pct" => sigma_p_pct, "sigma_q_pct" => sigma_q_pct, "gen_truth_source" => truth_source_raw, "gen_flow_ends" => flow_ends_raw, "gen_passive_sigma" => passive_sigma, "gen_passive_as_zi" => passive_as_zi, "gen_seed" => gen_seed); case_directory = directory)
  catch err
    record_webui_operation!(operation_log, "case_settings_save_failed"; route = "/stateestimation/generate-measurements", method = "POST", user_action = true, message = sprint(showerror, err))
  end
  sigma_note = ", sigma U=$(sigma_u_pct)% P=$(sigma_p_pct)% Q=$(sigma_q_pct)%$(include_i ? " I=$(sigma_i_pct)%" : "") of reading"
  return redirectq("generated $(out_name) ($(gen.rows) rows, $(gen.noisy ? "noisy" : "noise-free")$(sigma_note)$(gen.truth_note)$(gen.flow_note)$(gen.passive_note)$(gen.gross_note)$(gen.tap_note)$(gen.island_note))")
end

"""
    handle_se_add_noise(form; output_root, application_root, case_directory, operation_log)

Perturb an existing measurement set with Gaussian noise and store the result
as a set bound to the selected case. No power flow is computed and no truth
state is needed: the values that are there are perturbed with the sigma each
row declares.

This exists for the set a case file delivers, which is typically noise-free.
An estimation on ideal values returns `J = 0` by construction, and until now
the only way to a realistic run was regenerating the whole set from a fresh
solve - which also replaces the values themselves. Adding noise keeps the
operating point and only makes the readings imperfect.
"""
function handle_se_add_noise(form::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, operation_log::AbstractString = output_root)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory = case_directory, application_root = application_root, output_root = output_root)
  casefile = strip(String(something(_webui_form_value(form, "casefile", ""), "")))
  back(msg) = _webui_se_redirect(casefile, msg)
  isempty(casefile) && return back("Select a case before adding noise to its measurements.")
  (basename(casefile) == casefile && !occursin(r"[\\/]", casefile)) || return back("Invalid case name.")
  case_path = joinpath(directory, casefile)
  isfile(case_path) || return back("Case file not found: $(casefile)")
  source = strip(String(something(_webui_form_value(form, "measurement_file", ""), "")))
  seed = something(tryparse(Int, strip(String(something(_webui_form_value(form, "noise_seed", ""), "")))), 42)
  out_name = string(first(splitext(basename(casefile))), ".noisy.measurements.csv")
  out_path = joinpath(directory, out_name)
  rows = try
    config = load_sparlectra_config(DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
    net = _se_import_case_net(case_path, config; requested_format = _webui_case_format_hint(case_path))
    if !isempty(source)
      empty!(net.measurements)
      readMeasurementsCSV!(net; file = joinpath(directory, source), replace = true)
    end
    isempty(net.measurements) && return back("Neither the case file nor a selected set carries measurements to perturb.")
    addMeasurementNoise!(net; rng = Random.MersenneTwister(seed))
    writeMeasurementsCSV(net; file = out_path,
      headerComments = String["case: $(casefile)", "generator: noise-only (existing values perturbed, no power flow)",
                              "seed: $(seed)", "noise: gaussian (each row's own sigma)",
                              "source: $(isempty(source) ? "measurements carried by the case file" : source)"],
      busReference = _detect_case_format(case_path) === :cgmes ? :mrid : :name)
    length(net.measurements)
  catch err
    record_webui_operation!(operation_log, "se_measurements_noised_failed"; route = "/stateestimation/add-noise", method = "POST", user_action = true, casefile = casefile, message = sprint(showerror, err))
    return back("Could not add noise: $(first(split(sprint(showerror, err), '\n')))")
  end
  record_webui_operation!(operation_log, "se_measurements_noised"; route = "/stateestimation/add-noise", method = "POST", user_action = true, casefile = casefile, measurement_file = out_name, rows = rows, seed = seed)
  return back("Wrote $(out_name): $(rows) rows, each perturbed with its own sigma (seed $(seed)). No power flow was computed.")
end

"""
POST /stateestimation/reset-settings: deletes the saved case settings
(the case configuration file, plus any leftover legacy sidecar) of the
selected case, so the generator and run forms fall back to their defaults
on the next render. The reset covers the whole per-case profile; the
PowerFlow page shares the same file. The redirect intentionally carries NO sticky `g_*`
values, so stale form values do not survive the reset.
"""
function handle_se_reset_settings(form::AbstractDict; output_root::AbstractString, application_root::AbstractString = _webui_application_root(), case_directory = nothing, operation_log::AbstractString = output_root)
  directory = _webui_case_directory(; case_directory = case_directory, application_root = application_root, output_root = output_root)
  casefile = String(_webui_form_value(form, "casefile", ""))
  redirectq(msg) = _webui_se_redirect(casefile, msg)
  isempty(casefile) && return redirectq("no case selected")
  case_path = joinpath(directory, basename(casefile))
  isfile(case_path) || return redirectq("case not found: $(casefile)")
  removed_any = false
  try
    path = _webui_case_settings_path(output_root, case_path; case_directory = directory)
    if isfile(path)
      rm(path)
      removed_any = true
    end
    legacy = _webui_legacy_case_settings_path(output_root, case_path; case_directory = directory)
    if isfile(legacy)
      rm(legacy; force = true)
      removed_any = true
    end
  catch err
    return redirectq("settings reset failed: $(sprint(showerror, err))")
  end
  removed_any || return redirectq("no saved settings for $(basename(casefile)); the forms already use the defaults")
  record_webui_operation!(operation_log, "se_case_settings_reset"; route = "/stateestimation/reset-settings", method = "POST", user_action = true, casefile = casefile)
  return redirectq("saved settings for $(basename(casefile)) deleted; the forms are back to the defaults")
end

"""
    handle_webui_sysimage(; output_root, message) -> SparlectraWebUIResponse

Render the sysimage page (state of the image, a running build, the refresh
button). Pure rendering, no side effects, so the auto-refresh poll while a
build runs costs nothing but a file read.
"""
function handle_webui_sysimage(; output_root::AbstractString, message::AbstractString = "")::SparlectraWebUIResponse
  return _webui_html(render_webui_sysimage_page(; output_root, message))
end

"""
    handle_webui_sysimage_rebuild(; output_root, operation_log) -> SparlectraWebUIResponse

Start a sysimage rebuild in the background and go back to the sysimage page
with the outcome as a message. The redirect matters: a POST that renders its
own page turns the browser's reload button into a second build.
"""
function handle_webui_sysimage_rebuild(; output_root::AbstractString, operation_log::AbstractString)::SparlectraWebUIResponse
  result = start_sysimage_rebuild!(; output_root)
  record_webui_operation!(operation_log, "sysimage_rebuild_requested"; route = "/webui/sysimage/rebuild", method = "POST", user_action = true, status = result.started ? "started" : "rejected", message = result.message)
  return _webui_redirect("/webui/sysimage?message=$(_webui_urlencode(result.message))")
end

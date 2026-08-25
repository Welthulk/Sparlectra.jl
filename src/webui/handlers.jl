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
      role = ext == ".dat" ? _webui_dat_role_label(_webui_classify_dat_content(destination)) : ext == ".zip" ? _webui_cgmes_upload_role(destination, directory) : "matpower_case"
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
  query = isempty(selected) ? "" : "?casefile=$(_webui_urlencode(selected))"
  display_imported = ["$(name) ($(get(imported_roles, name, "unknown")))" for name in imported]
  message = _webui_urlencode(_webui_case_import_message(display_imported, rejected))
  separator = isempty(query) ? "?" : "&"
  return _webui_redirect("/powerflow$(query)$(separator)import_message=$(message)")
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
    return _webui_redirect("/powerflow?import_message=$(message)")
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
    return _webui_redirect("/powerflow?casefile=$(_webui_urlencode(resolved_name))&import_message=$(message)")
  end
  record_webui_operation!(operation_log, "case_resolve_failed"; route = "/powerflow/resolve-case", method = "POST", user_action = true, requested = manual_value, message = error_text)
  message = _webui_urlencode("Could not resolve case '$(manual_value)': $(error_text)")
  return _webui_redirect("/powerflow?import_message=$(message)")
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
    return _webui_redirect("/powerflow?import_message=$(message)")
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
        rm(target)
        sidecar = _webui_case_settings_path(output_root, requested; case_directory = directory)
        isfile(sidecar) && rm(sidecar; force = true)
        # the per-case N-1 weight list travels with the case (issue #331 Phase 5
        # follow-up); remove it too so a deleted case leaves no orphan behind
        weights = _webui_case_weights_path(requested; case_directory = directory)
        isfile(weights) && rm(weights; force = true)
      catch err
        error_text = sprint(showerror, err)
      end
    end
  end
  if isempty(error_text)
    record_webui_operation!(operation_log, "case_delete_completed"; route = "/powerflow/delete-case", method = "POST", user_action = true, casefile = requested)
    message = _webui_urlencode("Deleted case: $(requested)")
    return _webui_redirect("/powerflow?import_message=$(message)")
  end
  record_webui_operation!(operation_log, "case_delete_failed"; route = "/powerflow/delete-case", method = "POST", user_action = true, casefile = requested, message = error_text)
  message = _webui_urlencode("Could not delete case '$(requested)': $(error_text)")
  return _webui_redirect("/powerflow?import_message=$(message)")
end

"""
Delete a case's saved settings sidecar without touching the case file. Saved
settings outrank the configuration for their keys, so a stale sidecar can pin
a case to settings the user cannot see or override in the form (measured: a
delivery stuck on `power_flow_solver: dc` from an old run). This is the way
back to the plain configuration defaults.
"""
function handle_powerflow_case_settings_reset(form::AbstractDict; output_root::AbstractString = "results/powerflow_service", application_root::AbstractString = _webui_application_root(), case_directory::Union{Nothing,AbstractString} = nothing, operation_log::AbstractString = output_root)::SparlectraWebUIResponse
  directory = _webui_case_directory(; case_directory, application_root, output_root)
  requested = strip(String(something(_webui_form_value(form, "casefile", ""), "")))
  if isempty(requested)
    return _webui_redirect("/powerflow?import_message=$(_webui_urlencode("Select a case before resetting its saved settings."))")
  end
  if basename(requested) != requested || occursin(r"[\\/]", requested)
    record_webui_operation!(operation_log, "case_settings_reset_failed"; route = "/powerflow/case-settings/reset", method = "POST", user_action = true, casefile = requested, message = "invalid case name")
    return _webui_redirect("/powerflow?import_message=$(_webui_urlencode("Could not reset settings for '$(requested)': invalid case name"))")
  end
  sidecar = _webui_case_settings_path(output_root, requested; case_directory = directory)
  if !isfile(sidecar)
    record_webui_operation!(operation_log, "case_settings_reset_noop"; route = "/powerflow/case-settings/reset", method = "POST", user_action = true, casefile = requested, status = "no_sidecar")
    return _webui_redirect("/powerflow?casefile=$(_webui_urlencode(requested))&import_message=$(_webui_urlencode("No saved settings for '$(requested)' — the form already uses the configuration defaults."))")
  end
  try
    rm(sidecar; force = true)
  catch err
    record_webui_operation!(operation_log, "case_settings_reset_failed"; route = "/powerflow/case-settings/reset", method = "POST", user_action = true, casefile = requested, message = sprint(showerror, err))
    return _webui_redirect("/powerflow?casefile=$(_webui_urlencode(requested))&import_message=$(_webui_urlencode("Could not reset settings for '$(requested)': $(sprint(showerror, err))"))")
  end
  record_webui_operation!(operation_log, "case_settings_reset_completed"; route = "/powerflow/case-settings/reset", method = "POST", user_action = true, casefile = requested)
  return _webui_redirect("/powerflow?casefile=$(_webui_urlencode(requested))&import_message=$(_webui_urlencode("Saved settings for '$(requested)' deleted — the form now uses the configuration defaults."))")
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

function _webui_case_settings_saved_html(path::AbstractString, casefile::AbstractString, count::Integer, successful::Bool, override::Bool)::String
  status = successful ? "successful/converged" : (override ? "non-successful, saved via explicit override" : "non-successful")
  body = """
<section class=\"panel case-settings-saved\"><h2>Case settings saved</h2>
<p class=\"alert info\"><strong>Saved sidecar:</strong> <code>$(_webui_escape(path))</code></p>
<ul>
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
    return _webui_html(_webui_case_settings_saved_html(path, case_settings_source, length(settings), successful, !successful && override))
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

"""Serve the fast-start (sysimage) status page."""
handle_webui_fast_start(output_root::AbstractString)::SparlectraWebUIResponse = _webui_html(render_webui_fast_start(output_root))

"""Serve the sysimage build-log page."""
handle_webui_fast_start_log(output_root::AbstractString)::SparlectraWebUIResponse = _webui_html(render_webui_fast_start_log(output_root))

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

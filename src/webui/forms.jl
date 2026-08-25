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

# file: src/webui/forms.jl
# purpose: Web UI input resolution: application root and case directory
#          lookup, case/config file selectors, upload classification, and
#          form value parsing
const _WEBUI_PACKAGE_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

"""
    _webui_application_root([start_dir]) -> String

Resolve the Sparlectra application root used by the local Web UI. The lookup
supports starting Julia from the repository itself, from a parent directory
containing `Sparlectra` or `Sparlectra.jl`, and from the installed package root.
"""
function _webui_application_root(start_dir::AbstractString = pwd())::String
  start_root = abspath(start_dir)
  candidates = unique((
    start_root,
    joinpath(start_root, "Sparlectra"),
    joinpath(start_root, "Sparlectra.jl"),
    _WEBUI_PACKAGE_ROOT,
  ))
  for candidate in candidates
    isdir(joinpath(candidate, "data", "mpower")) && isdir(joinpath(candidate, "examples")) && return normpath(candidate)
  end
  return _WEBUI_PACKAGE_ROOT
end

"""
    _webui_casefile_options(application_root) -> Vector{String}

Return sorted user-selectable case filenames from the Web UI application's
`data/mpower` directory. The selector stays conservative: MATPOWER `.m` files
and copied internal DTF `.DAT` candidates are shown, while generated
Julia cache artifacts, warm-up cases, result artifacts, and sidecar profiles
stay hidden. Missing or empty directories produce an empty list.
"""
function _webui_casefile_options(application_root::AbstractString)::Vector{String}
  return _webui_casefile_options_in_directory(_webui_case_directory(; application_root))
end

"""Return the effective Web UI case directory used for selection and imports."""
function _webui_case_directory(; case_directory::Union{Nothing,AbstractString} = nothing, application_root::AbstractString = _webui_application_root(), output_root::AbstractString = default_webui_output_root())::String
  case_directory === nothing || return normpath(String(case_directory))
  package_case_directory = joinpath(application_root, "data", "mpower")
  try
    mkpath(package_case_directory)
    test_path = tempname(package_case_directory)
    open(test_path, "w") do io
      write(io, "")
    end
    rm(test_path; force = true)
    return normpath(package_case_directory)
  catch
    return normpath(default_webui_case_cache_dir(output_root))
  end
end

# .zip covers CGMES deliveries (single ZIP, base case plus boundary, or
# ZIP-in-ZIP) — the importer's container layer opens them in memory.
_webui_supported_upload_case_extension(name::AbstractString)::Bool = lowercase(splitext(basename(String(name)))[2]) in (".m", ".dat", ".zip")

function _webui_classify_dat_content(path::AbstractString)::Symbol
  lowercase(splitext(path)[2]) == ".dat" || return :not_dat
  text = try
    read(path, String)
  catch
    return :unknown_dat
  end
  try
    case = DTFImporter.read_dtf(path; strict = false)
    if case.size.NGES > 0 && length(case.buses) == case.size.NGES && length(case.branches) == case.size.LGES
      return isempty(case.outages) ? :dtf_network_case : :dtf_network_case_with_outages
    end
  catch
  end
  has_outages = occursin(r"(?im)^\s*AUSFALL\s*$", text) && occursin(r"(?im)^\s*ENDE\s*$", text)
  if has_outages
    return :dtf_outage_file
  elseif _is_for002_reference_dat(path) || occursin(r"(?i)(for002|reference|vergleich|result)", text)
    return :dtf_outage_or_reference
  end
  return :unknown_dat
end

_webui_dat_role_label(role::Symbol)::String = replace(String(role), '_' => ' ')

"""True when a CGMES delivery contains boundary profiles only (no equipment)."""
function _webui_is_cgmes_boundary_only(path::AbstractString)::Bool
  try
    profiles = union((f.profiles for f in CGMESImporter.summarizeCGMES(path = path).files)...)
    return isempty(intersect(profiles, (:EQ, :EQ_OP, :EQ_SC, :TP, :SSH))) && !isempty(intersect(profiles, (:EQ_BD, :EQ_BD_OP, :TP_BD)))
  catch
    return false
  end
end

"""
Pre-analysis of an uploaded CGMES ZIP, shown in the import summary so the
user sees immediately what the delivery contains and whether it is ready to
compute: version, profiles, element counts, and the boundary situation.
"""
function _webui_cgmes_upload_role(path::AbstractString, directory::AbstractString)::String
  summary = try
    CGMESImporter.summarizeCGMES(path = path)
  catch err
    return "⚠️ unreadable CGMES delivery (" * first(sprint(showerror, err), 80) * ")"
  end
  profiles = union((f.profiles for f in summary.files)...)
  counts = Dict(summary.class_histogram)
  cnt(cls) = get(counts, cls, 0)
  version = isempty(summary.version) ? "CGMES" : "CGMES " * summary.version
  proflist = join(sort(String.(collect(profiles))), ",")

  # a boundary-only delivery has no equipment of its own
  if isempty(intersect(profiles, (:EQ, :EQ_OP, :EQ_SC, :TP, :SSH))) && !isempty(intersect(profiles, (:EQ_BD, :EQ_BD_OP, :TP_BD)))
    return "🔗 $(version) boundary set · $(proflist) · $(cnt(:TopologicalNode)) boundary nodes — pair it with a base case"
  end

  parts = String[]
  push!(parts, "$(cnt(:TopologicalNode)) nodes")
  cnt(:ACLineSegment) > 0 && push!(parts, "$(cnt(:ACLineSegment)) lines")
  cnt(:PowerTransformer) > 0 && push!(parts, "$(cnt(:PowerTransformer)) transformers")
  total_loads = cnt(:EnergyConsumer) + cnt(:ConformLoad) + cnt(:NonConformLoad)
  total_loads > 0 && push!(parts, "$(total_loads) loads")
  cnt(:SynchronousMachine) > 0 && push!(parts, "$(cnt(:SynchronousMachine)) machines")
  cnt(:LinearShuntCompensator) > 0 && push!(parts, "$(cnt(:LinearShuntCompensator)) shunts")
  cnt(:SvVoltage) > 0 && push!(parts, "SV reference present")
  content = join(parts, ", ")

  missing_profiles = [p for p in (:EQ, :SSH, :TP) if !(p in profiles)]
  isempty(missing_profiles) || return "❌ $(version) · $(proflist) · $(content) — incomplete: missing profile(s) $(join(String.(missing_profiles), ", "))$(_webui_upload_analysis_note(path))"

  if !summary.boundary_missing_hint
    return "✅ $(version) · $(proflist) · $(content) · self-contained, ready to compute"
  end
  # a boundary delivery next to it completes this one
  neighbours = try
    filter(n -> occursin(r"(?i)boundary|_bd_|_bd\.", n) && joinpath(directory, n) != path, readdir(directory))
  catch
    String[]
  end
  for n in neighbours
    try
      s2 = CGMESImporter.summarizeCGMES(path = [path, joinpath(directory, n)])
      s2.boundary_missing_hint || return "✅ $(version) · $(proflist) · $(content) · boundary found in $(n), ready to compute"
    catch
    end
  end
  # Last resort before bothering the user: the ENTSO-E test configurations are
  # published as base case plus a separate boundary package. If the extracted
  # local test-set cache holds a boundary that completes this delivery, copy
  # it next to the case instead of asking for a second upload.
  fetched = _webui_try_supply_cgmes_boundary(path, directory)
  isempty(fetched) || return "✅ $(version) · $(proflist) · $(content) · boundary supplied automatically ($(fetched)), ready to compute"
  return "⚠️ $(version) · $(proflist) · $(content) · $(summary.unresolved_count) unresolved references — boundary set missing: upload it as well, or type cgmes:&lt;alias&gt;$(_webui_upload_analysis_note(path))"
end

"""
Upload-time import analysis for a CGMES delivery that stayed incomplete
(missing profiles or missing boundary): writes the full
[`importFailureAnalysis`](@ref Sparlectra.CGMESImporter.importFailureAnalysis)
report next to the case as `<case>.import_analysis.txt` and returns a short
message suffix naming the missing declared dependencies. Best-effort — an
unreadable delivery yields an empty suffix and no file.
"""
function _webui_upload_analysis_note(path::AbstractString)::String
  try
    store = CGMESImporter.loadCGMES(path)
    stats = CGMESImporter.importabilityStats(store)
    report_file = string(splitext(path)[1], ".import_analysis.txt")
    open(report_file, "w") do io
      print(io, CGMESImporter.importFailureAnalysis(store))
    end
    if !isempty(stats.missing_dependencies)
      shown = join(Iterators.take(stats.missing_dependencies, 3), ", ")
      length(stats.missing_dependencies) > 3 && (shown *= ", …")
      return " · missing declared dependencies: $(shown) — full analysis in $(basename(report_file))"
    end
    return " — full analysis in $(basename(report_file))"
  catch
    # keep the role summary useful even when the analysis cannot be built
    return ""
  end
end

"""
Look for a boundary delivery in the local ENTSO-E test-set cache that resolves
`path`, copy it next to the case, and return its filename (empty when none
fits). Only uses what is already extracted locally — no download is triggered
from an upload.
"""
function _webui_try_supply_cgmes_boundary(path::AbstractString, directory::AbstractString)::String
  extracted = joinpath(CGMESImporter.cgmesTestSetCacheDir(), "extracted")
  isdir(extracted) || return ""
  candidates = String[]
  for (root, dirs, files) in walkdir(extracted)
    for d in dirs
      occursin(r"(?i)boundary|_bd_v|_bd$", d) && push!(candidates, joinpath(root, d))
    end
  end
  for cand in candidates
    ok = try
      !CGMESImporter.summarizeCGMES(path = [path, cand]).boundary_missing_hint
    catch
      false
    end
    ok || continue
    name = basename(cand) * ".zip"
    dest = joinpath(directory, name)
    ispath(dest) && return name
    try
      open(dest, "w") do io
        CGMESImporter.ZipArchives.ZipWriter(io) do w
          for f in sort(readdir(cand))
            endswith(lowercase(f), ".xml") || continue
            CGMESImporter.ZipArchives.zip_newfile(w, f)
            write(w, read(joinpath(cand, f)))
          end
        end
      end
      return name
    catch
      isfile(dest) && rm(dest; force = true)
      return ""
    end
  end
  return ""
end


_webui_is_runnable_dat_role(role::Symbol)::Bool = role in (:dtf_network_case, :dtf_network_case_with_outages)

"""
    _webui_is_user_selectable_case(name) -> Bool

Return whether `name` should be shown in the normal Web UI case selector.
Internal warm-up cases use the reserved `warmup_` prefix and remain available
to startup warm-up code by explicit path only. Generated Julia cache artifacts
are also hidden from the selector; users can still enter an explicit path in the
manual case field when they intentionally want to run such a file.
"""
# Memoization for the per-file content checks of the case selector: the ZIP
# boundary-set detection and the DAT role classification read file content,
# and the form re-scans the whole case directory on every render (measured
# ~10 s per render over a 59-file cache with large CGMES ZIPs). Results are
# keyed by (absolute path, mtime, size); a replaced or edited file is
# re-classified, deleted files age out via the key check on the next hit.
const _WEBUI_CASE_SCAN_LOCK = ReentrantLock()
const _WEBUI_CASE_SCAN_CACHE = Dict{String,Tuple{Float64,Int64,Any}}()

function _webui_file_scan_memo(f::Function, path::AbstractString)
  st = try
    stat(path)
  catch
    return f(path)
  end
  key = string(abspath(path), "|", nameof(f))
  mt = Float64(Base.Filesystem.mtime(st))
  sz = Int64(st.size)
  return lock(_WEBUI_CASE_SCAN_LOCK) do
    hit = get(_WEBUI_CASE_SCAN_CACHE, key, nothing)
    if hit !== nothing && hit[1] == mt && hit[2] == sz
      hit[3]
    else
      value = f(path)
      _WEBUI_CASE_SCAN_CACHE[key] = (mt, sz, value)
      value
    end
  end
end

_webui_classify_dat_content_cached(path::AbstractString)::Symbol = _webui_file_scan_memo(_webui_classify_dat_content, path)
_webui_is_cgmes_boundary_only_cached(path::AbstractString)::Bool = _webui_file_scan_memo(_webui_is_cgmes_boundary_only, path)

function _webui_is_user_selectable_case(name::AbstractString)::Bool
  lowered_name = lowercase(basename(name))
  _, extension = splitext(lowered_name)
  endswith(lowered_name, ".sparlectra-webui.yaml") && return false
  endswith(lowered_name, ".contingency-weights.csv") && return false
  startswith(lowered_name, "warmup_") && extension in (".m", ".jl") && return false
  extension == ".jl" && return false
  extension in (".m", ".dat", ".zip") || return false
  if extension == ".zip"
    # A boundary set is a companion delivery, not a runnable case — hide it
    # the same way FOR002 reference files are hidden.
    isfile(name) || return true
    return !_webui_is_cgmes_boundary_only_cached(name)
  end
  if extension == ".dat"
    isfile(name) || return !_is_for002_reference_dat(name)
    return _webui_is_runnable_dat_role(_webui_classify_dat_content_cached(name))
  end
  return true
end

function _webui_casefile_options_in_directory(directory::AbstractString)::Vector{String}
  isdir(directory) || return String[]
  files = filter(readdir(directory)) do name
    path = joinpath(directory, name)
    return isfile(path) && _webui_is_user_selectable_case(path)
  end
  return sort!(files; by = lowercase)
end

function _webui_for002_reference_options_in_directory(directory::AbstractString)::Vector{String}
  isdir(directory) || return String[]
  files = filter(readdir(directory)) do name
    path = joinpath(directory, name)
    return isfile(path) && lowercase(splitext(name)[2]) == ".dat" && _webui_classify_dat_content_cached(path) === :dtf_outage_or_reference
  end
  return sort!(files; by = lowercase)
end

function _webui_is_config_file(path::AbstractString)::Bool
  name = lowercase(basename(path))
  return endswith(name, ".yaml") || endswith(name, ".yml") || endswith(name, ".yaml.example") || endswith(name, ".yml.example")
end

"""
    _webui_config_file_options(application_root) -> Vector{String}

Return sorted YAML configuration files and YAML example templates from the Web
UI application's `examples` directory. A local `configuration.yaml` or
`configuration.yaml.example` is ordered first when present.
"""
function _webui_config_file_options(application_root::AbstractString)::Vector{String}
  directory = joinpath(application_root, "examples")
  isdir(directory) || return String[]
  files = filter(path -> isfile(path) && _webui_is_config_file(path), readdir(directory; join = true))
  priority(path) = lowercase(basename(path)) in ("configuration.yaml", "configuration.yaml.example") ? 0 : 1
  return sort!(normpath.(files); by = path -> (priority(path), lowercase(basename(path))))
end

function _webui_form_string(value)::String
  value === nothing && return ""
  value === missing && return ""
  value isa AbstractString && return String(value)
  value isa Symbol && return String(value)
  value isa Bool && return value ? "true" : "false"
  value isa Integer && return string(value)
  value isa AbstractFloat && return isfinite(value) ? string(value) : throw(ArgumentError("Web UI form value must be finite."))
  throw(ArgumentError("Unsupported Web UI form value type $(typeof(value))."))
end

function _webui_form_bool(value)::Bool
  value isa AbstractVector && return any(_webui_form_bool, value)
  value isa Bool && return value
  value === nothing && return false
  value === missing && return false
  value isa AbstractString && return lowercase(strip(value)) in ("1", "true", "yes", "on")
  value isa Integer && return value != 0
  throw(ArgumentError("Unsupported Web UI checkbox value type $(typeof(value))."))
end

function _webui_form_number_string(value)::String
  value isa Bool && throw(ArgumentError("Boolean is not a numeric Web UI form value."))
  value isa Integer && return string(value)
  value isa AbstractFloat && return isfinite(value) ? string(value) : throw(ArgumentError("Web UI numeric form value must be finite."))
  value isa AbstractString && begin
    text = strip(value)
    isempty(text) && throw(ArgumentError("Web UI numeric form value must not be empty."))
    parsed = tryparse(Float64, text)
    parsed === nothing && throw(ArgumentError("Invalid numeric Web UI form value $(repr(text))."))
    isfinite(parsed) || throw(ArgumentError("Web UI numeric form value must be finite."))
    return text
  end
  throw(ArgumentError("Unsupported Web UI numeric value type $(typeof(value))."))
end

function _webui_case_settings_filename(casefile::AbstractString)::String
  stem = splitext(basename(strip(String(casefile))))[1]
  isempty(stem) && throw(ArgumentError("Case-settings profile requires a case filename."))
  return string(stem, ".sparlectra-webui.yaml")
end

function _webui_normalized_case_key(casefile::AbstractString)::String
  stem = replace(lowercase(splitext(basename(strip(String(casefile))))[1]), r"[^a-z0-9]+" => "_")
  stem = strip(stem, '_')
  return isempty(stem) ? "case" : stem
end

function _webui_resolve_case_profile_source(casefile::AbstractString; case_directory::Union{Nothing,AbstractString} = nothing, application_root::AbstractString = _webui_application_root())::String
  raw = strip(String(casefile))
  isempty(raw) && return ""
  if isabspath(raw)
    return normpath(raw)
  end
  ".." in splitpath(raw) && return ""
  base = case_directory === nothing ? joinpath(application_root, "data", "mpower") : String(case_directory)
  return normpath(joinpath(base, raw))
end

function _webui_case_settings_path(output_root::AbstractString, casefile::AbstractString; case_directory::Union{Nothing,AbstractString} = nothing)::String
  source = isabspath(strip(String(casefile))) ? normpath(String(casefile)) : _webui_resolve_case_profile_source(casefile; case_directory)
  isempty(source) && throw(ArgumentError("Unsafe MATPOWER case path for case-settings profile."))
  return joinpath(dirname(source), _webui_case_settings_filename(source))
end

# Per-case N-1 contingency weight list (issue #331 Phase 5 follow-up): stored
# next to the case as `<stem>.contingency-weights.csv`, the exact two-column
# format `readContingencyWeightsCSV` parses. Mirrors the sidecar-YAML naming so
# it is human-editable and travels with the case (list-excluded, delete-cascaded).
function _webui_case_weights_filename(casefile::AbstractString)::String
  stem = splitext(basename(strip(String(casefile))))[1]
  isempty(stem) && throw(ArgumentError("Contingency weights file requires a case filename."))
  return string(stem, ".contingency-weights.csv")
end

# No output_root parameter (unlike the sidecar helper it mirrors): the file
# lives next to the case, output_root is irrelevant, so it is left out.
function _webui_case_weights_path(casefile::AbstractString; case_directory::Union{Nothing,AbstractString} = nothing)::String
  source = isabspath(strip(String(casefile))) ? normpath(String(casefile)) : _webui_resolve_case_profile_source(casefile; case_directory)
  isempty(source) && throw(ArgumentError("Unsafe MATPOWER case path for contingency weights."))
  return joinpath(dirname(source), _webui_case_weights_filename(source))
end

function _webui_log_case_settings_load(output_root::AbstractString, event::AbstractString; fields...)
  try
    record_webui_operation!(output_root, event; route = "/powerflow", method = "GET", user_action = true, fields...)
  catch
  end
  return nothing
end

function _webui_normalize_case_profile_form_value(field::AbstractString, value)
  value === nothing && return nothing
  value === missing && return nothing
  value isa AbstractVector && throw(ArgumentError("Case-settings field $(field) does not support vector values for form rendering."))
  type = get(_WEBUI_CASE_PROFILE_FIELD_TYPES, String(field), String)
  allowed = get(_WEBUI_CASE_PROFILE_SELECT_VALUES, String(field), nothing)
  if allowed !== nothing && value isa Bool
    !value && "off" in allowed && return "off"
    value && "on" in allowed && return "on"
  end
  normalized = if type === Bool
    _webui_form_bool(value)
  elseif type <: Integer
    value isa Bool && throw(ArgumentError("Case-settings field $(field) must be an integer."))
    if value isa Integer
      Int(value)
    elseif value isa AbstractFloat && isinteger(value) && isfinite(value)
      Int(value)
    else
      text = value isa AbstractString ? strip(value) : _webui_form_string(value)
      parsed = tryparse(type, text)
      parsed === nothing && throw(ArgumentError("Case-settings field $(field) has invalid integer value $(repr(text))."))
      parsed
    end
  elseif type <: AbstractFloat
    text = _webui_form_number_string(value)
    parsed = tryparse(type, text)
    parsed === nothing && throw(ArgumentError("Case-settings field $(field) has invalid numeric value $(repr(text))."))
    parsed
  else
    _webui_form_string(value)
  end
  if allowed !== nothing && !(_webui_form_string(normalized) in allowed)
    throw(ArgumentError("Case-settings field $(field) has unsupported value $(repr(normalized))."))
  end
  return normalized
end


function _webui_config_field_values(config_file::AbstractString)::Dict{String,Any}
  values = Dict{String,Any}()
  path = strip(String(config_file))
  isempty(path) && return values
  isfile(path) || return values
  try
    raw = load_yaml_dict(path)
    for (config_key, field, _) in _WEBUI_FORM_CONFIG_FIELDS
      value = _dotted_config_value(raw, config_key)
      value === nothing && continue
      values[field] = _webui_normalize_case_profile_form_value(field, value)
    end
    return values
  catch
    return Dict{String,Any}()
  end
end

function webui_form_state(; selected_casefile::AbstractString = "", selected_config_file::AbstractString = "", sidecar_profile = nothing, submitted_form = nothing)
  config_path = isempty(selected_config_file) ? DEFAULT_SPARLECTRA_CONFIG_PATH : selected_config_file
  values = Dict{String,Any}(spec.field => spec.default for spec in WEBUI_OPTION_SPECS)
  config_values = _webui_config_field_values(config_path)
  merge!(values, config_values)
  values["casefile"] = selected_casefile
  values["casefile_manual"] = ""
  values["config_file"] = config_path
  if sidecar_profile isa AbstractDict
    # Last edit wins: when the configuration file is NEWER than the saved
    # case settings, the config's own keys take precedence — editing the
    # YAML must show up on the next page load without deleting the sidecar
    # (fields the YAML does not set keep their saved case value).
    config_beats_sidecar = false
    deliberate_config_values = Dict{String,Any}()
    if haskey(sidecar_profile, "_profile_path")
      ppath = String(sidecar_profile["_profile_path"])
      config_beats_sidecar = isfile(config_path) && isfile(ppath) && mtime(config_path) > mtime(ppath)
      if config_beats_sidecar
        # Only DELIBERATE config values may outrank a saved case setup. A
        # value equal to the shipped template carries no user intent — the
        # configuration file is rewritten on startup/migration, and letting
        # its untouched defaults win would silently discard every stored
        # case setting (measured: case_SyntheticUSA lost its solver setup
        # and diverged after an unrelated config write).
        template_values = _webui_config_field_values(DEFAULT_SPARLECTRA_CONFIG_PATH)
        for (field, value) in config_values
          get(template_values, field, nothing) == value && continue
          deliberate_config_values[field] = value
        end
      end
    end
    applied_config_over_sidecar = false
    # Human-visible record of the replaced values: without it, a newer
    # configuration file silently flips saved case settings and the first
    # symptom is a diverged run (cgmes_realgrid lost start_values=auto to a
    # config carrying flat). Only fields whose values actually differ are
    # listed; entries are (config key, saved value, config value).
    config_over_sidecar_details = Vector{Tuple{String,String,String}}()
    for (field, value) in sidecar_profile
      field == "_profile_path" && continue
      haskey(_WEBUI_OPTION_BY_FIELD, String(field)) || continue
      if config_beats_sidecar && haskey(deliberate_config_values, String(field))
        applied_config_over_sidecar = true
        sidecar_value = _webui_normalize_case_profile_form_value(String(field), value)
        config_value = deliberate_config_values[String(field)]
        if sidecar_value != config_value
          spec = _WEBUI_OPTION_BY_FIELD[String(field)]
          display_key = spec.config_key === nothing ? String(field) : String(spec.config_key)
          push!(config_over_sidecar_details, (display_key, _webui_form_string(sidecar_value), _webui_form_string(config_value)))
        end
        continue
      end
      values[String(field)] = _webui_normalize_case_profile_form_value(String(field), value)
    end
    haskey(sidecar_profile, "_profile_path") && (values["_profile_path"] = sidecar_profile["_profile_path"])
    applied_config_over_sidecar && (values["_config_newer_than_profile"] = true)
    isempty(config_over_sidecar_details) || (values["_config_over_profile_details"] = sort!(config_over_sidecar_details))
  end
  if submitted_form isa AbstractDict
    for spec in WEBUI_OPTION_SPECS
      if spec.control == :checkbox
        values[spec.field] = _webui_form_bool(_webui_form_value(submitted_form, spec.field, false))
      elseif _webui_form_value(submitted_form, spec.field, nothing) !== nothing
        raw = _webui_form_value(submitted_form, spec.field)
        values[spec.field] = try
          _webui_normalize_case_profile_form_value(spec.field, raw)
        catch
          _webui_form_string(raw)
        end
      end
    end
    for field in ("casefile", "casefile_manual", "config_file")
      raw = _webui_form_value(submitted_form, field, nothing)
      raw === nothing || (values[field] = strip(_webui_form_string(raw)))
    end
  end
  return values
end

function _webui_load_case_settings(output_root::AbstractString, casefile::AbstractString; case_directory::Union{Nothing,AbstractString} = nothing)
  path = try
    _webui_case_settings_path(output_root, casefile; case_directory)
  catch err
    _webui_log_case_settings_load(output_root, "case_settings_load_failed"; casefile, status = "rejected", message = sprint(showerror, err))
    return nothing
  end
  if !isfile(path)
    _webui_log_case_settings_load(output_root, "case_settings_not_found"; casefile, profile_path = path, status = "missing")
    return nothing
  end
  try
    data = load_yaml_dict(path)
    if get(data, "schema_version", nothing) != 1
      _webui_log_case_settings_load(output_root, "case_settings_load_failed"; casefile, profile_path = path, status = "rejected", message = "unsupported schema_version")
      return nothing
    end
    if get(data, "profile_kind", "") != "webui_case_settings"
      _webui_log_case_settings_load(output_root, "case_settings_load_failed"; casefile, profile_path = path, status = "rejected", message = "unsupported profile_kind")
      return nothing
    end
    settings = get(data, "settings", nothing)
    if !(settings isa AbstractDict)
      _webui_log_case_settings_load(output_root, "case_settings_load_failed"; casefile, profile_path = path, status = "rejected", message = "settings must be a dictionary")
      return nothing
    end
    profile = Dict{String,Any}()
    for (key, value) in settings
      field = String(key)
      if !(field in _WEBUI_CASE_PROFILE_FIELDS)
        _webui_log_case_settings_load(output_root, "case_settings_field_ignored"; casefile, profile_path = path, status = "ignored", field, message = "unknown field")
        continue
      end
      try
        normalized = _webui_normalize_case_profile_form_value(field, value)
        normalized === nothing && continue
        profile[field] = normalized
      catch err
        _webui_log_case_settings_load(output_root, "case_settings_field_ignored"; casefile, profile_path = path, status = "ignored", field, message = sprint(showerror, err))
      end
    end
    profile["_profile_path"] = path
    _webui_log_case_settings_load(output_root, "case_settings_loaded"; casefile, profile_path = path, status = "loaded", setting_count = length(profile) - 1)
    return profile
  catch err
    _webui_log_case_settings_load(output_root, "case_settings_load_failed"; casefile, profile_path = path, status = "failed", message = sprint(showerror, err))
    return nothing
  end
end

function _webui_input_value(values::AbstractDict, field::AbstractString, default)::String
  return _webui_escape(_webui_form_string(get(values, field, default)))
end

function _webui_checked(values::AbstractDict, field::AbstractString, default::Bool)::String
  return _webui_parse_bool(get(values, field, default)) ? " checked" : ""
end

function _webui_selected(values::AbstractDict, field::AbstractString, default)
  return get(values, field, default)
end

function _webui_form_value(form::AbstractDict, key::String, default = nothing)
  haskey(form, key) && return form[key]
  symbol_key = Symbol(key)
  haskey(form, symbol_key) && return form[symbol_key]
  return default
end

_webui_is_dat_casefile(casefile::AbstractString)::Bool = lowercase(splitext(strip(String(casefile)))[2]) == ".dat"

# "Short circuit" button gating (the button is only
# selectable when the data is actually there): cheap byte scan over the
# delivery's XML contents for short-circuit source markers — synchronous
# machines, feeder short-circuit currents, or equivalent-injection
# impedances. In-memory via collectCGMESFiles (nested ZIPs included), cached
# per path+mtime so a 6209-bus delivery is scanned once, not on every form
# render. A scan error must NOT lock the button: the service run reports
# missing data with an explicit failure reason either way.
const _WEBUI_SC_DATA_CACHE = Dict{String,Tuple{Float64,Bool}}()
const _WEBUI_SC_DATA_CACHE_LOCK = ReentrantLock()

function _webui_case_has_short_circuit_data(casefile::AbstractString)::Bool
  path = strip(String(casefile))
  (isfile(path) || isdir(path)) || return true
  stamp = try
    mtime(path)
  catch
    0.0
  end
  return lock(_WEBUI_SC_DATA_CACHE_LOCK) do
    hit = get(_WEBUI_SC_DATA_CACHE, path, nothing)
    hit !== nothing && hit[1] == stamp && return hit[2]
    result = try
      files = CGMESImporter.collectCGMESFiles(path)
      any(occursin("SynchronousMachine", f.content) || occursin("maxInitialSymShCCurrent", f.content) || occursin("EquivalentInjection.x", f.content) for f in files)
    catch
      true
    end
    _WEBUI_SC_DATA_CACHE[path] = (stamp, result)
    result
  end
end

# CGMES deliveries arrive either as a ZIP (typical upload) or as an unpacked
# folder; both are recognised by the API's auto-detection, this only
# pre-selects the matching form option.
function _webui_is_cgmes_casefile(casefile::AbstractString)::Bool
  value = strip(String(casefile))
  isempty(value) && return false
  lowercase(splitext(value)[2]) == ".zip" && return true
  return isempty(splitext(value)[2]) && isdir(value)
end

function _webui_parse_bool(value)::Bool
  return _webui_form_bool(value)
end

function _webui_parse_form_value(value, ::Type{Bool}, field::String)
  return _webui_parse_bool(value)
end

function _webui_parse_form_value(value, type::Type{<:Number}, field::String)
  text = value === nothing ? "" : strip(string(value))
  isempty(text) && throw(ArgumentError("Web UI field $(field) must not be empty."))
  try
    return parse(type, text)
  catch
    throw(ArgumentError("Web UI field $(field) has invalid value $(repr(text))."))
  end
end

_webui_parse_form_value(value, ::Type{String}, field::String) = strip(String(something(value, "")))

"""
    powerflow_webui_request(form; default_output_root="results/powerflow_service")

Convert browser form values into the dictionary accepted by
[`start_powerflow_run`](@ref). Only keys from
[`GUI_EDITABLE_CONFIG_KEYS`](@ref) are emitted as configuration overrides.
"""
function powerflow_webui_request(form::AbstractDict; default_output_root::AbstractString = "results/powerflow_service")::Dict{String,Any}
  existing_casefile = strip(String(something(_webui_form_value(form, "casefile", ""), "")))
  manual_casefile = strip(String(something(_webui_form_value(form, "casefile_manual", ""), "")))
  casefile = isempty(manual_casefile) ? existing_casefile : manual_casefile
  isempty(casefile) && throw(ArgumentError("Select an existing MATPOWER case or type a case name."))
  config_file = strip(String(something(_webui_form_value(form, "config_file", ""), "")))
  output_root = String(default_output_root)
  ignore_webui_settings = if _webui_form_value(form, "ignore_webui_settings", nothing) !== nothing
    _webui_parse_form_value(_webui_form_value(form, "ignore_webui_settings", "false"), Bool, "ignore_webui_settings")
  else
    !_webui_parse_form_value(_webui_form_value(form, "apply_webui_runtime_overrides", "true"), Bool, "apply_webui_runtime_overrides")
  end
  apply_runtime_overrides = !ignore_webui_settings
  overrides = Dict{String,Any}()
  if apply_runtime_overrides
    for (config_key, field, type) in _WEBUI_FORM_CONFIG_FIELDS
      config_key in GUI_EDITABLE_CONFIG_KEYS || error("Web UI field $(field) is not GUI-editable.")
      spec = _webui_option_spec(field)
      _webui_form_value(form, field, nothing) === nothing && continue
      raw = _webui_form_value(form, field)
      overrides[config_key] = _webui_parse_form_value(raw, type, field)
    end
  end
  request_options = Dict{String,Any}()
  for field in _WEBUI_CASE_PROFILE_EXTRA_FIELDS
    spec = _webui_option_spec(field)
    raw = _webui_form_value(form, field, spec.default)
    request_options[field] = _webui_parse_form_value(raw, spec.value_type, field)
  end
  case_format = strip(String(something(_webui_form_value(form, "case_format", "auto"), "auto")))
  for002_reference_file = strip(String(something(_webui_form_value(form, "for002_reference_file", ""), "")))
  outage_mode = strip(String(something(_webui_form_value(form, "dtf_outage_selection_mode", "none"), "none")))
  raw_selection = _webui_form_value(form, "dtf_outage_selection", String[])
  selection = raw_selection isa AbstractVector ? String.(raw_selection) : (isempty(strip(String(raw_selection))) ? String[] : [String(raw_selection)])
  return Dict{String,Any}(
    "casefile" => casefile,
    "case_format" => case_format,
    "for002_reference_file" => isempty(for002_reference_file) ? nothing : for002_reference_file,
    "run_dtf_outages" => outage_mode != "none",
    "dtf_outage_selection" => selection,
    "dtf_outage_selection_mode" => outage_mode,
    "compare_for002_outages" => !isempty(for002_reference_file) && outage_mode != "none",
    "write_outage_artifacts" => _webui_parse_form_value(_webui_form_value(form, "write_outage_artifacts", "true"), Bool, "write_outage_artifacts"),
    "write_outage_matpower_exports" => _webui_parse_form_value(_webui_form_value(form, "write_outage_matpower_exports", "false"), Bool, "write_outage_matpower_exports"),
    "matpower_export_requested" => _webui_parse_form_value(_webui_form_value(form, "matpower_export_requested", "false"), Bool, "matpower_export_requested"),
    "config_file" => config_file,
    "output_root" => output_root,
    "config_overrides" => overrides,
    "config_override_source" => apply_runtime_overrides ? "webui_form_runtime" : "user_yaml",
    "performance_timing" => request_options["performance_timing"],
    # Normal "Start PowerFlow run" submissions never write diagnose.log; only
    # the dedicated "Diagnose" self-check does (diagnose_mode forces
    # run_diagnostics = true server-side in start_powerflow_run regardless of
    # this value). There is no "Run diagnostics" checkbox for a normal run.
    "run_diagnostics" => false,
    "diagnose_mode" => _webui_parse_form_value(_webui_form_value(form, "diagnose_mode", "false"), Bool, "diagnose_mode"),
    "short_circuit_mode" => _webui_parse_form_value(_webui_form_value(form, "short_circuit_mode", "false"), Bool, "short_circuit_mode"),
    "import_analysis_mode" => _webui_parse_form_value(_webui_form_value(form, "import_analysis_mode", "false"), Bool, "import_analysis_mode"),
    "contingency_mode" => _webui_parse_form_value(_webui_form_value(form, "contingency_mode", "false"), Bool, "contingency_mode"),
    # N-1 outage kind is a RUN parameter (branch / generator), read as a plain
    # request key like dtf_outage_selection, not a config override
    "contingency_kind" => strip(String(something(_webui_form_value(form, "contingency_kind", "branch"), "branch"))),
    "detailed_result_csv" => request_options["detailed_result_csv"],
    "detailed_result_csv_format" => request_options["detailed_result_csv_format"],
    "export_cgmes" => request_options["export_cgmes"],
  )
end

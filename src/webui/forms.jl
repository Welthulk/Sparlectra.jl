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

Return sorted user-selectable MATPOWER `.m` case names from the Web UI
application's `data/mpower` directory. Generated `.jl` cache artifacts are
internal and are not shown in the selector. Missing or empty directories produce
an empty list.
"""
function _webui_casefile_options(application_root::AbstractString)::Vector{String}
  return _webui_casefile_options_in_directory(joinpath(application_root, "data", "mpower"))
end

function _webui_casefile_options_in_directory(directory::AbstractString)::Vector{String}
  isdir(directory) || return String[]
  files = filter(readdir(directory)) do name
    extension = lowercase(splitext(name)[2])
    return isfile(joinpath(directory, name)) && extension == ".m"
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

function webui_form_state(; selected_casefile::AbstractString = "", selected_config_file::AbstractString = "", sidecar_profile = nothing, submitted_form = nothing)
  values = Dict{String,Any}(spec.field => spec.default for spec in WEBUI_OPTION_SPECS)
  values["casefile"] = selected_casefile
  values["casefile_manual"] = ""
  values["config_file"] = isempty(selected_config_file) ? DEFAULT_SPARLECTRA_CONFIG_PATH : selected_config_file
  if sidecar_profile isa AbstractDict
    for (field, value) in sidecar_profile
      field == "_profile_path" && continue
      haskey(_WEBUI_OPTION_BY_FIELD, String(field)) || continue
      values[String(field)] = _webui_normalize_case_profile_form_value(String(field), value)
    end
    haskey(sidecar_profile, "_profile_path") && (values["_profile_path"] = sidecar_profile["_profile_path"])
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
  overrides = Dict{String,Any}()
  for (config_key, field, type) in _WEBUI_FORM_CONFIG_FIELDS
    config_key in GUI_EDITABLE_CONFIG_KEYS || error("Web UI field $(field) is not GUI-editable.")
    spec = _webui_option_spec(field)
    _webui_form_value(form, field, nothing) === nothing && continue
    raw = _webui_form_value(form, field)
    overrides[config_key] = _webui_parse_form_value(raw, type, field)
  end
  request_options = Dict{String,Any}()
  for field in _WEBUI_CASE_PROFILE_EXTRA_FIELDS
    spec = _webui_option_spec(field)
    raw = _webui_form_value(form, field, spec.default)
    request_options[field] = _webui_parse_form_value(raw, spec.value_type, field)
  end
  return Dict{String,Any}(
    "casefile" => casefile,
    "config_file" => config_file,
    "output_root" => output_root,
    "config_overrides" => overrides,
    "performance_timing" => request_options["performance_timing"],
    "run_diagnostics" => request_options["run_diagnostics"],
    "detailed_result_csv" => request_options["detailed_result_csv"],
    "detailed_result_csv_format" => request_options["detailed_result_csv_format"],
  )
end

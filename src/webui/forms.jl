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

const _WEBUI_FORM_CONFIG_FIELDS = (
  ("power_flow.tol", "power_flow_tol", Float64),
  ("power_flow.max_iter", "power_flow_max_iter", Int),
  ("power_flow.autodamp", "power_flow_autodamp", Bool),
  ("power_flow.autodamp_min", "power_flow_autodamp_min", Float64),
  ("power_flow.qlimits.enabled", "power_flow_qlimits_enabled", Bool),
  ("power_flow.wrong_branch_detection", "power_flow_wrong_branch_detection", String),
  ("power_flow.start_mode.angle_mode", "power_flow_start_angle_mode", String),
  ("power_flow.start_mode.voltage_mode", "power_flow_start_voltage_mode", String),
  ("matpower_import.auto_profile", "matpower_import_auto_profile", String),
  ("matpower_import.ratio", "matpower_import_ratio", String),
  ("matpower_import.shift_sign", "matpower_import_shift_sign", Float64),
  ("matpower_import.shift_unit", "matpower_import_shift_unit", String),
  ("matpower_import.bus_shunt_model", "matpower_import_bus_shunt_model", String),
  ("matpower_import.pv_voltage_source", "matpower_import_pv_voltage_source", String),
  ("matpower_import.compare_voltage_reference", "matpower_import_compare_voltage_reference", String),
  ("output.logfile_results", "output_logfile_results", String),
  ("benchmark.enabled", "benchmark_enabled", Bool),
  ("benchmark.samples", "benchmark_samples", Int),
  ("benchmark.seconds", "benchmark_seconds", Float64),
)

const _WEBUI_PERFORMANCE_TIMING_VALUES = WEBUI_PERFORMANCE_TIMING_VALUES

function _webui_form_value(form::AbstractDict, key::String, default = nothing)
  haskey(form, key) && return form[key]
  symbol_key = Symbol(key)
  haskey(form, symbol_key) && return form[symbol_key]
  return default
end

function _webui_parse_bool(value)::Bool
  value isa Bool && return value
  value === nothing && return false
  return lowercase(strip(String(value))) in ("1", "true", "yes", "on")
end

function _webui_parse_form_value(value, ::Type{Bool}, field::String)
  return _webui_parse_bool(value)
end

function _webui_parse_form_value(value, type::Type{<:Number}, field::String)
  text = value === nothing ? "" : strip(String(value))
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
    overrides[config_key] = _webui_parse_form_value(_webui_form_value(form, field), type, field)
  end
  return Dict{String,Any}(
    "casefile" => casefile,
    "config_file" => config_file,
    "output_root" => output_root,
    "config_overrides" => overrides,
    "performance_timing" => _webui_parse_form_value(_webui_form_value(form, "performance_timing", "off"), String, "performance_timing"),
    "run_diagnostics" => _webui_parse_bool(_webui_form_value(form, "run_diagnostics")),
    "detailed_result_csv" => _webui_parse_bool(_webui_form_value(form, "detailed_result_csv")),
    "detailed_result_csv_format" => _webui_parse_form_value(_webui_form_value(form, "detailed_result_csv_format", "technical"), String, "detailed_result_csv_format"),
  )
end

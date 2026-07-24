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

"""Dotted configuration keys accepted from GUI/API override input."""
const GUI_EDITABLE_CONFIG_KEYS = Set([
  "power_flow.method",
  "power_flow.tol",
  "power_flow.max_iter",
  "power_flow.autodamp",
  "power_flow.autodamp_min",
  "power_flow.qlimits.enabled",
  "power_flow.qlimits.enforcement_mode",
  "power_flow.solver",
  "power_flow.apslf.order",
  "power_flow.apslf.use_pade",
  "power_flow.apslf.nr_polish",
  "power_flow.apslf_start.enabled",
  "power_flow.apslf_start.order",
  "power_flow.wrong_branch_detection",
  "power_flow.start_mode.angle_mode",
  "power_flow.start_mode.voltage_mode",
  "power_flow.start_mode.dc_seed_unconditional",
  "power_flow.start_current_iteration.enabled",
  "power_flow.start_current_iteration.max_iter",
  "power_flow.start_current_iteration.tol",
  "power_flow.start_current_iteration.damping",
  "power_flow.start_current_iteration.accept_only_if_improved",
  "power_flow.start_current_iteration.min_improvement_factor",
  "power_flow.start_current_iteration.vm_min_pu",
  "power_flow.start_current_iteration.vm_max_pu",
  "power_flow.start_current_iteration.max_angle_step_deg",
  "power_flow.start_current_iteration.only_for_large_cases",
  "power_flow.merit.enabled",
  "power_flow.merit.armijo_c1",
  "power_flow.merit.fallback_max_mismatch",
  "power_flow.trust_region.enabled",
  "power_flow.trust_region.initial_radius",
  "power_flow.trust_region.eta_accept",
  "power_flow.trust_region.step_mode",
  "power_flow.islands.enabled",
  "power_flow.islands.mode",
  "power_flow.islands.reference_policy",
  "power_flow.islands.diagnostic_continue_after_failure",
  "matpower_import.auto_profile",
  "matpower_import.ratio",
  "matpower_import.shift_sign",
  "matpower_import.shift_unit",
  "matpower_import.bus_shunt_model",
  "matpower_import.pv_voltage_source",
  "matpower_import.compare_voltage_reference",
  "matpower_import.apply_bus_names",
  "matpower_import.apply_branch_names",
  "matpower_import.apply_branch_kind",
  "matpower_import.import_for001_contingencies",
  "matpower_import.matpower_dcline_mode",
  "matpower_export.write_solution",
  "transformer.tap_changer_model",
  "output.logfile_results",
  "output.detailed_result_csv_write_mode",
  "output.detailed_result_csv_exporter",
  "output.detailed_result_csv_direct_threshold_buses",
  "output.detailed_result_csv_buffer_initial_bytes",
  "output.detailed_result_csv_buffer_max_bytes",
  "output.detailed_result_csv_streaming_threshold_rows",
  "benchmark.enabled",
  "benchmark.samples",
  "benchmark.seconds",
])

function _flatten_config_keys!(keys_out::Set{String}, raw::AbstractDict, prefix::String = "")
  for (key, value) in raw
    path = isempty(prefix) ? String(key) : string(prefix, ".", key)
    value isa AbstractDict ? _flatten_config_keys!(keys_out, value, path) : push!(keys_out, path)
  end
  return keys_out
end

function _validate_override_type(key::String, value, expected::Type)
  valid = expected === Bool ? value isa Bool : expected === Int ? value isa Integer && !(value isa Bool) : value isa Real && !(value isa Bool)
  valid || throw(ArgumentError("Override $(key) has invalid type $(typeof(value)); expected $(expected)."))
  return value
end

function _validate_gui_override_value(key::String, value)
  if key in ("power_flow.autodamp", "power_flow.qlimits.enabled", "power_flow.start_current_iteration.enabled", "power_flow.start_current_iteration.accept_only_if_improved", "power_flow.start_current_iteration.only_for_large_cases", "power_flow.merit.enabled", "power_flow.merit.fallback_max_mismatch", "power_flow.trust_region.enabled", "power_flow.apslf.use_pade", "power_flow.apslf.nr_polish", "power_flow.apslf_start.enabled", "power_flow.islands.enabled", "power_flow.islands.diagnostic_continue_after_failure", "benchmark.enabled", "matpower_import.apply_bus_names", "matpower_import.apply_branch_names", "matpower_import.apply_branch_kind", "matpower_import.import_for001_contingencies", "matpower_export.write_solution")
    _validate_override_type(key, value, Bool)
  elseif key in ("power_flow.max_iter", "power_flow.start_current_iteration.max_iter", "power_flow.apslf.order", "power_flow.apslf_start.order", "benchmark.samples", "output.detailed_result_csv_direct_threshold_buses", "output.detailed_result_csv_buffer_initial_bytes", "output.detailed_result_csv_buffer_max_bytes", "output.detailed_result_csv_streaming_threshold_rows")
    _validate_override_type(key, value, Int)
    if key == "power_flow.start_current_iteration.max_iter"
      value >= 0 || throw(ArgumentError("Override $(key) must be non-negative; got $(value)."))
    elseif key == "output.detailed_result_csv_buffer_initial_bytes"
      value >= 0 || throw(ArgumentError("Override $(key) must be non-negative; got $(value)."))
    else
      value > 0 || throw(ArgumentError("Override $(key) must be positive; got $(value)."))
    end
  elseif key in ("power_flow.tol", "power_flow.autodamp_min", "power_flow.start_current_iteration.tol", "power_flow.start_current_iteration.damping", "power_flow.start_current_iteration.min_improvement_factor", "power_flow.start_current_iteration.vm_min_pu", "power_flow.start_current_iteration.vm_max_pu", "power_flow.start_current_iteration.max_angle_step_deg", "power_flow.merit.armijo_c1", "power_flow.trust_region.initial_radius", "power_flow.trust_region.eta_accept", "benchmark.seconds", "matpower_import.shift_sign")
    _validate_override_type(key, value, Float64)
    if key == "matpower_import.shift_sign"
      isfinite(value) && value in (-1.0, 1.0) || throw(ArgumentError("Override $(key) must be -1.0 or 1.0; got $(value)."))
    else
      isfinite(value) && value > 0 || throw(ArgumentError("Override $(key) must be finite and positive; got $(value)."))
    end
    key == "power_flow.autodamp_min" && value > 1 && throw(ArgumentError("Override power_flow.autodamp_min must be <= 1; got $(value)."))
    key == "power_flow.start_current_iteration.damping" && value > 1 && throw(ArgumentError("Override $(key) must be <= 1; got $(value)."))
    key == "power_flow.merit.armijo_c1" && value >= 0.5 && throw(ArgumentError("Override power_flow.merit.armijo_c1 must be < 0.5; got $(value)."))
    key == "power_flow.trust_region.initial_radius" && value > 10.0 && throw(ArgumentError("Override power_flow.trust_region.initial_radius must be <= 10.0 (the default trust_region.max_radius, not GUI-editable); got $(value)."))
    key == "power_flow.trust_region.eta_accept" && value >= 1.0 && throw(ArgumentError("Override power_flow.trust_region.eta_accept must be < 1.0; got $(value)."))
  elseif key == "power_flow.method"
    method = _as_symbol_cfg(value)
    method === :rectangular || throw(ArgumentError(unsupported_powerflow_method_message(method)))
  elseif key == "power_flow.solver"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), POWERFLOW_SOLVER_VALUES)
  elseif key == "power_flow.qlimits.enforcement_mode"
    _canonical_qlimit_enforcement_mode(_as_symbol_cfg(value))
  elseif key == "power_flow.wrong_branch_detection"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), WRONG_BRANCH_DETECTION_VALUES)
  elseif key == "power_flow.start_mode.angle_mode"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), POWERFLOW_START_ANGLE_MODE_VALUES)
  elseif key == "power_flow.start_mode.voltage_mode"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), POWERFLOW_START_VOLTAGE_MODE_VALUES)
  elseif key == "power_flow.trust_region.step_mode"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), TRUST_REGION_STEP_MODE_VALUES)
  elseif key == "power_flow.islands.mode"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), POWERFLOW_ISLAND_MODE_VALUES)
  elseif key == "power_flow.islands.reference_policy"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), POWERFLOW_ISLAND_REFERENCE_POLICY_VALUES)
  elseif key == "matpower_import.auto_profile"
    _validate_allowed_symbol(key, _as_auto_profile_symbol_cfg(value), MATPOWER_AUTO_PROFILE_VALUES)
  elseif key == "matpower_import.ratio"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), MATPOWER_RATIO_VALUES)
  elseif key == "matpower_import.shift_unit"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), MATPOWER_SHIFT_UNIT_VALUES)
  elseif key == "matpower_import.bus_shunt_model"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), MATPOWER_BUS_SHUNT_MODEL_VALUES)
  elseif key == "matpower_import.pv_voltage_source"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), MATPOWER_PV_VOLTAGE_SOURCE_VALUES)
  elseif key == "matpower_import.compare_voltage_reference"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), MATPOWER_COMPARE_VOLTAGE_REFERENCE_VALUES)
  elseif key == "matpower_import.matpower_dcline_mode"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), MATPOWER_DCLINE_MODE_VALUES)
  elseif key == "transformer.tap_changer_model"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), TRANSFORMER_TAP_CHANGER_MODEL_VALUES)
  elseif key == "output.logfile_results"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), OUTPUT_LOGFILE_RESULTS_VALUES)
  elseif key == "output.detailed_result_csv_write_mode"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), OUTPUT_DETAILED_RESULT_CSV_WRITE_MODE_VALUES)
  elseif key == "output.detailed_result_csv_exporter"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), OUTPUT_DETAILED_RESULT_CSV_EXPORTER_VALUES)
  end
  return nothing
end

function _set_dotted_override!(overrides::Dict{String,Any}, key::String, value)
  parts = split(key, '.')
  current = overrides
  for part in parts[1:(end-1)]
    child = get!(current, part, Dict{String,Any}())
    child isa Dict{String,Any} || throw(ArgumentError("Override path $(key) conflicts with another override."))
    current = child
  end
  current[parts[end]] = value
  return overrides
end

"""
    validate_gui_config_overrides(config_overrides) -> Dict{String,Any}

Validate GUI-supplied dotted configuration keys and return the equivalent
nested configuration dictionary. Only [`GUI_EDITABLE_CONFIG_KEYS`](@ref) are
accepted. Invalid keys, types, enum values, and ranges throw `ArgumentError`
before power-flow execution.
"""
function validate_gui_config_overrides(config_overrides::AbstractDict)::Dict{String,Any}
  known_keys = _flatten_config_keys!(Set{String}(), load_yaml_dict(DEFAULT_SPARLECTRA_CONFIG_PATH))
  nested = Dict{String,Any}()
  for (raw_key, value) in config_overrides
    key = String(raw_key)
    key in known_keys || throw(ArgumentError("Unknown Sparlectra configuration override key: $(key)."))
    key in GUI_EDITABLE_CONFIG_KEYS || throw(ArgumentError("Sparlectra configuration key $(key) is not allowed for GUI editing."))
    _validate_gui_override_value(key, value)
    _set_dotted_override!(nested, key, value)
  end
  return nested
end

function _load_api_config(config_file::String, nested_overrides::Dict{String,Any})
  raw, _ = _load_and_validate_config(DEFAULT_SPARLECTRA_CONFIG_PATH, config_file; cli_overrides = Dict{String,Any}(), overrides = nested_overrides)
  return SparlectraConfig(raw), raw
end

const CONFIG_OVERRIDE_REPORT_KEYS = (
  "matpower_import.auto_profile",
  "matpower_import.compare_voltage_reference",
  "matpower_import.matpower_dcline_mode",
  "transformer.tap_changer_model",
  "power_flow.tol",
  "power_flow.max_iter",
  "power_flow.autodamp",
  "power_flow.autodamp_min",
  "power_flow.start_mode.angle_mode",
  "power_flow.start_mode.voltage_mode",
  "power_flow.qlimits.enabled",
  "power_flow.qlimits.enforcement_mode",
  "power_flow.start_current_iteration.enabled",
  "power_flow.merit.enabled",
  "power_flow.trust_region.enabled",
  "power_flow.islands.enabled",
  "power_flow.islands.mode",
  "power_flow.islands.reference_policy",
  "power_flow.islands.diagnostic_continue_after_failure",
)

function _dotted_config_value(raw::AbstractDict, key::AbstractString)
  current = raw
  for part in split(String(key), '.')
    current isa AbstractDict || return nothing
    haskey(current, part) || return nothing
    current = current[part]
  end
  return current
end

function _set_dotted_metadata!(raw::Dict{String,Any}, key::AbstractString, value)
  parts = split(String(key), '.')
  current = raw
  for part in parts[1:(end - 1)]
    child = get!(current, part, Dict{String,Any}())
    child isa Dict{String,Any} || (child = current[part] = Dict{String,Any}())
    current = child
  end
  current[parts[end]] = value
  return raw
end

function _config_source_report(config_file::String, nested_overrides::Dict{String,Any}, effective_raw::AbstractDict; override_source::AbstractString = "explicit_api_request", auto_profile_result = nothing)
  user = isfile(config_file) ? load_yaml_dict(config_file) : Dict{String,Any}()
  report = Dict{String,Any}()
  for key in CONFIG_OVERRIDE_REPORT_KEYS
    source = _dotted_config_value(nested_overrides, key) !== nothing ? String(override_source) :
             _dotted_config_value(user, key) !== nothing ? "user_yaml" :
             "default"
    _set_dotted_metadata!(report, key, Dict{String,Any}(
      "value" => _dotted_config_value(effective_raw, key),
      "source" => source,
      "precedence" => "default < user_yaml < case_sidecar < webui_form_runtime < explicit_api_request",
    ))
  end
  if auto_profile_result !== nothing
    for pair in auto_profile_result.applied
      option = string("matpower_import.", first(pair))
      option in CONFIG_OVERRIDE_REPORT_KEYS || continue
      entry = _dotted_config_value(report, option)
      entry isa AbstractDict || continue
      entry["requested_value"] = _auto_profile_row_value(auto_profile_result.rows, option, :current)
      entry["value_before_auto_profile"] = _auto_profile_row_value(auto_profile_result.rows, option, :current)
      entry["applied_value"] = string(last(pair))
      entry["value"] = _dotted_config_value(effective_raw, option)
      entry["source"] = "matpower_auto_profile_apply"
      entry["reason"] = _auto_profile_row_value(auto_profile_result.rows, option, :reason)
      entry["replaced_explicit_or_runtime_value"] = _dotted_config_value(nested_overrides, option) !== nothing || _dotted_config_value(user, option) !== nothing
    end
  end
  return report
end

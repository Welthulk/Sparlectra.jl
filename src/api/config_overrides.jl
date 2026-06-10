"""Dotted configuration keys accepted from GUI/API override input."""
const GUI_EDITABLE_CONFIG_KEYS = Set([
  "power_flow.method",
  "power_flow.tol",
  "power_flow.max_iter",
  "power_flow.autodamp",
  "power_flow.autodamp_min",
  "power_flow.qlimits.enabled",
  "power_flow.wrong_branch_detection",
  "power_flow.start_mode.angle_mode",
  "power_flow.start_mode.voltage_mode",
  "output.logfile_results",
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
  if key in ("power_flow.autodamp", "power_flow.qlimits.enabled", "benchmark.enabled")
    _validate_override_type(key, value, Bool)
  elseif key in ("power_flow.max_iter", "benchmark.samples")
    _validate_override_type(key, value, Int)
    value > 0 || throw(ArgumentError("Override $(key) must be positive; got $(value)."))
  elseif key in ("power_flow.tol", "power_flow.autodamp_min", "benchmark.seconds")
    _validate_override_type(key, value, Float64)
    isfinite(value) && value > 0 || throw(ArgumentError("Override $(key) must be finite and positive; got $(value)."))
    key == "power_flow.autodamp_min" && value > 1 && throw(ArgumentError("Override power_flow.autodamp_min must be <= 1; got $(value)."))
  elseif key == "power_flow.method"
    method = _as_symbol_cfg(value)
    method === :rectangular || throw(ArgumentError(unsupported_powerflow_method_message(method)))
  elseif key == "power_flow.wrong_branch_detection"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), WRONG_BRANCH_DETECTION_VALUES)
  elseif key == "power_flow.start_mode.angle_mode"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), POWERFLOW_START_ANGLE_MODE_VALUES)
  elseif key == "power_flow.start_mode.voltage_mode"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), POWERFLOW_START_VOLTAGE_MODE_VALUES)
  elseif key == "output.logfile_results"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), OUTPUT_LOGFILE_RESULTS_VALUES)
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

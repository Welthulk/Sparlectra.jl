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

# file: src/config/config_overrides.jl
# purpose: GUI/API configuration override handling: the editable-key
#          allowlist, override validation (validate_gui_config_overrides),
#          and merging overrides into the loaded API configuration
"""Dotted configuration keys accepted from GUI/API override input."""
const GUI_EDITABLE_CONFIG_KEYS = Set([
  "power_flow.method",
  "power_flow.mode",
  "power_flow.tol",
  "power_flow.tol_MW",
  "power_flow.max_iter",
  "power_flow.autodamp",
  "power_flow.autodamp_min",
  "power_flow.qlimits.enabled",
  "power_flow.qlimits.enforcement_mode",
  "power_flow.solver",
  "power_flow.linear_solver",
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
  "power_flow.distributed_slack.enabled",
  "power_flow.distributed_slack.p_mode",
  "power_flow.external_grid.enabled",
  "power_flow.external_grid.source",
  "power_flow.external_grid.sk_MVA",
  "power_flow.external_grid.rx",
  "power_flow.rescue",
  "power_flow.dc.fallback",
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
  "cgmes_import.start_values",
  "cgmes_import.require_boundary",
  "cgmes_import.infer_base_voltages",
  "model.auto_profile",
  "matpower_import.ratio",
  "matpower_import.shift_sign",
  "matpower_import.shift_unit",
  "model.bus_shunt_model",
  "matpower_import.pv_voltage_source",
  "matpower_import.compare_voltage_reference",
  "matpower_import.apply_bus_names",
  "matpower_import.apply_branch_names",
  "matpower_import.apply_branch_kind",
  "matpower_import.import_for001_contingencies",
  "matpower_import.matpower_dcline_mode",
  "model.net_cache_enabled",
  "cgmes_import.hvdc_mode",
  "matpower_export.write_solution",
  "model.tap_changer_model",
  "short_circuit.sweep_method",
  "output.logfile_results",
  "output.console_summary",
  "output.console_auto_profile",
  "output.console_diagnostics",
  "output.console_q_limit_events",
  "output.console_max_rows",
  "output.result_table_max_rows",
  "output.result_table_large_case_threshold_buses",
  "output.result_table_large_case_mode",
  "output.logfile_diagnostics",
  "output.logfile_performance",
  "output.logfile_warnings",
  "output.startup_latency_hint",
  "output.console_live",
  "output.detailed_result_csv_write_mode",
  "output.detailed_result_csv_exporter",
  "output.detailed_result_csv_direct_threshold_buses",
  "output.detailed_result_csv_buffer_initial_bytes",
  "output.detailed_result_csv_buffer_max_bytes",
  "output.detailed_result_csv_streaming_threshold_rows",
  "benchmark.enabled",
  "benchmark.samples",
  "benchmark.seconds",
  "runtime.parallel.enabled",
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
  if key in ("power_flow.autodamp", "power_flow.qlimits.enabled", "power_flow.start_current_iteration.enabled", "power_flow.start_current_iteration.accept_only_if_improved", "power_flow.start_current_iteration.only_for_large_cases", "power_flow.merit.enabled", "power_flow.merit.fallback_max_mismatch", "power_flow.trust_region.enabled", "power_flow.apslf.use_pade", "power_flow.apslf.nr_polish", "power_flow.apslf_start.enabled", "power_flow.islands.enabled", "power_flow.islands.diagnostic_continue_after_failure", "power_flow.rescue", "power_flow.dc.fallback", "cgmes_import.require_boundary", "cgmes_import.infer_base_voltages", "benchmark.enabled", "matpower_import.apply_bus_names", "matpower_import.apply_branch_names", "matpower_import.apply_branch_kind", "matpower_import.import_for001_contingencies", "model.net_cache_enabled", "matpower_export.write_solution", "output.console_live", "output.console_summary", "output.startup_latency_hint")
    _validate_override_type(key, value, Bool)
  elseif key in ("power_flow.max_iter", "power_flow.start_current_iteration.max_iter", "power_flow.apslf.order", "power_flow.apslf_start.order", "benchmark.samples", "output.detailed_result_csv_direct_threshold_buses", "output.detailed_result_csv_buffer_initial_bytes", "output.detailed_result_csv_buffer_max_bytes", "output.detailed_result_csv_streaming_threshold_rows", "output.console_max_rows", "output.result_table_max_rows", "output.result_table_large_case_threshold_buses")
    _validate_override_type(key, value, Int)
    if key == "power_flow.start_current_iteration.max_iter"
      value >= 0 || throw(ArgumentError("Override $(key) must be non-negative; got $(value)."))
    elseif key == "output.detailed_result_csv_buffer_initial_bytes"
      value >= 0 || throw(ArgumentError("Override $(key) must be non-negative; got $(value)."))
    else
      value > 0 || throw(ArgumentError("Override $(key) must be positive; got $(value)."))
    end
  elseif key in ("power_flow.tol", "power_flow.tol_MW", "power_flow.autodamp_min", "power_flow.start_current_iteration.tol", "power_flow.start_current_iteration.damping", "power_flow.start_current_iteration.min_improvement_factor", "power_flow.start_current_iteration.vm_min_pu", "power_flow.start_current_iteration.vm_max_pu", "power_flow.start_current_iteration.max_angle_step_deg", "power_flow.merit.armijo_c1", "power_flow.trust_region.initial_radius", "power_flow.trust_region.eta_accept", "benchmark.seconds", "matpower_import.shift_sign")
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
  elseif key == "power_flow.mode"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), POWERFLOW_MODE_VALUES)
  elseif key == "power_flow.method"
    method = _as_symbol_cfg(value)
    method === :rectangular || throw(ArgumentError(unsupported_powerflow_method_message(method)))
  elseif key == "power_flow.solver"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), POWERFLOW_SOLVER_VALUES)
  elseif key == "power_flow.linear_solver"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), POWERFLOW_LINEAR_SOLVER_VALUES)
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
  elseif key == "model.auto_profile"
    _validate_allowed_symbol(key, _as_auto_profile_symbol_cfg(value), MATPOWER_AUTO_PROFILE_VALUES)
  elseif key == "matpower_import.ratio"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), MATPOWER_RATIO_VALUES)
  elseif key == "matpower_import.shift_unit"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), MATPOWER_SHIFT_UNIT_VALUES)
  elseif key == "model.bus_shunt_model"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), MATPOWER_BUS_SHUNT_MODEL_VALUES)
  elseif key == "matpower_import.pv_voltage_source"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), MATPOWER_PV_VOLTAGE_SOURCE_VALUES)
  elseif key == "matpower_import.compare_voltage_reference"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), MATPOWER_COMPARE_VOLTAGE_REFERENCE_VALUES)
  elseif key == "matpower_import.matpower_dcline_mode"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), MATPOWER_DCLINE_MODE_VALUES)
  elseif key == "cgmes_import.hvdc_mode"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), CGMES_HVDC_MODE_VALUES)
  elseif key == "short_circuit.sweep_method"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), SHORT_CIRCUIT_SWEEP_METHOD_VALUES)
  elseif key == "model.tap_changer_model"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), TRANSFORMER_TAP_CHANGER_MODEL_VALUES)
  elseif key == "output.logfile_results"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), OUTPUT_LOGFILE_RESULTS_VALUES)
  elseif key == "output.console_auto_profile"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), OUTPUT_CONSOLE_AUTO_PROFILE_VALUES)
  elseif key == "output.console_diagnostics"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), OUTPUT_CONSOLE_DIAGNOSTICS_VALUES)
  elseif key == "output.console_q_limit_events"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), OUTPUT_CONSOLE_Q_LIMIT_EVENTS_VALUES)
  elseif key == "output.result_table_large_case_mode"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), OUTPUT_RESULT_TABLE_LARGE_CASE_MODE_VALUES)
  elseif key == "output.logfile_diagnostics"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), OUTPUT_LOGFILE_DIAGNOSTICS_VALUES)
  elseif key == "output.logfile_performance"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), OUTPUT_LOGFILE_PERFORMANCE_VALUES)
  elseif key == "output.logfile_warnings"
    _validate_allowed_symbol(key, _as_symbol_cfg(value), OUTPUT_LOGFILE_WARNINGS_VALUES)
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

# AbstractString, not String: the service entry points take AbstractString,
# and a path that arrives as a SubString (a split request field) must not
# fail here with a MethodError.
function _load_api_config(config_file::AbstractString, nested_overrides::Dict{String,Any}; case_scope_from_defaults::Bool = false)
  user_set = Set{String}()
  raw, _ = _load_and_validate_config(DEFAULT_SPARLECTRA_CONFIG_PATH, String(config_file); cli_overrides = Dict{String,Any}(), overrides = nested_overrides, user_set_out = user_set, case_scope_from_defaults)
  return _copy_sparlectra_with_user_keys(SparlectraConfig(raw), user_set), raw
end

"""
    case_config_path(case_path) -> String

The case configuration file that belongs to one case: `<stem>.config.yaml`
next to the case file. For the canonical `.scf.json` double extension the
stem strips both parts, so `case57.scf.json` binds `case57.config.yaml`.
"""
function case_config_path(case_path::AbstractString)::String
  dir = dirname(abspath(String(case_path)))
  stem, ext = splitext(basename(String(case_path)))
  if lowercase(ext) == ".json" && endswith(lowercase(stem), ".scf")
    stem = stem[1:end-4]
  end
  return joinpath(dir, string(stem, ".config.yaml"))
end

"""
    load_case_config(case_path) -> Dict{String,Any}

The flat dotted case-scope configuration of one case, read from its
[`case_config_path`](@ref) file; empty when no such file exists. The file
must open with the header keys `config_version`, `scope: case`, and
`case: <case file basename>`; a `case` value naming a different file is the
hard error `case_config_mismatch`, so a copied config cannot silently steer
the wrong case. Keys outside the case scope are refused with the same
wording as the in-file `sparlectra.config` check.
"""
function load_case_config(case_path::AbstractString)::Dict{String,Any}
  path = case_config_path(case_path)
  isfile(path) || return Dict{String,Any}()
  raw = load_yaml_dict(path)
  version = _config_file_version(raw, path)
  version < CONFIG_VERSION_CURRENT && _apply_config_aliases!(raw, version, path)
  scope = lowercase(strip(string(get(raw, "scope", ""))))
  scope == "case" || throw(ArgumentError("Case configuration $(path) must declare scope: case (found $(repr(scope)))."))
  declared_case = strip(string(get(raw, "case", "")))
  expected_case = basename(String(case_path))
  declared_case == expected_case || throw(ArgumentError("case_config_mismatch: $(path) declares case $(repr(declared_case)) but is loaded for $(repr(expected_case)). The case header must name the case file it configures."))
  for header_key in ("config_version", "scope", "case")
    delete!(raw, header_key)
  end
  # the Web UI persists its per-case form defaults under `form`; those are
  # request fields, not configuration keys, and are read separately
  delete!(raw, "form")
  out = _flatten_config_values!(Dict{String,Any}(), raw)
  bad = sort!(String[k for k in keys(out) if !scf_is_case_config_key(k)])
  isempty(bad) || throw(ArgumentError("Case configuration $(path): config key(s) $(join(bad, ", ")) are not case scope (logging, benchmarking, parallelism, Web UI and export settings belong in the configuration file). Remove them from the case configuration."))
  return out
end

"""
    write_case_config(case_file, config) -> String

Write the case-scope keys of `config` (flat dotted keys) as the case
configuration file of `case_file` ([`case_config_path`](@ref)), with the D8
header (`config_version`, `scope: case`, `case:`). Keys outside the case
scope are dropped with one warning naming them; an empty case scope removes
an existing file instead of leaving a stale one. Returns the file path.
"""
function write_case_config(case_file::AbstractString, config::AbstractDict; form::AbstractDict = Dict{String,Any}())::String
  path = case_config_path(case_file)
  keep = Dict{String,Any}(String(k) => v for (k, v) in config if scf_is_case_config_key(String(k)))
  dropped = sort!(String[String(k) for k in keys(config) if !scf_is_case_config_key(String(k))])
  isempty(dropped) || @warn "write_case_config: $(length(dropped)) configuration key(s) are not case scope and were not written (they belong in the configuration file)" keys = dropped
  if isempty(keep) && isempty(form)
    isfile(path) && rm(path)
    return path
  end
  tree = Dict{String,Any}()
  for (k, v) in keep
    _dotted_config_set!(tree, k, v)
  end
  # `form` carries Web UI request defaults for this case (SE and generator
  # options); they are form fields, not configuration keys, and never enter
  # the configuration resolution
  isempty(form) || (tree["form"] = Dict{String,Any}(String(k) => v for (k, v) in form))
  open(path, "w") do io
    println(io, "config_version: ", CONFIG_VERSION_CURRENT)
    println(io, "scope: case")
    println(io, "case: ", basename(String(case_file)))
    _write_yaml_dict(io, tree)
  end
  return path
end

"""
    ConfigResolveError

Wraps a failure of one resolution level of [`resolve_config`](@ref) with the
stable service failure reason for that level (`invalid_case_file` for the
deprecated in-case block, `invalid_case_config` for the case configuration
file, `invalid_config_override` for override validation,
`invalid_configuration` for the general file), so every service maps
failures identically without parsing messages.
"""
struct ConfigResolveError <: Exception
  reason::String
  err::Exception
end

Base.showerror(io::IO, e::ConfigResolveError) = showerror(io, e.err)

_config_resolve_reason(err) = err isa ConfigResolveError ? err.reason : "invalid_configuration"

"""
    resolve_config(config_file, case_path, overrides) -> NamedTuple

The one configuration precedence of the run path (design decision D5),
highest first: explicit API/CLI `overrides`, the case configuration file,
`sparlectra.config` inside an SCF case (deprecated level, one warning naming
the case configuration file as the new place), the general configuration
file, packaged defaults. A key not set on one level falls through to the
next; no mtime logic anywhere.

Exception, whenever a case configuration FILE exists (issue #1 point 1,
decided for `defaults`; format independent per review point 2): the
case is then self-contained, so CASE-scope keys skip the general-file
level and fall through from the case levels directly to the packaged
defaults, and the same case-plus-config pair computes the same numbers
on every installation. Machine-scope keys (output, benchmark, runtime,
webui, matpower_export) still come from the general file. A case
without a config file, including a bare `.scf.json`, keeps the full
chain, and so does the deprecated in-file block alone (legacy level,
legacy semantics).

Returns the effective typed `config`, the `effective_raw` dictionary for
the effective-config artifact, the merged flat and nested override chains,
and the per-level dictionaries for provenance. Level errors propagate; the
caller owns the failure mapping.
"""
function resolve_config(config_file::AbstractString, case_path::AbstractString, overrides::AbstractDict = Dict{String,Any}(); auto_profile_overrides::AbstractDict = Dict{String,Any}())
  scf_level = Dict{String,Any}()
  is_scf_case = lowercase(splitext(String(case_path))[2]) == ".json" && isfile(case_path)
  if is_scf_case
    scf_level = try
      scf_case_config(String(case_path))
    catch err
      throw(ConfigResolveError("invalid_case_file", err))
    end
    # maxlog per case file: repeating the same deprecation for every run of
    # the same case says nothing new and buried the rest of the output
    isempty(scf_level) || @warn "sparlectra.config inside $(basename(String(case_path))) is deprecated; move these settings to $(basename(case_config_path(case_path))) next to the case file. The block still applies, directly below that file in precedence." maxlog = 1 _id = Symbol("scf_cfg_", basename(String(case_path)))
  end
  case_level = try
    load_case_config(case_path)
  catch err
    throw(ConfigResolveError("invalid_case_config", err))
  end
  merged = merge(scf_level, case_level, Dict{String,Any}(String(k) => v for (k, v) in overrides))
  nested = try
    validate_gui_config_overrides(merged)
  catch err
    throw(ConfigResolveError("invalid_config_override", err))
  end
  config, effective_raw = _load_api_config(String(config_file), nested; case_scope_from_defaults = isfile(case_config_path(case_path)))
  # D11: auto-profile recommendations are their own level, the weakest
  # user-facing one; a recommendation applies only where neither the user
  # YAML (user_set_keys of the first pass) nor a case level nor an
  # explicit override set the key
  auto_level = Dict{String,Any}()
  if !isempty(auto_profile_overrides)
    for (k, v) in auto_profile_overrides
      key = String(k)
      (haskey(scf_level, key) || haskey(case_level, key) || haskey(merged, key)) && continue
      if key in config.user_set_keys
        # same rule as apply_auto_profile_level: presence in a fully
        # written YAML is not a choice, a NON-DEFAULT value is
        field = findfirst(==(key), _AUTO_PROFILE_FIELD_KEYS)
        field !== nothing && _auto_profile_current_value(config, field) != _auto_profile_template_default(field) && continue
      end
      auto_level[key] = v
    end
    if !isempty(auto_level)
      merged_auto = merge(auto_level, merged)
      nested = validate_gui_config_overrides(merged_auto)
      config, effective_raw = _load_api_config(String(config_file), nested; case_scope_from_defaults = isfile(case_config_path(case_path)))
    end
  end
  return (config = config, effective_raw = effective_raw, merged_overrides = merged, nested_overrides = nested, scf_config = scf_level, case_config = case_level, auto_profile_config = auto_level)
end

# D11 (stage 5 straggler): the auto-profile recommendations form their own
# precedence level, the WEAKEST user-facing one. The mapping names the
# dotted key of every field the MATPOWER auto profile may apply; the level
# helper enforces that an explicitly set key (YAML, case file, sidecar, or
# override; the typed config carries them in user_set_keys) always wins.
const _AUTO_PROFILE_FIELD_KEYS = Dict{Symbol,String}(
  :shift_sign => "matpower_import.shift_sign",
  :shift_unit => "matpower_import.shift_unit",
  :ratio => "matpower_import.ratio",
  :bus_shunt_model => "model.bus_shunt_model",
  :compare_voltage_reference => "matpower_import.compare_voltage_reference",
)

"""
    apply_auto_profile_level(cfg, applied_pairs) -> (config, applied, skipped)

Applies auto-profile recommendations as their own precedence level (D11):
each pair lands only when the user did not set its key explicitly;
explicitly set keys are returned in `skipped` so the caller can report
the yield. The copy helpers resolve at call time (they live with the
MATPOWER engine), so this stays the single application site without a
definition-time dependency on the adapters block.
"""
# the five level keys with their template defaults, read once from a
# defaults-only config: a refresh-written YAML carries EVERY key, so
# "explicitly set" cannot mean mere presence; the level yields when the
# user's value actually DIFFERS from the template default (a key set to
# its default is not a choice the auto profile must respect), and always
# when the key arrives through a case level or an explicit override
# (resolve_config handles those levels)
function _auto_profile_current_value(cfg::SparlectraConfig, field::Symbol)
  field in fieldnames(MatpowerImportConfig) && return getfield(cfg.matpower, field)
  return getfield(cfg.model, field)
end
const _AUTO_PROFILE_TEMPLATE_DEFAULTS = Dict{Symbol,Any}()
function _auto_profile_template_default(field::Symbol)
  isempty(_AUTO_PROFILE_TEMPLATE_DEFAULTS) && for f in keys(_AUTO_PROFILE_FIELD_KEYS)
    _AUTO_PROFILE_TEMPLATE_DEFAULTS[f] = _auto_profile_current_value(SparlectraConfig(Dict{String,Any}()), f)
  end
  return _AUTO_PROFILE_TEMPLATE_DEFAULTS[field]
end

function apply_auto_profile_level(cfg::SparlectraConfig, applied_pairs)
  applied = Pair{Symbol,Any}[]
  skipped = Pair{Symbol,Any}[]
  for pair in applied_pairs
    field = first(pair)
    key = get(_AUTO_PROFILE_FIELD_KEYS, field, nothing)
    key === nothing && throw(ArgumentError("auto profile applied unknown field $(field); add it to _AUTO_PROFILE_FIELD_KEYS."))
    explicit = key in cfg.user_set_keys && _auto_profile_current_value(cfg, field) != _auto_profile_template_default(field)
    explicit ? push!(skipped, pair) : push!(applied, pair)
  end
  isempty(applied) && return (config = cfg, applied = applied, skipped = skipped)
  mat_updates = Pair{Symbol,Any}[pair for pair in applied if first(pair) in fieldnames(MatpowerImportConfig)]
  model_updates = Pair{Symbol,Any}[pair for pair in applied if first(pair) in fieldnames(ModelConfig)]
  mat2 = isempty(mat_updates) ? cfg.matpower : _copy_matpower_with(cfg.matpower; mat_updates...)
  model2 = isempty(model_updates) ? cfg.model : _copy_model_with(cfg.model; model_updates...)
  return (config = _copy_config_with(cfg; matpower = mat2, model = model2), applied = applied, skipped = skipped)
end

const CONFIG_OVERRIDE_REPORT_KEYS = (
  "model.auto_profile",
  "matpower_import.compare_voltage_reference",
  "matpower_import.matpower_dcline_mode",
  "model.tap_changer_model",
  "power_flow.tol",
  "power_flow.tol_MW",
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
  "power_flow.rescue",
  "power_flow.dc.fallback",
  "cgmes_import.require_boundary",
  "cgmes_import.infer_base_voltages",
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

function _config_source_report(config_file::String, nested_overrides::Dict{String,Any}, effective_raw::AbstractDict; override_source::AbstractString = "explicit_api_request", auto_profile_result = nothing, explicit_overrides::AbstractDict = Dict{String,Any}(), case_config::AbstractDict = Dict{String,Any}(), scf_config::AbstractDict = Dict{String,Any}())
  user = isfile(config_file) ? load_yaml_dict(config_file) : Dict{String,Any}()
  report = Dict{String,Any}()
  # resolve_config merges the case levels INTO the override structure
  # (functionally right), so the report needs the per-level dictionaries to
  # label a value's true origin: before this, a case-sidecar value was
  # reported under the caller's override label and the case_sidecar source
  # named in the precedence line could never appear (found by the demo-case
  # pinning test, 2026-09-03). Callers that pass no levels keep the old
  # merged behavior.
  explicit_flat = Dict{String,Any}(String(k) => v for (k, v) in explicit_overrides)
  for key in CONFIG_OVERRIDE_REPORT_KEYS
    source = if haskey(explicit_flat, key)
      String(override_source)
    elseif haskey(case_config, key)
      "case_sidecar"
    elseif haskey(scf_config, key)
      "case_file_block"
    elseif _dotted_config_value(nested_overrides, key) !== nothing
      # merged fallback for callers without per-level dictionaries
      String(override_source)
    elseif _dotted_config_value(user, key) !== nothing
      "user_yaml"
    else
      "default"
    end
    _set_dotted_metadata!(report, key, Dict{String,Any}(
      "value" => _dotted_config_value(effective_raw, key),
      "source" => source,
      "precedence" => "default < user_yaml < case_file_block < case_sidecar < webui_form_runtime < explicit_api_request",
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

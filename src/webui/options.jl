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

"""Central metadata for one Web UI run-option form field."""
struct WebUIOptionSpec
  config_key::Union{Nothing,String}
  field::String
  value_type::Type
  control::Symbol
  default::Any
  allowed_values::Any
  section::Symbol
  save_in_case_sidecar::Bool
end

const _WEBUI_QLIMIT_ENFORCEMENT_MODE_VALUES = (:active_set, :classic_simultaneous, :classic_one_at_a_time)

const _WEBUI_PERFORMANCE_TIMING_VALUES = WEBUI_PERFORMANCE_TIMING_VALUES

const WEBUI_OPTION_SPECS = (
  WebUIOptionSpec("power_flow.tol", "power_flow_tol", Float64, :number, "1e-8", (), :basic, true),
  WebUIOptionSpec("power_flow.max_iter", "power_flow_max_iter", Int, :number, 80, (), :basic, true),
  WebUIOptionSpec("power_flow.autodamp", "power_flow_autodamp", Bool, :checkbox, true, (), :basic, true),
  WebUIOptionSpec("power_flow.autodamp_min", "power_flow_autodamp_min", Float64, :number, 0.05, (), :basic, true),
  WebUIOptionSpec("power_flow.qlimits.enabled", "power_flow_qlimits_enabled", Bool, :checkbox, true, (), :basic, true),
  WebUIOptionSpec("power_flow.qlimits.enforcement_mode", "power_flow_qlimits_enforcement_mode", String, :select, "active_set", _WEBUI_QLIMIT_ENFORCEMENT_MODE_VALUES, :basic, true),
  WebUIOptionSpec("power_flow.solver", "power_flow_solver", String, :select, "rectangular", POWERFLOW_SOLVER_VALUES, :basic, true),
  WebUIOptionSpec("power_flow.apslf.order", "power_flow_apslf_order", Int, :number, 40, (), :basic, true),
  WebUIOptionSpec("power_flow.apslf.use_pade", "power_flow_apslf_use_pade", Bool, :checkbox, true, (), :basic, true),
  WebUIOptionSpec("power_flow.apslf.nr_polish", "power_flow_apslf_nr_polish", Bool, :checkbox, true, (), :basic, true),
  WebUIOptionSpec("power_flow.apslf_start.enabled", "power_flow_apslf_start_enabled", Bool, :checkbox, false, (), :basic, true),
  WebUIOptionSpec("power_flow.apslf_start.order", "power_flow_apslf_start_order", Int, :number, 40, (), :basic, true),
  WebUIOptionSpec("power_flow.wrong_branch_detection", "power_flow_wrong_branch_detection", String, :select, "warn", WRONG_BRANCH_DETECTION_VALUES, :basic, true),
  WebUIOptionSpec("power_flow.start_mode.angle_mode", "power_flow_start_angle_mode", String, :select, "dc", POWERFLOW_START_ANGLE_MODE_VALUES, :basic, true),
  WebUIOptionSpec("power_flow.start_mode.voltage_mode", "power_flow_start_voltage_mode", String, :select, "profile_blend", POWERFLOW_START_VOLTAGE_MODE_VALUES, :basic, true),
  WebUIOptionSpec("power_flow.start_current_iteration.enabled", "power_flow_start_current_iteration_enabled", Bool, :checkbox, false, (), :basic, true),
  WebUIOptionSpec("power_flow.start_current_iteration.max_iter", "power_flow_start_current_iteration_max_iter", Int, :number, 10, (), :basic, true),
  WebUIOptionSpec("power_flow.start_current_iteration.tol", "power_flow_start_current_iteration_tol", Float64, :number, "1e-3", (), :basic, true),
  WebUIOptionSpec("power_flow.start_current_iteration.damping", "power_flow_start_current_iteration_damping", Float64, :number, 0.5, (), :basic, true),
  WebUIOptionSpec("power_flow.start_current_iteration.accept_only_if_improved", "power_flow_start_current_iteration_accept_only_if_improved", Bool, :checkbox, true, (), :basic, true),
  WebUIOptionSpec("power_flow.start_current_iteration.min_improvement_factor", "power_flow_start_current_iteration_min_improvement_factor", Float64, :number, 0.98, (), :basic, true),
  WebUIOptionSpec("power_flow.start_current_iteration.vm_min_pu", "power_flow_start_current_iteration_vm_min_pu", Float64, :number, 0.5, (), :basic, true),
  WebUIOptionSpec("power_flow.start_current_iteration.vm_max_pu", "power_flow_start_current_iteration_vm_max_pu", Float64, :number, 1.5, (), :basic, true),
  WebUIOptionSpec("power_flow.start_current_iteration.max_angle_step_deg", "power_flow_start_current_iteration_max_angle_step_deg", Float64, :number, 30.0, (), :basic, true),
  WebUIOptionSpec("power_flow.start_current_iteration.only_for_large_cases", "power_flow_start_current_iteration_only_for_large_cases", Bool, :checkbox, false, (), :basic, true),
  WebUIOptionSpec("power_flow.merit.enabled", "power_flow_merit_enabled", Bool, :checkbox, false, (), :basic, true),
  WebUIOptionSpec("power_flow.merit.armijo_c1", "power_flow_merit_armijo_c1", Float64, :number, "1e-4", (), :basic, true),
  WebUIOptionSpec("power_flow.merit.fallback_max_mismatch", "power_flow_merit_fallback_max_mismatch", Bool, :checkbox, true, (), :basic, true),
  WebUIOptionSpec("power_flow.trust_region.enabled", "power_flow_trust_region_enabled", Bool, :checkbox, false, (), :basic, true),
  WebUIOptionSpec("power_flow.trust_region.initial_radius", "power_flow_trust_region_initial_radius", Float64, :number, 1.0, (), :basic, true),
  WebUIOptionSpec("power_flow.trust_region.eta_accept", "power_flow_trust_region_eta_accept", Float64, :number, 0.1, (), :basic, true),
  WebUIOptionSpec("matpower_import.auto_profile", "matpower_import_auto_profile", String, :select, "off", MATPOWER_AUTO_PROFILE_VALUES, :expert, true),
  WebUIOptionSpec("matpower_import.ratio", "matpower_import_ratio", String, :select, "normal", MATPOWER_RATIO_VALUES, :expert, true),
  WebUIOptionSpec("matpower_import.shift_sign", "matpower_import_shift_sign", Float64, :number, 1.0, (), :expert, true),
  WebUIOptionSpec("matpower_import.shift_unit", "matpower_import_shift_unit", String, :select, "deg", MATPOWER_SHIFT_UNIT_VALUES, :expert, true),
  WebUIOptionSpec("matpower_import.bus_shunt_model", "matpower_import_bus_shunt_model", String, :select, "admittance", MATPOWER_BUS_SHUNT_MODEL_VALUES, :expert, true),
  WebUIOptionSpec("matpower_import.pv_voltage_source", "matpower_import_pv_voltage_source", String, :select, "gen_vg", MATPOWER_PV_VOLTAGE_SOURCE_VALUES, :expert, true),
  WebUIOptionSpec("matpower_import.compare_voltage_reference", "matpower_import_compare_voltage_reference", String, :select, "imported_setpoint", MATPOWER_COMPARE_VOLTAGE_REFERENCE_VALUES, :expert, true),
  WebUIOptionSpec("transformer.tap_changer_model", "transformer_tap_changer_model", String, :select, "ideal", TRANSFORMER_TAP_CHANGER_MODEL_VALUES, :expert, true),
  WebUIOptionSpec("matpower_export.write_solution", "matpower_export_write_solution", Bool, :checkbox, true, (), :expert, true),
  WebUIOptionSpec("output.logfile_results", "output_logfile_results", String, :select, "compact", OUTPUT_LOGFILE_RESULTS_VALUES, :basic, true),
  WebUIOptionSpec("benchmark.enabled", "benchmark_enabled", Bool, :checkbox, false, (), :expert, true),
  WebUIOptionSpec("benchmark.samples", "benchmark_samples", Int, :number, 10, (), :expert, true),
  WebUIOptionSpec("benchmark.seconds", "benchmark_seconds", Float64, :number, 1.0, (), :expert, true),
  WebUIOptionSpec(nothing, "performance_timing", String, :select, "compact", _WEBUI_PERFORMANCE_TIMING_VALUES, :basic, true),
  WebUIOptionSpec(nothing, "run_diagnostics", Bool, :checkbox, false, (), :expert, true),
  WebUIOptionSpec(nothing, "detailed_result_csv", Bool, :checkbox, true, (), :basic, true),
  WebUIOptionSpec(nothing, "detailed_result_csv_format", String, :select, "excel_us", ("technical", "excel_de", "excel_us"), :basic, true),
)

const _WEBUI_OPTION_BY_FIELD = Dict(spec.field => spec for spec in WEBUI_OPTION_SPECS)
const _WEBUI_FORM_CONFIG_FIELDS = Tuple((spec.config_key, spec.field, spec.value_type) for spec in WEBUI_OPTION_SPECS if spec.config_key !== nothing)
const _WEBUI_CASE_PROFILE_EXTRA_FIELDS = Tuple(spec.field for spec in WEBUI_OPTION_SPECS if spec.config_key === nothing && spec.save_in_case_sidecar)
const _WEBUI_CASE_PROFILE_FIELDS = Tuple(spec.field for spec in WEBUI_OPTION_SPECS if spec.save_in_case_sidecar)
const _WEBUI_CASE_PROFILE_FIELD_TYPES = Dict{String,Type}(
  spec.field => spec.value_type for spec in WEBUI_OPTION_SPECS if spec.save_in_case_sidecar
)

const _WEBUI_CASE_PROFILE_SELECT_VALUES = Dict{String,Set{String}}(
  spec.field => Set(string.(collect(spec.allowed_values))) for spec in WEBUI_OPTION_SPECS if spec.control == :select
)

_webui_option_spec(field::AbstractString)::WebUIOptionSpec = _WEBUI_OPTION_BY_FIELD[String(field)]
_webui_option_default(field::AbstractString) = _webui_option_spec(field).default
_webui_option_allowed_values(field::AbstractString) = _webui_option_spec(field).allowed_values

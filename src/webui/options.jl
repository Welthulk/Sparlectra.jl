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

# file: src/webui/options.jl
# purpose: central Web UI run-option metadata (WEBUI_OPTION_SPECS) mapping
#          form fields to configuration keys, defaults, and allowed values
"""Central metadata for one Web UI run-option form field."""
struct WebUIOptionSpec
  config_key::Union{Nothing,String}
  field::String
  value_type::Type
  control::Symbol
  default::Any
  allowed_values::Any
  section::Symbol
  # placement scope (stage 4B): :adapter = format-bound, generated on the
  # Case page from options_type(adapter); :session = machine scope (D8:
  # output/benchmark/runtime/webui prefixes), never written to a case
  # configuration file; :case = everything a run of one case may carry.
  # Visibility stays in `section` (:basic / :expert).
  scope::Symbol
  save_in_case_sidecar::Bool
end

const _WEBUI_QLIMIT_ENFORCEMENT_MODE_VALUES = (:active_set, :classic_simultaneous, :classic_one_at_a_time)

const _WEBUI_PERFORMANCE_TIMING_VALUES = WEBUI_PERFORMANCE_TIMING_VALUES

# Visibility (stage 4B, criterion decided 2026-09-02): :basic exactly for
# the config keys that appear in at least one workshop under docs/lit or
# in configuration.yaml.example without a default; every other config-key
# spec is :expert. Request-only fields (config_key nothing) keep their own
# visibility, the criterion covers the 63 config-backed form options.
const WEBUI_OPTION_SPECS = (
  WebUIOptionSpec("power_flow.mode", "power_flow_mode", String, :select, "manual", ("manual", "auto"), :expert, :case, true),
  WebUIOptionSpec("power_flow.tol", "power_flow_tol", Float64, :number, "1e-8", (), :basic, :case, true),
  # The UNIT of the tolerance value above. Not a configuration key of its
  # own: it decides which key the value becomes, power_flow.tol (pu) or
  # power_flow.tol_MW (megawatt, converted with the case base at run time).
  # One value field with a unit beats two fields side by side, which is how
  # this started and read like two competing tolerances (maintainer,
  # 2026-09-06).
  WebUIOptionSpec(nothing, "power_flow_tol_unit", String, :select, "pu", ("pu", "MW"), :basic, :case, true),
  WebUIOptionSpec("power_flow.max_iter", "power_flow_max_iter", Int, :number, 80, (), :basic, :case, true),
  WebUIOptionSpec("power_flow.autodamp", "power_flow_autodamp", Bool, :checkbox, true, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.autodamp_min", "power_flow_autodamp_min", Float64, :number, 0.05, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.qlimits.enabled", "power_flow_qlimits_enabled", Bool, :checkbox, true, (), :basic, :case, true),
  WebUIOptionSpec("power_flow.qlimits.enforcement_mode", "power_flow_qlimits_enforcement_mode", String, :select, "active_set", _WEBUI_QLIMIT_ENFORCEMENT_MODE_VALUES, :expert, :case, true),
  WebUIOptionSpec("power_flow.solver", "power_flow_solver", String, :select, "rectangular", POWERFLOW_SOLVER_VALUES, :expert, :case, true),
  WebUIOptionSpec("power_flow.apslf.order", "power_flow_apslf_order", Int, :number, 40, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.apslf.use_pade", "power_flow_apslf_use_pade", Bool, :checkbox, true, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.apslf.nr_polish", "power_flow_apslf_nr_polish", Bool, :checkbox, true, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.apslf_start.enabled", "power_flow_apslf_start_enabled", Bool, :checkbox, false, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.apslf_start.order", "power_flow_apslf_start_order", Int, :number, 40, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.wrong_branch_detection", "power_flow_wrong_branch_detection", String, :select, "warn", WRONG_BRANCH_DETECTION_VALUES, :expert, :case, true),
  WebUIOptionSpec("power_flow.start_mode.angle_mode", "power_flow_start_angle_mode", String, :select, "dc", POWERFLOW_START_ANGLE_MODE_VALUES, :expert, :case, true),
  WebUIOptionSpec("power_flow.start_mode.voltage_mode", "power_flow_start_voltage_mode", String, :select, "profile_blend", POWERFLOW_START_VOLTAGE_MODE_VALUES, :expert, :case, true),
  WebUIOptionSpec("power_flow.start_mode.dc_seed_unconditional", "power_flow_dc_seed_unconditional", Bool, :checkbox, false, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.start_current_iteration.enabled", "power_flow_start_current_iteration_enabled", Bool, :checkbox, false, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.start_current_iteration.max_iter", "power_flow_start_current_iteration_max_iter", Int, :number, 10, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.start_current_iteration.tol", "power_flow_start_current_iteration_tol", Float64, :number, "1e-3", (), :expert, :case, true),
  WebUIOptionSpec("power_flow.start_current_iteration.damping", "power_flow_start_current_iteration_damping", Float64, :number, 0.5, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.start_current_iteration.accept_only_if_improved", "power_flow_start_current_iteration_accept_only_if_improved", Bool, :checkbox, true, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.start_current_iteration.min_improvement_factor", "power_flow_start_current_iteration_min_improvement_factor", Float64, :number, 0.98, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.start_current_iteration.vm_min_pu", "power_flow_start_current_iteration_vm_min_pu", Float64, :number, 0.5, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.start_current_iteration.vm_max_pu", "power_flow_start_current_iteration_vm_max_pu", Float64, :number, 1.5, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.start_current_iteration.max_angle_step_deg", "power_flow_start_current_iteration_max_angle_step_deg", Float64, :number, 30.0, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.start_current_iteration.only_for_large_cases", "power_flow_start_current_iteration_only_for_large_cases", Bool, :checkbox, false, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.merit.enabled", "power_flow_merit_enabled", Bool, :checkbox, false, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.merit.armijo_c1", "power_flow_merit_armijo_c1", Float64, :number, "1e-4", (), :expert, :case, true),
  WebUIOptionSpec("power_flow.merit.fallback_max_mismatch", "power_flow_merit_fallback_max_mismatch", Bool, :checkbox, true, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.trust_region.enabled", "power_flow_trust_region_enabled", Bool, :checkbox, false, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.trust_region.initial_radius", "power_flow_trust_region_initial_radius", Float64, :number, 1.0, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.trust_region.eta_accept", "power_flow_trust_region_eta_accept", Float64, :number, 0.1, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.trust_region.step_mode", "power_flow_trust_region_step_mode", String, :select, "scaled", TRUST_REGION_STEP_MODE_VALUES, :expert, :case, true),
  # Startup option: unlike the run options around it, this one takes effect at
  # the NEXT Web UI start, so it is not stored per case.
  WebUIOptionSpec("power_flow.rescue", "power_flow_rescue", Bool, :checkbox, true, (), :expert, :case, true),
  WebUIOptionSpec("runtime.parallel.enabled", "runtime_parallel_enabled", Bool, :checkbox, true, (), :expert, :session, true),
  WebUIOptionSpec("power_flow.dc.fallback", "power_flow_dc_fallback", Bool, :checkbox, false, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.distributed_slack.enabled", "power_flow_distributed_slack_enabled", Bool, :checkbox, false, (), :basic, :case, true),
  WebUIOptionSpec("power_flow.distributed_slack.p_mode", "power_flow_distributed_slack_p_mode", String, :select, "pg_weighted", DISTRIBUTED_SLACK_P_MODE_VALUES, :basic, :case, true),
  WebUIOptionSpec("power_flow.external_grid.enabled", "power_flow_external_grid_enabled", Bool, :checkbox, false, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.external_grid.source", "power_flow_external_grid_source", String, :select, "auto", EXTERNAL_GRID_SOURCE_VALUES, :expert, :case, true),
  WebUIOptionSpec("power_flow.external_grid.sk_MVA", "power_flow_external_grid_sk_mva", Float64, :number, 2000.0, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.external_grid.rx", "power_flow_external_grid_rx", Float64, :number, 0.1, (), :expert, :case, true),
  WebUIOptionSpec("power_flow.linear_solver", "power_flow_linear_solver", String, :select, "umfpack", POWERFLOW_LINEAR_SOLVER_VALUES, :expert, :case, true),
  WebUIOptionSpec("cgmes_import.start_values", "cgmes_start_values", String, :select, "auto", CGMES_START_VALUES_VALUES, :expert, :adapter, true),
  WebUIOptionSpec("cgmes_import.require_boundary", "cgmes_require_boundary", Bool, :checkbox, true, (), :basic, :adapter, true),
  WebUIOptionSpec("cgmes_import.infer_base_voltages", "cgmes_infer_base_voltages", Bool, :checkbox, false, (), :expert, :adapter, true),
  WebUIOptionSpec("cgmes_import.hvdc_mode", "cgmes_hvdc_mode", String, :select, "injections", CGMES_HVDC_MODE_VALUES, :basic, :adapter, true),
  WebUIOptionSpec("matpower_import.matpower_dcline_mode", "matpower_import_dcline_mode", String, :select, "pf_injections", MATPOWER_DCLINE_MODE_VALUES, :basic, :adapter, true),
  WebUIOptionSpec("model.auto_profile", "matpower_import_auto_profile", String, :select, "off", MATPOWER_AUTO_PROFILE_VALUES, :expert, :adapter, true),
  WebUIOptionSpec("matpower_import.ratio", "matpower_import_ratio", String, :select, "normal", MATPOWER_RATIO_VALUES, :basic, :adapter, true),
  WebUIOptionSpec("matpower_import.shift_sign", "matpower_import_shift_sign", Float64, :number, 1.0, (), :basic, :adapter, true),
  WebUIOptionSpec("matpower_import.shift_unit", "matpower_import_shift_unit", String, :select, "deg", MATPOWER_SHIFT_UNIT_VALUES, :basic, :adapter, true),
  WebUIOptionSpec("model.bus_shunt_model", "matpower_import_bus_shunt_model", String, :select, "admittance", MATPOWER_BUS_SHUNT_MODEL_VALUES, :basic, :adapter, true),
  WebUIOptionSpec("matpower_import.pv_voltage_source", "matpower_import_pv_voltage_source", String, :select, "gen_vg", MATPOWER_PV_VOLTAGE_SOURCE_VALUES, :expert, :adapter, true),
  WebUIOptionSpec("matpower_import.compare_voltage_reference", "matpower_import_compare_voltage_reference", String, :select, "imported_setpoint", MATPOWER_COMPARE_VOLTAGE_REFERENCE_VALUES, :expert, :adapter, true),
  WebUIOptionSpec("matpower_import.apply_bus_names", "matpower_import_apply_bus_names", Bool, :checkbox, false, (), :expert, :adapter, true),
  WebUIOptionSpec("model.tap_changer_model", "transformer_tap_changer_model", String, :select, "ideal", TRANSFORMER_TAP_CHANGER_MODEL_VALUES, :basic, :adapter, true),
  WebUIOptionSpec("matpower_export.write_solution", "matpower_export_write_solution", Bool, :checkbox, true, (), :expert, :adapter, true),
  WebUIOptionSpec("output.logfile_results", "output_logfile_results", String, :select, "compact", OUTPUT_LOGFILE_RESULTS_VALUES, :expert, :session, true),
  WebUIOptionSpec("benchmark.enabled", "benchmark_enabled", Bool, :checkbox, false, (), :expert, :session, true),
  WebUIOptionSpec("benchmark.samples", "benchmark_samples", Int, :number, 10, (), :expert, :session, true),
  WebUIOptionSpec("benchmark.seconds", "benchmark_seconds", Float64, :number, 1.0, (), :expert, :session, true),
  # state-estimation run options and measurement-generator options: request
  # keys (not config overrides), persisted in the form block of the case
  # configuration file so a case reload restores what the user set on the
  # SE page.
  #
  # These defaults MUST match `state_estimation` in the configuration: they
  # are what a form shows before the user touches anything, and a value that
  # differs here silently outranks the configured one. That is exactly how a
  # 25000-bus run kept estimating with k_suppress 6.0 while the
  # configuration said 4.0 (task_se_bad_data_v0100). A test walks this list
  # against StateEstimationConfig.
  WebUIOptionSpec(nothing, "se_flatstart", Bool, :checkbox, true, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "se_tol", Float64, :number, 1.0e-6, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "se_max_iter", Int, :number, 50, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "se_robust_mode", String, :select, "off", ("off", "staged", "replacement"), :basic, :case, true),
  WebUIOptionSpec(nothing, "se_k_eliminate", Float64, :number, 3.0, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "se_robust_k1", Float64, :number, 3.0, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "se_robust_k2", Float64, :number, 6.0, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "se_k_suppress", Float64, :number, 4.0, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "se_suppression_sigma", Float64, :number, 2000.0, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "se_max_eliminations", Int, :number, 3, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "se_update_shunts", Bool, :checkbox, false, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "se_tap_estimation", Bool, :checkbox, false, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "se_report_correlation", Bool, :checkbox, false, (), :basic, :case, true),
  # default ON: a noise-free set puts J near 0 instead of near dof, which
  # reads like a broken statistic (the low-band note explains it, but the
  # recommended demo profile is noisy)
  WebUIOptionSpec(nothing, "noise", Bool, :checkbox, true, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "gross_error_k", Float64, :number, 0.0, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "gross_error_count", Int, :number, 1, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "tap_error_steps", Float64, :number, 0.0, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "tap_error_count", Int, :number, 1, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "sigma_u_pct", Float64, :number, 0.5, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "include_currents", Bool, :checkbox, false, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "sigma_i_pct", Float64, :number, 1.0, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "sigma_ia_deg", Float64, :number, 0.0, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "sigma_p_pct", Float64, :number, 1.0, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "sigma_q_pct", Float64, :number, 1.0, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "gen_truth_source", String, :select, "fresh_solve", ("fresh_solve", "from_run"), :basic, :case, true),
  WebUIOptionSpec(nothing, "gen_flow_ends", String, :select, "both", ("both", "one_balance_aware"), :basic, :case, true),
  WebUIOptionSpec(nothing, "gen_passive_sigma", Float64, :number, 0.05, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "gen_passive_as_zi", Bool, :checkbox, false, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "gen_seed", Int, :number, 42, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "performance_timing", String, :select, "compact", _WEBUI_PERFORMANCE_TIMING_VALUES, :basic, :case, true),
  WebUIOptionSpec(nothing, "detailed_result_csv", Bool, :checkbox, true, (), :basic, :case, true),
  WebUIOptionSpec(nothing, "detailed_result_csv_format", String, :select, "excel_us", ("technical", "excel_de", "excel_us"), :basic, :case, true),
  WebUIOptionSpec(nothing, "export_cgmes", Bool, :checkbox, false, (), :basic, :case, true),
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

# hygiene assert (stage 4B): the scope column must stay consistent with the
# central case-scope predicate. Every :session key is machine scope, and no
# case-config key may claim :session; a machine-scope key MAY sit in
# :adapter when its placement is format-bound (matpower_export.*: the
# checkbox lives on the Case page, the save drops it from case files
# through the same predicate).
for _spec in WEBUI_OPTION_SPECS
  _spec.config_key === nothing && continue
  is_case_key = scf_is_case_config_key(String(_spec.config_key))
  _spec.scope == :session && @assert !is_case_key "case-config key $(_spec.config_key) must not be :session scope"
  is_case_key && @assert _spec.scope != :session "machine-scope column mismatch for $(_spec.config_key)"
  !is_case_key && @assert _spec.scope in (:session, :adapter) "machine-scope key $(_spec.config_key) must be :session or :adapter"
end

_webui_option_spec(field::AbstractString)::WebUIOptionSpec = _WEBUI_OPTION_BY_FIELD[String(field)]
_webui_option_default(field::AbstractString) = _webui_option_spec(field).default
_webui_option_allowed_values(field::AbstractString) = _webui_option_spec(field).allowed_values

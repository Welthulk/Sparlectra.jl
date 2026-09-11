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

# file: src/webui/docs.jl
# purpose: Web UI in-app documentation: help-topic registry, markdown page
#          loading and section extraction, and doc-link rewriting
const _WEBUI_DOCS_ROOT = normpath(joinpath(@__DIR__, "..", "..", "docs", "src"))

const WEBUI_HELP_TOPICS = Dict(
  # state-estimation page (0.10.0): estimator options come from the config
  # table (one row per key), workflow topics from state_estimation.md
  "webui.se_measurement_file" => (label = "Measurement set (CSV v1)", page = "state_estimation", heading = "Measurement files and the SE chain", selector = ""),
  "webui.se_noise_seed" => (label = "Noise seed", page = "state_estimation", heading = "Case files carry their own measurements", selector = ""),
  "state_estimation.flatstart" => (label = "SE flat start", page = "state_estimation_configuration", heading = "State-Estimation Configuration", selector = "`state_estimation.flatstart`"),
  "state_estimation.tol" => (label = "SE convergence tolerance", page = "state_estimation_configuration", heading = "State-Estimation Configuration", selector = "`state_estimation.tol`"),
  "state_estimation.max_iter" => (label = "SE maximum iterations", page = "state_estimation_configuration", heading = "State-Estimation Configuration", selector = "`state_estimation.max_iter`"),
  "state_estimation.robust" => (label = "Robust R modification", page = "state_estimation", heading = "Bad-data thresholds and robust modes", selector = ""),
  "state_estimation.robust_mode" => (label = "Down-weighting (off / staged / replacement)", page = "state_estimation", heading = "Bad-data thresholds and robust modes", selector = ""),
  "state_estimation.k_eliminate" => (label = "Elimination limit (normalized residual)", page = "state_estimation", heading = "Bad-data thresholds and robust modes", selector = ""),
  "state_estimation.suppression" => (label = "Down-weight limit and replacement sigma", page = "state_estimation", heading = "Bad-data thresholds and robust modes", selector = ""),
  "webui.se_generator_truth" => (label = "Generator: truth state", page = "state_estimation", heading = "Measurement generator v2", selector = ""),
  "webui.se_generator_flow_ends" => (label = "Generator: flows per branch", page = "state_estimation", heading = "Measurement generator v2", selector = ""),
  "webui.se_generator_passive" => (label = "Generator: passive nodes", page = "state_estimation", heading = "Measurement generator v2", selector = ""),
  "state_estimation.update_shunts" => (label = "Write estimated shunts back", page = "state_estimation_configuration", heading = "State-Estimation Configuration", selector = "`state_estimation.update_shunts`"),
  "state_estimation.report_residual_correlation" => (label = "Residual-correlation (K-matrix) report", page = "state_estimation", heading = "Residual correlations (optional K-matrix report)", selector = ""),
  "webui.se_max_eliminations" => (label = "Sequential elimination budget", page = "state_estimation", heading = "Localizability: the residual sensitivity `wii`", selector = ""),
  "webui.se_generator_noise" => (label = "Generator: seeded noise", page = "state_estimation", heading = "Measurement files and the SE chain", selector = ""),
  "webui.se_generator_gross_error" => (label = "Generator: bad data (k·sigma)", page = "state_estimation", heading = "Measurement files and the SE chain", selector = ""),
  "webui.se_generator_tap_error" => (label = "Generator: tap deviation", page = "state_estimation", heading = "Measurement files and the SE chain", selector = ""),
  "webui.power_flow_mode" => (label = "Auto mode", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.mode`"),
  "webui.se_generator_gross_count" => (label = "Generator: bad data rows", page = "state_estimation", heading = "Measurement files and the SE chain", selector = ""),
  "webui.se_generator_tap_count" => (label = "Generator: tap deviation transformers", page = "state_estimation", heading = "Measurement files and the SE chain", selector = ""),
  "webui.se_tap_estimation" => (label = "Estimate transformer taps", page = "state_estimation", heading = "Transformer tap estimation", selector = ""),
  "webui.se_generator_sigmas" => (label = "Generator: measurement sigmas (U, I, P, Q)", page = "state_estimation", heading = "Measurement files and the SE chain", selector = ""),
  "webui.casefile" => (label = "MATPOWER case file", page = "webui", heading = "PowerFlow input paths", selector = "`webui.casefile`"),
  "webui.config_file" => (label = "Configuration template file", page = "webui", heading = "PowerFlow input paths", selector = "`webui.config_file`"),
  "webui.import_case_files" => (label = "Import case files", page = "webui", heading = "Importing case files through the Web UI", selector = ""),
  # the case-export controls: what the two files contain and what a plain PGM
  # export deliberately leaves out
  "webui.scf_export" => (label = "Case export (SCF / plain PGM)", page = "scf", heading = "Writing a case file", selector = ""),
  "webui.case_format" => (label = "Case input format", page = "webui", heading = "Starting a PowerFlow run", selector = ""),
  "webui.for002_reference_file" => (label = "Optional FOR002 reference file", page = "webui", heading = "Starting a PowerFlow run", selector = ""),
  "webui.dtf_outage_selection" => (label = "Selected DTF outage labels/indices", page = "webui", heading = "Starting a PowerFlow run", selector = ""),
  "webui.config_maintenance" => (label = "Configuration maintenance", page = "webui", heading = "Configuration check and refresh", selector = ""),
  "webui.ignore_webui_settings" => (label = "Ignore Web UI settings and use configuration defaults", page = "webui", heading = "Configuration precedence and artifact downloads", selector = ""),
  "power_flow.tol" => (label = "PowerFlow tolerance", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.tol`"),
  "power_flow.tol_MW" => (label = "PowerFlow tolerance in MW", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.tol_MW`"),
  "power_flow.max_iter" => (label = "Maximum iterations", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.max_iter`"),
  "power_flow.autodamp" => (label = "Autodamping enabled", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.autodamp`"),
  "power_flow.autodamp_min" => (label = "Autodamping minimum", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.autodamp_min`"),
  "power_flow.qlimits.enabled" => (label = "Q-limit handling enabled", page = "powerflow_configuration", heading = "Q-limit options and guard", selector = "`power_flow.qlimits.enabled`"),
  "power_flow.qlimits.enforcement_mode" => (label = "Q-limit enforcement mode", page = "powerflow_configuration", heading = "Q-limit options and guard", selector = "`power_flow.qlimits.enforcement_mode`"),
  "power_flow.solver" => (label = "Solver", page = "powerflow_configuration", heading = "Solver selection (rectangular vs. APSLF)", selector = "`power_flow.solver`"),
  "power_flow.linear_solver" => (label = "Linear solver backend", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.linear_solver`"),
  "power_flow.apslf.order" => (label = "APSLF highest coefficient (order)", page = "powerflow_configuration", heading = "Solver selection (rectangular vs. APSLF)", selector = "`power_flow.apslf.order`"),
  "power_flow.apslf.use_pade" => (label = "APSLF Padé evaluation", page = "powerflow_configuration", heading = "Solver selection (rectangular vs. APSLF)", selector = "`power_flow.apslf.use_pade`"),
  "power_flow.apslf.nr_polish" => (label = "APSLF NR polish", page = "powerflow_configuration", heading = "Solver selection (rectangular vs. APSLF)", selector = "`power_flow.apslf.nr_polish`"),
  "power_flow.apslf_start.enabled" => (label = "Use APSLF start values", page = "powerflow_configuration", heading = "Solver selection (rectangular vs. APSLF)", selector = "`power_flow.apslf_start.enabled`"),
  "power_flow.apslf_start.order" => (label = "APSLF start highest coefficient (order)", page = "powerflow_configuration", heading = "Solver selection (rectangular vs. APSLF)", selector = "`power_flow.apslf_start.order`"),
  "power_flow.wrong_branch_detection" => (label = "Wrong-branch detection", page = "configuration", heading = "Wrong-branch detection semantics (rectangular PF)", selector = ""),
  "power_flow.start_mode.angle_mode" => (label = "Start angle mode", page = "powerflow_configuration", heading = "Start mode options", selector = "`power_flow.start_mode.angle_mode`"),
  "power_flow.start_mode.voltage_mode" => (label = "Start voltage mode", page = "powerflow_configuration", heading = "Start mode options", selector = "`power_flow.start_mode.voltage_mode`"),
  "power_flow.start_mode.dc_seed_unconditional" => (label = "DC start values", page = "powerflow_configuration", heading = "Start mode options", selector = "`power_flow.start_mode.dc_seed_unconditional`"),
  "power_flow.start_current_iteration.enabled" => (label = "Enable current-iteration pre-solve", page = "configuration", heading = "Complete default-key index", selector = ""),
  "power_flow.start_current_iteration.max_iter" => (label = "Current-iteration max iterations", page = "configuration", heading = "Complete default-key index", selector = ""),
  "power_flow.start_current_iteration.tol" => (label = "Current-iteration tolerance", page = "configuration", heading = "Complete default-key index", selector = ""),
  "power_flow.start_current_iteration.damping" => (label = "Current-iteration damping", page = "configuration", heading = "Complete default-key index", selector = ""),
  "power_flow.start_current_iteration.accept_only_if_improved" => (label = "Accept only if improved", page = "configuration", heading = "Complete default-key index", selector = ""),
  "power_flow.start_current_iteration.min_improvement_factor" => (label = "Minimum improvement factor", page = "configuration", heading = "Complete default-key index", selector = ""),
  "power_flow.start_current_iteration.vm_min_pu" => (label = "Minimum voltage guard [pu]", page = "configuration", heading = "Complete default-key index", selector = ""),
  "power_flow.start_current_iteration.vm_max_pu" => (label = "Maximum voltage guard [pu]", page = "configuration", heading = "Complete default-key index", selector = ""),
  "power_flow.start_current_iteration.max_angle_step_deg" => (label = "Maximum angle-step guard [deg]", page = "configuration", heading = "Complete default-key index", selector = ""),
  "power_flow.start_current_iteration.only_for_large_cases" => (label = "Only for large cases", page = "configuration", heading = "Complete default-key index", selector = ""),
  "power_flow.distributed_slack.enabled" => (label = "Distributed active-power slack enabled", page = "powerflow_configuration", heading = "Distributed active-power slack", selector = "`power_flow.distributed_slack.enabled`"),
  "power_flow.distributed_slack.p_mode" => (label = "Distributed-slack weight mode", page = "powerflow_configuration", heading = "Distributed active-power slack", selector = "`power_flow.distributed_slack.p_mode`"),
  "power_flow.external_grid.enabled" => (label = "External grid source enabled", page = "powerflow_configuration", heading = "External grid source", selector = "`power_flow.external_grid.enabled`"),
  "power_flow.external_grid.source" => (label = "External grid Sk''/R-X source", page = "powerflow_configuration", heading = "External grid source", selector = "`power_flow.external_grid.source`"),
  "power_flow.external_grid.sk_MVA" => (label = "External grid short-circuit power Sk'' [MVA]", page = "powerflow_configuration", heading = "External grid source", selector = "`power_flow.external_grid.sk_MVA`"),
  "power_flow.external_grid.rx" => (label = "External grid R/X ratio", page = "powerflow_configuration", heading = "External grid source", selector = "`power_flow.external_grid.rx`"),
  "power_flow.merit.enabled" => (label = "Enable Armijo merit-function line search", page = "powerflow_configuration", heading = "Merit-function line search options", selector = "`power_flow.merit.enabled`"),
  "power_flow.merit.armijo_c1" => (label = "Armijo sufficient-decrease constant", page = "powerflow_configuration", heading = "Merit-function line search options", selector = "`power_flow.merit.armijo_c1`"),
  "power_flow.merit.fallback_max_mismatch" => (label = "Merit fallback behavior", page = "powerflow_configuration", heading = "Merit-function line search options", selector = "`power_flow.merit.fallback_max_mismatch`"),
  "power_flow.trust_region.enabled" => (label = "Enable trust-region step control", page = "powerflow_configuration", heading = "Trust-region step control options", selector = "`power_flow.trust_region.enabled`"),
  "power_flow.trust_region.initial_radius" => (label = "Initial trust-region radius", page = "powerflow_configuration", heading = "Trust-region step control options", selector = "`power_flow.trust_region.initial_radius`"),
  "power_flow.trust_region.eta_accept" => (label = "Trust-region acceptance ratio (eta)", page = "powerflow_configuration", heading = "Trust-region step control options", selector = "`power_flow.trust_region.eta_accept`"),
  "power_flow.trust_region.step_mode" => (label = "Trust-region step mode", page = "powerflow_configuration", heading = "Trust-region step control options", selector = "`power_flow.trust_region.step_mode`"),
  "cgmes_import.start_values" => (label = "CGMES start values", page = "cgmes_import", heading = "Configuration (`cgmes_import`)", selector = "`cgmes_import.start_values`"),
  "cgmes_import.hvdc_mode" => (label = "HVDC converters", page = "cgmes_import", heading = "Configuration (`cgmes_import`)", selector = "`cgmes_import.hvdc_mode`"),
  "matpower_import.matpower_dcline_mode" => (label = "DC-line mode", page = "matpower_import", heading = "Option reference", selector = "`matpower_import.matpower_dcline_mode`"),
  "cgmes_import.require_boundary" => (label = "Require boundary set", page = "cgmes_import", heading = "Configuration (`cgmes_import`)", selector = "`cgmes_import.require_boundary`"),
  "cgmes_import.infer_base_voltages" => (label = "Infer missing base voltages", page = "cgmes_import", heading = "Configuration (`cgmes_import`)", selector = "`cgmes_import.infer_base_voltages`"),
  "power_flow.rescue" => (label = "Rescue ladder for failed AC solves", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.rescue`"),
  "runtime.parallel.enabled" => (label = "Parallel execution of independent work items", page = "performance_profiling", heading = "Runtime", selector = "`runtime.parallel.enabled`"),
  "power_flow.dc.fallback" => (label = "Standalone-DC fallback", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.dc.fallback`"),
  "model.auto_profile" => (label = "MATPOWER auto-profile", page = "matpower_import", heading = "Option reference", selector = "`model.auto_profile`"),
  "matpower_import.ratio" => (label = "Transformer ratio convention", page = "matpower_import", heading = "Option reference", selector = "`matpower_import.ratio`"),
  "matpower_import.apply_bus_names" => (label = "Apply bus names", page = "matpower_import", heading = "Option reference", selector = "`matpower_import.apply_bus_names`"),
  "matpower_import.shift_sign" => (label = "Phase-shift sign", page = "matpower_import", heading = "Option reference", selector = "`matpower_import.shift_sign`"),
  "matpower_import.shift_unit" => (label = "Phase-shift unit", page = "matpower_import", heading = "Option reference", selector = "`matpower_import.shift_unit`"),
  "model.bus_shunt_model" => (label = "Bus-shunt model", page = "matpower_import", heading = "Option reference", selector = "`model.bus_shunt_model`"),
  "matpower_import.pv_voltage_source" => (label = "PV voltage source", page = "matpower_import", heading = "Option reference", selector = "`matpower_import.pv_voltage_source`"),
  "matpower_import.compare_voltage_reference" => (label = "Voltage reference comparison", page = "matpower_import", heading = "Option reference", selector = "`matpower_import.compare_voltage_reference`"),
  "model.tap_changer_model" => (label = "Tap-changer model", page = "matpower_import", heading = "Option reference", selector = "`model.tap_changer_model`"),
  "matpower_export.write_solution" => (label = "Write solution into MATPOWER export", page = "matpower_import", heading = "Option reference", selector = "`matpower_export.write_solution`"),
  "output.logfile_results" => (label = "Logfile output mode", page = "performance_profiling", heading = "Output configuration", selector = "`output.logfile_results`"),
  "benchmark.enabled" => (label = "Enable benchmark measurements", page = "performance_profiling", heading = "Benchmark configuration", selector = "`benchmark.enabled`"),
  "benchmark.samples" => (label = "Benchmark samples (max. repeated measurements)", page = "performance_profiling", heading = "Benchmark configuration", selector = "`benchmark.samples`"),
  "benchmark.seconds" => (label = "Benchmark max. time budget [s]", page = "performance_profiling", heading = "Benchmark configuration", selector = "`benchmark.seconds`"),
  "webui.performance_timing" => (label = "Performance timing", page = "webui", heading = "Run artifacts and output modes", selector = ""),
  "webui.detailed_result_csv" => (label = "Detailed result CSV export", page = "webui", heading = "Run artifacts and output modes", selector = ""),
  "webui.detailed_result_csv_format" => (label = "Detailed CSV format", page = "webui", heading = "Run artifacts and output modes", selector = ""),
  "webui.export_cgmes" => (label = "CGMES export artifact", page = "cgmes_export", heading = "Export from the Web UI", selector = ""),
)

const WEBUI_FORM_HELP_TOPICS = Dict(
  # labels the case-export actions row, not a single input
  "scf_export" => "webui.scf_export",
  "measurement_file" => "webui.se_measurement_file",
  "noise_seed" => "webui.se_noise_seed",
  "se_flatstart" => "state_estimation.flatstart",
  "se_tol" => "state_estimation.tol",
  "se_max_iter" => "state_estimation.max_iter",
  "se_robust" => "state_estimation.robust",
  "se_robust_mode" => "state_estimation.robust_mode",
  "se_k_eliminate" => "state_estimation.k_eliminate",
  "se_robust_k1" => "state_estimation.robust_mode",
  "se_robust_k2" => "state_estimation.robust_mode",
  "se_k_suppress" => "state_estimation.suppression",
  "se_suppression_sigma" => "state_estimation.suppression",
  "se_generator_truth" => "webui.se_generator_truth",
  "se_generator_truth_run" => "webui.se_generator_truth",
  "se_generator_flow_ends" => "webui.se_generator_flow_ends",
  "se_generator_passive_sigma" => "webui.se_generator_passive",
  "se_generator_passive_zi" => "webui.se_generator_passive",
  "se_max_eliminations" => "webui.se_max_eliminations",
  "se_update_shunts" => "state_estimation.update_shunts",
  "se_report_correlation" => "state_estimation.report_residual_correlation",
  "se_generator_noise" => "webui.se_generator_noise",
  "se_generator_gross_error" => "webui.se_generator_gross_error",
  "se_generator_tap_error" => "webui.se_generator_tap_error",
  "power_flow_mode" => "webui.power_flow_mode",
  "se_generator_gross_count" => "webui.se_generator_gross_count",
  "se_generator_tap_count" => "webui.se_generator_tap_count",
  "se_tap_estimation" => "webui.se_tap_estimation",
  "se_generator_sigmas" => "webui.se_generator_sigmas",
  "casefile" => "webui.casefile",
  "config_file" => "webui.config_file",
  "casefiles" => "webui.import_case_files",
  "case_format" => "webui.case_format",
  "for002_reference_file" => "webui.for002_reference_file",
  "dtf_outage_selection" => "webui.dtf_outage_selection",
  "config_maintenance" => "webui.config_maintenance",
  "ignore_webui_settings" => "webui.ignore_webui_settings",
  "power_flow_tol" => "power_flow.tol",
  # its own topic: the field states the bound in MW, and pointing at the
  # per-unit entry sent the reader to the wrong row of the table
  "power_flow_tol_unit" => "power_flow.tol_MW",
  "power_flow_max_iter" => "power_flow.max_iter",
  "power_flow_autodamp" => "power_flow.autodamp",
  "power_flow_autodamp_min" => "power_flow.autodamp_min",
  "power_flow_qlimits_enabled" => "power_flow.qlimits.enabled",
  "power_flow_qlimits_enforcement_mode" => "power_flow.qlimits.enforcement_mode",
  "power_flow_solver" => "power_flow.solver",
  "power_flow_linear_solver" => "power_flow.linear_solver",
  "power_flow_apslf_order" => "power_flow.apslf.order",
  "power_flow_apslf_use_pade" => "power_flow.apslf.use_pade",
  "power_flow_apslf_nr_polish" => "power_flow.apslf.nr_polish",
  "power_flow_apslf_start_enabled" => "power_flow.apslf_start.enabled",
  "power_flow_apslf_start_order" => "power_flow.apslf_start.order",
  "power_flow_wrong_branch_detection" => "power_flow.wrong_branch_detection",
  "power_flow_start_angle_mode" => "power_flow.start_mode.angle_mode",
  "power_flow_start_voltage_mode" => "power_flow.start_mode.voltage_mode",
  "power_flow_dc_seed_unconditional" => "power_flow.start_mode.dc_seed_unconditional",
  "power_flow_start_current_iteration_enabled" => "power_flow.start_current_iteration.enabled",
  "power_flow_start_current_iteration_max_iter" => "power_flow.start_current_iteration.max_iter",
  "power_flow_start_current_iteration_tol" => "power_flow.start_current_iteration.tol",
  "power_flow_start_current_iteration_damping" => "power_flow.start_current_iteration.damping",
  "power_flow_start_current_iteration_accept_only_if_improved" => "power_flow.start_current_iteration.accept_only_if_improved",
  "power_flow_start_current_iteration_min_improvement_factor" => "power_flow.start_current_iteration.min_improvement_factor",
  "power_flow_start_current_iteration_vm_min_pu" => "power_flow.start_current_iteration.vm_min_pu",
  "power_flow_start_current_iteration_vm_max_pu" => "power_flow.start_current_iteration.vm_max_pu",
  "power_flow_start_current_iteration_max_angle_step_deg" => "power_flow.start_current_iteration.max_angle_step_deg",
  "power_flow_start_current_iteration_only_for_large_cases" => "power_flow.start_current_iteration.only_for_large_cases",
  "power_flow_distributed_slack_enabled" => "power_flow.distributed_slack.enabled",
  "power_flow_distributed_slack_p_mode" => "power_flow.distributed_slack.p_mode",
  "power_flow_external_grid_enabled" => "power_flow.external_grid.enabled",
  "power_flow_external_grid_source" => "power_flow.external_grid.source",
  "power_flow_external_grid_sk_mva" => "power_flow.external_grid.sk_MVA",
  "power_flow_external_grid_rx" => "power_flow.external_grid.rx",
  "power_flow_merit_enabled" => "power_flow.merit.enabled",
  "power_flow_merit_armijo_c1" => "power_flow.merit.armijo_c1",
  "power_flow_merit_fallback_max_mismatch" => "power_flow.merit.fallback_max_mismatch",
  "power_flow_trust_region_enabled" => "power_flow.trust_region.enabled",
  "power_flow_trust_region_initial_radius" => "power_flow.trust_region.initial_radius",
  "power_flow_trust_region_eta_accept" => "power_flow.trust_region.eta_accept",
  "power_flow_trust_region_step_mode" => "power_flow.trust_region.step_mode",
  "cgmes_start_values" => "cgmes_import.start_values",
  "cgmes_hvdc_mode" => "cgmes_import.hvdc_mode",
  "matpower_import_dcline_mode" => "matpower_import.matpower_dcline_mode",
  "cgmes_require_boundary" => "cgmes_import.require_boundary",
  "cgmes_infer_base_voltages" => "cgmes_import.infer_base_voltages",
  "power_flow_rescue" => "power_flow.rescue",
  "runtime_parallel_enabled" => "runtime.parallel.enabled",
  "power_flow_dc_fallback" => "power_flow.dc.fallback",
  "matpower_import_auto_profile" => "model.auto_profile",
  "matpower_import_ratio" => "matpower_import.ratio",
  "matpower_import_shift_sign" => "matpower_import.shift_sign",
  "matpower_import_shift_unit" => "matpower_import.shift_unit",
  "matpower_import_bus_shunt_model" => "model.bus_shunt_model",
  "matpower_import_pv_voltage_source" => "matpower_import.pv_voltage_source",
  "matpower_import_compare_voltage_reference" => "matpower_import.compare_voltage_reference",
  "matpower_import_apply_bus_names" => "matpower_import.apply_bus_names",
  "transformer_tap_changer_model" => "model.tap_changer_model",
  "matpower_export_write_solution" => "matpower_export.write_solution",
  "output_logfile_results" => "output.logfile_results",
  "benchmark_enabled" => "benchmark.enabled",
  "benchmark_samples" => "benchmark.samples",
  "benchmark_seconds" => "benchmark.seconds",
  "performance_timing" => "webui.performance_timing",
  "detailed_result_csv" => "webui.detailed_result_csv",
  "detailed_result_csv_format" => "webui.detailed_result_csv_format",
  "export_cgmes" => "webui.export_cgmes",
)

const WEBUI_HELP_EXCERPT_OVERRIDES = Dict(
  "webui.se_max_eliminations" => """
## Sequential elimination budget

Upper bound for the sequential bad-data elimination: while the chi-square band test fails AND a suspicious measurement (normalized residual at or above 3) is localizable (`wii` > 0.3), the worst row is deactivated and the estimation reruns, up to this many times. `0` disables elimination (diagnostics only).

Protected rows are never eliminated: zero-injection pseudo-measurements (`ZI`), derived shunt rows (`SHDERIV`), and `LINKAGG` cluster aggregates encode model knowledge, not telemetry.

The elimination trace (which row, normalized residual before, objective drop) lands in `se_diagnostics.md`.

Configuration key `state_estimation.max_eliminations` (default 3).
""",
  "webui.se_generator_noise" => """
## Generator: seeded noise

Adds Gaussian noise at the per-quantity sigmas to every generated value. The random generator is seeded, so regenerating with the same inputs reproduces the identical file.

Without noise the values are exact solutions of the power flow: the estimation then reports J close to 0 and the Wilson-Hilferty band test flags `:low` (residuals implausibly small against the declared sigmas). That is expected for a noise-free synthetic set, not an error.
""",
  "webui.se_generator_gross_error" => """
## Generator: bad data (k·sigma)

Corrupts ONE measurement, the first active-power flow row, by k times its sigma. `0` = off; `10` is a good value: one clearly detectable bad measurement (for example a stuck transducer).

Use it to exercise the bad-data workflow: the diagnostics should rank exactly this row first (largest normalized residual, localizable), the sequential elimination should remove it, and the robust option should suppress it without removal. The corrupted row is named in the confirmation message.
""",
  "webui.se_generator_tap_error" => """
## Generator: tap deviation

Generates the measurements from a network state whose FIRST in-service transformer runs the given number of MECHANICAL tap steps off the model position (0 = off; WHOLE steps only, a tap changer has no half positions). The deviation lives on the same tap-fraction grid the estimator's fixation uses, so it is exactly recoverable: tap estimation finds the step and the fixation run ends near J = 0. The model file is never modified; only the throwaway generation state is shifted. Machine (generator step-up) transformers are never chosen as the deviation target: the mass release skips them, so a deviation there could not be resolved.

The resulting measurement set is consistent in itself but disagrees with the model around that transformer: the band test reports `:high` with suspicious measurements clustered at the transformer. This is the test vector for bad-data localization and for the tap estimation (the "estimate taps" option of the estimator run resolves exactly this discrepancy). The affected transformer branch is named in the confirmation message.
""",
  "webui.se_tap_estimation" => """
## Estimate transformer taps

Releases the tap of every in-service transformer that carries a ratio tap changer as an additional estimation state (the model tap stays untouched). After the estimation converges, each tap is FIXED to its nearest mechanical step and one more run is solved in which the tap is no state variable any more.

The result page reports, per transformer, the continuous electrical step and the fixed mechanical step, plus J before versus after the fixation: a small J before with a large J after means the true position sits between mechanical steps; both small means the fixed step explains the measurements. The full table lands in `se_tap_estimates.csv`.

Estimating a tap needs redundancy around the transformer (a meshed path or measurements on both sides); a radial transformer with a single flow measurement makes the tap and the downstream voltage indistinguishable.

**Release guards.** Machine (generator step-up) transformers are never mass-released: their terminal voltage is set by the machine's AVR, not observed, so the tap state would absorb it (release such a transformer explicitly via `setTapEstimation!` when you really mean it). A bridge transformer whose cut-off side carries no voltage measurement is frozen at its current position (`radial_no_voltage_pin`), and any remaining tap state the measurement set cannot pin numerically is frozen too (`not_observable`). Frozen taps stay in the result table with their reason instead of vanishing silently; a fully frozen release behaves exactly like no release.
""",
  "webui.se_generator_sigmas" => """
## Generator: measurement sigmas (U, I, P, Q)

Measurement accuracy per quantity in PERCENT OF THE MEASURED VALUE (like a transducer accuracy class), written into the sigma column of every generated row and used for the noise, when enabled:

- **sigma U** (%): all voltage-magnitude rows. `0.5` corresponds to a class 0.5 device.
- **currents (I)** checkbox plus **sigma I** (%): branch current-magnitude rows at both branch ends, generated only when the checkbox is set. Currents are auxiliary in the estimator: gated at flat start and below 3 sigma, excluded from observability.
- **sigma P** (%): active-power injections AND branch flows.
- **sigma Q** (%): reactive-power injections AND branch flows.

Percent of reading keeps one setting meaningful across voltage levels: 1 percent of a 400 MW flow and of a 4 MW flow both get a class-appropriate sigma, where an absolute 2 MW sigma would be tight at 400 kV and absurd at 30 kV. Near-zero readings get a small per-type floor instead of a near-zero sigma (`measurementSigmaFloors`: 1e-4 pu, 0.05 MW/MVar, 0.1 A), modeling the range term of the transducer class, so a dead branch never enters with weight infinity and zero-injection buses do not stiffen the solve.

For a healthy noisy set the estimation should land near J = dof; sigmas that are too pessimistic against the enabled noise push the band test to `:low`, too optimistic ones to `:high`.
""",
  "power_flow.start_current_iteration.enabled" => """
## Enable current-iteration pre-solve

`power_flow.start_current_iteration.enabled` enables a guarded current-injection/current-iteration pre-solve before the Newton-Raphson power-flow solver starts.

This is not a separate power-flow solver and it does not replace Newton-Raphson. It is a start-value preconditioner: Sparlectra first builds the initial voltage profile from Start Voltage Mode and Start Angle Mode, then optionally tries a few current-iteration steps to improve that initial profile.

The improved voltage profile is accepted only if it passes the voltage and angle guards and improves the existing Sparlectra mismatch metric. If it does not improve the start, the original start values are restored and Newton-Raphson starts normally.

Default: disabled. Enable this only for difficult cases where the normal start profile or DC/profile-blend start is not robust enough.

Diagnostic artifact: `current_iteration_start.log`
""",
  "power_flow.start_current_iteration.max_iter" => """
## Current-iteration max iterations

`power_flow.start_current_iteration.max_iter` sets the maximum number of current-iteration pre-solve steps before Newton-Raphson starts.

A higher value gives the pre-solve more chances to reduce the initial mismatch, but it also costs extra time and may move the start profile too far away from the original initialization. The result is still guarded and will be rejected if it becomes implausible or does not improve the mismatch.

Default: 10. Keep this small. Increase it only when diagnostics show that the mismatch keeps improving but the pre-solve stops too early.
""",
  "power_flow.start_current_iteration.tol" => """
## Current-iteration tolerance

`power_flow.start_current_iteration.tol` sets the stopping tolerance for the current-iteration pre-solve.

If the current-iteration mismatch or update criterion falls below this tolerance, the pre-solve can stop before reaching the maximum number of iterations. This tolerance only controls the start-value pre-solve. It is not the final Newton-Raphson power-flow tolerance.

Default: 1.0e-3. Use a relatively loose value; the purpose is to improve the starting point, not to solve the final power flow.
""",
  "power_flow.start_current_iteration.damping" => """
## Current-iteration damping

`power_flow.start_current_iteration.damping` sets the damping factor for the current-iteration voltage update.

A value of 1.0 applies the full current-iteration update. Smaller values blend the update with the previous voltage and make the pre-solve more conservative. This can help avoid large voltage or angle jumps in difficult cases.

Default: 0.5. Lower it if the pre-solve is rejected by voltage or angle guards. Increase it only if the pre-solve is stable but improves too slowly.
""",
  "power_flow.start_current_iteration.accept_only_if_improved" => """
## Accept only if improved

`power_flow.start_current_iteration.accept_only_if_improved` controls whether the current-iteration candidate is accepted only when it improves the existing Sparlectra mismatch metric.

When enabled, the pre-solve is conservative: if the candidate start profile is not better than the original start profile, Sparlectra restores the original start values before Newton-Raphson starts.

Default: enabled. This should normally stay enabled. Disabling it is only useful for expert experiments because it can allow a worse start profile to enter Newton-Raphson.
""",
  "power_flow.start_current_iteration.min_improvement_factor" => """
## Minimum improvement factor

`power_flow.start_current_iteration.min_improvement_factor` sets the required improvement ratio for accepting the current-iteration candidate when **Accept only if improved** is enabled.

For example, 0.98 means the candidate mismatch must be at least about 2% lower than the original mismatch. Smaller values require a stronger improvement; values closer to 1.0 accept smaller improvements.

Default: 0.98. Keep this close to 1.0 for a conservative pre-solve. Lower it only when tiny improvements are not useful and you want to accept only clearly better starts.
""",
  "power_flow.start_current_iteration.vm_min_pu" => """
## Minimum voltage guard

`power_flow.start_current_iteration.vm_min_pu` sets the lower voltage-magnitude guard for accepting a current-iteration candidate.

If any candidate bus voltage falls below this value, the candidate is rejected and the original start values are restored. This prevents the pre-solve from sending Newton-Raphson into an implausible low-voltage start region.

Default: 0.5 pu. Lowering this value makes the guard more permissive; increasing it makes the pre-solve more conservative. Use `current_iteration_start.log` to see candidate voltage minima before changing this value.
""",
  "power_flow.start_current_iteration.vm_max_pu" => """
## Maximum voltage guard

`power_flow.start_current_iteration.vm_max_pu` sets the upper voltage-magnitude guard for accepting a current-iteration candidate.

If any candidate bus voltage exceeds this value, the candidate is rejected and the original start values are restored. This prevents unrealistic over-voltage start profiles from entering Newton-Raphson.

Default: 1.5 pu. Lowering this value makes the guard stricter; increasing it allows larger candidate voltages. Use `current_iteration_start.log` to see candidate voltage maxima before changing this value.
""",
  "power_flow.start_current_iteration.max_angle_step_deg" => """
## Maximum angle-step guard

`power_flow.start_current_iteration.max_angle_step_deg` sets the maximum allowed angle change during the current-iteration pre-solve.

If the candidate introduces an angle jump larger than this limit, the candidate is rejected and the original start values are restored. This guard is intended to prevent unstable or wrong-branch start profiles.

Default: 30 degrees. Lower it for a more conservative pre-solve. Increase it only when diagnostics show that otherwise plausible candidates are rejected solely by the angle-step guard.
""",
  "power_flow.start_current_iteration.only_for_large_cases" => """
## Only for large cases

`power_flow.start_current_iteration.only_for_large_cases` runs the current-iteration pre-solve only for cases that Sparlectra classifies as large enough for this extra start-value preparation.

This avoids spending time on small cases where normal start values usually work well and where the pre-solve is not needed. The exact large-case threshold follows the existing Sparlectra configuration logic.

Default: disabled. Enable this if you want current iteration available for difficult large MATPOWER cases without changing behavior for small examples.
""",
  "power_flow.apslf_start.enabled" => """
## Use APSLF start values

`power_flow.apslf_start.enabled` uses the AnalyticLoadFlow.jl-backed APSLF solver as a guarded start-value generator ahead of the rectangular Newton-Raphson solve, the same insertion point and accept/reject guard style as the current-iteration pre-solve: the candidate is only adopted when it strictly improves the rectangular mismatch, otherwise the original start values are restored.

This mode always runs with **no NR polish** and **no Q-limit enforcement**, and neither is configurable here:

- NR polish is always off internally (`nr_polish=false`) because the downstream rectangular Newton-Raphson solve performs that polishing step itself.
- Q-limits are always unconstrained during this pre-solve, independent of `power_flow.qlimits.enabled` or any other Q-limit setting. `power_flow.qlimits.*` only governs the rectangular NR solve that follows; the generator's only job is producing a better starting voltage profile, not enforcing reactive limits.

Requires AnalyticLoadFlow.jl to be loaded; mutually exclusive with `power_flow.solver = apslf` (rejected at configuration time — the start-value generator only makes sense ahead of the NR solve).

Default: disabled. Diagnostic artifact: `apslf_start.log`.
""",
  "power_flow.apslf_start.order" => """
## APSLF start highest coefficient (order)

`power_flow.apslf_start.order` sets the highest power-series coefficient used by the APSLF start-value generator (see **Use APSLF start values**). Same considerations as `power_flow.apslf.order`: higher orders can improve the series approximation but cost more before the candidate is even evaluated for acceptance.

This option has no effect unless `power_flow.apslf_start.enabled = true`.

Default: 40.
""",
)

const WEBUI_DOC_PAGES = Dict(
  "configuration" => (title = "Configuration", file = "configuration.md"),
  "powerflow_configuration" => (title = "Power-Flow Configuration", file = "powerflow_configuration.md"),
  "powerflow_service" => (title = "Local PowerFlow Service", file = "powerflow_service.md"),
  "q_limit_switching_strategy" => (title = "Q-limit Switching Strategy", file = "q_limit_switching_strategy.md"),
  "performance_profiling" => (title = "Performance and Profiling Configuration", file = "performance_profiling.md"),
  "matpower_format" => (title = "MATPOWER format", file = "matpower_format.md"),
  "dtf_format" => (title = "DTF legacy input format", file = "dtf_format.md"),
  "scf" => (title = "Sparlectra Case Format (SCF/PGM JSON)", file = "scf.md"),
  "matpower_import" => (title = "MATPOWER Import", file = "matpower_import.md"),
  # Reachable from the docs reader; per-option help topics follow once the
  # cgmes_import options get their own Web UI form fields (issue #294).
  "cgmes_import" => (title = "CGMES Import", file = "cgmes_import.md"),
  "cgmes_export" => (title = "CGMES Export", file = "cgmes_export.md"),
  "webui" => (title = "Local PowerFlow Web UI", file = "webui.md"),
  "state_estimation" => (title = "State Estimation", file = "state_estimation.md"),
  "state_estimation_configuration" => (title = "State-Estimation Configuration", file = "state_estimation_configuration.md"),
  "feature_matrix" => (title = "Feature Matrix", file = "feature_matrix.md"),
  "solver" => (title = "Solver", file = "solver.md"),
)

"""Resolve an allowlisted Web UI help topic to its documentation metadata."""
resolve_webui_help_topic(topic::AbstractString) = get(WEBUI_HELP_TOPICS, String(topic), nothing)

"""Resolve an allowlisted documentation page to its title and Markdown file."""
resolve_webui_doc_page(page::AbstractString) = get(WEBUI_DOC_PAGES, String(page), nothing)

function _webui_doc_path(page_metadata)::String
  path = normpath(joinpath(_WEBUI_DOCS_ROOT, page_metadata.file))
  dirname(path) == _WEBUI_DOCS_ROOT || throw(ArgumentError("Documentation path is outside the allowlisted documentation root."))
  return path
end

"""Load one allowlisted Markdown document used by the local Web UI."""
function load_webui_markdown_document(page::AbstractString)::Union{String,Nothing}
  metadata = resolve_webui_doc_page(page)
  metadata === nothing && return nothing
  path = _webui_doc_path(metadata)
  return isfile(path) ? read(path, String) : nothing
end

function _webui_markdown_heading(line::AbstractString)
  matched = match(r"^(#{1,6})\s+(.+?)\s*#*\s*$", line)
  matched === nothing && return nothing
  return (level = length(matched.captures[1]), text = strip(matched.captures[2]))
end

"""Extract a Markdown heading and its content through the next peer or parent heading."""
function extract_webui_markdown_section(markdown_text::AbstractString, heading::AbstractString)::Union{String,Nothing}
  lines = split(String(markdown_text), '\n'; keepempty = true)
  start_index = nothing
  heading_level = 0
  for index in eachindex(lines)
    parsed = _webui_markdown_heading(lines[index])
    if parsed !== nothing && parsed.text == String(heading)
      start_index = index
      heading_level = parsed.level
      break
    end
  end
  start_index === nothing && return nothing

  stop_index = lastindex(lines)
  for index in (start_index + 1):lastindex(lines)
    parsed = _webui_markdown_heading(lines[index])
    if parsed !== nothing && parsed.level <= heading_level
      stop_index = index - 1
      break
    end
  end
  return strip(join(lines[start_index:stop_index], "\n"))
end

function _webui_extract_markdown_table_row(section::AbstractString, selector::AbstractString)::Union{String,Nothing}
  lines = split(String(section), '\n'; keepempty = true)
  row_index = findfirst(line -> startswith(strip(line), "|") && occursin(selector, line), lines)
  row_index === nothing && return nothing
  header_indices = findall(line -> startswith(strip(line), "|"), lines[begin:(row_index - 1)])
  length(header_indices) >= 2 || return strip(lines[row_index])
  heading_index = findfirst(line -> _webui_markdown_heading(line) !== nothing, lines)
  excerpt = String[]
  heading_index !== nothing && push!(excerpt, lines[heading_index], "")
  append!(excerpt, (lines[header_indices[1]], lines[header_indices[2]], lines[row_index]))
  return join(excerpt, "\n")
end

"""Load the Markdown excerpt configured for an allowlisted Web UI help topic."""
function load_webui_help_excerpt(topic::AbstractString)::Union{String,Nothing}
  override = get(WEBUI_HELP_EXCERPT_OVERRIDES, String(topic), nothing)
  override !== nothing && return override
  metadata = resolve_webui_help_topic(topic)
  metadata === nothing && return nothing
  markdown_text = load_webui_markdown_document(metadata.page)
  markdown_text === nothing && return nothing
  section = extract_webui_markdown_section(markdown_text, metadata.heading)
  section === nothing && return nothing
  return isempty(metadata.selector) ? section : _webui_extract_markdown_table_row(section, metadata.selector)
end

function _webui_heading_slug(heading_html::AbstractString)::String
  text = replace(String(heading_html), r"<[^>]+>" => "")
  # Julia's Markdown writes parentheses and similar characters as numeric
  # entities (`&#40;`); decoded first, they separate words like any other
  # punctuation instead of leaving their code points in the slug
  text = replace(text, r"&#(\d+);" => m -> string(Char(parse(Int, m[3:prevind(m, lastindex(m))]))))
  text = lowercase(replace(text, "&amp;" => "and", "&quot;" => "", "'" => ""))
  return strip(replace(text, r"[^a-z0-9]+" => "-"), '-')
end

function _webui_rewritten_doc_href(target::AbstractString; current_page::Union{Nothing,String} = nothing)::Union{String,Nothing}
  href = String(target)
  startswith(href, "https://") && return href
  startswith(href, "http://") && return href
  # heading ids are lowercase here (see rewrite_webui_doc_links), while the
  # pages link Documenter-style anchors such as `#Configuration-precedence`
  startswith(href, "#") && return current_page === nothing ? nothing : lowercase(href)
  (startswith(href, "/") || occursin('\\', href) || occursin(':', href)) && return nothing

  relative = startswith(href, "./") ? href[3:end] : href
  matched = match(r"^([A-Za-z0-9_-]+)\.md(#[A-Za-z0-9._~:%-]+)?$", relative)
  matched === nothing && return nothing
  page = matched.captures[1]
  metadata = resolve_webui_doc_page(page)
  metadata === nothing && return nothing
  metadata.file == "$(page).md" || return nothing
  fragment = lowercase(something(matched.captures[2], ""))
  return "/docs/$(page)$(fragment)"
end

function _webui_doc_link_attributes(href::AbstractString)::String
  return startswith(href, "https://") || startswith(href, "http://") ? " target=\"_blank\" rel=\"noopener noreferrer\"" : ""
end

"""Rewrite rendered Markdown links to safe, allowlisted local documentation routes."""
function rewrite_webui_doc_links(rendered_html::AbstractString; current_page::Union{Nothing,String} = nothing)::String
  heading_pattern = r"<h([1-6])>(.*?)</h[1-6]>"s
  html = replace(String(rendered_html), heading_pattern => matched_text -> begin
    matched = match(heading_pattern, String(matched_text))
    level = matched.captures[1]
    contents = matched.captures[2]
    slug = _webui_heading_slug(contents)
    isempty(slug) ? String(matched_text) : "<h$(level) id=\"$(slug)\">$(contents)</h$(level)>"
  end)
  href_pattern = Regex("href=\"([^\"]+)\"")
  return replace(html, href_pattern => matched_text -> begin
    matched = match(href_pattern, String(matched_text))
    rewritten = _webui_rewritten_doc_href(matched.captures[1]; current_page = current_page)
    rewritten === nothing ? "aria-disabled=\"true\"" : "href=\"$(rewritten)\"$(_webui_doc_link_attributes(rewritten))"
  end)
end

## LaTeX math for the local documentation viewer (maintainer 2026-09-04:
## the formulas showed as raw markup). Julia's Markdown standard library
## escapes `$...$` and ```math blocks into literal text, so nothing in the
## HTML tells a browser that this is mathematics. The Documenter website
## renders it with KaTeX; the local viewer has no such asset and must not
## grow a CDN dependency, so the common notation is translated into
## Unicode BEFORE parsing and handed over as inline code. Anything not in
## the table survives verbatim, which is still better than a dollar sign.
const _WEBUI_MATH_SYMBOLS = [
  # multi-character names first, so \varepsilon is not eaten by \var
  "\\varepsilon" => "ε", "\\epsilon" => "ε", "\\vartheta" => "ϑ", "\\theta" => "θ",
  "\\alpha" => "α", "\\beta" => "β", "\\gamma" => "γ", "\\delta" => "δ",
  "\\zeta" => "ζ", "\\eta" => "η", "\\iota" => "ι", "\\kappa" => "κ",
  "\\lambda" => "λ", "\\mu" => "μ", "\\nu" => "ν", "\\xi" => "ξ",
  "\\rho" => "ρ", "\\sigma" => "σ", "\\tau" => "τ", "\\upsilon" => "υ",
  "\\phi" => "φ", "\\chi" => "χ", "\\psi" => "ψ", "\\omega" => "ω",
  "\\Gamma" => "Γ", "\\Delta" => "Δ", "\\Theta" => "Θ", "\\Lambda" => "Λ",
  "\\Xi" => "Ξ", "\\Pi" => "Π", "\\Sigma" => "Σ", "\\Upsilon" => "Υ",
  "\\Phi" => "Φ", "\\Psi" => "Ψ", "\\Omega" => "Ω", "\\pi" => "π",
  "\\cdot" => "·", "\\times" => "×", "\\pm" => "±", "\\mp" => "∓",
  "\\leq" => "≤", "\\le" => "≤", "\\geq" => "≥", "\\ge" => "≥",
  "\\neq" => "≠", "\\ne" => "≠", "\\approx" => "≈", "\\equiv" => "≡",
  "\\infty" => "∞", "\\partial" => "∂", "\\nabla" => "∇", "\\sum" => "Σ",
  "\\prod" => "Π", "\\int" => "∫", "\\in" => "∈", "\\notin" => "∉",
  "\\subset" => "⊂", "\\subseteq" => "⊆", "\\cup" => "∪", "\\cap" => "∩",
  "\\rightarrow" => "→", "\\to" => "→", "\\leftarrow" => "←",
  "\\Rightarrow" => "⇒", "\\Leftrightarrow" => "⇔", "\\mapsto" => "↦",
  "\\ldots" => "…", "\\dots" => "…", "\\quad" => " ", "\\," => " ",
  "\\;" => " ", "\\!" => "", "\\left" => "", "\\right" => "",
]

const _WEBUI_MATH_SUPERSCRIPT = Dict('0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴', '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹', '+' => '⁺', '-' => '⁻', 'n' => 'ⁿ', 'i' => 'ⁱ', 'T' => 'ᵀ')
const _WEBUI_MATH_SUBSCRIPT = Dict('0' => '₀', '1' => '₁', '2' => '₂', '3' => '₃', '4' => '₄', '5' => '₅', '6' => '₆', '7' => '₇', '8' => '₈', '9' => '₉', '+' => '₊', '-' => '₋', 'i' => 'ᵢ', 'j' => 'ⱼ', 'k' => 'ₖ', 'n' => 'ₙ', 'm' => 'ₘ', 'a' => 'ₐ', 'x' => 'ₓ')

## One math expression to readable Unicode. Scripts translate only when
## EVERY character has a Unicode form; otherwise the plain ^/_ notation
## stays, which reads fine in a monospace span.
function _webui_math_to_unicode(expr::AbstractString)::String
  s = String(expr)
  for (from, to) in _WEBUI_MATH_SYMBOLS
    s = replace(s, from => to)
  end
  # \frac{a}{b} -> (a)/(b), \sqrt{a} -> √(a), \text{a} -> a
  # one level of nesting is allowed, so \frac{P_{se}}{V} resolves too;
  # deeper nesting stays as written rather than being mangled
  s = replace(s, r"\\frac\{((?:[^{}]|\{[^{}]*\})*)\}\{((?:[^{}]|\{[^{}]*\})*)\}" => s"(\1)/(\2)")
  s = replace(s, r"\\sqrt\{((?:[^{}]|\{[^{}]*\})*)\}" => s"√(\1)")
  s = replace(s, r"\\(?:text|mathrm|mathit|operatorname)\{([^{}]*)\}" => s"\1")
  # `replace` hands the matched TEXT to the function, not a RegexMatch, so
  # the body is cut out here instead of read from a capture group
  script = (matched, table, braced) -> begin
    body = braced ? matched[nextind(matched, firstindex(matched), 2):prevind(matched, lastindex(matched))] : matched[nextind(matched, firstindex(matched)):end]
    all(c -> haskey(table, c), body) ? join(table[c] for c in body) : matched
  end
  s = replace(s, r"\^\{[^{}]*\}" => m -> script(m, _WEBUI_MATH_SUPERSCRIPT, true))
  s = replace(s, r"_\{[^{}]*\}" => m -> script(m, _WEBUI_MATH_SUBSCRIPT, true))
  s = replace(s, r"\^\w" => m -> script(m, _WEBUI_MATH_SUPERSCRIPT, false))
  s = replace(s, r"_\w" => m -> script(m, _WEBUI_MATH_SUBSCRIPT, false))
  # leftover grouping braces disappear only when no unresolved command is
  # left; otherwise the expression stays readable as written
  occursin("\\", s) || (s = replace(s, r"[{}]" => ""))
  return strip(s)
end

## Replace ```math blocks and $...$ spans in a Markdown source with inline
## code carrying the Unicode form, so the standard Markdown renderer emits
## something readable instead of escaped dollars.
function _webui_render_math(markdown_text::AbstractString)::String
  s = String(markdown_text)
  s = replace(s, r"```math\n.*?```"s => m -> string("`", _webui_math_to_unicode(strip(m[8:prevind(m, lastindex(m), 3)])), "`\n"))
  # inline: no newline inside, and not an escaped dollar
  s = replace(s, r"(?<!\\)\$[^\$\n]+?(?<!\\)\$" => m -> string("`", _webui_math_to_unicode(m[nextind(m, firstindex(m)):prevind(m, lastindex(m))]), "`"))
  return s
end

## Documenter cross references (maintainer 2026-09-11: the CGMES page
## showed its "Node-breaker deliveries without a TP profile" heading as a
## dead link). `[text](@id name)` labels a heading and `[text](@ref name)`
## points at it, possibly from another page; Julia's Markdown renders both
## as ordinary links whose target the href rewriter then disables. The
## label is dropped from the heading text, and a reference becomes the
## page-local `#slug` or `page.md#slug` link the rewriter already handles.
## References without a known label (docstring refs such as
## [`addACLine!`](@ref), or ids on pages the viewer does not serve) stay as
## they are and render disabled, as before.
const _WEBUI_DOCUMENTER_ID_PATTERN = r"\[([^\]]+)\]\(@id\s+([A-Za-z0-9_.-]+)\)"
const _WEBUI_DOCUMENTER_REF_PATTERN = r"\[([^\]]+)\]\(@ref\s+([A-Za-z0-9_.-]+)\)"

"""Drop a Documenter `[text](@id name)` label from heading text, keeping the text."""
_webui_documenter_heading_text(text::AbstractString) = replace(String(text), _WEBUI_DOCUMENTER_ID_PATTERN => s"\1")

"""Anchor slug of a Markdown heading, the same one the rendered heading gets."""
function _webui_markdown_heading_slug(text::AbstractString)::String
  plain = _webui_documenter_heading_text(text)
  # the rendered heading carries the code spans as tags (stripped by the
  # slug) and the ampersand escaped (turned into "and")
  return _webui_heading_slug(replace(plain, "`" => "", "&" => "&amp;"))
end

"""Map every Documenter `@id` label on the served pages to `(page, slug)`."""
function _webui_documenter_id_index()::Dict{String,Tuple{String,String}}
  index = Dict{String,Tuple{String,String}}()
  for (page, metadata) in WEBUI_DOC_PAGES
    path = _webui_doc_path(metadata)
    isfile(path) || continue
    for line in eachline(path)
      heading = _webui_markdown_heading(line)
      heading === nothing && continue
      for matched in eachmatch(_WEBUI_DOCUMENTER_ID_PATTERN, heading.text)
        index[String(matched.captures[2])] = (page, _webui_markdown_heading_slug(heading.text))
      end
    end
  end
  return index
end

"""Resolve Documenter `@id`/`@ref` syntax into links the viewer's href rewriter understands."""
function _webui_resolve_documenter_refs(markdown_text::AbstractString; current_page::Union{Nothing,String} = nothing)::String
  s = replace(String(markdown_text), _WEBUI_DOCUMENTER_ID_PATTERN => s"\1")
  occursin("(@ref ", s) || return s
  index = _webui_documenter_id_index()
  return replace(s, _WEBUI_DOCUMENTER_REF_PATTERN => matched_text -> begin
    matched = match(_WEBUI_DOCUMENTER_REF_PATTERN, String(matched_text))
    target = get(index, String(matched.captures[2]), nothing)
    target === nothing && return String(matched_text)
    page, slug = target
    href = page == current_page ? "#$(slug)" : "$(page).md#$(slug)"
    "[$(matched.captures[1])]($(href))"
  end)
end

"""Render trusted repository Markdown as HTML using Julia's Markdown standard library."""
function render_webui_markdown(markdown_text::AbstractString; current_page::Union{Nothing,String} = nothing)::String
  io = IOBuffer()
  resolved = _webui_resolve_documenter_refs(String(markdown_text); current_page = current_page)
  show(io, MIME"text/html"(), Markdown.parse(_webui_render_math(resolved)))
  return rewrite_webui_doc_links(String(take!(io)); current_page = current_page)
end

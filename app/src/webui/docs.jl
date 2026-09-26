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
# purpose: Web UI help registry: the hover hint and the documentation link
#          of every form control. The Web UI carries no documentation of
#          its own (maintainer decision 2026-09-25): the hint is the
#          operating text, the `doc` target opens the published
#          documentation in a new tab. Every `doc` anchor is a Documenter
#          `@id` label in docs/src, checked by tools/check_webui_doc_links.jl
#          in the docs gate.

# Registry of help topics. `hint` is plain text, at most 160 characters,
# English, no backticks; `doc` is "<page>/#<anchor>" relative to the
# documentation base URL, or "" for a control with nothing technical behind it.
const WEBUI_HELP_TOPICS = Dict(
  # state-estimation page (0.10.0): estimator options come from the config
  # table (one row per key), workflow topics from state_estimation.md
  "webui.se_measurement_file" => (label = "Measurement set (CSV v1)", hint = "Measurement CSV v1 set the estimation runs on. The set bound to the selected case is preselected; sets from other cases are labeled.", doc = "state_estimation_measurements/#se-measurement-files"),
  "webui.se_noise_seed" => (label = "Noise seed", hint = "Seed of the noise draw. The same seed reproduces the same noisy set, a different seed draws a fresh realization.", doc = "state_estimation_measurements/#se-case-file-measurements"),
  "state_estimation.flatstart" => (label = "SE flat start", hint = "Start the estimator from a flat voltage profile instead of the model's current state. Uncheck when a trusted start state exists.", doc = "state_estimation_configuration/#se-config"),
  "state_estimation.tol" => (label = "SE convergence tolerance", hint = "WLS convergence tolerance on the state step. The default 1e-6 is the noise floor of the finite-difference Jacobian; a tighter value is raised with a log line.", doc = "state_estimation_configuration/#se-config"),
  "state_estimation.max_iter" => (label = "SE maximum iterations", hint = "Iteration cap of the estimator, the same limit the service and the configuration use. A solve with released transformer taps can need close to 40 iterations.", doc = "state_estimation_configuration/#se-config"),
  "state_estimation.robust" => (label = "Robust R modification", hint = "Deprecated alias for the down-weighting mode staged, read only while that mode is off. Set the down-weighting mode instead, never both.", doc = "state_estimation/#se-bad-data"),
  "state_estimation.robust_mode" => (label = "Down-weighting (off / staged / replacement)", hint = "Treatment of a large residual: off is plain WLS, staged down-weights the row gradually between two knees, replacement gives it a fixed large sigma.", doc = "state_estimation/#se-bad-data"),
  "state_estimation.k_eliminate" => (label = "Elimination limit (normalized residual)", hint = "Normalized residual from which a measurement is removed from the estimate, not just down-weighted. Applies while the band test is high, within the budget.", doc = "state_estimation/#se-bad-data"),
  "state_estimation.suppression" => (label = "Down-weight limit and replacement sigma", hint = "Replacement mode: from this normalized residual a row is solved with the replacement sigma and loses its influence, but is not removed. Sigma in the row's unit.", doc = "state_estimation/#se-bad-data"),
  "webui.se_generator_truth" => (label = "Generator: truth state", hint = "Where the truth state comes from: fresh solve solves the case now, from run adopts the solved voltages of a successful run of this case without re-solving.", doc = "state_estimation_measurements/#se-generator-v2"),
  "webui.se_generator_flow_ends" => (label = "Generator: flows per branch", hint = "Flows per branch: both ends writes P and Q at from and to; one end keeps one flow group per branch, preferring the end whose bus has injection telemetry.", doc = "state_estimation_measurements/#se-generator-v2"),
  "webui.se_generator_passive" => (label = "Generator: passive nodes", hint = "How a passive node enters the set: on writes a protected zero-injection constraint, off an ordinary zero balance row. The sigma applies to both, floor 1 kW.", doc = "state_estimation_measurements/#se-generator-v2"),
  "state_estimation.update_shunts" => (label = "Write estimated shunts back", hint = "Write the estimated shunt susceptances back into the model after a converged run. Off keeps the model data authoritative.", doc = "state_estimation_configuration/#se-config"),
  "state_estimation.report_residual_correlation" => (label = "Residual-correlation (K-matrix) report", hint = "Add residual-correlation columns to the bad-data report: rows correlated above 0.707 form a group where a gross error cannot be told apart. Reporting only.", doc = "state_estimation/#se-k-report"),
  "webui.se_max_eliminations" => (label = "Sequential elimination budget", hint = "Upper bound of the sequential bad-data elimination: the worst localizable row is removed and the estimation rerun, up to this many times. 0 disables it.", doc = "state_estimation/#se-max-eliminations"),
  "webui.se_generator_noise" => (label = "Generator: seeded noise", hint = "Add seeded Gaussian noise at the sigmas below. Without noise the measurements match the model exactly and J lands near 0 instead of near dof.", doc = "state_estimation_measurements/#se-generator-noise"),
  "webui.se_generator_gross_error" => (label = "Generator: bad data (k·sigma)", hint = "Corrupt seed-randomly drawn measurements by k times their sigma. 0 is off; 10 gives clearly detectable bad data for the elimination and robust workflow.", doc = "state_estimation_measurements/#se-generator-gross-error"),
  "webui.se_generator_tap_error" => (label = "Generator: tap deviation", hint = "Generate from a state whose seed-selected transformers run this many mechanical tap steps off the model position; 0 is off. Needs truth state fresh solve.", doc = "state_estimation_measurements/#se-generator-tap-error"),
  "webui.power_flow_mode" => (label = "Auto mode", hint = "Pick start-value, step-control and Q-limit strategy from the network, retrying with stronger strategies on non-convergence. Explicit options below always win.", doc = "powerflow_configuration/#pf-solver-core"),
  "webui.se_generator_gross_count" => (label = "Generator: bad data rows", hint = "How many measurement rows get the gross error; the seed decides which rows. Used only when bad data k is above 0.", doc = "state_estimation_measurements/#se-generator-gross-error"),
  "webui.se_generator_tap_count" => (label = "Generator: tap deviation transformers", hint = "At most this many transformers get the tap deviation; the seed picks them among transformers the tap estimation can absorb. Unused when the deviation is 0.", doc = "state_estimation_measurements/#se-generator-tap-error"),
  "webui.se_tap_estimation" => (label = "Estimate transformer taps", hint = "Release every in-service ratio tap changer as an extra estimation state, fix it to the nearest mechanical step and rerun. Machine transformers stay untouched.", doc = "state_estimation_extensions/#se-tap-estimation"),
  "webui.se_generator_sigmas" => (label = "Generator: measurement sigmas (U, I, P, Q)", hint = "Measurement accuracy per quantity in percent of the measured value, like a transducer class; written into the sigma column of every row and used for the noise.", doc = "state_estimation_measurements/#se-generator-sigmas"),
  "webui.casefile" => (label = "MATPOWER case file", hint = "Choose a case from the list, or type a bare case name or a local path and press Enter to download or copy it into the case directory.", doc = "webui/#webui-form-options"),
  "webui.config_file" => (label = "Configuration template file", hint = "Configuration YAML the run starts from. Form values become per-run overrides on top of it; the file itself stays unchanged.", doc = "webui/#webui-form-options"),
  "webui.import_case_files" => (label = "Import case files", hint = "Copy MATPOWER, DTF, CGMES (ZIP or profile XML files) or SCF/PGM JSON files into the case directory. Nothing is run; an existing file is never overwritten.", doc = "webui/#webui-import-case-files"),
  # the case-export controls: what the two files contain and what a plain PGM
  # export deliberately leaves out
  "webui.scf_export" => (label = "Case export (SCF / plain PGM)", hint = "Write the selected case into the case directory as an SCF file (.scf.json) or as a plain power-grid-model dataset (.pgm.json) without the sparlectra block.", doc = "scf/#scf-writing"),
  "webui.case_format" => (label = "Case input format", hint = "Auto tells MATPOWER, CGMES, SCF/PGM JSON and DTF apart from the file content. Pick a format only when the file does not say what it is.", doc = "webui/#webui-form-options"),
  "webui.for002_reference_file" => (label = "Optional FOR002 reference file", hint = "Optional FOR002 file for the legacy reference comparison of a DTF run: absolute path, path in the case cache, or an offered candidate. Not a primary case.", doc = "webui/#webui-form-options"),
  "webui.dtf_outage_selection" => (label = "Selected DTF outage labels/indices", hint = "One outage label or index of the DTF case, used with the run mode selected. The result page shows a compact outage summary; the rows stay in the artifacts.", doc = "webui/#webui-form-options"),
  "webui.config_maintenance" => (label = "Configuration maintenance", hint = "Check compares the YAML with the template and previews the refresh without writing; Refresh adds missing keys after a backup; the Editor edits the file.", doc = "webui/#webui-configuration"),
  "webui.ignore_webui_settings" => (label = "Ignore Web UI settings and use configuration defaults", hint = "Run with the configuration file's values only: the form values and the saved case settings are ignored for this run.", doc = "webui/#webui-configuration"),
  "power_flow.tol" => (label = "PowerFlow tolerance", hint = "Convergence bound for the largest single bus mismatch, active and reactive alike. 1e-8 pu equals 1 W at a 100 MVA base.", doc = "powerflow_configuration/#pf-solver-core"),
  "power_flow.tol_MW" => (label = "PowerFlow tolerance in MW", hint = "Unit of the tolerance value. MW states the bound physically, converted with the case base at run time (1 MW at 100 MVA is 1e-2 pu); pu is the classic form.", doc = "powerflow_configuration/#pf-solver-core"),
  "power_flow.max_iter" => (label = "Maximum iterations", hint = "Iteration cap of the Newton-Raphson solve; a run that reaches it ends as non-converged. Very low values cut off hard cases.", doc = "powerflow_configuration/#pf-solver-core"),
  "power_flow.autodamp" => (label = "Autodamping enabled", hint = "Adaptive damping of the Newton step. Helps difficult convergence at a small overhead; switch it off only for strict algorithm comparisons.", doc = "powerflow_configuration/#pf-solver-core"),
  "power_flow.autodamp_min" => (label = "Autodamping minimum", hint = "Minimum damping factor while autodamping is active. Lower values stabilize hard cases but can increase the iteration count.", doc = "powerflow_configuration/#pf-solver-core"),
  "power_flow.qlimits.enabled" => (label = "Q-limit handling enabled", hint = "Master switch for generator Q-limit enforcement (PV to PQ switching). Off gives an unconstrained power flow; on can add switching iterations.", doc = "powerflow_configuration/#pf-qlimits"),
  "power_flow.qlimits.enforcement_mode" => (label = "Q-limit enforcement mode", hint = "Selects active-set or classical Q-limit switching; the classical modes may need repeated solves. off disables the handling like the checkbox above.", doc = "powerflow_configuration/#pf-qlimits"),
  "power_flow.solver" => (label = "Solver", hint = "Selects the executing solver: AC Newton-Raphson, the APSLF power series (needs AnalyticLoadFlow.jl), or a DC linear model without voltages and reactive results.", doc = "powerflow_configuration/#pf-solver-selection"),
  "power_flow.linear_solver" => (label = "Linear solver backend", hint = "Sparse backend of the Newton step. umfpack_reuse keeps the symbolic analysis across iterations, faster on large cases; umfpack analyzes every iteration anew.", doc = "powerflow_configuration/#pf-solver-core"),
  "power_flow.apslf.order" => (label = "APSLF highest coefficient (order)", hint = "Highest power-series coefficient computed by the APSLF solver. Higher orders help stressed cases but increase the solve cost.", doc = "powerflow_configuration/#pf-solver-selection"),
  "power_flow.apslf.use_pade" => (label = "APSLF Padé evaluation", hint = "Evaluate the voltage series via Padé approximants instead of direct Taylor summation. Usually better accuracy per order at a small overhead.", doc = "powerflow_configuration/#pf-solver-selection"),
  "power_flow.apslf.nr_polish" => (label = "APSLF NR polish", hint = "Run a Newton-Raphson polishing step on the APSLF series result. Default off: the series alone is a load-flow solution.", doc = "powerflow_configuration/#pf-solver-selection"),
  "power_flow.apslf.convergence_radius" => (label = "APSLF convergence radius", hint = "Evaluate the APSLF convergence radius, the distance of the nearest Padé pole to the evaluation point. Costs about as much as the solve; shown with the result.", doc = "powerflow_configuration/#pf-solver-selection"),
  "power_flow.flatstart" => (label = "Flat start", hint = "Start every bus at 1.0 pu and 0 degrees and ignore the imported start voltages. While on, the APSLF and DC start values and the pre-solve are switched off.", doc = "webui/#webui-flat-start"),
  "power_flow.apslf_start.enabled" => (label = "Use APSLF start values", hint = "APSLF solver as a guarded start-value generator ahead of Newton-Raphson: the candidate is kept only if it improves the mismatch. Not with solver APSLF.", doc = "powerflow_configuration/#pf-apslf-start"),
  "power_flow.apslf_start.order" => (label = "APSLF start highest coefficient (order)", hint = "Highest series coefficient of the APSLF start-value generator; higher orders cost more before the candidate is judged. No effect unless the generator is on.", doc = "powerflow_configuration/#pf-apslf-start"),
  "power_flow.wrong_branch_detection" => (label = "Wrong-branch detection", hint = "Plausibility check of the converged solution on the highest voltage level: warn reports a suspicious result, fail treats it as non-converged, off skips it.", doc = "configuration/#config-wrong-branch"),
  "power_flow.start_mode.angle_mode" => (label = "Start angle mode", hint = "Source of the start angles: classic zero angles, a DC pre-pass (default), or the imported angles. Ignored while the flat start is on.", doc = "powerflow_configuration/#pf-start-mode"),
  "power_flow.start_mode.voltage_mode" => (label = "Start voltage mode", hint = "Source of the start voltage magnitudes: classic flat, generator or bus setpoints, or a blend of the imported profile. Ignored while the flat start is on.", doc = "powerflow_configuration/#pf-start-mode"),
  "power_flow.start_mode.dc_seed_unconditional" => (label = "DC start values", hint = "Always run a full DC power flow first and seed the Newton-Raphson start angles from it, bypassing the start angle mode and its quality gate.", doc = "powerflow_configuration/#pf-start-mode"),
  "power_flow.start_current_iteration.enabled" => (label = "Enable current-iteration pre-solve", hint = "Guarded current-iteration pre-solve improving the start profile before Newton-Raphson; the result is kept only if it passes the guards and lowers the mismatch.", doc = "powerflow_configuration/#pf-current-iteration-start"),
  "power_flow.start_current_iteration.max_iter" => (label = "Current-iteration max iterations", hint = "Maximum number of pre-solve steps before Newton-Raphson starts. Keep it small; raise it only when the mismatch keeps improving but the pre-solve stops early.", doc = "powerflow_configuration/#pf-current-iteration-start"),
  "power_flow.start_current_iteration.tol" => (label = "Current-iteration tolerance", hint = "Stopping tolerance of the pre-solve, not the final Newton-Raphson tolerance. A loose value is enough; the goal is a better start, not a solved power flow.", doc = "powerflow_configuration/#pf-current-iteration-start"),
  "power_flow.start_current_iteration.damping" => (label = "Current-iteration damping", hint = "Damping factor of the pre-solve voltage update; 1.0 applies the full update. Lower it if the candidate is rejected by the voltage or angle guards.", doc = "powerflow_configuration/#pf-current-iteration-start"),
  "power_flow.start_current_iteration.accept_only_if_improved" => (label = "Accept only if improved", hint = "Accept the pre-solve candidate only when it improves the mismatch; otherwise the original start values are restored. Keep it on outside expert experiments.", doc = "powerflow_configuration/#pf-current-iteration-start"),
  "power_flow.start_current_iteration.min_improvement_factor" => (label = "Minimum improvement factor", hint = "Required improvement ratio of the candidate mismatch; 0.98 means at least about 2 percent better. Smaller values demand a stronger improvement.", doc = "powerflow_configuration/#pf-current-iteration-start"),
  "power_flow.start_current_iteration.vm_min_pu" => (label = "Minimum voltage guard [pu]", hint = "Lower voltage guard: a candidate with any bus voltage below it is rejected and the original start values are restored.", doc = "powerflow_configuration/#pf-current-iteration-start"),
  "power_flow.start_current_iteration.vm_max_pu" => (label = "Maximum voltage guard [pu]", hint = "Upper voltage guard: a candidate with any bus voltage above it is rejected and the original start values are restored.", doc = "powerflow_configuration/#pf-current-iteration-start"),
  "power_flow.start_current_iteration.max_angle_step_deg" => (label = "Maximum angle-step guard [deg]", hint = "Largest allowed angle change of one pre-solve update; a larger jump rejects the candidate. Lower it for a more conservative pre-solve.", doc = "powerflow_configuration/#pf-current-iteration-start"),
  "power_flow.start_current_iteration.only_for_large_cases" => (label = "Only for large cases", hint = "Run the pre-solve only for cases classified as large. Small cases keep their normal start values.", doc = "powerflow_configuration/#pf-current-iteration-start"),
  "power_flow.distributed_slack.enabled" => (label = "Distributed active-power slack enabled", hint = "Spread the active-power imbalance over participating generators instead of the single slack bus. Off is bit-identical to the classical single-slack solver.", doc = "powerflow_configuration/#pf-distributed-slack"),
  "power_flow.distributed_slack.p_mode" => (label = "Distributed-slack weight mode", hint = "Weight source of the participation factors: scheduled Pg, maximum P, remaining headroom, the imported factors, or an explicit weights table from the YAML.", doc = "powerflow_configuration/#pf-distributed-slack"),
  "power_flow.external_grid.enabled" => (label = "External grid source enabled", hint = "Model the slack as a non-ideal source behind a feeder impedance, so its voltage droops under load. Off keeps the ideal slack; not with the distributed slack.", doc = "powerflow_configuration/#pf-external-grid"),
  "power_flow.external_grid.source" => (label = "External grid Sk''/R-X source", hint = "Where Sk'' and R/X come from: auto prefers values the case declares on the slack bus, else the numbers below; config always uses the numbers below.", doc = "powerflow_configuration/#pf-external-grid"),
  "power_flow.external_grid.sk_MVA" => (label = "External grid short-circuit power Sk'' [MVA]", hint = "Initial symmetrical short-circuit power of the feeder in MVA. It sets the series impedance of the external grid source (z = baseMVA/Sk'').", doc = "powerflow_configuration/#pf-external-grid"),
  "power_flow.external_grid.rx" => (label = "External grid R/X ratio", hint = "R/X ratio of the feeder impedance of the external grid source.", doc = "powerflow_configuration/#pf-external-grid"),
  "power_flow.merit.enabled" => (label = "Enable Armijo merit-function line search", hint = "Armijo merit-function line search for step acceptance on difficult flat-start cases. Requires autodamping; adds one residual-norm evaluation per trial.", doc = "powerflow_configuration/#pf-merit"),
  "power_flow.merit.armijo_c1" => (label = "Armijo sufficient-decrease constant", hint = "Sufficient-decrease constant of the Armijo condition. Larger values reject more trial steps and increase the backtracking.", doc = "powerflow_configuration/#pf-merit"),
  "power_flow.merit.fallback_max_mismatch" => (label = "Merit fallback behavior", hint = "When no trial satisfies the Armijo condition, fall back to the max-mismatch criterion (first improving trial). Off takes the most conservative finite trial.", doc = "powerflow_configuration/#pf-merit"),
  "power_flow.trust_region.enabled" => (label = "Enable trust-region step control", hint = "Scaled-Newton trust-region step control as an alternative to autodamping; enabling one disables the other. Suited to difficult flat-start cases.", doc = "powerflow_configuration/#pf-trust-region"),
  "power_flow.trust_region.initial_radius" => (label = "Initial trust-region radius", hint = "Starting trust-region radius in per-unit state-vector norm. Larger values risk more rejected or shrunk first steps on hard cases.", doc = "powerflow_configuration/#pf-trust-region"),
  "power_flow.trust_region.eta_accept" => (label = "Trust-region acceptance ratio (eta)", hint = "Minimum actual-to-predicted reduction ratio to accept a trial step. Higher values reject more trials and add shrink-and-retry iterations.", doc = "powerflow_configuration/#pf-trust-region"),
  "power_flow.trust_region.step_mode" => (label = "Trust-region step mode", hint = "Trial-step construction: scaled rescales the full Newton direction to the radius; dogleg blends it with a steepest-descent step when the radius shrinks.", doc = "powerflow_configuration/#pf-trust-region"),
  "cgmes_import.start_values" => (label = "CGMES start values", hint = "auto starts from the delivery's SvVoltage state when it carries one, else flat; sv always uses the imported state, flat always starts flat.", doc = "cgmes_import/#cgmes-import-config"),
  "cgmes_import.hvdc_mode" => (label = "HVDC converters", hint = "HVDC converter model: fixed injections reproduce the delivery snapshot; a paired controller couples both converters of a link and makes the transfer steerable.", doc = "cgmes_import/#cgmes-import-config"),
  "matpower_import.matpower_dcline_mode" => (label = "DC-line mode", hint = "Model of active mpc.dcline rows: pf_injections adds two fixed terminal injections per row; paired_control also couples each pair as a steerable HVDC link.", doc = "matpower/#matpower-options"),
  "cgmes_import.require_boundary" => (label = "Require boundary set", hint = "Fail the CGMES import when topology references stay unresolved (boundary set missing). Uncheck to import an incomplete delivery anyway.", doc = "cgmes_import/#cgmes-import-config"),
  "cgmes_import.infer_base_voltages" => (label = "Infer missing base voltages", hint = "Reconstruct missing nominal voltages from SV voltages and transformer ratings when the delivery has no BaseVoltage catalog. Pair with an unchecked boundary set.", doc = "cgmes_import/#cgmes-import-config"),
  "power_flow.rescue" => (label = "Rescue ladder for failed AC solves", hint = "After a non-converged AC solve, retry from the original start through a fixed ladder: alternate start, autodamp, DC seed, settled Q-limits. First success wins.", doc = "powerflow_configuration/#pf-solver-core"),
  "runtime.parallel.enabled" => (label = "Parallel execution of independent work items", hint = "Use Julia threads for independent work items: island solves, short-circuit sweeps, contingency batches. Off forces every site onto the serial path.", doc = "performance_profiling/#perf-runtime"),
  "power_flow.dc.fallback" => (label = "Standalone-DC fallback", hint = "When AC and the rescue ladder fail, keep a standalone DC result: angles and branch P flows at 1 pu, no reactive results. The AC status stays non-converged.", doc = "powerflow_configuration/#pf-solver-core"),
  "model.auto_profile" => (label = "MATPOWER auto-profile", hint = "MATPOWER pre-run profile: off disables it, recommend logs import-convention recommendations, apply changes only safe conventions with clear evidence.", doc = "matpower/#matpower-options"),
  "matpower_import.ratio" => (label = "Transformer ratio convention", hint = "Interpretation of the branch ratio column: normal is the standard MATPOWER convention, reciprocal inverts it. Change it only for a known alternate convention.", doc = "matpower/#matpower-options"),
  "matpower_import.apply_bus_names" => (label = "Apply bus names", hint = "Use the case file's mpc.bus_name entries as bus names in results and logs instead of numeric BUS_I ids. Requires a bus_name block matching the bus count.", doc = "matpower/#matpower-options"),
  "matpower_import.shift_sign" => (label = "Phase-shift sign", hint = "Sign convention of the phase-shift column, 1 or -1. Flip it only to align with another tool's convention.", doc = "matpower/#matpower-options"),
  "matpower_import.shift_unit" => (label = "Phase-shift unit", hint = "Unit of the phase-shift column, degrees or radians. A wrong declaration shifts every phase shifter.", doc = "matpower/#matpower-options"),
  "model.bus_shunt_model" => (label = "Bus-shunt model", hint = "Interpretation of bus shunts: constant admittance (default) or a voltage-dependent injection. Change it only with residual evidence.", doc = "matpower/#matpower-options"),
  "matpower_import.pv_voltage_source" => (label = "PV voltage source", hint = "Source of the PV voltage setpoint: the generator VG (standard MATPOWER) or the bus VM; auto and strict_check use VG and warn when VG and VM differ.", doc = "matpower/#matpower-options"),
  "matpower_import.compare_voltage_reference" => (label = "Voltage reference comparison", hint = "Voltage reference used when results are compared with the case data: bus VM, generator VG, the imported setpoint, or hybrid when VM and VG disagree.", doc = "matpower/#matpower-options"),
  "model.tap_changer_model" => (label = "Tap-changer model", hint = "Tap-changer model applied to all transformers after import: ideal, or impedance correction that scales R and X with the tap position.", doc = "matpower/#matpower-options"),
  "matpower_export.write_solution" => (label = "Write solution into MATPOWER export", hint = "Write the solved bus VM/VA state and the branch flow columns into the MATPOWER export. Off exports a pure model file with flat voltages.", doc = "matpower/#matpower-options"),
  "output.logfile_results" => (label = "Logfile output mode", hint = "Detail of the solved result tables in run.log: off, compact, classic (result report plus timing summary) or full (adds the effective configuration).", doc = "performance_profiling/#perf-output"),
  "benchmark.enabled" => (label = "Enable benchmark measurements", hint = "Measure repeated solves and report their median instead of one timing. Bounded by the sample count and the time budget.", doc = "performance_profiling/#perf-benchmark"),
  "benchmark.samples" => (label = "Benchmark samples (max. repeated measurements)", hint = "Maximum number of repeated benchmark measurements per method. Collection stops earlier when the time budget is used up first.", doc = "performance_profiling/#perf-benchmark"),
  "benchmark.seconds" => (label = "Benchmark max. time budget [s]", hint = "Maximum time budget of the benchmark in seconds. Not a solver timeout: a running sample is never interrupted.", doc = "performance_profiling/#perf-benchmark"),
  "webui.performance_timing" => (label = "Performance timing", hint = "Write performance.log with the phases of one request (parsing, case loading, solve, artifacts); full adds internal profile entries, off writes nothing.", doc = "webui/#webui-output-modes"),
  "webui.detailed_result_csv" => (label = "Bus/branch CSV files", hint = "Write bus_voltages_complex.csv and branch_flows.csv with per-bus voltages and per-branch flows. Off by default because large networks produce large files.", doc = "webui/#webui-output-modes"),
  "webui.detailed_result_csv_format" => (label = "CSV format (every CSV file of a run)", hint = "Delimiter and decimal separator of every CSV a run writes: technical (comma, point), excel_de (semicolon, decimal comma) or excel_us. A machine-wide setting.", doc = "webui/#webui-output-modes"),
  "webui.export_cgmes" => (label = "CGMES export artifact", hint = "Write the case as one re-importable CGMES delivery (EQ, TP, SSH, SV in a ZIP) into the run's artifacts, for every case format and also on non-converged runs.", doc = "cgmes_export/#cgmes-export-webui"),
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
  "power_flow_apslf_convergence_radius" => "power_flow.apslf.convergence_radius",
  "power_flow_flatstart" => "power_flow.flatstart",
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

"""Resolve an allowlisted Web UI help topic to its label, hint and documentation link."""
resolve_webui_help_topic(topic::AbstractString) = get(WEBUI_HELP_TOPICS, String(topic), nothing)

# The documentation base URL the help icons open. Set from the running
# server's configuration (`webui.docs_base_url`) on every request; the
# environment variable wins for headless renders and tests. A link only,
# nothing is ever fetched.
const _WEBUI_DOCS_BASE_URL = Ref{String}("https://welthulk.github.io/Sparlectra.jl/")

function _webui_docs_base_url()::String
  base = get(ENV, "SPARLECTRA_WEBUI_DOCS_BASE_URL", _WEBUI_DOCS_BASE_URL[])
  return endswith(base, "/") ? base : base * "/"
end

# Refresh the base URL from a runtime's configuration file; a file that
# cannot be loaded keeps the previous value (the link is not worth failing a
# page for).
function _webui_refresh_docs_base_url!(runtime)
  runtime === nothing && return nothing
  config_file = String(something(getproperty(runtime, :config_file), ""))
  isempty(config_file) && return nothing
  isfile(config_file) || return nothing
  try
    _WEBUI_DOCS_BASE_URL[] = load_sparlectra_config(config_file).webui.docs_base_url
  catch err
    # a broken configuration file is reported by the pages that load it;
    # the help link falls back to the packaged default
    @debug "webui.docs_base_url not read" config_file exception = err
  end
  return nothing
end

"""Absolute documentation URL of a help topic, or "" when the topic has no `doc` target."""
function webui_help_doc_url(topic::AbstractString)::String
  metadata = resolve_webui_help_topic(topic)
  (metadata === nothing || isempty(metadata.doc)) && return ""
  return _webui_docs_base_url() * String(metadata.doc)
end

"""
    webui_help_page_url(topic) -> String

The in-app help page of a topic (`/help/<topic>`), or "" when the topic has
no `doc` target. The page shows the hint, the documentation section shipped
with the application (`WEBUI_HELP_EXCERPTS`) and the link to the same
section online, so the help does not depend on the published site carrying
the anchors of the version that is running.
"""
function webui_help_page_url(topic::AbstractString)::String
  metadata = resolve_webui_help_topic(topic)
  (metadata === nothing || isempty(metadata.doc)) && return ""
  return "/help/" * _webui_urlencode(String(topic))
end

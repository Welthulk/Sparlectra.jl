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

const _WEBUI_DOCS_ROOT = normpath(joinpath(@__DIR__, "..", "..", "docs", "src"))

const WEBUI_HELP_TOPICS = Dict(
  "webui.casefile" => (label = "MATPOWER case file", page = "webui", heading = "PowerFlow input paths", selector = "`webui.casefile`"),
  "webui.config_file" => (label = "Configuration template file", page = "webui", heading = "PowerFlow input paths", selector = "`webui.config_file`"),
  "webui.import_case_files" => (label = "Import case files", page = "webui", heading = "Importing case files through the Web UI", selector = ""),
  "webui.case_format" => (label = "Case input format", page = "webui", heading = "Starting a PowerFlow run", selector = ""),
  "webui.for002_reference_file" => (label = "Optional FOR002 reference file", page = "webui", heading = "Starting a PowerFlow run", selector = ""),
  "webui.dtf_outage_selection" => (label = "Selected DTF outage labels/indices", page = "webui", heading = "Starting a PowerFlow run", selector = ""),
  "webui.config_maintenance" => (label = "Configuration maintenance", page = "webui", heading = "Configuration check and refresh", selector = ""),
  "webui.ignore_webui_settings" => (label = "Ignore Web UI settings and use configuration defaults", page = "webui", heading = "Configuration precedence and artifact downloads", selector = ""),
  "power_flow.tol" => (label = "PowerFlow tolerance", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.tol`"),
  "power_flow.max_iter" => (label = "Maximum iterations", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.max_iter`"),
  "power_flow.autodamp" => (label = "Autodamping enabled", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.autodamp`"),
  "power_flow.autodamp_min" => (label = "Autodamping minimum", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.autodamp_min`"),
  "power_flow.qlimits.enabled" => (label = "Q-limit handling enabled", page = "powerflow_configuration", heading = "Q-limit options and guard", selector = "`power_flow.qlimits.enabled`"),
  "power_flow.qlimits.enforcement_mode" => (label = "Q-limit enforcement mode", page = "powerflow_configuration", heading = "Q-limit options and guard", selector = "`power_flow.qlimits.enforcement_mode`"),
  "power_flow.solver" => (label = "Solver", page = "powerflow_configuration", heading = "Solver selection (rectangular vs. APSLF)", selector = "`power_flow.solver`"),
  "power_flow.apslf.order" => (label = "APSLF highest coefficient (order)", page = "powerflow_configuration", heading = "Solver selection (rectangular vs. APSLF)", selector = "`power_flow.apslf.order`"),
  "power_flow.apslf.use_pade" => (label = "APSLF Padé evaluation", page = "powerflow_configuration", heading = "Solver selection (rectangular vs. APSLF)", selector = "`power_flow.apslf.use_pade`"),
  "power_flow.apslf.nr_polish" => (label = "APSLF NR polish", page = "powerflow_configuration", heading = "Solver selection (rectangular vs. APSLF)", selector = "`power_flow.apslf.nr_polish`"),
  "power_flow.apslf_start.enabled" => (label = "Use APSLF start values", page = "powerflow_configuration", heading = "Solver selection (rectangular vs. APSLF)", selector = "`power_flow.apslf_start.enabled`"),
  "power_flow.apslf_start.order" => (label = "APSLF start highest coefficient (order)", page = "powerflow_configuration", heading = "Solver selection (rectangular vs. APSLF)", selector = "`power_flow.apslf_start.order`"),
  "power_flow.wrong_branch_detection" => (label = "Wrong-branch detection", page = "configuration", heading = "Wrong-branch detection semantics (rectangular PF)", selector = ""),
  "power_flow.start_mode.angle_mode" => (label = "Start angle mode", page = "powerflow_configuration", heading = "Start mode options", selector = "`power_flow.start_mode.angle_mode`"),
  "power_flow.start_mode.voltage_mode" => (label = "Start voltage mode", page = "powerflow_configuration", heading = "Start mode options", selector = "`power_flow.start_mode.voltage_mode`"),
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
  "power_flow.merit.enabled" => (label = "Enable Armijo merit-function line search", page = "powerflow_configuration", heading = "Merit-function line search options", selector = "`power_flow.merit.enabled`"),
  "power_flow.merit.armijo_c1" => (label = "Armijo sufficient-decrease constant", page = "powerflow_configuration", heading = "Merit-function line search options", selector = "`power_flow.merit.armijo_c1`"),
  "power_flow.merit.fallback_max_mismatch" => (label = "Merit fallback behavior", page = "powerflow_configuration", heading = "Merit-function line search options", selector = "`power_flow.merit.fallback_max_mismatch`"),
  "power_flow.trust_region.enabled" => (label = "Enable trust-region step control", page = "powerflow_configuration", heading = "Trust-region step control options", selector = "`power_flow.trust_region.enabled`"),
  "power_flow.trust_region.initial_radius" => (label = "Initial trust-region radius", page = "powerflow_configuration", heading = "Trust-region step control options", selector = "`power_flow.trust_region.initial_radius`"),
  "power_flow.trust_region.eta_accept" => (label = "Trust-region acceptance ratio (eta)", page = "powerflow_configuration", heading = "Trust-region step control options", selector = "`power_flow.trust_region.eta_accept`"),
  "power_flow.trust_region.step_mode" => (label = "Trust-region step mode", page = "powerflow_configuration", heading = "Trust-region step control options", selector = "`power_flow.trust_region.step_mode`"),
  "matpower_import.auto_profile" => (label = "MATPOWER auto-profile", page = "matpower_import", heading = "Option reference", selector = "`matpower_import.auto_profile`"),
  "matpower_import.ratio" => (label = "Transformer ratio convention", page = "matpower_import", heading = "Option reference", selector = "`matpower_import.ratio`"),
  "matpower_import.shift_sign" => (label = "Phase-shift sign", page = "matpower_import", heading = "Option reference", selector = "`matpower_import.shift_sign`"),
  "matpower_import.shift_unit" => (label = "Phase-shift unit", page = "matpower_import", heading = "Option reference", selector = "`matpower_import.shift_unit`"),
  "matpower_import.bus_shunt_model" => (label = "Bus-shunt model", page = "matpower_import", heading = "Option reference", selector = "`matpower_import.bus_shunt_model`"),
  "matpower_import.pv_voltage_source" => (label = "PV voltage source", page = "matpower_import", heading = "Option reference", selector = "`matpower_import.pv_voltage_source`"),
  "matpower_import.compare_voltage_reference" => (label = "Voltage reference comparison", page = "matpower_import", heading = "Option reference", selector = "`matpower_import.compare_voltage_reference`"),
  "transformer.tap_changer_model" => (label = "Tap-changer model", page = "matpower_import", heading = "Option reference", selector = "`transformer.tap_changer_model`"),
  "matpower_export.write_solution" => (label = "Write solution into MATPOWER export", page = "matpower_import", heading = "Option reference", selector = "`matpower_export.write_solution`"),
  "output.logfile_results" => (label = "Logfile output mode", page = "performance_profiling", heading = "Output configuration", selector = "`output.logfile_results`"),
  "benchmark.enabled" => (label = "Enable benchmark measurements", page = "performance_profiling", heading = "Benchmark configuration", selector = "`benchmark.enabled`"),
  "benchmark.samples" => (label = "Benchmark samples (max. repeated measurements)", page = "performance_profiling", heading = "Benchmark configuration", selector = "`benchmark.samples`"),
  "benchmark.seconds" => (label = "Benchmark max. time budget [s]", page = "performance_profiling", heading = "Benchmark configuration", selector = "`benchmark.seconds`"),
  "webui.performance_timing" => (label = "Performance timing", page = "webui", heading = "Run artifacts and output modes", selector = ""),
  "webui.detailed_result_csv" => (label = "Detailed result CSV export", page = "webui", heading = "Run artifacts and output modes", selector = ""),
  "webui.detailed_result_csv_format" => (label = "Detailed CSV format", page = "webui", heading = "Run artifacts and output modes", selector = ""),
)

const WEBUI_FORM_HELP_TOPICS = Dict(
  "casefile" => "webui.casefile",
  "casefile_manual" => "webui.casefile",
  "config_file" => "webui.config_file",
  "casefiles" => "webui.import_case_files",
  "case_format" => "webui.case_format",
  "for002_reference_file" => "webui.for002_reference_file",
  "dtf_outage_selection" => "webui.dtf_outage_selection",
  "config_maintenance" => "webui.config_maintenance",
  "ignore_webui_settings" => "webui.ignore_webui_settings",
  "power_flow_tol" => "power_flow.tol",
  "power_flow_max_iter" => "power_flow.max_iter",
  "power_flow_autodamp" => "power_flow.autodamp",
  "power_flow_autodamp_min" => "power_flow.autodamp_min",
  "power_flow_qlimits_enabled" => "power_flow.qlimits.enabled",
  "power_flow_qlimits_enforcement_mode" => "power_flow.qlimits.enforcement_mode",
  "power_flow_solver" => "power_flow.solver",
  "power_flow_apslf_order" => "power_flow.apslf.order",
  "power_flow_apslf_use_pade" => "power_flow.apslf.use_pade",
  "power_flow_apslf_nr_polish" => "power_flow.apslf.nr_polish",
  "power_flow_apslf_start_enabled" => "power_flow.apslf_start.enabled",
  "power_flow_apslf_start_order" => "power_flow.apslf_start.order",
  "power_flow_wrong_branch_detection" => "power_flow.wrong_branch_detection",
  "power_flow_start_angle_mode" => "power_flow.start_mode.angle_mode",
  "power_flow_start_voltage_mode" => "power_flow.start_mode.voltage_mode",
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
  "power_flow_merit_enabled" => "power_flow.merit.enabled",
  "power_flow_merit_armijo_c1" => "power_flow.merit.armijo_c1",
  "power_flow_merit_fallback_max_mismatch" => "power_flow.merit.fallback_max_mismatch",
  "power_flow_trust_region_enabled" => "power_flow.trust_region.enabled",
  "power_flow_trust_region_initial_radius" => "power_flow.trust_region.initial_radius",
  "power_flow_trust_region_eta_accept" => "power_flow.trust_region.eta_accept",
  "power_flow_trust_region_step_mode" => "power_flow.trust_region.step_mode",
  "matpower_import_auto_profile" => "matpower_import.auto_profile",
  "matpower_import_ratio" => "matpower_import.ratio",
  "matpower_import_shift_sign" => "matpower_import.shift_sign",
  "matpower_import_shift_unit" => "matpower_import.shift_unit",
  "matpower_import_bus_shunt_model" => "matpower_import.bus_shunt_model",
  "matpower_import_pv_voltage_source" => "matpower_import.pv_voltage_source",
  "matpower_import_compare_voltage_reference" => "matpower_import.compare_voltage_reference",
  "transformer_tap_changer_model" => "transformer.tap_changer_model",
  "matpower_export_write_solution" => "matpower_export.write_solution",
  "output_logfile_results" => "output.logfile_results",
  "benchmark_enabled" => "benchmark.enabled",
  "benchmark_samples" => "benchmark.samples",
  "benchmark_seconds" => "benchmark.seconds",
  "performance_timing" => "webui.performance_timing",
  "detailed_result_csv" => "webui.detailed_result_csv",
  "detailed_result_csv_format" => "webui.detailed_result_csv_format",
)

const WEBUI_HELP_EXCERPT_OVERRIDES = Dict(
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
  "matpower_import" => (title = "MATPOWER Import", file = "matpower_import.md"),
  "matpower_case_matrix" => (title = "MATPOWER Case Diagnostics Matrix", file = "sparlectra_matpower_case_matrix.md"),
  "webui" => (title = "Local PowerFlow Web UI", file = "webui.md"),
  "feature_matrix" => (title = "Feature Matrix", file = "feature_matrix.md"),
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
  text = lowercase(replace(text, "&amp;" => "and", "&quot;" => "", "&#39;" => ""))
  return strip(replace(text, r"[^a-z0-9]+" => "-"), '-')
end

function _webui_rewritten_doc_href(target::AbstractString; current_page::Union{Nothing,String} = nothing)::Union{String,Nothing}
  href = String(target)
  startswith(href, "https://") && return href
  startswith(href, "http://") && return href
  startswith(href, "#") && return current_page === nothing ? nothing : href
  (startswith(href, "/") || occursin('\\', href) || occursin(':', href)) && return nothing

  relative = startswith(href, "./") ? href[3:end] : href
  matched = match(r"^([A-Za-z0-9_-]+)\.md(#[A-Za-z0-9._~:%-]+)?$", relative)
  matched === nothing && return nothing
  page = matched.captures[1]
  metadata = resolve_webui_doc_page(page)
  metadata === nothing && return nothing
  metadata.file == "$(page).md" || return nothing
  fragment = something(matched.captures[2], "")
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

"""Render trusted repository Markdown as HTML using Julia's Markdown standard library."""
function render_webui_markdown(markdown_text::AbstractString; current_page::Union{Nothing,String} = nothing)::String
  io = IOBuffer()
  show(io, MIME"text/html"(), Markdown.parse(String(markdown_text)))
  return rewrite_webui_doc_links(String(take!(io)); current_page = current_page)
end

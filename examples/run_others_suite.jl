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

# Date: 2026-07-26
# file: examples/run_others_suite.jl
# purpose: suite runner that executes the remaining example programs (transformer/tap demos, exports, network analysis, diagnostics) in fresh subprocesses and reports a summary

include(joinpath(@__DIR__, "internal", "example_suite_runner.jl"))

# Not registered here: the DTF validation examples
# (dtf_validation_report.jl, for002_matpower_metadata_validation.jl) are
# covered by examples/run_val_dtf_suite.jl; library files
# (experimental/qlimit_large_case_comparison.jl, internal/*) are not
# standalone programs.
const SUITE_SPECS = ExampleSpec[
  ExampleSpec(name = "3wt_phase_taps", file = "others/exp_3wt_phase_taps.jl", purpose = "3WT in three tap configurations (OLTC, PST/Schraegregler, combined) solved with runpf!"),
  ExampleSpec(name = "transformer_tap_changer_model", file = "others/exp_transformer_tap_changer_model.jl", purpose = "compares :ideal vs :impedance_correction tap_changer_model on an off-nominal-tap transformer"),
  ExampleSpec(name = "transformer_loss_extension", file = "others/exp_transformer_loss_extension.jl", purpose = "MATPOWER transformer-loss extension export/reimport round trip"),
  ExampleSpec(name = "tap_control_demo_grid", file = "others/tap_control_demo_grid.jl", purpose = "three-controller demo (OLTC + PST + Schraegregler) via run_sparlectra(net=...)"),
  ExampleSpec(name = "tap_control_schraeg_two_controllers", file = "others/tap_control_schraeg_two_controllers.jl", purpose = "split Schraegregelung: independent voltage (ratio tap) and active-power (phase tap) controllers on one transformer"),
  ExampleSpec(name = "pst_reactance_coupling", file = "others/exp_pst_reactance_coupling.jl", purpose = "PST outer-loop control with tap-dependent series reactance X(alpha) vs. a static-reactance run"),
  ExampleSpec(name = "machine_remote_voltage_control", file = "others/machine_remote_voltage_control.jl", purpose = "remote voltage control via machine reactive power (MachineVoltageControl + run_control!), incl. the honest at_limit outcome"),
  ExampleSpec(name = "svc_shunt_voltage_control", file = "others/exp_svc_shunt_voltage_control.jl", purpose = "SVC-style variable-shunt voltage control (in-range regulation, honest at_limit clamping) plus the generic controllable-element view"),
  ExampleSpec(name = "tcsc_series_reactance_control", file = "others/exp_tcsc_series_reactance_control.jl", purpose = "TCSC series-reactance control steering a loop-network flow split onto a branch target (issue #297), incl. honest at_limit and the controllable-element view"),
  ExampleSpec(name = "hvdc_b2b_pairing", file = "others/exp_hvdc_b2b_pairing.jl", purpose = "back-to-back HVDC pairing controller: Stage-0 snapshot vs steerable transfer on a two-island net (issue #297 Draft B)"),
  ExampleSpec(name = "auto_slack_selection", file = "others/exp_auto_slack_selection.jl", purpose = "automatic slack selection (power_flow.auto_slack / ensureSlack!) on a case whose data registers no voltage reference"),
  ExampleSpec(name = "ac_rescue_dc_fallback", file = "others/exp_ac_rescue_dc_fallback.jl", purpose = "non-convergence handling: power_flow.rescue strategy ladder plus power_flow.dc.fallback standalone-DC result"),
  ExampleSpec(name = "cgmes_import_analysis", file = "others/exp_cgmes_import_analysis.jl", purpose = "analyzeCGMES report naming the missing declared dependency of an incomplete CGMES delivery"),
  ExampleSpec(name = "cgmes_infer_base_voltages", file = "others/exp_cgmes_infer_base_voltages.jl", purpose = "cgmes_import.infer_base_voltages: reconstruct missing nominal voltages from the SV state and solve"),
  ExampleSpec(name = "export_solution", file = "others/export_solution.jl", args = ["case9.m"], purpose = "exports a solver-agnostic PFModel/PFSolution for case9"),
  ExampleSpec(name = "network_analyzer", file = "others/network_analyzer.jl", purpose = "topology analysis of a small network before and after removing a branch"),
  ExampleSpec(name = "using_links", file = "others/using_links.jl", purpose = "busbar coupler modeled as a bus link, demonstrating open/close link behavior"),
  ExampleSpec(name = "diagnose_self_check", file = "others/exp_diagnose_self_check.jl", purpose = "run_fixed_reference_self_check and the narrative diagnose.log report"),
  ExampleSpec(name = "cgmes_export_demo", file = joinpath("experimental", "cgmes_export_demo.jl"), optional = true, purpose = "CGMES 2.4.15 export (writeCGMESFiles) on a small net"),
]

const SUITE_NOTES = [
  "DTF validation examples are covered by `examples/run_val_dtf_suite.jl` and are not part of this suite.",
  "Short-circuit examples are covered by `examples/run_short_circuit_suite.jl` and are not part of this suite.",
]

# SPARLECTRA_EXAMPLE_SUITE_NO_MAIN=1 is a test-only escape hatch: it lets the
# test suite include this file to inspect the registry without running it.
if get(ENV, "SPARLECTRA_EXAMPLE_SUITE_NO_MAIN", "0") != "1"
  run_example(run_example_suite, "others_suite", "Others example suite (transformer, export, analysis)", SUITE_SPECS; notes = SUITE_NOTES)
end
nothing

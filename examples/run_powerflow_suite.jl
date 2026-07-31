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
# file: examples/run_powerflow_suite.jl
# purpose: suite runner that executes the power-flow example programs (AC/DC solves, MATPOWER import, Q-limit and control demos) in fresh subprocesses and reports a summary

include(joinpath(@__DIR__, "internal", "example_suite_runner.jl"))

# Not registered: current_iteration_start.jl is a legacy include-alias of
# exp_current_iteration_start.jl and would run the same demo twice.
const SUITE_SPECS = ExampleSpec[
  ExampleSpec(name = "dc_powerflow", file = "powerflow/exp_dc_powerflow.jl", purpose = "standalone DC power flow, optionally seeding the AC Newton-Raphson solve"),
  ExampleSpec(name = "programmatic_api", file = "powerflow/exp_programmatic_api.jl", purpose = "runs one MATPOWER case through the GUI-ready run_sparlectra_api contract"),
  ExampleSpec(name = "powerflow_service", file = "powerflow/exp_powerflow_service.jl", purpose = "local PowerFlow service run with result lookup by run ID (no HTTP server)"),
  ExampleSpec(name = "current_iteration_start", file = "powerflow/exp_current_iteration_start.jl", purpose = "guarded current-iteration start pre-solve via API configuration overrides"),
  ExampleSpec(name = "q_limit_voltage_adjustment", file = "powerflow/example_q_limit_voltage_adjustment.jl", purpose = "compares Q-limit :adjust_vset outcomes across three PV->PQ scenarios"),
  ExampleSpec(name = "voltage_dependent_control", file = "powerflow/example_voltage_dependent_control_rectangular.jl", args = ["--no-plot"], purpose = "voltage-dependent control (Q(U)/P(U)) demo without plotting"),
  # No requires_config gate: the example runs on the resolved configuration
  # chain, whose shipped default (src/configuration.yaml.example, case14.m)
  # always exists — a user configuration.yaml is optional, not required.
  ExampleSpec(name = "matpower_import", file = "powerflow/matpower_import.jl", purpose = "CLI MATPOWER import via run_matpower_case using the resolved configuration"),
  ExampleSpec(name = "matpower_import_multi_config", file = "powerflow/matpower_import_multi_config.jl", purpose = "compares one MATPOWER case across configuration files"),
  ExampleSpec(name = "configured_matpower_cases", file = "powerflow/exp_configured_matpower_cases.jl", purpose = "runs ordered matpower_import.cases config entries via run_sparlectra_cases"),
  ExampleSpec(name = "mc_probabilistic_powerflow", file = "powerflow/mc_probabilistic_powerflow.jl", purpose = "Monte-Carlo probabilistic power flow on case14"),
  ExampleSpec(name = "synthetic_tiled_grid_pf_perf", file = "powerflow/exp_synthetic_tiled_grid_pf_perf.jl", heavy = true, timeout_s = 1800, purpose = "synthetic tiled-grid power-flow performance benchmark"),
  ExampleSpec(name = "qlimit_large_case_mode_comparison", file = "powerflow/qlimit_large_case_mode_comparison.jl", heavy = true, timeout_s = 1800, purpose = "Q-limit start-profile/enforcement-mode comparison on very large MATPOWER cases"),
  ExampleSpec(name = "apslf_demo", file = "powerflow/apslf_demo.jl", optional = true, requires_package = "AnalyticLoadFlow", purpose = "compares the internal rectangular NR solver against the APSLF solver on case30"),
]

# SPARLECTRA_EXAMPLE_SUITE_NO_MAIN=1 is a test-only escape hatch: it lets the
# test suite include this file to inspect the registry without running it.
if get(ENV, "SPARLECTRA_EXAMPLE_SUITE_NO_MAIN", "0") != "1"
  run_example(run_example_suite, "powerflow_suite", "Power-flow example suite", SUITE_SPECS)
end
nothing

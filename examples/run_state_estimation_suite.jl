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
# file: examples/run_state_estimation_suite.jl
# purpose: suite runner that executes the state-estimation example programs (WLS, observability, diagnostics, Monte-Carlo study) in fresh subprocesses and reports a summary

include(joinpath(@__DIR__, "internal", "example_suite_runner.jl"))

# state_estimation_observability.jl internally includes
# state_estimation_wls.jl; that duplication is contained in its own
# subprocess and log.
const SUITE_SPECS = ExampleSpec[
  ExampleSpec(name = "wls", file = "state_estimation/state_estimation_wls.jl", purpose = "classical WLS state-estimation workflow with synthetic measurements"),
  ExampleSpec(name = "manual_measurements", file = "state_estimation/state_estimation_manual_measurements.jl", purpose = "measurement set built via the manual measurement helpers (addVmMeasurement!, ...)"),
  ExampleSpec(name = "observability", file = "state_estimation/state_estimation_observability.jl", purpose = "progressively deactivates branch-flow measurements and logs observability-redundancy metrics"),
  ExampleSpec(name = "passive_bus_zib_comparison", file = "state_estimation/state_estimation_passive_bus_zib_comparison.jl", purpose = "WLS state estimation with and without zero-injection (ZIB) measurements"),
  ExampleSpec(name = "pmu_angles", file = "state_estimation/state_estimation_pmu_angles.jl", purpose = "PMU voltage-angle measurements with the reference-angle offset state alpha (aligned vs. shifted PMU time base)"),
  ExampleSpec(name = "diagnostics", file = "state_estimation/usage_state_estimation_diagnostics.jl", purpose = "bad-data diagnostics: inject a bad measurement, validate, deactivate and rerun"),
  ExampleSpec(name = "h_matrix_observability", file = "state_estimation/h_matrix_observability_demo.jl", purpose = "small measurement Jacobians H and the public observability helpers"),
  ExampleSpec(name = "mc_study", file = "state_estimation/mc_state_estimation_study.jl", timeout_s = 1200, purpose = "Monte-Carlo WLS state-estimation error study on the 7-bus workshop net"),
]

# SPARLECTRA_EXAMPLE_SUITE_NO_MAIN=1 is a test-only escape hatch: it lets the
# test suite include this file to inspect the registry without running it.
if get(ENV, "SPARLECTRA_EXAMPLE_SUITE_NO_MAIN", "0") != "1"
  run_example(run_example_suite, "state_estimation_suite", "State-estimation example suite", SUITE_SPECS)
end
nothing

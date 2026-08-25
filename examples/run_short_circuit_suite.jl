# Copyright 2023-2026 Udo Schmitz
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

# Date: 2026-07-30
# file: examples/run_short_circuit_suite.jl
# purpose: suite runner for the IEC 60909 short-circuit examples (issue #277) in fresh subprocesses with a summary report

include(joinpath(@__DIR__, "internal", "example_suite_runner.jl"))

const SUITE_SPECS = ExampleSpec[
  ExampleSpec(name = "short_circuit", file = "others/exp_short_circuit.jl", purpose = "runShortCircuit! on a hand-built feeder+machine net — Ik'' max/min and the safety flag on defaulted data"),
  ExampleSpec(name = "short_circuit_reference", file = "others/exp_short_circuit_reference.jl", purpose = "PASS/FAIL check against the analytic IEC 60909-0 reference values"),
  ExampleSpec(name = "short_circuit_cgmes", file = "others/exp_short_circuit_cgmes.jl", purpose = "runShortCircuit! on the ENTSO-E MicroGrid BE delivery (local test-set cache; explains the fetch and skips cleanly when absent)"),
]

const SUITE_NOTES = ["Unbalanced faults (single line-to-earth etc.) are not supported yet — see docs/src/short_circuit.md."]

# SPARLECTRA_EXAMPLE_SUITE_NO_MAIN=1 is a test-only escape hatch: it lets the
# test suite include this file to inspect the registry without running it.
if get(ENV, "SPARLECTRA_EXAMPLE_SUITE_NO_MAIN", "0") != "1"
  run_example(run_example_suite, "short_circuit_suite", "Short-circuit example suite (IEC 60909)", SUITE_SPECS; notes = SUITE_NOTES)
end
nothing

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

# Date: 2026-08-21
# file: examples/run_parallel_suite.jl
# purpose: suite runner for the parallel-vs-serial demos of the multi-core
#          work: every example solves the same workload serially and on
#          Julia threads, prints both wall clocks side by side, and FAILS
#          if the results differ. Subprocesses get JULIA_NUM_THREADS=auto,
#          so a plain `julia --project=. examples/run_parallel_suite.jl`
#          already shows the speedups.

# The suite runner spawns each example via Base.julia_cmd(), which does NOT
# carry a --threads flag; the environment variable does reach the
# subprocesses. Respect an explicit user setting, default to auto.
haskey(ENV, "JULIA_NUM_THREADS") || (ENV["JULIA_NUM_THREADS"] = "auto")

include(joinpath(@__DIR__, "internal", "example_suite_runner.jl"))

const SUITE_SPECS = ExampleSpec[
  ExampleSpec(name = "parallel_islands", file = "powerflow/exp_parallel_islands.jl", purpose = "power_flow.islands.mode solve_parallel vs solve_independent on an 8-island net, bitwise-identical voltages"),
  ExampleSpec(name = "parallel_sc_sweep", file = "others/exp_parallel_sc_sweep.jl", purpose = "IEC 60909-0 all-bus sweep serial vs threaded chunks (8000 fault locations), row-identical results"),
  ExampleSpec(name = "contingency_n1", file = "others/exp_contingency_n1.jl", purpose = "full branch N-1 on case1354pegase serial vs parallel plus top-10 worst contingencies (skips politely without the cached case)"),
]

const SUITE_NOTES = [
  "Speedups depend on the machine's core count; the subprocesses run with JULIA_NUM_THREADS=auto (an explicit environment setting wins).",
  "Every example asserts serial/parallel result identity and fails loudly on a mismatch.",
]

# SPARLECTRA_EXAMPLE_SUITE_NO_MAIN=1 is a test-only escape hatch: it lets the
# test suite include this file to inspect the registry without running it.
if get(ENV, "SPARLECTRA_EXAMPLE_SUITE_NO_MAIN", "0") != "1"
  run_example(run_example_suite, "parallel_suite", "Parallel-vs-serial suite (islands, SC sweep, N-1)", SUITE_SPECS; notes = SUITE_NOTES)
end
nothing

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

# Date: 2026-06-19
# file: examples/exp_current_iteration_start.jl
# purpose: demonstrates enabling the guarded current-iteration start pre-solve via API configuration overrides and prints its metadata/artifact status

using Sparlectra

include(joinpath(@__DIR__, "internal", "example_header.jl"))

function main(; casefile::AbstractString = "case5.m", output_dir::AbstractString = joinpath(@__DIR__, "_out", "current_iteration_start"))
  print_example_banner("examples/exp_current_iteration_start.jl", "demonstrates enabling the guarded current-iteration start pre-solve via API configuration overrides and prints its metadata/artifact status")
  case_path = ensure_casefile(casefile)
  result = run_sparlectra_api(
    casefile = case_path,
    config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
    output_dir = output_dir,
    config_overrides = Dict(
      "power_flow.start_current_iteration.enabled" => true,
      "power_flow.start_current_iteration.max_iter" => 10,
      "power_flow.start_current_iteration.tol" => 1.0e-3,
      "power_flow.start_current_iteration.damping" => 0.5,
      "power_flow.start_current_iteration.accept_only_if_improved" => true,
      "power_flow.start_current_iteration.min_improvement_factor" => 0.98,
      "power_flow.start_current_iteration.vm_min_pu" => 0.5,
      "power_flow.start_current_iteration.vm_max_pu" => 1.5,
      "power_flow.start_current_iteration.max_angle_step_deg" => 30.0,
      "power_flow.start_current_iteration.only_for_large_cases" => false,
      "output.logfile_results" => "compact",
      "benchmark.enabled" => false,
    ),
  )

  println("run ID: ", result.run_id)
  println("status: ", result.status)
  println("current iteration enabled: ", get(result.metadata, "current_iteration_enabled", false))
  println("current iteration attempted: ", get(result.metadata, "current_iteration_attempted", false))
  println("current iteration accepted: ", get(result.metadata, "current_iteration_accepted", false))
  println("current iteration iterations: ", get(result.metadata, "current_iteration_iterations", 0))
  println("current iteration initial mismatch: ", get(result.metadata, "current_iteration_initial_mismatch", "n/a"))
  println("current iteration final mismatch: ", get(result.metadata, "current_iteration_final_mismatch", "n/a"))
  println("current iteration reason: ", get(result.metadata, "current_iteration_reason", ""))
  for artifact in result.artifacts
    artifact.name == "current_iteration_start.log" && println("diagnostic artifact: ", artifact.path)
  end
  return result
end

run_example(main)

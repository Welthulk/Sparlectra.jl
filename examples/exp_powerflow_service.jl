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

# Date: 2026-06-10
# file: examples/exp_powerflow_service.jl
# purpose: starts a local PowerFlow service run, looks up its serialized result by run ID, and lists its artifacts without an HTTP server

using Sparlectra

include(joinpath(@__DIR__, "internal", "example_header.jl"))

function main(;
  casefile::AbstractString = "case5.m",
  output_root::AbstractString = joinpath(@__DIR__, "_out", "powerflow_service"),
)
  print_example_banner("examples/exp_powerflow_service.jl", "starts a local PowerFlow service run, looks up its serialized result by run ID, and lists its artifacts without an HTTP server")
  request = Dict(
    "casefile" => ensure_casefile(casefile),
    "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
    "output_root" => output_root,
    "config_overrides" => Dict(
      "power_flow.tol" => 1.0e-8,
      "power_flow.max_iter" => 80,
      "benchmark.enabled" => false,
    ),
  )

  result = start_powerflow_run(request)
  result["success"] || return result

  run_id = result["run_id"]

  # This can be called in a later Julia process to recover indexed runs.
  Base.invokelatest(refresh_powerflow_run_registry!, output_root)
  runs = list_powerflow_runs(output_root)
  stored_result = get_powerflow_result(run_id)
  artifacts = list_powerflow_artifacts(run_id)

  println("indexed runs: ", length(runs))
  println("run ID: ", run_id)
  println("stored status: ", stored_result["status"])
  for artifact in artifacts
    println(artifact["name"], " (", artifact["mime_type"], ")")
  end
  return result
end

run_example(main)

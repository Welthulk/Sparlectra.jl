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
# file: examples/powerflow/exp_programmatic_api.jl
# purpose: runs one shipped case through the GUI-ready run_sparlectra_api contract and lists explicit artifacts

using Sparlectra
using SparlectraApp
include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

# the shipped 5-bus case (data/scf); examples never download a case
# (another case works too: fetch one yourself, e.g. with
# Sparlectra.ensure_casefile("case14.m"), and use its path here)
function main(; casefile::AbstractString = joinpath(pkgdir(Sparlectra), "data", "scf", "sp_case5.scf.json"), output_dir::AbstractString = joinpath(@__DIR__, "..", "_out", "api_run"))
  print_example_banner("examples/powerflow/exp_programmatic_api.jl", "runs one shipped case through the GUI-ready run_sparlectra_api contract and lists explicit artifacts")
  case_path = casefile
  result = run_sparlectra_api(
    casefile = case_path,
    config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
    output_dir = output_dir,
    config_overrides = Dict(
      "power_flow.tol" => 1.0e-8,
      "power_flow.max_iter" => 80,
      "power_flow.autodamp" => true,
      "output.logfile_results" => "compact",
      "benchmark.enabled" => false,
    ),
  )

  println("run ID: ", result.run_id)
  println("schema version: ", result.schema_version)
  println("status: ", result.status)
  println("solution available: ", result.solution_available)
  for artifact in result.artifacts
    println(artifact.kind, ": ", artifact.path)
  end
  return result
end

run_example(main)

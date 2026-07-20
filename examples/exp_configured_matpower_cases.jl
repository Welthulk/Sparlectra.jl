#!/usr/bin/env julia
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

# Date: 2026-05-30
# file: examples/exp_configured_matpower_cases.jl
# purpose: runs ordered matpower_import.cases entries sequentially through run_sparlectra_cases and prints each case's outcome

using Sparlectra

include(joinpath(@__DIR__, "internal", "example_header.jl"))

function main(args = ARGS)
  print_example_banner("examples/exp_configured_matpower_cases.jl", "runs ordered matpower_import.cases entries sequentially through run_sparlectra_cases and prints each case's outcome")
  config_file = Sparlectra.configuration_path_from_inputs(
    env_var = "SPARLECTRA_CONFIGURATION_YAML",
    fallback_paths = [Sparlectra.USER_SPARLECTRA_CONFIG_PATH, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH],
  )
  cfg = Sparlectra.load_sparlectra_config(config_file)
  results = Sparlectra.run_sparlectra_cases(config = cfg)
  for result in results
    println(result.net.name, ": ", result.outcome)
  end
  return results
end

run_example(main, ARGS)


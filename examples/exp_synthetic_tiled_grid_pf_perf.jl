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

# Date: 2026-05-18
# file: examples/exp_synthetic_tiled_grid_pf_perf.jl
# purpose: runs a synthetic tiled-grid power-flow performance benchmark via run_synthetic_tiled_grid_pf_perf

using Sparlectra

include(joinpath(@__DIR__, "internal", "example_header.jl"))

function main(args = ARGS)
  print_example_banner("examples/exp_synthetic_tiled_grid_pf_perf.jl", "runs a synthetic tiled-grid power-flow performance benchmark via run_synthetic_tiled_grid_pf_perf")
  config_file = Sparlectra.configuration_path_from_inputs(
    env_var = "SPARLECTRA_CONFIGURATION_YAML",
    fallback_paths = [Sparlectra.USER_SPARLECTRA_CONFIG_PATH],
  )
  return Sparlectra.run_synthetic_tiled_grid_pf_perf(; config_file = config_file, args = collect(String, args))
end

run_example(main, ARGS)


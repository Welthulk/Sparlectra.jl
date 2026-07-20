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
# file: examples/example_voltage_dependent_control_rectangular.jl
# purpose: runs the voltage-dependent control (Q(U)/P(U)) demo via run_voltage_dependent_control_demo, optionally plotting the control curve

using Sparlectra

include(joinpath(@__DIR__, "internal", "example_header.jl"))

function main(args = ARGS)
  print_example_banner("examples/example_voltage_dependent_control_rectangular.jl", "runs the voltage-dependent control (Q(U)/P(U)) demo via run_voltage_dependent_control_demo, optionally plotting the control curve")
  config_file = Sparlectra.configuration_path_from_inputs(
    env_var = "SPARLECTRA_CONFIGURATION_YAML",
    fallback_paths = [Sparlectra.USER_SPARLECTRA_CONFIG_PATH],
  )
  plot_curve = any(==("--plot"), args) ? true : (any(==("--no-plot"), args) ? false : nothing)
  return Sparlectra.run_voltage_dependent_control_demo(; config_file = config_file, plot_curve = plot_curve)
end

run_example(main, ARGS)


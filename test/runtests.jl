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

# file: test/runtests.jl

using Sparlectra
using Test
using Logging
using Printf
using LinearAlgebra

# keep logs quiet unless there's a warning or error
global_logger(ConsoleLogger(stderr, Logging.Info))

include("testgrid.jl")
include("testremove.jl")
include("test_solver_interface.jl")
include("test_state_estimation.jl")
include("test_voltage_dependent_control.jl")

@testset "Sparlectra.jl" begin
  run_grid_tests()
  run_remove_tests()
  run_solver_interface_tests()
  run_state_estimation_tests()
  run_voltage_dependent_control_tests()
end

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
global_logger(ConsoleLogger(stderr, Logging.Warn))

# Suppress one-time deprecation warnings in test output. Deprecated solver
# behavior is covered by dedicated tests and does not need to spam full-suite logs.
Sparlectra._warned_full_solver_deprecated[] = true
Sparlectra._warned_classic_solver_deprecated[] = true

include("testgrid.jl")
include("testremove.jl")
include("test_solver_interface.jl")
include("test_state_estimation.jl")
include("test_voltage_dependent_control.jl")
include("test_transformer_phase_shift.jl")
include("test_tap_controller.jl")

function _print_test_progress(step::Int, total::Int, label::String)
  width = 30
  filled = round(Int, width * step / total)
  bar = repeat("█", filled) * repeat("░", width - filled)
  @printf("\rTest progress [%s] %3d%% (%d/%d) %s", bar, round(Int, 100 * step / total), step, total, label)
  flush(stdout)
  if step == total
    println()
  end
end

@testset "Sparlectra.jl" begin
  test_steps = [
    ("grid", run_grid_tests),
    ("remove", run_remove_tests),
    ("solver_interface", run_solver_interface_tests),
    ("state_estimation", run_state_estimation_tests),
    ("voltage_dependent_control", run_voltage_dependent_control_tests),
    ("transformer_phase_shift", run_transformer_phase_shift_tests),
    ("tap_controller", run_tap_controller_tests),
  ]

  total = length(test_steps)
  _print_test_progress(0, total, "starting")
  for (idx, (name, runner)) in enumerate(test_steps)
    _print_test_progress(idx - 1, total, "running $(name)")
    runner()
    _print_test_progress(idx, total, "finished $(name)")
  end
  return nothing
end

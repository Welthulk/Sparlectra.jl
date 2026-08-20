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

# file: test/runtests.jl
# purpose: test suite entry point: selects the fast, extended, or all profile
#          from ARGS/SPARLECTRA_TEST_PROFILE, includes the per-topic test
#          files, and runs the grouped testsets with quiet output capture
using Sparlectra
using Test
using Logging
using Printf
using LinearAlgebra
using SparseArrays

global_logger(ConsoleLogger(stderr, Logging.Warn))

include("test_runner_helpers.jl")

const TEST_PROFILE = selected_test_profile()

function print_test_progress_header(profile::Symbol)
  println("Test framework: ", profile)
end

function print_group_progress(i::Int, total::Int, name::AbstractString) end

function include_fast_tests()
  include("testgrid.jl")
  include("test_piline_g.jl")
  include("test_state_estimation.jl")
  include("test_voltage_dependent_control.jl")
  include("test_transformer_phase_shift.jl")
  include("test_tap_controller.jl")
  include("test_series_reactance_control.jl")
  include("test_hvdc_pair_control.jl")
  include("test_terminal_status.jl")
  include("test_tap_changer_model.jl")
  include("test_phase_tap_changer_model.jl")
  include("test_phase_tap_table.jl")  
  include("test_configuration_coverage.jl")
  include("test_matpower_metadata.jl")
  include("test_api.jl")
  include("test_webui.jl")
  include("test_dc_powerflow.jl")
  include("test_distributed_slack.jl")
  include("test_island_diagnostics.jl")
  include("test_short_circuit.jl")
  include("test_external_grid.jl")
end

function include_extended_tests()
  # testgrid.jl defines run_grid_extended_tests (and shared helpers); in the
  # :all profile it is already loaded by include_fast_tests, so only include
  # it when the :extended profile runs standalone.
  isdefined(@__MODULE__, :run_grid_extended_tests) || include("testgrid.jl")
  include("test_3wt_phase_taps.jl")
  include("test_solver_interface.jl")
  include("test_klu_linear_solver.jl")
  include("test_api_extended.jl")
  include("test_webui_extended.jl")
  include("testremove.jl")
  include("test_pv_voltage_residuals.jl")
  include("test_matpower_example.jl")
  include("test_net_cache.jl")
  include("test_synthetic_grids.jl")
  include("test_configuration_docs.jl")
  include("test_repository_hygiene.jl")
  include("test_apslf.jl")
  include("test_cgmes_importer.jl")
  include("test_cgmes_export.jl")
  include("extended/test_dtf_importer.jl")
  include("extended/test_dtf_for002_validation_example.jl")
  include("extended/test_dtf_for002_outage_validation_example.jl")
  include("extended/test_dtf_matpower_export_validation_example.jl")
  include("extended/test_dtf_api_webui_integration.jl")
  # Experimental large-case comparison tooling is excluded from normal profiles.
end

function run_fast_profile_tests()
  function run_entry(name::Symbol)
    runner = Base.invokelatest(getfield, @__MODULE__, name)
    return Base.invokelatest(runner)
  end
  groups = [
    ("core_model", () -> begin
      run_entry(:run_grid_fast_tests)
      run_entry(:run_piline_g_tests)
    end),
    ("terminal_status", () -> run_entry(:run_terminal_status_tests)),
    ("configuration", () -> run_entry(:run_configuration_coverage_tests)),
    ("matpower_metadata", () -> run_entry(:run_matpower_metadata_tests)),
    ("programmatic_api", () -> run_entry(:run_api_fast_tests)),
    ("webui", () -> run_entry(:run_webui_fast_tests)),
    ("state_estimation", () -> run_entry(:run_state_estimation_tests)),
    ("dc_powerflow", () -> run_entry(:run_dc_powerflow_tests)),
    ("distributed_slack", () -> run_entry(:run_distributed_slack_tests)),
    ("island_diagnostics", () -> run_entry(:run_island_diagnostics_tests)),
    ("short_circuit", () -> run_entry(:run_short_circuit_tests)),
    ("external_grid", () -> run_entry(:run_external_grid_tests)),
    ("controls", () -> begin
      run_entry(:run_voltage_dependent_control_tests)
      run_entry(:run_transformer_phase_shift_tests)
      run_entry(:run_tap_controller_tests)
      run_entry(:run_series_reactance_control_tests)
      run_entry(:run_hvdc_pair_control_tests)
      run_entry(:run_tap_changer_model_tests)
      run_entry(:run_phase_tap_changer_model_tests)
      run_entry(:run_phase_tap_table_tests)
    end),
  ]
  @testset "Sparlectra.jl fast profile" begin
    total = length(groups)
    for (i, (name, runner)) in enumerate(groups)
      run_profile_group(i, total, name, runner)
    end
  end
end

function run_extended_profile_tests()
  function run_entry(name::Symbol)
    runner = Base.invokelatest(getfield, @__MODULE__, name)
    return Base.invokelatest(runner)
  end
  groups = [
    ("legacy/remove", () -> run_entry(:run_remove_tests)),
    ("core_model_extended", () -> run_entry(:run_grid_extended_tests)),
    ("powerflow_rectangular", () -> run_entry(:run_solver_interface_tests)),
    ("klu_linear_solver", () -> run_entry(:run_klu_linear_solver_tests)),
    ("3wt_phase_taps", () -> run_entry(:run_3wt_phase_taps_tests)),
    ("programmatic_api_extended", () -> run_entry(:run_api_extended_tests)),
    ("webui_extended", () -> run_entry(:run_webui_extended_tests)),
    ("pv_voltage_residuals", () -> run_entry(:run_pv_voltage_residual_tests)),
    ("matpower_examples", () -> run_entry(:run_matpower_example_tests)),
    ("net_cache", () -> run_entry(:run_net_cache_tests)),
    ("synthetic_grids", () -> run_entry(:run_synthetic_grid_tests)),
    ("configuration_docs", () -> run_entry(:run_configuration_docs_tests)),
    ("repository_hygiene", () -> run_entry(:run_repository_hygiene_tests)),
    ("apslf", () -> run_entry(:run_apslf_tests)),
    ("cgmes_importer", () -> run_entry(:run_cgmes_importer_tests)),
    ("cgmes_export", () -> run_entry(:run_cgmes_export_tests)),
    ("dtf_extended", () -> begin
      run_entry(:run_dtf_importer_tests)
      run_entry(:run_dtf_for002_validation_example_tests)
      run_entry(:run_dtf_for002_outage_validation_example_tests)
      run_entry(:run_dtf_matpower_export_validation_example_tests)
      run_entry(:run_dtf_api_webui_integration_tests)
    end),
  ]
  @testset "Sparlectra.jl extended profile" begin
    total = length(groups)
    for (i, (name, runner)) in enumerate(groups)
      run_profile_group(i, total, name, runner)
    end
  end
end

function main()
  if TEST_PROFILE === :fast
    print_test_progress_header(:fast)
    include_fast_tests()
    run_fast_profile_tests()
  elseif TEST_PROFILE === :extended
    print_test_progress_header(:extended)
    include_extended_tests()
    run_extended_profile_tests()
  elseif TEST_PROFILE === :all
    print_test_progress_header(:all)
    include_fast_tests()
    include_extended_tests()
    run_fast_profile_tests()
    run_extended_profile_tests()
  else
    error("Unknown test profile=$(TEST_PROFILE). Allowed: fast, extended, all. Selection precedence: CLI arg, SPARLECTRA_TEST_PROFILE, default fast.")
  end
end

Base.invokelatest(main)
return nothing

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
  include("test_solver_interface.jl")
  include("test_state_estimation.jl")
  include("test_voltage_dependent_control.jl")
  include("test_transformer_phase_shift.jl")
  include("test_tap_controller.jl")
  include("test_configuration_coverage.jl")
  include("test_matpower_metadata.jl")
  include("test_api.jl")
  include("test_webui.jl")
end

function include_extended_tests()
  include("testgrid.jl")
  include("test_api_extended.jl")
  include("test_webui_extended.jl")
  include("testremove.jl")
  include("test_pv_voltage_residuals.jl")
  include("test_matpower_example.jl")
  include("test_synthetic_grids.jl")
  include("test_configuration_docs.jl")
  include("test_repository_hygiene.jl")
  # Experimental large-case comparison tooling is excluded from normal profiles.
end

function run_fast_profile_tests()
  function run_entry(name::Symbol)
    runner = Base.invokelatest(getfield, @__MODULE__, name)
    return Base.invokelatest(runner)
  end
  groups = [
    ("core_model", () -> run_entry(:run_grid_fast_tests)),
    ("powerflow_rectangular", () -> run_entry(:run_solver_interface_tests)),
    ("configuration", () -> run_entry(:run_configuration_coverage_tests)),
    ("matpower_metadata", () -> run_entry(:run_matpower_metadata_tests)),
    ("programmatic_api", () -> run_entry(:run_api_fast_tests)),
    ("webui", () -> run_entry(:run_webui_fast_tests)),
    ("state_estimation", () -> run_entry(:run_state_estimation_tests)),
    ("controls", () -> begin
      run_entry(:run_voltage_dependent_control_tests)
      run_entry(:run_transformer_phase_shift_tests)
      run_entry(:run_tap_controller_tests)
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
    ("programmatic_api_extended", () -> run_entry(:run_api_extended_tests)),
    ("webui_extended", () -> run_entry(:run_webui_extended_tests)),
    ("pv_voltage_residuals", () -> run_entry(:run_pv_voltage_residual_tests)),
    ("matpower_examples", () -> run_entry(:run_matpower_example_tests)),
    ("synthetic_grids", () -> run_entry(:run_synthetic_grid_tests)),
    ("configuration_docs", () -> run_entry(:run_configuration_docs_tests)),
    ("repository_hygiene", () -> run_entry(:run_repository_hygiene_tests)),
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

using Sparlectra
using Test
using Logging
using Printf
using LinearAlgebra
using SparseArrays

global_logger(ConsoleLogger(stderr, Logging.Warn))
const TEST_PROFILE = Symbol(get(ENV, "SPARLECTRA_TEST_PROFILE", "fast"))

function include_fast_tests()
  include("testgrid.jl")
  include("test_solver_interface.jl")
  include("test_state_estimation.jl")
  include("test_voltage_dependent_control.jl")
  include("test_transformer_phase_shift.jl")
  include("test_tap_controller.jl")
  include("test_configuration_coverage.jl")
end

function include_extended_tests()
  include("testremove.jl")
  include("test_pv_voltage_residuals.jl")
  include("test_matpower_example.jl")
  include("test_synthetic_grids.jl")
  include("test_configuration_docs.jl")
end

function run_fast_profile_tests()
  function run_entry(name::Symbol)
    runner = Base.invokelatest(getfield, @__MODULE__, name)
    return Base.invokelatest(runner)
  end
  groups = [
    ("core_model", () -> run_entry(:run_grid_tests)),
    ("powerflow_rectangular", () -> run_entry(:run_solver_interface_tests)),
    ("configuration", () -> run_entry(:run_configuration_coverage_tests)),
    ("state_estimation", () -> run_entry(:run_state_estimation_tests)),
    ("controls", () -> begin
      run_entry(:run_voltage_dependent_control_tests)
      run_entry(:run_transformer_phase_shift_tests)
      run_entry(:run_tap_controller_tests)
    end),
  ]
  @testset "Sparlectra.jl fast profile" begin
    for (_, runner) in groups
      Base.invokelatest(runner)
    end
  end
end

function run_extended_profile_tests()
  function run_entry(name::Symbol)
    runner = Base.invokelatest(getfield, @__MODULE__, name)
    return Base.invokelatest(runner)
  end
  @testset "Sparlectra.jl extended profile" begin
    run_entry(:run_remove_tests)
    run_entry(:run_pv_voltage_residual_tests)
    run_entry(:run_matpower_example_tests)
    run_entry(:run_synthetic_grid_tests)
    run_entry(:run_configuration_docs_tests)
  end
end

if TEST_PROFILE === :fast
  include_fast_tests()
  run_fast_profile_tests()
elseif TEST_PROFILE === :extended
  include_fast_tests()
  include_extended_tests()
  run_fast_profile_tests()
  run_extended_profile_tests()
elseif TEST_PROFILE === :all
  include_fast_tests()
  include_extended_tests()
  # At the moment, `all` is an alias for `extended`.
  # Keep the profile for future all-only suites and CI matrix clarity.
  run_fast_profile_tests()
  run_extended_profile_tests()
else
  error("Unknown SPARLECTRA_TEST_PROFILE=$(TEST_PROFILE). Allowed: fast, extended, all")
end
return nothing

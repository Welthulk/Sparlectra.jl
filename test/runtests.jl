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
# purpose: test suite entry point: selects a profile (fast, pf, se, config,
#          webui, extd, extended, all) from ARGS/SPARLECTRA_TEST_PROFILE,
#          includes the test files of its groups, and runs the grouped
#          testsets with quiet output capture
using Sparlectra
# the application package (service layer, Web UI) lives in app/ with its own
# environment; it is loaded from there for the groups that test it and for
# the smoke test of the fast profile
pushfirst!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "app")))
using SparlectraApp
using Test
using Logging
using Printf
using LinearAlgebra
using SparseArrays

# Warnings the package raises ON PURPOSE while a test exercises the behaviour
# behind them (a deprecated configuration value, a version-less YAML file,
# non-case keys handed to write_case_config, an advisory topology precheck)
# are part of what is tested, not of the report. They are dropped here unless
# SPARLECTRA_TEST_SHOW_WARNINGS=1. Errors, and anything logged from outside a
# Sparlectra module, still reach stderr. A testset that asserts a warning
# installs its own logger through @test_logs and is not affected.
struct QuietPackageWarnings <: Logging.AbstractLogger
    inner::Logging.AbstractLogger
end
_logged_from_sparlectra(m) = m isa Module && string(Base.moduleroot(m)) == "Sparlectra"
Logging.min_enabled_level(l::QuietPackageWarnings) = Logging.min_enabled_level(l.inner)
Logging.catch_exceptions(l::QuietPackageWarnings) = Logging.catch_exceptions(l.inner)
function Logging.shouldlog(l::QuietPackageWarnings, level, _module, group, id)
    level == Logging.Warn && _logged_from_sparlectra(_module) && return false
    return Logging.shouldlog(l.inner, level, _module, group, id)
end
Logging.handle_message(l::QuietPackageWarnings, args...; kwargs...) = Logging.handle_message(l.inner, args...; kwargs...)

const SHOW_PACKAGE_WARNINGS = get(ENV, "SPARLECTRA_TEST_SHOW_WARNINGS", "0") == "1"
global_logger(SHOW_PACKAGE_WARNINGS ? ConsoleLogger(stderr, Logging.Warn) : QuietPackageWarnings(ConsoleLogger(stderr, Logging.Warn)))

include("test_runner_helpers.jl")

const TEST_PROFILE = selected_test_profile()

# Build workloads run this suite as a PackageCompiler precompile trace and
# can leave out groups that do not belong in the target image via
# SPARLECTRA_TEST_SKIP_GROUPS. The skip is printed explicitly so it can
# never be mistaken for a pass. Default: empty, every group runs.
const SKIP_GROUPS = Set{String}(String(strip(g)) for g in split(get(ENV, "SPARLECTRA_TEST_SKIP_GROUPS", ""), ','; keepempty=false))

function print_test_progress_header(profile::Symbol)
    println("Test framework: ", profile)
end

function print_group_progress(i::Int, total::Int, name::AbstractString) end

# One row per test group: the files that define its runners (included once,
# in this order, shared files only once across groups) and the runner
# functions called in order. The profiles below are lists of group names, so
# a group belongs to exactly one profile and `all` is the union.
struct TestGroup
    name::String
    files::Vector{String}
    runners::Vector{Symbol}
end

const TEST_GROUPS = TestGroup[
    # --- fast: the pull-request gate, model and Newton core, seconds not minutes
    TestGroup("core_model", ["testgrid.jl", "test_piline_g.jl"], [:run_grid_fast_tests, :run_piline_g_tests]),
    TestGroup("terminal_status", ["test_terminal_status.jl"], [:run_terminal_status_tests]),
    TestGroup("powerflow_rectangular", ["test_solver_interface.jl"], [:run_solver_interface_tests]),
    TestGroup("factorized_linear_solver", ["test_factorized_linear_solver.jl"], [:run_factorized_linear_solver_tests]),
    TestGroup("pv_voltage_residuals", ["test_pv_voltage_residuals.jl"], [:run_pv_voltage_residual_tests]),
    TestGroup("3wt_phase_taps", ["test_3wt_phase_taps.jl"], [:run_3wt_phase_taps_tests]),
    TestGroup("dc_powerflow", ["test_dc_powerflow.jl"], [:run_dc_powerflow_tests]),
    TestGroup("distributed_slack", ["test_distributed_slack.jl"], [:run_distributed_slack_tests]),
    TestGroup("island_diagnostics", ["test_island_diagnostics.jl"], [:run_island_diagnostics_tests]),
    TestGroup("external_grid", ["test_external_grid.jl"], [:run_external_grid_tests]),
    TestGroup("matpower_metadata", ["test_matpower_metadata.jl"], [:run_matpower_metadata_tests]),
    TestGroup("programmatic_api", ["test_api.jl"], [:run_api_fast_tests]),
    # --- pf: the rest of the power-flow surface, controllers, contingencies, scenarios
    TestGroup("auto_powerflow", ["test_auto_powerflow.jl"], [:run_auto_powerflow_tests]),
    TestGroup("short_circuit", ["test_short_circuit.jl"], [:run_short_circuit_tests]),
    TestGroup("parallel_foundation", ["test_parallel_foundation.jl"], [:run_parallel_foundation_tests]),
    TestGroup("contingency", ["test_contingency.jl"], [:run_contingency_tests]),
    TestGroup("scenarios", ["test_scenarios.jl"], [:run_scenario_patch_tests]),
    TestGroup("controls", ["test_voltage_dependent_control.jl", "test_transformer_phase_shift.jl", "test_tap_controller.jl", "test_series_reactance_control.jl", "test_upfc_control.jl", "test_hvdc_pair_control.jl", "test_tap_changer_model.jl", "test_phase_tap_changer_model.jl", "test_phase_tap_table.jl"],
        [:run_voltage_dependent_control_tests, :run_transformer_phase_shift_tests, :run_tap_controller_tests, :run_series_reactance_control_tests, :run_upfc_control_tests, :run_hvdc_pair_control_tests, :run_tap_changer_model_tests, :run_phase_tap_changer_model_tests, :run_phase_tap_table_tests]),
    TestGroup("core_model_extended", ["testgrid.jl"], [:run_grid_extended_tests]),
    TestGroup("contingency_extended", ["test_contingency.jl"], [:run_contingency_extended_tests]),
    TestGroup("scenario_engine", ["test_scenarios.jl"], [:run_scenario_engine_extended_tests]),
    TestGroup("apslf", ["test_apslf.jl"], [:run_apslf_tests]),
    # --- se: state estimation, observability, topology
    TestGroup("state_estimation", ["test_state_estimation.jl"], [:run_state_estimation_tests]),
    TestGroup("observability", ["test_observability.jl"], [:run_observability_tests]),
    TestGroup("topology_validation", ["test_topology_validation.jl"], [:run_topology_validation_tests]),
    # --- config: configuration surface, its documentation, repository hygiene
    TestGroup("configuration", ["test_configuration_coverage.jl"], [:run_configuration_coverage_tests]),
    TestGroup("configuration_docs", ["test_configuration_docs.jl"], [:run_configuration_docs_tests]),
    TestGroup("repository_hygiene", ["test_repository_hygiene.jl"], [:run_repository_hygiene_tests]),
    # --- webui: the local browser UI
    TestGroup("webui", ["test_webui.jl"], [:run_webui_fast_tests]),
    TestGroup("webui_extended", ["test_webui_extended.jl"], [:run_webui_extended_tests]),
    # --- extd: shipped cases, formats, service layer, examples, fixtures
    TestGroup("demo_cases", ["test_demo_cases.jl"], [:run_demo_case_tests]),
    TestGroup("scf", ["test_scf.jl"], [:run_scf_tests]),
    TestGroup("programmatic_api_extended", ["test_api_extended.jl"], [:run_api_extended_tests]),
    TestGroup("matpower_examples", ["test_matpower_example.jl"], [:run_matpower_example_tests]),
    TestGroup("example_infra", ["test_example_suites.jl"], [:run_example_suite_infra_tests]),
    TestGroup("net_cache", ["test_net_cache.jl"], [:run_net_cache_tests]),
    TestGroup("synthetic_grids", ["test_synthetic_grids.jl"], [:run_synthetic_grid_tests]),
    TestGroup("cgmes_importer", ["test_cgmes_importer.jl"], [:run_cgmes_importer_tests]),
    TestGroup("cgmes_export", ["test_cgmes_export.jl"], [:run_cgmes_export_tests]),
    TestGroup("dtf_extended", ["extended/test_dtf_importer.jl", "extended/test_dtf_for002_validation_example.jl", "extended/test_dtf_for002_outage_validation_example.jl", "extended/test_dtf_matpower_export_validation_example.jl", "extended/test_dtf_api_webui_integration.jl"],
        [:run_dtf_importer_tests, :run_dtf_for002_validation_example_tests, :run_dtf_for002_outage_validation_example_tests, :run_dtf_matpower_export_validation_example_tests, :run_dtf_api_webui_integration_tests]),
]
# Experimental large-case comparison tooling is excluded from every profile.

# The profiles. `fast` is the pull-request gate and stays short on purpose;
# a change to the solver or the estimator runs its own profile on top;
# `extended` is everything that is not fast (the former second profile);
# `all` is both. The documentation build is a gate of its own
# (tools/run_gates.sh docs), not a test profile.
const TEST_PROFILES = Dict{Symbol,Vector{String}}(
    :fast => ["core_model", "terminal_status", "powerflow_rectangular", "factorized_linear_solver", "pv_voltage_residuals", "3wt_phase_taps", "dc_powerflow", "distributed_slack", "island_diagnostics", "external_grid", "matpower_metadata", "programmatic_api"],
    :pf => ["auto_powerflow", "short_circuit", "parallel_foundation", "contingency", "scenarios", "controls", "core_model_extended", "contingency_extended", "scenario_engine", "apslf"],
    :se => ["state_estimation", "observability", "topology_validation"],
    :config => ["configuration", "configuration_docs", "repository_hygiene"],
    :webui => ["webui", "webui_extended"],
    :extd => ["demo_cases", "scf", "programmatic_api_extended", "matpower_examples", "example_infra", "net_cache", "synthetic_grids", "cgmes_importer", "cgmes_export", "dtf_extended"],
)
TEST_PROFILES[:extended] = vcat(TEST_PROFILES[:pf], TEST_PROFILES[:se], TEST_PROFILES[:config], TEST_PROFILES[:webui], TEST_PROFILES[:extd])
TEST_PROFILES[:all] = vcat(TEST_PROFILES[:fast], TEST_PROFILES[:extended])

# every group sits in exactly one of the six base profiles
let seen = String[]
    for key in (:fast, :pf, :se, :config, :webui, :extd), name in TEST_PROFILES[key]
        name in seen && error("test group $(name) is listed in two profiles")
        any(g -> g.name == name, TEST_GROUPS) || error("profile $(key) names an unknown test group $(name)")
        push!(seen, name)
    end
    missing_groups = [g.name for g in TEST_GROUPS if !(g.name in seen)]
    isempty(missing_groups) || error("test group(s) in no profile: " * join(missing_groups, ", "))
end

profile_groups(profile::Symbol) = [g for name in TEST_PROFILES[profile] for g in TEST_GROUPS if g.name == name]

"""
    include_group_files(groups)

Include the test files of `groups` once each, in group order. A file shared
by two groups (testgrid.jl, test_contingency.jl, test_scenarios.jl) is
included on its first use only.
"""
function include_group_files(groups::Vector{TestGroup})
    done = Set{String}()
    for g in groups, f in g.files
        f in done && continue
        include(f)
        push!(done, f)
    end
    return nothing
end

function run_profile_groups(profile::Symbol, groups::Vector{TestGroup})
    function run_entry(name::Symbol)
        runner = Base.invokelatest(getfield, @__MODULE__, name)
        return Base.invokelatest(runner)
    end
    skipped = [g.name for g in groups if g.name in SKIP_GROUPS]
    isempty(skipped) || println("Skipped group(s) via SPARLECTRA_TEST_SKIP_GROUPS: ", join(skipped, ", "))
    selected = [g for g in groups if !(g.name in SKIP_GROUPS)]
    @testset "Sparlectra.jl $(profile) profile" begin
        total = length(selected)
        for (i, g) in enumerate(selected)
            run_profile_group(i, total, g.name, () -> foreach(run_entry, g.runners))
        end
    end
end

"""
    timed_include(label, f)

Run the include phase of a profile and print what it cost.

The per-group `PASS` lines account for everything the GROUPS do, and for
nothing that happens before them: including the test files parses and
compiles every testset body in the profile, and that time is attributed
nowhere. Measuring it is the difference between "the tests are slow" and
"loading the tests is slow", which is not the same lever.
"""
function timed_include(label::AbstractString, f::Function)
    timed = @timed f()
    @printf("include %s: %.3f s (%.3f s compile, %.3f s recompile), %.1f MiB allocated\n",
        label, timed.time, timed.compile_time, timed.recompile_time, timed.bytes / 1024.0^2)
    return timed.value
end

function main()
    haskey(TEST_PROFILES, TEST_PROFILE) || error("Unknown test profile=$(TEST_PROFILE). Allowed: " * join(string.(sort(collect(keys(TEST_PROFILES)))), ", ") * ". Selection precedence: CLI arg, SPARLECTRA_TEST_PROFILE, default fast.")
    print_test_progress_header(TEST_PROFILE)
    groups = profile_groups(TEST_PROFILE)
    timed_include(string(TEST_PROFILE), () -> include_group_files(groups))
    run_profile_groups(TEST_PROFILE, groups)
end

Base.invokelatest(main)
return nothing
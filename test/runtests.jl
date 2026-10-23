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
#          webui, extd, cgmes, formats, legacy, install, workshops, extended,
#          all, threads)
#          from ARGS/SPARLECTRA_TEST_PROFILE,
#          includes the test files of its groups, and runs the grouped
#          testsets with quiet output capture
using Sparlectra
# output.csv_format defaults to `auto` (0.32.0), which follows the regional
# settings of the machine: on a German desktop every CSV of the suite would
# carry semicolons, on CI commas. The suite pins the portable locale, which
# `auto` maps to technical on every system (a POSIX locale variable wins on
# Windows too); child processes the tests start inherit it. The test of
# `auto` itself sets its locales with withenv.
ENV["LC_ALL"] = "C.UTF-8"
# the application package (service layer, Web UI) lives in app/ with its own
# environment; it is loaded from there for the groups that test it and for
# the smoke test of the fast profile. A fresh checkout (Pkg.test on CI, a
# registered install) has no app/Manifest.toml yet: the environment is set
# up here once, otherwise SparlectraApp cannot even find the library.
const _SPARLECTRA_APP_DIR = normpath(joinpath(@__DIR__, "..", "app"))
if !isfile(joinpath(_SPARLECTRA_APP_DIR, "Manifest.toml"))
  import Pkg
  Pkg.activate(_SPARLECTRA_APP_DIR) do
    Pkg.instantiate()
  end
end
pushfirst!(LOAD_PATH, _SPARLECTRA_APP_DIR)
# BenchmarkTools is a weak dependency of the library (the extension behind
# the benchmark mode of run_matpower_case). Pkg.test brings it through the
# test target; a run as `julia --project=. test/runtests.jl`, the way the
# gates run, does not, since an [extras] entry is not loadable there. The
# runner then keeps a small environment with BenchmarkTools alone in the
# depot, set up once, and stacks it onto the load path.
if Base.find_package("BenchmarkTools") === nothing
  _benchmark_env = joinpath(first(DEPOT_PATH), "environments", "sparlectra-test-benchmarktools")
  if !isfile(joinpath(_benchmark_env, "Manifest.toml"))
    import Pkg
    Pkg.activate(_benchmark_env) do
      Pkg.add("BenchmarkTools"; io = devnull)
    end
  end
  push!(LOAD_PATH, _benchmark_env)
end
# KLU the same way (0.30.1): the extension behind power mode's sparse LU;
# the factorized-solver test loads it and reports whether it ran
if Base.find_package("KLU") === nothing
  _klu_env = joinpath(first(DEPOT_PATH), "environments", "sparlectra-test-klu")
  if !isfile(joinpath(_klu_env, "Manifest.toml"))
    import Pkg
    Pkg.activate(_klu_env) do
      Pkg.add("KLU"; io = devnull)
    end
  end
  push!(LOAD_PATH, _klu_env)
end
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
    # --- fast: the pull-request gate, results against an independent truth
    TestGroup("docstrings", ["test_docstrings.jl"], [:run_docstring_tests]),
    TestGroup("api_smoke", ["fast/ref_api_smoke.jl"], [:run_ref_api_smoke_tests]),
    TestGroup("hand_calculations", ["fast/ref_hand_calculations.jl"], [:run_ref_hand_calculation_tests]),
    TestGroup("matpower_header", ["fast/ref_matpower_header.jl"], [:run_ref_matpower_header_tests]),
    TestGroup("entsoe_sv", ["fast/ref_entsoe_sv.jl"], [:run_ref_entsoe_sv_tests]),
    TestGroup("openloadflow", ["fast/ref_openloadflow.jl"], [:run_ref_openloadflow_tests]),
    # --- pf: the model and Newton core (formerly fast), solver options, island handling
    TestGroup("core_model", ["testgrid.jl", "test_piline_g.jl"], [:run_grid_fast_tests, :run_piline_g_tests]),
    TestGroup("terminal_status", ["test_terminal_status.jl"], [:run_terminal_status_tests]),
    TestGroup("powerflow_rectangular", ["test_solver_interface.jl"], [:run_solver_interface_tests]),
    TestGroup("factorized_linear_solver", ["test_factorized_linear_solver.jl"], [:run_factorized_linear_solver_tests]),
    TestGroup("3wt_phase_taps", ["test_3wt_phase_taps.jl"], [:run_3wt_phase_taps_tests]),
    TestGroup("dc_powerflow", ["test_dc_powerflow.jl"], [:run_dc_powerflow_tests]),
    TestGroup("distributed_slack", ["test_distributed_slack.jl"], [:run_distributed_slack_tests]),
    TestGroup("island_diagnostics", ["test_island_diagnostics.jl"], [:run_island_diagnostics_tests]),
    TestGroup("external_grid", ["test_external_grid.jl"], [:run_external_grid_tests]),
    TestGroup("matpower_metadata", ["testgrid.jl", "test_matpower_metadata.jl"], [:run_matpower_metadata_tests]),
    TestGroup("programmatic_api", ["test_api.jl"], [:run_api_fast_tests]),
    # --- pf: the rest of the power-flow surface, controllers, contingencies, scenarios
    TestGroup("auto_powerflow", ["test_auto_powerflow.jl"], [:run_auto_powerflow_tests, :run_auto_powerflow_ladder_tests]),
    TestGroup("short_circuit", ["test_short_circuit.jl"], [:run_short_circuit_tests]),
    TestGroup("parallel_foundation", ["test_parallel_foundation.jl"], [:run_parallel_foundation_tests]),
    TestGroup("contingency", ["test_contingency.jl"], [:run_contingency_tests, :run_contingency_case_tests]),
    TestGroup("scenarios", ["test_scenarios.jl"], [:run_scenario_patch_tests]),
    TestGroup("controls", ["test_voltage_dependent_control.jl", "test_transformer_phase_shift.jl", "test_tap_controller.jl", "test_series_reactance_control.jl", "test_upfc_control.jl", "test_hvdc_pair_control.jl", "test_tap_changer_model.jl", "test_phase_tap_changer_model.jl", "test_phase_tap_table.jl"],
        [:run_voltage_dependent_control_tests, :run_transformer_phase_shift_tests, :run_tap_controller_tests, :run_series_reactance_control_tests, :run_upfc_control_tests, :run_hvdc_pair_control_tests, :run_tap_changer_model_tests, :run_phase_tap_changer_model_tests, :run_phase_tap_table_tests]),
    TestGroup("core_model_extended", ["testgrid.jl"], [:run_grid_extended_tests]),
    TestGroup("contingency_extended", ["test_contingency.jl"], [:run_contingency_extended_tests]),
    TestGroup("scenario_engine", ["test_scenarios.jl"], [:run_scenario_engine_extended_tests]),
    TestGroup("apslf", ["test_apslf.jl"], [:run_apslf_tests]),
    # --- se: state estimation, observability, topology
    TestGroup("state_estimation", ["test_state_estimation.jl"], [:run_state_estimation_tests]),
    TestGroup("observability", ["test_observability.jl"], [:run_observability_tests, :run_observability_rank_tests]),
    TestGroup("topology_validation", ["test_legacy_topology_validation.jl"], [:run_legacy_topology_validation_tests]),
    TestGroup("measurement_table", ["test_measurement_table.jl"], [:run_measurement_table_tests]),
    # --- config: configuration surface, its documentation, repository hygiene
    TestGroup("configuration", ["test_configuration_coverage.jl"], [:run_configuration_coverage_tests]),
    TestGroup("configuration_docs", ["test_configuration_docs.jl"], [:run_configuration_docs_tests]),
    TestGroup("repository_hygiene", ["test_repository_hygiene.jl"], [:run_repository_hygiene_tests]),
    # --- webui: the local browser UI
    TestGroup("webui", ["test_webui.jl"], [:run_webui_fast_tests]),
    TestGroup("webui_extended", ["test_webui_extended.jl"], [:run_webui_extended_tests]),
    TestGroup("webui_state_estimation", ["test_webui.jl", "test_webui_state_estimation.jl"], [:run_webui_state_estimation_tests]),
    # --- extd: shipped cases, formats, service layer, examples, fixtures
    TestGroup("demo_cases", ["test_demo_cases.jl"], [:run_demo_case_tests]),
    TestGroup("scf", ["test_scf.jl"], [:run_scf_tests]),
    TestGroup("scf_series", ["test_scf_series.jl"], [:run_scf_series_tests]),
    TestGroup("programmatic_api_extended", ["test_api_extended.jl"], [:run_api_extended_tests]),
    TestGroup("matpower_examples", ["test_matpower_example.jl"], [:run_matpower_example_tests]),
    TestGroup("example_infra", ["test_example_suites.jl"], [:run_example_suite_infra_tests]),
    TestGroup("synthetic_grids", ["test_synthetic_grids.jl"], [:run_synthetic_grid_tests]),
    TestGroup("cgmes_importer", ["test_cgmes_importer.jl"], [:run_cgmes_importer_tests]),
    TestGroup("cgmes_export", ["test_cgmes_export.jl"], [:run_cgmes_export_tests]),
    TestGroup("cgmes_series", ["test_cgmes_series.jl"], [:run_cgmes_series_tests, :run_cgmes_union_tests]),
    TestGroup("cgmes_ssh_series", ["test_cgmes_ssh_series.jl"], [:run_cgmes_ssh_series_tests]),
    TestGroup("powsybl_importer", ["test_powsybl_importer.jl"], [:run_powsybl_importer_tests]),
    # --- install: the installation path of a fresh checkout, before a release
    TestGroup("install", ["test_install.jl"], [:run_install_tests]),
    TestGroup("workshops", ["test_workshops.jl"], [:run_workshop_tests]),
    TestGroup("dtf_extended", ["extended/test_dtf_importer.jl", "extended/test_dtf_for002_validation_example.jl", "extended/test_dtf_for002_outage_validation_example.jl", "extended/test_dtf_matpower_export_validation_example.jl", "extended/test_dtf_api_webui_integration.jl"],
        [:run_dtf_importer_tests, :run_dtf_for002_validation_example_tests, :run_dtf_for002_outage_validation_example_tests, :run_dtf_matpower_export_validation_example_tests, :run_dtf_api_webui_integration_tests]),
    # --- legacy: one group per test file that lost testsets to the legacy profile;
    # it includes that file for its helpers, then the test_legacy_ file next to it
    TestGroup("legacy_api", ["test_api.jl", "test_legacy_api.jl"], [:run_legacy_api_tests]),
    TestGroup("legacy_apslf", ["test_apslf.jl", "test_legacy_apslf.jl"], [:run_legacy_apslf_tests]),
    TestGroup("legacy_cgmes_export", ["test_cgmes_export.jl", "test_legacy_cgmes_export.jl"], [:run_legacy_cgmes_export_tests]),
    TestGroup("legacy_cgmes_importer", ["test_cgmes_importer.jl", "test_legacy_cgmes_importer.jl"], [:run_legacy_cgmes_importer_tests]),
    TestGroup("legacy_dc_powerflow", ["test_dc_powerflow.jl", "test_legacy_dc_powerflow.jl"], [:run_legacy_dc_powerflow_tests]),
    TestGroup("legacy_demo_cases", ["test_demo_cases.jl", "test_legacy_demo_cases.jl"], [:run_legacy_demo_cases_tests]),
    TestGroup("legacy_external_grid", ["test_external_grid.jl", "test_legacy_external_grid.jl"], [:run_legacy_external_grid_tests]),
    TestGroup("legacy_factorized_linear_solver", ["test_factorized_linear_solver.jl", "test_legacy_factorized_linear_solver.jl"], [:run_legacy_factorized_linear_solver_tests]),
    TestGroup("legacy_hvdc_pair_control", ["test_hvdc_pair_control.jl", "test_legacy_hvdc_pair_control.jl"], [:run_legacy_hvdc_pair_control_tests]),
    TestGroup("legacy_island_diagnostics", ["test_island_diagnostics.jl", "test_legacy_island_diagnostics.jl"], [:run_legacy_island_diagnostics_tests]),
    TestGroup("legacy_matpower_metadata", ["test_matpower_metadata.jl", "test_legacy_matpower_metadata.jl"], [:run_legacy_matpower_metadata_tests]),
    TestGroup("legacy_powsybl_importer", ["test_powsybl_importer.jl", "test_legacy_powsybl_importer.jl"], [:run_legacy_powsybl_importer_tests]),
    TestGroup("legacy_pv_voltage_residuals", ["testgrid.jl", "test_legacy_pv_voltage_residuals.jl"], [:run_legacy_pv_voltage_residuals_tests]),
    TestGroup("legacy_scenarios", ["test_scenarios.jl", "test_legacy_scenarios.jl"], [:run_legacy_scenarios_tests]),
    TestGroup("legacy_scf", ["test_scf.jl", "test_legacy_scf.jl"], [:run_legacy_scf_tests]),
    TestGroup("legacy_series_reactance_control", ["test_series_reactance_control.jl", "test_legacy_series_reactance_control.jl"], [:run_legacy_series_reactance_control_tests]),
    TestGroup("legacy_solver_interface", ["test_solver_interface.jl", "test_legacy_solver_interface.jl"], [:run_legacy_solver_interface_tests]),
    TestGroup("legacy_tap_controller", ["test_tap_controller.jl", "test_legacy_tap_controller.jl"], [:run_legacy_tap_controller_tests]),
    TestGroup("legacy_upfc_control", ["test_upfc_control.jl", "test_legacy_upfc_control.jl"], [:run_legacy_upfc_control_tests]),
    TestGroup("legacy_voltage_dependent_control", ["test_voltage_dependent_control.jl", "test_legacy_voltage_dependent_control.jl"], [:run_legacy_voltage_dependent_control_tests]),
    TestGroup("legacy_webui", ["test_webui.jl", "test_legacy_webui.jl"], [:run_legacy_webui_tests]),
    TestGroup("legacy_webui_extended", ["test_webui_extended.jl", "test_legacy_webui_extended.jl"], [:run_legacy_webui_extended_tests]),
    TestGroup("legacy_dtf_importer", ["extended/test_dtf_importer.jl", "extended/test_legacy_dtf_importer.jl"], [:run_legacy_dtf_importer_tests]),
    TestGroup("legacy_dtf_matpower_export_validation_example", ["extended/test_dtf_matpower_export_validation_example.jl", "extended/test_legacy_dtf_matpower_export_validation_example.jl"], [:run_legacy_dtf_matpower_export_validation_example_tests]),
    TestGroup("legacy_state_estimation", ["test_state_estimation.jl", "test_legacy_state_estimation.jl"], [:run_legacy_state_estimation_tests]),
    TestGroup("legacy_measurement_table", ["test_measurement_table.jl", "test_legacy_measurement_table.jl"], [:run_legacy_measurement_table_tests]),
    # --- legacy testsets (appended)
]

# The profiles. `fast` is the pull-request gate: results against an
# independent truth, nothing else, and it stays short on purpose;
# a change to the solver or the estimator runs its own profile on top;
# `extended` is everything that is not fast, legacy included;
# `install` is the installation path of a fresh checkout, due before a
# release; `all` is all of them. The documentation build is a gate of its own
# (tools/run_gates.sh docs), not a test profile.
const TEST_PROFILES = Dict{Symbol,Vector{String}}(
    :fast => ["docstrings", "api_smoke", "hand_calculations", "matpower_header", "entsoe_sv", "openloadflow"],
    :pf => ["core_model", "terminal_status", "powerflow_rectangular", "factorized_linear_solver", "3wt_phase_taps", "dc_powerflow", "distributed_slack", "island_diagnostics", "external_grid", "auto_powerflow", "short_circuit", "parallel_foundation", "contingency", "scenarios", "controls", "core_model_extended", "contingency_extended", "scenario_engine", "apslf"],
    :se => ["state_estimation", "observability", "topology_validation", "measurement_table"],
    :config => ["configuration", "configuration_docs", "repository_hygiene"],
    :webui => ["webui", "webui_extended", "webui_state_estimation"],
    # --- extd: the service layer and the shipped demo cases
    :extd => ["programmatic_api", "programmatic_api_extended", "demo_cases"],
    # --- cgmes and formats: the former adapters profile, by importer family
    :cgmes => ["cgmes_importer", "cgmes_export", "cgmes_series", "cgmes_ssh_series"],
    :formats => ["matpower_metadata", "scf", "scf_series", "matpower_examples", "synthetic_grids", "dtf_extended", "powsybl_importer"],
    # --- install: the installation path of a fresh checkout and the tests that start child processes
    :install => ["install", "example_infra"],
    :workshops => ["workshops"],
    # --- legacy: testsets that pin an earlier change; before a release and on a change to a protected file
    :legacy => String[g.name for g in TEST_GROUPS if startswith(g.name, "legacy_")],
)
TEST_PROFILES[:extended] = vcat(TEST_PROFILES[:pf], TEST_PROFILES[:se], TEST_PROFILES[:config], TEST_PROFILES[:webui], TEST_PROFILES[:extd], TEST_PROFILES[:cgmes], TEST_PROFILES[:formats], TEST_PROFILES[:legacy])
TEST_PROFILES[:all] = vcat(TEST_PROFILES[:fast], TEST_PROFILES[:extended], TEST_PROFILES[:install], TEST_PROFILES[:workshops])
# The two strands of extended (`tools/run_gates.sh extended --strands=2`):
# two processes at once, together exactly extended. Split by area, so each
# strand compiles the code of its areas once; the split is chosen so that
# both strands take about the same time.
TEST_PROFILES[:extended_a] = vcat(TEST_PROFILES[:cgmes], TEST_PROFILES[:se], TEST_PROFILES[:config], TEST_PROFILES[:legacy])
TEST_PROFILES[:extended_b] = vcat(TEST_PROFILES[:pf], TEST_PROFILES[:webui], TEST_PROFILES[:extd], TEST_PROFILES[:formats])
sort(vcat(TEST_PROFILES[:extended_a], TEST_PROFILES[:extended_b])) == sort(TEST_PROFILES[:extended]) || error("the strands extended_a and extended_b do not add up to extended")
# The three test strands of the full gate set (`tools/run_gates.sh set`):
# set_a and extended_b run at once with set_c, which takes fast and legacy
# (the docs build and the workshops follow it in the same strand); the
# threads profile runs alone afterwards. Together exactly fast + extended.
TEST_PROFILES[:set_a] = vcat(TEST_PROFILES[:cgmes], TEST_PROFILES[:se], TEST_PROFILES[:config])
TEST_PROFILES[:set_c] = vcat(TEST_PROFILES[:fast], TEST_PROFILES[:legacy])
sort(vcat(TEST_PROFILES[:set_a], TEST_PROFILES[:extended_b], TEST_PROFILES[:set_c])) == sort(vcat(TEST_PROFILES[:fast], TEST_PROFILES[:extended])) || error("the strands set_a, extended_b and set_c do not add up to fast and extended")

# every group sits in exactly one of the base profiles
let seen = String[]
    for key in (:fast, :pf, :se, :config, :webui, :extd, :cgmes, :formats, :legacy, :install, :workshops), name in TEST_PROFILES[key]
        name in seen && error("test group $(name) is listed in two profiles")
        any(g -> g.name == name, TEST_GROUPS) || error("profile $(key) names an unknown test group $(name)")
        push!(seen, name)
    end
    missing_groups = [g.name for g in TEST_GROUPS if !(g.name in seen)]
    isempty(missing_groups) || error("test group(s) in no profile: " * join(missing_groups, ", "))
end

# The threads profile: the parallel-vs-serial comparisons of every
# parallel site (short-circuit sweep, island power flow, N-1 batch,
# island-wise SE, QSTS), repeated THREAD_REPEATS times in ONE process that
# must run with several threads (tools/run_gates.sh threads starts it with
# --threads=4). In a single-threaded process every parallel site takes its
# serial loop, so these comparisons check the serial path against itself;
# that is how the base profiles run them. A single threaded run can pass by
# luck, hence the repeats. These runners are also called inside their area
# groups, so they are not base-profile groups and stay out of TEST_GROUPS.
const THREAD_GROUPS = TestGroup[
    TestGroup("short_circuit_parallel", ["test_short_circuit.jl"], [:run_short_circuit_parallel_tests]),
    TestGroup("island_powerflow_parallel", ["test_parallel_foundation.jl"], [:run_parallel_island_tests]),
    TestGroup("contingency_parallel", ["test_contingency.jl"], [:run_contingency_parallel_tests, :run_contingency_extended_tests]),
    TestGroup("state_estimation_parallel", ["test_state_estimation.jl"], [:run_state_estimation_parallel_tests]),
    TestGroup("qsts_parallel", ["test_scf_series.jl"], [:run_scf_series_parallel_tests]),
    TestGroup("cgmes_ssh_series_parallel", ["test_cgmes_ssh_series.jl"], [:run_cgmes_ssh_series_parallel_tests]),
]
const THREAD_REPEATS = 5

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

function run_entry(name::Symbol)
    runner = Base.invokelatest(getfield, @__MODULE__, name)
    return Base.invokelatest(runner)
end

function run_profile_groups(profile::Symbol, groups::Vector{TestGroup})
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

# the threads profile: every group of THREAD_GROUPS once per repeat, one
# result line each; the comparison report lines (RAN with n threads) show
# under every line
function run_thread_groups()
    Threads.nthreads() > 1 || error("the threads profile needs a multi-threaded process (julia --threads=4 ...; tools/run_gates.sh threads starts it so); with one thread every comparison checks the serial path against itself")
    println("threads: ", Threads.nthreads(), ", ", THREAD_REPEATS, " repeats")
    @testset "Sparlectra.jl threads profile" begin
        total = THREAD_REPEATS * length(THREAD_GROUPS)
        for rep in 1:THREAD_REPEATS, (k, g) in enumerate(THREAD_GROUPS)
            i = (rep - 1) * length(THREAD_GROUPS) + k
            run_profile_group(i, total, string(g.name, " (run ", rep, "/", THREAD_REPEATS, ")"), () -> foreach(run_entry, g.runners))
        end
    end
end

function main()
    if TEST_PROFILE === :threads
        print_test_progress_header(TEST_PROFILE)
        timed_include("threads", () -> include_group_files(THREAD_GROUPS))
        return run_thread_groups()
    end
    haskey(TEST_PROFILES, TEST_PROFILE) || error("Unknown test profile=$(TEST_PROFILE). Allowed: " * join(string.(sort(vcat(collect(keys(TEST_PROFILES)), :threads))), ", ") * ". Selection precedence: CLI arg, SPARLECTRA_TEST_PROFILE, default fast.")
    print_test_progress_header(TEST_PROFILE)
    groups = profile_groups(TEST_PROFILE)
    timed_include(string(TEST_PROFILE), () -> include_group_files(groups))
    run_profile_groups(TEST_PROFILE, groups)
end

Base.invokelatest(main)
return nothing
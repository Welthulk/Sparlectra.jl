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

# file: test/test_example_suites.jl
# purpose: tests the example suite runners and their shared backend: CLI
#          parsing, skip/filter logic, CSV and Markdown escaping, registry
#          integrity, and a --help subprocess smoke test
#
# Tests the example suite runners (examples/run_powerflow_suite.jl,
# run_state_estimation_suite.jl, run_others_suite.jl) and their shared
# backend examples/internal/example_suite_runner.jl: CLI parsing, skip/filter
# logic, CSV escaping, registry integrity, and a --help subprocess smoke
# test. The registered examples themselves are not executed here; their
# coverage lives in the dedicated example tests.
#
# Deliberately excluded from the fast/extended test profiles — the example
# suites are user-facing demos, not part of the library test surface. Run
# manually with:
#   julia --project=. -e 'using Test; include("test/test_example_suites.jl")'

const _EXAMPLES_DIR = normpath(joinpath(@__DIR__, "..", "examples"))
const _SUITE_SCRIPTS = ["run_powerflow_suite.jl", "run_state_estimation_suite.jl", "run_others_suite.jl", "run_short_circuit_suite.jl"]

function _sandbox_module()
  m = Module()
  Core.eval(m, :(include(path) = Base.include($m, path)))
  return m
end

function _include_with_no_main(path::AbstractString)
  previous = get(ENV, "SPARLECTRA_EXAMPLE_SUITE_NO_MAIN", nothing)
  ENV["SPARLECTRA_EXAMPLE_SUITE_NO_MAIN"] = "1"
  try
    m = _sandbox_module()
    Base.include(m, path)
    return m
  finally
    previous === nothing ? delete!(ENV, "SPARLECTRA_EXAMPLE_SUITE_NO_MAIN") : (ENV["SPARLECTRA_EXAMPLE_SUITE_NO_MAIN"] = previous)
  end
end

@testset "Example suites" begin
  runner = _sandbox_module()
  Base.include(runner, joinpath(_EXAMPLES_DIR, "internal", "example_suite_runner.jl"))

  @testset "CLI parsing" begin
    specs = [runner.ExampleSpec(name = "alpha", file = "alpha.jl", purpose = "p"), runner.ExampleSpec(name = "beta", file = "beta.jl", purpose = "p", heavy = true)]
    notes = String[]

    help_result = redirect_stdout(devnull) do
      runner.parse_example_suite_cli("x_suite", "X suite", specs, notes, ["--help"])
    end
    @test help_result === nothing

    opt = runner.parse_example_suite_cli("x_suite", "X suite", specs, notes, String[])
    @test opt["continue-on-error"] === true
    @test opt["include-heavy"] === false
    @test endswith(normpath(opt["output-dir"]), joinpath("examples", "_out", "x_suite"))

    opt = runner.parse_example_suite_cli("x_suite", "X suite", specs, notes, ["--include-heavy", "--strict", "--quiet", "--only=alpha", "--timeout=42", "--continue-on-error=false"])
    @test opt["include-heavy"] === true
    @test opt["strict"] === true
    @test opt["quiet"] === true
    @test opt["only"] == "alpha"
    @test opt["timeout"] == 42
    @test opt["continue-on-error"] === false

    @test_throws ArgumentError runner.parse_example_suite_cli("x_suite", "X suite", specs, notes, ["--nonsense=1"])
    @test_throws ArgumentError runner.parse_example_suite_cli("x_suite", "X suite", specs, notes, ["bare-arg"])
    @test_throws ArgumentError runner._parse_bool("maybe")
    @test runner._parse_bool("Yes") === true
    @test runner._parse_bool("0") === false
  end

  @testset "skip and filter logic" begin
    plain = runner.ExampleSpec(name = "plain", file = "plain.jl", purpose = "p")
    heavy = runner.ExampleSpec(name = "heavy", file = "heavy.jl", purpose = "p", heavy = true)
    optional = runner.ExampleSpec(name = "optional", file = "optional.jl", purpose = "p", optional = true)
    missing_dep = runner.ExampleSpec(name = "missing_dep", file = "missing_dep.jl", purpose = "p", optional = true, requires_package = "DefinitelyNotInstalledPackage123")
    needs_config = runner.ExampleSpec(name = "needs_config", file = "needs_config.jl", purpose = "p", requires_config = true)
    specs = [plain, heavy, optional, missing_dep, needs_config]

    opt = runner._default_suite_options("x_suite")
    no_filter = nothing
    no_skip = Set{String}()

    @test runner._spec_skip_status(plain, opt, no_filter, no_skip) === nothing
    @test runner._spec_skip_status(heavy, opt, no_filter, no_skip)[1] == "skipped_heavy"
    @test runner._spec_skip_status(optional, opt, no_filter, no_skip)[1] == "skipped_optional"

    opt["include-heavy"] = true
    opt["include-optional"] = true
    @test runner._spec_skip_status(heavy, opt, no_filter, no_skip) === nothing
    @test runner._spec_skip_status(optional, opt, no_filter, no_skip) === nothing
    @test runner._spec_skip_status(missing_dep, opt, no_filter, no_skip)[1] == "skipped_missing_dependency"
    @test runner._package_available("Test") === true

    config_expected = runner._user_config_available() ? nothing : "skipped_missing_config"
    skip = runner._spec_skip_status(needs_config, opt, no_filter, no_skip)
    @test (skip === nothing ? nothing : skip[1]) == config_expected

    only_names = runner._split_name_list("heavy, plain", specs, "only")
    @test only_names == Set(["heavy", "plain"])
    @test runner._spec_skip_status(optional, opt, only_names, no_skip)[1] == "skipped_by_filter"
    @test runner._spec_skip_status(plain, opt, only_names, Set(["plain"]))[1] == "skipped_by_filter"
    @test_throws ArgumentError runner._split_name_list("unknown_name", specs, "only")
  end

  @testset "CSV and Markdown escaping" begin
    @test runner._csv_cell("plain") == "plain"
    @test runner._csv_cell("a,b") == "\"a,b\""
    @test runner._csv_cell("say \"hi\"") == "\"say \"\"hi\"\"\""
    @test runner._csv_cell(missing) == ""
    @test runner._md("a|b\nc") == "a\\|b c"
    @test runner._is_failure_status("failed")
    @test runner._is_failure_status("timeout")
    @test !runner._is_failure_status("ok")
    @test !runner._is_failure_status("skipped_heavy")
  end

  @testset "registry integrity" begin
    all_names = String[]
    excluded = ("current_iteration_start.jl", "dtf_validation_report.jl", "for002_matpower_metadata_validation.jl", "run_val_dtf_suite.jl")
    for script in _SUITE_SCRIPTS
      m = _include_with_no_main(joinpath(_EXAMPLES_DIR, script))
      specs = Base.invokelatest(getglobal, m, :SUITE_SPECS)
      @test !isempty(specs)
      for spec in specs
        @test isfile(joinpath(_EXAMPLES_DIR, spec.file))
        @test !isempty(spec.purpose)
        @test spec.timeout_s > 0
        @test !(basename(spec.file) in excluded)
        push!(all_names, spec.name)
      end
    end
    @test allunique(all_names)
  end

  @testset "--help smoke test" begin
    for script in _SUITE_SCRIPTS
      cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $(joinpath(_EXAMPLES_DIR, script)) --help`
      proc = run(pipeline(Cmd(cmd; ignorestatus = true); stdout = devnull, stderr = devnull))
      @test success(proc)
    end
  end
end

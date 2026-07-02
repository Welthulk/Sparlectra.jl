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

using Test

function _csv_data_rows(path::AbstractString)
  lines = readlines(path)
  return max(length(lines) - 1, 0)
end

function _finite_metric_values(path::AbstractString)
  lines = readlines(path)
  @test length(lines) == 2
  headers = split(lines[1], ',')
  values = split(lines[2], ',')
  numeric = Float64[]
  for (h, v) in zip(headers, values)
    h in ("converged",) && continue
    push!(numeric, parse(Float64, v))
  end
  return numeric
end

function run_dtf_for002_validation_example_tests()
  @testset "native DTF/FOR002 validation example" begin
    repo = normpath(joinpath(@__DIR__, "..", ".."))
    script = joinpath(repo, "examples", "validate_dtf_for002_testnetz13.jl")
    source = read(script, String)
    @test occursin("DTFImporter.read_dtf", source)
    @test occursin("DTFImporter.build_net", source)
    @test !occursin("createNetFromMatPower", source)
    @test !occursin("build_for001_testnetz13", source)
    extended_runner = read(joinpath(repo, "test", "extended", "runtests_extended.jl"), String)
    @test !occursin("test/runtests.jl", extended_runner)
    @test !occursin("runtests.jl", replace(extended_runner, "runtests_extended.jl" => ""))

    outdir = mktempdir()
    Base.include(Main, joinpath(repo, "examples", "validate_dtf_for002_testnetz13.jl"))
    runner = Base.invokelatest(getfield, Main, :run_validation)
    result = Base.invokelatest(
      runner,
      [
        "--dtf-file=$(joinpath(repo, "test", "fixtures", "dtf", "FOR001.DAT"))",
        "--for002-file=$(joinpath(repo, "examples", "FOR002.DAT"))",
        "--output-dir=$outdir",
        "--write-csv=true",
        "--write-markdown=true",
        "--quiet=true",
      ],
    )
    @test result.converged
    @test nameof(typeof(result)) == :DTFFor002ValidationResult
    result_fields = fieldnames(typeof(result))
    for field in (:case, :net, :residuals, :generator_rows, :kcl_rows, :q_diagnostics)
      @test field ∉ result_fields
      @test !hasproperty(result, field)
    end

    detailed = Base.invokelatest(
      runner,
      [
        "--dtf-file=$(joinpath(repo, "test", "fixtures", "dtf", "FOR001.DAT"))",
        "--for002-file=$(joinpath(repo, "examples", "FOR002.DAT"))",
        "--output-dir=$(mktempdir())",
        "--write-csv=true",
        "--write-markdown=true",
        "--quiet=true",
      ];
      return_details = true,
    )
    @test detailed.converged
    @test hasproperty(detailed, :generator_rows)
    @test hasproperty(detailed, :kcl_rows)
    @test hasproperty(detailed, :residuals)
    @test hasproperty(detailed, :q_diagnostics)

    expected = [
      "dtf_for002_validation_summary.md",
      "dtf_import_summary.csv",
      "dtf_bus_comparison.csv",
      "dtf_generator_comparison.csv",
      "dtf_bus_kcl_comparison.csv",
      "dtf_q_semantics_diagnostics.csv",
      "dtf_branch_comparison.csv",
      "dtf_state_residual.csv",
      "dtf_validation_metrics.csv",
    ]
    for name in expected
      @test isfile(joinpath(outdir, name))
    end
    @test _csv_data_rows(joinpath(outdir, "dtf_bus_comparison.csv")) == 13
    @test _csv_data_rows(joinpath(outdir, "dtf_generator_comparison.csv")) >= count(r -> r.has_generator, detailed.generator_rows)
    @test _csv_data_rows(joinpath(outdir, "dtf_bus_kcl_comparison.csv")) == 13
    @test _csv_data_rows(joinpath(outdir, "dtf_q_semantics_diagnostics.csv")) > 0
    @test _csv_data_rows(joinpath(outdir, "dtf_branch_comparison.csv")) == 27
    @test all(isfinite, _finite_metric_values(joinpath(outdir, "dtf_validation_metrics.csv")))
    summary = read(joinpath(outdir, "dtf_for002_validation_summary.md"), String)
    @test occursin("Sparlectra slack bus: BSTADTS1", summary)
    @test occursin("## What are state residuals?", summary)

    cli_outdir = mktempdir()
    cli_output = read(
      `julia --project=$repo $script --dtf-file=$(joinpath(repo, "test", "fixtures", "dtf", "FOR001.DAT")) --for002-file=$(joinpath(repo, "examples", "FOR002.DAT")) --output-dir=$cli_outdir --write-csv=true --write-markdown=true`,
      String,
    )
    @test occursin("Native DTF/FOR002 validation", cli_output)
    @test occursin("Output directory: $cli_outdir", cli_output)
    @test occursin("Converged: true", cli_output)
    @test occursin("Iterations:", cli_output)
    @test occursin("Final mismatch:", cli_output)
    @test occursin("Max |dV|:", cli_output)
    @test occursin("Max branch |dP|:", cli_output)
    @test occursin("Written files:", cli_output)
    @test occursin("  - dtf_validation_metrics.csv", cli_output)
    @test !occursin("DTFCase(", cli_output)
    @test !occursin("Net:", cli_output)
    @test !occursin("generator_rows =", cli_output)
    @test !occursin("kcl_rows =", cli_output)
    for name in expected
      @test isfile(joinpath(cli_outdir, name))
    end
  end
end

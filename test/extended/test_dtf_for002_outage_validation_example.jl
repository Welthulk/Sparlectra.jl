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

function run_dtf_for002_outage_validation_example_tests()
  @testset "native DTF/FOR002 outage validation example" begin
    repo = normpath(joinpath(@__DIR__, "..", ".."))
    script = joinpath(repo, "examples", "validate_dtf_for002_outages_testnetz13.jl")
    source = read(script, String)
    @test occursin("DTFImporter.read_dtf", source)
    @test occursin("DTFImporter.build_net", source)
    @test !occursin("createNetFromMatPower", source)
    @test !occursin("build_for001_testnetz13", source)
    extended_runner = read(joinpath(repo, "test", "extended", "runtests_extended.jl"), String)
    @test !occursin("test/runtests.jl", extended_runner)
    @test !occursin("runtests.jl", replace(extended_runner, "runtests_extended.jl" => ""))

    dtf_file = joinpath(repo, "data", "DTF", "FOR001.DAT")
    for002_file = joinpath(repo, "data", "DTF", "FOR002.DAT")
    if !(isfile(dtf_file) && isfile(for002_file))
       @info "Skipping full DTF/FOR002 outage validation example; local external FOR001/FOR002 files are absent." dtf_file for002_file
      @test_skip "local FOR001/FOR002 validation data not available"
      return
    end

    outdir = mktempdir()
    Base.include(Main, script)
    runner = Base.invokelatest(getfield, Main, :run_validation)
    detailed = Base.invokelatest(runner, ["--dtf-file=$dtf_file", "--for002-file=$for002_file", "--output-dir=$outdir", "--write-csv=true", "--write-markdown=true", "--quiet=true"]; return_details=true)
    @test hasproperty(detailed, :outage_results)
    @test length(detailed.outage_results) == 2
    expected = ["dtf_outage_validation_summary.md", "dtf_outage_matching.csv", "dtf_outage_metrics.csv", "dtf_outage_bus_comparison.csv", "dtf_outage_generator_comparison.csv", "dtf_outage_branch_comparison.csv", "dtf_outage_bus_kcl_comparison.csv", "dtf_outage_state_residual.csv"]
    for name in expected
      @test isfile(joinpath(outdir, name))
    end
    @test _csv_data_rows(joinpath(outdir, "dtf_outage_matching.csv")) == 2
    @test all(r -> r.match_status == "matched" && r.matched_branch_index !== missing, detailed.matching_rows)
    @test _csv_data_rows(joinpath(outdir, "dtf_outage_metrics.csv")) == 2
    @test _csv_data_rows(joinpath(outdir, "dtf_outage_bus_comparison.csv")) == 26
    @test _csv_data_rows(joinpath(outdir, "dtf_outage_branch_comparison.csv")) == 54
    summary = read(joinpath(outdir, "dtf_outage_validation_summary.md"), String)
    @test occursin("## What are state residuals?", summary)
    for r in detailed.metrics_rows
      @test all(isfinite, Float64[r.final_mismatch, r.max_abs_d_vm_kV, r.max_abs_d_vm_pu, r.max_abs_d_va_deg, r.max_abs_branch_d_p_MW, r.max_abs_branch_d_q_MVar])
    end

    cli_outdir = mktempdir()
    cli_output = read(`julia --project=$repo $script --dtf-file=$dtf_file --for002-file=$for002_file --output-dir=$cli_outdir --write-csv=true --write-markdown=true`, String)
    @test occursin("Native DTF/FOR002 outage validation", cli_output)
    @test occursin("Parsed outages: 2", cli_output)
    @test occursin("Matched FOR002 outage blocks: 2", cli_output)
    @test occursin("[1] L1 ALPHA S1 -> BETA1 S1", cli_output)
    @test occursin("[2] T1A BETA1 S1 -> BETA2 S1", cli_output)
    @test !occursin("DTFCase(", cli_output)
    @test !occursin("Net:", cli_output)
    @test !occursin("generator_rows =", cli_output)
    @test !occursin("kcl_rows =", cli_output)
    @test length(split(cli_output, '\n')) < 80
  end
end

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

function run_dtf_matpower_export_validation_example_tests()
  repo = dirname(dirname(@__DIR__))
  script = joinpath(repo, "examples", "internal", "dtf_validation_matpower.jl")
  @testset "DTF MATPOWER export validation example" begin
    source = read(script, String)
    @test occursin("DTFImporter.read_dtf", source)
    @test occursin("DTFImporter.build_net", source)
    @test occursin("writeMatpowerCasefile", source)
    @test !occursin("DTF -> MATPOWER", source)
    @test occursin("createNetFromMatPowerFile", source)

    dtf_file = joinpath(repo, "data", "DTF", "FOR001.DAT")
    for002_file = joinpath(repo, "data", "DTF", "FOR002.DAT")
    if !(isfile(dtf_file) && isfile(for002_file))
       @info "Skipping full DTF MATPOWER export validation example; local external FOR001/FOR002 files are absent." dtf_file for002_file
      @test_skip "local FOR001/FOR002 validation data not available"
      return
    end

    outdir = mktempdir()
    cmd = `$(Base.julia_cmd()) --project=$repo $script --dtf-file=$dtf_file --for002-file=$for002_file --output-dir=$outdir --write-csv=true --write-markdown=true --write-matpower=true --run-outages=true`
    output = read(cmd, String)
    @test occursin("Native DTF -> existing MATPOWER export validation", output)
    @test !occursin("DTFCase", output)
    @test !occursin("Net(", output)
    @test isfile(joinpath(outdir, "dtf_matpower_export_summary.md"))
    @test isfile(joinpath(outdir, "dtf_matpower_export_metrics.csv"))
    @test isfile(joinpath(outdir, "dtf_matpower_export_bus_comparison.csv"))
    @test isfile(joinpath(outdir, "dtf_matpower_export_branch_comparison.csv"))
    @test isfile(joinpath(outdir, "dtf_matpower_export_generator_comparison.csv"))
    @test isfile(joinpath(outdir, "base.m"))
    @test isfile(joinpath(outdir, "outage_1.m"))
    @test isfile(joinpath(outdir, "outage_2.m"))
    mpc = Sparlectra.MatpowerIO.read_case_m(joinpath(outdir, "base.m"); legacy_compat = false)
    @test mpc.bus_name !== nothing
    @test length(mpc.bus_name) == 13
    @test mpc.branch_name !== nothing
    @test length(mpc.branch_name) == 27
    @test mpc.branch_kind !== nothing
    @test length(mpc.branch_kind) == 27
    @test Set(mpc.branch_kind) == Set(["L", "T"])
    @test all(mpc.branch[i, 9] == 0.0 for i in axes(mpc.branch, 1) if mpc.branch_kind[i] == "L")
    @test all(mpc.branch[i, 9] != 0.0 for i in axes(mpc.branch, 1) if mpc.branch_kind[i] == "T")
    @test mpc.for001_contingencies !== nothing
    @test length(mpc.for001_contingencies) == 2
    @test "L1 ALPHA S1 -> BETA1 S1" in mpc.for001_contingencies
    @test "T1A BETA1 S1 -> BETA2 S1" in mpc.for001_contingencies
    metrics = readlines(joinpath(outdir, "dtf_matpower_export_metrics.csv"))
    @test length(metrics) == 4
    @test any(line -> startswith(line, "\"base\","), metrics)
    @test count(line -> startswith(line, "\"outage_"), metrics) == 2
    @test length(readlines(joinpath(outdir, "dtf_matpower_export_bus_comparison.csv"))) == 1 + 13*3
    @test length(readlines(joinpath(outdir, "dtf_matpower_export_branch_comparison.csv"))) == 1 + 27*3

    mod = Module(:DTFMatpowerExportValidationExampleTest)
    Core.eval(mod, :(using Sparlectra))
    Core.eval(mod, :(include(path) = Base.include($mod, path)))
    Base.include(mod, script)
    runner = Base.invokelatest(() -> getfield(getfield(mod, :MatpowerRoundtripValidation), :run_validation))
    details = Base.invokelatest(runner, ["--dtf-file=$dtf_file", "--for002-file=$for002_file", "--output-dir=$(mktempdir())", "--quiet=true", "--write-matpower=true", "--run-outages=true"]; return_details = true)
    @test hasproperty(details, :scenario_results)
    @test length(details.scenario_results) == 3
    @test all(r -> r.metric.roundtrip_converged, details.scenario_results)
    # Testnetz13's BETA1S1<->BETA2S1 and DELTA1S1<->DELTA2S1 transformers carry
    # real nonzero no-load conductance (DTF G parameter); the roundtrip must
    # preserve their Ybus contribution via equivalent bus shunts (checked
    # below), not force branch g_pu to zero.
    @test all(r -> r.metric.count_nonzero_branch_g_pu_native == 5, details.scenario_results)
    @test all(r -> r.metric.max_abs_branch_g_pu_native > 0.0, details.scenario_results)
    @test all(r -> isapprox(r.metric.native_total_bus_Gs, r.metric.roundtrip_total_bus_Gs; atol = 1e-12), details.scenario_results)
    @test all(r -> length(r.roundtrip.linesAC) == 22, details.scenario_results)
    @test all(r -> length(r.roundtrip.trafos) == 5, details.scenario_results)

    # Regression: with tap_changer_model=impedance_correction, the exporter's
    # mpc.sparlectra.tap_changer_model marker must prevent the reimport from
    # reapplying calcTapCorrectedRX on top of the already-corrected r/x (see
    # matpower_import.md "Tap-impedance correction and reimport"). Transformer
    # branch impedances must be bit-identical between native and roundtrip.
    impedance_details = Base.invokelatest(runner, ["--dtf-file=$dtf_file", "--for002-file=$for002_file", "--output-dir=$(mktempdir())", "--quiet=true", "--write-matpower=true", "--run-outages=false", "--tap-changer-model=impedance_correction"]; return_details = true)
    @test all(r -> r.metric.roundtrip_converged, impedance_details.scenario_results)
    transformer_rows = [row for row in impedance_details.branch_rows if row.branch_kind_native == "T"]
    @test !isempty(transformer_rows)
    @test all(row -> row.x_pu_native == row.x_pu_roundtrip, transformer_rows)
    @test all(row -> row.r_pu_native == row.r_pu_roundtrip, transformer_rows)
  end
end

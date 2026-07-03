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
using Sparlectra

function run_dtf_api_webui_integration_tests()
  @testset "DTF/FOR001 API and Web UI integration" begin
    mktempdir() do tmp
      dtf = joinpath(@__DIR__, "..", "fixtures", "dtf", "FOR001.DAT")
      for002 = joinpath(@__DIR__, "..", "..", "examples", "FOR002.DAT")
      direct_output = joinpath(@__DIR__, "..", "..", "examples", "_out", "dtf_api_direct_smoke")
      isdir(direct_output) && rm(direct_output; recursive = true, force = true)
      result = Sparlectra.run_sparlectra_api(
        casefile = dtf,
        output_dir = direct_output,
        case_format = :dtf_for001,
        for002_reference_file = for002,
        detailed_result_csv = true,
        run_dtf_outages = true,
        dtf_outage_selection_mode = :selected,
        dtf_outage_selection = ["1"],
        write_outage_artifacts = true,
        matpower_export_requested = true,
        config_overrides = Dict("benchmark.enabled" => false, "output.logfile_results" => "classic"),
      )
      @test result.status == :succeeded
      @test result.metadata["input_format_detected"] == "dtf_for001"
      @test result.metadata["native_dtf_import_used"] == true
      @test result.metadata["dtf_bus_count"] == 13
      @test result.metadata["dtf_branch_count"] == 27
      @test result.metadata["dtf_outage_count"] == 2
      @test length(result.metadata["dtf_outage_results"]) >= 1
      artifact_names = Set(a.name for a in result.artifacts)
      bus_csv = read(joinpath(direct_output, "bus_voltages_complex.csv"), String)
      branch_csv = read(joinpath(direct_output, "branch_flows.csv"), String)
      bus_header = split(first(eachline(IOBuffer(bus_csv))), ',')
      branch_header = split(first(eachline(IOBuffer(branch_csv))), ',')
      @test "original_bus_name" in bus_header
      @test all(name -> name in branch_header, ("branch_name", "original_branch_name", "from_bus_name", "to_bus_name", "original_from_bus_name", "original_to_bus_name", "branch_kind"))
      @test occursin("ALPHA", bus_csv)
      @test occursin("ALPHA", branch_csv)
      @test occursin("L1 ALPHA S1 -> BETA1 S1", branch_csv)
      @test "dtf_import_summary.md" in artifact_names
      @test "dtf_import_summary.csv" in artifact_names
      @test "dtf_native_matpower_export.m" in artifact_names
      @test "dtf_for002_base_comparison.md" in artifact_names
      @test any(name -> occursin("dtf_outage_1_metrics.csv", name), artifact_names)
      @test any(name -> endswith(name, ".csv"), artifact_names)

      outage = Sparlectra.run_sparlectra_api(
        casefile = dtf,
        output_dir = joinpath(tmp, "dtf-outage"),
        case_format = :dtf_for001,
        run_dtf_outages = true,
        dtf_outage_selection_mode = :selected,
        dtf_outage_selection = ["1"],
        config_overrides = Dict("benchmark.enabled" => false, "output.logfile_results" => "classic"),
      )
      @test outage.status == :succeeded
      @test length(outage.metadata["dtf_outage_results"]) == 1
      @test outage.metadata["dtf_outage_results"][1]["converged"] == true
      @test any(a -> occursin("dtf_outage_1_metrics.csv", a.name), outage.artifacts)

      form_html = Sparlectra.render_powerflow_form(output_root = tmp, case_directory = dirname(dtf), selected_casefile = basename(dtf))
      @test occursin("Input format", form_html)
      @test occursin("DTF/FOR001 diagnostics (experimental/internal)", form_html)
      @test occursin("MATPOWER", form_html)
      @test occursin("option value=\"auto\">Auto", form_html)
      @test !occursin("New: full DTF/FOR001 support", form_html)

      request = Sparlectra.powerflow_webui_request(Dict(
        "casefile" => basename(dtf),
        "casefile_manual" => dtf,
        "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
        "case_format" => "dtf_for001",
        "for002_reference_file" => for002,
        "dtf_outage_selection_mode" => "selected",
        "dtf_outage_selection" => "1",
        "write_outage_artifacts" => "true",
        "write_outage_matpower_exports" => "true",
        "matpower_export_requested" => "true",
        "detailed_result_csv" => "true",
      ); default_output_root = tmp)
      @test request["case_format"] == "dtf_for001"
      @test request["run_dtf_outages"] == true
      @test request["dtf_outage_selection"] == ["1"]
      @test request["for002_reference_file"] == for002
      @test request["write_outage_artifacts"] == true
      @test request["write_outage_matpower_exports"] == true
      @test request["matpower_export_requested"] == true
      @test request["detailed_result_csv"] == true

      dc_case = joinpath(tmp, "FOR001_HVDC.DAT")
      write(dc_case, read(dtf, String) * "\nHVDC unsupported marker\n")
      dc_result = Sparlectra.run_sparlectra_api(
        casefile = dc_case,
        output_dir = joinpath(tmp, "dtf-dc"),
        case_format = :dtf_for001,
        config_overrides = Dict("benchmark.enabled" => false, "output.logfile_results" => "classic"),
      )
      @test dc_result.status == :failed
      @test dc_result.reason == "unsupported_dtf_dc_line"
      @test dc_result.metadata["unsupported_dcline_status"] == "unsupported_dtf_dc_line"
      @test occursin("DC lines are currently not supported by the native DTF/MATPOWER power-flow path.", dc_result.message)
    end
  end
end

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

using Sparlectra
using Test

include("test_webui_support.jl")

function run_webui_fast_tests()
  @testset "Web UI fast parsing, selectors, and rendering" begin
    @test Sparlectra._webui_parse_bool("on") === true
    @test Sparlectra._webui_parse_bool("false") === false
    @test Sparlectra.route_sparlectra_webui("GET", "/powerflow/artifact/..%2Fbad/result.json").status in (400, 404)

    mktempdir() do root
      case_directory = joinpath(root, "cases")
      mkpath(case_directory)
      write(joinpath(case_directory, "case14.m"), "function mpc = case14\nend\n")
      write(joinpath(case_directory, "FOR001.DAT"), _dtf_network_fixture())
      write(joinpath(case_directory, "FOR001_OUTAGES.DAT"), _dtf_network_with_outage_fixture())
      write(joinpath(case_directory, "OUTAGE.DAT"), _dtf_outage_fixture())
      write(joinpath(case_directory, "FOR002.DAT"), _dtf_reference_fixture())
      write(joinpath(case_directory, "UNKNOWN.DAT"), "plain unsupported data\n")
      @test Sparlectra._webui_classify_dat_content(joinpath(case_directory, "FOR001.DAT")) === :dtf_network_case
      @test Sparlectra._webui_classify_dat_content(joinpath(case_directory, "FOR001_OUTAGES.DAT")) === :dtf_network_case_with_outages
      @test Sparlectra._webui_classify_dat_content(joinpath(case_directory, "OUTAGE.DAT")) === :dtf_outage_file
      @test Sparlectra._webui_classify_dat_content(joinpath(case_directory, "FOR002.DAT")) === :dtf_outage_or_reference
      @test Sparlectra._webui_classify_dat_content(joinpath(case_directory, "UNKNOWN.DAT")) === :unknown_dat
      primary = Sparlectra._webui_casefile_options_in_directory(case_directory)
      reference = Sparlectra._webui_for002_reference_options_in_directory(case_directory)
      @test "FOR001.DAT" in primary
      @test "FOR001_OUTAGES.DAT" in primary
      @test !("OUTAGE.DAT" in primary)
      @test !("UNKNOWN.DAT" in primary)
      @test reference == ["FOR002.DAT"]

      active_html = Sparlectra.render_powerflow_result(Dict("run_id" => "active", "status" => "running", "elapsed_seconds" => 1.25))
      @test occursin("Run status", active_html)
      @test occursin("Elapsed time", active_html)
      @test !occursin("Solver time", active_html)
      @test !occursin("Wall time", active_html)

      terminal = Dict("run_id" => "ok", "status" => "completed", "solver_elapsed_s" => 0.125, "service_phase_timings" => [Dict("phase" => "total_service", "elapsed_seconds" => 1.5)])
      terminal_html = Sparlectra.render_powerflow_result(terminal)
      @test occursin("Solver time", terminal_html)
      @test occursin("Total time", terminal_html)
      @test !occursin("Elapsed time", terminal_html)
      @test !occursin("Wall time", terminal_html)
      @test !occursin("Solver time: n/a", terminal_html)

      failed_solver_html = Sparlectra.render_powerflow_result(Dict("run_id" => "failed", "status" => "failed", "metadata" => Dict("solver_elapsed_s" => 0.25), "service_phase_timings" => [Dict("phase" => "total_service", "elapsed_seconds" => 2.0)]))
      @test occursin("Solver time", failed_solver_html)
      presolver_html = Sparlectra.render_powerflow_result(Dict("run_id" => "bad", "status" => "failed", "service_phase_timings" => [Dict("phase" => "total_service", "elapsed_seconds" => 0.1)]))
      @test !occursin("Solver time", presolver_html)
      @test occursin("Total time", presolver_html)
      legacy_html = Sparlectra.render_powerflow_result(Dict("run_id" => "old", "status" => "completed", "elapsed_seconds" => 3.0))
      @test occursin("Total time", legacy_html)
    end
  end
end

run_webui_tests() = run_webui_fast_tests()

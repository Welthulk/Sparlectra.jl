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

function _fake_api_result(; status::Symbol, success::Bool, converged::Bool, output_dir::String, mode::String)
  artifacts = [
    Sparlectra.SparlectraApiArtifact("q_limit.log", :q_limit_log, joinpath(output_dir, "q_limit.log"), "text/plain", true, 12, "Q-limit diagnostic log"),
    Sparlectra.SparlectraApiArtifact("q_limit_events.csv", :q_limit_events, joinpath(output_dir, "q_limit_events.csv"), "text/csv", true, 12, "Q-limit events"),
  ]
  metadata = Dict{String,Any}(
    "run_status" => success ? "completed" : "completed_nonconverged",
    "numerical_status" => converged ? "converged" : "not_converged",
    "solver_elapsed_s" => 0.25,
    "q_limit_pv_to_pq_events" => 2,
    "q_limit_pv_to_pq_buses" => "2;5",
    "q_limit_active_set_events" => mode == "active_set" ? 2 : 0,
    "q_limit_reenable_events" => 1,
    "q_limit_classic_outer_loop_passes" => startswith(mode, "classic_") ? 1 : 0,
  )
  return Sparlectra.SparlectraApiResult("fake-$(mode)", "1", status, success, converged, true, 7, success ? 1.0e-10 : 0.01, success ? "none" : "nr_mismatch_not_converged", "fake result", "case.m", "config.yaml", output_dir, joinpath(output_dir, "run.log"), joinpath(output_dir, "result.json"), artifacts, Dict{String,Any}[], metadata, nothing)
end

function run_qlimit_large_case_comparison_tests()
  @testset "Q-limit large-case comparison orchestration" begin
    mktempdir() do dir
      calls = String[]
      resolver(case; outdir) = case == "missing.m" ?
                               Dict{String,Any}("available" => false, "path" => "", "status" => "skipped", "source" => "unavailable", "reason" => "download failed in test") :
                               Dict{String,Any}("available" => true, "path" => joinpath(outdir, case), "status" => case == "present.m" ? "cached" : "downloaded", "source" => case == "present.m" ? "cache" : "download")
      runner(; casefile, config_file, output_dir, config_overrides, kwargs...) = begin
        mode = String(config_overrides["power_flow.qlimits.enforcement_mode"])
        voltage = String(config_overrides["power_flow.start_mode.voltage_mode"])
        angle = String(config_overrides["power_flow.start_mode.angle_mode"])
        push!(calls, string(basename(casefile), ":", voltage, ":", angle, ":", mode, ":", config_overrides["output.logfile_results"]))
        mode == "active_set" && voltage == "profile_blend" && error("profile blend unavailable in test")
        mkpath(output_dir)
        return _fake_api_result(status = mode == "classic_one_at_a_time" ? :not_converged : :succeeded, success = mode != "classic_one_at_a_time", converged = mode != "classic_one_at_a_time", output_dir = output_dir, mode = mode)
      end
      progress = IOBuffer()
      result = Sparlectra.compare_qlimit_large_case_modes(
        cases = ["present.m", "downloaded.m", "missing.m"],
        output_root = joinpath(dir, "out"),
        case_cache_dir = joinpath(dir, "cache"),
        resolver = resolver,
        runner = runner,
        io = progress,
      )
      @test isfile(result.csv_path)
      @test isfile(result.json_path)
      @test length(result.rows) == 18
      @test count(row -> row["case_resolve_status"] == "skipped", result.rows) == 6
      @test count(row -> row["run_status"] == "failed", result.rows) == 2
      @test count(row -> row["numerical_status"] == "converged", result.rows) == 6
      @test count(row -> row["numerical_status"] == "not_converged", result.rows) == 4
      @test all(row -> haskey(row, "start_profile") && haskey(row, "case_path") && haskey(row, "api_completed"), result.rows)
      @test Set((row["case"], row["start_profile"], row["mode"]) for row in result.rows) == Set((case, String(profile), String(mode)) for case in ["present.m", "downloaded.m", "missing.m"] for profile in Sparlectra.QLIMIT_LARGE_CASE_START_PROFILES for mode in Sparlectra.QLIMIT_LARGE_CASE_MODES)
      @test all(endswith(call, ":compact") for call in calls)
      skipped = filter(row -> row["case"] == "missing.m", result.rows)
      @test length(skipped) == 6
      @test all(row -> occursin("download failed in test", row["error_message"]) && row["run_status"] == "skipped", skipped)
      text = String(take!(progress))
      @test occursin("Q-limit large-case mode comparison", text)
      @test occursin("resolving present.m", text)
      @test occursin("present.m | classic_start | active_set", text)
      @test occursin("summary csv", text)
      csv_text = read(result.csv_path, String)
      @test startswith(csv_text, "case,case_path,case_resolve_status")
      @test occursin("start_profile", csv_text)
      @test occursin("robust_dc_start", csv_text)
      @test occursin("profile blend unavailable in test", csv_text)
      json_text = read(result.json_path, String)
      @test occursin("\"start_profile\"", json_text)
    end
  end
end

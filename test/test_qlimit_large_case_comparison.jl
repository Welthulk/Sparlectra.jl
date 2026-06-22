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
  ]
  metadata = Dict{String,Any}(
    "run_status" => success ? "completed" : "completed_nonconverged",
    "numerical_status" => converged ? "converged" : "not_converged",
    "q_limit_pv_to_pq_events" => 2,
    "q_limit_active_set_events" => mode == "active_set" ? 2 : 0,
    "q_limit_classic_outer_loop_passes" => startswith(mode, "classic_") ? 1 : 0,
  )
  return Sparlectra.SparlectraApiResult("fake-$(mode)", "1", status, success, converged, true, 7, success ? 1.0e-10 : 0.01, success ? "none" : "nr_mismatch_not_converged", "fake result", "case.m", "config.yaml", output_dir, joinpath(output_dir, "run.log"), joinpath(output_dir, "result.json"), artifacts, Dict{String,Any}[], metadata, nothing)
end

function run_qlimit_large_case_comparison_tests()
  @testset "Q-limit large-case comparison orchestration" begin
    mktempdir() do dir
      calls = String[]
      resolver(case; outdir) = case == "missing.m" ?
                               Dict{String,Any}("available" => false, "path" => "", "source" => "unavailable", "reason" => "download failed in test") :
                               Dict{String,Any}("available" => true, "path" => joinpath(outdir, case), "source" => case == "present.m" ? "cache" : "download")
      runner(; casefile, config_file, output_dir, config_overrides, kwargs...) = begin
        mode = String(config_overrides["power_flow.qlimits.enforcement_mode"])
        push!(calls, string(basename(casefile), ":", mode))
        mkpath(output_dir)
        return _fake_api_result(status = mode == "classic_one_at_a_time" ? :not_converged : :succeeded, success = mode != "classic_one_at_a_time", converged = mode != "classic_one_at_a_time", output_dir = output_dir, mode = mode)
      end
      result = Sparlectra.compare_qlimit_large_case_modes(
        cases = ["present.m", "downloaded.m", "missing.m"],
        output_root = joinpath(dir, "out"),
        case_cache_dir = joinpath(dir, "cache"),
        resolver = resolver,
        runner = runner,
        io = IOBuffer(),
      )
      @test isfile(result.csv_path)
      @test isfile(result.json_path)
      @test length(result.rows) == 7
      @test count(row -> row["case_status"] == "skipped", result.rows) == 1
      @test count(row -> row["numerical_status"] == "converged", result.rows) == 4
      @test count(row -> row["numerical_status"] == "not_converged", result.rows) == 2
      @test Set(calls) == Set([
        "present.m:active_set", "present.m:classic_simultaneous", "present.m:classic_one_at_a_time",
        "downloaded.m:active_set", "downloaded.m:classic_simultaneous", "downloaded.m:classic_one_at_a_time",
      ])
      skipped = only(filter(row -> row["case"] == "missing.m", result.rows))
      @test skipped["reason"] == "download failed in test"
      @test skipped["api_status"] == "skipped"
    end
  end
end

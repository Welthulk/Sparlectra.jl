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

include("test_api_support.jl")

function run_api_fast_tests()
  @testset "API fast smoke and timing contracts" begin
    @testset "profiling wrapper records successes and exceptions" begin
      profile = Dict{Symbol,Any}(:enabled => true)
      @test Sparlectra._perf_profile_time!(profile, :solver_total) do
        sleep(0.001)
        42
      end == 42
      @test Sparlectra._solver_elapsed_from_profile(profile) > 0.0

      thrown_profile = Dict{Symbol,Any}(:enabled => true)
      @test_throws ErrorException Sparlectra._perf_profile_time!(thrown_profile, :solver_total) do
        sleep(0.001)
        error("intentional solver failure")
      end
      @test Sparlectra._solver_elapsed_from_profile(thrown_profile) > 0.0

      nested_profile = Dict{Symbol,Any}(:enabled => true)
      @test_throws ErrorException Sparlectra._perf_profile_time!(nested_profile, :solver_total) do
        Sparlectra._perf_profile_time!(nested_profile, :newton_step_linear_solve) do
          sleep(0.001)
          error("nested failure")
        end
      end
      solver = Sparlectra._solver_elapsed_from_profile(nested_profile)
      nested = nested_profile[:timings][:newton_step_linear_solve].elapsed_s
      @test solver >= nested
      @test nested_profile[:timings][:solver_total].calls == 1
    end

    mktempdir() do tmpdir
      casefile = _write_api_test_case(joinpath(tmpdir, "case_api.m"))
      template = joinpath(tmpdir, "configuration.yaml")
      write(template, "power_flow:\n  max_iter: 20\noutput:\n  logfile_results: compact\n")
      successful_transport = Sparlectra._api_result(run_id = "synthetic-success", status = :completed, success = true, reason = "none", message = "ok", casefile = casefile, config_file = template, output_dir = joinpath(tmpdir, "success"), logfile = joinpath(tmpdir, "success", "run.log"), result_file = joinpath(tmpdir, "success", "result.json"), raw_result = (solver_elapsed_s = 0.002,))
      @test Sparlectra.to_dict(successful_transport)["solver_elapsed_s"] > 0.0

      presolver_transport = Sparlectra._api_result(run_id = "synthetic-presolver", status = :failed, success = false, reason = "invalid_case", message = "missing case", casefile = joinpath(tmpdir, "missing.m"), config_file = template, output_dir = joinpath(tmpdir, "missing"), logfile = nothing, result_file = nothing)
      @test !haskey(Sparlectra.to_dict(presolver_transport), "solver_elapsed_s")

      @test Sparlectra.MatpowerIO.matpower_dcline_diagnostics((; dcline = [1 2 1 10 9 0 0 1 1 0 100]))["matpower_dcline_active_count"] == 1

      mkpath(joinpath(tmpdir, "synthetic-failure"))
      failed_transport = Sparlectra._api_failure("execution_error", "numerical failure"; casefile = casefile, config_file = template, output_dir = joinpath(tmpdir, "synthetic-failure"), logfile = joinpath(tmpdir, "synthetic-failure", "run.log"), result_file = joinpath(tmpdir, "synthetic-failure", "result.json"), metadata = Dict("solver_elapsed_s" => 0.001))
      @test Sparlectra.to_dict(failed_transport)["solver_elapsed_s"] > 0.0
    end
  end
end

run_api_tests() = run_api_fast_tests()

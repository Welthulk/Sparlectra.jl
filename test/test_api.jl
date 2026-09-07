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

# file: test/test_api.jl
# purpose: fast API smoke tests: version contract against Project.toml, the
#          profiling wrapper's success and exception recording, and run-index
#          handling of symlinked output roots
using Sparlectra
using Test

include("test_api_support.jl")

function run_api_fast_tests()
  @testset "API fast smoke and timing contracts" begin
    @testset "net parameters stamped exactly once per importer (task_import_direct D12)" begin
      # D12, corrected by task_import_direct: the stamping happens exactly
      # once PER IMPORTER, at the place each importer finishes; this test
      # is the guard against pulling the four call sites back together
      cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
      repo = dirname(@__DIR__)
      stamps() = Sparlectra._NET_PARAM_STAMP_COUNT[]
      assert_one = function (path, label; requested_format = :auto)
        before = stamps()
        imported = Sparlectra.import_case(path, cfg; requested_format = requested_format)
        @test stamps() - before == 1
        @test imported.net.cooldown_iters == cfg.powerflow.qlimits.cooldown_iters
        println("      import stamp: ", label, " ran (exactly once)")
        return imported
      end
      assert_one(joinpath(repo, "data", "mpower", "warmup_casePST.m"), "matpower")
      assert_one(joinpath(repo, "data", "scf", "sp_case5.scf.json"), "scf")
      dtf = joinpath(repo, "data", "DTF", "FOR001.DAT")
      if isfile(dtf)
        assert_one(dtf, "dtf"; requested_format = :dtf_for001)
      else
        println("      import stamp: dtf SKIPPED (data/DTF/FOR001.DAT not present)")
      end
      cgmes_dir = joinpath(repo, "data", "CGMES", "cases")
      cgmes_zip = ""
      if isdir(cgmes_dir)
        zips = sort(filter(f -> endswith(f, ".zip"), readdir(cgmes_dir; join = true)); by = filesize)
        isempty(zips) || (cgmes_zip = first(zips))
      end
      if isempty(cgmes_zip)
        println("      import stamp: cgmes SKIPPED (no cached delivery under data/CGMES/cases)")
      else
        assert_one(cgmes_zip, string("cgmes (", basename(cgmes_zip), ")"); requested_format = :cgmes)
      end
    end

    # Version-independent on purpose: assert only that the precompile-baked
    # version() matches the current Project.toml, so a version bump alone
    # cannot break the test.
    @test Sparlectra.version() == Sparlectra._read_project_version()

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

    @testset "run index accepts symlinked output roots" begin
      # Two environments can reach the same physical output root under
      # different names (a Flatpak XDG_STATE_HOME symlinked onto
      # ~/.local/state). Index entries store the writer's absolute paths,
      # so validation must resolve symlinks instead of comparing lexically,
      # in both directions: entry written under the real name and read
      # through the alias, and the reverse.
      if Sys.iswindows()
        @info "run-index symlink test SKIPPED on Windows (symlink creation may need privileges)"
      else
        mktempdir() do tmpdir
          real_root = joinpath(tmpdir, "real_root")
          run_id = "11111111-2222-3333-4444-555555555555"
          mkpath(joinpath(real_root, run_id))
          write(joinpath(real_root, run_id, "result.json"), "{}")
          alias_root = joinpath(tmpdir, "alias_root")
          symlink(real_root, alias_root)
          real_entry = Dict{String,Any}("run_id" => run_id, "output_dir" => joinpath(real_root, run_id), "result_file" => joinpath(real_root, run_id, "result.json"))
          alias_entry = Dict{String,Any}("run_id" => run_id, "output_dir" => joinpath(alias_root, run_id), "result_file" => joinpath(alias_root, run_id, "result.json"))
          @test Sparlectra._indexed_run_paths(real_entry, alias_root).valid
          @test Sparlectra._indexed_run_paths(alias_entry, real_root).valid
          # a symlink must still not smuggle a foreign directory into the root
          outside = joinpath(tmpdir, "outside")
          mkpath(outside)
          write(joinpath(outside, "result.json"), "{}")
          foreign_entry = Dict{String,Any}("run_id" => run_id, "output_dir" => outside, "result_file" => joinpath(outside, "result.json"))
          @test !Sparlectra._indexed_run_paths(foreign_entry, real_root).valid
        end
      end
    end
  end
end

run_api_tests() = run_api_fast_tests()

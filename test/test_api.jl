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
  @testset "API fast smoke and timing contracts" begin (function ()
    @testset "application layer smoke: one power flow and one state estimation" begin (function ()
      # the service API is the second way to every result; this stays in the
      # pull-request gate so a change to the library cannot break it unseen
      scf5 = joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json")
      meas5 = joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.measurements.csv")
      mktempdir() do tmpdir
        pf = run_sparlectra_api(casefile = scf5, config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, output_dir = joinpath(tmpdir, "pf"))
        @test pf.success
        @test pf.converged
        @test isfile(pf.result_file)
        # the Q-V verdict travels with every run: metadata line, count, run.log
        @test pf.metadata["qv_non_physical_states"] == 0
        @test pf.metadata["qv_characteristic_line"] == "Q-V characteristic: no non-physical generator states."
        @test occursin("Q-V characteristic: no non-physical", read(pf.logfile, String))
        # and names the buses of a non-physical end point: the Zeng/Chiang
        # case with the voltage-side release switched off (margin 1.0) ends
        # with three machines at Qmax above their setpoints
        zeng = joinpath(dirname(@__DIR__), "data", "scf", "case14_zeng_p306_activeSet_A.scf.json")
        no_release = joinpath(tmpdir, "configuration.yaml")
        write(no_release, replace(read(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, String), r"reenable_v_hyst_pu:\s*[0-9.e+-]+" => "reenable_v_hyst_pu: 1.0"))
        nonphys = run_sparlectra_api(casefile = zeng, config_file = no_release, output_dir = joinpath(tmpdir, "nonphys"))
        @test nonphys.converged
        @test nonphys.metadata["qv_non_physical_states"] == 3
        @test nonphys.metadata["qv_non_physical_buses"] == "2;3;6"
        @test occursin("3 non-physical generator state(s)", read(nonphys.logfile, String))
        se = start_powerflow_run(Dict{String,Any}("casefile" => scf5, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(tmpdir, "runs"), "se_mode" => true, "measurement_file" => meas5))
        @test se["status"] == "succeeded"
        @test se["metadata"]["run_mode"] == "se"
        # the measurements.csv artifact of an SE run follows the run's CSV
        # format like the other artifacts, and the technical set given as
        # input reads under that setting as before
        se_de = start_powerflow_run(Dict{String,Any}("casefile" => scf5, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(tmpdir, "runs"), "se_mode" => true, "measurement_file" => meas5, "config_overrides" => Dict{String,Any}("output.csv_format" => "excel_de")))
        @test se_de["status"] == "succeeded"
        meas_de = readlines(joinpath(se_de["output_dir"], "measurements.csv"))
        @test meas_de[1] == "# sparlectra-measurements v1"
        header_at = findfirst(l -> startswith(l, "type"), meas_de)
        @test header_at !== nothing && startswith(meas_de[header_at], "type;bus;")
        @test any(l -> occursin(r";\d+,\d+;", l), meas_de)
        # the artifact still reads back as a measurement set
        probe = importSCF(scf5)
        empty!(probe.measurements)
        @test readMeasurementsCSV!(probe; file = joinpath(se_de["output_dir"], "measurements.csv")).total == length(readlines(meas5)) - count(l -> startswith(l, "#") || startswith(l, "type"), readlines(meas5))
      end
    end)() end

    @testset "net parameters stamped exactly once per importer" begin (function ()
      # the stamping happens exactly
      # once PER IMPORTER, at the place each importer finishes; this test
      # is the guard against pulling the four call sites back together
      # (the CGMES leg runs on the checked-in delivery under
      # data/cgmes_demo).
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
      assert_one(cgmes_fixture_dir("sp_casePST"), "cgmes (sp_casePST fixture)"; requested_format = :cgmes)
    end)() end

    # Version-independent on purpose: assert only that the precompile-baked
    # version() matches the current Project.toml, so a version bump alone
    # cannot break the test.
    @test Sparlectra.version() == Sparlectra._read_project_version()

    @testset "profiling wrapper records successes and exceptions" begin (function ()
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
    end)() end

    mktempdir() do tmpdir
      casefile = _write_api_test_case(joinpath(tmpdir, "case_api.m"))
      template = joinpath(tmpdir, "configuration.yaml")
      write(template, "power_flow:\n  max_iter: 20\noutput:\n  logfile_results: compact\n")
      successful_transport = SparlectraApp._api_result(run_id = "synthetic-success", status = :completed, success = true, reason = "none", message = "ok", casefile = casefile, config_file = template, output_dir = joinpath(tmpdir, "success"), logfile = joinpath(tmpdir, "success", "run.log"), result_file = joinpath(tmpdir, "success", "result.json"), raw_result = (solver_elapsed_s = 0.002,))
      @test SparlectraApp.to_dict(successful_transport)["solver_elapsed_s"] > 0.0

      presolver_transport = SparlectraApp._api_result(run_id = "synthetic-presolver", status = :failed, success = false, reason = "invalid_case", message = "missing case", casefile = joinpath(tmpdir, "missing.m"), config_file = template, output_dir = joinpath(tmpdir, "missing"), logfile = nothing, result_file = nothing)
      @test !haskey(SparlectraApp.to_dict(presolver_transport), "solver_elapsed_s")

      @test Sparlectra.MatpowerIO.matpower_dcline_diagnostics((; dcline = [1 2 1 10 9 0 0 1 1 0 100]))["matpower_dcline_active_count"] == 1

      mkpath(joinpath(tmpdir, "synthetic-failure"))
      failed_transport = SparlectraApp._api_failure("execution_error", "numerical failure"; casefile = casefile, config_file = template, output_dir = joinpath(tmpdir, "synthetic-failure"), logfile = joinpath(tmpdir, "synthetic-failure", "run.log"), result_file = joinpath(tmpdir, "synthetic-failure", "result.json"), metadata = Dict("solver_elapsed_s" => 0.001))
      @test SparlectraApp.to_dict(failed_transport)["solver_elapsed_s"] > 0.0
    end

    @testset "run index accepts symlinked output roots" begin (function ()
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
          @test SparlectraApp._indexed_run_paths(real_entry, alias_root).valid
          @test SparlectraApp._indexed_run_paths(alias_entry, real_root).valid
          # a symlink must still not smuggle a foreign directory into the root
          outside = joinpath(tmpdir, "outside")
          mkpath(outside)
          write(joinpath(outside, "result.json"), "{}")
          foreign_entry = Dict{String,Any}("run_id" => run_id, "output_dir" => outside, "result_file" => joinpath(outside, "result.json"))
          @test !SparlectraApp._indexed_run_paths(foreign_entry, real_root).valid
        end
      end
    end)() end
  end)() end
end

run_api_tests() = run_api_fast_tests()

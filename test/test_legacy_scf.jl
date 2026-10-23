# Copyright 2023-2026 Udo Schmitz
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

# file: test/test_legacy_scf.jl
# purpose: legacy profile: the FACTS base-impedance guard of the SCF writer
#          and the Web UI save-as route (issue #378). The helpers they use
#          are defined in test_scf.jl, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_scf_tests()
    @testset "FACTS base-impedance guard is a hard error" begin (function ()
      net = _scf_test_net()
      d = mktempdir()
      # a compensated operating point must never be written as equipment data
      net.branchVec[1].x_pu = net.branchVec[1].x_base_pu * 0.5
      @test_throws ArgumentError exportSCF(net; file = joinpath(d, "bad.scf.json"))
      net.branchVec[1].x_pu = net.branchVec[1].x_base_pu
      @test isfile(exportSCF(net; file = joinpath(d, "good.scf.json")))
    end)() end

    @testset "Web UI save case as route (issue #378)" begin (function ()
      root = mktempdir()
      cases = joinpath(root, "cases")
      mkpath(cases)
      cp(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")), joinpath(cases, "sp_case14.scf.json"))
      cp(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.measurements.csv")), joinpath(cases, "sp_case14.measurements.csv"))
      rt = (; case_directory = cases, config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log = SparlectraApp.webui_operation_log_path(root), startup_config_error = nothing, runner = SparlectraApp.start_powerflow_run)

      # test 1 (issue text): save with a bound measurement set - three files
      # appear, the measurement file carries the rewritten binding, the new
      # case lists exactly this set
      resp = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/case/save-as", Dict{String,Any}("casefile" => "sp_case14.scf.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "save_as_name" => "sp_case14_WLS_B"); output_root = root, runtime = rt)
      @test resp.status == 303
      loc = only(h.second for h in resp.headers if h.first == "Location")
      @test occursin("casefile=sp_case14_WLS_B.scf.json", loc)
      @test isfile(joinpath(cases, "sp_case14_WLS_B.scf.json"))
      @test isfile(joinpath(cases, "sp_case14_WLS_B.measurements.csv"))
      @test any(l -> l == "# case: sp_case14_WLS_B.scf.json", eachline(joinpath(cases, "sp_case14_WLS_B.measurements.csv")))
      copy_case = Sparlectra.read_scf_json(joinpath(cases, "sp_case14_WLS_B.scf.json"))
      @test copy_case.sparlectra.meta["case_name"] == "sp_case14_WLS_B"
      @test copy_case.sparlectra.meta["source_reference"] == "copied from sp_case14.scf.json"
      se_state = SparlectraApp._webui_se_form_state(Dict{String,Any}("case" => "sp_case14_WLS_B.scf.json"); output_root = root, case_directory = cases)
      bound_to_copy = [name for (name, c) in se_state.meas_case if c == "sp_case14_WLS_B.scf.json"]
      @test bound_to_copy == ["sp_case14_WLS_B.measurements.csv"]

      # test 2: a MATPOWER case is saved as SCF, the .m source is untouched,
      # and the copy runs to the same result
      mfile = joinpath(cases, "case_api.m")
      write(mfile, """
      function mpc = case_api
      mpc.version = '2';
      mpc.baseMVA = 100;
      mpc.bus = [
      1 3 0 0 0 0 1 1.0 0 110 1 1.1 0.9;
      2 1 40 15 0 0 1 1.0 0 110 1 1.1 0.9;
      ];
      mpc.gen = [
      1 100 0 300 -300 1.02 100 1 300 0;
      ];
      mpc.branch = [
      1 2 0.01 0.05 0.0 999 999 999 0 0 1 -360 360;
      ];
      """)
      before_m = read(mfile, String)
      resp_m = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/case/save-as", Dict{String,Any}("casefile" => "case_api.m", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "save_as_name" => "case_api_saved"); output_root = root, runtime = rt)
      @test resp_m.status == 303
      @test read(mfile, String) == before_m
      @test isfile(joinpath(cases, "case_api_saved.scf.json"))
      original_run = SparlectraApp.start_powerflow_run(Dict("casefile" => "case_api.m", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(root, "runs_original")); case_directory = cases)
      saved_run = SparlectraApp.start_powerflow_run(Dict("casefile" => "case_api_saved.scf.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(root, "runs_saved")); case_directory = cases)
      @test original_run["status"] == "succeeded" && saved_run["status"] == "succeeded"
      @test isapprox(original_run["final_mismatch"], saved_run["final_mismatch"]; atol = 1e-12)

      # "start from the solved state" writes a real solved state, silently
      resolved_state_resp = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/case/save-as", Dict{String,Any}("casefile" => "case_api.m", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "save_as_name" => "case_api_start_state", "save_as_start_state" => "true"); output_root = root, runtime = rt)
      @test resolved_state_resp.status == 303
      start_state_case = Sparlectra.read_scf_json(joinpath(cases, "case_api_start_state.scf.json"))
      @test start_state_case.sparlectra.start_state.source == "solved_power_flow"
      # the saved copy states the frequency of the run configuration
      @test start_state_case.sparlectra.meta["f_nom"] == 50.0
      resp_f60 = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/case/save-as", Dict{String,Any}("casefile" => "case_api.m", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "save_as_name" => "case_api_60hz", "model_nominal_frequency_hz" => "60.0"); output_root = root, runtime = rt)
      @test resp_f60.status == 303
      @test Sparlectra.read_scf_json(joinpath(cases, "case_api_60hz.scf.json")).sparlectra.meta["f_nom"] == 60.0

      # test 3: an unsaved (not yet on-disk) solver change lands in the
      # COPY's sidecar; the source case has no sidecar at all and stays that way
      resp_solver = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/case/save-as", Dict{String,Any}("casefile" => "sp_case14.scf.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "save_as_name" => "sp_case14_solver_change", "power_flow_solver" => "apslf"); output_root = root, runtime = rt)
      @test resp_solver.status == 303
      solver_cfg = read(joinpath(cases, "sp_case14_solver_change.config.yaml"), String)
      @test occursin("case: sp_case14_solver_change.scf.json", solver_cfg)
      @test occursin("solver: apslf", solver_cfg)
      @test !isfile(joinpath(cases, "sp_case14.config.yaml"))
      # a pre-existing source sidecar survives untouched, and merges into the copy
      write(joinpath(cases, "sp_case14.config.yaml"), "config_version: 1\nscope: case\ncase: sp_case14.scf.json\npower_flow:\n  tol: 1.0e-8\n")
      before_source_cfg = read(joinpath(cases, "sp_case14.config.yaml"), String)
      resp_merge = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/case/save-as", Dict{String,Any}("casefile" => "sp_case14.scf.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "save_as_name" => "sp_case14_merge_test"); output_root = root, runtime = rt)
      @test resp_merge.status == 303
      @test read(joinpath(cases, "sp_case14.config.yaml"), String) == before_source_cfg
      merged_cfg = read(joinpath(cases, "sp_case14_merge_test.config.yaml"), String)
      @test occursin("tol: 1.0e-8", merged_cfg)
      @test occursin("case: sp_case14_merge_test.scf.json", merged_cfg)

      # test 4: refuse an existing name without the overwrite flag, accept with it
      dup_form = Dict{String,Any}("casefile" => "sp_case14.scf.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "save_as_name" => "sp_case14_WLS_B")
      resp_dup = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/case/save-as", dup_form; output_root = root, runtime = rt)
      loc_dup = only(h.second for h in resp_dup.headers if h.first == "Location")
      @test occursin("already%20exists", loc_dup) || occursin("already+exists", loc_dup)
      resp_dup_ow = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/case/save-as", merge(dup_form, Dict{String,Any}("save_as_overwrite" => "true")); output_root = root, runtime = rt)
      @test resp_dup_ow.status == 303

      # an invalid name (path separator) is refused before anything is written
      resp_invalid = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/case/save-as", Dict{String,Any}("casefile" => "sp_case14.scf.json", "save_as_name" => "sub/dir"); output_root = root, runtime = rt)
      loc_invalid = only(h.second for h in resp_invalid.headers if h.first == "Location")
      @test occursin("Invalid%20name", loc_invalid) || occursin("Invalid+name", loc_invalid)
      @test !isfile(joinpath(cases, "sub"))
    end)() end

  return true
end

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
# purpose: legacy profile: the testsets of test/test_scf.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_scf.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_scf_tests()
  # cut from testset: scf case format (writer)
    @testset "legacy key names in the in-file config block (0.13.0)" begin (function ()
      # a 0.10.x export carried model keys under their old prefixes in the
      # deprecated sparlectra.config block; the run refused them as "not
      # case scope" although they only have a new name (Web UI run
      # 58b55f73, case14.scf.json exported by 0.10.0)
      d = mktempdir()
      src = joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")
      root = Sparlectra.scf_json_parse(read(src, String))
      root["sparlectra"]["config"] = Dict{String,Any}("matpower_import.auto_profile" => "recommend", "transformer.tap_changer_model" => "ideal", "matpower_import.bus_shunt_model" => "admittance", "power_flow.tol" => 1.0e-7)
      legacy = joinpath(d, "sp_case14.scf.json")
      write(legacy, Sparlectra.scf_json_string(root))
      resolved = run_with_expected_warnings(() -> Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, legacy), (r"legacy configuration key name", r"is deprecated"))
      @test resolved.merged_overrides["model.auto_profile"] == "recommend"
      @test resolved.merged_overrides["model.tap_changer_model"] == "ideal"
      @test resolved.merged_overrides["power_flow.tol"] == 1.0e-7
      @test !haskey(resolved.merged_overrides, "matpower_import.auto_profile")
      @test resolved.config.powerflow.tol == 1.0e-7
      # a key that is not case scope under ANY name is still refused
      root["sparlectra"]["config"] = Dict{String,Any}("benchmark.enabled" => true)
      write(legacy, Sparlectra.scf_json_string(root))
      @test_throws Sparlectra.ConfigResolveError Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, legacy)
      # a file written before the scope split (created_by says so) keeps
      # loading: the machine keys are dropped, named, and never applied
      root["sparlectra"]["config"] = Dict{String,Any}("benchmark.enabled" => true, "output.console_summary" => false, "power_flow.tol" => 1.0e-7)
      root["sparlectra"]["meta"]["created_by"] = "Sparlectra 0.10.0"
      write(legacy, Sparlectra.scf_json_string(root))
      old = run_with_expected_warnings(() -> Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, legacy), (r"outside the case scope", r"is deprecated"))
      @test old.scf_dropped_keys == ["benchmark.enabled", "output.console_summary"]
      @test old.config.powerflow.tol == 1.0e-7
      @test old.config.benchmark.enabled == Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true).benchmark.enabled
      dropped = String[]
      @test !haskey(Sparlectra.scf_case_config(legacy; dropped), "benchmark.enabled")
      @test dropped == ["benchmark.enabled", "output.console_summary"]
      # the same keys from a current writer stay a hard error
      root["sparlectra"]["meta"]["created_by"] = string("Sparlectra ", Sparlectra.version())
      write(legacy, Sparlectra.scf_json_string(root))
      @test_throws Sparlectra.ConfigResolveError Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, legacy)
    end)() end

    @testset "FACTS base-impedance guard is a hard error" begin (function ()
      net = _scf_test_net()
      d = mktempdir()
      # a compensated operating point must never be written as equipment data
      net.branchVec[1].x_pu = net.branchVec[1].x_base_pu * 0.5
      @test_throws ArgumentError exportSCF(net; file = joinpath(d, "bad.scf.json"))
      net.branchVec[1].x_pu = net.branchVec[1].x_base_pu
      @test isfile(exportSCF(net; file = joinpath(d, "good.scf.json")))
    end)() end

    @testset "permanent v1 fixture keeps loading" begin (function ()
      # the regression guard from the format specification: this tracked file
      # must stay readable for every future reader version
      fixture = abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_casePST.scf.json"))
      @test isfile(fixture)
      net = importSCF(fixture)
      @test length(net.nodeVec) == 9
      @test length(net.branchVec) == 8
      @test length(net.measurements) == 59
      root = Sparlectra.scf_json_parse(read(fixture, String))
      @test root["sparlectra"]["format_version"] == Sparlectra.SCF_FORMAT_VERSION
      # it carries the tap cascade including the additional-voltage shifter
      taps = root["sparlectra"]["components"]["tap_changer"]
      @test any(t -> any(c -> c["index"] == 2 && haskey(c, "step"), t["controllers"]), taps)
      # and it still solves
      _, erg = runpf!(net, 60, 1e-8, 0)
      @test erg == 0
      # re-exporting the fixture reproduces it byte for byte (from a FRESH
      # read: the solve above moved the start voltages of `net`)
      d = mktempdir()
      pristine = importSCF(fixture)
      again = exportSCF(pristine; file = joinpath(d, "again.scf.json"), case_name = "warmup_casePST", source_format = "matpower", source_reference = "warmup_casePST.m", intended_calculations = String["power_flow", "state_estimation"], include_start_state = true, notes = root["sparlectra"]["meta"]["notes"])
      # This is the PERMANENT v1 format fixture: it keeps the version stamp of
      # the release that wrote it, on purpose, so a re-export cannot match it
      # byte for byte across a release. A plain equality broke on the
      # 0.10.0 -> 0.11.0 bump for exactly that reason, with nothing wrong in
      # the round trip (see test_scf_support.jl).
      verdict = scf_matches_fixture(again, fixture)
      @test verdict.stamp_is_current
      @test verdict.rest_identical
    end)() end

    @testset "real MATPOWER cases survive the round trip" begin (function ()
      # warmup_casePST alone did not exercise numerically named buses, shunts,
      # unlimited branch ratings or a neutral tap outside its own band; every
      # assertion here failed at least once before the round trip was fixed
      d = mktempdir()
      base_cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
      # case57 runs with the OTHER shunt model, so the reflection comparison
      # exercises both: the per-shunt models travel in the shunt_state block,
      # and the network-level default is stamped from the reader's own
      # configuration like the switching parameters
      vdi_cfg = Sparlectra.SparlectraConfig(Dict("model" => Dict("bus_shunt_model" => "voltage_dependent_injection")))
      # load_fixture_net: the real-MATPOWER quirks these legs need (numeric
      # bus names, unlimited ratings, a neutral tap outside its band) have
      # no shipped equivalent, so the legs run from the LOCAL cache only; a
      # fresh install skips them loudly instead of downloading
      real_legs = [(c, cfg) for (c, cfg) in (("case14.m", base_cfg), ("case57.m", vdi_cfg)) if large_case_path(c) !== nothing]
      length(real_legs) < 2 && println("      scf round trip: real-MATPOWER legs SKIPPED for ", join(setdiff(["case14.m", "case57.m"], first.(real_legs)), ", "), " (not in the large-case directory)")
      for (case, cfg) in real_legs
        cached = large_case_path(case)
        src = joinpath(d, case)
        cp(cached, src)
        net = SparlectraApp._se_import_case_net(src, cfg)
        a = exportSCF(net; file = joinpath(d, "a_$(case).json"), source_format = "matpower")
        # case57 declares a neutral tap outside its own regulation band; the
        # reader restores the band as stated and says so
        back = run_with_expected_warnings(() -> importSCF(a), (r"neutral tap position lies outside",))
        b = exportSCF(back; file = joinpath(d, "b_$(case).json"), source_format = "matpower")
        @test read(a, String) == read(b, String)
        # the generated component names embed the SOURCE bus number, so losing
        # it renamed everything and reshuffled the ids
        @test [getCompName(n.comp) for n in net.nodeVec] == [getCompName(n.comp) for n in back.nodeVec]
        @test [getCompName(br.comp) for br in net.branchVec] == [getCompName(br.comp) for br in back.branchVec]
        # bus ORDER is identity: ids are name-sorted, and "10" sorts before "2"
        @test [Sparlectra._scf_bus_name(net, i) for i in eachindex(net.nodeVec)] == [Sparlectra._scf_bus_name(back, i) for i in eachindex(back.nodeVec)]
        # tap ratios, neutral ratios and shunt admittances (the
        # capacitor-as-reactor regression) are covered field by field by the
        # reflection comparison below
        @test isapprox(Matrix(createYBUS(net = net)), Matrix(createYBUS(net = back)); rtol = 1e-12)
        # The Y-bus proves nothing about the OPERATING point: the reactive
        # band, the voltage setpoint and the regulation flag sit outside it,
        # and losing them moved the solved voltages by 0.06 pu on case14 while
        # every structural assertion above still passed. These four vectors
        # name the missing field directly and need no solver run.
        # EVERY field, by reflection: bus types, regulation, reactive bands,
        # per-bus voltage limits and the network's own Q-limit switching
        # parameters are all in here, and so is the next field somebody adds.
        # The network is compared as the service builds it, because the
        # switching parameters are stamped on that path.
        # with the start state, because carrying it is opt-in and the
        # comparison must see the same operating point on both sides
        full = exportSCF(net; file = joinpath(d, "full_$(case).json"), source_format = "matpower", include_start_state = true)
        # a tap-band advisory on some cases is the import's business, not a
        # finding of this round trip: captured, anything else fails
        via_service = run_with_expected_warnings(() -> Sparlectra._import_sparlectra_net(full, nothing, cfg), (r"neutral tap position",))
        diffs = scf_roundtrip_field_diffs(net, via_service)
        @test (case, diffs) == (case, String[])
        # A setpoint is written only where it ACTS. The check belongs on the
        # FILE, not on the network: the prosumer constructor gives a load its
        # own 1.0 either way, and that value has no effect on the calculation
        # - it is writing it that made the reader promote the machine.
        extra_written = Sparlectra.scf_json_parse(read(a, String))["sparlectra"]["extra"]
        @test !any(v isa AbstractDict && haskey(v, "vm_pu") && get(v, "regulated", false) !== true for (_, v) in extra_written)
        @test !any(v isa AbstractDict && haskey(v, "max_q_mvar") && !isa(v["max_q_mvar"], Real) for (_, v) in extra_written)
        # an unlimited rating (MATPOWER rateA = 0) has no PGM representation,
        # so it travels in the namespaced block and comes back as Inf
        unlimited = [i for i in eachindex(net.branchVec) if net.branchVec[i].sn_MVA !== nothing && !isfinite(net.branchVec[i].sn_MVA)]
        if !isempty(unlimited)
          root = Sparlectra.scf_json_parse(read(a, String))
          rows = vcat(get(root["data"], "line", []), get(root["data"], "generic_branch", []))
          @test !any(haskey(r, "i_n") || haskey(r, "sn") for r in rows if r["id"] in [0])   # no null slipped into the dataset
          @test !occursin("\"i_n\": null", read(a, String))
          @test all(back.branchVec[i].sn_MVA == Inf for i in unlimited)
        end
      end
      # T2: the solve comparison stays as the last line of defence, but on the
      # SMALL case only - the structural vectors above are the detector now.
      # Same cache-only premise as the loop (force: the loop already staged
      # this file into the same tempdir)
      if large_case_path("case14.m") !== nothing
        case = "case14.m"
        cached = large_case_path(case)
        src = joinpath(d, case)
        cp(cached, src; force = true)
        net = SparlectraApp._se_import_case_net(src, base_cfg)
        withstart = exportSCF(net; file = joinpath(d, "s_$(case).json"), source_format = "matpower", include_start_state = true)
        solved = importSCF(withstart)
        it_a, _ = runpf!(net, 40, 1e-10, 0)
        it_b, _ = runpf!(solved, 40, 1e-10, 0)
        @test it_a == it_b
        @test maximum(abs(getNodeVm(net.nodeVec[i]) - getNodeVm(solved.nodeVec[i])) for i in eachindex(net.nodeVec)) < 1e-12
        @test maximum(abs(net.nodeVec[i]._va_deg - solved.nodeVec[i]._va_deg) for i in eachindex(net.nodeVec)) < 1e-10
      end

      # the revision decides whether this reader may read the file at all
      # (load_fixture_net: probed on an export of the shipped sp_case5, so
      # this part runs on every install regardless of the cache gate above)
      shipped_export = exportSCF(Sparlectra.importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json"))); file = joinpath(d, "a_sp_case5.json"))
      stale = Sparlectra.scf_json_parse(read(shipped_export, String))
      stale["sparlectra"]["format_version"] = "0.9"
      staleF = joinpath(d, "stale.scf.json")
      write(staleF, Sparlectra.scf_json_string(stale))
      err_rev = try
        importSCF(staleF)
        ""
      catch e
        sprint(showerror, e)
      end
      @test occursin("format revision 0.9 found", err_rev) && occursin(Sparlectra.SCF_FORMAT_VERSION, err_rev)
      @test occursin("re-export", lowercase(err_rev))

      # "1.1" is the name this same structure carried before the release
      # (1.0 was never published), so a file written in that window must
      # still read: rejecting it would cost users their exported cases over
      # a renaming they never saw. Seen on a locally kept sp_case188 copy.
      pre = deepcopy(stale)
      pre["sparlectra"]["format_version"] = "1.1"
      preF = joinpath(d, "prerelease.scf.json")
      write(preF, Sparlectra.scf_json_string(pre))
      net_pre = importSCF(preF)
      @test length(net_pre.nodeVec) > 0

      # a setpoint without the regulation flag is a contradiction, not
      # something the reader resolves by promoting the machine
      contra = Sparlectra.scf_json_parse(read(shipped_export, String))
      target = first(k for (k, v) in contra["sparlectra"]["extra"] if v isa AbstractDict && haskey(v, "prosumption_type") && v["prosumption_type"] == "Load")
      contra["sparlectra"]["extra"][target]["vm_pu"] = 1.0
      contraF = joinpath(d, "contra.scf.json")
      write(contraF, Sparlectra.scf_json_string(contra))
      err_contra = try
        importSCF(contraF)
        ""
      catch e
        sprint(showerror, e)
      end
      @test occursin("voltage setpoint", err_contra) && occursin("regulated", err_contra)

      # a non-finite number never reaches the file as null: it is either
      # encoded or refused by name
      @test Sparlectra.scf_infinity_sentinel(Inf) == "inf"
      @test Sparlectra.scf_infinity_sentinel(-Inf) == "-inf"
      @test Sparlectra.scf_infinity_sentinel(2.5) == 2.5
      @test Sparlectra.scf_number_or_sentinel("inf", "x") == Inf
      @test_throws ArgumentError Sparlectra.scf_infinity_sentinel(NaN)
      @test_throws ArgumentError Sparlectra.scf_json_string(Dict{String,Any}("data" => Dict{String,Any}("bad" => Inf)))
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
      # a measurement CSV v1 sidecar (reader removed in 0.50.0) does not
      # stop save as or the SCF export: the case is written without it, the
      # set is not copied, and the message names it and the converter
      cp(mfile, joinpath(cases, "case_v1.m"))
      write(joinpath(cases, "case_v1.measurements.csv"), "﻿# sparlectra-measurements v1\n# case: case_v1.m\ntype,bus,from_bus,to_bus,branch_nr,link_nr,direction,value,sigma,active,id\n")
      for (route, form) in (("/powerflow/case/save-as", Dict{String,Any}("casefile" => "case_v1.m", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "save_as_name" => "case_v1_saved")), ("/powerflow/export-scf", Dict{String,Any}("casefile" => "case_v1.m", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH)))
        r_v1 = SparlectraApp.route_sparlectra_webui("POST", route, form; output_root = root, runtime = rt)
        msg_v1 = SparlectraApp._webui_urldecode(only(h.second for h in r_v1.headers if h.first == "Location"))
        @test r_v1.status == 303 && occursin("case_v1.measurements.csv", msg_v1) && occursin("convertMeasurementsV1", msg_v1) && !occursin("failed", msg_v1)
      end
      @test isfile(joinpath(cases, "case_v1_saved.scf.json")) && isfile(joinpath(cases, "case_v1.scf.json"))
      @test !isfile(joinpath(cases, "case_v1_saved.measurements.csv"))
      @test isempty(Sparlectra.importSCF(joinpath(cases, "case_v1_saved.scf.json")).measurements)
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

      # the action lives on the Case page only (one page per function):
      # no shortcut link from the SE
      # section, no separate route elsewhere
      case_page_html = String(copy(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case"; output_root = root, runtime = rt).body))
      @test occursin("/powerflow/case/save-as", case_page_html)
      se_page_html = String(copy(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=sp_case14.scf.json"; output_root = root, runtime = rt).body))
      @test !occursin("Save case as", se_page_html)
    end)() end

  return true
end

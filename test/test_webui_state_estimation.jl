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

# file: test/test_webui_state_estimation.jl
# purpose: the state estimation surface of the Web UI on warmup_casePST:
#          measurement upload and set selection, the generator, the sigma
#          table, the value editor, service SE runs with topology, tap and
#          bad-data checks, and the chain to power flow and N-1. One case
#          directory runs through all sections in order; names a later
#          section reads are declared local up front.
using Sparlectra
using Test

function run_webui_state_estimation_tests()
  @testset "state estimation page, measurement upload, and chain (SE phase 5)" begin (function ()
    root = mktempdir()
    cases = joinpath(root, "cases")
    mkpath(cases)
    case_path = joinpath(cases, "warmup_casePST.m")
    cp(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m"), case_path)
    rt = (; case_directory=cases, config_file=Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log=SparlectraApp.webui_operation_log_path(root), startup_config_error=nothing, runner=SparlectraApp.start_powerflow_run)
    bytes = s -> Vector{UInt8}(codeunits(s))
    # names one section sets and a later one reads
    local bad, gen, id1, loc, mfile, mname, n, other
    local p, page, page2, page3, pageInfo, r, regenerate
    local rejected, row, runs, saved, taps, v1, written_set

    @testset "measurement upload and set selection" begin (function ()
      # before the first set there is no run form, but the sigma table
      # is there, in a save form of its own (the generator reads it)
      noset_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test !occursin("name=\"se_mode\"", noset_page)
      @test occursin("name=\"se_sigma_lvl_n\"", noset_page) && occursin("name=\"se_sigma_lvl_1_vn_kv\"", noset_page)

      # upload classification: v1 CSV -> measurement_set; CSV without the
      # version comment -> unknown, retained but never offered
      v1 = "# sparlectra-measurements v1\n# case: warmup_casePST.m\ntype,bus,from_bus,to_bus,link_nr,direction,value,sigma,active,id\n"
      up = SparlectraApp.handle_powerflow_case_import(
          Dict{String,Any}("casefiles" => [SparlectraApp.WebUICaseUpload("meas_v1.csv", bytes(v1)), SparlectraApp.WebUICaseUpload("notes.csv", bytes("a,b\n1,2\n"))]);
          output_root=root, case_directory=cases, operation_log=root)
      @test up.status == 303
      loc = first(p for (k, p) in up.headers if k == "Location")
      @test occursin("measurement_set", SparlectraApp._webui_urldecode(loc))
      @test isfile(joinpath(cases, "meas_v1.csv"))
      @test isfile(joinpath(cases, "notes.csv"))   # retained
      offered = SparlectraApp._webui_measurement_options_in_directory(cases)
      @test "meas_v1.csv" in offered
      @test !("notes.csv" in offered)
      # the run picker lists only the sets of the selected case: the
      # `# case:` line decides, a file without one counts by its name
      # (case stem plus a dot; for x.scf.json the stem x as well)
      listdir = mktempdir()
      sethead = c -> string("# sparlectra-measurement-table v1\n", isempty(c) ? "" : "# case: $(c)\n", "timestamp,measurement_mrid,terminal_mrid,node_mrid,type,value,sigma,source,valid\n")
      for (fname, c) in (("case14.scf.measurements.csv", "case14.scf.json"), ("sp_case14.measurements.csv", "sp_case14.scf.json"), ("demo_clean.csv", "sp_case14.scf.json"), ("sp_case14.old.csv", ""), ("other.csv", ""), ("sp_case140.measurements.csv", ""))
          write(joinpath(listdir, fname), sethead(c))
      end
      # the header reader removes a leading BOM, skips blank lines and
      # takes CRLF: without that a BOM set lost its binding, so the
      # user's own set vanished and a foreign one matched by its name
      write(joinpath(listdir, "bom_set.csv"), "\ufeff" * sethead("sp_case14.scf.json"))
      write(joinpath(listdir, "sp_case14.foreign.csv"), "\ufeff" * sethead("case14.scf.json"))
      write(joinpath(listdir, "crlf_set.csv"), replace(sethead("sp_case14.scf.json"), "\n" => "\r\n"))
      write(joinpath(listdir, "blank_set.csv"), replace(sethead("sp_case14.scf.json"), "v1\n" => "v1\n\n"; count = 1))
      list_for = (c; kw...) -> SparlectraApp._webui_measurement_sets_for_case(listdir, c; kw...)
      @test list_for("sp_case14.scf.json") == [(name = n, case = c, model = "") for (n, c) in (("blank_set.csv", "sp_case14.scf.json"), ("bom_set.csv", "sp_case14.scf.json"), ("crlf_set.csv", "sp_case14.scf.json"), ("demo_clean.csv", "sp_case14.scf.json"), ("sp_case14.measurements.csv", "sp_case14.scf.json"), ("sp_case14.old.csv", ""))]
      @test list_for("case14.scf.json") == [(name = "case14.scf.measurements.csv", case = "case14.scf.json", model = ""), (name = "sp_case14.foreign.csv", case = "case14.scf.json", model = "")]
      @test isempty(list_for("case57.m"))
      # an SCF case exported from case14.scf.json accepts that source's sets
      @test [s.name for s in list_for("case14_copy.scf.json"; also_bound_to = ["case14.scf.json"])] == ["case14.scf.measurements.csv", "sp_case14.foreign.csv"]
      # the SE section of the page takes its list from there
      write(joinpath(cases, "foreign.csv"), sethead("case14.m"))
      @test SparlectraApp._webui_se_form_state(Dict{String,Any}("case" => "warmup_casePST.m"); output_root=root, case_directory=cases).measurements == ["meas_v1.csv"]
      rm(joinpath(cases, "foreign.csv"))
      # a measurement table (version line plus header) is recognized on
      # upload, offered, shown with its counts per type and handed out;
      # the CSV v1 editors stay off for it
      tbl = "# sparlectra-measurement-table v1\n# case: warmup_casePST.m\ntimestamp,measurement_mrid,terminal_mrid,node_mrid,type,value,sigma,source,valid\n,,,B1,U,110.0,0.2,test,true\n,,line_1:from,,P,10.0,0.1,test,true\n,,line_1:from,,Q,2.0,0.1,test,true\n"
      up_tbl = SparlectraApp.handle_powerflow_case_import(Dict{String,Any}("casefiles" => [SparlectraApp.WebUICaseUpload("meas_table.csv", bytes(tbl))]); output_root=root, case_directory=cases, operation_log=root)
      @test occursin("measurement_table", SparlectraApp._webui_urldecode(first(p for (k, p) in up_tbl.headers if k == "Location")))
      @test "meas_table.csv" in SparlectraApp._webui_measurement_options_in_directory(cases)
      tbl_state = SparlectraApp._webui_se_form_state(Dict{String,Any}("case" => "warmup_casePST.m", "measurement_file" => "meas_table.csv"); output_root=root, case_directory=cases)
      @test tbl_state.selected_measurement == "meas_table.csv"
      @test tbl_state.set_counts == [("P", 1), ("Q", 1), ("U", 1)]
      # a table opens in the text editor and the value editor as well
      @test length(tbl_state.set_rows) == 3 && startswith(tbl_state.set_text, "# sparlectra-measurement-table v1")
      @test SparlectraApp.handle_se_measurement_download(Dict{String,Any}("file" => "meas_table.csv"); output_root=root, case_directory=cases).status == 200
      # a table of another version is still offered (the run names the
      # version), but shows no counts
      write(joinpath(cases, "meas_table_v2.csv"), replace(tbl, "table v1" => "table v2"))
      @test "meas_table_v2.csv" in SparlectraApp._webui_measurement_options_in_directory(cases)
      @test isempty(SparlectraApp._webui_measurement_set_counts(joinpath(cases, "meas_table_v2.csv")))
      rm(joinpath(cases, "meas_table_v2.csv"))
      rm(joinpath(cases, "meas_table.csv"))
      # a table without its version line is named as such on upload,
      # not "unknown"
      up_unv = SparlectraApp.handle_powerflow_case_import(Dict{String,Any}("casefiles" => [SparlectraApp.WebUICaseUpload("meas_table_unv.csv", bytes(replace(tbl, "# sparlectra-measurement-table v1\n" => "")))]); output_root=root, case_directory=cases, operation_log=root)
      @test occursin("meas_table_unv.csv (measurement_table_without_version_line)", SparlectraApp._webui_urldecode(first(p for (k, p) in up_unv.headers if k == "Location")))
      rm(joinpath(cases, "meas_table_unv.csv"))
      # the structured value editor reads a table through its reader: a
      # quoted cell with commas stays one cell, the row one editor row
      qtbl = joinpath(cases, "meas_table_q.csv")
      write(qtbl, "# sparlectra-measurement-table v1\ntimestamp,measurement_mrid,terminal_mrid,node_mrid,type,value,sigma,source,valid\n,,,B1,U,110.0,0.2,\"a,b,c\",true\n")
      qrows = SparlectraApp._webui_measurement_set_rows(qtbl)
      @test length(qrows) == 1 && only(qrows).value == "110.0" && only(qrows).location == "node B1"
      rm(qtbl)
      # neither CSV appears in the case selector
      @test !any(endswith(name, ".csv") for name in SparlectraApp._webui_casefile_options_in_directory(cases))
      # an uploaded SCF time series zip is named a series, not read as a
      # CGMES delivery (it gave "unreadable CGMES delivery" before), and
      # it is offered in the case selector
      series_zip = joinpath(dirname(@__DIR__), "data", "qsts", "sp_case14_day.qsts.zip")
      upz = SparlectraApp.handle_powerflow_case_import(Dict{String,Any}("casefiles" => [SparlectraApp.WebUICaseUpload("day.qsts.zip", read(series_zip))]); output_root=root, case_directory=cases, operation_log=root)
      locz = SparlectraApp._webui_urldecode(first(p for (k, p) in upz.headers if k == "Location"))
      @test occursin("scf_time_series", locz)
      @test !occursin("CGMES", locz)
      @test "day.qsts.zip" in SparlectraApp._webui_casefile_options_in_directory(cases)
      # traversal names still rejected by the shared upload checks
      bad = SparlectraApp.handle_powerflow_case_import(Dict{String,Any}("casefiles" => [SparlectraApp.WebUICaseUpload("../evil.csv", bytes(v1))]); output_root=root, case_directory=cases, operation_log=root)
      @test bad.status == 303   # redirect with the rejection message
      @test !isfile(joinpath(dirname(cases), "evil.csv"))
    end)() end

    @testset "measurement generator and set storage" begin (function ()
      # SE page opens headless, demo generator writes an offered set,
      # and the run form carries no onsubmit button-disabling
      # /stateestimation is a real redirect (a rendering alias would be
      # two ways to one surface, the divergence the page split removed);
      # the SE section lives on the Runs page under its anchor
      page = SparlectraApp.route_sparlectra_webui("GET", "/stateestimation", Dict{String,String}(); output_root=root, runtime=rt)
      @test page.status == 303
      @test endswith(Dict(page.headers)["Location"], "#state-estimation")
      main_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test !occursin("href=\"/stateestimation\"", main_page)
      @test occursin("id=\"state-estimation\"", main_page)
      # two writes within one second get distinct names (a counter),
      # and the stamped form is what the generator and the noise action use
      stem_dir = mktempdir()
      @test occursin(r"^x\.measurements\.\d{8}-\d{6}\.csv$", SparlectraApp._webui_stamped_measurement_name(stem_dir, "x.measurements"))
      touch(joinpath(stem_dir, "x.measurements.20260927-120000.csv"))
      @test SparlectraApp._webui_stamped_measurement_name(stem_dir, "x.measurements"; stamp = "20260927-120000") == "x.measurements.20260927-120000-2.csv"
      touch(joinpath(stem_dir, "x.measurements.20260927-120000-2.csv"))
      @test SparlectraApp._webui_stamped_measurement_name(stem_dir, "x.measurements"; stamp = "20260927-120000") == "x.measurements.20260927-120000-3.csv"
      # the written set's name comes back in the redirect: every write is a
      # new, time-stamped file, so the tests take the name
      # from there instead of assuming it
      written_set(response) = SparlectraApp._webui_urldecode(String(match(r"measurement_file=([^&#]+)", Dict(response.headers)["Location"]).captures[1]))
      # a regeneration that the steps below then read: the written file
      # replaces the previous one of the test (which is removed)
      regenerate(form) = begin
          r = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", form; output_root=root, runtime=rt)
          isfile(mfile) && rm(mfile)
          n = written_set(r)
          (joinpath(cases, n), n)
      end
      gen = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root=root, runtime=rt)
      @test gen.status == 303
      mname = written_set(gen)
      @test occursin(r"^warmup_casePST\.measurements\.\d{8}-\d{6}(-\d+)?\.csv$", mname)
      mfile = joinpath(cases, mname)
      @test isfile(mfile)
      @test isMeasurementTable(mfile)
      # a rejected generation writes no file and says why in the
      # operation log (it used to leave no trace there at all)
      rejected = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "sigma_u_pct" => "-1"); output_root=root, runtime=rt)
      @test rejected.status == 303
      @test !occursin("measurement_file=", Dict(rejected.headers)["Location"])
      @test occursin(r"\"event\":\"se_measurements_rejected\".*sigma U must be a positive percentage", read(rt.operation_log, String))
      # the written set is armed for the next run: the redirect names it
      # and the page selects it (a noisy copy next to its
      # noise-free original was otherwise never the armed one)
      noised = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/add-noise", Dict{String,Any}("casefile" => "warmup_casePST.m", "measurement_file" => mname); output_root=root, runtime=rt)
      @test noised.status == 303
      noisy_location = Dict(noised.headers)["Location"]
      nname = written_set(noised)
      @test occursin(r"^warmup_casePST\.noisy\.\d{8}-\d{6}(-\d+)?\.csv$", nname)
      @test isfile(joinpath(cases, nname))
      # a browser never sends the #fragment of a redirect target
      noisy_page = String(SparlectraApp.route_sparlectra_webui("GET", first(split(noisy_location, '#')), Dict{String,String}(); output_root=root, runtime=rt).body)
      @test occursin("<option value=\"$(nname)\" selected>", noisy_page)
      @test occursin("<option value=\"$(mname)\">", noisy_page)
      @test occursin("/stateestimation/measurements/delete", noisy_page)
      # deleting a set: measurement files only, plain names only
      del_bad = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/measurements/delete", Dict{String,Any}("casefile" => "warmup_casePST.m", "measurement_file" => "../evil.csv"); output_root=root, runtime=rt)
      @test occursin("Invalid", SparlectraApp._webui_urldecode(Dict(del_bad.headers)["Location"]))
      del_case = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/measurements/delete", Dict{String,Any}("casefile" => "warmup_casePST.m", "measurement_file" => "warmup_casePST.m"); output_root=root, runtime=rt)
      @test isfile(joinpath(cases, "warmup_casePST.m"))
      @test occursin("Invalid", SparlectraApp._webui_urldecode(Dict(del_case.headers)["Location"]))
      del_ok = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/measurements/delete", Dict{String,Any}("casefile" => "warmup_casePST.m", "measurement_file" => nname); output_root=root, runtime=rt)
      @test del_ok.status == 303
      @test occursin("Deleted $(nname)", SparlectraApp._webui_urldecode(Dict(del_ok.headers)["Location"]))
      @test !isfile(joinpath(cases, nname))
      # "Save settings for this case" from the SE section keeps the chosen
      # set: the redirect names it and lands on the section, and the
      # saved choice outranks a newer bound set on a plain page load
      # (the choice was written but never read back)
      saved = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "warmup_casePST.m", "settings_target" => "this_case", "return_to" => "runs", "measurement_file" => mname, "se_flatstart" => "true"); output_root=root, runtime=rt)
      @test saved.status == 303
      saved_location = Dict(saved.headers)["Location"]
      @test occursin("measurement_file=$(SparlectraApp._webui_urlencode(mname))", saved_location)
      @test endswith(saved_location, "#state-estimation")
      newer = joinpath(cases, "newer_bound.measurements.csv")
      cp(mfile, newer)
      plain_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test occursin("<option value=\"$(mname)\" selected>", plain_page)
      @test occursin("<option value=\"newer_bound.measurements.csv\">", plain_page)
      rm(newer)
      # a saved choice whose file is gone falls back to the next rule and
      # says so; then the choice is set back to the set the steps below use
      SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "warmup_casePST.m", "settings_target" => "this_case", "return_to" => "runs", "measurement_file" => "gone.measurements.csv"); output_root=root, runtime=rt)
      gone_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test occursin("gone.measurements.csv is no longer in the case directory", gone_page)
      @test occursin("<option value=\"$(mname)\" selected>", gone_page)
      SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "warmup_casePST.m", "settings_target" => "this_case", "return_to" => "runs", "measurement_file" => mname); output_root=root, runtime=rt)
      page2 = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test occursin(mname, page2)
      @test occursin("name=\"se_mode\"", page2)
      @test !occursin("onsubmit", page2)
    end)() end

    @testset "sigma by voltage level in the form" begin (function ()
      # the sigma table posts one row per nominal voltage of the case;
      # Save settings writes it into the case configuration, the page
      # shows the saved values again, a run carries it as an override,
      # an empty table saves an empty mapping and a bad cell is refused
      case_vn = sort!(unique([nd.comp.cVN for nd in Sparlectra.import_case(case_path, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH)).net.nodeVec]); rev = true)
      posted_vn(page) = [parse(Float64, m.captures[1]) for m in eachmatch(r"name=\"se_sigma_lvl_\d+_vn_kv\" value=\"([^\"]+)\"", page)]
      @test posted_vn(page2) == case_vn
      sigma_key = "state_estimation.sigma_by_voltage_level"
      sigma_form = Dict{String,Any}("casefile" => "warmup_casePST.m", "settings_target" => "this_case", "return_to" => "runs", "measurement_file" => mname, "se_sigma_lvl_n" => string(length(case_vn)))
      for (i, v) in enumerate(case_vn)
          sigma_form["se_sigma_lvl_$(i)_vn_kv"] = string(v)
          sigma_form["se_sigma_lvl_$(i)_p_mw"] = i == 1 ? "1.5" : ""
          sigma_form["se_sigma_lvl_$(i)_u_kv"] = i == 1 ? "0.4" : ""
      end
      @test SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", sigma_form; output_root=root, runtime=rt).status == 303
      saved_levels = Sparlectra._parse_sigma_by_voltage_level(Sparlectra.load_case_config(case_path)[sigma_key])
      @test length(saved_levels) == 1 && saved_levels[1].vn_kv == case_vn[1] && saved_levels[1].p_mw == 1.5 && saved_levels[1].u_kv == 0.4 && saved_levels[1].q_mvar === nothing
      @test Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH).state_estimation.sigma_by_voltage_level == Sparlectra.SESigmaLevel[]
      page3 = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test occursin("name=\"se_sigma_lvl_1_p_mw\" value=\"1.5\"", page3)
      @test posted_vn(page3) == case_vn
      # round trip: the file keys the row by the nominal voltage, and
      # the form the reloaded page posts gives the saved table again
      @test occursin("\n    230:\n", read(Sparlectra.case_config_path(case_path), String))
      page_form(page) = Dict{String,Any}(String(m.captures[1]) => String(m.captures[2]) for m in eachmatch(r"name=\"(se_sigma_lvl_[^\"]+)\" value=\"([^\"]*)\"", page))
      full_form = merge(sigma_form, Dict{String,Any}("se_sigma_lvl_1_q_mvar" => "2.5", "se_sigma_lvl_1_i_a" => "30", "se_sigma_lvl_1_rel_pct" => "0.5"))
      @test SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", full_form; output_root=root, runtime=rt).status == 303
      full_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
      stored_table = Sparlectra.load_case_config(case_path)[sigma_key]
      @test stored_table == SparlectraApp._webui_sigma_levels_from_form(full_form) == SparlectraApp._webui_sigma_levels_from_form(page_form(full_page))
      @test collect(keys(stored_table)) == ["230"]
      # a level with a decimal point (15.5 kV) is saved like any other:
      # the file keys the row 15.5: and loads back to the same level
      @test Sparlectra._sigma_level_key(case_vn[2]) == "15.5"
      frac_save = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", merge(full_form, Dict{String,Any}("se_sigma_lvl_2_p_mw" => "0.5")); output_root=root, runtime=rt)
      @test frac_save.status == 303
      @test !occursin("bad_vn_kv", SparlectraApp._webui_urldecode(Dict(frac_save.headers)["Location"]))
      @test occursin("\n    15.5:\n      p_mw: 0.5\n", read(Sparlectra.case_config_path(case_path), String))
      frac_levels = Sparlectra._parse_sigma_by_voltage_level(Sparlectra.load_case_config(case_path)[sigma_key])
      @test [r.vn_kv for r in frac_levels] == [230.0, 15.5] && frac_levels[2].p_mw == 0.5 && frac_levels[1] == Sparlectra._parse_sigma_by_voltage_level(stored_table)[1]
      @test SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", sigma_form; output_root=root, runtime=rt).status == 303
      run_request = SparlectraApp.powerflow_webui_request(merge(sigma_form, Dict{String,Any}("se_mode" => "true")); case_directory = cases)
      @test Sparlectra._parse_sigma_by_voltage_level(run_request["config_overrides"][sigma_key]) == saved_levels
      # a request without the table leaves the configuration alone
      @test !haskey(SparlectraApp.powerflow_webui_request(Dict{String,Any}("casefile" => "warmup_casePST.m", "se_mode" => "true"); case_directory = cases)["config_overrides"], sigma_key)
      bad_form = merge(sigma_form, Dict{String,Any}("se_sigma_lvl_1_p_mw" => "0"))
      bad_save = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", bad_form; output_root=root, runtime=rt)
      @test occursin("bad_value", SparlectraApp._webui_urldecode(Dict(bad_save.headers)["Location"]))
      @test Sparlectra._parse_sigma_by_voltage_level(Sparlectra.load_case_config(case_path)[sigma_key]) == saved_levels
      # clearing the cells of a stored table clears it (the page marks
      # the table as stored); an empty table that was empty on the page
      # posts nothing, so it changes no run
      @test occursin("name=\"se_sigma_lvl_stored\"", page3)
      empty_form = merge(sigma_form, Dict{String,Any}("se_sigma_lvl_1_p_mw" => "", "se_sigma_lvl_1_u_kv" => "", "se_sigma_lvl_stored" => "1"))
      @test SparlectraApp._webui_sigma_levels_from_form(delete!(copy(empty_form), "se_sigma_lvl_stored")) === nothing
      SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", empty_form; output_root=root, runtime=rt)
      @test isempty(Sparlectra._parse_sigma_by_voltage_level(Sparlectra.load_case_config(case_path)[sigma_key]))
      # a stored table the parser refuses (flow style): the page names
      # the error and posts no table, so a run is not silently reset to
      # empty cells and stops with the error like an API run
      ccfg = Sparlectra.case_config_path(case_path)
      ccfg_text = read(ccfg, String)
      flow_text = replace(ccfg_text, r"\n  sigma_by_voltage_level:[^\n]*(\n    [^\n]*)*" => "\n  sigma_by_voltage_level: {p_mw: 2.0}")
      @test flow_text != ccfg_text
      write(ccfg, flow_text)
      flow_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test occursin("flow_style", flow_page)
      @test !occursin("name=\"se_sigma_lvl_n\"", flow_page)
      @test occursin("name=\"se_mode\"", flow_page)
      write(ccfg, ccfg_text)
      # the cells show the case configuration only, otherwise they are
      # empty: a table in the configuration file is not shown
      mktempdir() do d
          cp(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"), joinpath(d, "sp_case14.scf.json"))
          ucfg = joinpath(d, "user.yaml")
          write(ucfg, replace(read(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, String), "\n  sigma_by_voltage_level: {}" => "\n  sigma_by_voltage_level:\n    110:\n      p_mw: 3.0"))
          p110(t) = only(filter(r -> r.vn_kv == 110.0, t.rows)).values["p_mw"]
          no_case_cfg = redirect_stdout(() -> SparlectraApp._webui_sigma_level_table("sp_case14.scf.json", d, ucfg), devnull)
          @test no_case_cfg.editable && p110(no_case_cfg) == "" && !no_case_cfg.stored
          write(joinpath(d, "sp_case14.config.yaml"), "config_version: 1\nscope: case\ncase: sp_case14.scf.json\n")
          with_case_cfg = redirect_stdout(() -> SparlectraApp._webui_sigma_level_table("sp_case14.scf.json", d, ucfg), devnull)
          @test with_case_cfg.editable && p110(with_case_cfg) == "" && !with_case_cfg.stored
          write(joinpath(d, "sp_case14.config.yaml"), "config_version: 1\nscope: case\ncase: sp_case14.scf.json\nstate_estimation:\n  sigma_by_voltage_level:\n    110:\n      p_mw: 2.0\n")
          stored_cfg = redirect_stdout(() -> SparlectraApp._webui_sigma_level_table("sp_case14.scf.json", d, ucfg), devnull)
          @test stored_cfg.editable && p110(stored_cfg) == "2.0" && stored_cfg.stored
      end
    end)() end

    @testset "generator options, sigmas, currents and tap deviation" begin (function ()
      # generator options: noise + gross error produce a valid, different set
      plain = read(mfile, String)
      gen2 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "noise" => "true", "gross_error_k" => "8"); output_root=root, runtime=rt)
      @test gen2.status == 303
      @test occursin("bad%20data", Dict(gen2.headers)["Location"])
      gen2_file = joinpath(cases, written_set(gen2))
      @test gen2_file != mfile   # a new file, the clean set stays
      @test isMeasurementTable(gen2_file)
      @test read(gen2_file, String) != plain   # noise + gross error changed values
      @test read(mfile, String) == plain
      rm(gen2_file)
      # regenerate the clean set: another new file with the same content;
      # the resolution fields at 0 (the default, no rounding) leave every
      # byte as it is
      gen3 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_resolution_u_kv" => "0", "gen_resolution_p_mw" => "0", "gen_resolution_q_mvar" => "", "gen_resolution_i_a" => "0", "gen_resolution_angle_deg" => "0"); output_root=root, runtime=rt)
      gen3_file = joinpath(cases, written_set(gen3))
      @test read(gen3_file, String) == plain   # seeded generator is reproducible
      rm(gen3_file)
      # with a resolution every telemetry value (after noise and bad data)
      # is a multiple of its step, in table units, and the header names the
      # steps; zero-injection rows keep their value
      steps = Dict("U" => 0.1, "P" => 0.5, "Q" => 0.2, "I" => 1.0, "Ia" => 0.01)
      genr = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "noise" => "true", "gross_error_k" => "8", "include_currents" => "true", "sigma_ia_deg" => "0.1", "gen_resolution_u_kv" => "0.1", "gen_resolution_p_mw" => "0.5", "gen_resolution_q_mvar" => "0.2", "gen_resolution_i_a" => "1", "gen_resolution_angle_deg" => "0.01"); output_root=root, runtime=rt)
      @test genr.status == 303
      genr_file = joinpath(cases, written_set(genr))
      rtxt = read(genr_file, String)
      @test occursin("# resolution: U 0.1 kV, P 0.5 MW, Q 0.2 Mvar, I 1.0 A, angle 0.01 deg", rtxt)
      rrows = [split(l, ",") for l in split(rtxt, "\n") if !isempty(l) && !startswith(l, "#") && !startswith(l, "timestamp")]
      tele = [r for r in rrows if !occursin("zero_injection", r[8])]
      @test count(r -> r[5] == "Ia", tele) > 0 && count(r -> r[5] == "I", tele) > 0
      @test all(r -> (k = parse(Float64, r[6]) / steps[r[5]]; abs(k - round(k)) < 1e-6), tele)
      rm(genr_file)
      genrbad = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_resolution_p_mw" => "-0.1"); output_root=root, runtime=rt)
      @test occursin("resolution%20P", Dict(genrbad.headers)["Location"])

      # per-quantity sigmas are PERCENT OF THE MEASURED VALUE with the
      # per-type floors (voltage-level independent); the currents checkbox
      # adds current-magnitude rows. Noise off, so value == truth and the
      # per-row sigma law is exactly reproducible.
      gen3 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "sigma_u_pct" => "1.0", "include_currents" => "true", "sigma_i_pct" => "2.0", "sigma_p_pct" => "2.0", "sigma_q_pct" => "1.0"); output_root=root, runtime=rt)
      @test gen3.status == 303
      # every generation is a new file: the checks below read the one
      # just written, and it becomes the set of the steps that follow
      rm(mfile)
      mname = written_set(gen3)
      mfile = joinpath(cases, mname)
      netchk = Sparlectra._import_sparlectra_net(case_path, nothing, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true))
      readMeasurementTable!(netchk; file=mfile, zero_injection=false)
      fl = measurementSigmaFloors()
      frac(t) = t == Sparlectra.VmMeas ? 0.01 : (t in (Sparlectra.QinjMeas, Sparlectra.QflowMeas) ? 0.01 : 0.02)
      @test any(m -> m.typ == Sparlectra.ImagMeas, netchk.measurements)
      for m in netchk.measurements
          @test m.sigma == max(frac(m.typ) * abs(m.value), fl[m.typ])
      end
      # a 30 kV feeder flow and a 400 kV corridor flow both get a sigma
      # proportional to their own reading, never one absolute MW value
      psig = sort([m.sigma for m in netchk.measurements if m.typ == Sparlectra.PflowMeas])
      @test last(psig) / first(psig) > 2.0
      # invalid sigma is rejected with a clear message
      genbad = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "sigma_u_pct" => "-1"); output_root=root, runtime=rt)
      @test occursin("sigma%20U", Dict(genbad.headers)["Location"])

      # PMU current-angle rows via the sigma Ia field (absolute degrees)
      genia = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "include_currents" => "true", "sigma_ia_deg" => "0.1"); output_root=root, runtime=rt)
      @test genia.status == 303
      rm(mfile)
      mname = written_set(genia)
      mfile = joinpath(cases, mname)
      netia = Sparlectra._import_sparlectra_net(case_path, nothing, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true))
      readMeasurementTable!(netia; file=mfile, zero_injection=false)
      iarows = [m for m in netia.measurements if m.typ == Sparlectra.IaMeas]
      @test !isempty(iarows) && all(m.sigma == 0.1 for m in iarows)

      # tap deviation: measurements from a shifted-tap state differ from the
      # clean set and the message names the transformer branch
      plain2 = read(mfile, String)
      gentap = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "tap_error_steps" => "3"); output_root=root, runtime=rt)
      @test gentap.status == 303
      # on the PST warmup case the injected deviation is a Delta-u PHASE
      # step, and the message names the PST branch
      @test occursin("phase%20deviation", Dict(gentap.headers)["Location"])
      @test occursin("PST%20branch", Dict(gentap.headers)["Location"])
      rm(mfile)
      mname = written_set(gentap)
      mfile = joinpath(cases, mname)
      @test read(mfile, String) != plain2
      @test isMeasurementTable(mfile)
      # the file records the tap positions the set was generated from as a
      # structured table (electrical/fixed/transferred step columns, the
      # deviated transformer carries its percentage), the SE page offers
      # the download
      taptxt = read(mfile, String)
      @test occursin("# sparlectra-taps v1", taptxt)
      @test occursin("electrical_step,fixed_step,transferred_step,generation_deviation_steps", taptxt)
      @test occursin(",3.0", taptxt)   # the deviated transformer row
      pageInfo = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test occursin("/stateestimation/measurements/download?file=$(SparlectraApp._webui_urlencode(mname))", pageInfo)
      dlr = SparlectraApp.route_sparlectra_webui("GET", "/stateestimation/measurements/download?file=$(mname)", Dict{String,String}("file" => mname); output_root=root, runtime=rt)
      @test dlr.status == 200
      @test any(k == "Content-Disposition" for (k, _) in dlr.headers)
      dlbad = SparlectraApp.route_sparlectra_webui("GET", "/stateestimation/measurements/download?file=../evil.csv", Dict{String,String}("file" => "../evil.csv"); output_root=root, runtime=rt)
      @test dlbad.status in (400, 404)
      # the commented file still parses and round-trips
      netc = Sparlectra._import_sparlectra_net(case_path, nothing, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true))
      rc = readMeasurementTable!(netc; file=mfile, zero_injection=false)
      @test !isempty(rc.measurements) && isempty(rc.problems)
      # out-of-range percentage is rejected
      gentapbad = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "tap_error_steps" => "35"); output_root=root, runtime=rt)
      @test occursin("between%20-16%20and%2016", Dict(gentapbad.headers)["Location"])
      # half steps are no longer settable (a tap changer has no half positions)
      genthalf = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "tap_error_steps" => "1.5"); output_root=root, runtime=rt)
      @test occursin("whole%20number", Dict(genthalf.headers)["Location"])
      # restore the default set for the runs below (a new file; the
      # deviated one goes)
      gen5 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root=root, runtime=rt)
      rm(mfile)
      mname = written_set(gen5)
      mfile = joinpath(cases, mname)
      @test read(mfile, String) == plain
    end)() end

    @testset "page state and the value editor" begin (function ()
      # the SE page without a query FOLLOWS the shared selected-case
      # memory; a truly
      # fresh state (no memory under a fresh output root) still shows no
      # set info out of thin air
      SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m"; output_root=root, runtime=rt)
      pageRemembered = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test occursin("warmup_casePST.m", pageRemembered)
      pageFresh = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow", Dict{String,String}(); output_root=mktempdir(), runtime=rt).body)
      @test !occursin("Measurement set info", pageFresh)
      @test !occursin("Case binding:", pageFresh)

      # sticky generator inputs: the generate redirect carries the values
      # back and the re-rendered form keeps them instead of the defaults
      gsticky = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "sigma_u_pct" => "0.7", "noise" => "true"); output_root=root, runtime=rt)
      loc = Dict(gsticky.headers)["Location"]
      @test occursin("g_sigma_u_pct=0.7", loc)
      @test occursin("g_noise=true", loc)
      pageSticky = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m&g_sigma_u_pct=0.7&g_noise=true", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test occursin("name=\"sigma_u_pct\" value=\"0.7\"", pageSticky)
      @test occursin("name=\"noise\" value=\"true\" checked", pageSticky)
      mfile, mname = regenerate(Dict{String,Any}("casefile" => "warmup_casePST.m"))

      # structured value editor: every measurement kind is editable in the
      # table; the update handler rewrites only value/sigma/active and a
      # single invalid entry rejects the whole save
      pageTab = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test occursin("/stateestimation/measurements/update-values", pageTab)
      rowsTab = SparlectraApp._webui_measurement_set_rows(mfile)
      @test !isempty(rowsTab)
      # the generated set is a measurement table: types U, P, Q at nodes
      # and terminals
      @test any(r -> r.typ == "U", rowsTab) && any(r -> r.typ == "P" && startswith(r.location, "terminal "), rowsTab) && any(r -> r.typ == "P" && startswith(r.location, "node "), rowsTab)
      vrow = first(r for r in rowsTab if r.typ == "U")
      rup = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/measurements/update-values", Dict{String,Any}("file" => mname, "case" => "warmup_casePST.m", "v_$(vrow.line)" => "1.0777", "s_$(vrow.line)" => vrow.sigma, "a_$(vrow.line)" => "true"); output_root=root, runtime=rt)
      @test occursin("updated%201", Dict(rup.headers)["Location"])
      @test occursin("1.0777", read(mfile, String))
      # the edited file still reads as a table, with the new value
      @test any(r -> r.line == vrow.line && r.value == 1.0777 && r.problem === :none, readMeasurementTable(mfile))
      rbadv = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/measurements/update-values", Dict{String,Any}("file" => mname, "case" => "warmup_casePST.m", "v_$(vrow.line)" => "1.05", "s_$(vrow.line)" => "-1", "a_$(vrow.line)" => "true"); output_root=root, runtime=rt)
      @test occursin("nothing%20was%20saved", Dict(rbadv.headers)["Location"])
      @test occursin("1.0777", read(mfile, String))   # rejected save left the file alone
      # a fresh clean set for the runs below; the edited one goes
      gen4 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root=root, runtime=rt)
      rm(mfile)
      mname = written_set(gen4)
      mfile = joinpath(cases, mname)
    end)() end

    @testset "run summary, topology panel and case binding" begin (function ()
      # the run summary shows the chi-square band verdict with expected value
      fake = Dict("metadata" => Dict{String,Any}("run_mode" => "se", "se_observability_quality" => "good", "se_iterations" => 4, "se_objective" => 12.5, "se_dof" => 14, "se_band_reason" => "ok", "se_j_within_3sigma" => true))
      stxt = SparlectraApp._webui_se_summary(fake)
      # J/dof leads: a bare J grows with the row count, so the same healthy set
      # reads as an alarm on a larger case (J = 104 at dof 95 is a ratio of 1.1)
      @test occursin("J/dof = 0.89", stxt) && occursin("dof = 14", stxt)
      # a run whose final solve needed fewer iterations than its maximum
      # pass shows both, so the count is not read as the reported solve's
      fake2 = Dict("metadata" => Dict{String,Any}("run_mode" => "se", "se_observability_quality" => "good", "se_iterations" => 20, "se_iterations_last_solve" => 12, "se_objective" => 12.5, "se_dof" => 14, "se_band_reason" => "ok", "se_suspicious_count" => 0, "se_eliminations" => 0))
      @test occursin("20 iteration(s) (max over the passes; final solve 12)", SparlectraApp._webui_se_summary(fake2))
      @test !occursin("max over the passes", stxt)
      @test findfirst("J/dof", stxt).start < findfirst("(J = ", stxt).start
      @test occursin("within 3", stxt)
      # a doubled set is named where the J is shown, or the page shows an
      # alarming number with no cause
      fake["metadata"]["se_duplicate_rows"] = 82
      @test occursin("82 measurement(s) repeat an already measured quantity", SparlectraApp._webui_se_summary(fake))
      delete!(fake["metadata"], "se_duplicate_rows")
      @test !occursin("repeat an already measured", SparlectraApp._webui_se_summary(fake))
      fake["metadata"]["se_j_within_3sigma"] = false
      fake["metadata"]["se_band_reason"] = "high"
      @test occursin("OUTSIDE", SparlectraApp._webui_se_summary(fake))
      # :low reads as what it is (sigmas overstate the errors), never as an
      # alarm: the case57 confusion where "OUTSIDE (low)" was read as J too big
      fake["metadata"]["se_band_reason"] = "low"
      slow = SparlectraApp._webui_se_summary(fake)
      @test occursin("far BELOW", slow)
      @test occursin("not an alarm", slow)
      @test !occursin("OUTSIDE", slow)
      fake["metadata"]["se_band_reason"] = "high"

      # topology panel renders findings and the explicit hypothesis button
      fakeT = Dict("run_id" => "t1", "success" => true, "metadata" => Dict{String,Any}("run_mode" => "se", "se_topology_findings" => [Dict{String,Any}("stage" => "precheck", "kind" => "open_element_with_flow", "location" => "branch 2 (A-B, open)", "evidence" => "Pflow x at 12 sigma", "severity" => "strong")], "se_topology_station_findings" => [Dict{String,Any}("location" => "B2 (+1 linked)", "evidence" => "Pinj_B2 (|rn| 9.1)", "notes" => ["precheck_agreement"])]))
      tsecT = SparlectraApp._webui_se_topology_section(fakeT)
      @test occursin("open_element_with_flow", tsecT)
      @test occursin("Test topology hypotheses", tsecT)

      # tap estimation surfaces: run-form checkbox, summary line, result
      # table. The release guards are in, so the label is plain again and
      # the tooltip names the guard behavior instead of an experimental flag
      @test occursin("name=\"se_tap_estimation\"", pageInfo)
      @test !occursin("experimental", pageInfo)
      @test occursin("Machine and unobservable taps are skipped", pageInfo)
      # the box is the config key: the posted false (hidden input before the
      # box) becomes the override that switches the SE paths off
      treq = SparlectraApp.powerflow_webui_request(Dict{String,Any}("casefile" => case_path, "se_mode" => "true", "se_tap_estimation" => "false", "se_tap_fallback" => "false"))
      @test treq["config_overrides"]["state_estimation.tap_estimation.enabled"] == false && treq["config_overrides"]["state_estimation.tap_estimation.fallback"] == false
      fake["metadata"]["se_tap_count"] = 2
      fake["metadata"]["se_tap_fixed"] = true
      fake["metadata"]["se_tap_j_before"] = 1211.0
      fake["metadata"]["se_tap_j_after"] = 0.02
      fake["metadata"]["se_tap_dof_before"] = 28
      fake["metadata"]["se_tap_dof_after"] = 29
      stap = SparlectraApp._webui_se_summary(fake)
      @test occursin("tap estimation: 2 transformer(s)", stap)
      @test occursin("before fixation", stap)
      @test !occursin("off-grid tap residual", stap)
      fake["metadata"]["se_tap_offgrid_residual"] = true
      @test occursin("NOT from bad data", SparlectraApp._webui_se_summary(fake))
      fake["metadata"]["se_tap_offgrid_residual"] = false
      fake["metadata"]["se_tap_estimates"] = [
          Dict{String,Any}("branch" => 3, "name" => "T1", "mrid" => "", "mode" => "ratio", "electrical_step" => 2.03, "fixed_step" => 2, "electrical_shift_step" => 0.0, "fixed_shift_step" => 0, "out_of_range" => false, "fixed" => true),
          Dict{String,Any}("branch" => 4, "name" => "T2", "mrid" => "", "mode" => "pst", "electrical_step" => 0.0, "fixed_step" => 0, "electrical_shift_step" => -0.98, "fixed_shift_step" => -1, "out_of_range" => false, "fixed" => true),
      ]
      tsec = SparlectraApp._webui_se_tap_section(fake)
      @test occursin("Transformer tap estimates", tsec)
      @test occursin("J before fixation", tsec)
      @test occursin("<td>2</td>", tsec)              # fixed step, not raw r
      @test occursin("Shift step (fixed)", tsec)      # pst columns present
      @test !occursin("mRID", tsec)                   # no mRIDs on a MATPOWER-style row set
      @test occursin("se_tap_estimates.csv", tsec)
      # runs from 0.60.0 on: one row per regulator; a transformer with a
      # ratio and a phase changer adds the regulator count to the summary,
      # a corrected step and a note show in the status
      fake["metadata"]["se_tap_regulators"] = 3
      @test occursin("tap estimation: 2 transformer(s), 3 regulator(s)", SparlectraApp._webui_se_summary(fake))
      fake["metadata"]["se_tap_estimates"] = [
          Dict{String,Any}("branch" => 3, "regulator" => 1, "name" => "T1", "changer" => "ratio", "mrid" => "", "kind" => "ratio", "alpha_deg" => nothing, "step_estimated_raw" => 2.03, "du_est" => 0.0127, "ratio_est" => 0.975, "angle_est_deg" => 0.0, "step_transferred" => 0, "step_estimated" => 2, "deviation_steps" => 2, "out_of_range" => false, "fixed" => true, "status" => "estimated", "reason" => "none", "note" => "none", "notice" => ""),
          Dict{String,Any}("branch" => 4, "regulator" => 2, "name" => "T2", "changer" => "phase", "mrid" => "", "kind" => "asymmetrical", "alpha_deg" => 90.0, "step_estimated_raw" => -0.98, "du_est" => -0.0098, "ratio_est" => 1.0, "angle_est_deg" => -0.56, "step_transferred" => -1, "step_estimated" => -1, "deviation_steps" => 0, "out_of_range" => false, "fixed" => true, "status" => "estimated", "reason" => "none", "note" => "relative_only", "notice" => "parallel member: estimated only relative to the other member"),
      ]
      tsec2 = SparlectraApp._webui_se_tap_section(fake)
      @test occursin("fixed (corrected)", tsec2)
      @test occursin("fixed (transferred step confirmed), relative only", tsec2)
      @test occursin("<td>+2</td>", tsec2)              # deviation in steps

      # case binding lives IN the file: the generated set records its case,
      # the page shows the binding, the case selector stars bound cases, and
      # a COPY under a foreign name still says which case it belongs to
      @test occursin("# case: warmup_casePST.m", read(mfile, String))
      @test SparlectraApp._webui_measurement_set_case(mfile) == "warmup_casePST.m"
      other = joinpath(cases, "case9.measurements.csv")
      cp(mfile, other)
      page3 = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
      # both sets are bound to the case; the NEWEST is armed (the copy, by
      # its modification time), the other one is offered plain
      @test occursin("value=\"case9.measurements.csv\" selected", page3)
      @test occursin("value=\"$(mname)\">", page3)
      # block 4 dropped the SE-own case selector (shared selection), and its
      # "*" bound-set marker went with it; the binding stays visible through
      # the Case binding line and the per-set labels asserted here
      # under warmup_casePST the copy is bound to warmup_casePST and labeled with
      # that case (the file decides, not its name); ANOTHER case does not offer it
      @test occursin(">case9.measurements.csv (case warmup_casePST.m)</option>", page3)
      # suite migration (no-download rule): the foreign case only has to be
      # importable under ANOTHER name; the shipped sp_case5 does that
      case9_path = joinpath(cases, "sp_case5.scf.json")
      cp(joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json"), case9_path)
      page9 = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=sp_case5.scf.json", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test !occursin("value=\"case9.measurements.csv\"", page9)
      rm(other)

      # the SE service refuses a set bound to a different case up front,
      # naming both cases (the d27fbc77 confusion: warmup_casePST set on case85)
      rGate = start_powerflow_run(Dict{String,Any}("casefile" => case9_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
      @test rGate["status"] != "succeeded"
      @test occursin("bound to case warmup_casePST.m", rGate["message"])
      @test occursin("sp_case5.scf.json", rGate["message"])

      # inline editor: page carries the file verbatim, the save handler
      # writes it back atomically, containment rejects foreign names and
      # non-v1 content
      pageEd = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test occursin("/stateestimation/measurements/save", pageEd)
      content0 = read(mfile, String)
      edited = replace(content0, "# case: warmup_casePST.m" => "# case: warmup_casePST.m\n# note: edited inline"; count=1)
      rsave = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/measurements/save", Dict{String,Any}("file" => mname, "case" => "warmup_casePST.m", "content" => edited); output_root=root, runtime=rt)
      @test rsave.status == 303
      @test occursin("saved", Dict(rsave.headers)["Location"])
      @test occursin("# note: edited inline", read(mfile, String))
      rbad1 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/measurements/save", Dict{String,Any}("file" => "../evil.csv", "case" => "warmup_casePST.m", "content" => edited); output_root=root, runtime=rt)
      @test occursin("invalid", Dict(rbad1.headers)["Location"])
      rbad2 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/measurements/save", Dict{String,Any}("file" => mname, "case" => "warmup_casePST.m", "content" => "a,b\n1,2\n"); output_root=root, runtime=rt)
      @test occursin("rejected", Dict(rbad2.headers)["Location"])
      @test occursin("# note: edited inline", read(mfile, String))   # rejected save left the file alone
      write(mfile, content0)   # restore the clean set for the runs below
    end)() end

    @testset "service SE runs: topology, taps and bad data" begin (function ()
      # service-level SE run: artifacts present, history row kind se
      r1 = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
      @test r1["status"] == "succeeded"
      id1 = r1["run_id"]
      for artifact in ("se_state.csv", "se_diagnostics.md", "se_view.md", "measurements.csv", "measurements_verdict.csv")
          @test isfile(joinpath(root, id1, artifact))
      end
      # the bus/branch CSV checkbox on Web UI jobs: an SE run writes no
      # such files and logs the export skipped; a power flow whose CSV
      # export really fails (a directory takes the place of
      # bus_voltages_complex.csv) still logs failed
      csv_events = Tuple{String,Dict{Symbol,Any}}[]
      csv_record = (event; fields...) -> (startswith(String(event), "detailed_result_csv_export") && push!(csv_events, (String(event), Dict{Symbol,Any}(fields))); nothing)
      se_job = SparlectraApp.start_webui_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile, "detailed_result_csv" => true); event_callback=csv_record)
      wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[se_job["run_id"]]["task"])
      @test SparlectraApp._POWERFLOW_WEBUI_JOBS[se_job["run_id"]]["status"] == "success"
      @test [(e, get(f, :reason, nothing)) for (e, f) in csv_events if e != "detailed_result_csv_export_enabled"] == [("detailed_result_csv_export_skipped", "se_run")]
      empty!(csv_events)
      blocked_runner = function (request; case_directory=nothing)
          mkpath(joinpath(request["output_root"], request["run_id"], "bus_voltages_complex.csv"))
          return start_powerflow_run(request; case_directory)
      end
      pf_job = SparlectraApp.start_webui_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "detailed_result_csv" => true); runner=blocked_runner, event_callback=csv_record)
      wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[pf_job["run_id"]]["task"])
      @test SparlectraApp._POWERFLOW_WEBUI_JOBS[pf_job["run_id"]]["status"] == "success"
      @test any(((e, f),) -> e == "detailed_result_csv_export_failed" && get(f, :message, "") == "Detailed result CSV artifacts were not generated.", csv_events)
      @test !any(((e, _),) -> e == "detailed_result_csv_export_skipped", csv_events)
      runs = SparlectraApp.list_powerflow_runs(root)
      row = only(r for r in runs if string(get(r, "run_id", "")) == id1)
      @test get(row, "run_mode", "") == "se"
      # the SE result page carries the chain action
      resPage = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(id1)", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test occursin("se_start_run_id", resPage)

      # topology validation on the service run: the precheck RESULT is
      # logged on every run (the clean case says so explicitly), and the
      # result page offers the explicit hypothesis-test button
      @test occursin("topology precheck: no findings", read(joinpath(root, id1, "run.log"), String))
      @test occursin("Test topology hypotheses", resPage)
      rhyp = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/topology-hypotheses", Dict{String,Any}("run_id" => id1); output_root=root, runtime=rt)
      @test rhyp.status == 303
      @test isfile(joinpath(root, id1, "topology_hypotheses.md"))
      @test occursin("Recommendations only: NOTHING has been switched", read(joinpath(root, id1, "topology_hypotheses.md"), String))
      @test occursin("topology hypothesis test:", read(joinpath(root, id1, "run.log"), String))
      resPage2 = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(id1)", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test occursin("Topology validation (advisory)", resPage2)
      # the hypothesis test judges with the settings of the run it tests:
      # with a case configuration file its case-scope SE keys hold, not
      # the general file (k_eliminate also recorded by the run, and
      # zero_injection, which only the configuration carries)
      case_cfg = case_path * ".config.yaml"
      write(case_cfg, "config_version: 1\nscope: case\ncase: warmup_casePST.m\nstate_estimation:\n  k_eliminate: 5.5\n  zero_injection: off\n")
      rcc = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
      @test rcc["status"] == "succeeded"
      idc = rcc["run_id"]
      @test get(get_powerflow_result(idc)["metadata"], "se_k_eliminate", nothing) == 5.5
      SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/topology-hypotheses", Dict{String,Any}("run_id" => idc); output_root=root, runtime=rt)
      hyp_md = read(joinpath(root, idc, "topology_hypotheses.md"), String)
      @test occursin("k_eliminate 5.5,", hyp_md) && occursin("zero_injection off.", hyp_md)
      rm(case_cfg)
      # tap-estimation SE run: a set generated with a tap deviation
      # disagrees with the model around one transformer; releasing the taps
      # absorbs the discrepancy and the fixation lands on a mechanical step
      gentap2 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "tap_error_steps" => "3"); output_root=root, runtime=rt)
      rm(mfile)
      mname = written_set(gentap2)
      mfile = joinpath(cases, mname)
      rtap = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
      @test rtap["status"] == "succeeded"
      idtap = rtap["run_id"]
      tapcsv = joinpath(root, idtap, "se_tap_estimates.csv")
      @test isfile(tapcsv)
      taps = run_csv(tapcsv)
      # one row per tap regulator (0.60.0): the identification first, then
      # the transferred, estimated and raw estimated step; every step cell
      # is a number (a state stands in status, reason and notice)
      @test taps.header == SparlectraApp._SE_TAP_CSV_COLUMNS
      @test taps.header[1:9] == ["branch", "regulator", "name", "changer", "mrid", "changer_mrid", "step_transferred", "step_estimated", "step_estimated_raw"]
      @test all(f -> all(c -> taps.number(f[c]) isa Real, 7:9), taps.rows)
      # warmup_casePST: two transformers, the PST (one phase regulator)
      # estimated, the machine transformer calculated
      @test length(taps.rows) == 2
      @test length(unique(f[1] for f in taps.rows)) == 2
      status_col = findfirst(==("status"), taps.header)
      @test count(f -> f[status_col] == "calculated", taps.rows) == 1
      md = rtap["metadata"]
      # only the PST is estimable, the machine transformer stays calculated
      @test md["se_tap_count"] == 1
      @test md["se_tap_fixed"] == true
      # whole-step injection (half steps are no longer settable): the
      # release absorbs the deviation and the fixation lands exactly on the
      # mechanical step, so J stays numerically zero on both sides
      @test md["se_tap_j_before"] < 1e-10
      @test md["se_tap_j_after"] < 1e-10
      # the deviated PST lands on a nonzero mechanical step; on this case
      # the Delta-u deviation is a step of the phase regulator
      header_cols = taps.header
      fixed_col = findfirst(==("step_estimated"), header_cols)
      kind_col = findfirst(==("kind"), header_cols)
      @test any(f -> f[kind_col] == "asymmetrical" && taps.number(f[fixed_col]) != 0, taps.rows)
      # the injected deviation is a DEVIATION against the model step:
      # the tap CSV says so, the result table marks the row as corrected
      dev_col = findfirst(==("deviation_steps"), header_cols)
      @test any(f -> f[dev_col] != "" && taps.number(f[dev_col]) != 0, taps.rows)
      tapPage = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(idtap)", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test occursin("Transformer tap estimates", tapPage)
      @test occursin("fixed (corrected)", tapPage)
      @test occursin("<th>Transferred step</th>", tapPage)
      @test occursin("tap estimation: 1 transformer(s)", tapPage)

      # only declared tap changers are released (0.60.0): the same case
      # without its mpc.sparlectra.tap_changers block keeps both
      # transformers with their ratio only, so the same set gets no tap
      # state and no tap row, and run.log names the documented deviation
      # on branch 8 as one the run cannot release
      plain_dir = mkpath(joinpath(root, "plain"))
      plain_case = joinpath(plain_dir, "warmup_casePST.m")
      write(plain_case, replace(read(case_path, String), r"mpc\.sparlectra\.tap_changers = \[.*?\];"s => ""))
      @test !occursin("tap_changers", read(plain_case, String))
      rplain = start_powerflow_run(Dict{String,Any}("casefile" => plain_case, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
      @test rplain["status"] == "succeeded"
      @test get(rplain["metadata"], "se_tap_count", 0) == 0 && get(rplain["metadata"], "se_tap_regulators", 0) == 0
      plaincsv = joinpath(root, rplain["run_id"], "se_tap_estimates.csv")
      @test !isfile(plaincsv) || isempty(run_csv(plaincsv).rows)
      @test occursin("deviation of 3.0 step(s) on branch 8, which has no tap changer the run can release", read(joinpath(root, rplain["run_id"], "run.log"), String))

      # "estimate taps" off (the key state_estimation.tap_estimation.enabled,
      # here as the override the form posts): every changer keeps its
      # transferred step, the CSV says disabled, and a set that documents a
      # tap deviation gets the hint in the status line and run.log
      rhint = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile, "config_overrides" => Dict{String,Any}("state_estimation.tap_estimation.enabled" => false)))
      @test rhint["status"] == "succeeded"
      @test rhint["metadata"]["se_set_tap_deviation"] == true && rhint["metadata"]["se_tap_estimation"] == false
      @test occursin("tap estimation is switched off; enable 'estimate taps'", rhint["message"])
      @test occursin("but tap estimation is switched off", read(joinpath(root, rhint["run_id"], "run.log"), String))
      offcsv = run_csv(joinpath(root, rhint["run_id"], "se_tap_estimates.csv"))
      @test length(offcsv.rows) == 2 && all(f -> f[status_col] == "not_estimated" && f[status_col + 1] == "disabled", offcsv.rows)

      # tap estimation per transformer on the Case page: one "estimate" box
      # per transformer with a tap changer, the machine transformer
      # unchecked. The save writes the unchecked PST into
      # state_estimation.tap_estimation.exclude of the case configuration,
      # the reloaded page posts the same choice, and the next SE run keeps
      # the PST out (reason excluded, one run.log line). A checked machine
      # transformer sets include_machine_transformers; the defaults leave
      # both keys out of the file again.
      browser_form = function (page, prefix)
        # what a browser posts: unchecked boxes and disabled fields are
        # left out, the last value of a name wins
        out = Dict{String,Any}()
        for m in eachmatch(r"<input [^>]*>", page)
          tag = m.match
          nm = match(r"name=\"([^\"]+)\"", tag)
          (nm === nothing || !startswith(nm.captures[1], prefix) || occursin(" disabled", tag)) && continue
          (occursin("type=\"checkbox\"", tag) && !occursin(" checked", tag)) && continue
          v = match(r"value=\"([^\"]*)\"", tag)
          out[String(nm.captures[1])] = v === nothing ? "" : String(v.captures[1])
        end
        return out
      end
      case_page() = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
      excl_key = "state_estimation.tap_estimation.exclude"
      mach_key = "state_estimation.tap_estimation.include_machine_transformers"
      tf = browser_form(case_page(), "se_tap_")
      @test tf["se_tap_trafo_n"] == "2"
      i_pst = only(i for i in 1:2 if tf["se_tap_trafo_$(i)_machine"] == "false")
      i_gsu = 3 - i_pst
      @test tf["se_tap_trafo_$(i_pst)_estimate"] == "true" && tf["se_tap_trafo_$(i_gsu)_estimate"] == "false"
      pst_name = tf["se_tap_trafo_$(i_pst)_key"]
      save_base = Dict{String,Any}("casefile" => "warmup_casePST.m", "settings_target" => "this_case", "return_to" => "case")
      @test SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", merge(save_base, tf, Dict{String,Any}("se_tap_trafo_$(i_pst)_estimate" => "false")); output_root=root, runtime=rt).status == 303
      excl_cfg = Sparlectra.load_case_config(case_path)
      @test excl_cfg[excl_key] == [pst_name] && !haskey(excl_cfg, mach_key)
      tf2 = browser_form(case_page(), "se_tap_")
      @test tf2["se_tap_trafo_$(i_pst)_estimate"] == "false"
      @test SparlectraApp._webui_tap_exclusions_from_form(tf2) == (exclude = [pst_name], include_machine = false)
      rex = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
      @test rex["status"] == "succeeded"
      pst_rows = [t for t in rex["metadata"]["se_tap_estimates"] if t["name"] == pst_name]
      @test !isempty(pst_rows) && all(t -> t["status"] == "not_estimated" && t["reason"] == "excluded", pst_rows)
      @test rex["metadata"]["se_tap_excluded"] == 1 && rex["metadata"]["se_tap_machine"] == 1 && rex["metadata"]["se_tap_count"] == 0
      @test occursin("transformer $(pst_name) (branch 8) excluded (state_estimation.tap_estimation.exclude)", read(joinpath(root, rex["run_id"], "run.log"), String))
      @test SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", merge(save_base, tf2, Dict{String,Any}("se_tap_trafo_1_estimate" => "true", "se_tap_trafo_2_estimate" => "true")); output_root=root, runtime=rt).status == 303
      inc_cfg = Sparlectra.load_case_config(case_path)
      @test inc_cfg[mach_key] == true && !haskey(inc_cfg, excl_key)
      @test SparlectraApp._webui_tap_exclusions_from_form(browser_form(case_page(), "se_tap_")) == (exclude = String[], include_machine = true)
      rinc = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
      @test rinc["status"] == "succeeded"
      @test all(t -> t["reason"] != "machine_transformer" && t["reason"] != "excluded", rinc["metadata"]["se_tap_estimates"])
      @test SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", merge(save_base, tf); output_root=root, runtime=rt).status == 303
      def_cfg = Sparlectra.load_case_config(case_path)
      @test !haskey(def_cfg, excl_key) && !haskey(def_cfg, mach_key)
      # a machine transformer the case configuration names in exclude: the
      # page shows it unchecked, a save without touching the table keeps the
      # entry, and the run keeps it out without calculating it (excluded)
      gsu_name = tf["se_tap_trafo_$(i_gsu)_key"]
      Sparlectra.write_case_config(case_path, merge(def_cfg, Dict{String,Any}(excl_key => [gsu_name])))
      tf3 = browser_form(case_page(), "se_tap_")
      @test tf3["se_tap_trafo_$(i_gsu)_estimate"] == "false" && tf3["se_tap_trafo_$(i_gsu)_excluded"] == "true"
      @test SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", merge(save_base, tf3); output_root=root, runtime=rt).status == 303
      @test Sparlectra.load_case_config(case_path)[excl_key] == [gsu_name]
      rgx = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
      @test rgx["status"] == "succeeded"
      gsu_rows = [t for t in rgx["metadata"]["se_tap_estimates"] if t["name"] == gsu_name]
      @test !isempty(gsu_rows) && all(t -> t["status"] == "not_estimated" && t["reason"] == "excluded", gsu_rows)
      @test SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", merge(save_base, tf); output_root=root, runtime=rt).status == 303
      @test !haskey(Sparlectra.load_case_config(case_path), excl_key)
      # an exclude entry that names no transformer: one warning in run.log,
      # no effect on the rows
      rnk = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile, "config_overrides" => Dict{String,Any}(excl_key => ["no_such_trafo"])))
      @test rnk["status"] == "succeeded"
      @test occursin("WARNING: state_estimation.tap_estimation.exclude names 'no_such_trafo', which is no transformer of this case", read(joinpath(root, rnk["run_id"], "run.log"), String))
      @test !any(t -> t["reason"] == "excluded", rnk["metadata"]["se_tap_estimates"])

      # bad data lands findable: gross error in the set -> the verdict
      # file flags the corrupted row of measurements.csv as eliminated,
      # joined by its line number; the row is the one whose value is
      # 10 sigma off its truth value (truth_value comment of the set)
      mfile, mname = regenerate(Dict{String,Any}("casefile" => "warmup_casePST.m", "gross_error_k" => "10"))
      rbd = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
      @test rbd["status"] == "succeeded"
      vrows = run_csv(joinpath(root, rbd["run_id"], "measurements_verdict.csv")).rows
      art_by_line = Dict(r.line => r for r in readMeasurementTable(joinpath(root, rbd["run_id"], "measurements.csv")))
      elim_rows = filter(v -> v[9] == "true", vrows)
      @test length(unique(v[4] for v in elim_rows)) == rbd["metadata"]["se_eliminations"] >= 1
      @test all(v -> v[3] in ("P", "Q"), elim_rows)
      truth_bd = Dict{String,Float64}()
      for l in readlines(joinpath(root, rbd["run_id"], "measurements.csv"))
          startswith(l, "# truth_value,") || continue
          p = split(chopprefix(l, "# truth_value,"), ",")
          truth_bd[join(p[1:(end-1)], ",")] = parse(Float64, p[end])
      end
      @test all(v -> (r = art_by_line[parse(Int, v[1])]; abs(r.value - truth_bd[v[4]]) > 5 * r.sigma), elim_rows)
      @test occursin("measurements_verdict.csv", read(joinpath(root, rbd["run_id"], "run.log"), String))
      # the topology hypothesis test starts from the set the run ended
      # with: measurements.csv keeps the eliminated row valid, the
      # verdict file names it, and the base J of the test is the J of
      # the run (with the row back in, it carried the gross error)
      rhbd = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/topology-hypotheses", Dict{String,Any}("run_id" => rbd["run_id"]); output_root=root, runtime=rt)
      @test rhbd.status == 303
      mbase = match(r"Base run: J = ([0-9.eE+-]+) \(dof (\d+),", read(joinpath(root, rbd["run_id"], "topology_hypotheses.md"), String))
      @test mbase !== nothing
      # the set is noise-free apart from the gross error: J of the run
      # is round-off (about 1e-18), with the row back in it was 84 at
      # dof 42 instead of 41; the file prints four digits
      @test isapprox(parse(Float64, mbase.captures[1]), rbd["metadata"]["se_objective"]; rtol=1e-3, atol=1e-6)
      @test parse(Int, mbase.captures[2]) == rbd["metadata"]["se_dof"]
      # bad data is one click away on the result page (the direct download)
      bdPage = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(rbd["run_id"])", Dict{String,String}(); output_root=root, runtime=rt).body)
      @test occursin("/powerflow/artifact/$(rbd["run_id"])/measurements_verdict.csv", bdPage)
    end)() end

    @testset "generator variants, zero injection and truth state" begin (function ()
      # restore the clean default set for anything below
      mfile, mname = regenerate(Dict{String,Any}("casefile" => "warmup_casePST.m"))

      # --- measurement generator v2: truth state, flow ends, passive nodes,
      # delta comments, and the bad-data threshold surface of the run form
      genpage = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
      for needle in ("name=\"gen_truth_source\"", "name=\"gen_truth_run_id\"", "name=\"gen_flow_ends\"", "name=\"gen_passive_sigma\"", "name=\"gen_passive_as_zi\"", "name=\"se_robust_mode\"", "name=\"se_k_eliminate\"", "name=\"se_k_suppress\"", "name=\"se_suppression_sigma\"", "se-threshold-warning", "gen-truth-source")
          @test occursin(needle, genpage)
      end

      # one balance-aware flow end: deterministic (two generates produce the
      # identical file), one flow group per branch, choice documented
      mfile, mname = regenerate(Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_flow_ends" => "one_balance_aware"))
      one1 = read(mfile, String)
      mfile, mname = regenerate(Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_flow_ends" => "one_balance_aware"))
      one2 = read(mfile, String)
      @test one1 == one2
      @test occursin("# seed: 42", one1)
      # generating persists the generator options in the form block of the
      # case configuration file so a case reload restores them
      @test occursin("gen_flow_ends: one_balance_aware", read(joinpath(cases, "warmup_casePST.m.config.yaml"), String))
      @test occursin("# flow_ends: one_balance_aware", one1)
      @test occursin("# flow_end,", one1)
      @test occursin("# truth_value,U@", one1)
      # critical measurements on request: the set is thinned until two rows
      # are critical, stays observable, and names the rows in its comments
      mfile, mname = regenerate(Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_critical_count" => "2"))
      crit_text = read(mfile, String)
      @test occursin("# critical_target: 2 reached: ", crit_text)
      @test occursin("# critical_rows: ", crit_text)
      crit_net = Sparlectra.import_case(joinpath(cases, "warmup_casePST.m"), Sparlectra.load_sparlectra_config()).net
      Sparlectra.readMeasurementTable!(crit_net; file=mfile, zero_injection=false)
      crit_obs = evaluate_global_observability(crit_net)
      @test crit_obs.quality != :not_observable
      @test length(crit_obs.numerical_critical_measurement_indices) >= 2
      @test occursin("gen_critical_count: 2", read(joinpath(cases, "warmup_casePST.m.config.yaml"), String))
      # a case a synchronous action still holds refuses a second action and
      # a run, with a message instead of a half-written file
      @test SparlectraApp._webui_case_claim!("warmup_casePST.m", "generating measurements")
      busy_resp = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root=root, runtime=rt)
      @test occursin("is busy", SparlectraApp._webui_urldecode(Dict(busy_resp.headers)["Location"]))
      busy_run = SparlectraApp.start_webui_powerflow_run(Dict{String,Any}("casefile" => "warmup_casePST.m", "output_root" => root, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH); case_directory=cases)
      @test busy_run["status"] == "failed" && occursin("busy", String(busy_run["message"]))
      SparlectraApp._webui_case_release!("warmup_casePST.m")
      @test SparlectraApp._webui_case_busy("warmup_casePST.m") === nothing
      # the full default set measures injections at every bus, so every
      # branch keeps its from end and no to-direction flow rows remain
      @test all(l -> !(startswith(l, "PflowMeas") && length(split(l, ",")) >= 7 && split(l, ",")[7] == "to"), split(one1, "\n"))
      rone = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
      @test rone["status"] == "succeeded"

      # passive nodes as protected zero-injection constraints: ZI rows
      # written, no duplicate plain injection rows, elimination stays away
      # the checkbox is state_estimation.zero_injection of the case:
      # unchecked, the set has no ZI rows and the case configuration
      # says off; checked, it says on again (read by every SE run)
      mfile, mname = regenerate(Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_passive_as_zi" => "false"))
      @test !occursin("zero_injection", read(mfile, String))
      @test Sparlectra._as_symbol_cfg(Sparlectra.load_case_config(case_path)["state_estimation.zero_injection"]) === :off
      mfile, mname = regenerate(Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_passive_as_zi" => "true"))
      @test Sparlectra._as_symbol_cfg(Sparlectra.load_case_config(case_path)["state_estimation.zero_injection"]) === :on
      zi1 = read(mfile, String)
      # the table marks them by the source token zero_injection on a
      # bus balance row
      zilines = [l for l in split(zi1, "\n") if occursin(",P,", l) && occursin("zero_injection", l)]
      @test !isempty(zilines)
      zibus = split(first(zilines), ",")[4]
      @test !any(l -> startswith(l, ",,,$(zibus),P,") && !occursin("zero_injection", l), split(zi1, "\n"))
      rzi = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
      @test rzi["status"] == "succeeded"
      @test rzi["metadata"]["se_eliminations"] == 0
      # the numerics of the run under on: the gain matrix with its
      # condition number in the metadata and in the run.log summary line,
      # no multipliers; the form field is the configuration key
      @test (rzi["metadata"]["se_zero_injection"], rzi["metadata"]["se_condition_matrix"], rzi["metadata"]["se_zero_injection_multipliers"]) == ("on", "gain", 0)
      @test rzi["metadata"]["se_condition_number"] > 1.0
      @test occursin("numerics: zero_injection on, condition number of the gain matrix", read(joinpath(root, rzi["run_id"], "run.log"), String))
      @test !isfile(joinpath(root, rzi["run_id"], "se_zero_injection_multipliers.csv"))
      zreq = SparlectraApp.powerflow_webui_request(Dict{String,Any}("casefile" => case_path, "se_mode" => "true", "se_zero_injection_sigma_mva" => "0.01"))
      @test zreq["config_overrides"]["state_estimation.zero_injection_sigma_mva"] == 0.01
      # the SE flat start behind phase shifters is a configuration key as well
      rreq = SparlectraApp.powerflow_webui_request(Dict{String,Any}("casefile" => case_path, "se_mode" => "true", "se_rotate_behind_pst" => "true"))
      @test rreq["config_overrides"]["state_estimation.start_mode.rotate_behind_pst"] === true
      # the run's zero-injection sigma reaches the rows the service adds
      # (the process default stays 1e-3): the set without its ZI rows, run
      # with the key at 0.01, gets them back at 0.01 in measurements.csv
      nozi = joinpath(root, "warmup_nozi.measurements.csv")
      write(nozi, join([l for l in split(zi1, "\n") if startswith(l, "#") || !occursin("zero_injection", l)], "\n"))
      rzs = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => nozi, "config_overrides" => Dict{String,Any}("state_estimation.zero_injection_sigma_mva" => 0.01)))
      @test rzs["status"] == "succeeded"
      @test rzs["metadata"]["se_zero_injection_rows"] > 0
      @test occursin(r"zero injection: \d+ row\(s\) added .*sigma 0\.01 MW/Mvar per row", read(joinpath(root, rzs["run_id"], "run.log"), String))
      zsig = [parse(Float64, split(l, ",")[7]) for l in eachline(joinpath(root, rzs["run_id"], "measurements.csv")) if !startswith(l, "#") && occursin("zero_injection", l)]
      @test length(zsig) == rzs["metadata"]["se_zero_injection_rows"] && all(==(0.01), zsig)
      # a case on constraint (equality-constrained zero injection): the
      # box shows checked, a generation with the box checked keeps
      # constraint, unchecked writes off; checked again it is on
      @test [SparlectraApp._webui_normalize_case_profile_form_value("gen_passive_as_zi", v) for v in ("constraint", "off", "on")] == [true, false, true]
      SparlectraApp._webui_merge_case_settings!(root, case_path, Dict{String,Any}("gen_passive_as_zi" => "constraint"); config_file=Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, case_directory=cases)
      @test Sparlectra._as_symbol_cfg(Sparlectra.load_case_config(case_path)["state_estimation.zero_injection"]) === :constraint
      # a run of the case on constraint: the KKT matrix, and the
      # multipliers per bus in se_zero_injection_multipliers.csv (the
      # passive bus of warmup_casePST sits in a closed-link cluster and
      # enters as its LINKAGG pair)
      rcon = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
      @test rcon["status"] == "succeeded"
      @test (rcon["metadata"]["se_zero_injection"], rcon["metadata"]["se_condition_matrix"]) == ("constraint", "kkt")
      nzi = rcon["metadata"]["se_zero_injection_multipliers"]
      @test nzi > 0
      zcsv = run_csv(joinpath(root, rcon["run_id"], "se_zero_injection_multipliers.csv"))
      @test length(zcsv.rows) == nzi
      @test all(row -> row[2] in ("zero_injection", "link_aggregate") && isfinite(zcsv.number(row[3])) && isfinite(zcsv.number(row[4])), zcsv.rows)
      mfile, mname = regenerate(Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_passive_as_zi" => "true"))
      @test Sparlectra._as_symbol_cfg(Sparlectra.load_case_config(case_path)["state_estimation.zero_injection"]) === :constraint
      mfile, mname = regenerate(Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_passive_as_zi" => "false"))
      @test Sparlectra._as_symbol_cfg(Sparlectra.load_case_config(case_path)["state_estimation.zero_injection"]) === :off
      mfile, mname = regenerate(Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_passive_as_zi" => "true"))
      @test Sparlectra._as_symbol_cfg(Sparlectra.load_case_config(case_path)["state_estimation.zero_injection"]) === :on

      # truth state from run: adopts the SE run's state (bit-exact against
      # se_state.csv), documents the source, and the delta file appears
      mfile, mname = regenerate(Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_truth_source" => "from_run", "gen_truth_run_id" => id1))
      fr1 = read(mfile, String)
      @test occursin("# truth: run $(id1) (se,", fr1)
      sestate = run_csv(joinpath(root, id1, "se_state.csv"))
      strow = only([r for r in sestate.rows if r[1] == "2"])
      st_vm = sestate.number(strow[2])
      # the table carries U in kV; read back it is the pu value of the
      # state file bit for bit
      frnet = Sparlectra._import_sparlectra_net(case_path, nothing, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true))
      readMeasurementTable!(frnet; file=mfile, zero_injection=false)
      @test only(m for m in frnet.measurements if m.typ == Sparlectra.VmMeas && m.busIdx == 2).value == st_vm
      # rejections: tap deviation locked, unknown run id named
      rrej = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_truth_source" => "from_run", "gen_truth_run_id" => id1, "tap_error_steps" => "2"); output_root=root, runtime=rt)
      @test occursin("requires%20truth%20state", Dict(rrej.headers)["Location"])
      rrej2 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_truth_source" => "from_run", "gen_truth_run_id" => "nope"); output_root=root, runtime=rt)
      @test occursin("not%20found%20in%20the%20run%20history", Dict(rrej2.headers)["Location"])
    end)() end

    @testset "thresholds, request options, eliminations and reset" begin (function ()
      # restore the clean default set and exercise the threshold surface:
      # the staged service path equals the legacy Bool bitwise, the delta
      # artifact exists for generated sets, invalid modes reject
      mfile, mname = regenerate(Dict{String,Any}("casefile" => "warmup_casePST.m"))
      rst = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile, "se_robust_mode" => "staged", "se_robust_k1" => 3.0, "se_robust_k2" => 6.0))
      @test rst["status"] == "succeeded"
      @test rst["metadata"]["se_robust_mode"] == "staged"
      # knees out of order: the staged run refuses them naming both form
      # fields; without staged down-weighting the pair is not read, the
      # run goes through on the configured pair and run.log says so
      knee_request = Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile, "se_robust_k1" => 10.0, "se_robust_k2" => 6.0)
      knee_staged = start_powerflow_run(merge(knee_request, Dict{String,Any}("se_robust_mode" => "staged")))
      @test knee_staged["status"] == "failed"
      @test occursin("down-weight from |r|/sigma", knee_staged["message"]) && occursin("k1=10.0, k2=6.0", knee_staged["message"])
      knee_off = start_powerflow_run(merge(knee_request, Dict{String,Any}("se_robust_mode" => "off")))
      @test knee_off["status"] == "succeeded"
      @test occursin("are out of order; they are not used with down-weighting off", read(joinpath(knee_off["output_dir"], "run.log"), String))
      # the request builder records the SE options as sidecar-persistable
      # settings, so the browser flow's "save settings" keeps them.
      # se_robust_mode carries a real config_key since issue #377 (case
      # scope, like power_flow.solver), so it now travels as the dotted
      # config override "state_estimation.robust_mode", not as the bare
      # form field; se_robust_k2 has no config key and is unaffected.
      reqrec = SparlectraApp._webui_request_settings_for_profile(SparlectraApp.powerflow_webui_request(Dict{String,Any}("se_mode" => "true", "casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "measurement_file" => mfile, "se_robust_mode" => "staged", "se_robust_k2" => "6.0"); default_output_root=root))
      @test reqrec["state_estimation.robust_mode"] == "staged"
      @test reqrec["se_robust_k2"] == 6.0
      rleg = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile, "se_robust" => true))
      @test rleg["metadata"]["se_robust_mode"] == "staged"
      @test rst["metadata"]["se_objective"] == rleg["metadata"]["se_objective"]
      @test isfile(joinpath(root, rst["run_id"], "se_deltas.csv"))
      @test run_csv(joinpath(root, rst["run_id"], "se_deltas.csv")).header[1:9] == ["kind", "id", "type", "bus", "from_bus", "to_bus", "measured", "truth", "estimated"]
      rbadmode = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile, "se_robust_mode" => "bogus"))
      @test rbadmode["status"] != "succeeded"
      @test occursin("se_robust_mode", rbadmode["message"])

      # Delta-u PST angle estimation on the tracked demo case: the
      # generator targets the PST's additional-voltage stepper (nameplate
      # block mpc.sparlectra.tap_changers), the mass release estimates it
      # (:pst along the nameplate psi), and the fixation lands exactly on
      # the injected Delta-u step while the machine trafo stays calculated
      pstcase = joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m")
      pstout = joinpath(root, "warmup_casePST.measurements.csv")
      gpst = SparlectraApp._se_generate_measurement_set(pstcase, pstout, SparlectraApp.MeasurementGeneratorOptions(; noise=true, gross_k=0.0, tap_steps=2.0, include_i=false, sigma_u_pct=0.5, sigma_i_pct=1.0, sigma_p_pct=1.0, sigma_q_pct=1.0, sigma_ia_deg=0.0))
      @test occursin("Delta-u step", gpst.tap_note)
      @test occursin("PST branch 8", gpst.tap_note)

      # random multi-selection: the user chooses HOW MANY bad-data rows and
      # at most how many transformers, the seed decides WHERE; the same
      # seed reproduces the identical file, a different seed moves the picks.
      # load_fixture_net: the shipped sp_case60 offers several eligible
      # transformers, so a max of 2 actually draws 2 (derived)
      d14 = joinpath(dirname(@__DIR__), "data", "scf", "sp_case60.scf.json")
      go1 = joinpath(root, "gen_multi_a.csv")
      go2 = joinpath(root, "gen_multi_b.csv")
      go3 = joinpath(root, "gen_multi_c.csv")
      genmulti(out, seed) = SparlectraApp._se_generate_measurement_set(d14, out, SparlectraApp.MeasurementGeneratorOptions(; noise=false, gross_k=10.0, gross_count=3, tap_steps=2.0, tap_count=2, include_i=false, sigma_u_pct=0.5, sigma_i_pct=1.0, sigma_p_pct=1.0, sigma_q_pct=1.0, sigma_ia_deg=0.0, seed=seed))
      gm1 = genmulti(go1, 42)
      @test occursin("on 3 row(s):", gm1.gross_note)
      @test length(collect(eachmatch(r"on transformer branch", gm1.tap_note))) == 2
      gm2 = genmulti(go2, 42)
      @test gm2.gross_note == gm1.gross_note
      @test read(go1, String) == read(go2, String)
      gm3 = genmulti(go3, 7)
      @test gm3.gross_note != gm1.gross_note
      # the tap draw never spills onto a machine transformer while an
      # estimable one exists: on the 2-trafo PST case a max of 2 must hit
      # ONLY the PST and say so, or the set carries a deviation the mass
      # release skips by design (seen: J ~ 318 instead of ~ dof)
      gopst = joinpath(root, "gen_multi_pst.csv")
      gmp = SparlectraApp._se_generate_measurement_set(pstcase, gopst, SparlectraApp.MeasurementGeneratorOptions(; noise=false, gross_k=0.0, tap_steps=3.0, tap_count=2, include_i=false, sigma_u_pct=0.5, sigma_i_pct=1.0, sigma_p_pct=1.0, sigma_q_pct=1.0, sigma_ia_deg=0.0, seed=42))
      @test occursin("PST branch 8", gmp.tap_note)
      @test !occursin("branch 7", gmp.tap_note)
      @test occursin("limited to 1 eligible transformer(s), max 2 requested", gmp.tap_note)
      pstlines = readlines(gopst)
      @test any(l -> startswith(l, "# 7,") && endswith(l, ",0.0"), pstlines)
      @test any(l -> startswith(l, "# 8,") && endswith(l, ",3.0"), pstlines)
      @test_throws ArgumentError SparlectraApp.MeasurementGeneratorOptions(; noise=false, gross_k=10.0, gross_count=0, tap_steps=0.0, include_i=false, sigma_u_pct=0.5, sigma_i_pct=1.0, sigma_p_pct=1.0, sigma_q_pct=1.0, sigma_ia_deg=0.0)
      @test_throws ArgumentError SparlectraApp.MeasurementGeneratorOptions(; noise=false, gross_k=0.0, tap_steps=1.0, tap_count=0, include_i=false, sigma_u_pct=0.5, sigma_i_pct=1.0, sigma_p_pct=1.0, sigma_q_pct=1.0, sigma_ia_deg=0.0)

      # J_active: replacement suppression removes corrupted rows from the
      # STATE; the reported pair (honest J with original sigmas, J_active
      # over the trusted rows) makes that visible, and the band verdict
      # stays on the honest J (eliminations off so the rows STAY suppressed)
      gja = joinpath(root, "gen_jactive.csv")
      SparlectraApp._se_generate_measurement_set(case_path, gja, SparlectraApp.MeasurementGeneratorOptions(; noise=false, gross_k=12.0, gross_count=2, tap_steps=0.0, include_i=false, sigma_u_pct=0.5, sigma_i_pct=1.0, sigma_p_pct=1.0, sigma_q_pct=1.0, sigma_ia_deg=0.0, seed=42))
      rja = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => gja, "se_robust_mode" => "replacement", "se_max_eliminations" => 0))
      @test rja["status"] == "succeeded"
      @test rja["metadata"]["se_suppressed_rows"] >= 1
      @test rja["metadata"]["se_objective_active"] < rja["metadata"]["se_objective"]
      @test rja["metadata"]["se_dof_active"] == rja["metadata"]["se_dof"] - rja["metadata"]["se_suppressed_rows"]
      @test occursin("J_active", read(joinpath(root, rja["run_id"], "run.log"), String))

      # with eliminations enabled the HEADLINE describes the state AFTER
      # the elimination workflow: the injected gross error is eliminated
      # and J returns to ~dof (regression for the run where an eliminated
      # 10-sigma row still pushed the headline to J = 149 at dof 42)
      gjb = joinpath(root, "gen_jelim.csv")
      SparlectraApp._se_generate_measurement_set(case_path, gjb, SparlectraApp.MeasurementGeneratorOptions(; noise=true, gross_k=10.0, gross_count=1, tap_steps=0.0, include_i=false, sigma_u_pct=0.5, sigma_i_pct=1.0, sigma_p_pct=1.0, sigma_q_pct=1.0, sigma_ia_deg=0.0, seed=42))
      rje = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => gjb, "se_robust_mode" => "replacement", "se_k_suppress" => 4.0))
      @test rje["status"] == "succeeded"
      @test rje["metadata"]["se_eliminations"] == 1
      @test rje["metadata"]["se_band_reason"] == "ok"
      @test rje["metadata"]["se_objective"] < 2.0 * rje["metadata"]["se_dof"]

      # reset-settings deletes the per-case sidecar profile; a second
      # reset reports that the defaults are already active
      SparlectraApp._webui_merge_case_settings!(root, case_path, Dict{String,Any}("gen_seed" => 99); config_file=Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, case_directory=dirname(case_path))
      spath = SparlectraApp._webui_case_settings_path(root, case_path; case_directory=dirname(case_path))
      @test isfile(spath)
      rrst = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/reset-settings", Dict{String,Any}("casefile" => basename(case_path)); output_root=root, runtime=rt)
      @test !isfile(spath)
      @test occursin("deleted", string(rrst))
      rrst2 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/reset-settings", Dict{String,Any}("casefile" => basename(case_path)); output_root=root, runtime=rt)
      @test occursin("already", string(rrst2))
      @test occursin("reset-settings", SparlectraApp.render_se_form(; cases=[basename(case_path)], selected_case=basename(case_path)))

      # noise defaults ON: a fresh form (no sidecar, no stickies) checks
      # the box (a noise-free set puts J near 0 instead of near dof, which
      # reads like a broken statistic); an explicit false stays unchecked
      @test occursin("name=\"noise\" value=\"true\" checked", SparlectraApp.render_se_form(; cases=["case14.m"], selected_case="case14.m"))
      @test !occursin("name=\"noise\" value=\"true\" checked", SparlectraApp.render_se_form(; cases=["case14.m"], selected_case="case14.m", gen_values=Dict{String,String}("noise" => "false")))
      rpst = start_powerflow_run(Dict{String,Any}("casefile" => pstcase, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => pstout))
      @test rpst["status"] == "succeeded"
      @test rpst["metadata"]["se_band_reason"] == "ok"
      # the phase regulator of the PST (a nameplate Delta-u changer,
      # neutral step 0, so s and the step coincide)
      pstrow = only(t for t in rpst["metadata"]["se_tap_estimates"] if t["branch"] == 8 && t["role"] == "phase")
      @test pstrow["kind"] == "asymmetrical"
      @test pstrow["step_estimated"] == 2
      @test abs(pstrow["step_estimated_raw"] - 2.0) < 0.3
      gsurow = only(t for t in rpst["metadata"]["se_tap_estimates"] if t["branch"] == 7)
      @test gsurow["status"] == "calculated" && gsurow["reason"] == "machine_transformer"
    end)() end

    @testset "SE chain to power flow and N-1" begin (function ()
      # chain PF (:se_snapshot): immediate convergence, slack pickup recorded
      r2 = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_start_run_id" => id1, "se_start_mode" => "se_snapshot"))
      @test r2["status"] == "succeeded"
      @test r2["iterations"] <= 1
      @test abs(r2["metadata"]["slack_pickup_mw"]) < 1e-6
      @test r2["metadata"]["run_mode"] == "powerflow_se_start"
      # deviation report PF vs SE state: metadata extremes plus the per-bus
      # artifact (snapshot start on a consistent set: deviation tiny)
      @test isfinite(r2["metadata"]["se_pf_max_dvm_pu"])
      @test r2["metadata"]["se_pf_max_dvm_pu"] < 1e-6
      @test isfile(joinpath(r2["output_dir"], "se_pf_deviation.csv"))
      @test run_csv(joinpath(r2["output_dir"], "se_pf_deviation.csv")).header[1:9] == ["bus", "name", "mrid", "vm_se_pu", "vm_pf_pu", "dvm_pu", "va_se_deg", "va_pf_deg", "dva_deg"]

      # chain N-1: consumes the SE-started base case, notes the SE run id,
      # and matches the manual run on the same case
      # screening off on BOTH runs: this parity check compares the SE-started
      # base against the manual base through FULL solves (which re-converge to
      # 1e-9 agreement); screening estimates are first-order functions of the
      # base state itself and legitimately differ at the solver-tolerance
      # level, so they are not part of the chain-parity contract
      r3 = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "contingency_mode" => true, "contingency_kind" => "branch", "screening_mode" => "off", "se_start_run_id" => id1))
      @test r3["status"] == "succeeded"
      @test r3["metadata"]["se_run_id"] == id1
      r4 = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "contingency_mode" => true, "contingency_kind" => "branch", "screening_mode" => "off"))
      @test r4["status"] == "succeeded"
      # identical results: structural fields exact, numeric fields to 1e-9
      # (the SE-started base solve walks a different iteration path than the
      # flat start, so the last ULPs of the voltages differ legitimately)
      lines3 = readlines(joinpath(root, r3["run_id"], "contingency_n1.csv"))
      lines4 = readlines(joinpath(root, r4["run_id"], "contingency_n1.csv"))
      @test length(lines3) == length(lines4)
      # writeContingencyResultsCSV's default delimiter is "technical" (comma)
      # since issue #376, not the old hardcoded semicolon; a hardcoded ";"
      # here found no delimiter at all and silently degraded every row to a
      # single-field exact-string comparison, hiding the intended per-field
      # isapprox tolerance below.
      csv_delim = occursin(';', first(lines3)) ? ';' : ','
      # the note column is quoted and carries commas ("islanded: bus 8
      # cut off, 0.0 MW load ..."); a plain split at the delimiter
      # shifted the columns after it and made the number of checks
      # depend on the note text (webui 3015 -> 3014 when an N-1 note
      # changed). The quote-aware reader of the Web UI keeps one field
      # per column.
      ncols = length(SparlectraApp._webui_compare_split_csv(first(lines3), csv_delim))
      # excel_de (a direct call under a German locale) writes a decimal
      # comma and thousands dots; without this every number would fall
      # back to the exact-string comparison
      csv_number(s) = tryparse(Float64, csv_delim == ';' ? replace(s, "." => "", "," => ".") : s)
      for (l3, l4) in zip(lines3, lines4)
          f3 = SparlectraApp._webui_compare_split_csv(l3, csv_delim)
          f4 = SparlectraApp._webui_compare_split_csv(l4, csv_delim)
          @test length(f3) == length(f4) == ncols
          for (a, b) in zip(f3, f4)
              na = csv_number(a)
              nb = csv_number(b)
              if na !== nothing && nb !== nothing
                  @test (isnan(na) && isnan(nb)) || isapprox(na, nb; atol=1e-9)
              else
                  @test a == b
              end
          end
      end

      # chain without a preceding SE run errors clearly (both entry points)
      r5 = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_start_run_id" => "does-not-exist"))
      @test r5["status"] == "failed"
      # se_mode without a measurement file is rejected up front
      r6 = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true))
      @test r6["status"] == "failed"
    end)() end
  end)() end
  return nothing
end

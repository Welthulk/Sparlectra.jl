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

# file: test/test_webui.jl
# purpose: fast Web UI tests without a live server: form parsing, selectors
#          and rendering, commit sha resolution, saved case settings, and
#          the sysimage launcher decision
using Sparlectra
using Test
using Dates
using SHA

include("test_webui_support.jl")

function run_webui_fast_tests()
    @testset "Web UI fast parsing, selectors, and rendering" begin (function ()
        @test SparlectraApp._webui_parse_bool("on") === true
        @test SparlectraApp._webui_parse_bool("false") === false
        @test SparlectraApp.route_sparlectra_webui("GET", "/powerflow/artifact/..%2Fbad/result.json").status in (400, 404)

        mktempdir() do root
            case_directory = joinpath(root, "cases")
            mkpath(case_directory)
            write(joinpath(case_directory, "case14.m"), "function mpc = case14\nend\n")
            write(joinpath(case_directory, "FOR001.DAT"), _dtf_network_fixture())
            write(joinpath(case_directory, "FOR001_OUTAGES.DAT"), _dtf_network_with_outage_fixture())
            write(joinpath(case_directory, "OUTAGE.DAT"), _dtf_outage_fixture())
            write(joinpath(case_directory, "FOR002.DAT"), _dtf_reference_fixture())
            write(joinpath(case_directory, "UNKNOWN.DAT"), "plain unsupported data\n")
            @test SparlectraApp._webui_classify_dat_content(joinpath(case_directory, "FOR001.DAT")) === :dtf_network_case
            @test SparlectraApp._webui_classify_dat_content(joinpath(case_directory, "FOR001_OUTAGES.DAT")) === :dtf_network_case_with_outages
            @test SparlectraApp._webui_classify_dat_content(joinpath(case_directory, "OUTAGE.DAT")) === :dtf_outage_file
            @test SparlectraApp._webui_classify_dat_content(joinpath(case_directory, "FOR002.DAT")) === :dtf_outage_or_reference
            @test SparlectraApp._webui_classify_dat_content(joinpath(case_directory, "UNKNOWN.DAT")) === :unknown_dat
            primary = SparlectraApp._webui_casefile_options_in_directory(case_directory)
            reference = SparlectraApp._webui_for002_reference_options_in_directory(case_directory)
            @test "FOR001.DAT" in primary
            @test "FOR001_OUTAGES.DAT" in primary
            @test !("OUTAGE.DAT" in primary)
            @test !("UNKNOWN.DAT" in primary)
            @test reference == ["FOR002.DAT"]
        end

        @testset "request builder falls back to the case form block" begin (function ()
            # fields the run page no longer renders must reach the run from the
            # case configuration file's form block; a POSTed field always wins
            dir = mktempdir()
            case_path = joinpath(dir, "case_formblock.m")
            write(case_path, "function mpc = case_formblock\nend\n")
            Sparlectra.write_case_config(case_path, Dict{String,Any}(); form=Dict{String,Any}("export_cgmes" => true, "performance_timing" => "full", "case_format" => "matpower", "gen_seed" => 7))
            req = SparlectraApp.powerflow_webui_request(Dict("casefile" => case_path))
            @test req["export_cgmes"] === true
            @test req["performance_timing"] == "full"
            @test req["case_format"] == "matpower"
            req_posted = SparlectraApp.powerflow_webui_request(Dict("casefile" => case_path, "export_cgmes" => "false", "case_format" => "auto"))
            @test req_posted["export_cgmes"] === false
            @test req_posted["case_format"] == "auto"
            # relative case name resolves through case_directory
            req_rel = SparlectraApp.powerflow_webui_request(Dict("casefile" => "case_formblock.m"); case_directory=dir)
            @test req_rel["export_cgmes"] === true
            # without a stored format the hint chain answers: a free-typed .DAT
            # yields the explicit DTF format, everything else stays auto
            req_dat = SparlectraApp.powerflow_webui_request(Dict("casefile" => "no_such_FOR001.DAT"); case_directory=dir)
            @test req_dat["case_format"] == "dtf_for001"
            @test SparlectraApp.powerflow_webui_request(Dict("casefile" => "case.m"))["case_format"] == "auto"
            # the result-page save must not wipe the stored non-spec field: the
            # preservation loop carries case_format like the spec fields
            stored = SparlectraApp._webui_case_form_defaults(case_path, nothing)
            @test stored["case_format"] == "matpower"
            @test stored["gen_seed"] == 7

            # export_cgmes parsing: an unchecked checkbox submits only the
            # hidden false field, an absent field takes the spec default (false)
            @test SparlectraApp.powerflow_webui_request(Dict("casefile" => "case.m", "export_cgmes" => "true"))["export_cgmes"] === true
            @test SparlectraApp.powerflow_webui_request(Dict("casefile" => "case.m", "export_cgmes" => "false"))["export_cgmes"] === false
            @test SparlectraApp.powerflow_webui_request(Dict("casefile" => "case.m"))["export_cgmes"] === false

            # the export row never passes a partial export off as a full one:
            # notices and the zero-sequence line count show, an error is carried
            # into a failed row, and a run without export has no row
            completed_html = SparlectraApp.render_powerflow_result(
                Dict(
                    "run_id" => "x",
                    "status" => "completed",
                    "metadata" => Dict("cgmes_export_status" => "completed", "cgmes_export_files" => "n_EQ.xml, n_TP.xml, n_SSH.xml", "cgmes_export_notices" => "transformer T1: phase shift 2.0° not exported (fixed ratio only)", "cgmes_export_sc_lines" => 3),
                ),
            )
            @test occursin("n_SSH.xml", completed_html)
            @test occursin("phase shift 2.0° not exported", completed_html)
            @test occursin("zero-sequence data on 3 line(s)", completed_html)
            failed_html = SparlectraApp.render_powerflow_result(Dict("run_id" => "y", "status" => "completed", "metadata" => Dict("cgmes_export_status" => "failed", "cgmes_export_error" => "boom")))
            @test occursin("failed: boom", failed_html)
            plain_html = SparlectraApp.render_powerflow_result(Dict("run_id" => "z", "status" => "completed"))
            @test !occursin("CGMES export", plain_html)
        end)() end

        @testset "shared selected-case memory" begin (function ()
            # choosing on the Case page must reach the run page and the SE page
            # through plain nav links (no query): a live run hit "no case
            # selected" exactly this way
            root = mktempdir()
            cases = joinpath(root, "cases")
            mkpath(cases)
            # suite migration (no-download rule): this set only exercises the
            # shared-memory routing by NAME, so a fixture file suffices
            write(joinpath(cases, "case14.m"), "% case fixture\n")
            rt = (; case_directory=cases, config_file=Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log=SparlectraApp.webui_operation_log_path(root), startup_config_error=nothing, runner=SparlectraApp.start_powerflow_run)
            # nothing remembered yet: the run page renders without a case
            fresh = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow", Dict{String,String}(); output_root=root, runtime=rt).body)
            @test occursin("no case selected", fresh)
            # opening the Case page with an explicit case remembers it (GET
            # queries travel IN the target, the third argument is a POST body)
            SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case?casefile=case14.m"; output_root=root, runtime=rt)
            @test SparlectraApp._webui_recall_selected_case(root) == "case14.m"
            # the plain-nav run page now carries the remembered case as its hidden field
            run_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow", Dict{String,String}(); output_root=root, runtime=rt).body)
            @test occursin("type=\"hidden\" name=\"casefile\" value=\"case14.m\"", run_page)
            # /stateestimation is a real redirect onto the Runs page, whose SE
            # section consumes the same memory
            se_redirect = SparlectraApp.route_sparlectra_webui("GET", "/stateestimation", Dict{String,String}(); output_root=root, runtime=rt)
            @test se_redirect.status == 303
            @test Dict(se_redirect.headers)["Location"] == "/powerflow#state-estimation"
            @test occursin("id=\"state-estimation\"", run_page)
            # an empty run POST now points at the Case page, format-neutral
            err = try
                SparlectraApp.powerflow_webui_request(Dict{String,Any}())
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("Case page", sprint(showerror, err))
            @test !occursin("MATPOWER", sprint(showerror, err))
        end)() end

        @testset "settings save reaches the run without POST fields" begin (function ()
            # The riskiest single assumption of the page split: a value saved on the Settings page must reach a run whose
            # POST no longer carries the field, via resolve_config, and the
            # effective-config artifact must name case_sidecar as its source.
            root = mktempdir()
            cache = joinpath(root, "cases")
            app_root = normpath(joinpath(dirname(@__DIR__)))
            SparlectraApp._webui_stage_bundled_case!(app_root, cache, "sp_case14.scf.json")
            # the shipped CGMES deliveries are bundled as one ZIP per case,
            # packed into the case cache on first use and importable as is
            @test "sp_case14_cgmes.zip" in SparlectraApp._webui_bundled_scf_options(app_root)
            demo_zip = SparlectraApp._webui_stage_bundled_case!(app_root, cache, "sp_case14_cgmes.zip")
            @test demo_zip == joinpath(cache, "sp_case14_cgmes.zip") && isfile(demo_zip)
            @test length(importCGMES(path=demo_zip, name="sp_case14_cgmes").net.nodeVec) == 14
            @test SparlectraApp._webui_stage_bundled_case!(app_root, cache, "no_such_cgmes.zip") === nothing
            # 0.50.0: at Web UI start a copy that equals an earlier committed
            # version follows the package, a copy the user changed stays, a
            # current copy is not touched. The hash list carries the current
            # content of what staging wrote, the packed CGMES ZIP included
            # (the builder packs as the Web UI does); text files hash with
            # CRLF read as LF, so none of this depends on the line ends of
            # the checkout. The earlier versions are given here, the real
            # ones exist only in git.
            shipped = SparlectraApp._webui_shipped_hashes(app_root)
            for f in ("sp_case14.scf.json", "sp_case14.measurements.csv", "sp_case14.config.yaml", "sp_case14_cgmes.zip")
                @test SparlectraApp._webui_shipped_sha256(joinpath(cache, f)) in shipped[f]
            end
            fresh = joinpath(root, "refresh_cases")
            SparlectraApp._webui_stage_bundled_case!(app_root, fresh, "sp_case14.scf.json")
            SparlectraApp._webui_stage_bundled_case!(app_root, fresh, "sp_case14_cgmes.zip")
            # an earlier version as a Windows checkout (core.autocrlf) wrote it
            old_meas = "# an earlier shipped measurement set\nV;1;1.0\n"
            old_meas_crlf = replace(old_meas, "\n" => "\r\n")
            write(joinpath(fresh, "sp_case14.measurements.csv"), old_meas_crlf)
            write(joinpath(fresh, "sp_case14_cgmes.zip"), "an earlier shipped ZIP")
            write(joinpath(fresh, "sp_case14.measurements.baddata.csv"), "edited by the user\n")
            # the current configuration file in CRLF counts as current
            cfg_copy = joinpath(fresh, "sp_case14.config.yaml")
            cfg_crlf = replace(replace(read(cfg_copy, String), "\r\n" => "\n"), "\n" => "\r\n")
            write(cfg_copy, cfg_crlf)
            case_mtime = mtime(joinpath(fresh, "sp_case14.scf.json"))
            cfg_mtime = mtime(cfg_copy)
            known = Dict(name => copy(hashes) for (name, hashes) in shipped)
            push!(known["sp_case14.measurements.csv"], bytes2hex(sha256(old_meas)))
            push!(known["sp_case14_cgmes.zip"], bytes2hex(sha256("an earlier shipped ZIP")))
            refresh_log = joinpath(root, "refresh_log")
            day = Date(2026, 10, 18)
            replaced = SparlectraApp._webui_refresh_shipped_copies!(app_root, fresh, refresh_log; known, today = day)
            @test sort([r.file for r in replaced]) == ["sp_case14.measurements.csv", "sp_case14_cgmes.zip"]
            @test read(joinpath(fresh, "sp_case14.measurements.csv")) == read(joinpath(app_root, "data", "scf", "sp_case14.measurements.csv"))
            @test read(joinpath(fresh, "sp_case14.measurements.csv.replaced-2026-10-18"), String) == old_meas_crlf
            @test read(joinpath(fresh, "sp_case14_cgmes.zip")) == read(joinpath(cache, "sp_case14_cgmes.zip"))
            @test read(joinpath(fresh, "sp_case14.measurements.baddata.csv"), String) == "edited by the user\n"
            @test mtime(joinpath(fresh, "sp_case14.scf.json")) == case_mtime
            @test mtime(cfg_copy) == cfg_mtime && read(cfg_copy, String) == cfg_crlf
            @test count(f -> occursin(".replaced-", f), readdir(fresh)) == 2
            log_text = read(SparlectraApp.webui_operation_log_path(refresh_log), String)
            @test count("shipped_case_file_replaced", log_text) == 2
            @test count("shipped_case_file_replacing", log_text) == 2
            @test occursin(bytes2hex(sha256(old_meas)), log_text)
            # nothing left to do on the next start
            @test isempty(SparlectraApp._webui_refresh_shipped_copies!(app_root, fresh, refresh_log; known, today = day))
            # a replacement that fails while writing (full disk) leaves the
            # copy as it was and no temporary or kept file behind; an
            # unexpected failure is rethrown, with the same cleanup
            victim = joinpath(fresh, "sp_case5.measurements.csv")
            write(victim, old_meas)
            names_before = sort(readdir(fresh))
            disk_full = d -> (write(d, "half a fi"); throw(Base.IOError("write: no space left on device (ENOSPC)", -28)))
            failed = SparlectraApp._webui_stage_shipped_file!(disk_full, () -> "new", victim; known = Set([bytes2hex(sha256(old_meas))]), operation_log = refresh_log, today = day)
            @test failed.status === :failed && failed.kept_as == ""
            @test read(victim, String) == old_meas
            @test sort(readdir(fresh)) == names_before
            @test count("shipped_case_file_replace_failed", read(SparlectraApp.webui_operation_log_path(refresh_log), String)) == 1
            @test_throws ErrorException SparlectraApp._webui_stage_shipped_file!(d -> error("bug"), () -> "new", victim; known = Set([bytes2hex(sha256(old_meas))]), operation_log = refresh_log, today = day)
            @test read(victim, String) == old_meas && sort(readdir(fresh)) == names_before
            # a damaged list (wrong value type, no TOML) is named as such
            bad_root = joinpath(root, "bad_root")
            mkpath(joinpath(bad_root, "data"))
            write(joinpath(bad_root, "data", "shipped_hashes.toml"), "[files]\n\"a.csv\" = 3\n")
            @test_throws SparlectraApp._WebUIShippedHashesDamaged SparlectraApp._webui_shipped_hashes(bad_root)
            write(joinpath(bad_root, "data", "shipped_hashes.toml"), "[files\n")
            @test_throws SparlectraApp._WebUIShippedHashesDamaged SparlectraApp._webui_shipped_hashes(bad_root)
            # the runtime gets its own configuration copy: a save with the
            # case target writes machine-scope keys to the configuration file,
            # and the packaged template must never be that file
            cfg_rt = joinpath(root, "rt.configuration.yaml")
            cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, cfg_rt)
            rt = (; case_directory=cache, config_file=cfg_rt, operation_log=SparlectraApp.webui_operation_log_path(root), startup_config_error=nothing, runner=SparlectraApp.start_powerflow_run)
            resp = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "power_flow_max_iter" => "44"); output_root=root, runtime=rt)
            @test resp.status in (302, 303)
            @test Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json"))["power_flow.max_iter"] == 44
            # the block-3 run form posts NO override fields; the bare request
            # must carry none either
            req = SparlectraApp.powerflow_webui_request(Dict("casefile" => "sp_case14.scf.json"); case_directory=cache)
            @test isempty(req["config_overrides"])
            # an unchecked APSLF checkbox posts an explicit "false" (hidden plus
            # checkbox pair) and arrives as false; a field absent from the form
            # (a stale client) is left out instead of defaulting to false
            apslf_req = SparlectraApp.powerflow_webui_request(Dict{String,Any}("casefile" => "sp_case14.scf.json", "power_flow_apslf_nr_polish" => "false"); case_directory=cache)
            @test apslf_req["config_overrides"]["power_flow.apslf.nr_polish"] === false
            @test !haskey(apslf_req["config_overrides"], "power_flow.apslf_start.enabled")
            run = SparlectraApp.start_powerflow_run(merge(req, Dict("output_root" => joinpath(root, "runs"), "config_file" => cfg_rt)); case_directory=cache)
            @test run["status"] == "succeeded"
            eff = read(joinpath(String(run["output_dir"]), "effective_config.yaml"), String)
            seg = eff[first(findfirst("  power_flow:", eff)):end]
            mi = seg[first(findfirst("    max_iter:", seg)):(first(findfirst("    max_iter:", seg))+220)]
            @test occursin("value: 44", mi)
            @test occursin("source: case_sidecar", mi)
            # the solver choice wins over a generator toggle saved earlier: the
            # disabled toggle is not posted, so without this the sidecar kept
            # apslf_start.enabled = true under solver = apslf and the next run
            # failed (a Web UI run)
            resp_gen = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "power_flow_apslf_start_enabled" => "true"); output_root=root, runtime=rt)
            @test resp_gen.status in (302, 303)
            @test Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json"))["power_flow.apslf_start.enabled"] === true
            resp_solver = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "power_flow_solver" => "apslf"); output_root=root, runtime=rt)
            @test resp_solver.status in (302, 303)
            saved_cfg = Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json"))
            @test saved_cfg["power_flow.solver"] == "apslf"
            # the forced false is the inherited value, so the stored true is removed
            @test get(saved_cfg, "power_flow.apslf_start.enabled", false) === false
            @test Sparlectra.resolve_config(cfg_rt, joinpath(cache, "sp_case14.scf.json")).config.powerflow.apslf_start.enabled === false
            # an incompatible pair in one save is refused at save time, not at the run
            resp_bad = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "power_flow_solver" => "rectangular", "power_flow_apslf_start_enabled" => "true", "power_flow_dc_seed_unconditional" => "true"); output_root=root, runtime=rt)
            @test occursin("Could not save settings for this case", SparlectraApp._webui_urldecode(Dict(resp_bad.headers)["Location"]))
            @test Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json"))["power_flow.solver"] == "apslf"
            # the staged knees out of order (start above end) are refused at
            # save time as well, in the words of the form, and nothing is written
            resp_knees = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "se_robust_k1" => "10", "se_robust_k2" => "6"); output_root=root, runtime=rt)
            knees_note = SparlectraApp._webui_urldecode(Dict(resp_knees.headers)["Location"])
            @test occursin("Could not save settings", knees_note)
            @test occursin("down-weight from |r|/sigma", knees_note)
            @test !occursin("se_robust_k1: 10", read(Sparlectra.case_config_path(joinpath(cache, "sp_case14.scf.json")), String))
            SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "power_flow_solver" => "rectangular"); output_root=root, runtime=rt)
            # the flat start is a flat start: the saved start
            # settings stay as posted, the run switches the projection, the
            # DC seed and the two start modes off (named in run.log), the
            # current-iteration and APSLF starts follow their own keys;
            # unchecking gives the modes back. The ratio-profile switch
            # travels from the form into the case file.
            resp_flat = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "power_flow_flatstart" => "true", "power_flow_ratio_profile" => "false", "power_flow_start_current_iteration_enabled" => "true", "power_flow_start_angle_mode" => "dc", "power_flow_start_voltage_mode" => "profile_blend"); output_root=root, runtime=rt)
            @test resp_flat.status in (302, 303)
            flat_cfg = Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json"))
            @test flat_cfg["power_flow.flatstart"] === true
            @test flat_cfg["power_flow.start_mode.ratio_profile"] === false
            @test flat_cfg["power_flow.start_current_iteration.enabled"] === true
            run_flat = SparlectraApp.start_powerflow_run(Dict("casefile" => "sp_case14.scf.json", "config_file" => cfg_rt, "output_root" => root); case_directory=cache)
            @test run_flat["status"] == "succeeded"
            run_log = read(joinpath(String(run_flat["output_dir"]), "run.log"), String)
            @test occursin("power_flow.start_mode.angle_mode=classic, power_flow.start_mode.voltage_mode=classic", run_log)
            @test occursin(r"Flatstart\s+:\s+Yes", run_log)
            # the pre-solve was not switched off by the flat start
            @test run_flat["metadata"]["current_iteration_enabled"] === true
            profile_cfg, profile_keys = Sparlectra._flatstart_profile_config(Sparlectra.load_sparlectra_config(cfg_rt; reload = true, overrides = Dict{String,Any}("power_flow" => Dict{String,Any}("flatstart" => true, "apslf_start" => Dict{String,Any}("enabled" => true), "start_mode" => Dict{String,Any}("start_projection" => true, "angle_mode" => "dc")))))
            @test profile_cfg.powerflow.start_mode.start_projection === false
            @test profile_cfg.powerflow.apslf_start.enabled === true
            @test profile_cfg.powerflow.start_mode.angle_mode === :classic
            @test "power_flow.start_mode.angle_mode=classic" in profile_keys
            @test "power_flow.start_mode.start_projection=false" in profile_keys
            @test !any(k -> occursin("apslf_start", k), profile_keys)
            SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "power_flow_flatstart" => "false"); output_root=root, runtime=rt)
            # both are back at their inherited values (L20: such a key is not
            # kept in the case file), and the run reads them that way
            off_pf = Sparlectra.resolve_config(cfg_rt, joinpath(cache, "sp_case14.scf.json")).config.powerflow
            @test off_pf.start_mode.flatstart === false
            @test off_pf.start_mode.angle_mode === :dc
            @test !haskey(Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json")), "power_flow.flatstart")
            # machine-scope keys are named and kept out of the case file
            resp2 = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "output_logfile_results" => "compact"); output_root=root, runtime=rt)
            @test occursin("output.logfile_results", SparlectraApp._webui_urldecode(Dict(resp2.headers)["Location"]))
            @test !haskey(Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json")), "output.logfile_results")
            # the general target merges into the YAML with a backup
            cfg = joinpath(root, "configuration.yaml")
            cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, cfg)
            resp3 = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("settings_target" => "general", "config_file" => cfg, "power_flow_autodamp_min" => "0.09"); output_root=root)
            @test resp3.status in (302, 303)
            @test occursin("autodamp_min: 0.09", read(cfg, String))
            @test isfile(cfg * ".settings-save.bak")

            # the Q-limit block (0.20.2): values saved per case and in the
            # configuration file come back in the form and reach the run, and
            # a changed guard value shows on the comparison page
            eff_entry(dir, key) = begin
              eff_q = read(joinpath(String(dir), "effective_config.yaml"), String)
              at = findfirst(string(key, ":"), eff_q)
              at === nothing ? "" : eff_q[first(at):min(end, first(at) + 400)]
            end
            form_value(html, field) = (m = match(Regex("name=\"$(field)\"[^>]*value=\"([^\"]*)\""), html); m === nothing ? nothing : m.captures[1])
            resp_q = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}(
                "casefile" => "sp_case14.scf.json", "settings_target" => "this_case",
                "power_flow_qlimits_guard_max_switches" => "5",
                "power_flow_qlimits_guard_freeze_after_repeated_switching" => "false",
                "power_flow_qlimits_start_iter" => "4",
                "power_flow_qlimits_hysteresis_pu" => "2e-2",
                "power_flow_qlimits_trace_buses" => "1, 2",
                "power_flow_qlimits_guard_violation_mode" => "delayed_switch",
                "power_flow_qlimits_final_q_accept_pu" => "0.05"); output_root=root, runtime=rt)
            @test resp_q.status in (302, 303)
            q_cfg = Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json"))
            @test q_cfg["power_flow.qlimits.guard.max_switches"] == 5
            @test q_cfg["power_flow.qlimits.guard.freeze_after_repeated_switching"] === false
            @test q_cfg["power_flow.qlimits.start_iter"] == 4
            @test q_cfg["power_flow.qlimits.hysteresis_pu"] == 0.02
            @test q_cfg["power_flow.qlimits.trace_buses"] == [1, 2]
            @test q_cfg["power_flow.qlimits.guard.violation_mode"] == "delayed_switch"
            # the final check bound, read by every mode (the classic ones too):
            # a number is stored as a number, the word auto as the word
            @test q_cfg["power_flow.qlimits.final_q_accept_pu"] == 0.05
            @test SparlectraApp._webui_parse_form_value("auto", Union{Float64,String}, "f") == "auto"
            @test_throws ArgumentError SparlectraApp._webui_parse_form_value("wide", Union{Float64,String}, "f")
            q_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/settings?casefile=sp_case14.scf.json"; output_root=root, runtime=rt).body)
            @test form_value(q_page, "power_flow_qlimits_guard_max_switches") == "5"
            @test form_value(q_page, "power_flow_qlimits_start_iter") == "4"
            @test form_value(q_page, "power_flow_qlimits_trace_buses") == "1, 2"
            @test form_value(q_page, "power_flow_qlimits_final_q_accept_pu") == "0.05"
            @test !occursin(r"name=\"power_flow_qlimits_guard_freeze_after_repeated_switching\" type=\"checkbox\"[^>]*checked", q_page)
            run_q = SparlectraApp.start_powerflow_run(Dict("casefile" => "sp_case14.scf.json", "config_file" => cfg_rt, "output_root" => root); case_directory=cache)
            @test run_q["status"] == "succeeded"
            @test occursin("value: 5", eff_entry(run_q["output_dir"], "max_switches"))
            @test occursin("source: case_sidecar", eff_entry(run_q["output_dir"], "max_switches"))
            @test occursin("value: false", eff_entry(run_q["output_dir"], "freeze_after_repeated_switching"))
            # the configuration file: a case without a case configuration file
            # takes the general value (with one, the case file wins)
            cp(joinpath(Sparlectra.SPARLECTRA_ROOT, "data", "mpower", "sp_case9.m"), joinpath(cache, "sp_case9.m"))
            resp_g = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("settings_target" => "general", "config_file" => cfg, "power_flow_qlimits_guard_max_switches" => "7"); output_root=root)
            @test resp_g.status in (302, 303)
            @test Sparlectra.load_sparlectra_config(cfg; reload = true).powerflow.qlimits.guard_max_switches == 7
            rt_g = (; case_directory=cache, config_file=cfg, operation_log=SparlectraApp.webui_operation_log_path(root), startup_config_error=nothing, runner=SparlectraApp.start_powerflow_run)
            g_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/settings?casefile=sp_case9.m"; output_root=root, runtime=rt_g).body)
            @test form_value(g_page, "power_flow_qlimits_guard_max_switches") == "7"
            run_g = SparlectraApp.start_powerflow_run(Dict("casefile" => "sp_case9.m", "config_file" => cfg, "output_root" => root); case_directory=cache)
            @test run_g["status"] == "succeeded"
            @test occursin("value: 7", eff_entry(run_g["output_dir"], "max_switches"))
            @test ("power_flow.qlimits.guard.max_switches", "5", "7") in SparlectraApp._webui_compare_config_diff(run_q["output_dir"], run_g["output_dir"])
        end)() end

        @testset "the provisioned configuration follows changed template defaults" begin (function ()
            # a copy of the template provisioned by an older release still
            # carries the old default of a key the user never touched; it
            # follows the new template value, a value the user set stays
            mktempdir() do dir
                cfg = joinpath(dir, "configuration.yaml")
                text = read(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, String)
                @test occursin("linear_solver: umfpack_reuse", text)
                # an old copy: the old defaults of the solver and of the CSV
                # format (technical before 0.32.0), a user choice for the tolerance
                write(cfg, replace(text, "linear_solver: umfpack_reuse" => "linear_solver: umfpack", "csv_format: auto" => "csv_format: technical", r"^  tol: [^\n]*"m => "  tol: 1.0e-7"))
                changed = SparlectraApp._webui_follow_template_defaults!(cfg)
                # the report names key, old and new value for the start message
                @test any(c -> c.key == "power_flow.linear_solver" && c.old == "umfpack" && c.new == "umfpack_reuse", changed)
                @test any(c -> c.key == "output.csv_format" && c.old == "technical" && c.new == "auto", changed)
                @test Sparlectra.load_sparlectra_config(cfg; reload = true).powerflow.linear_solver === :umfpack_reuse
                @test Sparlectra.load_sparlectra_config(cfg; reload = true).powerflow.tol == 1.0e-7
                @test isfile(joinpath(dir, "configuration.template.yaml"))
                @test isfile(string(cfg, ".template-follow.bak"))
                # aligned: a second pass changes nothing
                @test isempty(SparlectraApp._webui_follow_template_defaults!(cfg))
                # with the template copy in place a changed template default follows,
                # a user value that differs from the old template stays
                template_copy = joinpath(dir, "configuration.template.yaml")
                write(template_copy, replace(read(template_copy, String), "linear_solver: umfpack_reuse" => "linear_solver: umfpack"))
                write(cfg, replace(read(cfg, String), "linear_solver: umfpack_reuse" => "linear_solver: umfpack"))
                @test [c.key for c in SparlectraApp._webui_follow_template_defaults!(cfg)] == ["power_flow.linear_solver"]
                write(template_copy, replace(read(template_copy, String), "linear_solver: umfpack_reuse" => "linear_solver: umfpack"))
                @test isempty(SparlectraApp._webui_follow_template_defaults!(cfg))   # user value equals the new template already
                @test Sparlectra.load_sparlectra_config(cfg; reload = true).powerflow.tol == 1.0e-7
                # a key the Web UI saved is the user's even when its value equals
                # the old template default: it is not followed
                write(template_copy, replace(read(template_copy, String), "linear_solver: umfpack_reuse" => "linear_solver: umfpack"))
                write(cfg, replace(read(cfg, String), "linear_solver: umfpack_reuse" => "linear_solver: umfpack"))
                SparlectraApp._webui_record_user_keys!(cfg, ("power_flow.linear_solver",))
                @test isempty(SparlectraApp._webui_follow_template_defaults!(cfg))
                @test Sparlectra.load_sparlectra_config(cfg; reload = true).powerflow.linear_solver === :umfpack
                # a settings save records its keys, so a saved value survives the next start
                general = SparlectraApp._webui_write_general_settings!(cfg, Dict{String,Any}("power_flow.max_iter" => 77))
                @test general.ok
                @test "power_flow.max_iter" in SparlectraApp._webui_user_keys(cfg)
            end
        end)() end

        @testset "the settings page shows the case's saved values; the configuration file view is a link" begin (function ()
            # Reversing the earlier default: with a case
            # selected every page shows the values a run of this case will use,
            # so the values just saved are visible right after "Saved settings
            # for this case" (before, the page showed the configuration file's
            # old values next to that message). The file view stays reachable
            # as a link (?config_view=1, also case_settings=0).
            dir = mktempdir()
            cache = joinpath(dir, "cases")
            mkpath(cache)
            root = joinpath(dir, "out")
            mkpath(root)
            cp(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")), joinpath(cache, "sp_case14.scf.json"))
            cfg = joinpath(root, "configuration.yaml")
            cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, cfg)
            rt = (; case_directory = cache, config_file = cfg, operation_log = SparlectraApp.webui_operation_log_path(root), startup_config_error = nothing, runner = SparlectraApp.start_powerflow_run)
            page(q) = String(copy(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/settings?casefile=sp_case14.scf.json" * q; output_root = root, runtime = rt).body))
            max_iter(html) = (m = match(r"name=\"power_flow_max_iter\"[^>]*value=\"(\d+)\"", html); m === nothing ? "" : String(m.captures[1]))
            saved = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "config_file" => cfg, "settings_target" => "this_case", "power_flow_max_iter" => "33"); output_root = root, runtime = rt)
            @test saved.status in (302, 303)
            @test Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json"))["power_flow.max_iter"] == 33
            cfg_max = string(Sparlectra.load_sparlectra_config(cfg; reload = true).powerflow.max_iter)
            default_view = page("")
            @test max_iter(default_view) == "33"
            @test occursin("case-settings-notice", default_view)
            @test occursin("config_view=1", default_view)
            for q in ("&config_view=1", "&case_settings=0")
                file_view = page(q)
                @test max_iter(file_view) == cfg_max
                @test !occursin("case-settings-notice", file_view)
                @test occursin("case_settings=1", file_view)
            end
            # the save message shows in both views
            @test occursin("Saved settings probe", page("&save_message=Saved%20settings%20probe"))
            @test occursin("Saved settings probe", page("&config_view=1&save_message=Saved%20settings%20probe"))
            # a save to the configuration file: the file view shows the new value,
            # the default view keeps the case's value and the message says why
            general = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "config_file" => cfg, "settings_target" => "general", "power_flow_max_iter" => "44"); output_root = root, runtime = rt)
            @test general.status in (302, 303)
            @test occursin("This case overrides 1 of these setting(s)", SparlectraApp._webui_urldecode(Dict(general.headers)["Location"]))
            @test Sparlectra.load_sparlectra_config(cfg; reload = true).powerflow.max_iter == 44
            @test max_iter(page("")) == "33"
            @test max_iter(page("&config_view=1")) == "44"
            # L20: the page posts every field it renders. Saving it unchanged
            # leaves the case file and the configuration file as they were
            # (it wrote about 90 keys into the case file before).
            case_path = joinpath(cache, "sp_case14.scf.json")
            case_cfg = Sparlectra.case_config_path(case_path)
            case_before = read(case_cfg, String)
            cfg_before = read(cfg, String)
            post(form) = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", form; output_root = root, runtime = rt)
            unchanged = _webui_browser_form(page(""), "settings-form")
            @test length(unchanged) > 50
            @test post(unchanged).status == 303
            @test read(case_cfg, String) == case_before
            @test read(cfg, String) == cfg_before
            # a changed field is written, and only that one
            @test post(merge(unchanged, Dict{String,Any}("power_flow_tol" => "1.0e-7"))).status == 303
            @test Sparlectra.load_case_config(case_path) == Dict{String,Any}("power_flow.max_iter" => 33, "power_flow.tol" => 1.0e-7)
            @test read(cfg, String) == cfg_before
            # set back to the inherited value the key leaves the case file; with
            # a case file the inherited value is the packaged default, not the
            # configuration file's 44, and the page shows what the run uses
            packaged_max = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true).powerflow.max_iter
            @test packaged_max != 44
            @test post(merge(_webui_browser_form(page(""), "settings-form"), Dict{String,Any}("power_flow_max_iter" => string(packaged_max)))).status == 303
            @test Sparlectra.load_case_config(case_path) == Dict{String,Any}("power_flow.tol" => 1.0e-7)
            @test Sparlectra.resolve_config(cfg, case_path).config.powerflow.max_iter == packaged_max
            @test max_iter(page("")) == string(packaged_max)
            # the last key set back would delete the file and hand the case to
            # the configuration file (44): the posted value stays instead
            packaged_tol = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true).powerflow.tol
            @test post(merge(_webui_browser_form(page(""), "settings-form"), Dict{String,Any}("power_flow_tol" => string(packaged_tol)))).status == 303
            @test Sparlectra.resolve_config(cfg, case_path).config.powerflow.max_iter == packaged_max
            @test Sparlectra.load_case_config(case_path) == Dict{String,Any}("power_flow.max_iter" => packaged_max)
            # A case WITHOUT a case file inherits the configuration file (44,
            # not the packaged default). Unchanged saves from the Settings
            # page, the Case page and the SE section create no case file; it
            # held max_iter 44 before and the case stopped following the
            # configuration file.
            for f in ("sp_case5.scf.json", "sp_case5.measurements.csv")
                cp(abspath(joinpath(dirname(@__DIR__), "data", "scf", f)), joinpath(cache, f))
            end
            case5 = joinpath(cache, "sp_case5.scf.json")
            get5(path) = String(copy(SparlectraApp.route_sparlectra_webui("GET", path * "?casefile=sp_case5.scf.json"; output_root = root, runtime = rt).body))
            @test max_iter(get5("/powerflow/settings")) == "44"
            for (path, id, marker) in (("/powerflow/settings", "settings-form", ""), ("/powerflow/case", "case-options-form", ""), ("/powerflow", "", "name=\"return_to\" value=\"runs\""))
                form5 = _webui_browser_form(get5(path), id; containing = marker)
                @test haskey(form5, "measurement_file") == (path == "/powerflow")
                @test post(form5).status == 303
                @test !isfile(Sparlectra.case_config_path(case5))
            end
            # a changed field creates the case file with that key only
            @test post(merge(_webui_browser_form(get5("/powerflow/settings"), "settings-form"), Dict{String,Any}("power_flow_max_iter" => "45"))).status == 303
            @test Sparlectra.load_case_config(case5) == Dict{String,Any}("power_flow.max_iter" => 45)
            @test Sparlectra.resolve_config(cfg, case5).config.powerflow.max_iter == 45
            # a measurement set chosen in the SE section survives a Settings save
            # (that save posts no set and dropped the stored choice)
            cp(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.measurements.baddata.csv")), joinpath(cache, "sp_case5.measurements.baddata.csv"))
            se5 = _webui_browser_form(get5("/powerflow"), ""; containing = "name=\"return_to\" value=\"runs\"")
            @test post(merge(se5, Dict{String,Any}("measurement_file" => "sp_case5.measurements.baddata.csv"))).status == 303
            @test post(merge(_webui_browser_form(get5("/powerflow/settings"), "settings-form"), Dict{String,Any}("power_flow_max_iter" => "46"))).status == 303
            @test occursin("<option value=\"sp_case5.measurements.baddata.csv\" selected>", get5("/powerflow"))
            @test Sparlectra.load_case_config(case5) == Dict{String,Any}("power_flow.max_iter" => 46)
            # The save that creates a case
            # file also takes over the case-scope keys in which the
            # configuration file differs from the packaged defaults. A case file
            # makes resolve_config skip the configuration file's case-scope
            # keys, so a save of tol alone moved max_iter from the
            # installation's 30 to the packaged 80. Machine keys stay out.
            rm(Sparlectra.case_config_path(case5))
            @test packaged_max != 30
            @test SparlectraApp._webui_write_general_settings!(cfg, Dict{String,Any}("power_flow.max_iter" => 30, "output.console_max_rows" => 7)).ok
            @test Sparlectra.resolve_config(cfg, case5).config.powerflow.max_iter == 30
            unchanged5 = _webui_browser_form(get5("/powerflow/settings"), "settings-form")
            @test post(unchanged5).status == 303
            @test !isfile(Sparlectra.case_config_path(case5))
            created = post(merge(unchanged5, Dict{String,Any}("power_flow_tol" => "1.0e-7")))
            @test created.status == 303
            @test occursin("1 value(s) taken from configuration.yaml", SparlectraApp._webui_urldecode(Dict(created.headers)["Location"]))
            @test Sparlectra.load_case_config(case5) == Dict{String,Any}("power_flow.tol" => 1.0e-7, "power_flow.max_iter" => 30)
            @test Sparlectra.resolve_config(cfg, case5).config.powerflow.max_iter == 30
            run5 = SparlectraApp.start_powerflow_run(Dict{String,Any}("casefile" => "sp_case5.scf.json", "config_file" => cfg, "output_root" => joinpath(root, "runs5")); case_directory = cache)
            @test run5["status"] == "succeeded"
            eff5 = read(joinpath(String(run5["output_dir"]), "effective_config.yaml"), String)
            seg5 = eff5[first(findfirst("  power_flow:", eff5)):end]
            @test occursin("value: 30", seg5[first(findfirst("    max_iter:", seg5)):(first(findfirst("    max_iter:", seg5))+220)])
            # an existing case file is not filled up: a later change of the
            # configuration file stays out of it
            @test SparlectraApp._webui_write_general_settings!(cfg, Dict{String,Any}("power_flow.autodamp_min" => 0.07)).ok
            @test post(merge(_webui_browser_form(get5("/powerflow/settings"), "settings-form"), Dict{String,Any}("power_flow_max_iter" => "31"))).status == 303
            @test Sparlectra.load_case_config(case5) == Dict{String,Any}("power_flow.tol" => 1.0e-7, "power_flow.max_iter" => 31)
            # the measurement generator creates the case file too (it stores
            # its options) and takes over the same values
            rm(Sparlectra.case_config_path(case5))
            gen5 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "sp_case5.scf.json", "config_file" => cfg); output_root = root, runtime = rt)
            @test gen5.status == 303
            @test occursin("2 value(s) taken from configuration.yaml", SparlectraApp._webui_urldecode(Dict(gen5.headers)["Location"]))
            stored5 = Sparlectra.load_case_config(case5)
            @test (stored5["power_flow.max_iter"], stored5["power_flow.autodamp_min"]) == (30, 0.07)
            @test !haskey(stored5, "output.console_max_rows")
            @test Sparlectra.resolve_config(cfg, case5).config.powerflow.max_iter == 30
            # without a case (fresh root, nothing recalled) the page shows the
            # configuration file and offers no file-view link
            root2 = joinpath(dir, "out2")
            mkpath(root2)
            rt2 = (; case_directory = cache, config_file = cfg, operation_log = SparlectraApp.webui_operation_log_path(root2), startup_config_error = nothing, runner = SparlectraApp.start_powerflow_run)
            nocase = String(copy(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/settings"; output_root = root2, runtime = rt2).body))
            @test max_iter(nocase) == "30"   # the configuration file value set above
            @test !occursin("config_view=1", nocase)
            # the form state takes a field the case sets from the case, however
            # new the configuration file is, and the rest from the configuration
            rank_cfg = joinpath(dir, "conf_rank.yaml")
            write(rank_cfg, "config_version: 1\nscope: general\npower_flow:\n  max_iter: 99\n")
            @test SparlectraApp.webui_form_state(selected_config_file=rank_cfg, sidecar_profile=Dict{String,Any}("power_flow_max_iter" => 55, "_profile_path" => joinpath(dir, "case57.config.yaml")))["power_flow_max_iter"] == 55
            @test SparlectraApp.webui_form_state(selected_config_file=rank_cfg)["power_flow_max_iter"] == 99
        end)() end

        @testset "no nested forms anywhere on the run page" begin (function ()
            # HTML forbids nested <form>. A browser closes the outer form where the
            # inner one starts, so every control AFTER it — including "Start
            # PowerFlow run" — silently falls out of the form and does nothing when
            # clicked. This bit exactly once; the depth check keeps it from
            # returning through any future in-form button.
            dir = mktempdir()
            prof = joinpath(dir, "c.sparlectra-webui.yaml")
            write(prof, "placeholder")
            for html in (
                SparlectraApp.render_powerflow_form(output_root=mktempdir()),
                SparlectraApp.render_powerflow_form(output_root=mktempdir(), selected_casefile="c.m", case_profile=Dict{String,Any}("power_flow_solver" => "dc", "_profile_path" => prof)),
            )
                depth = 0
                maxdepth = 0
                for m in eachmatch(r"</?form\b"i, html)
                    depth += startswith(lowercase(m.match), "</") ? -1 : 1
                    maxdepth = max(maxdepth, depth)
                end
                @test maxdepth == 1
                @test depth == 0
                # and the submit button must still sit inside the run form
                run_form_start = first(findfirst("<form id=\"powerflow-run-form\"", html))
                run_form_end = first(findnext("</form>", html, run_form_start))
                @test occursin("Start PowerFlow run", html[run_form_start:run_form_end])
            end
        end)() end

        @testset "saved case settings can be reset from the form" begin (function ()
            # Saved settings outrank the configuration for their keys, so a stale
            # sidecar can pin a case to a setting the user cannot override in the
            # form (measured: a delivery stuck on power_flow_solver: dc). The reset
            # path must be reachable independently of the dismissible notice.
            dir = mktempdir()
            prof = joinpath(dir, "c.sparlectra-webui.yaml")
            write(prof, "placeholder")
            with_sidecar = SparlectraApp.render_settings_page(output_root=mktempdir(), selected_casefile="c.m", case_profile=Dict{String,Any}("power_flow_solver" => "dc", "_profile_path" => prof))
            @test occursin("/powerflow/case-settings/reset", with_sidecar)
            # The solver radios must not sit inside a container the page grays
            # under the DC solver (data-ac-only-field, data-nr-only-field): the
            # script disables every control of such a container, a disabled
            # control is dropped from the POST, and the save then kept the
            # previously saved solver (a Web UI run: DC chosen, APSLF solved;
            # the delivery "stuck on dc" above was the same trap the other way
            # round). Walked through the tag stack, not grepped, because the
            # nesting is what decides.
            ancestors = let html = with_sidecar
                target = findfirst("name=\"power_flow_solver\" value=\"dc\"", html)
                voids = Set(["input", "br", "img", "meta", "link", "hr"])
                stack = String[]
                for m in eachmatch(r"<(/?)([a-zA-Z0-9]+)([^>]*)>", html)
                    m.offset >= first(target) && break
                    closing, tag, attrs = m.captures
                    tag = lowercase(tag)
                    tag in voids && continue
                    if closing == "/"
                        i = findlast(x -> startswith(x, tag * " "), stack)
                        i === nothing || deleteat!(stack, i:length(stack))
                    else
                        push!(stack, tag * " " * attrs)
                    end
                end
                stack
            end
            @test !isempty(ancestors)
            @test !any(a -> occursin("data-ac-only-field", a) || occursin("data-nr-only-field", a), ancestors)
            # without saved settings there is nothing to reset
            @test !occursin("/powerflow/case-settings/reset", SparlectraApp.render_settings_page(output_root=mktempdir()))

            # The handler deletes the sidecar and keeps the case file.
            root = joinpath(dir, "runs")
            cases = joinpath(dir, "cases")
            mkpath(root)
            mkpath(cases)
            write(joinpath(cases, "c.m"), "function mpc = c\nend\n")
            sc = SparlectraApp._webui_case_settings_path(root, "c.m"; case_directory=cases)
            mkpath(dirname(sc))
            write(sc, "values:\n  power_flow_solver: dc\n")
            response = SparlectraApp.handle_powerflow_case_settings_reset(Dict("casefile" => "c.m"); output_root=root, case_directory=cases, operation_log=root)
            @test response.status == 303
            @test !isfile(sc)
            @test isfile(joinpath(cases, "c.m"))
            # idempotent: a second reset is a no-op, not an error
            @test SparlectraApp.handle_powerflow_case_settings_reset(Dict("casefile" => "c.m"); output_root=root, case_directory=cases, operation_log=root).status == 303
            # path traversal is rejected
            @test SparlectraApp.handle_powerflow_case_settings_reset(Dict("casefile" => "../evil.m"); output_root=root, case_directory=cases, operation_log=root).status == 303
        end)() end

        @testset "operation log: clear from the page" begin (function ()
            root = mktempdir()
            log = SparlectraApp.webui_operation_log_path(root)
            rt = (; case_directory=mktempdir(), config_file=Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log=log, startup_config_error=nothing, runner=SparlectraApp.start_powerflow_run)
            SparlectraApp.record_webui_operation!(log, "probe"; route="/x", method="GET")
            SparlectraApp.record_webui_operation!(log, "probe"; route="/y", method="GET")
            page = String(copy(SparlectraApp.route_sparlectra_webui("GET", "/webui/operation-log"; output_root=root, runtime=rt).body))
            @test occursin("/webui/operation-log/clear", page)
            @test occursin("entries,", page)                      # the size is stated on the page
            before = length(readlines(log))          ## opening the page logs one entry too
            resp = SparlectraApp.route_sparlectra_webui("POST", "/webui/operation-log/clear"; output_root=root, runtime=rt)
            @test resp.status == 303
            # the file is emptied but keeps ONE entry recording the deletion, so the
            # log never becomes silently empty
            remaining = readlines(log)
            @test length(remaining) == 1
            @test occursin("operation_log_cleared", remaining[1])
            @test occursin("\"removed_entries\":$(before)", replace(remaining[1], " " => ""))
        end)() end

        @testset "case download follows a symlinked case directory" begin (function ()
            # the Web UI state directory is reachable through a symlink (a Flatpak
            # app data path pointing at ~/.local/state), so the form can carry the
            # linked spelling while the runtime holds the resolved one
            root = mktempdir()
            cases = joinpath(root, "cases")
            mkpath(cases)
            cp(abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m")), joinpath(cases, "warmup_casePST.m"))
            linked = joinpath(root, "linked")
            symlink(cases, linked)
            rt = (; case_directory=cases, config_file=Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log=SparlectraApp.webui_operation_log_path(root), startup_config_error=nothing, runner=SparlectraApp.start_powerflow_run)
            via_link = joinpath(linked, "warmup_casePST.m")
            dl = SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case/download?case=$(SparlectraApp._webui_urlencode(via_link))"; output_root=root, runtime=rt)
            @test dl.status == 200
            @test String(copy(dl.body)) == read(joinpath(cases, "warmup_casePST.m"), String)
            # resolving symlinks must not open the rest of the file system
            outside = SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case/download?case=$(SparlectraApp._webui_urlencode("/etc/passwd"))"; output_root=root, runtime=rt)
            @test outside.status == 303
            @test isempty(outside.body)
        end)() end

        @testset "operation-log retention" begin (function ()
            # the log is pruned by AGE at every start, and the retention is a
            # configuration key now (it used to be an environment variable only,
            # and a log written next to the runs aged without any limit at all)
            d = mktempdir()
            log = joinpath(d, "webui_operations.jsonl")
            stamp(days) = Dates.format(Dates.now(Dates.UTC) - Dates.Day(days), dateformat"yyyy-mm-ddTHH:MM:SS.sss")
            open(log, "w") do io
                println(io, "{\"timestamp\":\"$(stamp(40))Z\",\"event\":\"ancient\"}")
                println(io, "{\"timestamp\":\"$(stamp(2))Z\",\"event\":\"recent\"}")
            end
            SparlectraApp._prune_webui_operation_log!(log; SparlectraApp._webui_operation_log_options(; retention_days=10)...)
            kept = readlines(log)
            @test length(kept) == 1
            @test occursin("recent", kept[1])
            # a shorter retention drops more, which is the knob for an unwieldy log
            SparlectraApp._prune_webui_operation_log!(log; SparlectraApp._webui_operation_log_options(; retention_days=1)...)
            @test isempty(readlines(log))
            # the configuration carries it
            lowered = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true,
                overrides=Dict{String,Any}("webui" => Dict{String,Any}("operation_log_retention_days" => 3)))
            @test lowered.webui.operation_log_retention_days == 3
            @test SparlectraApp._webui_operation_log_options(; retention_days=3).retention_days == 3
            @test_throws ArgumentError Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true,
                overrides=Dict{String,Any}("webui" => Dict{String,Any}("operation_log_retention_days" => -1)))
        end)() end

        @testset "format-bound Case page options: rendered where the case reads them, saved values kept" begin (function ()
            root = mktempdir()
            cases = mktempdir()
            # a save with the case target may write machine-scope keys to the
            # configuration file: never the packaged template
            cfg_rt = joinpath(root, "rt.configuration.yaml")
            cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, cfg_rt)
            rt = (; case_directory=cases, config_file=cfg_rt, operation_log=SparlectraApp.webui_operation_log_path(root), startup_config_error=nothing, runner=SparlectraApp.start_powerflow_run)
            # a configuration that sets the frequency key away from the default
            cfg60 = joinpath(root, "conf60.yaml")
            write(cfg60, "config_version: 1\nscope: general\nmodel:\n  nominal_frequency_hz: 60.0\n")
            rt60 = merge(rt, (; config_file=cfg60))
            freq_disabled(body) = occursin(r"<select id=\"model_nominal_frequency_hz\"[^>]* disabled", body)
            selected(body) = (m = match(r"<select id=\"model_nominal_frequency_hz\"[^>]*>(.*?)</select>"s, body); m === nothing ? nothing : (s = match(r"<option value=\"([^\"]*)\" selected>", m.captures[1]); s === nothing ? nothing : s.captures[1]))
            case_page(name; runtime=rt) = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case?casefile=" * name; output_root=root, runtime).body)
            rendered(body, field) = occursin("name=\"" * field * "\"", body)
            # what the browser posts from the Case page options form: every
            # rendered input, the first value per name (the hidden false of a
            # checkbox pair), the selected option of a select that is not
            # disabled
            function posted_fields(body)
                form = match(r"<form id=\"case-options-form\".*?</form>"s, body).match
                fields = Dict{String,Any}()
                for m in eachmatch(r"<input[^>]*name=\"([^\"]+)\"[^>]*value=\"([^\"]*)\""s, form)
                    haskey(fields, m.captures[1]) || (fields[m.captures[1]] = m.captures[2])
                end
                for m in eachmatch(r"<select([^>]*)>(.*?)</select>"s, form)
                    occursin(" disabled", m.captures[1]) && continue
                    name = match(r"name=\"([^\"]+)\"", m.captures[1])
                    sel = match(r"<option value=\"([^\"]*)\" selected", m.captures[2])
                    (name === nothing || sel === nothing) || (fields[name.captures[1]] = sel.captures[1])
                end
                return fields
            end
            # a MATPOWER case: the MATPOWER options render; the system
            # frequency renders disabled (line charging is given in per
            # unit) with the configured key as its value
            write(joinpath(cases, "freq_probe.m"), "function mpc = freq_probe\nmpc.baseMVA = 100;\n")
            body = case_page("freq_probe.m"; runtime=rt60)
            @test rendered(body, "matpower_import_ratio")
            @test rendered(body, "transformer_tap_changer_model")
            @test selected(body) == "60.0"
            @test freq_disabled(body)
            # the key cannot act for a MATPOWER case: a POST that carries it
            # anyway does not store it for the case
            freq_stored(name) = haskey(Sparlectra.load_case_config(joinpath(cases, name)), "model.nominal_frequency_hz")
            SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "freq_probe.m", "settings_target" => "this_case", "model_nominal_frequency_hz" => "60.0", "matpower_import_ratio" => "reciprocal"); output_root=root, runtime=rt)
            @test Sparlectra.load_case_config(joinpath(cases, "freq_probe.m"))["matpower_import.ratio"] == "reciprocal"
            @test !freq_stored("freq_probe.m")
            # an SCF time series reads its base case (sp_case14: meta.f_nom
            # 50, every line c1 = 0): disabled at the file's 50, and a POST
            # carrying 60 does not store it
            cp(joinpath(Sparlectra.SPARLECTRA_ROOT, "data", "qsts", "sp_case14_day.qsts.zip"), joinpath(cases, "sp_case14_day.qsts.zip"))
            body_q = case_page("sp_case14_day.qsts.zip"; runtime=rt60)
            @test freq_disabled(body_q)
            @test selected(body_q) == "50.0"
            @test !haskey(posted_fields(body_q), "model_nominal_frequency_hz")
            SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14_day.qsts.zip", "settings_target" => "this_case", "model_nominal_frequency_hz" => "60.0"); output_root=root, runtime=rt)
            @test !freq_stored("sp_case14_day.qsts.zip")
            # an SCF file without meta.f_nom and without line c1: nothing
            # the frequency could change, the select is disabled and a posted
            # value is not stored
            noc = Sparlectra.scf_json_parse(read(joinpath(Sparlectra.SPARLECTRA_ROOT, "data", "scf", "sp_case14.scf.json"), String))
            delete!(noc["sparlectra"]["meta"], "f_nom")
            write(joinpath(cases, "noc1.scf.json"), Sparlectra.scf_json_string(noc))
            @test freq_disabled(case_page("noc1.scf.json"))
            SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "noc1.scf.json", "settings_target" => "this_case", "model_nominal_frequency_hz" => "60.0", "matpower_import_ratio" => "reciprocal"); output_root=root, runtime=rt)
            @test !freq_stored("noc1.scf.json")
            @test Sparlectra.load_case_config(joinpath(cases, "noc1.scf.json"))["matpower_import.ratio"] == "reciprocal"
            # an SCF file stating 60 Hz: the MATPOWER options do not render,
            # Export Solution (every format) does; the frequency select is
            # disabled at the file's 60 (the key is 50)
            scf = Sparlectra.scf_json_parse(read(joinpath(Sparlectra.SPARLECTRA_ROOT, "data", "scf", "sp_case14.scf.json"), String))
            scf["sparlectra"]["meta"]["f_nom"] = 60.0
            write(joinpath(cases, "freq60.scf.json"), Sparlectra.scf_json_string(scf))
            body60 = case_page("freq60.scf.json")
            @test !rendered(body60, "matpower_import_ratio")
            @test !rendered(body60, "matpower_import_bus_shunt_model")
            @test rendered(body60, "matpower_export_write_solution")
            @test selected(body60) == "60.0"
            @test freq_disabled(body60)
            @test occursin("data-hidden-options=\"11\"", body60)
            @test SparlectraApp.webui_form_state(selected_casefile="freq60.scf.json", case_directory=cases)["model_nominal_frequency_hz"] == 60.0
            # an SCF file without f_nom: the select is enabled and shows the
            # configured key
            pst = Sparlectra.scf_json_parse(read(joinpath(Sparlectra.SPARLECTRA_ROOT, "data", "scf", "sp_casePST.scf.json"), String))
            delete!(pst["sparlectra"]["meta"], "f_nom")
            write(joinpath(cases, "nofnom.scf.json"), Sparlectra.scf_json_string(pst))
            body_nof = case_page("nofnom.scf.json"; runtime=rt60)
            @test selected(body_nof) == "60.0"
            @test !freq_disabled(body_nof)
            # a save from the SCF case page keeps a MATPOWER value set before
            # and the key: hidden and disabled fields are not posted, and the
            # save writes posted fields only
            pre = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "freq60.scf.json", "settings_target" => "this_case", "matpower_import_ratio" => "reciprocal"); output_root=root, runtime=rt)
            @test pre.status in (302, 303)
            post_form = posted_fields(case_page("freq60.scf.json"))
            @test !haskey(post_form, "matpower_import_ratio")
            @test !haskey(post_form, "model_nominal_frequency_hz")
            @test post_form["settings_target"] == "this_case"
            saved = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", post_form; output_root=root, runtime=rt)
            @test saved.status in (302, 303)
            stored = Sparlectra.load_case_config(joinpath(cases, "freq60.scf.json"))
            @test stored["matpower_import.ratio"] == "reciprocal"
            @test !haskey(stored, "model.nominal_frequency_hz")
            @test Sparlectra.load_sparlectra_config(cfg_rt).model.nominal_frequency_hz == 50.0
            # configuration -> rendered select
            v = SparlectraApp.webui_form_state(selected_config_file=cfg60)
            @test v["model_nominal_frequency_hz"] == 60.0
            @test selected(SparlectraApp._webui_system_frequency_html(v)) == "60.0"
            # a bundled SCF case the combobox offers, not yet copied into the
            # case directory: its meta.f_nom (50) seeds the form state, not
            # the configured 60
            @test !isfile(joinpath(cases, "sp_case14.scf.json"))
            @test SparlectraApp.webui_form_state(selected_casefile="sp_case14.scf.json", selected_config_file=cfg60, case_directory=cases)["model_nominal_frequency_hz"] == 50.0
            # form -> override -> configuration
            req = SparlectraApp.powerflow_webui_request(Dict("casefile" => "freq_probe.m", "model_nominal_frequency_hz" => "60.0"); case_directory=cases, operation_log=root)
            @test req["config_overrides"]["model.nominal_frequency_hz"] == 60.0
            cfg, _ = Sparlectra._load_api_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, Sparlectra.validate_gui_config_overrides(req["config_overrides"]))
            @test cfg.model.nominal_frequency_hz == 60.0
            @test "model.nominal_frequency_hz" in cfg.user_set_keys
            # saved, then run (the run form posts no frequency): the
            # Case page save above stored no frequency, so no frequency
            # warning in run.log. A posted 50 is not saved against the
            # file's 60; a 50 stored before the lock still reaches the run
            # and gives exactly one warning
            freq_warnings(run) = count(l -> occursin("meta.f_nom = 60.0 Hz", l), eachline(joinpath(String(run["output_dir"]), "run.log")))
            freq_line(run) = only(filter(l -> startswith(l, "System frequency: "), collect(eachline(joinpath(String(run["output_dir"]), "run.log")))))
            for (posted, expected) in ((nothing, 0), ("50.0", 1))
                if posted !== nothing
                    save = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "freq60.scf.json", "settings_target" => "this_case", "model_nominal_frequency_hz" => posted); output_root=root, runtime=rt)
                    @test save.status in (302, 303)
                    @test !freq_stored("freq60.scf.json")
                    Sparlectra.write_case_config(joinpath(cases, "freq60.scf.json"), merge(Sparlectra.load_case_config(joinpath(cases, "freq60.scf.json")), Dict{String,Any}("model.nominal_frequency_hz" => parse(Float64, posted))))
                end
                patterns = expected == 0 ? () : (r"meta\.f_nom = 60\.0 Hz.*model\.nominal_frequency_hz = 50\.0 Hz",)
                run = run_with_expected_warnings(patterns) do
                    SparlectraApp.start_powerflow_run(Dict("casefile" => "freq60.scf.json", "config_file" => cfg_rt, "output_root" => joinpath(root, "runs_" * something(posted, "none"))); case_directory=cases)
                end
                @test run["status"] == "succeeded"
                @test freq_warnings(run) == expected
                # result.json names the frequency the network used: the
                # file's 60 also when the saved key says 50
                res = Sparlectra.scf_json_parse(read(joinpath(String(run["output_dir"]), "result.json"), String))
                @test res["metadata"]["network_nominal_frequency_hz"] == 60.0
                # sp_case14 carries c1 = 0 on every line: the frequency never
                # entered the model, and run.log says so instead of 60 Hz
                @test res["metadata"]["frequency_source"] == "not_used"
                @test freq_line(run) == "System frequency: not used (no capacitance in farad); case file states 60 Hz"
            end
            # the next save of the case options removes the stored 50
            @test SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", posted_fields(case_page("freq60.scf.json")); output_root=root, runtime=rt).status in (302, 303)
            @test !freq_stored("freq60.scf.json")
            # a case with line capacitance in farad (two_islands_prio, c1 > 0):
            # the stated 50 Hz is used; without meta.f_nom the configured
            # frequency is, and 60 Hz carries 1.2 times the line charging
            cap = Sparlectra.scf_json_parse(read(joinpath(Sparlectra.SPARLECTRA_ROOT, "data", "scf", "two_islands_prio.scf.json"), String))
            write(joinpath(cases, "cap_fnom.scf.json"), Sparlectra.scf_json_string(cap))
            delete!(cap["sparlectra"]["meta"], "f_nom")
            write(joinpath(cases, "cap_nof.scf.json"), Sparlectra.scf_json_string(cap))
            balance(run) = only(filter(l -> occursin("total network power balance", l), collect(eachline(joinpath(String(run["output_dir"]), "run.log")))))
            cap_runs = Dict{String,Any}()
            for (case, cfg, source, line) in (("cap_fnom.scf.json", cfg_rt, "case_file", "50 Hz (case file)"),
                                              ("cap_nof.scf.json", cfg_rt, "configuration", "50 Hz (configuration)"),
                                              ("cap_nof.scf.json", cfg60, "configuration", "60 Hz (configuration)"))
                run = SparlectraApp.start_powerflow_run(Dict("casefile" => case, "config_file" => cfg, "output_root" => joinpath(root, "cap_runs")); case_directory=cases)
                @test run["status"] == "succeeded"
                md = Sparlectra.scf_json_parse(read(joinpath(String(run["output_dir"]), "result.json"), String))["metadata"]
                @test md["frequency_source"] == source
                @test md["frequency_line"] == line
                @test freq_line(run) == "System frequency: " * line
                cap_runs[line] = run
            end
            # same file, same 50 Hz: the same solution; 60 Hz changes it
            @test balance(cap_runs["50 Hz (case file)"]) == balance(cap_runs["50 Hz (configuration)"])
            @test balance(cap_runs["60 Hz (configuration)"]) != balance(cap_runs["50 Hz (configuration)"])
            # cgmes_import.read_dynamics per case (default off): a CGMES
            # delivery with a DY profile (the SmallGrid fixture plus one DY
            # file with two machine models). The switch renders on the CGMES
            # Case page only; checked, the save writes true into the case
            # file and the run reads the DY profile; unchecked again, the
            # hidden false equals the default, the key leaves the case file
            # and the run skips the profile.
            @test !rendered(body60, "cgmes_read_dynamics")
            sg_zips = sort(filter(f -> endswith(f, ".zip"), readdir(joinpath(@__DIR__, "fixtures", "cgmes_smallgrid"); join = true)))
            sg_store = Sparlectra.CGMESImporter.loadCGMES(sg_zips)
            sm_mrid(name) = only(o.mrid for o in Sparlectra.CGMESImporter.objectsOf(sg_store, :SynchronousMachine) if Sparlectra.CGMESImporter.str(o, :name) == name)
            dy_model(id, mrid) = """
              <cim:SynchronousMachineTimeConstantReactance rdf:ID="$(id)">
                <cim:DynamicsFunctionBlock.enabled>true</cim:DynamicsFunctionBlock.enabled>
                <cim:RotatingMachineDynamics.inertia>4.25</cim:RotatingMachineDynamics.inertia>
                <cim:SynchronousMachineDynamics.SynchronousMachine rdf:resource="#$(mrid)"/>
              </cim:SynchronousMachineTimeConstantReactance>
            """
            dy_xml = """<?xml version="1.0" encoding="UTF-8"?>
            <rdf:RDF xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
              <md:FullModel rdf:about="urn:uuid:5f0e4c1a-0000-4000-8000-000000000507">
                <md:Model.profile>http://entsoe.eu/CIM/Dynamics/3/1</md:Model.profile>
              </md:FullModel>
            $(dy_model("_wdy_m1", sm_mrid("Claytor SM")))$(dy_model("_wdy_m2", sm_mrid("Torrey SM")))</rdf:RDF>
            """
            open(joinpath(cases, "sg_dy.zip"), "w") do io
                Sparlectra.CGMESImporter.ZipArchives.ZipWriter(io) do w
                    for z in sg_zips
                        Sparlectra.CGMESImporter.ZipArchives.zip_newfile(w, basename(z))
                        write(w, read(z))
                    end
                    Sparlectra.CGMESImporter.ZipArchives.zip_newfile(w, "sg_DY.xml")
                    write(w, dy_xml)
                end
            end
            body_cg = case_page("sg_dy.zip")
            @test rendered(body_cg, "cgmes_read_dynamics")
            post_cg = posted_fields(body_cg)
            @test post_cg["cgmes_read_dynamics"] == "false"
            # the DY section of the import report, up to the next heading;
            # the report escapes Markdown characters (cgmes\_import), so
            # the backslashes go before the comparison
            function dy_section(run)
                t = read(joinpath(String(run["output_dir"]), "cgmes_import_report.md"), String)
                i = findfirst("## Dynamic data (DY profile)\n", t)
                i === nothing && return ""
                j = findnext("\n## ", t, last(i))
                return replace(t[last(i)+1:(j === nothing ? lastindex(t) : first(j))], r"\\(.)" => s"\1")
            end
            for (checked, stored, section) in ((true, true, "DY profile read (cgmes_import.read_dynamics): 2 unit(s) with dynamic data"), (false, nothing, "DY profile present (1 file(s)), not read: cgmes_import.read_dynamics = false"))
                # the browser posts the hidden false first and, when checked,
                # the true of the checkbox; the save parser keeps the last
                save_cg = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", merge(post_cg, Dict{String,Any}("cgmes_read_dynamics" => checked ? "true" : "false")); output_root=root, runtime=rt)
                @test save_cg.status in (302, 303)
                @test get(Sparlectra.load_case_config(joinpath(cases, "sg_dy.zip")), "cgmes_import.read_dynamics", nothing) === stored
                run_cg = SparlectraApp.start_powerflow_run(Dict("casefile" => "sg_dy.zip", "config_file" => cfg_rt, "output_root" => joinpath(root, "runs_dy_" * string(checked))); case_directory=cases)
                @test run_cg["status"] == "succeeded"
                @test occursin(section, dy_section(run_cg))
            end
            # cgmes_import.equivalent_injection_as_source per case (default
            # off), same pattern: checked, the case file stores true and the
            # run imports the boundary EquivalentInjections of SmallGrid as
            # sources (one importer notice counts them); unchecked, the key
            # leaves the case file and no injection becomes a source
            @test !rendered(body60, "cgmes_equivalent_injection_as_source")
            @test rendered(body_cg, "cgmes_equivalent_injection_as_source")
            @test post_cg["cgmes_equivalent_injection_as_source"] == "false"
            for (checked, stored) in ((true, true), (false, nothing))
                save_ei = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", merge(post_cg, Dict{String,Any}("cgmes_equivalent_injection_as_source" => checked ? "true" : "false")); output_root=root, runtime=rt)
                @test save_ei.status in (302, 303)
                @test get(Sparlectra.load_case_config(joinpath(cases, "sg_dy.zip")), "cgmes_import.equivalent_injection_as_source", nothing) === stored
                run_ei = SparlectraApp.start_powerflow_run(Dict("casefile" => "sg_dy.zip", "config_file" => cfg_rt, "output_root" => joinpath(root, "runs_ei_" * string(checked))); case_directory=cases)
                @test run_ei["status"] == "succeeded"
                report = replace(read(joinpath(String(run_ei["output_dir"]), "cgmes_import_report.md"), String), r"\\(.)" => s"\1")
                @test occursin(r"notice: \d+ EquivalentInjection\(s\) imported as source", report) == checked
            end
        end)() end

        @testset "browser opening falls back to the system default" begin (function ()
            # Reported from Windows 11 with Edge uninstalled and only
            # Firefox present: nothing happened on start, the user had to type
            # 127.0.0.1:8080 by hand. Cause: the app-window attempt knows only the
            # Chromium family (they alone support the chromeless --app= window),
            # and the generic fallback returned early unless the platform was
            # Linux. Windows and macOS therefore had NO fallback at all.
            url = "http://127.0.0.1:8080/powerflow"
            no_exe(_) = nothing
            no_path(_) = false
            no_env = Dict{String,String}()

            # the reported case: Windows, not one Chromium browser anywhere
            win = SparlectraApp._webui_browser_open_command(url; platform=:windows,
                executable_lookup=no_exe, path_exists=no_path, environment=no_env)
            @test win !== nothing
            @test win[2] == :windows_start
            # cmd builtin, and the empty "" is the window TITLE: without it cmd
            # takes the quoted URL as a title and opens a console window
            @test occursin("cmd", string(win[1]))
            @test occursin(url, string(win[1]))

            # macOS keeps `open` as its fallback
            mac_open(name) = name == "open" ? "/usr/bin/open" : nothing
            mac = SparlectraApp._webui_browser_open_command(url; platform=:macos,
                executable_lookup=mac_open, path_exists=no_path, environment=no_env)
            @test mac !== nothing
            @test mac[2] == :macos_open

            # Linux is unchanged
            lin_xdg(name) = name == "xdg-open" ? "/usr/bin/xdg-open" : nothing
            lin = SparlectraApp._webui_browser_open_command(url; platform=:linux,
                executable_lookup=lin_xdg, path_exists=no_path, environment=no_env)
            @test lin !== nothing
            @test lin[2] == :xdg_open

            # and the precedence still holds: where a Chromium browser exists it
            # wins, because the app window is the nicer result
            chrome(name) = name == "chrome.exe" ? "C:\\chrome.exe" : nothing
            win_chrome = SparlectraApp._webui_browser_open_command(url; platform=:windows,
                executable_lookup=chrome, path_exists=no_path, environment=no_env)
            @test win_chrome[2] == :app_window
        end)() end

        @testset "app window geometry: configured size, saved by Stop Web UI" begin (function ()
            # The app window opens with webui.window_width/height and, once a
            # Stop Web UI click saved one, at webui.window_x/y. Stop sends the
            # geometry as hidden fields; a minimized window (the hidden flag,
            # or the -32000 Chromium parks it at on Windows) and implausible
            # values are not saved.
            url = "http://127.0.0.1:8080/powerflow"
            chrome(name) = name == "google-chrome" ? "/usr/bin/google-chrome" : nothing
            no_path(_) = false
            no_env = Dict{String,String}()
            default_cmd = SparlectraApp._webui_app_window_command(url; geometry=SparlectraApp._WEBUI_DEFAULT_WINDOW_GEOMETRY, profile_dir="", platform=:linux,
                executable_lookup=chrome, path_exists=no_path, environment=no_env)
            @test "--window-size=1500,950" in default_cmd.exec
            @test !any(startswith("--window-position="), default_cmd.exec)
            @test !any(startswith("--user-data-dir="), default_cmd.exec)
            # the flags act only in a new browser process (own profile), and the
            # position only through X11 under Wayland
            placed = SparlectraApp._webui_app_window_command(url; geometry=(width=1200, height=800, x=-10, y=20), profile_dir="/state/browser", platform=:linux,
                executable_lookup=chrome, path_exists=no_path, environment=Dict("WAYLAND_DISPLAY" => "wayland-0"))
            @test "--window-size=1200,800" in placed.exec
            @test "--window-position=-10,20" in placed.exec
            @test "--user-data-dir=/state/browser" in placed.exec
            @test "--ozone-platform=x11" in placed.exec
            @test !("--ozone-platform=x11" in default_cmd.exec)
            mac_cmd = SparlectraApp._webui_app_window_command(url; geometry=(width=1200, height=800, x=5, y=6), profile_dir="/state/browser", platform=:macos,
                executable_lookup=no_path, path_exists=(p -> p == "/Applications/Google Chrome.app"), environment=Dict("WAYLAND_DISPLAY" => "wayland-0"))
            @test "--window-position=5,6" in mac_cmd.exec
            @test "--user-data-dir=/state/browser" in mac_cmd.exec
            @test !("--ozone-platform=x11" in mac_cmd.exec)

            # a saved position on a monitor that is gone is dropped, so the
            # window opens on the main monitor; unknown layout keeps it
            screens = SparlectraApp._webui_parse_powershell_screens("0 0 1920 1080\r\n1920 0 2560 1440\r\n")
            @test screens == [(0, 0, 1920, 1080), (1920, 0, 2560, 1440)]
            @test SparlectraApp._webui_parse_xrandr_monitors("Monitors: 2\n 0: +*DP-1 2560/597x1440/336+0+0  DP-1\n 1: +HDMI-1 1920/510x1080/290+-1920+0  HDMI-1\n") == [(0, 0, 2560, 1440), (-1920, 0, 1920, 1080)]
            on_second = (width=1200, height=800, x=2000, y=100)
            @test SparlectraApp._webui_visible_geometry(on_second, screens) == on_second
            @test (@test_logs (:info, r"on no connected monitor") SparlectraApp._webui_visible_geometry(on_second, screens[1:1])) == (width=1200, height=800, x=nothing, y=nothing)
            @test SparlectraApp._webui_visible_geometry(on_second, nothing) == on_second
            @test SparlectraApp._webui_monitor_rects(:linux; executable_lookup=(_ -> nothing)) === nothing

            # the template carries the defaults, the loader validates
            cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true)
            @test (cfg.webui.window_width, cfg.webui.window_height, cfg.webui.window_x, cfg.webui.window_y) == (1500, 950, nothing, nothing)
            @test_throws ArgumentError Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true,
                overrides=Dict{String,Any}("webui" => Dict{String,Any}("window_width" => 0)))

            mktempdir() do dir
                file = joinpath(dir, "configuration.yaml")
                cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, file)
                @test SparlectraApp._webui_window_geometry(file) == SparlectraApp._WEBUI_DEFAULT_WINDOW_GEOMETRY
                form(; kw...) = Dict{String,String}("window_width" => "1234", "window_height" => "876", "window_x" => "12", "window_y" => "-34", "window_hidden" => "0",
                    (String(k) => v for (k, v) in kw)...)
                saved = SparlectraApp._webui_save_window_geometry!(file, form())
                @test saved.status == "saved"
                @test SparlectraApp._webui_window_geometry(file) == (width=1234, height=876, x=12, y=-34)
                @test Sparlectra.load_sparlectra_config(file; reload=true).webui.window_x == 12
                @test issubset(["webui.window_width", "webui.window_height", "webui.window_x", "webui.window_y"], SparlectraApp._webui_user_keys(file))
                # a minimized window, a parked one, nonsense sizes and a POST
                # without the fields (curl, an older page) keep the file
                for (fields, reason) in (
                        ((window_hidden="1",), "window_minimized"),
                        ((window_x="-32000", window_y="-32000"), "position_out_of_range"),
                        ((window_width="0",), "size_out_of_range"),
                        ((window_height="abc",), "not_a_number"),
                        ((window_width="",), "no_geometry"),
                    )
                    skipped = SparlectraApp._webui_save_window_geometry!(file, form(; fields...))
                    @test skipped.status == "skipped"
                    @test skipped.reason == reason
                end
                @test SparlectraApp._webui_save_window_geometry!(file, Dict{String,String}()).reason == "no_geometry"
                @test SparlectraApp._webui_window_geometry(file) == (width=1234, height=876, x=12, y=-34)
                # the packaged template is never written
                @test SparlectraApp._webui_save_window_geometry!(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, form()).reason == "packaged_template"
            end
            # the Stop form carries the fields the handler reads
            page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow"; output_root=mktempdir()).body)
            for name in ("window_width", "window_height", "window_x", "window_y", "window_hidden")
                @test occursin("name=\"$(name)\"", page)
            end
        end)() end

        @testset "sysimage launcher decision" begin (function ()
            # The launcher decides whether the Web UI starts from the image or
            # compiles on first use. It runs on plain Base BEFORE the package is
            # loaded (tools/sysimage_launcher.jl), so it is tested through the
            # module, not through Sparlectra. The last case is the expensive one:
            # the metadata pins the Julia version and the Manifest, i.e. the
            # DEPENDENCIES, so an image built before a src/ edit still looked
            # fresh and the Web UI silently served old code.
            launcher = Module(:LauncherUnderTest)
            Base.include(launcher, joinpath(Sparlectra.SPARLECTRA_ROOT, "tools", "sysimage_launcher.jl"))
            # Base.invokelatest also for the binding: Julia 1.12+ warns (and a later
            # version errors) on reading a global defined in this same top-level
            # expression by the include above
            SL = Base.invokelatest(getfield, launcher, :SysimageLauncher)
            # Julia 1.12: a module included at run time defines its methods in a
            # NEWER world than this function, so a direct call raises "method too
            # new". Binding lookup and call both go through invokelatest.
            sl(name, args...) = Base.invokelatest(Base.invokelatest(getfield, SL, name), args...)
            slval(name) = Base.invokelatest(getfield, SL, name)
            mktempdir() do tmp
                proj = joinpath(tmp, "proj")
                mkpath(joinpath(proj, "src"))
                manifest = joinpath(proj, "Manifest.toml")
                write(manifest, "# manifest fixture\n")
                write(joinpath(proj, "src", "Fixture.jl"), "module Fixture end\n")
                imgdir = joinpath(tmp, "image")
                mkpath(imgdir)
                img = joinpath(imgdir, "sparlectra.so")
                meta = joinpath(imgdir, "sysimage_meta.toml")
                sha = bytes2hex(open(SHA.sha256, manifest))
                full_meta = string("julia_version = \"", VERSION, "\"\nmanifest_sha256 = \"", sha, "\"\n")

                # SparlectraApp.webui_sysimage_problem answers the same question for
                # the sysimage page; the launcher cannot call it (it runs before the
                # package is loaded), so every state below checks that both agree,
                # or the page would one day call an image fine that the launcher
                # refuses to start from
                verdict() = (v = sl(:sysimage_problem, img, proj); @test v == SparlectraApp.webui_sysimage_problem(image_path=img, project_dir=proj); v)
                @test verdict() == "no sysimage found"
                write(img, "not a real image")
                # image without metadata is unusable: the contract is the pair
                @test verdict() == "no sysimage found"
                write(meta, "julia_version = \"9.9.9\"\nmanifest_sha256 = \"$(sha)\"\n")
                @test occursin("built for Julia 9.9.9", verdict())
                write(meta, string("julia_version = \"", VERSION, "\"\nmanifest_sha256 = \"deadbeef\"\n"))
                @test verdict() == "the sysimage does not match the current Manifest.toml"
                write(meta, full_meta)
                touch(img)
                @test verdict() === nothing
                # a source file newer than the image disables it (mtime resolution
                # can be one second, hence the wait)
                sleep(1.1)
                touch(joinpath(proj, "src", "Fixture.jl"))
                @test occursin("older than", verdict())
                write(meta, "kaputt = [[[")
                @test verdict() == "the sysimage metadata is unreadable"

                # the owner key the package writes is the one the launcher reads:
                # the image belongs to proj, whatever the trailing separator
                written = SparlectraApp.write_sysimage_meta(meta; manifest_path=manifest, project_dir=proj * "/")
                @test written["project_dir"] == rstrip(realpath(proj), ['/', '\\'])
                @test sl(:foreign_sysimage_owner, img, proj) === nothing
                @test sl(:foreign_sysimage_owner, img, tmp) == written["project_dir"]

                # the decision itself, through its seams: a fixed answer, a recorded
                # relaunch, the image this process runs on
                slk(name, args...; kw...) = Base.invokelatest(Base.invokelatest(getfield, SL, name), args...; kw...)
                decide(args; ask, current = "") = begin
                    relaunches = String[]
                    prompts = String[]
                    # stdout can be redirected into a file, not into a buffer
                    capture = joinpath(tmp, "decision.txt")
                    open(capture, "w") do out
                        redirect_stdout(out) do
                            slk(:handle_sysimage, args, "start_webui.jl", proj; image = img, ask = (p, _...) -> (push!(prompts, String(p)); ask), relaunch = (image, _...) -> push!(relaunches, String(image)), current_image = current)
                        end
                    end
                    (text = read(capture, String), relaunches = relaunches, prompts = prompts)
                end
                # an outdated image and the answer no: image and metadata are gone,
                # the start goes on without an image
                write(img, "not a real image")
                write(meta, string("julia_version = \"", VERSION, "\"\nmanifest_sha256 = \"deadbeef\"\n"))
                r = decide(String[]; ask = false)
                @test occursin("Sysimage is out of date and was removed.", r.text)
                @test occursin("Starting without a sysimage", r.text)
                @test isempty(r.relaunches)
                @test !isfile(img) && !isfile(meta)
                # a current image is left alone and used
                write(img, "not a real image")
                write(meta, full_meta)
                touch(img)
                r = decide(String[]; ask = false)
                @test r.relaunches == [img]
                @test isfile(img) && isfile(meta)
                # an outdated image while the process runs on an image given from
                # outside: a warning, nothing removed, nothing built
                write(meta, string("julia_version = \"", VERSION, "\"\nmanifest_sha256 = \"deadbeef\"\n"))
                r = decide(String[]; ask = false, current = joinpath(tmp, "someone_elses.so"))
                @test occursin("given from outside", r.text)
                @test isempty(r.relaunches)
                @test isfile(img) && isfile(meta)
                # one image per user: the image of another checkout (a release
                # clone next to the dev checkout) is neither used nor removed,
                # even when it would pass every check for this one; the build
                # question names the image it replaces. full_meta above has no
                # project_dir (an image built before the key) and keeps the old
                # behaviour, as the cases above show.
                other = joinpath(tmp, "other_checkout", "app")
                mkpath(other)
                write(meta, string(full_meta, "project_dir = \"", escape_string(other), "\"\n"))
                touch(img)
                r = decide(String[]; ask = false)
                @test occursin("belongs to $(other)", r.text)
                @test isempty(r.relaunches)
                @test isfile(img) && isfile(meta)
                @test length(r.prompts) == 1 && occursin("replaces the sysimage of $(other)", r.prompts[1])
                # same with an outdated Manifest: still not this checkout's image
                write(meta, string("julia_version = \"", VERSION, "\"\nmanifest_sha256 = \"deadbeef\"\nproject_dir = \"", escape_string(other), "\"\n"))
                r = decide(String[]; ask = false)
                @test !occursin("was removed", r.text)
                @test isfile(img) && isfile(meta)
                # this checkout recorded (with a trailing separator, as
                # SPARLECTRA_APP_ROOT writes it before normalization): used when
                # current, removed when outdated, as before
                write(meta, string(full_meta, "project_dir = \"", escape_string(proj * "/"), "\"\n"))
                touch(img)
                r = decide(String[]; ask = false)
                @test r.relaunches == [img]
                write(meta, string("julia_version = \"", VERSION, "\"\nmanifest_sha256 = \"deadbeef\"\nproject_dir = \"", escape_string(proj), "\"\n"))
                r = decide(String[]; ask = false)
                @test occursin("Sysimage is out of date and was removed.", r.text)
                @test !isfile(img) && !isfile(meta)
                @test slk(:external_sysimage, img; current = "") === nothing
                @test slk(:external_sysimage, img; current = joinpath(Sys.BINDIR, "..", "lib", "julia", "sys.so")) === nothing
                @test slk(:remove_stale_sysimage, joinpath(tmp, "nowhere.so")) === :absent
            end
            # no terminal in the test process, so the question answers itself with
            # no: an unattended start goes without an image, a build is opt-in
            @test sl(:ask_build, "never shown") == false
            # a reader left pending by a question that timed out is reused by
            # the next question on the same stream: a second readline would
            # race with it and the older task would take the answer
            console = Base.BufferStream()
            ask_console(seconds) = Base.invokelatest(Base.invokelatest(getfield, SL, :console_answer), "?"; io = devnull, input = console, seconds)
            @test ask_console(0.2) == ""
            first_reader = slval(:_PENDING_CONSOLE_READER)[][2]
            @test !istaskdone(first_reader)
            write(console, "y\n")
            @test ask_console(5.0) == "y"
            @test istaskdone(first_reader) && fetch(first_reader) == "y"
            @test slval(:_PENDING_CONSOLE_READER)[] === nothing
            # the end of the input is an empty answer
            @test Base.invokelatest(Base.invokelatest(getfield, SL, :console_answer), "?"; io = devnull, input = IOBuffer(""), seconds = 5.0) == ""
            @test slval(:REBUILD_FLAG) == "--rebuild-sysimage"
            @test slval(:NO_IMAGE_FLAG) == "--no-sysimage"
            # platform path contract, still owned by the package for everything
            # that runs after the load
            @test endswith(sl(:sysimage_path), SparlectraApp.webui_sysimage_ext())
            @test sl(:sysimage_path) == SparlectraApp.webui_sysimage_path()
        end)() end

        # Two runs of the same case under different settings are the normal way to
        # look at a Q-limit mode or a solver choice. The history could list them but
        # not put them side by side, so the comparison had to happen in the head or
        # in two browser tabs.
          # The box should appear only where the two runs
          # are comparable, and NOT decided by network size, because the same
          # case with an external-grid source and with a slack is a pair one
          # wants side by side. The criterion is the run kind: the page reads
          # power-flow artifacts, so only power-flow kinds qualify.
          @testset "only finished power-flow runs offer the box" begin (function ()
              @test SparlectraApp._webui_comparable_kind("")
              @test SparlectraApp._webui_comparable_kind("powerflow")
              @test SparlectraApp._webui_comparable_kind("diagnose")
              @test SparlectraApp._webui_comparable_kind("powerflow_se_start")
              for kind in ("se", "short_circuit", "contingency", "import_analysis")
                  @test !SparlectraApp._webui_comparable_kind(kind)
              end

              pf = Dict{String,Any}("run_id" => "pf-1", "available" => true, "status" => "succeeded", "run_mode" => "",
                  "timestamp" => "2026-09-09 01:00:00", "casefile" => "/very/long/path/to/cases/case118.m", "config_file" => "/etc/sparlectra/configuration.yaml")
              se = merge(pf, Dict{String,Any}("run_id" => "se-1", "run_mode" => "se"))
              running = merge(pf, Dict{String,Any}("run_id" => "pf-2", "status" => "running"))
              gone = merge(pf, Dict{String,Any}("run_id" => "pf-3", "available" => false))
              @test SparlectraApp._webui_run_comparable(pf)
              @test !SparlectraApp._webui_run_comparable(se)
              @test !SparlectraApp._webui_run_comparable(running)
              @test !SparlectraApp._webui_run_comparable(gone)

              html = SparlectraApp.render_powerflow_history([pf, se], mktempdir())
              @test occursin("name=\"run\" value=\"pf-1\"", html)
              @test !occursin("name=\"run\" value=\"se-1\"", html)
          end)() end

          # The result page names the case and the phase up top, where a reader
          # looks first.
          @testset "result page: case, phase and controllers up top" begin (function ()
              # the value of a summary card, nothing when the card is absent
              card(html, label) = (m = match(Regex("summary-label\">" * label * "</span><code[^>]*>([^<]*)</code>"), html); m === nothing ? nothing : String(m.captures[1]))
              probe = Dict{String,Any}("run_id" => "r-1", "status" => "succeeded", "success" => true, "converged" => true,
                  "casefile" => "/some/where/sp_case14.scf.json", "resolved_casefile" => "/some/where/sp_case14.scf.json",
                  "current_phase" => "finished", "final_outcome" => Dict{String,Any}("solver" => "rectangular"), "artifacts" => Any[])
              html = SparlectraApp.render_powerflow_result(probe)
              @test card(html, "Case") == "sp_case14.scf.json"
              @test card(html, "Phase") == "finished"
              # an active run carries the same two cards
              active_html = SparlectraApp.render_powerflow_result(merge(probe, Dict{String,Any}("status" => "running", "current_phase" => "linear_solve")))
              @test card(active_html, "Case") == "sp_case14.scf.json"
              @test card(active_html, "Phase") == "linear_solve"
              # A live snapshot carries `nothing` where the run has not produced
              # the value yet; the card once printed that word for the whole run
              # and named the case only at the end.
              live = Dict{String,Any}("run_id" => "r-2", "status" => "running", "success" => false,
                  "casefile" => "/some/where/case14.m", "resolved_casefile" => nothing, "current_phase" => nothing, "last_phase" => nothing, "artifacts" => Any[])
              live_html = SparlectraApp.render_powerflow_result(live)
              @test card(live_html, "Case") == "case14.m"
              @test card(live_html, "Phase") == "n/a"
              @test !occursin(">nothing<", live_html)

              # the controller counts from the run metadata reach the page (a
              # Q(U) machine read from a case file is named there, not only in
              # the Control column of the result print); a run without
              # controllers has no card
              @test card(html, "Controllers") === nothing
              with_ctrl = merge(probe, Dict{String,Any}("metadata" => Dict{String,Any}("controllers" => Dict{String,Any}("tap" => 1, "qu" => 1, "pu" => 0))))
              @test card(SparlectraApp.render_powerflow_result(with_ctrl), "Controllers") == "Q(U) 1 · tap 1"
              @test SparlectraApp._webui_control_summary(Dict{String,Any}("metadata" => Dict{String,Any}("controllers" => Dict{String,Any}("tap" => 0, "qu" => 0, "pu" => 0)))) === nothing
          end)() end

          @testset "the route insists on exactly two runs" begin (function ()
              root = mktempdir()
              for target in ("/powerflow/compare", "/powerflow/compare?run=only-one")
                  response = SparlectraApp.route_sparlectra_webui("GET", target; output_root=root)
                  @test (target, response.status) == (target, 400)
              end
              missing_response = SparlectraApp.route_sparlectra_webui("GET", "/powerflow/compare?run=nope-a&run=nope-b"; output_root=root)
              @test missing_response.status == 404
          end)() end

          @testset "repeated query keys survive the parser" begin (function ()
              # `_webui_parse_pairs` returns a Dict and keeps only the last value;
              # a set of checkboxes with one name needs all of them
              @test SparlectraApp._webui_query_values("/powerflow/compare?run=a&run=b", "run") == ["a", "b"]
              @test SparlectraApp._webui_query_values("/powerflow/compare?other=x", "run") == String[]
              @test SparlectraApp._webui_query_values("/powerflow/compare", "run") == String[]
          end)() end

          @testset "the page reads what the runs wrote" begin (function ()
              dir_a = mktempdir()
              dir_b = mktempdir()
              write(joinpath(dir_a, "effective_config.yaml"), "runtime:\n  case_name: sp_case_a\n  casefile: /state/data/sp_case_a.scf.json\npower_flow:\n  tol: 1.0e-8\n  qlimits:\n    enforcement_mode: active_set\n")
              write(joinpath(dir_b, "effective_config.yaml"), "runtime:\n  case_name: sp_case_b\n  casefile: /state/data/sp_case_b.scf.json\npower_flow:\n  tol: 1.0e-8\n  qlimits:\n    enforcement_mode: classic_simultaneous\n")
              diff = SparlectraApp._webui_compare_config_diff(dir_a, dir_b)
              # the nested key keeps its real path; a "remember the last section"
              # reader concatenated every section it had ever seen
              @test diff == [("power_flow.qlimits.enforcement_mode", "active_set", "classic_simultaneous"), ("runtime.case_name", "sp_case_a", "sp_case_b"), ("runtime.casefile", "/state/data/sp_case_a.scf.json", "/state/data/sp_case_b.scf.json")]

              write(joinpath(dir_a, "q_limit_events.csv"), "iteration,bus,side\n2,19,min\n2,32,min\n")
              write(joinpath(dir_b, "q_limit_events.csv"), "iteration,bus,side\n1,19,min\n")
              @test SparlectraApp._webui_compare_qlimit_buses(dir_a) == [19, 32]
              @test SparlectraApp._webui_compare_qlimit_buses(dir_b) == [19]
              @test SparlectraApp._webui_compare_qlimit_buses(mktempdir()) === nothing

              # The detailed export writes one of three formats, and A carries the
              # German one: ';' delimiter, decimal comma. A comma-splitting reader
              # found no columns there and the page then claimed the two runs had
              # nothing in common ("The two runs share no bus names", as
              # reported). The fixtures below are the real header of
              # bus_voltages_complex.csv, both formats, so the parser is judged on
              # what the runs actually write.
              head_de = "bus;bus_name;type;vm_pu;va_deg;vn_kV;q_limit_hit;original_bus_name"
              head_us = "bus,bus_name,type,vm_pu,va_deg,vn_kV,q_limit_hit,original_bus_name"
              write(joinpath(dir_a, "bus_voltages_complex.csv"),
                  "$(head_de)\n1;11001;PQ;1;0;138;false;NEWBERRY 1\n2;11002;PQ;0,99;-1,5;69;false;NEWBERRY 2\n")
              write(joinpath(dir_b, "bus_voltages_complex.csv"),
                  "$(head_us)\n1,11001,PQ,1.0,0.0,138,false,NEWBERRY 1\n2,11002,PQ,0.98,-1.0,69,false,NEWBERRY 2\n")
              va = SparlectraApp._webui_compare_voltages(dir_a)
              vb = SparlectraApp._webui_compare_voltages(dir_b)
              @test va["2"].vm == 0.99
              @test va["2"].va == -1.5           # decimal comma, not a second field
              @test va["2"].name == "NEWBERRY 2"  # the name a reader recognizes
              @test vb["2"].vm == 0.98
              @test SparlectraApp._webui_compare_voltages(mktempdir()) === nothing

              # grouped thousands, and a quoted cell because the group separator IS
              # the delimiter in the excel_us format
              @test SparlectraApp._webui_compare_split_csv("a,\"1,234.5\",b", ',') == ["a", "1,234.5", "b"]
              @test SparlectraApp._webui_compare_number("1,234.5", '.', ',') == 1234.5
              @test SparlectraApp._webui_compare_number("1.234,5", ',', '.') == 1234.5

              # losses come from the branch table the run wrote, summed
              write(joinpath(dir_a, "branch_flows.csv"),
                  "branch;from_bus;to_bus;p_loss_MW;q_loss_MVar\nB1;1;2;0,5;1,25\nB2;2;3;0,25;0,75\n")
              write(joinpath(dir_b, "branch_flows.csv"),
                  "branch,from_bus,to_bus,p_loss_MW,q_loss_MVar\nB1,1,2,0.4,1.25\nB2,2,3,0.25,0.75\n")
              la = SparlectraApp._webui_compare_losses(dir_a)
              @test la.p_MW ≈ 0.75
              @test la.q_MVAr ≈ 2.0
              @test la.branches == 2
              @test SparlectraApp._webui_compare_losses(mktempdir()) === nothing

              a = Dict{String,Any}("run_id" => "run-a", "output_dir" => dir_a, "casefile" => "case118.m",
                  "status" => "succeeded", "converged" => true, "iterations" => 6, "final_mismatch" => 1.0e-13)
              b = Dict{String,Any}("run_id" => "run-b", "output_dir" => dir_b, "casefile" => "case118.m",
                  "status" => "succeeded", "converged" => true, "iterations" => 23, "final_mismatch" => 7.4e-12)
              html = SparlectraApp.render_powerflow_compare(a, b)
              # the case of each run stands in the table on top; its path is no
              # configuration difference of its own
              @test !occursin("/state/data/sp_case_a.scf.json", html)
              # a path value shows its file name, and the parent directory
              # as well where both runs name the same file in different places
              @test SparlectraApp._webui_compare_value_cell("/x/data/case9.m", "/x/data/case14.m") == "<span title=\"/x/data/case9.m\">case9.m</span>"
              @test SparlectraApp._webui_compare_value_cell("/x/runA/configuration.yaml", "/x/runB/configuration.yaml") == "<span title=\"/x/runA/configuration.yaml\">runA/configuration.yaml</span>"
              @test SparlectraApp._webui_compare_value_cell("active_set", "classic_simultaneous") == "active_set"
              @test occursin("power_flow.qlimits.enforcement_mode", html)
              @test occursin("The runs clamped different buses", html)
              @test occursin("max |dVm|", html)
              @test occursin("0.01", html)          # 0.99 gegen 0.98
              @test occursin("NEWBERRY 2", html)
              @test occursin("Losses", html)
              @test occursin("0.1 MW", html) || occursin("0.09999", html)  # 0.75 gegen 0.65
              @test occursin("/powerflow/result/run-a", html)
          end)() end

        # Reported: "Case input format" offered
        # MATPOWER, DTF and CGMES but neither SCF nor power-grid-model, although
        # the API has accepted `scf` all along and a PGM `input.json` is read by
        # that very importer.
        @testset "case input format offers SCF and power-grid-model" begin (function ()
            @test SparlectraApp._normalize_case_format(:scf) === :scf
            # pgm is a spelling of scf, not a second reader
            @test SparlectraApp._normalize_case_format(:pgm) === :scf
            @test SparlectraApp._normalize_case_format("pgm") === :scf
            @test_throws ArgumentError SparlectraApp._normalize_case_format("nonsense")

            root = mktempdir()
            cases = mktempdir()
            write(joinpath(cases, "fmt_probe.m"), "function mpc = fmt_probe\nmpc.baseMVA = 100;\n")
            rt = (; case_directory=cases, config_file=Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log=SparlectraApp.webui_operation_log_path(root), startup_config_error=nothing, runner=SparlectraApp.start_powerflow_run)
            body = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case"; output_root=root, runtime=rt).body)
            @test occursin("value=\"scf\"", body)
            @test occursin("value=\"pgm\"", body)

            # the two whitelists that persist the choice (save, read back) both have
            # to know the new values, or the selection is silently dropped to auto
            write(joinpath(cases, "fmt_probe_json.json"), "{\"version\": \"1.0\", \"type\": \"input\", \"data\": {}}\n")
            SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save",
                Dict{String,Any}("casefile" => "fmt_probe_json.json", "case_format" => "pgm", "return_to" => "case"); output_root=root, runtime=rt)
            @test SparlectraApp._webui_case_form_defaults("fmt_probe_json.json", cases)["case_format"] == "pgm"
            SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save",
                Dict{String,Any}("casefile" => "fmt_probe.m", "case_format" => "matpower", "return_to" => "case"); output_root=root, runtime=rt)
            @test SparlectraApp._webui_case_form_defaults("fmt_probe.m", cases)["case_format"] == "matpower"
            # the selector offers PowSyBl as well; its choice survives the save
            cp(joinpath(Sparlectra.SPARLECTRA_ROOT, "data", "powsybl", "ieee14.xiidm"), joinpath(cases, "fmt_probe.xiidm"))
            SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save",
                Dict{String,Any}("casefile" => "fmt_probe.xiidm", "case_format" => "powsybl", "return_to" => "case"); output_root=root, runtime=rt)
            @test get(SparlectraApp._webui_case_form_defaults("fmt_probe.xiidm", cases), "case_format", nothing) == "powsybl"

            # Seen on Windows: `scf` stored for a MATPOWER case sent
            # the .m file into the JSON reader (`invalid integer ""`), while the
            # Case page showed nothing wrong. A stored format that contradicts
            # the file CONTENT is not loaded: the run falls back to auto, the
            # operation log and the Case page say so, and the next save writes
            # what the selector shows.
            Sparlectra.write_case_config(joinpath(cases, "fmt_probe.m"), Dict{String,Any}(); form=Dict{String,Any}("case_format" => "scf"))
            stored = SparlectraApp._webui_case_form_defaults("fmt_probe.m", cases)
            @test !haskey(stored, "case_format")
            @test occursin("fmt_probe.m", stored["_case_format_conflict"])
            req = SparlectraApp.powerflow_webui_request(Dict("casefile" => "fmt_probe.m"); case_directory=cases, operation_log=root)
            @test req["case_format"] == "auto"
            log_lines = readlines(SparlectraApp.webui_operation_log_path(root))
            ignored = filter(l -> occursin("case_format_ignored", l), log_lines)
            @test length(ignored) == 1
            @test occursin("\"stored_format\":\"scf\"", ignored[1]) || occursin("\"stored_format\": \"scf\"", ignored[1])
            body = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case?casefile=fmt_probe.m"; output_root=root, runtime=rt).body)
            @test occursin("case-format-notice", body)
            @test occursin("does not fit this file", body)
            # a value that fits the file is loaded and shows no notice
            Sparlectra.write_case_config(joinpath(cases, "fmt_probe.m"), Dict{String,Any}(); form=Dict{String,Any}("case_format" => "matpower"))
            @test SparlectraApp._webui_case_form_defaults("fmt_probe.m", cases)["case_format"] == "matpower"
            body_ok = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case?casefile=fmt_probe.m"; output_root=root, runtime=rt).body)
            @test !occursin("case-format-notice", body_ok)
            # the explicit wrong format posted by a form fails with the
            # readable message, not with a JSON parse error
            req_scf = SparlectraApp.powerflow_webui_request(Dict("casefile" => "fmt_probe.m", "case_format" => "scf"); case_directory=cases, operation_log=root)
            @test req_scf["case_format"] == "scf"
            failed = SparlectraApp.run_sparlectra_api(casefile = joinpath(cases, "fmt_probe.m"), case_format = :scf, output_dir = joinpath(root, "scf_on_m"))
            @test failed.status === :failed
            @test failed.reason == "case_format_mismatch"
            @test occursin("auto", String(failed.message))
            @test !occursin("invalid integer", String(failed.message))
        end)() end

        # Also reported live: the mode was set and nothing happened, because the
        # enable checkbox sat elsewhere in the form and stayed off. Both controls
        # now live in one block, and the mode itself can say "off".
        @testset "Q-limit handling reads as one setting" begin (function ()
            root = mktempdir()
            rt = (; case_directory=mktempdir(), config_file=Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log=SparlectraApp.webui_operation_log_path(root), startup_config_error=nothing, runner=SparlectraApp.start_powerflow_run)
            body = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/settings"; output_root=root, runtime=rt).body)
            block_start = findfirst("Q-limit handling", body)
            @test block_start !== nothing
            block = body[block_start[1]:min(lastindex(body), block_start[1]+1200)]
            @test occursin("power_flow_qlimits_enabled", block)
            @test occursin("value=\"off\"", block)

            # "off" disables the handling instead of being sent as a mode the
            # configuration would reject
            form = Dict{String,String}(
                "casefile" => "case14.m",
                "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
                "power_flow_qlimits_enforcement_mode" => "off",
                "power_flow_qlimits_enabled" => "true",
            )
            overrides = get(SparlectraApp.powerflow_webui_request(form), "config_overrides", Dict{String,Any}())
            @test !haskey(overrides, "power_flow.qlimits.enforcement_mode")
            @test overrides["power_flow.qlimits.enabled"] === false

            # a real mode is passed through untouched
            form["power_flow_qlimits_enforcement_mode"] = "classic_simultaneous"
            overrides2 = get(SparlectraApp.powerflow_webui_request(form), "config_overrides", Dict{String,Any}())
            @test overrides2["power_flow.qlimits.enforcement_mode"] == "classic_simultaneous"

            # The same word saved from the settings page ("Q-Limit
            # an, aber mode auf off" ended in an ArgumentError on the next run):
            # the case file and the configuration file get `enabled: false` and
            # no mode key.
            cache = rt.case_directory
            SparlectraApp._webui_stage_bundled_case!(normpath(dirname(@__DIR__)), cache, "sp_case14.scf.json")
            off_form = Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "power_flow_qlimits_enabled" => "true", "power_flow_qlimits_enforcement_mode" => "off")
            resp_case = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", off_form; output_root=root, runtime=rt)
            @test resp_case.status in (302, 303)
            saved = Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json"))
            @test saved["power_flow.qlimits.enabled"] === false
            @test !haskey(saved, "power_flow.qlimits.enforcement_mode")
            cfg_general = joinpath(root, "configuration.yaml")
            cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, cfg_general)
            resp_general = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("settings_target" => "general", "config_file" => cfg_general, "power_flow_qlimits_enabled" => "true", "power_flow_qlimits_enforcement_mode" => "off"); output_root=root)
            @test resp_general.status in (302, 303)
            general_cfg = Sparlectra.load_sparlectra_config(cfg_general; reload=true)
            @test general_cfg.powerflow.qlimits.ignore_q_limits
            @test general_cfg.powerflow.qlimits.enforcement_mode === :active_set
            # the case-options save follows the same rule
            opts_form = Dict{String,Any}("casefile" => "sp_case14.scf.json", "power_flow_qlimits_enforcement_mode" => "off", "return_to" => "case")
            SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", opts_form; output_root=root, runtime=rt)
            @test !haskey(Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json")), "power_flow.qlimits.enforcement_mode")
        end)() end

        # A diagnose run takes ONE step from the case's own voltages, so it never
        # converges by design. Reported as a failed power flow it was unusable:
        # every diagnosis looked like a crash.
        @testset "a completed diagnose run reads as a diagnosis" begin (function ()
            probe = Dict{String,Any}("run_mode" => "diagnose", "status" => "not_converged", "success" => false)
            @test SparlectraApp._webui_is_completed_diagnose(probe)
            @test SparlectraApp.webui_status_class(probe) == "status-info"

            # a result page has no run_mode; the self-check configuration the
            # diagnose flow writes is the marker there
            by_artifact = Dict{String,Any}("status" => "not_converged",
                "artifacts" => [Dict{String,Any}("name" => "diagnose_self_check_config.yaml")])
            @test SparlectraApp._webui_is_completed_diagnose(by_artifact)

            # a diagnose run that could NOT run keeps the failure vocabulary
            broken = Dict{String,Any}("run_mode" => "diagnose", "status" => "failed", "success" => false)
            @test !SparlectraApp._webui_is_completed_diagnose(broken)
            @test SparlectraApp.webui_status_class(broken) == "status-error"

            # an ordinary non-converged run stays an error
            plain = Dict{String,Any}("run_mode" => "", "status" => "not_converged", "success" => false)
            @test !SparlectraApp._webui_is_completed_diagnose(plain)
            @test SparlectraApp.webui_status_class(plain) == "status-error"
        end)() end

        @testset "sysimage page tells a native session how to use the image" begin (function ()
            # A Web UI started from the REPL runs without the image and cannot
            # switch to it, not even after a build from its own Sysimage page;
            # the page has to say how the image on disk is used.
            mktempdir() do tmp
                out = joinpath(tmp, "runs")
                mkpath(out)
                img = SparlectraApp.webui_sysimage_path(out)
                mkpath(dirname(img))
                # without an image on disk there is no hint; the status says why instead
                @test SparlectraApp.sysimage_use_hint(img; flavor_kind=:native, problem="no sysimage found") === nothing
                @test !occursin("sysimage-native-hint", SparlectraApp.render_webui_sysimage_page(output_root=out))
                # a valid image for THIS checkout: metadata matching the running
                # Julia and the application manifest, newer than every source file
                manifest = joinpath(SparlectraApp.SPARLECTRA_APP_ROOT, "Manifest.toml")
                sha = bytes2hex(open(SHA.sha256, manifest))
                write(joinpath(dirname(img), "sysimage_meta.toml"), string("julia_version = \"", VERSION, "\"\nmanifest_sha256 = \"", sha, "\"\n"))
                write(img, "not a real image")
                @test SparlectraApp.webui_sysimage_problem(image_path=img) === nothing
                hint = SparlectraApp.sysimage_use_hint(img; flavor_kind=:native, problem=nothing)
                @test hint !== nothing
                # the line mirrors what the launcher runs on a relaunch
                @test occursin("-J \"$(img)\"", hint.command)
                @test occursin("--startup-file=no", hint.command)
                @test occursin("--project=\"$(SparlectraApp.SPARLECTRA_APP_ROOT)\"", hint.command)
                @test occursin("start_sparlectra_webui", hint.repl)
                @test hint.start_script == normpath(joinpath(SparlectraApp.SPARLECTRA_APP_ROOT, ".."))
                @test isfile(joinpath(hint.start_script, "start_webui.jl"))
                # a session on the image or in the standalone app gets no hint
                @test SparlectraApp.sysimage_use_hint(img; flavor_kind=:sysimage, problem=nothing) === nothing
                @test SparlectraApp.sysimage_use_hint(img; flavor_kind=:app, problem=nothing) === nothing
                # the page renders it for a native session (this test process)
                if SparlectraApp.webui_runtime_flavor().kind === :native
                    page = SparlectraApp.render_webui_sysimage_page(output_root=out)
                    @test occursin("sysimage-native-hint", page)
                    @test occursin(SparlectraApp._webui_escape(img), page)
                    println("  sysimage native hint: page check ran")
                else
                    println("  sysimage native hint: page check skipped (this process runs on an image)")
                end
            end
        end)() end

        @testset "sysimage build progress and rebuild page" begin (function ()
            mktempdir() do tmp
                # webui_sysimage_dir puts the image NEXT TO the runs root, so an
                # output_root inside the temp directory keeps the whole fixture there
                out = joinpath(tmp, "runs")
                mkpath(out)
                imgdir = SparlectraApp.webui_sysimage_dir(out)
                mkpath(imgdir)
                progress = SparlectraApp.webui_sysimage_progress_path(out)

                # no build has ever run here
                @test SparlectraApp.read_sysimage_build_progress(output_root=out) === nothing
                @test SparlectraApp.sysimage_build_active(output_root=out) == false

                _progress_file(state, updated) = write(progress, string(
                    "state = \"", state, "\"\nstep = 3\nsteps = 4\n",
                    "phase = \"compiling the system image\"\ndetail = \"\"\nmessage = \"\"\n",
                    "elapsed_seconds = 42.0\nupdated_at = ", updated, "\npid = 1\n",
                    "log = \"", replace(SparlectraApp.webui_sysimage_build_log_path(out), "\\" => "\\\\"), "\"\n"))

                _progress_file("running", round(time(); digits=1))
                read_back = SparlectraApp.read_sysimage_build_progress(output_root=out)
                @test read_back !== nothing
                @test read_back["phase"] == "compiling the system image"
                @test SparlectraApp.sysimage_build_active(output_root=out) == true

                # a builder that stopped reporting is gone, not busy: a crashed build
                # must never lock the refresh button for the rest of the session
                _progress_file("running", round(time() - 3600; digits=1))
                @test SparlectraApp.sysimage_build_active(output_root=out) == false
                _progress_file("failed", round(time(); digits=1))
                @test SparlectraApp.sysimage_build_active(output_root=out) == false

                # the page renders in every state and offers the button when idle
                _progress_file("done", round(time(); digits=1))
                idle_page = SparlectraApp.render_webui_sysimage_page(output_root=out)
                @test occursin("/webui/sysimage/rebuild", idle_page)
                # the ATTRIBUTE with a value, not the bare name: the page skeleton
                # always carries `main[data-refresh-url]` inside its polling script,
                # so the bare name matches every page and proves nothing
                @test !occursin("data-refresh-url=", idle_page)
                _progress_file("running", round(time(); digits=1))
                busy_page = SparlectraApp.render_webui_sysimage_page(output_root=out)
                # while a build runs the page polls instead of offering a second build
                @test occursin("data-refresh-url=\"/webui/sysimage?autorefresh=1\"", busy_page)
                @test !occursin("/webui/sysimage/rebuild", busy_page)
                @test occursin("compiling the system image", busy_page)

                # a failed build shows the reason and the tail of the build log
                write(SparlectraApp.webui_sysimage_build_log_path(out), "line one\nERROR: linker died\n")
                _progress_file("failed", round(time(); digits=1))
                failed_page = SparlectraApp.render_webui_sysimage_page(output_root=out)
                @test occursin("ERROR: linker died", failed_page)

                # A progress file that exists but does not parse must not read as
                # "no build has ever run here": the reader returns nothing for both,
                # and only the page can tell them apart.
                write(progress, "kaputt = [[[")
                @test SparlectraApp.read_sysimage_build_progress(output_root=out) === nothing
                @test occursin("cannot be read", SparlectraApp.render_webui_sysimage_page(output_root=out))
                rm(progress)
                @test !occursin("cannot be read", SparlectraApp.render_webui_sysimage_page(output_root=out))
                _progress_file("failed", round(time(); digits=1))

                # The seconds between the button and the build process reporting for
                # the first time: the page only polls while a build is
                # active, and "active" used to require the CHILD's first entry. Julia
                # needs seconds to boot before it can write one, tens of seconds on a
                # cold Windows, and the page sat there looking dead meanwhile.
                rm(progress; force=true)
                SparlectraApp._write_sysimage_build_progress(out; state="starting", phase="starting the build process")
                @test SparlectraApp.sysimage_build_active(output_root=out) == true
                starting_page = SparlectraApp.render_webui_sysimage_page(output_root=out)
                @test occursin("data-refresh-url=\"/webui/sysimage?autorefresh=1\"", starting_page)
                @test occursin("starting the build process", starting_page)
                @test !occursin("/webui/sysimage/rebuild", starting_page)
                # The state must be `starting`, never `running`: tools/build_sysimage.jl
                # refuses to start when it finds a fresh `running` entry, so writing
                # `running` here would make the Web UI block the very build it just
                # launched.
                @test SparlectraApp.read_sysimage_build_progress(output_root=out)["state"] == "starting"
                # a start that never produced a process is gone after the wider window
                _progress_file("starting", round(time() - 3 * 3600; digits=1))
                @test SparlectraApp.sysimage_build_active(output_root=out) == false

                # route smoke: the page is reachable
                response = SparlectraApp.route_sparlectra_webui("GET", "/webui/sysimage"; output_root=out)
                @test response.status == 200
                @test occursin("Sysimage", String(response.body))

                # The POST is exercised against the "already running" guard ON
                # PURPOSE: start_sysimage_rebuild! spawns the real build script,
                # which writes to the user's own sysimage directory and works for
                # minutes. A test may not do that, so what is asserted here is the
                # wiring (route reaches the handler, handler redirects back with a
                # message) and the guard that stops a second build. The spawn itself
                # is covered by running a build, not by the suite.
                # the button starts the repository build script; 0.50.0 looked
                # for it inside app/ and failed
                @test isfile(SparlectraApp._sysimage_build_script())
                _progress_file("running", round(time(); digits=1))
                rebuild = SparlectraApp.route_sparlectra_webui("POST", "/webui/sysimage/rebuild"; output_root=out)
                @test rebuild.status == 303
                @test any(h -> first(h) == "Location" && startswith(last(h), "/webui/sysimage"), rebuild.headers)
                @test any(h -> first(h) == "Location" && occursin("already", last(h)), rebuild.headers)
            end
        end)() end

        @testset "update check and install" begin (function ()
            A = SparlectraApp
            # every fetch is a stub of the same shape as the Downloads one, so
            # no test reaches the network
            # the fetch follows the redirect of the latest-release page; `url`
            # is the address it ended at
            answer(status, final_url = ""; headers = Dict{String,String}()) = (url; timeout) -> (status = status, headers = headers, url = final_url, error = "", timed_out = false)
            failure(msg; timed_out = false) = (url; timeout) -> (status = 0, headers = Dict{String,String}(), url = "", error = msg, timed_out = timed_out)
            release(tag) = "https://github.com/Welthulk/Sparlectra.jl/releases/tag/$(tag)"
            cur = v"0.50.0"
            check(fetch) = A.webui_check_for_update(; fetch, current = cur)
            @test check(answer(200, release("v0.50.0"))).state === :up_to_date
            # a checkout ahead of the latest release is up to date as well
            @test check(answer(200, release("v0.49.3"))).state === :up_to_date
            newer = check(answer(200, release("v0.51.0")))
            @test (newer.state, newer.latest, newer.tag) == (:newer, v"0.51.0", "v0.51.0")
            @test newer.release_url == "https://github.com/Welthulk/Sparlectra.jl/releases/tag/v0.51.0"
            # the check asks the release page, never the REST API with its
            # anonymous quota of 60 requests per hour per IP
            asked = String[]
            A.webui_check_for_update(; fetch = (url; timeout) -> (push!(asked, url); answer(200, release("v0.50.0"))(url; timeout)), current = cur)
            @test asked == ["https://github.com/Welthulk/Sparlectra.jl/releases/latest"]
            offline = check(failure("Could not resolve host: github.com"))
            @test (offline.state, offline.failure) == (:failed, :network)
            @test occursin("Could not resolve host", offline.message)
            @test check(failure("Operation timed out"; timed_out = true)).failure === :timeout
            @test check(answer(429)).failure === :rate_limit
            @test check(answer(500)).failure === :http
            @test check(answer(404)).failure === :http
            # an address that is not a release tag of the repository: no
            # exception, a named failure (a repository without a release ends
            # on /releases, a proxy on its login page)
            for final in ("", "https://github.com/Welthulk/Sparlectra.jl/releases", "https://proxy.example.com/login", "https://github.com/Welthulk/Sparlectra.jl/releases/tag/latest", "https://github.com/Welthulk/Sparlectra.jl/releases/tag/v0.51", "https://github.com/Other/Sparlectra.jl/releases/tag/v0.51.0")
                @test (final, check(answer(200, final)).failure) == (final, :malformed)
            end

            mktempdir() do tmp
                # --- layout detection
                checkout = joinpath(tmp, "dev", "Sparlectra")
                mkpath(joinpath(checkout, ".git"))
                @test A.webui_install_layout(checkout; flavor_kind = :native) === :git_checkout
                worktree = joinpath(tmp, "wt")
                mkpath(worktree)
                write(joinpath(worktree, ".git"), "gitdir: /elsewhere/.git/worktrees/wt\n")
                @test A.webui_install_layout(worktree; flavor_kind = :native) === :git_checkout
                inst = joinpath(tmp, "apps", "Sparlectra")
                mkpath(joinpath(inst, "tools"))
                for f in ("Project.toml", "start_webui.jl", "start_webui.sh", "start_webui.bat", joinpath("tools", "install_webui.sh"), joinpath("tools", "install_webui.bat"))
                    write(joinpath(inst, f), "# fixture\n")
                end
                # SPARLECTRA_ROOT carries a trailing separator
                @test A.webui_install_layout(inst * "/"; flavor_kind = :sysimage) === :installed_copy
                @test A.webui_install_layout(inst; flavor_kind = :app) === :unknown
                renamed = joinpath(tmp, "apps", "Sparlectra.jl-0.50.0")
                cp(inst, renamed)
                @test A.webui_install_layout(renamed; flavor_kind = :native) === :unknown
                no_script = joinpath(tmp, "other", "Sparlectra")
                mkpath(dirname(no_script))
                cp(inst, no_script)
                rm(joinpath(no_script, "tools"); recursive = true)
                @test A.webui_install_layout(no_script; flavor_kind = :native) === :unknown

                # --- the unattended answers follow the current state
                home = joinpath(tmp, "home")
                mkpath(home)
                img = joinpath(tmp, "state", "sysimage", "sparlectra.so")
                none_yet = A._webui_update_answers(; flavor_kind = :native, sysimage_path = img, app_root = joinpath(inst, "app"), home, windows = false, apple = false)
                @test none_yet == (build_sysimage = false, create_shortcut = false)
                mkpath(joinpath(home, ".local", "share", "applications"))
                write(joinpath(home, ".local", "share", "applications", "sparlectra-webui.desktop"), "[Desktop Entry]\n")
                mkpath(dirname(img))
                write(img, "image")
                mkpath(joinpath(inst, "app"))
                write(joinpath(dirname(img), "sysimage_meta.toml"), "project_dir = \"$(escape_string(joinpath(inst, "app")))\"\n")
                @test A._webui_update_answers(; flavor_kind = :native, sysimage_path = img, app_root = joinpath(inst, "app"), home, windows = false, apple = false) == (build_sysimage = true, create_shortcut = true)
                # the image of another checkout is not this copy's: no build
                @test A._webui_update_answers(; flavor_kind = :native, sysimage_path = img, app_root = joinpath(checkout, "app"), home, windows = false, apple = false).build_sysimage == false
                @test A._webui_update_answers(; flavor_kind = :sysimage, sysimage_path = joinpath(tmp, "none.so"), app_root = joinpath(inst, "app"), home, windows = false, apple = false).build_sysimage == true

                # --- the check page offers what the layout allows
                logs = joinpath(tmp, "logs")
                out = joinpath(tmp, "runs")
                mkpath(out)
                ans = (build_sysimage = true, create_shortcut = false)
                page(layout, fetch; root = inst, started = true) = String(A.handle_webui_update_check(; operation_log = logs, fetch, layout, root, current = cur, output_root = out, started_by_script = started, answers = ans).body)
                git_page = page(:git_checkout, answer(200, release("v0.51.0")); root = checkout)
                @test occursin("git fetch --tags", git_page) && occursin("git checkout v0.51.0", git_page)
                @test !occursin("/webui/update/install", git_page)
                unknown_page = page(:unknown, answer(200, release("v0.51.0")); root = renamed)
                @test occursin(A.WEBUI_UPDATE_ONE_LINE_SH, unknown_page) || occursin(A._webui_escape(A.WEBUI_UPDATE_ONE_LINE_WINDOWS), unknown_page)
                @test !occursin("/webui/update/install", unknown_page)
                inst_page = page(:installed_copy, answer(200, release("v0.51.0")))
                @test occursin("action=\"/webui/update/install\"", inst_page)
                @test occursin(A._webui_escape(A.webui_update_log_path(logs)), inst_page)
                @test occursin(A._webui_escape(A._webui_update_token()), inst_page)
                # a session that does not end on stop cannot update itself
                @test !occursin("/webui/update/install", page(:installed_copy, answer(200, release("v0.51.0")); started = false))
                @test !occursin("/webui/update/install", page(:installed_copy, answer(200, release("v0.50.0"))))
                @test !occursin("/webui/update/install", page(:installed_copy, failure("Could not resolve host: api.github.com")))
                @test !isdir(A.webui_update_dir(logs))
                # the Info panel carries the button, and the route reaches the check
                @test occursin("action=\"/webui/update/check\"", A.render_webui_error(404, "x"))
                saved_fetch = A._WEBUI_UPDATE_FETCH[]
                try
                    A._WEBUI_UPDATE_FETCH[] = answer(200, release("v0.51.0"))
                    @test A.route_sparlectra_webui("POST", "/webui/update/check"; output_root = out).status == 200
                finally
                    A._WEBUI_UPDATE_FETCH[] = saved_fetch
                end

                # --- the install action
                spawned = Cmd[]
                spawn = cmd -> (push!(spawned, cmd); nothing)
                stopped = Ref(0)
                good = Dict{String,String}("confirm_update" => "yes", "token" => A._webui_update_token(), "tag" => "v0.51.0")
                no_job() = nothing
                busy_job() = Dict{String,Any}("status" => "running")
                install(form; layout = :installed_copy, started = true, active_job = no_job) = A.handle_webui_update_install(form; operation_log = logs, shutdown = () -> (stopped[] += 1), spawn, layout, root = inst, output_root = out, started_by_script = started, active_job, answers = ans, windows = false)
                # a git checkout and an unknown layout never start a process
                @test install(good; layout = :git_checkout).status == 409
                @test install(good; layout = :unknown).status == 409
                @test install(good; started = false).status == 409
                @test install(delete!(copy(good), "confirm_update")).status == 400
                @test install(merge(good, Dict("token" => "stale"))).status == 403
                @test install(good; active_job = busy_job).status == 409
                @test isempty(spawned) && stopped[] == 0
                @test !isdir(A.webui_update_dir(logs))
                # through the route: this checkout is a git checkout, so nothing starts
                @test A.route_sparlectra_webui("POST", "/webui/update/install", good; output_root = out).status == 409
                @test !isdir(A.webui_update_dir(out))

                accepted = install(good)
                @test accepted.status == 200
                @test length(spawned) == 1
                timedwait(() -> stopped[] == 1, 5.0)
                @test stopped[] == 1
                # the command that would run, built and handed to the stub only
                update_dir = A.webui_update_dir(logs)
                cmd = only(spawned)
                updater = cmd.exec[2]
                @test cmd.exec[1] == "sh"
                @test dirname(updater) == update_dir
                @test occursin(r"^sparlectra_update_\d{8}_\d{6}\.sh$", basename(updater))
                @test cmd.dir == update_dir
                script = read(updater, String)
                @test occursin("kill -0 $(getpid())", script)
                for line in ("SPARLECTRA_UPDATE=1", "SPARLECTRA_UPDATE_TAG='v0.51.0'", "SPARLECTRA_BUILD_SYSIMAGE=1", "SPARLECTRA_CREATE_SHORTCUT=0", "SPARLECTRA_START_WEBUI=0", "unset SPARLECTRA_SYSIMAGE_CHECKED", "SPARLECTRA_INSTALL_DIR='$(dirname(inst))'", "sh '$(joinpath(update_dir, "install_webui.sh"))'", "exec sh '$(joinpath(inst, "start_webui.sh"))'")
                    @test (line, occursin(line, script)) == (line, true)
                end
                @test read(joinpath(update_dir, "install_webui.sh"), String) == "# fixture\n"
                @test startswith(read(joinpath(update_dir, "update.log"), String), "Sparlectra update requested")
                if Sys.iswindows()
                    println("  update: sh -n of the updater skipped (Windows)")
                else
                    @test success(`sh -n $(updater)`)
                    println("  update: sh -n of the updater ran")
                end
                @test occursin("\"status\":\"started\"", read(A.webui_operation_log_path(logs), String))

                # --- the question of start_webui.jl: only an installed copy
                # is asked, and only y or yes starts the updater
                clogs = joinpath(tmp, "clogs")
                cspawned = Cmd[]
                asked = String[]
                fetched = Ref(0)
                counted(f) = (url; timeout) -> (fetched[] += 1; f(url; timeout))
                function offer(; reply = "y", fetch = answer(200, release("v0.51.0")), layout = :installed_copy, root = inst, interactive = true, args = String[], env = Dict{String,String}(), port = false, building = false, spawn = cmd -> (push!(cspawned, cmd); nothing), ask = p -> (push!(asked, p); reply))
                    return A.webui_console_update_offer(; args, env, interactive, fetch = counted(fetch), current = cur, layout, root, ask, port_busy = () -> port, sysimage_busy = () -> building, operation_log = clogs, answers = ans, spawn, io = devnull, windows = false)
                end
                # (a) a newer release, an installed copy and a yes
                @test offer(reply = "y") === :started
                @test length(cspawned) == 1 && length(asked) == 1
                @test cspawned[1].exec[1] == "sh" && dirname(cspawned[1].exec[2]) == A.webui_update_dir(clogs)
                @test occursin("kill -0 $(getpid())", read(cspawned[1].exec[2], String))
                @test offer(reply = " YES\n") === :started
                @test length(cspawned) == 2
                # (b) every other answer, an empty line, the end of the input
                # and no answer in time are a no
                for reply in ("n", "no", "", "  ", "ja", "yy")
                    @test (reply, offer(; reply)) == (reply, :declined)
                end
                eof_ask = p -> A._webui_console_read_answer(p; io = devnull, input = IOBuffer(""))
                @test offer(ask = eof_ask) === :declined
                silent = Base.BufferStream()
                @test A._webui_console_read_answer("?"; io = devnull, input = silent, seconds = 0.2) == ""
                typed = Base.BufferStream()
                write(typed, "y\n")
                @test A._webui_console_read_answer("?"; io = devnull, input = typed, seconds = 5.0) == "y"
                @test length(cspawned) == 2
                n_asked = length(asked)
                # (c) up to date and (d) a failed check ask nothing
                @test offer(fetch = answer(200, release("v0.50.0"))) === :up_to_date
                @test offer(fetch = failure("Could not resolve host: api.github.com")) === :failed
                @test offer(fetch = failure("Operation timed out"; timed_out = true)) === :failed
                @test offer(fetch = answer(500, "")) === :failed
                # (e) a git checkout or an unknown copy is never asked
                @test offer(layout = :git_checkout, root = checkout) === :manual
                @test offer(layout = :unknown, root = renamed) === :manual
                # a Web UI on the port or a sysimage build: no question
                @test offer(port = true) === :busy
                @test offer(building = true) === :busy
                @test length(asked) == n_asked && length(cspawned) == 2
                # no terminal and the opt-outs: not even the request
                fetched[] = 0
                @test offer(interactive = false) === :not_interactive
                @test offer(args = ["--no-update-check"]) === :disabled
                @test offer(env = Dict("SPARLECTRA_UPDATE_CHECK" => "0")) === :disabled
                @test offer(env = Dict("JULIA_PKG_OFFLINE" => "true")) === :disabled
                @test fetched[] == 0 && length(asked) == n_asked && length(cspawned) == 2
                # a shell that cannot start: nothing runs, the start goes on
                @test offer(spawn = cmd -> throw(Base.IOError("could not spawn sh", -2))) === :start_failed
                @test occursin("\"method\":\"CONSOLE\"", read(A.webui_operation_log_path(clogs), String))

                # the Windows plan, built on any platform: start opens the console
                # window, the batch file has CRLF line ends and the same answers
                wplan = A.webui_update_plan(; root = inst, current = cur, tag = "v0.51.0", operation_log = joinpath(tmp, "wlogs"), pid = 4242, answers = (build_sysimage = false, create_shortcut = true), windows = true)
                @test wplan.command.exec[2:end] == ["/c", "start", "\"Sparlectra update\"", "\"$(wplan.updater)\""]
                @test wplan.command.flags & Base.UV_PROCESS_WINDOWS_VERBATIM_ARGUMENTS != 0
                bat = read(wplan.updater, String)
                @test count("\r\n", bat) == count("\n", bat)
                for line in ("tasklist /FI \"PID eq 4242\"", "set \"SPARLECTRA_UPDATE_TAG=v0.51.0\"", "set \"SPARLECTRA_BUILD_SYSIMAGE=0\"", "set \"SPARLECTRA_CREATE_SHORTCUT=1\"", "set \"SPARLECTRA_INSTALL_DIR=$(dirname(inst))\"", "call \"$(joinpath(inst, "start_webui.bat"))\"")
                    @test (line, occursin(line, bat)) == (line, true)
                end
                # the restart runs in a console with a keyboard: the update
                # question is switched off before start_webui.bat is called
                restart = findfirst("call \"$(joinpath(inst, "start_webui.bat"))\"", bat)
                optout = findlast("set \"SPARLECTRA_UPDATE_CHECK=0\"", bat)
                @test optout !== nothing && restart !== nothing && last(optout) < first(restart)
                # the Windows updater records the result of the install script
                # and stops with the reason before it restarts the old copy
                for line in ("set \"RC=%ERRORLEVEL%\"", "findstr /b /c:\"Update FAILED:\"", "echo Update FAILED:%REASON%", ":update_done")
                    @test (line, occursin(line, bat)) == (line, true)
                end

                # --- a failed update stays visible: the result file of the
                # updater, read back with CRLF as the Windows updater writes it
                flog = joinpath(tmp, "flogs", "operations.jsonl")
                @test A.webui_last_update(flog) === nothing
                @test A._webui_info_last_update(flog) == ""
                mkpath(A.webui_update_dir(flog))
                write(A.webui_update_status_path(flog), "failed\r\nv0.51.0\r\n download of v0.51.0 failed\r\n2026-10-22 10:00\r\n")
                lastup = A.webui_last_update(flog)
                @test (lastup.ok, lastup.tag, lastup.reason, lastup.finished) == (false, "v0.51.0", "download of v0.51.0 failed", "2026-10-22 10:00")
                @test occursin("The last update to v0.51.0 failed (2026-10-22 10:00): download of v0.51.0 failed.", A._webui_info_last_update(flog))
                @test occursin("update-last-failed", String(A.handle_webui_update_check(; operation_log = flog, fetch = answer(200, release("v0.50.0")), layout = :installed_copy, root = inst, current = cur, output_root = out).body))
                # the start in a terminal names it even with the check switched
                # off, as after the restart of an update
                said = IOBuffer()
                @test A.webui_console_update_offer(; env = Dict("SPARLECTRA_UPDATE_CHECK" => "0"), operation_log = flog, io = said) === :disabled
                @test occursin("The last update to v0.51.0 failed", String(take!(said)))
                write(A.webui_update_status_path(flog), "ok\nv0.51.0\n\n2026-10-22 10:05\n")
                @test A.webui_last_update(flog).ok && A._webui_info_last_update(flog) == ""

                if Sys.iswindows()
                    println("  update: sh installer and updater runs skipped (Windows)")
                else
                    # --- the sh updater: a failing install script leaves the
                    # reason in last_update.txt, a good one "ok"
                    ulog = joinpath(tmp, "ulogs", "operations.jsonl")
                    uplan = A.webui_update_plan(; root = inst, current = cur, tag = "v0.51.0", operation_log = ulog, pid = 999999, answers = (build_sysimage = false, create_shortcut = false), windows = false)
                    write(uplan.installer, "echo \"Update FAILED: download of v0.51.0 failed\"\nexit 1\n")
                    write(joinpath(inst, "start_webui.sh"), "exit 0\n")
                    run(Cmd(`sh $(uplan.updater)`; dir = dirname(uplan.updater)))
                    ufail = A.webui_last_update(ulog)
                    @test (ufail.ok, ufail.tag, ufail.reason) == (false, "v0.51.0", "download of v0.51.0 failed")
                    write(uplan.installer, "exit 0\n")
                    run(Cmd(`sh $(uplan.updater)`; dir = dirname(uplan.updater)))
                    @test A.webui_last_update(ulog).ok

                    # --- the sh installer: an update that was asked for and
                    # cannot finish stops with the reason and exit code 1, the
                    # existing copy stays as it is (no network: curl fails)
                    sandbox = joinpath(tmp, "installer")
                    copy_dir = joinpath(sandbox, "Sparlectra")
                    mkpath(copy_dir)
                    write(joinpath(copy_dir, "Project.toml"), "name = \"Sparlectra\"\nversion = \"0.1.0\"\n")
                    fakebin = joinpath(tmp, "fakebin")
                    mkpath(fakebin)
                    write(joinpath(fakebin, "curl"), "#!/bin/sh\nexit 22\n")
                    chmod(joinpath(fakebin, "curl"), 0o755)
                    script = joinpath(A.SPARLECTRA_ROOT, "tools", "install_webui.sh")
                    # local names: a plain assignment in this closure would
                    # rebind the variables of the enclosing testset function
                    function run_installer(env)
                        local buf = IOBuffer()
                        local base = Dict("PATH" => fakebin * ":" * ENV["PATH"], "HOME" => tmp, "SPARLECTRA_INSTALL_DIR" => sandbox, "SPARLECTRA_START_WEBUI" => "0", "SPARLECTRA_CREATE_SHORTCUT" => "0")
                        local proc = run(pipeline(ignorestatus(setenv(`sh $(script)`, merge(base, env))); stdout = buf, stderr = buf))
                        return (proc.exitcode, String(take!(buf)))
                    end
                    code, text = run_installer(Dict("SPARLECTRA_UPDATE" => "1", "SPARLECTRA_UPDATE_TAG" => "v9.9.9"))
                    @test code == 1
                    @test occursin("Update FAILED: download of v9.9.9 failed", text)
                    code, text = run_installer(Dict("SPARLECTRA_UPDATE" => "1"))
                    @test code == 1
                    @test occursin("Update FAILED: could not determine the latest release", text)
                    @test read(joinpath(copy_dir, "Project.toml"), String) == "name = \"Sparlectra\"\nversion = \"0.1.0\"\n"
                    @test !isdir(joinpath(sandbox, "Sparlectra.old"))
                    println("  update: sh installer and updater runs ran")
                end
            end
        end)() end

        @testset "measurement tables bound by CGMES model id and by the assign switch" begin (function ()
            # The SmallGrid fixture is the CGMES case. The full RDF reader
            # (summarizeCGMES) gives the EQ model id the header-only reader
            # has to find.
            C = Sparlectra.CGMESImporter
            sg = joinpath(@__DIR__, "fixtures", "cgmes_smallgrid", "CGMES_v2.4.15_SmallGridTestConfiguration_BaseCase_Complete_v3.0.0.zip")
            eq = [C.cgmesModelKey(f.model_id) for f in C.summarizeCGMES(path = sg).files if :EQ in f.profiles]
            @test length(eq) == 1 && C.cgmesEquipmentModelIds(sg) == eq
            root = mktempdir()
            cases = joinpath(root, "cases")
            mkpath(cases)
            cp(sg, joinpath(cases, "smallgrid.zip"))
            # a second scenario on the same EQ lists the same tables
            cp(sg, joinpath(cases, "smallgrid_scenario2.zip"))
            # a .zip that is no ZIP has no model id and lists no model-bound table
            write(joinpath(cases, "broken.zip"), "no zip")
            head = "# sparlectra-measurement-table v1\n"
            body = string(join(Sparlectra.MEASUREMENT_TABLE_HEADER, ","), "\n,,,n1,U,110.0,,measured,true\n")
            # the id as a converter may write it: urn:uuid prefix, upper case
            write(joinpath(cases, "by_model.csv"), string(head, "# cgmes_model: urn:uuid:", uppercase(eq[1]), "\n", body))
            write(joinpath(cases, "other_model.csv"), string(head, "# cgmes_model: 00000000-0000-0000-0000-000000000001\n", body))
            # a `# case:` line outranks the model line
            write(joinpath(cases, "case_line.csv"), string(head, "# case: broken.zip\n# cgmes_model: ", eq[1], "\n", body))
            listed = c -> [s.name for s in SparlectraApp._webui_measurement_sets_for_case(cases, c)]
            @test listed("smallgrid.zip") == ["by_model.csv"]
            @test listed("smallgrid_scenario2.zip") == ["by_model.csv"]
            @test listed("broken.zip") == ["case_line.csv"]
            @test only(SparlectraApp._webui_measurement_sets_for_case(cases, "smallgrid.zip")).model == eq[1]

            # assign switch of the SE upload: the stored copy is bound to the
            # case of the form, whatever the file says; CRLF stays CRLF
            upload = (name, text) -> SparlectraApp.WebUICaseUpload(name, Vector{UInt8}(codeunits(text)))
            assign_form = (case, files) -> Dict{String,Any}("return_to" => "stateestimation", "return_case" => case, "assign_to_case" => "true", "casefiles" => files)
            foreign = replace(string(head, "# case: broken.zip\n# cgmes_model: 00000000-0000-0000-0000-000000000001\n", body), "\n" => "\r\n")
            r = SparlectraApp.handle_powerflow_case_import(assign_form("smallgrid.zip", [upload("assigned.csv", foreign)]); output_root = root, case_directory = cases, operation_log = root)
            loc = SparlectraApp._webui_urldecode(Dict(r.headers)["Location"])
            # the assigned file is armed for the next run
            @test occursin("measurement_file=assigned.csv", loc)
            stored = read(joinpath(cases, "assigned.csv"), String)
            @test startswith(stored, "# sparlectra-measurement-table v1\r\n# case: smallgrid.zip\r\n") && count("# case:", stored) == 1
            @test Sparlectra.readMeasurementTable(joinpath(cases, "assigned.csv"))[1].value == 110.0
            @test listed("smallgrid.zip") == ["assigned.csv", "by_model.csv"]
            @test !("assigned.csv" in listed("broken.zip"))
            # assigning again replaces the measurement file of that name
            SparlectraApp.handle_powerflow_case_import(assign_form("smallgrid_scenario2.zip", [upload("assigned.csv", foreign)]); output_root = root, case_directory = cases, operation_log = root)
            @test "assigned.csv" in listed("smallgrid_scenario2.zip") && !("assigned.csv" in listed("smallgrid.zip"))
            # without the switch the file keeps its own binding and an
            # existing name is refused
            plain = Dict{String,Any}("return_to" => "stateestimation", "return_case" => "smallgrid.zip", "casefiles" => [upload("assigned.csv", foreign), upload("plain.csv", foreign)])
            SparlectraApp.handle_powerflow_case_import(plain; output_root = root, case_directory = cases, operation_log = root)
            @test "assigned.csv" in listed("smallgrid_scenario2.zip") && listed("broken.zip") == ["case_line.csv", "plain.csv"]
            # the switch binds measurement files only: a CSV without the
            # version line is refused and not stored
            SparlectraApp.handle_powerflow_case_import(assign_form("smallgrid.zip", [upload("notes.csv", body)]); output_root = root, case_directory = cases, operation_log = root)
            @test !isfile(joinpath(cases, "notes.csv"))
        end)() end

        @testset "case-scan memo cache" begin (function ()
            mktempdir() do dir
                p = joinpath(dir, "SCAN.DAT")
                write(p, _dtf_network_fixture())
                @test SparlectraApp._webui_classify_dat_content_cached(p) === :dtf_network_case
                # second call answers from the (path, mtime, size) memo
                @test SparlectraApp._webui_classify_dat_content_cached(p) === :dtf_network_case
                # a replaced file re-classifies (mtime resolution can be one second)
                sleep(1.1)
                write(p, _dtf_outage_fixture())
                @test SparlectraApp._webui_classify_dat_content_cached(p) === :dtf_outage_file
                # the selector honors the refreshed classification
                @test !SparlectraApp._webui_is_user_selectable_case(p)
            end
        end)() end

        @testset "buildSysimage one-call API (dry run)" begin (function ()
            # the exported one-liner plans against the package project in a child
            # process; dry run must resolve the target paths without building
            r = SparlectraApp.buildSysimage(dry_run=true, quiet=true)
            @test r.built == false
            @test endswith(r.sysimage_path, SparlectraApp.webui_sysimage_ext()) || occursin("sparlectra.", basename(r.sysimage_path))
            @test occursin("sysimage_meta", basename(r.meta_path))
            # the relocatable executable is a checkout tool, not a package function
            @test isfile(joinpath(pkgdir(Sparlectra), "tools", "build_app.jl"))
            # the generated CLI module must parse and carry the core commands
            # (run/se/n1, config flags); a child process keeps ARGS handling and
            # the include of tools/build_app.jl out of the test session
            buildapp_script = joinpath(pkgdir(Sparlectra), "tools", "build_app.jl")
            gen_code = """
              push!(empty!(ARGS), "--dry-run")
              include(raw\"$(buildapp_script)\")
              for fl in (:full, :runtime)
                d = mktempdir()
                Main._write_app_package(d; flavor = fl, script = "")
                src = read(joinpath(d, "src", "SparlectraExe.jl"), String)
                ex = Meta.parseall(src)
                any(a -> a isa Expr && a.head == :error, ex.args) && error("generated module does not parse: " * string(fl))
                for needle in ("run <casefile>", "se <casefile>", "n1 <casefile>", "convert-measurements <casefile> <v1.csv>", "--config=", "--set key=value", "working directory no longer exists", "APP_BUILT_AT")
                  occursin(needle, src) || error("missing CLI needle " * needle * " in flavor " * string(fl))
                end
              end
              println("GEN_OK")
              """
            gen_out = read(`$(Base.julia_cmd()) --startup-file=no -e $(gen_code)`, String)
            @test occursin("GEN_OK", gen_out)
            # header start-flavor: native in the test session, :app when the
            # standalone executable stamped its build time into the environment
            @test SparlectraApp.webui_runtime_flavor().kind in (:native, :sysimage)
            withenv("SPARLECTRA_APP_BUILT" => "2026-08-29T11:11:11") do
                fa = SparlectraApp.webui_runtime_flavor()
                @test fa.kind === :app
                @test fa.built == "2026-08-29T11:11:11"
            end
        end)() end

        @testset "Web UI help pages and search" begin (function ()
            # a control renders its hover hint on the label and, with a doc
            # target, the ? to its in-app help page; that page carries the
            # section shipped with the application and the link into the
            # published documentation (base URL from the configuration)
            withenv("SPARLECTRA_WEBUI_DOCS_BASE_URL" => "https://docs.example/v1/") do
                lbl = SparlectraApp._webui_field_label("power_flow_tol", "Tolerance")
                @test occursin("<span class=\"field-label\" title=\"", lbl)
                @test occursin("class=\"help-link\" href=\"/help/power_flow.tol\"", lbl)
                @test occursin("target=\"_blank\" rel=\"noopener noreferrer\"", lbl)
                @test SparlectraApp.webui_help_doc_url("power_flow.tol") == "https://docs.example/v1/powerflow_configuration/#pf-solver-core"
                help = SparlectraApp.route_sparlectra_webui("GET", "/help/power_flow.tol")
                @test help.status == 200
                help_html = String(help.body)
                @test occursin("class=\"help-hint\"", help_html)
                @test !occursin("<table>", help_html)
                # a library section shows its lead paragraph and says the rest is online
                @test occursin("class=\"field-help help-more\"", help_html)
                @test occursin("href=\"https://docs.example/v1/powerflow_configuration/#pf-solver-core\"", help_html)
                # a relative page link of the section points at the same site
                flat = String(SparlectraApp.route_sparlectra_webui("GET", "/help/power_flow.flatstart").body)
                @test occursin("href=\"https://docs.example/v1/powerflow_configuration/#pf-solver-core\"", flat)
                @test !occursin("SPARLECTRADOCS", flat)
                @test !occursin("(@ref", flat)
            end
            @test SparlectraApp.route_sparlectra_webui("GET", "/help/no.such.topic").status == 404
            # the help ships with the application as two pages, the user
            # help and the reference: rendered with a contents list,
            # labelled headings keep the ids of the published site, links
            # between the two pages stay inside the Web UI
            manual = SparlectraApp.route_sparlectra_webui("GET", "/help")
            @test manual.status == 200
            manual_html = String(manual.body)
            reference = SparlectraApp.route_sparlectra_webui("GET", "/help/reference")
            @test reference.status == 200
            reference_html = String(reference.body)
            @test occursin("class=\"help-toc\"", manual_html)
            @test occursin("<h2 id=\"webui-case-formats\">", manual_html)
            @test occursin("<h2 id=\"start\">", manual_html)
            @test occursin("<h3 id=\"webui-form-options\">", reference_html)
            @test occursin("href=\"#webui-form-options\"", reference_html)
            @test occursin("href=\"/help/reference#webui_powsybl\"", manual_html)
            for page_html in (manual_html, reference_html)
                @test !occursin("SPARLECTRADOCS", page_html)
                @test !occursin("(@id", page_html)
                @test !occursin("(@ref", page_html)
                @test count("<h1", page_html) == 1
                @test occursin("href=\"/help\">Help</a>", page_html)
            end
            # a topic links to its section on the page that carries it
            @test occursin("href=\"/help/reference#webui-flat-start\"", String(SparlectraApp.route_sparlectra_webui("GET", "/help/power_flow.flatstart").body))
            # every target of the search index exists: a heading of its
            # page, or the help page of a control
            page_ids(html) = Set(String(m.captures[1]) for m in eachmatch(r" id=\"([^\"]+)\"", html))
            manual_ids, reference_ids = page_ids(manual_html), page_ids(reference_html)
            index_target_exists(url) = startswith(url, "/help/reference#") ? url[17:end] in reference_ids : startswith(url, "/help#") ? url[7:end] in manual_ids : SparlectraApp.route_sparlectra_webui("GET", url).status == 200
            @test [entry.url for entry in SparlectraApp.WEBUI_HELP_SEARCH_INDEX if !index_target_exists(entry.url)] == String[]
            @test length(SparlectraApp.WEBUI_HELP_SEARCH_INDEX) > 100
            # the in-app links of the two pages land on a heading as well,
            # and so do the help links the form pages carry (the hint of the
            # PowSyBl options pointed at the page its section had left)
            form_pages = join(String(SparlectraApp.route_sparlectra_webui("GET", path; output_root = mktempdir()).body) for path in ("/powerflow/case", "/powerflow", "/powerflow/settings"))
            inapp_links = [String(m.captures[1]) for m in eachmatch(r"href=\"(/help(?:/reference)?#[^\"]+)\"", manual_html * reference_html * form_pages)]
            @test "/help/reference#webui_powsybl" in inapp_links
            @test [url for url in inapp_links if !index_target_exists(url)] == String[]
            # the search: every word must occur, a part of a word is enough,
            # a hit in a title leads, the words are marked in the snippet,
            # the query comes back escaped, no hit names the
            # published documentation
            search_html(q; fragment = false) = String(SparlectraApp.route_sparlectra_webui("GET", "/help/search?q=" * q * (fragment ? "&fragment=1" : "")).body)
            help_hits(q) = [(String(m.captures[1]), String(m.captures[2])) for m in eachmatch(r"<li><a href=\"([^\"]+)\">([^<]*)</a> <span class=\"help-search-kind\">", search_html(q))]
            @test isempty(help_hits(""))
            @test first(help_hits("xiidm")) == ("/help/reference#webui_powsybl", "PowSyBl IIDM files")
            @test first(help_hits("short+circuit")) in (("/help/reference#short-circuit", "Short circuit"), ("/help#short-circuit", "Short circuit"))
            @test ("/help/reference#webui_cgmes_test_configurations", "CGMES test configurations") in help_hits("cgmes:")
            @test ("/help#webui-case-formats", "Case formats at a glance") in help_hits("cgmes:")
            @test first(help_hits("slack+generator")) == ("/help/powsybl_import.slack_ids", "Slack generator ids (PowSyBl)")
            @test ("/help/reference#webui_powsybl", "PowSyBl IIDM files") in help_hits("generatorShortCircuit")
            @test all(hit -> hit in help_hits("slack"), help_hits("slack+generator"))
            @test length(help_hits("slack+generator")) < length(help_hits("slack"))
            @test occursin("<mark>", search_html("xiidm"))
            @test occursin("<mark>xiidm</mark>", lowercase(search_html("XIIDM")))
            @test isempty(help_hits("zzzqqq"))
            withenv("SPARLECTRA_WEBUI_DOCS_BASE_URL" => "https://docs.example/v1/") do
                @test occursin("href=\"https://docs.example/v1/\"", search_html("zzzqqq"; fragment = true))
            end
            @test SparlectraApp.webui_help_search("a") == SparlectraApp.webui_help_search("")
            for hostile in (search_html("%3Cscript%3Ealert(1)%3C%2Fscript%3E"), search_html("%22%3E%3Cimg+src%3Dx%3E"; fragment = true), search_html("a%26b+%3Cb%3E"))
                @test !occursin("<script>alert(1)", hostile)
                @test !occursin("<img src=x>", hostile)
                @test !occursin("<b>\"", hostile)
            end
            @test occursin("&lt;script&gt;alert(1)", search_html("%3Cscript%3Ealert(1)%3C%2Fscript%3E"))
            # the hit list alone is what the search while typing puts into
            # the page; the form of the header is on every page
            @test !occursin("<html", search_html("xiidm"; fragment = true))
            @test occursin("class=\"header-help-search\" method=\"get\" action=\"/help/search\"", String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/history"; output_root = mktempdir()).body))
            # a query on the help page itself answers as the search does
            @test occursin("/help/reference#webui_powsybl", String(SparlectraApp.route_sparlectra_webui("GET", "/help?q=xiidm").body))
            for path in ("/docs", "/docs/webui", "/static/katex/katex.min.js")
                @test SparlectraApp.route_sparlectra_webui("GET", path).status == 404
            end
            # the save-target explanations are small print, not label text
            stg = SparlectraApp.render_settings_page(output_root=mktempdir(), selected_casefile="sp_case14.scf.json")
            @test occursin("settings-target-hint", stg)
            @test !occursin("value=\"this_case\" checked>this case (<code>sp_case14.scf.json</code>): case-scope", stg)
        end)() end

        @testset "sysimage build script dry run parses" begin (function ()
            script = joinpath(Sparlectra.SPARLECTRA_ROOT, "tools", "build_sysimage.jl")
            @test isfile(script)
            # the dry-run flag makes the script side-effect free: no build
            # environment changes, no PackageCompiler; safe to include in-process
            ENV["SPARLECTRA_SYSIMAGE_DRY_RUN"] = "1"
            try
                sandbox = Module()
                text = mktemp() do path, io
                    redirect_stdout(io) do
                        Base.include(sandbox, script)
                    end
                    flush(io)
                    read(path, String)
                end
                @test occursin("dry run", text)
                @test occursin("dry run finished", text)
            finally
                delete!(ENV, "SPARLECTRA_SYSIMAGE_DRY_RUN")
            end
        end)() end

        @testset "SBOM script dry run parses" begin (function ()
            # same mechanism as the sysimage smoke: the dry-run flag keeps the
            # script side-effect free (no build environment, no PkgToSoftwareBOM,
            # no network); the full generation runs only in the release workflow
            script = joinpath(Sparlectra.SPARLECTRA_ROOT, "tools", "generate_sbom.jl")
            @test isfile(script)
            ENV["SPARLECTRA_SBOM_DRY_RUN"] = "1"
            # the script reads the process-global ARGS for its output path; when
            # the test runner itself was invoked with a CLI profile argument
            # (`runtests.jl all`), that argument would leak in as the outfile
            # ("would write: .../all") and break the filename assertion below
            prev_args = copy(ARGS)
            empty!(ARGS)
            try
                sandbox = Module()
                text = mktemp() do path, io
                    redirect_stdout(io) do
                        Base.include(sandbox, script)
                    end
                    flush(io)
                    read(path, String)
                end
                @test occursin("dry run", text)
                @test occursin("Sparlectra.spdx.json", text)
                @test occursin("dry run finished", text)
            finally
                append!(ARGS, prev_args)
                delete!(ENV, "SPARLECTRA_SBOM_DRY_RUN")
            end
        end)() end
    end)() end
end

run_webui_tests() = run_webui_fast_tests()

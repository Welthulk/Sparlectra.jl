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

        @testset "git commit sha resolution (issue #290)" begin (function ()
            sha_a = "a" ^ 40
            sha_b = "b" ^ 40

            mktempdir() do root
                @test SparlectraApp._git_head_commit_sha(root) === nothing

                # Plain checkout: .git directory with a loose ref.
                gitdir = joinpath(root, ".git")
                mkpath(joinpath(gitdir, "refs", "heads"))
                write(joinpath(gitdir, "HEAD"), "ref: refs/heads/main\n")
                write(joinpath(gitdir, "refs", "heads", "main"), sha_a * "\n")
                @test SparlectraApp._git_head_commit_sha(root) == sha_a

                # Detached HEAD: plain SHA in HEAD.
                write(joinpath(gitdir, "HEAD"), sha_b * "\n")
                @test SparlectraApp._git_head_commit_sha(root) == sha_b

                # Packed ref: no loose ref file, SHA only in packed-refs.
                write(joinpath(gitdir, "HEAD"), "ref: refs/heads/packed\n")
                write(joinpath(gitdir, "packed-refs"), "# pack-refs with: peeled fully-peeled sorted\n$(sha_b) refs/heads/packed\n^$(sha_a)\n")
                @test SparlectraApp._git_head_commit_sha(root) == sha_b

                # Regression: the packed-refs handle must be closed on the early-return
                # path above. With a leaked handle this rm fails on Windows (EBUSY,
                # seen as a mktempdir-cleanup error); on Linux it always succeeds, so
                # the guard is deliberately GC- and OS-timing independent.
                rm(joinpath(gitdir, "packed-refs"))
                @test !isfile(joinpath(gitdir, "packed-refs"))
            end

            # Worktree checkout: .git is a "gitdir: <path>" pointer file and shared
            # refs live in the common git directory.
            mktempdir() do dir
                main_git = joinpath(dir, "main", ".git")
                mkpath(joinpath(main_git, "refs", "heads"))
                write(joinpath(main_git, "refs", "heads", "feature"), sha_a * "\n")
                worktree_gitdir = joinpath(main_git, "worktrees", "wt")
                mkpath(worktree_gitdir)
                write(joinpath(worktree_gitdir, "HEAD"), "ref: refs/heads/feature\n")
                write(joinpath(worktree_gitdir, "commondir"), "../..\n")
                worktree_root = joinpath(dir, "wt")
                mkpath(worktree_root)
                write(joinpath(worktree_root, ".git"), "gitdir: $(worktree_gitdir)\n")
                @test SparlectraApp._git_head_commit_sha(worktree_root) == sha_a

                write(joinpath(worktree_root, ".git"), "not a git pointer\n")
                @test SparlectraApp._git_head_commit_sha(worktree_root) === nothing
            end

            # In a git checkout of Sparlectra itself (plain or worktree) the Web UI
            # banner SHA must resolve; a registry install (no .git) yields nothing.
            package_root = normpath(joinpath(dirname(pathof(Sparlectra)), ".."))
            if ispath(joinpath(package_root, ".git"))
                sha = SparlectraApp._sparlectra_git_commit_sha()
                @test sha isa String
                @test occursin(r"^[0-9a-f]{40}$", sha)
            end
        end)() end

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

            active_html = SparlectraApp.render_powerflow_result(Dict("run_id" => "active", "status" => "running", "elapsed_seconds" => 1.25))
            @test occursin("Run status", active_html)
            @test occursin("Elapsed time", active_html)
            @test !occursin("Solver time", active_html)
            @test !occursin("Wall time", active_html)

            terminal = Dict("run_id" => "ok", "status" => "completed", "solver_elapsed_s" => 0.125, "service_phase_timings" => [Dict("phase" => "total_service", "elapsed_seconds" => 1.5)])
            terminal_html = SparlectraApp.render_powerflow_result(terminal)
            @test occursin("Solver time", terminal_html)
            @test occursin("Total time", terminal_html)
            @test !occursin("Elapsed time", terminal_html)
            @test !occursin("Wall time", terminal_html)
            @test !occursin("Solver time: n/a", terminal_html)

            failed_solver_html = SparlectraApp.render_powerflow_result(Dict("run_id" => "failed", "status" => "failed", "metadata" => Dict("solver_elapsed_s" => 0.25), "service_phase_timings" => [Dict("phase" => "total_service", "elapsed_seconds" => 2.0)]))
            @test occursin("Solver time", failed_solver_html)
            presolver_html = SparlectraApp.render_powerflow_result(Dict("run_id" => "bad", "status" => "failed", "service_phase_timings" => [Dict("phase" => "total_service", "elapsed_seconds" => 0.1)]))
            @test !occursin("Solver time", presolver_html)
            @test occursin("Total time", presolver_html)
            legacy_html = SparlectraApp.render_powerflow_result(Dict("run_id" => "old", "status" => "completed", "elapsed_seconds" => 3.0))
            @test occursin("Total time", legacy_html)
        end

        @testset "CGMES export checkbox and result row" begin (function ()
            # the export-CGMES checkbox renders on Settings
            form_html = SparlectraApp.render_settings_page(output_root=mktempdir())
            @test occursin("name=\"export_cgmes\"", form_html)
            @test SparlectraApp.resolve_webui_help_topic("webui.export_cgmes") !== nothing
            # excerpt loading needs the cgmes_export page in WEBUI_DOC_PAGES and the
            # "Export from the Web UI" heading in docs/src/cgmes_export.md
            @test SparlectraApp.load_webui_help_excerpt("webui.export_cgmes") !== nothing

            req = SparlectraApp.powerflow_webui_request(Dict("casefile" => "case.m", "export_cgmes" => "true"))
            @test req["export_cgmes"] === true
            # unchecked checkbox: browsers submit only the hidden false field
            req = SparlectraApp.powerflow_webui_request(Dict("casefile" => "case.m", "export_cgmes" => "false"))
            @test req["export_cgmes"] === false
            # absent field falls back to the spec default (false)
            req = SparlectraApp.powerflow_webui_request(Dict("casefile" => "case.m"))
            @test req["export_cgmes"] === false

            completed_html = SparlectraApp.render_powerflow_result(
                Dict(
                    "run_id" => "x",
                    "status" => "completed",
                    "metadata" => Dict("cgmes_export_status" => "completed", "cgmes_export_files" => "n_EQ.xml, n_TP.xml, n_SSH.xml", "cgmes_export_notices" => "transformer T1: phase shift 2.0° not exported (fixed ratio only)", "cgmes_export_sc_lines" => 3),
                ),
            )
            @test occursin("CGMES export", completed_html)
            @test occursin("n_SSH.xml", completed_html)
            @test occursin("phase shift 2.0° not exported", completed_html)
            @test occursin("zero-sequence data on 3 line(s)", completed_html)

            failed_html = SparlectraApp.render_powerflow_result(Dict("run_id" => "y", "status" => "completed", "metadata" => Dict("cgmes_export_status" => "failed", "cgmes_export_error" => "boom")))
            @test occursin("failed — boom", failed_html)

            plain_html = SparlectraApp.render_powerflow_result(Dict("run_id" => "z", "status" => "completed"))
            @test !occursin("CGMES export", plain_html)
        end)() end

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
            # failed (run f63b75c5)
            resp_gen = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "power_flow_apslf_start_enabled" => "true"); output_root=root, runtime=rt)
            @test resp_gen.status in (302, 303)
            @test Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json"))["power_flow.apslf_start.enabled"] === true
            resp_solver = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "power_flow_solver" => "apslf"); output_root=root, runtime=rt)
            @test resp_solver.status in (302, 303)
            saved_cfg = Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json"))
            @test saved_cfg["power_flow.solver"] == "apslf"
            @test saved_cfg["power_flow.apslf_start.enabled"] === false
            # an incompatible pair in one save is refused at save time, not at the run
            resp_bad = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "power_flow_solver" => "rectangular", "power_flow_apslf_start_enabled" => "true", "power_flow_dc_seed_unconditional" => "true"); output_root=root, runtime=rt)
            @test occursin("Could not save settings for this case", SparlectraApp._webui_urldecode(Dict(resp_bad.headers)["Location"]))
            @test Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json"))["power_flow.solver"] == "apslf"
            SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "power_flow_solver" => "rectangular"); output_root=root, runtime=rt)
            # the flat start is the one start switch: the saved start settings
            # stay as posted, the run switches them off while the flat start is
            # on and names them in run.log; unchecking gives them back
            resp_flat = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "power_flow_flatstart" => "true", "power_flow_apslf_start_enabled" => "true", "power_flow_start_current_iteration_enabled" => "true", "power_flow_start_angle_mode" => "dc", "power_flow_start_voltage_mode" => "profile_blend"); output_root=root, runtime=rt)
            @test resp_flat.status in (302, 303)
            flat_cfg = Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json"))
            @test flat_cfg["power_flow.flatstart"] === true
            @test flat_cfg["power_flow.apslf_start.enabled"] === true
            @test flat_cfg["power_flow.start_current_iteration.enabled"] === true
            run_flat = SparlectraApp.start_powerflow_run(Dict("casefile" => "sp_case14.scf.json", "config_file" => cfg_rt, "output_root" => root); case_directory=cache)
            @test run_flat["status"] == "succeeded"
            run_log = read(joinpath(String(run_flat["output_dir"]), "run.log"), String)
            @test occursin("Flat start: start-value machines forced off for this run: power_flow.apslf_start.enabled=false, power_flow.start_current_iteration.enabled=false, power_flow.start_mode.angle_mode=classic, power_flow.start_mode.voltage_mode=classic", run_log)
            @test occursin(r"Flatstart\s+:\s+Yes", run_log)
            @test run_flat["metadata"]["current_iteration_enabled"] === false
            SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "power_flow_flatstart" => "false"); output_root=root, runtime=rt)
            off_cfg = Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json"))
            @test off_cfg["power_flow.flatstart"] === false
            @test off_cfg["power_flow.apslf_start.enabled"] === true
            # machine-scope keys are named and kept out of the case file
            resp2 = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "benchmark_samples" => "5"); output_root=root, runtime=rt)
            @test occursin("benchmark.samples", SparlectraApp._webui_urldecode(Dict(resp2.headers)["Location"]))
            @test !haskey(Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json")), "benchmark.samples")
            # the general target merges into the YAML with a backup
            cfg = joinpath(root, "configuration.yaml")
            cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, cfg)
            resp3 = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("settings_target" => "general", "config_file" => cfg, "power_flow_autodamp_min" => "0.09"); output_root=root)
            @test resp3.status in (302, 303)
            @test occursin("autodamp_min: 0.09", read(cfg, String))
            @test isfile(cfg * ".settings-save.bak")
        end)() end

        @testset "the provisioned configuration follows changed template defaults" begin (function ()
            # a copy of the template provisioned by an older release still
            # carries the old default of a key the user never touched; it
            # follows the new template value, a value the user set stays
            mktempdir() do dir
                cfg = joinpath(dir, "configuration.yaml")
                text = read(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, String)
                @test occursin("linear_solver: umfpack_reuse", text)
                # an old copy: the old default of the solver, a user choice for the tolerance
                write(cfg, replace(replace(text, "linear_solver: umfpack_reuse" => "linear_solver: umfpack"), r"^  tol: [^\n]*"m => "  tol: 1.0e-7"))
                changed = SparlectraApp._webui_follow_template_defaults!(cfg)
                # the report names key, old and new value for the start message
                @test any(c -> c.key == "power_flow.linear_solver" && c.old == "umfpack" && c.new == "umfpack_reuse", changed)
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

        @testset "the settings page shows the configuration file; the case view is a switch" begin (function ()
            # Reported from the browser: after a save for one case the Settings
            # page kept showing that case's values, and a user could not tell
            # what the configuration file said. The page shows the
            # configuration file's values; ?case_settings=1 overlays the case.
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
            default_view = page("")
            @test max_iter(default_view) == string(Sparlectra.load_sparlectra_config(cfg; reload = true).powerflow.max_iter)
            @test !occursin("case-settings-notice", default_view)
            # the reset of the saved settings stays reachable in both views
            @test occursin("Reset saved settings for this case", default_view)
            case_view = page("&case_settings=1")
            @test max_iter(case_view) == "33"
            @test occursin("case-settings-notice", case_view)
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
            @test occursin("Reset saved settings for this case", with_sidecar)
            @test occursin("/powerflow/case-settings/reset", with_sidecar)
            without = SparlectraApp.render_settings_page(output_root=mktempdir())
            @test !occursin("Reset saved settings for this case", without)
            # switching the case reloads the page server-side; the wait must be
            # visible WHERE the switching happens, which is the Case page's
            # chooser
            @test occursin("case-loading-banner", SparlectraApp.render_case_page(output_root=mktempdir()))

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
            # the configuration carries it, with the documented default
            cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true)
            @test cfg.webui.operation_log_retention_days == 10
            lowered = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true,
                overrides=Dict{String,Any}("webui" => Dict{String,Any}("operation_log_retention_days" => 3)))
            @test lowered.webui.operation_log_retention_days == 3
            @test SparlectraApp._webui_operation_log_options(; retention_days=3).retention_days == 3
            @test_throws ArgumentError Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true,
                overrides=Dict{String,Any}("webui" => Dict{String,Any}("operation_log_retention_days" => -1)))
        end)() end

        @testset "contingency weights editor and storage (#331 Phase 5 follow-up)" begin (function ()
            dir = mktempdir()
            cases = joinpath(dir, "cases")
            mkpath(cases)
            root = joinpath(dir, "runs")
            mkpath(root)
            # load_fixture_net: the tracked PST warmup case gives the editor real
            # element names without any download
            cp(joinpath(pkgdir(Sparlectra), "data", "mpower", "warmup_casePST.m"), joinpath(cases, "warmup_casePST.m"))
            wf = SparlectraApp._webui_case_weights_path("warmup_casePST.m"; case_directory=cases)

            # the weight file lives next to the case as <stem>.contingency-weights.csv
            @test basename(wf) == "warmup_casePST.contingency-weights.csv"
            @test dirname(wf) == normpath(cases)

            # list exclusion: a weights file must not be offered as a selectable case
            touch(wf)
            opts = SparlectraApp._webui_casefile_options_in_directory(cases)
            @test "warmup_casePST.m" in opts
            @test !any(occursin("contingency-weights", o) for o in opts)
            rm(wf)

            # real element names for the fixtures
            net = redirect_stdout(devnull) do
                Sparlectra._import_sparlectra_net(joinpath(cases, "warmup_casePST.m"), nothing, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true))
            end
            bnames = [c.name for c in generateN1Branches(net)]
            bytes = s -> Vector{UInt8}(codeunits(s))
            up = (fname, data) -> SparlectraApp.handle_contingency_weights_upload(Dict{String,Any}("casefile" => "warmup_casePST.m", "casefiles" => [SparlectraApp.WebUICaseUpload(fname, data)]); output_root=root, case_directory=cases, operation_log=root)
            loc = r -> first(p for (k, p) in r.headers if k == "Location")

            # upload a valid weight file
            r = up("w.csv", bytes("name;weight\n$(bnames[1]);3.0\n"))
            @test r.status == 303
            @test isfile(wf)
            # a malformed CSV is rejected with the parser's line-numbered message and
            # the existing file is left untouched
            before = read(wf, String)
            rbad = up("bad.csv", bytes("name;weight\n$(bnames[2]);nope\n"))
            @test occursin("rejected", loc(rbad))
            @test occursin("line", loc(rbad))
            @test read(wf, String) == before
            # non-csv extension, oversized file, and a path-separator name all rejected
            @test occursin("rejected", loc(up("w.txt", UInt8[])))
            @test occursin("rejected", loc(up("w.csv", zeros(UInt8, SparlectraApp.WEBUI_CONTINGENCY_WEIGHTS_MAX_BYTES + 1))))
            @test occursin("rejected", loc(SparlectraApp.handle_contingency_weights_upload(Dict{String,Any}("casefile" => "../evil", "casefiles" => [SparlectraApp.WebUICaseUpload("w.csv", UInt8[])]); output_root=root, case_directory=cases, operation_log=root)))
            # uploading again replaces the file and says so
            @test occursin("replaced", loc(up("w.csv", bytes("name;weight\n$(bnames[1]);2.0\n"))))

            # the editor page seeds the case's real element names plus a raw-CSV editor
            page = redirect_stdout(devnull) do
                SparlectraApp.handle_contingency_weights_page(Dict{String,Any}("case" => "warmup_casePST.m"); output_root=root, case_directory=cases, operation_log=root)
            end
            body = String(page.body)
            @test page.status == 200
            @test occursin(bnames[1], body)
            @test occursin("Raw CSV", body)

            # saving from the seeded table omits rows left at exactly 1.0
            SparlectraApp.handle_contingency_weights_save(Dict{String,Any}("casefile" => "warmup_casePST.m", "element" => [bnames[1], bnames[2]], "weight" => ["2.5", "1.0"]); output_root=root, case_directory=cases, operation_log=root)
            saved = read(wf, String)
            @test occursin(bnames[1], saved)
            @test !occursin(bnames[2], saved)

            # download serves the stored file as an attachment
            dl = SparlectraApp.handle_contingency_weights_download(Dict{String,Any}("case" => "warmup_casePST.m"); output_root=root, case_directory=cases)
            @test dl.status == 200
            @test any(k == "Content-Disposition" for (k, _) in dl.headers)
            @test !isempty(dl.body)

            # reset deletes the weight file
            @test SparlectraApp.handle_contingency_weights_reset(Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root=root, case_directory=cases, operation_log=root).status == 303
            @test !isfile(wf)

            # deleting the case cascades to its weight file
            touch(wf)
            SparlectraApp.handle_powerflow_case_delete(Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root=root, case_directory=cases, operation_log=root)
            @test !isfile(joinpath(cases, "warmup_casePST.m"))
            @test !isfile(wf)
        end)() end

        @testset "state estimation page, measurement upload, and chain (SE phase 5)" begin (function ()
            root = mktempdir()
            cases = joinpath(root, "cases")
            mkpath(cases)
            case_path = joinpath(cases, "warmup_casePST.m")
            cp(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m"), case_path)
            rt = (; case_directory=cases, config_file=Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log=SparlectraApp.webui_operation_log_path(root), startup_config_error=nothing, runner=SparlectraApp.start_powerflow_run)
            bytes = s -> Vector{UInt8}(codeunits(s))

            # upload classification: v1 CSV -> measurement_set; CSV without the
            # version comment -> unknown, retained but never offered
            v1 = "# sparlectra-measurements v1\ntype,bus,from_bus,to_bus,link_nr,direction,value,sigma,active,id\n"
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
            # neither CSV appears in the case selector
            @test !any(endswith(name, ".csv") for name in SparlectraApp._webui_casefile_options_in_directory(cases))
            # traversal names still rejected by the shared upload checks
            bad = SparlectraApp.handle_powerflow_case_import(Dict{String,Any}("casefiles" => [SparlectraApp.WebUICaseUpload("../evil.csv", bytes(v1))]); output_root=root, case_directory=cases, operation_log=root)
            @test bad.status == 303   # redirect with the rejection message
            @test !isfile(joinpath(dirname(cases), "evil.csv"))

            # SE page opens headless, demo generator writes an offered v1 set,
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
            @test occursin(">Runs<", main_page)
            @test !occursin(">Network analysis<", main_page)
            @test !occursin(">New run<", main_page)
            gen = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root=root, runtime=rt)
            @test gen.status == 303
            mfile = joinpath(cases, "warmup_casePST.measurements.csv")
            @test isfile(mfile)
            @test SparlectraApp._webui_is_measurement_csv(mfile)
            page2 = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
            @test occursin("warmup_casePST.measurements.csv", page2)
            @test occursin("Run state estimation", page2)
            @test !occursin("onsubmit", page2)

            # every SE input parameter carries a PF-style help link, and the
            # topics resolve (config-table rows, doc headings, and overrides)
            @test count("help-link", page2) >= 12
            for topic in ("state_estimation.tol", "state_estimation.flatstart", "state_estimation.max_iter", "state_estimation.robust", "state_estimation.update_shunts", "state_estimation.report_residual_correlation", "webui.se_max_eliminations", "webui.se_measurement_file", "webui.se_generator_noise", "webui.se_generator_gross_error", "webui.se_generator_tap_error", "webui.se_generator_sigmas")
                @test SparlectraApp.handle_webui_help(topic).status == 200
            end
            @test occursin("MECHANICAL tap steps", String(SparlectraApp.handle_webui_help("webui.se_generator_tap_error").body))

            # generator options: noise + gross error produce a valid, different set
            plain = read(mfile, String)
            gen2 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "noise" => "true", "gross_error_k" => "8"); output_root=root, runtime=rt)
            @test gen2.status == 303
            @test occursin("bad%20data", Dict(gen2.headers)["Location"])
            @test SparlectraApp._webui_is_measurement_csv(mfile)
            @test read(mfile, String) != plain   # noise + gross error changed values
            # regenerate the clean set for the runs below
            SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root=root, runtime=rt)
            @test read(mfile, String) == plain   # seeded generator is reproducible

            # per-quantity sigmas are PERCENT OF THE MEASURED VALUE with the
            # per-type floors (voltage-level independent); the currents checkbox
            # adds current-magnitude rows. Noise off, so value == truth and the
            # per-row sigma law is exactly reproducible.
            gen3 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "sigma_u_pct" => "1.0", "include_currents" => "true", "sigma_i_pct" => "2.0", "sigma_p_pct" => "2.0", "sigma_q_pct" => "1.0"); output_root=root, runtime=rt)
            @test gen3.status == 303
            netchk = Sparlectra._import_sparlectra_net(case_path, nothing, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true))
            readMeasurementsCSV!(netchk; file=mfile)
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
            netia = Sparlectra._import_sparlectra_net(case_path, nothing, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true))
            readMeasurementsCSV!(netia; file=mfile)
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
            @test read(mfile, String) != plain2
            @test SparlectraApp._webui_is_measurement_csv(mfile)
            # the file records the tap positions the set was generated from as a
            # structured table (electrical/fixed/transferred step columns, the
            # deviated transformer carries its percentage), the SE page renders it
            # as a table and offers the download
            taptxt = read(mfile, String)
            @test occursin("# sparlectra-taps v1", taptxt)
            @test occursin("electrical_step,fixed_step,transferred_step,generation_deviation_steps", taptxt)
            @test occursin(",3.0", taptxt)   # the deviated transformer row
            pageInfo = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
            @test occursin("Transformer taps at generation", pageInfo)
            @test occursin("<th>fixed_step</th>", pageInfo)
            @test occursin("Measured values in this set:", pageInfo)   # per-type row counts
            @test occursin("Vm ×", pageInfo)
            @test occursin("/stateestimation/measurements/download?file=warmup_casePST.measurements.csv", pageInfo)
            dlr = SparlectraApp.route_sparlectra_webui("GET", "/stateestimation/measurements/download?file=warmup_casePST.measurements.csv", Dict{String,String}("file" => "warmup_casePST.measurements.csv"); output_root=root, runtime=rt)
            @test dlr.status == 200
            @test any(k == "Content-Disposition" for (k, _) in dlr.headers)
            dlbad = SparlectraApp.route_sparlectra_webui("GET", "/stateestimation/measurements/download?file=../evil.csv", Dict{String,String}("file" => "../evil.csv"); output_root=root, runtime=rt)
            @test dlbad.status in (400, 404)
            # the commented file still parses and round-trips
            netc = Sparlectra._import_sparlectra_net(case_path, nothing, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true))
            rc = readMeasurementsCSV!(netc; file=mfile)
            @test rc.total > 0
            # out-of-range percentage is rejected
            gentapbad = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "tap_error_steps" => "35"); output_root=root, runtime=rt)
            @test occursin("between%20-16%20and%2016", Dict(gentapbad.headers)["Location"])
            # half steps are no longer settable (a tap changer has no half positions)
            genthalf = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "tap_error_steps" => "1.5"); output_root=root, runtime=rt)
            @test occursin("whole%20number", Dict(genthalf.headers)["Location"])
            # restore the default set for the runs below
            SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root=root, runtime=rt)
            @test read(mfile, String) == plain

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
            SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root=root, runtime=rt)

            # structured value editor: every measurement kind is editable in the
            # table; the update handler rewrites only value/sigma/active and a
            # single invalid entry rejects the whole save
            pageTab = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
            @test occursin("Edit measurement values (table)", pageTab)
            @test occursin("/stateestimation/measurements/update-values", pageTab)
            rowsTab = SparlectraApp._webui_measurement_set_rows(mfile)
            @test !isempty(rowsTab)
            @test any(r -> r.typ == "VmMeas", rowsTab) && any(r -> r.typ == "PflowMeas", rowsTab) && any(r -> r.typ == "PinjMeas", rowsTab)
            vrow = first(r for r in rowsTab if r.typ == "VmMeas")
            rup = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/measurements/update-values", Dict{String,Any}("file" => "warmup_casePST.measurements.csv", "case" => "warmup_casePST.m", "v_$(vrow.line)" => "1.0777", "s_$(vrow.line)" => vrow.sigma, "a_$(vrow.line)" => "true"); output_root=root, runtime=rt)
            @test occursin("updated%201", Dict(rup.headers)["Location"])
            @test occursin("1.0777", read(mfile, String))
            rbadv = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/measurements/update-values", Dict{String,Any}("file" => "warmup_casePST.measurements.csv", "case" => "warmup_casePST.m", "v_$(vrow.line)" => "1.05", "s_$(vrow.line)" => "-1", "a_$(vrow.line)" => "true"); output_root=root, runtime=rt)
            @test occursin("nothing%20was%20saved", Dict(rbadv.headers)["Location"])
            @test occursin("1.0777", read(mfile, String))   # rejected save left the file alone
            SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root=root, runtime=rt)

            # the run summary shows the chi-square band verdict with expected value
            fake = Dict("metadata" => Dict{String,Any}("run_mode" => "se", "se_observability_quality" => "good", "se_iterations" => 4, "se_objective" => 12.5, "se_dof" => 14, "se_band_reason" => "ok", "se_j_within_3sigma" => true))
            stxt = SparlectraApp._webui_se_summary(fake)
            # J/dof leads: a bare J grows with the row count, so the same healthy set
            # reads as an alarm on a larger case (J = 104 at dof 95 is a ratio of 1.1)
            @test occursin("J/dof = 0.89", stxt) && occursin("dof = 14", stxt)
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
            @test occursin("Topology validation (advisory)", tsecT)
            @test occursin("open_element_with_flow", tsecT)
            @test occursin("Suspected topology error", tsecT)
            @test occursin("Test topology hypotheses", tsecT)
            @test occursin("switches NOTHING", tsecT)

            # tap estimation surfaces: run-form checkbox, summary line, result
            # table. The release guards are in, so the label is plain again and
            # the tooltip names the guard behavior instead of an experimental flag
            @test occursin("name=\"se_tap_estimation\"", pageInfo)
            @test !occursin("experimental", pageInfo)
            @test occursin("generator step-up", pageInfo)
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

            # case binding lives IN the file: the generated set records its case,
            # the page shows the binding, the case selector stars bound cases, and
            # a COPY under a foreign name still says which case it belongs to
            @test occursin("# case: warmup_casePST.m", read(mfile, String))
            @test SparlectraApp._webui_measurement_set_case(mfile) == "warmup_casePST.m"
            other = joinpath(cases, "case9.measurements.csv")
            cp(mfile, other)
            page3 = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
            @test occursin("value=\"warmup_casePST.measurements.csv\" selected", page3)
            # block 4 dropped the SE-own case selector (shared selection), and its
            # "*" bound-set marker went with it; the binding stays visible through
            # the Case binding line and the per-set labels asserted here
            @test occursin("Case binding:", page3)
            # under warmup_casePST the copy is bound to warmup_casePST and offered PLAIN (the
            # file decides, not its name); under ANOTHER case it carries the label
            @test occursin(">case9.measurements.csv</option>", page3)
            # suite migration (no-download rule): the foreign case only has to be
            # importable under ANOTHER name; the shipped sp_case5 does that
            case9_path = joinpath(cases, "sp_case5.scf.json")
            cp(joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json"), case9_path)
            page9 = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=sp_case5.scf.json", Dict{String,String}(); output_root=root, runtime=rt).body)
            @test occursin("case9.measurements.csv (for warmup_casePST.m)", page9)
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
            @test occursin("Edit measurement file (inline)", pageEd)
            @test occursin("/stateestimation/measurements/save", pageEd)
            content0 = read(mfile, String)
            edited = replace(content0, "# case: warmup_casePST.m" => "# case: warmup_casePST.m\n# note: edited inline"; count=1)
            rsave = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/measurements/save", Dict{String,Any}("file" => "warmup_casePST.measurements.csv", "case" => "warmup_casePST.m", "content" => edited); output_root=root, runtime=rt)
            @test rsave.status == 303
            @test occursin("saved", Dict(rsave.headers)["Location"])
            @test occursin("# note: edited inline", read(mfile, String))
            rbad1 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/measurements/save", Dict{String,Any}("file" => "../evil.csv", "case" => "warmup_casePST.m", "content" => edited); output_root=root, runtime=rt)
            @test occursin("invalid", Dict(rbad1.headers)["Location"])
            rbad2 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/measurements/save", Dict{String,Any}("file" => "warmup_casePST.measurements.csv", "case" => "warmup_casePST.m", "content" => "a,b\n1,2\n"); output_root=root, runtime=rt)
            @test occursin("rejected", Dict(rbad2.headers)["Location"])
            @test occursin("# note: edited inline", read(mfile, String))   # rejected save left the file alone
            write(mfile, content0)   # restore the clean set for the runs below

            # service-level SE run: artifacts present, history row kind se
            r1 = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
            @test r1["status"] == "succeeded"
            id1 = r1["run_id"]
            for artifact in ("se_state.csv", "se_diagnostics.md", "se_view.md", "measurements.csv")
                @test isfile(joinpath(root, id1, artifact))
            end
            runs = SparlectraApp.list_powerflow_runs(root)
            row = only(r for r in runs if string(get(r, "run_id", "")) == id1)
            @test get(row, "run_mode", "") == "se"
            # the SE result page carries the chain action
            resPage = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(id1)", Dict{String,String}(); output_root=root, runtime=rt).body)
            @test occursin("Run power flow from this estimate", resPage)
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

            # tap-estimation SE run: a set generated with a tap deviation
            # disagrees with the model around one transformer; releasing the taps
            # absorbs the discrepancy and the fixation lands on a mechanical step
            SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "tap_error_steps" => "3"); output_root=root, runtime=rt)
            rtap = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile, "se_tap_estimation" => true))
            @test rtap["status"] == "succeeded"
            idtap = rtap["run_id"]
            tapcsv = joinpath(root, idtap, "se_tap_estimates.csv")
            @test isfile(tapcsv)
            taplines = readlines(tapcsv)
            @test startswith(taplines[1], "branch,name,mrid,mode,alpha_deg,electrical_step,fixed_step")
            @test length(taplines) == 3                     # warmup_casePST: two transformers
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
            # the Delta-u deviation is a SHIFT step, so either step column counts
            header_cols = split(taplines[1], ",")
            step_cols = [findfirst(==(c), header_cols) for c in ("fixed_step", "fixed_shift_step")]
            @test all(i -> i !== nothing, step_cols)
            @test any(l -> (f=split(l, ","); any(parse(Int, f[i]) != 0 for i in step_cols)), taplines[2:end])
            tapPage = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(idtap)", Dict{String,String}(); output_root=root, runtime=rt).body)
            @test occursin("Transformer tap estimates", tapPage)
            @test occursin("tap estimation: 1 transformer(s)", tapPage)

            # A set that DOCUMENTS its tap deviations releases exactly those
            # transformers by itself, without the user asking for it. As a mere
            # hint this produced a J nobody could explain, and the hint only
            # appeared after a converged run (if transformers were changed, that option has to
            # be on by itself).
            rhint = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
            @test rhint["status"] == "succeeded"
            @test rhint["metadata"]["se_set_tap_deviation"] == true
            @test !isempty(rhint["metadata"]["se_auto_released_taps"])
            @test occursin("released automatically", rhint["message"])
            @test occursin("tap estimation released automatically", read(joinpath(root, rhint["run_id"], "run.log"), String))

            # bad data lands findable: gross error in the set -> se_bad_data.csv
            # with the measurement, its location, and the eliminated flag
            SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gross_error_k" => "10"); output_root=root, runtime=rt)
            rbd = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
            @test rbd["status"] == "succeeded"
            bdcsv = joinpath(root, rbd["run_id"], "se_bad_data.csv")
            @test isfile(bdcsv)
            bdtxt = read(bdcsv, String)
            @test startswith(bdtxt, "measurement_index,id,type,bus,from_bus,to_bus,normalized_residual,wii,localizable,eliminated,suppressed,downweighted")
            @test occursin(r"[PQ](flow|inj)_", bdtxt) # the corrupted power row is listed (P or Q, the seed decides)
            @test occursin("true", bdtxt)            # localizable/eliminated flags present
            @test occursin("se_bad_data.csv", read(joinpath(root, rbd["run_id"], "run.log"), String))
            # bad data is one click away on the result page (the direct download)
            bdPage = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(rbd["run_id"])", Dict{String,String}(); output_root=root, runtime=rt).body)
            @test occursin("Bad data (se_bad_data.csv)", bdPage)
            @test occursin("/powerflow/artifact/$(rbd["run_id"])/se_bad_data.csv", bdPage)

            # restore the clean default set for anything below
            SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root=root, runtime=rt)

            # --- measurement generator v2: truth state, flow ends, passive nodes,
            # delta comments, and the bad-data threshold surface of the run form
            genpage = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root=root, runtime=rt).body)
            for needle in ("name=\"gen_truth_source\"", "name=\"gen_truth_run_id\"", "name=\"gen_flow_ends\"", "name=\"gen_passive_sigma\"", "name=\"gen_passive_as_zi\"", "name=\"se_robust_mode\"", "name=\"se_k_eliminate\"", "name=\"se_k_suppress\"", "name=\"se_suppression_sigma\"", "se-threshold-warning", "gen-truth-source")
                @test occursin(needle, genpage)
            end

            # one balance-aware flow end: deterministic (two generates produce the
            # identical file), one flow group per branch, choice documented
            SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_flow_ends" => "one_balance_aware"); output_root=root, runtime=rt)
            one1 = read(mfile, String)
            SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_flow_ends" => "one_balance_aware"); output_root=root, runtime=rt)
            one2 = read(mfile, String)
            @test one1 == one2
            @test occursin("# seed: 42", one1)
            # generating persists the generator options in the form block of the
            # case configuration file so a case reload restores them
            @test occursin("gen_flow_ends: one_balance_aware", read(joinpath(cases, "warmup_casePST.m.config.yaml"), String))
            @test occursin("# flow_ends: one_balance_aware", one1)
            @test occursin("# flow_end,", one1)
            @test occursin("# truth_value,Vm_", one1)
            # critical measurements on request: the set is thinned until two rows
            # are critical, stays observable, and names the rows in its comments
            SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_critical_count" => "2"); output_root=root, runtime=rt)
            crit_text = read(mfile, String)
            @test occursin("# critical_target: 2 reached: ", crit_text)
            @test occursin("# critical_rows: ", crit_text)
            crit_net = Sparlectra.import_case(joinpath(cases, "warmup_casePST.m"), Sparlectra.load_sparlectra_config()).net
            Sparlectra.readMeasurementsCSV!(crit_net; file=mfile)
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
            SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_passive_as_zi" => "true"); output_root=root, runtime=rt)
            zi1 = read(mfile, String)
            @test occursin("ZI_PINJ_", zi1)
            zilines = [l for l in split(zi1, "\n") if startswith(l, "PinjMeas") && occursin("ZI_PINJ_bus_", l)]
            @test !isempty(zilines)
            ziline = first(zilines)
            zibus = split(ziline, ",")[2]
            @test !any(l -> startswith(l, "PinjMeas,$(zibus),") && !occursin("ZI_", l), split(zi1, "\n"))
            rzi = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
            @test rzi["status"] == "succeeded"
            @test rzi["metadata"]["se_eliminations"] == 0

            # truth state from run: adopts the SE run's state (bit-exact against
            # se_state.csv), documents the source, and the delta file appears
            SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_truth_source" => "from_run", "gen_truth_run_id" => id1); output_root=root, runtime=rt)
            fr1 = read(mfile, String)
            @test occursin("# truth: run $(id1) (se,", fr1)
            sestate = readlines(joinpath(root, id1, "se_state.csv"))
            strow = only([l for l in sestate if startswith(l, "2,")])
            st_vm = split(strow, ",")[2]
            vmrow = only([l for l in split(fr1, "\n") if startswith(l, "VmMeas,2,")])
            @test split(vmrow, ",")[8] == st_vm
            # rejections: tap deviation locked, unknown run id named
            rrej = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_truth_source" => "from_run", "gen_truth_run_id" => id1, "tap_error_steps" => "2"); output_root=root, runtime=rt)
            @test occursin("requires%20truth%20state", Dict(rrej.headers)["Location"])
            rrej2 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_truth_source" => "from_run", "gen_truth_run_id" => "nope"); output_root=root, runtime=rt)
            @test occursin("not%20found%20in%20the%20run%20history", Dict(rrej2.headers)["Location"])

            # restore the clean default set and exercise the threshold surface:
            # the staged service path equals the legacy Bool bitwise, the delta
            # artifact exists for generated sets, invalid modes reject
            SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root=root, runtime=rt)
            rst = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile, "se_robust_mode" => "staged", "se_robust_k1" => 3.0, "se_robust_k2" => 6.0))
            @test rst["status"] == "succeeded"
            @test rst["metadata"]["se_robust_mode"] == "staged"
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
            @test occursin("kind,id,type,bus,from_bus,to_bus,measured,truth,estimated", read(joinpath(root, rst["run_id"], "se_deltas.csv"), String))
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
            # transformers, so a max of 2 actually draws 2 (derived 2026-09-04)
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
            SparlectraApp._webui_merge_case_settings!(root, case_path, Dict{String,Any}("gen_seed" => 99); case_directory=dirname(case_path))
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
            rpst = start_powerflow_run(Dict{String,Any}("casefile" => pstcase, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => pstout, "se_tap_estimation" => true))
            @test rpst["status"] == "succeeded"
            @test rpst["metadata"]["se_band_reason"] == "ok"
            pstrow = only(t for t in rpst["metadata"]["se_tap_estimates"] if t["branch"] == 8)
            @test pstrow["mode"] == "pst"
            @test pstrow["fixed_shift_step"] == 2
            @test abs(pstrow["electrical_shift_step"] - 2.0) < 0.3
            gsurow = only(t for t in rpst["metadata"]["se_tap_estimates"] if t["branch"] == 7)
            @test gsurow["source"] == "calculated"

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
            devtxt = read(joinpath(r2["output_dir"], "se_pf_deviation.csv"), String)
            @test startswith(devtxt, "bus,name,mrid,vm_se_pu,vm_pf_pu,dvm_pu,va_se_deg,va_pf_deg,dva_deg")

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
            for (l3, l4) in zip(lines3, lines4)
                f3 = split(l3, csv_delim)
                f4 = split(l4, csv_delim)
                @test length(f3) == length(f4)
                for (a, b) in zip(f3, f4)
                    na = tryparse(Float64, a)
                    nb = tryparse(Float64, b)
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

        @testset "saved case settings outrank the configuration on the form" begin (function ()
            # No mtime logic anywhere: the case levels always win over the
            # general configuration for the fields they set, however new the
            # configuration file is.
            dir = mktempdir()
            cfg = joinpath(dir, "conf.yaml")
            write(cfg, "config_version: 1\nscope: general\npower_flow:\n  max_iter: 99\n")
            profile = Dict{String,Any}("power_flow_max_iter" => 55, "power_flow_autodamp_min" => 0.07, "_profile_path" => joinpath(dir, "case57.config.yaml"))
            v = SparlectraApp.webui_form_state(selected_config_file=cfg, sidecar_profile=profile)
            @test v["power_flow_max_iter"] == 55
            @test v["power_flow_autodamp_min"] == 0.07
            @test !haskey(v, "_config_newer_than_profile")
            # a field the case does not set comes from the configuration
            v_cfg_only = SparlectraApp.webui_form_state(selected_config_file=cfg)
            @test v_cfg_only["power_flow_max_iter"] == 99
        end)() end

        @testset "state estimation form settings are case scope (issue #377)" begin (function ()
            dir = mktempdir()
            casedir = joinpath(dir, "cases")
            mkpath(casedir)
            root = joinpath(dir, "out")
            mkpath(root)
            cp(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")), joinpath(casedir, "sp_case14.scf.json"))
            cp(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.measurements.csv")), joinpath(casedir, "sp_case14.measurements.csv"))
            cp(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case60.scf.json")), joinpath(casedir, "sp_case60.scf.json"))
            rt = (; case_directory=casedir, config_file=Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log=SparlectraApp.webui_operation_log_path(root), startup_config_error=nothing, runner=SparlectraApp.start_powerflow_run)

            # test: configuration.yaml with k_eliminate 3.5 and max_eliminations 0,
            # page opened, fields show the file's values without user input, and a
            # run without form input (case sidecar with max_eliminations: 0) reports
            # se_eliminations = 0 - the issue's own regression: it eliminated a row
            # although the file said 0.
            config_path = joinpath(dir, "config.yaml")
            cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, config_path)
            open(config_path, "a") do io
                println(io, "state_estimation:")
                println(io, "  k_eliminate: 3.5")
                println(io, "  max_eliminations: 0")
            end
            rt_yaml = (; case_directory=casedir, config_file=config_path, operation_log=SparlectraApp.webui_operation_log_path(root), startup_config_error=nothing, runner=SparlectraApp.start_powerflow_run)
            page = String(copy(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=sp_case14.scf.json"; output_root=root, runtime=rt_yaml).body))
            @test occursin("id=\"se-k-eliminate\" value=\"3.5\"", page)
            @test occursin("name=\"se_max_eliminations\" value=\"0\"", page)
            run_yaml = SparlectraApp.start_powerflow_run(Dict("casefile" => "sp_case14.scf.json", "config_file" => config_path, "output_root" => root, "se_mode" => true); case_directory=casedir)
            @test run_yaml["status"] == "succeeded"
            @test run_yaml["metadata"]["se_eliminations"] == 0
            @test run_yaml["metadata"]["se_k_eliminate"] == 3.5

            # submitting the SAME run form unchanged (as the browser would with the
            # page correctly pre-filled above) must not silently outrank the file
            submit_form = Dict{String,Any}(
                "casefile" => "sp_case14.scf.json", "config_file" => config_path, "se_mode" => "true",
                "measurement_file" => "sp_case14.measurements.csv", "se_flatstart" => "true", "se_tol" => "1e-6",
                "se_max_iter" => "50", "se_robust_mode" => "off", "se_k_eliminate" => "3.5", "se_robust_k1" => "3.0",
                "se_robust_k2" => "6.0", "se_k_suppress" => "4.0", "se_suppression_sigma" => "2000", "se_max_eliminations" => "0",
            )
            resp = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/run", submit_form; output_root=root, runtime=rt_yaml)
            run_id = basename(only(h.second for h in resp.headers if h.first == "Location"))
            wait(SparlectraApp._POWERFLOW_WEBUI_JOBS[run_id]["task"])
            submitted_result = get_powerflow_result(run_id)
            @test submitted_result["metadata"]["se_eliminations"] == 0
            @test submitted_result["metadata"]["se_k_eliminate"] == 3.5

            # test: form value changed and Save settings pressed -> value appears
            # in the case sidecar; next run without form input uses it
            save_form = copy(submit_form)
            save_form["se_k_eliminate"] = "5.0"
            save_form["settings_target"] = "this_case"
            save_resp = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/settings/save", save_form; output_root=root, runtime=rt)
            @test save_resp.status == 303
            sidecar = read(joinpath(casedir, "sp_case14.config.yaml"), String)
            @test occursin("k_eliminate: 5.0", sidecar)
            run_after_save = SparlectraApp.start_powerflow_run(Dict("casefile" => "sp_case14.scf.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true); case_directory=casedir)
            @test run_after_save["metadata"]["se_k_eliminate"] == 5.0

            # test: switching cases restores each case's own values - sp_case60 has
            # no sidecar, so it must show the STRUCT default (3.0), not case14's 5.0
            page14_after = String(copy(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=sp_case14.scf.json"; output_root=root, runtime=rt).body))
            @test occursin("id=\"se-k-eliminate\" value=\"5.0\"", page14_after)
            page60 = String(copy(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=sp_case60.scf.json"; output_root=root, runtime=rt).body))
            @test occursin("id=\"se-k-eliminate\" value=\"3.0\"", page60)

            # regression (already true, kept per the issue): power_flow.solver in
            # the case sidecar reaches the run without touching Settings. The
            # shipped sp_case* demos all carry active Q(U)/P(U) controllers,
            # which apslf refuses outright (a real, unrelated constraint), so the
            # check is that the sidecar's solver choice reaches the run at all
            # (visible in effective_config.yaml / the failure naming apslf by
            # name), not that apslf converges on one of these networks.
            write(joinpath(casedir, "sp_case14.config.yaml"), "config_version: 1\nscope: case\ncase: sp_case14.scf.json\npower_flow:\n  solver: apslf\n")
            run_solver = SparlectraApp.start_powerflow_run(Dict("casefile" => "sp_case14.scf.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root); case_directory=casedir)
            @test occursin("power_flow.solver=apslf", string(get(run_solver, "message", "")))
        end)() end

        @testset "browser opening falls back to the system default" begin (function ()
            # Reported 2026-09-07 from Windows 11 with Edge uninstalled and only
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

                @test sl(:sysimage_problem, img, proj) == "no sysimage found"
                write(img, "not a real image")
                # image without metadata is unusable: the contract is the pair
                @test sl(:sysimage_problem, img, proj) == "no sysimage found"
                write(meta, "julia_version = \"9.9.9\"\nmanifest_sha256 = \"$(sha)\"\n")
                @test occursin("built for Julia 9.9.9", sl(:sysimage_problem, img, proj))
                write(meta, string("julia_version = \"", VERSION, "\"\nmanifest_sha256 = \"deadbeef\"\n"))
                @test sl(:sysimage_problem, img, proj) == "the sysimage does not match the current Manifest.toml"
                write(meta, full_meta)
                touch(img)
                @test sl(:sysimage_problem, img, proj) === nothing
                # a source file newer than the image disables it (mtime resolution
                # can be one second, hence the wait)
                sleep(1.1)
                touch(joinpath(proj, "src", "Fixture.jl"))
                @test occursin("older than", sl(:sysimage_problem, img, proj))
                write(meta, "kaputt = [[[")
                @test sl(:sysimage_problem, img, proj) == "the sysimage metadata is unreadable"

                # the decision itself, through its seams: a fixed answer, a recorded
                # relaunch, the image this process runs on
                slk(name, args...; kw...) = Base.invokelatest(Base.invokelatest(getfield, SL, name), args...; kw...)
                decide(args; ask, current = "") = begin
                    relaunches = String[]
                    # stdout can be redirected into a file, not into a buffer
                    capture = joinpath(tmp, "decision.txt")
                    open(capture, "w") do out
                        redirect_stdout(out) do
                            slk(:handle_sysimage, args, "start_webui.jl", proj; image = img, ask = (_...) -> ask, relaunch = (image, _...) -> push!(relaunches, String(image)), current_image = current)
                        end
                    end
                    (text = read(capture, String), relaunches = relaunches)
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
                @test slk(:external_sysimage, img; current = "") === nothing
                @test slk(:external_sysimage, img; current = joinpath(Sys.BINDIR, "..", "lib", "julia", "sys.so")) === nothing
                @test slk(:remove_stale_sysimage, joinpath(tmp, "nowhere.so")) === :absent
            end
            # no terminal in the test process, so the question answers itself with
            # no: an unattended start goes without an image, a build is opt-in
            @test sl(:ask_build, "never shown") == false
            @test slval(:REBUILD_FLAG) == "--rebuild-sysimage"
            @test slval(:NO_IMAGE_FLAG) == "--no-sysimage"
            # platform path contract, still owned by the package for everything
            # that runs after the load
            @test endswith(sl(:sysimage_path), SparlectraApp.webui_sysimage_ext())
            @test sl(:sysimage_path) == SparlectraApp.webui_sysimage_path()
        end)() end

        @testset "the Info menu stays in the header on every page" begin (function ()
            # Regression (2026-09-07): the Info control was passed to the layout by
            # the Case, Settings and Runs renderers only, so it disappeared as soon
            # as the user clicked Operation Log, Run history, Last errors, Docs or a
            # result page. Everything it shows describes the running server, not the
            # page, so the header builds it itself now.
            root = SparlectraApp.default_webui_output_root()
            for path in ("/powerflow", "/powerflow/case", "/powerflow/settings", "/powerflow/history",
                "/webui/operation-log", "/webui/last-errors", "/docs", "/webui/sysimage")
                response = SparlectraApp.route_sparlectra_webui("GET", path; output_root=root)
                body = String(response.body)
                # the tuple form names the offending page in the failure output
                @test (path, occursin("topbar-info-menu", body)) == (path, true)
                # the box names the Julia the server runs on: with
                # 1.12 and 1.13 both in use, a screenshot has to say which
                @test occursin("<dt>Julia</dt><dd><code>$(VERSION)</code></dd>", body)
                # and the AnalyticLoadFlow it solves with, so an outdated
                # environment is visible on screen
                @test occursin("<dt>AnalyticLoadFlow</dt><dd><code>v$(pkgversion(Sparlectra.AnalyticLoadFlow))</code> <a href=\"https://github.com/Welthulk/AnalyticLoadFlow.jl\"", body)
            end
            # error pages carry it as well: they are where a user looks for the
            # output root and the operation log in the first place
            @test occursin("topbar-info-menu", SparlectraApp.render_webui_error(404, "not found"))
        end)() end

        # Two runs of the same case under different settings are the normal way to
        # look at a Q-limit mode or a solver choice. The history could list them but
        # not put them side by side, so the comparison had to happen in the head or
        # in two browser tabs.
        @testset "run comparison" begin (function ()
            @testset "history offers a selection and a compare button" begin (function ()
                root = SparlectraApp.default_webui_output_root()
                body = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/history"; output_root=root).body)
                @test occursin("action=\"/powerflow/compare\"", body)
                @test occursin("Compare selected runs", body)
                @test occursin("<th>Compare</th>", body)
            end)() end

            # Maintainer 2026-09-09: the box should appear only where the two runs
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
                # the paths no longer fit the row: name in the cell, path in the tooltip
                @test occursin("<td title=\"/very/long/path/to/cases/case118.m\">case118.m</td>", html)
                @test occursin("<td title=\"/etc/sparlectra/configuration.yaml\">configuration.yaml</td>", html)
                @test !occursin("<td>/very/long/path/to/cases/case118.m</td>", html)
            end)() end

            # The result page names the case and the phase up top, where a reader
            # looks first; the rows that carried the path twice and the
            # final_outcome row that said nothing are gone.
            @testset "result page: case and phase up top, path rows gone" begin (function ()
                probe = Dict{String,Any}("run_id" => "r-1", "status" => "succeeded", "success" => true, "converged" => true,
                    "casefile" => "/some/where/sp_case14.scf.json", "resolved_casefile" => "/some/where/sp_case14.scf.json",
                    "current_phase" => "finished", "final_outcome" => Dict{String,Any}("solver" => "rectangular"), "artifacts" => Any[])
                html = SparlectraApp.render_powerflow_result(probe)
                @test occursin("<span class=\"summary-label\">Case</span><code title=\"/some/where/sp_case14.scf.json\">sp_case14.scf.json</code>", html)
                @test occursin("<span class=\"summary-label\">Phase</span><code>finished</code>", html)
                for gone in ("<th>casefile</th>", "<th>resolved_casefile</th>", "<th>current_phase</th>", "<th>final_outcome</th>")
                    @test !occursin(gone, html)
                end
                # an active run carries the same two cards
                active = merge(probe, Dict{String,Any}("status" => "running", "current_phase" => "linear_solve"))
                active_html = SparlectraApp.render_powerflow_result(active)
                @test occursin("<span class=\"summary-label\">Phase</span><code>linear_solve</code>", active_html)
                @test occursin(">sp_case14.scf.json</code>", active_html)
                # A live snapshot carries `nothing` where the run has not produced
                # the value yet; the card once printed that word for the whole run
                # and named the case only at the end.
                live = Dict{String,Any}("run_id" => "r-2", "status" => "running", "success" => false,
                    "casefile" => "/some/where/case14.m", "resolved_casefile" => nothing, "current_phase" => nothing, "last_phase" => nothing, "artifacts" => Any[])
                live_html = SparlectraApp.render_powerflow_result(live)
                @test occursin("<span class=\"summary-label\">Case</span><code title=\"/some/where/case14.m\">case14.m</code>", live_html)
                @test occursin("<span class=\"summary-label\">Phase</span><code>n/a</code>", live_html)
                @test !occursin(">nothing<", live_html)

                # a Q(U) machine read from a case file is named on the page, not only
                # in the Control column of the result print (task qu_scf, 2026-09-11);
                # a run without controllers keeps its summary short
                @test !occursin("summary-label\">Controllers<", html)
                with_ctrl = merge(probe, Dict{String,Any}("metadata" => Dict{String,Any}("controllers" => Dict{String,Any}("tap" => 1, "qu" => 1, "pu" => 0))))
                @test occursin("<span class=\"summary-label\">Controllers</span><code>Q(U) 1 · tap 1</code>", SparlectraApp.render_powerflow_result(with_ctrl))
                @test SparlectraApp._webui_control_summary(Dict{String,Any}("metadata" => Dict{String,Any}("controllers" => Dict{String,Any}("tap" => 0, "qu" => 0, "pu" => 0)))) === nothing
            end)() end

            # A native file input speaks the BROWSER's language ("Durchsuchen",
            # "Keine Datei ausgewählt" on a German browser) on an English page.
            # The input stays in the form, off screen; a
            # label is the button and a span names the selection.
            @testset "file pickers speak the page's language" begin (function ()
                one = SparlectraApp._webui_file_input("weights_file"; accept=".csv", required=true)
                many = SparlectraApp._webui_file_input("casefiles"; accept=".m,.zip", multiple=true)
                alias = SparlectraApp._webui_file_input("casefiles"; accept=".csv", id="measurements")
                @test occursin("type=\"file\" name=\"weights_file\" accept=\".csv\" required>", one)
                @test occursin("<label for=\"file-field-weights_file\" class=\"button secondary-button file-field-button\">Choose file</label>", one)
                @test occursin("No file selected", one)
                @test occursin("multiple>", many)
                @test occursin(">Choose files</label>", many)
                # two pickers with one field name on one page need distinct ids
                @test occursin("id=\"file-field-measurements\"", alias)
                @test occursin("name=\"casefiles\"", alias)
                # the script that keeps the span current is part of every page
                page = SparlectraApp._webui_layout("t", "<p>x</p>")
                @test occursin("data-file-field-input", page)
                @test occursin("'No file selected'", page)
            end)() end

            @testset "the route insists on exactly two runs" begin (function ()
                root = SparlectraApp.default_webui_output_root()
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
                write(joinpath(dir_a, "effective_config.yaml"), "power_flow:\n  tol: 1.0e-8\n  qlimits:\n    enforcement_mode: active_set\n")
                write(joinpath(dir_b, "effective_config.yaml"), "power_flow:\n  tol: 1.0e-8\n  qlimits:\n    enforcement_mode: classic_simultaneous\n")
                diff = SparlectraApp._webui_compare_config_diff(dir_a, dir_b)
                # the nested key keeps its real path; a "remember the last section"
                # reader concatenated every section it had ever seen
                @test diff == [("power_flow.qlimits.enforcement_mode", "active_set", "classic_simultaneous")]

                write(joinpath(dir_a, "q_limit_events.csv"), "iteration,bus,side\n2,19,min\n2,32,min\n")
                write(joinpath(dir_b, "q_limit_events.csv"), "iteration,bus,side\n1,19,min\n")
                @test SparlectraApp._webui_compare_qlimit_buses(dir_a) == [19, 32]
                @test SparlectraApp._webui_compare_qlimit_buses(dir_b) == [19]
                @test SparlectraApp._webui_compare_qlimit_buses(mktempdir()) === nothing

                # The detailed export writes one of three formats, and A carries the
                # German one: ';' delimiter, decimal comma. A comma-splitting reader
                # found no columns there and the page then claimed the two runs had
                # nothing in common ("The two runs share no bus names", reported
                # 2026-09-08). The fixtures below are the real header of
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
                @test occursin("Configuration differences (1)", html)
                @test occursin("power_flow.qlimits.enforcement_mode", html)
                @test occursin("The runs clamped different buses", html)
                @test occursin("max |dVm|", html)
                @test occursin("0.01", html)          # 0.99 gegen 0.98
                @test occursin("NEWBERRY 2", html)
                @test occursin("Losses", html)
                @test occursin("0.1 MW", html) || occursin("0.09999", html)  # 0.75 gegen 0.65
                @test occursin("/powerflow/result/run-a", html)

                # two identical runs say so in one line instead of a table of zeros
                same = SparlectraApp.render_powerflow_compare(a, a)
                @test occursin("Both runs end on the same voltages", same)
                @test occursin("Both runs end on the same losses", same)
            end)() end
        end)() end

        # Reported from a live session 2026-09-08: "Case input format" offered
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

            # Seen on Windows 2026-09-24: `scf` stored for a MATPOWER case sent
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

            # The same word saved from the settings page (2026-09-11: "Q-Limit
            # an, aber mode auf off" ended in an ArgumentError on the next run):
            # the case file and the configuration file get `enabled: false` and
            # no mode key, and a file that already carries `off` still loads.
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
            # a file written by 0.11.1 with the word inside runs as "disabled"
            stale = joinpath(root, "stale.yaml")
            write(stale, "power_flow:\n  qlimits:\n    enforcement_mode: off\n")
            stale_cfg = Sparlectra.load_sparlectra_config(stale; reload=true)
            @test stale_cfg.powerflow.qlimits.ignore_q_limits
            @test stale_cfg.powerflow.qlimits.enforcement_mode === :active_set
            stale_run = SparlectraApp.start_powerflow_run(Dict{String,Any}("casefile" => "sp_case14.scf.json", "config_file" => stale, "output_root" => joinpath(root, "stale_runs"), "config_overrides" => Dict{String,Any}("power_flow.qlimits.enforcement_mode" => "off")); case_directory=cache)
            @test stale_run["status"] == "succeeded"
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

            # The tooltip carries the raw status only where the label says something
            # else. Putting it on every badge changed the markup of every ordinary
            # run and broke a live-page assertion in the extended profile
            # (test_webui_extended.jl, 2026-09-08); repeating "running" on hover
            # tells a reader nothing anyway.
            @test SparlectraApp._webui_status_badge("status-running", "running", "running") ==
                "<span class=\"status-badge status-running\">running</span>"
            @test SparlectraApp._webui_status_badge("status-info", "diagnosed", "not_converged") ==
                "<span class=\"status-badge status-info\" title=\"not_converged\">diagnosed</span>"
        end)() end

        @testset "the environment is checked before the sysimage question" begin (function ()
            # Regression 2026-09-07 (Windows 11). A checkout carried a Manifest.toml
            # from before AnalyticLoadFlow became a required dependency. The start
            # asked whether to build a sysimage FIRST, the user declined, and only
            # then did `using Sparlectra` fail with a KeyError deep in
            # Base.Precompilation. Answering yes would have been worse: the build
            # works in the separate @sparlectra-sysimage-build environment and would
            # have spent eleven minutes before the checkout's own manifest turned
            # out to be the problem.
            #
            # unresolved_dependencies answers that from two TOML reads, with no
            # package load, which is what lets it run before the question.
            launcher = Module(:LauncherEnvCheck)
            Base.include(launcher, joinpath(Sparlectra.SPARLECTRA_ROOT, "tools", "sysimage_launcher.jl"))
            # Base.invokelatest also for the binding: Julia 1.12+ warns (and a later
            # version errors) on reading a global defined in this same top-level
            # expression by the include above
            SL = Base.invokelatest(getfield, launcher, :SysimageLauncher)
            unresolved(dir) = Base.invokelatest(Base.invokelatest(getfield, SL, :unresolved_dependencies), dir)

            # this very checkout is resolvable, so the check stays silent
            @test isempty(unresolved(Sparlectra.SPARLECTRA_ROOT))

            mktempdir() do dir
                json = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
                toml = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
                # a manifest that predates a dependency: exactly the reported case
                stale = joinpath(dir, "stale")
                mkpath(stale)
                write(joinpath(stale, "Project.toml"), "[deps]\nJSON = \"$(json)\"\nTOML = \"$(toml)\"\n")
                write(joinpath(stale, "Manifest.toml"), "julia_version = \"$(VERSION)\"\nmanifest_format = \"2.0\"\n\n[deps]\n[[deps.TOML]]\nuuid = \"$(toml)\"\n")
                @test unresolved(stale) == ["JSON"]

                # a fresh checkout has no manifest at all
                fresh = joinpath(dir, "fresh")
                mkpath(fresh)
                write(joinpath(fresh, "Project.toml"), "[deps]\nJSON = \"$(json)\"\n")
                @test unresolved(fresh) == ["<no Manifest.toml>"]

                # a complete manifest is silent
                good = joinpath(dir, "good")
                mkpath(good)
                write(joinpath(good, "Project.toml"), "[deps]\nTOML = \"$(toml)\"\n")
                write(joinpath(good, "Manifest.toml"), "julia_version = \"$(VERSION)\"\nmanifest_format = \"2.0\"\n\n[deps]\n[[deps.TOML]]\nuuid = \"$(toml)\"\n")
                @test isempty(unresolved(good))

                # a manifest version below the compat bound: AnalyticLoadFlow
                # 0.9.14 against compat "0.9.15"
                alf = "19ecf91d-4042-490c-adc1-7fb0fdb392c0"
                old_version = joinpath(dir, "old_version")
                mkpath(old_version)
                write(joinpath(old_version, "Project.toml"), "[deps]\nAnalyticLoadFlow = \"$(alf)\"\nTOML = \"$(toml)\"\n\n[compat]\nAnalyticLoadFlow = \"0.9.15\"\n")
                write(joinpath(old_version, "Manifest.toml"), "julia_version = \"$(VERSION)\"\nmanifest_format = \"2.0\"\n\n[deps]\n[[deps.AnalyticLoadFlow]]\nuuid = \"$(alf)\"\nversion = \"0.9.14\"\n\n[[deps.TOML]]\nuuid = \"$(toml)\"\n")
                @test unresolved(old_version) == ["AnalyticLoadFlow 0.9.14 < 0.9.15"]
                write(joinpath(old_version, "Manifest.toml"), "julia_version = \"$(VERSION)\"\nmanifest_format = \"2.0\"\n\n[deps]\n[[deps.AnalyticLoadFlow]]\nuuid = \"$(alf)\"\nversion = \"0.9.15\"\n\n[[deps.TOML]]\nuuid = \"$(toml)\"\n")
                @test isempty(unresolved(old_version))
                bound(spec) = Base.invokelatest(Base.invokelatest(getfield, SL, :compat_lower_bound), spec)
                @test bound("0.9.15") == v"0.9.15"
                @test bound("^0.9.15") == v"0.9.15"
                @test bound("~1.2") == v"1.2.0"
                @test bound("=0.9.15") == v"0.9.15"
                @test bound("0.9.15, 0.10") == v"0.9.15"
                @test bound("") === nothing
                @test bound(">= 0.9") === nothing
                @test bound("0.9 - 0.10") === nothing

                # an unreadable manifest counts as unresolved, not as fine: the
                # environment cannot be trusted either way
                broken = joinpath(dir, "broken")
                mkpath(broken)
                write(joinpath(broken, "Project.toml"), "[deps]\n")
                write(joinpath(broken, "Manifest.toml"), "kaputt = [[[")
                broken_result = redirect_stdout(devnull) do
                    unresolved(broken)
                end
                @test length(broken_result) == 1
                @test occursin("unreadable", broken_result[1])
            end
        end)() end

        @testset "the repair resolves before it instantiates" begin (function ()
            # Regression 2026-09-07 (Windows 11). A checkout carried a Manifest.toml
            # from before AnalyticLoadFlow became a required dependency. `using
            # Sparlectra` failed with a KeyError deep in Base.Precompilation, the
            # recovery block caught it and ran Pkg.instantiate, and instantiate
            # cannot add a package the manifest never mentioned: it failed with
            # "AnalyticLoadFlow is a direct dependency, but does not appear in the
            # manifest ... run Pkg.resolve()". One unhelpful error became another.
            #
            # Verified in a throwaway environment: instantiate alone reproduces that
            # message, resolve followed by instantiate succeeds. The order is what
            # this test guards, and a text check is the honest tool for it: the
            # recovery lives at the top level of a script, which cannot be called
            # without starting a Web UI.
            script = read(joinpath(Sparlectra.SPARLECTRA_ROOT, "tools", "sysimage_launcher.jl"), String)
            resolve_at = findfirst("pkgm.resolve", script)
            instantiate_at = findfirst("pkgm.instantiate", script)
            @test resolve_at !== nothing
            @test instantiate_at !== nothing
            @test first(resolve_at) < first(instantiate_at)
            # and the advice printed on failure has to name the same two calls, or a
            # user who runs it by hand hits exactly the error we just fixed
            @test occursin("Pkg.resolve(); Pkg.instantiate()", script)
            # a manifest version below the compat bound is updated BEFORE the
            # resolve: resolve keeps manifest versions as explicit requirements
            # and fails as unsatisfiable on such a manifest
            update_at = findfirst("pkgm.update", script)
            @test update_at !== nothing
            @test first(update_at) < first(resolve_at)
            launcher = Module(:LauncherRepairCheck)
            Base.include(launcher, joinpath(Sparlectra.SPARLECTRA_ROOT, "tools", "sysimage_launcher.jl"))
            SL = Base.invokelatest(getfield, launcher, :SysimageLauncher)
            outdated(list) = Base.invokelatest(Base.invokelatest(getfield, SL, :outdated_dependencies), list)
            @test outdated(["AnalyticLoadFlow 0.9.14 < 0.9.15", "JSON", "<no Manifest.toml>"]) == ["AnalyticLoadFlow"]
            @test isempty(outdated(String[]))
        end)() end

        @testset "sysimage validity: launcher and package agree" begin (function ()
            # SparlectraApp.webui_sysimage_problem answers the same question for the
            # Web UI's sysimage page that SysimageLauncher.sysimage_problem answers
            # for the start decision. The launcher cannot call the package one (it
            # runs before the package is loaded), so the two implementations are
            # held together HERE: every fixture state below has to produce the same
            # verdict on both sides, or the sysimage page will one day claim an
            # image is fine that the launcher refuses to start from.
            launcher = Module(:LauncherParity)
            Base.include(launcher, joinpath(Sparlectra.SPARLECTRA_ROOT, "tools", "sysimage_launcher.jl"))
            # Base.invokelatest also for the binding: Julia 1.12+ warns (and a later
            # version errors) on reading a global defined in this same top-level
            # expression by the include above
            SL = Base.invokelatest(getfield, launcher, :SysimageLauncher)
            sl(name, args...) = Base.invokelatest(Base.invokelatest(getfield, SL, name), args...)
            both(img, proj) = (sl(:sysimage_problem, img, proj), SparlectraApp.webui_sysimage_problem(image_path=img, project_dir=proj))
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
                states = Pair{String,Function}[
                    "missing image"=>()->nothing,
                    "image without metadata"=>()->write(img, "not a real image"),
                    "wrong julia version"=>()->write(meta, "julia_version = \"9.9.9\"\nmanifest_sha256 = \"$(sha)\"\n"),
                    "manifest mismatch"=>()->write(meta, string("julia_version = \"", VERSION, "\"\nmanifest_sha256 = \"deadbeef\"\n")),
                    "valid"=>()->(write(meta, string("julia_version = \"", VERSION, "\"\nmanifest_sha256 = \"", sha, "\"\n")); touch(img)),
                    "source newer than image"=>()->(sleep(1.1); touch(joinpath(proj, "src", "Fixture.jl"))),
                    "unreadable metadata"=>()->write(meta, "kaputt = [[["),
                ]
                for (name, mutate) in states
                    mutate()
                    launcher_verdict, package_verdict = both(img, proj)
                    # the tuple form puts the disagreeing verdicts into the failure output
                    @test (name, launcher_verdict) == (name, package_verdict)
                end
                # and the verdicts are the expected ones, not merely equal
                @test SparlectraApp.webui_sysimage_problem(image_path=img, project_dir=proj) == "the sysimage metadata is unreadable"
            end
        end)() end

        @testset "sysimage page tells a native session how to use the image" begin (function ()
            # A Web UI started from the REPL runs without the image and cannot
            # switch to it, not even after a build from its own Sysimage page;
            # the page has to say how the image on disk is used (2026-09-24).
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
                # the first time (2026-09-07): the page only polls while a build is
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
                _progress_file("running", round(time(); digits=1))
                rebuild = SparlectraApp.route_sparlectra_webui("POST", "/webui/sysimage/rebuild"; output_root=out)
                @test rebuild.status == 303
                @test any(h -> first(h) == "Location" && startswith(last(h), "/webui/sysimage"), rebuild.headers)
                @test any(h -> first(h) == "Location" && occursin("already", last(h)), rebuild.headers)
            end
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

        @testset "warmup prefix hides only .jl workloads" begin (function ()
            # the reserved warmup workloads are .jl; a MATPOWER case carrying the
            # prefix (warmup_casePST.m) must stay selectable (regression: the demo
            # case vanished from both selectors after its rename)
            @test !SparlectraApp._webui_is_user_selectable_case("warmup_case3.jl")
            @test !SparlectraApp._webui_is_user_selectable_case("warmup_case118.jl")
            @test SparlectraApp._webui_is_user_selectable_case("warmup_casePST.m")
            @test !SparlectraApp._webui_is_user_selectable_case("warmup_casePST.sparlectra-webui.yaml")
        end)() end

        @testset "buildSysimage one-call API (dry run)" begin (function ()
            # the exported one-liner plans against the package project in a child
            # process; dry run must resolve the target paths without building
            r = SparlectraApp.buildSysimage(dry_run=true, quiet=true)
            @test r.built == false
            @test endswith(r.sysimage_path, SparlectraApp.webui_sysimage_ext()) || occursin("sparlectra.", basename(r.sysimage_path))
            @test occursin("sysimage_meta", basename(r.meta_path))
            # the relocatable executable is a CHECKOUT TOOL since 2026-09-04, not
            # a package function: Sparlectra runs on an installed Julia, so
            # buildApp left the module and the Web UI offers the sysimage only
            @test !isdefined(Sparlectra, :buildApp)
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
                for needle in ("run <casefile>", "se <casefile>", "n1 <casefile>", "--config=", "--set key=value", "working directory no longer exists", "APP_BUILT_AT")
                  occursin(needle, src) || error("missing CLI needle " * needle * " in flavor " * string(fl))
                end
              end
              println("GEN_OK")
              """
            gen_out = read(`$(Base.julia_cmd()) --startup-file=no -e $(gen_code)`, String)
            @test occursin("GEN_OK", gen_out)
            # local documentation viewer: LaTeX must not
            # reach the browser as raw markup, and a thousand-line reference page
            # needs a section index instead of a scrollbar
            math_html = SparlectraApp.render_webui_markdown(raw"Text with $\Omega_{ii} = w_i \cdot \sigma^2$ inline.")
            @test !occursin("\\Omega", math_html)
            @test !occursin("&#36;", math_html)
            @test occursin("Ω", math_html)
            @test occursin("σ", math_html)
            block_html = SparlectraApp.render_webui_markdown("a\n\n```math\n\\frac{P_{se}}{V}\n```\n")
            @test !occursin("frac", block_html)
            @test occursin("/", block_html)
            long_md = read(joinpath(pkgdir(Sparlectra), "docs", "src", "state_estimation.md"), String)
            toc = SparlectraApp._webui_doc_page_toc(long_md)
            @test occursin("On this page", toc)
            @test length(collect(eachmatch(r"<li>", toc))) > 8
            # a short page gets no index (it would be noise)
            @test SparlectraApp._webui_doc_page_toc("# T\n\n## One\n\ntext\n") == ""
            # Documenter cross references (the CGMES page showed its labelled
            # heading as a dead link): the `(@id ...)` label
            # leaves the heading text, `(@ref ...)` becomes the page-local anchor,
            # a label on another served page its /docs route, an unknown target
            # stays disabled; the section index uses the rendered heading ids
            cgmes_md = read(joinpath(pkgdir(Sparlectra), "docs", "src", "cgmes_import.md"), String)
            cgmes_html = SparlectraApp.render_webui_markdown(cgmes_md; current_page="cgmes_import")
            @test !occursin("@id", cgmes_html)
            @test !occursin("@ref topology_processor", cgmes_html)
            @test occursin("id=\"node-breaker-deliveries-without-a-tp-profile\"", cgmes_html)
            @test count("href=\"#node-breaker-deliveries-without-a-tp-profile\"", cgmes_html) == 2
            cgmes_ids = Set(m.captures[1] for m in eachmatch(r"<h[1-6] id=\"([^\"]*)\"", cgmes_html))
            cgmes_toc = SparlectraApp._webui_doc_page_toc(cgmes_md)
            @test occursin(">Node-breaker deliveries without a TP profile<", cgmes_toc)
            @test all(m.captures[1] in cgmes_ids for m in eachmatch(r"href=\"#([^\"]*)\"", cgmes_toc))
            ref_html = SparlectraApp.render_webui_markdown("[taps](@ref transformer-support) [unknown](@ref no_such_label) [fn](@ref)"; current_page="webui")
            @test occursin("href=\"/docs/feature_matrix#transformer-support\"", ref_html)
            @test count("aria-disabled=\"true\"", ref_html) == 2
            # a Documenter-style capitalized anchor reaches the lowercase heading id
            @test occursin("href=\"/docs/scf#configuration-precedence\"", SparlectraApp.render_webui_markdown("[p](scf.md#Configuration-precedence)"; current_page="configuration"))
            # the save-target explanations are small print, not label text
            stg = SparlectraApp.render_settings_page(output_root=mktempdir(), selected_casefile="sp_case14.scf.json")
            @test occursin("settings-target-hint", stg)
            @test !occursin("value=\"this_case\" checked>this case (<code>sp_case14.scf.json</code>): case-scope", stg)

            # header start-flavor: native in the test session, :app when the
            # standalone executable stamped its build time into the environment
            @test SparlectraApp.webui_runtime_flavor().kind in (:native, :sysimage)
            withenv("SPARLECTRA_APP_BUILT" => "2026-08-29T11:11:11") do
                fa = SparlectraApp.webui_runtime_flavor()
                @test fa.kind === :app
                @test fa.built == "2026-08-29T11:11:11"
            end
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
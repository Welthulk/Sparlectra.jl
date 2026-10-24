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

# file: test/test_legacy_webui.jl
# purpose: legacy profile: the testsets of test/test_webui.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here from test_webui.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_webui_tests()
  # cut from testset: Web UI fast parsing, selectors, and rendering
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

        @testset "YAML writers and reader give back every string (0.40.0)" begin (function ()
            # Both writers (the app's, used by the Settings save, and the
            # library's, used by case configuration files) doubled every
            # backslash and the reader kept them: a Windows path, a regular
            # expression and a string with quotes came back changed after one
            # save. One quoting rule now, and the reader resolves exactly its
            # two escapes.
            strings = ["C:\\data\\case.m", raw"[ _#(-]?\d{1,2}\)?$", "say \"hi\" # not a comment"]
            mktempdir() do dir
                for (label, writer) in (("app", (p, d) -> SparlectraApp._write_yaml_file(p, d)), ("library", (p, d) -> write(p, Sparlectra._yaml_dict_text(d))))
                    for v in strings
                        p = joinpath(dir, "rt.yaml")
                        writer(p, Dict{String,Any}("k" => v, "list" => Any[v, "plain"]))
                        back = Sparlectra.load_yaml_dict(p)
                        @test (label, back["k"]) == (label, v)
                        @test (label, back["list"]) == (label, Any[v, "plain"])
                    end
                    # a key that is no bare key is written in quotes, so a
                    # weight of "Bus 1", "a:b", Bus "N" or C:\ reads back
                    # after a save, also next to a value with a "#"; below
                    # the top level a key may hold a dot; a key the reader
                    # cannot give back (both quote kinds) is refused
                    keyed = Dict{String,Any}("Bus 1" => 2.0, "a:b" => 1, "Bus \"N\"" => "x # y", "C:\\" => "x # y", "N.1" => 1.5, "10.5" => Dict{String,Any}("p_mw" => 0.5), "max_iter" => 3)
                    p = joinpath(dir, "keys.yaml")
                    writer(p, Dict{String,Any}("weights" => keyed))
                    @test (label, Sparlectra.load_yaml_dict(p)["weights"]) == (label, keyed)
                    @test_throws "YAML key \"say \\\"it's\\\"\" cannot be written" writer(p, Dict{String,Any}("say \"it's\"" => 1))
                end
                # a run writes its effective configuration with the app
                # writer: a weight key in single quotes (Bus "N", C:\)
                # aborted the run there
                case5 = joinpath(dir, "sp_case5.scf.json")
                cp(joinpath(pkgdir(Sparlectra), "data", "scf", "sp_case5.scf.json"), case5)
                wcfg = joinpath(dir, "weights.yaml")
                write(wcfg, "config_version: 1\npower_flow:\n  distributed_slack:\n    weights:\n      'Bus \"N\"': 2.0\n      'C:\\': 1.0\n")
                out = mkpath(joinpath(dir, "run"))
                res = SparlectraApp.run_sparlectra_api(; casefile = case5, config_file = wcfg, output_dir = out)
                @test res.success
                @test Sparlectra.load_yaml_dict(joinpath(out, "effective_config.yaml"))["power_flow"]["distributed_slack"]["weights"] == Dict{String,Any}("Bus \"N\"" => 2.0, "C:\\" => 1.0)
                # the Settings save rewrites the whole file: a stem pattern
                # with backslashes (machine scope, not on the form) survives a
                # save of another key, and hand-written Windows paths and
                # regular expressions read as written
                cfg = joinpath(dir, "configuration.yaml")
                pattern = raw"[ _#(-]?\d{1,2}\)?$"
                write(cfg, replace(read(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, String), "name_stem_pattern: \"[ _#(-]?[0-9]{1,2}[)]?\$\"" => "name_stem_pattern: \"$(pattern)\""))
                @test Sparlectra.load_sparlectra_config(cfg; reload = true).contingency.events.name_stem_pattern == pattern
                saved = SparlectraApp._webui_write_general_settings!(cfg, Dict{String,Any}("contingency.events.enabled" => true))
                @test saved.ok
                again = Sparlectra.load_sparlectra_config(cfg; reload = true).contingency.events
                @test again.enabled && again.name_stem_pattern == pattern
                # both pattern forms strip the same suffixes
                stems(form) = [replace(n, Regex(form) => "") for n in ("Ostheim_Neuwiese_1", "Neubach_Waldsee", "Feld (2)", "Leitung_220")]
                @test stems(pattern) == stems(Sparlectra.ContingencyEventsConfig().name_stem_pattern)
                @test stems(pattern)[1:3] == ["Ostheim_Neuwiese", "Neubach_Waldsee", "Feld "]
            end
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

            # the editor page seeds the case's real element names
            page = redirect_stdout(devnull) do
                SparlectraApp.handle_contingency_weights_page(Dict{String,Any}("case" => "warmup_casePST.m"); output_root=root, case_directory=cases, operation_log=root)
            end
            body = String(page.body)
            @test page.status == 200
            @test occursin(bnames[1], body)

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
            # shipped sp_case* demos carry active controllers, which an apslf
            # run leaves inactive with one warning that names apslf (0.31.0);
            # that warning in run.log shows the sidecar's solver choice reached
            # the run.
            write(joinpath(casedir, "sp_case14.config.yaml"), "config_version: 1\nscope: case\ncase: sp_case14.scf.json\npower_flow:\n  solver: apslf\n")
            run_solver = SparlectraApp.start_powerflow_run(Dict("casefile" => "sp_case14.scf.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root); case_directory=casedir)
            @test get(run_solver, "success", false) === true
            @test occursin("power_flow.solver=apslf runs without controllers", read(joinpath(run_solver["output_dir"], "run.log"), String))

            # "estimate taps" is the key state_estimation.tap_estimation.enabled,
            # checked by its default true on a case with tap changers
            # (sp_case14) and on one without (feeder3_hardpv_A, lines
            # only, where it has nothing to release); a value saved for the
            # case wins over that default
            f3 = "feeder3_hardpv_A.scf.json"
            cp(abspath(joinpath(dirname(@__DIR__), "data", "scf", f3)), joinpath(casedir, f3))
            SparlectraApp._se_generate_measurement_set(joinpath(casedir, f3), joinpath(casedir, "feeder3_hardpv_A.measurements.csv"), SparlectraApp.MeasurementGeneratorOptions(; noise=false, gross_k=0.0, tap_steps=0.0, include_i=false, sigma_u_pct=0.5, sigma_i_pct=1.0, sigma_p_pct=1.0, sigma_q_pct=1.0, sigma_ia_deg=0.0))
            tapbox(case) = match(r"name=\"se_tap_estimation\" value=\"true\"( checked)?>", String(copy(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=$(case)"; output_root=root, runtime=rt).body)))
            @test tapbox("sp_case14.scf.json") !== nothing && tapbox("sp_case14.scf.json")[1] !== nothing
            @test tapbox(f3) !== nothing && tapbox(f3)[1] !== nothing
            SparlectraApp._webui_merge_case_settings!(root, joinpath(casedir, "sp_case14.scf.json"), Dict{String,Any}("se_tap_estimation" => false); config_file=Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, case_directory=casedir)
            @test tapbox("sp_case14.scf.json")[1] === nothing
        end)() end

        @testset "the Info menu stays in the header on every page" begin (function ()
            # Regression: the Info control was passed to the layout by
            # the Case, Settings and Runs renderers only, so it disappeared as soon
            # as the user clicked Operation Log, Run history, Last errors, Docs or a
            # result page. Everything it shows describes the running server, not the
            # page, so the header builds it itself now.
            root = mktempdir()   # never the user's real state directory: the pages log their visits there
            for path in ("/powerflow", "/powerflow/case", "/powerflow/settings", "/powerflow/history",
                "/webui/operation-log", "/webui/last-errors", "/webui/sysimage")
                response = SparlectraApp.route_sparlectra_webui("GET", path; output_root=root)
                body = String(response.body)
                # the tuple form names the offending page in the failure output
                @test (path, occursin("topbar-info-menu", body)) == (path, true)
                # the box names the Julia the server runs on: with
                # 1.12 and 1.13 both in use, a screenshot has to say which
                @test occursin("<dt>Julia</dt><dd><code>$(VERSION)</code></dd>", body)
                # and the AnalyticLoadFlow it solves with, so an outdated
                # environment is visible on screen
                @test occursin("<dt>AnalyticLoadFlow</dt><dd><code>v$(pkgversion(Sparlectra.AnalyticLoadFlow))</code> <a href=\"https://github.com/SOPTIM/AnalyticLoadFlow.jl\"", body)
                # and the application package, which has its own version;
                # the header bar keeps naming the library
                @test occursin("<dt>Application</dt><dd><code>SparlectraApp v$(pkgversion(SparlectraApp))</code></dd>", body)
                @test occursin("Sparlectra.jl v$(Sparlectra.version())", body)
            end
            # error pages carry it as well: they are where a user looks for the
            # output root and the operation log in the first place
            @test occursin("topbar-info-menu", SparlectraApp.render_webui_error(404, "not found"))
        end)() end

        @testset "the environment is checked before the sysimage question" begin (function ()
            # Regression (Windows 11). A checkout carried a Manifest.toml
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
            # Regression (Windows 11). A checkout carried a Manifest.toml
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

  return true
end

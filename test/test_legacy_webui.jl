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
#          They moved here unchanged from test_webui.jl; the helpers they use
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
            @test list_for("sp_case14.scf.json") == [(name = n, case = c) for (n, c) in (("blank_set.csv", "sp_case14.scf.json"), ("bom_set.csv", "sp_case14.scf.json"), ("crlf_set.csv", "sp_case14.scf.json"), ("demo_clean.csv", "sp_case14.scf.json"), ("sp_case14.measurements.csv", "sp_case14.scf.json"), ("sp_case14.old.csv", ""))]
            @test list_for("case14.scf.json") == [(name = "case14.scf.measurements.csv", case = "case14.scf.json"), (name = "sp_case14.foreign.csv", case = "case14.scf.json")]
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
            # a CSV v1 (reader removed in 0.50.0) carries the removal notice
            # naming the converter at the file, a table not
            v1_state = SparlectraApp._webui_se_form_state(Dict{String,Any}("case" => "warmup_casePST.m", "measurement_file" => "meas_v1.csv"); output_root=root, case_directory=cases)
            @test v1_state.set_v1 && !tbl_state.set_v1
            v1_page = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m&measurement_file=meas_v1.csv", Dict{String,String}(); output_root=root, runtime=rt).body)
            @test occursin("measurement-v1-notice", v1_page) && occursin("convertMeasurementsV1", v1_page)
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
            @test occursin("Run state estimation", page2)
            @test !occursin("onsubmit", page2)

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
            @test occursin("Run state estimation", flow_page)
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

            # every SE input parameter carries a hover hint, and the technical
            # ones a documentation link
            @test count("help-link", page2) >= 12
            for topic in ("state_estimation.tol", "state_estimation.flatstart", "state_estimation.max_iter", "state_estimation.robust", "state_estimation.update_shunts", "state_estimation.report_residual_correlation", "webui.se_max_eliminations", "webui.se_measurement_file", "webui.se_generator_noise", "webui.se_generator_gross_error", "webui.se_generator_tap_error", "webui.se_generator_sigmas")
                se_topic = SparlectraApp.resolve_webui_help_topic(topic)
                @test se_topic !== nothing && 0 < length(se_topic.hint) <= 160
            end
            @test occursin("state_estimation_measurements/#se-generator-tap-error", SparlectraApp.webui_help_doc_url("webui.se_generator_tap_error"))

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
            # regenerate the clean set: another new file with the same content
            gen3 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root=root, runtime=rt)
            gen3_file = joinpath(cases, written_set(gen3))
            @test read(gen3_file, String) == plain   # seeded generator is reproducible
            rm(gen3_file)

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
            @test occursin("U ×", pageInfo)
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
            @test occursin("Edit measurement values (table)", pageTab)
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
            @test occursin("Machine and unobservable taps are skipped", pageInfo)
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
            # both sets are bound to the case; the NEWEST is armed (the copy, by
            # its modification time), the other one is offered plain
            @test occursin("value=\"case9.measurements.csv\" selected", page3)
            @test occursin("value=\"$(mname)\">", page3)
            # block 4 dropped the SE-own case selector (shared selection), and its
            # "*" bound-set marker went with it; the binding stays visible through
            # the Case binding line and the per-set labels asserted here
            @test occursin("Case binding:", page3)
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
            @test occursin("Edit measurement file (inline)", pageEd)
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
            # a run directory from before 0.50.0 keeps its artifact in the
            # measurement CSV v1 layout: the page still renders, and the test
            # answers with a note instead of an error
            art1 = joinpath(root, id1, "measurements.csv")
            art1_table = read(art1, String)
            write(art1, "# sparlectra-measurements v1\ntype,bus,from_bus,to_bus,branch_nr,link_nr,direction,value,sigma,active,id,bad_data,status\nVmMeas,1,,,,,,1.0,0.01,true,v1,,\n")
            rold = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/topology-hypotheses", Dict{String,Any}("run_id" => id1); output_root=root, runtime=rt)
            @test rold.status == 303 && occursin("before%200.50.0", Dict(rold.headers)["Location"])
            @test SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(id1)", Dict{String,String}(); output_root=root, runtime=rt).status == 200
            write(art1, art1_table)

            # tap-estimation SE run: a set generated with a tap deviation
            # disagrees with the model around one transformer; releasing the taps
            # absorbs the discrepancy and the fixation lands on a mechanical step
            gentap2 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "tap_error_steps" => "3"); output_root=root, runtime=rt)
            rm(mfile)
            mname = written_set(gentap2)
            mfile = joinpath(cases, mname)
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
            # the injected deviation is a CHANGE against the model position:
            # the tap CSV says so, the result table marks the row as corrected
            change_col = findfirst(==("change"), header_cols)
            @test change_col !== nothing
            @test any(l -> (f=split(l, ","); f[change_col] != "" && f[change_col] != "0"), taplines[2:end])
            tapPage = String(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/result/$(idtap)", Dict{String,String}(); output_root=root, runtime=rt).body)
            @test occursin("Transformer tap estimates", tapPage)
            @test occursin("fixed (corrected)", tapPage)
            @test occursin("<th>Model step</th>", tapPage)
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

            # bad data lands findable: gross error in the set -> the verdict
            # file flags the corrupted row of measurements.csv as eliminated,
            # joined by its line number; the row is the one whose value is
            # 10 sigma off its truth value (truth_value comment of the set)
            mfile, mname = regenerate(Dict{String,Any}("casefile" => "warmup_casePST.m", "gross_error_k" => "10"))
            rbd = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
            @test rbd["status"] == "succeeded"
            vrows = [split(l, ",") for l in readlines(joinpath(root, rbd["run_id"], "measurements_verdict.csv"))[2:end]]
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
            # a case on constraint (equality-constrained zero injection): the
            # box shows checked, a generation with the box checked keeps
            # constraint, unchecked writes off; checked again it is on
            @test [SparlectraApp._webui_normalize_case_profile_form_value("gen_passive_as_zi", v) for v in ("constraint", "off", "on")] == [true, false, true]
            SparlectraApp._webui_merge_case_settings!(root, case_path, Dict{String,Any}("gen_passive_as_zi" => "constraint"); case_directory=cases)
            @test Sparlectra._as_symbol_cfg(Sparlectra.load_case_config(case_path)["state_estimation.zero_injection"]) === :constraint
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
            sestate = readlines(joinpath(root, id1, "se_state.csv"))
            strow = only([l for l in sestate if startswith(l, "2,")])
            st_vm = split(strow, ",")[2]
            # the table carries U in kV; read back it is the pu value of the
            # state file bit for bit
            frnet = Sparlectra._import_sparlectra_net(case_path, nothing, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload=true))
            readMeasurementTable!(frnet; file=mfile, zero_injection=false)
            @test only(m for m in frnet.measurements if m.typ == Sparlectra.VmMeas && m.busIdx == 2).value == parse(Float64, st_vm)
            # rejections: tap deviation locked, unknown run id named
            rrej = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_truth_source" => "from_run", "gen_truth_run_id" => id1, "tap_error_steps" => "2"); output_root=root, runtime=rt)
            @test occursin("requires%20truth%20state", Dict(rrej.headers)["Location"])
            rrej2 = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_truth_source" => "from_run", "gen_truth_run_id" => "nope"); output_root=root, runtime=rt)
            @test occursin("not%20found%20in%20the%20run%20history", Dict(rrej2.headers)["Location"])

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
            # the note column is quoted and carries commas ("islanded: bus 8
            # cut off, 0.0 MW load ..."); a plain split at the delimiter
            # shifted the columns after it and made the number of checks
            # depend on the note text (webui 3015 -> 3014 when an N-1 note
            # changed). The quote-aware reader of the Web UI keeps one field
            # per column.
            ncols = length(SparlectraApp._webui_compare_split_csv(first(lines3), csv_delim))
            for (l3, l4) in zip(lines3, lines4)
                f3 = SparlectraApp._webui_compare_split_csv(l3, csv_delim)
                f4 = SparlectraApp._webui_compare_split_csv(l4, csv_delim)
                @test length(f3) == length(f4) == ncols
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
                # the citation notes sit side by side, MATPOWER's and Sparlectra's
                @test occursin("matpower-citation-note", body)
                @test occursin("sparlectra-citation-note", body)
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

        @testset "warmup prefix hides only .jl workloads" begin (function ()
            # the reserved warmup workloads are .jl; a MATPOWER case carrying the
            # prefix (warmup_casePST.m) must stay selectable (regression: the demo
            # case vanished from both selectors after its rename)
            @test !SparlectraApp._webui_is_user_selectable_case("warmup_case3.jl")
            @test !SparlectraApp._webui_is_user_selectable_case("warmup_case118.jl")
            @test SparlectraApp._webui_is_user_selectable_case("warmup_casePST.m")
            @test !SparlectraApp._webui_is_user_selectable_case("warmup_casePST.sparlectra-webui.yaml")
        end)() end

  return true
end

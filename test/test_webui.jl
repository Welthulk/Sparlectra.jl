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
  @testset "Web UI fast parsing, selectors, and rendering" begin
    @test Sparlectra._webui_parse_bool("on") === true
    @test Sparlectra._webui_parse_bool("false") === false
    @test Sparlectra.route_sparlectra_webui("GET", "/powerflow/artifact/..%2Fbad/result.json").status in (400, 404)

    @testset "git commit sha resolution (issue #290)" begin
      sha_a = "a" ^ 40
      sha_b = "b" ^ 40

      mktempdir() do root
        @test Sparlectra._git_head_commit_sha(root) === nothing

        # Plain checkout: .git directory with a loose ref.
        gitdir = joinpath(root, ".git")
        mkpath(joinpath(gitdir, "refs", "heads"))
        write(joinpath(gitdir, "HEAD"), "ref: refs/heads/main\n")
        write(joinpath(gitdir, "refs", "heads", "main"), sha_a * "\n")
        @test Sparlectra._git_head_commit_sha(root) == sha_a

        # Detached HEAD: plain SHA in HEAD.
        write(joinpath(gitdir, "HEAD"), sha_b * "\n")
        @test Sparlectra._git_head_commit_sha(root) == sha_b

        # Packed ref: no loose ref file, SHA only in packed-refs.
        write(joinpath(gitdir, "HEAD"), "ref: refs/heads/packed\n")
        write(joinpath(gitdir, "packed-refs"), "# pack-refs with: peeled fully-peeled sorted\n$(sha_b) refs/heads/packed\n^$(sha_a)\n")
        @test Sparlectra._git_head_commit_sha(root) == sha_b

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
        @test Sparlectra._git_head_commit_sha(worktree_root) == sha_a

        write(joinpath(worktree_root, ".git"), "not a git pointer\n")
        @test Sparlectra._git_head_commit_sha(worktree_root) === nothing
      end

      # In a git checkout of Sparlectra itself (plain or worktree) the Web UI
      # banner SHA must resolve; a registry install (no .git) yields nothing.
      package_root = normpath(joinpath(dirname(pathof(Sparlectra)), ".."))
      if ispath(joinpath(package_root, ".git"))
        sha = Sparlectra._sparlectra_git_commit_sha()
        @test sha isa String
        @test occursin(r"^[0-9a-f]{40}$", sha)
      end
    end

    mktempdir() do root
      case_directory = joinpath(root, "cases")
      mkpath(case_directory)
      write(joinpath(case_directory, "case14.m"), "function mpc = case14\nend\n")
      write(joinpath(case_directory, "FOR001.DAT"), _dtf_network_fixture())
      write(joinpath(case_directory, "FOR001_OUTAGES.DAT"), _dtf_network_with_outage_fixture())
      write(joinpath(case_directory, "OUTAGE.DAT"), _dtf_outage_fixture())
      write(joinpath(case_directory, "FOR002.DAT"), _dtf_reference_fixture())
      write(joinpath(case_directory, "UNKNOWN.DAT"), "plain unsupported data\n")
      @test Sparlectra._webui_classify_dat_content(joinpath(case_directory, "FOR001.DAT")) === :dtf_network_case
      @test Sparlectra._webui_classify_dat_content(joinpath(case_directory, "FOR001_OUTAGES.DAT")) === :dtf_network_case_with_outages
      @test Sparlectra._webui_classify_dat_content(joinpath(case_directory, "OUTAGE.DAT")) === :dtf_outage_file
      @test Sparlectra._webui_classify_dat_content(joinpath(case_directory, "FOR002.DAT")) === :dtf_outage_or_reference
      @test Sparlectra._webui_classify_dat_content(joinpath(case_directory, "UNKNOWN.DAT")) === :unknown_dat
      primary = Sparlectra._webui_casefile_options_in_directory(case_directory)
      reference = Sparlectra._webui_for002_reference_options_in_directory(case_directory)
      @test "FOR001.DAT" in primary
      @test "FOR001_OUTAGES.DAT" in primary
      @test !("OUTAGE.DAT" in primary)
      @test !("UNKNOWN.DAT" in primary)
      @test reference == ["FOR002.DAT"]

      active_html = Sparlectra.render_powerflow_result(Dict("run_id" => "active", "status" => "running", "elapsed_seconds" => 1.25))
      @test occursin("Run status", active_html)
      @test occursin("Elapsed time", active_html)
      @test !occursin("Solver time", active_html)
      @test !occursin("Wall time", active_html)

      terminal = Dict("run_id" => "ok", "status" => "completed", "solver_elapsed_s" => 0.125, "service_phase_timings" => [Dict("phase" => "total_service", "elapsed_seconds" => 1.5)])
      terminal_html = Sparlectra.render_powerflow_result(terminal)
      @test occursin("Solver time", terminal_html)
      @test occursin("Total time", terminal_html)
      @test !occursin("Elapsed time", terminal_html)
      @test !occursin("Wall time", terminal_html)
      @test !occursin("Solver time: n/a", terminal_html)

      failed_solver_html = Sparlectra.render_powerflow_result(Dict("run_id" => "failed", "status" => "failed", "metadata" => Dict("solver_elapsed_s" => 0.25), "service_phase_timings" => [Dict("phase" => "total_service", "elapsed_seconds" => 2.0)]))
      @test occursin("Solver time", failed_solver_html)
      presolver_html = Sparlectra.render_powerflow_result(Dict("run_id" => "bad", "status" => "failed", "service_phase_timings" => [Dict("phase" => "total_service", "elapsed_seconds" => 0.1)]))
      @test !occursin("Solver time", presolver_html)
      @test occursin("Total time", presolver_html)
      legacy_html = Sparlectra.render_powerflow_result(Dict("run_id" => "old", "status" => "completed", "elapsed_seconds" => 3.0))
      @test occursin("Total time", legacy_html)
    end

    @testset "CGMES export checkbox and result row" begin
      # stage 4A block 3: the export-CGMES checkbox renders on Settings
      form_html = Sparlectra.render_settings_page(output_root = mktempdir())
      @test occursin("name=\"export_cgmes\"", form_html)
      @test Sparlectra.resolve_webui_help_topic("webui.export_cgmes") !== nothing
      # excerpt loading needs the cgmes_export page in WEBUI_DOC_PAGES and the
      # "Export from the Web UI" heading in docs/src/cgmes_export.md
      @test Sparlectra.load_webui_help_excerpt("webui.export_cgmes") !== nothing

      req = Sparlectra.powerflow_webui_request(Dict("casefile" => "case.m", "export_cgmes" => "true"))
      @test req["export_cgmes"] === true
      # unchecked checkbox: browsers submit only the hidden false field
      req = Sparlectra.powerflow_webui_request(Dict("casefile" => "case.m", "export_cgmes" => "false"))
      @test req["export_cgmes"] === false
      # absent field falls back to the spec default (false)
      req = Sparlectra.powerflow_webui_request(Dict("casefile" => "case.m"))
      @test req["export_cgmes"] === false

      completed_html = Sparlectra.render_powerflow_result(
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

      failed_html = Sparlectra.render_powerflow_result(Dict("run_id" => "y", "status" => "completed", "metadata" => Dict("cgmes_export_status" => "failed", "cgmes_export_error" => "boom")))
      @test occursin("failed — boom", failed_html)

      plain_html = Sparlectra.render_powerflow_result(Dict("run_id" => "z", "status" => "completed"))
      @test !occursin("CGMES export", plain_html)
    end

    @testset "request builder falls back to the case form block (stage 4A)" begin
      # fields the run page no longer renders must reach the run from the
      # case configuration file's form block; a POSTed field always wins
      dir = mktempdir()
      case_path = joinpath(dir, "case_formblock.m")
      write(case_path, "function mpc = case_formblock\nend\n")
      Sparlectra.write_case_config(case_path, Dict{String,Any}(); form = Dict{String,Any}("export_cgmes" => true, "performance_timing" => "full", "case_format" => "matpower", "gen_seed" => 7))
      req = Sparlectra.powerflow_webui_request(Dict("casefile" => case_path))
      @test req["export_cgmes"] === true
      @test req["performance_timing"] == "full"
      @test req["case_format"] == "matpower"
      req_posted = Sparlectra.powerflow_webui_request(Dict("casefile" => case_path, "export_cgmes" => "false", "case_format" => "auto"))
      @test req_posted["export_cgmes"] === false
      @test req_posted["case_format"] == "auto"
      # relative case name resolves through case_directory
      req_rel = Sparlectra.powerflow_webui_request(Dict("casefile" => "case_formblock.m"); case_directory = dir)
      @test req_rel["export_cgmes"] === true
      # without a stored format the hint chain answers: a free-typed .DAT
      # yields the explicit DTF format, everything else stays auto
      req_dat = Sparlectra.powerflow_webui_request(Dict("casefile" => "no_such_FOR001.DAT"); case_directory = dir)
      @test req_dat["case_format"] == "dtf_for001"
      @test Sparlectra.powerflow_webui_request(Dict("casefile" => "case.m"))["case_format"] == "auto"
      # the result-page save must not wipe the stored non-spec field: the
      # preservation loop carries case_format like the spec fields
      stored = Sparlectra._webui_case_form_defaults(case_path, nothing)
      @test stored["case_format"] == "matpower"
      @test stored["gen_seed"] == 7
    end

    @testset "shared selected-case memory (stage 4A harmonization)" begin
      # choosing on the Case page must reach the run page and the SE page
      # through plain nav links (no query): the maintainer's live run hit
      # "no case selected" exactly this way
      root = mktempdir()
      cases = joinpath(root, "cases")
      mkpath(cases)
      # suite migration (no-download rule): this set only exercises the
      # shared-memory routing by NAME, so a fixture file suffices
      write(joinpath(cases, "case14.m"), "% case fixture\n")
      rt = (; case_directory = cases, config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log = Sparlectra.webui_operation_log_path(root), startup_config_error = nothing, runner = Sparlectra.start_powerflow_run)
      # nothing remembered yet: the run page renders without a case
      fresh = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow", Dict{String,String}(); output_root = root, runtime = rt).body)
      @test occursin("no case selected", fresh)
      # opening the Case page with an explicit case remembers it (GET
      # queries travel IN the target, the third argument is a POST body)
      Sparlectra.route_sparlectra_webui("GET", "/powerflow/case?casefile=case14.m"; output_root = root, runtime = rt)
      @test Sparlectra._webui_recall_selected_case(root) == "case14.m"
      # the plain-nav run page now carries the remembered case as its hidden field
      run_page = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow", Dict{String,String}(); output_root = root, runtime = rt).body)
      @test occursin("type=\"hidden\" name=\"casefile\" value=\"case14.m\"", run_page)
      # stage 4A block 4: /stateestimation is a real redirect onto the Runs
      # page, whose SE section consumes the same memory
      se_redirect = Sparlectra.route_sparlectra_webui("GET", "/stateestimation", Dict{String,String}(); output_root = root, runtime = rt)
      @test se_redirect.status == 303
      @test Dict(se_redirect.headers)["Location"] == "/powerflow#state-estimation"
      @test occursin("id=\"state-estimation\"", run_page)
      # an empty run POST now points at the Case page, format-neutral
      err = try
        Sparlectra.powerflow_webui_request(Dict{String,Any}())
        nothing
      catch e
        e
      end
      @test err isa ArgumentError
      @test occursin("Case page", sprint(showerror, err))
      @test !occursin("MATPOWER", sprint(showerror, err))
    end

    @testset "settings save reaches the run without POST fields (stage 4A block 3)" begin
      # THE stage assumption (colleague: riskiest single assumption of the
      # stage): a value saved on the Settings page must reach a run whose
      # POST no longer carries the field, via resolve_config, and the
      # effective-config artifact must name case_sidecar as its source.
      root = mktempdir()
      cache = joinpath(root, "cases")
      app_root = normpath(joinpath(dirname(@__DIR__)))
      Sparlectra._webui_stage_bundled_case!(app_root, cache, "sp_case14.scf.json")
      rt = (; case_directory = cache, config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log = Sparlectra.webui_operation_log_path(root), startup_config_error = nothing, runner = Sparlectra.start_powerflow_run)
      resp = Sparlectra.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "power_flow_max_iter" => "44"); output_root = root, runtime = rt)
      @test resp.status in (302, 303)
      @test Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json"))["power_flow.max_iter"] == 44
      # the block-3 run form posts NO override fields; the bare request
      # must carry none either
      req = Sparlectra.powerflow_webui_request(Dict("casefile" => "sp_case14.scf.json"); case_directory = cache)
      @test isempty(req["config_overrides"])
      run = Sparlectra.start_powerflow_run(merge(req, Dict("output_root" => joinpath(root, "runs"), "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH)); case_directory = cache)
      @test run["status"] == "succeeded"
      eff = read(joinpath(String(run["output_dir"]), "effective_config.yaml"), String)
      seg = eff[first(findfirst("  power_flow:", eff)):end]
      mi = seg[first(findfirst("    max_iter:", seg)):first(findfirst("    max_iter:", seg)) + 220]
      @test occursin("value: 44", mi)
      @test occursin("source: case_sidecar", mi)
      # machine-scope keys are named and kept out of the case file
      resp2 = Sparlectra.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("casefile" => "sp_case14.scf.json", "settings_target" => "this_case", "benchmark_samples" => "5"); output_root = root, runtime = rt)
      @test occursin("benchmark.samples", Sparlectra._webui_urldecode(Dict(resp2.headers)["Location"]))
      @test !haskey(Sparlectra.load_case_config(joinpath(cache, "sp_case14.scf.json")), "benchmark.samples")
      # the general target merges into the YAML with a backup
      cfg = joinpath(root, "configuration.yaml")
      cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, cfg)
      resp3 = Sparlectra.route_sparlectra_webui("POST", "/powerflow/settings/save", Dict{String,Any}("settings_target" => "general", "config_file" => cfg, "power_flow_autodamp_min" => "0.09"); output_root = root)
      @test resp3.status in (302, 303)
      @test occursin("autodamp_min: 0.09", read(cfg, String))
      @test isfile(cfg * ".settings-save.bak")
    end

    @testset "no nested forms anywhere on the run page" begin
      # HTML forbids nested <form>. A browser closes the outer form where the
      # inner one starts, so every control AFTER it — including "Start
      # PowerFlow run" — silently falls out of the form and does nothing when
      # clicked. This bit exactly once; the depth check keeps it from
      # returning through any future in-form button.
      dir = mktempdir()
      prof = joinpath(dir, "c.sparlectra-webui.yaml")
      write(prof, "placeholder")
      for html in (
        Sparlectra.render_powerflow_form(output_root = mktempdir()),
        Sparlectra.render_powerflow_form(output_root = mktempdir(), selected_casefile = "c.m", case_profile = Dict{String,Any}("power_flow_solver" => "dc", "_profile_path" => prof)),
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
    end

    @testset "saved case settings can be reset from the form" begin
      # Saved settings outrank the configuration for their keys, so a stale
      # sidecar can pin a case to a setting the user cannot override in the
      # form (measured: a delivery stuck on power_flow_solver: dc). The reset
      # path must be reachable independently of the dismissible notice.
      dir = mktempdir()
      prof = joinpath(dir, "c.sparlectra-webui.yaml")
      write(prof, "placeholder")
      with_sidecar = Sparlectra.render_settings_page(output_root = mktempdir(), selected_casefile = "c.m", case_profile = Dict{String,Any}("power_flow_solver" => "dc", "_profile_path" => prof))
      @test occursin("Reset saved settings for this case", with_sidecar)
      @test occursin("/powerflow/case-settings/reset", with_sidecar)
      without = Sparlectra.render_settings_page(output_root = mktempdir())
      @test !occursin("Reset saved settings for this case", without)
      # switching the case reloads the page server-side; the wait must be
      # visible WHERE the switching happens, which since stage 4A is the
      # Case page's chooser
      @test occursin("case-loading-banner", Sparlectra.render_case_page(output_root = mktempdir()))

      # The handler deletes the sidecar and keeps the case file.
      root = joinpath(dir, "runs")
      cases = joinpath(dir, "cases")
      mkpath(root)
      mkpath(cases)
      write(joinpath(cases, "c.m"), "function mpc = c\nend\n")
      sc = Sparlectra._webui_case_settings_path(root, "c.m"; case_directory = cases)
      mkpath(dirname(sc))
      write(sc, "values:\n  power_flow_solver: dc\n")
      response = Sparlectra.handle_powerflow_case_settings_reset(Dict("casefile" => "c.m"); output_root = root, case_directory = cases, operation_log = root)
      @test response.status == 303
      @test !isfile(sc)
      @test isfile(joinpath(cases, "c.m"))
      # idempotent: a second reset is a no-op, not an error
      @test Sparlectra.handle_powerflow_case_settings_reset(Dict("casefile" => "c.m"); output_root = root, case_directory = cases, operation_log = root).status == 303
      # path traversal is rejected
      @test Sparlectra.handle_powerflow_case_settings_reset(Dict("casefile" => "../evil.m"); output_root = root, case_directory = cases, operation_log = root).status == 303
    end

    @testset "operation log: clear from the page" begin
      root = mktempdir()
      log = Sparlectra.webui_operation_log_path(root)
      rt = (; case_directory = mktempdir(), config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log = log, startup_config_error = nothing, runner = Sparlectra.start_powerflow_run)
      Sparlectra.record_webui_operation!(log, "probe"; route = "/x", method = "GET")
      Sparlectra.record_webui_operation!(log, "probe"; route = "/y", method = "GET")
      page = String(copy(Sparlectra.route_sparlectra_webui("GET", "/webui/operation-log"; output_root = root, runtime = rt).body))
      @test occursin("/webui/operation-log/clear", page)
      @test occursin("entries,", page)                      # the size is stated on the page
      before = length(readlines(log))          ## opening the page logs one entry too
      resp = Sparlectra.route_sparlectra_webui("POST", "/webui/operation-log/clear"; output_root = root, runtime = rt)
      @test resp.status == 303
      # the file is emptied but keeps ONE entry recording the deletion, so the
      # log never becomes silently empty
      remaining = readlines(log)
      @test length(remaining) == 1
      @test occursin("operation_log_cleared", remaining[1])
      @test occursin("\"removed_entries\":$(before)", replace(remaining[1], " " => ""))
    end

    @testset "case download follows a symlinked case directory" begin
      # the Web UI state directory is reachable through a symlink (a Flatpak
      # app data path pointing at ~/.local/state), so the form can carry the
      # linked spelling while the runtime holds the resolved one
      root = mktempdir()
      cases = joinpath(root, "cases")
      mkpath(cases)
      cp(abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m")), joinpath(cases, "warmup_casePST.m"))
      linked = joinpath(root, "linked")
      symlink(cases, linked)
      rt = (; case_directory = cases, config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log = Sparlectra.webui_operation_log_path(root), startup_config_error = nothing, runner = Sparlectra.start_powerflow_run)
      via_link = joinpath(linked, "warmup_casePST.m")
      dl = Sparlectra.route_sparlectra_webui("GET", "/powerflow/case/download?case=$(Sparlectra._webui_urlencode(via_link))"; output_root = root, runtime = rt)
      @test dl.status == 200
      @test String(copy(dl.body)) == read(joinpath(cases, "warmup_casePST.m"), String)
      # resolving symlinks must not open the rest of the file system
      outside = Sparlectra.route_sparlectra_webui("GET", "/powerflow/case/download?case=$(Sparlectra._webui_urlencode("/etc/passwd"))"; output_root = root, runtime = rt)
      @test outside.status == 303
      @test isempty(outside.body)
    end

    @testset "operation-log retention" begin
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
      Sparlectra._prune_webui_operation_log!(log; Sparlectra._webui_operation_log_options(; retention_days = 10)...)
      kept = readlines(log)
      @test length(kept) == 1
      @test occursin("recent", kept[1])
      # a shorter retention drops more, which is the knob for an unwieldy log
      Sparlectra._prune_webui_operation_log!(log; Sparlectra._webui_operation_log_options(; retention_days = 1)...)
      @test isempty(readlines(log))
      # the configuration carries it, with the documented default
      cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
      @test cfg.webui.operation_log_retention_days == 10
      lowered = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true,
        overrides = Dict{String,Any}("webui" => Dict{String,Any}("operation_log_retention_days" => 3)))
      @test lowered.webui.operation_log_retention_days == 3
      @test Sparlectra._webui_operation_log_options(; retention_days = 3).retention_days == 3
      @test_throws ArgumentError Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true,
        overrides = Dict{String,Any}("webui" => Dict{String,Any}("operation_log_retention_days" => -1)))
    end

    @testset "contingency weights editor and storage (#331 Phase 5 follow-up)" begin
      dir = mktempdir()
      cases = joinpath(dir, "cases")
      mkpath(cases)
      root = joinpath(dir, "runs")
      mkpath(root)
      # load_fixture_net: the tracked PST warmup case gives the editor real
      # element names without any download
      cp(joinpath(pkgdir(Sparlectra), "data", "mpower", "warmup_casePST.m"), joinpath(cases, "warmup_casePST.m"))
      wf = Sparlectra._webui_case_weights_path("warmup_casePST.m"; case_directory = cases)

      # the weight file lives next to the case as <stem>.contingency-weights.csv
      @test basename(wf) == "warmup_casePST.contingency-weights.csv"
      @test dirname(wf) == normpath(cases)

      # list exclusion: a weights file must not be offered as a selectable case
      touch(wf)
      opts = Sparlectra._webui_casefile_options_in_directory(cases)
      @test "warmup_casePST.m" in opts
      @test !any(occursin("contingency-weights", o) for o in opts)
      rm(wf)

      # real element names for the fixtures
      net = redirect_stdout(devnull) do
        Sparlectra._import_sparlectra_net(joinpath(cases, "warmup_casePST.m"), nothing, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true))
      end
      bnames = [c.name for c in generateN1Branches(net)]
      bytes = s -> Vector{UInt8}(codeunits(s))
      up = (fname, data) -> Sparlectra.handle_contingency_weights_upload(Dict{String,Any}("casefile" => "warmup_casePST.m", "casefiles" => [Sparlectra.WebUICaseUpload(fname, data)]); output_root = root, case_directory = cases, operation_log = root)
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
      @test occursin("rejected", loc(up("w.csv", zeros(UInt8, Sparlectra.WEBUI_CONTINGENCY_WEIGHTS_MAX_BYTES + 1))))
      @test occursin("rejected", loc(Sparlectra.handle_contingency_weights_upload(Dict{String,Any}("casefile" => "../evil", "casefiles" => [Sparlectra.WebUICaseUpload("w.csv", UInt8[])]); output_root = root, case_directory = cases, operation_log = root)))
      # uploading again replaces the file and says so
      @test occursin("replaced", loc(up("w.csv", bytes("name;weight\n$(bnames[1]);2.0\n"))))

      # the editor page seeds the case's real element names plus a raw-CSV editor
      page = redirect_stdout(devnull) do
        Sparlectra.handle_contingency_weights_page(Dict{String,Any}("case" => "warmup_casePST.m"); output_root = root, case_directory = cases, operation_log = root)
      end
      body = String(page.body)
      @test page.status == 200
      @test occursin(bnames[1], body)
      @test occursin("Raw CSV", body)

      # saving from the seeded table omits rows left at exactly 1.0
      Sparlectra.handle_contingency_weights_save(Dict{String,Any}("casefile" => "warmup_casePST.m", "element" => [bnames[1], bnames[2]], "weight" => ["2.5", "1.0"]); output_root = root, case_directory = cases, operation_log = root)
      saved = read(wf, String)
      @test occursin(bnames[1], saved)
      @test !occursin(bnames[2], saved)

      # download serves the stored file as an attachment
      dl = Sparlectra.handle_contingency_weights_download(Dict{String,Any}("case" => "warmup_casePST.m"); output_root = root, case_directory = cases)
      @test dl.status == 200
      @test any(k == "Content-Disposition" for (k, _) in dl.headers)
      @test !isempty(dl.body)

      # reset deletes the weight file
      @test Sparlectra.handle_contingency_weights_reset(Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root = root, case_directory = cases, operation_log = root).status == 303
      @test !isfile(wf)

      # deleting the case cascades to its weight file
      touch(wf)
      Sparlectra.handle_powerflow_case_delete(Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root = root, case_directory = cases, operation_log = root)
      @test !isfile(joinpath(cases, "warmup_casePST.m"))
      @test !isfile(wf)
    end

    @testset "state estimation page, measurement upload, and chain (SE phase 5)" begin
      root = mktempdir()
      cases = joinpath(root, "cases")
      mkpath(cases)
      case_path = joinpath(cases, "warmup_casePST.m")
      cp(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m"), case_path)
      rt = (; case_directory = cases, config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log = Sparlectra.webui_operation_log_path(root), startup_config_error = nothing, runner = Sparlectra.start_powerflow_run)
      bytes = s -> Vector{UInt8}(codeunits(s))

      # upload classification: v1 CSV -> measurement_set; CSV without the
      # version comment -> unknown, retained but never offered
      v1 = "# sparlectra-measurements v1\ntype,bus,from_bus,to_bus,link_nr,direction,value,sigma,active,id\n"
      up = Sparlectra.handle_powerflow_case_import(
        Dict{String,Any}("casefiles" => [Sparlectra.WebUICaseUpload("meas_v1.csv", bytes(v1)), Sparlectra.WebUICaseUpload("notes.csv", bytes("a,b\n1,2\n"))]);
        output_root = root, case_directory = cases, operation_log = root)
      @test up.status == 303
      loc = first(p for (k, p) in up.headers if k == "Location")
      @test occursin("measurement_set", Sparlectra._webui_urldecode(loc))
      @test isfile(joinpath(cases, "meas_v1.csv"))
      @test isfile(joinpath(cases, "notes.csv"))   # retained
      offered = Sparlectra._webui_measurement_options_in_directory(cases)
      @test "meas_v1.csv" in offered
      @test !("notes.csv" in offered)
      # neither CSV appears in the case selector
      @test !any(endswith(name, ".csv") for name in Sparlectra._webui_casefile_options_in_directory(cases))
      # traversal names still rejected by the shared upload checks
      bad = Sparlectra.handle_powerflow_case_import(Dict{String,Any}("casefiles" => [Sparlectra.WebUICaseUpload("../evil.csv", bytes(v1))]); output_root = root, case_directory = cases, operation_log = root)
      @test bad.status == 303   # redirect with the rejection message
      @test !isfile(joinpath(dirname(cases), "evil.csv"))

      # SE page opens headless, demo generator writes an offered v1 set,
      # and the run form carries no onsubmit button-disabling
      # stage 4A block 4: /stateestimation is a real redirect (a rendering
      # alias would be two ways to one surface, the divergence stage 4
      # removes); the SE section lives on the Runs page under its anchor
      page = Sparlectra.route_sparlectra_webui("GET", "/stateestimation", Dict{String,String}(); output_root = root, runtime = rt)
      @test page.status == 303
      @test endswith(Dict(page.headers)["Location"], "#state-estimation")
      main_page = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow", Dict{String,String}(); output_root = root, runtime = rt).body)
      @test !occursin("href=\"/stateestimation\"", main_page)
      @test occursin("id=\"state-estimation\"", main_page)
      @test occursin(">Runs<", main_page)
      @test !occursin(">Network analysis<", main_page)
      @test !occursin(">New run<", main_page)
      gen = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root = root, runtime = rt)
      @test gen.status == 303
      mfile = joinpath(cases, "warmup_casePST.measurements.csv")
      @test isfile(mfile)
      @test Sparlectra._webui_is_measurement_csv(mfile)
      page2 = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root = root, runtime = rt).body)
      @test occursin("warmup_casePST.measurements.csv", page2)
      @test occursin("Run state estimation", page2)
      @test !occursin("onsubmit", page2)

      # every SE input parameter carries a PF-style help link, and the
      # topics resolve (config-table rows, doc headings, and overrides)
      @test count("help-link", page2) >= 12
      for topic in ("state_estimation.tol", "state_estimation.flatstart", "state_estimation.max_iter", "state_estimation.robust", "state_estimation.update_shunts", "state_estimation.report_residual_correlation", "webui.se_max_eliminations", "webui.se_measurement_file", "webui.se_generator_noise", "webui.se_generator_gross_error", "webui.se_generator_tap_error", "webui.se_generator_sigmas")
        @test Sparlectra.handle_webui_help(topic).status == 200
      end
      @test occursin("MECHANICAL tap steps", String(Sparlectra.handle_webui_help("webui.se_generator_tap_error").body))

      # generator options: noise + gross error produce a valid, different set
      plain = read(mfile, String)
      gen2 = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "noise" => "true", "gross_error_k" => "8"); output_root = root, runtime = rt)
      @test gen2.status == 303
      @test occursin("bad%20data", Dict(gen2.headers)["Location"])
      @test Sparlectra._webui_is_measurement_csv(mfile)
      @test read(mfile, String) != plain   # noise + gross error changed values
      # regenerate the clean set for the runs below
      Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root = root, runtime = rt)
      @test read(mfile, String) == plain   # seeded generator is reproducible

      # per-quantity sigmas are PERCENT OF THE MEASURED VALUE with the
      # per-type floors (voltage-level independent); the currents checkbox
      # adds current-magnitude rows. Noise off, so value == truth and the
      # per-row sigma law is exactly reproducible.
      gen3 = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "sigma_u_pct" => "1.0", "include_currents" => "true", "sigma_i_pct" => "2.0", "sigma_p_pct" => "2.0", "sigma_q_pct" => "1.0"); output_root = root, runtime = rt)
      @test gen3.status == 303
      netchk = Sparlectra._import_sparlectra_net(case_path, nothing, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true))
      readMeasurementsCSV!(netchk; file = mfile)
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
      genbad = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "sigma_u_pct" => "-1"); output_root = root, runtime = rt)
      @test occursin("sigma%20U", Dict(genbad.headers)["Location"])

      # PMU current-angle rows via the sigma Ia field (absolute degrees)
      genia = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "include_currents" => "true", "sigma_ia_deg" => "0.1"); output_root = root, runtime = rt)
      @test genia.status == 303
      netia = Sparlectra._import_sparlectra_net(case_path, nothing, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true))
      readMeasurementsCSV!(netia; file = mfile)
      iarows = [m for m in netia.measurements if m.typ == Sparlectra.IaMeas]
      @test !isempty(iarows) && all(m.sigma == 0.1 for m in iarows)

      # tap deviation: measurements from a shifted-tap state differ from the
      # clean set and the message names the transformer branch
      plain2 = read(mfile, String)
      gentap = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "tap_error_steps" => "3"); output_root = root, runtime = rt)
      @test gentap.status == 303
      # on the PST warmup case the injected deviation is a Delta-u PHASE
      # step, and the message names the PST branch
      @test occursin("phase%20deviation", Dict(gentap.headers)["Location"])
      @test occursin("PST%20branch", Dict(gentap.headers)["Location"])
      @test read(mfile, String) != plain2
      @test Sparlectra._webui_is_measurement_csv(mfile)
      # the file records the tap positions the set was generated from as a
      # structured table (electrical/fixed/transferred step columns, the
      # deviated transformer carries its percentage), the SE page renders it
      # as a table and offers the download
      taptxt = read(mfile, String)
      @test occursin("# sparlectra-taps v1", taptxt)
      @test occursin("electrical_step,fixed_step,transferred_step,generation_deviation_steps", taptxt)
      @test occursin(",3.0", taptxt)   # the deviated transformer row
      pageInfo = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root = root, runtime = rt).body)
      @test occursin("Transformer taps at generation", pageInfo)
      @test occursin("<th>fixed_step</th>", pageInfo)
      @test occursin("Measured values in this set:", pageInfo)   # per-type row counts
      @test occursin("Vm ×", pageInfo)
      @test occursin("/stateestimation/measurements/download?file=warmup_casePST.measurements.csv", pageInfo)
      dlr = Sparlectra.route_sparlectra_webui("GET", "/stateestimation/measurements/download?file=warmup_casePST.measurements.csv", Dict{String,String}("file" => "warmup_casePST.measurements.csv"); output_root = root, runtime = rt)
      @test dlr.status == 200
      @test any(k == "Content-Disposition" for (k, _) in dlr.headers)
      dlbad = Sparlectra.route_sparlectra_webui("GET", "/stateestimation/measurements/download?file=../evil.csv", Dict{String,String}("file" => "../evil.csv"); output_root = root, runtime = rt)
      @test dlbad.status in (400, 404)
      # the commented file still parses and round-trips
      netc = Sparlectra._import_sparlectra_net(case_path, nothing, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true))
      rc = readMeasurementsCSV!(netc; file = mfile)
      @test rc.total > 0
      # out-of-range percentage is rejected
      gentapbad = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "tap_error_steps" => "35"); output_root = root, runtime = rt)
      @test occursin("between%20-16%20and%2016", Dict(gentapbad.headers)["Location"])
      # half steps are no longer settable (a tap changer has no half positions)
      genthalf = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "tap_error_steps" => "1.5"); output_root = root, runtime = rt)
      @test occursin("whole%20number", Dict(genthalf.headers)["Location"])
      # restore the default set for the runs below
      Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root = root, runtime = rt)
      @test read(mfile, String) == plain

      # since stage 4A the SE page without a query FOLLOWS the shared
      # selected-case memory (the maintainer's harmonization); a truly
      # fresh state (no memory under a fresh output root) still shows no
      # set info out of thin air
      Sparlectra.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m"; output_root = root, runtime = rt)
      pageRemembered = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow", Dict{String,String}(); output_root = root, runtime = rt).body)
      @test occursin("warmup_casePST.m", pageRemembered)
      pageFresh = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow", Dict{String,String}(); output_root = mktempdir(), runtime = rt).body)
      @test !occursin("Measurement set info", pageFresh)
      @test !occursin("Case binding:", pageFresh)

      # sticky generator inputs: the generate redirect carries the values
      # back and the re-rendered form keeps them instead of the defaults
      gsticky = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "sigma_u_pct" => "0.7", "noise" => "true"); output_root = root, runtime = rt)
      loc = Dict(gsticky.headers)["Location"]
      @test occursin("g_sigma_u_pct=0.7", loc)
      @test occursin("g_noise=true", loc)
      pageSticky = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m&g_sigma_u_pct=0.7&g_noise=true", Dict{String,String}(); output_root = root, runtime = rt).body)
      @test occursin("name=\"sigma_u_pct\" value=\"0.7\"", pageSticky)
      @test occursin("name=\"noise\" value=\"true\" checked", pageSticky)
      Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root = root, runtime = rt)

      # structured value editor: every measurement kind is editable in the
      # table; the update handler rewrites only value/sigma/active and a
      # single invalid entry rejects the whole save
      pageTab = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root = root, runtime = rt).body)
      @test occursin("Edit measurement values (table)", pageTab)
      @test occursin("/stateestimation/measurements/update-values", pageTab)
      rowsTab = Sparlectra._webui_measurement_set_rows(mfile)
      @test !isempty(rowsTab)
      @test any(r -> r.typ == "VmMeas", rowsTab) && any(r -> r.typ == "PflowMeas", rowsTab) && any(r -> r.typ == "PinjMeas", rowsTab)
      vrow = first(r for r in rowsTab if r.typ == "VmMeas")
      rup = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/measurements/update-values", Dict{String,Any}("file" => "warmup_casePST.measurements.csv", "case" => "warmup_casePST.m", "v_$(vrow.line)" => "1.0777", "s_$(vrow.line)" => vrow.sigma, "a_$(vrow.line)" => "true"); output_root = root, runtime = rt)
      @test occursin("updated%201", Dict(rup.headers)["Location"])
      @test occursin("1.0777", read(mfile, String))
      rbadv = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/measurements/update-values", Dict{String,Any}("file" => "warmup_casePST.measurements.csv", "case" => "warmup_casePST.m", "v_$(vrow.line)" => "1.05", "s_$(vrow.line)" => "-1", "a_$(vrow.line)" => "true"); output_root = root, runtime = rt)
      @test occursin("nothing%20was%20saved", Dict(rbadv.headers)["Location"])
      @test occursin("1.0777", read(mfile, String))   # rejected save left the file alone
      Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root = root, runtime = rt)

      # the run summary shows the chi-square band verdict with expected value
      fake = Dict("metadata" => Dict{String,Any}("run_mode" => "se", "se_observability_quality" => "good", "se_iterations" => 4, "se_objective" => 12.5, "se_dof" => 14, "se_band_reason" => "ok", "se_j_within_3sigma" => true))
      stxt = Sparlectra._webui_se_summary(fake)
      # J/dof leads: a bare J grows with the row count, so the same healthy set
      # reads as an alarm on a larger case (J = 104 at dof 95 is a ratio of 1.1)
      @test occursin("J/dof = 0.89", stxt) && occursin("dof = 14", stxt)
      @test findfirst("J/dof", stxt).start < findfirst("(J = ", stxt).start
      @test occursin("within 3", stxt)
      # a doubled set is named where the J is shown, or the page shows an
      # alarming number with no cause
      fake["metadata"]["se_duplicate_rows"] = 82
      @test occursin("82 measurement(s) repeat an already measured quantity", Sparlectra._webui_se_summary(fake))
      delete!(fake["metadata"], "se_duplicate_rows")
      @test !occursin("repeat an already measured", Sparlectra._webui_se_summary(fake))
      fake["metadata"]["se_j_within_3sigma"] = false
      fake["metadata"]["se_band_reason"] = "high"
      @test occursin("OUTSIDE", Sparlectra._webui_se_summary(fake))
      # :low reads as what it is (sigmas overstate the errors), never as an
      # alarm: the case57 confusion where "OUTSIDE (low)" was read as J too big
      fake["metadata"]["se_band_reason"] = "low"
      slow = Sparlectra._webui_se_summary(fake)
      @test occursin("far BELOW", slow)
      @test occursin("not an alarm", slow)
      @test !occursin("OUTSIDE", slow)
      fake["metadata"]["se_band_reason"] = "high"

      # topology panel renders findings and the explicit hypothesis button
      fakeT = Dict("run_id" => "t1", "success" => true, "metadata" => Dict{String,Any}("run_mode" => "se", "se_topology_findings" => [Dict{String,Any}("stage" => "precheck", "kind" => "open_element_with_flow", "location" => "branch 2 (A-B, open)", "evidence" => "Pflow x at 12 sigma", "severity" => "strong")], "se_topology_station_findings" => [Dict{String,Any}("location" => "B2 (+1 linked)", "evidence" => "Pinj_B2 (|rn| 9.1)", "notes" => ["precheck_agreement"])]))
      tsecT = Sparlectra._webui_se_topology_section(fakeT)
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
      stap = Sparlectra._webui_se_summary(fake)
      @test occursin("tap estimation: 2 transformer(s)", stap)
      @test occursin("before fixation", stap)
      @test !occursin("off-grid tap residual", stap)
      fake["metadata"]["se_tap_offgrid_residual"] = true
      @test occursin("NOT from bad data", Sparlectra._webui_se_summary(fake))
      fake["metadata"]["se_tap_offgrid_residual"] = false
      fake["metadata"]["se_tap_estimates"] = [
        Dict{String,Any}("branch" => 3, "name" => "T1", "mrid" => "", "mode" => "ratio", "electrical_step" => 2.03, "fixed_step" => 2, "electrical_shift_step" => 0.0, "fixed_shift_step" => 0, "out_of_range" => false, "fixed" => true),
        Dict{String,Any}("branch" => 4, "name" => "T2", "mrid" => "", "mode" => "pst", "electrical_step" => 0.0, "fixed_step" => 0, "electrical_shift_step" => -0.98, "fixed_shift_step" => -1, "out_of_range" => false, "fixed" => true),
      ]
      tsec = Sparlectra._webui_se_tap_section(fake)
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
      @test Sparlectra._webui_measurement_set_case(mfile) == "warmup_casePST.m"
      other = joinpath(cases, "case9.measurements.csv")
      cp(mfile, other)
      page3 = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root = root, runtime = rt).body)
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
      page9 = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow?casefile=sp_case5.scf.json", Dict{String,String}(); output_root = root, runtime = rt).body)
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
      pageEd = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root = root, runtime = rt).body)
      @test occursin("Edit measurement file (inline)", pageEd)
      @test occursin("/stateestimation/measurements/save", pageEd)
      content0 = read(mfile, String)
      edited = replace(content0, "# case: warmup_casePST.m" => "# case: warmup_casePST.m\n# note: edited inline"; count = 1)
      rsave = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/measurements/save", Dict{String,Any}("file" => "warmup_casePST.measurements.csv", "case" => "warmup_casePST.m", "content" => edited); output_root = root, runtime = rt)
      @test rsave.status == 303
      @test occursin("saved", Dict(rsave.headers)["Location"])
      @test occursin("# note: edited inline", read(mfile, String))
      rbad1 = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/measurements/save", Dict{String,Any}("file" => "../evil.csv", "case" => "warmup_casePST.m", "content" => edited); output_root = root, runtime = rt)
      @test occursin("invalid", Dict(rbad1.headers)["Location"])
      rbad2 = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/measurements/save", Dict{String,Any}("file" => "warmup_casePST.measurements.csv", "case" => "warmup_casePST.m", "content" => "a,b\n1,2\n"); output_root = root, runtime = rt)
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
      runs = Sparlectra.list_powerflow_runs(root)
      row = only(r for r in runs if string(get(r, "run_id", "")) == id1)
      @test get(row, "run_mode", "") == "se"
      # the SE result page carries the chain action
      resPage = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow/result/$(id1)", Dict{String,String}(); output_root = root, runtime = rt).body)
      @test occursin("Run power flow from this estimate", resPage)
      @test occursin("se_start_run_id", resPage)

      # topology validation on the service run: the precheck RESULT is
      # logged on every run (the clean case says so explicitly), and the
      # result page offers the explicit hypothesis-test button
      @test occursin("topology precheck: no findings", read(joinpath(root, id1, "run.log"), String))
      @test occursin("Test topology hypotheses", resPage)
      rhyp = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/topology-hypotheses", Dict{String,Any}("run_id" => id1); output_root = root, runtime = rt)
      @test rhyp.status == 303
      @test isfile(joinpath(root, id1, "topology_hypotheses.md"))
      @test occursin("Recommendations only: NOTHING has been switched", read(joinpath(root, id1, "topology_hypotheses.md"), String))
      @test occursin("topology hypothesis test:", read(joinpath(root, id1, "run.log"), String))
      resPage2 = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow/result/$(id1)", Dict{String,String}(); output_root = root, runtime = rt).body)
      @test occursin("Topology validation (advisory)", resPage2)

      # tap-estimation SE run: a set generated with a tap deviation
      # disagrees with the model around one transformer; releasing the taps
      # absorbs the discrepancy and the fixation lands on a mechanical step
      Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "tap_error_steps" => "3"); output_root = root, runtime = rt)
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
      @test any(l -> (f = split(l, ","); any(parse(Int, f[i]) != 0 for i in step_cols)), taplines[2:end])
      tapPage = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow/result/$(idtap)", Dict{String,String}(); output_root = root, runtime = rt).body)
      @test occursin("Transformer tap estimates", tapPage)
      @test occursin("tap estimation: 1 transformer(s)", tapPage)

      # A set that DOCUMENTS its tap deviations releases exactly those
      # transformers by itself, without the user asking for it. As a mere
      # hint this produced a J nobody could explain, and the hint only
      # appeared after a converged run (maintainer 2026-09-04: "if
      # transformers were changed, that option has to be on by itself").
      rhint = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
      @test rhint["status"] == "succeeded"
      @test rhint["metadata"]["se_set_tap_deviation"] == true
      @test !isempty(rhint["metadata"]["se_auto_released_taps"])
      @test occursin("released automatically", rhint["message"])
      @test occursin("tap estimation released automatically", read(joinpath(root, rhint["run_id"], "run.log"), String))

      # bad data lands findable: gross error in the set -> se_bad_data.csv
      # with the measurement, its location, and the eliminated flag
      Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gross_error_k" => "10"); output_root = root, runtime = rt)
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
      bdPage = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow/result/$(rbd["run_id"])", Dict{String,String}(); output_root = root, runtime = rt).body)
      @test occursin("Bad data (se_bad_data.csv)", bdPage)
      @test occursin("/powerflow/artifact/$(rbd["run_id"])/se_bad_data.csv", bdPage)

      # restore the clean default set for anything below
      Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root = root, runtime = rt)

      # --- measurement generator v2: truth state, flow ends, passive nodes,
      # delta comments, and the bad-data threshold surface of the run form
      genpage = String(Sparlectra.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.m", Dict{String,String}(); output_root = root, runtime = rt).body)
      for needle in ("name=\"gen_truth_source\"", "name=\"gen_truth_run_id\"", "name=\"gen_flow_ends\"", "name=\"gen_passive_sigma\"", "name=\"gen_passive_as_zi\"", "name=\"se_robust_mode\"", "name=\"se_k_eliminate\"", "name=\"se_k_suppress\"", "name=\"se_suppression_sigma\"", "se-threshold-warning", "gen-truth-source")
        @test occursin(needle, genpage)
      end

      # one balance-aware flow end: deterministic (two generates produce the
      # identical file), one flow group per branch, choice documented
      Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_flow_ends" => "one_balance_aware"); output_root = root, runtime = rt)
      one1 = read(mfile, String)
      Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_flow_ends" => "one_balance_aware"); output_root = root, runtime = rt)
      one2 = read(mfile, String)
      @test one1 == one2
      @test occursin("# seed: 42", one1)
      # generating persists the generator options in the form block of the
      # case configuration file so a case reload restores them
      @test occursin("gen_flow_ends: one_balance_aware", read(joinpath(cases, "warmup_casePST.config.yaml"), String))
      @test occursin("# flow_ends: one_balance_aware", one1)
      @test occursin("# flow_end,", one1)
      @test occursin("# truth_value,Vm_", one1)
      # the full default set measures injections at every bus, so every
      # branch keeps its from end and no to-direction flow rows remain
      @test all(l -> !(startswith(l, "PflowMeas") && length(split(l, ",")) >= 7 && split(l, ",")[7] == "to"), split(one1, "\n"))
      rone = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile))
      @test rone["status"] == "succeeded"

      # passive nodes as protected zero-injection constraints: ZI rows
      # written, no duplicate plain injection rows, elimination stays away
      Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_passive_as_zi" => "true"); output_root = root, runtime = rt)
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
      Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_truth_source" => "from_run", "gen_truth_run_id" => id1); output_root = root, runtime = rt)
      fr1 = read(mfile, String)
      @test occursin("# truth: run $(id1) (se,", fr1)
      sestate = readlines(joinpath(root, id1, "se_state.csv"))
      strow = only([l for l in sestate if startswith(l, "2,")])
      st_vm = split(strow, ",")[2]
      vmrow = only([l for l in split(fr1, "\n") if startswith(l, "VmMeas,2,")])
      @test split(vmrow, ",")[8] == st_vm
      # rejections: tap deviation locked, unknown run id named
      rrej = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_truth_source" => "from_run", "gen_truth_run_id" => id1, "tap_error_steps" => "2"); output_root = root, runtime = rt)
      @test occursin("requires%20truth%20state", Dict(rrej.headers)["Location"])
      rrej2 = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m", "gen_truth_source" => "from_run", "gen_truth_run_id" => "nope"); output_root = root, runtime = rt)
      @test occursin("not%20found%20in%20the%20run%20history", Dict(rrej2.headers)["Location"])

      # restore the clean default set and exercise the threshold surface:
      # the staged service path equals the legacy Bool bitwise, the delta
      # artifact exists for generated sets, invalid modes reject
      Sparlectra.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.m"); output_root = root, runtime = rt)
      rst = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => mfile, "se_robust_mode" => "staged", "se_robust_k1" => 3.0, "se_robust_k2" => 6.0))
      @test rst["status"] == "succeeded"
      @test rst["metadata"]["se_robust_mode"] == "staged"
      # the request builder records the SE options as sidecar-persistable
      # settings, so the browser flow's "save settings" keeps them
      reqrec = Sparlectra._webui_request_settings_for_profile(Sparlectra.powerflow_webui_request(Dict{String,Any}("se_mode" => "true", "casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "measurement_file" => mfile, "se_robust_mode" => "staged", "se_robust_k2" => "6.0"); default_output_root = root))
      @test reqrec["se_robust_mode"] == "staged"
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
      gpst = Sparlectra._se_generate_measurement_set(pstcase, pstout; noise = true, gross_k = 0.0, tap_steps = 2.0, include_i = false, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0, sigma_ia_deg = 0.0)
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
      genmulti(out, seed) = Sparlectra._se_generate_measurement_set(d14, out; noise = false, gross_k = 10.0, gross_count = 3, tap_steps = 2.0, tap_count = 2, include_i = false, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0, sigma_ia_deg = 0.0, seed = seed)
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
      gmp = Sparlectra._se_generate_measurement_set(pstcase, gopst; noise = false, gross_k = 0.0, tap_steps = 3.0, tap_count = 2, include_i = false, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0, sigma_ia_deg = 0.0, seed = 42)
      @test occursin("PST branch 8", gmp.tap_note)
      @test !occursin("branch 7", gmp.tap_note)
      @test occursin("limited to 1 eligible transformer(s), max 2 requested", gmp.tap_note)
      pstlines = readlines(gopst)
      @test any(l -> startswith(l, "# 7,") && endswith(l, ",0.0"), pstlines)
      @test any(l -> startswith(l, "# 8,") && endswith(l, ",3.0"), pstlines)
      @test_throws ArgumentError Sparlectra._se_generate_measurement_set(d14, go3; noise = false, gross_k = 10.0, gross_count = 0, tap_steps = 0.0, include_i = false, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0, sigma_ia_deg = 0.0)
      @test_throws ArgumentError Sparlectra._se_generate_measurement_set(d14, go3; noise = false, gross_k = 0.0, tap_steps = 1.0, tap_count = 0, include_i = false, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0, sigma_ia_deg = 0.0)

      # J_active: replacement suppression removes corrupted rows from the
      # STATE; the reported pair (honest J with original sigmas, J_active
      # over the trusted rows) makes that visible, and the band verdict
      # stays on the honest J (eliminations off so the rows STAY suppressed)
      gja = joinpath(root, "gen_jactive.csv")
      Sparlectra._se_generate_measurement_set(case_path, gja; noise = false, gross_k = 12.0, gross_count = 2, tap_steps = 0.0, include_i = false, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0, sigma_ia_deg = 0.0, seed = 42)
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
      Sparlectra._se_generate_measurement_set(case_path, gjb; noise = true, gross_k = 10.0, gross_count = 1, tap_steps = 0.0, include_i = false, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0, sigma_ia_deg = 0.0, seed = 42)
      rje = start_powerflow_run(Dict{String,Any}("casefile" => case_path, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "se_mode" => true, "measurement_file" => gjb, "se_robust_mode" => "replacement", "se_k_suppress" => 4.0))
      @test rje["status"] == "succeeded"
      @test rje["metadata"]["se_eliminations"] == 1
      @test rje["metadata"]["se_band_reason"] == "ok"
      @test rje["metadata"]["se_objective"] < 2.0 * rje["metadata"]["se_dof"]

      # reset-settings deletes the per-case sidecar profile; a second
      # reset reports that the defaults are already active
      Sparlectra._webui_merge_case_settings!(root, case_path, Dict{String,Any}("gen_seed" => 99); case_directory = dirname(case_path))
      spath = Sparlectra._webui_case_settings_path(root, case_path; case_directory = dirname(case_path))
      @test isfile(spath)
      rrst = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/reset-settings", Dict{String,Any}("casefile" => basename(case_path)); output_root = root, runtime = rt)
      @test !isfile(spath)
      @test occursin("deleted", string(rrst))
      rrst2 = Sparlectra.route_sparlectra_webui("POST", "/stateestimation/reset-settings", Dict{String,Any}("casefile" => basename(case_path)); output_root = root, runtime = rt)
      @test occursin("already", string(rrst2))
      @test occursin("reset-settings", Sparlectra.render_se_form(; cases = [basename(case_path)], selected_case = basename(case_path)))

      # noise defaults ON: a fresh form (no sidecar, no stickies) checks
      # the box (a noise-free set puts J near 0 instead of near dof, which
      # reads like a broken statistic); an explicit false stays unchecked
      @test occursin("name=\"noise\" value=\"true\" checked", Sparlectra.render_se_form(; cases = ["case14.m"], selected_case = "case14.m"))
      @test !occursin("name=\"noise\" value=\"true\" checked", Sparlectra.render_se_form(; cases = ["case14.m"], selected_case = "case14.m", gen_values = Dict{String,String}("noise" => "false")))
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
      for (l3, l4) in zip(lines3, lines4)
        f3 = split(l3, ";")
        f4 = split(l4, ";")
        @test length(f3) == length(f4)
        for (a, b) in zip(f3, f4)
          na = tryparse(Float64, a)
          nb = tryparse(Float64, b)
          if na !== nothing && nb !== nothing
            @test (isnan(na) && isnan(nb)) || isapprox(na, nb; atol = 1e-9)
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
    end

    @testset "saved case settings outrank the configuration on the form" begin
      # No mtime logic anywhere (D5): the case levels always win over the
      # general configuration for the fields they set, however new the
      # configuration file is.
      dir = mktempdir()
      cfg = joinpath(dir, "conf.yaml")
      write(cfg, "config_version: 1\nscope: general\npower_flow:\n  max_iter: 99\n")
      profile = Dict{String,Any}("power_flow_max_iter" => 55, "power_flow_autodamp_min" => 0.07, "_profile_path" => joinpath(dir, "case57.config.yaml"))
      v = Sparlectra.webui_form_state(selected_config_file = cfg, sidecar_profile = profile)
      @test v["power_flow_max_iter"] == 55
      @test v["power_flow_autodamp_min"] == 0.07
      @test !haskey(v, "_config_newer_than_profile")
      # a field the case does not set comes from the configuration
      v_cfg_only = Sparlectra.webui_form_state(selected_config_file = cfg)
      @test v_cfg_only["power_flow_max_iter"] == 99
    end

    @testset "sysimage launcher decision" begin
      # The launcher decides whether the Web UI starts from the image or
      # compiles on first use. It runs on plain Base BEFORE the package is
      # loaded (tools/sysimage_launcher.jl), so it is tested through the
      # module, not through Sparlectra. The last case is the expensive one:
      # the metadata pins the Julia version and the Manifest, i.e. the
      # DEPENDENCIES, so an image built before a src/ edit still looked
      # fresh and the Web UI silently served old code.
      launcher = Module(:LauncherUnderTest)
      Base.include(launcher, joinpath(Sparlectra.SPARLECTRA_ROOT, "tools", "sysimage_launcher.jl"))
      SL = getfield(launcher, :SysimageLauncher)
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
      end
      # no terminal in the test process, so the question answers itself with
      # the default: an unattended start ends up WITH an image
      @test sl(:ask_build, "never shown") == true
      @test slval(:REBUILD_FLAG) == "--rebuild-sysimage"
      @test slval(:NO_IMAGE_FLAG) == "--no-sysimage"
      # platform path contract, still owned by the package for everything
      # that runs after the load
      @test endswith(sl(:sysimage_path), Sparlectra.webui_sysimage_ext())
      @test sl(:sysimage_path) == Sparlectra.webui_sysimage_path()
    end

    @testset "case-scan memo cache" begin
      mktempdir() do dir
        p = joinpath(dir, "SCAN.DAT")
        write(p, _dtf_network_fixture())
        @test Sparlectra._webui_classify_dat_content_cached(p) === :dtf_network_case
        # second call answers from the (path, mtime, size) memo
        @test Sparlectra._webui_classify_dat_content_cached(p) === :dtf_network_case
        # a replaced file re-classifies (mtime resolution can be one second)
        sleep(1.1)
        write(p, _dtf_outage_fixture())
        @test Sparlectra._webui_classify_dat_content_cached(p) === :dtf_outage_file
        # the selector honors the refreshed classification
        @test !Sparlectra._webui_is_user_selectable_case(p)
      end
    end

    @testset "warmup prefix hides only .jl workloads" begin
      # the reserved warmup workloads are .jl; a MATPOWER case carrying the
      # prefix (warmup_casePST.m) must stay selectable (regression: the demo
      # case vanished from both selectors after its rename)
      @test !Sparlectra._webui_is_user_selectable_case("warmup_case3.jl")
      @test !Sparlectra._webui_is_user_selectable_case("warmup_case118.jl")
      @test Sparlectra._webui_is_user_selectable_case("warmup_casePST.m")
      @test !Sparlectra._webui_is_user_selectable_case("warmup_casePST.sparlectra-webui.yaml")
    end

    @testset "buildSysimage one-call API (dry run)" begin
      # the exported one-liner plans against the package project in a child
      # process; dry run must resolve the target paths without building
      r = Sparlectra.buildSysimage(dry_run = true, quiet = true)
      @test r.built == false
      @test endswith(r.sysimage_path, Sparlectra.webui_sysimage_ext()) || occursin("sparlectra.", basename(r.sysimage_path))
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
          src = read(joinpath(d, "src", "SparlectraApp.jl"), String)
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
      # local documentation viewer (maintainer 2026-09-04): LaTeX must not
      # reach the browser as raw markup, and a thousand-line reference page
      # needs a section index instead of a scrollbar
      math_html = Sparlectra.render_webui_markdown(raw"Text with $\Omega_{ii} = w_i \cdot \sigma^2$ inline.")
      @test !occursin("\\Omega", math_html)
      @test !occursin("&#36;", math_html)
      @test occursin("Ω", math_html)
      @test occursin("σ", math_html)
      block_html = Sparlectra.render_webui_markdown("a\n\n```math\n\\frac{P_{se}}{V}\n```\n")
      @test !occursin("frac", block_html)
      @test occursin("/", block_html)
      long_md = read(joinpath(pkgdir(Sparlectra), "docs", "src", "state_estimation.md"), String)
      toc = Sparlectra._webui_doc_page_toc(long_md)
      @test occursin("On this page", toc)
      @test length(collect(eachmatch(r"<li>", toc))) > 8
      # a short page gets no index (it would be noise)
      @test Sparlectra._webui_doc_page_toc("# T\n\n## One\n\ntext\n") == ""
      # the save-target explanations are small print, not label text
      stg = Sparlectra.render_settings_page(output_root = mktempdir(), selected_casefile = "sp_case14.scf.json")
      @test occursin("settings-target-hint", stg)
      @test !occursin("value=\"this_case\" checked>this case (<code>sp_case14.scf.json</code>): case-scope", stg)

      # header start-flavor: native in the test session, :app when the
      # standalone executable stamped its build time into the environment
      @test Sparlectra.webui_runtime_flavor().kind in (:native, :sysimage)
      withenv("SPARLECTRA_APP_BUILT" => "2026-08-29T11:11:11") do
        fa = Sparlectra.webui_runtime_flavor()
        @test fa.kind === :app
        @test fa.built == "2026-08-29T11:11:11"
      end
    end

    @testset "sysimage build script dry run parses" begin
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
    end

    @testset "SBOM script dry run parses" begin
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
    end
  end
end

run_webui_tests() = run_webui_fast_tests()

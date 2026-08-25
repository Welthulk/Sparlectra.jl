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
#          the fast-start sysimage metadata contract
using Sparlectra
using Test

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
      form_html = Sparlectra.render_powerflow_form(output_root = mktempdir())
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
      with_sidecar = Sparlectra.render_powerflow_form(output_root = mktempdir(), selected_casefile = "c.m", case_profile = Dict{String,Any}("power_flow_solver" => "dc", "_profile_path" => prof))
      @test occursin("Reset saved settings for this case", with_sidecar)
      @test occursin("/powerflow/case-settings/reset", with_sidecar)
      without = Sparlectra.render_powerflow_form(output_root = mktempdir())
      @test !occursin("Reset saved settings for this case", without)
      # switching the case reloads the form server-side; the wait must be visible
      @test occursin("case-loading-banner", without)

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

    @testset "contingency weights editor and storage (#331 Phase 5 follow-up)" begin
      dir = mktempdir()
      cases = joinpath(dir, "cases")
      mkpath(cases)
      root = joinpath(dir, "runs")
      mkpath(root)
      real14 = joinpath(pkgdir(Sparlectra), "test", "testdata", "mpower", "case14.m")
      isfile(real14) || (real14 = joinpath(pkgdir(Sparlectra), "data", "mpower", "case14.m"))
      cp(real14, joinpath(cases, "case14.m"))
      wf = Sparlectra._webui_case_weights_path("case14.m"; case_directory = cases)

      # the weight file lives next to the case as <stem>.contingency-weights.csv
      @test basename(wf) == "case14.contingency-weights.csv"
      @test dirname(wf) == normpath(cases)

      # list exclusion: a weights file must not be offered as a selectable case
      touch(wf)
      opts = Sparlectra._webui_casefile_options_in_directory(cases)
      @test "case14.m" in opts
      @test !any(occursin("contingency-weights", o) for o in opts)
      rm(wf)

      # real element names for the fixtures
      net = redirect_stdout(devnull) do
        Sparlectra._import_sparlectra_net(joinpath(cases, "case14.m"), nothing, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true))
      end
      bnames = [c.name for c in generateN1Branches(net)]
      bytes = s -> Vector{UInt8}(codeunits(s))
      up = (fname, data) -> Sparlectra.handle_contingency_weights_upload(Dict{String,Any}("casefile" => "case14.m", "casefiles" => [Sparlectra.WebUICaseUpload(fname, data)]); output_root = root, case_directory = cases, operation_log = root)
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
        Sparlectra.handle_contingency_weights_page(Dict{String,Any}("case" => "case14.m"); output_root = root, case_directory = cases, operation_log = root)
      end
      body = String(page.body)
      @test page.status == 200
      @test occursin(bnames[1], body)
      @test occursin("Raw CSV", body)

      # saving from the seeded table omits rows left at exactly 1.0
      Sparlectra.handle_contingency_weights_save(Dict{String,Any}("casefile" => "case14.m", "element" => [bnames[1], bnames[2]], "weight" => ["2.5", "1.0"]); output_root = root, case_directory = cases, operation_log = root)
      saved = read(wf, String)
      @test occursin(bnames[1], saved)
      @test !occursin(bnames[2], saved)

      # download serves the stored file as an attachment
      dl = Sparlectra.handle_contingency_weights_download(Dict{String,Any}("case" => "case14.m"); output_root = root, case_directory = cases)
      @test dl.status == 200
      @test any(k == "Content-Disposition" for (k, _) in dl.headers)
      @test !isempty(dl.body)

      # reset deletes the weight file
      @test Sparlectra.handle_contingency_weights_reset(Dict{String,Any}("casefile" => "case14.m"); output_root = root, case_directory = cases, operation_log = root).status == 303
      @test !isfile(wf)

      # deleting the case cascades to its weight file
      touch(wf)
      Sparlectra.handle_powerflow_case_delete(Dict{String,Any}("casefile" => "case14.m"); output_root = root, case_directory = cases, operation_log = root)
      @test !isfile(joinpath(cases, "case14.m"))
      @test !isfile(wf)
    end

    @testset "config edits win over older case settings" begin
      dir = mktempdir()
      cfg = joinpath(dir, "conf.yaml")
      prof = joinpath(dir, "profile.yaml")
      write(prof, "placeholder")
      sleep(1.1)
      write(cfg, "power_flow:\n  max_iter: 99\n")
      sidecar = Dict{String,Any}("power_flow_max_iter" => 55, "power_flow_autodamp_min" => 0.07, "_profile_path" => prof)
      v = Sparlectra.webui_form_state(selected_config_file = cfg, sidecar_profile = sidecar)
      # config file is newer: its keys win; fields the YAML does not set
      # keep their saved case value
      @test v["power_flow_max_iter"] == 99
      @test v["power_flow_autodamp_min"] == 0.07
      @test v["_config_newer_than_profile"] === true
      # the replaced values are named explicitly (config key, saved, config):
      # a silent flip cost a converged case its start_values before
      @test v["_config_over_profile_details"] == [("power_flow.max_iter", "55", "99")]
      rendered = Sparlectra.render_powerflow_form(output_root = mktempdir(), selected_config_file = cfg, case_profile = sidecar)
      @test occursin("These saved values were replaced", rendered)
      @test occursin("power_flow.max_iter", rendered)
      # sidecar newer again: saved case settings win as before
      sleep(1.1)
      write(prof, "placeholder2")
      v2 = Sparlectra.webui_form_state(selected_config_file = cfg, sidecar_profile = sidecar)
      @test v2["power_flow_max_iter"] == 55
      @test !haskey(v2, "_config_newer_than_profile")

      # Regression: a NEWER config that only repeats the shipped template
      # defaults carries no user intent (startup/migration rewrites touch the
      # file) and must NOT discard the saved case setup — case_SyntheticUSA
      # lost its solver settings this way and diverged.
      template_only = joinpath(dir, "template_only.yaml")
      template_defaults = Sparlectra._webui_config_field_values(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH)
      write(template_only, "power_flow:\n  qlimits:\n    enabled: $(get(template_defaults, "power_flow_qlimits_enabled", true))\n")
      sleep(1.1)
      touch(template_only)
      case_setup = Dict{String,Any}("power_flow_qlimits_enabled" => false, "power_flow_merit_enabled" => true, "_profile_path" => prof)
      v3 = Sparlectra.webui_form_state(selected_config_file = template_only, sidecar_profile = case_setup)
      @test v3["power_flow_qlimits_enabled"] === false
      @test v3["power_flow_merit_enabled"] === true
      @test !haskey(v3, "_config_newer_than_profile")
    end

    @testset "fast-start sysimage metadata contract" begin
      mktempdir() do tmp
        root = joinpath(tmp, "runs")
        mkpath(root)
        manifest = joinpath(tmp, "Manifest.toml")
        write(manifest, "# manifest fixture\n[[deps.Example]]\nuuid = \"0\"\n")
        img_path = Sparlectra.webui_sysimage_path(root)
        meta_path = Sparlectra.webui_sysimage_meta_path(root)
        mkpath(dirname(img_path))
        write(img_path, "not a real image, size is all that matters here")

        # fresh metadata next to an existing image: valid, no reasons
        Sparlectra.write_sysimage_meta(meta_path; manifest_path = manifest)
        st = Sparlectra.sysimage_status(output_root = root, manifest_path = manifest)
        @test st.present
        @test st.meta_present
        @test st.valid
        @test isempty(st.reasons)
        @test st.size_bytes > 0
        @test st.julia_version == string(VERSION)
        @test st.sparlectra_version == string(Sparlectra.SparlectraVersion)

        # wrong Julia version invalidates
        meta = Sparlectra.read_sysimage_meta(meta_path)
        meta["julia_version"] = "1.0.0"
        open(meta_path, "w") do io
          Sparlectra.TOML.print(io, meta)
        end
        st_julia = Sparlectra.sysimage_status(output_root = root, manifest_path = manifest)
        @test !st_julia.valid
        @test any(occursin("Julia version", r) for r in st_julia.reasons)

        # wrong manifest hash invalidates (package update case)
        Sparlectra.write_sysimage_meta(meta_path; manifest_path = manifest)
        write(manifest, "# changed manifest\n")
        st_sha = Sparlectra.sysimage_status(output_root = root, manifest_path = manifest)
        @test !st_sha.valid
        @test any(occursin("Manifest.toml changed", r) for r in st_sha.reasons)

        # missing image file invalidates
        Sparlectra.write_sysimage_meta(meta_path; manifest_path = manifest)
        rm(img_path)
        st_img = Sparlectra.sysimage_status(output_root = root, manifest_path = manifest)
        @test !st_img.present
        @test !st_img.valid
        @test any(occursin("sysimage file missing", r) for r in st_img.reasons)

        # missing metadata invalidates
        write(img_path, "image again")
        rm(meta_path)
        st_meta = Sparlectra.sysimage_status(output_root = root, manifest_path = manifest)
        @test st_meta.present
        @test !st_meta.meta_present
        @test !st_meta.valid
      end
    end

    @testset "fast-start build job and routes" begin
      mktempdir() do tmp
        root = joinpath(tmp, "runs")
        mkpath(root)
        # the suite must never run a real PackageCompiler build; the job
        # machinery is exercised with substitute commands
        @test Sparlectra.sysimage_build_state().status in (:idle, :completed, :failed)
        started = Sparlectra.start_sysimage_build!(root; _test_command = Cmd(["sleep", "1"]))
        @test started.ok
        @test Sparlectra.sysimage_build_state().running
        # one build at a time: a second request is rejected, not queued
        rejected = Sparlectra.start_sysimage_build!(root; _test_command = Cmd(["sleep", "1"]))
        @test !rejected.ok
        @test occursin("already running", rejected.message)
        # POST route while running: redirect plus rejected operation event
        response = Sparlectra.route_sparlectra_webui("POST", "/webui/fast-start/build"; output_root = root)
        @test response.status == 303
        @test Dict(response.headers)["Location"] == "/webui/fast-start"
        # the status page renders the running state with a visible spinner
        # and the elapsed time
        page = String(Sparlectra.route_sparlectra_webui("GET", "/webui/fast-start"; output_root = root).body)
        @test occursin("Fast start", page)
        @test occursin("build running", page)
        @test occursin("disabled", page)
        @test occursin("warmup-spinner", page)
        @test occursin("running for", page)
        @test Sparlectra.sysimage_build_state().elapsed_seconds >= 0.0
        # wait for the fake build; the job finishes as completed
        ok = timedwait(() -> !Sparlectra.sysimage_build_state().running, 30.0) === :ok
        @test ok
        @test Sparlectra.sysimage_build_state().status === :completed
        # failure path: a failing command ends as failed
        failed = Sparlectra.start_sysimage_build!(root; _test_command = Cmd(["false"]))
        @test failed.ok
        @test timedwait(() -> !Sparlectra.sysimage_build_state().running, 30.0) === :ok
        @test Sparlectra.sysimage_build_state().status === :failed
        # operation log carries the build events
        oplog = read(Sparlectra.webui_operation_log_path(root), String)
        @test occursin("sysimage_build_started", oplog)
        @test occursin("sysimage_build_finished", oplog)
        @test occursin("sysimage_build_failed", oplog)
        # idle page on a root without an image reports "not built"
        page_idle = String(Sparlectra.route_sparlectra_webui("GET", "/webui/fast-start"; output_root = root).body)
        @test occursin("not built", page_idle)
        # log page renders the captured command output artifact
        page_log = String(Sparlectra.route_sparlectra_webui("GET", "/webui/fast-start/log"; output_root = root).body)
        @test occursin("Fast-start build log", page_log)
        # leave the registry clean for other testsets
        lock(Sparlectra._SYSIMAGE_BUILD_LOCK) do
          Sparlectra._SYSIMAGE_BUILD_JOB[] = nothing
        end
      end
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

    @testset "fast-start sysimage detection" begin
      # the test process runs on the default system image, never on the
      # fast-start image
      @test Sparlectra.webui_sysimage_active() == false
    end

    @testset "fast-start build script dry run parses" begin
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

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
  end
end

run_webui_tests() = run_webui_fast_tests()

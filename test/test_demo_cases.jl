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

# file: test/test_demo_cases.jl
# purpose: the shipped demo cases (data/scf/sp_case*) are read through
#          import_case like any other SCF file; power flow, state estimation
#          on the measurements the file carries, and the declared controllers
#          must converge; sp_case6tap holds its three tap control targets.
#          The file checks (provenance, case config, size budget, bad-data
#          twin) live in test_repository_hygiene.jl.

const _DEMO_CASE_DIR = normpath(joinpath(dirname(@__DIR__), "data", "scf"))
const _DEMO_CASE_NAMES = ("sp_case5", "sp_case6tap", "sp_case14", "sp_case60", "sp_case188")
# the hand-written tap control case carries no measurement rows in its file;
# its measurement table next to it is read instead
const _DEMO_CASES_TABLE_ONLY = ("sp_case6tap",)
# scf_with_foreign_stamp for the legacy group of this file
isdefined(@__MODULE__, :scf_roundtrip_field_diffs) || include("test_scf_support.jl")
# the cases are machine-neutral: they run against the PACKAGE defaults and
# never resolve through this machine's configuration.yaml
_DEMO_REF_CONFIG() = Sparlectra.SparlectraConfig(Dict())

function run_demo_case_tests()
  # one outer set, so a failing case still runs the sets after it
  @testset "Shipped demo cases" begin (function ()
    for name in _DEMO_CASE_NAMES
      @testset "$(name)" begin (function ()
        scf_path = joinpath(_DEMO_CASE_DIR, "$(name).scf.json")

        # 1) power flow through the standard import path (start_state applies)
        net = Sparlectra.import_case(scf_path, _DEMO_REF_CONFIG()).net
        _, erg = runpf!(net; verbose = 0)
        @test erg == 0

        # 2) state estimation on the measurements the file itself carries;
        # the topology precheck advises on two of the shipped cases by design,
        # captured, anything else that warns fails
        meas_file = joinpath(_DEMO_CASE_DIR, "$(name).measurements.csv")
        @test isfile(meas_file)
        name in _DEMO_CASES_TABLE_ONLY && readMeasurementTable!(net; file = meas_file)
        se = run_with_expected_warnings(() -> runse!(net), (r"topology precheck reported",))
        @test se.converged

        # 3) the declared controllers of a shipped case must settle: a
        # comparison with a recorded run once froze a non-converging SSSC as
        # the reference. sp_case5 is the one case without a controller.
        ctrl_net = Sparlectra.import_case(scf_path, _DEMO_REF_CONFIG()).net
        controllers = Sparlectra.collect_outer_controllers(ctrl_net)
        @test isempty(controllers) == (name == "sp_case5")
        if !isempty(controllers)
          cres = Sparlectra.run_control!(ctrl_net; verbose = 0)
          @test cres.converged
        end
      end)() end
    end

    # sp_case6tap, the tap control case (diagram
    # data/scf/sp_case6tap_single_line.svg): its names, the ratio changer of
    # T2 that no controller drives, and the three targets of the case file
    # met inside their deadbands after the outer control loop
    @testset "sp_case6tap names, T2 changer and control targets" begin (function ()
      net = Sparlectra.import_case(joinpath(_DEMO_CASE_DIR, "sp_case6tap.scf.json"), _DEMO_REF_CONFIG()).net
      @test [n.comp.cName for n in net.nodeVec] == ["Nordwerk_220", "Suedfeld_220", "Kreuzau_110", "Mittelbach_110", "Ostheim_110", "Lindenau_20"]
      @test [b.comp.cName for b in net.branchVec] == ["Nordwerk_Suedfeld", "Kreuzau_Mittelbach", "Mittelbach_Ostheim", "Kreuzau_Ostheim", "T1_Nordwerk_Kreuzau", "T2_Suedfeld_Ostheim", "T3_Mittelbach_Lindenau"]
      # T2: declared ratio changer, +-12 steps of 1.25 percent around the
      # neutral ratio 1.0, at position 0, no phase changer; no controller
      # names it (the controllers sit on branch 5, T1, and 7, T3)
      t2 = net.branchVec[6]
      @test t2.has_ratio_tap && t2.ratio_tap_declared && !t2.has_phase_tap
      @test t2.tap_step == 0.0125 && t2.tap_ratio == 1.0
      @test Sparlectra._tap_position(t2.tap_ratio, 1.0, t2.tap_step) == 0
      @test t2.tap_min ≈ 1 / (1 + 12 * 0.0125) && t2.tap_max ≈ 1 / (1 - 12 * 0.0125)
      controllers = Sparlectra.collect_outer_controllers(net)
      @test sort([(c.trafo, c.mode) for c in controllers]) == [("5", :branch_active_power), ("5", :voltage), ("7", :voltage)]
      cres = Sparlectra.run_control!(net; verbose = 0)
      @test cres.converged
      row(n) = only(r for r in cres.controllers if r.controller_name == n)
      @test abs(row("T1_Nordwerk_Kreuzau OLTC").achieved_vm_pu - 1.0) <= 0.01
      @test abs(row("T3_Mittelbach_Lindenau OLTC").achieved_vm_pu - 1.0) <= 0.01
      @test abs(row("T1_Nordwerk_Kreuzau PST").achieved_p_mw - 60.0) <= 4.0
      # the result rows read the solved net: the bus voltages and the T1
      # flow agree with them
      vm(b) = only(n for n in net.nodeVec if n.comp.cName == b)._vm_pu
      @test row("T1_Nordwerk_Kreuzau OLTC").achieved_vm_pu ≈ vm("Kreuzau_110")
      @test row("T3_Mittelbach_Lindenau OLTC").achieved_vm_pu ≈ vm("Lindenau_20")
      @test row("T1_Nordwerk_Kreuzau PST").achieved_p_mw ≈ Sparlectra.get_branch_p_from_to_mw(net, "Nordwerk_220", "Kreuzau_110")
    end)() end

    # the shipped cases load through the Web
    # UI. The chooser offers them without any cache copy, and the run path
    # stages a bundled case into the cache WITH its sidecars (the per-case
    # config carries the machine-neutrality pin, so losing it on the copy
    # would silently undo that guarantee).
    @testset "shipped cases load through the Web UI" begin (function ()
      app_root = normpath(joinpath(dirname(@__DIR__)))
      offered = SparlectraApp._webui_bundled_scf_options(app_root)
      @test "sp_case14.scf.json" in offered
      @test "sp_case6tap.scf.json" in offered
      @test "sp_casePST.scf.json" in offered
      # 0.12.5: the feeder/Q-limit example cases and the plain PGM files are
      # offered too, each with its case configuration sidecar
      @test "feeder3_hardpv_A.scf.json" in offered
      @test "feeder3_hardpv_pgm.json" in offered
      @test !("README.md" in offered)
      cache = mktempdir()
      ctx = SparlectraApp._webui_case_context(; application_root = app_root, case_directory = cache)
      @test "sp_case14.scf.json" in ctx.casefiles
      staged = SparlectraApp._webui_stage_bundled_case!(app_root, cache, "sp_case14.scf.json")
      @test staged == joinpath(cache, "sp_case14.scf.json")
      @test isfile(staged)
      @test isfile(joinpath(cache, "sp_case14.config.yaml"))
      staged_pgm = SparlectraApp._webui_stage_bundled_case!(app_root, cache, "feeder3_hardpv_pgm.json")
      @test staged_pgm == joinpath(cache, "feeder3_hardpv_pgm.json")
      @test isfile(joinpath(cache, "feeder3_hardpv_pgm.config.yaml"))
      @test isfile(joinpath(cache, "sp_case14.measurements.csv"))
      # BOTH measurement variants travel: without the baddata twin the
      # diagnostics demo would die with a file error instead of a message
      @test isfile(joinpath(cache, "sp_case14.measurements.baddata.csv"))
      # a user-modified sidecar in the cache survives a second staging
      write(joinpath(cache, "sp_case14.config.yaml"), "config_version: 1\nscope: case\ncase: sp_case14.scf.json\npower_flow:\n  tol: 1.0e-8\n  max_iter: 44\n")
      SparlectraApp._webui_stage_bundled_case!(app_root, cache, "sp_case14.scf.json")
      @test occursin("max_iter: 44", read(joinpath(cache, "sp_case14.config.yaml"), String))
      # the staged case runs through the service front door, and the run
      # PROVES the pinning arrived (succeeded only confirms the copy, not
      # the config): the effective-config artifact must name the
      # case configuration file as a resolution source
      run = SparlectraApp.start_powerflow_run(Dict("casefile" => "sp_case14.scf.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(cache, "runs")); case_directory = cache)
      @test run["status"] == "succeeded"
      eff_path = joinpath(String(run["output_dir"]), "effective_config.yaml")
      @test isfile(eff_path)
      eff = read(eff_path, String)
      pf_seg = eff[first(findfirst("  power_flow:", eff)):end]
      tol_seg = pf_seg[first(findfirst("    tol:", pf_seg)):min(end, first(findfirst("    tol:", pf_seg)) + 250)]
      # the pinned tolerance is REPORTED as case_sidecar (this assert found
      # the source-report labeling bug: the case level used to surface
      # under the caller's override label)
      @test occursin("source: case_sidecar", tol_seg)
      @test occursin("value: 1.0e-8", tol_seg)
      # an unknown name stages nothing
      @test SparlectraApp._webui_stage_bundled_case!(app_root, cache, "sp_nope.scf.json") === nothing
    end)() end

    # REPLACES the former testset "hostile general config cannot move the
    # shipped fixtures" . That test could not fail:
    # it handed a hostile configuration to `import_case` and then called
    # `runpf!(net; verbose = 0)`, and that call form solved under the GLOBAL
    # configuration, so the hostile values never reached the solver. Its
    # premise was wrong twice over: `import_case` does not merge the case
    # sidecar at all (that happens in `resolve_config`, on the service path),
    # and a sidecar that pins only `power_flow.tol` cannot neutralize a
    # hostile `distributed_slack` anyway. Measured: with the
    # hostile configuration actually applied, sp_case14 throws
    # "distributed slack: no valid participant ... p_mode=pmax_weighted".
    #
    # What holds instead, and what this testset guards: the configuration a
    # net was imported with is the one it is solved under.
    @testset "the imported configuration reaches the solver" begin (function ()
      scf_path = joinpath(_DEMO_CASE_DIR, "sp_case14.scf.json")

      # structural: the net carries what it was built with
      strict = Sparlectra.SparlectraConfig(Dict{String,Any}("power_flow" => Dict{String,Any}("max_iter" => 1)))
      carried = Sparlectra.import_case(scf_path, strict).net._import_config
      @test carried isa Sparlectra.SparlectraConfig
      @test carried.powerflow.max_iter == 1

      # behavioural, and the actual regression: a cold start cannot converge
      # in a single iteration. Under the old behaviour this run silently used
      # the global max_iter and converged, which is exactly what hid the bug.
      cold_net = Sparlectra.import_case(scf_path, strict).net
      for nd in cold_net.nodeVec
        Sparlectra.setVmVa!(node = nd, vm_pu = 1.0, va_deg = 0.0)
      end
      _, erg_strict = runpf!(cold_net; verbose = 0)
      @test erg_strict != 0

      # the same cold start converges under the package defaults, so the
      # failure above comes from the configuration and not from the case
      ref_net = Sparlectra.import_case(scf_path, _DEMO_REF_CONFIG()).net
      for nd in ref_net.nodeVec
        Sparlectra.setVmVa!(node = nd, vm_pu = 1.0, va_deg = 0.0)
      end
      ite_ref, erg_ref = runpf!(ref_net; verbose = 0)
      @test erg_ref == 0
      @test ite_ref > 1

      # an explicit `config` still wins over the carried one
      loose = Sparlectra.SparlectraConfig(Dict{String,Any}("power_flow" => Dict{String,Any}("max_iter" => 50)))
      win_net = Sparlectra.import_case(scf_path, strict).net
      for nd in win_net.nodeVec
        Sparlectra.setVmVa!(node = nd, vm_pu = 1.0, va_deg = 0.0)
      end
      _, erg_win = runpf!(win_net, loose; verbose = 0)
      @test erg_win == 0
    end)() end
  end)() end
end

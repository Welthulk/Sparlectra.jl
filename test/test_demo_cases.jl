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
# purpose: the shipped demo cases (data/scf/sp_case*)
#          are regression fixtures, not decoration: each case is read through
#          import_case like any other SCF file, all five run kinds are
#          executed (power flow, N-1, explicit scenarios, state estimation,
#          short circuit), and the results are compared against the tracked
#          reference fixture tracked next to this file. The fixtures are
#          part of the repository: they are not generated during a test run,
#          they are the recorded reference the run is compared against.

const _DEMO_CASE_DIR = normpath(joinpath(dirname(@__DIR__), "data", "scf"))
const _DEMO_CASE_FIXTURES = normpath(joinpath(@__DIR__, "fixtures", "demo_cases"))
const _DEMO_CASE_NAMES = ("sp_case5", "sp_case14", "sp_case60", "sp_case188")
# same rule for the provenance stamp as the SCF tests use
include("test_scf_support.jl")
# the fixtures are machine-neutral: they were taken against the PACKAGE
# defaults, so the comparison must never
# resolve through this machine's configuration.yaml
_DEMO_REF_CONFIG() = Sparlectra.SparlectraConfig(Dict())

function _demo_fixture(name::String)::Dict{String,Any}
  path = joinpath(_DEMO_CASE_FIXTURES, "$(name).json")
  # the fixture is TRACKED, so a missing one means the checkout is
  # incomplete, not that something has to be generated first
  isfile(path) || error("demo case fixture missing: $(path); the reference fixtures are tracked in test/fixtures/demo_cases, so this checkout is incomplete")
  return Sparlectra.scf_json_parse(read(path, String))
end

function run_demo_case_tests()
  @testset "Shipped demo cases (regression fixtures)" begin (function ()
    for name in _DEMO_CASE_NAMES
      @testset "$(name)" begin (function ()
        scf_path = joinpath(_DEMO_CASE_DIR, "$(name).scf.json")
        fixture = _demo_fixture(name)
        @test isfile(scf_path)
        # provenance: the self-built statement travels in the file itself
        @test occursin("Self-built Sparlectra demo network", read(scf_path, String))
        # the per-case config file exists and pins the solve tolerance: its
        # presence pins the resolution chain to package defaults, so the
        # 1e-6 pu fixtures never resolve against a machine configuration
        cfg_path = Sparlectra.case_config_path(scf_path)
        @test isfile(cfg_path)
        @test occursin("scope: case", read(cfg_path, String))
        @test Sparlectra.load_case_config(scf_path)["power_flow.tol"] == 1.0e-8

        # 1) power flow through the standard import path (start_state applies)
        res = Sparlectra.import_case(scf_path, _DEMO_REF_CONFIG())
        net = res.net
        ite, erg = runpf!(net; verbose = 0)
        @test erg == 0
        @test ite == fixture["power_flow"]["iterations"]
        calcNetLosses!(net)
        p_losses, _ = Sparlectra.getTotalLosses(net = net)
        @test isapprox(p_losses, fixture["power_flow"]["losses_mw"]; atol = 1e-5)
        # the reported flows belong to the solved state: at both terminals of
        # every closed branch they equal the flow of the model the Ybus
        # stamps (live tap ratio and shift). sp_case14 has a nameplate tap off
        # its neutral position, where calcNetLosses! read the neutral ratio
        V = ComplexF64[something(nd._vm_pu, 1.0) * cis(deg2rad(something(nd._va_deg, 0.0))) for nd in net.nodeVec]
        flow_gap = 0.0
        for br in net.branchVec
          (Sparlectra._branch_terminal_state(br) === :closed && !(br.fromBus in net.isoNodes) && !(br.toBus in net.isoNodes)) || continue
          sf = Sparlectra._sv_branch_terminal_flow(br, V, :from, net.baseMVA)
          st = Sparlectra._sv_branch_terminal_flow(br, V, :to, net.baseMVA)
          flow_gap = max(flow_gap, abs(complex(br.fBranchFlow.pFlow, br.fBranchFlow.qFlow) - sf), abs(complex(br.tBranchFlow.pFlow, br.tBranchFlow.qFlow) - st))
        end
        @test flow_gap <= 1e-9
        ids = Sparlectra.scf_id_map(net)
        bus_state = fixture["power_flow"]["bus_state"]
        @test length(bus_state) == length(net.nodeVec)
        for i in eachindex(net.nodeVec)
          ref = bus_state[string(ids.node[i])]
          @test isapprox(net.nodeVec[i]._vm_pu, ref["vm_pu"]; atol = 1e-6)
          @test isapprox(net.nodeVec[i]._va_deg, ref["va_deg"]; atol = 1e-5)
        end

        # 2) state estimation on the measurements the file itself carries
        # the topology precheck advises on two of the shipped cases by design;
        # captured, anything else that warns fails
        se = run_with_expected_warnings(() -> runse!(net), (r"topology precheck reported",))
        @test se.converged
        @test se.dof == fixture["state_estimation"]["dof"]
        @test isapprox(se.objectiveJ, fixture["state_estimation"]["objective"]; atol = 1e-4)

        # 3) short circuit at the fixture's buses: the feeder with the data
        # of the file, the machines (the file names no data for them) with
        # the engine's default, which is why every row is flagged and the
        # engine warns once per machine
        sc_buses = sort(collect(keys(fixture["short_circuit"])))
        sc = run_with_expected_warnings(() -> runShortCircuit!(net, net.sc_sources; buses = sc_buses, case = :max), ("has no usable x''_d", "has no usable ratedS", "has neither x''_d nor ratedS"))
        @test length(sc.rows) == length(sc_buses)
        @test length(net.sc_sources.external_network_injections) == 1
        @test all(m -> m.satDirectSubtransX_pu === nothing, net.sc_sources.synchronous_machines)
        @test all(row -> row.contains_defaulted_data, sc.rows)
        for row in sc.rows
          ref = fixture["short_circuit"][String(row.bus)]
          @test isapprox(row.ik_kA, ref["ik_kA"]; atol = 1e-5)
          @test isapprox(row.sk_MVA, ref["sk_MVA"]; atol = 1e-3)
        end

        # 4) N-1 over the generated branch list on a fresh import
        n1net = Sparlectra.import_case(scf_path, _DEMO_REF_CONFIG()).net
        n1 = runContingencies!(n1net, generateN1Branches(n1net))
        @test length(n1) == fixture["n1"]["total"]
        @test count(r -> r.converged, n1) == fixture["n1"]["converged"]
        @test count(r -> r.converged && r.max_branch_loading_pct !== nothing && r.max_branch_loading_pct > 100.0, n1) == fixture["n1"]["overloads"]

        # 5) the file's own explicit scenarios through runScenarios!
        icase = read_scf_json(scf_path)
        iset = Sparlectra.scf_case_scenarios(icase)
        @test iset !== nothing
        @test String(iset.mode) == fixture["scenarios"]["mode"]
        explicit_ref = fixture["scenarios"]["explicit"]
        @test sort(collect(keys(explicit_ref))) == sort([s.name for s in iset.scenarios])
        scen_net = Sparlectra.import_case(scf_path, _DEMO_REF_CONFIG()).net
        scen_res = runScenarios!(scen_net, Sparlectra.ScenarioSet(scenarios = iset.scenarios, mode = :explicit); index = Sparlectra.ScenarioIndex(icase))
        for r in scen_res
          ref = explicit_ref[r.name]
          @test r.converged == ref["converged"]
          # islanding scenarios carry no finite loading; the fixture then
          # omits the key
          ref_loading = get(ref, "max_branch_loading_pct", nothing)
          ref_loading === nothing || @test isapprox(r.max_branch_loading_pct, ref_loading; atol = 1e-5)
        end

        # 6) cases with a declared controller: the parameters in the file
        # must match the fixture (a changed controller default then shows
        # up as a parameter diff next to the symptom), and the control run
        # reproduces its reference behavior
        if haskey(fixture, "control")
          ctrl_ref = fixture["control"]
          # a shipped case's controllers MUST settle; equality with the
          # fixture alone once froze a non-converging SSSC as the
          # reference (the first live run found it)
          @test ctrl_ref["converged"] === true
          file_params = icase.sparlectra.components.controllers
          @test Dict{String,Any}(String(k) => v for (k, v) in file_params) == ctrl_ref["parameters"]
          ctrl_net = Sparlectra.import_case(scf_path, _DEMO_REF_CONFIG()).net
          cres = Sparlectra.run_control!(ctrl_net; verbose = 0)
          @test cres.converged == ctrl_ref["converged"]
          @test cres.outer_iterations == ctrl_ref["outer_iterations"]
          @test cres.powerflow_solves == ctrl_ref["powerflow_solves"]
        end

        # size budget: no single bundle runs away; the BINDING cap is the
        # task's total (all cases together well under a megabyte, checked
        # after the loop). 420k since the place names: the
        # benchmark bundle measured 412k, up 3.2 percent from 399k, because
        # measurements and extra carry reference NAMES (Birkloh_110, not
        # C9); the names are already compact (max 9 chars a place). The
        # planned measurement id-mapping option (post-release) is the lever
        # that shrinks this below the old level
        bundle = filter(f -> startswith(f, name), readdir(_DEMO_CASE_DIR))
        @test sum(filesize(joinpath(_DEMO_CASE_DIR, f)) for f in bundle) < 420_000
        # the bad-data variant exists and differs in exactly the marked way
        @test isfile(joinpath(_DEMO_CASE_DIR, "$(name).measurements.baddata.csv"))
        @test occursin("baddata", read(joinpath(_DEMO_CASE_DIR, "$(name).measurements.baddata.csv"), String))
      end)() end
    end
    all_files = filter(f -> startswith(f, "sp_case"), readdir(_DEMO_CASE_DIR))
    @test sum(filesize(joinpath(_DEMO_CASE_DIR, f)) for f in all_files) < 900_000

    # the shipped cases load through the Web
    # UI. The chooser offers them without any cache copy, and the run path
    # stages a bundled case into the cache WITH its sidecars (the per-case
    # config carries the machine-neutrality pin, so losing it on the copy
    # would silently undo that guarantee).
    @testset "shipped cases load through the Web UI" begin (function ()
      app_root = normpath(joinpath(dirname(@__DIR__)))
      offered = SparlectraApp._webui_bundled_scf_options(app_root)
      @test "sp_case14.scf.json" in offered
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

      # and it still lands on the shipped fixture
      fixture = _demo_fixture("sp_case14")
      bus_state = fixture["power_flow"]["bus_state"]
      ids = Sparlectra.scf_id_map(ref_net)
      for i in eachindex(ref_net.nodeVec)
        @test isapprox(ref_net.nodeVec[i]._vm_pu, bus_state[string(ids.node[i])]["vm_pu"]; atol = 1e-6)
      end

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

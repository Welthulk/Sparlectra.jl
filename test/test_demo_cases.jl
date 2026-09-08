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
# purpose: the shipped demo cases (data/scf/sp_case*, task_demo_cases_v0100)
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
# defaults (colleague review 2026-09-03), so the comparison must never
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
  @testset "Shipped demo cases (regression fixtures)" begin
    for name in _DEMO_CASE_NAMES
      @testset "$(name)" begin
        scf_path = joinpath(_DEMO_CASE_DIR, "$(name).scf.json")
        fixture = _demo_fixture(name)
        @test isfile(scf_path)
        # provenance: the self-built statement travels in the file itself
        @test occursin("Self-built Sparlectra demo network", read(scf_path, String))
        # the per-case config file exists and pins the solve tolerance: its
        # presence pins the resolution chain to package defaults, so the
        # 1e-6 pu fixtures never resolve against a machine configuration
        # (colleague review 2026-09-03)
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
        ids = Sparlectra.scf_id_map(net)
        bus_state = fixture["power_flow"]["bus_state"]
        @test length(bus_state) == length(net.nodeVec)
        for i in eachindex(net.nodeVec)
          ref = bus_state[string(ids.node[i])]
          @test isapprox(net.nodeVec[i]._vm_pu, ref["vm_pu"]; atol = 1e-6)
          @test isapprox(net.nodeVec[i]._va_deg, ref["va_deg"]; atol = 1e-5)
        end

        # 2) state estimation on the measurements the file itself carries
        se = runse!(net, Vector{Sparlectra.Measurement}(net.measurements), Sparlectra.StateEstimationConfig())
        @test se.converged
        @test se.dof == fixture["state_estimation"]["dof"]
        @test isapprox(se.objectiveJ, fixture["state_estimation"]["objective"]; atol = 1e-4)

        # 3) short circuit at the fixture's buses (feeder data from the file)
        sc_buses = sort(collect(keys(fixture["short_circuit"])))
        sc = runShortCircuit!(net, net.sc_sources; buses = sc_buses, case = :max)
        @test length(sc.rows) == length(sc_buses)
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
          # reference (maintainer's first live run found it)
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
        # after the loop). 420k since the place names (2026-09-04): the
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
      end
    end
    all_files = filter(f -> startswith(f, "sp_case"), readdir(_DEMO_CASE_DIR))
    @test sum(filesize(joinpath(_DEMO_CASE_DIR, f)) for f in all_files) < 900_000

    # Why there is no byte round trip for these four, in case someone reaches
    # for one here: importing a shipped case and exporting it again does NOT
    # reproduce the file. `scenarios`, the `short_circuit` block with its PGM
    # `fault` row, the `extra` names and the measurement provenance are EXPORT
    # ARGUMENTS, not network state (measured 2026-09-08: sp_case5 comes back as
    # 15354 characters against the shipped 16643). Handing those blocks back to
    # the exporter from the file under test would compare it against its own
    # input. These files are built by a generator that passes those arguments
    # and is their provenance record, not by a re-export; the property such a
    # comparison would be reaching for is checked directly below instead.
    #
    # Colleague review 2026-09-08, after the 0.10.0 -> 0.11.0 bump turned a
    # fixture comparison red: a test that goes red on every version bump
    # trains the reflex to regenerate the file, and the next time the diff
    # may be more than the stamp line, with the regeneration hiding it. The
    # shipped cases must therefore not depend on the release recorded in
    # them at all. A copy stamped with a version that does not exist has to
    # import to the same network and solve to the same result; nothing may
    # read `created_by` except a human looking for provenance.
    @testset "the provenance stamp does not reach behavior" begin
      d = mktempdir()
      for name in _DEMO_CASE_NAMES
        original = joinpath(_DEMO_CASE_DIR, "$(name).scf.json")
        # the case config next to the file pins the tolerance; it has to
        # travel with the copy or the comparison resolves differently
        cp(Sparlectra.case_config_path(original), joinpath(d, basename(Sparlectra.case_config_path(original))); force = true)
        foreign = scf_with_foreign_stamp(original, d)
        @test occursin("Sparlectra 0.0.0-nonexistent", read(foreign, String))
        # and it is the ONLY difference
        @test scf_without_stamp(read(foreign, String)) == scf_without_stamp(read(original, String))

        ref = Sparlectra.import_case(original, _DEMO_REF_CONFIG()).net
        alt = Sparlectra.import_case(foreign, _DEMO_REF_CONFIG()).net
        @test length(alt.nodeVec) == length(ref.nodeVec)
        @test length(alt.branchVec) == length(ref.branchVec)
        @test length(alt.measurements) == length(ref.measurements)
        ite_ref, erg_ref = runpf!(ref; verbose = 0)
        ite_alt, erg_alt = runpf!(alt; verbose = 0)
        @test (ite_alt, erg_alt) == (ite_ref, erg_ref)
        calcNetLosses!(ref)
        calcNetLosses!(alt)
        p_ref, q_ref = Sparlectra.getTotalLosses(net = ref)
        p_alt, q_alt = Sparlectra.getTotalLosses(net = alt)
        @test p_alt == p_ref && q_alt == q_ref
      end
    end

    # maintainer request 2026-09-03: the shipped cases load through the Web
    # UI. The chooser offers them without any cache copy, and the run path
    # stages a bundled case into the cache WITH its sidecars (the per-case
    # config carries the machine-neutrality pin, so losing it on the copy
    # would silently undo that guarantee).
    @testset "shipped cases load through the Web UI" begin
      app_root = normpath(joinpath(dirname(@__DIR__)))
      offered = Sparlectra._webui_bundled_scf_options(app_root)
      @test "sp_case14.scf.json" in offered
      @test "sp_casePST.scf.json" in offered
      cache = mktempdir()
      ctx = Sparlectra._webui_case_context(; application_root = app_root, case_directory = cache)
      @test "sp_case14.scf.json" in ctx.casefiles
      staged = Sparlectra._webui_stage_bundled_case!(app_root, cache, "sp_case14.scf.json")
      @test staged == joinpath(cache, "sp_case14.scf.json")
      @test isfile(staged)
      @test isfile(joinpath(cache, "sp_case14.config.yaml"))
      @test isfile(joinpath(cache, "sp_case14.measurements.csv"))
      # BOTH measurement variants travel: without the baddata twin the
      # diagnostics demo would die with a file error instead of a message
      @test isfile(joinpath(cache, "sp_case14.measurements.baddata.csv"))
      # a user-modified sidecar in the cache survives a second staging
      write(joinpath(cache, "sp_case14.config.yaml"), "config_version: 1\nscope: case\ncase: sp_case14.scf.json\npower_flow:\n  tol: 1.0e-8\n  max_iter: 44\n")
      Sparlectra._webui_stage_bundled_case!(app_root, cache, "sp_case14.scf.json")
      @test occursin("max_iter: 44", read(joinpath(cache, "sp_case14.config.yaml"), String))
      # the staged case runs through the service front door, and the run
      # PROVES the pinning arrived (colleague: succeeded only confirms the
      # copy, not the config): the effective-config artifact must name the
      # case configuration file as a resolution source
      run = Sparlectra.start_powerflow_run(Dict("casefile" => "sp_case14.scf.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(cache, "runs")); case_directory = cache)
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
      @test Sparlectra._webui_stage_bundled_case!(app_root, cache, "sp_nope.scf.json") === nothing
    end

    # REPLACES the former testset "hostile general config cannot move the
    # shipped fixtures" (2026-09-03 to 2026-09-08). That test could not fail:
    # it handed a hostile configuration to `import_case` and then called
    # `runpf!(net; verbose = 0)`, and that call form solved under the GLOBAL
    # configuration, so the hostile values never reached the solver. Its
    # premise was wrong twice over: `import_case` does not merge the case
    # sidecar at all (that happens in `resolve_config`, on the service path),
    # and a sidecar that pins only `power_flow.tol` cannot neutralize a
    # hostile `distributed_slack` anyway. Measured on 2026-09-08: with the
    # hostile configuration actually applied, sp_case14 throws
    # "distributed slack: no valid participant ... p_mode=pmax_weighted".
    #
    # What holds instead, and what this testset guards: the configuration a
    # net was imported with is the one it is solved under.
    @testset "the imported configuration reaches the solver" begin
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
    end
  end
end

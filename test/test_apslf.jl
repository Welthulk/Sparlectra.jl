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

# file: test/test_apslf.jl
#
# Tests the APSLF (AnalyticLoadFlow.jl) integration: config validation,
# controller rejection, WebUI form parsing, and the "not installed" error
# path always run. The adapter-mapping, standalone-solve, and start-value
# generator tests additionally require AnalyticLoadFlow.jl to be resolvable
# in the active Julia environment (it is only a weak dependency of
# Sparlectra, loaded via ext/SparlectraAnalyticLoadFlowExt.jl) — they are
# skipped with an informational message otherwise, matching the weak-
# dependency design: production code never requires AnalyticLoadFlow.jl to
# be installed. To exercise the full test surface locally, `Pkg.add`
# AnalyticLoadFlow.jl into a Julia environment stacking this checkout and
# `using AnalyticLoadFlow` before running the extended test profile.

using Sparlectra
using Test
using Random

const _APSLF_AVAILABLE = Base.find_package("AnalyticLoadFlow") !== nothing
_APSLF_AVAILABLE && @eval import AnalyticLoadFlow

function run_apslf_tests()
  @testset "APSLF (AnalyticLoadFlow.jl) integration" begin
    @testset "Config validation: power_flow.solver / apslf / apslf_start" begin
      default_cfg = Sparlectra.SparlectraConfig()
      @test default_cfg.powerflow.solver === :rectangular
      @test default_cfg.powerflow.apslf.order == 40
      @test default_cfg.powerflow.apslf.use_pade === true
      @test default_cfg.powerflow.apslf.nr_polish === true
      @test default_cfg.powerflow.apslf_start.enabled === false
      @test default_cfg.powerflow.apslf_start.order == 40

      @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("solver" => "polar")))
      cfg_apslf = Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("solver" => "apslf")))
      @test cfg_apslf.powerflow.solver === :apslf

      @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("apslf" => Dict("order" => 0))))
      @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("apslf_start" => Dict("order" => -1))))

      # solver=apslf together with apslf_start.enabled=true is rejected (the
      # start-value generator only makes sense ahead of the NR solve).
      @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("solver" => "apslf", "apslf_start" => Dict("enabled" => true))))

      # apslf_start.enabled=true alone (solver stays rectangular) is fine.
      cfg_hybrid = Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("apslf_start" => Dict("enabled" => true, "order" => 12))))
      @test cfg_hybrid.powerflow.apslf_start.enabled
      @test cfg_hybrid.powerflow.apslf_start.order == 12
      @test cfg_hybrid.powerflow.solver === :rectangular

      # Unknown keys under the new sections are rejected like any other section.
      bad_key_file = tempname() * ".yaml"
      write(bad_key_file, "power_flow:\n  apslf:\n    bogus_key: 1\n")
      @test_throws ArgumentError Sparlectra.load_sparlectra_config(bad_key_file; reload = true)
    end

    @testset "Controller + APSLF solver rejection" begin
      mpc = Sparlectra.MatpowerIO.read_case(ensure_casefile("case14.m"))
      net = Sparlectra.createNetFromMatPowerCase(mpc = mpc, flatstart = true)
      addPowerTransformerControl!(net; trafo = "B_2WT_1_4_7", mode = :voltage, target_bus = "7", target_vm_pu = 1.0, control_ratio = true)
      @test length(Sparlectra.collect_outer_controllers(net)) == 1

      cfg = Sparlectra.SparlectraConfig(powerflow = Sparlectra.PowerFlowConfig(solver = :apslf), output = OutputConfig(logfile_results = :off))
      # No AnalyticLoadFlow.jl is required here: the controller check runs
      # before the solver is even constructed, so this rejects identically
      # regardless of whether the extension is loaded.
      @test_throws ArgumentError run_sparlectra(net = net, config = cfg)
    end

    @testset "WebUI form parsing for APSLF fields -> effective config" begin
      form = Dict{String,Any}(
        "casefile" => "case14.m",
        "power_flow_solver" => "apslf",
        "power_flow_apslf_order" => "25",
        "power_flow_apslf_use_pade" => "false",
        "power_flow_apslf_nr_polish" => "true",
        "power_flow_apslf_start_enabled" => "false",
        "power_flow_apslf_start_order" => "40",
      )
      request = Sparlectra.powerflow_webui_request(form)
      overrides = request["config_overrides"]
      @test overrides["power_flow.solver"] == "apslf"
      @test overrides["power_flow.apslf.order"] === 25
      @test overrides["power_flow.apslf.use_pade"] === false
      @test overrides["power_flow.apslf.nr_polish"] === true
      @test overrides["power_flow.apslf_start.enabled"] === false
      @test overrides["power_flow.apslf_start.order"] === 40

      nested = Sparlectra.validate_gui_config_overrides(overrides)
      cfg, _ = Sparlectra._load_api_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, nested)
      @test cfg.powerflow.solver === :apslf
      @test cfg.powerflow.apslf.order == 25
      @test cfg.powerflow.apslf.use_pade === false
      @test cfg.powerflow.apslf.nr_polish === true

      # Unchecked checkbox explicitly submits "false" (hidden+checkbox HTML
      # pairing), the same semantics as the Q-limit checkbox pattern: the key
      # is present with an explicit false, not omitted.
      disabled_form = copy(form)
      disabled_form["power_flow_apslf_use_pade"] = "false"
      @test Sparlectra.powerflow_webui_request(disabled_form)["config_overrides"]["power_flow.apslf.use_pade"] === false

      # A field genuinely absent from the submitted form (e.g. a stale
      # client) is skipped rather than defaulted to false.
      missing_field_form = copy(form)
      delete!(missing_field_form, "power_flow_apslf_start_enabled")
      @test !haskey(Sparlectra.powerflow_webui_request(missing_field_form)["config_overrides"], "power_flow.apslf_start.enabled")

      # GUI-layer validation rejects the same invalid values as the core config.
      @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.solver" => "polar"))
      @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.apslf.order" => 0))

      # Downstream config assembly still rejects the solver+start conflict
      # even when both values arrive through the GUI override path.
      conflict_overrides = Dict{String,Any}("power_flow.solver" => "apslf", "power_flow.apslf_start.enabled" => true)
      conflict_nested = Sparlectra.validate_gui_config_overrides(conflict_overrides)
      @test_throws ArgumentError Sparlectra._load_api_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, conflict_nested)
    end

    @testset "AnalyticLoadFlow.jl \"not installed\" error path" begin
      if _APSLF_AVAILABLE
        @info "AnalyticLoadFlow.jl is resolvable in this test session; skipping the not-installed error-path check (would require an actually uninstalled environment)."
      else
        err = try
          Sparlectra.apslf_solver()
          nothing
        catch caught
          caught
        end
        @test err isa ErrorException
        message = sprint(showerror, err)
        @test occursin("AnalyticLoadFlow.jl", message)
        @test occursin("nicht installiert", message)
      end
    end

    if _APSLF_AVAILABLE
      @testset "Adapter mapping (PFModel -> AnalyticLoadFlow spec, PF ordering)" begin
        ext = Base.get_extension(Sparlectra, :SparlectraAnalyticLoadFlowExt)
        @test ext !== nothing

        net = createTest3BusNet()
        model = buildPfModel(net; flatstart = true, include_limits = false)
        n = length(model.busIdx_net)

        spec = ext._apslf_spec_from_model(model)
        @test spec.Y === model.Ybus
        @test spec.bustype == model.busType
        @test spec.Pspec ≈ real.(model.Sspec)
        @test spec.Qspec ≈ imag.(model.Sspec)
        @test spec.Vm == model.Vset
        @test spec.Qmin == fill(-Inf, n)
        @test spec.Qmax == fill(Inf, n)
        @test spec.slack == model.slack_idx

        model_q = buildPfModel(net; flatstart = true, include_limits = true)
        spec_q = ext._apslf_spec_from_model(model_q)
        @test spec_q.Qmin == model_q.qmin_pu
        @test spec_q.Qmax == model_q.qmax_pu

        solver = ext.ApslfSolver(order = 20, use_pade = true, nr_polish = true)
        sol = solvePf(solver, model)
        @test sol isa PFSolution
        @test length(sol.V) == n
        @test sol.meta.solver === :apslf
        @test sol.meta.mode isa Symbol
      end

      @testset "Standalone APSLF run on case14 (convergence + NR agreement)" begin
        mpc = Sparlectra.MatpowerIO.read_case(ensure_casefile("case14.m"))

        net_nr = Sparlectra.createNetFromMatPowerCase(mpc = mpc, flatstart = true)
        iters_nr, status_nr = runpf!(net_nr, 30, 1e-8, 0)
        @test status_nr == 0

        net_apslf = Sparlectra.createNetFromMatPowerCase(mpc = mpc, flatstart = true)
        solver = apslf_solver(order = 40, use_pade = true, nr_polish = true)
        iters_apslf, status_apslf, sol = runpf_external!(net_apslf, solver; tol = 1e-8)
        @test status_apslf == 0
        @test sol.converged

        max_dVm = maximum(abs.(getfield.(net_apslf.nodeVec, :_vm_pu) .- getfield.(net_nr.nodeVec, :_vm_pu)))
        max_dVa = maximum(abs.(getfield.(net_apslf.nodeVec, :_va_deg) .- getfield.(net_nr.nodeVec, :_va_deg)))
        @test max_dVm < 1e-6
        @test max_dVa < 1e-4

        # Task 3 wiring: routing through the framework run path yields the
        # same solver identity and result contract as the direct external-
        # solver call above.
        cfg = Sparlectra.SparlectraConfig(powerflow = Sparlectra.PowerFlowConfig(solver = :apslf), output = OutputConfig(logfile_results = :off))
        net_framework = Sparlectra.createNetFromMatPowerCase(mpc = mpc, flatstart = true)
        result = run_sparlectra(net = net_framework, config = cfg)
        @test result.final_converged
        @test result.diagnostics.solver === :apslf
      end

      @testset "APSLF start-value generator guard and nr_polish=false in start mode" begin
        mpc = Sparlectra.MatpowerIO.read_case(ensure_casefile("case14.m"))
        net = Sparlectra.createNetFromMatPowerCase(mpc = mpc, flatstart = true)
        model = buildPfModel(net; flatstart = true, include_limits = false)

        # Disabled: no-op, restores the exact incoming vector.
        Vraw = copy(model.V0)
        V_disabled, summary_disabled = Sparlectra._run_guarded_apslf_start(model.Ybus, Vraw, model.Sspec, model.busType, model.Vset, model.slack_idx; enabled = false, order = 40, baseMVA = model.baseMVA)
        @test !summary_disabled.apslf_start_enabled
        @test !summary_disabled.apslf_start_attempted
        @test V_disabled == Vraw

        # case14 flat start: the raw (nr_polish=false) APSLF series result is
        # a worse start than the flat profile for this case (verified
        # separately: order 40 with nr_polish=true converges to ~1e-15,
        # nr_polish=false plateaus around 0.25 mismatch here) — the guard
        # must reject it and restore the original start. This also serves as
        # the nr_polish=false contract check: if the generator used
        # nr_polish=true internally, the candidate would trivially win.
        V_flat, summary_flat = Sparlectra._run_guarded_apslf_start(model.Ybus, Vraw, model.Sspec, model.busType, model.Vset, model.slack_idx; enabled = true, order = 40, baseMVA = model.baseMVA)
        @test summary_flat.apslf_start_enabled
        @test summary_flat.apslf_start_attempted
        @test summary_flat.apslf_start_accepted === false
        @test summary_flat.apslf_start_reason === :not_improved
        @test V_flat == Vraw

        # A deliberately perturbed, worse-than-APSLF start must be accepted.
        rng = Random.MersenneTwister(42)
        bad = model.V0 .* cis.(0.6 .* (rand(rng, length(model.V0)) .- 0.5))
        bad[model.slack_idx] = model.V0[model.slack_idx]
        V_bad, summary_bad = Sparlectra._run_guarded_apslf_start(model.Ybus, bad, model.Sspec, model.busType, model.Vset, model.slack_idx; enabled = true, order = 40, baseMVA = model.baseMVA)
        @test summary_bad.apslf_start_accepted === true
        @test summary_bad.apslf_start_reason === :improved
        @test V_bad != bad
        @test all(v -> isfinite(real(v)) && isfinite(imag(v)), V_bad)
      end
    else
      @info "AnalyticLoadFlow.jl is not resolvable in this Julia environment; skipping APSLF adapter-mapping, standalone-run, and start-generator tests. Add it to a stacked environment and rerun the extended profile to exercise them."
    end
  end
  return true
end

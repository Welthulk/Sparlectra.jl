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

# file: test/test_legacy_auto_powerflow.jl
# purpose: legacy profile: the testsets of test/test_auto_powerflow.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_auto_powerflow.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_auto_powerflow_tests()
  # cut from testset: auto powerflow
    @testset "escalation ladder keeps the step-control exclusions (issue #5)" begin (function ()
      # A hard_case profile starts with autodamp AND merit, a consistent pair.
      # A ladder stage used to switch autodamp off while leaving merit on, and
      # runpf_rectangular! rightly refused with
      #   ArgumentError: merit_enabled=true requires autodamp=true
      # The profiles the suite exercised before (case14: profile_assisted,
      # small_default) never produced that pair, which is why nothing caught
      # it. warmup_casePST is the case that classifies as hard_case, and it is
      # tracked, so this runs everywhere.
      cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
      case = joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m")
      net = redirect_stdout(devnull) do
        Sparlectra._import_sparlectra_net(case, nothing, cfg)
      end
      decision = Sparlectra.select_auto_pf_strategy(Sparlectra.collect_auto_pf_features(net))
      # the premise of the regression: if this case stops classifying as
      # hard_case the test silently stops testing anything
      @test decision.profile === :hard_case
      @test _autopf_step_control_ok(decision.options)

      pf = Sparlectra._apply_auto_pf_options(cfg.powerflow, decision.options, Set{String}())[1]
      @test _autopf_config_step_control_ok(pf)
      evidence = Sparlectra._auto_pf_qlimit_evidence(net)
      for stage in Sparlectra._auto_pf_escalation_stages(pf)
        open_gate, _ = stage.gate(evidence, pf)
        open_gate || continue
        pf = stage.transform(pf, evidence)[1]
        # the tuple form names the offending stage in the failure output
        @test (stage.id, _autopf_config_step_control_ok(pf)) == (stage.id, true)
      end
    end)() end

    @testset "DC fallback: honest labeling, off by default" begin (function ()
      # the service path cannot override power_flow.flatstart (not a GUI
      # key), so forcing real flat-start non-convergence at max_iter=1
      # (plus an unreachable tol, see above) needs a start_state-free copy
      # of the shipped case
      case = joinpath(mktempdir(), "sp_case14_flat.scf.json")
      raw = Sparlectra.scf_json_parse(read(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")), String))
      delete!(raw["sparlectra"], "start_state")
      write(case, Sparlectra.scf_json_string(raw))
      root = mktempdir()
      # enabled: the ladder ends in the labeled DC approximation
      rdc = start_powerflow_run(Dict{String,Any}("casefile" => case, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "config_overrides" => Dict{String,Any}("power_flow.mode" => "auto", "power_flow.max_iter" => 1, "power_flow.tol" => 1e-15, "power_flow.dc.fallback" => true)))
      @test rdc["status"] != "succeeded"
      @test rdc["converged"] === false
      @test rdc["metadata"]["dc_fallback_solution"] === true
      @test occursin("DC fallback", rdc["message"])
      # default: clean failure without a DC result
      rnf = start_powerflow_run(Dict{String,Any}("casefile" => case, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root, "config_overrides" => Dict{String,Any}("power_flow.mode" => "auto", "power_flow.max_iter" => 1, "power_flow.tol" => 1e-15)))
      @test rnf["metadata"]["dc_fallback_solution"] === false
      @test occursin("skipped_dc_fallback_disabled", read(joinpath(rnf["output_dir"], "auto_mode_decision.log"), String))
    end)() end

  return true
end

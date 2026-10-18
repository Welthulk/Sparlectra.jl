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

# file: test/test_legacy_dc_powerflow.jl
# purpose: legacy profile: the testsets of test/test_dc_powerflow.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_dc_powerflow.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_dc_powerflow_tests()
  # cut from testset: DC power flow (rundcpf!)
    @testset "Property: net injection balances to zero (lossless)" begin (function ()
      # load_fixture_net: the lossless-balance property is case-agnostic, so
      # it runs on tracked/shipped cases instead of downloaded IEEE cases (a
      # fresh install downloads nothing); warmup_casePST doubles as the
      # closed-link regression case below
      nets = (("warmup_casePST", load_fixture_net("warmup_casePST")), ("sp_case14", load_fixture_net("sp_case14")), ("sp_case5", load_fixture_net("sp_case5")))
      for (label, net) in nets
        report = rundcpf!(net)
        total = sum(row.p_gen_MW - row.p_load_MW for row in report.nodes)
        @test (label, isapprox(total, 0.0; atol = 1e-6)) == (label, true)
      end
    end)() end

    @testset "closed busbar couplers are contracted (link merge regression)" begin (function ()
      # earlier the DC path skipped _merged_pf_net: the two
      # sections of a CLOSED coupler were separate DC nodes, which either
      # split off a false reference-less island (warmup_casePST, its second
      # section hangs on the link alone) or silently solved the sections
      # 3.1/13.6 deg apart (sp_case60/sp_case188). The contraction makes
      # them one node: identical angles by construction.
      net60 = load_fixture_net("sp_case60")
      rep60 = rundcpf!(net60)
      va = Dict(row.bus_name => row.va_deg for row in rep60.nodes)
      @test va["Neubach_110"] == va["Neubach_110b"]
      # warmup_casePST solving AT ALL is the reference-less-island
      # regression; the balance property above already asserts its numbers
      @test rundcpf!(load_fixture_net("warmup_casePST")).metadata.converged
    end)() end

    @testset "dc_seed_unconditional: config-driven DC-seeded AC Newton-Raphson solve" begin (function ()
      cfg = Sparlectra.SparlectraConfig(powerflow = Sparlectra.PowerFlowConfig(start_mode = Sparlectra.StartModeConfig(dc_seed_unconditional = true)), output = OutputConfig(logfile_results = :off))
      result = run_sparlectra(casefile = "sp_case5.scf.json", path = joinpath(dirname(@__DIR__), "data", "scf"), config = cfg)
      @test result.method === :rectangular
      @test result.final_converged
      @test result.diagnostics.solver === :rectangular
      # net holds the AC-converged solution: voltage magnitudes are not flattened
      # to 1.0 everywhere the way a pure DC solve would leave them.
      @test any(!isapprox(n._vm_pu, 1.0; atol = 1e-9) for n in result.net.nodeVec)
      # The DC pre-seed step itself still ran and is recorded (mirrors
      # rundcpf!(seed_ac_start=true)'s own dc_pf_status/rectangular_pf_status
      # coexistence) -- it just isn't the run's final answer.
      @test dc_pf_status(result.net) !== nothing
      @test dc_pf_status(result.net).solver === :dc
      @test Sparlectra.rectangular_pf_status(result.net) !== nothing

      # The seed under a flat start (#463: the flat start no longer switches
      # it off): the solve must start from the seeded angles, not rebuild
      # 1.0 pu / 0 degrees over them. The seed is a flat-magnitude start
      # either way, so the first mismatch equals the seed's from the stored
      # profile and differs from the bare flat start's.
      seed_start(flat, seed) = (r = run_sparlectra(casefile = "sp_case14.scf.json", path = joinpath(dirname(@__DIR__), "data", "scf"), config = Sparlectra.SparlectraConfig(powerflow = Sparlectra.PowerFlowConfig(start_mode = Sparlectra.StartModeConfig(flatstart = flat, dc_seed_unconditional = seed)), output = OutputConfig(logfile_results = :off))); (r.final_converged, Sparlectra.rectangular_pf_status(r.net).initial_mismatch))
      conv_flat_seed, mis_flat_seed = seed_start(true, true)
      conv_stored_seed, mis_stored_seed = seed_start(false, true)
      conv_flat, mis_flat = seed_start(true, false)
      @test conv_flat_seed && conv_stored_seed && conv_flat
      @test mis_flat_seed == mis_stored_seed
      @test !isapprox(mis_flat_seed, mis_flat; rtol = 1e-3)

      # Regression: default (dc_seed_unconditional=false) still converges the same case.
      cfg_default = Sparlectra.SparlectraConfig(output = OutputConfig(logfile_results = :off))
      result_default = run_sparlectra(casefile = "sp_case5.scf.json", path = joinpath(dirname(@__DIR__), "data", "scf"), config = cfg_default)
      @test result_default.final_converged
      @test result_default.diagnostics.solver === :rectangular
    end)() end

    @testset "One specification source for AC and DC (#323)" begin (function ()
      # the DC injections come from the prosumers (buildComplexSVec), the
      # same source the rectangular AC solver reads. Two directions:
      # a node-sum edit moves NEITHER solver, a prosumer edit moves BOTH.
      base = _dc3_net()
      rep0 = rundcpf!(base)
      @test rep0.metadata.converged
      pf0 = [br.fBranchFlow.pFlow for br in base.branchVec]

      # node-sum helper: report layer only, DC result identical
      nsum = _dc3_net()
      addBusLoadPower!(net = nsum, busName = "B3", p = 30.0, q = 5.0)
      rep1 = rundcpf!(nsum)
      @test rep1.metadata.converged
      pf1 = [br.fBranchFlow.pFlow for br in nsum.branchVec]
      @test pf1 == pf0

      # prosumer edit: both solvers see the heavier load
      pros = _dc3_net()
      addProsumer!(net = pros, busName = "B3", type = "ENERGYCONSUMER", p = 30.0, q = 5.0)
      refreshBusTypesFromProsumers!(pros)
      rep2 = rundcpf!(pros)
      @test rep2.metadata.converged
      pf2 = [br.fBranchFlow.pFlow for br in pros.branchVec]
      @test pf2 != pf0
      # the DC spec equals the real part of the AC specification vector
      S = Sparlectra.buildComplexSVec(pros)
      p_b3 = real(S[geNetBusIdx(net = pros, busName = "B3")]) * pros.baseMVA
      @test p_b3 ≈ -(45.0 + 30.0) atol = 1e-9
      pros_ac = _dc3_net()
      addProsumer!(net = pros_ac, busName = "B3", type = "ENERGYCONSUMER", p = 30.0, q = 5.0)
      refreshBusTypesFromProsumers!(pros_ac)
      _, erg_ac = runpf!(pros_ac, 25, 1e-8, 0)
      @test erg_ac == 0
      @test real(Sparlectra.buildComplexSVec(pros_ac)) == real(S)
    end)() end

  return true
end

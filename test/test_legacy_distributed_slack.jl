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

# file: test/test_legacy_distributed_slack.jl
# purpose: legacy profile: the testsets of test/test_distributed_slack.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_distributed_slack.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_distributed_slack_tests()
  # cut from testset: Distributed slack (#192)
  (function ()
      # Classical baseline shared by several sub-tests.
      net0 = _dslack_testnet()
      _, erg0 = runpf_rectangular!(net0, 30, 1e-8, 0)
      vm0, va0 = _dslack_vmva(net0)

      @testset "off-path regression: disabled is bit-identical" begin (function ()
        net1 = _dslack_testnet()
        _, erg1 = runpf_rectangular!(net1, 30, 1e-8, 0; distributed_slack_enabled = false)
        @test erg1 == 0
        vm1, va1 = _dslack_vmva(net1)
        @test vm1 == vm0
        @test va1 == va0
        st = Sparlectra.rectangular_pf_status(net1)
        @test st.distributed_slack_active == false
        @test st.distributed_slack_lambda_p_mw == 0.0
      end)() end

      @testset "Q-limit interaction: PV→PQ switch keeps participation" begin (function ()
        # Tight Q limit forces B2 PV→PQ during the solve; its alpha must stay.
        net7 = _dslack_testnet(qmax_b2 = 2.0)
        _, erg7 = runpf_rectangular!(net7, 40, 1e-8, 0; distributed_slack_enabled = true, distributed_slack_p_mode = :pg_weighted)
        @test erg7 == 0
        st = Sparlectra.rectangular_pf_status(net7)
        @test st.distributed_slack_active == true
        @test st.distributed_slack_participants == 2
        @test isapprox(st.distributed_slack_alpha_sum, 1.0; atol = 1e-12)
        resid = _dslack_p_residual(net7)
        lambda_pu = st.distributed_slack_lambda_p_pu
        # B2 kept its 0.6 share although it no longer holds its voltage.
        @test isapprox(resid[2], 0.6 * lambda_pu; atol = 1e-7)
        @test isapprox(resid[3], 0.4 * lambda_pu; atol = 1e-7)
        # And the switch really happened: B2 logged a Q-limit hit.
        @test Sparlectra.qlimit_switch_count(net7, 2) >= 1
      end)() end

      @testset "config: default placeholder and block-style weights load" begin (function ()
        # Regression: the minimal YAML reader parses the default `weights: {}`
        # placeholder as the literal string "{}"; the config constructor must
        # treat it (and nothing/"") as an empty table instead of throwing —
        # otherwise EVERY default-merged load_sparlectra_config call fails.
        p0 = test_scratch_path(".yaml")
        write(p0, "config_version: 1\npower_flow:\n  tol: 1.0e-8\n")
        cfg = Sparlectra.load_sparlectra_config(p0; reload = true)
        @test cfg.powerflow.distributed_slack.enabled == false
        @test isempty(cfg.powerflow.distributed_slack.weights)
        # Real weight tables use block style; their child keys are user data
        # (bus names) and must pass unknown-key validation.
        p1 = test_scratch_path(".yaml")
        write(p1, """
  config_version: 1
  power_flow:
    distributed_slack:
      enabled: true
      p_mode: explicit
      weights:
        B1: 2.0
        B2: 1.0
  """)
        cfg1 = Sparlectra.load_sparlectra_config(p1; reload = true)
        ds = cfg1.powerflow.distributed_slack
        @test ds.enabled && ds.p_mode == :explicit
        @test ds.weights == Dict("B1" => 2.0, "B2" => 1.0)
      end)() end

      @testset "two islands get independent lambda_P" begin (function ()
        net8 = _dslack_two_island_net()
        profile = Dict{Symbol,Any}()
        _, erg8 = runpf!(
          net8, 40, 1e-8, 0;
          islands_enabled = true,
          distributed_slack_enabled = true,
          distributed_slack_p_mode = :pg_weighted,
          performance_profile = profile,
        )
        @test erg8 == 0
        statuses = profile[:ac_island_solver_statuses]
        @test length(statuses) == 2
        lambdas = sort([st.distributed_slack_lambda_p_mw for st in values(statuses)])
        @test all(l -> l > 0.0, lambdas)
        # Island A misses ~20 MW + losses, island C ~10 MW + losses: distinct.
        @test lambdas[2] - lambdas[1] > 5.0
        for st in values(statuses)
          @test st.distributed_slack_active == true
          @test st.distributed_slack_participants == 1
        end
        # The top-level status aggregates over the islands (regression: it used
        # to carry an arbitrary single island's lambda_P instead).
        agg = Sparlectra.rectangular_pf_status(net8)
        @test agg.distributed_slack_active == true
        @test agg.distributed_slack_participants == 2
        @test isapprox(agg.distributed_slack_lambda_p_mw, sum(st.distributed_slack_lambda_p_mw for st in values(statuses)); atol = 1e-9)
        # The aggregate participation holds both islands in the numbering of
        # the net (it held one island's rows on buses 1 and 2). A2 and C2
        # deliver schedule plus share in the report rows.
        part = agg.distributed_slack_participation
        @test sort([(r.bus_idx, r.gen_index) for r in part]) == [(2, 2), (5, 5)]
        rows = buildACPFlowReport(net8).nodes
        for r in part
          @test isapprox(rows[r.bus_idx].p_gen_MW, r.pg_mw + r.dp_mw; atol = 1e-9)
        end
      end)() end

  end)()
  return true
end

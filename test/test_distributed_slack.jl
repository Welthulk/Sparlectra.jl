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
#
# file: test/test_distributed_slack.jl
# purpose: regression and feature tests for the distributed active-power
# slack (issue #192): off-path identity, conservation, REF equivalence,
# weight modes (pg/explicit/imported via MATPOWER APF), fallback behavior,
# Q-limit interaction, and independent per-island lambda_P.

# Small 4-bus test net: REF (external injection) + two PV generators + one
# PQ load. The deliberate P imbalance (70 MW load vs 50 MW scheduled
# generation plus losses) is what lambda_P must absorb.
function _dslack_testnet(; qmax_b2::Float64 = 50.0)
  net = Net(name = "dslack_test", baseMVA = 100.0)
  for b in ("B1", "B2", "B3", "B4")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.02, x_pu = 0.10, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.01, x_pu = 0.06, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B4", r_pu = 0.03, x_pu = 0.12, b_pu = 0.02, status = 1)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  addProsumer!(net = net, busName = "B2", type = "GENERATOR", p = 30.0, vm_pu = 1.01, qMin = -50.0, qMax = qmax_b2, pMax = 100.0, pMin = 0.0)
  addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 20.0, vm_pu = 1.0, qMin = -50.0, qMax = 50.0, pMax = 100.0, pMin = 0.0)
  addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 70.0, q = 20.0)
  return net
end

_dslack_vmva(net) = ([n._vm_pu for n in net.nodeVec], [n._va_deg for n in net.nodeVec])

# Per-bus P residual P_calc(V) - P_spec at the solved voltages. With
# distributed slack this must equal alpha .* lambda_P; classically it is the
# whole imbalance sitting on the REF bus.
function _dslack_p_residual(net)
  Y = Sparlectra.createYBUS(net = net, sparse = true, printYBUS = false)
  V = [n._vm_pu * cis(deg2rad(n._va_deg)) for n in net.nodeVec]
  return real.(Sparlectra.calc_injections(Y, V)) .- real.(Sparlectra.buildComplexSVec(net))
end

# MATPOWER fixture with the optional 21st gen column (APF). `apf2`/`apf3` are
# the participation factors of the two PV generators; the slack unit stays 0.
function _dslack_write_matpower_case(dir::String, name::String; apf2::Float64, apf3::Float64)
  path = joinpath(dir, name * ".m")
  write(path, """
function mpc = $(name)
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
1 3 0 0 0 0 1 1.00 0 110 1 1.1 0.9;
2 2 0 0 0 0 1 1.01 0 110 1 1.1 0.9;
3 2 0 0 0 0 1 1.00 0 110 1 1.1 0.9;
4 1 70 20 0 0 1 1.00 0 110 1 1.1 0.9;
];
mpc.gen = [
1 0 0 300 -300 1.00 100 1 250 -250 0 0 0 0 0 0 0 0 0 0 $(0.0);
2 30 0 50 -50 1.01 100 1 100 0 0 0 0 0 0 0 0 0 0 0 $(apf2);
3 20 0 50 -50 1.00 100 1 100 0 0 0 0 0 0 0 0 0 0 0 $(apf3);
];
mpc.branch = [
1 2 0.01 0.08 0.02 0 0 0 0 0 1 -360 360;
2 3 0.02 0.10 0.02 0 0 0 0 0 1 -360 360;
3 4 0.01 0.06 0.02 0 0 0 0 0 1 -360 360;
1 4 0.03 0.12 0.02 0 0 0 0 0 1 -360 360;
];
""")
  return path
end

# Two electrically disconnected islands, each with its own reference and one
# PV participant, so the per-island lambda_P values must be independent.
function _dslack_two_island_net()
  net = Net(name = "dslack_islands", baseMVA = 100.0)
  for b in ("A1", "A2", "A3", "C1", "C2", "C3")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "A2", toBus = "A3", r_pu = 0.02, x_pu = 0.10, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "C1", toBus = "C2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "C2", toBus = "C3", r_pu = 0.02, x_pu = 0.10, b_pu = 0.02, status = 1)
  addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "A1")
  addProsumer!(net = net, busName = "A2", type = "GENERATOR", p = 30.0, vm_pu = 1.0, qMin = -50.0, qMax = 50.0)
  addProsumer!(net = net, busName = "A3", type = "ENERGYCONSUMER", p = 50.0, q = 10.0)
  addProsumer!(net = net, busName = "C1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "C1")
  addProsumer!(net = net, busName = "C2", type = "GENERATOR", p = 10.0, vm_pu = 1.0, qMin = -50.0, qMax = 50.0)
  addProsumer!(net = net, busName = "C3", type = "ENERGYCONSUMER", p = 20.0, q = 5.0)
  return net
end

function run_distributed_slack_tests()
  @testset "Distributed slack (#192)" begin
    # Classical baseline shared by several sub-tests.
    net0 = _dslack_testnet()
    _, erg0 = runpf_rectangular!(net0, 30, 1e-8, 0)
    @test erg0 == 0
    vm0, va0 = _dslack_vmva(net0)
    resid0 = _dslack_p_residual(net0)
    # Classically the REF bus absorbs the whole imbalance alone.
    @test abs(sum(resid0) - resid0[1]) < 1e-7

    @testset "off-path regression: disabled is bit-identical" begin
      net1 = _dslack_testnet()
      _, erg1 = runpf_rectangular!(net1, 30, 1e-8, 0; distributed_slack_enabled = false)
      @test erg1 == 0
      vm1, va1 = _dslack_vmva(net1)
      @test vm1 == vm0
      @test va1 == va0
      st = Sparlectra.rectangular_pf_status(net1)
      @test st.distributed_slack_active == false
      @test st.distributed_slack_lambda_p_mw == 0.0
    end

    @testset "conservation: residual pattern equals alpha*lambda" begin
      net2 = _dslack_testnet()
      _, erg2 = runpf_rectangular!(net2, 30, 1e-8, 0; distributed_slack_enabled = true, distributed_slack_p_mode = :pg_weighted)
      @test erg2 == 0
      st = Sparlectra.rectangular_pf_status(net2)
      @test st.distributed_slack_active == true
      @test st.distributed_slack_mode == :pg_weighted
      @test isapprox(st.distributed_slack_alpha_sum, 1.0; atol = 1e-12)
      lambda_pu = st.distributed_slack_lambda_p_pu
      @test lambda_pu > 0.0
      resid = _dslack_p_residual(net2)
      # alpha = (0.6, 0.4) from pg = (30, 20); REF and load keep their spec.
      @test isapprox(resid[2], 0.6 * lambda_pu; atol = 1e-7)
      @test isapprox(resid[3], 0.4 * lambda_pu; atol = 1e-7)
      @test abs(resid[1]) < 1e-7
      @test abs(resid[4]) < 1e-7
      # Total distributed correction matches the classically REF-absorbed
      # imbalance up to the (small) loss shift of the changed dispatch.
      @test isapprox(sum(resid), resid0[1]; atol = 2e-3)

      # Reporting surface: the participant table is persisted in the status
      # and rendered by the classical result print (bus, alpha, dP).
      part = st.distributed_slack_participation
      @test length(part) == 2
      @test isapprox(sort([r.alpha for r in part]), [0.4, 0.6]; atol = 1e-12)
      @test isapprox(sum(r.dp_mw for r in part), st.distributed_slack_lambda_p_mw; atol = 1e-9)
      buf = IOBuffer()
      Sparlectra._print_distributed_slack_participation(buf, net2)
      out = String(take!(buf))
      @test occursin("Distributed slack: mode pg_weighted", out)
      @test occursin("B2", out)
      @test occursin("Pg eff [MW]", out)
    end

    @testset "equivalence: all weight on the REF bus reproduces classical" begin
      net3 = _dslack_testnet()
      _, erg3 = runpf_rectangular!(
        net3, 30, 1e-10, 0;
        distributed_slack_enabled = true,
        distributed_slack_p_mode = :explicit,
        distributed_slack_weights = Dict("B1" => 1.0),
      )
      @test erg3 == 0
      vm3, va3 = _dslack_vmva(net3)
      @test isapprox(vm3, vm0; atol = 1e-8)
      @test isapprox(va3, va0; atol = 1e-8)
      st = Sparlectra.rectangular_pf_status(net3)
      # lambda_P is exactly the classical REF-absorbed power.
      @test isapprox(st.distributed_slack_lambda_p_pu, resid0[1]; atol = 1e-7)
    end

    @testset "modes: pg_weighted == matching explicit weights" begin
      net4a = _dslack_testnet()
      _, erg4a = runpf_rectangular!(net4a, 30, 1e-8, 0; distributed_slack_enabled = true, distributed_slack_p_mode = :pg_weighted)
      net4b = _dslack_testnet()
      _, erg4b = runpf_rectangular!(
        net4b, 30, 1e-8, 0;
        distributed_slack_enabled = true,
        distributed_slack_p_mode = :explicit,
        distributed_slack_weights = Dict("B2" => 30.0, "B3" => 20.0),
      )
      @test erg4a == 0 && erg4b == 0
      @test isapprox(first(_dslack_vmva(net4a)), first(_dslack_vmva(net4b)); atol = 1e-10)
      @test isapprox(last(_dslack_vmva(net4a)), last(_dslack_vmva(net4b)); atol = 1e-10)
    end

    @testset "modes: imported via MATPOWER APF column" begin
      dir = mktempdir()
      # APF present, solver disabled: identical to a zero-APF twin (carrying
      # the factor must never change results on its own).
      path_apf = _dslack_write_matpower_case(dir, "case4_apf", apf2 = 0.75, apf3 = 0.25)
      path_zero = _dslack_write_matpower_case(dir, "case4_apf_zero", apf2 = 0.0, apf3 = 0.0)
      net_apf = Sparlectra.createNetFromMatPowerFile(filename = path_apf)
      net_zero = Sparlectra.createNetFromMatPowerFile(filename = path_zero)
      # import-side arrival: positive APF lands, zero maps to unknown
      pf_by_gen = [ps.participationFactor for ps in net_apf.prosumpsVec if Sparlectra.isGenerator(ps)]
      @test any(==(0.75), skipmissing(replace(pf_by_gen, nothing => missing)))
      @test any(==(0.25), skipmissing(replace(pf_by_gen, nothing => missing)))
      @test all(pf === nothing for pf in (ps.participationFactor for ps in net_zero.prosumpsVec if Sparlectra.isGenerator(ps)))
      _, erga = runpf_rectangular!(net_apf, 30, 1e-8, 0)
      _, ergz = runpf_rectangular!(net_zero, 30, 1e-8, 0)
      @test erga == 0 && ergz == 0
      @test first(_dslack_vmva(net_apf)) == first(_dslack_vmva(net_zero))
      @test last(_dslack_vmva(net_apf)) == last(_dslack_vmva(net_zero))
      # enabled with :imported — shares proportional to APF (0.75 : 0.25)
      net_imp = Sparlectra.createNetFromMatPowerFile(filename = path_apf)
      _, ergi = runpf_rectangular!(net_imp, 30, 1e-8, 0; distributed_slack_enabled = true, distributed_slack_p_mode = :imported)
      @test ergi == 0
      resid = _dslack_p_residual(net_imp)
      @test isapprox(resid[2] / resid[3], 3.0; atol = 1e-4)
    end

    @testset "fallback: error throws, ref_only matches classical" begin
      # :imported with no factors anywhere -> no valid participant
      net5 = _dslack_testnet()
      @test_throws ArgumentError runpf_rectangular!(
        net5, 30, 1e-8, 0;
        distributed_slack_enabled = true,
        distributed_slack_p_mode = :imported,
        distributed_slack_fallback = :error,
      )
      net6 = _dslack_testnet()
      it6 = erg6 = nothing
      @test_logs (:warn, r"no valid participant") match_mode = :any begin
        it6, erg6 = runpf_rectangular!(
          net6, 30, 1e-8, 0;
          distributed_slack_enabled = true,
          distributed_slack_p_mode = :imported,
          distributed_slack_fallback = :ref_only,
        )
      end
      @test erg6 == 0
      @test first(_dslack_vmva(net6)) == vm0
      @test last(_dslack_vmva(net6)) == va0
      st = Sparlectra.rectangular_pf_status(net6)
      @test st.distributed_slack_active == false
      @test isempty(st.distributed_slack_participation)
      # inactive run: the print block stays silent
      buf = IOBuffer()
      Sparlectra._print_distributed_slack_participation(buf, net6)
      @test isempty(String(take!(buf)))
    end

    @testset "Q-limit interaction: PV→PQ switch keeps participation" begin
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
    end

    @testset "config: default placeholder and block-style weights load" begin
      # Regression: the minimal YAML reader parses the default `weights: {}`
      # placeholder as the literal string "{}"; the config constructor must
      # treat it (and nothing/"") as an empty table instead of throwing —
      # otherwise EVERY default-merged load_sparlectra_config call fails.
      p0 = tempname() * ".yaml"
      write(p0, "power_flow:\n  tol: 1.0e-8\n")
      cfg = Sparlectra.load_sparlectra_config(p0; reload = true)
      @test cfg.powerflow.distributed_slack.enabled == false
      @test isempty(cfg.powerflow.distributed_slack.weights)
      # Real weight tables use block style; their child keys are user data
      # (bus names) and must pass unknown-key validation.
      p1 = tempname() * ".yaml"
      write(p1, """
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
    end

    @testset "two islands get independent lambda_P" begin
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
    end
  end
end

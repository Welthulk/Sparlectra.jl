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

# Hand case with lossless lines: reference A scheduled 30 MW, PV B 20 MW,
# load C 80 MW. pg_weighted gives alpha 0.6/0.4 and lambda_P = 30 MW, so A
# delivers 48 MW and B 32 MW. With `link = true` B sits behind a closed bus
# link at A (bus "A2"), so the solver sees both units at one bus.
function _dslack_hand_net(; link::Bool = false)
  b = link ? "A2" : "B"
  net = Net(name = link ? "dslack_hand_link" : "dslack_hand", baseMVA = 100.0)
  for name in ("A", b, "C")
    addBus!(net = net, busName = name, vn_kV = 110.0)
  end
  addPIModelACLine!(net = net, fromBus = "A", toBus = "C", r_pu = 0.0, x_pu = 0.1, b_pu = 0.0, status = 1)
  if link
    addLink!(net = net, fromBus = "A", toBus = b, status = 1)
  else
    addPIModelACLine!(net = net, fromBus = b, toBus = "C", r_pu = 0.0, x_pu = 0.1, b_pu = 0.0, status = 1)
  end
  addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", p = 30.0, vm_pu = 1.0, va_deg = 0.0, referencePri = "A")
  addProsumer!(net = net, busName = b, type = "SYNCHRONOUSMACHINE", p = 20.0, vm_pu = 1.0)
  addProsumer!(net = net, busName = "C", type = "ENERGYCONSUMER", p = 80.0, q = 0.0)
  return net
end

# Pg and Pg eff cells (MW) of one bus row of the printed bus table, NaN for
# an empty cell.
function _dslack_table_pg(rendered::String, bus::String)
  lines = split(rendered, '\n')
  header = strip.(split(lines[findfirst(l -> occursin("Pg eff MW", l), lines)], '|'))
  row = strip.(split(lines[findfirst(l -> length(split(l, '|')) > 3 && strip(split(l, '|')[3]) == bus, lines)], '|'))
  cell(name) = something(tryparse(Float64, row[findfirst(==(name), header)]), NaN)
  return cell("Pg [MW]"), cell("Pg eff MW")
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
  @testset "Distributed slack (#192)" begin (function ()
    # Classical baseline shared by several sub-tests.
    net0 = _dslack_testnet()
    _, erg0 = runpf_rectangular!(net0, 30, 1e-8, 0)
    @test erg0 == 0
    vm0, va0 = _dslack_vmva(net0)
    resid0 = _dslack_p_residual(net0)
    # Classically the REF bus absorbs the whole imbalance alone.
    @test abs(sum(resid0) - resid0[1]) < 1e-7

    @testset "conservation: residual pattern equals alpha*lambda" begin (function ()
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
      # and gives the per-bus shares of the bus table.
      part = st.distributed_slack_participation
      @test length(part) == 2
      @test isapprox(sort([r.alpha for r in part]), [0.4, 0.6]; atol = 1e-12)
      @test isapprox(sum(r.dp_mw for r in part), st.distributed_slack_lambda_p_mw; atol = 1e-9)
      active, shares = Sparlectra._distributed_slack_bus_shares(net2)
      @test active
      @test length(shares) == 2
      @test isapprox(sort([a for (a, _) in values(shares)]), [0.4, 0.6]; atol = 1e-12)

      # Hand case: Pg is the schedule and Pg eff the delivered generation at
      # every participating bus, the reference included (the share was added
      # to the solved slack value, 66 MW); the report rows and the units
      # carry the delivered values (B had 20 MW, the units summed to 68).
      for link in (false, true)
        net_h = _dslack_hand_net(link = link)
        b = link ? "A2" : "B"
        _, erg_h = runpf!(net_h, 30, 1e-10, 0; distributed_slack_enabled = true, distributed_slack_p_mode = :pg_weighted)
        @test erg_h == 0
        @test isapprox(Sparlectra.rectangular_pf_status(net_h).distributed_slack_lambda_p_mw, 30.0; atol = 1e-6)
        tmp_h = mktempdir()
        printACPFlowResults(net_h, 0.0, 1, 1e-10, true, tmp_h)
        rendered_h = read(joinpath(tmp_h, "result_$(net_h.name).txt"), String)
        pg_a, eff_a = _dslack_table_pg(rendered_h, "A")
        pg_b, eff_b = _dslack_table_pg(rendered_h, b)
        @test (pg_a, eff_a, pg_b, eff_b) == (30.0, 48.0, 20.0, 32.0)
        rows = buildACPFlowReport(net_h).nodes
        @test isapprox([r.p_gen_MW for r in rows], [48.0, 32.0, 0.0]; atol = 1e-6)
        distributeBusResults!(net_h)
        @test isapprox([ps.pRes for ps in net_h.prosumpsVec], [48.0, 32.0, -80.0]; atol = 1e-6)
      end
    end)() end

    @testset "explicit weights reproduce the state of the matching mode" begin (function ()
      # pg_weighted reference: Pg 30/20 MW gives the shares 0.6/0.4
      net_pg = _dslack_testnet()
      _, erg_pg = runpf_rectangular!(net_pg, 30, 1e-8, 0; distributed_slack_enabled = true, distributed_slack_p_mode = :pg_weighted)
      @test erg_pg == 0
      vm_pg, va_pg = _dslack_vmva(net_pg)
      lambda_pg = Sparlectra.rectangular_pf_status(net_pg).distributed_slack_lambda_p_pu
      # weights -> state and lambda_P they must reproduce: all weight on the
      # REF bus gives the classical state, and lambda_P is the classically
      # REF-absorbed power; weights 30/20 give the pg_weighted state
      for (weights, tol, vm_ref, va_ref, lambda_ref, atol) in (
        (Dict("B1" => 1.0), 1e-10, vm0, va0, resid0[1], 1e-8),
        (Dict("B2" => 30.0, "B3" => 20.0), 1e-8, vm_pg, va_pg, lambda_pg, 1e-10),
      )
        net = _dslack_testnet()
        _, erg = runpf_rectangular!(
          net, 30, tol, 0;
          distributed_slack_enabled = true,
          distributed_slack_p_mode = :explicit,
          distributed_slack_weights = weights,
        )
        @test erg == 0
        vm, va = _dslack_vmva(net)
        @test isapprox(vm, vm_ref; atol = atol)
        @test isapprox(va, va_ref; atol = atol)
        @test isapprox(Sparlectra.rectangular_pf_status(net).distributed_slack_lambda_p_pu, lambda_ref; atol = 1e-7)
      end
    end)() end

    @testset "modes: imported via MATPOWER APF column" begin (function ()
      dir = mktempdir()
      path_apf = _dslack_write_matpower_case(dir, "case4_apf", apf2 = 0.75, apf3 = 0.25)
      path_zero = _dslack_write_matpower_case(dir, "case4_apf_zero", apf2 = 0.0, apf3 = 0.0)
      net_apf = Sparlectra.createNetFromMatPowerFile(filename = path_apf)
      net_zero = Sparlectra.createNetFromMatPowerFile(filename = path_zero)
      # import-side arrival: positive APF lands, zero maps to unknown
      pf_by_gen = [ps.participationFactor for ps in net_apf.prosumpsVec if Sparlectra.isGenerator(ps)]
      @test any(==(0.75), skipmissing(replace(pf_by_gen, nothing => missing)))
      @test any(==(0.25), skipmissing(replace(pf_by_gen, nothing => missing)))
      @test all(pf === nothing for pf in (ps.participationFactor for ps in net_zero.prosumpsVec if Sparlectra.isGenerator(ps)))
      # enabled with :imported: shares proportional to APF (0.75 : 0.25)
      net_imp = Sparlectra.createNetFromMatPowerFile(filename = path_apf)
      _, ergi = runpf_rectangular!(net_imp, 30, 1e-8, 0; distributed_slack_enabled = true, distributed_slack_p_mode = :imported)
      @test ergi == 0
      resid = _dslack_p_residual(net_imp)
      @test isapprox(resid[2] / resid[3], 3.0; atol = 1e-4)
    end)() end

    @testset "fallback: error throws, ref_only matches classical" begin (function ()
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
      # inactive run: summary line silent, no participation columns
      buf = IOBuffer()
      Sparlectra._print_distributed_slack_summary_line(buf, net6)
      @test isempty(String(take!(buf)))
      active, shares = Sparlectra._distributed_slack_bus_shares(net6)
      @test !active
      @test isempty(shares)
    end)() end

    @testset "Q-limit interaction: PV->PQ switch keeps participation" begin (function ()
      # Tight Q limit forces B2 PV->PQ during the solve; its alpha must stay.
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
      # The top-level status aggregates over the islands (it once carried
      # an arbitrary single island's lambda_P instead).
      agg = Sparlectra.rectangular_pf_status(net8)
      @test agg.distributed_slack_active == true
      @test agg.distributed_slack_participants == 2
      @test isapprox(agg.distributed_slack_lambda_p_mw, sum(st.distributed_slack_lambda_p_mw for st in values(statuses)); atol = 1e-9)
      # The aggregate participation holds both islands in the numbering of
      # the net. A2 and C2 deliver schedule plus share in the report rows.
      part = agg.distributed_slack_participation
      @test sort([(r.bus_idx, r.gen_index) for r in part]) == [(2, 2), (5, 5)]
      rows = buildACPFlowReport(net8).nodes
      for r in part
        @test isapprox(rows[r.bus_idx].p_gen_MW, r.pg_mw + r.dp_mw; atol = 1e-9)
      end
    end)() end

  end)() end
end

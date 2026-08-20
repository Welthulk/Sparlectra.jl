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

# file: test/test_hvdc_pair_control.jl
# purpose: HVDC back-to-back pairing controller (#297 Draft B): registration
#          validation, exact pairing invariant P_to = transfer - loss on a
#          two-island fixture, rating clamp with honest at_limit, per-side
#          voltage-target secant, bit-identical baseline with a disabled
#          controller, and the element-row vocabulary.

# Two areas coupled ONLY by the converter pair: island 1 = B1 (slack) - B2,
# island 2 = B3 (slack) - B4. The converters sit at B2 (from) and B4 (to)
# as PQ generator injections seeded with a Stage-0-like fixed operating
# point (equality clamps included, so registration must unclamp them).
function _build_hvdc_two_area_net(name::String; seed_from_p::Float64 = -80.0, seed_to_p::Float64 = 76.0)
  net = Net(name = name, baseMVA = 100.0)
  for b in ("B1", "B2", "B3", "B4")
    addBus!(net = net, busName = b, vn_kV = 380.0)
  end
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "B3", type = "EXTERNALNETWORKINJECTION", referencePri = "B3", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
  addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 50.0, q = 12.0)
  # converter injections (from side exports into the link, to side receives)
  addProsumer!(net = net, busName = "B2", type = "GENERATOR", p = seed_from_p, q = 0.0, pMin = seed_from_p, pMax = seed_from_p, qMin = 0.0, qMax = 0.0)
  from_idx = length(net.prosumpsVec)
  addProsumer!(net = net, busName = "B4", type = "GENERATOR", p = seed_to_p, q = 0.0, pMin = seed_to_p, pMax = seed_to_p, qMin = 0.0, qMax = 0.0)
  to_idx = length(net.prosumpsVec)
  return net, from_idx, to_idx
end

# island_feed fixture: island A keeps a classical reference, island C is fed
# ONLY by the grid-forming converter at C2 (modeled as the island reference);
# the island load hangs at the far end C1 so the PCC line carries real flow
function _build_hvdc_island_feed_net(name::String; sending_mw::Float64 = 0.0)
  net = Net(name = name, baseMVA = 100.0)
  for b in ("A1", "A2", "C1", "C2")
    addBus!(net = net, busName = b, vn_kV = 380.0)
  end
  addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "C2", toBus = "C1", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", referencePri = "A1", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
  addProsumer!(net = net, busName = "C2", type = "EXTERNALNETWORKINJECTION", referencePri = "C2", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "C1", type = "ENERGYCONSUMER", p = 50.0, q = 12.0)
  addProsumer!(net = net, busName = "A2", type = "GENERATOR", p = sending_mw, q = 0.0)
  return net, length(net.prosumpsVec)
end

function _hvdc_pf_cfg()
  return PowerFlowConfig(method = :rectangular, max_iter = 40, tol = 1e-9)
end

function _hvdc_run!(net)
  return run_control!(net; controllers = collect_outer_controllers(net), pf_config = _hvdc_pf_cfg(), control_config = ControlConfig(max_outer_iterations = 15, trace = true), verbose = 0)
end

function run_hvdc_pair_control_tests()
  @testset "HVDC pair controller (#297 Draft B)" begin
    @testset "registration validation" begin
      net, f_idx, t_idx = _build_hvdc_two_area_net("hvdc_val")
      @test_throws ErrorException addHvdcPairControl!(net; from_bus = "B2", to_bus = "B2", p_transfer_mw = 10.0)
      @test_throws ErrorException addHvdcPairControl!(net; from_bus = "B2", to_bus = "B4", p_transfer_mw = 10.0, loss_mw = -1.0)
      @test_throws ErrorException addHvdcPairControl!(net; from_bus = "B2", to_bus = "B4", p_transfer_mw = 10.0, loss_fraction = 1.0)
      @test_throws ErrorException addHvdcPairControl!(net; from_bus = "B2", to_bus = "B4", p_transfer_mw = 10.0, from_q_mvar = 1.0, from_vset_pu = 1.0)
      @test_throws ErrorException addHvdcPairControl!(net; from_bus = "NoSuchBus", to_bus = "B4", p_transfer_mw = 10.0)
      # the reference injection balances its island and cannot be paired
      @test_throws ErrorException addHvdcPairControl!(net; from_bus = "B1", to_bus = "B4", p_transfer_mw = 10.0)
      # voltage target without any usable reactive range (Stage-0 equality
      # clamps were lifted at resolution, leaving no range at all)
      @test_throws ErrorException addHvdcPairControl!(net; from_bus = "B2", to_bus = "B4", p_transfer_mw = 10.0, to_vset_pu = 1.0)
      ctrl = addHvdcPairControl!(net; from_bus = "B2", to_bus = "B4", p_transfer_mw = 10.0)
      @test ctrl isa HvdcPairControl
      # exclusivity: the pair's converters accept no second controller
      @test_throws ErrorException addHvdcPairControl!(net; from_bus = "B2", to_bus = "B4", p_transfer_mw = 20.0)
      clearHvdcPairControllers!(net)
      @test isempty(Sparlectra._hvdc_pair_controllers(net))
      # registration lifted the Stage-0 equality clamps on both converters
      @test net.prosumpsVec[f_idx].minP === nothing && net.prosumpsVec[f_idx].maxP === nothing
      @test net.prosumpsVec[t_idx].minP === nothing && net.prosumpsVec[t_idx].maxP === nothing
    end

    @testset "transfer target and exact pairing invariant" begin
      net, f_idx, t_idx = _build_hvdc_two_area_net("hvdc_transfer")
      addHvdcPairControl!(net; from_bus = "B2", to_bus = "B4", p_transfer_mw = 100.0, loss_mw = 4.0, loss_fraction = 0.05)
      result = _hvdc_run!(net)
      @test result.converged == true
      @test result.status == :converged
      # invariant exact: from = -transfer, to = transfer - (loss_mw + 5%)
      @test net.prosumpsVec[f_idx].pVal == -100.0
      @test net.prosumpsVec[t_idx].pVal == 100.0 - (4.0 + 0.05 * 100.0)
      ctrl = only(Sparlectra._hvdc_pair_controllers(net))
      @test ctrl.converged
      @test !ctrl.at_limit
      rows = Sparlectra.control_report_rows(ctrl, net, Sparlectra.NoControlState(), (outer_iteration = ctrl.outer_iters,))
      @test only(rows).p_transfer_mw == 100.0
      @test only(rows).loss_mw == 9.0
      @test only(rows).to_p_mw == 91.0
      # both islands solved with their own references
      @test result.powerflow_solves >= 1
    end

    @testset "reversed transfer" begin
      net, f_idx, t_idx = _build_hvdc_two_area_net("hvdc_reverse"; seed_from_p = 0.0, seed_to_p = 0.0)
      addHvdcPairControl!(net; from_bus = "B2", to_bus = "B4", p_transfer_mw = -60.0, loss_mw = 2.0)
      result = _hvdc_run!(net)
      @test result.converged == true
      @test net.prosumpsVec[f_idx].pVal == 60.0
      @test net.prosumpsVec[t_idx].pVal == -62.0
    end

    @testset "rating clamp is honest at_limit" begin
      net, f_idx, t_idx = _build_hvdc_two_area_net("hvdc_rating")
      addHvdcPairControl!(net; from_bus = "B2", to_bus = "B4", p_transfer_mw = 200.0, p_rating_mw = 100.0, loss_mw = 4.0)
      result = _hvdc_run!(net)
      ctrl = only(Sparlectra._hvdc_pair_controllers(net))
      @test net.prosumpsVec[f_idx].pVal == -100.0
      @test net.prosumpsVec[t_idx].pVal == 96.0
      @test ctrl.at_limit
      @test ctrl.status == :at_limit
      # framework semantics: a blocked controller counts as done, the run
      # result reports the finished loop; honesty lives on ctrl.at_limit
      # (same convention as the TCSC at_limit case)
      @test result.converged == true
      @test ctrl.converged == false
    end

    @testset "voltage-target side reaches the setpoint" begin
      net, f_idx, t_idx = _build_hvdc_two_area_net("hvdc_vset")
      addHvdcPairControl!(net; from_bus = "B2", to_bus = "B4", p_transfer_mw = 80.0, loss_mw = 3.0, to_vset_pu = 1.01, to_qmin_mvar = -60.0, to_qmax_mvar = 60.0)
      result = _hvdc_run!(net)
      ctrl = only(Sparlectra._hvdc_pair_controllers(net))
      @test result.converged == true
      @test ctrl.converged
      vm = get_bus_vm_pu(net, "B4")
      @test vm !== nothing && abs(vm - 1.01) <= ctrl.deadband_vm_pu
      @test -60.0 <= ctrl.to_q_now <= 60.0
      @test net.prosumpsVec[t_idx].qVal == ctrl.to_q_now
    end

    @testset "unreachable voltage target clamps honestly" begin
      net, _, _ = _build_hvdc_two_area_net("hvdc_vset_limit")
      addHvdcPairControl!(net; from_bus = "B2", to_bus = "B4", p_transfer_mw = 80.0, loss_mw = 3.0, to_vset_pu = 1.2, to_qmin_mvar = -5.0, to_qmax_mvar = 5.0)
      _hvdc_run!(net)
      ctrl = only(Sparlectra._hvdc_pair_controllers(net))
      @test !ctrl.converged
      @test ctrl.at_limit
      @test isapprox(ctrl.to_q_now, 5.0; atol = 1e-9)
    end

    @testset "disabled controller keeps the baseline bit-identical" begin
      net_a, _, _ = _build_hvdc_two_area_net("hvdc_base_a")
      net_b, _, _ = _build_hvdc_two_area_net("hvdc_base_b")
      addHvdcPairControl!(net_b; from_bus = "B2", to_bus = "B4", p_transfer_mw = 999.0, enabled = false)
      # the fixture is two islands by design; the raw legacy runpf! call
      # needs island solving enabled explicitly (PowerFlowConfig defaults
      # to enabled, the positional API does not)
      ite_a, erg_a = runpf!(net_a, 40, 1e-9, 0; method = :rectangular, islands_enabled = true)
      ite_b, erg_b = runpf!(net_b, 40, 1e-9, 0; method = :rectangular, islands_enabled = true)
      @test erg_a == 0 && erg_b == 0
      @test ite_a == ite_b
      vm_a = [something(n._vm_pu, NaN) for n in net_a.nodeVec]
      vm_b = [something(n._vm_pu, NaN) for n in net_b.nodeVec]
      @test vm_a == vm_b
    end

    @testset "YAML declaration round trip (#305 schema)" begin
      # declared via control.controllers, run through the real config path,
      # compared against the programmatic twin
      net_prog, f_idx, t_idx = _build_hvdc_two_area_net("hvdc_yaml_prog")
      addHvdcPairControl!(net_prog; from_bus = "B2", to_bus = "B4", p_transfer_mw = 90.0, loss_mw = 5.0)
      _hvdc_run!(net_prog)
      net_yaml, yf_idx, yt_idx = _build_hvdc_two_area_net("hvdc_yaml_decl")
      cfg = ControlConfig(controllers = Any[Dict{String,Any}("type" => "hvdc_pair", "name" => "b2b_main", "from_bus" => "B2", "to_bus" => "B4", "p_transfer_mw" => 90.0, "loss_mw" => 5.0)])
      @test applyConfiguredControllers!(net_yaml, cfg) == 1
      ctrl = only(Sparlectra._hvdc_pair_controllers(net_yaml))
      @test ctrl.name == "b2b_main"
      _hvdc_run!(net_yaml)
      @test net_yaml.prosumpsVec[yf_idx].pVal == net_prog.prosumpsVec[f_idx].pVal == -90.0
      @test net_yaml.prosumpsVec[yt_idx].pVal == net_prog.prosumpsVec[t_idx].pVal == 85.0
      # idempotent re-apply
      @test applyConfiguredControllers!(net_yaml, cfg) == 0
      # structural validation catches an unknown key at load time
      @test_throws ArgumentError Sparlectra._validate_controller_entries([Dict{String,Any}("type" => "hvdc_pair", "name" => "x", "from_bus" => "B2", "to_bus" => "B4", "p_transfer_mw" => 1.0, "typo_key" => 1)])
    end

    @testset "MATPOWER dcline paired_control mode" begin
      # two islands joined by one dcline row (17 columns, so the effective
      # PT is recomputed from LOSS0/LOSS1: 80 - (3 + 0.0125 * 80) = 76)
      case = mktempdir()
      path = joinpath(case, "case_b2b.m")
      write(
        path,
        """
function mpc = case_b2b
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
1 3 0 0 0 0 1 1.0 0 380 1 1.1 0.9;
2 1 40 10 0 0 1 1.0 0 380 1 1.1 0.9;
3 2 0 0 0 0 1 1.0 0 380 1 1.1 0.9;
4 1 50 12 0 0 1 1.0 0 380 1 1.1 0.9;
];
mpc.gen = [
1 60 0 300 -300 1.0 100 1 300 0;
3 40 0 300 -300 1.0 100 1 300 0;
];
mpc.branch = [
1 2 0.01 0.08 0 250 250 250 0 0 1 -360 360;
3 4 0.01 0.08 0 250 250 250 0 0 1 -360 360;
];
mpc.dcline = [
2 4 1 80 76 0 0 1.0 1.0 0 100 -50 50 -50 50 3 0.0125;
];
""",
      )
      net_paired = createNetFromMatPowerFile(filename = path, matpower_dcline_mode = :paired_control)
      ctrl = only(Sparlectra._hvdc_pair_controllers(net_paired))
      @test ctrl.name == "DCLINE_1"
      # r0.9.9: the persistent link record exists and carries the controller
      @test length(net_paired.hvdcLinks) == 1
      @test net_paired.hvdcLinks[1].controller_name == "DCLINE_1"
      @test net_paired.hvdcLinks[1].source == :matpower && net_paired.hvdcLinks[1].kind == :p2p
      @test ctrl.p_transfer_mw == 80.0
      @test ctrl.loss_mw == 3.0
      @test ctrl.loss_fraction == 0.0125
      result = _hvdc_run!(net_paired)
      @test result.converged == true
      # invariant matches the MATPOWER loss model exactly
      m = only(net_paired.matpowerDclineMetadata)
      @test net_paired.prosumpsVec[m.from_prosumer].pVal == -80.0
      @test net_paired.prosumpsVec[m.to_prosumer].pVal == 80.0 - (3.0 + 0.0125 * 80.0)
      # consistency anchor: with the target equal to the row's PF, the
      # controlled solve reproduces the pf_injections solve
      net_fixed = createNetFromMatPowerFile(filename = path, matpower_dcline_mode = :pf_injections)
      # r0.9.9: Stage-0 import records the link without a controller
      @test length(net_fixed.hvdcLinks) == 1
      @test net_fixed.hvdcLinks[1].controller_name === nothing
      _, erg = runpf!(net_fixed, 40, 1e-9, 0; method = :rectangular, islands_enabled = true)
      @test erg == 0
      vm_fixed = [something(n._vm_pu, NaN) for n in net_fixed.nodeVec]
      vm_paired = [something(n._vm_pu, NaN) for n in net_paired.nodeVec]
      @test isapprox(vm_fixed, vm_paired; atol = 1e-9)
      # a changed transfer target actually moves the flow
      ctrl.p_transfer_mw = 40.0
      ctrl.p_applied = false
      result2 = _hvdc_run!(net_paired)
      @test result2.converged == true
      @test net_paired.prosumpsVec[m.from_prosumer].pVal == -40.0
      @test net_paired.prosumpsVec[m.to_prosumer].pVal == 40.0 - (3.0 + 0.0125 * 40.0)
    end

    @testset "island_feed mode (grid-forming to side)" begin
      # mode-dependent validation rules
      net_val, _ = _build_hvdc_island_feed_net("hvdc_if_val")
      @test_throws ErrorException addHvdcPairControl!(net_val; from_bus = "A2", to_bus = "C2", mode = :island_feed, p_transfer_mw = 10.0)
      @test_throws ErrorException addHvdcPairControl!(net_val; from_bus = "A2", to_bus = "C2", mode = :island_feed, to_q_mvar = 5.0)
      @test_throws ErrorException addHvdcPairControl!(net_val; from_bus = "A2", to_bus = "C2", mode = :bogus)
      net_val2, _, _ = _build_hvdc_two_area_net("hvdc_if_val2")
      # island_feed needs the to side AS the island reference ...
      @test_throws ErrorException addHvdcPairControl!(net_val2; from_bus = "B2", to_bus = "B4", mode = :island_feed)
      # ... and setpoint mode still needs its transfer
      @test_throws ErrorException addHvdcPairControl!(net_val2; from_bus = "B2", to_bus = "B4")

      # happy path: the mirror settles on island draw + loss without touching
      # the reference injection
      net, s_idx = _build_hvdc_island_feed_net("hvdc_if_run")
      ctrl = addHvdcPairControl!(net; from_bus = "A2", to_bus = "C2", mode = :island_feed, loss_mw = 4.0)
      result = _hvdc_run!(net)
      @test result.converged
      @test ctrl.converged && !ctrl.at_limit && ctrl.status == :converged
      # island draw = 50 MW load plus a small line loss
      draw = ctrl.p_transfer_mw - 4.0
      @test 50.0 < draw < 51.0
      # the sending side mirrors draw + loss exactly
      @test isapprox(net.prosumpsVec[s_idx].pVal, -ctrl.p_transfer_mw; atol = 1e-9)
      rows = Sparlectra.control_report_rows(ctrl, net, Sparlectra.NoControlState(), (outer_iteration = ctrl.outer_iters,))
      @test rows[1].mode == :island_feed
      elems = controllableElements(net)
      row = only(filter(e -> e.actuator == :hvdc_p_transfer_mw, elems))
      @test row.device == "back-to-back HVDC pair (grid-forming)"
      # the bus table marks the converter roles in the Control column: the
      # steered sending injection as B2B, the grid-forming reference as
      # "B2B src" (its Type stays SLACK: the model is an ideal reference at
      # the PCC; SOURCE is reserved for the external-grid feeder element)
      mktempdir() do d
        cd(d) do
          printACPFlowResults(net, 0.0, 1, 1e-8, true)
          txt = read("result_$(net.name).txt", String)
          @test occursin("B2B src", txt)
        end
      end

      # rating clamp: the island draws more than the link can carry; honest
      # at_limit, sending side pinned at the rating (the island's slack still
      # balances in the model, the flag marks the violated rating)
      net_cl, s_cl = _build_hvdc_island_feed_net("hvdc_if_clamp")
      ctrl_cl = addHvdcPairControl!(net_cl; from_bus = "A2", to_bus = "C2", mode = :island_feed, loss_mw = 4.0, p_rating_mw = 40.0)
      result_cl = _hvdc_run!(net_cl)
      @test result_cl.converged == true    # blocked counts as done (TCSC convention)
      @test !ctrl_cl.converged
      @test ctrl_cl.at_limit && ctrl_cl.status == :at_limit
      @test isapprox(net_cl.prosumpsVec[s_cl].pVal, -40.0; atol = 1e-9)

      # YAML declaration with mode: island_feed instantiates and converges
      net_y, sy_idx = _build_hvdc_island_feed_net("hvdc_if_yaml")
      cfg = ControlConfig(controllers = Any[Dict{String,Any}("type" => "hvdc_pair", "name" => "if_link", "from_bus" => "A2", "to_bus" => "C2", "mode" => "island_feed", "loss_mw" => 4.0)])
      @test applyConfiguredControllers!(net_y, cfg) == 1
      ctrl_y = only(Sparlectra._hvdc_pair_controllers(net_y))
      @test ctrl_y.mode == :island_feed
      _hvdc_run!(net_y)
      @test isapprox(net_y.prosumpsVec[sy_idx].pVal, net.prosumpsVec[s_idx].pVal; atol = 1e-9)
    end

    @testset "meshed AC tie: one reference per synchronous island" begin
      # two formerly asynchronous areas tied by a new AC branch: both old
      # references survive -> actionable multi-reference error, both APIs
      net, _, _ = _build_hvdc_two_area_net("hvdc_mesh_two_refs")
      addPIModelACLine!(net = net, fromBus = "B1", toBus = "B3", r_pu = 0.02, x_pu = 0.16, b_pu = 0.0, status = 1)
      err = try
        runpf!(net, 40, 1e-9, 0; method = :rectangular, islands_enabled = true)
        nothing
      catch e
        e
      end
      @test err isa ErrorException
      @test occursin("angle references", sprint(showerror, err))
      @test occursin("B1", sprint(showerror, err)) && occursin("B3", sprint(showerror, err))

      # demote the second reference to a plain PV generator: the setpoint
      # pair keeps its invariant as a parallel PQ path next to the AC tie
      net2 = Net(name = "hvdc_mesh_pv", baseMVA = 100.0)
      for b in ("B1", "B2", "B3", "B4")
        addBus!(net = net2, busName = b, vn_kV = 380.0)
      end
      addPIModelACLine!(net = net2, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net2, fromBus = "B3", toBus = "B4", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net2, fromBus = "B1", toBus = "B3", r_pu = 0.02, x_pu = 0.16, b_pu = 0.0, status = 1)
      addProsumer!(net = net2, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net2, busName = "B3", type = "GENERATOR", p = 20.0, q = 0.0, vm_pu = 1.0, isRegulated = true)
      addProsumer!(net = net2, busName = "B2", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
      addProsumer!(net = net2, busName = "B4", type = "ENERGYCONSUMER", p = 50.0, q = 12.0)
      addProsumer!(net = net2, busName = "B2", type = "GENERATOR", p = -80.0, q = 0.0)
      f2 = length(net2.prosumpsVec)
      addProsumer!(net = net2, busName = "B4", type = "GENERATOR", p = 76.0, q = 0.0)
      t2 = length(net2.prosumpsVec)
      addHvdcPairControl!(net2; from_bus = "B2", to_bus = "B4", p_transfer_mw = 100.0, loss_mw = 4.0, from_prosumer = f2, to_prosumer = t2)
      result = _hvdc_run!(net2)
      @test result.converged
      @test net2.prosumpsVec[f2].pVal == -100.0
      @test net2.prosumpsVec[t2].pVal == 96.0

      # island_feed with an AC tie: the grid-forming reference makes the
      # synchronous island carry two references -> same early error
      net3, _ = _build_hvdc_island_feed_net("hvdc_mesh_if")
      addPIModelACLine!(net = net3, fromBus = "A1", toBus = "C1", r_pu = 0.02, x_pu = 0.16, b_pu = 0.0, status = 1)
      addHvdcPairControl!(net3; from_bus = "A2", to_bus = "C2", mode = :island_feed, loss_mw = 4.0)
      err3 = try
        _hvdc_run!(net3)
        nothing
      catch e
        e
      end
      @test err3 isa ErrorException
      @test occursin("angle references", sprint(showerror, err3))

      # demoting the grid-forming reference afterwards invalidates the
      # island_feed semantics: honest invalid_topology, no silent mode switch
      net4, s4 = _build_hvdc_island_feed_net("hvdc_mesh_if_demoted")
      addPIModelACLine!(net = net4, fromBus = "A1", toBus = "C1", r_pu = 0.02, x_pu = 0.16, b_pu = 0.0, status = 1)
      ctrl4 = addHvdcPairControl!(net4; from_bus = "A2", to_bus = "C2", mode = :island_feed, loss_mw = 4.0)
      # demote: the reference injection at C2 loses its reference role
      for ps in net4.prosumpsVec
        if Sparlectra.getPosumerBusIndex(ps) == net4.busDict["C2"] && ps.referencePri !== nothing
          ps.referencePri = nothing
        end
      end
      refreshBusTypesFromProsumers!(net4)
      result4 = _hvdc_run!(net4)
      @test ctrl4.status == :invalid_topology
      @test !ctrl4.converged
      # the mirror never fired: the sending injection kept its build value
      @test something(net4.prosumpsVec[s4].pVal, 0.0) == 0.0
    end

    @testset "result table alignment across all modes and statuses (Phase 6)" begin
      # two links, one per mode, and the LONGEST status value forced on the
      # island_feed controller; the regression guard is structural: in both
      # tables every row must carry the same number of column separators
      # and the same line length as its header, for every mode/status combo
      net = Net(name = "hvdc_align", baseMVA = 100.0)
      for b in ("B1", "B2", "B3", "B4", "A1", "A2", "C1", "C2")
        addBus!(net = net, busName = b, vn_kV = 380.0)
      end
      addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "C2", toBus = "C1", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net, busName = "B3", type = "EXTERNALNETWORKINJECTION", referencePri = "B3", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", referencePri = "A1", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net, busName = "C2", type = "EXTERNALNETWORKINJECTION", referencePri = "C2", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
      addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 50.0, q = 12.0)
      addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
      addProsumer!(net = net, busName = "C1", type = "ENERGYCONSUMER", p = 50.0, q = 12.0)
      addProsumer!(net = net, busName = "B2", type = "GENERATOR", p = -80.0, q = 0.0)
      addProsumer!(net = net, busName = "B4", type = "GENERATOR", p = 76.0, q = 0.0)
      addProsumer!(net = net, busName = "A2", type = "GENERATOR", p = 0.0, q = 0.0)
      ctrl1 = addHvdcPairControl!(net; from_bus = "B2", to_bus = "B4", p_transfer_mw = 100.0, loss_mw = 4.0, p_rating_mw = 150.0)
      ctrl2 = addHvdcPairControl!(net; from_bus = "A2", to_bus = "C2", mode = :island_feed, loss_mw = 4.0)
      result = _hvdc_run!(net)
      @test result.converged
      calcNetLosses!(net)
      # cycle through EVERY status value defined in hvdc_pair_control.jl,
      # including the 16-character :invalid_topology
      for forced in Sparlectra.HVDC_PAIR_STATUS_VALUES
        ctrl2.status = forced
        txt = mktempdir() do d
          cd(d) do
            printACPFlowResults(net, 0.0, 1, 1e-8, true)
            read("result_$(net.name).txt", String)
          end
        end
        lines = split(txt, "\n")
        # branch table: from its header to the separator after the rows
        bh = findfirst(l -> startswith(l, "| Branch"), lines)
        @test bh !== nothing
        width = length(lines[bh])
        pipes = count(==('|'), lines[bh])
        i = bh + 2
        while i <= length(lines) && startswith(lines[i], "|")
          @test length(lines[i]) == width
          @test count(==('|'), lines[i]) == pipes
          i += 1
        end
        # HVDC Link Flows table
        hh = findfirst(l -> startswith(l, "| Nr"), lines)
        @test hh !== nothing
        hwidth = length(lines[hh])
        hpipes = count(==('|'), lines[hh])
        j = hh + 2
        while j <= length(lines) && startswith(lines[j], "|")
          @test length(lines[j]) == hwidth
          @test count(==('|'), lines[j]) == hpipes
          j += 1
        end
      end
      ctrl2.status = :converged
      # P_tgt in the Link row is the controller setpoint, not the flow
      txt2 = mktempdir() do d
        cd(d) do
          printACPFlowResults(net, 0.0, 1, 1e-8, true)
          read("result_$(net.name).txt", String)
        end
      end
      link_row = only(filter(l -> occursin("B2B_B2_B4", l) && occursin("Link", l), split(txt2, "\n")))
      @test occursin("100.000", link_row)
    end

    @testset "element rows and summary output" begin
      net, _, _ = _build_hvdc_two_area_net("hvdc_rows")
      addHvdcPairControl!(net; from_bus = "B2", to_bus = "B4", p_transfer_mw = 80.0, loss_mw = 3.0, p_rating_mw = 150.0)
      _hvdc_run!(net)
      elems = controllableElements(net)
      row = only(filter(e -> e.actuator == :hvdc_p_transfer_mw, elems))
      @test row.device == "back-to-back HVDC pair"
      @test row.quantity == :hvdc_transfer
      @test row.actuator_min == -150.0 && row.actuator_max == 150.0
      buf = IOBuffer()
      printHvdcPairControllerSummary(buf, net)
      out = String(take!(buf))
      @test occursin("HVDC Pair Control Summary", out)
      @test occursin("link B2 -> B4", out)
      @test occursin("77.000 MW", out)
      # r0.9.9 result surface: link table, header count, report rows, and
      # the fixed-mode row of a hand-registered Stage-0 link
      report = buildACPFlowReport(net; ite = 1, converged = true)
      @test length(report.hvdc_links) == 1
      row = only(report.hvdc_links)
      @test row.mode == :setpoint && row.p_from_MW == 80.0 && row.loss_MW == 3.0
      net_s0, _, _ = _build_hvdc_two_area_net("hvdc_rows_stage0")
      link = addHvdcLink!(net_s0; from_bus = "B2", to_bus = "B4")
      @test link.source == :api && link.controller_name === nothing
      row0 = only(Sparlectra._hvdc_link_flow_rows(net_s0))
      @test row0.mode == :fixed
      @test row0.p_from_MW == 80.0 && row0.p_to_MW == 76.0 && row0.loss_MW == 4.0
      mktempdir() do d
        cd(d) do
          printACPFlowResults(net, 0.0, 1, 1e-8, true)
          txt = read("result_$(net.name).txt", String)
          @test occursin("HVDC links     :         1", txt)
          @test occursin("HVDC Link Flows", txt)
          @test occursin("converter losses (HVDC)", txt)
          # the link also appears as a marked row in the branch table
          @test occursin("| Link   |", txt)
          @test occursin("HVDC, not a branch", txt)
        end
      end
      # regression (multi-slack header): a two-island net reports both
      # references in the node-count line instead of a hardcoded 1
      mktempdir() do d
        cd(d) do
          printACPFlowResults(net, 0.0, 1, 1e-8, true)
          txt = read("result_$(net.name).txt", String)
          @test occursin("Slack: 2", txt)
        end
      end
      # regression (island writeback): the solver-written slack outcomes
      # reach the parent net, so the bus table shows Pg/Qg per reference
      for n in net.nodeVec
        if Sparlectra.getNodeType(n) == Sparlectra.Slack
          @test n._pƩGen !== nothing
          @test n._qƩGen !== nothing
        end
      end
    end
  end
  return true
end

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

# file: test/test_legacy_hvdc_pair_control.jl
# purpose: legacy profile: the testsets of test/test_hvdc_pair_control.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_hvdc_pair_control.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_hvdc_pair_control_tests()
  # cut from testset: HVDC pair controller (#297 Draft B)
    @testset "rating clamp is honest at_limit" begin (function ()
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
    end)() end

    @testset "unreachable voltage target clamps honestly" begin (function ()
      net, _, _ = _build_hvdc_two_area_net("hvdc_vset_limit")
      addHvdcPairControl!(net; from_bus = "B2", to_bus = "B4", p_transfer_mw = 80.0, loss_mw = 3.0, to_vset_pu = 1.2, to_qmin_mvar = -5.0, to_qmax_mvar = 5.0)
      _hvdc_run!(net)
      ctrl = only(Sparlectra._hvdc_pair_controllers(net))
      @test !ctrl.converged
      @test ctrl.at_limit
      @test isapprox(ctrl.to_q_now, 5.0; atol = 1e-9)
    end)() end

    @testset "disabled controller keeps the baseline bit-identical" begin (function ()
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
    end)() end

    @testset "YAML declaration round trip (#305 schema)" begin (function ()
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
    end)() end

    @testset "MATPOWER dcline paired_control mode" begin (function ()
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
      # #465: MATPOWER carries one reference per synchronous area. With bus 3
      # typed 3 as well, each island keeps the reference of the file (before
      # 0.30.2 only the first type-3 bus became Slack and bus 3 lost its
      # voltage control); joined by an AC branch, the second reference of the
      # same island stays a voltage-regulated PV bus, with a warning
      two_refs = replace(read(path, String), "3 2 0 0 0 0 1 1.0 0 380 1 1.1 0.9;" => "3 3 0 0 0 0 1 1.0 0 380 1 1.1 0.9;")
      path_refs = joinpath(case, "case_b2b_two_refs.m")
      write(path_refs, two_refs)
      net_refs = createNetFromMatPowerFile(filename = path_refs, matpower_dcline_mode = :pf_injections)
      @test Sparlectra.getNodeType(net_refs.nodeVec[1]) == Sparlectra.Slack
      @test Sparlectra.getNodeType(net_refs.nodeVec[3]) == Sparlectra.Slack
      _, erg_refs = runpf!(net_refs, 40, 1e-9, 0; method = :rectangular, islands_enabled = true)
      @test erg_refs == 0
      # the SCF conversion of the case carries the same two references
      mpc_refs = Sparlectra.MatpowerIO.read_case(path_refs; legacy_compat = true)
      scf_refs = Sparlectra.convert_case(Sparlectra.MatpowerAdapter(), mpc_refs, Sparlectra.matpower_adapter_options(Sparlectra.active_sparlectra_config()))
      @test length(scf_refs.sparlectra.roles.slack.nodes) == 2
      path_tied = joinpath(case, "case_b2b_tied.m")
      write(path_tied, replace(two_refs, "3 4 0.01 0.08 0 250 250 250 0 0 1 -360 360;" => "3 4 0.01 0.08 0 250 250 250 0 0 1 -360 360;\n2 3 0.01 0.08 0 250 250 250 0 0 1 -360 360;"))
      net_tied = @test_logs (:warn, r"second reference bus 3 in the AC island of reference bus 1") match_mode = :any createNetFromMatPowerFile(filename = path_tied, matpower_dcline_mode = :pf_injections)
      @test Sparlectra.getNodeType(net_tied.nodeVec[1]) == Sparlectra.Slack
      @test Sparlectra.getNodeType(net_tied.nodeVec[3]) == Sparlectra.PV
      # type 3 is reference priority 1 for every unit of the bus, the demoted
      # second reference included: it is the first to take over its island
      # (0.30.2); the SCF conversion writes the same priorities
      @test [ps.referencePriority for ps in net_tied.prosumpsVec if Sparlectra.getPosumerBusIndex(ps) in (1, 3)] == [1, 1]
      @test count(e -> e isa AbstractDict && get(e, "reference_priority", 0) == 1, values(scf_refs.sparlectra.extra)) == 2
    end)() end

    @testset "result table alignment across all modes and statuses (Phase 6)" begin (function ()
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
    end)() end

    @testset "element rows and summary output" begin (function ()
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
    end)() end

  return true
end

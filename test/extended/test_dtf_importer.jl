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

const DTF_FIXTURE = joinpath(@__DIR__, "..", "fixtures", "dtf", "FOR001.DAT")

function _multiplicity(branches, from::String, to::String)
  return count(b -> b.from == from && b.to == to, branches)
end

function _branch(case, from::String, to::String; kind::Char = 'L')
  matches = [b for b in case.branches if b.from == from && b.to == to && b.kind == kind]
  return first(matches)
end

function _branch(case, from::String, to::String, parallel_id::String; kind::Char = 'L')
  matches = [b for b in case.branches if b.from == from && b.to == to && b.kind == kind && b.parallel_id == parallel_id]
  return first(matches)
end

function _control(case, from::String, to::String, parallel_id::String)
  matches = [c for c in case.transformer_controls if c.from == from && c.to == to && c.parallel_id == parallel_id]
  return first(matches)
end

function _net_multiplicity(net, from::String, to::String)
  from_idx = Sparlectra.geNetBusIdx(net = net, busName = from)
  to_idx = Sparlectra.geNetBusIdx(net = net, busName = to)
  return count(br -> br.fromBus == from_idx && br.toBus == to_idx, net.branchVec)
end

function run_dtf_importer_tests()
  @testset "native DTF/FOR001 importer MVP" begin
    case = Sparlectra.DTFImporter.read_dtf(DTF_FIXTURE)

    @test case.nominal_voltages_kv == [231.0, 230.0, 400.0]
    @test case.size.NGES == 13
    @test case.size.LGES == 27
    @test case.size.IKOMP == 0
    @test case.size.IRETRA == 5
    @test case.size.slack == "BSTADTS1"

    @test length(case.branches) == 27
    @test count(b -> b.kind == 'L', case.branches) == 22
    @test count(b -> b.kind == 'T', case.branches) == 5

    expected_parallel = Dict(
      ("BETA2 S1", "BURG  S1") => 3,
      ("BETA1 S1", "BETA2 S1") => 3,
      ("DELTA1S1", "DELTA2S1") => 2,
      ("BETA1 S1", "DELTA1S1") => 2,
      ("BSTADTS1", "SUED  S1") => 2,
    )
    for ((from, to), expected) in expected_parallel
      @test _multiplicity(case.branches, from, to) == expected
    end

    line = _branch(case, "ALPHA S1", "BETA1 S1"; kind = 'L')
    line_pu = Sparlectra.DTFImporter._branch_pu(case, line)
    zbase_231 = 231.0^2 / 100.0
    @test line.r_ohm == 0.603
    @test line.x_ohm == 5.3
    @test line.b_s == 0.000805
    @test line.voltage_level_index == 1
    @test line_pu.u_ref_kv == 231.0
    @test line_pu.r ≈ 0.603 / zbase_231
    @test line_pu.x ≈ 5.3 / zbase_231
    @test line_pu.b ≈ 0.000805 * zbase_231
    @test !(line_pu.b ≈ 0.000805 * (400.0^2 / 100.0))

    trafo = _branch(case, "BETA1 S1", "BETA2 S1"; kind = 'T')
    trafo_pu = Sparlectra.DTFImporter._branch_pu(case, trafo)
    @test trafo.x_ohm == 8.22
    @test trafo.b_s == -1.84e-5
    @test trafo.voltage_level_index == 1
    @test trafo_pu.u_ref_kv == 231.0
    @test trafo_pu.x ≈ 8.22 / zbase_231
    @test trafo_pu.b ≈ -1.84e-5 * zbase_231
    @test !(trafo_pu.x ≈ 0.0051375)

    beta_control = _control(case, "BETA1 S1", "BETA2 S1", "C")
    @test beta_control.regulated_side == "R"
    @test beta_control.unregulated_side == ""
    @test beta_control.parallel_id == "C"
    @test beta_control.phase_shifter_flag == ""
    @test beta_control.nominal_unregulated_kv == 400.0
    @test beta_control.nominal_regulated_kv == 231.0
    @test beta_control.longitudinal_range_percent == 12.5
    @test beta_control.added_voltage_angle_deg == 0.0
    @test beta_control.max_tap_step == 9
    @test beta_control.actual_tap_step == 0
    @test beta_control.quadrature_range_percent == 0.0
    @test beta_control.quadrature_max_steps == 0
    @test beta_control.quadrature_actual_step == 0

    delta_control = _control(case, "DELTA1S1", "DELTA2S1", "B")
    @test delta_control.nominal_unregulated_kv == 400.0
    @test delta_control.nominal_regulated_kv == 231.0
    @test delta_control.longitudinal_range_percent == 18.0
    @test delta_control.max_tap_step == 13
    @test delta_control.actual_tap_step == 0

    case_e = Sparlectra.DTFImporter.read_dtf(joinpath(@__DIR__, "..", "..", "data", "DTF", "FOR001E.DAT"))
    delta_e_control = _control(case_e, "DELTA1S1", "DELTA2S1", "B")
    delta_e_branch = _branch(case_e, "DELTA1S1", "DELTA2S1", "B"; kind = 'T')
    bus_by_name_e = Dict(b.name => b for b in case_e.buses)
    delta_e_tap = Sparlectra.DTFImporter._dtf_effective_transformer_tap(
      case_e,
      delta_e_branch,
      delta_e_control,
      bus_by_name_e[delta_e_branch.from],
      bus_by_name_e[delta_e_branch.to],
    )
    expected_fraction = 0.18 * 7 / 13
    expected_complex = 1 + expected_fraction * cis(deg2rad(60.0))
    @test delta_e_tap.model == :skew_angle
    @test delta_e_tap.tap_fraction ≈ expected_fraction
    @test delta_e_tap.skew_angle_deg == 60.0
    @test delta_e_tap.effective_complex ≈ expected_complex
    @test delta_e_tap.ratio ≈ (230.0 / 231.0) / abs(expected_complex)
    @test delta_e_tap.shift_deg ≈ -rad2deg(angle(expected_complex))
    @test delta_e_tap.convention == :dtf_regulating_vector_reciprocal

    net = Sparlectra.createNetFromDTFFile(DTF_FIXTURE)
    @test net isa Sparlectra.Net
    ok, msg = Sparlectra.validate!(net = net)
    @test ok
    @test length(net.nodeVec) == 13
    @test length(net.branchVec) == 27
    for ((from, to), expected) in expected_parallel
      @test _net_multiplicity(net, from, to) == expected
    end

    beta_branch = _branch(case, "BETA1 S1", "BETA2 S1", "C"; kind = 'T')
    delta_branch = _branch(case, "DELTA1S1", "DELTA2S1", "B"; kind = 'T')
    @test net.branchVec[beta_branch.index].ratio ≈ 230.0 / 231.0
    @test !(net.branchVec[beta_branch.index].ratio ≈ 400.0 / 231.0)
    @test net.branchVec[delta_branch.index].ratio ≈ 230.0 / 231.0

    for bus in case.buses
      idx = Sparlectra.geNetBusIdx(net = net, busName = bus.name)
      @test Sparlectra.getNodeVn(net.nodeVec[idx]) == case.nominal_voltages_kv[bus.voltage_level_index]
      @test net.nodeVec[idx]._vm_pu == 1.0
    end

    for bus_name in ("BETA2 S1", "NORD  S1", "WEILERS1", "SUED  S1")
      idx = Sparlectra.geNetBusIdx(net = net, busName = bus_name)
      gens = [ps for ps in net.prosumpsVec if Sparlectra.getPosumerBusIndex(ps) == idx && Sparlectra.isGenerator(ps.proSumptionType)]
      @test !isempty(gens)
      @test all(ps -> !ps.isRegulated, gens)
      @test Sparlectra.isPQNode(net.nodeVec[idx])
    end

    slack_idx = Sparlectra.geNetBusIdx(net = net, busName = "BSTADTS1")
    @test slack_idx in net.slackVec
    slack_gens = [ps for ps in net.prosumpsVec if Sparlectra.getPosumerBusIndex(ps) == slack_idx && ps.referencePri == slack_idx]
    @test !isempty(slack_gens)

    synthetic = Sparlectra.DTFImporter.DTFCase(
      "synthetic",
      100.0,
      Sparlectra.DTFImporter.DTFParams("", Float64[]),
      ["synthetic"],
      [230.0],
      Sparlectra.DTFImporter.DTFSize("", 2, 1, 0, 0, "SLACK"),
      [Sparlectra.DTFImporter.DTFBranch("", 1, 'L', 1, "", "PV", "SLACK", 0.01, 0.1, 0.0, 0.0, nothing)],
      Sparlectra.DTFImporter.DTFCompensation[],
      Sparlectra.DTFImporter.DTFTransformerControl[],
      [
        Sparlectra.DTFImporter.DTFBus("", 1, 1, 1, "PV", 230.0, 0.0, 0.0, 0.0, 10.0, 2.0, -5.0, 5.0),
        Sparlectra.DTFImporter.DTFBus("", 2, 2, 1, "SLACK", 230.0, 0.0, 0.0, 0.0, 20.0, 3.0, -10.0, 10.0),
      ],
      Sparlectra.DTFImporter.DTFOutage[],
      Sparlectra.DTFImporter.DTFTrailingRecord[],
    )
    synthetic_net = Sparlectra.DTFImporter.build_net(synthetic)
    pv_idx = Sparlectra.geNetBusIdx(net = synthetic_net, busName = "PV")
    pv_gens = [ps for ps in synthetic_net.prosumpsVec if Sparlectra.getPosumerBusIndex(ps) == pv_idx && Sparlectra.isGenerator(ps.proSumptionType)]
    @test !isempty(pv_gens)
    @test any(ps -> ps.isRegulated && ps.vm_pu == 1.0, pv_gens)

    @test case.outages isa Vector{Sparlectra.DTFImporter.DTFOutage}
    @test length(case.outages) == 2
    @test case.outages[1].from == "ALPHA S1"
    @test case.outages[1].to == "BETA1 S1"
    @test case.outages[2].parallel_id == "A"
    @test case.outages[2].from == "BETA1 S1"
    @test case.outages[2].to == "BETA2 S1"
  end
end

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

function _synthetic_dtf_case(; kind::Char = 'T', g_s::Float64 = 4.0e-5, b_s::Float64 = -1.0e-5, path::String = "synthetic")
  return Sparlectra.DTFImporter.DTFCase(
    path,
    100.0,
    Sparlectra.DTFImporter.DTFParams("", Float64[]),
    ["synthetic"],
    [110.0],
    Sparlectra.DTFImporter.DTFSize("", 2, 1, 0, 0, "SLACK"),
    [Sparlectra.DTFImporter.DTFBranch("", 1, kind, 1, "A", "PV", "SLACK", 1.21, 6.05, g_s, b_s, nothing)],
    Sparlectra.DTFImporter.DTFCompensation[],
    Sparlectra.DTFImporter.DTFTransformerControl[],
    [
      Sparlectra.DTFImporter.DTFBus("", 1, 1, 1, "PV", 110.0, 0.0, 0.0, 0.0, 10.0, 2.0, -5.0, 5.0),
      Sparlectra.DTFImporter.DTFBus("", 2, 2, 1, "SLACK", 110.0, 0.0, 0.0, 0.0, 20.0, 3.0, -10.0, 10.0),
    ],
    [Sparlectra.DTFImporter.DTFOutage("", 1, 'T', 1, "A", "PV", "SLACK")],
    [Sparlectra.DTFImporter.DTFTrailingRecord("trailing branch echo", 1, :branch, :echo, nothing)],
  )
end

function _synthetic_dtf_case_with_tap_control(; longitudinal_range_percent::Float64, actual_tap_step::Int, max_tap_step::Int, added_voltage_angle_deg::Float64 = 0.0)
  control = Sparlectra.DTFImporter.DTFTransformerControl(
    "", 1, "PV", "SLACK", "A", "", "PV", "SLACK",
    110.0, 110.0, longitudinal_range_percent, added_voltage_angle_deg, max_tap_step, actual_tap_step,
    nothing, nothing, nothing,
  )
  return Sparlectra.DTFImporter.DTFCase(
    "synthetic_tap",
    100.0,
    Sparlectra.DTFImporter.DTFParams("", Float64[]),
    ["synthetic_tap"],
    [110.0],
    Sparlectra.DTFImporter.DTFSize("", 2, 1, 0, 1, "SLACK"),
    [Sparlectra.DTFImporter.DTFBranch("", 1, 'T', 1, "A", "PV", "SLACK", 1.21, 6.05, 0.0, 0.0, nothing)],
    Sparlectra.DTFImporter.DTFCompensation[],
    [control],
    [
      Sparlectra.DTFImporter.DTFBus("", 1, 1, 1, "PV", 110.0, 0.0, 0.0, 0.0, 10.0, 2.0, -5.0, 5.0),
      Sparlectra.DTFImporter.DTFBus("", 2, 2, 1, "SLACK", 110.0, 0.0, 0.0, 0.0, 20.0, 3.0, -10.0, 10.0),
    ],
    Sparlectra.DTFImporter.DTFOutage[],
    Sparlectra.DTFImporter.DTFTrailingRecord[],
  )
end

function run_dtf_importer_tests()
  @testset "native DTF importer focused synthetic checks" begin
    case = _synthetic_dtf_case()
    branch = only(case.branches)
    pu = Sparlectra.DTFImporter._branch_pu(case, branch)
    zbase = 110.0^2 / 100.0
    @test pu.r ≈ 1.21 / zbase
    @test pu.x ≈ 6.05 / zbase
    @test pu.g ≈ 4.0e-5 * zbase
    @test pu.b ≈ -1.0e-5 * zbase

    net = Sparlectra.DTFImporter.build_net(case)
    @test net isa Sparlectra.Net
    @test length(net.trafos) == 1
    @test length(net.branchVec) == 1
    winding = only(net.trafos).side1
    @test winding.g ≈ pu.g
    @test winding.b ≈ pu.b
    @test net.branchVec[1].g_pu ≈ pu.g
    @test net.branchVec[1].b_pu ≈ pu.b
    @test Sparlectra.getTrafoRXBG(winding) == (winding.r, winding.x, winding.b, winding.g)
    r_pu, x_pu, b_pu, g_pu = Sparlectra.getTrafoRXBG_pu(winding, 110.0, 100.0)
    @test r_pu ≈ pu.r
    @test x_pu ≈ pu.x
    @test b_pu ≈ pu.b
    @test g_pu ≈ pu.g
    @test Sparlectra.bus_shunt_totals_pu(net).total_g_pu ≈ 0.0
    @test get(net.matpower_branch_metadata, 1, nothing).transformer_loss_allocation == :native_branch_pi

    zero_g_net = Sparlectra.DTFImporter.build_net(_synthetic_dtf_case(g_s = 0.0))
    @test zero_g_net.branchVec[1].g_pu == 0.0
    @test Sparlectra.bus_shunt_totals_pu(zero_g_net).total_g_pu ≈ 0.0

    @test length(case.outages) == 1
    @test only(case.outages).from == "PV"
    @test length(case.trailing_records) == 1

    # The classic result print shows `net.name` as "Case"; it must reflect the
    # originating file name, not a free-text description card from the file.
    named_net = Sparlectra.DTFImporter.build_net(_synthetic_dtf_case(path = "/some/dir/FOR001B.DAT"))
    @test named_net.name == "FOR001B.DAT"
    @test net.name == "synthetic"
  end

  @testset "native DTF importer honors configured tap-changer model" begin
    case = _synthetic_dtf_case_with_tap_control(longitudinal_range_percent = 10.0, actual_tap_step = 7, max_tap_step = 10)

    net_ideal = Sparlectra.DTFImporter.build_net(case; tap_changer_model = :ideal)
    net_corrected = Sparlectra.DTFImporter.build_net(case; tap_changer_model = :impedance_correction)

    branch_ideal = only(net_ideal.branchVec)
    branch_corrected = only(net_corrected.branchVec)

    tap_fraction = (10.0 / 100.0) * 7 / 10
    expected_factor = (1.0 + tap_fraction)^2
    @test isapprox(branch_corrected.r_pu, branch_ideal.r_pu * expected_factor; atol = 1e-12)
    @test isapprox(branch_corrected.x_pu, branch_ideal.x_pu * expected_factor; atol = 1e-12)
    @test get(net_corrected.matpower_branch_metadata, 1, nothing).tap_changer_model == :impedance_correction
    @test isapprox(get(net_corrected.matpower_branch_metadata, 1, nothing).tap_impedance_correction_factor, expected_factor; atol = 1e-12)
    @test get(net_ideal.matpower_branch_metadata, 1, nothing).tap_impedance_correction_factor == 1.0
  end

  @testset "native DTF full local fixture" begin
    if !isfile(DTF_FIXTURE)
      @info "Skipping full FOR001 fixture validation; external reference network is not tracked. Place a local file at test/fixtures/dtf/FOR001.DAT or data/DTF/FOR001.DAT for manual validation."
      @test_skip "FOR001 full reference fixture not available"
      return
    end
    case = Sparlectra.DTFImporter.read_dtf(DTF_FIXTURE)
    @test case.size.NGES == 13
    @test length(case.branches) == 27
    @test count(b -> b.kind == 'T', case.branches) == 5
    net = Sparlectra.createNetFromDTFFile(DTF_FIXTURE)
    @test length(net.branchVec) == 27
    @test Sparlectra.bus_shunt_totals_pu(net).total_g_pu ≈ 0.0
  end
end

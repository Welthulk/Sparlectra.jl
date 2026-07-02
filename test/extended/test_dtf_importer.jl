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

    net = Sparlectra.createNetFromDTFFile(DTF_FIXTURE)
    @test net isa Sparlectra.Net
    ok, msg = Sparlectra.validate!(net = net)
    @test ok
    @test length(net.nodeVec) == 13
    @test length(net.branchVec) == 27
    for ((from, to), expected) in expected_parallel
      @test _net_multiplicity(net, from, to) == expected
    end

    @test case.outages isa Vector{Sparlectra.DTFImporter.DTFOutage}
    @test length(case.outages) == 2
    @test case.outages[1].from == "ALPHA S1"
    @test case.outages[1].to == "BETA1 S1"
    @test case.outages[2].parallel_id == "A"
    @test case.outages[2].from == "BETA1 S1"
    @test case.outages[2].to == "BETA2 S1"
  end
end

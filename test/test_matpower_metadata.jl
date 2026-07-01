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

function _metadata_case(; branch_kind = ["L"], dcline = nothing)
  return Sparlectra.MatpowerIO.MatpowerCase(
    "metadata_case",
    100.0,
    [2 1 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9;
     1 3 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9],
    [1 10.0 0.0 100.0 -100.0 1.0 100.0 1 100.0 0.0 0 0 0 0 0 0 0 0 0 0 0],
    [1 2 0.01 0.05 0.0 100.0 0.0 0.0 0.0 0.0 1 -60.0 60.0],
    nothing,
    ["LoadBus", "SlackBus"],
    ["L1 SlackBus -> LoadBus"],
    branch_kind,
    ["L1 SlackBus -> LoadBus"],
    dcline,
  )
end

function run_matpower_metadata_tests()
  @testset "MATPOWER metadata import" begin
    mktempdir() do dir
      path = joinpath(dir, "case_meta.m")
      write(path, """
function mpc = case_meta
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
2 1 0 0 0 0 1 1 0 110 1 1.1 0.9;
1 3 0 0 0 0 1 1 0 110 1 1.1 0.9;
];
mpc.gen = [
1 10 0 100 -100 1 100 1 100 0 0 0 0 0 0 0 0 0 0 0 0;
];
mpc.branch = [
1 2 0.01 0.05 0 100 0 0 0 0 1 -60 60;
];
mpc.bus_name = {
'LoadBus';
'SlackBus';
};
mpc.branch_name = {
'L1 SlackBus -> LoadBus';
};
mpc.branch_kind = {
'T';
};
mpc.for001_contingencies = {
'L1 SlackBus -> LoadBus';
};
mpc.dcline = [
1 2 1 100 0 0 0 1 1 0 200 -50 50 -50 50 2 0.01;
2 1 0 50 49 0 0 1 1 0 200;
];
""")
      mpc_unsorted = Sparlectra.MatpowerIO.read_case_m(path; legacy_compat = false)
      @test mpc_unsorted.bus_name == ["LoadBus", "SlackBus"]
      @test mpc_unsorted.branch_name == ["L1 SlackBus -> LoadBus"]
      @test mpc_unsorted.branch_kind == ["T"]
      @test mpc_unsorted.for001_contingencies == ["L1 SlackBus -> LoadBus"]
      @test size(mpc_unsorted.dcline) == (2, 17)
      mpc_sorted = Sparlectra.MatpowerIO.legacy_sort_bus(mpc_unsorted)
      @test mpc_sorted.bus[:, 1] == [1.0, 2.0]
      @test mpc_sorted.bus_name == ["SlackBus", "LoadBus"]
      @test Sparlectra.MatpowerIO.for001_contingency_branch_indices(mpc_sorted) == [1]
    end

    mpc = Sparlectra.MatpowerIO.legacy_sort_bus(_metadata_case())
    net_default = Sparlectra.createNetFromMatPowerCase(mpc = mpc)
    @test haskey(net_default.busDict, "1")
    @test haskey(net_default.busDict, "2")

    net_named = Sparlectra.createNetFromMatPowerCase(mpc = mpc, apply_bus_names = true, apply_branch_names = true)
    @test haskey(net_named.busDict, "SlackBus")
    @test haskey(net_named.busDict, "LoadBus")
    @test net_named.busOrigIdxDict[net_named.busDict["SlackBus"]] == 1
    @test first(values(net_named.matpower_branch_metadata)).orig_name == "L1 SlackBus -> LoadBus"
    @test net_named.for001Contingencies == ["L1 SlackBus -> LoadBus"]

    net_heuristic = Sparlectra.createNetFromMatPowerCase(mpc = mpc, apply_bus_names = true, apply_branch_kind = false)
    @test length(net_heuristic.linesAC) == 1
    @test length(net_heuristic.trafos) == 0
    net_kind = Sparlectra.createNetFromMatPowerCase(mpc = _metadata_case(branch_kind = ["T"]), apply_bus_names = true, apply_branch_kind = true)
    @test length(net_kind.linesAC) == 0
    @test length(net_kind.trafos) == 1

    dcline = [1.0 2.0 1.0 100.0 0.0 4.0 5.0 1.0 1.0 0.0 200.0 -50.0 50.0 -50.0 50.0 2.0 0.01]
    mpc_dcline = Sparlectra.MatpowerIO.legacy_sort_bus(_metadata_case(dcline = dcline))
    @test_throws Sparlectra.MatpowerIO.UnsupportedMatpowerDclineError Sparlectra.createNetFromMatPowerCase(mpc = mpc_dcline)
    net_dcline = Sparlectra.createNetFromMatPowerCase(mpc = mpc_dcline, apply_bus_names = true, matpower_dcline_mode = :pf_injections)
    meta = only(net_dcline.matpowerDclineMetadata)
    @test meta.pf_mw == 100.0
    @test meta.pt_mw == 97.0
    @test meta.loss_mw == 3.0
    @test net_dcline.prosumpsVec[meta.from_prosumer].pVal == -100.0
    @test net_dcline.prosumpsVec[meta.to_prosumer].pVal == 97.0

    inactive = [1.0 2.0 0.0 100.0 99.0]
    net_inactive = Sparlectra.createNetFromMatPowerCase(mpc = Sparlectra.MatpowerIO.legacy_sort_bus(_metadata_case(dcline = inactive)), matpower_dcline_mode = :reject_active)
    @test isempty(net_inactive.matpowerDclineMetadata)
  end
end

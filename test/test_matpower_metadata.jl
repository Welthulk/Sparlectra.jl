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
    # Regression: a local loop variable used to build bus_original_name_by_orig
    # from mpc.bus_name shared its name with the outer case-name variable, so
    # net.name ended up as the last bus's original name (e.g. "LoadBus")
    # instead of the case name, whenever mpc.bus_name metadata was present
    # (independent of apply_bus_names).
    @test net_default.name == "metadata_case"

    net_named = Sparlectra.createNetFromMatPowerCase(mpc = mpc, apply_bus_names = true, apply_branch_names = true)
    @test net_named.name == "metadata_case"
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

    @testset "Sparlectra transformer-loss extension" begin
      net_loss = Sparlectra.Net(name = "loss_extension", baseMVA = 100.0)
      Sparlectra.addBus!(net = net_loss, busName = "FROM", vn_kV = 110.0, oBusIdx = 1)
      Sparlectra.addBus!(net = net_loss, busName = "TO", vn_kV = 110.0, oBusIdx = 2)
      Sparlectra.addProsumer!(net = net_loss, busName = "FROM", type = "GENERATOR", p = 0.0, q = 0.0, referencePri = "FROM", vm_pu = 1.0)
      Sparlectra._addPIModelTrafo_by_idx!(net = net_loss, from = 1, to = 2, r_pu = 0.01, x_pu = 0.05, b_pu = -0.002, g_pu = 0.005, status = 1, ratedS = 100.0, ratio = 1.0, shift_deg = 0.0)
      net_loss.matpower_branch_metadata[1] = (
        orig_name = "T1A FROM -> TO",
        source_label = "T1A",
        orig_kind = :transformer,
        orig_index = 1,
        dtf_kind = 'T',
        u_ref_kV = 110.0,
        from_bus_vn_kV = 110.0,
        to_bus_vn_kV = 110.0,
        r_ohm = 1.21,
        x_ohm = 6.05,
        g_s = 4.132231404958678e-5,
        b_s = -1.652892561983471e-5,
        g_pu = 0.005,
        b_pu = -0.002,
        active_no_load_g_pu = 0.005,
        transformer_loss_allocation = :native_branch_pi,
        tap_ratio = 1.0,
        actual_tap_step = 0,
        max_tap_step = 9,
        phase_shift_deg = 0.0,
      )
      @test net_loss.branchVec[1].g_pu ≈ 0.005
      @test only(net_loss.trafos).side1.g ≈ 0.005
      @test Sparlectra.getTrafoRXBG(only(net_loss.trafos).side1)[4] ≈ 0.005
      @test Sparlectra.getTrafoRXBG_pu(only(net_loss.trafos).side1, 110.0, 100.0)[4] ≈ 0.005
      @test Sparlectra.bus_shunt_totals_pu(net_loss).total_g_pu ≈ 0.0
      mktempdir() do dir
        path = joinpath(dir, "loss_extension.m")
        Sparlectra.writeMatpowerCasefile(net_loss, path)
        txt = read(path, String)
        @test occursin("SPARLECTRA EXTENSION WARNING", txt)
        @test occursin("mpc.sparlectra.format_version = 1;", txt)
        @test occursin("mpc.sparlectra.transformer_losses", txt)
        mpc_loss = Sparlectra.MatpowerIO.read_case_m(path; legacy_compat = false)
        @test mpc_loss.sparlectra !== nothing
        @test length(mpc_loss.sparlectra.transformer_losses) == 1
        net_roundtrip = Sparlectra.createNetFromMatPowerCase(mpc = mpc_loss, apply_bus_names = true, apply_branch_names = true, apply_branch_kind = true)
        @test net_roundtrip.branchVec[1].g_pu ≈ 0.005
        @test only(net_roundtrip.trafos).side1.g ≈ 0.005
        @test Sparlectra.bus_shunt_totals_pu(net_roundtrip).total_g_pu ≈ 0.0
        @test only(values(net_roundtrip.matpower_branch_metadata)).transformer_loss.g_pu ≈ 0.005
        path2 = joinpath(dir, "loss_extension_roundtrip.m")
        Sparlectra.writeMatpowerCasefile(net_roundtrip, path2)
        mpc_loss2 = Sparlectra.MatpowerIO.read_case_m(path2; legacy_compat = false)
        @test length(mpc_loss2.sparlectra.transformer_losses) == 1
        net_roundtrip2 = Sparlectra.createNetFromMatPowerCase(mpc = mpc_loss2, apply_bus_names = true, apply_branch_names = true, apply_branch_kind = true)
        @test net_roundtrip2.branchVec[1].g_pu ≈ 0.005
        @test only(net_roundtrip2.trafos).side1.g ≈ 0.005
        @test Sparlectra.bus_shunt_totals_pu(net_roundtrip2).total_g_pu ≈ 0.0
      end
    end

    dcline = [1.0 2.0 1.0 100.0 0.0 4.0 5.0 1.0 1.0 0.0 200.0 -50.0 50.0 -50.0 50.0 2.0 0.01]
    mpc_dcline = Sparlectra.MatpowerIO.legacy_sort_bus(_metadata_case(dcline = dcline))
    @test_throws Sparlectra.MatpowerIO.UnsupportedMatpowerDclineError Sparlectra.createNetFromMatPowerCase(mpc = mpc_dcline, matpower_dcline_mode = :reject_active)
    net_dcline = Sparlectra.createNetFromMatPowerCase(mpc = mpc_dcline, apply_bus_names = true, matpower_dcline_mode = :pf_injections)
    meta = only(net_dcline.matpowerDclineMetadata)
    @test meta.pf_mw == 100.0
    @test meta.input_pt_mw == 0.0
    @test meta.effective_pt_mw == 97.0
    @test meta.pt_mw == 97.0
    @test meta.loss_mw == 3.0
    @test meta.qf_mvar == 4.0
    @test meta.qt_mvar == 5.0
    @test meta.vf_pu == 1.0
    @test meta.vt_pu == 1.0
    @test meta.qminf_mvar == -50.0
    @test meta.qmaxf_mvar == 50.0
    @test meta.qmint_mvar == -50.0
    @test meta.qmaxt_mvar == 50.0
    @test meta.from_voltage_controlled === false
    @test meta.to_voltage_controlled === true
    @test net_dcline.prosumpsVec[meta.from_prosumer].pVal == -100.0
    @test net_dcline.prosumpsVec[meta.to_prosumer].pVal == 97.0
    @test net_dcline.prosumpsVec[meta.from_prosumer].qVal == 4.0
    @test net_dcline.prosumpsVec[meta.to_prosumer].qVal == 5.0
    @test net_dcline.prosumpsVec[meta.from_prosumer].minQ == -50.0
    @test net_dcline.prosumpsVec[meta.to_prosumer].maxQ == 50.0
    @test Sparlectra.getNodeType(net_dcline.nodeVec[net_dcline.busDict["SlackBus"]]) == Sparlectra.Slack
    @test Sparlectra.getNodeType(net_dcline.nodeVec[net_dcline.busDict["LoadBus"]]) == Sparlectra.PV

    inactive = [1.0 2.0 0.0 100.0 99.0]
    net_inactive = Sparlectra.createNetFromMatPowerCase(mpc = Sparlectra.MatpowerIO.legacy_sort_bus(_metadata_case(dcline = inactive)), matpower_dcline_mode = :reject_active)
    @test isempty(net_inactive.matpowerDclineMetadata)

    short_dcline = [1.0 2.0 1.0 25.0 24.0 1.5 2.5 1.01 1.02 0.0 100.0 -8.0 9.0 -6.0 7.0]
    net_short = Sparlectra.createNetFromMatPowerCase(mpc = Sparlectra.MatpowerIO.legacy_sort_bus(_metadata_case(dcline = short_dcline)), apply_bus_names = true, matpower_dcline_mode = :pf_injections)
    short_meta = only(net_short.matpowerDclineMetadata)
    @test short_meta.input_pt_mw == 24.0
    @test short_meta.effective_pt_mw == 24.0
    @test short_meta.loss0_mw == 0.0
    @test short_meta.loss1 == 0.0
    @test net_short.prosumpsVec[short_meta.from_prosumer].pVal == -25.0
    @test net_short.prosumpsVec[short_meta.to_prosumer].pVal == 24.0

    unknown_dcline = [1.0 999.0 1.0 25.0 24.0]
    @test_throws ArgumentError Sparlectra.createNetFromMatPowerCase(mpc = Sparlectra.MatpowerIO.legacy_sort_bus(_metadata_case(dcline = unknown_dcline)), matpower_dcline_mode = :pf_injections)

    isolated_case = Sparlectra.MatpowerIO.MatpowerCase(
      "isolated_dcline_case",
      100.0,
      [1 3 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9;
       2 4 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9],
      [1 10.0 0.0 100.0 -100.0 1.0 100.0 1 100.0 0.0 0 0 0 0 0 0 0 0 0 0 0],
      zeros(0, 13),
      nothing,
      nothing,
      nothing,
      nothing,
      nothing,
      [1.0 2.0 1.0 5.0 4.0 0.0 0.0 1.0 1.03],
    )
    net_isolated = Sparlectra.createNetFromMatPowerCase(mpc = isolated_case, matpower_dcline_mode = :pf_injections)
    isolated_meta = only(net_isolated.matpowerDclineMetadata)
    @test isolated_meta.to_voltage_controlled === false
    @test Sparlectra.getNodeType(net_isolated.nodeVec[net_isolated.busDict["2"]]) != Sparlectra.PV

    @testset "MATPOWER export write_solution" begin
      net_ws = Sparlectra.Net(name = "write_solution_case", baseMVA = 100.0)
      Sparlectra.addBus!(net = net_ws, busName = "SLACK", vn_kV = 110.0, oBusIdx = 1)
      Sparlectra.addBus!(net = net_ws, busName = "LOAD", vn_kV = 110.0, oBusIdx = 2)
      Sparlectra.addProsumer!(net = net_ws, busName = "SLACK", type = "GENERATOR", p = 0.0, q = 0.0, referencePri = "SLACK", vm_pu = 1.0)
      Sparlectra.addProsumer!(net = net_ws, busName = "LOAD", type = "ENERGYCONSUMER", p = 10.0, q = 2.0)
      Sparlectra._addPIModelACLine_by_idx!(net = net_ws, from = 1, to = 2, r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1, ratedS = 100.0)
      ite, erg = Sparlectra.runpf!(net_ws, 20, 1e-8, 0)
      @test erg == 0
      Sparlectra.calcNetLosses!(net_ws)

      mktempdir() do dir
        path_true = joinpath(dir, "solved_true.m")
        Sparlectra.writeMatpowerCasefile(net_ws, path_true; write_solution = true)
        txt_true = read(path_true, String)
        @test occursin("%fbus\ttbus\tr\tx\tb\trateA\trateB\trateC\tratio\tangle\tstatus\tangmin\tangmax\tPF\tQF\tPT\tQT\n", txt_true)
        @test occursin("mpc.sparlectra.solution_written = 1;", txt_true)
        mpc_true = Sparlectra.MatpowerIO.read_case_m(path_true; legacy_compat = false)
        @test mpc_true.sparlectra !== nothing
        @test mpc_true.sparlectra.solution_written === true

        path_false = joinpath(dir, "solved_false.m")
        Sparlectra.writeMatpowerCasefile(net_ws, path_false; write_solution = false)
        txt_false = read(path_false, String)
        @test occursin("%fbus\ttbus\tr\tx\tb\trateA\trateB\trateC\tratio\tangle\tstatus\tangmin\tangmax\n", txt_false)
        @test !occursin("\tPF\tQF\tPT\tQT", txt_false)
        @test !occursin("mpc.sparlectra.solution_written", txt_false)
        mpc_false = Sparlectra.MatpowerIO.read_case_m(path_false; legacy_compat = false)
        @test mpc_false.sparlectra === nothing || !mpc_false.sparlectra.solution_written
        load_row = findfirst(==(2.0), mpc_false.bus[:, 1])
        @test mpc_false.bus[load_row, 8] == 1.0 # Vm reset for non-slack/PV bus
        @test mpc_false.bus[load_row, 9] == 0.0 # Va reset for non-slack/PV bus
        slack_row = findfirst(==(1.0), mpc_false.bus[:, 1])
        @test mpc_false.bus[slack_row, 8] == net_ws.nodeVec[1]._vm_pu # slack setpoint preserved
      end
    end

    @testset "Tap-changer-model marker prevents double correction on reimport" begin
      tap_bus = [1 3 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9; 2 1 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9]
      tap_gen = [1 10.0 0.0 100.0 -100.0 1.0 100.0 1 100.0 0.0 0 0 0 0 0 0 0 0 0 0 0]
      tap_branch = [1 2 0.01 0.05 0.0 100.0 0.0 0.0 0.9 0.0 1 -60.0 60.0]
      mpc_tap = Sparlectra.MatpowerIO.MatpowerCase("tap_marker_case", 100.0, tap_bus, tap_gen, tap_branch, nothing, nothing)

      net_native = Sparlectra.createNetFromMatPowerCase(mpc = mpc_tap, tap_changer_model = :impedance_correction)
      @test length(net_native.trafos) == 1
      x_pu_corrected = net_native.branchVec[1].x_pu
      @test x_pu_corrected != 0.05 # correction actually applied

      mktempdir() do dir
        path = joinpath(dir, "tap_marker.m")
        Sparlectra.writeMatpowerCasefile(net_native, path; write_solution = false)
        txt = read(path, String)
        @test occursin("mpc.sparlectra.tap_changer_model = 'impedance_correction';", txt)
        mpc_reimport = Sparlectra.MatpowerIO.read_case_m(path; legacy_compat = false)
        @test mpc_reimport.sparlectra !== nothing
        @test mpc_reimport.sparlectra.tap_changer_model == "impedance_correction"

        # Reimport with impedance_correction active in the config path too: the
        # marker must make the importer skip a second correction, so x_pu/r_pu
        # match the already-corrected native values exactly (no double factor).
        net_roundtrip = Sparlectra.createNetFromMatPowerCase(mpc = mpc_reimport, tap_changer_model = :impedance_correction)
        @test net_roundtrip.branchVec[1].x_pu ≈ x_pu_corrected
        @test net_roundtrip.branchVec[1].r_pu ≈ net_native.branchVec[1].r_pu
      end
    end

    @testset "Standard MATPOWER case without tap-changer-model marker is unaffected" begin
      no_marker_bus = [1 3 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9; 2 1 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9]
      no_marker_gen = [1 10.0 0.0 100.0 -100.0 1.0 100.0 1 100.0 0.0 0 0 0 0 0 0 0 0 0 0 0]
      no_marker_branch = [1 2 0.01 0.05 0.0 100.0 0.0 0.0 0.9 0.0 1 -60.0 60.0]
      mpc_no_marker = Sparlectra.MatpowerIO.MatpowerCase("no_marker_case", 100.0, no_marker_bus, no_marker_gen, no_marker_branch, nothing, nothing)
      @test mpc_no_marker.sparlectra === nothing

      net_ideal = Sparlectra.createNetFromMatPowerCase(mpc = mpc_no_marker, tap_changer_model = :ideal)
      net_corrected = Sparlectra.createNetFromMatPowerCase(mpc = mpc_no_marker, tap_changer_model = :impedance_correction)
      @test net_ideal.branchVec[1].x_pu ≈ 0.05
      @test net_corrected.branchVec[1].x_pu != net_ideal.branchVec[1].x_pu
    end
  end
end

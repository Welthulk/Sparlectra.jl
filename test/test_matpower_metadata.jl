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

# file: test/test_matpower_metadata.jl
# purpose: tests MATPOWER parsing and metadata import/export: tokenizer and
#          comment-stripping edge cases, the Sparlectra transformer-loss
#          extension, write_solution, and the tap-changer-model marker
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

# Y-bus of the mpc.branch rows of a written file with plain MATPOWER
# semantics: TAP (column 9) is the ratio, SHIFT (column 10) the shift in
# degrees, TAP 0 means 1; only the transformer conductance comes from the
# Sparlectra loss block, MATPOWER has no column for it (#513)
function _mp_rows_ybus(path::AbstractString)
  mpc = Sparlectra.MatpowerIO.read_case_m(path; legacy_compat = false)
  g = Dict{Int,Float64}()
  if mpc.sparlectra !== nothing
    for row in mpc.sparlectra.transformer_losses
      g[Int(row.branch_row)] = Float64(row.g_pu)
    end
  end
  n = size(mpc.bus, 1)
  Y = zeros(ComplexF64, n, n)
  for b in eachrow(mpc.bus)
    Y[Int(b[1]), Int(b[1])] += (b[5] + im * b[6]) / mpc.baseMVA
  end
  for (r, l) in enumerate(eachrow(mpc.branch))
    l[11] == 0.0 && continue
    f, t = Int(l[1]), Int(l[2])
    ys = 1 / (l[3] + im * l[4])
    ysh = (get(g, r, 0.0) + im * l[5]) / 2
    tap = (l[9] == 0.0 ? 1.0 : l[9]) * cis(deg2rad(l[10]))
    Y[f, f] += (ys + ysh) / abs2(tap)
    Y[t, t] += ys + ysh
    Y[f, t] -= ys / conj(tap)
    Y[t, f] -= ys / tap
  end
  return Y
end

function run_matpower_metadata_tests()
  @testset "Newline-separated matrix rows (RTS-GMLC style)" begin (function ()
    # MATLAB accepts bare newlines as matrix row separators; official MATPOWER
    # cases terminate rows with ';' but e.g. RTS-GMLC does not. Splitting on
    # ';' alone silently collapsed such matrices into a single row.
    mktempdir() do dir
      path = joinpath(dir, "case_newline_rows.m")
      write(path, """
function mpc = case_newline_rows
mpc.version = '2';
mpc.baseMVA = 100.0;
mpc.bus = [
\t1\t3\t0.0\t0.0\t0.0\t0.0\t1\t1.0\t0.0\t110.0\t1\t1.05\t0.95
\t2\t1\t10.0\t3.0\t0.0\t0.0\t1\t1.0\t0.0\t110.0\t1\t1.05\t0.95
];
mpc.gen = [
\t1\t10.0\t3.0\t99\t-99\t1.0\t100.0\t1\t99\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0
];
mpc.branch = [
\t1\t2\t0.01\t0.1\t0.0\t100\t100\t100\t0.0\t0.0\t1\t-360\t360
];
""")
      mpc = Sparlectra.MatpowerIO.read_case_m(path; legacy_compat = false)
      @test size(mpc.bus) == (2, 13)
      @test size(mpc.gen) == (1, 21)
      @test size(mpc.branch) == (1, 13)
      @test mpc.bus[2, 1] == 2.0
      @test mpc.bus[2, 3] == 10.0
      # Mixed style (';' plus newlines) keeps working unchanged.
      net = createNetFromMatPowerFile(filename = path)
      @test length(net.nodeVec) == 2
    end

    # MATLAB comments: '%' up to the line end, also in CRLF files and after
    # multibyte characters
    strip_comments = Sparlectra.MatpowerIO._strip_matlab_comments
    @test strip_comments("no comments\nat all\n") === "no comments\nat all\n"
    @test strip_comments("abc%x\ndef\n") == "abc\ndef\n"
    @test strip_comments("abc%x\r\ndef\n") == "abc\ndef\n"
    @test strip_comments("%full line\nkeep\n") == "\nkeep\n"
    @test strip_comments("tail%no newline") == "tail"
    @test strip_comments("a%%b\nc%d") == "a\nc"
    @test strip_comments("ümläut%kommentar\nzeile") == "ümläut\nzeile"

    mktempdir() do dir
      path = joinpath(dir, "case_tokenizer_edges.m")
      write(path, """
function mpc = case_tokenizer_edges
mpc.version = '2';
mpc.baseMVA = 100.0; % inline comment after statement
% full-line comment
mpc.bus = [
\t1\t3\t0.0\t0.0\t0.0\t0.0\t1\t1.0\t0.0\t110.0\t1\t1.05\t0.95; % row comment

\t2\t1\t10.0\t3.0\t0.0\t0.0\t1\t1.0\t0.0\t110.0\t1\t1.05\t0.95
\t3\t1\t5.0\t1.0\t0.0\t0.0\t1\t1.0\t0.0\t110.0\t1\t1.05\t0.95\t99.0\t88.0
\t4\t1\t2.0
];
mpc.gen = [
\t1\t10.0\t3.0\t99\t-99\t1.0\t100.0\t1\t99\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0
];
mpc.branch = [
\t1\t2\t0.01\t0.1\t0.0\t100\t100\t100\t0.0\t0.0\t1\t-360\t360;
\t2\t3\t0.01\t0.1\t0.0\t100\t100\t100\t0.0\t0.0\t1\t-360\t360;
\t3\t4\t0.01\t0.1\t0.0\t100\t100\t100\t0.0\t0.0\t1\t-360\t360
];
mpc.gencost = [
\t2\t0\t0\t3\t0.1\t5.0\t0.0;
\t2\t0\t0\t2\t4.0\t0.0
];
""")
      mpc = Sparlectra.MatpowerIO.read_case_m(path; legacy_compat = false)
      # blank line + comment-only rows ignored; short row zero-padded;
      # overlong row truncated to the fixed 13 columns
      @test size(mpc.bus) == (4, 13)
      @test mpc.bus[3, 13] == 0.95
      @test mpc.bus[4, 3] == 2.0
      @test mpc.bus[4, 4] == 0.0
      # auto-width block: ragged rows padded with 0.0 to maxcols
      @test size(mpc.gencost) == (2, 7)
      @test mpc.gencost[2, 5] == 4.0
      @test mpc.gencost[2, 7] == 0.0
    end
  end)() end

  @testset "MATPOWER metadata import" begin (function ()
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

    # missing bus_name with apply_bus_names (0.60.0): the fallback warns once
    # per case and process, a repeated import of the same case logs it at
    # debug level only (the Web UI imports a case for the start page and for
    # every run), a different case warns again. Case names unique to this
    # check, so no earlier import has used up the warning.
    nobus(case) = Sparlectra.MatpowerIO.MatpowerCase(case, mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, nothing, nothing)
    fallback = (:warn, r"bus_name metadata missing")
    @test_logs fallback match_mode = :any Sparlectra.createNetFromMatPowerCase(mpc = nobus("bus_name_once_a"), apply_bus_names = true)
    @test_logs min_level = Logging.Warn Sparlectra.createNetFromMatPowerCase(mpc = nobus("bus_name_once_a"), apply_bus_names = true)
    @test_logs (:debug, r"bus_name metadata missing") match_mode = :any min_level = Logging.Debug Sparlectra.createNetFromMatPowerCase(mpc = nobus("bus_name_once_a"), apply_bus_names = true)
    @test_logs fallback match_mode = :any Sparlectra.createNetFromMatPowerCase(mpc = nobus("bus_name_once_b"), apply_bus_names = true)

    # line or transformer without a branch kind (0.60.0): TAP 0 or exactly 1,
    # SHIFT 0 and one voltage level is a line (case6495rte writes TAP 1 on
    # every line); TAP 1 between 110 and 20 kV is a transformer with ratio 1,
    # TAP 1 with a shift and any other TAP are transformers. Both readers
    # agree, an explicit kind still wins, and the Y-bus equals the one of the
    # rows with plain MATPOWER semantics (the kind changes no admittance).
    mktempdir() do dir
      path_kind = joinpath(dir, "kind_rule.m")
      write(path_kind, """
      function mpc = kind_rule
      mpc.version = '2';
      mpc.baseMVA = 100;
      mpc.bus = [
      1 3 0 0 0 0 1 1.0 0 110 1 1.1 0.9;
      2 1 20 5 0 0 1 1.0 0 110 1 1.1 0.9;
      3 1 10 2 0 0 1 1.0 0 20 1 1.1 0.9;
      ];
      mpc.gen = [
      1 30 0 100 -100 1.0 100 1 100 0 0 0 0 0 0 0 0 0 0 0 0;
      ];
      mpc.branch = [
      1 2 0.01 0.05 0.02 100 0 0 1 0 1 -360 360;
      2 3 0.005 0.08 0 100 0 0 1 0 1 -360 360;
      1 2 0.01 0.06 0 100 0 0 1 2 1 -360 360;
      1 2 0.02 0.07 0.01 100 0 0 0 0 1 -360 360;
      1 2 0.01 0.05 0 100 0 0 0.98 0 1 -360 360;
      ];
      """)
      cfg_kind = Sparlectra.active_sparlectra_config()
      net_dir = Sparlectra.createNetFromMatPowerFile(filename = path_kind)
      net_conv = Sparlectra.build_net(Sparlectra.convert_case(Sparlectra.MatpowerAdapter(), Sparlectra.MatpowerIO.read_case(path_kind; legacy_compat = true), Sparlectra.matpower_adapter_options(cfg_kind)); config = cfg_kind)
      for n in (net_dir, net_conv)
        @test [br.ratio != 0.0 for br in n.branchVec] == [false, true, true, false, true]
        @test n.branchVec[2].ratio == 1.0
        @test maximum(abs.(Matrix(createYBUS(net = n, sparse = true)) .- _mp_rows_ybus(path_kind))) <= 1e-12
      end
      # an explicit kind overrides the rule in both directions
      mpc_kind = Sparlectra.MatpowerIO.read_case(path_kind)
      mpc_kind = Sparlectra.MatpowerIO.MatpowerCase(mpc_kind.name, mpc_kind.baseMVA, mpc_kind.bus, mpc_kind.gen, mpc_kind.branch, nothing, nothing, nothing, ["T", "T", "T", "L", "T"], nothing, nothing)
      net_k = Sparlectra.createNetFromMatPowerCase(mpc = mpc_kind, apply_branch_kind = true)
      @test [br.ratio != 0.0 for br in net_k.branchVec] == [true, true, true, false, true]
    end

    @testset "Sparlectra transformer-loss extension" begin (function ()
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
        Sparlectra.writeMatpowerCasefile(net_loss, path; write_solution = false)
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
        Sparlectra.writeMatpowerCasefile(net_roundtrip, path2; write_solution = false)
        mpc_loss2 = Sparlectra.MatpowerIO.read_case_m(path2; legacy_compat = false)
        @test length(mpc_loss2.sparlectra.transformer_losses) == 1
        net_roundtrip2 = Sparlectra.createNetFromMatPowerCase(mpc = mpc_loss2, apply_bus_names = true, apply_branch_names = true, apply_branch_kind = true)
        @test net_roundtrip2.branchVec[1].g_pu ≈ 0.005
        @test only(net_roundtrip2.trafos).side1.g ≈ 0.005
        @test Sparlectra.bus_shunt_totals_pu(net_roundtrip2).total_g_pu ≈ 0.0
      end
    end)() end

    @testset "Sparlectra links extension (busbar coupler)" begin (function ()
      mktempdir() do dir
        # Bus 2 carries the load but has no branch: it is fed only through
        # the closed impedance-less link onto busbar bus 3, so the test
        # covers the galvanic link contraction end to end.
        path = joinpath(dir, "case_link.m")
        write(path, """
function mpc = case_link
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
1 3 0 0 0 0 1 1 0 110 1 1.1 0.9;
2 1 40 10 0 0 1 1 0 110 1 1.1 0.9;
3 1 0 0 0 0 1 1 0 110 1 1.1 0.9;
];
mpc.gen = [
1 50 0 100 -100 1 100 1 100 0 0 0 0 0 0 0 0 0 0 0 0;
];
mpc.branch = [
1 3 0.01 0.05 0 100 0 0 0 0 1 -360 360;
];
mpc.sparlectra = struct();
mpc.sparlectra.format_version = 1;
mpc.sparlectra.links = [
2 3 1;
];
""")
        mpc_link = Sparlectra.MatpowerIO.read_case_m(path; legacy_compat = false)
        @test mpc_link.sparlectra !== nothing
        @test mpc_link.sparlectra.links == [2.0 3.0 1.0]

        net_link = Sparlectra.createNetFromMatPowerFile(filename = path)
        @test length(net_link.linkVec) == 1
        l = only(net_link.linkVec)
        @test l.status == 1
        @test Set([Int(l.fromBus), Int(l.toBus)]) == Set([net_link.busDict["2"], net_link.busDict["3"]])

        _, erg = runpf!(net_link, 30, 1e-8, 0)
        @test erg == 0
        # closed link => the two buses are one electrical node
        @test Sparlectra.get_bus_vm_pu(net_link, "2") ≈ Sparlectra.get_bus_vm_pu(net_link, "3") atol = 1e-12

        path2 = joinpath(dir, "case_link_roundtrip.m")
        Sparlectra.writeMatpowerCasefile(net_link, path2; write_solution = false)
        txt = read(path2, String)
        @test occursin("mpc.sparlectra.links", txt)
        net_link2 = Sparlectra.createNetFromMatPowerFile(filename = path2)
        @test length(net_link2.linkVec) == 1
        @test only(net_link2.linkVec).status == 1

        # malformed block: wrong column count must be rejected, not guessed
        bad = joinpath(dir, "case_link_bad.m")
        write(bad, replace(read(path, String), "2 3 1;" => "2 3;"))
        @test_throws ArgumentError Sparlectra.MatpowerIO.read_case_m(bad; legacy_compat = false)
      end
    end)() end

    @testset "Sparlectra tap-changers extension (nameplate data)" begin (function ()
      mktempdir() do dir
        # branch 2 is a transformer at neutral ratio 1.0; the extension
        # declares a +-2 x 2.5 percent changer standing at step 1 and a
        # pure phase shifter on branch 3 (tap_step 0 = no ratio changer)
        path = joinpath(dir, "case_taps.m")
        write(path, """
function mpc = case_taps
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
1 3 0 0 0 0 1 1 0 110 1 1.1 0.9;
2 1 40 10 0 0 1 1 0 20 1 1.1 0.9;
3 1 20 5 0 0 1 1 0 20 1 1.1 0.9;
];
mpc.gen = [
1 70 0 100 -100 1 100 1 100 0 0 0 0 0 0 0 0 0 0 0 0;
];
mpc.branch = [
1 2 0.002 0.06 0 100 0 0 1.0 0 1 -360 360;
1 3 0.002 0.06 0 100 0 0 1.0 -2.0 1 -360 360;
2 3 0.02 0.1 0 100 0 0 0 0 1 -360 360;
];
mpc.sparlectra = struct();
mpc.sparlectra.format_version = 1;
mpc.sparlectra.tap_changers = [
1 0.025 -2 2 1 0 0 0 0;
2 0 0 0 0 0.5 -10 10 -2;
];
""")
        net = Sparlectra.createNetFromMatPowerFile(filename = path)
        b1 = net.branchVec[1]
        @test b1.has_ratio_tap
        @test b1.tap_step == 0.025
        @test b1.ratio == 1.0                       # neutral stays the TAP column
        @test b1.tap_ratio ≈ 1.0 / 1.025            # standing at step 1
        @test b1.tap_min ≈ 1.0 / 1.05
        @test b1.tap_max ≈ 1.0 / 0.95
        b2 = net.branchVec[2]
        @test !b2.has_ratio_tap                     # tap_step 0: declared ratio-less
        @test b2.has_phase_tap
        @test b2.phase_step_deg == 0.5
        @test b2.angle == -2.0                      # neutral shift stays the SHIFT column
        @test b2.phase_shift_deg ≈ -3.0             # -2 deg neutral + (-2 steps) * 0.5 deg
        # the band is stored RELATIVE to neutral (the fixation compares the
        # regulating-vector angle, which passes 0 at the neutral position)
        @test b2.phase_min_deg ≈ -5.0 && b2.phase_max_deg ≈ 5.0
        @test !b1.ratio_step_assumed && !b2.phase_step_assumed

        # a NaN current step is no step: the direct import
        # and the adapter put the changer on its middle position, the
        # neutral position 0 (TAP and SHIFT of the branch row), mark it and
        # name it in one info line per transformer
        pathn = joinpath(dir, "case_taps_nostep.m")
        write(pathn, replace(read(path, String), "1 0.025 -2 2 1 0 0 0 0;" => "1 0.025 -2 2 NaN 0 0 0 0;", "2 0 0 0 0 0.5 -10 10 -2;" => "2 0 0 0 0 0.5 -10 10 NaN;"))
        netn = @test_logs (:info, "MATPOWER: tap changer branch 1: step not transferred, middle position 0 assumed") (:info, "MATPOWER: tap changer branch 2: step not transferred, middle position 0 assumed") match_mode = :any Sparlectra.createNetFromMatPowerFile(filename = pathn)
        cfgn = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
        line_n = r"^SCF: tap changer branch \d+: step not transferred, middle position 0 assumed$"
        convn = @test_logs (:info, line_n) (:info, line_n) match_mode = :any Sparlectra.build_net(Sparlectra.convert_case(MatpowerAdapter(), Sparlectra.MatpowerIO.read_case(pathn), Sparlectra.matpower_adapter_options(cfgn)); config = cfgn)
        for back in (netn, convn)
          n1 = only(b for b in back.branchVec if b.has_ratio_tap && b.ratio_tap_declared)
          n2 = only(b for b in back.branchVec if b.has_phase_tap)
          @test (n1.tap_ratio, n1.ratio_step_assumed, n1.phase_step_assumed) == (1.0, true, false)
          @test (n2.phase_shift_deg, n2.phase_step_assumed, n2.ratio_step_assumed) == (-2.0, true, false)
        end
        # MATPOWER has no mark for it: the export writes the middle position
        # as the step with one warning per transformer, and the reimport
        # reads it as a transferred step with the same taps
        pathw = joinpath(dir, "case_taps_nostep_out.m")
        assumed_warn(br, which) = "MATPOWER export: transformer $(br.comp.cName): $(which) step not transferred, the export writes the assumed middle position"
        w1 = only(b for b in netn.branchVec if b.ratio_step_assumed)
        w2 = only(b for b in netn.branchVec if b.phase_step_assumed)
        wlogs, _ = Test.collect_test_logs() do
          Sparlectra.writeMatpowerCasefile(netn, pathw; write_solution = false)
        end
        wmsgs = [string(l.message) for l in wlogs if l.level == Base.CoreLogging.Warn && occursin("step not transferred", string(l.message))]
        @test sort(wmsgs) == sort([assumed_warn(w1, "ratio"), assumed_warn(w2, "phase")])
        netw = Sparlectra.createNetFromMatPowerFile(filename = pathw)
        @test !any(b -> b.ratio_step_assumed || b.phase_step_assumed, netw.branchVec)
        @test [(b.tap_ratio, b.phase_shift_deg) for b in netw.branchVec] == [(b.tap_ratio, b.phase_shift_deg) for b in netn.branchVec]
        # a tap controller that moves the assumed ratio changer takes its
        # mark off: the step it reaches is a step, only the phase changer
        # is still warned about
        netc = deepcopy(netn)
        wc1 = only(b for b in netc.branchVec if b.ratio_step_assumed)
        Sparlectra.addPowerTransformerControl!(netc; trafo = string(wc1.branchIdx), mode = :voltage, target_bus = only(k for (k, v) in netc.busDict if v == Int(wc1.toBus)), target_vm_pu = 1.04, deadband_vm_pu = 2e-3)
        @test Sparlectra.run_sparlectra(net = netc, config = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-9), output = OutputConfig(logfile_results = :off, startup_latency_hint = false))).numerical_converged
        @test wc1.tap_ratio != 1.0 && !wc1.ratio_step_assumed
        clogs, _ = Test.collect_test_logs() do
          Sparlectra.writeMatpowerCasefile(netc, joinpath(dir, "case_taps_ctl_out.m"); write_solution = false)
        end
        @test [string(l.message) for l in clogs if l.level == Base.CoreLogging.Warn && occursin("step not transferred", string(l.message))] == [assumed_warn(w2, "phase")]

        # the live positions enter the power flow and the generated flows
        # stay consistent with the injections (branchFlow_pu live-tap rule)
        _, terg = runpf!(net, 40, 1e-10, 0)
        @test terg == 0
        tmeas = generateMeasurementsFromPF(net; noise = false)
        for bus in eachindex(net.nodeVec)
          pinj = only(m.value for m in tmeas if m.typ == Sparlectra.PinjMeas && m.busIdx == bus)
          pfl = sum(m.value for m in tmeas if m.typ == Sparlectra.PflowMeas && ((m.direction == :from && Int(net.branchVec[m.branchIdx].fromBus) == bus) || (m.direction == :to && Int(net.branchVec[m.branchIdx].toBus) == bus)); init = 0.0)
          @test isapprox(pinj, pfl; atol = 1e-8)
        end

        # export writes the nameplate rows back; the reimport restores the
        # grid and the live positions
        path2 = joinpath(dir, "case_taps_roundtrip.m")
        Sparlectra.writeMatpowerCasefile(net, path2; write_solution = false)
        @test occursin("mpc.sparlectra.tap_changers", read(path2, String))
        net2 = Sparlectra.createNetFromMatPowerFile(filename = path2)
        @test net2.branchVec[1].tap_step == 0.025
        @test net2.branchVec[1].tap_ratio ≈ b1.tap_ratio
        @test !net2.branchVec[2].has_ratio_tap
        @test net2.branchVec[2].phase_shift_deg ≈ -3.0
        # a transformer at ratio 1 between two buses of one level, without
        # MATPOWER metadata (as from a CGMES import, the 380 kV legs of a
        # three-winding star), goes to mpc.branch with TAP 1 and "T", comes
        # back as a transformer and keeps its nameplate row (until the
        # ratio-based kind it was written as a line with TAP 0 and no row)
        path110 = joinpath(dir, "case_taps_110.m")
        write(path110, replace(read(path, String), "2 1 40 10 0 0 1 1 0 20 1 1.1 0.9;" => "2 1 40 10 0 0 1 1 0 110 1 1.1 0.9;"))
        net3 = Sparlectra.createNetFromMatPowerFile(filename = path110)
        delete!(net3.matpower_branch_metadata, 1)
        net3.branchVec[1].tap_ratio = 1.0
        path3 = joinpath(dir, "case_taps_ratio1.m")
        Sparlectra.writeMatpowerCasefile(net3, path3; write_solution = false)
        mpc_r1 = Sparlectra.MatpowerIO.read_case_m(path3; legacy_compat = false)
        @test (mpc_r1.branch[1, 9], mpc_r1.branch_kind[1]) == (1.0, "T")
        net3b = Sparlectra.createNetFromMatPowerFile(filename = path3)
        @test net3b.branchVec[1].ratio == 1.0
        @test net3b.branchVec[1].tap_step == 0.025
        # 14 columns since 0.50.0: the neutral position in 13 and 14 (#513)
        @test size(Sparlectra._matpower_sparlectra_tap_changers(Sparlectra.MatpowerIO.read_case(path2)), 2) == 14

        # a CGMES changer (tap_ratio_grid :linear) writes its positions on
        # its own grid with ratio_grid 1 in column 12, and both readers bring
        # the grid back. Before, the 3W CGMES fixture 8 steps of 1 % off its
        # import ratio wrote position -7.41 on the fraction grid and came
        # back as a :fraction branch.
        netl = Sparlectra.createNetFromMatPowerFile(filename = path)
        Sparlectra.applyTapNameplate!(netl.branchVec[1]; tap_step = 0.01, tap_min_step = -10.0, tap_max_step = 10.0, tap_current_step = 8.0, tap_ratio_grid = :linear)
        @test netl.branchVec[1].tap_ratio == 1.08
        pathl = joinpath(dir, "case_taps_linear.m")
        Sparlectra.writeMatpowerCasefile(netl, pathl; write_solution = false)
        mpl = Sparlectra.MatpowerIO.read_case(pathl)
        tcl = Sparlectra._matpower_sparlectra_tap_changers(mpl)
        @test size(tcl, 2) == 14
        @test (tcl[1, 3], tcl[1, 4], tcl[1, 5], tcl[1, 12]) == (-10.0, 10.0, 8.0, 1.0)
        @test tcl[2, 12] == 0.0
        cfgl = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
        for back in (Sparlectra.createNetFromMatPowerFile(filename = pathl), Sparlectra.build_net(Sparlectra.convert_case(MatpowerAdapter(), mpl, Sparlectra.matpower_adapter_options(cfgl)); config = cfgl))
          bl = only(b for b in back.branchVec if b.has_ratio_tap)
          @test bl.tap_ratio_grid === :linear
          @test (bl.tap_ratio, bl.tap_min, bl.tap_max, bl.tap_step) == (1.08, 0.9, 1.1, 0.01)
          @test only(b for b in back.branchVec if b.has_phase_tap).tap_ratio_grid === :fraction
        end

        # Delta-u PST (11th column): the mechanical grid is the
        # additional-voltage amplitude, the shift angle follows via atan
        pathdu = joinpath(dir, "case_taps_du.m")
        write(pathdu, replace(read(path, String), "2 0 0 0 0 0.5 -10 10 -2;" => "2 0 0 0 0 0 -10 10 -2 90 0.01;"))
        netdu = Sparlectra.createNetFromMatPowerFile(filename = pathdu)
        bdu = netdu.branchVec[2]
        @test bdu.has_phase_tap && !bdu.has_ratio_tap
        @test bdu.phase_du_step == 0.01
        @test bdu.tap_est_alpha_deg == 90.0
        # live polar position from the cascade at r2 = -2 * 0.01
        tb = 1.0 * cis(deg2rad(-2.0))
        tl = tb / (1.0 + (-2 * 0.01) * cis(deg2rad(90.0)))
        @test bdu.tap_ratio ≈ abs(tl)
        @test bdu.phase_shift_deg ≈ rad2deg(angle(tl))
        # roundtrip keeps the Delta-u grid and the live position
        pathdu2 = joinpath(dir, "case_taps_du_rt.m")
        Sparlectra.writeMatpowerCasefile(netdu, pathdu2; write_solution = false)
        netdu2 = Sparlectra.createNetFromMatPowerFile(filename = pathdu2)
        @test netdu2.branchVec[2].phase_du_step == 0.01
        @test netdu2.branchVec[2].tap_ratio ≈ bdu.tap_ratio
        @test netdu2.branchVec[2].phase_shift_deg ≈ bdu.phase_shift_deg
        # a release and a disable keep psi 90; a psi a caller sets to 0 on
        # this 90 model at step -2: the writer splits along 90 and the reader
        # builds 90, so the round trip is exact and the writer gives no psi-0
        # warning
        netdc = Sparlectra.createNetFromMatPowerFile(filename = pathdu)
        Sparlectra.setTapEstimation!(netdc; trafo = 2, mode = :pst, alpha_deg = 90.0)
        Sparlectra.setTapEstimation!(netdc; trafo = 2, enabled = false)
        @test netdc.branchVec[2].tap_est_alpha_deg == 90.0
        netdc.branchVec[2].tap_est_alpha_deg = 0.0
        pathdc = joinpath(dir, "case_taps_du_clr.m")
        @test_logs min_level = Logging.Warn Sparlectra.writeMatpowerCasefile(netdc, pathdc; write_solution = false)
        bdc = Sparlectra.createNetFromMatPowerFile(filename = pathdc).branchVec[2]
        @test (bdc.phase_du_step, bdc.tap_est_alpha_deg) == (0.01, 90.0)
        @test bdc.tap_ratio ≈ bdu.tap_ratio
        @test bdc.phase_shift_deg ≈ bdu.phase_shift_deg
        # SCF with alpha_deg 0 stated for this changer: 0 means 0, a Delta-u
        # changer in phase with the voltage. Hand calculation: step -2 of
        # 0.01 at psi 0 gives t = cis(-2 deg) / (1 - 0.02). The SCF writer
        # gives step -2 and alpha_deg 0 back. Until 0.60.0 the import built
        # the model along 90 and stored 0, and the writer wrote step 0. A
        # missing alpha_deg key still reads as 90, as in the shipped cases.
        scfdu = read(exportSCF(netdu; file = joinpath(dir, "case_taps_du.scf.json")), String)
        @test count(r"\"alpha_deg\":\s*90\.0", scfdu) == 1
        scfdu0 = joinpath(dir, "case_taps_du_0.scf.json")
        write(scfdu0, replace(scfdu, r"\"alpha_deg\":\s*90\.0" => "\"alpha_deg\": 0.0"))
        netdu0 = importSCF(scfdu0)
        bdu0 = only(b for b in netdu0.branchVec if b.has_phase_tap && b.phase_du_step > 0.0)
        tl0 = cis(deg2rad(-2.0)) / (1.0 - 0.02)
        @test bdu0.tap_est_alpha_deg == 0.0
        @test bdu0.tap_ratio ≈ abs(tl0)
        @test bdu0.phase_shift_deg ≈ rad2deg(angle(tl0))
        tcdu0 = Sparlectra.net_to_scf(netdu0)["sparlectra"]["components"]["tap_changer"]
        cdu0 = only(c for t in tcdu0 for c in t["controllers"] if c["index"] == 2)
        @test (Int(cdu0["pos"]), cdu0["alpha_deg"]) == (-2, 0.0)
        bdu0rt = only(b for b in importSCF(exportSCF(netdu0; file = joinpath(dir, "case_taps_du_0_rt.scf.json"))).branchVec if b.has_phase_tap && b.phase_du_step > 0.0)
        @test (bdu0rt.tap_est_alpha_deg, bdu0rt.phase_du_step) == (0.0, 0.01)
        @test bdu0rt.tap_ratio ≈ bdu0.tap_ratio
        @test bdu0rt.phase_shift_deg ≈ bdu0.phase_shift_deg
        # MATPOWER reads psi 0 as 90 and cannot carry this changer: the
        # writer refuses it instead of writing step 0
        @test_throws ArgumentError Sparlectra.writeMatpowerCasefile(netdu0, joinpath(dir, "case_taps_du_0.m"); write_solution = false)
        scfdum = joinpath(dir, "case_taps_du_nokey.scf.json")
        write(scfdum, replace(scfdu, r"\"alpha_deg\":\s*90\.0,\s*" => ""))
        bdum = only(b for b in importSCF(scfdum).branchVec if b.has_phase_tap && b.phase_du_step > 0.0)
        @test bdum.tap_est_alpha_deg == 90.0
        @test bdum.tap_ratio ≈ bdu.tap_ratio
        @test bdum.phase_shift_deg ≈ bdu.phase_shift_deg
        # the same stated 0 at position 0: t is the neutral tap, the MATPOWER
        # writer passes it, but its psi column knows only 0 = not given and
        # the reader builds 90. The writer warns; the direction comes back
        # as 90 (known gap, no format change).
        scfdz = joinpath(dir, "case_taps_du_0_pos0.scf.json")
        rdz = r"\"alpha_deg\":\s*90\.0,(\s*\"index\":\s*2,\s*\"pos\":\s*)-2"
        @test count(rdz, scfdu) == 1
        write(scfdz, replace(scfdu, rdz => s"\"alpha_deg\": 0.0,\g<1>0"))
        netdz = importSCF(scfdz)
        bdz = only(b for b in netdz.branchVec if b.has_phase_tap && b.phase_du_step > 0.0)
        @test bdz.tap_est_alpha_deg == 0.0
        @test bdz.tap_ratio ≈ 1.0 && bdz.phase_shift_deg ≈ -2.0
        pathdz = joinpath(dir, "case_taps_du_0_pos0.m")
        @test_logs (:warn, r"Delta-u changer with psi 0") match_mode = :any Sparlectra.writeMatpowerCasefile(netdz, pathdz; write_solution = false)
        bdz2 = only(b for b in Sparlectra.createNetFromMatPowerFile(filename = pathdz).branchVec if b.has_phase_tap && b.phase_du_step > 0.0)
        @test bdz2.tap_est_alpha_deg == 90.0
        @test bdz2.tap_ratio ≈ bdz.tap_ratio
        @test bdz2.phase_shift_deg ≈ bdz.phase_shift_deg
        # a degree-form phase changer as the CGMES import leaves it: psi
        # not set, tap_est_alpha_deg at the constructor 0 (_setPhaseTapRange!
        # does not touch it), which is also what a stated SCF alpha_deg 0
        # stores. The SCF writer writes the degree form with alpha_deg 0 and
        # the reader keeps kind, neutral shift, positions and the 0. MATPOWER
        # writes psi 0 and the reader builds 90, so the changer comes back
        # declared (known gap until S5).
        netg = Sparlectra.createNetFromMatPowerFile(filename = path)
        bg = netg.branchVec[2]
        bg.tap_est_alpha_deg = 0.0
        tcg = only(c for t in Sparlectra.net_to_scf(netg)["sparlectra"]["components"]["tap_changer"] for c in t["controllers"] if c["index"] == 2)
        @test (tcg["alpha_deg"], tcg["step_deg"], Int(tcg["pos"]), Int(tcg["pos_min"]), Int(tcg["pos_max"]), haskey(tcg, "step")) == (0.0, 0.5, -2, -10, 10, false)
        bg2 = only(b for b in importSCF(exportSCF(netg; file = joinpath(dir, "case_taps_deg0.scf.json"))).branchVec if b.has_phase_tap && b.ratio != 0.0 && b.phase_step_deg == 0.5)
        @test (bg2.phase_du_step, bg2.tap_est_alpha_deg, bg2.angle) == (0.0, 0.0, -2.0)
        @test bg2.phase_shift_deg ≈ -3.0
        @test bg2.phase_min_deg ≈ -5.0 && bg2.phase_max_deg ≈ 5.0
        pathg = joinpath(dir, "case_taps_deg0.m")
        Sparlectra.writeMatpowerCasefile(netg, pathg; write_solution = false)
        bg3 = Sparlectra.createNetFromMatPowerFile(filename = pathg).branchVec[2]
        @test (bg3.phase_step_deg, bg3.angle, bg3.tap_est_alpha_deg) == (0.5, -2.0, 90.0)
        @test bg3.phase_shift_deg ≈ -3.0
        # both phase grids at once are contradictory nameplate data
        badgrid = joinpath(dir, "case_taps_badgrid.m")
        write(badgrid, replace(read(path, String), "2 0 0 0 0 0.5 -10 10 -2;" => "2 0 0 0 0 0.5 -10 10 -2 90 0.01;"))
        @test_throws ArgumentError Sparlectra.createNetFromMatPowerFile(filename = badgrid)

        # ratio changer plus Delta-u changer on one branch: the Delta-u factor
        # divides the LIVE ratio. Hand calculation: ratio step 2 of 2.5 %
        # gives 1 / 1.05, Delta-u step -10 of 0.01 at psi 90 gives 1 - 0.1j,
        # so t = (1 / 1.05) * cis(-2 deg) / (1 - 0.1j). Until 0.60.0 the
        # import built t on the neutral ratio 1 and lost the ratio position;
        # the writers and the generator read both positions from t as if the
        # branch had one changer (the SCF writer gave Delta-u step -11).
        pathpair = joinpath(dir, "case_taps_pair.m")
        write(pathpair, replace(read(path, String), "2 0 0 0 0 0.5 -10 10 -2;" => "2 0.025 -2 2 2 0 -10 10 -10 90 0.01;"))
        netp = Sparlectra.createNetFromMatPowerFile(filename = pathpair)
        bp = netp.branchVec[2]
        tpair = (1.0 / 1.05) * cis(deg2rad(-2.0)) / (1.0 - 0.1im)
        @test bp.has_ratio_tap && bp.has_phase_tap
        @test bp.tap_ratio ≈ abs(tpair)
        @test bp.phase_shift_deg ≈ rad2deg(angle(tpair))
        # MATPOWER writer: both positions in their own columns, reimport
        # gives the same branch tap
        pathpair2 = joinpath(dir, "case_taps_pair_rt.m")
        Sparlectra.writeMatpowerCasefile(netp, pathpair2; write_solution = false)
        tcp = Sparlectra._matpower_sparlectra_tap_changers(Sparlectra.MatpowerIO.read_case(pathpair2))
        rowp = only(r for r in eachrow(tcp) if r[1] == 2.0)
        @test rowp[5] ≈ 2.0 && rowp[9] ≈ -10.0
        netp2 = Sparlectra.createNetFromMatPowerFile(filename = pathpair2)
        @test netp2.branchVec[2].tap_ratio ≈ bp.tap_ratio
        @test netp2.branchVec[2].phase_shift_deg ≈ bp.phase_shift_deg
        # SCF writer: controller 1 at position 2, controller 2 at -10, and
        # the reimport gives the same branch tap
        tcs = Sparlectra.net_to_scf(netp)["sparlectra"]["components"]["tap_changer"]
        ctrl = only(t["controllers"] for t in tcs if length(t["controllers"]) == 2)
        @test sort([(Int(c["index"]), Int(c["pos"])) for c in ctrl]) == [(1, 2), (2, -10)]
        netp3 = importSCF(exportSCF(netp; file = joinpath(dir, "case_taps_pair.scf.json")))
        bp3 = only(b for b in netp3.branchVec if b.has_ratio_tap && b.has_phase_tap)
        @test bp3.tap_ratio ≈ bp.tap_ratio
        @test bp3.phase_shift_deg ≈ bp.phase_shift_deg
        # SE generator: one ratio step down on both transformers. The pair
        # keeps its Delta-u position. Truth: the same net with branch 1 at
        # ratio 1 (step 0) and branch 2 at t = (1 / 1.025) * cis(-2 deg) /
        # (1 - 0.1j), set by hand; the branch flows have to agree.
        gen = SparlectraApp._se_generate_measurement_set(pathpair, joinpath(dir, "pair.measurements.csv"), SparlectraApp.MeasurementGeneratorOptions(; noise = false, gross_k = 0.0, tap_steps = -1.0, tap_count = 2, include_i = false, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0, sigma_ia_deg = 0.0))
        nett = Sparlectra.createNetFromMatPowerFile(filename = pathpair)
        tgen = (1.0 / 1.025) * cis(deg2rad(-2.0)) / (1.0 - 0.1im)
        nett.branchVec[1].tap_ratio = 1.0
        nett.branchVec[2].tap_ratio = abs(tgen)
        nett.branchVec[2].phase_shift_deg = rad2deg(angle(tgen))
        _, terr = runpf!(nett, 40, 1e-10, 0)
        @test terr == 0
        tmp = generateMeasurementsFromPF(nett; noise = false)
        for k in (1, 2), typ in (Sparlectra.PflowMeas, Sparlectra.QflowMeas)
          want = only(m.value for m in tmp if m.typ == typ && m.branchIdx == k && m.direction == :from)
          got = only(m.value for m in gen.measurements if m.typ == typ && m.branchIdx == k && m.direction == :from)
          @test isapprox(got, want; atol = 1e-6)
        end
        # psi 180 puts the Delta-u factor on the real axis: the two
        # positions cannot be told apart, so the import refuses the pair
        pathcol = joinpath(dir, "case_taps_pair_180.m")
        write(pathcol, replace(read(path, String), "2 0 0 0 0 0.5 -10 10 -2;" => "2 0.025 -2 2 1 0 -10 10 -2 180 0.01;"))
        @test_throws ArgumentError Sparlectra.createNetFromMatPowerFile(filename = pathcol)
        # SCF keeps a stated alpha_deg 0 (MATPOWER reads 0 as 90), which is
        # the same collinear pair: the SCF import refuses it as well
        scfpair = read(exportSCF(netp; file = joinpath(dir, "case_taps_pair_a.scf.json")), String)
        scf0 = joinpath(dir, "case_taps_pair_0.scf.json")
        write(scf0, replace(scfpair, r"\"alpha_deg\":\s*90\.0" => "\"alpha_deg\": 0.0"))
        @test_throws ArgumentError importSCF(scf0)
        # setTapEstimation! with enabled = false and mode = :ratio keeps the
        # nameplate psi (30 here), so both writers split the pair at either
        # Delta-u position and the reimport gives the same branch tap. Until
        # 0.60.0 both calls set psi to 0 and the writers refused the pair
        # off Delta-u position 0.
        # A psi a caller sets to 0 is gone. Off Delta-u position 0 no writer
        # can split t without it: both refuse instead of writing positions
        # along psi 0 (SCF) or a guessed 90 (MATPOWER, ratio step 2 and
        # Delta-u step -10 at psi 30 came out as other positions). At
        # Delta-u position 0 the factor is 1, the split is exact and the SCF
        # file with alpha_deg 0 reads back.
        for (dupos, writable) in ((-10, false), (0, true))
          pathclr = joinpath(dir, "case_taps_pair_clr.m")
          write(pathclr, replace(read(path, String), "2 0 0 0 0 0.5 -10 10 -2;" => "2 0.025 -2 2 2 0 -10 10 $(dupos) 30 0.01;"))
          netc = Sparlectra.createNetFromMatPowerFile(filename = pathclr)
          Sparlectra.setTapEstimation!(netc; trafo = 2, mode = :both, alpha_deg = 30.0)
          Sparlectra.setTapEstimation!(netc; trafo = 2, enabled = false)
          # :ratio is consistent with the branch state only at Delta-u position 0
          dupos == 0 && Sparlectra.setTapEstimation!(netc; trafo = 2, mode = :ratio)
          @test netc.branchVec[2].tap_est_alpha_deg == 30.0
          pathkeep = joinpath(dir, "case_taps_pair_keep_rt.m")
          Sparlectra.writeMatpowerCasefile(netc, pathkeep; write_solution = false)
          for bk in (Sparlectra.createNetFromMatPowerFile(filename = pathkeep).branchVec[2], only(b for b in importSCF(exportSCF(netc; file = joinpath(dir, "case_taps_pair_keep.scf.json"))).branchVec if b.has_ratio_tap && b.has_phase_tap))
            @test bk.tap_est_alpha_deg == 30.0
            @test bk.tap_ratio ≈ netc.branchVec[2].tap_ratio
            @test bk.phase_shift_deg ≈ netc.branchVec[2].phase_shift_deg
          end
          netc.branchVec[2].tap_est_alpha_deg = 0.0
          pathclr2 = joinpath(dir, "case_taps_pair_clr_rt.m")
          if writable
            # the psi 0 goes out as not given (90): the writer warns
            @test_logs (:warn, r"Delta-u changer with psi 0") match_mode = :any Sparlectra.writeMatpowerCasefile(netc, pathclr2; write_solution = false)
            rowc = only(r for r in eachrow(Sparlectra._matpower_sparlectra_tap_changers(Sparlectra.MatpowerIO.read_case(pathclr2))) if r[1] == 2.0)
            @test rowc[5] == 2.0 && rowc[9] == 0.0
            netc3 = importSCF(exportSCF(netc; file = joinpath(dir, "case_taps_pair_clr.scf.json")))
            bc3 = only(b for b in netc3.branchVec if b.has_ratio_tap && b.has_phase_tap)
            @test bc3.tap_ratio ≈ netc.branchVec[2].tap_ratio
            @test bc3.phase_shift_deg ≈ netc.branchVec[2].phase_shift_deg
          else
            @test_throws ArgumentError Sparlectra.writeMatpowerCasefile(netc, pathclr2; write_solution = false)
            @test_throws ArgumentError Sparlectra.net_to_scf(netc)
          end
        end
        # a 1.25 % ratio step at either band end: the band ends of the
        # fraction grid carry their own rounding noise (1.999999999999993
        # for +2), so the writer snaps them with the position, else the
        # reimport found the ratio step outside its band
        for (row, pos) in (("2 0.0125 -2 2 2 0 -10 10 -10 90 0.01;", 2.0), ("2 0.0125 -1 1 -1 0 -10 10 -10 90 0.01;", -1.0))
          pathend = joinpath(dir, "case_taps_pair_end.m")
          write(pathend, replace(read(path, String), "2 0 0 0 0 0.5 -10 10 -2;" => row))
          nete = Sparlectra.createNetFromMatPowerFile(filename = pathend)
          pathend2 = joinpath(dir, "case_taps_pair_end_rt.m")
          Sparlectra.writeMatpowerCasefile(nete, pathend2; write_solution = false)
          rowe = only(r for r in eachrow(Sparlectra._matpower_sparlectra_tap_changers(Sparlectra.MatpowerIO.read_case(pathend2))) if r[1] == 2.0)
          @test rowe[5] == pos && rowe[9] == -10.0
          nete2 = Sparlectra.createNetFromMatPowerFile(filename = pathend2)
          @test nete2.branchVec[2].tap_ratio ≈ nete.branchVec[2].tap_ratio
          @test nete2.branchVec[2].phase_shift_deg ≈ nete.branchVec[2].phase_shift_deg
        end

        # short rows are zero-padded like MATPOWER optional trailing
        # columns (the phase columns of a pure ratio changer may be
        # omitted); rejections: too many columns, unknown branch, line
        # target, and a current step outside the declared band
        shortr = joinpath(dir, "case_taps_short.m")
        write(shortr, replace(read(path, String), "1 0.025 -2 2 1 0 0 0 0;" => "1 0.025 -2 2 1;"))
        nshort = Sparlectra.createNetFromMatPowerFile(filename = shortr)
        @test nshort.branchVec[1].tap_ratio ≈ 1.0 / 1.025
        badc = joinpath(dir, "case_taps_badcols.m")
        write(badc, replace(read(path, String), "1 0.025 -2 2 1 0 0 0 0;" => "1 0.025 -2 2 1 0 0 0 0 7 7 7;"))
        @test_throws ArgumentError Sparlectra.MatpowerIO.read_case_m(badc; legacy_compat = false)
        badb = joinpath(dir, "case_taps_badbranch.m")
        write(badb, replace(read(path, String), "1 0.025 -2 2 1 0 0 0 0;" => "9 0.025 -2 2 1 0 0 0 0;"))
        @test_throws ArgumentError Sparlectra.createNetFromMatPowerFile(filename = badb)
        badl = joinpath(dir, "case_taps_online.m")
        write(badl, replace(read(path, String), "1 0.025 -2 2 1 0 0 0 0;" => "3 0.025 -2 2 1 0 0 0 0;"))
        @test_throws ArgumentError Sparlectra.createNetFromMatPowerFile(filename = badl)
        bado = joinpath(dir, "case_taps_outofband.m")
        write(bado, replace(read(path, String), "1 0.025 -2 2 1 0 0 0 0;" => "1 0.025 -2 2 5 0 0 0 0;"))
        @test_throws ArgumentError Sparlectra.createNetFromMatPowerFile(filename = bado)
        # 13 columns (a neutral ratio without its shift) and a neutral ratio
        # of 0 are rejected
        bad13 = joinpath(dir, "case_taps_13.m")
        write(bad13, replace(read(path, String), "1 0.025 -2 2 1 0 0 0 0;" => "1 0.025 -2 2 1 0 0 0 0 0 0 0 1.0;"))
        @test_throws ArgumentError Sparlectra.MatpowerIO.read_case_m(bad13; legacy_compat = false)
        bad0 = joinpath(dir, "case_taps_neutral0.m")
        write(bad0, replace(read(path, String), "1 0.025 -2 2 1 0 0 0 0;" => "1 0.025 -2 2 1 0 0 0 0 0 0 0 0 0;", "2 0 0 0 0 0.5 -10 10 -2;" => "2 0 0 0 0 0.5 -10 10 -2 0 0 0 1.0 -2.0;"))
        @test_throws ArgumentError Sparlectra.MatpowerIO.read_case_m(bad0; legacy_compat = false)

        # #513: TAP/SHIFT carry the ratio and shift in effect, so the rows
        # alone give the Y-bus Sparlectra solves; the PST of branch 2 off
        # neutral (-2 deg neutral, -3 deg in effect) and the Delta-u PST
        # (ratio and shift from the cascade) included
        for (n, p) in ((net, path2), (netdu, pathdu2))
          @test maximum(abs.(Matrix(createYBUS(net = n, sparse = true)) .- _mp_rows_ybus(p))) <= 1e-12
        end
        mp2 = Sparlectra.MatpowerIO.read_case_m(path2; legacy_compat = false)
        @test (mp2.branch[1, 9], mp2.branch[2, 10]) == (b1.tap_ratio, b2.phase_shift_deg)
        @test (mp2.sparlectra.tap_changers[1, 13], mp2.sparlectra.tap_changers[2, 14]) == (1.0, -2.0)

        # sp_case14, transformer Wehrden (branch 17): neutral 1.025, four
        # steps of 0.625 %, 1.0 in effect. The export writes 1.0 as TAP,
        # the neutral and the position into the block
        sp14 = importSCF(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"))
        bw = sp14.branchVec[17]
        @test (bw.ratio, bw.tap_ratio) == (1.025, 1.0)
        p14 = joinpath(dir, "sp_case14.m")
        Sparlectra.writeMatpowerCasefile(sp14, p14; write_solution = false)
        Y14 = Matrix(createYBUS(net = sp14, sparse = true))
        @test maximum(abs.(Y14 .- _mp_rows_ybus(p14))) <= 1e-12
        mp14 = Sparlectra.MatpowerIO.read_case_m(p14; legacy_compat = false)
        tc14 = mp14.sparlectra.tap_changers
        w = findfirst(==(17.0), tc14[:, 1])
        @test mp14.branch[17, 9] == 1.0
        @test (tc14[w, 13], round(tc14[w, 5]; digits = 9)) == (1.025, 4.0)
        # the Sparlectra re-import: same neutral, position, band and grid on
        # both readers, the same Y-bus bitwise; the power flow agrees to
        # 1e-12 (not bitwise: the MATPOWER import starts from other values)
        cfg14 = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
        back14 = Sparlectra.createNetFromMatPowerFile(filename = p14)
        conv14 = Sparlectra.build_net(Sparlectra.convert_case(MatpowerAdapter(), Sparlectra.MatpowerIO.read_case(p14), Sparlectra.matpower_adapter_options(cfg14)); config = cfg14)
        # transformer tap fields (the line rows of the SCF and MATPOWER
        # imports differ in fields a line never reads)
        tapfields(n) = [(b.ratio, b.angle, b.tap_ratio, b.phase_shift_deg, b.tap_min, b.tap_max, b.tap_step, b.tap_ratio_grid, b.has_ratio_tap, b.has_phase_tap, b.phase_min_deg, b.phase_max_deg) for b in n.branchVec if b.ratio != 0.0]
        for back in (back14, conv14)
          @test tapfields(back) == tapfields(sp14)
          @test Matrix(createYBUS(net = back, sparse = true)) == Y14
        end
        pf14 = deepcopy(sp14)
        pfb = deepcopy(back14)
        @test runpf!(pf14, 40, 1e-10, 0)[2] == 0
        @test runpf!(pfb, 40, 1e-10, 0)[2] == 0
        @test maximum(abs(a._vm_pu - b._vm_pu) for (a, b) in zip(pf14.nodeVec, pfb.nodeVec)) <= 1e-12
        @test maximum(abs(a._va_deg - b._va_deg) for (a, b) in zip(pf14.nodeVec, pfb.nodeVec)) <= 1e-12

        # TAP edited by hand in a new file: the tap position is kept and the
        # difference is named
        p14ed = joinpath(dir, "sp_case14_edited.m")
        txt14 = read(p14, String)
        txt_ed = replace(txt14, "\t25.0\t0\t0\t1.0\t0.0\t1\t" => "\t25.0\t0\t0\t1.01\t0.0\t1\t")
        @test count("1.01", txt_ed) == count("1.01", txt14) + 1
        write(p14ed, txt_ed)
        ed14 = @test_logs (:warn, r"differ from the ratio and shift of the tap position") match_mode = :any Sparlectra.createNetFromMatPowerFile(filename = p14ed)
        @test ed14.branchVec[17].tap_ratio == 1.0
      end
    end)() end

    # one active row (status 1 in column 3) counts as active
    @test Sparlectra.MatpowerIO.matpower_dcline_diagnostics((; dcline = [1 2 1 10 9 0 0 1 1 0 100]))["matpower_dcline_active_count"] == 1
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

    @testset "MATPOWER export write_solution" begin (function ()
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
        # a unit without an active-power limit is written with the MATPOWER
        # "unlimited" convention, not with +-baseMVA
        gen_row = findfirst(==(1.0), mpc_false.gen[:, 1])
        @test (mpc_false.gen[gen_row, 9], mpc_false.gen[gen_row, 10]) == (9999.0, -9999.0)

        # a hand-built net exported before any solve: its node voltages are
        # still 1.0, the setpoints (1.03 for every CIGRE unit) live on the
        # generators; Vg and the PV and slack bus Vm must carry them
        cigre = createCIGRE()
        path_cigre = joinpath(dir, "cigre_unsolved.m")
        Sparlectra.writeMatpowerCasefile(cigre, path_cigre; write_solution = false)
        mpc_cigre = Sparlectra.MatpowerIO.read_case_m(path_cigre; legacy_compat = false)
        @test size(mpc_cigre.gen, 1) == 4
        @test all(==(1.03), mpc_cigre.gen[:, 6])
        regulated = [r for r in axes(mpc_cigre.bus, 1) if mpc_cigre.bus[r, 2] in (2.0, 3.0)]
        @test length(regulated) == 4
        @test all(==(1.03), mpc_cigre.bus[regulated, 8])
        # the round trip keeps the result: the reimported net solves to the
        # same voltages as the net it was written from
        cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 50, tol = 1e-10, qlimits = Sparlectra.QLimitConfig(ignore_q_limits = true)), output = OutputConfig(logfile_results = :off, startup_latency_hint = false))
        back = Sparlectra.createNetFromMatPowerFile(filename = path_cigre)
        @test run_sparlectra(net = cigre, config = cfg).numerical_converged
        @test run_sparlectra(net = back, config = cfg).numerical_converged
        @test maximum(abs(a._vm_pu - b._vm_pu) for (a, b) in zip(cigre.nodeVec, back.nodeVec)) < 1e-8
        @test maximum(abs(a._va_deg - b._va_deg) for (a, b) in zip(cigre.nodeVec, back.nodeVec)) < 1e-6
      end

      # three-winding star points through a model-only round trip: two star
      # equivalents whose 110 kV windings (r 0.00044 pu, ratio 0.9649) close
      # a loop. At 1.0 pu the star points drive a large circulating current
      # and Newton does not converge; the solver's flat start puts them at
      # the ratio of that winding. Before the fix the file wrote 1.0 pu and
      # the reader built plain buses, so the re-read net failed both from the
      # file voltages and from a flat start (50 iterations, status 1, as the
      # CGMES MiniGrid did).
      n3 = Sparlectra.Net(name = "star_roundtrip", baseMVA = 100.0)
      for (b, vn, aux) in (("S10", 10.0, false), ("M110", 110.0, false), ("H380", 380.0, false), ("X1", 380.0, true), ("X2", 380.0, true), ("T1", 30.0, false), ("T2", 30.0, false), ("L10", 10.0, false), ("N110", 110.0, false))
        Sparlectra.addBus!(net = n3, busName = b, vn_kV = vn, isAux = aux)
      end
      Sparlectra.addProsumer!(net = n3, busName = "S10", type = "SYNCHRONOUSMACHINE", p = 0.0, q = 0.0, vm_pu = 1.0, va_deg = 0.0, referencePri = "S10")
      Sparlectra.addProsumer!(net = n3, busName = "L10", type = "ENERGYCONSUMER", p = 9.0, q = 5.0)
      Sparlectra.addPIModelTrafo!(net = n3, fromBus = "M110", toBus = "S10", r_pu = 0.0055, x_pu = 0.132, b_pu = 0.0, status = 1, ratio = 1.039)
      Sparlectra.addPIModelTrafo!(net = n3, fromBus = "M110", toBus = "L10", r_pu = 0.0175, x_pu = 0.42, b_pu = 0.0, status = 1, ratio = 0.9957)
      for (x, t) in (("X1", "T1"), ("X2", "T2"))
        Sparlectra.addPIModelTrafo!(net = n3, fromBus = "H380", toBus = x, r_pu = 0.00041, x_pu = 0.0665, b_pu = 0.0, status = 1, ratio = 1.0)
        Sparlectra.addPIModelTrafo!(net = n3, fromBus = x, toBus = "M110", r_pu = 0.00044, x_pu = 1.0e-5, b_pu = 0.0, status = 1, ratio = 0.9649)
        Sparlectra.addPIModelTrafo!(net = n3, fromBus = x, toBus = t, r_pu = 0.0028, x_pu = 0.14, b_pu = 0.0, status = 1, ratio = 1.0526)
      end
      # a line beside the star: the exporter must tell it from the ratio
      # 1.0 legs (a native line carries ratio 0, no metadata on either)
      Sparlectra.addPIModelACLine!(net = n3, fromBus = "M110", toBus = "N110", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
      # since 0.50.0 addBus! gives the node its busName, and an isAux bus is
      # auxiliary by its component type whatever it is called (X1, X2 carry
      # no Aux_ prefix)
      @test [Sparlectra.getCompName(n.comp) for n in n3.nodeVec] == ["S10", "M110", "H380", "X1", "X2", "T1", "T2", "L10", "N110"]
      aux3 = [Sparlectra._is_aux_node(n) for n in n3.nodeVec]
      @test aux3 == [false, false, false, true, true, false, false, false, false]
      @test Sparlectra.runpf!(n3, 50, 1e-8, 0; opt_flatstart = true) == (5, 0)
      # the re-read net from the file voltages and from a flat start
      solves_like_n3(net) = for flat in (false, true)
        nb = deepcopy(net)
        @test Sparlectra.runpf!(nb, 50, 1e-8, 0; opt_flatstart = flat) == (5, 0)
        @test maximum(abs(nb.nodeVec[i]._vm_pu - n3.nodeVec[i]._vm_pu) for i in eachindex(n3.nodeVec)) < 1e-8
      end
      cfg3 = Sparlectra.active_sparlectra_config()
      mktempdir() do dir
        path3 = joinpath(dir, "star.m")
        Sparlectra.writeMatpowerCasefile(n3, path3; write_solution = false)
        stars = [n3.busDict["X1"], n3.busDict["X2"]]
        mpc3 = Sparlectra.MatpowerIO.read_case_m(path3; legacy_compat = false)
        @test all(mpc3.bus[stars, 8] .≈ 0.9649)
        # the role travels as the bus list mpc.sparlectra.aux_buses
        @test vec(mpc3.sparlectra.aux_buses) == stars
        # mpc.branch_kind and TAP row by row: the eight transformers "T"
        # with their ratio (the 380 kV legs at 1.0), the line "L" with TAP
        # 0. Judging by ratio != 1 wrote the line as "T" and the legs as
        # lines with TAP 0, which the re-read then built as lines.
        @test mpc3.branch_kind == vcat(fill("T", 8), ["L"])
        @test mpc3.branch[:, 9] == [1.039, 0.9957, 1.0, 0.9649, 1.0526, 1.0, 0.9649, 1.0526, 0.0]
        back3 = Sparlectra.createNetFromMatPowerFile(filename = path3)
        @test (length(back3.trafos), length(back3.linesAC)) == (8, 1)
        @test [Sparlectra._is_aux_node(n) for n in back3.nodeVec] == aux3
        solves_like_n3(back3)
        # the SCF converter path reads the same role
        conv3 = Sparlectra.build_net(Sparlectra.convert_case(Sparlectra.MatpowerAdapter(), Sparlectra.MatpowerIO.read_case(path3; legacy_compat = true), Sparlectra.matpower_adapter_options(cfg3)); config = cfg3)
        @test [Sparlectra._is_aux_node(n) for n in conv3.nodeVec] == aux3
        # an SCF round trip keeps the role (roles.aux_nodes)
        scf3 = Sparlectra.scf_to_net(Sparlectra.net_to_scf(n3))
        @test [Sparlectra._is_aux_node(n) for n in scf3.nodeVec] == aux3
      end
    end)() end

    @testset "Tap-changer-model marker prevents double correction on reimport" begin (function ()
      tap_bus = [1 3 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9; 2 1 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9]
      tap_gen = [1 10.0 0.0 100.0 -100.0 1.0 100.0 1 100.0 0.0 0 0 0 0 0 0 0 0 0 0 0]
      tap_branch = [1 2 0.01 0.05 0.0 100.0 0.0 0.0 0.9 0.0 1 -60.0 60.0]
      mpc_tap = Sparlectra.MatpowerIO.MatpowerCase("tap_marker_case", 100.0, tap_bus, tap_gen, tap_branch, nothing, nothing)

      # a case without the marker: :ideal keeps x, the correction changes it
      @test mpc_tap.sparlectra === nothing
      @test Sparlectra.createNetFromMatPowerCase(mpc = mpc_tap, tap_changer_model = :ideal).branchVec[1].x_pu ≈ 0.05

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
    end)() end

    @testset "one reference per synchronous area (#465)" begin (function ()
      # each island keeps the type-3 bus of the file; a second type-3 bus in
      # the same AC island stays a voltage-regulated PV bus, with a warning
      dir = mktempdir()
      two_refs = """
function mpc = case_b2b
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
1 3 0 0 0 0 1 1.0 0 380 1 1.1 0.9;
2 1 40 10 0 0 1 1.0 0 380 1 1.1 0.9;
3 3 0 0 0 0 1 1.0 0 380 1 1.1 0.9;
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
"""
      path_refs = joinpath(dir, "case_b2b_two_refs.m")
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
      path_tied = joinpath(dir, "case_b2b_tied.m")
      write(path_tied, replace(two_refs, "3 4 0.01 0.08 0 250 250 250 0 0 1 -360 360;" => "3 4 0.01 0.08 0 250 250 250 0 0 1 -360 360;\n2 3 0.01 0.08 0 250 250 250 0 0 1 -360 360;"))
      net_tied = @test_logs (:warn, r"second reference bus 3 in the AC island of reference bus 1") match_mode = :any createNetFromMatPowerFile(filename = path_tied, matpower_dcline_mode = :pf_injections)
      @test Sparlectra.getNodeType(net_tied.nodeVec[1]) == Sparlectra.Slack
      @test Sparlectra.getNodeType(net_tied.nodeVec[3]) == Sparlectra.PV
      # type 3 is reference priority 1 for every unit of the bus, the demoted
      # second reference included; the SCF conversion writes the same
      @test [ps.referencePriority for ps in net_tied.prosumpsVec if Sparlectra.getPosumerBusIndex(ps) in (1, 3)] == [1, 1]
      @test count(e -> e isa AbstractDict && get(e, "reference_priority", 0) == 1, values(scf_refs.sparlectra.extra)) == 2
    end)() end
  end)() end
end

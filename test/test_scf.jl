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

# file: test/test_scf.jl
# purpose: Sparlectra Case Format writer (issue #342, stage 1): deterministic
#          ids and byte-identical repeat exports, SI unit conversion against
#          hand-computed references, PGM subset shape, the namespaced
#          sparlectra block (roles, tap cascade, extra names, measurements),
#          the FACTS base-impedance guard, and the Web UI export route

# the shared capture helper needs both stdlibs when this file runs standalone
using Logging
# the Web UI export test packs a shipped CGMES delivery into a zip, the form
# the case directory holds it in
import ZipArchives

# the fixture comparison lives with the demo-case tests, so both use one rule
isdefined(@__MODULE__, :scf_roundtrip_field_diffs) || include("test_scf_support.jl")
# the synthetic CGMES deliveries with tap changers (estimation round trip)
isdefined(@__MODULE__, :_cgmes_synth_tapside_dir) || include(joinpath(@__DIR__, "cgmes_synth3w_fixture.jl"))

# warmup_casePST.m is TRACKED and lives in the checkout; anything else has to
# come from SPARLECTRA_LARGE_CASES_DIR, so a caller passing another name gates
# on large_case_path first
function _scf_test_net(name::AbstractString = "warmup_casePST.m")
  cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
  path = name == "warmup_casePST.m" ? abspath(joinpath(dirname(@__DIR__), "data", "mpower", name)) : large_case_path(name)
  path === nothing && error(string("_scf_test_net: ", name, " is not in the large-case directory; gate on large_case_path before calling"))
  return SparlectraApp._se_import_case_net(path, cfg)
end





function run_scf_tests()
  @testset "units declaration" begin (function ()
    # absent means si: every existing file reads unchanged
    plain = Sparlectra.importSCF(joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json"))
    case_plain = Sparlectra.read_scf_json(joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json"))
    @test Sparlectra.scf_case_units(case_plain) === :si
    @test !haskey(case_plain.sparlectra.meta, "units")
    # explicit si is accepted
    case_plain.sparlectra.meta["units"] = "si"
    @test Sparlectra.scf_case_units(case_plain) === :si
    # pu with a complete base is accepted
    case_plain.sparlectra.meta["units"] = "pu"
    @test Sparlectra.scf_case_units(case_plain) === :pu
    # junk is a hard error, not a silent si
    case_plain.sparlectra.meta["units"] = "furlongs"
    @test_throws ArgumentError Sparlectra.scf_case_units(case_plain)
    case_plain.sparlectra.meta["units"] = 42
    @test_throws ArgumentError Sparlectra.scf_case_units(case_plain)
    # pu without its base names the missing field
    case_plain.sparlectra.meta["units"] = "pu"
    s_base_kept = pop!(case_plain.sparlectra.meta, "s_base")
    err = try
      Sparlectra.scf_case_units(case_plain)
      nothing
    catch e
      e
    end
    @test err isa ArgumentError
    @test occursin("s_base", sprint(showerror, err))
    case_plain.sparlectra.meta["s_base"] = s_base_kept
    # pu with a broken node base names the node
    u_kept = case_plain.data.node[1].u_rated
    case_plain.data.node[1].u_rated = 0.0
    err2 = try
      Sparlectra.scf_case_units(case_plain)
      nothing
    catch e
      e
    end
    @test err2 isa ArgumentError
    @test occursin("u_rated", sprint(showerror, err2))
    case_plain.data.node[1].u_rated = u_kept
    @test plain isa Sparlectra.Net

    # units = pu: build_net reads the pu block as it stands. The case is
    # annotated with pu values formed from its own SI view by the reader's
    # SI formula, so both builds must agree bit for bit (the SI regression),
    # and one planted value the SI path cannot express must arrive unchanged
    src = joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json")
    case_si = Sparlectra.read_scf_json(src)
    case_pu = Sparlectra.read_scf_json(src)
    s_base = Sparlectra._scf_num(case_pu.sparlectra.meta["s_base"], "s_base")
    u_by = Dict(n.id => n.u_rated for n in case_pu.data.node)
    zb(id) = u_by[id]^2 / s_base
    w_nom = 2.0 * pi * Sparlectra._scf_num(get(case_pu.sparlectra.meta, "f_nom", 50.0), "f_nom")
    magic = 0.015404508911002418
    magic_id = nothing
    for row in case_pu.data.line
      z = zb(row.to_node)
      e = get!(case_pu.sparlectra.extra, string(row.id), Dict{String,Any}())
      b_s = row.c1 === nothing ? 0.0 : row.c1 * w_nom
      tan1 = row.tan1 === nothing ? 0.0 : row.tan1
      pu = Dict{String,Any}("r_pu" => row.r1 / z, "x_pu" => row.x1 / z, "b_pu" => b_s * z, "g_pu" => row.g1 === nothing ? tan1 * b_s * z : row.g1 * z)
      if magic_id === nothing
        pu["x_pu"] = magic
        magic_id = row.id
      end
      e["pu"] = pu
    end
    for row in case_pu.data.generic_branch
      z = zb(row.to_node)
      e = get!(case_pu.sparlectra.extra, string(row.id), Dict{String,Any}())
      e["pu"] = Dict{String,Any}("r_pu" => row.r1 / z, "x_pu" => row.x1 / z, "b_pu" => row.b1 === nothing ? 0.0 : row.b1 * z, "g_pu" => row.g1 === nothing ? 0.0 : row.g1 * z)
    end
    for row in case_pu.data.shunt
      z = zb(row.node)
      e = get!(case_pu.sparlectra.extra, string(row.id), Dict{String,Any}())
      e["pu"] = Dict{String,Any}("g_pu" => something(row.g1, 0.0) * z, "b_pu" => something(row.b1, 0.0) * z)
    end
    case_pu.sparlectra.meta["units"] = "pu"
    net_si = Sparlectra.build_net(case_si)
    net_pu = Sparlectra.build_net(case_pu)
    @test length(net_si.branchVec) == length(net_pu.branchVec)
    @test length(net_si.shuntVec) == length(net_pu.shuntVec)
    planted = 0
    for i in eachindex(net_si.branchVec)
      a = net_si.branchVec[i]
      b = net_pu.branchVec[i]
      if b.x_pu === magic && a.x_pu !== magic
        planted += 1
      else
        @test a.x_pu === b.x_pu
      end
      @test a.r_pu === b.r_pu
      @test a.b_pu === b.b_pu
      @test a.g_pu === b.g_pu
    end
    @test planted == 1
    for i in eachindex(net_si.shuntVec)
      @test net_si.shuntVec[i].y_pu_shunt === net_pu.shuntVec[i].y_pu_shunt
    end
    # a pu file whose branch lacks its pu block is refused, naming the id
    broken = Sparlectra.read_scf_json(src)
    broken.sparlectra.meta["units"] = "pu"
    err3 = try
      Sparlectra.build_net(broken)
      nothing
    catch e
      e
    end
    @test err3 isa ArgumentError
    @test occursin("pu block", sprint(showerror, err3))

    # MATPOWER states its parameters per unit. Both ways into a case file,
    # the adapter (convert_case) and the export of the directly imported
    # net with units = :pu, must read back to the direct import's values
    # bit for bit; through SI, sp_case118 loses the last bit on 123 of its
    # 744 branch values (no exact preimage)
    mfile = joinpath(dirname(@__DIR__), "data", "mpower", "sp_case118.m")
    cfg = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, mfile).config
    mpc = Sparlectra.MatpowerIO.read_case(mfile; legacy_compat = true)
    direct = Sparlectra.import_net(MatpowerAdapter(), mpc, Sparlectra.matpower_adapter_options(cfg))
    d = mktempdir()
    adapter_file = Sparlectra.write_scf_json(Sparlectra.convert_case(MatpowerAdapter(), mpc, Sparlectra.matpower_adapter_options(cfg)), joinpath(d, "adapter.scf.json"))
    export_file = Sparlectra.exportSCF(direct; file = joinpath(d, "export.scf.json"), units = :pu)
    si_file = Sparlectra.exportSCF(direct; file = joinpath(d, "export_si.scf.json"))
    # the default export stays SI and declares nothing (shipped files unchanged)
    @test !haskey(Sparlectra.read_scf_json(si_file).sparlectra.meta, "units")
    for f in (adapter_file, export_file)
      case_m = Sparlectra.read_scf_json(f)
      @test Sparlectra.scf_case_units(case_m) === :pu
      back = Sparlectra.build_net(case_m; config = cfg)
      @test length(back.branchVec) == length(direct.branchVec)
      @test count(i -> any(fld -> !isequal(getfield(direct.branchVec[i], fld), getfield(back.branchVec[i], fld)), (:r_pu, :x_pu, :b_pu, :g_pu)), eachindex(direct.branchVec)) == 0
      @test length(back.shuntVec) == length(direct.shuntVec)
      @test all(i -> isequal(direct.shuntVec[i].y_pu_shunt, back.shuntVec[i].y_pu_shunt), eachindex(direct.shuntVec))
    end
    # the SI export of the same net is the case the pu block exists for
    back_si = Sparlectra.build_net(Sparlectra.read_scf_json(si_file); config = cfg)
    @test any(i -> !isequal(direct.branchVec[i].x_pu, back_si.branchVec[i].x_pu) || !isequal(direct.branchVec[i].r_pu, back_si.branchVec[i].r_pu) || !isequal(direct.branchVec[i].b_pu, back_si.branchVec[i].b_pu), eachindex(direct.branchVec))
    # the pu values live in the sparlectra block, which strict PGM drops
    @test_throws ArgumentError Sparlectra.exportSCF(direct; file = joinpath(d, "strict.pgm.json"), units = :pu, strict_pgm = true)
  end)() end

  @testset "own bus names vs external_id" begin (function ()
    # the shipped demo cases use the split deliberately: extra.name is the
    # invented place name (the reference name), external_id the short
    # technical handle. The import must apply the place name everywhere
    # AND restore the handle, so a re-export preserves both fields.
    net = Sparlectra.importSCF(joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json"))
    @test haskey(net.busDict, "Sandau_110")
    @test haskey(net.busDict, "Lindhof_20")
    slack = net.nodeVec[net.busDict["Sandau_110"]]
    @test slack.comp.cID == "HV1"
    @test Sparlectra.getCompName(slack.comp) == "Sandau_110"
    d = mktempdir()
    f = exportSCF(net; file = joinpath(d, "names.scf.json"))
    root = Sparlectra.scf_json_parse(read(f, String))
    node_extra = [v for v in values(root["sparlectra"]["extra"]) if haskey(v, "bus_index")]
    by_handle = Dict(String(v["external_id"]) => String(v["name"]) for v in node_extra)
    @test by_handle["HV1"] == "Sandau_110"
    @test by_handle["MV1"] == "Lindhof_20"
    # no component_name rows appear: reference name and component name
    # agree on purpose in the shipped files
    @test all(!haskey(v, "component_name") for v in node_extra)

    # branches, appliances and links keep the name of the file as well, not
    # the generated one (B_2WT_110_4_5 before 0.60.0), and the external_id
    # as their component id; the kind of a named branch still reads from its
    # build (Trafo, Line), and the re-export writes the names back
    tr = net.branchVec[5]
    @test tr.comp.cName == "Lindhof_110_20" && tr.comp.cID == "#B_2WT_110_4_5#5"
    @test Sparlectra._branch_kind_label(net, tr) == "Trafo"
    @test net.branchVec[1].comp.cName == "Sandau_Borntal" && Sparlectra._branch_kind_label(net, net.branchVec[1]) == "Line"
    br_extra = Dict(Int(v["branch_index"]) => String(v["name"]) for v in values(root["sparlectra"]["extra"]) if haskey(v, "branch_index"))
    @test br_extra[5] == "Lindhof_110_20" && br_extra[1] == "Sandau_Borntal"
    f6 = joinpath(dirname(@__DIR__), "data", "scf", "sp_case6tap.scf.json")
    net6 = Sparlectra.importSCF(f6)
    @test [b.comp.cName for b in net6.branchVec[5:7]] == ["T1_Nordwerk_Kreuzau", "T2_Suedfeld_Ostheim", "T3_Mittelbach_Lindenau"]
    @test [p.comp.cName for p in net6.prosumpsVec] == ["Gen_Nordwerk", "Gen_Suedfeld", "Ld_Mittelbach", "Ld_Ostheim", "Ld_Lindenau"]

    # a generic_branch without a tap_changer row has no changer: the
    # constructor's default grid (0.9 to 1.1, read as positions -14..17)
    # was kept before 0.60.0. T2 of sp_case6tap declares one (+-12 steps of
    # 1.25 %); with its row taken out of the file it has none
    t2 = net6.branchVec[6]
    @test t2.has_ratio_tap && !t2.has_phase_tap && t2.tap_step == 0.0125
    @test t2.tap_min ≈ 1 / 1.15 && t2.tap_max ≈ 1 / 0.85
    root6 = Sparlectra.scf_json_parse(read(f6, String))
    filter!(r -> r["branch"] != 12, root6["sparlectra"]["components"]["tap_changer"])
    bare = Sparlectra.scf_to_net(root6)
    @test !bare.branchVec[6].has_ratio_tap && !bare.branchVec[6].has_phase_tap
    @test bare.branchVec[7].has_ratio_tap && !bare.branchVec[7].has_phase_tap
  end)() end

  @testset "deterministic ids and byte-identical repeat export" begin (function ()
    net = _scf_test_net()
    d = mktempdir()
    a = exportSCF(net; file = joinpath(d, "a.scf.json"), source_format = "matpower")
    b = exportSCF(net; file = joinpath(d, "b.scf.json"), source_format = "matpower")
    @test read(a, String) == read(b, String)
    # a freshly imported net of the same case must produce the same bytes:
    # ids come from the component names, not from Dict iteration order
    net2 = _scf_test_net()
    c = exportSCF(net2; file = joinpath(d, "c.scf.json"), source_format = "matpower")
    @test read(a, String) == read(c, String)
    # the file ends with a newline and parses as the documented root shape
    txt = read(a, String)
    @test endswith(txt, "\n")
    @test startswith(txt, "{\n  \"version\": \"1.0\",\n  \"type\": \"input\",\n  \"is_batch\": false,")
  end)() end

  @testset "PGM subset shape and SI unit conversion" begin (function ()
    net = _scf_test_net()
    root = Sparlectra.net_to_scf(net; f_nom = 50.0)
    data = root["data"]
    @test root["version"] == "1.0"
    @test root["type"] == "input"
    @test root["is_batch"] === false
    @test haskey(data, "node") && haskey(data, "line") && haskey(data, "generic_branch")
    # only components of the documented Rev. 1 subset appear
    allowed = Set(["node", "line", "generic_branch", "link", "source", "sym_load", "sym_gen", "shunt", "voltage_regulator", "sym_voltage_sensor", "sym_power_sensor", "sym_current_sensor", "fault"])
    @test all(k -> k in allowed, keys(data))
    # ids are unique across every component of the file, including the
    # sparlectra-only ones
    ids = Int[]
    for (_, rows) in data
      append!(ids, Int[r["id"] for r in rows])
    end
    for (_, rows) in get(root["sparlectra"], "components", Dict{String,Any}())
      append!(ids, Int[r["id"] for r in rows])
    end
    @test length(ids) == length(unique(ids))
    # node voltages in V
    nodes = Dict(r["id"] => r for r in data["node"])
    names = root["sparlectra"]["extra"]
    bus1 = parse(Int, only(id for (id, e) in names if get(e, "bus_index", 0) == 1))
    @test nodes[bus1]["u_rated"] ≈ Sparlectra.getNodeVn(net.nodeVec[1]) * 1000.0
    # line impedance in ohm on the TO-side base (Sparlectra stamps the tap
    # on the from side, so r_pu/x_pu are to-side referenced like PGM's
    # generic_branch); check one line against the hand computation
    lrow = first(data["line"])
    lidx = only(get(e, "branch_index", 0) for (id, e) in names if string(lrow["id"]) == id)
    br = net.branchVec[lidx]
    vn = Sparlectra.getNodeVn(net.nodeVec[Int(br.toBus)])
    zb = (vn * 1000.0)^2 / (net.baseMVA * 1.0e6)
    @test lrow["r1"] ≈ br.r_base_pu * zb
    @test lrow["x1"] ≈ br.x_base_pu * zb
    # charging susceptance becomes a capacitance in farad
    @test lrow["c1"] ≈ (br.b_pu / zb) / (2 * pi * 50.0)
    # transformer: k and theta describe the LIVE tap position
    grow = first(data["generic_branch"])
    gidx = only(get(e, "branch_index", 0) for (id, e) in names if string(grow["id"]) == id)
    t = Sparlectra.calcBranchRatio(net.branchVec[gidx])
    @test grow["k"] ≈ abs(t)
    @test grow["theta"] ≈ angle(t)
    # appliance powers in SI: the internal MW/MVar become W/var
    if haskey(data, "sym_load")
      ids = Sparlectra.scf_id_map(net)
      byid = Dict(v => k for (k, v) in ids.prosumer)
      lo = first(data["sym_load"])
      ps = net.prosumpsVec[byid[lo["id"]]]
      @test lo["p_specified"] ≈ (ps.pVal === nothing ? 0.0 : ps.pVal) * 1.0e6
      @test lo["q_specified"] ≈ (ps.qVal === nothing ? 0.0 : ps.qVal) * 1.0e6
    end
  end)() end

  @testset "sparlectra block: roles, taps, names" begin (function ()
    net = _scf_test_net()
    root = Sparlectra.net_to_scf(net; intended_calculations = String["power_flow", "state_estimation"])
    sp = root["sparlectra"]
    @test sp["format_version"] == Sparlectra.SCF_FORMAT_VERSION
    @test sp["meta"]["s_base"] ≈ net.baseMVA * 1.0e6
    @test sp["meta"]["f_nom"] == 50.0
    @test "state_estimation" in sp["meta"]["intended_calculations"]
    # the slack is stated explicitly, never inferred by the reader
    @test sp["roles"]["slack"]["mode"] == "single"
    @test length(sp["roles"]["slack"]["nodes"]) == 1
    # every component carries its reference name, the resolution table the
    # controller and contingency sections rely on
    @test all(e -> haskey(e, "name"), values(sp["extra"]))
    # the tap cascade: the demo case has a ratio changer and a Delta-u PST
    taps = sp["components"]["tap_changer"]
    @test length(taps) == 2
    ratio_tc = only(t for t in taps if any(c -> c["index"] == 1, t["controllers"]))
    @test only(c for c in ratio_tc["controllers"] if c["index"] == 1)["step"] > 0.0
    pst_tc = only(t for t in taps if any(c -> c["index"] == 2, t["controllers"]))
    pstc = only(c for c in pst_tc["controllers"] if c["index"] == 2)
    @test pstc["alpha_deg"] == 90.0
    @test pstc["step"] == 0.01
    @test pstc["pos_min"] == -10 && pstc["pos_max"] == 10
    # start state is opt-in and never called a result
    @test !haskey(sp, "start_state")
    root2 = Sparlectra.net_to_scf(net; include_start_state = true)
    @test haskey(root2["sparlectra"], "start_state")
    @test root2["sparlectra"]["start_state"]["source"] == "solved_power_flow"
    @test !haskey(root2["sparlectra"], "results")
  end)() end

  @testset "three-winding grouping and reference names" begin (function ()
    # the star equivalent as Sparlectra builds it (and as PGM documents
    # it): three generic_branch legs plus one auxiliary star node
    net = Net(name = "t3w", baseMVA = 100.0)
    addBus!(net = net, busName = "HV", vn_kV = 380.0)
    addBus!(net = net, busName = "MV", vn_kV = 110.0)
    addBus!(net = net, busName = "LV", vn_kV = 30.0)
    addBus!(net = net, busName = "AUX3WT_T1", vn_kV = 380.0, isAux = true)
    for (term, r, x) in (("HV", 0.001, 0.05), ("MV", 0.002, 0.04), ("LV", 0.003, 0.03))
      addPIModelTrafo!(net = net, fromBus = "AUX3WT_T1", toBus = term, r_pu = r, x_pu = x, b_pu = 0.0, status = 1, ratio = 1.0, ratedS = 300.0)
    end
    addProsumer!(net = net, busName = "HV", type = "ExternalNetworkInjection", referencePri = "HV", vm_pu = 1.0, va_deg = 0.0)
    addProsumer!(net = net, busName = "MV", type = "ENERGYCONSUMER", p = 50.0, q = 10.0)
    root = Sparlectra.net_to_scf(net)
    sp = root["sparlectra"]
    groups = sp["components"]["transformer3w"]
    @test length(groups) == 1
    g = only(groups)
    # the star node is declared auxiliary, so reports and observability
    # statistics can exclude it
    @test g["star_node"] in sp["roles"]["aux_nodes"]
    @test length(g["ends"]) == 3
    @test [e["role"] for e in g["ends"]] == ["hv", "mv", "lv"]
    # roles follow the terminal voltage, highest first
    @test [e["u_rated"] for e in g["ends"]] == [380000.0, 110000.0, 30000.0]
    # the electrical values live exactly once, in the three legs
    legs = Set(e["branch"] for e in g["ends"])
    @test length(legs) == 3
    @test all(r -> r["id"] in legs, root["data"]["generic_branch"])
    @test g["leg_direction"] == "star_to_terminal"
    # every leg starts at the star node, so no reader has to guess
    for row in root["data"]["generic_branch"]
      @test row["from_node"] == g["star_node"]
    end
    # a leg turned to terminal → star (a tapped CGMES leg) is reported
    turned = deepcopy(net)
    br = turned.branchVec[2]
    br.fromBus, br.toBus = br.toBus, br.fromBus
    @test only(Sparlectra.net_to_scf(turned)["sparlectra"]["components"]["transformer3w"])["leg_direction"] == "mixed"
    # reference names (the busDict keys every other section resolves
    # against) are what `extra.name` carries; addBus! names the node after
    # its busName, so no separate component name travels
    byname = Dict(e["name"] => e for e in values(sp["extra"]))
    @test haskey(byname, "AUX3WT_T1")
    @test !haskey(byname["AUX3WT_T1"], "component_name")
    @test haskey(byname, "HV") && haskey(byname, "MV") && haskey(byname, "LV")
  end)() end

  @testset "round-trip identity: bytes, power flow, state estimation" begin (function ()
    # the format's core promise: build, write, read, write produces the
    # same bytes, and a run on the re-read network reproduces the original
    # numerically, not just approximately
    net = _scf_test_net()
    Sparlectra.readMeasurementTable!(net; file = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.measurements.csv")), zero_injection = false)
    d = mktempdir()
    kw = (source_format = "matpower", source_reference = "warmup_casePST.m", include_start_state = true)
    a = exportSCF(net; file = joinpath(d, "a.scf.json"), kw...)
    back = importSCF(a)
    b = exportSCF(back; file = joinpath(d, "b.scf.json"), kw...)
    @test read(a, String) == read(b, String)
    # structure survives element by element
    @test length(back.nodeVec) == length(net.nodeVec)
    @test length(back.branchVec) == length(net.branchVec)
    @test length(back.prosumpsVec) == length(net.prosumpsVec)
    @test length(back.shuntVec) == length(net.shuntVec)
    @test length(back.linkVec) == length(net.linkVec)
    @test length(back.measurements) == length(net.measurements)
    @test Set(keys(back.busDict)) == Set(keys(net.busDict))
    # the model the solver sees is bit-identical
    @test Sparlectra.createYBUS(net = net, sparse = true) == Sparlectra.createYBUS(net = back, sparse = true)
    @test [b_.tap_ratio for b_ in back.branchVec] == [b_.tap_ratio for b_ in net.branchVec]
    @test [b_.phase_shift_deg for b_ in back.branchVec] == [b_.phase_shift_deg for b_ in net.branchVec]
    @test [b_.has_phase_tap for b_ in back.branchVec] == [b_.has_phase_tap for b_ in net.branchVec]
    @test [b_.phase_du_step for b_ in back.branchVec] == [b_.phase_du_step for b_ in net.branchVec]
    # measurements survive with value, sigma, order, and id
    @test [m.id for m in back.measurements] == [m.id for m in net.measurements]
    @test [m.value for m in back.measurements] == [m.value for m in net.measurements]
    @test [m.sigma for m in back.measurements] == [m.sigma for m in net.measurements]
    # the measurements as PGM sensors: P and Q of one location pair into a
    # single power sensor, so there are fewer sensors than measurements, and
    # every sensor is named by a row
    adoc = Sparlectra.scf_json_parse(read(a, String))
    adata = adoc["data"]
    arows = adoc["sparlectra"]["measurements"]["rows"]
    sensor_kinds = ("sym_voltage_sensor", "sym_power_sensor", "sym_current_sensor")
    @test sum(length(get(adata, k, [])) for k in sensor_kinds) < length(net.measurements)
    @test Set(r["sensor"] for r in arows) == Set(s["id"] for k in sensor_kinds for s in get(adata, k, []))
    # units without a base: U and its sigma in kV of the node, voltage
    # sensors in volt, power sensors at a node or a branch end
    ku = findfirst(r -> r["type"] == "U", arows)
    @test arows[ku]["sigma"] ≈ net.measurements[ku].sigma * Sparlectra.getNodeVn(net.nodeVec[net.measurements[ku].busIdx])
    @test first(adata["sym_voltage_sensor"])["u_measured"] > 1000.0
    @test first(adata["sym_power_sensor"])["measured_terminal_type"] in ("node", "branch_from", "branch_to")
    # every measurement kind, shunt-bay rows included: the export used to
    # drop ShuntQMeas and the bus-referenced ImagMeas without a word
    full = importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")))
    runpf!(full, 60, 1e-8, 0)
    full_meas = Sparlectra.generateMeasurementsFromPF(full; includeImag = true, includeIa = true, includeShuntQ = true, includeVa = true)
    empty!(full.measurements)
    append!(full.measurements, full_meas)
    bay = first(sh.busIdx for sh in full.shuntVec if sh.status == 1)
    Sparlectra.addMeasurement!(full.measurements; typ = Sparlectra.ImagMeas, value = 12.5, sigma = 0.5, busIdx = bay, id = "Imag_shunt_probe")
    Sparlectra.addMeasurement!(full.measurements; typ = Sparlectra.IaMeas, value = 91.0, sigma = 0.05, busIdx = bay, active = false, id = "Ia_shunt_probe")
    @test any(m -> m.typ == Sparlectra.ShuntQMeas, full.measurements)
    # redundant readings without an explicit id share the default id; one of
    # each pair goes first. The file position follows the row, not the id
    # (keyed by id, both rows got the last position and the reader moved the
    # first one to the back)
    dup = Sparlectra.Measurement[]
    for (typ, v) in ((Sparlectra.ShuntQMeas, 4.0), (Sparlectra.VmMeas, 1.01), (Sparlectra.ShuntQMeas, 4.2), (Sparlectra.VmMeas, 1.02))
      Sparlectra.addMeasurement!(dup; typ = typ, value = v, sigma = 0.5, busIdx = bay)
    end
    @test dup[1].id == dup[3].id && dup[2].id == dup[4].id
    pushfirst!(full.measurements, dup[1], dup[2])
    append!(full.measurements, dup[3:4])
    fa = exportSCF(full; file = joinpath(d, "full_a.scf.json"))
    full_back = importSCF(fa)
    @test length(full_back.measurements) == length(full.measurements)
    @test [(m.typ, m.busIdx, m.branchIdx, m.direction, m.id, m.value, m.sigma, m.active) for m in full_back.measurements] == [(m.typ, m.busIdx, m.branchIdx, m.direction, m.id, m.value, m.sigma, m.active) for m in full.measurements]
    @test read(exportSCF(full_back; file = joinpath(d, "full_b.scf.json")), String) == read(fa, String)
    mtuple(n) = [(m.typ, m.busIdx, m.branchIdx, m.direction, m.id, m.value, m.sigma, m.active) for m in n.measurements]
    # SCF since 1.2: one row per measurement with the table semantics; the shunt
    # bay is a row of its own with the value, current sensors carry ampere
    fdoc = Sparlectra.scf_json_parse(read(fa, String))
    @test fdoc["sparlectra"]["format_version"] == Sparlectra.SCF_FORMAT_VERSION == "1.3"
    rows = fdoc["sparlectra"]["measurements"]["rows"]
    @test !haskey(fdoc["sparlectra"]["measurements"], "shunt_bay_rows")
    @test length(rows) == length(full.measurements)
    @test String[r["type"] for r in rows] == String[Sparlectra._scf_table_type(m) for m in full.measurements]
    # an id is written only where it is not the default of its location
    # (the generated rows), and it reads back either way
    @test [get(r, "id", nothing) for r in rows] == [m.id == Sparlectra._default_measurement_id(m.typ; busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.branchIdx === nothing ? :none : m.direction) ? nothing : m.id for m in full.measurements]
    @test count(r -> haskey(r, "id"), rows) == 2
    bay_row = only(r for r in rows if get(r, "id", "") == "Imag_shunt_probe")
    @test haskey(bay_row["terminal"], "shunt_bay") && bay_row["type"] == "I" && bay_row["value"] == 12.5 && !haskey(bay_row, "sensor")
    @test only(r for r in rows if get(r, "id", "") == "Ia_shunt_probe")["valid"] == false
    @test all(r -> !haskey(r, "sensor") && haskey(r, "value"), [r for (r, m) in zip(rows, full.measurements) if m.typ == Sparlectra.ShuntQMeas])
    ki = findfirst(r -> r["type"] == "I" && haskey(r, "sensor"), rows)
    isensor = only(s for s in fdoc["data"]["sym_current_sensor"] if s["id"] == rows[ki]["sensor"])
    @test isensor["i_measured"] == full.measurements[ki].value
    # a U value no Float64 divides back to its per-unit value is written
    # with 18 digits and read at 256 bits, as the measurement table does
    vn = Sparlectra.getNodeVn(full.nodeVec[bay])
    pu18 = first(x for x in (0.9 + k * 1.0e-6 for k in 1:200_000) if Sparlectra._measurement_table_digits(Sparlectra._measurement_table_kv_text(x, vn * 1.0e3 / 1.0e3)) > 17)
    long = deepcopy(full)
    # (the sigma is always written in kV, the value only where the sensor
    # does not read back, so the sigma carries the long number here)
    Sparlectra.addMeasurement!(long.measurements; typ = Sparlectra.VmMeas, value = pu18, sigma = pu18, busIdx = bay, id = "Vm_long")
    fl = exportSCF(long; file = joinpath(d, "long.scf.json"))
    @test occursin("\"sigma\": " * Sparlectra._measurement_table_kv_text(pu18, vn * 1.0e3 / 1.0e3) * ",", read(fl, String))
    @test last(importSCF(fl).measurements).value === pu18 && last(importSCF(fl).measurements).sigma === pu18
    @test read(exportSCF(importSCF(fl); file = joinpath(d, "long_b.scf.json")), String) == read(fl, String)
    # a zero-injection row carries its flag in the row (written only when
    # true); in 1.2 the id decides nothing, so a row with the ZI id shape
    # but without the flag stays an ordinary measurement
    zinet = deepcopy(full)
    empty!(zinet.measurements)
    Sparlectra.addZeroInjectionMeasurements!(zinet.measurements; net = zinet, busIdxs = [bay])
    Sparlectra.addMeasurement!(zinet.measurements; typ = Sparlectra.PinjMeas, value = 0.0, sigma = 0.001, busIdx = bay, id = "ZI_PINJ_bus_99")
    Sparlectra.addMeasurement!(zinet.measurements; typ = Sparlectra.QinjMeas, value = 0.5, sigma = 0.001, busIdx = bay, id = "meter_q", zero_injection = true)
    fz = exportSCF(zinet; file = joinpath(d, "zi.scf.json"))
    @test [get(r, "zero_injection", false) for r in Sparlectra.scf_json_parse(read(fz, String))["sparlectra"]["measurements"]["rows"]] == [true, true, false, true]
    @test [m.zero_injection for m in importSCF(fz).measurements] == [true, true, false, true]
    @test read(exportSCF(importSCF(fz); file = joinpath(d, "zi_b.scf.json")), String) == read(fz, String)
    # only value and sigma of measurements.rows keep a long literal as
    # text; anywhere else it is a Float64 (parsed as before), and the writer
    # never refuses what the reader produced
    for lit in ("0.10000000000000000000", "0.123456789012345678", "110.000000000000000000", "1.0000000000000000001e-8", "12345678901234567890.5")
      doc = Sparlectra.scf_json_parse("{\"a\": $(lit), \"sparlectra\": {\"measurements\": {\"rows\": [{\"value\": $(lit), \"type\": \"P\"}]}}}")
      @test doc["a"] isa Float64 && doc["a"] === parse(Float64, lit)
      cell = doc["sparlectra"]["measurements"]["rows"][1]["value"]
      @test cell isa Sparlectra.ScfDecimalText && cell.text == lit && cell.value === parse(Float64, lit)
      out = Sparlectra.scf_json_string(doc)
      @test occursin("\"a\": $(repr(parse(Float64, lit)))", out) && occursin("\"value\": $(lit)", out)
    end
    # a 1.2 file whose U row carries such a literal (a hand edit) reads and
    # writes back unchanged through the typed reader and writer
    hand = replace(read(fl, String), r"\"type\": \"U\"\n( *)\}" => s"\"type\": \"U\",\n\1  \"value\": 110.000000000000000000\n\1}"; count = 1)
    @test occursin("110.000000000000000000", hand)
    fh = joinpath(d, "hand.scf.json")
    write(fh, hand)
    @test any(m -> m.typ == Sparlectra.VmMeas, importSCF(fh).measurements)
    fh2 = joinpath(d, "hand_b.scf.json")
    write_scf_json(read_scf_json(fh), fh2)
    @test read(fh2, String) == hand
    # reader validation: a wrong type at a node, an undefined node, a sensor
    # at another location and an unknown key are named
    # (row 2 is a U at a node with its sensor, row 1 a shunt-bay Q)
    @test rows[2]["type"] == "U" && haskey(rows[1]["terminal"], "shunt_bay")
    psens = first(r["sensor"] for r in rows if r["type"] == "P")
    for (k, f) in ((2, r -> r["zero_injection"] = true), (2, r -> r["type"] = "I"), (2, r -> r["node"] = 99999), (2, r -> r["sensor"] = psens), (2, r -> r["mystery"] = 1), (1, r -> r["terminal"] = Dict("shunt_bay" => 99999)), (1, r -> r["type"] = "P"))
      bad = deepcopy(fdoc)
      f(bad["sparlectra"]["measurements"]["rows"][k])
      @test_throws ArgumentError Sparlectra.scf_to_net(bad)
    end
    # a sensor no row names (added by hand) is read after the rows with its
    # default id; the rows keep their values
    dropped = deepcopy(fdoc)
    gone = popat!(dropped["sparlectra"]["measurements"]["rows"], 2)
    moved = Sparlectra.scf_to_net(dropped).measurements
    @test length(moved) == length(full.measurements)
    @test last(moved).typ == Sparlectra.VmMeas && isapprox(last(moved).value, full.measurements[2].value; rtol = 1e-14) && last(moved).id != full.measurements[2].id && !haskey(gone, "id")
    @test mtuple((measurements = moved[1:(end-1)],)) == mtuple((measurements = full.measurements[[1; 3:end]],))
    # power flow: identical iterations and identical voltages
    i1, e1 = runpf!(net, 60, 1e-8, 0)
    i2, e2 = runpf!(back, 60, 1e-8, 0)
    @test (i1, e1) == (i2, e2)
    @test [n._vm_pu for n in back.nodeVec] == [n._vm_pu for n in net.nodeVec]
    @test [n._va_deg for n in back.nodeVec] == [n._va_deg for n in net.nodeVec]
    # a row the file has no slot for (a link flow) is named, not dropped
    # silently
    @test !isempty(back.linkVec)
    Sparlectra.addMeasurement!(back.measurements; typ = Sparlectra.PflowMeas, value = 1.0, sigma = 0.1, linkIdx = 1)
    nf = @test_logs (:warn, r"link flow row\(s\) are not written") match_mode = :any exportSCF(back; file = joinpath(d, "noslot.scf.json"), kw...)
    @test length(importSCF(nf).measurements) == length(back.measurements) - 1
    # The shunt of a shunt voltage controller belongs to the controller: the
    # export leaves it out and the controller entry rebuilds it. Written as a
    # plain shunt too, sp_case60 gained one shunt per round trip, and after a
    # control run that row held the solved susceptance (twice the
    # compensation). The present susceptance travels as the start value.
    svc = importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case60.scf.json")))
    n_shunts = length(svc.shuntVec)
    sa = exportSCF(svc; file = joinpath(d, "svc_a.scf.json"), include_start_state = true)
    svc_back = importSCF(sa)
    svc_again = importSCF(exportSCF(svc_back; file = joinpath(d, "svc_b.scf.json"), include_start_state = true))
    @test (length(svc_back.shuntVec), length(svc_again.shuntVec)) == (n_shunts, n_shunts)
    @test read(joinpath(d, "svc_b.scf.json"), String) == read(sa, String)
    runpf!(svc, 60, 1e-8, 0)
    runpf!(svc_again, 60, 1e-8, 0)
    @test [n._vm_pu for n in svc_again.nodeVec] == [n._vm_pu for n in svc.nodeVec]
    @test Sparlectra.run_control!(svc; verbose = 0).converged
    svc_ctrl = only(c for c in collect_outer_controllers(svc) if c isa Sparlectra.ShuntVoltageControl)
    @test svc_ctrl.bs_mvar != 0.0
    solved = importSCF(exportSCF(svc; file = joinpath(d, "svc_solved.scf.json")))
    @test length(solved.shuntVec) == n_shunts
    @test only(c for c in collect_outer_controllers(solved) if c isa Sparlectra.ShuntVoltageControl).bs_mvar == svc_ctrl.bs_mvar
    Sparlectra.run_control!(solved; verbose = 0)
    @test maximum(abs(solved.nodeVec[i]._vm_pu - svc.nodeVec[i]._vm_pu) for i in eachindex(svc.nodeVec)) < 1e-8
    # the strict PGM file drops the controller, so there the shunt stays as
    # a plain shunt with its present susceptance
    strict_svc = @test_logs (:warn, r"strict_pgm = true") match_mode = :any exportSCF(svc; file = joinpath(d, "svc.pgm.json"), strict_pgm = true)
    @test length(Sparlectra.scf_json_parse(read(strict_svc, String))["data"]["shunt"]) == n_shunts
    # The state of a controller shunt travels in its controller entry: the
    # SVC shunt out of service, a switched bank (MSC) shunt released for
    # estimation. The bank sits one ulp off its block grid, as after a
    # solve, and a plain shunt added after both controllers gets another
    # index on import; all of it survives export, import, export bytewise.
    st_net = importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case60.scf.json")))
    st_svc = only(c for c in collect_outer_controllers(st_net) if c isa Sparlectra.ShuntVoltageControl)
    st_net.shuntVec[st_svc.shunt_idx].status = 0
    msc = addShuntVoltageControl!(st_net; bus = "Waldsee_110", target_vm_pu = 1.0, bs_min_mvar = -10.0, bs_max_mvar = 10.0, step_mvar = 5.0)
    st_net.shuntVec[msc.shunt_idx].estimate = true
    msc.bs_mvar = nextfloat(5.0)
    addShunt!(net = st_net, busName = "Lindeck_110", pShunt = 0.0, qShunt = 3.0)
    sta = exportSCF(st_net; file = joinpath(d, "st_a.scf.json"))
    st_back = importSCF(sta)
    stb = exportSCF(st_back; file = joinpath(d, "st_b.scf.json"))
    @test read(stb, String) == read(sta, String)
    ctrl_shunt(n, bus) = n.shuntVec[only(c for c in collect_outer_controllers(n) if c isa Sparlectra.ShuntVoltageControl && c.bus == bus).shunt_idx]
    @test (ctrl_shunt(st_back, st_svc.bus).status, ctrl_shunt(st_back, st_svc.bus).estimate) == (0, false)
    @test (ctrl_shunt(st_back, "Waldsee_110").status, ctrl_shunt(st_back, "Waldsee_110").estimate) == (1, true)
    @test length(st_back.shuntVec) == length(st_net.shuntVec)
    # the entry key is for the shunt of a shunt voltage controller only
    st_root = Sparlectra.scf_json_parse(read(sta, String))
    st_bad = deepcopy(st_root)
    st_bad["sparlectra"]["components"]["controllers"][st_svc.name]["shunt_state"]["position"] = 1
    @test_throws ArgumentError Sparlectra.scf_to_net(st_bad)
    # state estimation: identical objective and identical state
    se_src = _scf_test_net()
    Sparlectra.readMeasurementTable!(se_src; file = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.measurements.csv")), zero_injection = false)
    se_file = exportSCF(se_src; file = joinpath(d, "se.scf.json"), kw...)
    se_back = importSCF(se_file)
    expected_se = (r"topology precheck reported",)
    s1 = run_with_expected_warnings(() -> runse!(se_src), expected_se)
    s2 = run_with_expected_warnings(() -> runse!(se_back), expected_se)
    @test s1.objectiveJ == s2.objectiveJ
    @test s1.dof == s2.dof
    @test s1.iterations == s2.iterations
    @test s1.voltages == s2.voltages
  end)() end

  # The tap changer sets of the state estimation (Branch.se_tap_changers)
  # travel in the optional key `estimation` of a tap_changer_models entry
  # A hand-built net with a stored set of every kind
  # (ratio on the to side with the typed impedance correction, an additive
  # ratio with a linear phase changer, asymmetrical and symmetrical with
  # reactance curves, a phase table with r, x, g, b, a ratio table with its
  # outside increment, an empty set) and two typed windings off their
  # neutral steps, one with a stored set: write, read, write gives the same
  # bytes, the same live ratio and shift, and every set field by field
  # (before base_ratio and base_shift_deg each trip applied the steps once
  # more, and the writer refused an r or x set by the model). The same net without stored sets writes no key, reads with no
  # set and the same admittance matrix (the key does not touch the power
  # flow model). A set of origin typed or fields is refused on write and on
  # read. The CGMES deliveries of A4b (the own export of typed models, the
  # 2W RatioTapChangerTable with r, x, g, b, the 3W with two tapped legs)
  # keep their sets, and the first export is the fixed point of the writer
  # (with meta.cgmes_ids the re-read keeps the mRIDs). A tap estimation on
  # the re-read delivery gives the voltages and tap rows of the original.
  @testset "SCF tap_changer_models estimation round trip" begin (function ()
    S = Sparlectra
    # every field of two sets (tables and curves included)
    function same_set(a, b)
      (a === nothing || b === nothing) && return a === b
      (a.ratio_n, a.shift_n_deg, a.origin, length(a.changers)) == (b.ratio_n, b.shift_n_deg, b.origin, length(b.changers)) || return false
      return all(k -> all(f -> getfield(a.changers[k], f) == getfield(b.changers[k], f), fieldnames(S.SETapChanger)), eachindex(a.changers))
    end
    sets_of(n) = [br.se_tap_changers for br in n.branchVec]
    models_of(file) = S.scf_json_parse(read(file, String))["sparlectra"]["tap_changer_models"]
    d = mktempdir()

    # --- hand-built stored sets
    net = Net(name = "scf_sets", baseMVA = 100.0)
    for b in ("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
      addBus!(net = net, busName = b, vn_kV = 110.0)
    end
    addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
    addProsumer!(net = net, busName = "J", type = "ENERGYCONSUMER", p = 10.0, q = 2.0)
    # the typed windings sit off their neutral steps: branch 8 a ratio model
    # and a symmetrical changer with its reactance curve, branch 9 an
    # asymmetrical changer with the impedance correction (both set r or x
    # from the model)
    buses = ("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
    for k in 1:7
      addPIModelTrafo!(net = net, fromBus = buses[k], toBus = buses[k + 1], r_pu = 0.002, x_pu = 0.05, b_pu = 0.0, status = 1, ratio = 1.0, shift_deg = k == 7 ? -3.0 : 0.0)
    end
    addPIModelTrafo!(net = net, fromBus = "H", toBus = "I", r_pu = 0.002, x_pu = 0.05, b_pu = 0.0, status = 1, ratio = 1.0, shift_deg = 0.0, taps = PowerTransformerTaps(Vn_kV = 110.0, step = 2, lowStep = -5, highStep = 5, neutralStep = 0, voltageIncrement_kV = 1.1), phase_taps = PhaseTapChangerModel(kind = :symmetrical, step = 3, lowStep = -8, highStep = 8, neutralStep = 0, voltage_step_increment = 0.02, x_min = 0.05, x_max = 0.06))
    addPIModelTrafo!(net = net, fromBus = "I", toBus = "J", r_pu = 0.002, x_pu = 0.05, b_pu = 0.0, status = 1, ratio = 1.0, shift_deg = 0.0, phase_taps = PhaseTapChangerModel(kind = :asymmetrical, step = 3, lowStep = -10, highStep = 10, neutralStep = 0, voltage_step_increment = 0.01, winding_connection_angle_deg = 60.0, convention = :reciprocal_from_side), tap_changer_model = :impedance_correction)
    @test net.branchVec[8].x_pu != net.branchVec[8].x_base_pu && net.branchVec[9].r_pu != net.branchVec[9].r_base_pu
    pts(f) = [S.SETapPoint(step = k, ratio = f(k), angle_deg = 0.5 * k, r = 1.0 + 0.01 * k, x = 1.0 + 0.02 * k, g = 1.0 - 0.03 * k, b = 1.0 + 0.04 * k) for k in -2:2]
    C(; kw...) = S.SETapChanger(; kw...)
    stored = Dict(
      1 => S.SETapChangerSet(ratio_n = 1.02, shift_n_deg = 0.0, origin = :cgmes, changers = [C(name = "rtc", mrid = "_m1", kind = :ratio, role = :ratio, side = :to, du = 0.0125, neutral = 0, low = -10, high = 10, step = 2.5, rx_correction = true)]),
      2 => S.SETapChangerSet(ratio_n = 0.98, shift_n_deg = 1.5, origin = :powsybl, changers = [C(name = "lin", kind = :ratio, role = :ratio, side = :from, additive = true, step_abs = 0.01, neutral = 0, low = -5, high = 5, step = -1.0), C(name = "ptc", kind = :linear, role = :phase, side = :to, dphi_deg = 1.2, neutral = 2, low = 0, high = 4, step = 3.0, x_curve = (0.04, 0.05, 2.4), x_curve_form = :tan2, impedance = :x_curve)]),
      3 => S.SETapChangerSet(ratio_n = 1.0, origin = :cgmes, changers = [C(name = "asym", kind = :asymmetrical, role = :phase, side = :from, du = 0.015, alpha_deg = 60.0, neutral = 13, low = 1, high = 25, step = 16.0, x_curve = (0.05, 0.06, 18.2), impedance = :x_curve, rx_correction = true)]),
      4 => S.SETapChangerSet(ratio_n = 1.0, origin = :cgmes, changers = [C(name = "sym", kind = :symmetrical, role = :phase, side = :to, du = 0.02, neutral = 0, low = -8, high = 8, step = -3.0)]),
      5 => S.SETapChangerSet(ratio_n = 1.0, origin = :cgmes, changers = [C(name = "ptab", kind = :tabular, role = :phase, side = :from, neutral = 0, low = -2, high = 2, step = 1.0, table = pts(k -> 1.0 + 0.001 * k^2), impedance = :table)]),
      6 => S.SETapChangerSet(ratio_n = 1.01, origin = :cgmes, changers = [C(name = "rtab", kind = :tabular, role = :ratio, side = :from, neutral = 0, low = -10, high = 10, step = -4.0, table = pts(k -> 1.0 + 0.015 * k), outside_du = 0.015, impedance = :table)]),
      7 => S.SETapChangerSet(ratio_n = 1.0, shift_n_deg = -3.0, origin = :cgmes),
      9 => S.SETapChangerSet(ratio_n = 1.0, origin = :cgmes, changers = [C(name = "asym9", kind = :asymmetrical, role = :phase, side = :to, du = 0.01, alpha_deg = 90.0, neutral = 0, low = -10, high = 10, step = 3.0)]),
    )
    for (k, st) in stored
      net.branchVec[k].se_tap_changers = st
    end
    a = exportSCF(net; file = joinpath(d, "sets_a.scf.json"))
    back = importSCF(a)
    b = exportSCF(back; file = joinpath(d, "sets_b.scf.json"))
    @test read(a, String) == read(b, String)
    # phase_tap_declared: branches 1 to 7 carry the undeclared default phase
    # grid, which the writer marks (before, the re-read declared it)
    live(n) = [(br.tap_ratio, br.phase_shift_deg, br.r_pu, br.x_pu, br.phase_tap_declared) for br in n.branchVec]
    @test live(back) == live(net)
    @test all(k -> same_set(sets_of(net)[k], sets_of(back)[k]), eachindex(net.branchVec))
    # restoreBaseImpedances! puts the base r and x on the typed windings; the
    # file describes the same equipment (before, the guard refused it)
    rest = restoreBaseImpedances!(deepcopy(net))
    @test read(exportSCF(rest; file = joinpath(d, "sets_rest.scf.json")), String) == read(a, String)
    @test count(st -> st !== nothing, sets_of(back)) == 8
    # entry shape: only the key on a branch without typed models, both on
    # branch 9, none on the typed branch 8 without a stored set
    ids = Dict(get(e, "branch_index", 0) => k for (k, e) in S.scf_json_parse(read(a, String))["sparlectra"]["extra"] if e isa AbstractDict)
    m = models_of(a)
    @test all(k -> collect(keys(m[ids[k]])) == ["estimation"], 1:7)
    @test isempty(m[ids[7]]["estimation"]["changers"])
    @test !haskey(m[ids[8]], "estimation") && haskey(m[ids[8]], "phase")
    @test haskey(m[ids[9]], "estimation") && haskey(m[ids[9]], "phase")
    # a typed winding still resolves; without a stored set its set is derived
    @test back.branchVec[8].taps_derived && back.branchVec[9].taps_derived
    @test S._effective_tap_changers(back.branchVec[8]).origin === :typed
    @test S._effective_tap_changers(back.branchVec[9]).origin === :cgmes
    # an empty set stays the declaration "no tap changer"
    @test S._tap_release_mode(back.branchVec[7]) === nothing
    # without stored sets: no key, no set, the same model
    plain = deepcopy(net)
    foreach(br -> br.se_tap_changers = nothing, plain.branchVec)
    p1 = exportSCF(plain; file = joinpath(d, "plain_a.scf.json"))
    @test all(e -> !haskey(e, "estimation"), values(models_of(p1))) && length(models_of(p1)) == 2
    pback = importSCF(p1)
    @test all(st -> st === nothing, sets_of(pback))
    @test S.createYBUS(net = pback, sparse = true) == S.createYBUS(net = back, sparse = true)
    @test [br.tap_ratio for br in pback.branchVec] == [br.tap_ratio for br in back.branchVec]
    @test read(exportSCF(pback; file = joinpath(d, "plain_b.scf.json")), String) == read(p1, String)
    # the reader names a broken estimation object
    root = S.scf_json_parse(read(a, String))
    function refused(f)
      r = deepcopy(root)
      f(r["sparlectra"]["tap_changer_models"])
      try
        S.scf_to_net(r)
      catch err
        err isa ArgumentError && return err.msg
        rethrow()
      end
      return ""
    end
    # no step: a stored changer without `step` sits at its
    # neutral step, typed models without `step` at their neutral_step, a
    # nameplate controller without `pos` at its neutral position 0; each is
    # marked and named in one info line, the rest of the file is unchanged
    r0 = deepcopy(root)
    tm = r0["sparlectra"]["tap_changer_models"]
    delete!(tm[ids[3]]["estimation"]["changers"][1], "step")
    delete!(tm[ids[8]]["ratio"], "step")
    delete!(tm[ids[8]]["phase"], "step")
    row4 = only(t for t in r0["sparlectra"]["components"]["tap_changer"] if t["branch"] == parse(Int, ids[4]))
    foreach(c -> delete!(c, "pos"), row4["controllers"])
    nb = @test_logs (:info, "SCF: tap changer branch $(ids[3])/asym: step not transferred, middle position 13 assumed") (:info, "SCF: tap changer branch $(ids[8]): step not transferred, middle position 0 assumed") (:info, "SCF: tap changer branch $(ids[4]): step not transferred, middle position 0 assumed") match_mode = :any S.scf_to_net(r0)
    c3 = nb.branchVec[3].se_tap_changers.changers[1]
    @test (c3.step, c3.step_assumed) == (13.0, true)
    @test (nb.branchVec[8].tap_winding.taps.step, nb.branchVec[8].tap_winding.phase_taps.step) == (0, 0)
    @test nb.branchVec[8].ratio_step_assumed && nb.branchVec[8].phase_step_assumed
    @test nb.branchVec[4].ratio_step_assumed && nb.branchVec[4].tap_ratio == nb.branchVec[4].ratio
    @test !any(br -> br.ratio_step_assumed || br.phase_step_assumed, back.branchVec) && !back.branchVec[3].se_tap_changers.changers[1].step_assumed
    # the writer puts the middle position into the file with the mark
    # step_assumed next to it; the re-read marks the same changers at the
    # same steps, names them again, and a second write is the first byte
    # for byte
    ma = exportSCF(nb; file = joinpath(d, "assumed_a.scf.json"))
    mroot = S.scf_json_parse(read(ma, String))
    mm = mroot["sparlectra"]["tap_changer_models"]
    c3w = mm[ids[3]]["estimation"]["changers"][1]
    @test c3w["step_assumed"] === true && c3w["step"] == 13
    @test mm[ids[8]]["ratio"]["step_assumed"] === true && mm[ids[8]]["phase"]["step_assumed"] === true && mm[ids[8]]["ratio"]["step"] == 0
    mrow4 = only(t for t in mroot["sparlectra"]["components"]["tap_changer"] if t["branch"] == parse(Int, ids[4]))
    @test all(c -> get(c, "step_assumed", false) === true && c["pos"] == 0, mrow4["controllers"])
    @test !occursin("step_assumed", read(a, String))
    mb = @test_logs (:info, "SCF: tap changer branch $(ids[3])/asym: step not transferred, middle position 13 assumed") (:info, "SCF: tap changer branch $(ids[8]): step not transferred, middle position 0 assumed") (:info, "SCF: tap changer branch $(ids[4]): step not transferred, middle position 0 assumed") match_mode = :any importSCF(ma)
    marks(n) = [(br.ratio_step_assumed, br.phase_step_assumed, br.se_tap_changers === nothing ? Bool[] : [c.step_assumed for c in br.se_tap_changers.changers]) for br in n.branchVec]
    @test marks(mb) == marks(nb)
    @test live(mb) == live(nb)
    @test mb.branchVec[3].se_tap_changers.changers[1].step == 13.0
    @test read(exportSCF(mb; file = joinpath(d, "assumed_b.scf.json")), String) == read(ma, String)
    # the write path of a QSTS carry moves the typed ratio model off its
    # middle position: its mark goes, the phase model keeps it
    nc = deepcopy(nb)
    @test S._set_model_tap_steps!(nc.branchVec[8], 1, typemin(Int))
    @test !nc.branchVec[8].ratio_step_assumed && nc.branchVec[8].phase_step_assumed
    mc = S.scf_json_parse(read(exportSCF(nc; file = joinpath(d, "assumed_c.scf.json")), String))["sparlectra"]["tap_changer_models"][ids[8]]
    @test (mc["ratio"]["step"], haskey(mc["ratio"], "step_assumed"), get(mc["phase"], "step_assumed", false)) == (1, false, true)
    e1 = ids[1]
    @test occursin("tap_changer_models $(e1) estimation.changers[1] has the unknown key \"dU\"", refused(t -> t[e1]["estimation"]["changers"][1]["dU"] = 0.01))
    @test occursin("tap_changer_models $(e1) estimation.changers[1]: SETapChanger rtc: kind must be one of", refused(t -> t[e1]["estimation"]["changers"][1]["kind"] = "rotary"))
    @test occursin("tap_changer_models $(e1) estimation is missing the required key \"origin\"", refused(t -> delete!(t[e1]["estimation"], "origin")))
    for o in ("typed", "fields", "scf")
      @test occursin("tap_changer_models $(e1) estimation.origin must be \"cgmes\" or \"powsybl\", got \"$(o)\"", refused(t -> t[e1]["estimation"]["origin"] = o))
    end
    # the tap impedance correction is defined for the reciprocal convention:
    # branch 9 (asymmetrical, impedance_correction) in the direct convention
    # reads with exactly one warning naming its entry, the file as written
    # (reciprocal) reads without a warning
    e9 = ids[9]
    direct = deepcopy(root)
    direct["sparlectra"]["tap_changer_models"][e9]["phase"]["convention"] = "direct_regulating_vector"
    @test_logs (:warn, Regex("SCF: tap_changer_models $(e9): the tap impedance correction is defined for the convention reciprocal_from_side")) min_level = Logging.Warn S.scf_to_net(direct)
    @test_logs min_level = Logging.Warn S.scf_to_net(deepcopy(root))
    # a derived set put on a branch is no import record: the writer refuses it
    hand = deepcopy(net)
    hand.branchVec[8].se_tap_changers = S._effective_tap_changers(back.branchVec[8])
    @test hand.branchVec[8].se_tap_changers.origin === :typed
    @test_throws ArgumentError exportSCF(hand; file = joinpath(d, "hand.scf.json"))

    # --- the CGMES deliveries of A4b
    cn = Net(name = "sets", baseMVA = 100.0)
    for bn in ("A", "B", "C", "D", "E", "F", "G", "H")
      addBus!(net = cn, busName = bn, vn_kV = 110.0)
    end
    addProsumer!(net = cn, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
    addProsumer!(net = cn, busName = "H", type = "ENERGYCONSUMER", p = 10.0, q = 2.0)
    casym() = PhaseTapChangerModel(kind = :asymmetrical, step = 16, lowStep = 1, highStep = 25, neutralStep = 13, voltage_step_increment = 0.0125, winding_connection_angle_deg = 90.0, convention = :direct_regulating_vector)
    sym = PhaseTapChangerModel(kind = :symmetrical, step = 3, lowStep = -8, highStep = 8, neutralStep = 0, voltage_step_increment = 0.02, x_min = 0.05, x_max = 0.06, convention = :direct_regulating_vector)
    tab = PhaseTapChangerModel(kind = :tabular, step = 2, neutralStep = 0, convention = :direct_regulating_vector, table = [TapTablePoint(step = k, ratio = 1.0 + 0.001 * k^2, angle_deg = 2.0 * k) for k in -3:3])
    rtc = PowerTransformerTaps(Vn_kV = 110.0, step = 2, lowStep = -5, highStep = 5, neutralStep = 0, voltageIncrement_kV = 1.1)
    tr(f, t; shift = 0.0, phase_taps = nothing, taps = nothing) = addPIModelTrafo!(net = cn, fromBus = f, toBus = t, r_pu = 0.002, x_pu = 0.05, b_pu = 0.0, status = 1, ratio = 1.0, shift_deg = shift, phase_taps = phase_taps, taps = taps)
    tr("A", "B"; phase_taps = casym())
    tr("B", "C"; phase_taps = sym)
    tr("C", "D"; phase_taps = tab)
    tr("D", "E"; taps = rtc)
    tr("E", "F")
    tr("F", "G"; shift = -3.0)
    tr("G", "H"; phase_taps = casym())
    cn.trafos[end].side2.taps = PowerTransformerTaps(Vn_kV = 110.0, step = -1, lowStep = -4, highStep = 4, neutralStep = 0, voltageIncrement_kV = 1.65)
    cdir = mktempdir()
    writeCGMESFiles(cn; path = cdir)
    deliveries = [
      ("typed export", importCGMES(path = cdir, name = "sets", require_boundary = false).net),
      ("2W ratio table", importCGMES(path = _cgmes_synth_tapside_dir(:two_winding; table = (r = 4.0, x = 10.0, g = -20.0, b = 15.0), control_kv = 21.2), name = "tap2w_table").net),
      ("3W two tapped legs", importCGMES(path = _cgmes_synth3w_dir(end1_tap_step = 13), name = "synth3w_sets").net),
    ]
    for (label, src) in deliveries
      f1 = exportSCF(src; file = joinpath(d, "cg_1.scf.json"), source_format = "cgmes")
      n1 = importSCF(f1)
      f2 = exportSCF(n1; file = joinpath(d, "cg_2.scf.json"), source_format = "cgmes")
      f3 = exportSCF(importSCF(f2); file = joinpath(d, "cg_3.scf.json"), source_format = "cgmes")
      ok = read(f1, String) == read(f2, String) == read(f3, String)
      ok &= all(k -> same_set(sets_of(src)[k], sets_of(n1)[k]), eachindex(src.branchVec))
      ok &= count(st -> st !== nothing, sets_of(n1)) == count(br -> br.ratio != 0.0, src.branchVec)
      @test ok
      ok || println("      SCF estimation round trip of the ", label, " delivery: sets or bytes differ")
    end
    # tap estimation on the delivery and on its re-read: a phase and a
    # ratio changer released, exact measurements of the import state
    src = deliveries[1][2]
    reread = importSCF(exportSCF(src; file = joinpath(d, "cg_se.scf.json"), source_format = "cgmes"))
    truth = deepcopy(src)
    ite, erg = runpf!(truth, 40, 1e-12, 0; method = :rectangular)
    @test erg == 0
    meas = generateMeasurementsFromPF(truth; includeImag = true, includeVa = true)
    runs = map((src, reread)) do n
      setTapEstimation!(n; trafo = 1, mode = :pst)
      setTapEstimation!(n; trafo = 4, mode = :ratio)
      with_state_estimation_config(max_iter = 40, tol = 1e-8, update_net = false) do
        runse!(n, meas)
      end
    end
    @test runs[1].converged && runs[2].converged
    @test maximum(abs.(runs[1].voltages .- runs[2].voltages)) <= 1e-12
    @test [(r.branch, r.kind, r.step_fixed) for r in runs[1].tapEstimates] == [(r.branch, r.kind, r.step_fixed) for r in runs[2].tapEstimates] == [(1, :asymmetrical, 16.0), (4, :ratio, 2.0)]
    @test all(k -> abs(runs[1].tapEstimates[k].s_est - runs[2].tapEstimates[k].s_est) <= 1e-9, 1:2)
  end)() end

  @testset "line shunt conductance without capacitance (g1)" begin (function ()
    # tan1 = g/b cannot spell a conductance when b == 0 (a CGMES gch with
    # no bch, say): the writer falls back to the direct g1 field, the
    # reader prefers it, and the strict-PGM writer strips it BY NAME
    gnet = Net(name = "pure_g", baseMVA = 100.0)
    addBus!(net = gnet, busName = "A", vn_kV = 110.0)
    addBus!(net = gnet, busName = "B", vn_kV = 110.0)
    addPIModelACLine!(net = gnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, g_pu = 3.0e-4, status = 1)
    addProsumer!(net = gnet, busName = "A", type = "ExternalNetworkInjection", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
    addProsumer!(net = gnet, busName = "B", type = "ENERGYCONSUMER", p = 20.0, q = 5.0)
    dg = mktempdir()
    fg = exportSCF(gnet; file = joinpath(dg, "pure_g.scf.json"))
    row = only(Sparlectra.scf_json_parse(read(fg, String))["data"]["line"])
    @test row["tan1"] == 0.0
    @test haskey(row, "g1")
    gback = importSCF(fg)
    gb = only(gback.branchVec)
    @test (gb.g_pu, gb.b_pu) == (3.0e-4, 0.0)
    # the write-read-write fixed point holds with the extension field
    fg2 = exportSCF(gback; file = joinpath(dg, "pure_g2.scf.json"))
    @test read(fg, String) == read(fg2, String)
    fstrict = @test_logs (:warn, r"line shunt conductance") match_mode = :any exportSCF(gnet; file = joinpath(dg, "pure_g.pgm.json"), strict_pgm = true)
    srow = only(Sparlectra.scf_json_parse(read(fstrict, String))["data"]["line"])
    @test !haskey(srow, "g1")
  end)() end

  @testset "reader validation: unknown keys and broken references" begin (function ()
    fixture = abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_casePST.scf.json"))
    base = Sparlectra.scf_json_parse(read(fixture, String))
    # unknown keys are hard errors, at the root and inside the namespace
    bad_root = deepcopy(base)
    bad_root["unexpected"] = 1
    @test_throws ArgumentError Sparlectra.scf_to_net(bad_root)
    bad_ns = deepcopy(base)
    bad_ns["sparlectra"]["mystery"] = Dict{String,Any}()
    @test_throws ArgumentError Sparlectra.scf_to_net(bad_ns)
    bad_comp = deepcopy(base)
    bad_comp["data"]["asym_line"] = Any[]
    @test_throws ArgumentError Sparlectra.scf_to_net(bad_comp)
    # a dangling reference is named, not silently dropped
    bad_ref = deepcopy(base)
    bad_ref["data"]["line"][1]["to_node"] = 99999
    @test_throws ArgumentError Sparlectra.scf_to_net(bad_ref)
    # duplicate ids across component types are rejected
    dup = deepcopy(base)
    dup["data"]["line"][1]["id"] = dup["data"]["node"][1]["id"]
    @test_throws ArgumentError Sparlectra.scf_to_net(dup)
    # batch datasets are not an input form
    batch = deepcopy(base)
    batch["is_batch"] = true
    @test_throws ArgumentError Sparlectra.scf_to_net(batch)
    # the JSON parser reports the position of a malformed file
    @test_throws ArgumentError Sparlectra.scf_json_parse("{\"a\": }")
    @test_throws ArgumentError Sparlectra.scf_json_parse("[1, 2]")
    # the case configuration is readable without building the network
    cfgkeys = scf_case_config(fixture)
    @test cfgkeys isa Dict{String,Any}
    # the revision decides whether this reader may read the file at all, and
    # a setpoint without the regulation flag is a contradiction the reader
    # refuses instead of promoting the machine (moved from test_legacy_scf.jl)
    dv = mktempdir()
    shipped_export = exportSCF(importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json"))); file = joinpath(dv, "a_sp_case5.json"))
    stale = Sparlectra.scf_json_parse(read(shipped_export, String))
    stale["sparlectra"]["format_version"] = "0.9"
    staleF = joinpath(dv, "stale.scf.json")
    write(staleF, Sparlectra.scf_json_string(stale))
    @test_throws ["format revision 0.9 found", Sparlectra.SCF_FORMAT_VERSION] importSCF(staleF)
    @test_throws r"re-export"i importSCF(staleF)
    contra = Sparlectra.scf_json_parse(read(shipped_export, String))
    target = first(k for (k, v) in contra["sparlectra"]["extra"] if v isa AbstractDict && haskey(v, "prosumption_type") && v["prosumption_type"] == "Load")
    contra["sparlectra"]["extra"][target]["vm_pu"] = 1.0
    contraF = joinpath(dv, "contra.scf.json")
    write(contraF, Sparlectra.scf_json_string(contra))
    @test_throws ["voltage setpoint", "regulated"] importSCF(contraF)
    # a non-finite number never reaches the file as null: it is either
    # encoded or refused by name
    @test Sparlectra.scf_infinity_sentinel(Inf) == "inf"
    @test Sparlectra.scf_infinity_sentinel(-Inf) == "-inf"
    @test Sparlectra.scf_infinity_sentinel(2.5) == 2.5
    @test Sparlectra.scf_number_or_sentinel("inf", "x") == Inf
    @test_throws ArgumentError Sparlectra.scf_infinity_sentinel(NaN)
    @test_throws ArgumentError Sparlectra.scf_json_string(Dict{String,Any}("data" => Dict{String,Any}("bad" => Inf)))
  end)() end

  # The MATPOWER quirks the round trip once failed on, on a small net of its
  # own. They came from case14 and case57 in the large-case directory, which
  # a fresh install does not have: numeric bus names ("10" sorts before "2"),
  # an unlimited rating (rateA = 0), a neutral ratio 0.895 outside the 0.9 to
  # 1.1 band (the case57 transformer), and a bus shunt, read once with each
  # shunt model so the per-shunt model travels in shunt_state.
  @testset "round trip of MATPOWER quirks on a small net" begin (function ()
    d = mktempdir()
    src = joinpath(d, "scf_quirks.m")
    write(src, """
    function mpc = scf_quirks
    mpc.version = '2';
    mpc.baseMVA = 100;
    mpc.bus = [
    2 3 0 0 0 0 1 1.02 0 110 1 1.1 0.9;
    10 1 40 15 0 19 1 1.0 0 110 1 1.1 0.9;
    7 2 20 5 0 0 1 1.01 0 20 1 1.1 0.9;
    ];
    mpc.gen = [
    2 40 0 300 -300 1.02 100 1 300 0;
    7 20 0 20 -20 1.01 100 1 50 0;
    ];
    mpc.branch = [
    2 10 0.01 0.05 0.02 0 0 0 0 0 1 -360 360;
    10 7 0.005 0.08 0.0 100 100 100 0.895 0 1 -360 360;
    ];
    """)
    base_cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
    vdi_cfg = Sparlectra.SparlectraConfig(Dict("model" => Dict("bus_shunt_model" => "voltage_dependent_injection")))
    for (label, cfg) in (("admittance", base_cfg), ("voltage_dependent_injection", vdi_cfg))
      net = SparlectraApp._se_import_case_net(src, cfg)
      a = exportSCF(net; file = joinpath(d, "a_$(label).json"), source_format = "matpower")
      # the reader restores the band as stated in the file and says so
      back = @test_logs (:warn, r"neutral tap position lies outside") match_mode = :any importSCF(a)
      b = exportSCF(back; file = joinpath(d, "b_$(label).json"), source_format = "matpower")
      @test read(a, String) == read(b, String)
      # the generated component names embed the source bus number, so losing
      # it renamed everything and reshuffled the ids
      @test [getCompName(n.comp) for n in net.nodeVec] == [getCompName(n.comp) for n in back.nodeVec]
      @test [getCompName(br.comp) for br in net.branchVec] == [getCompName(br.comp) for br in back.branchVec]
      # bus order is identity: ids are name-sorted
      @test [Sparlectra._scf_bus_name(net, i) for i in eachindex(net.nodeVec)] == [Sparlectra._scf_bus_name(back, i) for i in eachindex(back.nodeVec)]
      @test isapprox(Matrix(createYBUS(net = net)), Matrix(createYBUS(net = back)); rtol = 1e-12)
      # every field by reflection, on the network as the service builds it
      # (the switching parameters are stamped on that path), with the start
      # state on both sides
      full = exportSCF(net; file = joinpath(d, "full_$(label).json"), source_format = "matpower", include_start_state = true)
      via_service = run_with_expected_warnings(() -> Sparlectra._import_sparlectra_net(full, nothing, cfg), (r"neutral tap position",))
      @test (label, scf_roundtrip_field_diffs(net, via_service)) == (label, String[])
      # a setpoint and a reactive band are written only for regulated units
      extra_written = Sparlectra.scf_json_parse(read(a, String))["sparlectra"]["extra"]
      @test !any(v isa AbstractDict && haskey(v, "vm_pu") && get(v, "regulated", false) !== true for (_, v) in extra_written)
      @test !any(v isa AbstractDict && haskey(v, "max_q_mvar") && !isa(v["max_q_mvar"], Real) for (_, v) in extra_written)
      # an unlimited rating has no PGM representation: it travels in the
      # namespaced block, no null is written, and it comes back as Inf
      unlimited = [i for i in eachindex(net.branchVec) if net.branchVec[i].sn_MVA !== nothing && !isfinite(net.branchVec[i].sn_MVA)]
      @test !isempty(unlimited)
      @test !occursin("\"i_n\": null", read(a, String))
      @test all(back.branchVec[i].sn_MVA == Inf for i in unlimited)
    end
    # the copy with the start state solves to the same point
    net = SparlectraApp._se_import_case_net(src, base_cfg)
    withstart = exportSCF(net; file = joinpath(d, "s_quirks.json"), source_format = "matpower", include_start_state = true)
    solved = run_with_expected_warnings(() -> importSCF(withstart), (r"neutral tap position",))
    it_a, _ = runpf!(net, 40, 1e-10, 0)
    it_b, _ = runpf!(solved, 40, 1e-10, 0)
    @test it_a == it_b
    @test maximum(abs(getNodeVm(net.nodeVec[i]) - getNodeVm(solved.nodeVec[i])) for i in eachindex(net.nodeVec)) < 1e-12
    @test maximum(abs(net.nodeVec[i]._va_deg - solved.nodeVec[i]._va_deg) for i in eachindex(net.nodeVec)) < 1e-10
  end)() end

  # DynamicData on the machines (#506): defaults and rules of the record,
  # and the block extra.<machine>.dynamics of revision 1.3. A machine with
  # every field, one with a subset and one without go through write and read
  # bit for bit; a case without the block writes none, and the same case as
  # a revision 1.2 file reads to the same network
  @testset "dynamic unit data and the SCF dynamics block" begin (function ()
    @test DynamicData().unit_type == GridFollowingUnit
    @test DynamicData(inertia_s = 4).unit_type == SynchronousUnit
    local dflt = DynamicData()
    @test (dflt.inertia_s, dflt.damping_pu, dflt.rating_mva, dflt.droop_pu, dflt.deadband_hz, dflt.primary_reserve_mw, dflt.governor_time_s, dflt.agc_participation, dflt.secondary_reserve_mw, dflt.virtual_inertia_s, dflt.black_start, dflt.min_stable_mw) == (nothing, 0.0, nothing, nothing, 0.0, 0.0, nothing, 0.0, 0.0, nothing, false, nothing)
    @test DynamicData(damping_pu = -0.0).damping_pu === 0.0
    # each rule fails with a DynamicDataError naming the field
    for (field, value) in ((:inertia_s, -1.0), (:virtual_inertia_s, -0.5), (:governor_time_s, -0.1), (:primary_reserve_mw, -1.0), (:secondary_reserve_mw, -2.0), (:deadband_hz, -0.01),
                           (:droop_pu, 0.0), (:droop_pu, -0.05), (:rating_mva, 0.0), (:agc_participation, 1.5), (:agc_participation, -0.1), (:damping_pu, NaN), (:min_stable_mw, Inf))
      local kw = NamedTuple{(field,)}((value,))
      @test_throws DynamicDataError DynamicData(; kw...)
      @test_throws Regex("^DynamicDataError: $(field) = ") DynamicData(; kw...)
    end

    local d = mktempdir()
    # the shipped case written once without start state, so the network
    # read back from it is the reference of every comparison below
    local plain = exportSCF(importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"))); file = joinpath(d, "plain.scf.json"))
    # a case without dynamic data writes no block
    @test !occursin("\"dynamics\"", read(plain, String))
    local net = importSCF(plain)
    local gens = [ps for ps in net.prosumpsVec if Sparlectra.isGenerator(ps)]
    @test length(gens) >= 3
    # 0.1 + 0.2 and 1/23 have no short decimal form: they test bit identity
    local full = DynamicData(inertia_s = 0.1 + 0.2, damping_pu = 0.7, rating_mva = 247.5, droop_pu = 0.045, deadband_hz = 0.01, primary_reserve_mw = 12.5, governor_time_s = 0.35,
                             agc_participation = 0.2, secondary_reserve_mw = 30.0, unit_type = GridFormingUnit, virtual_inertia_s = 2.1, black_start = true, min_stable_mw = 40.0)
    local part = DynamicData(inertia_s = 4.25, droop_pu = 1 / 23)
    gens[1].dynamics = full
    gens[2].dynamics = part
    # machines are told apart by their source id (several share a name)
    local id(ps) = ps.comp.cID
    local f = exportSCF(net; file = joinpath(d, "dyn.scf.json"))
    local doc = Sparlectra.scf_json_parse(read(f, String))
    @test doc["sparlectra"]["format_version"] == "1.3"
    local blocks = Dict(e["external_id"] => e["dynamics"] for e in values(doc["sparlectra"]["extra"]) if e isa AbstractDict && haskey(e, "dynamics"))
    @test sort(collect(keys(blocks))) == sort([id(gens[1]), id(gens[2])])
    @test sort(collect(keys(blocks[id(gens[1])]))) == sort(String.(collect(fieldnames(DynamicData))))
    @test sort(collect(keys(blocks[id(gens[2])]))) == ["agc_participation", "black_start", "damping_pu", "deadband_hz", "droop_pu", "inertia_s", "primary_reserve_mw", "secondary_reserve_mw", "unit_type"]
    @test (blocks[id(gens[1])]["unit_type"], blocks[id(gens[2])]["unit_type"]) == ("grid_forming", "synchronous")
    local back = importSCF(f)
    local by_id = Dict(id(ps) => ps for ps in back.prosumpsVec)
    @test length(by_id) == length(back.prosumpsVec)
    # === compares every field bit for bit
    @test by_id[id(gens[1])].dynamics === full
    @test by_id[id(gens[2])].dynamics === part
    @test by_id[id(gens[3])].dynamics === nothing
    @test count(ps -> ps.dynamics !== nothing, back.prosumpsVec) == 2
    @test isempty(scf_roundtrip_field_diffs(net, back))
    @test read(exportSCF(back; file = joinpath(d, "dyn_again.scf.json")), String) == read(f, String)
    # the round-trip comparison sees the field: a lost block is a difference
    local lost = importSCF(plain)
    @test any(m -> startswith(m, "prosumer.dynamics"), scf_roundtrip_field_diffs(net, lost))
    # a broken block is named: unknown key, a value against the rules, an
    # unknown unit type, a flag that is no Bool, a number given as text
    local gid = only(k for (k, e) in doc["sparlectra"]["extra"] if e isa AbstractDict && get(e, "external_id", "") == id(gens[2]))
    for (key, value, msg) in (("inertia", 4.0, "unknown key \"inertia\""), ("inertia_s", -1.0, "inertia_s = -1.0 must be >= 0"),
                              ("unit_type", "synchron", "unit_type must be"), ("black_start", 1, "black_start"), ("droop_pu", "0.05", "droop_pu"))
      local bad = deepcopy(doc)
      bad["sparlectra"]["extra"][gid]["dynamics"][key] = value
      @test_throws ArgumentError Sparlectra.scf_to_net(bad)
      @test_throws msg Sparlectra.scf_to_net(bad)
    end
  end)() end

  @testset "system frequency: meta.f_nom, the key, and the MATPOWER path" begin (function ()
    d = mktempdir()
    key_cfg(f; user_set = false) = SparlectraConfig(model = ModelConfig(nominal_frequency_hz = f), user_set_keys = user_set ? Set(["model.nominal_frequency_hz"]) : Set{String}())
    sp14 = abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"))
    vstate(n) = ([nd._vm_pu for nd in n.nodeVec], [nd._va_deg for nd in n.nodeVec])

    # round trip at 60 Hz: export, import, the net carries 60, and the
    # power flow is bit identical to the same export at 50 Hz (b in pu does
    # not depend on f; c1 = b/(2*pi*f) and back). sp_case14 has no line
    # charging; the Zeng case14 has six lines with c1 > 0.
    zeng = abspath(joinpath(dirname(@__DIR__), "data", "scf", "case14_zeng_p306_activeSet_A.scf.json"))
    f60 = ""
    net60 = nothing
    for (tag, file) in (("sp14", sp14), ("zeng", zeng))
      orig = importSCF(file)
      @test orig.nominal_frequency_hz == 50.0
      f50 = exportSCF(orig; file = joinpath(d, "$(tag)_50.scf.json"))
      src60 = importSCF(file)
      src60.nominal_frequency_hz = 60.0
      f60 = exportSCF(src60; file = joinpath(d, "$(tag)_60.scf.json"))
      net60 = src60
      r50 = Sparlectra.scf_json_parse(read(f50, String))
      r60 = Sparlectra.scf_json_parse(read(f60, String))
      @test (r50["sparlectra"]["meta"]["f_nom"], r60["sparlectra"]["meta"]["f_nom"]) == (50.0, 60.0)
      c50 = [row["c1"] for row in get(r50["data"], "line", Any[])]
      c60 = [row["c1"] for row in get(r60["data"], "line", Any[])]
      tag == "zeng" && @test count(>(0.0), c50) == 6
      @test all(isapprox.(c60, c50 .* (5 / 6); rtol = 1e-14))
      back50 = importSCF(f50)
      back60 = importSCF(f60)
      @test back60.nominal_frequency_hz == 60.0
      @test Sparlectra.createYBUS(net = back60, sparse = true) == Sparlectra.createYBUS(net = orig, sparse = true)
      @test runpf!(back60, 30, 1e-10, 0) == runpf!(back50, 30, 1e-10, 0)
      @test vstate(back60) == vstate(back50)
    end
    root60 = Sparlectra.scf_json_parse(read(f60, String))

    # the file outranks the key: a configured 60 does not move a 50 Hz file,
    # and a key nobody set explicitly gives no warning
    at50 = @test_logs min_level = Logging.Warn Sparlectra.build_net(read_scf_json(sp14); config = key_cfg(60.0))
    @test at50.nominal_frequency_hz == 50.0
    # a file with 60 and an explicitly set key 50: the file wins, with
    # exactly one warning that names both values
    case60 = read_scf_json(f60)
    won = @test_logs min_level = Logging.Warn (:warn, r"meta\.f_nom = 60\.0 Hz.*model\.nominal_frequency_hz = 50\.0 Hz") Sparlectra.build_net(case60; config = key_cfg(50.0; user_set = true))
    @test won.nominal_frequency_hz == 60.0
    # a file without f_nom takes the key, silently
    delete!(root60["sparlectra"]["meta"], "f_nom")
    nofreq = joinpath(d, "nofreq.scf.json")
    write(nofreq, Sparlectra.scf_json_string(root60))
    taken = @test_logs min_level = Logging.Warn Sparlectra.build_net(read_scf_json(nofreq); config = key_cfg(60.0; user_set = true))
    @test taken.nominal_frequency_hz == 60.0
    @test Sparlectra.build_net(read_scf_json(nofreq); config = key_cfg(50.0)).nominal_frequency_hz == 50.0

    # a value that is not a positive number is a named error on load
    for bad in (0.0, 0, -50.0, true, "50", nothing)
      broken = Sparlectra.scf_json_parse(read(sp14, String))
      broken["sparlectra"]["meta"]["f_nom"] = bad
      err = try
        Sparlectra.scf_to_net(broken)
        nothing
      catch e
        e
      end
      @test err isa ArgumentError
      @test err isa ArgumentError && occursin("meta.f_nom", err.msg)
    end

    # every shipped case states 50 Hz and imports with it under any key
    shipped = sort!(filter(f -> endswith(f, ".scf.json"), readdir(dirname(sp14); join = true)))
    @test length(shipped) >= 10
    for f in shipped
      root = Sparlectra.scf_json_parse(read(f, String))
      @test root["sparlectra"]["meta"]["f_nom"] == 50.0
      n50 = Sparlectra.build_net(read_scf_json(f); config = key_cfg(50.0))
      n60 = Sparlectra.build_net(read_scf_json(f); config = key_cfg(60.0))
      @test n50.nominal_frequency_hz == n60.nominal_frequency_hz == 50.0
      @test Sparlectra.createYBUS(net = n60, sparse = true) == Sparlectra.createYBUS(net = n50, sparse = true)
    end

    # a plain PGM file (no sparlectra block) follows the key: c1 > 0 gives
    # 6/5 times b at 60 Hz
    pgm = joinpath(d, "plain_c1.json")
    write(pgm, """
{"version": "1.0", "type": "input", "is_batch": false, "attributes": {},
 "data": {
  "node": [{"id": 1, "u_rated": 110000.0}, {"id": 2, "u_rated": 110000.0}],
  "line": [{"id": 3, "from_node": 1, "to_node": 2, "from_status": 1, "to_status": 1, "r1": 1.21, "x1": 9.68, "c1": 2.0e-7, "tan1": 0.0}],
  "source": [{"id": 4, "node": 1, "status": 1, "u_ref": 1.0, "u_ref_angle": 0.0}],
  "sym_load": [{"id": 5, "node": 2, "status": 1, "type": 0, "p_specified": 20000000.0, "q_specified": 5000000.0}]
 }
}
""")
    p50 = Sparlectra.build_net(read_scf_json(pgm); config = key_cfg(50.0))
    p60 = Sparlectra.build_net(read_scf_json(pgm); config = key_cfg(60.0; user_set = true))
    @test (p50.nominal_frequency_hz, p60.nominal_frequency_hz) == (50.0, 60.0)
    @test isapprox(only(p60.branchVec).b_pu / only(p50.branchVec).b_pu, 6 / 5; rtol = 1e-14)
    @test only(p60.branchVec).x_pu == only(p50.branchVec).x_pu

    # the export refuses a frequency that is not positive (the keyword, or
    # an adapter option such as DTFAdapterOptions(f_nom = 0.0)), instead
    # of writing c1 as NaN or negative
    for bad in (0.0, -60.0)
      err = try
        Sparlectra.net_to_scf(net60; f_nom = bad)
        nothing
      catch e
        e
      end
      @test err isa ArgumentError && occursin("f_nom", err.msg)
    end

    # strict PGM drops meta.f_nom and says so with the value
    @test_logs (:warn, r"system frequency \(f_nom = 60\.0 Hz\)") match_mode = :any exportSCF(net60; file = joinpath(d, "sp14_60.pgm.json"), strict_pgm = true)

    # MATPOWER carries no frequency: the key decides, the power flow is
    # bit identical at 50 and 60 Hz, and the conversion to SCF writes the
    # key value as meta.f_nom (c1 needs one)
    mfile = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m"))
    m50 = Sparlectra._import_sparlectra_net(mfile, nothing, key_cfg(50.0))
    m60 = Sparlectra._import_sparlectra_net(mfile, nothing, key_cfg(60.0))
    @test (m50.nominal_frequency_hz, m60.nominal_frequency_hz) == (50.0, 60.0)
    @test runpf!(m60, 30, 1e-10, 0) == runpf!(m50, 30, 1e-10, 0)
    @test vstate(m60) == vstate(m50)
    mpc = Sparlectra.MatpowerIO.read_case(mfile; legacy_compat = true)
    conv = Sparlectra.convert_case(MatpowerAdapter(), mpc, Sparlectra.matpower_adapter_options(key_cfg(60.0)))
    @test conv.sparlectra.meta["f_nom"] == 60.0
    @test Sparlectra.build_net(conv; config = key_cfg(60.0)).nominal_frequency_hz == 60.0
  end)() end

  @testset "controllers round-trip through the declarative schema" begin (function ()
    # FACTS is what PGM cannot express at all; the case format carries the
    # controllers in the control.controllers schema verbatim, and the
    # reader instantiates them through applyConfiguredControllers!
    net = _scf_test_net()
    addSeriesReactanceControl!(net; fromBus = "3", toBus = "4", p_target_mw = 35.0, x_min_pu = 0.003, x_max_pu = 0.02, name = "tcsc_3_4", priority = 2)
    restoreBaseImpedances!(net)
    d = mktempdir()
    f = exportSCF(net; file = joinpath(d, "ctrl.scf.json"))
    entry = Sparlectra.scf_json_parse(read(f, String))["sparlectra"]["components"]["controllers"]["tcsc_3_4"]
    # the entry uses the same type name and keyword names as the YAML form
    @test entry["type"] == "series_reactance"
    @test entry["from_bus"] == "3" && entry["to_bus"] == "4"
    @test entry["p_target_mw"] == 35.0
    @test entry["x_min_pu"] == 0.003 && entry["x_max_pu"] == 0.02
    # the outer-loop rank is written when set (rank 0 stays out of the file)
    @test entry["priority"] == 2
    back = importSCF(f)
    c1 = only(collect_outer_controllers(net))
    c2 = only(collect_outer_controllers(back))
    @test c2.name == c1.name
    @test c2.p_target_mw == c1.p_target_mw
    @test (c2.x_min_pu, c2.x_max_pu) == (c1.x_min_pu, c1.x_max_pu)
    @test c2.limit_mode == c1.limit_mode
    @test c2.priority == c1.priority == 2
    restoreBaseImpedances!(back)
    @test read(exportSCF(back; file = joinpath(d, "ctrl2.scf.json")), String) == read(f, String)

    # the control range of a tap controller travels with it: a
    # configuration narrows the range of the Moorau controller of the
    # shipped sp_case14, export writes tap_min/tap_max, import restores
    # them, and a controller without a range writes neither key
    tnet = importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")))
    tctrl = only(Sparlectra._tap_controllers(tnet))
    tf0 = exportSCF(tnet; file = joinpath(d, "tap0.scf.json"))
    c0 = Sparlectra.scf_json_parse(read(tf0, String))["sparlectra"]["components"]["controllers"]
    @test !any(v -> haskey(v, "tap_min") || haskey(v, "tap_max"), values(c0))
    range_entry = Dict{String,Any}("type" => "power_transformer", "name" => "tap_range", "trafo" => tctrl.trafo, "mode" => "voltage", "tap_min" => 1.0, "tap_max" => 1.05)
    @test applyConfiguredControllers!(tnet, ControlConfig(enabled = true, controllers = Any[range_entry])) == 0
    tf = exportSCF(tnet; file = joinpath(d, "tap.scf.json"))
    tentry = only(v for v in values(Sparlectra.scf_json_parse(read(tf, String))["sparlectra"]["components"]["controllers"]) if v["type"] == "power_transformer")
    @test (tentry["tap_min"], tentry["tap_max"]) == (1.0, 1.05)
    tback = only(Sparlectra._tap_controllers(importSCF(tf)))
    @test (tback.tap_min, tback.tap_max) == (1.0, 1.05)

    # Voltage-dependent control travels in extra.<machine> as its points,
    # its interpolation mode and its limits. Before
    # that, exporting a network with a Q(U) machine and reading it back lost
    # the control silently. sp_case14 is the base because warmup_casePST
    # carries active links, which the solver refuses next to voltage-
    # dependent injections.
    qnet = importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")))
    qu = QUController(make_characteristic([(0.95, 0.30), (1.00, 0.0), (1.05, -0.20)]; interpolation = :spline); qmin_MVAr = -50.0, qmax_MVAr = 50.0, sbase_MVA = qnet.baseMVA)
    pu = PUController(make_characteristic([(0.90, 0.10), (1.10, 0.10)]); pmin_MW = 0.0, pmax_MW = 20.0, sbase_MVA = qnet.baseMVA)
    addProsumer!(net = qnet, busName = "Ilmrode_110", type = "SYNCHRONOUSMACHINE", p = 10.0, q = 0.0, qu_controller = qu, pu_controller = pu)
    @test validate!(net = qnet)[1]
    qf = exportSCF(qnet; file = joinpath(d, "qu.scf.json"))
    qroot = Sparlectra.scf_json_parse(read(qf, String))
    qentry = only(v for v in values(qroot["sparlectra"]["extra"]) if haskey(v, "qu_control"))
    @test qentry["qu_control"]["points"] == [[0.95, 0.30], [1.0, 0.0], [1.05, -0.20]]
    @test qentry["qu_control"]["interpolation"] == "spline"
    @test qentry["qu_control"]["qmin_mvar"] == -50.0 && qentry["qu_control"]["qmax_mvar"] == 50.0
    # linear is the default and is not written, like every other default
    @test !haskey(qentry["pu_control"], "interpolation")
    @test qentry["pu_control"]["pmin_mw"] == 0.0 && qentry["pu_control"]["pmax_mw"] == 20.0
    @test !haskey(qentry, "pq_gen_controller")
    qback = importSCF(qf)
    qps = only(filter(has_qu_controller, qback.prosumpsVec))
    @test qps.quController.characteristic.points == qu.characteristic.points
    @test qps.quController.characteristic.interpolation === :spline
    @test (qps.quController.qmin_pu, qps.quController.qmax_pu) == (qu.qmin_pu, qu.qmax_pu)
    @test qps.puController.characteristic.points == pu.characteristic.points
    @test (qps.puController.pmin_pu, qps.puController.pmax_pu) == (pu.pmin_pu, pu.pmax_pu)
    runpf!(qnet, 30, 1e-10, 0)
    runpf!(qback, 30, 1e-10, 0)
    @test maximum(abs(getNodeVm(qnet.nodeVec[i]) - getNodeVm(qback.nodeVec[i])) for i in eachindex(qnet.nodeVec)) < 1e-9
    @test read(exportSCF(qback; file = joinpath(d, "qu2.scf.json")), String) == read(qf, String)

    # the file is validated when it is read, by component name, not when
    # it is solved
    text = read(qf, String)
    # the writer puts one point per line; the block is matched as a whole
    qu_points = r"\"points\": \[\s*\[0\.95, 0\.3\],\s*\[1\.0, 0\.0\],\s*\[1\.05, -0\.2\]\s*\]"
    @test occursin(qu_points, text)
    broken = (
      ("one point", replace(text, qu_points => "\"points\": [[1.0, 0.0]]")),
      ("unknown mode", replace(text, "\"interpolation\": \"spline\"" => "\"interpolation\": \"cubic\"")),
      ("limits crossed", replace(text, "\"qmin_mvar\": -50.0" => "\"qmin_mvar\": 60.0")),
      ("voltages not increasing", replace(text, qu_points => "\"points\": [[1.0, 0.3], [0.95, 0.0], [1.05, -0.2]]")),
    )
    for (label, content) in broken
      @test content != text
      bf = joinpath(d, "broken.scf.json")
      write(bf, content)
      @test_throws ArgumentError importSCF(bf)
      err = try; importSCF(bf); nothing; catch e; e; end
      @test occursin("SCF: ", sprint(showerror, err))
    end

    # The MATPOWER converter's constant controllers (a PQ-bus generator's
    # limits) keep their short form: the flag, not an object, so every
    # existing case file written from MATPOWER stays byte for byte. The
    # shipped MATPOWER cases in the checkout have no PQ-bus generator, so
    # the case is built here.
    write(joinpath(d, "case_pq.m"), """
function mpc = case_pq
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
1 3 0 0 0 0 1 1.0 0 110 1 1.1 0.9;
2 1 40 15 0 0 1 1.0 0 110 1 1.1 0.9;
3 1 30 10 0 0 1 1.0 0 110 1 1.1 0.9;
];
mpc.gen = [
1 100 0 300 -300 1.02 100 1 300 0;
2 20 5 30 -30 1.0 100 1 50 0;
];
mpc.branch = [
1 2 0.01 0.05 0.0 999 999 999 0 0 1 -360 360;
2 3 0.01 0.05 0.0 999 999 999 0 0 1 -360 360;
];
""")
    mnet = SparlectraApp._se_import_case_net(joinpath(d, "case_pq.m"), Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true))
    @test count(has_qu_controller, mnet.prosumpsVec) == 1 && count(has_pu_controller, mnet.prosumpsVec) == 1
    mf = exportSCF(mnet; file = joinpath(d, "pq.scf.json"))
    mextra = values(Sparlectra.scf_json_parse(read(mf, String))["sparlectra"]["extra"])
    @test count(v -> get(v, "pq_gen_controller", false) === true, mextra) == 1
    @test count(v -> haskey(v, "qu_control") || haskey(v, "pu_control"), mextra) == 0
    mback = importSCF(mf)
    @test count(has_qu_controller, mback.prosumpsVec) == 1
    @test read(exportSCF(mback; file = joinpath(d, "pq2.scf.json")), String) == read(mf, String)

    # two entries on one transformer with disjoint actuators (ratio for U,
    # phase for P) are two controllers: the existence check of the
    # declarative path took the second for a repeat of the first. sp_case6tap
    # holds U Kreuzau_110 and P over T1 with T1 and U Lindenau_20 with T3; a
    # second apply adds nothing, the outer loop settles inside the deadbands
    # and the controllers survive a write and read
    c6 = importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case6tap.scf.json")))
    @test length(Sparlectra._tap_controllers(c6)) == 3
    entries6 = Sparlectra.scf_json_parse(read(joinpath(dirname(@__DIR__), "data", "scf", "sp_case6tap.scf.json"), String))["sparlectra"]["components"]["controllers"]
    @test Sparlectra._scf_apply_controllers!(c6, entries6; branch_by_name = Dict("T1_Nordwerk_Kreuzau" => 5, "T3_Mittelbach_Lindenau" => 7)) == 0
    f6c = exportSCF(c6; file = joinpath(d, "case6tap.scf.json"))
    c6b = importSCF(f6c)
    @test length(Sparlectra._tap_controllers(c6b)) == 3
    @test read(exportSCF(c6b; file = joinpath(d, "case6tap2.scf.json")), String) == read(f6c, String)
    res6 = redirect_stdout(() -> run_control!(c6), devnull)
    @test res6.status === :converged
    row(n) = only(r for r in res6.controllers if r.controller_name == n)
    @test abs(row("T1_Nordwerk_Kreuzau OLTC").achieved_vm_pu - 1.0) <= 0.01
    @test abs(row("T3_Mittelbach_Lindenau OLTC").achieved_vm_pu - 1.0) <= 0.01
    @test abs(row("T1_Nordwerk_Kreuzau PST").achieved_p_mw - 60.0) <= 4.0
  end)() end

  @testset "power-grid-model interoperability" begin (function ()
    # A PGM `source` is the reference by definition: its u_ref is the slack
    # voltage whether or not a hand-written sparlectra block marks the
    # machine `regulated`. Found on a delivered file: a
    # hand-written feeder with u_ref 1.02 and an extra entry without the
    # flag solved with the slack at 1.0 pu, every other bus 0.02 pu low,
    # while the plain PGM twin of the same network solved at 1.02.
    let hand = mktempdir()
      text = """
{"version": "1.0", "type": "input", "is_batch": false, "attributes": {},
 "data": {
  "node": [{"id": 1, "u_rated": 110000.0}, {"id": 2, "u_rated": 110000.0}],
  "line": [{"id": 3, "from_node": 1, "to_node": 2, "from_status": 1, "to_status": 1, "r1": 1.21, "x1": 9.68, "c1": 0.0, "tan1": 0.0}],
  "source": [{"id": 4, "node": 1, "status": 1, "u_ref": 1.02, "u_ref_angle": 0.0}],
  "sym_load": [{"id": 5, "node": 2, "status": 1, "type": 0, "p_specified": 20000000.0, "q_specified": 5000000.0}]
 },
 "sparlectra": {"format_version": "1.0", "meta": {"case_name": "hand", "s_base": 100000000.0, "f_nom": 50.0, "intended_calculations": ["power_flow"]},
  "roles": {"slack": {"mode": "single", "nodes": [1]}},
  "extra": {"1": {"name": "B1"}, "2": {"name": "B2"}, "4": {"name": "Grid", "reference_pri": true}}}
}
"""
      hf = joinpath(hand, "hand.scf.json")
      write(hf, text)
      hnet = importSCF(hf)
      # the setpoint reaches the node when the bus types are resolved for
      # the solve; what the reader must get right is the machine itself
      hsrc = only(filter(p -> p.referencePri !== nothing, hnet.prosumpsVec))
      @test hsrc.vm_pu == 1.02 && hsrc.isRegulated
      runpf!(hnet, 30, 1e-10, 0)
      @test getNodeVm(hnet.nodeVec[1]) == 1.02
      # the same network without the sparlectra block is a plain PGM file and
      # has always solved at u_ref; both must agree
      plain = replace(text, r",\s*\"sparlectra\": \{.*\}\}\s*\}\s*$"s => "}")
      pf = joinpath(hand, "plain.json")
      write(pf, plain)
      pnet = importSCF(pf)
      runpf!(pnet, 30, 1e-10, 0)
      @test abs(getNodeVm(pnet.nodeVec[2]) - getNodeVm(hnet.nodeVec[2])) < 1e-12
    end
    # A file written by PGM itself: no namespaced block at all, `source`
    # instead of a slack flag, and a generic_branch (the component PGM users
    # have the fewest examples for).
    file = abspath(joinpath(dirname(@__DIR__), "data", "scf", "pgm_interop.json"))
    root = Sparlectra.scf_json_parse(read(file, String))
    @test !haskey(root, "sparlectra")
    net = importSCF(file)
    @test length(net.nodeVec) == 3
    # PGM has no slack flag: an in-service source IS the reference
    @test net.slackVec == [1]
    @test getNodeType(net.nodeVec[1]) === Sparlectra.Slack
    trafo = net.branchVec[2]
    @test trafo.ratio == 1.025                       # k is the from-side ratio
    @test isapprox(trafo.x_pu, 0.11 / (10.5e3^2 / 100.0e6); rtol = 1e-15)   # to-side referenced
    it, _ = runpf!(net, 40, 1e-12, 0)
    @test it <= 10
    # independent model of the same network: ideal ratio k on the from side,
    # impedance on the to side, constant-power load. No Sparlectra internals.
    zl = (15.0 + 45.0im) / (150.0e3^2 / 100.0e6)
    zt = (0.0 + 0.11im) / (10.5e3^2 / 100.0e6)
    k = 1.025
    s_load = (5.0 + 1.0im) / 100.0
    v3 = 1.0 + 0.0im
    v2 = 1.0 + 0.0im
    for _ in 1:5000
      il = conj(s_load / v3)
      v2 = 1.0 - zl * (il / k)
      v3 = v2 / k - zt * il
    end
    @test isapprox(getNodeVm(net.nodeVec[2]), abs(v2); atol = 1e-12)
    @test isapprox(getNodeVm(net.nodeVec[3]), abs(v3); atol = 1e-12)
    @test isapprox(net.nodeVec[3]._va_deg, rad2deg(angle(v3)); atol = 1e-10)
    # PGM's `source` IS Sparlectra's external network injection: the case
    # must come back as the same model, not as a generator (which is what
    # made the exported file a DIFFERENT network than the one read in)
    @test net.prosumpsVec[1].comp.cTyp === Sparlectra.ExternalNetworkInjection
    d0 = mktempdir()
    re = Sparlectra.scf_json_parse(read(exportSCF(net; file = joinpath(d0, "src.scf.json")), String))
    @test haskey(re["data"], "source") && length(re["data"]["source"]) == 1
    @test !haskey(re["data"], "sym_gen")
    @test re["data"]["source"][1]["u_ref"] == 1.0
    # a source we WRITE is complete in the same sense: sk and rx_ratio
    @test re["data"]["source"][1]["sk"] == 2.0e9
    @test re["data"]["source"][1]["rx_ratio"] == 0.1
    # A complete PGM source states its short-circuit power and R/X ratio:
    # PGM models it as a voltage source BEHIND that impedance, so the
    # numbers belong to the model. They arrive as declared case data.
    feeder = only(net.sc_sources.external_network_injections)
    @test String(feeder.bus) == Sparlectra._scf_bus_name(net, 1)
    @test isapprox(sqrt(3.0) * 150.0 * feeder.maxInitialSymShCCurrent_A / 1000.0, 2000.0; rtol = 1e-9)
    @test feeder.maxR1ToX1Ratio == 0.1
    @test Sparlectra._declared_slack_feeder(net, net.sc_sources).sk_MVA ≈ 2000.0
    # ... so the case is short-circuit ready as it stands
    sc = runShortCircuit!(importSCF(file); case = :max)
    @test count(r -> r.status === :ok, sc.rows) == 3
    # ... and power_flow.external_grid computes with the FILE's numbers
    cfg_eg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true,
      overrides = Dict{String,Any}("power_flow" => Dict{String,Any}("external_grid" => Dict{String,Any}("enabled" => true, "source" => "auto"))))
    nonideal = importSCF(file)
    note = Sparlectra._apply_external_grid_config!(nonideal, cfg_eg.powerflow; declared = Sparlectra._declared_slack_feeder(nonideal, nonideal.sc_sources))
    @test occursin("declared by the case data", note)
    @test occursin("2000.0 MVA", note)
    runpf!(nonideal, 40, 1e-12, 0)
    # the same independent model, now with the source impedance in front
    zs = 100.0 / 2000.0
    xs = zs / sqrt(1.0 + 0.1^2)
    rs = 0.1 * xs
    w3 = 1.0 + 0.0im
    w1 = 1.0 + 0.0im
    for _ in 1:5000
      il = conj(s_load / w3)
      w1 = 1.0 - (rs + xs * im) * (il / k)
      w3 = (w1 - zl * (il / k)) / k - zt * il
    end
    @test isapprox(getNodeVm(nonideal.nodeVec[1]), abs(w1); atol = 1e-10)
    @test isapprox(getNodeVm(nonideal.nodeVec[3]), abs(w3); atol = 1e-10)
    # the ideal reading (the default) is the one WITHOUT that impedance
    @test getNodeVm(net.nodeVec[1]) == 1.0

    # exporting it adds the namespaced block; the result still solves the same
    d = mktempdir()
    again = importSCF(exportSCF(net; file = joinpath(d, "pgm_again.scf.json"), source_format = "pgm"))
    it2, _ = runpf!(again, 40, 1e-12, 0)
    @test it2 == it
    @test isapprox(getNodeVm(again.nodeVec[3]), getNodeVm(net.nodeVec[3]); atol = 1e-12)

    # PGM's own JSON stores enums as integers (MeasuredTerminalType node = 9,
    # FaultType three_phase = 0); the reader aborted with a MethodError on
    # String(::Int). The integer form must read as the same case as the name.
    case14 = Sparlectra.scf_json_parse(read(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"), String))
    terminals = Sparlectra._SCF_PGM_ENUMS["measured_terminal_type"]
    as_int = deepcopy(case14)
    for s in as_int["data"]["sym_power_sensor"]
      s["measured_terminal_type"] = findfirst(==(s["measured_terminal_type"]), terminals) - 1
    end
    as_int["data"]["fault"][1]["fault_type"] = 0
    mkey(m) = (m.typ, m.value, m.sigma, m.busIdx, m.branchIdx, m.direction, m.id)
    @test map(mkey, Sparlectra.scf_to_net(as_int).measurements) == map(mkey, Sparlectra.scf_to_net(case14).measurements)
    # a value outside the enum is refused with the field name
    off = deepcopy(case14)
    off["data"]["sym_power_sensor"][1]["measured_terminal_type"] = 42
    @test_throws ArgumentError Sparlectra.scf_to_net(off)
    @test_throws r"measured_terminal_type.*got 42" Sparlectra.scf_to_net(off)
  end)() end

  @testset "strict PGM writer mode" begin (function ()
    net = _scf_test_net()
    d = mktempdir()
    it0, _ = runpf!(net, 30, 1e-10, 0)
    strict = @test_logs (:warn,) match_mode = :any exportSCF(net; file = joinpath(d, "strict.json"), strict_pgm = true)
    root = Sparlectra.scf_json_parse(read(strict, String))
    # the plain dataset: no namespaced block at all
    @test !haskey(root, "sparlectra")
    @test sort(collect(keys(root))) == ["attributes", "data", "is_batch", "type", "version"]
    # PGM knows no slack flag, so a slack generator becomes a source, and
    # its sym_gen/voltage_regulator rows are gone (or the node would inject twice)
    @test haskey(root["data"], "source") && !isempty(root["data"]["source"])
    src_nodes = Set(Int(r["node"]) for r in root["data"]["source"])
    @test all(!(Int(r["node"]) in src_nodes) for r in get(root["data"], "sym_gen", []))
    # ... and the file still solves, to the same voltages
    back = importSCF(strict)
    it, _ = runpf!(back, 30, 1e-10, 0)
    @test it == it0
    @test maximum(abs(getNodeVm(net.nodeVec[i]) - getNodeVm(back.nodeVec[i])) for i in eachindex(net.nodeVec)) < 1e-12
    # the full form is unaffected and still carries everything
    full = Sparlectra.scf_json_parse(read(exportSCF(net; file = joinpath(d, "full.scf.json")), String))
    @test haskey(full, "sparlectra") && haskey(full["sparlectra"], "extra")
    @test filesize(joinpath(d, "strict.json")) < filesize(joinpath(d, "full.scf.json"))

    # a voltage-dependent controller has no PGM counterpart; the strict
    # writer names it among what the plain dataset does not carry
    qnet = importSCF(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14_qu.scf.json")))
    @test any(has_qu_controller, qnet.prosumpsVec)
    @test_logs (:warn, r"voltage-dependent Q\(U\)/P\(U\) control") match_mode = :any exportSCF(qnet; file = joinpath(d, "strict_qu.json"), strict_pgm = true)

    # an unregulated external network injection (a CGMES boundary or area
    # equivalent, built the way the CGMES importer builds it) is a PQ
    # injection, not a voltage source: the full file keeps it unregulated
    # and the solution unchanged, the strict file has one source only.
    # Before 0.31.3 both wrote it as a PGM source: it came back regulated
    # at 1.0 pu (CGMES MiniGrid: 0.039 pu off), and a strict file had one
    # reference per injection
    eni = Net(name = "scf_unregulated_eni", baseMVA = 100.0)
    addBus!(net = eni, busName = "A", vn_kV = 110.0)
    addBus!(net = eni, busName = "B", vn_kV = 110.0)
    addPIModelACLine!(net = eni, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.1, b_pu = 0.0, status = 1)
    addProsumer!(net = eni, busName = "A", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "A")
    addProsumer!(net = eni, busName = "B", type = "EXTERNALNETWORKINJECTION", p = -5.0, q = -1.0, isRegulated = false)
    eni.prosumpsVec[end].comp.cTyp = Sparlectra.ExternalNetworkInjection
    addProsumer!(net = eni, busName = "B", type = "ENERGYCONSUMER", p = 20.0, q = 4.0)
    runpf!(eni, 30, 1e-10, 0)
    vm_b = getNodeVm(eni.nodeVec[2])
    back_eni = importSCF(exportSCF(eni; file = joinpath(d, "eni.scf.json")))
    kept = only(ps for ps in back_eni.prosumpsVec if Int(ps.comp.cFrom_bus) == 2 && ps.comp.cTyp == Sparlectra.ExternalNetworkInjection)
    @test kept.isRegulated === false
    runpf!(back_eni, 30, 1e-10, 0)
    @test abs(getNodeVm(back_eni.nodeVec[2]) - vm_b) < 1e-12
    strict_eni_file = @test_logs (:warn,) match_mode = :any exportSCF(eni; file = joinpath(d, "eni.pgm.json"), strict_pgm = true)
    strict_eni = Sparlectra.scf_json_parse(read(strict_eni_file, String))
    @test length(strict_eni["data"]["source"]) == 1
  end)() end

  @testset "three-winding transformer from a nameplate" begin (function ()
    # the convenience form: one block in CGMES PowerTransformerEnd terms
    # instead of three generic_branch legs plus a star node
    d = mktempdir()
    plate = Dict{String,Any}("ends" => Any[
      Dict{String,Any}("role" => "hv", "rated_u" => 220000.0, "rated_s" => 1.0e8, "r" => 0.6, "x" => 30.0, "b" => 0.0),
      Dict{String,Any}("role" => "mv", "rated_u" => 110000.0, "rated_s" => 1.0e8, "r" => 0.4, "x" => 15.0, "b" => 0.0),
      Dict{String,Any}("role" => "lv", "rated_u" => 10500.0, "rated_s" => 4.0e7, "r" => 0.2, "x" => 8.0, "b" => 0.0),
    ])
    root = Dict{String,Any}(
      "version" => Sparlectra.SCF_PGM_VERSION, "type" => "input", "is_batch" => false, "attributes" => Dict{String,Any}(),
      "data" => Dict{String,Any}(
        "node" => Any[Dict{String,Any}("id" => 1, "u_rated" => 220000.0), Dict{String,Any}("id" => 2, "u_rated" => 110000.0), Dict{String,Any}("id" => 3, "u_rated" => 10500.0)],
        "source" => Any[Dict{String,Any}("id" => 4, "node" => 1, "status" => 1, "u_ref" => 1.0)],
        "sym_load" => Any[Dict{String,Any}("id" => 5, "node" => 2, "status" => 1, "type" => 0, "p_specified" => 3.0e7, "q_specified" => 8.0e6)],
      ),
      "sparlectra" => Dict{String,Any}(
        "format_version" => Sparlectra.SCF_FORMAT_VERSION,
        "meta" => Dict{String,Any}("case_name" => "t3w", "s_base" => 1.0e8, "f_nom" => 50.0),
        "transformer_types" => Dict{String,Any}("cgmes_220_110_10" => plate),
        "components" => Dict{String,Any}("transformer3w" => Any[Dict{String,Any}("id" => 7, "type" => "cgmes_220_110_10",
          "nodes" => Any[Dict{String,Any}("role" => "hv", "node" => 1), Dict{String,Any}("role" => "mv", "node" => 2), Dict{String,Any}("role" => "lv", "node" => 3)])]),
        "extra" => Dict{String,Any}("1" => Dict{String,Any}("name" => "HV", "bus_index" => 1), "2" => Dict{String,Any}("name" => "MV", "bus_index" => 2),
          "3" => Dict{String,Any}("name" => "LV", "bus_index" => 3), "4" => Dict{String,Any}("name" => "Feeder", "element_index" => 1, "prosumption_type" => "ExternalNetworkInjection"),
          "5" => Dict{String,Any}("name" => "LoadMV", "element_index" => 2, "prosumption_type" => "EnergyConsumer")),
        "roles" => Dict{String,Any}("slack" => Dict{String,Any}("nodes" => Any[1])),
      ),
    )
    file = joinpath(d, "t3w.json")
    write(file, Sparlectra.scf_json_string(root))
    net = importSCF(file)
    # the star equivalent: an auxiliary node plus three legs, built through
    # the SAME constructor a hand-written network would use
    @test length(net.nodeVec) == 4
    @test length(net.branchVec) == 3
    it, _ = runpf!(net, 40, 1e-10, 0)
    @test it <= 10
    # an export writes the explicit legs again (the nameplate is an INPUT
    # shorthand, not a second representation), and that file runs the same
    again = importSCF(exportSCF(net; file = joinpath(d, "t3w_out.scf.json")))
    it2, _ = runpf!(again, 40, 1e-10, 0)
    @test it2 == it
    @test maximum(abs(getNodeVm(net.nodeVec[i]) - getNodeVm(again.nodeVec[i])) for i in eachindex(net.nodeVec)) == 0.0
    # the two forms are exclusive, and a broken nameplate fails on load
    bad = deepcopy(root)
    bad["sparlectra"]["components"]["transformer3w"][1]["ends"] = Any[]
    @test_throws ArgumentError Sparlectra.scf_to_net(bad)
    unknown = deepcopy(root)
    unknown["sparlectra"]["components"]["transformer3w"][1]["type"] = "does_not_exist"
    @test_throws ArgumentError Sparlectra.scf_to_net(unknown)
    missing_role = deepcopy(root)
    missing_role["sparlectra"]["transformer_types"]["cgmes_220_110_10"]["ends"][3]["role"] = "hv"
    @test_throws ArgumentError Sparlectra.scf_to_net(missing_role)
  end)() end

  @testset "shunt state and PGM fault rows" begin (function ()
    d = mktempdir()
    # suite migration: the shipped sp_case14 carries a shunt too, so the
    # shunt-state export runs off the bundle without any download
    net = Sparlectra.importSCF(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"))
    net.shuntVec[1].status = 0
    setShuntEstimation!(net; busName = Sparlectra._scf_bus_name(net, net.shuntVec[1].busIdx), enabled = true)
    ids = Sparlectra.scf_id_map(net)
    sc = Dict{String,Any}("case" => "max", "sweep" => "explicit", "buses" => Any[ids.node[2]])
    f = exportSCF(net; file = joinpath(d, "state.scf.json"), short_circuit = sc)
    root = Sparlectra.scf_json_parse(read(f, String))
    # what PGM's g1/b1 cannot say travels in the namespaced block, and only
    # for shunts that deviate from the defaults
    state = root["sparlectra"]["components"]["shunt_state"]
    @test length(state) == 1
    @test state[1]["status"] == 0 && state[1]["estimate"] == true
    # the bus selection also exists in PGM's own vocabulary
    fault = root["data"]["fault"]
    @test length(fault) == 1
    @test fault[1]["fault_type"] == "three_phase" && fault[1]["r_f"] == 0.0
    @test fault[1]["fault_object"] == ids.node[2]
    @test scf_fault_nodes(f) == [ids.node[2]]
    back = importSCF(f)
    @test back.shuntVec[1].status == 0
    @test back.shuntVec[1].estimate
    @test read(f, String) == read(exportSCF(back; file = joinpath(d, "state2.scf.json"), short_circuit = sc), String)
    # an unbalanced or non-bolted fault is refused rather than ignored
    bad = deepcopy(root)
    bad["data"]["fault"][1]["fault_type"] = "single_phase_to_ground"
    @test_throws ArgumentError Sparlectra.scf_to_net(bad)
    bad2 = deepcopy(root)
    bad2["data"]["fault"][1]["x_f"] = 0.5
    @test_throws ArgumentError Sparlectra.scf_to_net(bad2)
  end)() end

  @testset "studies run from the case file" begin (function ()
    # the point of carrying a study is that a run executes it: the case list
    # and the short-circuit study come out of the file, not out of the request
    net = _scf_test_net()
    d = mktempdir()
    cfg = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH
    # 0.40.0: the study is the case file's scenarios block, run through the
    # file_block source; a two-outage scenario runs as one event instead of
    # being refused
    base_file = exportSCF(net; file = joinpath(d, "cont.scf.json"))
    scase = read_scf_json(base_file)
    idx = ScenarioIndex(scase)
    ids = Sparlectra.scf_id_map(net)
    b1, b2 = ids.branch[1], ids.branch[2]
    n1, n2 = getCompName(net.branchVec[1].comp), getCompName(net.branchVec[2].comp)
    out(id) = PatchOp(op = :status, target = idx.kind_by_id[id], id = id, value = 0.0)
    set = ScenarioSet(scenarios = [
      Scenario(name = "loss of $(n1)", weight = 2.5, ops = [out(b1)]),
      Scenario(name = "loss of $(n2)", ops = [out(b2)]),
      Scenario(name = "common mode", ops = [out(b1), out(b2)]),
    ])
    scase.sparlectra.scenarios = scenario_set_dict(set)
    f = write_scf_json(scase, base_file)
    res = redirect_stdout(devnull) do
      SparlectraApp._run_contingency_service(f, cfg, joinpath(d, "run_explicit"), "scf_ct", "branch"; scenario_source = "file_block")
    end
    md = SparlectraApp.to_dict(res)["metadata"]
    @test SparlectraApp.to_dict(res)["status"] == "succeeded"
    @test md["contingency_cases_source"] == "scenario_file_block"
    @test md["contingency_cases"] == 3                       # the file's list, not every branch
    # N-1 names the system frequency like a power-flow run: the export wrote
    # the line charging as c1 in farad together with meta.f_nom
    @test md["frequency_source"] == "case_file"
    @test count(==("System frequency: " * md["frequency_line"]), eachline(joinpath(d, "run_explicit", "run.log"))) == 1
    csv_text = read(joinpath(d, "run_explicit", "contingency_n1.csv"), String)
    @test occursin("loss of $(n1)", csv_text)
    @test occursin("common mode", csv_text)

    # A case file that carries its own configuration is loaded through the
    # API config path, and the service entry points take AbstractString: a
    # config path arriving as a SubString (a split request field) must not
    # end the run with a MethodError, which is how an N-1 on an exported
    # case failed in the Web UI.
    with_cfg = exportSCF(net; file = joinpath(d, "withcfg.scf.json"))
    Sparlectra.write_case_config(with_cfg, Dict{String,Any}("power_flow.mode" => "auto"))
    sub_cfg = split(string(cfg, "|marker"), "|")[1]
    @test sub_cfg isa SubString{String}
    res_sub = redirect_stdout(devnull) do
      SparlectraApp._run_contingency_service(with_cfg, sub_cfg, joinpath(d, "run_subcfg"), "scf_sub", "branch")
    end
    d_sub = SparlectraApp.to_dict(res_sub)
    @test d_sub["status"] == "succeeded"
    @test d_sub["metadata"]["contingency_cases"] > 0
    # the case configuration file's settings reached the run (the file is
    # nested YAML with the case-scope header)
    cc_text = read(Sparlectra.case_config_path(with_cfg), String)
    @test occursin("scope: case", cc_text)
    @test occursin("mode: auto", cc_text)
    res_sub_sc = redirect_stdout(devnull) do
      SparlectraApp._run_short_circuit_service(with_cfg, sub_cfg, joinpath(d, "run_subsc"), "scf_subsc")
    end
    # no sources in a MATPOWER-derived case: the reason must be the DATA,
    # never a type error
    @test SparlectraApp.to_dict(res_sub_sc)["reason"] == "short_circuit_data_missing"

    # State estimation runs on a case file, and the measurements come WITH
    # the model: no separate CSV has to exist for it (that combination was
    # rejected outright before, with "needs a MATPOWER or CGMES case").
    se_net = _scf_test_net()
    Sparlectra.readMeasurementTable!(se_net; file = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.measurements.csv")), zero_injection = false)
    se_case = exportSCF(se_net; file = joinpath(d, "se_case.scf.json"), intended_calculations = String["power_flow", "state_estimation"])
    @test !isempty(importSCF(se_case).measurements)
    # the service's SE run advises on the topology of these small sets by
    # design; the advisory is captured, anything else that warns fails.
    # Tap estimation is on by default (state_estimation.tap_estimation.enabled)
    # and releases the PST of this case; its fixation notice warns when the
    # rise exceeds the 95 % quantile (five runs of this testset, J 877 to
    # 1066 before the fixation)
    se_quiet(f) = run_with_expected_warnings(() -> redirect_stdout(f, devnull), (r"topology precheck reported", r"runse!: tap fixation J"))
    res_se = se_quiet() do
      SparlectraApp._run_state_estimation_service(se_case, cfg, joinpath(d, "run_se"), "scf_se", joinpath(d, "no_such_set.csv"))
    end
    d_se = SparlectraApp.to_dict(res_se)
    @test d_se["status"] == "succeeded"
    @test d_se["metadata"]["se_measurement_rows"] == length(importSCF(se_case).measurements)
    # the run writes its own measurement artifact, so it stays reproducible
    @test isfile(joinpath(d, "run_se", "measurements.csv"))
    @test occursin("from the case file", read(joinpath(d, "run_se", "run.log"), String))
    # phase timing, same file the power-flow path writes: without it a slow
    # estimation can only be guessed at (the topology precheck once cost
    # 45 s of a 55 s run on a 25000-bus net and nothing in the run
    # directory showed it). The precheck needs its own line.
    se_perf = joinpath(d, "run_se", "performance.log")
    @test isfile(se_perf)
    se_perf_text = read(se_perf, String)
    @test occursin("topology_precheck:", se_perf_text)
    @test occursin("state_estimation:", se_perf_text)
    @test occursin("total_service:", se_perf_text)
    @test length(d_se["metadata"]["service_phase_timings"]) >= 6
    # The set's provenance travels WITH the case: whether the values carry
    # noise, the per-row truth values, and the tap deviations the generator
    # applied. Losing those made an SCF run silently drop the tap warning
    # and the measured-vs-truth delta file that a CSV run produces.
    prov_set = joinpath(d, "prov.measurements.csv")
    writeMeasurementTable(prov_set, Sparlectra._measurement_table_from_measurements(se_net, se_net.measurements, nameMeasurementKeys(se_net)).rows; comments = String[
      "case: se_case.scf.json", "generator: v2", "seed: 7",
      "noise: gaussian (sigma U 0.5%, P 1%, Q 1% of reading)",
      "sparlectra-taps v1",
      "branch,mrid,name,from_bus,to_bus,neutral_ratio,step,electrical_step,fixed_step,transferred_step,generation_deviation_steps",
      "2,,T1,B1,B2,1.0,0.01,1.5,2,0,2.0",
      string("truth_value,", se_net.measurements[1].id, ",", repr(se_net.measurements[1].value + 0.01)),
    ])
    prov = SparlectraApp._webui_measurement_set_provenance(prov_set)
    @test prov["noise"] == true
    @test prov["seed"] == 7
    @test length(prov["tap_deviations"]) == 1
    @test prov["tap_deviations"][1]["branch"] == 2
    @test prov["tap_deviations"][1]["steps"] == 2.0
    @test length(prov["truth_values"]) == 1
    prov_case = exportSCF(se_net; file = joinpath(d, "prov_case.scf.json"), source_reference = "warmup_casePST.m",
      intended_calculations = String["power_flow", "state_estimation"], measurement_provenance = prov)
    back = Sparlectra._scf_measurement_provenance(prov_case)
    @test back["noise"] == true
    @test back["tap_deviations"][1]["steps"] == 2.0
    @test length(back["truth_values"]) == 1
    @test Sparlectra._scf_source_reference(prov_case) == "warmup_casePST.m"
    res_prov = se_quiet() do
      SparlectraApp._run_state_estimation_service(prov_case, cfg, joinpath(d, "run_prov"), "scf_prov", joinpath(d, "no_such_set.csv"))
    end
    @test SparlectraApp.to_dict(res_prov)["status"] == "succeeded"
    # the deviation reaches the run, and the truth value reaches the deltas
    @test occursin("documents a tap deviation", read(joinpath(d, "run_prov", "run.log"), String))
    @test isfile(joinpath(d, "run_prov", "se_deltas.csv"))
    @test count(row -> first(row) == "meas", run_csv(joinpath(d, "run_prov", "se_deltas.csv")).rows) == 1
    # the headline leads with J/dof: a bare J grows with the row count and
    # reads as an alarm where J/dof says the set is healthy
    @test occursin("J/dof = ", SparlectraApp.to_dict(res_prov)["message"])

    # A set that measures the same quantity twice is the signature of a
    # corrupted file: it inflates J without any single residual looking
    # wrong, so the run has to name it instead of estimating from it.
    second_reading = [Sparlectra.Measurement(typ = m.typ, value = m.value * 1.05, sigma = m.sigma, busIdx = m.busIdx,
      branchIdx = m.branchIdx, direction = m.direction, linkIdx = m.linkIdx, id = m.id * "_b") for m in se_net.measurements]
    doubled = vcat(collect(se_net.measurements), second_reading)
    @test sum(values(Sparlectra._duplicate_measured_quantities(doubled))) == length(se_net.measurements)
    @test isempty(Sparlectra._duplicate_measured_quantities(se_net.measurements))
    dbl_net = _scf_test_net()
    append!(dbl_net.measurements, doubled)
    dbl_case = exportSCF(dbl_net; file = joinpath(d, "dbl.scf.json"), intended_calculations = String["power_flow", "state_estimation"])
    # a second reading at one location is legal (redundant transducers) and
    # must survive the file: keying sensors by location alone silently kept
    # only the last reading, which lost half the set without a word
    dbl_back = importSCF(dbl_case)
    @test length(dbl_back.measurements) == length(doubled)
    @test sort([m.id for m in dbl_back.measurements]) == sort([m.id for m in doubled])
    @test sort([m.value for m in dbl_back.measurements]) == sort([m.value for m in doubled])
    res_dbl = se_quiet() do
      SparlectraApp._run_state_estimation_service(dbl_case, cfg, joinpath(d, "run_dbl"), "scf_dbl", joinpath(d, "no_such_set.csv"))
    end
    d_dbl = SparlectraApp.to_dict(res_dbl)
    @test d_dbl["metadata"]["se_duplicate_rows"] == length(se_net.measurements)
    @test occursin("repeat an already measured quantity", d_dbl["message"])

    # a case file WITHOUT measurements says what is missing instead of
    # failing on a file that was never there
    plain_case = exportSCF(_scf_test_net(); file = joinpath(d, "plain_case.scf.json"))
    res_none = se_quiet() do
      SparlectraApp._run_state_estimation_service(plain_case, cfg, joinpath(d, "run_se_none"), "scf_se_none", joinpath(d, "no_such_set.csv"))
    end
    d_none = SparlectraApp.to_dict(res_none)
    @test d_none["reason"] == "invalid_measurements"
    @test occursin("carries no measurements", d_none["message"])

    # short circuit: a case file WITHOUT sources says so instead of
    # reporting a zero-current network
    f0 = exportSCF(net; file = joinpath(d, "nosrc.scf.json"), short_circuit = Dict{String,Any}("case" => "max"))
    res0 = redirect_stdout(devnull) do
      SparlectraApp._run_short_circuit_service(f0, cfg, joinpath(d, "run_sc0"), "scf_sc0")
    end
    @test SparlectraApp.to_dict(res0)["reason"] == "short_circuit_data_missing"

    # ... and one WITH sources runs the study the file defines
    scnet = Net(name = "scstudy", baseMVA = 100.0)
    addBus!(net = scnet, busName = "A", vn_kV = 110.0)
    addBus!(net = scnet, busName = "B", vn_kV = 110.0)
    addPIModelACLine!(net = scnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
    addProsumer!(net = scnet, busName = "A", type = "ExternalNetworkInjection", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
    addProsumer!(net = scnet, busName = "B", type = "ENERGYCONSUMER", p = 20.0, q = 5.0)
    addExternalGrid!(net = scnet, busName = "A", sk_max_MVA = 5000.0, sk_min_MVA = 3000.0, rx_max = 0.1, rx_min = 0.1)
    scids = Sparlectra.scf_id_map(scnet)
    study = Dict{String,Any}("case" => "min", "c_factor" => 0.95, "sweep" => "explicit", "buses" => Any[scids.node[2]])
    fsc = exportSCF(scnet; file = joinpath(d, "sc.scf.json"), short_circuit = study)
    # the source data survives the round trip, which is what makes the run possible
    @test length(importSCF(fsc).sc_sources.external_network_injections) == 1
    # the technical side of the format comparison below: `auto` follows the
    # locale, so this run pins the portable one (a German desktop would give
    # excel_de and the comparison would test excel_de against itself)
    res_sc = withenv("LC_ALL" => "C.UTF-8") do
      redirect_stdout(devnull) do
        SparlectraApp._run_short_circuit_service(fsc, cfg, joinpath(d, "run_sc"), "scf_sc")
      end
    end
    dsc = SparlectraApp.to_dict(res_sc)
    @test dsc["status"] == "succeeded"
    msc = dsc["metadata"]
    @test msc["input_format_detected"] == "scf"
    @test msc["sc_study_from_case_file"] == true
    @test msc["sc_case_selected"] == "min"                   # the file picks the headline case
    @test msc["sc_c_factor"] == 0.95                         # ... and the c-factor, config being at its default
    @test msc["sc_case_rows"] == 1                           # explicit sweep: only the listed bus
    @test msc["sc_worst_bus"] == "B"
    @test isfile(joinpath(d, "run_sc", "short_circuit_max.csv"))
    @test isfile(joinpath(d, "run_sc", "short_circuit_min.csv"))
    sclog = read(joinpath(d, "run_sc", "run.log"), String)
    @test occursin("from the case file", sclog)
    # issue #376 follow-up: the short-circuit CSVs follow output.csv_format
    # like every other CSV artifact of a run. Default (technical): comma
    # delimiter, dot decimal. excel_de (machine-scope key, so it comes from
    # the general configuration file): semicolon delimiter, decimal comma.
    sc_default_lines = readlines(joinpath(d, "run_sc", "short_circuit_max.csv"))
    @test sc_default_lines[1] == "bus,vn_kV,island,status,c,zk_ohm,rx_ratio,ik_kA,sk_MVA,kappa,ip_kA,flagged,reasons,lower_bound"
    @test !occursin(';', sc_default_lines[2])
    cfg_de = joinpath(d, "config_excel_de.yaml")
    cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, cfg_de)
    open(cfg_de, "a") do io
      println(io, "output:")
      println(io, "  csv_format: excel_de")
    end
    res_sc_de = redirect_stdout(devnull) do
      SparlectraApp._run_short_circuit_service(fsc, cfg_de, joinpath(d, "run_sc_de"), "scf_sc_de")
    end
    @test SparlectraApp.to_dict(res_sc_de)["status"] == "succeeded"
    for name in ("short_circuit_max.csv", "short_circuit_min.csv")
      sc_de_lines = readlines(joinpath(d, "run_sc_de", name))
      @test sc_de_lines[1] == "bus;vn_kV;island;status;c;zk_ohm;rx_ratio;ik_kA;sk_MVA;kappa;ip_kA;flagged;reasons;lower_bound"
      @test length(sc_de_lines) == length(sc_default_lines)
      # the c column (fifth field, 1.1 for max / 0.95 for min at HV) always
      # carries a fraction, unlike vn_kV, which a whole kV value would
      # print without any separator: dot decimal by default, decimal comma
      # under excel_de
      @test occursin('.', split(sc_default_lines[2], ',')[5])
      @test !occursin('.', split(sc_de_lines[2], ';')[5])
      @test occursin(',', split(sc_de_lines[2], ';')[5])
    end
    # issue #386: one writer for every result CSV. The SE state CSV follows
    # the format too and reads back in every format; the AC island report
    # and the generic writer as well.
    cfg_obj_de = Sparlectra.load_sparlectra_config(cfg_de; reload = true)
    @test Sparlectra.result_csv_format() == "auto"
    @test with_sparlectra_config(Sparlectra.result_csv_format, cfg_obj_de) == "excel_de"
    gen_path = joinpath(d, "generic.csv")
    Sparlectra.write_result_csv(gen_path, ("name", "value", "flag"), (("a;b", 1234.5, true), ("c", 2, false)); format = "excel_de")
    gen_lines = readlines(gen_path)
    @test gen_lines[1] == "name;value;flag"
    @test gen_lines[2] == "\"a;b\";1.234,5;true"
    withenv("LC_ALL" => "C.UTF-8") do
      Sparlectra.write_result_csv(gen_path, ("name", "value"), (("x", 0.5),))
    end
    @test readlines(gen_path) == ["name,value", "x,0.5"]
    # SE state round trip under excel_de (decimal comma, semicolon)
    scf_dir = joinpath(dirname(@__DIR__), "data", "scf")
    se_net = Sparlectra.importSCF(joinpath(scf_dir, "sp_case5.scf.json"))
    Sparlectra.readMeasurementTable!(se_net; file = joinpath(scf_dir, "sp_case5.measurements.csv"), zero_injection = false)
    se_res = run_with_expected_warnings(() -> runse!(se_net), (r"topology precheck reported",))
    @test se_res.converged
    se_path = joinpath(d, "se_state_de.csv")
    with_sparlectra_config(() -> writeSEStateCSV(se_net; file = se_path), cfg_obj_de)
    se_lines = readlines(se_path)
    @test se_lines[1] == "# sparlectra-se-state v1"
    @test se_lines[2] == "bus;vm_pu;va_deg;pinj_MW;qinj_MVar"
    @test occursin(',', split(se_lines[3], ';')[2])
    se_back = Sparlectra.importSCF(joinpath(scf_dir, "sp_case5.scf.json"))
    readSEStateCSV!(se_back; file = se_path)
    @test all(isapprox(se_back.nodeVec[i]._vm_pu, se_net.nodeVec[i]._vm_pu; atol = 1e-12) for i in eachindex(se_net.nodeVec))
    @test all(isapprox(se_back.nodeVec[i]._va_deg, se_net.nodeVec[i]._va_deg; atol = 1e-10) for i in eachindex(se_net.nodeVec))
    # every CSV artifact of one power-flow run carries the declared
    # delimiter (the AC island report and the mismatch history included)
    pf_de = start_powerflow_run(Dict{String,Any}("casefile" => fsc, "config_file" => cfg_de, "output_root" => joinpath(d, "pf_de_root")))
    @test pf_de["status"] == "succeeded"
    pf_csvs = filter(f -> endswith(f, ".csv"), readdir(pf_de["output_dir"]))
    @test !isempty(pf_csvs)
    for f in pf_csvs
      data_lines = filter(l -> !startswith(l, "#"), readlines(joinpath(pf_de["output_dir"], f)))
      first_line = isempty(data_lines) ? "" : data_lines[1]
      @test occursin(';', first_line)
    end
    # the Web UI form field (request detailed_result_csv_format) with a
    # technical configuration file: the request value is folded into the
    # run's configuration, so EVERY CSV of the run carries the semicolon,
    # not only the two detailed exports (a Web UI run)
    pf_req = start_powerflow_run(Dict{String,Any}("casefile" => fsc, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(d, "pf_req_root"), "detailed_result_csv" => true, "detailed_result_csv_format" => "excel_de"))
    @test pf_req["status"] == "succeeded"
    req_csvs = filter(f -> endswith(f, ".csv"), readdir(pf_req["output_dir"]))
    @test "bus_powers.csv" in req_csvs
    for f in req_csvs
      data_lines = filter(l -> !startswith(l, "#"), readlines(joinpath(pf_req["output_dir"], f)))
      @test occursin(';', isempty(data_lines) ? "" : data_lines[1])
    end
    @test occursin("csv_format: excel_de", read(joinpath(pf_req["output_dir"], "effective_config.yaml"), String))
    # the registry is the template again after the run
    @test Sparlectra.result_csv_format() == "auto"
    # the same (older) request field reaches the state-estimation service
    # through the one override output.csv_format (a Web UI run kept
    # the comma in every SE artifact); the services themselves take no
    # request-level format since 0.16.0
    se_req = start_powerflow_run(Dict{String,Any}("casefile" => joinpath(scf_dir, "sp_case5.scf.json"), "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(d, "se_req_root"), "se_mode" => true, "measurement_file" => joinpath(scf_dir, "sp_case5.measurements.csv"), "detailed_result_csv_format" => "excel_de"))
    @test se_req["status"] == "succeeded"
    se_dir = se_req["output_dir"]
    # the estimator names the frequency too; sp_case5 has c1 = 0 on every line
    se_md = Sparlectra.scf_json_parse(read(joinpath(se_dir, "result.json"), String))["metadata"]
    @test se_md["frequency_source"] == "not_used"
    @test count(==("System frequency: not used (no capacitance in farad); case file states 50 Hz"), eachline(joinpath(se_dir, "run.log"))) == 1
    for f in ("measurements_verdict.csv", "se_state.csv", "se_deltas.csv")
      isfile(joinpath(se_dir, f)) || continue
      se_lines = filter(l -> !startswith(l, "#"), readlines(joinpath(se_dir, f)))
      @test occursin(';', isempty(se_lines) ? "" : se_lines[1])
    end
    # an unknown node id in the sweep fails on the run, not with an empty table
    root = Sparlectra.scf_json_parse(read(fsc, String))
    root["sparlectra"]["short_circuit"]["buses"] = Any[999999]
    write(joinpath(d, "badbus.scf.json"), Sparlectra.scf_json_string(root))
    res_bad = redirect_stdout(devnull) do
      SparlectraApp._run_short_circuit_service(joinpath(d, "badbus.scf.json"), cfg, joinpath(d, "run_scbad"), "scf_scbad")
    end
    @test SparlectraApp.to_dict(res_bad)["status"] == "failed"

    # the Web UI offers the button for a case file carrying sources
    @test SparlectraApp._webui_case_has_short_circuit_data(fsc)
    @test !SparlectraApp._webui_case_has_short_circuit_data(f0)

    # a shipped case whose machines carry no source data: the feeder has
    # data, the machine counts on the default, so the run is a warning
    # that counts; the substituted defaults stand in run.log and none of
    # them reaches the console (an empty list fails on any warning)
    res_mixed = run_with_expected_warnings(String[]) do
      redirect_stdout(devnull) do
        SparlectraApp._run_short_circuit_service(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json")), cfg, joinpath(d, "run_sc_mixed"), "scf_sc_mixed")
      end
    end
    dmixed = SparlectraApp.to_dict(res_mixed)
    @test (dmixed["status"], dmixed["reason"]) == ("warning", "short_circuit_partial_defaults")
    @test (dmixed["metadata"]["sc_sources"], dmixed["metadata"]["sc_sources_with_data"]) == (2, 1)
    mixed_log = read(joinpath(d, "run_sc_mixed", "run.log"), String)
    @test occursin("Substituted defaults:", mixed_log)
    # the machine states its maximum active power, the default stands on it
    @test occursin("default 0.2 pu on its maximum active power 20.0 MW", mixed_log)
  end)() end

  @testset "study definitions are carried and validated" begin (function ()
    net = _scf_test_net()
    d = mktempdir()
    sc = Dict{String,Any}("case" => "max", "sweep" => "all_buses")
    f = exportSCF(net; file = joinpath(d, "study.scf.json"), short_circuit = sc)
    studies = scf_case_studies(f)
    @test studies.short_circuit["case"] == "max"
    # a case without study blocks reports empty, never a guessed default
    plain = scf_case_studies(abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_casePST.scf.json")))
    @test isempty(plain.short_circuit)
    # broken study definitions fail on load, not at the start of a sweep
    root = Sparlectra.scf_json_parse(read(f, String))
    bad_case = deepcopy(root)
    bad_case["sparlectra"]["short_circuit"]["case"] = "bogus"
    @test_throws ArgumentError Sparlectra.scf_to_net(bad_case)
    bad_cf = deepcopy(root)
    bad_cf["sparlectra"]["short_circuit"]["c_factor"] = 3.0
    @test_throws ArgumentError Sparlectra.scf_to_net(bad_cf)
  end)() end

  @testset "SCF cases run through the framework and the service" begin (function ()
    fixture = abspath(joinpath(dirname(@__DIR__), "data", "scf", "sp_casePST.scf.json"))
    # the case format is a first-class case format: detected, selectable,
    # and runnable through the same paths as every other source
    @test Sparlectra._detect_case_format(fixture) === :scf
    # An explicit format is judged against the file CONTENT, not the
    # extension (`scf` saved for a MATPOWER case sent it into
    # the JSON reader, which blamed a JSON file nobody was reading). The
    # message names the file and the way out.
    mfile = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "sp_case9.m"))
    fmt_err = try
      Sparlectra._detect_case_format(mfile; requested = :scf)
      nothing
    catch err
      err
    end
    @test fmt_err isa ArgumentError
    @test occursin("sp_case9.m", sprint(showerror, fmt_err))
    @test occursin("auto", sprint(showerror, fmt_err))
    @test_throws ArgumentError Sparlectra._detect_case_format(fixture; requested = :matpower)
    @test Sparlectra._detect_case_format(mfile; requested = :matpower) === :matpower
    # cgmes and dtf keep the explicit choice without a content verdict
    @test Sparlectra._detect_case_format(mfile; requested = :cgmes) === :cgmes
    oddities = mktempdir()
    # a JSON case under a foreign name is still recognized by its content
    renamed = joinpath(oddities, "renamed_case.txt")
    cp(fixture, renamed)
    @test Sparlectra._detect_case_format(renamed) === :scf
    # a UTF-8 byte order mark (Windows editors write one) is skipped by the
    # reader and by the detection; the case then runs like the plain file
    withbom = joinpath(oddities, "bom.scf.json")
    write(withbom, vcat(UInt8[0xEF, 0xBB, 0xBF], read(fixture)))
    @test Sparlectra.scf_json_parse(read(withbom, String)) == Sparlectra.scf_json_parse(read(fixture, String))
    @test Sparlectra._detect_case_format(withbom; requested = :scf) === :scf
    @test SparlectraApp._webui_is_user_selectable_case("warmup_casePST.scf.json")
    @test !SparlectraApp._webui_is_user_selectable_case("run_metadata.json")
    cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
    res = run_sparlectra(casefile = fixture, config = cfg)
    @test res.final_converged
    res_bom = run_sparlectra(casefile = withbom, config = cfg)
    @test res_bom.final_converged
    @test res_bom.iterations == res.iterations
    root_dir = mktempdir()
    rs = start_powerflow_run(Dict{String,Any}("casefile" => fixture, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root_dir))
    @test rs["status"] == "succeeded"
    # the case file's own configuration applies, and an explicit override
    # still wins over it
    net = _scf_test_net()
    d = mktempdir()
    withcfg = exportSCF(net; file = joinpath(d, "cc.scf.json"))
    # the writer never emits a byte order mark
    @test first(read(withcfg)) == UInt8('{')
    Sparlectra.write_case_config(withcfg, Dict{String,Any}("power_flow.max_iter" => 44))
    r1 = start_powerflow_run(Dict{String,Any}("casefile" => withcfg, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root_dir))
    @test occursin("max_iter: 44", read(joinpath(r1["output_dir"], "effective_config.yaml"), String))
    r2 = start_powerflow_run(Dict{String,Any}("casefile" => withcfg, "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => root_dir, "config_overrides" => Dict{String,Any}("power_flow.max_iter" => 55)))
    @test occursin("max_iter: 55", read(joinpath(r2["output_dir"], "effective_config.yaml"), String))
  end)() end

  @testset "Web UI export route" begin (function ()
    root = mktempdir()
    cases = joinpath(root, "cases")
    mkpath(cases)
    src = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m"))
    cp(src, joinpath(cases, "warmup_casePST.m"))
    cp(string(first(splitext(src)), ".measurements.csv"), joinpath(cases, "warmup_casePST.measurements.csv"))
    rt = (; case_directory = cases, config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, operation_log = SparlectraApp.webui_operation_log_path(root), startup_config_error = nothing, runner = SparlectraApp.start_powerflow_run)
    resp = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/export-scf", Dict{String,Any}("casefile" => "warmup_casePST.m", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "power_flow_mode" => "auto"); output_root = root, runtime = rt)
    out = joinpath(cases, "warmup_casePST.scf.json")
    @test occursin("Exported", string(resp))
    @test isfile(out)
    txt = read(out, String)
    # the form's configuration travels into the case configuration file
    # next to the export (the writer emits no config block any more);
    # measurements found next to the case are exported with the file
    @test !occursin("\"config\"", txt)
    # the export is an SCF case, so its file is <stem>.config.yaml
    export_cc = read(joinpath(cases, "warmup_casePST.config.yaml"), String)
    @test occursin("scope: case", export_cc)
    @test occursin("mode: auto", export_cc)
    @test occursin("sym_voltage_sensor", txt)
    # a MATPOWER source is exported in its own unit system (per unit)
    @test Sparlectra.scf_case_units(Sparlectra.read_scf_json(out)) === :pu
    @test occursin("state_estimation", txt)
    # an SCF case that declares units = pu is re-exported per unit with the
    # same pu values; before the fix it went out as SI and lost its pu blocks
    pu_before = Sparlectra.scf_json_parse(txt)["sparlectra"]["extra"]
    @test count(k -> haskey(pu_before[k], "pu"), collect(keys(pu_before))) > 0
    again = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/export-scf", Dict{String,Any}("casefile" => "warmup_casePST.scf.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH); output_root = root, runtime = rt)
    @test occursin("Exported", string(again))
    @test Sparlectra.scf_case_units(Sparlectra.read_scf_json(out)) === :pu
    pu_after = Sparlectra.scf_json_parse(read(out, String))["sparlectra"]["extra"]
    @test all(k -> !haskey(pu_before[k], "pu") || get(get(pu_after, k, Dict()), "pu", nothing) == pu_before[k]["pu"], collect(keys(pu_before)))
    # unknown case and empty selection fail loudly, without writing
    bad = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/export-scf", Dict{String,Any}("casefile" => "does_not_exist.m"); output_root = root, runtime = rt)
    @test occursin("not%20found", string(bad)) || occursin("not+found", string(bad))
    empty_sel = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/export-scf", Dict{String,Any}(); output_root = root, runtime = rt)
    @test occursin("Select%20a%20case", string(empty_sel)) || occursin("Select+a+case", string(empty_sel))
    # the plain PGM variant writes its own file
    pgm = run_with_expected_warnings((r"strict_pgm = true",)) do
      SparlectraApp.route_sparlectra_webui("POST", "/powerflow/export-scf", Dict{String,Any}("casefile" => "warmup_casePST.m", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "scf_strict_pgm" => "true"); output_root = root, runtime = rt)
    end
    @test isfile(joinpath(cases, "warmup_casePST.pgm.json"))
    @test !haskey(Sparlectra.scf_json_parse(read(joinpath(cases, "warmup_casePST.pgm.json"), String)), "sparlectra")

    # a CGMES case (a zip in the case directory) exports as well, through the
    # importer the power-flow run uses; before 0.31.3 the export read cases
    # with the file-format importer and stopped with "file extension .zip is
    # not supported". The shipped sp_case14 delivery, packed at test time
    zip14 = joinpath(cases, "sp_case14_cgmes.zip")
    ZipArchives.ZipWriter(zip14) do w
      for f in sort(readdir(cgmes_fixture_dir("sp_case14"); join = true))
        ZipArchives.zip_newfile(w, basename(f))
        write(w, read(f))
      end
    end
    cg = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/export-scf", Dict{String,Any}("casefile" => "sp_case14_cgmes.zip", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH); output_root = root, runtime = rt)
    @test occursin("Exported", string(cg))
    @test isfile(joinpath(cases, "sp_case14_cgmes.scf.json"))
    # a CGMES source states SI values and stays SI
    @test Sparlectra.scf_case_units(Sparlectra.read_scf_json(joinpath(cases, "sp_case14_cgmes.scf.json"))) === :si
    cg_pgm = run_with_expected_warnings((r"strict_pgm = true",)) do
      SparlectraApp.route_sparlectra_webui("POST", "/powerflow/export-scf", Dict{String,Any}("casefile" => "sp_case14_cgmes.zip", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "scf_strict_pgm" => "true"); output_root = root, runtime = rt)
    end
    @test occursin("Exported", string(cg_pgm))
    @test isfile(joinpath(cases, "sp_case14_cgmes.pgm.json"))

    # the exported file can be taken back out of the case directory: the
    # export names it in the redirect, the page offers it, and the route
    # serves it as a download
    @test occursin("&download=warmup_casePST.scf.json", string(resp))
    page = SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case?download=warmup_casePST.scf.json&import_message=Exported"; output_root = root, runtime = rt)
    page_html = String(copy(page.body))
    @test occursin("/powerflow/case/download?case=warmup_casePST.scf.json", page_html)
    dl = SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case/download?case=warmup_casePST.scf.json"; output_root = root, runtime = rt)
    @test dl.status == 200
    @test any(k == "Content-Disposition" && occursin("warmup_casePST.scf.json", v) for (k, v) in dl.headers)
    @test any(k == "Content-Type" && occursin("application/json", v) for (k, v) in dl.headers)
    @test String(copy(dl.body)) == read(joinpath(cases, "warmup_casePST.scf.json"), String)
    # a name that escapes the case directory, or names nothing, is refused
    for bad_name in ("..%2F..%2Fetc%2Fpasswd", "nope.json")
      refused = SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case/download?case=$(bad_name)"; output_root = root, runtime = rt)
      @test refused.status == 303
      @test isempty(refused.body)
    end

    # the plain PGM export must be findable again: it is a runnable case,
    # so it belongs in the selector, and it downloads and runs like any
    # other case (it was written but invisible, which is how it got lost)
    @test isfile(joinpath(cases, "warmup_casePST.pgm.json"))
    form_html = String(copy(SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case"; output_root = root, runtime = rt).body))
    @test occursin("warmup_casePST.pgm.json", form_html)
    @test SparlectraApp._webui_is_user_selectable_case("warmup_casePST.pgm.json")
    pgm_run = SparlectraApp.start_powerflow_run(Dict("casefile" => "warmup_casePST.pgm.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(root, "runs_pgm")); case_directory = cases)
    @test pgm_run["status"] == "succeeded"
    # the run form carries the FULL path back after saving case settings, so
    # the download has to accept one that points into the case directory
    full_path = joinpath(cases, "warmup_casePST.pgm.json")
    dl_full = SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case/download?case=$(SparlectraApp._webui_urlencode(full_path))"; output_root = root, runtime = rt)
    @test dl_full.status == 200
    @test String(copy(dl_full.body)) == read(full_path, String)
    # ... without opening a path outside it
    outside = SparlectraApp.route_sparlectra_webui("GET", "/powerflow/case/download?case=$(SparlectraApp._webui_urlencode("/etc/passwd"))"; output_root = root, runtime = rt)
    @test outside.status == 303
    @test isempty(outside.body)

    # The state-estimation page must never arm a FOREIGN measurement set.
    # It used to preselect the first file in the directory when no bound set
    # existed, so a case14 run started with case118's set and failed the
    # binding check on a set the user never picked.
    cp(abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.measurements.csv")), joinpath(cases, "case118.measurements.csv"); force = true)
    se_page = String(copy(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.pgm.json"; output_root = root, runtime = rt).body))
    # since 0.50.0 the picker offers only the sets bound to the case, so the
    # foreign set is not even listed and the page says no set belongs here
    @test !occursin("value=\"case118.measurements.csv\"", se_page)
    @test occursin("No measurement set in the case directory belongs to warmup_casePST.pgm.json", se_page)
    # a case file that CARRIES its measurements says so, with the number the
    # estimator actually reads, and offers them as the first choice
    carrying = _scf_test_net()
    Sparlectra.readMeasurementTable!(carrying; file = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.measurements.csv")), zero_injection = false)
    exportSCF(carrying; file = joinpath(cases, "carrying.scf.json"), intended_calculations = String["power_flow", "state_estimation"])
    carry_page = String(copy(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=carrying.scf.json"; output_root = root, runtime = rt).body))
    @test !occursin("case118.measurements.csv\" selected", carry_page)
    # generating a set for a case file works and binds it to THAT case, so
    # the page preselects it afterwards
    # every generation writes a NEW, time-stamped file; its name travels in
    # the redirect (measurement_file=...), which is how the page arms it
    written_set(response) = SparlectraApp._webui_urldecode(String(match(r"measurement_file=([^&#]+)", Dict(response.headers)["Location"]).captures[1]))
    gen_resp = SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "warmup_casePST.pgm.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH); output_root = root, runtime = rt)
    gen_name = written_set(gen_resp)
    @test startswith(gen_name, "warmup_casePST.pgm.measurements.")
    generated = joinpath(cases, gen_name)
    @test isfile(generated)
    @test any(line -> line == "# case: warmup_casePST.pgm.json", eachline(generated))
    after_page = String(copy(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=warmup_casePST.pgm.json"; output_root = root, runtime = rt).body))
    @test occursin("$(gen_name)\" selected", after_page)
    # ... and the run with it goes through
    se_run = SparlectraApp.start_powerflow_run(Dict("casefile" => "warmup_casePST.pgm.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(root, "runs_se_gen"), "se_mode" => true, "measurement_file" => gen_name); case_directory = cases)
    @test se_run["status"] == "succeeded"
    # the ONE CSV setting (0.15.1, a Web UI run: the bad-data extract kept
    # the comma): a state-estimation request WITHOUT any CSV field under a
    # configuration file that says excel_de writes every CSV of the run
    # with the semicolon; nothing in the request may fall back to a
    # request default
    cfg_de = joinpath(root, "configuration_excel_de.yaml")
    write(cfg_de, replace(read(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, String), r"^(\s*csv_format:)\s*\S+"m => s"\1 excel_de"))
    @test occursin("csv_format: excel_de", read(cfg_de, String))
    se_de = SparlectraApp.start_powerflow_run(Dict("casefile" => "warmup_casePST.pgm.json", "config_file" => cfg_de, "output_root" => joinpath(root, "runs_se_de"), "se_mode" => true, "measurement_file" => gen_name); case_directory = cases)
    @test se_de["status"] == "succeeded"
    # measurements.csv is the measurement set itself (a measurement table,
    # an input format with its own fixed layout), not a result artifact;
    # its verdict file is one
    se_de_csvs = filter(f -> endswith(f, ".csv") && f != "measurements.csv", readdir(se_de["output_dir"]))
    @test "se_state.csv" in se_de_csvs
    @test "measurements_verdict.csv" in se_de_csvs
    for f in se_de_csvs
      data_lines = filter(l -> !startswith(l, "#"), readlines(joinpath(se_de["output_dir"], f)))
      @test occursin(';', isempty(data_lines) ? "" : data_lines[1])
    end
    # exporting a case FILE again must not grow format suffixes
    @test_logs (:warn, r"strict_pgm = true") match_mode = :any SparlectraApp.route_sparlectra_webui("POST", "/powerflow/export-scf", Dict{String,Any}("casefile" => "warmup_casePST.scf.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "scf_strict_pgm" => "true"); output_root = root, runtime = rt)
    @test !isfile(joinpath(cases, "warmup_casePST.scf.pgm.json"))
    @test isfile(joinpath(cases, "warmup_casePST.pgm.json"))

    # A GENERATED set wins over the one a case file carries: generating is a
    # deliberate act, and the newer data. Preferring the carried set made a
    # freshly generated one invisible (the run reported the old row count).
    carry_gen = _scf_test_net()
    Sparlectra.readMeasurementTable!(carry_gen; file = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.measurements.csv")), zero_injection = false)
    exportSCF(carry_gen; file = joinpath(cases, "wins.scf.json"), intended_calculations = String["power_flow", "state_estimation"])
    before_gen = String(copy(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=wins.scf.json"; output_root = root, runtime = rt).body))
    @test occursin("(from the case file: $(length(carry_gen.measurements)) rows)", before_gen)
    wins_name = written_set(SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/generate-measurements", Dict{String,Any}("casefile" => "wins.scf.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH); output_root = root, runtime = rt))
    @test startswith(wins_name, "wins.scf.measurements.")
    after_gen = String(copy(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=wins.scf.json"; output_root = root, runtime = rt).body))
    @test occursin("$(wins_name)\" selected", after_gen)
    # ... and the generated set REPLACES the carried rows instead of being
    # appended to them (59 carried + 75 generated = 134 rows in a file that
    # should have had 75)
    probe = importSCF(joinpath(cases, "wins.scf.json"))
    empty!(probe.measurements)
    Sparlectra.readMeasurementTable!(probe; file = joinpath(cases, wins_name), zero_injection = false)
    @test length(probe.measurements) < length(carry_gen.measurements) + 10
    # the set records WHETHER it carries noise, and the case file keeps that
    # statement: a noise-free set has J = 0 by construction
    @test any(l -> startswith(l, "# noise:"), eachline(joinpath(cases, wins_name)))
    ideal_net = _scf_test_net()
    runpf!(ideal_net, 30, 1e-10, 0)
    calcNetLosses!(ideal_net)
    setMeasurementsFromPF!(ideal_net; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false)
    ideal_case = exportSCF(ideal_net; file = joinpath(cases, "ideal.scf.json"), measurement_provenance = Dict{String,Any}("noise" => false, "generator" => "v2"))
    @test Sparlectra.scf_json_parse(read(ideal_case, String))["sparlectra"]["measurements"]["provenance"]["noise"] === false
    @test SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=ideal.scf.json"; output_root = root, runtime = rt).status == 200
    # noising needs NO power flow and keeps the sigmas; J moves from 0 to a
    # real value on the very same operating point
    noisy_name = written_set(SparlectraApp.route_sparlectra_webui("POST", "/stateestimation/add-noise", Dict{String,Any}("casefile" => "ideal.scf.json", "noise_seed" => "7"); output_root = root, runtime = rt))
    @test startswith(noisy_name, "ideal.scf.noisy.")
    noisy_set = joinpath(cases, noisy_name)
    @test isfile(noisy_set)
    @test any(l -> l == "# case: ideal.scf.json", eachline(noisy_set))
    ideal_run = SparlectraApp.start_powerflow_run(Dict("casefile" => "ideal.scf.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(root, "runs_ideal"), "se_mode" => true, "measurement_file" => ""); case_directory = cases)
    noisy_run = SparlectraApp.start_powerflow_run(Dict("casefile" => "ideal.scf.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(root, "runs_noisy"), "se_mode" => true, "measurement_file" => noisy_name); case_directory = cases)
    @test ideal_run["status"] == "succeeded" && noisy_run["status"] == "succeeded"
    @test occursin("J = 0.0", ideal_run["message"])
    @test !occursin("J = 0.0", noisy_run["message"])

    # a case file is selectable AND resolvable by its bare name: the service
    # resolver rejected .json outright, so an exported case could not be run
    resolved = SparlectraApp._resolve_powerflow_casefile("warmup_casePST.scf.json", cases)
    @test resolved == abspath(joinpath(cases, "warmup_casePST.scf.json"))
    @test_throws ArgumentError SparlectraApp._resolve_powerflow_casefile("no_such_case.scf.json", cases)
    run_scf = SparlectraApp.start_powerflow_run(Dict("casefile" => "warmup_casePST.scf.json", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "output_root" => joinpath(root, "runs_scf")); case_directory = cases)
    @test run_scf["status"] == "succeeded"
    @test run_scf["converged"] === true

    # uploading a case file: accepted by extension AND validated before it
    # is stored, so a foreign .json never lands in the case directory
    up = Dict{String,Any}("casefiles" => [SparlectraApp.WebUICaseUpload("uploaded.scf.json", Vector{UInt8}(read(joinpath(cases, "warmup_casePST.scf.json"))))])
    SparlectraApp.route_sparlectra_webui("POST", "/powerflow/import-cases", up; output_root = root, runtime = rt)
    @test isfile(joinpath(cases, "uploaded.scf.json"))
    foreign = Dict{String,Any}("casefiles" => [SparlectraApp.WebUICaseUpload("foreign.json", Vector{UInt8}(codeunits("{\"hello\": 1}")))])
    resp_foreign = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/import-cases", foreign; output_root = root, runtime = rt)
    @test !isfile(joinpath(cases, "foreign.json"))
    @test SparlectraApp._webui_scf_upload_reason(Vector{UInt8}(codeunits("{\"hello\": 1}"))) !== nothing
    @test SparlectraApp._webui_scf_upload_reason(Vector{UInt8}(read(joinpath(cases, "warmup_casePST.scf.json")))) === nothing

    # meta.f_nom is the frequency of the run configuration: 50 Hz without a
    # posted value, 60 Hz when the form posts 60 (MATPOWER states none)
    fnom(f) = Sparlectra.scf_json_parse(read(joinpath(cases, f), String))["sparlectra"]["meta"]["f_nom"]
    @test fnom("warmup_casePST.scf.json") == 50.0
    resp60 = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/export-scf", Dict{String,Any}("casefile" => "warmup_casePST.m", "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, "model_nominal_frequency_hz" => "60.0"); output_root = root, runtime = rt)
    @test occursin("Exported", string(resp60))
    @test fnom("warmup_casePST.scf.json") == 60.0
  end)() end

  return nothing
end

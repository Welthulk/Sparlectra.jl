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
    # the writer emits no config block any more; case settings live in the
    # case configuration file (write_case_config)
    @test !haskey(sp, "config")
    # start state is opt-in and never called a result
    @test !haskey(sp, "start_state")
    root2 = Sparlectra.net_to_scf(net; include_start_state = true)
    @test haskey(root2["sparlectra"], "start_state")
    @test root2["sparlectra"]["start_state"]["source"] == "solved_power_flow"
    @test !haskey(root2["sparlectra"], "results")
  end)() end

  @testset "measurements become PGM sensors" begin (function ()
    net = _scf_test_net()
    Sparlectra.readMeasurementTable!(net; file = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.measurements.csv")), zero_injection = false)
    root = Sparlectra.net_to_scf(net)
    data = root["data"]
    @test haskey(data, "sym_voltage_sensor")
    @test haskey(data, "sym_power_sensor")
    # P and Q of one location pair into a single PGM power sensor, so the
    # sensor count stays below the raw row count
    rows = root["sparlectra"]["measurements"]["rows"]
    nsensors = length(data["sym_voltage_sensor"]) + length(data["sym_power_sensor"]) + length(get(data, "sym_current_sensor", []))
    @test nsensors == length(rows)
    @test nsensors < length(net.measurements)
    # every Sparlectra row is accounted for exactly once
    allids = String[]
    for r in rows
      append!(allids, String[String(x) for x in r["sparlectra_ids"]])
    end
    @test sort(allids) == sort(String[m.id for m in net.measurements])
    # voltage sensors carry SI volts, power sensors watt
    v = first(data["sym_voltage_sensor"])
    @test v["u_measured"] > 1000.0
    p = first(data["sym_power_sensor"])
    @test p["measured_terminal_type"] in ("node", "branch_from", "branch_to")
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
    # against) are what `extra.name` carries; the internal component name
    # travels separately when it differs
    byname = Dict(e["name"] => e for e in values(sp["extra"]))
    @test haskey(byname, "AUX3WT_T1")
    @test byname["AUX3WT_T1"]["component_name"] != "AUX3WT_T1"
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
    # reader validation of the shunt-bay rows: an unknown type and an
    # undefined node are named, a row without a position goes to the end
    root = Sparlectra.scf_json_parse(read(fa, String))
    bay_rows = root["sparlectra"]["measurements"]["shunt_bay_rows"]
    bad_type = deepcopy(root)
    bad_type["sparlectra"]["measurements"]["shunt_bay_rows"][1]["type"] = "Pinj"
    @test_throws ArgumentError Sparlectra.scf_to_net(bad_type)
    bad_node = deepcopy(root)
    bad_node["sparlectra"]["measurements"]["shunt_bay_rows"][1]["node"] = 99999
    @test_throws ArgumentError Sparlectra.scf_to_net(bad_node)
    no_pos = deepcopy(root)
    delete!(no_pos["sparlectra"]["measurements"]["shunt_bay_rows"][1], "position")
    moved = Sparlectra.scf_to_net(no_pos).measurements
    @test length(moved) == length(full.measurements)
    @test (last(moved).id, last(moved).value) == (bay_rows[1]["id"], bay_rows[1]["value"])
    @test last(full.measurements).value != bay_rows[1]["value"]
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
    nf = @test_logs (:warn, r"no slot in the case file") match_mode = :any exportSCF(back; file = joinpath(d, "noslot.scf.json"), kw...)
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
    addSeriesReactanceControl!(net; fromBus = "3", toBus = "4", p_target_mw = 35.0, x_min_pu = 0.003, x_max_pu = 0.02, name = "tcsc_3_4")
    restoreBaseImpedances!(net)
    d = mktempdir()
    f = exportSCF(net; file = joinpath(d, "ctrl.scf.json"))
    entry = Sparlectra.scf_json_parse(read(f, String))["sparlectra"]["components"]["controllers"]["tcsc_3_4"]
    # the entry uses the same type name and keyword names as the YAML form
    @test entry["type"] == "series_reactance"
    @test entry["from_bus"] == "3" && entry["to_bus"] == "4"
    @test entry["p_target_mw"] == 35.0
    @test entry["x_min_pu"] == 0.003 && entry["x_max_pu"] == 0.02
    back = importSCF(f)
    c1 = only(collect_outer_controllers(net))
    c2 = only(collect_outer_controllers(back))
    @test c2.name == c1.name
    @test c2.p_target_mw == c1.p_target_mw
    @test (c2.x_min_pu, c2.x_max_pu) == (c1.x_min_pu, c1.x_max_pu)
    @test c2.limit_mode == c1.limit_mode
    restoreBaseImpedances!(back)
    @test read(exportSCF(back; file = joinpath(d, "ctrl2.scf.json")), String) == read(f, String)

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
    csv_text = read(joinpath(d, "run_explicit", "contingency_n1.csv"), String)
    @test occursin("loss of $(n1)", csv_text)
    @test occursin("common mode", csv_text)

    # a file that still carries the removed contingencies block fails with
    # the way out, on the kind path as well
    legacy = Sparlectra.scf_json_parse(read(f, String))
    legacy["sparlectra"]["contingencies"] = Dict{String,Any}("mode" => "explicit", "cases" => Any[Dict{String,Any}("name" => "old", "outages" => Any[Dict{String,Any}("component" => b1)])])
    fl = joinpath(d, "legacy.scf.json")
    write(fl, Sparlectra.scf_json_string(legacy))
    res_l = redirect_stdout(devnull) do
      SparlectraApp._run_contingency_service(fl, cfg, joinpath(d, "run_legacy"), "scf_ctl", "branch")
    end
    dl = SparlectraApp.to_dict(res_l)
    @test dl["status"] == "failed"
    @test occursin("sparlectra.scenarios", dl["message"])

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
    # design; the advisory is captured, anything else that warns fails
    se_quiet(f) = run_with_expected_warnings(() -> redirect_stdout(f, devnull), (r"topology precheck reported",))
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
    @test count(l -> startswith(l, "meas,"), readlines(joinpath(d, "run_prov", "se_deltas.csv"))) == 1
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
    res_sc = redirect_stdout(devnull) do
      SparlectraApp._run_short_circuit_service(fsc, cfg, joinpath(d, "run_sc"), "scf_sc")
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
    Sparlectra.write_result_csv(gen_path, ("name", "value"), (("x", 0.5),))
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
    for f in ("se_bad_data.csv", "se_state.csv", "se_deltas.csv")
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
    legacy = deepcopy(root)
    legacy["sparlectra"]["contingencies"] = Dict{String,Any}("mode" => "explicit", "cases" => Any[])
    @test_throws ArgumentError Sparlectra.scf_to_net(legacy)
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

  @testset "every output option is case-file and GUI editable" begin (function ()
    # the case file's config block uses the one existing allowlist, so the
    # logging surface has to be in it
    for key in ("output.console_summary", "output.console_diagnostics", "output.console_q_limit_events", "output.logfile_diagnostics", "output.logfile_performance", "output.logfile_warnings", "output.result_table_max_rows", "output.startup_latency_hint")
      @test key in Sparlectra.GUI_EDITABLE_CONFIG_KEYS
    end
    ok = Sparlectra.validate_gui_config_overrides(Dict{String,Any}("output.logfile_warnings" => "table", "output.console_summary" => false, "output.result_table_max_rows" => 50, "output.startup_latency_hint" => false))
    @test haskey(ok, "output")
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("output.logfile_warnings" => "bogus"))
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
    @test occursin("state_estimation", txt)
    # unknown case and empty selection fail loudly, without writing
    bad = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/export-scf", Dict{String,Any}("casefile" => "does_not_exist.m"); output_root = root, runtime = rt)
    @test occursin("not%20found", string(bad)) || occursin("not+found", string(bad))
    empty_sel = SparlectraApp.route_sparlectra_webui("POST", "/powerflow/export-scf", Dict{String,Any}(); output_root = root, runtime = rt)
    @test occursin("Select%20a%20case", string(empty_sel)) || occursin("Select+a+case", string(empty_sel))
    # the button is on the Case page for every case, not only
    # for cases that already carry saved settings
    @test occursin("/powerflow/export-scf", SparlectraApp.render_case_page())
    # the plain PGM variant has its own button and writes its own file
    @test occursin("scf_strict_pgm", SparlectraApp.render_case_page())
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
    @test occursin("Download selected case", page_html)
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
    @test occursin("carries <strong>$(length(carrying.measurements)) measurements</strong>", carry_page)
    @test occursin("(from the case file: $(length(carrying.measurements)) rows)", carry_page)
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
    # the ONE CSV setting (0.15.1, a Web UI run: se_bad_data.csv kept
    # the comma): a state-estimation request WITHOUT any CSV field under a
    # configuration file that says excel_de writes every CSV of the run
    # with the semicolon; nothing in the request may fall back to a
    # request default
    cfg_de = joinpath(root, "configuration_excel_de.yaml")
    write(cfg_de, replace(read(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, String), r"^(\s*csv_format:)\s*\S+"m => s"\1 excel_de"))
    @test occursin("csv_format: excel_de", read(cfg_de, String))
    se_de = SparlectraApp.start_powerflow_run(Dict("casefile" => "warmup_casePST.pgm.json", "config_file" => cfg_de, "output_root" => joinpath(root, "runs_se_de"), "se_mode" => true, "measurement_file" => gen_name); case_directory = cases)
    @test se_de["status"] == "succeeded"
    # measurements.csv is the measurement set itself (measurement CSV v1,
    # an input format with its own fixed header), not a result artifact
    se_de_csvs = filter(f -> endswith(f, ".csv") && f != "measurements.csv", readdir(se_de["output_dir"]))
    @test "se_state.csv" in se_de_csvs
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
    ideal_page = String(copy(SparlectraApp.route_sparlectra_webui("GET", "/powerflow?casefile=ideal.scf.json"; output_root = root, runtime = rt).body))
    @test occursin("ideal (noise-free)", ideal_page)
    @test occursin("/stateestimation/add-noise", ideal_page)
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

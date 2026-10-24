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
#
# file: test/test_cgmes_export.jl
# purpose: CGMES export identity and roundtrip contract — structural keys in
#          net.cgmes_ids, parallel-line numbering, deterministic uuid5
#          minting, byte-identical re-export, the duplicate-mRID guard, the
#          power-flow-identical self-roundtrip on EQ+TP+SSH+SV (PI lines,
#          transformers incl. phase shift and ratio-tap
#          machinery, machines, SVC, loads, shunts, bus links), SSH/SV
#          content, ZIP packaging and tool
#          provenance on nets built in memory, and the export-import-export
#          identity of the checked-in deliveries under data/cgmes_demo
#          (object-wise, with the fields the importer does not carry named).

using Test
using Sparlectra
using Dates
using UUIDs
using Logging
isdefined(@__MODULE__, :_cgmes_synth_eb_dir) || include(joinpath(@__DIR__, "cgmes_synth_eb_fixture.jl"))

# Fixed header stamp: with `created` pinned the exported files are
# byte-reproducible, which the determinism tests rely on.
const _EXPORT_STAMP = DateTime(2026, 1, 1, 12, 0, 0)

isdefined(@__MODULE__, :_cgmes_synth_tapside_dir) || include(joinpath(@__DIR__, "cgmes_synth3w_fixture.jl"))

# bare-call defaults are too tight for these fixtures — always solve through
# the explicit signature (maxIte, tol, verbose)
_solve!(net) = runpf!(net, 30, 1e-8, 0)

function _export_test_net()
  net = Net(name = "exportnet", baseMVA = 100.0)
  addBus!(net = net, busName = "A", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "B", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "C", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  # Two parallel lines on the A-B corridor; the second is declared B->A to
  # exercise the lexicographic bus-pair normalization of the structural key.
  addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.05, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B", toBus = "A", r_pu = 0.02, x_pu = 0.06, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B", toBus = "C", r_pu = 0.015, x_pu = 0.04, b_pu = 0.01, status = 1)
  return net
end

# every stage-2/3 class in one net: two voltage levels, transformer with
# off-nominal ratio + phase shift + ratio-tap machinery, PV machine, slack
# injection, SVC, load, shunt, bus link
function _roundtrip_net()
  net = Net(name = "rt", baseMVA = 100.0)
  addBus!(net = net, busName = "HV1", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "HV1B", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "HV2", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "MV1", vn_kV = 20.0, vm_pu = 1.0, va_deg = 0.0)
  addPIModelACLine!(net = net, fromBus = "HV1", toBus = "HV2", r_pu = 0.01, x_pu = 0.05, b_pu = 0.02, status = 1)
  # the two addACLine! lines of 0.42 (per unit, length 10 took r, x, b as
  # totals, length 0.5 halved them) as PI-model lines with the same values,
  # so the network stays the same after the addACLine! rebuild
  addPIModelACLine!(net = net, fromBus = "HV1B", toBus = "HV2", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
  addPIModelACLine!(net = net, fromBus = "HV1", toBus = "HV2", r_pu = 0.03 * 0.5, x_pu = 0.2 * 0.5, b_pu = 0.02 * 0.5, status = 1)
  addLink!(net = net, fromBus = "HV1", toBus = "HV1B", status = 1)
  addPIModelTrafo!(net = net, fromBus = "HV2", toBus = "MV1", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, status = 1, ratio = 1.04, shift_deg = 2.0)
  addProsumer!(net = net, busName = "HV1", type = "EXTERNALNETWORKINJECTION", referencePri = "HV1", vm_pu = 1.02, va_deg = 0.0, pMin = -30.0, pMax = 120.0)
  # the class marker of the importer, so the slack infeed leaves as
  # ExternalNetworkInjection with its asymmetric P band
  net.prosumpsVec[end].comp.cTyp = Sparlectra.ExternalNetworkInjection
  addProsumer!(net = net, busName = "HV2", type = "SYNCHRONOUSMACHINE", p = 15.0, q = 0.0, vm_pu = 1.01, isRegulated = true, qMin = -50.0, qMax = 50.0, pMin = 0.0, pMax = 80.0)
  addProsumer!(net = net, busName = "MV1", type = "ENERGYCONSUMER", p = 20.0, q = 5.0)
  addShuntMatpower!(net = net, busName = "MV1", Gs = 0.0, Bs = 4.0)
  addProsumer!(net = net, busName = "HV2", type = "STATICVARCOMPENSATOR", p = 0.0, q = 0.0, qMin = -121.0, qMax = 121.0)
  # the concrete class marker the CGMES importer maintains; addProsumer!
  # collapses the type string to the generic Generator component
  net.prosumpsVec[end].comp.cTyp = Sparlectra.StaticVarCompensator
  # ratio-tap machinery on the transformer winding (range export)
  net.trafos[1].side1.taps = Sparlectra.PowerTransformerTaps(Vn_kV = 110.0, step = 2, lowStep = -9, highStep = 9, neutralStep = 0, voltageIncrement_kV = 1.1)
  return net
end

function _compare_solved(n1, n2; atol = 1e-6)
  for (bus, i1) in n1.busDict
    i2 = n2.busDict[bus]
    @test isapprox(n1.nodeVec[i1]._vm_pu, n2.nodeVec[i2]._vm_pu; atol = atol)
    @test isapprox(n1.nodeVec[i1]._va_deg, n2.nodeVec[i2]._va_deg; atol = atol)
  end
end

# --- profile-file comparison for the export-import-export identity ----------
#
# The exporter writes objects in the order of the net it reads, and an
# imported net orders its buses by the topology walk, so a byte comparison
# of a re-export against the fixture fails on order alone. The comparison
# below is object-wise: every top-level RDF object (two-space indented
# `<cim:Class rdf:ID|about=...>` up to its close tag) is keyed by its mRID,
# the header lines are compared verbatim, and attribute values are compared
# numerically where both parse as numbers (the pu <-> physical unit
# conversions leave 1e-15 relative noise).
function _cgmes_profile_objects(text::AbstractString)
  header = String[]
  objects = Dict{String,Tuple{String,Vector{Pair{String,String}}}}()
  cls = ""
  id = ""
  attrs = Pair{String,String}[]
  for line in split(text, '\n')
    m = match(r"^  <(cim:[A-Za-z0-9]+) rdf:(?:ID|about)=\"#?_?([^\"]+)\">$", line)
    if m !== nothing
      cls = String(m.captures[1])
      id = String(m.captures[2])
      attrs = Pair{String,String}[]
    elseif !isempty(cls) && line == string("  </", cls, ">")
      haskey(objects, id) && error("duplicate object id in profile: ", id)
      objects[id] = (cls, attrs)
      cls = ""
    elseif !isempty(cls)
      a = match(r"^    <([^ >]+)(?: rdf:resource=\"([^\"]*)\"/>|>(.*)</[^>]+>)$", line)
      a === nothing && error("unparsed attribute line in ", cls, " ", id, ": ", line)
      push!(attrs, String(a.captures[1]) => String(something(a.captures[2], a.captures[3])))
    else
      # the exporter stamps its own version into the header comment and the
      # model description; a release bump must not read as a changed export
      # (the fixtures were written by 0.17.0, the test ran green on 0.17.0
      # only), so the stamp is neutralized before the verbatim comparison
      push!(header, replace(String(line), r"Sparlectra\.jl v\d+\.\d+\.\d+" => "Sparlectra.jl v<version>"))
    end
  end
  return header, objects
end

# SvPowerFlow rows are evaluated from the voltage state, so two solves that
# stop at the same tolerance agree on them only to that tolerance (1e-8
# MW-level residuals); every other value is model data and compares tight.
function _cgmes_values_equal(tag::AbstractString, a::AbstractString, b::AbstractString)::Bool
  a == b && return true
  fa = tryparse(Float64, a)
  fb = tryparse(Float64, b)
  (fa === nothing || fb === nothing) && return false
  tag in ("cim:SvPowerFlow.p", "cim:SvPowerFlow.q") && return isapprox(fa, fb; atol = 1e-6)
  return isapprox(fa, fb; rtol = 1e-9, atol = 1e-12)
end

# Compare one re-exported profile file against its fixture. `drop_classes`
# names the object classes the fixture carries and the re-export cannot;
# `drop_attr(cls, id, tag)` names the attributes that legitimately differ.
# One assertion per file: the mismatch list must be empty, and its first
# entries (class, mRID, attribute, fixture value, re-export value) land in
# the failure output. Returns the number of objects compared.
# `per_unit_ends` names PowerTransformerEnd mRIDs whose ratedU legitimately
# differs (it is then dropped by `drop_attr`): their g and b are in siemens
# on that ratedU, so they are compared as g * ratedU^2 and b * ratedU^2 of
# each side, the per-unit admittance, with the same tolerance.
function _cgmes_end_per_unit(attrs::Vector{Pair{String,String}})::Vector{Pair{String,String}}
  iu = findfirst(p -> p.first == "cim:PowerTransformerEnd.ratedU", attrs)
  iu === nothing && return attrs
  u2 = parse(Float64, attrs[iu].second)^2
  return [p.first in ("cim:PowerTransformerEnd.g", "cim:PowerTransformerEnd.b") ? (p.first => string(parse(Float64, p.second) * u2)) : p for p in attrs]
end

function _compare_cgmes_profile(fixture_text::AbstractString, reexport_text::AbstractString; drop_classes = (), drop_attr = (cls, id, tag) -> false, per_unit_ends = Set{String}())::Int
  # a checkout with CRLF conversion (Windows before the .gitattributes rule
  # covered XML) must compare equal to the LF re-export
  fixture_text = replace(String(fixture_text), "\r\n" => "\n")
  reexport_text = replace(String(reexport_text), "\r\n" => "\n")
  hf, of = _cgmes_profile_objects(fixture_text)
  hr, or = _cgmes_profile_objects(reexport_text)
  @test hf == hr
  expected_ids = Set(id for (id, (cls, _)) in of if !(cls in drop_classes))
  @test Set(keys(or)) == expected_ids
  mismatches = String[]
  compared = 0
  for id in intersect(expected_ids, Set(keys(or)))
    cf, af = of[id]
    cr, ar = or[id]
    if id in per_unit_ends && cf == cr == "cim:PowerTransformerEnd"
      af = _cgmes_end_per_unit(af)
      ar = _cgmes_end_per_unit(ar)
    end
    cf == cr || push!(mismatches, string(cf, " ", id, ": class ", cr, " in the re-export"))
    keep_f = [p for p in af if !drop_attr(cf, id, p.first)]
    keep_r = [p for p in ar if !drop_attr(cr, id, p.first)]
    if first.(keep_f) != first.(keep_r)
      push!(mismatches, string(cf, " ", id, ": attributes ", first.(keep_f), " vs ", first.(keep_r)))
    else
      for (pf, pr) in zip(keep_f, keep_r)
        _cgmes_values_equal(pf.first, pf.second, pr.second) || push!(mismatches, string(cf, " ", id, " ", pf.first, ": ", pf.second, " vs ", pr.second))
      end
    end
    compared += 1
  end
  @test (length(mismatches), first(mismatches, 8)) == (0, String[])
  return compared
end

function run_cgmes_export_tests()
  @testset "structural keys and parallel lines" begin (function ()
    net = _export_test_net()
    dir = mktempdir()
    files = writeCGMESFiles(net; path = dir, created = _EXPORT_STAMP)
    @test length(files) == 4
    @test all(isfile, files)
    @test endswith(files[3], "_SSH.xml")
    @test endswith(files[4], "_SV.xml")
    @test haskey(net.cgmes_ids, "MODEL|SSH")
    @test haskey(net.cgmes_ids, "MODEL|SV")
    # tool provenance: version + export date in every file, stamp = created
    for f in files
      head = read(f, String)
      @test occursin("Generated by Sparlectra.jl v$(Sparlectra.version()) on 2026-01-01T12:00:00Z", head)
      @test occursin("<md:Model.description>Sparlectra.jl v$(Sparlectra.version()) export</md:Model.description>", head)
    end
    @test haskey(net.cgmes_ids, "ACL|A|B|1")
    @test haskey(net.cgmes_ids, "ACL|A|B|2")
    @test haskey(net.cgmes_ids, "ACL|B|C|1")
    @test haskey(net.cgmes_ids, "TN|A")
    @test haskey(net.cgmes_ids, "BV|110.0")
    @test haskey(net.cgmes_ids, "MODEL|EQ")
    @test haskey(net.cgmes_ids, "MODEL|TP")
    # parallel lines share the corridor but never an mRID
    @test net.cgmes_ids["ACL|A|B|1"] != net.cgmes_ids["ACL|A|B|2"]
    eq = read(files[1], String)
    tp = read(files[2], String)
    for key in ("ACL|A|B|1", "ACL|A|B|2", "ACL|B|C|1")
      @test occursin("rdf:ID=\"_$(net.cgmes_ids[key])\"", eq)
    end
    @test occursin("rdf:ID=\"_$(net.cgmes_ids["TN|A"])\"", tp)
    # terminal keys follow the equipment sequence (T1 = from side)
    @test haskey(net.cgmes_ids, "ACL|A|B|1|T1")
    @test haskey(net.cgmes_ids, "ACL|A|B|1|T2")
    # the SV profile carries one voltage per bus
    sv = read(files[4], String)
    @test count("<cim:SvVoltage rdf:ID", sv) == 3
  end)() end

  @testset "minted ids are uuid5 over the key" begin (function ()
    net = _export_test_net()
    writeCGMESFiles(net; path = mktempdir(), created = _EXPORT_STAMP)
    ns = Sparlectra.CGMESImporter.CGMES_UUID_NAMESPACE
    @test net.cgmes_ids["TN|A"] == string(UUIDs.uuid5(ns, "TN|A"))
    @test net.cgmes_ids["ACL|A|B|2"] == string(UUIDs.uuid5(ns, "ACL|A|B|2"))
    # the FullModel ids are minted over the key plus the net name: two
    # networks never share an EQ model id (L8, every export carried the same
    # one and the Web UI binds measurement tables by it)
    @test net.cgmes_ids["MODEL|EQ"] == string(UUIDs.uuid5(ns, "MODEL|EQ|exportnet"))
    other = _export_test_net()
    other.name = "othernet"
    dir_a = mktempdir()
    dir_b = mktempdir()
    writeCGMESFiles(_export_test_net(); path = dir_a, created = _EXPORT_STAMP)
    writeCGMESFiles(other; path = dir_b, created = _EXPORT_STAMP)
    eq_a = Sparlectra.CGMESImporter.cgmesEquipmentModelIds(dir_a)
    eq_b = Sparlectra.CGMESImporter.cgmesEquipmentModelIds(dir_b)
    @test length(eq_a) == 1 && length(eq_b) == 1 && eq_a != eq_b
  end)() end

  @testset "independent builds identical, names carry no identity" begin (function ()
    n1 = _export_test_net()
    n2 = _export_test_net()
    f1 = writeCGMESFiles(n1; path = mktempdir(), created = _EXPORT_STAMP)
    f2 = writeCGMESFiles(n2; path = mktempdir(), created = _EXPORT_STAMP)
    # both nets start with an empty cgmes_ids dict -> same minted ids
    @test n1.cgmes_ids == n2.cgmes_ids
    for i in eachindex(f1)
      @test read(f1[i]) == read(f2[i])
    end
    # a second export of the same net, its ids already minted, is byte-equal
    f1b = writeCGMESFiles(n1; path = mktempdir(), created = _EXPORT_STAMP)
    @test all(read(f1[i]) == read(f1b[i]) for i in eachindex(f1))
    # renaming a component changes only its display name, never its mRID
    n2.linesAC[1].comp.cName = "renamed_line"
    f3 = writeCGMESFiles(n2; path = mktempdir(), created = _EXPORT_STAMP)
    @test n2.cgmes_ids == n1.cgmes_ids
    eq3 = read(f3[1], String)
    @test occursin("rdf:ID=\"_$(n1.cgmes_ids["ACL|A|B|1"])\"", eq3)
    @test occursin(">renamed_line<", eq3)
    # a line pairs with its branch by bus pair and order, not by the branch
    # name: renaming the branch leaves TP and SV unchanged, with the open
    # terminal and two flows per line; a line without a branch stops the
    # export before anything is written
    net = _export_test_net()
    net.branchVec[2].from_status = 0
    plain = mktempdir()
    writeCGMESFiles(net; path = plain, created = _EXPORT_STAMP)
    renamed = deepcopy(net)
    renamed.branchVec[1].comp.cName = "Corridor_north"
    again = mktempdir()
    writeCGMESFiles(renamed; path = again, created = _EXPORT_STAMP)
    for profile in ("TP", "SV")
      @test read(joinpath(plain, "exportnet_$(profile).xml"), String) == read(joinpath(again, "exportnet_$(profile).xml"), String)
    end
    @test count("<cim:ACDCTerminal.connected>false", read(joinpath(again, "exportnet_TP.xml"), String)) == 1
    @test count("<cim:SvPowerFlow ", read(joinpath(again, "exportnet_SV.xml"), String)) == 2 * length(renamed.linesAC)
    broken = _export_test_net()
    push!(broken.linesAC, broken.linesAC[1])
    dir = mktempdir()
    err = @test_throws ErrorException writeCGMESFiles(broken; path = dir, created = _EXPORT_STAMP)
    @test occursin("has no branch in the network", sprint(showerror, err.value))
    @test isempty(readdir(dir))
  end)() end

  @testset "self-roundtrip is power-flow-identical" begin (function ()
    original = _roundtrip_net()
    @test _solve!(original)[2] == 0
    exported = _roundtrip_net()
    # 0.30.2 reference priorities: the slack infeed states 1, the PV
    # machine 3; both travel through SSH and come back on the prosumers
    exported.prosumpsVec[1].referencePriority = 1
    exported.prosumpsVec[2].referencePriority = 3
    @test _solve!(exported)[2] == 0
    dir = mktempdir()
    notes = String[]
    files = writeCGMESFiles(exported; path = dir, created = _EXPORT_STAMP, notices = notes)
    @test length(files) == 4
    # every class in the fixture is representable — nothing may be dropped
    @test isempty(notes)
    eq = read(files[1], String)
    @test occursin("<cim:Breaker", eq)
    @test occursin("<cim:StaticVarCompensator", eq)
    @test occursin("<cim:RatioTapChanger", eq)
    @test occursin("<cim:PhaseTapChangerLinear", eq)
    res = importCGMES(path = dir, name = "rt_back")
    net2 = res.net
    @test res.slack_bus == "HV1"
    prio_by_bus(n) = sort!([(only(k for (k, v) in n.busDict if v == Int(ps.comp.cFrom_bus)), ps.referencePriority) for ps in n.prosumpsVec if ps.referencePriority > 0])
    @test prio_by_bus(net2) == prio_by_bus(exported) == [("HV1", 1), ("HV2", 3)]
    @test length(net2.nodeVec) == length(original.nodeVec)
    @test length(net2.branchVec) == length(original.branchVec)
    @test length(net2.prosumpsVec) == length(original.prosumpsVec)
    @test length(net2.shuntVec) == length(original.shuntVec)
    @test length(net2.linkVec) == length(original.linkVec)
    # the reconstructed transformer carries ratio AND shift exactly
    tbr = [br for br in net2.branchVec if occursin("_2WT_", br.comp.cName)]
    @test length(tbr) == 1
    @test isapprox(tbr[1].ratio, 1.04; atol = 1e-9)
    @test isapprox(tbr[1].angle, 2.0; atol = 1e-9)
    # SVC survives as a regulated unit with its rating-derived Q limits
    svc = [p for p in net2.prosumpsVec if p.comp.cTyp == Sparlectra.StaticVarCompensator]
    @test length(svc) == 1
    @test svc[1].minQ == -121.0
    @test svc[1].maxQ == 121.0
    # the asymmetric P band of the network injection comes back as
    # written, in generator sign, not as the hull (-120, 120)
    eni = only(p for p in net2.prosumpsVec if p.comp.cTyp == Sparlectra.ExternalNetworkInjection)
    @test (eni.minP, eni.maxP) == (-30.0, 120.0)
    # the re-imported model reproduces the SV flows at the SV state
    @test compareWithSV(res).flows.max_dp < 1e-6
    # the exported SV state starts the re-import at the solution
    its, erg = _solve!(net2)
    @test erg == 0
    @test its <= 2
    _compare_solved(original, net2; atol = 1e-8)
    # a ratio step that was not transferred: CGMES has no mark for it, the
    # SSH step is the assumed middle position and one warning names the
    # transformer
    tidx = only(i for (i, br) in enumerate(exported.branchVec) if occursin("_2WT_", br.comp.cName))
    exported.branchVec[tidx].ratio_step_assumed = true
    assumed_note = "transformer $(exported.trafos[1].comp.cName): ratio step not transferred, the export writes the assumed middle position"
    notes2 = String[]
    files2 = writeCGMESFiles(exported; path = mktempdir(), created = _EXPORT_STAMP, notices = notes2)
    @test notes2 == [assumed_note]
    @test occursin("<cim:TapChanger.step>2</cim:TapChanger.step>", read(only(f for f in files2 if endswith(f, "_SSH.xml")), String))
    # read back, the written step counts as transferred
    @test !any(m -> occursin("step not transferred", m), importCGMES(path = dirname(files2[1]), name = "rt_assumed").messages)
    @test_logs (:warn, "CGMES export: " * assumed_note) match_mode = :any writeCGMESFiles(exported; path = mktempdir(), created = _EXPORT_STAMP)
  end)() end

  @testset "SSH and SV profiles carry the operating point" begin (function ()
    net = Net(name = "sshnet", baseMVA = 100.0)
    addBus!(net = net, busName = "A", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
    addBus!(net = net, busName = "B", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
    addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
    addProsumer!(net = net, busName = "A", type = "SYNCHRONOUSMACHINE", p = 25.0, q = 5.0, referencePri = "A", vm_pu = 1.02, isRegulated = true)
    addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 25.0, q = 5.0)
    addShuntMatpower!(net = net, busName = "B", Gs = 0.0, Bs = 4.0)
    # a second line, so that no single outage cuts B off (N-1 check below)
    addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.02, x_pu = 0.1, b_pu = 0.0, status = 1)
    # C hangs on B through a closed switch alone, with a unit, a load and a
    # shunt of its own: the re-import marks it isolated (branch rule), the
    # solver merges it with B, and the SV comparison compares it
    addBus!(net = net, busName = "C", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
    addLink!(net = net, fromBus = "B", toBus = "C", status = 1)
    addProsumer!(net = net, busName = "C", type = "SYNCHRONOUSMACHINE", p = 10.0, q = 2.0)
    addProsumer!(net = net, busName = "C", type = "ENERGYCONSUMER", p = 4.0, q = 1.0)
    addShuntMatpower!(net = net, busName = "C", Gs = 0.5, Bs = 2.0)
    @test _solve!(net)[2] == 0
    files = writeCGMESFiles(net; path = mktempdir(), created = _EXPORT_STAMP)
    ssh = read(files[3], String)
    # machine convention: injection-positive storage flips sign in SSH
    @test occursin("<cim:RotatingMachine.p>-25.0</cim:RotatingMachine.p>", ssh)
    @test occursin("<cim:RotatingMachine.q>-5.0</cim:RotatingMachine.q>", ssh)
    @test occursin("<cim:SynchronousMachine.referencePriority>1</cim:SynchronousMachine.referencePriority>", ssh)
    @test occursin("<cim:EnergyConsumer.p>25.0</cim:EnergyConsumer.p>", ssh)
    @test occursin("<cim:EnergyConsumer.q>5.0</cim:EnergyConsumer.q>", ssh)
    # voltage target in kV of the regulated bus
    @test occursin("<cim:RegulatingControl.targetValue>$(1.02 * 110.0)</cim:RegulatingControl.targetValue>", ssh)
    @test occursin("<cim:ShuntCompensator.sections>1</cim:ShuntCompensator.sections>", ssh)
    eq = read(files[1], String)
    @test occursin("RegulatingControlModeKind.voltage", eq)
    @test occursin("<cim:LinearShuntCompensator.bPerSection>$(4.0 / 110.0^2)</cim:LinearShuntCompensator.bPerSection>", eq)
    # SV: the re-import validates cleanly against the exported state
    res = importCGMES(path = dirname(files[1]), name = "sv_back")
    @test isempty(res.no_sv_buses)
    cmp = Sparlectra.CGMESImporter.compareWithSV(res)
    @test maximum(abs(r.dvm) for r in cmp.rows) < 1e-12
    @test maximum(abs(r.dva) for r in cmp.rows) < 1e-12
    @test !isempty(cmp.flows.rows)
    @test maximum(abs(r.dp) for r in cmp.flows.rows) < 1e-9
    @test maximum(abs(r.dq) for r in cmp.flows.rows) < 1e-9
    # C has no branch, so the import marks it isolated, but the closed
    # switch joins it to B: its voltage and its shunt are compared at the
    # voltage of the switch group, and the units of B and C are balanced
    # over the group (C alone would miss the power the switch carries)
    @test res.net.busDict["C"] in res.net.isoNodes
    @test any(r -> r.bus == "C", cmp.rows)
    @test count(r -> r.kind == :shunt, cmp.flows.rows) == 2
    @test count(r -> r.kind == :units, cmp.flows.rows) == 2
    # the machine on C delivers 10 MW and 2 Mvar (load sign: -10, -2); a
    # balance of C alone would give its load and shunt instead
    a_name = res.net.nodeVec[res.net.busDict["A"]].comp.cName
    c_units = only(r for r in cmp.flows.rows if r.kind == :units && r.bus != a_name)
    @test isapprox(c_units.sv_p, -10.0; atol = 1e-9) && isapprox(c_units.sv_q, -2.0; atol = 1e-9)
    @test isapprox(c_units.p, -10.0; atol = 1e-9) && isapprox(c_units.q, -2.0; atol = 1e-9)
    # the solved state matches as well
    @test _solve!(res.net)[2] == 0
    solved = Sparlectra.CGMESImporter.compareWithSV(res)
    @test solved.n == 3
    @test solved.max_dvm < 1e-9
    @test maximum(abs(r.dp) for r in solved.flows.rows) < 1e-6
    # the SV export of the import keeps the voltage and the rows of C
    files2 = writeCGMESFiles(res.net; path = mktempdir(), created = _EXPORT_STAMP)
    res2 = importCGMES(path = dirname(files2[1]), name = "sv_back2")
    cmp2 = Sparlectra.CGMESImporter.compareWithSV(res2)
    @test cmp2.n == 3
    @test cmp2.max_dvm < 1e-9
    @test length(cmp2.flows.rows) == length(cmp.flows.rows)
    @test maximum(abs(r.dp) for r in cmp2.flows.rows) < 1e-6
    # N-1 applies the same branch rule after each outage: C stays marked
    # and is neither cut off nor counted as shed load
    n1 = runContingencies!(res.net, generateN1Branches(res.net); warm_note = Ref(""))
    @test length(n1) == 2
    for r in n1
      @test r.shed_load_mw == 0.0
      @test r.error === nothing
    end
  end)() end

  @testset "a PowerElectronicsConnection goes out as EquivalentInjection" begin (function ()
    # CGMES 2.4.15 has no PowerElectronicsConnection class. One PQ and
    # one voltage-regulating converter, marked as the CGMES importer
    # marks them, leave as EquivalentInjection and come back with the same
    # injection, Q limits, regulation and converter type; without them the
    # same net exports for the short-circuit comparison. `bind` gives the
    # regulating one a Q limit it cannot hold its target with.
    function pec_net(with_pec::Bool; bind::Bool = false)
      net = Net(name = with_pec ? "pecnet" : "nopecnet", baseMVA = 100.0)
      for b in ("A", "B", "C")
        addBus!(net = net, busName = b, vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
      end
      addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.05, b_pu = 0.02, status = 1)
      addPIModelACLine!(net = net, fromBus = "B", toBus = "C", r_pu = 0.015, x_pu = 0.04, b_pu = 0.01, status = 1)
      addProsumer!(net = net, busName = "A", type = "SYNCHRONOUSMACHINE", referencePri = "A", vm_pu = 1.02, va_deg = 0.0)
      addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
      addProsumer!(net = net, busName = "C", type = "ENERGYCONSUMER", p = 20.0, q = bind ? 25.0 : 5.0)
      if with_pec
        addProsumer!(net = net, busName = "B", type = "GENERATOR", p = 12.5, q = 3.0, pMin = -5.0, pMax = 15.0, qMin = -8.0, qMax = 8.0)
        net.prosumpsVec[end].comp.cTyp = Sparlectra.PowerElectronicsConnection
        qlim = bind ? 5.0 : 30.0
        addProsumer!(net = net, busName = "C", type = "GENERATOR", p = 20.0, q = 0.0, vm_pu = bind ? 1.05 : 1.01, isRegulated = true, qMin = -qlim, qMax = qlim)
        net.prosumpsVec[end].comp.cTyp = Sparlectra.PowerElectronicsConnection
      end
      return net
    end
    C = Sparlectra.CGMESImporter
    # the re-imported converters by bus name
    back_of(res) = Dict(only(k for (k, v) in res.net.busDict if v == Int(ps.comp.cFrom_bus)) => ps for ps in (res.net.prosumpsVec[last(res.element_by_mrid[o.mrid])] for o in C.objectsOf(res.store, :EquivalentInjection)))
    net = pec_net(true)
    @test _solve!(net)[2] == 0
    dir = mktempdir()
    notes = String[]
    files = writeCGMESFiles(net; path = dir, created = _EXPORT_STAMP, notices = notes)
    # the class change is the one line of the export report
    @test notes == ["2 PowerElectronicsConnection prosumer(s) written as EquivalentInjection with description PowerElectronicsConnection (CGMES 2.4.15 has no class for them; a re-import keeps them as converters without the fault current maxIFault and ratedS, so Ik''max in their island is flagged as a lower bound)"]
    eq = read(files[1], String)
    ssh = read(files[3], String)
    @test count("<cim:EquivalentInjection rdf:ID", eq) == 2
    @test count("<cim:IdentifiedObject.description>PowerElectronicsConnection</cim:IdentifiedObject.description>", eq) == 2
    @test count("<cim:SynchronousMachine rdf:ID", eq) == 1
    # limits and regulation capability in EQ, regulation in SSH (kV)
    for tag in ("minP>-5.0", "maxP>15.0", "minQ>-8.0", "maxQ>8.0", "minQ>-30.0", "maxQ>30.0")
      @test occursin("<cim:EquivalentInjection.$(tag)</cim:EquivalentInjection.", eq)
    end
    @test count("<cim:EquivalentInjection.regulationCapability>true<", eq) == 1
    @test count("<cim:EquivalentInjection.regulationStatus>true<", ssh) == 1
    @test occursin("<cim:EquivalentInjection.regulationTarget>$(1.01 * 110.0)</cim:EquivalentInjection.regulationTarget>", ssh)
    # the one RegulatingControl is the slack machine's
    @test count("<cim:RegulatingControl rdf:ID", eq) == 1

    res = importCGMES(path = dir, name = "pec_back")
    eis = C.objectsOf(res.store, :EquivalentInjection)
    @test length(eis) == 2 && length(C.objectsOf(res.store, :SynchronousMachine)) == 1
    @test all(first(res.element_by_mrid[o.mrid]) == :external_grid for o in eis)
    back = back_of(res)
    orig = Dict(b => net.prosumpsVec[k] for (b, k) in (("B", 4), ("C", 5)))
    for b in ("B", "C")
      @test (back[b].pVal, back[b].qVal) == (orig[b].pVal, orig[b].qVal)
      @test (back[b].minQ, back[b].maxQ) == (orig[b].minQ, orig[b].maxQ)
      # minP/maxP in generator sign come back as written, (-5, 15) at B
      @test (back[b].minP, back[b].maxP) == (orig[b].minP, orig[b].maxP)
      @test Sparlectra.isRegulating(back[b]) == Sparlectra.isRegulating(orig[b])
      # the description brings the converter type back
      @test back[b].comp.cTyp == Sparlectra.PowerElectronicsConnection
    end
    @test back["C"].vm_pu ≈ 1.01
    @test !Sparlectra.isRegulating(back["B"])
    # the SV of the export validates, and the re-import solves to the
    # same state from it
    cmp = C.compareWithSV(res)
    @test maximum(abs(r.dvm) for r in cmp.rows) < 1e-12
    @test count(r -> r.kind == :units, cmp.flows.rows) == 3
    @test maximum(abs(r.dp) for r in cmp.flows.rows) < 1e-9
    its, erg = _solve!(res.net)
    @test erg == 0 && its <= 2
    _compare_solved(net, res.net)
    # a second export writes them as EquivalentInjection again
    notes2 = String[]
    writeCGMESFiles(res.net; path = mktempdir(), created = _EXPORT_STAMP, notices = notes2)
    @test notes2 == notes

    # a Q limit that binds: C cannot hold 1.05 pu with 5 MVAr and ends PQ
    # at the limit, before the export and after the re-import
    bnet = pec_net(true; bind = true)
    @test _solve!(bnet)[2] == 0
    @test bnet.nodeVec[bnet.busDict["C"]]._nodeType == Sparlectra.PQ
    bdir = mktempdir()
    writeCGMESFiles(bnet; path = bdir, created = _EXPORT_STAMP)
    bres = importCGMES(path = bdir, name = "pec_bind_back")
    @test (back_of(bres)["C"].minQ, back_of(bres)["C"].maxQ) == (-5.0, 5.0)
    @test _solve!(bres.net)[2] == 0
    @test bres.net.nodeVec[bres.net.busDict["C"]]._nodeType == Sparlectra.PQ
    _compare_solved(bnet, bres.net)

    # no machine record for them in the short-circuit harvest: Ik'' as
    # the network without them, and every row flagged as a lower bound
    # for the converter fault current the export could not carry
    @test all(m -> m.bus == "A", res.shortcircuit.synchronous_machines)
    @test isempty(res.shortcircuit.equivalent_injections) && isempty(res.shortcircuit.power_electronics_connections)
    ref_dir = mktempdir()
    ref_net = pec_net(false)
    @test _solve!(ref_net)[2] == 0
    writeCGMESFiles(ref_net; path = ref_dir, created = _EXPORT_STAMP)
    sc_rows(r) = Dict(row.bus => row for row in run_with_expected_warnings(() -> runShortCircuit!(r; case = :max), ["short-circuit: SynchronousMachine", "short-circuit: PowerElectronicsConnection"]).rows)
    sc_pec, sc_ref = sc_rows(res), sc_rows(importCGMES(path = ref_dir, name = "nopec_back"))
    @test sort!(collect(keys(sc_pec))) == sort!(collect(keys(sc_ref))) == ["A", "B", "C"]
    @test all(isapprox(sc_pec[b].ik_kA, sc_ref[b].ik_kA; rtol = 1e-12) for b in keys(sc_ref))
    for (b, row) in sc_pec
      @test row.contains_defaulted_data
      @test count(r -> occursin("PowerElectronicsConnection at bus", r) && occursin("Ik''max in its island is a lower bound", r), row.reasons) == 2
      @test !any(r -> occursin("PowerElectronicsConnection", r), sc_ref[b].reasons)
    end
  end)() end

  @testset "three-winding transformer with a tapped leg round-trips" begin (function ()
    # the importer builds the tapped end-1 leg side bus -> star; the
    # export must still write one PowerTransformer with three ends
    res = importCGMES(path = _cgmes_synth_tapside_dir(:three_winding), name = "tap3w")
    @test _solve!(res.net)[2] == 0
    notes = String[]
    files = writeCGMESFiles(res.net; path = mktempdir(), created = _EXPORT_STAMP, notices = notes)
    @test isempty(notes)
    eq = read(files[1], String)
    @test count("<cim:PowerTransformer rdf:ID", eq) == 1
    @test count("<cim:PowerTransformerEnd rdf:ID", eq) == 3
    back = importCGMES(path = dirname(files[1]), name = "tap3w_back")
    @test _solve!(back.net)[2] == 0
    for bus in ("BUS1", "BUS2", "BUS3")
      n1 = res.net.nodeVec[res.net.busDict[bus]]
      n2 = back.net.nodeVec[back.net.busDict[bus]]
      @test isapprox(n1._vm_pu, n2._vm_pu; atol = 1e-9)
      @test isapprox(n1._va_deg, n2._va_deg; atol = 1e-7)
    end
    # the exported SV flows of the three ends are the solved leg flows
    @test compareWithSV(back).flows.max_dp < 1e-6
    # ends 1 and 3 carry a magnetizing admittance: before and after the
    # round trip it sits at the star node (the arm seen at the star), the
    # arm at the side bus is zero
    function star_arms(net)
      s = net.busDict["AUX3WT_T3W"]
      arms = Tuple{ComplexF64,ComplexF64}[]
      for br in net.branchVec
        (br.fromBus == s || br.toBus == s) || continue
        y0f = Sparlectra._branch_y0_from(br) / abs2(Sparlectra.calcBranchRatio(br))
        y0t = Sparlectra._branch_y0_to(br)
        push!(arms, br.fromBus == s ? (y0f, y0t) : (y0t, y0f))
      end
      return arms
    end
    for net in (res.net, back.net)
      arms = star_arms(net)
      @test length(arms) == 3
      @test all(a -> iszero(last(a)), arms)
      @test count(a -> !iszero(first(a)), arms) == 2
    end
  end)() end

  @testset "zip packaging re-imports directly" begin (function ()
    net = _export_test_net()
    addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
    files = writeCGMESFiles(net; path = mktempdir(), created = _EXPORT_STAMP, zip = true)
    @test length(files) == 5
    @test endswith(files[5], "_CGMES.zip")
    res = importCGMES(path = files[5], name = "zip_back")
    @test length(res.net.nodeVec) == 3
    @test length(res.net.linesAC) == 3
  end)() end

  # A branch built from an EquivalentBranch goes out under its class and
  # mRID; written as ACLineSegment or PowerTransformer it would come back
  # as another class with another branch model (the boundary-line ratio
  # across two levels). The branch model of an asymmetric one uses the
  # geometric mean of |z12| and |z21|; the export writes the delivered r, x,
  # r21, x21 from net.equivalent_branch_metadata, so the re-import has the
  # delivery's values and builds the same branch bit for bit. The record
  # survives the typed case (SCF) the service path builds the net from.
  @testset "an imported EquivalentBranch exports as EquivalentBranch" begin (function ()
    num = Sparlectra.CGMESImporter.num
    fields = (:r_pu, :x_pu, :g_from_pu, :b_from_pu, :g_to_pu, :b_to_pu, :ratio, :tap_ratio, :angle, :phase_shift_deg)
    for mag in (1.0, 1.01)
      res = importCGMES(path = _cgmes_synth_eb_dir(mag = mag), name = "synth_eb")
      @test _solve!(res.net)[2] == 0
      @test length(res.net.equivalent_branch_metadata) == 3
      # the typed case keeps the record per branch
      case = Sparlectra.cgmes_enrich_case!(Sparlectra.net_to_scfcase(res.net; source_format = "cgmes", source_reference = "synth_eb", include_start_state = true), res.net)
      @test Sparlectra.build_net(case).equivalent_branch_metadata == res.net.equivalent_branch_metadata
      for (label, net) in (("import", res.net), ("typed case", Sparlectra.build_net(case)))
        @test _solve!(net)[2] == 0
        dir = mktempdir()
        notes = String[]
        files = writeCGMESFiles(net; path = dir, created = _EXPORT_STAMP, notices = notes)
        # the export names the objects whose model and export values differ
        @test notes == (mag == 1.0 ? String[] : ["2 EquivalentBranch object(s) with |z21| other than |z12|: r, x, r21, x21 exported as delivered, the branch model uses the geometric mean of |z12| and |z21|"])
        eq = read(files[1], String)
        @test count("<cim:EquivalentBranch ", eq) == 3
        @test all(m -> occursin("<cim:EquivalentBranch rdf:ID=\"_$(m)\">", eq), ("eb_same", "eb_asym", "eb_cross"))
        @test count("<cim:EquivalentBranch.r21>", eq) == 2 && count("<cim:EquivalentBranch.x21>", eq) == 2
        @test count("<cim:EquivalentNetwork ", eq) == 1
        back = importCGMES(path = dir, name = "eb_back")
        @test count(m -> startswith(m, "notice: EquivalentBranch: 3 mapped (2 same voltage level, 1 across voltage levels)") && occursin("2 asymmetric", m), back.messages) == 1
        name_of(n) = Dict(v => k for (k, v) in n.busDict)
        z21s = _cgmes_synth_eb_z21(mag)
        for (m, z12, z21) in (("_eb_same", complex(1.21, 6.05), complex(1.21, 6.05)), ("_eb_asym", complex(1.5, 7.0), z21s[1]), ("_eb_cross", complex(0.6, 4.84), z21s[2]))
          o = back.store.objects[m]
          @test o.class == :EquivalentBranch
          # the delivered values, exactly
          @test complex(num(o, :r), num(o, :x)) == z12
          if m == "_eb_same"
            @test num(o, :r21) === nothing && num(o, :x21) === nothing
          else
            @test complex(num(o, :r21), num(o, :x21)) == z21
          end
          b1 = res.net.branchVec[res.element_by_mrid[m][2]]
          b2 = back.net.branchVec[back.element_by_mrid[m][2]]
          @test (name_of(back.net)[b2.fromBus], name_of(back.net)[b2.toBus]) == (name_of(res.net)[b1.fromBus], name_of(res.net)[b1.toBus])
          # the same branch, bit for bit
          @test all(f -> getfield(b2, f) === getfield(b1, f), fields)
        end
        @test back.net.equivalent_branch_metadata == res.net.equivalent_branch_metadata
        # the exported SV is the solution: the re-import starts there
        its, erg = _solve!(back.net)
        @test erg == 0 && its <= 2
        _compare_solved(res.net, back.net)
      end
    end
  end)() end

  # Export-import-export on the checked-in deliveries (data/cgmes_demo,
  # written by tools/gen_cgmes_fixtures.jl with the same header stamp):
  # the re-export of an imported delivery reproduces every object of the
  # fixture with its mRID and every attribute value, except what the
  # importer does not carry into the net. Named exactly:
  #   1. RatioTapChanger and TapChangerControl objects (EQ and SSH): the
  #      importer folds the ratio-tap step into the branch ratio and keeps
  #      the range only as branch nameplate data, the exporter writes tap
  #      machinery from the winding record alone, so the re-export has
  #      none; the effective ratio itself is reproduced (the SV profile
  #      and the transformer parameters compare equal).
  #   2. PowerTransformerEnd.ratedU of an end that carried a ratio tap
  #      changer: the fixture absorbed the step correction into ratedU so
  #      the live ratio survives the re-import, the re-export writes the
  #      live ratio without a step. The g and b of such an end are in
  #      siemens on that ratedU and are compared in per unit
  #      (g * ratedU^2, b * ratedU^2), not dropped.
  #   3. SynchronousMachine.minQ/maxQ: the importer reads them as the hull
  #      of both sign readings (sp_case118 has asymmetric pairs) and
  #      substitutes wide symmetric limits where the fixture has none
  #      (sp_case14 Gen_110).
  #   4. IdentifiedObject.name of a LinearShuntCompensator and of its
  #      Terminal: the importer names shunts by its own bus index.
  # Everything else (topology, lines, transformers, loads, machines,
  # regulating controls, breakers, the SSH operating point, the SV
  # voltages and flows) must compare equal object by object.
  @testset "export-import-export identity on the checked-in deliveries" begin (function ()
    # each delivery carries its own EQ model id (L8: all three shared one,
    # so a measurement table bound to one demo case showed up at all three)
    demo_eq = [Sparlectra.CGMESImporter.cgmesEquipmentModelIds(cgmes_fixture_dir(c)) for c in ("sp_case14", "sp_case118", "sp_casePST")]
    @test all(ids -> length(ids) == 1, demo_eq)
    @test allunique(first.(demo_eq))
    for case in ("sp_case14", "sp_case118", "sp_casePST")
      dir = cgmes_fixture_dir(case)
      res = importCGMES(path = dir, name = case)
      # the SV of the re-export is the solved state: solve from the SV
      # start with Q-limits off, as the fixture was solved (a unit that
      # switched to a limit would legitimately change SSH q and SV)
      @test runpf!(res.net, 60, 1e-8, 0; method = :rectangular, qlimits_enabled = false)[2] == 0
      out = mktempdir()
      notes = String[]
      files = writeCGMESFiles(res.net; path = out, created = _EXPORT_STAMP, notices = notes)
      @test isempty(notes)
      @test [basename(f) for f in files] == [string(case, "_", p, ".xml") for p in ("EQ", "TP", "SSH", "SV")]
      eq_fixture = read(joinpath(dir, basename(files[1])), String)
      _, eq_objects = _cgmes_profile_objects(eq_fixture)
      ref_id(v) = String(last(split(v, "_"; limit = 2)))
      rtc_ends = Set(ref_id(v) for (_, (cls, attrs)) in eq_objects if cls == "cim:RatioTapChanger" for (tag, v) in attrs if tag == "cim:RatioTapChanger.TransformerEnd")
      shunt_ids = Set(id for (id, (cls, _)) in eq_objects if cls == "cim:LinearShuntCompensator")
      shunt_terminals = Set(id for (id, (cls, attrs)) in eq_objects if cls == "cim:Terminal" && any(tag == "cim:Terminal.ConductingEquipment" && ref_id(v) in shunt_ids for (tag, v) in attrs))
      drop_attr = (cls, id, tag) -> begin
        (cls == "cim:PowerTransformerEnd" && tag == "cim:PowerTransformerEnd.ratedU" && id in rtc_ends) ||
          (cls == "cim:SynchronousMachine" && tag in ("cim:SynchronousMachine.minQ", "cim:SynchronousMachine.maxQ")) ||
          (cls == "cim:LinearShuntCompensator" && tag == "cim:IdentifiedObject.name") ||
          (cls == "cim:Terminal" && tag == "cim:IdentifiedObject.name" && id in shunt_terminals)
      end
      counts = Int[]
      for f in files
        # per_unit_ends: the ends under exception 2 carry magnetizing since
        # the demo cases were rebuilt (0.50.0, transformer magnetizing from
        # nameplate since 0.32.1); their g and b refer to the differing
        # ratedU and are compared in per unit
        compared = _compare_cgmes_profile(read(joinpath(dir, basename(f)), String), read(f, String); drop_classes = ("cim:RatioTapChanger", "cim:TapChangerControl"), drop_attr = drop_attr, per_unit_ends = rtc_ends)
        push!(counts, compared)
      end
      println("      ", case, ": objects compared EQ/TP/SSH/SV = ", join(counts, "/"), ", ratio-tap ends excluded: ", length(rtc_ends))
      @test all(>(0), counts)
      @test !isempty(rtc_ends)
      # the excluded tap machinery is exactly what the fixture carries
      @test count(o -> o[2][1] == "cim:RatioTapChanger", collect(eq_objects)) == length(rtc_ends)
    end
  end)() end

  @testset "lines in ohm and siemens export unchanged with their length" begin (function ()
    # addACLine! and addACLinePerKm! store ohm and siemens totals; the
    # exporter writes them as they are (no per-unit round trip) and takes
    # Conductor.length from length_km
    net = Net(name = "ohmnet", baseMVA = 100.0)
    addBus!(net = net, busName = "A", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
    addBus!(net = net, busName = "B", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
    addBus!(net = net, busName = "C", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
    addACLine!(net = net, fromBus = "A", toBus = "B", r_ohm = 1.3, x_ohm = 4.7, b_s = 3.1e-5, g_s = 1.7e-7, length_km = 12.5)
    addACLinePerKm!(net = net, fromBus = "B", toBus = "C", length_km = 8.0, r_ohm_per_km = 0.121, x_ohm_per_km = 0.391, c_nf_per_km = 9.5, g_us_per_km = 0.02)
    addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
    addProsumer!(net = net, busName = "C", type = "ENERGYCONSUMER", p = 20.0, q = 5.0)
    @test _solve!(net)[2] == 0
    f_hz = Sparlectra.model_config().nominal_frequency_hz
    expected = [(1.3, 4.7, 3.1e-5, 1.7e-7, 12.5), (0.121 * 8.0, 0.391 * 8.0, 2.0 * pi * f_hz * 9.5 * 1.0e-9 * 8.0, 0.02 * 1.0e-6 * 8.0, 8.0)]
    dir = mktempdir()
    files = writeCGMESFiles(net; path = dir, created = _EXPORT_STAMP)
    eq = read(files[1], String)
    segs = [m.match for m in eachmatch(r"<cim:ACLineSegment rdf:ID=.*?</cim:ACLineSegment>"s, eq)]
    @test length(segs) == 2
    val(seg, tag) = parse(Float64, match(Regex("<cim:$(tag)>([^<]*)</cim:$(tag)>"), seg).captures[1])
    for (seg, (r, x, b, g, len)) in zip(segs, expected)
      @test val(seg, "ACLineSegment.r") == r
      @test val(seg, "ACLineSegment.x") == x
      @test val(seg, "ACLineSegment.bch") == b
      @test val(seg, "ACLineSegment.gch") == g
      @test val(seg, "Conductor.length") == len
    end
    # the re-import builds the same admittance matrix
    net2 = importCGMES(path = dir, name = "ohm_back").net
    buses = ["A", "B", "C"]
    Y1 = Matrix(createYBUS(net = net))[[net.busDict[b] for b in buses], [net.busDict[b] for b in buses]]
    Y2 = Matrix(createYBUS(net = net2))[[net2.busDict[b] for b in buses], [net2.busDict[b] for b in buses]]
    @test isapprox(Y1, Y2; rtol = 1e-12)
  end)() end

  @testset "duplicate mRID aborts before writing" begin (function ()
    net = _export_test_net()
    net.cgmes_ids["TN|A"] = "deadbeef-0000-0000-0000-000000000001"
    net.cgmes_ids["TN|B"] = "deadbeef-0000-0000-0000-000000000001"
    dir = mktempdir()
    err = try
      writeCGMESFiles(net; path = dir, created = _EXPORT_STAMP)
      nothing
    catch e
      e
    end
    @test err isa ErrorException
    msg = sprint(showerror, err)
    @test occursin("TN|A", msg)
    @test occursin("TN|B", msg)
    # the guard fires before any file is opened
    @test isempty(readdir(dir))
  end)() end
end

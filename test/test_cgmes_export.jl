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
#          power-flow-identical self-roundtrip on EQ+TP+SSH+SV (transformers
#          incl. phase shift and ratio-tap machinery, machines, SVC, loads,
#          shunts, bus links), SSH/SV content, ZIP packaging and tool
#          provenance on nets built in memory, and the export-import-export
#          identity of the checked-in deliveries under data/cgmes_demo
#          (object-wise, with the fields the importer does not carry named).

using Test
using Sparlectra
using Dates
using UUIDs
using Logging

# Fixed header stamp: with `created` pinned the exported files are
# byte-reproducible, which the determinism tests rely on.
const _EXPORT_STAMP = DateTime(2026, 1, 1, 12, 0, 0)

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
  addLink!(net = net, fromBus = "HV1", toBus = "HV1B", status = 1)
  addPIModelTrafo!(net = net, fromBus = "HV2", toBus = "MV1", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, status = 1, ratio = 1.04, shift_deg = 2.0)
  addProsumer!(net = net, busName = "HV1", type = "EXTERNALNETWORKINJECTION", referencePri = "HV1", vm_pu = 1.02, va_deg = 0.0)
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
      push!(header, String(line))
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
function _compare_cgmes_profile(fixture_text::AbstractString, reexport_text::AbstractString; drop_classes = (), drop_attr = (cls, id, tag) -> false)::Int
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
  @testset "CGMES export identity" begin
    @testset "structural keys and parallel lines" begin
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
    end

    @testset "minted ids are uuid5 over the key" begin
      net = _export_test_net()
      writeCGMESFiles(net; path = mktempdir(), created = _EXPORT_STAMP)
      ns = Sparlectra.CGMESImporter.CGMES_UUID_NAMESPACE
      @test net.cgmes_ids["TN|A"] == string(UUIDs.uuid5(ns, "TN|A"))
      @test net.cgmes_ids["ACL|A|B|2"] == string(UUIDs.uuid5(ns, "ACL|A|B|2"))
    end

    @testset "re-export is byte-identical" begin
      net = _export_test_net()
      f1 = writeCGMESFiles(net; path = mktempdir(), created = _EXPORT_STAMP)
      f2 = writeCGMESFiles(net; path = mktempdir(), created = _EXPORT_STAMP)
      for i in eachindex(f1)
        @test read(f1[i]) == read(f2[i])
      end
    end

    @testset "independent builds identical, names carry no identity" begin
      n1 = _export_test_net()
      n2 = _export_test_net()
      f1 = writeCGMESFiles(n1; path = mktempdir(), created = _EXPORT_STAMP)
      f2 = writeCGMESFiles(n2; path = mktempdir(), created = _EXPORT_STAMP)
      # both nets start with an empty cgmes_ids dict -> same minted ids
      @test n1.cgmes_ids == n2.cgmes_ids
      for i in eachindex(f1)
        @test read(f1[i]) == read(f2[i])
      end
      # renaming a component changes only its display name, never its mRID
      n2.linesAC[1].comp.cName = "renamed_line"
      f3 = writeCGMESFiles(n2; path = mktempdir(), created = _EXPORT_STAMP)
      @test n2.cgmes_ids == n1.cgmes_ids
      eq3 = read(f3[1], String)
      @test occursin("rdf:ID=\"_$(n1.cgmes_ids["ACL|A|B|1"])\"", eq3)
      @test occursin(">renamed_line<", eq3)
    end

    @testset "self-roundtrip is power-flow-identical" begin
      original = _roundtrip_net()
      @test _solve!(original)[2] == 0
      exported = _roundtrip_net()
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
      # the exported SV state starts the re-import at the solution
      its, erg = _solve!(net2)
      @test erg == 0
      @test its <= 2
      _compare_solved(original, net2)
    end

    @testset "SSH and SV profiles carry the operating point" begin
      net = Net(name = "sshnet", baseMVA = 100.0)
      addBus!(net = net, busName = "A", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
      addBus!(net = net, busName = "B", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
      addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.05, b_pu = 0.0, status = 1)
      addProsumer!(net = net, busName = "A", type = "SYNCHRONOUSMACHINE", p = 25.0, q = 5.0, referencePri = "A", vm_pu = 1.02, isRegulated = true)
      addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 25.0, q = 5.0)
      addShuntMatpower!(net = net, busName = "B", Gs = 0.0, Bs = 4.0)
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
    end

    @testset "regulated tap group exports one shared TapChangerControl" begin
      # #322 export half: master and follower reference the SAME control
      # and both carry controlEnabled, so a reimport regroups them instead
      # of seeing independent (fighting) controllers
      n = Net(name = "tccpar", baseMVA = 100.0)
      for (b, vn) in (("H1", 110.0), ("H2", 110.0), ("L1", 20.0), ("L2", 20.0))
        addBus!(net = n, busName = b, vn_kV = vn)
      end
      addProsumer!(net = n, busName = "H1", type = "EXTERNALNETWORKINJECTION", referencePri = "H1", vm_pu = 1.02, va_deg = 0.0)
      addProsumer!(net = n, busName = "L1", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
      addProsumer!(net = n, busName = "L2", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
      addPIModelACLine!(net = n, fromBus = "H1", toBus = "H2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = n, fromBus = "L1", toBus = "L2", r_pu = 0.02, x_pu = 0.1, b_pu = 0.0, status = 1)
      add2WTrafo!(net = n, fromBus = "H1", toBus = "L1", sn_mva = 63.0, vk_percent = 12.0, vkr_percent = 0.4, pfe_kw = 0.0, i0_percent = 0.0)
      add2WTrafo!(net = n, fromBus = "H1", toBus = "L1", sn_mva = 63.0, vk_percent = 12.0, vkr_percent = 0.4, pfe_kw = 0.0, i0_percent = 0.0)
      for tf in n.trafos
        tf.side1.taps = PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = -9, highStep = 9, neutralStep = 0, voltageIncrement_kV = 110.0 * 0.00625)
      end
      for br in n.branchVec
        br.ratio == 0.0 && continue
        br.has_ratio_tap = true
        br.tap_step = 0.00625
      end
      m = string(n.branchVec[3].branchIdx)
      fl = string(n.branchVec[4].branchIdx)
      addPowerTransformerControl!(n; trafo = m, followers = [fl], mode = :voltage, target_bus = "L1", target_vm_pu = 1.01, deadband_vm_pu = 0.004)
      dir = mktempdir()
      files = writeCGMESFiles(n; path = dir, created = _EXPORT_STAMP)
      eq = read(files[1], String)
      ssh = read(files[3], String)
      @test count("cim:TapChangerControl rdf:ID", eq) == 1
      @test count("TapChanger.TapChangerControl rdf:resource", eq) == 2
      @test count("TapChanger.controlEnabled>true", ssh) == 2
      @test occursin("RegulatingControl.targetDeadband", ssh)
      res = importCGMES(path = dir, name = "tcc_back", tap_control = true)
      ctrls = Sparlectra._tap_controllers(res.net)
      @test length(ctrls) == 1
      c = only(ctrls)
      @test length(c.followers) == 1
      @test isapprox(something(c.target_vm_pu, NaN), 1.01; atol = 1e-9)
      @test isapprox(c.deadband_vm_pu, 0.004; atol = 1e-9)
    end

    @testset "zip packaging re-imports directly" begin
      net = _export_test_net()
      addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      files = writeCGMESFiles(net; path = mktempdir(), created = _EXPORT_STAMP, zip = true)
      @test length(files) == 5
      @test endswith(files[5], "_CGMES.zip")
      res = importCGMES(path = files[5], name = "zip_back")
      @test length(res.net.nodeVec) == 3
      @test length(res.net.linesAC) == 3
    end

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
    #      live ratio without a step.
    #   3. SynchronousMachine.minQ/maxQ: the importer reads them as the hull
    #      of both sign readings (sp_case118 has asymmetric pairs) and
    #      substitutes wide symmetric limits where the fixture has none
    #      (sp_case14 Gen_110).
    #   4. IdentifiedObject.name of a LinearShuntCompensator and of its
    #      Terminal: the importer names shunts by its own bus index.
    # Everything else (topology, lines, transformers, loads, machines,
    # regulating controls, breakers, the SSH operating point, the SV
    # voltages and flows) must compare equal object by object.
    @testset "export-import-export identity on the checked-in deliveries" begin
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
          compared = _compare_cgmes_profile(read(joinpath(dir, basename(f)), String), read(f, String); drop_classes = ("cim:RatioTapChanger", "cim:TapChangerControl"), drop_attr = drop_attr)
          push!(counts, compared)
        end
        println("      ", case, ": objects compared EQ/TP/SSH/SV = ", join(counts, "/"), ", ratio-tap ends excluded: ", length(rtc_ends))
        @test all(>(0), counts)
        @test !isempty(rtc_ends)
        # the excluded tap machinery is exactly what the fixture carries
        @test count(o -> o[2][1] == "cim:RatioTapChanger", collect(eq_objects)) == length(rtc_ends)
      end
    end

    @testset "duplicate mRID aborts before writing" begin
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
    end
  end
end

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

# file: tools/gen_cgmes_fixtures.jl
# purpose: build the checked-in CGMES test deliveries under
#          data/cgmes_demo/<case>/ from the shipped cases: import the
#          case, solve it with the rectangular solver, export EQ+TP+SSH+SV
#          with writeCGMESFiles under a FIXED header stamp. The output is
#          byte-reproducible (deterministic uuid5 ids, pinned `created`), so
#          rerunning the script on an unchanged case leaves git clean. No
#          zip is written; the tests pack one at run time when they need it.
#
# usage: julia --project=. tools/gen_cgmes_fixtures.jl [sp_case14 ...]

using Sparlectra
using Dates
using Printf

# One stamp for every fixture: `created` is the only time-dependent value
# in the exported profiles, so pinning it is what makes the regeneration
# byte-identical. The tests that re-export an imported fixture use the same
# value (test/test_cgmes_export.jl, _EXPORT_STAMP).
const CGMES_FIXTURE_STAMP = DateTime(2026, 1, 1, 12, 0, 0)

const _DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))
const _FIXTURE_ROOT = normpath(joinpath(@__DIR__, "..", "data", "cgmes_demo"))

# The MATPOWER import options of the workshops and of load_fixture_net in
# the test suite (test/test_runner_helpers.jl); the fixture must be the same
# network the tests compare against.
function _load_matpower(path::AbstractString)
  return createNetFromMatPowerFile(
    filename = path,
    flatstart = false,
    enable_pq_gen_controllers = true,
    bus_shunt_model = :admittance,
    matpower_shift_sign = 1.0,
    matpower_shift_unit = :deg,
    matpower_ratio = :normal,
    tap_changer_model = :ideal,
  )
end

# case name => loader. sp_case14 carries the OLTC tap controller, sp_case118
# the Q-limited units and synchronous condensers, sp_casePST the
# phase-shifting transformer and a bus link.
const CGMES_FIXTURE_CASES = [
  "sp_case14" => () -> importSCF(joinpath(_DATA_DIR, "scf", "sp_case14.scf.json")),
  "sp_case118" => () -> _load_matpower(joinpath(_DATA_DIR, "mpower", "sp_case118.m")),
  "sp_casePST" => () -> importSCF(joinpath(_DATA_DIR, "scf", "sp_casePST.scf.json")),
]

# The exporter writes ratio-tap machinery (RatioTapChanger range, neutral,
# current step) only from `PowerTransformerTaps` on a transformer winding,
# and a voltage tap controller only when every group member carries that
# machinery. The SCF and MATPOWER importers fill the branch nameplate
# instead (has_ratio_tap, tap_step, tap_min, tap_max, tap_ratio on the
# cascade grid tap = neutral / (1 + n * step)); this maps the nameplate back
# onto the winding so the OLTC controller of sp_case14 reaches the delivery
# as a TapChangerControl. Transformer branches and net.trafos pair
# positionally, the same pairing the exporter relies on.
function _winding_taps_from_nameplate!(net::Net)
  # by build kind, not by the B_2WT_ name marker: SCF cases restore their
  # own branch names, so the marker is missing on sp_case14's transformers
  trafo_branches = [br for br in net.branchVec if Sparlectra._branch_build_kind(br) === :transformer]
  length(trafo_branches) == length(net.trafos) || return 0
  applied = 0
  for (tf, br) in zip(net.trafos, trafo_branches)
    br.has_ratio_tap && br.tap_step > 0.0 || continue
    tf.side1.taps === nothing && tf.side2.taps === nothing || continue
    base = br.ratio
    step_of(ratio) = round(Int, (base / ratio - 1.0) / br.tap_step)
    vn1 = Sparlectra.getNodeVn(net.nodeVec[br.fromBus])
    tf.side1.taps = Sparlectra.PowerTransformerTaps(
      Vn_kV = vn1,
      step = step_of(br.tap_ratio),
      lowStep = step_of(br.tap_max),
      highStep = step_of(br.tap_min),
      neutralStep = 0,
      voltageIncrement_kV = vn1 * br.tap_step,
    )
    applied += 1
  end
  return applied
end

function build_cgmes_fixture(case::AbstractString, loader)
  net = loader()
  windings = _winding_taps_from_nameplate!(net)
  # the exporter names the profile files after the net, so the fixture
  # folder and the file prefix agree regardless of what the case file calls
  # itself
  net.name = String(case)
  # Q-limits stay off in the fixture solve. A delivery carries voltage
  # targets and the SSH operating point, not the solver's active set: a
  # unit that ended Q-pinned in the source solve re-imports as a regulated
  # PV unit at its target, the importer reads minQ/maxQ as the hull of both
  # sign readings, and the switching path then decides which units pin
  # again (measured on sp_case118: 14 pinned units, 0.012 pu against the SV
  # profile even with the source limits restored). With every regulated
  # unit holding its target, the SV profile is the solution of the
  # delivery itself, and the tests reproduce it to floating-point noise.
  ite, erg = runpf!(net, 60, 1e-8, 0; method = :rectangular, qlimits_enabled = false)
  erg == 0 || error(string(case, ": power flow did not converge (erg = ", erg, ")"))
  outdir = joinpath(_FIXTURE_ROOT, case)
  mkpath(outdir)
  # start from an empty folder so a renamed or dropped profile never lingers
  for f in readdir(outdir; join = true)
    rm(f)
  end
  notices = String[]
  files = writeCGMESFiles(net; path = outdir, created = CGMES_FIXTURE_STAMP, notices = notices)
  total = 0
  println(case, ": ", length(net.nodeVec), " buses, ", length(net.branchVec), " branches, solved in ", ite, " iterations, ", windings, " winding(s) with ratio-tap machinery from the branch nameplate")
  for f in files
    sz = filesize(f)
    total += sz
    @printf("  %-22s %9d bytes\n", basename(f), sz)
  end
  @printf("  %-22s %9d bytes\n", "total", total)
  for n in notices
    println("  exporter notice: ", n)
  end
  return total
end

# --- sp_case6_igm: two operator models and a boundary set -------------------
#
# A hand-written CGM, not an export: the exporter writes one model, and the
# area of a bus in a CGMES import comes from the operator model that defines
# it (cgmes_import.area_source). Six 110 kV buses:
#
#   Nord (http://nord.example/sparlectra): Holm (ExternalNetworkInjection,
#     slack), Rabe (load); sub-regions Nord-West (Holm) and Nord-Ost (Rabe)
#     of the region Nord; control area CA_Nord
#   Sued (http://sued.example/sparlectra): Marbach and Nesselau (loads);
#     sub-region Sued of the region Sued; control area CA_Sued
#   boundary set: the X nodes XN_Rabe_Mar and XN_Holm_Nes
#
# Each tie line is two halves, one in each operator model, joined at its X
# node (L_Rabe_X1 + L_X1_Mar, L_Holm_X2 + L_X2_Nes); the TieFlow objects sit
# on the X-side terminals of the halves. Lines 0.12 + j0.4 ohm/km, 9.5 nF/km.
# The files carry no version stamp, so a release does not change them.

const _IGM_NS = """xmlns:cim="http://iec.ch/TC57/2013/CIM-schema-cim16#" xmlns:md="http://iec.ch/TC57/61970-552/ModelDescription/1#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#\""""
const _IGM_PROFILE = Dict(
  :EQ => "http://entsoe.eu/CIM/EquipmentCore/3/1",
  :TP => "http://entsoe.eu/CIM/Topology/4/1",
  :SSH => "http://entsoe.eu/CIM/SteadyStateHypothesis/1/1",
  :EQ_BD => "http://entsoe.eu/CIM/EquipmentBoundary/3/1",
  :TP_BD => "http://entsoe.eu/CIM/TopologyBoundary/3/1",
)
const _IGM_AUTHORITY = Dict("Nord" => "http://nord.example/sparlectra", "Sued" => "http://sued.example/sparlectra", "BD" => "http://boundary.example/sparlectra")

# name, operator, sub-region
const _IGM_BUSES = [("Holm", "Nord", "Nord-West"), ("Rabe", "Nord", "Nord-Ost"), ("Marbach", "Sued", "Sued"), ("Nesselau", "Sued", "Sued")]
const _IGM_XNODES = ["XN_Rabe_Mar", "XN_Holm_Nes"]
# name, operator, from, to, km; an X node as end makes the line a tie half
const _IGM_LINES = [
  ("L_Holm_Rabe", "Nord", "Holm", "Rabe", 20.0),
  ("L_Rabe_X1", "Nord", "Rabe", "XN_Rabe_Mar", 15.0),
  ("L_Holm_X2", "Nord", "Holm", "XN_Holm_Nes", 25.0),
  ("L_X1_Mar", "Sued", "XN_Rabe_Mar", "Marbach", 10.0),
  ("L_X2_Nes", "Sued", "XN_Holm_Nes", "Nesselau", 20.0),
  ("L_Mar_Nes", "Sued", "Marbach", "Nesselau", 12.0),
]
# bus => (P MW, Q Mvar), load sign
const _IGM_LOADS = Dict("Rabe" => (30.0, 10.0), "Marbach" => (25.0, 8.0), "Nesselau" => (15.0, 5.0))

_igm_id(parts...) = string("_", join(parts, "_"))
_igm_fmt(x::Float64) = @sprintf("%.6g", x)

function _igm_head(io, profile::Symbol, op::AbstractString)
  println(io, """<?xml version="1.0" encoding="UTF-8"?>""")
  println(io, "<rdf:RDF ", _IGM_NS, ">")
  println(io, "  <md:FullModel rdf:about=\"urn:uuid:sp-case6-igm-", lowercase(op), "-", lowercase(String(profile)), "\">")
  println(io, "    <md:Model.scenarioTime>2026-01-01T12:00:00Z</md:Model.scenarioTime>")
  println(io, "    <md:Model.description>sp_case6_igm, operator model ", op, "</md:Model.description>")
  println(io, "    <md:Model.profile>", _IGM_PROFILE[profile], "</md:Model.profile>")
  println(io, "    <md:Model.modelingAuthoritySet>", _IGM_AUTHORITY[op], "</md:Model.modelingAuthoritySet>")
  println(io, "  </md:FullModel>")
end

_igm_name(io, name) = println(io, "    <cim:IdentifiedObject.name>", name, "</cim:IdentifiedObject.name>")

function _igm_terminal(io, eq::AbstractString, seq::Int)
  println(io, "  <cim:Terminal rdf:ID=\"", _igm_id("t", eq, seq), "\">")
  println(io, "    <cim:Terminal.ConductingEquipment rdf:resource=\"#", _igm_id(eq), "\"/>")
  println(io, "    <cim:ACDCTerminal.sequenceNumber>", seq, "</cim:ACDCTerminal.sequenceNumber>")
  println(io, "  </cim:Terminal>")
end

_igm_tn(bus) = _igm_id("tn", bus)

function _igm_eq(op::AbstractString)::String
  io = IOBuffer()
  _igm_head(io, :EQ, op)
  regions = op == "Nord" ? ["Nord-West", "Nord-Ost"] : ["Sued"]
  println(io, "  <cim:GeographicalRegion rdf:ID=\"", _igm_id("gr", op), "\">")
  _igm_name(io, op)
  println(io, "  </cim:GeographicalRegion>")
  for r in regions
    println(io, "  <cim:SubGeographicalRegion rdf:ID=\"", _igm_id("sgr", r), "\">")
    _igm_name(io, r)
    println(io, "    <cim:SubGeographicalRegion.Region rdf:resource=\"#", _igm_id("gr", op), "\"/>")
    println(io, "  </cim:SubGeographicalRegion>")
  end
  for (bus, o, r) in _IGM_BUSES
    o == op || continue
    println(io, "  <cim:Substation rdf:ID=\"", _igm_id("ss", bus), "\">")
    _igm_name(io, bus)
    println(io, "    <cim:Substation.Region rdf:resource=\"#", _igm_id("sgr", r), "\"/>")
    println(io, "  </cim:Substation>")
    println(io, "  <cim:VoltageLevel rdf:ID=\"", _igm_id("vl", bus), "\">")
    _igm_name(io, string(bus, "_110"))
    println(io, "    <cim:VoltageLevel.Substation rdf:resource=\"#", _igm_id("ss", bus), "\"/>")
    println(io, "    <cim:VoltageLevel.BaseVoltage rdf:resource=\"#_bv110\"/>")
    println(io, "  </cim:VoltageLevel>")
  end
  if op == "Nord"
    println(io, "  <cim:ExternalNetworkInjection rdf:ID=\"_eni_Holm\">")
    _igm_name(io, "Grid_Holm")
    println(io, "    <cim:ExternalNetworkInjection.referencePriority>1</cim:ExternalNetworkInjection.referencePriority>")
    println(io, "  </cim:ExternalNetworkInjection>")
    _igm_terminal(io, "eni_Holm", 1)
  end
  for (bus, _) in sort!(collect(_IGM_LOADS))
    any(b -> b[1] == bus && b[2] == op, _IGM_BUSES) || continue
    println(io, "  <cim:EnergyConsumer rdf:ID=\"", _igm_id("load", bus), "\">")
    _igm_name(io, string("Load_", bus))
    println(io, "  </cim:EnergyConsumer>")
    _igm_terminal(io, string("load_", bus), 1)
  end
  bch = 2pi * 50.0 * 9.5e-9
  for (name, o, _, _, km) in _IGM_LINES
    o == op || continue
    println(io, "  <cim:ACLineSegment rdf:ID=\"", _igm_id(name), "\">")
    _igm_name(io, name)
    println(io, "    <cim:ACLineSegment.r>", _igm_fmt(0.12 * km), "</cim:ACLineSegment.r>")
    println(io, "    <cim:ACLineSegment.x>", _igm_fmt(0.4 * km), "</cim:ACLineSegment.x>")
    println(io, "    <cim:ACLineSegment.bch>", _igm_fmt(bch * km), "</cim:ACLineSegment.bch>")
    println(io, "    <cim:ConductingEquipment.BaseVoltage rdf:resource=\"#_bv110\"/>")
    println(io, "  </cim:ACLineSegment>")
    _igm_terminal(io, name, 1)
    _igm_terminal(io, name, 2)
  end
  println(io, "  <cim:ControlArea rdf:ID=\"", _igm_id("ca", op), "\">")
  _igm_name(io, string("CA_", op))
  println(io, "    <cim:ControlArea.type rdf:resource=\"http://iec.ch/TC57/2013/CIM-schema-cim16#ControlAreaTypeKind.Interchange\"/>")
  println(io, "  </cim:ControlArea>")
  for (name, o, f, t, _) in _IGM_LINES
    o == op || continue
    for (seq, bus) in ((1, f), (2, t))
      bus in _IGM_XNODES || continue
      println(io, "  <cim:TieFlow rdf:ID=\"", _igm_id("tf", name), "\">")
      println(io, "    <cim:TieFlow.ControlArea rdf:resource=\"#", _igm_id("ca", op), "\"/>")
      println(io, "    <cim:TieFlow.positiveFlowIn>true</cim:TieFlow.positiveFlowIn>")
      println(io, "    <cim:TieFlow.Terminal rdf:resource=\"#", _igm_id("t", name, seq), "\"/>")
      println(io, "  </cim:TieFlow>")
    end
  end
  println(io, "</rdf:RDF>")
  return String(take!(io))
end

function _igm_tp(op::AbstractString)::String
  io = IOBuffer()
  _igm_head(io, :TP, op)
  for (bus, o, _) in _IGM_BUSES
    o == op || continue
    println(io, "  <cim:TopologicalNode rdf:ID=\"", _igm_tn(bus), "\">")
    _igm_name(io, bus)
    println(io, "    <cim:TopologicalNode.BaseVoltage rdf:resource=\"#_bv110\"/>")
    println(io, "    <cim:TopologicalNode.ConnectivityNodeContainer rdf:resource=\"#", _igm_id("vl", bus), "\"/>")
    println(io, "  </cim:TopologicalNode>")
  end
  term(eq, seq, bus) = println(io, "  <cim:Terminal rdf:about=\"#", _igm_id("t", eq, seq), "\">\n    <cim:Terminal.TopologicalNode rdf:resource=\"#", _igm_tn(bus), "\"/>\n  </cim:Terminal>")
  op == "Nord" && term("eni_Holm", 1, "Holm")
  for (bus, _) in sort!(collect(_IGM_LOADS))
    any(b -> b[1] == bus && b[2] == op, _IGM_BUSES) && term(string("load_", bus), 1, bus)
  end
  for (name, o, f, t, _) in _IGM_LINES
    o == op || continue
    term(name, 1, f)
    term(name, 2, t)
  end
  println(io, "</rdf:RDF>")
  return String(take!(io))
end

function _igm_ssh(op::AbstractString)::String
  io = IOBuffer()
  _igm_head(io, :SSH, op)
  if op == "Nord"
    println(io, "  <cim:ExternalNetworkInjection rdf:about=\"#_eni_Holm\">")
    println(io, "    <cim:ExternalNetworkInjection.p>-70</cim:ExternalNetworkInjection.p>")
    println(io, "    <cim:ExternalNetworkInjection.q>-23</cim:ExternalNetworkInjection.q>")
    println(io, "  </cim:ExternalNetworkInjection>")
  end
  for (bus, (p, q)) in sort!(collect(_IGM_LOADS))
    any(b -> b[1] == bus && b[2] == op, _IGM_BUSES) || continue
    println(io, "  <cim:EnergyConsumer rdf:about=\"#", _igm_id("load", bus), "\">")
    println(io, "    <cim:EnergyConsumer.p>", _igm_fmt(p), "</cim:EnergyConsumer.p>")
    println(io, "    <cim:EnergyConsumer.q>", _igm_fmt(q), "</cim:EnergyConsumer.q>")
    println(io, "  </cim:EnergyConsumer>")
  end
  eqs = String[]
  op == "Nord" && push!(eqs, "eni_Holm")
  append!(eqs, [string("load_", b) for (b, o, _) in _IGM_BUSES if o == op && haskey(_IGM_LOADS, b)])
  for e in eqs
    println(io, "  <cim:Terminal rdf:about=\"#", _igm_id("t", e, 1), "\">\n    <cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected>\n  </cim:Terminal>")
  end
  for (name, o, _, _, _) in _IGM_LINES
    o == op || continue
    for seq in 1:2
      println(io, "  <cim:Terminal rdf:about=\"#", _igm_id("t", name, seq), "\">\n    <cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected>\n  </cim:Terminal>")
    end
  end
  println(io, "</rdf:RDF>")
  return String(take!(io))
end

function _igm_eq_bd()::String
  io = IOBuffer()
  _igm_head(io, :EQ_BD, "BD")
  println(io, "  <cim:BaseVoltage rdf:ID=\"_bv110\">\n    <cim:BaseVoltage.nominalVoltage>110</cim:BaseVoltage.nominalVoltage>\n  </cim:BaseVoltage>")
  for x in _IGM_XNODES
    println(io, "  <cim:Line rdf:ID=\"", _igm_id("bdline", x), "\">")
    _igm_name(io, string("Line_", x))
    println(io, "  </cim:Line>")
    println(io, "  <cim:ConnectivityNode rdf:ID=\"", _igm_id("cn", x), "\">")
    _igm_name(io, x)
    println(io, "    <cim:ConnectivityNode.ConnectivityNodeContainer rdf:resource=\"#", _igm_id("bdline", x), "\"/>")
    println(io, "  </cim:ConnectivityNode>")
  end
  println(io, "</rdf:RDF>")
  return String(take!(io))
end

function _igm_tp_bd()::String
  io = IOBuffer()
  _igm_head(io, :TP_BD, "BD")
  for x in _IGM_XNODES
    println(io, "  <cim:TopologicalNode rdf:ID=\"", _igm_tn(x), "\">")
    _igm_name(io, x)
    println(io, "    <cim:TopologicalNode.BaseVoltage rdf:resource=\"#_bv110\"/>")
    println(io, "    <cim:TopologicalNode.ConnectivityNodeContainer rdf:resource=\"#", _igm_id("bdline", x), "\"/>")
    println(io, "  </cim:TopologicalNode>")
  end
  println(io, "</rdf:RDF>")
  return String(take!(io))
end

function build_igm_fixture(case::AbstractString)
  outdir = joinpath(_FIXTURE_ROOT, case)
  mkpath(outdir)
  for f in readdir(outdir; join = true)
    rm(f)
  end
  files = [
    "$(case)_Nord_EQ.xml" => _igm_eq("Nord"), "$(case)_Nord_TP.xml" => _igm_tp("Nord"), "$(case)_Nord_SSH.xml" => _igm_ssh("Nord"),
    "$(case)_Sued_EQ.xml" => _igm_eq("Sued"), "$(case)_Sued_TP.xml" => _igm_tp("Sued"), "$(case)_Sued_SSH.xml" => _igm_ssh("Sued"),
    "$(case)_BD_EQ.xml" => _igm_eq_bd(), "$(case)_BD_TP.xml" => _igm_tp_bd(),
  ]
  total = 0
  println(case, ": hand-written CGM, two operator models and a boundary set")
  for (name, text) in files
    write(joinpath(outdir, name), text)
    total += sizeof(text)
    @printf("  %-26s %9d bytes\n", name, sizeof(text))
  end
  @printf("  %-26s %9d bytes\n", "total", total)
  return total
end

const CGMES_HANDWRITTEN_CASES = ["sp_case6_igm" => build_igm_fixture]

function main(args)
  known = vcat(first.(CGMES_FIXTURE_CASES), first.(CGMES_HANDWRITTEN_CASES))
  selected = isempty(args) ? known : String.(args)
  for (case, loader) in CGMES_FIXTURE_CASES
    case in selected || continue
    build_cgmes_fixture(case, loader)
  end
  for (case, builder) in CGMES_HANDWRITTEN_CASES
    case in selected && builder(case)
  end
  unknown = setdiff(selected, known)
  isempty(unknown) || error(string("unknown fixture case(s): ", join(unknown, ", ")))
  return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
  Base.invokelatest(main, ARGS)
end

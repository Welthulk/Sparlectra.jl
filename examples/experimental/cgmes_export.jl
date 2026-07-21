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

# file: examples/experimental/cgmes_export.jl
# purpose: experimental Stage-1 CGMES 2.4.15 (RDF/XML) exporter for a Sparlectra Net.
#
# Experimental example tooling. Not part of the stable Sparlectra API, not
# wired into the Web UI, and may change without notice.
#
# Stage 1 scope: EQ + TP profiles for buses (TopologicalNode) and
# ACLineSegments, including optional zero-sequence/short-circuit attributes
# (EquipmentShortCircuit).
#
# Deliberately not included (later stages):
#   - Loads/generators/shunts (EnergyConsumer, SynchronousMachine, ...) + SSH
#   - Transformers (PowerTransformerEnd, tap changer)
#   - SV profile (SvVoltage/SvPowerFlow from the power-flow result)
#   - ZIP container
#
# Accuracy note: class/attribute names follow CGMES 2.4.15 (cim16). Validate
# against the ENTSO-E test models (e.g. MicroGrid) before any production use;
# individual profiles may require additional mandatory fields.

import UUIDs
import Dates

# Fixed namespace for deterministic mRIDs (uuid5). Same network + same
# component names => identical UUIDs on every export (roundtrip-/diff-stable).
const CGMES_UUID_NAMESPACE = UUIDs.UUID("6ba7b810-9dad-11d1-80b4-00c04fd430c8")

const CIM_NS = "http://iec.ch/TC57/2013/CIM-schema-cim16#"
const RDF_NS = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
const MD_NS = "http://iec.ch/TC57/61970-552/ModelDescription/1#"

const EQ_PROFILES = [
  "http://entsoe.eu/CIM/EquipmentCore/3/1",
  "http://entsoe.eu/CIM/EquipmentShortCircuit/3/1",
  "http://entsoe.eu/CIM/EquipmentOperation/3/1",
]
const TP_PROFILES = ["http://entsoe.eu/CIM/Topology/4/1"]

"""
    CGMESLineShortCircuit

Zero-sequence/short-circuit data for one line, for the EquipmentShortCircuit
profile. Physical units (Ohm, Siemens, °C). Only written when supplied — no
zero-sequence values are invented from positive-sequence data.
"""
Base.@kwdef struct CGMESLineShortCircuit
  r0_ohm::Float64
  x0_ohm::Float64
  b0ch_S::Float64 = 0.0
  g0ch_S::Float64 = 0.0
  endTemperature_C::Union{Nothing,Float64} = nothing
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

cgmesUUID(kind::AbstractString, name::AbstractString) = string(UUIDs.uuid5(CGMES_UUID_NAMESPACE, string(kind, "|", name)))

function xmlEscape(s::AbstractString)
  s = replace(s, "&" => "&amp;")
  s = replace(s, "<" => "&lt;")
  s = replace(s, ">" => "&gt;")
  s = replace(s, "\"" => "&quot;")
  return s
end

fmtVal(x::Real) = string(Float64(x))

function writeXmlHeader(io::IO)
  println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
  println(io, "<rdf:RDF xmlns:cim=\"$(CIM_NS)\" xmlns:md=\"$(MD_NS)\" xmlns:rdf=\"$(RDF_NS)\">")
end

function writeFullModel(io::IO, modelId::AbstractString, profiles::Vector{String}, dependentOn::Vector{String})
  scenarioTime = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SSZ")
  println(io, "  <md:FullModel rdf:about=\"urn:uuid:$(modelId)\">")
  println(io, "    <md:Model.scenarioTime>$(scenarioTime)</md:Model.scenarioTime>")
  println(io, "    <md:Model.created>$(scenarioTime)</md:Model.created>")
  println(io, "    <md:Model.version>1</md:Model.version>")
  for p in profiles
    println(io, "    <md:Model.profile>$(p)</md:Model.profile>")
  end
  for d in dependentOn
    println(io, "    <md:Model.DependentOn rdf:resource=\"urn:uuid:$(d)\"/>")
  end
  println(io, "    <md:Model.modelingAuthoritySet>http://sparlectra.local/</md:Model.modelingAuthoritySet>")
  println(io, "  </md:FullModel>")
end

# Shared collection context so the EQ and TP profiles use the same mRIDs.
struct CGMESContext
  netName::String
  idx2busName::Dict{Int,String}
  baseVoltageIds::Dict{Float64,String}    # vn_kV -> mRID BaseVoltage
  voltageLevelIds::Dict{Float64,String}   # vn_kV -> mRID VoltageLevel
  substationId::String
  regionId::String
  subRegionId::String
  topoNodeIds::Dict{Int,String}           # busIdx -> mRID TopologicalNode
  lineTerminalIds::Vector{NTuple{2,String}} # per line (Terminal1, Terminal2)
  lineIds::Vector{String}                 # mRID ACLineSegment per net.linesAC index
  lineContainerIds::Vector{String}        # mRID cim:Line container per line
end

function buildContext(net)
  idx2busName = Dict{Int,String}(v => k for (k, v) in net.busDict)

  baseVoltageIds = Dict{Float64,String}()
  voltageLevelIds = Dict{Float64,String}()
  for node in net.nodeVec
    vn = node.comp.cVN
    if !haskey(baseVoltageIds, vn)
      baseVoltageIds[vn] = cgmesUUID("BaseVoltage", string(vn))
      voltageLevelIds[vn] = cgmesUUID("VoltageLevel", string(net.name, "|", vn))
    end
  end

  topoNodeIds = Dict{Int,String}()
  for node in net.nodeVec
    busName = get(idx2busName, node.busIdx, string("Bus_", node.busIdx))
    topoNodeIds[node.busIdx] = cgmesUUID("TopologicalNode", busName)
  end

  lineTerminalIds = NTuple{2,String}[]
  lineIds = String[]
  lineContainerIds = String[]
  for line in net.linesAC
    lname = line.comp.cName
    push!(lineIds, cgmesUUID("ACLineSegment", lname))
    push!(lineContainerIds, cgmesUUID("Line", lname))
    push!(lineTerminalIds, (cgmesUUID("Terminal", string(lname, "|T1")), cgmesUUID("Terminal", string(lname, "|T2"))))
  end

  return CGMESContext(
    net.name,
    idx2busName,
    baseVoltageIds,
    voltageLevelIds,
    cgmesUUID("Substation", net.name),
    cgmesUUID("GeographicalRegion", net.name),
    cgmesUUID("SubGeographicalRegion", net.name),
    topoNodeIds,
    lineTerminalIds,
    lineIds,
    lineContainerIds,
  )
end

# Always deliver line parameters in physical units (Ohm/S). PI-model lines
# carry p.u. values and are converted back.
function lineParamsOhm(net, line)
  if line._isPIModel
    vn = line.comp.cVN
    g_pu = isnothing(line.g) ? 0.0 : line.g
    b_pu = isnothing(line.b) ? 0.0 : line.b
    r, x, g, b = fromPU_RXBG(r_pu = line.r, x_pu = line.x, g_pu = g_pu, b_pu = b_pu, v_kv = vn, baseMVA = net.baseMVA)
    return r, x, b, g
  else
    r, x, b, g = getLineRXBG(line)
    return r, x, (isnothing(b) ? 0.0 : b), (isnothing(g) ? 0.0 : g)
  end
end

# ---------------------------------------------------------------------------
# EQ profile
# ---------------------------------------------------------------------------

function writeEQFile(net, ctx::CGMESContext, path::AbstractString, eqModelId::AbstractString; sc_line_data::Dict{Int,CGMESLineShortCircuit} = Dict{Int,CGMESLineShortCircuit}())
  open(path, "w") do io
    writeXmlHeader(io)
    writeFullModel(io, eqModelId, EQ_PROFILES, String[])

    # Regions / substation / voltage levels (minimal container hierarchy).
    println(io, "  <cim:GeographicalRegion rdf:ID=\"_$(ctx.regionId)\">")
    println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(ctx.netName))_Region</cim:IdentifiedObject.name>")
    println(io, "  </cim:GeographicalRegion>")

    println(io, "  <cim:SubGeographicalRegion rdf:ID=\"_$(ctx.subRegionId)\">")
    println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(ctx.netName))_SubRegion</cim:IdentifiedObject.name>")
    println(io, "    <cim:SubGeographicalRegion.Region rdf:resource=\"#_$(ctx.regionId)\"/>")
    println(io, "  </cim:SubGeographicalRegion>")

    println(io, "  <cim:Substation rdf:ID=\"_$(ctx.substationId)\">")
    println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(ctx.netName))_Substation</cim:IdentifiedObject.name>")
    println(io, "    <cim:Substation.Region rdf:resource=\"#_$(ctx.subRegionId)\"/>")
    println(io, "  </cim:Substation>")

    for (vn, bvId) in sort(collect(ctx.baseVoltageIds); by = first)
      println(io, "  <cim:BaseVoltage rdf:ID=\"_$(bvId)\">")
      println(io, "    <cim:IdentifiedObject.name>BV_$(vn)_kV</cim:IdentifiedObject.name>")
      println(io, "    <cim:BaseVoltage.nominalVoltage>$(fmtVal(vn))</cim:BaseVoltage.nominalVoltage>")
      println(io, "  </cim:BaseVoltage>")

      vlId = ctx.voltageLevelIds[vn]
      println(io, "  <cim:VoltageLevel rdf:ID=\"_$(vlId)\">")
      println(io, "    <cim:IdentifiedObject.name>VL_$(vn)_kV</cim:IdentifiedObject.name>")
      println(io, "    <cim:VoltageLevel.Substation rdf:resource=\"#_$(ctx.substationId)\"/>")
      println(io, "    <cim:VoltageLevel.BaseVoltage rdf:resource=\"#_$(bvId)\"/>")
      println(io, "  </cim:VoltageLevel>")
    end

    # ACLineSegments + terminals.
    for (i, line) in enumerate(net.linesAC)
      lname = line.comp.cName
      vn = line.comp.cVN
      r, x, b, g = lineParamsOhm(net, line)

      println(io, "  <cim:Line rdf:ID=\"_$(ctx.lineContainerIds[i])\">")
      println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(lname))_Line</cim:IdentifiedObject.name>")
      println(io, "    <cim:Line.Region rdf:resource=\"#_$(ctx.subRegionId)\"/>")
      println(io, "  </cim:Line>")

      println(io, "  <cim:ACLineSegment rdf:ID=\"_$(ctx.lineIds[i])\">")
      println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(lname))</cim:IdentifiedObject.name>")
      println(io, "    <cim:Equipment.EquipmentContainer rdf:resource=\"#_$(ctx.lineContainerIds[i])\"/>")
      println(io, "    <cim:ConductingEquipment.BaseVoltage rdf:resource=\"#_$(ctx.baseVoltageIds[vn])\"/>")
      println(io, "    <cim:Conductor.length>$(fmtVal(line.length))</cim:Conductor.length>")
      println(io, "    <cim:ACLineSegment.r>$(fmtVal(r))</cim:ACLineSegment.r>")
      println(io, "    <cim:ACLineSegment.x>$(fmtVal(x))</cim:ACLineSegment.x>")
      println(io, "    <cim:ACLineSegment.bch>$(fmtVal(b))</cim:ACLineSegment.bch>")
      println(io, "    <cim:ACLineSegment.gch>$(fmtVal(g))</cim:ACLineSegment.gch>")

      # Short-circuit/zero-sequence data (EquipmentShortCircuit) — only when supplied.
      if haskey(sc_line_data, i)
        sc = sc_line_data[i]
        println(io, "    <cim:ACLineSegment.r0>$(fmtVal(sc.r0_ohm))</cim:ACLineSegment.r0>")
        println(io, "    <cim:ACLineSegment.x0>$(fmtVal(sc.x0_ohm))</cim:ACLineSegment.x0>")
        println(io, "    <cim:ACLineSegment.b0ch>$(fmtVal(sc.b0ch_S))</cim:ACLineSegment.b0ch>")
        println(io, "    <cim:ACLineSegment.g0ch>$(fmtVal(sc.g0ch_S))</cim:ACLineSegment.g0ch>")
        if !isnothing(sc.endTemperature_C)
          println(io, "    <cim:ACLineSegment.shortCircuitEndTemperature>$(fmtVal(sc.endTemperature_C))</cim:ACLineSegment.shortCircuitEndTemperature>")
        end
      end
      println(io, "  </cim:ACLineSegment>")

      # Terminals (EQ side: equipment association + sequence order).
      for (seq, termId) in enumerate(ctx.lineTerminalIds[i])
        println(io, "  <cim:Terminal rdf:ID=\"_$(termId)\">")
        println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(lname))_T$(seq)</cim:IdentifiedObject.name>")
        println(io, "    <cim:Terminal.ConductingEquipment rdf:resource=\"#_$(ctx.lineIds[i])\"/>")
        println(io, "    <cim:ACDCTerminal.sequenceNumber>$(seq)</cim:ACDCTerminal.sequenceNumber>")
        println(io, "  </cim:Terminal>")
      end
    end

    println(io, "</rdf:RDF>")
  end
  return path
end

# ---------------------------------------------------------------------------
# TP profile
# ---------------------------------------------------------------------------

function writeTPFile(net, ctx::CGMESContext, path::AbstractString, tpModelId::AbstractString, eqModelId::AbstractString)
  open(path, "w") do io
    writeXmlHeader(io)
    writeFullModel(io, tpModelId, TP_PROFILES, [eqModelId])

    # TopologicalNodes.
    for node in net.nodeVec
      busName = get(ctx.idx2busName, node.busIdx, string("Bus_", node.busIdx))
      vn = node.comp.cVN
      tnId = ctx.topoNodeIds[node.busIdx]
      println(io, "  <cim:TopologicalNode rdf:ID=\"_$(tnId)\">")
      println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(busName))</cim:IdentifiedObject.name>")
      println(io, "    <cim:TopologicalNode.BaseVoltage rdf:resource=\"#_$(ctx.baseVoltageIds[vn])\"/>")
      println(io, "    <cim:TopologicalNode.ConnectivityNodeContainer rdf:resource=\"#_$(ctx.voltageLevelIds[vn])\"/>")
      println(io, "  </cim:TopologicalNode>")
    end

    # Terminal -> TopologicalNode association.
    for (i, line) in enumerate(net.linesAC)
      fromIdx = line.comp.cFrom_bus
      toIdx = line.comp.cTo_bus
      busIdxs = (fromIdx, toIdx)
      for (seq, termId) in enumerate(ctx.lineTerminalIds[i])
        busIdx = busIdxs[seq]
        if isnothing(busIdx) || !haskey(ctx.topoNodeIds, busIdx)
          @warn "CGMES-TP: line $(line.comp.cName) terminal $(seq) has no valid bus index — skipped"
          continue
        end
        println(io, "  <cim:Terminal rdf:about=\"#_$(termId)\">")
        println(io, "    <cim:Terminal.TopologicalNode rdf:resource=\"#_$(ctx.topoNodeIds[busIdx])\"/>")
        println(io, "    <cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected>")
        println(io, "  </cim:Terminal>")
      end
    end

    println(io, "</rdf:RDF>")
  end
  return path
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

"""
    writeCGMESFiles(net; path::AbstractString = pwd(),
                    sc_line_data::Dict{Int,CGMESLineShortCircuit} = Dict())

Exports a Sparlectra `Net` as CGMES-2.4.15 profile files (Stage 1: EQ + TP,
buses and ACLineSegments). `sc_line_data` supplies optional zero-sequence/
short-circuit data per line index (order of `net.linesAC`).

Works for any `Net` — built programmatically, or imported from MATPOWER or
DTF — since the export reads directly from the network model.

Returns a vector of the written file paths.

Not covered (later stages): transformers, prosumers/SSH, SV, ZIP container.
Existing transformers/prosumers/shunts in the network are skipped with a
warning.
"""
function writeCGMESFiles(net; path::AbstractString = pwd(), sc_line_data::Dict{Int,CGMESLineShortCircuit} = Dict{Int,CGMESLineShortCircuit}())
  isempty(net.trafos) || @warn "CGMES export Stage 1: $(length(net.trafos)) transformer(s) are not yet exported"
  isempty(net.prosumpsVec) || @warn "CGMES export Stage 1: $(length(net.prosumpsVec)) prosumer(s) are not yet exported (SSH follows in a later stage)"
  isempty(net.shuntVec) || @warn "CGMES export Stage 1: $(length(net.shuntVec)) shunt(s) are not yet exported"

  ctx = buildContext(net)
  eqModelId = cgmesUUID("Model", string(net.name, "|EQ"))
  tpModelId = cgmesUUID("Model", string(net.name, "|TP"))

  eqPath = joinpath(path, string(net.name, "_EQ.xml"))
  tpPath = joinpath(path, string(net.name, "_TP.xml"))

  writeEQFile(net, ctx, eqPath, eqModelId; sc_line_data = sc_line_data)
  writeTPFile(net, ctx, tpPath, tpModelId, eqModelId)

  return [eqPath, tpPath]
end

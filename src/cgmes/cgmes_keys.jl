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

# file: src/cgmes/cgmes_keys.jl
# purpose: structural identity keys for `net.cgmes_ids`, shared by the CGMES
#          importer (mRID capture) and the CGMES exporter (mRID reuse/minting).
#
# A structural key describes WHAT an object is in the network — not what it is
# called in a source file. The importer records original mRIDs under these
# keys; the exporter looks the same keys up and reuses the recorded mRID, so a
# CGMES → Net → CGMES roundtrip preserves object identity. Keys that are not
# recorded (programmatic nets, MATPOWER imports) are minted deterministically
# via uuid5 over the key string itself.
#
# Key syntax (one line per object class):
#   TN|<busname>                     TopologicalNode of a bus
#   BV|<vn_kV>                       BaseVoltage per nominal voltage
#   VL|<vn_kV>                       VoltageLevel per nominal voltage
#   ACL|<busA>|<busB>|<k>            ACLineSegment; bus pair sorted
#                                    lexicographically, k = 1-based first-seen
#                                    index for parallel lines on the same pair
#   LNC|<busA>|<busB>|<k>            cim:Line container of that segment
#   PT|<busA>|<busB>|<k>             PowerTransformer (same pair/k rules,
#                                    counted over transformer branches only)
#   <ptkey>|E1 / |E2                 PowerTransformerEnd by end number
#   <ptkey>|PTC                      single-step linear phase tap changer
#                                    carrying the branch phase shift
#   EC|<bus>|<k>                     EnergyConsumer; k = 1-based first-seen
#                                    index per class and bus
#   SM|<bus>|<k>                     SynchronousMachine (+ "|GU", "|RC")
#   ENI|<bus>|<k>                    ExternalNetworkInjection (+ "|RC")
#   ASM|<bus>|<k>                    AsynchronousMachine
#   SH|<bus>|<k>                     LinearShuntCompensator
#   <parentkey>|T1 / |T2             Terminal by equipment sequence number
#   MODEL|EQ / MODEL|TP / MODEL|SSH  md:FullModel header ids per profile
#   REGION / SUBREGION / SUBSTATION  the single container hierarchy objects

# A single leading underscore in an mRID is RDF/XML id syntax
# (rdf:ID="_<uuid>"), not part of the model identity. Captured mRIDs are
# stored without it; the exporter re-adds the underscore when writing rdf:ID.
cgmesCanonicalMrid(mrid::AbstractString)::String = startswith(mrid, "_") ? String(mrid[2:end]) : String(mrid)

cgmesKeyTopologicalNode(bus::AbstractString)::String = string("TN|", bus)
cgmesKeyBaseVoltage(vn_kV::Real)::String = string("BV|", Float64(vn_kV))
cgmesKeyVoltageLevel(vn_kV::Real)::String = string("VL|", Float64(vn_kV))
cgmesKeyModel(profile::AbstractString)::String = string("MODEL|", profile)
cgmesKeyTerminal(parentKey::AbstractString, seq::Integer)::String = string(parentKey, "|T", seq)

function _cgmesOrderedPair(busA::AbstractString, busB::AbstractString)
  return busA <= busB ? (String(busA), String(busB)) : (String(busB), String(busA))
end

function cgmesKeyACLineSegment(busA::AbstractString, busB::AbstractString, k::Integer)::String
  a, b = _cgmesOrderedPair(busA, busB)
  return string("ACL|", a, "|", b, "|", k)
end

function cgmesKeyLineContainer(busA::AbstractString, busB::AbstractString, k::Integer)::String
  a, b = _cgmesOrderedPair(busA, busB)
  return string("LNC|", a, "|", b, "|", k)
end

function cgmesKeyPowerTransformer(busA::AbstractString, busB::AbstractString, k::Integer)::String
  a, b = _cgmesOrderedPair(busA, busB)
  return string("PT|", a, "|", b, "|", k)
end

# normalized pair as the "bus" part of the shared (kind, bus) counter, so PT
# counting works through cgmesNextBusEquipmentIndex! like the bus classes
function _cgmesTrafoPairKey(busA::AbstractString, busB::AbstractString)::String
  a, b = _cgmesOrderedPair(busA, busB)
  return string(a, "|", b)
end

# three-winding transformer: key over the sorted side buses (the star bus is
# an artifact of the equivalent, not part of the identity); end/terminal
# children use the source end numbering
function cgmesKeyPowerTransformer3W(busA::AbstractString, busB::AbstractString, busC::AbstractString, k::Integer)::String
  s = sort([String(busA), String(busB), String(busC)])
  return string("PT3|", s[1], "|", s[2], "|", s[3], "|", k)
end

cgmesKeyTransformerEnd(ptKey::AbstractString, endNumber::Integer)::String = string(ptKey, "|E", endNumber)

# single-terminal equipment on a bus: EC (EnergyConsumer), SM
# (SynchronousMachine), ENI (ExternalNetworkInjection), ASM
# (AsynchronousMachine), SH (LinearShuntCompensator)
cgmesKeyBusEquipment(kind::AbstractString, bus::AbstractString, k::Integer)::String = string(kind, "|", bus, "|", k)

# First-seen 1-based index per (class kind, bus). Same consistency argument as
# the line-pair counter: importer capture and exporter walk the equipment of
# one class in the same creation order.
function cgmesNextBusEquipmentIndex!(counter::Dict{Tuple{String,String},Int}, kind::AbstractString, bus::AbstractString)::Int
  key = (String(kind), String(bus))
  k = get(counter, key, 0) + 1
  counter[key] = k
  return k
end

# First-seen 1-based parallel index for a bus pair. Importer and exporter walk
# the lines in the same order (import order == net.linesAC order), so both
# sides assign identical <k> values without storing the counter anywhere.
function cgmesNextParallelIndex!(counter::Dict{Tuple{String,String},Int}, busA::AbstractString, busB::AbstractString)::Int
  pair = _cgmesOrderedPair(busA, busB)
  k = get(counter, pair, 0) + 1
  counter[pair] = k
  return k
end

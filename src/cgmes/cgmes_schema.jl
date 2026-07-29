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

# file: src/cgmes/cgmes_schema.jl
# purpose: shared CGMES schema tables (namespaces, profile URIs, version
# detection) used by both the importer and the export path. Class and
# attribute names exist exactly once, here.

# --- namespaces -------------------------------------------------------------

const RDF_NS = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
const MD_NS = "http://iec.ch/TC57/61970-552/ModelDescription/1#"
const DM_NS = "http://iec.ch/TC57/61970-552/DifferenceModel/1#"

# CIM data namespace per CGMES version.
const CIM16_NS = "http://iec.ch/TC57/2013/CIM-schema-cim16#"   # CGMES 2.4.15
const CIM100_NS = "http://iec.ch/TC57/CIM100#"                 # CGMES 3.0

"""CGMES version string derived from the `cim` namespace URI; `""` if unknown."""
function cgmesVersionFromNamespace(cim_uri::AbstractString)::String
  cim_uri == CIM16_NS && return "2.4.15"
  cim_uri == CIM100_NS && return "3.0"
  return ""
end

# --- profile classification -------------------------------------------------

# `md:Model.profile` URI → profile tag. One file may declare several URIs
# (A2 finding: MicroGrid EQ declares EquipmentCore AND EquipmentShortCircuit),
# so classification always yields a Set of tags per file.
#
# CGMES 2.4.15 profile URIs carry entsoe.eu hosts; CGMES 3.0 uses
# iec.ch/ns URIs with "-EU" suffixes. Substring keys keep the table short and
# version-agnostic where safe.
const PROFILE_URI_TAGS = [
  ("EquipmentBoundaryOperation", :EQ_BD_OP),
  ("EquipmentBoundary", :EQ_BD),
  ("TopologyBoundary", :TP_BD),
  ("EquipmentShortCircuit", :EQ_SC),
  ("EquipmentOperation", :EQ_OP),
  ("CoreEquipment", :EQ),          # CGMES 3.0 naming
  ("EquipmentCore", :EQ),          # CGMES 2.4.15 naming
  ("SteadyStateHypothesis", :SSH),
  ("StateVariables", :SV),
  ("Topology", :TP),               # after TopologyBoundary!
  ("DiagramLayout", :DL),
  ("GeographicalLocation", :GL),
  ("Dynamics", :DY),
]

"""Map one `md:Model.profile` URI to a profile tag (`:EQ`, `:TP`, …); `:UNKNOWN` if unmatched."""
function profileTagFromURI(uri::AbstractString)::Symbol
  for (needle, tag) in PROFILE_URI_TAGS
    occursin(needle, uri) && return tag
  end
  return :UNKNOWN
end

# CGMES 3.0 / Network Codes deliveries may carry a `dcat:Dataset` document
# header instead of `md:FullModel` (seen in ENTSO-E's ReliCapGrid boundary
# files). There the profile is declared as a short code in `dcat:keyword`,
# not as a Model.profile URI. Keys are the codes ENTSO-E uses in that header.
const PROFILE_KEYWORD_TAGS = Dict(
  "EQ" => :EQ,
  "OP" => :EQ_OP,
  "SC" => :EQ_SC,
  "SSH" => :SSH,
  "TP" => :TP,
  "SV" => :SV,
  "BD" => :EQ_BD,      # ReliCapGrid: one combined boundary file per border
  "EQBD" => :EQ_BD,
  "TPBD" => :TP_BD,
  "DL" => :DL,
  "GL" => :GL,
  "DY" => :DY,
)

"""Map one `dcat:keyword` short code to a profile tag; `:UNKNOWN` if unmatched."""
profileTagFromKeyword(kw::AbstractString)::Symbol = get(PROFILE_KEYWORD_TAGS, uppercase(strip(String(kw))), :UNKNOWN)

# Profile tags that carry power-flow relevant data (Stage 0 reads only these;
# DL/GL/DY are counted and skipped).
const IMPORT_PROFILE_TAGS = Set([:EQ, :EQ_OP, :EQ_SC, :SSH, :TP, :SV, :EQ_BD, :EQ_BD_OP, :TP_BD])

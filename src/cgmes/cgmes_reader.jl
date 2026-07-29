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

# file: src/cgmes/cgmes_reader.jl
# purpose: layer 2 of the CGMES importer — generic RDF/XML property-bag
# reader. Every top-level `<cim:Class rdf:ID|about>` becomes/updates one
# CIMObject; no per-class code. `rdf:ID` creates, `rdf:about` overlays
# (SSH/TP/SV on EQ). Verified against the ENTSO-E test sets (A2): flat
# structure, enum values as rdf:resource fragments, inherited attributes
# with foreign class prefixes.

"""
One CIM instance as a property bag. `attrs` holds literal values and enum
fragments, `refs` holds mRID targets of `rdf:resource="#…"` references.
Keys are the property local names without their class prefix
(`ACLineSegment.r` → `:r`, `IdentifiedObject.name` → `:name`); when two
different property names share a suffix inside one object, the full dotted
name is kept as an additional key.

`source` is the name of the file whose `rdf:ID` created the object (overlays
via `rdf:about` do not change it). In a multi-area delivery this is the only
way to tell which area contributed a piece of equipment — the sides of an
assembled border are otherwise indistinguishable.
"""
struct CIMObject
  mrid::String
  class::Symbol
  attrs::Dict{Symbol,String}
  refs::Dict{Symbol,String}
  source::String
end

"""Per-file metadata from the md:Model header (layer-1/2 handshake)."""
struct CGMESFileInfo
  name::String
  header::Symbol                 # :FullModel | :DifferenceModel | :Dataset (dcat) | :none
  profile_uris::Vector{String}
  profiles::Set{Symbol}          # mapped tags, see cgmes_schema.jl
  version::String                # "2.4.15" | "3.0" | ""
  object_count::Int
  skipped::Bool
  skip_reason::String
end

_strip_mrid(v::AbstractString) = startswith(v, "#") ? String(v[2:end]) : String(v)

function _property_key!(attrs_or_refs::Dict{Symbol,String}, pname::AbstractString, value::String)
  key = Symbol(last(split(pname, '.')))
  if haskey(attrs_or_refs, key) && attrs_or_refs[key] != value
    # suffix collision inside one object → keep the full dotted name too
    attrs_or_refs[Symbol(pname)] = value
  else
    attrs_or_refs[key] = value
  end
  return nothing
end

"""
    readCGMESFile!(objects, byclass, file; profile_filter=IMPORT_PROFILE_TAGS) -> CGMESFileInfo

Parse one CGMES XML payload into the shared object store. DifferenceModel
files and files outside `profile_filter` are classified but not parsed
(`skipped = true`). Returns the file metadata.
"""
function readCGMESFile!(objects::Dict{String,CIMObject}, byclass::Dict{Symbol,Vector{String}}, file::CGMESFile; profile_filter::Set{Symbol} = IMPORT_PROFILE_TAGS)::CGMESFileInfo
  doc = EzXML.parsexml(file.content)
  r = EzXML.root(doc)
  cim_uri = ""
  for (p, uri) in EzXML.namespaces(r)
    p == "cim" && (cim_uri = uri)
  end
  version = cgmesVersionFromNamespace(cim_uri)

  header = :none
  profile_uris = String[]
  keywords = String[]
  for el in EzXML.eachelement(r)
    n = EzXML.nodename(el)
    if n == "FullModel" || n == "DifferenceModel"
      header = Symbol(n)
      for p in EzXML.eachelement(el)
        EzXML.nodename(p) == "Model.profile" && push!(profile_uris, EzXML.nodecontent(p))
      end
      break
    elseif n == "Dataset"
      # CGMES 3.0 / NC document header (dcat:Dataset) instead of md:FullModel —
      # ReliCapGrid boundary files ship like this. The profile is declared as a
      # short code in dcat:keyword ("BD", "EQ", …), there is no Model.profile.
      # Without this branch such a file gets an empty profile set, is parsed
      # anyway, but never classified — a boundary file then never reaches
      # store.boundary and the assembled-border equivalent-injection rule in
      # the mapping cannot fire.
      header = :Dataset
      for p in EzXML.eachelement(el)
        EzXML.nodename(p) == "keyword" && push!(keywords, EzXML.nodecontent(p))
      end
      break
    end
  end
  profiles = Set{Symbol}(profileTagFromURI(u) for u in profile_uris)
  for kw in keywords
    tag = profileTagFromKeyword(kw)
    tag == :UNKNOWN || push!(profiles, tag)
  end

  if header == :DifferenceModel
    return CGMESFileInfo(file.name, header, profile_uris, profiles, version, 0, true, "difference models are not supported")
  end
  if !isempty(profiles) && isempty(intersect(profiles, profile_filter))
    return CGMESFileInfo(file.name, header, profile_uris, profiles, version, 0, true, "profile not in import set")
  end

  count = 0
  for el in EzXML.eachelement(r)
    cls = Symbol(EzXML.nodename(el))
    (cls == :FullModel || cls == :DifferenceModel) && continue

    mrid = ""
    is_about = false
    for a in EzXML.eachattribute(el)
      an = EzXML.nodename(a)
      if an == "ID"
        mrid = EzXML.nodecontent(a)
      elseif an == "about"
        mrid = _strip_mrid(EzXML.nodecontent(a))
        is_about = true
      end
    end
    isempty(mrid) && continue

    obj = get(objects, mrid, nothing)
    if obj === nothing
      obj = CIMObject(mrid, cls, Dict{Symbol,String}(), Dict{Symbol,String}(), file.name)
      objects[mrid] = obj
      push!(get!(Vector{String}, byclass, cls), mrid)
    elseif !is_about && obj.class != cls
      # duplicate rdf:ID with a different class — keep first, count as object anyway
      @debug "CGMES reader: duplicate rdf:ID $(mrid) ($(obj.class) vs $(cls))"
    end
    count += 1

    for p in EzXML.eachelement(el)
      pname = EzXML.nodename(p)
      res = nothing
      for a in EzXML.eachattribute(p)
        EzXML.nodename(a) == "resource" && (res = EzXML.nodecontent(a))
      end
      if res === nothing
        _property_key!(obj.attrs, pname, EzXML.nodecontent(p))
      elseif startswith(res, "#")
        _property_key!(obj.refs, pname, _strip_mrid(res))
      elseif occursin('#', res)
        # enum value: keep the fragment tail (e.g. WindingConnection.D → "D"
        # stays qualified as "WindingConnection.D" for unambiguous matching)
        _property_key!(obj.attrs, pname, String(last(split(res, '#'))))
      else
        # absolute URI reference (e.g. boundary mRIDs in urn:uuid form)
        _property_key!(obj.refs, pname, _strip_mrid(res))
      end
    end
  end
  return CGMESFileInfo(file.name, header, profile_uris, profiles, version, count, false, "")
end

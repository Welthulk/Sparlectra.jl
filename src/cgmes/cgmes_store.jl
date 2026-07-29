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

# file: src/cgmes/cgmes_store.jl
# purpose: layer 3 of the CGMES importer — merged profile store with by-class
# index and deliberately dumb typed accessors. No query language. The API is
# shaped so that a later column-store swap stays invisible to layers 4/5.

"""
Merged CGMES data set: all profiles of one delivery in a single object store.
`boundary` holds the mRIDs of objects defined in boundary-profile files
(EQ_BD/TP_BD) — needed to recognize X-nodes in assembled models.
"""
struct CGMESStore
  objects::Dict{String,CIMObject}
  byclass::Dict{Symbol,Vector{String}}
  files::Vector{CGMESFileInfo}
  version::String                # "2.4.15" | "3.0" | "" (mixed/unknown)
  boundary::Set{String}
end

"""
    loadCGMES(path; profile_filter=IMPORT_PROFILE_TAGS) -> CGMESStore

Load a CGMES delivery (folder / ZIP / vector of paths, see
`collectCGMESFiles`) into a merged store. EQ is read before the overlay
profiles so that `rdf:about` updates hit existing objects; DifferenceModel
files and out-of-filter profiles are recorded as skipped.
"""
function loadCGMES(path; profile_filter::Set{Symbol} = IMPORT_PROFILE_TAGS)::CGMESStore
  files = collectCGMESFiles(path)
  objects = Dict{String,CIMObject}()
  byclass = Dict{Symbol,Vector{String}}()
  infos = Vector{CGMESFileInfo}(undef, length(files))

  # Order: EQ-family first, then TP, then SSH/SV overlays. Classification
  # needs the header, so read it cheaply once per file up front.
  order = collect(eachindex(files))
  rank = Dict(:EQ => 0, :EQ_OP => 0, :EQ_SC => 0, :EQ_BD => 1, :EQ_BD_OP => 1, :TP_BD => 2, :TP => 2, :SSH => 3, :SV => 4)
  headers = [readCGMESHeader(f) for f in files]
  sort!(order; by = i -> minimum(get(rank, t, 9) for t in headers[i].profiles; init = 9))

  boundary = Set{String}()
  for i in order
    known = Set(keys(objects))
    infos[i] = readCGMESFile!(objects, byclass, files[i]; profile_filter = profile_filter)
    if !isempty(intersect(infos[i].profiles, (:EQ_BD, :EQ_BD_OP, :TP_BD)))
      for m in keys(objects)
        m in known || push!(boundary, m)
      end
    end
  end

  versions = unique(filter(!isempty, [fi.version for fi in infos]))
  version = length(versions) == 1 ? versions[1] : ""
  return CGMESStore(objects, byclass, infos, version, boundary)
end

"""Cheap header-only pass (no object parsing) used for read ordering."""
function readCGMESHeader(file::CGMESFile)::CGMESFileInfo
  empty_objects = Dict{String,CIMObject}()
  empty_byclass = Dict{Symbol,Vector{String}}()
  return readCGMESFile!(empty_objects, empty_byclass, file; profile_filter = Set{Symbol}())
end

# --- accessors (kept deliberately dumb) -------------------------------------

"""All objects of `class`, in file order."""
objectsOf(store::CGMESStore, class::Symbol) = (store.objects[m] for m in get(store.byclass, class, String[]))

"""Number of objects of `class`."""
countOf(store::CGMESStore, class::Symbol)::Int = length(get(store.byclass, class, String[]))

"""Literal attribute as `Float64`, or `default` when absent/unparsable."""
function num(obj::CIMObject, key::Symbol, default::Union{Nothing,Float64} = nothing)
  v = get(obj.attrs, key, nothing)
  v === nothing && return default
  parsed = tryparse(Float64, v)
  return parsed === nothing ? default : parsed
end

"""Literal attribute as `String`, or `default` when absent."""
str(obj::CIMObject, key::Symbol, default::Union{Nothing,String} = nothing) = get(obj.attrs, key, default)

"""Literal attribute as `Bool` (`"true"`/`"false"`), or `default` when absent."""
function boolval(obj::CIMObject, key::Symbol, default::Union{Nothing,Bool} = nothing)
  v = get(obj.attrs, key, nothing)
  v === nothing && return default
  return lowercase(v) == "true"
end

"""Enum attribute value (stored as `Kind.value` fragment), or `default`."""
enumval(obj::CIMObject, key::Symbol, default::Union{Nothing,String} = nothing) = get(obj.attrs, key, default)

"""Follow the reference `key` of `obj`; returns the target `CIMObject` or `nothing`."""
function ref(store::CGMESStore, obj::CIMObject, key::Symbol)::Union{Nothing,CIMObject}
  target = get(obj.refs, key, nothing)
  target === nothing && return nothing
  return get(store.objects, target, nothing)
end

"""
    unresolvedReferences(store) -> Vector{@NamedTuple{mrid, class, key, target}}

All `rdf:resource` targets that do not exist in the store — the concept §7.4
"boundary set missing" detection signal.
"""
function unresolvedReferences(store::CGMESStore)
  out = NamedTuple{(:mrid, :class, :key, :target),Tuple{String,Symbol,Symbol,String}}[]
  for (mrid, obj) in store.objects
    for (k, target) in obj.refs
      haskey(store.objects, target) || push!(out, (mrid = mrid, class = obj.class, key = k, target = target))
    end
  end
  return out
end

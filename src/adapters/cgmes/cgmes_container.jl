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

# file: src/adapters/cgmes/cgmes_container.jl
# purpose: layer 1 of the CGMES importer, turn a folder, a ZIP, several ZIPs
# or ZIP-in-ZIP containers into named XML byte streams. No XML parsing here;
# profile classification happens in the reader (from the md:Model header,
# never from filenames).

"""
One XML payload from a CGMES delivery: `name` is the entry path inside its
`source` container (or the file path for plain folders), `content` the raw
XML. Stored as `String` (not `Vector{UInt8}`) on purpose: the payload is
parsed more than once (header pass + full pass), and `String(bytes)` would
steal the byte buffer on first use.
"""
struct CGMESFile
  name::String
  source::String
  content::String
end

const _MAX_ZIP_NESTING = 4

function _collect_from_zip!(out::Vector{CGMESFile}, bytes::Vector{UInt8}, source::String, depth::Int)
  depth > _MAX_ZIP_NESTING && error("CGMES container: ZIP nesting deeper than $(_MAX_ZIP_NESTING) in $(source)")
  reader = ZipArchives.ZipReader(bytes)
  for entry in ZipArchives.zip_names(reader)
    lname = lowercase(entry)
    if endswith(lname, ".xml")
      push!(out, CGMESFile(entry, source, String(ZipArchives.zip_readentry(reader, entry))))
    elseif endswith(lname, ".zip")
      _collect_from_zip!(out, ZipArchives.zip_readentry(reader, entry), source * "!" * entry, depth + 1)
    end
  end
  return out
end

"""
    collectCGMESFiles(path) -> Vector{CGMESFile}
    collectCGMESFiles(paths::AbstractVector) -> Vector{CGMESFile}

Collect all XML payloads from a CGMES delivery. Accepts a folder (searched
recursively), a `.zip` (nested ZIPs are opened in-memory, up to depth
$( _MAX_ZIP_NESTING)), a single `.xml`, or a vector of any of these
(e.g. base case + boundary ZIP).
"""
function collectCGMESFiles(path::AbstractString)::Vector{CGMESFile}
  out = CGMESFile[]
  if isdir(path)
    for (root, _, files) in walkdir(path)
      for f in sort(files)
        p = joinpath(root, f)
        lf = lowercase(f)
        if endswith(lf, ".xml")
          push!(out, CGMESFile(relpath(p, path), path, read(p, String)))
        elseif endswith(lf, ".zip")
          _collect_from_zip!(out, read(p), p, 1)
        end
      end
    end
  elseif isfile(path)
    lf = lowercase(path)
    if endswith(lf, ".zip")
      _collect_from_zip!(out, read(path), path, 1)
    elseif endswith(lf, ".xml")
      push!(out, CGMESFile(basename(path), dirname(path), read(path, String)))
    else
      error("CGMES container: unsupported file type: $path (expected folder, .zip or .xml)")
    end
  else
    error("CGMES container: path not found: $path")
  end
  isempty(out) && error("CGMES container: no XML files found in $path")
  return out
end

collectCGMESFiles(paths::AbstractVector)::Vector{CGMESFile} = reduce(vcat, (collectCGMESFiles(String(p)) for p in paths); init = CGMESFile[])

# --- model identity without a full read ---------------------------------------

# Bytes read from the start of each XML file to find its model header. The
# md:FullModel block is the first element of a CGMES file and a few KB long;
# 256 KiB leaves room for long profile lists and comments in front of it.
const _CGMES_HEAD_BYTES = 262_144

"""
    cgmesModelKey(id) -> String

Canonical form of a CGMES model id for comparisons: without `urn:uuid:`, a
leading `#` or `_`, surrounding blanks, in lower case. The measurement table
writes the id as found in the EQ header; both sides go through this before
they are compared.
"""
function cgmesModelKey(id::AbstractString)::String
  k = strip(String(id))
  k = startswith(k, "#") ? k[2:end] : k
  k = startswith(lowercase(k), "urn:uuid:") ? k[10:end] : k
  k = startswith(k, "_") ? k[2:end] : k
  return lowercase(String(k))
end

# model id and profile tags of one XML payload, from its first bytes only;
# nothing when the head carries no md:FullModel or dcat:Dataset element
function _cgmes_head_model(bytes::AbstractVector{UInt8})
  # cut at the last '>' so a multi-byte character split by the byte limit
  # does not end up in the text the regular expressions run on
  stop = findlast(==(UInt8('>')), bytes)
  stop === nothing && return nothing
  head = String(bytes[1:stop])
  m = match(r"<(?:[\w\-]+:)?(FullModel|Dataset)\b[^>]*?\brdf:about\s*=\s*\"([^\"]+)\"", head)
  m === nothing && return nothing
  close_tag = findnext(Regex("</(?:[\\w\\-]+:)?" * m.captures[1] * "\\s*>"), head, m.offset)
  block = close_tag === nothing ? head[m.offset:end] : head[m.offset:last(close_tag)]
  profiles = Set{Symbol}(profileTagFromURI(u.captures[1]) for u in eachmatch(r"Model\.profile\s*>\s*([^<]*?)\s*<", block))
  for k in eachmatch(r"keyword\s*>\s*([^<]*?)\s*<", block)
    tag = profileTagFromKeyword(k.captures[1])
    tag == :UNKNOWN || push!(profiles, tag)
  end
  return (id = String(m.captures[2]), profiles = profiles)
end

function _cgmes_head_models_zip!(out::Vector, bytes::Vector{UInt8}, source::String, depth::Int)
  depth > _MAX_ZIP_NESTING && error("CGMES container: ZIP nesting deeper than $(_MAX_ZIP_NESTING) in $(source)")
  reader = ZipArchives.ZipReader(bytes)
  for entry in ZipArchives.zip_names(reader)
    lname = lowercase(entry)
    if endswith(lname, ".xml")
      head = ZipArchives.zip_openentry(io -> read(io, _CGMES_HEAD_BYTES), reader, entry)
      h = _cgmes_head_model(head)
      h === nothing || push!(out, h)
    elseif endswith(lname, ".zip")
      _cgmes_head_models_zip!(out, ZipArchives.zip_readentry(reader, entry), source * "!" * entry, depth + 1)
    end
  end
  return out
end

"""
    cgmesEquipmentModelIds(path) -> Vector{String}

The model ids of the EQ profile files of a CGMES delivery (folder, ZIP with
nested ZIPs, or one XML file): the `rdf:about` of each md:FullModel (or
dcat:Dataset) header that declares the equipment profile, in the form of
[`cgmesModelKey`](@ref), sorted and unique. Boundary EQ files do not count.
Only the first bytes of each file are read and no object is parsed, so the
call is cheap enough for a page that lists measurement tables per case.

Empty when no file declares the EQ profile. A missing path or an unsupported
file type raises the same error as [`collectCGMESFiles`](@ref); a damaged ZIP
raises the error of the ZIP reader.
"""
function cgmesEquipmentModelIds(path::AbstractString)::Vector{String}
  heads = NamedTuple{(:id, :profiles),Tuple{String,Set{Symbol}}}[]
  if isdir(path)
    for (root, _, files) in walkdir(path)
      for f in sort(files)
        p = joinpath(root, f)
        lf = lowercase(f)
        if endswith(lf, ".xml")
          h = _cgmes_head_model(open(io -> read(io, _CGMES_HEAD_BYTES), p))
          h === nothing || push!(heads, h)
        elseif endswith(lf, ".zip")
          _cgmes_head_models_zip!(heads, read(p), p, 1)
        end
      end
    end
  elseif isfile(path)
    lf = lowercase(path)
    if endswith(lf, ".zip")
      _cgmes_head_models_zip!(heads, read(path), String(path), 1)
    elseif endswith(lf, ".xml")
      h = _cgmes_head_model(open(io -> read(io, _CGMES_HEAD_BYTES), path))
      h === nothing || push!(heads, h)
    else
      error("CGMES container: unsupported file type: $path (expected folder, .zip or .xml)")
    end
  else
    error("CGMES container: path not found: $path")
  end
  return sort!(unique([cgmesModelKey(h.id) for h in heads if :EQ in h.profiles]))
end

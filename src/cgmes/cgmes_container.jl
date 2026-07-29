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

# file: src/cgmes/cgmes_container.jl
# purpose: layer 1 of the CGMES importer — turn a folder, a ZIP, several ZIPs
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

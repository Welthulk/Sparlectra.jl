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

# file: examples/experimental/cgmes_fetch_testsets.jl
# purpose: fetch the ENTSO-E/ReliCapGrid CGMES test sets into the local cache (see docs/src/cgmes_import.md).
# Fetch the ENTSO-E CGMES conformity test configurations on demand, unpack the
# ZIP-in-ZIP layout recursively, and print an inventory of grids, profiles and
# CGMES versions found. Analysis tooling only — not part of the package API and
# free of non-stdlib dependencies (Downloads, SHA, p7zip_jll are stdlibs).
#
# The archive lands in the repository's data directory but is never
# committed (`data/CGMES/` is gitignored):
#   default   <repo>/data/CGMES
#   override  ENV["SPARLECTRA_CGMES_CACHE"]
#
# Usage:
#   julia examples/experimental/cgmes_fetch_testsets.jl
#   julia examples/experimental/cgmes_fetch_testsets.jl --overwrite
#
# If the download fails (URL moved), download the package manually from the
# ENTSO-E "CIM Conformity and Interoperability" page and place the ZIP at
# <cache>/TestConfigurations_packageCASv2.0.zip, then re-run this script.

using Downloads
using SHA
using Printf
using p7zip_jll

# Primary URL redirects to eepublicdownloads.entsoe.eu (verified 2026-07-27,
# ~22 MB). The legacy docstore.entsoe.eu host no longer resolves.
const CGMES_TESTSET_URLS = [
  "https://www.entsoe.eu/Documents/CIM_documents/Grid_Model_CIM/TestConfigurations_packageCASv2.0.zip",
  "https://eepublicdownloads.entsoe.eu/clean-documents/CIM_documents/Grid_Model_CIM/TestConfigurations_packageCASv2.0.zip",
]

cgmes_cache_dir() = get(ENV, "SPARLECTRA_CGMES_CACHE", joinpath(dirname(dirname(@__DIR__)), "data", "CGMES"))

function fetch_package(; cache::AbstractString = cgmes_cache_dir(), overwrite::Bool = false)::String
  mkpath(cache)
  dst = joinpath(cache, split(CGMES_TESTSET_URLS[1], "/")[end])
  if isfile(dst) && !overwrite
    @printf("exists: %s (%.1f MB)\n", dst, filesize(dst) / 1e6)
    return dst
  end
  for (i, url) in enumerate(CGMES_TESTSET_URLS)
    try
      @printf("downloading (%d/%d): %s\n", i, length(CGMES_TESTSET_URLS), url)
      Downloads.download(url, dst)
      @printf("downloaded %s %s\n", basename(dst), bytes2hex(sha256(read(dst))))
      return dst
    catch err
      @printf("  failed: %s\n", sprint(showerror, err))
    end
  end
  error("""
        Could not download the ENTSO-E test configurations from any known URL.
        Download manually from the ENTSO-E "CIM Conformity and Interoperability" page
        (https://www.entsoe.eu/data/cim/cim-conformity-and-interoperability/) and place
        the file at: $dst
        """)
end

# Recursively extract: top-level package plus any nested ZIPs (the ENTSO-E
# package wraps each test configuration in its own ZIP, sometimes two deep).
function extract_recursive(zip::AbstractString, outdir::AbstractString; depth::Int = 0, maxdepth::Int = 4)
  depth > maxdepth && error("ZIP nesting deeper than $maxdepth at $zip — aborting (unexpected layout)")
  mkpath(outdir)
  run(pipeline(`$(p7zip_jll.p7zip()) x -y -o$outdir $zip`; stdout = devnull))
  for (root, _, files) in walkdir(outdir)
    for f in files
      if endswith(lowercase(f), ".zip")
        inner = joinpath(root, f)
        innerdir = joinpath(root, first(splitext(f)))
        if !isdir(innerdir)
          extract_recursive(inner, innerdir; depth = depth + 1, maxdepth = maxdepth)
        end
      end
    end
  end
end

# Filename-token based classification — good enough for an inventory. The
# authoritative header-based profile detection is work package A2.
function classify_profile(fname::AbstractString)::String
  n = uppercase(fname)
  for token in ("EQ_BD", "TP_BD", "EQBD", "TPBD")
    occursin(token, n) && return "BD"
  end
  for token in ("_EQ", "_TP", "_SSH", "_SV", "_DL", "_GL", "_DY")
    occursin(token, n) && return strip(token, '_') |> String
  end
  return "?"
end

function inventory(extracted::AbstractString)::Vector{NamedTuple}
  rows = NamedTuple[]
  for (root, _, files) in walkdir(extracted)
    for f in files
      endswith(lowercase(f), ".xml") || continue
      rel = relpath(joinpath(root, f), extracted)
      grid = first(splitpath(rel))
      push!(rows, (grid = grid, profile = classify_profile(f), file = rel, kB = round(Int, filesize(joinpath(root, f)) / 1024)))
    end
  end
  return rows
end

function print_inventory(rows::Vector{NamedTuple}, outfile::AbstractString)
  open(outfile, "w") do io
    for sink in (stdout, io)
      println(sink, "# CGMES test-set inventory (", length(rows), " XML files)\n")
      grids = sort(unique(r.grid for r in rows))
      println(sink, "| Top-level entry | XML files | Profiles found |")
      println(sink, "|---|---:|---|")
      for g in grids
        sub = filter(r -> r.grid == g, rows)
        profs = sort(unique(r.profile for r in sub))
        println(sink, "| ", g, " | ", length(sub), " | ", join(profs, ", "), " |")
      end
      println(sink)
    end
    println(io, "\n## All files\n\n| File | Profile | kB |\n|---|---|---:|")
    for r in sort(rows; by = r -> r.file)
      println(io, "| ", r.file, " | ", r.profile, " | ", r.kB, " |")
    end
  end
  println("full inventory written to: ", outfile)
end

function main(args::Vector{String} = ARGS)
  overwrite = "--overwrite" in args
  cache = cgmes_cache_dir()
  pkg = fetch_package(; cache = cache, overwrite = overwrite)
  extracted = joinpath(cache, "extracted")
  if !isdir(extracted) || overwrite
    isdir(extracted) && rm(extracted; recursive = true)
    println("extracting (recursive) ...")
    extract_recursive(pkg, extracted)
  end
  rows = inventory(extracted)
  print_inventory(rows, joinpath(cache, "inventory.md"))
  println("cache dir: ", cache)
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end

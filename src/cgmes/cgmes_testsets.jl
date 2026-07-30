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

# file: src/cgmes/cgmes_testsets.jl
# purpose: on-demand access to the ENTSO-E CGMES conformity test
# configurations. Downloads the official package once into the gitignored
# `data/CGMES` cache, extracts it (pure Julia, nested ZIPs included), and
# packs a named test set — base case plus boundary where needed — into a
# single ZIP that any CGMES entry point (importCGMES, the API, the Web UI
# upload flow) accepts.

# Primary URL redirects to eepublicdownloads.entsoe.eu (verified 2026-07-27);
# the legacy docstore.entsoe.eu host no longer resolves.
const CGMES_TESTSET_URLS = [
  "https://www.entsoe.eu/Documents/CIM_documents/Grid_Model_CIM/TestConfigurations_packageCASv2.0.zip",
  "https://eepublicdownloads.entsoe.eu/clean-documents/CIM_documents/Grid_Model_CIM/TestConfigurations_packageCASv2.0.zip",
]

"""Cache directory for the ENTSO-E test configurations (`data/CGMES`, override via `SPARLECTRA_CGMES_CACHE`)."""
cgmesTestSetCacheDir() = get(ENV, "SPARLECTRA_CGMES_CACHE", joinpath(dirname(dirname(@__DIR__)), "data", "CGMES"))

const _MG = joinpath("MicroGrid", "BaseCase_BC")
const _SGB = joinpath("SmallGrid", "BusBranch")
const _SGN = joinpath("SmallGrid", "NodeBreaker")

"""
Known test-set aliases for `fetchCGMESTestSet` (and the Web UI `cgmes:` case
entry): alias → subdirectories of the extracted package that form the
delivery (base case plus boundary set where one exists).
"""
const CGMES_TESTSET_ALIASES = Dict{String,Vector{String}}(
  # convenience shorthand: "microgrid" is the BE base configuration — the
  # reference demonstrator used throughout the docs and tests
  "microgrid" => [joinpath(_MG, "CGMES_v2.4.15_MicroGridTestConfiguration_BC_BE_v2"), joinpath(_MG, "CGMES_v2.4.15_MicroGridTestConfiguration_BD_v2")],
  "microgrid_be" => [joinpath(_MG, "CGMES_v2.4.15_MicroGridTestConfiguration_BC_BE_v2"), joinpath(_MG, "CGMES_v2.4.15_MicroGridTestConfiguration_BD_v2")],
  "microgrid_nl" => [joinpath(_MG, "CGMES_v2.4.15_MicroGridTestConfiguration_BC_NL_v2"), joinpath(_MG, "CGMES_v2.4.15_MicroGridTestConfiguration_BD_v2")],
  "microgrid_assembled" => [joinpath(_MG, "CGMES_v2.4.15_MicroGridTestConfiguration_BC_Assembled_v2"), joinpath(_MG, "CGMES_v2.4.15_MicroGridTestConfiguration_BD_v2")],
  "smallgrid" => [joinpath(_SGB, "CGMES_v2.4.15_SmallGridTestConfiguration_BaseCase_Complete_v3.0.0"), joinpath(_SGB, "CGMES_v2.4.15_SmallGridTestConfiguration_Boundary_v3.0.0")],
  "smallgrid_nb" => [joinpath(_SGN, "CGMES_v2.4.15_SmallGridTestConfiguration_BaseCase_Complete_v3.0.0"), joinpath(_SGN, "CGMES_v2.4.15_SmallGridTestConfiguration_Boundary_v3.0.0")],
  "fullgrid" => [joinpath("FullGrid", "CGMES_v2.4.15_FullGridTestConfiguration_BB_BE_v2"), joinpath("FullGrid", "CGMES_v2.4.15_FullGridTestConfiguration_BD_v1")],
  "fullgrid_nb" => [joinpath("FullGrid", "CGMES_v2.4.15_FullGridTestConfiguration_NB_BE_v4"), joinpath("FullGrid", "CGMES_v2.4.15_FullGridTestConfiguration_BD_v1")],
  "realgrid" => [joinpath("RealGrid", "CGMES_v2.4.15_RealGridTestConfiguration_v2")],
  # same directories the CGMES example suite validates (bus-branch base case
  # vs. node-breaker complete case, each with the MiniGrid boundary set)
  "minigrid" => [joinpath("MiniGrid", "BusBranch", "CGMES_v2.4.15_MiniGridTestConfiguration_BaseCase_v3"), joinpath("MiniGrid", "BusBranch", "CGMES_v2.4.15_MiniGridTestConfiguration_Boundary_v3")],
  "minigrid_nb" => [joinpath("MiniGrid", "NodeBreaker", "CGMES_v2.4.15_MiniGridTestConfiguration_BaseCase_Complete_v3"), joinpath("MiniGrid", "NodeBreaker", "CGMES_v2.4.15_MiniGridTestConfiguration_Boundary_v3")],
)

# --- second source: ENTSO-E ReliCapGrid (GitHub, CGMES 3.0) ----------------
#
# A different publication route than the conformity package: individual
# CIM/XML files in a repository, one folder per synthetic model, with the
# boundary set split per border. Files are fetched over raw.githubusercontent
# and packed into one ZIP, exactly like the conformity aliases.

const RELICAPGRID_REPO = "entsoe/relicapgrid"
const RELICAPGRID_BRANCH = "cgmes-3.0_ncp-2.5_tc-2.0"
const _RCG_GRID = "Instance/%s/Grid/cimxml"
const _RCG_BD = "Instance/boundaryData/Grid/cimxml"
const _RCG_COMMON = "Instance/commonData/Grid/cimxml/Grid_CommonData_CGM-CD.xml"

"""
ReliCapGrid aliases: alias → (model folder, grid file names, boundary file
names). The models are CGMES **3.0**; their boundary set is published per
border, so a model may need more than one boundary file.
"""
const RELICAPGRID_ALIASES = Dict{String,NamedTuple{(:model, :grid, :boundary),Tuple{String,Vector{String},Vector{String}}}}(
  "svedala" => (
    model = "Svedala",
    grid = ["20220615T2230Z__Svedala_EQ_1.xml", "20220615T2230Z_2D_Svedala_SSH_1.xml", "20220615T2230Z_2D_Svedala_TP_1.xml", "20220615T2230Z_2D_Svedala_SV_1.xml"],
    boundary = ["Boundary_Border-Svedala-Belgovia.xml", "Boundary_Border-Svedala-Espheim.xml"],
  ),
  "espheim" => (
    model = "Espheim",
    grid = ["20220615T2230Z__Espheim_EQ_1.xml", "20220615T2230Z_2D_Espheim_SSH_1.xml", "20220615T2230Z_2D_Espheim_TP_1.xml", "20220615T2230Z_2D_Espheim_SV_1.xml"],
    boundary = ["Boundary_Border-Svedala-Espheim.xml", "Boundary_Border-Espheim-Portheim.xml"],
  ),
  "belgovia" => (
    model = "Belgovia",
    grid = ["20220615T2230Z__Belgovia_EQ_1.xml", "20220615T2230Z_2D_Belgovia_SSH_1.xml", "20220615T2230Z_2D_Belgovia_TP_1.xml", "20220615T2230Z_2D_Belgovia_SV_1.xml"],
    boundary = ["Boundary_Border-Svedala-Belgovia.xml", "Boundary_Border-Galia-Belgovia.xml"],
  ),
  "galia" => (
    model = "Galia",
    grid = ["20220615T2230Z__Galia_EQ_1.xml", "20220615T2230Z_2D_Galia_SSH_1.xml", "20220615T2230Z_2D_Galia_TP_1.xml", "20220615T2230Z_2D_Galia_SV_1.xml"],
    boundary = ["Boundary_Border-Galia-Belgovia.xml", "Boundary_Border-Galia-Britheim.xml", "Boundary_Border-Galia-Nordheim.xml"],
  ),
  "nordheim" => (
    model = "Nordheim",
    grid = ["20220615T2230Z__Nordheim_EQ_1.xml", "20220615T2230Z_2D_Nordheim_SSH_1.xml", "20220615T2230Z_2D_Nordheim_TP_1.xml", "20220615T2230Z_2D_Nordheim_SV_1.xml"],
    boundary = ["Boundary_Border-Galia-Nordheim.xml"],
  ),
  "britheim" => (
    model = "Britheim",
    grid = ["20220615T2230Z__Britheim_EQ_1.xml", "20220615T2230Z_2D_Britheim_SSH_1.xml", "20220615T2230Z_2D_Britheim_TP_1.xml", "20220615T2230Z_2D_Britheim_SV_1.xml"],
    boundary = ["Boundary_Border-Galia-Britheim.xml"],
  ),
  # Portheim is Espheim's other neighbour. Note its files carry a different
  # timestamp (2024-12-23) than the rest of the family (2022-06-15) — it was
  # published later, so a combined model mixes two snapshot times.
  "portheim" => (
    model = "Portheim",
    grid = ["20241223T0642Z_2D_Portheim_EQ_1.xml", "20241223T0642Z_2D_Portheim_SSH_1.xml", "20241223T0642Z_2D_Portheim_TP_1.xml", "20241223T0642Z_2D_Portheim_SV_1.xml"],
    boundary = ["Boundary_Border-Espheim-Portheim.xml"],
  ),
)

"""
Combined ReliCapGrid deliveries: alias → member aliases of
[`RELICAPGRID_ALIASES`](@ref), fetched and packed into ONE delivery.

A single ReliCapGrid model is one area of a multi-area system. Imported alone,
the nodes on its borders hang free — there is nothing behind them — so the
power flow has no chance of producing a meaningful result no matter how good
the solver is. Combining the areas across a shared boundary file closes those
borders.

The family has exactly six borders:

    Svedala — Espheim — Portheim
       |
    Belgovia — Galia — Britheim
                  |
               Nordheim

`relicapgrid_cgm` therefore takes all seven models and is the only fully
closed combination. `svedala_neighbours` is the cheap variant that closes
Svedala's own two borders but leaves Espheim—Portheim and Belgovia—Galia open.

(`Jotunheim` exists in the repository but ships only TP and SV — no EQ/SSH and
no border — so it cannot be imported and is not offered.)
"""
const RELICAPGRID_COMBINED = Dict{String,Vector{String}}(
  "relicapgrid_cgm" => ["svedala", "espheim", "portheim", "belgovia", "galia", "britheim", "nordheim"],
  "svedala_neighbours" => ["svedala", "espheim", "belgovia"],
)

_rcgRawURL(path::AbstractString) = "https://raw.githubusercontent.com/" * RELICAPGRID_REPO * "/" * RELICAPGRID_BRANCH * "/" * path

"""All member aliases of `key`: itself for a single model, the member list for a combined one."""
_rcgMembers(key::AbstractString)::Vector{String} = get(RELICAPGRID_COMBINED, key, [String(key)])

"""
Files one ReliCapGrid member contributes, as `(zip_entry, repo_path, cache_dir)`.

Grid files are cached per model folder (mirroring the repository layout);
boundary and commonData files are shared between models and therefore cached
once in a common folder, so a combined fetch does not download them per member.
"""
function _rcgMemberFiles(member::AbstractString)
  spec = RELICAPGRID_ALIASES[member]
  root = joinpath(cgmesTestSetCacheDir(), "relicapgrid", RELICAPGRID_BRANCH)
  model_cache = joinpath(root, spec.model)
  shared_cache = joinpath(root, "_shared")
  out = Tuple{String,String,String}[]
  for f in spec.grid
    push!(out, (f, replace(_RCG_GRID, "%s" => spec.model) * "/" * f, model_cache))
  end
  for f in spec.boundary
    push!(out, (f, _RCG_BD * "/" * f, shared_cache))
  end
  # shared common data (base voltages, geographical regions, …) that every
  # ReliCapGrid model references
  push!(out, (basename(_RCG_COMMON), _RCG_COMMON, shared_cache))
  return out
end

"""
    fetchReliCapGridSet(alias; outdir) -> String

Fetch a ReliCapGrid delivery (CGMES 3.0) — grid profiles plus the boundary
files of its borders — from GitHub and pack it into `<outdir>/cgmes_<alias>.zip`.
Downloaded files are cached under the CGMES cache directory, so repeated calls
do not hit the network.

`alias` is either a single model (`"svedala"`) or a combined delivery
(`"relicapgrid_cgm"`, see [`RELICAPGRID_COMBINED`](@ref)); a combined delivery
packs all its members plus the shared boundary files into one ZIP.
"""
function fetchReliCapGridSet(alias::AbstractString; outdir::AbstractString)::String
  key = lowercase(strip(String(alias)))
  known = union(keys(RELICAPGRID_ALIASES), keys(RELICAPGRID_COMBINED))
  key in known || error("unknown ReliCapGrid set '$(alias)'. Available: $(join(sort(collect(known)), ", "))")
  mkpath(outdir)
  dest = joinpath(outdir, "cgmes_" * key * ".zip")
  isfile(dest) && return dest

  # Collect over all members and de-duplicate: neighbouring areas share the
  # boundary file of their common border, and every model references the same
  # commonData file. Writing a ZIP entry twice would produce a duplicate name.
  wanted = Tuple{String,String,String}[]
  seen = Set{String}()
  for member in _rcgMembers(key)
    for entry in _rcgMemberFiles(member)
      entry[1] in seen && continue
      push!(seen, entry[1])
      push!(wanted, entry)
    end
  end

  for (name, path, cachedir) in wanted
    mkpath(cachedir)
    local_file = joinpath(cachedir, name)
    isfile(local_file) && continue
    try
      Downloads.download(_rcgRawURL(path), local_file)
    catch err
      isfile(local_file) && rm(local_file; force = true)
      error("""
        Could not download '$(name)' from the ReliCapGrid repository ($(sprint(showerror, err))).
        Source: https://github.com/$(RELICAPGRID_REPO)/tree/$(RELICAPGRID_BRANCH)
        Place the files manually in: $(cachedir)
        """)
    end
  end
  tmp = dest * ".tmp"
  open(tmp, "w") do io
    ZipArchives.ZipWriter(io) do w
      for (name, _, cachedir) in wanted
        ZipArchives.zip_newfile(w, name)
        write(w, read(joinpath(cachedir, name)))
      end
    end
  end
  mv(tmp, dest; force = true)
  return dest
end

"""All known test-set aliases across both sources (conformity package and ReliCapGrid)."""
allCGMESTestSetAliases() = sort(union(collect(keys(CGMES_TESTSET_ALIASES)), collect(keys(RELICAPGRID_ALIASES)), collect(keys(RELICAPGRID_COMBINED))))

# recursive pure-Julia extraction; nested ZIPs land in a folder named after
# the ZIP (matching the layout the analysis phase established)
function _extractZipTree(zipbytes::Vector{UInt8}, outdir::AbstractString; depth::Int = 0)
  depth > 4 && error("CGMES test sets: ZIP nesting deeper than 4 — unexpected layout")
  reader = ZipArchives.ZipReader(zipbytes)
  for entry in ZipArchives.zip_names(reader)
    endswith(entry, "/") && continue
    data = ZipArchives.zip_readentry(reader, entry)
    if endswith(lowercase(entry), ".zip")
      _extractZipTree(data, joinpath(outdir, dirname(entry), first(splitext(basename(entry)))); depth = depth + 1)
    else
      dest = joinpath(outdir, entry)
      mkpath(dirname(dest))
      write(dest, data)
    end
  end
  return nothing
end

"""
    ensureCGMESTestConfigurations(; cache = cgmesTestSetCacheDir()) -> String

Make sure the ENTSO-E test-configuration package is downloaded and extracted
under `cache`; returns the `extracted` directory. Downloads once (~22 MB);
if every URL fails, the error explains where to place the ZIP manually.
"""
function ensureCGMESTestConfigurations(; cache::AbstractString = cgmesTestSetCacheDir())::String
  extracted = joinpath(cache, "extracted")
  isdir(extracted) && return extracted
  mkpath(cache)
  pkg = joinpath(cache, "TestConfigurations_packageCASv2.0.zip")
  if !isfile(pkg)
    lasterr = ""
    ok = false
    for url in CGMES_TESTSET_URLS
      try
        Downloads.download(url, pkg)
        ok = true
        break
      catch err
        lasterr = sprint(showerror, err)
      end
    end
    ok || error("""
      Could not download the ENTSO-E CGMES test configurations ($(lasterr)).
      Download the package manually from the ENTSO-E "CIM Conformity and
      Interoperability" page and place it at: $(pkg)
      """)
  end
  _extractZipTree(read(pkg), extracted)
  return extracted
end

"""
    fetchCGMESTestSet(alias; outdir) -> String

Resolve a test-set alias (see [`CGMES_TESTSET_ALIASES`](@ref), e.g.
`"microgrid_be"`, `"smallgrid"`, `"realgrid"`) and pack the delivery — base
case plus boundary where one exists — into `<outdir>/cgmes_<alias>.zip`.
Returns the ZIP path; the file is reused when it already exists. Downloads
and extracts the ENTSO-E package on first use.
"""
function fetchCGMESTestSet(alias::AbstractString; outdir::AbstractString)::String
  key = lowercase(strip(String(alias)))
  # ReliCapGrid is a second source with its own fetch route (GitHub raw files
  # instead of the conformity package); single models and combined deliveries
  # both dispatch there.
  (haskey(RELICAPGRID_ALIASES, key) || haskey(RELICAPGRID_COMBINED, key)) && return fetchReliCapGridSet(key; outdir = outdir)
  haskey(CGMES_TESTSET_ALIASES, key) || error("unknown CGMES test set '$(alias)'. Available: $(join(allCGMESTestSetAliases(), ", "))")
  mkpath(outdir)
  dest = joinpath(outdir, "cgmes_" * key * ".zip")
  isfile(dest) && return dest
  extracted = ensureCGMESTestConfigurations()
  dirs = [joinpath(extracted, d) for d in CGMES_TESTSET_ALIASES[key]]
  for d in dirs
    isdir(d) || error("CGMES test set '$(key)': expected directory missing after extraction: $(d)")
  end
  tmp = dest * ".tmp"
  open(tmp, "w") do io
    ZipArchives.ZipWriter(io) do w
      # Multi-directory aliases (e.g. microgrid_assembled = BE + NL + BD) ship
      # the identical boundary files in every side directory; ZipWriter rejects
      # duplicate entry names, so each file name is packed exactly once.
      seen = Set{String}()
      for d in dirs, f in sort(readdir(d))
        endswith(lowercase(f), ".xml") || continue
        f in seen && continue
        push!(seen, f)
        ZipArchives.zip_newfile(w, f)
        write(w, read(joinpath(d, f)))
      end
    end
  end
  mv(tmp, dest; force = true)
  return dest
end

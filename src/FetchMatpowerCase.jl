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
#
# file: src/FetchMatpowerCase.jl
#

using Printf

"""
Utilities for downloading MATPOWER-compatible case files on demand.

Note:
- Downloaded files are stored locally and are not part of the Sparlectra
  source distribution.
- Users are responsible for complying with the license terms of the
  respective upstream sources.
"""
module FetchMatpowerCase


export ensure_matpower_case, fetch_matpower_case, emit_julia_case, ensure_casefile, main

using Downloads
using SHA
using Printf
using ..MatpowerIO


basename(path::AbstractString) = splitpath(path)[end]

"""
    fetch_matpower_case(url, outdir; overwrite=false) -> mfile::String

Downloads a MATPOWER *.m file from `url` into `outdir` (unless already present).
Returns the local path to the downloaded file.
"""
function fetch_matpower_case(url::AbstractString, outdir::AbstractString; overwrite::Bool=false)::String
    mkpath(outdir)
    fname = split(url, "/")[end]
    isempty(fname) && error("Could not derive filename from URL: $url")
    dst = joinpath(outdir, fname)

    if isfile(dst) && !overwrite
        return dst
    end

    Downloads.download(url, dst)
    return dst
end

"""
    emit_julia_case(mfile, outdir; legacy_compat=true, overwrite=false) -> jlfile::String

Parses a MATPOWER *.m case and writes a Julia NamedTuple (*.jl) compatible with
MatpowerIO.read_case_julia().
"""
function emit_julia_case(mfile::AbstractString, outdir::AbstractString;
                         legacy_compat::Bool=true, overwrite::Bool=false)::String
    mkpath(outdir)
    mpc = MatpowerIO.read_case(mfile; legacy_compat=legacy_compat)

    jlfile = joinpath(outdir, "$(mpc.name).jl")
    if isfile(jlfile) && !overwrite
        return jlfile
    end

    open(jlfile, "w") do io
        println(io, "# Julia MATPOWER-like case (NamedTuple) for Sparlectra / MatpowerIO.read_case_julia()")
        println(io, "# Auto-generated from: $(basename(mfile))")
        println(io, "")
        println(io, "(")
        println(io, "  name = \"$(mpc.name)\",")
        println(io, "  baseMVA = $(mpc.baseMVA),")
        println(io, "")

        print(io, "  bus = ")
        write_matrix_literal(io, mpc.bus; indent="    ")
        println(io, ",")

        println(io, "")
        print(io, "  gen = ")
        write_matrix_literal(io, mpc.gen; indent="    ")
        println(io, ",")

        println(io, "")
        print(io, "  branch = ")
        write_matrix_literal(io, mpc.branch; indent="    ")
        println(io, ",")

        println(io, "")
        if mpc.gencost === nothing
            println(io, "  gencost = nothing,")
        else
            print(io, "  gencost = ")
            write_matrix_literal(io, mpc.gencost; indent="    ")
            println(io, ",")
        end

        if mpc.bus_name === nothing
            println(io, "  bus_name = nothing,")
        else
            # bus_name is rare; emit as a literal vector of strings
            println(io, "  bus_name = [")
            for s in mpc.bus_name
                println(io, "    \"$(replace(s, "\"" => "\\\""))\",")
            end
            println(io, "  ],")
        end

        println(io, ")")
    end

    return jlfile
end

# Write a Matrix{Float64} as a Julia literal in MATPOWER table style.
function write_matrix_literal(io::IO, M::AbstractMatrix{<:Real}; indent::String="    ")
    println(io, "[")
    maxrows = size(M, 1)
    for i in 1:maxrows
        print(io, indent)
        maxcols = size(M, 2)
        for j in 1:maxcols
            # Use %.12g: compact, stable, enough digits for cases
            if j == 1
                @printf(io, "%.12g", Float64(M[i, j]))
            else
                @printf(io, "  %.12g", Float64(M[i, j]))
            end
        end
        println(io, ";")
    end
    print(io, rstrip(indent))  # align closing bracket
    println(io, "]")
end


"""
    ensure_matpower_case(; url, outdir, to_jl=true, overwrite=false, legacy_compat=true, verbose=true)
        -> (mfile::String, jlfile::Union{Nothing,String})

Ensures the MATPOWER *.m file exists in `outdir` by downloading it if missing.
Optionally emits a Julia NamedTuple case (*.jl).

This is intended to be called from examples and other code.
"""
function ensure_matpower_case(; url::AbstractString,
                               outdir::AbstractString,
                               to_jl::Bool=true,
                               overwrite::Bool=false,
                               legacy_compat::Bool=true,
                               verbose::Bool=true)

    mfile = joinpath(outdir, split(url, "/")[end])
    existed = isfile(mfile)
    if !existed || overwrite
        mfile = fetch_matpower_case(url, outdir; overwrite=overwrite)
        if verbose
            bytes = read(mfile)
            @printf("downloaded %s %s\n", basename(mfile), bytes2hex(sha256(bytes)))
        end
    else
        verbose && @printf("exists: %s\n", basename(mfile))
    end

    jlfile = nothing
    if to_jl
        jlfile = emit_julia_case(mfile, outdir; legacy_compat=legacy_compat, overwrite=overwrite)
        verbose && @printf("ready: %s\n", basename(jlfile))
    end

    return (mfile=mfile, jlfile=jlfile)
end

# ---------------------------
# Optional CLI wrapper
# ---------------------------

function print_help(io::IO=stdout)
    println(io, """
fetch_matpower_case.jl

Download a MATPOWER *.m case file and optionally emit a Julia NamedTuple case.

Options:
  --url <RAW_URL>           Direct URL to a *.m file (raw content).
  --outdir <DIR>            Output directory (default: data/).
  --to-jl                   Also generate a *.jl NamedTuple case.
  --overwrite               Overwrite existing files.
  --legacy-compat           Keep legacy bus sorting (default).
  --no-legacy-compat        Disable legacy bus sorting.
  -h, --help                Show this help.

Note:
  If --url is omitted, this script prints help and exits without error.
""")
end

function parse_args(args::Vector{String})
    opts = Dict{String,Any}(
        "url" => nothing,
        "outdir" => normpath(joinpath(@__DIR__, "..", "data")),
        "to_jl" => false,
        "overwrite" => false,
        "legacy_compat" => true,
    )

    want_help = false
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--url"
            i += 1
            i > length(args) && error("Missing value after --url")
            opts["url"] = args[i]
        elseif a == "--outdir"
            i += 1
            i > length(args) && error("Missing value after --outdir")
            opts["outdir"] = normpath(args[i])
        elseif a == "--to-jl"
            opts["to_jl"] = true
        elseif a == "--overwrite"
            opts["overwrite"] = true
        elseif a == "--legacy-compat"
            opts["legacy_compat"] = true
        elseif a == "--no-legacy-compat"
            opts["legacy_compat"] = false
        elseif a in ("-h", "--help")
            want_help = true
        else
            error("Unknown argument: $a (use --help)")
        end
        i += 1
    end

    return opts, want_help
end

"""
    main(args=ARGS) -> nothing

CLI entry point. If no --url is provided, it prints help and returns without error.
"""
function main(args=ARGS)
    opts, want_help = parse_args(collect(args))
    url = opts["url"]

    if want_help || url === nothing
        print_help()
        return nothing
    end

    ensure_matpower_case(
        url=url,
        outdir=opts["outdir"],
        to_jl=opts["to_jl"],
        overwrite=opts["overwrite"],
        legacy_compat=opts["legacy_compat"],
        verbose=true,
    )

    return nothing
end

"""
    ensure_casefile(casefile; outdir=nothing, overwrite=false, to_jl=true) -> String

Ensure a MATPOWER-compatible case file exists locally.

- If `casefile` is an existing path, it is returned unchanged.
- If `casefile` is a bare filename (e.g. `case14.m`), it is downloaded into `outdir`.
- If `casefile` ends with `.jl` and is missing, the corresponding `.m` is downloaded and `.jl` is generated.

Returns the local path to the requested case file.
"""
function ensure_casefile(casefile::AbstractString;
                         outdir::Union{Nothing,AbstractString}=nothing,
                         overwrite::Bool=false,
                         to_jl::Bool=true)::String

  # 1) If user passed an existing file path, just use it.
  if isfile(casefile)
    return normpath(casefile)
  end

  # 2) Decide output directory.
  # Prefer repo root: <repo>/data/mpower
  if outdir === nothing
    # Find package root from Sparlectra source file location:
    # <pkg>/src/Sparlectra.jl -> <pkg>
    pkgroot = normpath(joinpath(@__DIR__, ".."))
    outdir = normpath(joinpath(pkgroot, "data", "mpower"))
  end
  mkpath(outdir)

  # 3) If looks like a path but does not exist, fail explicitly.
  if occursin(r"[\\/]", casefile)
    error("Case file not found: $casefile")
  end

  lcase = lowercase(casefile)

  if endswith(lcase, ".m")
    url = "https://raw.githubusercontent.com/MATPOWER/matpower/master/data/$(casefile)"
    ensure_matpower_case(url=url, outdir=outdir, to_jl=to_jl, overwrite=overwrite)
    return joinpath(outdir, casefile)
  elseif endswith(lcase, ".jl")
    mname = casefile[1:end-3] * ".m"
    url = "https://raw.githubusercontent.com/MATPOWER/matpower/master/data/$(mname)"
    ensure_matpower_case(url=url, outdir=outdir, to_jl=true, overwrite=overwrite)
    return joinpath(outdir, casefile)
  else
    error("Unsupported casefile extension: $casefile (expected .m or .jl)")
  end
end


end # module


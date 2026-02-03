#!/usr/bin/env julia
#------------------------------------------------------------------------------
# examples/fetch_matpower_case.jl
#
# Download a MATPOWER *.m case file and optionally emit a Julia NamedTuple case
# compatible with MatpowerIO.read_case_julia().
#
# Usage:
#   julia --project=. examples/fetch_matpower_case.jl --url <RAW_URL> [--outdir data] [--to-jl] [--overwrite]
#
# Examples:
#   julia --project=. examples/fetch_matpower_case.jl --url https://raw.githubusercontent.com/MATPOWER/matpower/master/data/case14.m --to-jl
#   julia --project=. examples/fetch_matpower_case.jl --url https://raw.githubusercontent.com/MATPOWER/matpower/master/data/case118.m --outdir data/cases --to-jl
#------------------------------------------------------------------------------

using Downloads
using SHA

# Load your local MatpowerIO module (adjust if your path differs)
include(joinpath(@__DIR__, "..", "..", "src", "MatpowerIO.jl"))
using .MatpowerIO

# ---------------------------
# Minimal CLI parsing
# ---------------------------
function parse_args(args::Vector{String})
    opts = Dict{String,Any}(
        "url" => nothing,
        "outdir" => normpath(joinpath(@__DIR__, "..", "data")),
        "to_jl" => false,
        "overwrite" => false,
        "legacy_compat" => true,
    )

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
            # default true, allow explicit flag to keep consistent style
            opts["legacy_compat"] = true
        elseif a == "--no-legacy-compat"
            opts["legacy_compat"] = false
        elseif a in ("-h", "--help")
            return opts, true
        else
            error("Unknown argument: $a (use --help)")
        end
        i += 1
    end

    return opts, false
end

function print_help()
    println("""
fetch_matpower_case.jl

Download a MATPOWER *.m case file and optionally emit a Julia NamedTuple case.

Options:
  --url <RAW_URL>           Required. Direct URL to a *.m file (raw content).
  --outdir <DIR>            Output directory (default: data/).
  --to-jl                   Also generate a *.jl NamedTuple case.
  --overwrite               Overwrite existing files.
  --legacy-compat           Keep legacy bus sorting (default).
  --no-legacy-compat        Disable legacy bus sorting.
  -h, --help                Show this help.

Example:
  julia --project=. examples/fetch_matpower_case.jl --url https://raw.githubusercontent.com/MATPOWER/matpower/master/data/case14.m --to-jl
""")
end

# ---------------------------
# Download + hashing
# ---------------------------
function fetch_matpower_case(url::AbstractString, outdir::AbstractString; overwrite::Bool=false)
    mkpath(outdir)
    fname = split(url, "/")[end]
    isempty(fname) && error("Could not derive filename from URL: $url")
    dst = joinpath(outdir, fname)

    if isfile(dst) && !overwrite
        @info "skip (exists)" file=fname path=dst
        return dst
    end

    @info "download" url=url dst=dst
    Downloads.download(url, dst)

    bytes = read(dst)
    @info "sha256" file=fname sha256=bytes2hex(sha256(bytes))
    return dst
end

# ---------------------------
# Emit Julia case from parsed mpc
# ---------------------------
function emit_julia_case(mfile::AbstractString, outdir::AbstractString;
                         legacy_compat::Bool=true, overwrite::Bool=false)

    mpc = MatpowerIO.read_case(mfile; legacy_compat=legacy_compat)

    jlfile = joinpath(outdir, "$(mpc.name).jl")
    if isfile(jlfile) && !overwrite
        @info "skip (exists)" file=basename(jlfile) path=jlfile
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

        # Print matrices in a readable literal form
        print(io, "  bus = ")
        show(io, "text/plain", mpc.bus); println(io, ",\n")

        print(io, "  gen = ")
        show(io, "text/plain", mpc.gen); println(io, ",\n")

        print(io, "  branch = ")
        show(io, "text/plain", mpc.branch); println(io, ",\n")

        if mpc.gencost === nothing
            println(io, "  gencost = nothing,")
        else
            print(io, "  gencost = ")
            show(io, "text/plain", mpc.gencost); println(io, ",")
        end

        if mpc.bus_name === nothing
            println(io, "  bus_name = nothing,")
        else
            print(io, "  bus_name = ")
            show(io, "text/plain", mpc.bus_name); println(io, ",")
        end

        println(io, ")")
    end

    @info "wrote" file=basename(jlfile) path=jlfile
    return jlfile
end

# Helpers
basename(path::AbstractString) = splitpath(path)[end]

# ---------------------------
# Main
# ---------------------------
function main()
    opts, want_help = parse_args(collect(ARGS))
    if want_help
        print_help()
        return
    end

    url = opts["url"]
    url === nothing && (print_help(); error("Missing required --url"))

    outdir = opts["outdir"]
    overwrite = opts["overwrite"]
    to_jl = opts["to_jl"]
    legacy_compat = opts["legacy_compat"]

    mfile = fetch_matpower_case(url, outdir; overwrite=overwrite)

    if to_jl
        emit_julia_case(mfile, outdir; legacy_compat=legacy_compat, overwrite=overwrite)
    end
end

main()

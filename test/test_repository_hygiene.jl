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

# file: test/test_repository_hygiene.jl
# purpose: scans tracked source, docs, test, and example files for forbidden
#          terminology tokens and fails with a bounded hit report; also
#          asserts that every src/ Julia file is listed on exactly one
#          reference page or on the explicit exclusion list below
using Test
using Unicode
import SHA

const _HYGIENE_ROOTS = ("docs/src", "src", "test", "examples")
const _HYGIENE_EXCLUDED_PARTS = Set(["docs/build", "examples/_out", "results", "coverage", ".git", ".julia", "node_modules"])
const _HYGIENE_EXTENSIONS = Set([".jl", ".md", ".toml", ".yaml", ".yml", ".json"])

function _repository_root()
    return normpath(joinpath(@__DIR__, ".."))
end

function _is_hygiene_excluded(rel::AbstractString)::Bool
    rel_norm = replace(rel, '\\' => '/')
    for part in _HYGIENE_EXCLUDED_PARTS
        rel_norm == part && return true
        startswith(rel_norm, part * "/") && return true
    end
    occursin(r"(^|/)(tmp|temp|\.tmp)(/|$)"i, rel_norm) && return true
    return false
end

function _tracked_hygiene_files(repo::AbstractString)::Vector{String}
    files = String[]
    git_dir = joinpath(repo, ".git")
    if isdir(git_dir)
        out = read(`git -C $repo ls-files`, String)
        append!(files, split(chomp(out), '\n'))
    else
        for root in _HYGIENE_ROOTS
            absroot = joinpath(repo, root)
            isdir(absroot) || continue
            for (dir, _, names) in walkdir(absroot)
                for name in names
                    push!(files, relpath(joinpath(dir, name), repo))
                end
            end
        end
    end
    return sort!(unique(filter(files) do rel
        rel = replace(rel, '\\' => '/')
        any(root -> rel == root || startswith(rel, root * "/"), _HYGIENE_ROOTS) || return false
        _is_hygiene_excluded(rel) && return false
        lowercase(splitext(rel)[2]) in _HYGIENE_EXTENSIONS
    end))
end

function _hygiene_normalize(value::AbstractString)::String
    return lowercase(Unicode.normalize(String(value), :NFC))
end

# A name that must not appear in the repository, stored only as SHA-256 of
# its normalized spellings so this file does not carry it either. Every
# window of 7 or 8 characters that starts with "sch" is hashed and compared.
const _HYGIENE_FORBIDDEN_SHA256 = Set([
    "f894e2aef608cf5635b2559c69844e9ee620a270e0df1ae9a25f6017416b7d9a",
    "2c43dcd8d59a7558af96a78da634bfb2dbd93fef82124b7109a7a48d73f8fa5c",
    "9bf59152a9eaa11a7b5fe14456b279a33fdf0b1277e8b190de15622242d09f20",
])

function _hygiene_forbidden_hit(text::AbstractString)::Bool
    chars = collect(text)
    n = length(chars)
    for i in 1:(n - 6)
        (chars[i] == 's' && chars[i + 1] == 'c' && chars[i + 2] == 'h') || continue
        for len in (7, 8)
            i + len - 1 <= n || continue
            bytes2hex(SHA.sha256(String(chars[i:(i + len - 1)]))) in _HYGIENE_FORBIDDEN_SHA256 && return true
        end
    end
    return false
end

function _bounded_hygiene_failure(hits::Vector{String}; limit::Int=20)::String
    shown = first(hits, min(limit, length(hits)))
    lines = ["repository hygiene found $(length(hits)) forbidden terminology hit(s):"]
    append!(lines, shown)
    extra = length(hits) - length(shown)
    extra > 0 && push!(lines, "... $(extra) additional hits omitted")
    msg = join(lines, '\n')
    return sizeof(msg) > 16 * 1024 ? String(take!(IOBuffer(codeunits(msg)[1:(16*1024)]))) : msg
end

# Source files that deliberately appear on no reference page: the Web UI
# server and the build tools have no pages of their own; their exported entry
# points render on reference_api.md (Application entry points section).
const _REFERENCE_PAGE_EXCLUDED_SRC = Set([
    "src/build/precompile.jl",
    "app/src/precompile.jl",
    "app/src/SparlectraApp.jl",
    "app/src/sysimage_builder.jl",
    "app/src/webui/docs.jl",
    "app/src/webui/forms.jl",
    "app/src/webui/handlers.jl",
    "app/src/webui/help_excerpts.jl",  # generated from docs/src by tools/generate_webui_help_excerpts.jl
    "app/src/webui/help_math.jl",
    "app/src/webui/operations.jl",
    "app/src/webui/options.jl",
    "app/src/webui/routes.jl",
    "app/src/webui/sysimage.jl",
    "app/src/webui/update.jl",
    "app/src/webui/views.jl",
    "app/src/webui/webui.jl",
])

# Collect every src/ file named in a reference page Pages list, mapped to the
# pages naming it. A page lists each file twice (Public API and Internals
# blocks share one Pages list), so pages are deduplicated per file.
function _reference_page_src_entries(repo::AbstractString)::Dict{String,Vector{String}}
    docsdir = joinpath(repo, "docs", "src")
    out = Dict{String,Vector{String}}()
    for name in sort!(filter(n -> startswith(n, "reference_") && endswith(n, ".md"), readdir(docsdir)))
        # CRLF-normalized: on a Windows checkout with autocrlf the page ends its
        # lines with "\r\n", and `$` in multiline mode does not match before the
        # "\r"; every Pages entry then went unseen and every tracked file was
        # reported as "on no reference page"
        text = replace(read(joinpath(docsdir, name), String), "\r\n" => "\n")
        for m in eachmatch(r"^\s*\"((?:app/)?src/[^\"]+\.jl)\",$"m, text)
            pages = get!(out, String(m.captures[1]), String[])
            name in pages || push!(pages, name)
        end
    end
    return out
end

# Every tracked src/ Julia file must be on exactly one reference page or on
# the exclusion list; stale Pages entries and stale exclusions fail too, so
# a moved or deleted file cannot silently keep a dead reference.
function _reference_page_coverage_violations(repo::AbstractString)::Vector{String}
    entries = _reference_page_src_entries(repo)
    tracked = filter(f -> endswith(f, ".jl"), split(chomp(read(`git -C $repo ls-files src app/src`, String)), "\n"))
    violations = String[]
    for rel in tracked
        listed = get(entries, rel, String[])
        excluded = rel in _REFERENCE_PAGE_EXCLUDED_SRC
        if excluded && !isempty(listed)
            push!(violations, "$(rel): excluded from reference pages but listed on $(join(listed, ", "))")
        elseif !excluded && isempty(listed)
            push!(violations, "$(rel): on no reference page and not on the exclusion list")
        elseif length(listed) > 1
            push!(violations, "$(rel): listed on more than one reference page ($(join(listed, ", ")))")
        end
    end
    for rel in sort!(collect(keys(entries)))
        rel in tracked || push!(violations, "$(rel): reference page entry has no tracked source file")
    end
    for rel in sort!(collect(_REFERENCE_PAGE_EXCLUDED_SRC))
        rel in tracked || push!(violations, "$(rel): exclusion list entry has no tracked source file")
    end
    return violations
end

# The detached-docstring check moved to test/test_docstrings.jl (fast
# profile, parser based, all of src/, app/src/, ext/, test/, examples/).

function run_repository_hygiene_tests()
    @testset "repository hygiene" begin (function ()
        repo = _repository_root()
        hits = String[]
        for rel in _tracked_hygiene_files(repo)
            rel == "test/test_repository_hygiene.jl" && continue
            normalized_rel = _hygiene_normalize(replace(rel, '\\' => '/'))
            # a tracked path can be absent from the worktree (staged rename or
            # delete); an absent file has no content to scan, its path is checked
            text = isfile(joinpath(repo, rel)) ? _hygiene_normalize(read(joinpath(repo, rel), String)) : ""
            _hygiene_forbidden_hit(normalized_rel) && push!(hits, "$(rel): path contains forbidden token")
            _hygiene_forbidden_hit(text) && push!(hits, "$(rel): content contains forbidden token")
        end
        if !isempty(hits)
            error(_bounded_hygiene_failure(hits))
        end
        # the library must not know the application package: no SparlectraApp
        # symbol below src/, the dependency points one way only
        boundary = String[]
        for f in filter(f -> endswith(f, ".jl"), split(chomp(read(`git -C $repo ls-files src`, String)), "\n"))
            occursin("SparlectraApp", read(joinpath(repo, f), String)) && push!(boundary, f)
        end
        isempty(boundary) || error(join(["library files that name SparlectraApp:"; boundary], "\n"))
        coverage = _reference_page_coverage_violations(repo)
        if !isempty(coverage)
            error(join(["reference page coverage violated:"; coverage], "\n"))
        end
        # every file the Web UI serves from its static directory is tracked:
        # an asset that exists only on the machine it was created on (an
        # ignore rule such as *.png catches it without a word) is served
        # here and answers 404 on every other checkout
        static_rel = "app/src/webui/static"
        static_tracked = split(chomp(read(`git -C $repo ls-files $static_rel`, String)), "\n")
        static_missing = String[string(static_rel, "/", name) for name in readdir(joinpath(repo, static_rel)) if !(string(static_rel, "/", name) in static_tracked)]
        isempty(static_missing) || error(join(["Web UI static files that are not tracked (check .gitignore):"; static_missing], "\n"))
        # additional check files next to this one (name pattern private_*.jl)
        # are included when present and named in the group report
        extra = sort(filter(f -> startswith(f, "private_") && endswith(f, ".jl"), readdir(@__DIR__)))
        if !isempty(extra)
            println("      additional check file(s) ran: ", join(extra, ", "))
            for f in extra
                include(joinpath(@__DIR__, f))
            end
            # the include above defines the function in a NEWER world than this
            # frame: both the binding lookup and the call go through
            # invokelatest (the same pattern as run_entry in runtests.jl), or
            # Julia 1.12+ warns about the binding access and a future version
            # errors
            checker = Base.invokelatest(getfield, @__MODULE__, :run_private_boundary_check)
            @test Base.invokelatest(checker, repo) == true
        end
        @test true
    end)() end
    # an own testset, so an error in the checks above does not skip these
    @testset "Web UI layering and demo case files" begin (function ()
        repo = _repository_root()
        # layering: the Web UI reaches the solver only through the service
        # layer, so no file under app/src/webui calls a power-flow entry point
        webui_src = filter(f -> endswith(f, ".jl"), split(chomp(read(`git -C $repo ls-files app/src/webui`, String)), "\n"))
        solver_calls = ["run" * "pf!", "run" * "pf_rectangular!", "run_complex_nr_" * "rectangular"]
        layering = String[f for f in webui_src if any(name -> occursin(name, read(joinpath(repo, f), String)), solver_calls)]
        @test isempty(layering)
        # the shipped demo cases: provenance in the file, a case config that
        # pins the tolerance (machine-neutral), the bad-data twin, and the
        # size budget (420 kB per bundle since the place names, 900 kB total)
        scf_dir = joinpath(repo, "data", "scf")
        for name in ("sp_case5", "sp_case14", "sp_case60", "sp_case188")
          scf_path = joinpath(scf_dir, "$(name).scf.json")
          @test isfile(scf_path)
          @test occursin("Self-built Sparlectra demo network", read(scf_path, String))
          @test occursin("scope: case", read(Sparlectra.case_config_path(scf_path), String))
          @test Sparlectra.load_case_config(scf_path)["power_flow.tol"] == 1.0e-8
          bundle = filter(f -> startswith(f, name), readdir(scf_dir))
          @test sum(filesize(joinpath(scf_dir, f)) for f in bundle) < 420_000
          baddata = joinpath(scf_dir, "$(name).measurements.baddata.csv")
          @test isfile(baddata) && occursin("baddata", read(baddata, String))
        end
        @test sum(filesize(joinpath(scf_dir, f)) for f in filter(f -> startswith(f, "sp_case"), readdir(scf_dir))) < 900_000
    end)() end
end
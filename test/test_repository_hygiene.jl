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

function _bounded_hygiene_failure(hits::Vector{String}; limit::Int = 20)::String
  shown = first(hits, min(limit, length(hits)))
  lines = ["repository hygiene found $(length(hits)) forbidden terminology hit(s):"]
  append!(lines, shown)
  extra = length(hits) - length(shown)
  extra > 0 && push!(lines, "... $(extra) additional hits omitted")
  msg = join(lines, '\n')
  return sizeof(msg) > 16 * 1024 ? String(take!(IOBuffer(codeunits(msg)[1:16 * 1024]))) : msg
end

# Source files that deliberately appear on no reference page: the Web UI
# server and the build tools have no pages of their own; their exported entry
# points render on reference_api.md (Application entry points section).
const _REFERENCE_PAGE_EXCLUDED_SRC = Set([
  "src/build/precompile.jl",
  "src/build/sysimage_builder.jl",
  "src/webui/docs.jl",
  "src/webui/forms.jl",
  "src/webui/handlers.jl",
  "src/webui/operations.jl",
  "src/webui/options.jl",
  "src/webui/routes.jl",
  "src/webui/sysimage.jl",
  "src/webui/views.jl",
  "src/webui/webui.jl",
])

# Collect every src/ file named in a reference page Pages list, mapped to the
# pages naming it. A page lists each file twice (Public API and Internals
# blocks share one Pages list), so pages are deduplicated per file.
function _reference_page_src_entries(repo::AbstractString)::Dict{String,Vector{String}}
  docsdir = joinpath(repo, "docs", "src")
  out = Dict{String,Vector{String}}()
  for name in sort!(filter(n -> startswith(n, "reference_") && endswith(n, ".md"), readdir(docsdir)))
    for m in eachmatch(r"^\s*\"(src/[^\"]+\.jl)\",$"m, read(joinpath(docsdir, name), String))
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
  tracked = filter(f -> endswith(f, ".jl"), split(chomp(read(`git -C $repo ls-files src`, String)), "\n"))
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

# A docstring is attached to the definition that FOLLOWS it. Put a blank
# line in between and Julia keeps the string as a standalone expression: the
# binding stays undocumented, without an error and without a warning. Nothing
# in a test run sees it; only the documentation build does, and only when a
# page happens to link the binding with @ref. On 2026-09-07 exactly that took
# the docs build down (`ensure_casefile`, linked from three pages), and the
# same blank line had quietly eaten the docstrings of `addBranch!`,
# `addZeroInjectionMeasurements!` and the keyword method of
# `runpf_rectangular!` months earlier.
#
# The scan tracks triple-quote delimiters and only treats a block as a
# docstring when its OPENING delimiter stands alone on its line: that
# separates a docstring from a multi-line string constant such as
# `const TEXT = """`, where a blank line before the next definition is
# perfectly normal.
const _DEFINITION_LINE = r"^\s*(?:@\w+[^\s]*\s+)*(?:function|macro|module|baremodule|const|abstract\s+type|primitive\s+type|mutable\s+struct|struct)\b"

function _detached_docstring_violations(repo::AbstractString)::Vector{String}
  violations = String[]
  tracked = filter(f -> endswith(f, ".jl"), split(chomp(read(`git -C $repo ls-files src test examples`, String)), "\n"))
  for rel in tracked
    # this file documents the pattern in its own comments and fixtures
    rel == "test/test_repository_hygiene.jl" && continue
    path = joinpath(repo, rel)
    isfile(path) || continue
    lines = readlines(path)
    inside = false          # inside a triple-quoted block
    bare_opener = false     # ... and that block opened with a lone delimiter
    for (i, line) in pairs(lines)
      delimiters = length(collect(eachmatch(r"\"\"\"", line)))
      closed = false
      if isodd(delimiters)
        if inside
          inside = false
          closed = bare_opener
        else
          inside = true
          # `"""`, `raw"""` or `md"""` with nothing else on the line
          bare_opener = occursin(r"^\s*(?:[A-Za-z_]\w*)?\"\"\"\s*$", line)
        end
      end
      # a lone single-line string (`"Short doc."`) detaches the same way
      if !inside && !closed && occursin(r"^\s*\"[^\"].*[^\\]\"\s*$", line)
        closed = true
      end
      closed || continue
      j = i + 1
      while j <= length(lines) && isempty(strip(lines[j]))
        j += 1
      end
      j > length(lines) && continue
      j == i + 1 && continue  # no blank line: correctly attached
      if occursin(_DEFINITION_LINE, lines[j]) || occursin(r"^\s*[A-Za-z_]\w*[!?]?\s*\(.*\)\s*=", lines[j])
        push!(violations, "$(rel):$(i): blank line between the docstring and the definition on line $(j) detaches it")
      end
    end
  end
  return violations
end

function run_repository_hygiene_tests()
  @testset "repository hygiene" begin
    repo = _repository_root()
    forbidden_terms = _hygiene_normalize.(String[
      "sch" * "ae" * "fer",
      "sch" * "ä" * "fer",
      "sch" * "a" * "fer",
    ])
    hits = String[]
    for rel in _tracked_hygiene_files(repo)
      rel == "test/test_repository_hygiene.jl" && continue
      normalized_rel = _hygiene_normalize(replace(rel, '\\' => '/'))
      # a tracked path can be absent from the worktree (staged rename or
      # delete); an absent file has no content to scan, its path is checked
      text = isfile(joinpath(repo, rel)) ? _hygiene_normalize(read(joinpath(repo, rel), String)) : ""
      for term in forbidden_terms
        occursin(term, normalized_rel) && push!(hits, "$(rel): path contains forbidden token")
        occursin(term, text) && push!(hits, "$(rel): content contains forbidden token")
      end
    end
    if !isempty(hits)
      error(_bounded_hygiene_failure(hits))
    end
    coverage = _reference_page_coverage_violations(repo)
    if !isempty(coverage)
      error(join(["reference page coverage violated:"; coverage], "\n"))
    end
    detached = _detached_docstring_violations(repo)
    if !isempty(detached)
      error(join(["docstrings detached from their definition:"; detached], "\n"))
    end
    # The maintainer's working repository carries further checks here that a
    # published checkout has no subject for, the private-to-public boundary
    # among them. They live in their own files, found by convention rather
    # than by name, so this file never points at something a public reader
    # cannot have. Their absence is STATED, not skipped quietly: a silent
    # skip reads as a pass and hides a coverage gap, the defect class this
    # repository removed twice in the week of 2026-09-06.
    extra = sort(filter(f -> startswith(f, "private_") && endswith(f, ".jl"), readdir(@__DIR__)))
    if isempty(extra)
      # the word SKIPPED is load-bearing: the group runner surfaces exactly
      # those lines through its output capture, so the absence reaches the
      # report instead of being swallowed
      println("      SKIPPED maintainer-only repository checks: none present in this checkout. ",
              "They verify the private-to-public boundary, which has no subject in a published ",
              "tree, so nothing beyond the checks above was verified here.")
    else
      for f in extra
        include(joinpath(@__DIR__, f))
      end
      # Base.invokelatest: the include above defines the function in a NEWER
      # world than this frame, so a direct call raises "the applicable method
      # may be too new" on Julia 1.12
      @test Base.invokelatest(run_private_boundary_check, repo) == true
    end
    @test true
  end
end

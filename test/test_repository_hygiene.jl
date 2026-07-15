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
      text = _hygiene_normalize(read(joinpath(repo, rel), String))
      for term in forbidden_terms
        occursin(term, normalized_rel) && push!(hits, "$(rel): path contains forbidden token")
        occursin(term, text) && push!(hits, "$(rel): content contains forbidden token")
      end
    end
    if !isempty(hits)
      error(_bounded_hygiene_failure(hits))
    end
    @test true
  end
end

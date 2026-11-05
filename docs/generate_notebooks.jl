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

# file: docs/generate_notebooks.jl
# purpose: regenerate the committed Literate.jl outputs from docs/lit/*.jl:
#          a Documenter page under docs/src/generated/ and a notebook under
#          notebooks/. Run manually after editing a source:
#          `julia --project=docs docs/generate_notebooks.jl`. Both outputs
#          are committed; docs/make.jl does NOT call this script. As a safety
#          net, the Documentation workflow (jekyll-gh-pages.yml) reruns this
#          script on every push to main and commits stale artifacts back.

# Self-contained environment setup, same pattern as docs/make.jl: activate the
# docs project regardless of how this script is started.
using Pkg
Pkg.activate(@__DIR__)
# resolve keeps every manifest version as an explicit requirement, so a docs
# manifest that predates a compat bump of the package fails as unsatisfiable
# instead of moving; the fallback updates this untracked docs environment and
# resolves again
try
    Pkg.resolve()
catch err
    println("Pkg.resolve failed (", first(sprint(showerror, err), 120), "); updating the docs environment.")
    Pkg.update()
    Pkg.resolve()
end
Pkg.instantiate()

using Literate
using JSON3

const LIT_DIR       = joinpath(@__DIR__, "lit")
const GENERATED_DIR = joinpath(@__DIR__, "src", "generated")
const NOTEBOOK_DIR  = normpath(joinpath(@__DIR__, "..", "notebooks"))

"""
    patch_notebook_metadata!(path::AbstractString)

Rewrite the generated notebook's `metadata.kernelspec` to the generic
`{name: "julia", display_name: "Julia", language: "julia"}`. Literate stamps
the kernelspec (and `language_info.version`) with the local Julia version,
which would churn the committed file on every local Julia upgrade;
`language_info` is therefore normalized to a version-free form as well. The
workshop installer (`tools/install_workshops.jl`) replaces the generic
kernel name in its copies by the kernel it installs. The rewritten file is
re-parsed to prove it is still valid JSON.
"""
function patch_notebook_metadata!(path::AbstractString)
  nb = copy(JSON3.read(read(path, String)))
  metadata = get!(nb, :metadata, Dict{Symbol,Any}())
  metadata[:kernelspec] = Dict(:name => "julia", :display_name => "Julia", :language => "julia")
  metadata[:language_info] = Dict(:name => "julia", :file_extension => ".jl", :mimetype => "application/julia")
  open(path, "w") do io
    JSON3.pretty(io, nb)
    println(io)
  end
  reparsed = JSON3.read(read(path, String))
  ks = reparsed.metadata.kernelspec
  @assert ks.name == "julia" && ks.language == "julia" && ks.display_name == "Julia"
  return nothing
end

"""
    check_outputs(mdpath, nbpath)

Post-generation sanity checks for the `#nb` line filtering: the setup cell
(the check that Sparlectra is installed in the workshop environment) must be
present in the notebook and absent from the Documenter page, and no
notebook may install packages itself (the workshop installer does that).
"""
function check_outputs(mdpath::AbstractString, nbpath::AbstractString)
  nb = JSON3.read(read(nbpath, String))
  sources = [join(cell.source) for cell in nb.cells if cell.cell_type == "code"]
  @assert any(occursin("Base.find_package(\"Sparlectra\")", s) for s in sources) "setup cell missing from $(nbpath)"
  @assert !any(occursin(r"Pkg\.(add|develop|activate)\(", s) for s in sources) "a code cell of $(nbpath) installs packages; the workshop installer does that"
  @assert !occursin("find_package", read(mdpath, String)) "setup cell leaked into $(mdpath)"
  return nothing
end

function generate()
  for source in filter(endswith(".jl"), sort(readdir(LIT_DIR)))
    srcpath = joinpath(LIT_DIR, source)
    stem    = first(splitext(source))
    Literate.markdown(srcpath, GENERATED_DIR; credit = false)
    Literate.notebook(srcpath, NOTEBOOK_DIR; execute = false, credit = false)
    nbpath = joinpath(NOTEBOOK_DIR, stem * ".ipynb")
    mdpath = joinpath(GENERATED_DIR, stem * ".md")
    patch_notebook_metadata!(nbpath)
    check_outputs(mdpath, nbpath)
    println("generated: $(relpath(mdpath)) and $(relpath(nbpath))")
  end
end

Base.invokelatest(generate)

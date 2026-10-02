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

# file: test/test_docstrings.jl
# purpose: fast-profile check that every docstring sits directly on the
#          definition it documents (task_docstrings_0302). Replaces the
#          line-based scan of the config-profile hygiene group, which covered
#          src/test/examples only and missed one-line docstrings, a helper
#          placed directly below a docstring and doubled docstrings.

# Measured on Julia 1.13.1 (task_docstrings_0302, item 4): a docstring is
# attached only when the definition follows on the very next line. A blank
# line or a comment line in between detaches it (the text is discarded and
# `Meta.parseall` shows a bare String at block level); indentation of the
# definition does not matter. Any other code in between (a helper, a const)
# takes the docstring itself, so the documented function stays undocumented
# and the helper carries text that names another function.
#
# The check therefore parses every tracked .jl file below src/, app/src/,
# ext/, test/ and examples/ and reports
#   1. a bare string at block level (module body, top level, begin/if
#      blocks): a detached docstring or a stray string;
#   2. a docstring whose signature line (its first indented line, the
#      Documenter convention "    name(args)") names a different binding
#      than the definition it is attached to: a docstring that slid onto a
#      helper.

_doc_is_docmacro(e) = e isa Expr && e.head === :macrocall && e.args[1] == GlobalRef(Core, Symbol("@doc"))
_doc_is_bare_string(e) = e isa String || (e isa Expr && e.head === :string) ||
                         (e isa Expr && e.head === :macrocall && e.args[1] isa Symbol && endswith(String(e.args[1]), "_str"))

# the name a definition binds (function, short form, const, struct, ...)
function _doc_defname(e)
  e isa Symbol && return e
  e isa QuoteNode && return e.value
  e isa Expr || return nothing
  e.head in (:function, :macro, :(=), :call, :where, :curly, :<:) && return _doc_defname(e.args[1])
  e.head === :(::) && return length(e.args) == 2 ? _doc_defname(e.args[1]) : nothing
  e.head === :. && return _doc_defname(e.args[end])
  e.head in (:const, :global) && return _doc_defname(e.args[1])
  e.head === :struct && return _doc_defname(e.args[2])
  e.head in (:abstract, :primitive) && return _doc_defname(e.args[1])
  return nothing
end

# the binding a docstring's first indented line names ("    foo(x) -> y",
# "    Module.foo!", "    Foo"); nothing when the docstring has no such line
function _doc_signature_name(s::AbstractString)
  for l in split(s, '\n')
    m = match(r"^    ([A-Za-z_][A-Za-z0-9_!]*(?:\.[A-Za-z_][A-Za-z0-9_!]*)*)", l)
    m === nothing || return Symbol(split(m.captures[1], '.')[end])
    isempty(strip(l)) || return nothing
  end
  return nothing
end

function _doc_walk!(violations, args, rel, line)
  ln = line
  for a in args
    if a isa LineNumberNode
      ln = a.line
    elseif _doc_is_bare_string(a)
      push!(violations, "$(rel):$(ln): a bare string at block level: a docstring detached from its definition (blank or comment line in between) or a stray string")
    elseif _doc_is_docmacro(a)
      text = a.args[3] isa String ? a.args[3] : nothing
      sig = text === nothing ? nothing : _doc_signature_name(text)
      name = _doc_defname(a.args[end])
      if sig !== nothing && name !== nothing && sig != name
        push!(violations, "$(rel):$(ln): the docstring names $(sig) but is attached to $(name) (code between the docstring and its definition)")
      end
    elseif a isa Expr && a.head in (:module, :baremodule)
      _doc_walk!(violations, a.args[3].args, rel, ln)
    elseif a isa Expr && a.head in (:block, :toplevel)
      _doc_walk!(violations, a.args, rel, ln)
    elseif a isa Expr && a.head in (:if, :elseif)
      for b in a.args[2:end]
        b isa Expr || continue
        b.head === :block && _doc_walk!(violations, b.args, rel, ln)
        b.head === :elseif && _doc_walk!(violations, [b], rel, ln)
      end
    end
  end
  return violations
end

"""
    docstring_violations(repo; dirs) -> Vector{String}

Every detached or misplaced docstring below `dirs` (tracked `.jl` files of
the repository at `repo`), one line each with file and line. Empty when all
docstrings sit directly on the definition they document.
"""
function docstring_violations(repo::AbstractString; dirs = ["src", "app/src", "ext", "test", "examples"])::Vector{String}
  violations = String[]
  files = filter(f -> endswith(f, ".jl"), split(chomp(read(`git -C $repo ls-files $dirs`, String)), "\n"))
  for rel in files
    path = joinpath(repo, rel)
    isfile(path) || continue
    _doc_walk!(violations, Meta.parseall(read(path, String); filename = path).args, rel, 0)
  end
  return violations
end

function run_docstring_tests()
  @testset "docstrings sit on their definitions" begin (function ()
    repo = normpath(joinpath(@__DIR__, ".."))
    violations = docstring_violations(repo)
    isempty(violations) || println("      detached or misplaced docstrings:\n        ", join(violations, "\n        "))
    @test isempty(violations)
    # the check itself: the three layouts the Julia 1.13 experiment found
    # broken are reported, the attached one is not
    probe(code) = _doc_walk!(String[], Meta.parseall(code).args, "probe.jl", 0)
    @test isempty(probe("\"\"\"\n    f(x)\n\ndoc\n\"\"\"\nf(x) = x\n"))
    @test length(probe("\"\"\"\n    f(x)\n\ndoc\n\"\"\"\n\nf(x) = x\n")) == 1
    @test length(probe("\"\"\"\n    f(x)\n\ndoc\n\"\"\"\n# note\nf(x) = x\n")) == 1
    @test length(probe("\"\"\"\n    f(x)\n\ndoc\n\"\"\"\n_h() = 1\nf(x) = x\n")) == 1
  end)() end
  return true
end

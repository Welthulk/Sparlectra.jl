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

# file: test/test_workshops.jl
# purpose: run every Literate workshop (docs/lit/workshop_*.jl) top to
#          bottom in a fresh module, so the notebooks and the docs pages
#          generated from them cannot drift from the library. Each workshop
#          carries an @assert next to every printed number; a failing
#          assert fails this test. The workshops use shipped cases only.

# a fresh module that knows `include` and `eval` the way Main of a notebook
# kernel and Documenter's sandbox do (a plain Module() has neither): the
# first cell of every workshop includes its warm-up file
function _workshop_module(name::Symbol)::Module
  sym = gensym(name)
  return Core.eval(Main, :(baremodule $sym
    using Base
    eval(x) = Core.eval($sym, x)
    include(x) = Base.include($sym, abspath(x))
  end))
end

function run_workshop_tests()
  @testset "Workshops run with their assertions" begin (function ()
    lit_dir = abspath(joinpath(dirname(@__DIR__), "docs", "lit"))
    workshops = sort!(filter(f -> startswith(f, "workshop_") && endswith(f, ".jl"), readdir(lit_dir)))
    @test !isempty(workshops)
    ran = String[]
    for f in workshops
      path = joinpath(lit_dir, f)
      # a fresh module per workshop: no name leaks between them, exactly
      # like a notebook kernel started for that one workshop
      mod = _workshop_module(Symbol("Workshop_", splitext(f)[1]))
      # the workshops print their results and log their decisions; the test
      # run keeps both quiet (a failing assert still throws)
      ok = try
        with_logger(NullLogger()) do
          redirect_stdout(devnull) do
            Base.include(mod, path)
          end
        end
        true
      catch err
        @error "workshop failed" workshop = f exception = (err, catch_backtrace())
        false
      end
      @test ok
      ok && push!(ran, f)
    end
    println("      workshops: RAN ", length(ran), " of ", length(workshops), " (", join(replace.(ran, "workshop_" => "", ".jl" => ""), ", "), ")")
  end)() end
  @testset "Workshop uses lines are current" begin (function ()
    # tools/workshop_uses_comments.jl keeps the "## uses:" line of every
    # code cell; a cell right under a markdown line once sent its rewrite
    # into an endless loop (26 GB before it was killed)
    tool = Module(:WorkshopUsesTool)
    Base.include(tool, joinpath(dirname(@__DIR__), "tools", "workshop_uses_comments.jl"))
    rewrite = getfield(tool, :rewrite)
    tight = ["# ## Head", "# text", "x = 1", "", "# more", "y = x + 1"]
    @test Base.invokelatest(rewrite, tight) == ["# ## Head", "# text", "x = 1", "", "# more", "## uses: x (Head)", "@isdefined(x) || error(\"Run the section \\\"Head\\\" first: it sets up x.\")", "y = x + 1"]
    # a name from the warm-up cell (the one with `using Sparlectra`) gets
    # the uses line but no guard: the warm-up is the one cell never skipped
    warm = ["# ## Warm-up", "using Sparlectra", "h(n) = n", "# **Example 1**", "z = h(2)"]
    @test Base.invokelatest(rewrite, warm)[end-1:end] == ["## uses: h (warm-up cell)", "z = h(2)"]
    @test Base.invokelatest(rewrite, Base.invokelatest(rewrite, tight)) == Base.invokelatest(rewrite, tight)
    # a keyword name is not a use: `f(x = 2)` does not read x
    @test Base.invokelatest(rewrite, ["# a", "x = 1", "# b", "f(; x) = x", "# c", "f(x = 2)"])[end] == "f(x = 2)"
    lit_dir = joinpath(dirname(@__DIR__), "docs", "lit")
    for f in sort!(filter(f -> startswith(f, "workshop_") && endswith(f, ".jl"), readdir(lit_dir)))
      lines = readlines(joinpath(lit_dir, f))
      @test (f, Base.invokelatest(rewrite, lines) == lines) == (f, true)
    end
  end)() end
  @testset "Workshop notebooks run cell by cell" begin (function ()
    # what a participant runs is the notebook, not the Literate script:
    # every code cell of every committed notebook, in order, in a fresh
    # module, like one kernel per notebook. The install cell is skipped (the
    # package under test is already loaded); a failing cell names the
    # notebook and the cell number of the Colab view (1-based)
    nb_dir = abspath(joinpath(dirname(@__DIR__), "notebooks"))
    notebooks = sort!(filter(f -> startswith(f, "workshop_") && endswith(f, ".ipynb"), readdir(nb_dir)))
    @test !isempty(notebooks)
    ran = String[]
    for f in notebooks
      nb = Sparlectra.scf_json_parse(read(joinpath(nb_dir, f), String))
      mod = _workshop_module(Symbol("Notebook_", splitext(f)[1]))
      failed_cell = 0
      for (k, cell) in enumerate(nb["cells"])
        cell["cell_type"] == "code" || continue
        code = join(cell["source"])
        # the install cell: only its `using` lines run (a notebook may load
        # its packages there, the APSLF workshop did)
        if occursin("Pkg.add(", code)
          code = join(filter(l -> startswith(l, "using ") && !startswith(l, "using Pkg"), split(code, '\n')), "\n")
        end
        ok = try
          with_logger(NullLogger()) do
            redirect_stdout(devnull) do
              include_string(mod, code, string(f, " cell ", k))
            end
          end
          true
        catch err
          @error "notebook cell failed" notebook = f cell = k exception = (err, catch_backtrace())
          false
        end
        ok || (failed_cell = k; break)
      end
      @test (f, failed_cell) == (f, 0)
      failed_cell == 0 && push!(ran, f)
    end
    println("      notebooks: RAN ", length(ran), " of ", length(notebooks))
  end)() end
  return true
end

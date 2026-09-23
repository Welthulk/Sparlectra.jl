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
      mod = Module(Symbol("Workshop_", splitext(f)[1]))
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
  return true
end

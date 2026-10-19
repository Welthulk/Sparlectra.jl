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

# file: test/test_legacy_matpower_example.jl
# purpose: legacy profile: the testsets of test/test_matpower_example.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_matpower_example.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_matpower_example_tests()
  # cut from the runner function
  @testset "MATPOWER example cleanup guard" begin (function ()
    source = read(joinpath(@__DIR__, "..", "examples", "powerflow", "matpower_import.jl"), String)
    @test !occursin("load_yaml_dict", source)
    @test !occursin("_normalize_matpower_example_config", source)
    @test !occursin("_copy_nested_yaml!", source)
    @test !occursin("_nested_yaml_get", source)
    @test !occursin("Dict{String,Any}", source)
    @test occursin("Sparlectra.run_matpower_case", source)
    @test occursin("run_example(main)", source)
    @test occursin("Base.invokelatest(getfield(@__MODULE__, :_ensure_julia_threads_for_script)", source)
    @test occursin("Base.invokelatest(getfield(@__MODULE__, :_resolve_julia_threads_request)", source)
    @test occursin("Base.invokelatest(", source)

    example_path = joinpath(@__DIR__, "..", "examples", "powerflow", "matpower_import.jl")
    old_no_main = get(ENV, "SPARLECTRA_MATPOWER_IMPORT_NO_MAIN", nothing)
    ENV["SPARLECTRA_MATPOWER_IMPORT_NO_MAIN"] = "1"
    mod = Module(:MatpowerImportExampleSmoke)
    Core.eval(mod, :(using Sparlectra))
    Core.eval(mod, :(include(path) = Base.include($mod, path)))
    try
      Base.include(mod, example_path)
      @test isdefined(mod, :main)
    finally
      if isnothing(old_no_main)
        delete!(ENV, "SPARLECTRA_MATPOWER_IMPORT_NO_MAIN")
      else
        ENV["SPARLECTRA_MATPOWER_IMPORT_NO_MAIN"] = old_no_main
      end
    end
  end)() end

  return true
end

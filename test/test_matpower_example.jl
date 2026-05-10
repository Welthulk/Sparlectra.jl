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

# file: test/test_matpower_example.jl

function run_matpower_example_tests()
  @testset "MATPOWER example keyword forwarding" begin
    source = read(joinpath(@__DIR__, "..", "src", "examples", "matpower_import.jl"), String)
    signature_match = match(r"function bench_run_acpflow\(;.*\)", source)
    @test !isnothing(signature_match)
    signature = signature_match.match

    expected_keywords = [
      "autodamp",
      "autodamp_min",
      "start_projection",
      "start_projection_try_dc_start",
      "start_projection_try_blend_scan",
      "start_projection_blend_lambdas",
      "start_projection_dc_angle_limit_deg",
    ]

    for keyword in expected_keywords
      @test occursin(keyword, signature)
      @test occursin(keyword * " = " * keyword, source)
    end

    @test occursin("Base.invokelatest(() -> bench_run_acpflow(;", source)
    @test occursin("Base.invokelatest(() -> main())", source)
    @test occursin("function _enable_pq_gen_controllers_for_method(method::Symbol, requested::Bool)::Bool", source)
    @test occursin("return requested && method === :rectangular", source)
    @test occursin("!enable_pq_gen_controllers && m === :rectangular", source)
  end
  return true
end

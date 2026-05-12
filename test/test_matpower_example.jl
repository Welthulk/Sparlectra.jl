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

using Sparlectra
using Test

function run_matpower_example_tests()
  @testset "MATPOWER example keyword forwarding" begin
    source = read(joinpath(@__DIR__, "..", "src", "examples", "matpower_import.jl"), String)
    signature_match = match(r"(?s)function bench_run_acpflow\(;.*?\n\)", source)
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

    @test occursin("Base.invokelatest(", source)
    @test occursin("() -> getfield(@__MODULE__, :bench_run_acpflow)(;", source)
    @test occursin("Base.invokelatest(() -> getfield(@__MODULE__, :main)())", source)
    @test occursin("SPARLECTRA_MATPOWER_IMPORT_NO_MAIN", source)
    @test occursin("function _enable_pq_gen_controllers_for_method(method::Symbol, requested::Bool)::Bool", source)
    @test occursin("return requested && method === :rectangular", source)
    @test occursin("!enable_pq_gen_controllers && m === :rectangular", source)

    example_path = joinpath(@__DIR__, "..", "src", "examples", "matpower_import.jl")
    old_no_main = get(ENV, "SPARLECTRA_MATPOWER_IMPORT_NO_MAIN", nothing)
    ENV["SPARLECTRA_MATPOWER_IMPORT_NO_MAIN"] = "1"
    mod = Module(:MatpowerImportExampleSmoke)
    try
      Base.include(mod, example_path)
      cfg = Base.invokelatest(() -> getfield(mod, :bench_config_for_case)("case14.m", Dict{String,Any}(
        "autodamp" => true,
        "autodamp_min" => 0.002,
        "start_projection" => true,
        "start_projection_try_dc_start" => false,
        "start_projection_try_blend_scan" => false,
        "start_projection_blend_lambdas" => [0.1, 0.9],
        "start_projection_dc_angle_limit_deg" => 45.0,
      )))
      @test cfg.autodamp === true
      @test cfg.autodamp_min == 0.002
      @test cfg.start_projection === true
      @test cfg.start_projection_try_dc_start === false
      @test cfg.start_projection_try_blend_scan === false
      @test cfg.start_projection_blend_lambdas == [0.1, 0.9]
      @test cfg.start_projection_dc_angle_limit_deg == 45.0
      @test Base.invokelatest(() -> getfield(mod, :_enable_pq_gen_controllers_for_method)(:rectangular, true)) === true
      @test Base.invokelatest(() -> getfield(mod, :_enable_pq_gen_controllers_for_method)(:polar, true)) === false
      mpc_seeded = (; bus = hcat(collect(1.0:3.0), fill(1.0, 3), zeros(3, 5), [1.02, 1.01, 0.99], [0.0, -1.0, -2.0]))
      @test Base.invokelatest(() -> getfield(mod, :_warn_if_flatstart_uses_only_voltage_setpoints)("case1951rte.m", (; opt_flatstart = true), mpc_seeded)) === nothing
      @test_logs (:info, r"opt_flatstart=false uses stored MATPOWER voltage magnitudes and angles") Base.invokelatest(() -> getfield(mod, :_warn_if_flatstart_uses_only_voltage_setpoints)("case1951rte.m", (; opt_flatstart = false), mpc_seeded))

      logfile, io = mktemp()
      close(io)
      result = Base.invokelatest(() -> getfield(mod, :bench_run_acpflow)(;
        casefile = "case14.m",
        methods = Symbol[],
        mpc = nothing,
        logfile = logfile,
        show_diff = false,
        tol_vm = 0.02,
        tol_va = 2.0,
        autodamp = true,
        autodamp_min = 0.002,
        start_projection = true,
        start_projection_try_dc_start = false,
        start_projection_try_blend_scan = false,
        start_projection_blend_lambdas = [0.1, 0.9],
        start_projection_dc_angle_limit_deg = 45.0,
        show_once = false,
        benchmark = false,
      ))
      @test result == Dict{Symbol,Any}()
      rm(logfile; force = true)
    finally
      if isnothing(old_no_main)
        delete!(ENV, "SPARLECTRA_MATPOWER_IMPORT_NO_MAIN")
      else
        ENV["SPARLECTRA_MATPOWER_IMPORT_NO_MAIN"] = old_no_main
      end
    end
  end

  @testset "MATPOWER local Julia case resolution" begin
    mktempdir() do dir
      case_path = joinpath(dir, "local_case.jl")
      write(case_path, """
      (
        name = \"local_case\",
        baseMVA = 100.0,
        bus = [1 3 0 0 0 0 1 1.0 0.0 110 1 1.1 0.9],
        gen = [1 0 0 10 -10 1.0 100 1 10 0],
        branch = zeros(0, 13),
        gencost = nothing,
        bus_name = nothing,
      )
      """)

      resolved = Sparlectra.FetchMatpowerCase.ensure_casefile("local_case.jl"; outdir = dir, to_jl = true, overwrite = false)
      @test resolved == abspath(case_path)

      old_pwd = pwd()
      try
        cd(dirname(@__DIR__))
        relative_path = relpath(case_path, pwd())
        mpc = Sparlectra.MatpowerIO.read_case(relative_path; legacy_compat = false)
        @test mpc.name == "local_case"
        @test size(mpc.bus) == (1, 13)
      finally
        cd(old_pwd)
      end
    end
  end
  return true
end

if abspath(PROGRAM_FILE) == @__FILE__
  Base.invokelatest(run_matpower_example_tests)
end

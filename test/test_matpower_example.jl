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

function _is_gitlab_ci()::Bool
  return lowercase(get(ENV, "GITLAB_CI", "false")) == "true"
end

function _ensure_matpower_case_for_example_tests(case_name::AbstractString)
  local_case = joinpath(Sparlectra.MPOWER_DIR, case_name)
  isfile(local_case) && return abspath(local_case)
  _is_gitlab_ci() && return nothing
  try
    return Sparlectra.FetchMatpowerCase.ensure_casefile(case_name; outdir = Sparlectra.MPOWER_DIR, to_jl = false, overwrite = false)
  catch
    return nothing
  end
end

function run_matpower_example_tests()
  @testset "MATPOWER example keyword forwarding" begin
    source = read(joinpath(@__DIR__, "..", "src", "examples", "matpower_import.jl"), String)
    normalized_source = replace(source, "\r\n" => "\n")
    signature_match = match(r"(?s)function bench_run_acpflow\(;.*?\n\)", normalized_source)
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
      "matpower_ratio",
    ]

    for keyword in expected_keywords
      @test occursin(keyword, signature)
      @test occursin(keyword * " = " * keyword, normalized_source)
    end

    @test occursin("Base.invokelatest(", normalized_source)
    @test occursin(r"Base\.invokelatest\(\s*getfield\(@__MODULE__, :bench_run_acpflow\);", normalized_source)
    @test occursin("Base.invokelatest(getfield(@__MODULE__, :main))", normalized_source)
    @test occursin("SPARLECTRA_MATPOWER_IMPORT_NO_MAIN", normalized_source)
    @test occursin("function _enable_pq_gen_controllers_for_method(method::Symbol, requested::Bool)::Bool", normalized_source)
    @test occursin("return requested && method === :rectangular", normalized_source)
    @test occursin("!enable_pq_gen_controllers && m === :rectangular", normalized_source)

    example_path = joinpath(@__DIR__, "..", "src", "examples", "matpower_import.jl")
    old_no_main = get(ENV, "SPARLECTRA_MATPOWER_IMPORT_NO_MAIN", nothing)
    ENV["SPARLECTRA_MATPOWER_IMPORT_NO_MAIN"] = "1"
    mod = Module(:MatpowerImportExampleSmoke)
    try
      Base.include(mod, example_path)

      example_yaml = joinpath(@__DIR__, "..", "src", "examples", "matpower_import.yaml.example")
      @test isfile(example_yaml)
      yaml_cfg = Base.invokelatest(() -> getfield(mod, :load_yaml_config)(example_yaml))
      example_case = String(yaml_cfg["case"])
      @test example_case == "case118.m"
      example_case_path = _ensure_matpower_case_for_example_tests(example_case)
      if isnothing(example_case_path)
        @test_skip "MATPOWER example case missing; dataset download skipped or unavailable"
      else
        @test isfile(example_case_path)
      end
      old_yaml = get(ENV, "SPARLECTRA_MATPOWER_IMPORT_YAML", nothing)
      try
        ENV["SPARLECTRA_MATPOWER_IMPORT_YAML"] = example_yaml
        @test Base.invokelatest(() -> getfield(mod, :_yaml_path_from_inputs)()) == example_yaml
      finally
        if isnothing(old_yaml)
          delete!(ENV, "SPARLECTRA_MATPOWER_IMPORT_YAML")
        else
          ENV["SPARLECTRA_MATPOWER_IMPORT_YAML"] = old_yaml
        end
      end

      cfg = Base.invokelatest(() -> getfield(mod, :bench_config_for_case)("case14.m", Dict{String,Any}(
        "autodamp" => true,
        "autodamp_min" => 0.002,
        "start_projection" => true,
        "start_projection_try_dc_start" => false,
        "start_projection_try_blend_scan" => false,
        "start_projection_blend_lambdas" => [0.1, 0.9],
        "start_projection_dc_angle_limit_deg" => 45.0,
        "diagnose_matpower_reference" => true,
        "diagnose_branch_shift_conventions" => true,
        "diagnose_maxlines" => 3,
        "matpower_shift_unit" => "rad",
        "matpower_shift_sign" => -1.0,
        "matpower_ratio" => "reciprocal",
        "log_effective_config" => true,
      )))
      @test cfg.autodamp === true
      @test cfg.autodamp_min == 0.002
      @test cfg.start_projection === true
      @test cfg.start_projection_try_dc_start === false
      @test cfg.start_projection_try_blend_scan === false
      @test cfg.start_projection_blend_lambdas == [0.1, 0.9]
      @test cfg.start_projection_dc_angle_limit_deg == 45.0
      @test cfg.diagnose_matpower_reference === true
      @test cfg.diagnose_branch_shift_conventions === true
      @test cfg.diagnose_maxlines == 3
      @test cfg.matpower_shift_unit == "rad"
      @test cfg.matpower_shift_sign == -1.0
      @test cfg.matpower_ratio == "reciprocal"
      @test cfg.log_effective_config === true
      @test occursin("matpower_shift_sign", normalized_source)
      @test occursin("matpower_ratio", normalized_source)
      @test occursin("_print_effective_config", normalized_source)
      @test occursin("function _print_matpower_reference_diagnostics", normalized_source)
      @test occursin("Branch-shift convention scan", normalized_source)
      @test Base.invokelatest(() -> getfield(mod, :_enable_pq_gen_controllers_for_method)(:rectangular, true)) === true
      @test Base.invokelatest(() -> getfield(mod, :_enable_pq_gen_controllers_for_method)(:polar, true)) === false
      skipped_compare_summary = Base.invokelatest(() -> getfield(mod, :_show_once_summary_row)(:rectangular, (; converged = true), nothing, false; compare_available = false))
      @test skipped_compare_summary.converged === true
      @test isnan(skipped_compare_summary.max_dvm)
      @test isnan(skipped_compare_summary.max_dva)
      @test skipped_compare_summary.cmp_ok === false
      compared_summary = Base.invokelatest(() -> getfield(mod, :_show_once_summary_row)(:rectangular, (; converged = true), (; max_dvm = 0.01, max_dva = 0.2), true; compare_available = true))
      @test compared_summary.converged === true
      @test compared_summary.max_dvm == 0.01
      @test compared_summary.max_dva == 0.2
      @test compared_summary.cmp_ok === true
      mpc_seeded = (; bus = hcat(collect(1.0:3.0), fill(1.0, 3), zeros(3, 5), [1.02, 1.01, 0.99], [0.0, -1.0, -2.0]))
      @test Base.invokelatest(() -> getfield(mod, :_warn_if_flatstart_uses_only_voltage_setpoints)("case1951rte.m", (; opt_flatstart = true), mpc_seeded)) === nothing
      @test_logs (:info, r"opt_flatstart=false uses stored MATPOWER voltage magnitudes and angles") Base.invokelatest(() -> getfield(mod, :_warn_if_flatstart_uses_only_voltage_setpoints)("case1951rte.m", (; opt_flatstart = false), mpc_seeded))

      diagnostic_mpc = (;
        baseMVA = 100.0,
        bus = [1.0 3.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 110.0 1.0 1.1 0.9;
               2.0 1.0 0.0 0.0 0.0 0.0 1.0 1.0 -0.1 110.0 1.0 1.1 0.9],
        gen = [1.0 0.0 0.0 0.0 0.0 1.0 100.0 1.0 0.0 0.0 zeros(1, 11)...],
        branch = [1.0 2.0 0.01 0.10 0.02 0.0 0.0 0.0 1.0 0.01 1.0 0.0 0.0],
      )
      diag_io = IOBuffer()
      @test Base.invokelatest(() -> getfield(mod, :_print_matpower_reference_diagnostics)(diag_io, diagnostic_mpc; matpower_shift_sign = -1.0, matpower_shift_unit = "rad", diagnose_branch_shift_conventions = true, maxlines = 2)) === nothing
      @test occursin("Branch-shift convention scan", String(take!(diag_io)))

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
        matpower_shift_unit = "rad",
        matpower_shift_sign = -1.0,
        matpower_ratio = "reciprocal",
        log_effective_config = true,
        effective_config = cfg,
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

  @testset "MATPOWER SHIFT convention options" begin
    bus = [1.0 3.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 110.0 1.0 1.1 0.9;
           2.0 1.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 110.0 1.0 1.1 0.9]
    branch_rad = [1.0 2.0 0.01 0.10 0.02 0.0 0.0 0.0 1.0 (pi / 6) 1.0 -360.0 360.0]
    branch_deg = copy(branch_rad)
    branch_deg[1, 10] = -30.0
    y_rad = Sparlectra.MatpowerIO.build_ybus_matpower(bus, branch_rad, 100.0; matpower_shift_unit = "rad", matpower_shift_sign = -1.0)
    y_deg = Sparlectra.MatpowerIO.build_ybus_matpower(bus, branch_deg, 100.0)
    @test isapprox(y_rad, y_deg; atol = 1e-12, rtol = 0.0)
  end

  @testset "MATPOWER TAP ratio convention options" begin
    bus = [1.0 3.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 110.0 1.0 1.1 0.9;
           2.0 1.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 110.0 1.0 1.1 0.9]
    branch = [1.0 2.0 0.01 0.10 0.02 0.0 0.0 0.0 2.0 0.0 1.0 0.0 0.0]
    branch_recip = copy(branch)
    branch_recip[1, 9] = 0.5
    y_recip_option = Sparlectra.MatpowerIO.build_ybus_matpower(bus, branch, 100.0; matpower_ratio = "reciprocal")
    y_recip_data = Sparlectra.MatpowerIO.build_ybus_matpower(bus, branch_recip, 100.0)
    @test isapprox(y_recip_option, y_recip_data; atol = 1e-12, rtol = 0.0)

    mpc = (;
      name = "ratio_convention",
      baseMVA = 100.0,
      bus = bus,
      gen = [1.0 0.0 0.0 10.0 -10.0 1.0 100.0 1.0 10.0 0.0 zeros(1, 11)...],
      branch = branch,
      gencost = nothing,
      bus_name = nothing,
    )
    net_normal = Sparlectra.createNetFromMatPowerCase(mpc = mpc, matpower_ratio = "normal")
    net_recip = Sparlectra.createNetFromMatPowerCase(mpc = mpc, matpower_ratio = "reciprocal")
    @test net_normal.branchVec[1].ratio == 2.0
    @test net_recip.branchVec[1].ratio == 0.5
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

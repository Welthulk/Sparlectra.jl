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

# file: test/test_matpower_example.jl
# purpose: tests the MATPOWER runner: the standard path with warmup runs, the
#          compact summary, the timing coverage of the performance profile,
#          the benchmark mode with and without BenchmarkTools, and the
#          runtime thread status
using Sparlectra
using LinearAlgebra
using Test
# the benchmark mode of the runner lives in the package extension
# SparlectraBenchmarkToolsExt; loading BenchmarkTools here loads it (under
# Pkg.test the test target brings the package, otherwise the test runner
# stacks a small environment with it onto the load path)
using BenchmarkTools

function run_matpower_example_tests()
  # the shipped case, so the runner never fetches anything
  case9 = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "sp_case9.m"))

  # The only run of the standard (non-benchmark) path of run_matpower_case,
  # with the warmup loop switched on: two extra solves before the
  # representative one.
  @testset "MATPOWER runner standard path with warmup runs" begin (function ()
    test_cfg = test_scratch_path(".yaml")
    write(test_cfg, """
config_version: 1
scope: general
benchmark:
  enabled: false
performance:
  enabled: true
  representative_warmup_runs: 2
power_flow:
  flatstart: true
  start_mode:
    angle_mode: dc
    voltage_mode: pv_gen_vg
""")
    captured = Sparlectra.run_with_output_capture() do
      Sparlectra.run_matpower_case(; casefile = case9, config_file = test_cfg)
    end
    result = captured.result
    @test isfile(result.logfile)
    @test result.status_ref !== nothing
    @test occursin("summary method=", read(result.logfile, String))
    profile = result.performance_profile
    @test profile[:representative_warmup_runs] == 2
    # set only inside the warmup loop
    @test haskey(profile, :last_warmup_elapsed_s)
    @test profile[:timing_mode] === :warm_steady_state
  end)() end

  @testset "Compact summary reports the canonical status fields" begin (function ()
    base = (method = :rectangular, final_converged = false, iterations = 30, elapsed_s = 1.25)
    # status fields, fragments the summary carries, fragments it must not carry
    rows = [
      # outcome, solution and limit check come from the canonical fields
      ((outcome = :not_converged, numerical_converged = false, solution_available = false, limit_validation_status = :skip,
        final_mismatch = 2.575667589, reason_text = "NR mismatch did not converge"),
       ["outcome=not_converged", "numerical_solution=FAIL", "solution_available=false", "limit_validation=SKIP", "reason=\"NR mismatch did not converge\""],
       String[]),
      # a converged solve that fails the limit check is not a numerical failure
      ((outcome = :converged_limits_failed, numerical_converged = true, solution_available = true, limit_validation_status = :fail,
        final_mismatch = 1.2e-7, reason_text = "remaining PV Q-limit violations", iterations = 8, elapsed_s = 0.5),
       ["outcome=converged_limits_failed", "numerical_solution=OK", "solution_available=true", "limit_validation=FAIL"],
       String[]),
      # the Q-limit counters are passed through
      ((outcome = :not_converged, numerical_converged = false, solution_available = false, limit_validation_status = :skip,
        final_mismatch = 1.0, iterations = 80, pv2pq_events = 0, pv2pq_buses = 0, qlimit_active_set_changes = 76, qlimit_reenable_events = 63,
        reason_text = "NR mismatch did not converge; Q-limit active set changed repeatedly"),
       ["pv2pq_events=0", "qlimit_active_set_changes=76", "qlimit_reenable_events=63"],
       String[]),
      # representative and solver time are not swapped, the benchmark median is in ms
      ((outcome = :converged, numerical_converged = true, solution_available = true, limit_validation_status = :ok, final_converged = true,
        final_mismatch = 1e-8, iterations = 5, elapsed_s = 14.522494, solver_elapsed_s = 8.948731, benchmark_median_s = 0.018227,
        pv2pq_events = 0, pv2pq_buses = 0, reason_text = "none"),
       ["representative_time=14.522494 s", "solver_time=8.948731 s", "benchmark_median=18.227 ms"],
       ["solver_time=14.522494 s"]),
    ]
    for (status, present, absent) in rows
      summary = Sparlectra._compact_run_summary(merge(base, status))
      for fragment in present
        @test occursin(fragment, summary)
      end
      for fragment in absent
        @test !occursin(fragment, summary)
      end
    end
  end)() end

  @testset "Q-limit stop reasons at the iteration limit have their own text" begin (function ()
    # the three reasons of a run at max_iter (they replace
    # :nr_mismatch_not_converged_active_set_unstable) have their own text,
    # not the underscore fallback
    for reason in Sparlectra.QLIMIT_MAX_ITER_STOP_REASONS
      @test Sparlectra._rectangular_rejection_reason_text(reason) != replace(String(reason), "_" => " ")
    end
  end)() end

  @testset "MATPOWER performance profile timing coverage" begin (function ()
    profile = Dict{Symbol,Any}(
      :enabled => true,
      :level => :full,
      :show_allocations => true,
      :show_iteration_table => true,
      :representative_elapsed_s => 10.0,
      :events => [(label = :benchmark_rectangular, seconds = 0.025)],
      :timings => Dict(
        :network_construction => (calls = 1, elapsed_s = 2.0, bytes = 100),
        :start_projection => (calls = 1, elapsed_s = 1.0, bytes = 200),
        :solver_total => (calls = 1, elapsed_s = 6.0, bytes = 300),
        :newton_step_linear_solve => (calls = 4, elapsed_s = 0.5, bytes = 400),
      ),
      :iterations => [(iter = 1, mismatch = 1.0, elapsed_s = 0.01)],
      :timing_mode => :cold_representative,
      :representative_warmup_runs => 0,
    )
    # the number on the one line that starts with the label
    reported(txt, label) = parse(Float64, first(split(strip(split(only(filter(l -> startswith(l, label), split(txt, '\n'))), ':'; limit = 2)[2]))))

    io_full = IOBuffer()
    Sparlectra.print_performance_profile(io_full, profile; level = :full, max_rows = 25)
    txt_full = String(take!(io_full))
    @test reported(txt_full, "Representative wall time") == 10.0
    # top level: network construction 2 + start projection 1 + solver 6;
    # the nested linear solve and the benchmark event stay out
    @test reported(txt_full, "Top-level measured time") == 9.0
    @test reported(txt_full, "Sum of all recorded phase timings") == 9.5
    @test reported(txt_full, "Benchmark events total") == 0.025
    # the full level carries the phase table and the iteration block
    @test occursin("network_construction", txt_full)
    @test occursin("Iteration diagnostics", txt_full)

    # compact logging prints the coverage alone: no phase table, no
    # iteration block, the same numbers
    compact_profile = merge(profile, Dict{Symbol,Any}(:compact_logging => true))
    output = Sparlectra.run_with_output_capture() do
      Sparlectra.emit_performance_summary(compact_profile; print_to_console = true, write_to_logfile = false, console_level = :compact)
    end
    txt_compact = output.stdout
    @test !occursin("network_construction", txt_compact)
    @test !occursin("Iteration diagnostics", txt_compact)
    @test reported(txt_compact, "Top-level measured time") == 9.0
  end)() end

  # Without BenchmarkTools the benchmark mode stops before the case is read,
  # with the way to get it. A process of its own, in the library environment
  # alone: in this process the extension is loaded already. With the
  # extension (this process) the benchmark block is written once, one row
  # per method.
  @testset "benchmark mode with and without BenchmarkTools" begin (function ()
    cfg_on = test_scratch_path(".yaml")
    write(cfg_on, "config_version: 1\nscope: general\nbenchmark:\n  enabled: true\n  methods: [rectangular]\n  seconds: 0.05\n  samples: 2\n")
    script = test_scratch_path(".jl")
    write(script, """
    using Sparlectra
    println("extension loaded: ", Base.get_extension(Sparlectra, :SparlectraBenchmarkToolsExt) !== nothing)
    try
      Sparlectra.run_matpower_case(; casefile = "no_such_case_is_read.m", config_file = $(repr(cfg_on)))
      println("no error")
    catch err
      # the expected failure is the ArgumentError naming BenchmarkTools;
      # anything else is printed too and fails the assertions below
      println("error: ", sprint(showerror, err))
    end
    """)
    library = normpath(joinpath(@__DIR__, ".."))
    # stderr goes into the same text, so a warning of the child process
    # would show in the assertions instead of in the test log
    buffer = IOBuffer()
    run(pipeline(Cmd(`$(Base.julia_cmd()) --startup-file=no --project=$(library) $(script)`; ignorestatus = true); stdout = buffer, stderr = buffer))
    out = String(take!(buffer))
    @test !occursin("Warning", out)
    @test occursin("extension loaded: false", out)
    @test occursin("benchmark.enabled needs the package BenchmarkTools", out)
    @test occursin("using BenchmarkTools", out)
    # the check comes first: the case of the run is never looked for
    @test !occursin("no_such_case_is_read", out)

    @test Base.get_extension(Sparlectra, :SparlectraBenchmarkToolsExt) !== nothing
    captured = Sparlectra.run_with_output_capture() do
      Sparlectra.run_matpower_case(; casefile = case9, config_file = cfg_on)
    end
    lines = split(captured.stdout, '\n')
    @test count(l -> occursin("Benchmark results", l), lines) == 1
    @test count(l -> startswith(l, "benchmark method=rectangular"), lines) == 1
  end)() end

  @testset "Runtime thread status reporting" begin (function ()
    # defaults: nothing requested, the parallel sites are on with one task
    # per Julia thread
    status_default = Sparlectra.runtime_thread_status(Sparlectra.RuntimeConfig(Dict{String,Any}()))
    @test status_default.parallel_enabled === true
    @test status_default.parallel_max_tasks == Threads.nthreads()
    @test occursin("parallel: enabled=true max_tasks=$(Threads.nthreads())", sprint(Sparlectra.print_runtime_thread_config, status_default))

    cfg_keep = Sparlectra.RuntimeConfig(; print_thread_config = true, julia_threads = "keep", blas_threads = "keep")
    status_keep = Sparlectra.runtime_thread_status(cfg_keep)
    @test status_keep.julia_applied === true
    @test status_keep.julia_startup_required === false

    # a Julia thread count can only be set at process start
    requested = Threads.nthreads() == 1 ? 2 : 1
    cfg_julia = Sparlectra.RuntimeConfig(; print_thread_config = true, julia_threads = string(requested), blas_threads = "keep")
    status_julia = Sparlectra.runtime_thread_status(cfg_julia)
    @test status_julia.requested_julia_threads_resolved == requested
    @test status_julia.julia_applied === false
    @test status_julia.julia_startup_required === true
    @test occursin("julia --threads=$(requested)", sprint(Sparlectra.print_runtime_thread_config, status_julia))

    # a BLAS thread count applies at run time
    cfg_blas = Sparlectra.RuntimeConfig(; print_thread_config = true, julia_threads = "keep", blas_threads = string(LinearAlgebra.BLAS.get_num_threads()))
    status_blas = Sparlectra.runtime_thread_status(cfg_blas)
    @test status_blas.blas_applied === true
  end)() end
end

# Executed from test/runtests.jl to avoid duplicate execution and
# Julia 1.12 world-age warnings when this file is included into Main.

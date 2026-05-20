using Sparlectra
using Test

function run_matpower_example_tests()
  @testset "Central Sparlectra configuration" begin
    cfg = Sparlectra.load_sparlectra_config(; reload = true)
    @test cfg.powerflow.method === :rectangular
    @test cfg.powerflow.sparse === true

    bad_method_cfg = tempname() * ".yaml"
    write(bad_method_cfg, "power_flow:\n  method: polar\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(bad_method_cfg; reload = true)

    bad_sparse_cfg = tempname() * ".yaml"
    write(bad_sparse_cfg, "power_flow:\n  sparse: false\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(bad_sparse_cfg; reload = true)

    bad_unknown_cfg = tempname() * ".yaml"
    write(bad_unknown_cfg, "power_flow:\n  typo_tol: 1.0e-4\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(bad_unknown_cfg; reload = true)

    removed_key_cfg = tempname() * ".yaml"
    write(removed_key_cfg, "matpower_import:\n  benchmark: true\n")
    err_removed = try
      Sparlectra.load_sparlectra_config(removed_key_cfg; reload = true)
      nothing
    catch err
      err
    end
    @test err_removed isa ArgumentError
    @test occursin("matpower_import.benchmark", sprint(showerror, err_removed))
    @test occursin("benchmark.enabled", sprint(showerror, err_removed))

    bench_cfg = tempname() * ".yaml"
    write(bench_cfg, "benchmark:\n  enabled: false\n  methods: [rectangular]\n  seconds: 0.1\n  samples: 2\n  show_once: true\n")
    cfg_bench = Sparlectra.load_sparlectra_config(bench_cfg; reload = true)
    @test cfg_bench.benchmark.enabled === false
    @test cfg_bench.benchmark.methods == [:rectangular]

    bad_bench_method_cfg = tempname() * ".yaml"
    write(bad_bench_method_cfg, "benchmark:\n  methods: [polar]\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(bad_bench_method_cfg; reload = true)

    startmode_cfg = tempname() * ".yaml"
    write(startmode_cfg, """
power_flow:
  start_mode:
    angle_mode: dc
    voltage_mode: pv_gen_vg
""")
    cfg_startmode = Sparlectra.load_sparlectra_config(startmode_cfg; reload = true)
    @test cfg_startmode.powerflow.start_mode.angle_mode === :dc
    @test cfg_startmode.powerflow.start_mode.voltage_mode === :pv_gen_vg

    for mode in ("classic", "dc", "bus_va_blend", "matpower_va")
      cfg_mode = tempname() * ".yaml"
      write(cfg_mode, "power_flow:\n  start_mode:\n    angle_mode: $(mode)\n")
      cfg_loaded = Sparlectra.load_sparlectra_config(cfg_mode; reload = true)
      @test cfg_loaded.powerflow.start_mode.angle_mode === Symbol(mode)
    end

    for mode in ("classic", "pv_gen_vg", "pv_bus_vm", "all_bus_vm", "bus_vm_va_blend")
      cfg_mode = tempname() * ".yaml"
      write(cfg_mode, "power_flow:\n  start_mode:\n    voltage_mode: $(mode)\n")
      cfg_loaded = Sparlectra.load_sparlectra_config(cfg_mode; reload = true)
      @test cfg_loaded.powerflow.start_mode.voltage_mode === Symbol(mode)
    end

    bad_startmode_cfg = tempname() * ".yaml"
    write(bad_startmode_cfg, "power_flow:\n  start_mode:\n    voltage_mode: nonsense\n")
    err_startmode = try
      Sparlectra.load_sparlectra_config(bad_startmode_cfg; reload = true)
      nothing
    catch err
      err
    end
    @test err_startmode isa ArgumentError
    @test occursin("power_flow.start_mode.voltage_mode", sprint(showerror, err_startmode))
    @test occursin("pv_gen_vg", sprint(showerror, err_startmode))

    cfg_roundtrip = tempname() * ".yaml"
    write(cfg_roundtrip, """
power_flow:
  qlimits:
    enabled: false
    hysteresis_pu: 0.123
    cooldown_iters: 7
matpower_import:
  enable_pq_gen_controllers: false
  ratio: reciprocal
  bus_shunt_model: voltage_dependent_injection
state_estimation:
  method: wls
""")
    cfg_loaded = Sparlectra.load_sparlectra_config(cfg_roundtrip; reload = true)
    @test cfg_loaded.powerflow.qlimits.ignore_q_limits === true
    @test cfg_loaded.powerflow.qlimits.hysteresis_pu == 0.123
    @test cfg_loaded.powerflow.qlimits.cooldown_iters == 7
    @test cfg_loaded.matpower.enable_pq_gen_controllers === false
    @test cfg_loaded.matpower.ratio === :reciprocal
    @test cfg_loaded.matpower.bus_shunt_model === :voltage_dependent_injection
    @test cfg_loaded.state_estimation.method === :wls

    for mode in ("gen_vg", "bus_vm", "auto", "strict_check")
      cfg_mode = tempname() * ".yaml"
      write(cfg_mode, "matpower_import:\n  pv_voltage_source: $(mode)\n")
      @test Sparlectra.load_sparlectra_config(cfg_mode; reload = true).matpower.pv_voltage_source === Symbol(mode)
    end
    for mode in ("bus_vm", "gen_vg", "imported_setpoint", "hybrid")
      cfg_mode = tempname() * ".yaml"
      write(cfg_mode, "matpower_import:\n  compare_voltage_reference: $(mode)\n")
      @test Sparlectra.load_sparlectra_config(cfg_mode; reload = true).matpower.compare_voltage_reference === Symbol(mode)
    end
    for mode in ("off", "on", "auto")
      cfg_mode = tempname() * ".yaml"
      write(cfg_mode, "power_flow:\n  rectangular_preallocate_workspace: $(mode)\n")
      @test Sparlectra.load_sparlectra_config(cfg_mode; reload = true).powerflow.rectangular_preallocate_workspace === Symbol(mode)
    end
    bad_ws_cfg = tempname() * ".yaml"
    write(bad_ws_cfg, "power_flow:\n  rectangular_preallocate_workspace: invalid\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(bad_ws_cfg; reload = true)
  end

  @testset "MATPOWER example cleanup guard" begin
    source = read(joinpath(@__DIR__, "..", "examples", "matpower_import.jl"), String)
    @test !occursin("load_yaml_dict", source)
    @test !occursin("_normalize_matpower_example_config", source)
    @test !occursin("_copy_nested_yaml!", source)
    @test !occursin("_nested_yaml_get", source)
    @test !occursin("Dict{String,Any}", source)
    @test occursin("Sparlectra.run_matpower_case", source)
    @test occursin("Base.invokelatest(getfield(@__MODULE__, :main))", source)
    @test occursin("Base.invokelatest(", source)

    example_path = joinpath(@__DIR__, "..", "examples", "matpower_import.jl")
    old_no_main = get(ENV, "SPARLECTRA_MATPOWER_IMPORT_NO_MAIN", nothing)
    ENV["SPARLECTRA_MATPOWER_IMPORT_NO_MAIN"] = "1"
    mod = Module(:MatpowerImportExampleSmoke)
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
  end

  @testset "MATPOWER runner operational output" begin
    test_cfg = tempname() * ".yaml"
    write(test_cfg, """
benchmark:
  enabled: false
power_flow:
  flatstart: true
  start_mode:
    angle_mode: dc
    voltage_mode: pv_gen_vg
matpower_import:
  case: case14.m
""")
    try
      result = Sparlectra.run_matpower_case(; casefile = "case14.m", config_file = test_cfg)
      @test isfile(result.logfile)
      @test result.status_ref !== nothing
      logtxt = read(result.logfile, String)
      @test occursin("Sparlectra version", logtxt)
      @test occursin("config default", logtxt)
      @test occursin("config user", logtxt)
      @test occursin("casefile", logtxt)
      @test occursin("logfile", logtxt)
      @test occursin("performance", logtxt)
      @test occursin("summary method=", logtxt)
      @test occursin("summary", logtxt)
    catch err
      # Network fetch can fail in CI/offline runs depending on HTTP stack versions.
      et = string(typeof(err))
      em = sprint(showerror, err)
      @test occursin("RequestError", et) || occursin("FieldError", et) || occursin("RequestError", em)
    end
  end

  @testset "Compact summary uses canonical outcome fields" begin
    summary_not_conv = Sparlectra._compact_run_summary((
      method = :rectangular,
      outcome = :not_converged,
      numerical_converged = false,
      solution_available = false,
      limit_validation_status = :skip,
      final_converged = false,
      final_mismatch = 2.575667589,
      reason_text = "NR mismatch did not converge",
      iterations = 30,
      elapsed_s = 1.25,
    ))
    @test occursin("outcome=not_converged", summary_not_conv)
    @test occursin("numerical_solution=FAIL", summary_not_conv)
    @test occursin("solution_available=false", summary_not_conv)
    @test occursin("limit_validation=SKIP", summary_not_conv)

    summary_limits = Sparlectra._compact_run_summary((
      method = :rectangular,
      outcome = :converged_limits_failed,
      numerical_converged = true,
      solution_available = true,
      limit_validation_status = :fail,
      final_converged = false,
      final_mismatch = 1.2e-7,
      reason_text = "remaining PV Q-limit violations",
      iterations = 8,
      elapsed_s = 0.5,
    ))
    @test occursin("outcome=converged_limits_failed", summary_limits)
    @test occursin("numerical_solution=OK", summary_limits)
    @test occursin("solution_available=true", summary_limits)
    @test occursin("limit_validation=FAIL", summary_limits)
  end

  @testset "MATPOWER performance profile output coverage" begin
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
    )
    io_full = IOBuffer()
    Sparlectra.print_performance_profile(io_full, profile; level = :full, max_rows = 25)
    txt_full = String(take!(io_full))
    @test occursin("network_construction", txt_full)
    @test occursin("benchmark_rectangular", txt_full)
    @test occursin("Timing coverage", txt_full)
    @test occursin("Representative wall time", txt_full)
    @test occursin("Top-level measured time", txt_full)
    @test occursin("Sum of all recorded phase timings", txt_full)
    @test occursin("Benchmark events total", txt_full)
    @test occursin("Coverage status", txt_full)
    @test occursin("Nested timing note", txt_full)
    @test occursin("10.0", txt_full)
    @test occursin("9.0", txt_full)
    @test !occursin("9.025", txt_full)
    @test occursin("Bytes", txt_full)
    @test occursin("Bytes/Call", txt_full)
    @test occursin("Iteration diagnostics", txt_full)

    io_compact = IOBuffer()
    Sparlectra.print_performance_profile(io_compact, profile; level = :compact, max_rows = 25)
    txt_compact = String(take!(io_compact))
    @test occursin("network_construction", txt_compact)
    @test occursin("benchmark_rectangular", txt_compact)
    @test occursin("Timing coverage", txt_compact)
  end

  @testset "MATPOWER benchmark output routing" begin
    test_cfg = tempname() * ".yaml"
    write(test_cfg, """
benchmark:
  enabled: true
  methods: [rectangular]
  seconds: 0.05
  samples: 3
output:
  console_summary: true
  console_diagnostics: compact
  console_q_limit_events: summary
  logfile_diagnostics: full
performance:
  enabled: true
  write_to_logfile: true
diagnostics:
  log_effective_config: true
matpower_import:
  case: case14.m
""")
    result = try
      captured = Sparlectra.run_with_output_capture() do
        Sparlectra.run_matpower_case(; casefile = "case14.m", config_file = test_cfg)
      end
      captured
    catch err
      et = string(typeof(err))
      em = sprint(showerror, err)
      if occursin("RequestError", et) || occursin("FieldError", et) || occursin("RequestError", em)
        @test true
        nothing
      else
        rethrow(err)
      end
    end
    isnothing(result) && return
    txt = result.stdout
    @test count(occursin("Benchmark results", line) for line in split(txt, '\n')) == 1
    @test count(occursin("benchmark method=rectangular", line) for line in split(txt, '\n')) == 1
    @test occursin("representative=", txt)
    @test count(occursin("Q-Limit Active-Set Summary", line) for line in split(txt, '\n')) <= 1
    logtxt = read(result.result.logfile, String)
    @test occursin("Effective Sparlectra Configuration", logtxt)
    @test occursin("summary method=", logtxt)
    @test occursin("representative_time=", logtxt)
    @test occursin("solver_time=", logtxt)
    @test occursin("benchmark_median=", logtxt)
    @test occursin("Top-level measured time", logtxt)
    @test occursin("Top-level excluding output", logtxt)
    @test occursin("Result output time", logtxt)
    @test !occursin("q_limit_active_set=OK final_converged=false", logtxt)
    @test occursin("Representative solve diagnostics", logtxt) || occursin("Solve diagnostics (compact)", logtxt)
  end

  @testset "MATPOWER compact timing summary uses solver_total semantics" begin
    summary = Sparlectra._compact_run_summary((
      method = :rectangular,
      outcome = :converged,
      numerical_converged = true,
      solution_available = true,
      limit_validation_status = :ok,
      final_converged = true,
      final_mismatch = 1e-8,
      iterations = 5,
      elapsed_s = 14.522494,
      solver_elapsed_s = 8.948731,
      benchmark_median_s = 0.018227,
      pv2pq_events = 0,
      pv2pq_buses = 0,
      reason_text = "none",
    ))
    @test occursin("representative_time=14.522494 s", summary)
    @test occursin("solver_time=8.948731 s", summary)
    @test occursin("benchmark_median=18.227 ms", summary)
    @test !occursin("solver_time=14.522494 s", summary)
  end



  @testset "Runtime config accepts numeric thread values" begin
    cfg_file = tempname() * ".yaml"
    write(cfg_file, "runtime:
  julia_threads: 4
  blas_threads: 16
")
    cfg_runtime = Sparlectra.load_sparlectra_config(cfg_file; reload = true)
    @test cfg_runtime.runtime.julia_threads == "4"
    @test cfg_runtime.runtime.blas_threads == "16"
  end

  @testset "Runtime thread status reporting" begin
    cfg_keep = Sparlectra.RuntimeConfig(; print_thread_config = true, julia_threads = "keep", blas_threads = "keep")
    status_keep = Sparlectra.runtime_thread_status(cfg_keep)
    @test status_keep.julia_applied === true
    @test status_keep.julia_startup_required === false

    requested = Threads.nthreads() == 1 ? 2 : 1
    cfg_julia = Sparlectra.RuntimeConfig(; print_thread_config = true, julia_threads = string(requested), blas_threads = "keep")
    status_julia = Sparlectra.runtime_thread_status(cfg_julia)
    @test status_julia.requested_julia_threads_resolved == requested
    @test status_julia.julia_applied === false
    @test status_julia.julia_startup_required === true
    io = IOBuffer()
    Sparlectra.print_runtime_thread_config(io, status_julia)
    txt = String(take!(io))
    @test occursin("Request not applied: Julia threads must be set at process startup.", txt)
    @test occursin("julia --threads=$(requested)", txt)

    cfg_blas = Sparlectra.RuntimeConfig(; print_thread_config = true, julia_threads = "keep", blas_threads = string(LinearAlgebra.BLAS.get_num_threads()))
    status_blas = Sparlectra.runtime_thread_status(cfg_blas)
    @test status_blas.blas_applied === true
  end

  @testset "Compact performance summary output" begin
    profile = Dict{Symbol,Any}(
      :compact_logging => true,
      :representative_elapsed_s => 10.0,
      :events => [(label = :benchmark_rectangular, seconds = 0.02)],
      :timings => Dict(
        :network_construction => (calls = 1, elapsed_s = 2.5, bytes = 10),
        :result_output => (calls = 1, elapsed_s = 1.0, bytes = 10),
        :solver_total => (calls = 1, elapsed_s = 5.0, bytes = 10),
      ),
      :iterations => NamedTuple[],
    )
    output = Sparlectra.run_with_output_capture() do
      Sparlectra.emit_performance_summary(profile; print_to_console = true, write_to_logfile = false, console_level = :compact)
    end
    txt = output.stdout
    @test occursin("Performance Summary (compact)", txt)
    @test occursin("Coverage status", txt)
    @test !occursin("Phase", txt)
    @test !occursin("Iteration diagnostics", txt)
  end
end

# Executed from test/runtests.jl to avoid duplicate execution and
# Julia 1.12 world-age warnings when this file is included into Main.

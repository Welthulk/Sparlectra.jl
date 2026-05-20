"""
    new_performance_profile(cfg::PerformanceConfig)

Create a lightweight mutable performance profile dictionary used by runner APIs.
"""
function new_performance_profile(cfg::PerformanceConfig)
  return Dict{Symbol,Any}(
    :enabled => cfg.enabled,
    :level => cfg.level,
    :show_allocations => cfg.show_allocations,
    :show_iteration_table => cfg.show_iteration_table,
    :compact_logging => cfg.compact_logging,
    :max_diagnostic_rows => cfg.max_diagnostic_rows,
    :started_at => time(),
    :events => NamedTuple[],
    :timings => Dict{Symbol,Any}(),
    :iterations => NamedTuple[],
  )
end

function _record_perf!(profile, label::Symbol, seconds::Real)
  profile isa AbstractDict || return profile
  events = get!(profile, :events, NamedTuple[])
  push!(events, (label = label, seconds = Float64(seconds)))
  return profile
end

_perf_timing_seconds(row)::Float64 = row isa NamedTuple && hasproperty(row, :elapsed_s) ? Float64(getproperty(row, :elapsed_s)) : 0.0
_perf_timing_bytes(row)::Int = row isa NamedTuple && hasproperty(row, :bytes) ? Int(getproperty(row, :bytes)) : 0

function _is_top_level_perf_phase(phase::Symbol)::Bool
  return phase in (:network_construction, :start_projection, :solver_total, :result_output)
end

function _sum_perf_timings(profile; top_level_only::Bool = false)::Float64
  profile isa AbstractDict || return 0.0
  timings = get(profile, :timings, Dict{Symbol,Any}())
  total = 0.0
  for (phase, row) in timings
    top_level_only && !_is_top_level_perf_phase(phase) && continue
    total += _perf_timing_seconds(row)
  end
  return total
end

function _sum_perf_events(profile)::Float64
  profile isa AbstractDict || return 0.0
  events = get(profile, :events, NamedTuple[])
  total = 0.0
  for row in events
    total += row isa NamedTuple && hasproperty(row, :seconds) ? Float64(getproperty(row, :seconds)) : 0.0
  end
  return total
end

function _print_timing_coverage(io::IO, profile; level::Symbol = :compact)
  profile isa AbstractDict || return nothing
  representative_elapsed_s = Float64(get(profile, :representative_elapsed_s, 0.0))
  sum_top_level = _sum_perf_timings(profile; top_level_only = true)
  result_output_s = haskey(get(profile, :timings, Dict{Symbol,Any}()), :result_output) ? _perf_timing_seconds(profile[:timings][:result_output]) : 0.0
  sum_top_level_ex_output = sum_top_level - result_output_s
  sum_all = _sum_perf_timings(profile; top_level_only = false)
  benchmark_total = _sum_perf_events(profile)
  coverage = representative_elapsed_s > 0.0 ? sum_top_level / representative_elapsed_s : 0.0
  delta_s = representative_elapsed_s - sum_top_level
  status = coverage >= 0.9 ? "OK" : coverage >= 0.7 ? "PARTIAL" : "WARNING"
  status_detail = representative_elapsed_s > 0.0 ? "$(status), top-level timings explain $(round(coverage * 100.0; digits = 1))% of wall time" : "$(status), representative wall time unavailable"

  println(io)
  println(io, "Timing coverage")
  println(io, "---------------")
  timing_mode = get(profile, :timing_mode, :cold_representative)
  println(io, "Timing mode                    : ", timing_mode == :warm_steady_state ? "warm / steady-state" : "cold / representative")
  timing_mode == :cold_representative && println(io, "Timing note                    : first representative run may include JIT, sparse/BLAS init, and cache construction")
  println(io, "Representative wall time       : ", round(representative_elapsed_s; digits = 6), " s")
  println(io, "Top-level measured time        : ", round(sum_top_level; digits = 6), " s")
  println(io, "Top-level excluding output     : ", round(sum_top_level_ex_output; digits = 6), " s")
  println(io, "Result output time             : ", round(result_output_s; digits = 6), " s")
  println(io, "Sum of all recorded phase timings : ", round(sum_all; digits = 6), " s")
  println(io, "Benchmark events total         : ", round(benchmark_total; digits = 6), " s")
  if delta_s >= 0.0
    println(io, "Unaccounted / overhead time    : ", round(delta_s; digits = 6), " s")
  else
    println(io, "Nested/overlapping excess      : ", round(abs(delta_s); digits = 6), " s")
  end
  println(io, "Coverage status                : ", status_detail)
  println(io, "Nested timing note             : recorded phases include nested timings; their sum may exceed wall time")
  level === :compact && println(io, "  (compact mode: iteration diagnostics omitted)")
  println(io)
  println(io, "Rectangular workspace")
  println(io, "---------------------")
  println(io, "reuse         : ", get(profile, :rectangular_workspace_reuse, false))
  println(io, "preallocated  : ", get(profile, :rectangular_workspace_preallocated, false))
  println(io, "reason        : ", get(profile, :rectangular_workspace_reason, :unsupported_solver))
  haskey(profile, :rectangular_workspace_nbus) && println(io, "nbus          : ", get(profile, :rectangular_workspace_nbus, 0))
  haskey(profile, :rectangular_workspace_nstate) && println(io, "nstate        : ", get(profile, :rectangular_workspace_nstate, 0))
  return nothing
end

function print_performance_profile(io::IO, profile; title::AbstractString = "Performance Summary", max_rows::Int = 25, level::Union{Symbol,Nothing} = nothing)
  println(io, title)
  println(io, "-"^length(title))
  if !(profile isa AbstractDict)
    println(io, "(no recorded events)")
    return nothing
  end
  effective_level = isnothing(level) ? get(profile, :level, :compact) : level
  show_allocations = Bool(get(profile, :show_allocations, false))
  show_iteration_table = Bool(get(profile, :show_iteration_table, false))
  timings = get(profile, :timings, Dict{Symbol,Any}())
  events = get(profile, :events, NamedTuple[])
  iterations = get(profile, :iterations, NamedTuple[])

  has_timing_rows = false
  for row in values(timings)
    if row isa NamedTuple && hasproperty(row, :calls) && hasproperty(row, :elapsed_s)
      has_timing_rows = true
      break
    end
  end

  if has_timing_rows || !isempty(events)
    phase_w = 42
    println(io, rpad("Phase", phase_w), lpad("Calls", 8), lpad("Time [s]", 14), show_allocations ? lpad("Bytes", 14) * lpad("Bytes/Call", 14) : "")
    shown = 0
    for phase in sort!(collect(keys(timings)); by = string)
      shown >= max_rows && break
      row = timings[phase]
      row isa NamedTuple || continue
      hasproperty(row, :calls) || continue
      hasproperty(row, :elapsed_s) || continue
      calls = getproperty(row, :calls)
      elapsed_s = getproperty(row, :elapsed_s)
      bytes = _perf_timing_bytes(row)
      bytes_per_call = calls > 0 ? Int(fld(bytes, calls)) : 0
      bytes_str = show_allocations ? lpad(string(bytes), 14) * lpad(string(bytes_per_call), 14) : ""
      println(io, rpad(string(phase), phase_w), lpad(string(calls), 8), lpad(string(round(elapsed_s; digits = 6)), 14), bytes_str)
      shown += 1
    end
    for row in events
      shown >= max_rows && break
      phase = hasproperty(row, :label) ? getproperty(row, :label) : :event
      elapsed_s = hasproperty(row, :seconds) ? getproperty(row, :seconds) : 0.0
      bytes_str = show_allocations ? lpad("-", 14) * lpad("-", 14) : ""
      println(io, rpad(string(phase), phase_w), lpad("1", 8), lpad(string(round(elapsed_s; digits = 6)), 14), bytes_str)
      shown += 1
    end
    total_rows = length(timings) + length(events)
    if total_rows > shown
      println(io, "  ... ", total_rows - shown, " additional rows")
    end
  else
    println(io, "(no recorded events)")
  end

  if effective_level === :full && !isempty(events)
    println(io)
    println(io, "Benchmark events")
    println(io, "----------------")
    rows = min(length(events), max_rows)
    for i in 1:rows
      row = events[i]
      label = hasproperty(row, :label) ? getproperty(row, :label) : :event
      seconds = hasproperty(row, :seconds) ? getproperty(row, :seconds) : 0.0
      println(io, "  ", label, ": ", round(seconds; digits = 6), " s")
    end
  end

  if effective_level === :full && show_iteration_table && !isempty(iterations)
    println(io)
    println(io, "Iteration diagnostics")
    println(io, "---------------------")
    rows = min(length(iterations), max_rows)
    for i in 1:rows
      row = iterations[i]
      println(io, "  ", row)
    end
    if length(iterations) > rows
      println(io, "  ... ", length(iterations) - rows, " additional iterations")
    end
  end
  _print_timing_coverage(io, profile; level = effective_level)
  return nothing
end

function emit_performance_summary(profile; logfile::AbstractString = "", print_to_console::Bool = true, write_to_logfile::Bool = true, max_rows::Int = 25, console_level::Symbol = :compact, logfile_level::Union{Symbol,Nothing} = nothing)
  if print_to_console
    if console_level === :compact && Bool(get(profile, :compact_logging, false))
      println(stdout, "Performance Summary (compact)")
      println(stdout, "-----------------------------")
      _print_timing_coverage(stdout, profile; level = :compact)
    else
      print_performance_profile(stdout, profile; max_rows = max_rows, level = console_level)
    end
  end
  if write_to_logfile && !isempty(logfile)
    open(logfile, "a") do io
      print_performance_profile(io, profile; max_rows = max_rows, level = logfile_level)
      println(io)
    end
  end
  return nothing
end

function parse_runtime_threads_request(raw)::Union{Int,Nothing}
  sval = lowercase(strip(string(raw)))
  sval in ("", "keep", "default", "false", "off") && return nothing
  sval == "auto" && return Sys.CPU_THREADS
  iv = tryparse(Int, sval)
  isnothing(iv) && return nothing
  return iv > 0 ? iv : nothing
end

function runtime_thread_status(cfg::RuntimeConfig)
  requested_julia = parse_runtime_threads_request(cfg.julia_threads)
  active_julia = Threads.nthreads()
  julia_applied = isnothing(requested_julia) || requested_julia == active_julia
  julia_startup_required = !isnothing(requested_julia) && requested_julia != active_julia

  req_blas = lowercase(strip(cfg.blas_threads))
  requested_blas = parse_runtime_threads_request(cfg.blas_threads)
  before = BLAS.get_num_threads()
  if req_blas != "" && req_blas != "keep" && !isnothing(requested_blas)
    BLAS.set_num_threads(requested_blas)
  end
  active_blas = BLAS.get_num_threads()
  blas_applied = isnothing(requested_blas) || requested_blas == active_blas
  return (; 
    cpu_threads = Sys.CPU_THREADS,
    julia_threads = active_julia,
    blas_threads_before = before,
    blas_threads = active_blas,
    requested_blas_threads = cfg.blas_threads,
    requested_julia_threads = cfg.julia_threads,
    requested_julia_threads_resolved = requested_julia,
    requested_blas_threads_resolved = requested_blas,
    julia_applied = julia_applied,
    julia_startup_required = julia_startup_required,
    blas_applied = blas_applied,
  )
end

function print_runtime_thread_config(io::IO, status)
  println(io, "Runtime thread configuration")
  println(io, "CPU threads   : ", status.cpu_threads)
  println(io, "Julia threads : requested=", isnothing(status.requested_julia_threads_resolved) ? "keep" : status.requested_julia_threads_resolved, " active=", status.julia_threads, " applied=", status.julia_applied, " startup_required=", status.julia_startup_required)
  println(io, "BLAS threads  : requested=", isnothing(status.requested_blas_threads_resolved) ? "keep" : status.requested_blas_threads_resolved, " active=", status.blas_threads, " applied=", status.blas_applied)
  if status.julia_startup_required
    println(io, "Request not applied: Julia threads must be set at process startup.")
    println(io, "Start with:")
    println(io, "  julia --threads=", status.requested_julia_threads_resolved, " --project=. examples/matpower_import.jl")
  end
  return nothing
end

function matpower_run_logfile_path(casefile::AbstractString, cfg::OutputConfig)
  log_dir = joinpath(pkgdir(Sparlectra), "examples", "_out")
  mkpath(log_dir)
  stem = splitext(basename(casefile))[1]
  return joinpath(log_dir, "run_$(stem)_$(Dates.format(now(), "yyyymmdd_HHMMSS")).log")
end

function print_matpower_runner_header(io::IO; default_config::AbstractString, user_config::AbstractString, casefile::AbstractString, logfile::AbstractString, methods::AbstractVector{Symbol}, benchmark::BenchmarkConfig, performance::PerformanceConfig)
  println(io, "-"^94)
  println(io, "Sparlectra version: ", pkgversion(Sparlectra))
  println(io, "config default    : ", default_config)
  println(io, "config user       : ", isempty(user_config) ? "not found" : user_config)
  println(io, "casefile          : ", casefile)
  println(io, "logfile           : ", logfile)
  println(io, "method(s)         : ", join(string.(methods), ", "))
  println(io, "benchmark         : ", benchmark.enabled ? "enabled" : "disabled")
  println(io, "performance       : ", performance.enabled ? "enabled" : "disabled", ", level=", performance.level)
  println(io, "-"^94)
  return nothing
end

function run_with_output_capture(f::Function; capture_stdout::Bool = true, capture_stderr::Bool = true)
  stdout_text = ""
  stderr_text = ""
  result = nothing
  mktemp() do out_path, out_io
    mktemp() do err_path, err_io
      out_target = capture_stdout ? out_io : stdout
      err_target = capture_stderr ? err_io : stderr
      result = redirect_stdio(stdout = out_target, stderr = err_target) do
        f()
      end
      flush(out_io)
      flush(err_io)
      capture_stdout && (stdout_text = read(out_path, String))
      capture_stderr && (stderr_text = read(err_path, String))
    end
  end
  return (result = result, stdout = stdout_text, stderr = stderr_text)
end

run_silent_for_benchmark(f::Function) = run_with_output_capture(() -> Base.invokelatest(f); capture_stdout = true, capture_stderr = true).result

function append_captured_output_to_logfile(logfile::AbstractString, captured; section_title::AbstractString = "")
  chunks = String[]
  !isempty(captured.stdout) && push!(chunks, captured.stdout)
  !isempty(captured.stderr) && push!(chunks, captured.stderr)
  isempty(chunks) && return nothing
  open(logfile, "a") do io
    !isempty(section_title) && println(io, section_title)
    for chunk in chunks
      print(io, chunk)
      endswith(chunk, '\n') || println(io)
    end
    println(io)
  end
  return nothing
end

function _compact_run_summary(status)::String
  isnothing(status) && return "summary status=unavailable"
  representative_elapsed_s = hasproperty(status, :elapsed_s) ? getproperty(status, :elapsed_s) : 0.0
  solver_elapsed_s = hasproperty(status, :solver_elapsed_s) ? getproperty(status, :solver_elapsed_s) : nothing
  benchmark_median_s = hasproperty(status, :benchmark_median_s) ? getproperty(status, :benchmark_median_s) : nothing
  method = hasproperty(status, :method) ? getproperty(status, :method) : :rectangular
  outcome = hasproperty(status, :outcome) ? getproperty(status, :outcome) : :not_converged
  numerical_solution = hasproperty(status, :numerical_converged) ? (getproperty(status, :numerical_converged) ? "OK" : "FAIL") : "FAIL"
  solution_available = hasproperty(status, :solution_available) ? Bool(getproperty(status, :solution_available)) : (numerical_solution == "OK")
  limit_validation = hasproperty(status, :limit_validation_status) ? uppercase(String(getproperty(status, :limit_validation_status))) : (hasproperty(status, :q_limit_active_set_ok) ? (getproperty(status, :q_limit_active_set_ok) ? "OK" : "FAIL") : "SKIP")
  final_converged = hasproperty(status, :final_converged) ? getproperty(status, :final_converged) : true
  final_mismatch = hasproperty(status, :final_mismatch) ? getproperty(status, :final_mismatch) : NaN
  reason_text = hasproperty(status, :reason_text) ? getproperty(status, :reason_text) : "none"
  iterations = hasproperty(status, :iterations) ? getproperty(status, :iterations) : -1
  pv2pq_events = hasproperty(status, :pv2pq_events) ? getproperty(status, :pv2pq_events) : 0
  pv2pq_buses = hasproperty(status, :pv2pq_buses) ? getproperty(status, :pv2pq_buses) : 0
  benchmark_fragment = isnothing(benchmark_median_s) ? "" : " benchmark_median=$(round(benchmark_median_s * 1000.0; digits = 6)) ms"
  solver_fragment = isnothing(solver_elapsed_s) ? " solver_time=unavailable" : " solver_time=$(round(Float64(solver_elapsed_s); digits = 6)) s"
  result_output_s = hasproperty(status, :result_output_s) ? getproperty(status, :result_output_s) : nothing
  result_output_fragment = isnothing(result_output_s) ? "" : " result_output_time=$(round(Float64(result_output_s); digits = 6)) s"
  return "summary method=$(method) outcome=$(outcome) numerical_solution=$(numerical_solution) solution_available=$(solution_available) limit_validation=$(limit_validation) final_converged=$(final_converged) final_mismatch=$(round(final_mismatch; digits = 9)) iterations=$(iterations) representative_time=$(round(representative_elapsed_s; digits = 6)) s$(solver_fragment)$(result_output_fragment)$(benchmark_fragment) pv2pq_events=$(pv2pq_events) pv2pq_buses=$(pv2pq_buses) reason=\"$(reason_text)\""
end

function _run_matpower_single(local_case::AbstractString, cfg::SparlectraConfig, profile, status_ref)
  t = @elapsed run_acpflow(; casefile = local_case, config = cfg, performance_profile = profile, status_ref = status_ref, show_compact_result = false)
  status = status_ref[]
  if isnothing(status)
    status_ref[] = (method = cfg.powerflow.method, elapsed_s = t)
  else
    timings = get(profile, :timings, Dict{Symbol,Any}())
    result_output_s = haskey(timings, :result_output) ? _perf_timing_seconds(timings[:result_output]) : nothing
    status_ref[] = merge(status, (method = cfg.powerflow.method, elapsed_s = t, result_output_s = result_output_s))
  end
  return t
end

function _run_matpower_single_routed(local_case::AbstractString, cfg::SparlectraConfig, profile, status_ref; to_console::Bool = false, capture_for_log::Bool = false)
  if to_console && !capture_for_log
    return _run_matpower_single(local_case, cfg, profile, status_ref)
  end
  captured = run_with_output_capture(capture_stdout = true, capture_stderr = true) do
    _run_matpower_single(local_case, cfg, profile, status_ref)
  end
  return (elapsed = captured.result, captured = captured)
end

function with_powerflow_method(cfg::SparlectraConfig, method::Symbol)::SparlectraConfig
  pf = cfg.powerflow
  pf2 = PowerFlowConfig(;
    method = method,
    tol = pf.tol,
    max_iter = pf.max_iter,
    sparse = pf.sparse,
    autodamp = pf.autodamp,
    autodamp_min = pf.autodamp_min,
    start_mode = pf.start_mode,
    qlimits = pf.qlimits,
  )
  return SparlectraConfig(;
    powerflow = pf2,
    state_estimation = cfg.state_estimation,
    matpower = cfg.matpower,
    performance = cfg.performance,
    benchmark = cfg.benchmark,
    runtime = cfg.runtime,
    diagnostics = cfg.diagnostics,
    output = cfg.output,
  )
end

"""
    run_matpower_case(; casefile::AbstractString = "", config_file::AbstractString = "")

Run MATPOWER import and power-flow with central typed configuration plus
runner-level operational output and logging.
"""
function run_matpower_case(; casefile::AbstractString = "", config_file::AbstractString = "")
  if !isempty(strip(config_file))
    load_sparlectra_config!(String(config_file); reload = true)
  else
    load_sparlectra_config!()
  end
  cfg = active_sparlectra_config()
  resolved_case = isempty(strip(casefile)) ? cfg.matpower.case : String(casefile)
  isempty(strip(resolved_case)) && throw(ArgumentError("No MATPOWER case selected. Set matpower_import.case in the active configuration or pass casefile."))
  local_case = FetchMatpowerCase.ensure_casefile(resolved_case)
  logfile = matpower_run_logfile_path(local_case, cfg.output)
  print_matpower_runner_header(stdout; default_config = DEFAULT_SPARLECTRA_CONFIG_PATH, user_config = isfile(USER_SPARLECTRA_CONFIG_PATH) ? USER_SPARLECTRA_CONFIG_PATH : "", casefile = local_case, logfile = logfile, methods = cfg.benchmark.methods, benchmark = cfg.benchmark, performance = cfg.performance)
  open(logfile, "w") do io
    print_matpower_runner_header(io; default_config = DEFAULT_SPARLECTRA_CONFIG_PATH, user_config = isfile(USER_SPARLECTRA_CONFIG_PATH) ? USER_SPARLECTRA_CONFIG_PATH : "", casefile = local_case, logfile = logfile, methods = cfg.benchmark.methods, benchmark = cfg.benchmark, performance = cfg.performance)
  end
  if cfg.runtime.print_thread_config
    status = runtime_thread_status(cfg.runtime)
    print_runtime_thread_config(stdout, status)
    open(logfile, "a") do io
      print_runtime_thread_config(io, status)
    end
  end
  cfg.diagnostics.log_effective_config && open(logfile, "a") do io; print_effective_config(io, cfg); end

  profile = new_performance_profile(cfg.performance)
  status_ref = Ref{Any}(nothing)

  if cfg.benchmark.enabled
    open(logfile, "a") do io
      println(io, "Benchmark results")
    end
    println("Benchmark results")
    for method in cfg.benchmark.methods
      method_cfg = with_powerflow_method(cfg, method)
      warm_run = _run_matpower_single_routed(
        local_case,
        method_cfg,
        profile,
        status_ref;
        to_console = cfg.output.console_diagnostics === :full,
        capture_for_log = cfg.output.logfile_diagnostics !== :off || cfg.output.logfile_warnings !== :off,
      )
      warm = warm_run isa NamedTuple ? warm_run.elapsed : warm_run
      profile[:representative_elapsed_s] = warm
      if warm_run isa NamedTuple && cfg.output.logfile_diagnostics === :full
        append_captured_output_to_logfile(logfile, warm_run.captured; section_title = "Representative solve diagnostics (method=$(method))")
      end
      bench_status_ref = Ref{Any}(nothing)
      bench = BenchmarkTools.@benchmarkable run_silent_for_benchmark() do
        Base.invokelatest(
          run_acpflow;
          casefile = $local_case,
          config = $method_cfg,
          performance_profile = nothing,
          status_ref = $bench_status_ref,
          show_compact_result = false,
        )
      end
      trial = Base.invokelatest(BenchmarkTools.run, bench; samples = max(1, cfg.benchmark.samples), seconds = max(0.01, cfg.benchmark.seconds))
      med = BenchmarkTools.median(trial).time / 1e9
      mn = BenchmarkTools.minimum(trial).time / 1e9
      line = "benchmark method=$(method) representative=$(round(warm;digits=6))s median=$(round(med * 1000.0;digits=6)) ms min=$(round(mn * 1000.0;digits=6)) ms samples=$(cfg.benchmark.samples) seconds=$(cfg.benchmark.seconds)"
      println(line)
      open(logfile, "a") do io
        println(io, line)
      end
      _record_perf!(profile, Symbol("benchmark_$(method)"), med)
      if !isnothing(status_ref[])
        status_ref[] = merge(status_ref[], (benchmark_median_s = med, benchmark_min_s = mn, representative_elapsed_s = warm))
      end
    end
  else
    cold_solver_time = nothing
    warm_solver_time = nothing
    warmup_runs = max(0, cfg.performance.representative_warmup_runs)
    if warmup_runs > 0
      for _ in 1:warmup_runs
        warm_status_ref = Ref{Any}(nothing)
        _run_matpower_single(local_case, cfg, Dict{Symbol,Any}(), warm_status_ref)
      end
      profile[:timing_mode] = :warm_steady_state
    else
      profile[:timing_mode] = :cold_representative
    end
    run_single = _run_matpower_single_routed(
      local_case,
      cfg,
      profile,
      status_ref;
      to_console = cfg.output.console_diagnostics === :full,
      capture_for_log = cfg.output.logfile_diagnostics !== :off || cfg.output.logfile_warnings !== :off,
    )
    if run_single isa NamedTuple
      profile[:representative_elapsed_s] = run_single.elapsed
      st = status_ref[]
      if !isnothing(st) && hasproperty(st, :solver_elapsed_s)
        warm_solver_time = getproperty(st, :solver_elapsed_s)
      end
      if cfg.output.logfile_diagnostics === :full
        append_captured_output_to_logfile(logfile, run_single.captured; section_title = "Solve diagnostics")
      elseif cfg.output.logfile_diagnostics === :compact
        compact = _compact_run_summary(status_ref[])
        open(logfile, "a") do io
          println(io, "Solve diagnostics (compact)")
          println(io, compact)
          println(io)
        end
      end
    else
      profile[:representative_elapsed_s] = run_single
    end
    if cfg.performance.compare_cold_warm && warmup_runs > 0
      cold_status_ref = Ref{Any}(nothing)
      cold_profile = new_performance_profile(cfg.performance)
      _run_matpower_single(local_case, cfg, cold_profile, cold_status_ref)
      st = cold_status_ref[]
      if !isnothing(st) && hasproperty(st, :solver_elapsed_s)
        cold_solver_time = getproperty(st, :solver_elapsed_s)
      end
      if !(isnothing(cold_solver_time) || isnothing(warm_solver_time) || warm_solver_time == 0)
        speedup = cold_solver_time / warm_solver_time
        println("Cold run solver time : ", round(cold_solver_time; digits = 6), " s")
        println("Warm run solver time : ", round(warm_solver_time; digits = 6), " s")
        println("Warmup speedup       : ", round(speedup; digits = 3), "x")
      end
    end
  end

  summary_line = _compact_run_summary(status_ref[])
  cfg.output.console_summary && println(summary_line)
  open(logfile, "a") do io
    println(io, summary_line)
  end

  write_perf_log = cfg.performance.write_to_logfile && cfg.output.logfile_performance !== :off
  try
    emit_performance_summary(profile; logfile = logfile, print_to_console = cfg.performance.print_to_console, write_to_logfile = write_perf_log, max_rows = cfg.performance.max_diagnostic_rows, console_level = cfg.performance.compact_logging ? :compact : cfg.performance.level, logfile_level = cfg.output.logfile_performance)
    println("Wrote log file: ", logfile)
  catch err
    if err isa InterruptException
      println(stderr, "Interrupted by user while writing console output.")
      println(stderr, "Log file already written: ", logfile)
      rethrow()
    end
    rethrow()
  end
  return (status_ref = status_ref[], logfile = logfile)
end

function run_synthetic_tiled_grid_pf_perf(; config_file::AbstractString = "", args = String[], outdir::AbstractString = "")
  !isempty(strip(config_file)) && load_sparlectra_config!(String(config_file); reload = true)
  limits = isempty(args) ? [20] : Int[parse(Int, a) for a in args if !startswith(a, "--")]
  isempty(limits) && (limits = [20])
  dir = isempty(strip(outdir)) ? joinpath(pkgdir(Sparlectra), "examples", "_out") : String(outdir)
  mkpath(dir)
  logfile = joinpath(dir, "synthetic_tiled_grid_pf_perf_$(Dates.format(now(), "yyyymmdd_HHMMSS")).log")
  rows = NamedTuple[]
  for limit in limits
    net, meta = build_synthetic_tiled_grid_net(limit)
    iterations, erg, et = run_net_acpflow(net = net, show_results = false)
    push!(rows, (limit = limit, nbus = meta.actual_buses, converged = erg == 0, iterations = iterations, solve_ms = et * 1000.0))
  end
  open(logfile, "w") do io
    println(io, "limit nbus converged iterations solve_ms")
    for r in rows
      println(io, "$(r.limit) $(r.nbus) $(r.converged) $(r.iterations) $(round(r.solve_ms;digits=3))")
    end
  end
  println("Synthetic tiled-grid benchmark")
  for r in rows
    println("$(r.limit) $(r.nbus) $(r.converged) $(r.iterations) $(round(r.solve_ms;digits=3))")
  end
  println("Wrote log file: ", logfile)
  return (rows = rows, logfile = logfile)
end

function run_voltage_dependent_control_demo(; config_file::AbstractString = "", plot_curve::Union{Bool,Nothing} = nothing)
  !isempty(strip(config_file)) && load_sparlectra_config!(String(config_file); reload = true)
  qu_points = [(103.4, 35.0), (110.0, 0.0), (116.6, -25.0)]
  pu_points = [(103.4, 25.0), (110.0, 10.0), (116.6, 0.0)]
  net = Net(name = "voltage_dependent_control", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "Prosumer", vn_kV = 110.0)
  addACLine!(net = net, fromBus = "Slack", toBus = "Prosumer", length = 30.0, r = 0.05, x = 0.45)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
  qu_curve = make_characteristic(qu_points; voltage_unit = :kV, value_unit = :MVAr, vn_kV = 110.0, sbase_MVA = net.baseMVA)
  pu_curve = make_characteristic(pu_points; voltage_unit = :kV, value_unit = :MW, vn_kV = 110.0, sbase_MVA = net.baseMVA)
  addProsumer!(net = net, busName = "Prosumer", type = "SYNCHRONOUSMACHINE", p = 10.0, q = 0.0, qu_controller = QUController(qu_curve; qmin_MVAr = -50.0, qmax_MVAr = 50.0, sbase_MVA = net.baseMVA), pu_controller = PUController(pu_curve; pmin_MW = 0.0, pmax_MW = 50.0, sbase_MVA = net.baseMVA))
  addProsumer!(net = net, busName = "Prosumer", type = "ENERGYCONSUMER", p = 45.0, q = 18.0)
  ite, erg, _ = run_net_acpflow(net = net, max_ite = 40, tol = 1e-9, verbose = 1, opt_sparse = true, method = :rectangular, show_results = false)
  println("Voltage-dependent control demo converged=", erg == 0, " iterations=", ite, " plot_curve=", isnothing(plot_curve) ? "default" : string(plot_curve))
  return net
end

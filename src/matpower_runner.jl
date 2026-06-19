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

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 20.5.2026
# file: src/matpower_runner.jl

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

"""
    _record_perf!(profile, label, seconds)

Append a benchmark/performance event to `profile[:events]` when a mutable
profile dictionary is available.
"""
function _record_perf!(profile, label::Symbol, seconds::Real)
  profile isa AbstractDict || return profile
  events = get!(profile, :events, NamedTuple[])
  push!(events, (label = label, seconds = Float64(seconds)))
  return profile
end

_perf_timing_seconds(row)::Float64 = row isa NamedTuple && hasproperty(row, :elapsed_s) ? Float64(getproperty(row, :elapsed_s)) : 0.0
_perf_timing_bytes(row)::Int = row isa NamedTuple && hasproperty(row, :bytes) ? Int(getproperty(row, :bytes)) : 0

const _SESSION_REPRESENTATIVE_TIMING_SEEN = Set{Tuple{Symbol,Symbol}}()

function _auto_timing_mode(method::Symbol, warmup_runs::Int)::Symbol
  warmup_runs > 0 && return :warm_steady_state

  active_project = Base.active_project()
  project_key = isnothing(active_project) ? "__no_project__" : active_project
  key = (method, project_key)

  if key in _SESSION_REPRESENTATIVE_TIMING_SEEN
    return :auto_warm_session_reuse
  end

  push!(_SESSION_REPRESENTATIVE_TIMING_SEEN, key)
  return :auto_cold_session_first
end

function _timing_mode_label(mode::Symbol)::String
  mode === :warm_steady_state && return "warm / steady-state"
  mode === :auto_warm_session_reuse && return "auto: warm session-reuse (already executed in this Julia session)"
  mode === :auto_cold_session_first && return "auto: cold first session-run (first execution in this Julia session)"
  return "cold / representative"
end

function _is_top_level_perf_phase(phase::Symbol)::Bool
  return phase in (:network_construction, :start_projection, :solver_total, :result_output)
end

"""
    _sum_perf_timings(profile; top_level_only=false)

Return accumulated timing seconds from `profile[:timings]`.

When `top_level_only=true`, only major runner phases are included so coverage
can be compared against representative wall time.
"""
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

"""
    _sum_perf_events(profile)

Return total seconds across benchmark events recorded in `profile[:events]`.
"""
function _sum_perf_events(profile)::Float64
  profile isa AbstractDict || return 0.0
  events = get(profile, :events, NamedTuple[])
  total = 0.0
  for row in events
    total += row isa NamedTuple && hasproperty(row, :seconds) ? Float64(getproperty(row, :seconds)) : 0.0
  end
  return total
end

"""
    _print_timing_coverage(io, profile; level=:compact)

Print coverage diagnostics that compare representative wall time against
recorded phase timings and benchmark-event totals.
"""
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
  warmup_runs = Int(get(profile, :representative_warmup_runs, 0))
  warmup_override_reason = haskey(profile, :warmup_override_reason) ? String(get(profile, :warmup_override_reason, "")) : ""
  println(io, "Timing mode                    : ", _timing_mode_label(timing_mode))
  println(io, "Representative warmup runs     : ", warmup_runs)
  timing_mode in (:cold_representative, :auto_cold_session_first) && println(io, "Timing note                    : first representative run may include JIT, sparse/BLAS init, and cache construction")
  if warmup_runs == 0
    if !isempty(warmup_override_reason)
      println(io, "Warmup override reason         : ", warmup_override_reason)
    else
      println(io, "Warmup guidance                : This timing is a cold first-run measurement. To measure steady-state runtime, set performance.representative_warmup_runs >= 1.")
    end
  end
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
  if haskey(profile, :last_warmup_elapsed_s)
    println(io)
    println(io, "Warmup")
    println(io, "------")
    println(io, "runs              : ", warmup_runs)
    println(io, "last_warmup_time  : ", round(Float64(get(profile, :last_warmup_elapsed_s, 0.0)); digits = 6), " s")
  end
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

"""
    print_performance_profile(io, profile; ...)

Render a human-readable performance report from the MATPOWER runner profile.

Supports compact and full output levels, optional allocation columns, and
optional iteration diagnostics.
"""
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
    for i = 1:rows
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
    for i = 1:rows
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

"""
    emit_performance_summary(profile; ...)

Emit performance output to console and/or logfile according to runner output
configuration.
"""
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

"""
    parse_runtime_threads_request(raw)

Parse runtime-thread configuration tokens (`"keep"`, `"auto"`, integer-like
strings) into `Int` thread counts or `nothing`.
"""
function parse_runtime_threads_request(raw)::Union{Int,Nothing}
  sval = lowercase(strip(string(raw)))
  sval in ("", "keep", "default", "false", "off") && return nothing
  sval == "auto" && return Sys.CPU_THREADS
  iv = tryparse(Int, sval)
  isnothing(iv) && return nothing
  return iv > 0 ? iv : nothing
end

"""
    runtime_thread_status(cfg)

Resolve requested Julia/BLAS thread settings, apply BLAS runtime changes, and
return a status tuple for reporting.
"""
function runtime_thread_status(cfg::RuntimeConfig)
  requested_julia = parse_runtime_threads_request(cfg.julia_threads)
  active_julia = Threads.nthreads()
  julia_applied = isnothing(requested_julia) || requested_julia == active_julia
  julia_startup_required = !isnothing(requested_julia) && requested_julia != active_julia

  req_blas = lowercase(strip(cfg.blas_threads))
  requested_blas = parse_runtime_threads_request(cfg.blas_threads)
  before = BLAS.get_num_threads()
  # BLAS threads are runtime-mutable; Julia threads are startup-only.
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

"""
    print_runtime_thread_config(io, status)

Print runtime thread request/active/apply status, including startup guidance
when Julia thread requests cannot be applied in-process.
"""
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

"""
    matpower_run_logfile_path(casefile, cfg)

Build the per-run MATPOWER logfile path under `examples/_out/`.
"""
function matpower_run_logfile_path(casefile::AbstractString, cfg::OutputConfig)
  log_dir = joinpath(pkgdir(Sparlectra), "examples", "_out")
  mkpath(log_dir)
  stem = splitext(basename(casefile))[1]
  return joinpath(log_dir, "run_$(stem)_$(Dates.format(now(), "yyyymmdd_HHMMSS")).log")
end

"""
    print_matpower_runner_header(io; ...)

Print the standard MATPOWER runner banner with resolved config, case, method,
and logging/performance settings.
"""
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

"""
    run_with_output_capture(f; capture_stdout=true, capture_stderr=true)

Run `f()` with redirected stdio and return `(result, stdout, stderr)`.

Used by the runner to isolate diagnostic output for logfile routing.
"""
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

"""
    run_silent_for_benchmark(f)

Execute benchmark payload `f` via `invokelatest` while suppressing captured
stdio output.
"""
run_silent_for_benchmark(f::Function) = run_with_output_capture(() -> Base.invokelatest(f); capture_stdout = true, capture_stderr = true).result

"""
    append_captured_output_to_logfile(logfile, captured; section_title="")

Append captured stdout/stderr chunks to the MATPOWER logfile.
"""
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

"""
    _compact_run_summary(status)

Build a single-line compact summary string from a MATPOWER run status object.
"""
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
  qlimit_active_set_changes = hasproperty(status, :qlimit_active_set_changes) ? getproperty(status, :qlimit_active_set_changes) : 0
  qlimit_reenable_events = hasproperty(status, :qlimit_reenable_events) ? getproperty(status, :qlimit_reenable_events) : 0
  benchmark_fragment = isnothing(benchmark_median_s) ? "" : " benchmark_median=$(round(benchmark_median_s * 1000.0; digits = 6)) ms"
  solver_fragment = isnothing(solver_elapsed_s) ? " solver_time=unavailable" : " solver_time=$(round(Float64(solver_elapsed_s); digits = 6)) s"
  result_output_s = hasproperty(status, :result_output_s) ? getproperty(status, :result_output_s) : nothing
  result_output_fragment = isnothing(result_output_s) ? "" : " result_output_time=$(round(Float64(result_output_s); digits = 6)) s"
  return "summary method=$(method) outcome=$(outcome) numerical_solution=$(numerical_solution) solution_available=$(solution_available) limit_validation=$(limit_validation) final_converged=$(final_converged) final_mismatch=$(round(final_mismatch; digits = 9)) iterations=$(iterations) representative_time=$(round(representative_elapsed_s; digits = 6)) s$(solver_fragment)$(result_output_fragment)$(benchmark_fragment) pv2pq_events=$(pv2pq_events) pv2pq_buses=$(pv2pq_buses) qlimit_active_set_changes=$(qlimit_active_set_changes) qlimit_reenable_events=$(qlimit_reenable_events) reason=\"$(reason_text)\""
end

"""
    _run_matpower_single(local_case, cfg, profile, status_ref)

Run one MATPOWER power-flow solve and update `status_ref` with elapsed-time and
result-output timing metadata.
"""
function _run_matpower_single(local_case::AbstractString, cfg::SparlectraConfig, profile, status_ref; show_results::Bool = true)
  result = Ref{SparlectraRunResult}()
  t = @elapsed result[] = _run_sparlectra(; casefile = String(local_case), config = cfg, performance_profile = profile, emit_output = show_results)
  status_ref[] = _sparlectra_status_row(result[])
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

"""
    _run_matpower_single_routed(local_case, cfg, profile, status_ref; ...)

Run one solve either directly to console or through captured stdio for later
logfile emission.
"""
function _run_matpower_single_routed(local_case::AbstractString, cfg::SparlectraConfig, profile, status_ref; to_console::Bool = false, capture_for_log::Bool = false, show_results::Bool = true)
  if to_console && !capture_for_log
    return _run_matpower_single(local_case, cfg, profile, status_ref; show_results = show_results)
  end
  captured = run_with_output_capture(capture_stdout = true, capture_stderr = true) do
    _run_matpower_single(local_case, cfg, profile, status_ref; show_results = show_results)
  end
  return (elapsed = captured.result, captured = captured)
end

"""
    with_powerflow_method(cfg, method)

Create a configuration copy with `powerflow.method` replaced by `method`.
"""
function with_powerflow_method(cfg::SparlectraConfig, method::Symbol)::SparlectraConfig
  pf = cfg.powerflow
  pf_kwargs = NamedTuple{fieldnames(PowerFlowConfig)}(getfield.(Ref(pf), fieldnames(PowerFlowConfig)))
  pf2 = PowerFlowConfig(; pf_kwargs..., method = method)
  return SparlectraConfig(; powerflow = pf2, state_estimation = cfg.state_estimation, matpower = cfg.matpower, performance = cfg.performance, benchmark = cfg.benchmark, runtime = cfg.runtime, diagnostics = cfg.diagnostics, output = cfg.output, control = cfg.control)
end

_auto_profile_value_string(x)::String = x isa AbstractFloat ? string(round(x; digits = 8)) : string(x)

function _matpower_import_auto_profile_residual_score(stats)::Float64
  stats isa NamedTuple || return Inf
  get(stats, :ok, false) || return Inf
  return max(Float64(get(stats, :max_p_mis_pu, Inf)), Float64(get(stats, :max_q_mis_pu, Inf)))
end

function _matpower_import_auto_profile_convention_scan(mpc)
  rows = NamedTuple[]
  for shift_unit in (:deg, :rad), shift_sign in (1.0, -1.0), ratio in (:normal, :reciprocal)
    stats = MatpowerIO.vmva_power_mismatch_stats(mpc; matpower_shift_unit = shift_unit, matpower_shift_sign = shift_sign, matpower_ratio = ratio)
    push!(rows, (shift_unit = shift_unit, shift_sign = shift_sign, ratio = ratio, stats = stats, score = _matpower_import_auto_profile_residual_score(stats)))
  end
  sort!(rows; by = row -> row.score)
  return rows
end

function _matpower_import_auto_profile_shunt_score(mpc, include_shunts::Bool; shift_unit, shift_sign, ratio)::Float64
  bus = copy(mpc.bus)
  if !include_shunts && size(bus, 2) >= 6
    bus[:, 5] .= 0.0
    bus[:, 6] .= 0.0
  end
  mpc2 = MatpowerIO.MatpowerCase(mpc.name, mpc.baseMVA, bus, mpc.gen, mpc.branch, mpc.gencost, mpc.bus_name)
  return _matpower_import_auto_profile_residual_score(MatpowerIO.vmva_power_mismatch_stats(mpc2; matpower_shift_unit = shift_unit, matpower_shift_sign = shift_sign, matpower_ratio = ratio))
end

function _matpower_import_auto_profile_action(mode::Symbol, current, recommended; safe_to_apply::Bool)::Symbol
  current == recommended && return :keep
  mode === :apply && safe_to_apply && return :applied
  return :recommend
end

function _push_auto_profile_row!(rows, mode::Symbol, option::AbstractString, current, recommended; safe_to_apply::Bool, reason::AbstractString, evidence::AbstractString = "")
  push!(
    rows,
    (
      option = String(option),
      current = _auto_profile_value_string(current),
      recommended = _auto_profile_value_string(recommended),
      action = _matpower_import_auto_profile_action(mode, current, recommended; safe_to_apply = safe_to_apply),
      reason = String(reason),
      evidence = String(evidence),
      raw_current = current,
      raw_recommended = recommended,
    ),
  )
  return rows
end

function _matpower_import_auto_profile_has_vg_mismatch(mpc; tol::Float64)
  rows = MatpowerIO.pv_voltage_reference_rows(mpc; matpower_pv_voltage_source = :gen_vg, tol = tol, warn = false)
  mismatches = [row for row in rows if isfinite(row.dvg_bus) && abs(row.dvg_bus) > tol]
  return (rows = rows, mismatches = mismatches)
end

function _copy_matpower_with(mat::MatpowerImportConfig; kwargs...)
  fields = NamedTuple{fieldnames(MatpowerImportConfig)}(getfield.(Ref(mat), fieldnames(MatpowerImportConfig)))
  return MatpowerImportConfig(; fields..., kwargs...)
end

function _copy_config_with(cfg::SparlectraConfig; powerflow = cfg.powerflow, matpower = cfg.matpower)
  return SparlectraConfig(; powerflow = powerflow, state_estimation = cfg.state_estimation, matpower = matpower, performance = cfg.performance, benchmark = cfg.benchmark, runtime = cfg.runtime, diagnostics = cfg.diagnostics, output = cfg.output, control = cfg.control)
end

function matpower_import_auto_profile(mpc, cfg::SparlectraConfig; mode::Symbol = cfg.matpower.auto_profile)
  mode in (:off, :recommend, :apply) || throw(ArgumentError("matpower auto-profile mode must be off, recommend, or apply."))
  mode === :off && return (config = cfg, rows = NamedTuple[], applied = NamedTuple[])

  mat = cfg.matpower
  pf = cfg.powerflow
  rows = NamedTuple[]
  applied_pairs = Pair{Symbol,Any}[]

  convention_rows = _matpower_import_auto_profile_convention_scan(mpc)
  current_row = only([row for row in convention_rows if row.shift_unit === mat.shift_unit && row.shift_sign == mat.shift_sign && row.ratio === mat.ratio])
  best = first(convention_rows)
  current_score = current_row.score
  best_score = best.score
  strong_convention = isfinite(current_score) && isfinite(best_score) && best_score + 1.0e-8 < min(current_score * 0.25, current_score - 1.0e-6)
  convention_reason = strong_convention ? "residual scan found a clearly better MATPOWER branch convention" : "residual scan is ambiguous or current convention is best"
  convention_evidence = "current=$(round(current_score; digits = 8)) best=$(round(best_score; digits = 8))"
  rec_shift_unit = strong_convention ? best.shift_unit : mat.shift_unit
  rec_shift_sign = strong_convention ? best.shift_sign : mat.shift_sign
  rec_ratio = strong_convention ? best.ratio : mat.ratio
  _push_auto_profile_row!(rows, mode, "matpower_import.shift_unit", mat.shift_unit, rec_shift_unit; safe_to_apply = strong_convention, reason = convention_reason, evidence = convention_evidence)
  _push_auto_profile_row!(rows, mode, "matpower_import.shift_sign", mat.shift_sign, rec_shift_sign; safe_to_apply = strong_convention, reason = convention_reason, evidence = convention_evidence)
  _push_auto_profile_row!(rows, mode, "matpower_import.ratio", mat.ratio, rec_ratio; safe_to_apply = strong_convention, reason = convention_reason, evidence = convention_evidence)
  if mode === :apply && strong_convention
    push!(applied_pairs, :shift_unit => rec_shift_unit, :shift_sign => rec_shift_sign, :ratio => rec_ratio)
  end

  shunt_adm_score = _matpower_import_auto_profile_shunt_score(mpc, true; shift_unit = rec_shift_unit, shift_sign = rec_shift_sign, ratio = rec_ratio)
  shunt_vdi_score = _matpower_import_auto_profile_shunt_score(mpc, false; shift_unit = rec_shift_unit, shift_sign = rec_shift_sign, ratio = rec_ratio)
  current_shunt_score = mat.bus_shunt_model === :admittance ? shunt_adm_score : shunt_vdi_score
  best_shunt_model = shunt_adm_score <= shunt_vdi_score ? :admittance : :voltage_dependent_injection
  best_shunt_score = min(shunt_adm_score, shunt_vdi_score)
  strong_shunt = isfinite(current_shunt_score) && isfinite(best_shunt_score) && best_shunt_score + 1.0e-8 < min(current_shunt_score * 0.25, current_shunt_score - 1.0e-6)
  rec_shunt = strong_shunt ? best_shunt_model : mat.bus_shunt_model
  _push_auto_profile_row!(
    rows,
    mode,
    "matpower_import.bus_shunt_model",
    mat.bus_shunt_model,
    rec_shunt;
    safe_to_apply = strong_shunt,
    reason = strong_shunt ? "bus-shunt residual scan found a clearly better model" : "bus-shunt residual scan is ambiguous or current model is best",
    evidence = "admittance=$(round(shunt_adm_score; digits = 8)) voltage_dependent_injection=$(round(shunt_vdi_score; digits = 8))",
  )
  mode === :apply && strong_shunt && push!(applied_pairs, :bus_shunt_model => rec_shunt)

  vg = _matpower_import_auto_profile_has_vg_mismatch(mpc; tol = mat.pv_voltage_mismatch_tol_pu)
  has_online_vg = any(row -> !isempty(row.gen_vgs), vg.rows)
  rec_pv_source = has_online_vg ? :gen_vg : mat.pv_voltage_source
  pv_reason = has_online_vg ? "online GEN.VG values are available for PV/REF buses" : "no online GEN.VG values available for PV/REF buses"
  _push_auto_profile_row!(rows, mode, "matpower_import.pv_voltage_source", mat.pv_voltage_source, rec_pv_source; safe_to_apply = false, reason = pv_reason, evidence = "pv_ref_buses=$(length(vg.rows)) vg_bus_mismatches=$(length(vg.mismatches))")
  rec_compare = isempty(vg.mismatches) ? mat.compare_voltage_reference : :hybrid
  _push_auto_profile_row!(
    rows,
    mode,
    "matpower_import.compare_voltage_reference",
    mat.compare_voltage_reference,
    rec_compare;
    safe_to_apply = false,
    reason = isempty(vg.mismatches) ? "BUS.VM and GEN.VG are consistent within tolerance" : "BUS.VM / GEN.VG mismatch detected; hybrid comparison keeps PQ buses on BUS.VM",
    evidence = "vg_bus_mismatches=$(length(vg.mismatches)) tol=$(mat.pv_voltage_mismatch_tol_pu)",
  )

  nbus = size(mpc.bus, 1)
  npv = count(r -> Int(mpc.bus[r, 2]) == 2, axes(mpc.bus, 1))
  fragile_start = nbus >= 1000 || npv >= 100 || (isfinite(best_score) && best_score > 1.0e-3)
  if fragile_start
    robust_start = (angle_mode = :dc, voltage_mode = :bus_vm_va_blend, start_projection = true, try_dc_start = true, try_blend_scan = true, branch_guard = true, measure_candidates = true, reuse_import_data = true)
    for (field, rec) in pairs(robust_start)
      cur = getfield(pf.start_mode, field)
      _push_auto_profile_row!(rows, mode, "power_flow.start_mode.$(field)", cur, rec; safe_to_apply = false, reason = "large or fragile MATPOWER case; robust imported/DC start is safer", evidence = "nbus=$(nbus) npv=$(npv) best_residual=$(round(best_score; digits = 8))")
    end
  end

  q_ranges = Float64[]
  if size(mpc.gen, 2) >= 5
    for g in axes(mpc.gen, 1)
      size(mpc.gen, 2) >= 8 && mpc.gen[g, 8] <= 0.0 && continue
      push!(q_ranges, abs(Float64(mpc.gen[g, 4] - mpc.gen[g, 5])) / mpc.baseMVA)
    end
  end
  narrow_q = count(x -> x <= 0.02, q_ranges)
  qlimit_risk = npv >= 100 || narrow_q >= 1
  if qlimit_risk
    q_recs = (
      start_iter = 3,
      start_mode = :iteration_or_auto,
      hysteresis_pu = 0.01,
      cooldown_iters = 1,
      guard = true,
      guard_min_q_range_pu = 0.02,
      guard_narrow_range_mode = :lock_pq,
      guard_zero_range_mode = :lock_pq,
      guard_violation_mode = :lock_pq,
      guard_max_switches = 10,
      guard_max_remaining_violations = 0,
    )
    for (field, rec) in pairs(q_recs)
      cur = getfield(pf.qlimits, field)
      option = field in (:guard_min_q_range_pu, :guard_narrow_range_mode, :guard_zero_range_mode, :guard_violation_mode, :guard_max_switches, :guard_max_remaining_violations) ? "power_flow.qlimits.guard.$(replace(String(field), "guard_" => ""))" : "power_flow.qlimits.$(field)"
      _push_auto_profile_row!(rows, mode, option, cur, rec; safe_to_apply = false, reason = "many PV buses or narrow Q ranges; conservative Q-limit guard recommended", evidence = "npv=$(npv) narrow_q_ranges=$(narrow_q)")
    end
  end

  mat2 = isempty(applied_pairs) ? mat : _copy_matpower_with(mat; applied_pairs...)
  cfg2 = mat2 === mat ? cfg : _copy_config_with(cfg; matpower = mat2)
  return (config = cfg2, rows = rows, applied = Tuple(applied_pairs))
end

function print_matpower_import_auto_profile(io::IO, rows)
  isempty(rows) && return nothing
  println(io, "MATPOWER auto-profile recommendations")
  println(io, "-------------------------------------")
  @printf(io, "%-46s %-14s %-14s %-10s %s\n", "option", "current", "recommended", "action", "reason")
  for row in rows
    reason = isempty(row.evidence) ? row.reason : string(row.reason, " (", row.evidence, ")")
    @printf(io, "%-46s %-14s %-14s %-10s %s\n", row.option, row.current, row.recommended, string(row.action), reason)
  end
  println(io)
  return nothing
end

function _auto_profile_row_value(rows, option::AbstractString, field::Symbol)
  for row in rows
    row.option == option && return getfield(row, field)
  end
  return ""
end

function print_matpower_import_auto_profile_summary(io::IO, auto_profile_result, cfg::SparlectraConfig)
  rows = auto_profile_result.rows
  isempty(rows) && return nothing
  println(io, "Original MATPOWER import options:")
  println(io, "  auto_profile = ", cfg.matpower.auto_profile)
  println(io, "  ratio        = ", _auto_profile_row_value(rows, "matpower_import.ratio", :current))
  println(io, "  shift_unit   = ", _auto_profile_row_value(rows, "matpower_import.shift_unit", :current))
  println(io, "  shift_sign   = ", _auto_profile_row_value(rows, "matpower_import.shift_sign", :current))
  println(io)
  println(io, "Auto-profile recommendation:")
  println(io, "  ratio        = ", _auto_profile_row_value(rows, "matpower_import.ratio", :recommended))
  println(io, "  shift_unit   = ", _auto_profile_row_value(rows, "matpower_import.shift_unit", :recommended))
  println(io, "  shift_sign   = ", _auto_profile_row_value(rows, "matpower_import.shift_sign", :recommended))
  println(io)
  println(io, "Final effective MATPOWER import options:")
  println(io, "  auto_profile = ", cfg.matpower.auto_profile)
  println(io, "  ratio        = ", cfg.matpower.ratio)
  println(io, "  shift_unit   = ", cfg.matpower.shift_unit)
  println(io, "  shift_sign   = ", cfg.matpower.shift_sign)
  if cfg.matpower.auto_profile === :apply && isempty(auto_profile_result.applied)
    println(io, "  applied      = none (all safe recommendations already matched the runtime options or the scan was ambiguous)")
  elseif cfg.matpower.auto_profile === :apply
    println(io, "  applied      = ", join((string(first(pair), "=", last(pair)) for pair in auto_profile_result.applied), ", "))
  end
  println(io)
  return nothing
end

function print_matpower_import_auto_profile_effective_options(io::IO, cfg::SparlectraConfig)
  println(io, "Final effective MATPOWER auto-profile options")
  println(io, "---------------------------------------------")
  println(io, "matpower_import.auto_profile: ", cfg.matpower.auto_profile)
  println(io, "matpower_import.shift_unit: ", cfg.matpower.shift_unit)
  println(io, "matpower_import.shift_sign: ", cfg.matpower.shift_sign)
  println(io, "matpower_import.ratio: ", cfg.matpower.ratio)
  println(io, "matpower_import.bus_shunt_model: ", cfg.matpower.bus_shunt_model)
  println(io, "matpower_import.pv_voltage_source: ", cfg.matpower.pv_voltage_source)
  println(io, "matpower_import.compare_voltage_reference: ", cfg.matpower.compare_voltage_reference)
  println(io, "power_flow.start_mode.angle_mode: ", cfg.powerflow.start_mode.angle_mode)
  println(io, "power_flow.start_mode.voltage_mode: ", cfg.powerflow.start_mode.voltage_mode)
  println(io, "power_flow.start_mode.start_projection: ", cfg.powerflow.start_mode.start_projection)
  println(io, "power_flow.start_mode.try_dc_start: ", cfg.powerflow.start_mode.try_dc_start)
  println(io, "power_flow.start_mode.try_blend_scan: ", cfg.powerflow.start_mode.try_blend_scan)
  println(io, "power_flow.start_mode.branch_guard: ", cfg.powerflow.start_mode.branch_guard)
  println(io, "power_flow.start_mode.measure_candidates: ", cfg.powerflow.start_mode.measure_candidates)
  println(io, "power_flow.start_mode.reuse_import_data: ", cfg.powerflow.start_mode.reuse_import_data)
  println(io, "power_flow.qlimits.start_iter: ", cfg.powerflow.qlimits.start_iter)
  println(io, "power_flow.qlimits.start_mode: ", cfg.powerflow.qlimits.start_mode)
  println(io, "power_flow.qlimits.hysteresis_pu: ", cfg.powerflow.qlimits.hysteresis_pu)
  println(io, "power_flow.qlimits.cooldown_iters: ", cfg.powerflow.qlimits.cooldown_iters)
  println(io, "power_flow.qlimits.guard.enabled: ", cfg.powerflow.qlimits.guard)
  println(io)
  return nothing
end

function run_matpower_import_auto_profile(mpc, cfg::SparlectraConfig)
  cfg.matpower.auto_profile === :off && return (config = cfg, rows = NamedTuple[], applied = NamedTuple[])
  return matpower_import_auto_profile(mpc, cfg; mode = cfg.matpower.auto_profile)
end

function write_matpower_import_auto_profile(io::IO, auto_profile_result, cfg::SparlectraConfig; casefile::AbstractString = "")
  !isempty(casefile) && println(io, "Runtime casefile: ", casefile, "\n")
  print_matpower_import_auto_profile_summary(io, auto_profile_result, cfg)
  print_matpower_import_auto_profile(io, auto_profile_result.rows)
  print_matpower_import_auto_profile_effective_options(io, cfg)
  return nothing
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
  auto_profile_result = nothing
  if cfg.matpower.auto_profile !== :off
    mpc_for_profile = MatpowerIO.read_case(local_case; legacy_compat = true)
    auto_profile_result = run_matpower_import_auto_profile(mpc_for_profile, cfg)
    cfg = auto_profile_result.config
  end
  print_matpower_runner_header(
    stdout;
    default_config = DEFAULT_SPARLECTRA_CONFIG_PATH,
    user_config = isfile(USER_SPARLECTRA_CONFIG_PATH) ? USER_SPARLECTRA_CONFIG_PATH : "",
    casefile = local_case,
    logfile = logfile,
    methods = cfg.benchmark.methods,
    benchmark = cfg.benchmark,
    performance = cfg.performance,
  )
  open(logfile, "w") do io
    print_matpower_runner_header(
      io;
      default_config = DEFAULT_SPARLECTRA_CONFIG_PATH,
      user_config = isfile(USER_SPARLECTRA_CONFIG_PATH) ? USER_SPARLECTRA_CONFIG_PATH : "",
      casefile = local_case,
      logfile = logfile,
      methods = cfg.benchmark.methods,
      benchmark = cfg.benchmark,
      performance = cfg.performance,
    )
  end
  if auto_profile_result !== nothing && cfg.matpower.auto_profile_log
    write_matpower_import_auto_profile(stdout, auto_profile_result, cfg; casefile = local_case)
    open(logfile, "a") do io
      write_matpower_import_auto_profile(io, auto_profile_result, cfg; casefile = local_case)
      print_effective_config(io, cfg)
    end
  end
  if cfg.runtime.print_thread_config
    status = runtime_thread_status(cfg.runtime)
    print_runtime_thread_config(stdout, status)
    open(logfile, "a") do io
      print_runtime_thread_config(io, status)
    end
  end
  cfg.diagnostics.log_effective_config && open(logfile, "a") do io
    print_effective_config(io, cfg)
  end

  profile = new_performance_profile(cfg.performance)
  status_ref = Ref{Any}(nothing)

  # Benchmark mode runs one representative solve per method and then measures
  # steady repeated solves separately so benchmark statistics are not polluted by
  # runner/logging output.
  if cfg.benchmark.enabled
    profile[:representative_warmup_runs] = max(0, cfg.performance.representative_warmup_runs)
    if cfg.performance.representative_warmup_runs > 0
      profile[:warmup_override_reason] = "Benchmark mode uses BenchmarkTools repeated samples; dedicated representative warmup loops are not run."
    end
    open(logfile, "a") do io
      println(io, "Benchmark results")
    end
    println("Benchmark results")
    for method in cfg.benchmark.methods
      method_cfg = with_powerflow_method(cfg, method)
      warm_run = _run_matpower_single_routed(local_case, method_cfg, profile, status_ref; to_console = cfg.output.console_diagnostics === :full, capture_for_log = cfg.output.logfile_diagnostics !== :off || cfg.output.logfile_warnings !== :off)
      warm = warm_run isa NamedTuple ? warm_run.elapsed : warm_run
      profile[:representative_elapsed_s] = warm
      if warm_run isa NamedTuple && cfg.output.logfile_diagnostics === :full
        append_captured_output_to_logfile(logfile, warm_run.captured; section_title = "Representative solve diagnostics (method=$(method))")
      end
      bench_status_ref = Ref{Any}(nothing)
      bench = BenchmarkTools.@benchmarkable run_silent_for_benchmark() do
        Base.invokelatest(_run_sparlectra; casefile = String($local_case), config = $method_cfg, performance_profile = nothing, emit_output = false)
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
    # Non-benchmark mode keeps one representative solve and optional warmup
    # preparation for cold-vs-warm comparison diagnostics.
    cold_solver_time = nothing
    warm_solver_time = nothing
    warmup_runs = max(0, cfg.performance.representative_warmup_runs)
    cold_profile = nothing
    profile[:representative_warmup_runs] = warmup_runs
    if warmup_runs > 0
      last_warmup_elapsed_s = 0.0
      for _ = 1:warmup_runs
        warm_status_ref = Ref{Any}(nothing)
        cold_profile = Dict{Symbol,Any}()
        warmup_elapsed_s = _run_matpower_single(local_case, cfg, cold_profile, warm_status_ref; show_results = false)
        last_warmup_elapsed_s = warmup_elapsed_s
        if cfg.performance.compare_cold_warm
          st = warm_status_ref[]
          if !isnothing(st) && hasproperty(st, :solver_elapsed_s) && isnothing(cold_solver_time)
            cold_solver_time = getproperty(st, :solver_elapsed_s)
          end
        end
      end
      profile[:last_warmup_elapsed_s] = last_warmup_elapsed_s
    end
    profile[:timing_mode] = _auto_timing_mode(cfg.powerflow.method, warmup_runs)
    run_single = _run_matpower_single_routed(local_case, cfg, profile, status_ref; to_console = cfg.output.console_diagnostics === :full, capture_for_log = cfg.output.logfile_diagnostics !== :off || cfg.output.logfile_warnings !== :off)
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
      _run_matpower_single(local_case, cfg, cold_profile, cold_status_ref; show_results = false)
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
    emit_performance_summary(
      profile;
      logfile = logfile,
      print_to_console = cfg.performance.print_to_console,
      write_to_logfile = write_perf_log,
      max_rows = cfg.performance.max_diagnostic_rows,
      console_level = cfg.performance.compact_logging ? :compact : cfg.performance.level,
      logfile_level = cfg.output.logfile_performance,
    )
    println("Wrote log file: ", logfile)
  catch err
    if err isa InterruptException
      println(stderr, "Interrupted by user while writing console output.")
      println(stderr, "Log file already written: ", logfile)
      rethrow()
    end
    rethrow()
  end
  return (status_ref = status_ref[], logfile = logfile, performance_profile = profile)
end

function _synthetic_pf_perf_row(limit, meta, result)
  return (
    limit = limit,
    nbus = meta.actual_buses,
    converged = result.final_converged,
    numerical_converged = result.numerical_converged,
    solution_available = result.solution_available,
    outcome = result.outcome,
    reason = result.reason,
    iterations = result.iterations,
    solve_ms = result.elapsed_s * 1000.0,
  )
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
    result = _run_sparlectra(net = net; emit_output = false)
    push!(rows, _synthetic_pf_perf_row(limit, meta, result))
  end
  open(logfile, "w") do io
    println(io, "limit nbus converged numerical_converged solution_available outcome reason iterations solve_ms")
    for r in rows
      println(io, "$(r.limit) $(r.nbus) $(r.converged) $(r.numerical_converged) $(r.solution_available) $(r.outcome) $(r.reason) $(r.iterations) $(round(r.solve_ms;digits=3))")
    end
  end
  println("Synthetic tiled-grid benchmark")
  println("limit nbus converged numerical_converged solution_available outcome reason iterations solve_ms")
  for r in rows
    println("$(r.limit) $(r.nbus) $(r.converged) $(r.numerical_converged) $(r.solution_available) $(r.outcome) $(r.reason) $(r.iterations) $(round(r.solve_ms;digits=3))")
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
  cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 40, tol = 1e-9), output = OutputConfig(logfile_results = :off))
  result = run_sparlectra(net = net, config = cfg)
  println("Voltage-dependent control demo converged=", result.numerical_converged, " iterations=", result.iterations, " plot_curve=", isnothing(plot_curve) ? "default" : string(plot_curve))
  return net
end

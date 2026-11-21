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
#
# file: benchmark/gpu/batch_study_gpu.jl
# purpose: GPU study, stage 2: throughput of the injection batch
#          (runInjectionBatch) on one imported network, import once per
#          case and not timed. Workloads: seeded random load scaling per bus
#          (band +-20 percent) with N = 1, 8, 64, 256, 1024, and the uniform
#          ramp of 96 steps; START = net (the net's start iterate) or base
#          (warm start from the solved base case). Variants: C1 batch
#          (default configuration) and C2 batch (power mode, power_mode_lu
#          auto, KLU loaded, workers warmed by a base solve), both over
#          Julia tasks with BLAS at 1 thread inside, at the task count that
#          was fastest in a scan at N = 64 (TASKS); GPU with the strategies
#          sequential (deterministic, the default), sequential_fast,
#          sequential_f32 (float32_mixed, deterministic), nonuniform
#          (deterministic), uniform (fast, ubatch_mask on) and
#          uniform_nomask; the two-phase variants g5, g5r and g6
#          (gpu_two_phase, switch at SWITCH). Per run the wall time of the whole call (setup
#          included, the GPU synchronized before and after); the first call
#          is reported on its own, the median is over REPEATS further calls.
#          A GPU variant that runs out of device memory at some N gets its
#          limit searched by bisection (single calls); every N above the
#          limit of the default is skipped for all variants.
#          Checks: the C1 batch equals single runpf! solves of its first
#          scenarios bit for bit; every run against the C1 batch (max abs
#          dVm, max abs dVa, convergence, iterations, Q-limit counts); the
#          repeats of a variant against each other (SHA-1 of vm and va).
#          Per-scenario residuals (Float64) of the last call of every row go
#          to OUT_SCEN, for the two-phase variants with the steps per phase,
#          the switch reason and the reduction rate rho per step.
#
# Run: julia -t <threads> --project=benchmark/gpu benchmark/gpu/batch_study_gpu.jl [case ...]
#   default: the nine cases of the study. Environment: NS (default
#   "1,8,64,256,1024"), REPEATS (default 5), SEED (default 1), START (net or
#   base, default net), STRATEGIES (default all six GPU variants), TASKS
#   (task counts of the CPU scan, default "4,8,16,24"), TASKS_FROM (a CSV
#   of an earlier run: take the fastest scan count per case and variant from
#   it and skip the scan), RAMP (default 1; 0
#   leaves the ramp out), SWITCH (phase-1 switch of the two-phase variants,
#   default 1e-5), EXECUTOR (the CPU executor of C1, C2 and the CPU phase
#   of g6, default threads, so the task count is a thread count), CPU_TIMED (default 1; 0 runs C1 once untimed as the
#   reference of the GPU rows and leaves C2 out, for a series of further GPU
#   variants), OUT and OUT_SCEN (CSV paths; default stdout only, no
#   scenario file). One CSV row per case, workload and variant on stdout,
#   behind a header with commit, Julia, threads, device, cuDSS version.

using Sparlectra, CUDA, CUDSS, KLU, LinearAlgebra, Printf, Dates, Statistics, SHA

const QUIET = OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)
const CASES = ["case14", "case118", "case1354pegase", "case2869pegase", "case6495rte", "case9241pegase", "case_ACTIVSg10k", "case_ACTIVSg25k", "case_ACTIVSg70k"]
const NS = parse.(Int, split(get(ENV, "NS", "1,8,64,256,1024"), ","))
const REPEATS = parse(Int, get(ENV, "REPEATS", "5"))
const SEED = parse(Int, get(ENV, "SEED", "1"))
const START = Symbol(get(ENV, "START", "net"))
const STRATEGIES = Symbol.(split(get(ENV, "STRATEGIES", "sequential,sequential_fast,sequential_f32,nonuniform,uniform,uniform_nomask"), ","))
const TASKS = parse.(Int, split(get(ENV, "TASKS", "4,8,16,24"), ","))
const SCAN_N = 64
const TASKS_FROM = get(ENV, "TASKS_FROM", "")
const SWITCH = parse(Float64, get(ENV, "SWITCH", "1e-5"))
const CPU_TIMED = get(ENV, "CPU_TIMED", "1") == "1"
const EXECUTOR = Symbol(get(ENV, "EXECUTOR", "threads"))
const TWO_PHASE = (:g5, :g5r, :g6)
# scenarios of the C1 batch checked against single runpf! solves
const SINGLE_CHECKS = 2

# a GPU run that the device cannot take: out of device memory, or a cuDSS
# status (allocation failure) reported as an error; the GPU variants run
# with gpu_fallback_to_cpu off, so these reach the script
gpu_refused(e) = e isa CUDA.OutOfGPUMemoryError || e isa CUDSS.CUDSSError
oom(e) = e isa CUDA.OutOfGPUMemoryError || (e isa CUDSS.CUDSSError && occursin("ALLOC", sprint(showerror, e)))

vhash(r) = bytes2hex(sha1(vcat(reinterpret(UInt8, vec(r.vm)), reinterpret(UInt8, vec(r.va_deg)))))[1:12]

function base_config(net)
  eff = Sparlectra._runpf_effective_config(net, nothing)
  return eff === nothing ? Sparlectra.powerflow_config() : (eff isa SparlectraConfig ? eff.powerflow : eff)
end

gpu_cb(; kw...) = Sparlectra.ComputeBackendConfig(; mode = :gpu, gpu_fallback_to_cpu = false, kw...)

function variant_config(pf, var::Symbol)
  var === :C2 && return Sparlectra._copy_powerflow_with(pf; power_mode = true, power_mode_lu = :auto)
  var === :C1 && return pf
  var === :sequential_fast && return Sparlectra._copy_powerflow_with(pf; compute_backend = gpu_cb(gpu_fast_mode = true))
  var === :sequential_f32 && return Sparlectra._copy_powerflow_with(pf; compute_backend = gpu_cb(gpu_precision = :float32_mixed))
  return Sparlectra._copy_powerflow_with(pf; compute_backend = gpu_cb())
end

gpu_strategy(var) = var in (:sequential, :sequential_fast, :sequential_f32) || var in TWO_PHASE ? :sequential : var === :uniform_nomask ? :uniform : var

# tasks: the CPU task count (C1, C2, and the CPU phase of g6)
function run_variant(net, sc, pf, var::Symbol, tasks::Int)
  cfg = variant_config(pf, var)
  if var in (:C1, :C2)
    return runInjectionBatch(net, sc; config = cfg, backend = :cpu, start = START, blas_threads = 1, executor = EXECUTOR, parallel_enabled = true, parallel_max_tasks = tasks, parallel_min_work_items = 1)
  end
  CUDA.synchronize()
  r = runInjectionBatch(net, sc; config = cfg, backend = :gpu, start = START, gpu_strategy = gpu_strategy(var), gpu_uniform_mask = var !== :uniform_nomask,
    gpu_two_phase = var in TWO_PHASE ? var : :off, gpu_switch_mismatch = SWITCH, blas_threads = 1, executor = EXECUTOR, parallel_enabled = true, parallel_max_tasks = tasks, parallel_min_work_items = 1)
  CUDA.synchronize()
  return r
end

# the first call, then REPEATS timed calls; keeps the last result and the
# SHA-1 of every result (the matrices of N = 1024 are too large to keep all)
function timed_runs(net, sc, pf, var, tasks)
  times = Float64[]
  hashes = String[]
  stats = Dict{Symbol,Any}[]
  first_s = NaN
  r = nothing
  for k = 0:REPEATS
    GC.gc()
    t0 = time_ns()
    r = run_variant(net, sc, pf, var, tasks)
    t = (time_ns() - t0) / 1.0e9
    push!(hashes, vhash(r))
    k == 0 ? (first_s = t) : (push!(times, t); push!(stats, r.stats))
  end
  return (times = times, first_s = first_s, r = r, hashes = hashes, stats = stats)
end

function compare(r, ref)
  ref === nothing && return (dvm = NaN, dva = NaN, conv_equal = "", iter_equal = "", qviol_equal = "")
  both = r.converged .& ref.converged
  dvm = any(both) ? maximum(maximum(abs.(r.vm[:, s] .- ref.vm[:, s])) for s in findall(both)) : NaN
  dva = any(both) ? maximum(maximum(abs.(r.va_deg[:, s] .- ref.va_deg[:, s])) for s in findall(both)) : NaN
  return (dvm = dvm, dva = dva, conv_equal = r.converged == ref.converged, iter_equal = r.iterations == ref.iterations, qviol_equal = r.qlimit_violations == ref.qlimit_violations)
end

# the first scenarios of the C1 batch against single runpf! solves on fresh
# copies of the batch template with the batch configuration, BLAS at 1
# thread as inside the batch (UMFPACK gives other bits with another BLAS
# thread count)
function single_check(net, sc, ref)
  cfg = Sparlectra._injection_batch_config(net, nothing)
  template, cfg = START === :base ? Sparlectra._injection_warm_template(net, cfg)[1:2] : (net, cfg)
  blas = BLAS.get_num_threads()
  BLAS.set_num_threads(1)
  try
    return _single_check(template, sc, ref, cfg)
  finally
    BLAS.set_num_threads(blas)
  end
end

function _single_check(template, sc, ref, cfg)
  ok = true
  for s = 1:min(SINGLE_CHECKS, length(sc))
    w = deepcopy(template)
    Sparlectra._apply_load_scale!(w, sc.load_scale[:, s])
    runpf!(w; config = cfg)
    vm = [nd._vm_pu for nd in w.nodeVec]
    va = [nd._va_deg for nd in w.nodeVec]
    ok &= vm == ref.vm[:, s] && va == ref.va_deg[:, s] && length(Sparlectra.rectangular_pf_status(w).mismatch_history) == ref.iterations[s]
  end
  return ok
end

const HEADER = "case,nbus,start,workload,N,variant,tasks,status,repeats,first_s,median_s,min_s,max_s,scen_per_s,setup_median_s,worker_copy_median_s,worker_warm_median_s,steady_median_s,gpu_loop_median_s,converged,iter_min,iter_max,final_mismatch_max,qviol_share,qviol_mean,qviol_max,max_dvm_vs_c1,max_dva_deg_vs_c1,conv_equal_c1,iter_equal_c1,qviol_equal_c1,repeat_bitwise,distinct_results,single_solves_bitwise,precision,cudss_mode,reordering_alg,matching_alg,ubatch_mask,analyses,refactorizations,used_gpu_gb,peak_gpu_gb,switch_tol,p1_steps_mean,p1_steps_max,p2_steps_mean,p2_steps_max,switch_by_tol,switch_by_stagnation,phase1_failed,p2_refactorizations,to_float64,phase1_median_s,phase2_median_s,cpu_phase2_median_s,max_abs_y"
const SCEN_HEADER = "case,start,workload,N,variant,scenario,converged,iterations,final_mismatch,qlimit_violations,phase1_steps,phase2_steps,switch_reason,rho"

fmt(x) = x isa AbstractFloat ? @sprintf("%.6g", x) : string(x)

function emit(io, fields...)
  line = join(map(fmt, fields), ",")
  println(line)
  io === nothing || (println(io, line); flush(io))
  flush(stdout)
end

function emit_skip(io, name, n, wl, N, var, why)
  emit(io, name, n, START, wl, N, var, "", why, 0, ntuple(_ -> "", length(split(HEADER, ",")) - 9)...)
end

function emit_scenarios(sio, name, wl, N, var, r)
  sio === nothing && return
  tp = haskey(r.stats, :tp_switch_reason)
  for s = 1:length(r.converged)
    extra = tp ? (r.stats[:tp_phase1_steps][s], r.stats[:tp_phase2_steps][s], r.stats[:tp_switch_reason][s], join([@sprintf("%.3g", x) for x in r.stats[:tp_rho][s]], ";")) : ("", "", "", "")
    println(sio, join(fmt.((name, START, wl, N, var, s, r.converged[s], r.iterations[s], r.final_mismatch[s], r.qlimit_violations[s], extra...)), ","))
  end
  flush(sio)
end

function emit_row(io, sio, name, n, wl, N, var, tasks, tr, ref, single)
  r = tr.r
  cmp = compare(r, ref)
  med = median(tr.times)
  conv = count(r.converged)
  qv = r.qlimit_violations[r.converged]
  st = r.stats
  gpu = r.backend === :gpu
  steady = gpu ? NaN : median([x[:solve_s] - x[:worker_copy_s] - x[:worker_warm_s] for x in tr.stats])
  status = gpu || var in (:C1, :C2) ? "ok" : "cpu fallback: " * replace(r.note, "," => ";")
  emit(io, name, n, START, wl, N, var, gpu ? "" : tasks, status, REPEATS, tr.first_s, med, minimum(tr.times), maximum(tr.times), N / med,
    median([x[:setup_s] for x in tr.stats]), gpu ? NaN : median([x[:worker_copy_s] for x in tr.stats]), gpu ? NaN : median([x[:worker_warm_s] for x in tr.stats]), steady,
    gpu ? median([x[:gpu_loop_s] for x in tr.stats]) : NaN,
    conv, minimum(r.iterations), maximum(r.iterations), conv == 0 ? NaN : maximum(r.final_mismatch[r.converged]),
    conv == 0 ? NaN : count(>(0), qv) / conv, conv == 0 ? NaN : sum(qv) / conv, conv == 0 ? "" : maximum(qv),
    cmp.dvm, cmp.dva, cmp.conv_equal, cmp.iter_equal, cmp.qviol_equal,
    length(unique(tr.hashes)) == 1, length(unique(tr.hashes)), single,
    gpu ? st[:gpu_precision] : "float64", gpu ? st[:gpu_cudss_mode] : "", gpu ? st[:gpu_reordering_alg] : "", gpu ? st[:gpu_matching_alg] : "", gpu ? st[:gpu_uniform_mask] : "",
    gpu ? st[:gpu_analyses] : "", gpu ? st[:gpu_refactorizations] : "",
    gpu ? @sprintf("%.2f", st[:gpu_used_bytes_after_analysis] / 2^30) : "", gpu ? @sprintf("%.2f", maximum(x[:gpu_used_bytes_peak] for x in tr.stats) / 2^30) : "",
    two_phase_fields(r, tr)...)
  emit_scenarios(sio, name, wl, N, var, r)
  @printf(stderr, "  %s %s %-15s t%-3s median %.4f s  %.1f scen/s  conv %d/%d  dvm %.1e  distinct %d\n", Dates.format(now(), "HH:MM:SS"), wl, var, gpu ? "" : string(tasks), med, N / med, conv, N, cmp.dvm, length(unique(tr.hashes)))
end

# the two-phase columns: steps per phase, switch reasons, phase-2
# refactorizations, scenarios that went to Float64, median phase times
function two_phase_fields(r, tr)
  st = r.stats
  haskey(st, :tp_switch_reason) || return ntuple(_ -> "", 14)
  p1 = st[:tp_phase1_steps]
  p2 = st[:tp_phase2_steps]
  why = st[:tp_switch_reason]
  return (st[:gpu_switch_tol], mean(p1), maximum(p1), mean(p2), maximum(p2), count(==(:tol), why), count(==(:stagnation), why), count(in((:failed, :maxiter)), why),
    sum(get(st, :tp_phase2_refactorizations, [0])), count(get(st, :tp_float64, Bool[])),
    median([x[:gpu_phase1_s] for x in tr.stats]), median([x[:gpu_phase2_s] for x in tr.stats]), median([get(x, :cpu_phase2_s, NaN) for x in tr.stats]), st[:gpu_max_abs_y])
end

# largest N in (lo, hi) the variant fits on the device, by bisection on
# single calls, to a resolution of hi/16 (lo fits, hi does not)
function memory_limit(net, pf, var, lo, hi)
  while hi - lo > max(1, hi ÷ 16)
    mid = (lo + hi) ÷ 2
    ok = try
      GC.gc()
      CUDA.reclaim()
      run_variant(net, random_load_scenarios(net, mid; band = 0.2, seed = SEED), pf, var, 1)
      true
    catch e
      oom(e) || rethrow(e)
      false
    end
    ok ? (lo = mid) : (hi = mid)
  end
  return lo, hi
end

# the fastest task count per CPU variant in the scan rows of `path`
function tasks_from(path, name)
  cols = SubString{String}[]
  best = Dict{Symbol,Tuple{Float64,Int}}()
  for line in eachline(path)
    # the header of the file itself (older runs have fewer columns)
    startswith(line, "case,") && (cols = split(line, ","); continue)
    ic(c) = findfirst(==(c), cols)
    f = split(line, ",")
    (length(f) == length(cols) && f[1] == name && f[ic("workload")] == "scan64") || continue
    var = Symbol(f[ic("variant")])
    med = parse(Float64, f[ic("median_s")])
    t = parse(Int, f[ic("tasks")])
    (!haskey(best, var) || med < best[var][1]) && (best[var] = (med, t))
  end
  haskey(best, :C1) && haskey(best, :C2) || error("no scan rows for $(name) in $(path)")
  return Dict(k => v[2] for (k, v) in best)
end

function study_case(io, sio, name)
  p = Sparlectra.ensure_casefile(name * ".m"; to_jl = false)
  t0 = time_ns()
  net = Sparlectra.import_case(p, SparlectraConfig(output = QUIET)).net
  import_s = (time_ns() - t0) / 1.0e9
  n = length(net.nodeVec)
  pf = base_config(net)
  @printf(stderr, "%s %s: %d buses, import %.2f s, start %s\n", Dates.format(now(), "HH:MM:SS"), name, n, import_s, START)

  # CPU task-count scan at N = 64: the fastest count per variant is used for
  # every other N (the N = 64 runs of that count are the random64 rows)
  scan_sc = random_load_scenarios(net, SCAN_N; band = 0.2, seed = SEED)
  best = isempty(TASKS_FROM) ? Dict{Symbol,Int}() : tasks_from(TASKS_FROM, name)
  scan = Dict{Symbol,Any}()
  for var in (isempty(TASKS_FROM) ? (:C1, :C2) : ())
    for t in filter(<=(Threads.nthreads()), TASKS)
      tr = timed_runs(net, scan_sc, pf, var, t)
      emit_row(io, nothing, name, n, "scan64", SCAN_N, var, t, tr, nothing, "")
      if !haskey(scan, var) || median(tr.times) < median(scan[var].times)
        scan[var] = tr
        best[var] = t
      end
    end
  end
  @printf(stderr, "  task counts: C1 %d, C2 %d\n", best[:C1], best[:C2])

  workloads = [(string("random", N), N) for N in NS]
  get(ENV, "RAMP", "1") == "1" && push!(workloads, ("ramp96", 96))
  limit = Dict{Symbol,Int}()
  for (wl, N) in workloads
    sc = startswith(wl, "ramp") ? load_ramp_scenarios(net, N; low = 0.8, high = 1.2) : random_load_scenarios(net, N; band = 0.2, seed = SEED)
    if haskey(limit, :sequential) && N > limit[:sequential]
      emit_skip(io, name, n, wl, N, "all", string("skipped: above the GPU memory limit of the default (", limit[:sequential], ")"))
      continue
    end
    # the GPU default first (its memory limit decides whether the CPU runs)
    pending = Any[]
    gpu_vars = STRATEGIES
    for var in gpu_vars
      if haskey(limit, var) && N > limit[var]
        emit_skip(io, name, n, wl, N, var, string("skipped: above its GPU memory limit (", limit[var], ")"))
        continue
      end
      tr = try
        timed_runs(net, sc, pf, var, best[:C1])
      catch e
        gpu_refused(e) || rethrow(e)
        msg = replace(first(sprint(showerror, e), 120), "," => ";", "\n" => " ")
        emit_skip(io, name, n, wl, N, var, "refused: " * msg)
        @printf(stderr, "  %s %s refused: %s\n", wl, var, msg)
        if oom(e)
          lower = maximum(filter(<(N), NS); init = 0)
          lo, hi = memory_limit(net, pf, var, lower, N)
          limit[var] = lo
          emit_skip(io, name, n, wl, N, var, string("memory limit: fits N = ", lo, "; out of memory at N = ", hi))
          @printf(stderr, "  %s memory limit: fits %d, fails %d\n", var, lo, hi)
        end
        continue
      end
      push!(pending, (var, tr))
    end
    if haskey(limit, :sequential) && N > limit[:sequential]
      for (var, tr) in pending
        emit_row(io, sio, name, n, wl, N, var, 0, tr, nothing, "")
      end
      continue
    end
    ref = nothing
    if !CPU_TIMED
      # the reference only: C1 once, not timed and not written
      ref = run_variant(net, sc, pf, :C1, best[:C1])
    end
    for var in (CPU_TIMED ? (:C1, :C2) : ())
      tr = haskey(scan, var) && wl == "random$(SCAN_N)" ? scan[var] : timed_runs(net, sc, pf, var, best[var])
      var === :C1 && (ref = tr.r)
      single = var === :C1 ? string(single_check(net, sc, tr.r)) : ""
      emit_row(io, sio, name, n, wl, N, var, best[var], tr, ref, single)
    end
    for (var, tr) in pending
      emit_row(io, sio, name, n, wl, N, var, 0, tr, ref, "")
    end
    pending = nothing
    GC.gc()
    CUDA.reclaim()
  end
  net = nothing
  GC.gc()
  return nothing
end

function main()
  cases = isempty(ARGS) ? CASES : ARGS
  commit = try
    readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`)
  catch e
    e isa ProcessFailedException || e isa Base.IOError || rethrow(e)
    "unknown"
  end
  out = get(ENV, "OUT", "")
  out_scen = get(ENV, "OUT_SCEN", "")
  io = isempty(out) ? nothing : open(out, "a")
  sio = isempty(out_scen) ? nothing : open(out_scen, "a")
  head = string("# batch_study_gpu started ", now(), " commit ", commit, " Julia ", VERSION, " threads ", Threads.nthreads(), " BLAS threads outside the batch ", BLAS.get_num_threads(), " (1 inside the CPU batch) CUDA.jl ", pkgversion(CUDA), " CUDSS.jl ", pkgversion(CUDSS), " libcudss ", CUDSS.version(), " device ", CUDA.name(CUDA.device()), " repeats ", REPEATS, " after one first call, seed ", SEED, " start ", START, " cases ", join(cases, " "))
  println(head)
  println(HEADER)
  if io !== nothing
    println(io, head)
    println(io, HEADER)
  end
  if sio !== nothing
    println(sio, head)
    println(sio, SCEN_HEADER)
  end
  try
    for name in cases
      study_case(io, sio, name)
    end
  finally
    io === nothing || close(io)
    sio === nothing || close(sio)
  end
  println("# finished ", now())
end

Base.invokelatest(main)

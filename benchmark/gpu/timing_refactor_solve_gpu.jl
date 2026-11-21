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
# file: benchmark/gpu/timing_refactor_solve_gpu.jl
# purpose: GPU study, stage 1b: the time of one cuDSS refactorization plus
#          solve, deterministic and fast mode, Float64 and Float32, next to
#          the end-to-end time of the same solves and of C1, so the per-step
#          numbers can be checked against the total.
#          1. Isolated: the device Jacobian at the C1 solution, one solver
#             per mode and precision; MICRO repetitions of refactorization
#             plus solve, each with CUDA.synchronize() before and after on
#             the host clock (and refactorization and solve alone).
#          2. In the solver: REPEATS device Newton solves (G3, G4) after one
#             warm-up, each wrapped in CUDA.synchronize() on the host clock
#             (end to end, run_sparlectra), with the loop's own phase times
#             (CUDA events) from the solver status; C1 the same way.
#          Medians throughout. The phase times of the loop must add up to
#          its core time up to the host work between them; the rest of the
#          end-to-end time is import, Ybus and result work outside the loop.
#
# Run: julia --project=benchmark/gpu benchmark/gpu/timing_refactor_solve_gpu.jl [case ...]
#   default cases case_ACTIVSg25k case_ACTIVSg70k; REPEATS (default 5),
#   MICRO (default 20). CSV: $GPU_STUDY_RUNS/timing_refactor_solve_gpu.csv
#   (default benchmark/gpu/results).

using Sparlectra, CUDA, CUDSS, SparseArrays, LinearAlgebra, Statistics, Printf, Dates, Logging

const EXT = Base.get_extension(Sparlectra, :SparlectraCUDAExt)
EXT === nothing && error("SparlectraCUDAExt not loaded")

const QUIET = OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)
const REPEATS = parse(Int, get(ENV, "REPEATS", "5"))
const MICRO = parse(Int, get(ENV, "MICRO", "20"))

cfg(cb) = SparlectraConfig(powerflow = PowerFlowConfig(compute_backend = cb, max_iter = 30), output = QUIET)
bus_type_symbol(t) = t == Sparlectra.Slack ? :Slack : t == Sparlectra.PV ? :PV : :PQ

function out_dir()
  d = get(ENV, "GPU_STUDY_RUNS", "")
  isempty(d) || return d
  d = joinpath(@__DIR__, "results")
  mkpath(d)
  return d
end

# seconds of f() with the device idle before and finished after
function synced(f)
  CUDA.synchronize()
  t0 = time_ns()
  f()
  CUDA.synchronize()
  return (time_ns() - t0) / 1.0e9
end

# 1. refactorization plus solve on the device Jacobian at the solution
function isolated(net, T, fast::Bool)
  n = length(net.nodeVec)
  Ybus = Sparlectra.createYBUS(net = net, sparse = true)
  _, slack = Sparlectra.initialVrect(net; flatstart = true)
  S = Sparlectra.buildComplexSVec(net)
  Vset = Sparlectra._bus_voltage_setpoints_from_prosumers(net)
  V = Sparlectra.buildVoltageVector(net)
  types = [bus_type_symbol(Sparlectra.getNodeType(nd)) for nd in net.nodeVec]
  g = EXT._device_jacobian_mismatch(Ybus, slack, reshape(V, n, 1), [types], reshape(S, n, 1), reshape(Vset, n, 1), zeros(n, 1), zeros(n, 1); T = T)
  nvar = size(g.F, 1)
  nz = CuVector{T}(vec(g.J[:, 1]))
  solver = CudssSolver(CuVector{Int32}(Int32.(g.jrowptr)), CuVector{Int32}(Int32.(g.jcol)), nz, "G", 'F')
  fast || cudss_set(solver, "deterministic_mode", 1)
  cudss_set(solver, "pivot_epsilon", zero(T))
  x = CUDA.zeros(T, nvar)
  b = CuVector{T}(T.(-g.F[:, 1]))
  X = CudssMatrix(x)
  B = CudssMatrix(b)
  cudss("analysis", solver, X, B)
  cudss("factorization", solver, X, B)
  CUDA.synchronize()
  both = Float64[]
  refac = Float64[]
  solve = Float64[]
  for _ = 1:(MICRO+2)
    push!(both, synced(() -> (cudss("refactorization", solver, X, B); cudss("solve", solver, X, B))))
    push!(refac, synced(() -> cudss("refactorization", solver, X, B)))
    push!(solve, synced(() -> cudss("solve", solver, X, B)))
  end
  # the first two of each are warm-up
  return (nvar = nvar, nnzJ = length(nz), both = median(both[3:end]), refac = median(refac[3:end]), solve = median(solve[3:end]))
end

# 2. one solve end to end, host clock around synchronized work
function one_solve(p, cb, device::Bool)
  r = nothing
  t = synced() do
    r = with_logger(NullLogger()) do
      device ? Sparlectra._with_device_newton(() -> run_sparlectra(casefile = p, config = cfg(cb))) : run_sparlectra(casefile = p, config = cfg(cb))
    end
  end
  return r, t
end

const PHASES = (:gpu_newton_build_s, :gpu_newton_mismatch_s, :gpu_newton_refactor_s, :gpu_newton_solve_s, :gpu_newton_update_s, :gpu_newton_analysis_s, :gpu_newton_factor_first_s, :gpu_newton_transfer_s, :gpu_newton_host_qlimit_s)

function main()
  cases = isempty(ARGS) ? ["case_ACTIVSg25k", "case_ACTIVSg70k"] : ARGS
  commit = try
    readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`)
  catch e
    e isa ProcessFailedException || e isa Base.IOError || rethrow(e)
    "unknown"
  end
  started = now()
  println("timing_refactor_solve_gpu  started ", started, "  commit ", commit, "  Julia ", VERSION, "  threads ", Threads.nthreads(), "  BLAS threads ", BLAS.get_num_threads(), "  CUDA.jl ", pkgversion(CUDA), "  CUDSS.jl ", pkgversion(CUDSS), "  device ", CUDA.name(CUDA.device()), "  repeats ", REPEATS, " micro ", MICRO, ", medians")
  csv = String["case,variant,mode,iso_refac_solve_ms,iso_refac_ms,iso_solve_ms,loop_refac_plus_solve_ms,refactors,solves,loop_core_s,loop_core_min_s,loop_core_max_s,analysis_s,phase_sum_s,end_to_end_s,c1_end_to_end_s"]
  for name in cases
    p = Sparlectra.ensure_casefile(name * ".m"; to_jl = false)
    # C1: warm-up, then REPEATS
    rc, _ = one_solve(p, Sparlectra.ComputeBackendConfig(), false)
    c1_t = median([one_solve(p, Sparlectra.ComputeBackendConfig(), false)[2] for _ = 1:REPEATS])
    @printf("\n%s: n %d, C1 end to end %.3f s (converged %s, %d entries)\n", name, length(rc.net.nodeVec), c1_t, rc.final_converged, length(Sparlectra.rectangular_pf_status(rc.net).mismatch_history))
    for (var, prec) in (("G3", :float64), ("G4", :float32_mixed)), fast in (false, true)
      T = prec === :float64 ? Float64 : Float32
      iso = isolated(rc.net, T, fast)
      cb = Sparlectra.ComputeBackendConfig(mode = :gpu, gpu_precision = prec, gpu_fast_mode = fast)
      one_solve(p, cb, true)
      runs = [one_solve(p, cb, true) for _ = 1:REPEATS]
      sts = [Sparlectra.rectangular_pf_status(r.net) for (r, _) in runs]
      all(st -> get(st, :gpu_newton, false) && st.linear_solver_fallback_count == 0, sts) || @warn "a run did not stay on the device loop" name var
      e2e = median([t for (_, t) in runs])
      refs = sts[1].linear_solver_refactor_count
      steps = get(sts[1], :jacobian_reuse_refactorisations, length(sts[1].mismatch_history) - 1)
      # one refactorization plus one solve in the loop, as the isolated
      # number: the loop's refactorization time over its refactorizations
      # plus its solve time over its solves (one solve per Newton step; the
      # first step uses the factorization of the pre-flight, no refactorization)
      per_step = median([st.gpu_newton_refactor_s / max(refs, 1) + st.gpu_newton_solve_s / max(steps, 1) for st in sts]) * 1e3
      cores = [st.gpu_newton_core_s for st in sts]
      core = median(cores)
      ana = median([st.gpu_newton_analysis_s for st in sts])
      psum = median([sum(getproperty(st, ph) for ph in PHASES) for st in sts])
      mode = fast ? "fast" : "deterministic"
      @printf("  %-3s %-13s | isolated, synchronized: refac+solve %7.3f ms (refac %7.3f, solve %7.3f) | in the loop, events: refac+solve %7.3f ms (%d refactorizations, %d solves) | loop core %.3f s (min %.3f, max %.3f), analysis %.3f s, phases %.3f s | end to end %.3f s (C1 %.3f s)\n",
        var, mode, iso.both * 1e3, iso.refac * 1e3, iso.solve * 1e3, per_step, refs, steps, core, minimum(cores), maximum(cores), ana, psum, e2e, c1_t)
      println("      phases (median of the runs, s): ", join([@sprintf("%s %.4f", replace(String(ph), "gpu_newton_" => "", "_s" => ""), median([getproperty(st, ph) for st in sts])) for ph in PHASES], ", "))
      push!(csv, @sprintf("%s,%s,%s,%.6f,%.6f,%.6f,%.6f,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f", name, var, mode, iso.both * 1e3, iso.refac * 1e3, iso.solve * 1e3, per_step, refs, steps, core, minimum(cores), maximum(cores), ana, psum, e2e, c1_t))
      flush(stdout)
    end
  end
  path = joinpath(out_dir(), "timing_refactor_solve_gpu.csv")
  open(path, "w") do io
    println(io, "# timing_refactor_solve_gpu.jl commit ", commit, " started ", started, " finished ", now(), " threads ", Threads.nthreads(), " repeats ", REPEATS, " micro ", MICRO)
    foreach(l -> println(io, l), csv)
  end
  println("CSV ", path, "\nfinished ", now())
end

Base.invokelatest(main)

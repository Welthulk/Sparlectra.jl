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
# file: benchmark/gpu/history_bound_gpu.jl
# purpose: GPU study, stage 1b: the data behind the bounds of the mismatch
#          history comparison in check_device_newton_gpu.jl. Per case the C1
#          history, cond(J, 1) of the Jacobian at the solution (1-norm times
#          the estimated 1-norm of the inverse), and the deviation of every
#          later history entry of G1, G3 (Float64) and G2, G4 (float32_mixed)
#          from C1, in the deterministic mode (default) and in the fast mode.
#
#          Model: a backward stable solve of J dx = F in precision u gives
#          a step whose residual is at most about u ||J|| ||dx||, so the
#          next mismatch moves by about cond(J, 1) u ||F_{k-1}||. Below that
#          sits the rounding of the mismatch itself, eps(Float64) times
#          sigma = max_i |V_i| sum_j |Y_ij| |V_j| at the solution. The
#          script reports per entry the ratio
#          dev_k / (cond(J, 1) u c_{k-1} + eps(Float64) sigma)
#          with c the C1 history and u = eps(Float64) or eps(Float32); a
#          ratio up to about 10 fits the model.
#
#          Rule for the bounds, first form (written before the first data
#          run): Float64: HIST_REL = 10 times the largest |g_k - c_k| / |c_k|
#          over the entries with |c_k| >= 1e-6, HIST_ABS = 10 times the
#          largest |g_k - c_k| over the entries with |c_k| < 1e-6 (the
#          rounding floor), each rounded up to a power of ten, over all
#          cases, both variants and both modes. float32_mixed, intermediate
#          entries (k = 2 .. L-1, L the shorter history): HIST32 = 10 times
#          the largest |g_k - c_k| / c_{k-1}, rounded up to a power of ten.
#          The last entry of a float32_mixed run is judged by the tolerance,
#          not by C1.
#          Second form (set after the data was read, on review): the check
#          adds HIST_ABS and HIST_REL |c_k|, so HIST_REL is the smallest
#          power of ten with 10 |g_k - c_k| <= HIST_ABS + HIST_REL |c_k| on
#          every entry, with HIST_ABS as above; float32_mixed also gets
#          HIST32_CUR = 10 times the largest |g_k - c_k| / c_k, rounded up;
#          all over the deterministic runs, the mode the check runs in. The
#          fast runs are reported against the same constants.
#
# Run: julia --project=benchmark/gpu benchmark/gpu/history_bound_gpu.jl [case ...]
#   RUNS_DET (default 20) and RUNS_FAST (default 10) runs per variant and
#   mode. CSV: $GPU_STUDY_RUNS/history_bound_gpu.csv (default
#   benchmark/gpu/results).

using Sparlectra, CUDA, CUDSS, SparseArrays, LinearAlgebra, Printf, Dates, Logging

const EXT = Base.get_extension(Sparlectra, :SparlectraCUDAExt)
EXT === nothing && error("SparlectraCUDAExt not loaded")

const STUDY_CASES = ["case14", "case118", "case1354pegase", "case2869pegase", "case6495rte", "case9241pegase", "case_ACTIVSg10k", "case_ACTIVSg25k", "case_ACTIVSg70k"]
const QUIET = OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)
const RUNS_DET = parse(Int, get(ENV, "RUNS_DET", "20"))
const RUNS_FAST = parse(Int, get(ENV, "RUNS_FAST", "10"))
const SPLIT = 1e-6

cfg(cb) = SparlectraConfig(powerflow = PowerFlowConfig(compute_backend = cb, max_iter = 30), output = QUIET)
gpu_cb(prec, fast) = Sparlectra.ComputeBackendConfig(mode = :gpu, gpu_precision = prec, gpu_fast_mode = fast)
hist(r) = copy(Sparlectra.rectangular_pf_status(r.net).mismatch_history)
bus_type_symbol(t) = t == Sparlectra.Slack ? :Slack : t == Sparlectra.PV ? :PV : :PQ

function out_dir()
  d = get(ENV, "GPU_STUDY_RUNS", "")
  isempty(d) || return d
  d = joinpath(@__DIR__, "results")
  mkpath(d)
  return d
end

# cond(J, 1) of the device Jacobian (Float64) at the C1 solution with the
# final bus types
function cond1_at_solution(net)
  n = length(net.nodeVec)
  Ybus = Sparlectra.createYBUS(net = net, sparse = true)
  _, slack = Sparlectra.initialVrect(net; flatstart = true)
  S = Sparlectra.buildComplexSVec(net)
  Vset = Sparlectra._bus_voltage_setpoints_from_prosumers(net)
  V = Sparlectra.buildVoltageVector(net)
  types = [bus_type_symbol(Sparlectra.getNodeType(nd)) for nd in net.nodeVec]
  g = EXT._device_jacobian_mismatch(Ybus, slack, reshape(V, n, 1), [types], reshape(S, n, 1), reshape(Vset, n, 1), zeros(n, 1), zeros(n, 1); T = Float64)
  nvar = size(g.F, 1)
  J = SparseMatrixCSC(transpose(SparseMatrixCSC(nvar, nvar, Int.(g.jrowptr), Int.(g.jcol), vec(g.J[:, 1]))))
  sigma = maximum(abs.(V) .* (abs.(Ybus) * abs.(V)))
  return opnorm(J, 1) * SparseArrays.opnormestinv(J), sigma
end

function run_variant(p, var::Symbol, fast::Bool)
  prec = var in (:G2, :G4) ? :float32_mixed : :float64
  device = var in (:G3, :G4)
  r = with_logger(NullLogger()) do
    device ? Sparlectra._with_device_newton(() -> run_sparlectra(casefile = p, config = cfg(gpu_cb(prec, fast)))) : run_sparlectra(casefile = p, config = cfg(gpu_cb(prec, fast)))
  end
  return (hist = hist(r), converged = r.final_converged, mode = Sparlectra.rectangular_pf_status(r.net).gpu_cudss_mode)
end

# deviations of one GPU history g from the C1 history c
function deviations(g, c, cond1, sigma, f32::Bool)
  u = f32 ? Float64(eps(Float32)) : eps(Float64)
  L = min(length(g), length(c))
  last = f32 ? L - 1 : L
  rows = NamedTuple[]
  for k = 2:last
    d = abs(g[k] - c[k])
    push!(rows, (k = k, c = c[k], d = d, rel = d / abs(c[k]), rel_prev = d / c[k-1], model = d / (cond1 * u * c[k-1] + eps(Float64) * sigma)))
  end
  return (first_bitwise = g[1] === c[1], same_length = length(g) == length(c), rows = rows)
end

pow10up(x) = x <= 0 ? 0.0 : 10.0^ceil(log10(x))

function main()
  cases = isempty(ARGS) ? STUDY_CASES : ARGS
  commit = try
    readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`)
  catch e
    e isa ProcessFailedException || e isa Base.IOError || rethrow(e)
    "unknown"
  end
  started = now()
  println("history_bound_gpu  started ", started, "  commit ", commit, "  Julia ", VERSION, "  CUDA.jl ", pkgversion(CUDA), "  CUDSS.jl ", pkgversion(CUDSS), "  device ", CUDA.name(CUDA.device()), "  runs deterministic ", RUNS_DET, " fast ", RUNS_FAST)
  csv = String["case,var,mode,run,k,c1,dev,rel,rel_prev,model_ratio,cond1,same_length,first_bitwise"]
  # global maxima for the rule: (value, case, var, mode, run, k)
  best_rel = (0.0, "", :none, :none, 0, 0)
  best_abs = (0.0, "", :none, :none, 0, 0)
  best_32 = (0.0, "", :none, :none, 0, 0)
  cond_max = 0.0
  for name in cases
    p = Sparlectra.ensure_casefile(name * ".m"; to_jl = false)
    rc = with_logger(() -> run_sparlectra(casefile = p, config = cfg(Sparlectra.ComputeBackendConfig())), NullLogger())
    c = hist(rc)
    cond1, sigma = cond1_at_solution(rc.net)
    cond_max = max(cond_max, cond1)
    @printf("\n%s: C1 %d entries, last %.2e, cond(J,1) %.2e, cond*eps64 %.2e, cond*eps32 %.2e, mismatch rounding eps*sigma %.2e\n", name, length(c), c[end], cond1, cond1 * eps(Float64), cond1 * eps(Float32), eps(Float64) * sigma)
    println("  C1 history: ", join([@sprintf("%.3e", x) for x in c], " "))
    for var in (:G1, :G3, :G2, :G4), fast in (false, true)
      f32 = var in (:G2, :G4)
      nruns = fast ? RUNS_FAST : RUNS_DET
      # per run: the largest deviations; per variant and mode: their maxima
      worst_rel = (0.0, 0, 0)
      worst_abs = (0.0, 0, 0)
      worst_32 = (0.0, 0, 0)
      worst_model = 0.0
      lengths = Int[]
      firsts = Bool[]
      patterns = Set{Vector{UInt64}}()
      for run = 1:nruns
        g = run_variant(p, var, fast)
        push!(patterns, reinterpret(UInt64, g.hist))
        dv = deviations(g.hist, c, cond1, sigma, f32)
        push!(lengths, length(g.hist))
        push!(firsts, dv.first_bitwise)
        for row in dv.rows
          push!(csv, @sprintf("%s,%s,%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%s,%s", name, var, fast ? "fast" : "deterministic", run, row.k, row.c, row.d, row.rel, row.rel_prev, row.model, cond1, dv.same_length, dv.first_bitwise))
          worst_model = max(worst_model, row.model)
          if f32
            row.rel_prev > worst_32[1] && (worst_32 = (row.rel_prev, run, row.k))
          elseif row.c >= SPLIT
            row.rel > worst_rel[1] && (worst_rel = (row.rel, run, row.k))
          else
            row.d > worst_abs[1] && (worst_abs = (row.d, run, row.k))
          end
        end
      end
      mode = fast ? :fast : :deterministic
      if f32
        worst_32[1] > best_32[1] && (best_32 = (worst_32[1], name, var, mode, worst_32[2], worst_32[3]))
        @printf("  %-3s %-13s runs %2d, distinct histories %2d, lengths %s, entry 1 bitwise %s | intermediate max |g-c|/c_{k-1} %.2e (run %d, k %d) | max model ratio %.2e\n", var, mode, nruns, length(patterns), string(sort(unique(lengths))), all(firsts), worst_32..., worst_model)
      else
        worst_rel[1] > best_rel[1] && (best_rel = (worst_rel[1], name, var, mode, worst_rel[2], worst_rel[3]))
        worst_abs[1] > best_abs[1] && (best_abs = (worst_abs[1], name, var, mode, worst_abs[2], worst_abs[3]))
        @printf("  %-3s %-13s runs %2d, distinct histories %2d, lengths %s, entry 1 bitwise %s | max rel (|c| >= 1e-6) %.2e (run %d, k %d) | max abs (|c| < 1e-6) %.2e (run %d, k %d) | max model ratio %.2e\n", var, mode, nruns, length(patterns), string(sort(unique(lengths))), all(firsts), worst_rel..., worst_abs..., worst_model)
      end
      flush(stdout)
    end
  end
  path = joinpath(out_dir(), "history_bound_gpu.csv")
  open(path, "w") do io
    println(io, "# history_bound_gpu.jl commit ", commit, " started ", started, " finished ", now())
    foreach(l -> println(io, l), csv)
  end
  println("\nlargest Float64 relative deviation (|c| >= 1e-6): ", best_rel)
  println("largest Float64 absolute deviation (|c| < 1e-6): ", best_abs)
  println("largest float32_mixed intermediate |g-c|/c_{k-1}: ", best_32)
  @printf("rule: HIST_REL = %.0e, HIST_ABS = %.0e, HIST32 = %.0e; model bound for float32_mixed max cond(J,1) eps32 = %.2e (%s)\n", pow10up(10 * best_rel[1]), pow10up(10 * best_abs[1]), pow10up(10 * best_32[1]), cond_max * eps(Float32), pow10up(10 * best_32[1]) <= cond_max * eps(Float32) ? "HIST32 below it" : "HIST32 ABOVE it")
  println("CSV ", path, "\nfinished ", now())
end

Base.invokelatest(main)

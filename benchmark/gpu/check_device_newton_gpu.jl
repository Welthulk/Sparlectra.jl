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
# file: benchmark/gpu/check_device_newton_gpu.jl
# purpose: check of the device Newton loop (GPU study, stage 1b) as a solver:
#          1. cuDSS marks a singular matrix only through "npivots" or a
#             non-finite solution, never through "info"; with pivot_epsilon
#             0 (the device loop) always through a non-finite solution;
#          2. G1, G2 (stage 1) and G3, G4 (device loop) against C1 on the
#             study cases: convergence, mismatch history (length and
#             values), Newton steps, Q-limit counters and log, final active
#             set, max |dV|, fallbacks, host copies per iteration. Float32
#             runs may take other step counts (reported); everything else
#             must match. Per case also G3 against C1 with Q limits
#             without cooldown (a bus can be clamped and released in one
#             iteration);
#          3. every recovery path of the device loop by injected faults
#             (GPU errors, non-finite steps, a failed step check, a
#             singular CPU redo, fallback off), against C1;
#          4. a kernel fault that kills the CUDA context, in a child process
#             per fault point: nothing can be read from the device, the
#             CPU runs the whole solve again from its start and must give
#             the path and the result of C1.
#          Exit code 0 only if every check passed.
#
# Run: julia --project=benchmark/gpu benchmark/gpu/check_device_newton_gpu.jl [case ...]
#   CASES default: the nine cases of the study. FAULT_CASE (default
#   case2869pegase) and VDEP_CASE (default case6495rte) for parts 3 and 4.
#   PARTS (default 1,2,3,4) selects the parts.
#   Internal: STICKY=<point>:<iteration> runs one part-4 child.

using Sparlectra, CUDA, CUDSS, SparseArrays, LinearAlgebra, Printf, Dates, Logging
using CUDA.CUSPARSE: CuSparseMatrixCSR

const EXT = Base.get_extension(Sparlectra, :SparlectraCUDAExt)
EXT === nothing && error("SparlectraCUDAExt not loaded")

const STUDY_CASES = ["case14", "case118", "case1354pegase", "case2869pegase", "case6495rte", "case9241pegase", "case_ACTIVSg10k", "case_ACTIVSg25k", "case_ACTIVSg70k"]
const QUIET = OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)
# G3/G4 against C1: the update kernel is not bitwise (design 2.4), so V
# differs at rounding level; this bound only catches a wrong solution
const DV_BOUND = 1e-6

# rescue off in parts 3 and 4: a non-converged solve (singular step, max_iter)
# would otherwise be retried with other start values and report the retry
cfg(cb; max_iter = 30, rescue = true, qlimits = QLimitConfig()) = SparlectraConfig(powerflow = PowerFlowConfig(compute_backend = cb, max_iter = max_iter, rescue = rescue, qlimits = qlimits), output = QUIET)
gpu_cb(prec; fallback = true) = Sparlectra.ComputeBackendConfig(mode = :gpu, gpu_precision = prec, gpu_fallback_to_cpu = fallback)

casefile(name) = Sparlectra.ensure_casefile(name * ".m"; to_jl = false)

# the final active set: node types and the Q-limit clamps (bus => side)
active_set(net) = (types = [Sparlectra.getNodeType(nd) for nd in net.nodeVec], clamps = sort!(collect(net.qLimitEvents)))

# what is compared between two runs
function summary_of(r)
  st = Sparlectra.rectangular_pf_status(r.net)
  return (
    converged = r.final_converged,
    reason = st.reason,
    history = length(st.mismatch_history),
    hist = copy(st.mismatch_history),
    qlog = [(e.iter, e.bus, e.side, e.kind) for e in r.net.qLimitLog],
    steps = get(st, :jacobian_reuse_refactorisations, -1),
    qchanges = st.qlimit_active_set_changes,
    qreenable = st.qlimit_reenable_events,
    qswitch_its = st.type_switch_iterations,
    V = Sparlectra.buildVoltageVector(r.net),
    active = active_set(r.net),
    st = st,
  )
end

# The mismatch history by value. Entry 1 belongs to the start iterate both
# runs share, so it must be bitwise equal (the device mismatch repeats the
# CPU arithmetic, check_jacobian_gpu.jl). The later iterates differ at
# rounding level (cuDSS against UMFPACK, the update kernel), so a later entry
# passes if |gpu - cpu| <= HIST_ABS + HIST_REL * |cpu|. Both constants come
# from history_bound_gpu.jl by the rule in its header, in this combined form
# and over the deterministic runs (the mode the check runs in): HIST_ABS is
# 10 times the largest deviation on the entries below 1e-6 pu (3.5e-12,
# case_ACTIVSg70k), rounded up to a power of ten; HIST_REL the smallest
# power of ten with 10 |gpu - cpu| <= HIST_ABS + HIST_REL |cpu| on every
# entry (largest need 2.5e-11, case_ACTIVSg70k entry 3, where |cpu| is 6.6).
# The fast-mode runs fit the same constants (need 2.9e-11). HIST_ABS stays
# two orders below the tolerance 1e-8, so both runs take the same
# convergence decisions; a wrong step or another active set changes an
# entry by many orders more. `ratio` is the largest |gpu - cpu| / bound (1
# passes).
const HIST_ABS = 1e-10
const HIST_REL = 1e-10
function hist_check(g, c)
  length(g.hist) == length(c.hist) || return (ok = false, first_bitwise = false, ratio = Inf)
  first_bitwise = isempty(c.hist) || g.hist[1] === c.hist[1]
  ratio = maximum((abs(g.hist[k] - c.hist[k]) / (HIST_ABS + HIST_REL * abs(c.hist[k])) for k in eachindex(c.hist)); init = 0.0)
  return (ok = first_bitwise && ratio <= 1.0, first_bitwise = first_bitwise, ratio = ratio)
end

# float32_mixed: the Float32 factorization makes each step an inexact Newton
# step. The step differs from the Float64 step by an error proportional to
# the step, so the next mismatch differs from C1 by an amount proportional
# to the previous mismatch: an intermediate entry k = 2 .. L-1 (L the
# shorter history) passes if |gpu - cpu| <= HIST32 * cpu_{k-1}. Near
# quadratic convergence cpu_k is far below cpu_{k-1}, and that bound alone
# would allow a deviation far above cpu_k; so the entry must also stay
# within HIST32_CUR * cpu_k. The last entry is judged by the tolerance
# (convergence), and a run may take another step count (task section 9), so
# the entries behind L are not compared. Entry 1 must be bitwise equal as
# for Float64. Both constants by the same rule over the deterministic runs
# (history_bound_gpu.jl): largest |gpu - cpu| / cpu_{k-1} 1.9e-4 (G4,
# case_ACTIVSg70k), largest |gpu - cpu| / cpu_k 0.11 (G4, case9241pegase,
# entry 7 at 1.2e-8 pu). The fast mode reaches 7.2e-3 and 0.12, inside both
# bounds but without the factor 10. cond(J, 1) eps(Float32) is above 1 from
# case2869pegase on, so the bounds come from the data, not from that model.
const HIST32 = 1e-2
const HIST32_CUR = 10.0
function hist_check32(g, c)
  L = min(length(g.hist), length(c.hist))
  L == 0 && return (ok = false, first_bitwise = false, ratio = Inf)
  first_bitwise = g.hist[1] === c.hist[1]
  bound(k) = min(HIST32 * c.hist[k-1], HIST32_CUR * c.hist[k])
  ratio = maximum((abs(g.hist[k] - c.hist[k]) / bound(k) for k = 2:(L-1)); init = 0.0)
  return (ok = first_bitwise && ratio <= 1.0, first_bitwise = first_bitwise, ratio = ratio)
end

run_c1(p; max_iter = 30, rescue = true, qlimits = QLimitConfig()) = summary_of(run_sparlectra(casefile = p, config = cfg(Sparlectra.ComputeBackendConfig(); max_iter = max_iter, rescue = rescue, qlimits = qlimits)))
run_stage1(p, prec) = summary_of(run_sparlectra(casefile = p, config = cfg(gpu_cb(prec))))
run_device(p, prec; fallback = true, max_iter = 30, rescue = true, qlimits = QLimitConfig()) = summary_of(Sparlectra._with_device_newton(() -> run_sparlectra(casefile = p, config = cfg(gpu_cb(prec; fallback = fallback); max_iter = max_iter, rescue = rescue, qlimits = qlimits))))

const FAILS = String[]
check!(ok::Bool, what) = (ok || push!(FAILS, what); ok)

# --- 1. cuDSS on singular matrices -------------------------------------------

function cudss_pivots(A::SparseMatrixCSC{Float64}, T::Type; pivot_epsilon = nothing, deterministic = false)
  At = copy(transpose(A))
  d = CuSparseMatrixCSR{T,Int32}(CuVector{Int32}(Int32.(At.colptr)), CuVector{Int32}(Int32.(At.rowval)), CuVector{T}(T.(At.nzval)), size(A))
  n = size(A, 1)
  x = CUDA.zeros(T, n)
  b = CuVector{T}(T.(collect(1.0:n)))
  s = CudssSolver(d, "G", 'F')
  pivot_epsilon === nothing || cudss_set(s, "pivot_epsilon", T(pivot_epsilon))
  deterministic && cudss_set(s, "deterministic_mode", 1)
  cudss("analysis", s, x, b)
  cudss("factorization", s, x, b)
  cudss("solve", s, x, b)
  xs = Array(x)
  return (info = Int(cudss_get(s, "info")), npivots = Int(cudss_get(s, "npivots")), maxx = Float64(maximum(abs.(xs))), finite = all(isfinite, xs), x = xs)
end

function part1()
  println("\n== 1. cuDSS: info and npivots on regular and singular matrices, default pivot_epsilon (stage 1) and 0 (device loop)")
  n = 200
  A = spdiagm(-1 => fill(-1.0, n - 1), 0 => fill(4.0, n), 1 => fill(-1.0, n - 1))
  Z = copy(A)
  Z[5, :] .= 0.0
  # the same row with only a stored zero on the diagonal
  Zd = dropzeros(Z)
  Zd[5, 5] = 0.0
  E = copy(A)
  E[7, :] = E[8, :]
  W = sparse([0.0 1.0; 1.0 0.0])
  for (label, M, singular) in (("regular tridiagonal", A, false), ("zero row, stored zeros", Z, true), ("zero row, zero diagonal", Zd, true), ("two equal rows", E, true), ("2x2 swap, zero diagonal", W, false))
    for T in (Float64, Float32)
      r = cudss_pivots(M, T)
      r0 = cudss_pivots(M, T; pivot_epsilon = 0.0)
      # default: a singular matrix shows as a perturbed pivot or a non-finite
      # solution (the stage 1 checks), info is never set; pivot_epsilon 0: a
      # singular matrix gives a non-finite solution (the device loop check)
      # and a regular one the bits of the default (compared in deterministic
      # mode: without it two solves of one matrix differ in the last bits)
      ok1 = r.info == 0 && (singular ? (r.npivots > 0 || !r.finite) : (r.npivots == 0 && r.finite))
      same_bits = singular || (rd = cudss_pivots(M, T; deterministic = true); rd0 = cudss_pivots(M, T; pivot_epsilon = 0.0, deterministic = true); all(k -> rd0.x[k] === rd.x[k], eachindex(rd.x)))
      ok0 = r0.info == 0 && r0.npivots == 0 && (singular ? !r0.finite : (r0.finite && same_bits))
      ok = check!(ok1 && ok0, "cuDSS $label $T")
      @printf("  %-26s %-8s default: info %d npivots %d max|x| %.3e finite %s | epsilon 0: info %d npivots %d finite %s%s  %s\n", label, T, r.info, r.npivots, r.maxx, r.finite, r0.info, r0.npivots, r0.finite, singular ? "" : (same_bits ? ", bitwise as default" : ", DIFFERS from default"), ok ? "ok" : "FAIL")
    end
  end
  # repeatability: ten solves of one matrix, distinct bit patterns of x
  for T in (Float64, Float32), det in (false, true)
    xs = [cudss_pivots(A, T; deterministic = det).x for _ = 1:10]
    nd = length(unique(map(x -> reinterpret(UInt8, x), xs)))
    ok = !det || check!(nd == 1, "cuDSS deterministic_mode repeatable $T")
    @printf("  repeatability %-8s deterministic_mode %d: %2d distinct results of 10  %s\n", T, det, nd, det ? (ok ? "ok" : "FAIL") : "(reported)")
  end
end

# --- 2. solver runs against C1 -----------------------------------------------

function part2(cases)
  println("\n== 2. G1, G2 (stage 1) and G3, G4 (device loop) against C1")
  println("  dhist: largest difference of the mismatch history to C1 over its bound (1 passes; float32_mixed: HIST32 on the intermediate entries), entry 1 bitwise (b) or not (B!)")
  @printf("  %-16s %-4s %-5s %4s %4s %5s %5s %8s %10s %4s %4s %9s %6s %4s %6s %8s %s\n", "case", "var", "conv", "hist", "C1", "steps", "C1", "max|dV|", "dhist", "an", "ref", "qchanges", "clamps", "fb", "d2h/it", "B/it", "result")
  for name in cases
    p = casefile(name)
    c1 = run_c1(p)
    @printf("  %-16s %-4s %-5s %4d %4s %5d %5s %8s %10s %4s %4s %4d/%-4d %6d\n", name, "C1", c1.converged, c1.history, "", c1.steps, "", "", "", "", "", c1.qchanges, c1.qreenable, length(c1.active.clamps))
    for (var, prec, device) in (("G1", :float64, false), ("G2", :float32_mixed, false), ("G3", :float64, true), ("G4", :float32_mixed, true))
      g = device ? run_device(p, prec) : run_stage1(p, prec)
      dv = maximum(abs.(g.V .- c1.V))
      st = g.st
      fb = st.linear_solver_fallback_count
      f32 = prec === :float32_mixed
      # history length and Newton steps must equal C1; float32_mixed may take
      # other step counts (task section 9: reported, not forbidden), its
      # Q-limit switch iterations then shift, so only their number counts
      same_path = g.history == c1.history && g.steps == c1.steps
      same_q = g.qchanges == c1.qchanges && g.qreenable == c1.qreenable && (f32 ? length(g.qswitch_its) == length(c1.qswitch_its) : g.qswitch_its == c1.qswitch_its)
      same_active = g.active == c1.active
      # float64: the mismatch history by value and the Q-limit log entry by
      # entry; float32_mixed: the intermediate history entries by HIST32,
      # the log only reported
      hc = f32 ? hist_check32(g, c1) : hist_check(g, c1)
      same_log = g.qlog == c1.qlog
      ok = check!(g.converged == c1.converged && fb == 0 && dv <= DV_BOUND && same_q && same_active && hc.ok && (f32 || (same_path && same_log)) && (!device || get(st, :gpu_newton, false)), "$name $var")
      d2h_it = device ? get(st, :gpu_newton_d2h_count, 0) / g.history : NaN
      bytes_it = device ? get(st, :gpu_newton_d2h_bytes, 0) / g.history : NaN
      notes = String[]
      same_path || push!(notes, f32 ? "steps differ from C1 (float32_mixed, reported)" : "history or steps differ from C1")
      same_q || push!(notes, "Q-limit counters differ from C1")
      same_active || push!(notes, "final active set differs from C1")
      hc.ok || push!(notes, "mismatch history values differ from C1")
      same_log || push!(notes, f32 ? "Q-limit log differs from C1 (float32_mixed, reported)" : "Q-limit log differs from C1")
      @printf("  %-16s %-4s %-5s %4d %4d %5d %5d %8.1e %8.1e%-2s %4d %4d %4d/%-4d %6d %4d %6.2f %8.1f %s%s\n", name, var, g.converged, g.history, c1.history, g.steps, c1.steps, dv, hc.ratio, hc.first_bitwise ? "b" : "B!", st.linear_solver_analyze_count, st.linear_solver_refactor_count, g.qchanges, g.qreenable, length(g.active.clamps), fb, d2h_it, bytes_it, ok ? "ok" : "FAIL", isempty(notes) ? "" : "  (" * join(notes, "; ") * ")")
      fb > 0 && println("      fallback: ", first(st.compute_backend_reason, 200))
      if !same_active
        dtypes = findall(k -> g.active.types[k] != c1.active.types[k], eachindex(c1.active.types))
        println("      active set: node types differ at buses ", first(dtypes, 20), ", clamps only in ", var, " ", first(setdiff(g.active.clamps, c1.active.clamps), 20), ", only in C1 ", first(setdiff(c1.active.clamps, g.active.clamps), 20))
      end
      if dv > DV_BOUND
        k = argmax(abs.(g.V .- c1.V))
        @printf("      max |dV| at bus %d: |V| C1 %.6f, %s %.6f, angle difference %.3e deg; final mismatch C1 %.2e, %s %.2e; buses above the bound %d\n", k, abs(c1.V[k]), var, abs(g.V[k]), rad2deg(angle(g.V[k] / c1.V[k])), c1.st.mismatch_history[end], var, st.mismatch_history[end], count(>(DV_BOUND), abs.(g.V .- c1.V)))
      end
    end
    cooldown0(name, p)
  end
end

# Q limits without cooldown, rescue off: the handler reads abs(V) of a bus in
# the call that clamped it, so the sample must carry abs(V) of every bus that
# is PV in the data (a value it does not carry is an error). G3 against C1
# with the same settings, the Q-limit log entry by entry. The voltage margin
# stays at its default: with margin 0 the release of a bus that was PV in
# the measured iterate compares abs(V) with its set point at rounding level,
# and already stage 1 (V on the host) takes other releases than C1.
function cooldown0(name, p)
  qc = QLimitConfig(cooldown_iters = 0)
  c = run_c1(p; rescue = false, qlimits = qc)
  g = run_device(p, :float64; rescue = false, qlimits = qc)
  st = g.st
  dv = maximum(abs.(g.V .- c.V))
  fb = st.linear_solver_fallback_count
  hc = hist_check(g, c)
  same_log = g.qlog == c.qlog
  nrel(s) = count(e -> e[4] === :release, s.qlog)
  ok = check!(g.converged == c.converged && fb == 0 && g.history == c.history && g.steps == c.steps && hc.ok && same_log && g.qchanges == c.qchanges && g.qreenable == c.qreenable && g.qswitch_its == c.qswitch_its && g.active == c.active && (!c.converged || dv <= DV_BOUND) && get(st, :gpu_newton, false), "$name G3 cooldown 0")
  @printf("  %-16s %-4s cooldown 0: conv %s (C1 %s) hist %d (C1 %d) dhist %.1e%s q %d/%d (C1 %d/%d) log %d (C1 %d) releases %d (C1 %d) same log %s max|dV| %.1e fb %d  %s\n", name, "G3", g.converged, c.converged, g.history, c.history, hc.ratio, hc.first_bitwise ? "b" : "B!", g.qchanges, g.qreenable, c.qchanges, c.qreenable, length(g.qlog), length(c.qlog), nrel(g), nrel(c), same_log, dv, fb, ok ? "ok" : "FAIL")
  return ok
end

# --- 3. recovery paths by injected faults ------------------------------------

# A hook for EXT._FAULT_HOOK: acts once at (point, iteration).
mutable struct Fault
  point::Symbol
  it::Int
  action::Symbol      # :throw, :nan, :nanv, :flag (take the failure branch), :sticky, :sticky_async
  fired::Bool
end

function (f::Fault)(point::Symbol, it::Int, st)
  (f.fired || point !== f.point || it != f.it) && return false
  f.fired = true
  f.action === :throw && throw(CUDA.CuError(CUDA.CUDACore.CUDA_ERROR_LAUNCH_FAILED))
  f.action === :nan && (fill!(st.sol, NaN); return false)
  f.action === :nanv && (fill!(st.V, ComplexF64(NaN)); return false)
  f.action === :sticky && (sticky_fault!(); return false)
  f.action === :sticky_async && (sticky_fault!(; sync = false); return false)
  return true
end

function _k_bad_store!(a)
  # far outside any allocation: illegal address, the context is dead after it
  @inbounds a[1<<40] = 1.0
  return nothing
end

# sync = false: the fault surfaces at a later synchronizing call of the
# loop, as an asynchronous error from the device copies would
function sticky_fault!(; sync = true)
  a = CUDA.zeros(Float64, 1)
  @cuda threads = 1 _k_bad_store!(a)
  sync && CUDA.synchronize()
  return nothing
end

function with_fault(f, fault)
  EXT._FAULT_HOOK[] = fault
  try
    return with_logger(f, NullLogger())
  finally
    EXT._FAULT_HOOK[] = nothing
  end
end

# compare a recovered run with C1; `expect` is :c1 (same path, result and
# active set as C1) or (:singular, k) (the CPU redo of step k was singular);
# `values`: the mismatch history by value and the Q-limit log too (float64);
# false: the float32_mixed bound on the history, no log comparison
function judge(label, g, c1, fault, expect, reason_has; values = true)
  st = g.st
  dv = maximum(abs.(g.V .- c1.V))
  fb = st.linear_solver_fallback_count
  reason = st.compute_backend_reason
  # with expect (:singular, k) the first k entries of the C1 history
  cref = expect === :c1 ? c1 : (hist = c1.hist[1:min(expect[2], end)],)
  hc = values ? hist_check(g, cref) : hist_check32(g, cref)
  same_values = hc.ok && (!values || expect !== :c1 || g.qlog == c1.qlog)
  if expect === :c1
    ok = same_values && fault.fired && g.converged == c1.converged && g.history == c1.history && g.steps == c1.steps && g.qchanges == c1.qchanges && g.qreenable == c1.qreenable && g.qswitch_its == c1.qswitch_its && g.active == c1.active && fb == 1 && occursin(reason_has, reason) && (!c1.converged || dv <= DV_BOUND)
  else
    k = expect[2]
    # the CPU loop: a singular step in iteration k ends the solve with k
    # mismatch entries and k - 1 Newton steps
    ok = same_values && fault.fired && !g.converged && g.reason === :singular_newton_step && g.history == k && g.steps == k - 1 && fb == 1 && occursin(reason_has, reason)
  end
  check!(ok, label)
  @printf("  %-44s fired %-5s conv %-5s hist %3d (C1 %3d) steps %3d (C1 %3d) dhist %.1e%-2s q %d/%d (C1 %d/%d) max|dV| %.1e fb %d  %s\n", label, fault.fired, g.converged, g.history, c1.history, g.steps, c1.steps, hc.ratio, hc.first_bitwise ? "b" : "B!", g.qchanges, g.qreenable, c1.qchanges, c1.qreenable, dv, fb, ok ? "ok" : "FAIL")
  println("      reason: ", first(reason, 160), "  status reason ", g.reason)
  return ok
end

function part3(fault_case, vdep_case)
  println("\n== 3. recovery paths by injected faults ($fault_case, $vdep_case)")
  p = casefile(fault_case)
  c1 = run_c1(p; rescue = false)
  k = isempty(c1.qswitch_its) ? 2 : first(c1.qswitch_its)
  println("  C1: converged ", c1.converged, ", ", c1.history, " mismatch entries, ", c1.steps, " steps, first Q-limit switch in iteration ", k)
  cases = [
    # (label, fault, precision, expected, text in the reason)
    ("preflight error", Fault(:preflight, 0, :throw, false), :float64, :c1, "before the device Newton loop"),
    ("error at mismatch, iteration 1 (no step)", Fault(:mismatch, 1, :throw, false), :float64, :c1, "in iteration 1"),
    ("error at mismatch, step 2 unchecked", Fault(:mismatch, 3, :throw, false), :float64, :c1, "in iteration 2"),
    ("error while reading step 2", Fault(:check, 3, :throw, false), :float64, :c1, "in iteration 2"),
    ("failed check of step 2", Fault(:check, 3, :flag, false), :float64, :c1, "non-finite iterate"),
    ("non-finite step 2", Fault(:update, 2, :nan, false), :float64, :c1, "non-finite iterate"),
    ("non-finite step 2, float32_mixed", Fault(:update, 2, :nan, false), :float32_mixed, :c1, "non-finite iterate"),
    ("error after Q-limit switch, iteration $k", Fault(:switch, k, :throw, false), :float64, :c1, "in iteration $k"),
    # the norm after the switch of iteration k comes with the download of
    # k + 1; an error before that download computes it on the CPU
    ("error at mismatch after the switch of $k", Fault(:mismatch, k + 1, :throw, false), :float64, :c1, "in iteration $k"),
    # the device copies of the measured iterate (Vlast, Vbest) run before
    # the history entry is pushed: the CPU loop takes the iteration, the
    # entry must not stand twice
    ("error in the iterate copies, iteration 2", Fault(:records, 2, :throw, false), :float64, :c1, "in iteration 2"),
    ("error in the iterate copies, iteration $k", Fault(:records, k, :throw, false), :float64, :c1, "in iteration $k"),
    ("error in the step of iteration 2", Fault(:step, 2, :throw, false), :float64, :c1, "in iteration 2"),
    ("error in the step of iteration $k (switch)", Fault(:step, k, :throw, false), :float64, :c1, "in iteration $k"),
    # after the converged iteration only V and the records are read: no
    # step and no history entry may be added
    ("error reading V after convergence", Fault(:final, c1.history, :throw, false), :float64, :c1, "after the device loop ended"),
  ]
  for (label, fault, prec, expect, has) in cases
    g = with_fault(() -> run_device(p, prec; rescue = false), fault)
    judge(label, g, c1, fault, expect, has; values = prec === :float64)
  end
  # a singular CPU redo of step 2: the solve ends there as the CPU loop does
  let fault = Fault(:update, 2, :nan, false), cpu = Fault(:cpu_step, 2, :flag, false)
    hook = (pt, it, st) -> fault(pt, it, st) | cpu(pt, it, st)
    g = with_fault(() -> run_device(p, :float64; rescue = false), hook)
    judge("non-finite step 2, singular CPU redo", g, c1, cpu, (:singular, 2), "non-finite iterate")
  end
  # a non-finite norm in iteration 1 (no step pending) ends the loop; an
  # error reading V after it must give the run without that error
  let g0 = with_fault(() -> run_device(p, :float64; rescue = false), Fault(:mismatch, 1, :nanv, false))
    nan1 = Fault(:mismatch, 1, :nanv, false)
    fin = Fault(:final, 1, :throw, false)
    g = with_fault(() -> run_device(p, :float64; rescue = false), (pt, it, st) -> nan1(pt, it, st) | fin(pt, it, st))
    fb = g.st.linear_solver_fallback_count
    ok = check!(fin.fired && g0.reason === :nr_nonfinite && g.reason === g0.reason && g.history == g0.history == 1 && isequal(g.hist, g0.hist) && g.steps == g0.steps && g.converged == g0.converged && fb == 1 && occursin("after the device loop ended", g.st.compute_backend_reason), "non-finite norm 1, error reading V")
    @printf("  %-44s fired %-5s conv %-5s hist %3d (without the error %3d) steps %3d (%3d) status reason %s (%s) fb %d  %s\n", "non-finite norm 1, error reading V", fin.fired, g.converged, g.history, g0.history, g.steps, g0.steps, g.reason, g0.reason, fb, ok ? "ok" : "FAIL")
    println("      reason: ", first(g.st.compute_backend_reason, 160))
  end
  # the last step of a run that uses up max_iter is checked after the loop
  let c13 = run_c1(p; max_iter = 3, rescue = false), fault = Fault(:check, 4, :flag, false)
    g = with_fault(() -> run_device(p, :float64; max_iter = 3, rescue = false), fault)
    judge("max_iter 3, failed check of step 3", g, c13, fault, :c1, "in iteration 3")
    # step 3 passed its check and is confirmed; a GPU error after that, in
    # the host records or in the outcome, must neither redo the step nor
    # add a history entry
    for (label, pt) in (("max_iter 3, error in the records after step 3", :host_records), ("max_iter 3, error in the outcome after step 3", :outcome))
      f = Fault(pt, 3, :throw, false)
      g = with_fault(() -> run_device(p, :float64; max_iter = 3, rescue = false), f)
      judge(label, g, c13, f, :c1, "after the device loop ended")
    end
  end
  # fallback off: the failure goes up
  let fault = Fault(:update, 2, :nan, false)
    threw = try
      with_fault(() -> run_device(p, :float64; fallback = false, rescue = false), fault)
      false
    catch e
      e isa ErrorException && occursin("non-finite iterate", e.msg)
    end
    check!(threw, "fallback off throws")
    @printf("  %-44s throws %s  %s\n", "non-finite step 2, fallback off", threw, threw ? "ok" : "FAIL")
  end
  # voltage-dependent injections: V goes to the host each iteration and the
  # step is checked on it before the rebuild
  pv = casefile(vdep_case)
  cv = run_c1(pv; rescue = false)
  println("  C1 $vdep_case: converged ", cv.converged, ", ", cv.history, " mismatch entries, ", cv.steps, " steps")
  for (label, fault, has) in (("vdep: error in the uploads, iteration 3", Fault(:vdep, 3, :throw, false), "in iteration 3"), ("vdep: failed check of step 2", Fault(:check, 3, :flag, false), "non-finite iterate"), ("vdep: error at mismatch, iteration 3", Fault(:mismatch, 3, :throw, false), "in iteration 3"))
    g = with_fault(() -> run_device(pv, :float64; rescue = false), fault)
    judge(label, g, cv, fault, :c1, has)
  end
end

# --- 4. a kernel fault that kills the context --------------------------------

function sticky_child(spec)
  parts = split(spec, ":")
  point, it = parts[1], parts[2]
  action = length(parts) > 2 && parts[3] == "async" ? :sticky_async : :sticky
  p = casefile(get(ENV, "FAULT_CASE", "case2869pegase"))
  c1 = run_c1(p; rescue = false)
  fault = Fault(Symbol(point), parse(Int, it), action, false)
  g = with_fault(() -> run_device(p, :float64; rescue = false), fault)
  # in both cases the CPU runs the whole solve: same path as C1
  ok = judge("kernel fault at $spec", g, c1, fault, :c1, point == "preflight_kernels" ? "before the device Newton loop" : "context is lost")
  dead = try
    CUDA.zeros(Float64, 1)
    CUDA.synchronize()
    false
  catch e
    e isa CUDA.CuError
  end
  println("  context dead after the fault: ", dead)
  println("STICKY_RESULT ", ok && dead)
  flush(stdout)
  exit(ok && dead ? 0 : 1)
end

function part4(fault_case)
  println("\n== 4. kernel fault (illegal address) that kills the context, one child process each ($fault_case)")
  c1 = run_c1(casefile(fault_case); rescue = false)
  k = isempty(c1.qswitch_its) ? 2 : first(c1.qswitch_its)
  # records:2:async: the fault is launched in the device copies of the
  # measured iterate and not waited for; it surfaces at a later call, and
  # the solve must not be redone from the possibly bad Vlast
  for spec in ("preflight_kernels:0", "mismatch:3", "step:2", "switch:$k", "records:2:async")
    cmd = setenv(`$(Base.julia_cmd()) --project=$(Base.active_project()) $(@__FILE__)`, merge(ENV, Dict("STICKY" => spec, "FAULT_CASE" => fault_case)))
    out = read(ignorestatus(cmd), String)
    lines = filter(l -> occursin("kernel fault", l) || occursin("reason:", l) || occursin("context dead", l) || occursin("STICKY_RESULT", l), split(out, '\n'))
    foreach(l -> println(l), lines)
    check!(occursin("STICKY_RESULT true", out), "sticky $spec")
  end
end

# --- main ----------------------------------------------------------------------

function main()
  haskey(ENV, "STICKY") && return sticky_child(ENV["STICKY"])
  cases = isempty(ARGS) ? STUDY_CASES : ARGS
  started = now()
  commit = try
    readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`)
  catch e
    e isa ProcessFailedException || e isa Base.IOError || rethrow(e)
    "unknown"
  end
  parts = parse.(Int, split(get(ENV, "PARTS", "1,2,3,4"), ","))
  println("check_device_newton_gpu  started ", started, "  commit ", commit, "  Julia ", VERSION, "  CUDA.jl ", pkgversion(CUDA), "  CUDSS.jl ", pkgversion(CUDSS), "  device ", CUDA.name(CUDA.device()), "  parts ", join(parts, ","))
  1 in parts && part1()
  2 in parts && part2(cases)
  3 in parts && part3(get(ENV, "FAULT_CASE", "case2869pegase"), get(ENV, "VDEP_CASE", "case6495rte"))
  4 in parts && part4(get(ENV, "FAULT_CASE", "case2869pegase"))
  println("\nfinished ", now(), "  failed checks: ", length(FAILS))
  foreach(f -> println("  FAIL ", f), FAILS)
  exit(isempty(FAILS) ? 0 : 1)
end

Base.invokelatest(main)

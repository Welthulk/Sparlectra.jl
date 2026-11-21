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
# file: benchmark/gpu/repeat_fresh_gpu.jl
# purpose: GPU study, stage 1b: repeatability of one GPU variant across
#          fresh processes. C1 once, then RUNS child processes per variant,
#          each one solve; compared with C1: convergence, history length,
#          final active set (node types and Q-limit clamps), max |dV|, and a
#          SHA-1 of the bits of V. The cuDSS mode is the configuration
#          default (deterministic) unless FAST=1.
#
# Run: julia --project=benchmark/gpu benchmark/gpu/repeat_fresh_gpu.jl [case [variants]]
#   default case_ACTIVSg70k G2,G4; RUNS (default 10); FAST=1 for
#   gpu_fast_mode. Exit code 0 only if every run of a variant gives the same
#   V bits and the active set of C1.
#   Internal: CHILD=<variant> OUT=<file> runs one solve.

using Sparlectra, CUDA, CUDSS, Printf, Dates, Logging, Serialization, SHA

const QUIET = OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)
const FAST = get(ENV, "FAST", "0") == "1"

cfg(cb) = SparlectraConfig(powerflow = PowerFlowConfig(compute_backend = cb), output = QUIET)

function solve(p, var::String)
  var == "C1" && return run_sparlectra(casefile = p, config = cfg(Sparlectra.ComputeBackendConfig()))
  prec = var in ("G2", "G4") ? :float32_mixed : :float64
  cb = Sparlectra.ComputeBackendConfig(mode = :gpu, gpu_precision = prec, gpu_fast_mode = FAST)
  var in ("G3", "G4") && return Sparlectra._with_device_newton(() -> run_sparlectra(casefile = p, config = cfg(cb)))
  return run_sparlectra(casefile = p, config = cfg(cb))
end

function summary_of(r)
  st = Sparlectra.rectangular_pf_status(r.net)
  V = Sparlectra.buildVoltageVector(r.net)
  return (
    converged = r.final_converged,
    hist = length(st.mismatch_history),
    types = [Sparlectra.getNodeType(nd) for nd in r.net.nodeVec],
    clamps = sort!(collect(r.net.qLimitEvents)),
    V = V,
    hash = bytes2hex(sha1(reinterpret(UInt8, V)))[1:12],
    mode = get(st, :gpu_cudss_mode, :none),
    gpu_newton = get(st, :gpu_newton, false),
    fallback = get(st, :linear_solver_fallback_count, 0),
  )
end

function child(var, out)
  p = Sparlectra.ensure_casefile(ARGS_CASE() * ".m"; to_jl = false)
  t0 = time()
  s = with_logger(() -> summary_of(solve(p, var)), NullLogger())
  serialize(out, merge(s, (elapsed = time() - t0,)))
  return nothing
end

ARGS_CASE() = isempty(ARGS) ? "case_ACTIVSg70k" : ARGS[1]

function main()
  haskey(ENV, "CHILD") && return child(ENV["CHILD"], ENV["OUT"])
  name = ARGS_CASE()
  variants = length(ARGS) >= 2 ? split(ARGS[2], ",") : ["G2", "G4"]
  runs = parse(Int, get(ENV, "RUNS", "10"))
  commit = readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`)
  println("repeat_fresh_gpu  started ", now(), "  commit ", commit, "  case ", name, "  variants ", join(variants, ","), "  runs ", runs, "  cuDSS mode ", FAST ? "fast" : "deterministic (default)")
  tmp = mktempdir()
  spawn(var, k) = begin
    out = joinpath(tmp, "$(var)_$(k).jls")
    cmd = setenv(`$(Base.julia_cmd()) --project=$(Base.active_project()) $(@__FILE__) $name`, merge(ENV, Dict("CHILD" => var, "OUT" => out)))
    run(pipeline(cmd; stdout = devnull, stderr = devnull))
    deserialize(out)
  end
  c1 = spawn("C1", 0)
  @printf("C1 converged %s hist %d clamps %d V %s (%.1f s)\n", c1.converged, c1.hist, length(c1.clamps), c1.hash, c1.elapsed)
  ok = true
  for var in variants
    hashes = String[]
    for k = 1:runs
      g = spawn(var, k)
      same_active = g.types == c1.types && g.clamps == c1.clamps
      dv = maximum(abs.(g.V .- c1.V))
      push!(hashes, g.hash)
      @printf("%2d %s conv %s hist %d clamps %d active %s max|dV| %.2e V %s mode %s device loop %s fallback %d (%.1f s)\n", k, var, g.converged, g.hist, length(g.clamps), same_active ? "same as C1" : "DIFFERS", dv, g.hash, g.mode, g.gpu_newton, g.fallback, g.elapsed)
      if !same_active
        d = findall(i -> g.types[i] != c1.types[i], eachindex(g.types))
        println("   node types differ at buses ", first(d, 20), "; clamps only here ", first(setdiff(g.clamps, c1.clamps), 20), ", only in C1 ", first(setdiff(c1.clamps, g.clamps), 20))
      end
      ok &= same_active && g.converged
      flush(stdout)
    end
    nd = length(unique(hashes))
    ok &= nd == 1
    println(var, ": ", nd, " distinct V of ", runs, " fresh processes")
  end
  println("finished ", now(), "  ", ok ? "all runs equal and with the active set of C1" : "NOT all runs equal or with the active set of C1")
  exit(ok ? 0 : 1)
end

Base.invokelatest(main)

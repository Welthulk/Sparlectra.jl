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
# file: benchmark/gpu/probe_cudss_sync_gpu.jl
# purpose: GPU study, stage 1b: does a cuDSS refactorization or solve block
#          the host until the device work queued before it is done? The
#          CUPTI trace of the device Newton loop shows a cudaMalloc and a
#          4 B host-to-device cudaMemcpyAsync inside every refactorization;
#          either can synchronize. Probe: a kernel that sleeps on the device
#          for SPIN_MS is queued on the stream, then the call under test is
#          made without a synchronize, and the host clock reads when the
#          call returns. A call that returns after about SPIN_MS waited for
#          the device; one that returns at once did not. Controls: an 8 B
#          device-to-host copy (blocks) and a device-to-device copy (does
#          not). Device Jacobian at the C1 solution as in
#          timing_refactor_solve_gpu.jl, pivot_epsilon 0, both modes,
#          Float64 and Float32.
#
# Run: julia --project=benchmark/gpu benchmark/gpu/probe_cudss_sync_gpu.jl [case ...]
#   default cases case118 case_ACTIVSg25k; SPIN_MS (default 50), REPEATS
#   (default 5, median of the host return times).

using Sparlectra, CUDA, CUDSS, SparseArrays, LinearAlgebra, Statistics, Printf, Dates, Logging

const EXT = Base.get_extension(Sparlectra, :SparlectraCUDAExt)
EXT === nothing && error("SparlectraCUDAExt not loaded")

const QUIET = OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)
const SPIN_MS = parse(Int, get(ENV, "SPIN_MS", "50"))
const REPEATS = parse(Int, get(ENV, "REPEATS", "5"))

bus_type_symbol(t) = t == Sparlectra.Slack ? :Slack : t == Sparlectra.PV ? :PV : :PQ

# one thread sleeps on the device in steps of 1 ms
function _k_spin!(ms)
  for _ = 1:ms
    CUDA.nanosleep(UInt32(1_000_000))
  end
  return nothing
end

spin!() = @cuda threads = 1 _k_spin!(Int32(SPIN_MS))

# host milliseconds until f() returns, with the spin kernel queued before it;
# then the device milliseconds until everything is done
function probe(f)
  ret = Float64[]
  tot = Float64[]
  for _ = 1:(REPEATS+1)
    CUDA.synchronize()
    t0 = time_ns()
    spin!()
    f()
    t1 = time_ns()
    CUDA.synchronize()
    t2 = time_ns()
    push!(ret, (t1 - t0) / 1.0e6)
    push!(tot, (t2 - t0) / 1.0e6)
  end
  # the first is warm-up
  return (ret = median(ret[2:end]), tot = median(tot[2:end]))
end

function setup(net, T, fast::Bool)
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
  EXT._cudss_configure!(solver, fast)
  cudss_set(solver, "pivot_epsilon", zero(T))
  X = CudssMatrix(CUDA.zeros(T, nvar))
  B = CudssMatrix(CuVector{T}(T.(-g.F[:, 1])))
  cudss("analysis", solver, X, B)
  cudss("factorization", solver, X, B)
  CUDA.synchronize()
  return (solver = solver, X = X, B = B, nvar = nvar)
end

function main()
  cases = isempty(ARGS) ? ["case118", "case_ACTIVSg25k"] : ARGS
  commit = try
    readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`)
  catch e
    e isa ProcessFailedException || e isa Base.IOError || rethrow(e)
    "unknown"
  end
  println("probe_cudss_sync_gpu  started ", now(), "  commit ", commit, "  Julia ", VERSION, "  threads ", Threads.nthreads(), "  CUDA.jl ", pkgversion(CUDA), "  CUDSS.jl ", pkgversion(CUDSS), "  device ", CUDA.name(CUDA.device()), "  spin ", SPIN_MS, " ms, repeats ", REPEATS, ", medians")
  a = CUDA.zeros(Float64, 1)
  b = CUDA.zeros(Float64, 1)
  h = zeros(Float64, 1)
  sp = probe(() -> nothing)
  dh = probe(() -> copyto!(h, a))
  dd = probe(() -> copyto!(b, a))
  @printf("  spin alone: returns after %.2f ms, done after %.2f ms\n", sp.ret, sp.tot)
  @printf("  control 8 B device-to-host copy: returns after %.2f ms (blocks)\n", dh.ret)
  @printf("  control device-to-device copy: returns after %.2f ms (does not block)\n", dd.ret)
  for name in cases
    p = Sparlectra.ensure_casefile(name * ".m"; to_jl = false)
    r = with_logger(() -> run_sparlectra(casefile = p, config = SparlectraConfig(output = QUIET)), NullLogger())
    println("\n", name, ": n ", length(r.net.nodeVec))
    for T in (Float64, Float32), fast in (false, true)
      s = setup(r.net, T, fast)
      rf = probe(() -> cudss("refactorization", s.solver, s.X, s.B))
      so = probe(() -> cudss("solve", s.solver, s.X, s.B))
      both = probe(() -> (cudss("refactorization", s.solver, s.X, s.B); cudss("solve", s.solver, s.X, s.B)))
      verdict(x) = x.ret >= 0.5 * SPIN_MS ? "waits" : "does not wait"
      @printf("  %-7s %-13s nvar %7d | refactorization returns after %7.2f ms (%s) | solve %7.2f ms (%s) | both %7.2f ms (%s), done after %7.2f ms\n",
        T, fast ? "fast" : "deterministic", s.nvar, rf.ret, verdict(rf), so.ret, verdict(so), both.ret, verdict(both), both.tot)
      flush(stdout)
    end
  end
  println("finished ", now())
end

Base.invokelatest(main)

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

# file: examples/matpower_import_case14.jl
using Sparlectra
import Sparlectra: MatpowerIO
using BenchmarkTools
using Printf

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

case = "case30.m"      # or "case14.m"
flatstart = true
verbose = 1  # 0=quiet, 1=normal, 2=debug
cooldown_iters = 3  # number of iterations to continue after convergence (for logging)
q_hyst_pu = 0.01  # hysteresis for Q-limit enforcement (pu of max Q)

println("Sparlectra version: ", Sparlectra.version(), "\n")
println("Importing MATPOWER case file: $case\n")
println("DEMO: Matpower_import_case14.jl\n  ------------------------------\n")


# -----------------------------------------------------------------------------
# Ensure case is available locally (download on demand)
# -----------------------------------------------------------------------------

# This will:
# - download case14.m into data/mpower/ if missing
# - generate case14.jl if requested and missing
local_case = Sparlectra.FetchMatpowerCase.ensure_casefile(
    case;
    outdir = Sparlectra.MPOWER_DIR,
    to_jl = true,
    overwrite = false,
)

# -----------------------------------------------------------------------------
# Read case (Julia or MATPOWER)
# -----------------------------------------------------------------------------

mpc = MatpowerIO.read_case(local_case)

net = createNetFromMatPowerFile(filename = local_case, log = (verbose > 1), flatstart = flatstart, cooldown = cooldown_iters, q_hyst_pu = q_hyst_pu)

# -----------------------------------------------------------------------------
# Benchmark helper: benchmark exactly run_acpflow(...)
# -----------------------------------------------------------------------------
function bench_run_acpflow(; casefile::String, methods::Vector{Symbol},
    opt_fd::Bool=true, opt_sparse::Bool=true, opt_flatstart::Bool=true,
    verbose::Int=0, cooldown_iters::Int=0, q_hyst_pu::Float64=0.0,
    seconds::Float64=2.0, samples::Int=50,
    show_once::Bool=false)

  results = Dict{Symbol,Any}()

  # Warmup (compile) once per method with minimal output
  for m in methods
    run_acpflow(
      casefile = casefile,
      opt_fd = opt_fd,
      opt_sparse = opt_sparse,
      method = m,
      opt_flatstart = opt_flatstart,
      show_results = false,
      verbose = 0,
      cooldown_iters = cooldown_iters,
      q_hyst_pu = q_hyst_pu,
    )
  end

  println("\n==================== Benchmark run_acpflow ====================")
  println("casefile        = ", casefile)
  println("opt_fd          = ", opt_fd, "   opt_sparse = ", opt_sparse, "   flatstart = ", opt_flatstart)
  println("cooldown_iters  = ", cooldown_iters, "   q_hyst_pu = ", q_hyst_pu)
  println("seconds/method  = ", seconds, "   samples = ", samples)
  println("===============================================================\n")

    for m in methods
    benchable = @benchmarkable run_acpflow(
        casefile = casefile_,
        opt_fd = opt_fd_,
        opt_sparse = opt_sparse_,
        method = method_,
        opt_flatstart = opt_flatstart_,
        show_results = false,
        verbose = 0,
        cooldown_iters = cooldown_iters_,
        q_hyst_pu = q_hyst_pu_,
      ) setup=(
        casefile_ = $casefile;
        opt_fd_ = $opt_fd;
        opt_sparse_ = $opt_sparse;
        opt_flatstart_ = $opt_flatstart;
        cooldown_iters_ = $cooldown_iters;
        q_hyst_pu_ = $q_hyst_pu;
        method_ = $m;
      )

    b = run(benchable; seconds=seconds, samples=samples)
    results[m] = b

    tmed = BenchmarkTools.median(b).time / 1e6
    tmin = BenchmarkTools.minimum(b).time / 1e6
    tmax = BenchmarkTools.maximum(b).time / 1e6
    alloc = BenchmarkTools.median(b).memory
    allocs = BenchmarkTools.median(b).allocs

    println("method = ", m)
    println("  time:  median = ", round(tmed, digits=4), " ms",
            "   min = ", round(tmin, digits=4), " ms",
            "   max = ", round(tmax, digits=4), " ms")
    println("  alloc: ", alloc, " bytes   allocs: ", allocs, "\n")
  end

  # Optional: show results once (not benchmarked)
  if show_once
    println("\n--- One non-benchmarked run with show_results=true (for verification) ---\n")
    for m in methods
      println("Running power flow with method: $m\n")
      run_acpflow(
        casefile = casefile,
        opt_fd = opt_fd,
        opt_sparse = opt_sparse,
        method = m,
        opt_flatstart = opt_flatstart,
        show_results = true,
        verbose = verbose,
        cooldown_iters = cooldown_iters,
        q_hyst_pu = q_hyst_pu,
      )
    end
  end

  return results
end

methods = [:polar_full, :rectangular, :classic]

bench = bench_run_acpflow(
  casefile = basename(local_case),
  methods = methods,
  opt_fd = true,
  opt_sparse = true,
  opt_flatstart = flatstart,
  verbose = verbose,
  cooldown_iters = cooldown_iters,
  q_hyst_pu = q_hyst_pu,
  seconds = 2.0,
  samples = 10,
  show_once = true, 
)

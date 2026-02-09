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
using Dates
# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

#case = "case30.m"      # or "case14.m"
case = "case118.m"      # or "case14.m"
#case = "case39.m"      # or "case30.m"
#case = "case57.m"     # or "case39.m"
#case = "case300.m"     # or "case39.m"
#case = "case18.m"     # or "case14.m"
#case = "case141.m"     # or "case14.m"
#case = "case69.m"     # or "case14.m"
#case = "case85.m"     # or "case14.m"

flatstart = false
verbose = 1  # 0=quiet, 1=normal, 2=debug
cooldown_iters = 10  # number of iterations to continue after convergence (for logging)
q_hyst_pu = 0.1  # hysteresis for Q-limit enforcement (pu of max Q)
print("\e[2J\e[H") # clear screen and move cursor to home position
println("----------------------------------------------------------------------------------------------")
println("Sparlectra version: ", Sparlectra.version(), "\n")
println("Importing MATPOWER case file: $case\n")

# -----------------------------------------------------------------------------
# Output redirection (write verbose output to git-ignored file)
# -----------------------------------------------------------------------------
const OUTDIR = joinpath(@__DIR__, "_out")
mkpath(OUTDIR)

timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
logfile   = joinpath(OUTDIR, "run_$(case)_$(timestamp).log")

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

# -----------------------------------------------------------------------------
# Compare against MATPOWER reference (if present)
# -----------------------------------------------------------------------------
show_diff = true       # <- dein Flag
tol_vm    = 2e-2
tol_va    = 5e-1

function mp_has_vm_va(mpc)
    size(mpc.bus, 2) >= 9 || return false
    vm = mpc.bus[:, 8]
    va = mpc.bus[:, 9]

    # need finite everywhere (not just "any"), otherwise it's not a trustworthy reference
    all(isfinite, vm) || return false
    all(isfinite, va) || return false

    # reject near-flat placeholders
    vm_dev = maximum(abs.(vm .- 1.0))
    va_dev = maximum(abs.(va .- 0.0))
    return (vm_dev > 1e-3) || (va_dev > 1e-2)
end
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
     @printf("bench  method=%-12s  med=%8.4f ms  min=%8.4f ms  alloc=%9d B  allocs=%d\n", String(m), tmed, tmin, alloc, allocs)

  end

  # Optional: show results once (not benchmarked)
  # Optional: show results once (not benchmarked)
  if show_once
    # write full solver output to logfile; keep terminal clean
    local summaries = Vector{NamedTuple{(:method,:converged,:max_dvm,:max_dva,:cmp_ok),Tuple{Symbol,Bool,Float64,Float64,Bool}}}()

    open(logfile, "w") do io
      redirect_stdout(io) do
        redirect_stderr(io) do
          println("Sparlectra version: ", Sparlectra.version())
          println("casefile: ", casefile)
          println("timestamp: ", Dates.now())
          println()

          for m in methods
            println("=================================================================")
            println("RUN method = ", m)
            println("=================================================================\n")

            net_res = run_acpflow(
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

            # compare (still computed, but details go to logfile)
            if mp_has_vm_va(mpc)
              ok, stats = MatpowerIO.compare_vm_va(
                net_res, mpc;
                show_diff = show_diff,
                tol_vm    = tol_vm,
                tol_va    = tol_va,
                maxlines  = 20
              )
              push!(summaries, (method=m, converged=true,
                                max_dvm=Float64(get(stats, :max_dvm, NaN)),
                                max_dva=Float64(get(stats, :max_dva, NaN)),
                                cmp_ok=ok))
            else
              push!(summaries, (method=m, converged=true, max_dvm=NaN, max_dva=NaN, cmp_ok=false))
              println("Compare skipped: no solution-like VM/VA in mpc.bus(:,8:9)")
            end

            println()
          end
        end
      end
    end

    # print short summary to terminal
    println("\n==================== Summary ====================")
    println("logfile: ", logfile)
    for s in summaries
      if isnan(s.max_dvm) || isnan(s.max_dva)
        @printf("method=%-12s  compare=SKIP\n", String(s.method))
      else
        @printf("method=%-12s  compare=%s  max|dVm|=%8.5f pu  max|dVa|=%7.4f deg\n",
                String(s.method),
                s.cmp_ok ? "OK " : "FAIL",
                s.max_dvm, s.max_dva)
      end
    end
    println("=================================================\n")
  end
  return results
end

methods = [:polar_full, :rectangular, :classic]


function main()
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
        samples = 50,
        show_once = true,
    )
    return bench
end

main()
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

# file: examples/matpower_import.jl
using Sparlectra
import Sparlectra: MatpowerIO
using BenchmarkTools
using Printf
using Dates

# -----------------------------------------------------------------------------
# YAML config helpers (simple subset)
# -----------------------------------------------------------------------------
function _parse_yaml_scalar(raw::AbstractString)
  s = strip(raw)
  isempty(s) && return nothing

  if (startswith(s, "\"") && endswith(s, "\"")) || (startswith(s, "'") && endswith(s, "'"))
    return s[2:end-1]
  end
  ls = lowercase(s)
  ls == "true" && return true
  ls == "false" && return false
  ls == "null" && return nothing

  iv = tryparse(Int, s)
  !isnothing(iv) && return iv
  fv = tryparse(Float64, s)
  !isnothing(fv) && return fv
  return s
end

function _parse_yaml_list(raw::AbstractString)
  inner = strip(raw)[2:end-1]
  isempty(strip(inner)) && return Any[]
  return [_parse_yaml_scalar(part) for part in split(inner, ",")]
end

function load_yaml_config(path::AbstractString)
  isempty(path) && return Dict{String,Any}()
  isfile(path) || error("YAML config file not found: $path")

  cfg = Dict{String,Any}()
  for line in eachline(path)
    stripped = strip(line)
    isempty(stripped) && continue
    startswith(stripped, "#") && continue
    occursin(":", stripped) || continue

    key, value_raw = split(stripped, ":"; limit = 2)
    key = strip(key)
    value_raw = strip(split(value_raw, "#"; limit = 2)[1]) # remove inline comments

    if startswith(value_raw, "[") && endswith(value_raw, "]")
      cfg[key] = _parse_yaml_list(value_raw)
    else
      cfg[key] = _parse_yaml_scalar(value_raw)
    end
  end
  return cfg
end

function _as_symbol_vec(v)
  v isa AbstractVector || return Symbol[]
  return Symbol.(String.(v))
end

function _as_int_vec(v)
  v isa AbstractVector || return Int[]
  return Int[x for x in v]
end

function _yaml_path_from_inputs()
  !isempty(ARGS) && return ARGS[1]
  env_path = get(ENV, "SPARLECTRA_MATPOWER_IMPORT_YAML", "")
  !isempty(env_path) && return env_path

  local_default = joinpath(@__DIR__, "matpower_import.yaml")
  isfile(local_default) && return local_default
  return ""
end
# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
const DEFAULT_CASE = "case141.m"
const DEFAULT_METHODS = [:polar_full, :rectangular, :classic]
const METHODS = [:polar_full, :rectangular, :classic]

# -----------------------------------------------------------------------------
# Output redirection (write verbose output to git-ignored file)
# -----------------------------------------------------------------------------
const OUTDIR = joinpath(@__DIR__, "_out")
mkpath(OUTDIR)

# -----------------------------------------------------------------------------
# Compare against MATPOWER reference (if present)
# -----------------------------------------------------------------------------
show_diff = true
tol_vm = 2e-2
tol_va = 5e-1

# Keep scenario-specific benchmark options centralized here.
# If additional case-specific options are introduced (e.g. PV->PQ lock lists),
# they can be added in this table instead of using local `if case == ...` blocks.
const CASE_BENCH_OVERRIDES = Dict{String,NamedTuple}(
  # Example:
  # "case1951rte.m" => (; lock_pv_to_pq_buses = [44]),
)

function bench_config_for_case(case_name::AbstractString, yaml_cfg::Dict{String,Any})
  base = (;
    opt_fd = false,
    opt_sparse = true,
    opt_flatstart = false,
    verbose = 1,
    cooldown_iters = 0,
    q_hyst_pu = 0.0,
    lock_pv_to_pq_buses = Int[],
    seconds = 2.0,
    samples = 50,
    show_once = true,
  )
  case_override = get(CASE_BENCH_OVERRIDES, String(case_name), (;))
  yaml_override = (;)
  if !isempty(yaml_cfg)
    yaml_override = (;
      opt_fd = Bool(get(yaml_cfg, "opt_fd", base.opt_fd)),
      opt_sparse = Bool(get(yaml_cfg, "opt_sparse", base.opt_sparse)),
      opt_flatstart = Bool(get(yaml_cfg, "opt_flatstart", base.opt_flatstart)),
      verbose = Int(get(yaml_cfg, "verbose", base.verbose)),
      cooldown_iters = Int(get(yaml_cfg, "cooldown_iters", base.cooldown_iters)),
      q_hyst_pu = Float64(get(yaml_cfg, "q_hyst_pu", base.q_hyst_pu)),
      lock_pv_to_pq_buses = _as_int_vec(get(yaml_cfg, "lock_pv_to_pq_buses", base.lock_pv_to_pq_buses)),
      seconds = Float64(get(yaml_cfg, "seconds", base.seconds)),
      samples = Int(get(yaml_cfg, "samples", base.samples)),
      show_once = Bool(get(yaml_cfg, "show_once", base.show_once)),
    )
  end
  return merge(base, case_override, yaml_override)
end

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
function bench_run_acpflow(; casefile::String, methods::Vector{Symbol}, mpc, logfile::String, show_diff::Bool, tol_vm::Float64, tol_va::Float64, opt_fd::Bool = true, opt_sparse::Bool = true, opt_flatstart::Bool = true, verbose::Int = 0, cooldown_iters::Int = 0, q_hyst_pu::Float64 = 0.0, lock_pv_to_pq_buses::AbstractVector{Int} = Int[], seconds::Float64 = 2.0, samples::Int = 50, show_once::Bool = false)
  results = Dict{Symbol,Any}()

  # Warmup (compile) once per method with minimal output
  for m in methods
    run_acpflow(casefile = casefile, opt_fd = opt_fd, opt_sparse = opt_sparse, method = m, opt_flatstart = opt_flatstart, show_results = false, verbose = 0, cooldown_iters = cooldown_iters, q_hyst_pu = q_hyst_pu, lock_pv_to_pq_buses = lock_pv_to_pq_buses)
  end

  println("\n==================== Benchmark run_acpflow ====================")
  println("casefile        = ", casefile)
  println("opt_fd          = ", opt_fd, "   opt_sparse = ", opt_sparse, "   flatstart = ", opt_flatstart)
  println("cooldown_iters  = ", cooldown_iters, "   q_hyst_pu = ", q_hyst_pu)
  println("lock PV->PQ     = ", collect(lock_pv_to_pq_buses))
  println("seconds/method  = ", seconds, "   samples = ", samples)
  println("===============================================================\n")

  for m in methods
    benchable = @benchmarkable run_acpflow(casefile = casefile_, opt_fd = opt_fd_, opt_sparse = opt_sparse_, method = method_, opt_flatstart = opt_flatstart_, show_results = false, verbose = 0, cooldown_iters = cooldown_iters_, q_hyst_pu = q_hyst_pu_, lock_pv_to_pq_buses = lock_pv_to_pq_buses_) setup = (casefile_ = $casefile;
    opt_fd_ = $opt_fd;
    opt_sparse_ = $opt_sparse;
    opt_flatstart_ = $opt_flatstart;
    cooldown_iters_ = $cooldown_iters;
    q_hyst_pu_ = $q_hyst_pu;
    pv_table_rows_ = $pv_table_rows;
    check_q_limit_signs_ = $check_q_limit_signs;
    autocorrect_q_limit_signs_ = $autocorrect_q_limit_signs;
    validate_limits_after_pf_ = $validate_limits_after_pf;
    q_limit_violation_headroom_ = $q_limit_violation_headroom;
    lock_pv_to_pq_buses_ = $lock_pv_to_pq_buses;
    method_ = $m)

    b = run(benchable; seconds = seconds, samples = samples)
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
    local summaries = Vector{NamedTuple{(:method, :converged, :max_dvm, :max_dva, :cmp_ok),Tuple{Symbol,Bool,Float64,Float64,Bool}}}()

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

            net_res = run_acpflow(casefile = casefile, opt_fd = opt_fd, opt_sparse = opt_sparse, method = m, opt_flatstart = opt_flatstart, show_results = true, verbose = verbose, cooldown_iters = cooldown_iters, q_hyst_pu = q_hyst_pu, lock_pv_to_pq_buses = lock_pv_to_pq_buses)

            # compare (still computed, but details go to logfile)
            if mp_has_vm_va(mpc)
              ok, stats = MatpowerIO.compare_vm_va(net_res, mpc; show_diff = show_diff, tol_vm = tol_vm, tol_va = tol_va, maxlines = 20)
              push!(summaries, (method = m, converged = true, max_dvm = Float64(get(stats, :max_dvm, NaN)), max_dva = Float64(get(stats, :max_dva, NaN)), cmp_ok = ok))
            else
              push!(summaries, (method = m, converged = true, max_dvm = NaN, max_dva = NaN, cmp_ok = false))
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
        @printf("method=%-12s  compare=%s  max|dVm|=%8.5f pu  max|dVa|=%7.4f deg\n", String(s.method), s.cmp_ok ? "OK " : "FAIL", s.max_dvm, s.max_dva)
      end
    end
    println("=================================================\n")
  end
  return results
end

function main()
  yaml_path = _yaml_path_from_inputs()
  yaml_cfg = load_yaml_config(yaml_path)
  case = String(get(yaml_cfg, "case", DEFAULT_CASE))
  methods_cfg = haskey(yaml_cfg, "methods") ? _as_symbol_vec(yaml_cfg["methods"]) : Symbol[]
  methods = isempty(methods_cfg) ? DEFAULT_METHODS : methods_cfg

  print("\e[2J\e[H") # clear screen and move cursor to home position
  println("----------------------------------------------------------------------------------------------")
  println("Sparlectra version: ", Sparlectra.version(), "\n")
  println("Importing MATPOWER case file: $case\n")
  if !isempty(yaml_path)
    println("Using YAML config: $yaml_path\n")
  end

  timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
  logfile = joinpath(OUTDIR, "run_$(case)_$(timestamp).log")

  # This will:
  # - download case14.m into data/mpower/ if missing
  # - generate case14.jl if requested and missing
  local_case = Sparlectra.FetchMatpowerCase.ensure_casefile(case; outdir = Sparlectra.MPOWER_DIR, to_jl = true, overwrite = false)
  mpc = MatpowerIO.read_case(local_case)

  vmva_chk = MatpowerIO.vmva_power_mismatch_stats(mpc)
  if get(vmva_chk, :ok, false)
    println("==================== MATPOWER VM/VA self-check ====================")
    println("checked eqns     : P(PQ+PV)=", vmva_chk.n_p, "  Q(PQ)=", vmva_chk.n_q)
    println("max |ΔP| (PQ+PV) : ", vmva_chk.max_p_mis_pu, " pu  (", vmva_chk.max_p_mis_MW, " MW)")
    println("max |ΔQ| (PQ)    : ", vmva_chk.max_q_mis_pu, " pu  (", vmva_chk.max_q_mis_MVar, " MVar)")
    println("===================================================================\n")
  else
    println("MATPOWER VM/VA self-check skipped: ", get(vmva_chk, :msg, "unknown reason"))
  end

  cfg = bench_config_for_case(case, yaml_cfg)
  bench = bench_run_acpflow(
    casefile = basename(local_case),
    methods = methods,
    mpc = mpc,
    logfile = logfile,
    show_diff = show_diff,
    tol_vm = tol_vm,
    tol_va = tol_va,
    opt_fd = cfg.opt_fd,
    opt_sparse = cfg.opt_sparse,
    opt_flatstart = cfg.opt_flatstart,
    verbose = cfg.verbose,
    cooldown_iters = cfg.cooldown_iters,
    q_hyst_pu = cfg.q_hyst_pu,
    lock_pv_to_pq_buses = cfg.lock_pv_to_pq_buses,
    seconds = cfg.seconds,
    samples = cfg.samples,
    show_once = cfg.show_once,
  )
  return bench
end

main()

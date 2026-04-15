# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""
Compare load-flow behavior (convergence + runtime) for one or multiple MATPOWER cases.

Intended workflow to compare versions/states (e.g. r0.7.1 vs current):
1) Run this script in checkout A.
2) Run the same command in checkout B.
3) Compare the printed table.

Usage:
  julia --project=. src/examples/compare_stands_loadflow.jl
  julia --project=. src/examples/compare_stands_loadflow.jl case14.m case118.m

Env overrides:
  SPARLECTRA_COMPARE_REPEATS     (default: 3)
  SPARLECTRA_COMPARE_MAX_ITE     (default: 40)
  SPARLECTRA_COMPARE_TOL         (default: 1e-8)
  SPARLECTRA_COMPARE_METHOD      (default: rectangular)
  SPARLECTRA_COMPARE_OPT_SPARSE  (default: true)
  SPARLECTRA_COMPARE_OPT_FD      (default: false)
  SPARLECTRA_COMPARE_ENABLE_PQ_CONTROLLERS (default: true)
"""

using Sparlectra
using Statistics
using Printf

const DEFAULT_CASES = ["case14.m"]

_bool_env(name::String, default::Bool) = lowercase(get(ENV, name, string(default))) in ("1", "true", "yes", "on")

function run_case_benchmark(casefile::String;
                            repeats::Int,
                            max_ite::Int,
                            tol::Float64,
                            method::Symbol,
                            opt_sparse::Bool,
                            opt_fd::Bool,
                            enable_pq_gen_controllers::Bool)
  local_case = Sparlectra.FetchMatpowerCase.ensure_casefile(casefile; outdir = Sparlectra.MPOWER_DIR, to_jl = false, overwrite = false)

  elapsed_s = Float64[]
  iterations = Int[]
  converged = Bool[]

  for _ in 1:repeats
    status_ref = Ref{Any}(nothing)
    run_acpflow(
      casefile = basename(local_case),
      path = dirname(local_case),
      max_ite = max_ite,
      tol = tol,
      method = method,
      opt_sparse = opt_sparse,
      opt_fd = opt_fd,
      verbose = 0,
      show_results = false,
      show_compact_result = false,
      status_ref = status_ref,
      enable_pq_gen_controllers = enable_pq_gen_controllers,
    )

    st = status_ref[]
    push!(elapsed_s, st.elapsed_s)
    push!(iterations, st.iterations)
    push!(converged, st.converged)
  end

  return (
    casefile = casefile,
    runs = repeats,
    converged_runs = count(identity, converged),
    median_elapsed_s = median(elapsed_s),
    min_elapsed_s = minimum(elapsed_s),
    max_elapsed_s = maximum(elapsed_s),
    median_iterations = Int(round(median(iterations))),
    max_iterations = maximum(iterations),
  )
end

function main()
  cases = isempty(ARGS) ? DEFAULT_CASES : collect(ARGS)

  repeats = parse(Int, get(ENV, "SPARLECTRA_COMPARE_REPEATS", "3"))
  max_ite = parse(Int, get(ENV, "SPARLECTRA_COMPARE_MAX_ITE", "40"))
  tol = parse(Float64, get(ENV, "SPARLECTRA_COMPARE_TOL", "1e-8"))
  method = Symbol(get(ENV, "SPARLECTRA_COMPARE_METHOD", "rectangular"))
  opt_sparse = _bool_env("SPARLECTRA_COMPARE_OPT_SPARSE", true)
  opt_fd = _bool_env("SPARLECTRA_COMPARE_OPT_FD", false)
  enable_pq_gen_controllers = _bool_env("SPARLECTRA_COMPARE_ENABLE_PQ_CONTROLLERS", true)

  println("=== Sparlectra load-flow comparison run ===")
  println("Version: ", Sparlectra.version())
  println("Settings: repeats=", repeats,
          ", method=", method,
          ", max_ite=", max_ite,
          ", tol=", tol,
          ", opt_sparse=", opt_sparse,
          ", opt_fd=", opt_fd,
          ", enable_pq_gen_controllers=", enable_pq_gen_controllers)
  println("")

  results = NamedTuple[]
  for casefile in cases
    try
      push!(results, run_case_benchmark(casefile;
                                        repeats = repeats,
                                        max_ite = max_ite,
                                        tol = tol,
                                        method = method,
                                        opt_sparse = opt_sparse,
                                        opt_fd = opt_fd,
                                        enable_pq_gen_controllers = enable_pq_gen_controllers))
    catch err
      @warn "Case benchmark failed" casefile exception = (err, catch_backtrace())
    end
  end

  println("casefile | converged/runs | median_iter | max_iter | median_s | min_s | max_s")
  println("---------|----------------|-------------|----------|----------|-------|------")
  for r in results
    @printf("%s | %d/%d | %d | %d | %.4f | %.4f | %.4f\n",
            r.casefile, r.converged_runs, r.runs, r.median_iterations, r.max_iterations,
            r.median_elapsed_s, r.min_elapsed_s, r.max_elapsed_s)
  end
end

Base.invokelatest(main)

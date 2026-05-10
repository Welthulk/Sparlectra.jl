#!/usr/bin/env julia
# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

using Dates
using LinearAlgebra
using Printf
using Sparlectra

const DEFAULT_CONFIG = Dict{String,Any}(
  "bus_limits" => [100, 300, 500, 10000, 20000, 30000],
  "solver" => Dict{String,Any}("max_iter" => 25, "tol" => 1e-8, "verbose" => 0, "flatstart" => false, "method" => :rectangular, "opt_sparse" => true),
  "performance" => Dict{String,Any}("blas_threads" => 0, "show_thread_info" => true),
  "synthetic_network" => Dict{String,Any}("aspect_ratio" => 1.0, "base_mva" => 100.0, "r" => 0.01, "x" => 0.05, "g" => 0.0, "b" => 0.0, "load_mw_per_right_corner" => 50.0, "load_mvar_per_right_corner" => 15.0, "generation_balance" => 0.995, "vm_slack" => 1.0, "vm_flat" => 1.0),
)

function _copy_config_dict(x::Dict{String,Any})::Dict{String,Any}
  out = Dict{String,Any}()
  for (key, value) in x
    out[key] = value isa Dict{String,Any} ? _copy_config_dict(value) : value
  end
  return out
end

_get(d::Dict{String,Any}, key::String, default) = haskey(d, key) ? d[key] : default
_as_float(x) = x isa Real ? Float64(x) : parse(Float64, String(x))
_as_int(x) = x isa Integer ? Int(x) : parse(Int, String(x))
_as_symbol(x) = x isa Symbol ? x : Symbol(String(x))

function _parse_args(args::Vector{String})
  config_path = nothing
  limits = Int[]
  max_override = nothing
  for arg in args
    if startswith(arg, "--max-buses=")
      max_override = parse(Int, split(arg, "=", limit = 2)[2])
    elseif _looks_like_yaml_path(arg)
      config_path = arg
    else
      push!(limits, parse(Int, arg))
    end
  end
  return config_path, limits, max_override
end

function _looks_like_yaml_path(arg::AbstractString)::Bool
  lower = lowercase(arg)
  return endswith(lower, ".yaml") || endswith(lower, ".yml") || endswith(lower, ".yaml.example") || endswith(lower, ".yml.example")
end

function _fallback_example_path(path::AbstractString)::String
  lower = lowercase(path)
  if endswith(lower, ".yaml") || endswith(lower, ".yml")
    return path * ".example"
  end
  return path
end

function _resolve_config_path(config_path)::Union{Nothing,String}
  if isnothing(config_path)
    println("No configuration file specified; using built-in defaults.")
    return nothing
  end

  if isfile(config_path)
    return String(config_path)
  end
  fallback_path = _fallback_example_path(config_path)
  if fallback_path != config_path && isfile(fallback_path)
    println("Configuration file not found: ", config_path)
    println("Using example configuration fallback: ", fallback_path)
    return fallback_path
  end

  println("Configuration file not found: ", config_path)
  println("No fallback configuration file found; using built-in defaults.")
  return nothing
end

function _configure_blas!(cfg::Dict{String,Any})
  perf = _get(cfg, "performance", Dict{String,Any}())
  blas_threads = _as_int(_get(perf, "blas_threads", 0))
  if blas_threads > 0
    BLAS.set_num_threads(blas_threads)
  end
  if as_bool(_get(perf, "show_thread_info", true))
    println("BLAS threads: ", BLAS.get_num_threads())
    println("Julia threads: ", Threads.nthreads())
  end
end

function _max_abs_mismatch(net::Net)
  V = buildVoltageVector(net)
  Y = createYBUS(net = net, sparse = false, printYBUS = false)
  Sspec = buildComplexSVec(net)
  Scalc = V .* conj(Y * V)
  mismatch = Sspec - Scalc
  for slack in net.slackVec
    mismatch[slack] = 0.0 + 0.0im
  end
  return maximum(abs.(real.(mismatch))), maximum(abs.(imag.(mismatch)))
end

function _network_kwargs(cfg::Dict{String,Any})
  sn = _get(cfg, "synthetic_network", Dict{String,Any}())
  return (
    aspect_ratio = _as_float(_get(sn, "aspect_ratio", 1.0)),
    base_mva = _as_float(_get(sn, "base_mva", 100.0)),
    r = _as_float(_get(sn, "r", 0.01)),
    x = _as_float(_get(sn, "x", 0.05)),
    g = _as_float(_get(sn, "g", 0.0)),
    b = _as_float(_get(sn, "b", 0.0)),
    load_mw_per_right_corner = _as_float(_get(sn, "load_mw_per_right_corner", 50.0)),
    load_mvar_per_right_corner = _as_float(_get(sn, "load_mvar_per_right_corner", 15.0)),
    generation_balance = _as_float(_get(sn, "generation_balance", 0.995)),
    vm_slack = _as_float(_get(sn, "vm_slack", 1.0)),
    vm_flat = _as_float(_get(sn, "vm_flat", 1.0)),
  )
end

function _print_ascii_plot(rows)
  isempty(rows) && return nothing
  max_ms = maximum(row.solve_ms for row in rows)
  max_ms <= 0.0 && return nothing
  println("\nASCII plot: nbus versus solve_ms")
  for row in rows
    width = max(1, round(Int, 40 * row.solve_ms / max_ms))
    @printf("%8d | %s %.3f ms\n", row.nbus, repeat("█", width), row.solve_ms)
  end
  return nothing
end

function main(args = ARGS)
  config_path, cli_limits, max_override = _parse_args(collect(args))
  cfg = _copy_config_dict(DEFAULT_CONFIG)
  resolved_config_path = _resolve_config_path(config_path)
  if !isnothing(resolved_config_path)
    println("Loading configuration: ", resolved_config_path)
    merge_yaml_dict!(cfg, load_yaml_dict(resolved_config_path))
  end
  if !isempty(cli_limits)
    cfg["bus_limits"] = cli_limits
  end
  if !isnothing(max_override)
    cfg["bus_limits"] = [max_override]
  end

  _configure_blas!(cfg)
  bus_limits = as_int_vector(cfg["bus_limits"])
  solver = _get(cfg, "solver", Dict{String,Any}())
  max_iter = _as_int(_get(solver, "max_iter", 25))
  tol = _as_float(_get(solver, "tol", 1e-8))
  verbose = _as_int(_get(solver, "verbose", 0))
  flatstart = as_bool(_get(solver, "flatstart", false))
  method = _as_symbol(_get(solver, "method", :rectangular))
  opt_sparse = as_bool(_get(solver, "opt_sparse", true))

  outdir = joinpath(@__DIR__, "_out")
  mkpath(outdir)
  logfile = joinpath(outdir, "synthetic_tiled_grid_pf_perf_" * Dates.format(now(), "yyyymmdd_HHMMSS") * ".log")

  rows = NamedTuple[]
  header = "Limit nbus rows cols branch_count converged iterations solve_ms max|ΔP| max|ΔQ| system_build_ms system_build_allocated_bytes total_case_runtime_ms"
  open(logfile, "w") do io
    println(io, "# Sparlectra synthetic tiled-grid PF benchmark")
    println(io, "# timestamp = ", now())
    println(io, header)

    println("\n", header)
    for limit in bus_limits
      build_alloc = 0
      net = nothing
      meta = nothing
      build_time = @elapsed begin
        build_alloc = @allocated begin
          net, meta = build_synthetic_tiled_grid_net(limit; _network_kwargs(cfg)...)
        end
      end
      solve_time = @elapsed begin
        iterations, erg, solver_elapsed = run_net_acpflow(net = net, max_ite = max_iter, tol = tol, verbose = verbose, opt_sparse = opt_sparse, method = method, show_results = false, opt_flatstart = flatstart)
      end
      max_dp, max_dq = _max_abs_mismatch(net)
      total_runtime_ms = 1000.0 * (build_time + solve_time)
      row = (
        limit = limit,
        nbus = meta.actual_buses,
        rows = meta.rows,
        cols = meta.cols,
        branch_count = meta.branch_count,
        converged = (erg == 0),
        iterations = iterations,
        solve_ms = 1000.0 * solver_elapsed,
        max_dp = max_dp,
        max_dq = max_dq,
        system_build_ms = 1000.0 * build_time,
        system_build_allocated_bytes = build_alloc,
        total_case_runtime_ms = total_runtime_ms,
      )
      push!(rows, row)
      @printf(
        "%5d %5d %4d %4d %12d %9s %10d %8.3f %.3e %.3e %15.3f %28d %21.3f\n",
        row.limit,
        row.nbus,
        row.rows,
        row.cols,
        row.branch_count,
        string(row.converged),
        row.iterations,
        row.solve_ms,
        row.max_dp,
        row.max_dq,
        row.system_build_ms,
        row.system_build_allocated_bytes,
        row.total_case_runtime_ms
      )
      @printf(io, "%5d %5d %4d %4d %12d %9s %10d %8.3f %.3e %.3e %15.3f %28d %21.3f\n", row.limit, row.nbus, row.rows, row.cols, row.branch_count, row.converged, row.iterations, row.solve_ms, row.max_dp, row.max_dq, row.system_build_ms, row.system_build_allocated_bytes, row.total_case_runtime_ms)
    end
  end

  _print_ascii_plot(rows)
  println("\nWrote log file: ", logfile)
  return rows
end

# IMPORTANT for Julia 1.12 / Revise world-age safety:
# Default: run immediately (also when included from REPL), like a script.
if get(ENV, "SPARLECTRA_SUITE_NO_AUTORUN", "0") != "1"
  main_fn = getfield(@__MODULE__, :main)
  Base.invokelatest(main_fn, ARGS)
  nothing
end

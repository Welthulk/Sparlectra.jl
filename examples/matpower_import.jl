# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

using Sparlectra

function _parse_cli_args(args::Vector{String})
  casefile = ""
  julia_threads_override = ""
  passthrough = String[]
  for arg in args
    if startswith(arg, "--julia-threads=")
      julia_threads_override = split(arg, "=", limit = 2)[2]
    elseif startswith(arg, "--")
      push!(passthrough, arg)
    elseif isempty(casefile)
      casefile = arg
    else
      push!(passthrough, arg)
    end
  end
  return (; casefile = casefile, julia_threads_override = julia_threads_override, passthrough = passthrough)
end

function _resolve_julia_threads_request(config_file::AbstractString, cli_override::AbstractString)
  !isempty(strip(cli_override)) && return Sparlectra.parse_runtime_threads_request(cli_override)
  env_override = get(ENV, "SPARLECTRA_JULIA_THREADS", "")
  !isempty(strip(env_override)) && return Sparlectra.parse_runtime_threads_request(env_override)
  cfg = Sparlectra.load_sparlectra_config(config_file; reload = true)
  return Sparlectra.parse_runtime_threads_request(cfg.runtime.julia_threads)
end

function _ensure_julia_threads_for_script(config_file::AbstractString, cli_override::AbstractString)
  get(ENV, "SPARLECTRA_THREADS_REEXEC", "0") == "1" && return nothing
  requested = _resolve_julia_threads_request(config_file, cli_override)
  isnothing(requested) && return nothing
  requested == Threads.nthreads() && return nothing
  program_file = abspath(PROGRAM_FILE)
  isempty(program_file) && return nothing
  # Only re-exec when this file is the actual script entry point.
  # In VS Code / REPL workflows PROGRAM_FILE can point to terminalserver.jl.
  program_file == abspath(@__FILE__) || return nothing
  cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) --threads=$(requested) $(program_file) $(ARGS...)`
  withenv("SPARLECTRA_THREADS_REEXEC" => "1") do
    run(cmd)
  end
  exit(0)
end

function main()
  parsed = _parse_cli_args(copy(ARGS))
  config_file = Sparlectra.configuration_path_from_inputs(
    env_var = "SPARLECTRA_CONFIGURATION_YAML",
    fallback_paths = [Sparlectra.USER_SPARLECTRA_CONFIG_PATH],
  )
  _ensure_julia_threads_for_script(config_file, parsed.julia_threads_override)
  return Sparlectra.run_matpower_case(; config_file = config_file, casefile = parsed.casefile)
end

if get(ENV, "SPARLECTRA_MATPOWER_IMPORT_NO_MAIN", "0") != "1"
  try
    Base.invokelatest(getfield(@__MODULE__, :main))
  catch err
    if err isa InterruptException
      println(stderr, "\nInterrupted by user.")
      exit(130)
    end
    rethrow()
  end
end

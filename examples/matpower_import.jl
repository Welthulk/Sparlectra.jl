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

# Date: 2026-05-18
# file: examples/matpower_import.jl
# purpose: CLI runner that imports a MATPOWER case via run_matpower_case using the resolved Sparlectra configuration, with optional Julia-thread re-exec

using Sparlectra

include(joinpath(@__DIR__, "internal", "example_header.jl"))

"""
    _parse_cli_args(args)

Parse MATPOWER-runner CLI arguments into case-file, Julia-thread override, and
passthrough flags.
"""
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

"""
    _resolve_julia_threads_request(config_file, cli_override)

Resolve requested Julia thread count from CLI override, environment, or YAML
configuration (in that priority order).
"""
function _resolve_julia_threads_request(config_file::AbstractString, cli_override::AbstractString)
  !isempty(strip(cli_override)) && return Sparlectra.parse_runtime_threads_request(cli_override)
  env_override = get(ENV, "SPARLECTRA_JULIA_THREADS", "")
  !isempty(strip(env_override)) && return Sparlectra.parse_runtime_threads_request(env_override)
  cfg = Sparlectra.load_sparlectra_config(config_file; reload = true)
  return Sparlectra.parse_runtime_threads_request(cfg.runtime.julia_threads)
end

"""
    _ensure_julia_threads_for_script(config_file, cli_override)

Re-execute the current script with the requested Julia thread count when needed.

Returns `nothing` when no re-exec is required. This mutates process control flow
and exits the current process after handing off to the re-executed process.
"""
function _ensure_julia_threads_for_script(config_file::AbstractString, cli_override::AbstractString)
  get(ENV, "SPARLECTRA_THREADS_REEXEC", "0") == "1" && return nothing
  requested = Base.invokelatest(getfield(@__MODULE__, :_resolve_julia_threads_request), config_file, cli_override)
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

"""
    main()

Entry point for the MATPOWER import runner.
"""
function main()
  print_example_banner("examples/matpower_import.jl", "CLI runner that imports a MATPOWER case via run_matpower_case using the resolved Sparlectra configuration, with optional Julia-thread re-exec")
  parsed = _parse_cli_args(copy(ARGS))
  config_file = Sparlectra.configuration_path_from_inputs(
    env_var = "SPARLECTRA_CONFIGURATION_YAML",
    fallback_paths = [Sparlectra.USER_SPARLECTRA_CONFIG_PATH],
  )
  Base.invokelatest(getfield(@__MODULE__, :_ensure_julia_threads_for_script), config_file, parsed.julia_threads_override)
  return Sparlectra.run_matpower_case(; config_file = config_file, casefile = parsed.casefile)
end

# SPARLECTRA_MATPOWER_IMPORT_NO_MAIN=1 is a test-only escape hatch: it lets
# test/test_matpower_example.jl load this file's function definitions into a
# throwaway module (to smoke-test its structure) without actually running the
# CLI import. It has no effect on normal console/interactive use — those
# always run main() unconditionally, exactly like every other examples/*.jl
# program (see examples/exp_transformer_tap_changer_model.jl for why no
# `abspath(PROGRAM_FILE) == @__FILE__` guard is used here).
if get(ENV, "SPARLECTRA_MATPOWER_IMPORT_NO_MAIN", "0") != "1"
  try
    _ = run_example(main)
    nothing
  catch err
    if err isa InterruptException
      println(stderr, "\nInterrupted by user.")
      exit(130)
    end
    rethrow()
  end
end


# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

using Sparlectra

# ---------------------------------------------------------------------------
# VS Code defaults
# ---------------------------------------------------------------------------
#
# These defaults are used when the script is started without command-line
# arguments. They make it convenient to press "Run" in VS Code.
#
# For normal command-line usage, pass the case and configs explicitly:
#
#   julia --project=. examples/matpower_import_multi_config.jl \
#     data/mpower/case_ACTIVSg10k.m \
#     --config=examples/configuration_activsg10k_wrong_branch_warn.yaml \
#     --config=examples/configuration_activsg10k_wrong_branch_branch20_warn.yaml \
#     --config=examples/configuration_activsg10k_wrong_branch_branch20_fail.yaml \
#     --status-only
#
#=
bash:

cd ~/.julia/dev/Sparlectra

julia --project=. examples/matpower_import_multi_config.jl \
  data/mpower/case_ACTIVSg10k.m \
  --config=examples/configuration_activsg10k_wrong_branch_warn.yaml \
  --config=examples/configuration_activsg10k_wrong_branch_branch20_warn.yaml \
  --config=examples/configuration_activsg10k_wrong_branch_branch20_fail.yaml \
  --status-only

  =#

const VS_CODE_CASEFILE = "data/mpower/case_ACTIVSg10k.m"

const VS_CODE_CONFIG_FILES = String["examples/configuration_activsg10k_wrong_branch_warn.yaml", "examples/configuration_activsg10k_wrong_branch_branch20_warn.yaml", "examples/configuration_activsg10k_wrong_branch_branch20_fail.yaml"]

# :runner      -> use Sparlectra.run_matpower_case, like the original script
# :status_only -> run the PF path and print rectangular wrong-branch status fields
const VS_CODE_MODE = :status_only

"""
    ParsedArgs

Command-line argument container for the MATPOWER multi-configuration runner.
"""
struct ParsedArgs
  casefile::String
  config_files::Vector{String}
  julia_threads_override::String
  mode::Symbol
  show_help::Bool
end

"""
    _print_help()

Print supported command-line usage.
"""
function _print_help()
  println("""
MATPOWER import runner with multi-configuration support.

Usage:
  julia --project=. examples/matpower_import_multi_config.jl CASEFILE [options]

Options:
  --config=PATH              Add one YAML configuration file. May be repeated.
  --configs=A,B,C            Add multiple YAML configuration files.
  --status-only              Run a compact status check and print wrong-branch fields.
  --runner                   Use Sparlectra.run_matpower_case output mode.
  --julia-threads=N          Request Julia thread count, with script re-exec when possible.
  --help                     Show this help text.

Examples:
  julia --project=. examples/matpower_import_multi_config.jl \\
    data/mpower/case_ACTIVSg10k.m \\
    --config=examples/configuration_activsg10k_wrong_branch_warn.yaml \\
    --config=examples/configuration_activsg10k_wrong_branch_branch20_warn.yaml \\
    --config=examples/configuration_activsg10k_wrong_branch_branch20_fail.yaml \\
    --status-only
""")
  return nothing
end

"""
    _split_config_list(raw)

Split a comma- or semicolon-separated config list into paths.
"""
function _split_config_list(raw::AbstractString)::Vector{String}
  parts = split(raw, r"[,;]")
  return String[strip(p) for p in parts if !isempty(strip(p))]
end

"""
    _parse_cli_args(args)

Parse MATPOWER-runner CLI arguments into case-file, config files, mode, and
Julia-thread override.
"""
function _parse_cli_args(args::Vector{String})
  casefile = ""
  config_files = String[]
  julia_threads_override = ""
  mode = :runner
  show_help = false

  for arg in args
    if arg == "--help" || arg == "-h"
      show_help = true
    elseif startswith(arg, "--julia-threads=")
      julia_threads_override = split(arg, "=", limit = 2)[2]
    elseif startswith(arg, "--config=")
      cfg = strip(split(arg, "=", limit = 2)[2])
      isempty(cfg) && throw(ArgumentError("`--config=` requires a non-empty path."))
      push!(config_files, cfg)
    elseif startswith(arg, "--configs=")
      cfgs = strip(split(arg, "=", limit = 2)[2])
      append!(config_files, _split_config_list(cfgs))
    elseif arg == "--status-only"
      mode = :status_only
    elseif arg == "--runner"
      mode = :runner
    elseif startswith(arg, "--case=")
      casefile = strip(split(arg, "=", limit = 2)[2])
    elseif startswith(arg, "--")
      throw(ArgumentError("Unsupported flag: $(arg)"))
    elseif isempty(casefile)
      casefile = arg
    else
      throw(ArgumentError("Only one positional case file is supported. Extra argument: $(arg)"))
    end
  end

  return ParsedArgs(casefile, config_files, julia_threads_override, mode, show_help)
end

"""
    _default_config_file()

Resolve the default configuration file using the standard Sparlectra mechanism.
"""
function _default_config_file()::String
  return Sparlectra.configuration_path_from_inputs(env_var = "SPARLECTRA_CONFIGURATION_YAML", fallback_paths = [Sparlectra.USER_SPARLECTRA_CONFIG_PATH])
end

"""
    _resolve_casefile(parsed)

Resolve case-file from CLI arguments or VS Code defaults.
"""
function _resolve_casefile(parsed::ParsedArgs)::String
  !isempty(strip(parsed.casefile)) && return parsed.casefile
  !isempty(strip(VS_CODE_CASEFILE)) && return VS_CODE_CASEFILE
  throw(ArgumentError("No MATPOWER case provided. Pass a casefile or set VS_CODE_CASEFILE."))
end

"""
    _resolve_config_files(parsed)

Resolve configuration files from CLI arguments, environment/default lookup, or
VS Code defaults.
"""
function _resolve_config_files(parsed::ParsedArgs)::Vector{String}
  !isempty(parsed.config_files) && return parsed.config_files

  if !isempty(VS_CODE_CONFIG_FILES)
    return copy(VS_CODE_CONFIG_FILES)
  end

  cfg = _default_config_file()
  return isempty(cfg) ? String[""] : String[cfg]
end

"""
    _resolve_mode(parsed)

Resolve run mode from CLI arguments or VS Code default.
"""
function _resolve_mode(parsed::ParsedArgs)::Symbol
  isempty(ARGS) && return VS_CODE_MODE
  return parsed.mode
end

"""
    _resolve_julia_threads_request(config_file, cli_override)

Resolve requested Julia thread count from CLI override, environment, or YAML
configuration, in that priority order.
"""
function _resolve_julia_threads_request(config_file::AbstractString, cli_override::AbstractString)
  !isempty(strip(cli_override)) && return Sparlectra.parse_runtime_threads_request(cli_override)
  env_override = get(ENV, "SPARLECTRA_JULIA_THREADS", "")
  !isempty(strip(env_override)) && return Sparlectra.parse_runtime_threads_request(env_override)
  isempty(strip(config_file)) && return nothing
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
    _split_casefile_for_run_acpflow(casefile)

Convert a possibly path-qualified case file into `(casefile, path)` arguments
for `run_acpflow`.
"""
function _split_casefile_for_run_acpflow(casefile::AbstractString)
  dir = dirname(casefile)
  file = basename(casefile)

  if isempty(dir) || dir == "."
    return (casefile = file, path = nothing)
  end

  return (casefile = file, path = dir)
end

"""
    _print_config_header(index, total, config_file, casefile, mode)

Print a compact per-run header.
"""
function _print_config_header(index::Int, total::Int, config_file::AbstractString, casefile::AbstractString, mode::Symbol)
  println()
  println("="^96)
  println("MATPOWER configuration run ", index, " / ", total)
  println("mode        : ", mode)
  println("casefile    : ", casefile)
  println("config file : ", isempty(config_file) ? "default/none" : config_file)
  println("="^96)
  return nothing
end

"""
    _print_rectangular_status(st)

Print the rectangular solver status fields relevant for wrong-branch diagnostics.
"""
function _print_rectangular_status(st)
  if isnothing(st)
    println("rectangular status             = unavailable")
    return nothing
  end

  println("status                         = ", st.status)
  println("numerical_converged            = ", st.numerical_converged)
  println("final_converged                = ", st.final_converged)
  println("reason                         = ", st.reason)
  println("final_mismatch                 = ", st.final_mismatch)
  println("wrong_branch_detection         = ", st.wrong_branch_detection)
  println("wrong_branch_status            = ", st.wrong_branch_status)
  println("wrong_branch_reason            = ", st.wrong_branch_reason)
  println("wrong_branch_low_vm_count      = ", st.wrong_branch_low_vm_count)
  println("wrong_branch_high_vm_count     = ", st.wrong_branch_high_vm_count)
  println("wrong_branch_angle_spread_deg  = ", st.wrong_branch_angle_spread_deg)
  println("wrong_branch_branch_violations = ", st.wrong_branch_branch_angle_violation_count)
  println("wrong_branch_worst_branch_deg  = ", st.wrong_branch_worst_branch_angle_deg)
  println("branch_quality_metrics         = ", st.branch_quality_metrics)
  return nothing
end

"""
    _run_status_only(config_file, casefile)

Run the public AC power-flow path and print the rectangular wrong-branch status.
"""
function _run_status_only(config_file::AbstractString, casefile::AbstractString)
  cfg = isempty(strip(config_file)) ? Sparlectra.load_sparlectra_config(; reload = true) : Sparlectra.load_sparlectra_config(config_file; reload = true)

  println("configured detection           = ", cfg.powerflow.wrong_branch_detection)
  println("configured branch angle limit  = ", cfg.powerflow.wrong_branch_max_branch_angle_deg)

  case_args = _split_casefile_for_run_acpflow(casefile)

  net = run_acpflow(casefile = case_args.casefile, path = case_args.path, config = cfg, show_results = false, printResultToFile = false, verbose = 0)

  st = Sparlectra.rectangular_pf_status(net)
  _print_rectangular_status(st)
  return st
end

"""
    _run_runner(config_file, casefile)

Run the standard Sparlectra MATPOWER runner.
"""
function _run_runner(config_file::AbstractString, casefile::AbstractString)
  return Sparlectra.run_matpower_case(; config_file = config_file, casefile = casefile)
end

"""
    main()

Entry point for the MATPOWER multi-configuration runner.
"""
function main()
  parsed = _parse_cli_args(copy(ARGS))

  if parsed.show_help
    _print_help()
    return nothing
  end

  casefile = _resolve_casefile(parsed)
  config_files = _resolve_config_files(parsed)
  mode = _resolve_mode(parsed)

  isempty(config_files) && throw(ArgumentError("No configuration files resolved."))

  # Use the first config for runtime-thread resolution. Running multiple configs
  # with different Julia thread requests in one process is intentionally avoided.
  Base.invokelatest(getfield(@__MODULE__, :_ensure_julia_threads_for_script), first(config_files), parsed.julia_threads_override)

  for (i, cfgfile) in enumerate(config_files)
    _print_config_header(i, length(config_files), cfgfile, casefile, mode)

    if mode === :runner
      _run_runner(cfgfile, casefile)
    elseif mode === :status_only
      _run_status_only(cfgfile, casefile)
    else
      throw(ArgumentError("Unsupported mode: $(mode). Use :runner or :status_only."))
    end
  end

  return nothing
end

if get(ENV, "SPARLECTRA_MATPOWER_IMPORT_NO_MAIN", "0") != "1"
  try
    # Use invokelatest for script-style execution under Revise/Julia 1.12.
    _ = Base.invokelatest(getfield(@__MODULE__, :main))
    nothing
  catch err
    if err isa InterruptException
      println(stderr, "\nInterrupted by user.")
      exit(130)
    end
    rethrow()
  end
end

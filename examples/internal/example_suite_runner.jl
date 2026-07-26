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

# Date: 2026-07-26
# file: examples/internal/example_suite_runner.jl
# purpose: shared infrastructure for the example suite runners (registry struct, CLI parsing, subprocess execution, CSV/Markdown summaries)

# Shared backend for examples/run_powerflow_suite.jl,
# examples/run_state_estimation_suite.jl and examples/run_others_suite.jl.
#
# Each registered example runs in a fresh Julia subprocess with the same
# project, reproducing the documented invocation
# `julia --project=. examples/<file>.jl [args...]`. A subprocess (instead of
# `include` into the suite session) is required because several examples call
# `exit(...)` (matpower_import.jl even re-execs itself with a different
# thread count), read the process-global `ARGS`, and define top-level
# `const`s.
#
# This file only defines functions and types; it has no side effects when
# included. It deliberately does not load Sparlectra so that `--help` and
# `--list` respond without package load time.

using Dates

include(joinpath(@__DIR__, "example_header.jl"))

"""
    ExampleSpec

One registered example program of a suite.

- `name`: registry key used for `--only=`/`--skip=` and log/summary rows.
- `file`: path relative to `examples/`.
- `purpose`: one line, mirroring the example file's `# purpose:` header.
- `args`: extra command-line arguments passed to the subprocess.
- `heavy`: skipped unless `--include-heavy` (long runtime or very large cases).
- `optional`: skipped unless `--include-optional` (optional dependency or
  experimental API).
- `requires_package`: name of a package that must be installed; when it is not
  resolvable the example is reported as `skipped_missing_dependency`.
- `requires_config`: the example needs a local user configuration
  (`examples/configuration.yaml` or `SPARLECTRA_CONFIGURATION_YAML`); without
  one it is reported as `skipped_missing_config`.
- `timeout_s`: per-example wall-clock timeout (overridable via `--timeout=`).
"""
Base.@kwdef struct ExampleSpec
  name::String
  file::String
  purpose::String
  args::Vector{String} = String[]
  heavy::Bool = false
  optional::Bool = false
  requires_package::Union{Nothing,String} = nothing
  requires_config::Bool = false
  timeout_s::Int = 600
end

function _parse_bool(value::AbstractString)
  normalized = lowercase(strip(value))
  normalized in ("true", "1", "yes", "on") && return true
  normalized in ("false", "0", "no", "off") && return false
  throw(ArgumentError("invalid boolean value: $value"))
end

function _default_suite_options(suite_name::AbstractString)
  return Dict{String,Any}(
    "only" => "",
    "skip" => "",
    "include-heavy" => false,
    "include-optional" => false,
    "continue-on-error" => true,
    "strict" => false,
    "quiet" => false,
    "list" => false,
    "timeout" => 0,
    "write-csv" => true,
    "write-markdown" => true,
    "output-dir" => joinpath(@__DIR__, "..", "_out", suite_name),
  )
end

_spec_tags(spec::ExampleSpec) = join(filter(!isempty, [spec.heavy ? "heavy" : "", spec.optional ? "optional" : "", spec.requires_package === nothing ? "" : "requires:$(spec.requires_package)", spec.requires_config ? "requires:user-config" : ""]), " ")

function _print_suite_help(suite_name::AbstractString, title::AbstractString, specs::Vector{ExampleSpec}, notes::Vector{String})
  println("""
$title

Runs each registered example in a fresh Julia subprocess and writes a
summary (CSV + Markdown) plus per-example logs to the output directory.

Usage:
  julia --project=. examples/run_$(suite_name).jl [options]

Options:
  --list
      Print the registry (names, tags, purposes) and exit.
  --only=<name,...>       Run only the named examples.
  --skip=<name,...>       Skip the named examples.
  --include-heavy         Also run examples tagged heavy (long runtime / very large cases).
  --include-optional      Also run examples tagged optional (optional dependency / experimental API).
  --timeout=<seconds>     Global per-example timeout override (0 keeps per-example defaults).
  --continue-on-error=true|false   (default true)
  --strict                Error at the end if any example failed or timed out.
  --write-csv=true|false           (default true)
  --write-markdown=true|false      (default true)
  --output-dir=<path>     (default examples/_out/$(suite_name)/)
  --quiet
  --help
""")
  println("Registered examples:")
  for spec in specs
    tags = _spec_tags(spec)
    println("  ", rpad(spec.name, 34), isempty(tags) ? "" : "[$tags] ", spec.purpose)
  end
  for note in notes
    println("\n", note)
  end
  return nothing
end

function _print_suite_list(specs::Vector{ExampleSpec})
  println("| name | file | tags | timeout_s | purpose |")
  println("|---|---|---|---:|---|")
  for spec in specs
    println("| ", spec.name, " | ", spec.file, " | ", _spec_tags(spec), " | ", spec.timeout_s, " | ", _md(spec.purpose), " |")
  end
  return nothing
end

function parse_example_suite_cli(suite_name::AbstractString, title::AbstractString, specs::Vector{ExampleSpec}, notes::Vector{String}, args = ARGS)
  opt = _default_suite_options(suite_name)
  boolean_keys = Set(["continue-on-error", "strict", "quiet", "include-heavy", "include-optional", "list", "write-csv", "write-markdown"])
  for arg in args
    if arg in ("--help", "-h")
      _print_suite_help(suite_name, title, specs, notes)
      return nothing
    elseif arg in ("--list", "--quiet", "--strict", "--include-heavy", "--include-optional")
      opt[arg[3:end]] = true
      continue
    elseif arg == "--no-continue-on-error"
      opt["continue-on-error"] = false
      continue
    end

    startswith(arg, "--") || throw(ArgumentError("unsupported argument: $arg"))
    parts = split(arg[3:end], "="; limit = 2)
    length(parts) == 2 || throw(ArgumentError("expected --key=value, got: $arg"))
    key, value = parts
    haskey(opt, key) || throw(ArgumentError("unknown option --$key"))

    if key == "timeout"
      opt[key] = parse(Int, value)
    elseif key in boolean_keys
      opt[key] = _parse_bool(value)
    else
      opt[key] = String(value)
    end
  end
  return opt
end

function _split_name_list(value::AbstractString, specs::Vector{ExampleSpec}, option::AbstractString)
  names = [strip(part) for part in split(value, ",") if !isempty(strip(part))]
  known = Set(spec.name for spec in specs)
  for name in names
    name in known || throw(ArgumentError("unknown example name '$name' in --$option (see --list)"))
  end
  return Set(String.(names))
end

function _package_available(name::AbstractString)
  pkgid = Base.identify_package(name)
  pkgid === nothing && return false
  return Base.locate_package(pkgid) !== nothing
end

function _user_config_available()
  isempty(strip(get(ENV, "SPARLECTRA_CONFIGURATION_YAML", ""))) || return true
  return isfile(normpath(joinpath(@__DIR__, "..", "configuration.yaml")))
end

# Keep in sync with the identical helpers in examples/validate_dtf_suite.jl.
function _csv_cell(value)
  value isa Missing && return ""
  value === nothing && return ""
  text = string(value)
  if occursin(',', text) || occursin('"', text) || occursin('\n', text) || occursin('\r', text)
    return "\"" * replace(text, "\"" => "\"\"") * "\""
  end
  return text
end

function _write_csv(path::AbstractString, rows)
  mkpath(dirname(path))
  open(path, "w") do io
    isempty(rows) && return
    columns = propertynames(first(rows))
    println(io, join(string.(columns), ","))
    for row in rows
      println(io, join((_csv_cell(getproperty(row, column)) for column in columns), ","))
    end
  end
  return path
end

function _md(value)
  value isa Missing && return ""
  value === nothing && return ""
  return replace(string(value), "|" => "\\|", "\n" => " ", "\r" => " ")
end

const _SUITE_SKIP_STATUSES = ("skipped_heavy", "skipped_optional", "skipped_missing_dependency", "skipped_missing_config", "skipped_by_filter")

_is_failure_status(status::AbstractString) = !(status == "ok" || status in _SUITE_SKIP_STATUSES)

function _summary_row(spec::ExampleSpec; status, exit_code = missing, duration_s = missing, log_file = missing, message = "")
  return (
    name = spec.name,
    file = spec.file,
    status = status,
    exit_code = exit_code,
    duration_s = duration_s isa Missing ? missing : round(duration_s; digits = 1),
    tags = _spec_tags(spec),
    log_file = log_file,
    message = message,
  )
end

"""
    _spec_skip_status(spec, opt, only_names, skip_names) -> Union{Nothing,Tuple{String,String}}

Returns `nothing` when the example should run, otherwise `(status, message)`.
"""
function _spec_skip_status(spec::ExampleSpec, opt, only_names, skip_names)
  only_names !== nothing && !(spec.name in only_names) && return ("skipped_by_filter", "not selected by --only")
  spec.name in skip_names && return ("skipped_by_filter", "excluded by --skip")
  spec.heavy && !opt["include-heavy"] && return ("skipped_heavy", "run with --include-heavy")
  spec.optional && !opt["include-optional"] && return ("skipped_optional", "run with --include-optional")
  if spec.requires_package !== nothing && !_package_available(spec.requires_package)
    return ("skipped_missing_dependency", "package $(spec.requires_package) is not installed")
  end
  if spec.requires_config && !_user_config_available()
    return ("skipped_missing_config", "create examples/configuration.yaml (see examples/README.md) or set SPARLECTRA_CONFIGURATION_YAML")
  end
  return nothing
end

function _run_spec_subprocess(spec::ExampleSpec, opt, examples_dir::AbstractString, logs_dir::AbstractString)
  file_path = normpath(joinpath(examples_dir, spec.file))
  log_path = joinpath(logs_dir, spec.name * ".log")
  timeout_s = opt["timeout"] > 0 ? opt["timeout"] : spec.timeout_s
  cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $(file_path) $(spec.args...)`
  # Examples are documented to run from the repository root
  # (`julia --project=. examples/<name>.jl`); several resolve data files and
  # configuration relative to the working directory.
  cmd = Cmd(cmd; dir = normpath(joinpath(examples_dir, "..")))

  started = time()
  timed_out = false
  proc = open(log_path, "w") do io
    println(io, "# command: ", join(cmd.exec, " "))
    println(io, "# started: ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"), "\n")
    flush(io)
    p = run(pipeline(cmd; stdout = io, stderr = io); wait = false)
    if timedwait(() -> process_exited(p), float(timeout_s); pollint = 0.5) == :timed_out
      timed_out = true
      kill(p)
      if timedwait(() -> process_exited(p), 10.0; pollint = 0.5) == :timed_out && Sys.isunix()
        kill(p, Base.SIGKILL)
      end
    end
    wait(p)
    p
  end
  duration = time() - started

  if timed_out
    return _summary_row(spec; status = "timeout", exit_code = proc.exitcode, duration_s = duration, log_file = log_path, message = "killed after $(timeout_s) s")
  end
  status = success(proc) ? "ok" : "failed"
  message = status == "ok" ? "" : "exit code $(proc.exitcode), see log"
  return _summary_row(spec; status = status, exit_code = proc.exitcode, duration_s = duration, log_file = log_path, message = message)
end

function _write_suite_markdown(path::AbstractString, title::AbstractString, opt, rows, notes::Vector{String})
  mkpath(dirname(path))
  open(path, "w") do io
    println(io, "# ", title, "\n")
    println(io, "- generated: ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
    println(io, "- output directory: `", opt["output-dir"], "`")
    println(io, "- include heavy: `", opt["include-heavy"], "`; include optional: `", opt["include-optional"], "`\n")

    println(io, "## Results\n")
    println(io, "| Name | File | Status | Exit code | Duration s | Tags |")
    println(io, "|---|---|---|---:|---:|---|")
    for row in rows
      println(io, "| ", _md(row.name), " | ", _md(row.file), " | ", _md(row.status), " | ", _md(row.exit_code), " | ", _md(row.duration_s), " | ", _md(row.tags), " |")
    end

    failures = [row for row in rows if _is_failure_status(row.status)]
    if !isempty(failures)
      println(io, "\n## Diagnostics\n")
      for row in failures
        println(io, "- **", row.name, " / ", row.status, "**: ", isempty(row.message) ? "no additional message" : row.message, " (log: `", row.log_file, "`)")
      end
    end
    for note in notes
      println(io, "\n", note)
    end
  end
  return path
end

"""
    run_example_suite(suite_name, title, specs; notes = String[], args = ARGS)

Runs every registered example of the suite in a fresh subprocess (subject to
the CLI filters) and writes summary CSV/Markdown plus per-example logs to the
output directory. Returns a NamedTuple
`(output_dir, summary_rows, written_files, successful)` or `nothing` when
`--help`/`--list` was requested.
"""
function run_example_suite(suite_name::AbstractString, title::AbstractString, specs::Vector{ExampleSpec}; notes::Vector{String} = String[], args = ARGS)
  print_example_banner("examples/run_$(suite_name).jl", title)
  opt = parse_example_suite_cli(suite_name, title, specs, notes, args)
  opt === nothing && return nothing
  if opt["list"]
    _print_suite_list(specs)
    return nothing
  end

  duplicate_names = length(Set(spec.name for spec in specs)) == length(specs)
  duplicate_names || throw(ArgumentError("duplicate example names in the $suite_name registry"))
  only_names = isempty(opt["only"]) ? nothing : _split_name_list(opt["only"], specs, "only")
  skip_names = isempty(opt["skip"]) ? Set{String}() : _split_name_list(opt["skip"], specs, "skip")

  output_dir = normpath(opt["output-dir"])
  opt["output-dir"] = output_dir
  logs_dir = joinpath(output_dir, "logs")
  mkpath(logs_dir)
  examples_dir = normpath(joinpath(@__DIR__, ".."))

  rows = NamedTuple[]
  name_width = maximum(length(spec.name) for spec in specs)
  for (index, spec) in enumerate(specs)
    opt["quiet"] || print("[", lpad(index, 2), "/", length(specs), "] ", rpad(spec.name, name_width), " ... ")

    skip = _spec_skip_status(spec, opt, only_names, skip_names)
    if skip !== nothing
      status, message = skip
      push!(rows, _summary_row(spec; status = status, message = message))
      opt["quiet"] || println("SKIPPED (", replace(status, "skipped_" => ""), "): ", message)
      continue
    end

    row = _run_spec_subprocess(spec, opt, examples_dir, logs_dir)
    push!(rows, row)
    if !opt["quiet"]
      println(uppercase(row.status), " (", row.duration_s, " s)")
      _is_failure_status(row.status) && println(stderr, "  ", row.message, " — ", row.log_file)
    end
    if _is_failure_status(row.status) && !opt["continue-on-error"]
      break
    end
  end

  written_files = String[]
  opt["write-csv"] && push!(written_files, _write_csv(joinpath(output_dir, "$(suite_name)_summary.csv"), rows))
  opt["write-markdown"] && push!(written_files, _write_suite_markdown(joinpath(output_dir, "$(suite_name)_summary.md"), title, opt, rows, notes))

  failed_rows = [row for row in rows if _is_failure_status(row.status)]
  skipped_rows = [row for row in rows if row.status in _SUITE_SKIP_STATUSES]
  if !opt["quiet"]
    println("\n", title)
    println("Output directory: ", output_dir)
    println("Examples: ", length(rows), " (ok: ", length(rows) - length(failed_rows) - length(skipped_rows), ", failed: ", length(failed_rows), ", skipped: ", length(skipped_rows), ")")
    for path in written_files
      println("  - ", path)
    end
    for note in notes
      println(note)
    end
  end

  result = (output_dir = output_dir, summary_rows = rows, written_files = written_files, successful = isempty(failed_rows))

  if opt["strict"] && !isempty(failed_rows)
    details = join(("$(row.name): $(row.status)" for row in failed_rows), "; ")
    error("$title strict checks failed: $details")
  end

  return result
end

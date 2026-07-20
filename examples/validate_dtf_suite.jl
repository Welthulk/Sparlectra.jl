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

# Date: 2026-07-16
# file: examples/validate_dtf_suite.jl
# purpose: unified CLI runner for DTF import audit, native DTF/FOR002 base/outage validation, and DTF->MATPOWER->Sparlectra roundtrip validation

# Unified DTF validation suite.
#
# Thin CLI runner over the shared modules in examples/internal/
# (dtf_validation_base/outages/matpower/audit plus
# dtf_for002_validation_utils.jl). It combines:
#   - DTF import audit
#   - native DTF/FOR002 base-case validation
#   - native DTF/FOR002 outage validation
#   - DTF -> MATPOWER -> Sparlectra roundtrip validation
#
# Default input:  <repository>/data/DTF/FOR001*.DAT and FOR002*.DAT
# Default output: examples/_out/dtf_validation_suite/

using Dates
using Printf
using Sparlectra



include(joinpath(@__DIR__, "internal", "dtf_validation_base.jl"))
include(joinpath(@__DIR__, "internal", "dtf_validation_outages.jl"))
include(joinpath(@__DIR__, "internal", "dtf_validation_matpower.jl"))
include(joinpath(@__DIR__, "internal", "dtf_validation_audit.jl"))
include(joinpath(@__DIR__, "internal", "example_header.jl"))

const _SUITE_MODES = ("audit", "base", "outages", "matpower")

function _parse_bool(value::AbstractString)
  normalized = lowercase(strip(value))
  normalized in ("true", "1", "yes", "on") && return true
  normalized in ("false", "0", "no", "off") && return false
  throw(ArgumentError("invalid boolean value: $value"))
end

function _default_suite_options()
  repo = normpath(joinpath(@__DIR__, ".."))
  return Dict{String,Any}(
    "mode" => "all",
    "case" => "all",
    "data-dir" => joinpath(repo, "data", "DTF"),
    "output-dir" => joinpath(@__DIR__, "_out", "dtf_validation_suite"),
    "tol" => 1e-8,
    "max-iter" => 50,
    "method" => "rectangular",
    "write-csv" => true,
    "write-markdown" => true,
    "write-matpower" => true,
    "run-outages" => true,
    "quiet" => false,
    "strict" => false,
    "parser-strict" => false,
    "continue-on-error" => true,
    "legacy-voltage-level-collapse-230kv" => false,
    "transformer-ratio-mode" => "neutral_one",
    "tap-changer-model" => "impedance_correction",
  )
end

function _print_help()
  println("""
Unified DTF validation suite

Usage:
  julia --project=. examples/validate_dtf_suite.jl [options]

Main options:
  --mode=all|audit|base|outages|matpower
      Multiple modes may be comma-separated, for example --mode=audit,base.
  --case=all|A|B|C|D|E
      Multiple cases may be comma-separated. Case A maps to FOR001.DAT/FOR002.DAT.
  --data-dir=<path>
  --output-dir=<path>

Solver options:
  --tol=1e-8
  --max-iter=50
  --method=rectangular|polar|default
  --transformer-ratio-mode=neutral_one|winding_over_network
  --tap-changer-model=ideal|impedance_correction
      ideal (default) keeps prior Sparlectra behavior where the tap changer
      does not affect R/X. impedance_correction re-refers R/X through the
      tapped winding (|1 + f*e^(j*phi)|^2); this removes the residual
      deviations seen in tap-stepped cases (for example C/E) against FOR002.
  --legacy-voltage-level-collapse-230kv=true|false

Output and control:
  --write-csv=true|false
  --write-markdown=true|false
  --write-matpower=true|false
  --run-outages=true|false
  --parser-strict=true|false
  --strict=true|false
  --continue-on-error=true|false
  --quiet
  --help

File pairing:
  FOR001.DAT  <-> FOR002.DAT   (case A)
  FOR001B.DAT <-> FOR002B.DAT  (case B)
  Additional matching suffixes are detected automatically.
""")
end

function parse_suite_cli(args = ARGS)
  opt = _default_suite_options()
  boolean_keys = Set([
    "write-csv",
    "write-markdown",
    "write-matpower",
    "run-outages",
    "quiet",
    "strict",
    "parser-strict",
    "continue-on-error",
    "legacy-voltage-level-collapse-230kv",
  ])
  for arg in args
    if arg in ("--help", "-h")
      _print_help()
      return nothing
    elseif arg == "--quiet"
      opt["quiet"] = true
      continue
    elseif arg == "--strict"
      opt["strict"] = true
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

    if key == "tol"
      opt[key] = parse(Float64, value)
    elseif key == "max-iter"
      opt[key] = parse(Int, value)
    elseif key in boolean_keys
      opt[key] = _parse_bool(value)
    else
      opt[key] = value
    end
  end

  opt["data-dir"] = normpath(abspath(opt["data-dir"]))
  opt["output-dir"] = normpath(abspath(opt["output-dir"]))
  opt["tol"] > 0 || throw(ArgumentError("--tol must be positive"))
  opt["max-iter"] > 0 || throw(ArgumentError("--max-iter must be positive"))
  return opt
end

function _selected_modes(value::AbstractString)
  raw = lowercase.(strip.(split(value, ",")))
  "all" in raw && return collect(_SUITE_MODES)
  aliases = Dict("native" => "base", "outage" => "outages", "roundtrip" => "matpower")
  modes = String[]
  for token in raw
    isempty(token) && continue
    mode = get(aliases, token, token)
    mode in _SUITE_MODES || throw(ArgumentError("unknown mode '$token'; allowed: all, " * join(_SUITE_MODES, ", ")))
    mode in modes || push!(modes, mode)
  end
  isempty(modes) && throw(ArgumentError("no validation mode selected"))
  return modes
end

function _selected_cases(value::AbstractString)
  tokens = uppercase.(strip.(split(value, ",")))
  "ALL" in tokens && return nothing
  selected = Set(filter(token -> !isempty(token), tokens))
  isempty(selected) && throw(ArgumentError("no case selected"))
  return selected
end

function _discover_for_files(data_dir::AbstractString)
  isdir(data_dir) || throw(ArgumentError("DTF data directory does not exist: $data_dir"))
  for001 = Dict{String,String}()
  for002 = Dict{String,String}()

  for name in readdir(data_dir)
    path = joinpath(data_dir, name)
    isfile(path) || continue

    m1 = match(r"^FOR001(.*)\.DAT$"i, name)
    if m1 !== nothing
      suffix = uppercase(String(something(m1.captures[1], "")))
      haskey(for001, suffix) && throw(ArgumentError("duplicate FOR001 suffix '$suffix': $(for001[suffix]) and $path"))
      for001[suffix] = path
      continue
    end

    m2 = match(r"^FOR002(.*)\.DAT$"i, name)
    if m2 !== nothing
      suffix = uppercase(String(something(m2.captures[1], "")))
      haskey(for002, suffix) && throw(ArgumentError("duplicate FOR002 suffix '$suffix': $(for002[suffix]) and $path"))
      for002[suffix] = path
    end
  end

  suffixes = sort!(collect(union(keys(for001), keys(for002))); by = s -> (isempty(s) ? 0 : 1, s))
  return [(suffix = s, case_id = isempty(s) ? "A" : s, dtf_file = get(for001, s, missing), for002_file = get(for002, s, missing)) for s in suffixes]
end

_safe_case_dir(case_id::AbstractString) = "case_" * replace(case_id, r"[^A-Za-z0-9_.-]+" => "_")

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

function _finite_max(values)
  usable = Float64[]
  for value in values
    (value isa Missing || value === nothing) && continue
    number = try
      Float64(value)
    catch
      continue
    end
    isfinite(number) && push!(usable, abs(number))
  end
  return isempty(usable) ? missing : maximum(usable)
end

# Keep in sync with the identical helper in examples/internal/dtf_for002_validation_utils.jl.
function _fmt_num(x::AbstractFloat)::String
  isfinite(x) || return string(x)
  x == 0.0 && return "0.0"
  ax = abs(x)
  (ax < 1.0e-4 || ax >= 1.0e6) && return @sprintf("%.3e", x)
  ax < 0.01 && return string(round(x; digits = 1 - floor(Int, log10(ax))))
  return string(round(x; digits = 3))
end
_fmt_num(x::Real) = _fmt_num(Float64(x))

function _md(value)
  value isa Missing && return ""
  value === nothing && return ""
  value isa AbstractFloat && return _fmt_num(value)
  return replace(string(value), "|" => "\\|", "\n" => " ", "\r" => " ")
end

function _summary_row(;
  case_id,
  suffix,
  mode,
  status,
  dtf_file = missing,
  for002_file = missing,
  output_dir = missing,
  converged = missing,
  scenario_count = missing,
  max_abs_d_vm_kV = missing,
  max_abs_d_va_deg = missing,
  max_abs_d_p_MW = missing,
  max_abs_d_q_MVar = missing,
  message = "",
)
  return (
    case_id = case_id,
    suffix = suffix,
    mode = mode,
    status = status,
    converged = converged,
    scenario_count = scenario_count,
    max_abs_d_vm_kV = max_abs_d_vm_kV,
    max_abs_d_va_deg = max_abs_d_va_deg,
    max_abs_d_p_MW = max_abs_d_p_MW,
    max_abs_d_q_MVar = max_abs_d_q_MVar,
    dtf_file = dtf_file,
    for002_file = for002_file,
    output_dir = output_dir,
    message = message,
  )
end

function _common_validation_args(opt, dtf_file, for002_file, output_dir)
  return String[
    "--dtf-file=$(dtf_file)",
    "--for002-file=$(for002_file)",
    "--output-dir=$(output_dir)",
    "--tol=$(opt["tol"])",
    "--max-iter=$(opt["max-iter"])",
    "--method=$(opt["method"])",
    "--write-csv=$(opt["write-csv"])",
    "--write-markdown=$(opt["write-markdown"])",
    "--quiet=true",
    "--print-summary=false",
    "--tap-changer-model=$(opt["tap-changer-model"])",
  ]
end

function _run_import_audit(case_id, dtf_file, output_dir, opt)
  mkpath(output_dir)
  path = joinpath(output_dir, "dtf_import_audit.md")
  open(path, "w") do io
    println(io, "# DTF import audit\n")
    println(io, "- generated: ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
    println(io, "- parser strict mode: ", opt["parser-strict"], "\n")
    ImportAudit.audit_case(io, case_id, dtf_file; strict = opt["parser-strict"])
  end
  return path
end

function _run_base_validation(dtf_file, for002_file, output_dir, opt)
  args = _common_validation_args(opt, dtf_file, for002_file, output_dir)
  push!(args, "--strict=$(opt["parser-strict"])")
  push!(args, "--legacy-voltage-level-collapse-230kv=$(opt["legacy-voltage-level-collapse-230kv"])")
  push!(args, "--transformer-ratio-mode=$(opt["transformer-ratio-mode"])")
  return NativeBaseValidation.run_validation(args; return_details = false)
end

function _run_outage_validation(dtf_file, for002_file, output_dir, opt)
  args = _common_validation_args(opt, dtf_file, for002_file, output_dir)
  push!(args, "--strict=false")
  return NativeOutageValidation.run_validation(args; return_details = false)
end

function _run_matpower_validation(dtf_file, for002_file, output_dir, opt)
  opt["write-matpower"] || throw(ArgumentError("MATPOWER roundtrip validation requires --write-matpower=true"))
  args = _common_validation_args(opt, dtf_file, for002_file, output_dir)
  append!(args, [
    "--write-matpower=true",
    "--run-outages=$(opt["run-outages"])",
    "--strict=false",
    "--details=true",
  ])
  return MatpowerRoundtripValidation.run_validation(args; return_details = true)
end

function _write_suite_markdown(path, opt, inventory, rows, modes)
  open(path, "w") do io
    println(io, "# Unified DTF validation suite\n")
    println(io, "- generated: ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
    println(io, "- data directory: `", opt["data-dir"], "`")
    println(io, "- selected modes: `", join(modes, ", "), "`")
    println(io, "- selected cases: `", opt["case"], "`")
    println(io, "- method: `", opt["method"], "`; tolerance: `", opt["tol"], "`; max iterations: `", opt["max-iter"], "`")
    println(io, "- transformer ratio mode: `", opt["transformer-ratio-mode"], "`")
    println(io, "- tap-changer model: `", opt["tap-changer-model"], "`")
    println(io, "- parser strict mode: `", opt["parser-strict"], "`\n")

    println(io, "## Input inventory\n")
    println(io, "| Case | Suffix | FOR001 | FOR002 |")
    println(io, "|---|---:|---|---|")
    for item in inventory
      println(io, "| ", _md(item.case_id), " | ", _md(item.suffix), " | ", _md(item.dtf_file), " | ", _md(item.for002_file), " |")
    end

    println(io, "\n## Results\n")
    println(io, "| Case | Mode | Status | Converged | Scenarios | max abs(dV) kV | max abs(dVa) deg | max abs(dP) MW | max abs(dQ) MVar |")
    println(io, "|---|---|---|---:|---:|---:|---:|---:|---:|")
    for row in rows
      println(
        io,
        "| ", _md(row.case_id),
        " | ", _md(row.mode),
        " | ", _md(row.status),
        " | ", _md(row.converged),
        " | ", _md(row.scenario_count),
        " | ", _md(row.max_abs_d_vm_kV),
        " | ", _md(row.max_abs_d_va_deg),
        " | ", _md(row.max_abs_d_p_MW),
        " | ", _md(row.max_abs_d_q_MVar),
        " |",
      )
    end

    failures = [row for row in rows if !(row.status in ("ok", "ok_no_outages"))]
    if !isempty(failures)
      println(io, "\n## Diagnostics\n")
      for row in failures
        println(io, "- **", row.case_id, " / ", row.mode, " / ", row.status, "**: ", isempty(row.message) ? "no additional message" : row.message)
      end
    end
  end
  return path
end

function _handle_error!(rows, item, mode, mode_output, err, opt)
  message = sprint(showerror, err)
  push!(
    rows,
    _summary_row(
      case_id = item.case_id,
      suffix = item.suffix,
      mode = mode,
      status = "error",
      dtf_file = item.dtf_file,
      for002_file = item.for002_file,
      output_dir = mode_output,
      message = message,
    ),
  )
  opt["quiet"] || println(stderr, "[", item.case_id, "/", mode, "] ERROR: ", message)
  opt["continue-on-error"] || throw(err)
  return nothing
end

function run_suite(args = ARGS)
  print_example_banner("examples/validate_dtf_suite.jl", "unified CLI runner for DTF import audit, native DTF/FOR002 base/outage validation, and DTF->MATPOWER->Sparlectra roundtrip validation")
  opt = parse_suite_cli(args)
  opt === nothing && return nothing

  modes = _selected_modes(opt["mode"])
  selected_cases = _selected_cases(opt["case"])
  inventory = _discover_for_files(opt["data-dir"])
  isempty(inventory) && throw(ArgumentError("no FOR001*.DAT or FOR002*.DAT files found in $(opt["data-dir"])"))

  if selected_cases !== nothing
    inventory = [item for item in inventory if uppercase(item.case_id) in selected_cases]
    isempty(inventory) && throw(ArgumentError("none of the requested cases were found: $(join(sort!(collect(selected_cases)), ", "))"))
  end

  mkpath(opt["output-dir"])
  rows = NamedTuple[]

  for item in inventory
    case_output = joinpath(opt["output-dir"], _safe_case_dir(item.case_id))
    opt["quiet"] || println("\nCase ", item.case_id, ": ", item.dtf_file, " / ", item.for002_file)

    if item.dtf_file isa Missing
      for mode in modes
        push!(
          rows,
          _summary_row(
            case_id = item.case_id,
            suffix = item.suffix,
            mode = mode,
            status = "skipped_missing_for001",
            for002_file = item.for002_file,
            output_dir = joinpath(case_output, mode),
            message = "No matching FOR001 file was found.",
          ),
        )
      end
      continue
    end

    for mode in modes
      mode_output = joinpath(case_output, mode)
      opt["quiet"] || print("  ", rpad(mode, 10), " ... ")

      if mode in ("base", "outages") && item.for002_file isa Missing
        push!(
          rows,
          _summary_row(
            case_id = item.case_id,
            suffix = item.suffix,
            mode = mode,
            status = "skipped_missing_for002",
            dtf_file = item.dtf_file,
            output_dir = mode_output,
            message = "No matching FOR002 reference file was found.",
          ),
        )
        opt["quiet"] || println("SKIPPED (missing FOR002)")
        continue
      end

      try
        if mode == "audit"
          path = _run_import_audit(item.case_id, item.dtf_file, mode_output, opt)
          push!(
            rows,
            _summary_row(
              case_id = item.case_id,
              suffix = item.suffix,
              mode = mode,
              status = "ok",
              dtf_file = item.dtf_file,
              for002_file = item.for002_file,
              output_dir = mode_output,
              scenario_count = 1,
              message = path,
            ),
          )

        elseif mode == "base"
          result = _run_base_validation(item.dtf_file, item.for002_file, mode_output, opt)
          metrics = result.metrics
          status = result.converged ? "ok" : "not_converged"
          push!(
            rows,
            _summary_row(
              case_id = item.case_id,
              suffix = item.suffix,
              mode = mode,
              status = status,
              dtf_file = item.dtf_file,
              for002_file = item.for002_file,
              output_dir = mode_output,
              converged = result.converged,
              scenario_count = 1,
              max_abs_d_vm_kV = metrics.max_abs_d_vm_kV,
              max_abs_d_va_deg = metrics.max_abs_d_va_deg,
              max_abs_d_p_MW = result.max_branch_d_p_MW,
              max_abs_d_q_MVar = result.max_branch_d_q_MVar,
              message = "iterations=$(result.iterations), final_mismatch=$(result.final_mismatch)",
            ),
          )

        elseif mode == "outages"
          result = _run_outage_validation(item.dtf_file, item.for002_file, mode_output, opt)
          metrics = result.metrics_rows
          all_matched = all(row -> row.match_status == "matched", result.matching_rows)
          all_converged = all(row -> row.converged, metrics)
          status = isempty(metrics) ? "ok_no_outages" : (all_matched && all_converged ? "ok" : (!all_matched ? "outage_match_failed" : "not_converged"))
          push!(
            rows,
            _summary_row(
              case_id = item.case_id,
              suffix = item.suffix,
              mode = mode,
              status = status,
              dtf_file = item.dtf_file,
              for002_file = item.for002_file,
              output_dir = mode_output,
              converged = isempty(metrics) ? missing : all_converged,
              scenario_count = length(metrics),
              max_abs_d_vm_kV = _finite_max([m.max_abs_d_vm_kV for m in metrics]),
              max_abs_d_va_deg = _finite_max([m.max_abs_d_va_deg for m in metrics]),
              max_abs_d_p_MW = _finite_max([m.max_abs_branch_d_p_MW for m in metrics]),
              max_abs_d_q_MVar = _finite_max([m.max_abs_branch_d_q_MVar for m in metrics]),
              message = "parsed_outages=$(result.parsed_outages), parsed_for002_blocks=$(result.parsed_for002_outage_blocks)",
            ),
          )

        elseif mode == "matpower"
          for002_arg = item.for002_file isa Missing ? joinpath(opt["data-dir"], "FOR002.DAT") : item.for002_file
          result = _run_matpower_validation(item.dtf_file, for002_arg, mode_output, opt)
          metrics = result.metrics_rows
          all_converged = all(row -> row.native_converged && row.roundtrip_converged, metrics)
          all_status = all(row -> row.status_match, metrics)
          status = all_converged && all_status ? "ok" : (!all_converged ? "not_converged" : "status_mismatch")
          push!(
            rows,
            _summary_row(
              case_id = item.case_id,
              suffix = item.suffix,
              mode = mode,
              status = status,
              dtf_file = item.dtf_file,
              for002_file = item.for002_file,
              output_dir = mode_output,
              converged = all_converged,
              scenario_count = length(metrics),
              max_abs_d_vm_kV = _finite_max([m.max_abs_d_vm_kV for m in metrics]),
              max_abs_d_va_deg = _finite_max([m.max_abs_d_va_deg for m in metrics]),
              max_abs_d_p_MW = _finite_max([m.max_abs_branch_d_p_MW for m in metrics]),
              max_abs_d_q_MVar = _finite_max([m.max_abs_branch_d_q_MVar for m in metrics]),
              message = "roundtrip branch status match=$(all_status)",
            ),
          )
        end
        opt["quiet"] || println(last(rows).status)
      catch err
        _handle_error!(rows, item, mode, mode_output, err, opt)
      end
    end
  end

  written_files = String[]
  if opt["write-csv"]
    push!(written_files, _write_csv(joinpath(opt["output-dir"], "dtf_validation_suite_summary.csv"), rows))
    inventory_rows = [
      (
        case_id = item.case_id,
        suffix = item.suffix,
        dtf_file = item.dtf_file,
        for002_file = item.for002_file,
        has_for001 = !(item.dtf_file isa Missing),
        has_for002 = !(item.for002_file isa Missing),
      ) for item in inventory
    ]
    push!(written_files, _write_csv(joinpath(opt["output-dir"], "dtf_validation_suite_inventory.csv"), inventory_rows))
  end
  if opt["write-markdown"]
    push!(
      written_files,
      _write_suite_markdown(
        joinpath(opt["output-dir"], "dtf_validation_suite_summary.md"),
        opt,
        inventory,
        rows,
        modes,
      ),
    )
  end

  failed_rows = [row for row in rows if !(row.status in ("ok", "ok_no_outages"))]
  if !opt["quiet"]
    println("\nUnified DTF validation suite")
    println("Output directory: ", opt["output-dir"])
    println("Checks: ", length(rows))
    println("Successful: ", length(rows) - length(failed_rows))
    println("Non-successful: ", length(failed_rows))
    for path in written_files
      println("  - ", path)
    end
  end

  result = (
    output_dir = opt["output-dir"],
    inventory = inventory,
    summary_rows = rows,
    written_files = written_files,
    successful = isempty(failed_rows),
  )

  if opt["strict"] && !isempty(failed_rows)
    details = join(("$(row.case_id)/$(row.mode): $(row.status)" for row in failed_rows), "; ")
    error("DTF validation suite strict checks failed: $details")
  end

  return result
end

run_example(run_suite, ARGS)
nothing

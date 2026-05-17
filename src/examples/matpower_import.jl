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
import Sparlectra: MatpowerIO, _perf_profile_time!
using BenchmarkTools
using Printf
using Dates
using Logging

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

function _strip_yaml_inline_comment(raw::AbstractString)
  buf = IOBuffer()
  qchar = '\0'
  for ch in raw
    if qchar != '\0'
      ch == qchar && (qchar = '\0')
      print(buf, ch)
    elseif ch == '\'' || ch == '"'
      qchar = ch
      print(buf, ch)
    elseif ch == '#'
      break
    else
      print(buf, ch)
    end
  end
  return String(take!(buf))
end

function load_yaml_config(path::AbstractString)
  isempty(path) && return Dict{String,Any}()
  isfile(path) || error("YAML config file not found: $path")

  cfg = Dict{String,Any}()
  current_section = nothing
  pending_list_key = nothing
  pending_list_section = nothing
  for rawline in eachline(path)
    uncommented = rstrip(_strip_yaml_inline_comment(rawline))
    stripped = strip(uncommented)
    isempty(stripped) && continue
    startswith(stripped, "#") && continue

    indent = firstindex(uncommented)
    while indent <= lastindex(uncommented) && uncommented[indent] == ' '
      indent = nextind(uncommented, indent)
    end
    is_nested = indent > firstindex(uncommented)

    if startswith(stripped, "-") && !isnothing(pending_list_key)
      item_raw = strip(stripped[2:end])
      target = isnothing(pending_list_section) ? cfg : cfg[pending_list_section]
      target[pending_list_key] isa AbstractVector || (target[pending_list_key] = Any[])
      push!(target[pending_list_key], _parse_yaml_scalar(item_raw))
      continue
    end

    occursin(":", stripped) || continue
    key, value_raw = split(stripped, ":"; limit = 2)
    key = strip(key)
    value_raw = strip(value_raw)
    target = cfg
    if is_nested && !isnothing(current_section) && cfg[current_section] isa Dict{String,Any}
      target = cfg[current_section]
    elseif !is_nested
      current_section = nothing
    end

    if isempty(value_raw)
      if is_nested
        target[key] = Any[]
        pending_list_key = key
        pending_list_section = current_section
      else
        cfg[key] = Dict{String,Any}()
        current_section = key
        pending_list_key = key
        pending_list_section = nothing
      end
    elseif startswith(value_raw, "[") && endswith(value_raw, "]")
      target[key] = _parse_yaml_list(value_raw)
      pending_list_key = nothing
      pending_list_section = nothing
    else
      target[key] = _parse_yaml_scalar(value_raw)
      pending_list_key = nothing
      pending_list_section = nothing
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

function _as_output_mode(v)::Symbol
  s = lowercase(String(v))
  if s == "classic"
    return :classic
  elseif s == "dataframe"
    return :dataframe
  else
    @warn "Unknown show_once_output; using :classic" value = v
    return :classic
  end
end

function _as_performance_level(v)::Symbol
  s = lowercase(String(v))
  s in ("false", "off", "none", "no", "0") && return :off
  s in ("summary", "detailed", "iteration") && return Symbol(s)
  @warn "Unknown performance.level; using :summary" value = v
  return :summary
end

function _yaml_section_get(cfg::Dict{String,Any}, section::AbstractString, key::AbstractString, default)
  flat_key = string(section, "_", key)
  haskey(cfg, flat_key) && return cfg[flat_key]
  sec = get(cfg, section, nothing)
  sec isa Dict{String,Any} && haskey(sec, key) && return sec[key]
  return default
end

function _as_console_mode(v)::Symbol
  v === false && return :off
  v === true && return :compact
  s = lowercase(String(v))
  s in ("false", "off", "none", "no", "0") && return :off
  s in ("compact", "summary", "full") && return Symbol(s)
  @warn "Unknown console/logfile verbosity mode; using :compact" value = v
  return :compact
end

function _yaml_path_from_inputs()
  !isempty(ARGS) && return ARGS[1]
  env_path = get(ENV, "SPARLECTRA_MATPOWER_IMPORT_YAML", "")
  !isempty(env_path) && return env_path

  local_default = joinpath(@__DIR__, "matpower_import.yaml")
  isfile(local_default) && return local_default
  local_example = joinpath(@__DIR__, "matpower_import.yaml.example")
  isfile(local_example) && return local_example
  return ""
end
# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
const DEFAULT_CASE = "case141.m"
const DEFAULT_METHODS = [:rectangular]
const METHODS = [:rectangular]

# -----------------------------------------------------------------------------
# Output redirection (write verbose output to git-ignored file)
# -----------------------------------------------------------------------------
const OUTDIR = joinpath(@__DIR__, "_out")
mkpath(OUTDIR)

mutable struct _MatpowerImportLogStatus
  warnings::Int
  errors::Int
end

_MatpowerImportLogStatus() = _MatpowerImportLogStatus(0, 0)

mutable struct _MatpowerImportTableLogger <: Logging.AbstractLogger
  io::IO
  min_level::Logging.LogLevel
  status::_MatpowerImportLogStatus
  header_printed::Bool
end

Logging.min_enabled_level(logger::_MatpowerImportTableLogger) = logger.min_level
Logging.catch_exceptions(logger::_MatpowerImportTableLogger) = true
Logging.shouldlog(logger::_MatpowerImportTableLogger, level, _module, group, id) = level >= logger.min_level

function _log_level_label(level)::String
  level >= Logging.Error && return "ERROR"
  level >= Logging.Warn && return "WARN"
  return string(level)
end

function _log_kv_string(kwargs)::String
  isempty(kwargs) && return ""
  parts = String[]
  for (key, value) in kwargs
    push!(parts, string(key, "=", repr(value)))
  end
  return join(parts, "; ")
end

function _log_source(_module, file, line)::String
  file_label = isnothing(file) ? "unknown" : basename(String(file))
  return string(_module, ":", file_label, ":", line)
end

function Logging.handle_message(logger::_MatpowerImportTableLogger, level, message, _module, group, id, file, line; kwargs...)
  level >= Logging.Error && (logger.status.errors += 1)
  Logging.Warn <= level < Logging.Error && (logger.status.warnings += 1)
  if !logger.header_printed
    println(logger.io, "")
    println(logger.io, "==================== Reported warnings/errors ====================")
    @printf(logger.io, "%-23s  %-5s  %-28s  %-24s  %s\n", "timestamp", "level", "source", "message", "details")
    @printf(logger.io, "%-23s  %-5s  %-28s  %-24s  %s\n", repeat("-", 23), repeat("-", 5), repeat("-", 28), repeat("-", 24), repeat("-", 40))
    logger.header_printed = true
  end
  source = _log_source(_module, file, line)
  @printf(logger.io, "%-23s  %-5s  %-28s  %-24s  %s\n", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"), _log_level_label(level), source, string(message), _log_kv_string(kwargs))
  flush(logger.io)
  return nothing
end

function _with_log_table(f::Function, logfile::AbstractString, status::_MatpowerImportLogStatus)
  open(logfile, "a") do io
    logger = _MatpowerImportTableLogger(io, Logging.Warn, status, false)
    return Logging.with_logger(logger) do
      f()
    end
  end
end

function _print_log_status(io::IO, status::_MatpowerImportLogStatus)
  println(io, "reported warnings/errors: warnings=", status.warnings, " errors=", status.errors)
  return nothing
end

function _new_performance_profile(cfg)
  return Dict{Symbol,Any}(
    :enabled => cfg.performance_enabled && cfg.performance_level != :off,
    :level => cfg.performance_level,
    :show_allocations => cfg.performance_show_allocations,
    :show_iteration_table => cfg.performance_show_iteration_table,
  )
end

function _print_performance_profile(io::IO, profile; title::AbstractString = "Performance Summary", max_rows::Int = 20)
  profile isa AbstractDict || return nothing
  Bool(get(profile, :enabled, false)) || return nothing
  timings = get(profile, :timings, Dict{Symbol,Any}())
  println(io, "\n==================== ", title, " ====================")
  if isempty(timings)
    println(io, "No timing samples recorded.")
  else
    rows = sort(collect(timings); by = pair -> pair.second.elapsed_s, rev = true)
    nshow = min(length(rows), max(0, max_rows))
    show_alloc = Bool(get(profile, :show_allocations, false))
    if show_alloc
      @printf(io, "%-40s %8s %12s %14s\n", "phase", "calls", "seconds", "allocated")
    else
      @printf(io, "%-40s %8s %12s\n", "phase", "calls", "seconds")
    end
    for (phase, row) in rows[1:nshow]
      if show_alloc
        @printf(io, "%-40s %8d %12.6f %14d\n", String(phase), row.calls, row.elapsed_s, row.bytes)
      else
        @printf(io, "%-40s %8d %12.6f\n", String(phase), row.calls, row.elapsed_s)
      end
    end
    if nshow < length(rows)
      @printf(io, "... truncated performance rows: showing %d/%d\n", nshow, length(rows))
    end
    sp = get(profile, :start_projection_summary, nothing)
    if sp !== nothing
      @printf(io, "start_projection: selected=%s, candidates=%d, best_mismatch=%.6e, time=%.6f s\n", String(sp.selected), sp.candidates, sp.best_mismatch, sp.elapsed_s)
    end
    if haskey(timings, :solver_total)
      println(io, "Note: solver_total is inclusive; nested phase rows explain visible subwork, and any remainder is solver control flow, status bookkeeping, or other uninstrumented overhead.")
    end
  end
  iterations = get(profile, :iterations, NamedTuple[])
  if Bool(get(profile, :show_iteration_table, true)) && get(profile, :level, :summary) === :iteration && !isempty(iterations)
    println(io, "\nNewton iteration table:")
    @printf(io, "%8s %16s %16s %16s\n", "iter", "max_mismatch", "qlimit_changed", "qlimit_reenabled")
    nshow = min(length(iterations), max(0, max_rows))
    for row in iterations[1:nshow]
      @printf(io, "%8d %16.6e %16s %16s\n", row.iteration, row.max_mismatch, string(row.qlimit_changed), string(row.qlimit_reenabled))
    end
    if nshow < length(iterations)
      @printf(io, "... truncated iteration rows: showing %d/%d\n", nshow, length(iterations))
    end
  end
  println(io, "===============================================================")
  return nothing
end

function _emit_performance_summary(profile; logfile::AbstractString = "", print_to_console::Bool = true, write_to_logfile::Bool = true, max_rows::Int = 20)
  profile isa AbstractDict || return nothing
  Bool(get(profile, :enabled, false)) || return nothing
  if print_to_console
    _print_performance_profile(stdout, profile; title = "Performance Summary", max_rows = max_rows)
  end
  if write_to_logfile && !isempty(logfile)
    open(logfile, "a") do io
      _print_performance_profile(io, profile; title = "Performance Summary", max_rows = max_rows)
    end
  end
  return nothing
end

# -----------------------------------------------------------------------------
# Compare against MATPOWER reference (if present)
# -----------------------------------------------------------------------------
const DEFAULT_SHOW_DIFF = true
const DEFAULT_TOL_VM = 2e-2
const DEFAULT_TOL_VA = 5e-1

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
    autodamp = false,
    autodamp_min = 1e-3,
    start_projection = false,
    start_projection_try_dc_start = true,
    start_projection_try_blend_scan = true,
    start_projection_branch_guard = true,
    start_projection_measure_candidates = true,
    start_projection_reuse_import_data = true,
    start_projection_blend_lambdas = [0.25, 0.5, 0.75],
    start_projection_dc_angle_limit_deg = 60.0,
    qlimit_start_iter = 2,
    qlimit_start_mode = :iteration,
    qlimit_auto_q_delta_pu = 1e-4,
    qlimit_guard = false,
    qlimit_guard_min_q_range_pu = 1e-4,
    qlimit_guard_zero_range_mode = :lock_pq,
    qlimit_guard_narrow_range_mode = :prefer_pq,
    qlimit_guard_max_switches = 10,
    qlimit_guard_freeze_after_repeated_switching = true,
    qlimit_guard_accept_bounded_violations = false,
    qlimit_guard_max_remaining_violations = 0,
    qlimit_guard_violation_mode = :delayed_switch,
    qlimit_guard_violation_threshold_pu = 1e-4,
    qlimit_guard_log = true,
    verbose = 1,
    cooldown_iters = 0,
    q_hyst_pu = 0.0,
    qlimit_trace_buses = Int[],
    lock_pv_to_pq_buses = Int[],
    ignore_q_limits = false,
    max_ite = 30,
    tol = 1e-6,
    show_diff = DEFAULT_SHOW_DIFF,
    tol_vm = DEFAULT_TOL_VM,
    tol_va = DEFAULT_TOL_VA,
    seconds = 2.0,
    samples = 50,
    show_once = true,
    show_once_output = :classic,
    show_once_max_nodes = 0,
    benchmark = true,
    matpower_shift_sign = 1.0,
    matpower_shift_unit = "deg",
    matpower_ratio = "normal",
    reference_override = false,
    reference_vm_pu = 1.0,
    reference_va_deg = 0.0,
    diagnose_matpower_reference = false,
    diagnose_branch_shift_conventions = false,
    diagnose_branch_neighborhood = false,
    diagnose_residual_clusters = false,
    diagnose_residual_cluster_threshold_mw = 1.0,
    diagnose_residual_cluster_maxlines = 20,
    diagnose_branch_neighborhood_buses = Int[],
    diagnose_branch_neighborhood_depth = 1,
    diagnose_branch_neighborhood_maxlines = 200,
    diagnose_nodal_balance_breakdown = false,
    diagnose_nodal_balance_buses = Int[],
    diagnose_nodal_balance_maxlines = 200,
    diagnose_nodal_balance_include_branches = true,
    diagnose_nodal_balance_include_generators = true,
    diagnose_nodal_balance_include_shunts = true,
    diagnose_negative_branch_impedance = true,
    diagnose_negative_branch_impedance_maxlines = 100,
    diagnose_negative_branch_impedance_fail_on_negative_r = false,
    diagnose_negative_branch_impedance_fail_on_negative_x = false,
    diagnose_negative_branch_impedance_warn_threshold_abs_r = 0.0,
    diagnose_negative_branch_impedance_warn_threshold_abs_x = 0.0,
    diagnose_maxlines = 12,
    log_effective_config = false,
    enable_pq_gen_controllers = true,
    bus_shunt_model = "admittance",
    trace_legacy_bus_type_warnings = false,
    matpower_pv_voltage_source = :gen_vg,
    matpower_pv_voltage_mismatch_tol_pu = 1e-4,
    compare_voltage_reference = :bus_vm,
    diagnose_pv_voltage_references = false,
    diagnose_pv_voltage_maxlines = 30,
    flatstart_voltage_mode = :classic,
    flatstart_angle_mode = :classic,
    flatstart_branch_guard = true,
    flatstart_min_vm_pu = 0.5,
    flatstart_max_angle_step_deg = 20.0,
    flatstart_max_vm_step_pu = 0.1,
    wrong_branch_detection = true,
    wrong_branch_min_vm_pu = 0.4,
    wrong_branch_max_angle_spread_deg = 120.0,
    wrong_branch_rescue = false,
    wrong_branch_rescue_modes = [:dc, :bus_vm_va_blend, :matpower_va],
    matpower_auto_profile = :off,
    matpower_auto_profile_log = true,
    console_summary = true,
    console_auto_profile = :compact,
    console_diagnostics = :compact,
    console_q_limit_events = :summary,
    console_max_rows = 20,
    logfile_diagnostics = :full,
    performance_enabled = false,
    performance_level = :summary,
    performance_print_to_console = true,
    performance_write_to_logfile = true,
    performance_show_allocations = false,
    performance_show_iteration_table = true,
    performance_compact_logging = true,
    performance_skip_reference_comparison = false,
    performance_skip_expensive_diagnostics = false,
    performance_skip_branch_neighborhood_report = true,
    performance_max_diagnostic_rows = 20,
  )
  case_override = get(CASE_BENCH_OVERRIDES, String(case_name), (;))
  yaml_override = (;)
  if !isempty(yaml_cfg)
    yaml_override = (;
      opt_fd = Bool(get(yaml_cfg, "opt_fd", base.opt_fd)),
      opt_sparse = Bool(get(yaml_cfg, "opt_sparse", base.opt_sparse)),
      opt_flatstart = Bool(get(yaml_cfg, "opt_flatstart", base.opt_flatstart)),
      autodamp = Bool(get(yaml_cfg, "autodamp", base.autodamp)),
      autodamp_min = Float64(get(yaml_cfg, "autodamp_min", base.autodamp_min)),
      start_projection = Bool(get(yaml_cfg, "start_projection", base.start_projection)),
      start_projection_try_dc_start = Bool(get(yaml_cfg, "start_projection_try_dc_start", base.start_projection_try_dc_start)),
      start_projection_try_blend_scan = Bool(get(yaml_cfg, "start_projection_try_blend_scan", base.start_projection_try_blend_scan)),
      start_projection_branch_guard = Bool(get(yaml_cfg, "start_projection_branch_guard", base.start_projection_branch_guard)),
      start_projection_measure_candidates = Bool(get(yaml_cfg, "start_projection_measure_candidates", base.start_projection_measure_candidates)),
      start_projection_reuse_import_data = Bool(get(yaml_cfg, "start_projection_reuse_import_data", base.start_projection_reuse_import_data)),
      start_projection_blend_lambdas = Float64[x for x in get(yaml_cfg, "start_projection_blend_lambdas", base.start_projection_blend_lambdas)],
      start_projection_dc_angle_limit_deg = Float64(get(yaml_cfg, "start_projection_dc_angle_limit_deg", base.start_projection_dc_angle_limit_deg)),
      qlimit_start_iter = Int(get(yaml_cfg, "qlimit_start_iter", base.qlimit_start_iter)),
      qlimit_start_mode = Symbol(get(yaml_cfg, "qlimit_start_mode", base.qlimit_start_mode)),
      qlimit_auto_q_delta_pu = Float64(get(yaml_cfg, "qlimit_auto_q_delta_pu", base.qlimit_auto_q_delta_pu)),
      qlimit_guard = Bool(get(yaml_cfg, "qlimit_guard", base.qlimit_guard)),
      qlimit_guard_min_q_range_pu = Float64(get(yaml_cfg, "qlimit_guard_min_q_range_pu", base.qlimit_guard_min_q_range_pu)),
      qlimit_guard_zero_range_mode = Symbol(get(yaml_cfg, "qlimit_guard_zero_range_mode", base.qlimit_guard_zero_range_mode)),
      qlimit_guard_narrow_range_mode = Symbol(get(yaml_cfg, "qlimit_guard_narrow_range_mode", base.qlimit_guard_narrow_range_mode)),
      qlimit_guard_max_switches = Int(get(yaml_cfg, "qlimit_guard_max_switches", base.qlimit_guard_max_switches)),
      qlimit_guard_freeze_after_repeated_switching = Bool(get(yaml_cfg, "qlimit_guard_freeze_after_repeated_switching", base.qlimit_guard_freeze_after_repeated_switching)),
      qlimit_guard_accept_bounded_violations = Bool(get(yaml_cfg, "qlimit_guard_accept_bounded_violations", base.qlimit_guard_accept_bounded_violations)),
      qlimit_guard_max_remaining_violations = Int(get(yaml_cfg, "qlimit_guard_max_remaining_violations", base.qlimit_guard_max_remaining_violations)),
      qlimit_guard_violation_mode = Symbol(get(yaml_cfg, "qlimit_guard_violation_mode", base.qlimit_guard_violation_mode)),
      qlimit_guard_violation_threshold_pu = Float64(get(yaml_cfg, "qlimit_guard_violation_threshold_pu", base.qlimit_guard_violation_threshold_pu)),
      qlimit_guard_log = Bool(get(yaml_cfg, "qlimit_guard_log", base.qlimit_guard_log)),
      verbose = Int(get(yaml_cfg, "verbose", base.verbose)),
      cooldown_iters = Int(get(yaml_cfg, "cooldown_iters", base.cooldown_iters)),
      q_hyst_pu = Float64(get(yaml_cfg, "q_hyst_pu", base.q_hyst_pu)),
      qlimit_trace_buses = _as_int_vec(get(yaml_cfg, "qlimit_trace_buses", base.qlimit_trace_buses)),
      lock_pv_to_pq_buses = _as_int_vec(get(yaml_cfg, "lock_pv_to_pq_buses", base.lock_pv_to_pq_buses)),
      ignore_q_limits = Bool(get(yaml_cfg, "ignore_q_limits", base.ignore_q_limits)),
      max_ite = Int(get(yaml_cfg, "max_ite", base.max_ite)),
      tol = Float64(get(yaml_cfg, "tol", base.tol)),
      show_diff = Bool(get(yaml_cfg, "show_diff", base.show_diff)),
      tol_vm = Float64(get(yaml_cfg, "tol_vm", base.tol_vm)),
      tol_va = Float64(get(yaml_cfg, "tol_va", base.tol_va)),
      seconds = Float64(get(yaml_cfg, "seconds", base.seconds)),
      samples = Int(get(yaml_cfg, "samples", base.samples)),
      show_once = Bool(get(yaml_cfg, "show_once", base.show_once)),
      show_once_output = _as_output_mode(get(yaml_cfg, "show_once_output", base.show_once_output)),
      show_once_max_nodes = Int(get(yaml_cfg, "show_once_max_nodes", base.show_once_max_nodes)),
      benchmark = Bool(get(yaml_cfg, "benchmark", base.benchmark)),
      matpower_shift_sign = Float64(get(yaml_cfg, "matpower_shift_sign", base.matpower_shift_sign)),
      matpower_shift_unit = String(get(yaml_cfg, "matpower_shift_unit", base.matpower_shift_unit)),
      matpower_ratio = String(get(yaml_cfg, "matpower_ratio", base.matpower_ratio)),
      reference_override = Bool(get(yaml_cfg, "reference_override", base.reference_override)),
      reference_vm_pu = Float64(get(yaml_cfg, "reference_vm_pu", base.reference_vm_pu)),
      reference_va_deg = Float64(get(yaml_cfg, "reference_va_deg", base.reference_va_deg)),
      diagnose_matpower_reference = Bool(get(yaml_cfg, "diagnose_matpower_reference", base.diagnose_matpower_reference)),
      diagnose_branch_shift_conventions = Bool(get(yaml_cfg, "diagnose_branch_shift_conventions", base.diagnose_branch_shift_conventions)),
      diagnose_branch_neighborhood = Bool(get(yaml_cfg, "diagnose_branch_neighborhood", base.diagnose_branch_neighborhood)),
      diagnose_residual_clusters = Bool(get(yaml_cfg, "diagnose_residual_clusters", base.diagnose_residual_clusters)),
      diagnose_residual_cluster_threshold_mw = Float64(get(yaml_cfg, "diagnose_residual_cluster_threshold_mw", base.diagnose_residual_cluster_threshold_mw)),
      diagnose_residual_cluster_maxlines = Int(get(yaml_cfg, "diagnose_residual_cluster_maxlines", base.diagnose_residual_cluster_maxlines)),
      diagnose_branch_neighborhood_buses = _as_int_vec(get(yaml_cfg, "diagnose_branch_neighborhood_buses", base.diagnose_branch_neighborhood_buses)),
      diagnose_branch_neighborhood_depth = Int(get(yaml_cfg, "diagnose_branch_neighborhood_depth", base.diagnose_branch_neighborhood_depth)),
      diagnose_branch_neighborhood_maxlines = Int(get(yaml_cfg, "diagnose_branch_neighborhood_maxlines", base.diagnose_branch_neighborhood_maxlines)),
      diagnose_nodal_balance_breakdown = Bool(get(yaml_cfg, "diagnose_nodal_balance_breakdown", base.diagnose_nodal_balance_breakdown)),
      diagnose_nodal_balance_buses = _as_int_vec(get(yaml_cfg, "diagnose_nodal_balance_buses", base.diagnose_nodal_balance_buses)),
      diagnose_nodal_balance_maxlines = Int(get(yaml_cfg, "diagnose_nodal_balance_maxlines", base.diagnose_nodal_balance_maxlines)),
      diagnose_nodal_balance_include_branches = Bool(get(yaml_cfg, "diagnose_nodal_balance_include_branches", base.diagnose_nodal_balance_include_branches)),
      diagnose_nodal_balance_include_generators = Bool(get(yaml_cfg, "diagnose_nodal_balance_include_generators", base.diagnose_nodal_balance_include_generators)),
      diagnose_nodal_balance_include_shunts = Bool(get(yaml_cfg, "diagnose_nodal_balance_include_shunts", base.diagnose_nodal_balance_include_shunts)),
      diagnose_negative_branch_impedance = Bool(get(yaml_cfg, "diagnose_negative_branch_impedance", base.diagnose_negative_branch_impedance)),
      diagnose_negative_branch_impedance_maxlines = Int(get(yaml_cfg, "diagnose_negative_branch_impedance_maxlines", base.diagnose_negative_branch_impedance_maxlines)),
      diagnose_negative_branch_impedance_fail_on_negative_r = Bool(get(yaml_cfg, "diagnose_negative_branch_impedance_fail_on_negative_r", base.diagnose_negative_branch_impedance_fail_on_negative_r)),
      diagnose_negative_branch_impedance_fail_on_negative_x = Bool(get(yaml_cfg, "diagnose_negative_branch_impedance_fail_on_negative_x", base.diagnose_negative_branch_impedance_fail_on_negative_x)),
      diagnose_negative_branch_impedance_warn_threshold_abs_r = Float64(get(yaml_cfg, "diagnose_negative_branch_impedance_warn_threshold_abs_r", base.diagnose_negative_branch_impedance_warn_threshold_abs_r)),
      diagnose_negative_branch_impedance_warn_threshold_abs_x = Float64(get(yaml_cfg, "diagnose_negative_branch_impedance_warn_threshold_abs_x", base.diagnose_negative_branch_impedance_warn_threshold_abs_x)),
      diagnose_maxlines = Int(get(yaml_cfg, "diagnose_maxlines", base.diagnose_maxlines)),
      log_effective_config = Bool(get(yaml_cfg, "log_effective_config", base.log_effective_config)),
      enable_pq_gen_controllers = Bool(get(yaml_cfg, "enable_pq_gen_controllers", base.enable_pq_gen_controllers)),
      bus_shunt_model = String(get(yaml_cfg, "bus_shunt_model", base.bus_shunt_model)),
      trace_legacy_bus_type_warnings = Bool(get(yaml_cfg, "trace_legacy_bus_type_warnings", base.trace_legacy_bus_type_warnings)),
      matpower_pv_voltage_source = Symbol(get(yaml_cfg, "matpower_pv_voltage_source", base.matpower_pv_voltage_source)),
      matpower_pv_voltage_mismatch_tol_pu = Float64(get(yaml_cfg, "matpower_pv_voltage_mismatch_tol_pu", base.matpower_pv_voltage_mismatch_tol_pu)),
      compare_voltage_reference = Symbol(get(yaml_cfg, "compare_voltage_reference", base.compare_voltage_reference)),
      diagnose_pv_voltage_references = Bool(get(yaml_cfg, "diagnose_pv_voltage_references", base.diagnose_pv_voltage_references)),
      diagnose_pv_voltage_maxlines = Int(get(yaml_cfg, "diagnose_pv_voltage_maxlines", base.diagnose_pv_voltage_maxlines)),
      flatstart_voltage_mode = Symbol(get(yaml_cfg, "flatstart_voltage_mode", base.flatstart_voltage_mode)),
      flatstart_angle_mode = Symbol(get(yaml_cfg, "flatstart_angle_mode", base.flatstart_angle_mode)),
      flatstart_branch_guard = Bool(get(yaml_cfg, "flatstart_branch_guard", base.flatstart_branch_guard)),
      flatstart_min_vm_pu = Float64(get(yaml_cfg, "flatstart_min_vm_pu", base.flatstart_min_vm_pu)),
      flatstart_max_angle_step_deg = Float64(get(yaml_cfg, "flatstart_max_angle_step_deg", base.flatstart_max_angle_step_deg)),
      flatstart_max_vm_step_pu = Float64(get(yaml_cfg, "flatstart_max_vm_step_pu", base.flatstart_max_vm_step_pu)),
      wrong_branch_detection = Bool(get(yaml_cfg, "wrong_branch_detection", base.wrong_branch_detection)),
      wrong_branch_min_vm_pu = Float64(get(yaml_cfg, "wrong_branch_min_vm_pu", base.wrong_branch_min_vm_pu)),
      wrong_branch_max_angle_spread_deg = Float64(get(yaml_cfg, "wrong_branch_max_angle_spread_deg", base.wrong_branch_max_angle_spread_deg)),
      wrong_branch_rescue = Bool(get(yaml_cfg, "wrong_branch_rescue", base.wrong_branch_rescue)),
      wrong_branch_rescue_modes = _as_symbol_vec(get(yaml_cfg, "wrong_branch_rescue_modes", base.wrong_branch_rescue_modes)),
      matpower_auto_profile = _as_auto_profile_mode(get(yaml_cfg, "matpower_auto_profile", base.matpower_auto_profile)),
      matpower_auto_profile_log = Bool(get(yaml_cfg, "matpower_auto_profile_log", base.matpower_auto_profile_log)),
      console_summary = Bool(get(yaml_cfg, "console_summary", base.console_summary)),
      console_auto_profile = _as_console_mode(get(yaml_cfg, "console_auto_profile", base.console_auto_profile)),
      console_diagnostics = _as_console_mode(get(yaml_cfg, "console_diagnostics", base.console_diagnostics)),
      console_q_limit_events = _as_console_mode(get(yaml_cfg, "console_q_limit_events", base.console_q_limit_events)),
      console_max_rows = Int(get(yaml_cfg, "console_max_rows", base.console_max_rows)),
      logfile_diagnostics = _as_console_mode(get(yaml_cfg, "logfile_diagnostics", base.logfile_diagnostics)),
      performance_enabled = Bool(_yaml_section_get(yaml_cfg, "performance", "enabled", base.performance_enabled)),
      performance_level = _as_performance_level(_yaml_section_get(yaml_cfg, "performance", "level", base.performance_level)),
      performance_print_to_console = Bool(_yaml_section_get(yaml_cfg, "performance", "print_to_console", base.performance_print_to_console)),
      performance_write_to_logfile = Bool(_yaml_section_get(yaml_cfg, "performance", "write_to_logfile", base.performance_write_to_logfile)),
      performance_show_allocations = Bool(_yaml_section_get(yaml_cfg, "performance", "show_allocations", base.performance_show_allocations)),
      performance_show_iteration_table = Bool(_yaml_section_get(yaml_cfg, "performance", "show_iteration_table", base.performance_show_iteration_table)),
      performance_compact_logging = Bool(_yaml_section_get(yaml_cfg, "performance", "compact_logging", base.performance_compact_logging)),
      performance_skip_reference_comparison = Bool(_yaml_section_get(yaml_cfg, "performance", "skip_reference_comparison", base.performance_skip_reference_comparison)),
      performance_skip_expensive_diagnostics = Bool(_yaml_section_get(yaml_cfg, "performance", "skip_expensive_diagnostics", base.performance_skip_expensive_diagnostics)),
      performance_skip_branch_neighborhood_report = Bool(_yaml_section_get(yaml_cfg, "performance", "skip_branch_neighborhood_report", base.performance_skip_branch_neighborhood_report)),
      performance_max_diagnostic_rows = Int(_yaml_section_get(yaml_cfg, "performance", "max_diagnostic_rows", base.performance_max_diagnostic_rows)),
    )
  end
  return merge(base, case_override, yaml_override)
end

function _warn_if_flatstart_uses_only_voltage_setpoints(case_name::AbstractString, cfg, mpc)
  cfg.opt_flatstart && return nothing
  mp_has_vm_va(mpc) || return nothing

  @info "opt_flatstart=false uses stored MATPOWER voltage magnitudes and angles as the initial solve state. Set opt_flatstart=true to start from scratch with 1.0 pu / 0° on PQ buses and slack/PV voltage setpoints." case = case_name
  return nothing
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

function _matpower_reference_residuals(mpc; matpower_shift_sign::Real = 1.0, matpower_shift_unit = "deg", matpower_ratio = "normal", keep_shunts::Bool = true)
  bus = keep_shunts ? mpc.bus : copy(mpc.bus)
  if !keep_shunts && size(bus, 2) >= 6
    bus[:, 5] .= 0.0
    bus[:, 6] .= 0.0
  end
  branch = mpc.branch

  busrow = Dict{Int,Int}()
  sizehint!(busrow, size(bus, 1))
  for r in axes(bus, 1)
    busrow[Int(bus[r, 1])] = r
  end

  nbus = size(bus, 1)
  Pinj = zeros(Float64, nbus)
  Qinj = zeros(Float64, nbus)
  for r in axes(bus, 1)
    Pinj[r] -= bus[r, 3] / mpc.baseMVA
    Qinj[r] -= bus[r, 4] / mpc.baseMVA
  end
  for g in axes(mpc.gen, 1)
    if size(mpc.gen, 2) >= 8 && mpc.gen[g, 8] <= 0.0
      continue
    end
    busI = Int(mpc.gen[g, 1])
    haskey(busrow, busI) || continue
    r = busrow[busI]
    Pinj[r] += mpc.gen[g, 2] / mpc.baseMVA
    Qinj[r] += mpc.gen[g, 3] / mpc.baseMVA
  end

  Y = MatpowerIO.build_ybus_matpower(bus, branch, mpc.baseMVA; matpower_shift_sign = matpower_shift_sign, matpower_shift_unit = matpower_shift_unit, matpower_ratio = matpower_ratio)
  V = bus[:, 8] .* cis.(deg2rad.(bus[:, 9]))
  mis = V .* conj.(Y * V) .- ComplexF64.(Pinj, Qinj)
  p_rows = Int[]
  q_rows = Int[]
  for r in axes(bus, 1)
    btype = Int(bus[r, 2])
    if btype == 1
      push!(p_rows, r)
      push!(q_rows, r)
    elseif btype == 2
      push!(p_rows, r)
    end
  end
  return (; bus, mis, p_rows, q_rows)
end

function _print_top_residuals(io::IO, label::AbstractString, diag, baseMVA::Float64; maxlines::Int = 12)
  bus = diag.bus
  mis = diag.mis
  p_order = sort(diag.p_rows; by = r -> abs(real(mis[r])), rev = true)
  q_order = sort(diag.q_rows; by = r -> abs(imag(mis[r])), rev = true)
  n_p = min(maxlines, length(p_order))
  n_q = min(maxlines, length(q_order))

  println(io, "\n==================== MATPOWER reference residual diagnostics: ", label, " ====================")
  println(io, "Top active-power residuals on enforced PQ/PV buses:")
  println(io, " rank  bus      type        dP_pu         dP_MW       Vm_ref     Va_ref_deg")
  for rank = 1:n_p
    r = p_order[rank]
    @printf(io, "%5d %6d %5d %13.6f %13.3f %10.6f %12.6f\n", rank, Int(bus[r, 1]), Int(bus[r, 2]), real(mis[r]), real(mis[r]) * baseMVA, bus[r, 8], bus[r, 9])
  end
  println(io, "Top reactive-power residuals on enforced PQ buses:")
  println(io, " rank  bus      type        dQ_pu       dQ_MVAr      Vm_ref     Va_ref_deg")
  for rank = 1:n_q
    r = q_order[rank]
    @printf(io, "%5d %6d %5d %13.6f %13.3f %10.6f %12.6f\n", rank, Int(bus[r, 1]), Int(bus[r, 2]), imag(mis[r]), imag(mis[r]) * baseMVA, bus[r, 8], bus[r, 9])
  end
  println(io, "================================================================================\n")
end

function _matpower_busrow_map(bus)
  busrow = Dict{Int,Int}()
  sizehint!(busrow, size(bus, 1))
  for r in axes(bus, 1)
    busrow[Int(bus[r, 1])] = r
  end
  return busrow
end

function _matpower_branch_stamp(row; matpower_shift_sign::Real = 1.0, matpower_shift_unit = "deg", matpower_ratio = "normal")
  shift_unit = MatpowerIO._normalize_matpower_shift_unit(matpower_shift_unit)
  shift_sign = Float64(matpower_shift_sign)
  ratio_mode = MatpowerIO._normalize_matpower_ratio_mode(matpower_ratio)

  r = Float64(row[3])
  x = Float64(row[4])
  b = Float64(row[5])
  y = inv(complex(r, x))
  ysh = 1im * b / 2
  tap_raw = length(row) >= 9 ? Float64(row[9]) : 1.0
  shift_raw = length(row) >= 10 ? Float64(row[10]) : 0.0
  ratio = MatpowerIO._matpower_import_ratio(tap_raw; mode = ratio_mode)
  shift_rad = MatpowerIO._matpower_shift_radians(shift_raw; sign = shift_sign, unit = shift_unit)
  tap = ratio * cis(shift_rad)

  return (;
    br_r = r,
    br_x = x,
    br_b = b,
    tap_raw = tap_raw,
    shift_raw = shift_raw,
    ratio = ratio,
    shift_rad = shift_rad,
    shift_deg = rad2deg(shift_rad),
    tap = tap,
    Yff = (y + ysh) / (tap * conj(tap)),
    Yft = -y / conj(tap),
    Ytf = -y / tap,
    Ytt = y + ysh,
  )
end

function _matpower_branch_flow_from_stamp(mpc, busrow::Dict{Int,Int}, e::Int, stamp)
  row = view(mpc.branch, e, :)
  f_bus = Int(row[1])
  t_bus = Int(row[2])
  f = busrow[f_bus]
  t = busrow[t_bus]
  Vf = mpc.bus[f, 8] * cis(deg2rad(mpc.bus[f, 9]))
  Vt = mpc.bus[t, 8] * cis(deg2rad(mpc.bus[t, 9]))
  If = stamp.Yff * Vf + stamp.Yft * Vt
  It = stamp.Ytf * Vf + stamp.Ytt * Vt
  return (; f_bus, t_bus, f, t, Vf, Vt, Sf = Vf * conj(If), St = Vt * conj(It))
end

function _matpower_high_residual_bus_set(diag; maxlines::Int = 12)
  high = Set{Int}()
  p_order = sort(diag.p_rows; by = r -> abs(real(diag.mis[r])), rev = true)
  q_order = sort(diag.q_rows; by = r -> abs(imag(diag.mis[r])), rev = true)
  for r in p_order[1:min(max(maxlines, 0), length(p_order))]
    push!(high, Int(diag.bus[r, 1]))
  end
  for r in q_order[1:min(max(maxlines, 0), length(q_order))]
    push!(high, Int(diag.bus[r, 1]))
  end
  return high
end

function _matpower_residual_bus_rows(diag, baseMVA::Float64; threshold_mw::Real = 1.0, maxlines::Int = 20)
  rows = Set{Int}()
  threshold_pu = abs(Float64(threshold_mw)) / baseMVA
  for r in diag.p_rows
    abs(real(diag.mis[r])) >= threshold_pu && push!(rows, r)
  end
  for r in diag.q_rows
    abs(imag(diag.mis[r])) >= threshold_pu && push!(rows, r)
  end
  if isempty(rows)
    p_order = sort(diag.p_rows; by = r -> abs(real(diag.mis[r])), rev = true)
    q_order = sort(diag.q_rows; by = r -> abs(imag(diag.mis[r])), rev = true)
    for r in p_order[1:min(max(maxlines, 0), length(p_order))]
      push!(rows, r)
    end
    for r in q_order[1:min(max(maxlines, 0), length(q_order))]
      push!(rows, r)
    end
  end
  return rows
end

function _matpower_residual_bus_clusters(mpc, diag; threshold_mw::Real = 1.0, maxlines::Int = 20)
  selected_rows = _matpower_residual_bus_rows(diag, Float64(mpc.baseMVA); threshold_mw = threshold_mw, maxlines = maxlines)
  isempty(selected_rows) && return NamedTuple[]
  busrow = _matpower_busrow_map(diag.bus)
  selected_buses = Set(Int(diag.bus[r, 1]) for r in selected_rows)
  adjacency = Dict{Int,Set{Int}}(busI => Set{Int}() for busI in selected_buses)
  for e in axes(mpc.branch, 1)
    status = size(mpc.branch, 2) >= 11 ? Float64(mpc.branch[e, 11]) : 1.0
    status == 0.0 && continue
    f_bus = Int(mpc.branch[e, 1])
    t_bus = Int(mpc.branch[e, 2])
    if f_bus in selected_buses && t_bus in selected_buses
      push!(adjacency[f_bus], t_bus)
      push!(adjacency[t_bus], f_bus)
    end
  end

  clusters = NamedTuple[]
  seen = Set{Int}()
  for seed in sort!(collect(selected_buses))
    seed in seen && continue
    queue = [seed]
    push!(seen, seed)
    buses = Int[]
    while !isempty(queue)
      busI = popfirst!(queue)
      push!(buses, busI)
      for nb in adjacency[busI]
        nb in seen && continue
        push!(seen, nb)
        push!(queue, nb)
      end
    end
    sort!(buses)
    row_ids = [busrow[busI] for busI in buses if haskey(busrow, busI)]
    internal_branch_rows = Int[]
    boundary_branch_rows = Int[]
    bus_set = Set(buses)
    for e in axes(mpc.branch, 1)
      status = size(mpc.branch, 2) >= 11 ? Float64(mpc.branch[e, 11]) : 1.0
      status == 0.0 && continue
      f_bus = Int(mpc.branch[e, 1])
      t_bus = Int(mpc.branch[e, 2])
      f_in = f_bus in bus_set
      t_in = t_bus in bus_set
      if f_in && t_in
        push!(internal_branch_rows, e)
      elseif f_in || t_in
        push!(boundary_branch_rows, e)
      end
    end
    max_p_mw = isempty(row_ids) ? 0.0 : maximum(abs(real(diag.mis[r])) for r in row_ids) * mpc.baseMVA
    max_q_mvar = isempty(row_ids) ? 0.0 : maximum(abs(imag(diag.mis[r])) for r in row_ids) * mpc.baseMVA
    sum_mis = isempty(row_ids) ? 0.0 + 0.0im : sum(diag.mis[row_ids])
    push!(clusters, (;
      buses = buses,
      rows = row_ids,
      internal_branch_rows = internal_branch_rows,
      boundary_branch_rows = boundary_branch_rows,
      max_p_mw = max_p_mw,
      max_q_mvar = max_q_mvar,
      sum_mis = sum_mis,
    ))
  end
  return sort!(clusters; by = c -> max(c.max_p_mw, c.max_q_mvar), rev = true)
end

function _print_residual_cluster_diagnostics(io::IO, mpc, diag; threshold_mw::Real = 1.0, maxlines::Int = 20, matpower_shift_sign::Real = 1.0, matpower_shift_unit = "deg", matpower_ratio = "normal")
  clusters = _matpower_residual_bus_clusters(mpc, diag; threshold_mw = threshold_mw, maxlines = maxlines)
  isempty(clusters) && return nothing
  nshow = min(max(maxlines, 0), length(clusters))
  println(io, "==================== MATPOWER fixed-reference residual clusters ====================")
  @printf(io, "threshold: %.6g MW/MVAr; clusters: %d; showing: %d\n", Float64(threshold_mw), length(clusters), nshow)
  println(io, "Clusters connect buses whose fixed-reference P/Q residuals exceed the threshold through online MATPOWER branches.")
  println(io, "Use this to separate local data/model inconsistencies from one global residual list.")
  for rank in 1:nshow
    cluster = clusters[rank]
    @printf(io, "\ncluster %d: buses=%d internal_branches=%d boundary_branches=%d max|dP|=%.3f MW max|dQ|=%.3f MVAr sum_dS=% .6g%+ .6gi pu\n", rank, length(cluster.buses), length(cluster.internal_branch_rows), length(cluster.boundary_branch_rows), cluster.max_p_mw, cluster.max_q_mvar, real(cluster.sum_mis), imag(cluster.sum_mis))
    println(io, "  buses: ", join(cluster.buses[1:min(length(cluster.buses), maxlines)], ", "), length(cluster.buses) > maxlines ? " ..." : "")
    bus_order = sort(cluster.rows; by = r -> max(abs(real(diag.mis[r])), abs(imag(diag.mis[r]))), rev = true)
    println(io, "  top residual buses: BUS_I type dP_MW dQ_MVAr Vm_ref Va_ref_deg")
    for r in bus_order[1:min(length(bus_order), maxlines)]
      @printf(io, "    %6d %4d % .6f % .6f % .6f % .6f\n", Int(diag.bus[r, 1]), Int(diag.bus[r, 2]), real(diag.mis[r]) * mpc.baseMVA, imag(diag.mis[r]) * mpc.baseMVA, diag.bus[r, 8], diag.bus[r, 9])
    end
    if !isempty(cluster.internal_branch_rows)
      println(io, "  internal branch rows: ", join(cluster.internal_branch_rows[1:min(length(cluster.internal_branch_rows), maxlines)], ", "), length(cluster.internal_branch_rows) > maxlines ? " ..." : "")
      busrow = _matpower_busrow_map(mpc.bus)
      for e in cluster.internal_branch_rows[1:min(length(cluster.internal_branch_rows), maxlines)]
        row = view(mpc.branch, e, :)
        stamp = _matpower_branch_stamp(row; matpower_shift_sign = matpower_shift_sign, matpower_shift_unit = matpower_shift_unit, matpower_ratio = matpower_ratio)
        flow = _matpower_branch_flow_from_stamp(mpc, busrow, e, stamp)
        @printf(io, "    row %d %d->%d R=% .6g X=% .6g TAP=% .6g SHIFT=% .6g Sf=% .6g%+ .6gi pu St=% .6g%+ .6gi pu\n", e, flow.f_bus, flow.t_bus, stamp.br_r, stamp.br_x, stamp.tap_raw, stamp.shift_raw, real(flow.Sf), imag(flow.Sf), real(flow.St), imag(flow.St))
      end
    end
  end
  if nshow < length(clusters)
    @printf(io, "\n... truncated residual cluster diagnostics: showing %d/%d clusters (diagnose_residual_cluster_maxlines)\n", nshow, length(clusters))
  end
  println(io, "==================================================================================\n")
  return nothing
end

function _negative_branch_markers(row, selected_set::Set{Int}; high_residual_buses::Set{Int} = Set{Int}())
  f_bus = Int(row[1])
  t_bus = Int(row[2])
  markers = String[]
  Float64(row[3]) < 0.0 && push!(markers, "<negative BR_R>")
  Float64(row[4]) < 0.0 && push!(markers, "<negative BR_X>")
  (f_bus in selected_set && t_bus in selected_set) && push!(markers, "<connects two selected high-residual buses>")
  (f_bus in high_residual_buses && t_bus in high_residual_buses) && push!(markers, "<connects two top-residual buses>")
  return isempty(markers) ? "" : "  " * join(markers, " ")
end

function _negative_branch_rows(branch; warn_threshold_abs_r::Real = 0.0, warn_threshold_abs_x::Real = 0.0)
  rows = Int[]
  for e in axes(branch, 1)
    r = Float64(branch[e, 3])
    x = Float64(branch[e, 4])
    if r < -abs(warn_threshold_abs_r) || x < -abs(warn_threshold_abs_x)
      push!(rows, e)
    end
  end
  return rows
end

function _print_negative_branch_impedance_diagnostics(io::IO, mpc; matpower_shift_sign::Real = 1.0, matpower_shift_unit = "deg", matpower_ratio = "normal", maxlines::Int = 100, fail_on_negative_r::Bool = false, fail_on_negative_x::Bool = false, warn_threshold_abs_r::Real = 0.0, warn_threshold_abs_x::Real = 0.0, diagnose_branch_neighborhood_buses::AbstractVector{Int} = Int[], diagnose_nodal_balance_buses::AbstractVector{Int} = Int[], high_residual_buses::Set{Int} = Set{Int}())
  branch = mpc.branch
  nbranch = size(branch, 1)
  neg_r = [e for e in axes(branch, 1) if Float64(branch[e, 3]) < 0.0]
  neg_x = [e for e in axes(branch, 1) if Float64(branch[e, 4]) < 0.0]
  both_neg = [e for e in axes(branch, 1) if Float64(branch[e, 3]) < 0.0 && Float64(branch[e, 4]) < 0.0]
  zero_z = [e for e in axes(branch, 1) if Float64(branch[e, 3]) == 0.0 && Float64(branch[e, 4]) == 0.0]
  tap_zero = [e for e in axes(branch, 1) if (size(branch, 2) < 9 || Float64(branch[e, 9]) == 0.0)]
  tap_nonzero = [e for e in axes(branch, 1) if size(branch, 2) >= 9 && Float64(branch[e, 9]) != 0.0]
  shift_nonzero = [e for e in axes(branch, 1) if size(branch, 2) >= 10 && Float64(branch[e, 10]) != 0.0]
  neg_r_shift = [e for e in neg_r if size(branch, 2) >= 10 && Float64(branch[e, 10]) != 0.0]
  neg_r_nonunity_tap = [e for e in neg_r if size(branch, 2) >= 9 && Float64(branch[e, 9]) != 0.0 && Float64(branch[e, 9]) != 1.0]

  if fail_on_negative_r && !isempty(neg_r)
    e = neg_r[1]
    error("Negative MATPOWER branch resistance detected at branch row $(e) (f_bus=$(Int(branch[e, 1])), t_bus=$(Int(branch[e, 2])), BR_R=$(branch[e, 3])). Set diagnose_negative_branch_impedance_fail_on_negative_r=false to allow this model-critical data.")
  end
  if fail_on_negative_x && !isempty(neg_x)
    e = neg_x[1]
    error("Negative MATPOWER branch reactance detected at branch row $(e) (f_bus=$(Int(branch[e, 1])), t_bus=$(Int(branch[e, 2])), BR_X=$(branch[e, 4])). Set diagnose_negative_branch_impedance_fail_on_negative_x=false to allow this model-critical data.")
  end

  rows_to_print = _negative_branch_rows(branch; warn_threshold_abs_r = warn_threshold_abs_r, warn_threshold_abs_x = warn_threshold_abs_x)
  isempty(rows_to_print) && isempty(zero_z) && return nothing

  min_r_row = nbranch == 0 ? 0 : argmin(branch[:, 3])
  min_x_row = nbranch == 0 ? 0 : argmin(branch[:, 4])
  neighborhood_set = Set(Int[x for x in diagnose_branch_neighborhood_buses])
  nodal_set = Set(Int[x for x in diagnose_nodal_balance_buses])
  busrow = _matpower_busrow_map(mpc.bus)
  nshow = min(max(maxlines, 0), length(rows_to_print))

  println(io, "==================== MATPOWER negative branch impedance diagnostics ====================")
  if !isempty(neg_r) || !isempty(neg_x)
    println(io, "Negative branch resistance/reactance detected. This can occur in reduced or equivalent network models, including some PEGASE-style cases, but it is numerically and physically sensitive. Sparlectra preserves the signed impedance values. No clipping or correction is applied.")
  end
  println(io, "branch rows: ", nbranch)
  println(io, "BR_R < 0: ", length(neg_r), "   BR_X < 0: ", length(neg_x), "   both negative: ", length(both_neg))
  println(io, "BR_R == 0 and BR_X == 0: ", length(zero_z))
  println(io, "TAP == 0: ", length(tap_zero), "   nonzero TAP: ", length(tap_nonzero), "   nonzero SHIFT: ", length(shift_nonzero))
  println(io, "negative R with nonzero SHIFT: ", length(neg_r_shift), "   negative R with non-unity TAP: ", length(neg_r_nonunity_tap))
  if nbranch > 0
    @printf(io, "minimum BR_R: row %d  f_bus=%d  t_bus=%d  BR_R=% .10g\n", min_r_row, Int(branch[min_r_row, 1]), Int(branch[min_r_row, 2]), branch[min_r_row, 3])
    @printf(io, "minimum BR_X: row %d  f_bus=%d  t_bus=%d  BR_X=% .10g\n", min_x_row, Int(branch[min_x_row, 1]), Int(branch[min_x_row, 2]), branch[min_x_row, 4])
  end

  for pos in 1:nshow
    e = rows_to_print[pos]
    row = view(branch, e, :)
    f_bus = Int(row[1])
    t_bus = Int(row[2])
    stamp = _matpower_branch_stamp(row; matpower_shift_sign = matpower_shift_sign, matpower_shift_unit = matpower_shift_unit, matpower_ratio = matpower_ratio)
    flags = String[]
    (f_bus in neighborhood_set || t_bus in neighborhood_set) && push!(flags, "branch-neighborhood endpoint")
    (f_bus in nodal_set || t_bus in nodal_set) && push!(flags, "nodal-balance endpoint")
    (f_bus in high_residual_buses && t_bus in high_residual_buses) && push!(flags, "both endpoints high-residual")
    marker = _negative_branch_markers(row, neighborhood_set; high_residual_buses = high_residual_buses)
    @printf(io, "\nbranch row %d: f_bus=%d  t_bus=%d%s\n", e, f_bus, t_bus, marker)
    @printf(io, "  raw branch: BR_R=% .10g  BR_X=% .10g  BR_B=% .10g  TAP=% .10g  SHIFT=% .10g\n", stamp.br_r, stamp.br_x, stamp.br_b, stamp.tap_raw, stamp.shift_raw)
    @printf(io, "  interpreted: ratio=% .10g  shift=% .10g rad (% .10g deg)  t=% .10g%+ .10gi\n", stamp.ratio, stamp.shift_rad, stamp.shift_deg, real(stamp.tap), imag(stamp.tap))
    @printf(io, "  y=1/(R+jX)=% .10g%+ .10gi\n", real(inv(complex(stamp.br_r, stamp.br_x))), imag(inv(complex(stamp.br_r, stamp.br_x))))
    @printf(io, "  Yff=% .10g%+ .10gi  Yft=% .10g%+ .10gi\n", real(stamp.Yff), imag(stamp.Yff), real(stamp.Yft), imag(stamp.Yft))
    @printf(io, "  Ytf=% .10g%+ .10gi  Ytt=% .10g%+ .10gi\n", real(stamp.Ytf), imag(stamp.Ytf), real(stamp.Ytt), imag(stamp.Ytt))
    isempty(flags) || println(io, "  flags: ", join(flags, "; "))
    if haskey(busrow, f_bus) && haskey(busrow, t_bus)
      flow = _matpower_branch_flow_from_stamp(mpc, busrow, e, stamp)
      @printf(io, "  fixed-reference flow Sf=% .10g%+ .10gi pu  St=% .10g%+ .10gi pu\n", real(flow.Sf), imag(flow.Sf), real(flow.St), imag(flow.St))
    end
  end
  if nshow < length(rows_to_print)
    @printf(io, "\n... truncated negative branch impedance diagnostics: showing %d/%d branches (diagnose_negative_branch_impedance_maxlines)\n", nshow, length(rows_to_print))
  end
  println(io, "=======================================================================================\n")
  return nothing
end

function _matpower_branch_neighborhood(branch, seed_buses::AbstractVector{Int}; depth::Int = 1)
  seed_set = Set(seed_buses)
  isempty(seed_set) && return Int[]
  max_depth = max(depth, 1)
  seen_buses = copy(seed_set)
  frontier = copy(seed_set)
  selected = Set{Int}()

  for _depth in 1:max_depth
    next_frontier = Set{Int}()
    for e in axes(branch, 1)
      status = size(branch, 2) >= 11 ? branch[e, 11] : 1.0
      status == 0.0 && continue
      f_bus = Int(branch[e, 1])
      t_bus = Int(branch[e, 2])
      if f_bus in frontier || t_bus in frontier
        push!(selected, e)
        if !(f_bus in seen_buses)
          push!(seen_buses, f_bus)
          push!(next_frontier, f_bus)
        end
        if !(t_bus in seen_buses)
          push!(seen_buses, t_bus)
          push!(next_frontier, t_bus)
        end
      end
    end
    isempty(next_frontier) && break
    frontier = next_frontier
  end
  return sort!(collect(selected))
end

function _print_matpower_branch_neighborhood_diagnostics(io::IO, mpc, diag; buses::AbstractVector{Int}, depth::Int = 1, maxlines::Int = 200, matpower_shift_sign::Real = 1.0, matpower_shift_unit = "deg", matpower_ratio = "normal")
  mp_has_vm_va(mpc) || return nothing
  selected_buses = unique(Int[x for x in buses])
  isempty(selected_buses) && return nothing
  busrow = _matpower_busrow_map(mpc.bus)
  branch_rows = _matpower_branch_neighborhood(mpc.branch, selected_buses; depth = depth)
  selected_set = Set(selected_buses)
  high_residual_buses = _matpower_high_residual_bus_set(diag; maxlines = maxlines)
  nshow = min(max(maxlines, 0), length(branch_rows))

  println(io, "==================== MATPOWER branch-neighborhood fixed-reference diagnostics ====================")
  println(io, "seed buses: ", selected_buses)
  println(io, "depth: ", max(depth, 1), "   selected branches: ", length(branch_rows), "   showing: ", nshow)
  println(io, "convention: SHIFT sign=", matpower_shift_sign, " unit=", matpower_shift_unit, " ratio=", matpower_ratio)

  for pos in 1:nshow
    e = branch_rows[pos]
    row = view(mpc.branch, e, :)
    f_bus = Int(row[1])
    t_bus = Int(row[2])
    haskey(busrow, f_bus) || continue
    haskey(busrow, t_bus) || continue
    f = busrow[f_bus]
    t = busrow[t_bus]
    stamp = _matpower_branch_stamp(row; matpower_shift_sign = matpower_shift_sign, matpower_shift_unit = matpower_shift_unit, matpower_ratio = matpower_ratio)
    Vf = mpc.bus[f, 8] * cis(deg2rad(mpc.bus[f, 9]))
    Vt = mpc.bus[t, 8] * cis(deg2rad(mpc.bus[t, 9]))
    If = stamp.Yff * Vf + stamp.Yft * Vt
    It = stamp.Ytf * Vf + stamp.Ytt * Vt
    Sf = Vf * conj(If)
    St = Vt * conj(It)
    marker = _negative_branch_markers(row, selected_set; high_residual_buses = high_residual_buses)

    @printf(io, "\nbranch row %d: f_bus=%d  t_bus=%d%s\n", e, f_bus, t_bus, marker)
    @printf(io, "  raw branch: BR_R=% .10g  BR_X=% .10g  BR_B=% .10g  TAP=% .10g  SHIFT=% .10g\n", stamp.br_r, stamp.br_x, stamp.br_b, stamp.tap_raw, stamp.shift_raw)
    @printf(io, "  interpreted: ratio=% .10g  shift=% .10g rad (% .10g deg)  t=% .10g%+ .10gi\n", stamp.ratio, stamp.shift_rad, stamp.shift_deg, real(stamp.tap), imag(stamp.tap))
    @printf(io, "  Yff=% .10g%+ .10gi  Yft=% .10g%+ .10gi\n", real(stamp.Yff), imag(stamp.Yff), real(stamp.Yft), imag(stamp.Yft))
    @printf(io, "  Ytf=% .10g%+ .10gi  Ytt=% .10g%+ .10gi\n", real(stamp.Ytf), imag(stamp.Ytf), real(stamp.Ytt), imag(stamp.Ytt))
    @printf(io, "  Vref f: Vm=% .10g pu  Va=% .10g deg;  Vref t: Vm=% .10g pu  Va=% .10g deg\n", mpc.bus[f, 8], mpc.bus[f, 9], mpc.bus[t, 8], mpc.bus[t, 9])
    @printf(io, "  fixed-reference flow Sf=% .10g%+ .10gi pu (% .6f MW,% .6f MVAr)\n", real(Sf), imag(Sf), real(Sf) * mpc.baseMVA, imag(Sf) * mpc.baseMVA)
    @printf(io, "  fixed-reference flow St=% .10g%+ .10gi pu (% .6f MW,% .6f MVAr)\n", real(St), imag(St), real(St) * mpc.baseMVA, imag(St) * mpc.baseMVA)
    @printf(io, "  residual contribution at f_bus: dP=% .10g pu (% .6f MW)  dQ=% .10g pu (% .6f MVAr); total residual dP=% .10g pu dQ=% .10g pu\n", real(Sf), real(Sf) * mpc.baseMVA, imag(Sf), imag(Sf) * mpc.baseMVA, real(diag.mis[f]), imag(diag.mis[f]))
    @printf(io, "  residual contribution at t_bus: dP=% .10g pu (% .6f MW)  dQ=% .10g pu (% .6f MVAr); total residual dP=% .10g pu dQ=% .10g pu\n", real(St), real(St) * mpc.baseMVA, imag(St), imag(St) * mpc.baseMVA, real(diag.mis[t]), imag(diag.mis[t]))
  end
  if nshow < length(branch_rows)
    @printf(io, "\n... truncated branch-neighborhood diagnostics: showing %d/%d branches (diagnose_branch_neighborhood_maxlines)\n", nshow, length(branch_rows))
  end
  println(io, "===================================================================================================\n")
  return nothing
end

function _matpower_bus_type_label(code::Int)::String
  code == 1 && return "PQ"
  code == 2 && return "PV"
  code == 3 && return "REF"
  code == 4 && return "isolated"
  return string("unknown(", code, ")")
end

function _matpower_voltage_setpoint_source(mpc, busI::Int, busrow::Int; matpower_pv_voltage_source = :gen_vg)::String
  if matpower_pv_voltage_source == :bus_vm
    return "BUS.VM"
  end
  for g in axes(mpc.gen, 1)
    if Int(mpc.gen[g, 1]) == busI && (size(mpc.gen, 2) < 8 || mpc.gen[g, 8] > 0.0)
      return "GEN.VG"
    end
  end
  return "BUS.VM fallback"
end

function _nodal_balance_dominant_source(Sbranch::ComplexF64, Sspec_no_shunt::ComplexF64, Sshunt::ComplexF64, residual::ComplexF64; bus_type::Int)::String
  mag_res = abs(residual)
  mag_res < 1e-9 && return "mixed"
  parts = (branch_sum = abs(Sbranch), load_or_gen_spec = abs(Sspec_no_shunt), bus_shunt = abs(Sshunt))
  max_part = maximum(values(parts))
  max_part < 1e-12 && return bus_type in (2, 3) ? "q_treatment" : "unknown"
  winners = [String(k) for (k, v) in pairs(parts) if v >= 0.67 * max_part]
  length(winners) == 1 && return winners[1]
  return bus_type in (2, 3) && abs(imag(residual)) > abs(real(residual)) ? "q_treatment" : "mixed"
end

function _print_nodal_balance_breakdown(io::IO, mpc, diag; buses::AbstractVector{Int}, maxlines::Int = 200, include_branches::Bool = true, include_generators::Bool = true, include_shunts::Bool = true, matpower_shift_sign::Real = 1.0, matpower_shift_unit = "deg", matpower_ratio = "normal", matpower_pv_voltage_source = :gen_vg)
  mp_has_vm_va(mpc) || return nothing
  selected_buses = unique(Int[x for x in buses])
  isempty(selected_buses) && return nothing
  busrow = _matpower_busrow_map(mpc.bus)
  nshow = min(max(maxlines, 0), length(selected_buses))
  branch_sums = zeros(ComplexF64, size(mpc.bus, 1))
  incident = Dict{Int,Vector{Tuple{Int,Symbol,ComplexF64}}}()
  for e in axes(mpc.branch, 1)
    status = size(mpc.branch, 2) >= 11 ? mpc.branch[e, 11] : 1.0
    status == 0.0 && continue
    f_bus = Int(mpc.branch[e, 1])
    t_bus = Int(mpc.branch[e, 2])
    haskey(busrow, f_bus) || continue
    haskey(busrow, t_bus) || continue
    stamp = _matpower_branch_stamp(view(mpc.branch, e, :); matpower_shift_sign = matpower_shift_sign, matpower_shift_unit = matpower_shift_unit, matpower_ratio = matpower_ratio)
    flow = _matpower_branch_flow_from_stamp(mpc, busrow, e, stamp)
    branch_sums[flow.f] += flow.Sf
    branch_sums[flow.t] += flow.St
    push!(get!(incident, f_bus, Tuple{Int,Symbol,ComplexF64}[]), (e, :from, flow.Sf))
    push!(get!(incident, t_bus, Tuple{Int,Symbol,ComplexF64}[]), (e, :to, flow.St))
  end

  println(io, "==================== MATPOWER fixed-reference nodal balance breakdown ====================")
  println(io, "formula: S_spec = S_gen_online - S_load - S_shunt; S_branch_sum = sum incident fixed-reference branch injections; residual dS = S_branch_sum - S_spec")
  println(io, "units: p.u. on baseMVA=", mpc.baseMVA, " (MW/MVAr values are p.u. times baseMVA)")
  for busI in selected_buses[1:nshow]
    if !haskey(busrow, busI)
      println(io, "\nBUS_I=", busI, " not found in MATPOWER bus table")
      continue
    end
    r = busrow[busI]
    btype = Int(mpc.bus[r, 2])
    Pd = mpc.bus[r, 3] / mpc.baseMVA
    Qd = mpc.bus[r, 4] / mpc.baseMVA
    Vm = mpc.bus[r, 8]
    Va = mpc.bus[r, 9]
    Sload = complex(Pd, Qd)
    Sshunt = (Vm^2) * complex(mpc.bus[r, 5], -mpc.bus[r, 6]) / mpc.baseMVA
    online_gen = 0 + 0im
    offline_gen = 0 + 0im
    online_count = 0
    offline_count = 0
    for g in axes(mpc.gen, 1)
      Int(mpc.gen[g, 1]) == busI || continue
      Sg = complex(mpc.gen[g, 2], mpc.gen[g, 3]) / mpc.baseMVA
      if size(mpc.gen, 2) >= 8 && mpc.gen[g, 8] <= 0.0
        offline_gen += Sg
        offline_count += 1
      else
        online_gen += Sg
        online_count += 1
      end
    end
    Sspec_no_shunt = online_gen - Sload
    Sspec = Sspec_no_shunt - Sshunt
    Sbranch = branch_sums[r]
    residual = Sbranch - Sspec
    total_residual = diag.mis[r]
    source = _nodal_balance_dominant_source(Sbranch, Sspec_no_shunt, Sshunt, residual; bus_type = btype)

    @printf(io, "\nBUS_I=%d  internal_index=%d  type=%s  Vm_ref=% .10g  Va_ref=% .10g deg\n", busI, r, _matpower_bus_type_label(btype), Vm, Va)
    @printf(io, "  Pd=% .10g pu (% .6f MW)  Qd=% .10g pu (% .6f MVAr)\n", Pd, Pd * mpc.baseMVA, Qd, Qd * mpc.baseMVA)
    @printf(io, "  Gs=% .10g MW at V=1  Bs=% .10g MVAr at V=1  S_shunt=% .10g%+ .10gi pu\n", mpc.bus[r, 5], mpc.bus[r, 6], real(Sshunt), imag(Sshunt))
    if include_generators
      @printf(io, "  online generators: count=%d  sum Pg=% .10g pu (% .6f MW)  sum Qg=% .10g pu (% .6f MVAr)\n", online_count, real(online_gen), real(online_gen) * mpc.baseMVA, imag(online_gen), imag(online_gen) * mpc.baseMVA)
      if offline_count > 0 || abs(offline_gen) > 0.0
        @printf(io, "  offline generators: count=%d  sum Pg=% .10g pu (% .6f MW)  sum Qg=% .10g pu (% .6f MVAr)\n", offline_count, real(offline_gen), real(offline_gen) * mpc.baseMVA, imag(offline_gen), imag(offline_gen) * mpc.baseMVA)
      end
    end
    if btype in (2, 3)
      println(io, "  PV/REF treatment: P is enforced in the MATPOWER power-flow equations; Q is diagnostic/calculated here and is not part of the checked residual unless the bus is PQ.")
      println(io, "  voltage setpoint source: ", _matpower_voltage_setpoint_source(mpc, busI, r; matpower_pv_voltage_source = matpower_pv_voltage_source))
    end
    include_shunts || println(io, "  shunt details suppressed by diagnose_nodal_balance_include_shunts=false; S_shunt is still included in the formula above.")
    if include_branches
      entries = get(incident, busI, Tuple{Int,Symbol,ComplexF64}[])
      println(io, "  incident online branches: ", length(entries))
      for (e, side, S) in entries[1:min(length(entries), max(0, maxlines))]
        @printf(io, "    branch row %d (%s): S=% .10g%+ .10gi pu (% .6f MW,% .6f MVAr)\n", e, String(side), real(S), imag(S), real(S) * mpc.baseMVA, imag(S) * mpc.baseMVA)
      end
    end
    @printf(io, "  S_spec=% .10g%+ .10gi pu  (S_gen_online - S_load - S_shunt)\n", real(Sspec), imag(Sspec))
    @printf(io, "  S_branch_sum=% .10g%+ .10gi pu\n", real(Sbranch), imag(Sbranch))
    @printf(io, "  residual dS=S_branch_sum-S_spec=% .10g%+ .10gi pu (% .6f MW,% .6f MVAr)\n", real(residual), imag(residual), real(residual) * mpc.baseMVA, imag(residual) * mpc.baseMVA)
    @printf(io, "  reported total residual comparison: diag.mis=% .10g%+ .10gi pu  delta=% .3e%+ .3ei pu\n", real(total_residual), imag(total_residual), real(residual - total_residual), imag(residual - total_residual))
    println(io, "  dominant source classification: ", source)
  end
  if nshow < length(selected_buses)
    @printf(io, "\n... truncated nodal balance diagnostics: showing %d/%d buses (diagnose_nodal_balance_maxlines)\n", nshow, length(selected_buses))
  end
  println(io, "========================================================================================\n")
  return nothing
end

function _print_matpower_reference_diagnostics(io::IO, mpc; matpower_shift_sign::Real = 1.0, matpower_shift_unit = "deg", matpower_ratio = "normal", diagnose_branch_shift_conventions::Bool = false, diagnose_branch_neighborhood::Bool = false, diagnose_residual_clusters::Bool = false, diagnose_residual_cluster_threshold_mw::Float64 = 1.0, diagnose_residual_cluster_maxlines::Int = 20, diagnose_branch_neighborhood_buses::AbstractVector{Int} = Int[], diagnose_branch_neighborhood_depth::Int = 1, diagnose_branch_neighborhood_maxlines::Int = 200, diagnose_nodal_balance_breakdown::Bool = false, diagnose_nodal_balance_buses::AbstractVector{Int} = Int[], diagnose_nodal_balance_maxlines::Int = 200, diagnose_nodal_balance_include_branches::Bool = true, diagnose_nodal_balance_include_generators::Bool = true, diagnose_nodal_balance_include_shunts::Bool = true, diagnose_negative_branch_impedance::Bool = true, diagnose_negative_branch_impedance_maxlines::Int = 100, diagnose_negative_branch_impedance_fail_on_negative_r::Bool = false, diagnose_negative_branch_impedance_fail_on_negative_x::Bool = false, diagnose_negative_branch_impedance_warn_threshold_abs_r::Float64 = 0.0, diagnose_negative_branch_impedance_warn_threshold_abs_x::Float64 = 0.0, matpower_pv_voltage_source = :gen_vg, maxlines::Int = 12)
  if !mp_has_vm_va(mpc)
    if diagnose_negative_branch_impedance || diagnose_branch_neighborhood
      Base.invokelatest(
        getfield(@__MODULE__, :_print_negative_branch_impedance_diagnostics),
        io,
        mpc;
        matpower_shift_sign = matpower_shift_sign,
        matpower_shift_unit = matpower_shift_unit,
        matpower_ratio = matpower_ratio,
        maxlines = diagnose_negative_branch_impedance_maxlines,
        fail_on_negative_r = diagnose_negative_branch_impedance_fail_on_negative_r,
        fail_on_negative_x = diagnose_negative_branch_impedance_fail_on_negative_x,
        warn_threshold_abs_r = diagnose_negative_branch_impedance_warn_threshold_abs_r,
        warn_threshold_abs_x = diagnose_negative_branch_impedance_warn_threshold_abs_x,
        diagnose_branch_neighborhood_buses = diagnose_branch_neighborhood_buses,
        diagnose_nodal_balance_buses = diagnose_nodal_balance_buses,
      )
    end
    return nothing
  end
  base_diag = Base.invokelatest(getfield(@__MODULE__, :_matpower_reference_residuals), mpc; matpower_shift_sign = matpower_shift_sign, matpower_shift_unit = matpower_shift_unit, matpower_ratio = matpower_ratio)
  label = "configured SHIFT sign=$(matpower_shift_sign), unit=$(matpower_shift_unit), ratio=$(matpower_ratio)"
  Base.invokelatest(getfield(@__MODULE__, :_print_top_residuals), io, label, base_diag, mpc.baseMVA; maxlines = maxlines)

  if diagnose_negative_branch_impedance || diagnose_branch_neighborhood
    Base.invokelatest(
      getfield(@__MODULE__, :_print_negative_branch_impedance_diagnostics),
      io,
      mpc;
      matpower_shift_sign = matpower_shift_sign,
      matpower_shift_unit = matpower_shift_unit,
      matpower_ratio = matpower_ratio,
      maxlines = diagnose_negative_branch_impedance_maxlines,
      fail_on_negative_r = diagnose_negative_branch_impedance_fail_on_negative_r,
      fail_on_negative_x = diagnose_negative_branch_impedance_fail_on_negative_x,
      warn_threshold_abs_r = diagnose_negative_branch_impedance_warn_threshold_abs_r,
      warn_threshold_abs_x = diagnose_negative_branch_impedance_warn_threshold_abs_x,
      diagnose_branch_neighborhood_buses = diagnose_branch_neighborhood_buses,
      diagnose_nodal_balance_buses = diagnose_nodal_balance_buses,
      high_residual_buses = _matpower_high_residual_bus_set(base_diag; maxlines = maxlines),
    )
  end

  if diagnose_residual_clusters
    Base.invokelatest(
      getfield(@__MODULE__, :_print_residual_cluster_diagnostics),
      io,
      mpc,
      base_diag;
      threshold_mw = diagnose_residual_cluster_threshold_mw,
      maxlines = diagnose_residual_cluster_maxlines,
      matpower_shift_sign = matpower_shift_sign,
      matpower_shift_unit = matpower_shift_unit,
      matpower_ratio = matpower_ratio,
    )
  end

  if diagnose_branch_neighborhood
    Base.invokelatest(
      getfield(@__MODULE__, :_print_matpower_branch_neighborhood_diagnostics),
      io,
      mpc,
      base_diag;
      buses = diagnose_branch_neighborhood_buses,
      depth = diagnose_branch_neighborhood_depth,
      maxlines = diagnose_branch_neighborhood_maxlines,
      matpower_shift_sign = matpower_shift_sign,
      matpower_shift_unit = matpower_shift_unit,
      matpower_ratio = matpower_ratio,
    )
  end

  if diagnose_nodal_balance_breakdown
    Base.invokelatest(
      getfield(@__MODULE__, :_print_nodal_balance_breakdown),
      io,
      mpc,
      base_diag;
      buses = diagnose_nodal_balance_buses,
      maxlines = diagnose_nodal_balance_maxlines,
      include_branches = diagnose_nodal_balance_include_branches,
      include_generators = diagnose_nodal_balance_include_generators,
      include_shunts = diagnose_nodal_balance_include_shunts,
      matpower_shift_sign = matpower_shift_sign,
      matpower_shift_unit = matpower_shift_unit,
      matpower_ratio = matpower_ratio,
      matpower_pv_voltage_source = matpower_pv_voltage_source,
    )
  end

  max_p_pu = isempty(base_diag.p_rows) ? 0.0 : maximum(abs.(real.(base_diag.mis[base_diag.p_rows])))
  max_q_pu = isempty(base_diag.q_rows) ? 0.0 : maximum(abs.(imag.(base_diag.mis[base_diag.q_rows])))
  if max(max_p_pu, max_q_pu) >= 1.0
    println(io, "The dominant residual is already present in the raw MATPOWER-style fixed-reference check.")
    println(io, "This indicates that the stored VM/VA reference values are not power-balanced for the imported branch data.")
    println(io, "This is not necessarily a Sparlectra solver or import error.")
    println(io)
  end

  if diagnose_branch_shift_conventions
    if size(mpc.branch, 2) >= 10
      shifted = findall(e -> mpc.branch[e, 10] != 0.0, axes(mpc.branch, 1))
      println(io, "Non-zero branch SHIFT entries: ", length(shifted))
      if !isempty(shifted)
        println(io, " first entries:  row   f_bus  t_bus   ratio       shift_raw")
        for e in shifted[1:min(maxlines, length(shifted))]
          @printf(io, "              %5d %7d %6d %8.5f %14.6f\n", e, Int(mpc.branch[e, 1]), Int(mpc.branch[e, 2]), mpc.branch[e, 9], mpc.branch[e, 10])
        end
      end
    end
    println(io, "Branch-shift convention scan (using MATPOWER VM/VA as a fixed reference):")
    println(io, " convention                         max|dP|_MW     max|dQ|_MVAr")
    variants = (
      ("MATPOWER sign, degrees", 1.0, "deg", "normal", true),
      ("opposite sign, degrees", -1.0, "deg", "normal", true),
      ("MATPOWER sign, radians", 1.0, "rad", "normal", true),
      ("opposite sign, radians", -1.0, "rad", "normal", true),
      ("configured, reciprocal ratio", matpower_shift_sign, matpower_shift_unit, "reciprocal", true),
      ("configured, bus shunts disabled", matpower_shift_sign, matpower_shift_unit, matpower_ratio, false),
    )
    for (variant_label, shift_sign, shift_unit, ratio_mode, keep_shunts) in variants
      d = Base.invokelatest(getfield(@__MODULE__, :_matpower_reference_residuals), mpc; matpower_shift_sign = shift_sign, matpower_shift_unit = shift_unit, matpower_ratio = ratio_mode, keep_shunts = keep_shunts)
      max_p = isempty(d.p_rows) ? NaN : maximum(abs.(real.(d.mis[d.p_rows]))) * mpc.baseMVA
      max_q = isempty(d.q_rows) ? NaN : maximum(abs.(imag.(d.mis[d.q_rows]))) * mpc.baseMVA
      @printf(io, " %-34s %13.3f %15.3f\n", variant_label, max_p, max_q)
    end
    println(io)
  end
  return nothing
end

function _print_effective_config(io::IO, cfg; yaml_path::String = "", case_name::String = "", methods = Symbol[])
  println(io, "==================== Effective YAML / MATPOWER import config ====================")
  !isempty(yaml_path) && println(io, "yaml_path: ", yaml_path)
  !isempty(case_name) && println(io, "case: ", case_name)
  !isempty(methods) && println(io, "methods: ", collect(methods))
  for key in sort!(collect(keys(pairs(cfg))); by = String)
    println(io, key, ": ", getproperty(cfg, key))
  end
  println(io, "===============================================================================\n")
  return nothing
end
function _as_auto_profile_mode(v)::Symbol
  v === false && return :off
  v === true && return :apply
  s = lowercase(String(v))
  s in ("false", "off", "none", "no", "0") && return :off
  s == "recommend" && return :recommend
  s == "apply" && return :apply
  @warn "Unknown matpower_auto_profile; using false/off" value = v
  return :off
end

function _namedtuple_from_symbol_dict(d::Dict{Symbol,Any})
  isempty(d) && return (;)
  ordered = sort!(collect(keys(d)); by = String)
  return NamedTuple{Tuple(ordered)}(Tuple(d[k] for k in ordered))
end

function _max_reference_residual_mw(diag, baseMVA::Real)
  max_p = isempty(diag.p_rows) ? 0.0 : maximum(abs.(real.(diag.mis[diag.p_rows]))) * baseMVA
  max_q = isempty(diag.q_rows) ? 0.0 : maximum(abs.(imag.(diag.mis[diag.q_rows]))) * baseMVA
  return (; max_p, max_q, score = max(max_p, max_q))
end

function _matpower_auto_add!(recommendations::Dict{Symbol,Any}, reasons::Dict{Symbol,String}, option::Symbol, value, reason::AbstractString)
  recommendations[option] = value
  reasons[option] = String(reason)
  return nothing
end

function _matpower_auto_profile_shift_scan!(recommendations::Dict{Symbol,Any}, reasons::Dict{Symbol,String}, evidence::Vector{String}, mpc, cfg)
  mp_has_vm_va(mpc) || (push!(evidence, "branch-shift scan skipped: no usable MATPOWER VM/VA reference"); return nothing)
  variants = [
    (; label = "standard sign/degrees + normal ratio", matpower_shift_sign = 1.0, matpower_shift_unit = "deg", matpower_ratio = "normal"),
    (; label = "opposite sign/radians + normal ratio", matpower_shift_sign = -1.0, matpower_shift_unit = "rad", matpower_ratio = "normal"),
    (; label = "standard sign/degrees + reciprocal ratio", matpower_shift_sign = 1.0, matpower_shift_unit = "deg", matpower_ratio = "reciprocal"),
    (; label = "opposite sign/radians + reciprocal ratio", matpower_shift_sign = -1.0, matpower_shift_unit = "rad", matpower_ratio = "reciprocal"),
  ]
  rows = NamedTuple[]
  current = nothing
  for v in variants
    diag = Base.invokelatest(getfield(@__MODULE__, :_matpower_reference_residuals), mpc; matpower_shift_sign = v.matpower_shift_sign, matpower_shift_unit = v.matpower_shift_unit, matpower_ratio = v.matpower_ratio)
    res = _max_reference_residual_mw(diag, mpc.baseMVA)
    row = merge(v, res)
    push!(rows, row)
    if isapprox(v.matpower_shift_sign, cfg.matpower_shift_sign; atol = 0.0, rtol = 0.0) && v.matpower_shift_unit == cfg.matpower_shift_unit && v.matpower_ratio == cfg.matpower_ratio
      current = row
    end
  end
  best = rows[argmin([r.score for r in rows])]
  current = isnothing(current) ? best : current
  push!(evidence, @sprintf("branch-shift scan best='%s' score=%.3f MW/MVAr; current='%s' score=%.3f", best.label, best.score, current.label, current.score))
  if best.score + 1e-9 < 0.5 * max(current.score, 1e-9)
    reason = @sprintf("branch-shift scan strongly improves fixed-reference residual score from %.3f to %.3f MW/MVAr", current.score, best.score)
    _matpower_auto_add!(recommendations, reasons, :matpower_shift_sign, best.matpower_shift_sign, reason)
    _matpower_auto_add!(recommendations, reasons, :matpower_shift_unit, best.matpower_shift_unit, reason)
    _matpower_auto_add!(recommendations, reasons, :matpower_ratio, best.matpower_ratio, reason)
  end
  return nothing
end

function _matpower_auto_profile_shunt_scan!(recommendations::Dict{Symbol,Any}, reasons::Dict{Symbol,String}, evidence::Vector{String}, mpc, cfg)
  mp_has_vm_va(mpc) || (push!(evidence, "bus-shunt scan skipped: no usable MATPOWER VM/VA reference"); return nothing)
  keep = Base.invokelatest(getfield(@__MODULE__, :_matpower_reference_residuals), mpc; matpower_shift_sign = cfg.matpower_shift_sign, matpower_shift_unit = cfg.matpower_shift_unit, matpower_ratio = cfg.matpower_ratio, keep_shunts = true)
  drop = Base.invokelatest(getfield(@__MODULE__, :_matpower_reference_residuals), mpc; matpower_shift_sign = cfg.matpower_shift_sign, matpower_shift_unit = cfg.matpower_shift_unit, matpower_ratio = cfg.matpower_ratio, keep_shunts = false)
  keep_res = _max_reference_residual_mw(keep, mpc.baseMVA)
  drop_res = _max_reference_residual_mw(drop, mpc.baseMVA)
  push!(evidence, @sprintf("bus-shunt scan keep max|dP|=%.3f MW max|dQ|=%.3f MVAr; drop max|dP|=%.3f MW max|dQ|=%.3f MVAr", keep_res.max_p, keep_res.max_q, drop_res.max_p, drop_res.max_q))
  if drop_res.max_q > 1.05 * max(keep_res.max_q, 1e-9) || drop_res.max_p >= 0.95 * max(keep_res.max_p, 1e-9)
    _matpower_auto_add!(recommendations, reasons, :bus_shunt_model, "admittance", "disabling MATPOWER bus shunts does not improve active residuals or worsens reactive residuals")
  end
  return nothing
end

function _matpower_auto_profile_voltage_source!(recommendations::Dict{Symbol,Any}, reasons::Dict{Symbol,String}, evidence::Vector{String}, mpc)
  rows = MatpowerIO.pv_voltage_reference_rows(mpc; warn = false)
  isempty(rows) && (push!(evidence, "PV/REF voltage-source scan found no PV/REF buses"); return nothing)
  with_online = count(r -> !isempty(r.gen_vgs), rows)
  mismatched = count(r -> !isempty(r.gen_vgs) && isfinite(r.dvg_bus) && abs(r.dvg_bus) > 1e-4, rows)
  push!(evidence, "PV/REF voltage-source scan: $(with_online)/$(length(rows)) buses have online GEN.VG targets; $(mismatched) differ from BUS.VM by more than 1e-4 pu")
  if with_online > 0
    _matpower_auto_add!(recommendations, reasons, :matpower_pv_voltage_source, :gen_vg, "PV/REF buses have online generator voltage targets")
    _matpower_auto_add!(recommendations, reasons, :compare_voltage_reference, :hybrid, "hybrid comparison uses active setpoints for PV/REF buses while preserving BUS.VM for PQ and switched buses")
  end
  return nothing
end

function _matpower_auto_profile_flatstart!(recommendations::Dict{Symbol,Any}, reasons::Dict{Symbol,String}, evidence::Vector{String}, mpc)
  nbus = size(mpc.bus, 1)
  nbranch = size(mpc.branch, 1)
  nonzero_shift = size(mpc.branch, 2) >= 10 ? count(!iszero, mpc.branch[:, 10]) : 0
  large_or_shifted = nbus >= 1000 || nbranch >= 2000 || nonzero_shift >= 10
  push!(evidence, "flatstart scan: nbus=$(nbus), nbranch=$(nbranch), nonzero branch SHIFT=$(nonzero_shift)")
  if large_or_shifted
    reason = "large or phase-shifted MATPOWER case benefits from DC-angle and blended-voltage flat-start seeds"
    _matpower_auto_add!(recommendations, reasons, :opt_flatstart, true, reason)
    _matpower_auto_add!(recommendations, reasons, :flatstart_angle_mode, :dc, reason)
    _matpower_auto_add!(recommendations, reasons, :flatstart_voltage_mode, :bus_vm_va_blend, reason)
    _matpower_auto_add!(recommendations, reasons, :start_projection, true, reason)
    _matpower_auto_add!(recommendations, reasons, :start_projection_try_dc_start, true, reason)
    _matpower_auto_add!(recommendations, reasons, :start_projection_try_blend_scan, false, "DC start plus bus VM/VA blend is the compact default pre-run profile for large cases")
    _matpower_auto_add!(recommendations, reasons, :start_projection_measure_candidates, false, "large-case fast path skips candidate mismatch scans when the DC start is explicitly requested")
    _matpower_auto_add!(recommendations, reasons, :start_projection_reuse_import_data, true, "reuse the already parsed MATPOWER case for import and flat-start lookup")
  end
  return nothing
end

function _matpower_auto_profile_qlimits!(recommendations::Dict{Symbol,Any}, reasons::Dict{Symbol,String}, evidence::Vector{String}, mpc)
  size(mpc.gen, 2) >= 8 || (push!(evidence, "Q-limit scan skipped: MATPOWER gen table has no status column"); return nothing)
  active = [g for g in axes(mpc.gen, 1) if mpc.gen[g, 8] > 0.0]
  isempty(active) && (push!(evidence, "Q-limit scan found no online generators"); return nothing)
  narrow = 0
  for g in active
    qmax = Float64(mpc.gen[g, 4])
    qmin = Float64(mpc.gen[g, 5])
    width = qmax - qmin
    if abs(width) <= 1e-6 || max(abs(qmax), abs(qmin)) <= 1e-6 || width <= 5.0
      narrow += 1
    end
  end
  share = narrow / length(active)
  push!(evidence, @sprintf("Q-limit scan: %d/%d online generators have zero or narrow Q range", narrow, length(active)))
  if narrow >= 5 && share >= 0.25
    _matpower_auto_add!(recommendations, reasons, :qlimit_start_iter, 3, "many generators have zero or narrow reactive-power limits; delay switching slightly")
    _matpower_auto_add!(recommendations, reasons, :qlimit_start_mode, :iteration_or_auto, "many narrow Q limits benefit from conservative iteration-or-stability start")
    _matpower_auto_add!(recommendations, reasons, :q_hyst_pu, 0.01, "many narrow Q limits benefit from a small hysteresis deadband")
    _matpower_auto_add!(recommendations, reasons, :cooldown_iters, 1, "many narrow Q limits benefit from a short PV/PQ cooldown")
  end
  return nothing
end


function _print_vmva_self_check(io::IO, mpc; matpower_shift_sign::Real = 1.0, matpower_shift_unit = "deg", matpower_ratio = "normal")
  vmva_chk = MatpowerIO.vmva_power_mismatch_stats(mpc; matpower_shift_sign = matpower_shift_sign, matpower_shift_unit = matpower_shift_unit, matpower_ratio = matpower_ratio)
  if get(vmva_chk, :ok, false)
    println(io, "==================== MATPOWER VM/VA self-check ====================")
    println(io, "SHIFT convention : sign=", matpower_shift_sign, " unit=", matpower_shift_unit, " ratio=", matpower_ratio)
    println(io, "checked eqns     : P(PQ+PV)=", vmva_chk.n_p, "  Q(PQ)=", vmva_chk.n_q)
    println(io, "max |ΔP| (PQ+PV) : ", vmva_chk.max_p_mis_pu, " pu  (", vmva_chk.max_p_mis_MW, " MW)")
    println(io, "max |ΔQ| (PQ)    : ", vmva_chk.max_q_mis_pu, " pu  (", vmva_chk.max_q_mis_MVar, " MVar)")
    println(io, "===================================================================\n")
  else
    println(io, "MATPOWER VM/VA self-check skipped: ", get(vmva_chk, :msg, "unknown reason"))
  end
  return nothing
end

function _mpc_pv_bus_ids(mpc)
  size(mpc.bus, 2) >= 2 || return Int[]
  pv_mask = mpc.bus[:, 2] .== 2
  return sort!(unique!(Int.(round.(mpc.bus[pv_mask, 1]))))
end

function _run_iterations(status)::Int
  return Int(get(status, :iterations, -1))
end

function _run_elapsed_s(status)::Float64
  return Float64(get(status, :elapsed_s, NaN))
end

function _print_converged_loss_summary(io::IO, method::Symbol, status, net::Sparlectra.Net)
  iterations = _run_iterations(status)
  elapsed_s = _run_elapsed_s(status)
  converged = Bool(get(status, :converged, false))
  if converged
    p_loss, q_loss = getTotalLosses(net = net)
    @printf(io, "summary method=%-12s  converged=yes  iterations=%d  time=%8.6f s  losses P=%10.6f MW  Q=%10.6f MVar\n", String(method), iterations, elapsed_s, p_loss, q_loss)
  else
    @printf(io, "summary method=%-12s  converged=no   iterations=%d  time=%8.6f s  losses=SKIP\n", String(method), iterations, elapsed_s)
  end
  return nothing
end

function _show_once_summary_row(method::Symbol, status, stats, cmp_ok::Bool; compare_available::Bool, net = nothing)
  converged = Bool(get(status, :converged, false))
  iterations = _run_iterations(status)
  elapsed_s = _run_elapsed_s(status)
  numerical_solution = Bool(get(status, :numerical_converged, converged)) ? :ok : :fail
  q_limit_active_set = Bool(get(status, :q_limit_active_set_ok, converged)) ? :ok : :fail
  final_converged = Bool(get(status, :final_converged, converged))
  reason_text = String(get(status, :reason_text, final_converged ? "none" : "not converged"))
  status_text = String(get(status, :status, final_converged ? :converged : :not_converged))
  qcounts = _qlimit_summary_counts(net)
  if compare_available
    return (
      method = method,
      converged = converged,
      iterations = iterations,
      elapsed_s = elapsed_s,
      max_dvm = Float64(get(stats, :max_dvm, NaN)),
      max_dva = Float64(get(stats, :max_dva, NaN)),
      slack_delta_va = Float64(get(stats, :slack_delta_va, NaN)),
      angle_alignment = get(stats, :angle_alignment, :none),
      cmp_ok = cmp_ok,
      compare_status = Symbol(get(stats, :compare_status, cmp_ok ? :ok : :fail)),
      numerical_solution = numerical_solution,
      q_limit_active_set = q_limit_active_set,
      final_converged = final_converged,
      reason_text = reason_text,
      status = status_text,
      pv2pq_events = qcounts.pv2pq_events,
      pv2pq_buses = qcounts.pv2pq_buses,
      guarded_pv_buses = qcounts.guarded_pv_buses,
      oscillating_buses = qcounts.oscillating_buses,
    )
  end
  return (
    method = method,
    converged = converged,
    iterations = iterations,
    elapsed_s = elapsed_s,
    max_dvm = NaN,
    max_dva = NaN,
    slack_delta_va = NaN,
    angle_alignment = :none,
    cmp_ok = false,
    compare_status = :skip,
    numerical_solution = numerical_solution,
    q_limit_active_set = q_limit_active_set,
    final_converged = final_converged,
    reason_text = reason_text,
    status = status_text,
    pv2pq_events = qcounts.pv2pq_events,
    pv2pq_buses = qcounts.pv2pq_buses,
    guarded_pv_buses = qcounts.guarded_pv_buses,
    oscillating_buses = qcounts.oscillating_buses,
  )
end

function _print_pv_voltage_reference_diagnostics(io::IO, mpc, net; matpower_pv_voltage_source = :gen_vg, compare_voltage_reference = :bus_vm, tol::Float64 = 1e-4, maxlines::Int = 30)
  rows = MatpowerIO.pv_voltage_reference_rows(mpc; net = net, matpower_pv_voltage_source = matpower_pv_voltage_source, tol = tol, warn = true)
  qmin_pu, qmax_pu = getQLimits_pu(net)
  enriched_rows = NamedTuple[]
  for row in rows
    node_idx = findfirst(k -> (haskey(net.busOrigIdxDict, k) ? net.busOrigIdxDict[k] : k) == row.busI, eachindex(net.nodeVec))
    if node_idx === nothing
      continue
    end
    node = net.nodeVec[node_idx]
    final_type = getNodeType(node)
    switched_to_pq = haskey(net.qLimitEvents, node_idx)
    active_regulating = (final_type == Sparlectra.Slack || final_type == Sparlectra.PV) && !switched_to_pq
    qcalc = node._qƩGen === nothing ? NaN : Float64(node._qƩGen)
    push!(
      enriched_rows,
      (
        row...,
        node_idx = node_idx,
        final_type = final_type == Sparlectra.Slack ? "REF" : final_type == Sparlectra.PV ? "PV" : final_type == Sparlectra.PQ ? "PQ" : string(final_type),
        switched_to_pq = switched_to_pq,
        active_regulating = active_regulating,
        qcalc = qcalc,
        qmin = qmin_pu[node_idx] * net.baseMVA,
        qmax = qmax_pu[node_idx] * net.baseMVA,
      ),
    )
  end
  println(io, "==================== PV voltage reference diagnostics ====================")
  println(io, "source option: matpower_pv_voltage_source = ", matpower_pv_voltage_source)
  println(io, "compare option: compare_voltage_reference = ", compare_voltage_reference)
  println(io, "tolerance: ", tol, " pu")
  if isempty(rows)
    println(io, "No MATPOWER PV/REF buses found.")
    println(io, "==========================================================================")
    return nothing
  end
  ord_mis = sortperm(1:length(rows); by = i -> (isfinite(rows[i].dvg_bus) ? 0 : 1, isfinite(rows[i].dvg_bus) ? -abs(rows[i].dvg_bus) : 0.0))
  println(io, "\nTop PV/REF BUS.VM vs GEN.VG mismatches (online GEN.VG rows first):")
  println(io, " rank  BUS_I  type  BUS.VM   GEN.VG_values        imported_Vset  dVG_BUS")
  for rank = 1:min(maxlines, length(ord_mis))
    row = rows[ord_mis[rank]]
    @printf(io, " %4d  %5d   %3s  %7.5f  %-20s %11.5f  %+8.5f\n", rank, row.busI, row.bus_type == 3 ? "REF" : "PV", row.bus_vm, string(row.gen_vgs), row.imported_vset, row.dvg_bus)
  end
  ord_final = sortperm(1:length(rows); by = i -> max(abs(rows[i].dvm_bus), abs(rows[i].dvm_vset)), rev = true)
  println(io, "\nTop final Vm deviations explained by selected setpoint:")
  println(io, " rank  BUS_I  type  BUS.VM   GEN.VG  imported_Vset  Vm_calc  dVm_BUS  dVm_Vset")
  for rank = 1:min(maxlines, length(ord_final))
    row = rows[ord_final[rank]]
    gen_vg = isempty(row.gen_vgs) ? NaN : row.gen_vgs[1]
    @printf(io, " %4d  %5d   %3s  %7.5f  %7.5f  %12.5f  %7.5f  %+8.5f  %+9.5f\n", rank, row.busI, row.bus_type == 3 ? "REF" : "PV", row.bus_vm, gen_vg, row.imported_vset, row.vm_calc, row.dvm_bus, row.dvm_vset)
  end
  active_rows = filter(row -> row.active_regulating, enriched_rows)
  max_vset_residual = isempty(active_rows) ? 0.0 : maximum(abs(row.dvm_vset) for row in active_rows)
  @printf(io, "\nPost-solve active PV/REF setpoint residual: max |Vm_calc - imported_Vset| = %.8g pu over %d buses\n", max_vset_residual, length(active_rows))
  ord_active = sortperm(1:length(active_rows); by = i -> abs(active_rows[i].dvm_vset), rev = true)
  println(io, "Worst active PV/REF buses by imported setpoint residual:")
  println(io, " rank  BUS_I  BUS.VM   GEN.VG  imported_Vset  Vm_calc  dVm_BUS  dVm_Vset  final  Qcalc     Qmin      Qmax")
  for rank = 1:min(maxlines, length(ord_active))
    row = active_rows[ord_active[rank]]
    gen_vg = isempty(row.gen_vgs) ? NaN : row.gen_vgs[1]
    @printf(io, " %4d  %5d  %7.5f  %7.5f  %12.5f  %7.5f  %+8.5f  %+9.5f  %-5s  %8.3f  %8.3f  %8.3f\n", rank, row.busI, row.bus_vm, gen_vg, row.imported_vset, row.vm_calc, row.dvm_bus, row.dvm_vset, row.final_type, row.qcalc, row.qmin, row.qmax)
  end
  fallback_count = count(row -> row.imported_kind == :BUS_VM_FALLBACK, enriched_rows)
  println(io, "PV/REF buses without online generators using BUS.VM fallback: ", fallback_count)
  println(io, "\nInterpretation:")
  println(io, "- Large dVm_BUS but small dVm_Vset means Sparlectra is following the imported PV setpoint.")
  println(io, "- dVm_BUS is the deviation from MATPOWER BUS.VM; dVm_Vset is the deviation from imported GEN.VG/imported_Vset.")
  println(io, "- BUS.VM fallback rows are PV/REF buses without an online generator voltage setpoint.")
  println(io, "- If compare_voltage_reference=bus_vm, this can produce a compare FAIL even when the solver is correct.")
  println(io, "==========================================================================")
  return nothing
end

function _wrong_branch_stats(net, mpc; compare_voltage_reference = :bus_vm, matpower_pv_voltage_source = :gen_vg, tol::Float64 = 1e-4)
  vm = filter(isfinite, [n._vm_pu === nothing ? NaN : Float64(n._vm_pu) for n in net.nodeVec if !isIsolated(n)])
  cmp_ok, cmp_stats = MatpowerIO.compare_vm_va(net, mpc; show_diff = false, tol_vm = 1e9, tol_va = 1e9, compare_voltage_reference = compare_voltage_reference, matpower_pv_voltage_source = matpower_pv_voltage_source, matpower_pv_voltage_mismatch_tol_pu = tol)
  return (min_vm = isempty(vm) ? NaN : minimum(vm), max_vm = isempty(vm) ? NaN : maximum(vm), max_dva = Float64(get(cmp_stats, :max_dva, NaN)), max_dvm = Float64(get(cmp_stats, :max_dvm, NaN)))
end

function _print_wrong_branch_warning(io::IO, net, mpc; wrong_branch_detection::Bool, wrong_branch_min_vm_pu::Float64, wrong_branch_max_angle_spread_deg::Float64, compare_voltage_reference = :bus_vm, matpower_pv_voltage_source = :gen_vg, tol::Float64 = 1e-4)
  wrong_branch_detection || return nothing
  stats = _wrong_branch_stats(net, mpc; compare_voltage_reference = compare_voltage_reference, matpower_pv_voltage_source = matpower_pv_voltage_source, tol = tol)
  suspicious = stats.min_vm < wrong_branch_min_vm_pu || stats.max_dva > wrong_branch_max_angle_spread_deg
  if suspicious
    println(io, "WARNING: Solver converged to a suspicious low-voltage / wrong-branch solution.")
    @printf(io, "  minVm=%.5f pu, max|dVm|=%.5f pu, max|dVa|=%.5f deg\n", stats.min_vm, stats.max_dvm, stats.max_dva)
    println(io, "  Consider flatstart_angle_mode=dc or flatstart_voltage_mode=bus_vm_va_blend.")
  end
  return suspicious
end

function _print_dataframe_nodes(io::IO, net::Sparlectra.Net; max_nodes::Int = 0)
  bus_name_by_idx = Dict{Int,String}()
  for (name, idx) in net.busDict
    bus_name_by_idx[idx] = name
  end
  nodes = sort(net.nodeVec, by = x -> x.busIdx)
  n_total = length(nodes)
  n_show = max_nodes > 0 ? min(max_nodes, n_total) : n_total

  println(io, "DataFrame-style node output:")
  println(io, " row  bus  bus_name              type     vm_pu     va_deg")
  println(io, "-----------------------------------------------------------")

  shown = 0
  for i in eachindex(nodes)
    shown >= n_show && break
    n = nodes[i]
    bname = get(bus_name_by_idx, n.busIdx, string(n.busIdx))
    vm = isnothing(n._vm_pu) ? NaN : Float64(n._vm_pu)
    va = isnothing(n._va_deg) ? NaN : Float64(n._va_deg)
    @printf(io, " %4d %4d  %-20s %-6s %8.5f %10.4f\n", i, n.busIdx, bname, String(toString(n._nodeType)), vm, va)
    shown += 1
  end

  if n_show < n_total
    @printf(io, "... truncated: showing %d/%d nodes (show_once_max_nodes)\n", n_show, n_total)
  else
    @printf(io, "shown nodes: %d/%d\n", n_show, n_total)
  end
  println(io)
end

function _enable_pq_gen_controllers_for_method(method::Symbol, requested::Bool)::Bool
  # PQ-generator voltage-dependent controllers are supported only by the rectangular solver.
  return requested && method === :rectangular
end

function _matpower_slack_reference(mpc)
  if isnothing(mpc) || !hasproperty(mpc, :bus) || size(mpc.bus, 2) < 9
    return (bus = missing, vm_pu = missing, va_deg = missing)
  end
  slack_rows = findall(r -> Int(round(mpc.bus[r, 2])) == 3, axes(mpc.bus, 1))
  isempty(slack_rows) && return (bus = missing, vm_pu = missing, va_deg = missing)
  row = first(slack_rows)
  return (bus = Int(round(mpc.bus[row, 1])), vm_pu = Float64(mpc.bus[row, 8]), va_deg = Float64(mpc.bus[row, 9]))
end

function _print_reference_override_status(io::IO, mpc; reference_override::Bool, reference_vm_pu::Float64, reference_va_deg::Float64)
  slack_ref = _matpower_slack_reference(mpc)
  println(io, "MATPOWER slack reference bus: ", slack_ref.bus)
  println(io, "MATPOWER slack reference read: ", slack_ref.vm_pu, " pu / ", slack_ref.va_deg, " deg")
  if reference_override
    println(io, "MATPOWER slack reference override: enabled -> ", reference_vm_pu, " pu / ", reference_va_deg, " deg")
  else
    println(io, "MATPOWER slack reference override: disabled; using imported MATPOWER slack reference")
    if reference_vm_pu != 1.0 || reference_va_deg != 0.0
      println(io, "MATPOWER slack reference override values ignored while disabled: ", reference_vm_pu, " pu / ", reference_va_deg, " deg")
    end
  end
  return nothing
end
# -----------------------------------------------------------------------------
# Benchmark helper: benchmark exactly run_acpflow(...)
# -----------------------------------------------------------------------------
function bench_run_acpflow(;
  casefile::String,
  methods::Vector{Symbol},
  mpc,
  logfile::String,
  show_diff::Bool,
  tol_vm::Float64,
  tol_va::Float64,
  max_ite::Int = 30,
  tol::Float64 = 1e-6,
  opt_fd::Bool = true,
  opt_sparse::Bool = true,
  opt_flatstart::Bool = true,
  autodamp::Bool = false,
  autodamp_min::Float64 = 1e-3,
  start_projection::Bool = false,
  start_projection_try_dc_start::Bool = true,
  start_projection_try_blend_scan::Bool = true,
  start_projection_branch_guard::Bool = true,
  start_projection_measure_candidates::Bool = true,
  start_projection_reuse_import_data::Bool = true,
  start_projection_blend_lambdas::AbstractVector{<:Real} = [0.25, 0.5, 0.75],
  start_projection_dc_angle_limit_deg::Float64 = 60.0,
  qlimit_start_iter::Int = 2,
  qlimit_start_mode::Symbol = :iteration,
  qlimit_auto_q_delta_pu::Float64 = 1e-4,
  qlimit_guard::Bool = false,
  qlimit_guard_min_q_range_pu::Float64 = 1e-4,
  qlimit_guard_zero_range_mode::Symbol = :lock_pq,
  qlimit_guard_narrow_range_mode::Symbol = :prefer_pq,
  qlimit_guard_max_switches::Int = 10,
  qlimit_guard_freeze_after_repeated_switching::Bool = true,
  qlimit_guard_accept_bounded_violations::Bool = false,
  qlimit_guard_max_remaining_violations::Int = 0,
  qlimit_guard_violation_mode::Symbol = :delayed_switch,
  qlimit_guard_violation_threshold_pu::Float64 = 1e-4,
  qlimit_guard_log::Bool = true,
  verbose::Int = 0,
  cooldown_iters::Int = 0,
  q_hyst_pu::Float64 = 0.0,
  qlimit_trace_buses::AbstractVector{Int} = Int[],
  qlimit_lock_reason::Symbol = :manual,
  lock_pv_to_pq_buses::AbstractVector{Int} = Int[],
  seconds::Float64 = 2.0,
  samples::Int = 50,
  show_once::Bool = false,
  show_once_output::Symbol = :classic,
  show_once_max_nodes::Int = 0,
  benchmark::Bool = true,
  matpower_shift_sign::Float64 = 1.0,
  matpower_shift_unit::String = "deg",
  matpower_ratio::String = "normal",
  reference_override::Bool = false,
  reference_vm_pu::Float64 = 1.0,
  reference_va_deg::Float64 = 0.0,
  diagnose_matpower_reference::Bool = false,
  diagnose_branch_shift_conventions::Bool = false,
  diagnose_branch_neighborhood::Bool = false,
  diagnose_residual_clusters::Bool = false,
  diagnose_residual_cluster_threshold_mw::Float64 = 1.0,
  diagnose_residual_cluster_maxlines::Int = 20,
  diagnose_branch_neighborhood_buses::AbstractVector{Int} = Int[],
  diagnose_branch_neighborhood_depth::Int = 1,
  diagnose_branch_neighborhood_maxlines::Int = 200,
  diagnose_nodal_balance_breakdown::Bool = false,
  diagnose_nodal_balance_buses::AbstractVector{Int} = Int[],
  diagnose_nodal_balance_maxlines::Int = 200,
  diagnose_nodal_balance_include_branches::Bool = true,
  diagnose_nodal_balance_include_generators::Bool = true,
  diagnose_nodal_balance_include_shunts::Bool = true,
  diagnose_negative_branch_impedance::Bool = true,
  diagnose_negative_branch_impedance_maxlines::Int = 100,
  diagnose_negative_branch_impedance_fail_on_negative_r::Bool = false,
  diagnose_negative_branch_impedance_fail_on_negative_x::Bool = false,
  diagnose_negative_branch_impedance_warn_threshold_abs_r::Float64 = 0.0,
  diagnose_negative_branch_impedance_warn_threshold_abs_x::Float64 = 0.0,
  diagnose_maxlines::Int = 12,
  log_effective_config::Bool = false,
  yaml_path::String = "",
  effective_config = nothing,
  enable_pq_gen_controllers::Bool = true,
  bus_shunt_model::String = "admittance",
  matpower_pv_voltage_source::Symbol = :gen_vg,
  matpower_pv_voltage_mismatch_tol_pu::Float64 = 1e-4,
  compare_voltage_reference::Symbol = :bus_vm,
  diagnose_pv_voltage_references::Bool = false,
  diagnose_pv_voltage_maxlines::Int = 30,
  flatstart_voltage_mode::Symbol = :classic,
  flatstart_angle_mode::Symbol = :classic,
  wrong_branch_detection::Bool = true,
  wrong_branch_min_vm_pu::Float64 = 0.4,
  wrong_branch_max_angle_spread_deg::Float64 = 120.0,
  wrong_branch_rescue::Bool = false,
  wrong_branch_rescue_modes::AbstractVector{Symbol} = [:dc, :bus_vm_va_blend, :matpower_va],
  console_summary::Bool = true,
  console_auto_profile::Symbol = :compact,
  console_diagnostics::Symbol = :compact,
  console_q_limit_events::Symbol = :summary,
  console_max_rows::Int = 20,
  logfile_diagnostics::Symbol = :full,
  performance_profile = nothing,
  performance_level::Symbol = :summary,
  performance_print_to_console::Bool = true,
  performance_write_to_logfile::Bool = true,
  performance_show_iteration_table::Bool = true,
  performance_skip_reference_comparison::Bool = false,
  performance_skip_expensive_diagnostics::Bool = false,
  performance_skip_branch_neighborhood_report::Bool = false,
  performance_max_diagnostic_rows::Int = 20,
  log_status::_MatpowerImportLogStatus = _MatpowerImportLogStatus(),
)
  t0 = time()
  results = Dict{Symbol,Any}()
  effective_reference_vm_pu = reference_override ? reference_vm_pu : nothing
  effective_reference_va_deg = reference_override ? reference_va_deg : nothing
  effective_diagnose_matpower_reference = diagnose_matpower_reference && !performance_skip_expensive_diagnostics
  effective_diagnose_pv_voltage_references = diagnose_pv_voltage_references && !performance_skip_expensive_diagnostics
  effective_diagnose_branch_neighborhood = diagnose_branch_neighborhood && !performance_skip_branch_neighborhood_report && !performance_skip_expensive_diagnostics
  effective_diagnose_residual_clusters = diagnose_residual_clusters && !performance_skip_expensive_diagnostics
  effective_diagnose_negative_branch_impedance = diagnose_negative_branch_impedance && !performance_skip_expensive_diagnostics
  effective_diagnose_maxlines = min(diagnose_maxlines, performance_max_diagnostic_rows)


  # Append to the logfile because main() may already have captured import warnings there.
  open(logfile, "a") do io
    println(io, "Sparlectra version: ", Sparlectra.version())
    println(io, "casefile: ", casefile)
    println(io, "timestamp: ", Dates.now())
    Base.invokelatest(getfield(@__MODULE__, :_print_reference_override_status), io, mpc; reference_override = reference_override, reference_vm_pu = reference_vm_pu, reference_va_deg = reference_va_deg)
    println(io)
    if log_effective_config && !isnothing(effective_config)
      Base.invokelatest(getfield(@__MODULE__, :_print_effective_config), io, effective_config; yaml_path = yaml_path, case_name = casefile, methods = methods)
    end
    if !isnothing(mpc)
      Base.invokelatest(getfield(@__MODULE__, :_print_vmva_self_check), io, mpc; matpower_shift_sign = matpower_shift_sign, matpower_shift_unit = matpower_shift_unit, matpower_ratio = matpower_ratio)
      if effective_diagnose_matpower_reference
        Base.invokelatest(
          getfield(@__MODULE__, :_print_matpower_reference_diagnostics),
          io,
          mpc;
          matpower_shift_sign = matpower_shift_sign,
          matpower_shift_unit = matpower_shift_unit,
          matpower_ratio = matpower_ratio,
          diagnose_branch_shift_conventions = diagnose_branch_shift_conventions,
          diagnose_branch_neighborhood = effective_diagnose_branch_neighborhood,
          diagnose_residual_clusters = effective_diagnose_residual_clusters,
          diagnose_residual_cluster_threshold_mw = diagnose_residual_cluster_threshold_mw,
          diagnose_residual_cluster_maxlines = diagnose_residual_cluster_maxlines,
          diagnose_branch_neighborhood_buses = diagnose_branch_neighborhood_buses,
          diagnose_branch_neighborhood_depth = diagnose_branch_neighborhood_depth,
          diagnose_branch_neighborhood_maxlines = diagnose_branch_neighborhood_maxlines,
          diagnose_nodal_balance_breakdown = diagnose_nodal_balance_breakdown,
          diagnose_nodal_balance_buses = diagnose_nodal_balance_buses,
          diagnose_nodal_balance_maxlines = diagnose_nodal_balance_maxlines,
          diagnose_nodal_balance_include_branches = diagnose_nodal_balance_include_branches,
          diagnose_nodal_balance_include_generators = diagnose_nodal_balance_include_generators,
          diagnose_nodal_balance_include_shunts = diagnose_nodal_balance_include_shunts,
          diagnose_negative_branch_impedance = effective_diagnose_negative_branch_impedance,
          diagnose_negative_branch_impedance_maxlines = diagnose_negative_branch_impedance_maxlines,
          diagnose_negative_branch_impedance_fail_on_negative_r = diagnose_negative_branch_impedance_fail_on_negative_r,
          diagnose_negative_branch_impedance_fail_on_negative_x = diagnose_negative_branch_impedance_fail_on_negative_x,
          diagnose_negative_branch_impedance_warn_threshold_abs_r = diagnose_negative_branch_impedance_warn_threshold_abs_r,
          diagnose_negative_branch_impedance_warn_threshold_abs_x = diagnose_negative_branch_impedance_warn_threshold_abs_x,
          matpower_pv_voltage_source = matpower_pv_voltage_source,
          maxlines = effective_diagnose_maxlines,
        )
      elseif effective_diagnose_negative_branch_impedance || effective_diagnose_branch_neighborhood
        Base.invokelatest(
          getfield(@__MODULE__, :_print_negative_branch_impedance_diagnostics),
          io,
          mpc;
          matpower_shift_sign = matpower_shift_sign,
          matpower_shift_unit = matpower_shift_unit,
          matpower_ratio = matpower_ratio,
          maxlines = diagnose_negative_branch_impedance_maxlines,
          fail_on_negative_r = diagnose_negative_branch_impedance_fail_on_negative_r,
          fail_on_negative_x = diagnose_negative_branch_impedance_fail_on_negative_x,
          warn_threshold_abs_r = diagnose_negative_branch_impedance_warn_threshold_abs_r,
          warn_threshold_abs_x = diagnose_negative_branch_impedance_warn_threshold_abs_x,
          diagnose_branch_neighborhood_buses = diagnose_branch_neighborhood_buses,
          diagnose_nodal_balance_buses = diagnose_nodal_balance_buses,
        )
      end
    end
  end

  # Optional: show results once (not benchmarked)
  if show_once
    local summaries = NamedTuple[]
    show_classic = (show_once_output == :classic)
    for m in methods
      open(logfile, "a") do io
        redirect_stdout(io) do
          redirect_stderr(io) do
            println("=================================================================")
            println("RUN method = ", m)
            println("=================================================================\n")
            status_ref = Ref{Any}(nothing)
            net_res = _with_log_table(logfile, log_status) do
              run_acpflow(
              casefile = casefile,
              max_ite = max_ite,
              tol = tol,
              opt_fd = opt_fd,
              opt_sparse = opt_sparse,
              method = m,
              autodamp = autodamp,
              autodamp_min = autodamp_min,
              start_projection = start_projection,
              start_projection_try_dc_start = start_projection_try_dc_start,
              start_projection_try_blend_scan = start_projection_try_blend_scan,
              start_projection_branch_guard = start_projection_branch_guard,
              start_projection_measure_candidates = start_projection_measure_candidates,
              start_projection_blend_lambdas = start_projection_blend_lambdas,
              start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
              qlimit_start_iter = qlimit_start_iter,
              qlimit_start_mode = qlimit_start_mode,
              qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu,
              qlimit_guard = qlimit_guard,
              qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu,
              qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode,
              qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode,
              qlimit_guard_max_switches = qlimit_guard_max_switches,
              qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching,
              qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations,
              qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations,
              qlimit_guard_violation_mode = qlimit_guard_violation_mode,
              qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu,
              qlimit_guard_log = qlimit_guard_log,
              pv_table_rows = console_max_rows,
              opt_flatstart = opt_flatstart,
              show_results = show_classic,
              show_compact_result = true,
              status_ref = status_ref,
              verbose = verbose,
              cooldown_iters = cooldown_iters,
              q_hyst_pu = q_hyst_pu,
              qlimit_trace_buses = qlimit_trace_buses,
              qlimit_lock_reason = qlimit_lock_reason,
              lock_pv_to_pq_buses = lock_pv_to_pq_buses,
              enable_pq_gen_controllers = _enable_pq_gen_controllers_for_method(m, enable_pq_gen_controllers),
              bus_shunt_model = bus_shunt_model,
              matpower_shift_sign = matpower_shift_sign,
              matpower_shift_unit = matpower_shift_unit,
              matpower_ratio = matpower_ratio,
              reference_vm_pu = effective_reference_vm_pu,
              reference_va_deg = effective_reference_va_deg,
              matpower_pv_voltage_source = matpower_pv_voltage_source,
              matpower_pv_voltage_mismatch_tol_pu = matpower_pv_voltage_mismatch_tol_pu,
              flatstart_voltage_mode = flatstart_voltage_mode,
              flatstart_angle_mode = flatstart_angle_mode,
                imported_matpower_case = start_projection_reuse_import_data ? mpc : nothing,
                performance_profile = performance_profile,
              )
            end
            status = status_ref[]
            _print_converged_loss_summary(io, m, status, net_res)
            if !show_classic
              _print_dataframe_nodes(io, net_res; max_nodes = show_once_max_nodes)
            end
            if mp_has_vm_va(mpc) && !performance_skip_reference_comparison
              if effective_diagnose_pv_voltage_references
                Base.invokelatest(getfield(@__MODULE__, :_print_pv_voltage_reference_diagnostics), io, mpc, net_res; matpower_pv_voltage_source = matpower_pv_voltage_source, compare_voltage_reference = compare_voltage_reference, tol = matpower_pv_voltage_mismatch_tol_pu, maxlines = diagnose_pv_voltage_maxlines)
              end
              ok, stats = _perf_profile_time!(performance_profile, :reference_comparison) do
                MatpowerIO.compare_vm_va(
                  net_res,
                  mpc;
                  show_diff = show_diff,
                  tol_vm = tol_vm,
                  tol_va = tol_va,
                  maxlines = 20,
                  compare_voltage_reference = compare_voltage_reference,
                  matpower_pv_voltage_source = matpower_pv_voltage_source,
                  matpower_pv_voltage_mismatch_tol_pu = matpower_pv_voltage_mismatch_tol_pu,
                )
              end
              _print_wrong_branch_warning(
                io,
                net_res,
                mpc;
                wrong_branch_detection = wrong_branch_detection,
                wrong_branch_min_vm_pu = wrong_branch_min_vm_pu,
                wrong_branch_max_angle_spread_deg = wrong_branch_max_angle_spread_deg,
                compare_voltage_reference = compare_voltage_reference,
                matpower_pv_voltage_source = matpower_pv_voltage_source,
                tol = matpower_pv_voltage_mismatch_tol_pu,
              )
              push!(summaries, _show_once_summary_row(m, status, stats, ok; compare_available = true, net = net_res))
            else
  push!(summaries, _show_once_summary_row(m, status, nothing, false; compare_available = false, net = net_res))
  if performance_skip_reference_comparison
    println("Compare skipped: performance.skip_reference_comparison=true")
  else
    println("Compare skipped: no solution-like VM/VA in mpc.bus(:,8:9)")
  end
end
            if !Bool(get(status, :converged, false)) && !enable_pq_gen_controllers && m === :rectangular
              println("Fallback diagnostic: rerun with enable_pq_gen_controllers=true")
              fb_ref = Ref{Any}(nothing)
              _with_log_table(logfile, log_status) do
                run_acpflow(
                  casefile = casefile,
                max_ite = max_ite,
                tol = tol,
                opt_fd = opt_fd,
                opt_sparse = opt_sparse,
                method = m,
                autodamp = autodamp,
                autodamp_min = autodamp_min,
                start_projection = start_projection,
                start_projection_try_dc_start = start_projection_try_dc_start,
                start_projection_try_blend_scan = start_projection_try_blend_scan,
                start_projection_branch_guard = start_projection_branch_guard,
                start_projection_measure_candidates = start_projection_measure_candidates,
                start_projection_blend_lambdas = start_projection_blend_lambdas,
                start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
                qlimit_start_iter = qlimit_start_iter,
                qlimit_start_mode = qlimit_start_mode,
                qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu,
                qlimit_guard = qlimit_guard,
                qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu,
                qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode,
                qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode,
                qlimit_guard_max_switches = qlimit_guard_max_switches,
                qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching,
                qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations,
                qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations,
                qlimit_guard_violation_mode = qlimit_guard_violation_mode,
                qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu,
                qlimit_guard_log = qlimit_guard_log,
                pv_table_rows = console_max_rows,
                opt_flatstart = opt_flatstart,
                show_results = false,
                show_compact_result = true,
                status_ref = fb_ref,
                verbose = 0,
                cooldown_iters = cooldown_iters,
                q_hyst_pu = q_hyst_pu,
                qlimit_trace_buses = qlimit_trace_buses,
                qlimit_lock_reason = qlimit_lock_reason,
                lock_pv_to_pq_buses = lock_pv_to_pq_buses,
                enable_pq_gen_controllers = true,
                bus_shunt_model = bus_shunt_model,
                matpower_shift_sign = matpower_shift_sign,
                matpower_shift_unit = matpower_shift_unit,
                matpower_ratio = matpower_ratio,
                reference_vm_pu = effective_reference_vm_pu,
                reference_va_deg = effective_reference_va_deg,
                matpower_pv_voltage_source = matpower_pv_voltage_source,
                matpower_pv_voltage_mismatch_tol_pu = matpower_pv_voltage_mismatch_tol_pu,
                flatstart_voltage_mode = flatstart_voltage_mode,
                  flatstart_angle_mode = flatstart_angle_mode,
                imported_matpower_case = start_projection_reuse_import_data ? mpc : nothing,
                performance_profile = performance_profile,
                )
              end
            end
            @printf(io, "show_once method=%s output=%s\n\n", String(m), String(show_once_output))
          end
        end
      end
    end

    # print compact operational summary to terminal
    if console_summary
      println()
      _print_matpower_run_summary(stdout, summaries; logfile = logfile, q_limit_mode = console_q_limit_events)
      println()
    end
  end

  if !show_once
    println("\nInitial run (convergence + iterations):")
    for m in methods
      status_ref = Ref{Any}(nothing)
      net_res = _with_log_table(logfile, log_status) do
        run_acpflow(
        casefile = casefile,
        max_ite = max_ite,
        tol = tol,
        opt_fd = opt_fd,
        opt_sparse = opt_sparse,
        method = m,
        autodamp = autodamp,
        autodamp_min = autodamp_min,
        start_projection = start_projection,
        start_projection_try_dc_start = start_projection_try_dc_start,
        start_projection_try_blend_scan = start_projection_try_blend_scan,
        start_projection_branch_guard = start_projection_branch_guard,
        start_projection_measure_candidates = start_projection_measure_candidates,
        start_projection_blend_lambdas = start_projection_blend_lambdas,
        start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
        qlimit_start_iter = qlimit_start_iter,
        qlimit_start_mode = qlimit_start_mode,
        qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu,
        qlimit_guard = qlimit_guard,
        qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu,
        qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode,
        qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode,
        qlimit_guard_max_switches = qlimit_guard_max_switches,
        qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching,
        qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations,
        qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations,
        qlimit_guard_violation_mode = qlimit_guard_violation_mode,
        qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu,
        qlimit_guard_log = qlimit_guard_log,
        pv_table_rows = console_max_rows,
        opt_flatstart = opt_flatstart,
        show_results = false,
        show_compact_result = true,
        status_ref = status_ref,
        verbose = 0,
        cooldown_iters = cooldown_iters,
        q_hyst_pu = q_hyst_pu,
        qlimit_trace_buses = qlimit_trace_buses,
        qlimit_lock_reason = qlimit_lock_reason,
        lock_pv_to_pq_buses = lock_pv_to_pq_buses,
        enable_pq_gen_controllers = _enable_pq_gen_controllers_for_method(m, enable_pq_gen_controllers),
        bus_shunt_model = bus_shunt_model,
        matpower_shift_sign = matpower_shift_sign,
        matpower_shift_unit = matpower_shift_unit,
        matpower_ratio = matpower_ratio,
        reference_vm_pu = effective_reference_vm_pu,
        reference_va_deg = effective_reference_va_deg,
        matpower_pv_voltage_source = matpower_pv_voltage_source,
        matpower_pv_voltage_mismatch_tol_pu = matpower_pv_voltage_mismatch_tol_pu,
        flatstart_voltage_mode = flatstart_voltage_mode,
        flatstart_angle_mode = flatstart_angle_mode,
                imported_matpower_case = start_projection_reuse_import_data ? mpc : nothing,
                performance_profile = performance_profile,
        )
      end
      _print_converged_loss_summary(stdout, m, status_ref[], net_res)
      open(logfile, "a") do io
        _print_converged_loss_summary(io, m, status_ref[], net_res)
      end
    end
  end

  if !benchmark
    total_s = time() - t0
    @printf("total runtime        : %.3f s\n", total_s)
    _perf_profile_time!(performance_profile, :logging_diagnostics) do
      open(logfile, "a") do io
        @printf(io, "total runtime: %.3f s\n", total_s)
        _print_log_status(io, log_status)
      end
      _print_log_status(stdout, log_status)
    end
    _emit_performance_summary(performance_profile; logfile = logfile, print_to_console = performance_print_to_console, write_to_logfile = performance_write_to_logfile, max_rows = performance_max_diagnostic_rows)
    return results
  end

  # Warmup (compile) once per method with minimal output
  println("Warmup run:")
  for m in methods
    _with_log_table(logfile, log_status) do
      run_acpflow(
        casefile = casefile,
      max_ite = max_ite,
      tol = tol,
      opt_fd = opt_fd,
      opt_sparse = opt_sparse,
      method = m,
      autodamp = autodamp,
      autodamp_min = autodamp_min,
      start_projection = start_projection,
      start_projection_try_dc_start = start_projection_try_dc_start,
      start_projection_try_blend_scan = start_projection_try_blend_scan,
      start_projection_branch_guard = start_projection_branch_guard,
      start_projection_measure_candidates = start_projection_measure_candidates,
      start_projection_blend_lambdas = start_projection_blend_lambdas,
      start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
      qlimit_start_iter = qlimit_start_iter,
      qlimit_start_mode = qlimit_start_mode,
      qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu,
      qlimit_guard = qlimit_guard,
      qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu,
      qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode,
      qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode,
      qlimit_guard_max_switches = qlimit_guard_max_switches,
      qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching,
      qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations,
      qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations,
      qlimit_guard_violation_mode = qlimit_guard_violation_mode,
      qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu,
      qlimit_guard_log = qlimit_guard_log,
      pv_table_rows = console_max_rows,
      opt_flatstart = opt_flatstart,
      show_results = false,
      verbose = 0,
      cooldown_iters = cooldown_iters,
      q_hyst_pu = q_hyst_pu,
      qlimit_trace_buses = qlimit_trace_buses,
      qlimit_lock_reason = qlimit_lock_reason,
      lock_pv_to_pq_buses = lock_pv_to_pq_buses,
      enable_pq_gen_controllers = _enable_pq_gen_controllers_for_method(m, enable_pq_gen_controllers),
      bus_shunt_model = bus_shunt_model,
      matpower_shift_sign = matpower_shift_sign,
      matpower_shift_unit = matpower_shift_unit,
      matpower_ratio = matpower_ratio,
      reference_vm_pu = effective_reference_vm_pu,
      reference_va_deg = effective_reference_va_deg,
      matpower_pv_voltage_source = matpower_pv_voltage_source,
      matpower_pv_voltage_mismatch_tol_pu = matpower_pv_voltage_mismatch_tol_pu,
      flatstart_voltage_mode = flatstart_voltage_mode,
      flatstart_angle_mode = flatstart_angle_mode,
                imported_matpower_case = start_projection_reuse_import_data ? mpc : nothing,
                performance_profile = performance_profile,
    )
  end
  end

  println("\n==================== Benchmark run_acpflow ====================")
  println("casefile        = ", casefile)
  println("opt_fd          = ", opt_fd, "   opt_sparse = ", opt_sparse, "   flatstart = ", opt_flatstart)
  println("autodamp        = ", autodamp, "   autodamp_min = ", autodamp_min)
  println("start_projection= ", start_projection, "   dc_angle_limit_deg = ", start_projection_dc_angle_limit_deg)
  println("matpower_ratio  = ", matpower_ratio)
  println("reference read  = ", _matpower_slack_reference(mpc).vm_pu, " pu / ", _matpower_slack_reference(mpc).va_deg, " deg")
  println("reference override = ", reference_override ? "enabled" : "disabled", reference_override ? " -> " * string(reference_vm_pu) * " pu / " * string(reference_va_deg) * " deg" : "")
  println("cooldown_iters  = ", cooldown_iters, "   q_hyst_pu = ", q_hyst_pu)
  println("lock PV->PQ     = ", collect(lock_pv_to_pq_buses))
  println("tol             = ", tol, "   max_ite = ", max_ite)
  println("seconds/method  = ", seconds, "   samples = ", samples)
  println("===============================================================\n")

  for m in methods
    method_enable_pq_gen_controllers = _enable_pq_gen_controllers_for_method(m, enable_pq_gen_controllers)
    benchable = @benchmarkable _with_log_table(logfile_, log_status_) do
      run_acpflow(
        casefile = casefile_,
      max_ite = max_ite_,
      tol = tol_,
      opt_fd = opt_fd_,
      opt_sparse = opt_sparse_,
      method = method_,
      autodamp = autodamp_,
      autodamp_min = autodamp_min_,
      start_projection = start_projection_,
      start_projection_try_dc_start = start_projection_try_dc_start_,
      start_projection_try_blend_scan = start_projection_try_blend_scan_,
      start_projection_blend_lambdas = start_projection_blend_lambdas_,
      start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg_,
      qlimit_start_iter = qlimit_start_iter_,
      qlimit_start_mode = qlimit_start_mode_,
      qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu_,
      qlimit_guard = qlimit_guard_,
      qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu_,
      qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode_,
      qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode_,
      qlimit_guard_max_switches = qlimit_guard_max_switches_,
      qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching_,
      qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations_,
      qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations_,
      qlimit_guard_violation_mode = qlimit_guard_violation_mode_,
      qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu_,
      qlimit_guard_log = qlimit_guard_log_,
      pv_table_rows = pv_table_rows_,
      opt_flatstart = opt_flatstart_,
      show_results = false,
      verbose = 0,
      cooldown_iters = cooldown_iters_,
      q_hyst_pu = q_hyst_pu_,
      qlimit_trace_buses = qlimit_trace_buses_,
      qlimit_lock_reason = qlimit_lock_reason_,
      lock_pv_to_pq_buses = lock_pv_to_pq_buses_,
      enable_pq_gen_controllers = enable_pq_gen_controllers_,
      bus_shunt_model = bus_shunt_model_,
      matpower_shift_sign = matpower_shift_sign_,
      matpower_shift_unit = matpower_shift_unit_,
      matpower_ratio = matpower_ratio_,
      reference_vm_pu = reference_vm_pu_,
      reference_va_deg = reference_va_deg_,
      matpower_pv_voltage_source = matpower_pv_voltage_source_,
      matpower_pv_voltage_mismatch_tol_pu = matpower_pv_voltage_mismatch_tol_pu_,
      flatstart_voltage_mode = flatstart_voltage_mode_,
        flatstart_angle_mode = flatstart_angle_mode_,
      )
    end setup = (casefile_ = $casefile;
    max_ite_ = $max_ite;
    tol_ = $tol;
    opt_fd_ = $opt_fd;
    opt_sparse_ = $opt_sparse;
    autodamp_ = $autodamp;
    autodamp_min_ = $autodamp_min;
    start_projection_ = $start_projection;
    start_projection_try_dc_start_ = $start_projection_try_dc_start;
    start_projection_try_blend_scan_ = $start_projection_try_blend_scan;
    start_projection_blend_lambdas_ = $start_projection_blend_lambdas;
    start_projection_dc_angle_limit_deg_ = $start_projection_dc_angle_limit_deg;
    qlimit_start_iter_ = $qlimit_start_iter;
    qlimit_start_mode_ = $qlimit_start_mode;
    qlimit_auto_q_delta_pu_ = $qlimit_auto_q_delta_pu;
    qlimit_guard_ = $qlimit_guard;
    qlimit_guard_min_q_range_pu_ = $qlimit_guard_min_q_range_pu;
    qlimit_guard_zero_range_mode_ = $qlimit_guard_zero_range_mode;
    qlimit_guard_narrow_range_mode_ = $qlimit_guard_narrow_range_mode;
    qlimit_guard_max_switches_ = $qlimit_guard_max_switches;
    qlimit_guard_freeze_after_repeated_switching_ = $qlimit_guard_freeze_after_repeated_switching;
    qlimit_guard_accept_bounded_violations_ = $qlimit_guard_accept_bounded_violations;
    qlimit_guard_max_remaining_violations_ = $qlimit_guard_max_remaining_violations;
    qlimit_guard_violation_mode_ = $qlimit_guard_violation_mode;
    qlimit_guard_violation_threshold_pu_ = $qlimit_guard_violation_threshold_pu;
    qlimit_guard_log_ = $qlimit_guard_log;
    pv_table_rows_ = $console_max_rows;
    opt_flatstart_ = $opt_flatstart;
    cooldown_iters_ = $cooldown_iters;
    q_hyst_pu_ = $q_hyst_pu;
    qlimit_trace_buses_ = $qlimit_trace_buses;
    qlimit_lock_reason_ = $qlimit_lock_reason;
    lock_pv_to_pq_buses_ = $lock_pv_to_pq_buses;
    enable_pq_gen_controllers_ = $method_enable_pq_gen_controllers;
    bus_shunt_model_ = $bus_shunt_model;
    matpower_shift_sign_ = $matpower_shift_sign;
    matpower_shift_unit_ = $matpower_shift_unit;
    matpower_ratio_ = $matpower_ratio;
    reference_vm_pu_ = $effective_reference_vm_pu;
    reference_va_deg_ = $effective_reference_va_deg;
    matpower_pv_voltage_source_ = $matpower_pv_voltage_source;
    matpower_pv_voltage_mismatch_tol_pu_ = $matpower_pv_voltage_mismatch_tol_pu;
    flatstart_voltage_mode_ = $flatstart_voltage_mode;
    flatstart_angle_mode_ = $flatstart_angle_mode;
    method_ = $m;
    logfile_ = $logfile;
    log_status_ = $log_status)

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
  total_s = time() - t0
  @printf("\ntotal runtime        : %.3f s\n", total_s)
  open(logfile, "a") do io
    @printf(io, "total runtime: %.3f s\n", total_s)
    _print_log_status(io, log_status)
  end
  _print_log_status(stdout, log_status)
  _emit_performance_summary(performance_profile; logfile = logfile, print_to_console = performance_print_to_console, write_to_logfile = performance_write_to_logfile, max_rows = performance_max_diagnostic_rows)
  return results
end

"""
    _matpower_auto_profile(mpc, cfg, yaml_cfg)

Build the optional MATPOWER example auto-profile for an imported case. The
profile inspects fixed-reference residuals, shunt impact, PV/REF voltage targets,
large-case flat-start risk, and narrow Q-limit patterns, then returns the
original configuration plus recommendations and any safely applied changes.
Explicit YAML keys are preserved in `apply` mode so user configuration remains
reproducible.
"""
function _matpower_auto_profile(mpc, cfg, yaml_cfg::Dict{String,Any})
  mode = _as_auto_profile_mode(cfg.matpower_auto_profile)
  explicit_keys = Set(Symbol(k) for k in keys(yaml_cfg))
  recommendations = Dict{Symbol,Any}()
  reasons = Dict{Symbol,String}()
  evidence = String[]
  applied = Dict{Symbol,Any}()
  preserved = Symbol[]

  if mode != :off
    _matpower_auto_profile_shift_scan!(recommendations, reasons, evidence, mpc, cfg)
    _matpower_auto_profile_shunt_scan!(recommendations, reasons, evidence, mpc, cfg)
    _matpower_auto_profile_voltage_source!(recommendations, reasons, evidence, mpc)
    _matpower_auto_profile_flatstart!(recommendations, reasons, evidence, mpc)
    _matpower_auto_profile_qlimits!(recommendations, reasons, evidence, mpc)

    if mode == :apply
      for option in sort!(collect(keys(recommendations)); by = String)
        if option in explicit_keys
          push!(preserved, option)
        else
          current = hasproperty(cfg, option) ? getproperty(cfg, option) : nothing
          recommended = recommendations[option]
          current == recommended || (applied[option] = recommended)
        end
      end
    end
  end

  cfg_out = isempty(applied) ? cfg : merge(cfg, _namedtuple_from_symbol_dict(applied))
  return (; mode, cfg = cfg_out, recommendations = _namedtuple_from_symbol_dict(recommendations), applied = _namedtuple_from_symbol_dict(applied), preserved = sort!(preserved; by = String), reasons, evidence)
end

function _print_matpower_auto_profile(io::IO, result)
  result.mode == :off && return nothing
  println(io, "==================== MATPOWER auto-profile ====================")
  println(io, "mode: ", result.mode)
  if isempty(result.evidence)
    println(io, "diagnostics: no evidence collected")
  else
    println(io, "diagnostics:")
    for item in result.evidence
      println(io, "  - ", item)
    end
  end
  rec_keys = sort!(collect(keys(pairs(result.recommendations))); by = String)
  if isempty(rec_keys)
    println(io, "recommendations: none; keeping user/default configuration")
  else
    println(io, "recommendations:")
    for option in rec_keys
      value = getproperty(result.recommendations, option)
      action = hasproperty(result.applied, option) ? "applied" : (option in result.preserved ? "preserved explicit value" : "recommended only")
      println(io, "  - ", option, " => ", value, " [", action, "]")
      println(io, "    reason: ", get(result.reasons, option, "diagnostic recommendation"))
    end
  end
  println(io, "================================================================\n")
  return nothing
end

function _auto_profile_recommendation_text(result, option::Symbol, yes::AbstractString, no::AbstractString = "OK")::String
  return hasproperty(result.recommendations, option) ? yes : no
end

function _auto_profile_evidence_with_prefix(result, prefix::AbstractString)::String
  for item in result.evidence
    startswith(item, prefix) && return item
  end
  return ""
end

function _print_matpower_auto_profile_compact(io::IO, result)
  result.mode == :off && return nothing
  println(io, "MATPOWER auto-profile: ", result.mode)
  branch_evidence = _auto_profile_evidence_with_prefix(result, "branch-shift scan")
  branch_match = match(r"branch-shift scan best='([^']+)'", branch_evidence)
  branch_text = isnothing(branch_match) ? "OK" : string("OK, ", branch_match.captures[1])
  shunt_text = _auto_profile_recommendation_text(result, :bus_shunt_model, "admittance recommended", "keep/admittance")
  pv_text = _auto_profile_recommendation_text(result, :matpower_pv_voltage_source, "GEN.VG recommended", "OK")
  qlimit_evidence = _auto_profile_evidence_with_prefix(result, "Q-limit scan")
  qlimit_text = isempty(qlimit_evidence) ? "OK" : replace(qlimit_evidence, "online generators have zero or narrow Q range" => "online generators have zero/narrow Q range")
  start_text = hasproperty(result.recommendations, :start_projection) || hasproperty(result.recommendations, :flatstart_angle_mode) ? "DC angle + blended voltage recommended" : "OK"
  @printf(io, "  branch shift      : %s\n", branch_text)
  @printf(io, "  bus shunts        : %s\n", shunt_text)
  @printf(io, "  PV voltage source : %s\n", pv_text)
  @printf(io, "  Q-limit scan      : %s\n", qlimit_text)
  @printf(io, "  start profile     : %s\n\n", start_text)
  return nothing
end

function _print_matpower_console_diagnostics_compact(io::IO, mpc; matpower_shift_sign::Real = 1.0, matpower_shift_unit = "deg", matpower_ratio = "normal")
  isnothing(mpc) && return nothing
  println(io, "Diagnostics:")
  vmva_chk = MatpowerIO.vmva_power_mismatch_stats(mpc; matpower_shift_sign = matpower_shift_sign, matpower_shift_unit = matpower_shift_unit, matpower_ratio = matpower_ratio)
  if get(vmva_chk, :ok, false)
    @printf(io, "  VM/VA self-check : max |dP|=%.4g MW, max |dQ|=%.4g MVAr\n", Float64(get(vmva_chk, :max_p_mis_MW, NaN)), Float64(get(vmva_chk, :max_q_mis_MVar, NaN)))
  else
    println(io, "  VM/VA self-check : skipped")
  end
  if hasproperty(mpc, :branch) && size(mpc.branch, 2) >= 4
    neg_r = count(<(0.0), mpc.branch[:, 3])
    neg_x = count(<(0.0), mpc.branch[:, 4])
    @printf(io, "  negative branches: BR_R<0=%d, BR_X<0=%d\n", neg_r, neg_x)
  end
  println(io)
  return nothing
end

function _qlimit_summary_counts(net)
  isnothing(net) && return (; pv2pq_events = 0, pv2pq_buses = 0, guarded_pv_buses = 0, oscillating_buses = 0)
  counts = Dict{Int,Int}()
  guarded = Set{Int}()
  for ev in net.qLimitLog
    counts[ev.bus] = get(counts, ev.bus, 0) + 1
    ev.iter == 0 && push!(guarded, ev.bus)
  end
  return (; pv2pq_events = length(net.qLimitLog), pv2pq_buses = length(net.qLimitEvents), guarded_pv_buses = length(guarded), oscillating_buses = count(>(1), values(counts)))
end

function _print_matpower_run_summary(io::IO, summaries; logfile::AbstractString = "", q_limit_mode::Symbol = :summary)
  println(io, "Run summary:")
  if isempty(summaries)
    println(io, "  status              : no methods executed")
  end
  for s in summaries
    compare_text = s.compare_status == :ok ? "OK" : s.compare_status == :warn ? "WARN" : s.compare_status == :skip ? "SKIP" : "FAIL"
    println(io, "  method              : ", s.method)
    println(io, "  numerical_solution  : ", s.numerical_solution == :ok ? "OK" : "FAIL")
    println(io, "  q_limit_active_set  : ", s.q_limit_active_set == :ok ? "OK" : "FAIL")
    println(io, "  final_converged     : ", s.final_converged)
    println(io, "  iterations          : ", s.iterations)
    @printf(io, "  time                : %.6f s\n", s.elapsed_s)
    println(io, "  compare             : ", compare_text)
    if !isnan(s.max_dvm)
      @printf(io, "  max |dVm|           : %.5g pu\n", s.max_dvm)
    end
    if !isnan(s.max_dva)
      @printf(io, "  max |dVa| aligned   : %.5g deg\n", s.max_dva)
    end
    if q_limit_mode != :off
      println(io, "  pv2pq_events        : ", s.pv2pq_events)
      println(io, "  pv2pq_buses         : ", s.pv2pq_buses)
      println(io, "  guarded_pv_buses    : ", s.guarded_pv_buses)
      println(io, "  oscillating_buses   : ", s.oscillating_buses)
    end
    if hasproperty(s, :status)
      println(io, "  status              : ", s.status)
    end
    println(io, "  reason              : ", s.reason_text)
  end
  !isempty(logfile) && println(io, "\nDetails written to logfile.")
  return nothing
end

function main()
  yaml_path = _yaml_path_from_inputs()
  yaml_cfg = load_yaml_config(yaml_path)
  case = String(get(yaml_cfg, "case", DEFAULT_CASE))
  methods_cfg = haskey(yaml_cfg, "methods") ? _as_symbol_vec(yaml_cfg["methods"]) : Symbol[]
  methods = isempty(methods_cfg) ? DEFAULT_METHODS : methods_cfg

  print("\e[2J\e[H") # clear screen and move cursor to home position
  println("----------------------------------------------------------------------------------------------")

  timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
  log_case = replace(basename(case), r"[^A-Za-z0-9_.-]" => "_")
  logfile = joinpath(OUTDIR, "run_$(log_case)_$(timestamp).log")

  # This will:
  # - download case14.m into data/mpower/ if missing
  # - generate case14.jl if requested and missing
  local_case = Sparlectra.FetchMatpowerCase.ensure_casefile(case; outdir = Sparlectra.MPOWER_DIR, to_jl = true, overwrite = false)
  println("Sparlectra version: ", Sparlectra.version())
  println("casefile: ", local_case)
  !isempty(yaml_path) && println("yaml: ", yaml_path)
  println("logfile: ", logfile, "\n")
  log_status = _MatpowerImportLogStatus()
  cfg = bench_config_for_case(case, yaml_cfg)
  performance_profile = _new_performance_profile(cfg)
  open(logfile, "w") do io
    println(io, "Sparlectra version: ", Sparlectra.version())
    println(io, "casefile: ", local_case)
    println(io, "timestamp: ", Dates.now())
  end
  mpc = _perf_profile_time!(performance_profile, :matpower_parse) do
    _with_log_table(logfile, log_status) do
      MatpowerIO.read_case(local_case)
    end
  end

  auto_profile = _matpower_auto_profile(mpc, cfg, yaml_cfg)
  cfg = auto_profile.cfg
  performance_profile[:enabled] = cfg.performance_enabled && cfg.performance_level != :off
  performance_profile[:level] = cfg.performance_level
  performance_profile[:show_allocations] = cfg.performance_show_allocations
  performance_profile[:show_iteration_table] = cfg.performance_show_iteration_table
  if cfg.performance_enabled && cfg.performance_compact_logging
    cfg = merge(cfg, (; console_diagnostics = :compact, console_q_limit_events = :summary, console_max_rows = min(cfg.console_max_rows, cfg.performance_max_diagnostic_rows)))
  end
  if cfg.matpower_auto_profile_log && auto_profile.mode != :off
    if cfg.console_auto_profile == :full
      _print_matpower_auto_profile(stdout, auto_profile)
    elseif cfg.console_auto_profile != :off
      _print_matpower_auto_profile_compact(stdout, auto_profile)
    end
    open(logfile, "a") do io
      _print_matpower_auto_profile(io, auto_profile)
    end
  end
  if cfg.console_diagnostics == :compact || cfg.console_diagnostics == :summary
    _print_matpower_console_diagnostics_compact(stdout, mpc; matpower_shift_sign = cfg.matpower_shift_sign, matpower_shift_unit = cfg.matpower_shift_unit, matpower_ratio = cfg.matpower_ratio)
  end
  _warn_if_flatstart_uses_only_voltage_setpoints(case, cfg, mpc)
  if cfg.trace_legacy_bus_type_warnings
    ENV["SPARLECTRA_TRACE_LEGACY_BUSTYPE"] = "1"
    println("legacy trace     = enabled (SPARLECTRA_TRACE_LEGACY_BUSTYPE=1)")
  else
    delete!(ENV, "SPARLECTRA_TRACE_LEGACY_BUSTYPE")
  end
  lock_pv_to_pq_buses = cfg.lock_pv_to_pq_buses
  qlimit_lock_reason = :manual
  if cfg.ignore_q_limits
    lock_pv_to_pq_buses = _mpc_pv_bus_ids(mpc)
    qlimit_lock_reason = :ignore_q_limits
  end

  bench = Base.invokelatest(
    getfield(@__MODULE__, :bench_run_acpflow);
    casefile = local_case,
    methods = methods,
    mpc = mpc,
    logfile = logfile,
    show_diff = cfg.show_diff,
    tol_vm = cfg.tol_vm,
    tol_va = cfg.tol_va,
    max_ite = cfg.max_ite,
    tol = cfg.tol,
    opt_fd = cfg.opt_fd,
    opt_sparse = cfg.opt_sparse,
    opt_flatstart = cfg.opt_flatstart,
    autodamp = cfg.autodamp,
    autodamp_min = cfg.autodamp_min,
    start_projection = cfg.start_projection,
    start_projection_try_dc_start = cfg.start_projection_try_dc_start,
    start_projection_try_blend_scan = cfg.start_projection_try_blend_scan,
    start_projection_branch_guard = cfg.start_projection_branch_guard,
    start_projection_measure_candidates = cfg.start_projection_measure_candidates,
    start_projection_reuse_import_data = cfg.start_projection_reuse_import_data,
    start_projection_blend_lambdas = cfg.start_projection_blend_lambdas,
    start_projection_dc_angle_limit_deg = cfg.start_projection_dc_angle_limit_deg,
    qlimit_start_iter = cfg.qlimit_start_iter,
    qlimit_start_mode = cfg.qlimit_start_mode,
    qlimit_auto_q_delta_pu = cfg.qlimit_auto_q_delta_pu,
    qlimit_guard = cfg.qlimit_guard,
    qlimit_guard_min_q_range_pu = cfg.qlimit_guard_min_q_range_pu,
    qlimit_guard_zero_range_mode = cfg.qlimit_guard_zero_range_mode,
    qlimit_guard_narrow_range_mode = cfg.qlimit_guard_narrow_range_mode,
    qlimit_guard_max_switches = cfg.qlimit_guard_max_switches,
    qlimit_guard_freeze_after_repeated_switching = cfg.qlimit_guard_freeze_after_repeated_switching,
    qlimit_guard_accept_bounded_violations = cfg.qlimit_guard_accept_bounded_violations,
    qlimit_guard_max_remaining_violations = cfg.qlimit_guard_max_remaining_violations,
    qlimit_guard_violation_mode = cfg.qlimit_guard_violation_mode,
    qlimit_guard_violation_threshold_pu = cfg.qlimit_guard_violation_threshold_pu,
    qlimit_guard_log = cfg.qlimit_guard_log,
    verbose = cfg.verbose,
    cooldown_iters = cfg.cooldown_iters,
    q_hyst_pu = cfg.q_hyst_pu,
    qlimit_trace_buses = cfg.qlimit_trace_buses,
    qlimit_lock_reason = qlimit_lock_reason,
    lock_pv_to_pq_buses = lock_pv_to_pq_buses,
    seconds = cfg.seconds,
    samples = cfg.samples,
    show_once = cfg.show_once,
    show_once_output = cfg.show_once_output,
    show_once_max_nodes = cfg.show_once_max_nodes,
    matpower_shift_sign = cfg.matpower_shift_sign,
    matpower_shift_unit = cfg.matpower_shift_unit,
    matpower_ratio = cfg.matpower_ratio,
    reference_override = cfg.reference_override,
    reference_vm_pu = cfg.reference_vm_pu,
    reference_va_deg = cfg.reference_va_deg,
    diagnose_matpower_reference = cfg.diagnose_matpower_reference,
    diagnose_branch_shift_conventions = cfg.diagnose_branch_shift_conventions,
    diagnose_branch_neighborhood = cfg.diagnose_branch_neighborhood,
    diagnose_residual_clusters = cfg.diagnose_residual_clusters,
    diagnose_residual_cluster_threshold_mw = cfg.diagnose_residual_cluster_threshold_mw,
    diagnose_residual_cluster_maxlines = cfg.diagnose_residual_cluster_maxlines,
    diagnose_branch_neighborhood_buses = cfg.diagnose_branch_neighborhood_buses,
    diagnose_branch_neighborhood_depth = cfg.diagnose_branch_neighborhood_depth,
    diagnose_branch_neighborhood_maxlines = cfg.diagnose_branch_neighborhood_maxlines,
    diagnose_nodal_balance_breakdown = cfg.diagnose_nodal_balance_breakdown,
    diagnose_nodal_balance_buses = cfg.diagnose_nodal_balance_buses,
    diagnose_nodal_balance_maxlines = cfg.diagnose_nodal_balance_maxlines,
    diagnose_nodal_balance_include_branches = cfg.diagnose_nodal_balance_include_branches,
    diagnose_nodal_balance_include_generators = cfg.diagnose_nodal_balance_include_generators,
    diagnose_nodal_balance_include_shunts = cfg.diagnose_nodal_balance_include_shunts,
    diagnose_negative_branch_impedance = cfg.diagnose_negative_branch_impedance,
    diagnose_negative_branch_impedance_maxlines = cfg.diagnose_negative_branch_impedance_maxlines,
    diagnose_negative_branch_impedance_fail_on_negative_r = cfg.diagnose_negative_branch_impedance_fail_on_negative_r,
    diagnose_negative_branch_impedance_fail_on_negative_x = cfg.diagnose_negative_branch_impedance_fail_on_negative_x,
    diagnose_negative_branch_impedance_warn_threshold_abs_r = cfg.diagnose_negative_branch_impedance_warn_threshold_abs_r,
    diagnose_negative_branch_impedance_warn_threshold_abs_x = cfg.diagnose_negative_branch_impedance_warn_threshold_abs_x,
    diagnose_maxlines = cfg.diagnose_maxlines,
    log_effective_config = cfg.log_effective_config,
    yaml_path = yaml_path,
    effective_config = cfg,
    benchmark = cfg.benchmark,
    enable_pq_gen_controllers = cfg.enable_pq_gen_controllers,
    bus_shunt_model = cfg.bus_shunt_model,
    matpower_pv_voltage_source = cfg.matpower_pv_voltage_source,
    matpower_pv_voltage_mismatch_tol_pu = cfg.matpower_pv_voltage_mismatch_tol_pu,
    compare_voltage_reference = cfg.compare_voltage_reference,
    diagnose_pv_voltage_references = cfg.diagnose_pv_voltage_references,
    diagnose_pv_voltage_maxlines = cfg.diagnose_pv_voltage_maxlines,
    flatstart_voltage_mode = cfg.flatstart_voltage_mode,
    flatstart_angle_mode = cfg.flatstart_angle_mode,
    wrong_branch_detection = cfg.wrong_branch_detection,
    wrong_branch_min_vm_pu = cfg.wrong_branch_min_vm_pu,
    wrong_branch_max_angle_spread_deg = cfg.wrong_branch_max_angle_spread_deg,
    wrong_branch_rescue = cfg.wrong_branch_rescue,
    wrong_branch_rescue_modes = cfg.wrong_branch_rescue_modes,
    console_summary = cfg.console_summary,
    console_auto_profile = cfg.console_auto_profile,
    console_diagnostics = cfg.console_diagnostics,
    console_q_limit_events = cfg.console_q_limit_events,
    console_max_rows = cfg.console_max_rows,
    logfile_diagnostics = cfg.logfile_diagnostics,
    performance_profile = performance_profile,
    performance_level = cfg.performance_level,
    performance_print_to_console = cfg.performance_print_to_console,
    performance_write_to_logfile = cfg.performance_write_to_logfile,
    performance_show_iteration_table = cfg.performance_show_iteration_table,
    performance_skip_reference_comparison = cfg.performance_skip_reference_comparison,
    performance_skip_expensive_diagnostics = cfg.performance_skip_expensive_diagnostics,
    performance_skip_branch_neighborhood_report = cfg.performance_skip_branch_neighborhood_report,
    performance_max_diagnostic_rows = cfg.performance_max_diagnostic_rows,
    log_status = log_status,
  )
  return bench
end

if get(ENV, "SPARLECTRA_MATPOWER_IMPORT_NO_MAIN", "") != "1"
  Base.invokelatest(getfield(@__MODULE__, :main))
end

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
    start_projection_blend_lambdas = [0.25, 0.5, 0.75],
    start_projection_dc_angle_limit_deg = 60.0,
    qlimit_start_iter = 2,
    qlimit_start_mode = :iteration,
    qlimit_auto_q_delta_pu = 1e-4,
    verbose = 1,
    cooldown_iters = 0,
    q_hyst_pu = 0.0,
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
    reference_vm_pu = 1.0,
    reference_va_deg = 0.0,
    diagnose_matpower_reference = false,
    diagnose_branch_shift_conventions = false,
    diagnose_maxlines = 12,
    log_effective_config = false,
    enable_pq_gen_controllers = true,
    bus_shunt_model = "admittance",
    trace_legacy_bus_type_warnings = false,
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
      start_projection_blend_lambdas = Float64[x for x in get(yaml_cfg, "start_projection_blend_lambdas", base.start_projection_blend_lambdas)],
      start_projection_dc_angle_limit_deg = Float64(get(yaml_cfg, "start_projection_dc_angle_limit_deg", base.start_projection_dc_angle_limit_deg)),
      qlimit_start_iter = Int(get(yaml_cfg, "qlimit_start_iter", base.qlimit_start_iter)),
      qlimit_start_mode = Symbol(get(yaml_cfg, "qlimit_start_mode", base.qlimit_start_mode)),
      qlimit_auto_q_delta_pu = Float64(get(yaml_cfg, "qlimit_auto_q_delta_pu", base.qlimit_auto_q_delta_pu)),
      verbose = Int(get(yaml_cfg, "verbose", base.verbose)),
      cooldown_iters = Int(get(yaml_cfg, "cooldown_iters", base.cooldown_iters)),
      q_hyst_pu = Float64(get(yaml_cfg, "q_hyst_pu", base.q_hyst_pu)),
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
      reference_vm_pu = Float64(get(yaml_cfg, "reference_vm_pu", base.reference_vm_pu)),
      reference_va_deg = Float64(get(yaml_cfg, "reference_va_deg", base.reference_va_deg)),
      diagnose_matpower_reference = Bool(get(yaml_cfg, "diagnose_matpower_reference", base.diagnose_matpower_reference)),
      diagnose_branch_shift_conventions = Bool(get(yaml_cfg, "diagnose_branch_shift_conventions", base.diagnose_branch_shift_conventions)),
      diagnose_maxlines = Int(get(yaml_cfg, "diagnose_maxlines", base.diagnose_maxlines)),
      log_effective_config = Bool(get(yaml_cfg, "log_effective_config", base.log_effective_config)),
      enable_pq_gen_controllers = Bool(get(yaml_cfg, "enable_pq_gen_controllers", base.enable_pq_gen_controllers)),
      bus_shunt_model = String(get(yaml_cfg, "bus_shunt_model", base.bus_shunt_model)),
      trace_legacy_bus_type_warnings = Bool(get(yaml_cfg, "trace_legacy_bus_type_warnings", base.trace_legacy_bus_type_warnings)),
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
  for rank in 1:n_p
    r = p_order[rank]
    @printf(io, "%5d %6d %5d %13.6f %13.3f %10.6f %12.6f\n", rank, Int(bus[r, 1]), Int(bus[r, 2]), real(mis[r]), real(mis[r]) * baseMVA, bus[r, 8], bus[r, 9])
  end
  println(io, "Top reactive-power residuals on enforced PQ buses:")
  println(io, " rank  bus      type        dQ_pu       dQ_MVAr      Vm_ref     Va_ref_deg")
  for rank in 1:n_q
    r = q_order[rank]
    @printf(io, "%5d %6d %5d %13.6f %13.3f %10.6f %12.6f\n", rank, Int(bus[r, 1]), Int(bus[r, 2]), imag(mis[r]), imag(mis[r]) * baseMVA, bus[r, 8], bus[r, 9])
  end
  println(io, "================================================================================\n")
end

function _print_matpower_reference_diagnostics(io::IO, mpc; matpower_shift_sign::Real = 1.0, matpower_shift_unit = "deg", matpower_ratio = "normal", diagnose_branch_shift_conventions::Bool = false, maxlines::Int = 12)
  mp_has_vm_va(mpc) || return nothing
  base_diag = Base.invokelatest(getfield(@__MODULE__, :_matpower_reference_residuals), mpc; matpower_shift_sign = matpower_shift_sign, matpower_shift_unit = matpower_shift_unit, matpower_ratio = matpower_ratio)
  label = "configured SHIFT sign=$(matpower_shift_sign), unit=$(matpower_shift_unit), ratio=$(matpower_ratio)"
  Base.invokelatest(getfield(@__MODULE__, :_print_top_residuals), io, label, base_diag, mpc.baseMVA; maxlines = maxlines)

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

function _show_once_summary_row(method::Symbol, status, stats, cmp_ok::Bool; compare_available::Bool)
  converged = Bool(get(status, :converged, false))
  iterations = _run_iterations(status)
  elapsed_s = _run_elapsed_s(status)
  if compare_available
    return (method = method, converged = converged, iterations = iterations, elapsed_s = elapsed_s, max_dvm = Float64(get(stats, :max_dvm, NaN)), max_dva = Float64(get(stats, :max_dva, NaN)), cmp_ok = cmp_ok)
  end
  return (method = method, converged = converged, iterations = iterations, elapsed_s = elapsed_s, max_dvm = NaN, max_dva = NaN, cmp_ok = false)
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
  start_projection_blend_lambdas::AbstractVector{<:Real} = [0.25, 0.5, 0.75],
  start_projection_dc_angle_limit_deg::Float64 = 60.0,
  qlimit_start_iter::Int = 2,
  qlimit_start_mode::Symbol = :iteration,
  qlimit_auto_q_delta_pu::Float64 = 1e-4,
  verbose::Int = 0,
  cooldown_iters::Int = 0,
  q_hyst_pu::Float64 = 0.0,
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
  reference_vm_pu::Float64 = 1.0,
  reference_va_deg::Float64 = 0.0,
  diagnose_matpower_reference::Bool = false,
  diagnose_branch_shift_conventions::Bool = false,
  diagnose_maxlines::Int = 12,
  log_effective_config::Bool = false,
  yaml_path::String = "",
  effective_config = nothing,
  enable_pq_gen_controllers::Bool = true,
  bus_shunt_model::String = "admittance",
)
  t0 = time()
  results = Dict{Symbol,Any}()

  # Create/seed logfile up-front so users can see progress even if benchmarks are long.
  open(logfile, "w") do io
    println(io, "Sparlectra version: ", Sparlectra.version())
    println(io, "casefile: ", casefile)
    println(io, "timestamp: ", Dates.now())
    println(io)
    if log_effective_config && !isnothing(effective_config)
      Base.invokelatest(getfield(@__MODULE__, :_print_effective_config), io, effective_config; yaml_path = yaml_path, case_name = casefile, methods = methods)
    end
    if !isnothing(mpc)
      Base.invokelatest(getfield(@__MODULE__, :_print_vmva_self_check), io, mpc; matpower_shift_sign = matpower_shift_sign, matpower_shift_unit = matpower_shift_unit, matpower_ratio = matpower_ratio)
      if diagnose_matpower_reference
        Base.invokelatest(
          getfield(@__MODULE__, :_print_matpower_reference_diagnostics),
          io,
          mpc;
          matpower_shift_sign = matpower_shift_sign,
          matpower_shift_unit = matpower_shift_unit,
          matpower_ratio = matpower_ratio,
          diagnose_branch_shift_conventions = diagnose_branch_shift_conventions,
          maxlines = diagnose_maxlines,
        )
      end
    end
  end

  # Optional: show results once (not benchmarked)
  if show_once
    local summaries = Vector{NamedTuple{(:method, :converged, :iterations, :elapsed_s, :max_dvm, :max_dva, :cmp_ok),Tuple{Symbol,Bool,Int,Float64,Float64,Float64,Bool}}}()
    show_classic = (show_once_output == :classic)
    println("show_once output is written to logfile: ", logfile)
    for m in methods
      open(logfile, "a") do io
        redirect_stdout(io) do
          redirect_stderr(io) do
            println("=================================================================")
            println("RUN method = ", m)
            println("=================================================================\n")
            status_ref = Ref{Any}(nothing)
            net_res = run_acpflow(
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
              start_projection_blend_lambdas = start_projection_blend_lambdas,
              start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
              qlimit_start_iter = qlimit_start_iter,
              qlimit_start_mode = qlimit_start_mode,
              qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu,
              opt_flatstart = opt_flatstart,
              show_results = show_classic,
              show_compact_result = true,
              status_ref = status_ref,
              verbose = verbose,
              cooldown_iters = cooldown_iters,
              q_hyst_pu = q_hyst_pu,
              lock_pv_to_pq_buses = lock_pv_to_pq_buses,
              enable_pq_gen_controllers = _enable_pq_gen_controllers_for_method(m, enable_pq_gen_controllers),
              bus_shunt_model = bus_shunt_model,
              matpower_shift_sign = matpower_shift_sign,
              matpower_shift_unit = matpower_shift_unit,
              matpower_ratio = matpower_ratio,
              reference_vm_pu = reference_vm_pu,
              reference_va_deg = reference_va_deg,
            )
            status = status_ref[]
            _print_converged_loss_summary(io, m, status, net_res)
            if !show_classic
              _print_dataframe_nodes(io, net_res; max_nodes = show_once_max_nodes)
            end
            if mp_has_vm_va(mpc)
              ok, stats = MatpowerIO.compare_vm_va(net_res, mpc; show_diff = show_diff, tol_vm = tol_vm, tol_va = tol_va, maxlines = 20)
              push!(summaries, _show_once_summary_row(m, status, stats, ok; compare_available = true))
            else
              push!(summaries, _show_once_summary_row(m, status, nothing, false; compare_available = false))
              println("Compare skipped: no solution-like VM/VA in mpc.bus(:,8:9)")
            end
            if !Bool(get(status, :converged, false)) && !enable_pq_gen_controllers && m === :rectangular
              println("Fallback diagnostic: rerun with enable_pq_gen_controllers=true")
              fb_ref = Ref{Any}(nothing)
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
                start_projection_blend_lambdas = start_projection_blend_lambdas,
                start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
                qlimit_start_iter = qlimit_start_iter,
                qlimit_start_mode = qlimit_start_mode,
                qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu,
                opt_flatstart = opt_flatstart,
                show_results = false,
                show_compact_result = true,
                status_ref = fb_ref,
                verbose = 0,
                cooldown_iters = cooldown_iters,
                q_hyst_pu = q_hyst_pu,
                lock_pv_to_pq_buses = lock_pv_to_pq_buses,
                enable_pq_gen_controllers = true,
                bus_shunt_model = bus_shunt_model,
                matpower_shift_sign = matpower_shift_sign,
                matpower_shift_unit = matpower_shift_unit,
                matpower_ratio = matpower_ratio,
                reference_vm_pu = reference_vm_pu,
                reference_va_deg = reference_va_deg,
              )
            end
            @printf(io, "show_once method=%s output=%s\n\n", String(m), String(show_once_output))
          end
        end
      end
    end

    # print short summary to terminal
    println("\n==================== Summary ====================")
    println("logfile: ", logfile)
    for s in summaries
      if isnan(s.max_dvm) || isnan(s.max_dva)
        @printf("method=%-12s  converged=%s  iterations=%d  time=%8.6f s  compare=SKIP\n", String(s.method), s.converged ? "yes" : "no", s.iterations, s.elapsed_s)
      else
        @printf("method=%-12s  converged=%s  iterations=%d  time=%8.6f s  compare=%s  max|dVm|=%8.5f pu  max|dVa|=%7.4f deg\n", String(s.method), s.converged ? "yes" : "no", s.iterations, s.elapsed_s, s.cmp_ok ? "OK " : "FAIL", s.max_dvm, s.max_dva)
      end
    end
    println("=================================================\n")
  end

  if !show_once
    println("\nInitial run (convergence + iterations):")
    for m in methods
      status_ref = Ref{Any}(nothing)
      net_res = run_acpflow(
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
        start_projection_blend_lambdas = start_projection_blend_lambdas,
        start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
        qlimit_start_iter = qlimit_start_iter,
        qlimit_start_mode = qlimit_start_mode,
        qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu,
        opt_flatstart = opt_flatstart,
        show_results = false,
        show_compact_result = true,
        status_ref = status_ref,
        verbose = 0,
        cooldown_iters = cooldown_iters,
        q_hyst_pu = q_hyst_pu,
        lock_pv_to_pq_buses = lock_pv_to_pq_buses,
        enable_pq_gen_controllers = _enable_pq_gen_controllers_for_method(m, enable_pq_gen_controllers),
        bus_shunt_model = bus_shunt_model,
        matpower_shift_sign = matpower_shift_sign,
        matpower_shift_unit = matpower_shift_unit,
        matpower_ratio = matpower_ratio,
        reference_vm_pu = reference_vm_pu,
        reference_va_deg = reference_va_deg,
      )
      _print_converged_loss_summary(stdout, m, status_ref[], net_res)
      open(logfile, "a") do io
        _print_converged_loss_summary(io, m, status_ref[], net_res)
      end
    end
  end

  if !benchmark
    total_s = time() - t0
    println("\nBenchmarkTools loop skipped (`benchmark: false`).")
    @printf("total runtime    = %.3f s\n", total_s)
    println("logfile: ", logfile)
    open(logfile, "a") do io
      @printf(io, "total runtime: %.3f s\n", total_s)
    end
    return results
  end

  # Warmup (compile) once per method with minimal output
  println("Warmup run:")
  for m in methods
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
      start_projection_blend_lambdas = start_projection_blend_lambdas,
      start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
      qlimit_start_iter = qlimit_start_iter,
      qlimit_start_mode = qlimit_start_mode,
      qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu,
      opt_flatstart = opt_flatstart,
      show_results = false,
      verbose = 0,
      cooldown_iters = cooldown_iters,
      q_hyst_pu = q_hyst_pu,
      lock_pv_to_pq_buses = lock_pv_to_pq_buses,
      enable_pq_gen_controllers = _enable_pq_gen_controllers_for_method(m, enable_pq_gen_controllers),
      bus_shunt_model = bus_shunt_model,
      matpower_shift_sign = matpower_shift_sign,
      matpower_shift_unit = matpower_shift_unit,
      matpower_ratio = matpower_ratio,
      reference_vm_pu = reference_vm_pu,
      reference_va_deg = reference_va_deg,
    )
  end

  println("\n==================== Benchmark run_acpflow ====================")
  println("casefile        = ", casefile)
  println("opt_fd          = ", opt_fd, "   opt_sparse = ", opt_sparse, "   flatstart = ", opt_flatstart)
  println("autodamp        = ", autodamp, "   autodamp_min = ", autodamp_min)
  println("start_projection= ", start_projection, "   dc_angle_limit_deg = ", start_projection_dc_angle_limit_deg)
  println("matpower_ratio  = ", matpower_ratio)
  println("reference       = ", reference_vm_pu, " pu / ", reference_va_deg, " deg")
  println("cooldown_iters  = ", cooldown_iters, "   q_hyst_pu = ", q_hyst_pu)
  println("lock PV->PQ     = ", collect(lock_pv_to_pq_buses))
  println("tol             = ", tol, "   max_ite = ", max_ite)
  println("seconds/method  = ", seconds, "   samples = ", samples)
  println("===============================================================\n")

  for m in methods
    method_enable_pq_gen_controllers = _enable_pq_gen_controllers_for_method(m, enable_pq_gen_controllers)
    benchable = @benchmarkable run_acpflow(
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
      opt_flatstart = opt_flatstart_,
      show_results = false,
      verbose = 0,
      cooldown_iters = cooldown_iters_,
      q_hyst_pu = q_hyst_pu_,
      lock_pv_to_pq_buses = lock_pv_to_pq_buses_,
      enable_pq_gen_controllers = enable_pq_gen_controllers_,
      bus_shunt_model = bus_shunt_model_,
      matpower_shift_sign = matpower_shift_sign_,
      matpower_shift_unit = matpower_shift_unit_,
      matpower_ratio = matpower_ratio_,
      reference_vm_pu = reference_vm_pu_,
      reference_va_deg = reference_va_deg_,
    ) setup = (casefile_ = $casefile;
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
    opt_flatstart_ = $opt_flatstart;
    cooldown_iters_ = $cooldown_iters;
    q_hyst_pu_ = $q_hyst_pu;
    lock_pv_to_pq_buses_ = $lock_pv_to_pq_buses;
    enable_pq_gen_controllers_ = $method_enable_pq_gen_controllers;
    bus_shunt_model_ = $bus_shunt_model;
    matpower_shift_sign_ = $matpower_shift_sign;
    matpower_shift_unit_ = $matpower_shift_unit;
    matpower_ratio_ = $matpower_ratio;
    reference_vm_pu_ = $reference_vm_pu;
    reference_va_deg_ = $reference_va_deg;
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
  total_s = time() - t0
  @printf("\ntotal runtime    = %.3f s\n", total_s)
  println("logfile: ", logfile)
  open(logfile, "a") do io
    @printf(io, "total runtime: %.3f s\n", total_s)
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
  log_case = replace(basename(case), r"[^A-Za-z0-9_.-]" => "_")
  logfile = joinpath(OUTDIR, "run_$(log_case)_$(timestamp).log")

  # This will:
  # - download case14.m into data/mpower/ if missing
  # - generate case14.jl if requested and missing
  local_case = Sparlectra.FetchMatpowerCase.ensure_casefile(case; outdir = Sparlectra.MPOWER_DIR, to_jl = true, overwrite = false)
  mpc = MatpowerIO.read_case(local_case)

  cfg = bench_config_for_case(case, yaml_cfg)
  _warn_if_flatstart_uses_only_voltage_setpoints(case, cfg, mpc)
  if cfg.trace_legacy_bus_type_warnings
    ENV["SPARLECTRA_TRACE_LEGACY_BUSTYPE"] = "1"
    println("legacy trace     = enabled (SPARLECTRA_TRACE_LEGACY_BUSTYPE=1)")
  else
    delete!(ENV, "SPARLECTRA_TRACE_LEGACY_BUSTYPE")
  end
  lock_pv_to_pq_buses = cfg.lock_pv_to_pq_buses
  if cfg.ignore_q_limits
    lock_pv_to_pq_buses = _mpc_pv_bus_ids(mpc)
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
      start_projection_blend_lambdas = cfg.start_projection_blend_lambdas,
      start_projection_dc_angle_limit_deg = cfg.start_projection_dc_angle_limit_deg,
      qlimit_start_iter = cfg.qlimit_start_iter,
      qlimit_start_mode = cfg.qlimit_start_mode,
      qlimit_auto_q_delta_pu = cfg.qlimit_auto_q_delta_pu,
      verbose = cfg.verbose,
      cooldown_iters = cfg.cooldown_iters,
      q_hyst_pu = cfg.q_hyst_pu,
      lock_pv_to_pq_buses = lock_pv_to_pq_buses,
      seconds = cfg.seconds,
      samples = cfg.samples,
      show_once = cfg.show_once,
      show_once_output = cfg.show_once_output,
      show_once_max_nodes = cfg.show_once_max_nodes,
      matpower_shift_sign = cfg.matpower_shift_sign,
      matpower_shift_unit = cfg.matpower_shift_unit,
      matpower_ratio = cfg.matpower_ratio,
      reference_vm_pu = cfg.reference_vm_pu,
      reference_va_deg = cfg.reference_va_deg,
      diagnose_matpower_reference = cfg.diagnose_matpower_reference,
      diagnose_branch_shift_conventions = cfg.diagnose_branch_shift_conventions,
      diagnose_maxlines = cfg.diagnose_maxlines,
      log_effective_config = cfg.log_effective_config,
      yaml_path = yaml_path,
      effective_config = cfg,
      benchmark = cfg.benchmark,
      enable_pq_gen_controllers = cfg.enable_pq_gen_controllers,
      bus_shunt_model = cfg.bus_shunt_model,
  )
  return bench
end

if get(ENV, "SPARLECTRA_MATPOWER_IMPORT_NO_MAIN", "") != "1"
  Base.invokelatest(getfield(@__MODULE__, :main))
end

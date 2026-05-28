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

#
# file: src/jacobian_complex.jl
#

# jacobian_complex.jl — Complex-State Newton-Raphson Power Flow Formulation
#
# This module implements a Newton-Raphson power flow solver using complex voltages
# in rectangular coordinates (Vr + jVi) as state variables, as an alternative to
# the conventional polar formulation (Vm, θ).
#
# Features:
# - Rectangular complex-state Newton-Raphson with PQ and PV bus handling
# - Wirtinger calculus-based Jacobian construction for complex power equations
# - Active-set Q-limit management for PV→PQ switching with optional re-enable
# - Sparse analytic Jacobian construction
# - Hysteresis and cooldown mechanisms for robust PV bus management
# - Direct integration with Sparlectra.jl network data structures
#
# Mathematical Foundation:
# - State vector: x = [Vr(non-slack); Vi(non-slack)] ∈ ℝ^(2(n-1))
# - Complex power: S = V .* conj(Y * V) where Y is the bus admittance matrix  
# - PQ buses: ΔP = Re(S_calc - S_spec), ΔQ = Im(S_calc - S_spec)
# - PV buses: ΔP = Re(S_calc - S_spec), ΔV = |V| - V_set
# - Slack bus voltage is held constant throughout the iteration
#
# Key Functions:
# - runpf_rectangular!(): Main solver interface (network-integrated rectangular NR)
# - build_complex_jacobian(): Wirtinger-based Jacobian block construction
# - mismatch_rectangular(): Residual function for PQ/PV bus constraints
#
# Note:
# - The power-flow core uses sparse Y-bus and Jacobian matrices by default.
# - Dense builders remain for small helper-level diagnostics, but unsupported PF
#   entry options are rejected before the solver core is entered.
#
# References:
# - Wirtinger calculus for complex derivatives

using LinearAlgebra
using SparseArrays
using Printf

_is_rectangular_linear_step_failure(e) = e isa LinearAlgebra.SingularException || e isa LinearAlgebra.LAPACKException

function _print_rectangular_qlimit_summary(
  io::IO,
  net::Net,
  V::Vector{ComplexF64},
  Sbus_pu::Vector{ComplexF64},
  bus_types::Vector{Symbol},
  qmin_pu::AbstractVector,
  qmax_pu::AbstractVector,
  Qload_pu::Vector{Float64};
  q_hyst_pu::Float64,
  tolerance_pu::Float64 = 0.0,
  max_rows::Int = 30,
  max_console_rows::Union{Nothing,Int} = nothing,
)
  checked_pv = 0
  checked_ref = 0
  rows = NamedTuple[]
  for bus in eachindex(bus_types)
    bus_type = bus_types[bus]
    bus_type in (:PV, :Slack) || continue
    (bus <= length(qmin_pu)) && (bus <= length(qmax_pu)) || continue
    isfinite(qmin_pu[bus]) && isfinite(qmax_pu[bus]) || continue

    if bus_type == :PV
      checked_pv += 1
    else
      checked_ref += 1
    end

    qcalc = imag(Sbus_pu[bus]) + Qload_pu[bus]
    qmin = qmin_pu[bus]
    qmax = qmax_pu[bus]

    side = :none
    amount = 0.0
    if qcalc < qmin - tolerance_pu
      side = :low
      amount = qmin - qcalc
    elseif qcalc > qmax + tolerance_pu
      side = :high
      amount = qcalc - qmax
    end
    side == :none && continue

    last_it = lastQLimitIter(net, bus)
    lo_return, hi_return = q_limit_band(qmin_pu, qmax_pu, bus, q_hyst_pu)
    in_return_band = (qcalc > lo_return) && (qcalc < hi_return)
    push!(rows, (bus = bus, busI = _qlimit_original_bus_id(net, bus), type = bus_type, qcalc = qcalc, qmin = qmin, qmax = qmax, side = side, amount = amount, switched = !isnothing(last_it), last_it = last_it, in_return_band = in_return_band))
  end

  sort!(rows; by = r -> -abs(r.amount * net.baseMVA))
  pv_violations = count(r -> r.type == :PV, rows)
  ref_violations = count(r -> r.type == :Slack, rows)
  # REF/slack buses anchor the angle reference and are not moved through the
  # PV→PQ active-set path; report their Q-limit residuals as diagnostics unless
  # a future slack-control mode explicitly handles them.
  status_text = pv_violations > 0 ? "FAIL" : (ref_violations > 0 ? "WARN" : "OK")

  println(io, "Final PV Q-limit active-set check: ", status_text == "FAIL" ? "FAIL" : "OK")
  @printf(io, "  active PV buses checked with finite Qmin/Qmax: %d\n", checked_pv)
  @printf(io, "  PV violations: %d (tolerance %.6g pu / %.6g MVAr)\n", pv_violations, tolerance_pu, tolerance_pu * net.baseMVA)
  println(io, "Final REF Q-limit diagnostic: ", ref_violations > 0 ? "WARN" : "OK")
  @printf(io, "  active REF buses checked with finite Qmin/Qmax: %d\n", checked_ref)
  @printf(io, "  REF violations: %d (tolerance %.6g pu / %.6g MVAr)\n", ref_violations, tolerance_pu, tolerance_pu * net.baseMVA)

  if isempty(rows)
    println(io, "Final PV/REF Q-limit check: OK")
  else
    if pv_violations > 0
      println(io, "Final PV/REF Q-limit check: FAIL")
    else
      println(io, "Final PV/REF Q-limit check: WARN (REF/slack violations are reported but do not fail the active PV set by default)")
    end
    @printf(io, "  active PV/REF buses checked with finite Qmin/Qmax: %d\n", checked_pv + checked_ref)
    @printf(io, "  remaining violations: %d (PV %d, REF %d)\n", length(rows), pv_violations, ref_violations)
    println(io, "  BUS_I │ type │ Qcalc pu │ Qcalc MVAr │ Qmin pu │ Qmin MVAr │ Qmax pu │ Qmax MVAr │ dir │ viol pu │ viol MVAr │ switched │ last it │ return band")
    row_limit = isnothing(max_console_rows) ? max_rows : max_console_rows
    shown = row_limit < 0 ? length(rows) : min(row_limit, length(rows))
    for row in Iterators.take(rows, shown)
      last_text = isnothing(row.last_it) ? "-" : string(row.last_it)
      band_text = row.in_return_band ? "inside" : "outside"
      type_text = row.type == :Slack ? "REF" : String(row.type)
      @printf(
        io,
        " %6d │ %4s │ %8.5f │ %10.3f │ %7.5f │ %9.3f │ %7.5f │ %9.3f │ %4s │ %7.5f │ %9.3f │ %8s │ %7s │ %s\n",
        row.busI,
        type_text,
        row.qcalc,
        row.qcalc * net.baseMVA,
        row.qmin,
        row.qmin * net.baseMVA,
        row.qmax,
        row.qmax * net.baseMVA,
        String(row.side),
        row.amount,
        row.amount * net.baseMVA,
        string(row.switched),
        last_text,
        band_text
      )
    end
    if shown < length(rows)
      @printf(io, "  (%d more violation rows omitted; increase max_rows/max_console_rows for full output)\n", length(rows) - shown)
    end
  end
  return (
    checked = checked_pv + checked_ref,
    checked_pv = checked_pv,
    checked_ref = checked_ref,
    violating = length(rows),
    pv_violations = pv_violations,
    ref_violations = ref_violations,
    violating_hyst = count(r -> !r.in_return_band, rows),
    switches = length(net.qLimitLog),
    status = Symbol(lowercase(status_text)),
  )
end

function _expand_ybus_for_isolated_nodes(Yred, n::Int, iso_nodes::Vector{Int})
  isempty(iso_nodes) && return Yred

  iso_mask = falses(n)
  for bus in iso_nodes
    if 1 <= bus <= n
      iso_mask[bus] = true
    end
  end

  active = Int[]
  for bus in eachindex(iso_mask)
    iso_mask[bus] || push!(active, bus)
  end

  size(Yred, 1) == length(active) || error("_expand_ybus_for_isolated_nodes: size mismatch between reduced Ybus and active buses.")
  size(Yred, 2) == length(active) || error("_expand_ybus_for_isolated_nodes: Ybus is not square in active-bus space.")

  Yfull = issparse(Yred) ? spzeros(ComplexF64, n, n) : zeros(ComplexF64, n, n)
  Yfull[active, active] = Yred
  return Yfull
end

@inline function _has_vset_adjust_config(ps::ProSumer)::Bool
  return !isnothing(ps.vset_adjust) || !(isnothing(ps.vstep_pu) && isnothing(ps.tap_steps_down) && isnothing(ps.tap_steps_up))
end

@inline function _bus_label(net::Net, bus::Int)::String
  return getCompName(net.nodeVec[bus].comp)
end

@inline function _resolve_vset_adjust_config(ps::ProSumer)::Union{Nothing,VoltageAdjustConfig}
  if !isnothing(ps.vset_adjust)
    return ps.vset_adjust
  end
  has_any = _has_vset_adjust_config(ps)
  has_any || return nothing
  all_defined = !isnothing(ps.vstep_pu) && !isnothing(ps.tap_steps_down) && !isnothing(ps.tap_steps_up)
  all_defined || return nothing
  return VoltageAdjustConfig(Float64(ps.vstep_pu), Int(ps.tap_steps_down), Int(ps.tap_steps_up))
end

function _try_adjust_vset_on_q_limit!(net::Net, bus::Int, side::Symbol, it::Int, controllers::Dict{Int,NamedTuple{(:prosumer_idx, :config),Tuple{Int,VoltageAdjustConfig}}}, base_vset::Vector{Float64}, Vset::Vector{Float64}, adjust_counter::Vector{Int}, qlimit_max_outer::Int, verbose::Int)::Bool
  cname = _bus_label(net, bus)
  if !haskey(controllers, bus)
    if verbose > 0
      @info "Bus $cname: no voltage adjustment controller -> fallback PV→PQ (it=$it)"
    end
    return false
  end

  ctrl = controllers[bus].config
  vm_base = base_vset[bus]
  vm_min = vm_base - ctrl.tap_steps_down * ctrl.vstep_pu
  vm_max = vm_base + ctrl.tap_steps_up * ctrl.vstep_pu
  vm_old = Vset[bus]
  vm_new = side == :max ? vm_old - ctrl.vstep_pu : vm_old + ctrl.vstep_pu
  can_step = (vm_new >= vm_min - 1e-12) && (vm_new <= vm_max + 1e-12)

  if can_step && adjust_counter[bus] < qlimit_max_outer
    vm_new = clamp(vm_new, vm_min, vm_max)
    Vset[bus] = vm_new
    net.nodeVec[bus]._vm_pu = vm_new
    adjust_counter[bus] += 1
    if verbose > 0
      event = side == :max ? "Qmax violated" : "Qmin violated"
      @info "Bus $cname: $event -> voltage adjusted from $vm_old to $vm_new (it=$it, step=$(adjust_counter[bus]))"
    end
    return true
  end

  if verbose > 0
    @info "Bus $cname: no further voltage steps (vm_min=$vm_min, vm_max=$vm_max, vm=$vm_old) -> fallback PV→PQ (it=$it)"
  end
  return false
end

"""
    runpf_rectangular!(net; maxiter=20, tol=1e-8, damp=0.2, verbose=0)

Run a complex-state Newton-Raphson power flow in rectangular coordinates on a Sparlectra network.

# Arguments
- `net::Net`: Network object containing bus, branch, and generation data
- `maxiter::Int=20`: Maximum number of Newton-Raphson iterations
- `tol::Float64=1e-8`: Convergence tolerance for maximum mismatch
- `damp::Float64=0.2`: Damping factor for Newton step (0 < damp ≤ 1)
- `verbose::Int=0`: Verbosity level (0=quiet, 1=basic info, 2=detailed)

# Returns
- `Tuple{Int, Int}`: (iterations_used, error_code)
  - error_code: 0=converged, 1=max_iterations_reached

# Details
This function implements a complete power flow solver using complex voltages in rectangular 
coordinates (Vr + jVi) as state variables, providing an alternative to conventional 
polar formulations.

## Mathematical Foundation
- **State vector**: x = [Vr(non-slack); Vi(non-slack)] ∈ ℝ^(2(n-1))
- **Complex power**: S = V .* conj(Y * V) where Y is the bus admittance matrix
- **PQ buses**: ΔP = Re(S_calc - S_spec), ΔQ = Im(S_calc - S_spec)
- **PV buses**: ΔP = Re(S_calc - S_spec), ΔV = |V| - V_set
- **Slack bus**: Voltage held constant throughout iterations

## Active Set Q-Limit Management
- **PV→PQ switching**: When reactive power demand violates generator Q-limits
- **Optional PQ→PV re-enable**: With hysteresis band and cooldown mechanisms
- **Robust handling**: Guards against inappropriate switching of non-generator buses

## Algorithm Steps
1. **Initialization**: Extract voltages, build Y-bus, classify bus types
2. **Power specification**: Build S = P + jQ from network loads/generation
3. **Iterative solution**: Newton-Raphson with mismatch function for PQ/PV constraints
4. **Q-limit enforcement**: Active-set management during iterations
5. **Result update**: Write final voltages and computed powers back to network

## Network Integration
- **Input**: Uses `net.nodeVec` for bus data, `net.baseMVA` for per-unit conversion
- **Output**: Updates `node._vm_pu`, `node._va_deg`, `node._pƩGen`, `node._qƩGen`
- **Q-limits**: Integrates with `net.qLimitEvents` and `net.qLimitLog` for tracking
- **Compatibility**: Maintains same interface as `runpf!()` for easy substitution

# See Also
- `runpf_rectangular!(net, maxIte, tolerance, verbose; ...)`: compatibility overload matching `runpf!()` signature
- `mismatch_rectangular()`: Core mismatch function for PQ/PV constraints
- `build_rectangular_jacobian_pq_pv()`: Analytic Jacobian construction
"""

function runpf_rectangular!(
  net::Net;
  maxiter::Int = 20,
  tol::Float64 = 1e-8,
  damp::Float64 = 0.2,
  verbose::Int = 0,
  autodamp::Bool = false,
  autodamp_min::Float64 = 0.05,
  opt_flatstart::Bool = net.flatstart,
  pv_table_rows::Int = 30,
  lock_pv_to_pq_buses::AbstractVector{Int} = Int[],
  qlimit_mode::Symbol = :switch_to_pq,
  qlimit_max_outer::Int = 30,
  start_projection::Bool = false,
  start_projection_try_dc_start::Bool = true,
  start_projection_try_blend_scan::Bool = true,
  start_projection_branch_guard::Bool = true,
  start_projection_measure_candidates::Bool = true,
  start_projection_accept_unmeasured_dc_start::Bool = false,
  start_projection_blend_lambdas::AbstractVector{<:Real} = [0.25, 0.5, 0.75],
  start_projection_dc_angle_limit_deg::Float64 = 60.0,
  qlimit_start_iter::Int = 2,
  qlimit_start_mode::Symbol = :iteration,
  qlimit_auto_q_delta_pu::Float64 = 1e-4,
  qlimit_trace_buses::AbstractVector{Int} = Int[],
  qlimit_lock_reason::Symbol = :manual,
  qlimit_guard::Bool = false,
  qlimit_guard_min_q_range_pu::Float64 = 1e-4,
  qlimit_guard_zero_range_mode::Symbol = :lock_pq,
  qlimit_guard_narrow_range_mode::Symbol = :prefer_pq,
  qlimit_guard_log::Bool = true,
  qlimit_guard_max_switches::Int = 10,
  qlimit_guard_accept_bounded_violations::Bool = false,
  qlimit_guard_max_remaining_violations::Int = 0,
  qlimit_guard_freeze_after_repeated_switching::Bool = true,
  qlimit_guard_violation_mode::Symbol = :delayed_switch,
  qlimit_guard_violation_threshold_pu::Float64 = 1e-4,
  wrong_branch_detection::Symbol = :warn,
  wrong_branch_rescue::Bool = false,
  wrong_branch_min_vm_pu::Float64 = 0.70,
  wrong_branch_max_vm_pu::Float64 = 1.30,
  wrong_branch_max_angle_spread_deg::Float64 = 180.0,
  wrong_branch_max_branch_angle_deg::Float64 = 90.0,
  wrong_branch_min_low_vm_count::Int = 1,
  wrong_branch_rescue_max_attempts::Int = 2,
  performance_profile = nothing,
  rectangular_workspace_reuse::Bool = true,
  rectangular_preallocate_workspace::Symbol = :auto,
  rectangular_workspace_min_buses::Int = 1000,
)
  _validate_rectangular_powerflow_options(method = :rectangular, sparse = true)
  if verbose > 1
    @info "Running complex rectangular NR power flow..."
  end

  nodes = net.nodeVec
  n = length(nodes)
  Sbase = net.baseMVA
  Yred = _perf_profile_time!(performance_profile, :ybus_assembly) do
    createYBUS(net = net, sparse = true, printYBUS = (verbose > 1))
  end
  Ybus = _perf_profile_time!(performance_profile, :ybus_expand_isolated) do
    (size(Yred, 1) == n) ? Yred : _expand_ybus_for_isolated_nodes(Yred, n, net.isoNodes)
  end

  # 1) Initial complex voltages V0 and slack index
  V0, slack_idx = _perf_profile_time!(performance_profile, :solver_initial_voltage) do
    initialVrect(net; flatstart = opt_flatstart)
  end

  # 2) Specified complex power injections S (p.u.). For Q(U)/P(U) controllers
  # this vector becomes state-dependent and is re-evaluated per Newton iteration.
  S = _perf_profile_time!(performance_profile, :solver_initial_injections) do
    buildComplexSVec(net)
  end
  dPinj_dVm = zeros(Float64, n)
  dQinj_dVm = zeros(Float64, n)
  has_vdep_control = has_voltage_dependent_control(net)

  # 3) Bus types from Node data, and PV setpoints from regulating prosumers.
  # Node voltages may be temporary start guesses for MATPOWER flat-start modes.
  bus_types = Vector{Symbol}(undef, n)
  Vset = _perf_profile_time!(performance_profile, :solver_voltage_setpoint_lookup) do
    _bus_voltage_setpoints_from_prosumers(net; performance_profile = performance_profile)
  end

  _perf_profile_time!(performance_profile, :solver_bus_type_scan) do
    @inbounds for (k, node) in enumerate(nodes)
      BusType = getNodeType(node)
      if BusType == Slack
        bus_types[k] = :Slack
      elseif BusType == PV
        bus_types[k] = :PV
      elseif BusType == PQ
        bus_types[k] = :PQ
      elseif BusType == Isolated
        # Keep isolated buses in the rectangular state vector as neutral PQ rows.
        # Their injections are forced to zero so they do not affect the solved grid.
        bus_types[k] = :PQ
        S[k] = 0.0 + 0.0im
      else
        error("runpf_rectangular!: unsupported bus type at bus $k, given: $(BusType)")
      end
    end
  end

  # The slack/reference voltage is fixed rather than solved by a residual row.
  # Keep its magnitude at the regulating prosumer setpoint even if a MATPOWER
  # flat-start mode temporarily changed node._vm_pu as an initial guess.
  _perf_profile_time!(performance_profile, :solver_slack_voltage_fix) do
    V0[slack_idx] = ComplexF64(Vset[slack_idx] * cos(angle(V0[slack_idx])), Vset[slack_idx] * sin(angle(V0[slack_idx])))
  end

  V0 = _perf_profile_time!(performance_profile, :start_projection) do
    project_rectangular_start(
      Ybus,
      V0,
      S,
      bus_types,
      Vset,
      slack_idx;
      enabled = start_projection,
      try_dc_start = start_projection_try_dc_start,
      try_blend_scan = start_projection_try_blend_scan,
      branch_guard = start_projection_branch_guard,
      measure_candidates = start_projection_measure_candidates,
      accept_unmeasured_dc_start = start_projection_accept_unmeasured_dc_start,
      blend_lambdas = start_projection_blend_lambdas,
      dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
      verbose = verbose,
      performance_profile = performance_profile,
    )
  end

  # 4) Q-limit data 
  qmin_pu, qmax_pu = _perf_profile_time!(performance_profile, :solver_qlimit_extraction) do
    getQLimits_pu(net)
  end
  # Start fresh each PF run before guard pre-processing records locked buses.
  _perf_profile_time!(performance_profile, :solver_qlimit_log_reset) do
    resetQLimitLog!(net)
  end
  if verbose > 1
    printPVQLimitsTable(net; max_rows = typemax(Int))
  elseif verbose > 0
    printPVQLimitsTable(net; max_rows = pv_table_rows)
  end

  guarded_qlimit_buses = Int[]
  if qlimit_guard
    guarded_qlimit_buses = _perf_profile_time!(performance_profile, :qlimit_guard_preprocess) do
      _apply_qlimit_guard_to_rectangular_active_set!(net, bus_types, S, build_qload_pu(net), qmin_pu, qmax_pu; min_q_range_pu = qlimit_guard_min_q_range_pu, zero_range_mode = qlimit_guard_zero_range_mode, narrow_range_mode = qlimit_guard_narrow_range_mode, log = qlimit_guard_log, verbose = verbose)
    end
  end

  # --- Active-set bookkeeping (rectangular solver) ------------------------
  nb = n  # number of buses

  # PV origin mask (guards for PQ->PV re-enable)
  pv_orig_mask = _perf_profile_time!(performance_profile, :solver_active_set_origin_mask) do
    mask = falses(nb)
    @inbounds for k in eachindex(mask)
      mask[k] = (bus_types[k] == :PV)
    end
    mask
  end

  cooldown_iters = hasfield(typeof(net), :cooldown_iters) ? net.cooldown_iters : 0
  q_hyst_pu      = hasfield(typeof(net), :q_hyst_pu) ? net.q_hyst_pu : 0.0
  allow_reenable = (cooldown_iters > 0) || (q_hyst_pu > 0.0)
  qlimit_mode in (:switch_to_pq, :adjust_vset) || error("Unsupported qlimit_mode=$(qlimit_mode). Supported: :switch_to_pq, :adjust_vset.")
  qlimit_max_outer > 0 || error("qlimit_max_outer must be > 0 (got $(qlimit_max_outer)).")
  qlimit_start_iter > 0 || error("qlimit_start_iter must be > 0 (got $(qlimit_start_iter)).")
  qlimit_start_mode in (:iteration, :auto_q_delta, :iteration_or_auto) || error("Unsupported qlimit_start_mode=$(qlimit_start_mode). Supported: :iteration, :auto_q_delta, :iteration_or_auto.")
  qlimit_auto_q_delta_pu >= 0.0 || error("qlimit_auto_q_delta_pu must be >= 0 (got $(qlimit_auto_q_delta_pu)).")
  qlimit_guard_violation_mode in (:delayed_switch, :lock_pq, :ignore) || error("Unsupported qlimit_guard_violation_mode=$(qlimit_guard_violation_mode). Supported: :delayed_switch, :lock_pq, :ignore.")
  qlimit_guard_violation_threshold_pu >= 0.0 || error("qlimit_guard_violation_threshold_pu must be >= 0 (got $(qlimit_guard_violation_threshold_pu)).")

  controllers, base_vset, adjust_counter, qlimit_trace_internal = _perf_profile_time!(performance_profile, :solver_active_set_setup) do
    controllers_ = qlimit_mode == :adjust_vset ? _build_vset_adjust_controllers(net) : Dict{Int,NamedTuple{(:prosumer_idx, :config),Tuple{Int,VoltageAdjustConfig}}}()
    (controllers_, copy(Vset), zeros(Int, nb), _resolve_qlimit_trace_buses(net, qlimit_trace_buses))
  end
  qlimit_trace_enabled = !isempty(qlimit_trace_internal)
  if qlimit_trace_enabled
    missing = setdiff(collect(qlimit_trace_buses), [_qlimit_original_bus_id(net, bus) for bus in qlimit_trace_internal])
    isempty(missing) || @warn "qlimit_trace_buses entries not found in network" missing = missing
    println("Q-limit trace enabled for BUS_I values: ", [_qlimit_original_bus_id(net, bus) for bus in qlimit_trace_internal])
  end

  # 5) NR-Loop
  V = copy(V0)
  history = Float64[]
  converged = false
  iters = 0
  rejection_reason = :nr_mismatch_not_converged
  rectangular_workspace_reason = :disabled
  rectangular_workspace_preallocated = false
  if !rectangular_workspace_reuse || rectangular_preallocate_workspace == :off
    rectangular_workspace_reason = rectangular_workspace_reuse ? :off : :reuse_disabled
  elseif rectangular_preallocate_workspace == :on
    rectangular_workspace_preallocated = true
    rectangular_workspace_reason = :forced_on
  elseif rectangular_preallocate_workspace == :auto
    rectangular_workspace_preallocated = nb >= rectangular_workspace_min_buses
    rectangular_workspace_reason = rectangular_workspace_preallocated ? :auto_threshold : :auto_below_threshold
  else
    error("Unsupported rectangular_preallocate_workspace=$(rectangular_preallocate_workspace). Supported: :off, :on, :auto.")
  end
  workspace = RectangularIterationWorkspace(nb)
  if performance_profile !== nothing
    performance_profile[:rectangular_workspace_reuse] = rectangular_workspace_reuse
    performance_profile[:rectangular_workspace_preallocated] = rectangular_workspace_preallocated
    performance_profile[:rectangular_workspace_reason] = rectangular_workspace_reason
    performance_profile[:rectangular_workspace_nbus] = nb
    performance_profile[:rectangular_workspace_nstate] = 2 * max(nb - 1, 0)
  end

  if verbose > 1
    @info "Starting rectangular complex NR power flow..."
    @info "Initial complex voltages V0:" V0
    @info "Slack bus index:" slack_idx
    @info "maxiter = $maxiter, tol = $tol, damp = $damp, autodamp = $autodamp, autodamp_min = $autodamp_min, start_projection = $start_projection"
  end

  qlimit_active_set_changes = 0
  qlimit_reenable_events = 0
  for it = 1:maxiter
    iters = it

    if has_vdep_control
      S, dPinj_dVm, dQinj_dVm = _perf_profile_time!(performance_profile, :iteration_controlled_injections) do
        buildControlledSVec(net, V)
      end
    end

    # Mismatch with current bus_types and (possibly) voltage-dependent S.
    F = _perf_profile_time!(performance_profile, :iteration_mismatch) do
      mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
    end
    max_mis = maximum(abs.(F))
    push!(history, max_mis)
    if rectangular_workspace_preallocated
      resize!(workspace.rhs_vector, length(F))
      copyto!(workspace.rhs_vector, F)
    end

    (verbose > 1) && @debug "Rectangular NR iteration" iter = it max_mismatch = max_mis

    # --- Q-Limit Active Set: PV -> PQ, optional PQ -> PV (rectangular) ------
    changed   = false
    reenabled = false

    _perf_profile_time!(performance_profile, :iteration_qload) do
      copyto!(workspace.qload_pu, build_qload_pu(net))
    end

    Scalc_pu = _perf_profile_time!(performance_profile, :iteration_calc_injections) do
      calc_injections(Ybus, V)
    end
    current_pv_qreq_pu = _perf_profile_time!(performance_profile, :iteration_qreq_vector) do
      fill!(workspace.current_pv_qreq_pu, NaN)
      @inbounds for bus in eachindex(workspace.current_pv_qreq_pu)
        if bus_types[bus] == :PV
          workspace.current_pv_qreq_pu[bus] = imag(Scalc_pu[bus]) + workspace.qload_pu[bus]
        end
      end
      workspace.current_pv_qreq_pu
    end

    qlimit_iter_ready, qlimit_auto_ready, qlimit_ready, converged_this_iter, qlimit_check_active = _perf_profile_time!(performance_profile, :iteration_control_bookkeeping) do
      qlimit_iter_ready_ = it >= qlimit_start_iter
      qlimit_auto_ready_ = false
      if qlimit_start_mode in (:auto_q_delta, :iteration_or_auto)
        max_q_delta = 0.0
        compared = false
        @inbounds for bus in eachindex(current_pv_qreq_pu)
          if isfinite(current_pv_qreq_pu[bus]) && isfinite(workspace.prev_pv_qreq_pu[bus])
            max_q_delta = max(max_q_delta, abs(current_pv_qreq_pu[bus] - workspace.prev_pv_qreq_pu[bus]))
            compared = true
          end
        end
        qlimit_auto_ready_ = compared && (max_q_delta <= qlimit_auto_q_delta_pu)
      end
      qlimit_ready_ = qlimit_start_mode == :iteration ? qlimit_iter_ready_ : qlimit_start_mode == :auto_q_delta ? qlimit_auto_ready_ : (qlimit_iter_ready_ || qlimit_auto_ready_)
      converged_this_iter_ = max_mis <= tol
      violation_guard_active = qlimit_guard_violation_mode == :lock_pq
      qlimit_check_active_ = qlimit_ready_ && (!converged_this_iter_ || violation_guard_active)
      (qlimit_iter_ready_, qlimit_auto_ready_, qlimit_ready_, converged_this_iter_, qlimit_check_active_)
    end

    if qlimit_trace_enabled
      _perf_profile_time!(performance_profile, :iteration_qlimit_trace) do
        lock_mask_trace = falses(nb)
        for bus in lock_pv_to_pq_buses
          if 1 <= bus <= nb
            lock_mask_trace[bus] = true
          end
        end
        for bus in qlimit_trace_internal
          bus <= nb || continue
          qreq = bus_types[bus] in (:PV, :Slack) ? imag(Scalc_pu[bus]) + workspace.qload_pu[bus] : NaN
          _print_rectangular_qlimit_trace(
            stdout,
            net,
            it,
            bus,
            bus_types[bus],
            V,
            Vset,
            qreq,
            qmin_pu,
            qmax_pu;
            q_hyst_pu = q_hyst_pu,
            qlimit_start_iter = qlimit_start_iter,
            qlimit_start_mode = qlimit_start_mode,
            qlimit_iter_ready = qlimit_iter_ready,
            qlimit_auto_ready = qlimit_auto_ready,
            qlimit_ready = qlimit_ready,
            qlimit_check_active = qlimit_check_active,
            converged_before_switching = converged_this_iter,
            qlimit_mode = qlimit_mode,
            lock_mask = lock_mask_trace,
            qlimit_lock_reason = qlimit_lock_reason,
          )
        end
      end
    end

    if qlimit_check_active
      changed, reenabled = _perf_profile_time!(performance_profile, :iteration_qlimit_active_set) do
        active_set_q_limits!(
          net,
          it,
          nb;
          qmin_pu = qmin_pu,
          qmax_pu = qmax_pu,
          pv_orig_mask = pv_orig_mask,
          allow_reenable = qlimit_mode == :switch_to_pq ? allow_reenable : false,
          q_hyst_pu = q_hyst_pu,
          cooldown_iters = cooldown_iters,
          lock_pv_to_pq_buses = lock_pv_to_pq_buses,
          on_violation! = qlimit_mode == :adjust_vset ? ((bus, qreq, side, qclamp) -> _try_adjust_vset_on_q_limit!(net, bus, side, it, controllers, base_vset, Vset, adjust_counter, qlimit_max_outer, verbose)) : nothing,
          verbose = verbose,
          get_qreq_pu = bus -> begin
            (bus_types[bus] == :Slack) && return 0.0
            return imag(Scalc_pu[bus]) + workspace.qload_pu[bus]
          end,
          is_pv = bus -> (bus_types[bus] == :PV),
          make_pq! = (bus, qclamp_gen_pu, side) -> begin
            bus_types[bus] = :PQ
            qinj_pu = qclamp_gen_pu - workspace.qload_pu[bus]
            S[bus] = ComplexF64(real(S[bus]), qinj_pu)
            net.nodeVec[bus]._qƩGen = qclamp_gen_pu * net.baseMVA
          end,
          make_pv! = (bus) -> begin
            bus_types[bus] = :PV
          end,
          qlimit_guard_max_switches = qlimit_guard_max_switches,
          qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching,
          qlimit_guard_violation_mode = qlimit_guard_violation_mode,
          qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu,
          lock_mask_buffer = workspace.lock_mask,
        )
      end

      # If bus_types/spec changed, mismatch definition changed (ΔQ ↔ ΔV) => rebuild F
      if changed || reenabled
        changed && (qlimit_active_set_changes += 1)
        reenabled && (qlimit_reenable_events += 1)
        F = _perf_profile_time!(performance_profile, :iteration_mismatch_after_qlimit) do
          mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
        end
        max_mis = maximum(abs.(F))
        history[end] = max_mis  # optional: overwrite last stored value for this iteration
      end
    end
    if converged_this_iter && !(changed || reenabled)
      converged = true
      rejection_reason = :none
      break
    end
    copyto!(workspace.prev_pv_qreq_pu, current_pv_qreq_pu)
    # --- Newton step (FD or analytic) -----------------------------------
    try
      V = _perf_profile_time!(performance_profile, :iteration_newton_step) do
        complex_newton_step_rectangular(Ybus, V, S; slack_idx = slack_idx, damp = damp, autodamp = autodamp, autodamp_min = autodamp_min, bus_types = bus_types, Vset = Vset, dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm, performance_profile = performance_profile)
      end
    catch step_error
      if _is_rectangular_linear_step_failure(step_error)
        verbose > 0 && @warn "Rectangular Newton step failed because the linear Jacobian solve was singular; returning non-convergence." iteration = it max_mismatch = max_mis exception = (typeof(step_error), sprint(showerror, step_error))
        converged = false
        rejection_reason = :singular_newton_step
        break
      end
      rethrow(step_error)
    end

    # Keep slack voltage fixed (safety belt)
    _perf_profile_time!(performance_profile, :iteration_state_update) do
      V[slack_idx] = V0[slack_idx]
      _perf_profile_push_iteration!(performance_profile, (iteration = it, max_mismatch = max_mis, qlimit_changed = changed, qlimit_reenabled = reenabled))
    end
  end

  # 6) Update voltages back to network
  # --- mirror bus_types back into Net/node types (PV->PQ switching) ---
  _perf_profile_time!(performance_profile, :solver_final_active_set_sync) do
    _sync_rectangular_bus_types_to_net!(net, bus_types)
  end
  numerical_converged = converged
  final_pv_voltage_residual = _perf_profile_time!(performance_profile, :solver_final_status_checks) do
    _max_rectangular_pv_voltage_residual(V, Vset, bus_types, net.qLimitEvents)
  end
  if numerical_converged && final_pv_voltage_residual > tol
    verbose > 0 && @warn "Rectangular NR convergence rejected because active PV voltage setpoint residual exceeds tolerance." max_pv_voltage_residual = final_pv_voltage_residual tolerance = tol
    converged = false
    rejection_reason = :active_pv_voltage_residual
  end

  _perf_profile_time!(performance_profile, :solver_voltage_writeback) do
    _finalize_rectangular_voltage_writeback!(net, V)
  end

  # 7) Compute bus injections from final voltages  
  Sbus_pu, Sbus_MVA = _perf_profile_time!(performance_profile, :solver_final_injection_vectors) do
    _compute_rectangular_final_injections(Ybus, V, Sbase)
  end

  @debug "Final Voltages Mag = " [abs.(V)...]
  @debug "Final Voltages Ang = " [angle.(V) .* (180.0 / π)...]

  _perf_profile_time!(performance_profile, :solver_result_bus_writeback) do
    _write_rectangular_bus_power_results!(net, nodes, Sbus_MVA)
  end

  # 8) Update total bus power (sum of complex injections in p.u.)
  p, q = _write_rectangular_total_bus_power!(net, Sbus_pu, Sbase, verbose, performance_profile)

  qlimit_summary = nothing
  branch_quality = _wrong_branch_not_checked_result()
  wrong_branch_rescue_attempted = false
  wrong_branch_rescue_reason = :disabled
  if numerical_converged
    qlimit_status = _perf_profile_time!(performance_profile, :solver_final_qlimit_summary) do
      _finalize_rectangular_qlimit_summary(
        net,
        V,
        Sbus_pu,
        bus_types,
        qmin_pu,
        qmax_pu,
        converged,
        rejection_reason;
        verbose = verbose,
        qlimit_trace_enabled = qlimit_trace_enabled,
        q_hyst_pu = q_hyst_pu,
        tol = tol,
        pv_table_rows = pv_table_rows,
        qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations,
        qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations,
      )
    end
    qlimit_summary = qlimit_status.qlimit_summary
    converged = qlimit_status.converged
    rejection_reason = qlimit_status.rejection_reason
    wrong_branch_status = _finalize_rectangular_wrong_branch_diagnostics(
      V,
      bus_types,
      Vset,
      slack_idx,
      converged,
      rejection_reason;
      wrong_branch_detection = wrong_branch_detection,
      wrong_branch_min_vm_pu = wrong_branch_min_vm_pu,
      wrong_branch_max_vm_pu = wrong_branch_max_vm_pu,
      wrong_branch_max_angle_spread_deg = wrong_branch_max_angle_spread_deg,
      wrong_branch_max_branch_angle_deg = wrong_branch_max_branch_angle_deg,
      wrong_branch_min_low_vm_count = wrong_branch_min_low_vm_count,
      net = net,
    )
    branch_quality = wrong_branch_status.branch_quality
    converged = wrong_branch_status.converged
    rejection_reason = wrong_branch_status.rejection_reason
    wrong_branch_rescue_attempted = wrong_branch_status.wrong_branch_rescue_attempted
    wrong_branch_rescue_reason = wrong_branch_status.wrong_branch_rescue_reason
  end

  switch_counts, oscillating_buses, max_switching_exceeded, q_limit_active_set_ok, converged, rejection_reason, final_reason, final_status, status = _perf_profile_time!(performance_profile, :solver_status_bookkeeping) do
    switch_counts_ = qlimit_switch_counts(net)
    oscillating_buses_ = count(>=(max(qlimit_guard_max_switches, 1)), values(switch_counts_))
    max_switching_exceeded_ = qlimit_guard_freeze_after_repeated_switching && oscillating_buses_ > 0
    q_limit_active_set_ok_ = numerical_converged && final_pv_voltage_residual <= tol && (isnothing(qlimit_summary) || qlimit_summary.pv_violations == 0 || (qlimit_guard_accept_bounded_violations && qlimit_summary.pv_violations <= qlimit_guard_max_remaining_violations)) && !max_switching_exceeded_
    converged_ = converged
    rejection_reason_ = rejection_reason
    if numerical_converged && max_switching_exceeded_ && !q_limit_active_set_ok_
      rejection_reason_ = :max_switching_exceeded
      converged_ = false
    end
    status_build_ = _build_rectangular_final_status(
      net,
      numerical_converged,
      q_limit_active_set_ok_,
      converged_,
      rejection_reason_,
      qlimit_summary,
      final_pv_voltage_residual,
      history,
      qlimit_active_set_changes,
      qlimit_reenable_events,
      oscillating_buses_,
      guarded_qlimit_buses,
      branch_quality,
      wrong_branch_detection,
      wrong_branch_rescue_attempted,
      wrong_branch_rescue_reason,
    )
    final_reason_ = status_build_.final_reason
    final_status_ = status_build_.final_status
    status_ = _store_and_print_rectangular_final_status!(net, status_build_.status, verbose)
    (switch_counts_, oscillating_buses_, max_switching_exceeded_, q_limit_active_set_ok_, converged_, rejection_reason_, final_reason_, final_status_, status_)
  end
  if verbose > 0
    _print_qlimit_active_set_summary(stdout, status)
    println(stdout, "Wrong-branch check:")
    @printf(stdout, "  status           = %s\n", uppercase(String(status.branch_quality_status)))
    @printf(stdout, "  detection_mode   = %s\n", String(status.wrong_branch_detection))
    @printf(stdout, "  reason           = %s\n", String(status.branch_quality_reason))
    @printf(stdout, "  min_vm_pu        = %.6f\n", status.branch_quality_metrics.min_vm_pu)
    @printf(stdout, "  max_vm_pu        = %.6f\n", status.branch_quality_metrics.max_vm_pu)
    @printf(stdout, "  low_vm_count     = %d\n", status.branch_quality_metrics.low_vm_count)
    @printf(stdout, "  angle_spread_deg = %.6f\n", status.branch_quality_metrics.angle_spread_deg)
    @printf(stdout, "  max_branch_angle = %.6f\n", status.branch_quality_metrics.max_branch_angle_deg)
    @printf(stdout, "  violation_count  = %d\n", status.branch_quality_metrics.branch_angle_violation_count)
    !isnothing(status.branch_quality_metrics.worst_branch) && @printf(stdout, "  worst_branch     = %s\n", string(status.branch_quality_metrics.worst_branch))
    @printf(stdout, "  rescue_requested = %s\n", string(wrong_branch_rescue || wrong_branch_detection == :rescue))
    @printf(stdout, "  rescue_attempted = %s\n", string(status.wrong_branch_rescue_attempted))
    @printf(stdout, "  rescue_reason    = %s\n", String(status.wrong_branch_rescue_reason))
  end

  return iters, converged ? 0 : 1
end

"""
    runpf!(net, maxIte, tolerance=1e-6, verbose=0; method=:rectangular)

Unified AC power flow interface.

Arguments:
- `net::Net`: network
- `maxIte::Int`: maximum iterations
- `tolerance::Float64`: mismatch tolerance
- `verbose::Int`: verbosity level
- `method::Symbol`: must be `:rectangular`
- `autodamp::Bool`: enable residual-based backtracking for rectangular Newton steps
- `autodamp_min::Float64`: minimum automatic damping factor when `autodamp = true`
- `qlimit_start_iter::Int`: first Newton iteration where PV→PQ Q-limit switching may run in `:iteration` mode
- `qlimit_start_mode::Symbol`: `:iteration`, `:auto_q_delta`, or `:iteration_or_auto` start criterion for PV→PQ switching
- `qlimit_auto_q_delta_pu::Float64`: PV reactive-power request change threshold for automatic switching start

Notes:
- Link-flow recovery (`calcLinkFlowsKCL!`) is method-agnostic and uses solved PF results.
- If active-link merges create internal isolated buses, the rectangular sparse
  solver remains the only supported PF path; there is no polar fallback.

Returns:
    (iterations::Int, status::Int)
where `status == 0` indicates convergence.
"""
function runpf_rectangular!(
  net::Net,
  maxIte::Int,
  tolerance::Float64 = 1e-6,
  verbose::Int = 0;
  damp = 1.0,
  autodamp::Bool = false,
  autodamp_min::Float64 = 0.05,
  opt_flatstart::Bool = net.flatstart,
  pv_table_rows::Int = 30,
  validate_limits_after_pf::Bool = false,
  q_limit_violation_headroom::Float64 = 0.0,
  lock_pv_to_pq_buses::AbstractVector{Int} = Int[],
  qlimit_mode::Symbol = :switch_to_pq,
  qlimit_max_outer::Int = 30,
  start_projection::Bool = false,
  start_projection_try_dc_start::Bool = true,
  start_projection_try_blend_scan::Bool = true,
  start_projection_branch_guard::Bool = true,
  start_projection_measure_candidates::Bool = true,
  start_projection_accept_unmeasured_dc_start::Bool = false,
  start_projection_blend_lambdas::AbstractVector{<:Real} = [0.25, 0.5, 0.75],
  start_projection_dc_angle_limit_deg::Float64 = 60.0,
  qlimit_start_iter::Int = 2,
  qlimit_start_mode::Symbol = :iteration,
  qlimit_auto_q_delta_pu::Float64 = 1e-4,
  qlimit_trace_buses::AbstractVector{Int} = Int[],
  qlimit_lock_reason::Symbol = :manual,
  qlimit_guard::Bool = false,
  qlimit_guard_min_q_range_pu::Float64 = 1e-4,
  qlimit_guard_zero_range_mode::Symbol = :lock_pq,
  qlimit_guard_narrow_range_mode::Symbol = :prefer_pq,
  qlimit_guard_log::Bool = true,
  qlimit_guard_max_switches::Int = 10,
  qlimit_guard_accept_bounded_violations::Bool = false,
  qlimit_guard_max_remaining_violations::Int = 0,
  qlimit_guard_freeze_after_repeated_switching::Bool = true,
  qlimit_guard_violation_mode::Symbol = :delayed_switch,
  qlimit_guard_violation_threshold_pu::Float64 = 1e-4,
  wrong_branch_detection::Symbol = :warn,
  wrong_branch_rescue::Bool = false,
  wrong_branch_min_vm_pu::Float64 = 0.70,
  wrong_branch_max_vm_pu::Float64 = 1.30,
  wrong_branch_max_angle_spread_deg::Float64 = 180.0,
  wrong_branch_max_branch_angle_deg::Float64 = 90.0,
  wrong_branch_min_low_vm_count::Int = 1,
  wrong_branch_rescue_max_attempts::Int = 2,
  rectangular_workspace_reuse::Bool = true,
  rectangular_preallocate_workspace::Symbol = :auto,
  rectangular_workspace_min_buses::Int = 1000,
  performance_profile = nothing,
)
  iters, erg = runpf_rectangular!(
    net;
    maxiter = maxIte,
    tol = tolerance,
    damp = damp,
    autodamp = autodamp,
    autodamp_min = autodamp_min,
    verbose = verbose,
    opt_flatstart = opt_flatstart,
    pv_table_rows = pv_table_rows,
    lock_pv_to_pq_buses = lock_pv_to_pq_buses,
    qlimit_mode = qlimit_mode,
    qlimit_max_outer = qlimit_max_outer,
    start_projection = start_projection,
    start_projection_try_dc_start = start_projection_try_dc_start,
    start_projection_try_blend_scan = start_projection_try_blend_scan,
    start_projection_branch_guard = start_projection_branch_guard,
    start_projection_measure_candidates = start_projection_measure_candidates,
    start_projection_accept_unmeasured_dc_start = start_projection_accept_unmeasured_dc_start,
    start_projection_blend_lambdas = start_projection_blend_lambdas,
    start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
    qlimit_start_iter = qlimit_start_iter,
    qlimit_start_mode = qlimit_start_mode,
    qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu,
    qlimit_trace_buses = qlimit_trace_buses,
    qlimit_lock_reason = qlimit_lock_reason,
    qlimit_guard = qlimit_guard,
    qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu,
    qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode,
    qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode,
    qlimit_guard_log = qlimit_guard_log,
    qlimit_guard_max_switches = qlimit_guard_max_switches,
    qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations,
    qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations,
    qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching,
    qlimit_guard_violation_mode = qlimit_guard_violation_mode,
    qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu,
    wrong_branch_detection = wrong_branch_detection,
    wrong_branch_rescue = wrong_branch_rescue,
    wrong_branch_min_vm_pu = wrong_branch_min_vm_pu,
    wrong_branch_max_vm_pu = wrong_branch_max_vm_pu,
    wrong_branch_max_angle_spread_deg = wrong_branch_max_angle_spread_deg,
    wrong_branch_max_branch_angle_deg = wrong_branch_max_branch_angle_deg,
    wrong_branch_min_low_vm_count = wrong_branch_min_low_vm_count,
    wrong_branch_rescue_max_attempts = wrong_branch_rescue_max_attempts,
    rectangular_workspace_reuse = rectangular_workspace_reuse,
    rectangular_preallocate_workspace = rectangular_preallocate_workspace,
    rectangular_workspace_min_buses = rectangular_workspace_min_buses,
    performance_profile = performance_profile,
  )
  return iters, erg
end

function _runpf_with_config!(net::Net, config::PowerFlowConfig; verbose::Int = 0, damp = 1.0, pv_table_rows::Int = 30, validate_limits_after_pf::Bool = false, q_limit_violation_headroom::Float64 = 0.0, qlimit_lock_reason::Symbol = :manual, performance_profile = nothing)
  start = config.start_mode
  qlim = config.qlimits
  return runpf!(
    net,
    config.max_iter,
    config.tol,
    verbose;
    method = config.method,
    opt_flatstart = start.flatstart,
    damp = damp,
    autodamp = config.autodamp,
    autodamp_min = config.autodamp_min,
    wrong_branch_detection = config.wrong_branch_detection,
    wrong_branch_rescue = config.wrong_branch_rescue,
    wrong_branch_min_vm_pu = config.wrong_branch_min_vm_pu,
    wrong_branch_max_vm_pu = config.wrong_branch_max_vm_pu,
    wrong_branch_max_angle_spread_deg = config.wrong_branch_max_angle_spread_deg,
    wrong_branch_max_branch_angle_deg = config.wrong_branch_max_branch_angle_deg,
    wrong_branch_min_low_vm_count = config.wrong_branch_min_low_vm_count,
    wrong_branch_rescue_max_attempts = config.wrong_branch_rescue_max_attempts,
    pv_table_rows = pv_table_rows,
    validate_limits_after_pf = validate_limits_after_pf,
    q_limit_violation_headroom = q_limit_violation_headroom,
    lock_pv_to_pq_buses = qlim.lock_pv_to_pq_buses,
    start_projection = start.start_projection,
    start_projection_try_dc_start = start.try_dc_start,
    start_projection_try_blend_scan = start.try_blend_scan,
    start_projection_branch_guard = start.branch_guard,
    start_projection_measure_candidates = start.measure_candidates,
    start_projection_accept_unmeasured_dc_start = start.accept_unmeasured_dc_start,
    start_projection_blend_lambdas = start.blend_lambdas,
    start_projection_dc_angle_limit_deg = start.dc_angle_limit_deg,
    qlimit_start_iter = qlim.start_iter,
    qlimit_start_mode = qlim.start_mode,
    qlimit_auto_q_delta_pu = qlim.auto_q_delta_pu,
    qlimit_trace_buses = qlim.trace_buses,
    qlimit_lock_reason = qlimit_lock_reason,
    qlimit_guard = qlim.guard,
    qlimit_guard_min_q_range_pu = qlim.guard_min_q_range_pu,
    qlimit_guard_zero_range_mode = qlim.guard_zero_range_mode,
    qlimit_guard_narrow_range_mode = qlim.guard_narrow_range_mode,
    qlimit_guard_log = qlim.guard_log,
    qlimit_guard_max_switches = qlim.guard_max_switches,
    qlimit_guard_accept_bounded_violations = qlim.guard_accept_bounded_violations,
    qlimit_guard_max_remaining_violations = qlim.guard_max_remaining_violations,
    qlimit_guard_freeze_after_repeated_switching = qlim.guard_freeze_after_repeated_switching,
    qlimit_guard_violation_mode = qlim.guard_violation_mode,
    qlimit_guard_violation_threshold_pu = qlim.guard_violation_threshold_pu,
    rectangular_workspace_reuse = config.rectangular_workspace_reuse,
    rectangular_preallocate_workspace = config.rectangular_preallocate_workspace,
    rectangular_workspace_min_buses = config.rectangular_workspace_min_buses,
    performance_profile = performance_profile,
  )
end

runpf!(net::Net, config::PowerFlowConfig; kwargs...) = _runpf_with_config!(net, config; kwargs...)
runpf!(net::Net, config::SparlectraConfig; kwargs...) = _runpf_with_config!(net, config.powerflow; kwargs...)

function runpf!(net::Net; config::Union{Nothing,PowerFlowConfig,SparlectraConfig} = nothing, kwargs...)
  runtime_keys = Set((:verbose, :damp, :pv_table_rows, :validate_limits_after_pf, :q_limit_violation_headroom, :qlimit_lock_reason, :performance_profile))
  cfg0 = config === nothing ? powerflow_config() : (config isa SparlectraConfig ? config.powerflow : config)
  if !isempty(kwargs)
    raw = Dict{String,Any}(String(k) => v for (k, v) in pairs(kwargs))
    if all(k -> k in runtime_keys, keys(kwargs))
      return _runpf_with_config!(
        net,
        cfg0;
        verbose = Int(get(raw, "verbose", 0)),
        damp = get(raw, "damp", 1.0),
        pv_table_rows = Int(get(raw, "pv_table_rows", 30)),
        validate_limits_after_pf = Bool(get(raw, "validate_limits_after_pf", false)),
        q_limit_violation_headroom = Float64(get(raw, "q_limit_violation_headroom", 0.0)),
        qlimit_lock_reason = Symbol(get(raw, "qlimit_lock_reason", :manual)),
        performance_profile = get(raw, "performance_profile", nothing),
      )
    else
      throw(ArgumentError("runpf!: solver options must be supplied through PowerFlowConfig or set_sparlectra_config!; only runtime keywords $(collect(runtime_keys)) are accepted."))
    end
  end
  return _runpf_with_config!(net, cfg0)
end

function _active_link_representative_map(net::Net)
  n = length(net.nodeVec)
  parent = collect(1:n)

  function find_root(i::Int)
    while parent[i] != i
      parent[i] = parent[parent[i]]
      i = parent[i]
    end
    return i
  end

  function union_set(a::Int, b::Int)
    ra = find_root(a)
    rb = find_root(b)
    ra == rb && return
    if ra < rb
      parent[rb] = ra
    else
      parent[ra] = rb
    end
  end

  for l in net.linkVec
    l.status == 1 || continue
    union_set(Int(l.fromBus), Int(l.toBus))
  end

  return [find_root(i) for i in 1:n]
end

function _merged_pf_net(net::Net)
  reps = _active_link_representative_map(net)
  all(reps[i] == i for i in eachindex(reps)) && return net, reps, false

  wnet = deepcopy(net)
  n = length(wnet.nodeVec)
  cluster_members = [Int[] for _ in 1:n]
  for bus in 1:n
    push!(cluster_members[reps[bus]], bus)
  end

  for rep in 1:n
    members = cluster_members[rep]
    isempty(members) && continue

    ref = wnet.nodeVec[rep]
    p_load = 0.0
    q_load = 0.0
    p_gen = 0.0
    q_gen = 0.0
    p_sh = 0.0
    q_sh = 0.0
    has_slack = false
    has_pv = false

    for b in members
      nref = wnet.nodeVec[b]
      p_load += isnothing(nref._pƩLoad) ? 0.0 : nref._pƩLoad
      q_load += isnothing(nref._qƩLoad) ? 0.0 : nref._qƩLoad
      p_gen += isnothing(nref._pƩGen) ? 0.0 : nref._pƩGen
      q_gen += isnothing(nref._qƩGen) ? 0.0 : nref._qƩGen
      p_sh += isnothing(nref._pShunt) ? 0.0 : nref._pShunt
      q_sh += isnothing(nref._qShunt) ? 0.0 : nref._qShunt
      has_slack |= (nref._nodeType == Slack)
      has_pv |= (nref._nodeType == PV)
    end

    ref._pƩLoad = p_load
    ref._qƩLoad = q_load
    ref._pƩGen = p_gen
    ref._qƩGen = q_gen
    ref._pShunt = p_sh
    ref._qShunt = q_sh
    ref._nodeType = has_slack ? Slack : (has_pv ? PV : PQ)

    for b in members
      b == rep && continue
      nref = wnet.nodeVec[b]
      nref._pƩLoad = 0.0
      nref._qƩLoad = 0.0
      nref._pƩGen = 0.0
      nref._qƩGen = 0.0
      nref._pShunt = 0.0
      nref._qShunt = 0.0
      nref._nodeType = Isolated
    end
  end

  for br in wnet.branchVec
    f = reps[Int(br.fromBus)]
    t = reps[Int(br.toBus)]
    br.fromBus = f
    br.toBus = t
    if f == t
      br.status = 0
    end
  end

  for sh in wnet.shuntVec
    sh.busIdx = reps[Int(sh.busIdx)]
    if hasproperty(sh.comp, :cFrom_bus) && getfield(sh.comp, :cFrom_bus) !== nothing
      setfield!(sh.comp, :cFrom_bus, sh.busIdx)
    end
  end

  for ps in wnet.prosumpsVec
    if hasproperty(ps.comp, :cFrom_bus) && getfield(ps.comp, :cFrom_bus) !== nothing
      new_bus = reps[Int(getfield(ps.comp, :cFrom_bus))]
      setfield!(ps.comp, :cFrom_bus, new_bus)
    end
    if hasproperty(ps.comp, :cTo_bus) && getfield(ps.comp, :cTo_bus) !== nothing
      new_bus = reps[Int(getfield(ps.comp, :cTo_bus))]
      setfield!(ps.comp, :cTo_bus, new_bus)
    end
  end

  empty!(wnet.isoNodes)
  markIsolatedBuses!(net = wnet, log = false)
  return wnet, reps, true
end

function runpf!(
  net::Net,
  maxIte::Int,
  tolerance::Float64 = 1e-6,
  verbose::Int = 0;
  method::Symbol = :rectangular,
  damp = 1.0,
  autodamp::Bool = false,
  autodamp_min::Float64 = 0.05,
  opt_flatstart::Bool = net.flatstart,
  pv_table_rows::Int = 30,
  validate_limits_after_pf::Bool = false,
  q_limit_violation_headroom::Float64 = 0.0,
  lock_pv_to_pq_buses::AbstractVector{Int} = Int[],
  qlimit_mode::Symbol = :switch_to_pq,
  qlimit_max_outer::Int = 30,
  start_projection::Bool = false,
  start_projection_try_dc_start::Bool = true,
  start_projection_try_blend_scan::Bool = true,
  start_projection_branch_guard::Bool = true,
  start_projection_measure_candidates::Bool = true,
  start_projection_accept_unmeasured_dc_start::Bool = false,
  start_projection_blend_lambdas::AbstractVector{<:Real} = [0.25, 0.5, 0.75],
  start_projection_dc_angle_limit_deg::Float64 = 60.0,
  qlimit_start_iter::Int = 2,
  qlimit_start_mode::Symbol = :iteration,
  qlimit_auto_q_delta_pu::Float64 = 1e-4,
  qlimit_trace_buses::AbstractVector{Int} = Int[],
  qlimit_lock_reason::Symbol = :manual,
  qlimit_guard::Bool = false,
  qlimit_guard_min_q_range_pu::Float64 = 1e-4,
  qlimit_guard_zero_range_mode::Symbol = :lock_pq,
  qlimit_guard_narrow_range_mode::Symbol = :prefer_pq,
  qlimit_guard_log::Bool = true,
  qlimit_guard_max_switches::Int = 10,
  qlimit_guard_accept_bounded_violations::Bool = false,
  qlimit_guard_max_remaining_violations::Int = 0,
  qlimit_guard_freeze_after_repeated_switching::Bool = true,
  qlimit_guard_violation_mode::Symbol = :delayed_switch,
  qlimit_guard_violation_threshold_pu::Float64 = 1e-4,
  wrong_branch_detection::Symbol = :warn,
  wrong_branch_rescue::Bool = false,
  wrong_branch_min_vm_pu::Float64 = 0.70,
  wrong_branch_max_vm_pu::Float64 = 1.30,
  wrong_branch_max_angle_spread_deg::Float64 = 180.0,
  wrong_branch_max_branch_angle_deg::Float64 = 90.0,
  wrong_branch_min_low_vm_count::Int = 1,
  wrong_branch_rescue_max_attempts::Int = 2,
  rectangular_workspace_reuse::Bool = true,
  rectangular_preallocate_workspace::Symbol = :auto,
  rectangular_workspace_min_buses::Int = 1000,
  performance_profile = nothing,
)
  _validate_rectangular_powerflow_options(method = :rectangular, sparse = true)
  wnet, reps, has_merges = _merged_pf_net(net)
  refreshBusTypesFromProsumers!(wnet)
  has_vdep_control = has_voltage_dependent_control(wnet)

  function _sync_merged_results_to_original!()
    for i in eachindex(net.nodeVec)
      src = wnet.nodeVec[reps[i]]
      net.nodeVec[i]._vm_pu = src._vm_pu
      net.nodeVec[i]._va_deg = src._va_deg
    end
    updateShuntPowers!(net = net)
  end

  if method === :rectangular
    if has_merges
      has_vdep_control && error("runpf!: voltage-dependent injections, including P(U)/Q(U) controllers and bus_shunt_model=voltage_dependent_injection, are not supported with active-link merge handling in rectangular mode. Disable merges or use a topology without internal isolated buses.")
      iters, erg = runpf_rectangular!(
        wnet,
        maxIte,
        tolerance,
        verbose;
        damp = damp,
        autodamp = autodamp,
        autodamp_min = autodamp_min,
        wrong_branch_detection = wrong_branch_detection,
        wrong_branch_rescue = wrong_branch_rescue,
        wrong_branch_min_vm_pu = wrong_branch_min_vm_pu,
        wrong_branch_max_vm_pu = wrong_branch_max_vm_pu,
        wrong_branch_max_angle_spread_deg = wrong_branch_max_angle_spread_deg,
        wrong_branch_max_branch_angle_deg = wrong_branch_max_branch_angle_deg,
        wrong_branch_min_low_vm_count = wrong_branch_min_low_vm_count,
        wrong_branch_rescue_max_attempts = wrong_branch_rescue_max_attempts,
        opt_flatstart = opt_flatstart,
        pv_table_rows = pv_table_rows,
        lock_pv_to_pq_buses = lock_pv_to_pq_buses,
        qlimit_mode = qlimit_mode,
        qlimit_max_outer = qlimit_max_outer,
        start_projection = start_projection,
        start_projection_try_dc_start = start_projection_try_dc_start,
        start_projection_try_blend_scan = start_projection_try_blend_scan,
        start_projection_branch_guard = start_projection_branch_guard,
        start_projection_measure_candidates = start_projection_measure_candidates,
        start_projection_accept_unmeasured_dc_start = start_projection_accept_unmeasured_dc_start,
        start_projection_blend_lambdas = start_projection_blend_lambdas,
        start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
        qlimit_start_iter = qlimit_start_iter,
        qlimit_start_mode = qlimit_start_mode,
        qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu,
        qlimit_trace_buses = qlimit_trace_buses,
        qlimit_lock_reason = qlimit_lock_reason,
        qlimit_guard = qlimit_guard,
        qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu,
        qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode,
        qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode,
        qlimit_guard_log = qlimit_guard_log,
        qlimit_guard_max_switches = qlimit_guard_max_switches,
        qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations,
        qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations,
        qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching,
        qlimit_guard_violation_mode = qlimit_guard_violation_mode,
        qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu,
        rectangular_workspace_reuse = rectangular_workspace_reuse,
        rectangular_preallocate_workspace = rectangular_preallocate_workspace,
        rectangular_workspace_min_buses = rectangular_workspace_min_buses,
        performance_profile = performance_profile,
      )
    else
      iters, erg = runpf_rectangular!(
        wnet,
        maxIte,
        tolerance,
        verbose;
        damp = damp,
        autodamp = autodamp,
        autodamp_min = autodamp_min,
        wrong_branch_detection = wrong_branch_detection,
        wrong_branch_rescue = wrong_branch_rescue,
        wrong_branch_min_vm_pu = wrong_branch_min_vm_pu,
        wrong_branch_max_vm_pu = wrong_branch_max_vm_pu,
        wrong_branch_max_angle_spread_deg = wrong_branch_max_angle_spread_deg,
        wrong_branch_max_branch_angle_deg = wrong_branch_max_branch_angle_deg,
        wrong_branch_min_low_vm_count = wrong_branch_min_low_vm_count,
        wrong_branch_rescue_max_attempts = wrong_branch_rescue_max_attempts,
        opt_flatstart = opt_flatstart,
        pv_table_rows = pv_table_rows,
        lock_pv_to_pq_buses = lock_pv_to_pq_buses,
        qlimit_mode = qlimit_mode,
        qlimit_max_outer = qlimit_max_outer,
        start_projection = start_projection,
        start_projection_try_dc_start = start_projection_try_dc_start,
        start_projection_try_blend_scan = start_projection_try_blend_scan,
        start_projection_branch_guard = start_projection_branch_guard,
        start_projection_measure_candidates = start_projection_measure_candidates,
        start_projection_accept_unmeasured_dc_start = start_projection_accept_unmeasured_dc_start,
        start_projection_blend_lambdas = start_projection_blend_lambdas,
        start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
        qlimit_start_iter = qlimit_start_iter,
        qlimit_start_mode = qlimit_start_mode,
        qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu,
        qlimit_trace_buses = qlimit_trace_buses,
        qlimit_lock_reason = qlimit_lock_reason,
        qlimit_guard = qlimit_guard,
        qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu,
        qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode,
        qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode,
        qlimit_guard_log = qlimit_guard_log,
        qlimit_guard_max_switches = qlimit_guard_max_switches,
        qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations,
        qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations,
        qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching,
        qlimit_guard_violation_mode = qlimit_guard_violation_mode,
        qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu,
        rectangular_workspace_reuse = rectangular_workspace_reuse,
        rectangular_preallocate_workspace = rectangular_preallocate_workspace,
        rectangular_workspace_min_buses = rectangular_workspace_min_buses,
        performance_profile = performance_profile,
      )
    end
    rect_status = rectangular_pf_status(wnet)
    if rect_status !== nothing
      _set_rectangular_pf_status!(net, rect_status)
    end
    if erg == 0 && has_merges
      _sync_merged_results_to_original!()
    end
    if validate_limits_after_pf && (verbose > 0)
      printFinalLimitValidation(has_merges ? net : wnet; q_headroom = q_limit_violation_headroom)
    end
    return iters, erg
  else
    throw(ArgumentError(unsupported_powerflow_method_message(method)))
  end
end

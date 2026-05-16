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
# - Both analytic and finite-difference Jacobian options
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
# - run_complex_nr_rectangular_for_net!(): Main solver interface
# - runpf_rectangular!(): Convenience wrapper matching runpf!() signature
# - build_complex_jacobian(): Wirtinger-based Jacobian block construction
# - mismatch_rectangular(): Residual function for PQ/PV bus constraints
#
# Note:
# - The FD Jacobian is mathematically dense; sparse storage does not bring much benefit.
# - The analytic Jacobian currently uses a dense rectangular build, even though the
#   underlying structure is sparse (Ybus-like). A true sparse implementation would
#   require a dedicated builder similar to `calcJacobian(...; sparse=true)`.
#
# References:
# - Wirtinger calculus for complex derivatives

using LinearAlgebra
using SparseArrays
using Printf

_is_rectangular_linear_step_failure(e) = e isa LinearAlgebra.SingularException || e isa LinearAlgebra.LAPACKException

"""
    build_complex_jacobian(Ybus, V)

Builds the 2n × 2n Wirtinger-type Jacobian blocks for the complex-state
Newton–Raphson formulation.

Given:
    I = Ybus * V
    S = V .* conj.(I)

We construct the blocks:
    J11 = ∂S/∂V
    J12 = ∂S/∂V*
    J21 = ∂conj(S)/∂V
    J22 = ∂conj(S)/∂V*

Returns:
    J11, J12, J21, J22  (all full matrices, not Diagonal)
"""
function build_complex_jacobian(Ybus, V)
  I = Ybus * V
  n = length(V)

  # J11 = diag(conj(I))
  J11 = Matrix(Diagonal(conj.(I)))

  # J12 = diag(V) * conj(Ybus)
  J12 = Matrix(Diagonal(V)) * conj.(Ybus)

  # J21 = conj(J12)
  J21 = conj.(J12)

  # J22 = conj(J11)
  J22 = conj.(J11)

  return J11, J12, J21, J22
end

"""
    mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx) -> F::Vector{Float64}

Compute the real-valued mismatch vector F(V) for the rectangular
complex-state formulation with PQ and PV buses.

For each non-slack bus i:
- if bus_types[i] == :PQ:
      ΔP_i = Re(S_calc[i]) - Re(S_spec[i])
      ΔQ_i = Im(S_calc[i]) - Im(S_spec[i])

- if bus_types[i] == :PV:
      ΔP_i = Re(S_calc[i]) - Re(S_spec[i])
      ΔV_i = |V[i]| - Vset[i]

F is stacked as [ΔP_2, ΔQ/ΔV_2, ..., ΔP_n, ΔQ/ΔV_n] over all non-slack buses.
"""
function mismatch_rectangular(Ybus, V::Vector{ComplexF64}, S::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int)
  # Residual/update convention for the rectangular solver:
  # F = calc - spec, including the PV row ΔV = |V| - Vset, and Newton solves
  # J * δx = -F before applying rectangular absolute updates ΔVr, ΔVi.
  # The positive PV voltage Jacobian row below is therefore consistent with
  # the opposite residual sign used by the polar spec-minus-calc solvers.
  n = length(V)
  @assert length(S) == n
  @assert length(bus_types) == n
  @assert length(Vset) == n
  # Network-based injections for the current state
  S_calc = calc_injections(Ybus, V)

  # F has 2*(n-1) entries: for each non-slack bus two residuals
  # PQ:  ΔP_i, ΔQ_i
  # PV:  ΔP_i, ΔV_i
  F = zeros(Float64, 2 * (n - 1))

  row = 1
  @inbounds for i = 1:n
    if i == slack_idx
      continue
    end

    S_ci = S_calc[i]
    S_si = S[i]

    if bus_types[i] == :PQ
      ΔP = real(S_ci) - real(S_si)
      ΔQ = imag(S_ci) - imag(S_si)

      F[row]   = ΔP
      F[row+1] = ΔQ

    elseif bus_types[i] == :PV
      ΔP = real(S_ci) - real(S_si)
      ΔV = abs(V[i]) - Vset[i]

      F[row]   = ΔP
      F[row+1] = ΔV

    else
      error("mismatch_rectangular: unsupported bus type $(bus_types[i]) at bus $i")
    end

    row += 2
  end

  return F
end

function _bus_voltage_setpoint_from_prosumers(net::Net, bus::Int, fallback::Float64)::Float64
  setpoint = fallback
  found = false
  use_slack_prosumer_setpoint = net.flatstart
  for ps in net.prosumpsVec
    getPosumerBusIndex(ps) == bus || continue
    ((use_slack_prosumer_setpoint && isSlack(ps)) || (isGenerator(ps) && isRegulating(ps))) || continue
    isnothing(ps.vm_pu) && continue
    if found && abs(setpoint - Float64(ps.vm_pu)) > 1e-8
      @debug "Multiple voltage-regulating prosumers with different setpoints at bus $(bus); keeping first setpoint." kept = setpoint ignored = ps.vm_pu
      continue
    end
    setpoint = Float64(ps.vm_pu)
    found = true
  end
  return setpoint
end

function _bus_voltage_setpoints_from_prosumers(net::Net)::Vector{Float64}
  Vset = Vector{Float64}(undef, length(net.nodeVec))
  @inbounds for k in eachindex(net.nodeVec)
    node = net.nodeVec[k]
    fallback = isnothing(node._vm_pu) ? 1.0 : Float64(node._vm_pu)
    Vset[k] = _bus_voltage_setpoint_from_prosumers(net, k, fallback)
  end
  return Vset
end

function _max_rectangular_pv_voltage_residual(V::Vector{ComplexF64}, Vset::Vector{Float64}, bus_types::Vector{Symbol}, q_limit_events::Dict{Int,Symbol})::Float64
  max_residual = 0.0
  @inbounds for k in eachindex(V)
    bus_types[k] == :PV || continue
    haskey(q_limit_events, k) && continue
    max_residual = max(max_residual, abs(abs(V[k]) - Vset[k]))
  end
  return max_residual
end

_qlimit_original_bus_id(net::Net, bus::Int)::Int = haskey(net.busOrigIdxDict, bus) ? net.busOrigIdxDict[bus] : bus

function _resolve_qlimit_trace_buses(net::Net, requested::AbstractVector{Int})::Vector{Int}
  isempty(requested) && return Int[]
  orig_to_net = Dict{Int,Int}()
  sizehint!(orig_to_net, length(net.busOrigIdxDict))
  for (net_idx, orig_idx) in net.busOrigIdxDict
    orig_to_net[orig_idx] = net_idx
  end
  resolved = Int[]
  for bus in requested
    if haskey(orig_to_net, bus)
      push!(resolved, orig_to_net[bus])
    elseif 1 <= bus <= length(net.nodeVec)
      push!(resolved, bus)
    end
  end
  unique!(resolved)
  sort!(resolved)
  return resolved
end

function _bus_has_online_voltage_regulator(net::Net, bus::Int)::Bool
  for ps in net.prosumpsVec
    getPosumerBusIndex(ps) == bus || continue
    (isSlack(ps) || (isGenerator(ps) && isRegulating(ps))) && return true
  end
  return false
end

function _qlimit_violation(qreq::Float64, qmin::Float64, qmax::Float64; hyst::Float64 = 0.0)
  if isfinite(qmax) && qreq > qmax + hyst
    return (:high, qreq - qmax)
  elseif isfinite(qmin) && qreq < qmin - hyst
    return (:low, qmin - qreq)
  end
  return (:none, 0.0)
end

function _print_rectangular_qlimit_trace(
  io::IO,
  net::Net,
  it::Int,
  bus::Int,
  bus_type::Symbol,
  V::Vector{ComplexF64},
  Vset::Vector{Float64},
  qreq_pu::Float64,
  qmin_pu::AbstractVector,
  qmax_pu::AbstractVector;
  q_hyst_pu::Float64,
  qlimit_start_iter::Int,
  qlimit_start_mode::Symbol,
  qlimit_iter_ready::Bool,
  qlimit_auto_ready::Bool,
  qlimit_ready::Bool,
  qlimit_check_active::Bool,
  converged_before_switching::Bool,
  qlimit_mode::Symbol,
  lock_mask::AbstractVector{Bool},
  qlimit_lock_reason::Symbol = :manual,
)
  busI = _qlimit_original_bus_id(net, bus)
  qmin = (bus <= length(qmin_pu)) ? qmin_pu[bus] : -Inf
  qmax = (bus <= length(qmax_pu)) ? qmax_pu[bus] : Inf
  side, amount_pu = _qlimit_violation(qreq_pu, qmin, qmax; hyst = 0.0)
  side_hyst, _ = _qlimit_violation(qreq_pu, qmin, qmax; hyst = q_hyst_pu)
  has_limits = has_q_limits(qmin_pu, qmax_pu, bus)
  locked = bus <= length(lock_mask) && lock_mask[bus]
  has_regulator = _bus_has_online_voltage_regulator(net, bus)

  decision = "keep PV"
  reason = "active PV/REF bus checked; no physical Q-limit violation"
  if bus_type == :Slack
    decision = "not applicable"
    reason = "slack/reference bus is excluded from PV->PQ switching"
  elseif bus_type != :PV
    decision = bus_type == :PQ ? "keep PQ" : "not applicable"
    reason = bus_type == :PQ ? "bus is no longer an active PV/REF bus" : "bus is not an active PV/REF bus in this iteration"
  elseif locked
    decision = "not applicable"
    reason = qlimit_lock_reason == :ignore_q_limits ? "ignore_q_limits=true disables Q-limit switching for this bus" : "bus is locked/excluded from PV->PQ switching"
  elseif !has_regulator
    decision = "not applicable"
    reason = "missing voltage setpoint or online voltage-regulating generator"
  elseif !has_limits
    decision = "not applicable"
    reason = "bus has no finite Q limit"
  elseif side == :none
    reason = "active PV/REF bus checked; no physical Q-limit violation"
  elseif converged_before_switching
    reason = "active PV/REF bus checked, but solver convergence was reached before rectangular active-set switching in this iteration"
  elseif !qlimit_iter_ready && qlimit_start_mode in (:iteration, :iteration_or_auto)
    reason = "active PV/REF bus checked, but qlimit_start_iter has not been reached"
  elseif !qlimit_ready
    reason = "active PV/REF bus checked, but qlimit_start_mode suppresses switching in this iteration"
  elseif side_hyst == :none
    reason = "active PV/REF bus checked, but physical violation is inside q_hyst_pu switching deadband"
  elseif qlimit_mode == :switch_to_pq
    decision = "switch to PQ"
    reason = "active PV/REF bus checked; Q-limit enforcement is active and violation exceeds hysteresis"
  elseif qlimit_mode == :adjust_vset
    decision = "keep PV"
    reason = "active PV/REF bus checked; qlimit_mode=:adjust_vset handles the violation by voltage-setpoint adjustment instead of PV->PQ switching"
  else
    decision = "not applicable"
    reason = "rectangular solver path does not enforce this Q-limit mode"
  end

  @printf(io, "Q-limit trace it=%d BUS_I=%d type=%s Vm_calc=%.6f imported_Vset=%.6f Qcalc=%.6f pu (%.3f MVAr) Qmin=%.6f pu (%.3f MVAr) Qmax=%.6f pu (%.3f MVAr) q_hyst_pu=%.6f qlimit_start_iter=%d qlimit_start_mode=%s active=%s violation=%s amount=%.6f pu (%.3f MVAr) decision=%s reason=%s\n",
          it, busI, String(bus_type), abs(V[bus]), Vset[bus], qreq_pu, qreq_pu * net.baseMVA,
          qmin, qmin * net.baseMVA, qmax, qmax * net.baseMVA, q_hyst_pu, qlimit_start_iter,
          String(qlimit_start_mode), string(qlimit_check_active), String(side), amount_pu, amount_pu * net.baseMVA,
          decision, reason)
  return nothing
end

mutable struct _RectangularPFStatusTable
  entries::Vector{Tuple{UInt,WeakRef,Any}}
end

const _RECTANGULAR_PF_STATUS = _RectangularPFStatusTable(Tuple{UInt,WeakRef,Any}[])

function _prune_rectangular_pf_status!(table::_RectangularPFStatusTable)
  filter!(entry -> entry[2].value !== nothing, table.entries)
  return table
end

function _set_rectangular_pf_status!(net::Net, status)
  _prune_rectangular_pf_status!(_RECTANGULAR_PF_STATUS)
  key = objectid(net)
  for i in eachindex(_RECTANGULAR_PF_STATUS.entries)
    entry = _RECTANGULAR_PF_STATUS.entries[i]
    if entry[1] == key && entry[2].value === net
      _RECTANGULAR_PF_STATUS.entries[i] = (key, WeakRef(net), status)
      return status
    end
  end
  push!(_RECTANGULAR_PF_STATUS.entries, (key, WeakRef(net), status))
  return status
end

function rectangular_pf_status(net::Net)
  _prune_rectangular_pf_status!(_RECTANGULAR_PF_STATUS)
  key = objectid(net)
  for entry in _RECTANGULAR_PF_STATUS.entries
    entry[1] == key && entry[2].value === net && return entry[3]
  end
  return nothing
end

function _rectangular_rejection_reason_text(reason::Symbol)
  reason == :none && return "none"
  reason == :remaining_pv_q_limit_violations && return "remaining PV Q-limit violations"
  reason == :active_pv_voltage_residual && return "active PV voltage residual exceeds tolerance"
  reason == :singular_newton_step && return "singular Newton step"
  reason == :nr_mismatch_not_converged && return "NR mismatch did not converge"
  return replace(String(reason), "_" => " ")
end

function _print_rectangular_convergence_summary(io::IO, status)
  numerical_text = status.numerical_converged ? "OK" : "FAIL"
  active_text = status.q_limit_active_set_ok ? "OK" : "FAIL"
  final_text = status.final_converged ? "true" : "false"
  @printf(io, "rectangular convergence: numerical_solution=%s  q_limit_active_set=%s  final_converged=%s  reason=%s\n", numerical_text, active_text, final_text, status.reason_text)
  return nothing
end

function _rectangular_solver_status_symbol(numerical_converged::Bool, active_set_ok::Bool, final_converged::Bool, reason::Symbol)::Symbol
  final_converged && return :converged
  numerical_converged || return reason == :singular_newton_step ? :singular_jacobian : :nr_mismatch_not_converged
  active_set_ok && return :converged_with_active_set_instability
  reason == :max_switching_exceeded && return :max_switching_exceeded
  reason == :bounded_q_limit_violations_accepted && return :converged_with_bounded_q_limit_violations
  reason == :remaining_pv_q_limit_violations && return :numerically_converged_but_active_set_failed
  return :qlimit_chatter
end

function _print_qlimit_active_set_summary(io::IO, status)
  println(io, "==================== Q-Limit Active-Set Summary ====================")
  println(io)
  @printf(io, "NR convergence             : %s\n", status.numerical_converged ? "yes" : "no")
  @printf(io, "Final mismatch             : %.6g\n", status.final_mismatch)
  @printf(io, "Active-set convergence     : %s\n", status.q_limit_active_set_ok ? "yes" : "no")
  @printf(io, "PV→PQ switching events     : %d\n", status.pv_pq_switching_events)
  @printf(io, "Oscillating buses          : %d\n", status.oscillating_buses)
  @printf(io, "Guarded narrow-Q PV buses  : %d\n", status.guarded_narrow_q_pv_buses)
  @printf(io, "Final status               : %s\n", String(status.status))
  println(io)
  println(io, "===================================================================")
  return nothing
end

function _apply_qlimit_guard_to_rectangular_active_set!(net::Net, bus_types::Vector{Symbol}, S::Vector{ComplexF64}, Qload_pu::Vector{Float64}, qmin_pu::AbstractVector, qmax_pu::AbstractVector; min_q_range_pu::Float64, zero_range_mode::Symbol, narrow_range_mode::Symbol, log::Bool, verbose::Int)
  min_q_range_pu >= 0.0 || error("qlimit_guard_min_q_range_pu must be >= 0 (got $(min_q_range_pu)).")
  zero_range_mode in (:lock_pq, :prefer_pq, :delayed_switch, :ignore) || error("Unsupported qlimit_guard_zero_range_mode=$(zero_range_mode). Supported: :lock_pq, :prefer_pq, :delayed_switch, :ignore.")
  narrow_range_mode in (:lock_pq, :prefer_pq, :delayed_switch, :ignore) || error("Unsupported qlimit_guard_narrow_range_mode=$(narrow_range_mode). Supported: :lock_pq, :prefer_pq, :delayed_switch, :ignore.")

  guarded = Int[]
  @inbounds for bus in eachindex(bus_types)
    bus_types[bus] == :PV || continue
    bus <= length(qmin_pu) && bus <= length(qmax_pu) || continue
    qmin = qmin_pu[bus]
    qmax = qmax_pu[bus]
    isfinite(qmin) && isfinite(qmax) || continue
    qrange = abs(qmax - qmin)
    qrange < min_q_range_pu || continue
    mode = qrange <= eps(Float64) ? zero_range_mode : narrow_range_mode
    mode in (:lock_pq, :prefer_pq) || continue

    qclamp = 0.5 * (qmin + qmax)
    bus_types[bus] = :PQ
    S[bus] = ComplexF64(real(S[bus]), qclamp - Qload_pu[bus])
    net.nodeVec[bus]._qƩGen = qclamp * net.baseMVA
    logQLimitHit!(net, 0, bus, qclamp >= 0.0 ? :max : :min)
    push!(guarded, bus)
  end

  if log && verbose > 0 && !isempty(guarded)
    @printf("Q-limit guard: locked %d narrow-range PV bus(es) as PQ before rectangular NR.\n", length(guarded))
  end
  return guarded
end

function _print_rectangular_qlimit_summary(io::IO, net::Net, V::Vector{ComplexF64}, Sbus_pu::Vector{ComplexF64}, bus_types::Vector{Symbol}, qmin_pu::AbstractVector, qmax_pu::AbstractVector, Qload_pu::Vector{Float64}; q_hyst_pu::Float64, tolerance_pu::Float64 = 0.0, max_rows::Int = 30, max_console_rows::Union{Nothing,Int} = nothing)
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
    push!(rows, (
      bus = bus,
      busI = _qlimit_original_bus_id(net, bus),
      type = bus_type,
      qcalc = qcalc,
      qmin = qmin,
      qmax = qmax,
      side = side,
      amount = amount,
      switched = !isnothing(last_it),
      last_it = last_it,
      in_return_band = in_return_band,
    ))
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
      @printf(io, " %6d │ %4s │ %8.5f │ %10.3f │ %7.5f │ %9.3f │ %7.5f │ %9.3f │ %4s │ %7.5f │ %9.3f │ %8s │ %7s │ %s\n",
              row.busI, type_text, row.qcalc, row.qcalc * net.baseMVA,
              row.qmin, row.qmin * net.baseMVA, row.qmax, row.qmax * net.baseMVA,
              String(row.side), row.amount, row.amount * net.baseMVA,
              string(row.switched), last_text, band_text)
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

function run_complex_nr_rectangular(
  Ybus,
  V0,
  S;
  slack_idx::Int = 1,
  maxiter::Int = 20,
  tol::Float64 = 1e-8,
  verbose::Bool = false,
  damp::Float64 = 1.0,
  autodamp::Bool = false,
  autodamp_min::Float64 = 1e-3,
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  use_fd::Bool = false,
  use_sparse::Bool = false,
  dPinj_dVm::Vector{Float64} = zeros(Float64, length(V0)),
  dQinj_dVm::Vector{Float64} = zeros(Float64, length(V0)),
)
  V = copy(V0)
  history = Float64[]

  for iter = 1:maxiter
    F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
    max_mis = maximum(abs.(F))
    push!(history, max_mis)

    if max_mis <= tol
      return V, true, iter, history
    end

    if use_fd
      V = complex_newton_step_rectangular_fd(Ybus, V, S; slack_idx = slack_idx, damp = damp, autodamp = autodamp, autodamp_min = autodamp_min, h = 1e-6, bus_types = bus_types, Vset = Vset, dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm)
    else
      V = complex_newton_step_rectangular(Ybus, V, S; slack_idx = slack_idx, damp = damp, autodamp = autodamp, autodamp_min = autodamp_min, bus_types = bus_types, Vset = Vset, use_sparse = use_sparse, dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm)
    end
  end

  return V, false, maxiter, history
end

function _sanitize_rectangular_start(V::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int)
  Vs = copy(V)
  @inbounds for k in eachindex(Vs)
    Vk = Vs[k]
    vm = abs(Vk)
    if !isfinite(real(Vk)) || !isfinite(imag(Vk)) || vm <= 0.0
      vm = (bus_types[k] in (:Slack, :PV) && isfinite(Vset[k]) && Vset[k] > 0.0) ? Vset[k] : 1.0
      Vs[k] = ComplexF64(vm, 0.0)
    end
  end
  Vs[slack_idx] = V[slack_idx]
  return Vs
end

function _voltage_magnitude_for_projection(Vraw::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, k::Int)
  vm_raw = abs(Vraw[k])
  if bus_types[k] in (:Slack, :PV) && isfinite(Vset[k]) && Vset[k] > 0.0
    return Vset[k]
  elseif isfinite(vm_raw) && vm_raw > 0.0
    return vm_raw
  else
    return 1.0
  end
end

function _dc_angle_start_rectangular(
  Ybus,
  Vraw::Vector{ComplexF64},
  S::Vector{ComplexF64},
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  slack_idx::Int;
  dc_angle_limit_deg::Float64 = 60.0,
)
  n = length(Vraw)
  non_slack = non_slack_indices(n, slack_idx)
  nred = length(non_slack)
  nred == 0 && return copy(Vraw)
  pos = build_pos_map(non_slack, n)

  B = zeros(Float64, nred, nred)
  if Ybus isa SparseMatrixCSC
    rv = rowvals(Ybus)
    nz = nonzeros(Ybus)
    @inbounds for j in 1:n
      for ptr in nzrange(Ybus, j)
        i = rv[ptr]
        i == j && continue
        ri = pos[i]
        ri == 0 && continue
        bij = imag(nz[ptr])
        bij == 0.0 && continue
        if j != slack_idx
          cj = pos[j]
          if cj != 0
            B[ri, cj] -= bij
          end
        end
        B[ri, ri] += bij
      end
    end
  else
    @inbounds for i in non_slack
      ri = pos[i]
      for j in 1:n
        j == i && continue
        bij = imag(Ybus[i, j])
        bij == 0.0 && continue
        if j != slack_idx
          cj = pos[j]
          if cj != 0
            B[ri, cj] -= bij
          end
        end
        B[ri, ri] += bij
      end
    end
  end

  P = real.(S[non_slack])
  θred = solve_linear(B, P; allow_pinv = true)
  limit = deg2rad(dc_angle_limit_deg)
  θslack = angle(Vraw[slack_idx])
  Vdc = similar(Vraw)

  @inbounds for k in 1:n
    vm = _voltage_magnitude_for_projection(Vraw, bus_types, Vset, k)
    if k == slack_idx
      Vdc[k] = Vraw[k]
    else
      θ = θslack + clamp(θred[pos[k]], -limit, limit)
      Vdc[k] = ComplexF64(vm * cos(θ), vm * sin(θ))
    end
  end
  return Vdc
end

function _blend_voltage_starts(Vraw::Vector{ComplexF64}, Vdc::Vector{ComplexF64}, λ::Float64, slack_idx::Int)
  0.0 <= λ <= 1.0 || error("blend lambda must satisfy 0 ≤ λ ≤ 1 (got $(λ)).")
  V = similar(Vraw)
  @inbounds for k in eachindex(Vraw)
    if k == slack_idx
      V[k] = Vraw[k]
      continue
    end
    vm = (1.0 - λ) * abs(Vraw[k]) + λ * abs(Vdc[k])
    θ = (1.0 - λ) * angle(Vraw[k]) + λ * angle(Vdc[k])
    V[k] = ComplexF64(vm * cos(θ), vm * sin(θ))
  end
  return V
end

"""
    project_rectangular_start(Ybus, Vraw, S, bus_types, Vset, slack_idx; ...)

Build a projected initial voltage for the rectangular power-flow solver.

When enabled, the projection sanitizes the raw seed, optionally computes a
DC-angle start from active-power injections and the Y-bus off-diagonal
susceptances, and optionally scans convex blends between the raw and DC starts.
The candidate with the lowest rectangular mismatch is returned.
"""
function project_rectangular_start(
  Ybus,
  Vraw::Vector{ComplexF64},
  S::Vector{ComplexF64},
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  slack_idx::Int;
  enabled::Bool = false,
  try_dc_start::Bool = true,
  try_blend_scan::Bool = true,
  blend_lambdas::AbstractVector{<:Real} = [0.25, 0.5, 0.75],
  dc_angle_limit_deg::Float64 = 60.0,
  verbose::Int = 0,
)
  enabled || return Vraw
  dc_angle_limit_deg > 0.0 || error("dc_angle_limit_deg must be > 0 (got $(dc_angle_limit_deg)).")

  raw = _sanitize_rectangular_start(Vraw, bus_types, Vset, slack_idx)
  best = raw
  best_name = :raw
  best_mis = _max_rectangular_mismatch(Ybus, raw, S, bus_types, Vset, slack_idx)
  Vdc = nothing

  if try_dc_start
    Vdc = _dc_angle_start_rectangular(Ybus, raw, S, bus_types, Vset, slack_idx; dc_angle_limit_deg = dc_angle_limit_deg)
    dc_mis = _max_rectangular_mismatch(Ybus, Vdc, S, bus_types, Vset, slack_idx)
    if isfinite(dc_mis) && dc_mis < best_mis
      best = Vdc
      best_name = :dc_start
      best_mis = dc_mis
    end
  end

  if try_blend_scan && Vdc !== nothing
    for λ_raw in blend_lambdas
      λ = Float64(λ_raw)
      Vblend = _blend_voltage_starts(raw, Vdc, λ, slack_idx)
      blend_mis = _max_rectangular_mismatch(Ybus, Vblend, S, bus_types, Vset, slack_idx)
      if isfinite(blend_mis) && blend_mis < best_mis
        best = Vblend
        best_name = Symbol("blend_", λ)
        best_mis = blend_mis
      end
    end
  end

  if verbose > 0
    raw_mis = _max_rectangular_mismatch(Ybus, raw, S, bus_types, Vset, slack_idx)
    @info "start projection selected $(best_name)" raw_mismatch = raw_mis projected_mismatch = best_mis
  end
  return best
end

function _validate_rectangular_damping(damp::Float64, autodamp_min::Float64)
  isfinite(damp) || error("damp must be finite (got $(damp)).")
  isfinite(autodamp_min) || error("autodamp_min must be finite (got $(autodamp_min)).")
  0.0 < damp <= 1.0 || error("damp must satisfy 0 < damp ≤ 1 (got $(damp)).")
  0.0 < autodamp_min <= damp || error("autodamp_min must satisfy 0 < autodamp_min ≤ damp (got autodamp_min=$(autodamp_min), damp=$(damp)).")
  return nothing
end

function _apply_rectangular_delta(V::Vector{ComplexF64}, δx::Vector{Float64}, slack_idx::Int, non_slack::Vector{Int}, alpha::Float64)
  n = length(V)
  Vr = real.(V)
  Vi = imag.(V)
  Vr_new = copy(Vr)
  Vi_new = copy(Vi)

  @inbounds for (idx, bus) in enumerate(non_slack)
    Vr_new[bus] += alpha * δx[idx]
    Vi_new[bus] += alpha * δx[(n - 1) + idx]
  end

  Vr_new[slack_idx] = Vr[slack_idx]
  Vi_new[slack_idx] = Vi[slack_idx]
  return ComplexF64.(Vr_new, Vi_new)
end

function _max_rectangular_mismatch(Ybus, V::Vector{ComplexF64}, S::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int)
  F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
  return maximum(abs.(F))
end

"""
    choose_rectangular_autodamp(Ybus, V, S, δx, F0; slack_idx, damp, autodamp_min, bus_types, Vset)

Select a Newton step length for the rectangular power-flow solver by backtracking
from `damp` toward `autodamp_min`. The first trial step that reduces the maximum
absolute mismatch is accepted. If no trial reduces the mismatch, the smallest
finite trial is returned so the solver can continue safely with a conservative
step.

Returns `(alpha, Vtrial, trial_mismatch)`.
"""
function choose_rectangular_autodamp(
  Ybus,
  V::Vector{ComplexF64},
  S::Vector{ComplexF64},
  δx::Vector{Float64},
  F0::Vector{Float64};
  slack_idx::Int,
  damp::Float64 = 1.0,
  autodamp_min::Float64 = 1e-3,
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
)
  _validate_rectangular_damping(damp, autodamp_min)
  non_slack = non_slack_indices(length(V), slack_idx)
  current_mismatch = maximum(abs.(F0))
  best_alpha = autodamp_min
  best_V = _apply_rectangular_delta(V, δx, slack_idx, non_slack, autodamp_min)
  best_mismatch = _max_rectangular_mismatch(Ybus, best_V, S, bus_types, Vset, slack_idx)

  alpha = damp
  while alpha >= autodamp_min
    Vtrial = _apply_rectangular_delta(V, δx, slack_idx, non_slack, alpha)
    trial_mismatch = _max_rectangular_mismatch(Ybus, Vtrial, S, bus_types, Vset, slack_idx)
    if isfinite(trial_mismatch) && trial_mismatch < current_mismatch
      return alpha, Vtrial, trial_mismatch
    end
    if isfinite(trial_mismatch) && trial_mismatch < best_mismatch
      best_alpha = alpha
      best_V = Vtrial
      best_mismatch = trial_mismatch
    end
    alpha *= 0.5
  end

  return best_alpha, best_V, best_mismatch
end

"""
    update_net_voltages_from_complex!(net, V)

Update the bus voltage magnitudes and angles in the network from the
final complex voltages V (in per-unit).
"""
function update_net_voltages_from_complex!(net::Net, V::Vector{ComplexF64})
  nodes = net.nodeVec
  n = length(nodes)
  @assert length(V) == n

  for (k, node) in enumerate(nodes)
    Vk = V[k]
    vm = abs(Vk)
    va_rad = angle(Vk)
    va_deg = rad2deg(va_rad)
    node._vm_pu = vm
    node._va_deg = va_deg
  end
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

function _build_vset_adjust_controllers(net::Net)
  controllers = Dict{Int,NamedTuple{(:prosumer_idx, :config),Tuple{Int,VoltageAdjustConfig}}}()

  for (ps_idx, ps) in enumerate(net.prosumpsVec)
    isGenerator(ps) || continue
    bus = getPosumerBusIndex(ps)

    has_any = _has_vset_adjust_config(ps)
    if has_any
      cfg = _resolve_vset_adjust_config(ps)
      all_defined = !isnothing(cfg)
      all_defined || error("Bus $(_bus_label(net, bus)): invalid voltage adjustment config at prosumer $ps_idx. vstep_pu, tap_steps_down and tap_steps_up must be provided together.")
      cfg.vstep_pu > 0.0 || error("Bus $(_bus_label(net, bus)): invalid vstep_pu=$(cfg.vstep_pu). Must be > 0.")
      cfg.tap_steps_down >= 0 || error("Bus $(_bus_label(net, bus)): invalid tap_steps_down=$(cfg.tap_steps_down). Must be ≥ 0.")
      cfg.tap_steps_up >= 0 || error("Bus $(_bus_label(net, bus)): invalid tap_steps_up=$(cfg.tap_steps_up). Must be ≥ 0.")
      haskey(controllers, bus) && error("Bus $(_bus_label(net, bus)): multiple prosumers define voltage adjustment data. Only one controller per bus is allowed.")

      controllers[bus] = (prosumer_idx = ps_idx, config = cfg)
    else
      partially_set = !isnothing(ps.vstep_pu) || !isnothing(ps.tap_steps_down) || !isnothing(ps.tap_steps_up)
      partially_set && error("Bus $(_bus_label(net, bus)): incomplete voltage adjustment config at prosumer $ps_idx. vstep_pu, tap_steps_down and tap_steps_up must be provided together.")
    end
  end

  return controllers
end

function _try_adjust_vset_on_q_limit!(
  net::Net,
  bus::Int,
  side::Symbol,
  it::Int,
  controllers::Dict{Int,NamedTuple{(:prosumer_idx, :config),Tuple{Int,VoltageAdjustConfig}}},
  base_vset::Vector{Float64},
  Vset::Vector{Float64},
  adjust_counter::Vector{Int},
  qlimit_max_outer::Int,
  verbose::Int,
)::Bool
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
    build_rectangular_jacobian_pq_pv_sparse(
        Ybus,
        V,
        bus_types,
        Vset,
        slack_idx,
    ) -> SparseMatrixCSC{Float64}

Builds the analytic rectangular Jacobian corresponding to `mismatch_rectangular`
using the sparsity pattern of `Ybus`.

State vector:
    x = [Vr(non-slack); Vi(non-slack)] ∈ ℝ^(2(n-1))

Residual F(V):
    - PQ buses: ΔP_i, ΔQ_i
    - PV buses: ΔP_i, ΔV_i
    - Slack bus: no equations

Jacobian entries are derived from
    S_i(V) = V_i * conj( (Ybus * V)_i )

Wirtinger-based identities:
    ∂S/∂V   = diag(conj(I)) + diag(V) * conj(Ybus)
    ∂S/∂V*  = diag(V) * conj(Ybus)

Chain rule to rectangular:
    ∂S/∂Vr = ∂S/∂V + ∂S/∂V*
    ∂S/∂Vi = j(∂S/∂V - ∂S/∂V*)

With ΔP_i = Re(ΔS_i), ΔQ_i = Im(ΔS_i), ΔV_i = |V_i| - Vset[i].

Returns:
    J :: SparseMatrixCSC{Float64} with size (2(n-1)) × (2(n-1)).
"""
function build_rectangular_jacobian_pq_pv_sparse(
  Ybus::SparseMatrixCSC{ComplexF64},
  V::Vector{ComplexF64},
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  slack_idx::Int;
  dPinj_dVm::Vector{Float64} = zeros(Float64, length(V)),
  dQinj_dVm::Vector{Float64} = zeros(Float64, length(V)),
  vm_eps::Float64 = 1e-9,
)
  n = length(V)
  @assert length(bus_types) == n
  @assert length(Vset) == n
  @assert length(dPinj_dVm) == n
  @assert length(dQinj_dVm) == n

  I = Ybus * V

  non_slack = non_slack_indices(n, slack_idx)
  pos_non_slack = build_pos_map(non_slack, n)
  nvar = 2 * (n - 1)   # [Vr(non-slack); Vi(non-slack)]
  m = nvar          # F has 2 equations per non-slack bus

  # Row blocks: for each non-slack bus i
  #   row_block[i]   = index of ΔP row for bus i
  #   row_block[i]+1 = index of ΔQ / ΔV row for bus i
  row_block = zeros(Int, n)
  row = 0
  for i = 1:n
    if i == slack_idx
      continue
    end
    row += 2
    row_block[i] = row - 1  # ΔP row
  end

  # Triplet storage
  Iidx = Int[]
  Jidx = Int[]
  Vals = Float64[]
  # Rough capacity hint (purely heuristic)
  sizehint!(Iidx, 16 * nnz(Ybus))
  sizehint!(Jidx, 16 * nnz(Ybus))
  sizehint!(Vals, 16 * nnz(Ybus))

  # Helpers for sparse access
  rv    = rowvals(Ybus)
  nzval = nonzeros(Ybus)

  # --- 1) Contributions from complex power equations (ΔP, ΔQ) ---------------
  #
  # For bus i, bus j:
  #   ∂S_i/∂Vr_j = conj(I_i) * δ_ij + V_i * conj(Y_ij)
  #   ∂S_i/∂Vi_j = j*(conj(I_i) * δ_ij - V_i * conj(Y_ij))
  #
  # Then:
  #   ∂P_i/∂Vr_j = Re(∂S_i/∂Vr_j)
  #   ∂P_i/∂Vi_j = Re(∂S_i/∂Vi_j)
  #   ∂Q_i/∂Vr_j = Im(∂S_i/∂Vr_j)
  #   ∂Q_i/∂Vi_j = Im(∂S_i/∂Vi_j)
  #
  # For PV buses, only the ΔP row uses these derivatives; the second row is ΔV.

  for j = 1:n
    col_pos = pos_non_slack[j]
    if col_pos == 0
      # Slack bus column -> no state variable
      continue
    end

    colVr = col_pos
    colVi = (n - 1) + col_pos

    for ptr in nzrange(Ybus, j)
      i = rv[ptr]

      # Slack bus has no equations
      if i == slack_idx
        continue
      end

      rb = row_block[i]
      if rb == 0
        continue
      end

      rowP = rb          # ΔP row for bus i
      rowQ = rb + 1      # ΔQ/ΔV row for bus i

      Yij = nzval[ptr]

      # Base contributions from V_i * conj(Y_ij)
      S_dVr = V[i] * conj(Yij)
      S_dVi = im * (-V[i] * conj(Yij))

      # Diagonal term from conj(I_i) * δ_ij
      if i == j
        Ii = I[i]
        S_dVr += conj(Ii)
        S_dVi += im * conj(Ii)
      end

      # Real / imaginary parts
      ∂P_Vr = real(S_dVr)
      ∂P_Vi = real(S_dVi)
      ∂Q_Vr = imag(S_dVr)
      ∂Q_Vi = imag(S_dVi)

      # First equation: always ΔP_i for PQ and PV
      if abs(∂P_Vr) > 0.0
        push!(Iidx, rowP)
        push!(Jidx, colVr)
        push!(Vals, ∂P_Vr)
      end
      if abs(∂P_Vi) > 0.0
        push!(Iidx, rowP)
        push!(Jidx, colVi)
        push!(Vals, ∂P_Vi)
      end

      # Second equation:
      #   - PQ: ΔQ_i -> uses Q derivatives
      #   - PV: ΔV_i -> no contribution from S, handled separately
      bt = bus_types[i]
      if bt == :PQ
        if abs(∂Q_Vr) > 0.0
          push!(Iidx, rowQ)
          push!(Jidx, colVr)
          push!(Vals, ∂Q_Vr)
        end
        if abs(∂Q_Vi) > 0.0
          push!(Iidx, rowQ)
          push!(Jidx, colVi)
          push!(Vals, ∂Q_Vi)
        end
      elseif bt == :PV
        # nothing here, ΔV row added below
      else
        error("build_rectangular_jacobian_pq_pv_sparse: unsupported bus type $(bt) at bus $i")
      end
    end
  end

  # --- 2) Contributions for ΔV_i = |V_i| - Vset[i] on PV buses --------------
  #
  # Only depends on local Vr_i, Vi_i:
  #   |V_i| = sqrt(Vr_i^2 + Vi_i^2)
  #   ∂|V_i|/∂Vr_i = Vr_i / |V_i|
  #   ∂|V_i|/∂Vi_i = Vi_i / |V_i|

  for i = 1:n
    if i == slack_idx || bus_types[i] != :PV
      continue
    end

    rb   = row_block[i]
    rowV = rb + 1  # second row for that bus

    pos = pos_non_slack[i]
    if pos == 0
      continue
    end

    vm = abs(V[i])
    if vm == 0.0
      continue
    end

    dVr = real(V[i]) / vm
    dVi = imag(V[i]) / vm

    colVr = pos
    colVi = (n - 1) + pos

    if abs(dVr) > 0.0
      push!(Iidx, rowV)
      push!(Jidx, colVr)
      push!(Vals, dVr)
    end
    if abs(dVi) > 0.0
      push!(Iidx, rowV)
      push!(Jidx, colVi)
      push!(Vals, dVi)
    end
  end

  # Local chain-rule terms for voltage-dependent specified injections.
  # For ΔP = Pcalc - Pspec(|V|): subtract dPspec/d|V| * d|V|/dVr and d|V|/dVi.
  # For PQ second row ΔQ = Qcalc - Qspec(|V|): analogous subtraction.
  for i in non_slack
    pos = pos_non_slack[i]
    rb = row_block[i]
    pos == 0 && continue
    rb == 0 && continue

    vm = abs(V[i])
    vm_safe = vm > vm_eps ? vm : vm_eps
    dvm_dvr = real(V[i]) / vm_safe
    dvm_dvi = imag(V[i]) / vm_safe
    colVr = pos
    colVi = (n - 1) + pos

    dP = dPinj_dVm[i]
    if dP != 0.0
      push!(Iidx, rb)
      push!(Jidx, colVr)
      push!(Vals, -dP * dvm_dvr)
      push!(Iidx, rb)
      push!(Jidx, colVi)
      push!(Vals, -dP * dvm_dvi)
    end

    if bus_types[i] == :PQ
      dQ = dQinj_dVm[i]
      if dQ != 0.0
        push!(Iidx, rb + 1)
        push!(Jidx, colVr)
        push!(Vals, -dQ * dvm_dvr)
        push!(Iidx, rb + 1)
        push!(Jidx, colVi)
        push!(Vals, -dQ * dvm_dvi)
      end
    end
  end

  return sparse(Iidx, Jidx, Vals, m, nvar)
end

"""
    build_rectangular_jacobian_pq_pv_dense(
        Ybus, V, bus_types, Vset, slack_idx
    ) -> J::Matrix{Float64}

Build the analytic rectangular Jacobian for the mismatch vector `F(V)`
defined in `mismatch_rectangular`.

- State vector: x = [Vr(non-slack); Vi(non-slack)]
- Rows: for each non-slack bus i
    * PQ: [ΔP_i; ΔQ_i]
    * PV: [ΔP_i; ΔV_i]  with ΔV_i = |V_i| - Vset[i]

`bus_types` and `Vset` must be consistent with `mismatch_rectangular`.
"""
function build_rectangular_jacobian_pq_pv_dense(Ybus, V::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int; dPinj_dVm::Vector{Float64} = zeros(Float64, length(V)), dQinj_dVm::Vector{Float64} = zeros(Float64, length(V)), vm_eps::Float64 = 1e-9)
  n = length(V)
  @assert length(bus_types) == n
  @assert length(Vset) == n
  @assert length(dPinj_dVm) == n
  @assert length(dQinj_dVm) == n

  # --- 1) Wirtinger blocks for S(V) = V .* conj(Ybus * V)
  J11, J12, J21, J22 = build_complex_jacobian(Ybus, V)

  # --- 2) Full 2n×2n rectangular J for ΔP/ΔQ wrt Vr/Vi (all buses)
  # Rows: [ΔP_1..ΔP_n; ΔQ_1..ΔQ_n]
  # Cols: [Vr_1..Vr_n; Vi_1..Vi_n]
  Jrect_full = zeros(Float64, 2n, 2n)

  @inbounds for j = 1:n
    col_sum  = J11[:, j] .+ J12[:, j]  # corresponds to dS/dVr_j
    col_diff = J11[:, j] .- J12[:, j]  # used for dS/dVi_j

    # dP/dVr_j, dQ/dVr_j
    @views Jrect_full[1:n, j]      .= real.(col_sum)
    @views Jrect_full[(n+1):2n, j] .= imag.(col_sum)

    # dP/dVi_j, dQ/dVi_j
    @views Jrect_full[1:n, n+j]      .= -imag.(col_diff)
    @views Jrect_full[(n+1):2n, n+j] .= real.(col_diff)
  end

  # --- 3) Reduce to non-slack variables and rows matching mismatch_rectangular

  non_slack = collect(1:n)
  deleteat!(non_slack, slack_idx)

  nvar = 2 * (n - 1)
  m    = 2 * (n - 1)
  @assert nvar == m

  pos_map = build_pos_map(non_slack, n)

  # Column indices in the full rectangular J that correspond to
  # [Vr(non_slack); Vi(non_slack)]
  col_idx_full = vcat(non_slack, n .+ non_slack)

  J = zeros(Float64, m, nvar)

  row = 1
  @inbounds for i = 1:n
    if i == slack_idx
      continue
    end

    # First row for this bus: ΔP_i
    rowP_full = i                  # P row index in full J
    @views J[row, :] .= Jrect_full[rowP_full, col_idx_full]

    # Second row: ΔQ_i (PQ) or ΔV_i (PV)
    if bus_types[i] == :PQ
      rowQ_full = n + i          # Q row index in full J
      @views J[row+1, :] .= Jrect_full[rowQ_full, col_idx_full]

    elseif bus_types[i] == :PV
      # ΔV_i = |V_i| - Vset[i]
      J[row+1, :] .= 0.0

      pos = pos_map[i]
      if pos != 0
        vm = abs(V[i])
        if vm > 0.0
          dVr = real(V[i]) / vm
          dVi = imag(V[i]) / vm

          # Columns in reduced J:
          #   Vr_i -> index pos
          #   Vi_i -> index (n-1) + pos
          J[row+1, pos]       = dVr
          J[row+1, (n-1)+pos] = dVi
        end
      end
    else
      error("build_rectangular_jacobian_pq_pv: unsupported bus type $(bus_types[i]) at bus $i")
    end

    row += 2
  end

  for i in non_slack
    rowP = 2 * pos_map[i] - 1
    pos = pos_map[i]
    vm = abs(V[i])
    vm_safe = vm > vm_eps ? vm : vm_eps
    dvm_dvr = real(V[i]) / vm_safe
    dvm_dvi = imag(V[i]) / vm_safe

    J[rowP, pos] -= dPinj_dVm[i] * dvm_dvr
    J[rowP, (n-1)+pos] -= dPinj_dVm[i] * dvm_dvi

    if bus_types[i] == :PQ
      rowQ = rowP + 1
      J[rowQ, pos] -= dQinj_dVm[i] * dvm_dvr
      J[rowQ, (n-1)+pos] -= dQinj_dVm[i] * dvm_dvi
    end
  end

  return J
end

"""
    build_rectangular_jacobian_pq_pv(
        Ybus,
        V,
        bus_types,
        Vset,
        slack_idx;
        use_sparse::Bool = false,
    )

Dispatches to either the dense or sparse rectangular Jacobian builder matching
`mismatch_rectangular`.

- If `use_sparse == true` and `Ybus` is a `SparseMatrixCSC{ComplexF64}`, the
  sparse builder is used.
- Otherwise, the dense builder is used.
"""
function build_rectangular_jacobian_pq_pv(Ybus, V::Vector{ComplexF64}, bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int; use_sparse::Bool = false, dPinj_dVm::Vector{Float64} = zeros(Float64, length(V)), dQinj_dVm::Vector{Float64} = zeros(Float64, length(V)), vm_eps::Float64 = 1e-9)
  if use_sparse && Ybus isa SparseMatrixCSC{ComplexF64}
    return build_rectangular_jacobian_pq_pv_sparse(Ybus, V, bus_types, Vset, slack_idx; dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm, vm_eps = vm_eps)
  else
    return build_rectangular_jacobian_pq_pv_dense(Ybus, V, bus_types, Vset, slack_idx; dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm, vm_eps = vm_eps)
  end
end

"""
    complex_newton_step_rectangular(
        Ybus,
        V,
        S;
        slack_idx,
        damp,
        bus_types,
        Vset,
        use_sparse=false,
    )

Performs one Newton–Raphson step in rectangular coordinates using the analytic
Jacobian that matches `mismatch_rectangular`.

- State: x = [Vr(non-slack); Vi(non-slack)]
- Residual: F(x) = mismatch_rectangular(...)
"""
function complex_newton_step_rectangular(
  Ybus,
  V::Vector{ComplexF64},
  S::Vector{ComplexF64};
  slack_idx::Int,
  damp::Float64 = 1.0,
  autodamp::Bool = false,
  autodamp_min::Float64 = 1e-3,
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  use_sparse::Bool = false,
  dPinj_dVm::Vector{Float64} = zeros(Float64, length(V)),
  dQinj_dVm::Vector{Float64} = zeros(Float64, length(V)),
)
  n = length(V)
  @assert length(S) == n
  @assert length(bus_types) == n
  @assert length(Vset) == n
  @assert length(dPinj_dVm) == n
  @assert length(dQinj_dVm) == n

  # Non-slack indices
  non_slack = non_slack_indices(n, slack_idx)
  # Residual matching the FD variant
  F0 = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
  m = length(F0)
  nvar = 2 * (n - 1)
  @assert m == nvar "complex_newton_step_rectangular: mismatch and state dimension differ"

  # Analytic Jacobian (dense or sparse)
  J = build_rectangular_jacobian_pq_pv(Ybus, V, bus_types, Vset, slack_idx; use_sparse = use_sparse, dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm)

  # Solve J * δx = -F
  δx = solve_linear(J, -F0; allow_pinv = true)
  if autodamp
    _, Vtrial, _ = choose_rectangular_autodamp(Ybus, V, S, δx, F0; slack_idx = slack_idx, damp = damp, autodamp_min = autodamp_min, bus_types = bus_types, Vset = Vset)
    return Vtrial
  end

  _validate_rectangular_damping(damp, min(autodamp_min, damp))
  return _apply_rectangular_delta(V, δx, slack_idx, non_slack, damp)
end

"""
    run_complex_nr_rectangular_for_net!(net; maxiter=20, tol=1e-8, damp=0.2, verbose=0, use_fd=false)

Run a complex-state Newton-Raphson power flow in rectangular coordinates on a Sparlectra network.

# Arguments
- `net::Net`: Network object containing bus, branch, and generation data
- `maxiter::Int=20`: Maximum number of Newton-Raphson iterations
- `tol::Float64=1e-8`: Convergence tolerance for maximum mismatch
- `damp::Float64=0.2`: Damping factor for Newton step (0 < damp ≤ 1)
- `verbose::Int=0`: Verbosity level (0=quiet, 1=basic info, 2=detailed)
- `use_fd::Bool=false`: Use finite-difference Jacobian instead of analytic

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
- `runpf_rectangular!()`: Convenience wrapper matching `runpf!()` signature
- `mismatch_rectangular()`: Core mismatch function for PQ/PV constraints
- `build_rectangular_jacobian_pq_pv()`: Analytic Jacobian construction
"""

function run_complex_nr_rectangular_for_net!(
  net::Net;
  maxiter::Int = 20,
  tol::Float64 = 1e-8,
  damp::Float64 = 0.2,
  verbose::Int = 0,
  use_fd::Bool = false,
  autodamp::Bool = false,
  autodamp_min::Float64 = 1e-3,
  opt_sparse::Bool = true,
  opt_flatstart::Bool = net.flatstart,
  pv_table_rows::Int = 30,
  lock_pv_to_pq_buses::AbstractVector{Int} = Int[],
  qlimit_mode::Symbol = :switch_to_pq,
  qlimit_max_outer::Int = 30,
  start_projection::Bool = false,
  start_projection_try_dc_start::Bool = true,
  start_projection_try_blend_scan::Bool = true,
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
)
  if verbose > 1
    @info "Running complex rectangular NR power flow... use_fd=$use_fd, opt_sparse=$opt_sparse"
  end

  nodes = net.nodeVec
  n     = length(nodes)
  Sbase = net.baseMVA
  if !opt_sparse
    sparse = n > 1000
  else
    sparse = opt_sparse
  end
  Yred = createYBUS(net = net, sparse = sparse, printYBUS = (verbose > 1))
  Ybus = (size(Yred, 1) == n) ? Yred : _expand_ybus_for_isolated_nodes(Yred, n, net.isoNodes)

  # 1) Initial complex voltages V0 and slack index
  V0, slack_idx = initialVrect(net; flatstart = opt_flatstart)

  # 2) Specified complex power injections S (p.u.). For Q(U)/P(U) controllers
  # this vector becomes state-dependent and is re-evaluated per Newton iteration.
  S = buildComplexSVec(net)
  dPinj_dVm = zeros(Float64, n)
  dQinj_dVm = zeros(Float64, n)
  has_vdep_control = has_voltage_dependent_control(net)

  # 3) Bus types from Node data, and PV setpoints from regulating prosumers.
  # Node voltages may be temporary start guesses for MATPOWER flat-start modes.
  bus_types = Vector{Symbol}(undef, n)
  Vset      = _bus_voltage_setpoints_from_prosumers(net)

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
      error("run_complex_nr_rectangular_for_net!: unsupported bus type at bus $k, given: $(BusType)")
    end

  end

  # The slack/reference voltage is fixed rather than solved by a residual row.
  # Keep its magnitude at the regulating prosumer setpoint even if a MATPOWER
  # flat-start mode temporarily changed node._vm_pu as an initial guess.
  V0[slack_idx] = ComplexF64(Vset[slack_idx] * cos(angle(V0[slack_idx])), Vset[slack_idx] * sin(angle(V0[slack_idx])))

  V0 = project_rectangular_start(
    Ybus,
    V0,
    S,
    bus_types,
    Vset,
    slack_idx;
    enabled = start_projection,
    try_dc_start = start_projection_try_dc_start,
    try_blend_scan = start_projection_try_blend_scan,
    blend_lambdas = start_projection_blend_lambdas,
    dc_angle_limit_deg = start_projection_dc_angle_limit_deg,
    verbose = verbose,
  )

  # 4) Q-limit data 
  qmin_pu, qmax_pu = getQLimits_pu(net)
  # Start fresh each PF run before guard pre-processing records locked buses.
  resetQLimitLog!(net)
  if verbose > 1
    printPVQLimitsTable(net; max_rows = typemax(Int))
  elseif verbose > 0
    printPVQLimitsTable(net; max_rows = pv_table_rows)
  end

  guarded_qlimit_buses = Int[]
  if qlimit_guard
    guarded_qlimit_buses = _apply_qlimit_guard_to_rectangular_active_set!(
      net,
      bus_types,
      S,
      build_qload_pu(net),
      qmin_pu,
      qmax_pu;
      min_q_range_pu = qlimit_guard_min_q_range_pu,
      zero_range_mode = qlimit_guard_zero_range_mode,
      narrow_range_mode = qlimit_guard_narrow_range_mode,
      log = qlimit_guard_log,
      verbose = verbose,
    )
  end

  # --- Active-set bookkeeping (rectangular solver) ------------------------
  nb = n  # number of buses

  # PV origin mask (guards for PQ->PV re-enable)
  pv_orig_mask = falses(nb)
  @inbounds for k = 1:nb
    pv_orig_mask[k] = (bus_types[k] == :PV)
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

  controllers = qlimit_mode == :adjust_vset ? _build_vset_adjust_controllers(net) : Dict{Int,NamedTuple{(:prosumer_idx, :config),Tuple{Int,VoltageAdjustConfig}}}()
  base_vset = copy(Vset)
  adjust_counter = zeros(Int, nb)
  qlimit_trace_internal = _resolve_qlimit_trace_buses(net, qlimit_trace_buses)
  qlimit_trace_enabled = !isempty(qlimit_trace_internal)
  if qlimit_trace_enabled
    missing = setdiff(collect(qlimit_trace_buses), [_qlimit_original_bus_id(net, bus) for bus in qlimit_trace_internal])
    isempty(missing) || @warn "qlimit_trace_buses entries not found in network" missing = missing
    println("Q-limit trace enabled for BUS_I values: ", [_qlimit_original_bus_id(net, bus) for bus in qlimit_trace_internal])
  end

  # 5) NR-Loop
  V         = copy(V0)
  history   = Float64[]
  converged = false
  iters     = 0
  rejection_reason = :nr_mismatch_not_converged
  prev_pv_qreq_pu = fill(NaN, nb)

  if verbose > 1
    @info "Starting rectangular complex NR power flow..."
    @info "Initial complex voltages V0:" V0
    @info "Slack bus index:" slack_idx
    @info "maxiter = $maxiter, tol = $tol, damp = $damp, autodamp = $autodamp, autodamp_min = $autodamp_min, start_projection = $start_projection"
  end

  for it = 1:maxiter
    iters = it

    if has_vdep_control
      S, dPinj_dVm, dQinj_dVm = buildControlledSVec(net, V)
    end

    # Mismatch with current bus_types and (possibly) voltage-dependent S.
    F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
    max_mis = maximum(abs.(F))
    push!(history, max_mis)

    (verbose > 1) && @debug "Rectangular NR iteration" iter = it max_mismatch = max_mis

    # --- Q-Limit Active Set: PV -> PQ, optional PQ -> PV (rectangular) ------
    changed   = false
    reenabled = false

    Qload_pu = build_qload_pu(net)

    Scalc_pu = calc_injections(Ybus, V)
    current_pv_qreq_pu = fill(NaN, nb)
    @inbounds for bus in eachindex(current_pv_qreq_pu)
      if bus_types[bus] == :PV
        current_pv_qreq_pu[bus] = imag(Scalc_pu[bus]) + Qload_pu[bus]
      end
    end

    qlimit_iter_ready = it >= qlimit_start_iter
    qlimit_auto_ready = false
    if qlimit_start_mode in (:auto_q_delta, :iteration_or_auto)
      max_q_delta = 0.0
      compared = false
      @inbounds for bus in eachindex(current_pv_qreq_pu)
        if isfinite(current_pv_qreq_pu[bus]) && isfinite(prev_pv_qreq_pu[bus])
          max_q_delta = max(max_q_delta, abs(current_pv_qreq_pu[bus] - prev_pv_qreq_pu[bus]))
          compared = true
        end
      end
      qlimit_auto_ready = compared && (max_q_delta <= qlimit_auto_q_delta_pu)
    end
    qlimit_ready = qlimit_start_mode == :iteration ? qlimit_iter_ready :
                   qlimit_start_mode == :auto_q_delta ? qlimit_auto_ready :
                   (qlimit_iter_ready || qlimit_auto_ready)
    converged_this_iter = max_mis <= tol
    violation_guard_active = qlimit_guard_violation_mode == :lock_pq
    qlimit_check_active = qlimit_ready && (!converged_this_iter || violation_guard_active)

    if qlimit_trace_enabled
      lock_mask_trace = falses(nb)
      for bus in lock_pv_to_pq_buses
        if 1 <= bus <= nb
          lock_mask_trace[bus] = true
        end
      end
      for bus in qlimit_trace_internal
        bus <= nb || continue
        qreq = bus_types[bus] in (:PV, :Slack) ? imag(Scalc_pu[bus]) + Qload_pu[bus] : NaN
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

    if qlimit_check_active
      changed, reenabled = active_set_q_limits!(
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
          return imag(Scalc_pu[bus]) + Qload_pu[bus]
        end,
        is_pv = bus -> (bus_types[bus] == :PV),
        make_pq! = (bus, qclamp_gen_pu, side) -> begin
          bus_types[bus] = :PQ
          qinj_pu = qclamp_gen_pu - Qload_pu[bus]
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
      )

      # If bus_types/spec changed, mismatch definition changed (ΔQ ↔ ΔV) => rebuild F
      if changed || reenabled
        F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
        max_mis = maximum(abs.(F))
        history[end] = max_mis  # optional: overwrite last stored value for this iteration
      end
    end
    if converged_this_iter && !(changed || reenabled)
      converged = true
      rejection_reason = :none
      break
    end
    prev_pv_qreq_pu = current_pv_qreq_pu
    # --- Newton step (FD or analytic) -----------------------------------
    try
      if use_fd
        V = complex_newton_step_rectangular_fd(Ybus, V, S; slack_idx = slack_idx, damp = damp, autodamp = autodamp, autodamp_min = autodamp_min, h = 1e-6, bus_types = bus_types, Vset = Vset, dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm)
      else
        V = complex_newton_step_rectangular(Ybus, V, S; slack_idx = slack_idx, damp = damp, autodamp = autodamp, autodamp_min = autodamp_min, bus_types = bus_types, Vset = Vset, use_sparse = sparse, dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm)
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
    V[slack_idx] = V0[slack_idx]
  end

  # 6) Update voltages back to network
  # --- mirror bus_types back into Net/node types (PV->PQ switching) ---
  @inbounds for k = 1:n
    if bus_types[k] == :PQ
      setNodeType!(net.nodeVec[k], "PQ")
    elseif bus_types[k] == :PV
      setNodeType!(net.nodeVec[k], "PV")
    end
  end
  numerical_converged = converged
  final_pv_voltage_residual = _max_rectangular_pv_voltage_residual(V, Vset, bus_types, net.qLimitEvents)
  if numerical_converged && final_pv_voltage_residual > tol
    verbose > 0 && @warn "Rectangular NR convergence rejected because active PV voltage setpoint residual exceeds tolerance." max_pv_voltage_residual = final_pv_voltage_residual tolerance = tol
    converged = false
    rejection_reason = :active_pv_voltage_residual
  end

  update_net_voltages_from_complex!(net, V)

  # 7) Compute bus injections from final voltages  
  Sbus_pu = calc_injections(Ybus, V)
  Sbus_MVA = Sbus_pu .* Sbase

  @debug "Final Voltages Mag = " [abs.(V)...]
  @debug "Final Voltages Ang = " [angle.(V) .* (180.0 / π)...]

  isoNodes = net.isoNodes

  @inbounds for (k, node) in enumerate(nodes)
    # Skip isolated nodes (if any)
    if k in isoNodes
      continue
    end

    Sbus      = Sbus_MVA[k]
    Pbus_MW   = real(Sbus)
    Qbus_MVar = imag(Sbus)

    # Slack bus: always write P and Q generation from the solved state
    if node._nodeType == Sparlectra.Slack
      node._pƩGen = Pbus_MW
      node._qƩGen = Qbus_MVar

    elseif node._nodeType == Sparlectra.PV
      # PV bus: Qgen is a result (unless it later switches)
      node._qƩGen = Qbus_MVar

    elseif node._nodeType == Sparlectra.PQ
      # If this bus was forced PV->PQ by Q-limit, keep the clamped generator Q
      # (already written in make_pq! as MVar).
      if haskey(net.qLimitEvents, k)
        @debug "Bus $(k) is PQ due to Q-limit; keeping clamped Qgen = $(node._qƩGen) MVar."
      end
    end
    # PQ buses / pure loads: do not touch _pƩLoad / _pƩGen here.
    # The original load/generation specification remains intact.
  end

  # 8) Update total bus power (sum of complex injections in p.u.)
  p = (sum(real.(Sbus_pu))) * Sbase
  q = (sum(imag.(Sbus_pu))) * Sbase

  if verbose > 1
    @info "Set total bus power to p = $p MW and q = $q MVar"
  end

  setTotalBusPower!(net = net, p = p, q = q)
  updateShuntPowers!(net = net)

  qlimit_summary = nothing
  if numerical_converged
    final_Qload_pu = build_qload_pu(net)
    qlimit_summary_io = (verbose > 0 || qlimit_trace_enabled) ? stdout : devnull
    qlimit_summary = _print_rectangular_qlimit_summary(qlimit_summary_io, net, V, Sbus_pu, bus_types, qmin_pu, qmax_pu, final_Qload_pu; q_hyst_pu = q_hyst_pu, tolerance_pu = tol, max_rows = pv_table_rows, max_console_rows = pv_table_rows)
    remaining_pv_violations = qlimit_summary.pv_violations
    bounded_ok = qlimit_guard_accept_bounded_violations && remaining_pv_violations <= qlimit_guard_max_remaining_violations
    if remaining_pv_violations > 0 && !bounded_ok
      verbose > 0 && @warn "Rectangular NR active-set failed because active PV Q-limit violations remain after the numerical solve." pv_violations = qlimit_summary.pv_violations ref_violations = qlimit_summary.ref_violations
      converged = false
      rejection_reason = :remaining_pv_q_limit_violations
    elseif remaining_pv_violations > 0 && bounded_ok
      rejection_reason = :bounded_q_limit_violations_accepted
    end
  end

  switch_counts = qlimit_switch_counts(net)
  oscillating_buses = count(>=(max(qlimit_guard_max_switches, 1)), values(switch_counts))
  max_switching_exceeded = qlimit_guard_freeze_after_repeated_switching && oscillating_buses > 0
  q_limit_active_set_ok = numerical_converged && final_pv_voltage_residual <= tol && (isnothing(qlimit_summary) || qlimit_summary.pv_violations == 0 || (qlimit_guard_accept_bounded_violations && qlimit_summary.pv_violations <= qlimit_guard_max_remaining_violations)) && !max_switching_exceeded
  if numerical_converged && max_switching_exceeded && !q_limit_active_set_ok
    rejection_reason = :max_switching_exceeded
    converged = false
  end
  final_reason = converged ? :none : rejection_reason
  final_status = _rectangular_solver_status_symbol(numerical_converged, q_limit_active_set_ok, converged, final_reason)
  status = _set_rectangular_pf_status!(net, (
    numerical_converged = numerical_converged,
    nr_converged = numerical_converged,
    active_set_converged = q_limit_active_set_ok,
    q_limit_active_set_ok = q_limit_active_set_ok,
    final_converged = converged,
    status = final_status,
    reason = final_reason,
    reason_text = _rectangular_rejection_reason_text(final_reason),
    pv_q_limit_violations = isnothing(qlimit_summary) ? 0 : qlimit_summary.pv_violations,
    ref_q_limit_violations = isnothing(qlimit_summary) ? 0 : qlimit_summary.ref_violations,
    final_pv_voltage_residual = final_pv_voltage_residual,
    final_mismatch = isempty(history) ? Inf : history[end],
    pv_pq_switching_events = length(net.qLimitLog),
    oscillating_buses = oscillating_buses,
    guarded_narrow_q_pv_buses = length(guarded_qlimit_buses),
  ))
  if verbose > 0 || (numerical_converged && !q_limit_active_set_ok)
    _print_rectangular_convergence_summary(stdout, status)
    _print_qlimit_active_set_summary(stdout, status)
  end

  return iters, converged ? 0 : 1
end

"""
    runpf_rectangular!(net, maxIte, tolerance=1e-6, verbose=0)

Runs a rectangular complex-state Newton–Raphson power flow on `net::Net`.

Returns:
    (iterations::Int, status::Int)
where `status == 0` indicates convergence.
"""
function runpf_rectangular!(
  net::Net,
  maxIte::Int,
  tolerance::Float64 = 1e-6,
  verbose::Int = 0;
  opt_fd::Bool = false,
  opt_sparse::Bool = true,
  damp = 1.0,
  autodamp::Bool = false,
  autodamp_min::Float64 = 1e-3,
  opt_flatstart::Bool = net.flatstart,
  pv_table_rows::Int = 30,
  lock_pv_to_pq_buses::AbstractVector{Int} = Int[],
  qlimit_mode::Symbol = :switch_to_pq,
  qlimit_max_outer::Int = 30,
  start_projection::Bool = false,
  start_projection_try_dc_start::Bool = true,
  start_projection_try_blend_scan::Bool = true,
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
)
  iters, erg = run_complex_nr_rectangular_for_net!(
    net;
    maxiter = maxIte,
    tol = tolerance,
    damp = damp,
    autodamp = autodamp,
    autodamp_min = autodamp_min,
    verbose = verbose,
    use_fd = opt_fd,
    opt_sparse = opt_sparse,
    opt_flatstart = opt_flatstart,
    pv_table_rows = pv_table_rows,
    lock_pv_to_pq_buses = lock_pv_to_pq_buses,
    qlimit_mode = qlimit_mode,
    qlimit_max_outer = qlimit_max_outer,
    start_projection = start_projection,
    start_projection_try_dc_start = start_projection_try_dc_start,
    start_projection_try_blend_scan = start_projection_try_blend_scan,
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
  )
  return iters, erg
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

  return [find_root(i) for i = 1:n]
end

function _merged_pf_net(net::Net)
  reps = _active_link_representative_map(net)
  all(reps[i] == i for i in eachindex(reps)) && return net, reps, false

  wnet = deepcopy(net)
  n = length(wnet.nodeVec)
  cluster_members = [Int[] for _ = 1:n]
  for bus = 1:n
    push!(cluster_members[reps[bus]], bus)
  end

  for rep = 1:n
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

"""
    runpf!(net, maxIte, tolerance=1e-6, verbose=0; method=:rectangular)

Unified AC power flow interface.

Arguments:
- `net::Net`: network
- `maxIte::Int`: maximum iterations
- `tolerance::Float64`: mismatch tolerance
- `verbose::Int`: verbosity level
- `method::Symbol`: `:rectangular` (recommended), `:polar_full` (deprecated), or `:classic` (deprecated)
- `autodamp::Bool`: enable residual-based backtracking for rectangular Newton steps
- `autodamp_min::Float64`: minimum automatic damping factor when `autodamp = true`
- `qlimit_start_iter::Int`: first Newton iteration where PV→PQ Q-limit switching may run in `:iteration` mode
- `qlimit_start_mode::Symbol`: `:iteration`, `:auto_q_delta`, or `:iteration_or_auto` start criterion for PV→PQ switching
- `qlimit_auto_q_delta_pu::Float64`: PV reactive-power request change threshold for automatic switching start

Notes:
- Link-flow recovery (`calcLinkFlowsKCL!`) is method-agnostic and uses solved PF results.
- If active-link merges create internal isolated buses, `:rectangular` currently falls
  back to `:polar_full` for robustness.

Returns:
    (iterations::Int, status::Int)

where `status == 0` indicates convergence.
"""
function runpf!(
  net::Net,
  maxIte::Int,
  tolerance::Float64 = 1e-6,
  verbose::Int = 0;
  method::Symbol = :rectangular,
  opt_fd::Bool = false,
  opt_sparse::Bool = true,
  opt_flatstart::Bool = net.flatstart,
  damp = 1.0,
  autodamp::Bool = false,
  autodamp_min::Float64 = 1e-3,
  pv_table_rows::Int = 30,
  validate_limits_after_pf::Bool = false,
  q_limit_violation_headroom::Float64 = 0.0,
  lock_pv_to_pq_buses::AbstractVector{Int} = Int[],
  qlimit_mode::Symbol = :switch_to_pq,
  qlimit_max_outer::Int = 30,
  start_projection::Bool = false,
  start_projection_try_dc_start::Bool = true,
  start_projection_try_blend_scan::Bool = true,
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
)
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

  #@info "Running AC Power Flow using method: $(method)"
  if method === :polar_full
    has_vdep_control && error("runpf!: voltage-dependent injections, including P(U)/Q(U) controllers and bus_shunt_model=voltage_dependent_injection, are currently supported only for method=:rectangular.")
    if qlimit_mode != :switch_to_pq
      @warn "runpf!: qlimit_mode=$(qlimit_mode) is only supported for method=:rectangular. Falling back to :switch_to_pq behavior."
    end
    iters, erg = runpf_full!(wnet, maxIte, tolerance, verbose; opt_sparse = opt_sparse, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, lock_pv_to_pq_buses = lock_pv_to_pq_buses, warn_deprecated = true)
    if erg == 0 && has_merges
      _sync_merged_results_to_original!()
    end
    if validate_limits_after_pf && (verbose > 0)
      printFinalLimitValidation(has_merges ? net : wnet; q_headroom = q_limit_violation_headroom)
    end
    return iters, erg
  elseif method === :rectangular
    if has_merges
      has_vdep_control && error("runpf!: voltage-dependent injections, including P(U)/Q(U) controllers and bus_shunt_model=voltage_dependent_injection, are not supported with active-link merge handling in rectangular mode. Disable merges or use a topology without internal isolated buses.")
      if verbose > 0
        @warn "runpf!: rectangular solver detected internal Isolated buses from active-link merges; using rectangular FD fallback instead of :polar_full"
      end
      iters, erg = runpf_rectangular!(wnet, maxIte, tolerance, verbose; opt_fd = true, opt_sparse = opt_sparse, damp = damp, autodamp = autodamp, autodamp_min = autodamp_min, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, lock_pv_to_pq_buses = lock_pv_to_pq_buses, qlimit_mode = qlimit_mode, qlimit_max_outer = qlimit_max_outer, start_projection = start_projection, start_projection_try_dc_start = start_projection_try_dc_start, start_projection_try_blend_scan = start_projection_try_blend_scan, start_projection_blend_lambdas = start_projection_blend_lambdas, start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg, qlimit_start_iter = qlimit_start_iter, qlimit_start_mode = qlimit_start_mode, qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu, qlimit_trace_buses = qlimit_trace_buses, qlimit_lock_reason = qlimit_lock_reason, qlimit_guard = qlimit_guard, qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu, qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode, qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode, qlimit_guard_log = qlimit_guard_log, qlimit_guard_max_switches = qlimit_guard_max_switches, qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations, qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations, qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching, qlimit_guard_violation_mode = qlimit_guard_violation_mode, qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu)
    else
      iters, erg = runpf_rectangular!(wnet, maxIte, tolerance, verbose; opt_fd = opt_fd, opt_sparse = opt_sparse, damp = damp, autodamp = autodamp, autodamp_min = autodamp_min, opt_flatstart = opt_flatstart, pv_table_rows = pv_table_rows, lock_pv_to_pq_buses = lock_pv_to_pq_buses, qlimit_mode = qlimit_mode, qlimit_max_outer = qlimit_max_outer, start_projection = start_projection, start_projection_try_dc_start = start_projection_try_dc_start, start_projection_try_blend_scan = start_projection_try_blend_scan, start_projection_blend_lambdas = start_projection_blend_lambdas, start_projection_dc_angle_limit_deg = start_projection_dc_angle_limit_deg, qlimit_start_iter = qlimit_start_iter, qlimit_start_mode = qlimit_start_mode, qlimit_auto_q_delta_pu = qlimit_auto_q_delta_pu, qlimit_trace_buses = qlimit_trace_buses, qlimit_lock_reason = qlimit_lock_reason, qlimit_guard = qlimit_guard, qlimit_guard_min_q_range_pu = qlimit_guard_min_q_range_pu, qlimit_guard_zero_range_mode = qlimit_guard_zero_range_mode, qlimit_guard_narrow_range_mode = qlimit_guard_narrow_range_mode, qlimit_guard_log = qlimit_guard_log, qlimit_guard_max_switches = qlimit_guard_max_switches, qlimit_guard_accept_bounded_violations = qlimit_guard_accept_bounded_violations, qlimit_guard_max_remaining_violations = qlimit_guard_max_remaining_violations, qlimit_guard_freeze_after_repeated_switching = qlimit_guard_freeze_after_repeated_switching, qlimit_guard_violation_mode = qlimit_guard_violation_mode, qlimit_guard_violation_threshold_pu = qlimit_guard_violation_threshold_pu)
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
  elseif method === :classic
    has_vdep_control && error("runpf!: voltage-dependent injections, including P(U)/Q(U) controllers and bus_shunt_model=voltage_dependent_injection, are currently supported only for method=:rectangular.")
    if qlimit_mode != :switch_to_pq
      @warn "runpf!: qlimit_mode=$(qlimit_mode) is only supported for method=:rectangular. Falling back to :switch_to_pq behavior."
    end
    iters, erg = runpf_classic!(wnet, maxIte, tolerance, verbose, opt_sparse, opt_flatstart)
    if erg == 0 && has_merges
      _sync_merged_results_to_original!()
    end
    if validate_limits_after_pf && (verbose > 0)
      printFinalLimitValidation(has_merges ? net : wnet; q_headroom = q_limit_violation_headroom)
    end
    return iters, erg
  else
    error("runpf!: unknown method $(method). Use :rectangular (recommended), :polar_full (deprecated), or :classic (deprecated).")
  end
end

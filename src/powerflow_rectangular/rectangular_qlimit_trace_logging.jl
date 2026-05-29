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
# This file is included inside module Sparlectra. Do not add a module wrapper here.
#
# Rectangular power-flow Q-limit trace logging helpers.

# Date: 29.5.2026
# file: src/powerflow_rectangular/rectangular_qlimit_trace_logging.jl

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

  @printf(
    io,
    "Q-limit trace it=%d BUS_I=%d type=%s Vm_calc=%.6f imported_Vset=%.6f Qcalc=%.6f pu (%.3f MVAr) Qmin=%.6f pu (%.3f MVAr) Qmax=%.6f pu (%.3f MVAr) q_hyst_pu=%.6f qlimit_start_iter=%d qlimit_start_mode=%s active=%s violation=%s amount=%.6f pu (%.3f MVAr) decision=%s reason=%s\n",
    it,
    busI,
    String(bus_type),
    abs(V[bus]),
    Vset[bus],
    qreq_pu,
    qreq_pu * net.baseMVA,
    qmin,
    qmin * net.baseMVA,
    qmax,
    qmax * net.baseMVA,
    q_hyst_pu,
    qlimit_start_iter,
    String(qlimit_start_mode),
    string(qlimit_check_active),
    String(side),
    amount_pu,
    amount_pu * net.baseMVA,
    decision,
    reason
  )
  return nothing
end

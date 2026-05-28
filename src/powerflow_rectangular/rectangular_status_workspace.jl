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
# Rectangular power-flow status registry, iteration workspace, and status reporting helpers.

mutable struct _RectangularPFStatusTable
  entries::Vector{Tuple{UInt,WeakRef,Any}}
end

const _RECTANGULAR_PF_STATUS = _RectangularPFStatusTable(Tuple{UInt,WeakRef,Any}[])

mutable struct RectangularIterationWorkspace
  qload_pu::Vector{Float64}
  current_pv_qreq_pu::Vector{Float64}
  prev_pv_qreq_pu::Vector{Float64}
  lock_mask::BitVector
  rhs_vector::Vector{Float64}
end

function RectangularIterationWorkspace(nb::Int)
  return RectangularIterationWorkspace(zeros(Float64, nb), fill(NaN, nb), fill(NaN, nb), falses(nb), zeros(Float64, 2 * max(nb - 1, 0)))
end

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
  reason == :nr_mismatch_not_converged_active_set_unstable && return "NR mismatch did not converge; Q-limit active set changed repeatedly"
  reason == :wrong_branch_detected && return "wrong-branch solution detected"
  reason == :wrong_branch_rescue_not_implemented && return "wrong-branch rescue requested but not implemented"
  reason == :rescue_requested_but_not_available && return "wrong-branch rescue requested but not available"
  reason == :nonfinite_voltage && return "non-finite voltage state"
  reason == :nonfinite_branch_angle && return "non-finite branch angle state"
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
  numerical_converged || return reason == :singular_newton_step ? :singular_jacobian : :not_converged
  reason == :wrong_branch_detected && return :wrong_branch_detected
  reason == :wrong_branch_rescue_not_implemented && return :wrong_branch_rescue_not_implemented
  reason in (:nonfinite_voltage, :nonfinite_branch_angle) && return :wrong_branch_detected
  active_set_ok && return :converged_with_limit_warnings
  reason == :max_switching_exceeded && return :converged_limits_failed
  reason == :bounded_q_limit_violations_accepted && return :converged_with_limit_warnings
  reason == :remaining_pv_q_limit_violations && return :converged_limits_failed
  return :converged_limits_failed
end

function _print_qlimit_active_set_summary(io::IO, status)
  println(io, "==================== Q-Limit Active-Set Summary ====================")
  println(io)
  @printf(io, "NR convergence             : %s\n", status.numerical_converged ? "yes" : "no")
  @printf(io, "Final mismatch             : %.6g\n", status.final_mismatch)
  @printf(io, "Active-set convergence     : %s\n", status.q_limit_active_set_ok ? "yes" : "no")
  @printf(io, "PV→PQ switching events     : %d\n", status.pv_pq_switching_events)
  @printf(io, "Q-limit active-set changes : %d\n", status.qlimit_active_set_changes)
  @printf(io, "Q-limit re-enable events   : %d\n", status.qlimit_reenable_events)
  @printf(io, "Oscillating buses          : %d\n", status.oscillating_buses)
  @printf(io, "Guarded narrow-Q PV buses  : %d\n", status.guarded_narrow_q_pv_buses)
  @printf(io, "Final status               : %s\n", String(status.status))
  println(io)
  println(io, "===================================================================")
  return nothing
end

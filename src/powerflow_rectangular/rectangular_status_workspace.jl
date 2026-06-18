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

# Date: 29.5.2026
# file: src/powerflow_rectangular/rectangular_status_workspace.jl 

mutable struct _RectangularPFStatusTable
  # Weak-reference registry keeps per-Net rectangular status without preventing GC.
  # Allows solver to track state across multiple solves on the same Net object
  # without creating reference cycles or leaking memory when the Net is freed.
  entries::Vector{Tuple{UInt,WeakRef,Any}}
end

const _RECTANGULAR_PF_STATUS = _RectangularPFStatusTable(Tuple{UInt,WeakRef,Any}[])
# Global weak-reference registry for rectangular solver status.
# Indexed by (objectid, WeakRef) to maintain per-Net status across Julia sessions
# without preventing garbage collection of discarded networks.

mutable struct RectangularIterationWorkspace
  # Reusable vectors avoid per-iteration allocations in rectangular NR runs.
  # Allocated once at solve setup; cleared and reused in each outer (Q-limit) iteration.
  # This workspace is optional and is only created when rectangular_preallocate_workspace
  # is enabled or when matrix size exceeds the auto-threshold.
  qload_pu::Vector{Float64}
  current_pv_qreq_pu::Vector{Float64}
  prev_pv_qreq_pu::Vector{Float64}
  lock_mask::BitVector
  rhs_vector::Vector{Float64}
end

function RectangularIterationWorkspace(nb::Int)
  # Pre-allocate workspace buffers for a system with nb buses.
  # Vectors are sized to match per-bus data (nb elements) or reduced-state RHS (2(n-1) for non-slack).
  # Initialization with NaN helps detect accidental use of stale Q-requirement values.
  return RectangularIterationWorkspace(
    zeros(Float64, nb),           # qload_pu: injected/absorbed reactive load per bus
    fill(NaN, nb),                # current_pv_qreq_pu: PV Q request in this iteration
    fill(NaN, nb),                # prev_pv_qreq_pu: PV Q request in previous iteration (for delta tracking)
    falses(nb),                   # lock_mask: which buses are locked in active-set switching
    zeros(Float64, 2 * max(nb - 1, 0)),  # rhs_vector: RHS for reduced (n-1)-bus Newton system
  )
end

function _prune_rectangular_pf_status!(table::_RectangularPFStatusTable)
  # Dead weakrefs are removed opportunistically before any lookup/insert.
  # This keeps the global registry bounded even if many Net objects are created and freed.
  # Pruning is cheap (O(n) scan) and is performed before each status operation.
  filter!(entry -> entry[2].value !== nothing, table.entries)
  return table
end

function _set_rectangular_pf_status!(net::Net, status)
  # Store or update solver status for a specific Net object.
  # Uses objectid(net) + identity check (===) to distinguish between different Net objects,
  # even if they happen to have the same memory address after GC.
  # This prevents accidental status cross-contamination in stress tests or repeated solves.
  _prune_rectangular_pf_status!(_RECTANGULAR_PF_STATUS)
  key = objectid(net)
  for i in eachindex(_RECTANGULAR_PF_STATUS.entries)
    entry = _RECTANGULAR_PF_STATUS.entries[i]
    if entry[1] == key && entry[2].value === net
      # Found existing entry for this Net; update its status in-place.
      _RECTANGULAR_PF_STATUS.entries[i] = (key, WeakRef(net), status)
      return status
    end
  end
  # No existing entry; append new one to registry.
  push!(_RECTANGULAR_PF_STATUS.entries, (key, WeakRef(net), status))
  return status
end

"""
    rectangular_pf_status(net::Net) -> Any

Retrieve the most recent rectangular power-flow solver status for a network.

If no solver has run on this network, or the network has been garbage-collected,
returns `nothing`. Otherwise returns a status struct (type depends on solver
implementation) containing convergence flags, iteration counts, and diagnosis data.

# Purpose
Track solver outcome per Net object across multiple solve attempts without
creating reference cycles or preventing garbage collection. This enables
post-solve diagnostics, convergence checks, and iteration-count queries
without modifying the Net data structure.

# Arguments
- `net::Net`: Network object to query.

# Returns
- Status struct or `nothing` if no status is available.
"""
function rectangular_pf_status(net::Net)
  _prune_rectangular_pf_status!(_RECTANGULAR_PF_STATUS)
  key = objectid(net)
  for entry in _RECTANGULAR_PF_STATUS.entries
    entry[1] == key && entry[2].value === net && return entry[3]
  end
  return nothing
end

function _rectangular_rejection_reason_text(reason::Symbol)
  # Convert symbolic rejection reason to human-readable text.
  # Used in convergence summaries and solver status reports.
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
  # Fallback for unforeseen reasons: convert underscore to space.
  return replace(String(reason), "_" => " ")
end

"""
    _print_rectangular_convergence_summary(io::IO, status)

Write a one-line convergence summary to an output stream.

Displays:
- Numerical convergence (NR mismatch) status
- Q-limit active-set convergence status
- Final overall convergence flag
- Symbolic rejection/stop reason

# Purpose
Provide quick console/log summary of what halted the solver (convergence, max iter,
singular Jacobian, wrong-branch, or Q-limit oscillation). Intended for solver
report output, not for interactive diagnostics.

# Arguments
- `io::IO`: Output stream (typically `stdout` or a log file)
- `status`: Status struct with fields `numerical_converged`, `q_limit_active_set_ok`,
  `final_converged`, and `reason_text`.

# Side effects
- Writes formatted text to `io`.
"""
function _print_rectangular_convergence_summary(io::IO, status)
  # Compact one-line summary of convergence outcome.
  # Format: "rectangular convergence: numerical_solution=X  q_limit_active_set=Y  final_converged=Z  reason=..."
  numerical_text = status.numerical_converged ? "OK" : "FAIL"
  active_text = status.q_limit_active_set_ok ? "OK" : "FAIL"
  final_text = status.final_converged ? "true" : "false"
  @printf(io, "rectangular convergence: numerical_solution=%s  q_limit_active_set=%s  final_converged=%s  reason=%s\n", numerical_text, active_text, final_text, status.reason_text)
  return nothing
end

function _rectangular_solver_status_symbol(numerical_converged::Bool, active_set_ok::Bool, final_converged::Bool, reason::Symbol)::Symbol
  # Map detailed status flags to a single symbolic outcome for downstream reporting.
  # Keep symbol mapping stable: downstream status semantics/log text depend on it.
  # This avoids exposing internal flag details to logging, UI, and user APIs.
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

"""
    _print_qlimit_active_set_summary(io::IO, status)

Write a detailed Q-limit active-set convergence summary to an output stream.

Displays:
- Numerical NR convergence
- Final mismatch value
- Active-set convergence status
- Switching event counters (PV→PQ, active-set changes, re-enable events)
- Count of oscillating buses and guarded narrow-Q buses
- Final solver status symbol

# Purpose
Provide a structured diagnostic report for Q-limit iterations when the
rectangular solver has run. Used in verbose output and detailed solver
convergence analysis. Helps users understand whether Q-limit switching
was the limiting factor or whether the NR mismatch itself failed to converge.

# Arguments
- `io::IO`: Output stream (typically `stdout` or a log file)
- `status`: Status struct with fields for convergence flags, mismatch, event counts,
  and final status symbol.

# Side effects
- Writes formatted table to `io`.
"""
function _print_qlimit_active_set_summary(io::IO, status)
  # Detailed multi-line table of Q-limit iteration diagnostics.
  # Helps distinguish between "NR solved but Q-limits unsatisfied" vs.
  # "active-set oscillation prevented convergence" vs. "fully converged".
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

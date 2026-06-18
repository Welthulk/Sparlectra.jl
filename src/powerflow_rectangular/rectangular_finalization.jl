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
# Rectangular power-flow post-iteration finalization helpers.
# Date: 29.5.2026
# file: src/powerflow_rectangular/rectangular_finalization.jl

function _sync_rectangular_bus_types_to_net!(net::Net, bus_types)
  @inbounds for k in eachindex(bus_types)
    # Mirror solver-internal symbolic bus states back to Net node types.
    # Only PV/PQ are synchronized here; Slack/Isolated are managed in their
    # dedicated control paths to avoid accidental role changes in finalization.
    if bus_types[k] == :PQ
      setNodeType!(net.nodeVec[k], "PQ")
    elseif bus_types[k] == :PV
      setNodeType!(net.nodeVec[k], "PV")
    end
  end
  return nothing
end

function _finalize_rectangular_voltage_writeback!(net::Net, V::Vector{ComplexF64})
  # Single writeback point for solved complex voltages to keep downstream
  # post-processing aligned with the final NR state.
  update_net_voltages_from_complex!(net, V)
  return nothing
end

function _compute_rectangular_final_injections(Ybus, V::Vector{ComplexF64}, Sbase::Float64)
  # Compute injections in p.u. first, then scale once to physical units (MVA).
  # This keeps numerical reduction and reporting paths explicit and consistent.
  Sbus_pu = calc_injections(Ybus, V)
  Sbus_MVA = Sbus_pu .* Sbase
  return Sbus_pu, Sbus_MVA
end

function _write_rectangular_bus_power_results!(net::Net, nodes, Sbus_MVA)
  isoNodes = net.isoNodes
  @inbounds for (k, node) in enumerate(nodes)
    # Isolated nodes are intentionally skipped: no physically meaningful
    # network injection should be written back for disconnected buses.
    if k in isoNodes
      continue
    end

    Sbus = Sbus_MVA[k]
    Pbus_MW = real(Sbus)
    Qbus_MVar = imag(Sbus)

    if node._nodeType == Sparlectra.Slack
      # Slack absorbs the system residual; both P and Q are solver outcomes.
      node._pƩGen = Pbus_MW
      node._qƩGen = Qbus_MVar
    elseif node._nodeType == Sparlectra.PV
      # PV keeps scheduled active power; only reactive output is finalized here.
      node._qƩGen = Qbus_MVar
    elseif node._nodeType == Sparlectra.PQ
      if haskey(net.qLimitEvents, k)
        # Bus was converted PV->PQ due to reactive limit enforcement.
        # Keep previously clamped generator Q instead of overwriting it.
        @debug "Bus $(k) is PQ due to Q-limit; keeping clamped Qgen = $(node._qƩGen) MVar."
      end
    end
  end
  return nothing
end

function _write_rectangular_total_bus_power!(net::Net, Sbus_pu, Sbase::Float64, verbose::Int, performance_profile)
  p, q = _perf_profile_time!(performance_profile, :solver_total_power_reduction) do
    # Reduce from final per-bus injections (p.u.) and scale once to MW/MVAr.
    # This ensures totals are derived from exactly the same solved state as
    # per-bus writeback and avoids duplicate recomputation paths.
    ((sum(real.(Sbus_pu))) * Sbase, (sum(imag.(Sbus_pu))) * Sbase)
  end

  if verbose > 1
    @info "Set total bus power to p = $p MW and q = $q MVar"
  end

  _perf_profile_time!(performance_profile, :solver_total_power_writeback) do
    # Keep total power and shunt updates in one profiled block, because both
    # belong to the same post-solve consistency step.
    setTotalBusPower!(net = net, p = p, q = q)
    updateShuntPowers!(net = net)
  end

  return p, q
end

# Copyright 2023-2026 Udo Schmitz
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
# purpose: post-iteration finalization: bus-type sync back to the Net, solved
#          voltage write-back, final injection computation, and bus power
#          result write-back (P(U)/Q(U) units at the solved voltage included)

function _sync_rectangular_bus_types_to_net!(net::Net, bus_types)
  @inbounds for k in eachindex(bus_types)
    # Mirror solver-internal symbolic bus states back to Net node types.
    # Only PV/PQ are synchronized here; Slack/Isolated are managed in their
    # dedicated control paths to avoid accidental role changes in finalization.
    # An Isolated bus enters the solver as a neutral PQ row (zero injection)
    # and keeps its type: written back as PQ, the next solve on the same net
    # (a controller pass, a second runpf!) sees a bus without a branch and
    # without a source as an island with no reference and stops.
    getNodeType(net.nodeVec[k]) == Isolated && continue
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
  # The solve gives the bus injection V conj(Y V), generation minus what the
  # bus draws itself. The nodes keep the GENERATION (every reader prints the
  # bus load next to it, the Q-limit checks compare it with the generator
  # limits), so the bus's own consumption at the solved voltage goes back on
  # top. Before 0.41.0 the injection itself was written: a slack or PV bus
  # with a load showed its generation reduced by that load.
  p_cons, q_cons = _bus_consumption_MVA(net, length(nodes))
  @inbounds for (k, node) in enumerate(nodes)
    # Isolated nodes are intentionally skipped: no physically meaningful
    # network injection should be written back for disconnected buses.
    if k in isoNodes
      continue
    end

    Sbus = Sbus_MVA[k]
    Pgen_MW = real(Sbus) + p_cons[k]
    Qgen_MVar = imag(Sbus) + q_cons[k]

    if node._nodeType == Sparlectra.Slack
      # Slack absorbs the system residual; both P and Q are solver outcomes.
      node._pƩGen = Pgen_MW
      node._qƩGen = Qgen_MVar
    elseif node._nodeType == Sparlectra.PV
      # PV keeps scheduled active power; only reactive output is finalized here.
      node._qƩGen = Qgen_MVar
    elseif node._nodeType == Sparlectra.PQ
      if haskey(net.qLimitEvents, k)
        # Bus was converted PV->PQ due to reactive limit enforcement.
        # Keep previously clamped generator Q instead of overwriting it.
        @debug "Bus $(k) is PQ due to Q-limit; keeping clamped Qgen = $(node._qƩGen) MVar."
      end
    end
  end
  # a clamped bus with a Q(U) unit that follows its characteristic: the
  # clamp record holds that unit at its Q at the held voltage, the bus
  # reports its physical generator Q instead
  _clamped_bus_physical_qgen!(net, Sbus_MVA, q_cons)
  # P(U)/Q(U) units are written by _write_controlled_generation!, which the
  # solver calls only for an accepted solution
  return nothing
end

# Generation of the buses with a P(U) or Q(U) generator unit at the solved
# voltage magnitude (MW, MVAr): the sum over the generator units of the
# bus, a unit with a characteristic at that voltage, every other unit at
# its schedule pVal/qVal (the sum of _bus_scheduled_pgen and
# _bus_power_component_cache). `has_p`/`has_q` mark the buses with such a
# unit; `p`/`q` hold the sum on every bus with a mark of either kind, the
# other entries stay 0.0. Returns (p, q, has_p, has_q).
function _controlled_bus_generation(net::Net)
  nodes = net.nodeVec
  nb = length(nodes)
  has_p = falses(nb)
  has_q = falses(nb)
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    bus = getPosumerBusIndex(ps)
    1 <= bus <= nb || continue
    has_pu_controller(ps) && (has_p[bus] = true)
    has_qu_controller(ps) && (has_q[bus] = true)
  end
  p = zeros(Float64, nb)
  q = zeros(Float64, nb)
  (any(has_p) || any(has_q)) || return p, q, has_p, has_q
  base = net.baseMVA
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    bus = getPosumerBusIndex(ps)
    1 <= bus <= nb || continue
    (has_p[bus] || has_q[bus]) || continue
    vm = _safe_vm(Float64(something(nodes[bus]._vm_pu, 1.0)))
    p[bus] += has_pu_controller(ps) ? first(evaluate_controller(ps.puController, vm)) * base : something(ps.pVal, 0.0)
    q[bus] += has_qu_controller(ps) ? first(evaluate_controller(ps.quController, vm)) * base : something(ps.qVal, 0.0)
  end
  return p, q, has_p, has_q
end

# Write the generation of the buses with a P(U) or Q(U) generator unit,
# evaluated at the solved voltage magnitude, to `node._pƩGen` and
# `node._qƩGen` (MW, MVAr), as _controlled_bus_generation sums it. The
# solve takes such a unit's P and Q from its characteristic at the voltage
# of each iteration (buildControlledSVec), while the node kept the sum of
# the scheduled pVal/qVal from the import: sp_case118 with buses 4, 8 and
# 10 as PQ with a Q(U) characteristic (+30 MVAr at 0.95 pu to -30 MVAr at
# 1.05 pu) and a 15 MVAr schedule kept 15 MVAr on the node where the
# machines deliver -1.42, -5.67 and -10.78 MVAr.
#
# Written for P: every bus but the slack (slack P is the solve's own
# result, PV and PQ P is a schedule). For Q: a PQ bus without a Q-limit
# event. Slack and PV Q are the solve's own result, and a clamped bus keeps
# its clamp record or the physical Q of _clamped_bus_physical_qgen!, which
# the active set and the Q-V check read as the enforced value. The Q-V
# check is not affected: _qv_effective_qgen takes the characteristic, not
# the node, at a bus with a Q(U) unit, and qvCharacteristicViolations only
# looks at buses with a Q-limit event, which this function leaves alone.
# Isolated buses are skipped like in the write-back above. Called by
# runpf_rectangular! only for a converged and accepted solve: a diverged or
# rejected run keeps the schedule on the node. A run on a link-merged copy
# writes the representative here; _split_merged_bus_powers! hands the
# values of the cluster members to the original net. Returns the number of
# buses written.
function _write_controlled_generation!(net::Net)
  p, q, has_p, has_q = _controlled_bus_generation(net)
  (any(has_p) || any(has_q)) || return 0
  nodes = net.nodeVec
  iso = Set{Int}(net.isoNodes)
  written = 0
  for bus in eachindex(nodes)
    (has_p[bus] || has_q[bus]) || continue
    bus in iso && continue
    t = getNodeType(nodes[bus])
    (t == Slack || t == Isolated) && continue
    write_q = has_q[bus] && t == PQ && !haskey(net.qLimitEvents, bus)
    has_p[bus] && (nodes[bus]._pƩGen = p[bus])
    write_q && (nodes[bus]._qƩGen = q[bus])
    (has_p[bus] || write_q) && (written += 1)
  end
  return written
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

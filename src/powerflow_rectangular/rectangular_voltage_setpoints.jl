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
# Rectangular power-flow voltage-setpoint lookup helpers.
#
# This file is included inside module Sparlectra. Do not add a module wrapper here.

# Date: 29.5.2026
# file: src/powerflow_rectangular/rectangular_voltage_setpoints.jl

"""
Build per-bus voltage-magnitude setpoints for rectangular power-flow initialization.

Setpoints are resolved in explicit phases (mapping, regulator scan, fallback, fill)
so timing/profiling can attribute costs to each step and diagnostics can explain
where a bus setpoint came from.

Resolution precedence per bus:
1. Voltage-regulating prosumer setpoint (`ps.vm_pu`) if available.
   - Slack prosumer setpoints are considered only when `net.flatstart == true`.
2. Fallback to node voltage (`node._vm_pu`) from imported/current network state.
3. Final fallback `1.0` p.u. when node voltage is missing.

If multiple regulating prosumers define different setpoints for the same bus,
the first encountered value is kept and a debug message is emitted.

# Arguments
- `net::Net`: Network model containing buses and prosumers.

# Keywords
- `performance_profile=nothing`: Optional profiling dictionary used by
  `_perf_profile_time!` for phase timing and counters.

# Returns
- `Vector{Float64}`: Voltage setpoint vector in solver bus order.
"""
function _bus_voltage_setpoints_from_prosumers(net::Net; performance_profile = nothing)::Vector{Float64}
  # Keep resolution in explicit phases: easier to profile and to diagnose
  # whether setpoints come from regulators or fallback node metadata.
  nodes = net.nodeVec
  nbus = length(nodes)

  bus_idx_to_position, matpower_bus_i_to_position = _perf_profile_time!(performance_profile, :vset_map_build) do
    # Build index maps once to avoid repeated dictionary lookups/scans in later phases.
    # `bus_idx_to_position`: internal bus index -> solver position
    # `matpower_bus_i_to_position`: original MATPOWER BUS_I -> solver position
    bus_idx_to_position_ = Dict{Int,Int}()
    matpower_bus_i_to_position_ = Dict{Int,Int}()
    sizehint!(bus_idx_to_position_, nbus)
    sizehint!(matpower_bus_i_to_position_, nbus)
    @inbounds for k in eachindex(nodes)
      bus_idx_to_position_[k] = k
      matpower_bus_i_to_position_[get(net.busOrigIdxDict, k, k)] = k
    end
    (bus_idx_to_position_, matpower_bus_i_to_position_)
  end

  generator_bus_position_to_vg, has_online_voltage_regulator = _perf_profile_time!(performance_profile, :vset_generator_scan) do
    # Regulating setpoints override fallback node Vm because they encode the
    # operational control target used by the solver.
    generator_bus_position_to_vg_ = Dict{Int,Float64}()
    has_online_voltage_regulator_ = falses(nbus)
    sizehint!(generator_bus_position_to_vg_, min(length(net.prosumpsVec), nbus))
    use_slack_prosumer_setpoint = net.flatstart
    @inbounds for ps in net.prosumpsVec
      bus_idx = getPosumerBusIndex(ps)
      bus_pos = get(bus_idx_to_position, bus_idx, 0)
      bus_pos == 0 && continue

      # Slack setpoint is considered only for flat-start seeds.
      # Generator regulation is always considered if prosumer is regulating.
      ((use_slack_prosumer_setpoint && isSlack(ps)) || (isGenerator(ps) && isRegulating(ps))) || continue

      has_online_voltage_regulator_[bus_pos] = true
      isnothing(ps.vm_pu) && continue
      vm_pu = Float64(ps.vm_pu)

      if haskey(generator_bus_position_to_vg_, bus_pos)
        # Deterministic precedence: keep the first discovered setpoint to avoid
        # non-reproducible startup behavior from iteration/order side effects.
        kept = generator_bus_position_to_vg_[bus_pos]
        if abs(kept - vm_pu) > 1e-8
          @debug "Multiple voltage-regulating prosumers with different setpoints at bus $(bus_idx); keeping first setpoint." kept = kept ignored = vm_pu
        end
        continue
      end
      generator_bus_position_to_vg_[bus_pos] = vm_pu
    end
    (generator_bus_position_to_vg_, has_online_voltage_regulator_)
  end

  # Fallback captures imported/current per-bus voltage magnitude when no explicit
  # regulating setpoint exists; 1.0 pu is the safe neutral default.
  fallback_bus_vm = _perf_profile_time!(performance_profile, :vset_fallback_bus_vm) do
    fallback_bus_vm_ = Vector{Float64}(undef, nbus)
    @inbounds for k in eachindex(nodes)
      node = nodes[k]
      fallback_bus_vm_[k] = isnothing(node._vm_pu) ? 1.0 : Float64(node._vm_pu)
    end
    fallback_bus_vm_
  end

  _perf_profile_time!(performance_profile, :vset_node_metadata_scan) do
    # Keep mapping present in profile output for traceability between solver order
    # and original case numbering in logs/diagnostics.
    length(matpower_bus_i_to_position)
  end

  Vset = _perf_profile_time!(performance_profile, :vset_vector_fill) do
    # Two-stage fill avoids conditional branching in the solver:
    # initialize all buses with fallback, then override regulated buses.
    Vset_ = copy(fallback_bus_vm)
    @inbounds for (bus_pos, vm_pu) in generator_bus_position_to_vg
      Vset_[bus_pos] = vm_pu
    end
    Vset_
  end

  _perf_profile_time!(performance_profile, :vset_missing_online_gen_summary) do
    # Diagnostic counter: PV/Slack buses without active regulator setpoint.
    # Useful to explain when fallback Vm determined startup for controlled buses.
    missing = 0
    @inbounds for k in eachindex(nodes)
      nt = getNodeType(nodes[k])
      if (nt == PV || nt == Slack) && !has_online_voltage_regulator[k]
        missing += 1
      end
    end
    missing
  end

  return Vset
end

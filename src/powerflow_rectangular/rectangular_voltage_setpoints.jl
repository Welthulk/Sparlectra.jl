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

function _bus_voltage_setpoints_from_prosumers(net::Net; performance_profile = nothing)::Vector{Float64}
  nodes = net.nodeVec
  nbus = length(nodes)

  bus_idx_to_position, matpower_bus_i_to_position = _perf_profile_time!(performance_profile, :vset_map_build) do
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
    # Generator/slack regulator setpoints take precedence over fallback node Vm.
    generator_bus_position_to_vg_ = Dict{Int,Float64}()
    has_online_voltage_regulator_ = falses(nbus)
    sizehint!(generator_bus_position_to_vg_, min(length(net.prosumpsVec), nbus))
    use_slack_prosumer_setpoint = net.flatstart
    @inbounds for ps in net.prosumpsVec
      bus_idx = getPosumerBusIndex(ps)
      bus_pos = get(bus_idx_to_position, bus_idx, 0)
      bus_pos == 0 && continue
      ((use_slack_prosumer_setpoint && isSlack(ps)) || (isGenerator(ps) && isRegulating(ps))) || continue
      has_online_voltage_regulator_[bus_pos] = true
      isnothing(ps.vm_pu) && continue
      vm_pu = Float64(ps.vm_pu)
      if haskey(generator_bus_position_to_vg_, bus_pos)
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

  # Fallback is imported/current node voltage when no regulating setpoint is found.
  fallback_bus_vm = _perf_profile_time!(performance_profile, :vset_fallback_bus_vm) do
    fallback_bus_vm_ = Vector{Float64}(undef, nbus)
    @inbounds for k in eachindex(nodes)
      node = nodes[k]
      fallback_bus_vm_[k] = isnothing(node._vm_pu) ? 1.0 : Float64(node._vm_pu)
    end
    fallback_bus_vm_
  end

  _perf_profile_time!(performance_profile, :vset_node_metadata_scan) do
    # Keep MATPOWER original BUS_I to solver-position mapping explicit in profiles.
    # The setpoint fill below uses solver order, while diagnostics can still relate
    # a position to the imported MATPOWER bus number without rescanning nodes.
    length(matpower_bus_i_to_position)
  end

  Vset = _perf_profile_time!(performance_profile, :vset_vector_fill) do
    Vset_ = copy(fallback_bus_vm)
    @inbounds for (bus_pos, vm_pu) in generator_bus_position_to_vg
      Vset_[bus_pos] = vm_pu
    end
    Vset_
  end

  _perf_profile_time!(performance_profile, :vset_missing_online_gen_summary) do
    # Count PV/slack buses without online regulators for diagnostics/profile summaries.
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

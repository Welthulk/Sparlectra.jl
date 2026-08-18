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
# Rectangular power-flow Q-limit guard preprocessing helper.

# Date: 29.5.2026
# file: src/powerflow_rectangular/rectangular_qlimit_guard.jl
# purpose: Q-limit guard preprocessing for the rectangular active set:
#          reclassifies PV buses with zero or narrow reactive range as PQ
#          before the Newton iteration starts

"""
    _lock_zero_range_pv_buses!(net, bus_types, S, Qload_pu, qmin_pu, qmax_pu; log, verbose) -> Vector{Int}

Reclassify PV buses whose reactive limits coincide (`qmin == qmax`) as PQ,
with their reactive injection fixed at that value. This is modelling, not a
heuristic: without reactive headroom a machine cannot hold a voltage, so it
is not a voltage-controlled bus. Runs unconditionally whenever Q limits are
enforced — unlike the narrow-range guard, there is no threshold to choose
and no case where the PV treatment would be correct.

Returns the reclassified bus indices.
"""
function _lock_zero_range_pv_buses!(net::Net, bus_types::Vector{Symbol}, S::Vector{ComplexF64}, Qload_pu::Vector{Float64}, qmin_pu::AbstractVector, qmax_pu::AbstractVector; log::Bool = true, verbose::Int = 0)
  locked = Int[]
  @inbounds for bus in eachindex(bus_types)
    bus_types[bus] == :PV || continue
    bus <= length(qmin_pu) && bus <= length(qmax_pu) || continue
    qmin = qmin_pu[bus]
    qmax = qmax_pu[bus]
    (isfinite(qmin) && isfinite(qmax)) || continue
    qmax - qmin == 0.0 || continue
    bus_types[bus] = :PQ
    S[bus] = ComplexF64(real(S[bus]), qmin - Qload_pu[bus])
    net.nodeVec[bus]._qƩGen = qmin * net.baseMVA
    logQLimitHit!(net, 0, bus, qmin >= 0.0 ? :max : :min)
    push!(locked, bus)
  end
  if log && verbose > 0 && !isempty(locked)
    @printf("Q-limit: %d PV bus(es) have no reactive headroom (qmin == qmax) and start as PQ.\n", length(locked))
  end
  return locked
end

function _apply_qlimit_guard_to_rectangular_active_set!(net::Net, bus_types::Vector{Symbol}, S::Vector{ComplexF64}, Qload_pu::Vector{Float64}, qmin_pu::AbstractVector, qmax_pu::AbstractVector; min_q_range_pu::Float64, zero_range_mode::Symbol, narrow_range_mode::Symbol, log::Bool, verbose::Int)
  # Guard runs once before NR to lock obviously problematic narrow-Q PV buses.
  min_q_range_pu >= 0.0 || error("qlimit_guard_min_q_range_pu must be >= 0 (got $(min_q_range_pu)).")
  zero_range_mode in (:lock_pq, :prefer_pq, :delayed_switch, :ignore) || error("Unsupported qlimit_guard_zero_range_mode=$(zero_range_mode). Supported: :lock_pq, :prefer_pq, :delayed_switch, :ignore.")
  narrow_range_mode in (:lock_pq, :prefer_pq, :delayed_switch, :ignore) || error("Unsupported qlimit_guard_narrow_range_mode=$(narrow_range_mode). Supported: :lock_pq, :prefer_pq, :delayed_switch, :ignore.")

  # Preprocess suspiciously narrow finite Q-ranges before active-set iteration.
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

  # The remaining network-integrated solver loop still performs PV→PQ switching.
  if log && verbose > 0 && !isempty(guarded)
    @printf("Q-limit guard: locked %d narrow-range PV bus(es) as PQ before rectangular NR.\n", length(guarded))
  end
  return guarded
end

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
#          before the Newton iteration starts and gives them a start
#          voltage consistent with their fixed reactive injection

"""
    _start_guarded_buses_from_fixed_injection!(V0, Ybus, S, buses; max_sweeps = 10, tol = 1e-10) -> Int

Start voltage of the buses the Q-limit guard ran as PQ (`buses`): each
one gets the voltage its own bus equation gives with the fixed injection
`S[bus]` against its neighbours' start voltages, instead of the setpoint
of a machine that cannot hold it. A machine whose reactive range is zero
or narrower than `min_q_range_pu` is a fixed-Q injection, so its setpoint
says nothing about the voltage the bus will settle at; starting the bus
there leaves a reactive mismatch of the size the machine would need to
hold the setpoint (measured on case_SyntheticUSA island 1: 145 pu on one
zero-range 69 kV bus whose setpoint is 1.04 pu and whose bus voltage in
the case file is 1.009 pu), and Newton did not converge from that start.

The equations of the guarded buses are solved together by Jacobi sweeps
with every other bus held at its start value: per sweep,
`V[b] = (conj(S[b]) / conj(V[b]) - sum(Y[b,j] * V[j], j != b)) / Y[b,b]`
for each guarded bus, until the largest voltage change is below `tol` or
`max_sweeps` sweeps have run. Only `V0` entries of `buses` change. A bus
with a zero diagonal or a non-finite update keeps its previous value.

Returns the number of sweeps run (0 when `buses` is empty).
"""
function _start_guarded_buses_from_fixed_injection!(V0::AbstractVector{ComplexF64}, Ybus::AbstractMatrix, S::AbstractVector{ComplexF64}, buses::Vector{Int}; max_sweeps::Int = 10, tol::Float64 = 1e-10)
  isempty(buses) && return 0
  Vprev = Vector{ComplexF64}(undef, length(buses))
  sweeps = 0
  for _ = 1:max_sweeps
    sweeps += 1
    # Jacobi: the network current of this sweep is taken from the voltages
    # of the previous sweep for every bus, guarded neighbours included
    I = Ybus * V0
    @inbounds for (k, b) in enumerate(buses)
      Vprev[k] = V0[b]
    end
    max_change = 0.0
    @inbounds for (k, b) in enumerate(buses)
      ybb = Ybus[b, b]
      ybb == 0 && continue
      others = I[b] - ybb * Vprev[k]
      Vn = (conj(S[b]) / conj(Vprev[k]) - others) / ybb
      (isfinite(real(Vn)) && isfinite(imag(Vn)) && abs(Vn) > 0.0) || continue
      V0[b] = Vn
      max_change = max(max_change, abs(Vn - Vprev[k]))
    end
    max_change < tol && break
  end
  return sweeps
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

"""
    _lock_listed_pv_buses_as_pq!(net, bus_types, S, qmin_pu, qmax_pu, buses; verbose) -> Vector{Int}

Run the PV buses listed in `power_flow.qlimits.lock_pv_to_pq_buses` (internal
bus positions) as PQ from the start of the solve. Each listed bus that is PV
in `bus_types` becomes PQ with its scheduled reactive generation (the sum of
the generator `qVal` at the bus) clamped into the bus's `[qmin, qmax]`, so the
run cannot end with a reactive violation at a bus it was told not to regulate.
Listed buses that are not PV (slack, PQ, already reclassified by the guard)
and positions outside the network are left alone.

The conversion is logged at iteration 0 like every other setup-time
reclassification (zero-range and narrow-range guard, warm active set): the
side is the limit the scheduled value was clamped to, otherwise the sign of
the run value (the guard's convention). The caller runs this BEFORE the PV
origin mask is taken, so the PQ->PV release rule never releases a listed
bus: the list is a lock, not a start state.

Returns the converted bus positions.
"""
function _lock_listed_pv_buses_as_pq!(net::Net, bus_types::Vector{Symbol}, S::Vector{ComplexF64}, qmin_pu::AbstractVector, qmax_pu::AbstractVector, buses::AbstractVector{Int}; verbose::Int = 0)
  locked = Int[]
  isempty(buses) && return locked
  nb = length(bus_types)
  # scheduled reactive generation per bus (pu), the same sum buildComplexSVec
  # puts into S; only the generator part is replaced below, so loads and any
  # voltage-dependent shunt injection already in S stay untouched
  qgen_sched = zeros(Float64, nb)
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    bus = getPosumerBusIndex(ps)
    (1 <= bus <= nb) || continue
    qgen_sched[bus] += something(ps.qVal, 0.0) / net.baseMVA
  end
  for bus in unique(buses)
    (1 <= bus <= nb && bus_types[bus] == :PV) || continue
    qmin = bus <= length(qmin_pu) ? Float64(qmin_pu[bus]) : -Inf
    qmax = bus <= length(qmax_pu) ? Float64(qmax_pu[bus]) : Inf
    q = qgen_sched[bus]
    side = q >= 0.0 ? :max : :min
    if isfinite(qmax) && q > qmax
      q, side = qmax, :max
    elseif isfinite(qmin) && q < qmin
      q, side = qmin, :min
    end
    bus_types[bus] = :PQ
    S[bus] = ComplexF64(real(S[bus]), imag(S[bus]) - qgen_sched[bus] + q)
    net.nodeVec[bus]._qƩGen = q * net.baseMVA
    logQLimitHit!(net, 0, bus, side)
    push!(locked, bus)
  end
  if verbose > 0 && !isempty(locked)
    @printf("Q-limit: %d PV bus(es) of power_flow.qlimits.lock_pv_to_pq_buses run as PQ from the start (bus %s).\n", length(locked), join((_qlimit_original_bus_id(net, b) for b in locked), ", "))
  end
  return locked
end

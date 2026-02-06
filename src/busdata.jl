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


# file: src/busdata.jl — Data type for NR / power flow

mutable struct BusData
  idx::Int          # # PF index (non-iso compact), Bus index (after sorting)
  nodeIdx::Int      # index in net.nodeVec (original)
  vm_pu::Float64    # Voltage in p.u.
  va_rad::Float64   # Angle in rad
  pƩ::Float64       # Sum active power (p.u.)
  qƩ::Float64       # Sum reactive power (p.u.)
  _pRes::Float64    # calculated active power (p.u.)
  _qRes::Float64    # calculated reactive power (p.u.)
  type::NodeType
  isPVactive::Union{Nothing,Bool}

  function BusData(idx::Int, nodeIdx::Int, vm_pu::Float64, va_rad::Float64, sumP::Float64, sumQ::Float64, type::NodeType)
    if type == PV
      isPV = true
    else
      isPV = nothing
    end

    new(idx, nodeIdx, vm_pu, va_rad, sumP, sumQ, 0.0, 0.0, type, isPV)
  end
end

function Base.show(io::IO, bus::BusData)
  va_deg = round(rad2deg(bus.va_rad), digits = 3)
  print(io, "BusData($(bus.idx), $(bus.vm_pu), $(va_deg)°), $(bus.pƩ), $(bus.qƩ), $(bus._pRes), $(bus._qRes), $(bus.type))")
end

function getBusData(nodes::Vector{Node}, Sbase_MVA::Float64, flatStart)
  busVec = Vector{BusData}()

  slackIdx = 0

  sumLoad_p = 0.0
  sumLoad_q = 0.0
  sumGen_p = 0.0
  sumGen_q = 0.0
  idx = 0
  for (i, n) in enumerate(nodes)
    isIsolated(n) && continue
    idx += 1
    nodeIdx = i
    if isSlack(n)
      slackIdx = idx
    end

    type = n._nodeType

    p = 0
    q = 0
    p += n._pƩLoad === nothing ? 0.0 : n._pƩLoad * -1.0
    q += n._qƩLoad === nothing ? 0.0 : n._qƩLoad * -1.0

    sumLoad_p += n._pƩLoad === nothing ? 0.0 : n._pƩLoad * -1.0
    sumLoad_q += n._qƩLoad === nothing ? 0.0 : n._qƩLoad * -1.0
    # Note: Shunts are considered in Y-Bus Matrix

    p += n._pƩGen === nothing ? 0.0 : n._pƩGen
    q += n._qƩGen === nothing ? 0.0 : n._qƩGen

    sumGen_p += n._pƩGen === nothing ? 0.0 : n._pƩGen
    sumGen_q += n._qƩGen === nothing ? 0.0 : n._qƩGen

    p = p / Sbase_MVA
    q = q / Sbase_MVA

    va_rad = 0.0
    if flatStart
      if type == PQ
        vm_pu = 1.0
        va_rad = 0.0
      elseif type == PV
        vm_pu = n._vm_pu === nothing ? 1.0 : n._vm_pu
        va_rad = 0.0
      elseif type == Slack
        vm_pu = n._vm_pu === nothing ? 1.0 : n._vm_pu
        #va_rad = n._va_deg === nothing ? 0.0 : deg2rad(n._va_deg)
        va_rad = 0.0
      end
    else
      vm_pu = n._vm_pu === nothing ? 1.0 : n._vm_pu
      va_rad = n._va_deg === nothing ? 0.0 : deg2rad(n._va_deg)

      if angle_limit && type == PQ
        if abs(n._va_deg) > deg2rad(30.0)
          va_rad = sign(n._va_deg) * deg2rad(30.0)
          @warn "getBusData: bus $(idx) type=$(type), va_deg=$(n._va_deg) exceeds angle limit, set to 30°"
        end
      end
    end
    busIdx = idx
    b = BusData(busIdx, nodeIdx, vm_pu, va_rad, p, q, type)    
    push!(busVec, b)
  end

  if slackIdx == 0
    throw("No slack node found")
  end

  sort!(busVec, by = x -> x.idx)

  sumLoad_p = sumLoad_p * -1.0
  sumLoad_q = sumLoad_q * -1.0
  delta_p = round((sumGen_p - sumLoad_p), digits = 3)
  delta_q = round((sumGen_q - sumLoad_q), digits = 3)

  if debug
    println("\n∑Load: [$(sumLoad_p), $(sumLoad_q)], ∑Gen [$(sumGen_p), $(sumGen_q)]  Δp, Δq: [$(delta_p), $(delta_q)]")
  end

  if debug
    println("\nslack bus: $(slackIdx)")
    for b in busVec
      println("$(b)")
    end
  end

  return busVec, slackIdx
end # getBusData

# helper function to count number of nodes of type = value [PQ, PV]
function getBusTypeVec(busVec::Vector{BusData})
  busTypeVec = Vector{NodeType}()
  slackIdx = 0
  idx = 0
  for (i, bus) in enumerate(busVec) #1:n    
    idx += 1
    if bus.type == Slack
      slackIdx = idx
    end
    push!(busTypeVec, bus.type)
  end
  if debug
    for (i, bus) in enumerate(busTypeVec)
      print("busVec[$(i)]: $(bus), ")
    end
  end
  return busTypeVec, slackIdx
end

# count number of nodes of type = value [PQ, PV]
function countNodes(busTypeVec::Vector{NodeType}, pos, value::NodeType)
  sum = 0

  for (index, bus_type) in enumerate(busTypeVec)
    if index >= pos
      break
    end

    if bus_type == value
      sum += 1
    end
  end

  return sum
end

"""
    map_NR_voltage_to_net!(V_nr, busVec, net) -> V_net

Maps the NR-solver voltage vector `V_nr` (in busVec order)
into the original bus index order used by `net.nodeVec`.

Returns a Vector{ComplexF64} such that:
    V_net[busIdx] == V_nr[k]
where `k` is the index in busVec with `busVec[k].idx == busIdx`.
"""
function map_NR_voltage_to_net!(V_nr::Vector{ComplexF64}, busVec::Vector{BusData}, net::Net)
  V_net = Vector{ComplexF64}(undef, length(net.nodeVec))

  @inbounds for k in eachindex(busVec)
    V_net[busVec[k].nodeIdx] = V_nr[k]
  end

  return V_net
end

"""
    buildVoltageVector_from_busVec(busVec::Vector{BusData}) -> Vector{ComplexF64}

Builds the complex voltage vector in NR ordering:
    V_nr[k] = busVec[k].vm_pu * exp(j * busVec[k].va_rad)

This is the canonical source of V_nr for mapping and post-processing.
"""
function buildVoltageVector_from_busVec(busVec::Vector{BusData})
  V_nr = Vector{ComplexF64}(undef, length(busVec))
  @inbounds for k in eachindex(busVec)
    bus = busVec[k]
    V_nr[k] = bus.vm_pu * exp(im * bus.va_rad)
  end
  return V_nr
end

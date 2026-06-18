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
# Rectangular power-flow Q-limit bus-ID and trace helpers.

# Preserve original MATPOWER BUS_I identifiers for diagnostics output.

# Date: 29.5.2026
# file: src/powerflow_rectangular/rectangular_qlimit_trace.jl

_qlimit_original_bus_id(net::Net, bus::Int)::Int = haskey(net.busOrigIdxDict, bus) ? net.busOrigIdxDict[bus] : bus

function _resolve_qlimit_trace_buses(net::Net, requested::AbstractVector{Int})::Vector{Int}
  isempty(requested) && return Int[]
  orig_to_net = Dict{Int,Int}()
  sizehint!(orig_to_net, length(net.busOrigIdxDict))
  for (net_idx, orig_idx) in net.busOrigIdxDict
    orig_to_net[orig_idx] = net_idx
  end
  resolved = Int[]
  # Trace requests may be provided in original numbering and are mapped to internal bus positions.
  for bus in requested
    if haskey(orig_to_net, bus)
      push!(resolved, orig_to_net[bus])
    elseif 1 <= bus <= length(net.nodeVec)
      push!(resolved, bus)
    end
  end
  unique!(resolved)
  # Stable sorting keeps logs deterministic across runs.
  sort!(resolved)
  return resolved
end

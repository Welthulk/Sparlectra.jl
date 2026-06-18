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
# Rectangular power-flow result write-back helpers.
#
# This file is included inside module Sparlectra. Do not add a module wrapper here.

# Date: 29.5.2026
# file: src/powerflow_rectangular/rectangular_result_updates.jl

"""
    update_net_voltages_from_complex!(net, V)

Update the bus voltage magnitudes and angles in the network from the
final complex voltages V (in per-unit).
"""
function update_net_voltages_from_complex!(net::Net, V::Vector{ComplexF64})
  nodes = net.nodeVec
  n = length(nodes)
  @assert length(V) == n

  # Node order is solver order; no bus-index remapping is applied here.
  for (k, node) in enumerate(nodes)
    # Write-back uses solver-order complex bus voltages directly.
    Vk = V[k]
    # Persist magnitude/angle view expected by downstream reporting code.
    vm = abs(Vk)
    va_rad = angle(Vk)
    va_deg = rad2deg(va_rad)
    node._vm_pu = vm
    node._va_deg = va_deg
  end
end

# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0.
#
# Rectangular power-flow result write-back helpers.
#
# This file is included inside module Sparlectra. Do not add a module wrapper here.

"""
    update_net_voltages_from_complex!(net, V)

Update the bus voltage magnitudes and angles in the network from the
final complex voltages V (in per-unit).
"""
function update_net_voltages_from_complex!(net::Net, V::Vector{ComplexF64})
  nodes = net.nodeVec
  n = length(nodes)
  @assert length(V) == n

  for (k, node) in enumerate(nodes)
    Vk = V[k]
    vm = abs(Vk)
    va_rad = angle(Vk)
    va_deg = rad2deg(va_rad)
    node._vm_pu = vm
    node._va_deg = va_deg
  end
end

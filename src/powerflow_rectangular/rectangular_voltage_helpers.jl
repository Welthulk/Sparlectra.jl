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
# Shared rectangular power-flow voltage helpers.
#
# This file is included inside module Sparlectra. Do not add a module wrapper here.

# Date: 11.6.2026
# file: src/powerflow_rectangular/rectangular_voltage_helpers.jl

function _apply_voltage_magnitude_preserving_angle(V::ComplexF64, vm::Float64)::ComplexF64
  if iszero(V)
    # A zero phasor has no defined angle, so use the explicit zero-angle fallback.
    return ComplexF64(vm, 0.0)
  end
  return ComplexF64(vm * cis(angle(V)))
end


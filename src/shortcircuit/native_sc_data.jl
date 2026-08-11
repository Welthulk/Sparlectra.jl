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
# This file is included inside module Sparlectra. Do not add a module wrapper
# here (same convention as shortcircuit/short_circuit.jl).
#
# file: src/shortcircuit/native_sc_data.jl
# purpose: native short-circuit source container for networks built without a
#          CGMES delivery (issue #299). Field-identical to
#          `CGMESImporter.CGMESShortCircuitData` — the short-circuit engine is
#          duck-typed over that record contract and consumes either container
#          unchanged. Included BEFORE network.jl because `Net` carries a typed
#          `sc_sources::NativeShortCircuitData` field.

"""
    NativeShortCircuitData

Short-circuit source records of a natively built network (`addExternalGrid!`),
consumed by [`runShortCircuit!`](@ref) exactly like the CGMES harvest.

Deliberately **field-identical** to `CGMESImporter.CGMESShortCircuitData` and
**not** related to it through a common supertype: the short-circuit engine is
duck-typed over the record contract (issue #299 task decision) — the six
`Vector{NamedTuple}` fields with their tuple keys ARE the interface, and a
type hierarchy would suggest a coupling that does not exist. Keep both
structs in sync field-by-field when the contract evolves.

All-empty by default; `addExternalGrid!` fills
`external_network_injections`. The remaining vectors exist so the engine's
machine/equivalent loops iterate empty collections instead of needing
per-container special cases.
"""
struct NativeShortCircuitData
  external_network_injections::Vector{NamedTuple}
  synchronous_machines::Vector{NamedTuple}
  ac_line_segments::Vector{NamedTuple}
  transformer_ends::Vector{NamedTuple}
  equivalent_injections::Vector{NamedTuple}
  asynchronous_machines::Vector{NamedTuple}
end

NativeShortCircuitData() = NativeShortCircuitData(NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[])

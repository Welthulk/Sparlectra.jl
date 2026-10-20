# Copyright 2023-2026 Udo Schmitz
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

# file: src/sv_compare_core.jl
# purpose: the parts of a comparison against an SV state that do not depend
#          on the CGMES importer: the model flow at one branch terminal and
#          the angle alignment per island. compareWithSV uses both, and so
#          does every other comparison against an SV state, so the numbers
#          agree by construction.

"""
    _sv_branch_terminal_flow(br, V, side, baseMVA) -> ComplexF64

Complex power in MW/MVAr flowing from the bus into branch `br` at its
terminal `side` (`:from`, anything else is the to side), from the complex
bus voltages `V` in pu indexed by bus number. The model is the solver's:
series admittance, the shunt arm of each terminal (the from arm behind the
tap), and the off-nominal ratio `tap_ratio` with `phase_shift_deg` on the
from side; a branch with `ratio == 0` has no tap. `calcNetLosses!` uses
the same ratio and shift.
"""
function _sv_branch_terminal_flow(br::Branch, V::AbstractVector{ComplexF64}, side::Symbol, baseMVA::Float64)::ComplexF64
  ys = inv(br.r_pu + im * br.x_pu)
  # per-terminal shunt arms (0.20.0), the from arm behind the tap
  ysh_from = _branch_y0_from(br)
  ysh_to = _branch_y0_to(br)
  tr = br.ratio == 0.0 ? 1.0 + 0im : br.tap_ratio * cis(deg2rad(br.phase_shift_deg))
  Vf, Vt = V[br.fromBus], V[br.toBus]
  return side == :from ? Vf * conj(((ys + ysh_from) / abs2(tr)) * Vf - (ys / conj(tr)) * Vt) * baseMVA : Vt * conj((ys + ysh_to) * Vt - (ys / tr) * Vf) * baseMVA
end

# Median of the values; 0.0 for none. Median rather than mean: one bus
# whose SV entry is off does not move the offset of its island.
function _sv_median(v::Vector{Float64})::Float64
  sorted = sort(v)
  n = length(sorted)
  n == 0 && return 0.0
  return isodd(n) ? sorted[(n + 1) ÷ 2] : (sorted[n ÷ 2] + sorted[n ÷ 2 + 1]) / 2.0
end

"""
    _sv_align_angles(dva, island_of) -> (aligned = Vector{Float64}, offset_deg = Float64, offsets = Dict{Int,Float64})

Remove one angle offset per island from the angle differences `dva`
(degrees, model minus SV, one entry per compared bus). `island_of[i]` is
the island of entry `i` (0 for none); the caller builds it from the
islands the power flow solves (`_power_flow_island_of_bus`: closed links
connect, also through a bus the branch-only marking isolates), never from
the branch-only view, where a part that hangs on the rest through a
closed switch alone becomes an island of its own and its median hides its
own deviation.

The offset of an island is the median of its finite differences, and
`aligned[i]` is `dva[i]` minus the offset of its island. `offset_deg` is
the offset of the island with the most finite differences, the reference
offset of the delivery as a whole; an entry whose island has no finite
difference is aligned with it. `offsets` maps every island with a finite
difference to its offset, so a caller can align further values of the
same islands (`get(offsets, island, offset_deg)`) with exactly the offsets
the entries got.
"""
function _sv_align_angles(dva::Vector{Float64}, island_of::Vector{Int})
  length(dva) == length(island_of) || throw(DimensionMismatch("_sv_align_angles: $(length(dva)) angle differences but $(length(island_of)) island numbers"))
  per_island = Dict{Int,Vector{Float64}}()
  for (d, island) in zip(dva, island_of)
    isfinite(d) || continue
    push!(get!(per_island, island, Float64[]), d)
  end
  offset_of = Dict(id => _sv_median(v) for (id, v) in per_island)
  main_island = isempty(per_island) ? 0 : argmax(Dict(id => length(v) for (id, v) in per_island))
  offset_deg = get(offset_of, main_island, 0.0)
  aligned = Float64[d - get(offset_of, island, offset_deg) for (d, island) in zip(dva, island_of)]
  return (aligned = aligned, offset_deg = offset_deg, offsets = offset_of)
end

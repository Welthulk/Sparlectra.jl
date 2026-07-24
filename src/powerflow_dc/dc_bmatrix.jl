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
# DC power-flow B' matrix assembly.
#
# This file is included inside module Sparlectra. Do not add a module wrapper here.

# file: src/powerflow_dc/dc_bmatrix.jl

"""
    _DcBranchTerm

One in-service branch's contribution to the DC power-flow matrices:
`branch_idx`/`from`/`to` identify it, `b` is its DC series susceptance
(`1/(x_pu*tap)`), and `shift_rad` is its phase-shift angle in radians (0 for
a plain line or an off-nominal-tap-only transformer). `Pf = b*(theta[from] -
theta[to]) - b*shift_rad` reproduces MATPOWER's `Bf*theta + Pfinj` per branch
without assembling a separate sparse `Bf` matrix.
"""
struct _DcBranchTerm
  branch_idx::Int
  from::Int
  to::Int
  b::Float64
  shift_rad::Float64
end

"""
    assemble_dc_bbus(net::Net) -> (Bbus, Pbusinj, terms::Vector{_DcBranchTerm})

DC power-flow susceptance matrix and phase-shift injection vector, ported
from MATPOWER's `makeBdc.m` (BSD-3-Clause; algorithm only, no code copied).
Series reactance only — resistance, shunt conductance/susceptance, and line
charging are all ignored, per the lossless DC model. Tap magnitude divides
the series susceptance; phase shift enters only as an injection term, never
as part of `Bbus` itself (MATPOWER's DC convention).

Reads `Branch.tap_ratio`/`Branch.phase_shift_deg` — the branch's *live*
effective tap this network was last stamped with — gated by
`Branch.ratio != 0.0` marking a transformer (mirrors `calcBranchRatio`'s
"line vs. transformer" convention in `src/branch.jl`). Only in-service
(`status == 1`) branches contribute.

Returns the full `n×n` sparse `Bbus` (before any slack elimination — reduce
with `reduce_susceptance_to_nonslack` using `value_of = x -> -x`,
since that helper expects positive off-diagonal branch weights and `Bbus`'s
off-diagonals are the standard-convention `-b`; see [`solve_dc_powerflow`](@ref)'s
docstring for why), the `n`-length bus-level shift injection `Pbusinj`
(per-unit), and `terms`, the per-branch data needed to compute branch flows
after `theta` is solved.
"""
function assemble_dc_bbus(net::Net)
  n = length(net.nodeVec)
  Pbusinj = zeros(Float64, n)
  terms = _DcBranchTerm[]
  sizehint!(terms, length(net.branchVec))
  I = Int[]
  J = Int[]
  V = Float64[]
  @inbounds for br in net.branchVec
    br.status == 1 || continue
    f = br.fromBus
    t = br.toBus
    is_transformer = br.ratio != 0.0
    tap = is_transformer && br.tap_ratio != 0.0 ? br.tap_ratio : 1.0
    shift_rad = is_transformer ? deg2rad(br.phase_shift_deg) : 0.0
    x = br.x_pu
    x == 0.0 && throw(ArgumentError("DC power flow requires a nonzero series reactance; branch $(br.branchIdx) ($(hasproperty(br.comp, :cName) ? br.comp.cName : "unnamed")) has x_pu == 0."))
    b = 1.0 / (x * tap)
    push!(I, f); push!(J, f); push!(V, b)
    push!(I, t); push!(J, t); push!(V, b)
    push!(I, f); push!(J, t); push!(V, -b)
    push!(I, t); push!(J, f); push!(V, -b)
    pfinj = -b * shift_rad
    Pbusinj[f] += pfinj
    Pbusinj[t] -= pfinj
    push!(terms, _DcBranchTerm(br.branchIdx, f, t, b, shift_rad))
  end
  Bbus = sparse(I, J, V, n, n)
  return Bbus, Pbusinj, terms
end

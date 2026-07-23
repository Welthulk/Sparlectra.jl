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
# DC power-flow linear solve and branch-flow computation.
#
# This file is included inside module Sparlectra. Do not add a module wrapper here.

# file: src/powerflow_dc/dc_solve.jl

_dc_bus_injection_pu(node::Node, baseMVA::Float64)::Float64 = (_bus_power_value(node._pƩGen) - _bus_power_value(node._pƩLoad)) / baseMVA

"""
    solve_dc_powerflow(net::Net, slack_idx::Int; angle_reference_rad=0.0, performance_profile=nothing)
        -> (theta::Vector{Float64}, Pf::Vector{Float64}, Pt::Vector{Float64}, slack_p_mw::Float64, terms::Vector{_DcBranchTerm})

Solve `Bbus * theta = Pbus - Pbusinj` for bus angles (radians), with the
slack bus fixed. `Pbus[k] = (node._pƩGen - node._pƩLoad)/baseMVA` — MATPOWER's
`Pbus`, built from specified generation/load, independent of any AC solve
(shunts are deliberately excluded: the lossless DC model has no shunt term).

The reduced linear system is built assuming a zero slack angle (matching
MATPOWER's `dcpf.m` convention and reusing `reduce_susceptance_to_nonslack`
and `_solve_dc_angle_system` unmodified, shared with the rectangular
solver's DC-angle start projection); `angle_reference_rad` is then added as a
uniform shift to every bus angle afterward. This is exact, not approximate:
`Bbus`'s rows sum to zero by construction (a weighted graph Laplacian), so
`Bbus*(theta .+ c) == Bbus*theta` for any constant `c` — fixing the slack
angle to a nonzero reference is mathematically equivalent to solving with a
zero reference and shifting the whole result afterward.

Note on `reduce_susceptance_to_nonslack`'s `value_of`: the helper
reconstructs the reduced Laplacian purely from `Bsrc`'s *off-diagonal*
entries (its own diagonal is ignored), and expects those off-diagonal
entries to already be the *positive* per-branch weight — exactly what
`imag(Ybus)` off-diagonals are for the legacy AC-start heuristic. `Bbus`
(from [`assemble_dc_bbus`](@ref)) is a real MATPOWER-style Laplacian with
*negative* off-diagonals (`-b`, matching `makeBdc`'s actual `Bbus`, which
must stay that way since it's a documented public return value), so
`value_of = x -> -x` is passed here to flip it back to the positive
convention the helper expects — `value_of = identity` would silently
double-negate every angle and flow.

`Pf[l] = b_l*(theta[from_l] - theta[to_l]) - b_l*shift_rad_l`, in MW after
scaling by `baseMVA`; `Pt = -Pf` (lossless, not an independent unknown).
`slack_p_mw` is the slack bus's injection implied by the solved flows (its
net outflow across incident branches), not read from `Pbus` at the slack bus
(which is not part of the reduced solve).
"""
function solve_dc_powerflow(net::Net, slack_idx::Int; angle_reference_rad::Float64 = 0.0, performance_profile = nothing)
  n = length(net.nodeVec)
  Bbus, Pbusinj, terms = assemble_dc_bbus(net)

  non_slack = non_slack_indices(n, slack_idx)
  nred = length(non_slack)
  theta = zeros(Float64, n)
  if nred > 0
    pos = build_pos_map(non_slack, n)
    Bred = reduce_susceptance_to_nonslack(Bbus, non_slack, pos, slack_idx, nred; value_of = x -> -x)
    Pbus = [_dc_bus_injection_pu(net.nodeVec[i], net.baseMVA) for i in non_slack]
    P = Pbus .- Pbusinj[non_slack]
    θred = _solve_dc_angle_system(Bred, P, performance_profile, nred)
    @inbounds for (k, i) in enumerate(non_slack)
      theta[i] = θred[k]
    end
  end
  if angle_reference_rad != 0.0
    theta .+= angle_reference_rad
  end

  nbr = length(terms)
  Pf = Vector{Float64}(undef, nbr)
  Pt = Vector{Float64}(undef, nbr)
  net_outflow_pu = zeros(Float64, n)
  @inbounds for (l, term) in enumerate(terms)
    pf_pu = term.b * (theta[term.from] - theta[term.to]) - term.b * term.shift_rad
    Pf[l] = pf_pu * net.baseMVA
    Pt[l] = -Pf[l]
    net_outflow_pu[term.from] += pf_pu
    net_outflow_pu[term.to] -= pf_pu
  end
  slack_p_mw = net_outflow_pu[slack_idx] * net.baseMVA

  return theta, Pf, Pt, slack_p_mw, terms
end

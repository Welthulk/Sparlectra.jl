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

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 04.09.2023
# file: src/losses.jl
# purpose: post-solve branch-flow and network-loss calculation (calcNetLosses!)
#          and KCL-based flow allocation for bus links (calcLinkFlowsKCL!)

"""
    calcNetLosses!(net::Net)

Calculates branch flows and network losses for the given network.

This default method builds the complex voltage vector internally and forwards
to `calcNetLosses!(net, V)`. If the NR solver already has V available, it can
call the two-argument variant directly to avoid recomputing V.
"""
function calcNetLosses!(net::Net)
  V = buildVoltageVector(net)
  calcNetLosses!(net, V)
end

"""
    _closed_branch_flow_pu(V, from, to, br, tapSide) -> ComplexF64

Complex branch power S_ij in per unit from bus `from` to bus `to` of a
fully CLOSED branch, evaluated on the complex voltage vector `V`. This is
the flow formula `calcNetLosses!` uses (extracted so the contingency
screening can estimate loadings from a trial voltage vector without
writing the net); `tapSide` says which terminal carries the tap (1 =
`from`, 2 = `to`).
"""
function _closed_branch_flow_pu(V::Vector{ComplexF64}, from::Int, to::Int, br::Branch, tapSide::Int)
  @assert tapSide == 1 || tapSide == 2
  ui = V[from]
  uj = V[to]
  # the LIVE tap (tap_ratio, phase_shift_deg), the complex ratio the Ybus
  # stamps. The neutral ratio and angle differ from it where a nameplate tap
  # sits off its neutral position (an SCF tap changer at pos 4): the flows
  # then missed the solved state, 5.9 MVA of KCL error at the transformer
  # bus on sp_case14.
  tap = calcBranchRatio(br)
  if tapSide == 1
    ui /= tap
  elseif tapSide == 2
    uj /= tap
  end
  Yik = inv(br.r_pu + im * br.x_pu)
  # the shunt arm of the end the flow leaves (0.20.0): tapSide 1 is the from
  # end (arm behind the tap), tapSide 2 the to end
  Y0ik = tapSide == 1 ? _branch_y0_from(br) : _branch_y0_to(br)
  return abs(ui)^2 * conj(Y0ik + Yik) - ui * conj(uj) * conj(Yik)
end

"""
    calcNetLosses!(net::Net, V::Vector{ComplexF64})

Calculates branch flows and network losses using an externally provided complex
voltage vector `V` (typically from the final NR residual).
"""
function calcNetLosses!(net::Net, V::Vector{ComplexF64})
  nodes     = net.nodeVec
  branchVec = net.branchVec
  Sbase_MVA = net.baseMVA

  # Safety: V must be large enough to index all buses in `nodes`
  @assert length(V) >= maximum(n.busIdx for n in nodes)

  # -------------------------------------------------------------------------
  # Local helper: complex branch power S_ij (per unit) from bus `from` to `to`
  # -------------------------------------------------------------------------
  # the main loop dispatches on _branch_terminal_state and calls this helper
  # only for fully closed branches; the formula lives in the module-level
  # _closed_branch_flow_pu so the contingency screening shares it
  calcBranchFlow(from::Int, to::Int, br::Branch, tapSide::Int) = _closed_branch_flow_pu(V, from, to, br, tapSide)

  # -------------------------------------------------------------------------
  # Main loop over all branches
  # -------------------------------------------------------------------------
  ∑pfrom = 0.0
  ∑qfrom = 0.0
  ∑pto   = 0.0
  ∑qto   = 0.0
  ∑pv    = 0.0
  ∑qv    = 0.0

  for br in branchVec
    state = _branch_terminal_state(br)
    if state == :open
      @debug "Branch: $(br.comp.cName) is out of service"
      setBranchLosses!(br, 0.0, 0.0)
      br.open_end_vm_pu = nothing
      br.open_end_va_deg = nothing
      continue
    end

    if state == :open_to || state == :open_from
      closed_bus = state == :open_to ? Int(br.fromBus) : Int(br.toBus)
      if closed_bus in net.isoNodes
        # dead closed end (bus without source, see markIsolatedBuses!):
        # nothing is energized, same bookkeeping as fully open
        setBranchLosses!(br, 0.0, 0.0)
        setBranchFlow!(br, BranchFlow(nothing, nothing, 0.0, 0.0), BranchFlow(nothing, nothing, 0.0, 0.0))
        br.open_end_vm_pu = nothing
        br.open_end_va_deg = nothing
        continue
      end
      # one-sided open branch (r0.9.10): the closed terminal carries
      # S = |U|^2 * conj(Y_in) (the full charging plus the r loss of the
      # charging current), the open terminal carries zero by definition.
      # Everything entering is dissipated or stored in the branch, so the
      # branch loss equals S at the closed terminal. The open-end voltage
      # follows from the pi divider and is stored as a result on the branch
      # (Ferranti rise), the open bus itself stays whatever it is.
      closed = state == :open_to ? Int(br.fromBus) : Int(br.toBus)
      u_closed = V[closed]
      S_closed = abs2(u_closed) * conj(_open_terminal_yin(br)) * Sbase_MVA
      zeroFlow = BranchFlow(nothing, nothing, 0.0, 0.0)
      closedFlow = BranchFlow(nodes[closed]._vm_pu, nodes[closed]._va_deg, real(S_closed), imag(S_closed))
      if state == :open_to
        setBranchFlow!(br, zeroFlow, closedFlow)
      else
        setBranchFlow!(br, closedFlow, zeroFlow)
      end
      u_open = _open_end_voltage(br, u_closed)
      br.open_end_vm_pu = abs(u_open)
      br.open_end_va_deg = rad2deg(angle(u_open))
      setBranchLosses!(br, real(S_closed), imag(S_closed))
      ∑pv += real(S_closed)
      ∑qv += imag(S_closed)
      continue
    end
    br.open_end_vm_pu = nothing
    br.open_end_va_deg = nothing

    # From-side flow (from -> to)
    from = br.fromBus
    to = br.toBus
    S_from = calcBranchFlow(from, to, br, 1) * Sbase_MVA
    brFromFlow = BranchFlow(nodes[from]._vm_pu, nodes[from]._va_deg, real(S_from), imag(S_from))

    # To-side flow (to -> from)
    S_to = calcBranchFlow(to, from, br, 2) * Sbase_MVA
    brToFlow = BranchFlow(nodes[to]._vm_pu, nodes[to]._va_deg, real(S_to), imag(S_to))

    # Store flows in branch (keep existing order: to-side, from-side)
    setBranchFlow!(br, brToFlow, brFromFlow)

    # Optional accumulators (if needed later)
    # ∑pfrom += brFromFlow.pFlow
    # ∑qfrom += brFromFlow.qFlow
    # ∑pto   += brToFlow.pFlow
    # ∑qto   += brToFlow.qFlow

    # Define branch losses as sum of powers at both ends
    S_loss = S_from + S_to
    p_loss = real(S_loss)
    q_loss = imag(S_loss)

    setBranchLosses!(br, p_loss, q_loss)
    ∑pv += p_loss
    ∑qv += q_loss
  end

  setTotalLosses!(net = net, pLosses = ∑pv, qLosses = ∑qv)
end

"""
    calcLinkFlowsKCL!(net::Net; tol::Float64 = 1e-6)

Compute bus-link active/reactive flows from nodal KCL after a power-flow run,
without introducing links into the YBUS matrix. Link direction uses the
`fromBus -> toBus` sign convention.

For each bus i:
    sum(P_link,out - P_link,in) = P_inj(i) - P_branch,out(i)
(and analog for Q).

Algorithm overview:
1. Build net nodal injections `P_inj/Q_inj` from static generation minus load.
2. Subtract the solved shunt draw (`node._pShunt/_qShunt`, consumption
   positive, `updateShuntPowers!`).
3. Subtract outgoing terminal branch powers to get the link right-hand side `b`.
4. Build oriented incidence matrix `A` for active links.
5. Solve per connected link component (BFS) using `pinv(A_component) * b_component`.
6. If `sum(b_component)` is not near zero, distribute the residual uniformly
   before solving so the component system becomes consistent.
7. Write resulting link P/Q flows and derive terminal currents from |S| and V_LL.

Notes:
- `tol` is only used for the component residual-balancing step.
- For meshed/singular components (e.g., rings), `pinv` yields the minimum-norm
  least-squares solution consistent with KCL.
"""
function calcLinkFlowsKCL!(net::Net; tol::Float64 = 1e-6)
  _link_flow_allocation!(net; tol = tol)
  return nothing
end

"""
    calcLinkFlowsSE!(net; tol=1e-6) -> Union{Nothing,Vector{NamedTuple}}

Link flows after a state estimation. The nodal balances of every link
component fix its flows up to the circulations of its rings (the KCL
allocation of `calcLinkFlowsKCL!`, minimum-norm on a ring).

After `runse!` with `state_estimation.update_net = true` on this net, and
as long as the node voltages are still exactly those that run wrote, the
flows are that KCL split on the estimated state plus the ring circulations
the estimator estimated. Link flow measurements are rows of the WLS there
(see State Estimation, link measurements): a measured tree link has its
flow from the balance of the estimated state, a measured ring adds the
estimated circulation, and a link measurement the bad-data elimination
removed does not enter. `source` is `:estimated` for a link with an active
flow measurement or an estimated circulation, `:kcl` otherwise. A released
shunt enters with the susceptance the net holds: the estimate with
`state_estimation.update_shunts = true`, the model value otherwise.

Without such a run, or once the state has changed since (a `runpf!` after
the estimate), the active link flow measurements on `net.measurements`
split the flow per link component by weighted least squares under the nodal
balances as equality constraints: the measurements `f_j = z_j` with weight
`1/sigma^2`, the balances `A f = b` exact, solved on the same KKT system as
`state_estimation.zero_injection = constraint`. A measurement moves only the
circulations of a ring, a flow the balance fixes keeps its KCL value, and a
circulation no measurement sees stays at the minimum-norm split;
`source = :measured_ls` for a measured link. Without link measurements the
result equals `calcLinkFlowsKCL!` exactly (shared core).

The balance construction reads solved branch flows (`fBranchFlow`/
`tBranchFlow`, filled by `calcNetLosses!`, which `runse!(updateNet = true)`
calls) and the shunt draw (`node._pShunt/_qShunt`, refreshed by the SE
write-back).

Returns one row per active link: `(linkIdx, name, pFlow_MW, qFlow_MVar,
source, p_meas_residual, q_meas_residual)` with the residuals `f - z` of the
active link measurements on `net.measurements` (NaN when unmeasured).
`nothing` when the net has no active link.
"""
function calcLinkFlowsSE!(net::Net; tol::Float64 = 1e-6)
  pMeas = [(linkIdx = m.linkIdx, value = m.value, sigma = m.sigma) for m in net.measurements if m.active && m.linkIdx !== nothing && m.typ == PflowMeas]
  qMeas = [(linkIdx = m.linkIdx, value = m.value, sigma = m.sigma) for m in net.measurements if m.active && m.linkIdx !== nothing && m.typ == QflowMeas]
  # the circulations of the last runse! with update_net on this net, only
  # while the net still holds the state of that estimate
  circ = _se_link_circulation_current(net)
  return _link_flow_allocation!(net; tol = tol, pMeas = pMeas, qMeas = qMeas, circulation = circ)
end

# Shared allocation core. Without measurement rows this is EXACTLY the
# original KCL algorithm (same operations, bitwise-identical results); the
# measurement branches only run when a matching row exists in the component.
# With `circulation` (linkIdx => (P, Q), from a state estimation) the
# measurements do not split: the flows are the KCL split plus these
# circulations, the measurements only give the residuals.
function _link_flow_allocation!(net::Net; tol::Float64 = 1e-6, pMeas = nothing, qMeas = nothing, circulation::Union{Nothing,AbstractDict} = nothing)
  isempty(net.linkVec) && return nothing

  active = [l for l in net.linkVec if l.status == 1]
  isempty(active) && return nothing

  # link measurements by position in `active` (the incidence columns)
  activePos = Dict{Int,Int}(l.linkIdx => j for (j, l) in enumerate(active))
  pByCol = Dict{Int,Tuple{Float64,Float64}}()   # col -> (value, sigma)
  qByCol = Dict{Int,Tuple{Float64,Float64}}()
  for (meas, target) in ((pMeas, pByCol), (qMeas, qByCol))
    meas === nothing && continue
    for mm in meas
      j = get(activePos, mm.linkIdx, 0)
      if j == 0
        @warn "link-flow allocation: measurement on link $(mm.linkIdx) ignored (link not active)"
        continue
      end
      haskey(target, j) && @warn "link-flow allocation: several measurements on link $(mm.linkIdx); the last one wins"
      target[j] = (mm.value, mm.sigma)
    end
  end

  nbus = length(net.nodeVec)

  # 1) Build nodal injections from static specs in MW/MVar
  #    (same sign convention as the PF setup).
  P_inj = zeros(Float64, nbus)
  Q_inj = zeros(Float64, nbus)
  for (k, node) in enumerate(net.nodeVec)
    Pgen = isnothing(node._pƩGen) ? 0.0 : node._pƩGen
    Qgen = isnothing(node._qƩGen) ? 0.0 : node._qƩGen
    Pload = isnothing(node._pƩLoad) ? 0.0 : node._pƩLoad
    Qload = isnothing(node._qƩLoad) ? 0.0 : node._qƩLoad
    P_inj[k] = Pgen - Pload
    Q_inj[k] = Qgen - Qload
  end

  # 2) Collect outgoing branch terminal powers per bus in MW/MVar.
  #    These are already available after PF and represent non-link network exchange.
  P_branch_out = zeros(Float64, nbus)
  Q_branch_out = zeros(Float64, nbus)
  for br in net.branchVec
    # only fully open branches carry no terminal power; a one-sided open
    # branch (r0.9.10) contributes its closed-terminal charging draw here
    # (the stored open-side flow is zero, adding it is a no-op)
    _branch_terminal_state(br) == :open && continue
    if !isnothing(br.fBranchFlow)
      P_branch_out[br.fromBus] += br.fBranchFlow.pFlow
      Q_branch_out[br.fromBus] += br.fBranchFlow.qFlow
    end
    if !isnothing(br.tBranchFlow)
      P_branch_out[br.toBus] += br.tBranchFlow.pFlow
      Q_branch_out[br.toBus] += br.tBranchFlow.qFlow
    end
  end

  # 3) Shunt draw at the bus (results from PF or SE) if available.
  #    node._pShunt/node._qShunt are the power the shunts consume
  #    (updateShuntPowers!: |V|^2 conj(y)), so they leave the balance.
  #    Until 0.60.0 they were added: a reactor at a link bus moved its
  #    Mvar to the wrong side (single link, reactor 8 Mvar at 1 pu on the to
  #    bus: 16.20 Mvar instead of the 23.49 the feeding line delivers).
  P_shunt = zeros(Float64, nbus)
  Q_shunt = zeros(Float64, nbus)
  for (k, node) in enumerate(net.nodeVec)
    P_shunt[k] = isnothing(node._pShunt) ? 0.0 : node._pShunt
    Q_shunt[k] = isnothing(node._qShunt) ? 0.0 : node._qShunt
  end

  # 4) Remaining nodal balance that must be carried by links.
  bP = P_inj .- P_shunt .- P_branch_out
  bQ = Q_inj .- Q_shunt .- Q_branch_out

  m = length(active)
  A = zeros(Float64, nbus, m)
  for (j, l) in enumerate(active)
    A[l.fromBus, j] += 1.0
    A[l.toBus, j] -= 1.0
  end

  # 5) Split by connected components in the link graph for robust solving.
  #    Each component can be solved independently.
  bus_to_links = [Int[] for _ in 1:nbus]
  for (j, l) in enumerate(active)
    push!(bus_to_links[l.fromBus], j)
    push!(bus_to_links[l.toBus], j)
  end

  visited = falses(nbus)
  flowsP = zeros(Float64, m)
  flowsQ = zeros(Float64, m)

  for b0 in 1:nbus
    if visited[b0] || isempty(bus_to_links[b0])
      continue
    end

    # BFS buses in this link component and collect all incident links.
    q = [b0]
    comp_buses = Int[]
    comp_links_set = Set{Int}()
    visited[b0] = true

    while !isempty(q)
      b = popfirst!(q)
      push!(comp_buses, b)
      for lj in bus_to_links[b]
        push!(comp_links_set, lj)
        l = active[lj]
        nb = (l.fromBus == b) ? l.toBus : l.fromBus
        if !visited[nb]
          visited[nb] = true
          push!(q, nb)
        end
      end
    end

    comp_links = sort!(collect(comp_links_set))
    Ab = A[comp_buses, comp_links]
    bbP = copy(bP[comp_buses])
    bbQ = copy(bQ[comp_buses])

    # 6) Enforce component solvability condition sum(b)=0.
    #    If a small mismatch exists (numerics/modeling residual), distribute it
    #    uniformly across component buses.
    sumP = sum(bbP)
    sumQ = sum(bbQ)
    if abs(sumP) > tol
      bbP .-= sumP / length(bbP)
    end
    if abs(sumQ) > tol
      bbQ .-= sumQ / length(bbQ)
    end

    # 7) Least-squares solve. `pinv` also handles singular meshed components
    #    (e.g., rings) by returning the minimum-norm solution. With link
    #    measurements in this component the split is the equality-constrained
    #    least squares of `_constrained_link_solve`, on b AFTER the sum-zero
    #    correction (the balances exact, minimum norm where no measurement
    #    reaches).
    # with circulations from a state estimation the split is that estimate,
    # the measurements only give the residuals below
    compP = circulation === nothing ? [(k, pByCol[gj]) for (k, gj) in enumerate(comp_links) if haskey(pByCol, gj)] : Tuple{Int,Tuple{Float64,Float64}}[]
    compQ = circulation === nothing ? [(k, qByCol[gj]) for (k, gj) in enumerate(comp_links) if haskey(qByCol, gj)] : Tuple{Int,Tuple{Float64,Float64}}[]
    fp = if isempty(compP)
      pinv(Ab) * bbP
    else
      _constrained_link_solve(Ab, bbP, compP)
    end
    fq = if isempty(compQ)
      pinv(Ab) * bbQ
    else
      _constrained_link_solve(Ab, bbQ, compQ)
    end

    for (k, gj) in enumerate(comp_links)
      flowsP[gj] = fp[k]
      flowsQ[gj] = fq[k]
    end
  end

  # 8) Reset all link outputs, then write back solved values for active links.
  for l in net.linkVec
    setLinkFlow!(l, 0.0, 0.0)
    setLinkCurrent!(l, 0.0, 0.0)
  end
  rows = NamedTuple[]
  for (j, l) in enumerate(active)
    cpq = circulation === nothing ? nothing : get(circulation, l.linkIdx, nothing)
    p = cpq === nothing ? flowsP[j] : flowsP[j] + cpq[1]
    q = cpq === nothing ? flowsQ[j] : flowsQ[j] + cpq[2]
    setLinkFlow!(l, p, q)

    # 9) Three-phase current magnitudes per terminal:
    # |I|[kA] = |S|[MVA] / (sqrt(3) * V_LL[kV])
    s_abs = hypot(p, q)
    vf_kV = net.nodeVec[l.fromBus]._vm_pu * net.nodeVec[l.fromBus].comp.cVN
    vt_kV = net.nodeVec[l.toBus]._vm_pu * net.nodeVec[l.toBus].comp.cVN

    ifrom = (vf_kV > 1e-12) ? (s_abs / (Wurzel3 * vf_kV)) : NaN
    ito = (vt_kV > 1e-12) ? (s_abs / (Wurzel3 * vt_kV)) : NaN
    setLinkCurrent!(l, ifrom, ito)

    # report row: allocation source and, where measured, the residual f - z
    measured = haskey(pByCol, j) || haskey(qByCol, j)
    source = circulation === nothing ? (measured ? :measured_ls : :kcl) : ((measured || cpq !== nothing) ? :estimated : :kcl)
    push!(
      rows,
      (
        linkIdx = l.linkIdx,
        name = l.cName,
        pFlow_MW = p,
        qFlow_MVar = q,
        source = source,
        p_meas_residual = haskey(pByCol, j) ? p - pByCol[j][1] : NaN,
        q_meas_residual = haskey(qByCol, j) ? q - qByCol[j][1] : NaN,
      ),
    )
  end
  return rows
end

# Cut for the singular values in `_link_measured_cycles`. N has orthonormal
# columns, so the singular values of N[measured, :] lie in [0, 1] and the
# cut is absolute: a link outside every loop (a pendant coupler) has a zero
# row in exact arithmetic, nullspace returns rounding instead (4.8e-18 on a
# ring of three links with one pendant link), and a cut relative to the
# largest singular value counted that as a loop.
const _LINK_CYCLE_TOL = 1.0e-8

# The ring circulations a set of measured links sees. `N` is nullspace(A)
# of a link component (one orthonormal column per independent loop),
# `measured` the positions of the measured links in A. Returns (Nc, Nfree):
# Nc spans the loop directions that move at least one measured flow, Nfree
# those that move none. Entries at rounding level (|a| <= 1e-12) are set to
# zero, so a link outside the loops gets no circulation term. Used by the
# split below and by the link rows of the state estimation
# (`_se_link_table`), one rank rule for both.
function _link_measured_cycles(N::AbstractMatrix{Float64}, measured::Vector{Int})
  n = size(N, 2)
  (n == 0 || isempty(measured)) && return (zeros(Float64, size(N, 1), 0), Matrix{Float64}(N))
  F = svd(N[measured, :]; full = true)
  r = count(s -> s > _LINK_CYCLE_TOL, F.S)
  Nc = N * F.V[:, 1:r]
  Nfree = N * F.V[:, (r+1):n]
  for M in (Nc, Nfree)
    M[abs.(M).<=1e-12] .= 0.0
  end
  return Nc, Nfree
end

# Link split of one component with measurements: the flows f minimize
# sum_j w_j (f_j - z_j)^2 with w_j = 1/sigma_j^2 under the nodal balances
# A f = b as equality constraints, on the KKT system of
# state_estimation.zero_injection = constraint (`_se_solve_constrained!`;
# the problem is linear, one step from f = 0 is the solution). Constraint
# rows: the balances without the last bus, on b minus its mean (the part of
# b outside the range of A, which pinv ignores as well; the caller removes
# it only above `tol`). The component is connected, so the last row then
# follows from the others and would make the KKT matrix singular. Plus
# Nfree' f = 0 for the loop directions no measurement sees, so they keep
# the minimum-norm split that pinv gives without measurements. A singular
# KKT matrix can
# only come from a construction error here; the solver then falls back to
# least squares and the split warns.
function _constrained_link_solve(Ab::Matrix{Float64}, bb::Vector{Float64}, meas::Vector{<:Tuple{Int,Tuple{Float64,Float64}}})
  nrow, ncol = size(Ab)
  G = zeros(Float64, ncol, ncol)
  g = zeros(Float64, ncol)
  for (col, vz) in meas
    w = vz[2] > 0.0 ? inv(vz[2]^2) : 1.0
    G[col, col] += w
    g[col] += w * vz[1]
  end
  _, Nfree = _link_measured_cycles(nullspace(Ab), sort!(unique(Int[col for (col, _) in meas])))
  C = vcat(Ab[1:(nrow-1), :], permutedims(Nfree))
  bc = bb .- sum(bb) / nrow
  rc = vcat(bc[1:(nrow-1)], zeros(Float64, size(Nfree, 2)))
  f, _, fell = _se_solve_constrained!(_SENormalEquationsSolver(), sparse(G), g, sparse(C), rc)
  fell && @warn "calcLinkFlowsSE!: the KKT matrix of a link split was singular (a bug in the split); it was solved by least squares"
  return f
end

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
  # tap handling (magnitude + angle)
  ratio = (br.ratio != 0.0) ? br.ratio : 1.0
  angle = (br.ratio != 0.0) ? br.angle : 0.0
  tap = calcComplexRatio(tapRatio = ratio, angleInDegrees = angle)
  if tapSide == 1
    ui /= tap
  elseif tapSide == 2
    uj /= tap
  end
  Yik = inv(br.r_pu + im * br.x_pu)
  Y0ik = 0.5 * (br.g_pu + im * br.b_pu)
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
2. Add solved shunt injections (`node._pShunt/_qShunt`) to nodal injections.
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

W2 link-flow allocation with measurements (SE phase 3): the KCL allocation of
`calcLinkFlowsKCL!` extended to a weighted least-squares split per link
component. Active link flow measurements on `net.measurements` (`PflowMeas`/
`QflowMeas` with a `linkIdx`) enter as extra rows `f_j = z_j` with weight
`1/sigma^2` (P rows into the P solve, Q rows into the Q solve); the nodal
balance rows keep weight 1. With no link measurements the result equals
`calcLinkFlowsKCL!` exactly (shared core, identical code path).

Link measurements constrain only the flow split, never the system state: the
estimator excludes them from the WLS (the link is not in the Ybus).

The balance construction reads solved branch flows (`fBranchFlow`/
`tBranchFlow`, filled by `calcNetLosses!`, which `runse!(updateNet = true)`
calls) and the shunt injections (`node._pShunt/_qShunt`, refreshed by the
SE write-back). Run after `runse!` with `updateNet = true`.

Returns one row per active link: `(linkIdx, name, pFlow_MW, qFlow_MVar,
source, p_meas_residual, q_meas_residual)` with `source = :kcl` (no
measurement on that link) or `:measured_ls`, and the measurement residuals
`f_est - z` (NaN when unmeasured). `nothing` when the net has no active
link.
"""
function calcLinkFlowsSE!(net::Net; tol::Float64 = 1e-6)
  pMeas = [(linkIdx = m.linkIdx, value = m.value, sigma = m.sigma) for m in net.measurements if m.active && m.linkIdx !== nothing && m.typ == PflowMeas]
  qMeas = [(linkIdx = m.linkIdx, value = m.value, sigma = m.sigma) for m in net.measurements if m.active && m.linkIdx !== nothing && m.typ == QflowMeas]
  return _link_flow_allocation!(net; tol = tol, pMeas = pMeas, qMeas = qMeas)
end

# Shared allocation core. Without measurement rows this is EXACTLY the
# original KCL algorithm (same operations, bitwise-identical results); the
# measurement branches only run when a matching row exists in the component.
function _link_flow_allocation!(net::Net; tol::Float64 = 1e-6, pMeas = nothing, qMeas = nothing)
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

  # 3) Include shunt bus injections (results from PF) if available.
  #    node._pShunt/node._qShunt are signed injections at the bus.
  P_shunt = zeros(Float64, nbus)
  Q_shunt = zeros(Float64, nbus)
  for (k, node) in enumerate(net.nodeVec)
    P_shunt[k] = isnothing(node._pShunt) ? 0.0 : node._pShunt
    Q_shunt[k] = isnothing(node._qShunt) ? 0.0 : node._qShunt
  end

  # 4) Remaining nodal balance that must be carried by links.
  bP = P_inj .+ P_shunt .- P_branch_out
  bQ = Q_inj .+ Q_shunt .- Q_branch_out

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
    #    measurements in this component (W2, SE phase 3) one weighted row
    #    `f_j = z_j` per measurement is appended AFTER the sum-zero
    #    correction; the balance rows keep weight 1 and the minimum-norm
    #    property survives through the pinv of the row-scaled system.
    compP = [(k, pByCol[gj]) for (k, gj) in enumerate(comp_links) if haskey(pByCol, gj)]
    compQ = [(k, qByCol[gj]) for (k, gj) in enumerate(comp_links) if haskey(qByCol, gj)]
    fp = if isempty(compP)
      pinv(Ab) * bbP
    else
      _weighted_link_solve(Ab, bbP, compP)
    end
    fq = if isempty(compQ)
      pinv(Ab) * bbQ
    else
      _weighted_link_solve(Ab, bbQ, compQ)
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
    p = flowsP[j]
    q = flowsQ[j]
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
    push!(
      rows,
      (
        linkIdx = l.linkIdx,
        name = l.cName,
        pFlow_MW = p,
        qFlow_MVar = q,
        source = measured ? :measured_ls : :kcl,
        p_meas_residual = haskey(pByCol, j) ? p - pByCol[j][1] : NaN,
        q_meas_residual = haskey(qByCol, j) ? q - qByCol[j][1] : NaN,
      ),
    )
  end
  return rows
end

# Weighted least-squares split for one component: balance rows (weight 1) plus
# one row per measured link, scaled by sqrt(1/sigma^2). pinv keeps minimum
# norm in the directions no row constrains.
function _weighted_link_solve(Ab::Matrix{Float64}, bb::Vector{Float64}, meas::Vector{<:Tuple{Int,Tuple{Float64,Float64}}})
  nrow, ncol = size(Ab)
  M = zeros(Float64, nrow + length(meas), ncol)
  rhs = zeros(Float64, nrow + length(meas))
  M[1:nrow, :] .= Ab
  rhs[1:nrow] .= bb
  for (r, (col, vz)) in enumerate(meas)
    z, σ = vz
    sw = σ > 0.0 ? 1.0 / σ : 1.0
    M[nrow+r, col] = sw
    rhs[nrow+r] = sw * z
  end
  return pinv(M) * rhs
end

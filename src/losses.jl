# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 04.09.2023
# file losses.jl

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
  function calcBranchFlow(from::Int, to::Int, br::Branch, tapSide::Int)
    @assert tapSide == 1 || tapSide == 2
    if br.status < 0
      return 0.0 + 0.0im
    end

    # Base voltages taken directly from V (already includes vm and va)
    ui = V[from]
    uj = V[to]

    # Tap handling (magnitude + angle)
    ratio = (br.ratio != 0.0) ? br.ratio : 1.0
    angle = (br.ratio != 0.0) ? br.angle : 0.0
    tap   = calcComplexRatio(tapRatio = ratio, angleInDegrees = angle)

    if tapSide == 1
      ui /= tap
    elseif tapSide == 2
      uj /= tap
    end

    rpu = br.r_pu
    xpu = br.x_pu
    bpu = br.b_pu
    gpu = br.g_pu

    Yik  = inv(rpu + im * xpu)
    Y0ik = 0.5 * (gpu + im * bpu)

    # Complex power S_ij in per unit
    s = abs(ui)^2 * conj(Y0ik + Yik) - ui * conj(uj) * conj(Yik)

    return s
  end # calcBranchFlow

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
    if br.status == 0
      @debug "Branch: $(br.comp.cName) is out of service"
      setBranchLosses!(br, 0.0, 0.0)
      continue
    end

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

Computes active bus-link flows after power flow via nodal KCL, without adding
links to YBUS. Link directions follow `fromBus -> toBus` sign convention.

For each bus i:
    sum(P_link,out - P_link,in) = P_inj(i) - P_branch,out(i)
(and analog for Q).

If a connected link component has a small non-zero residual sum, the residual is
uniformly distributed across buses in that component to enforce solvability.
"""
function calcLinkFlowsKCL!(net::Net; tol::Float64 = 1e-6)
  isempty(net.linkVec) && return

  active = [l for l in net.linkVec if l.status == 1]
  isempty(active) && return

  nbus = length(net.nodeVec)

  # Bus injections from static specs in MW/MVar (same sign convention as PF setup)
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

  # Outgoing branch terminal powers per bus in MW/MVar
  P_branch_out = zeros(Float64, nbus)
  Q_branch_out = zeros(Float64, nbus)
  for br in net.branchVec
    br.status == 0 && continue
    if !isnothing(br.fBranchFlow)
      P_branch_out[br.fromBus] += br.fBranchFlow.pFlow
      Q_branch_out[br.fromBus] += br.fBranchFlow.qFlow
    end
    if !isnothing(br.tBranchFlow)
      P_branch_out[br.toBus] += br.tBranchFlow.pFlow
      Q_branch_out[br.toBus] += br.tBranchFlow.qFlow
    end
  end

  bP = P_inj .- P_branch_out
  bQ = Q_inj .- Q_branch_out

  m = length(active)
  A = zeros(Float64, nbus, m)
  for (j, l) in enumerate(active)
    A[l.fromBus, j] += 1.0
    A[l.toBus, j] -= 1.0
  end

  # Split by connected components in link graph for robust solving
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

    # BFS buses in this link component
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

    # enforce sum(b)=0 for each connected component
    sumP = sum(bbP)
    sumQ = sum(bbQ)
    if abs(sumP) > tol
      bbP .-= sumP / length(bbP)
    end
    if abs(sumQ) > tol
      bbQ .-= sumQ / length(bbQ)
    end

    # least-squares (handles tree and meshed components)
    fp = Ab \ bbP
    fq = Ab \ bbQ

    for (k, gj) in enumerate(comp_links)
      flowsP[gj] = fp[k]
      flowsQ[gj] = fq[k]
    end
  end

  # reset + write back to original link objects
  for l in net.linkVec
    setLinkFlow!(l, 0.0, 0.0)
    setLinkCurrent!(l, 0.0, 0.0)
  end
  for (j, l) in enumerate(active)
    p = flowsP[j]
    q = flowsQ[j]
    setLinkFlow!(l, p, q)

    # Three-phase current magnitudes per terminal:
    # |I|[kA] = |S|[MVA] / (sqrt(3) * V_LL[kV])
    s_abs = hypot(p, q)
    vf_kV = net.nodeVec[l.fromBus]._vm_pu * net.nodeVec[l.fromBus].comp.cVN
    vt_kV = net.nodeVec[l.toBus]._vm_pu * net.nodeVec[l.toBus].comp.cVN

    ifrom = (vf_kV > 1e-12) ? (s_abs / (Wurzel3 * vf_kV)) : NaN
    ito = (vt_kV > 1e-12) ? (s_abs / (Wurzel3 * vt_kV)) : NaN
    setLinkCurrent!(l, ifrom, ito)
  end
end

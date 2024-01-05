# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 04.09.2023
# include-file losses.jl

function calcNetLosses!(nodes::Vector{ResDataTypes.Node}, branchVec::Vector{ResDataTypes.Branch}, Sbase_MVA::Float64, log::Bool = false)
  if debug
    println("\ncalcNetworkLosses (BaseMVA=$(Sbase_MVA))\n")
  end

  # Sij = vi*exp(j*phi_i)*( (vi*exp(j*phi_i) - vk*exp(j*phi_k)*Y_ik +  vi*exp(j*phi_i)*Y0ik)*                     
  # Y0ik: Queradmittanz
  function calcBranchFlow(from::Int, to::Int, br::ResDataTypes.Branch, tapSide::Int)
    argi = deg2rad(nodes[from]._va_deg)
    argj = deg2rad(nodes[to]._va_deg)

    tap = (br.ratio != 0.0) ? br.ratio : 1.0

    ui = nodes[from]._vm_pu * exp(im * argi)
    if tapSide == 1
      ui = ui / tap
    end

    uj = nodes[to]._vm_pu * exp(im * argj)
    if tapSide == 2
      uj = uj / tap
    end

    u_diff = ui - uj
    Yik = inv((br.r_pu + im * br.x_pu))
    Y0ik = 0.5 * (br.g_pu + im * br.b_pu)
    s = ui * conj(u_diff * Yik + ui * Y0ik)
    return (s)
  end

  n = length(nodes)
  Pvk = zeros(n)
  Qvk = zeros(n)

  p_total = 0.0
  q_total = 0.0
  if debug
    println("Branch Power Flows:")
  end

  # Parallel branches, only 2 branches per tuple!
  branchTupleSet = Set{Tuple}()
  for br in branchVec
    if br.isParallel
      @debug "Parallel branch detected: $(br.fromBus) -> $(br.toBus)"
      tupple = (br.fromBus, br.toBus)
      push!(branchTupleSet, tupple)
    end
  end

  for br in branchVec
    tupple = (br.fromBus, br.toBus)
    pCorr = 1.0
    if tupple in branchTupleSet
      @debug "taking parallel branch $(br.fromBus) -> $(br.toBus) into account"
      pCorr = 0.5
    end

    from = br.fromBus
    to = br.toBus
    S = calcBranchFlow(from, to, br, 1) * Sbase_MVA
    P = real(S) * pCorr
    Q = imag(S) * pCorr

    Pvk[from] += P
    Qvk[from] += Q

    brFromFlow = BranchFlow(nodes[br.fromBus]._vm_pu, nodes[br.fromBus]._va_deg, P, Q)
    if debug
      p_val = round(P, digits = 2)
      q_val = round(Q, digits = 2)
      println("Bus_$(from) -> Bus_$(to): P = $(p_val), Q = $(q_val)")
    end
    p_total += P
    q_total += Q

    from = br.toBus
    to = br.fromBus
    S = calcBranchFlow(from, to, br, 2) * Sbase_MVA
    P = real(S) * pCorr
    Q = imag(S) * pCorr

    Pvk[from] += P
    Qvk[from] += Q

    p_total += P
    q_total += Q
    if debug
      p_val = round(P, digits = 2)
      q_val = round(Q, digits = 2)
      println("Bus_$(from) -> Bus_$(to): P = $(p_val), Q = $(q_val)\n")
    end
    brToFlow = BranchFlow(nodes[br.toBus]._vm_pu, nodes[br.toBus]._va_deg, P, Q)
    setBranchFlow!(brToFlow, brFromFlow, br)
  end
  if debug
    p_val = round(p_total, digits = 2)
    q_val = round(q_total, digits = 2)
    println("Losses:")
    println("∑Pv = $(p_val), ∑Qv = $(q_val)\n")
    println("Bus Powers:")
    for i = 1:n
      p_val = round(Pvk[i], digits = 3)
      q_val = round(Qvk[i], digits = 3)
      println("Bus_[$(i)]: P = $(p_val), Q = $(q_val) ")
    end
  end

  #TODO: Node PQ vs. GenPower
  for n in nodes
    if Pvk[n.busIdx] <= 0.0
      setNodePQ!(n, -Pvk[n.busIdx], -Qvk[n.busIdx])
    elseif Pvk[n.busIdx] > 0.0
      setGenPower!(n, Pvk[n.busIdx], Qvk[n.busIdx])
    end
  end
end

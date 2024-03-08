# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 04.09.2023
# include-file losses.jl

function calcNetLosses!(nodes::Vector{ResDataTypes.Node}, branchVec::Vector{ResDataTypes.Branch}, Sbase_MVA::Float64)
  @debug "\ncalcNetworkLosses (BaseMVA=$(Sbase_MVA))\n"
  # Sij = vi*exp(j*phi_i)*( (vi*exp(j*phi_i) - vk*exp(j*phi_k)*Y_ik +  vi*exp(j*phi_i)*Y0ik)*  
  # Sij = vi^2*conj(Y0ik+Yik)-vi*conj(vk)*conj(Yik)                   
  # Y0ik: Queradmittanz
  function calcBranchFlow(from::Int, to::Int, br::ResDataTypes.Branch, tapSide::Int)
    @assert tapSide == 1 || tapSide == 2
    if br.status < 0
      return (0.0 + 0.0im)
    end

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
    rpu = br.r_pu
    xpu = br.x_pu
    bpu = br.b_pu
    gpu = br.g_pu

    Yik = inv((rpu + im * xpu))      
    Y0ik = 0.5*(gpu + im * bpu)
    s = abs(ui)^2*conj(Y0ik+Yik)-ui*conj(uj)*conj(Yik) 
    
    return s
  end # calcBranchFlow

  n = length(nodes)
  
  #∑pfrom=∑qfrom=∑pto=∑qto=0.0
  for br in branchVec    
    from = br.fromBus
    to = br.toBus
    
    S = calcBranchFlow(from, to, br, 1) * Sbase_MVA
    


    brFromFlow = BranchFlow(nodes[br.fromBus]._vm_pu, nodes[br.fromBus]._va_deg, real(S), imag(S))

    from = br.toBus
    to = br.fromBus
    S = calcBranchFlow(from, to, br, 2) * Sbase_MVA    
    

    brToFlow = BranchFlow(nodes[br.toBus]._vm_pu, nodes[br.toBus]._va_deg, real(S), imag(S))
    setBranchFlow!(brToFlow, brFromFlow, br)

    #=
    ∑pfrom += brFromFlow.pFlow
    ∑qfrom += brFromFlow.qFlow
    ∑pto += brToFlow.pFlow
    ∑qto += brToFlow.qFlow
    =#
  end

  #=
  ∑pv = ∑pfrom+∑pto
  ∑qv = ∑qfrom+∑qto
  ∑pl = ∑qg = ∑pg = ∑ql = 0.0
  idx = 0
  for n in nodes
    if n._nodeType == ResDataTypes.Slack
      idx = n.busIdx
      continue
    end
    ∑pl += isnothing(n._pƩLoad) ? 0.0 : n._pƩLoad
    ∑ql += isnothing(n._qƩLoad) ? 0.0 : n._qƩLoad
    ∑pg += isnothing(n._pƩGen) ? 0.0 : n._pƩGen
    ∑qg += isnothing(n._qƩGen) ? 0.0 : n._qƩGen
  end
  
  pSlack = abs(∑pg-∑pl) + ∑pv
  qSlack = abs(∑qg-∑ql) + ∑qv
  if pSlack > 1e-12
    @debug "Slack-Generator is feeding power into the network (P=$pSlack MW, Q=$qSlack MVar)"
    setGenPower!(nodes[idx], pSlack, qSlack)
  elseif pSlack < -1e-12
    @debug "Slack-Generator is consuming power from the network (P=$pSlack MW, Q=$qSlack MVar)"
    setNodePQ!(nodes[idx], pSlack, qSlack)
  else
    @debug "Slack-Generator is not feeding or consuming power (P=$pSlack MW, Q=$qSlack MVar)"
    setNodePQ!(nodes[idx], 0.0, 0.0)
    setGenPower!(nodes[idx], 0.0, 0.0)
  end  
  =#
  
end

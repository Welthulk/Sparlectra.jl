# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 04.09.2023
# include-file losses.jl

"""
    calcNetLosses!(net::Net)

Calculates the network losses for the given network.

# Arguments
- `net::Net`: The network.

# Example
```julia
calcNetLosses!(net = network)
"""
function calcNetLosses!(net::Net)
  @debug "\ncalcNetworkLosses\n"
  # Sij = vi*exp(j*phi_i)*( (vi*exp(j*phi_i) - vk*exp(j*phi_k)*Y_ik +  vi*exp(j*phi_i)*Y0ik)*  
  # Sij = vi^2*conj(Y0ik+Yik)-vi*conj(vk)*conj(Yik)                   
  # Y0ik: Queradmittanz
  nodes = net.nodeVec
  branchVec = net.branchVec
  Sbase_MVA = net.baseMVA

  function calcBranchFlow(from::Int, to::Int, br::Branch, tapSide::Int)
    @assert tapSide == 1 || tapSide == 2
    if br.status < 0
      return (0.0 + 0.0im)
    end

    argi = deg2rad(nodes[from]._va_deg)
    argj = deg2rad(nodes[to]._va_deg)

    ratio = (br.ratio != 0.0) ? br.ratio : 1.0
    angle = (br.ratio != 0.0) ? br.angle : 0.0

    tap = calcComplexRatio(tapRatio = ratio, angleInDegrees = angle)

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
    Y0ik = 0.5 * (gpu + im * bpu)
    s = abs(ui)^2 * conj(Y0ik + Yik) - ui * conj(uj) * conj(Yik)

    return s
  end # calcBranchFlow

  function calcLosses(br::Branch)
    VmFrom = br.fBranchFlow.vm_pu
    VaFrom = br.fBranchFlow.va_deg
    VmTo = br.tBranchFlow.vm_pu
    VaTo = br.tBranchFlow.va_deg

    v1 = VmFrom * exp(im * deg2rad(VaFrom))
    v2 = VmTo * exp(im * deg2rad(VaTo))
    # Losses: abs( Vf / tau - Vt ) ^ 2 / (Rs - j Xs)
    ratio = (br.ratio != 0.0) ? br.ratio : 1.0
    angle = (br.ratio != 0.0) ? br.angle : 0.0
    tap = calcComplexRatio(tapRatio = ratio, angleInDegrees = angle)

    Vdiff = v1 / tap - v2

    y = inv((br.r_pu + br.x_pu * im)) #+ 0.5*(br.g_pu + br.b_pu*im)

    Sdiff = Vdiff * conj(Vdiff * y) * net.baseMVA
    pv = real(Sdiff)
    qv = imag(Sdiff)
    return pv, qv
  end

  n = length(nodes)

  ∑pfrom = ∑qfrom = ∑pv = ∑qv = ∑pto = ∑qto = 0.0
  for br in branchVec
    from = br.fromBus
    to = br.toBus
    if br.status == 0
      @debug "Branch: $(br.comp.cName) is out of service"
      setBranchLosses!(br, 0.0, 0.0)
      continue
    end
    S = calcBranchFlow(from, to, br, 1) * Sbase_MVA
    brFromFlow = BranchFlow(nodes[br.fromBus]._vm_pu, nodes[br.fromBus]._va_deg, real(S), imag(S))

    from = br.toBus
    to = br.fromBus
    S = calcBranchFlow(from, to, br, 2) * Sbase_MVA
    brToFlow = BranchFlow(nodes[br.toBus]._vm_pu, nodes[br.toBus]._va_deg, real(S), imag(S))

    setBranchFlow!(br, brToFlow, brFromFlow)

    #= for later use...
    ∑pfrom += brFromFlow.pFlow
    ∑qfrom += brFromFlow.qFlow
    ∑pto += brToFlow.pFlow
    ∑qto += brToFlow.qFlow
    =#

    pv, qv = calcLosses(br)
    setBranchLosses!(br, abs(pv), abs(qv))
    ∑pv += abs(pv)
    ∑qv += abs(qv)
  end
  setTotalLosses!(net = net, pLosses = ∑pv, qLosses = ∑qv)
end

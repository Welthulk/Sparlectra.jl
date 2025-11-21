# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 04.09.2023
# include-file losses.jl
# Sij = vi*exp(j*phi_i)*( (vi*exp(j*phi_i) - vk*exp(j*phi_k)*Y_ik +  vi*exp(j*phi_i)*Y0ik)*
# Sij = vi^2*conj(Y0ik+Yik)-vi*conj(vk)*conj(Yik)
# Y0ik: shunt admittance

"""
    buildVoltageVector(net::Net) -> Vector{ComplexF64}

Builds the complex bus voltage vector V[k] = vm_pu[k] * exp(j * va_rad[k])
using the current nodal state stored in `net.nodeVec`.
"""
function buildVoltageVector(net::Net)
    V = Vector{ComplexF64}(undef, length(net.nodeVec))
    @inbounds for n in net.nodeVec
        V[n.busIdx] = n._vm_pu * cis(deg2rad(n._va_deg))
    end
    return V
end

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
    # Local helper: branch losses based on bus voltages V and series admittance
    # -------------------------------------------------------------------------
    function calcLosses(br::Branch)
        # Use complex bus voltages directly from V (no recomputation from vm/va)
        v1 = V[br.fromBus]
        v2 = V[br.toBus]

        ratio = (br.ratio != 0.0) ? br.ratio : 1.0
        angle = (br.ratio != 0.0) ? br.angle : 0.0
        tap   = calcComplexRatio(tapRatio = ratio, angleInDegrees = angle)

        # Losses: | V1/tap - V2 |^2 / (Rs - j Xs)
        Vdiff = v1 / tap - v2

        y = inv(br.r_pu + br.x_pu * im)  # + 0.5*(br.g_pu + br.b_pu*im) if shunt is included
        Sdiff = Vdiff * conj(Vdiff * y) * Sbase_MVA

        pv = real(Sdiff)
        qv = imag(Sdiff)
        return pv, qv
    end

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

        # From-side flow
        from = br.fromBus
        to   = br.toBus
        S_from = calcBranchFlow(from, to, br, 1) * Sbase_MVA
        brFromFlow = BranchFlow(
            nodes[from]._vm_pu,
            nodes[from]._va_deg,
            real(S_from),
            imag(S_from),
        )

        # To-side flow
        from = br.toBus
        to   = br.fromBus
        S_to = calcBranchFlow(from, to, br, 2) * Sbase_MVA
        brToFlow = BranchFlow(
            nodes[br.toBus]._vm_pu,
            nodes[br.toBus]._va_deg,
            real(S_to),
            imag(S_to),
        )

        # Store flows in branch (keep existing order)
        setBranchFlow!(br, brToFlow, brFromFlow)

        # (accumulators kept for potential later use)
        # ∑pfrom += brFromFlow.pFlow
        # ∑qfrom += brFromFlow.qFlow
        # ∑pto   += brToFlow.pFlow
        # ∑qto   += brToFlow.qFlow

        pv, qv = calcLosses(br)
        # Keep behavior: store absolute values of P/Q losses
        setBranchLosses!(br, abs(pv), abs(qv))
        ∑pv += abs(pv)
        ∑qv += abs(qv)
    end

    setTotalLosses!(net = net, pLosses = ∑pv, qLosses = ∑qv)
end

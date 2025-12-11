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
        to   = br.toBus
        S_from = calcBranchFlow(from, to, br, 1) * Sbase_MVA
        brFromFlow = BranchFlow(
            nodes[from]._vm_pu,
            nodes[from]._va_deg,
            real(S_from),
            imag(S_from),
        )

        # To-side flow (to -> from)
        S_to = calcBranchFlow(to, from, br, 2) * Sbase_MVA
        brToFlow = BranchFlow(
            nodes[to]._vm_pu,
            nodes[to]._va_deg,
            real(S_to),
            imag(S_to),
        )

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

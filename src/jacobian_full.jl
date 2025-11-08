# jacobian_full.jl — Full-system Newton-Raphson including PV identity rows (minimal invasive)

using LinearAlgebra
using SparseArrays
using Printf

# ------------------------------
# Power feed vector for the full system
# (P for all buses, Q only for PQ; PV gets placeholder rows)
# ------------------------------
function getPowerFeeds_full(busVec::Vector{BusData}, n_pq::Int, n_pv::Int)::Vector{Float64}
    n = n_pq + n_pv
    size = 2*n
    pf = zeros(Float64, size)
    vindex = 0
    @inbounds for bus in busVec
        if bus.type == PQ
            vindex += 1; pf[vindex] = bus.pƩ   # P
            vindex += 1; pf[vindex] = bus.qƩ   # Q
        elseif bus.type == PV
            vindex += 1; pf[vindex] = bus.pƩ   # P
            vindex += 1; pf[vindex] = 0.0      # rPV placeholder
        end
    end
    return pf
end

# ------------------------------
# Residual for the full system: PV’s second equation becomes rPV = V − Vset
# ------------------------------
function residuum_full_withPV(
    Y::AbstractMatrix{ComplexF64},
    busVec::Vector{BusData},
    Vset::Vector{Float64},           # length = number of buses (only relevant for PV)
    n_pq::Int, n_pv::Int,
    log::Bool)::Vector{Float64}

    V = [bus.vm_pu * exp(im * bus.va_rad) for bus in busVec]
    S = Diagonal(V) * conj(Y * V)

    n = n_pq + n_pv
    Δ = zeros(Float64, 2*n)

    vindex = 0
    @inbounds for k in eachindex(busVec)
        bus = busVec[k]
        if bus.type == PQ
            vindex += 1; Δ[vindex] = bus.pƩ - real(S[k])        # ΔP
            vindex += 1; Δ[vindex] = bus.qƩ - imag(S[k])        # ΔQ
        elseif bus.type == PV
            vindex += 1; Δ[vindex] = bus.pƩ - real(S[k])        # ΔP
            vindex += 1; Δ[vindex] = bus.vm_pu - Vset[k]        # rPV
        end
        bus._pRes = real(S[k]); bus._qRes = imag(S[k])
    end

    if log
        println("residuum_full_withPV: ‖Δ‖ = ", norm(Δ))
    end
    return Δ
end

# ------------------------------
# Utility: set a full row of a CSC matrix to zero
# ------------------------------
function zero_row!(J::SparseMatrixCSC, r::Int)
    @inbounds for c in axes(J, 2)
        for k in nzrange(J, c)
            if rowvals(J)[k] == r
                J.nzval[k] = 0.0
            end
        end
    end
    return J
end


function zero_row!(J::AbstractMatrix, r::Int)
    @views J[r, :] .= zero(eltype(J))
    return J
end

# ------------------------------
# Full Jacobian matrix including PV voltage columns and PV identity rows
# ------------------------------
function calcJacobian_withPVIdentity(
    Y::AbstractMatrix{ComplexF64},
    busVec::Vector{BusData},
    adjBranch::Vector{Vector{Int}},
    slackIdx::Int,
    n_pq::Int, n_pv::Int;
    log::Bool=false, sparse::Bool=true)

    # Remove only the slack bus from θ/V as usual; do not eliminate anything else
    shiftIJ(idx::Int) = (idx >= slackIdx) ? idx - 1 : idx

    n = n_pq + n_pv
    size = 2*n
    J = sparse ? spzeros(Float64, size, size) : zeros(Float64, size, size)

    @inbounds for i in eachindex(busVec)
        # SKIP: never build Jacobian rows for the slack bus
        if busVec[i].type == Sparlectra.Slack
            continue
        end
        
        vm_i = busVec[i].vm_pu
        va_i = busVec[i].va_rad
        type_i = busVec[i].type

        i1 = 2*shiftIJ(i)-1   # P-row i
        i2 = 2*shiftIJ(i)     # Q/constraint row i

        for j in adjBranch[i]
            # SKIP: never build Jacobian columns for the slack bus
            if busVec[j].type == Sparlectra.Slack
                continue
            end
            
            vm_j = busVec[j].vm_pu
            va_j = busVec[j].va_rad

            j1 = 2*shiftIJ(j)-1   # θ-column j
            j2 = 2*shiftIJ(j)     # V-column j (now also for PV buses)

            if i == j
                Yii = abs(Y[i,i]); αii = angle(Y[i,i]); arg = -αii
                # P-θ
                J[i1, j1] = vm_i*(Yii*vm_i*sin(arg)) - vm_i*sum(abs(Y[i,k])*busVec[k].vm_pu *
                               sin(va_i - busVec[k].va_rad - angle(Y[i,k])) for k in adjBranch[i])
                # Q-θ
                J[i2, j1] = -vm_i*(Yii*vm_i*cos(arg)) + vm_i*sum(abs(Y[i,k])*busVec[k].vm_pu *
                               cos(va_i - busVec[k].va_rad - angle(Y[i,k])) for k in adjBranch[i])
                # P-V (relative)
                J[i1, j2] = vm_i*(Yii*vm_i*cos(arg)) + vm_i*sum(abs(Y[i,k])*busVec[k].vm_pu *
                               cos(va_i - busVec[k].va_rad - angle(Y[i,k])) for k in adjBranch[i])
                # Q-V (relative)
                J[i2, j2] = vm_i*(Yii*vm_i*sin(arg)) + vm_i*sum(abs(Y[i,k])*busVec[k].vm_pu *
                               sin(va_i - busVec[k].va_rad - angle(Y[i,k])) for k in adjBranch[i])
            else
                Yij = abs(Y[i,j]); arg = va_i - va_j - angle(Y[i,j])
                # P-θ_j
                J[i1, j1] = vm_i * Yij * vm_j * sin(arg)
                # Q-θ_j
                J[i2, j1] = -vm_i * Yij * vm_j * cos(arg)
                # P-V_j (relative)
                J[i1, j2] = vm_i * Yij * vm_j * cos(arg)
                # Q-V_j (relative)
                J[i2, j2] = vm_i * Yij * vm_j * sin(arg)
            end
        end

        # PV constraint: row i2 as (V - Vset) = 0 => ∂/∂ΔV_rel_i = +V_i
        if type_i == Sparlectra.PV
            if issparse(J)
                zero_row!(J, i2)               # sparse path
            else
                @views J[i2, :] .= 0.0         # dense path
            end
            J[i2, 2*shiftIJ(i)] = vm_i
        end
    end

    if log
        println("calcJacobian_withPVIdentity: size = ", size, " × ", size)
    end
    return J
end

# ------------------------------
# Full Newton-Raphson including PV identity rows (separate from the reduced version)
# ------------------------------
function calcNewtonRaphson_withPVIdentity!(
    net::Net, Y::AbstractMatrix{ComplexF64}, maxIte::Int;
    tolerance::Float64=1e-6, verbose::Int=0, sparse::Bool=false,
    flatStart::Bool=false, angle_limit::Bool=false, debug::Bool=false)

    setJacobianAngleLimit(angle_limit)
    setJacobianDebug(debug)

    # Q-Limits (p.u.) 
    qmin_pu, qmax_pu = get_Q_limits_pu(net)

    nodes = net.nodeVec
    Sbase_MVA = net.baseMVA

    busVec, slackNum = getBusData(nodes, Sbase_MVA, flatStart)
    busTypeVec, slackIdx = getBusTypeVec(busVec)

    n_pv = count(b->b.type==PV, busVec)
    n_pq = count(b->b.type==PQ, busVec)
    n    = n_pq + n_pv

    # Q-Limits vorbereiten (falls Zahl der Busse nicht passt → (re)build)
    if length(net.qmin_pu) != length(busVec) || length(net.qmax_pu) != length(busVec)
        buildQLimits!(net)
    end
    # Logging-Reset am Start
    resetQLimitLog!(net)


    # Vset: actual setpoints for PV, 1.0 (dummy) otherwise
    Vset = [ (b.type==Sparlectra.PV) ?
         (isnothing(net.nodeVec[b.idx]._vm_pu) ? b.vm_pu : net.nodeVec[b.idx]._vm_pu) :
         1.0
       for b in busVec ]


    adjBranch = adjacentBranches(Y, debug)

    it = 0; erg = 1
    feeders = getPowerFeeds_full(busVec, n_pq, n_pv)

    while it <= maxIte
        Δ = residuum_full_withPV(Y, busVec, Vset, n_pq, n_pv, (verbose>2))

        # --- Active-Set: check Q-Limits of PV buses and switch PV->PQ if necessary ---
        changed = false
        @inbounds for k in eachindex(busVec)
            b = busVec[k]
            if b.type == PV
                qreq = b._qRes                      # aktueller Q-Bedarf (p.u.)
                busIdx = b.idx

                # -----------------------------------------------------------
                # [OPTIONAL: Cooldown benutzen, wenn du irgendwann PQ->PV
                #  re-enablen willst. Hier NICHTS tun: wir schalten nur PV->PQ.
                #  Den Cooldown brauchst du beim Rückschalten (siehe unten).]
                # -----------------------------------------------------------

                if qreq > net.qmax_pu[k]
                    # trifft obere Grenze
                    b.type = PQ
                    b.qƩ  = net.qmax_pu[k]
                    changed = true
                    net.qLimitEvents[busIdx] = :max
                    Sparlectra.logQLimitHit!(net, it, busIdx, :max)   # <— strukturierter Log
                    (verbose>0) && @printf "  PV->PQ Bus %d: Q=%.4f > Qmax=%.4f\n" busIdx qreq net.qmax_pu[k]

                elseif qreq < net.qmin_pu[k]
                    # trifft untere Grenze
                    b.type = PQ
                    b.qƩ  = net.qmin_pu[k]
                    changed = true
                    net.qLimitEvents[busIdx] = :min
                    Sparlectra.logQLimitHit!(net, it, busIdx, :min)   # <— strukturierter Log
                    (verbose>0) && @printf "  PV->PQ Bus %d: Q=%.4f < Qmin=%.4f\n" busIdx qreq net.qmin_pu[k]
                end
            end
        end
        if changed
            Δ = residuum_full_withPV(Y, busVec, Vset, n_pq, n_pv, (verbose>2))
        end
        nrm = norm(Δ)
        (verbose>0) && @printf " norm %e, tol %e, ite %d\n" nrm tolerance it

        if nrm < tolerance
            (verbose>0) && println("Convergence after $(it) iterations")
            erg = 0; break
        end

        J = calcJacobian_withPVIdentity(Y, busVec, adjBranch, slackIdx, n_pq, n_pv;
                                        log=(verbose>=3), sparse=sparse)

        # Solve
        Δx = nothing
        try
            Δx = J \ Δ
            if verbose > 3
                println("max |Δ| = ", maximum(abs.(Δx)))
            end
        catch err
            println("Linear solve failed: ", err)
            erg = 2; break
        end

        # Update (PV buses get ΔV_rel from the identity row)
        kSlack = 0
        for i in eachindex(busVec)
            if busVec[i].type == Slack
                kSlack += 1; continue
            end
            # Full index without elimination (only shift for slack)
            idx1 = 2*((i >= slackIdx) ? (i-1) : i) - 1
            idx2 = idx1 + 1

            if angle_limit
                busVec[i].va_rad += min(Δx[idx1], 0.7)
            else
                busVec[i].va_rad += Δx[idx1]
            end
            busVec[i].vm_pu += busVec[i].vm_pu * Δx[idx2]

            # Reprojection (keep polar form consistent)
            Vc = busVec[i].vm_pu * exp(im*busVec[i].va_rad)
            busVec[i].vm_pu = abs(Vc)
            busVec[i].va_rad = angle(Vc)
        end

        it += 1
    end

    # Write results back to the network object
    isoNodes = net.isoNodes
    for bus in busVec
        idx = bus.idx
        if idx in isoNodes; continue; end
        idx += count(i -> i < idx, isoNodes)
        setVmVa!(node=nodes[idx], vm_pu=bus.vm_pu, va_deg=rad2deg(bus.va_rad))
        if bus.type == PV
            nodes[idx]._qƩGen = bus._qRes * Sbase_MVA
        elseif bus.type == Slack
            nodes[idx]._pƩGen = bus._pRes * Sbase_MVA
            nodes[idx]._qƩGen = bus._qRes * Sbase_MVA
        end
    end

    setTotalBusPower!(net=net,
        p=sum(b->b._pRes, busVec), q=sum(b->b._qRes, busVec))


    if !isempty(net.qLimitEvents) && verbose>0
        println("Q-limit events (BusIdx => side): ", net.qLimitEvents)
    end

    return it, erg
end

# ------------------------------
# Convenience wrapper similar to runpf!, but for the full system
# ------------------------------
function runpf_full!(net::Net, maxIte::Int, tolerance::Float64=1e-6, verbose::Int=0)
    sparse = length(net.nodeVec) > 60
    printYBus = (length(net.nodeVec) < 20) && (verbose > 1)
    Y = createYBUS(net=net, sparse=sparse, printYBUS=printYBus)
    return calcNewtonRaphson_withPVIdentity!(net, Y, maxIte;
        tolerance=tolerance, verbose=verbose, sparse=sparse,
        flatStart=false, angle_limit=false, debug=false)
end

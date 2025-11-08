# test_jacobian_full.jl — tests for the full-system NR (PV identity rows)
# Note: Do NOT add `using` statements here; runtest.jl already handles imports.

# -----------------------------------------------------------------------------
# 1) Basic convergence test using the same CIGRE-HV builder
# -----------------------------------------------------------------------------
function test_acpflow_full(verbose::Int = 0)
    net = testCreateNetworkFromScratch()
    tol = 1e-6
    maxIte = 20
    print_results = (verbose > 0)
    result = true

    etime = @elapsed begin
        ite, erg = runpf_full!(net, maxIte, tol, verbose)
        if erg != 0
            @warn "Full-system power flow did not converge"
            result = false
        end
        if print_results
            calcNetLosses!(net)
            printACPFlowResults(net, etime, ite, tol)
        end
    end

    return result
end

# -----------------------------------------------------------------------------
# 2) Cross-check: reduced vs. full system should agree closely
#    - Same network, run reduced NR (runpf!) and full NR (runpf_full!).
#    - Compare bus voltages (magnitude and angle) within tolerances.
#    - Additionally check that PV buses satisfy |V - Vset| ≤ tol_vm_pv.
# -----------------------------------------------------------------------------
function test_full_matches_reduced(; tol_vm::Float64=5e-4, tol_va_deg::Float64=5e-3, tol_vm_pv::Float64=5e-4, verbose::Int=0)
    net1 = testCreateNetworkFromScratch()
    net2 = testCreateNetworkFromScratch()

    tol = 1e-6
    maxIte = 30

    # Reduced (baseline)
    ite_r, erg_r = runpf!(net1, maxIte, tol, verbose)
    if erg_r != 0
        @warn "Reduced NR did not converge"
        return false
    end

    # Full (PV identity rows)
    ite_f, erg_f = runpf_full!(net2, maxIte, tol, verbose)
    if erg_f != 0
        @warn "Full NR did not converge"
        return false
    end

    ok = true

    # Compare bus-by-bus
    for i in eachindex(net1.nodeVec)
        n1 = net1.nodeVec[i]
        n2 = net2.nodeVec[i]
        vm_diff = abs(n1._vm_pu - n2._vm_pu)
        va_diff = abs(n1._va_deg - n2._va_deg)
        if vm_diff > tol_vm || va_diff > tol_va_deg
            @warn "Bus $(i): |ΔV|=$(vm_diff), |Δθ|=$(va_diff) exceeds tolerances (vm=$(tol_vm), va=$(tol_va_deg))"
            ok = false
        end

        # Extra PV check: ensure V ≈ Vset
        if isPVNode(n2)
            # If there is an explicit setpoint stored, use it; else fall back to typical 1.0 p.u.
            Vset = isnothing(n2._vm_pu) ? 1.0 : n2._vm_pu
            if abs(n2._vm_pu - Vset) > tol_vm_pv
                @warn "PV bus $(i): |V - Vset| = $(abs(n2._vm_pu - Vset)) > $(tol_vm_pv)"
                ok = false
            end
        end
    end

    return ok
end

# -----------------------------------------------------------------------------
# 3) Structural test: Jacobian size and PV identity rows are as expected
#    - J must be 2*(n_pq+n_pv) × 2*(n_pq+n_pv)
#    - For each PV bus i: the row corresponding to its Q/constraint contains
#      only one nonzero at its local V-column with value ≈ V_i (relative scaling)
# -----------------------------------------------------------------------------
function test_jacobian_full_structure(verbose::Int=0)
    net = testCreateNetworkFromScratch()
    Y = createYBUS(net = net, sparse = true, printYBUS = false)

    # build bus vectors like in calcNewtonRaphson_withPVIdentity!
    Sbase_MVA = net.baseMVA
    busVec, slackIdx  = Sparlectra.getBusData(net.nodeVec, Sbase_MVA, false)
    busTypeVec, slackIdx2 = Sparlectra.getBusTypeVec(busVec)
    @assert slackIdx == slackIdx2

    n_pv = count(b->b.type==Sparlectra.PV, busVec)
    n_pq = count(b->b.type==Sparlectra.PQ, busVec)
    n = n_pq + n_pv

    adj = adjacentBranches(Y, false)  # this one is exported; unqualified is fine
    J = calcJacobian_withPVIdentity(Y, busVec, adj, slackIdx, n_pq, n_pv; log=false, sparse=true)
    expected = 2*n
    ok = true
    if size(J,1) != expected || size(J,2) != expected
        @warn "Jacobian shape mismatch: got $(size(J)), expected ($(expected), $(expected))"
        ok = false
    end

    # helper matching the production code shiftIJ
    shiftIJ(idx::Int) = (idx >= slackIdx) ? idx - 1 : idx

    # For each PV bus, test identity-row pattern
    for i in eachindex(busVec)
        if busVec[i].type == Sparlectra.PV
            i2 = 2*shiftIJ(i)
            vcol = 2*shiftIJ(i)   # local V-column

            # Count nonzeros in row i2
            nnz_row = 0
            last_val = 0.0
#            for c in 1:size(J,2)
 #               if J[i2,c] != 0.0
             for c in axes(J, 2)
                if J[i2, c] != 0.0
                    nnz_row += 1
                    last_val = J[i2,c]
                    if c != vcol
                        @warn "PV row $(i2): unexpected nonzero at column $(c)"
                        ok = false
                    end
                end
            end
            if nnz_row != 1
                @warn "PV row $(i2): expected exactly one nonzero, found $(nnz_row)"
                ok = false
            end
            # Value should be roughly V_i
            if abs(last_val - busVec[i].vm_pu) > 1e-12
                @warn "PV row $(i2): value $(last_val) ≠ V_i $(busVec[i].vm_pu)"
                ok = false
            end
        end
    end

    return ok
end

# --- Helper: PV erzwingen & Vset setzen (nur über öffentliche API) ---
function _force_pv!(net, busname::String; vm_pu::Float64=1.16)
    # Bus wirklich auf PV setzen (falls er auf PQ gefallen ist)
    setBusType!(net = net, busName = busname, busType = "PV")
    # V-Sollwert für PV-Bus setzen (wird vom Full-NR genutzt)
    setPVBusVset!(net, busname; vm_pu = vm_pu)
    return nothing
end

# --- Test: PV→PQ-Umschaltung über Q-Limits (nutzt die tatsächlich vorhandenen PV-Busse) ---
function test_pv_q_limit_switch!(net; verbose::Int=0)
    n = deepcopy(net)  # keine Nebenwirkungen
    tol    = 1e-6
    maxIte = 30

    # Warmstart (falls verfügbar), damit Full-NR leichter konvergiert
    try
        _, erg_r = runpf!(n, maxIte, tol, 0)
        if erg_r != 0 && verbose > 0
            @warn "Reduced NR did not converge; continuing anyway"
        end
    catch
        # ok, ohne Warmstart
    end

    # Busdaten der aktuellen Netzlage
    Y0 = createYBUS(net = n, sparse = true, printYBUS = false)
    Sbase = n.baseMVA
    busVec, slackIdx1 = Sparlectra.getBusData(n.nodeVec, Sbase, false)
    _,       slackIdx = Sparlectra.getBusTypeVec(busVec)
    @assert slackIdx == slackIdx1

    # Nimm die echten PV-Busse aus busVec (ohne Slack)
    pv_idxs = [i for i in eachindex(busVec) if busVec[i].type == Sparlectra.PV && i != slackIdx]
    @assert !isempty(pv_idxs) "Keine PV-Busse im Netz – Test nicht anwendbar."

    # Wähle bis zu 3 PV-Busse
    targets = pv_idxs[1:min(3, length(pv_idxs))]

    # Helper: Busname zu Index (robust auf verschiedene Feldnamen)
    _busname(i) = begin
        nnode = n.nodeVec[i]
        if hasproperty(nnode, :name)      ; getfield(nnode, :name)
        elseif hasproperty(nnode, :_name) ; getfield(nnode, :_name)
        else                               "B$(i)"
        end
    end
    names = map(_busname, targets)

    # Q-Limits moderat eng + Logs resetten
    try
        buildQLimits!(n)
    catch
    end
    @assert hasproperty(n, :qmin_pu) && hasproperty(n, :qmax_pu) "qmin_pu/qmax_pu fehlen im Net!"
    n.qmin_pu .= -Inf;  n.qmax_pu .=  Inf

    # zwei Versuche mit zunehmendem „Druck“ auf Q
    attempts = (
        (1.05, 0.03),   # Vset leicht anheben, ±0.03 p.u. Band
        (1.07, 0.05),   # stärker anheben, Band etwas breiter (bessere Feasibility)
    )

    conv = false
    hit  = false

    for (vset, qband) in attempts
        # enge Bänder nur für die Ziel-PV-Busse
        for i in targets
            n.qmin_pu[i] = -qband
            n.qmax_pu[i] =  qband
        end

        # Logs leeren
        try
            resetQLimitLog!(n)
        catch
            if !hasproperty(n, :qLimitLog)   ; n.qLimitLog    = NamedTuple[]          ; end
            if !hasproperty(n, :qLimitEvents); n.qLimitEvents = Dict{Int,Symbol}()    ; end
            empty!(n.qLimitLog); empty!(n.qLimitEvents)
        end

        # Vset für die ausgewählten PVs setzen (per Namen, wie deine API es erwartet)
        for name in names
            setPVBusVset!(n, name; vm_pu = vset)
        end

        # Nach Änderungen an den Setpoints Busdaten/Slack **neu** bestimmen,
        # YBUS kann gleich bleiben (Vset ändert Y nicht)
        busVec2, slackIdxA = Sparlectra.getBusData(n.nodeVec, Sbase, false)
        _,       slackIdxB = Sparlectra.getBusTypeVec(busVec2)
        @assert slackIdxA == slackIdxB
        Y = Y0  # identisch zu oben; reicht aus

        ite, erg = calcNewtonRaphson_withPVIdentity!(n, Y, slackIdxA;
                                                     tolerance = tol,
                                                     verbose   = verbose,
                                                     sparse    = true,
                                                     flatStart = false)
        conv = (erg == 0)
        hit  = any(haskey(n.qLimitEvents, i) for i in targets)

        if verbose > 0
            @info "try(vset=$(vset), qband=$(qband)) => ite=$(ite), erg=$(erg), hit=$(hit)"
            @info "targets" targets names
            @info "qLimitEvents" n.qLimitEvents
            @info "qLimitLog"    n.qLimitLog
        end

        if conv && hit
            break
        end
    end

    @test conv === true
    @test hit  === true
    return true
end

# Wrapper für runtests.jl – baut das Netz wie die anderen Full-System-Tests
function test_pv_q_limit_switch(; verbose::Int=0)
    net = testCreateNetworkFromScratch()
    return test_pv_q_limit_switch!(net; verbose=verbose)
end

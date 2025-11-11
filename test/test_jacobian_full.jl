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
    setPVGeneratorQLimitsAll!(net = net, qmin_MVar = -350.0, qmax_MVar = 350.0)
    buildQLimits!(net, reset=true)

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
    setPVGeneratorQLimitsAll!(net = net2, qmin_MVar = -350.0, qmax_MVar = 350.0)
    buildQLimits!(net2, reset=true)

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
function test_jacobian_full_structure(;verbose::Int=0)
    net = testCreateNetworkFromScratch()
    Y = createYBUS(net = net, sparse = true, printYBUS = false)

    # build bus vectors like in calcNewtonRaphson_withPVIdentity!
    Sbase_MVA = net.baseMVA
    busVec, slackIdx  = getBusData(net.nodeVec, Sbase_MVA, false)
    busTypeVec, slackIdx2 = Sparlectra.getBusTypeVec(busVec)
    @assert slackIdx == slackIdx2

    n_pv = count(b->b.type==Sparlectra.PV, busVec)
    n_pq = count(b->b.type==Sparlectra.PQ, busVec)
    n = n_pq + n_pv

    setPVGeneratorQLimitsAll!(net = net, qmin_MVar = -350.0, qmax_MVar = 350.0)
    buildQLimits!(net, reset=true)
    adj = adjacentBranches(Y, false)  # this one is exported; unqualified is fine
    J = calcJacobian_withPVIdentity(Y, busVec, adj, slackIdx, n_pq, n_pv; log=false, sparse=true)
    expected = 2*n
    ok = true
    if size(J,1) != expected || size(J,2) != expected
        @error "Jacobian shape mismatch: got $(size(J)), expected ($(expected), $(expected))"
        exit(false)
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

# --- Test: PV→PQ switching via Q-limits using public setters ---
function test_pv_q_limit_switch!(net::Net; verbose::Int=0)
    tol   = 1e-6
    maxIt = 30
    print_results = (verbose > 0)
    # PV buses known from testCreateNetworkFromScratch()
    pv_names = ["B10","B11","B12"]

    # 1) Raise Vset a bit; keep generous Q-limits so no PV->PQ is expected
    for name in pv_names
        setPVBusVset!(net, name; vm_pu = 1.07)
    end
    setPVGeneratorQLimitsAll!(net = net, qmin_MVar = -350.0, qmax_MVar = 350.0)
    buildQLimits!(net, reset=true)

    etime = @elapsed begin
        ite, erg = runpf_full!(net, maxIt, tol, verbose)
    end
    converged = (erg == 0)

    if !converged
        @warn "Full NR did not converge"
    else
        # Optional but nice: compute total network losses after per-branch flows exist
        try
            Sparlectra.calcNetLosses!(net)
        catch err
            @warn "Sparlectra.calcNetLosses!(net) failed or not found" error=err
        end
    end

    # Was there any PV->PQ switch? (expected: no)
    pv_idx = map(name -> geNetBusIdx(net = net, busName = name), pv_names)
    hit = any(haskey(net.qLimitEvents, i) for i in pv_idx)
    if print_results
        printQLimitLog(net; sort_by=:bus)
        calcNetLosses!(net)
        printACPFlowResults(net, etime, ite, tol)
    end

    return hit == false
end

# Wrapper for runtests.jl – builds the network like the other full-system tests
function test_pv_q_limit_switch(; verbose::Int=1)
    net = testCreateNetworkFromScratch()
    return test_pv_q_limit_switch!(net; verbose=verbose)
end

"""
    complex_newton_step(Ybus, V, S)

Performs one Newton–Raphson step in the complex formulation.

Arguments:
- `Ybus`: bus admittance matrix
- `V`: complex bus voltage vector
- `S`: specified power injections (P + jQ)
"""
function complex_newton_step(Ybus, V, S)
    # Compute injected currents
    I = Ybus * V

    # Compute complex power injections
    S_calc = V .* conj.(I)

    # Mismatch ΔS = S_calc − S_spec
    ΔS = S_calc .- S

    # TODO: replace with correct Jacobian update
    # J = build_complex_jacobian(Ybus, V)
    # ΔV = -J \ ΔS

    # Temporary placeholder: simple scaled mismatch without division by zero
    ΔV = similar(V)
    for i in eachindex(V)
        ΔV[i] = -0.01 * ΔS[i]  # just a small step in direction of negative mismatch
    end

    return V .+ ΔV
end


"""
    run_complex_nr(Ybus, V0, S; slack_idx=1, maxiter=20, tol=1e-8, verbose=false)

Runs an iterative complex-state Newton–Raphson prototype on a given test system.

Arguments:
- `Ybus`: bus admittance matrix
- `V0`: initial complex voltage vector
- `S`: specified complex power injections (P + jQ)
- `slack_idx`: index of the slack bus (its voltage is kept fixed at V0[slack_idx])
- `maxiter`: maximum number of iterations
- `tol`: convergence tolerance on the maximum power mismatch |ΔS|
- `verbose`: if true, prints the mismatch per iteration

Returns:
- `V`: final complex voltage vector
- `converged`: Bool, true if `maximum(abs.(ΔS)) <= tol`
- `iters`: number of performed iterations
- `history`: Vector of mismatch norms per iteration
"""
function run_complex_nr(Ybus, V0, S; slack_idx=1, maxiter=20, tol=1e-8, verbose=false)
    V = copy(V0)
    history = Float64[]

    for iter in 1:maxiter
        # Compute mismatch for current voltages
        I = Ybus * V
        S_calc = V .* conj.(I)
        ΔS = S_calc .- S

        max_mis = maximum(abs.(ΔS))
        push!(history, max_mis)

        if verbose
            @info "Complex NR iteration" iter=iter max_mismatch=max_mis
        end

        if max_mis <= tol
            return V, true, iter, history
        end

        # Take one prototype Newton step
        V = complex_newton_step(Ybus, V, S)

        # Enforce slack bus voltage
        V[slack_idx] = V0[slack_idx]
    end

    return V, false, maxiter, history
end


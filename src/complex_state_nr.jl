"""
    complex_newton_step(Ybus, V, S)

Performs one Newton–Raphson step in the complex formulation using
a Wirtinger-type Jacobian on the complex power injections.

Arguments:
- `Ybus`: bus admittance matrix
- `V`: current complex bus voltage vector
- `S`: specified power injections (P + jQ)

Returns:
- updated complex voltage vector `V_new`
"""
function complex_newton_step(Ybus, V, S)
    # Compute currents and complex power
    I = Ybus * V
    S_calc = V .* conj.(I)

    # Complex power mismatch ΔS = S_calc − S_spec
    ΔS = S_calc .- S

    # Build Jacobian blocks
    J11, J12, J21, J22 = build_complex_jacobian(Ybus, V)

    n = length(V)

    # Assemble full 2n × 2n Jacobian
    J_top = hcat(J11, J12)
    J_bot = hcat(J21, J22)
    J = vcat(J_top, J_bot)

    # Right-hand side: -[ΔS; conj(ΔS)]
    rhs = vcat(-ΔS, -conj.(ΔS))

    # Solve for [ΔV; ΔV*] in a robust way
    sol = nothing
    try
        sol = J \ rhs
    catch e
        if e isa LinearAlgebra.SingularException
            # Fallback: use pseudoinverse (minimum-norm solution)
            sol = pinv(J) * rhs
        else
            rethrow(e)
        end
    end

    ΔV = sol[1:n]  # upper block is the actual voltage correction

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


"""
    complex_newton_step_rectangular(Ybus, V, S; damp=1.0)

Performs one Newton–Raphson step in rectangular coordinates using
the complex bus voltages V = Vr + j*Vi as state. The mismatch is

    F(Vr, Vi) = [ real(S_calc - S_spec);
                  imag(S_calc - S_spec) ]

and the Jacobian is computed w.r.t. [Vr; Vi].

Arguments:
- `Ybus`: bus admittance matrix (n×n, Complex)
- `V`: current complex bus voltage vector (length n)
- `S`: specified complex power injections P + jQ (length n)
- `damp`: scalar damping factor for the Newton step (0 < damp ≤ 1)

Returns:
- updated complex voltage vector `V_new`
"""
function complex_newton_step_rectangular(Ybus, V, S; damp=1.0)
    n = length(V)

    # Currents and complex power
    I = Ybus * V
    S_calc = V .* conj.(I)
    ΔS = S_calc .- S

    # Real state representation
    Vr = real.(V)
    Vi = imag.(V)

    # Complex Jacobian Jc: n × 2n, mapping [dVr; dVi] → dS (complex)
    Jc = Matrix{ComplexF64}(undef, n, 2n)

    # Derivative wrt Vr: perturb dVr = e_k, dVi = 0
    for k in 1:n
        dVr = zeros(Float64, n)
        dVi = zeros(Float64, n)
        dVr[k] = 1.0
        dV = ComplexF64.(dVr, dVi)

        dI = Ybus * dV
        dS = dV .* conj.(I) .+ V .* conj.(dI)

        Jc[:, k] = dS
    end

    # Derivative wrt Vi: perturb dVi = e_k, dVr = 0
    for k in 1:n
        dVr = zeros(Float64, n)
        dVi = zeros(Float64, n)
        dVi[k] = 1.0
        dV = ComplexF64.(dVr, dVi)

        dI = Ybus * dV
        dS = dV .* conj.(I) .+ V .* conj.(dI)

        Jc[:, n + k] = dS
    end

    # Build real 2n × 2n Jacobian for F = [real(ΔS); imag(ΔS)]
    J = zeros(Float64, 2n, 2n)
    for k in 1:2n
        col = Jc[:, k]
        J[1:n, k] .= real.(col)
        J[n+1:2n, k] .= imag.(col)
    end

    # Real mismatch vector
    F = vcat(real.(ΔS), imag.(ΔS))

    # Solve J * δx = -F for δx = [δVr; δVi] (robust against singularity)
    δx = nothing
    try
        δx = J \ (-F)
    catch e
        if e isa LinearAlgebra.SingularException
            δx = pinv(J) * (-F)
        else
            rethrow(e)
        end
    end

    # Apply damping
    δx .*= damp

    δVr = δx[1:n]
    δVi = δx[n+1:2n]

    V_new = ComplexF64.(Vr .+ δVr, Vi .+ δVi)
    return V_new
end

"""
    run_complex_nr_rectangular(Ybus, V0, S; slack_idx=1, maxiter=20, tol=1e-8, verbose=false)

Runs an iterative Newton–Raphson power flow in rectangular coordinates
using complex bus voltages as state and the function
    complex_newton_step_rectangular.

Arguments:
- `Ybus`: bus admittance matrix
- `V0`: initial complex voltages
- `S`: specified power injections (P + jQ)
- `slack_idx`: index of the slack bus (its voltage is kept fixed at V0[slack_idx])
- `maxiter`: maximum number of iterations
- `tol`: convergence tolerance on max |ΔS|
- `verbose`: print mismatch per iteration if true

Returns:
- `V`: final complex voltages
- `converged`: Bool
- `iters`: performed iterations
- `history`: Vector of mismatch norms per iteration
"""
function run_complex_nr_rectangular(Ybus, V0, S;
                                    slack_idx=1,
                                    maxiter=20,
                                    tol=1e-8,
                                    verbose=false,
                                    damp=0.2)
    V = copy(V0)
    history = Float64[]

    for iter in 1:maxiter
        I = Ybus * V
        S_calc = V .* conj.(I)
        ΔS = S_calc .- S

        max_mis = maximum(abs.(ΔS))
        push!(history, max_mis)

        if verbose
            @info "Rectangular NR iteration" iter=iter max_mismatch=max_mis
        end

        if max_mis <= tol
            return V, true, iter, history
        end

        # Newton step in rectangular coordinates
        V = complex_newton_step_rectangular(Ybus, V, S; damp=damp)

        # Enforce slack bus voltage
        V[slack_idx] = V0[slack_idx]
    end

    return V, false, maxiter, history
end

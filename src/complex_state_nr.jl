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

function complex_newton_step_rectangular(Ybus, V, S; slack_idx::Int, damp::Float64=1.0)
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
        J[1:n, k]       .= real.(col)
        J[n+1:2n, k]    .= imag.(col)
    end

    # Real mismatch vector F = [Re(ΔS); Im(ΔS)]
    F = vcat(real.(ΔS), imag.(ΔS))

    # --- Slack-Elimination: nur Nicht-Slack-Busse als Unbekannte/ Gleichungen ---

    # Bus-Indizes ohne Slack
    non_slack = collect(1:n)
    deleteat!(non_slack, slack_idx)

    # Zeilen/Gleichungen: P/Q für alle Nicht-Slack-Busse
    row_idx = vcat(non_slack, n .+ non_slack)  # erst P, dann Q

    # Spalten/Variablen: Vr/Vi für alle Nicht-Slack-Busse
    col_idx = vcat(non_slack, n .+ non_slack)

    Jred = J[row_idx, col_idx]
    Fred = F[row_idx]

    # Solve Jred * δx_red = -Fred (robust gegen Singularität)
    δx_red = nothing
    try
        δx_red = Jred \ (-Fred)
    catch e
        if e isa LinearAlgebra.SingularException
            δx_red = pinv(Jred) * (-Fred)
        else
            rethrow(e)
        end
    end

    # Auf vollen Zustandsvektor zurückheben (Slack hat Δ=0)
    δx = zeros(Float64, 2n)
    δx[col_idx] .= δx_red

    # Dämpfung
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
                                    slack_idx,
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
        V = complex_newton_step_rectangular_fd(Ybus, V, S; slack_idx=slack_idx, damp=damp)

        # Enforce slack bus voltage
        V[slack_idx] = V0[slack_idx]
    end

    return V, false, maxiter, history
end


"""
    run_complex_nr_rectangular_for_net!(net; maxiter=20, tol=1e-8, damp=0.2, verbose=false)

Run a complex-state Newton–Raphson power flow in rectangular coordinates
on the given `net::Net` object.

The function:
- builds Ybus, the complex voltage start vector V0 and the power injections S
  from the network data,
- runs `run_complex_nr_rectangular` on this data,
- writes the resulting voltages back into the network (vm_pu / va_deg),
- and returns (iters, erg), similar to `runpf!`.

Returns:
- `iters::Int`: number of iterations performed
- `erg::Int`: 0 if converged, nonzero otherwise
"""
function run_complex_nr_rectangular_for_net!(net::Net;
                                             maxiter::Int = 20,
                                             tol::Float64 = 1e-8,
                                             damp::Float64 = 0.2,
                                             verbose::Int = 0)
    
    sparse    = length(net.nodeVec) > 60
    Ybus = createYBUS(net=net, sparse=sparse, printYBUS = (verbose > 1))    
    

    # 2) Build initial complex voltages V0 from bus data
    V0, slack_idx = initial_Vrect_from_net(net)

    # 3) Build specified complex power injections S = P + jQ from bus data
    S = build_S_from_net(net)
    
    if verbose > 0
        @info "Starting rectangular complex NR power flow..."
        @info "Initial complex voltages V0:" V0
        @info "Specified complex power injections S (pu):" S
        @info "Slack bus index:" slack_idx
        @info "max iter:" maxiter "tolerance:" tol "damping:" damp
    end

    # 4) Run rectangular NR
    V, converged, iters, history = run_complex_nr_rectangular(
        Ybus, V0, S;
        slack_idx = slack_idx,
        maxiter   = maxiter,
        tol       = tol,
        verbose   = (verbose>0),
        damp      = damp,
    )

    # 5) Write back voltages to network (vm_pu, va_deg)
    update_net_voltages_from_complex!(net, V)

    erg = converged ? 0 : 1
    return iters, erg
end


"""
    initial_Vrect_from_net(net) -> (V0, slack_idx)

Build the initial complex voltage vector V0 from the network bus data
(Vm, Va), and detect the slack bus index.

Returns:
- V0::Vector{ComplexF64}
- slack_idx::Int
"""
function initial_Vrect_from_net(net::Net)
    nodes = net.nodeVec  
    n = length(nodes)
    V0 = Vector{ComplexF64}(undef, n)

    slack_idx = 0

    for (k, node) in enumerate(nodes)
        
        vm = node._vm_pu
        va_deg = node._va_deg
        va_rad = deg2rad(va_deg)

        V0[k] = vm * cis(va_rad)

        if isSlack(node)
            @assert slack_idx == 0 "Multiple slack buses detected"
            slack_idx = k
        end
    end

    if slack_idx == 0
        error("No slack bus found in network")
    end

    return V0, slack_idx
end

"""
    build_S_from_net(net) -> S::Vector{ComplexF64}

Build the specified complex power injection vector S = P + jQ in per-unit
for each bus, based on the net's bus load / generation / shunt data.
Positive P/Q means net injection into the bus (generation),
negative means net consumption (load).

"""
function build_S_from_net(net::Net)
    nodes = net.nodeVec
    n = length(nodes)
    S = Vector{ComplexF64}(undef, n)

    baseMVA = net.baseMVA  # laut Doku vorhanden :contentReference[oaicite:3]{index=3}

    for (k, node) in enumerate(nodes)
        Pgen_MW   = isnothing(node._pƩGen) ? 0.0 : node._pƩGen
        Qgen_MVar = isnothing(node._qƩGen) ? 0.0 : node._qƩGen
        Pload_MW  = isnothing(node._pƩLoad) ? 0.0 : node._pƩLoad
        Qload_MVar= isnothing(node._qƩLoad) ? 0.0 : node._qƩLoad
        Psh_MW    = isnothing(node._pShunt) ? 0.0 : node._pShunt
        Qsh_MVar  = isnothing(node._qShunt) ? 0.0 : node._qShunt

        Pinj_MW   = Pgen_MW - Pload_MW - Psh_MW
        Qinj_MVar = Qgen_MVar - Qload_MVar - Qsh_MVar

        Ppu = Pinj_MW   / baseMVA
        Qpu = Qinj_MVar / baseMVA

        S[k] = ComplexF64(Ppu, Qpu)
    end

    return S
end


"""
    update_net_voltages_from_complex!(net, V)

Update the bus voltage magnitudes and angles in the network from the
final complex voltages V (in per-unit).
"""
function update_net_voltages_from_complex!(net::Net, V::Vector{ComplexF64})
    nodes = net.nodeVec
    n = length(nodes)
    @assert length(V) == n

    for (k, node) in enumerate(nodes)
        Vk = V[k]
        vm = abs(Vk)
        va_rad = angle(Vk)
        va_deg = rad2deg(va_rad)

        # TODO: adjust field names to your Node type
        node._vm_pu = vm
        node._va_deg = va_deg
    end

    return nothing
end

"""
    mismatch_rectangular(Ybus, V, S, slack_idx) -> F::Vector{Float64}

Compute the real-valued mismatch vector F(V) for the rectangular
complex-state formulation:

    F(V) = [ real(ΔS_non_slack);
             imag(ΔS_non_slack) ]

where ΔS = S_calc(V) - S_spec, and entries corresponding to the
slack bus are removed (no equations at the slack bus).
"""
function mismatch_rectangular(Ybus, V, S, slack_idx::Int)
    n = length(V)

    I = Ybus * V
    S_calc = V .* conj.(I)
    ΔS = S_calc .- S

    # remove slack bus from equations
    non_slack = collect(1:n)
    deleteat!(non_slack, slack_idx)

    ΔS_ns = ΔS[non_slack]

    F = vcat(real.(ΔS_ns), imag.(ΔS_ns))
    return F
end

"""
    complex_newton_step_rectangular_fd(Ybus, V, S; slack_idx=1, damp=1.0, h=1e-6)

Performs one Newton–Raphson step in rectangular coordinates using a
finite-difference Jacobian on the mismatch

    F(V) = [ real(ΔS_non_slack); imag(ΔS_non_slack) ]

with ΔS = S_calc(V) - S_spec.

Arguments:
- `Ybus`: bus admittance matrix (n×n, Complex)
- `V`: current complex bus voltage vector (length n)
- `S`: specified complex power injections P + jQ (length n)
- `slack_idx`: index of the slack bus
- `damp`: scalar damping factor for the Newton step (0 < damp ≤ 1)
- `h`: perturbation step for finite differences

Returns:
- updated complex voltage vector `V_new`
"""
function complex_newton_step_rectangular_fd(Ybus, V, S;
                                            slack_idx::Int=1,
                                            damp::Float64=1.0,
                                            h::Float64=1e-6)
    n = length(V)

    # --- Basis-Residuum F(V)
    F0 = mismatch_rectangular(Ybus, V, S, slack_idx)
    m = length(F0)  # = 2 * (n-1)

    # Nicht-Slack-Busse
    non_slack = collect(1:n)
    deleteat!(non_slack, slack_idx)

    # Variablen: Vr(non_slack) und Vi(non_slack) => 2*(n-1) Unbekannte
    nvar = 2 * (n - 1)

    J = zeros(Float64, m, nvar)

    # Hilfsfunktion: baue V aus Vr/Vi-Änderung
    Vr = real.(V)
    Vi = imag.(V)

    # Spaltenweise finite Differenzen
    # Reihenfolge der Variablen:
    #   1..(n-1):   Vr(non_slack[k])
    #   (n)..2(n-1): Vi(non_slack[k])
    for (col_idx, bus) in enumerate(non_slack)
        # --- Perturbation in Vr(bus)
        V_pert = copy(V)
        V_pert[bus] = ComplexF64(Vr[bus] + h, Vi[bus])

        Fp = mismatch_rectangular(Ybus, V_pert, S, slack_idx)
        J[:, col_idx] .= (Fp .- F0) ./ h
    end

    for (offset, bus) in enumerate(non_slack)
        col_idx = (n - 1) + offset
        # --- Perturbation in Vi(bus)
        V_pert = copy(V)
        V_pert[bus] = ComplexF64(Vr[bus], Vi[bus] + h)

        Fp = mismatch_rectangular(Ybus, V_pert, S, slack_idx)
        J[:, col_idx] .= (Fp .- F0) ./ h
    end

    # --- Newton-Gleichung: J * δx = -F0
    δx = nothing
    try
        δx = J \ (-F0)
    catch e
        if e isa LinearAlgebra.SingularException
            δx = pinv(J) * (-F0)
        else
            rethrow(e)
        end
    end

    # Dämpfung
    δx .*= damp

    # Auf VR/VI der Nicht-Slack-Busse verteilen
    Vr_new = copy(Vr)
    Vi_new = copy(Vi)

    for (idx, bus) in enumerate(non_slack)
        Vr_new[bus] += δx[idx]
        Vi_new[bus] += δx[(n - 1) + idx]
    end

    # Slack bleibt unverändert
    Vr_new[slack_idx] = Vr[slack_idx]
    Vi_new[slack_idx] = Vi[slack_idx]

    V_new = ComplexF64.(Vr_new, Vi_new)
    return V_new
end

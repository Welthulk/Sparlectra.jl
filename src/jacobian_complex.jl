# jacobian_complex.jl - Complex-State Newton–Raphson Formulation
"""
    build_complex_jacobian(Ybus, V)

Builds the 2n × 2n Wirtinger-type Jacobian blocks for the complex-state
Newton–Raphson formulation.

Given:
    I = Ybus * V
    S = V .* conj.(I)

We construct the blocks:
    J11 = ∂S/∂V
    J12 = ∂S/∂V*
    J21 = ∂conj(S)/∂V
    J22 = ∂conj(S)/∂V*

Returns:
    J11, J12, J21, J22  (all full matrices, not Diagonal)
"""
function build_complex_jacobian(Ybus, V)
    I = Ybus * V
    n = length(V)

    # J11 = diag(conj(I))
    J11 = Matrix(Diagonal(conj.(I)))

    # J12 = diag(V) * conj(Ybus)
    J12 = Matrix(Diagonal(V)) * conj.(Ybus)

    # J21 = conj(J12)
    J21 = conj.(J12)

    # J22 = conj(J11)
    J22 = conj.(J11)

    return J11, J12, J21, J22
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

function run_complex_nr_rectangular(Ybus, V0, S;
                                    slack_idx::Int=1,
                                    maxiter::Int=20,
                                    tol::Float64=1e-8,
                                    verbose::Bool=false,
                                    damp::Float64=0.2,
                                    bus_types::Vector{Symbol},
                                    Vset::Vector{Float64},
                                    use_fd::Bool=true)
    V = copy(V0)
    history = Float64[]

    n = length(V)
    non_slack = collect(1:n)
    deleteat!(non_slack, slack_idx)

    for iter in 1:maxiter
        F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
        max_mis = maximum(abs.(F))
        push!(history, max_mis)

        if verbose
            @info "Rectangular NR iteration" iter=iter max_mismatch=max_mis
        end

        if max_mis <= tol
            return V, true, iter, history
        end

        if use_fd
            V = complex_newton_step_rectangular_fd(
                Ybus, V, S;
                slack_idx = slack_idx,
                damp      = damp,
                h         = 1e-6,
                bus_types = bus_types,
                Vset      = Vset,
            )
        else
            V = complex_newton_step_rectangular(
                Ybus, V, S;
                slack_idx = slack_idx,
                damp      = damp,
                bus_types = bus_types,
                Vset      = Vset,
            )
        end

        V[slack_idx] = V0[slack_idx]
    end

    return V, false, maxiter, history
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
    mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx) -> F::Vector{Float64}

Compute the real-valued mismatch vector F(V) for the rectangular
complex-state formulation with PQ and PV buses.

For each non-slack bus i:
- if bus_types[i] == :PQ:
      ΔP_i = Re(S_calc[i]) - Re(S_spec[i])
      ΔQ_i = Im(S_calc[i]) - Im(S_spec[i])

- if bus_types[i] == :PV:
      ΔP_i = Re(S_calc[i]) - Re(S_spec[i])
      ΔV_i = |V[i]| - Vset[i]

F is stacked as [ΔP_2, ΔQ/ΔV_2, ..., ΔP_n, ΔQ/ΔV_n] over all non-slack buses.
"""
function mismatch_rectangular(Ybus,
                              V::Vector{ComplexF64},
                              S::Vector{ComplexF64},
                              bus_types::Vector{Symbol},
                              Vset::Vector{Float64},
                              slack_idx::Int)

    n = length(V)
    @assert length(S)         == n
    @assert length(bus_types) == n
    @assert length(Vset)      == n

    # Network-based injections for the current state
    I      = Ybus * V
    S_calc = V .* conj.(I)

    # F has 2*(n-1) entries: for each non-slack bus two residuals
    # PQ:  ΔP_i, ΔQ_i
    # PV:  ΔP_i, ΔV_i
    F = zeros(Float64, 2*(n-1))

    row = 1
    @inbounds for i in 1:n
        if i == slack_idx
            continue
        end

        S_ci = S_calc[i]
        S_si = S[i]

        if bus_types[i] == :PQ
            ΔP = real(S_ci) - real(S_si)
            ΔQ = imag(S_ci) - imag(S_si)

            F[row]   = ΔP
            F[row+1] = ΔQ

        elseif bus_types[i] == :PV
            ΔP = real(S_ci) - real(S_si)
            ΔV = abs(V[i]) - Vset[i]

            F[row]   = ΔP
            F[row+1] = ΔV

        else
            error("mismatch_rectangular: unsupported bus type $(bus_types[i]) at bus $i")
        end

        row += 2
    end

    return F
end

"""
    complex_newton_step_rectangular_fd(Ybus, V, S;
                                       slack_idx=1,
                                       damp=1.0,
                                       h=1e-6)

Performs one Newton–Raphson step in rectangular coordinates using a
finite-difference Jacobian on the mismatch

    F(V) = [ real(ΔS_non_slack); imag(ΔS_non_slack) ]

with ΔS = S_calc(V) - S_spec.

Arguments:
- `Ybus`: bus admittance matrix (n×n, Complex)
- `V`: current complex bus voltage vector (length n)
- `S`: specified complex power injections P + jQ (length n)
- `slack_idx`: index of the slack bus (no equations / no variables)
- `damp`: scalar damping factor for the Newton step (0 < damp ≤ 1)
- `h`: perturbation step for finite differences

Returns:
- updated complex voltage vector `V_new`
"""
function complex_newton_step_rectangular_fd(Ybus, V, S;
                                            slack_idx::Int=1,
                                            damp::Float64=1.0,
                                            h::Float64=1e-6,
                                            bus_types::Vector{Symbol},
                                            Vset::Vector{Float64})
    n = length(V)

    # Base mismatch F(V)
    F0 = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
    m = length(F0)  # = 2 * (n-1)

    # Non-slack buses
    non_slack = collect(1:n)
    deleteat!(non_slack, slack_idx)

    # Variables: Vr(non_slack) and Vi(non_slack)
    nvar = 2 * (n - 1)
    @assert nvar == m "Rectangular FD-Newton: nvar and m should both equal 2*(n-1)"

    J = zeros(Float64, m, nvar)

    Vr = real.(V)
    Vi = imag.(V)

    # Columns 1..(n-1): perturb Vr(non_slack[k])
    for (col_idx, bus) in enumerate(non_slack)
        V_pert = copy(V)
        V_pert[bus] = ComplexF64(Vr[bus] + h, Vi[bus])

        Fp = mismatch_rectangular(Ybus, V_pert, S, bus_types, Vset, slack_idx)
        J[:, col_idx] .= (Fp .- F0) ./ h
    end

    # Columns (n)..2(n-1): perturb Vi(non_slack[k])
    for (offset, bus) in enumerate(non_slack)
        col_idx = (n - 1) + offset
        V_pert = copy(V)
        V_pert[bus] = ComplexF64(Vr[bus], Vi[bus] + h)

        Fp = mismatch_rectangular(Ybus, V_pert, S, bus_types, Vset, slack_idx)
        J[:, col_idx] .= (Fp .- F0) ./ h
    end

    # Solve J * δx = -F0
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

    # Damping
    δx .*= damp

    Vr_new = copy(Vr)
    Vi_new = copy(Vi)

    # Apply update to non-slack buses
    for (idx, bus) in enumerate(non_slack)
        Vr_new[bus] += δx[idx]
        Vi_new[bus] += δx[(n - 1) + idx]
    end

    # Keep slack bus fixed
    Vr_new[slack_idx] = Vr[slack_idx]
    Vi_new[slack_idx] = Vi[slack_idx]

    V_new = ComplexF64.(Vr_new, Vi_new)
    return V_new
end

"""
    build_rectangular_jacobian_pq_pv(
        Ybus, V, bus_types, Vset, slack_idx
    ) -> J::Matrix{Float64}

Build the analytic rectangular Jacobian for the mismatch vector `F(V)`
defined in `mismatch_rectangular`.

- State vector: x = [Vr(non-slack); Vi(non-slack)]
- Rows: for each non-slack bus i
    * PQ: [ΔP_i; ΔQ_i]
    * PV: [ΔP_i; ΔV_i]  with ΔV_i = |V_i| - Vset[i]

`bus_types` and `Vset` must be consistent with `mismatch_rectangular`.
"""
function build_rectangular_jacobian_pq_pv(
    Ybus,
    V::Vector{ComplexF64},
    bus_types::Vector{Symbol},
    Vset::Vector{Float64},
    slack_idx::Int
)
    n = length(V)
    @assert length(bus_types) == n
    @assert length(Vset)      == n

    # --- 1) Wirtinger blocks for S(V) = V .* conj(Ybus * V)
    J11, J12, J21, J22 = build_complex_jacobian(Ybus, V)

    # --- 2) Full 2n×2n rectangular J for ΔP/ΔQ wrt Vr/Vi (all buses)
    # Rows: [ΔP_1..ΔP_n; ΔQ_1..ΔQ_n]
    # Cols: [Vr_1..Vr_n; Vi_1..Vi_n]
    Jrect_full = zeros(Float64, 2n, 2n)

    @inbounds for j in 1:n
        col_sum  = J11[:, j] .+ J12[:, j]  # corresponds to dS/dVr_j
        col_diff = J11[:, j] .- J12[:, j]  # used for dS/dVi_j

        # dP/dVr_j, dQ/dVr_j
        @views Jrect_full[1:n,        j] .= real.(col_sum)
        @views Jrect_full[n+1:2n,     j] .= imag.(col_sum)

        # dP/dVi_j, dQ/dVi_j
        @views Jrect_full[1:n,    n + j] .= -imag.(col_diff)
        @views Jrect_full[n+1:2n, n + j] .=  real.(col_diff)
    end

    # --- 3) Reduce to non-slack variables and rows matching mismatch_rectangular

    non_slack = collect(1:n)
    deleteat!(non_slack, slack_idx)

    nvar = 2 * (n - 1)
    m    = 2 * (n - 1)
    @assert nvar == m

    # Column indices in the full rectangular J that correspond to
    # [Vr(non_slack); Vi(non_slack)]
    col_idx_full = vcat(non_slack, n .+ non_slack)

    J = zeros(Float64, m, nvar)

    row = 1
    @inbounds for i in 1:n
        if i == slack_idx
            continue
        end

        # First row for this bus: ΔP_i
        rowP_full = i                  # P row index in full J
        @views J[row, :] .= Jrect_full[rowP_full, col_idx_full]

        # Second row: ΔQ_i (PQ) or ΔV_i (PV)
        if bus_types[i] == :PQ
            rowQ_full = n + i          # Q row index in full J
            @views J[row + 1, :] .= Jrect_full[rowQ_full, col_idx_full]

        elseif bus_types[i] == :PV
            # ΔV_i = |V_i| - Vset[i]
            J[row + 1, :] .= 0.0

            pos = findfirst(==(i), non_slack)
            if pos !== nothing
                vm = abs(V[i])
                if vm > 0.0
                    dVr = real(V[i]) / vm
                    dVi = imag(V[i]) / vm

                    # Columns in reduced J:
                    #   Vr_i -> index pos
                    #   Vi_i -> index (n-1) + pos
                    J[row + 1, pos]           = dVr
                    J[row + 1, (n - 1) + pos] = dVi
                end
            end
        else
            error("build_rectangular_jacobian_pq_pv: unsupported bus type $(bus_types[i]) at bus $i")
        end

        row += 2
    end

    return J
end

"""
    complex_newton_step_rectangular(
        Ybus, V, S;
        slack_idx=1,
        damp=1.0,
        bus_types,
        Vset
    )

Perform one Newton–Raphson step in rectangular coordinates using an
analytic Jacobian that matches `mismatch_rectangular`.

The mismatch vector F(V) is defined as in `mismatch_rectangular`:
- PQ buses: ΔP_i, ΔQ_i
- PV buses: ΔP_i, ΔV_i with ΔV_i = |V_i| - Vset[i]

Only non-slack buses are treated as variables (Vr/Vi), the slack bus
voltage is kept fixed.
"""
function complex_newton_step_rectangular(
    Ybus,
    V::Vector{ComplexF64},
    S::Vector{ComplexF64};
    slack_idx::Int = 1,
    damp::Float64  = 1.0,
    bus_types::Vector{Symbol},
    Vset::Vector{Float64}
)
    n = length(V)
    @assert length(S)         == n
    @assert length(bus_types) == n
    @assert length(Vset)      == n

    # Non-slack bus indices
    non_slack = collect(1:n)
    deleteat!(non_slack, slack_idx)

    # Base mismatch F(V)
    F0 = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)

    # Analytic Jacobian for the same F(V)
    J  = build_rectangular_jacobian_pq_pv(Ybus, V, bus_types, Vset, slack_idx)

    # Solve J * δx = -F0
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

    # Damping
    δx .*= damp

    Vr = real.(V)
    Vi = imag.(V)

    # Apply corrections only to non-slack buses
    @inbounds for (k, bus) in enumerate(non_slack)
        Vr[bus] += δx[k]
        Vi[bus] += δx[(n - 1) + k]
    end

    V_new = ComplexF64.(Vr, Vi)
    # keep slack fixed (defensive, sollte schon passen)
    V_new[slack_idx] = V[slack_idx]

    return V_new
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
#=
function run_complex_nr_rectangular_for_net!(net::Net;
                                             maxiter::Int = 20,
                                             tol::Float64 = 1e-8,
                                             damp::Float64 = 0.2,
                                             verbose::Int = 0)
    
    sparse = length(net.nodeVec) > 60
    Ybus   = createYBUS(net=net, sparse=sparse, printYBUS = (verbose > 1))    
    

    # 2) Build initial complex voltages V0 from bus data
    V0, slack_idx = initial_Vrect_from_net(net)

    # 3) Build specified complex power injections S = P + jQ from bus data
    S = build_S_from_net(net)
        # 3a) Build bus type vector and PV voltage setpoints
    nodes = net.nodeVec
    n = length(nodes)

    bus_types = Vector{Symbol}(undef, n)
    Vset      = Vector{Float64}(undef, n)

    for (k, node) in enumerate(nodes)
        BusType = getNodeType(node)        
        if BusType == Slack
            bus_types[k] = :Slack
        elseif BusType == PV
            bus_types[k] = :PV
        elseif BusType == PQ
            bus_types[k] = :PQ
        else
            error("Unsupported bus type for complex NR rectangular PF")    
        end
        Vset[k] = node._vm_pu  
    end

    
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
        bus_types = bus_types,
        Vset      = Vset,
    )


        # 5) Write back voltages to network (vm_pu, va_deg)
    update_net_voltages_from_complex!(net, V)

    # 6) Recompute bus powers from final voltages and Ybus
    nodes   = net.nodeVec
    Sbase   = net.baseMVA
    Ibus    = Ybus * V
    Sbus_pu = V .* conj.(Ibus)   # S in pu

    #TODO: write back to net bus data structure
    for (k, node) in enumerate(nodes)   
        if isSlack(node)
           Sbus = Sbus_pu[k] * Sbase
           Pbus_MW = real(Sbus)
           Qbus_MVar = imag(Sbus)

           if Pbus_MW < 0.0
               node._pƩLoad = Pbus_MW
               node._qƩLoad = Qbus_MVar
           else
               node._pƩGen = Pbus_MW
               node._qƩGen = Qbus_MVar
           end           
        end        
    end

    # Total bus power in pu (wie im klassischen Solver)
    setTotalBusPower!(net = net,p = sum(real.(Sbus_pu)),  q = sum(imag.(Sbus_pu)))

    erg = converged ? 0 : 1
    return iters, erg
end
=#
function run_complex_nr_rectangular_for_net!(net::Net;
                                             maxiter::Int = 20,
                                             tol::Float64 = 1e-8,
                                             damp::Float64 = 0.2,
                                             verbose::Int = 0,
                                             use_fd::Bool = false)

    nodes  = net.nodeVec
    n      = length(nodes)
    Sbase  = net.baseMVA
    sparse = n > 60
    Ybus   = createYBUS(net=net, sparse=sparse, printYBUS = (verbose > 1))

    # 1) Initial complex voltages V0 and slack index
    V0, slack_idx = initial_Vrect_from_net(net)

    # 2) Specified complex power injections S (p.u.), sign convention wie im Polar-NR
    S = build_S_from_net(net)

    # 3) Bus types and PV setpoints from Node data
    bus_types = Vector{Symbol}(undef, n)
    Vset      = Vector{Float64}(undef, n)

    @inbounds for (k, node) in enumerate(nodes)
        BusType = getNodeType(node)
        if BusType == Slack
            bus_types[k] = :Slack
        elseif BusType == PV
            bus_types[k] = :PV
        elseif BusType == PQ
            bus_types[k] = :PQ
        else
            error("run_complex_nr_rectangular_for_net!: unsupported bus type at bus $k")
        end
        Vset[k] = isnothing(node._vm_pu) ? 1.0 : node._vm_pu
    end

    # 4) Q-limit data (analog zu calcNewtonRaphson_withPVIdentity!)
    qmin_pu, qmax_pu = getQLimits_pu(net)
    resetQLimitLog!(net)

    cooldown_iters = hasfield(typeof(net), :cooldown_iters) ? net.cooldown_iters : 0
    q_hyst_pu      = hasfield(typeof(net), :q_hyst_pu)      ? net.q_hyst_pu      : 0.0

    # PV-Busse am Anfang (Guard für PQ->PV-Reenable)
    pv_orig = [k for k in 1:n if bus_types[k] == :PV]

    # 5) NR-Schleife
    V        = copy(V0)
    history  = Float64[]
    converged = false
    iters     = 0

    if verbose > 0
        @info "Starting rectangular complex NR power flow..."
        @info "Initial complex voltages V0:" V0
        @info "Slack bus index:" slack_idx
        @info "maxiter = $maxiter, tol = $tol, damp = $damp"
    end

    for it in 1:maxiter
        iters = it

        # Mismatch mit aktuellem bus_types & S
        F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
        max_mis = maximum(abs.(F))
        push!(history, max_mis)

        (verbose > 0) && @info "Rectangular NR iteration" iter=it max_mismatch=max_mis

        if max_mis <= tol
            converged = true
            break
        end

        # --- Q-Limit Active Set: PV -> PQ, optional PQ -> PV -------------------
        changed   = false
        reenabled = false

        if it > 2    # kleine Verzögerung wie im Polar-NR
            # aktuelles S(V) aus Netz
            I      = Ybus * V
            S_calc = V .* conj.(I)

            # 4a) PV -> PQ
            @inbounds for k in 1:n
                if k == slack_idx
                    continue
                end
                if bus_types[k] != :PV
                    continue
                end

                qreq   = imag(S_calc[k])    # aktuelle Q-Anforderung (p.u.)
                busIdx = k                  # hier ist Busindex = Position in nodeVec

                if qreq > qmax_pu[k]
                    # Bus wird PQ, Q auf Qmax geklemmt
                    bus_types[k] = :PQ
                    S[k] = complex(real(S[k]), net.qmax_pu[k])
                    changed = true

                    net.qLimitEvents[busIdx] = :max
                    logQLimitHit!(net, it, busIdx, :max)
                    (verbose>0) && @printf "PV->PQ Bus %d: Q=%.4f > Qmax=%.4f (it=%d)\n" busIdx qreq qmax_pu[k] it

                elseif qreq < qmin_pu[k]
                    bus_types[k] = :PQ
                    S[k] = complex(real(S[k]), net.qmin_pu[k])
                    changed = true

                    net.qLimitEvents[busIdx] = :min
                    logQLimitHit!(net, it, busIdx, :min)
                    (verbose>0) && @printf "PV->PQ Bus %d: Q=%.4f < Qmin=%.4f (it=%d)\n" busIdx qreq qmin_pu[k] it
                end
            end

            if changed
                F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
                max_mis = maximum(abs.(F))
            end

            # 4b) Optional PQ -> PV Re-enable (Hysterese + Cooldown)
            #     Wenn du das erstmal nicht willst: kompletten Block auskommentieren.
            if q_hyst_pu > 0.0 || cooldown_iters > 0
                @inbounds for k in 1:n
                    # Guards:
                    # (1) aktuell PQ
                    # (2) war am Anfang PV
                    # (3) hat einen Q-Limit-Event in dieser Rechnung
                    if bus_types[k] != :PQ
                        continue
                    end
                    if !(k in pv_orig)
                        continue
                    end
                    if !haskey(net.qLimitEvents, k)
                        continue
                    end

                    qreq = imag(S_calc[k])
                    lo   = net.qmin_pu[k] + q_hyst_pu
                    hi   = net.qmax_pu[k] - q_hyst_pu
                    ready = (qreq > lo) && (qreq < hi)

                    if ready && (cooldown_iters > 0)
                        last_it = lastQLimitIter(net, k)
                        if !isnothing(last_it) && (it - last_it) < cooldown_iters
                            ready = false
                        end
                    end

                    if ready
                        bus_types[k] = :PV
                        delete!(net.qLimitEvents, k)
                        reenabled = true
                        (verbose>0) && @printf "PQ->PV Bus %d: Q=%.4f in [%.4f, %.4f] (cooldown=%d)\n" k qreq lo hi cooldown_iters
                    end
                end

                if reenabled
                    F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
                    max_mis = maximum(abs.(F))
                end
            end
        end

        # --- Newton-Schritt (FD oder analytisch) -----------------------------
        if use_fd
            V = complex_newton_step_rectangular_fd(
                Ybus, V, S;
                slack_idx = slack_idx,
                damp      = damp,
                h         = 1e-6,
                bus_types = bus_types,
                Vset      = Vset,
            )
        else
            V = complex_newton_step_rectangular(
                Ybus, V, S;
                slack_idx = slack_idx,
                damp      = damp,
                bus_types = bus_types,
                Vset      = Vset,
            )
        end

        # Slackspannung festhalten (Sicherheitsgurt)
        V[slack_idx] = V0[slack_idx]
    end

    # 6) Spannungen zurück ins Netz schreiben
    update_net_voltages_from_complex!(net, V)

    # 7) Busleistungen (wie im klassischen Solver) zurückrechnen
    Ibus     = Ybus * V
    Sbus_pu  = V .* conj.(Ibus)
    Sbus_MVA = Sbus_pu .* Sbase

    @inbounds for (k, node) in enumerate(nodes)
        Sbus      = Sbus_MVA[k]
        Pbus_MW   = real(Sbus)
        Qbus_MVar = imag(Sbus)

        # sehr einfache Heuristik: P<0 → Load, P>=0 → Gen
        if Pbus_MW < 0.0
            node._pƩLoad = Pbus_MW
            node._qƩLoad = Qbus_MVar
        else
            node._pƩGen = Pbus_MW
            node._qƩGen = Qbus_MVar
        end
    end

    # 8) Gesamtbusleistung wie im Polar-NR aktualisieren    
    setTotalBusPower!(net = net,p = sum(real.(Sbus_pu)),  q = sum(imag.(Sbus_pu)))
    return iters, converged ? 0 : 1
end

"""
    runpf_rectangular!(net, maxIte, tolerance=1e-6, verbose=0)

Runs a rectangular complex-state Newton–Raphson power flow on `net::Net`.

Returns:
    (iterations::Int, status::Int)
where `status == 0` indicates convergence.
"""
function runpf_rectangular!(net::Net, maxIte::Int, tolerance::Float64=1e-6, verbose::Int=0)
    iters, erg = run_complex_nr_rectangular_for_net!(
        net;
        maxiter = maxIte,
        tol     = tolerance,
        damp    = 1.0,
        verbose = verbose,
    )
    return iters, erg
end

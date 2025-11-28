using Sparlectra
using Test
using Logging
using Printf
using LinearAlgebra

# keep logs quiet unless there's a warning or error
global_logger(ConsoleLogger(stderr, Logging.Warn))

include("testgrid.jl")
include("testremove.jl")
include("test_jacobian_full.jl")

@testset "Sparlectra.jl" begin
    @test testNetwork() == true
    @test test_NBI_MDO() == true
    @test test_acpflow(0) == true
    @test testISOBusses() == true
    @test testRemoveFunctions() == true    
    @testset "Full-system NR (PV identity rows)" begin
        @test test_acpflow_full(0) == true
        @test test_full_matches_reduced() == true
        @test test_jacobian_full_structure() == true        
        @test test_5BusNet(0,15.0) == true
        @test test_5BusNet(0,30.0) == false
    end
    @testset "Rectangular complex NR vs. runpf_full!" begin
        net1 = createTest5BusNet(pq_only=false)
        net2 = createTest5BusNet(pq_only=false)

        # Klassischer Solver
        it_polar, erg_polar = runpf_full!(net1, 20, 1e-8, 0)
        @test erg_polar == 0

        # Rechteck-Solver
        it_rect, erg_rect = run_complex_nr_rectangular_for_net!(net2; maxiter=20, tol=1e-6, damp=1.0, verbose=0)
        @test erg_rect == 0

        V_polar = [n._vm_pu * cis(deg2rad(n._va_deg)) for n in net1.nodeVec]
        V_rect  = [n._vm_pu * cis(deg2rad(n._va_deg)) for n in net2.nodeVec]

        @test maximum(abs.(abs.(V_polar) .- abs.(V_rect))) < 1e-4
        @test maximum(abs.(angle.(V_polar) .- angle.(V_rect))) < 1e-3
    end

@testset "Rectangular analytic Jacobian matches FD (3-bus example)" begin
    # 3-bus Ybus matrix (same as in example_complex_nr_3bus.jl)
    Ydata = ComplexF64[
        10 - 5im,
       -10 + 5im,
        0 + 0im,
       -10 + 5im,
        20 - 10im,
       -10 + 5im,
        0 + 0im,
       -10 + 5im,
        10 - 5im,
    ]
    Ybus = reshape(Ydata, 3, 3)

    # Initial voltages: slack at bus 1 = 1∠0, flat start at others
    V = ComplexF64[1.0 + 0im, 1.0 + 0im, 1.0 + 0im]

    # Specified power injections S = P + jQ
    S = ComplexF64[
        1.0 - 0.2im,    # slack / generator at bus 1
       -0.5 + 0.1im,    # load at bus 2
       -0.5 + 0.1im,    # load at bus 3
    ]

    slack_idx = 1
    n = length(V)

    # Bus types and voltage setpoints
    bus_types = Symbol[:Slack, :PQ, :PQ]
    Vset      = fill(1.0, n)

    # Base mismatch (already only non-slack buses)
    F0 = Sparlectra.mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
    m  = length(F0)
    @test m == 2*(n-1)

    # Build FD Jacobian in rectangular coordinates for F(V)
    non_slack = collect(1:n)
    deleteat!(non_slack, slack_idx)

    nvar = 2*(n-1)
    J_fd = zeros(Float64, m, nvar)

    h = 1e-6

    # Perturb real parts Vr(non_slack)
    Vr = real.(V)
    Vi = imag.(V)

    for (col_idx, bus) in enumerate(non_slack)
        Vp = copy(V)
        Vp[bus] = ComplexF64(Vr[bus] + h, Vi[bus])

        Fp = Sparlectra.mismatch_rectangular(Ybus, Vp, S, bus_types, Vset, slack_idx)
        J_fd[:, col_idx] .= (Fp .- F0) ./ h
    end

    # Perturb imaginary parts Vi(non_slack)
    for (offset, bus) in enumerate(non_slack)
        col_idx = (n - 1) + offset
        Vp = copy(V)
        Vp[bus] = ComplexF64(Vr[bus], Vi[bus] + h)

        Fp = Sparlectra.mismatch_rectangular(Ybus, Vp, S, bus_types, Vset, slack_idx)
        J_fd[:, col_idx] .= (Fp .- F0) ./ h
    end

    # Analytic Jacobian
    J_an = Sparlectra.build_rectangular_jacobian_pq_pv(Ybus, V, bus_types, Vset, slack_idx)

    @test size(J_an) == size(J_fd)

    diff    = J_an .- J_fd
    rel_err = norm(diff) / max(norm(J_fd), eps())

    @info "Rectangular analytic vs FD Jacobian: rel_err = $rel_err"
    @test rel_err < 1e-6
end
@testset "Rectangular analytic vs FD step" begin
    bus_types = Symbol[:Slack, :PQ, :PQ]
    Vset      = fill(1.0, 3)
    slack_idx = 1

    V0 = ComplexF64[1.0 + 0im, 1.0 + 0im, 1.0 + 0im]

    V_fd = Sparlectra.complex_newton_step_rectangular_fd(
        Ybus, V0, S;
        slack_idx = slack_idx,
        damp      = 1.0,
        h         = 1e-6,
        bus_types = bus_types,
        Vset      = Vset,
    )

    V_an = Sparlectra.complex_newton_step_rectangular(
        Ybus, V0, S;
        slack_idx = slack_idx,
        damp      = 1.0,
        bus_types = bus_types,
        Vset      = Vset,
    )

    @test isapprox(V_fd, V_an; rtol=1e-6, atol=1e-8)
end

end

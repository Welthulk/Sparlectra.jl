using Sparlectra
using Test
using Logging
using Printf

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

end

using Sparlectra
using Test
using Logging

# keep logs quiet unless there's a warning or error
global_logger(ConsoleLogger(stderr, Logging.Warn))

# existing tests
include("testgrid.jl")
include("testremove.jl")

# NEW: full Jacobian / PV-identity tests
include("test_jacobian_full.jl")

@testset "Sparlectra.jl" begin
    @test testNetwork() == true
    @test test_NBI_MDO() == true
    @test test_acpflow(0) == true
    @test testISOBusses() == true
    @test testRemoveFunctions() == true
    # @test testImportMatpower() == true  # needs file in repo

    @testset "Full-system NR (PV identity rows)" begin
        @test test_acpflow_full(0) == true
        @test test_full_matches_reduced() == true
        @test test_jacobian_full_structure() == true
        @test test_pv_q_limit_switch(verbose=1) == true
    end
end

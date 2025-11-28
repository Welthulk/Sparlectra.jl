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
    @test test_5BusNet(0, 15.0) == true
    @test test_5BusNet(0, 30.0) == false
  end
end

using Sparlectra
using Test
using Logging
using Printf
using LinearAlgebra

# keep logs quiet unless there's a warning or error
global_logger(ConsoleLogger(stderr, Logging.Warn))

include("testgrid.jl")
include("testremove.jl")

@testset "Sparlectra.jl" begin
  @test test_2WTPITrafo() == true
  @test test_3WTPITrafo() == true
  @test testNetwork() == true
  @test test_NBI_MDO() == true
  @test test_acpflow(0; lLine_6a6b = 0.01, damp = 1.0, method = :rectangular, opt_sparse = true) == true
  @test test_acpflow(0; lLine_6a6b = 0.01, damp = 1.0, method = :rectangular, opt_sparse = false) == true
  @test test_acpflow(0; lLine_6a6b = 0.01, damp = 1.0, method = :polar_full, opt_sparse = true) == true
  @test test_acpflow(0; lLine_6a6b = 0.01, damp = 1.0, method = :classic, opt_sparse = true) == true
  @test testISOBusses() == true
  @test testRemoveFunctions() == true
  @test test_5BusNet(0, 10.0) == true
  @test test_3BusNet(0, 150.0, :rectangular, false, false) == true
  @test test_3BusNet(0, 150.0, :polar_full, false, false) == true
end

using Sparlectra
using Test
using Logging

global_logger(ConsoleLogger(stderr, Logging.Warn))
include("testgrid.jl")
include("testremove.jl")
@testset "Sparlectra.jl" begin
  @test testNetwork() == true
  @test test_NBI_MDO() == true
  @test test_acpflow(0) == true
  @test testISOBusses() == true
  @test testRemoveFunctions() == true
  #@test testImportMatpower() == true # This test is not working because the file is not found in github
end

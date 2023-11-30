using Sparlectra
using Test
using Logging

global_logger(ConsoleLogger(stderr, Logging.Warn))
include("testnetworks.jl")
@testset "Sparlectra.jl" begin
  @test testNetwork() == true
  @test test_NBI_MDO() == true
  @test test_acpflow() == true
end

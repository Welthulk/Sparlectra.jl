using Sparlectra
using Test
using Logging

global_logger(ConsoleLogger(stderr, Logging.Debug))

include("testgrid.jl")
include("testpst.jl")
include("testABCD.jl")

@testset "Sparlectra.jl" begin    
  #@test testNetwork() == true
  #@test test_NBI_MDO() == true
  #@test test_acpflow(0) == true
  #@test testOpenBranches(0) == true
  #@test testParallelBranches(0) == true
  #@test testISOBusses() == true
  #@test testImportMatpower() == true
  @test test_phaseshifters(1) == true
  #@test test_ABCD(1) == true
  
end

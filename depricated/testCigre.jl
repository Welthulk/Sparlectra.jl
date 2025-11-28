using Sparlectra
using Test
using Logging
using Printf

# keep logs quiet unless there's a warning or error
global_logger(ConsoleLogger(stderr, Logging.Info))

include("testgrid.jl")



  net = testCreateNetworkFromScratch()
  @info net
  print(net.branchVec)
  tol = 1e-6
  maxIte = 10
  verbose = 1
  result = true
  etime = @elapsed begin
    ite, erg = runpf!(net, maxIte, tol, verbose)
  end
  if erg != 0
    @warn "Power flow did not converge"
    result = false
  end


  calcNetLosses!(net)
  printACPFlowResults(net, etime, ite, tol)





using Test
using Sparlectra
using Sparlectra.ResDataTypes
using Sparlectra.SparlectraExport
using Sparlectra.SparlectraImport
using Sparlectra.SparlectraNet
using Sparlectra.SparlectraResult
using BenchmarkTools
using Logging

global_logger(ConsoleLogger(stderr, Logging.Info))

function acpflow(casefile::String, iterations::Int=10, verbose::Int=0, mdo::Bool = false, tol::Float64 = 1e-6, base_MVA::Float64 = 0.0, printResultToFile::Bool = false, printResultAnyCase::Bool = false)
  
  ext = splitext(casefile)[2]
  myNet = nothing
  if ext == ".m"    
    path = "mpower"
    jpath = joinpath(pwd(), "data", path, strip(casefile))
    if !isfile(jpath)
      @error "File $(jpath) not found"
      return
    end
    myNet = createNetFromMatPowerFile(jpath, base_MVA, (verbose >= 1), mdo)
  else
    @error "File extension $(ext) not supported"
    return
  end

  
  # run power flow
  ite = 0
  etime = @elapsed begin
    ite, erg = runpf!(myNet, iterations, tol, verbose) 
  end

  if erg == 0 || printResultAnyCase
    calcNetLosses!(myNet.nodeVec, myNet.branchVec, myNet.baseMVA)
    jpath = printResultToFile ? joinpath(pwd(), "data", path) : ""
    printACPFlowResults(myNet, etime, ite, tol, printResultToFile, jpath)
  elseif erg == 1
    @warn "Newton-Raphson did not converge"
  else
    @error "error during calculation of Newton-Raphson"
  end
end

file = "case2.m"
#file = "case7.m"
#file = "case7a.m"
#file = "case7b.m"
#file = "case7c.m"
#file = "case_ieee30.m"
printResultToFile = false
printResultAnyCase = false
mdo = false
tol = 1e-8
ite = 25
verbose = 2 # 0: no output, 1: iteration norm, 2: + Y-Bus, 3: + Jacobian, 4: + Power Flow
base_MVA = 0.0 # # Set a value to a non-zero value to override the corresponding value in the case file.
@time acpflow(file, ite, verbose, mdo, tol, base_MVA, printResultToFile, printResultAnyCase)
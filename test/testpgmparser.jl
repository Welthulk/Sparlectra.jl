using Test
using Sparlectra
using Sparlectra.ResDataTypes
using Sparlectra.SparlectraTools
using Sparlectra.SparlectraExport
using Sparlectra.SparlectraImport
using Sparlectra.SparlectraNet
using Sparlectra.SparlectraResult
using BenchmarkTools
using Logging

global_logger(ConsoleLogger(stderr, Logging.Info))

function acpflow(casefile::String, writeCase::Bool, iterations::Int, verbose::Int)
  sparse=true  
  tol = 1e-6     
  log = verbose >= 1

  # load + create net
  jpath=joinpath(pwd(),"data","pgmodel",casefile)
  if !isfile(jpath)
    @error "File $(jpath) not found"
    return
  end
   
  myNet =createNetFromPGM(jpath)
  filename=joinpath(pwd(),"data","pgmodel",myNet.name * "_epx")
  exportPGM(myNet,filename)
  # matpower export  
  if writeCase       
    case = myNet.name
    file =joinpath(pwd(),"data","pgmodel",case * ".m")    
    writeMatpowerCasefile(myNet, file, case)    
  end

  # Y-Bus Matrix  
  Y = createYBUS(myNet.branchVec, myNet.shuntVec, (verbose>2), sparse)
  ite = 0
  etime = @elapsed begin
    ite,erg = calcNewtonRaphson!(Y, myNet.nodeVec, myNet.baseMVA, iterations, tol, verbose, sparse)      
  end      
  
  if erg == 0
   calcNetLosses!(myNet.nodeVec, myNet.branchVec, myNet.baseMVA, log)    
   printACPFlowResults(myNet, etime, ite)
  elseif erg == 1
    @warn "Newton-Raphson did not converge"
  else
    @error "error during calculation of Newton-Raphson"  
  end
  
end

@time acpflow("input.json", true,  8, 1)
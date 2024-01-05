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
  sparse=false  
  tol = 1e-6   
  base_MVA = 100.0
  log = verbose >= 1

  # load + create net
  jpath=joinpath(pwd(),"data","ppower",casefile)
  if !isfile(jpath)
    @error "File $(jpath) not found"
    return
  end
   
  myNet = createNetFromFile(jpath,base_MVA, (verbose>2))  
  case = myNet.name
  
   
  # matpower export
  file =joinpath(pwd(),"Sparlectra","data","ppower",case * ".m")
  if writeCase       
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

@time acpflow("bsp7.json", false,  6, 0)
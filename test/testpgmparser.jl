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

function acpflow(casefile::String, writeCase::Bool, writePGM::Bool, iterations::Int, verbose::Int, tol::Float64=1e-6)
  sparse=false  
  # load + create net
  jpath=joinpath(pwd(),"data","pgmodel",casefile)
  if !isfile(jpath)
    @error "File $(jpath) not found"
    return
  end
  myNet =createNetFromPGM(jpath)
  if length(myNet.nodeVec) == 0
    @error "No nodes found in $(jpath)"
    return
  elseif length(myNet.nodeVec) > 30
    sparse = true
  end
  # matpower export  
  if writeCase       
    case = myNet.name
    file =joinpath(pwd(),"data","pgmodel",case * ".m")    
    writeMatpowerCasefile(myNet, file, case)    
  end

  # Y-Bus Matrix  
  Y = createYBUS(myNet.branchVec, myNet.shuntVec, sparse, (verbose>2))
  ite = 0
  etime = @elapsed begin
    ite,erg = calcNewtonRaphson!(Y, myNet.nodeVec, myNet.baseMVA, iterations, tol, verbose, sparse)      
  end      
  
  if erg == 0
   calcNetLosses!(myNet.nodeVec, myNet.branchVec, myNet.baseMVA)    
   printACPFlowResults(myNet, etime, ite)
   if writePGM
    convertPVtoPQ!(myNet)
    filename=joinpath(pwd(),"data","pgmodel",myNet.name * "_epx")
    exportPGM(net=myNet,filename=filename,useMVATrafoModell=false,exportSlackGen=true)
   end
  elseif erg == 1
    @warn "Newton-Raphson did not converge"
  else
    @error "error during calculation of Newton-Raphson"  
  end
  
end

#@time acpflow("input.json", true, true, 10, 1)
#]@time acpflow("MiniGrid_v4.json", true, true, 10, 1)
@time acpflow("SmallGrid.json", true, true, 20, 1, 1e-6)
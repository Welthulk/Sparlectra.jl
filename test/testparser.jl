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

function acpflow(casefile::String, writeCase::Bool, writePGM::Bool, iterations::Int, verbose::Int, mdo::Bool = false, tol::Float64 = 1e-6, base_MVA::Float64 = 0.0, printResultToFile::Bool = false, printResultAnyCase::Bool = false)
  sparse = false
  ext = splitext(casefile)[2]
  myNet = nothing
  path = ""
  pgmpath = "pgmodel"
  if ext == ".json"
    # load + create net    
    writePGM = false
    path = "pgmodel"
    jpath = joinpath(pwd(), "data", path, casefile)
    path = joinpath(pwd(), "data", "pgmodel")
    if !isfile(jpath)
      @error "File $(jpath) not found"
      return
    end

    myNet = createNetFromPGM(jpath)
  elseif ext == ".m"
    writeCase = false
    path = "mpower"
    jpath = joinpath(pwd(), "data", path, casefile)
    if !isfile(jpath)
      @error "File $(jpath) not found"
      return
    end
    myNet = createNetFromMatPowerFile(jpath, base_MVA, (verbose >= 1), mdo)
  else
    @error "File extension $(ext) not supported"
    return
  end

  if length(myNet.nodeVec) == 0
    @error "No nodes found in $(jpath)"
    return
  elseif length(myNet.nodeVec) > 30
    sparse = true
  end
  # matpower export  
  if writeCase
    case = myNet.name
    file = joinpath(pwd(), "data", path, case * ".m")
    writeMatpowerCasefile(myNet, file, case)
  end

  # Y-Bus Matrix  
  Y = createYBUS(myNet.branchVec, myNet.shuntVec, sparse, (verbose > 0))
  ite = 0
  etime = @elapsed begin
    ite, erg = calcNewtonRaphson!(Y, myNet.nodeVec, myNet.baseMVA, iterations, tol, verbose, sparse)
  end

  if erg == 0 || printResultAnyCase
    calcNetLosses!(myNet.nodeVec, myNet.branchVec, myNet.baseMVA)
    if printResultToFile
      jpath = joinpath(pwd(), "data", path)
      @info "...export results to $(jpath)"
    else
      jpath = ""
    end
    printACPFlowResults(myNet, etime, ite, tol, printResultToFile, jpath)
    if writePGM
      convertPVtoPQ!(myNet)
      filename = joinpath(pwd(), "data", pgmpath, strip(myNet.name) * "_epx")
      exportPGM(net = myNet, filename = filename, useMVATrafoModell = false, exportSlackGen = true)
    end
  elseif erg == 1
    @warn "Newton-Raphson did not converge"
  else
    @error "error during calculation of Newton-Raphson"
  end
end

#file = "MiniGrid_v4.json"
#file = "SmallGrid.json"
file = "case7.m"
#file = "case30.m"
printResultToFile = false
printResultAnyCase = false
writeCase = false
writePGM = true
mdo = false
tol = 1e-6
ite = 20
verbose = 0
base_MVA = 0.0 # # Set a value to a non-zero value to override the corresponding value in the case file.
@time acpflow(file, writeCase, writePGM, ite, verbose, mdo, 1e-6, base_MVA, printResultToFile, printResultAnyCase)
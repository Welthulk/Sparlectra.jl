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

function acpflow(casefile::String, iterations::Int, verbose::Int, mdo::Bool, exportToPGM::Bool=false, printResultAnyCase::Bool=false, printResultToFile::Bool=false)
  sparse = true
  tol = 1e-6
  base_MVA = 0.0 # # Set a value to a non-zero value to override the corresponding value in the case file.
  # load + create net
  jpath = joinpath(pwd(), "data", "mpower", casefile)
  if !isfile(jpath)
    @error "File $(jpath) not found"
    return
  end

  @time myNet = SparlectraNet.createNetFromMatPowerFile(jpath, base_MVA, (verbose >= 1), mdo)

  if exportToPGM    
    filename = joinpath(pwd(), "data", "pgmodel", myNet.name * "_pgm")
    @info "...export to PGM-Files, Filename: ($filename)"
    exportPGM(net=myNet, filename=filename, useMVATrafoModell=false)
  end

  # Y-Bus Matrix
  @time Y = SparlectraNet.createYBUS(myNet.branchVec, myNet.shuntVec, sparse, (verbose > 2))

  ite = 0
  etime = @elapsed begin
    ite, erg = calcNewtonRaphson!(Y, myNet.nodeVec, myNet.baseMVA, iterations, tol, verbose, sparse)
  end

  if erg == 0 || printResultAnyCase    
    calcNetLosses!(myNet.nodeVec, myNet.branchVec, myNet.baseMVA)
    if printResultToFile
      jpath = joinpath(pwd(), "data", "mpower")
      @info "...export results to $(jpath)"
    else
      jpath = ""
    end
    printACPFlowResults(myNet, etime, ite, printResultToFile, jpath)
  elseif erg == 1
    @warn "Newton-Raphson did not converge"
  else
    @error "error during calculation of Newton-Raphson"
  end

end

verbose = 0
ite = 6
mdo = false
pgm_export = true
showAnyResult = false
writeResultToFile = false
#case = "case_ieee118.m"
#case = "case_ieee30.m"
case="ieee_30bus.m"
@time acpflow(case, ite, verbose, mdo, pgm_export, showAnyResult, writeResultToFile)



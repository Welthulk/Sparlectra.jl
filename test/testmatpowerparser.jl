using Test
using Sparlectra
using Sparlectra.ResDataTypes
using Sparlectra.SparlectraExport
using Sparlectra.SparlectraImport
using Sparlectra.SparlectraNet
using Sparlectra.SparlectraResult
using Sparlectra.SparlectraTools
using BenchmarkTools
using Logging


global_logger(ConsoleLogger(stderr, Logging.Info))

function acpflow(casefile::String, iterations::Int, verbose::Int, mdo::Bool, printResultAnyCase::Bool=false, printResultToFile::Bool=false)
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

  # Y-Bus Matrix
  @time Y = SparlectraNet.createYBUS(myNet.branchVec, myNet.shuntVec, sparse, (verbose > 2))

  ite = 0
  etime = @elapsed begin
    ite, erg = calcNewtonRaphson!(Y, myNet.nodeVec, myNet.baseMVA, iterations, tol, verbose, sparse)
  end

  if erg == 0 || printResultAnyCase
    calcNetLosses!(myNet.nodeVec, myNet.branchVec, myNet.baseMVA, (verbose > 1))
    if printResultToFile
      jpath = joinpath(pwd(), "data", "mpower")
      @info "...export results to $(jpath)"
    else
      jpath = ""
    end
    printACPFlowResults(myNet, etime, ite, printResultToFile, jpath)
    if printResultToFile
      @info "...done"
    end
  elseif erg == 1
    @warn "Newton-Raphson did not converge"
  else
    @error "error during calculation of Newton-Raphson"
  end

end

verbose = 1
ite = 6
mdo = false
showAnyResult = false
writeResultToFile = true
case = "case_ieee118.m"

@time acpflow(case, ite, verbose, mdo, showAnyResult, writeResultToFile)



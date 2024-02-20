using Test
using Sparlectra
using Sparlectra.ResDataTypes
using Sparlectra.SparqlQueryCGMES
using Sparlectra.SparlectraTools
using Sparlectra.SparlectraExport
using Sparlectra.SparlectraImport
using Sparlectra.SparlectraNet
using Sparlectra.SparlectraResult
using BenchmarkTools
using Logging

global_logger(ConsoleLogger(stderr, Logging.Warn))
# endpoint = Jena Fuseki-Server
const endpoint = "http://localhost:3030/explore2/query"
sparse = false
iterations = 12
tol = 1e-6
log = false
slackbusName = "HG2" #HG2_Node-9"
Netname = "MiniGrid"
sbase_mva = 1.0
verbose = 0
writeCasePGM = true 
writeCaseMP = false
@time myNet = SparlectraNet.createNetFromTripleStore(endpoint, sbase_mva, Netname, slackbusName, log)
case = myNet.name

if writeCaseMP
  # matpower export
  file = joinpath(pwd(), "data", "MiniGrid", case * ".m")
  SparlectraExport.writeMatpowerCasefile(myNet, file, case)
end

if writeCasePGM  
    # pgm export
    file = joinpath(pwd(),"data","pgmodel",myNet.name * "_pgm")
    exportPGM(myNet,file)
end

# calculate Y-Bus Matrix
Y = SparlectraNet.createYBUS(myNet.branchVec, myNet.shuntVec, log, sparse)

# calculate Newton-Raphson
ite = 0
etime = @elapsed begin
  ite, erg = calcNewtonRaphson!(Y, myNet.nodeVec, myNet.baseMVA, iterations, tol, verbose, sparse)
end

if erg == 0
  calcNetLosses!(myNet.nodeVec, myNet.branchVec, myNet.baseMVA, log)
  printACPFlowResults(myNet, etime, ite)
elseif erg == 1
  @warn "Newton Raphson did not converge"
else
  @error "error during calculation of Newton-Raphson"
end

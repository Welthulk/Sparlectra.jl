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

global_logger(ConsoleLogger(stderr, Logging.Debug))

function acpflow(casefile::String, writeCase::Bool, iterations::Int, verbose::Int)
  sparse=true  
  tol = 1e-6   
  base_MVA = 100.0
  log = verbose >= 1

  # load + create net
  jpath=joinpath(pwd(),"data","pgmodel",casefile)
  if !isfile(jpath)
    @error "File $(jpath) not found"
    return
  end
   
  createNetFromPGM(jpath)
  
  #myNet = createNetFromFile(jpath,base_MVA, (verbose>2))  
  #case = myNet.name
  
  
  
  #= 
    
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
  =#
end

#=
# Funktion zur Berechnung von Yshunt und Umrechnung in Mikrosiemens
function calculate_yshunt(frequency::Float64, capacitance::Float64, tan_delta::Float64)::Float64
    omega = 2 * π * frequency
    y_shunt = omega * capacitance / (tan_delta + im)

    # Umrechnung von Siemens zu Mikrosiemens
    y_shunt_microsiemens = abs(y_shunt) * 1e6
    return y_shunt_microsiemens
end

# Beispielaufruf
frequency = 50.0  # Annahme: Frequenz in Hertz
capacitance = 1.10248e-07  # Annahme: Kapazität in Farad
tan_delta = 0.1  # Annahme: tan delta

yshunt_microsiemens = calculate_yshunt(frequency, capacitance, tan_delta)

println("Yshunt in Mikrosiemens: $yshunt_microsiemens µS")
=#


@time acpflow("input.json", false,  6, 1)
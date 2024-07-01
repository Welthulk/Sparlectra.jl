using Sparlectra
using BenchmarkTools
using Logging


function test_ABCD(verbose::Int = 0)
  
  # Netzwerk erstellen
  Sbase_MVA = 100.0  # Basis-Scheinleistung in MVA
  netName = "testnet"
  net = Net(name=netName, baseMVA=Sbase_MVA)

  # Busse hinzufügen
  addBus!(net = net, busName = "B1", busType = "Slack", vn_kV = 220.0)  
  addBus!(net= net,  busName=  "B2", busType = "PQ", vn_kV = 220.0)

  addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 200.0, q = 10.0)

  # Transformator hinzufügen
  #add2WTrafo!(net=net, fromBus="Bus1", toBus="Bus2", vnom_from_kV=110.0, vnom_to_kV=10.0, snom_MVA=50.0, vnom_percent=10.0, tap=1.0)

  # YBus-Matrix berechnen
  Y_bus_1 = createYBUS(net=net, sparse=false, printYBUS=(verbose > 0))
  Y_bus_2 = createYBUS_ABCD(net=net, sparse=false, printYBUS=(verbose > 0))
  
  br = getNetBranch(net = net, fromBus = "B1", toBus = "B2")
  if verbose > 0
    @show br
  end
  if isapprox(Y_bus_1, Y_bus_2, atol=1e-6)
    return true
  else
    return false
  end

end


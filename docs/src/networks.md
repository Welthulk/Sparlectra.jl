Network Module
=============
The `Net` module provides functionality for creating and manipulating power system network models in Julia. It includes features for defining buses, branches, transformers, prosumers, and shunts, as well as methods for running power flow analysis.

```julia
# Import the Net module
  using Sparlectra
  using BenchmarkTools
  net = Net(name = "case5", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", busType = "PQ", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "B2", busType = "PQ", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "B3", busType = "PQ", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "B4", busType = "PQ", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "B5", busType = "Slack", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)

  addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 25.0, r = 0.2, x = 0.39)
  addACLine!(net = net, fromBus = "B1", toBus = "B3", length = 25.0, r = 0.2, x = 0.39)
  addACLine!(net = net, fromBus = "B2", toBus = "B4", length = 25.0, r = 0.2, x = 0.39)
  addACLine!(net = net, fromBus = "B3", toBus = "B4", length = 25.0, r = 0.2, x = 0.39)
  addACLine!(net = net, fromBus = "B4", toBus = "B5", length = 25.0, r = 0.2, x = 0.39)

  addProsumer!(net = net, busName = "B1", type = "ENERGYCONSUMER", p = 1.0, q = 2.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 1.0, q = 2.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 1.0, q = 2.0)

  addProsumer!(net = net, busName = "B5", type = "SYNCHRONMASCHINE", referencePri = "B5", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "B1", type = "GENERATOR", p = 1.1, q = 2.0)

  # Run power flow
  tol = 1e-6
  maxIte = 10  
  etime = @elapsed begin
    ite, erg = runpf!(net, maxIte, tol, 0)
  end
  if erg != 0
    @warn "Power flow did not converge"        
  else
    calcNetLosses!(net)
    printACPFlowResults(net, etime, ite, tol)
  end
  
```

```@autodocs 
  Modules = [Sparlectra]   
  Pages = ["network.jl"]
  Order = [:type, :function]
```  

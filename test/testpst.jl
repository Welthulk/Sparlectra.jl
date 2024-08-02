using Sparlectra
using BenchmarkTools
using Logging


function test_phaseshifters(verbose::Int = 0)
  function run_acpflow(net::Net)
    
    if verbose < 1
      return
    end

    tol = 1e-3
    maxIte = 10
  
    result = true
    etime = @elapsed begin
      ite, erg = runpf!(net, maxIte, tol, verbose)
    end
    if erg != 0
      @warn "test_ABCD: Power flow did not converge"
      result = false
    end

    if verbose > 0
      calcNetLosses!(net)
      printACPFlowResults(net, etime, ite, tol)
    end
  end

  #=
  x_α_max = 8.0
  maxStep = 12
  neutralStep = 0
  step = 0
  δu = 0.1
  x_0 = 4.0

  myPST = SymmetricalPhaseShifter(vn=380.0, from=1, to=2, neutralStep=neutralStep, step=step, maxStep=maxStep, δu=δu, x_0=x_0, x_α_max=x_α_max)
  @debug myPST
  for i in neutralStep:maxStep
    setCurrentStep(myPST, i)
    @debug getCurrentAngle(myPST)
    @debug getCurrentX(myPST)
  end
  =#
  #=

       G -> 1-----------------2         4------------------5 <- G
                              |---PST---| 
       L <- 3-----------------2         4------------------6 -> L

  =#

  net = Net(name="psttest", baseMVA=1000.0)  
  vn_kV_b1 = 220.0
  vn_kV_b2 = 220.0
  addBus!(net = net, busName = "B1", busType = "Slack", vn_kV = vn_kV_b1)  
  addBus!(net= net,  busName=  "B2", busType = "PV", vn_kV = vn_kV_b2)
  addBus!(net= net,  busName=  "B3", busType = "PQ", vn_kV = vn_kV_b2)
  
  addBus!(net= net,  busName=  "B4", busType = "PV", vn_kV = vn_kV_b2)
  #addBus!(net= net,  busName=  "B5", busType = "Slack", vn_kV = vn_kV_b2)
  addBus!(net= net,  busName=  "B5", busType = "PV", vn_kV = vn_kV_b2)
  addBus!(net= net,  busName=  "B6", busType = "PQ", vn_kV = vn_kV_b2)
  

  
  addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 200.0, r = 0.02, x = 0.4, c_nf_per_km = 9.08, tanδ = 0.0)
  addACLine!(net = net, fromBus = "B3", toBus = "B2", length = 200.0, r = 0.02, x = 0.4, c_nf_per_km = 9.08, tanδ = 0.0)
  addACLine!(net = net, fromBus = "B5", toBus = "B4", length = 10.0, r = 0.02, x = 0.4, c_nf_per_km = 9.08, tanδ = 0.0)
  addACLine!(net = net, fromBus = "B4", toBus = "B6", length = 10.0, r = 0.02, x = 0.4, c_nf_per_km = 9.08, tanδ = 0.0)
  addPST!(net = net, fromBus = "B4", toBus = "B2", sn_mva=1000.0, neutralStep = 0, maxStep = 12, step = 0, δu = 0.1, x_0 = 4.0, x_α_max = 8.0, status=1)
  
  
  addProsumer!(net = net, busName = "B1", type = "SYNCHRONOUSMACHINE", p = 500.0, vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")  
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 50.0, q = 5.0)
  
  
  addProsumer!(net = net, busName = "B5", type = "SYNCHRONOUSMACHINE", p = 55.0, vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "B6", type = "ENERGYCONSUMER", p = 100.0, q = 5.0)

  
  
  
  createYBUS(net=net, sparse=false, printYBUS=(verbose > 0))
  
  run_acpflow(net)
  
  return true
end

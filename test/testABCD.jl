using Sparlectra
using BenchmarkTools
using Logging


function test_ABCD(verbose::Int = 0)
  function run_acpflow(net::Net, verbose::Int = 0)
    
    if verbose < 2
      return
    end

    tol = 1e-6
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

  function createNet(mode)
    if verbose > 0  @info "createNet: $mode" end
    Sbase_MVA = 100.0  # Basis-Scheinleistung in MVA
    netName = "testnet-ABCD"
    net = Net(name=netName, baseMVA=Sbase_MVA)  
    if mode == "line"
     vn_kV_b1 = 220.0
     vn_kV_b2 = 220.0

    elseif mode =="2WT"
      vn_kV_b1 = 110.0
      vn_kV_b2 = 10.0

    elseif mode =="combine"
     vn_kV_b1 = 220.0
     vn_kV_b2 = 220.0
     vn_kV_b3 = 110.0
    elseif mode == "longline"
     vn_kV_b1 = 220.0
     vn_kV_b2 = 220.0
    else
     @error "Unknown mode: $mode"
     return
    end
    
    # Busse hinzufügen
    addBus!(net = net, busName = "B1", busType = "Slack", vn_kV = vn_kV_b1)  
    addBus!(net= net,  busName=  "B2", busType = "PQ", vn_kV = vn_kV_b2)
    if mode == "combine"
      addBus!(net= net,  busName=  "B3", busType = "PQ", vn_kV = vn_kV_b3)
    end
    
    if mode == "line" 
      addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
    elseif mode == "longline"
      addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 5500.0, r = 0.05, x = 0.1, c_nf_per_km = 1000.0, tanδ = 5.0, isLongLine = true)
    elseif mode =="2WT" 
      add2WTrafo!(net=net, fromBus="B1", toBus="B2", sn_mva = Sbase_MVA, vk_percent = 10.0, vkr_percent = 0.1, pfe_kw = 0.0, i0_percent = 0.0)
    elseif mode == "combine"
      addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 100.0, r = 0.0653, x = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
      add2WTrafo!(net=net, fromBus="B2", toBus="B3", sn_mva = Sbase_MVA, vk_percent = 10.0, vkr_percent = 0.1, pfe_kw = 0.0, i0_percent = 0.0)
    end  

    if mode == "combine"
      addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
      addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 200.0, q = 10.0)
    elseif mode == "longline"
      addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
      addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 20.0, q = 1.0)
    else
      addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
      addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 200.0, q = 10.0)
    end


    return net
  end
  result = true  
  
  # --------------------------------------------------------------------------------
  net = createNet("line")  
  Y_bus_1 = createYBUS(net=net, sparse=false, printYBUS=(verbose > 0))
  Y_bus_2 = createYBUS_ABCD(net=net, sparse=false, printYBUS=(verbose > 0))  
  br = getNetBranch(net = net, fromBus = "B1", toBus = "B2")
  if verbose > 0
    @show br
  end
  if !isapprox(Y_bus_1, Y_bus_2, atol=1e-6)
    result =  false
  end  
  run_acpflow(net, verbose)

  # --------------------------------------------------------------------------------
  net = createNet("2WT")  
  Y_bus_1 = createYBUS(net=net, sparse=false, printYBUS=(verbose > 0))
  Y_bus_2 = createYBUS_ABCD(net=net, sparse=false, printYBUS=(verbose > 0))
  br = getNetBranch(net = net, fromBus = "B1", toBus = "B2")
  if verbose > 0
    @show br
  end
  if !isapprox(Y_bus_1, Y_bus_2, atol=1e-6)
    result =  false
  end
  run_acpflow(net, verbose)  
  
  # --------------------------------------------------------------------------------
  net = createNet("combine")  
  Y_bus_1 = createYBUS(net=net, sparse=false, printYBUS=(verbose > 0))
  Y_bus_2 = createYBUS_ABCD(net=net, sparse=false, printYBUS=(verbose > 0))
  if !isapprox(Y_bus_1, Y_bus_2, atol=1e-6)
    result =  false
  end
  run_acpflow(net, verbose)


  # --------------------------------------------------------------------------------
  net = createNet("longline")  
  Y_bus_1 = createYBUS(net=net, sparse=false, printYBUS=(verbose > 0))
  Y_bus_2 = createYBUS_ABCD(net=net, sparse=false, printYBUS=(verbose > 0))  
  br = getNetBranch(net = net, fromBus = "B1", toBus = "B2")
  #if !isapprox(Y_bus_1, Y_bus_2, atol=1e-6)
  #  result =  false
  #end  
  run_acpflow(net, 2)


  return result
end


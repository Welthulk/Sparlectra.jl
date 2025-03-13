# Network Analyzer for Sparlectra

using Sparlectra

"""
    analyze_network(net::Net)

Analyze a Sparlectra network structure and print detailed information about its components.
Helpful for debugging issues with removal operations.
"""
function analyze_network(net::Net)
  println("\n====================== NETWORK ANALYSIS ======================")
  println("Network: $(net.name) (BaseMVA: $(net.baseMVA))")

  # Bus analysis
  println("\n--- BUSES ($(length(net.nodeVec))) ---")
  for (i, node) in enumerate(net.nodeVec)
    bus_type = toString(node._nodeType)
    slack = isSlack(node) ? "SLACK" : ""
    println("$i: $(node.comp.cName) (busIdx: $(node.busIdx), type: $bus_type) $slack")
  end

  # Branch analysis
  println("\n--- BRANCHES ($(length(net.branchVec))) ---")
  for (i, branch) in enumerate(net.branchVec)
    println("$i: $(branch.comp.cName) (from: $(branch.fromBus), to: $(branch.toBus), status: $(branch.status))")
  end

  # Bus dictionary mapping
  println("\n--- BUS DICTIONARY ---")
  for (name, idx) in net.busDict
    println("$name => $idx")
  end

  # Shunt analysis
  println("\n--- SHUNTS ($(length(net.shuntVec))) ---")
  for (i, shunt) in enumerate(net.shuntVec)
    println("$i: $(shunt.comp.cName) (busIdx: $(shunt.busIdx))")
  end

  # Shunt dictionary mapping
  println("\n--- SHUNT DICTIONARY ---")
  for (busIdx, shuntIdx) in net.shuntDict
    bus_name = "Unknown"
    for (name, idx) in net.busDict
      if idx == busIdx
        bus_name = name
        break
      end
    end
    println("Bus $busIdx ($bus_name) => Shunt $shuntIdx")
  end

  # Prosumer analysis
  println("\n--- PROSUMERS ($(length(net.prosumpsVec))) ---")
  for (i, prosumer) in enumerate(net.prosumpsVec)
    p_type = toString(prosumer.proSumptionType)
    bus_idx = prosumer.comp.cFrom_bus
    bus_name = "Unknown"
    for (name, idx) in net.busDict
      if idx == bus_idx
        bus_name = name
        break
      end
    end
    println("$i: $(prosumer.comp.cName) (busIdx: $bus_idx ($bus_name), type: $p_type)")
  end

  # Isolated nodes
  println("\n--- ISOLATED NODES ($(length(net.isoNodes))) ---")
  for idx in net.isoNodes
    println("$idx")
  end

  println("\n================================================================")
end

"""
    test_debug()

Create a test network and analyze it with detailed debug information.
This helps understand the network structure before and after removal operations.
"""
function test_debug()
  # Create a test network as in testremove.jl
  net = Net(name = "remove_test", baseMVA = 100.0)

  # Add buses
  addBus!(net = net, busName = "B1", busType = "PQ", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", busType = "PQ", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", busType = "PQ", vn_kV = 110.0)
  addBus!(net = net, busName = "B4", busType = "PQ", vn_kV = 110.0)
  addBus!(net = net, busName = "B5", busType = "SLACK", vn_kV = 110.0)
  addBus!(net = net, busName = "B6", busType = "PV", vn_kV = 20.0)

  # Add lines
  addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 10.0, r = 0.01, x = 0.1, c_nf_per_km = 10.0, tanδ = 0.0, ratedS = 100.0)
  addACLine!(net = net, fromBus = "B2", toBus = "B3", length = 15.0, r = 0.01, x = 0.1, c_nf_per_km = 10.0, tanδ = 0.0, ratedS = 100.0)
  addACLine!(net = net, fromBus = "B3", toBus = "B4", length = 20.0, r = 0.01, x = 0.1, c_nf_per_km = 10.0, tanδ = 0.0, ratedS = 100.0)
  addACLine!(net = net, fromBus = "B4", toBus = "B5", length = 25.0, r = 0.01, x = 0.1, c_nf_per_km = 10.0, tanδ = 0.0, ratedS = 100.0)

  # Add transformer
  add2WTrafo!(net = net, fromBus = "B1", toBus = "B6", sn_mva = 100.0, vk_percent = 10.0, vkr_percent = 0.5, pfe_kw = 20.0, i0_percent = 0.1)

  # Add shunts
  addShunt!(net = net, busName = "B2", pShunt = 0.0, qShunt = 10.0)
  addShunt!(net = net, busName = "B4", pShunt = 0.0, qShunt = 15.0)

  # Add prosumers
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 50.0, q = 20.0)
  addProsumer!(net = net, busName = "B5", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.03, va_deg = 0.0, referencePri = "B5")
  addProsumer!(net = net, busName = "B6", type = "GENERATOR", p = 80.0, q = 30.0, vm_pu = 1.02)

  println("\n=================== INITIAL NETWORK STATE ===================")
  analyze_network(net)

  # Test removing a branch
  println("\n\n=================== AFTER BRANCH REMOVAL ===================")
  println("Removing branch 2...")
  removeBranch!(net = net, branchNr = 2)
  analyze_network(net)

  return net
end

# Call the test_debug function
net = test_debug()
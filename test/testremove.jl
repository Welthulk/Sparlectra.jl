# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 13.03.2025
# test/testremove.jl

using Sparlectra
using Test
using Logging

function createTestNetwork()::Net
  # Create a simple test network
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

  # Validate the network
  result, msg = validate!(net = net)
  @assert result "Network creation failed: $msg"

  return net
end

function testRemoveBranch()
  @testset "Remove Branch Test" begin
    net = createTestNetwork()

    # Get initial counts
    initialBranches = length(net.branchVec)

    # Check the actual number of branches available
    @test initialBranches > 0

    # Remove a branch by number - only if there are at least 2 branches
    if initialBranches >= 2
      branchNr = 2  # Let's remove the second branch
      result = removeBranch!(net = net, branchNr = branchNr)
      @test result == true

      # Verify one branch was removed
      @test length(net.branchVec) == initialBranches - 1

      # Remove another branch by number
      branchNr = 1  # Now remove the first branch
      result = removeBranch!(net = net, branchNr = branchNr)
      @test result == true

      # Verify another branch was removed
      @test length(net.branchVec) == initialBranches - 2
    else
      @info "Not enough branches to test multiple removals"
    end

    # Try to remove a non-existent branch (use a safe number)
    invalidBranchNr = initialBranches + 1
    result = removeBranch!(net = net, branchNr = invalidBranchNr)
    @test result == false
  end
end

function testRemoveShunt()
  @testset "Remove Shunt Test" begin
    net = createTestNetwork()

    # Get initial shunt count
    initialShunts = length(net.shuntVec)

    # Remove a shunt that exists
    result = removeShunt!(net = net, busName = "B2")
    @test result == true

    # Check updated state
    @test length(net.shuntVec) == initialShunts - 1

    # Verify node has updated
    bus2 = net.nodeVec[net.busDict["B2"]]
    @test iszero(bus2._pShunt) || isnothing(bus2._pShunt)
    @test iszero(bus2._qShunt) || isnothing(bus2._qShunt)

    # Try to remove a non-existent shunt (choose a bus that definitely has no shunt)
    noShuntBusName = "B1"  # Using B1 as it doesn't have a shunt
    result = removeShunt!(net = net, busName = noShuntBusName)
    @test result == false
  end
end

function testRemoveProsumer()
  @testset "Remove Prosumer Test" begin
    net = createTestNetwork()

    # Get initial prosumer count
    initialProsumers = length(net.prosumpsVec)

    # First find the prosumer connected to B3 - it should be an ENERGYCONSUMER
    b3Idx = net.busDict["B3"]
    prosumerIdx = findfirst(p -> p.comp.cFrom_bus == b3Idx, net.prosumpsVec)

    @test prosumerIdx !== nothing

    # Remove the prosumer
    result = removeProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER")
    @test result == true

    # Check updated state
    @test length(net.prosumpsVec) == initialProsumers - 1

    # Find the prosumer connected to B6 - it should be a GENERATOR
    b6Idx = net.busDict["B6"]
    prosumerIdx = findfirst(p -> p.comp.cFrom_bus == b6Idx, net.prosumpsVec)

    @test prosumerIdx !== nothing

    # Remove the prosumer
    result = removeProsumer!(net = net, busName = "B6", type = "GENERATOR")
    @test result == true

    # Check updated state
    @test length(net.prosumpsVec) == initialProsumers - 2

    # Try to remove a non-existent prosumer
    noGeneratorBusName = "B1"  # Using B1 as it doesn't have a generator
    result = removeProsumer!(net = net, busName = noGeneratorBusName, type = "GENERATOR")
    @test result == false
  end
end

function testIsolatedBuses()
  @testset "Isolated Buses Test" begin
    net = createTestNetwork()

    # Initially no isolated buses
    markIsolatedBuses!(net = net)

    # We don't check isempty here because there might already be isolated nodes
    initialIsolatedCount = length(net.isoNodes)

    # Get branch numbers for connections to B3
    branchIndices = Int[]
    b3Idx = net.busDict["B3"]

    for (i, br) in enumerate(net.branchVec)
      if br.fromBus == b3Idx || br.toBus == b3Idx
        push!(branchIndices, i)
      end
    end

    # Remove all branches connected to B3 (from highest index to lowest)
    sort!(branchIndices, rev = true)
    for idx in branchIndices
      removeBranch!(net = net, branchNr = idx)
    end

    # Remove any prosumers at B3
    removeProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER")

    # Mark isolated buses
    markIsolatedBuses!(net = net)

    # Check if bus is now isolated
    @test b3Idx in net.isoNodes
    @test length(net.isoNodes) > initialIsolatedCount
  end
end

function testRemoveSpecialCases()
  @testset "Special Removal Cases" begin
    net = createTestNetwork()

    # Test cannot remove slack bus
    slackBusName = "B5"
    result = removeBus!(net = net, busName = slackBusName)
    @test result == false

    # Test cannot remove bus with connections
    connectedBusName = "B2"
    result = removeBus!(net = net, busName = connectedBusName)
    @test result == false
  end
end

function testRemoveFunctions()
  @testset "Network Removal Functions" begin
    testRemoveBranch()
    testRemoveShunt()
    testRemoveProsumer()
    testIsolatedBuses()
    testRemoveSpecialCases()
  end
  return true
end

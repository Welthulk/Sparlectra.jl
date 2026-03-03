# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Network Analyzer for Sparlectra

using Sparlectra

"""
    analyze_network(net::Net; verbose_components::Bool=true, io::IO=stdout)

Analyze a `Sparlectra.Net` and print:
- a compact summary (number of buses/branches/shunts/prosumers),
- optional component-level details,
- dictionary mappings,
- consistency checks (invalid references/mismatches).

Useful for debugging network edits and remove operations.
"""
function analyze_network(net::Net; verbose_components::Bool=true, io::IO=stdout)
  bus_name_by_idx = Dict{Int, String}(idx => String(name) for (name, idx) in net.busDict)

  println(io, "\n====================== NETWORK ANALYSIS ======================")
  println(io, "Network: $(net.name) (BaseMVA: $(net.baseMVA))")
  println(io, "Summary: buses=$(length(net.nodeVec)), branches=$(length(net.branchVec)), shunts=$(length(net.shuntVec)), prosumers=$(length(net.prosumpsVec)), isolated=$(length(net.isoNodes))")

  if verbose_components
    println(io, "\n--- BUSES ($(length(net.nodeVec))) ---")
    for (i, node) in enumerate(net.nodeVec)
      bus_type = toString(node._nodeType)
      slack = isSlack(node) ? " [SLACK]" : ""
      println(io, "$i: $(node.comp.cName) (busIdx=$(node.busIdx), type=$bus_type)$slack")
    end

    println(io, "\n--- BRANCHES ($(length(net.branchVec))) ---")
    for (i, branch) in enumerate(net.branchVec)
      from_name = get(bus_name_by_idx, branch.fromBus, "Unknown")
      to_name = get(bus_name_by_idx, branch.toBus, "Unknown")
      println(io, "$i: $(branch.comp.cName) (from=$(branch.fromBus):$from_name, to=$(branch.toBus):$to_name, status=$(branch.status))")
    end

    println(io, "\n--- SHUNTS ($(length(net.shuntVec))) ---")
    for (i, shunt) in enumerate(net.shuntVec)
      bus_name = get(bus_name_by_idx, shunt.busIdx, "Unknown")
      println(io, "$i: $(shunt.comp.cName) (busIdx=$(shunt.busIdx):$bus_name)")
    end

    println(io, "\n--- PROSUMERS ($(length(net.prosumpsVec))) ---")
    for (i, prosumer) in enumerate(net.prosumpsVec)
      p_type = toString(prosumer.proSumptionType)
      bus_idx = prosumer.comp.cFrom_bus
      bus_name = get(bus_name_by_idx, bus_idx, "Unknown")
      println(io, "$i: $(prosumer.comp.cName) (busIdx=$bus_idx:$bus_name, type=$p_type)")
    end
  end

  println(io, "\n--- BUS DICTIONARY ---")
  for (name, idx) in sort(collect(net.busDict); by=x -> x[2])
    println(io, "$name => $idx")
  end

  println(io, "\n--- SHUNT DICTIONARY ---")
  for (busIdx, shuntIdx) in sort(collect(net.shuntDict); by=x -> x[1])
    bus_name = get(bus_name_by_idx, busIdx, "Unknown")
    println(io, "Bus $busIdx ($bus_name) => Shunt $shuntIdx")
  end

  println(io, "\n--- ISOLATED NODES ($(length(net.isoNodes))) ---")
  for idx in sort(collect(net.isoNodes))
    bus_name = get(bus_name_by_idx, idx, "Unknown")
    println(io, "$idx ($bus_name)")
  end

  _print_consistency_report(net, bus_name_by_idx; io)
  println(io, "\n================================================================")
end

"""
    build_demo_network() -> Net

Build a compact demo grid that can be used for analyzer debugging.
"""
function build_demo_network()::Net
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

  return net
end

"""
    run_demo() -> Net

Build a demo network, run analysis before and after removing branch 2.
Returns the modified network.
"""
function run_demo()::Net
  net = build_demo_network()

  println("\n=================== INITIAL NETWORK STATE ===================")
  analyze_network(net)

  println("\n\n=================== AFTER BRANCH REMOVAL ===================")
  println("Removing branch 2...")
  removeBranch!(net = net, branchNr = 2)
  analyze_network(net)

  return net
end

function _print_consistency_report(net::Net, bus_name_by_idx::Dict{Int, String}; io::IO=stdout)
  issues = String[]
  bus_count = length(net.nodeVec)
  shunt_count = length(net.shuntVec)

  # busDict indices should point to valid nodes
  for (name, idx) in net.busDict
    if !(1 <= idx <= bus_count)
      push!(issues, "busDict[$name] = $idx is out of bounds (1:$bus_count)")
    elseif String(name) != net.nodeVec[idx].comp.cName
      push!(issues, "busDict[$name] = $idx points to node $(net.nodeVec[idx].comp.cName)")
    end
  end

  # branch endpoints should reference valid buses
  for (i, branch) in enumerate(net.branchVec)
    if !(1 <= branch.fromBus <= bus_count)
      push!(issues, "branch[$i] has invalid fromBus=$(branch.fromBus)")
    end
    if !(1 <= branch.toBus <= bus_count)
      push!(issues, "branch[$i] has invalid toBus=$(branch.toBus)")
    end
  end

  # shunt vector and shuntDict should reference valid buses/shunt indices
  for (i, shunt) in enumerate(net.shuntVec)
    if !(1 <= shunt.busIdx <= bus_count)
      push!(issues, "shunt[$i] has invalid busIdx=$(shunt.busIdx)")
    end
  end
  for (busIdx, shuntIdx) in net.shuntDict
    if !(1 <= busIdx <= bus_count)
      push!(issues, "shuntDict has invalid busIdx=$busIdx")
    end
    if !(1 <= shuntIdx <= shunt_count)
      push!(issues, "shuntDict[$busIdx] has invalid shunt index=$shuntIdx")
    end
  end

  # prosumers must point to valid buses
  for (i, prosumer) in enumerate(net.prosumpsVec)
    bus_idx = prosumer.comp.cFrom_bus
    if !(1 <= bus_idx <= bus_count)
      push!(issues, "prosumer[$i] has invalid cFrom_bus=$bus_idx")
    end
  end

  # isolated node list should reference existing buses
  for idx in net.isoNodes
    if !(1 <= idx <= bus_count)
      push!(issues, "isoNodes contains invalid bus index=$idx")
    end
  end

  println(io, "\n--- CONSISTENCY CHECKS ---")
  if isempty(issues)
    println(io, "OK: no structural issues detected.")
  else
    println(io, "Found $(length(issues)) issue(s):")
    for msg in issues
      println(io, "- $msg")
    end
  end

  # quick reachability info for isolated nodes
  if !isempty(net.isoNodes)
    names = [get(bus_name_by_idx, idx, "Unknown") for idx in sort(collect(net.isoNodes))]
    println(io, "Isolated buses: ", join(names, ", "))
  end
end

if abspath(PROGRAM_FILE) == @__FILE__
  run_demo()
end

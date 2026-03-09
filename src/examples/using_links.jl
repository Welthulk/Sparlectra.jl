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

"""
using_links.jl

Example network with a busbar coupler (`Bus1` ↔ `Bus1a`) represented as a
Sparlectra bus link.

Requested topology:
- Branches: Bus1→Bus3, Bus2→Bus4, Bus3→Bus4, Bus1a→Bus4
- Additional topological link: Bus1↔Bus1a
- Injections at Bus1, Bus4 and Bus5 (Bus5 is slack)

The script runs two scenarios:
1) Link closed (`status = 1`)
2) Link open   (`status = 0`)
"""

using Sparlectra

function build_link_demo_net(; link_closed::Bool)
  net = Net(name = link_closed ? "using_links_closed" : "using_links_open", baseMVA = 100.0)

  # Buses: five base buses + additional Bus1a.
  addBus!(net = net, busName = "Bus1", busType = "PQ", vn_kV = 110.0)
  addBus!(net = net, busName = "Bus1a", busType = "PQ", vn_kV = 110.0)
  addBus!(net = net, busName = "Bus2", busType = "PQ", vn_kV = 110.0)
  addBus!(net = net, busName = "Bus3", busType = "PQ", vn_kV = 110.0)
  addBus!(net = net, busName = "Bus4", busType = "PQ", vn_kV = 110.0)
  addBus!(net = net, busName = "Bus5", busType = "Slack", vn_kV = 110.0)

  # Requested branches.
  addPIModelACLine!(net = net, fromBus = "Bus1", toBus = "Bus3", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "Bus2", toBus = "Bus4", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "Bus3", toBus = "Bus4", r_pu = 0.008, x_pu = 0.060, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "Bus1a", toBus = "Bus4", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)

  # Tie to slack area.
  addPIModelACLine!(net = net, fromBus = "Bus4", toBus = "Bus5", r_pu = 0.006, x_pu = 0.050, b_pu = 0.0, status = 1)

  # Busbar coupler / sectionalizer between Bus1 and Bus1a.
  addLink!(net = net, fromBus = "Bus1", toBus = "Bus1a", status = link_closed ? 1 : 0)

  # Injections as requested: Bus1, Bus4 and Bus5 (Slack).
  addProsumer!(net = net, busName = "Bus1", type = "GENERATOR", p = 45.0, q = 0.0, vm_pu = 1.01)
  addProsumer!(net = net, busName = "Bus4", type = "GENERATOR", p = 35.0, q = 0.0, vm_pu = 1.00)
  addProsumer!(net = net, busName = "Bus5", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Bus5")

  # Loads to create non-trivial power transfers.
  addProsumer!(net = net, busName = "Bus2", type = "LOAD", p = 20.0, q = 6.0)
  addProsumer!(net = net, busName = "Bus3", type = "LOAD", p = 30.0, q = 4.0)
  addProsumer!(net = net, busName = "Bus4", type = "LOAD", p = 25.0, q = 8.0)
  addProsumer!(net = net, busName = "Bus1a", type = "LOAD", p = 30.0, q = 10.0)

  return net
end

function run_link_demo(; link_closed::Bool)
  net = build_link_demo_net(link_closed = link_closed)
  iterations = 0
  status = -1
  elapsed_s = @elapsed begin
    iterations, status = runpf!(net, 25)
  end

  println("\n====================")
  println("Scenario: ", link_closed ? "Link CLOSED (Bus1-Bus1a)" : "Link OPEN (Bus1-Bus1a)")
  println("Converged: ", status == 0, " (status=", status, ", iterations=", iterations, ")")

  if status == 0
    # runpf! solves node voltages and injections.
    # Branch flow/loss tables in printACPFlowResults are populated by calcNetLosses!.
    calcNetLosses!(net)
    calcLinkFlowsKCL!(net)
    printACPFlowResults(net, elapsed_s, iterations, 1e-8)
  else
    println("Power flow did not converge.")
  end

  return net
end

# Run both variants for direct comparison.
run_link_demo(link_closed = true)
run_link_demo(link_closed = false)

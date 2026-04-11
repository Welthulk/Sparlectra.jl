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

# file: src/examples/using_netreports.jl

using Sparlectra
using Printf
using DataFrames

"""
Build a minimal network, run AC power flow and create `ACPFlowReport`.

This example demonstrates:
- creating a network in code,
- solving power flow with `run_net_acpflow`,
- reading structured results from `ACPFlowReport`,
- including bus-link flows (`report.links`),
- conversion to DataFrame.
"""
function main()
  net = Net(name = "netreport_demo", baseMVA = 100.0)

  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)
  addBus!(net = net, busName = "B4", vn_kV = 110.0)

  addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 5.0, r = 0.05, x = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "B3", toBus = "B4", length = 20.0, r = 0.05, x = 0.5, c_nf_per_km = 10.0, tanδ = 0.0)

  # Link between two PQ buses: not part of YBUS, but reported via KCL allocation.
  addLink!(net = net, fromBus = "B2", toBus = "B3", status = 1)

  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  addProsumer!(net = net, busName = "B2", type = "GENERATOR", p = 10.0, q = 1.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
  addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 25.0, q = 10.0)

  ite, erg, etime = run_net_acpflow(net = net, max_ite = 40, tol = 1e-10, method = :polar_full, opt_sparse = true, show_results = false)

  report = buildACPFlowReport(net; ct = etime, ite = ite, tol = 1e-10, converged = (erg == 0), solver = :polar_full)

  println("\n--- ACPFlowReport summary ---")
  println(report)
  @printf("Converged: %s, iterations: %d, elapsed: %.4f s\n", string(report.metadata.converged), report.metadata.iterations, report.metadata.elapsed_s)
  @printf("Total losses: P = %.4f MW, Q = %.4f MVar\n", report.metadata.total_p_loss_MW, report.metadata.total_q_loss_MVar)

  println("\nNode table rows:")
  for row in report.nodes
    println(row)
  end

  println("\nLink table rows:")
  for row in report.links
    println(row)
  end

  nodes_df = DataFrame(report.nodes)
  branches_df = DataFrame(report.branches)
  links_df = DataFrame(report.links)

  println("\nDataFrame(nodes):")
  show(nodes_df; allrows = true, allcols = true)
  println("\n\nDataFrame(branches):")
  show(branches_df; allrows = true, allcols = true)
  println("\n\nDataFrame(links):")
  show(links_df; allrows = true, allcols = true)
  println()

  return report
end

main()
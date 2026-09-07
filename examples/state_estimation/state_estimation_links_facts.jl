# Copyright 2023-2026 Udo Schmitz
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

# Date: 2026-08-26
# file: examples/state_estimation/state_estimation_links_facts.jl
# purpose: SE phase 3 walkthrough. A busbar coupler (closed link) splits one
#          electrical node into two modelled buses: the SE runs on the
#          contracted net, aggregates the cluster injections (LINKAGG), and
#          syncs the member voltages. A link P measurement then steers the W2
#          flow allocation (calcLinkFlowsSE!), and se_view reports the frozen
#          operating point.

using Sparlectra
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

function _create_demo_net()
  net = Net(name = "se_links_demo", baseMVA = 100.0)
  for b in ("S", "B1", "B1a", "B2")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addACLine!(net = net, fromBus = "S", toBus = "B1", length = 10.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "B1a", toBus = "B2", length = 8.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addProsumer!(net = net, busName = "S", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "S")
  addProsumer!(net = net, busName = "B1", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  addProsumer!(net = net, busName = "B1a", type = "ENERGYCONSUMER", p = 5.0, q = 2.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 20.0, q = 6.0)
  # the busbar coupler: B1 and B1a are one electrical node
  addLink!(net = net, fromBus = "B1", toBus = "B1a", status = 1)
  ok, msg = validate!(net = net)
  ok || error("demo net invalid: $msg")
  return net
end

function run_state_estimation_links_facts()
  print_example_banner("examples/state_estimation/state_estimation_links_facts.jl", "SE on the contracted net (busbar coupler), LINKAGG aggregation, W2 allocation with a link measurement, se_view")

  net = _create_demo_net()
  ite_pf, erg_pf = runpf!(net, 40, 1e-10, 0; method = :rectangular)
  erg_pf == 0 || error("Power flow did not converge")
  calcNetLosses!(net)
  @printf("PF iterations: %d\n\n", ite_pf)

  # 1) SE on the contracted net: the coupler fuses B1/B1a, injections on the
  # cluster aggregate to one LINKAGG measurement per kind
  std = measurementStdDevs(vm = 1e-3, pinj = 0.1, qinj = 0.1, pflow = 0.1, qflow = 0.1)
  meas = generateMeasurementsFromPF(net; noise = false, stddev = std)
  append!(net.measurements, meas)
  res = runse!(net; maxIte = 30, tol = 1e-10, updateNet = true)
  b1 = geNetBusIdx(net = net, busName = "B1")
  b1a = geNetBusIdx(net = net, busName = "B1a")
  println("SE on the contracted net:")
  @printf("  converged = %s; vm(B1) = %.6f pu == vm(B1a) = %.6f pu (members share the representative)\n\n", string(res.converged), net.nodeVec[b1]._vm_pu, net.nodeVec[b1a]._vm_pu)

  # 2) W2 allocation: without a measurement the coupler flow is the KCL
  # (minimum-norm) split; a link P measurement then pins it
  rep0 = calcLinkFlowsSE!(net)
  p0 = net.linkVec[1].pFlow_MW
  @printf("coupler flow (no measurement) : %.3f MW  source = %s\n", p0, string(rep0[1].source))
  addPflowMeasurement!(net; value = p0 + 0.5, sigma = 0.01, linkNr = 1)
  rep1 = calcLinkFlowsSE!(net)
  @printf("coupler flow (measured z)     : %.3f MW  source = %s  residual = %.4f MW\n\n", net.linkVec[1].pFlow_MW, string(rep1[1].source), rep1[1].p_meas_residual)

  # 3) the frozen operating point at a glance
  view = se_view(net)
  print_se_view(view)
  return nothing
end

run_example(run_state_estimation_links_facts)

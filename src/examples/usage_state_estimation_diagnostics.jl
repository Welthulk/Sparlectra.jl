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
#
# Example: Usage of bad-data diagnostics for state estimation.
#
# Workflow:
# 1) Build a small 3-bus demo network
# 2) Run PF and generate synthetic measurements
# 3) Inject one artificial bad measurement
# 4) Run validate_measurements(...) and inspect ranking
# 5) Run runse_diagnostics(...; deactivate_and_rerun=true)

using Sparlectra
using Random
using Printf

function _create_demo_net()
  net = Net(name = "se_diagnostics_demo", baseMVA = 100.0)

  addBus!(net = net, busName = "Slack", busType = "SLACK", vn_kV = 110.0)
  addBus!(net = net, busName = "LoadA", busType = "PQ", vn_kV = 110.0)
  addBus!(net = net, busName = "LoadB", busType = "PQ", vn_kV = 110.0)

  addACLine!(net = net, fromBus = "Slack", toBus = "LoadA", length = 12.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "LoadA", toBus = "LoadB", length = 9.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "Slack", toBus = "LoadB", length = 11.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "LoadA", type = "ENERGYCONSUMER", p = 35.0, q = 12.0)
  addProsumer!(net = net, busName = "LoadB", type = "ENERGYCONSUMER", p = 28.0, q = 9.0)

  return net
end

function run_usage_state_estimation_diagnostics()
  net = _create_demo_net()

  ite_pf, erg_pf = runpf!(net, 40, 1e-10, 0; method = :polar_full, opt_sparse = true)
  erg_pf == 0 || error("Power flow did not converge")

  std = measurementStdDevs(vm = 0.01, pinj = 1.0, qinj = 1.0, pflow = 1.0, qflow = 1.0)
  meas = generateMeasurementsFromPF(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = true, stddev = std, rng = MersenneTwister(1234))
  isempty(meas) && error("No measurements generated")

  vm_idx = findfirst(m -> m.typ == Sparlectra.VmMeas, meas)
  isnothing(vm_idx) && error("No Vm measurement found for outlier injection")
  idx = something(vm_idx, 0)
  m = meas[idx]
  meas[idx] = Measurement(typ = m.typ, value = m.value + 0.2, sigma = m.sigma, active = m.active, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id)

  report = validate_measurements(net, meas; maxIte = 12, tol = 1e-6, flatstart = true, jacEps = 1e-6, normalizedThreshold = 3.0)
  @printf("PF iterations: %d\n", ite_pf)

  print_se_diagnostics(report; topN = 8)

  #println("\n--- Same report as Markdown ---")
  #print_se_diagnostics(report; topN = 8, format = :markdown)

  diag = runse_diagnostics(net, meas; deactivate_and_rerun = true, maxIte = 12, tol = 1e-6, flatstart = true, jacEps = 1e-6, normalizedThreshold = 3.0)

  println("\n--- Diagnostics with deactivate-and-rerun ---")
  print_se_diagnostics(diag; topN = 8)

  return nothing #diag
end

run_usage_state_estimation_diagnostics()

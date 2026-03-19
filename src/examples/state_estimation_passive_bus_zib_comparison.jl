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
# Example: compare state estimation with and without ZIB (zero-injection bus)

# file: src/examples/state_estimation_passive_bus_zib_comparison.jl

using Sparlectra
using Printf
using Random

function build_passive_transit_example_net()
  net = Net(name = "se_passive_transit_example", baseMVA = 100.0)

  addBus!(net = net, busName = "Slack", busType = "SLACK", vn_kV = 110.0)
  addBus!(net = net, busName = "Transit", busType = "PQ", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", busType = "PQ", vn_kV = 110.0)

  addACLine!(net = net, fromBus = "Slack", toBus = "Transit", length = 10.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "Transit", toBus = "Load", length = 8.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = 40.0, q = 15.0)

  ok, msg = validate!(net = net)
  ok || error("Validation failed: $(msg)")
  return net
end

function _reset_non_slack_buses!(net::Net)
  for node in net.nodeVec
    if getNodeType(node) != :Slack
      node._vm_pu = 1.0
      node._va_deg = 0.0
    end
  end
  return nothing
end

function build_sparse_measurements_for_zib_demo(net::Net; rng::AbstractRNG = MersenneTwister(7))
  ite_pf, status_pf = runpf!(net, 40, 1e-10, 0; method = :polar_full, opt_sparse = true)
  status_pf == 0 || error("Power flow did not converge")

  std = measurementStdDevs(vm = 1e-5, pinj = 1e-5, qinj = 1e-5, pflow = 1e-5, qflow = 1e-5)
  all_meas = generateMeasurementsFromPF(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false, stddev = std, rng = rng)

  sparse = Measurement[]
  for m in all_meas
    keep = (m.typ == Sparlectra.VmMeas && m.busIdx == 1) || ((m.typ == Sparlectra.PflowMeas || m.typ == Sparlectra.QflowMeas) && m.branchIdx == 1 && m.direction == :from) || ((m.typ == Sparlectra.PflowMeas || m.typ == Sparlectra.QflowMeas) && m.branchIdx == 2 && m.direction == :to)
    keep && push!(sparse, m)
  end

  return ite_pf, sparse
end

function _measurement_target(m::Measurement)
  if !isnothing(m.busIdx)
    return "bus=$(m.busIdx)"
  elseif !isnothing(m.branchIdx)
    return "branch=$(m.branchIdx), dir=$(m.direction)"
  end
  return "-"
end

function _print_measurements(title::String, measurements::Vector{Measurement})
  println(title)
  println(repeat("-", length(title)))
  for m in measurements
    @printf("%-10s %-22s value=%11.6f sigma=%9.6f\n", string(m.typ), _measurement_target(m), m.value, m.sigma)
  end
  println()
  return nothing
end

function _run_scenario(label::String, base_net::Net, measurements::Vector{Measurement})
  net = deepcopy(base_net)
  _reset_non_slack_buses!(net)
  empty!(net.measurements)
  append!(net.measurements, measurements)

  obs = evaluate_global_observability(net; flatstart = false, jacEps = 1e-6)
  se = runse!(net; maxIte = 12, tol = 1e-8, flatstart = false, jacEps = 1e-6, updateNet = false)

  println(label)
  println(repeat("=", length(label)))
  @printf("Measurements: %d\n", length(measurements))
  @printf("Quality: %s\n", string(obs.quality))
  @printf("Redundancy: %d\n", obs.redundancy)
  @printf("Numerical critical measurements: %s\n", isempty(obs.numerical_critical_measurement_indices) ? "none" : string(obs.numerical_critical_measurement_indices))
  @printf("Structural critical measurements: %s\n", isempty(obs.structural_critical_measurement_indices) ? "none" : string(obs.structural_critical_measurement_indices))
  @printf("SE converged: %s in %d iterations\n", string(se.converged), se.iterations)
  @printf("Objective J: %.6e\n", se.objectiveJ)
  @printf("Residual norm: %.6e\n\n", se.residualNorm)

  return obs, se
end

function run_passive_bus_zib_comparison()
  net = build_passive_transit_example_net()
  ite_pf, sparse_meas = build_sparse_measurements_for_zib_demo(net)
  passive = findPassiveBuses(net)

  println("State estimation comparison with and without ZIB")
  println("===============================================")
  @printf("PF iterations for synthetic reference: %d\n", ite_pf)
  @printf("Detected passive buses: %s\n\n", string(passive))

  _print_measurements("Sparse measurement set without ZIB", sparse_meas)

  with_zib = copy(sparse_meas)
  added_zib = addZeroInjectionMeasurements!(with_zib; net = net, sigma = 1e-6)
  _print_measurements("Added zero-injection measurements (ZIB)", added_zib)

  _run_scenario("Scenario A: without ZIB", net, sparse_meas)
  _run_scenario("Scenario B: with ZIB", net, with_zib)

  return nothing
end

run_passive_bus_zib_comparison()

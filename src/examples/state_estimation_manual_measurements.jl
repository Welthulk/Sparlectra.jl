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
# Example: state estimation with a manually assembled measurement set
#
# This script demonstrates how to use the new measurement helper functions
# (`addVmMeasurement!`, `addPinjMeasurement!`, ...) to build a measurement
# vector without manual `Measurement(...)` construction.

# file: src/examples/state_estimation_manual_measurements.jl

using Sparlectra
using Printf
using Random

function build_manual_measurement_example_net()
  net = Net(name = "se_manual_measurements", baseMVA = 100.0)

  addBus!(net = net, busName = "B1", busType = "Slack", vn_kV = 110.0, vm_pu = 1.02, va_deg = 0.0)
  addBus!(net = net, busName = "B2", busType = "PQ", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "B3", busType = "PV", vn_kV = 110.0, vm_pu = 1.01, va_deg = 0.0)

  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B1", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)

  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)
  addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 55.0, q = 8.0, vm_pu = 1.01)
  addProsumer!(net = net, busName = "B2", type = "LOAD", p = 40.0, q = 15.0)
  addProsumer!(net = net, busName = "B3", type = "LOAD", p = 12.0, q = 4.0)

  ok, msg = validate!(net = net)
  ok || error("Validation failed: $(msg)")
  return net
end

function build_manual_measurements(net::Net; rng::AbstractRNG = MersenneTwister(42))
  ite_pf, status_pf = runpf!(net, 40, 1e-10, 0; method = :polar_full, opt_sparse = true)
  status_pf == 0 || error("Power flow did not converge")

  std = measurementStdDevs(vm = 1e-3, pinj = 0.8, qinj = 0.8, pflow = 0.5, qflow = 0.5)
  synthetic = generateMeasurementsFromPF(
    net;
    includeVm = true,
    includePinj = true,
    includeQinj = true,
    includePflow = true,
    includeQflow = true,
    noise = true,
    stddev = std,
    rng = rng,
  )

  function pick_measurement(typ::Sparlectra.MeasurementType; busIdx::Union{Nothing,Int} = nothing, branchIdx::Union{Nothing,Int} = nothing, direction::Symbol = :none)
    idx = findfirst(i -> begin
      m = synthetic[i]
      return m.typ == typ && m.busIdx == busIdx && m.branchIdx == branchIdx && m.direction == direction
    end, eachindex(synthetic))
    isnothing(idx) && error("Synthetic measurement not found for $(typ)")
    return synthetic[idx]
  end

  vm_b1 = pick_measurement(Sparlectra.VmMeas; busIdx = 1)
  vm_b2 = pick_measurement(Sparlectra.VmMeas; busIdx = 2)
  pinj_b2 = pick_measurement(Sparlectra.PinjMeas; busIdx = 2)
  qinj_b2 = pick_measurement(Sparlectra.QinjMeas; busIdx = 2)
  pflow_12 = pick_measurement(Sparlectra.PflowMeas; branchIdx = 1, direction = :from)
  qflow_12 = pick_measurement(Sparlectra.QflowMeas; branchIdx = 1, direction = :from)
  pflow_31 = pick_measurement(Sparlectra.PflowMeas; branchIdx = 3, direction = :from)

  empty!(net.measurements)
  addVmMeasurement!(net; busName = "B1", value = vm_b1.value, sigma = vm_b1.sigma, id = "VM_B1")
  addVmMeasurement!(net; busName = "B2", value = vm_b2.value, sigma = vm_b2.sigma, id = "VM_B2")
  addPinjMeasurement!(net; busName = "B2", value = pinj_b2.value, sigma = pinj_b2.sigma, id = "PINJ_B2")
  addQinjMeasurement!(net; busName = "B2", value = qinj_b2.value, sigma = qinj_b2.sigma, id = "QINJ_B2")
  addPflowMeasurement!(net; fromBus = "B1", toBus = "B2", value = pflow_12.value, sigma = pflow_12.sigma, direction = :from, id = "PFLOW_12_FROM")
  addQflowMeasurement!(net; branchNr = 1, value = qflow_12.value, sigma = qflow_12.sigma, direction = :from, id = "QFLOW_12_FROM")
  addPflowMeasurement!(net; fromBus = "B3", toBus = "B1", value = pflow_31.value, sigma = pflow_31.sigma, direction = :from, id = "PFLOW_31_FROM")

  return ite_pf
end

function reset_to_flatstart!(net::Net)
  for node in net.nodeVec
    if getNodeType(node) != :Slack
      node._vm_pu = 1.0
      node._va_deg = 0.0
    end
  end
  return nothing
end

function run_manual_measurement_example()
  net = build_manual_measurement_example_net()
  ite_pf = build_manual_measurements(net)
  reset_to_flatstart!(net)

  gobs = evaluate_global_observability(net; flatstart = false, jacEps = 1e-6)
  println("Manual measurement example")
  println("==========================")
  @printf("PF iterations used for synthetic data: %d\n", ite_pf)
  @printf("Measurements: %d, states: %d, quality: %s\n", gobs.n_measurements, gobs.n_states, string(gobs.quality))

  se = runse!(net; maxIte = 12, tol = 1e-6, flatstart = false, jacEps = 1e-6, updateNet = true)
  @printf("SE converged: %s in %d iterations\n", string(se.converged), se.iterations)
  @printf("Objective J: %.6e\n", se.objectiveJ)
  @printf("Residual norm: %.6e\n", se.residualNorm)

  println("\nMeasurements")
  println("------------")
  for m in net.measurements
    target = isnothing(m.busIdx) ? "branch=$(m.branchIdx), dir=$(m.direction)" : "bus=$(m.busIdx)"
    @printf("%-16s %-24s value=%10.5f  sigma=%8.5f\n", string(m.typ), target, m.value, m.sigma)
  end

  println("\nEstimated bus voltages")
  println("----------------------")
  for (i, node) in enumerate(net.nodeVec)
    @printf("%-4s Vm=%8.5f pu  Va=%9.4f deg\n", node.comp.cName, node._vm_pu, node._va_deg)
  end

  return se
end

run_manual_measurement_example()

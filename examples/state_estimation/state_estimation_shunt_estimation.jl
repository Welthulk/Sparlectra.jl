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
# file: examples/state_estimation/state_estimation_shunt_estimation.jl
# purpose: shunt estimation (SE phase 2). Case A: a reactor's model
#          susceptance is off by 20 percent; a direct bay Q measurement plus
#          a local Vm measurement recover the true B as an estimator state
#          (setShuntEstimation!), with write-back only behind updateShunts.
#          Case B: without a direct Q measurement, deriveShuntPseudoMeasurements!
#          derives a protected SHDERIV pseudo-measurement from the bay current
#          and the measured voltage.

using Sparlectra
using Random
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

function _create_demo_net()
  net = Net(name = "se_shunt_demo", baseMVA = 100.0)
  for b in ("Slack", "LoadA", "LoadB")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addACLine!(net = net, fromBus = "Slack", toBus = "LoadA", length = 12.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "LoadA", toBus = "LoadB", length = 9.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "Slack", toBus = "LoadB", length = 11.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "LoadA", type = "ENERGYCONSUMER", p = 35.0, q = 12.0)
  addProsumer!(net = net, busName = "LoadB", type = "ENERGYCONSUMER", p = 28.0, q = 9.0)
  # a 30-MVar (at 1 pu) compensation reactor at LoadB; the PF runs with this
  # TRUE susceptance, the estimation then starts from a stale model value
  addShunt!(net = net, busName = "LoadB", pShunt = 0.0, qShunt = 30.0)
  ok, msg = validate!(net = net)
  ok || error("demo net invalid: $msg")
  return net
end

function run_state_estimation_shunt_estimation()
  print_example_banner("examples/state_estimation/state_estimation_shunt_estimation.jl", "shunt estimation: recover a stale reactor susceptance (case A) and derive a bay pseudo-measurement (case B)")

  # --- case A: direct bay measurement, B as estimator state ---------------
  net = _create_demo_net()
  ite_pf, erg_pf = runpf!(net, 40, 1e-10, 0; method = :rectangular)
  erg_pf == 0 || error("Power flow did not converge")
  @printf("PF iterations: %d\n\n", ite_pf)

  std = measurementStdDevs(vm = 1e-3, pinj = 0.1, qinj = 0.1, pflow = 0.1, qflow = 0.1, shuntq = 0.1)
  meas = generateMeasurementsFromPF(net; includeShuntQ = true, noise = false, stddev = std)

  sh = net.shuntVec[1]
  b_true = imag(sh.y_pu_shunt)
  sh.y_pu_shunt = complex(real(sh.y_pu_shunt), 1.2 * b_true)   # stale model: +20 percent
  setShuntEstimation!(net; busName = "LoadB")

  res = runse!(net, meas; maxIte = 30, tol = 1e-10, updateNet = false)
  row = only(res.shuntEstimates)
  println("case A (direct ShuntQ measurement, B released as a state):")
  @printf("  B_model = %.4f pu (stale, +20%%)   B_est = %.4f pu   B_true = %.4f pu   frozen = %s\n", row.B_model, row.B_est, b_true, string(row.frozen))
  @printf("  model untouched without updateShunts: imag(y_pu_shunt) = %.4f pu\n", imag(sh.y_pu_shunt))
  runse!(net, meas; maxIte = 30, tol = 1e-10, updateNet = false, updateShunts = true)
  @printf("  after updateShunts = true          : imag(y_pu_shunt) = %.4f pu\n\n", imag(sh.y_pu_shunt))

  # --- case B: no direct Q measurement, derive it from bay current + Vm ---
  net2 = _create_demo_net()
  runpf!(net2, 40, 1e-10, 0; method = :rectangular)
  updateShuntPowers!(net = net2)
  sh2 = net2.shuntVec[1]
  V2 = buildVoltageVector(net2)
  i2 = sh2.busIdx
  i_A = 1000.0 * sqrt(sh2.p_shunt^2 + sh2.q_shunt^2) / (sqrt(3.0) * getNodeVn(net2.nodeVec[i2]) * abs(V2[i2]))
  addImagMeasurement!(net2; value = i_A, sigma = 2.0, busName = "LoadB")        # bay current
  addVmMeasurement!(net2; busName = "LoadB", value = abs(V2[i2]), sigma = 0.002) # the hard prerequisite
  derived = deriveShuntPseudoMeasurements!(net2)
  d = only(derived)
  println("case B (bay current + measured voltage, no direct Q measurement):")
  @printf("  derived %s = %.3f MVar (solved q_shunt = %.3f MVar), sigma = %.3f MVar\n", d.id, d.value, sh2.q_shunt, d.sigma)
  println("  the SHDERIV pseudo-measurement is protected from bad-data elimination.")
  return nothing
end

run_example(run_state_estimation_shunt_estimation)

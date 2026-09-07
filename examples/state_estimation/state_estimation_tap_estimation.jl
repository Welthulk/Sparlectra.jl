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

# Date: 2026-08-27
# file: examples/state_estimation/state_estimation_tap_estimation.jl
# purpose: transformer tap estimation (0.10.0): a transformer whose
#          actual tap disagrees with the model poisons the estimate around
#          it. setTapEstimation! releases the tap as an extra WLS state;
#          after convergence the estimator FIXES it to the nearest
#          mechanical step and reruns without the tap state, reporting J
#          before versus after the fixation. Shows the exact integer-step
#          recovery, the honest off-grid half-step case, the write-back
#          protection, and the machine-trafo back-calculation.

using Sparlectra
using Random
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

# two voltage levels meshed through two parallel transformer paths: the
# closed loop plus terminal voltages is what makes a tap observable
function _create_tap_net()
  net = Net(name = "se_tap_demo", baseMVA = 100.0)
  for (b, vn) in (("H1", 110.0), ("H2", 110.0), ("L1", 20.0), ("L2", 20.0))
    addBus!(net = net, busName = b, vn_kV = vn)
  end
  addProsumer!(net = net, busName = "H1", type = "EXTERNALNETWORKINJECTION", referencePri = "H1", vm_pu = 1.02, va_deg = 0.0)
  addProsumer!(net = net, busName = "L1", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
  addProsumer!(net = net, busName = "L2", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
  addPIModelACLine!(net = net, fromBus = "H1", toBus = "H2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "L1", toBus = "L2", r_pu = 0.02, x_pu = 0.1, b_pu = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "H1", toBus = "L1", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "H2", toBus = "L2", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  # branch 3 (H1-L1) carries the tap changer under study
  net.branchVec[3].has_ratio_tap = true
  net.branchVec[3].tap_step = 0.00625
  return net
end

# measurements from a TRUE state whose tap sits `steps` mechanical steps
# off the model position (fraction grid, the same grid the fixation uses)
function _measurements_with_true_tap(steps::Float64)
  tnet = _create_tap_net()
  br = tnet.branchVec[3]
  ct = Sparlectra._cascade_tap(br, steps * br.tap_step, 0.0, 0.0)
  br.tap_ratio = ct.tap_ratio
  br.phase_shift_deg = ct.phase_shift_deg
  ite, erg = runpf!(tnet, 40, 1e-12, 0; method = :rectangular)
  erg == 0 || error("power flow of the true state did not converge")
  return generateMeasurementsFromPF(tnet; includeImag = true, noise = false)
end

function run_state_estimation_tap_estimation()
  # 1) the problem: the true tap is 2 steps off, the model is neutral. A
  # plain estimation has nowhere to put the discrepancy: J explodes and
  # the suspicious measurements cluster at the transformer although every
  # telemetry row is healthy.
  meas = _measurements_with_true_tap(2.0)
  net = _create_tap_net()
  res0 = runse!(net, meas; maxIte = 40, tol = 1e-10, updateNet = false)
  @printf("without release : J = %10.3f (dof %d)  <- model discrepancy, NOT bad data\n", res0.objectiveJ, res0.dof)

  # 2) the fix: release the tap, estimate, and let the mandatory fixation
  # round it to the mechanical step. An integer deviation is recovered
  # exactly: J after the fixation is numerically zero again.
  net = _create_tap_net()
  setTapEstimation!(net; trafo = 3, mode = :ratio)
  res = runse!(net, meas; maxIte = 40, tol = 1e-10, updateNet = false)
  t = res.tapEstimates[1]
  tf = res.tapFixation
  @printf("with release    : electrical step %.3f -> fixed step %d\n", t.electrical_step_1, t.fixed_step_1)
  @printf("                  J before fixation %.3e (dof %d) -> after %.3e (dof %d)\n", tf.j_before, tf.dof_before, tf.j_after, tf.dof_after)

  # 3) the off-grid case: a HALF step cannot be a mechanical position. The
  # continuous estimate absorbs it, the fixation must round, and the
  # remaining J honestly says so (with noise the offgrid_residual note
  # marks the band flip as "not bad data").
  meas15 = _measurements_with_true_tap(1.5)
  net = _create_tap_net()
  setTapEstimation!(net; trafo = 3, mode = :ratio)
  res15 = runse!(net, meas15; maxIte = 40, tol = 1e-10, updateNet = false)
  t15 = res15.tapEstimates[1]
  tf15 = res15.tapFixation
  @printf("off-grid truth  : electrical step %.3f -> fixed step %d, J %.3e -> %.3f\n", t15.electrical_step_1, t15.fixed_step_1, tf15.j_before, tf15.j_after)

  # 4) write-back protection: the model is only touched on explicit
  # request (updateTaps), and only with the FIXED mechanical position
  net = _create_tap_net()
  setTapEstimation!(net; trafo = 3, mode = :ratio)
  before = net.branchVec[3].tap_ratio
  runse!(net, meas; maxIte = 40, tol = 1e-10)
  println("write-back      : without updateTaps the model tap stays bitwise untouched: ", net.branchVec[3].tap_ratio === before)
  runse!(net, meas; maxIte = 40, tol = 1e-10, updateTaps = true)
  @printf("                  with updateTaps = true the fixed position lands in the model: tap_ratio = %.6f\n", net.branchVec[3].tap_ratio)

  # 5) machine (GSU) transformers are the OTHER way around: never released
  # (their terminal voltage is AVR-set, not observed); after the SE their
  # tap is back-calculated from the AVR setpoint, the dispatch P, and the
  # MEASURED machine Q against the estimated network voltage.
  gnet = _create_tap_net()
  addBus!(net = gnet, busName = "G1", vn_kV = 10.0)
  addProsumer!(net = gnet, busName = "G1", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 5.0, vm_pu = 1.02)
  addPIModelTrafo!(net = gnet, fromBus = "H2", toBus = "G1", r_pu = 0.001, x_pu = 0.09, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  gnet.branchVec[5].has_ratio_tap = true
  gnet.branchVec[5].tap_step = 0.00625
  brg = gnet.branchVec[5]
  ctg = Sparlectra._cascade_tap(brg, 2 * brg.tap_step, 0.0, 0.0)   # true GSU position: 2 steps
  brg.tap_ratio = ctg.tap_ratio
  brg.phase_shift_deg = ctg.phase_shift_deg
  runpf!(gnet, 40, 1e-12, 0; method = :rectangular)
  V = buildVoltageVector(gnet)
  y = Sparlectra.calcAdmittance(brg, brg.comp.cVN, gnet.baseMVA)
  q_scada = imag(V[5] * conj(y[3] * V[2] + y[4] * V[5]) * gnet.baseMVA)   # the plant's Q telemetry
  append!(gnet.measurements, generateMeasurementsFromPF(gnet; includeImag = true, noise = false))
  runse!(gnet; maxIte = 40, tol = 1e-10, updateNet = true)
  bt = calcMachineTrafoTapFromSE(gnet; trafo = 5, v_machine_pu = 1.02, q_mvar = q_scada)
  @printf("machine trafo   : back-calculated electrical step %.3f -> fixed step %d (Q residual %.3f MVar)\n", bt.electrical_step, bt.fixed_step, bt.q_residual_mvar)
  return nothing
end

run_example(run_state_estimation_tap_estimation)

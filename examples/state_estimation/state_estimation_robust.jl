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
# file: examples/state_estimation/state_estimation_robust.jl
# purpose: SE phase 4 walkthrough. A 10 sigma gross error distorts the plain
#          WLS estimate; the two-stage robust R modification (robust = true)
#          suppresses it while the original-sigma diagnostics still name the
#          culprit first. The Wilson-Hilferty band test then shows both
#          failure directions: :high on the bad-data set and :low on a
#          noise-free synthetic set (sigmas overestimated).

using Sparlectra
using Random
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

function _create_demo_net()
  net = Net(name = "se_robust_demo", baseMVA = 100.0)
  for b in ("Slack", "LoadA", "LoadB")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addACLine!(net = net, fromBus = "Slack", toBus = "LoadA", length = 12.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "LoadA", toBus = "LoadB", length = 9.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "Slack", toBus = "LoadB", length = 11.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "LoadA", type = "ENERGYCONSUMER", p = 35.0, q = 12.0)
  addProsumer!(net = net, busName = "LoadB", type = "ENERGYCONSUMER", p = 28.0, q = 9.0)
  ok, msg = validate!(net = net)
  ok || error("demo net invalid: $msg")
  return net
end

function run_state_estimation_robust()
  print_example_banner("examples/state_estimation/state_estimation_robust.jl", "robust R modification suppresses a gross error; Wilson-Hilferty flags :high and :low")

  net = _create_demo_net()
  ite_pf, erg_pf = runpf!(net, 40, 1e-10, 0; method = :rectangular)
  erg_pf == 0 || error("Power flow did not converge")
  Vref = buildVoltageVector(net)
  @printf("PF iterations: %d\n\n", ite_pf)

  std = measurementStdDevs(vm = 0.002, pinj = 0.5, qinj = 0.5, pflow = 0.5, qflow = 0.5)
  meas = generateMeasurementsFromPF(net; noise = true, stddev = std, rng = MersenneTwister(42))
  badIdx = findfirst(m -> m.typ == Sparlectra.PflowMeas, meas)
  bm = meas[badIdx]
  meas[badIdx] = Measurement(typ = bm.typ, value = bm.value + 10.0 * bm.sigma, sigma = bm.sigma, active = bm.active, busIdx = bm.busIdx, branchIdx = bm.branchIdx, direction = bm.direction, id = bm.id)
  println("gross error: 10 sigma on ", bm.id, "\n")

  # 1) plain WLS versus robust: the modification pulls the estimate back to
  # the noise floor while the culprit reaches stage 2
  resPlain = runse!(net, meas; maxIte = 40, tol = 1e-10, updateNet = false)
  resRobust = runse!(net, meas; maxIte = 40, tol = 1e-10, updateNet = false, robust = true)
  @printf("max |V - V_true| : plain WLS = %.5f pu   robust = %.5f pu\n", maximum(abs.(resPlain.voltages .- Vref)), maximum(abs.(resRobust.voltages .- Vref)))
  rr = only(r for r in resRobust.robustRows if r.id == bm.id)
  @printf("robust stage of %s: stage %d, t = %.1f, sigma widened by %.2fx\n\n", rr.id, rr.stage, rr.t, rr.sigma_factor)

  # 2) solve/diagnosis separation: the original-sigma ranking still names the
  # culprit first (Wilson-Hilferty verdict :high)
  report = validate_measurements(net, meas; maxIte = 40, tol = 1e-8, robust = true)
  top = first(report.measurement_ranking)
  @printf("diagnosis (original sigmas): top suspect = %s (|rn| = %.1f), band test reason = %s\n\n", top.id, top.abs_normalized_residual, string(report.objective.reason))

  # 3) the other failure direction: a noise-free synthetic set is flagged
  # :low (J implausibly small, sigmas overestimated) instead of silently
  # passing like the old symmetric test
  measClean = generateMeasurementsFromPF(net; noise = false, stddev = std)
  repLow = validate_measurements(net, measClean; maxIte = 30, tol = 1e-8)
  @printf("noise-free synthetic set: reason = %s (z_wh = %.1f, nu = %d)\n", string(repLow.objective.reason), repLow.objective.z_wh, repLow.objective.dof)
  println("interpretation: ", summarize_se_diagnostics(repLow).reason)
  return nothing
end

run_example(run_state_estimation_robust)

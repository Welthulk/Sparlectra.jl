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
# file: examples/state_estimation/state_estimation_imag_bad_data.jl
# purpose: current-magnitude measurements (ImagMeas) and bad-data localization
#          stage 1 (SE phase 1): shows that added branch currents raise the
#          residual sensitivity wii of the power measurements without touching
#          observability, then locates a gross error via sequential elimination
#          and prints the trace. Also demonstrates the ImagMeas value gate.

using Sparlectra
using Random
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

function _create_demo_net()
  net = Net(name = "se_imag_demo", baseMVA = 100.0)

  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "LoadA", vn_kV = 110.0)
  addBus!(net = net, busName = "LoadB", vn_kV = 110.0)

  addACLine!(net = net, fromBus = "Slack", toBus = "LoadA", length = 12.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "LoadA", toBus = "LoadB", length = 9.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "Slack", toBus = "LoadB", length = 11.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "LoadA", type = "ENERGYCONSUMER", p = 35.0, q = 12.0)
  addProsumer!(net = net, busName = "LoadB", type = "ENERGYCONSUMER", p = 28.0, q = 9.0)

  return net
end

# inject a gross error of n_sigma standard deviations into measurement idx
function _inject_gross_error!(meas, idx, n_sigma)
  m = meas[idx]
  meas[idx] = Measurement(typ = m.typ, value = m.value + n_sigma * m.sigma, sigma = m.sigma, active = m.active, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id)
  return m.id
end

function run_state_estimation_imag_bad_data()
  print_example_banner("examples/state_estimation/state_estimation_imag_bad_data.jl", "ImagMeas currents raise bad-data localizability (wii); sequential elimination locates a gross error")
  net = _create_demo_net()

  ite_pf, erg_pf = runpf!(net, 40, 1e-10, 0; method = :rectangular)
  erg_pf == 0 || error("Power flow did not converge")
  @printf("PF iterations: %d\n\n", ite_pf)

  std = measurementStdDevs(vm = 0.005, pinj = 0.5, qinj = 0.5, pflow = 0.5, qflow = 0.5, imag = 2.0)

  # the same synthetic measurement set twice: without and with branch currents
  _fresh(withImag) = generateMeasurementsFromPF(net; includeImag = withImag, noise = true, stddev = std, rng = MersenneTwister(1234))

  # 1) currents do not carry observability: identical verdict
  obs_without = evaluate_global_observability(net, _fresh(false))
  obs_with = evaluate_global_observability(net, _fresh(true))
  println("observability : rank ", obs_without.numerical_rank, " / quality ", obs_without.quality, " without currents")
  println("                rank ", obs_with.numerical_rank, " / quality ", obs_with.quality, " with currents (ImagMeas rows excluded by design)")
  println()

  # 2) localizability: the wii of the faulted power measurement with and
  # without current measurements (current measurements raise wii)
  println("gross error: 10 sigma on the first Pflow measurement")
  for withImag in (false, true)
    meas = _fresh(withImag)
    bad_idx = findfirst(m -> m.typ == Sparlectra.PflowMeas, meas)
    bad_id = _inject_gross_error!(meas, bad_idx, 10.0)
    diag = runse_diagnostics(net, meas; max_eliminations = 3, maxIte = 20, tol = 1e-8)
    row = only(r for r in diag.diagnostics.measurement_ranking if r.measurement_index == bad_idx)
    @printf("%-18s wii(%s) = %.3f  localizable = %-5s  eliminated first: %s  stop: %s\n",
      withImag ? "with currents:" : "without currents:", bad_id, row.wii, string(row.localizable),
      string(diag.eliminations[1].id == bad_id), string(diag.stop_reason))
    if withImag
      println("\nsequential elimination trace (with currents):")
      print_se_diagnostics(diag; topN = 5)
    end
  end
  println()

  # 3) the value gate: a current reading below 3 sigma never enters the
  # estimation (derivative of |I| is discontinuous near zero current)
  meas = _fresh(true)
  addImagMeasurement!(meas; net = net, value = 1.0, sigma = 10.0, fromBus = "Slack", toBus = "LoadA", id = "Imag_low_current")
  res = runse!(net, meas; maxIte = 20, tol = 1e-8, updateNet = false)
  println("value gate    : added Imag_low_current with value 1 A < 3 sigma = 30 A;")
  println("                the estimator logs its exclusion and converged = ", res.converged)
  return nothing
end

run_example(run_state_estimation_imag_bad_data)

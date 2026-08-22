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

# Date: 2026-07-24
# file: examples/state_estimation/mc_state_estimation_study.jl
# purpose: Monte-Carlo state-estimation error study on the 7-bus workshop net —
#          repeated noisy measurement sets, WLS estimation, error statistics

using Sparlectra
using Random
using Statistics
using Printf

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# --- Study parameters --------------------------------------------------------
const M_RUNS = 500              # Monte-Carlo runs (one noisy measurement set each)
const SEED_BASE = 20260700      # run k uses MersenneTwister(SEED_BASE + k)
const SE_MAX_ITER = 12
const SE_TOL = 1e-6

"""Build the 7-bus ring network from the workshop tour (docs/lit/workshop_tour.jl, chapter 1)."""
function buildWorkshop7BusNet()::Net
  net = Net(name = "mc_se_7bus", baseMVA = 100.0)

  addBus!(net = net, busName = "B1", vn_kV = 110.0, vm_pu = 1.02, va_deg = 0.0)
  for i in 2:7
    addBus!(net = net, busName = "B$(i)", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  end

  # Ring + cross-connections
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B4", toBus = "B5", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B5", toBus = "B6", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B6", toBus = "B7", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B7", toBus = "B1", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B5", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B6", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)

  # Source / generation / loads
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)
  addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 60.0, q = 10.0)
  addProsumer!(net = net, busName = "B2", type = "LOAD", p = 35.0, q = 10.0)
  addProsumer!(net = net, busName = "B4", type = "LOAD", p = 45.0, q = 15.0)
  addProsumer!(net = net, busName = "B5", type = "LOAD", p = 25.0, q = 8.0)
  addProsumer!(net = net, busName = "B6", type = "LOAD", p = 30.0, q = 10.0)
  addProsumer!(net = net, busName = "B7", type = "LOAD", p = 20.0, q = 6.0)

  ok, msg = validate!(net = net)
  ok || error("Validation failed: $(msg)")
  return net
end

function main()
  print_example_banner("examples/state_estimation/mc_state_estimation_study.jl", "Monte-Carlo WLS state-estimation error study on the 7-bus workshop network")

  net = buildWorkshop7BusNet()
  nb = length(net.nodeVec)

  # Reference power flow: the "true" system state the estimator should recover.
  ite_pf, status_pf = runpf!(net, 40, 1e-10, 0)
  status_pf == 0 || error("Reference power flow did not converge.")
  vm_ref = [node._vm_pu for node in net.nodeVec]
  va_ref = [node._va_deg for node in net.nodeVec]
  println("Reference power flow solved in $(ite_pf) iterations.")

  # Measurement noise model (same standard deviations for every run).
  noise_std = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7)

  # One-time observability check on a noise-free measurement set.
  setMeasurementsFromPF!(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false, stddev = noise_std)
  gobs = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
  println("Global observability: quality=$(gobs.quality), measurements=$(gobs.n_measurements), states=$(gobs.n_states)")

  vm_rmse = Float64[]
  va_rmse = Float64[]
  vm_maxerr = Float64[]
  va_maxerr = Float64[]
  iterations = Int[]
  n_converged = 0
  n_failed = 0

  for k in 1:M_RUNS
    # Fresh noisy measurement set from the (unchanged) reference state.
    setMeasurementsFromPF!(
      net;
      includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true,
      noise = true, stddev = noise_std, rng = MersenneTwister(SEED_BASE + k),
    )

    se = runse!(net; maxIte = SE_MAX_ITER, tol = SE_TOL, flatstart = true, jacEps = 1e-6, updateNet = true)
    if !se.converged
      n_failed += 1
      continue
    end
    n_converged += 1
    push!(iterations, se.iterations)

    vm_err = [net.nodeVec[b]._vm_pu - vm_ref[b] for b in 1:nb]
    va_err = [net.nodeVec[b]._va_deg - va_ref[b] for b in 1:nb]
    push!(vm_rmse, sqrt(mean(abs2, vm_err)))
    push!(va_rmse, sqrt(mean(abs2, va_err)))
    push!(vm_maxerr, maximum(abs, vm_err))
    push!(va_maxerr, maximum(abs, va_err))
  end

  # Restore the reference state on the net (the last run left its estimate).
  for b in 1:nb
    net.nodeVec[b]._vm_pu = vm_ref[b]
    net.nodeVec[b]._va_deg = va_ref[b]
  end

  n_converged > 0 || error("No state-estimation run converged.")

  println()
  println("=== Monte-Carlo state-estimation study (7-bus workshop net) ===")
  @printf("runs             : %d\n", M_RUNS)
  @printf("converged        : %d (%.1f %%)\n", n_converged, 100.0 * n_converged / M_RUNS)
  @printf("failed           : %d\n", n_failed)
  @printf("mean iterations  : %.2f\n", mean(iterations))
  @printf("seed base        : %d\n", SEED_BASE)
  println()
  @printf("%-22s %12s %12s %12s\n", "metric", "mean", "max", "std")
  println("-"^62)
  @printf("%-22s %12.3e %12.3e %12.3e\n", "Vm RMSE [pu]", mean(vm_rmse), maximum(vm_rmse), std(vm_rmse))
  @printf("%-22s %12.3e %12.3e %12.3e\n", "Vm max error [pu]", mean(vm_maxerr), maximum(vm_maxerr), std(vm_maxerr))
  @printf("%-22s %12.3e %12.3e %12.3e\n", "Va RMSE [deg]", mean(va_rmse), maximum(va_rmse), std(va_rmse))
  @printf("%-22s %12.3e %12.3e %12.3e\n", "Va max error [deg]", mean(va_maxerr), maximum(va_maxerr), std(va_maxerr))

  # Note on parallelization: runs are independent, so this loop could run on
  # Threads.@threads with one deepcopy(net) per thread (never share a Net
  # across threads). Kept serial here for clarity and determinism.
  return nothing
end

main()

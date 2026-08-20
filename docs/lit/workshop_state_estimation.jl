# Copyright 2023-2026 Udo Schmitz                                             #src
#                                                                             #src
# Licensed under the Apache License, Version 2.0 (the "License");             #src
# you may not use this file except in compliance with the License.            #src
# You may obtain a copy of the License at                                     #src
#                                                                             #src
#     http://www.apache.org/licenses/LICENSE-2.0                              #src
#                                                                             #src
# Unless required by applicable law or agreed to in writing, software         #src
# distributed under the License is distributed on an "AS IS" BASIS,           #src
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    #src
# See the License for the specific language governing permissions and         #src
# limitations under the License.                                              #src
#                                                                             #src
# file: docs/lit/workshop_state_estimation.jl                                 #src
# purpose: Literate.jl source of the state-estimation workshop notebook and   #src
#          its Documenter page. Follow-up to workshop_intro.jl; regenerate    #src
#          the committed outputs with                                         #src
#          `julia --project=docs docs/generate_notebooks.jl`.                 #src

# # State estimation from noisy measurements
#
# [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb)
#
# A power flow computes the network state from an exact specification. A
# real control room has the opposite problem: it receives many **redundant,
# noisy measurements** (voltage magnitudes, injections, branch flows) and
# must reconstruct the most likely state from them. That is *state
# estimation*: a weighted-least-squares (WLS) fit of the bus voltages to
# the measurement set, where each measurement counts with the inverse of
# its variance.
#
# In this notebook you build a 7-bus network with
# [Sparlectra.jl](https://github.com/Welthulk/Sparlectra.jl), solve a
# reference power flow, derive a noisy synthetic measurement set from it,
# check observability, and run the estimator: the closed loop that lets
# you judge estimation quality against a known truth.
#
# The network is the workshop ring (B1 carries the grid connection, the
# diagonals are the cross-ties B2-B5 and B3-B6):
#
# ```text
#  (slack)
#    B1 ---- B2 ---- B3 ---- B4
#    |         \    /         |
#    |          \  /          |
#    |           \/           |
#    |           /\           |
#    B7 ---- B6 ---- B5 ------+
# ```
#
# > **Note:** On Google Colab the install cell takes a few minutes on a
# > fresh session (package download and precompilation). Colab's Julia
# > version may change over time; this notebook targets Julia ≥ 1.12.

#nb # ## Setup (Colab)
#nb # This cell installs Sparlectra from GitHub (branch `main`) into a fresh
#nb # temporary environment. The isolation matters on Colab: the shared
#nb # default environment ships many preinstalled packages, and installing
#nb # anything there triggers precompilation of that whole stack. Run this
#nb # cell first, once per session; it takes a few minutes.
#nb using Pkg
#nb Pkg.activate(temp = true)
#nb Pkg.add(url = "https://github.com/Welthulk/Sparlectra.jl", rev = "main")
#nb ## To test another branch, set rev to its name, e.g. rev = "dev/r0.9.8".
#nb ## For the latest registered release use: Pkg.add("Sparlectra")
#nb ## Switching versions in a running session? A "[loaded: ...]" note means
#nb ## the old version is still active — restart the runtime, then rerun
#nb ## this cell.

# ## Load the packages
#
# `Random` (standard library) provides the seeded generator that makes the
# synthetic measurement noise reproducible.

using Sparlectra
using Random

# ## Warm-up
#
# Julia compiles each function on first use. This cell warms the two paths
# the notebook exercises, the power-flow solver (which produces the
# reference state) and the WLS estimator, on a tiny throwaway network, so
# the real study runs at full speed.

wnet = Net(name = "warmup", baseMVA = 100.0)
addBus!(net = wnet, busName = "A", vn_kV = 110.0)
addBus!(net = wnet, busName = "B", vn_kV = 110.0)
addProsumer!(net = wnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
addPIModelACLine!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
t_pf = @elapsed runpf!(wnet, 10, 1e-8, 0)
setMeasurementsFromPF!(wnet; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false)
t_se = @elapsed runse!(wnet; maxIte = 8, tol = 1e-6, flatstart = true, jacEps = 1e-6, updateNet = false)
println("warm: power flow ", round(t_pf; digits = 2), " s, estimator ", round(t_se; digits = 2), " s (first calls compile)")

# ## Build the study network
#
# Seven 110 kV buses in a ring with two cross-connections, enough meshing
# that branch-flow measurements carry real information. An external network
# injection at `B1` is the slack, a generator feeds at `B3`, the other
# buses carry loads.

net = Net(name = "workshop_se_7bus", baseMVA = 100.0)

addBus!(net = net, busName = "B1", vn_kV = 110.0, vm_pu = 1.02, va_deg = 0.0)
for i in 2:7
  addBus!(net = net, busName = "B$(i)", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
end

addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B4", toBus = "B5", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B5", toBus = "B6", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B6", toBus = "B7", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B7", toBus = "B1", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B2", toBus = "B5", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B3", toBus = "B6", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)

addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)
addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 60.0, q = 10.0)
addProsumer!(net = net, busName = "B2", type = "LOAD", p = 35.0, q = 10.0)
addProsumer!(net = net, busName = "B4", type = "LOAD", p = 45.0, q = 15.0)
addProsumer!(net = net, busName = "B5", type = "LOAD", p = 25.0, q = 8.0)
addProsumer!(net = net, busName = "B6", type = "LOAD", p = 30.0, q = 10.0)
addProsumer!(net = net, busName = "B7", type = "LOAD", p = 20.0, q = 6.0)

ok, msg = validate!(net = net)
ok || error("Network validation failed: $msg")

# ## Solve the reference power flow
#
# This solved state plays the role of the (in reality unknown) *true*
# system state: the measurements are derived from it, and the estimator
# never sees it directly.

ite_pf, status_pf = runpf!(net, 40, 1e-10, 0)
status_pf == 0 || error("Power flow did not converge")
println("reference power flow converged in $ite_pf iterations")

# ## Derive a noisy measurement set
#
# `measurementStdDevs` defines one standard deviation per measurement kind
# (pu for voltages, MW/MVar for powers); `setMeasurementsFromPF!` then
# reads the solved state and creates the measurements, with Gaussian noise
# of exactly those standard deviations, seeded for reproducibility. This is
# the standard trick for exercising an estimator: the truth is known, so
# the residuals are meaningful.

std = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7)
setMeasurementsFromPF!(
  net;
  includeVm = true,
  includePinj = true,
  includeQinj = true,
  includePflow = true,
  includeQflow = true,
  noise = true,
  stddev = std,
  rng = MersenneTwister(42),
)
println(length(net.measurements), " measurements created")

# ## Check observability
#
# Estimation only works where the measurement set actually determines the
# state. The global check compares measurement count against state count
# and probes the numerical rank of the measurement Jacobian.

gobs = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
println("global observability quality: ", gobs.quality)
println("measurements: ", gobs.n_measurements, ", states: ", gobs.n_states)

# ## Run the estimator
#
# `runse!` iterates the WLS normal equations from a flat start. With
# `updateNet = true` the estimated voltages are written back into the
# network, so the usual result printers show the *estimated* state. The
# objective $J$ is the weighted sum of squared residuals; for healthy
# Gaussian noise it should land near the number of redundant measurements.

se = runse!(
  net;
  maxIte = 12,
  tol = 1e-6,
  flatstart = true,
  jacEps = 1e-6,
  updateNet = true,
)

println("SE converged: ", se.converged, " in ", se.iterations, " iterations")
println("objective J:  ", round(se.objectiveJ; digits = 2))

# ## Inspect the estimated state
#
# `SEResult.voltages` holds the estimated complex bus voltage per bus
# index. Because the measurements carried only mild noise, the estimate
# reproduces the reference power flow closely (magnitudes to a few
# 1e-3 pu), and the slack reference `B1` comes back at 1.02 pu / 0°. The
# degrees of freedom (`dof`) and the 3σ check on $J$ summarize whether the
# residuals are consistent with the declared measurement accuracies.

println("dof (redundant measurements): ", se.dof)
println("J within 3σ band:             ", se.jWithin3Sigma)
println()
for (name, idx) in sort(collect(net.busDict); by = last)
  v = se.voltages[idx]
  println(rpad(name, 4), "  Vm = ", round(abs(v); digits = 4), " pu   Va = ", round(rad2deg(angle(v)); digits = 3), "°")
end

# ## Building measurement sets manually
#
# Real measurement sets are not derived from a solved power flow; they
# arrive from SCADA one by one. The `add*Measurement!` helpers resolve bus
# and branch references for you, so assembling a set works just like
# building the network. A deliberately sparse set like this one is a good
# way to explore *when observability breaks down*:

empty!(net.measurements)

addVmMeasurement!(net; busName = "B1", value = 1.02, sigma = 0.002)
addPinjMeasurement!(net; busName = "B2", value = -35.0, sigma = 1.0)
addQinjMeasurement!(net; busName = "B2", value = -10.0, sigma = 1.0)
addPflowMeasurement!(net; fromBus = "B1", toBus = "B2", value = 22.0, sigma = 0.8, direction = :from)
addQflowMeasurement!(net; branchNr = 1, value = 7.0, sigma = 0.8, direction = :to)

obs = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
println("sparse manual set quality: ", obs.quality)
println("measurements: ", obs.n_measurements, ", states: ", obs.n_states)

# Five measurements against 13 states: the check reports the shortfall
# instead of letting the estimator run into a rank-deficient system.
#
# ## Where to go next
#
# - New to Sparlectra? The
#   [introduction notebook](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_intro.ipynb)
#   builds a network from scratch step by step, directly in Colab.
# - [Slack types and short-circuit currents](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb):
#   the grid-connection notebook: slack representations and IEC 60909-0
#   fault currents.
# - [State Estimation](https://welthulk.github.io/Sparlectra.jl/state_estimation/):
#   the full documentation: observability analysis, PMU angle measurements,
#   bad-data handling, and the measurement model.
# - [State-Estimation Configuration](https://welthulk.github.io/Sparlectra.jl/state_estimation_configuration/):
#   every `state_estimation.*` configuration key.

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
# file: docs/lit/workshop_intro.jl                                            #src
# purpose: Literate.jl source of the introductory workshop notebook and the   #src
#          "Try it in your Browser" Documenter page. One maintained file      #src
#          generates both outputs; regenerate them with                       #src
#          `julia --project=docs docs/generate_notebooks.jl`. Lines prefixed  #src
#          with `#nb` appear only in the notebook (Colab install cell);       #src
#          lines suffixed with `#src` appear in neither output.               #src

# # Your first power flow with Sparlectra
#
# [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_intro.ipynb)
#
# [Sparlectra.jl](https://github.com/Welthulk/Sparlectra.jl) is a Julia
# framework for AC power-flow and state-estimation studies. In this guided
# first tour you build a small 110 kV network from scratch (seven buses in a
# ring with two cross-connections), validate it, solve it with the built-in
# Newton-Raphson solver, and read the classical result tables. No input files
# are needed; everything is created programmatically.
#
# The network at a glance (B1 carries the grid connection, the diagonals
# are the two cross-ties B2-B5 and B3-B6):
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
# > **Note:** When running this notebook on Google Colab, the install cell
# > takes a few minutes on a fresh session (package download and
# > precompilation). Colab's Julia version may change over time; this
# > notebook targets Julia ≥ 1.12.

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

# ## Load the package

using Sparlectra

# ## Create an empty network
#
# Every Sparlectra model starts from a [`Net`](https://welthulk.github.io/Sparlectra.jl/reference_network/)
# object. `baseMVA` is the system base power used for all per-unit
# conversions.

net = Net(name = "workshop_intro", baseMVA = 100.0)

# ## Add buses
#
# `addBus!` creates the electrical nodes. `vn_kV` is the nominal voltage,
# `vm_pu` and `va_deg` provide the starting voltage for the solver. Note that
# the operational bus type (slack / PV / PQ) is *not* declared here: it is
# derived later from the devices attached to each bus.

addBus!(net = net, busName = "B1", vn_kV = 110.0, vm_pu = 1.02, va_deg = 0.0)
for i in 2:7
  addBus!(net = net, busName = "B$(i)", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
end

# ## Connect the buses
#
# `addPIModelACLine!` adds a line as a π-equivalent branch with per-unit
# parameters (`r_pu` resistance, `x_pu` reactance, `b_pu` total line
# charging). The seven buses form a ring, and two extra cross-connections
# give the power flow alternative paths.

addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B4", toBus = "B5", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B5", toBus = "B6", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B6", toBus = "B7", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B7", toBus = "B1", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B2", toBus = "B5", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B3", toBus = "B6", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)

# ## Attach loads and generation
#
# Devices that consume or produce power are added with `addProsumer!`. The
# external network injection at `B1` references its own bus as the voltage
# reference: that makes `B1` the slack bus, which balances the network and
# fixes voltage magnitude and angle. The generator at `B3` feeds in 60 MW,
# and the remaining buses carry loads (`p` in MW, `q` in MVar).

addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)
addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 60.0, q = 10.0)
addProsumer!(net = net, busName = "B2", type = "LOAD", p = 35.0, q = 10.0)
addProsumer!(net = net, busName = "B4", type = "LOAD", p = 45.0, q = 15.0)
addProsumer!(net = net, busName = "B5", type = "LOAD", p = 25.0, q = 8.0)
addProsumer!(net = net, busName = "B6", type = "LOAD", p = 30.0, q = 10.0)
addProsumer!(net = net, busName = "B7", type = "LOAD", p = 20.0, q = 6.0)

# ## Validate the network
#
# `validate!` checks the model for structural problems (unconnected buses,
# missing slack, inconsistent parameters) before any numerics run. Make this
# a habit after every round of model edits.

ok, msg = validate!(net = net)
ok || error("Network validation failed: $msg")

# ## Solve the power flow
#
# `runpf!` runs the rectangular Newton-Raphson solver. The arguments are the
# maximum number of iterations, the convergence tolerance, and a verbosity
# level; it returns the iteration count and a status (`0` means converged).

tol    = 1e-8
maxIte = 10

etime = @elapsed begin
  ite, erg = runpf!(net, maxIte, tol, 0)
end

# ## Inspect the results
#
# `calcNetLosses!` computes the branch flows and total network losses from
# the converged voltages; `printACPFlowResults` prints the classical result
# tables: bus voltages, branch flows, and losses.

if erg == 0
  calcNetLosses!(net)
  printACPFlowResults(net, etime, ite, tol)
else
  @warn "Power flow did not converge (status = $erg)"
end

# ## Where to go next
#
# You have built, validated, and solved a complete network in a few dozen
# lines. From here:
#
# - [Slack types and short-circuit currents](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb):
#   the follow-up notebook, directly in Colab: ideal slack vs. external-grid
#   source vs. distributed slack, plus IEC 60909-0 fault currents.
# - [Workshop](https://welthulk.github.io/Sparlectra.jl/workshop/): file
#   import and export, transformers, bus links, tap control, and Q-limits.
# - [State estimation from noisy measurements](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb):
#   the estimation notebook, directly in Colab: reconstruct the network
#   state from redundant, noisy measurements instead of a load
#   specification.
# - [Feature Matrix](https://welthulk.github.io/Sparlectra.jl/feature_matrix/):
#   what Sparlectra covers, at a glance.

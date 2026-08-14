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
# file: docs/lit/workshop_series_compensation.jl                              #src
# purpose: Literate.jl source of the series-compensation (TCSC) workshop      #src
#          notebook and its Documenter page: why a series reactance steers    #src
#          flow, and how the SeriesReactanceControl outer-loop controller     #src
#          moves a loop-network flow split onto a branch target. Regenerate   #src
#          the committed outputs with                                         #src
#          `julia --project=docs docs/generate_notebooks.jl`.                 #src

# # Series compensation: steering flow with a TCSC
#
# [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_series_compensation.ipynb)
#
# In a meshed AC network, power does not follow contracts, it follows
# impedance: parallel paths split the transfer in inverse proportion to
# their reactances. A TCSC (thyristor controlled series capacitor) exploits
# exactly that lever: a variable series reactance in one line steers how
# much flow that corridor carries. In this notebook you build a small loop
# network with [Sparlectra.jl](https://github.com/Welthulk/Sparlectra.jl),
# watch the natural flow split, and then let the
# `SeriesReactanceControl` outer-loop controller move the split onto a
# target, including the honest failure mode when the target is out of
# reach.
#
# > **Note:** On Google Colab the install cell takes a few minutes on a
# > fresh session (package download and precompilation). Colab's Julia
# > version may change over time; this notebook targets Julia ≥ 1.12.

#nb # ## Setup (Colab)
#nb # Installing Sparlectra takes a few minutes on a fresh Colab
#nb # session (package download + precompilation). Run this cell first.
#nb using Pkg
#nb Pkg.add("Sparlectra")

# ## Load the package

using Sparlectra

# ## Why a series reactance steers flow
#
# Every branch enters the power flow through its admittance matrix (see the
# [Branch Model](https://welthulk.github.io/Sparlectra.jl/branchmodel/)
# page for the derivation):
#
# ```math
# Y_{br} = \begin{bmatrix}
#   \frac{1}{\tau^2}\left(y_{ser} + \frac{y_{shunt}}{2}\right) &
#   -y_{ser}\,\frac{1}{\tau e^{-j\phi}} \\
#   -y_{ser}\,\frac{1}{\tau e^{j\phi}} &
#   y_{ser} + \frac{y_{shunt}}{2}
# \end{bmatrix},
# \qquad y_{ser} = \frac{1}{R + jX},
# ```
#
# with $N = 1$ for lines. The TCSC acts purely through $X$ inside
# $y_{ser}$: every accepted controller step changes one branch stamp and
# the outer loop re-stamps the Y-bus before the next solve.
#
# For a lossless line the transfer relation
#
# ```math
# P_{12} = \frac{V_1 V_2}{X}\,\sin(\delta_1 - \delta_2)
# ```
#
# says that a lower series reactance carries more power at a given angle
# difference. In a loop, flow redistributes between the parallel paths
# according to their reactance ratio, which is exactly what we are about
# to watch.

# ## A loop network with two corridors
#
# 80 MW travel from source `A` to sink `B` over two parallel corridors.
# The upper corridor (`A` to `M2` to `B`) has twice the reactance of the
# lower one, so it naturally carries only one third of the transfer.

function build_loop()
  net = Net(name = "tcsc_workshop", baseMVA = 100.0)
  for b in ("A", "M1", "M2", "B")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 80.0, q = 20.0)
  addPIModelACLine!(net = net, fromBus = "A", toBus = "M1", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "M1", toBus = "B", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "A", toBus = "M2", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "M2", toBus = "B", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
  ok, msg = validate!(net = net)
  ok || error("Network validation failed: $msg")
  return net
end

net = build_loop()
run_sparlectra(net = net)
println("natural split: corridor 1 (A->M1) = ", round(get_branch_p_from_to_mw(net, "A", "M1"); digits = 2), " MW")
println("               corridor 2 (A->M2) = ", round(get_branch_p_from_to_mw(net, "A", "M2"); digits = 2), " MW")

# The 2:1 reactance ratio produces the 2:1 flow split, independent of any
# thermal ratings: the low-reactance corridor attracts the flow.
#
# ## Attach the TCSC and steer the split
#
# `addSeriesReactanceControl!` registers the controller on the line from
# `A` to `M2`. The target of 35 MW needs a visible reactance move: the
# outer loop measures the branch flow after each converged solve, steps
# `x_pu` via secant iteration (the first step is a bounded probe, because
# the sign of $dP/dX$ depends on the network), and stops inside the
# 0.5 MW default deadband.

ctrl = addSeriesReactanceControl!(net; fromBus = "A", toBus = "M2", p_target_mw = 35.0, x_min_pu = 0.02, x_max_pu = 0.30)
run_sparlectra(net = net)
println("steered:  corridor 2 (A->M2) = ", round(get_branch_p_from_to_mw(net, "A", "M2"); digits = 2), " MW (target 35)")
println("          x_pu moved 0.20 -> ", round(ctrl.x_pu; digits = 4), ", status = ", ctrl.status)

# The controller row from the last control run and the generic
# controllable-element view carry the shared vocabulary (actuator, range,
# quantity, target) that all outer controllers report:

cr = latest_control_result(net)
println("outer loop: status = ", cr.status, ", outer iterations = ", cr.outer_iterations, ", pf solves = ", cr.powerflow_solves)
for row in cr.controllers
  println("controller: ", row.controller_name, " achieved ", round(row.achieved_p_mw; digits = 2), " MW of ", row.p_target_mw, " MW, x_pu = ", round(row.x_pu; digits = 4))
end
for e in controllableElements(net)
  println("element:    ", e.element, " | ", e.device, " | ", e.actuator, " in [", e.actuator_min, ", ", e.actuator_max, "] | ", e.quantity, " @ ", e.target)
end

# ## The honest limit
#
# Ask the same corridor for 70 MW and the range `[0.02, 0.30]` is not
# enough: the reactance clamps at the capacitive end, the branch behaves
# as a fixed compensated line, and the controller reports `at_limit`
# instead of pretending convergence. The power flow itself stays valid.

net2 = build_loop()
ctrl2 = addSeriesReactanceControl!(net2; fromBus = "A", toBus = "M2", p_target_mw = 70.0, x_min_pu = 0.02, x_max_pu = 0.30)
run_sparlectra(net = net2)
println("limited: corridor 2 (A->M2) = ", round(get_branch_p_from_to_mw(net2, "A", "M2"); digits = 2), " MW (target 70)")
println("         x_pu = ", round(ctrl2.x_pu; digits = 4), " (clamped), at_limit = ", ctrl2.at_limit, ", converged = ", ctrl2.converged)

# ## Where to go next
#
# - [Series Compensation (TCSC)](https://welthulk.github.io/Sparlectra.jl/series_compensation/):
#   the theory page, compensation degree, device versus model, and the
#   resonance guard.
# - [Control Framework](https://welthulk.github.io/Sparlectra.jl/control_framework/):
#   the outer loop all controllers share, and the uniform element view.
# - [Workshop tour](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb):
#   all workshop examples in one Colab session.

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
# See the License for the specific language governing permissions and        #src
# limitations under the License.                                              #src
#                                                                             #src
# file: docs/lit/workshop_tour.jl                                             #src
# purpose: Literate.jl source of the all-in-one workshop tour notebook and    #src
#          its Documenter page: one install, one warm-up, then compact        #src
#          chapters for power flow, slack types and short circuit, tap        #src
#          control, Q(U) control, remote voltage control, and state           #src
#          estimation. Regenerate the committed outputs with                  #src
#          `julia --project=docs docs/generate_notebooks.jl`.                 #src

# # The Sparlectra workshop tour
#
# [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb)
#
# All workshop examples in **one session**: install
# [Sparlectra.jl](https://github.com/Welthulk/Sparlectra.jl) once, warm the
# compiler up once, and then walk through the chapters without waiting
# again. The focused single-topic notebooks cover the same ground with more
# narrative; this tour is the fast lane.
#
# After the warm-up (compilation happens there, everything after is fast)
# the chapters are:
#
# 1. A first power flow
# 2. Slack types and short-circuit currents
# 3. Transformer tap control (OLTC)
# 4. Voltage-dependent reactive power, Q(U)
# 5. Remote voltage control by a machine
# 6. A steerable HVDC link (back-to-back pairing)
# 7. State estimation
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

# ## Warm-up
#
# Julia compiles each function on first use. The tiny network below triggers
# that compilation for the solver path once, so every later chapter runs at
# full speed. The two timings make the effect visible.

using Sparlectra

wnet = Net(name = "warmup", baseMVA = 100.0)
addBus!(net = wnet, busName = "A", vn_kV = 110.0)
addBus!(net = wnet, busName = "B", vn_kV = 110.0)
addProsumer!(net = wnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
addPIModelACLine!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)

t_first = @elapsed runpf!(wnet, 10, 1e-8, 0)
t_second = @elapsed runpf!(wnet, 10, 1e-8, 0)
println("first solve (includes compilation): ", round(t_first; digits = 2), " s")
println("second solve (compiled):            ", round(t_second * 1000; digits = 2), " ms")

# A solve helper used by all chapters (25 iterations maximum, tolerance
# $10^{-8}$):

function solve!(net; kwargs...)
  etime = @elapsed begin
    ite, erg = runpf!(net, 25, 1e-8, 0; kwargs...)
  end
  erg == 0 || error("Power flow did not converge (status = $erg)")
  calcNetLosses!(net)
  return etime, ite
end

# ## Chapter 1: a first power flow
#
# Seven 110 kV buses in a ring with two cross-connections, an external
# network injection as slack at `B1`, a generator at `B3`, loads elsewhere.
# The builder is a function because chapter 7 reuses the same network.
# The full guided version of this chapter is the
# [introduction notebook](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_intro.ipynb).
#
# ```text
#  (slack)
#    B1 ---- B2 ---- B3 ---- B4
#    |         \    /         |
#    |          \  /          |
#    |           \/           |     diagonals: B2-B5 and B3-B6
#    |           /\           |
#    B7 ---- B6 ---- B5 ------+
# ```

function build_ring7(name::String)
  net = Net(name = name, baseMVA = 100.0)
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
  return net
end

net1 = build_ring7("tour_first_pf")
etime, ite = solve!(net1)
printACPFlowResults(net1, etime, ite, 1e-8)

# Reading aid: the slack at `B1` covers the difference between 155 MW of
# load, 60 MW of scheduled generation, and the network losses. All bus
# voltages stay near 1.0 pu.

# ### How much can you trust these numbers?
#
# Every Newton iteration solves the linear system $J \, \Delta x = -F$ with
# the power-flow Jacobian $J$. The condition number $\kappa(J)$ measures how
# strongly that solve amplifies tiny perturbations: rounding, measurement
# noise in the input data, small parameter changes. The attainable relative
# accuracy in Float64 is roughly $\kappa \cdot 2 \cdot 10^{-16}$, so every
# power of ten in $\kappa$ costs one significant digit of the result.
# `condestJacobian(net)` estimates $\kappa_1$ at the operating point the net
# currently holds, on the same sparse Jacobian the solver factors:

println("ring network: kappa = ", round(condestJacobian(net1), sigdigits = 3))

# Around 45: excellent. Rule of thumb: below about $10^6$ well conditioned,
# around $10^{10}$ borderline, beyond $10^{14}$ numerically singular in
# Float64.
#
# The instructive part is how conditioning degrades when the physics
# degenerate, long before the solver visibly fails. Take a small feeder with
# a measurement stub at `B3` and make the stub line weaker in each round:

for x_weak in (0.08, 800.0, 8.0e6, 8.0e10)
  net = Net(name = "tour_cond", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 20.0, q = 5.0)
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = x_weak / 8, x_pu = x_weak, b_pu = 0.0, status = 1)
  etime, ite = solve!(net)
  vm3 = round(something(net.nodeVec[3]._vm_pu, NaN); digits = 4)
  println("x = ", lpad(x_weak, 8), " pu:  ", ite, " iterations,  Vm(B3) = ", vm3, ",  kappa = ", round(condestJacobian(net), sigdigits = 3))
end

# Reading aid: the solver converges in 4 iterations in every round, and
# `Vm(B3)` prints the same plausible 0.9938 each time. Nothing in the result
# table reveals that the last round sits at $\kappa \approx 10^{12}$, where
# only about 4 significant digits survive: the printed voltage is already at
# the edge of what the arithmetic can guarantee, and any sensitivity built
# on this Jacobian (voltage per tap step, voltage per MVar) is numerically
# meaningless. That is exactly what the estimate is for: the classic result
# log reports it as a `Jacobian cond.` line, and diagnose runs grade it with
# a plain-language verdict.

# ## Chapter 2: slack types and short-circuit currents
#
# One grid connection, modeled three ways, plus an IEC 60909-0 fault-current
# sweep from the declared feeder data. The detailed walk-through with full
# reading aids is the
# [slack-types notebook](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb).

function build_grid(mode::Symbol)
  net = Net(name = "tour_eg8_$(mode)", baseMVA = 100.0)
  for b in ("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.010, x_pu = 0.060, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.015, x_pu = 0.080, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.020, x_pu = 0.090, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B4", toBus = "B5", r_pu = 0.012, x_pu = 0.070, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B5", toBus = "B6", r_pu = 0.015, x_pu = 0.075, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B6", toBus = "B7", r_pu = 0.018, x_pu = 0.085, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B7", toBus = "B8", r_pu = 0.010, x_pu = 0.055, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B8", toBus = "B1", r_pu = 0.011, x_pu = 0.065, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B7", r_pu = 0.020, x_pu = 0.100, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B6", r_pu = 0.022, x_pu = 0.110, b_pu = 0.02, status = 1)
  addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 60.0, vm_pu = 1.01, qMin = -60.0, qMax = 60.0)
  addProsumer!(net = net, busName = "B6", type = "GENERATOR", p = 40.0, vm_pu = 1.00, qMin = -40.0, qMax = 40.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 45.0, q = 12.0)
  addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 50.0, q = 15.0)
  addProsumer!(net = net, busName = "B7", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
  addProsumer!(net = net, busName = "B8", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
  addExternalGrid!(net = net, busName = "B1", vm_pu = 1.02, sk_max_MVA = 2000.0, sk_min_MVA = 1500.0, rx_max = 0.1, internal_impedance = (mode === :source))
  ok, msg = validate!(net = net)
  ok || error("Network validation failed: $msg")
  return net
end

# **Ideal slack**: `B1` is pinned at exactly 1.02 pu / 0° and absorbs the
# whole imbalance (the `SLACK` row).

net_slack = build_grid(:slack)
etime, ite = solve!(net_slack)
printACPFlowResults(net_slack, etime, ite, 1e-8)

# **Non-ideal source**: with `internal_impedance = true` the setpoint moves
# to the hidden internal bus (last row, type `SOURCE`); the terminal `B1`
# in the first row droops below 1.02 pu.

net_source = build_grid(:source)
etime, ite = solve!(net_source)
printACPFlowResults(net_source, etime, ite, 1e-8)

# **Distributed slack**: the generators pick up the imbalance according to
# their scheduled output (0.6/0.4); the slack row keeps only the reactive
# balance.

net_dist = build_grid(:slack)
etime, ite = solve!(net_dist; distributed_slack_enabled = true, distributed_slack_p_mode = :pg_weighted)
printACPFlowResults(net_dist, etime, ite, 1e-8)

# The three connection models side by side. The losses differ because the
# flow pattern differs; a negative Q loss means the line charging produces
# more reactive power than the flows consume.

println(rpad("scenario", 20), lpad("Vm(B1) pu", 11), lpad("P loss MW", 11), lpad("Q loss MVAr", 13), "   balanced by")
for (label, net, by) in (
  ("ideal slack", net_slack, "slack bus B1"),
  ("non-ideal source", net_source, "hidden source bus"),
  ("distributed slack", net_dist, "B3 (α=0.6) + B6 (α=0.4)"),
)
  pl, ql = getTotalLosses(net = net)
  println(rpad(label, 20), lpad(string(round(get_bus_vm_pu(net, "B1"); digits = 4)), 11), lpad(string(round(pl; digits = 3)), 11), lpad(string(round(ql; digits = 3)), 13), "   ", by)
end

# **Short circuit**: the feeder data declared in `addExternalGrid!` feeds
# `runShortCircuit!` directly. $I_k''$ is largest at the connection bus and
# decays with electrical distance.

printShortCircuitResult(runShortCircuit!(net_slack; case = :max))
printShortCircuitResult(runShortCircuit!(net_slack; case = :min))

# ## Chapter 3: transformer tap control (OLTC)
#
# A transformer with a ratio tap changer holds the voltage at a remote load
# bus. The outer control loop moves the discrete tap until the target is
# inside the deadband; the power flow itself stays untouched. Details:
# [Control Framework](https://welthulk.github.io/Sparlectra.jl/control_framework/).

function build_oltc()
  net = Net(name = "tour_oltc", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "MV", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", referencePri = "Slack", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = 60.0, q = 20.0)
  addPIModelTrafo!(net = net, fromBus = "Slack", toBus = "MV", r_pu = 0.004, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "MV", toBus = "Load", r_pu = 0.02, x_pu = 0.10, b_pu = 0.01, status = 1)
  ## enable the ratio tap on the transformer branch and give it a name the
  ## controller can address
  t = getNetBranch(net = net, fromBus = "Slack", toBus = "MV")
  t.comp.cName = "T1"
  t.has_ratio_tap = true
  t.tap_min = 0.90
  t.tap_max = 1.10
  t.tap_step = 0.0125
  ok, msg = validate!(net = net)
  ok || error("Network validation failed: $msg")
  return net
end

net_oltc = build_oltc()
run_sparlectra(net = net_oltc)
println("uncontrolled: Vm(Load) = ", round(get_bus_vm_pu(net_oltc, "Load"); digits = 4), " pu")

# Now attach the controller (voltage mode, discrete steps) and rerun.
# `run_sparlectra` executes the outer control loop automatically when
# controllers are present.

addTapController!(
  net_oltc;
  trafo = "T1",
  mode = :voltage,
  target_bus = "Load",
  target_vm_pu = 1.0,
  control_ratio = true,
  control_phase = false,
  is_discrete = true,
  deadband_vm_pu = 5e-3,
)
run_sparlectra(net = net_oltc)
println("controlled:   Vm(Load) = ", round(get_bus_vm_pu(net_oltc, "Load"); digits = 4), " pu")
printTapControllerSummary(stdout, net_oltc)

# Reading aid: the summary shows the chosen tap position and the achieved
# voltage. With a discrete 0.0125 step the controller stops as soon as the
# target is inside the deadband, not at the exact setpoint.

# ## Chapter 4: voltage-dependent reactive power, Q(U)
#
# A machine can follow a Q(U) droop characteristic: absorb reactive power
# when its voltage is high, inject when it is low. Unlike the outer-loop
# controllers above, Q(U) is solved **inside** Newton-Raphson. Details:
# [Voltage Dependent Control](https://welthulk.github.io/Sparlectra.jl/voltage_dependent_control/).

function build_qu(p_load::Float64, q_load::Float64)
  net = Net(name = "tour_qu", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 110.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  qu = QUController(
    make_characteristic(
      [(104.5, 30.0), (107.0, 20.0), (110.0, 0.0), (112.0, -10.0), (115.5, -20.0)];
      voltage_unit = :kV,
      value_unit = :MVAr,
      vn_kV = 110.0,
      sbase_MVA = 100.0,
      interpolation = :polynomial,
    );
    qmin_MVAr = -50.0,
    qmax_MVAr = 50.0,
    sbase_MVA = 100.0,
  )
  addProsumer!(net = net, busName = "B2", type = "SYNCHRONOUSMACHINE", p = 10.0, q = 0.0, qu_controller = qu)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = p_load, q = q_load)
  ok, msg = validate!(net = net)
  ok || error("Network validation failed: $msg")
  return net
end

# Light load: the machine bus sits above 110 kV, so the characteristic asks
# the machine to **absorb** reactive power (negative Q).

net_qu_light = build_qu(5.0, 1.0)
etime, ite = solve!(net_qu_light)
printACPFlowResults(net_qu_light, etime, ite, 1e-8)

# Heavy load: the voltage sags below 110 kV and the same characteristic
# turns the machine into a reactive power **injector** (positive Q).

net_qu_heavy = build_qu(80.0, 25.0)
etime, ite = solve!(net_qu_heavy)
printACPFlowResults(net_qu_heavy, etime, ite, 1e-8)

# Reading aid: compare the `Qg` value and the `Control` column (`Q(U)`) of
# bus `B2` between the two tables; the sign flips with the voltage level,
# exactly along the declared characteristic.

# ## Chapter 5: remote voltage control by a machine
#
# A machine regulates the voltage at a **different** bus via its reactive
# output, the counterpart of a CGMES `RegulatingControl` at a foreign
# terminal. The outer loop drives the machine's Q until the remote target
# is met, and parks honestly `at_limit` when the reactive range is too
# small. Details:
# [Remote Voltage Control](https://welthulk.github.io/Sparlectra.jl/remote_voltage_control/).

function build_rvc(qmin_mvar::Float64, qmax_mvar::Float64)
  net = Net(name = "tour_rvc", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "GenBus", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "GenBus", type = "SYNCHRONOUSMACHINE", p = 30.0, q = 0.0, qMin = qmin_mvar, qMax = qmax_mvar)
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -20.0)
  addPIModelACLine!(net = net, fromBus = "Slack", toBus = "GenBus", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
  addPIModelACLine!(net = net, fromBus = "GenBus", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
  return net
end

pf = PowerFlowConfig(max_iter = 30, tol = 1e-9)

net_rvc = build_rvc(-50.0, 50.0)
runpf!(net_rvc; config = pf, verbose = 0)
println("uncontrolled: Vm(Load) = ", round(get_bus_vm_pu(net_rvc, "Load"); digits = 4), " pu")

addMachineVoltageControl!(net_rvc; bus = "GenBus", target_bus = "Load", target_vm_pu = 1.05, deadband_vm_pu = 5e-4)
result = run_control!(net_rvc; controllers = collect_outer_controllers(net_rvc), pf_config = pf, control_config = ControlConfig(max_outer_iterations = 15), verbose = 0)
println("controlled:   Vm(Load) = ", round(get_bus_vm_pu(net_rvc, "Load"); digits = 4), " pu (target 1.05)")
println("loop: status = ", result.status, ", outer iterations = ", result.outer_iterations, ", pf solves = ", result.powerflow_solves)
printMachineControllerSummary(stdout, net_rvc)

# The honest failure mode: cut the reactive range to ±2 MVAr and the same
# target is out of reach. The controller parks at its limit and says so
# instead of pretending convergence.

net_rvc2 = build_rvc(-2.0, 2.0)
runpf!(net_rvc2; config = pf, verbose = 0)
addMachineVoltageControl!(net_rvc2; bus = "GenBus", target_bus = "Load", target_vm_pu = 1.05, deadband_vm_pu = 5e-4)
result2 = run_control!(net_rvc2; controllers = collect_outer_controllers(net_rvc2), pf_config = pf, control_config = ControlConfig(max_outer_iterations = 15), verbose = 0)
println("limited:      Vm(Load) = ", round(get_bus_vm_pu(net_rvc2, "Load"); digits = 4), " pu (target 1.05)")
printMachineControllerSummary(stdout, net_rvc2)

# ## Chapter 6: a steerable HVDC link (back-to-back pairing)
#
# Two AC areas joined ONLY by an HVDC converter pair: no AC tie, no angle
# coupling, so the areas stay two separate electrical islands with their
# own references. The transfer through the link is a control setpoint, not
# the result of an angle difference. First the Stage-0 view: two fixed
# injections reproduce a snapshot of 80 MW transfer with 4 MW converter
# loss.
#
# Note that EACH island keeps a classical reference of its own (`A1` and
# `C1`); the converters at `A2`/`C2` are plain injections the controller
# will steer. The result header therefore reports `Slack: 2`, and both
# reference buses appear as `SLACK` rows in the bus table:
#
# ```text
#      island A                          island C
#   A1 -------- A2  ===== DC link =====  C2 -------- C1
# (slack)    load 40 MW              load 50 MW    (slack)
#            + converter             + converter
#            (from side)             (to side)
# ```

function build_b2b(name::String)
  net = Net(name = name, baseMVA = 100.0)
  for b in ("A1", "A2", "C1", "C2")
    addBus!(net = net, busName = b, vn_kV = 380.0)
  end
  addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "C1", toBus = "C2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", referencePri = "A1", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "C1", type = "EXTERNALNETWORKINJECTION", referencePri = "C1", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
  addProsumer!(net = net, busName = "C2", type = "ENERGYCONSUMER", p = 50.0, q = 12.0)
  addProsumer!(net = net, busName = "A2", type = "GENERATOR", p = -80.0, q = 0.0)  ## converter, exports
  addProsumer!(net = net, busName = "C2", type = "GENERATOR", p = 76.0, q = 0.0)   ## converter, receives 80 - 4
  return net
end

net6 = build_b2b("tour_b2b")
etime, ite = solve!(net6; islands_enabled = true)
println("two islands solved in ", ite, " iteration(s)")

# Reading aid: both areas balance on their own reference; the link carries
# whatever the snapshot says. To make it steerable, pair the two converter
# injections: `addHvdcPairControl!` enforces the invariant
# $P_\text{to} = P_\text{transfer} - P_\text{loss}$ exactly and lets you
# retarget the transfer.

addHvdcPairControl!(net6; from_bus = "A2", to_bus = "C2", p_transfer_mw = 120.0, loss_mw = 4.0, p_rating_mw = 150.0)
result6 = run_control!(net6; controllers = collect_outer_controllers(net6), pf_config = PowerFlowConfig(method = :rectangular, max_iter = 25, tol = 1e-8), control_config = ControlConfig(max_outer_iterations = 8, trace = false))
calcNetLosses!(net6)
printACPFlowResults(net6, etime, result6.last_pf_iterations, 1e-8)

# Reading aid: the HVDC pair reports inside the `Control` section of the
# classical result output, in the same aligned label/value layout as the
# transformer, machine, and TCSC summaries (`printHvdcPairControllerSummary`
# prints the same block standalone). The link now carries 120 MW instead
# of 80: in the branch table the line A1 -> A2 supplies 40 MW more (area
# A's reference generates the export), while C1 -> C2 turns around and
# carries the received power away from the converter bus.
#
# **Why does the solver still report two islands?** An AC voltage angle is
# only defined *within* one synchronous island, relative to that island's
# own reference. The link transfers power but no angle information (there
# is no branch, no admittance, no angle coupling between the areas), so
# each island keeps its own reference pinned at 0 degrees. The two-island
# report is the model telling you the areas are asynchronous; it would be
# wrong for it to disappear. The peek below makes that visible: both
# reference buses sit at exactly 0.0 deg, and comparing an A-side angle
# with a C-side angle carries no information, because each is measured
# against a different zero.

bus_va_deg(net, bus) = net.nodeVec[net.busDict[bus]]._va_deg  ## peek into the solved state
for (bus, role) in (("A1", "reference of island A"), ("A2", "converter, exports 120 MW"), ("C1", "reference of island C"), ("C2", "converter, receives 116 MW"))
  println(rpad(bus, 4), rpad(role, 27), ": Vm = ", round(get_bus_vm_pu(net6, bus); digits = 4), " pu, Va = ", round(bus_va_deg(net6, bus); digits = 3), " deg")
end

# Reading aid: within island A the angle falls toward A2 (the converter
# bus imports 120 MW plus the local load from the reference), within
# island C it rises toward C2 (the converter bus feeds the island). Each
# gradient is meaningful only against its own 0-degree reference. The same
# controller attaches automatically on import when a MATPOWER case sets
# `matpower_dcline_mode = paired_control` or a CGMES delivery is loaded
# with `hvdc_mode = paired_control`. Theory:
# [HVDC Back-to-Back](https://welthulk.github.io/Sparlectra.jl/hvdc_back_to_back/).

# ### HVDC as the island's source (grid-forming)
#
# Can the converter itself BE the reference (slack) of island C? Not as
# one side of the paired controller: the pairing treats the transfer as a
# *setpoint* on both sides, while a slack's power is the *outcome* of its
# island's balance (load plus losses). One injection cannot be both at
# once, so `addHvdcPairControl!` refuses a reference bus by design.
#
# The converse is a perfectly valid model though, called a grid-forming
# (Vf) converter: the receiving converter IS the island's source. It
# holds voltage and angle at its PCC, and its power output follows from
# whatever the island draws. Think of an offshore platform or an
# asynchronously supplied island grid. Island C below has NO source of its
# own; its reference moves onto the converter bus C2 (compare the sketch
# with the setpoint variant above):
#
# ```text
#      island A                          island C
#   A1 -------- A2  ===== DC link =====  C2 -------- C1
# (slack)    load 40 MW            grid-forming    load 50 MW
#            + sending             converter
#            converter             (= island C reference)
# ```

function build_b2b_source(name::String; sending_mw::Float64 = 0.0)
  net = Net(name = name, baseMVA = 100.0)
  for b in ("A1", "A2", "C1", "C2")
    addBus!(net = net, busName = b, vn_kV = 380.0)
  end
  addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "C2", toBus = "C1", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", referencePri = "A1", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
  ## island C: no classical slack. The receiving converter is the
  ## grid-forming source, holding 1.0 pu / 0 deg at its PCC bus C2.
  addProsumer!(net = net, busName = "C2", type = "EXTERNALNETWORKINJECTION", referencePri = "C2", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "C1", type = "ENERGYCONSUMER", p = 50.0, q = 12.0)
  ## sending-side converter; mirrored below once the island balance is known
  addProsumer!(net = net, busName = "A2", type = "GENERATOR", p = sending_mw, q = 0.0)
  return net
end

net7 = build_b2b_source("tour_b2b_source")
solve!(net7; islands_enabled = true)
p_island_c = get_branch_p_from_to_mw(net7, "C2", "C1")
println("grid-forming converter at C2 delivers ", round(p_island_c; digits = 3), " MW (island C load + line loss)")
for bus in ("C2", "C1")
  println("  ", bus, ": Vm = ", round(get_bus_vm_pu(net7, bus); digits = 4), " pu, Va = ", round(bus_va_deg(net7, bus); digits = 3), " deg")
end

# Reading aid: the transfer is no longer a setpoint. Island C decides how
# much it draws (load plus line loss), the converter delivers exactly
# that, and the island's reference sits at the converter PCC: 1.0 pu and
# 0 degrees at C2, while the load bus C1 hangs below it. The sending side
# must now *mirror* the island draw plus the converter loss:

net8 = build_b2b_source("tour_b2b_mirrored"; sending_mw = -(p_island_c + 4.0))
solve!(net8; islands_enabled = true)
println("sending side A1 -> A2 carries ", round(get_branch_p_from_to_mw(net8, "A1", "A2"); digits = 3), " MW (40 MW local load + ", round(p_island_c + 4.0; digits = 3), " MW export + line loss)")

# And the refusal from above, demonstrated: pairing the grid-forming
# converter is rejected, because its injection is the island balance, not
# a setpoint:

try
  addHvdcPairControl!(net8; from_bus = "A2", to_bus = "C2", p_transfer_mw = 50.0)
catch err
  println(sprint(showerror, err))
end

# The pairing controller automates exactly this mirror: `mode =
# :island_feed` reads the island balance after each solve and keeps the
# sending side matched, with an honest `at_limit` once the island draw
# exceeds `p_rating_mw`. No transfer setpoint is given, the island
# decides:

net9 = build_b2b_source("tour_b2b_grid_forming")
addHvdcPairControl!(net9; from_bus = "A2", to_bus = "C2", mode = :island_feed, loss_mw = 4.0, p_rating_mw = 150.0)
result9 = run_control!(net9; controllers = collect_outer_controllers(net9), pf_config = PowerFlowConfig(method = :rectangular, max_iter = 25, tol = 1e-8), control_config = ControlConfig(max_outer_iterations = 8, trace = false))
calcNetLosses!(net9)
printACPFlowResults(net9, etime, result9.last_pf_iterations, 1e-8)

# Reading aid, and the direct comparison with the setpoint variant above.
# Both versions report `Slack: 2` in the header (one reference per island,
# always), but WHERE the references sit is the whole difference:
#
# - **Setpoint variant:** references at `A1` and `C1`, each island brings
#   its own source. Both converters are plain PQ injections steered by the
#   controller (`-120 / +116 MW` in the bus table), and the transfer is an
#   order the link follows. Look at C1's `Pg`: it is negative, island C's
#   own slack absorbs the surplus the link pumps in.
# - **Grid-forming variant (this table):** references at `A1` and `C2`.
#   The receiving converter itself is island C's `SLACK` row; its `Pg`
#   column is the island balance outcome (load plus line loss, no
#   setpoint anywhere), `C1` carries only load, and the HVDC block in the
#   `Control` section shows the transfer as "mirrored from island draw".
#
# Try `p_rating_mw = 40.0`: the sending side pins at the rating with
# `at_limit = true` and `converged = false`. The island's reference still
# balances in the model (a power flow cannot show the collapse), so the
# honest flag is what marks the undeliverable draw.

# ## Chapter 7: state estimation
#
# Close the loop: solve a reference power flow on the chapter-1 network,
# derive a noisy synthetic measurement set from it, check observability,
# and let the WLS estimator reconstruct the state. The full narrative is
# the
# [state-estimation notebook](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb).

using Random

net_se = build_ring7("tour_se")
ite_pf, status_pf = runpf!(net_se, 40, 1e-10, 0)
status_pf == 0 || error("Power flow did not converge")

std = measurementStdDevs(vm = 1e-3, pinj = 1.0, qinj = 1.0, pflow = 0.7, qflow = 0.7)
setMeasurementsFromPF!(
  net_se;
  includeVm = true,
  includePinj = true,
  includeQinj = true,
  includePflow = true,
  includeQflow = true,
  noise = true,
  stddev = std,
  rng = MersenneTwister(42),
)

gobs = evaluate_global_observability(net_se; flatstart = true, jacEps = 1e-6)
println("observability: ", gobs.quality, " (", gobs.n_measurements, " measurements, ", gobs.n_states, " states)")

se = runse!(net_se; maxIte = 12, tol = 1e-6, flatstart = true, jacEps = 1e-6, updateNet = true)
println("SE converged: ", se.converged, " in ", se.iterations, " iterations")
println("objective J:  ", round(se.objectiveJ; digits = 2), " (dof ", se.dof, ", within 3σ: ", se.jWithin3Sigma, ")")
for (name, idx) in sort(collect(net_se.busDict); by = last)
  v = se.voltages[idx]
  println(rpad(name, 4), "  Vm = ", round(abs(v); digits = 4), " pu   Va = ", round(rad2deg(angle(v)); digits = 3), "°")
end

# Reading aid: with mild noise the estimate reproduces the reference state
# to a few 1e-3 pu, and $J$ lands near the degrees of freedom, the
# textbook health check for a WLS estimator.
#
# ## Where to go next
#
# The focused notebooks with the full narrative:
#
# - [Introduction](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_intro.ipynb)
# - [Slack types and short circuit](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb)
# - [State estimation](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb)
#
# And the documentation for the chapters that have no own notebook yet:
#
# - [Control Framework](https://welthulk.github.io/Sparlectra.jl/control_framework/) and
#   the [Workshop](https://welthulk.github.io/Sparlectra.jl/workshop/) tap-control section
# - [Voltage Dependent Control](https://welthulk.github.io/Sparlectra.jl/voltage_dependent_control/)
# - [Remote Voltage Control](https://welthulk.github.io/Sparlectra.jl/remote_voltage_control/)
# - [Feature Matrix](https://welthulk.github.io/Sparlectra.jl/feature_matrix/)

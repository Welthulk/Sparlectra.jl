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
# This tour is the FIRST HALF of the Sparlectra workshop: install
# [Sparlectra.jl](https://github.com/Welthulk/Sparlectra.jl) once, warm the
# compiler up once, and then climb from the very first bus to the solver's
# control features in one session. The
# [ADVANCED tour](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour_advanced.ipynb)
# continues with the Expert and Beyond tiers (remote voltage control,
# HVDC, state estimation, FACTS, N-1, threads); the focused single-topic
# notebooks go deeper on individual chapters.
#
# After the warm-up (compilation happens there, everything after is fast)
# the chapters climb three tiers:
#
# **Newcomer**
# 1. Your first network, built step by step
#
# **Beginner**
# 2. Working with the model: trust, switching, editing, Q-limits
#
# **Advanced**
# 3. Slack types and short-circuit currents
# 4. Transformer tap control (OLTC)
# 5. Voltage-dependent reactive power, Q(U)
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
#nb ## the old version is still active: restart the runtime, then rerun
#nb ## this cell.

# ## Warm-up and shared helpers
#
# Julia compiles each function on first use. This one cell warms the
# paths the chapters exercise, so nothing stalls mid-tour: the
# Newton-Raphson solver and the IEC 60909 short circuit (chapter 3). The
# `using` clauses and the small helpers of the whole tour live here too,
# collected up top so they cannot be missed.

using Sparlectra
using Random

## solve helper used by all chapters (25 iterations, tolerance 1e-8)
function solve!(net; kwargs...)
  etime = @elapsed begin
    ite, erg = runpf!(net, 25, 1e-8, 0; kwargs...)
  end
  erg == 0 || error("Power flow did not converge (status = $erg)")
  calcNetLosses!(net)
  return etime, ite
end

## tiny warm-up net: a grid connection WITH declared short-circuit data
wnet = Net(name = "warmup", baseMVA = 100.0)
addBus!(net = wnet, busName = "A", vn_kV = 110.0)
addBus!(net = wnet, busName = "B", vn_kV = 110.0)
addExternalGrid!(net = wnet, busName = "A", vm_pu = 1.0, sk_max_MVA = 2000.0, sk_min_MVA = 1500.0, rx_max = 0.1, internal_impedance = false)
addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
addPIModelACLine!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)

t_first = @elapsed runpf!(wnet, 10, 1e-8, 0)
t_second = @elapsed runpf!(wnet, 10, 1e-8, 0)
println("power flow     : first solve ", round(t_first; digits = 2), " s (compiles), second ", round(t_second * 1000; digits = 2), " ms")

t_sc = @elapsed runShortCircuit!(wnet; case = :max)
println("short circuit  : ", round(t_sc; digits = 2), " s, everything warm")

# ## Part I: Newcomer
#
# ## Chapter 1: your first network, built step by step
#
# No input files, no configuration: a complete 110 kV network from scratch,
# validated, solved, and read. The network is seven buses in a ring with
# two cross-connections; `B1` carries the grid connection:
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
#
# Every Sparlectra model starts from a `Net` object; `baseMVA` is the
# system base power for all per-unit conversions. `addBus!` creates the
# electrical nodes (`vn_kV` nominal voltage, `vm_pu`/`va_deg` the solver's
# starting voltage). Note what is NOT declared here: the operational bus
# type (slack / PV / PQ) is derived later from the devices attached to
# each bus.

net1 = Net(name = "tour_first_pf", baseMVA = 100.0)
addBus!(net = net1, busName = "B1", vn_kV = 110.0, vm_pu = 1.02, va_deg = 0.0)
for i in 2:7
  addBus!(net = net1, busName = "B$(i)", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
end

# `addPIModelACLine!` connects buses with a line as a pi-equivalent branch
# in per-unit (`r_pu`, `x_pu`, and `b_pu` for the total charging):

addPIModelACLine!(net = net1, fromBus = "B1", toBus = "B2", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net1, fromBus = "B2", toBus = "B3", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net1, fromBus = "B3", toBus = "B4", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net1, fromBus = "B4", toBus = "B5", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net1, fromBus = "B5", toBus = "B6", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net1, fromBus = "B6", toBus = "B7", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net1, fromBus = "B7", toBus = "B1", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net1, fromBus = "B2", toBus = "B5", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net1, fromBus = "B3", toBus = "B6", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)

# Devices that consume or produce power are `addProsumer!` calls. The
# external network injection at `B1` references its OWN bus as the voltage
# reference; that is what makes `B1` the slack bus. The generator at `B3`
# feeds in 60 MW, the remaining buses carry loads (`p` in MW, `q` in MVar):

addProsumer!(net = net1, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)
addProsumer!(net = net1, busName = "B3", type = "GENERATOR", p = 60.0, q = 10.0)
addProsumer!(net = net1, busName = "B2", type = "LOAD", p = 35.0, q = 10.0)
addProsumer!(net = net1, busName = "B4", type = "LOAD", p = 45.0, q = 15.0)
addProsumer!(net = net1, busName = "B5", type = "LOAD", p = 25.0, q = 8.0)
addProsumer!(net = net1, busName = "B6", type = "LOAD", p = 30.0, q = 10.0)
addProsumer!(net = net1, busName = "B7", type = "LOAD", p = 20.0, q = 6.0)

# `validate!` checks the model for structural problems (unconnected buses,
# missing slack, inconsistent parameters) BEFORE any numerics run; make it
# a habit after every round of model edits. Then `runpf!` runs the
# rectangular Newton-Raphson solver (max iterations, tolerance, verbosity;
# status `0` means converged), `calcNetLosses!` derives branch flows and
# losses from the converged voltages, and `printACPFlowResults` prints the
# classical result tables:

ok1, msg1 = validate!(net = net1)
ok1 || error("Network validation failed: $msg1")
etime, ite = solve!(net1)   ## solve! wraps exactly runpf! + calcNetLosses! (see warm-up)
printACPFlowResults(net1, etime, ite, 1e-8)

# Reading aid: the slack at `B1` covers the difference between 155 MW of
# load, 60 MW of scheduled generation, and the network losses; all bus
# voltages stay near 1.0 pu. Loading a case from a FILE instead is one
# call through the framework workflow:
# `run_sparlectra(casefile = "case14.m", path = ...)` after
# `ensure_casefile("case14.m")`, which downloads the case on demand; the
# result carries the solved net as `result.net`.
#
# The same construction, packed into a function: later chapters (state
# estimation, model editing) reuse this network.

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

# ## Part II: Beginner
#
# ## Chapter 2: working with the model
#
# Solving once is the easy part. This chapter covers what day-to-day work
# actually consists of: judging how much the numbers can be trusted,
# editing and switching the model, exporting it, and letting the solver
# enforce reactive-power limits.
#
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
  _, ite_weak = solve!(net)
  vm3 = round(something(net.nodeVec[3]._vm_pu, NaN); digits = 4)
  println("x = ", lpad(x_weak, 8), " pu:  ", ite_weak, " iterations,  Vm(B3) = ", vm3, ",  kappa = ", round(condestJacobian(net), sigdigits = 3))
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

# ### Editing, switching, exporting
#
# Model work is iterative: change a parameter, switch an element, remove
# one, validate, solve again. The dedicated helpers keep the bookkeeping
# consistent (branch indices, prosumer injections, isolated buses):

net_edit = build_ring7("tour_edit")
## stiffen the B1-B2 line (per-branch parameter update)
brVec = getNetBranchNumberVec(net = net_edit, fromBus = "B1", toBus = "B2")
updateBranchParameters!(net = net_edit, branchNr = brVec[1], branch = BranchModel(0.005, 0.040, 0.0, 0.0, 0.0, 0.0, 100.0))
## add 5 MW / 1 MVAr of load at B4. Loads and generators are PROSUMER
## objects and the AC solver reads its injections from them, so growing a
## load means adding (or editing) a prosumer. The node-sum helpers
## (addBusLoadPower!) only feed the report layer, not the AC solve (#323).
addProsumer!(net = net_edit, busName = "B4", type = "LOAD", p = 5.0, q = 1.0)
## switch the B3-B6 cross-tie out of service (aggregate switch; a single
## terminal would be setBranchTerminalStatus!, see the next section)
tie = getNetBranchNumberVec(net = net_edit, fromBus = "B3", toBus = "B6")
setNetBranchStatus!(net = net_edit, branchNr = tie[1], status = 0)
## remove the B2-B5 cross-tie outright, then re-validate and solve
removeACLine!(net = net_edit, fromBus = "B2", toBus = "B5")
markIsolatedBuses!(net = net_edit, log = false)
ok_e, msg_e = validate!(net = net_edit)
ok_e || error("edit round left the net invalid: $msg_e")
_, ite_e = solve!(net_edit)
println("edited net solves in ", ite_e, " iterations; branches now: ", length(net_edit.branchVec))
## the edited model exports as a MATPOWER case file
case_out = joinpath(mktempdir(), "tour_edit.m")
writeMatpowerCasefile(net_edit, case_out)
println("exported: ", case_out, " (", filesize(case_out), " bytes)")

# Reading aid: removal helpers (`removeACLine!`, `removeTrafo!`,
# `removeProsumer!`, ...) mutate the net and can leave isolated buses
# behind; `markIsolatedBuses!` flags them for the solver and
# `clearIsolatedBuses!` deletes the safe ones. `removeBus!` deliberately
# only CHECKS removability. Re-validate after every edit round.
#
# One trap worth demonstrating (issue #323): there are also NODE-level
# power helpers (`addBusLoadPower!`, `addBusGenPower!`). They edit report
# sums that NO solver reads; the solvers build their injections from the
# prosumer objects. So this "edit" changes nothing:

vm_before = get_bus_vm_pu(net_edit, "B4")
addBusLoadPower!(net = net_edit, busName = "B4", p = 25.0, q = 5.0)  ## report layer only!
solve!(net_edit)
println("after addBusLoadPower!(+25 MW): Vm(B4) = ", round(get_bus_vm_pu(net_edit, "B4"); digits = 4), " pu (before: ", round(vm_before; digits = 4), " pu, unchanged)")
addProsumer!(net = net_edit, busName = "B4", type = "LOAD", p = 25.0, q = 5.0)  ## THIS is a load
solve!(net_edit)
println("after addProsumer!(LOAD, 25 MW): Vm(B4) = ", round(get_bus_vm_pu(net_edit, "B4"); digits = 4), " pu (sags, the solver saw it)")
#
# ### Reactive-power limits (PV to PQ switching)
#
# A voltage-regulating machine holds its bus voltage only while its
# reactive power stays inside `[qMin, qMax]`. When the solver hits a
# limit, it switches the bus from PV to PQ at the violated bound within
# the Newton iteration (the active-set strategy) and reports every event:

net_q = Net(name = "tour_qlimits", baseMVA = 100.0)
for b in ("Q1", "Q2", "Q3")
  addBus!(net = net_q, busName = b, vn_kV = 110.0)
end
addProsumer!(net = net_q, busName = "Q1", type = "EXTERNALNETWORKINJECTION", referencePri = "Q1", vm_pu = 1.0, va_deg = 0.0)
addProsumer!(net = net_q, busName = "Q2", type = "GENERATOR", p = 20.0, vm_pu = 1.05, qMin = -5.0, qMax = 5.0)
addProsumer!(net = net_q, busName = "Q3", type = "ENERGYCONSUMER", p = 45.0, q = 20.0)
addPIModelACLine!(net = net_q, fromBus = "Q1", toBus = "Q2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net_q, fromBus = "Q2", toBus = "Q3", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
validate!(net = net_q)
solve!(net_q)
println("Vm(Q2) = ", round(get_bus_vm_pu(net_q, "Q2"); digits = 4), " pu (setpoint was 1.05)")
printQLimitLog(net_q)
distributeBusResults!(net_q)

# Reading aid: holding 1.05 pu at `Q2` would need more than the 5 MVar the
# machine may deliver, so the solver pins Q at `qMax` and lets the voltage
# float below the setpoint; the log names bus, iteration, and bound.
# `distributeBusResults!` pushes the solved bus totals back onto the
# individual prosumers; with several machines on one bus it redistributes
# water-filling style, so no unit leaves its Q range. The deep dive
# (enforcement modes, guards, oscillation handling) is on the
# [Q-limit strategy page](https://welthulk.github.io/Sparlectra.jl/q_limit_switching_strategy/).
#
# ### Opening one end of a line
#
# A breaker can open a single terminal while the other stays connected.
# The line then carries no through flow, but it does NOT disappear: seen
# from the closed bus it collapses to its exact pi reduction and keeps
# drawing its FULL charging (for realistic lines the two shunt arms act
# almost in parallel, so it is `g + jb`, not half of it), plus the small
# ohmic loss of the charging current. The voltage at the open end rises
# above the feeding bus, the classical Ferranti effect, and is reported as
# a branch result without adding a bus. Open one end with
# `setBranchTerminalStatus!` and compare:

net_open = Net(name = "tour_open_end", baseMVA = 100.0)
addBus!(net = net_open, busName = "A", vn_kV = 380.0)
addBus!(net = net_open, busName = "B", vn_kV = 380.0)
addProsumer!(net = net_open, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
addProsumer!(net = net_open, busName = "B", type = "ENERGYCONSUMER", p = 120.0, q = 30.0)
addPIModelACLine!(net = net_open, fromBus = "A", toBus = "B", r_pu = 0.02, x_pu = 0.16, b_pu = 0.9, g_pu = 0.004, status = 1)
validate!(net = net_open)
solve!(net_open)
println("closed : line carries ", round(net_open.branchVec[1].fBranchFlow.pFlow; digits = 1), " MW to the load")

setBranchTerminalStatus!(net_open.branchVec[1]; to = false)
markIsolatedBuses!(net = net_open, log = false)
solve!(net_open)
br_open = net_open.branchVec[1]
println("open@to: charging draw ", round(br_open.fBranchFlow.qFlow; digits = 1), " MVAr, active loss ", round(br_open.fBranchFlow.pFlow; digits = 3), " MW")
println("         voltage at the OPEN end: ", round(br_open.open_end_vm_pu; digits = 4), " pu, HIGHER than the ", round(get_bus_vm_pu(net_open, "A"); digits = 2), " pu at the feeding bus A")
println("         (Ferranti effect: the charging current flowing through the line reactance lifts the voltage toward the open end)")

# Reading aid: the branch table marks the row `open@to`, the header counts
# it under `Open terminals`, and bus `B` is reported isolated: the open end
# is a branch RESULT (`open_end_vm_pu`), not a solved bus. The full story,
# including the Schur reduction and why it is the full charging, is on the
# branch-model docs page under "One-sided open branches"; the runnable twin
# is `exp_open_terminal_line.jl`.

# ### Links: connections without impedance
#
# A link (`addLink!`) models a busbar coupler or sectionalizer: a closed
# switch between two buses. It is NOT a branch. It has no impedance, it is
# never stamped into the Y-bus, and it never appears in the branch table.
# Instead the solver contracts every cluster of buses joined by closed
# links onto one representative bus before the Y-bus is built, so all
# linked buses share one voltage by construction. (Do not confuse these
# links with the HVDC "Link" rows of chapter 7: a bus link is a switch,
# an HVDC link is a converter pair.)
#
# Because the link has no admittance, the power flow cannot tell how much
# power crosses it: that is reconstructed AFTER the solve, from Kirchhoff's
# current law, with `calcLinkFlowsKCL!`. The interesting case is a ring of
# links, a zero-impedance cycle:
#
# ```text
#      S (slack)
#      |
#      |  real line (r = 0.01, x = 0.08)
#      |
#      R1
#     /  \            R1, R2, R3 joined by three closed links:
#    /    \           an impedance-less ring, electrically ONE node
#   R3 --- R2
# (20 MW) (30 MW)
# ```
#
# In a zero-impedance loop the flow split is physically NOT unique: any
# circulating current can be added without changing a single voltage.
# Sparlectra returns the minimum-norm KCL solution (Moore-Penrose
# pseudoinverse), the unique split with zero artificial circulation, so
# the result is deterministic and reproducible:

net_ring = Net(name = "tour_link_ring", baseMVA = 100.0)
addBus!(net = net_ring, busName = "S", vn_kV = 110.0)
for b in ("R1", "R2", "R3")
  addBus!(net = net_ring, busName = b, vn_kV = 110.0)
end
addProsumer!(net = net_ring, busName = "S", type = "EXTERNALNETWORKINJECTION", referencePri = "S", vm_pu = 1.0, va_deg = 0.0)
addProsumer!(net = net_ring, busName = "R2", type = "ENERGYCONSUMER", p = 30.0, q = 8.0)
addProsumer!(net = net_ring, busName = "R3", type = "ENERGYCONSUMER", p = 20.0, q = 5.0)
addPIModelACLine!(net = net_ring, fromBus = "S", toBus = "R1", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
addLink!(net = net_ring, fromBus = "R1", toBus = "R2", status = 1)
addLink!(net = net_ring, fromBus = "R2", toBus = "R3", status = 1)
addLink!(net = net_ring, fromBus = "R3", toBus = "R1", status = 1)
validate!(net = net_ring)
solve!(net_ring)
calcLinkFlowsKCL!(net_ring)

println("one voltage for the whole ring: Vm(R1/R2/R3) = ", join((round(get_bus_vm_pu(net_ring, b); digits = 5) for b in ("R1", "R2", "R3")), " / "), " pu")
ring_bus_name = Dict(v => k for (k, v) in net_ring.busDict)   ## BusLink stores bus INDICES
for l in net_ring.linkVec
  println("link ", ring_bus_name[l.fromBus], " -> ", ring_bus_name[l.toBus], ": ", lpad(round(l.pFlow_MW; digits = 2), 7), " MW")
end

# Reading aid: 50 MW enter the ring at `R1`. The minimum-norm split sends
# 26.67 MW directly to the 30 MW load at `R2` and 23.33 MW the other way
# round to `R3`; the third coupler carries only the 3.33 MW that `R2` still
# needs. A negative sign just means the flow runs against the link's
# from-to direction. Any other split (say 30 and 20 with an idle third
# coupler) would satisfy KCL too, but only by adding a circulating
# component; the pseudoinverse is exactly the split without one. The
# links page of the docs has the math and the modeling guidelines (for
# example: never link the slack bus itself).

# ## Part III: Advanced
#
# ## Chapter 3: slack types and short-circuit currents
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

# ## Chapter 4: transformer tap control (OLTC)
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

# ## Chapter 5: voltage-dependent reactive power, Q(U)
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


# ## Where to go next
#
# The workshop continues in the
# [ADVANCED tour](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour_advanced.ipynb)
# (Expert and Beyond tiers): remote voltage control by a machine, a
# steerable HVDC link, state estimation, the FACTS devices and their
# limits, N-1 contingency analysis, and parallel sweeps on Julia threads.
#
# The focused notebooks with the full narrative of single topics:
#
# - [Slack types and short circuit](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb)
# - [Distributed slack](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_distributed_slack.ipynb)
# - [Transformer taps](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_transformers.ipynb)
# - [TCSC flow steering](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_series_compensation.ipynb)
# - [State estimation](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb)
#
# And the documentation for going further:
#
# - [Solver Guide](https://welthulk.github.io/Sparlectra.jl/solver/)
# - [Slack and External Grid Sources](https://welthulk.github.io/Sparlectra.jl/slack_vs_source/)
# - [Control Framework](https://welthulk.github.io/Sparlectra.jl/control_framework/)
# - [Voltage Dependent Control](https://welthulk.github.io/Sparlectra.jl/voltage_dependent_control/)
# - [Feature Matrix](https://welthulk.github.io/Sparlectra.jl/feature_matrix/)

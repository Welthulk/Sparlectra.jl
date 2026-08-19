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
# file: docs/lit/workshop_slack_short_circuit.jl                              #src
# purpose: Literate.jl source of the slack-representations and short-circuit  #src
#          workshop notebook and its Documenter page. Follow-up to            #src
#          workshop_intro.jl; regenerate the committed outputs with           #src
#          `julia --project=docs docs/generate_notebooks.jl`.                 #src

# # Slack types and short-circuit currents
#
# [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb)
#
# Every AC power flow needs one bus that balances the network: the *slack*.
# But the real grid connection is not an infinitely stiff busbar; it is the
# superordinate network behind an impedance, and it also determines how much
# short-circuit current arrives at your buses. In this notebook you model one
# and the same grid connection three ways with
# [Sparlectra.jl](https://github.com/Welthulk/Sparlectra.jl) and compare the
# results, then reuse the connection's declared short-circuit power for an
# IEC 60909-0 short-circuit calculation:
#
# 1. **Ideal slack**: the connection bus holds voltage magnitude and angle,
#    no matter what.
# 2. **Non-ideal external-grid source**: the reference voltage sits behind
#    the feeder impedance $Z_Q = U_n^2 / S_k''$, so the connection-bus
#    voltage droops under load.
# 3. **Distributed slack**: the active-power imbalance is picked up by the
#    generators according to participation factors, the primary-control
#    picture.
#
# The theory behind the comparison is on
# [Slack Bus and External Grid Sources](https://welthulk.github.io/Sparlectra.jl/slack_vs_source/);
# the short-circuit method is documented in the
# [Short-Circuit Compendium](https://welthulk.github.io/Sparlectra.jl/short_circuit/).
#
# The test network is an eight-bus ring with two cross-ties (the middle
# verticals B2-B7 and B3-B6); the grid connection under study sits at B1:
#
# ```text
#  (grid connection)
#    B1 ---- B2 ---- B3 ---- B4
#    |       |       |       |
#    B8 ---- B7 ---- B6 ---- B5
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

# ## Load the package

using Sparlectra

# ## The study network
#
# A meshed 110 kV ring `B1..B8` with two chords, two PV generators (60 MW at
# `B3`, 40 MW at `B6`) and 160 MW of load. Scheduled generation deliberately
# undershoots the load, so the grid connection at `B1` has to import a
# visible amount of power: that import is what makes the three slack
# representations distinguishable.
#
# `addExternalGrid!` models the connection as an IEC 60909-0 network feeder:
# it creates the load-flow side (ideal slack by default, non-ideal source
# with `internal_impedance = true`) **and** records the declared
# short-circuit power ($S_{k,\mathrm{max}}'' = 2000$ MVA,
# $S_{k,\mathrm{min}}'' = 1500$ MVA) for the short-circuit engine.

function build_grid(mode::Symbol)
  net = Net(name = "workshop_eg8_$(mode)", baseMVA = 100.0)
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

  addExternalGrid!(
    net = net,
    busName = "B1",
    vm_pu = 1.02,
    sk_max_MVA = 2000.0,
    sk_min_MVA = 1500.0,
    rx_max = 0.1,
    internal_impedance = (mode === :source),
  )

  ok, msg = validate!(net = net)
  ok || error("Network validation failed: $msg")
  return net
end

# A small solve helper so every scenario runs identically (25 iterations
# maximum, tolerance $10^{-8}$):

function solve!(net; kwargs...)
  etime = @elapsed begin
    ite, erg = runpf!(net, 25, 1e-8, 0; kwargs...)
  end
  erg == 0 || error("Power flow did not converge (status = $erg)")
  calcNetLosses!(net)
  return etime, ite
end

# ## Scenario 1: ideal slack
#
# The default: the connection bus `B1` becomes the reference (REF) bus and
# holds exactly 1.02 pu at 0° while absorbing whatever active and reactive
# power the network is missing.

net_slack = build_grid(:slack)
etime, ite = solve!(net_slack)
printACPFlowResults(net_slack, etime, ite, 1e-8)

# Reading aid: the `SLACK` row at `B1` imports the scheduled 60 MW imbalance
# plus all network losses, and `B1` sits at exactly 1.02 pu / 0.000°: the
# ideal, infinitely stiff busbar.

# ## Scenario 2: non-ideal external-grid source
#
# With `internal_impedance = true` the reference voltage moves to a hidden
# internal bus `B1__extgrid_int` behind the feeder impedance
# $z_{pu} = \mathrm{baseMVA} / S_k'' = 100/2000 = 0.05$ (split by the
# declared $R/X = 0.1$). The terminal bus `B1` becomes an ordinary solved
# bus.

net_source = build_grid(:source)
etime, ite = solve!(net_source)
printACPFlowResults(net_source, etime, ite, 1e-8)

# Reading aid: look at the **first row**. The terminal bus `B1` is now an
# ordinary `PQ` bus at about 1.009 pu and -1.6°, below the 1.02 pu
# setpoint. The setpoint itself is held by the hidden internal bus
# `B1__extgrid_int` in the **last row** (type `SOURCE`, exactly 1.020 pu at
# 0°, the actual angle reference): the import current drops the difference
# between those two rows across the feeder impedance. The stiffer the
# declared $S_k''$, the smaller the droop; for $S_k'' \to \infty$ this
# variant degenerates to Scenario 1.

# ## Scenario 3: distributed slack
#
# Scenarios 1 and 2 answer the question "who supplies the missing power?"
# the same way: one dedicated bus absorbs the whole imbalance plus all
# losses. That is an accounting fiction; in a real interconnection,
# primary control raises the output of many machines at once. The
# distributed slack models exactly that. The solver gains one unknown, the
# total correction `lambda_P`, and every participating generator covers
# its share `alpha * lambda_P`. The weights come from the scheduled output
# (`:pg_weighted`: 60 MW gives α = 0.6 for `B3`, 40 MW gives α = 0.4 for
# `B6`). The reference bus keeps the angle reference and the reactive
# balance, but its active import drops to zero.

net_dist = build_grid(:slack)
etime, ite = solve!(net_dist; distributed_slack_enabled = true, distributed_slack_p_mode = :pg_weighted)
printACPFlowResults(net_dist, etime, ite, 1e-8)

# Reading aid: the bus table now carries the participation directly in the
# columns **`dSl alpha`** and **`Pg eff MW`**: `B3` picks up α = 0.6 of the
# 61.5 MW correction (60 to 96.9 MW effective), `B6` the remaining α = 0.4
# (40 to 64.6 MW). The `Pg` column keeps showing the *schedule*, the header
# line above the table sums it up (mode and `lambda_P`), and the slack row
# at `B1` no longer imports active power; only the reactive balance stays
# with the reference bus. The branch flows confirm the pickup: summing them
# around `B3` yields exactly the effective 96.9 MW.
#
# ## The three scenarios side by side
#
# Same network, three grid-connection models. The losses differ because the
# flow pattern differs: in Scenario 3 the extra power from `B3` and `B6`
# travels longer paths through the ring, and in Scenario 2 the feeder
# branch itself dissipates a share. (A negative Q loss means the line
# charging produces more reactive power than the flows consume.)

println(rpad("scenario", 20), lpad("Vm(B1) pu", 11), lpad("P loss MW", 11), lpad("Q loss MVAr", 13), "   balanced by")
for (label, net, by) in (
  ("ideal slack", net_slack, "slack bus B1"),
  ("non-ideal source", net_source, "hidden source bus"),
  ("distributed slack", net_dist, "B3 (α=0.6) + B6 (α=0.4)"),
)
  pl, ql = getTotalLosses(net = net)
  println(rpad(label, 20), lpad(string(round(get_bus_vm_pu(net, "B1"); digits = 4)), 11), lpad(string(round(pl; digits = 3)), 11), lpad(string(round(ql; digits = 3)), 13), "   ", by)
end

# ## Short-circuit currents (IEC 60909-0)
#
# The external grid is more than a voltage boundary condition: its declared
# short-circuit power says how much fault current the superordinate network
# can deliver. `runShortCircuit!` replaces the operating state by the
# equivalent voltage source at the fault location and computes the initial
# symmetrical short-circuit current $I_k''$, power $S_k''$, and peak current
# $i_p$ per fault bus. Only sources with **declared** short-circuit data
# contribute; here that is the feeder at `B1`. The two generators carry no
# short-circuit attributes, so they are simply not short-circuit sources in
# this sweep; near those machines the real fault level would be somewhat
# higher than the feeder-only result below.

sc_max = runShortCircuit!(net_slack; case = :max)
printShortCircuitResult(sc_max)

# The minimum case (protection sensitivity) uses the declared
# $S_{k,\mathrm{min}}'' = 1500$ MVA and the lower IEC 60909-0 voltage
# factor $c_\mathrm{min}$:

sc_min = runShortCircuit!(net_slack; case = :min)
printShortCircuitResult(sc_min)

# Reading aid: $I_k''$ is largest at the connection bus `B1`, whose
# short-circuit level is exactly the declared feeder strength (2000 MVA
# resp. 1500 MVA), and decays with electrical distance as line impedance
# accumulates in the fault loop; `B4`, the electrically farthest bus, sees
# less than a third of the connection-bus current.
#
# ## Where to go next
#
# - New to Sparlectra? The
#   [introduction notebook](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_intro.ipynb)
#   builds a network from scratch step by step, directly in Colab.
# - [Slack Bus and External Grid Sources](https://welthulk.github.io/Sparlectra.jl/slack_vs_source/):
#   the full theory: why the load flow needs a slack, the source model, and
#   how the equation system changes.
# - [Short-Circuit Compendium](https://welthulk.github.io/Sparlectra.jl/short_circuit/):
#   method, c-factors, safety flagging, and CGMES-fed short circuits.
# - `examples/powerflow/exp_external_grid_comparison.jl` in the repository:
#   the same comparison as a script, tabulated bus by bus.

```@meta
EditURL = "../../lit/workshop_slack_short_circuit.jl"
```

# Slack types and short-circuit currents

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb)

Every AC power flow needs one bus that balances the network: the *slack*.
But the real grid connection is not an infinitely stiff busbar; it is the
superordinate network behind an impedance, and it also determines how much
short-circuit current arrives at your buses. In this notebook you model one
and the same grid connection three ways with
[Sparlectra.jl](https://github.com/Welthulk/Sparlectra.jl) and compare the
results, then reuse the connection's declared short-circuit power for an
IEC 60909-0 short-circuit calculation:

1. **Ideal slack**: the connection bus holds voltage magnitude and angle,
   no matter what.
2. **Non-ideal external-grid source**: the reference voltage sits behind
   the feeder impedance $Z_Q = U_n^2 / S_k''$, so the connection-bus
   voltage droops under load.
3. **Distributed slack**: the active-power imbalance is picked up by the
   generators according to participation factors, the primary-control
   picture.

The theory behind the comparison is on
[Slack Bus and External Grid Sources](https://welthulk.github.io/Sparlectra.jl/slack_vs_source/);
the short-circuit method is documented in the
[Short-Circuit Compendium](https://welthulk.github.io/Sparlectra.jl/short_circuit/).

> **Note:** On Google Colab the install cell takes a few minutes on a
> fresh session (package download and precompilation). Colab's Julia
> version may change over time; this notebook targets Julia ≥ 1.12.

## Load the package

````@example workshop_slack_short_circuit
using Sparlectra
````

## The study network

A meshed 110 kV ring `B1..B8` with two chords, two PV generators (60 MW at
`B3`, 40 MW at `B6`) and 160 MW of load. Scheduled generation deliberately
undershoots the load, so the grid connection at `B1` has to import a
visible amount of power: that import is what makes the three slack
representations distinguishable.

`addExternalGrid!` models the connection as an IEC 60909-0 network feeder:
it creates the load-flow side (ideal slack by default, non-ideal source
with `internal_impedance = true`) **and** records the declared
short-circuit power ($S_{k,\mathrm{max}}'' = 2000$ MVA,
$S_{k,\mathrm{min}}'' = 1500$ MVA) for the short-circuit engine.

````@example workshop_slack_short_circuit
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
````

A small solve helper so every scenario runs identically (25 iterations
maximum, tolerance $10^{-8}$):

````@example workshop_slack_short_circuit
function solve!(net; kwargs...)
  etime = @elapsed begin
    ite, erg = runpf!(net, 25, 1e-8, 0; kwargs...)
  end
  erg == 0 || error("Power flow did not converge (status = $erg)")
  calcNetLosses!(net)
  return etime, ite
end
````

## Scenario 1: ideal slack

The default: the connection bus `B1` becomes the reference (REF) bus and
holds exactly 1.02 pu at 0° while absorbing whatever active and reactive
power the network is missing.

````@example workshop_slack_short_circuit
net_slack = build_grid(:slack)
etime, ite = solve!(net_slack)
printACPFlowResults(net_slack, etime, ite, 1e-8)
````

Reading aid: the `SLACK` row at `B1` imports the scheduled 60 MW imbalance
plus all network losses, and `B1` sits at exactly 1.02 pu / 0.000°: the
ideal, infinitely stiff busbar.

## Scenario 2: non-ideal external-grid source

With `internal_impedance = true` the reference voltage moves to a hidden
internal bus `B1__extgrid_int` behind the feeder impedance
$z_{pu} = \mathrm{baseMVA} / S_k'' = 100/2000 = 0.05$ (split by the
declared $R/X = 0.1$). The terminal bus `B1` becomes an ordinary solved
bus.

````@example workshop_slack_short_circuit
net_source = build_grid(:source)
etime, ite = solve!(net_source)
printACPFlowResults(net_source, etime, ite, 1e-8)
````

Reading aid: look at the **first row**. The terminal bus `B1` is now an
ordinary `PQ` bus at about 1.009 pu and -1.6°, below the 1.02 pu
setpoint. The setpoint itself is held by the hidden internal bus
`B1__extgrid_int` in the **last row** (type `SOURCE`, exactly 1.020 pu at
0°, the actual angle reference): the import current drops the difference
between those two rows across the feeder impedance. The stiffer the
declared $S_k''$, the smaller the droop; for $S_k'' \to \infty$ this
variant degenerates to Scenario 1.

## Scenario 3: distributed slack

The single ideal slack is an accounting fiction: in reality, primary
control spreads an imbalance over many machines. With
`distributed_slack_enabled` the active-power mismatch is distributed over
the generators, here weighted by their scheduled output
(`:pg_weighted`: 60 MW ⇒ α = 0.6 for `B3`, 40 MW ⇒ α = 0.4 for `B6`).

````@example workshop_slack_short_circuit
net_dist = build_grid(:slack)
etime, ite = solve!(net_dist; distributed_slack_enabled = true, distributed_slack_p_mode = :pg_weighted)
printACPFlowResults(net_dist, etime, ite, 1e-8)
````

Reading aid: compare with Scenario 1. The slack row at `B1` no longer
imports any active power (empty `Pg`; only the *reactive* balance stays
with the reference bus). The ≈ 61.5 MW of imbalance plus losses moved
into the network: sum the printed branch flows around each generator bus
and you find `B3` injecting ≈ 97 MW and `B6` ≈ 65 MW, their schedule
plus the α-shares 0.6/0.4. Note that the `Pg` column still shows the
60/40 MW *schedule*; the applied correction is part of the solver's
structured status and of the compact summary printed at `verbose ≥ 1`.

## Short-circuit currents (IEC 60909-0)

The external grid is more than a voltage boundary condition: its declared
short-circuit power says how much fault current the superordinate network
can deliver. `runShortCircuit!` replaces the operating state by the
equivalent voltage source at the fault location and computes the initial
symmetrical short-circuit current $I_k''$, power $S_k''$, and peak current
$i_p$ per fault bus. Only sources with **declared** short-circuit data
contribute; here that is the feeder at `B1`. The two generators carry no
short-circuit attributes, so they are simply not short-circuit sources in
this sweep; near those machines the real fault level would be somewhat
higher than the feeder-only result below.

````@example workshop_slack_short_circuit
sc_max = runShortCircuit!(net_slack; case = :max)
printShortCircuitResult(sc_max)
````

The minimum case (protection sensitivity) uses the declared
$S_{k,\mathrm{min}}'' = 1500$ MVA and the lower IEC 60909-0 voltage
factor $c_\mathrm{min}$:

````@example workshop_slack_short_circuit
sc_min = runShortCircuit!(net_slack; case = :min)
printShortCircuitResult(sc_min)
````

Reading aid: $I_k''$ is largest at the connection bus `B1`, whose
short-circuit level is exactly the declared feeder strength (2000 MVA
resp. 1500 MVA), and decays with electrical distance as line impedance
accumulates in the fault loop; `B4`, the electrically farthest bus, sees
less than a third of the connection-bus current.

## Where to go next

- New to Sparlectra? The
  [introduction notebook](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_intro.ipynb)
  builds a network from scratch step by step, directly in Colab.
- [Slack Bus and External Grid Sources](https://welthulk.github.io/Sparlectra.jl/slack_vs_source/):
  the full theory: why the load flow needs a slack, the source model, and
  how the equation system changes.
- [Short-Circuit Compendium](https://welthulk.github.io/Sparlectra.jl/short_circuit/):
  method, c-factors, safety flagging, and CGMES-fed short circuits.
- `examples/powerflow/exp_external_grid_comparison.jl` in the repository:
  the same comparison as a script, tabulated bus by bus.


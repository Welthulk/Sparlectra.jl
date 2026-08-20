```@meta
EditURL = "../../lit/workshop_intro.jl"
```

# Your first power flow with Sparlectra

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_intro.ipynb)

[Sparlectra.jl](https://github.com/Welthulk/Sparlectra.jl) is a Julia
framework for AC power-flow and state-estimation studies. In this guided
first tour you build a small 110 kV network from scratch (seven buses in a
ring with two cross-connections), validate it, solve it with the built-in
Newton-Raphson solver, and read the classical result tables. No input files
are needed; everything is created programmatically.

The network at a glance (B1 carries the grid connection, the diagonals
are the two cross-ties B2-B5 and B3-B6):

```text
 (slack)
   B1 ---- B2 ---- B3 ---- B4
   |         \    /         |
   |          \  /          |
   |           \/           |
   |           /\           |
   B7 ---- B6 ---- B5 ------+
```

> **Note:** When running this notebook on Google Colab, the install cell
> takes a few minutes on a fresh session (package download and
> precompilation). Colab's Julia version may change over time; this
> notebook targets Julia ≥ 1.12.

## Warm-up and shared imports

Julia compiles each function on first use, so the very first solve of a
session carries the compilation cost. This cell loads the package and
warms the power-flow path on a tiny throwaway network; every solve in
the actual walkthrough then runs at full speed.

````@example workshop_intro
using Sparlectra

wnet = Net(name = "warmup", baseMVA = 100.0)
addBus!(net = wnet, busName = "A", vn_kV = 110.0)
addBus!(net = wnet, busName = "B", vn_kV = 110.0)
addProsumer!(net = wnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
addPIModelACLine!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
t_first = @elapsed runpf!(wnet, 10, 1e-8, 0)
t_second = @elapsed runpf!(wnet, 10, 1e-8, 0)
println("first solve (includes compilation): ", round(t_first; digits = 2), " s, second: ", round(t_second * 1000; digits = 2), " ms")
````

## Create an empty network

Every Sparlectra model starts from a [`Net`](https://welthulk.github.io/Sparlectra.jl/reference_network/)
object. `baseMVA` is the system base power used for all per-unit
conversions.

````@example workshop_intro
net = Net(name = "workshop_intro", baseMVA = 100.0)
````

## Add buses

`addBus!` creates the electrical nodes. `vn_kV` is the nominal voltage,
`vm_pu` and `va_deg` provide the starting voltage for the solver. Note that
the operational bus type (slack / PV / PQ) is *not* declared here: it is
derived later from the devices attached to each bus.

````@example workshop_intro
addBus!(net = net, busName = "B1", vn_kV = 110.0, vm_pu = 1.02, va_deg = 0.0)
for i in 2:7
  addBus!(net = net, busName = "B$(i)", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
end
````

## Connect the buses

`addPIModelACLine!` adds a line as a π-equivalent branch with per-unit
parameters (`r_pu` resistance, `x_pu` reactance, `b_pu` total line
charging). The seven buses form a ring, and two extra cross-connections
give the power flow alternative paths.

````@example workshop_intro
addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B4", toBus = "B5", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B5", toBus = "B6", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B6", toBus = "B7", r_pu = 0.012, x_pu = 0.090, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B7", toBus = "B1", r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B2", toBus = "B5", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)
addPIModelACLine!(net = net, fromBus = "B3", toBus = "B6", r_pu = 0.009, x_pu = 0.070, b_pu = 0.0, status = 1)
````

## Attach loads and generation

Devices that consume or produce power are added with `addProsumer!`. The
external network injection at `B1` references its own bus as the voltage
reference: that makes `B1` the slack bus, which balances the network and
fixes voltage magnitude and angle. The generator at `B3` feeds in 60 MW,
and the remaining buses carry loads (`p` in MW, `q` in MVar).

````@example workshop_intro
addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)
addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 60.0, q = 10.0)
addProsumer!(net = net, busName = "B2", type = "LOAD", p = 35.0, q = 10.0)
addProsumer!(net = net, busName = "B4", type = "LOAD", p = 45.0, q = 15.0)
addProsumer!(net = net, busName = "B5", type = "LOAD", p = 25.0, q = 8.0)
addProsumer!(net = net, busName = "B6", type = "LOAD", p = 30.0, q = 10.0)
addProsumer!(net = net, busName = "B7", type = "LOAD", p = 20.0, q = 6.0)
````

## Validate the network

`validate!` checks the model for structural problems (unconnected buses,
missing slack, inconsistent parameters) before any numerics run. Make this
a habit after every round of model edits.

````@example workshop_intro
ok, msg = validate!(net = net)
ok || error("Network validation failed: $msg")
````

## Solve the power flow

`runpf!` runs the rectangular Newton-Raphson solver. The arguments are the
maximum number of iterations, the convergence tolerance, and a verbosity
level; it returns the iteration count and a status (`0` means converged).

````@example workshop_intro
tol    = 1e-8
maxIte = 10

etime = @elapsed begin
  ite, erg = runpf!(net, maxIte, tol, 0)
end
````

## Inspect the results

`calcNetLosses!` computes the branch flows and total network losses from
the converged voltages; `printACPFlowResults` prints the classical result
tables: bus voltages, branch flows, and losses.

````@example workshop_intro
if erg == 0
  calcNetLosses!(net)
  printACPFlowResults(net, etime, ite, tol)
else
  @warn "Power flow did not converge (status = $erg)"
end
````

## Where to go next

You have built, validated, and solved a complete network in a few dozen
lines. From here:

- [Slack types and short-circuit currents](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb):
  the follow-up notebook, directly in Colab: ideal slack vs. external-grid
  source vs. distributed slack, plus IEC 60909-0 fault currents.
- [Workshop](https://welthulk.github.io/Sparlectra.jl/workshop/): file
  import and export, transformers, bus links, tap control, and Q-limits.
- [State estimation from noisy measurements](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb):
  the estimation notebook, directly in Colab: reconstruct the network
  state from redundant, noisy measurements instead of a load
  specification.
- [Feature Matrix](https://welthulk.github.io/Sparlectra.jl/feature_matrix/):
  what Sparlectra covers, at a glance.


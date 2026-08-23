```@meta
EditURL = "../../lit/workshop_series_compensation.jl"
```

# FACTS flow control: from the TCSC to the UPFC

> **Level: Expert**, companion of the advanced tour's FACTS chapter.
> Covers the series-reactance controller (TCSC) and builds up to the full
> **UPFC** (unified power flow controller), including the DC-link-coupled
> model that steers a line's active and reactive flow independently.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_series_compensation.ipynb)

> **Note:** This workshop was created with AI assistance and is reviewed
> and curated by the maintainer; it is not a fully machine-generated text.

In a meshed AC network, power does not follow contracts, it follows
impedance: parallel paths split the transfer in inverse proportion to
their reactances. A TCSC (thyristor controlled series capacitor) exploits
exactly that lever: a variable series reactance in one line steers how
much flow that corridor carries. In this notebook you build a small loop
network with [Sparlectra.jl](https://github.com/Welthulk/Sparlectra.jl),
watch the natural flow split, and then let the
`SeriesReactanceControl` outer-loop controller move the split onto a
target, including the honest failure mode when the target is out of
reach.

> **Note:** On Google Colab the install cell takes a few minutes on a
> fresh session (package download and precompilation). Colab's Julia
> version may change over time; this notebook targets Julia ≥ 1.12.

## Warm-up

Julia compiles each function on first use. This cell loads the package
and warms the two paths this notebook exercises, the power-flow solver
and the outer control loop with a series-reactance controller, on a tiny
throwaway corridor, so the real study runs at full speed.

````@example workshop_series_compensation
using Sparlectra

wnet = Net(name = "warmup", baseMVA = 100.0)
for b in ("A", "M", "B")
  addBus!(net = wnet, busName = b, vn_kV = 110.0)
end
addProsumer!(net = wnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
addPIModelACLine!(net = wnet, fromBus = "A", toBus = "M", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
addPIModelACLine!(net = wnet, fromBus = "M", toBus = "B", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
addPIModelACLine!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
addSeriesReactanceControl!(wnet; fromBus = "A", toBus = "M", p_target_mw = 6.0, x_min_pu = 0.05, x_max_pu = 0.2)
t_ctrl = @elapsed run_control!(wnet; controllers = collect_outer_controllers(wnet), pf_config = PowerFlowConfig(method = :rectangular, max_iter = 15, tol = 1e-8), control_config = ControlConfig(max_outer_iterations = 4, trace = false))
println("warm: power flow plus series-reactance control ", round(t_ctrl; digits = 2), " s (first calls compile)")
````

## Why a series reactance steers flow

Every branch enters the power flow through its admittance matrix (see the
[Branch Model](https://welthulk.github.io/Sparlectra.jl/branchmodel/)
page for the derivation):

```math
Y_{br} = \begin{bmatrix}
  \frac{1}{\tau^2}\left(y_{ser} + \frac{y_{shunt}}{2}\right) &
  -y_{ser}\,\frac{1}{\tau e^{-j\phi}} \\
  -y_{ser}\,\frac{1}{\tau e^{j\phi}} &
  y_{ser} + \frac{y_{shunt}}{2}
\end{bmatrix},
\qquad y_{ser} = \frac{1}{R + jX},
```

with $N = 1$ for lines. The TCSC acts purely through $X$ inside
$y_{ser}$: every accepted controller step changes one branch stamp and
the outer loop re-stamps the Y-bus before the next solve.

For a lossless line the transfer relation

```math
P_{12} = \frac{V_1 V_2}{X}\,\sin(\delta_1 - \delta_2)
```

says that a lower series reactance carries more power at a given angle
difference. In a loop, flow redistributes between the parallel paths
according to their reactance ratio, which is exactly what we are about
to watch.

## A loop network with two corridors

**Example 1: the natural flow split.** 80 MW travel from source `A` to
sink `B` over two parallel corridors.
The upper corridor (`A` to `M2` to `B`) has twice the reactance of the
lower one, so it naturally carries only one third of the transfer.

```text
         +---- M1 ----+      corridor 1: x = 0.10 per line (TCSC here)
         |            |
   A ----+            +---- B (load 80 MW)
 (slack) |            |
         +---- M2 ----+      corridor 2: x = 0.20 per line
```

````@example workshop_series_compensation
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
````

Reading aid (Example 1): the 2:1 reactance ratio produces the 2:1 flow
split, independent of any thermal ratings: the low-reactance corridor
attracts the flow.

## Attach the TCSC and steer the split

**Example 2: steering the split onto a target.**
`addSeriesReactanceControl!` registers the controller on the line from
`A` to `M2`, continuing on the Example 1 network (diagram above).
The target of 35 MW needs a visible reactance move: the
outer loop measures the branch flow after each converged solve, steps
`x_pu` via secant iteration (the first step is a bounded probe, because
the sign of $dP/dX$ depends on the network), and stops inside the
0.5 MW default deadband.

````@example workshop_series_compensation
ctrl = addSeriesReactanceControl!(net; fromBus = "A", toBus = "M2", p_target_mw = 35.0, x_min_pu = 0.02, x_max_pu = 0.30)
run_sparlectra(net = net)
println("steered:  corridor 2 (A->M2) = ", round(get_branch_p_from_to_mw(net, "A", "M2"); digits = 2), " MW (target 35)")
println("          x_pu moved 0.20 -> ", round(ctrl.x_pu; digits = 4), ", status = ", ctrl.status)
````

The controller row from the last control run (Example 2) and the generic
controllable-element view carry the shared vocabulary (actuator, range,
quantity, target) that all outer controllers report:

````@example workshop_series_compensation
cr = latest_control_result(net)
println("outer loop: status = ", cr.status, ", outer iterations = ", cr.outer_iterations, ", pf solves = ", cr.powerflow_solves)
for row in cr.controllers
  println("controller: ", row.controller_name, " achieved ", round(row.achieved_p_mw; digits = 2), " MW of ", row.p_target_mw, " MW, x_pu = ", round(row.x_pu; digits = 4))
end
for e in controllableElements(net)
  println("element:    ", e.element, " | ", e.device, " | ", e.actuator, " in [", e.actuator_min, ", ", e.actuator_max, "] | ", e.quantity, " @ ", e.target)
end
````

## The honest limit

**Example 3: the honest limit.** Ask the same corridor for 70 MW, on a
fresh copy of the Example 1 loop (diagram there), and the range
`[0.02, 0.30]` is not enough: the reactance clamps at the capacitive
end, the branch behaves
as a fixed compensated line, and the controller reports `at_limit`
instead of pretending convergence. The power flow itself stays valid.

````@example workshop_series_compensation
net2 = build_loop()
ctrl2 = addSeriesReactanceControl!(net2; fromBus = "A", toBus = "M2", p_target_mw = 70.0, x_min_pu = 0.02, x_max_pu = 0.30)
run_sparlectra(net = net2)
println("limited: corridor 2 (A->M2) = ", round(get_branch_p_from_to_mw(net2, "A", "M2"); digits = 2), " MW (target 70)")
println("         x_pu = ", round(ctrl2.x_pu; digits = 4), " (clamped), at_limit = ", ctrl2.at_limit, ", converged = ", ctrl2.converged)
````

## From series reactance to the UPFC

The TCSC above steers flow by changing a series REACTANCE: its injected
voltage is in quadrature with the line current, so it exchanges no active
power. A UPFC (unified power flow controller) removes that restriction. It
adds a SHUNT converter at one line end and couples the two converters
through a DC link, so the SERIES converter may inject a voltage of ARBITRARY
phase. The in-phase component now carries active power, balanced through the
DC link by the shunt, and that is the extra degree of freedom: the line can
hold INDEPENDENT active and reactive targets at once.

```text
   bus i (from)      series converter        bus j (to)
     V_i o----+---[ + V_se - ]---[ line ]---o  V_j
              |                                 controlled flow ->
        [ shunt conv ]  P_sh = -P_se (DC balance) + Q_sh
              |
             === DC link ===   (couples the two converters)
```

**Example 4: the quadrature composite (SSSC + STATCOM).** In the quadrature
limit `P_se = 0`, the DC link idles, and the UPFC is exactly an SSSC on the
branch plus a STATCOM at the bus. `addUpfcControl!` (default `model =
:quadrature`) registers that pair as one device, on the loop of Example 1
(diagram above) with a machine added at `M2` for the shunt converter:

````@example workshop_series_compensation
netu = build_loop()
addProsumer!(net = netu, busName = "M2", type = "GENERATOR", p = 0.0, q = 0.0)
upfc = addUpfcControl!(netu; fromBus = "A", toBus = "M2", shunt_bus = "M2",
                       target_bus = "B", target_vm_pu = 0.99, p_target_mw = 35.0,
                       v_inj_max_pu = 0.08, s_max_mva = 40.0)
run_control!(netu)
println("quadrature UPFC: two converter rows for one device:")
for row in controllableElements(netu)
  println("  ", rpad(row.device, 48), row.actuator, ", at_limit = ", row.at_limit)
end
````

Reading aid (Example 4): one call, one composite name, but honestly TWO
result rows, the series (SSSC) and the shunt (STATCOM), each with its own
`at_limit`. This is the whole device in the quadrature limit; it steers ONE
line quantity plus the shunt voltage, but not independent P and Q.

**Example 5: the full model, independent P and Q.** Lifting the quadrature
restriction (`model = :full`) lets the series converter inject an
arbitrary-phase voltage; the line then holds distinct P and Q targets, with
the active part balanced across the DC link. The full model needs the shunt
at the SENDING bus, so a small mesh (the parallel `S->L` path lets the flow
be steered):

```text
   S (slack) --- I ==[UPFC series]== J --- L (load)
                 |                            S ------------- L
           shunt converter at I               (parallel path)
```

````@example workshop_series_compensation
function build_upfc_mesh()
  m = Net(name = "upfc_mesh", baseMVA = 100.0)
  for b in ("S", "I", "J", "L")
    addBus!(net = m, busName = b, vn_kV = 110.0)
  end
  addProsumer!(net = m, busName = "S", type = "EXTERNALNETWORKINJECTION", referencePri = "S", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = m, busName = "I", type = "GENERATOR", p = 0.0, q = 0.0)   # shunt converter
  addProsumer!(net = m, busName = "L", type = "ENERGYCONSUMER", p = 90.0, q = 30.0)
  addPIModelACLine!(net = m, fromBus = "S", toBus = "I", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = m, fromBus = "I", toBus = "J", r_pu = 0.02, x_pu = 0.18, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = m, fromBus = "J", toBus = "L", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = m, fromBus = "S", toBus = "L", r_pu = 0.02, x_pu = 0.16, b_pu = 0.0, status = 1)
  ok, msg = validate!(net = m)
  ok || error("mesh net invalid: $msg")
  return m
end

netf = build_upfc_mesh()
full = addUpfcControl!(netf; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I",
                       p_target_mw = 40.0, q_target_mvar = 10.0, q_shunt_mvar = 0.0,
                       v_inj_max_pu = 0.30, s_max_mva = 120.0,
                       deadband_p_mw = 1e-2, deadband_q_mvar = 1e-2, max_outer_iters = 80)
run_control!(netf; control_config = ControlConfig(max_outer_iterations = 80))
f = full.upfc
println("full UPFC on I->J:")
println("  line P = ", round(f.achieved_p_mw; digits = 2), " MW (target 40) and Q = ", round(f.achieved_q_mvar; digits = 2), " MVAr (target 10), both at once")
println("  series V_se = ", round(abs(f.v_se_pu); digits = 4), " pu, P_se = ", round(f.p_se_mw; digits = 3), " MW")
println("  DC-link balance P_se + P_sh = ", round(f.p_se_mw + f.p_sh_mw; digits = 4), " MW")
````

The classical result tables carry the whole picture: the "Controllers" line
now counts the UPFC, and the "UPFC Control Summary" block reports the line
P/Q targets vs achieved, the series voltage, and the DC-link residual.

````@example workshop_series_compensation
calcNetLosses!(netf)
printACPFlowResults(netf, 0.0, 1, 1e-8)
````

Forcing the series phase back to quadrature collapses P_se to zero, back to
the Example 4 behaviour:

````@example workshop_series_compensation
netfq = build_upfc_mesh()
fq = addUpfcControl!(netfq; model = :full, series_phase = :quadrature, fromBus = "I", toBus = "J",
                     shunt_bus = "I", p_target_mw = 40.0, q_target_mvar = 0.0, q_shunt_mvar = 0.0,
                     v_inj_max_pu = 0.30, s_max_mva = 120.0, deadband_p_mw = 1e-2, max_outer_iters = 80)
run_control!(netfq; control_config = ControlConfig(max_outer_iterations = 80))
println("quadrature-forced: P_se = ", round(fq.upfc.p_se_mw; digits = 4), " MW (zero: no phase-shifter DOF)")
````

Reading aid (Example 5): the number that unlocks independent P and Q is
`P_se`, the active power the series converter pushes through the DC link
(nonzero here, exactly zero when forced to quadrature). First-cut honesty:
the shunt runs on a reactive setpoint (closed-loop shunt voltage is a
follow-up), and the model converges for feasible, moderate targets. The full
limitation list and the phasor picture are on the FACTS Devices page.

## Where to go next

- [FACTS Devices](https://welthulk.github.io/Sparlectra.jl/facts/):
  the device taxonomy (STATCOM, SVC, SSSC, TCSC, both UPFC models), the
  limit-characteristic comparison, and the UPFC phasor picture.
- [Series Compensation (TCSC)](https://welthulk.github.io/Sparlectra.jl/series_compensation/):
  the theory page, compensation degree, device versus model, and the
  resonance guard.
- [Control Framework](https://welthulk.github.io/Sparlectra.jl/control_framework/):
  the outer loop all controllers share, and the uniform element view.
- [Workshop tour](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb):
  all workshop examples in one Colab session.


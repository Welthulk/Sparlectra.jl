# Remote Voltage Control

A machine regulates the voltage magnitude at a bus that is not its own
connection point: its reactive output is the actuator, the voltage at a
foreign target bus the controlled variable. Sparlectra implements this as
the outer-loop controller `MachineVoltageControl` on the
[generic control framework](control_framework.md), like the transformer
tap controllers.

## Why this is not just a PV bus

A PV bus couples actuator (the machine's reactive output $Q_g$) and
controlled variable (the voltage $|V|$) at one node; Newton-Raphson fixes
$|V|$ there and drops $Q_g$ from the unknowns. With remote regulation the
roles separate: at the machine bus $m$ the reactive injection is free
(neither $|V_m|$ nor $Q_m$ is fixed a priori), at the target bus $t$ the
magnitude $|V_t|$ is prescribed without an adjustable injection.

Folding this into the Newton iteration (drop the $Q$ mismatch at $m$, add
$|V_t| - V^{\mathrm{set}}$) changes the Jacobian structure, interacts with
Q-limit switching and couples buses that share no branch. Sparlectra keeps
the inner solver untouched: the machine stays a PQ injection for every
inner solve, and its reactive setpoint moves between solves.

## The scalar control problem

Between two solves the controller sees the scalar map

```math
Q_m \mapsto V_t(Q_m),
```

with everything else (loads, other setpoints, taps) held by the power
flow. Around an operating point the map is close to linear; its slope is
the network's voltage sensitivity

```math
s = \frac{\partial V_t}{\partial Q_m} > 0,
```

positive for any working actuator. Its magnitude follows the electrical
distance between $m$ and $t$ (dominated by the reactance of the path) and
shrinks toward zero when the target bus is far away or held stiff by nearby
sources.

## Secant iteration

The controller solves $V_t(Q_m) = V^{\mathrm{set}}$ by secant iteration,
without Jacobian or probe solves:

1. **Bootstrap.** The first move is a bounded fraction (25 %) of the
   remaining reactive headroom in the expected direction: voltage too low,
   toward `qmax_mvar`; too high, toward `qmin_mvar`.
2. **Secant step.** Every following move uses the two previous operating
   points $(Q^{k-1}, V^{k-1})$ and $(Q^k, V^k)$:

   ```math
   Q^{k+1} = Q^k + \frac{V^{\mathrm{set}} - V^k}{s_k}, \qquad
   s_k = \frac{V^k - V^{k-1}}{Q^k - Q^{k-1}},
   ```

   clamped to $[Q_{\min}, Q_{\max}]$. This settles within the deadband in
   three to five outer iterations.
3. **Physical-sign guard.** A measured slope $s_k \le 0$ (the target barely
   responds, or other controllers moved the state) triggers the bootstrap
   step instead of a step toward the wrong bound.

Convergence: $|V_t - V^{\mathrm{set}}| \le$ `deadband_vm_pu`.
`ControlConfig.max_outer_iterations` caps the loop.

## Reactive limits: honest `at_limit`

The actuator range is the machine's reactive capability
$[Q_{\min}, Q_{\max}]$ (the imported `ReactiveCapabilityCurve` evaluated
at the scheduled P where one exists, else the scalar hull). When the
secant step is clamped at a bound and the target is still outside the
deadband, the controller parks with status `at_limit`, the outer-loop
analogue of a PV bus switching to PQ. The report row carries
`at_limit = true`, `converged = false` and the achieved voltage.

## STATCOM mode: current-based limit

The constant box models a synchronous machine. A STATCOM is a
voltage-source converter bounded by its current, so the deliverable
reactive power scales with the terminal voltage,

```math
Q_{lim}(V) = V \cdot S_{max}
```

with $S_{max}$ the converter rating at 1.0 pu (`s_max_mva`, alternatively
`i_max_ka` converted via $\sqrt{3}\,U_n I_{max}$ at registration). Only the
limit handling changes:

- the bounds $\pm V \cdot S_{max}$ are re-evaluated from the solved
  machine-bus voltage before every outer step (the element row shows the
  currently deliverable range, not the nameplate);
- an at-limit STATCOM whose bound still moves keeps adjusting, so the
  delivered Q tracks the voltage linearly; it parks `at_limit` once the
  bound has settled;
- the machine's own `minQ`/`maxQ` are ignored: the converter current is
  the limit.

The linear collapse ($Q \propto V$) is the STATCOM's advantage over the
SVC's quadratic one ($Q \propto V^2$); see [FACTS Devices](@ref facts_devices).
In range the mode behaves like the constant-Q controller.

## Interaction with the rest of the solver

- **Bus typing.** Machine bus and target bus stay PQ. A PV or slack target
  is rejected (two authorities for one voltage).
- **Q-limit machinery.** A remote-controlled machine is exempt by
  construction: active-set switching only considers `isRegulating`
  prosumers (an RVC machine has `isRegulated = false` and no
  voltage-adjust controller), and PV→PQ switching only acts on PV buses.
  The controller's clamping is the single limit instance.
- **Bookkeeping.** Each applied step updates the machine's `ProSumer.qVal`
  and its bus-level generation sum by the same delta, so per-machine and
  per-bus views stay coherent.
- **Several controllers.** Tap and machine controllers run in the same
  outer loop (`run_control!` evaluates all, applies all, re-solves). One
  machine controller per target bus is enforced; several machines at one
  bus targeting different buses pollute each other's measured
  sensitivities (more outer iterations).
- **Tap and machine controller on one target bus**: `addMachineVoltageControl!`
  warns but does not resolve it (a tap-regulated bus stays PQ); reconfigure
  one of them.
- **Coordinated Q-sharing** among several machines on one target
  (participation factors) is not implemented; the first machine claims
  the target.

## API

| Use | Where |
|---|---|
| Call | `addMachineVoltageControl!` (keywords in the example), `run_control!`, `collect_outer_controllers`, `printMachineControllerSummary` |
| Config key | `cgmes_import.machine_control` (`importCGMES(machine_control = true)`): attaches the controllers for machines whose voltage `RegulatingControl` points at a foreign bus ([CGMES Import](cgmes_import.md)) |
| Result field | report row: `converged`, `at_limit`, achieved voltage; STATCOM element row: the currently deliverable range |
| Example | `examples/others/machine_remote_voltage_control.jl` |

```julia
addMachineVoltageControl!(net;
  bus = "GenBus",            # machine's own bus (PQ machine required)
  target_bus = "Load",       # remote regulated bus (PQ required)
  target_vm_pu = 1.02,
  deadband_vm_pu = 1e-3,     # convergence band
  # qmin_mvar / qmax_mvar default to the machine's minQ/maxQ
)
# STATCOM variant: current-based limit instead of the constant box
addMachineVoltageControl!(net;
  bus = "StatcomBus", target_bus = "Load", target_vm_pu = 1.0,
  s_max_mva = 25.0,          # converter rating at 1.0 pu; Q_lim = V * S_max
)
result = run_control!(net; controllers = collect_outer_controllers(net))
printMachineControllerSummary(stdout, net)
```

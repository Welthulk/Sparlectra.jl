# Remote Voltage Control

Remote voltage control lets a machine regulate the voltage magnitude at a bus
that is **not** its own connection point: the machine's reactive output is the
actuator, the voltage at a foreign *target bus* is the controlled variable.
Sparlectra implements this as an outer-loop controller
(`MachineVoltageControl`) on top of the
[generic control framework](control_framework.md), the same architecture the
transformer tap controllers use.

## Why this is not just a PV bus

A classic PV bus couples two things at one node: the *actuator* (the
machine's reactive output $Q_g$) and the *controlled variable* (the voltage
magnitude $|V|$ of the same bus). The Newton-Raphson formulation exploits
that coupling — at a PV bus, $|V|$ is fixed and $Q_g$ drops out of the
unknowns.

With remote regulation the two roles separate:

- at the **machine bus** $m$, the reactive injection is a free control
  variable — neither $|V_m|$ nor $Q_m$ is fixed a priori;
- at the **target bus** $t$, the voltage magnitude $|V_t|$ is prescribed,
  while the bus itself has no adjustable injection.

Folding this into the inner Newton iteration is possible (the classic
formulation drops the $Q$ mismatch equation at $m$ and adds a
$|V_t| - V^{\mathrm{set}}$ equation instead), but it changes the Jacobian
structure, interacts with the Q-limit switching machinery, and couples buses
that share no branch. Sparlectra deliberately keeps the inner solver
untouched and treats remote regulation as an **outer loop**, exactly like tap
control: the machine stays an ordinary PQ injection for every inner solve,
and its reactive setpoint moves *between* solves.

## The scalar control problem

Between two power-flow solves, the controller sees a scalar map

```math
Q_m \mapsto V_t(Q_m),
```

the voltage magnitude at the target bus as a function of the machine's
reactive output, with everything else (loads, other setpoints, taps) held by
the power flow. Around an operating point this map is close to linear; its
slope is the network's voltage sensitivity

```math
s = \frac{\partial V_t}{\partial Q_m} > 0,
```

which is positive for any physically working actuator: injecting more
reactive power raises the surrounding voltage profile. The magnitude of $s$
depends on the electrical distance between $m$ and $t$ — dominated by the
reactance of the path — and shrinks toward zero when the target bus is
electrically far away or held stiff by nearby sources.

## Secant iteration

The controller solves $V_t(Q_m) = V^{\mathrm{set}}$ with a secant iteration
that never computes a Jacobian and needs no probe solves:

1. **Bootstrap.** With no measured sensitivity yet, the first move is a
   bounded fraction (25 %) of the remaining reactive headroom in the
   physically expected direction — voltage too low → toward `qmax_mvar`,
   too high → toward `qmin_mvar`. A deliberately short first step only costs
   one outer iteration; the secant update extrapolates past it immediately.
2. **Secant step.** Every following move uses the two previous operating
   points $(Q^{k-1}, V^{k-1})$ and $(Q^k, V^k)$:

   ```math
   Q^{k+1} = Q^k + \frac{V^{\mathrm{set}} - V^k}{s_k}, \qquad
   s_k = \frac{V^k - V^{k-1}}{Q^k - Q^{k-1}},
   ```

   clamped to $[Q_{\min}, Q_{\max}]$. Because $V_t(Q_m)$ is nearly linear,
   this typically settles within the deadband in three to five outer
   iterations.
3. **Physical-sign guard.** A measured slope $s_k \le 0$ contradicts the
   physics of a working actuator (it appears when the target barely responds,
   e.g. numerically, or when other controllers moved the state in between).
   The controller then falls back to the bootstrap step instead of stepping
   toward the wrong bound.

Convergence is voltage-based: $|V_t - V^{\mathrm{set}}| \le$
`deadband_vm_pu`. The framework's `ControlConfig.max_outer_iterations` caps
the loop.

## Reactive limits: honest `at_limit`

The actuator range is the machine's reactive capability
$[Q_{\min}, Q_{\max}]$ (from the import: the `ReactiveCapabilityCurve`
evaluated at the scheduled P where one exists, else the scalar hull). When
the secant step is clamped at a bound and the target is still outside the
deadband, the controller parks with status `at_limit` — the machine
physically cannot deliver the target. This is the exact outer-loop analogue
of a PV bus switching to PQ under Q-limit enforcement, and it is reported
honestly instead of iterating further: the report row carries
`at_limit = true`, `converged = false` and the achieved voltage.

## STATCOM mode: current-based limit (issue #297 Draft A)

The constant box above models a synchronous machine. A STATCOM is a
voltage-source converter, and its bound is the converter CURRENT: the
deliverable reactive power scales with the terminal voltage,

```math
Q_{lim}(V) = V \cdot S_{max}
```

with $S_{max}$ the converter rating at 1.0 pu (`s_max_mva`, alternatively
`i_max_ka` converted via $\sqrt{3}\,U_n I_{max}$ at registration). The
controller keeps the full secant machinery and replaces only the limit
handling:

- the symmetric bounds $\pm V \cdot S_{max}$ are re-evaluated from the
  solved machine-bus voltage before every outer step (LIVE bounds; the
  element row shows the currently deliverable range, not the nameplate);
- an at-limit STATCOM whose bound still moves keeps adjusting, so the
  delivered Q TRACKS the sagging or recovering voltage linearly; it parks
  `at_limit` only once the bound has settled;
- the machine's own `minQ`/`maxQ` are deliberately ignored in this mode:
  the converter current is the limit.

The linear collapse ($Q \propto V$) is the STATCOM's defining advantage
over the SVC's quadratic one ($Q \propto V^2$); the comparison table and
the device taxonomy live on the [FACTS Devices](@ref facts_devices) page.
In range, the mode behaves like the constant-Q controller and converges
into the same deadband.

## Interaction with the rest of the solver

- **Bus typing.** The machine bus stays PQ throughout; the target bus stays
  PQ as well (a PV or slack target is already voltage-held by another unit
  and is rejected — there would be two authorities for one voltage).
- **Q-limit machinery.** A remote-controlled machine is exempt from the
  native Q-limit path by construction, twice over: the active-set switching
  only considers `isRegulating` prosumers (an RVC machine has
  `isRegulated = false` and no voltage-adjust controller), and PV→PQ
  switching only acts on PV buses while the machine bus stays PQ throughout.
  The controller's own clamping is therefore the *single* limit instance —
  no double clamping, no fight between outer loop and active set.
- **Bookkeeping.** Each applied step updates both the machine's
  `ProSumer.qVal` and its bus-level generation sum by the same delta, so
  per-machine and per-bus views stay coherent when several injections share
  the bus.
- **Several controllers.** Tap controllers and machine controllers run in the
  same outer loop (`run_control!` evaluates all, then applies all, then
  re-solves). One machine controller per target bus is enforced; several
  machines at one bus targeting *different* buses are possible but their
  measured sensitivities pollute each other — expect more outer iterations.
  One cross-type case is **warned about but not resolved automatically**: a
  tap controller and a machine controller regulating the *same* target bus
  (the PQ check alone cannot catch it, because a tap-regulated bus stays
  PQ). `addMachineVoltageControl!` emits a warning when a transformer
  controller already regulates the target — the two would fight over one
  voltage, so reconfigure one of them; no cached ENTSO-E delivery exercises
  the pattern.
- **Coordinated Q-sharing** among several machines on one target (a power
  plant with n units, participation factors) is not implemented; the first
  machine claims the target, the others keep their scheduled reactive output.

## API

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

Runnable demo: `examples/others/machine_remote_voltage_control.jl` (reachable
target and the `at_limit` outcome). On CGMES deliveries the controllers are
attached by `importCGMES(machine_control = true)` — config key
`cgmes_import.machine_control` — for machines whose voltage
`RegulatingControl` points at a foreign bus; see
[CGMES Import](cgmes_import.md).

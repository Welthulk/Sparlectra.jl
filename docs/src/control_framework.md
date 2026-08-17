# Generic Control Framework

## Purpose

The inner numerical solver remains `runpf!`. The generic orchestration layer is `run_control!`, which runs controller-driven outer iterations around repeated PF solves.

`run_sparlectra` automatically dispatches through `run_control!` when `collect_outer_controllers(net)` returns at least one controller.

Built-in controllers: `PowerTransformerControl` (tap changers — voltage,
active power, combined; theory below) and `MachineVoltageControl`
(remote voltage regulation via machine reactive power — theory in
[Remote Voltage Control](remote_voltage_control.md)).

## Architecture

```text
run_sparlectra (public entry)
        |
        v
collect_outer_controllers(net)
        |
        v
run_control!
        |
        +--> runpf!
        +--> control_evaluate!
        +--> control_propose_update!
        +--> control_apply_update!
        +--> runpf! again
        |
        v
ControlRunResult stored on net.control_result
```

Preferred public entries:

```julia
run_sparlectra(; net = net, ...)
run_sparlectra(; casefile = "case14.m", path = "...", ...)
```

`run_sparlectra` is the preferred public framework entry point for both in-memory
and file-based workflows. For AC power-flow examples, `run_acpflow` is kept as a
thin alias with the same minimal configuration-driven signature. Both names
return `SparlectraRunResult`.

| Layer | Function | Purpose |
|---|---|---|
| Framework | `run_sparlectra` (`run_acpflow` alias) | Import/config/control/solve/output orchestration |
| Solver | `runpf!` | Solve an already built `Net` using `PowerFlowConfig` |
| Control | `run_control!` | Execute outer-loop controllers |
| Import | `createNetFromMatPowerFile` | Convert a MATPOWER file into a `Net` |

## Hook interface

- `AbstractOuterController`
- `AbstractControlState`
- `AbstractControlUpdate`
- `control_initialize!`
- `control_evaluate!`
- `control_propose_update!`
- `control_apply_update!`
- `control_is_converged`
- `control_is_blocked`
- `control_status`
- `control_report_rows`
- `control_trace_rows`
- `control_max_outer_iterations`

`control_max_outer_iterations` provides controller-specific outer-loop limits. The global outer-loop budget is `control.max_outer_iterations` and is combined with controller limits.

For framework runs, `SparlectraRunResult.numerical_converged` remains strictly
about the last numerical PF solve and `solution_available` remains true when
that solve produced a usable solution. `final_converged` is true only when the
numerical PF solve converged, limit validation did not fail, and
`control_status` is either `:none` or `:converged`. In particular, `:blocked`,
`:max_outer_iterations`, and other non-success control terminal states keep the
last usable solution but do not count as final framework convergence.

`control_max_outer_iterations` is currently treated as an internal extension hook (not exported). External custom controllers can still extend it via `Sparlectra.control_max_outer_iterations(::MyController)`.

## Result model

`ControlRunResult` contains:

- `status`
- `converged`
- `outer_iterations`
- `powerflow_solves`
- `last_pf_iterations`
- `last_pf_status`
- `controllers`
- `trace`
- `elements` — the generic controllable-element records at run end (see below)

Terminal statuses:

- `:no_controllers`
- `:disabled`
- `:no_active_controllers`
- `:pf_failed`
- `:converged`
- `:blocked`
- `:max_outer_iterations`

## Legacy status boundary (`erg`)

At the legacy API boundary, `erg` reflects inner numerical PF success/failure only.

- `:pf_failed` maps to failure (`erg = 1`).
- Control-loop outcomes such as `:blocked` or `:max_outer_iterations` are not inner numerical PF failures.

Inspect control-loop outcome via `latest_control_result(net)` or `net.control_result`.

## Latest-result access

Use:

```julia
latest_control_result(net)
net.control_result
```

These expose the latest control run associated with the `Net` instance.

## Declarative controllers in configuration (issue #305)

Controllers can be declared under `control.controllers` instead of being
attached programmatically: one named mapping per controller, `type`
selecting the device function, the remaining keys mirroring its keyword
arguments. The run pipeline applies the declarations to the net before the
outer loop starts; `applyConfiguredControllers!` does the same for
a programmatically built net. Schema, supported types, and validation
behavior are documented in the control section of
[Configuration](configuration.md).

Declarations are idempotent per net: an element that already carries a
controller of the declared type is skipped, so repeated
`run_sparlectra(net = ...)` calls do not stack duplicates. To change an
already-applied controller, rebuild the net or adjust it programmatically.

### Public API decision: no generic `addController!`

A generic `addController!(net, controller)` entry point is deliberately
**not** exposed. Two reasons:

- Attachment is device-specific: transformer controllers live on the
  winding (`side.controls`), machine/shunt/series controllers in
  `net.machineControls`. A generic attach would need per-type knowledge
  anyway and would bypass the reference resolution, exclusivity checks,
  and cross-controller warnings that live in the `add*Control!` functions.
- User-defined `AbstractOuterController` subtypes already have a supported
  path without touching the net:
  `run_control!(net; controllers = [my_controller, ...])` accepts any
  controller list explicitly.

The device-specific `add*Control!` functions therefore remain the only
public way to attach controllers to a `Net`.

## Controllable elements (generic view)

Every registered controller describes itself in one shared vocabulary —
what the element is, which actuator it moves within which limits, and which
quantity it steers onto which target:

```julia
controllableElements(net)   # -> Vector{NamedTuple}
```

Each record carries `name`, `element`, `device`, `actuator`,
`actuator_min`/`actuator_max`, `quantity`, `target`, `target_value`,
`discrete`, `enabled`, and the live `status`/`converged`/`at_limit` flags.
The same records are stored on `ControlRunResult.elements` at the end of a
control run. The view is derived on demand from the registered controllers —
purely reporting, no control behavior attached. Current devices:

| Device | Actuator | Quantity |
|---|---|---|
| OLTC transformer | `:tap_ratio` | `:bus_voltage` |
| Phase-shifting transformer | `:phase_shift_deg` | `:branch_active_power` |
| Combined regulation | `:tap_ratio_and_phase_shift` | both |
| Machine remote voltage control | `:machine_q_mvar` | `:bus_voltage` |
| SVC (variable shunt) | `:shunt_bs_mvar` | `:bus_voltage` |
| TCSC (series compensation) | `:series_x_pu` | `:branch_active_power` |

## SVC: variable-shunt voltage control

`addShuntVoltageControl!(net; bus, target_vm_pu, bs_min_mvar, bs_max_mvar,
…)` adds an SVC-style controller: its own shunt element whose susceptance
(MVAr at 1.0 p.u., capacitive positive) the outer loop moves via secant
iteration to hold the bus voltage. At a limit the susceptance stays clamped
and the reactive output follows the bus voltage squared through the Y-bus —
the constant-B region of a real SVC, reported honestly as `at_limit` (a
Q-limited machine would hold constant Q instead; that difference is the
point of the device model). The bus must be PQ; a second shunt controller
on the same bus is rejected, and a transformer tap controller regulating
the same bus voltage triggers the cross-type warning. `runShortCircuit!`
and the power flow see the SVC only through its shunt stamp — disabled or
absent controllers leave results untouched.

## TCSC: series-reactance flow control

`addSeriesReactanceControl!(net; fromBus, toBus, p_target_mw, x_min_pu,
x_max_pu, ...)` adds a TCSC-like controller on a line branch: the outer
loop moves the branch series reactance `x_pu` within its range via secant
iteration until the branch carries the active-power target (measured in
the registered from to to direction). Every accepted step changes one
branch stamp; the Y-bus is re-stamped before the next solve. At a range
end the branch behaves as a fixed compensated line, reported honestly as
`at_limit`. Transformer branches are rejected (taps own transformer
reactance), and ranges whose series impedance magnitude enters the
resonance guard `eps_z` are refused at registration. Theory and device
mapping: [Series Compensation (TCSC)](series_compensation.md).

## Trace rows (transformer control)

Minimal row schema:

- `outer_iteration`
- `controller_name`
- `controller_type`
- `transformer_id`
- `mode`
- `status`
- `converged`
- `at_limit`
- `achieved_vm_pu`
- `target_vm_pu`
- `achieved_p_mw`
- `target_p_mw`
- `tap_ratio`
- `phase_shift_deg`

## [Transformer regulation theory: OLTC, PST, and combined regulation](@id transformer_regulation_theory)

A regulated transformer inserts an additional voltage into the winding.
Three cases are distinguished by the phase of that added voltage relative to
the winding voltage:

* **In-phase regulation — OLTC** (German *Längsregelung*): the added voltage
  is parallel to the winding voltage. Only the magnitude of the complex
  turns ratio changes — this is the ordinary on-load tap changer (OLTC)
  ratio tap, and its dominant network effect is on voltage magnitudes and
  reactive-power flow.
* **Quadrature regulation — PST** (German *Querregelung*): the added voltage
  is perpendicular (90°). Mostly the angle of the turns ratio changes — this
  is the phase-shifting transformer (PST) or quadrature booster, and its
  dominant effect is on active-power flow through meshed paths.
* **Combined (oblique) regulation** (German *Schrägregelung*; no established
  English term): the added voltage has an intermediate angle, or — the case
  modeled here — the unit carries **both** an in-phase and a quadrature tap.
  The complex turns ratio becomes `n = ρ · e^{jα}` with independently
  switchable magnitude `ρ` (ratio tap) and angle `α` (phase tap).

In Sparlectra's branch model these are the independent branch fields
`tap_ratio` (with `tap_min/max/step`) and `phase_shift_deg`
(`phase_min/max/step_deg`); both enter the complex ratio of the transformer
branch. The idealization is that a ratio step does not move the angle and a
phase step does not move the magnitude (a real asymmetrical combined
regulator couples them; the CGMES-style `PhaseTapChangerModel` on the winding is
staged for that but not yet wired into the branch admittance).

### Control: why V→ratio and P→phase may be split

In transmission grids the sensitivities decouple well: voltage magnitudes
respond mainly to the ratio tap (V–Q coupling) and active-power flow to the
phase angle (P–θ coupling). Real combined-regulation units therefore
usually run **two separate controllers** — a voltage regulator
driving the in-phase stage and an active-power regulator driving the
quadrature stage. Note the practical convention: the "angle controller"
regulates a *power* setpoint; the phase angle is its actuator, not its
control target.

Sparlectra supports both realizations:

1. **One combined controller** — `mode = :voltage_and_branch_active_power`
   with `control_ratio = true` and `control_phase = true` (one status, one
   report row; see `examples/others/tap_control_demo_grid.jl`).
2. **Two independent controllers on the same transformer** (split combined
   regulation) — one `mode = :voltage` controller owning the ratio tap and
   one `mode = :branch_active_power` controller owning the phase tap:

   ```julia
   addTapController!(net; trafo = "T_SCHRAEG", mode = :voltage,
       target_bus = "Load_MV", target_vm_pu = 1.02,
       control_ratio = true, control_phase = false, deadband_vm_pu = 5e-3)

   addTapController!(net; trafo = "T_SCHRAEG", mode = :branch_active_power,
       target_branch = ("HV", "MV"), p_target_mw = 120.0,
       control_ratio = false, control_phase = true, deadband_p_mw = 2.0)
   ```

   Each channel keeps its own target, deadband, convergence status, and
   report/trace rows; the outer loop alternates both until each is inside
   its deadband. This mirrors the separately parameterized device
   controllers in the field. Demo:
   `examples/others/tap_control_schraeg_two_controllers.jl`.

### Per-actuator exclusivity

`addPowerTransformerControl!` enforces that each actuator (ratio tap, phase
tap) of a transformer is driven by **at most one** active controller. Two
controllers on one transformer are accepted exactly when their actuator
sets are disjoint; a second claim on an already-driven actuator raises an
error.

### Discrete-step sizing rule

With discrete taps the deadband must cover at least the effect of half a
tap step on the controlled quantity (e.g. one 0.5° phase step moving
≈2.5 MW requires `deadband_p_mw ≥ ~1.5`). A tighter deadband makes the
controller hunt around the target — alternating steps without ever
converging — until `max_outer_iterations` stops the loop.

# Generic Control Framework

`runpf!` is the inner numerical solver; `run_control!` runs
controller-driven outer iterations around repeated power-flow solves, and
`run_sparlectra` dispatches through it whenever
`collect_outer_controllers(net)` returns a controller. Built-in
controllers: `PowerTransformerControl` (tap changers: voltage, active
power, combined; theory below), `MachineVoltageControl`
([Remote Voltage Control](remote_voltage_control.md)),
`ShuntVoltageControl` (SVC, and MSC/MSR switched banks via `step_mvar`),
`SeriesReactanceControl` (TCSC, and SSSC via `v_inj_max_pu`) and
`HvdcPairControl` (back-to-back HVDC pairs).

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

Public entries (`run_acpflow` is an alias; both return
`SparlectraRunResult`):

```julia
run_sparlectra(; net = net, ...)
run_sparlectra(; casefile = "case14.m", path = "...", ...)
```

**Use**

| Layer | Function | Purpose |
|---|---|---|
| Framework | `run_sparlectra` (`run_acpflow` alias) | Import, configuration, control, solve, output |
| Solver | `runpf!` | Solve a built `Net` with a `PowerFlowConfig` |
| Control | `run_control!` | Execute the outer-loop controllers; `run_control!(net; controllers = [...])` runs user-defined ones |
| Configuration | `control.controllers`, `control.max_outer_iterations` | Declarative controllers, global outer-loop budget ([Configuration](configuration.md)) |
| Result | `latest_control_result(net)`, `net.control_result`, `controllableElements(net)` | The `ControlRunResult` of the last run and the generic element records |
| `ControlRunResult` fields | `status`, `converged`, `outer_iterations`, `powerflow_solves`, `last_pf_iterations`, `total_pf_iterations` (sum over all passes; the result header reports the last pass and this total), `last_pf_status`, `controllers`, `trace`, `elements` | Terminal `status`: `:no_controllers`, `:disabled`, `:no_active_controllers`, `:pf_failed`, `:converged`, `:blocked` (a pass in which no controller moved; one closing solve at the final controller state follows and counts in both totals), `:max_outer_iterations` |
| Trace rows (transformer control) | `outer_iteration`, `controller_name`, `controller_type`, `transformer_id`, `mode`, `status`, `converged`, `at_limit`, `achieved_vm_pu`, `target_vm_pu`, `achieved_p_mw`, `target_p_mw`, `tap_ratio`, `phase_shift_deg` | One row per controller and outer iteration in which it acted; a controller held back by its [priority](@ref controller_priority_ranks) writes none |
| Element record (`controllableElements`) | `name`, `element`, `device`, `actuator`, `actuator_min`/`actuator_max`, `quantity`, `target`, `target_value`, `discrete`, `enabled`, live `status`/`converged`/`at_limit` | `ControlRunResult.elements` stores the same records at run end |
| Framework result | `SparlectraRunResult.numerical_converged` (the last PF solve only), `solution_available` (that solve is usable), `final_converged` (converged PF solve, passed limit validation, `control_status` `:none` or `:converged`) | `:blocked`, `:max_outer_iterations` and other non-success states keep the last usable solution without counting as framework convergence |
| Legacy boundary | `erg` reflects inner PF success only: `:pf_failed` maps to `erg = 1` | `:blocked` or `:max_outer_iterations` are not PF failures |

A user-defined controller subtypes `AbstractOuterController` (with
`AbstractControlState` and `AbstractControlUpdate`) and implements the
hook methods of `run_control!` from the
[Controllers Reference](reference_controller.md); its own outer-loop
limit combines with the global budget `control.max_outer_iterations`.

## Declarative controllers in configuration

Under `control.controllers` each named mapping declares one controller:
`type` selects the device function, the remaining keys mirror its keyword
arguments. The run pipeline applies the declarations before the outer
loop; `applyConfiguredControllers!` does the same for a programmatically
built net. Schema, types and validation: [Configuration](configuration.md).
An element that already carries a controller of the declared type is
skipped, so repeated runs do not stack duplicates (a transformer counts per
actuator: entries on its ratio and on its phase tap are two controllers); to change it, rebuild
the net or adjust it programmatically. There is no generic
`addController!`: attachment is device-specific (transformer controllers
live on the winding in `side.controls`, the others in
`net.machineControls`), and the `add*Control!` functions carry the
reference resolution and exclusivity checks.

## [Controller priority](@id controller_priority_ranks)

Every controller takes an integer `priority` (keyword of the
`add*Control!` functions, key `priority` in `control.controllers`, default
0). In each pass the controllers act by rank, smallest number first. A
rank acts only if no controller with a smaller number moved its actuator in
that pass: each of them is inside its deadband, at its limit, or parked.
Controllers with the same number act together, so the default (all 0) is
the loop without ranks. Unlike `referencePriority` of a machine, 0 is a
rank here (the first), not "no preference".

A locked controller (the dead time of a sequential [QSTS](qsts.md) run)
stays in its rank and is evaluated, but moves nothing. While it stands
outside its band and not at its limit it holds the larger numbers; in its
band or at its limit they act.

A controller is at its limit when its actuator stands at the end of its
range:

| Controller | At its limit |
|---|---|
| Tap changer | ratio at the end of its control range: `tap_min`/`tap_max` of the controller where set, of the transformer otherwise (phase tap: `phase_min_deg`/`phase_max_deg`) |
| Switched bank (`step_mvar`) | on its outermost block, or parked (no whole block brings the voltage closer to the target) |
| SVC | susceptance at `bs_min_mvar`/`bs_max_mvar` |
| Machine, STATCOM | reactive output at its Q limit (STATCOM: `V * S_max`, which moves with the terminal voltage) |
| TCSC, SSSC, UPFC, HVDC pair | at the end of their own range or rating |

The common case is a switched bank behind the tap changer of its bus:

```yaml
control:
  controllers:
    bank_LV:
      type: shunt_voltage
      bus: LV
      target_vm_pu: 1.0
      bs_min_mvar: 0.0
      bs_max_mvar: 20.0
      step_mvar: 5.0
      deadband_vm_pu: 0.02
      priority: 1
```

With the tap controller at rank 0 the bank switches only when the tap
stands in its band or at a tap limit. At equal rank both move in the same
pass, and the bank switches as soon as the bus leaves the bank's band,
while the tap would have brought the voltage back alone. The cascade needs
targets that agree within both deadbands (a bank at 1.01 pu behind a tap
at 1.0 +/- 0.01 pu agrees): with targets further apart the lower rank
pulls the voltage out of the tap's band whenever the tap is in band, and
the two walk against each other at any ranks. The cross-type warning at
registration names that case.

A tap controller can get its own control range inside the physical one
(`tap_min`/`tap_max` on the `power_transformer` entry or as keywords of
`addPowerTransformerControl!`). The range sets where the bank behind it
takes over. A ratio outside the range goes onto its nearer end in the
first pass, whatever the voltage error; until then the controller does
not count as converged. The move logs a warning with the transformer, its
start ratio and step, the limit and the ratio it goes to, once per
controller and run; a QSTS run over its base case and all steps is one
run. The entry below sets the range of the controller the case file puts
on `T1`; on a transformer without a controller the entry adds one and
needs `target_bus` and `target_vm_pu` as well:

```yaml
control:
  controllers:
    tap_T1:
      type: power_transformer
      trafo: T1
      mode: voltage
      tap_min: 0.95
      tap_max: 1.05
```

A move along the limit does not hold the next rank: a controller that was
at its limit at the end of the previous pass and is there again after its
move (a STATCOM following `V * S_max`, a tap settling onto `tap_min` from
a ratio a rounding error above it) lets the next rank act in the same
pass. A last partial tap step onto a limit off the step grid is a move.

A waiting rank costs passes: after the last move of a higher rank one
solve runs before the lower rank acts, so a long chain of ranks may need a
larger `control.max_outer_iterations`. A controller held back proposes
nothing and writes no trace row for that pass. `control_priority(ctrl)`
returns the rank; a user-defined controller without a `priority` field
ranks at 0.

## Controllable elements (generic view)

Every registered controller describes itself in one vocabulary (record
fields in the Use table above):

```julia
controllableElements(net)   # -> Vector{NamedTuple}
```

| Device | Actuator | Quantity |
|---|---|---|
| OLTC transformer | `:tap_ratio` | `:bus_voltage` |
| Phase-shifting transformer | `:phase_shift_deg` | `:branch_active_power` |
| Combined regulation | `:tap_ratio_and_phase_shift` | both |
| Machine remote voltage control | `:machine_q_mvar` | `:bus_voltage` |
| SVC (variable shunt) | `:shunt_bs_mvar` | `:bus_voltage` |
| MSC/MSR (switched shunt bank) | `:shunt_bs_mvar` | `:bus_voltage` |
| TCSC (series compensation) | `:series_x_pu` | `:branch_active_power` |
| Back-to-back HVDC pair | `:hvdc_p_transfer_mw` | `:hvdc_transfer` |

## SVC: variable-shunt voltage control

`addShuntVoltageControl!(net; bus, target_vm_pu, bs_min_mvar, bs_max_mvar,
...)` adds an SVC controller: its own shunt element whose susceptance
(MVAr at 1.0 p.u., capacitive positive) the outer loop moves by secant
iteration to hold the bus voltage; at a limit the susceptance stays
clamped (`at_limit`). `step_mvar` switches to the discrete MSC/MSR mode;
the STATCOM alternative is `addMachineVoltageControl!` with `s_max_mva`.
Limit physics of all three: [FACTS Devices](@ref facts_devices).

The bus must be PQ; a second shunt controller on the same bus is rejected,
a tap controller regulating the same bus triggers the cross-type warning,
except for a shunt controller at another priority with the same target
within both deadbands (the cascade in
[Controller priority](@ref controller_priority_ranks)).
`runShortCircuit!` and the power flow see the SVC only through its shunt
stamp.

## TCSC: series-reactance flow control

`addSeriesReactanceControl!(net; fromBus, toBus, p_target_mw, x_min_pu,
x_max_pu, ...)` adds a TCSC controller on a line branch: the outer loop
moves the series reactance `x_pu` within its range by secant iteration
until the branch carries the active-power target (registered from-to
direction), re-stamping the Y-bus before each solve; at a range end the
branch is a fixed compensated line (`at_limit`). Transformer branches are
rejected (taps own transformer reactance), and ranges whose series
impedance magnitude enters the resonance guard `eps_z` are refused. Theory
and the SSSC mode (`v_inj_max_pu`):
[Series Compensation (TCSC)](series_compensation.md); the UPFC
(`addUpfcControl!`): [FACTS Devices](@ref facts_devices).

## Master/slave groups for parallel transformers

Two transformers in parallel between the same busbars regulate as a group;
unequal tap positions would drive a circulating reactive power around the
loop, and two secant loops on one target voltage oscillate against each
other:

```julia
addPowerTransformerControl!(net; trafo = "1", followers = ["2"],
                            mode = :voltage, target_bus = "LV",
                            target_vm_pu = 1.0, deadband_vm_pu = 5e-3)
```

The master runs the normal discrete voltage loop; every accepted master
move is mirrored onto the followers step-synchronously (whole steps of
each follower's own `tap_step`, clamped to its range). Declaratively the
group is the `followers` list on a `power_transformer` entry. A second
independent voltage controller on an already-regulated target bus
triggers a warning; a follower cannot carry its own ratio controller and
cannot follow two groups. Synchronized steps multiply the voltage effect
per master step by the group size, so the deadband must cover at least
half the aggregated step effect.

CGMES: `RatioTapChanger`s sharing one `RegulatingControl`
(`TapChanger.TapChangerControl`) import as a group, the first enabled one
as master (message `follows the group of ...`); `controlEnabled = false`
stays fixed. The exporter writes one shared `TapChangerControl` per active
voltage tap controller (mode voltage, regulated-bus terminal in EQ;
`enabled`, `targetValue`/`targetDeadband` in kV in SSH), referenced by the
master's and every follower's `RatioTapChanger` with
`controlEnabled = true`, so a group survives a CGMES roundtrip.

## [Transformer regulation theory: OLTC, PST, and combined regulation](@id transformer_regulation_theory)

A regulated transformer inserts an additional voltage into the winding;
its phase relative to the winding voltage distinguishes three cases:

* **In-phase regulation, OLTC** (German *Längsregelung*): the added voltage
  is parallel to the winding voltage, so only the magnitude of the complex
  turns ratio changes: the ordinary ratio tap, acting mainly on voltage
  magnitudes and reactive-power flow.
* **Quadrature regulation, PST** (German *Querregelung*): the added voltage
  is perpendicular (90°), so mostly the angle changes: the phase-shifting
  transformer or quadrature booster, acting mainly on active-power flow.
* **Combined (oblique) regulation** (German *Schrägregelung*): the added
  voltage has an intermediate angle, or, as modeled here, the unit carries
  both an in-phase and a quadrature tap: `n = ρ · e^{jα}` with
  independently switchable magnitude `ρ` (ratio tap) and angle `α` (phase
  tap).

In the branch model these are the independent fields `tap_ratio` (with
`tap_min/max/step`) and `phase_shift_deg` (`phase_min/max/step_deg`); both
enter the complex ratio of the transformer branch. On a transformer
without a typed tap model a ratio step does not move the angle and a phase
step does not move the magnitude; a typed `PhaseTapChangerModel` on the
winding couples them like a real asymmetrical regulator, and the
[branch model](branchmodel.md) derives the split of its move.

### Control: why V→ratio and P→phase may be split

In transmission grids voltage magnitudes respond mainly to the ratio tap
(V-Q coupling) and active-power flow to the phase angle (P-θ coupling), so
combined-regulation units usually run a voltage regulator on the in-phase
stage and an active-power regulator on the quadrature stage (the "angle
controller" regulates a *power* setpoint). Sparlectra supports both:

1. **One combined controller**: `mode = :voltage_and_branch_active_power`
   with `control_ratio = true` and `control_phase = true` (one status, one
   report row; `examples/others/tap_control_demo_grid.jl`).
2. **Two independent controllers on one transformer** (split combined
   regulation): a `mode = :voltage` controller owning the ratio tap and a
   `mode = :branch_active_power` controller owning the phase tap:

   ```julia
   addTapController!(net; trafo = "T_SCHRAEG", mode = :voltage,
       target_bus = "Load_MV", target_vm_pu = 1.02,
       control_ratio = true, control_phase = false, deadband_vm_pu = 5e-3)

   addTapController!(net; trafo = "T_SCHRAEG", mode = :branch_active_power,
       target_branch = ("HV", "MV"), p_target_mw = 120.0,
       control_ratio = false, control_phase = true, deadband_p_mw = 2.0)
   ```

   Each channel keeps its own target, deadband, status, and report/trace
   rows. Demo: `examples/others/tap_control_schraeg_two_controllers.jl`.

### Per-actuator exclusivity and step sizing

Each actuator (ratio tap, phase tap) is driven by at most one active
controller: two controllers on one transformer need disjoint actuator
sets, and a second claim on a driven actuator raises an error. With
discrete taps the deadband must cover at least the effect of half a tap
step on the controlled quantity (one 0.5° phase step moving about 2.5 MW
requires `deadband_p_mw ≥ ~1.5`); a tighter deadband makes the controller
hunt around the target until `max_outer_iterations` stops the loop.

On a transformer whose winding carries a typed tap model (`taps`,
`phase_taps`, see [Branch model](branchmodel.md)), the controller moves
the model's step and the resolver rewrites ratio, shift and the tap
dependent reactance; the step sizes above then are the model's own, which
for a symmetrical or asymmetrical phase shifter are not uniform in
degrees. The controller report rows show the model step as the position.

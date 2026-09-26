# Generic Control Framework

`runpf!` is the inner numerical solver; `run_control!` runs
controller-driven outer iterations around repeated power-flow solves, and
`run_sparlectra` dispatches through it whenever
`collect_outer_controllers(net)` returns a controller. Built-in
controllers: `PowerTransformerControl` (tap changers: voltage, active
power, combined; theory below), `MachineVoltageControl`
([Remote Voltage Control](remote_voltage_control.md)),
`ShuntVoltageControl` (SVC, and MSC/MSR switched banks via `step_mvar`),
`SeriesReactanceControl` (TCSC, and SSSC via `v_inj_max_pu`), and
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
| Framework | `run_sparlectra` (`run_acpflow` alias) | Import/config/control/solve/output orchestration |
| Solver | `runpf!` | Solve an already built `Net` using `PowerFlowConfig` |
| Control | `run_control!` | Execute outer-loop controllers; `run_control!(net; controllers = [...])` runs user-defined ones |
| Import | `createNetFromMatPowerFile` | Convert a MATPOWER file into a `Net` |
| Configuration | `control.controllers`, `control.max_outer_iterations` | Declarative controllers, global outer-loop budget ([Configuration](configuration.md)) |
| Result | `latest_control_result(net)`, `net.control_result`, `controllableElements(net)` | The `ControlRunResult` of the last run and the generic element records |
| `ControlRunResult` fields | `status`, `converged`, `outer_iterations`, `powerflow_solves`, `last_pf_iterations`, `total_pf_iterations` (sum over all passes; the result header reports the last pass and this total), `last_pf_status`, `controllers`, `trace`, `elements` | Terminal `status`: `:no_controllers`, `:disabled`, `:no_active_controllers`, `:pf_failed`, `:converged`, `:blocked`, `:max_outer_iterations` |
| Trace rows (transformer control) | `outer_iteration`, `controller_name`, `controller_type`, `transformer_id`, `mode`, `status`, `converged`, `at_limit`, `achieved_vm_pu`, `target_vm_pu`, `achieved_p_mw`, `target_p_mw`, `tap_ratio`, `phase_shift_deg` | One row per controller and outer iteration |
| Element record (`controllableElements`) | `name`, `element`, `device`, `actuator`, `actuator_min`/`actuator_max`, `quantity`, `target`, `target_value`, `discrete`, `enabled`, live `status`/`converged`/`at_limit` | `ControlRunResult.elements` stores the same records at run end |
| Framework result | `SparlectraRunResult.numerical_converged` (the last PF solve only), `solution_available` (that solve is usable), `final_converged` (converged PF solve, passed limit validation, `control_status` `:none` or `:converged`) | `:blocked`, `:max_outer_iterations` and other non-success states keep the last usable solution without counting as framework convergence |
| Legacy boundary | `erg` reflects inner PF success only: `:pf_failed` maps to `erg = 1` | `:blocked` or `:max_outer_iterations` are not PF failures |

A user-defined controller subtypes `AbstractOuterController` (with
`AbstractControlState` and `AbstractControlUpdate`) and implements the
hook methods of `run_control!` from the
[Controllers Reference](reference_controller.md); its own outer-loop
limit combines with the global budget `control.max_outer_iterations`.

## Declarative controllers in configuration

Controllers can be declared under `control.controllers`: one named mapping
per controller, `type` selecting the device function, the remaining keys
mirroring its keyword arguments. The run pipeline applies the declarations
before the outer loop; `applyConfiguredControllers!` does the same for a
programmatically built net. Schema, types and validation:
[Configuration](configuration.md). An element that already carries a
controller of the declared type is skipped, so repeated runs do not stack
duplicates; to change it, rebuild the net or adjust it programmatically.

!!! details "Why there is no generic `addController!`"
    Attachment is device-specific (transformer controllers live on the
    winding in `side.controls`, the others in `net.machineControls`), and
    the `add*Control!` functions carry the reference resolution,
    exclusivity checks and cross-controller warnings that
    `addController!(net, controller)` would have to duplicate.

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
...)` adds an SVC-style controller: its own shunt element whose
susceptance (MVAr at 1.0 p.u., capacitive positive) the outer loop moves
by secant iteration to hold the bus voltage; at a limit the susceptance
stays clamped and the reactive output follows $V^2$ (`at_limit`).
`step_mvar` switches to the discrete MSC/MSR mode; the STATCOM
alternative is `addMachineVoltageControl!` with `s_max_mva`
([FACTS Devices](@ref facts_devices)).

**Notes**

- The bus must be PQ; a second shunt controller on the same bus is
  rejected, a tap controller regulating the same bus triggers the
  cross-type warning.
- `runShortCircuit!` and the power flow see the SVC only through its
  shunt stamp.

## TCSC: series-reactance flow control

`addSeriesReactanceControl!(net; fromBus, toBus, p_target_mw, x_min_pu,
x_max_pu, ...)` adds a TCSC-like controller on a line branch: the outer
loop moves the series reactance `x_pu` within its range by secant
iteration until the branch carries the active-power target (registered
from-to direction), re-stamping the Y-bus before each solve; at a range
end the branch is a fixed compensated line (`at_limit`). Theory:
[Series Compensation (TCSC)](series_compensation.md). The same controller
carries the SSSC mode (`v_inj_max_pu` instead of the fixed window);
`addUpfcControl!` adds a UPFC as `model = :quadrature` (SSSC plus STATCOM
as one named device) or `model = :full` (one controller steering line P
and Q independently): [FACTS Devices](@ref facts_devices).

**Notes**

- Transformer branches are rejected (taps own transformer reactance).
- Ranges whose series impedance magnitude enters the resonance guard
  `eps_z` are refused.

## Master/slave groups for parallel transformers

Two transformers in parallel between the same busbars regulate as a
group:

```julia
addPowerTransformerControl!(net; trafo = "1", followers = ["2"],
                            mode = :voltage, target_bus = "LV",
                            target_vm_pu = 1.0, deadband_vm_pu = 5e-3)
```

The master runs the normal discrete voltage loop; every accepted master
move is mirrored onto the followers step-synchronously (whole steps of
each follower's own `tap_step`, clamped to its range). Declaratively the
group is the `followers` list on a `power_transformer` entry.

**Notes**

- A second independent voltage controller on an already-regulated target
  bus triggers a warning; a follower cannot carry its own ratio
  controller and cannot follow two groups.
- Synchronized steps multiply the voltage effect per master step by the
  group size, so the deadband must cover at least half the aggregated
  step effect.
- CGMES: `RatioTapChanger`s sharing one `RegulatingControl`
  (`TapChanger.TapChangerControl`) import as a group, the first enabled
  one as master, the others as followers (message
  `follows the group of ...`); `controlEnabled = false` stays fixed. The
  exporter writes one shared `TapChangerControl` per active voltage tap
  controller (mode voltage, regulated-bus terminal in EQ; `enabled`,
  `targetValue`/`targetDeadband` in kV in SSH), referenced by the
  master's and every follower's `RatioTapChanger` with
  `controlEnabled = true`, so a group survives a CGMES roundtrip.

!!! details "Why a group instead of two controllers"
    Unequal tap positions on parallel units drive a circulating reactive
    power around the loop (one step apart in opposite directions splits
    the flows of the example into -35 and +53 MVAr), and two secant
    loops on one target voltage oscillate against each other.

## [Transformer regulation theory: OLTC, PST, and combined regulation](@id transformer_regulation_theory)

A regulated transformer inserts an additional voltage into the winding; its
phase relative to the winding voltage distinguishes three cases:

* **In-phase regulation, OLTC** (German *Längsregelung*): the added voltage
  is parallel to the winding voltage, so only the magnitude of the complex
  turns ratio changes: the ordinary ratio tap, acting mainly on voltage
  magnitudes and reactive-power flow.
* **Quadrature regulation, PST** (German *Querregelung*): the added voltage
  is perpendicular (90°), so mostly the angle changes: the phase-shifting
  transformer or quadrature booster, acting mainly on active-power flow.
* **Combined (oblique) regulation** (German *Schrägregelung*): the added
  voltage has an intermediate angle, or, as modeled here, the unit carries
  **both** an in-phase and a quadrature tap: `n = ρ · e^{jα}` with
  independently switchable magnitude `ρ` (ratio tap) and angle `α` (phase
  tap).

In the branch model these are the independent fields `tap_ratio` (with
`tap_min/max/step`) and `phase_shift_deg` (`phase_min/max/step_deg`); both
enter the complex ratio of the transformer branch. A ratio step does not
move the angle and a phase step does not move the magnitude (a real
asymmetrical combined regulator couples them; the CGMES-style
`PhaseTapChangerModel` on the winding is staged for that but not wired
into the branch admittance). The [branch model](branchmodel.md) page
derives the split of a typed phase tap changer's move.

### Control: why V→ratio and P→phase may be split

In transmission grids the sensitivities decouple well: voltage magnitudes
respond mainly to the ratio tap (V-Q coupling), active-power flow to the
phase angle (P-θ coupling). Real combined-regulation units therefore
usually run a voltage regulator on the in-phase stage and an active-power
regulator on the quadrature stage (the "angle controller" regulates a
*power* setpoint). Sparlectra supports both:

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
dependent reactance; the step sizes above then are the model's own,
which for a symmetrical or asymmetrical phase shifter are not uniform in
degrees. The controller report rows show the model step as the position.

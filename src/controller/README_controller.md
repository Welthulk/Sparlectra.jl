# Outer-loop controllers (developer notes)

## Purpose

The generic outer-control framework and its controllers: transformer taps,
remote machine voltage control, variable shunts (SVC), series reactance
(TCSC), UPFC, and back-to-back HVDC pairing. Controllers run OUTSIDE the
Newton iteration: evaluate after a solved power flow, propose a discrete or
continuous update, apply it, and re-solve until settled.

## The controller contract

Every controller implements the `AbstractOuterController` interface from
`control_framework.jl`:

```text
control_evaluate!       read the solved state, decide whether to act
control_propose_update! compute the next actuator value
control_apply_update!   write it into the network
```

`run_control!` owns the outer loop: baseline solve, controller rounds,
convergence bookkeeping, and the `ControlRunResult`.

## Include order

`control_framework.jl` is included early in `src/Sparlectra.jl` (the config
structs reference its types); the concrete controllers follow after the
network types, and `controller_config.jl` last because it instantiates
controllers from configuration.

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `control_framework.jl` | Interface, outer loop, result types | `run_control!`, `ControlConfig`, `ControlRunResult` |
| `controller_config.jl` | Instantiate controllers declared in configuration | `applyConfiguredControllers!` |
| `tap_control.jl` | Transformer ratio/phase tap controllers | `addTapController!` |
| `machine_control.jl` | Remote voltage control for machines | `addMachineVoltageControl!` |
| `shunt_control.jl` | SVC-style variable-shunt voltage controller | `addShuntVoltageControl!` |
| `series_reactance_control.jl` | TCSC-like series-reactance controller | `addSeriesReactanceControl!` |
| `upfc_control.jl` | UPFC as stationary quadrature composite | `addUpfcControl!` |
| `upfc_full_control.jl` | Full UPFC: series injection plus shunt support | `addUpfcFullControl!` |
| `hvdc_pair_control.jl` | Back-to-back HVDC pairing controller | `addHvdcLink!`, `addHvdcPairControl!` |

## Conventions

- Controllers act on solved states only; nothing here may reach into the
  Newton loop or the active-set logic.
- Each controller family has `clear*Controllers!` and a private `_*_controllers`
  registry accessor; follow that pattern for new controllers.
- Cross-type interactions (for example a tap controller and a shunt controller
  regulating the same bus) are not guarded at construction time; see
  `docs/src/remote_voltage_control.md`.

## Reference

Rendered API documentation: `docs/src/reference_controller.md`. User-facing
documentation: `docs/src/control_framework.md` and `docs/src/facts.md`.

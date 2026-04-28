# Transformer Control (Complex Tap, Outer Loop)

## Concept

Sparlectra models transformer regulation in the existing branch PI model with a complex tap
`t = τ * exp(jφ)` and **without auxiliary nodes**. The ratio `τ` is used for voltage control,
and the phase shift `φ` for branch active power control (PST behavior).

A coupled Schrägregler is represented by enabling both controls on one transformer.

## Equations

With branch series admittance `y`, total shunt `y_sh`, and from-side complex tap `t`:

- `Y_ff = (y + 0.5*y_sh)/|t|^2`
- `Y_ft = -y/conj(t)`
- `Y_tf = -y/t`
- `Y_tt = y + 0.5*y_sh`

This is the exact Sparlectra sign/conjugate convention.

## Numerical Method

Current implementation uses an outer loop around `runpf!`:

1. Solve PF with current tap values.
2. Evaluate controller errors.
3. Update `tap_ratio` and/or `phase_shift_deg` (continuous or discrete).
4. Re-run PF.
5. Stop on deadband convergence, tap limits, or `max_outer_iters`.

The augmented Newton formulation with tap variables is intentionally deferred.

## API

```julia
addPowerTransformerControl!(net;
  trafo = "1",
  mode = :voltage,
  target_bus = "B5",
  target_vm_pu = 1.01,
  control_ratio = true,
  control_phase = false,
  is_discrete = true)
```

```julia
addPowerTransformerControl!(net;
  trafo = "1",
  mode = :branch_active_power,
  target_branch = ("B1", "B2"),
  p_target_mw = 250.0,
  control_ratio = false,
  control_phase = true,
  is_discrete = true)
```

```julia
addPowerTransformerControl!(net;
  trafo = "1",
  mode = :voltage_and_branch_active_power,
  target_bus = "B5",
  target_vm_pu = 1.01,
  target_branch = ("B1", "B2"),
  p_target_mw = 250.0,
  control_ratio = true,
  control_phase = true,
  is_discrete = true)
```

## Discrete Tap Behavior

- ratio: `tap_ratio_new = clamp(tap_ratio ± tap_step, tap_min, tap_max)`
- phase: `phase_shift_deg_new = clamp(phase_shift_deg ± phase_step_deg, phase_min_deg, phase_max_deg)`

## Limits / Scope (v1)

- No internal/auxiliary transformer control nodes.
- No coupled Jacobian with tap states yet.
- No multi-transformer coordinated remote control yet.

You can also pass controller definitions directly while creating a transformer:

```julia
ctrl = PowerTransformerControl(
  trafo = "",
  mode = :voltage,
  target_bus = "B5",
  target_vm_pu = 1.01,
  control_ratio = true,
  control_phase = false,
)

addPIModelTrafo!(
  net = net,
  fromBus = "B1",
  toBus = "B2",
  r_pu = 0.01,
  x_pu = 0.08,
  b_pu = 0.0,
  ratio = 1.0,
  shift_deg = 0.0,
  status = 1,
  controls = [ctrl],
)
```

See `src/examples/example_transformer_tap.jl` for a runnable example.

Internally, passed controller channels are also attached to the selected transformer winding
(`PowerTransformerWinding.controls`), so controller-side assignment is explicitly represented
on the transformer model itself.

# Transformer Control (Complex Tap, Outer Loop)

## Concept

Sparlectra models transformer regulation within the existing branch PI model using a complex tap
`t = τ * exp(jφ)` and **without auxiliary nodes**.

* `τ` → voltage control
* `φ` → active power flow control (PST behavior)

A combined transformer (“Schrägregler”) is represented by enabling both controls on one device.

---

## Equations

With series admittance `y`, total shunt `y_sh`, and from-side tap `t`:

```math
\begin{aligned}
Y_{ff} &= \frac{y + 0.5\,y_{sh}}{\lvert t \rvert^2} \\
Y_{ft} &= -\frac{y}{\overline{t}} \\
Y_{tf} &= -\frac{y}{t} \\
Y_{tt} &= y + 0.5\,y_{sh}
\end{aligned}
```

This matches the Sparlectra sign and conjugation convention.

---

## Numerical Method

Transformer control is implemented as an **outer loop around the power flow**:

1. Solve PF with current taps
2. Evaluate control error
3. Update tap(s) (continuous or discrete)
4. Re-run PF
5. Stop on convergence, limits, or iteration cap

No augmentation of the Newton system is performed.

---

## Controller Resolution

Transformer tap/phase control is implemented as an `AbstractOuterController`.
Controllers are collected by `collect_outer_controllers(net)` from
`PowerTransformerWinding.controls`.

They are:

* deduplicated by identity
* used consistently for execution and reporting
* evaluated to determine if outer-loop control is required

---

## Design Rationale

Tap control is intentionally **kept outside the Newton solver**:

* no extension of Jacobian/state vector
* simpler solver backends
* centralized control logic (deadbands, limits, discrete steps)
* deterministic post-processing between iterations

Taps are treated as **supervisory control updates**, not algebraic unknowns.

---

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

## Generic control framework integration

Preferred high-level solve path for already constructed networks:

```julia
run_acpflow(net = net; show_results = false)
```

When at least one controller is present (`collect_outer_controllers(net)`), the
high-level runner calls `run_control!` internally for outer-loop orchestration.
`runpf!` remains the inner numerical solver.

Advanced direct use:

```julia
result = run_control!(
  net;
  pf_config = powerflow_config(),
  control_config = control_config(),
)

result.status
result.trace
latest_control_result(net)
```

`run_control!` remains available for advanced direct use.
`run_net_acpflow` is retained only as a compatibility/secondary wrapper.
There is no public compatibility execution path for `run_tap_controllers_outer!`.

### Result inspection

```julia
run_acpflow(net = net; show_results = false)

result = latest_control_result(net)

println(result.status)
println(result.outer_iterations)
println(result.powerflow_solves)
println(result.controllers)
println(result.trace)
```

`result.status` is the outer control-loop terminal state.
It is separate from numerical PF success/failure.
`result.trace` is machine-readable and does not require parsing console output.

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

---

## Discrete Tap Behavior

* ratio:
  `tap_ratio_new = clamp(tap_ratio ± tap_step, tap_min, tap_max)`

* phase:
  `phase_shift_deg_new = clamp(phase_shift_deg ± phase_step_deg, phase_min_deg, phase_max_deg)`

---

## Limits / Scope

* No auxiliary transformer nodes
* No coupling of tap variables into NR
* No coordinated multi-transformer control

---

## Remote Voltage Control (Current vs. Complete)

### Current capability

Sparlectra already supports **basic remote voltage control**:

* measurement: `target_bus`
* actuator: one transformer tap
* objective: `target_vm_pu ± deadband`

This corresponds to a **single-controller remote regulation**.

---

### What is missing for full remote voltage control

A complete implementation (as used in real grid control systems or CGMES-based models) typically requires:

* **multiple transformers controlling the same remote bus**
* **coordination between controllers**, e.g.:

  * participation factors / weighting
  * priority rules
* **limit handling with redistribution**
  (when one transformer hits tap limits)
* **anti-hunting / stabilization mechanisms**
* **deterministic group convergence logic**

---

### Key gap in Sparlectra

Currently:

* controllers operate **independently**
* no grouping or coordination exists
* no shared control objective across multiple transformers

---

### Implication

Sparlectra implements:

→ **single-actuator remote control**

but not yet:

→ **coordinated multi-actuator remote voltage control**

---

## Example: Inline Controller Definition

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

---

## Runnable example

`examples/tap_control_demo_grid.jl`

This runnable example demonstrates:

* building a small grid
* adding tap/phase controllers through the public API
* running `run_acpflow(net = net; ...)`
* inspecting `latest_control_result(net)`
* printing controller rows and trace rows

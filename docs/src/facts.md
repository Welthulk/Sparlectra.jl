# [FACTS Devices](@id facts_devices)

FACTS devices (Flexible AC Transmission Systems) control bus voltage,
branch flow, and the sharing of power between parallel paths. In a
steady-state power flow each device is a controllable network parameter
(a reactive injection, a shunt susceptance, a series reactance, or a
converter power) driven toward a target by the
[outer control loop](control_framework.md) around `runpf!`. Hands-on:
chapter 4 of the [advanced workshop tour](generated/workshop_tour_advanced.md)
and `examples/others/exp_facts_limit_modes.jl`.

## Device family and where each one lives

| Device | Kind | Actuator | Sparlectra controller | Theory |
|---|---|---|---|---|
| SVC (thyristor-controlled reactor/capacitor) | shunt | susceptance $B$, continuous | `addShuntVoltageControl!` | this page, [Control Framework](control_framework.md) |
| MSC/MSR (mechanically switched capacitor/reactor bank) | shunt | susceptance $B$, whole blocks | `addShuntVoltageControl!` with `step_mvar` | this page |
| STATCOM (VSC shunt converter) | shunt | reactive current | `addMachineVoltageControl!` with `s_max_mva` | this page, [Remote Voltage Control](remote_voltage_control.md) |
| TCSC (thyristor-controlled series capacitor) | series | reactance $x$ in a fixed window | `addSeriesReactanceControl!` with `x_min_pu`/`x_max_pu` | [Series Compensation](series_compensation.md) |
| SSSC (VSC series converter) | series | reactance deviation, voltage-bounded | `addSeriesReactanceControl!` with `v_inj_max_pu` | this page, [Series Compensation](series_compensation.md) |
| PST / Schrägregler (phase-shifting transformer) | phase | tap angle | `addPowerTransformerControl!` (`mode = :branch_active_power`) | [Control Framework](control_framework.md) |
| HVDC back-to-back (paired VSC/LCC converters) | converter | paired P injections, Q or voltage per terminal | `addHvdcPairControl!` | [HVDC Back-to-Back](hvdc_back_to_back.md) |
| UPFC (combined shunt + series converter) | combined | series voltage (quadrature composite, or arbitrary-phase full model) + shunt | `addUpfcControl!` (`model = :quadrature` or `:full`) | this page |

All controllers report through `ControlRunResult`,
`controllableElements` (element, device, actuator with live range, target,
`status`/`converged`/`at_limit`), and the per-controller summary printers.

## The limit characteristic is the device

In range, every shunt compensator holds a voltage target by injecting
reactive power. The devices differ at their limit, which is where the
depressed-voltage conditions they were installed for put them.

**Synchronous machine (constant-Q box).** A fixed reactive capability
$Q \in [Q_{min}, Q_{max}]$, independent of the terminal voltage. At the limit
the machine is a fixed injection: the `MachineVoltageControl` default
(`limit_mode = :constant_q`), the outer-loop analogue of PV to PQ switching.

**SVC (constant-B limit).** A continuous susceptance $B$ within
$[B_{min}, B_{max}]$. At the limit the susceptance clamps and the delivered
reactive power follows the voltage quadratically through the Y-bus stamp:

```math
Q_{SVC} = V^2 \, B_{lim}
```

At $V = 0.9$ pu a fully switched-in SVC delivers 81 percent of its rating.
`ShuntVoltageControl` needs no extra modeling for this: the clamped
susceptance stays stamped in the Y-bus.

**STATCOM (constant-current limit).** A voltage-source converter bounded by
its converter current, so the deliverable reactive power scales linearly
with the terminal voltage:

```math
Q_{STATCOM} = V \, I_{max} \quad\text{, in Sparlectra: } Q_{lim} = V \cdot S_{max}
```

with $S_{max}$ the rating at 1.0 pu; at $V = 0.9$ pu the STATCOM still
delivers 90 percent. `addMachineVoltageControl!` with `s_max_mva` (or
`i_max_ka`, converted via $\sqrt{3}\,U_n I_{max}$) selects this mode. The
bound $\pm V \cdot S_{max}$ is re-evaluated from the solved terminal
voltage before every outer step, so an at-limit STATCOM tracks the
voltage.

Delivered reactive power at the capacitive limit, relative to the rating:

| Terminal voltage | Machine box | STATCOM ($\propto V$) | SVC ($\propto V^2$) |
|---:|---:|---:|---:|
| 1.00 pu | 100 % | 100 % | 100 % |
| 0.95 pu | 100 % | 95 % | 90 % |
| 0.90 pu | 100 % | 90 % | 81 % |
| 0.80 pu | 100 % | 80 % | 64 % |

The machine column is idealized (a real machine derates through its
capability curve).

## Discrete banks: MSC/MSR

A switched bank shares the shunt physics above and lives in the same
controller. With `step_mvar` the susceptance moves in whole blocks (e.g.
four times 10 MVAr):

- the secant proposal is truncated toward the target to whole steps, so
  the bank never overshoots and cannot hunt between two adjacent blocks;
- when no whole block improves the voltage further, the controller parks
  on the reached step (`status = :parked`, the last step before crossing,
  so under-compensation rather than overvoltage) and releases itself when
  another controller moves the operating point enough for a block to help;
- at the outermost block the constant-B limit applies: delivered Q follows
  $V^2$ with the last block connected.

Coordination of a bank with an OLTC on the same bus is not implemented.

## Series side: fixed window versus voltage-bounded window

**TCSC (fixed reactance window).** The branch reactance moves within a
hardware window $[x_{min}, x_{max}]$, independent of loading; at a window
end the branch is a fixed compensated line. The impedance-magnitude guard
`eps_z` excludes the resonance region (see
[Series Compensation](series_compensation.md)).

**SSSC (injected-voltage window).** The converter injects a voltage in
quadrature with the line current: in steady state a reactance deviation from
the natural line reactance $x_{base}$, bounded by the injectable voltage:

```math
|V_{inj}| = |I| \cdot |x - x_{base}| \le V_{inj,max}
\quad\Longleftrightarrow\quad
|x - x_{base}| \le \frac{V_{inj,max}}{|I|}
```

The window shrinks with loading: at high transfer the SSSC saturates,
while a TCSC keeps its full window. In `addSeriesReactanceControl!` the
mode is selected with `v_inj_max_pu`; the window
$x_{base} \pm V_{inj,max}/|I|$ is re-evaluated from the solved branch
current before every outer step, with a floor on $|I|$ (a currentless
branch is unconstrained) and the `eps_z` guard applied as a clamp.

## Live bounds in the outer loop

STATCOM and SSSC share the live bound: the actuator range follows the
solved operating point and is refreshed at the start of every outer
iteration, before the secant step is clamped against it. A parked
controller whose bound still moves is released and keeps adjusting; it
reports `at_limit` only once its bound has settled. `controllableElements`
and the report rows show the bounds of the last evaluated operating
point, so `actuator_min`/`actuator_max` are the currently deliverable
range, not the nameplate. In range, both modes behave like their
fixed-limit counterparts.

## Usage

Programmatic (full keyword sets in the docstrings):

```julia
# SVC: continuous susceptance, quadratic limit collapse
addShuntVoltageControl!(net; bus = "B", target_vm_pu = 1.0,
                        bs_min_mvar = -60.0, bs_max_mvar = 60.0)

# MSC/MSR: the same controller as a switched bank, whole 10-MVAr blocks
addShuntVoltageControl!(net; bus = "B", target_vm_pu = 1.0,
                        bs_min_mvar = -40.0, bs_max_mvar = 40.0,
                        step_mvar = 10.0)

# STATCOM: current-based limit, linear in V
addMachineVoltageControl!(net; bus = "B", target_bus = "C",
                          target_vm_pu = 1.0, s_max_mva = 25.0)

# TCSC: fixed reactance window
addSeriesReactanceControl!(net; fromBus = "A", toBus = "B",
                           p_target_mw = 35.0, x_min_pu = 0.02, x_max_pu = 0.30)

# SSSC: injected-voltage window, shrinks with loading
addSeriesReactanceControl!(net; fromBus = "A", toBus = "B",
                           p_target_mw = 35.0, v_inj_max_pu = 0.05)
```

Declarative, under `control.controllers` (see
[Control Framework](control_framework.md)):

```yaml
control:
  enabled: true
  controllers:
    statcom_c:
      type: machine_voltage
      bus: B
      target_bus: C
      target_vm_pu: 1.0
      s_max_mva: 25.0
    sssc_ab:
      type: series_reactance
      from_bus: A
      to_bus: B
      p_target_mw: 35.0
      v_inj_max_pu: 0.05
```

## UPFC: the stationary quadrature composite

The unified power flow controller combines a STATCOM (shunt side, bus
voltage) and an SSSC (series side, branch flow) behind one DC link.
Restricting the injected series voltage to quadrature with the line
current makes the series converter exchange (approximately) no active
power with the line, and what remains is an SSSC on the branch plus a
STATCOM at the bus. `addUpfcControl!` (YAML type `upfc`) registers the
pair as one named device:

- one call, one composite name; the series controller steers the branch
  active power inside `v_inj_max_pu`, the shunt controller holds a remote
  bus voltage inside `s_max_mva` (or `i_max_ka`);
- registration is all-or-nothing, and the composite equals the manually
  registered pair to machine precision;
- the result table keeps one row per actuator with `at_limit` per
  converter side; both rows carry the device string
  `UPFC series/shunt (VSC pair, stationary quadrature model)`.

The composite has no series active-power injection; independent P and Q
steering needs the full model below.

```yaml
control:
  enabled: true
  controllers:
    upfc_main:
      type: upfc                 # model: quadrature is the default
      from_bus: A
      to_bus: B
      shunt_bus: B
      target_bus: LOAD
      target_vm_pu: 1.0
      p_target_mw: 35.0
      v_inj_max_pu: 0.05
      s_max_mva: 25.0
```

!!! details "Why the quadrature restriction"
    The DC link couples the two converters through an active-power
    balance: a two-actuator controller with one coupling constraint,
    which the single-actuator secant pattern of the outer loop does not
    fit. Quadrature injection removes the coupling, so each side runs as
    the existing single-actuator controller.

## UPFC: the full DC-link-coupled model

`model = :full` (YAML `model: full`) adds the phase-shifter degree of
freedom: a series voltage
`V_se` of arbitrary phase, so the line carries independent active and
reactive targets at once. The active part of the series injection,
`P_se = Re(V_se·conj(I_s))`, flows through the DC link and is balanced by the
shunt converter (`P_sh = -P_se`); in quadrature it is zero and the device
collapses onto the composite above:

```text
        Im (quadrature to I_s: reactance, NO DC power)
         ^
         |      V_se
         |     /
         |    /  in-phase part -> P_se -> DC link -> shunt   (the UPFC DOF)
         |   /
         +--------------->  Re, aligned with the line current I_s
```

The series source is an equivalent series impedance `z_add = V_se / I_s`
on the branch (`Re(z_add) < 0` when the converter injects active power),
so the line stays an ordinary branch. With the terminal voltages frozen
each outer iteration the from-end flow is affine in `V_se`, so the series
step is an exact 2x2 solve; the coupled iteration uses an adaptive damping
line search.

- one call, one controller: `model = :full` steers the from-end line flow
  to `p_target_mw` and `q_target_mvar`;
- the shunt converter provides the DC-link balance plus a reactive setpoint
  `q_shunt_mvar`, inside the current-based rating with headroom
  `Q_max = sqrt((V·s_max)^2 - P_sh^2)`;
- the result row carries `V_se` magnitude and angle, `P_se`, `P_sh`, the
  shunt Q, and the DC-link residual `|P_se + P_sh|` (a convergence quantity:
  the balance holds by construction only at the frozen state);
- `series_phase = :quadrature` forces `P_se = 0` and reproduces the composite.

```yaml
control:
  enabled: true
  controllers:
    upfc_full:
      type: upfc
      model: full
      from_bus: I
      to_bus: J
      shunt_bus: I
      p_target_mw: 40.0
      q_target_mvar: 10.0
      q_shunt_mvar: 0.0
      v_inj_max_pu: 0.20
      s_max_mva: 120.0
```

**Limitations**

- **Stationary model.** No dynamics; IPFC (a shared DC bus across several
  lines) is out of scope.
- **Shunt reactive setpoint, not closed-loop voltage.** Coupling a
  shunt-voltage secant with the line reactive-flow control does not
  converge in the sequential outer loop; closed-loop shunt voltage
  regulation would need power-flow sensitivities or an augmented in-solver
  state.
- **No explicit series current limit.** Only `|V_se| <= v_inj_max_pu` is
  clamped, not the series-converter current `|I_s| <= i_max`.
- **The branch impedance is modified in place** (like SSSC/TCSC): `z_add`
  stays on the branch after the control run, with a negative resistance
  part for the full model. That is the power-flow construct, not the
  physical line, so every branch carries its physical base impedance
  (`r_base_pu`/`x_base_pu`) next to the live value, and `runShortCircuit!`
  and the CGMES/MATPOWER exports read the base without a manual reset.
  `restoreBaseImpedances!(net)` returns the live field to the base,
  `clearUpfcFullControllers!(net)` does the same and drops the controller.
- **Low line current.** `z_add = V_se / I_s` is floored at a minimum
  current, so a dead line stays finite and keeps its base impedance.
- **Convergence regime.** Feasible, moderate flow targets converge
  reliably; aggressive targets near the injectable-voltage limit may not.

## Examples

`examples/others/exp_facts_limit_modes.jl` shows the three limit
characteristics on one weak corridor plus the SSSC window on a loop
network; `examples/others/exp_facts_base_impedance.jl` runs a full UPFC on
a meshed corridor, then short circuit and MATPOWER export on the base
impedance. Chapter 4 of the advanced workshop tour walks the same
contrasts.

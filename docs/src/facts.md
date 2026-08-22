# [FACTS Devices](@id facts_devices)

FACTS (Flexible AC Transmission Systems) devices use power electronics to
control quantities that a classical grid can only influence indirectly:
bus voltage, branch flow, and the sharing of power between parallel paths.
In a steady-state power flow a FACTS device reduces to a controllable
network parameter, a reactive injection, a shunt susceptance, a series
reactance, or a converter power, driven toward a target by the
[generic outer control loop](control_framework.md) around `runpf!`.

This page collects the device family, the limit characteristics that
distinguish the devices from each other, and the modeling decisions behind
Sparlectra's implementations. Hands-on: chapter 9 of the
[workshop tour](generated/workshop_tour.md) and the example
`examples/others/exp_facts_limit_modes.jl`.

## Device family and where each one lives

| Device | Kind | Actuator | Sparlectra controller | Theory |
|---|---|---|---|---|
| SVC (thyristor-controlled reactor/capacitor) | shunt | susceptance $B$ | `addShuntVoltageControl!` | this page, [Control Framework](control_framework.md) |
| STATCOM (VSC shunt converter) | shunt | reactive current | `addMachineVoltageControl!` with `s_max_mva` | this page, [Remote Voltage Control](remote_voltage_control.md) |
| TCSC (thyristor-controlled series capacitor) | series | reactance $x$ in a fixed window | `addSeriesReactanceControl!` with `x_min_pu`/`x_max_pu` | [Series Compensation](series_compensation.md) |
| SSSC (VSC series converter) | series | reactance deviation, voltage-bounded | `addSeriesReactanceControl!` with `v_inj_max_pu` | this page, [Series Compensation](series_compensation.md) |
| PST / Schrägregler (phase-shifting transformer) | phase | tap angle | `addPowerTransformerControl!` (`mode = :branch_active_power`) | [Control Framework](control_framework.md) |
| HVDC back-to-back (paired VSC/LCC converters) | converter | paired P injections, Q or voltage per terminal | `addHvdcPairControl!` | [HVDC Back-to-Back](hvdc_back_to_back.md) |
| UPFC (combined shunt + series converter) | combined | two coupled actuators | not implemented, see the design note below | this page |

All controllers report through the same surfaces: `ControlRunResult`,
`controllableElements` (element, device, actuator with live range, target,
`status`/`converged`/`at_limit`), and the per-controller summary printers.

## The limit characteristic is the device

In range, every shunt compensator does the same thing: it holds a voltage
target by injecting reactive power. The devices differ at their LIMIT, and
the limit is precisely where the distinction matters, because a compensator
reaches it during the depressed-voltage conditions it was installed for.

**Synchronous machine (constant-Q box).** The classical limit is a fixed
reactive capability, $Q \in [Q_{min}, Q_{max}]$, independent of the
terminal voltage. At the limit the machine behaves like a fixed injection;
this is the `MachineVoltageControl` default (`limit_mode = :constant_q`)
and the exact outer-loop analogue of PV to PQ switching.

**SVC (constant-B limit).** An SVC regulates a continuous susceptance $B$
within $[B_{min}, B_{max}]$. At the limit the susceptance clamps and the
delivered reactive power follows the voltage QUADRATICALLY through the
Y-bus stamp:

```math
Q_{SVC} = V^2 \, B_{lim}
```

This is the well-known weakness of the SVC: as the voltage sags, the
support collapses with $V^2$. At $V = 0.9$ pu a fully switched-in SVC
delivers only 81 percent of its nominal rating. Sparlectra's
`ShuntVoltageControl` produces this behavior with no extra modeling: the
clamped susceptance stays stamped in the Y-bus and the $V^2$ law falls out
of the power flow itself.

**STATCOM (constant-current limit).** A STATCOM is a voltage-source
converter; its bound is the converter CURRENT. The deliverable reactive
power therefore scales LINEARLY with the terminal voltage:

```math
Q_{STATCOM} = V \, I_{max} \quad\text{, in Sparlectra: } Q_{lim} = V \cdot S_{max}
```

with $S_{max}$ the rating at 1.0 pu. At $V = 0.9$ pu the STATCOM still
delivers 90 percent of its rating, the decisive advantage over the SVC
under exactly the sag conditions that matter. `addMachineVoltageControl!`
with `s_max_mva` (or `i_max_ka`, converted via $\sqrt{3}\,U_n I_{max}$)
switches the machine controller into this mode: the symmetric bound
$\pm V \cdot S_{max}$ is re-evaluated from the solved terminal voltage
before every outer step, so an at-limit STATCOM TRACKS the sagging or
recovering voltage instead of freezing at a stale bound.

Summarized on one axis (delivered reactive power at the capacitive limit,
relative to the 1.0-pu rating):

| Terminal voltage | Machine box | STATCOM ($\propto V$) | SVC ($\propto V^2$) |
|---:|---:|---:|---:|
| 1.00 pu | 100 % | 100 % | 100 % |
| 0.95 pu | 100 % | 95 % | 90 % |
| 0.90 pu | 100 % | 90 % | 81 % |
| 0.80 pu | 100 % | 80 % | 64 % |

The machine column is idealized (a real machine derates too, through its
capability curve); the STATCOM/SVC columns are the model behavior and match
the device physics to first order.

## Series side: fixed window versus voltage-bounded window

**TCSC (fixed reactance window).** The thyristor-controlled series
capacitor changes the branch reactance within a hardware-defined window
$[x_{min}, x_{max}]$, independent of loading. At a window end the branch is
a fixed compensated line. The resonance region between capacitive and
inductive operation is excluded by the impedance-magnitude guard `eps_z`
(see [Series Compensation](series_compensation.md)).

**SSSC (injected-voltage window).** The static synchronous series
compensator injects a voltage in quadrature with the line current. In
steady state that is equivalent to a reactance DEVIATION from the natural
line reactance $x_{base}$, bounded by the injectable voltage magnitude:

```math
|V_{inj}| = |I| \cdot |x - x_{base}| \le V_{inj,max}
\quad\Longleftrightarrow\quad
|x - x_{base}| \le \frac{V_{inj,max}}{|I|}
```

The usable window is therefore CURRENT-dependent and shrinks with loading:
at high transfer, exactly when a large flow correction would need a large
reactance swing, the SSSC saturates, while a TCSC keeps its full window.
Conversely, at light loading the SSSC window is wide. In
`addSeriesReactanceControl!` the mode is selected with `v_inj_max_pu`; the
window $x_{base} \pm V_{inj,max}/|I|$ is re-evaluated from the solved
branch current before every outer step, with a floor on $|I|$ (a
currentless branch is physically unconstrained) and the same `eps_z`
resonance guard applied as a clamp.

## Live bounds in the outer loop

Both converter-based modes (STATCOM, SSSC) share one mechanism, the LIVE
BOUND: the actuator range is a function of the solved operating point and
is refreshed at the start of every outer iteration, before the secant step
is clamped against it. Two consequences:

- **At-limit tracking.** A parked controller whose bound still moves (the
  voltage keeps sagging, the current keeps rising) is released and keeps
  adjusting; it reports `at_limit` only once its bound has settled. The
  delivered quantity therefore follows the physical limit law across outer
  iterations instead of freezing at the first clamp.
- **Honest element rows.** `controllableElements` and the report rows show
  the bounds of the LAST evaluated operating point, so `actuator_min`/
  `actuator_max` are the currently deliverable range, not the nameplate.

In range, both modes behave like their fixed-limit counterparts: the same
secant iteration on the same scalar map, converging into the same deadband.

## Usage

Programmatic (see the docstrings for the full keyword sets):

```julia
# SVC: continuous susceptance, quadratic limit collapse
addShuntVoltageControl!(net; bus = "B", target_vm_pu = 1.0,
                        bs_min_mvar = -60.0, bs_max_mvar = 60.0)

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

## UPFC: design status

The unified power flow controller combines a STATCOM (shunt side, bus
voltage) and an SSSC (series side, branch P and Q) behind one DC link; the
link couples the two converters through an active-power balance. That makes
it a TWO-actuator controller with one coupling constraint, which does not
fit the single-actuator secant pattern the outer loop is built on: each
Sparlectra controller today owns one actuator and one target, and the loop
coordinates controllers only through the shared power flow.

The UPFC therefore stays UNIMPLEMENTED by decision (issue #297 Draft G).
The stationary approximation that IS expressible today: place a PST (angle),
a series controller (TCSC or SSSC mode), and a STATCOM on the same corridor.
This reproduces the three control channels of a UPFC (angle, series
compensation, shunt voltage) but NOT the DC-link power balance between the
shunt and series converters; the three controllers converge through the
outer loop as independent devices. For studies where the coupling matters,
the approximation overestimates the device's degrees of freedom, and the
honest answer is that a dedicated multi-actuator controller with a coupling
constraint is future work.

## Validation

- `test/test_tap_controller.jl`: STATCOM registration validation, limit
  tracking on a depressed-voltage corridor ($Q = V \cdot S_{max}$ across
  operating points), in-range equivalence with the constant-Q mode, and
  the SVC-versus-STATCOM contrast on one case ($V^2$ versus $V$ collapse,
  asserted numerically).
- `test/test_series_reactance_control.jl`: SSSC registration validation,
  converged operation inside the live window, pinned operation with the
  effective injected voltage at $V_{inj,max}$, and TCSC-mode regression.
- `examples/others/exp_facts_limit_modes.jl`: the three limit
  characteristics side by side on one weak corridor plus the SSSC window
  on a loop network.
- Workshop tour chapter 9 walks the same contrasts interactively.

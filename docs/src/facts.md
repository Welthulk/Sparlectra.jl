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
Sparlectra's implementations. Hands-on: chapter 4 of the
[advanced workshop tour](generated/workshop_tour_advanced.md) and the
example `examples/others/exp_facts_limit_modes.jl`.

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
| UPFC (combined shunt + series converter) | combined | series voltage (quadrature composite, or arbitrary-phase full model) + shunt | `addUpfcControl!` (`model = :quadrature` or `:full`, see the notes below) | this page |

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

## Discrete banks: MSC/MSR

Strictly speaking a switched capacitor/reactor bank is not power
electronics, but it is the working end of most voltage-control schemes and
shares the shunt physics above, so it lives in the same controller
(issue #324). With `step_mvar` the susceptance moves in whole switched
blocks (e.g. four times 10 MVAr):

- the secant proposal is TRUNCATED toward the target to whole steps, so
  the bank approaches the voltage target from one side and never
  overshoots; this is the anti-hunting guarantee, a bank that never
  crosses its target cannot oscillate between two adjacent blocks;
- when no whole block improves the voltage further, the controller PARKS
  on the reached step (`status = :parked`, deliberately the last step
  BEFORE crossing: conservative under-compensation instead of a possible
  overvoltage) and releases itself when another controller moves the
  operating point far enough that a block helps again;
- at the outermost admissible block the constant-B limit region applies
  unchanged, the delivered Q follows $V^2$ with the last block connected.

The classical coordination case, a switched bank plus an OLTC on one bus,
remains open follow-up work.

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
voltage) and an SSSC (series side, branch flow) behind one DC link; the
link couples the two converters through an active-power balance. That makes
the full device a TWO-actuator controller with one coupling constraint,
which does not fit the single-actuator secant pattern the outer loop is
built on: each Sparlectra controller owns one actuator and one target, and
the loop coordinates controllers only through the shared power flow. That
was the reason for deferring the UPFC in issue #297 Draft G, and it still
holds for the full device.

The way in is the QUADRATURE argument (issue #325): restrict the injected
series voltage to quadrature with the line current. Then the series
converter exchanges (approximately) no active power with the line, the DC
link carries about zero, the coupling constraint degenerates, and what
remains is exactly an SSSC on the branch plus a STATCOM at the bus. Both
controllers exist, so `addUpfcControl!` registers them together as one
named device:

- one call, one composite name; the series controller steers the branch
  active power inside the injected-voltage limit `v_inj_max_pu`, the shunt
  controller holds a remote bus voltage inside the current-based rating
  `s_max_mva` (or `i_max_ka`);
- registration is all-or-nothing (a rejected call leaves the net
  untouched), and the composite behaves exactly like the manually
  registered pair, to machine precision;
- the result table keeps one row per actuator with `at_limit` per converter
  side; both rows carry the device string
  `UPFC series/shunt (VSC pair, stationary quadrature model)`.

What the composite is NOT: it has no series ACTIVE-power injection. The
phase-shifter degree of freedom stays unavailable, and independent P and Q
steering of the line needs the full model below.

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

## UPFC: the full DC-link-coupled model

The full model (issue #326, `model = :full`) delivers the phase-shifter
degree of freedom: a series voltage `V_se` of ARBITRARY phase, so the line
carries INDEPENDENT active and reactive targets at once. The active part of
the series injection, `P_se = Re(V_se·conj(I_s))`, flows through the DC link
and is balanced by the shunt converter (`P_sh = -P_se`). In quadrature the
in-phase component is zero and the device collapses onto the composite
above; the picture is the split of the injected voltage relative to the line
current:

```text
        Im (quadrature to I_s: reactance, NO DC power)
         ^
         |      V_se
         |     /
         |    /  in-phase part -> P_se -> DC link -> shunt   (the UPFC DOF)
         |   /
         +--------------->  Re, aligned with the line current I_s
```

The series source is realised as an equivalent series impedance
`z_add = V_se / I_s` added to the branch (`Re(z_add) < 0` when the converter
injects active power), so the line stays an ordinary branch and no fictitious
injection is created. With the terminal voltages frozen each outer iteration
the from-end flow is affine in `V_se`, so the series step is an exact 2x2
solve; the coupled iteration is globalised with an adaptive damping line
search.

- one call, one controller (not a pair): `model = :full` steers the from-end
  line flow to `p_target_mw` AND `q_target_mvar`;
- the shunt converter provides the DC-link balance plus a reactive SETPOINT
  `q_shunt_mvar`, inside the current-based rating whose reactive headroom is
  coupled to the active load, `Q_max = sqrt((V·s_max)^2 - P_sh^2)`;
- the result row carries `V_se` magnitude and angle, `P_se`, `P_sh`, the
  shunt Q, and the DC-link residual `|P_se + P_sh|` (a genuine convergence
  quantity, since the balance holds by construction only at the frozen state);
- `series_phase = :quadrature` forces `P_se = 0` and reproduces the composite.

Honest limitations of the first cut:

- **Stationary model.** No dynamics or transients; IPFC (a shared DC bus
  across several lines) is out of scope.
- **Shunt reactive SETPOINT, not closed-loop voltage.** Coupling a
  shunt-voltage secant with the line reactive-flow control does not converge
  in the sequential outer loop (a known behaviour of injection-model UPFCs);
  closed-loop shunt voltage regulation is a follow-up that needs the AC
  power-flow sensitivity framework (issue #217) or an augmented in-solver
  state.
- **No explicit series current limit.** Only the injected-voltage magnitude
  `|V_se| <= v_inj_max_pu` is clamped, not the series-converter current
  `|I_s| <= i_max`; adding it is one more clamp on the same step.
- **The branch impedance is modified in place** (like the SSSC/TCSC): the
  equivalent series impedance `z_add` stays on the branch after the control
  run, and for the full model its resistance part goes NEGATIVE. That is the
  correct steady-state power-flow construct, but not the physical line, so a
  short circuit or an export run on the SAME net must use the base network:
  `runShortCircuit!` FAILS LOUDLY on a negative-resistance branch (the
  converter is bypassed under fault), and a CGMES/MATPOWER export would write
  the compensated impedance, not the physical one. Run the fault calculation
  or the export on the unmodified network, or before the control run.
- **Low line current.** `z_add = V_se / I_s` is floored at a minimum current
  (`|I_s|` guard) so a lightly loaded or dead line stays finite and keeps its
  base impedance; a UPFC on an essentially currentless line has nothing to
  steer.
- **Convergence regime.** The full model converges reliably for feasible,
  moderate flow targets (the realistic operating envelope of a UPFC). Very
  aggressive targets near the injectable-voltage limit may not converge in
  the outer loop; a robust envelope needs the sensitivity/Jacobian work above.

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

## Validation

- `test/test_tap_controller.jl`: STATCOM registration validation, limit
  tracking on a depressed-voltage corridor ($Q = V \cdot S_{max}$ across
  operating points), in-range equivalence with the constant-Q mode, and
  the SVC-versus-STATCOM contrast on one case ($V^2$ versus $V$ collapse,
  asserted numerically).
- `test/test_series_reactance_control.jl`: SSSC registration validation,
  converged operation inside the live window, pinned operation with the
  effective injected voltage at $V_{inj,max}$, and TCSC-mode regression.
- `test/test_upfc_control.jl`: the quadrature composite equals the manually
  registered SSSC+STATCOM pair to machine precision, both limit
  characteristics at their clamps, all-or-nothing registration, and the YAML
  type `upfc` with the double-apply no-op; the full model reaches independent
  P and Q on one line simultaneously with the DC-link balance closed, reduces
  to the SSSC when the series phase is forced to quadrature, and round-trips
  through `model: full` in YAML.
- `examples/others/exp_facts_limit_modes.jl`: the three limit
  characteristics side by side on one weak corridor plus the SSSC window
  on a loop network.
- Chapter 4 of the advanced workshop tour walks the same contrasts
  interactively.

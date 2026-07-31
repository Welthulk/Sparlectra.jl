# Analysis: #227 — Generic framework for controllable grid elements

Status: analysis, 2026-07-31. Question set: what does Sparlectra already
have, what would the generic layer add, what are FACTS devices, and how are
they modeled in a steady-state power flow?

> **Implementation status:** increments 1 (generic reporting view,
> `controllableElements`) and 2 (SVC susceptance mode,
> `ShuntVoltageControl`) are implemented. Ready-to-post issue texts for
> increments 3–6 live in
> `docs/dev/issue_drafts_controllable_elements.md`.

## 1. What #227 asks for

One internal abstraction for everything that follows the pattern
*actuator → controlled quantity → target → limits → status*, so that future
devices (SVC, STATCOM, TCSC, back-to-back HVDC, …) do not each become a
special case. First scope per the issue: internal types, a controller
registry, reporting, one simple device (SVC/STATCOM-like), outer-loop
coupling only.

## 2. What already exists (inventory)

The outer-loop control framework is in place and already covers most of the
issue's "suggested first scope":

| #227 first-scope item | Status today |
|---|---|
| 1. Internal types for actuator/target/limits/status | **Largely present, implicitly.** `AbstractOuterController` with the lifecycle `control_evaluate` → `control_apply_update!` → `control_report_rows`/`control_trace_rows`; targets, deadbands, limits, and an honest `at_limit`/`converged` status live on the two concrete controllers. What is missing is a *shared, explicit* data model (see gap analysis). |
| 2. Controller registry on the network | **Done.** Transformer controllers live on the windings (`PowerTransformerWinding.controls`), machine controllers on `net.machineControls`; `collect_outer_controllers(net)` is the single collection point. |
| 3. Reporting for configured/active controllers | **Done.** `ControlRunResult` / `latest_control_result(net)`, per-controller report and trace rows, WebUI/`run.log` integration. |
| 4. First simple device (SVC/STATCOM-like) | **The actual gap** — see below. |
| 5. Outer-loop coupling, no Newton extension | **Done.** `run_control!` solves → evaluates → applies → re-solves, with `calcNetLosses!` refreshes between iterations. Discrete and continuous actuators both work. |

Concrete controllable elements today:

- **OLTC transformer** (`PowerTransformerControl`, `mode = :voltage`):
  ratio-tap actuator, local or remote bus-voltage target, discrete or
  continuous steps, tap limits, probe-based direction estimate.
- **PST** (`mode = :branch_active_power`): phase-tap actuator, branch-P
  target; since #274 the series reactance follows the device characteristic
  X(α) on every accepted move (typed `PhaseTapChangerModel` on the winding,
  formula or tabular), and the probe perturbs angle *and* reactance
  consistently. Split "Schrägregelung" (two disjoint actuators on one
  transformer) is supported.
- **Machine remote voltage control** (`MachineVoltageControl`): machine-Q
  actuator via secant iteration onto a remote PQ bus target, honest
  `at_limit` on the reactive bounds, one controller per target bus enforced,
  cross-type warning when a tap controller regulates the same bus.
- **Static voltage holding** (not outer-loop, but part of the picture):
  PV/slack machines with Q limits (active-set or classic enforcement),
  voltage-dependent controls, and the imported **SVCs from CGMES**, which
  today are mapped as *static* PV injections with rating-derived Q limits —
  they hold their setpoint through the PV machinery, not through a
  controller object.
- **Stage-0 HVDC converters** from CGMES multi-area imports: fixed PCC
  injections per side (no coupling constraint, no controller).

So the user's observation is exactly right: the framework is there. What
#227 adds on top is (a) a *shared vocabulary* over the existing controllers
and (b) the missing device families.

## 3. What are FACTS?

**FACTS — Flexible AC Transmission Systems** — are power-electronics-based
devices (thyristor or voltage-source-converter technology) that make
quantities controllable which are fixed in a conventional AC grid: nodal
voltage, effective series impedance of a line, and the voltage angle across
it. Because transmitted active power over a corridor follows
`P ≈ V1·V2·sin(δ1−δ2)/X`, controlling `V`, `X`, or `δ` means controlling
the power flow itself — without building new lines. The classic families:

- **Shunt devices** (voltage support at a node):
  - *SVC* (Static Var Compensator): thyristor-switched/controlled reactors
    and capacitors — a **variable shunt susceptance** `B` within
    `[B_min, B_max]`.
  - *STATCOM* (Static Synchronous Compensator): a voltage-source converter
    behaving like a synchronous condenser — a **controllable reactive
    current source**.
- **Series devices** (flow control on a line):
  - *TCSC* (Thyristor-Controlled Series Capacitor): a **variable series
    reactance** `X_TCSC` within `[X_min, X_max]` inserted into the line
    (capacitive to reduce, inductive to increase the effective X; a
    resonance band in between is excluded).
  - *SSSC* (Static Synchronous Series Compensator): a converter injecting a
    **series voltage in quadrature** with the line current — compensation
    independent of loading.
- **Combined**: *UPFC* (Unified Power Flow Controller) = STATCOM + SSSC on
  a common DC link — controls the line's P and Q and the local voltage
  simultaneously (explicitly a non-goal in #227).
- **Related converter technology**: *HVDC*, especially **back-to-back**
  links (two converter stations without a DC line) that couple two AC areas
  with a fully controlled power transfer.

The phase-shifting transformer is the electromechanical ancestor of the
series family — which is why Sparlectra's PST controller is the right
template for TCSC-like devices.

## 4. How the modern devices are modeled in a steady-state power flow

The industry-standard steady-state treatments, per device:

- **SVC — two operating regions.** In range, the SVC holds its voltage
  setpoint: modeled as a **PV bus with Q limits** derived from the
  susceptance range (`Q = B·V²`, so `Q_lim = B_lim·V²`). At a limit the
  device becomes a **fixed susceptance**: its reactive output then varies
  with `V²` instead of staying constant. That second region is the real
  modeling difference from a generator — a Q-limited generator holds
  constant Q, a limited SVC does not. Sparlectra today implements the
  first region (PV + rating-derived Q limits, with the SV-consistency
  demotion guard); the constant-B limit region is not modeled.
- **STATCOM — current-limited, not power-limited.** In range identical to
  the SVC (PV bus). At its limit a STATCOM injects its maximum reactive
  **current**, so the Q limit scales with `V` (`Q_lim = V·I_max`) rather
  than `V²` (SVC) or constant (machine). In practice many tools
  approximate it as a machine with constant Q limits — that is exactly what
  Sparlectra's existing machine model already provides.
- **TCSC — adjustable series reactance.** The standard model is simply a
  branch whose `x` is an actuator within `[X_min, X_max]`, steered onto
  a branch-flow target; an alternative formulation replaces the device by
  equivalent nodal power injections. The outer-loop mechanics — "accepted
  move updates `x_pu`, next solve re-stamps the Y-bus" — are precisely what
  #274 just built for PST reactance tracking.
- **SSSC — series voltage source.** Steady state: an ideal series voltage
  `V_se` perpendicular to the line current, with magnitude limits, usually
  folded into an injection model. More involved than TCSC because it adds a
  genuine new variable.
- **UPFC** — combined shunt + series source with an internal active-power
  balance constraint between the two converters; the standard PF treatment
  either augments the Newton system or uses a two-source injection model.
  Rightly excluded from the first scope.
- **Back-to-back HVDC** — two paired injections with a coupling constraint
  `P1 + P2 + P_loss = 0` and per-terminal Q/V control. Sparlectra's
  Stage-0 CGMES import already creates the two fixed injections; a
  controller would enforce the pairing and steer the transfer.

## 5. What the generic layer buys (Vorteile)

1. **No new special cases.** SVC/STATCOM/TCSC each reduce to "existing
   lifecycle + one new actuator": variable shunt-B, machine-Q with a
   V-dependent limit, series-X on a line. The framework, registry,
   reporting, deadband/limit semantics, and WebUI surfaces are reused as-is.
2. **Uniform operational semantics.** `at_limit` honesty, one-actuator-one-
   controller exclusivity, cross-controller interaction warnings, and trace
   rows already exist once — every new device inherits them instead of
   re-inventing them (the tap/machine cross-type warning shows why that
   matters).
3. **A named vocabulary enables tooling.** An explicit
   `(element, actuator, quantity, target, limits, mode, status)` record —
   even as a thin reporting view over the existing controllers — is the
   precondition for the issue's sensitivity perspective:
   `dx/du = -J⁻¹ · ∂r/∂u` needs to know *which* u belongs to *which*
   controller. The Jacobian exists in the solver; the association layer is
   what is missing.
4. **Roadmap coherence.** TCSC via the #274 mechanics, STATCOM via the
   machine controller, SVC-B via a new shunt actuator, B2B-HVDC via paired
   injections — each is an incremental step on the same rail instead of
   four unrelated features.

## 6. Recommended increments (smallest first)

1. **Reporting view first (cheap):** derive the generic
   controllable-element record from the existing controllers and surface it
   in `ControlRunResult`/WebUI — the issue's items 1–3 become *explicit*
   without touching any control logic.
2. **First device: SVC in susceptance mode** (`ShuntVoltageControl`):
   actuator B on the existing shunt model, voltage target, at-limit
   behavior "freeze B" (constant-B region) — upgrades the imported CGMES
   SVCs from static PV mapping to honest device behavior at their limits.
3. **TCSC-like series-X controller** on a line branch — reuses the #274
   update path (`x_pu` assignment + Y-bus re-stamp) and the PST probe.
4. **STATCOM variant** of `MachineVoltageControl` with `Q_lim = V·I_max`
   evaluated per outer iteration.
5. **B2B pairing controller** over the Stage-0 converter injections
   (transfer target, per-terminal Q/V).
6. **Sensitivities** as a separate follow-up once the association layer
   exists (solver-side `J` access plus `∂r/∂u` per actuator kind).

Non-goals stay as in the issue: no UPFC first, no full AC/DC model, no
migration of the existing tap/PV logic.

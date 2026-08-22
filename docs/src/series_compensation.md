# Series Compensation (TCSC)

A TCSC (thyristor controlled series capacitor) is a variable series
reactance inserted into a line branch. Steered onto a branch active-power
target, it redistributes flow between the parallel paths of a meshed
network. Sparlectra models it as the outer-loop controller
`SeriesReactanceControl`, registered with `addSeriesReactanceControl!`
(issue #297). This page explains why a series reactance steers flow, how
the physical device maps onto the model, and how the controller iterates;
the outer-loop machinery itself is described in the
[Control Framework](control_framework.md).

## The branch as a four-terminal element

Every branch enters the power flow through its admittance matrix (see
[Branch Model](branchmodel.md) for the derivation and sign conventions):

```math
Y_{br} = \begin{bmatrix}
  \frac{1}{\tau^2}\left(y_{ser} + \frac{y_{shunt}}{2}\right) &
  -y_{ser}\,\frac{1}{\tau e^{-j\phi}} \\
  -y_{ser}\,\frac{1}{\tau e^{j\phi}} &
  y_{ser} + \frac{y_{shunt}}{2}
\end{bmatrix},
\qquad y_{ser} = \frac{1}{R + jX},
```

with the complex turns ratio $N = \tau e^{j\phi}$ and $N = 1$ for lines.
The TCSC acts purely through $X$ inside $y_{ser}$: the controller assigns a
new `x_pu` to its line branch, which changes exactly one branch stamp, and
the outer loop re-stamps the Y-bus before the next solve. No solver
internals are touched; this is the same actuator mechanism the
phase-shifting transformer uses for its tap-dependent reactance
$X(\alpha)$ in the [tap controller](control_framework.md).

## Why a series reactance steers flow

For a lossless line between buses 1 and 2 the transferred active power is

```math
P_{12} = \frac{V_1 V_2}{X}\,\sin(\delta_1 - \delta_2) .
```

Lowering the series reactance raises the power carried by that path at a
given angle difference. In a meshed network the effect appears as a flow
redistribution: two parallel corridors split the transfer in inverse
proportion to their reactances, so compensating one corridor pulls flow
onto it and off its neighbors. This loop-network mechanism is exactly what
the example and the tests exercise: a corridor with twice the reactance of
its parallel path carries one third of the transfer at baseline, and the
TCSC moves that share onto the target.

## Compensation degree

Series compensation is usually quoted as the compensation degree

```math
k = \frac{X_C}{X_{line}}, \qquad
X = X_{line} - X_C = (1 - k)\,X_{line},
```

with $X_C$ the inserted capacitive reactance and $X$ the net branch
reactance the power flow sees. Practical installations run at roughly
$k = 0.2$ to $0.7$; beyond that, subsynchronous resonance concerns and
protection coordination dominate. The model does not track $k$ explicitly:
the controller works directly on the net reactance `x_pu`, clamped to
`[x_min_pu, x_max_pu]`. A range for a line with $X_{line} = 0.20$ p.u. and
up to 70 percent compensation is therefore simply
`x_min_pu = 0.06, x_max_pu = 0.20`.

## Physical device versus model

A real TCSC is a firing-angle controlled parallel circuit of a capacitor
bank and a thyristor-switched reactor. Its apparent reactance as a function
of the firing angle passes through a resonance region between the
capacitive and the inductive operating range; a real device must not dwell
there, and vendor curves exclude it.

The model is deliberately simpler: one continuous, clamped reactance. Two
consequences follow:

- The impedance-magnitude guard `eps_z` is the documented stand-in for the
  resonance exclusion. Registration rejects any range in which the series
  impedance magnitude $|R + jX|$ falls below the guard, checked at both
  range ends and, for a sign-crossing range, at the crossing itself (there
  the magnitude bottoms out at $|R|$). Splitting the admissible range into
  two disjoint intervals (capacitive and inductive, with the resonance band
  excluded between them) is a possible later refinement, mirroring the
  issue text.
- Negative net reactance (a net capacitive branch) is admissible by design:
  the guard already protects the singular neighborhood of $X = 0$, so the
  model does not need a positivity restriction. This is a documented design
  choice; most practical targets are reachable well inside the inductive
  range.

At a range end the controller stops moving and the branch behaves as a
fixed compensated line; there is no special casing beyond the honest
`at_limit` report.

The measured quantity is the active power of the controlled branch itself,
in the registered `fromBus` to `toBus` direction (the tap controller's
`achieved_p_mw` convention). Steering a remote branch through a TCSC
elsewhere is a conceivable extension (an optional `target_branch`), but the
first cut keeps measurement and actuator on the same element.

## Numerical method

The controller runs inside the generic outer loop
([Control Framework](control_framework.md)): measure after a converged
solve, propose a step, apply it, re-solve.

**Secant iteration.** The scalar map from `x_pu` to the branch power is
smooth in the operating range, so the controller uses the same secant
update as the machine and shunt controllers: from the last two points
$(x_k, P_k)$ and $(x_{k-1}, P_{k-1})$ it estimates the sensitivity and
steps toward the target.

**Bootstrap step.** Before two points exist there is no measurable
sensitivity. The first move is a bounded probe: a fixed fraction of the
headroom toward the more distant range end. Unlike the shunt controller,
the sign of $dP/dX$ is not hard-coded; in a meshed network it depends on
where the branch sits relative to its parallel paths, so the probe
measures it and the secant uses the signed slope from then on.

**Sensitivity guard.** A slope magnitude below a minimum threshold would
turn the secant step into noise amplification; the controller then falls
back to the bounded probe instead, the same guard the machine controller
uses.

**Honest at limit.** Every proposed step is clamped to
`[x_min_pu, x_max_pu]`. `at_limit` is set only when the reactance is
clamped at a range end and the target is still outside the deadband:
the controller does not pretend convergence, and the solve itself remains
valid with the branch as a fixed compensated line.

## SSSC mode: injected-voltage limit (issue #297 Draft F)

The fixed window above models a TCSC, whose hardware defines an admissible
reactance range independent of loading. A static synchronous series
compensator (SSSC) is a voltage-source converter injecting a voltage in
quadrature with the line current; in steady state that acts as a reactance
deviation from the natural line reactance $x_{base}$, bounded by the
injectable voltage magnitude:

```math
|V_{inj}| = |I| \cdot |x - x_{base}| \le V_{inj,max}
\quad\Longleftrightarrow\quad
|x - x_{base}| \le \frac{V_{inj,max}}{|I|}
```

`addSeriesReactanceControl!` with `v_inj_max_pu` (instead of
`x_min_pu`/`x_max_pu`) switches the controller into this mode. The same
secant machinery runs on a LIVE window:

- the bounds $x_{base} \pm V_{inj,max}/|I|$ are re-evaluated from the
  solved branch current (measured at the registered from side,
  $|I| = |S|/V$) before every outer step;
- the window SHRINKS with loading: at high transfer, exactly when a large
  flow correction would need a large reactance swing, the SSSC saturates,
  while a TCSC keeps its full window (the defining contrast between the
  two devices, see [FACTS Devices](@ref facts_devices));
- a floor on $|I|$ keeps the window finite on a currentless branch
  (physically an SSSC injects no voltage without current), and the `eps_z`
  resonance guard is applied as a clamp on the window instead of an error,
  so a transient operating point can never abort the run;
- an at-limit SSSC whose window still moves keeps adjusting and parks
  `at_limit` only once the window has settled; at the limit the effective
  injected voltage $|I| \cdot |x - x_{base}|$ sits at $V_{inj,max}$.

The element row and the summary printer show the live window, the natural
reactance, the measured current, and the currently injected voltage.

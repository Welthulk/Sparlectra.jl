# Series Compensation (TCSC)

A TCSC (thyristor controlled series capacitor) is a variable series
reactance in a line branch. Steered onto a branch active-power target, it
redistributes flow between the parallel paths of a meshed network.
Sparlectra models it as the outer-loop controller `SeriesReactanceControl`,
registered with `addSeriesReactanceControl!` on the
[Control Framework](control_framework.md).

| Use | Where |
|---|---|
| Call | `addSeriesReactanceControl!` with `x_min_pu`/`x_max_pu` (TCSC, fixed window) or `v_inj_max_pu` (SSSC, current-dependent window); impedance guard `eps_z` |
| Result field | `achieved_p_mw` (branch power in the registered `fromBus` to `toBus` direction), `at_limit`; the element row and the summary printer show the live window, the natural reactance, the measured current and the injected voltage |

## The branch as a four-terminal element

The TCSC acts through $X$ inside the series admittance
$y_{ser} = 1/(R + jX)$ of the branch matrix ([Branch Model](branchmodel.md))
only: the controller assigns a new `x_pu` to its line branch, and the outer
loop re-stamps the Y-bus before the next solve. The tap-dependent
reactance $X(\alpha)$ of a phase-shifting transformer uses the same
mechanism ([tap controller](control_framework.md)).

## Why a series reactance steers flow

For a lossless line between buses 1 and 2

```math
P_{12} = \frac{V_1 V_2}{X}\,\sin(\delta_1 - \delta_2) .
```

A lower series reactance carries more power on that path at a given angle
difference. Parallel corridors split the transfer in inverse proportion to
their reactances, so compensating one corridor pulls flow onto it and off
its neighbors.

## Compensation degree

Series compensation is quoted as the compensation degree

```math
k = \frac{X_C}{X_{line}}, \qquad
X = X_{line} - X_C = (1 - k)\,X_{line},
```

with $X_C$ the inserted capacitive reactance and $X$ the net branch
reactance the power flow sees. Installations run at about $k = 0.2$ to
$0.7$; beyond that, subsynchronous resonance and protection coordination
dominate. The controller works on the net reactance `x_pu`, clamped to
`[x_min_pu, x_max_pu]`: a line with $X_{line} = 0.20$ p.u. and up to 70
percent compensation gets `x_min_pu = 0.06, x_max_pu = 0.20`.

## Physical device versus model

A real TCSC is a firing-angle controlled parallel circuit of a capacitor
bank and a thyristor-switched reactor; its apparent reactance passes
through a resonance region, between the capacitive and the inductive
range, where the device must not dwell. The model is one continuous,
clamped reactance. The impedance-magnitude guard `eps_z` stands in for the
resonance exclusion: registration rejects a range in which $|R + jX|$ falls
below the guard, checked at both range ends and, for a sign-crossing range,
at the crossing (where the magnitude bottoms out at $|R|$). A negative net
reactance (a net capacitive branch) is admissible; the guard protects the
singular neighborhood of $X = 0$. At a range end the controller stops and
the branch is a fixed compensated line, reported as `at_limit`.

The measured quantity is the active power of the controlled branch itself,
in the registered `fromBus` to `toBus` direction (the tap controller's
`achieved_p_mw` convention).

## Numerical method

The controller runs in the generic outer loop (measure after a converged
solve, propose a step, apply it, re-solve) with the secant update of the
machine and shunt controllers on the smooth map from `x_pu` to the branch
power. Before two points exist, the first move is a fixed fraction of the
headroom toward the more distant range end; the sign of $dP/dX$ is not
hard-coded (in a meshed network it depends on where the branch sits
relative to its parallel paths), the probe measures it, and a slope below
a minimum threshold falls back to the probe. Every step is clamped to
`[x_min_pu, x_max_pu]`; `at_limit` is set only when the reactance sits at
a range end and the target is still outside the deadband.

## SSSC mode: injected-voltage limit

The fixed window models a TCSC, whose admissible reactance range is
independent of loading. A static synchronous series compensator (SSSC)
injects a voltage in quadrature with the line current; in steady state
that is a reactance deviation from the natural line reactance $x_{base}$,
bounded by the injectable voltage:

```math
|V_{inj}| = |I| \cdot |x - x_{base}| \le V_{inj,max}
\quad\Longleftrightarrow\quad
|x - x_{base}| \le \frac{V_{inj,max}}{|I|}
```

`addSeriesReactanceControl!` with `v_inj_max_pu` (instead of
`x_min_pu`/`x_max_pu`) selects this mode. The window
$x_{base} \pm V_{inj,max}/|I|$ is re-evaluated from the solved branch
current (at the registered from side, $|I| = |S|/V$) before every outer
step, so it shrinks with loading: at high transfer the SSSC saturates while
a TCSC keeps its full window. A floor on $|I|$ keeps the window finite on a
currentless branch, and the `eps_z` guard clamps the window instead of
raising an error. At the limit $|I| \cdot |x - x_{base}|$ sits at
$V_{inj,max}$. The live-bound mechanics (adjusting while the window still
moves, `at_limit` once it has settled, the current window in the element
row) are shared with the STATCOM: [FACTS Devices](@ref facts_devices).

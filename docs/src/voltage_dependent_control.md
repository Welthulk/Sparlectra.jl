# Voltage-dependent prosumer control: Q(U) and P(U) (rectangular solver)

Voltage-dependent active/reactive power control for prosumers:

- `QUController`: `Q = f_Q(|V|)`
- `PUController`: `P = f_P(|V|)`

These are soft controls (state-dependent injections), not hard PV equality
constraints; the structural bus type (`Slack`, `PV`, `PQ`) stays unchanged.
The control is attached to the prosumer, not to bus typing: `Q(U)` droop
for volt/var support, `P(U)` for curtailed or voltage-sensitive active
injection.

Capability limits are not controls. A CGMES `ReactiveCapabilityCurve`
describes Q bounds as a function of active power Q(P): it is evaluated at
import time into ordinary `minQ`/`maxQ` limits and enforced by the
[Q-limit switching machinery](q_limit_switching_strategy.md), never through
a `QUController` (a bound as voltage-dependent setpoint would wander
with $|V|$). Both mechanisms can coexist on one machine.

## Rectangular state and voltage magnitude

```math
V_i = e_i + j f_i, \qquad |V_i| = \sqrt{e_i^2 + f_i^2}.
```

For `|V_i| > 0`,

```math
\frac{\partial |V_i|}{\partial e_i} = \frac{e_i}{|V_i|},
\qquad
\frac{\partial |V_i|}{\partial f_i} = \frac{f_i}{|V_i|}.
```

The implementation uses

```math
|V_i|_\varepsilon = \max(|V_i|, \varepsilon),\quad \varepsilon = 10^{-9}
```

in these derivatives to avoid division by zero near a collapsed voltage.

## Controlled specified injections

For a controlled prosumer at bus `i`:

```math
Q_i^{\mathrm{spec}} = f_Q(|V_i|),
\qquad
P_i^{\mathrm{spec}} = f_P(|V_i|).
```

Several prosumers on one bus sum their specified injections (generation
positive, load negative):

```math
P_i^{\mathrm{spec,tot}} = \sum_{k\in\mathcal{D}_i} P_{i,k}^{\mathrm{spec}},
\qquad
Q_i^{\mathrm{spec,tot}} = \sum_{k\in\mathcal{D}_i} Q_{i,k}^{\mathrm{spec}}.
```

## Mismatch equations in rectangular NR

For each non-slack bus:

- PQ bus:

```math
\Delta P_i = P_i^{\mathrm{calc}}(e,f) - P_i^{\mathrm{spec,tot}}(|V_i|),
\qquad
\Delta Q_i = Q_i^{\mathrm{calc}}(e,f) - Q_i^{\mathrm{spec,tot}}(|V_i|).
```

- PV bus:

```math
\Delta P_i = P_i^{\mathrm{calc}}(e,f) - P_i^{\mathrm{spec,tot}}(|V_i|),
\qquad
\Delta V_i = |V_i| - V_i^{\mathrm{set}}.
```

Only the specified injection part becomes state-dependent.

## Jacobian extension by chain rule (local terms)

`P_spec` / `Q_spec` depend only on the local `|V_i|`, so the extra Jacobian
terms are local (bus `i` row, bus `i` state columns).

Active-power mismatch row:

```math
\frac{\partial \Delta P_i}{\partial e_i}
= \frac{\partial P_i^{\mathrm{calc}}}{\partial e_i}
- \frac{dP_i^{\mathrm{spec}}}{d|V_i|}\,\frac{e_i}{|V_i|},
```

```math
\frac{\partial \Delta P_i}{\partial f_i}
= \frac{\partial P_i^{\mathrm{calc}}}{\partial f_i}
- \frac{dP_i^{\mathrm{spec}}}{d|V_i|}\,\frac{f_i}{|V_i|}.
```

Reactive-power mismatch row (PQ buses):

```math
\frac{\partial \Delta Q_i}{\partial e_i}
= \frac{\partial Q_i^{\mathrm{calc}}}{\partial e_i}
- \frac{dQ_i^{\mathrm{spec}}}{d|V_i|}\,\frac{e_i}{|V_i|},
```

```math
\frac{\partial \Delta Q_i}{\partial f_i}
= \frac{\partial Q_i^{\mathrm{calc}}}{\partial f_i}
- \frac{dQ_i^{\mathrm{spec}}}{d|V_i|}\,\frac{f_i}{|V_i|}.
```

PV second rows (`ΔV`) get no `Q(U)` term: that row is the voltage magnitude
mismatch.

## Characteristic interpolation and derivative

Controllers use ordered points `(u_k, y_k)` internally in p.u. Points are
given directly in p.u. (`voltage_unit=:pu`, `value_unit=:pu`) or in
physical units (`voltage_unit=:kV`, `value_unit=:MW`/`:MVAr`) with
conversion metadata (`vn_kV`, `sbase_MVA`).

`make_characteristic(...; interpolation = ...)` selects the interpolation:

- `:linear` (default): piecewise linear.
- `:spline`: natural cubic spline through all points.
- `:polynomial`: one global interpolating polynomial through all points.

With only two points, `:spline` and `:polynomial` reduce to a straight
line. For `:linear`, inside segment `[u_k, u_{k+1}]`:

```math
y(u) = y_k + \frac{y_{k+1}-y_k}{u_{k+1}-u_k}(u-u_k),
\qquad
\frac{dy}{du} = \frac{y_{k+1}-y_k}{u_{k+1}-u_k}.
```

Conventions for all modes:

- Below the first point: clamp to the first value, derivative `0`.
- Above the last point: clamp to the last value, derivative `0`.
- At breakpoints: segment-wise evaluation (continuous value; slope is
  side-dependent).
- At an explicit min/max saturation: output clamped, derivative `0`.

## API usage

```julia
using Sparlectra

qu = QUController(
    make_characteristic([(104.5, 30.0), (107.0, 20.0), (110.0, 0.0), (112.0, -10.0), (115.5, -20.0)];
                        voltage_unit = :kV, value_unit = :MVAr,
                        vn_kV = 110.0, sbase_MVA = 100.0,
                        interpolation = :polynomial);
    qmin_MVAr = -50.0,
    qmax_MVAr = 50.0,
    sbase_MVA = 100.0,
)

pu = PUController(
    make_characteristic([(104.5, 20.0), (108.0, 14.0), (110.0, 10.0), (113.0, 6.0), (115.5, 0.0)];
                        voltage_unit = :kV, value_unit = :MW,
                        vn_kV = 110.0, sbase_MVA = 100.0,
                        interpolation = :polynomial);
    pmin_MW = 0.0,
    pmax_MW = 50.0,
    sbase_MVA = 100.0,
)

addProsumer!(
    net = net,
    busName = "B2",
    type = "SYNCHRONOUSMACHINE",
    p = 10.0,
    q = 0.0,
    qu_controller = qu,
    pu_controller = pu,
)

runpf!(net, 30, 1e-8, 0)
```

## Where to configure

| Use | Where |
|---|---|
| Call | `make_characteristic` plus `QUController`/`PUController` on `addProsumer!` |
| Case file (SCF) | `extra.<machine>.qu_control` / `pu_control`: points in per unit, interpolation mode, limits in MVAr / MW; `exportSCF` writes them for every machine with a controller, a hand-written file is validated when read ([SCF](scf.md)) |
| MATPOWER | `Pmin/Pmax/Qmin/Qmax` of a generator on a `PQ` bus (`BUS_TYPE = 1`) become constant P(U)/Q(U) controllers (a fixed value with limits, no curve), logged as an import message |
| CGMES, YAML | CGMES carries no Q(U) characteristic; the YAML configuration has no entry either, a characteristic is network data, not a setting |

Whatever the source, the Jacobian uses the analytic derivative of the
interpolated curve, not a difference quotient; outside the point range and
at a limit it is zero.

## Solver support and limitation

Supported on the default rectangular solver path; legacy polar/classic
solver modes are deprecated.

## Result printout semantics (`Type` vs `Control`)

- `Type`: structural bus class (`Slack`, `PV`, `PQ`).
- `Control`: voltage-dependent behavior (`Q(U)`, `P(U)`, `Q(U), P(U)`, or
  `-`).
- `Pg/Qg/Pl/Ql`: effective solved bus power components; for prosumers with
  `Q(U)`/`P(U)` these are the controller-evaluated setpoints at the solved
  bus voltage, not the static `p`/`q` inputs.

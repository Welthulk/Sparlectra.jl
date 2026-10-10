# Voltage-dependent prosumer control: Q(U) and P(U) (rectangular solver)

`QUController` sets `Q = f_Q(|V|)` and `PUController` sets `P = f_P(|V|)`
for a prosumer. Both are state-dependent injections, not PV equality
constraints: the bus type (`Slack`, `PV`, `PQ`) stays unchanged. `Q(U)` is
the droop for volt/var support, `P(U)` a curtailed or voltage-sensitive
active injection.

Capability limits are not controls. A CGMES `ReactiveCapabilityCurve` gives
the Q bounds as a function of P; the import evaluates it into ordinary
`minQ`/`maxQ` limits, which the
[Q-limit switching](q_limit_switching_strategy.md) enforces, never a
`QUController` (a bound as voltage-dependent setpoint would wander with
$|V|$). Both can coexist on one machine.

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

These derivatives use

```math
|V_i|_\varepsilon = \max(|V_i|, \varepsilon),\quad \varepsilon = 10^{-9}
```

against division by zero at a collapsed voltage.

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

Only the specified injection becomes state-dependent.

## Jacobian extension by chain rule (local terms)

`P_spec` / `Q_spec` depend on the local `|V_i|` only, so the extra Jacobian
terms sit in the row of bus `i` and the state columns of bus `i`.

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

The PV second row (`ΔV`) is the voltage magnitude mismatch and gets no
`Q(U)` term.

## Characteristic interpolation and derivative

A controller holds ordered points `(u_k, y_k)` in p.u. They are given in
p.u. (`voltage_unit=:pu`, `value_unit=:pu`) or in physical units
(`voltage_unit=:kV`, `value_unit=:MW`/`:MVAr`) with the conversion data
`vn_kV` and `sbase_MVA`.

`make_characteristic(...; interpolation = ...)` selects:

- `:linear` (default): piecewise linear.
- `:spline`: natural cubic spline through all points.
- `:polynomial`: one interpolating polynomial through all points.

With two points, `:spline` and `:polynomial` are the straight line. For
`:linear`, inside segment `[u_k, u_{k+1}]`:

```math
y(u) = y_k + \frac{y_{k+1}-y_k}{u_{k+1}-u_k}(u-u_k),
\qquad
\frac{dy}{du} = \frac{y_{k+1}-y_k}{u_{k+1}-u_k}.
```

For all modes:

- Below the first point and above the last: value clamped to that end
  point, derivative `0`.
- At an interior breakpoint: continuous value; `:linear` takes the slope of
  the segment left of the breakpoint, `:spline` and `:polynomial` have a
  continuous slope there.
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
| Case file (SCF) | `extra.<machine>.qu_control` / `pu_control`; `exportSCF` writes them for every machine with a controller ([SCF](scf.md)) |
| MATPOWER | `Pmin/Pmax/Qmin/Qmax` of a generator on a `PQ` bus (`BUS_TYPE = 1`) become constant P(U)/Q(U) controllers (a fixed value with limits, no curve), logged as an import message |
| CGMES, YAML | CGMES carries no Q(U) characteristic, and the YAML configuration has no entry: a characteristic is network data, not a setting |

Whatever the source, the Jacobian uses the analytic derivative of the
interpolated curve, not a difference quotient; outside the point range and
at a limit it is zero.

## Result printout semantics (`Type` vs `Control`)

- `Type`: structural bus class (`Slack`, `PV`, `PQ`).
- `Control`: voltage-dependent behavior (`Q(U)`, `P(U)`, `Q(U), P(U)`, or
  `-`).
- `Pg/Qg/Pl/Ql`: solved bus power components; for prosumers with
  `Q(U)`/`P(U)` these are the setpoints evaluated at the solved bus
  voltage, not the static `p`/`q` inputs.

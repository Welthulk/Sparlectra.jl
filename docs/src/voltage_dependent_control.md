# Voltage-dependent prosumer control: Q(U) and P(U) (rectangular solver)

This page documents voltage-dependent active/reactive power control for prosumers:

- `QUController`: `Q = f_Q(|V|)`
- `PUController`: `P = f_P(|V|)`

These are **soft controls** (state-dependent injections), not hard PV equality constraints.
The structural bus type (`Slack`, `PV`, `PQ`) remains unchanged.

## Physical interpretation

A controlled prosumer changes its setpoint according to local voltage magnitude.
Examples:

- `Q(U)` droop-like behavior for volt/var support.
- `P(U)` behavior for curtailed generation or voltage-sensitive active injection.

In Sparlectra this is attached to the prosumer (injection element), not to bus typing.

## Rectangular state and voltage magnitude

The rectangular state is

```math
V_i = e_i + j f_i, \qquad |V_i| = \sqrt{e_i^2 + f_i^2}.
```

For `|V_i| > 0`,

```math
\frac{\partial |V_i|}{\partial e_i} = \frac{e_i}{|V_i|},
\qquad
\frac{\partial |V_i|}{\partial f_i} = \frac{f_i}{|V_i|}.
```

Numerical safeguard: the implementation uses

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

If multiple prosumers share one bus, specified injections are aggregated by summation
(with sign convention: generation positive, load negative):

```math
P_i^{\mathrm{spec,tot}} = \sum_{k\in\mathcal{D}_i} P_{i,k}^{\mathrm{spec}},
\qquad
Q_i^{\mathrm{spec,tot}} = \sum_{k\in\mathcal{D}_i} Q_{i,k}^{\mathrm{spec}}.
```

## Mismatch equations in rectangular NR

The rectangular solver uses (for each non-slack bus):

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

So only the specified injection part becomes state-dependent.

## Jacobian extension by chain rule (local terms)

Because `P_spec` / `Q_spec` depend only on local `|V_i|`, extra Jacobian terms are local
(bus `i` row, bus `i` state columns only).

For active-power mismatch row:

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

For reactive-power mismatch row (PQ buses):

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

For PV second rows (`ΔV`) no `Q(U)` term is added because that row is voltage magnitude mismatch.

## Piecewise linear characteristic and derivative

Controllers use ordered points `(u_k, y_k)` internally in p.u.

You can now provide characteristic points either:

- directly in p.u. (`voltage_unit=:pu`, `value_unit=:pu`), or
- in physical units (`voltage_unit=:kV`, `value_unit=:MW`/`:MVAr`) with
  conversion metadata (`vn_kV`, `sbase_MVA`).

Inside segment `[u_k, u_{k+1}]`:

```math
y(u) = y_k + \frac{y_{k+1}-y_k}{u_{k+1}-u_k}(u-u_k),
\qquad
\frac{dy}{du} = \frac{y_{k+1}-y_k}{u_{k+1}-u_k}.
```

Behavior conventions:

- Below first point: clamp to first value, derivative `0`.
- Above last point: clamp to last value, derivative `0`.
- At breakpoints: segment-wise evaluation (continuous value; slope is side-dependent).
- If explicit min/max saturation is hit, output is clamped and derivative forced to `0`.

## API usage

```julia
using Sparlectra

qu = QUController(
    make_characteristic([(104.5, 30.0), (110.0, 0.0), (115.5, -20.0)];
                        voltage_unit = :kV, value_unit = :MVAr,
                        vn_kV = 110.0, sbase_MVA = 100.0);
    qmin_MVAr = -50.0,
    qmax_MVAr = 50.0,
    sbase_MVA = 100.0,
)

pu = PUController(
    make_characteristic([(104.5, 20.0), (110.0, 10.0), (115.5, 0.0)];
                        voltage_unit = :kV, value_unit = :MW,
                        vn_kV = 110.0, sbase_MVA = 100.0);
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

runpf!(net, 30, 1e-8, 0; method = :rectangular)
```

Notes:

- Existing p.u.-based calls remain fully supported.
- Controller evaluation and Jacobian assembly still operate in p.u. internally.
- MATPOWER import: if a generator is connected to a `PQ` bus (`BUS_TYPE = 1`),
  Sparlectra maps `Pmin/Pmax/Qmin/Qmax` to constant `P(U)`/`Q(U)` controllers
  (fixed characteristic with limits) and emits an import log message.

## Solver support and limitation

- Supported: `method = :rectangular`.
- Unsupported: legacy/deprecated `:polar_full` and `:classic` (clear runtime error if used with these controllers).

## Result printout semantics (`Type` vs `Control`)

Detailed result output has separate columns:

- `Type`: structural bus class (`Slack`, `PV`, `PQ`), unchanged semantics.
- `Control`: additional voltage-dependent behavior (`Q(U)`, `P(U)`, `Q(U), P(U)`, or `-`).
- `Pg/Qg/Pl/Ql`: effective solved bus power components; for prosumers with
  `Q(U)`/`P(U)` this reflects the controller-evaluated setpoints at the solved
  bus voltage (not only initial static `p`/`q` inputs).

This avoids overloading bus type strings and keeps existing workflows stable.

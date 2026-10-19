# Component Types

## Network Components

### Net

The container of a power system network.

```
Sparlectra.Net
```

### Node

A bus of the power system.

```
Sparlectra.Node
```

### Branch Components

```
Sparlectra.Branch
Sparlectra.BranchFlow
Sparlectra.BranchModel
```

A `Branch` carries the per-terminal flags `from_status`/`to_status`
(1 = closed, 0 = open) next to the aggregate `status` (`1` iff both
terminals are closed); `setBranchStatus!` sets all three,
`setBranchTerminalStatus!(br; from =, to =)` toggles single terminals.
Branches open at one terminal: "One-sided open branches" in the
[branch model](branchmodel.md).

### Prosumer Components

```
Sparlectra.ProSumer
```

### Shunt Components

```
Sparlectra.Shunt
```

### Transformer Components

```
Sparlectra.PowerTransformer
Sparlectra.PowerTransformerWinding
Sparlectra.PowerTransformerTaps
```

### Line Components

```
Sparlectra.ACLineSegment
```

Three builders add a line, one per kind of input data:

| Builder | Input | Stored on the line |
|---|---|---|
| [`addACLine!`](@ref) | totals of the line: `r_ohm`, `x_ohm`, `b_s`, `g_s`, optional `length_km` | ohm and siemens |
| [`addACLinePerKm!`](@ref) | `length_km`, `r_ohm_per_km`, `x_ohm_per_km`, `c_nf_per_km`, and `g_us_per_km` or `tanδ` | ohm and siemens (values times length) |
| [`addPIModelACLine!`](@ref) | per unit: `r_pu`, `x_pu`, `b_pu`, `g_pu` | per unit |

The first two convert to per unit with [`toPU_RXBG`](@ref) on the bus
voltage and `net.baseMVA`. `length_km` of `addACLine!` is information
only (reports, CGMES `Conductor.length`) and does not scale the values.
`addACLinePerKm!` derives `b` from the capacitance with the system
frequency of the network (`net.nominal_frequency_hz`, set from
`model.nominal_frequency_hz`, default 50 Hz) and `g` from `g_us_per_km` or
`tanδ`. All importers build per-unit PI-model lines.

### [Migrating from 0.42](@id migrating_from_0_42)

Up to 0.42 `addACLine!` took per-unit values `r`, `x`, `b` and a `length`.
For `length > 1` the values were the totals of the line, for `length <= 1`
(also 1.0) they were multiplied by `length`. The keywords `length`, `r`,
`x`, `b`, `g`, `c_nf_per_km` and `tanδ` now raise an `ArgumentError` that
names the replacement.

Lines built with `c_nf_per_km` had their charging susceptance read as per
unit, Z_base times too small; results of such networks change. If your own
scripts build lines with `c_nf_per_km`, recompute them with the new form
and compare the results.

In the table `Zb = Vn_kV^2 / baseMVA` (121 ohm at 110 kV and 100 MVA) and
`L` is the old `length`.

| Old form (0.42) | Old behaviour | New form |
|---|---|---|
| `addACLine!(length = L, r, x)` with `L > 1` | `r`, `x` per-unit totals, `b = g = 0` | same numbers: `addPIModelACLine!(r_pu = r, x_pu = x, b_pu = 0.0, status = 1)`; in physical units: `addACLine!(r_ohm = r*Zb, x_ohm = x*Zb, length_km = L)` |
| same with `L <= 1` | `r`, `x`, `b` multiplied by `L` | `addPIModelACLine!(r_pu = r*L, x_pu = x*L, b_pu = b*L, status = 1)` |
| with `b` (per unit) | `b` taken, `c_nf_per_km` and `tanδ` ignored | `b_pu = b` (`b*L` for `L <= 1`); in physical units `b_s = b/Zb` |
| `c_nf_per_km = c, tanδ = t`, no `b` | siemens read as per unit: `b = 2π·50·c·L·1e-9`, `g = t·b` | correct: `addACLinePerKm!(length_km = L, r_ohm_per_km = r*Zb/L, x_ohm_per_km = x*Zb/L, c_nf_per_km = c, tanδ = t)` (`r*Zb`, `x*Zb` per km for `L <= 1`); the old numbers: `addPIModelACLine!(b_pu = 2π*50*c*L*1e-9, g_pu = t*b_pu, ...)` |
| `c_nf_per_km = c`, no `tanδ`, no `b` | `c` ignored, `b = 0`, warning since 0.41.0 | drop `c` (first row); `addACLinePerKm!` would give `b = 2π·f·C·L` |
| `c_nf_per_km = 0.0` | `b = g = 0` | first or second row |
| `ACLineSegment(...; length, paramsBasedOnLength, c_nf_per_km, tanδ)` | constructor with the length rule | `ACLineSegment(; vn_kv, from, to, r, x, b = 0.0, g = 0.0, length_km = nothing, ratedS = nothing, isPIModel)` |

Further changes:
- `ratedS`, `status`, `from_status` and `to_status` keep their names.
  `addPIModelACLine!` requires `status` (the old `addACLine!` default was 1).
- `paramsBasedOnLength` was never a keyword of `addACLine!`; the field and
  the constructor keyword are gone.
- Line fields: `line.length` is now `line.length_km` (`nothing` for PI-model
  lines, before 1.0); `c_nf_per_km`, `tanδ` and `paramsBasedOnLength` are
  gone; `b` and `g` are never `nothing`.
- [`getLineRXBG`](@ref) returns ohm and siemens for a line from
  `addACLine!` or `addACLinePerKm!` (before: per unit with the length rule);
  per unit comes from `getLineRXBG_pu`.
- A missing `r_ohm`, `x_ohm` or per-km value raises `UndefKeywordError`;
  `addACLinePerKm!` refuses `length_km <= 0` and `g_us_per_km` together with
  `tanδ`.
- The CGMES export writes `Conductor.length` from `length_km`, else 1.0.

The fixed 50 Hz of the capacitance path is replaced by the system
frequency of the network, `net.nominal_frequency_hz`: a new `Net` takes
`model.nominal_frequency_hz` (default 50.0, must be positive), an SCF file
sets it from `meta.f_nom`. `addACLinePerKm!` reads the field when it is
called, and the SCF writer states it. See
[System frequency](@ref system_frequency).

### Link Components

```
Sparlectra.BusLink
```

## Basic Components

```
Sparlectra.Component
Sparlectra.ImpPGMComp
Sparlectra.ImpPGMComp3WT
```

## Enumerations

```
Sparlectra.ComponentTyp
Sparlectra.TrafoTyp
Sparlectra.NodeType
Sparlectra.ProSumptionType
```

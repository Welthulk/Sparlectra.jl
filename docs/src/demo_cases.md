# Shipped Demo Cases

Sparlectra ships four self-built networks in the [Sparlectra Case
Format](scf.md), tracked in `data/scf`, so a fresh install can run power
flow, state estimation, short circuit, N-1 and scenarios without
downloading anything.

!!! note "Provenance"
    These networks were built for Sparlectra. They are not the IEEE cases
    of the same bus count, and no data was taken from MATPOWER, DTF or any
    other source. Every number was chosen at construction or derived from
    a solved run of that construction, with fixed seeds.

The number in the name is the bus count; the `sp_` prefix marks them as
Sparlectra's own. Components carry invented place names: busbars as
`<Ort>_<kV>` (`Ostheim_110`, `Moorau_20`; a second busbar section gets a
`b` suffix), lines by their endpoints (`Ostheim_Neuwiese_1`/`_2` for the
double circuit; sp_case60's deliberate same-name pair carries the bare
corridor name on both legs), transformers by station plus both voltage
levels (`Moorau_110_20`). The short technical handle (`H1`, `A12`, `MC9`)
stays in the file as `external_id` ([case format](scf.md)). Each case
ships with measurements (a noisy CSV with the seed in the header, plus a
bad-data variant with one gross reading), a scenarios block (explicit
mode; N-1 expansions are chosen on the Runs page), short-circuit study
buses and a solved start state.

| Case | Buses | What it demonstrates |
|---|---|---|
| `sp_case5` | 5 | The hand-checkable entry: two voltage levels, an external grid as slack and IEC 60909 feeder, one PV machine, ratings chosen so that line outages overload the surviving parallel path (`line_Borntal_Lindhof_out`: 122 percent) while the base keeps 32 percent reserve. |
| `sp_case14` | 14 | The regulated feeder: a meshed 110 kV level, a radial 20 kV feeder behind a transformer with a REAL discrete tap controller (deadband 0.01 pu, settles in 5 outer iterations), a fixed-tap transformer off neutral, and N-1 outages that island the load-only feeder. |
| `sp_case60` | 60 | The feature-dense operated grid: 220/110/20 kV, a PST holding its corridor flow, a combined ratio-and-phase Schraegregler, an SVC, a 40 percent series-compensated corridor, a same-name double circuit across two busbar sections joined by an impedance-less link, and a bridge branch whose outage island solves via PV promotion. No branch above 80 percent in the base case: the screening fixture. |
| `sp_case188` | 188 | The benchmark: a 380 kV double ring with a double busbar, a steerable HVDC pair (60 MW setpoint), four parametric 110/20 kV zones, a STATCOM with a remote target, an SSSC, a PST and a series-compensated corridor. Import plus power flow in roughly 20 ms. |

One more shipped case is a MATPOWER file: `data/mpower/sp_case118.m`, a
synthetic 118-bus case with the cardinalities of the IEEE 118-bus system
(118 buses, 186 branches of which 9 transformers, 54 generators with a
good third of synchronous condensers, 99 loads, three areas, an 11-bus
345 kV backbone over a 138 kV grid) and its own topology, parameters and
names, produced and validated by a generator script outside the package.
It stands in for the downloaded IEEE case in the APSLF and scenarios
workshops: several machines sit at a reactive limit, every branch keeps
25 percent headroom, and a handful of N-1 outages island a bus or do not
converge.

In the Web UI the shipped cases appear in the case chooser; picking one
stages it into the case cache with its sidecars (per-case configuration
file, measurement CSVs). Programmatically, load them like any other SCF
file:

```julia
using Sparlectra
path = joinpath(pkgdir(Sparlectra), "data", "scf", "sp_case14.scf.json")
res = Sparlectra.import_case(path, load_sparlectra_config())
runpf!(res.net)
```

Reference values (bus voltages, iterations, losses, N-1 summary counts,
state-estimation objective and degrees of freedom, short-circuit currents
at two buses) live in `test/fixtures/demo_cases/` and
`data/scf/README.md`; tests use these cases, only tests of the download
and cache machinery fetch cases.

## PowSyBl bundles

`data/powsybl/` ships three PowSyBl networks (`ieee14.powsybl`,
`four_substations.powsybl`, `micro_grid_be.powsybl`: the IEEE 14-bus
case, a node-breaker network with two synchronous components, an HVDC
link and a PST, and the CGMES MicroGrid BE with a three-winding
transformer and dangling lines). Each directory holds the `.xiidm` file,
which Sparlectra reads directly (`run_sparlectra(casefile =
"data/powsybl/ieee14.powsybl/ieee14.xiidm")`), and the table bundle
pypowsybl wrote from it with the OpenLoadFlow reference solution in
`reference_buses.csv`. The Web UI case selector offers the
bundles; a `.xiidm` file of your own is imported the same way, no Python
is needed for either. The import reproduces the OpenLoadFlow voltages
within the bands of the test suite. See [PowSyBl Import](powsybl_import.md).

## CGMES deliveries of the demo cases

`data/cgmes_demo/<case>/` holds CGMES 2.4.15 deliveries (EQ, TP, SSH, SV)
that `tools/gen_cgmes_fixtures.jl` exports from `sp_case14`, `sp_case118`
and `sp_casePST` with Sparlectra's own exporter; a regeneration is byte
identical. They appear in the Web UI case selector as `<case>_cgmes.zip`.
The exporter
writes bus-branch deliveries only: no node-breaker topology, boundary sets
or DifferenceModel.

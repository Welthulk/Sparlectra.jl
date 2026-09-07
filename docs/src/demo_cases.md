# Shipped Demo Cases

Sparlectra ships four self-built networks in the [Sparlectra Case
Format](scf.md), tracked in `data/scf`, so a fresh install can run power
flow, state estimation, short circuit, N-1 and scenarios without
downloading anything.

!!! note "Provenance"
    These networks were built for Sparlectra. They are not the IEEE cases
    of the same bus count, and no data was taken from MATPOWER, DTF or any
    other source; public test networks served only as orientation for what
    a plausible network contains. Every number was either chosen when the
    case was constructed or derived from a solved run of that
    construction, and the shipped files are the result of that
    construction with fixed seeds.

The number in the name is the bus count; the `sp_` prefix marks them as
Sparlectra's own. Components carry invented PLACE names the way a real
control room labels them: busbars as `<Ort>_<kV>` (`Ostheim_110`,
`Moorau_20`; a second busbar section gets a `b` suffix), lines by
their endpoints (`Ostheim_Neuwiese_1`/`_2` for the double circuit;
sp_case60's deliberate same-name pair carries the bare corridor name on
both legs), transformers by their station plus both voltage levels
(`Moorau_110_20`). The names double as documentation, and the short
technical handle (`H1`, `A12`, `MC9`) stays in the file as
`external_id`, so the shipped cases also demonstrate the
name-vs-external_id channel of the [case format](scf.md). Each case
ships with its measurements (a noisy CSV with
the seed in the header, plus a bad-data variant with one deliberately
gross reading for the diagnostics demo), a scenarios block (explicit
scenarios; the block ships in explicit mode so an editor shows them, and the N-1 expansions are chosen on the Runs page), short-circuit study buses, and a
solved start state. The fast test profile runs all five run kinds against
tracked reference fixtures, so the cases are regression fixtures, not
decoration.

| Case | Buses | What it demonstrates |
|---|---|---|
| `sp_case5` | 5 | The hand-checkable entry: two voltage levels, an external grid as slack and IEC 60909 feeder, one PV machine, ratings chosen so that line outages overload the surviving parallel path (`line_Borntal_Lindhof_out`: 122 percent) while the base keeps 32 percent reserve. |
| `sp_case14` | 14 | The regulated feeder: a meshed 110 kV level, a radial 20 kV feeder behind a transformer with a REAL discrete tap controller (deadband 0.01 pu, settles in 5 outer iterations), a fixed-tap transformer off neutral, and N-1 outages that island the load-only feeder. |
| `sp_case60` | 60 | The feature-dense operated grid: 220/110/20 kV, a PST holding its corridor flow, a combined ratio-and-phase Schraegregler, an SVC, a 40 percent series-compensated corridor, a same-name double circuit across two busbar sections joined by an impedance-less link, and a bridge branch whose outage island solves via PV promotion. No branch above 80 percent in the base case: the screening fixture. |
| `sp_case188` | 188 | The benchmark: a 380 kV double ring with a double busbar, a steerable HVDC pair (60 MW setpoint), four parametric 110/20 kV zones, a STATCOM with a remote target, an SSSC, a PST and a series-compensated corridor. Import plus power flow in roughly 20 ms. |

In the local Web UI the shipped cases appear in the case chooser out of
the box; picking one stages it into the case cache together with its
sidecars (the per-case configuration file and the measurement CSVs), so
every run kind works on it immediately. Programmatically, load them like
any other SCF file; no special path exists:

```julia
using Sparlectra
path = joinpath(pkgdir(Sparlectra), "data", "scf", "sp_case14.scf.json")
res = Sparlectra.import_case(path, load_sparlectra_config())
runpf!(res.net)
```

The per-case reference values (bus voltages, iterations, losses, N-1
summary counts, state-estimation objective and degrees of freedom,
short-circuit currents at two buses) live in
`test/fixtures/demo_cases/` and in the per-case sections of
`data/scf/README.md`. Tests should prefer these shipped cases over
downloaded ones; only tests of the download and cache machinery itself
still fetch cases.

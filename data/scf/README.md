# Shipped SCF cases

These networks were built for Sparlectra. They are NOT the IEEE cases of
the same bus count, and no data was taken from MATPOWER, DTF or any other
source; MATPOWER and the DTF test networks served only as orientation for
what a plausible network contains (voltage levels, transformer and line
data, shunts, load and generation mix, ratings). Every number was either
chosen when the case was constructed or derived from a solved run of that
construction; the construction uses fixed seeds and reproduces the shipped
files bit for bit.

Since 2026-09-04 the cases carry invented PLACE names as reference names
(maintainer scheme): busbars are `<Ort>_<kV>` (`Ostheim_110`,
`Moorau_20`; a second busbar section gets a `b` suffix), lines are
named by their endpoints (`Ostheim_Neuwiese_1`/`_2` for a double circuit;
sp_case60's deliberate same-name pair carries the bare corridor name on
BOTH legs), and transformers by their station plus both voltage levels
(`Moorau_110_20`, role suffix `_PST` where needed). The short
technical handle (`H1`, `A12`, `MC9`) stays as `external_id`, so every
shipped file demonstrates the format's name-vs-external_id channel. The
two small cases use hand-picked names, sp_case60 and sp_case188 draw
from the deterministic prefix-root pool in the builder; all names are
deliberately compact (5 to 9 characters a place) and ASCII only, because
they travel through ids, CSVs and JSON.

Each `sp_case*` ships as four files: the case itself (`.scf.json`, with
roles, human component names, a solved `start_state`, a scenarios block,
short-circuit study buses, and its measurements), a noisy measurement CSV
(seed in the file header), a bad-data variant with one deliberately gross
reading for the diagnostics demo, and a tracked reference fixture under
`test/fixtures/demo_cases/`. The fast-profile group `demo_cases` runs all
five run kinds (power flow, N-1, scenarios, state estimation, short
circuit) against those fixtures, so the cases are regression fixtures,
not decoration. A fresh install can run every run kind on them without
downloading anything; new tests should use them instead of downloaded
cases (maintainer rule 2026-09-03).

## sp_case5 (5 buses, the hand-checkable entry)

Two voltage levels: a 110 kV four-bus double path plus one 20 kV bus
behind the transformer. External grid as slack and IEC 60909 feeder
(2000 MVA, reproduced exactly at the connection bus by the SC chain),
one 12 MW PV machine, 32 MW load. The base case keeps reserve (max
loading 68 percent); three of the four line outages overload a surviving
parallel path (`line_Borntal_Lindhof_out`: 122 percent), the
transformer outage islands the 20 kV load. Reference: 0.94 MW losses,
SE J 37.1 at dof 26.

## sp_case14 (14 buses, the regulated feeder)

Meshed 110 kV level (ring, two chords, the double circuit
`Ostheim_Neuwiese_1`/`_2`), a radial 20 kV feeder behind the REGULATED
transformer `Moorau_110_20` (a real discrete tap controller
holding `Moorau_20` at 1.0 pu, deadband 0.01 pu, settles
in 5 outer iterations), a fixed-tap transformer one step off neutral,
and a capacitive shunt. Three N-1 outages island the load-only feeder
(14 of 17 converge). The fixture carries the controller PARAMETERS next
to the result numbers. Reference: 1.68 MW losses, SE J 78.4 at dof 83.

## sp_case60 (60 buses, the feature-dense operated grid)

Three voltage levels: a 220 kV double ring feeds the meshed 110 kV zone
A through two coupling transformers and the PST `Nordtal_220_110_PST`
(flow controller), the combined ratio-and-phase Schraegregler
`Steinhof_110_20` feeds the 20 kV group MA, an SVC regulates the
`Bornloh_110` pocket (handle A9), one corridor is 40 percent
series-compensated, the SAME-NAME double circuit `Neubach_Waldsee`
spans the two busbar sections of station Neubach (`Neubach_110` and
`Neubach_110b`, impedance-less coupler), and zone B
hangs on a single bridge branch with its own 40 MW unit (the
bridge-outage island solves via PV promotion). Base state keeps reserve
(no branch above 80 percent, worst converged N-1 loading 72 percent):
the intended screening fixture. 74 N-1 cases, 51 converge, all 23
failures are structural islands. Reference: 11.1 MW losses, SE J 158.8
at dof 209.

## sp_case188 (188 buses, the benchmark)

The full family at benchmark size: a 380 kV double ring (double busbar
`Falkwik_380`/`Falkwik_380b`, handles V4/V4B, on an impedance-less
coupler), a steerable HVDC pair bridging Birkloh-Holtau (handles
V2-V6) at a 60 MW setpoint, four parametric 110/20 kV zones, a STATCOM
(remote target), an SSSC, a PST, and a series-compensated corridor.
Import plus power flow in about 20 ms. 217 N-1 cases, 173 converge, all
44 failures are structural islands. Ships an injection-based measurement
set (no flow rows, dof 189) to hold the size budget; all four bundles
together stay around 0.72 MB, inside the 0.9 MB budget. Reference:
23.8 MW losses, SE J 201.8.

## sp_casePST.scf.json (format fixture, not a demo case)

Renamed from warmup_casePST.scf.json on 2026-09-03 (same bytes, maintainer request). The **permanent version 1 fixture** of the Sparlectra Case Format
(issue #342). Every future reader must keep loading this file unchanged;
it is the regression guard against accidental breaking changes to the
FORMAT (the `sp_case*` fixtures above guard BEHAVIOR instead). It is
generated from `data/mpower/warmup_casePST.m` plus its measurement set
and carries the tap cascade (a ratio changer and an additional-voltage
phase shifter), a busbar coupler, shunts, state-estimation measurements,
and a start state. Regenerate it only together with a documented format
version bump.

`pgm_interop.json` is the plain power-grid-model interoperability probe
used by the adapter detection tests.

# HVDC Back-to-Back and the Pairing Controller

How Sparlectra models HVDC links, why there is no angle coupling, and how
the pairing controller (#297 Draft B) makes a link steerable without
changing the solver.

## What a back-to-back link is

A back-to-back (B2B) HVDC station couples two AC areas on one site: two
converter stations, a common DC circuit of negligible length, no long DC
line. The areas may run asynchronously (different frequency regimes or
uncoordinated phase references), and the power exchanged between them is a
**control setpoint** of the converters, not the result of a voltage-angle
difference. A point-to-point link is the same arrangement with a real DC
line between the stations; for the power flow both look identical.

## Why there is no angle coupling

An AC tie line couples the voltage angles of its terminals: the flow
follows the angle difference, and both areas share one angle reference. An
HVDC link does not. The converters decouple the two AC sides completely,
so:

- two areas joined only through HVDC remain **separate electrical
  islands**, each with its own reference,
- the transfer does not react to angle changes; it is whatever the
  converter control commands,
- removing the link changes the area balances, not the island structure.

Sparlectra models exactly this: the DC circuit is never mapped into the
Y-bus, and island detection treats HVDC terminals as injections.

## The model ladder

**Stage 0 (default): fixed injections.** Each converter becomes a fixed
PQ injection at its AC connection point, carrying the snapshot operating
point (CGMES SSH values, MATPOWER `dcline` `PF`/`PT` columns). This is the
industry-standard load-flow treatment (MATPOWER `toggle_dcline` builds the
same two bounded dummy generators) and reproduces the delivered snapshot
exactly. The link is not steerable; the two injections know nothing about
each other.

**Paired control (opt-in): the steerable link.** The same two injections,
plus one `HvdcPairControl` outer-loop controller that knows they belong to
one link. It enforces the **pairing invariant**

```math
P_\text{to} = P_\text{transfer} - P_\text{loss}, \qquad
P_\text{from} = -P_\text{transfer}
```

in injection convention (the from side exports `P_transfer` into the
link), with the loss model

```math
P_\text{loss} = \text{loss}_\text{mw} + \text{loss}_\text{fraction} \cdot |P_\text{transfer}|
```

which maps one-to-one onto MATPOWER's `LOSS0`/`LOSS1` pair; for CGMES the
loss is derived from the difference of the two SSH operating points. Each
terminal holds either a fixed reactive injection or a voltage target
(per-side secant on the terminal Q within its reactive range, the
machine-controller method). An optional transfer rating clamps
`|P_transfer|` with honest `at_limit`.

**Deliberately not modeled:** DC-side electrics and dynamics, converter
internals, multi-terminal DC grids, embedded (in-Newton) converter
equations. The outer loop is the mechanism, as for every other controller.

## How the controller works

The transfer side needs no iteration: injections are solver inputs, so the
setpoints are applied once and the invariant holds exactly after every
apply step. Only voltage-target terminals iterate, one secant per side
against the measured terminal voltage. Honesty rules match the other
controllers: `converged` only when the (possibly clamped) transfer equals
the target and every voltage target sits inside its deadband; `at_limit`
when the rating clamps the transfer or a voltage side is stuck at its
reactive bound. The reference injection of an island can never be part of
a setpoint pair (its power balances the island), and a PV-regulated
converter terminal stays P-only.

## Grid-forming mode (`mode = :island_feed`)

A setpoint and a slack exclude each other: the pairing needs the transfer
as a given, while a reference's power is the outcome of its island's
balance. The controller therefore refuses to pair a reference bus. The
converse is its own valid model, a **grid-forming (Vf) converter** feeding
an island that has no other source, such as an offshore platform or an
asynchronously supplied island grid.

`addHvdcPairControl!(net; from_bus, to_bus, mode = :island_feed, ...)`
models exactly that: the receiving converter is declared as the reference
of its island (an `EXTERNALNETWORKINJECTION` with `referencePri` at the
PCC), it holds voltage and angle there, and its output is whatever the
island draws. The dependency of the pairing inverts: instead of
`P_to = transfer - loss` with a given transfer, each outer iteration reads
the island balance and mirrors it onto the sending side,
`P_from = -(P_island + loss)`. `p_transfer_mw` must be omitted, the to
side carries neither `q_mvar` nor `vset_pu` (the slack holds its own
voltage), and the mirror counts as settled when applied and derived
transfer agree within `deadband_p_mw` (default `1e-3` MW).

`p_rating_mw` keeps its honest semantics: once the island draw exceeds the
rating, the sending side is pinned at the rating with `at_limit = true`
and `converged = false`. Note the model's limit here: the power flow's
reference always balances its island, so the deficit does not appear as a
voltage collapse in the solution; the flag is what marks the violated
rating. Grid-forming links are attached programmatically or via YAML
(`mode: island_feed`); the importers attach setpoint pairs only.

## Data sources

- **Programmatic:** `addHvdcPairControl!(net; from_bus, to_bus,
  p_transfer_mw, loss_mw, loss_fraction, ...)`.
- **YAML:** controller type `hvdc_pair` under `control.controllers`
  (see [Configuration](configuration.md)).
- **MATPOWER:** `matpower_import.matpower_dcline_mode = paired_control`
  attaches one controller per active `mpc.dcline` row, seeded with `PF`,
  `LOSS0`/`LOSS1`, and the terminal reactive values
  (see [MATPOWER Import Configuration](matpower_import.md)).
- **CGMES:** `cgmes_import.hvdc_mode = paired_control`. The importer
  groups the converters through the DC topology classes
  (`ACDCConverterDCTerminal`, `DCNode`, `DCLineSegment`); a component with
  exactly two converters is a link, back-to-back when no line segment
  participates. Detection runs in every mode and names the pairs in the
  import messages (see [CGMES Import](cgmes_import.md)). Validated on the
  conformity FullGrid set (two links) and the ReliCapGrid combined model,
  whose real border crossing yields transfer 109.118 MW with 9.098 MW loss
  from the two SSH operating points.

## Choosing the mode

Stage 0 is right whenever the task is reproducing a snapshot: it is exact,
free of extra machinery, and every existing result stays unchanged. Paired
control answers what-if questions: change `p_transfer_mw` and resolve to
see how both areas redispatch, let a terminal hold its voltage, or study
transfer limits with an honest `at_limit`. In the Web UI both modes are
one select away (PowerFlow form: "HVDC converters" for CGMES cases,
"DC-line mode" for MATPOWER cases); the controller shows up in the
controller summary, the element table, and the trace like every other
outer-loop device.

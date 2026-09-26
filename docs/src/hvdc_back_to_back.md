# HVDC Back-to-Back and the Pairing Controller

An HVDC link enters the power flow as two converter injections without
angle coupling. The optional pairing controller makes the link steerable
without changing the solver.

## Why there is no angle coupling

A back-to-back (B2B) station couples two AC areas on one site through two
converters and a short common DC circuit; a point-to-point link adds a
real DC line, and for the power flow both are identical. Unlike an AC tie
line, the converters decouple the two AC sides completely: the exchanged
power is a converter setpoint, not the result of an angle difference, two
areas joined only through HVDC remain separate electrical islands with
their own references, and removing the link changes the area balances,
not the island structure. The DC circuit is never mapped into the Y-bus;
island detection treats HVDC terminals as injections.

## The model ladder

**Stage 0 (default): fixed injections.** Each converter is a fixed PQ
injection at its AC connection point with the snapshot operating point
(CGMES SSH values, MATPOWER `dcline` `PF`/`PT` columns; MATPOWER's
`toggle_dcline` does the same). The snapshot is reproduced exactly; the
link is not steerable.

**Paired control (opt-in): the steerable link.** The same two injections
plus one `HvdcPairControl` outer-loop controller enforcing the pairing
invariant

```math
P_\text{to} = P_\text{transfer} - P_\text{loss}, \qquad
P_\text{from} = -P_\text{transfer}
```

in injection convention (the from side exports `P_transfer` into the
link), with the loss model

```math
P_\text{loss} = \text{loss}_\text{mw} + \text{loss}_\text{fraction} \cdot |P_\text{transfer}|
```

which maps onto MATPOWER's `LOSS0`/`LOSS1`; for CGMES the loss is the
difference of the two SSH operating points. Each terminal holds either a
fixed reactive injection or a voltage target (per-side secant on the
terminal Q within its reactive range). An optional transfer rating clamps
`|P_transfer|` and reports `at_limit`.

**Not modeled:** DC-side electrics and dynamics, converter internals,
multi-terminal DC grids, embedded (in-Newton) converter equations.

## How the controller works

Transfer setpoints are solver inputs and hold exactly after every apply
step; only voltage-target terminals iterate, one secant per side. `converged`
requires the (possibly clamped) transfer to equal the target and every
voltage target inside its deadband; `at_limit` means the rating clamps
the transfer or a voltage side is stuck at its reactive bound. The
reference injection of an island can never be part of a setpoint pair; a
PV-regulated converter terminal stays P-only.

## Grid-forming mode (`mode = :island_feed`)

A setpoint and a slack exclude each other, so the controller refuses to
pair a reference bus. The converse is its own model: a grid-forming (Vf)
converter feeding an island without any other source (offshore platform,
asynchronous island grid).

`addHvdcPairControl!(net; from_bus, to_bus, mode = :island_feed, ...)`
declares the receiving converter as the reference of its island (an
`EXTERNALNETWORKINJECTION` with `referencePri` at the PCC) whose output is
whatever the island draws. Each outer
iteration reads the island balance and mirrors it onto the sending side,
`P_from = -(P_island + loss)`. `p_transfer_mw` must be omitted, the to side
carries neither `q_mvar` nor `vset_pu`, and the mirror is settled when
applied and derived transfer agree within `deadband_p_mw` (default `1e-3`
MW).

Once the island draw exceeds `p_rating_mw`, the sending side is pinned at
the rating with `at_limit = true` and `converged = false`; the reference
still balances its island, so the flag, not a voltage collapse, marks the
violated rating. Grid-forming links are attached programmatically or via
YAML (`mode: island_feed`); the importers attach setpoint pairs only.

## Data sources

| Use | Where |
|---|---|
| Call | `addHvdcPairControl!(net; from_bus, to_bus, p_transfer_mw, loss_mw, loss_fraction, ...)`; `mode = :island_feed` for the grid-forming variant |
| YAML | controller type `hvdc_pair` under `control.controllers` ([Configuration](configuration.md)) |
| MATPOWER | `matpower_import.matpower_dcline_mode = paired_control`: one controller per active `mpc.dcline` row, seeded with `PF`, `LOSS0`/`LOSS1` and the terminal reactive values ([MATPOWER cases](matpower.md)) |
| CGMES | `cgmes_import.hvdc_mode = paired_control`: converters grouped through the DC topology classes (`ACDCConverterDCTerminal`, `DCNode`, `DCLineSegment`), a component with exactly two converters is a link, back-to-back when no line segment participates; detection runs in every mode and names the pairs in the import messages ([CGMES Import](cgmes_import.md)) |
| Result field | controller row: `converged`, `at_limit`, `invalid_topology` |
| Web UI | PowerFlow form: "HVDC converters" for CGMES cases, "DC-line mode" for MATPOWER cases |

## Meshed operation (AC tie in parallel to the link)

Two areas also tied by an AC branch form one synchronous island with
exactly one angle reference. Keeping both former references fails fast
with a message naming the island and its reference buses
(`AC island 1 has 2 angle references (A1, C1). ...`); the solver never
demotes a reference on its own.

After demoting one reference to PV (an `ExternalNetworkInjection` without
`referencePri`, or a voltage-regulated generator), the setpoint pair works
unchanged as a parallel PQ path: the link carries its ordered transfer,
the AC tie the rest.
`mode = :island_feed` is out of scope in a synchronous grid; the controller
reports `invalid_topology` once its reference bus was demoted.
Walkthrough: `examples/others/exp_hvdc_meshed_ac_tie.jl` and "Meshed
operation" in chapter 2 of the advanced workshop tour.

## Choosing the mode

Stage 0 reproduces a snapshot exactly. Paired control answers what-if
questions: change `p_transfer_mw` and resolve, let a terminal hold its
voltage, or study transfer limits with `at_limit`. In the Web UI both
modes are one select away.

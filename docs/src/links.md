# Links (Bus Couplers / Sectionalizers)

## Concept

A link is an impedance-less connection between two buses: a busbar coupler,
a sectionalizer, a node split or merge from a CIM import. Links come from
`addLink!`, from retained CGMES switches, or from the MATPOWER block
`mpc.sparlectra.links` (one `fbus tbus status` row per coupler, written
back on export; see [MATPOWER cases](matpower.md)).

A link has no impedance and no Y-bus entry. A closed link between bus $i$
and $j$ imposes the voltage constraint

```math
[
V_i = V_j
]
```

## Relation to KCL

Kirchhoff's current law holds per bus after topology processing; link flows
are reconstructed after the solve from the bus power balances.

## Zero-impedance loops (critical case)

Links that form a loop create a zero-impedance cycle:

```
Bus1 ──link── Bus2
  │             │
 link          link
  │             │
Bus3 ───────────
```

Without a voltage drop the current distribution is underdetermined;
treated electrically, the system is singular.

## Resolution via Pseudoinverse

With $A$ the incidence matrix of the link graph, $f$ the link flows and $b$
the nodal power imbalance,

```math
[
A f = b
]
```

is rank-deficient in a loop. The flows are the minimum-norm solution

```math
[
f = A^{+} b
]
```

with $A^{+}$ the Moore-Penrose pseudoinverse: a consistent KCL solution,
uniform flows in symmetric loops, no artificial circulating currents. Loop
flows are not unique; the minimum-norm solution is deterministic.

## Modeling guidelines

* Do not connect links to slack buses
* Prefer identical bus types for linked buses
* Use links only for topology, not impedance modeling
* Avoid large link-only subgraphs without measurements (state estimation)

## Example

```julia
linkNr = addLink!(net = net, fromBus = "Bus1", toBus = "Bus1a", status = 1)
```

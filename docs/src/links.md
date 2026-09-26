# Links (Bus Couplers / Sectionalizers)

## Concept

Links are impedance-less topological connections between buses: busbar
couplers, sectionalizers, node splitting or merging in CIM imports. They
enter a network via `addLink!`, from retained CGMES switches, or from the
MATPOWER extension block `mpc.sparlectra.links` (one `fbus tbus status`
row per coupler, written back on MATPOWER export; see
[MATPOWER cases](matpower.md)).

A link is not a physical branch: no impedance, no Y-bus entry, a voltage
equality constraint instead.

## Mathematical interpretation

A closed link between bus *i* and *j* imposes

```math
[
V_i = V_j
]
```

a topological constraint, not an admittance.

## Relation to KCL

Kirchhoff's Current Law is enforced per bus after topology processing,
not via admittance equations at the link; link flows are reconstructed
after solving via power balancing.

## Zero-impedance loops (critical case)

Multiple links forming a loop create a zero-impedance cycle:

```
Bus1 ──link── Bus2
  │             │
 link          link
  │             │
Bus3 ───────────
```

No voltage drop in the loop, an underdetermined current distribution, a
singular system if treated electrically.

## Resolution via Pseudoinverse

Link flows in such loops are computed as a minimum-norm solution. With
$A$ the incidence matrix of the link graph, $f$ the unknown link flows and
$b$ the nodal power imbalance,

```math
[
A f = b
]
```

is rank-deficient, so

```math
[
f = A^{+} b
]
```

with $A^{+}$ the Moore-Penrose pseudoinverse: a consistent KCL solution
with the minimum 2-norm flow distribution, uniform flows in symmetric
loops, no artificial circulation currents. Link flows in loops are not
unique; the minimum-energy solution is deterministic.

## Modeling guidelines

* Do not connect links to slack buses
* Prefer identical bus types for linked buses
* Use links only for topology, not impedance modeling
* Avoid large link-only subgraphs without measurements (SE context)

## Example

```julia
linkNr = addLink!(net = net, fromBus = "Bus1", toBus = "Bus1a", status = 1)
```

## Summary

| Aspect            | Behavior              |
| ----------------- | --------------------- |
| Electrical model  | none (no Y-bus entry) |
| Constraint        | voltage equality      |
| Loop handling     | pseudoinverse         |
| Flow uniqueness   | not unique            |
| Returned solution | minimum-norm          |


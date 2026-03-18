# Links (Bus Couplers / Sectionalizers)

## Concept

Links represent **impedance-less topological connections** between buses.

Typical use cases:

* Busbar couplers
* Sectionalizers
* Node splitting / merging in CIM imports

A link is **not a physical branch**:

* It has **no impedance**
* It is **not stamped into the Y-bus**
* It enforces **voltage equality constraints**

---

## Mathematical interpretation

A closed link between bus *i* and *j* imposes:

```math
[
V_i = V_j
]
```

This introduces a **topological constraint**, not an admittance.

---

## Relation to KCL

Since links are not part of the Y-bus, Kirchhoff’s Current Law (KCL) is not enforced via admittance equations at the link.

Instead:

* KCL is enforced **per bus after topology processing**
* Link flows are reconstructed **after solving** via power balancing

---

## Zero-impedance loops (critical case)

If multiple links form a loop, the system contains a **zero-impedance cycle**.

Example:

```
Bus1 ──link── Bus2
  │             │
 link          link
  │             │
Bus3 ───────────
```

This leads to:

* No voltage drop in the loop
* Underdetermined current distribution
* Singular system if treated electrically

---

## Resolution via Pseudoinverse

To compute link flows in such cases, Sparlectra uses a **minimum-norm solution**.

Let:

*  $A$  = incidence matrix of link graph
*  $f$  = unknown link flows
*  $b$  = nodal power imbalance

We solve:

```math
[
A f = b
]
```
Since the system is rank-deficient:

```math
[
f = A^{+} b
]
```

where $A^{+}$ is the **Moore–Penrose pseudoinverse**.

---

## Interpretation of the solution

The pseudoinverse yields:

* A **consistent KCL solution**
* The **minimum 2-norm flow distribution**
* Physically equivalent to:

  * uniform distribution of flows in symmetric loops
  * no artificial circulation currents

---

## Practical implications

* Link flows in loops are **not unique**
* Sparlectra returns the **minimum-energy solution**
* Results are stable and deterministic

---

## Modeling guidelines

* Do not connect links to slack buses
* Prefer identical bus types for linked buses
* Use links only for topology, not impedance modeling
* Avoid large link-only subgraphs without measurements (SE context)

---

## Example

```julia
linkNr = addLink!(net = net, fromBus = "Bus1", toBus = "Bus1a", status = 1)
```

---

## Summary

| Aspect            | Behavior              |
| ----------------- | --------------------- |
| Electrical model  | none (no Y-bus entry) |
| Constraint        | voltage equality      |
| Loop handling     | pseudoinverse         |
| Flow uniqueness   | not unique            |
| Returned solution | minimum-norm          |

---

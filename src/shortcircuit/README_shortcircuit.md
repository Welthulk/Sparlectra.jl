# Short circuit (developer notes)

## Purpose

Balanced 3-phase initial symmetrical short-circuit current calculation
(IEC 60909 style) plus the native short-circuit data container for networks
built without a CGMES delivery.

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `short_circuit.jl` | c-factor handling, source admittances, per-island matrices, the calculation | `runShortCircuit!`, `printShortCircuitResult` |
| `native_sc_data.jl` | Short-circuit source data for natively built networks | `NativeShortCircuitData` |

## Conventions

- Fault sweeps can run in parallel over chunks; per-chunk factorization
  copies are created serially BEFORE spawning (see the repository threading
  rules; a `local` in the task body is load-bearing and commented).
- Source data comes either from a CGMES delivery
  (`cgmesLineShortCircuitData`) or from the native container; the calculation
  does not care which.

## Reference

Rendered API documentation: `docs/src/reference_shortcircuit.md`. User-facing
documentation: `docs/src/short_circuit.md`.

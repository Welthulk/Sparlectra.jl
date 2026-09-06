# Contingency analysis (developer notes)

## Purpose

The N-1 contingency batch API: case generation (branches, generators, FOR001
outage lists), batch execution over a base network, weighting, and reporting.

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `contingency.jl` | Everything: case types, generation, batch run, metrics, CSV/report output | `runContingencies!`, `generateN1Branches`, `generateN1Generators`, `buildContingencyReport` |

## Conventions

- Batch execution can run contingency cases in parallel; the parallel site
  follows the repository threading rules (chunk indexing, serial fallback via
  `runtime.parallel.*`, per-call overrides). The serial path is the same
  function, not a copy.
- Contingency studies can also come from an SCF case file
  (`scf_case_studies`); the service run consumes them via
  `src/api/run_contingency_service.jl`.

## Reference

Rendered API documentation: `docs/src/reference_contingency.md`. User-facing
documentation: `docs/src/contingency.md`.

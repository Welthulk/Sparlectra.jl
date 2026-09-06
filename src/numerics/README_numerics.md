# Shared numerics (developer notes)

## Purpose

Numerical utilities shared by more than one solver path: condition-number
diagnostics for Newton-Raphson Jacobians and the Takahashi/Erisman-Tinney
selected inverse on a UMFPACK factorization.

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `condition_number.jl` | Condition estimation and verdict reporting for Jacobians | `condestJacobian`, `reportCondition` |
| `takahashi.jl` | Selected inverse (diagonal or pattern-restricted) from an existing LU | `takahashi_diag`, `takahashi_selected_inverse` |

## Conventions

- The Takahashi routines reuse an existing UMFPACK factorization; never copy
  a shared `UmfpackLU` under threads (see the repository threading rules).
- Condition diagnostics are opt-in and report-only; they never change solver
  decisions.

## Reference

Rendered API documentation: `docs/src/reference_numerics.md`.

# DC power flow (developer notes)

## Purpose

The linear DC approximation: B-matrix assembly, reduced linear solve,
per-island orchestration, status registry, and reporting. Also the DC angle
seed used by the rectangular start projection.

## Include order

```julia
include("powerflow_dc/dc_bmatrix.jl")
include("powerflow_dc/dc_solve.jl")
include("powerflow_dc/dc_status_workspace.jl")
include("powerflow_dc/dc_network_solver.jl")
include("powerflow_dc/dc_report.jl")
```

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `dc_bmatrix.jl` | MATPOWER-style B matrix assembly | `assemble_dc_bbus` |
| `dc_solve.jl` | Reduced angle solve for one island | `solve_dc_powerflow` |
| `dc_network_solver.jl` | Per-island orchestration and Net write-back | `_run_dc_powerflow!`, `_write_dc_solution!` |
| `dc_status_workspace.jl` | Weak-ref-keyed DC status registry | `_set_dc_pf_status!` |
| `dc_report.jl` | Report types, printing, public entry, DC angle seed | `rundcpf!`, `buildDcPowerFlowReport` |

## Conventions

- The DC status registry is deliberately separate from the rectangular status
  table; the two solvers must not share result metadata.
- `_dc_seed_rectangular_angles!` is a start-value helper for the rectangular
  solver, not part of the DC result contract.

## Reference

Rendered API documentation: `docs/src/reference_powerflow_dc.md`.

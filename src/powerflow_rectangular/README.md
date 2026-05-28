# Rectangular power-flow implementation

## Purpose

This directory contains helper layers for the rectangular complex-state Newton–Raphson power-flow path.

The main solver loop, public wrappers, and remaining Q-limit active-set orchestration still live in `src/jacobian_complex.jl`.

## Include order

The includes in `src/Sparlectra.jl` must stay in this dependency-aware order:

```julia
include("powerflow_rectangular/rectangular_core_equations.jl")
include("powerflow_rectangular/rectangular_jacobian_builders.jl")
include("powerflow_rectangular/rectangular_newton_step.jl")
include("powerflow_rectangular/rectangular_wrong_branch.jl")
include("powerflow_rectangular/rectangular_start_projection.jl")
include("powerflow_rectangular/rectangular_result_updates.jl")
include("powerflow_rectangular/rectangular_voltage_setpoints.jl")
include("jacobian_complex.jl")
```

Why this matters:

- `rectangular_jacobian_builders.jl` depends on core mismatch/equation helpers.
- `rectangular_newton_step.jl` depends on mismatch and Jacobian builders.
- start/diagnostic/result/setpoint helpers must be available before `jacobian_complex.jl` consumes them.

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `rectangular_core_equations.jl` | Rectangular mismatch and core equation helpers | `build_complex_jacobian`, `mismatch_rectangular`, `_max_rectangular_pv_voltage_residual` |
| `rectangular_jacobian_builders.jl` | Analytic rectangular Jacobian assembly (sparse and dense diagnostic variants) | `build_rectangular_jacobian_pq_pv_sparse`, `build_rectangular_jacobian_pq_pv_dense`, `build_rectangular_jacobian_pq_pv` |
| `rectangular_newton_step.jl` | Newton step update and autodamping/backtracking | `_validate_rectangular_damping`, `_apply_rectangular_delta`, `choose_rectangular_autodamp`, `complex_newton_step_rectangular` |
| `rectangular_start_projection.jl` | Start sanitization, optional DC-angle seed, and blend projection | `_sanitize_rectangular_start`, `_dc_angle_start_rectangular`, `project_rectangular_start` |
| `rectangular_wrong_branch.jl` | Post-solve wrong-branch diagnostics and status classification | `_wrap_to_180_deg`, `_circular_angle_spread_deg`, `_check_wrong_branch_solution` |
| `rectangular_result_updates.jl` | Final complex-voltage write-back into network node fields | `update_net_voltages_from_complex!` |
| `rectangular_voltage_setpoints.jl` | Bus voltage setpoint lookup and fallback preparation for rectangular setup | `_bus_voltage_setpoints_from_prosumers` |
| `../jacobian_complex.jl` | Remaining rectangular solver loop/public wrappers and Q-limit active-set workflow | `run_complex_nr_rectangular_for_net!`, `runpf_rectangular!`, `runpf!`-related rectangular entry points |

## Execution flow

High-level rectangular PF flow in the current split:

1. Build or access Y-bus and initial complex voltage state.
2. Build specified injections (`S`) for the current network state.
3. Resolve bus types and voltage setpoints.
4. Optionally project/sanitize the rectangular start.
5. Run Newton iterations with rectangular mismatch and Jacobian builders.
6. Optionally apply autodamping backtracking.
7. Handle Q-limit active-set logic in the remaining solver loop in `jacobian_complex.jl`.
8. Write final complex voltages back to node magnitude/angle fields.
9. Run wrong-branch diagnostics/status reporting where configured.

## Refactoring rules

- Keep comments and docstrings in English.
- Do not add `module` wrappers in this directory.
- Do not change public APIs from helper files.
- Preserve include-order dependencies.
- Do not mix mechanical extraction with behavior changes.
- Add/update this README when adding a new rectangular helper file.

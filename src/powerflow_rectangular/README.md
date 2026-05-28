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
include("powerflow_rectangular/rectangular_standalone_solver.jl")
include("powerflow_rectangular/rectangular_wrong_branch.jl")
include("powerflow_rectangular/rectangular_start_projection.jl")
include("powerflow_rectangular/rectangular_result_updates.jl")
include("powerflow_rectangular/rectangular_voltage_setpoints.jl")
include("powerflow_rectangular/rectangular_qlimit_trace.jl")
include("powerflow_rectangular/rectangular_qlimit_trace_logging.jl")
include("powerflow_rectangular/rectangular_qlimit_vset_adjustment.jl")
include("powerflow_rectangular/rectangular_qlimit_guard.jl")
include("jacobian_complex.jl")
```

Why this matters:

- `rectangular_jacobian_builders.jl` depends on core mismatch/equation helpers.
- `rectangular_newton_step.jl` depends on mismatch and Jacobian builders.
- `rectangular_standalone_solver.jl` depends on mismatch and Newton-step helpers.
- start/diagnostic/result/setpoint helpers must be available before `jacobian_complex.jl` consumes them.
- Q-limit trace bus-ID mapping helpers must be loaded before the remaining Q-limit workflow in `jacobian_complex.jl`.
- Q-limit trace logging helpers must be loaded after trace mapping helpers and before the remaining Q-limit workflow in `jacobian_complex.jl`.
- Q-limit `:adjust_vset` helper construction must be loaded before the remaining Q-limit workflow in `jacobian_complex.jl`.
- Q-limit guard preprocessing must be loaded before the remaining rectangular Q-limit active-set loop in `jacobian_complex.jl`.

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `rectangular_core_equations.jl` | Rectangular mismatch and core equation helpers | `build_complex_jacobian`, `mismatch_rectangular`, `_max_rectangular_pv_voltage_residual` |
| `rectangular_jacobian_builders.jl` | Analytic rectangular Jacobian assembly (sparse and dense diagnostic variants) | `build_rectangular_jacobian_pq_pv_sparse`, `build_rectangular_jacobian_pq_pv_dense`, `build_rectangular_jacobian_pq_pv` |
| `rectangular_newton_step.jl` | Newton step update and autodamping/backtracking | `_validate_rectangular_damping`, `_apply_rectangular_delta`, `choose_rectangular_autodamp`, `complex_newton_step_rectangular` |
| `rectangular_standalone_solver.jl` | Standalone rectangular array-level NR driver | `run_complex_nr_rectangular` |
| `rectangular_start_projection.jl` | Start sanitization, optional DC-angle seed, and blend projection | `_sanitize_rectangular_start`, `_dc_angle_start_rectangular`, `project_rectangular_start` |
| `rectangular_wrong_branch.jl` | Post-solve wrong-branch diagnostics and status classification | `_wrap_to_180_deg`, `_circular_angle_spread_deg`, `_check_wrong_branch_solution` |
| `rectangular_result_updates.jl` | Final complex-voltage write-back into network node fields | `update_net_voltages_from_complex!` |
| `rectangular_voltage_setpoints.jl` | Bus voltage setpoint lookup and fallback preparation for rectangular setup | `_bus_voltage_setpoints_from_prosumers` |
| `rectangular_qlimit_trace.jl` | Q-limit trace bus-ID mapping helpers used by rectangular diagnostics | `_qlimit_original_bus_id`, `_resolve_qlimit_trace_buses` |
| `rectangular_qlimit_trace_logging.jl` | Q-limit trace/logging diagnostics for rectangular active-set decisions | `_bus_has_online_voltage_regulator`, `_qlimit_violation`, `_print_rectangular_qlimit_trace` |
| `rectangular_qlimit_vset_adjustment.jl` | Q-limit `:adjust_vset` controller construction for rectangular workflow | `_build_vset_adjust_controllers` |
| `rectangular_qlimit_guard.jl` | Q-limit guard preprocessing before rectangular active-set iterations | `_apply_qlimit_guard_to_rectangular_active_set!` |
| `../jacobian_complex.jl` | Remaining rectangular solver loop/public wrappers and Q-limit active-set workflow | `run_complex_nr_rectangular_for_net!`, `runpf_rectangular!`, `runpf!`-related rectangular entry points |

## Execution flow

High-level rectangular PF flow in the current split:

1. Build or access Y-bus and initial complex voltage state.
2. Build specified injections (`S`) for the current network state.
3. Resolve bus types and voltage setpoints.
4. Optionally project/sanitize the rectangular start.
5. Run Newton iterations via the standalone array solver (`run_complex_nr_rectangular`) using rectangular mismatch and Newton-step helpers.
6. Optionally apply autodamping backtracking.
7. Apply Q-limit guard preprocessing, then handle active PV→PQ switching logic in the remaining solver loop in `jacobian_complex.jl`.
8. Write final complex voltages back to node magnitude/angle fields.
9. Run wrong-branch diagnostics/status reporting where configured.

## Refactoring rules

- Keep comments and docstrings in English.
- Do not add `module` wrappers in this directory.
- Do not change public APIs from helper files.
- Preserve include-order dependencies.
- Do not mix mechanical extraction with behavior changes.
- Add/update this README when adding a new rectangular helper file.

# ACPFlow framework (developer notes)

## Purpose

This directory contains the configuration-driven AC power-flow framework: the
public entry points (`run_sparlectra`, `run_sparlectra_cases`, `run_acpflow`),
case resolution and import context, solver dispatch, the automatic power-flow
mode, result postprocessing, and status classification. It orchestrates; the
numerical solvers live in `src/powerflow_rectangular/` and `src/powerflow_dc/`.

## Include order

`acpflow.jl` is the include hub; `src/Sparlectra.jl` includes it once. The
order inside is dependency-aware:

```julia
include("start_modes.jl")
include("net_cache.jl")
include("import_context.jl")
include("execution.jl")
include("auto_powerflow.jl")
include("apslf_execution.jl")
include("status.jl")
include("output.jl")
include("entrypoint.jl")
```

Four files are included separately by `src/Sparlectra.jl` because they depend
on later layers: `island_diagnostics.jl` (needs the result types),
`solver_core.jl`, `solver_interface.jl` (needs the DC solver), and
`apslf_solver.jl` (needs the solver interface).

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `acpflow.jl` | Include hub for the framework layer | (includes only) |
| `entrypoint.jl` | Public framework entry points | `run_sparlectra`, `run_sparlectra_cases`, `run_acpflow` |
| `import_context.jl` | Casefile resolution, shared MATPOWER/SCF import context, config copy helpers | `_resolve_sparlectra_casefile`, `_import_sparlectra_context`, `_copy_sparlectra_with_powerflow` |
| `start_modes.jl` | Apply configured MATPOWER start modes to the imported net | `_apply_matpower_start_modes!` |
| `net_cache.jl` | Opt-in binary cache of parsed case and built Net, keyed by file hash and import options | `_net_cache_key`, `_net_cache_load`, `_net_cache_store` |
| `execution.jl` | Central solver dispatch for one configured run | `_execute_sparlectra_powerflow!` |
| `auto_powerflow.jl` | Automatic power-flow mode: feature collection, strategy selection, escalation ladder | `collect_auto_pf_features`, `select_auto_pf_strategy`, `auto_pf_hints` |
| `apslf_execution.jl` | Solver wiring for `power_flow.solver = :apslf` | `_run_apslf_powerflow!` |
| `apslf_solver.jl` | Bridge to the external AnalyticLoadFlow solver via the PFModel interface | `ApslfSolver`, `solvePf` |
| `solver_core.jl` | Shared numerical helpers: currents, injections, sparse solves | `calc_currents`, `calc_injections`, `solve_sparse_system`, `solve_linear` |
| `solver_interface.jl` | Solver-agnostic PF model bridge for external solvers | `buildPfModel`, `solvePf`, `applyPfSolution!`, `runpf_external!` |
| `status.jl` | `SparlectraRunResult` and framework status classification | `_compose_framework_status`, `_build_sparlectra_result` |
| `output.jl` | Result postprocessing: losses, KCL link flows, printing | `_postprocess_sparlectra_result!` |
| `island_diagnostics.jl` | AC island detection and per-island mismatch diagnostics | `_ac_island_components`, `_status_for_island` |

## Execution flow

```text
run_sparlectra
  -> import_case (src/import/case_import.jl) for case input
  -> _execute_sparlectra_powerflow!  (solver dispatch)
     -> runpf_rectangular! / rundcpf! / _run_apslf_powerflow!
  -> _postprocess_sparlectra_result!
  -> _build_sparlectra_result  (status classification)
```

The automatic mode (`power_flow.mode = auto`) wraps the dispatch in a
read-only feature analysis plus an escalation ladder; it rewrites only the run
configuration handed to the attempt, never the loaded template.

## Reference

Rendered API documentation: `docs/src/reference_acpflow.md` (the APSLF bridge
has its own page, `docs/src/reference_apslf.md`).

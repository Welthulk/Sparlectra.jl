# Central Configuration

Sparlectra uses one central YAML configuration flow.

| Item | Path / mechanism | Notes |
|---|---|---|
| Default config | `src/configuration.yaml.example` | Shipped defaults, version-controlled. |
| User override | `examples/configuration.yaml` | Local override, gitignored. |
| Explicit config file | `load_sparlectra_config("/path/to/file.yaml")` | Replaces default user override path. |
| Environment selector | `SPARLECTRA_CONFIGURATION_YAML` | Used by helper path selection in example workflows. |

## Precedence and validation

| Priority (low → high) | Source |
|---|---|
| 1 | `src/configuration.yaml.example` |
| 2 | `examples/configuration.yaml` or explicit `user_path` |
| 3 | `cli_overrides` |
| 4 | `overrides` |

Unknown keys are rejected. Removed keys like `matpower_import.benchmark` are rejected with a migration hint.

## Configuration section overview

| Section | YAML path | Typed config | Purpose | Typical user |
|---|---|---|---|---|
| Power flow | `power_flow` | `PowerFlowConfig` | Solver method, tolerances, start values, Q limits | PF users |
| MATPOWER import | `matpower_import` | `MatpowerImportConfig` | Case selection and MATPOWER interpretation | MATPOWER users |
| State estimation | `state_estimation` | `StateEstimationConfig` | WLS solver and observability options | SE users |
| Runtime | `runtime` | `RuntimeConfig` | Julia/BLAS thread handling | Performance users |
| Output | `output` | `OutputConfig` | Console/logfile output detail | All users |
| Diagnostics | `diagnostics` | `DiagnosticsConfig` | Expensive/compact diagnostics | Advanced users |
| Performance | `performance` | `PerformanceConfig` | Profiling and timing reports | Performance users |
| Benchmark | `benchmark` | `BenchmarkConfig` | Repeated timing runs | Benchmark users |
| Extensions | `extensions` | reserved | Compatibility placeholder | Internal/advanced |

Reserved extension leaf: `extensions.reserved` (kept for forward compatibility, not an active runtime control).

## Migration notes

| Old key | New key | Action |
|---|---|---|
| `matpower_import.benchmark` | `benchmark.enabled` | Removed from `matpower_import`; use top-level benchmark section. |
| `methods` (top-level legacy) | `benchmark.methods` | Keep methods under top-level `benchmark`. |
| `opt_sparse` | `power_flow.sparse` | Prefer structured power-flow key. |
| `max_ite` | `power_flow.max_iter` | Prefer structured power-flow key. |

## Detailed references

- [Power-Flow Configuration](powerflow_configuration.md)
- [MATPOWER Import Configuration](matpower_import.md)
- [State-Estimation Configuration](state_estimation_configuration.md)
- [Performance and Profiling Configuration](performance_profiling.md)

## Key index (selected paths)

- `power_flow.rectangular_workspace_reuse`
- `power_flow.rectangular_preallocate_workspace`
- `power_flow.rectangular_workspace_min_buses`
- `performance.representative_warmup_runs`
- `performance.compare_cold_warm`
- `output.result_table_max_rows`
- `output.result_table_large_case_threshold_buses`
- `output.result_table_large_case_mode`

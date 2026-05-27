# Central Configuration

Sparlectra uses one **central, typed configuration model**.

The main entry point is `Sparlectra.load_sparlectra_config(...)`, which:

1. loads YAML defaults,
2. optionally overlays a user YAML,
3. applies programmatic overrides,
4. validates key names and value domains,
5. builds a typed `SparlectraConfig` object used by runtime modules.

## Configuration files and selectors

| Item | Path / mechanism | Role |
|---|---|---|
| Default config | `src/configuration.yaml.example` | Version-controlled baseline for all options. |
| User override | `examples/configuration.yaml` | Local override file (project default user path). |
| Explicit config path | `load_sparlectra_config("/path/to/file.yaml")` | Replaces default user override path for that load call. |
| Environment-based path selection helper | `SPARLECTRA_CONFIGURATION_YAML` via `configuration_path_from_inputs(...)` | Used by script/example path resolution workflows. |

## Merge precedence

Configuration precedence (low → high):

1. `src/configuration.yaml.example`
2. `examples/configuration.yaml` (or explicit `user_path`)
3. `cli_overrides`
4. `overrides`

Unknown keys are rejected during validation. Removed keys are also rejected with migration hints (for example `matpower_import.benchmark` → `benchmark.enabled`).

## Typed central object

The merged YAML is converted into:

- `SparlectraConfig`
  - `powerflow::PowerFlowConfig`
  - `state_estimation::StateEstimationConfig`
  - `matpower::MatpowerImportConfig`
  - `performance::PerformanceConfig`
  - `benchmark::BenchmarkConfig`
  - `runtime::RuntimeConfig`
  - `diagnostics::DiagnosticsConfig`
  - `output::OutputConfig`
  - `control::ControlConfig`

This typed model is the canonical internal representation that should be consumed by power-flow, MATPOWER import, state estimation, output/reporting, performance profiling, benchmark runners, and future modules.

## YAML structure (section map)

| YAML section | Typed section | Purpose | Status |
|---|---|---|---|
| `power_flow` | `PowerFlowConfig` | Rectangular power-flow solver controls, start mode, Q-limits | Public / supported |
| `matpower_import` | `MatpowerImportConfig` | MATPOWER case path + import interpretation options | Public / supported |
| `state_estimation` | `StateEstimationConfig` | State-estimation runtime controls | Public / supported |
| `output` | `OutputConfig` | Console/logfile behavior and result table sizing | Public / supported |
| `performance` | `PerformanceConfig` | Profiling/reporting toggles and diagnostic volume controls | Public / supported |
| `benchmark` | `BenchmarkConfig` | Repeated benchmark-run controls | Public / supported |
| `control` | `ControlConfig` | Generic controller outer-loop orchestration controls (`control.controllers` reserved for future YAML definitions) | Public / supported |
| `runtime` | `RuntimeConfig` | Julia/BLAS thread control knobs for entry workflows | Public / supported |
| `diagnostics` | `DiagnosticsConfig` | Effective-config logging and diagnostics render controls | Public / supported |
| `extensions` | reserved (not mapped to typed runtime fields) | Future extension placeholder | Reserved |

### Supported power-flow path

The supported public power-flow solver path is **rectangular** (`power_flow.method: rectangular`).

## Option lifecycle and compatibility policy

Sparlectra docs distinguish option categories:

- **Public / supported**: keys in `src/configuration.yaml.example` and typed section constructors.
- **Reserved**: schema placeholders such as `extensions.reserved` for forward compatibility.
- **Deprecated compatibility aliases**: accepted for migration but not preferred in new YAML (for example `max_ite` and some start-projection alias keys).
- **Removed**: explicitly rejected keys with migration error guidance (for example `matpower_import.benchmark`).
- **Internal-only implementation details**: not documented as stable external user API.

Prefer canonical nested keys shown in the example YAML and module pages.

## Loader and validation behavior

### Key validation

- User YAML keys and override keys are validated against the default-schema tree from `src/configuration.yaml.example`.
- Unknown keys throw `ArgumentError`.
- Type/domain checks are applied while constructing typed config objects (for example Symbol allow-lists and positivity checks).

### Alias handling

Some legacy aliases are accepted in constructors for migration convenience, but canonical docs and examples use structured section keys.

### Caching

`load_sparlectra_config(...)` caches the typed result when loading from unchanged files without overrides. File hash/mtime changes invalidate cache reuse automatically.

## Module consumption model

The central typed config is intended to be consumed by each module from its own section:

- **Power flow**: reads `config.powerflow` (solver controls, start mode, Q-limit controls).
- **MATPOWER import**: reads `config.matpower` and combines with PF/output/performance sections as needed by example/runner paths.
- **State estimation**: reads `config.state_estimation` for method/tolerances/observability toggles.
- **Output/reporting**: reads `config.output` and `config.diagnostics`.
- **Performance profiling**: reads `config.performance` and `config.benchmark`.
- **Future modules**: should add a dedicated typed section and canonical YAML subtree.

## Minimal example YAML

```yaml
power_flow:
  method: rectangular
  tol: 1.0e-5
  max_iter: 80
  autodamp: true
  autodamp_min: 0.05

  start_mode:
    angle_mode: dc
    voltage_mode: profile_blend
    profile_source: matpower_reference
    start_projection: true
    try_dc_start: true
    try_blend_scan: true
    blend_lambdas: [0.25, 0.5, 0.75]
    dc_angle_limit_deg: 60.0

  qlimits:
    enabled: true
    start_iter: 3
    start_mode: iteration_or_auto

matpower_import:
  case: case14.m
  auto_profile: recommend
  pv_voltage_source: gen_vg

state_estimation:
  enabled: true
  method: wls

output:
  console_summary: true
  logfile_results: full

performance:
  enabled: true
  level: iteration

benchmark:
  enabled: true
  methods: [rectangular]

runtime:
  julia_threads: keep
  blas_threads: keep

extensions:
  reserved: true
```

For complete key references and allowed-value tables, see the module-specific pages below.

## Control configuration (generic outer loop)

```yaml
control:
  enabled: true
  max_outer_iterations: 20
  trace: true
  log_iterations: true
  stop_on_pf_failure: true
  controllers: []
```

| Key | Type | Default | Meaning |
|---|---:|---:|---|
| `enabled` | Bool | `true` | Enables the generic outer-loop control framework. |
| `max_outer_iterations` | Int | `20` | Global outer-loop cap. Does not control inner NR iterations. |
| `trace` | Bool | `true` | Collect machine-readable control trace rows. |
| `log_iterations` | Bool | `true` | Enables optional per-iteration control logging hooks. |
| `stop_on_pf_failure` | Bool | `true` | Stops control orchestration when inner PF fails. |
| `controllers` | Vector | `[]` | Reserved for future YAML controller definitions; leave empty for current programmatic controller setup. |

In Stage 1, controllers are typically attached programmatically via
`addTapController!` / `addPowerTransformerControl!`.

### Demo controller YAML vs. central `control.controllers`

The tap-control demo may read `examples/tap_control_demo_grid.yaml` for
example setpoints and transformer tap/phase parameters (`oltc`, `pst`,
`schraeg`). This is an example-specific
input file consumed by `examples/tap_control_demo_grid.jl`.

It does not define the central `control.controllers` schema. Today,
`control.controllers` remains reserved/future and should be left as `[]` in
central configuration files.

## Migration notes

| Legacy / old key | Canonical key | Notes |
|---|---|---|
| `matpower_import.benchmark` | `benchmark.enabled` | Removed from `matpower_import`; now top-level benchmark section. |
| `methods` (top-level legacy path) | `benchmark.methods` | Keep benchmark methods in `benchmark`. |
| `max_ite` | `power_flow.max_iter` | Legacy alias; prefer canonical nested key. |

## Detailed references

- [Power-Flow Configuration](powerflow_configuration.md)
- [MATPOWER Import Configuration](matpower_import.md)
- [State-Estimation Configuration](state_estimation_configuration.md)
- [Performance and Profiling Configuration](performance_profiling.md)


## Complete default-key index

The following canonical keys are currently present in `src/configuration.yaml.example`:

- `benchmark`
- `benchmark.enabled`
- `benchmark.methods`
- `benchmark.samples`
- `benchmark.seconds`
- `benchmark.show_once`
- `benchmark.show_once_max_nodes`
- `benchmark.show_once_output`
- `control`
- `control.controllers`
- `control.enabled`
- `control.log_iterations`
- `control.max_outer_iterations`
- `control.stop_on_pf_failure`
- `control.trace`
- `diagnostics`
- `diagnostics.console_auto_profile`
- `diagnostics.console_diagnostics`
- `diagnostics.console_max_rows`
- `diagnostics.console_q_limit_events`
- `diagnostics.console_summary`
- `diagnostics.log_effective_config`
- `diagnostics.logfile_diagnostics`
- `extensions`
- `extensions.reserved`
- `matpower_import`
- `matpower_import.auto_profile`
- `matpower_import.auto_profile_log`
- `matpower_import.bus_shunt_model`
- `matpower_import.case`
- `matpower_import.compare_voltage_reference`
- `matpower_import.enable_pq_gen_controllers`
- `matpower_import.preallocate_min_buses`
- `matpower_import.preallocate_network`
- `matpower_import.pv_voltage_mismatch_tol_pu`
- `matpower_import.pv_voltage_source`
- `matpower_import.ratio`
- `matpower_import.shift_sign`
- `matpower_import.shift_unit`
- `output`
- `output.console_auto_profile`
- `output.console_diagnostics`
- `output.console_max_rows`
- `output.console_q_limit_events`
- `output.console_summary`
- `output.logfile_diagnostics`
- `output.logfile_performance`
- `output.logfile_results`
- `output.logfile_warnings`
- `output.result_table_large_case_mode`
- `output.result_table_large_case_threshold_buses`
- `output.result_table_max_rows`
- `performance`
- `performance.compact_logging`
- `performance.compare_cold_warm`
- `performance.enabled`
- `performance.level`
- `performance.max_diagnostic_rows`
- `performance.print_to_console`
- `performance.representative_warmup_runs`
- `performance.show_allocations`
- `performance.show_iteration_table`
- `performance.skip_branch_neighborhood_report`
- `performance.skip_expensive_diagnostics`
- `performance.skip_reference_comparison`
- `performance.write_to_logfile`
- `power_flow`
- `power_flow.autodamp`
- `power_flow.autodamp_min`
- `power_flow.wrong_branch_detection`
- `power_flow.wrong_branch_rescue`
- `power_flow.wrong_branch_min_vm_pu`
- `power_flow.wrong_branch_max_vm_pu`
- `power_flow.wrong_branch_max_angle_spread_deg`
- `power_flow.wrong_branch_max_branch_angle_deg`
- `power_flow.wrong_branch_min_low_vm_count`
- `power_flow.wrong_branch_rescue_max_attempts`
- `power_flow.flatstart`
- `power_flow.max_iter`
- `power_flow.method`
- `power_flow.qlimits`
- `power_flow.qlimits.auto_q_delta_pu`
- `power_flow.qlimits.cooldown_iters`
- `power_flow.qlimits.enabled`
- `power_flow.qlimits.guard`
- `power_flow.qlimits.guard.accept_bounded_violations`
- `power_flow.qlimits.guard.enabled`
- `power_flow.qlimits.guard.freeze_after_repeated_switching`
- `power_flow.qlimits.guard.log`
- `power_flow.qlimits.guard.max_remaining_violations`
- `power_flow.qlimits.guard.max_switches`
- `power_flow.qlimits.guard.min_q_range_pu`
- `power_flow.qlimits.guard.narrow_range_mode`
- `power_flow.qlimits.guard.violation_mode`
- `power_flow.qlimits.guard.violation_threshold_pu`
- `power_flow.qlimits.guard.zero_range_mode`
- `power_flow.qlimits.hysteresis_pu`
- `power_flow.qlimits.lock_pv_to_pq_buses`
- `power_flow.qlimits.start_iter`
- `power_flow.qlimits.start_mode`
- `power_flow.qlimits.trace_buses`
- `power_flow.rectangular_preallocate_workspace`
- `power_flow.rectangular_workspace_min_buses`
- `power_flow.rectangular_workspace_reuse`
- `power_flow.start_mode`
- `power_flow.start_mode.accept_unmeasured_dc_start`
- `power_flow.start_mode.angle_mode`
- `power_flow.start_mode.blend_lambdas`
- `power_flow.start_mode.branch_guard`
- `power_flow.start_mode.dc_angle_limit_deg`
- `power_flow.start_mode.measure_candidates`
- `power_flow.start_mode.profile_source`
- `power_flow.start_mode.reuse_import_data`
- `power_flow.start_mode.start_projection`
- `power_flow.start_mode.try_blend_scan`
- `power_flow.start_mode.try_dc_start`
- `power_flow.start_mode.voltage_mode`
- `power_flow.tol`
- `runtime`
- `runtime.blas_threads`
- `runtime.julia_threads`
- `runtime.print_thread_config`
- `state_estimation`
- `state_estimation.enabled`
- `state_estimation.flatstart`
- `state_estimation.jac_eps`
- `state_estimation.max_iter`
- `state_estimation.method`
- `state_estimation.observability`
- `state_estimation.observability.enabled`- `state_estimation.tol`
- `state_estimation.update_net`

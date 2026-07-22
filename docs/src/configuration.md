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
  - `transformer::TransformerConfig`
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
| `transformer` | `TransformerConfig` | Transformer-modeling options shared by all importers (tap-changer model) | Public / supported |
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

### MATPOWER import metadata and DC-line options

`matpower_import.apply_bus_names`, `apply_branch_names`, and
`apply_branch_kind` default to `false` so existing imports keep numeric bus
names and heuristic branch classification. Enable them when MATPOWER cases
carry validation metadata such as standard `mpc.bus_name` and user-defined
`mpc.branch_name`/`mpc.branch_kind`. `import_for001_contingencies` defaults to
`true` and preserves user-defined `mpc.for001_contingencies` for validation
workflows. `matpower_import.matpower_dcline_mode` defaults to
`:pf_injections`; use explicit `:reject_active` when strict active-row rejection is required. The default uses simple MATPOWER power-flow
DC-line terminal injections are desired. OPF and `dclinecost` remain
unsupported.

### AC island solving

Disconnected AC islands are not tied together by MATPOWER `mpc.dcline`
terminal injections: those injections affect bus power balance but do not add
Ybus branches. `power_flow.islands.enabled: true` is the default, so these disconnected AC components are solved independently unless the user explicitly disables island solving:

```yaml
power_flow:
  islands:
    enabled: true
    mode: solve_independent
    reference_policy: matpower_like
```

The `matpower_like` policy keeps an existing island REF/Slack bus, otherwise
promotes the deterministic first PV/voltage-controlled bus as that island's
angle reference. Islands without REF/Slack or PV support fail before NR. When
multiple islands are detected, Sparlectra writes `ac_islands.csv` in the run
directory with bus, branch, generator/load, DC-line terminal, power-balance,
reference, and status diagnostics.

### Transformer tap-changer model

`transformer.tap_changer_model` selects, for **all** transformers of an
imported case, whether the tap changer is treated as electrically ideal or as
affecting the series impedance:

```yaml
transformer:
  # Allowed values: ideal, impedance_correction
  tap_changer_model: ideal
```

- `ideal` (default): tap steps only change the complex winding ratio; the
  transformer series impedance (R, X) keeps its neutral-position value. This
  preserves prior Sparlectra behavior.
- `impedance_correction`: tap steps additionally re-refer the series impedance
  through the tapped winding. R and X are scaled with the squared magnitude of
  the regulating vector, `|1 + f·e^(jφ)|²`, where `f` is the longitudinal
  regulating-voltage fraction and `φ` the skew angle (0° for a pure
  longitudinal/ratio tap changer).

The option is read by both the MATPOWER importer
(`createNetFromMatPowerFile`/`createNetFromMatPowerCase`) and the native DTF
importer (`DTFImporter.build_net`/`createNetFromDTFFile`), but the correction
math itself lives in a single place,
[`calcTapCorrectedRX`](@ref)/[`calcTapImpedanceCorrectionFactor`](@ref) in
`src/equicircuit.jl`, so importers stay free of duplicated tap-impedance
mathematics.

The corrected R/X values are what a subsequent [`writeMatpowerCasefile`](@ref)
export writes out (both the MATPOWER and the native DTF import path use the
same `Branch.r_pu`/`Branch.x_pu` fields the exporter reads). When such a case
was corrected on import, the exporter records a
`mpc.sparlectra.tap_changer_model = 'impedance_correction'` roundtrip marker so
a Sparlectra reimport does not reapply the correction a second time; see [Tap-impedance
correction and reimport](matpower_import.md#tap-impedance-correction-and-reimport).

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

  start_current_iteration:
    enabled: false
    max_iter: 10
    tol: 1.0e-3
    damping: 0.5
    accept_only_if_improved: true
    min_improvement_factor: 0.98
    vm_min_pu: 0.5
    vm_max_pu: 1.5
    max_angle_step_deg: 30.0
    only_for_large_cases: false

  qlimits:
    enabled: true
    enforcement_mode: active_set
    start_iter: 3
    start_mode: iteration_or_auto

matpower_import:
  case: case14.m
  # Non-empty cases take precedence for run_sparlectra_cases.
  cases: [case14.m, case118.m]
  auto_profile: recommend
  auto_profile_log: true
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
  casefile: ""
  case_name: ""
  case_source: ""
  configured_default_casefile: ""

extensions:
  reserved: true
```

`power_flow.qlimits.enforcement_mode` selects the reactive-limit algorithm.
`active_set` preserves Sparlectra's existing in-iteration PV→PQ active-set
behavior. `classic_simultaneous` and `classic_one_at_a_time` run a
classical/reference-style outer loop: solve the base AC power flow first with
Q-limit switching disabled, clamp violating generator Q to the violated limit
only after a successful solve, convert affected voltage-controlled buses to PQ,
and rerun without PQ→PV re-enable inside that enforcement loop. Legacy aliases
`matpower_simultaneous` and `matpower_one_at_a_time` are still accepted for old
YAML files and are normalized to the corresponding `classic_*` value.

`matpower_import.auto_profile` controls a MATPOWER pre-run profile. Use `off`
to disable it, `recommend` to print/log a recommendation table without changing
the active configuration, or `apply` to apply only unambiguous import-convention
recommendations. The runner also logs final effective MATPOWER options; it never
rewrites user YAML files.

For complete key references and allowed-value tables, see the module-specific pages below.

## Wrong-branch detection semantics (rectangular PF)

`power_flow.wrong_branch_detection` is a post-convergence plausibility check for rectangular PF results. It is heuristic and does **not** prove global branch correctness.

- `off`: checker disabled.
- `warn`: suspicious solutions are reported in rectangular status metadata; numerical convergence remains accepted.
- `fail`: suspicious solutions are treated as failed final convergence.
- `rescue`: reserved/request mode. If a suspicious solution is detected, status reports `rescue_requested_but_not_available`; active retry/rescue loops are not implemented yet.

The thresholds include voltage magnitude range, global angle spread, and active-branch angle-difference checks via `power_flow.wrong_branch_max_branch_angle_deg`.

### Where the result is visible

The check result — not just the setting — is surfaced in every output surface, so a suspicious solution is visible without reading console warnings:

- **`ACPFlowReport.metadata`** (`src/results.jl`, `buildACPFlowReport`): `wrong_branch_status` and `wrong_branch_reason`.
- **AC island diagnostics CSV** (`ac_island_solver_summary.csv`, one row per AC island): trailing `wrong_branch_status`/`wrong_branch_reason` columns, alongside the existing `wrong_branch_detection` *setting* column. The matching per-island `ac_island_<id>_solver.log` also lists both fields.
- **Console/log summary** (`printACPFlowResults`): a single `wrong-branch check: SUSPECT (...)` or `wrong-branch check: FAIL (...)` line is printed only when the result is neither `ok` nor `not_checked`; clean runs and `wrong_branch_detection = off` print nothing extra.
- **Web UI run result page**: a `status-badge` ("Wrong-branch check" row) using the same success/warning/error styling as the run-status badge; omitted entirely when the result is `not_checked`.
- **`run_sparlectra_api` result metadata**: `wrong_branch_status`, `wrong_branch_reason`, `wrong_branch_low_vm_count`, `wrong_branch_high_vm_count`, `wrong_branch_angle_spread_deg`, `wrong_branch_branch_angle_violation_count`.

`status` values: `ok` (checked, no finding), `warn` (suspicious but the numerical result was still accepted per `wrong_branch_detection = warn`), `fail` (suspicious and treated as non-convergence per `wrong_branch_detection = fail`), `wrong_branch_rescue_not_implemented` (the reserved `rescue` mode was requested; no retry loop runs — see below), or `not_checked` (`wrong_branch_detection = off`, or the check never ran, e.g. a non-finite solution short-circuits earlier reporting). `reason` values mirror the case listed above (`none`, `low_voltage_magnitude`, `high_voltage_magnitude`, `angle_spread_exceeded`, `branch_angle_exceeded`, `nonfinite_voltage`, `disabled`, `rescue_requested_but_not_available`).

**Maintainer decision:** the wrong-branch **rescue retry loop is intentionally out of scope** and will not be implemented as part of this detection work; `wrong_branch_detection = rescue` stays a reserved mode that reports `wrong_branch_rescue_not_implemented` rather than retrying. Detection with full output visibility (this section) plus the APSLF solver as an alternative start/solve path (see `docs/src/external_solvers.md`) are the supported mitigations for hard flat-start cases.

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
- `matpower_export`
- `matpower_export.write_solution`
- `matpower_import`
- `matpower_import.auto_profile`
- `matpower_import.auto_profile_log`
- `matpower_import.bus_shunt_model`
- `matpower_import.case`
- `matpower_import.cases`
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
- `power_flow.merit`
- `power_flow.merit.armijo_c1`
- `power_flow.merit.enabled`
- `power_flow.merit.fallback_max_mismatch`
- `power_flow.merit.scale_p`
- `power_flow.merit.scale_q`
- `power_flow.merit.scale_v`
- `power_flow.trust_region`
- `power_flow.trust_region.enabled`
- `power_flow.trust_region.initial_radius`
- `power_flow.trust_region.min_radius`
- `power_flow.trust_region.max_radius`
- `power_flow.trust_region.eta_accept`
- `power_flow.trust_region.shrink_factor`
- `power_flow.trust_region.expand_factor`
- `power_flow.trust_region.expand_threshold`
- `power_flow.trust_region.step_mode`
- `power_flow.method`
- `power_flow.qlimits`
- `power_flow.qlimits.auto_q_delta_pu`
- `power_flow.qlimits.cooldown_iters`
- `power_flow.qlimits.enabled`
- `power_flow.qlimits.enforcement_mode`
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
- `power_flow.start_current_iteration`
- `power_flow.start_current_iteration.accept_only_if_improved`
- `power_flow.start_current_iteration.damping`
- `power_flow.start_current_iteration.enabled`
- `power_flow.start_current_iteration.max_angle_step_deg`
- `power_flow.start_current_iteration.max_iter`
- `power_flow.start_current_iteration.min_improvement_factor`
- `power_flow.start_current_iteration.only_for_large_cases`
- `power_flow.start_current_iteration.tol`
- `power_flow.start_current_iteration.vm_max_pu`
- `power_flow.start_current_iteration.vm_min_pu`
- `power_flow.tol`
- `runtime`
- `runtime.blas_threads`
- `runtime.case_name`
- `runtime.case_source`
- `runtime.casefile`
- `runtime.configured_default_casefile`
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
- `transformer`
- `transformer.tap_changer_model`

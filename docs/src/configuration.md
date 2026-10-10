# Central Configuration

`Sparlectra.load_sparlectra_config(...)` loads the packaged defaults,
overlays a user YAML and programmatic overrides, validates keys and values,
and builds the typed `SparlectraConfig` every module reads.

## Configuration files and selectors

| Item | Path / mechanism | Role |
|---|---|---|
| Packaged defaults | `src/config/configuration.yaml.example` | Baseline for every key, commented. |
| User override | `examples/configuration.yaml` | Local override (default user path). |
| Explicit path | `load_sparlectra_config("/path/to/file.yaml")` | Replaces the default user path for that call. |
| Environment variable | `SPARLECTRA_CONFIGURATION_YAML` via `configuration_path_from_inputs(...)` | Path selection in scripts and examples. |

## Merge precedence

A later level wins:

| Level | Source | Notes |
|---|---|---|
| 1 | `src/config/configuration.yaml.example` | packaged defaults |
| 2 | user YAML: `examples/configuration.yaml` or an explicit `user_path` | |
| 3 | case configuration file next to the case | `<stem>.config.yaml` for an SCF case, else `<file name>.config.yaml` (`case118.m.config.yaml`), so a MATPOWER case and its SCF export keep separate files |
| 4 | Web UI form and runtime values of a service or Web UI run | |
| 5 | programmatic overrides | `cli_overrides`, then `overrides` (the API's `config_overrides`) |
| 6 | post-processing | such as `model.auto_profile: apply` |

An unknown or removed key in a file is warned about once (a removed key
with its migration hint) and ignored, so a file from an earlier release
keeps loading.

Case-scope keys are `matpower_import.*`, `cgmes_import.*`,
`powsybl_import.*`, `power_flow.*`, `state_estimation.*` and
`short_circuit.*` (`scf_is_case_config_key(key)` answers per key). When
a case ships its own configuration file, these keys skip level 2 and
resolve from the case file straight to the packaged defaults, so the case
computes the same numbers on every installation. Machine-scope keys (`output.*`, `benchmark.*`, `runtime.*`, `webui.*`,
`matpower_export.*`) still read the user file; a case file that states one
is refused by name. The in-file block `sparlectra.config` is deprecated and
applied below the case configuration file with a warning.
`effective_config.yaml` of each run records what took effect. See
[Configuration precedence](scf.md#Configuration-precedence).

## File version and scope

Every configuration file starts with its format version and scope:

```yaml
config_version: 1
scope: general
```

A file without `config_version` reads as version 0: it loads through the
version-0 aliases (`model.*` lived in `matpower_import` and `transformer`,
`runtime.case`/`runtime.cases` in `matpower_import`), and the loader names
the missing version and every translated key. A `config_version` newer
than the running Sparlectra is an error. `scope` is `general` for the
installation-wide file and `case` for a file next to its case.

`refresh_sparlectra_config_file(path; write = false)` checks a file against
the current template: user values are kept, missing keys come from the
template, duplicate keys are reported, `normalize_deprecated = true`
rewrites known aliases, and every key no loader reads any more (removed,
deprecated, unknown, a `form` block) is deleted and named. Nothing is
written without `write = true`; a write first saves a
`.bak-YYYYmmdd-HHMMSS` backup, and duplicate keys block it. The Web UI
button **Refresh configuration** and a [sysimage](sysimage.md) build call
it; nothing rewrites a user file at startup.

## YAML structure (section map)

The merged YAML becomes a `SparlectraConfig`; each module reads its own
typed section.

| YAML section | Typed section | Purpose |
|---|---|---|
| `power_flow` | `PowerFlowConfig` | Solver, start mode, Q limits, `power_flow.mode` (`manual` or `auto`): [Power-Flow Configuration](powerflow_configuration.md), [Integration Guide](integration.md). |
| `matpower_import` | `MatpowerImportConfig` | Import conventions, the metadata keys `matpower_import.apply_bus_names`, `apply_branch_names`, `apply_branch_kind`, `import_for001_contingencies` and `matpower_import.matpower_dcline_mode`: [Option reference](@ref matpower-options). |
| `cgmes_import` | `CGMESImportConfig` | [CGMES Import](cgmes_import.md). |
| `powsybl_import` | `PowsyblImportConfig` | [PowSyBl Import](powsybl_import.md). |
| `short_circuit` | `ShortCircuitConfig` | IEC 60909 short circuit: `short_circuit.c_factor`, `short_circuit.sweep_method` (`auto`, `solves`, `takahashi`), `short_circuit.takahashi_min_buses`; [Short-Circuit Analysis](short_circuit.md). |
| `model` | `ModelConfig` | Model construction for every importer: bus shunt model, tap-changer model (below), auto profile, net cache, preallocation. |
| `state_estimation` | `StateEstimationConfig` | [State-Estimation Configuration](state_estimation_configuration.md). |
| `output` | `OutputConfig` | Console and logfile output, result tables (keys below and in [Output configuration](@ref perf-output)); `output.startup_latency_hint: false` silences the one-per-process JIT warm-up note of runs without a sysimage or executable. |
| `performance` | `PerformanceConfig` | Profiling and diagnostic volume: [Performance and Profiling](performance_profiling.md). |
| `benchmark` | `BenchmarkConfig` | Repeated benchmark runs: [Benchmark configuration](@ref perf-benchmark). |
| `contingency` | `ContingencyConfig` | N-1 batch: `contingency.rescue_ladder`, `contingency.max_iter`, `contingency.screening.mode`, `contingency.screening.margin_pct`, `contingency.warm_active_set`, `contingency.warm_cold_check`, `contingency.warm_cold_check_margin_pu`; [N-1 Contingency Analysis](contingency.md). |
| `control` | `ControlConfig` | Controller outer loop and declarative controllers, below. |
| `runtime` | `RuntimeConfig` | Case selection (`runtime.case`, `runtime.cases`), thread control, `runtime.parallel.*` ([Parallel execution](parallel_execution.md)). |
| `diagnostics` | `DiagnosticsConfig` | `log_effective_config` only; the former `diagnostics.console_*`/`logfile_diagnostics` duplicate `output.*` and are ignored with a warning. |
| `webui` | `WebUIConfig` | Web UI preferences such as `webui.show_case_settings_notice`; keys below. |
| `extensions` | reserved | Placeholder, not mapped to typed fields. |

## Option lifecycle and compatibility policy

Prefer the nested keys of the template.

| Class | Meaning |
|---|---|
| Public / supported | keys in `src/config/configuration.yaml.example` and the typed constructors |
| Reserved | placeholders such as `extensions.reserved` |
| Deprecated aliases | accepted for migration, not for new files (`max_ite`, some start-projection aliases) |
| Removed | warn with a migration hint, ignored (`matpower_import.benchmark`) |
| Internal | not a stable API |

### Transformer tap-changer model

`model.tap_changer_model` decides for every transformer of an imported case
whether a tap step changes the winding ratio alone or also the series
impedance:

```yaml
model:
  # Allowed values: ideal, impedance_correction
  tap_changer_model: ideal
```

| Value | Effect |
|---|---|
| `ideal` (default) | Tap steps change the complex winding ratio; R and X keep their neutral-position values. |
| `impedance_correction` | Tap steps also re-refer R and X through the tapped winding; the correction math is in [`calcTapCorrectedRX`](@ref) / [`calcTapImpedanceCorrectionFactor`](@ref), see the [branch model](branchmodel.md). |

The MATPOWER and DTF importers read the key. [`writeMatpowerCasefile`](@ref)
writes the corrected `Branch.r_pu`/`Branch.x_pu` with a roundtrip marker,
so a reimport does not correct twice; see
[Tap-impedance correction and reimport](matpower.md#tap-impedance-correction-and-reimport).

## Loader and validation behavior

- Keys are validated against the schema tree of the template: an unknown
  key in a file is warned about and dropped, an unknown override key throws
  `ArgumentError`. The Web UI configuration editor rejects an unknown key
  while the text is edited.
- Type and domain checks (allowed symbols, positivity) run while the typed
  objects are built; some legacy aliases are accepted there.
- `load_sparlectra_config(...)` caches the typed result for unchanged files
  without overrides; a changed hash or mtime invalidates it.

## Minimal example YAML

```yaml
config_version: 1
scope: general

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
  qlimits:
    enabled: true
    enforcement_mode: active_set

model:
  auto_profile: off

matpower_import:
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

runtime:
  case: case14.m
  cases: [case14.m, case118.m]   # non-empty wins for run_sparlectra_cases
  julia_threads: keep
  blas_threads: keep
```

## [Wrong-branch detection semantics (rectangular PF)](@id config-wrong-branch)

`power_flow.wrong_branch_detection` is a heuristic plausibility check of
every numerically converged AC result, with or without reactive limits; it
proves nothing about the solution branch.

| Mode | Effect |
|---|---|
| `off` | no check |
| `warn` (default) | a suspicious solution is reported in the status metadata, the run stays converged |
| `fail` | a suspicious solution counts as non-convergence |
| `rescue` | reserved: reports `rescue_requested_but_not_available`, no retry runs |

The magnitude band (`wrong_branch_min_vm_pu` to `wrong_branch_max_vm_pu`)
and the level-share rule are judged on every energised bus whose nominal
voltage is at or above `wrong_branch_min_vn_kV`; when no level reaches the
floor, on the highest level alone. A bus below
`wrong_branch_collapse_vm_pu` is a finding on every level. The angle spread
and the branch-angle rule stay on the highest level (both ends of the branch
on that level). The angle across every in-service branch at the reference
bus is judged on every level, phase shift compensated
(`reference_branch_angle_exceeded`); it catches a solution rotated as a
whole against a reference machine behind its own transformer. The reported
magnitudes, lowest-bus list
and `wrong_branch_level_kV` with its counts refer to the judged levels. A
run with more than `wrong_branch_max_plain_steps` Newton steps without a
switching event is reported (`wrong_branch_plain_steps_exceeded`) without
changing the status; a well-posed case takes 7 to 14 steps.

!!! details "Why the angle rules stay on the highest voltage level"
    Sub-transmission levels run at 0.94 to 0.97 pu in healthy snapshots
    with their own angle spread. The magnitude band (0.70 pu) is far below
    that and is judged on every transmission level, because a collapse can
    sit on a lower level under a clean top level (a 13659-bus case ended
    with 130 buses of its 150 kV level at 0.33 pu under three clean 750 kV
    buses).

### Where the result is visible

| Surface | Fields |
|---|---|
| `ACPFlowReport.metadata` | `wrong_branch_status`, `wrong_branch_reason` |
| `run_sparlectra_api` result metadata | the two fields above plus `wrong_branch_low_vm_count`, `wrong_branch_high_vm_count`, `wrong_branch_level_kV`, `wrong_branch_level_low_vm_count`, `wrong_branch_level_bus_count`, `wrong_branch_max_bus_angle_deg`, `wrong_branch_plain_steps`, `wrong_branch_plain_steps_exceeded`, `wrong_branch_angle_spread_deg`, `wrong_branch_branch_angle_violation_count`, `wrong_branch_reference_branch_angle_deg` with `wrong_branch_reference_branch` (the largest angle across a branch at the reference bus and that branch in case bus numbers) |
| `ac_island_solver_summary.csv`, `ac_island_<id>_solver.log` | `wrong_branch_status`/`wrong_branch_reason` per island, next to the `wrong_branch_detection` setting column |
| `printACPFlowResults` | a `Wrong-branch   : SUSPECT (...)` or `FAIL (...)` line unless the result is `ok` or `not_checked` |
| Web UI run result page | a "Wrong-branch check" badge unless `not_checked` |

| `status` value | Meaning |
|---|---|
| `ok` | checked, no finding |
| `warn` | suspicious, accepted under `warn` |
| `fail` | suspicious, non-convergence under `fail` |
| `wrong_branch_rescue_not_implemented` | the reserved `rescue` mode was requested |
| `not_checked` | `off`, or the check never ran (non-finite solution) |

`reason` values: `none`, `voltage_collapse`, `low_voltage_level_share`,
`low_voltage_magnitude`, `high_voltage_magnitude`, `bus_angle_exceeded`,
`angle_spread_exceeded`, `branch_angle_exceeded`,
`reference_branch_angle_exceeded`, `nonfinite_voltage`, `disabled`,
`rescue_requested_but_not_available`.

The reserved `rescue` mode is not the rescue ladder for failed AC solves
(`power_flow.rescue`, [Solver core options](@ref pf-solver-core)). For hard
flat-start cases the mitigations are this check and the APSLF solver
([External Solvers](external_solvers.md)).

Tuning keys (all under `power_flow.`):

| Key | Default | Meaning |
|---|---|---|
| `power_flow.wrong_branch_min_vm_pu` | `0.70` | Lower edge of the plausibility band. |
| `power_flow.wrong_branch_max_vm_pu` | `1.30` | Upper edge of the plausibility band. |
| `power_flow.wrong_branch_min_low_vm_count` | `1` | Sub-band buses it takes to raise the finding. |
| `power_flow.wrong_branch_min_vn_kV` | `100.0` | Nominal-voltage floor of the judged levels (magnitude band and level share). |
| `power_flow.wrong_branch_low_vm_share` | `0.05` | Share of a level's buses below the band that raises `low_voltage_level_share`. |
| `power_flow.wrong_branch_collapse_vm_pu` | `0.5` | A bus below this magnitude on any level raises `voltage_collapse`; 0 switches the rule off. |
| `power_flow.wrong_branch_max_angle_spread_deg` | `180.0` | Maximum total angle spread of the solution. |
| `power_flow.wrong_branch_max_branch_angle_deg` | `90.0` | Bound on the angle across an active branch of the highest level, phase shift compensated. |
| `power_flow.wrong_branch_max_bus_angle_deg` | `120.0` | Bound on any judged bus angle relative to the slack (`bus_angle_exceeded`). |
| `power_flow.wrong_branch_max_reference_branch_angle_deg` | `90.0` | Bound on the angle across any in-service branch at the reference bus, every level, phase shift compensated (`reference_branch_angle_exceeded`, ranked before `bus_angle_exceeded`); 0 switches the rule off. |
| `power_flow.wrong_branch_max_plain_steps` | `20` | Newton steps without a switching event above which `wrong_branch_plain_steps_exceeded` is reported. |
| `power_flow.wrong_branch_rescue` | `false` | Reserved switch of the unimplemented rescue loop; reports instead of retrying. |
| `power_flow.wrong_branch_rescue_max_attempts` | `2` | Attempt bound of that reserved mode. |

## Control configuration (generic outer loop)

```yaml
control:
  enabled: true
  max_outer_iterations: 20
  trace: true
  log_iterations: true
  stop_on_pf_failure: true
  verbose_passes: false
  controllers: {}
```

| Key | Type | Default | Meaning |
|---|---:|---:|---|
| `control.enabled` | Bool | `true` | Run the outer control loop (`run_control!`) around the inner solver when controllers exist. |
| `control.max_outer_iterations` | Int | `20` | Outer-loop budget shared by all controllers; independent of `power_flow.max_iter`. |
| `control.trace` | Bool | `true` | Collect machine-readable rows in `ControlRunResult.trace`. |
| `control.log_iterations` | Bool | `true` | One console line per control pass when `output.console_diagnostics` is `full`. |
| `control.stop_on_pf_failure` | Bool | `true` | Abort when the inner power flow fails. |
| `control.verbose_passes` | Bool | `false` | Repeat the full inner-solver diagnostic blocks on every pass; off prints them once and summarizes later passes in one line. |
| `control.controllers` | Mapping | `{}` | Declarative controllers, applied to the net before the outer loop (below; mechanics in [Control Framework](control_framework.md)). |

### Declarative controllers (`control.controllers`)

One named entry per controller; `type` selects the device function, the
other keys mirror its keyword arguments (bus, branch and transformer
references by name). Block style only: the minimal YAML reader has no
`- item` sequences.

```yaml
control:
  controllers:
    tap_T1:
      type: power_transformer
      trafo: T1
      mode: voltage
      target_bus: B2
      target_vm_pu: 1.02
    tcsc_B1_B2:
      type: series_reactance
      from_bus: B1
      to_bus: B2
      p_target_mw: 80.0
      x_min_pu: 0.05
      x_max_pu: 0.30
```

| `type` | Device function | Notes |
|---|---|---|
| `power_transformer` | `addPowerTransformerControl!` | |
| `machine_voltage` | `addMachineVoltageControl!` | |
| `shunt_voltage` | `addShuntVoltageControl!` | |
| `series_reactance` | `addSeriesReactanceControl!` | |
| `hvdc_pair` | `addHvdcPairControl!` | |
| `upfc` | `addUpfcControl!` | `model: quadrature` registers the SSSC+STATCOM pair, `model: full` the DC-link-coupled independent-P/Q model, see [FACTS Devices](facts.md) |

Unknown types, unknown keys and missing required keys fail at configuration
load; unresolvable references and invalid limits fail at apply time, naming
the entry. An element that already carries a controller of the same type
is skipped. `applyConfiguredControllers!` applies the block to a
programmatically built net.

## Bookkeeping and console keys

| Key | Default | Meaning |
|---|---|---|
| `runtime.casefile`, `runtime.case_name`, `runtime.case_source`, `runtime.configured_default_casefile` | `""` | Written by the API and the Web UI into `effective_config.yaml`: which case a run used and where it came from. Not user inputs. |
| `output.console_live` | `false` | Mirror the captured run output live to the console during API runs; `run.log` stays identical. |
| `output.result_table_max_rows` | `200` | Row cap of the classical result tables. |
| `output.result_table_large_case_threshold_buses` | `1000` | Bus count from which a case counts as large for result rendering. |
| `output.result_table_large_case_mode` | `summary` | What large cases print instead of full tables: `summary`, `classic`, `full`. |
| `output.csv_format` | `auto` | Delimiter and decimal separator of every CSV file a run writes; only the measurement CSV Sparlectra reads back keeps its fixed layout. `auto` follows the regional settings of the machine (decimal comma gives `excel_de`, decimal point `excel_us`, the C/POSIX locale `technical`); `technical` is comma delimiter, dot decimal; `excel_de` semicolon delimiter, comma decimal, dot thousands separator; `excel_us` comma delimiter, dot decimal, comma thousands separator. The API keywords `detailed_result_csv_format`/`detailed_result_csv_semicolon` are a deprecated per-request override. |
| `webui.operation_log_retention_days` | `10` | Int `>= 0`: days the operation log reaches back. Every Web UI start drops older entries; `0` keeps only the current session. The environment variable `SPARLECTRA_WEBUI_OPERATION_LOG_RETENTION_DAYS` wins (headless runs without a configuration file). |
| `webui.docs_base_url` | `https://welthulk.github.io/Sparlectra.jl/` | Base URL of the published documentation that the help pages (the **?** next to a control) and the header link open; point it at a local build or a pinned version. A link only, nothing is fetched. |

## Migration notes

| Legacy / old key | Canonical key | Notes |
|---|---|---|
| `matpower_import.benchmark` | `benchmark.enabled` | removed, warns with this hint |
| `methods` (top-level legacy path) | `benchmark.methods` | |
| `max_ite` | `power_flow.max_iter` | deprecated alias |

## Detailed references

[Power-Flow Configuration](powerflow_configuration.md),
[State-Estimation Configuration](state_estimation_configuration.md),
[MATPOWER cases](matpower.md),
[Performance and Profiling Configuration](performance_profiling.md).
The canonical key set is
[`src/config/configuration.yaml.example`](https://github.com/Welthulk/Sparlectra.jl/blob/main/src/config/configuration.yaml.example);
`print_effective_config` prints a run's effective configuration.

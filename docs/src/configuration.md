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
  - `runtime::RuntimeConfig` (incl. `runtime.parallel.*` as `ParallelRuntimeConfig`)
  - `diagnostics::DiagnosticsConfig`
  - `output::OutputConfig`
  - `control::ControlConfig`
  - `cgmes::CGMESImportConfig`
  - `shortcircuit::ShortCircuitConfig`

This typed model is the canonical internal representation that should be consumed by power-flow, MATPOWER import, state estimation, output/reporting, performance profiling, benchmark runners, and future modules.

## YAML structure (section map)

| YAML section | Typed section | Purpose | Status |
|---|---|---|---|
| `power_flow` | `PowerFlowConfig` | Rectangular power-flow solver controls, start mode, Q-limits | Public / supported |
| `matpower_import` | `MatpowerImportConfig` | MATPOWER case path + import interpretation options | Public / supported |
| `cgmes_import` | `CGMESImportConfig` | CGMES delivery path + import options (see [CGMES Import](cgmes_import.md)) | Public / supported |
| `short_circuit` | `ShortCircuitConfig` | IEC 60909 short-circuit evaluation; `short_circuit.c_factor` overrides the Table-1 voltage factor (see [Short-Circuit Analysis](short_circuit.md)) | Public / supported |
| `transformer` | `TransformerConfig` | Transformer-modeling options shared by all importers (tap-changer model) | Public / supported |
| `state_estimation` | `StateEstimationConfig` | State-estimation runtime controls | Public / supported |
| `output` | `OutputConfig` | Console/logfile behavior and result table sizing | Public / supported |
| `performance` | `PerformanceConfig` | Profiling/reporting toggles and diagnostic volume controls | Public / supported |
| `benchmark` | `BenchmarkConfig` | Repeated benchmark-run controls | Public / supported |
| `control` | `ControlConfig` | Generic controller outer-loop orchestration controls; `control.controllers` holds declarative controller definitions (implemented, issue #305), see the controllers section below | Public / supported |
| `runtime` | `RuntimeConfig` | Julia/BLAS thread control knobs for entry workflows | Public / supported |
| `diagnostics` | `DiagnosticsConfig` | Effective-config logging (`log_effective_config` only — the former `diagnostics.console_*`/`logfile_diagnostics` duplicates of `output.*` are deprecated and ignored with a warning) | Public / supported |
| `webui` | `WebUIConfig` | Web UI presentation preferences (e.g. `webui.show_case_settings_notice`) | Public / supported |
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

All heuristics are evaluated **only on the network's highest nominal voltage
level** (branch-angle checks: both ends on that level). Sub-transmission
levels routinely run at 0.94–0.97 pu in healthy snapshots; judging them
against the transmission-level band produced false `SUSPECT` verdicts (seen
on a real 6209-bus CGMES delivery whose 45 kV feeders were flagged while the
380 kV level was clean). The reported `min_vm_pu`/`max_vm_pu` and the
lowest-bus list refer to the checked level.

### Where the result is visible

The check result — not just the setting — is surfaced in every output surface, so a suspicious solution is visible without reading console warnings:

- **`ACPFlowReport.metadata`** (`src/results.jl`, `buildACPFlowReport`): `wrong_branch_status` and `wrong_branch_reason`.
- **AC island diagnostics CSV** (`ac_island_solver_summary.csv`, one row per AC island): trailing `wrong_branch_status`/`wrong_branch_reason` columns, alongside the existing `wrong_branch_detection` *setting* column. The matching per-island `ac_island_<id>_solver.log` also lists both fields.
- **Console/log summary** (`printACPFlowResults`): a single `Wrong-branch   : SUSPECT (...)` or `Wrong-branch   : FAIL (...)` line (label aligned with the other classical header lines) is printed only when the result is neither `ok` nor `not_checked`; clean runs and `wrong_branch_detection = off` print nothing extra.
- **Web UI run result page**: a `status-badge` ("Wrong-branch check" row) using the same success/warning/error styling as the run-status badge; omitted entirely when the result is `not_checked`.
- **`run_sparlectra_api` result metadata**: `wrong_branch_status`, `wrong_branch_reason`, `wrong_branch_low_vm_count`, `wrong_branch_high_vm_count`, `wrong_branch_angle_spread_deg`, `wrong_branch_branch_angle_violation_count`.

`status` values: `ok` (checked, no finding), `warn` (suspicious but the numerical result was still accepted per `wrong_branch_detection = warn`), `fail` (suspicious and treated as non-convergence per `wrong_branch_detection = fail`), `wrong_branch_rescue_not_implemented` (the reserved `rescue` mode was requested; no retry loop runs — see below), or `not_checked` (`wrong_branch_detection = off`, or the check never ran, e.g. a non-finite solution short-circuits earlier reporting). `reason` values mirror the case listed above (`none`, `low_voltage_magnitude`, `high_voltage_magnitude`, `angle_spread_exceeded`, `branch_angle_exceeded`, `nonfinite_voltage`, `disabled`, `rescue_requested_but_not_available`).

The rescue retry loop is intentionally not implemented. `wrong_branch_detection = rescue` is a reserved mode that reports `wrong_branch_rescue_not_implemented` rather than retrying. Supported mitigations for hard flat-start cases are detection with full output visibility (this section) and the APSLF solver as an alternative start/solve path (see [External Solvers](external_solvers.md)).

Tuning keys of the detector (all under `power_flow.`):

| Key | Default | Meaning |
|---|---|---|
| `power_flow.wrong_branch_min_vm_pu` | `0.70` | Lower edge of the plausibility band; solved magnitudes below it count as suspicious. |
| `power_flow.wrong_branch_max_vm_pu` | `1.30` | Upper edge of the plausibility band. |
| `power_flow.wrong_branch_min_low_vm_count` | `1` | How many sub-band buses it takes to raise the finding. |
| `power_flow.wrong_branch_max_angle_spread_deg` | `180.0` | Maximum admissible total angle spread of the solution. |
| `power_flow.wrong_branch_rescue` | `false` | Reserved switch for the unimplemented rescue loop; reports instead of retrying. |
| `power_flow.wrong_branch_rescue_max_attempts` | `2` | Attempt bound for that reserved mode. |

## Control configuration (generic outer loop)

```yaml
control:
  enabled: true
  max_outer_iterations: 20
  trace: true
  log_iterations: true
  stop_on_pf_failure: true
  controllers: {}
```

| Key | Type | Default | Meaning |
|---|---:|---:|---|
| `enabled` | Bool | `true` | Enables the generic outer-loop control framework. |
| `max_outer_iterations` | Int | `20` | Global outer-loop cap. Does not control inner NR iterations. |
| `trace` | Bool | `true` | Collect machine-readable control trace rows. |
| `log_iterations` | Bool | `true` | Enables optional per-iteration control logging hooks. |
| `stop_on_pf_failure` | Bool | `true` | Stops control orchestration when inner PF fails. |
| `controllers` | Mapping | `{}` | Declarative controller instantiation (issue #305): one named mapping per controller, applied to the net before the outer loop. See [Control Framework](control_framework.md) for the schema. |

### Declarative controllers (`control.controllers`)

One named entry per controller; the `type` key selects the device function,
the remaining keys mirror its keyword arguments (bus/branch/transformer
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

Supported types: `power_transformer` (`addPowerTransformerControl!`),
`machine_voltage` (`addMachineVoltageControl!`), `shunt_voltage`
(`addShuntVoltageControl!`), `series_reactance`
(`addSeriesReactanceControl!`), `hvdc_pair`
(`addHvdcPairControl!`). Unknown types or keys and missing required
keys fail at configuration load; unknown bus/branch/transformer references
and invalid limits fail at apply time naming the entry. Entries whose
element already carries a controller of the same type are skipped, so
repeated runs on one net stay idempotent
(`applyConfiguredControllers!` applies a configuration to a
programmatically built net directly).

### Demo controller YAML vs. central `control.controllers`

The tap-control demo may read `examples/others/tap_control_demo_grid.yaml` for
example setpoints and transformer tap/phase parameters (`oltc`, `pst`,
`schraeg`). This is an example-specific
input file consumed by `examples/others/tap_control_demo_grid.jl` and does
not define the central `control.controllers` schema above.

The section's plain switches:

| Key | Default | Meaning |
|---|---|---|
| `control.enabled` | `true` | Run the outer control loop (`run_control!`) around the inner solver when controllers exist. |
| `control.max_outer_iterations` | `20` | Outer-loop budget shared by all active controllers. |
| `control.trace` | `true` | Collect machine-readable rows in `ControlRunResult.trace`. |
| `control.log_iterations` | `true` | Reserved for per-iteration control logging. |
| `control.stop_on_pf_failure` | `true` | Abort the orchestration when the inner power flow fails. |

## Bookkeeping and console keys

A few keys exist for run bookkeeping and console behavior rather than for
solver tuning:

| Key | Default | Meaning |
|---|---|---|
| `runtime.casefile`, `runtime.case_name`, `runtime.case_source`, `runtime.configured_default_casefile` | `""` | Populated by the API/Web UI into the `effective_config.yaml` run artifact; they record which case a run actually used and where it came from. Not user inputs. |
| `output.console_live` | `false` | Mirror captured run output live to the real console during API/service runs; the archived `run.log` stays identical. |
| `output.result_table_max_rows` | `200` | Row cap for the classical result tables. |
| `output.result_table_large_case_threshold_buses` | `1000` | Bus count from which a case counts as large for result rendering. |
| `output.result_table_large_case_mode` | `summary` | What large cases print instead of full tables (`summary`, `classic`, `full`). |
| `webui.warmup` | `true` | Compile the hot paths once at Web UI start (a few seconds) so the first run does not pay them; `false` for a minimal-footprint start. |

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

The canonical key set lives in ONE place,
[`src/configuration.yaml.example`](https://github.com/Welthulk/Sparlectra.jl/blob/main/src/configuration.yaml.example):
the version-controlled default configuration that the loader validates
every user YAML against. A hand-maintained copy of that list used to sit
here and drifted (it was missing about sixty keys); read the file itself,
it is commented per key and always current. The effective configuration
of a run, with every default and override resolved, prints via
`print_effective_config`.

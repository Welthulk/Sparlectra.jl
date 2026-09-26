# Central Configuration

Sparlectra uses one central, typed configuration model. The entry point
`Sparlectra.load_sparlectra_config(...)` loads the YAML defaults, overlays a
user YAML, applies programmatic overrides, validates key names and value
domains, and builds the typed `SparlectraConfig` the runtime modules read.

## Configuration files and selectors

| Item | Path / mechanism | Role |
|---|---|---|
| Default config | `src/config/configuration.yaml.example` | Version-controlled baseline for all options. |
| User override | `examples/configuration.yaml` | Local override file (project default user path). |
| Explicit config path | `load_sparlectra_config("/path/to/file.yaml")` | Replaces default user override path for that load call. |
| Environment-based path selection helper | `SPARLECTRA_CONFIGURATION_YAML` via `configuration_path_from_inputs(...)` | Used by script/example path resolution workflows. |

## Merge precedence

A later level wins:

| Level | Source | Notes |
|---|---|---|
| 1 | `src/config/configuration.yaml.example` | packaged defaults |
| 2 | user YAML: `examples/configuration.yaml` or an explicit `user_path` | |
| 3 | case configuration file next to the case | `<stem>.config.yaml` for an SCF case, `<file name>.config.yaml` (such as `case118.m.config.yaml`) for every other format, so a MATPOWER case and its SCF export keep separate files |
| 4 | Web UI form and runtime values of a service or Web UI run | |
| 5 | programmatic overrides | `cli_overrides`, then `overrides` (the API's `config_overrides`) |
| 6 | post-processing | such as `model.auto_profile: apply` |

Unknown keys are rejected during validation; removed keys are rejected with
a migration hint (for example `matpower_import.benchmark` to
`benchmark.enabled`).

When a case ships its own configuration file, case-scope keys skip level 2
and resolve from the case file straight to the packaged defaults, so the
same case-plus-config pair computes the same numbers on every installation;
only machine-scope keys (`output.*`, `benchmark.*`, `runtime.*`, `webui.*`,
`matpower_export.*`) still read the user file. See
[Configuration precedence](scf.md#Configuration-precedence).

The case file's own `sparlectra.config` block is deprecated: the writer
does not emit it, the reader still applies it below the case
configuration file (between levels 3 and 1 for case-scope keys) with a
warning naming the case configuration file as the new place. Case-scope
keys are `matpower_import.*`, `cgmes_import.*`, `powsybl_import.*`, `power_flow.*`,
`state_estimation.*` and `short_circuit.*`; `scf_is_case_config_key(key)`
answers it for a single key, and a case file that states a key outside
that scope is refused by name when it is read. `effective_config.yaml`
of each run records what took effect.

## File version and scope

Every configuration file declares its format as the first keys:

```yaml
config_version: 1
scope: general
```

A file without `config_version` reads as version 0: it still loads, the
version-0 aliases are applied (the `model.*` keys lived in
`matpower_import`/`transformer`, `runtime.case`/`runtime.cases` in
`matpower_import`), and the loader reports the missing version and every
translated legacy key. `refresh_sparlectra_config_file` rewrites such a file
to the current layout; a [sysimage](sysimage.md) build does it for the Web
UI configuration. A `config_version` newer than the running Sparlectra is
an error. `scope` is `general` for the installation-wide file, `case` for a
per-case file next to its case.

## Typed central object

The merged YAML is converted into `SparlectraConfig`; every module reads its
own section, and a new module adds its own typed section and YAML subtree.

| `SparlectraConfig` field | Read by |
|---|---|
| `powerflow::PowerFlowConfig` | power flow (`config.powerflow`) |
| `state_estimation::StateEstimationConfig` | state estimation (`config.state_estimation`) |
| `matpower::MatpowerImportConfig` | MATPOWER import (`config.matpower`) |
| `model::ModelConfig` | all importers |
| `performance::PerformanceConfig`, `benchmark::BenchmarkConfig` | profiling and benchmarks (`config.performance`, `config.benchmark`) |
| `contingency::ContingencyConfig` | N-1 batches |
| `runtime::RuntimeConfig` (incl. `runtime.parallel.*` as `ParallelRuntimeConfig`) | entry workflows, thread control |
| `diagnostics::DiagnosticsConfig`, `output::OutputConfig` | output and reporting (`config.output`, `config.diagnostics`) |
| `control::ControlConfig` | controller outer loop |
| `cgmes::CGMESImportConfig` | CGMES import |
| `shortcircuit::ShortCircuitConfig` | short-circuit evaluation |

## YAML structure (section map)

| YAML section | Typed section | Purpose | Status |
|---|---|---|---|
| `power_flow` | `PowerFlowConfig` | Rectangular power-flow solver controls, start mode, Q-limits; `power_flow.mode` switches between `manual` and the network-driven `auto` strategy (see [Power-Flow Configuration](powerflow_configuration.md) and the [Integration Guide](integration.md)) | Public / supported |
| `matpower_import` | `MatpowerImportConfig` | MATPOWER import interpretation options | Public / supported |
| `cgmes_import` | `CGMESImportConfig` | CGMES delivery path + import options (see [CGMES Import](cgmes_import.md)) | Public / supported |
| `powsybl_import` | `PowsyblImportConfig` | PowSyBl IIDM import options (see [PowSyBl Import](powsybl_import.md)) | Public / supported |
| `short_circuit` | `ShortCircuitConfig` | IEC 60909 short-circuit evaluation; `short_circuit.c_factor` overrides the Table-1 voltage factor, `short_circuit.sweep_method` (`auto`/`solves`/`takahashi`) selects the all-bus Thevenin sweep and `short_circuit.takahashi_min_buses` the island size from which the selected inverse is used (see [Short-Circuit Analysis](short_circuit.md)) | Public / supported |
| `model` | `ModelConfig` | Model construction shared by all importers (bus shunt model, tap-changer model, auto profile, net cache, preallocation) | Public / supported |
| `state_estimation` | `StateEstimationConfig` | State-estimation runtime controls | Public / supported |
| `output` | `OutputConfig` | Console/logfile behavior and result table sizing; `output.startup_latency_hint` silences the one-per-process note about JIT warm-up in sessions without a sysimage or executable (see [Sysimage](sysimage.md)) | Public / supported |
| `performance` | `PerformanceConfig` | Profiling/reporting toggles and diagnostic volume controls | Public / supported |
| `benchmark` | `BenchmarkConfig` | Repeated benchmark-run controls | Public / supported |
| `contingency` | `ContingencyConfig` | N-1 contingency batch controls; `contingency.rescue_ladder` is the per-case start-value ladder (subset of `warm`/`apslf`/`dc`/`flat`), `contingency.screening.mode` (`off`/`flag`/`only`, default `off`; screening is a deliberate opt-in) switches the base-factorization outage screening on the service path, and `contingency.screening.margin_pct` (default `10.0`) is its flagging margin; see [N-1 Contingency Analysis](contingency.md) | Public / supported |
| `control` | `ControlConfig` | Generic controller outer-loop orchestration controls; `control.controllers` holds declarative controller definitions, see the controllers section below | Public / supported |
| `runtime` | `RuntimeConfig` | Case selection (`runtime.case`/`runtime.cases`) and Julia/BLAS thread control knobs for entry workflows | Public / supported |
| `diagnostics` | `DiagnosticsConfig` | Effective-config logging (`log_effective_config` only; `diagnostics.console_*`/`logfile_diagnostics` duplicate `output.*`, are deprecated and ignored with a warning) | Public / supported |
| `webui` | `WebUIConfig` | Web UI presentation preferences (e.g. `webui.show_case_settings_notice`) | Public / supported |
| `extensions` | reserved (not mapped to typed runtime fields) | Future extension placeholder | Reserved |

The supported public power-flow solver path is rectangular
(`power_flow.method: rectangular`).

## Option lifecycle and compatibility policy

Prefer the canonical nested keys of the example YAML and the module pages.

| Class | Meaning |
|---|---|
| Public / supported | keys in `src/config/configuration.yaml.example` and typed section constructors |
| Reserved | schema placeholders such as `extensions.reserved` for forward compatibility |
| Deprecated compatibility aliases | accepted for migration but not preferred in new YAML (for example `max_ite` and some start-projection alias keys) |
| Removed | explicitly rejected keys with migration error guidance (for example `matpower_import.benchmark`) |
| Internal-only implementation details | not a stable external user API |

### MATPOWER import metadata and DC-line options

`matpower_import.apply_bus_names`, `apply_branch_names`, and
`apply_branch_kind` default to `false` (numeric bus names, heuristic branch
classification); enable them for cases that carry `mpc.bus_name` and
user-defined `mpc.branch_name`/`mpc.branch_kind`.
`import_for001_contingencies` (default `true`) preserves user-defined
`mpc.for001_contingencies` for validation workflows.
`matpower_import.matpower_dcline_mode` defaults to `:pf_injections`
(DC-line terminals as simple power-flow injections); `:reject_active`
rejects active DC-line rows. OPF and `dclinecost` are unsupported.

### AC island solving

MATPOWER `mpc.dcline` terminal injections add no Ybus branches, so they do
not tie disconnected AC islands together. `power_flow.islands.enabled: true`
is the default: disconnected AC components are solved independently.

```yaml
power_flow:
  islands:
    enabled: true
    mode: solve_independent
    reference_policy: matpower_like
```

The `matpower_like` policy keeps an existing island REF/Slack bus, otherwise
promotes the deterministic first PV/voltage-controlled bus as the island's
angle reference; islands without REF/Slack or PV support fail before NR.
With multiple islands, `ac_islands.csv` in the run directory lists bus,
branch, generator/load, DC-line terminal, power-balance, reference, and
status diagnostics.

### Transformer tap-changer model

`model.tap_changer_model` selects, for all transformers of an imported case,
whether the tap changer is electrically ideal or affects the series
impedance:

```yaml
model:
  # Allowed values: ideal, impedance_correction
  tap_changer_model: ideal
```

| Value | Effect |
|---|---|
| `ideal` (default) | Tap steps only change the complex winding ratio; the series impedance (R, X) keeps its neutral-position value. |
| `impedance_correction` | Tap steps additionally re-refer the series impedance through the tapped winding (R and X scaled, see below). |

With `impedance_correction`, R and X are scaled with the squared magnitude
of the regulating vector, `|1 + f·e^(jφ)|²`, where `f` is the longitudinal
regulating-voltage fraction and `φ` the skew angle (0° for a pure
longitudinal/ratio tap changer). Read by the MATPOWER importer
(`createNetFromMatPowerFile`/`createNetFromMatPowerCase`) and the native DTF
importer (`DTFImporter.build_net`/`createNetFromDTFFile`); the correction
math lives once, in [`calcTapCorrectedRX`](@ref) /
[`calcTapImpedanceCorrectionFactor`](@ref) (`src/equicircuit.jl`). A
subsequent [`writeMatpowerCasefile`](@ref) export writes the corrected
`Branch.r_pu`/`Branch.x_pu` values and records the roundtrip marker
`mpc.sparlectra.tap_changer_model = 'impedance_correction'` so a reimport
does not reapply the correction; see
[Tap-impedance correction and reimport](matpower.md#tap-impedance-correction-and-reimport).

## Loader and validation behavior

- User YAML and override keys are validated against the schema tree of
  `src/config/configuration.yaml.example`; unknown keys throw `ArgumentError`.
- Type and domain checks (Symbol allow-lists, positivity) run while
  constructing the typed objects; some legacy aliases are accepted there.
- `load_sparlectra_config(...)` caches the typed result for unchanged files
  without overrides; file hash/mtime changes invalidate the cache.

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

`power_flow.qlimits.enforcement_mode` selects the reactive-limit algorithm:
`active_set` is the in-iteration PV to PQ active-set switching;
`classic_simultaneous` and `classic_one_at_a_time` run a classical outer
loop (solve with switching disabled, clamp violating generator Q to the
limit, convert those buses to PQ, rerun without PQ to PV re-enable). The
legacy aliases `matpower_simultaneous` and `matpower_one_at_a_time` are
normalized to the `classic_*` value.

`model.auto_profile` controls the MATPOWER pre-run profile: `off`,
`recommend` (log a recommendation table, change nothing), `apply` (apply
only unambiguous import-convention recommendations). The runner logs the
effective MATPOWER options and never rewrites user YAML files.

## [Wrong-branch detection semantics (rectangular PF)](@id config-wrong-branch)

`power_flow.wrong_branch_detection` is a post-convergence plausibility check
for rectangular PF results; it is heuristic and does not prove global branch
correctness.

| Mode | Effect |
|---|---|
| `off` | checker disabled |
| `warn` | suspicious solutions are reported in rectangular status metadata; numerical convergence remains accepted |
| `fail` | suspicious solutions are treated as failed final convergence |
| `rescue` | reserved mode: a suspicious solution reports `rescue_requested_but_not_available`; no retry loop runs |

The thresholds cover the voltage magnitude range, the global angle spread,
and active-branch angle differences via
`power_flow.wrong_branch_max_branch_angle_deg`. All heuristics are evaluated
only on the network's highest nominal voltage level (branch-angle checks:
both ends on that level); the reported `min_vm_pu`/`max_vm_pu` and the
lowest-bus list refer to the checked level.

!!! details "Why only the highest voltage level is judged"
    Sub-transmission levels routinely run at 0.94 to 0.97 pu in healthy
    snapshots, and judging them against the transmission-level band
    produces false `SUSPECT` verdicts.

### Where the result is visible

| Surface | Fields |
|---|---|
| `ACPFlowReport.metadata` | `wrong_branch_status` and `wrong_branch_reason`. |
| AC island diagnostics CSV (`ac_island_solver_summary.csv`, one row per island) | trailing `wrong_branch_status`/`wrong_branch_reason` columns next to the `wrong_branch_detection` *setting* column; the per-island `ac_island_<id>_solver.log` lists both fields. |
| Console/log summary (`printACPFlowResults`) | a `Wrong-branch   : SUSPECT (...)` or `Wrong-branch   : FAIL (...)` line, printed only when the result is neither `ok` nor `not_checked`. |
| Web UI run result page | a "Wrong-branch check" badge with the run-status styling; omitted when the result is `not_checked`. |
| `run_sparlectra_api` result metadata | `wrong_branch_status`, `wrong_branch_reason`, `wrong_branch_low_vm_count`, `wrong_branch_high_vm_count`, `wrong_branch_angle_spread_deg`, `wrong_branch_branch_angle_violation_count`. |

| `status` value | Meaning |
|---|---|
| `ok` | checked, no finding |
| `warn` | suspicious, result accepted under `wrong_branch_detection = warn` |
| `fail` | suspicious, treated as non-convergence under `fail` |
| `wrong_branch_rescue_not_implemented` | the reserved `rescue` mode was requested |
| `not_checked` | `wrong_branch_detection = off`, or the check never ran (e.g. a non-finite solution) |

`reason` values: `none`, `low_voltage_magnitude`, `high_voltage_magnitude`,
`angle_spread_exceeded`, `branch_angle_exceeded`, `nonfinite_voltage`,
`disabled`, `rescue_requested_but_not_available`.

The wrong-branch retry loop (`wrong_branch_detection = rescue`) is a
reserved mode and not implemented; it reports
`rescue_requested_but_not_available`. It is distinct from the general
rescue ladder for failed AC solves (`power_flow.rescue`, see
[Power-Flow Configuration](powerflow_configuration.md#pf-solver-core)).
For hard flat-start cases the mitigations are this detection and the APSLF
solver as an alternative start/solve path (see
[External Solvers](external_solvers.md)).

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
  verbose_passes: false
  controllers: {}
```

| Key | Type | Default | Meaning |
|---|---:|---:|---|
| `control.enabled` | Bool | `true` | Run the outer control loop (`run_control!`) around the inner solver when controllers exist. |
| `control.max_outer_iterations` | Int | `20` | Outer-loop budget shared by all active controllers. Does not control inner NR iterations. |
| `control.trace` | Bool | `true` | Collect machine-readable rows in `ControlRunResult.trace`. |
| `control.log_iterations` | Bool | `true` | One console line per control pass (converged, mismatch, active-set changes) when `output.console_diagnostics` is `full`. |
| `control.stop_on_pf_failure` | Bool | `true` | Abort the orchestration when the inner power flow fails. |
| `control.verbose_passes` | Bool | `false` | Repeat the full inner-solver diagnostic blocks on every control pass. Off: the first pass prints them once, later passes get the one-line summary. |
| `control.controllers` | Mapping | `{}` | Declarative controller instantiation: one named mapping per controller, applied to the net before the outer loop. See [Control Framework](control_framework.md) for the schema. |

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

| `type` | Device function | Notes |
|---|---|---|
| `power_transformer` | `addPowerTransformerControl!` | |
| `machine_voltage` | `addMachineVoltageControl!` | |
| `shunt_voltage` | `addShuntVoltageControl!` | |
| `series_reactance` | `addSeriesReactanceControl!` | |
| `hvdc_pair` | `addHvdcPairControl!` | |
| `upfc` | `addUpfcControl!` | `model: quadrature` registers the SSSC+STATCOM pair, `model: full` the DC-link-coupled independent-P/Q model, see [FACTS Devices](facts.md) |

Unknown types or keys and missing required keys fail at configuration load;
unknown bus/branch/transformer references and invalid limits fail at apply
time naming the entry. An entry whose element already carries a controller
of the same type is skipped. `applyConfiguredControllers!` applies a
configuration to a programmatically built net.
`examples/others/tap_control_demo_grid.yaml` (setpoints and tap/phase
parameters `oltc`, `pst`, `schraeg` for
`examples/others/tap_control_demo_grid.jl`) is an example-specific input
file, not the `control.controllers` schema.

## Bookkeeping and console keys

| Key | Default | Meaning |
|---|---|---|
| `runtime.casefile`, `runtime.case_name`, `runtime.case_source`, `runtime.configured_default_casefile` | `""` | Populated by the API/Web UI into the `effective_config.yaml` run artifact; they record which case a run actually used and where it came from. Not user inputs. |
| `output.console_live` | `false` | Mirror captured run output live to the real console during API/service runs; the archived `run.log` stays identical. |
| `output.result_table_max_rows` | `200` | Row cap for the classical result tables. |
| `output.result_table_large_case_threshold_buses` | `1000` | Bus count from which a case counts as large for result rendering. |
| `output.result_table_large_case_mode` | `summary` | What large cases print instead of full tables (`summary`, `classic`, `full`). |
| `output.csv_format` | `technical` | Delimiter/decimal-separator format of every CSV file a run writes (`write_result_csv`): `bus_voltages_complex.csv`, `branch_flows.csv`, `bus_powers.csv`, `q_limit_*.csv`, the short-circuit, contingency and scenario tables, the SE diagnostic exports and `se_state.csv`, `ac_islands.csv`, the SV comparison and the DTF outage metrics. Only the measurement CSV that Sparlectra reads back keeps its fixed layout. Allowed values: `technical` (comma delimiter, dot decimal), `excel_de` (semicolon delimiter, comma decimal, dot thousands separator), `excel_us` (comma delimiter, dot decimal, comma thousands separator). The API's `detailed_result_csv_format`/`detailed_result_csv_semicolon` request keywords are a deprecated per-request override of this key. |
| `webui.operation_log_retention_days` | `10` | Int, `>= 0`. How far the operation log reaches back. Every Web UI start drops older entries from every operation log it knows; `0` keeps only the current session. Lower it when the log page grows unwieldy: its size comes from the number of entries, not from their age. The environment variable `SPARLECTRA_WEBUI_OPERATION_LOG_RETENTION_DAYS` still wins, for headless runs that read no configuration file. |

| `webui.docs_base_url` | `https://welthulk.github.io/Sparlectra.jl/` | String. Base URL of the published documentation the help pages (the **?** next to a control) and the header link open; set it to a local docs build or a pinned version folder. A link only, nothing is fetched; the help pages themselves ship with the application. |

## Migration notes

| Legacy / old key | Canonical key | Notes |
|---|---|---|
| `matpower_import.benchmark` | `benchmark.enabled` | Moved to the top-level `benchmark` section. |
| `methods` (top-level legacy path) | `benchmark.methods` | Keep benchmark methods in `benchmark`. |
| `max_ite` | `power_flow.max_iter` | Legacy alias; prefer the canonical nested key. |

## Detailed references

- [Power-Flow Configuration](powerflow_configuration.md)
- [MATPOWER cases](matpower.md)
- [State-Estimation Configuration](state_estimation_configuration.md)
- [Performance and Profiling Configuration](performance_profiling.md)

The canonical key set lives in one place,
[`src/config/configuration.yaml.example`](https://github.com/Welthulk/Sparlectra.jl/blob/main/src/config/configuration.yaml.example),
commented per key and always current. The effective configuration of a run,
with every default and override resolved, prints via
`print_effective_config`.

# MATPOWER Import Configuration

## MATPOWER reference terminology

| Term | Meaning in Sparlectra | Is it a solver start? | Is it a comparison reference? |
|---|---|---:|---:|
| `BUS.VM` / `BUS.VA` | Imported MATPOWER bus voltage columns | sometimes | yes |
| `GEN.VG` | Generator PV voltage setpoint column | yes (PV setpoint source) | yes |
| imported setpoint | Voltage setpoint chosen by import logic | sometimes | yes |
| historical value | Prior Sparlectra/SCADA/SE state | yes | no (unless explicitly configured elsewhere) |

## Option reference

| YAML path | Type | Default | Allowed values | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `matpower_import.case` | String | `case14.m` | case path/name | Compatible single-case selector and fallback when `cases` is empty. | Single case studies and benchmarks. | Invalid/missing paths. | Parse/solve scales with case size. | `matpower_import.cases`, runtime/profile. |
| `matpower_import.cases` | Vector{String} | `[case14.m, case118.m]` | non-empty case names | Ordered batch selector for `run_sparlectra_cases`; a non-empty list takes precedence over `case`. | Deterministic multi-case validation and release checks. | Empty case names or expecting `run_sparlectra` to return a vector. | Sequential parse/solve cost per case. | `matpower_import.case`, `run_sparlectra_cases`. |
| `matpower_import.auto_profile` | Symbol/String | `recommend` | `off`, `recommend`, `apply` | Run a MATPOWER pre-run profile. `off` disables it, `recommend` logs decisions without changing the active config, and `apply` changes only safe import-convention recommendations with clear evidence. Solver-start and Q-limit recommendations remain logged but skipped unless configured directly. | Development, large-case investigation, reproducible robust imports. | Expecting YAML files to be rewritten; applying ambiguous diagnostics. | Low; scans existing VM/VA residuals before the solve. | Output profile visibility options. |
| `matpower_import.auto_profile_log` | Bool | `true` | `true`, `false` | Print/log auto-profile reasoning and final effective options. | Debug import decisions and reproduce final settings. | Quiet high-volume runs. | Logging overhead only. | `output.console_auto_profile`, logfile settings. |
| `matpower_import.pv_voltage_source` | Symbol/String | `gen_vg` | `gen_vg`, `bus_vm`, `auto`, `strict_check` | PV voltage setpoint source policy. | Standard MATPOWER semantics. | Nonstandard conversion assumptions. | None. | `compare_voltage_reference`, PF starts. |
| `matpower_import.pv_voltage_mismatch_tol_pu` | Float64 | `1e-4` | nonnegative real | Tolerance for PV voltage mismatch checks. | Tight validation studies. | Overly strict noisy data. | Low. | `compare_voltage_reference`. |
| `matpower_import.compare_voltage_reference` | Symbol/String | `imported_setpoint` | `bus_vm`, `gen_vg`, `imported_setpoint`, `hybrid` | Voltage reference used for comparisons. Auto-profile recommends `hybrid` when BUS.VM / GEN.VG mismatches are detected. | MATPOWER comparison workflows. | When historical/SCADA ref should dominate. | Low. | `pv_voltage_source`, diagnostics. |
| `matpower_import.bus_shunt_model` | Symbol/String | `admittance` | `admittance`, `voltage_dependent_injection` | Bus shunt interpretation model. | Default import path. | Alternative modeling studies without residual evidence. | None. | Import convention diagnostics. |
| `matpower_import.shift_unit` | Symbol/String | `deg` | `deg`, `rad` | Phase-shift input unit. | Cases with radians metadata. | Wrong unit declaration. | Negligible. | `shift_sign`, branch shift diagnostics. |
| `matpower_import.shift_sign` | Float64 | `1.0` | real (typ. `1.0`, `-1.0`) | Phase-shift sign convention. | Cross-tool convention alignment. | Unnecessary flipping. | None. | `shift_unit`, branch-shift diagnostics. |
| `matpower_import.ratio` | Symbol/String | `normal` | `normal`, `reciprocal` | Branch ratio interpretation mode. | Standard MATPOWER import. | Unsupported alternate conventions. | None. | Transformer/tap interpretation. |
| `matpower_import.enable_pq_gen_controllers` | Bool | `true` | `true`, `false` | Enable controller behavior on imported PQ generators. | Realistic controlled studies. | Raw imported behavior reproduction. | Small control bookkeeping cost. | PF Q-limit behavior. |
| `matpower_import.preallocate_network` | Symbol/String | `auto` | `off`, `on`, `auto` | Controls import-time `sizehint!` preallocation for large MATPOWER network construction. | Large imports where construction allocations dominate runtime. | Tiny cases where tuning is unnecessary. | Can reduce import allocations/time; no model changes. | `matpower_import.preallocate_min_buses`. |
| `matpower_import.preallocate_min_buses` | Int | `1000` | positive integer | Bus-count threshold used when `preallocate_network = auto`. | Auto-tuning preallocation trigger for site-specific case sizes. | If fixed always-on/off behavior is preferred. | Threshold only; no model changes. | `matpower_import.preallocate_network`. |
| `matpower_import.apply_bus_names` | Bool | `false` | `true`, `false` | Use standard `mpc.bus_name` metadata for imported bus names. | FOR001/FOR002 validation and named-bus workflows. | Preserve historical numeric names. | None. | Fails on duplicate names. |
| `matpower_import.apply_branch_names` | Bool | `false` | `true`, `false` | Attach user-defined `mpc.branch_name` metadata to `net.matpower_branch_metadata`. | Outage and contingency mapping. | Cases without branch metadata. | None. | `import_for001_contingencies`. |
| `matpower_import.apply_branch_kind` | Bool | `false` | `true`, `false` | Use user-defined `mpc.branch_kind` to override line/transformer classification. | Conversion workflows that know row kinds. | Prefer electrical heuristic. | None. | Accepts `L`/`LINE`/`ACL` and `T`/`TRAFO`/`TRANSFORMER`/`2WT`. |
| `matpower_import.import_for001_contingencies` | Bool | `true` | `true`, `false` | Preserve user-defined `mpc.for001_contingencies`. | FOR001/FOR002 validation. | Ignore validation metadata. | None. | `mpc.branch_name` enables index mapping. |
| `matpower_import.matpower_dcline_mode` | Symbol/String | `reject_active` | `reject_active`, `ignore_inactive`, `pf_injections` | Controls active `mpc.dcline` rows. | Use `pf_injections` to emulate MATPOWER simple PF DC-line injections. | OPF/dclinecost studies. | Adds two fixed prosumers per active row in `pf_injections`. | Default keeps fail-fast behavior. |

## Auto-profile pre-run

`matpower_import.auto_profile` is evaluated by the MATPOWER runner before the
main solve. The pre-run reads the MATPOWER case, computes compact diagnostics,
and emits a table with the option path, current value, recommended value,
action, reason, and diagnostic evidence. The shipped default is
`auto_profile: recommend` with `auto_profile_log: true`, so normal MATPOWER
runner usage logs recommendation-only diagnostics by default without changing
the active configuration.

Modes:

- `off`: no pre-run is performed.
- `recommend`: diagnostics are logged, but the active run configuration is not
  changed.
- `apply`: Sparlectra applies only clearly safe convention changes. The first
  implementation may apply `shift_unit`, `shift_sign`, `ratio`, and
  `bus_shunt_model` when residual scans show a large, unambiguous improvement.
  PV/REF voltage, robust-start, and Q-limit guard settings are recommendation-only.

Diagnostics used by the pre-run:

- MATPOWER VM/VA power-balance residual scans over branch shift unit/sign and
  transformer ratio conventions.
- Bus-shunt residual comparison with and without MATPOWER bus shunt admittance.
- PV/REF `BUS.VM` versus online `GEN.VG` mismatch counts.
- Case-size, PV-bus count, and generator Q-range heuristics for robust-start and
  Q-limit recommendations.

Explicit YAML values remain visible. In `apply` mode, a changed option is shown
as `applied` in the auto-profile table and the final effective options block;
ambiguous diagnostics and conservative solver-start or Q-limit recommendations
are shown as `keep` or `skipped` and do not change the active run. Auto-profile
never rewrites user YAML files. To reproduce a run, copy the logged final
effective options (or enable
`diagnostics.log_effective_config`) into a tracked configuration file.

### Comparing multiple MATPOWER configurations

`examples/matpower_import_multi_config.jl` is a developer-oriented helper for
running the same MATPOWER case against several YAML configuration files. It is
useful for checking whether options such as `matpower_import.auto_profile`,
`power_flow.wrong_branch_detection`, or start-mode settings affect the final
rectangular solver status.

Example:

```bash
julia --project=. examples/matpower_import_multi_config.jl \
  data/mpower/case14.m \
  --config=path/to/config_a.yaml \
  --config=path/to/config_b.yaml \
  --status-only
```

The script does not create or rewrite YAML files. Repeated `--config=...`
arguments and comma- or semicolon-separated `--configs=A,B,C` lists are both
supported. In `--status-only` mode it prints the rectangular status and
wrong-branch diagnostic fields for each configuration; `--runner` delegates to
the standard `Sparlectra.run_matpower_case` runner output.

## Citation and case-file usage

Sparlectra references MATPOWER case names for diagnostics and comparison workflows,
but does **not** redistribute MATPOWER case files.

If you use MATPOWER software, data formats, or case files in your workflow, please
cite MATPOWER as recommended by the official guidance page:

- <https://matpower.org/citing/>

When Sparlectra documentation mentions case names (for example `case300.m`,
`case1354pegase.m`, `case1951rte.m`, or `case_ACTIVSg10k.m`), users should obtain
the original files from MATPOWER and/or the original data sources and follow the
applicable license, citation, and redistribution terms.

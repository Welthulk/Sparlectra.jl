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
| `matpower_import.case` | String | `case14.m` | case path/name | MATPOWER case selector. | Case studies and benchmarks. | Invalid/missing paths. | Parse/solve scales with case size. | `benchmark.*`, runtime/profile. |
| `matpower_import.auto_profile` | Symbol/String | `recommend` | `off`, `recommend`, `apply` | Run a MATPOWER pre-run profile. `off` disables it, `recommend` logs decisions without changing the active config, and `apply` changes only safe convention recommendations with clear residual evidence. | Development, large-case investigation, reproducible robust imports. | Expecting YAML files to be rewritten; applying ambiguous diagnostics. | Low; scans existing VM/VA residuals before the solve. | Output profile visibility options. |
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

## Auto-profile pre-run

`matpower_import.auto_profile` is evaluated by the MATPOWER runner before the
main solve. The pre-run reads the MATPOWER case, computes compact diagnostics,
and emits a table with the option path, current value, recommended value,
action, reason, and diagnostic evidence.

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
ambiguous diagnostics are shown as `keep` or `recommend` and do not change the
active run. Auto-profile never rewrites user YAML files. To reproduce a run,
copy the logged final effective options (or enable
`diagnostics.log_effective_config`) into a tracked configuration file.

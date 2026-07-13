# MATPOWER format

Sparlectra can read power-flow test cases in the MATPOWER case format. A MATPOWER case is usually a `.m` file that defines or returns a MATLAB/Octave structure named `mpc`. In the current case format, this structure contains `baseMVA`, `bus`, `gen`, `branch` and optionally `gencost`.[^matpower-data-format]

For a normal power-flow run in Sparlectra, the most relevant parts are `baseMVA`, `bus`, `gen` and `branch`. The optional `gencost` table is mainly relevant for optimal power flow and is usually not needed for a plain load-flow calculation.[^matpower-caseformat]

## Minimal structure

A typical case file has this shape:

```matlab
function mpc = case_example
mpc.version = '2';
mpc.baseMVA = 100;

mpc.bus = [
    % bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
      1     3   0  0  0  0  1    1  0  110    1    1.1  0.9;
      2     1  50 30  0  0  1    1  0  110    1    1.1  0.9;
];

mpc.gen = [
    % bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin
      1   80 0  100  -100 1.0 100   1      100  0;
];

mpc.branch = [
    % fbus tbus r     x      b      rateA rateB rateC ratio angle status angmin angmax
      1    2    0.02  0.06   0.03   100   100   100   0     0     1      -360   360;
];
```

The exact number and meaning of columns is fixed by MATPOWER. The official MATPOWER documentation is the authoritative reference for the complete column lists.[^matpower-data-format]

## `baseMVA`

`baseMVA` is the common power base of the case. MATPOWER uses it for per-unit quantities. Typical examples use values such as `100`, but the correct value depends on the case data.

## `bus`: network nodes

`mpc.bus` describes the electrical buses. Each row is one bus. The most important columns for users are:

| Column | Meaning |
|---|---|
| `BUS_I` | Bus number. Branches and generators refer to this number. |
| `BUS_TYPE` | Bus type: `1` = PQ, `2` = PV, `3` = reference/slack, `4` = isolated. |
| `PD`, `QD` | Active and reactive demand in MW and MVAr. |
| `VM`, `VA` | Voltage magnitude in p.u. and voltage angle in degrees. Often used as initial or reference values. |
| `BASE_KV` | Nominal voltage level in kV. |
| `VMAX`, `VMIN` | Upper and lower voltage magnitude limits in p.u. |

A PQ bus normally represents a load bus. A PV bus usually has generation with a specified active power and voltage magnitude. The reference bus supplies the angular reference for the island.

## `gen`: generators and voltage setpoints

`mpc.gen` describes generators or equivalent injections. Each row is one generator connected to a bus.

Important columns are:

| Column | Meaning |
|---|---|
| `GEN_BUS` | Bus number where the generator is connected. |
| `PG`, `QG` | Active and reactive generator output in MW and MVAr. |
| `QMAX`, `QMIN` | Reactive power limits. These are relevant for Q-limit handling. |
| `VG` | Voltage magnitude setpoint in p.u. |
| `GEN_STATUS` | `> 0` means in service, `<= 0` means out of service. |
| `PMAX`, `PMIN` | Active power limits in MW. |

For load-flow studies, `PG`, `VG`, `QMAX` and `QMIN` are usually the values users inspect first.

## `branch`: lines and transformers

`mpc.branch` describes network connections. Each row is one line, cable, transformer or equivalent branch.

Important columns are:

| Column | Meaning |
|---|---|
| `F_BUS`, `T_BUS` | From-bus and to-bus numbers. |
| `BR_R`, `BR_X` | Series resistance and reactance in p.u. |
| `BR_B` | Total line charging susceptance in p.u. |
| `RATE_A`, `RATE_B`, `RATE_C` | Thermal ratings. `0` is commonly used for unlimited. |
| `TAP` | Transformer off-nominal tap ratio. `0` indicates no off-nominal transformer tap and is mathematically treated like a ratio of `1`. |
| `SHIFT` | Transformer phase-shift angle in degrees. |
| `BR_STATUS` | `1` means in service, `0` means out of service. |
| `ANGMIN`, `ANGMAX` | Optional voltage angle difference limits. |

In MATPOWER, lines and transformers are represented in the same `branch` matrix. A non-zero tap ratio or phase shift indicates transformer behavior.

## `gencost`: optional cost data

`mpc.gencost` contains generator cost data for optimal power-flow studies. It can be present in many public MATPOWER cases, but it is usually not needed for a plain power-flow run.

## Practical checks before importing a case

Before running a case in Sparlectra, check the following points:

1. The file should define or return `mpc`.
2. `mpc.version` should normally be `'2'`.
3. Bus numbers in `branch` and `gen` must exist in `mpc.bus`.
4. At least one in-service reference bus is required per connected island.
5. Decimal values must use a decimal point, not a decimal comma.
6. Inactive generators or branches should have status `0`.
7. Per-unit branch impedances must be consistent with `baseMVA` and the voltage base of the case.

## What Sparlectra uses

For a standard Web UI power-flow run, Sparlectra primarily uses the topology and electrical data from `bus`, `gen` and `branch`. The result page and log files then report whether the numerical power flow converged, how many iterations were needed and what final mismatch remained.

When Q-limit handling is enabled, the generator reactive power limits from `QMAX` and `QMIN` become important. If a generator cannot maintain its voltage setpoint within those limits, the solver may have to treat the bus differently depending on the selected configuration.

## Sparlectra-specific optional metadata and import options

The standard MATPOWER case concepts are `mpc.bus`, `mpc.gen`, `mpc.branch` and, when present, `mpc.gencost`. The official MATPOWER documentation remains authoritative for those standard tables and their columns.

Sparlectra can also read a small set of optional metadata fields when a conversion or validation workflow provides them. These fields are Sparlectra-recognized extensions, not required MATPOWER fields. If they are absent, standard MATPOWER imports still work. If they are present but the corresponding `apply_*` import options remain `false`, Sparlectra keeps conservative/default naming and branch-classification behavior.

Names and source metadata are useful for CSV exports, FOR001/FOR002 validation, outage mapping, diagnostics and support workflows.

| Field | Purpose | Used when |
|---|---|---|
| `mpc.bus_name` | Optional readable bus names. | `matpower_import.apply_bus_names = true` |
| `mpc.branch_name` | Optional readable branch/source names. | `matpower_import.apply_branch_names = true` |
| `mpc.branch_kind` | Optional branch-kind hints such as line/transformer. | `matpower_import.apply_branch_kind = true` |
| `mpc.for001_contingencies` | Optional FOR001-derived outage/contingency metadata for validation workflows. | `matpower_import.import_for001_contingencies = true` |
| `mpc.dcline` | MATPOWER DC-line data. | Controlled by `matpower_import.matpower_dcline_mode`; active rows are rejected by default unless `pf_injections` is selected. |
| `mpc.sparlectra` | Versioned Sparlectra extension namespace for metadata that standard MATPOWER does not model. | Automatically recognized when present. |

### Sparlectra transformer-loss extension

Standard MATPOWER has no branch-local field equivalent to the active no-load
conductance found on some FOR/DTF transformer records. Sparlectra therefore
exports transformer active-loss metadata under `mpc.sparlectra` when such data
exists. The block is versioned with `mpc.sparlectra.format_version = 1` and
stores `mpc.sparlectra.transformer_losses` entries with the MATPOWER branch row,
original FOR/DTF identifier, terminal bus names, raw `R/X/G/B` values where
available, converted per-unit values, tap/phase-shift data, source marker, and
the allocation convention.

Sparlectra applies FOR/DTF transformer conductance as part of the native
transformer PI branch model (`Branch.g_pu`) rather than as synthetic terminal
bus shunts. Exported `.m` files include an explicit `SPARLECTRA EXTENSION
WARNING` comment because plain MATPOWER ignores the proprietary metadata and may
therefore compute different transformer active losses. When Sparlectra reimports
one of these files it restores the branch conductance metadata without adding
terminal bus-shunt approximations.

## Sparlectra MATPOWER import options overview

The options below are the most common import-convention controls users may need when a case was converted from another tool or carries optional metadata. See the full [MATPOWER import configuration reference](matpower_import.md) for allowed values, defaults and interactions.

| Option | User-facing purpose |
|---|---|
| `matpower_import.auto_profile` | Runs a pre-run profile that can log or safely apply import-convention recommendations. |
| `matpower_import.auto_profile_log` | Controls whether auto-profile reasoning and final effective options are printed/logged. |
| `matpower_import.pv_voltage_source` | Selects whether PV voltage setpoints come from generator `VG`, bus `VM`, or an automatic/strict policy. |
| `matpower_import.compare_voltage_reference` | Chooses the voltage reference used for MATPOWER comparison diagnostics. |
| `matpower_import.bus_shunt_model` | Selects how MATPOWER bus shunts are interpreted during import. |
| `matpower_import.shift_unit` | Declares whether branch phase-shift values are in degrees or radians. |
| `matpower_import.shift_sign` | Controls the phase-shift sign convention used for branch import. |
| `matpower_import.ratio` | Controls the transformer tap-ratio interpretation. |
| `matpower_import.enable_pq_gen_controllers` | Enables controller behavior for imported PQ generators when appropriate. |
| `matpower_import.apply_bus_names` | Applies optional `mpc.bus_name` metadata to imported bus names. |
| `matpower_import.apply_branch_names` | Preserves optional `mpc.branch_name` metadata for branch/source mapping. |
| `matpower_import.apply_branch_kind` | Uses optional `mpc.branch_kind` metadata to classify lines and transformers. |
| `matpower_import.import_for001_contingencies` | Imports optional `mpc.for001_contingencies` validation metadata. |
| `matpower_import.matpower_dcline_mode` | Controls whether active `mpc.dcline` rows are rejected, ignored when inactive, or mapped with the MATPOWER DC-line PF injections approximation. |

For example, a case-conversion or validation workflow might use:

```yaml
matpower_import:
  auto_profile: recommend
  pv_voltage_source: gen_vg
  compare_voltage_reference: imported_setpoint
  ratio: normal
  shift_unit: deg
  shift_sign: 1.0
  bus_shunt_model: admittance
  apply_bus_names: true
  apply_branch_names: true
  apply_branch_kind: true
  import_for001_contingencies: true
  matpower_dcline_mode: reject_active
```

This is an example configuration, not a requirement for every case.

## MATPOWER `mpc.dcline`

Sparlectra recognizes `mpc.dcline` tables. The default
`matpower_import.matpower_dcline_mode: reject_active` remains fail-fast for
active rows so that DC lines are not silently dropped. Empty and inactive-only
tables are tolerated.

For ordinary power-flow compatibility studies, opt in to
`matpower_import.matpower_dcline_mode: pf_injections`. This imports active rows
with a MATPOWER `toggle_dcline`-compatible PF approximation: each DC line is
represented by from/to generator-like terminal injections, from-side
`PG = -PF`, to-side received power from `PF - (LOSS0 + LOSS1 * PF)` when loss
columns are available (or input `PT` otherwise), row `QF`/`QT`, voltage
setpoints `VF`/`VT`, and terminal Q limits where present. Terminal buses with
voltage setpoints become voltage-controlled where appropriate, while existing
reference buses stay reference buses and isolated buses are not activated.
API/Web UI runs with imported active DC lines include `matpower_dcline.csv` as
a compact diagnostic artifact.

This option is a `toggle_dcline`-compatible PF approximation, not a full HVDC
converter model. OPF-style behavior, converter controls, `dclinecost`, and
detailed DC-grid modeling are not provided by Sparlectra yet.

The Web UI exposes selected MATPOWER import convention controls in the form and writes the effective configuration as a run artifact. It does not rewrite the user YAML automatically. With `auto_profile = recommend`, Sparlectra logs recommendations without changing the active configuration; with `auto_profile = apply`, it applies only safe convention recommendations. Use the result artifacts and logs to reproduce the final effective settings for a run.

## Full reference

This page is intentionally short. For the complete MATPOWER format, use the official MATPOWER documentation:

- MATPOWER Manual, Appendix “Data File Format”: <https://matpower.app/manual/matpower/DataFileFormat.html>
- MATPOWER `caseformat` reference: <https://matpower.org/documentation/ref-manual/legacy/functions/caseformat.html>
- MATPOWER named column constants: <https://matpower.org/docs/ref/matpower6.0/define_constants.html>

[^matpower-data-format]: MATPOWER Manual, “Data File Format”. It describes version 2 case files as an `mpc` structure with `baseMVA`, `bus`, `branch`, `gen` and optional `gencost`; it also defines the bus, generator, branch and cost columns. <https://matpower.app/manual/matpower/DataFileFormat.html>
[^matpower-caseformat]: MATPOWER `caseformat` reference. It defines a MATPOWER case file as an M-file or MAT-file that returns or defines an `mpc` case structure and documents version 2, optional fields and solver-added result columns. <https://matpower.org/documentation/ref-manual/legacy/functions/caseformat.html>

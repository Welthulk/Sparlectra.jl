# MATPOWER cases

Sparlectra reads MATPOWER `.m` case files (format version 2) and writes
networks back as `.m` files, including the solved state. The `mpc`
structure (`baseMVA`, `bus`, `gen`, `branch`, optional `gencost`) and
every column are defined by the MATPOWER Manual, Appendix "Data File
Format": <https://matpower.app/manual/matpower/DataFileFormat.html>
(`caseformat` reference:
<https://matpower.org/documentation/ref-manual/legacy/functions/caseformat.html>,
named column constants:
<https://matpower.org/docs/ref/matpower6.0/define_constants.html>). A
power flow needs `baseMVA`, `bus`, `gen` and `branch`; `gencost` serves
optimal power flow, which Sparlectra does not have: the cost rows stay with
the generators and are written back on export
([Generator costs](@ref matpower_gencost)).

## Import

`run_sparlectra` imports a case and runs the configured framework workflow,
returning one `SparlectraRunResult`:

```julia
using Sparlectra

case_path = ensure_casefile("case14.m")

# Import a MATPOWER case file and run the configured framework workflow.
result = run_sparlectra(
    casefile = basename(case_path),
    path = dirname(case_path),
)

net = result.net
```

`ensure_casefile` downloads a standard MATPOWER test case on demand into
the user case cache and returns the absolute path; a `.jl` name
(`ensure_casefile("case14.jl")`) requests a generated `.jl` companion
file. Pass `path` for local files and `config` for a loaded configuration
(`load_sparlectra_config("examples/configuration.yaml")`).

`run_sparlectra_cases` runs several configured cases in order; a non-empty
`runtime.cases` list wins over `runtime.case`, and performance profiles
are per `run_sparlectra` call, not per batch:

```yaml
matpower_import:
  cases: [case14.m, case118.m]
```

```julia
using Sparlectra

cfg = load_sparlectra_config("path/to/configuration.yaml")
results = run_sparlectra_cases(config = cfg)

for result in results
    println(result.net.name, ": ", result.outcome)
end
```

`createNetFromMatPowerFile` imports without a power flow. The transformer
conventions come from `SparlectraConfig.matpower` in framework runs or are
passed directly (`matpower_ratio = "normal"`, the default, uses the branch
`TAP` values as stored, `"reciprocal"` their reciprocal):

```julia
using Sparlectra

# Import the network configuration
net = createNetFromMatPowerFile(filename = "path/to/case5.m", log = false)

# Explicit import conventions
net = createNetFromMatPowerFile(
    filename = case_path,
    matpower_shift_unit = "rad",
    matpower_shift_sign = -1.0,
    matpower_ratio = "reciprocal",
)
```

Run the power flow on an imported network directly:

```julia
using Sparlectra

net = createNetFromMatPowerFile(filename = "case5.m", log = false)

tol = 1e-6
max_ite = 10
verbose = 0
ite, erg = runpf!(net, max_ite, tol, verbose)

if erg == 0
  calcNetLosses!(net)
  printACPFlowResults(net, 0.0, ite, tol)
else
  @warn "Power flow did not converge after \$ite iterations"
end
```

`casefileparser` returns the raw data arrays without building a network:

```julia
using Sparlectra

case_name, baseMVA, busData, genData, branchData = casefileparser("case9.m")
println("Number of buses: \$(size(busData, 1))")
```

`examples/powerflow/matpower_import.jl` resolves `runtime.julia_threads`
after the CLI (`--julia-threads=8`) and environment
(`SPARLECTRA_JULIA_THREADS`) overrides and re-executes itself once with
that `--threads` setting when it differs from `Threads.nthreads()`:

```bash
julia --threads=8 --project=. examples/powerflow/matpower_import.jl
julia --project=. examples/powerflow/matpower_import.jl --julia-threads=8
```

The Web UI imports `.m` files copy-only through **Import case files** and
shows the import-convention controls in its form
([Import case files](@ref webui-import-case-files)).

## Import options

The comparison options use these MATPOWER reference terms:

| Term | Meaning in Sparlectra | Is it a solver start? | Is it a comparison reference? |
|---|---|---:|---:|
| `BUS.VM` / `BUS.VA` | Imported MATPOWER bus voltage columns | sometimes | yes |
| `GEN.VG` | Generator PV voltage setpoint column | yes (PV setpoint source) | yes |
| imported setpoint | Voltage setpoint chosen by import logic | sometimes | yes |
| historical value | Prior Sparlectra/SCADA/SE state | yes | no (unless explicitly configured elsewhere) |

## [Option reference](@id matpower-options)

The keys of a MATPOWER run: the case selectors, the import profile that recommends or applies the conventions a case needs, the conversion switches (phase-shift unit and sign, tap ratio, PV voltage source, shunt and tap-changer model, generator controllers), the solution export and the network preallocation.

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `runtime.case` | String | `case14.m` | case path/name | Single-case selector and fallback when `cases` is empty. |
| `runtime.cases` | Vector{String} | `[case14.m, case118.m]` | non-empty case names | Ordered batch selector for `run_sparlectra_cases`; a non-empty list takes precedence over `case`. |
| `model.auto_profile` | Symbol/String | `off` | `off`, `recommend`, `apply` | Experimental convention scan of the stored voltages ([MATPOWER conventions (experimental)](@ref matpower_conventions_experimental)). `off` disables it, `recommend` logs decisions without changing the active config, `apply` changes only `shift_unit`, `shift_sign`, `ratio` and `bus_shunt_model`, and only on clear evidence; PV/REF voltage, solver-start and Q-limit recommendations are logged, never applied. YAML files are never rewritten. |
| `model.auto_profile_log` | Bool | `true` | `true`, `false` | Print/log auto-profile reasoning and final effective options. |
| `model.auto_profile_max_fit_pu` | Float64 | `0.1` | non-negative | Largest power mismatch (pu, worst bus) at which the stored `VM`/`VA` columns still count as a solved state. When the best convention reading scores above it, the scan recommends no convention change and keeps the configured `shift_unit`, `shift_sign`, `ratio` and `bus_shunt_model`. Raise only for files with a slightly stale but consistent solution; lower to make the scan stricter. |
| `matpower_import.pv_voltage_source` | Symbol/String | `gen_vg` | `gen_vg`, `bus_vm`, `auto`, `strict_check` | PV voltage setpoint source; `gen_vg` is the standard MATPOWER semantics. |
| `matpower_import.pv_voltage_mismatch_tol_pu` | Float64 | `1e-4` | nonnegative real | Tolerance for the PV voltage mismatch checks. |
| `matpower_import.compare_voltage_reference` | Symbol/String | `imported_setpoint` | `bus_vm`, `gen_vg`, `imported_setpoint`, `hybrid` | Voltage reference used for comparisons. Auto-profile recommends `hybrid` when `BUS.VM`/`GEN.VG` mismatches are detected. |
| `model.bus_shunt_model` | Symbol/String | `admittance` | `admittance`, `voltage_dependent_injection` | Bus shunt interpretation ([Branch model](branchmodel.md)); change it only with residual evidence. |
| `matpower_import.shift_unit` | Symbol/String | `deg` | `deg`, `rad` | Phase-shift input unit. |
| `matpower_import.shift_sign` | Float64 | `1.0` | real (typ. `1.0`, `-1.0`) | Phase-shift sign convention. |
| `matpower_import.ratio` | Symbol/String | `normal` | `normal`, `reciprocal` | Branch ratio interpretation: `TAP` as stored, or its reciprocal. |
| `model.tap_changer_model` | Symbol/String | `ideal` | `ideal`, `impedance_correction` | Tap-changer model applied to all transformers after import, shared with the DTF importer; see [Transformer tap-changer model](configuration.md#transformer-tap-changer-model). |
| `matpower_export.write_solution` | Bool | `true` | `true`, `false` | Whether [`writeMatpowerCasefile`](@ref) writes the solved state (see [Export](#Export)). |
| `matpower_import.enable_pq_gen_controllers` | Bool | `true` | `true`, `false` | Enable controller behavior on imported PQ generators; `false` reproduces the raw imported behavior. |
| `model.preallocate_network` | Symbol/String | `auto` | `off`, `on`, `auto` | Import-time `sizehint!` preallocation for large network construction; no model change. |
| `model.preallocate_min_buses` | Int | `500` | positive integer | Bus-count threshold used when `preallocate_network = auto`. |
| `matpower_import.apply_bus_names` | Bool | `false` | `true`, `false` | Use `mpc.bus_name` for the imported bus names; fails on duplicate names. |
| `matpower_import.apply_branch_names` | Bool | `false` | `true`, `false` | Attach `mpc.branch_name` to `net.matpower_branch_metadata` (outage and contingency mapping). |
| `matpower_import.apply_branch_kind` | Bool | `false` | `true`, `false` | Use `mpc.branch_kind` to override the line/transformer classification; accepts `L`/`LINE`/`ACL` and `T`/`TRAFO`/`TRANSFORMER`/`2WT`. |
| `matpower_import.import_for001_contingencies` | Bool | `true` | `true`, `false` | Preserve `mpc.for001_contingencies` (FOR001/FOR002 validation); `mpc.branch_name` enables index mapping. |
| `matpower_import.matpower_dcline_mode` | Symbol/String | `pf_injections` | `reject_active` (deprecated), `ignore_inactive`, `pf_injections`, `paired_control` | Handling of active `mpc.dcline` rows, see [DC lines and disconnected AC islands](#DC-lines-and-disconnected-AC-islands). A configuration file carrying the deprecated `reject_active` loads as `pf_injections` with a warning; the strict fail-fast check stays available programmatically: `createNetFromMatPowerCase(matpower_dcline_mode = :reject_active)`. |

Example configuration for a case-conversion or validation workflow:

```yaml
matpower_import:
  auto_profile: off
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
  matpower_dcline_mode: pf_injections
```

### [MATPOWER conventions (experimental)](@id matpower_conventions_experimental)

The convention overrides (`matpower_import.shift_unit`, `shift_sign`,
`ratio`, `model.bus_shunt_model`) and the convention scan
(`model.auto_profile`) are experimental, off by default and silent. The
standard MATPOWER reading (shift in degrees, sign +1, ratio as stored, bus
shunt as admittance) satisfies the power balance of the file on every case
checked; the alternative readings fail wherever a file has a phase shift
or an off-nominal tap. Change these only for a file you know to deviate.
At the standard values nothing about conventions is written to the run
log, the solver log, the result metadata or the configuration report, and
no `matpower_auto_profile.log` artifact is created; with the scan on, or a
convention set away from the standard reading, the lines and the artifact
appear, prefixed `experimental`, and the effective values appear in the
configuration report. In the Web UI the block sits under Advanced on the
Case page as "MATPOWER conventions (experimental)" and opens by itself
when a case configuration carries a non-standard value.
[Start strategies by case](start_strategies.md) explains why the stored
voltage columns of some files must not be used to change the reading.

### Auto-profile pre-run

With `model.auto_profile` on, the MATPOWER runner reads the case before
the main solve and prints a table with option path, current value,
recommended value, action, reason and evidence. Diagnostics: `VM`/`VA`
power-balance residual scans over the shift unit/sign and ratio
conventions, bus-shunt residuals with and without the bus shunt
admittance, PV/REF `BUS.VM` versus online `GEN.VG` mismatch counts, and
case size, PV-bus count and generator Q-range heuristics for the
robust-start and Q-limit recommendations. For large cases the pre-run
keeps start projection disabled, uses DC-angle and blended-voltage
flat-start seeds, recommends a practical validation tolerance, and
disables expensive diagnostics unless requested.

The convention ranking is evidence only when the stored `VM`/`VA` columns
are a solved state under some reading: when the best of the eight
readings still scores above `model.auto_profile_max_fit_pu` (default
`0.1` pu at the worst bus), the scan recommends no convention change,
`apply` changes nothing, and the table and the console line say so with
the best score ("stored VM/VA are not a solved state under any reading").
The PEGASE files ship such columns (best fit 0.5 to 3.9 pu); a solved
state scores below 0.1 pu. In `apply` mode a changed option is shown as
`applied` in the table and in the final effective options block; to
reproduce a run, copy the logged final effective options (or enable
`diagnostics.log_effective_config`) into a tracked configuration file.

`examples/powerflow/matpower_import_multi_config.jl` runs one case against
several YAML files (repeated `--config=...`, or `--configs=A,B,C`) to
check whether `model.auto_profile`, `power_flow.wrong_branch_detection` or
start-mode settings change the final solver status; `--status-only` prints
the status and wrong-branch fields per configuration, `--runner` delegates
to `Sparlectra.run_matpower_case`:

```bash
julia --project=. examples/powerflow/matpower_import_multi_config.jl \
  data/mpower/case14.m \
  --config=path/to/config_a.yaml \
  --config=path/to/config_b.yaml \
  --status-only
```

### DC lines and disconnected AC islands

`matpower_import.matpower_dcline_mode = pf_injections` (the default)
imports each active `mpc.dcline` row as MATPOWER's `toggle_dcline` power
flow does: two generator-like terminal prosumers, from-side `PG = -PF`,
to-side received power `PF - (LOSS0 + LOSS1 * PF)` when loss columns exist
(else the input `PT`), row `QF`/`QT`, voltage setpoints `VF`/`VT`, and
terminal Q limits where present. Terminal buses with voltage setpoints
become voltage-controlled where MATPOWER would make them PV; reference
buses stay reference buses and isolated buses are not activated. API and
Web UI runs write `matpower_dcline.csv` describing the mapping.

- `reject_active` (deprecated as a configuration value): an active row
  (`status != 0`) aborts before solving with
  `failure_reason = unsupported_matpower_dcline`.
- `ignore_inactive`: the same active-row rejection. Empty or
  inactive-only tables are tolerated in every mode.
- `paired_control`: one steerable `HvdcPairControl` per row, see
  [HVDC Back-to-Back](hvdc_back_to_back.md).

This is a power-flow approximation, not an HVDC converter or DC-grid
model: OPF constraints, converter controls, `dclinecost` and DC-line
optimization are unsupported. No AC branch or dummy admittance joins the
terminals, so a case such as `case_SyntheticUSA.m` consists of several AC
islands coupled only through DC-line injections; each island is solved on
its own and `ac_islands.csv` records the result per island
([Power-Flow Configuration](powerflow_configuration.md),
[Parallel Execution](parallel_execution.md)).

### Tap-impedance correction and reimport

A case built with `model.tap_changer_model = impedance_correction` exports
`BR_R`/`BR_X` values that already carry the correction
([Transformer tap-changer model](configuration.md#transformer-tap-changer-model));
`writeMatpowerCasefile` marks this with
`mpc.sparlectra.tap_changer_model = 'impedance_correction'`. On reimport,
`createNetFromMatPowerFile`/`createNetFromMatPowerCase` detect the marker
and skip `calcTapCorrectedRX` whatever `model.tap_changer_model` says, so
the correction is never applied twice. Cases without the marker, including
third-party MATPOWER cases, import unchanged.

## What Sparlectra reads

A power-flow run uses `bus`, `gen` and `branch`; with Q-limit handling
enabled, `QMAX` and `QMIN` decide whether a generator can hold its voltage
setpoint. MATLAB functions inside the `.m` file are not supported. The
fields below the standard matrices are Sparlectra extensions: when absent,
standard imports work unchanged; when present with the `apply_*` option at
`false`, naming and branch classification keep their defaults. The direct
importer takes the same switches as keywords (`apply_bus_names = true`,
`apply_branch_names = true`, `apply_branch_kind = true`,
`import_for001_contingencies = true`).

| Field | Used for | Notes |
|---|---|---|
| `mpc.baseMVA` | Common power base for per-unit quantities | Typically `100`. |
| `mpc.bus` | Network nodes: type (`1` PQ, `2` PV, `3` reference, `4` isolated), demand, shunts, `VM`/`VA` start or reference values, base voltage, limits | Lines and transformers refer to `BUS_I`. |
| `mpc.gen` | Generators: `PG`, `VG` setpoint, `QMAX`/`QMIN` for Q-limit handling, status | `PG`, `VG`, `QMAX` and `QMIN` are the values to inspect first. |
| `mpc.branch` | Lines and transformers: `BR_R`/`BR_X`, charging, ratings, `TAP`, `SHIFT`, status | A non-zero tap ratio or phase shift marks a transformer; `TAP = 0` is treated as ratio `1`. |
| `mpc.gencost` | Optimal power flow | Kept per generator ([`GenCost`](@ref)) and written back unchanged on export; no calculation reads it. |
| `mpc.bus_name` | Readable bus names | `matpower_import.apply_bus_names: true`; the original `BUS_I` values stay in `busOrigIdxDict`. |
| `mpc.branch_name` | Readable branch/source names, outage mapping | `matpower_import.apply_branch_names = true`; lands in `net.matpower_branch_metadata`. |
| `mpc.branch_kind` | Line/transformer classification override | `matpower_import.apply_branch_kind = true`; `L`, `LINE`, `ACL` force a line, `T`, `TRAFO`, `TRANSFORMER`, `2WT` a transformer. Missing or length-mismatched metadata falls back to the electrical heuristic. |
| `mpc.for001_contingencies` | FOR001-derived contingency names for validation workflows | `matpower_import.import_for001_contingencies = true`; lands in `net.for001Contingencies`. |
| `mpc.dcline` | DC-line data | `matpower_import.matpower_dcline_mode: pf_injections` by default, see above. |
| `mpc.sparlectra` | Versioned extension namespace (`mpc.sparlectra.format_version = 1`) for metadata that standard MATPOWER does not model | Recognized automatically when present; sub-blocks below. |
| `mpc.sparlectra.transformer_losses` | Transformer no-load conductance from FOR/DTF sources | One entry per transformer branch: MATPOWER branch row, original FOR/DTF identifier, terminal bus names, raw `R/X/G/B` values where available, converted per-unit values, tap/phase-shift data, source marker and allocation convention. |
| `mpc.sparlectra.links` | Busbar couplers, one `fbus tbus status` row per impedance-less coupler (`1` = closed) | Each row becomes a `BusLink`: a closed link merges its two buses into one electrical node (the same contraction CGMES switch imports use), an open link keeps them separate. Exports write the block from `net.linkVec`. Unknown bus numbers or a wrong column count are rejected. |
| `mpc.sparlectra.tap_changers` | Tap-changer nameplates, one row per transformer branch | Columns `branch tap_step tap_min_step tap_max_step tap_current_step phase_step_deg phase_min_step phase_max_step phase_current_step` plus optional `psi_deg` and `phase_du_step`; see below. |
| `mpc.sparlectra.tap_changer_model` | Roundtrip marker `'impedance_correction'` | See [Tap-impedance correction and reimport](#tap-impedance-correction-and-reimport). |
| `mpc.sparlectra.solution_written` | Marker `1` that columns 8/9 and 14-17 are a solution | Written by `writeMatpowerCasefile` with `write_solution = true`. |

**Notes**

- The file must define or return `mpc` with `mpc.version` `'2'`; decimal
  values use a decimal point.
- Bus numbers in `branch` and `gen` must exist in `mpc.bus`. Every type-3
  bus becomes the reference of its AC island; a second one in the same
  island stays a PV bus, with a warning. The units of every type-3 bus get
  reference priority 1, so they are the first to take over when their
  island loses its reference
  ([Reference priority](slack_vs_source.md#Reference-priority)). An island
  without a type-3 bus takes the bus of its best voltage-controlled unit;
  one without any generating unit cannot be solved.
- Inactive generators or branches carry status `0`; per-unit branch
  impedances must be consistent with `baseMVA` and the voltage base.
- Transformer conductance flows into `Branch.g_pu`, part of the
  transformer PI branch model and not a terminal bus shunt
  ([Branch model](branchmodel.md)); the `transformer_losses` block
  reimports without bus-shunt approximations.

### Tap-changer nameplates

Standard MATPOWER knows only the continuous `TAP`/`SHIFT` columns. With
`mpc.sparlectra.tap_changers` present, `TAP`/`SHIFT` are the neutral
position and the current steps move the live position off it (ratio grid
`tap = neutral / (1 + n * tap_step)`, the grid the tap estimation fixes to).
`psi_deg` is the PST regulating-vector direction (`0` = unspecified, default
the symmetric 90 degrees).

Two phase-changer flavors exist, and declaring both grids at once is
rejected: `phase_step_deg` steps the shift angle in degrees; `phase_du_step`
is an additional-voltage stepper (the common Delta-u PST): each step adds
`phase_du_step` per unit of additional voltage in the direction `psi_deg`,
the shift angle follows from the cascade (`atan`), and the step columns
count Delta-u steps. The tap estimation fixes a Delta-u PST linearly on that
grid.

A `tap_step` of `0` declares a pure phase shifter, which keeps it out of
the estimator's ratio mass release; declared phase changers join the mass
release in `:pst` mode along their nameplate direction. Current steps may
be fractional; short rows are zero-padded like MATPOWER optional trailing
columns. Exports write the block for every transformer off neutral or with
a non-default grid. Rows referencing non-transformer branches, unknown
branches, or positions outside the declared band are rejected at import.

## Export

```julia
using Sparlectra

# First create or import a network
net = Net(name = "export_example", baseMVA = 100.0)
# Add components to the network...

# Export the network to a Matpower case file
filepath = "path/to/output/export_example.m"
writeMatpowerCasefile(net, filepath)
```

`writeMatpowerCasefile` takes `write_solution::Union{Nothing,Bool}`; the
default reads `matpower_export.write_solution` from the active
configuration (default `true`):

```julia
writeMatpowerCasefile(net, filepath; write_solution = true)   # solved state
writeMatpowerCasefile(net, filepath; write_solution = false)  # model only
writeMatpowerCasefile(net, filepath)                          # configuration decides
```

| `write_solution` | Export |
|---|---|
| `true` (default) | `mpc.bus` `VM`/`VA` carry the solved node state and `mpc.branch` gains the standard MATPOWER result columns 14-17 (`PF`, `QF`, `PT`, `QT`, flow into the branch at each end) from the branch-flow report; the exporter does not recompute flows. `mpc.sparlectra.solution_written = 1` marks columns 8/9 and 14-17 as a solution. An unsolved network (no branch flows) warns and falls back to the 13-column model-only export. |
| `false` | A pure model file: `mpc.branch` keeps its 13 columns, `VM = 1.0`/`VA = 0.0` for all non-slack/non-PV buses (slack and PV setpoints are preserved). |

Sparlectra writes no OPF columns 18-21. Exported `.m` files with
transformer losses carry a `SPARLECTRA EXTENSION WARNING` comment, because
plain MATPOWER ignores the `mpc.sparlectra` block and computes different
transformer active losses.

A branch whose two terminal shunt arms differ (a PowSyBl transformer with
its magnetizing admittance on side 1, a line with unequal `b1`/`b2`, see
[Branch model](branchmodel.md)) has no place in MATPOWER's one `BR_B`.
The rule `asymmetric_shunts = :bus_shunt` (the only value of the keyword)
writes the symmetric part `2 * min(from, to)` on the branch row and the
excess of each terminal as a bus shunt at that bus (the from excess seen
through the tap), named in a `% branch-derived shunt of branch ...`
comment line and listed in `mpc.sparlectra.branch_shunts` (branch, bus,
Gs in MW, Bs in MVar). A MATPOWER solver reproduces the same Y-bus; the
Sparlectra reimport keeps these parts on the bus shunt as parts of their
branch, so an N-1 outage of the branch takes them away with it.

### [Generator costs](@id matpower_gencost)

Sparlectra has no optimal power flow, but a case it writes stays usable
with MATPOWER's `runopf`: the `mpc.gencost` rows are kept with their
generators and written back unchanged, in the model export and in the
`<case>_calc_<date>.m` file of a service run.

- Import: row `g` of `mpc.gencost` belongs to generator row `g`; a block
  with twice as many rows carries the reactive-power costs of the same
  units in its second half. Both cost models are kept, 1 (piecewise
  linear) and 2 (polynomial), as a [`GenCost`](@ref) in
  `ProSumer.gencost`. A block whose row count fits neither, or that holds
  an invalid row, is not imported, with a warning naming the reason; the
  network imports as before. Out-of-service generators are dropped with
  their cost rows.
- Export: the rows are written in the order of the `mpc.gen` rows, padded
  with zeros to the widest row. When only some generators carry costs (a
  unit added in Sparlectra, a slack generator the exporter adds), no block
  is written and a warning names the units without costs, because a
  partial block would assign the rows to the wrong generators.
  Reactive-power rows are written only when every unit has one.
- SCF stores the rows per machine as `extra.<machine>.gencost`
  ([Case Format](scf.md)).

## [Citation and case-file usage](@id matpower-citation)

Sparlectra references MATPOWER case names for diagnostics and comparisons
but does **not** redistribute MATPOWER case files. If you use MATPOWER
software, data formats or case files, cite MATPOWER as its guidance page
asks (<https://matpower.org/citing/>):

> R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas, "MATPOWER:
> Steady-State Operations, Planning and Analysis Tools for Power Systems
> Research and Education," IEEE Transactions on Power Systems, 26(1),
> 12-19, 2011. <https://doi.org/10.1109/TPWRS.2010.2051168>

Some case files (ACTIVSg, PEGASE, RTE) request additional case-specific
citations in their headers. Case names in this documentation (`case300.m`,
`case1354pegase.m`, `case1951rte.m`, `case_ACTIVSg10k.m`) refer to files
you obtain from MATPOWER or the original data sources under their license,
citation and redistribution terms.

## Binary case cache (`model.net_cache_enabled`)

`model.net_cache_enabled: true` is inert (every importer builds its
network directly from its own format) and logs a warning; the key stays
readable for old configuration files, and leftover `.sparlectra_net_cache`
directories can be deleted.

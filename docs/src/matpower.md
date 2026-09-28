# MATPOWER cases

Sparlectra reads MATPOWER `.m` case files (format version 2) and writes
networks back as `.m` files, including the solved state. The format itself,
the `mpc` structure with `baseMVA`, `bus`, `gen`, `branch` and the optional
`gencost`, and every column of those matrices are defined by the
MATPOWER Manual, Appendix "Data File Format":
<https://matpower.app/manual/matpower/DataFileFormat.html> (see also the
`caseformat` reference
<https://matpower.org/documentation/ref-manual/legacy/functions/caseformat.html>
and the named column constants
<https://matpower.org/docs/ref/matpower6.0/define_constants.html>). A plain
power flow needs `baseMVA`, `bus`, `gen` and `branch`; `gencost` serves
optimal power flow and is read but not used.

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

Standard MATPOWER test cases are downloaded on demand by `ensure_casefile`
(into `Sparlectra.data/mpower`, returning the absolute path); a `.jl` name
(`ensure_casefile("case14.jl")`) requests a generated `.jl` companion file.
Pass `path` for local fixtures or site-specific files and `config` for a
loaded configuration (`load_sparlectra_config("examples/configuration.yaml")`).

Several configured cases run sequentially in configured order with
`run_sparlectra_cases`; a non-empty `runtime.cases` list wins over
`runtime.case`, and batch-level performance-profile aggregation is not
supported (profile individual `run_sparlectra` calls):

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

Import without a power flow, with the transformer conventions either taken
from `SparlectraConfig.matpower` (framework runs) or passed directly (the
default `matpower_ratio = "normal"` uses branch `TAP` values as stored,
`"reciprocal"` imports their reciprocal):

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

The `examples/powerflow/matpower_import.jl` workflow also controls the Julia
thread count: `runtime.julia_threads` is resolved after the CLI and
environment overrides, and if the requested value differs from
`Threads.nthreads()` the script re-executes once with that `--threads`
setting.

```bash
julia --threads=8 --project=. examples/powerflow/matpower_import.jl
julia --project=. examples/powerflow/matpower_import.jl --julia-threads=8
```

```powershell
$env:SPARLECTRA_JULIA_THREADS = "8"
julia --project=. examples/powerflow/matpower_import.jl
```

The Web UI imports `.m` files copy-only through **Import case files** and
exposes the import-convention controls below in its form; see
[Local PowerFlow Web UI](webui.md).

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

| YAML path | Type | Default | Allowed values | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `runtime.case` | String | `case14.m` | case path/name | Compatible single-case selector and fallback when `cases` is empty. | Single case studies and benchmarks. | Invalid/missing paths. | Parse/solve scales with case size. | `runtime.cases`, runtime/profile. |
| `runtime.cases` | Vector{String} | `[case14.m, case118.m]` | non-empty case names | Ordered batch selector for `run_sparlectra_cases`; a non-empty list takes precedence over `case`. | Deterministic multi-case validation and release checks. | Empty case names or expecting `run_sparlectra` to return a vector. | Sequential parse/solve cost per case. | `runtime.case`, `run_sparlectra_cases`. |
| `model.auto_profile` | Symbol/String | `recommend` | `off`, `recommend`, `apply` | Run a MATPOWER pre-run profile. `off` disables it, `recommend` logs decisions without changing the active config, and `apply` changes only safe import-convention recommendations with clear evidence. Solver-start and Q-limit recommendations remain logged but skipped unless configured directly. | Development, large-case investigation, reproducible robust imports. | Expecting YAML files to be rewritten; applying ambiguous diagnostics. | Low; scans existing VM/VA residuals before the solve. | Output profile visibility options. |
| `model.auto_profile_log` | Bool | `true` | `true`, `false` | Print/log auto-profile reasoning and final effective options. | Debug import decisions and reproduce final settings. | Quiet high-volume runs. | Logging overhead only. | `output.console_auto_profile`, logfile settings. |
| `matpower_import.pv_voltage_source` | Symbol/String | `gen_vg` | `gen_vg`, `bus_vm`, `auto`, `strict_check` | PV voltage setpoint source policy. | Standard MATPOWER semantics. | Nonstandard conversion assumptions. | None. | `compare_voltage_reference`, PF starts. |
| `matpower_import.pv_voltage_mismatch_tol_pu` | Float64 | `1e-4` | nonnegative real | Tolerance for PV voltage mismatch checks. | Tight validation studies. | Overly strict noisy data. | Low. | `compare_voltage_reference`. |
| `matpower_import.compare_voltage_reference` | Symbol/String | `imported_setpoint` | `bus_vm`, `gen_vg`, `imported_setpoint`, `hybrid` | Voltage reference used for comparisons. Auto-profile recommends `hybrid` when BUS.VM / GEN.VG mismatches are detected. | MATPOWER comparison workflows. | When historical/SCADA ref should dominate. | Low. | `pv_voltage_source`, diagnostics. |
| `model.bus_shunt_model` | Symbol/String | `admittance` | `admittance`, `voltage_dependent_injection` | Bus shunt interpretation model. | Default import path. | Alternative modeling studies without residual evidence. | None. | Import convention diagnostics. |
| `matpower_import.shift_unit` | Symbol/String | `deg` | `deg`, `rad` | Phase-shift input unit. | Cases with radians metadata. | Wrong unit declaration. | Negligible. | `shift_sign`, branch shift diagnostics. |
| `matpower_import.shift_sign` | Float64 | `1.0` | real (typ. `1.0`, `-1.0`) | Phase-shift sign convention. | Cross-tool convention alignment. | Unnecessary flipping. | None. | `shift_unit`, branch-shift diagnostics. |
| `matpower_import.ratio` | Symbol/String | `normal` | `normal`, `reciprocal` | Branch ratio interpretation mode. | Standard MATPOWER import. | Unsupported alternate conventions. | None. | Transformer/tap interpretation. |
| `model.tap_changer_model` | Symbol/String | `ideal` | `ideal`, `impedance_correction` | Tap-changer model applied to all transformers after import; see [Transformer tap-changer model](configuration.md#transformer-tap-changer-model). | Cases where reference data corrects R/X with the tap position. | Standard ideal-tap-changer imports. | None (constant per-branch scale factor). | Shared with the native DTF importer; implemented centrally in `src/equicircuit.jl`. |
| `matpower_export.write_solution` | Bool | `true` | `true`, `false` | Whether [`writeMatpowerCasefile`](@ref) writes the solved bus `VM`/`VA` state and branch result columns 14 to 17 (`PF`/`QF`/`PT`/`QT`) into the export, marked with `mpc.sparlectra.solution_written`. When `false`, `mpc.branch` keeps 13 columns and `VM = 1.0`/`VA = 0.0` for all non-slack/non-PV buses (slack and PV setpoints are preserved). | Exporting a solved case for downstream tools or roundtrip validation. | Exporting a pure, state-independent model file. | None; sourced from the existing branch-flow report, no recomputation. | If the network is unsolved, falls back to a 13-column model-only export with a warning; interacts with `model.tap_changer_model` through the `mpc.sparlectra.tap_changer_model` roundtrip marker (see [Tap-impedance correction and reimport](#tap-impedance-correction-and-reimport)). |
| `matpower_import.enable_pq_gen_controllers` | Bool | `true` | `true`, `false` | Enable controller behavior on imported PQ generators. | Realistic controlled studies. | Raw imported behavior reproduction. | Small control bookkeeping cost. | PF Q-limit behavior. |
| `model.preallocate_network` | Symbol/String | `auto` | `off`, `on`, `auto` | Controls import-time `sizehint!` preallocation for large MATPOWER network construction. | Large imports where construction allocations dominate runtime. | Tiny cases where tuning is unnecessary. | Can reduce import allocations/time; no model changes. | `model.preallocate_min_buses`. |
| `model.preallocate_min_buses` | Int | `1000` | positive integer | Bus-count threshold used when `preallocate_network = auto`. | Auto-tuning preallocation trigger for site-specific case sizes. | If fixed always-on/off behavior is preferred. | Threshold only; no model changes. | `model.preallocate_network`. |
| `matpower_import.apply_bus_names` | Bool | `false` | `true`, `false` | Use standard `mpc.bus_name` metadata for imported bus names. | FOR001/FOR002 validation and named-bus workflows. | Preserve historical numeric names. | None. | Fails on duplicate names. |
| `matpower_import.apply_branch_names` | Bool | `false` | `true`, `false` | Attach user-defined `mpc.branch_name` metadata to `net.matpower_branch_metadata`. | Outage and contingency mapping. | Cases without branch metadata. | None. | `import_for001_contingencies`. |
| `matpower_import.apply_branch_kind` | Bool | `false` | `true`, `false` | Use user-defined `mpc.branch_kind` to override line/transformer classification. | Conversion workflows that know row kinds. | Prefer electrical heuristic. | None. | Accepts `L`/`LINE`/`ACL` and `T`/`TRAFO`/`TRANSFORMER`/`2WT`. |
| `matpower_import.import_for001_contingencies` | Bool | `true` | `true`, `false` | Preserve user-defined `mpc.for001_contingencies`. | FOR001/FOR002 validation. | Ignore validation metadata. | None. | `mpc.branch_name` enables index mapping. |
| `matpower_import.matpower_dcline_mode` | Symbol/String | `pf_injections` | `reject_active` (deprecated), `ignore_inactive`, `pf_injections`, `paired_control` | Controls active `mpc.dcline` rows. | Use `pf_injections` to emulate MATPOWER simple PF DC-line injections; `paired_control` additionally attaches one steerable `HvdcPairControl` per row (transfer `PF`, losses `LOSS0`/`LOSS1`, see [HVDC Back-to-Back](hvdc_back_to_back.md)). | OPF/dclinecost studies. | Adds two fixed prosumers per active row in `pf_injections`; `paired_control` adds one controller per row on top. | Default supports active DC-line rows as fixed terminal injections. Configuration files carrying the deprecated `reject_active` load as `pf_injections` with a warning. The strict fail-fast check stays available programmatically: `createNetFromMatPowerCase(matpower_dcline_mode = :reject_active)`. |

Example configuration for a case-conversion or validation workflow:

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
  matpower_dcline_mode: pf_injections
```

### Auto-profile pre-run

The MATPOWER runner evaluates `model.auto_profile` before the main solve: it
reads the case, computes compact diagnostics and prints a table with option
path, current value, recommended value, action, reason and evidence. The
shipped default `auto_profile: recommend` with `auto_profile_log: true` logs
recommendations without changing the configuration.

| Mode | Effect |
|---|---|
| `off` | No pre-run. |
| `recommend` | Diagnostics are logged, the run configuration is unchanged. |
| `apply` | Only clearly safe convention changes are applied: `shift_unit`, `shift_sign`, `ratio` and `bus_shunt_model` when residual scans show a large, unambiguous improvement. PV/REF voltage, robust-start and Q-limit guard settings stay recommendation-only (`keep` or `skipped` in the table). |

Diagnostics used: VM/VA power-balance residual scans over branch shift
unit/sign and transformer ratio conventions; bus-shunt residual comparison
with and without the MATPOWER bus shunt admittance; PV/REF `BUS.VM` versus
online `GEN.VG` mismatch counts; case size, PV-bus count and generator
Q-range heuristics for robust-start and Q-limit recommendations. For large
cases the pre-run keeps start projection disabled, uses DC-angle and
blended-voltage flat-start seeds, recommends a practical validation
tolerance, and disables expensive diagnostics unless requested.

Explicit YAML values stay visible; in `apply` mode a changed option is shown
as `applied` in the table and in the final effective options block.
Auto-profile never rewrites user YAML files. To reproduce a run, copy the
logged final effective options (or enable
`diagnostics.log_effective_config`) into a tracked configuration file.

`examples/powerflow/matpower_import_multi_config.jl` runs one MATPOWER case
against several YAML files (repeated `--config=...` or a comma- or
semicolon-separated `--configs=A,B,C` list) to check whether
`model.auto_profile`, `power_flow.wrong_branch_detection` or start-mode
settings change the final rectangular solver status; `--status-only` prints
the status and wrong-branch diagnostic fields per configuration, `--runner`
delegates to the standard `Sparlectra.run_matpower_case` output, and the
script never creates or rewrites YAML files:

```bash
julia --project=. examples/powerflow/matpower_import_multi_config.jl \
  data/mpower/case14.m \
  --config=path/to/config_a.yaml \
  --config=path/to/config_b.yaml \
  --status-only
```

### DC lines and disconnected AC islands

`matpower_import.matpower_dcline_mode = pf_injections` (the default) imports
each active `mpc.dcline` row with the MATPOWER `toggle_dcline`-compatible
approximation: two generator-like terminal prosumers, from-side `PG = -PF`,
to-side received power `PF - (LOSS0 + LOSS1 * PF)` when loss columns exist
(else the input `PT`), row `QF`/`QT`, voltage setpoints `VF`/`VT`, and
terminal Q limits where present. Terminal buses with voltage setpoints
become voltage-controlled where MATPOWER would make them PV; reference buses
stay reference buses and isolated buses are not activated. API and Web UI
runs write `matpower_dcline.csv` describing the mapping.

- `reject_active` (strict, deprecated as a configuration value): an active
  row (`status != 0`) aborts before solving with
  `failure_reason = unsupported_matpower_dcline`.
- `ignore_inactive`: the same active-row rejection, documenting the
  inactive-row ignore policy. Empty or inactive-only tables are tolerated in
  every mode.
- `paired_control`: one steerable `HvdcPairControl` per row, see
  [HVDC Back-to-Back](hvdc_back_to_back.md).
- This is a power-flow approximation, not an HVDC converter or DC-grid
  model: OPF constraints, converter controls, `dclinecost` and DC-line
  optimization are unsupported.

No AC branches, dummy admittances or tiny impedance bridges are created
between the terminals, so large cases such as `case_SyntheticUSA.m` consist
of several AC Ybus components coupled only through DC-line injections.
`power_flow.islands.enabled = true` is the default; with
`power_flow.islands.mode = solve_independent` each AC island is solved by the
rectangular NR solver and the states are merged into one result. The
artifact `ac_islands.csv` records per island the reference bus, status,
active DC-line terminal count, power totals and pre-slack active-power
imbalance.

### Tap-impedance correction and reimport

A case built with `model.tap_changer_model = impedance_correction` exports
`BR_R`/`BR_X` values that already carry the correction (see
[Transformer tap-changer model](configuration.md#transformer-tap-changer-model));
`writeMatpowerCasefile` marks this with
`mpc.sparlectra.tap_changer_model = 'impedance_correction'`. On reimport,
`createNetFromMatPowerFile`/`createNetFromMatPowerCase` detect the marker
and skip `calcTapCorrectedRX` whatever `model.tap_changer_model` says. Cases
without the marker, including third-party MATPOWER cases, import unchanged.

## What Sparlectra reads

A power-flow run uses `bus`, `gen` and `branch`; with Q-limit handling
enabled, `QMAX` and `QMIN` decide whether a generator can hold its voltage
setpoint, and the configured strategy decides how a limited bus is treated.
Additional MATLAB functions inside the `.m` file are not supported. The
optional fields below are Sparlectra extensions, not MATPOWER fields: when
absent, standard imports work unchanged; when present with the `apply_*`
option at `false`, naming and branch classification keep their defaults.
The direct importer takes the same switches as keywords
(`apply_bus_names = true`, `apply_branch_names = true`,
`apply_branch_kind = true`, `import_for001_contingencies = true`).

| Field | Used for | Notes |
|---|---|---|
| `mpc.baseMVA` | Common power base for per-unit quantities | Typically `100`. |
| `mpc.bus` | Network nodes: type (`1` PQ, `2` PV, `3` reference, `4` isolated), demand, shunts, `VM`/`VA` start or reference values, base voltage, limits | Lines and transformers refer to `BUS_I`. |
| `mpc.gen` | Generators: `PG`, `VG` setpoint, `QMAX`/`QMIN` for Q-limit handling, status | `PG`, `VG`, `QMAX` and `QMIN` are the values to inspect first. |
| `mpc.branch` | Lines and transformers: `BR_R`/`BR_X`, charging, ratings, `TAP`, `SHIFT`, status | A non-zero tap ratio or phase shift marks a transformer; `TAP = 0` is treated as ratio `1`. |
| `mpc.gencost` | Optimal power flow | Read, not used; never written. |
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
- Bus numbers in `branch` and `gen` must exist in `mpc.bus`, and each
  connected island needs at least one in-service reference bus.
- Inactive generators or branches carry status `0`; per-unit branch
  impedances must be consistent with `baseMVA` and the voltage base.
- Transformer conductance lives on `PowerTransformerWinding.g` and flows
  through `getTrafoRXBG`/`getTrafoRXBG_pu` into `Branch.g_pu`, part of the
  transformer PI branch model rather than a synthetic terminal bus shunt;
  Sparlectra reimports the `transformer_losses` block without adding
  bus-shunt approximations.
- Support for individual network data issues is beyond the scope of this
  project; users are encouraged to resolve such issues independently and
  share their results with the community.

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

A `tap_step` of `0` declares a pure phase shifter, which keeps it out
of the estimator's ratio mass release; declared phase changers join the mass
release in `:pst` mode along their nameplate direction. Current steps may be
fractional; short rows are zero-padded like MATPOWER optional trailing
columns. Exports write the block for every transformer off neutral or with a
non-default grid. Rows referencing non-transformer branches, unknown
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
| `true` (default) | `mpc.bus` `VM`/`VA` carry the solved node state and `mpc.branch` gains the standard MATPOWER result columns 14-17 (`PF`, `QF`, `PT`, `QT`, flow into the branch at each end) from the branch-flow report path; the exporter does not recompute flows. `mpc.sparlectra.solution_written = 1` marks columns 8/9 and 14-17 as a solution. An unsolved network (no branch flows) warns and falls back to the 13-column model-only export. |
| `false` | A pure model file: `mpc.branch` keeps its 13 columns, `VM = 1.0`/`VA = 0.0` for all non-slack/non-PV buses (slack and PV setpoints are preserved). |

Sparlectra writes no OPF columns 18-21 and no `mpc.gencost`. Exported `.m`
files with transformer losses carry a `SPARLECTRA EXTENSION WARNING` comment
because plain MATPOWER ignores the `mpc.sparlectra` block.

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

!!! details "Why it is built this way"
    Transformer no-load conductance has no branch-local field in standard
    MATPOWER, and mapping it onto a terminal bus shunt would change the
    active losses of the transformer. The `mpc.sparlectra.transformer_losses`
    block keeps the conductance on the branch across an export and reimport;
    plain MATPOWER ignores the block and computes different transformer
    active losses, which the comment in the exported file says.

    The `mpc.sparlectra.tap_changer_model` marker exists because a second
    tap-impedance correction on reimport would stack a differently derived
    factor on the already corrected values: the MATPOWER reimport uses
    `1/ratio²`, the native DTF importer the `tap_fraction`-based
    regulating-vector factor. Skipping `calcTapCorrectedRX` when the marker
    is present keeps the transformer impedances bit-identical between the
    native and the roundtrip case.

    DC lines are imported as fixed terminal injections rather than as an
    HVDC model because that is what MATPOWER's own `toggle_dcline` power
    flow does; the same approximation gives comparable results and needs no
    converter parameters the case does not carry. The price is that the AC
    components joined only by DC lines are separate islands, which the
    island solver handles instead of an artificial bridge branch.

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

`model.net_cache_enabled: true` is inert (every importer builds its network
directly from its own format) and logs a warning; the key stays readable for
old configuration files, and leftover `.sparlectra_net_cache` directories
can be deleted.

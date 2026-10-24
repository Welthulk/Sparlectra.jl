State Estimation Measurements
==============================

Measurement sets travel as measurement tables or inside a case file. This
page covers the files, the measurement generator, the case-file binding and
the power flow chained onto an estimate. The estimator itself is on
[State Estimation](state_estimation.md); the format of the table (columns,
keys, units, signs) on [Measurement table](state_estimation_measurement_table.md).

## [Measurement files and the SE chain](@id se-measurement-files)

The measurement table (first line `# sparlectra-measurement-table v1`, one
row per measured value, keyed by element name or CGMES mRID, values in MW,
MVAr, A, kV and degrees) is the measurement format of Sparlectra. The Web
UI generator and "add noise" write it, and the sets shipped with the demo
cases are tables. `readMeasurementTable!(net; file)` reads one set into
`net.measurements`; `writeMeasurementTable(file, rows)` writes rows, for
example from `generateMeasurementTableFromPF`. The measurement CSV v1 is
no longer read (see [Removed: measurement CSV v1](@ref se_measurement_csv_v1)).

**Run artifacts.** An estimation run writes the set it estimated as a
measurement table, `measurements.csv` in the run directory, and the
verdict per row in `measurements_verdict.csv` next to it. The name
`measurements.csv` is kept: the Web UI links and the history know it, and
it matches the `<case>.measurements.csv` sets. The table loads directly
with `readMeasurementTable!` (for a CGMES case with
`index = cgmesMeasurementKeys(result)`) or as the measurement file of a
new run, and gives the same measurement set: the run estimates exactly
the set this file reads back to.

- A table set keeps the rows that became a measurement, in file order,
  with their keys, value text, source and valid flag; rows that did not
  are listed in `measurement_import.txt`. A row without a sigma cell gets
  the sigma the run used and the source token `sigma_level` or
  `sigma_default`, so the file reads back the same under any
  configuration.
- Zero-injection rows the run added are bus balance rows with source
  `se_run;zero_injection`; a set carried by a case file is written with
  source `case_file`. A measurement the table cannot carry (a flow on a link
  without a key) is named as a comment line and stays in the estimate.
- The comments of the input table are carried over (case binding,
  generator provenance, tap table, truth values).

The verdict file has one line per data row of the table, joined by its
line number in `measurements.csv` (a key with its type is not unique: a
bus sum covers the rows of all its prosumers, and each of them gets the
verdict of the sum). Its CSV format follows `output.csv_format`.

| column | content |
|---|---|
| `line` | line number of the row in `measurements.csv` |
| `key`, `type` | the key cell and the type of the row |
| `measurement_id` | id of the measurement the row belongs to |
| `residual` | measured minus estimated, in the unit of the row (MW, MVAr, A, degrees, `U` in kV) |
| `normalized_residual` | residual over its standard deviation; 0 for a zero-injection constraint |
| `wii` | residual sensitivity of the row |
| `localizable` | `wii` above the localizability threshold |
| `eliminated` | removed by the bad-data elimination; residual and normalized residual are those of the solve before its removal |
| `constraint` | zero-injection row solved as a constraint (`state_estimation.zero_injection: constraint`) |
| `reason` | tokens separated by `;`: `inactive` (valid=false), `deenergized` (the row sits at a bus that is not energized, or at one of its branches, and left the run), `range_check_substituted` (a voltage reading outside the range, replaced by a pseudo measurement; residual and normalized residual are those of the substitute, see [Range check of the voltage measurements](@ref se-range-check)), `eliminated`, `suppressed`, `downweighted`, `suspicious`, `critical` |

Every other active row carries the residuals of the final diagnostics
(after every elimination). The suspicious count in `run.log`, on the
result page and in `se_suspicious` is taken from the same diagnostics;
the pass before the elimination lists its own suspects in
`se_diagnostics.md`.
A run directory written before 0.50.0 holds the CSV v1 artifact
(`bad_data` and `status` behind `id`); the Web UI shows such a run, and
the topology hypothesis test asks for a new run instead of reading it.

**Comment blocks.** A table carries them as `#` lines before its header.

| comment | written by | content |
|---|---|---|
| `# case: <name>` | generator | the case binding: the SE service refuses a set bound to a different case up front, a set without a binding still runs, with a log note |
| `sparlectra-taps v1` | generator | one row per transformer with a tap changer: the changer the deviation moves, its neutral ratio and step size, the absolute electrical, fixed and transferred step of the generation state and the deviation in steps |
| `# truth_value,...` | generator | the noise-free truth value of every row (source of `se_deltas.csv`, see [Measurement generator v2](@ref se-generator-v2)) |
| `# flow_end,...` | generator | the flow-end choice per branch |
| `# generator: v2` | generator | the set's provenance (truth source, flow-end choices, passive handling) |

**Bad data.** The verdict file names the suspicious, eliminated,
suppressed (replacement) and down-weighted (staged) rows. Until 0.50.0
`se_bad_data.csv` listed them; the verdict file replaces it. The tap steps
of a run are in `se_tap_estimates.csv`. Eliminated rows are inactive in
the reported run (its `J` is the `J` without them); `se_diagnostics.md` names the run, the case, the set and
the reported numbers, then the pass that found them ("Diagnostics pass").
An elimination counts only when the control solve without the row
converges, with released taps also without `runse!` freezing them again
where the run before kept them released; otherwise the row stays active, the attempt is named
("Reverted") and the sequence stops with `not_converged`. With released
taps the service eliminates in rounds on the end steps of the result run,
one row per round, see [Tap estimation switch and fallback](@ref se-tap-fallback);
every row of all rounds is `eliminated`, each with the residual of the
round that removed it. The
topology hypothesis test reads `measurements.csv` and leaves out the rows
the verdict file marks as eliminated.

The set tab of the Web UI (download, re-upload, inline editor, tap comment
table) is described in [Web UI](webui_reference.md#State-estimation).

## Measurement generator

The generator solves the case once (island-wise) and writes a new set
`<case>.measurements.<yyyymmdd-HHMMSS>.csv` next to the case, never over an
earlier one ("add noise" writes `<case>.noisy.<stamp>.csv`). A generated
set takes precedence over the one a case file carries. Set selection and
deletion: [Web UI](webui_reference.md#State-estimation).

| | |
|---|---|
| Call | `generateMeasurementsFromPF(...)` / `setMeasurementsFromPF!(net; ...)` with sigmas from `measurementStdDevs`, floors from `measurementSigmaFloors`, `relativeSigma = true` for percent sigmas; passive buses via `findPassiveBuses` and `addZeroInjectionMeasurements!` |
| Web UI | generator tab of the state-estimation section; "Reset saved settings" deletes the per-case settings sidecar |

### [Measurement generator controls](@id se-generator-controls)

#### [Seeded noise](@id se-generator-noise)

Adds Gaussian noise at the per-quantity sigmas to every generated value.
The random generator is seeded, so regenerating with the same inputs
reproduces the identical file.

Without noise the values are exact solutions of the power flow: the
estimation reports `J` close to 0 and the band test flags `:low`. That is
expected for a noise-free synthetic set, not an error.

#### [Bad data (k times sigma)](@id se-generator-gross-error)

Corrupts seed-randomly drawn telemetry rows by k times their sigma. `0` =
off; `10` gives a clearly detectable bad measurement (a stuck transducer).
**Bad data rows** says how many rows get the gross error; the seed decides
which, protected `ZI` rows are never drawn, and the confirmation message
names every corrupted row. The diagnostics should rank exactly these rows
first, the sequential elimination remove them, the robust option suppress
them without removal.

#### [Tap deviation](@id se-generator-tap-error)

Generates the measurements from a network state whose transformers run the
given number of mechanical tap steps off the model position (0 = off;
whole steps only); the model file stays untouched. **Transformers (max)**
bounds how many transformers deviate; the seed picks them from the
in-service, non-machine transformers with a tap changer the estimation
can release, and the message names the affected branches. The generator
does not read `state_estimation.tap_estimation.exclude`. When the set puts
a deviation on a transformer the run keeps out, run.log names it. The
deviation moves one changer per transformer, the first ratio changer,
else the first phase changer. A machine transformer is used only when the
case has nothing else. The tap estimation keeps machine transformers out
by default (`state_estimation.tap_estimation.include_machine_transformers = false`)
and calculates their step after the estimate, so a deviation there stays
in the estimate as a model error (`J` far above `dof`). run.log warns in
that case and names the key.

The deviation lies on a physical step of the changer (ratio, angle and
the impedance follow the changer data by the rule of the estimator, see
Impedance per step in [State estimation extensions](state_estimation_extensions.md)),
the step the estimator fixes on. When the set makes the step observable
(see [Observability per changer](@ref se-tap-observability)), tap
estimation recovers it exactly and the fixation run ends near `J = 0`. A
changer the set does not observe keeps its transferred step
(`status = not_estimated`), and the deviation stays as a model error.
Without tap estimation the band test reports `:high` with the suspects
clustered at that transformer. Off-grid positions are not a generator
setting; the estimator reports them through the off-grid note.

#### [Measurement sigmas (U, I, P, Q)](@id se-generator-sigmas)

Measurement accuracy per quantity in percent of the measured value (like a
transducer accuracy class), written into the sigma column of every
generated row and used for the noise, when enabled:

- **sigma U** (%): all voltage-magnitude rows. `0.5` corresponds to a class
  0.5 device.
- **currents (I)** checkbox plus **sigma I** (%): branch current-magnitude
  rows at both branch ends, generated only when the checkbox is set.
  Currents are auxiliary in the estimator, see
  [`ImagMeas`](state_estimation.md#Branch-current-magnitude-measurements-(ImagMeas)).
- **sigma P** (%): active-power injections and branch flows.
- **sigma Q** (%): reactive-power injections and branch flows.

Percent of reading keeps one setting meaningful across voltage levels.
Near-zero readings get a per-type floor instead of a near-zero sigma
(`measurementSigmaFloors`: 1e-4 pu, 0.05 MW/MVar, 0.1 A); the power floor
also keeps zero-injection buses from stiffening the flat-start solve. A
healthy noisy set lands near `J = dof`; sigmas too pessimistic against the
noise push the band test to `:low`, too optimistic ones to `:high`.

Where `state_estimation.sigma_by_voltage_level` has a value for the voltage
level of the measuring bus and the quantity, that sigma is written instead
of the percent field and the noise is drawn with it.

#### [Measurement resolution (U, P, Q, I, angle)](@id se-generator-resolution)

Rounds each generated value after noise and bad data to the step an online
system delivers. Default off (`0`, no rounding) for every quantity; steps in
U kV, P MW, Q Mvar, I A and angles in degrees (form fields
`gen_resolution_u_kv`, `gen_resolution_p_mw`, `gen_resolution_q_mvar`,
`gen_resolution_i_a`, `gen_resolution_angle_deg`; library keyword
`resolution` of `generateMeasurementsFromPF`). The written sigma stays the
configured one: the rounding is part of the measurement error (standard
deviation step / sqrt(12)), so a step coarser than about sigma times
sqrt(12) adds visible error. A rounded set names the steps in a
`# resolution:` comment line; zero-injection and passive balance rows are
not rounded.

### [Measurement generator v2](@id se-generator-v2)

Truth state, flow placement, passive-node handling and a target number of
critical measurements are configurable in the Web UI form.

**Options**

| Option | Choices | What it does |
|---|---|---|
| Truth state | `fresh solve` (default) | solves the case now (tolerance tightened to at most 1e-8, island-wise); the confirmation names iterations and tolerance |
| | `from run` | adopts the solved voltages of a successful run of the same case from the run history (SE runs through `se_state.csv`, power-flow runs through `bus_voltages_complex.csv`, any CSV export format), bit-exact against the artifact |
| Flow measurements per branch | `both ends` (default) | measures P/Q (and I) at from and to |
| | `one end (balance-aware)` | keeps one flow group per branch; the end at a bus with injection telemetry is preferred (`ZI` rows do not count), the from end wins when both or neither qualify; the choice is documented per branch in the set's `# flow_end,...` comments |
| Passive nodes | checkbox (`state_estimation.zero_injection`, default on; `constraint` shows as checked) | passive buses get zero-injection rows (prefix `ZI`, at the configured sigma, at least `state_estimation.zero_injection_sigma_mva`) and no duplicate injection rows; see [Zero injections and coupler balances](@ref se_zero_injection) |
| | off | explicit zero balance rows `Pinj/Qinj = 0` at the configured sigma (default 0.05 MW/MVar) |
| | both | the checkbox is the case's `state_estimation.zero_injection`: it is seeded from the configuration, and a generation writes it into the case configuration file, so the SE runs of the case (measurement CSV, case file, measurement table) follow it; see [State-Estimation Configuration](state_estimation_configuration.md) |
| `critical measurements` | count (0 = off) | thins the set until that many rows are critical; the set comments name the target, the rows that ended critical and the rows removed; when the target is out of reach, the message says how many were reached |

**Notes**

- `from run` rejects a missing artifact, a foreign case or an unknown run
  id up front, and warns (in the confirmation and in the file's `# truth:`
  line) when the source run is no solution of the model (state estimation
  outside the band, taps not estimated after a tap fallback, not converged): flows
  derived from such a state violate the node balances. The tap-deviation
  option is locked with `from run`.
- Critical-measurement thinning reads the criticality from `diag(Omega)`
  (see [Observability](observability.md)) and removes per step the
  strongest partner (largest normalized covariance
  `|Omega_ij| / sqrt(Omega_ii Omega_jj)`) of the telemetry row closest to
  critical; above the dense-path size caps it removes the redundant row
  with the smallest `wii`. Zero-injection and passive balance rows are
  never removed; a removal that fails the rank test is undone.
- **Delta file.** An SE run on a set with `# truth_value,...` comments
  writes `se_deltas.csv`: per measurement the measured, truth and estimated
  values with the three deltas, the sigma, the normalized residual and the
  elimination flag, plus one row per released tap changer comparing the
  documented deviation (0 for a changer the generator did not move) with
  the deviation of the fixed step.

## [Case files carry their own measurements](@id se-case-file-measurements)

A [Sparlectra Case Format](scf.md) case stores its measurement set inside
the file. An estimation on such a case needs no separate CSV; the run still
writes `measurements.csv` into its output directory for reproducibility.
Naming a set explicitly takes precedence, and a set generated for the case
replaces the one the file carries (the newer data). Every other format
requires a measurement file.

| | |
|---|---|
| Provenance | a set records whether it carries noise; a case file keeps that in `sparlectra.measurements.provenance`, together with the per-row truth values and the tap deviations the generator applied |
| Add noise | [`addMeasurementNoise!`](@ref) perturbs each value with the sigma its own row declares, without recomputing anything |
| Web UI | the page names the state of the carried set (ideal or measured) and offers "Add noise to this set", see [Web UI](webui_reference.md#State-estimation) |

A noise-free set returns $J = 0$ by construction, which says nothing about
the estimate. A run on a case file produces the same `se_deltas.csv` and
tap warning as a run from a measurement file, and a set bound to the case an SCF
file was exported from stays usable, because the file records its source.
A set in which every quantity appears twice (one generation appended onto
existing rows) inflates $J$ without any single residual looking wrong,
which reads as a topology error; the reader and the run report the
repeated rows, and regenerating the set is the fix.

## PF from the estimate

A power flow can start from the last estimate, either as start values only
or as the study base with the estimated nodal balances. The estimated state
persists in `se_state.csv`, which is how the Web UI hands an SE result to a
later power-flow or N-1 run. Both modes leave the persistent model
untouched ([Write-back policy](state_estimation.md#Write-back-policy)).

| | |
|---|---|
| Call | `runpf_from_se!(net, maxIte, tol, verbose; mode)` after `with_state_estimation_config(() -> runse!(...); update_net = true)` or a state restored via `readSEStateCSV!` |
| `mode = :se_state` | config `profile_source = state_estimation`: start from the estimated voltages; the model injections stay authoritative, the measurement/model difference goes into the (possibly distributed) slack |
| `mode = :se_snapshot` | config `profile_source = se_snapshot`: also take the nodal balances from the estimation, through temporary delta load prosumers removed after the run, the slack exempt |
| State file | `writeSEStateCSV`/`readSEStateCSV!` (`se_state.csv`); the N-1 report metadata records the SE run id |
| Service naming | request key `se_mode`; an SE run reports `run_mode = "se"`, the chained power flow `run_mode = "powerflow_se_start"` (with `se_run_id` and `se_start_mode`), an SE-started N-1 keeps `run_mode = "contingency"` plus the same two keys |
| Web UI | "Run power flow from this estimate" on the SE result page |

With consistent PV setpoints the snapshot power flow converges in 0 or 1
iterations and `slack_pickup_mw/mvar` stays below tolerance; active
outer-loop controllers (taps, machine control, Q-limits) still iterate and
legitimately move the operating point, which shows as a small residual
pickup.

## [Removed: measurement CSV v1](@id se_measurement_csv_v1)

The reader of the measurement CSV v1 (`readMeasurementsCSV!`) was removed
in 0.50.0. The SE service refuses a v1 file with
an error that names the converter, and so does "add noise" in the Web UI.
The Web UI export and "Save case as" write the case without a v1 set and
name the skipped file and the converter in their message. Convert a v1 file and its
case with `convertMeasurementsV1` or `sparlectra convert-measurements` (see
[Migrating from the v1 CSV](@ref se_table_migration)); the converter still
reads every v1 file.

`writeMeasurementsCSV(net; file)` stays and writes this layout, the input
of the converter. The first line is the version comment, the second the
header:

    # sparlectra-measurements v1
    type,bus,from_bus,to_bus,branch_nr,link_nr,direction,value,sigma,active,id

| column | meaning | rule |
|---|---|---|
| `type` | measurement type (`VmMeas`, `PinjMeas`, ...) | an unknown type name rejects the whole file |
| `bus` | bus of a bus-referenced row | location group 1; buses by name, matched whitespace-insensitively; on CGMES-sourced nets the writer can reference buses by their mRID (`busReference = :mrid`); the reader resolves names, component ids and mRIDs alike |
| `from_bus`, `to_bus` | branch endpoints | location group 2, together with `branch_nr` |
| `branch_nr` | branch index | disambiguates parallel branches; the reader verifies it against the named endpoints and rejects a mismatch. Files without a `branch_nr` column are still read; their branch rows resolve by bus pair and reject parallels (same rule as `addPflowMeasurement!`) |
| `link_nr` | link index | location group 3 |
| `direction` | branch end (`:from`/`:to`) | branch rows only |
| `value`, `sigma` | measured value and standard deviation in the type's unit (Vm in pu, powers in MW/MVAr, currents in A, angles in degrees) | numbers keep their shortest round-trip form |
| `active` | active flag | inactive rows are kept in the file |
| `id` | measurement id | prefixes `ZI` (zero injection) and `SHDERIV` (derived shunt Q) mark protected rows |
| `bad_data`, `status` | run artifact before 0.50.0 only | the converter ignores both columns |

Exactly one location group per row (`bus` |
`from_bus`+`to_bus`+`branch_nr` | `link_nr`). Delimiter and decimal
separator follow `output.csv_format` of the writing run or session
(`excel_de`: semicolon and decimal comma); the converter tells the format
from the header line and reports errors with file and line
(`file:line: reason`).

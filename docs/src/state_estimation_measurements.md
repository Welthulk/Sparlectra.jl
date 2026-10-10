State Estimation Measurements
==============================

Measurement sets travel as CSV files or inside a case file. This page
covers the file format, the measurement generator, the case-file binding
and the power flow chained onto an estimate. The estimator itself is on
[State Estimation](state_estimation.md).

## [Measurement files and the SE chain](@id se-measurement-files)

`writeMeasurementsCSV(net; file)` and
`readMeasurementsCSV!(net; file, replace = true)` persist measurement sets
as UTF-8 CSV. The first line is the version comment, the second the header:

    # sparlectra-measurements v1
    type,bus,from_bus,to_bus,branch_nr,link_nr,direction,value,sigma,active,id

### File format

| column | meaning | rule |
|---|---|---|
| `type` | measurement type (`VmMeas`, `PinjMeas`, ...) | an unknown type name rejects the whole file |
| `bus` | bus of a bus-referenced row | location group 1; buses by name, matched whitespace-insensitively. On CGMES-sourced nets the generator references buses by their preserved ENTSO-E mRID (`writeMeasurementsCSV(busReference = :mrid)`); the reader resolves names, component ids and mRIDs alike; MATPOWER and DTF nets keep bus numbers/names |
| `from_bus`, `to_bus` | branch endpoints | location group 2, together with `branch_nr` |
| `branch_nr` | branch index | disambiguates parallel branches; the reader verifies it against the named endpoints and rejects a mismatch. Files without a `branch_nr` column are still read; their branch rows resolve by bus pair and reject parallels (same rule as `addPflowMeasurement!`) |
| `link_nr` | link index | location group 3 |
| `direction` | branch end (`:from`/`:to`) | branch rows only |
| `value`, `sigma` | measured value and standard deviation in the type's unit | numbers keep their shortest round-trip form |
| `active` | active flag | inactive rows are kept in the file |
| `id` | measurement id | prefixes `ZI` (zero injection) and `SHDERIV` (derived shunt Q) mark protected rows |
| `bad_data`, `status` | run artifact only: the `measurements.csv` of an estimation run repeats the verdict per row behind `id`, `bad_data` holds `*` on every row the diagnostics flagged, `status` says `critical measurement` where the loss of the row would leave the set unobservable | the reader ignores both columns, so the artifact reads back as the same set |

Exactly one location group per row (`bus` |
`from_bus`+`to_bus`+`branch_nr` | `link_nr`). Delimiter and decimal
separator follow `output.csv_format` of the writing run or session
(`excel_de`: semicolon and decimal comma); the reader tells the format from
the header line, and a round trip reproduces every row exactly. The import
is atomic with line-precise errors (`file:line: reason`).

**Comment blocks**

| comment | written by | content |
|---|---|---|
| `# case: <name>` | generator | the case binding: the SE service refuses a set bound to a different case up front, a set without a binding still runs, with a log note |
| `sparlectra-taps v1` | generator | the transformer tap positions of the generation state (electrical, fixed and transferred step, deviation percentage) |
| `# truth_value,...` | generator | the noise-free truth value of every row (source of `se_deltas.csv`, see [Measurement generator v2](@ref se-generator-v2)) |
| `# flow_end,...` | generator | the flow-end choice per branch |
| `# generator: v2` | generator | the set's provenance (truth source, flow-end choices, passive handling) |
| `sparlectra-taps-estimated v1` | estimation run (`measurements.csv` artifact) | one line per transformer with model step, estimated step, fixed step and change, when taps were estimated |

**Bad-data artifact.** `se_bad_data.csv` lists the suspicious, eliminated,
suppressed (replacement) and down-weighted (staged) rows with their network
locations. Eliminated rows are inactive in the reported run (its `J` is the
`J` without them); `se_diagnostics.md` names the run, the case, the set and
the reported numbers, then the pass that found them ("Diagnostics pass").
An elimination counts only when the control solve without the row
converges; otherwise the row stays active, the attempt is named
("Reverted") and the sequence stops with `not_converged`. After a
tap-estimation fallback the diagnostics run again on the frozen taps.

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
transformers whose deviation the tap estimation can absorb (in service,
non-machine, with a declared changer), and the message names the affected
branches. A machine transformer is used only when the case has nothing
else, because the Web UI mass release skips it and the deviation would
leave an unexplainable model error (`J` far above `dof`); the estimator run
warns in that case.

The deviation lies on the tap-fraction grid of the estimator's fixation, so
tap estimation recovers it exactly and the fixation run ends near `J = 0`;
without tap estimation the band test reports `:high` with the suspects
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
| Passive nodes | checkbox (default on) | buses without generation, load and shunt (`findPassiveBuses`) are written as protected zero-injection constraints (`addZeroInjectionMeasurements!` rows, prefix `ZI`, excluded from elimination and robust modification, at the configured sigma but never tighter than `ZERO_INJECTION_SIGMA` = 0.001 MW) and no duplicate injection rows |
| | off | explicit zero balance rows `Pinj/Qinj = 0` at the configured sigma (default 0.05 MW/MVar) |
| `critical measurements` | count (0 = off) | thins the set until that many rows are critical; the set comments name the target, the rows that ended critical and the rows removed; when the target is out of reach, the message says how many were reached |

**Notes**

- `from run` rejects a missing artifact, a foreign case or an unknown run
  id up front, and warns (in the confirmation and in the file's `# truth:`
  line) when the source run is no solution of the model (state estimation
  outside the band, taps frozen at model positions, not converged): flows
  derived from such a state violate the node balances. The tap-deviation
  option is locked with `from run`.
- Critical-measurement thinning reads the criticality from `diag(Omega)`
  (see [Observability](observability.md)) and removes per step the
  strongest partner (largest normalized covariance
  `|Omega_ij| / sqrt(Omega_ii Omega_jj)`) of the telemetry row closest to
  critical; above the dense-path size caps it removes the redundant row
  with the smallest `wii`. Zero-injection and passive balance rows are
  never removed; a removal that fails the rank test is undone.
- The zero-injection sigma floor is 1 kW (`ZERO_INJECTION_SIGMA`); a
  tighter row, squared in `G = H' W H`, pushes the normal equations past
  double precision.
- **Delta file.** An SE run on a set with `# truth_value,...` comments
  writes `se_deltas.csv`: per measurement the measured, truth and estimated
  values with the three deltas, the sigma, the normalized residual and the
  elimination flag, plus one row per released tap comparing the documented
  deviation with the fixed step.

## [Case files carry their own measurements](@id se-case-file-measurements)

A [Sparlectra Case Format](scf.md) case stores its measurement set inside
the file. An estimation on such a case needs no separate CSV; the run still
writes `measurements.csv` into its output directory for reproducibility.
Naming a set explicitly takes precedence, and a set generated for the case
replaces the one the file carries (the newer data). Every other format
requires a CSV.

| | |
|---|---|
| Provenance | a set records whether it carries noise; a case file keeps that in `sparlectra.measurements.provenance`, together with the per-row truth values and the tap deviations the generator applied |
| Add noise | [`addMeasurementNoise!`](@ref) perturbs each value with the sigma its own row declares, without recomputing anything |
| Web UI | the page names the state of the carried set (ideal or measured) and offers "Add noise to this set", see [Web UI](webui_reference.md#State-estimation) |

A noise-free set returns $J = 0$ by construction, which says nothing about
the estimate. A run on a case file produces the same `se_deltas.csv` and
tap warning as a run from a CSV set, and a set bound to the case an SCF
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

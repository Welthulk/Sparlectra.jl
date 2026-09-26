State Estimation Measurements
==============================

Measurement sets travel as CSV files or inside a case file. Here are
the file format, the measurement generator, the case-file binding and the
power flow chained onto an estimate. The estimator itself is on
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
the header line, and writing and re-reading reproduces every row exactly,
link rows included. The import is atomic with line-precise errors
(`file:line: reason`); an unknown version comment or type name rejects the
whole file.

**Comment blocks**

| comment | written by | content |
|---|---|---|
| `# case: <name>` | generator | the case binding: the SE service refuses a set bound to a different case up front, a set without a binding still runs, with a log note |
| `sparlectra-taps v1` | generator | the transformer tap positions of the generation state (electrical, fixed and transferred step, deviation percentage) |
| `# truth_value,...` | generator | the noise-free truth value of every row (source of `se_deltas.csv`, see [Measurement generator v2](@ref se-generator-v2)) |
| `# flow_end,...` | generator | the flow-end choice per branch |
| `# generator: v2` | generator | the set's provenance (truth source, flow-end choices, passive handling) |
| `sparlectra-taps-estimated v1` | estimation run (`measurements.csv` artifact) | one line per transformer with model step, estimated step, fixed step and change, when taps were estimated |

**Bad-data artifact.** `se_bad_data.csv` lists the suspicious and eliminated
measurements with their network locations, plus per row whether the robust
solve suppressed it (replacement) or down-weighted it (staged); suppressed
rows below the elimination threshold are listed too.

**Web UI.** The set tab of the state-estimation section offers download,
re-upload and an inline editor and renders the tap comment table; see
[Web UI](webui.md#State-estimation).

## Measurement generator

The demo generator solves the case once (island-wise) and writes
`<case>.measurements.csv`; a generated set replaces the one a case file
carries.

**Use**

| | |
|---|---|
| Call | `generateMeasurementsFromPF(...)` / `setMeasurementsFromPF!(net; ...)` with sigmas from `measurementStdDevs`, floors from `measurementSigmaFloors`, `relativeSigma = true` for percent sigmas; passive buses via `findPassiveBuses` and `addZeroInjectionMeasurements!` |
| Web UI | generator tab of the state-estimation section; "Reset saved settings" under the generator deletes the per-case settings sidecar |

**Options**

| Option | What it does |
|---|---|
| Sigmas | entered per quantity in percent of the measured value (U, P, Q, and I with a checkbox for current rows), which keeps one accuracy setting meaningful across voltage levels (`generateMeasurementsFromPF(relativeSigma = true)` with the `measurementSigmaFloors` floors for near-zero readings) |
| Seeded noise | on by default; the same seed reproduces the identical set |
| Bad-data injection | k times sigma (10 is clearly detectable); "bad data rows" says how many telemetry rows get the gross error, the seed decides which (protected `ZI` rows are never corrupted), and the generation message names every corrupted row |
| Tap deviation | generates the measurements from a state whose transformers run a chosen whole number of mechanical steps off the model position (the model file stays untouched); "transformers (max)" caps how many transformers get it, the seed decides which, drawn only from transformers whose deviation the tap estimation can absorb (non-machine, with a declared changer); the message notes when fewer were eligible |

**Notes**

- Without noise `J` lands near 0 instead of near `dof` (band note `:low`,
  not an alarm).
- A machine transformer is used for the tap deviation only when the case has
  nothing else, because the mass release skips it and the deviation would
  leave an unexplainable model error (`J` far above `dof`); the estimator
  run warns in that case.
- Off-grid positions are not a generator setting; the estimator reports them
  through the off-grid note.

!!! details "Why it is built this way"
    The 0.05 MW/MVar power floor models the transducer range term and keeps
    zero-injection buses from stiffening the flat-start solve.

### [Measurement generator controls](@id se-generator-controls)

The four generator controls of the Web UI, one subsection each.

#### [Seeded noise](@id se-generator-noise)

Adds Gaussian noise at the per-quantity sigmas to every generated value.
The random generator is seeded, so regenerating with the same inputs
reproduces the identical file.

Without noise the values are exact solutions of the power flow: the
estimation then reports `J` close to 0 and the Wilson-Hilferty band test
flags `:low` (residuals implausibly small against the declared sigmas).
That is expected for a noise-free synthetic set, not an error.

#### [Bad data (k times sigma)](@id se-generator-gross-error)

Corrupts one measurement, the first active-power flow row, by k times its
sigma. `0` = off; `10` is a good value: one clearly detectable bad
measurement (for example a stuck transducer). **Bad data rows** says how
many telemetry rows get the gross error; the seed decides which.

Use it to exercise the bad-data workflow: the diagnostics should rank
exactly this row first (largest normalized residual, localizable), the
sequential elimination should remove it, and the robust option should
suppress it without removal. The corrupted row is named in the
confirmation message.

#### [Tap deviation](@id se-generator-tap-error)

Generates the measurements from a network state whose first in-service
transformer runs the given number of mechanical tap steps off the model
position (0 = off; whole steps only, a tap changer has no half positions).
**Transformers (max)** bounds how many transformers deviate. The deviation
lives on the same tap-fraction grid the estimator's fixation uses, so it is
exactly recoverable: tap estimation finds the step and the fixation run
ends near `J = 0`.

The resulting measurement set is consistent in itself but disagrees with
the model around that transformer: the band test reports `:high` with
suspicious measurements clustered at the transformer. This is the test
vector for bad-data localization and for the tap estimation (the
"estimate taps" option of the estimator run resolves exactly this
discrepancy). The affected transformer branch is named in the confirmation
message.

#### [Measurement sigmas (U, I, P, Q)](@id se-generator-sigmas)

Measurement accuracy per quantity in percent of the measured value (like a
transducer accuracy class), written into the sigma column of every
generated row and used for the noise, when enabled:

- **sigma U** (%): all voltage-magnitude rows. `0.5` corresponds to a class
  0.5 device.
- **currents (I)** checkbox plus **sigma I** (%): branch current-magnitude
  rows at both branch ends, generated only when the checkbox is set.
  Currents are auxiliary in the estimator: gated at flat start and below 3
  sigma, excluded from observability.
- **sigma P** (%): active-power injections and branch flows.
- **sigma Q** (%): reactive-power injections and branch flows.

Percent of reading keeps one setting meaningful across voltage levels: 1
percent of a 400 MW flow and of a 4 MW flow both get a class-appropriate
sigma, where an absolute 2 MW sigma would be tight at 400 kV and absurd at
30 kV. Near-zero readings get a small per-type floor instead of a
near-zero sigma (`measurementSigmaFloors`: 1e-4 pu, 0.05 MW/MVar, 0.1 A),
modeling the range term of the transducer.

For a healthy noisy set the estimation should land near `J = dof`; sigmas
that are too pessimistic against the enabled noise push the band test to
`:low`, too optimistic ones to `:high`.

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

- `from run` rejects a missing artifact, a foreign case or an unknown run id
  up front; the tap-deviation option is locked with `from run` (a run state
  is a finished snapshot) and the service enforces the lock.
- Critical-measurement thinning never removes zero-injection and passive
  balance rows; a removal that fails the rank test is undone and never
  retried, so the set stays observable.
- **Delta file.** An SE run on a set with `# truth_value,...` comments
  writes `se_deltas.csv`: per measurement the measured, truth and estimated
  values with the three deltas, the sigma, the normalized residual and the
  elimination flag, plus one row per released tap comparing the documented
  deviation with the fixed step.

!!! details "Why it is built this way"
    **The 1 kW floor is measured, not chosen.** At 1e-6 a zero-injection
    row weighs a million times a normal power measurement, and the squaring
    in `G = H' W H` pushes the normal equations past double precision (on
    the 25000-bus set: no convergence, `J/dof = 3e17`, zero-injection
    residuals of 182 GW).

    **Critical measurements by partner search.** Every thinning step reads
    the criticality from `diag(Omega)` (see
    [Observability](observability.md)) and removes one row by a partner
    search on the residual covariance. The telemetry row that needs the
    fewest partner removals to become critical is the anchor; its strongest
    partner, the row with the largest normalized covariance
    `|Omega_ij| / sqrt(Omega_ii Omega_jj)` (1 for a mutually redundant
    pair), is removed. Above the dense-path size caps the step removes the
    redundant row with the smallest `wii` instead.

## [Case files carry their own measurements](@id se-case-file-measurements)

A [Sparlectra Case Format](scf.md) case stores its measurement set inside
the file. An estimation on such a case needs no separate CSV; the run still
writes `measurements.csv` into its output directory for reproducibility.
Naming a set explicitly takes precedence, and a set generated for the case
replaces the one the file carries (the newer data). Every other format
requires a CSV.

**Use**

| | |
|---|---|
| Provenance | a set records whether it carries noise; a case file keeps that in `sparlectra.measurements.provenance`, together with the per-row truth values and the tap deviations the generator applied |
| Add noise | [`addMeasurementNoise!`](@ref) perturbs each value with the sigma its own row declares, without recomputing anything |
| Web UI | the page names the state of the carried set (ideal or measured) and offers "Add noise to this set"; a case file with its own measurements is offered first, see [Web UI](webui.md#State-estimation) |

**Notes**

- A noise-free set returns $J = 0$ by construction, which says nothing about
  the estimate (on the warm-up case $J$ moves from 0 to about 43 at 42
  degrees of freedom after adding noise).
- A run on a case file produces the same `se_deltas.csv` and tap warning as
  a run from a CSV set.
- A set bound to the case an SCF file was exported from stays usable,
  because the file records its source.
- A quantity measured twice: a set in which every quantity appears twice is
  a corrupted file, usually one generation appended onto existing rows. It
  inflates $J$ without any single residual looking wrong, which reads as a
  topology error. Reading a set reports how many active rows repeat an
  already measured quantity, and the run states it in the log and its
  message. Regenerating the set is the fix.

## PF from the estimate

A power flow can start from the last estimate, either as start values only
or as the study base with the estimated nodal balances. The estimated state
persists in `se_state.csv`, which is how the Web UI hands an SE result to a
later power-flow or N-1 run.

**Use**

| | |
|---|---|
| Call | `runpf_from_se!(net, maxIte, tol, verbose; mode)` after `with_state_estimation_config(() -> runse!(...); update_net = true)` or a state restored via `readSEStateCSV!` |
| `mode = :se_state` | config `profile_source = state_estimation`: start from the estimated voltages; the model injections stay authoritative, the measurement/model difference goes into the (possibly distributed) slack |
| `mode = :se_snapshot` | config `profile_source = se_snapshot`: also take the nodal balances from the estimation, through temporary delta load prosumers removed after the run, the slack exempt; the persistent model is never mutated |
| State file | `writeSEStateCSV`/`readSEStateCSV!` (`se_state.csv`); the N-1 report metadata records the SE run id |
| Service naming | request key `se_mode`; an SE run reports `run_mode = "se"`, the chained power flow `run_mode = "powerflow_se_start"` (with `se_run_id` and `se_start_mode`), an SE-started N-1 keeps `run_mode = "contingency"` plus the same two keys |
| Web UI | "Run power flow from this estimate" on the SE result page |

**Notes**

- With consistent PV setpoints the snapshot power flow converges in 0 or 1
  iterations and `slack_pickup_mw/mvar` stays below tolerance.
- Active outer-loop controllers (taps, machine control, Q-limits) still
  iterate and legitimately move the operating point, which shows as a small
  residual pickup.

!!! details "Why it is built this way"
    `:se_state` keeps the model authoritative, `:se_snapshot` makes the
    estimated operating point the study base. Both leave the persistent
    model untouched, in line with the
    [Write-back policy](state_estimation.md#Write-back-policy).

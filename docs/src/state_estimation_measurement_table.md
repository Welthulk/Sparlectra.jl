State Estimation Measurement Table
==================================

A measurement table is a CSV file with one row per measured value and one
measurement set per timestamp. It is addressed by CGMES mRID for a CGMES
case and by element name for every other case format. It is the
measurement format of Sparlectra: the Web UI generator writes it and the
shipped sets use it. The [measurement CSV v1](@ref se_measurement_csv_v1)
is no longer read since 0.50.0; the two formats are told apart by their
first line, and a v1 file is converted once (see
[Migrating from the v1 CSV](@ref se_table_migration)).

The format is new and may change before 1.0.

## File

- CSV, UTF-8 (a BOM at the start is ignored), comma as separator, decimal
  point, LF line ends (a CR before the LF is tolerated).
- The version line comes first, then optional comment lines, then the
  header:

      # sparlectra-measurement-table v1
      # free comments: origin, generator, seed
      timestamp,measurement_mrid,terminal_mrid,node_mrid,type,value,sigma,source,valid

- Lines starting with `#` before the header are comments; blank lines are
  skipped. After the header there are no comments.
- A field with a comma or a double quote is enclosed in double quotes,
  `""` inside is one quote (RFC 4180). In practice this concerns `source`
  only. A field cannot hold a line break: each row is one line, and
  `writeMeasurementTable` refuses a cell or comment with a line break.
- Blanks around a cell are not part of it: the reader strips them and the
  writer does not write them.
- `readMeasurementTable` refuses a file without the header, without the
  version line or with another version, with an `ArgumentError` that says
  what was expected. Blanks after the version are allowed; other text
  after it is refused with its own message. `isMeasurementTable` decides by content (version line
  plus header), not by the file name. The measurement CSV v1 starts with
  `# sparlectra-measurements v1` and is no measurement table.
- A file that starts with the table header but lacks the version line is
  still treated as a table by the SE service, so the run names the missing
  line; the Web UI lists it on upload as
  `measurement_table_without_version_line`.

## Columns

| column | required | content |
|---|---|---|
| `timestamp` | no | ISO 8601 in UTC, `YYYY-MM-DDThh:mm[:ss[.f]]Z`. Empty for a single snapshot. Each instant forms its own measurement set, in the order of first appearance; `12:00Z`, `12:00:00Z` and `12:00:00.000Z` are one instant, read as `12:00Z`. An impossible date or time (`2025-02-30`, `25:61`) is `bad_timestamp`. |
| `measurement_mrid` | one of three | mRID of a `cim:Analog` (CGMES only). Type and location come from the Analog, `type` stays empty. |
| `terminal_mrid` | one of three | key of a terminal: branch end, injection or shunt. `type` is `P`, `Q`, `I` or `Ia`. |
| `node_mrid` | one of three | key of a node. `type` is `U`, `P`, `Q` or `Va`. |
| `type` | see above | `P`, `Q`, `I`, `Ia` (current angle), `U`, `Va` (voltage angle) |
| `value` | yes | the reading in the unit of the type (see Units) |
| `sigma` | no | standard deviation in the same unit, positive. Empty: the sigma of the voltage level of the measuring bus when `state_estimation.sigma_by_voltage_level` has one (see [Measurement sigma by voltage level](@ref se-sigma-by-level)), else the default of the measurement type (`measurementStdDevs`). Below the floor, the smaller of 1e-3 (the default of `state_estimation.zero_injection_sigma_mva`, in the unit of the estimator, pu for `U`) and the value of that key, it is read at the floor and counted as `sigma_lifted`, as in the measurement CSV v1. |
| `source` | no | free text. The generators write `synthetic` or `synthetic;bad_data`, the pseudo generator `pseudo;ssh`, the converter `converted;v1`; the `measurements.csv` artifact of an SE run adds `se_run`, `case_file`, `sigma_level` and `sigma_default` (see [run artifacts](@ref se-measurement-files)). One token is evaluated: `zero_injection` (tokens separated by `;`) on a bus balance row marks a zero-injection constraint (see Signs and mapping). |
| `valid` | yes | `true` or `false`. A row with `false` becomes an inactive measurement (in a bus sum see Group sums). |

Exactly one of the three key columns is set per row.

## Addressing

The names of the key columns are historical. In both cases they hold keys;
which kind of key depends on the case.

**CGMES case** (`cgmesMeasurementKeys`):

- `terminal_mrid`: mRID of a Terminal.
- `node_mrid`: mRID of a TopologicalNode, ConnectivityNode or
  BusbarSection.
- `measurement_mrid`: mRID of an Analog. Its `measurementType` gives the
  type (`ThreePhaseActivePower` P, `ThreePhaseReactivePower` Q,
  `LineCurrent` I, `LineToLineVoltage` U), its terminal the location
  (without a terminal its PowerSystemResource, when that is a node or a
  busbar). Other measurement types are reported as `unknown_type`. A `U`
  Analog on a terminal gives the voltage of the terminal's node (also on
  a switch terminal); at a node that is not energized it is
  `not_modelled`.
- A terminal of a branch out of service or open at one terminal is
  `not_modelled` at both ends: the estimator does not take flow
  measurements on a partially open branch yet. The same holds for name
  keys.
- mRIDs are compared in canonical form: a leading `_` (rdf:ID in CGMES
  2.4.15) and a `urn:uuid:` prefix are removed.
- A bus or branch end the delivery has no mRID for (the star point of a
  three-winding transformer, the legs that end there, auxiliary buses of
  the import) is keyed by its name as for any other case (below). The
  mRID comes first: a name key exists only where there is no mRID, and
  it never replaces an mRID key.

**Any other case** (MATPOWER, SCF, DTF, PowSyBl, a net built in code;
`nameMeasurementKeys`):

- `node_mrid`: the bus name.
- `terminal_mrid`: `<branch name>:from` or `<branch name>:to` for a branch
  end; the prosumer name for a load or generator; the shunt name for a
  shunt.
- A branch name that several branches share (parallel lines) is qualified
  with the branch number, `<branch name>#<k>:from`, where `k` is the
  position in `net.branchVec` (the `branch_nr` of the CSV v1). A prosumer
  or shunt name that several of them share is qualified with its bus,
  `<name>@<bus name>`. When several at one bus share that as well (two
  generators of one voltage level carry the same generated name), each
  gets its order at the bus, `<name>@<bus name>#<k>` (prosumers first,
  then shunts, each in model order). A key that is still not unique
  resolves to `not_modelled` with the reason.
- Both `#<k>` forms are a way out, not a stable address: `k` follows the
  element order of the import and can differ after the next import or
  export. Element names should be unique; a table meant to be reused keys
  on unique names.
- Every prosumer and every in-service shunt at an energized bus is a member
  of the sum at its bus, whatever its key (see Group sums). An element
  without a name has no key; its bus sum is therefore always
  `partial_injection`.
- `measurement_mrid` does not apply; a row with it is reported as
  `unknown_key`.
- Names are compared as written, only surrounding blanks are ignored. The
  canonical mRID form does not apply: buses `_B9` and `B9` are two keys.
  Two bus names that differ only in surrounding blanks give one key; a
  `U` row on it is `not_modelled` (ambiguous), never one of the buses.

The names are the ones Sparlectra shows in its results (`comp.cName`, the
bus names of `net.busDict`). For imported MATPOWER cases the branch and
prosumer names are generated by the importer, for example
`B_ACL_110_1_2` or `Ld_110.0`; an SCF case uses the names of its
`extra` block (`Sandau_Borntal`).

## Units

| type | value and sigma in the file | in the estimator |
|---|---|---|
| `P` | MW | MW (`PflowMeas`, or the bus sum `PinjMeas`) |
| `Q` | MVAr | MVAr (`QflowMeas`, `QinjMeas`, `ShuntQMeas`) |
| `I` | A, magnitude at the branch end or the shunt bay | A (`ImagMeas`), value and sigma unchanged |
| `Ia` | degrees, current angle at the branch end | degrees (`IaMeas`) |
| `U` | kV, line-to-line | `VmMeas` in pu of the node nominal voltage |
| `Va` | degrees, voltage angle at the node | degrees (`VaMeas`) |

`P` and `Q` at a node are the bus balance in MW and MVAr (`PinjMeas`,
`QinjMeas`). Angles refer to the time base of the measuring device; the
estimator treats the offset to the slack as in the measurement CSV v1
(see [State Estimation](state_estimation.md)).

For `U` the import divides value and sigma by the nominal voltage
`Vn_kV` of the bus: `vm_pu = U_kV / Vn_kV`, `sigma_pu = sigma_kV / Vn_kV`.
No other quantity is converted: no power or current base is applied. A
`U` cell with more than 17 significant digits is divided from its decimal
text at 256 bits and rounded once. The converter
(`convertMeasurementsV1`) writes such a cell where no Float64 kV value
divides back to the pu value of the CSV v1, so the pu values come back bit
for bit.
An empty sigma takes the level sigma `u_kv / Vn_kV` when the level has
one, else the `VmMeas` default directly in pu. `I` stays in
ampere: `ImagMeas` is an ampere-valued measurement (see
[State Estimation](state_estimation.md)), the import applies no current
base.

## Signs and mapping

- `P` and `Q` at a terminal are positive into the equipment at that
  terminal (the CGMES `SvPowerFlow` convention).
- Branch end: `P`, `Q`, `I` become `PflowMeas`, `QflowMeas`, `ImagMeas` at
  that end. The Sparlectra flow direction from the bus into the branch is
  the same direction.
- Injections (loads, machines, external grids, equivalent injections,
  SVCs, mapped energy sources; for name keys every prosumer) at one bus:
  `P` and `Q` of all modelled members are summed into one `PinjMeas` and
  one `QinjMeas` with the sign turned (bus injection = minus the power
  into the equipment), `sigma = sqrt(sum sigma_i^2)`.
- Shunts at one bus: `Q` summed the same way into one `ShuntQMeas`, sign
  kept (power into the shunt, the `q_shunt` convention). The estimator
  predicts `ShuntQMeas` as the sum over all in-service shunts at the bus.
- `I` at a shunt: the bay current, `ImagMeas` at the bus, when the bus
  has one in-service shunt. With several shunts the rows are
  `unsupported`: the bay current is the magnitude of the summed currents,
  not the sum of the magnitudes.
- `Ia` at a branch end: `IaMeas` at that end.
- `P` at a shunt, `I` at an injection and `Ia` at an injection or shunt
  are `unsupported`.
- Node with `U`: `VmMeas` at the bus.
- Node with `P` or `Q`: the bus balance, generation minus load at the bus
  (injection positive), `PinjMeas` or `QinjMeas` with the value as written.
  It does not replace the prosumer rows of the same bus: where both
  exist, the balance and the prosumer sum enter the set as two redundant
  measurements.
- A bus balance row whose `source` carries the token `zero_injection` is
  a zero-injection row: it gets the id `ZI_PINJ_bus_<k>` or
  `ZI_QINJ_bus_<k>` and is protected from bad-data elimination like the
  rows `state_estimation.zero_injection` adds. Its sigma is the one of the
  row.
- Node with `Va`: `VaMeas` at the bus.
- A bus link (a closed switch or coupler joining two buses outside the
  Y-bus): `P` and `Q` on its terminals are link flow measurements, positive
  from the link's from bus to its to bus (a row at the to end changes
  sign). The estimator takes them as WLS rows (on a coupler ring they
  decide the circulation, two parallel couplers for example);
  `calcLinkFlowsSE!` gives the link flows after the estimate. The SE
  service runs it after every estimate on a case with a closed link and
  writes them to `se_link_flows.csv`. An open
  link is `not_modelled`, `I` on a link `unsupported`. Name keys:
  `<link name>:from`, `<link name>:to`; in a CGMES case the terminals of
  the switch.

## Group sums

A bus sum needs a value for every modelled member of the group at that
bus.

- All members have a row with `valid=true`: the sum is active. A
  `valid=false` row of a member that also has a valid row is not used and
  is reported as `invalid_unused`.
- Every member has a row, but at least one member only rows with
  `valid=false`: the sum is created **inactive** (the input says "do not
  use", as for a single branch row). The other members still use only
  their valid rows; an invalid row next to a valid one is
  `invalid_unused` here too.
- A member has no row at all: no sum is formed; its rows are counted as
  `partial_injection`.
- Several rows of one member (redundant meters on one terminal) are merged
  to their inverse-variance mean before the sum and counted as `merged`.

## Redundancy and duplicates

Several keys at one location are redundant measurements and all enter the
set: two Analogs on one terminal, or `U` at a busbar section and at its
TopologicalNode give two `VmMeas` at one bus. `duplicate` means the same
key with the same type twice in one set; the second row is reported.

## Row outcomes

No row is dropped silently. Each row ends as an active measurement, an
inactive measurement, or with a reason in the import report
(`measurementImportReport`, one line per problem row with its line number):

| reason | meaning |
|---|---|
| `field_count` | wrong number of fields (no set can take the row, listed as rejected) |
| `bad_timestamp` | timestamp not ISO 8601 UTC or not a valid date and time (rejected) |
| `key_count` | not exactly one of the three key columns set |
| `unknown_type` | type does not fit the key column, or an Analog of a type the table does not carry |
| `bad_value`, `bad_sigma`, `bad_valid` | value, sigma or valid not readable; also a negative `U` or `I` (magnitudes), from a row or an Analog |
| `unknown_key` | key not in the index of the case |
| `not_modelled` | key known, but not part of the solved network (switch, unmapped or disconnected equipment, de-energized node, open branch), or a name key that is still shared; the report names the reason |
| `unsupported` | quantity not modelled at this location |
| `partial_injection` | bus sum without a row for every member |
| `invalid_unused` | `valid=false` row of a group member that also has a valid row; the sum (active or inactive) uses the valid rows of that member |
| `duplicate` | same key and type twice in one set |

## Generators

Both generators return rows with the header fields; `writeMeasurementTable`
writes them with the version line. The same arguments, seed included, give
the same rows and therefore the same file.

`generateMeasurementTableFromSV(result)` works on a CGMES import result and
writes mRID keys: `U` at every energized TopologicalNode with an
`SvVoltage`, `P` and `Q` at both ends of every branch closed at both
terminals, computed
with the network model at the SV voltage state.

`generateMeasurementTableFromPF(net)` works on any net after `runpf!` and
writes name keys: `U` at the buses and `P` and `Q` at both branch ends,
computed with the network model at the solved voltages.

| keyword | default | SV | PF | meaning |
|---|---|---|---|---|
| `timestamp` | `""` | yes | yes | timestamp column of every row |
| `include_injections` | `false` | yes | yes | also `P` and `Q` of the injections and `Q` of the shunts (SV: the delivered `SvPowerFlow`; PF: the solved values per prosumer of `distributeBusResults!`, computed on a copy) |
| `coverage` | `1.0` | no | yes | share of the buses and of the branch ends that get rows, drawn by the seed |
| `noise` | `true` | yes | yes | Gaussian noise with the row sigma; `false` writes the model values |
| `sigma_pq_rel`, `sigma_pq_floor` | `0.01`, `0.1` | yes | yes | sigma of `P` and `Q`: `max(sigma_pq_rel * abs(value), sigma_pq_floor)` in MW/MVAr |
| `sigma_u_rel` | `0.002` | yes | yes | sigma of `U` relative to the reading (kV) |
| `sigma_by_voltage_level` | the configured table | yes | yes | where the level of the measuring bus has a row with a column for the type, its sigma replaces the three keywords above, for the written sigma and the noise (see [Measurement sigma by voltage level](@ref se-sigma-by-level)) |
| `bad_data_fraction`, `bad_data_sigmas` | `0.0`, `20.0` | yes | yes | share of rows with a gross error of that many sigma (random sign), marked `source = "synthetic;bad_data"` |
| `seed` | `1` | yes | yes | seed of the random stream (coverage draw, noise, gross errors) |

The sigma keywords are checked before anything is generated:
`sigma_pq_floor` and `sigma_u_rel` must be positive, `sigma_pq_rel` and
`bad_data_sigmas` finite and at least 0, otherwise an `ArgumentError`
names the keyword. Every written
sigma is then positive and passes the reader.

With `noise = false` an estimation on the table reproduces the state it was
generated from: on `sp_case14` the estimated voltages agree with the power
flow to better than 1e-15 pu (all buses and branch ends, with or without the
injection rows; the test asserts 1e-6 pu).

## Pseudo measurements

A real delivery rarely measures every load and generator.
`generatePseudoMeasurementTable(net, rows; sigma)` (name keys) and
`generatePseudoMeasurementTable(result, rows; sigma)` (CGMES import, mRID
keys) write rows for the injection members that `rows` does not observe:
`P` and `Q` of each load, generator, external grid, equivalent injection or
SVC, taken from the scheduled values of the case (for CGMES the SSH
setpoints), as power into the equipment. A member counts as observed for a
type when a well-formed row carries its terminal key, or an Analog on its
terminal, with that type; a `valid=false` row counts as well. Members are
filled one by one, so a partial bus sum becomes complete and keeps its real
rows. A unit without a Q setpoint gets no `Q` row; shunts get none.

| keyword | default | meaning |
|---|---|---|
| `sigma` | required | sigma of the `P` rows in MW |
| `sigma_q` | `sigma` | sigma of the `Q` rows in MVAr |
| `timestamp` | `nothing` | `nothing`: every timestamp of `rows`; a string: only that one |
| `sigma_by_voltage_level` | the configured table | a level row with `p_mw` (`q_mvar`) replaces `sigma` (`sigma_q`) at that level; `Sparlectra.SESigmaLevel[]` keeps the pseudo sigmas |

Both sigmas must be finite and positive (`ArgumentError` otherwise). Choose
them large against the real meters, so the real rows win. The rows carry
`source = "pseudo;ssh"` and `valid = true`; appended to the table they are
read and summed like any other row:

```julia
rows = readMeasurementTable("delivery.table.csv")
pseudo = generatePseudoMeasurementTable(result, rows; sigma = 5.0)
writeMeasurementTable("filled.table.csv", vcat(rows, pseudo))
```

On the SmallGrid fixture the second set observes every injection and gets
no pseudo rows; the first set has no injection rows and gets `P` and `Q`
for all 125 members.

## Incomplete measurement sets

A real delivery is never complete. Two kinds of rows close the gaps:

- Zero injection at passive buses: the run adds it where the set has no
  injection row, see
  [Zero injections and coupler balances](@ref se_zero_injection).
- Pseudo measurements from the schedule (the SSH setpoints, see above)
  carry a large sigma. They restore observability where the gaps are and
  lose against every real meter.

Missing flows are not derived from measured neighbours. Such a value is a
linear combination of measurements the estimator already has through the
network model, so it adds no information. It counts those measurements a
second time, and its error is correlated with theirs, while WLS treats
every row as independent. Observability does not improve; only the
residuals look better.

## SmallGrid fixture

`test/fixtures/measurements/smallgrid_synthetic.table.csv` is a table for
the ENTSO-E SmallGrid (NodeBreaker BaseCase in
`test/fixtures/cgmes_smallgrid`), built by
`tools/gen_smallgrid_measurements.jl` with `generateMeasurementTableFromSV`.
It holds two sets. The numbers below come from importing the fixture and
running `runse_diagnostics` with `state_estimation.max_eliminations = 10`.

| set | rows | measurements | before elimination | eliminated | after |
|---|---|---|---|---|---|
| `2026-01-01T12:00Z` | 865 | 862 (861 active, 1 inactive) | J = 1971.01, dof 626, J/dof 3.149 | 4 rows marked `synthetic;bad_data`: two `U` (normalized residual 20.03 and -19.87), two branch `Q` (-18.93 and 14.06) | J = 618.55, dof 622, J/dof 0.994, largest abs normalized residual 3.20, stop `:consistent` |
| `2026-01-01T13:00Z` | 1126 | 1092 (all active, injection rows in bus sums) | J = 915.01, dof 857, J/dof 1.068 | none | largest abs normalized residual 3.36, stop `:consistent` |

The first set also reports three handmade rows, one per reason
(`unknown_key`, `not_modelled` for a breaker terminal, `unknown_type`), and
one healthy branch row with `valid=false` that becomes an inactive
measurement. These numbers depend on the network model; a change of the
CGMES mapping that moves them is a change of the model, not of the table.

## Use in the SE service and the Web UI

The SE service reads a table while `state_estimation.measurement_table` is
`true` (the default); with the switch off a table is refused and the key
is named. The format of the case picks the index: mRIDs for a CGMES case,
names for every other format. The run estimates the first set of the file
and writes the import report of all sets to `measurement_import.txt`.
The run adds the zero-injection rows of
[Zero injections and coupler balances](@ref se_zero_injection) to the set;
`measurement_import.txt` and the run log give the count
(`zero-injection rows N`). `applyMeasurementSet!` does the same for a
table in a script and stores the count in `set.counts[:zero_injection]`;
its keyword `zero_injection` overrides the configured switch. The switch is
case scope: a case with its own `<case>.config.yaml` takes it from that file
(written by the generator checkbox), not from the general configuration file.
Before 0.50.0 a table run added no such rows, so a table that leaves passive
buses unmeasured estimates differently now; `off` restores the former
result. Redundant readings of one quantity (several keys at one location) are
named in the run log as a note; unlike a doubled measurement CSV v1 set
they are no duplicate finding of the run. See
[State-Estimation Configuration](state_estimation_configuration.md).

The Web UI generator and "add noise" write tables, keyed the way the
service reads them. The state-estimation section offers a table next to the
case (upload, selection, row counts per type, text and value editor,
download, delete); the run goes through the same service. See
[Web UI Reference](webui_reference.md).

### [Binding a table to its case](@id se_table_binding)

Two comment lines tie a table to a case:

    # case: sp_case14.scf.json
    # cgmes_model: 0b6f1c2e-5a7d-4e8f-9c3b-2d1e4f5a6b7c

`# case:` names the case file. The generator writes it, and the run
refuses a table whose line names another case. `# cgmes_model:` names the
model id (`rdf:about` of the md:FullModel header) of the EQ file the
table's mRIDs come from; a table may carry several. The Web UI offers such
a table for every CGMES case whose EQ header carries one of these ids, so
the binding does not depend on the file name, and several scenarios built
on one EQ all offer the table. The id is compared without `urn:uuid:`, a
leading `_` and case (`cgmesModelKey`); `cgmesEquipmentModelIds(path)`
reads the ids of a delivery from the file headers without importing it.
A `# case:` line outranks the model line. A table without either line is
offered for a case when its file name starts with the case name and a dot.
`# cgmes_description:` is free text for the reader.

The run does not refuse a table bound by the model id alone: it notes in
the run log whether the id is the EQ model of the case, and rows the case
cannot resolve are counted in `measurement_import.txt`.

In a script, `readMeasurementTable!(net; file)` does the read, the import
with `nameMeasurementKeys(net)` (or the `index` given, for example
`cgmesMeasurementKeys(result)`) and `applyMeasurementSet!` in one call and
returns the set.

The example `examples/state_estimation/state_estimation_measurement_table.jl`
runs the generator and the estimation on `sp_case14` and the bad-data
elimination on the SmallGrid fixture.

## [Migrating from the v1 CSV](@id se_table_migration)

`convertMeasurementsV1(v1_file, net; out)` writes the table of a
measurement CSV v1 file. The v1 file addresses buses and branches by name
and number, the table by keys, so the converter takes the net of the case;
for a CGMES case pass `index = cgmesMeasurementKeys(result)` so the table
carries mRIDs as the SE service expects.

```julia
net = importSCF("sp_case14.scf.json")
res = convertMeasurementsV1("sp_case14.measurements.v1.csv", net; out = "sp_case14.measurements.csv")
print(res.report)
```

From the command line: `sparlectra convert-measurements <case> <v1.csv>`
in the built application, or
`julia --project=app tools/convert_measurements.jl <case> <v1.csv> [<out.csv>]`
in a checkout; both key the table like the SE service.

- Every v1 row becomes one row with the same measurement. `Vm` becomes `U`
  in kV, written with up to 18 digits where needed so the import gives back
  the pu value bit for bit; P and Q stay in MW and MVAr, I in A, angles in
  degrees. Bus injections become node balances (`P`, `Q` at the node),
  zero-injection rows (`ZI` ids) keep their sigma and get the source token
  `zero_injection`, `active` becomes `valid`.
- The comments of the v1 file move in front of the header; the ids in the
  `truth_value` and critical-row lines are rewritten to the ids of the
  import, so tap estimation and `se_deltas.csv` work as before.
- Rows the table cannot carry (a flow on a link without a key, a shunt bay at a bus with
  several shunts, a negative `U` or `I`, a second row with the same key and
  type) are not written; the result and the comment lines of the table name
  each with its line.
- A row at an isolated bus is written with that bus key; the import
  reports it as `not_modelled`, as the SE service left such a v1 row out.

The estimate on the converted table equals the one on the v1 file: the
same J, state and eliminations on every set shipped with Sparlectra (shown
in 0.50.0, when both readers existed). The v1 reader `readMeasurementsCSV!`
was removed in 0.50.0; the SE service and the Web UI refuse a v1 file with
an error that names the converter.

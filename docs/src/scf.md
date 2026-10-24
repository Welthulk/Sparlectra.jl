# Sparlectra Case Format (SCF)

One self-describing case file (`.scf.json`): a power-grid-model (PGM)
input dataset that a PGM installation reads directly, plus a `sparlectra`
block with what PGM has no model for (slack strategy, tap changers,
component names and source ids, state-estimation measurements, study
definitions). The framework, the service and the Web UI run it like a
MATPOWER or CGMES case. The case's own settings travel in a case
configuration file next to it, `<stem>.config.yaml`.

**Use**

| | |
|---|---|
| Write | `exportSCF(net; file, ...)`; `write_case_config(file, Dict(...))` for the case configuration file |
| Read | `importSCF(file)`; without building the network `load_case_config(file)`, `scf_case_config(file)`, `scf_case_studies(file)` |
| Typed form | `SCFCase` with `read_scf_json`, `write_scf_json`, `build_net`, `net_to_scfcase`; dict entry points `scf_to_net`, `net_to_scf` |
| CLI | `sparlectra export`, with `--pgm` for the plain PGM dataset |
| Web UI | **Export as SCF case file**, **Export as plain PGM**, **Save case as**, **Download selected case** on the Case page ([Web UI](webui_reference.md#Export-and-download)) |
| Fixtures | `data/scf/sp_casePST.scf.json`, the permanent version 1 file every future reader must keep loading unchanged; `data/scf/pgm_interop.json`, a file as power-grid-model writes it |

## [Writing a case file](@id scf-writing)

`exportSCF` writes a network as one JSON file that reads back machine-exact, with its source, the intended calculations and, on request, the solved start state; the Web UI offers the same export from the Case page.

```julia
using Sparlectra

net = createNetFromMatPowerFile(filename = "case118.m")
exportSCF(net; file = "case118.scf.json",
          source_format = "matpower", source_reference = "case118.m",
          intended_calculations = ["power_flow", "state_estimation"])
write_case_config("case118.scf.json", Dict("power_flow.mode" => "auto"))
```

The Web UI export picks up a `<case>.measurements.csv` next to the case
and stores the form's overrides in the case configuration file;
**Save case as** writes the same set under a new name
([Save case as](webui_reference.md#Save-case-as)).

## Reading a case file

```julia
net = importSCF("case118.scf.json")
runpf!(net, 30, 1e-8, 0)
```

The reader builds the network through the public constructors (`addBus!`,
`addPIModelACLine!`, `addPIModelTrafo!`, `addProsumer!`, `addShunt!`,
`addLink!`) and ends with `validate!`. It refuses unknown top-level keys,
unknown keys inside `sparlectra` and unsupported component types by name,
and requires unique ids and resolvable `from_node`, `node`,
`measured_object`, `regulated_object`, `branch` and `star_node` references
of the right kind. `load_case_config(file)`, `scf_case_config(file)` and
`scf_case_studies(file)` return the case configuration file, the
deprecated in-file block and the study definitions without building the
network.

## Configuration precedence

A case that ships a `<stem>.config.yaml` computes the same numbers on every
installation: its case-scope keys resolve from explicit `config_overrides`
(API, CLI, Web UI form), then the case configuration file, then the
packaged defaults; the machine's YAML supplies machine-scope keys only
([Merge precedence](configuration.md#Merge-precedence)).
`effective_config.yaml` of each run records what took effect.

Case-scope keys are `matpower_import.*`, `cgmes_import.*`,
`powsybl_import.*`, `power_flow.*`, `state_estimation.*` and
`short_circuit.*`; `output.*`, `benchmark.*`, `runtime.*`, `webui.*` and
`matpower_export.*` stay in the machine configuration
(`scf_is_case_config_key(key)` answers per key). A key outside the case
scope is refused by name on read. The in-file block `sparlectra.config` is
deprecated: the writer does not emit it, the reader applies it below the
case configuration file with a warning naming the new place. The Web UI
form is seeded from the case file and posts every option as an explicit
override ([Case-specific settings](webui_reference.md#Case-specific-settings)).

## The reference: single, distributed, or a source

`roles.slack` states which reference the case was built with:

| `mode` | Meaning |
|---|---|
| `single` | one ideal slack per island, `nodes` lists them |
| `distributed` | the reference bus keeps the angle, the island's active-power imbalance is shared; `participation` carries one factor per generator |
| `source` | no slack flag at all: the in-service PGM `source` components are the reference (this is what a file written by power-grid-model looks like) |

```json
"roles": {
  "slack": {
    "mode": "distributed",
    "nodes": [1],
    "normalize": true,
    "participation": [{"object": 23, "factor": 2.0}]
  }
}
```

The factors are the network's own participation factors (MATPOWER `APF`,
CGMES `normalPF`). They travel with the model but do not switch the
feature on; the case has to ask for it in its configuration (these keys,
in `<stem>.config.yaml` or in the deprecated in-file block shown here):

```json
"config": {
  "power_flow.distributed_slack.enabled": true,
  "power_flow.distributed_slack.p_mode": "imported"
}
```

`imported` reads exactly those factors; `explicit` weights
(`power_flow.distributed_slack.weights`) are a run setting and stay in the
YAML configuration ([Power-Flow Configuration](powerflow_configuration.md)).

`extra.<appliance>.reference_pri` flags the appliance that holds the
reference. `extra.<appliance>.reference_priority` says which one takes it
over when an island has none or the slack is lost in N-1: an integer with
the CGMES semantics (1 is the strongest, larger is weaker), written only
where one is stated; a file without it reads with no priority on any unit
([Reference priority](slack_vs_source.md#Reference-priority)).

Generator costs from a MATPOWER case (`mpc.gencost`) sit on their machine
as `extra.<machine>.gencost`: `{"p": [...], "q": [...]}`, each a MATPOWER
cost row (`MODEL`, `STARTUP`, `SHUTDOWN`, `NCOST`, then points or
coefficients), `q` optional. No calculation reads them; a MATPOWER export
writes the block back unchanged ([Generator costs](@ref matpower_gencost)).

### [Dynamic unit data](@id scf_dynamics)

Since revision 1.3 a machine can carry its [`DynamicData`](@ref) as
`extra.<machine>.dynamics`, one key per field under the field name. A
field at `nothing` is left out, the fields with a numeric default are
written with their value, `unit_type` is `"synchronous"`,
`"grid_following"` or `"grid_forming"`:

```json
"dynamics": {"inertia_s": 4.25, "damping_pu": 0.0, "droop_pu": 0.05,
             "deadband_hz": 0.0, "primary_reserve_mw": 0.0,
             "agc_participation": 0.0, "secondary_reserve_mw": 0.0,
             "unit_type": "synchronous", "black_start": false}
```

A machine without dynamic data has no block, and a case without any writes
the bytes it wrote before. On read a missing key takes the default of
[Dynamic Unit Data](@ref dynamic_data), a missing `unit_type` follows the
inertia. An unknown key, a value of the wrong type or a value against the
rules of its field fails with an `ArgumentError` naming the machine and the
key. Write and read keep every value bit for bit. No calculation reads the
block.

## Importers and the typed case

Every importer builds its network directly (an SCF file through
`build_net`); converting a case into an SCF file is always an explicit
export, so imported values stay bit-exact.

`SCFCase` is the in-memory form of one document (the PGM `data` section
as typed component vectors, the `sparlectra` block as sub-objects).
`read_scf_json(path)` parses and validates a file into it,
`write_scf_json(case, path)` serializes it to the canonical bytes,
`build_net(case)` constructs the network and applies the run
configuration's net parameters exactly once, `net_to_scfcase(net)` is the
export direction (same keywords as `exportSCF`). `importSCF` and
`exportSCF` wrap these four; `scf_to_net`/`net_to_scf` keep the dict entry
points.

## The round trip

A network written by Sparlectra, read back and written again gives two
byte-identical files: bus and branch order, the source system's bus
numbers, tap ratios, shunt susceptances with their sign, controllers,
measurements and the power-flow result come back bit for bit, and a run
on the re-read network gives the same power-flow iterations and voltages
and the same state-estimation objective, degrees of freedom and state. A
network imported from MATPOWER, CGMES or DTF keeps its electrical model
and solution, not every detail of the source file (CGMES mRIDs are
recorded, MATPOWER comments and column layout are not).

**Notes**

- Per-unit values can land one ulp away when no exact SI preimage exists;
  the Y-bus deviation is at the 1e-16 level.
- Tap regulation bands snap to the step grid: the operating tap ratio is
  exact, a band edge can move by less than one step. A neutral ratio
  outside its own band is kept and named in a warning on load. On
  `case13659pegase` one phase-tap band moves by one step on the second
  write; the third write equals the second.
- An unlimited branch rating (MATPOWER `rateA = 0`, arriving as `Inf`) is
  absent from the PGM dataset (no `i_n` is no limit there) and travels as
  `"inf"` in the namespaced block; any other non-finite number is a hard
  error naming the key.
- An unlimited reactive band is stated by absence: `max_q_mvar` and
  `min_q_mvar` are written only when finite, in both blocks; a missing and
  an infinite limit give the same power flow.
- Impedances are written from the base (physical) values; exporting a
  FACTS-compensated operating point is a hard error naming
  `restoreBaseImpedances!`.

The operating point travels in `extra`: the reactive band of every
machine, its voltage setpoint and whether it regulates (`max_q_mvar`,
`min_q_mvar`, `vm_pu`, `regulated`), because a machine written as a PGM
`source` has no `voltage_regulator` row. So does a voltage-dependent
controller: `extra.<machine>.qu_control` and `pu_control` carry the
characteristic's points in per unit (voltage, power on `meta.s_base`),
its interpolation mode and its limits in MVAr / MW:

```json
"qu_control": {
  "points": [[0.95, 0.30], [1.00, 0.00], [1.05, -0.20]],
  "interpolation": "spline",
  "qmin_mvar": -50.0,
  "qmax_mvar": 50.0
}
```

`interpolation` is `linear`, `spline` or `polynomial`
([Voltage-dependent control](voltage_dependent_control.md)); the reader
validates at least two points, strictly increasing voltages, a known mode
and `qmin_mvar` not above `qmax_mvar`. The MATPOWER converter's constant
controllers keep their short form, the flag `pq_gen_controller`; an
explicit `qu_control`/`pu_control` object wins over it. A setpoint alone
does not promote a machine: the reader restores the `regulated` flag the
file states.

Ids are deterministic: the writer walks the components in a fixed type
order and, inside a type, in lexicographic name order, so the same
network, or a case rebuilt from its source, gets the same ids. `extra`
records each element's internal index, and `measurements.rows` lists the
measurements in their internal order, because internal order fixes
generated ids and summation order. The writer emits the SI value that
reads back into exactly the per-unit value it started from (both sides
derive the impedance base from `u_rated` and `s_base` with the same
formula); `measurements.rows` carries the sigma of each measurement (one
PGM sensor sigma cannot express a P and a Q sigma) and a value only where
the SI value would not read back exactly.
Floats use the shortest round-trip representation; the file is
pretty-printed with one key per line (lists of plain scalars inline), so
a diff points at the changed field.

### Fields that are not written

A field whose value equals its documented default is left out and the
reader puts it back; this applies to the namespaced block only, `data`
keeps every attribute PGM requires.

| Field | Default | Meaning when absent |
|---|---|---|
| `extra.<branch>.kind` | `"BranchC"` | an ordinary line |
| `extra.<machine>.vm_pu` | absent | the machine does not regulate (a regulating one keeps its setpoint in its `voltage_regulator` row, the reference in `source.u_ref`) |
| `extra.<machine>.max_q_mvar`, `min_q_mvar` | absent | no reactive limit |
| `extra.<machine>.max_p_mw`, `min_p_mw` | absent | no active limit |
| `extra.<machine>.qu_control.interpolation`, `pu_control.interpolation` | `"linear"` | piecewise linear characteristic |
| `extra.<machine>.qu_control.qmin_mvar`, `qmax_mvar`, `pu_control.pmin_mw`, `pmax_mw` | absent | no limit on the controlled value |
| `extra.<machine>.regulated`, `apu_node` | `false` | |
| `extra.<machine>.reference_priority` | `0` | no reference priority stated |
| `extra.<machine>.gencost` | absent | no generator cost data |
| `extra.<machine>.dynamics` | absent | no dynamic data (`ProSumer.dynamics = nothing`) |
| `extra.<branch>.meta.sn_mva` | absent | no rating recorded |
| `components.tap_changer[].tap_est_mode` | `"none"` | the tap is not released for estimation |
| `components.tap_changer[].controllers[].pos` | absent | the step was not transferred: the changer sits at its middle position, the neutral position 0, and the import names it in an info line |
| `tap_changer_models.<id>.ratio.step`, `.phase.step`, `.estimation.changers[].step` | absent | the changer sits at its middle position, `neutral_step` (`neutral` in `estimation`), named in an info line |
| `components.tap_changer[].controllers[].step_assumed`, `tap_changer_models.<id>.ratio.step_assumed`, `.phase.step_assumed`, `.estimation.changers[].step_assumed` | `false` | `pos` or `step` is a step of the source. The writer puts `true` on a changer whose step was not transferred: `pos` or `step` then holds the assumed middle position, and the reader marks the changer again and names it in an info line, so a written middle position never reads back as a real step. A tap controller or a QSTS carry that moves the changer takes the mark off, the step it reaches is a step |
| `components.tap_changer[].controllers[].grid` | `"fraction"` | ratio positions on the fraction grid `ratio_base / (1 + pos * step)`; `"linear"` for a CGMES changer |
| `components.tap_changer[].controllers[].declared` | `true` | the changer is nameplate data; a declared phase changer joins the mass release of the tap estimation. The writer puts `false` on a controller no nameplate row declared (the default ratio or phase grid of a transformer), so the re-read keeps that state |
| `components.transformer3w[].status` | `1` | in service |
| `roles.slack.normalize` | `true` | participation factors are normalized |

`area` and `zone` are not in this list (their constructor value is
"unknown", not `1`) and are written only when set.

## File layout

The root object is the PGM root object with exactly one added key:

```json
{
  "version": "1.0",
  "type": "input",
  "is_batch": false,
  "attributes": {},
  "data": { },
  "sparlectra": { }
}
```

`data` is a valid PGM input dataset in SI units (V, W, var, VA, ohm,
siemens, farad, ampere, radian); `sparlectra` holds everything PGM has no
model for. A plain PGM reader still gets a computable network; see
[Interoperability](#Interoperability) for what it loses.

### The PGM subset in `data`

| Component | Sparlectra element | Notes |
|---|---|---|
| `node` | bus | `u_rated` is the rated line-to-line voltage in V |
| `line` | line, cable | `r1`, `x1` in ohm, `c1` in farad, `tan1`, `i_n` in A; a shunt conductance rides in `tan1` (= g/b), and the one value that encoding cannot spell, g without capacitance (b = 0), is written as a direct `g1` in siemens, an SCF extension the strict PGM writer strips by name |
| `generic_branch` | transformer, phase shifter | `r1`, `x1`, `g1`, `b1` referenced to the **to** side, `k` the live off-nominal ratio, `theta` the live shift in radians, `sn` in VA |
| `link` | busbar coupler | zero-impedance connection |
| `source` | external network injection | `u_ref`, optional `u_ref_angle` |
| `sym_load`, `sym_gen` | load, generator | `p_specified` in W, `q_specified` in var |
| `shunt` | shunt | `g1`, `b1` in siemens |
| `voltage_regulator` | PV control of a generator | `u_ref` in per-unit, `q_min`/`q_max` in var (PGM 1.13 or newer) |
| `sym_voltage_sensor`, `sym_power_sensor`, `sym_current_sensor` | state-estimation measurements | values and sigmas in SI |

The Sparlectra Y-bus stamps the tap on the **from** side, so
`r_pu`/`x_pu` are referenced to the **to** side, like PGM's
`generic_branch`.

### [Unit system](@id scf_units)

The numbers in `data` are base SI units:

| Field group | Fields | Unit |
|---|---|---|
| line impedance | `line.r1`, `line.x1` | ohm |
| line capacitance | `line.c1` | farad |
| line loss factor | `line.tan1` (conductance over susceptance) | 1 |
| line conductance without capacitance | `line.g1` (optional) | siemens |
| conductance, susceptance | `g1`, `b1` of `generic_branch` and `shunt` | siemens |
| power | `p_specified`, `q_specified` | W, var |
| voltage base | `node.u_rated` | V (line to line) |
| current rating | `line.i_n` | A |
| power rating | `generic_branch.sn`, `meta.s_base` | VA |
| angle | `theta` | rad |

`sparlectra.meta.units` declares in which unit system the branch and
shunt parameters are stated. Two values exist:

- `"si"`: the values in `data` are the ones the network is built from.
  This is what a file without the key means, so every file written before
  the key existed reads unchanged.
- `"pu"`: the impedances and admittances are stated per unit in a `pu`
  block of the component's `extra` entry, and the network is built from
  them. The `data` section stays a complete SI view for PGM readers.

A `pu` block holds `r_pu`, `x_pu`, `b_pu`, `g_pu` for a `line` or
`generic_branch` (total shunt values, on the base of the **to** node) and
`g_pu`, `b_pu` for a `shunt`. The base is the one the file already
carries: `meta.s_base` and the `u_rated` of the node; a `"pu"` file without
`meta.s_base` or with a node without a positive `u_rated` is an
`ArgumentError` naming the field, and so is a branch or shunt without its
`pu` block. Powers, voltages, ratings and angles stay in SI in both
systems. Powers and voltages convert by a power of ten and read back
exactly. Ratings and phase-shift angles do not always: a line rating
travels as `i_n` in A, a transformer rating as `sn` in VA, a shift as
`theta` in rad, and a rating like 200 MVA can read back as
200.00000000000003 MVA, a shift like -0.06138 degrees as
-0.06137999999999999 degrees. The difference is in the last bit only.

The declaration exists because a per-unit value converted to ohm and back
does not always return the same number: in about one value in ten no SI
value maps back exactly. A case that comes from per-unit data therefore
keeps it per unit, and the only conversion happens when the network is
built.

An SI file (excerpt, node 2 rated 110 kV, 100 MVA base, 50 Hz):

```json
"data": { "line": [ { "id": 7, "from_node": 1, "to_node": 2, "r1": 1.2342, "x1": 3.0976, "c1": 1.1312e-7 } ] },
"sparlectra": { "meta": { "s_base": 100000000.0 } }
```

The same line in a `"pu"` file:

```json
"data": { "line": [ { "id": 7, "from_node": 1, "to_node": 2, "r1": 1.2342, "x1": 3.0976, "c1": 1.1312e-7 } ] },
"sparlectra": {
  "meta": { "s_base": 100000000.0, "units": "pu" },
  "extra": { "7": { "name": "ACL_1", "pu": { "r_pu": 0.0102, "x_pu": 0.0256, "b_pu": 0.0043, "g_pu": 0.0 } } }
}
```

`exportSCF` writes SI by default and `"pu"` with `units = :pu`; the
MATPOWER adapter (`convert_case`) and the Web UI export of a MATPOWER case
write `"pu"`. A strict PGM file (`strict_pgm = true`) has no `sparlectra`
block and is always SI.

### The `sparlectra` block

| Key | Content |
|---|---|
| `format_version` | schema version of the block (`1.3`) |
| `meta` | case name, `s_base` in VA, `units` (`"si"` or `"pu"`, see [Unit system](@ref scf_units)), `f_nom` (system frequency in Hz, optional; the writer always states the frequency of the net; on read it sets the frequency of the net and wins over `model.nominal_frequency_hz`, which applies only when `f_nom` is absent, default 50 Hz; a value that is not a positive number is an `ArgumentError` naming `meta.f_nom`; see [System frequency](@ref system_frequency)), source format and reference, producer, `intended_calculations` as the file's contract, and for a CGMES delivery `cgmes_ids`, the mRID registry of the import keyed by structural identity, which the reader restores so a re-export keeps the mRIDs |
| `roles` | what PGM cannot express: `slack.mode` (`single`, `distributed`, `source`) with its nodes or participation factors, auxiliary nodes, isolated nodes |
| `extra` | per-component annotations keyed by id: `name` (the reference name every other section resolves against), `external_id` (CGMES mRID, MATPOWER index), bus/branch index, node type, ratings, and in a `"pu"` file the `pu` block of each branch and shunt |
| `components.tap_changer` | the tap changers PGM has no model for: per controller the step grid, position, band, and nameplate angle, plus the neutral ratio and angle. Ratio positions count on the fraction grid `ratio_base / (1 + pos * step)` (nameplate data), or with `grid: "linear"` on the grid of a CGMES changer, `ratio_base + pos * step`. `step_assumed: true` marks a position that was not transferred and is the assumed middle position |
| `components.transformer3w` | three-winding transformer: either the grouping of three existing `generic_branch` legs (their star node, and which end is hv/mv/lv; the electrical values stay in the legs and are never duplicated) or a nameplate that describes the transformer directly, see below |
| `components.sc_source` | IEC 60909 source data (external network injections, machines) that PGM's `source` cannot carry beyond `sk`/`rx_ratio` |
| `components.shunt_state` | what PGM's `g1`/`b1` cannot say: whether a shunt is in service, whether it is a voltage-dependent injection rather than an admittance, and whether its susceptance is released as a state-estimation state. Only deviating shunts produce a row |
| `transformer_types` | named transformer nameplates in CGMES `PowerTransformerEnd` vocabulary, referenced by `components.transformer3w` entries |
| `branch_shunt_split` | per branch id, only for a branch whose two terminal shunt arms differ: `g_from_pu`, `b_from_pu`, `g_to_pu`, `b_to_pu` in pu on the branch base. The `data` section keeps the PGM totals (`c1`, `b1`); on load the row overrides the symmetric half. Absent block: symmetric |
| `tap_changer_models` | per branch id, the typed tap-changer models of the from-side winding: `ratio` (the `PowerTransformerTaps` inputs: `vn_kv`, `step`, `low_step`, `high_step`, `neutral_step`, `voltage_increment_kv`, `neutral_u_kv`, `convention`, `step_assumed`), `phase` (`kind` `symmetrical`, `asymmetrical` or `tabular`, the steps, `voltage_step_increment`, `step_phase_shift_increment`, `winding_connection_angle_deg`, `x_min`, `x_max`, `convention`, `step_assumed`, and `table` rows `step`, `ratio`, `angle_deg`, `x_pu` for `tabular`), `tap_changer_model` when the impedance correction is on, and `base_ratio`, `base_shift_deg`, the winding ratio and shift the models move from (`k` and `theta` of the data row carry the steps). On load the models are restored onto the winding and the resolver derives ratio, shift and grid from them; an entry without `base_ratio` or `base_shift_deg` takes `k` or `theta` as the base; `components.tap_changer` rows still describe the grid for readers that know no model. The older spelling `extra[<id>].phase_taps` is still read, with a warning. The optional key `estimation` carries the tap changer set the state estimation uses for the branch, as an import that knows the changers stored it (CGMES): `origin` (`cgmes` or `powsybl`, any other value is refused on write and read), `ratio_n`, `shift_n_deg` and the `changers` with `name`, `mrid`, `kind`, `role`, `side`, `neutral`, `low`, `high`, `step`, `impedance` and, where set, `additive`, `du`, `alpha_deg`, `dphi_deg`, `step_abs`, `outside_du`, `rx_correction`, `step_assumed`, `x_curve` (`x_min`, `x_max`, `a_max_deg`, `form`) and `table` rows (`step`, `ratio`, `angle_deg`, `r`, `x`, `g`, `b`). An empty `changers` list means the transformer has no tap changer. An entry with only this key runs no resolver; a file without it reads as before |
| `components.controllers` | FACTS and regulation, in the declarative `control.controllers` schema verbatim: same type names, same keyword names. The reader hands the entries to `applyConfiguredControllers!`, so there is one construction path. The shunt of a shunt voltage controller is not written as a `shunt` row; its entry carries the present susceptance as `bs_start_mvar` and, where they deviate, the shunt's `status`, `model` and `estimate` in a `shunt_state` object with the keys of a `components.shunt_state` row |
| `contingencies` | removed in 0.40.0. A file that still carries the block fails on load with a message naming the `scenarios` block and the Web UI scenario editor |
| `scenarios` | the scenario model: `mode` (`explicit`, `n1_branches`, `n1_generators`, `n1_all`), `exclusions`, `events` (correlated event classes generated next to the list: `double_circuits`, `busbars`, `substations`, `links`, `name_groups`, `pairs`), and the ordered scenario list, each scenario `name`, `weight`, an optional `kind` (written when it is not `explicit`) and its `ops` (`op` = `status`/`set`/`scale`, `target` component class, `id` the SCF component id, plus `value`, `field`, `factor` as the op needs). The N-1 modes and the event classes expand through the same generators `runContingencies!` uses ([N-1 Contingency Analysis](contingency.md)); a scenario of status-0 ops only runs as one outage event |
| `short_circuit` | IEC 60909 study definition: `case` (`max`/`min`), optional `c_factor`, `sweep` (`all_buses`/`explicit`) and, for an explicit sweep, the `buses` list of node ids |
| `measurements.rows` | one row per measured value with the semantics of the [measurement table](state_estimation_measurement_table.md), in the order of the measurement vector; see [Measurements](#Measurements) |
| `config` | case-specific dotted configuration keys (same allowlist the API and Web UI use); deprecated in favour of the case configuration file |
| `start_state` | opt-in bus voltages as start values, never as results |

The study blocks are validated on load (an unknown mode, an out-of-band
`c_factor` or an empty outage list fails on read, not at the start of a
long sweep).

The two blocks for one transformer branch with id 12:

```json
"branch_shunt_split": {
  "12": {"g_from_pu": 0.0004, "b_from_pu": -0.012, "g_to_pu": 0.0, "b_to_pu": 0.0}
},
"tap_changer_models": {
  "12": {
    "phase": {"kind": "asymmetrical", "step": 3, "low_step": -10, "high_step": 10,
              "neutral_step": 0, "voltage_step_increment": 0.01,
              "winding_connection_angle_deg": 90.0, "convention": "reciprocal_from_side"}
  }
}
```

The `data` section stays a valid power-grid-model dataset: the split and
the models live only in the `sparlectra` block, and `exportSCF(strict_pgm
= true)` names both among what the plain dataset does not carry.

### Names and source ids

`extra[<id>].name` is the reference name that measurements, controllers,
contingencies and reports resolve against (Sparlectra's bus dictionary
key); a differing internal component name travels as `component_name`,
the source-system id as `external_id` (for a CGMES delivery the mRID,
marked `source_id_kind: cgmes_mrid`, otherwise the internal component
id). With `"name": "Ostheim_110", "external_id": "H1"` the import uses
the place name everywhere, restores `external_id` and `component_name`
onto the network, and a re-export preserves both byte for byte, as the
[shipped demo cases](demo_cases.md) do. Branches, links, appliances and
shunts take `name` and `external_id` as their component name and id in
the same way, so results, measurement tables and contingency rows show
the names of the file (`Lindhof_110_20`, not the generated
`B_2WT_110_4_5`). A `generic_branch` without a `tap_changer` row has no
tap changer.

### Three-winding transformers without writing legs

Instead of three hand-written `generic_branch` legs plus a star node, a
`transformer3w` entry may describe the transformer as CGMES does: one
`PowerTransformerEnd` per winding with `rated_u`, `rated_s`, `r`, `x` and
`b`. The reader builds the star equivalent through `add3WTPiModelTrafo!`
and creates the auxiliary star node itself.

```json
"sparlectra": {
  "transformer_types": {
    "cgmes_220_110_10_100MVA": {
      "ends": [
        {"role": "hv", "rated_u": 220000.0, "rated_s": 100000000.0, "r": 0.6, "x": 30.0, "b": 0.0},
        {"role": "mv", "rated_u": 110000.0, "rated_s": 100000000.0, "r": 0.4, "x": 15.0, "b": 0.0},
        {"role": "lv", "rated_u": 10500.0,  "rated_s": 40000000.0,  "r": 0.2, "x": 8.0,  "b": 0.0}
      ]
    }
  },
  "components": {
    "transformer3w": [
      {"id": 7, "type": "cgmes_220_110_10_100MVA",
       "nodes": [{"role": "hv", "node": 1}, {"role": "mv", "node": 2}, {"role": "lv", "node": 3}]}
    ]
  }
}
```

A `type` entry is a named nameplate shared by transformers of one design;
an entry can also carry its `nameplate` inline. An entry either points at
legs or describes them; mixing the two fails on load.

#### Which CIM attributes the nameplate speaks

One nameplate end is one CIM `PowerTransformerEnd`:

| Nameplate key | CIM attribute | Unit here | |
|---|---|---|---|
| `role` | the winding's voltage level (CIM orders ends by `PowerTransformerEnd.endNumber`) | `hv`, `mv`, `lv` | required |
| `node` | the end's `Terminal` to `TopologicalNode` | node id | required, in the end or in `nodes` |
| `rated_u` | `PowerTransformerEnd.ratedU` | V | required |
| `rated_s` | `PowerTransformerEnd.ratedS` | VA | required |
| `r` | `PowerTransformerEnd.r` | ohm | required |
| `x` | `PowerTransformerEnd.x` | ohm | required |
| `b` | `PowerTransformerEnd.b` | siemens | optional, default 0 |
| `status` (on the group, not the end) | in service | 0 or 1 | optional, default 1 |

All three ends state `r`, `x` and `b` on the HV base, the star-equivalent
form the CGMES import produces; the reader converts to per unit with the
highest `rated_u` as the voltage base. Not taken from CIM: the magnetising
conductance `g` (the losses sit in `r`), `phaseAngleClock` and
`connectionKind` (the vector group, not evaluated by the balanced
positive-sequence model), and the tap changers. A `RatioTapChanger` or
`PhaseTapChanger` lives in its own `components.tap_changer` entry,
referenced by the group's `tap_changer`; controller `index` 1 is the
ratio changer, `index` 2 the phase changer, each with step size, live
position and position band. Any other nameplate key fails on load. The
nameplate is an input shorthand: an export always writes the explicit
legs plus the grouping, the form that round-trips.

#### What the star equivalent is called

The reader generates the names of the star point and the three legs, and
whatever addresses them later (a contingency case, a measurement, a tap
block) uses exactly these:

| Object | Name | Built from |
|---|---|---|
| star node | `Aux3WT_<hv>_<mv>_<lv>` | the three terminal bus names |
| leg | `B_2WT_<Vn>_<star>_<terminal>` | rated voltage in kV, then the two bus indices |

For the nameplate above, with terminals `HV`, `MV` and `LV`, the export
writes the star node `Aux3WT_HV_MV_LV` (listed under `roles.aux_nodes`)
and the legs `B_2WT_220_4_1`, `B_2WT_220_4_2` and `B_2WT_220_4_3`; the
grouping states which leg is which end:

```json
"components": {
  "transformer3w": [
    {"id": 13, "star_node": 1, "leg_direction": "star_to_terminal", "tap_changer": 10,
     "ends": [
       {"role": "hv", "branch": 5, "terminal_node": 2, "u_rated": 220000.0, "sn": 100000000.0},
       {"role": "mv", "branch": 6, "terminal_node": 4, "u_rated": 110000.0, "sn": 100000000.0},
       {"role": "lv", "branch": 7, "terminal_node": 3, "u_rated": 10500.0,  "sn": 40000000.0}
     ]}
  ]
}
```

`leg_direction` says that every leg runs from the star node to its
terminal (`terminal_to_star` for the reverse, `mixed` when the legs
differ, as for a CGMES leg with a tap changer), and the roles follow the
terminal voltages, highest first. A
CGMES delivery keeps its own names from the `PowerTransformer` and its
ends; the generated names appear only for a nameplate-built transformer.

### Running the studies

A study run on a case file executes what the file defines and applies the
file's own configuration.

| Study | File block | Entry | Result |
|---|---|---|---|
| N-1 and scenarios | `scenarios` (a scenario may combine several status, setpoint and scaling operations, so a double outage is one scenario) | `runScenarios!`; the service with `scenario_source = file_block`; the Web UI's scenario editor writes the block | the per-scenario `weight` is honored, `run.log` records where the case list came from |
| State estimation | `measurements` | the run needs no separate CSV; an explicitly picked CSV wins | `measurements.csv` is still written as the run's artifact |
| Short circuit | `short_circuit`: `case` (`max` or `min`), optional `c_factor`, the sweep | `sweep: explicit` takes a `buses` list of node ids resolved through `extra`; without an explicit sweep PGM's own `data.fault` rows name the faulted node | both cases are always evaluated (`short_circuit_max.csv`, `short_circuit_min.csv`); `case` only decides where the headline numbers come from |

**Notes**

- Scenarios are a feature of this format. Converting a case to MATPOWER,
  CGMES or PowSyBl does not carry the `scenarios` block.
- An unknown bus id in an explicit sweep fails the run; a `fault_type`
  other than `three_phase` or a non-zero `r_f`/`x_f` is refused on load
  (balanced bolted fault only). The file's `c_factor` applies at the
  configuration default and loses against an explicit
  `short_circuit.c_factor`.
- Source data is in `components.sc_source` (empty fields are written as
  `null`). A machine without a source row enters with the default
  reactance, the rows it feeds are flagged and the run ends with status
  `warning`; a file whose sources carry neither data nor a rated power is
  refused ([Short-Circuit Analysis](@ref short_circuit_source_data)). The
  reference unit of a bus with a feeder row is that feeder, not counted
  twice. An export writes the rows of the file only.

## Source formats

The writer serializes a `Net`, so every imported case exports the same
way, with the payload its importer put into the network.

| Source | What the file carries |
|---|---|
| MATPOWER | branch and shunt parameters per unit (`meta.units = "pu"`) when written by the Web UI export or `convert_case`, bus and branch names, the `mpc.sparlectra` extensions (tap-changer nameplates as `components.tap_changer`, busbar couplers as `data.link`), branch kinds, ratings |
| CGMES | mRIDs as `external_id`, the three-winding star grouping with its auxiliary node, IEC 60909 source data, the reference names of the topological nodes |
| DTF | station reference names, transformer tap data, ratings |

Three-winding transformers need no conversion: the star equivalent (three
legs plus an auxiliary star node) is what PGM documents for
`generic_branch`; the grouping block only records which legs belong
together.

## Measurements

The measurements of a case travel twice: as PGM sensors in `data`, in SI,
for a PGM reader, and as `sparlectra.measurements.rows` with the semantics
of the [measurement table](state_estimation_measurement_table.md), which
Sparlectra reads.

PGM sensors: a bus voltage (plus its angle if present) becomes one
`sym_voltage_sensor`, a bus balance P/Q pair one `sym_power_sensor` with
`measured_terminal_type = node`, a branch-end P/Q pair one branch power
sensor, a current magnitude/angle pair one `sym_current_sensor` (magnitude
and sigma in A); a missing half of a pair is `null`. The shunt bay has no
PGM sensor.

`measurements.rows` holds one entry per measured value, in the order of
the measurement vector (the WLS objective sums in that order):

| Key | Content |
|---|---|
| `type` | `U`, `P`, `Q`, `Va` at a node; `P`, `Q`, `I`, `Ia` at a terminal |
| `node` | node id: `U`, `Va`, and `P`/`Q` as the bus balance (injection positive) |
| `terminal` | `{"branch": id, "side": "from"}` (or `"to"`) for a branch end; `{"shunt_bay": node id}` for the shunt bay of a node: `Q` is the sum over its shunts, `I` and `Ia` the bay current |
| `sensor` | the PGM sensor that carries the reading; absent for a shunt-bay row |
| `id` | the Sparlectra measurement id; absent when it equals the default id of its location, which the reader restores: `Vm_bus_<k>`, `Pflow_branch_<k>_from` and so on, where `k` is the bus or branch number inside Sparlectra (the index recorded in `extra`), not the node or branch id of the file |
| `value` | in kV (`U`), MW, MVAr, A or degrees; written only where the sensor value would not read back bit for bit, and always for a shunt-bay row |
| `sigma` | in the same unit |
| `valid` | `false` for an inactive measurement; absent means `true` |
| `zero_injection` | `true` for a zero-injection row (a `P` or `Q` balance at a node, the source token `zero_injection` of the table); absent means `false`. From 1.2 on the id decides nothing; a 1.0 file marks such a row only by the id `ZI_PINJ_bus_<k>` / `ZI_QINJ_bus_<k>` with value 0, and the reader keeps that rule for 1.0 |

Units carry no base, as in the table: only `U` is divided by the rated
voltage of the node (`u_rated` / 1000) on reading. A `U` value or sigma
that no 64-bit float divides back to its per-unit value is written with up
to 18 significant digits. The reader keeps the text of `value` and `sigma`
in `measurements.rows` only, and divides such a `U` text at 256 bits, as
the table does; every other number of the file is read as a 64-bit float.
A row without `value` or `sigma` takes it from its sensor; a sensor no row
names is read after the rows with the id `Vm_<node id>` (`Pinj_<node id>`,
`Imag_<branch id>`, ... with the ids of the file). A bus link flow has no row:
the file has no sensor for a bus link, and the export names such rows in a
warning (the measurement table carries them on the terminals of the link).

Two places where the block differs from the table: it is keyed by the ids
of the file instead of names or mRIDs, and the shunt bay is keyed by its
node (the measurement is the sum over the shunts of the bus) and also
carries a current angle `Ia`, which the table has no place for. Provenance,
noise, file formats and generators:
[Measurements](state_estimation_measurements.md).

## Interoperability

| Level | What you get |
|---|---|
| **PGM-readable** | the `data` section alone is a valid PGM input dataset |
| **Sparlectra-complete** | the full file: slack strategy, tap changers, names, measurement detail, configuration |
| **Lossy toward PGM** | a PGM reader loses the tap changer detail (the live ratio survives in `k`/`theta`), the slack strategy (it sees `source` components), the names, and the configuration |
| **Strict PGM** | `exportSCF(net; strict_pgm = true)`, the "Export as plain PGM" button, or `sparlectra export --pgm` writes the dataset alone, for a consumer that reads PGM and nothing else |

Strict mode drops the namespaced block, names what it dropped in a
warning (always `f_nom`: PGM has no frequency field, the line `c1` is
written with the frequency of the net, and a reader of the file takes
`model.nominal_frequency_hz`), and rewrites a slack generator as a PGM
`source` (PGM's reference is a source; a generator left in place would
make the node inject twice).

A file written by power-grid-model has no `sparlectra` block; its
in-service `source` components become the reference, and since PGM's
`source` is Sparlectra's external network injection, a case that arrives
with a source leaves with a source. Its `sk` and `rx_ratio` are read as
feeder data, so a short circuit runs on the file as it stands, and
`power_flow.external_grid` with `source: auto`
([Power-Flow Configuration](powerflow_configuration.md)) moves the
reference behind $z = U_n^2 / S_k''$ and reproduces what PGM computes; the
default stays the ideal slack. An export writes `sk`/`rx_ratio` back
whenever the network carries them. Its `sym_current_sensor` rows are read
in A, as PGM defines them (`i_measured = 50.0` is 50 A); before 0.50.0
Sparlectra divided them by 1000.

## Time series

A time series on a case file is a zip with the unchanged case, a manifest
and one table of patch operations for all steps; see
[SCF time series](@ref scf_time_series). The case file has no series
block.

## Version and compatibility

`sparlectra.format_version` is always written with the current revision.
A reader accepts its own and the older revisions listed below. Backward
compatibility across revisions is not promised: a breaking change is named
in the changelog; a case file is an exchange and archive format for the
version that wrote it.

| SCF revision | Reader behaviour |
|---|---|
| `1.3` | the current revision: accepted, written by every export. Adds the optional block `extra.<machine>.dynamics` ([Dynamic unit data](@ref scf_dynamics)) |
| `1.2` | accepted and read as written: it differs from 1.3 only by the missing `dynamics` block. Current sensors carry A, `measurements.rows` follows the measurement table |
| `1.0` | accepted; its `i_measured` and `i_sigma` are 1000 times too large and are divided by 1000 on reading, giving the value that was exported; its `measurements.rows` (one entry per sensor with `sparlectra_ids`, `sigmas`, `positions`) and `measurements.shunt_bay_rows` are read as before |
| `1.1` | the pre-release name of `1.0`: read as `1.0` with an info message asking for a re-export |
| any other value | refused by name, stating what it found, what it expected, and that the case has to be re-exported |
| field absent | refused the same way; a dataset written by power-grid-model has no namespaced block and skips the check |

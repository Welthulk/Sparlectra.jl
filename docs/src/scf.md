# Sparlectra Case Format (SCF)

A single self-describing case file that a power-grid-model (PGM)
installation can read directly, while carrying everything Sparlectra needs
beyond the PGM model: slack strategy, the tap-changer cascade, component
names and source ids, and state-estimation measurements. The case's own
settings travel in a case configuration file next to it
(`<stem>.config.yaml`).

A `.scf.json` file is a **first-class case format**: it appears in the Web
UI case selector, runs through the framework and the service like a
MATPOWER or CGMES case, and brings its own configuration with it through
its case configuration file.

## Writing a case file

```julia
using Sparlectra

net = createNetFromMatPowerFile(filename = "case118.m")
exportSCF(net; file = "case118.scf.json",
          source_format = "matpower", source_reference = "case118.m",
          intended_calculations = ["power_flow", "state_estimation"])
write_case_config("case118.scf.json", Dict("power_flow.mode" => "auto"))
```

In the Web UI, the **Export as SCF case file** button under the case
selection on the PowerFlow page writes `<case>.scf.json` into the case
directory. It imports the selected case through the normal run import path,
picks up a `<case>.measurements.csv` sitting next to it, and stores the
form's own configuration overrides in the case configuration file next to
the export, so the case ships with the settings it was meant to run with.

## Reading a case file

```julia
net = importSCF("case118.scf.json")
runpf!(net, 30, 1e-8, 0)
```

The reader builds the network through the same public constructors every
other importer uses (`addBus!`, `addPIModelACLine!`, `addPIModelTrafo!`,
`addProsumer!`, `addShunt!`, `addLink!`, and the shared tap-nameplate
application) and finishes with `validate!`, so an SCF case behaves exactly
like a case from any other source. Validation runs in three stages:

1. **Schema**: unknown top-level keys, unknown keys inside `sparlectra`, and
   unsupported component types are hard errors naming the offender. Silent
   acceptance is the failure mode this format exists to prevent.
2. **Reference integrity**: every id is unique across the whole file, and
   every `from_node`, `node`, `measured_object`, `regulated_object`,
   `branch`, and `star_node` resolves to an object of the right kind.
3. **Model plausibility**: the checks `validate!` already performs.

`load_case_config(file)` returns the dotted configuration keys of a case's
case configuration file, `scf_case_config(file)` those of the deprecated
in-file block, and `scf_case_studies(file)` its study definitions, both without
building the network.

### Getting the file out of the Web UI

The export buttons write next to the case, because that is where a run
finds the file again. To take it to another machine, the page offers it
right after the export ("Download <name>"), and the actions row carries a
**Download selected case** link that hands over whatever the case selector
currently shows. The route serves plain files from the case directory only:
a name with a path separator, a name that escapes the directory, and a
CGMES delivery directory are refused.

## Configuration precedence

The case configuration file is OPTIONAL, and its PRESENCE is the switch,
for every input format alike: a case that ships a `<stem>.config.yaml`
is self-contained, so its case-scope keys resolve in two stages: what
the case levels state, then the packaged defaults. The machine's own
YAML configuration does not participate in those keys, and the same
case-plus-config pair computes the same numbers on every installation.
A case without a config file, a bare `.scf.json` included, keeps the
full chain with the machine's configuration between case levels and
defaults (the tracked fixture `warmup_casePST.scf.json` runs that way),
and so does the deprecated in-file block alone. Machine-scope keys
always come from the machine's configuration file.
`effective_config.yaml` of each run records what actually took effect.

A case file ships with the settings it was meant to run with. For a run
on a case that carries its configuration file, the precedence of a
case-scope key is, highest first:

1. explicit `config_overrides` from the API, CLI, or Web UI form,
2. the case configuration file `<stem>.config.yaml` next to the case,
3. the case file's own `sparlectra.config` (DEPRECATED: the writer no
   longer emits this block; the reader still applies it, with a warning
   naming the case configuration file as the new place),
4. the packaged template defaults.

The machine's YAML configuration file takes part only for machine-scope
keys of such a run; the full chain of a config-less run is documented in
[Configuration](configuration.md).

A case file decides how ITS network is computed, not how the installation
behaves. The permitted keys are therefore the GUI-editable allowlist minus
everything that describes the machine and the session: `output.*` (logging
and result tables), `benchmark.*`, `runtime.*` (parallelism), `webui.*` and
`matpower_export.*` stay in the configuration file, so a delivered case
cannot reconfigure someone else's reporting. What remains is the case
itself: import conventions (`matpower_import.*`, `cgmes_import.*`), the full
solver surface (`power_flow.*`, including accuracy, iteration limits, solver
choice, Q-limit strategy and distributed slack), `state_estimation.*` and
`short_circuit.*`. `scf_is_case_config_key(key)` answers it for a single
key. Every key that took effect appears in `effective_config.yaml` of the
run.

A key outside that scope is refused by name when the file is read, and the
writer does not produce one: a file that states a setting which quietly does
not apply is the silent acceptance this format exists to prevent. Case files
written before this rule carried the whole Web UI form and need one
re-export.

**In the Web UI the form is seeded from the case file.** The form posts a
value for every option it shows, and those are explicit overrides, level 1
above. Selecting a case file therefore moves its controls to the values the
file carries, the same way a saved settings profile does, and the page says
which settings came from the case. Without that seeding the form would post
its own defaults and silently outrank the file it had just loaded. A saved
settings profile for the case still wins over the file, and any control the
user edits wins over both.

## The reference: single, distributed, or a source

`roles.slack` states which reference the case was built with, and the reader
restores it:

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

The factors are the network's own participation factors, the same ones the
MATPOWER (`APF`) and CGMES (`normalPF`) importers produce. They travel with
the model, but they do not switch the feature on: the case still has to ask
for it in its `config` block, with

```json
"config": {
  "power_flow.distributed_slack.enabled": true,
  "power_flow.distributed_slack.p_mode": "imported"
}
```

`imported` is the mode that reads exactly those factors. `explicit` weights
(`power_flow.distributed_slack.weights`) are not part of the case-file
allowlist and stay in the YAML configuration, because a weight list keyed by
bus name is a run setting, not a property of the network. See
[Power-Flow Configuration](powerflow_configuration.md).

## One importer per format, conversion on request

Every input format has its OWN importer that builds the network
directly and hands it to the solver: MATPOWER through its native
constructor, DTF and CGMES through theirs, and an SCF file simply
through `build_net`. SCF is the preferred WORKING format, never a
mandatory way station: no import converts into it behind your back.
Converting a case INTO an SCF file is an explicit action with its own
result (the SCF export, the Case page's export button, the
shipped-case builders), performed by the format's converter. This
keeps imported values bit-exact: nothing is re-encoded on the way in,
so nothing is lost.

## The typed case in memory

`SCFCase` is the in-memory form of one SCF document: the PGM `data`
section as typed component vectors, the namespaced `sparlectra` block as
its document sub-objects. `read_scf_json(path)` parses and validates a
file into it, `write_scf_json(case, path)` serializes it back to the
canonical bytes, `build_net(case)` constructs the network and applies the
run configuration's net parameters exactly once, and `net_to_scfcase(net)`
is the export direction (same keywords as `exportSCF`). `importSCF` and
`exportSCF` are wrappers over these four; `scf_to_net`/`net_to_scf` keep
the historical dict entry points.

## What the round trip preserves

The contract is between a case file and the network built from it. For a
network written by Sparlectra, read back, and written again, the two files
are **byte identical**, and a run on the re-read network reproduces the
original numerically, not approximately: identical power-flow iterations
and identical bus voltages, identical state-estimation objective, degrees
of freedom, and state. The test suite checks this on the tracked fixture
and on case14, case57, case118 and case300.

The contract does not extend to other formats. A network imported from
MATPOWER, CGMES or DTF and exported as a case file keeps its electrical
model and its solution, but not every detail of the source file (CGMES
mRIDs are recorded where available, MATPOWER comments and column layout
are not), so reading the case file and writing the source format again is
a new export, not a reproduction of the original file.

Three details make the case-file contract possible, and each is a
deliberate choice:

- **Internal order is part of the identity.** The branch vector order, the
  prosumer order, and the measurement order all influence generated
  component ids or floating-point summation order, so `extra` records each
  element's internal index and `measurements.rows` records each row's
  position; the reader restores them.
- **The PGM sensor stays the authority for measured values**, but PGM
  carries one sigma per sensor in SI, which cannot express two different row
  sigmas (P and Q) exactly, and the SI conversion is not always
  bit-reversible. `measurements.rows` therefore carries the internal sigmas
  always, and a value only where reading it back would not reproduce the
  exact number.
- **Both sides compute the impedance base from the file's own numbers**
  (`u_rated` in V and `s_base` in VA) with the same formula, and the writer
  emits the SI value that the reader turns back into exactly the per-unit
  value it started from (it checks the neighbouring floats). Without that
  last step a file drifted by one ulp on EVERY cycle instead of being a
  fixed point.

A write, read, write cycle produces byte-identical files. That is a property
of the format, not of the numbers: unit conversion is not exactly invertible
in floating point, so the writer emits the SI value that the reader turns
back into exactly the per-unit value it started from (checking the
neighbouring floats), and falls back to a value that maps to itself. Without
this a file drifted by one ulp on every cycle. Verified on
`warmup_casePST`, case14, case57, case118 and case300.

What comes back bit-for-bit: bus and branch ORDER (they fix the Y-bus
numbering and the generated component names), the source system's bus
numbers, tap ratios, shunt susceptances including their sign, controllers,
measurements, and the network's power-flow result.

The operating point is part of that, and it does not live in the Y-bus: the
reactive band of every machine, its voltage setpoint, and whether it
regulates at all travel in `extra` (`max_q_mvar`, `min_q_mvar`, `vm_pu`,
`regulated`), because a machine written as a PGM `source` has no
`voltage_regulator` row to keep them in, and an unlimited band has no JSON
representation in the dataset at all. A setpoint alone must not promote a
machine either: the constructor treats any setpoint as regulation, so the
reader restores the flag the file states. Read the case file back, solve it,
and the voltages match the source network to the last bits.

Two things are restored within a documented tolerance rather than bit-exactly:

- **Per-unit values** can land one ulp away when no exact preimage exists
  (`0.05403` with a 0.01 impedance base has none). The resulting Y-bus
  deviation is at the 1e-16 level, far below any modelling accuracy.
- **Tap regulation bands** snap to the step grid, because a step-based
  format can only express band edges that lie on it. The operating tap
  ratio is exact; a band edge can move by less than one step. A source
  case whose neutral ratio lies outside its own band (MATPOWER case57)
  keeps that oddity and is named in a warning on load.

One case is known to break the byte-identical repeat export:
`case13659pegase`, where the band of a PHASE tap changer moves by one step
on the second write. The model is unaffected (identical Y-bus, identical
ratios, identical number of regulating machines), and the third write equals
the second. The ratio band and every other quantity are fixed points; this
one is not yet.

Values JSON cannot carry are encoded, never silently dropped: an unlimited
branch rating (MATPOWER `rateA = 0`, which arrives as `Inf`) is absent from
the PGM dataset, because "no `i_n`" IS "no limit" there, and travels as
`"inf"` in the namespaced block. A non-finite number anywhere else is a hard
error naming the key, not a `null` the reader would later reject.

**An unlimited reactive band is stated by absence.** `max_q_mvar` and
`min_q_mvar` are written only when the limit is finite; a missing limit means
there is none. The rule is the same in the dataset and in the namespaced
block, so there is no sentinel to interpret and no second spelling for the
same fact. It is verified that a missing and an infinite reactive limit
produce the same `qmin_pu`/`qmax_pu` and the same power flow.

### Fields that are not written

A field whose value equals its documented default is left out, and the reader
puts the default back. This is limited to the namespaced block: the `data`
section keeps every attribute power-grid-model requires, or it would stop
being a valid PGM dataset.

| Field | Default | Meaning when absent |
|---|---|---|
| `extra.<branch>.kind` | `"BranchC"` | an ordinary line |
| `extra.<machine>.vm_pu` | absent | the machine does not regulate (a regulating one keeps its setpoint in its `voltage_regulator` row, the reference in `source.u_ref`) |
| `extra.<machine>.max_q_mvar`, `min_q_mvar` | absent | no reactive limit |
| `extra.<machine>.max_p_mw`, `min_p_mw` | absent | no active limit |
| `extra.<machine>.regulated`, `apu_node` | `false` | |
| `extra.<branch>.meta.sn_mva` | absent | no rating recorded |
| `components.tap_changer[].tap_est_mode` | `"none"` | the tap is not released for estimation |
| `components.transformer3w[].status` | `1` | in service |
| `roles.slack.normalize` | `true` | participation factors are normalized |

`area` and `zone` are NOT in this list: their constructor value is "unknown",
not `1`, so writing them only when set is what keeps them faithful.

`data/scf/sp_casePST.scf.json` is the permanent version 1 fixture: a
tracked file that every future reader must keep loading unchanged.

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
siemens, farad, ampere, radian). `sparlectra` holds everything PGM has no
model for. A reader that knows nothing about Sparlectra still gets a
complete, computable network; what it loses is listed under
[Interoperability](#Interoperability).

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

The reference side of `generic_branch` is not a convention choice: the
Sparlectra Y-bus stamps the tap on the **from** side, so its `r_pu`/`x_pu`
are referenced to the **to** side, exactly like PGM. Impedances are always
written from the **base** (physical) values; exporting a network that
carries a FACTS-compensated operating point is a hard error naming
`restoreBaseImpedances!`, because such a file would silently misrepresent
the equipment.

### The `sparlectra` block

| Key | Content |
|---|---|
| `format_version` | schema version of the block (currently `1.0`) |
| `meta` | case name, `s_base` in VA, `f_nom` in Hz, source format and reference, producer, and `intended_calculations` as the file's contract |
| `roles` | what PGM cannot express: `slack.mode` (`single`, `distributed`, `source`) with its nodes or participation factors, auxiliary nodes, isolated nodes |
| `extra` | per-component annotations keyed by id: `name` (the reference name every other section resolves against), `external_id` (CGMES mRID, MATPOWER index), bus/branch index, node type, ratings |
| `components.tap_changer` | the tap cascade PGM has no model for: per controller the step grid, position, band, and nameplate angle, plus the neutral ratio and angle |
| `components.transformer3w` | three-winding transformer: either the GROUPING of three existing `generic_branch` legs (their star node, and which end is hv/mv/lv; the electrical values stay in the legs and are never duplicated) or a NAMEPLATE that describes the transformer directly, see below |
| `components.sc_source` | IEC 60909 source data (external network injections, machines) that PGM's `source` cannot carry beyond `sk`/`rx_ratio` |
| `components.shunt_state` | what PGM's `g1`/`b1` cannot say: whether a shunt is in service, whether it is a voltage-dependent injection rather than an admittance, and whether its susceptance is released as a state-estimation state. Only deviating shunts produce a row |
| `transformer_types` | named transformer nameplates in CGMES `PowerTransformerEnd` vocabulary, referenced by `components.transformer3w` entries |
| `components.controllers` | FACTS and regulation, in the declarative `control.controllers` schema verbatim: same type names, same keyword names. The reader hands the entries to `applyConfiguredControllers!`, so there is one construction path and no second vocabulary |
| `contingencies` | LEGACY N-1 study definition: `mode` (`explicit`, `all_branches`, `all_branches_plus`), the case list with their outages, exclusions. Still read, mapped onto the scenario model at load time; the writer emits `scenarios` only |
| `scenarios` | the scenario model (0.10.0): `mode` (`explicit`, `n1_branches`, `n1_generators`, `n1_all`), `exclusions`, and the ordered scenario list, each scenario `name`, `weight` and its `ops` (`op` = `status`/`set`/`scale`, `target` component class, `id` the SCF component id, plus `value`, `field`, `factor` as the op needs). The N-1 modes expand through the same generators `runContingencies!` uses, so N-1 is the special case of the model; see [N-1 Contingency Analysis](contingency.md) and the Web UI's scenario editor |
| `short_circuit` | IEC 60909 study definition: `case` (`max`/`min`), optional `c_factor`, `sweep` (`all_buses`/`explicit`) and, for an explicit sweep, the `buses` list of node ids |
| `measurements.rows` | what PGM sensors cannot carry: the Sparlectra row ids and sigmas behind each sensor, their positions in the measurement vector, active flags where they deviate, and a value only where the SI conversion is not bit-reversible |
| `config` | case-specific dotted configuration keys (same allowlist the API and Web UI use) |
| `start_state` | opt-in bus voltages as START values, never as results |

The study blocks state WHAT to compute; the result contract of the runs
themselves is unchanged. They are validated on load (an unknown mode, an
out-of-band `c_factor`, or an empty outage list fails when the file is read,
not at the start of a long sweep) and are readable without building the
network through `scf_case_studies(file)`.

### Three-winding transformers without writing legs

Writing three `generic_branch` legs plus a star node by hand is the part of
a case file people get wrong. A `transformer3w` entry may therefore describe
the transformer instead of pointing at legs, in the vocabulary CGMES uses
for it: one `PowerTransformerEnd` per winding with `rated_u`, `rated_s`,
`r`, `x` and `b`. The reader builds the star equivalent through the same
constructor a hand-written network uses (`add3WTPiModelTrafo!`), so there is
no second modelling path, and it creates the auxiliary star node itself.

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

A `type` entry is the same nameplate under a name, so several transformers
of one design are written once and referenced; an entry can also carry its
`nameplate` inline. The two forms are exclusive: an entry either points at
legs or describes them, and mixing them fails on load.

#### Which CIM attributes the nameplate speaks

One nameplate end is one CIM `PowerTransformerEnd` of a `PowerTransformer`.
These attributes are read:

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

All three ends state `r`, `x` and `b` **on the HV base**, which is the
star-equivalent form the CGMES import produces; the reader converts to per
unit with the highest `rated_u` of the three as the voltage base.

What the nameplate deliberately does NOT take from CIM: the magnetising
conductance `g` (a three-winding equivalent in this format carries its
losses in `r`), `phaseAngleClock` and `connectionKind` (the vector group,
which the balanced positive-sequence model does not evaluate), and the tap
changers. A `RatioTapChanger` or `PhaseTapChanger` on the transformer lives
in its own `components.tap_changer` entry, and the group points at it with
`tap_changer`; controller `index` 1 is the ratio changer, `index` 2 the
phase changer, each with its step size, live position and position band.
Anything else in the nameplate is an unknown key and fails on load rather
than being ignored.

This is an INPUT shorthand. An export always writes the explicit legs plus
the grouping, because that is the form that round-trips; a file that used
the nameplate therefore comes back in the explicit form, describing the same
network.

#### What the star equivalent is called

The reader creates the star point and the three legs itself, so their names
are generated, and everything that addresses a component by name later (a
contingency case, a measurement, a tap block) has to use exactly these:

| Object | Name | Built from |
|---|---|---|
| star node | `Aux3WT_<hv>_<mv>_<lv>` | the three terminal bus names |
| leg | `B_2WT_<Vn>_<star>_<terminal>` | rated voltage in kV, then the two bus indices |

For the nameplate above, with terminals named `HV`, `MV` and `LV`, the
export writes the star node as `Aux3WT_HV_MV_LV` and the three legs as
`B_2WT_220_4_1`, `B_2WT_220_4_2` and `B_2WT_220_4_3`. The star node is
listed under `roles.aux_nodes`, so a reader can tell it apart from a real
busbar, and the grouping states which leg is which end:

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

`leg_direction` says it outright: every leg runs from the star node to its
terminal, so no reader has to infer an orientation. The roles follow the
terminal voltages, highest first. A CGMES delivery keeps its own names
instead, because the importer takes them from the `PowerTransformer` and its
ends; the generated names above appear where a case file builds the
transformer from a nameplate.

### Running the studies

A study run on a case file executes what the file defines, so shipping a
case ships the study it was meant for.

**N-1 and scenarios.** The `scenarios` block is the study definition since
0.10.0: an ordered list of named, weighted patch scenarios (a scenario may
combine several status, setpoint and scaling operations, so a double
outage IS one scenario), plus the N-1 modes that expand at load time. The
service runs it with `scenario_source = file_block`, the Web UI's scenario
editor writes it, and `runScenarios!` is the programmatic entry. The
legacy `contingencies` block stays readable and is mapped onto the
scenario model at load time (`mode` and exclusions carry over, every
listed case becomes one scenario with a status patch per outage
component); the writer emits `scenarios` only, and no existing file needs
regeneration. On the legacy service path (no scenario source in the
request) the historical semantics hold: `explicit` runs exactly the
listed cases, `all_branches` generates the sweep and applies `exclude`,
`all_branches_plus` adds the listed extras, outages resolve through
`extra[<id>].name`, and a legacy case with more than one outage is
refused there (the scenario path is the one that runs multi-op cases).
The per-case `weight` is honored, and `run.log` records where the case
list came from.

**State estimation.** A case file carries its measurements with the model,
so an estimation runs on it without a separate measurement CSV: the set that
came with the case is used, and the run still writes `measurements.csv` as
its artifact, so it stays reproducible from its own output. Picking a CSV
explicitly still works and then wins. In the Web UI the measurement selector
of a case file offers "(from the case file)" as its first entry, and the run
is offered even when the case cache holds no CSV at all. A case file without
measurements says exactly that instead of complaining about a missing file.

**Short circuit.** The `short_circuit` block selects the headline case
(`max` or `min`), an optional `c_factor`, and the sweep. `sweep: explicit`
takes a `buses` list of node ids, again resolved through `extra`; an unknown
id fails the run instead of producing an empty table. The same selection
also exists in PGM's own vocabulary: `data.fault` rows name the faulted
node, and a run without an explicit sweep uses them, so a case written in
PGM terms is runnable as it stands. Sparlectra evaluates the balanced bolted
fault, so a `fault_type` other than `three_phase` or a non-zero `r_f`/`x_f`
is refused on load rather than silently ignored. Both cases are always
evaluated, so `short_circuit_max.csv` and `short_circuit_min.csv` exist for
every format and the file's `case` only decides where the headline numbers
come from. The `c_factor` in the file applies when the configuration is at
its default and loses against an explicitly set `short_circuit.c_factor`, in
line with the precedence for the rest of the case configuration.

Short circuit on a case file needs the sources in
`components.sc_source`; a file without them says so instead of reporting a
zero-current network. That is also what the Web UI's "Short circuit" button
checks before it enables itself. Empty source fields are written as `null`
rather than dropped, so a record keeps its shape across a round trip.

Both study runs also apply the case file's own `config` block, the same way
a power-flow run does.

## Deterministic ids and readable diffs

PGM requires integer ids. The writer assigns them by a documented rule:
components are walked in a fixed type order and, inside a type, in
lexicographic name order. Exporting the same network twice therefore
produces byte-identical files, and rebuilding a case from its source gives
the same ids again, so Git diffs stay reviewable. Floats use the shortest
round-trip representation, which is lossless.

The file is pretty-printed throughout: every object breaks one key per
line, and only a list of plain scalars (node ids, for instance) stays
inline. A diff therefore points at the changed FIELD. Version 1 files up to
0.10.0 wrote each component as a single compact line, which kept diffs
narrow but read like a data stream; the content is unchanged, so an older
file and a re-exported one differ in whitespace only.

`extra[<id>].name` is the REFERENCE name: the name measurements,
controllers, contingencies, and reports resolve against (Sparlectra's bus
dictionary key). Where the internal component name differs, it travels as
`component_name`. The source-system id goes into `external_id`: for a CGMES
delivery that is the mRID (marked `source_id_kind: cgmes_mrid`), for other
formats the internal component id.

This makes OWN bus names a first-class channel, independent of the source
id: a file may say `"name": "Ostheim_110", "external_id": "H1"` and the
import uses the place name everywhere while the technical handle stays
addressable. The import restores `external_id` and `component_name` onto
the built network, so re-exporting an imported case preserves both fields
byte for byte. The [shipped demo cases](demo_cases.md) use exactly this
split: invented place names as reference names, short handles as
`external_id`.

## Source formats

The writer serializes a `Net`, so it is format-agnostic by construction: a
case imported from MATPOWER, CGMES, or DTF exports the same way, and the
format-specific payload each importer put into the network travels along.

| Source | What the file carries |
|---|---|
| MATPOWER | bus and branch names, the `mpc.sparlectra` extensions (tap-changer nameplates as `components.tap_changer`, busbar couplers as `data.link`), branch kinds, ratings |
| CGMES | mRIDs as `external_id`, the three-winding star grouping with its auxiliary node, IEC 60909 source data, the reference names of the topological nodes |
| DTF | station reference names, transformer tap data, ratings |

Three-winding transformers need no conversion: Sparlectra already builds the
star equivalent PGM documents for `generic_branch` (three legs plus an
auxiliary star node), so the writer serializes it directly and the grouping
block only records which legs belong together.

## Measurements

Sparlectra measurement rows map onto PGM sensors: a bus `Vm` (plus a `Va`
if present) becomes one `sym_voltage_sensor`, a `Pinj`/`Qinj` pair one
`sym_power_sensor` with `measured_terminal_type = node`, a `Pflow`/`Qflow`
pair one branch power sensor, and current magnitude/angle one
`sym_current_sensor`. Where a location carries only one half of a pair, the
missing attribute is `null` and the Sparlectra rows behind the sensor are
listed in `measurements.rows`, so nothing is lost.

A case that carries measurements needs no separate CSV to be estimated: a
state-estimation run on it uses the set that came with the model and still
writes `measurements.csv` as its own artifact, so every run stays
reproducible from its output alone.

`measurements.provenance` records what the numbers are worth. It keeps the
generator and seed, whether the values carry noise (an ideal set has no noise
at all, and its $J$ is zero by construction, which proves nothing about the
estimator), the per-row truth values, and the tap deviations the generator
applied. Those last two are what a run needs to produce `se_deltas.csv` and
to warn when a documented tap deviation cannot be absorbed because tap
estimation is off. Without them a case-file run would report a large $J$ and
leave you guessing why.

An ideal set can be perturbed in place, without solving a power flow again:

```julia
addMeasurementNoise!(net)   # uses each row's own sigma
```

## Interoperability

| Level | What you get |
|---|---|
| **PGM-readable** | the `data` section alone is a valid PGM input dataset |
| **Sparlectra-complete** | the full file: slack strategy, tap cascade, names, measurement detail, configuration |
| **Lossy toward PGM** | a PGM reader loses the tap cascade detail (the live ratio survives in `k`/`theta`), the slack strategy (it sees `source` components), the names, and the configuration |
| **Strict PGM** | `exportSCF(net; strict_pgm = true)`, the "Export as plain PGM" button, or `sparlectra export --pgm` writes the dataset ALONE, for a consumer that reads PGM and nothing else |

Strict mode drops the namespaced block, names what it dropped in a warning,
and rewrites a slack generator as a PGM `source`: PGM has no slack flag, its
reference IS a source, and leaving the generator in place would make the node
inject twice. Reading such a file back therefore loses names, roles, tap
nameplates, measurements and configuration, which is the point of the mode,
not a defect.

Reading in the other direction needs nothing special: a file written by
power-grid-model has no `sparlectra` block at all, and its in-service
`source` components become the network's reference. PGM's `source` IS
Sparlectra's external network injection, so it is read and written as one:
a case that arrives with a source leaves with a source, not with a generator.

A complete PGM source also states `sk` and `rx_ratio`, because PGM models it
as a voltage source BEHIND that impedance. Those numbers are read as feeder
data the case DECLARES: a short circuit can run on the file as it stands,
and `power_flow.external_grid` (with `source: auto`, see
[Power-Flow Configuration](powerflow_configuration.md)) computes with the
file's own values instead of the configuration defaults. The default stays
the ideal slack, which is Sparlectra's documented behavior; enabling the
external grid moves the reference behind $z = U_n^2 / S_k''$ and reproduces
what PGM computes. An export writes `sk`/`rx_ratio` back whenever the
network carries the feeder data.
`data/scf/pgm_interop.json` is that case as a tracked fixture, including a
`generic_branch` (the component PGM users have the fewest examples for); the
test suite checks the solved voltages against an independently derived model
rather than against Sparlectra's own numbers.

## Version and compatibility

`sparlectra.format_version` is written from day one, and the writer always
emits the current revision. A reader accepts exactly its own revision and
refuses anything else BY NAME, stating what it found, what it expected, and
that the case has to be re-exported. A file without the field is refused the
same way. A dataset written by power-grid-model has no namespaced block at
all and is never subject to this check.

**The comparison rule is exact equality, and that is a deliberate choice for
the pre-release phase.** While the format still moves, a reader that guessed
at an older file would be more dangerous than one that refuses it. It also
means the next revision breaks every file of this one. Before the format is
published, the rule has to be replaced by an explicit compatibility
statement: a MINOR increment stays readable (only additions and omissions
whose default the reader knows), a MAJOR increment does not. The revision
is a single string today, so that change is itself a breaking one, which is
why it belongs before publication rather than after.

Backward compatibility is not promised. The format is young and still
moving: while it is, a version bump may change or drop what an older file
carries, and an older file may need to be re-exported rather than migrated.
A breaking change is named in the changelog. Treat a case file as an
exchange and archive format for the version that wrote it, not as a
guaranteed long-term container.

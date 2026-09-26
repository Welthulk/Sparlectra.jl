# [PowSyBl Import](@id powsybl_import)

IIDM is the native network format of [PowSyBl](https://www.powsybl.org),
the open-source power system framework of the LF Energy foundation; its
XML variant XIIDM (`.xiidm`, `.xml`) carries node-breaker topology, tap
changers, HVDC links and limits. Sparlectra reads the file itself, in
Julia, without Python: `import_case`, `run_sparlectra` and the Web UI take
a `.xiidm` file as they take a MATPOWER case, and the import reproduces the
OpenLoadFlow solution of the same network.

```julia
using Sparlectra
res = run_sparlectra(casefile = "network.xiidm")
```

## What the reader resolves

The reader (`read_iidm_tables`) turns the file into the same tables that
pypowsybl's getters deliver, so one builder serves both sources:

- the bus-breaker and bus views of a node-breaker substation: nodes merge
  through internal connections and closed switches, a retained switch
  stays in the bus-breaker view (as a link of the built network), a bus
  is valid on powsybl's rule (a busbar section with a feeder, or a branch
  with two feeders), and the buses are named `<voltage level>_<lowest
  node>` as powsybl names them;
- the connected and synchronous components, numbered by size;
- the current step of every ratio and phase tap changer: `rho`, `alpha`
  and the step corrections of `r`, `x`, `g`, `b`, on two- and
  three-winding transformers;
- the reactive limits of a generator or VSC station at its active
  setpoint, interpolated on a capability curve as powsybl interpolates;
- the operational limits per side and limits group, permanent and
  temporary;
- dangling and tie lines, HVDC lines with their converter stations,
  static var compensators, linear and non-linear shunts.

Extensions, areas and properties are ignored. Two constructs are refused
with a message naming them: a compressed file (`.xiidm.bz2`, `.gz`: unpack
it first) and the tie-line form of IIDM versions before 1.10 (the two
half lines inline). The bus state written in the file seeds the start
voltages, unless a tap changer or shunt was solved at another position
than the one the file carries (`solvedTapPosition`, `solvedSectionCount`);
the import then starts flat and says so in the report.

### IIDM versions

The reader takes the IIDM XML schema versions 1.0 to 1.17 (the namespace
`http://www.powsybl.org/schema/iidm/1_<n>`). It is verified against files
of version 1.17, which is what pypowsybl 1.16 writes; for older files it
knows the element forms that changed on the way: the linear shunt model
on the element itself (before 1.3), `danglingLine` with `ucteXnodeCode`
(before 1.17, `boundaryLine` with `pairingKey` since), tie lines by
`danglingLineId1`/`2` (1.10 to 1.16) and `boundaryLineId1`/`2` (1.17),
`currentLimits1`/`2` on the element (before 1.12) next to the limits
groups, a static var compensator without the `regulating` attribute
(mode `OFF` then means not regulating), `targetV` on a ratio tap changer
(before 1.12) next to `regulationValue`. A construct of a newer schema
that the reader does not know is ignored, never guessed; a construct it
cannot map is refused with its name. Compressed files are not read.

## Two sources, one builder

| Source | What runs | Needs |
|---|---|---|
| An **IIDM file** (`.xiidm`, `.xml`) | `read_iidm_tables`, then `build_net_from_powsybl` | nothing but Sparlectra |
| A **table bundle**, a directory `<case>.powsybl` with `manifest.json` and one CSV per pypowsybl table | `read_powsybl_bundle`, then the same builder | nothing but Sparlectra |

`import_case` and the Web UI detect both (a file whose first bytes carry
the IIDM namespace, a directory whose name ends in `.powsybl` with the
bundle manifest). The bundle is the reference form: pypowsybl writes it
with OpenLoadFlow's solution in the state columns and `reference_buses.csv`
next to it, which is what the test suite compares against. The five test
fixtures ship as bundles with their `.xiidm` files, and the reader is
checked against them column by column.

## Reference bundles (optional)

Nothing here is needed to import a network. A table bundle is reference
material: pypowsybl writes every table with all attributes, runs
OpenLoadFlow (default parameters, except a slack mismatch bound of 1e-4
MW and a Newton-Raphson tolerance of 1e-9 per equation so that a 1e-6 pu
comparison has a converged reference) and writes `reference_buses.csv`
and the manifest next to the CSV files. That is how the five shipped
bundles were produced (pypowsybl 1.16.1), and how a comparison of
Sparlectra against OpenLoadFlow on a network of your own is set up: the
tables of the bundle carry OpenLoadFlow's solution in their state
columns, the reader's tables carry the file's state, everything else is
equal. CSV rules of a bundle: comma separated, strings always quoted,
floats in full precision with `NaN` and `Inf` as literals, booleans
`true`/`false`, a missing integer as `-1`; `read_powsybl_bundle` and
`write_powsybl_bundle` read and write the form.

## Bus model

Every bus of PowSyBl's bus-breaker view becomes one Sparlectra bus, named
by its bus-breaker id, with the bus-view id in the external-id channel so
results join back to the OpenLoadFlow bus table. Retained switches become
links (closed or open), non-retained switches never appear (PowSyBl
already contracted them), and the nominal voltage comes from the voltage
level. Three-winding transformers get a star bus `<id>_star` at
`rated_u0`; a dangling line runs to a boundary bus `<id>_boundary` that
carries its `p0`/`q0` as a load, and the two dangling lines of a tie line
meet at `<pairing_key>_xnode`.

## Example files

Three of the bundles also ship under `data/powsybl` (`ieee14`,
`four_substations`, `micro_grid_be`), each with its `.xiidm` file, and
appear in the Web UI case selector.

| case | source | what it covers |
|---|---|---|
| `ieee14` | `pypowsybl.network.create_ieee14` | the IEEE 14-bus case: Y-bus identity with MATPOWER `case14` on every entry, including the two branches PowSyBl imports as lines between voltage levels (their one-sided shunts become bus shunts, see Conventions) |
| `ieee57` | `create_ieee57` | 57 buses, 17 transformers, two parallel transformer pairs, OpenLoadFlow's distributed slack |
| `four_substations` | `create_four_substations_node_breaker_network` | node-breaker topology with retained switches, a phase-shifting transformer, VSC and LCC HVDC links, an SVC, curve reactive limits, temporary limits, two synchronous components |
| `micro_grid_be` | `create_micro_grid_be_network` | CGMES origin, a three-winding transformer with a ratio tap changer, remote voltage regulation, dangling lines |
| `eurostag_tie_lines` | `create_eurostag_tutorial_example1_with_tie_lines_and_areas` | tie lines built from paired dangling lines |

Each fixture directory under `test/fixtures/powsybl/` holds the bundle,
the `.xiidm` file pypowsybl wrote and `reference_buses.csv`;
the bundles were written with pypowsybl 1.16.1. The test suite reads
every `.xiidm` with the native reader and requires its tables to equal
the bundle's in every column except the state columns, and the power flow
from the `.xiidm` to land in the same bands as from the bundle.

## Conventions

Settled against the OpenLoadFlow voltages of the example files:

- **Transformer ratio.** Sparlectra's ratio is the MATPOWER tap, the no-load
  ratio of the from-side to the to-side voltage in pu. IIDM's `rho` is the
  ideal ratio at side 1 in kV and already carries `rated_u2 / rated_u1`, so
  `ratio = vn_to / (rho * vn_from)`; the phase shift is `-alpha`. The
  `_at_current_tap` impedances are used, no tap tables reach the network.
- **Magnetizing admittance.** The whole `g`, `b` of a transformer sits on
  side 1, behind the ideal transformer, as OpenLoadFlow places it; three-
  winding legs alike. The branch carries it as its from-terminal arm
  (`g_from_pu`, `b_from_pu`), a line keeps `g1`, `b1` and `g2`, `b2` on
  their own terminals, and a dangling line its admittance on the network
  terminal; no bus shunt is created for any of them, so the shunt list of
  an imported network holds the shunt compensators only, an outage takes
  the admittance away with its branch, and the SCF export keeps the split
  (see [Branch model](branchmodel.md)).
- **Lines between voltage levels.** PowSyBl keeps a plain conductor and
  OpenLoadFlow's default line model is exactly that; Sparlectra builds it
  as the ratio branch `vn_to / vn_from` with the impedance on the to-side
  base, the same form the CGMES import uses.
- **Flow sign.** A PowSyBl branch column `p1` is positive when power flows
  from the side-1 bus into the branch; generator `target_p` is positive
  for generation, the result columns `p`, `q` of injections are in load
  convention.
- **HVDC.** Fixed injections by default: the rectifier station draws
  `target_p / (1 - loss_factor / 100)`, the inverter delivers
  `target_p * (1 - loss_factor / 100)`, as OpenLoadFlow reports them; a
  VSC that regulates voltage is PV at its bus, an LCC absorbs
  `P * tan(acos(power_factor))`.
- **Active power balance.** OpenLoadFlow distributes the mismatch over the
  generators with a nonzero target proportionally to `max_p`; the importer
  carries that rule as participation factors, so
  `power_flow.distributed_slack.enabled = true` with `p_mode = imported`
  reproduces OpenLoadFlow's balance; in that mode a unit takes part
  wherever it sits, also as PQ (a non-regulating generator, or a unit the
  `remote` regulation mode runs as PQ). Without it the slack generator of
  each component absorbs the whole mismatch.
- **Slack choice.** Per synchronous component the regulating unit with
  the largest `max_p`, a unit whose setpoint applies to its own bus before
  one that regulates a remote bus (the slack bus keeps its start voltage,
  so a remotely regulating slack would leave its own bus wherever the
  file's state put it); `slack_ids` overrides the choice.
- **Reactive limits.** OpenLoadFlow switches a unit at its limit without
  hysteresis; Sparlectra's `power_flow.qlimits.hysteresis_pu` (default
  0.01) keeps a unit PV inside the band. Set it small to match
  OpenLoadFlow bus by bus.
- **Remote voltage regulation.** `powsybl_import.remote_regulation =
  remote` attaches an outer-loop machine voltage control on the regulated
  bus (the machinery of `cgmes_import.machine_control`), which reproduces
  OpenLoadFlow's remote control; `hold_local` (the default) keeps the unit
  PV at its own bus.

## [Configuration (`powsybl_import`)](@id powsybl-import-config)

The `powsybl_import` scope (canonical section `powsybl`, alias
`powsybl_import`), see [Configuration](configuration.md):

| Key | Default | Meaning |
|---|---|---|
| `powsybl_import.base_mva` | `100.0` | System base in MVA. |
| `powsybl_import.hvdc_mode` | `fixed_injection` | `fixed_injection` or `paired_control`; the latter is not implemented for PowSyBl sources in this release and is rejected with the mode named. |
| `powsybl_import.slack_ids` | `[]` | Generator ids that override the slack choice of their synchronous component; a YAML list or one string with `;` between the ids. |
| `powsybl_import.multi_slack` | `true` | One slack per synchronous component; `false` keeps only the component of the first slack and reports the others. |
| `powsybl_import.remote_regulation` | `hold_local` | `hold_local`, `pq` (PQ with `target_q`) or `remote` (outer-loop machine voltage control on the regulated bus). |

The import report (counts per element type, every skipped element with
its reason, the slack decision per component, notices) is printed by
`format_powsybl_report` and travels with the imported case as
`provenance["powsybl_report"]`.

## Out of scope

Snapshots with missing injections and injection patching, time-series
batches, tap and voltage controllers from IIDM regulation data (the step
tables are read, not used), IIDM export, compressed IIDM files, a Web UI
directory picker for bundles (a bundle is selected by its path).

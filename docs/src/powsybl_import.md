# [PowSyBl Import](@id powsybl_import)

IIDM is the native network format of [PowSyBl](https://www.powsybl.org);
its XML variant XIIDM (`.xiidm`, `.xml`) carries node-breaker topology,
tap changers, HVDC links and limits. Sparlectra reads the file in Julia,
without Python: `import_case`, `run_sparlectra` and the Web UI take a
`.xiidm` file as they take a MATPOWER case, and the import reproduces the
OpenLoadFlow solution of the same network.

```julia
using Sparlectra
res = run_sparlectra(casefile = "network.xiidm")
```

## What the reader resolves

The reader (`read_iidm_tables`) turns the file into the tables PowSyBl's
own getters deliver, and the builder (`build_net_from_powsybl`) makes the
network from them. It resolves:

- the bus-breaker and bus views of a node-breaker substation: nodes merge
  through internal connections and closed switches, a retained switch
  stays in the bus-breaker view (a link of the built network), a bus is
  valid on powsybl's rule (a busbar section with a feeder, or a branch
  with two feeders), and buses are named `<voltage level>_<lowest node>`
  as powsybl names them;
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

Areas and properties are ignored; of the extensions only the
short-circuit data of a generator (`generatorShortCircuit`, see
[Run kinds](@ref powsybl_run_kinds)) is read. Two constructs are refused
with a message naming them: a compressed file (`.xiidm.bz2`, `.gz`: unpack
it first) and the tie-line form of IIDM versions before 1.10 (the two
half lines inline). The bus state written in the file seeds the start
voltages, unless a tap changer or shunt was solved at another position
than the one the file carries (`solvedTapPosition`, `solvedSectionCount`);
the import then starts flat and says so in the report.

### IIDM versions

The reader takes the IIDM XML schema versions 1.0 to 1.17 (the namespace
`http://www.powsybl.org/schema/iidm/1_<n>`). It is verified against files
of version 1.17 (what pypowsybl 1.16 writes); for older files it knows
the element forms that changed on the way: the linear shunt model on the
element itself (before 1.3), `danglingLine` with `ucteXnodeCode` (before
1.17, `boundaryLine` with `pairingKey` since), tie lines by
`danglingLineId1`/`2` (1.10 to 1.16) and `boundaryLineId1`/`2` (1.17),
`currentLimits1`/`2` on the element (before 1.12) next to the limits
groups, a static var compensator without the `regulating` attribute
(mode `OFF` then means not regulating), `targetV` on a ratio tap changer
(before 1.12) next to `regulationValue`. A construct of a newer schema
that the reader does not know is ignored, never guessed; a construct it
cannot map is refused with its name.

## Source

`import_case` and the Web UI detect an IIDM file by its first bytes (the
IIDM namespace), not by the extension.

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

`data/powsybl` ships three of PowSyBl's own example networks and a
short-circuit variant of the IEEE case; the Web UI case selector offers
them, and the import reproduces the OpenLoadFlow voltages of each.

| case | what it covers |
|---|---|
| `ieee14` | the IEEE 14-bus case: Y-bus identity with MATPOWER `case14` on every entry, including the two branches PowSyBl imports as lines between voltage levels (each keeps its one-sided shunt on its own terminal, see Conventions) |
| `four_substations` | node-breaker topology with retained switches, a phase-shifting transformer, VSC and LCC HVDC links, an SVC, curve reactive limits, temporary limits, two synchronous components |
| `micro_grid_be` | CGMES origin, a three-winding transformer with a ratio tap changer, remote voltage regulation, dangling lines |
| `ieee14_sc` | the IEEE case with short-circuit data on its generators, see [Run kinds](@ref powsybl_run_kinds) |

## Conventions

Settled against the OpenLoadFlow voltages of the example files:

- **Transformer ratio.** Sparlectra's ratio is the MATPOWER tap, the no-load
  ratio of the from-side to the to-side voltage in pu. IIDM's `rho` is the
  ideal ratio at side 1 in kV and already carries `rated_u2 / rated_u1`, so
  `ratio = vn_to / (rho * vn_from)`; the phase shift is `-alpha`. The
  `_at_current_tap` impedances are used, no tap tables reach the network.
- **Magnetizing admittance.** The whole `g`, `b` of a transformer sits on
  side 1, behind the ideal transformer, as OpenLoadFlow places it; the
  branch carries it as its from-terminal arm, a line keeps `g1`, `b1` and
  `g2`, `b2` on their own terminals, a dangling line its admittance on the
  network terminal; no bus shunt is created for any of them
  ([Where the magnetizing admittance sits](@ref magnetizing_placement)).
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
  reproduces OpenLoadFlow's balance (a unit takes part wherever it sits,
  also as PQ). Without it the slack generator of each component absorbs
  the whole mismatch.
- **Slack choice.** Per synchronous component a regulating unit, by the
  ranking every reference choice uses ([`reference_candidate_rank`](@ref)):
  a unit whose setpoint applies to its own bus before one that regulates a
  remote bus (a remotely regulating slack would leave its own bus at the
  file's start voltage), then the largest unit, by rated power where the
  file states one and by `max_p` otherwise; `slack_ids` overrides the
  choice.
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
| `powsybl_import.hvdc_mode` | `fixed_injection` | `fixed_injection` or `paired_control`; the latter is not implemented for PowSyBl sources and is rejected with the mode named. |
| `powsybl_import.slack_ids` | `[]` | Generator ids that override the slack choice of their synchronous component; a YAML list or one string with `;` between the ids. |
| `powsybl_import.multi_slack` | `true` | One slack per synchronous component; `false` keeps only the component of the first slack and reports the others. |
| `powsybl_import.remote_regulation` | `hold_local` | `hold_local`, `pq` (PQ with `target_q`) or `remote` (outer-loop machine voltage control on the regulated bus). |

The import report (counts per element type, every skipped element with
its reason, the slack decision per component, notices) is printed by
`format_powsybl_report` and travels with the imported case as
`provenance["powsybl_report"]`; a run of the service or the Web UI writes
it as the artifact `powsybl_import.log`.

## [Run kinds](@id powsybl_run_kinds)

A PowSyBl case runs through every run kind of the service and the Web UI.

| Run kind | What the case contributes |
|---|---|
| Power flow, Diagnose | the network as imported |
| State estimation | the network; measurements come from a generated or uploaded set |
| N-1 | branch and generator outages generated from the network, solved with the slack model of the run configuration; scenario files address component ids of a case file, so they need an export to the Sparlectra Case Format first |
| Short circuit | the connected generators as sources |
| Exports | the Sparlectra Case Format file (Case page), the MATPOWER case and the CGMES delivery as run artifacts |

For the short circuit the reader takes `directSubtransX`, `directTransX`
and `stepUpTransformerX` of the extension `generatorShortCircuit`
(reactances in ohm at the generator's voltage level; a step-up reactance
lies in series with the subtransient one). A generator without the
extension is evaluated with the default `x''_d` of
[Short-Circuit Analysis](short_circuit.md), and every row that depends on
it carries the flag and the reason; the run status follows the rule of
every format ([Source data and run status](@ref short_circuit_source_data)).
Of the shipped files `ieee14.xiidm` carries neither reactances nor rated
powers (the run is refused), `micro_grid_be.xiidm` carries rated powers
only (defaults with a warning), and `ieee14_sc.xiidm` carries complete
data: the same network with assumed ratings and `x''_d = 0.2` pu, built
by `tools/build_powsybl_sc_demo.jl`.

## Out of scope

Snapshots with missing injections and injection patching, time-series
batches, tap and voltage controllers from IIDM regulation data (the step
tables are read, not used), IIDM export.

# [PowSyBl Import](@id powsybl_import)

IIDM is the native network format of [PowSyBl](https://www.powsybl.org),
the open-source power system framework of the LF Energy foundation; its
XML variant XIIDM (`.xiidm`, `.xiidm.bz2`, `.xml`) carries node-breaker
topology, tap changers, HVDC links and limits. Sparlectra does not parse
IIDM: pypowsybl loads the file, resolves the topology and the tap
positions, and hands over plain tables; Sparlectra builds its network
from those tables and reproduces the OpenLoadFlow solution of the same
network.

## Two entry points

| Source | What runs | Needs |
|---|---|---|
| A **table bundle**, a directory `<case>.powsybl` with `manifest.json` and one CSV per pypowsybl table | `read_powsybl_bundle`, then `build_net_from_powsybl` | nothing but Sparlectra |
| An **IIDM file** read live | the PythonCall extension's `read_powsybl_network`, then the same builder | PythonCall and pypowsybl in the session |

The bundle is the primary path: fixtures, examples, tests and the
`run_sparlectra` case picker work on bundles. `import_case` and the Web UI
detect both (a directory whose name ends in `.powsybl` with the bundle
manifest, a file whose first bytes carry the IIDM namespace); a file
without the extension loaded fails with the instruction to create a
bundle.

## Creating a bundle

```sh
python3 -m venv ~/.venv-powsybl
~/.venv-powsybl/bin/pip install pypowsybl pandas
~/.venv-powsybl/bin/python tools/powsybl_dump.py network.xiidm cases/network.powsybl --case network
~/.venv-powsybl/bin/python tools/powsybl_dump.py builtin:ieee14 cases/ieee14.powsybl
```

`tools/powsybl_dump.py` writes every table with all attributes, runs
OpenLoadFlow for the reference solution (default parameters, except a
slack mismatch bound of 1e-4 MW and a Newton-Raphson tolerance of 1e-9
per equation so that a 1e-6 pu comparison has a converged reference),
writes `reference_buses.csv` and the manifest. `builtin:<factory>` names a
`pypowsybl.network.create_<factory>` example. CSV rules: comma separated,
strings always quoted, floats in full precision with `NaN` and `Inf` as
literals, booleans `true`/`false`, a missing integer as `-1`.

## Reading an IIDM file live

With `PythonCall` in the environment the extension loads on its own and
`import_case`, `run_sparlectra` and the Web UI accept `.xiidm` files.
Point PythonCall at a Python that has pypowsybl before Julia starts,
otherwise it installs its own through CondaPkg:

```sh
export JULIA_CONDAPKG_BACKEND=Null
export JULIA_PYTHONCALL_EXE=~/.venv-powsybl/bin/python
julia --project=. -e 'using Pkg; Pkg.add("PythonCall")'
```

```julia
using Sparlectra, PythonCall
tables = Sparlectra.read_powsybl_network("network.xiidm")   # OpenLoadFlow runs, tables with all attributes
Sparlectra.dump_powsybl_bundle("network.xiidm", "cases/network.powsybl")  # the same bundle the script writes
```

The live read runs the reference load flow with the tolerances of the
dump script, so its tables equal a bundle of the same file cell by cell.

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

| case | source | what it covers |
|---|---|---|
| `ieee14` | `pypowsybl.network.create_ieee14` | the IEEE 14-bus case: Y-bus identity with MATPOWER `case14` on every entry, including the two branches PowSyBl imports as lines between voltage levels (their one-sided shunts become bus shunts, see Conventions) |
| `ieee57` | `create_ieee57` | 57 buses, 17 transformers, two parallel transformer pairs, OpenLoadFlow's distributed slack |
| `four_substations` | `create_four_substations_node_breaker_network` | node-breaker topology with retained switches, a phase-shifting transformer, VSC and LCC HVDC links, an SVC, curve reactive limits, temporary limits, two synchronous components |
| `micro_grid_be` | `create_micro_grid_be_network` | CGMES origin, a three-winding transformer with a ratio tap changer, remote voltage regulation, dangling lines |
| `eurostag_tie_lines` | `create_eurostag_tutorial_example1_with_tie_lines_and_areas` | tie lines built from paired dangling lines |

Each fixture directory under `test/fixtures/powsybl/` holds the bundle,
the `.xiidm` file pypowsybl wrote (the example input for the live path)
and `reference_buses.csv`; `tools/gen_powsybl_fixtures.sh` regenerates
them.

## Conventions

Settled against the OpenLoadFlow voltages of the example files:

- **Transformer ratio.** Sparlectra's ratio is the MATPOWER tap, the no-load
  ratio of the from-side to the to-side voltage in pu. IIDM's `rho` is the
  ideal ratio at side 1 in kV and already carries `rated_u2 / rated_u1`, so
  `ratio = vn_to / (rho * vn_from)`; the phase shift is `-alpha`. The
  `_at_current_tap` impedances are used, no tap tables reach the network.
- **Magnetizing admittance.** The whole `g`, `b` of a transformer sits on
  the bus of side 1, behind the ideal transformer, as OpenLoadFlow places
  it; three-winding legs alike. Sparlectra's branch model carries one
  symmetric charging admittance, so this share, the one-sided excess of a
  line whose `b1` and `b2` differ, and the network-side admittance of a
  dangling line become **bus shunts** (`Sh_<bus>` entries next to the real
  shunt compensators). They are not tied to their branch: a branch outage
  in the contingency engine leaves them on the bus, an SCF or MATPOWER
  export writes them as bus shunts, and a state estimation that estimates
  shunts treats them as compensation. An outage study on an imported
  network should account for that; a branch-bound shunt is planned.
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
  reproduces OpenLoadFlow's balance. Without it the slack generator of
  each component (the regulating unit with the largest `max_p`, or a
  `slack_ids` override) absorbs the whole mismatch.
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
| `powsybl_import.python_exe` | `""` | Python executable for the PythonCall extension; empty means the PythonCall default. Read only when an IIDM file is imported live. |

The import report (counts per element type, every skipped element with
its reason, the slack decision per component, notices) is printed by
`format_powsybl_report` and travels with the imported case as
`provenance["powsybl_report"]`.

## Out of scope

Snapshots with missing injections and injection patching, time-series
batches, tap and voltage controllers from IIDM regulation data (the step
tables are dumped, not used), IIDM export, an IIDM XML reader in Julia, a
Web UI directory picker for bundles (a bundle is selected by its path).

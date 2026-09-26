# PowSyBl networks

Three example networks of the PowSyBl framework, each as the IIDM file
(`<case>.xiidm`) that Sparlectra reads directly and as the table bundle
pypowsybl wrote from it. No Python is needed for either:

```julia
using Sparlectra
run_sparlectra(casefile = "data/powsybl/ieee14.powsybl/ieee14.xiidm")   # the IIDM file
run_sparlectra(casefile = "data/powsybl/ieee14.powsybl")               # the bundle
```

| Directory | Network | What it exercises |
|---|---|---|
| `ieee14.powsybl` | IEEE 14-bus case (pypowsybl `create_ieee14`) | lines between voltage levels with unequal `b1`/`b2`, the Y-bus identity with MATPOWER `case14` |
| `four_substations.powsybl` | pypowsybl `create_four_substations_node_breaker_network` | node-breaker topology (retained switches as links), two synchronous components, an HVDC link, a phase-shifting transformer, an SVC |
| `micro_grid_be.powsybl` | CGMES MicroGrid BE (pypowsybl `create_micro_grid_be_network`) | a three-winding transformer, dangling lines, remote voltage regulation, a file state solved at other tap positions (the import starts flat and says so) |

## Your own `.xiidm` file

Select or upload it in the Web UI, or pass its path to `run_sparlectra`
or `import_case`. Sparlectra's IIDM reader (`read_iidm_tables`) resolves
what pypowsybl resolves: the bus-breaker and bus views of a node-breaker
substation, the components, the current tap step of every transformer,
the reactive limits on a capability curve, the operational limits. A
compressed file (`.xiidm.bz2`) is unpacked first (`bzip2 -d`). See the
documentation page [PowSyBl Import](../../docs/src/powsybl_import.md).

## What the bundle directory holds

- `<case>.xiidm`: the IIDM file, the source of everything else;
- `manifest.json`: the case name, the pypowsybl and IIDM versions, the
  list of tables with their row counts, and the OpenLoadFlow reference
  (iterations, reference bus per component);
- one CSV file per pypowsybl table (`buses.csv`, `lines.csv`,
  `2_windings_transformers.csv`, `generators.csv`, ... 25 tables, written
  with `all_attributes = true`), the state columns carrying OpenLoadFlow's
  solution;
- `reference_buses.csv`: the bus voltages OpenLoadFlow computed, which the
  test suite compares against (ieee14 within 1e-6 pu, the others within
  1e-5 pu), from the `.xiidm` file and from the bundle alike.

The bundle is the reference form: the test suite checks the reader against
it column by column. Nobody writes a CSV by hand: pypowsybl wrote the
bundles (every table with `all_attributes = true`, OpenLoadFlow with a
slack mismatch bound of 1e-4 MW and a Newton-Raphson tolerance of 1e-9
per equation for the reference). The three directories here are copies of
the test fixtures under `test/fixtures/powsybl`.

## Conventions the importer settled

The importer reproduces the OpenLoadFlow solution of every example: the
transformer ratio is `vn_to / (rho * vn_from)`, the phase shift `-alpha`,
a transformer's magnetizing admittance sits on its side 1 as the from
arm of the branch, a line keeps `b1` and `b2` on their own terminals, a
dangling line its admittance at the network terminal, HVDC stations are
fixed injections with the loss factors OpenLoadFlow reports, OpenLoadFlow's
slack distribution is carried as participation factors
(`power_flow.distributed_slack.p_mode = imported` reproduces it), and the
slack of a component is the largest locally regulating unit. The full list
with the measured alternatives is on the documentation page.

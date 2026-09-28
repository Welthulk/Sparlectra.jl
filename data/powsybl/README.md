# PowSyBl networks

Example networks of the PowSyBl framework as IIDM files (`.xiidm`).
Sparlectra reads the file itself, in Julia; no Python is needed:

```julia
using Sparlectra
run_sparlectra(casefile = "data/powsybl/ieee14.xiidm")
```

| File | Network | What it exercises |
|---|---|---|
| `ieee14.xiidm` | IEEE 14-bus case | lines between voltage levels with unequal `b1`/`b2`, the Y-bus identity with MATPOWER `case14` |
| `four_substations.xiidm` | four substations in node-breaker topology | retained switches as links, two synchronous components, an HVDC link, a phase-shifting transformer, an SVC |
| `micro_grid_be.xiidm` | CGMES MicroGrid BE | a three-winding transformer, dangling lines, remote voltage regulation, a file state solved at other tap positions (the import starts flat and says so) |
| `ieee14_sc.xiidm` | the IEEE 14-bus case with short-circuit data | the short-circuit run on a PowSyBl case: rated powers and the extension `generatorShortCircuit` on every generator |

`ieee14_sc.xiidm` is built by `tools/build_powsybl_sc_demo.jl` from
`ieee14.xiidm`, which stays as it is. The IEEE case defines neither rated
powers nor reactances; the values are assumed for demonstration
(`x''_d = 0.2` pu on the machine base; 300, 80, 60, 40 and 40 MVA).

The Web UI case selector offers the four files.

## Your own `.xiidm` file

Select or upload it in the Web UI, or pass its path to `run_sparlectra`
or `import_case`. The reader (`read_iidm_tables`) resolves the bus-breaker
and bus views of a node-breaker substation, the components, the current
tap step of every transformer, the reactive limits on a capability curve
and the operational limits. A compressed file (`.xiidm.bz2`) is unpacked
first (`bzip2 -d`). See the documentation page
[PowSyBl Import](../../docs/src/powsybl_import.md).

## Conventions the importer settled

The importer reproduces the OpenLoadFlow solution of the first three
files; the reference voltages are frozen in the test suite
(`test/fixtures/powsybl`). The transformer ratio is
`vn_to / (rho * vn_from)`, the phase shift `-alpha`, a transformer's
magnetizing admittance sits on its side 1 as the from arm of the branch, a
line keeps `b1` and `b2` on their own terminals, a dangling line its
admittance at the network terminal, HVDC stations are fixed injections
with their loss factors, the slack distribution is carried as
participation factors (`power_flow.distributed_slack.p_mode = imported`
reproduces it), and the slack of a component is the largest locally
regulating unit. The full list with the measured alternatives is on the
documentation page.

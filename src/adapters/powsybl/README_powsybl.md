# PowSyBl (IIDM) adapter

## Purpose

PowSyBl's native network format is IIDM (XML variant XIIDM). Sparlectra
does not parse IIDM: pypowsybl loads the file, resolves the node-breaker
topology and the tap positions, and hands over plain tables. This adapter
reads those tables and builds a `Net` from them. Two entry points share
one builder: a table bundle written by `tools/powsybl_dump.py` (no Python
at Sparlectra run time), and, when the PythonCall extension is loaded, an
`.xiidm` file read live.

## Bundle format

A bundle is a directory whose name ends in `.powsybl`, holding
`manifest.json` and one CSV per pypowsybl table (the getter name without
`get_`, plus `tie_lines`). The manifest names the case, the source file,
the pypowsybl and IIDM versions, every table with its file, row count and
index columns, and the OpenLoadFlow reference solution. The CSV rules and
the manifest fields are documented on the docs page
`docs/src/powsybl_import.md`; the reader is `read_powsybl_bundle`, the
writer `write_powsybl_bundle`, the column types come from
`POWSYBL_SCHEMA` and never from the data.

## Bus model

Every row of `get_bus_breaker_view_buses()` becomes one Sparlectra bus,
named by its bus-breaker id, with the bus-view id in the `external_id`
channel so results join back to the OpenLoadFlow bus table. Retained
switches become links (closed or open); non-retained switches never
appear, PowSyBl already contracted them. The nominal voltage comes from
the voltage level, the substation id is kept as bus metadata.

## Conventions

Tap changers are a fixed operating point: two-winding transformers use
the `_at_current_tap` impedances with `rho` and `alpha`, three-winding
legs the same columns per leg. No controllers and no tap tables reach the
`Net`; the step tables stay in the bundle for a later task.

The slack of a synchronous component is the regulating, connected
generator with the largest `max_p`, ties broken by id, unless
`powsybl_import.slack_ids` names one. Distributed slack remains a power
flow option.

Units follow pypowsybl (ohm, siemens, kV, degrees, MW, MVar, A) and are
converted to per unit on `powsybl_import.base_mva` and the nominal
voltage of the element's own side; transformers are referred through the
shared equivalent-circuit helpers. Flow sign: a PowSyBl branch column
`p1` is positive when power enters the branch at side 1. HVDC stations
are fixed injections by default: the rectifier draws
`target_p / (1 - loss_factor / 100)`, the inverter injects
`target_p * (1 - loss_factor / 100)`, as OpenLoadFlow reports them.

The two-winding transformer ratio convention was settled by the
acceptance test against the OpenLoadFlow voltages of the
`four_substations` case: `ratio = vn_to / (rho * vn_from)`, phase shift
`-alpha` (the convention comment at the head of `powsybl_mapping.jl`
carries the measured alternatives).

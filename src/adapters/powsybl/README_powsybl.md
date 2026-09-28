# PowSyBl (IIDM) adapter

## Purpose

PowSyBl's native network format is IIDM (XML variant XIIDM). Sparlectra
reads the file in Julia (`iidm_reader.jl`): it resolves the node-breaker
topology, the components and the tap positions the way pypowsybl does and
delivers the same plain tables; this adapter builds a `Net` from them. The
source is the `.xiidm` file. No Python is involved at any point.

## Reference bundles (test suite only)

The fixtures under `test/fixtures/powsybl` hold, next to each `.xiidm`
file, a table bundle: a directory whose name ends in `.powsybl`, with
`manifest.json` and one CSV per table (the getter name without `get_`,
plus `tie_lines`), the state columns carrying the OpenLoadFlow solution,
and `reference_buses.csv`. They are frozen references: the reader is
checked against them column by column, the power flow against their
voltages. `read_powsybl_bundle` and `write_powsybl_bundle` read and write
the form; the column types come from `POWSYBL_SCHEMA` and never from the
data.

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
`Net`; the step tables stay in the tables for a later task.

The slack of a synchronous component follows `reference_candidate_rank`:
a connected unit that regulates its own bus before one that regulates a
remote bus, then the largest (rated power where the file states one, else
`max_p`), unless `powsybl_import.slack_ids` names one. Distributed slack remains a power
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

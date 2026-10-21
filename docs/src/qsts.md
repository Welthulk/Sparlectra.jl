# Quasi-Static Time Series (QSTS)

QSTS solves a network over a sequence of operating states, one power flow per
step. No dynamics are modelled between steps.

| In scope | Out of scope |
|---|---|
| steady-state power flow per step | RMS dynamics, swing equations |
| comparison against delivered state variables (SV) | EMT, switching transients |
| SCF and CGMES snapshot series | time-coupled OPF, storage scheduling |
| | interpolation between snapshots |

## Terms

| Term | Meaning |
|---|---|
| snapshot | one complete operating state (CGMES EQ, SSH, TP, SV, or an SCF case) |
| series | snapshots of one model, consecutive in time, passing the series check |
| union network | one network containing every element of every snapshot of a series |
| delta | list of patch operations turning the union network into snapshot k |
| step | one entry of a series; solved once per run |

## Modes

| `qsts.mode` | Controller state | Steps | Use |
|---|---|---|---|
| `off` | | | default; `runQsts` refuses to run |
| `snapshots` | taken from each snapshot (SSH taps and sections) | independent, parallel | reference run with the set points of the source system |
| `sequential` | carried from step to step, the series values as reference | in order; several series in parallel | controller behaviour over a day: switching, hunting, violations |

`snapshots` runs without controllers unless `qsts.snapshot_controls = true`.
`sequential` runs the controllers of the `control` section.

## Workflow, CGMES

```julia
using Sparlectra

inputs = [CGMESSeriesInput(p) for p in sort(readdir("deliveries"; join = true))]
report = checkCGMESSeries(inputs)
printCGMESSeriesReport(stdout, report)

config = SparlectraConfig(Dict(
  "qsts" => Dict("mode" => "snapshots"),
  "cgmes_import" => Dict("tap_control" => true, "map_energy_source" => true),
))
union  = buildCGMESSeriesUnion(inputs; config = config)
result = runQsts(union; config = config)
writeQstsResults("out", result)
```

`tap_control = true` reads the SSH tap steps; with the default `false` the
taps are the `SvTapStep` positions of each snapshot. Deliveries that carry
generation as `EnergySource` need `map_energy_source = true`, otherwise that
generation is missing. `runQsts` refuses a configuration whose
`cgmes_import` keys differ from those the union was built with.

## Workflow, SCF

```julia
config = SparlectraConfig(Dict("qsts" => Dict("mode" => "snapshots")))
series = read_scf_series("data/qsts/sp_case14_day.qsts.zip")
result = runQsts(series; config = config)
```

The same series in order with the example controllers:

```julia
cfg = load_sparlectra_config("data/qsts/sp_case14_day.control.yaml"; overrides = Dict("qsts" => Dict("mode" => "sequential")))
result = runQsts(series; config = cfg)
qstsSwitching(result)          # step, element_id, old, new, steps, trigger
qstsMetrics(result)            # element_id, metric, value
```

## Series check

[`checkCGMESSeries`](@ref) groups snapshots into series. Each input is a
[`CGMESSeriesInput`](@ref): a delivery (zip, folder, files) or a loaded
`CGMESStore`, with an optional model key, scenario time and label. Without
them, the model key is the `modelingAuthoritySet` of the EQ header and the
time is the `scenarioTime` of the SSH header, else of the EQ header.
`writeCGMESSeriesReport("series.csv", report)` writes the report as CSV.

Computed per pair of consecutive snapshots of one model.

| Metric | Classes | Effect |
|---|---|---|
| structural overlap | ACLineSegment, PowerTransformer, PowerTransformerEnd, BusbarSection, Breaker, Disconnector, ConnectivityNode | below `min_structural_overlap`: new series |
| injection churn | SynchronousMachine, EnergyConsumer, EquivalentInjection, EnergySource, GeneratingUnit, PowerElectronicsConnection and its units | above `max_injection_churn`: pair marked |
| parameter drift | common objects with changed EQ parameters (r, x, b, g, ratio) | reported |
| time step | | above `max_gap`: new series |

Overlap is common objects over union. Objects with a new mRID but the same
name and connections count as common. Topological nodes are reported, not
counted: exports may renumber them. Switch states, SSH set points and the
stored energy of a battery are no drift.

The report also lists, per snapshot, closed non-retained switches whose two
terminals lie on different topological nodes (SSH and TP from different
switch states). The import turns such a switch into a link, so a solve
stays valid.

## Union network

[`buildCGMESSeriesUnion`](@ref) builds one network that holds every element
of every snapshot of a series and expresses each snapshot as a delta
against it.

| Item | Rule |
|---|---|
| node identity | set of terminals on the same topological node in every snapshot |
| switch | terminal pair; a link with one state per step |
| `_TPL` coupling | closes where a step joins two union nodes without a switch between them |
| missing unit | parked: zero injection, position kept |
| missing branch or shunt | status 0 |
| three-winding star point | addressed as `<PowerTransformer mRID>#star` |
| class change of one mRID | second object addressed as `<mRID>#<class>` |

A delta sets per step: values, branch parameters, transformer ratio and
shift, tap changer data (range, step, table ratios, impedance per step),
Q limits, regulation, reference role, tap controllers, bus isolation and
start voltages.

The link ops of a step are not all switch operations.
`n_switch_link_changes[k]` counts the switch links whose state in step `k`
differs from the template, `n_tpl_link_changes[k]` the `_TPL` couplings.
A `_TPL` coupling changes when an export puts terminals on other
topological nodes, for example equivalents on other connectivity nodes.

A switch that the TP keeps inside one topological node in every step is no
link, so nothing can open it. `extra_link_switches` (switch mRIDs, default
none) splits the union node at such a switch along its connectivity nodes
and the other switches of the node; the switch becomes a link that is
closed where the TP joins its sides, and a link op can open it. The split
needs node-breaker connectivity, and a second switch that joins the same
two sides has to be named as well.

```julia
u = buildCGMESSeriesUnion(inputs; config = cfg, extra_link_switches = [switch_mrid])
```

```julia
net = deepcopy(union.template)                       # the base snapshot
apply!(net, vcat(union.deltas[5].ops, union.deltas[5].control_ops), union.index)
runpf!(net, 50, 1e-8, 0)                             # snapshot 5
```

[`compareSnapshotReproduction`](@ref) applies a delta to a copy of the union
network and compares it with the direct import of that snapshot by mRID.
Empty result: snapshot reproduced.

Memory: the union network plus one snapshot. Series with
`cgmes_import.machine_control = true` or `cgmes_import.hvdc_mode =
:paired_control` are refused, because those controllers bind to element
positions.

## [Snapshot run](@id qsts_mode_a)

[`runQsts`](@ref) with `qsts.mode = snapshots`:

1. Solves the base case: for an SCF series the base case of the file, for
   a CGMES union the snapshot of its base step (`qsts.base_step`).
2. Solves every step on a fresh copy of the unsolved template with the
   full delta of the step applied. A step of an SCF series starts from the
   bus voltages of the base solution. A CGMES step starts from the state of
   its own snapshot: `start_ops` of its delta carry the start of every bus
   in service, the SV voltage or, for a bus without one, the neighbour
   start of the importer. A CGMES step applies `ops` and `control_ops`:
   shunt banks take the SSH sections, taps the SSH steps (`tap_control =
   true`) or the `SvTapStep` positions.
3. Compares every CGMES step with the SV state of its snapshot.

A step of an SCF series isolates a part of the net without a source as a
sequential step does (see [Sequential run](@ref qsts_mode_b)). A snapshot
run writes no `n_isolated` or `shed_load_mw`; one note per such step names
the number of buses and the load left unsupplied.

| Setting in base and steps | Value |
|---|---|
| `control.enabled` | `false`, unless `qsts.snapshot_controls = true` |
| start | SCF: the base solution; CGMES: the SV of the step, no fallback; one attempt per step |
| start generators, rescue paths, power mode, island parallelism | off |
| factorization reuse between steps | none |

Every key that differs from the passed configuration is listed in the
result `notes`. With `power_flow.flatstart = true` the base case starts
flat as in a power-flow run, without start projection or DC seed; the notes
name the switched-off keys. A `tap_pos` op of an SCF series sets the tap of its step,
also on a transformer with a tap controller; with the controllers running
such an op is refused. The run is refused, with the key named, for
`power_flow.solver` other than `rectangular`, `power_flow.mode = auto` and
`power_flow.external_grid.enabled = true`.

The steps run in parallel over Julia threads under `runtime.parallel` (the
keywords `parallel_enabled`, `parallel_max_tasks` and
`parallel_min_work_items` override it). Every step starts from a fresh copy
and nothing is carried between steps, so a parallel run gives the same
values and the same file bytes as a serial one.

## [Sequential run](@id qsts_mode_b)

[`runQsts`](@ref) with `qsts.mode = sequential` solves the steps of an SCF
series or a CGMES union in order. Each step starts from a fresh copy of the template with the
full delta of the step, and the outer-loop controllers run in every step.
Tap controllers (masters and followers) and shunt bank controllers carry
their state to the next step; all other controllers start every step from
the template.

1. The base case is solved with the controllers. Its solution is the start
   of step 1.
2. Every step: start voltages, full delta, the taps and banks of the series
   read as reference, the carried taps and banks written, one solve.
3. Step 1 is always a resync step. With `qsts.resync_every = n` every step
   with `(k - 1) mod n = 0` is one as well. A resync step keeps the taps and
   banks of the series and clears dead time and hunting memory.
4. A step that does not converge holds taps, banks and start voltages of the
   last solved step. Its dead time runs on.

| Setting | Value |
|---|---|
| `control.enabled` | as configured (`false`: no controllers, taps and banks follow the series) |
| `qsts.snapshot_controls` | no effect |
| start | `qsts.sequential_start = previous`: last solved step; `snapshot`: the same step in a snapshot run without controllers; a bus without a value there takes the base solution |
| start generators, rescue paths, power mode | off, as in a snapshot run |
| island parallelism | `runtime.parallel` for one series, off for several |

The carried tap is written absolutely: the model step of a modelled
transformer, else ratio and shift. A follower takes its own carried
position. A controller only carries the actuators it moves: an OLTC keeps
the shift of the series, a PST the ratio. On a transformer with a tap table
the base impedance is put on the values the delta wrote for the step, so the
table rescales the step's impedance when the tap moves. A `tap_pos` op on a
regulated transformer is the reference of its step.

**Start on a CGMES union.** Between snapshots of different days the state
can jump far (other switching, other angles), and Newton started from the
step before may not converge. For a union of different days set
`qsts.sequential_start = snapshot`; it starts every step from its own
snapshot run.

**Controller targets.** A target taken from the state of one snapshot (its
SV voltage) does not fit the switching state of the next: buses join or
split, and two controllers can then hold one node at different targets. A
sequential run with tap controllers needs the setpoints of the tap
changers, not values of the state.

**Dead time.** `qsts.dead_time_steps = d` locks a controller in the d steps
after a counted move. A locked controller stays in its rank, is evaluated
and moves nothing. While it stands outside its band and not at its limit it
holds every controller with a larger priority number, so a bank ranked
behind a locked tap waits; in its band or at its limit the larger numbers
act ([Controller priority](@ref controller_priority_ranks)). With `d = 0`
(default) nothing is locked.

**Positions.** A tap position is the model step, else the position counted
from the ratio 1.0 (`buildTapControllerReportRows`), for a PST from the
shift 0. A switched bank counts blocks (`bs_mvar / step_mvar`). A continuous
SVC has no position: it is carried, but gets no switching rows.

**Switching list.** One row per step and controller whose position changed
from the end of the step before to the end of the step (`qsts_switching.csv`,
[`qstsSwitching`](@ref)). A search back and forth inside one controller run
is not listed; it shows in `outer_iterations`.

| Trigger | When | Counted |
|---|---|---|
| `control` | the controller moved | yes |
| `range` | the position stood outside the control range before the solve | no |
| `resync` | a resync step | no |

**Hunting.** Three counted moves of one controller with alternating sign
whose first and last step lie at most `qsts.hunting_window - 1` steps apart
are one pattern. Overlapping patterns count each (up, down, up, down: two).
A resync step breaks a pattern.

On a CGMES union the taps and banks of each snapshot (SSH steps and
sections, `control_ops`) are the reference of their step; a step compares
with its SV as in a snapshot run.

In a step of an SCF series or an event series (also in a snapshot run), a
part of the net without a source (no unit with P > 0 or voltage regulation; a load is never a source)
is isolated, as the CGMES import does: its injections are zeroed, its
branches go out of service, and its buses have no values in that step. The
step counts them in `n_isolated` and their load in `shed_load_mw`. A CGMES
delta carries its isolation itself.

Several series: `runQsts([s1, s2, s3]; config)` returns one result per
series (SCF series or CGMES unions). The series run in parallel over Julia
threads, the steps of one series in order; parallel and serial runs give
the same values and file bytes. `writeQstsResults(dir, results)` writes
each result into `dir/series_<k>`, its manifest with `series` (k) and
`n_series`.

## [Event list](@id qsts_events)

A sequential run can be driven by switching events instead of snapshot
deltas. [`readQstsEvents`](@ref) reads the list, [`qstsEventSeries`](@ref)
turns it into a series on an SCF case or a CGMES union:

```
step,element_id,action,value
3,26,open,
5,32,out,
9,26,close,
```

| Column | Rule |
|---|---|
| header | exactly `step,element_id,action,value` |
| `step` | whole number from 1 |
| `element_id` | SCF component id, or the mRID of a union element (a switch must be a link of the union, see `extra_link_switches`) |
| `action` | `open`, `close` on a branch, transformer or switch; `in`, `out` on a unit, load or shunt |
| `value` | empty |

An element takes one event per step and keeps its state until its next
event. Step k of the series is the full delta of the switching state after
all events up to k; a unit goes out by parking. Every error names the CSV
line.

```julia
events = readQstsEvents("events.csv")
series = qstsEventSeries(read_scf_json("data/scf/sp_case14.scf.json"), events)
config = SparlectraConfig(Dict("qsts" => Dict("mode" => "sequential")))
result = runQsts(series; config = config)
qstsValues(result, "", "n_islands")
```

An event series has no time reference and no SV. A snapshot run refuses it,
and so does `qsts.sequential_start = snapshot`. An SCF event series keeps
the source file and the row count in `provenance["event_list"]`, and the
manifest of the run records them. On a CGMES union `qstsEventSeries`
returns a [`CGMESEventSeries`](@ref); its base is the union template.

## Chart

[`qstsChart`](@ref) draws one quantity of one element over the steps as SVG
text: one point per step with a value and a tooltip, the line broken at a
step without a value, the band dashed (`band = :auto`: `qsts.vm_min_pu` and
`qsts.vm_max_pu` for `vm_pu`, `qsts.loading_max_pct` for `loading_pct`). The
SVG holds no script; a notebook or the documentation shows the returned
[`QstsChart`](@ref) directly.

```julia
c = qstsChart(result, "8", "vm_pu")
write("moorau_vm.svg", c.svg)
```

## Results

`writeQstsResults(dir, result)` writes:

| File | Content |
|---|---|
| `qsts_values.csv` | long format `step,element_id,quantity,value` |
| `qsts_steps.csv` | per step: label, timestamp, base flag, status, control status, iterations, outer iterations, residual, frozen buses, slack buses, message |
| `qsts_elements.csv` | `element_id, kind, name, from, to, vn_kv, rating_mva` |
| `qsts_manifest.json` | format version, mode, step count, base step, step length, solver settings, notes; `sequential`: the keys of mode B and `n_failed`, an event series its `event_list` (source, rows) |
| `qsts_switching.csv` | `sequential` only: `step,element_id,old,new,steps,trigger` |
| `qsts_metrics.csv` | `sequential` only: `element_id,metric,value` |

The `kind` of a branch in `qsts_elements.csv` is `line`, `transformer` or
`equivalent_branch` (a CGMES `EquivalentBranch`, whose `p_loss_mw` can be
negative and is reported as computed).

Long format, quantities:

| Element | Quantities | Unit |
|---|---|---|
| bus | `vm_pu`, `va_deg` | pu, deg |
| branch | `p_from_mw`, `q_from_mvar`, `p_to_mw`, `q_to_mvar`, `p_loss_mw`, `q_loss_mvar`, `loading_pct` | MW, MVAr, % |
| transformer with a tap controller | `tap_ratio`, `shift_deg` | 1, deg |
| step | `converged`, `iterations`, `outer_iterations`, `n_frozen`, `total_p_loss_mw`, `total_q_loss_mvar`, `n_low_vm` | |
| bus, branch, step | SV comparison, see below | pu, deg, MW, MVAr |
| step, `sequential` | `resync`, `n_moves`, `n_islands`, `n_isolated`, `shed_load_mw`, `max_dvm_step_pu`, `max_dp_step_mw`, `n_vm_low`, `n_vm_high`, `n_overload` | count, MW, pu |
| controller, `sequential` | `position`, `bs_mvar`, `ref_position`, `dpos_ref`, `moves`, `at_limit`, `locked` | steps, MVAr, 0/1 |

Step quantities of a sequential run: `n_moves` is the sum of the counted
steps, `n_islands` the islands of the step (a sourceless one included),
`n_isolated` the buses without a value because they are isolated and
`shed_load_mw` their load, `max_dvm_step_pu` and `max_dp_step_mw` the largest
change against the step before (only when that step converged),
`n_vm_low`, `n_vm_high` the buses outside `qsts.vm_min_pu` and
`qsts.vm_max_pu`, `n_overload` the branches above `qsts.loading_max_pct`.
`ref_position` is the position the series gives the step, `dpos_ref` the
simulated minus the reference position in steps, `locked` 1 when the dead
time locked the controller. A controller has the address `<actuator>#<kind>`: the
transformer, the bank's shunt when the series addresses it, else the bank's
bus, and `OLTC`, `PST`, `SCHRAEG`, `MSC` or `SVC`. In `qsts_elements.csv` it
has kind `controller`, the actuator in `from` and the regulated bus in `to`.

Metrics of a sequential run (`qsts_metrics.csv`, [`qstsMetrics`](@ref)):

| Element | Metrics |
|---|---|
| run (`element_id` empty) | `n_steps_solved`, `n_failed`, `switch_count_total`, `n_hunting_total`, `n_resync`, `loss_energy_mwh` |
| controller | `switch_count`, `n_hunting`, `steps_at_limit`, `steps_locked`, `max_abs_dpos_ref`, `rms_dpos_ref` |
| bus outside the band | `vm_low_episodes`, `vm_low_steps`, `vm_low_minutes`, `vm_low_max_pu` and the same for `vm_high` |
| branch above the limit | `overload_episodes`, `overload_steps`, `overload_minutes`, `overload_max_pct` |

An episode is a run of consecutive steps in violation; a step without a
value ends it. `_max_` is the largest excess beyond the limit. The step
length is `step_minutes`, else the distance to the next timestamp; without
either, the `_minutes` rows and `loss_energy_mwh` are left out.
`loss_energy_mwh` sums `total_p_loss_mw` times the step length over the
solved steps.

Rows exist only for solved values: no rows for isolated buses (a bus
that hangs on the network by a closed link or switch alone is solved with
its link group and has rows), buses in an island without slack and PV bus (counted in `n_frozen`), star points,
branches out of service, or non-converged steps (those keep only the step
quantities of the solve). `loading_pct` needs a rating. Units are part of
the quantity name. Format version 1.0; other columns, a renamed quantity,
another unit or another rule for addresses or gaps raise it to 2.0, new
quantities and element kinds do not.

`status` in `qsts_steps.csv` is `solved`, `not_converged` or
`solver_error`. `not_converged` covers the whole net and the AC islands
of a net with several; for islands `iterations` holds the Newton
iterations of the first island that did not converge, and `message` names
each such island with its status, iterations, reason and final mismatch.
`solver_error` is a solve that threw, with the first line of the error in
`message`. `residual` counts differences of a CGMES snapshot the delta
cannot express; the slack buses are separated by `|` (the slack can move
between steps). The CSV files follow `output.csv_format` with 15
significant digits.

A sequential run adds controller elements, step and controller quantities
and two files. Mode A writes none of them, and its files are unchanged.

Every converged solve gets the plausibility check of a single power flow
(`power_flow.plausibility.vmin_warn`, default 0.7 pu, see
[Configuration](powerflow_configuration.md)). `n_low_vm` counts the
energized buses below the threshold; it has no row when the step did not
converge or the threshold is 0. A base case or step with such a bus stays
`solved`; it gets one entry in `notes` and one warning in the log, both
with the step, the lowest bus and its voltage.

Element addresses: on a CGMES union the mRID (union topological node, line,
two-winding transformer, transformer end of a three-winding leg, the star
and class-change addresses above). Branch addresses are stable across runs;
node addresses follow the base step. On an SCF series the component id.

```julia
using DataFrames
df = DataFrame(qstsLongTable(result))
v  = qstsValues(result, "BUS_17", "vm_pu")   # NaN where no value
```

## SV comparison

For CGMES series each step is compared with its delivered SV.

| Quantity | Element | Step metrics |
|---|---|---|
| `sv_dvm_pu`, `sv_dva_deg` | bus | `sv_n_nodes`, max and RMS, `sv_va_offset_deg`, `sv_tn_spread_pu` |
| `sv_dp_from_mw` ... `sv_dq_to_mvar` | branch terminal | `sv_n_terminals`, max and RMS |

Angles are aligned per electrical island (median). The metrics count each
topological node once and equal those of [`compareWithSV`](@ref) on the
solved direct import; `sv_tn_spread_pu`, the largest voltage difference
between the union buses of one node, should be 0. Flows compare every
`SvPowerFlow` at a line or transformer terminal. Zero-impedance links
cannot reproduce voltage differences carried by the SV across closed
couplers; this sets a lower bound on achievable deviation.

## Measurement file

One format for state estimation input, independent of the source system.
CSV, UTF-8, comma, decimal point, LF.

```
timestamp,measurement_mrid,terminal_mrid,node_mrid,type,value,sigma,source,valid
2026-01-15T10:00:00Z,AN_101,,,,231.4,0.5,scada,true
2026-01-15T10:00:00Z,,T_55,,P,-42.1,,scada,true
2026-01-15T10:00:00Z,,,TN_12,U,412.0,1.0,scada,false
```

| Column | Rule |
|---|---|
| key | exactly one of `measurement_mrid`, `terminal_mrid`, `node_mrid` |
| `measurement_mrid` | mRID of a CGMES `Analog`; type and location from EQ, `type` empty |
| `terminal_mrid` | `type` one of `P`, `Q`, `I` |
| `node_mrid` | `type` `U`; topological or connectivity node |
| `timestamp` | ISO 8601 UTC; one measurement set per timestamp; empty for a single snapshot |
| `sigma` | SI unit of the type; empty: SE default |
| `valid` | `false` rows are counted, not used |

## Synthetic series

| Tool | Input | Output |
|---|---|---|
| `tools/cgmes_ssh_series.jl` | one CGMES snapshot, role file (TOML), steps, step length, seed, optional event list | one bundle per step and `manifest.toml`; EQ, TP, boundary unchanged, SSH varied |
| `tools/build_scf_series.jl` | SCF case | SCF series zip; without `--write` it checks the shipped file against a fresh build |

Roles: `load` (daily profile), `wind` (ramp), `thermal` (balance),
`fixed` (unchanged, default). Taps and sections stay fixed. A `[range]`
table sets the factor range of a role or a single class key and the bounds
of the thermal balance; where a bound binds, the rest stays at the slack
and the generator notes it per step. A machine can be keyed by the class
of its generating unit (`"SynchronousMachine:WindGeneratingUnit"`), a
`PowerElectronicsConnection` by the class of its units
(`"PowerElectronicsConnection:PhotoVoltaicUnit"`); a connection with units
of several classes takes the plain key and gets a note.

```toml
EnergyConsumer = "load"
EnergySource = "load"
EquivalentInjection = "wind"
"SynchronousMachine:ThermalGeneratingUnit" = "thermal"

[range]
EnergySource = [0.8, 1.2]
thermal = [0.5, 1.5]
```

Defaults: `load` `[0.6, 1.2]` (night to evening peak), `wind` `[0.0, 1.3]`
(ramp scaled onto the range), noise clipped to the range. Without a thermal
range a step that needs a negative thermal factor stops the run. If the
mapping has thermal units, the manifest gives `thermal_factor` and
`imbalance_mw` per step, and the command line prints the steps where a
bound binds. A range that reaches no object, for example `range.thermal`
without thermal units, has no effect and is named in the notes of every
bundle. A mapping without `[range]` writes the same bytes as before the
table existed.

```
julia --project=. tools/cgmes_ssh_series.jl base.zip boundary.zip out --roles roles.toml --steps 96 --minutes 15
```

An event list (CSV with the header `step,element_id,action,value`, action
`open` or `close` on a switch) opens and closes switches; the TP is not
recomputed, so the series check names those steps.

## [SCF series zip](@id scf_time_series)

An SCF series is a base case plus one delta per step, addressed by SCF
component ids. It needs no CGMES and no union network. The zip holds
exactly three entries; any other entry is an error.

| Entry | Content |
|---|---|
| `series.json` | format `sparlectra-scf-series`, `format_version` 1.0, `name`, `created_by`, `step_minutes` (`null` for irregular steps), `steps`, free `provenance`; a step has its number (1 to n, no gaps), an ISO 8601 UTC `timestamp` (on every step or on none) and an axis `label` |
| `base.scf.json` | an unchanged case file |
| `delta.csv` | one patch operation per row, sorted by step |

`delta.csv` is UTF-8 with comma, decimal point and LF; `#` comments are
allowed above the header only:

```
step,op,target,id,field,value
1,scale,load,32,,0.507
1,scale,generator,41,p,0.5178
57,status,branch,26,,0
```

The columns are those of a [`PatchOp`](@ref): `set` with a field and a
value, `scale` with the factor as value and an optional field `p` or `q`,
`status` with 0 or 1. Units are those of the scenario block (MW, MVAr, pu,
degrees). Each step is the full delta against the base case and `scale`
refers to the base value, so every step applies on its own. A step writes
each field of a component once. Node ids address buses; tap controllers
have no id in a case file, so a series cannot address them.

[`read_scf_series`](@ref) checks the whole file before the first step, and
every error names the entry and, for `delta.csv`, the line.
[`is_scf_series`](@ref) decides by content; a power flow on a series zip
(`run_sparlectra`, the services, the Web UI) runs the base case; the
**QSTS** button of the Web UI runs the series
([Web UI Reference](@ref webui_qsts)).
`write_scf_series("copy.qsts.zip", series)` writes a series back.

**Units out of service in the base case.** A `sym_load`, `sym_gen` or
`source` with `status: 0` in the base case is parked when the series is
loaded: `build_net(series)` turns it into a 0 MW load in place
(`set in_service 0`), and a step brings it back with `set in_service 1`.
A `set`, `scale` or `status 1` on a parked unit needs `set in_service 1`
in the same step; `status 1` alone does not bring it back. A unit that a
step parks with `set in_service 0` takes no other operation in that step.
A plain case file read still builds such a unit in service, and
`status 0` in a step or an N-1 scenario still removes the unit.

**Example.** `data/qsts/sp_case14_day.qsts.zip` holds 96 steps of 15
minutes on `sp_case14`. Loads follow a daily shape by voltage level, the
two machines other than the reference follow the load, and the line
`Ostheim_Neuwiese_1` is out from 14:00 to 15:45. Taps and the shunt stay
as in the base case.

For runs with controllers, `data/qsts/sp_case14_day.control.yaml` adds a
tap controller for `Wehrden_110_20` and makes the shunt at `Moorau_20` a
switched bank of 2.5 MVAr blocks (0 to 10 MVAr in total), next to the
tap controller of `Moorau_110_20` from the case file. It also gives that
controller the control range `tap_min: 1.0`, so the tap does not go below
the neutral ratio, and puts the bank behind it (`priority: 1`, target
1.01 pu, band 0.02 pu). The bank band of 0.99 to 1.03 pu contains the tap
band, so the bank switches only when the tap stands at its limit and the
voltage falls below 0.99 pu. In a snapshot run with
`qsts.snapshot_controls = true` that happens from 18:00 to 20:30: the Moorau
tap stands at 1.0 and the bank switches in one block of 2.5 MVAr. In all
other steps the bank stays off. In a sequential run the bank switches in at
18:00 and stays in to the end of the day, since the voltage stays inside the
bank band once the block is in.

## Configuration

| Key | Default | Meaning |
|---|---|---|
| `qsts.mode` | `off` | `off`, `snapshots`, `sequential` |
| `qsts.base_step` | 1 | base snapshot of a CGMES union |
| `qsts.snapshot_controls` | `false` | controllers per step in `snapshots` |
| `qsts.sequential_start` | `previous` | start of a step in `sequential`: `previous`, `snapshot` |
| `qsts.resync_every` | 0 | resync every n steps; 0: step 1 only |
| `qsts.dead_time_steps` | 0 | steps a controller stays locked after a counted move |
| `qsts.hunting_window` | 3 | steps of a hunting pattern, at least 3 |
| `qsts.vm_min_pu`, `qsts.vm_max_pu` | 0.9, 1.1 | voltage band of the metrics |
| `qsts.loading_max_pct` | 100 | loading limit of the metrics |
| `qsts.series.min_structural_overlap` | 0.95 | below: new series |
| `qsts.series.max_injection_churn` | none | above: pair marked |
| `qsts.series.max_gap` | none | minutes; above: new series |
| `cgmes_import.map_energy_source` | `false` | map EnergySource as PQ injection |
| `cgmes_import.map_power_electronics` | `false` | map PowerElectronicsConnection |

All keys accept case-level overrides.

## Limits

- Steps are stationary states. In `snapshots` they are independent; in
  `sequential` only the tap and bank positions and the dead time link them.
  Time between steps carries no physics.
- The CGMES import de-energizes islands without a source, also those whose
  units only draw power. An island with infeed but no reference gets a
  slack from the importer: among the units with infeed the strongest
  `referencePriority`, else the largest active power reserve, else the
  largest infeed, buses of degree greater than 1 first; the import messages
  name each choice.
- A slack bus keeps its start angle, and an island without slack keeps
  the start values of its buses. The deltas of the series union carry both,
  so every step reproduces the direct import.
- Converged is not plausible: where available, check the SV comparison.

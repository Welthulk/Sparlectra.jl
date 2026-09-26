# CGMES Import

Sparlectra reads ENTSO-E **CGMES 2.4.15** deliveries (Common Grid Model
Exchange Standard) and builds a bus-branch `Net` that solves with `runpf!`.
The importer reports everything it skips and validates against the
delivery's own solved state.

## Quick start

```julia
using Sparlectra

# 1. Diagnose first: works on incomplete or broken deliveries
summary = summarizeCGMES(path = ["grid_EQ_SSH_TP_SV.zip", "boundary.zip"])
show(stdout, MIME"text/plain"(), summary)

# 2. Import (Net only)
net = createNetFromCGMES(path = ["grid.zip", "boundary.zip"])

# 3. Import with full result: Net, store, topology, short-circuit data, notices
result = importCGMES(path = ["grid.zip", "boundary.zip"])
runpf!(result.net, 30, 1e-8, 0)

# 4. Validate against the delivery's own SV profile
cmp = compareWithSV(result)
@show cmp.max_dvm cmp.max_dva
```

`path` accepts a folder, a single ZIP, or a vector of either; nested ZIPs
are opened in memory. Profiles are classified from the `md:Model` header.

## What is imported

| CGMES source | Result |
|---|---|
| `TopologicalNode` (TP) | buses, nominal voltage from `BaseVoltage`; without a TP profile the buses are derived from `ConnectivityNode`s and switch states (see [the topology processor](@ref topology_processor)) |
| `ACLineSegment`, `SeriesCompensator` | π-model branches; lines spanning two nominal voltages become ratio branches (boundary lines) |
| `PowerTransformer` (2 and 3 windings) | π-model transformers; 3-winding as star equivalent with an AUX bus |
| `RatioTapChanger`, phase tap changers | fixed tap positions, or outer-loop controllers with `tap_control = true`; `PhaseTapChangerTabular` resolves its `PhaseTapChangerTable` row (ratio and angle) at the tap position. Tap angles fold with end-referral semantics: an end-2 (to-side) angle enters negated (`θ_eff = θ1 − θ2`), matching RealGrid's SV state; the ENTSO-E PSEI `PTE2` conformity toys expect the unflipped angle and deviate by about 0.3° |
| retained `Switch`/`Breaker`/… | zero-impedance bus links |
| `LinearShuntCompensator`, `NonlinearShuntCompensator` | shunts; the nonlinear characteristic sums its per-section points up to the active `sections` count (each point read as one switched group; CIM would also permit reading a point as the absolute value at that section count) |
| `EnergyConsumer`, `ConformLoad`, `NonConformLoad`, `StationSupply` | loads (SSH values) |
| `SynchronousMachine`, `ExternalNetworkInjection`, `EquivalentInjection` | injections, PV where a local voltage control is active; a machine whose voltage `RegulatingControl` points at a *different* bus is held PV at its own bus by default, or becomes an outer-loop remote voltage controller with `machine_control = true` (see [Remote Voltage Control](remote_voltage_control.md)); machine Q limits from the `ReactiveCapabilityCurve` (evaluated at the scheduled P) where one exists, else the scalar `minQ`/`maxQ` hull; a positive `GeneratingUnit.normalPF` arrives as `ProSumer.participationFactor` for the [distributed slack](powerflow_configuration.md) (`p_mode: imported`), zero/absent maps to unknown |
| `StaticVarCompensator` | P = 0 reactive injection, Q limits from the Ω ratings |
| `AsynchronousMachine` | fixed PQ operating point from the SSH `RotatingMachine.p`/`q` (load convention: a motor consumes with `p > 0`); no voltage regulation, no slack candidacy. The induction machine's voltage/slip dependence is a dynamics topic, not a steady-state one |
| `SvVoltage`, `SvPowerFlow`, `SvTapStep` | start values and validation reference |

Everything else is counted in the coverage report; `summarizeCGMES` shows
the full class histogram. Short-circuit source data is harvested on every
import (read, not evaluated); `shortCircuitCoverage` reports its
per-attribute completeness, printed in `cgmes.log`. Theory and data
table: [Short-Circuit Compendium](short_circuit.md).

## Conventions

| Convention | Sparlectra behaviour |
|---|---|
| Sign convention | CGMES uses the load convention (a machine with `p < 0` injects). Machine, external-injection and equivalent-injection values are inverted; a load with negative `p` passes through unchanged. |
| `Terminal.connected`, both ends open | Branch out of service. |
| `Terminal.connected`, exactly one end open | No longitudinal current; the charging arm of the closed terminal stays as a shunt at the closed bus (it belongs to the reactive balance). |
| `PowerTransformerEnd.g`, `.b` | Each end's magnetizing admittance stays on its own terminal of the branch (end 1 referred to the end-2 base as the from arm, end 2 as the to arm; a three-winding end on the leg's bus side). The export writes them back per end, so a round trip keeps the placement. A line's `bch`, `gch` is one value and splits in halves. |
| `Terminal.connected`, both ends closed | Normal branch. |
| Boundary sets | Cross-border deliveries reference topological nodes that live in the boundary files. Detected through unresolved references, not filenames; the import fails with an explicit message unless `require_boundary = false`. |
| Import analysis | On abort the importer first prints the supplied model files (profile, version, model id), every `md:Model.DependentOn` prerequisite matched against the supplied models (a missing boundary set is named by its model id), an unresolved-reference histogram by class and property, and a verdict. In Web UI runs it lands in `cgmes.log`; the **Analyze import** button runs the same check before a full import (see [Web UI](webui.md)); `analyzeCGMES(path = ...)` does it on demand. |
| Unresolved `TopologicalNode.BaseVoltage` | Cannot be repaired with `require_boundary = false`: real ENTSO-E deliveries keep their base-voltage catalog in the boundary EQ file (see `infer_base_voltages` below). |
| Slack selection | In order: `referencePriority ≥ 1` (a notice on ties), a single `ExternalNetworkInjection`, the largest external injection, the largest synchronous machine. |
| Sourceless parts | Components without any generator are de-energized (injections zeroed, buses isolated) with one message each, instead of aborting the power flow. |

## [Node-breaker deliveries without a TP profile](@id topology_processor)

Some EMS and substation-level exports ship EQ+SSH only, topology as
`ConnectivityNode`s plus switches. The importer then derives the bus
partition itself.

- **Trigger.** Only when the delivery contains `ConnectivityNode`s and no
  non-boundary `TopologicalNode`; a TP-carrying delivery is unaffected.
- **Aggregation.** Connectivity nodes merge across closed, non-retained
  switches (`Switch`, `Breaker`, `Disconnector`, `LoadBreakSwitch`,
  `Jumper`, `ProtectedSwitch`). A switch counts as open when SSH
  `Switch.open` says so (overriding EQ `normalOpen`) or when it is out of
  service. Retained switches are never merged; they become zero-impedance
  bus links.
- **Result.** Each connectivity group becomes a synthetic
  `TopologicalNode`, named after a busbar section in the group where one
  exists. Nominal voltage resolves through the container chain
  (`ConnectivityNode` to `VoltageLevel`, hopping over `Bay` containers) to
  the `BaseVoltage`. Boundary connectivity nodes adopt the boundary set's
  existing TP_BD nodes, so cross-border stitching works as with a shipped
  TP.
- **Visibility.** The import messages announce the processor
  (`topology processor: derived N topological node(s) from M connectivity
  node(s) ...`); without that message the shipped TP was used. A delivery
  with neither a TP profile nor connectivity nodes aborts with the import
  analysis.

On the ENTSO-E conformity sets the derived partition reproduces the
shipped TP bus for bus (MiniGrid, SmallGrid). Two data findings: FullGrid's
TP assigns the terminals of one connectivity node (the shared load node of
`BE-Load_1`/`BE_CL_1`/`BE_NC_1`) to two topological nodes, which no
processor can derive from the graph, so 10 MW of load sit one bus apart;
and dead network fragments (all equipment out of service) may be
partitioned differently than the sender's TP, showing up as extra isolated,
de-energized buses.

## [Configuration (`cgmes_import`)](@id cgmes-import-config)

The keys that steer a CGMES import: where the delivery is, the system base, whether unresolved boundary references fail the import, and the per-element handling described below. The Web UI carries the ones a delivery is usually adjusted with.

| Key | Default | Purpose |
|---|---|---|
| `cgmes_import.path` | `""` | Delivery location(s); `;`-separated for multi-part deliveries (base case plus boundary set). |
| `cgmes_import.base_mva` | `100.0` | System base in MVA (CGMES does not define one). |
| `cgmes_import.require_boundary` | `true` | Fail when topology references stay unresolved. Allowed values: `true`, `false`. |
| `cgmes_import.tap_control` | `false` | Start from the SSH tap positions and attach the CGMES-defined outer-loop tap controllers, instead of importing the solved `SvTapStep` positions as fixed taps. Allowed values: `true`, `false`. |
| `cgmes_import.machine_control` | `false` | Attach outer-loop remote voltage controllers (`MachineVoltageControl`) for machines whose voltage `RegulatingControl` points at a different bus, instead of holding those machines PV at their own bus. Allowed values: `true`, `false`. |
| `cgmes_import.ignore_connected` | `false` | Diagnostic override treating every terminal as connected, for snapshots whose SSH flags contradict their own SV state. Allowed values: `true`, `false`. |
| `cgmes_import.vset_min_pu` | `0.5` | Lower bound of the plausibility band for a voltage `RegulatingControl.targetValue`, in p.u. of the regulated bus's nominal voltage. |
| `cgmes_import.vset_max_pu` | `1.5` | Upper bound of that band. A target outside `[vset_min_pu, vset_max_pu]` is treated as a placeholder: it is ignored, the unit is held PV at the bus voltage derived from the nominal data, and the substitution is reported as a `warning:`. |
| `cgmes_import.multi_slack` | `true` | Give every electrical island its own SV-declared angle reference (at most one per island). Required for multi-island deliveries; `false` forces the single-reference behavior. Allowed values: `true`, `false`. |
| `cgmes_import.placeholder_guards` | `warn_skip` | Behavior of the placeholder guards (implausible shunt admittances, tap corrections outside 0.5 … 2.0). `warn_skip`: keep the filler value out of the solve with a warning. `strict`: abort the import with an error naming the object, for deliveries where dropped data must never go unnoticed. Allowed values: `warn_skip`, `strict`. |
| `cgmes_import.infer_base_voltages` | `false` | Reconstruct missing nominal voltages when the delivery ships without its `BaseVoltage` catalog (in real ENTSO-E deliveries the catalog lives in the boundary EQ): nodes are seeded from the SV voltages (kV, snapped to the standard level series) and transformer rated voltages, then the levels propagate across level-preserving equipment (anything but a transformer). All substitutions are summarized as one `warning:` message with per-source and per-level counts. Pair with `require_boundary: false`; nodes that stay unresolved still abort with the import analysis. Allowed values: `true`, `false`. |
| `cgmes_import.hvdc_mode` | `injections` | HVDC converter handling. `injections`: fixed PCC injections with the SSH operating point, the industry-standard load-flow treatment; areas joined only through HVDC stay separate islands. `paired_control`: the same injections plus one steerable `HvdcPairControl` per detected converter pair; the DC topology (`ACDCConverterDCTerminal`, `DCNode`, `DCLineSegment`) groups the converters into links, the transfer and loss are derived from the two SSH operating points, and the outer loop keeps the pairing invariant exact. Detection runs in both modes and names every pair in the import messages; a pair that cannot be attached (skipped converter, inconsistent snapshot) degrades to the fixed injections with a notice. Allowed values: `injections`, `paired_control`. |
| `cgmes_import.start_values` | `auto` | Newton-Raphson start state for CGMES runs. `auto`: use the delivery's `SvVoltage` state when it carries one, else the flat start (a real delivery is built around its own operating point, and starting elsewhere makes it diverge for no good reason). `sv`: always start from the imported state. `flat`: always use the synthetic flat start, for method studies. With an SV start the competing start-value machines (`start_projection`, `dc_seed_unconditional`, `start_current_iteration`, `apslf_start`) are forced off for the run. On CGMES runs an explicit `sv` or `flat` wins over `power_flow.flatstart`; under `auto` a set `power_flow.flatstart` chooses the flat start although the delivery carries an SV state (the Settings page offers that one switch). MATPOWER and DTF runs ignore this key. The resolved decision (including what `auto` chose and any overridden keys) is logged to `run.log` and `cgmes.log`; the SV comparison (`sv_compare.csv`) runs in every mode. Allowed values: `auto`, `flat`, `sv`. |

### Implausible voltage setpoints

A voltage `RegulatingControl.targetValue` outside the
`vset_min_pu`/`vset_max_pu` band is treated as a placeholder (see the
table). Real targets in the tested deliveries span 0.92 … 1.15 p.u. Widen
the band if a delivery legitimately regulates outside it, or set
`vset_min_pu: 0.0` with a large `vset_max_pu` to accept every value.

!!! details "Why the band exists"
    Some deliveries carry placeholder regulation targets: ReliCapGrid's
    Svedala model declares `targetValue = 0.001` kV (about `5e-5` p.u.)
    on 17 kV and 20 kV generator busbars. Taken literally, such a target
    turns the reference bus of its island into a bus at zero volts and
    the power flow converges onto the all-zero solution. The Svedala
    units also carry `Equipment.inService = false` and are skipped
    earlier; the band guards deliveries that park units without setting
    `inService`.

### Placeholder guards: shunt admittances and tap corrections

Two guards catch filler values; `cgmes_import.placeholder_guards: strict`
aborts with an error naming the object instead of skipping.

- a shunt whose admittance exceeds **10 × baseMVA** at nominal voltage is
  skipped with a `warning:` (not clamped: a placeholder carries no
  information to clamp to);
- a single tap correction factor outside **0.5 … 2.0** is ignored with a
  `warning:`; the transformer keeps its nominal ratio (real tap ranges
  stay within a few ten percent of neutral).

!!! details "Why the guards exist"
    Conformity sets such as FullGrid fill attributes with the `X.99`
    scheme (tabular PST row `ratio 9.99 / angle 0.99°`, a
    `NonlinearShuntCompensatorPoint` with `b = g = 0.99 S`, a 50-GW shunt
    at 225 kV, switch `ratedCurrent 999.99`, …). With the guards FullGrid
    solves from a flat start; its shipped SV profile is internally
    inconsistent (a 14.5° angle jump across a 0.3 Ω line), so the SV start
    and comparison are meaningless there.

### Machine Q limits: the `ReactiveCapabilityCurve` is Q(P), not Q(U)

A `ReactiveCapabilityCurve` gives reactive limits as a function of active
power: the `CurveData` x axis is the machine's own P in the CGMES machine
convention (it may span both signs), `y1`/`y2` are the Q limits at that
operating point. The importer evaluates the curve once, at the machine's
scheduled SSH P (linear interpolation, clamped to the curve's P domain),
passes the pair through the same sign-convention hull as the scalar
`minQ`/`maxQ`, and stores the result as ordinary Q limits for the solver's
native Q-limit machinery (PV to PQ switching with hysteresis, cooldown and
guard).

Priority: curve, then scalar hull, then wide symmetric fallback
(so MicroGrid BE-G1, whose scalars are the degenerate `0/0` pair, gets its
±210 MVAr at P = −90 MW).

!!! details "Why the curve is not a voltage-dependent controller"
    The `QUController`/`PUController` path models droop injections with
    dQ/d|V| terms in the Jacobian (see
    [Voltage Dependent Control](voltage_dependent_control.md)), and a Q(P)
    bound folded into them would wander with voltage. A machine may carry
    both a Q(U) characteristic and curve-derived limits. Since P is fixed
    for a PV/PQ machine, the one-time evaluation is exact, with one
    simplification: under [distributed slack](powerflow_configuration.md)
    the λ_P correction shifts machine P, which strictly moves the curve
    limits; they stay at the SSH operating point.

### Out-of-service equipment (`Equipment.inService`)

CGMES 3.0 carries the operational status on the equipment itself in the
SSH profile; 2.4.15 only has `Terminal.connected`. The importer treats
`inService = false` as out of service for every mapped class: injections,
loads and shunts are skipped, branches go out of service, switches count
as open, each skip reported. ReliCapGrid parks whole plants and hundreds
of switches this way; importing them would add phantom generation and
merge separately solved islands. The `ignore_connected` override also
revives out-of-service equipment.

### HVDC

The importer never maps the DC side (`DCLineSegment`, `DCNode`); each
converter station is a fixed injection at its AC connection point, and
areas joined only through HVDC stay separate islands. Model ladder,
pairing controller and the reason there is no angle coupling:
[HVDC Back-to-Back](hvdc_back_to_back.md). Two delivery patterns:

- **Explicit converters** (`VsConverter`, `CsConverter`): fixed injections
  with the SSH operating point (p, q); the setpoint difference between the
  two stations of a link is the DC loss.
- **Ends without an SSH power**: an end that declares `targetUdc > 0` and
  carries `ACDCConverter.p = 0` gets its power derived from the opposite
  end, reduced by both converters' `idleLoss`, with a message. A classic
  LCC link that declares `CsConverter.targetIdc` (`pPccControl = dcCurrent`)
  on one end and `ACDCConverter.targetUdc` (`pPccControl = dcVoltage`) on
  the other gets their product as the DC power, `CsConverter.operatingMode`
  saying which end draws (`rectifier`) and which delivers (`inverter`); the
  DC line resistance is not applied. A pair with neither target is reported
  as not recoverable rather than given an invented number.
- **DC border crossings in assembled multi-area models**: a boundary node
  whose equivalent-injection pair does not cancel is such a crossing (a
  cancelling pair declares the same AC exchange twice and is discarded).
  The node is split per side, each area keeping its tie line and its
  equivalent on its own bus.

!!! warning "Precondition for the per-side split"
    Side identity comes from the defining file (`CIMObject.source`), the
    only criterion available since both sides reference the same boundary
    node. This works for deliveries with separate files per area
    (ReliCapGrid, the conformity assemblies); a CGM merged into a single
    file loses the criterion, the split silently does not apply, and the
    crossing stays galvanically joined.

!!! details "Why a zero SSH power at one end is not taken literally"
    In CIM the active power of a DC-voltage-controlling converter is a
    result of the DC power balance, not a setpoint, so an SSH snapshot
    legitimately carries `ACDCConverter.p = 0` at that end; taken
    literally the link would swallow its whole throughput.

## Tap control

With `tap_control = true` the importer starts from the SSH tap positions
and attaches outer-loop controllers for every tap changer whose
`controlEnabled` and `TapChangerControl.enabled` flags are set: voltage
control on ratio tap changers, active-power control on phase shifters,
also on three-winding transformers (each star-equivalent leg is an
ordinary PI-model branch). The run then goes through `run_sparlectra`
instead of a plain `runpf!`. A voltage controller whose target bus is held
by a generator (slack or PV) is disabled with a notice; CGMES target
deadbands are often wide, so a controller may settle one step away from
the recorded position.

## Machine remote voltage control

With `machine_control = true` a machine whose voltage `RegulatingControl`
terminal sits at a foreign bus is imported as a PQ injection with the SSH
operating point and gets a
[`MachineVoltageControl`](remote_voltage_control.md) that moves its
reactive output within the imported Q limits until the remote bus reaches
the target. By default such machines are held PV at their own bus, with a
notice. A plan also falls back to held-PV, with a notice in
`result.messages`, when the target bus is already voltage-held
(PV/slack), isolated, not part of the built network, or claimed by
another machine controller (one per target bus). The target value passes
through the `vset_min_pu`/`vset_max_pu` band against the remote bus's
nominal voltage.

## Validation against the SV profile

`compareWithSV(result)` compares the solved state with the delivery's SV
profile:

- **Voltages**: per-bus Δvm/Δva against `SvVoltage`, with max and RMS.
  Angles are only defined up to one constant per island (an IGM keeps the
  CGM's angle reference, the local solve pins its own slack), so the
  comparison removes the median offset before judging the angles
  (`dva_aligned` drives `max_dva`/`rms_dva`) and reports it separately as
  `va_ref_offset_deg` (cgmes.log, run metadata, Web UI summary row); the
  raw `dva` column stays in `sv_compare.csv`. Secondary islands with their
  own reference may keep a residual offset.
- **Flows** (`.flows`): per-terminal comparison against `SvPowerFlow` in
  the CGMES sign convention (branch terminals, shunts at the solved
  voltage, loads as an SSH-SV consistency check, units aggregated per
  bus). The ENTSO-E conformity sets only ship injection terminals.

De-energized and isolated buses are excluded. A sweep over every cached or
fetchable test set is `examples/run_cgmes_suite.jl`.

## Limitations

- CGMES 3.0 deliveries are read (`dcat:Dataset` headers, per-border
  boundary files, SSH `Equipment.inService`), validated against the
  ReliCapGrid/Svedala 3.0 sets only.
- Node-breaker deliveries are imported at bus-branch granularity (shipped
  TP, or [the topology processor](@ref topology_processor)); individual
  switches do not become model elements (retained ones become bus links).
- Per-step `r`/`x`/`g`/`b` corrections of tabular phase-tap tables are not
  folded into the branch impedance (the table's ratio and angle are); such
  rows are flagged in the import messages.
- Multi-valued references: `ref()` reads the first occurrence of a repeated
  property, `refsAll` the full list; the import emits one notice per
  affected class/property (`TopologicalIsland.TopologicalNodes` is the
  typical case). No mapped path consumes a list-valued reference.
- Difference models (`dm:DifferenceModel`) are skipped with a report line.

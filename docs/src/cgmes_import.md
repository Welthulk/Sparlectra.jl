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
are opened in memory. `importCGMES` also takes a store already loaded with
`loadCGMES` and leaves it unchanged. Profiles are classified from the
`md:Model` header. `result.element_by_mrid` maps the mRID of every imported
object to the element it became (kind and index).

## What is imported

| CGMES source | Result |
|---|---|
| `TopologicalNode` (TP) | buses, nominal voltage from `BaseVoltage`; without a TP profile the buses are derived from `ConnectivityNode`s and switch states ([topology processor](@ref topology_processor)) |
| `ACLineSegment`, `SeriesCompensator` | π-model branches; lines spanning two nominal voltages become ratio branches (boundary lines) |
| `EquivalentBranch` | branches without charging: on one voltage level a π branch with `r`, `x` on that base, across two levels a branch with per-unit ratio 1 and the impedance on the base of terminal 1; one with both terminals on one topological node builds nothing and gets a notice; one with an impedance below `cgmes_import.equivalent_branch_min_z_pu` joins its two buses as a bus link, with a notice; `r21`/`x21` other than `r`/`x` give a phase-stamped transformer branch (see [Asymmetric EquivalentBranch](@ref cgmes_eb_phase_stamp)); one notice counts mapped, skipped and stamped objects |
| `PowerTransformer` (2 and 3 windings) | π-model transformers; 3-winding as star equivalent with an AUX bus. A tap changer acts at the terminal of its own end, the series impedance lies on the inner side: a 3-winding leg with a tap changer runs from its bus to the star point, a 2-winding transformer with tap changers only on end 2 runs from end 2 to end 1 |
| `RatioTapChanger`, phase tap changers | fixed tap positions, or outer-loop controllers with `tap_control = true`; a `RatioTapChanger` with a `RatioTapChangerTable` takes the ratio of the table point at the tap position, and its controller steps from table point to table point; without a table the linear step formula applies; `PhaseTapChangerTabular` resolves its `PhaseTapChangerTable` row (ratio and angle) at the tap position. The `r`, `x`, `g`, `b` of a `RatioTapChangerTable` or `PhaseTapChangerTable` row are deviations in percent: the row of the step scales the series impedance and the magnetizing admittance (of a 3-winding leg, that leg's). With `tap_control = true` the import starts at the SSH step, and a controller that moves the tap sets the impedance of the step it reaches. The ENTSO-E `ExplicitLoadFlowCalculation` conformity set gives its phase tap table an `x` of 13.4036 at every step, the neutral step included; its SV does not apply that deviation, so the import differs from it by about 1.2 MW on the transformer. An end-2 (to-side) tap angle enters negated (`θ_eff = θ1 − θ2`), matching RealGrid's SV state; the ENTSO-E PSEI `PTE2` conformity sets expect the unflipped angle and deviate by about 0.3° |
| retained `Switch`/`Breaker`/… | zero-impedance bus links |
| `LinearShuntCompensator`, `NonlinearShuntCompensator` | shunts; the nonlinear characteristic sums its per-section points up to the active `sections` count (each point read as one switched group) |
| `EnergyConsumer`, `ConformLoad`, `NonConformLoad`, `StationSupply` | loads (SSH values) |
| `SynchronousMachine`, `ExternalNetworkInjection`, `EquivalentInjection` | injections, PV where a local voltage control is active; a `SynchronousMachine` or `ExternalNetworkInjection` with SSH `controlEnabled = false` does not regulate and is a PQ injection at its SSH `p`/`q` (absent counts as `true`), it stays a slack candidate by `referencePriority` or SV angle reference and as slack holds its bus at the start voltage (the SV value where present); an `EquivalentInjection` has no `controlEnabled` and regulates only with EQ `regulationCapability = true`, SSH `regulationStatus = true` and a positive `regulationTarget` (absent counts as `false`, one notice counts the capable units without a status); a machine whose voltage `RegulatingControl` points at a *different* bus is held PV at its own bus by default, or becomes an outer-loop remote voltage controller with `machine_control = true`; Q limits from the `ReactiveCapabilityCurve` (evaluated at the scheduled P) where one exists, else the scalar `minQ`/`maxQ` hull; the `minP`/`maxP` of an `ExternalNetworkInjection` or `EquivalentInjection` are read as delivered, in generator sign (a pair 0/0 or with `maxP` below `minP` gives none); a P limit of 9999 MW and above in magnitude (also `GeneratingUnit.minOperatingP`/`maxOperatingP`) is a placeholder and arrives as no limit, so it ranks no unit in the slack choice and weights none in the distributed slack (one notice counts these units per class; Q limits keep their delivered values); a positive `GeneratingUnit.normalPF` arrives as `ProSumer.participationFactor` for the [distributed slack](powerflow_configuration.md) (`p_mode: imported`), zero/absent maps to unknown |
| `StaticVarCompensator` | P = 0 reactive injection, Q limits from the Ω ratings; PQ at its SSH `q` with `controlEnabled = false` |
| `AsynchronousMachine` | fixed PQ operating point from the SSH `RotatingMachine.p`/`q` (load convention: a motor consumes with `p > 0`); no voltage regulation, no slack candidacy |
| `EnergySource` | only with `map_energy_source = true`: fixed PQ injection from the SSH `activePower`/`reactivePower`; no voltage regulation, a slack only in an island without other reference ([EnergySource](@ref cgmes_energy_source)) |
| `SvVoltage`, `SvPowerFlow`, `SvTapStep` | start values and validation reference |

Everything else is counted in the coverage report; `summarizeCGMES` shows
the full class histogram. Short-circuit source data is harvested on every
import (read, not evaluated); `shortCircuitCoverage` reports its
per-attribute completeness, printed in `cgmes.log`
([Short-Circuit Compendium](short_circuit.md)). A class with terminals
or with SSH active power of which no object reached the network (for
example `EnergySource` with `map_energy_source` off, or
`PowerElectronicsConnection`) gets one
`class not mapped:` row in the import messages and in `cgmes.log`, with its
object count and its SSH P sum; a class without SSH active power has no P
sum. One header line above the rows states the sign convention: load sign,
positive means flow out of the node.

Not counted as missing are the objects the topology consumes:
`BusbarSection`, `Junction`, and every switch (`Switch`, `Breaker`,
`Disconnector`, `LoadBreakSwitch`, `Jumper`, `ProtectedSwitch`, `Fuse`,
`GroundDisconnector`, `Cut`, `DisconnectingCircuitBreaker`) that is open,
disconnected or inside one TopologicalNode. A closed switch between two
nodes is a bus link for the first six classes; for the others the coupling
is missing and the class gets a row. Out-of-service or disconnected objects
of a class the importer maps are skipped on purpose and do not count.

## Conventions

| Convention | Sparlectra behaviour |
|---|---|
| Sign convention | CGMES uses the load convention (a machine with `p < 0` injects). Machine, external-injection, equivalent-injection and energy-source values are inverted; a load with negative `p` passes through unchanged. |
| `Terminal.connected`, both ends open | Branch out of service. |
| `Terminal.connected`, exactly one end open | No longitudinal current; the charging arm of the closed terminal stays as a shunt at the closed bus. |
| `Terminal.connected`, both ends closed | Normal branch. |
| `Terminal.connected` of a three-winding end | An open end leaves its leg open at the side bus while another end of the unit is connected; the magnetizing admittance of that end stays at the star point. With every end open or the unit out of service all three legs are out of service. |
| `PowerTransformerEnd.g`, `.b` | The magnetizing admittance of both ends of a two-winding transformer is summed and put half on each terminal of the branch, the from half behind the tap, as in the SV of the ENTSO-E test configurations (since 0.41.0). A three-winding end puts its admittance at the star point, the leg's arm at its bus is zero (since 0.42.0). The export writes equal arms back per end; a branch with unequal arms (from PowSyBl) writes the smaller arm on both ends and the rest as a `LinearShuntCompensator` at the bus ([Where the magnetizing admittance sits](@ref magnetizing_placement)). |
| Boundary sets | Cross-border deliveries reference topological nodes that live in the boundary files. Detected through unresolved references, not filenames; the import fails with an explicit message unless `require_boundary = false`. |
| Import analysis | On abort the importer prints the supplied model files (profile, version, model id), every `md:Model.DependentOn` prerequisite matched against the supplied models (a missing boundary set is named by its model id), an unresolved-reference histogram by class and property, and a verdict. `analyzeCGMES(path = ...)` does it on demand; the Web UI runs it on every upload ([Web UI Reference](webui_reference.md)). |
| Unresolved `TopologicalNode.BaseVoltage` | Cannot be repaired with `require_boundary = false`: real ENTSO-E deliveries keep their base-voltage catalog in the boundary EQ file (see `infer_base_voltages`). |
| Slack selection | In order: `referencePriority ≥ 1` (a notice on ties), a single `ExternalNetworkInjection`, the largest external injection (by `maxP`, else by `|p|`), the largest synchronous machine by `ratedS` (a rating of 9999 MVA and above is a placeholder and ranks as none; the rating also breaks ties). The `referencePriority` of every machine and network injection also stays on its unit (`ProSumer.referencePriority`) for island references and N-1 ([Reference priority](slack_vs_source.md#Reference-priority)). |
| Sourceless parts | Components without a source are de-energized (injections zeroed, branches out of service, buses isolated at voltage 0) with one message each and one notice with their number, their buses and the load left unsupplied, instead of aborting the power flow. A source is a generator-type unit with P > 0 (generator sign) or with voltage regulation, also of a remote bus in its island with `machine_control = true` (see Machine remote voltage control below); a unit without regulation that has P = 0 or draws power is none, so an island whose units only draw power is de-energized. A load is never a source, also not one that feeds power in. The unsupplied load includes the power the drawing units drew; the power that loads fed in is reported separately as dropped load infeed. |
| Islands without a reference | With `multi_slack = true` (default), an energized electrical island with infeed (a unit with P > 0 in generator sign) but without slack and without a voltage-regulating unit gets a slack from the importer. Slack candidates are only units with infeed, `EnergySource` included, which the slack choice above passes over; a unit with P = 0 or drawing power is never one, also not by `referencePriority`. Among them: the strongest `referencePriority`, else the largest active power reserve (`maxP` minus P) on a bus of degree greater than 1, else the largest infeed on a bus of degree greater than 1, and only when no such bus carries infeed, the largest infeed on a bus of degree 1 (both buses of a two-bus island have degree 1). Its voltage setpoint is the start (SV) voltage of its bus. Each choice is a `warning:` in the import messages, one notice counts them. Without a slack the island-wise solve would hold a promoted unit at its start voltage, and the result there would depend on the start value only. |
| Bus on a closed switch only | A bus without an in-service branch is marked isolated, also when a closed switch joins it to a bus with one. The solver still merges such a bus with that bus and solves it, and the SV comparison compares it at the voltage of its switch group. The import message counts these buses apart. |

## [Asymmetric EquivalentBranch](@id cgmes_eb_phase_stamp)

`EquivalentBranch.r21` and `x21` are the impedance seen from terminal 2.
With $z_{12} = r + jx$ and $z_{21} = r_{21} + jx_{21}$ the object is the
two-port

```math
I_1 = \frac{V_1 - V_2}{z_{12}}, \qquad I_2 = \frac{V_2 - V_1}{z_{21}},
```

the reading that fits the `SvPowerFlow` of real deliveries best. Its
admittance matrix is not symmetric. The importer builds such an object as a
transformer branch with ratio 1 and the phase stamp

```math
\varphi = \tfrac{1}{2}\arg\frac{z_{21}}{z_{12}}, \quad
c = \frac{1}{z_{12}}\sqrt{\frac{|z_{12}|}{|z_{21}|}}, \quad
z_s = \frac{e^{j\varphi}}{c}, \quad
y_{0f} = c\,(1 - e^{-j\varphi}), \quad
y_{0t} = c\,e^{-j\varphi}(e^{-j\varphi} - 1),
```

with the phase shift $\varphi$, the series impedance $z_s$ and the shunt
arms $y_{0f}$, $y_{0t}$, both impedances on the base of terminal 1 (across
two voltage levels with per-unit ratio 1). $|z_s|$ is the geometric mean
$\sqrt{|z_{12}||z_{21}|}$, and the arms make both row sums of the two-port
zero: $Y_{11} = -Y_{12} = c$, $Y_{22} = -Y_{21} = c\,e^{-2j\varphi}$. For
$|z_{21}| = |z_{12}|$ (a pure rotation) this is the two-port exactly. A
branch always has $|Y_{12}| = |Y_{21}|$, so otherwise each entry is off in
magnitude by $\sqrt{|z_{21}|/|z_{12}|}$ or its inverse. The stamped branch
gets no tap changer and no controller.

An `EquivalentBranch` is a network equivalent, not equipment: a negative
`r` or `x`, normal in Ward and Kron reductions, is taken over unchanged.
The same holds for a leg of a three-winding transformer in its star
equivalent: a single leg of a real unit can have a negative `r` or `x`,
only the sum of two legs is physical. The CGMES and MATPOWER exports and
the short circuit accept both; the check for a negative series resistance
applies to lines and two-winding transformers.

Limit: the error per branch is of the order
$|V_1 - V_2| \cdot \delta / |z|$ with
$\delta = \left||z_{21}|/|z_{12}| - 1\right|$, up to 1.7 % in the real
deliveries checked so far. It scales with the voltage difference, like the
branch flow. A stamp that keeps $Y_{11}$, $Y_{12}$, $Y_{22}$ exact puts the
whole error on $Y_{21}$, where it acts on the whole voltage $V_1$ and is of
the order of the short-circuit current; on a test delivery with
$\delta = 1\,\%$ that moved the solved voltages by 0.069 pu instead of
0.0002 pu. The import notice counts the stamped objects, the ones with
$|z_{21}| \ne |z_{12}|$ and the largest $\delta$. The model symmetrizes the magnitudes, but
the delivered `r`, `x`, `r21`, `x21` stay with the branch
(`Net.equivalent_branch_metadata`), and the export writes them unchanged
([CGMES Export](cgmes_export.md)).

## [Node-breaker deliveries without a TP profile](@id topology_processor)

Some EMS and substation-level exports ship EQ+SSH only, topology as
`ConnectivityNode`s plus switches. The importer then derives the bus
partition itself.

- **Trigger.** Only when the delivery contains `ConnectivityNode`s and no
  non-boundary `TopologicalNode`; a TP-carrying delivery is unaffected.
- **Aggregation.** Connectivity nodes merge across closed, non-retained
  switches (`Switch`, `Breaker`, `Disconnector`, `LoadBreakSwitch`,
  `Jumper`, `ProtectedSwitch`); a switch counts as open when SSH
  `Switch.open` says so (overriding EQ `normalOpen`) or when it is out of
  service. Retained switches become zero-impedance bus links.
- **Result.** Each connectivity group becomes a synthetic
  `TopologicalNode`, named after a busbar section in the group where one
  exists; the nominal voltage resolves through `VoltageLevel` (hopping
  over `Bay` containers) to the `BaseVoltage`. Boundary connectivity nodes
  adopt the boundary set's TP_BD nodes, so cross-border stitching works as
  with a shipped TP.
- **Visibility.** The import messages announce the processor
  (`topology processor: derived N topological node(s) from M connectivity
  node(s) ...`); without that message the shipped TP was used. A delivery
  with neither a TP profile nor connectivity nodes aborts with the import
  analysis.

On MiniGrid and SmallGrid the derived partition reproduces the shipped TP
bus for bus. FullGrid's TP assigns the terminals of one connectivity node
(the shared load node of `BE-Load_1`/`BE_CL_1`/`BE_NC_1`) to two
topological nodes, which no processor can derive from the graph, so 10 MW
of load sit one bus apart; dead network fragments may be partitioned
differently than the sender's TP and show up as extra isolated,
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
| `cgmes_import.ignore_connected` | `false` | Diagnostic override treating every terminal as connected, for snapshots whose SSH flags contradict their own SV state; also revives out-of-service equipment. Allowed values: `true`, `false`. |
| `cgmes_import.map_energy_source` | `false` | Map `EnergySource` objects as fixed PQ injections ([EnergySource](@ref cgmes_energy_source)). Allowed values: `true`, `false`. |
| `cgmes_import.vset_min_pu` | `0.5` | Lower bound of the plausibility band for a voltage `RegulatingControl.targetValue`, in p.u. of the regulated bus's nominal voltage. |
| `cgmes_import.vset_max_pu` | `1.5` | Upper bound of that band. A target outside `[vset_min_pu, vset_max_pu]` is treated as a placeholder: it is ignored, the unit is held PV at the bus voltage derived from the nominal data, and the substitution is reported as a `warning:`. |
| `cgmes_import.multi_slack` | `true` | Give every electrical island its own SV-declared angle reference (at most one per island), and a slack to every island with infeed but without reference (see Islands without a reference above). Required for multi-island deliveries; `false` forces the single-reference behavior. Allowed values: `true`, `false`. |
| `cgmes_import.equivalent_branch_min_z_pu` | `1e-6` | An `EquivalentBranch` whose impedance in both directions (`r`, `x` and `r21`, `x21`, per unit on the base of terminal 1) lies below this value builds no branch: its two buses are joined as a bus link, with a notice. `0` switches it off. |
| `cgmes_import.placeholder_guards` | `warn_skip` | Behavior of the placeholder guards (implausible shunt admittances, tap corrections outside 0.5 … 2.0). `warn_skip`: keep the filler value out of the solve with a warning. `strict`: abort the import with an error naming the object. Allowed values: `warn_skip`, `strict`. |
| `cgmes_import.infer_base_voltages` | `false` | Reconstruct missing nominal voltages when the delivery ships without its `BaseVoltage` catalog: nodes are seeded from the SV voltages (kV, snapped to the standard level series) and transformer rated voltages, then the levels propagate across level-preserving equipment (anything but a transformer). All substitutions are summarized as one `warning:` with per-source and per-level counts. Pair with `require_boundary: false`; nodes that stay unresolved still abort with the import analysis. Allowed values: `true`, `false`. |
| `cgmes_import.hvdc_mode` | `injections` | HVDC converter handling. `injections`: fixed PCC injections with the SSH operating point; areas joined only through HVDC stay separate islands. `paired_control`: the same injections plus one steerable `HvdcPairControl` per detected converter pair; the DC topology (`ACDCConverterDCTerminal`, `DCNode`, `DCLineSegment`) groups the converters into links, the transfer and loss are derived from the two SSH operating points, and the outer loop keeps the pairing exact. Detection runs in both modes and names every pair in the import messages; a pair that cannot be attached (skipped converter, inconsistent snapshot) degrades to the fixed injections with a notice. Allowed values: `injections`, `paired_control`. |
| `cgmes_import.start_values` | `auto` | Newton-Raphson start state for CGMES runs. `auto`: the delivery's `SvVoltage` state when it carries one, else the flat start. `sv`: always the imported state. `flat`: always the flat start. With an SV start the other start-value machines (`start_projection`, `dc_seed_unconditional`, `start_current_iteration`, `apslf_start`) are forced off for the run. An explicit `sv` or `flat` wins over `power_flow.flatstart`; under `auto` a set `power_flow.flatstart` chooses the flat start although the delivery carries an SV state. MATPOWER and DTF runs ignore this key. The resolved decision (what `auto` chose, overridden keys) is logged to `run.log` and `cgmes.log`; the SV comparison (`sv_compare.csv`) runs in every mode. Allowed values: `auto`, `flat`, `sv`. |

### Implausible voltage setpoints

A voltage `RegulatingControl.targetValue` outside the
`vset_min_pu`/`vset_max_pu` band is treated as a placeholder (see the
table); some deliveries park units with targets such as `0.001` kV. Real
targets in the tested deliveries span 0.92 … 1.15 p.u. Widen the band if
a delivery legitimately regulates outside it, or set `vset_min_pu: 0.0`
with a large `vset_max_pu` to accept every value.

### Placeholder guards: shunt admittances and tap corrections

Two guards catch filler values such as FullGrid's `X.99` attributes;
`cgmes_import.placeholder_guards: strict` aborts with an error naming the
object instead of skipping.

- a shunt whose admittance exceeds **10 × baseMVA** at nominal voltage is
  skipped (not clamped) with a `warning:`;
- a single tap correction factor outside **0.5 … 2.0** is ignored with a
  `warning:`; the transformer keeps its nominal ratio.

With the guards FullGrid solves from a flat start; its shipped SV profile
is internally inconsistent, so the SV start and comparison are meaningless
there.

### Machine Q limits: the `ReactiveCapabilityCurve` is Q(P), not Q(U)

A `ReactiveCapabilityCurve` gives reactive limits as a function of active
power: the `CurveData` x axis is the machine's own P in the CGMES machine
convention (it may span both signs), `y1`/`y2` are the Q limits there.
The importer evaluates the curve once, at the scheduled SSH P (linear
interpolation, clamped to the curve's P domain), passes the pair through
the same sign-convention hull as the scalar `minQ`/`maxQ`, and stores the
result as ordinary Q limits. Priority: curve, then scalar hull, then a
wide symmetric fallback when the scalars are the degenerate `0/0` pair.
An `EquivalentInjection` follows the same order without the fallback: it
has no rating to derive one from, so without a usable pair it keeps no Q
limits. A
machine may carry both a
Q(U) characteristic ([Voltage Dependent Control](voltage_dependent_control.md))
and curve-derived limits. Under [distributed slack](powerflow_configuration.md)
the λ_P correction shifts machine P; the curve limits stay at the SSH
operating point.

### [EnergySource](@id cgmes_energy_source)

Some deliveries carry generation as `EnergySource` instead of
`SynchronousMachine`. With `cgmes_import.map_energy_source: true` every
`EnergySource` that is in service, has a connected terminal and sits on a
TopologicalNode becomes a fixed PQ injection with
`p = -activePower` and `q = -reactivePower` from the SSH. CGMES states both
in load sign convention (positive means flow out of the node into the
source), like `SvPowerFlow` at the source's terminal, so a consistent
delivery reproduces `p = -SvPowerFlow.p`. The injection has no voltage
regulation and no Q limits, and the importer's slack choice passes it over. It
is addressed as `:generator` in `element_by_mrid`, and `compareWithSV`
counts it with the units of its bus.

A mapped `EnergySource` with P > 0 (injecting) is a source, so an island
it feeds stays energized (with the key off the same island is de-energized
as sourceless). In an island whose only infeed is `EnergySource` the importer
makes the largest one the slack, by the rule for islands without a
reference, with a warning; that source's P then follows the island balance
instead of the SSH value. With `multi_slack: false` such an island has no
slack after the import: the island-wise solve (`power_flow.islands`) makes
its strongest `EnergySource` the island reference, with
`power_flow.islands.enabled: false` the solve rejects the island, and the
import names each such island in a notice.

The key is off by default, because it changes the imported network of
every delivery that carries such objects. Off, they stay out of the
network and are listed in `skipped_equipment`; one notice in the import
messages (and in `cgmes.log`) names the key that maps them, and their
count and P sum stand in the table of classes not mapped. On, one notice
counts the mapped and skipped
objects with the P and Q sums. A delivery whose generation sits in
`EnergySource` objects reproduces its SV state only with the key on.

### Out-of-service equipment (`Equipment.inService`)

CGMES 3.0 carries the operational status on the equipment itself in the
SSH profile; 2.4.15 only has `Terminal.connected`. The importer treats
`inService = false` as out of service for every mapped class: injections,
loads and shunts are skipped, branches go out of service, switches count
as open, each skip reported. The `ignore_connected` override revives such
equipment.

### HVDC

The importer never maps the DC side (`DCLineSegment`, `DCNode`); each
converter station is a fixed injection at its AC connection point, and
areas joined only through HVDC stay separate islands. Model ladder,
pairing controller and the reason there is no angle coupling:
[HVDC Back-to-Back](hvdc_back_to_back.md). Three delivery patterns:

- **Explicit converters** (`VsConverter`, `CsConverter`): fixed injections
  with the SSH operating point (p, q); the setpoint difference between the
  two stations of a link is the DC loss.
- **Ends without an SSH power**: a DC-voltage-controlling end may carry
  `ACDCConverter.p = 0` in the SSH (its power is a result in CIM, not a
  setpoint); with `targetUdc > 0` it gets its power from the opposite
  end, reduced by both converters' `idleLoss`, with a message. A classic
  LCC link with
  `CsConverter.targetIdc` (`pPccControl = dcCurrent`) on one end and
  `ACDCConverter.targetUdc` (`pPccControl = dcVoltage`) on the other gets
  their product as the DC power, `CsConverter.operatingMode` saying which
  end draws (`rectifier`) and which delivers (`inverter`); the DC line
  resistance is not applied. A pair with neither target is reported as
  not recoverable.
- **DC border crossings in assembled multi-area models**: a boundary node
  whose equivalent-injection pair does not cancel is such a crossing (a
  cancelling pair declares the same AC exchange twice and is discarded);
  the node is split per side, each area keeping its tie line and its
  equivalent on its own bus. Side identity comes from the defining file
  (`CIMObject.source`), so the split needs separate files per area; a
  CGM merged into a single file loses the criterion, the split silently
  does not apply, and the crossing stays galvanically joined.

## Tap control

With `tap_control = true` the importer starts from the SSH tap positions
and attaches outer-loop controllers for every tap changer whose
`controlEnabled` and `TapChangerControl.enabled` flags are set: voltage
control on ratio tap changers, active-power control on phase shifters,
also on three-winding transformers. The run then goes through
`run_sparlectra` instead of a plain `runpf!`. A voltage controller whose
target bus is held by a generator (slack or PV) is disabled with a
notice; CGMES target deadbands are often wide, so a controller may settle
one step away from the recorded position.

Every tap changer, regulated or not and with or without `tap_control`,
gives its branch the tap range, the ratios of its `RatioTapChangerTable`
and the impedance per step. A controller added later through the
`control` section then steps on the delivered table; the solve itself
uses the imported tap position.

## Machine remote voltage control

With `machine_control = true` a machine whose voltage `RegulatingControl`
terminal sits at a foreign bus is imported as a PQ injection with the SSH
operating point and gets a
[`MachineVoltageControl`](remote_voltage_control.md) that moves its
reactive output within the imported Q limits until the remote bus reaches
the target. A machine with `controlEnabled = false` gets no controller
and stays PQ. By default, and when the target bus is already voltage-held
(PV/slack), isolated, not part of the built network or claimed by another
machine controller (one per target bus), the machine is held PV at its
own bus with a notice in `result.messages`. The same holds when its
regulation is the only source of its island (no slack, no PV bus, no unit
with P > 0, for example a pump): only units with infeed are slack
candidates, so the island gets no slack from the importer, and the machine
is held PV at its own bus at its start voltage. The target value passes
through the `vset_min_pu`/`vset_max_pu` band against the remote bus's
nominal voltage.

## Validation against the SV profile

`compareWithSV(result)` compares the solved state with the delivery's SV
profile:

- **Voltages**: per-bus Δvm/Δva against `SvVoltage`, with max and RMS.
  Angles are only defined up to one constant per island, so the
  comparison removes the median offset per island (`dva_aligned` drives
  `max_dva`/`rms_dva`) and reports it as `va_ref_offset_deg` (cgmes.log,
  run metadata, Web UI summary row); the raw `dva` column stays in
  `sv_compare.csv`.
- **Flows** (`.flows`): per-terminal comparison against `SvPowerFlow` in
  the CGMES sign convention (branch terminals, shunts at the solved
  voltage, loads as an SSH-SV consistency check, units aggregated per
  group of buses joined by closed switches, since a switch carries power
  the balance of a single bus does not see). The ENTSO-E conformity sets
  only ship injection terminals.
- **Tap steps**: a two-winding transformer whose `SvPowerFlow` at both
  terminals fits another tap step clearly better than its `SvTapStep`
  (mismatch at most a tenth and at least 1 MVA smaller, at the SV
  voltages) was solved with a step the delivery does not state. The
  import counts such tap changers in one notice in `result.messages` and
  lists them in `result.sv_tap_step_fits` with both steps and mismatches;
  the fixed-reference self-check lists them in `self_check.log`.

De-energized and isolated buses are excluded, except a bus that a closed
switch joins to a bus with a branch. `examples/run_cgmes_suite.jl`
sweeps every cached or fetchable test set.

## Limitations

- CGMES 3.0 deliveries are read (`dcat:Dataset` headers, per-border
  boundary files, SSH `Equipment.inService`), validated against the
  ReliCapGrid/Svedala 3.0 sets only.
- Node-breaker deliveries are imported at bus-branch granularity (shipped
  TP, or [the topology processor](@ref topology_processor)); switches do
  not become model elements (retained ones become bus links).
- Multi-valued references: `ref()` reads the first occurrence of a repeated
  property, `refsAll` the full list; the import emits one notice per
  affected class/property (typically `TopologicalIsland.TopologicalNodes`).
  No mapped path consumes a list-valued reference.
- An `EquivalentBranch` with $|z_{21}| \ne |z_{12}|$ is not reproduced
  exactly: the branch model has $|Y_{12}| = |Y_{21}|$, so it is built with
  the geometric mean of the two magnitudes ([Asymmetric EquivalentBranch](@ref cgmes_eb_phase_stamp)).
- Difference models (`dm:DifferenceModel`) are skipped with a report line.

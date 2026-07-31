# CGMES Import

Sparlectra reads ENTSO-E **CGMES 2.4.15** deliveries (Common Grid Model
Exchange Standard) and builds a bus-branch `Net` that solves with `runpf!`.
The importer is deliberately lean: it reads the profiles it needs, reports
everything it skips, and validates the result against the solved state the
delivery ships with.

## Quick start

```julia
using Sparlectra

# 1. Diagnose first — works on incomplete or broken deliveries
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

`path` accepts a folder, a single ZIP, or a vector of either — nested ZIPs are
opened in memory. Profiles are classified from the `md:Model` header, never
from filenames.

## What is imported

| CGMES source | Result |
|---|---|
| `TopologicalNode` (TP) | buses, nominal voltage from `BaseVoltage` |
| `ACLineSegment`, `SeriesCompensator` | π-model branches; lines spanning two nominal voltages become ratio branches (boundary lines) |
| `PowerTransformer` (2 and 3 windings) | π-model transformers; 3-winding as star equivalent with an AUX bus |
| `RatioTapChanger`, phase tap changers | fixed tap positions, or outer-loop controllers with `tap_control = true`; `PhaseTapChangerTabular` resolves its `PhaseTapChangerTable` row (ratio and angle) at the tap position. Tap angles fold with end-referral semantics: an end-2 (to-side) angle enters negated (`θ_eff = θ1 − θ2`), pinned by RealGrid's SV state — the ENTSO-E PSEI `PTE2` conformity toys expect the unflipped angle and deviate by ≈0.3° by design of this choice |
| retained `Switch`/`Breaker`/… | zero-impedance bus links |
| `LinearShuntCompensator`, `NonlinearShuntCompensator` | shunts; the nonlinear characteristic sums its per-section points up to the active `sections` count (each point read as one switched group — an interpretation choice, CIM would also permit reading a point as the absolute value at that section count) |
| `EnergyConsumer`, `ConformLoad`, `NonConformLoad`, `StationSupply` | loads (SSH values) |
| `SynchronousMachine`, `ExternalNetworkInjection`, `EquivalentInjection` | injections, PV where a local voltage control is active; a machine whose voltage `RegulatingControl` points at a *different* bus is held PV at its own bus by default, or becomes an outer-loop remote voltage controller with `machine_control = true` (see [Remote Voltage Control](remote_voltage_control.md)); machine Q limits from the `ReactiveCapabilityCurve` (evaluated at the scheduled P) where one exists, else the scalar `minQ`/`maxQ` hull; a positive `GeneratingUnit.normalPF` arrives as `ProSumer.participationFactor` for the [distributed slack](powerflow_configuration.md) (`p_mode: imported`), zero/absent maps to unknown |
| `StaticVarCompensator` | P = 0 reactive injection, Q limits from the Ω ratings |
| `AsynchronousMachine` | Fixed PQ operating point from the SSH `RotatingMachine.p`/`q` (load convention — a motor consumes with `p > 0`); no voltage regulation, no slack candidacy. The induction machine's voltage/slip dependence is a dynamics topic, not a steady-state one |
| `SvVoltage`, `SvPowerFlow`, `SvTapStep` | start values and validation reference |

Everything else is counted in the coverage report rather than silently
dropped. `summarizeCGMES` shows the full class histogram.

Short-circuit source data is harvested on every import (read, not evaluated);
`shortCircuitCoverage` reports its per-attribute completeness, and `cgmes.log`
prints that view. Theory, data table and the staged evaluation concept:
[Short-Circuit Compendium](short_circuit.md).

## Conventions worth knowing

**Sign conventions.** CGMES uses the load convention: a machine with `p < 0`
injects. The importer inverts machine, external-injection and
equivalent-injection values; a load with negative `p` is passed through
unchanged (it genuinely is an injection).

**`Terminal.connected`.** Three cases, following real-snapshot semantics:

- both ends open → branch out of service,
- **exactly one end open** → the branch carries no longitudinal current, but
  its half charging admittance stays as a shunt at the closed bus (dropping
  it entirely would distort the reactive balance),
- both ends closed → normal branch.

**Boundary sets.** Cross-border deliveries reference topological nodes that
live in the boundary files. The importer detects this through actually
unresolved references, not through filenames, and fails with an explicit
message unless `require_boundary = false`.

**Import analysis.** When such an import aborts, the importer first prints a
full analysis: the supplied model files (profile, version, model id), every
`md:Model.DependentOn` prerequisite the file headers declare — matched
against the supplied models, so a missing boundary set is named by its exact
model id — an unresolved-reference histogram by class and property, and a
plain-language verdict. In Web UI runs the report lands in `cgmes.log`, and
the **Analyze import** button runs the same check before a full import (see
[Web UI](webui.md)). The report is also available on demand without
importing: `analyzeCGMES(path = ...)` accepts the same path forms as
`importCGMES`. Note that a delivery whose `TopologicalNode.BaseVoltage`
references stay unresolved cannot be repaired with `require_boundary =
false`: real ENTSO-E deliveries keep their base-voltage catalog in the
boundary EQ file, so without the matching boundary there is no voltage level
to build buses from.

**Slack selection.** In order: `referencePriority ≥ 1` (a notice is emitted
when several units tie), a single `ExternalNetworkInjection`, the largest
external injection, the largest synchronous machine.

**Sourceless parts.** Components without any generator are de-energized
deliberately (injections zeroed, buses isolated) with one message each,
instead of aborting the power flow with "island without reference".

## Configuration (`cgmes_import`)

| Key | Default | Purpose |
|---|---|---|
| `cgmes_import.path` | `""` | Delivery location(s); `;`-separated for multi-part deliveries (base case plus boundary set). |
| `cgmes_import.base_mva` | `100.0` | System base in MVA — CGMES does not define one. |
| `cgmes_import.require_boundary` | `true` | Fail when topology references stay unresolved. Allowed values: `true`, `false`. |
| `cgmes_import.tap_control` | `false` | Start from the SSH tap positions and attach the CGMES-defined outer-loop tap controllers, instead of importing the solved `SvTapStep` positions as fixed taps. Allowed values: `true`, `false`. |
| `cgmes_import.machine_control` | `false` | Attach outer-loop remote voltage controllers (`MachineVoltageControl`) for machines whose voltage `RegulatingControl` points at a different bus, instead of holding those machines PV at their own bus. Allowed values: `true`, `false`. |
| `cgmes_import.ignore_connected` | `false` | Diagnostic override treating every terminal as connected — for snapshots whose SSH flags contradict their own SV state. Allowed values: `true`, `false`. |
| `cgmes_import.vset_min_pu` | `0.5` | Lower bound of the plausibility band for a voltage `RegulatingControl.targetValue`, in p.u. of the regulated bus's nominal voltage. |
| `cgmes_import.vset_max_pu` | `1.5` | Upper bound of that band. A target outside `[vset_min_pu, vset_max_pu]` is treated as a placeholder: it is ignored, the unit is held PV at the bus voltage derived from the nominal data, and the substitution is reported as a `warning:`. |
| `cgmes_import.multi_slack` | `true` | Give every electrical island its own SV-declared angle reference (at most one per island). Required for multi-island deliveries; `false` forces the legacy single-reference behavior. Allowed values: `true`, `false`. |
| `cgmes_import.placeholder_guards` | `warn_skip` | Behavior of the placeholder guards (implausible shunt admittances, tap corrections outside 0.5 … 2.0). `warn_skip`: keep the filler value out of the solve with a warning. `strict`: abort the import with an error naming the object — for deliveries where dropped data must never go unnoticed. Allowed values: `warn_skip`, `strict`. |
| `cgmes_import.infer_base_voltages` | `false` | Reconstruct missing nominal voltages when the delivery ships without its `BaseVoltage` catalog (in real ENTSO-E deliveries the catalog lives in the boundary EQ): nodes are seeded from the SV voltages (kV, snapped to the standard level series) and transformer rated voltages, then the levels propagate across level-preserving equipment (anything but a transformer). All substitutions are summarized as one `warning:` message with per-source and per-level counts. Pair with `require_boundary: false`; nodes that stay unresolved still abort with the import analysis. Allowed values: `true`, `false`. |
| `cgmes_import.start_values` | `flat` | Newton-Raphson start state for CGMES runs. `flat`: synthetic flat start — the solver earns the solution itself. `sv`: start from the delivery's imported `SvVoltage` state; the competing start-value machines (`start_projection`, `dc_seed_unconditional`, `start_current_iteration`, `apslf_start`) are forced off for the run. On CGMES runs this key wins over `power_flow.flatstart` / `power_flow.start_mode.flatstart`; MATPOWER and DTF runs ignore it. The decision (including any overridden keys) is logged to `run.log` and `cgmes.log`, and the SV comparison (`sv_compare.csv`) runs in both modes. Allowed values: `flat`, `sv`. |

### Implausible voltage setpoints

Some deliveries carry placeholder regulation targets. ReliCapGrid's Svedala
model declares `targetValue = 0.001` kV on 17 kV and 20 kV generator busbars —
about `5e-5` p.u. — for exactly the five units that also carry no `SvVoltage`,
i.e. machines the exporting tool left out of its solved state. Taken literally
such a target turns the reference bus of its island into a bus at essentially
zero volts, and the power flow then converges onto the matching all-zero
solution: formally correct, physically meaningless.

The band exists to catch that. Across every ENTSO-E and ReliCapGrid delivery
tested, the real targets span 0.92 … 1.15 p.u., so the default leaves ample
margin. It is configurable because the number comes from observed data, not from
the standard: widen it if a delivery legitimately regulates outside the band, or
set `vset_min_pu: 0.0` with a large `vset_max_pu` to accept every value.

Note that the ReliCapGrid units in question also carry SSH
`Equipment.inService = false`, which the importer honors (see below) — they are
skipped before their setpoint is ever read. The band remains as a guard for
deliveries that park units without setting `inService`.

### Placeholder guards: shunt admittances and tap corrections

The same philosophy covers two more fields where placeholders otherwise
destroy the solve. FullGrid — the completeness configuration — systematically
fills attributes with the `X.99` scheme (tabular PST table row
`ratio 9.99 / angle 0.99°`, a `NonlinearShuntCompensatorPoint` with
`b = g = 0.99 S`, which at 225 kV is a 50-GW shunt, switch
`ratedCurrent 999.99`, …). Two guards catch these:

- a shunt whose admittance exceeds **10 × baseMVA** at nominal voltage is
  skipped with a `warning:` (skipped, not clamped — a placeholder carries no
  information to clamp to);
- a single tap correction factor outside **0.5 … 2.0** is ignored with a
  `warning:`; the transformer keeps its nominal ratio. Real tap ranges stay
  within a few ten percent of neutral (RealGrid tabular rows: 0.9 … 1.1).

With the guards in place FullGrid's network solves from a flat start (its
shipped SV profile remains internally inconsistent — a 14.5° angle jump
across a 0.3 Ω line — so the SV-based start and the SV comparison stay
meaningless for this set).

Both guards act globally with warn-and-skip semantics by default. When
silently dropping data is not acceptable — a productive delivery rather
than a conformity set — set `cgmes_import.placeholder_guards: strict`: a
suspected placeholder then aborts the import with an error naming the
offending object instead of skipping it.

### Machine Q limits: the `ReactiveCapabilityCurve` is Q(P), not Q(U)

A `ReactiveCapabilityCurve` is the machine's operating chart: **reactive
limits as a function of active power** — the `CurveData` x axis is the
machine's own P in the CGMES machine convention (it may legitimately span
both signs; MicroGrid BE-G1: −100 … +100 MW), `y1`/`y2` are the Q limits at
that operating point. It is *not* a voltage dependence.

The importer therefore evaluates the curve **once, at the machine's
scheduled SSH P** (linear interpolation, clamped to the curve's P domain),
passes the pair through the same sign-convention hull as the scalar
`minQ`/`maxQ`, and stores the result as the machine's ordinary Q limits.
From there the limits feed the rectangular solver's **native Q-limit
machinery** — the PV→PQ active-set switching with hysteresis, cooldown and
guard — exactly like scalar limits do. Priority chain: curve where one
exists → scalar hull → wide symmetric fallback (so a machine like BE-G1,
whose scalars are the degenerate `0/0` pair precisely because its real
limits live in the curve, gets its ±210 MVAr at P = −90 MW).

Deliberately **not** used for this: the `QUController`/`PUController`
voltage-dependent control path. Those model droop *injections* — a setpoint
as a function of the local bus voltage magnitude, re-evaluated every
iteration with dQ/d|V| terms in the Jacobian (see
[Voltage Dependent Control](voltage_dependent_control.md)). Folding a Q(P)
capability *bound* into that machinery would make the limits wander with
voltage and turn bounds into setpoints — both data-unfaithful. The two
mechanisms stay orthogonal: a machine may carry a Q(U) characteristic *and*
curve-derived limits.

Since P is fixed for a PV/PQ machine in the power flow, the one-time
evaluation is exact, with one documented simplification: under
[distributed slack](powerflow_configuration.md) the λ_P correction shifts
machine P, which strictly moves the curve limits; they currently stay at
the SSH operating point (negligible at measured corrections — RealGrid:
λ_P ≈ 5 MW over 368 participants, and zero curves in that delivery).

### Out-of-service equipment (`Equipment.inService`)

CGMES 3.0 carries the operational status on the equipment itself in the SSH
profile; 2.4.15 only has `Terminal.connected`. The importer treats
`inService = false` as out of service for every mapped class: injections,
loads and shunts are skipped, branches go out of service, switches count as
open. This matters — ReliCapGrid parks whole plants and hundreds of switches
with `connected = true, inService = false` and no SV state; importing them adds
phantom generation and falsely merges network parts the delivery solved as
separate islands. Each skip is reported. The `ignore_connected` diagnostic
override also revives out-of-service equipment, consistent with its "treat
everything as live" meaning.

### HVDC

The one property of HVDC that matters for a load flow is that it has **no
angle coupling**: the transfer is a control setpoint, not the result of an
angle difference. Two areas joined only through HVDC are electrically separate
islands with their own angle references. The importer therefore never maps the
DC side (`DCLineSegment`, `DCNode`) and represents each converter station as a
fixed injection at its AC connection point — the industry-standard treatment
(MATPOWER's `dcline` does the same with two bounded dummy generators). Two
delivery patterns are handled:

- **Explicit converters** (`VsConverter`, `CsConverter`): mapped as fixed
  injections with the SSH operating point (p, q). The setpoint difference
  between the two stations of a link is the DC loss.
- **DC border crossings in assembled multi-area models**: a boundary node
  whose equivalent-injection pair does *not* cancel is such a crossing (a
  cancelling pair is two declarations of the same AC exchange and is
  discarded; a non-cancelling pair is the two converter injections). The node
  is split per side — each area keeps its tie line and its equivalent on its
  own bus — instead of galvanically short-circuiting a coupling that does not
  exist. The side of a piece of equipment is the file that defined it.

!!! warning "Precondition for the per-side split"
    Side identity comes from the *defining file* (`CIMObject.source`) — the
    only criterion available, since both sides reference the same boundary
    node. This works for deliveries with separate files per area (ReliCapGrid,
    the conformity assemblies). A CGM merged into a **single file** loses the
    criterion, and the split silently does not apply: the crossing then stays
    galvanically joined and the convergence/SV symptoms described above
    return. This is a known limitation, not a bug.

## Tap control

With `tap_control = true` the importer starts from the SSH tap positions and
attaches Sparlectra's outer-loop controllers for every tap changer whose
`controlEnabled` and `TapChangerControl.enabled` flags are set — voltage
control on ratio tap changers, active-power control on phase shifters. This
includes tap changers on three-winding transformers: each star-equivalent leg
is an ordinary PI-model branch (AUX bus → side bus), so the leg carries the
controller like a two-winding transformer would. The
run then goes through `run_sparlectra` (control framework) instead of a plain
`runpf!`.

Two guards apply: a voltage controller whose target bus is held by a generator
(slack or PV) cannot regulate anything and is disabled with a notice; and the
CGMES target deadbands are often wide, so a controller may legitimately settle
one step away from the position the exporting tool recorded.

## Machine remote voltage control

With `machine_control = true` a machine whose voltage `RegulatingControl`
terminal sits at a *foreign* bus is imported as a PQ injection with the SSH
operating point and gets a [`MachineVoltageControl`](remote_voltage_control.md)
attached: the outer control loop moves the machine's reactive output within
its imported Q limits until the remote bus reaches the control target. Without
the option (the default) such machines keep the Stage-1 behavior — held PV at
their own bus, with a notice.

A plan falls back to held-PV, each time with a notice in `result.messages`,
when the target bus is already voltage-held (PV/slack), isolated, not part of
the built network, or already claimed by another machine controller (one
controller per target bus; further machines keep their SSH reactive output).
The target value passes through the same `vset_min_pu`/`vset_max_pu`
plausibility band as local setpoints, evaluated against the *remote* bus's
nominal voltage.

## Validation against the SV profile

`compareWithSV(result)` compares the solved state with the delivery's own
State Variables profile:

- **Voltages** — per-bus Δvm/Δva against `SvVoltage`, with max and RMS.
  Angles are only defined up to one constant per island: an IGM cut out of
  the continental CGM keeps the CGM's global angle reference, while the
  local solve pins its own slack — a uniform offset of tens of degrees that
  says nothing about the state. The comparison removes the median offset
  before judging the angles (`dva_aligned` drives `max_dva`/`rms_dva`) and
  reports it separately as `va_ref_offset_deg` (cgmes.log, run metadata,
  Web UI summary row); the raw `dva` column stays in `sv_compare.csv`.
  Secondary islands with their own reference may keep a residual offset.
- **Flows** (`.flows`) — per-terminal comparison against `SvPowerFlow` in the
  CGMES sign convention: branch terminals, shunts at the solved voltage, loads
  as an SSH↔SV consistency check, and units aggregated per bus. Note that the
  ENTSO-E conformity data sets only ship injection terminals; branch flows
  appear in real exchanges.

De-energized and isolated buses are excluded so the metrics describe the
solved grid.

A full sweep over every cached/fetchable test set — import, solve, SV
deviation, one table row per case incl. RealGrid and the ReliCapGrid/Svedala
3.0 family — is part of `examples/run_cgmes_suite.jl`; the measured state is
recorded in `docs/dev/cgmes_testset_overview.md`.

## Limitations

- CGMES 3.0 deliveries are read (`dcat:Dataset` headers, per-border boundary
  files, SSH `Equipment.inService`), validated against the ReliCapGrid/Svedala
  3.0 sets; conformity coverage beyond those deliveries is still limited.
- Node-breaker modelling is not implemented; it is not needed in practice
  because CGMES deliveries ship the TP profile alongside, which the importer
  reads as bus-branch.
- Per-step `r`/`x`/`g`/`b` corrections of tabular phase-tap tables are not
  folded into the branch impedance (the table's ratio and angle are); rows
  with such corrections are flagged in the import messages.
- Multi-valued references: `ref()` reads the first occurrence of a repeated
  property; the full document-order list is available via `refsAll`, and the
  import emits one notice per affected class/property
  (`TopologicalIsland.TopologicalNodes` is the typical case). No mapped path
  consumes a list-valued reference yet.
- Difference models (`dm:DifferenceModel`) are skipped with a report line.

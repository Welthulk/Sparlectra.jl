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
| `RatioTapChanger`, phase tap changers | fixed tap positions, or outer-loop controllers with `tap_control = true` |
| retained `Switch`/`Breaker`/… | zero-impedance bus links |
| `LinearShuntCompensator` | shunts |
| `EnergyConsumer`, `ConformLoad`, `NonConformLoad`, `StationSupply` | loads (SSH values) |
| `SynchronousMachine`, `ExternalNetworkInjection`, `EquivalentInjection` | injections, PV where a local voltage control is active |
| `StaticVarCompensator` | P = 0 reactive injection, Q limits from the Ω ratings |
| `SvVoltage`, `SvPowerFlow`, `SvTapStep` | start values and validation reference |

Everything else is counted in the coverage report rather than silently
dropped. `summarizeCGMES` shows the full class histogram.

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
| `cgmes_import.ignore_connected` | `false` | Diagnostic override treating every terminal as connected — for snapshots whose SSH flags contradict their own SV state. Allowed values: `true`, `false`. |
| `cgmes_import.vset_min_pu` | `0.5` | Lower bound of the plausibility band for a voltage `RegulatingControl.targetValue`, in p.u. of the regulated bus's nominal voltage. |
| `cgmes_import.vset_max_pu` | `1.5` | Upper bound of that band. A target outside `[vset_min_pu, vset_max_pu]` is treated as a placeholder: it is ignored, the unit is held PV at the bus voltage derived from the nominal data, and the substitution is reported as a `warning:`. |
| `cgmes_import.multi_slack` | `true` | Give every electrical island its own SV-declared angle reference (at most one per island). Required for multi-island deliveries; `false` forces the legacy single-reference behavior. Allowed values: `true`, `false`. |

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

### HVDC (Stage-0 model)

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

## Tap control (Stage 2)

With `tap_control = true` the importer starts from the SSH tap positions and
attaches Sparlectra's outer-loop controllers for every tap changer whose
`controlEnabled` and `TapChangerControl.enabled` flags are set — voltage
control on ratio tap changers, active-power control on phase shifters. The
run then goes through `run_sparlectra` (control framework) instead of a plain
`runpf!`.

Two guards apply: a voltage controller whose target bus is held by a generator
(slack or PV) cannot regulate anything and is disabled with a notice; and the
CGMES target deadbands are often wide, so a controller may legitimately settle
one step away from the position the exporting tool recorded.

## Validation against the SV profile

`compareWithSV(result)` compares the solved state with the delivery's own
State Variables profile:

- **Voltages** — per-bus Δvm/Δva against `SvVoltage`, with max and RMS.
- **Flows** (`.flows`) — per-terminal comparison against `SvPowerFlow` in the
  CGMES sign convention: branch terminals, shunts at the solved voltage, loads
  as an SSH↔SV consistency check, and units aggregated per bus. Note that the
  ENTSO-E conformity data sets only ship injection terminals; branch flows
  appear in real exchanges.

De-energized and isolated buses are excluded so the metrics describe the
solved grid.

## Limitations

- CGMES 3.0 is prepared in the schema layer but not yet validated against a
  3.0 delivery.
- Node-breaker modelling is not implemented; it is not needed in practice
  because CGMES deliveries ship the TP profile alongside, which the importer
  reads as bus-branch.
- Tabular phase tap changers are read but not applied to the solved branch.
- HVDC, `AsynchronousMachine` and nonlinear shunts are counted, not mapped.
- Difference models (`dm:DifferenceModel`) are skipped with a report line.

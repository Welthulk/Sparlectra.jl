State Estimation Extensions
============================

These features build on the WLS core of
[State Estimation](state_estimation.md), share its result structure
(`SEResult`) and follow its
[Write-back policy](state_estimation.md#Write-back-policy).

## FACTS in state estimation (`se_view`)

State estimation is a snapshot: `runse!` never invokes the outer control
loop, every controller is frozen at its current operating point. `se_view`
reports the frozen view without mutating the net.

| Device (PF view) | SE view | estimable |
|---|---|---|
| OLTC (in-phase regulator) | transformer with fixed tap, or released tap state | yes, [tap estimation](#Transformer-tap-estimation) with `mode = :ratio` |
| PST (phase-shifting regulator) | transformer with fixed tap, or released tap state | yes, [tap estimation](#Transformer-tap-estimation) with `mode = :pst` (API) |
| Schraegregler | transformer with fixed taps, or released cascade | yes, [tap estimation](#Transformer-tap-estimation) with `mode = :both` (API) |
| SVC / STATCOM | shunt with Q injection | `B` per case A of [shunt estimation](state_estimation.md#Shunt-estimation-and-back-calculation-(ShuntQMeas)) |
| Compensation reactor / capacitor bank | shunt | `B` per case A/B |
| Series compensation | fixed branch impedance | no |

| | |
|---|---|
| Call | `se_view(net)`: the registered controllers, the shunt-estimation releases and their classification, the link clusters the estimator will fuse, and every measurement the WLS will exclude or aggregate |
| Print | `print_se_view(view)` (`format = :markdown`) |
| Artifact | `se_view.md` |

## Links in state estimation

Bus links (impedance-less couplers, `addLink!`) are not part of the Ybus, so
SE runs on the contracted net, like the power flow: every closed-link
cluster is fused onto its representative bus, open links remain real
separations. Link flow measurements are allocation inputs, never WLS rows.

| | |
|---|---|
| Link flow measurement | `addPflowMeasurement!(net; linkNr = ...)` (positive from `link.fromBus` to `link.toBus`); constrains only the flow split |
| Allocation | `calcLinkFlowsSE!(net)` after `with_state_estimation_config(() -> runse!(...); update_net = true)`; the returned rows carry `source = :kcl` or `:measured_ls` plus the residual per link |
| Write-back | `updateNet = true` gives all cluster members the representative's estimated voltage |
| Diagnostics | `LINKAGG` rows appear with `measurement_index = 0` and are never eliminated |

On a fused cluster, `VmMeas`/`VaMeas` on a member remap exactly.
`PinjMeas`/`QinjMeas` aggregate only as a whole: when every member carries
an active injection measurement of the kind, one `LINKAGG` row with the
summed value and root-sum-square sigma replaces them, otherwise all member
injections of that kind are excluded with a warning (a partial sum would
bias the fused balance). `ShuntQMeas` and bus-referenced `ImagMeas` follow
their device to the representative; branch measurements keep their branch,
and a branch shorted inside a cluster loses its flow measurements with a
warning. `calcLinkFlowsSE!` is a weighted least-squares split per link
component (nodal balance rows at weight 1, each link measurement at
`1/sigma^2`; without measurements it equals `calcLinkFlowsKCL!`) and runs
after the write-back because it consumes the estimated branch flows and
shunt powers.

## Island-wise estimation

A net can split into several synchronous AC islands (multi-area MATPOWER
cases, CGMES deliveries with stub islands). Like the power flow, `runse!`
partitions the measurement set onto the islands and estimates every
measured island with its own angle reference; unmeasured islands are skipped
and reported.

| | |
|---|---|
| Call | `runse!` partitions automatically; `evaluate_global_observability` judges each measured island on its own subnet (quality = worst measured island, note `:island_partition`, per-island results in `islands`) |
| Result field | `SEResult.islands`: one row per island (estimated/skipped, iterations, `J`/`dof`, own band verdict `band_reason`, `z_wh`); `SEResult.voltages` indexed by the original bus numbering; `objectiveJ`/`dof` summed over the estimated islands |
| Artifact | the island table in `se_diagnostics.md` (the printed diagnostics show the same) |
| Chain | `runpf_from_se!` solves multi-island nets island-wise automatically; the slack pickup sums over the island references |

Closed links fuse islands into one estimation group; an open link remains
a real separation. An island with its own slack keeps it as reference; one
without gets the MATPOWER-like promotion of the power flow
(`detect_ac_islands`' reference bus), an unpowered stub without any
candidate an angle pin on its first bus. Unestimated islands keep their
start values in `SEResult.voltages`; `validate_measurements` and the
sequential elimination work across islands with the caller's measurement
indices. Independent chi-squares add, so the summed `J`/`dof` is a valid
overall statistic, while the per-island verdict keeps a local problem
visible.

## [Transformer tap estimation](@id se-tap-estimation)

A tap position that disagrees with the model (a stale SCADA position, a
local operation) poisons the estimate around the transformer: the band test
reports `:high` and the suspects cluster there although every telemetry row
is healthy. Releasing the tap as an estimation state resolves this. The
estimator solves with the tap as a continuous state, rounds it to the
nearest mechanical step, and solves once more with the tap fixed; the
reported result is that of the fixed run.

**Use**

| | |
|---|---|
| Release | `setTapEstimation!(net; trafo, mode = :ratio \| :pst \| :both, alpha_deg, enabled = true)`; `trafo` is a branch index or component name, `mode = :ratio` releases the longitudinal regulator $r_1$, `:pst` the skew regulator $r_2$ with its fixed nameplate direction `alpha_deg`, `:both` the pair |
| Write back to the model | `state_estimation.update_taps = true` (keyword `updateTaps`), off by default |
| Machine transformers | `calcMachineTrafoTapFromSE(net; trafo, v_machine_pu, p_mw, q_mvar)`: back-calculation after the run, reports the electrical and nearest mechanical step plus the reactive-power residual, writes nothing |
| Result | `SEResult.tapEstimates` (per released transformer: branch index and name, mRID for CGMES cases, branch and bus numbers for MATPOWER/DTF, model step, continuous electrical step, fixed mechanical step, change, range flag, freeze reason, island id on multi-island nets), `SEResult.tapFixation` (`j_before`/`dof_before` with the taps still states, `j_after`/`dof_after` with the taps fixed, `offgrid_residual`) |
| Artifact | `se_tap_estimates.csv` (model step next to the fixed step and the change between them) |
| Web UI | "estimate taps" (on by default) releases every in-service transformer with a ratio tap changer (mode `:ratio`) except machine transformers; PST releases (`:pst`, `:both`) are an API-level choice. Result table: [Web UI](webui_reference.md#State-estimation) |

**Notes**

- The regulator states start from the current branch position; a position
  contradicting the declared mode is an error. Fixation rounds regulator 1
  on the fraction grid `tap_step` and regulator 2 on the shift grid
  `phase_step_deg`; out-of-range positions are flagged and clamped.
  `SEResult.voltages`/`objectiveJ`/`dof` are those of the fixed run.
- Freeze reasons: `radial_no_voltage_pin` for a bridge transformer (its
  removal cuts the net) whose cut-off side carries no active voltage
  measurement, `not_observable` for a tap state the measurements cannot
  pin (a `:both` release can freeze partially). A frozen transformer solves
  bitwise identical to no release and still appears in `tapEstimates` with
  its reason. Release and back-calculation are mutually exclusive per
  transformer.
- The Web UI's "estimate taps" first solves once with the taps frozen and
  then estimates them from that state (from a flat start many tap states
  crawl to the solution). The estimate is kept only if its fixed positions
  fit the measurements at least as well as the model positions of that
  first solve. If it is rejected, or the run does not converge with
  released taps (a set collectively too thin for them), every released tap
  is frozen back to its model position and the estimation repeats once;
  the run log says so and the result metadata carries
  `se_tap_estimation_fallback = true` with
  `se_tap_estimation_fallback_reason` `not_converged` or `worse_fit`: the
  tap positions in such a result are model values, not estimates.
- A band failure that appears only through the fixation is not bad data:
  when the continuous fit passes or undershoots the band and the fixed run
  fails it `:high`, `tapFixation.offgrid_residual` is set (summary note
  `:offgrid_tap_residual`). The true position then sits between steps, or
  the step table (`tap_step`, neutral position) is wrong.

**Cascade model.** The released tap enters as one or two WLS states with

```math
t(r_1, r_2) = \frac{t_{base}}{(1 + r_1)\,(1 + r_2\, e^{j\alpha})}
```

in the `calcSkewAngleTap`/`calcAdmittance` convention of the importers;
every prediction re-adds the branch at the cascade position, so at the
initial position the released model reproduces the stamped power flow to
machine precision.

**Fixation and dof.** The fixation raises `dof`: a released tap is one
extra state, `dof_before = m - (voltage states + taps)`, and the fixation
removes it, `dof_after = m - voltage states = dof_before + released taps`.
`J` moves the other way, e.g. `J 58.5 (dof 41) -> J 62.5 (dof 42)`, both
inside their band; only a fixation that jumps out of the band is a finding.
On multi-island nets the sums aggregate like the island chi-squares.

**Machine transformers.** A generator step-up transformer's machine
terminal is AVR-held and not independently observed, so a released tap
would only absorb the voltage control; the Web UI mass release skips such
transformers (a generator bus hanging on the transformer alone), an
explicit `setTapEstimation!` stays possible. The back-calculation uses the
estimated network-side voltage, the AVR setpoint, the dispatch P and the
measured machine reactive power.

**Guards.** Tap states are outside the observability count (like the gated
current rows); the bridge test and the per-column numerical test check the
redundancy around the transformer (a meshed path or measurements on both
sides) directly. A frozen regulator keeps its exact model value inside a
written cascade.

## Topology validation

A wrong service state (a breaker recorded closed that is open, or the
reverse) poisons an estimation differently from bad telemetry: the errors
cluster around one station and no elimination cures them. Three advisory
stages validate the topology: findings and recommendations only, never a
mutation of a status, a measurement or the model.

| | |
|---|---|
| Stage 1, pre-checks | `validate_topology(net, measurements)` runs linearly before any estimation, automatic at the start of `runse!` with `state_estimation.topology_precheck = true` (the default); findings land in `SEResult.topologyFindings` and warn, the result is logged on every service run |
| Stage 2, classification | `runse_diagnostics` reports `:topology_error_suspected_at_station` when the sequential elimination exhausts its budget while the band test still fails `:high` and the surviving suspects cluster at one station (at least `state_estimation.topology_cluster_min`) |
| Stage 3, hypothesis test | `test_topology_hypotheses(net, measurements; candidates = :auto, max_candidates)` toggles each candidate's service state on a working copy, re-runs the estimation, and reports J and the band z-score before versus after, ranked |
| Config keys | thresholds are sigma multiples (`state_estimation.topology_open_flow_k` and friends, see [configuration](state_estimation_configuration.md)) |
| Web UI | the stage-3 test sits behind a button on the SE result page and writes `topology_hypotheses.md` |

Stage-1 findings: a measured flow over an open element
(`:open_element_with_flow`), a closed branch reading dead at both measured
ends while the model expects a clear flow (`:closed_element_without_flow`,
low severity), voltage measurements disagreeing across a closed link
(`:closed_link_voltage_mismatch`), and the node balance at completely
measured nodes (`:kcl_violation`; shunt term from the measured voltage,
partially measured nodes skipped). Out-of-service elements enter only the
status-contradiction check.

Stage 2 takes a station as a closed-link contraction cluster; stage-1
findings at the same station are noted as `precheck_agreement`, and with
correlation columns on a fully correlated cluster carries
`correlated_group`. A curable gross error never produces a topology
finding, because the elimination succeeds first.

Stage 3 supports a hypothesis when the `:high` failure disappears under the
toggle (inside or below the band); several supported hypotheses are flagged
ambiguous. Candidates come from the stage-2 stations and the stage-1 status
contradictions (`candidates = :auto`, capped by `max_candidates`, with
per-candidate timing), or from an explicit list. Nothing is switched
automatically.

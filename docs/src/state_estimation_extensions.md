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
| Schraegregler | transformer with fixed taps, or released tap states | yes, [tap estimation](#Transformer-tap-estimation) with `mode = :both` (API) |
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
| Link flow measurement | `addPflowMeasurement!(net; linkNr = ...)` (positive from `link.fromBus` to `link.toBus`), or a `P`/`Q` row on a link terminal of the measurement table; constrains only the flow split |
| Allocation | `calcLinkFlowsSE!(net)` after `with_state_estimation_config(() -> runse!(...); update_net = true)`; the returned rows carry `source = :kcl` or `:measured_ls` plus the residual per link. The SE service calls it and writes `se_link_flows.csv` |
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
component (the nodal balance rows as equations, each link measurement at
`1/sigma^2`, so a measurement moves only the circulations of a ring and a
flow the balance fixes keeps its value; without measurements it equals
`calcLinkFlowsKCL!`) and runs
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
is healthy. Releasing the tap changers as estimation states resolves this.
The estimator solves with one continuous step per tap changer, rounds each
to the nearest physical step, and solves once more with the steps fixed; the
reported result is that of the fixed run.

**Use**

| | |
|---|---|
| Release | `setTapEstimation!(net; trafo, mode = :ratio \| :pst \| :both, alpha_deg = nothing, enabled = true)`; `trafo` is a branch index or component name. `mode` selects the tap changers of the transformer by role: `:ratio` the ratio changers, `:pst` the phase changers, `:both` all. `:pst` and `:both` need a declared phase changer, `:ratio` a ratio changer (`ArgumentError` otherwise). `alpha_deg` is optional and, when given, has to equal the direction of the asymmetrical changer from the import. Returns the start steps `s0` |
| Write back to the model | `state_estimation.update_taps = true` (keyword `updateTaps`), off by default: ratio, shift and the stored steps of the fixed run, and r, x, g, b when an impedance source applied (see Impedance per step) |
| Machine transformers | `calcMachineTrafoTapFromSE(net; trafo, v_machine_pu, p_mw, q_mvar)`: back-calculation of the step of the first ratio changer after the run (`s_est`, the electrical and the nearest physical step, the reactive-power residual), writes nothing |
| Result | `SEResult.tapEstimates`, one row per released tap changer: branch, `regulator` (index in the changer set), transformer name and mRID, `changer` and `changer_mrid`, `kind`, `role`, `side`, `alpha_deg`, `du`, `neutral_step`, the estimate `s_est` with `du_est`, `ratio_est`, `angle_est_deg`, `s_fixed` (the position of the fixed run), `step_model`, `step_fixed`, `deviation_steps`, `out_of_range`, `fixed`, `frozen_reason`, `note`, `notice`, and `island` on multi-island nets. `SEResult.tapFixation`: `fixed` (false when the steps were not fixed), `not_converged_released` (the released run did not converge and the taps were frozen again, see Fixation), `j_before`/`dof_before` with the steps as states, `j_after`/`dof_after` with the steps fixed, the iterations of both runs, `offgrid_residual`, `k_removed`, `delta_j_notice`, and `j_rise_exceeds_quantile`, `sequential` and `j_simultaneous` for the sequential fixation (see Fixation) |
| Artifact | `se_tap_estimates.csv`, one row per regulator: changer, kind, side, `alpha_deg`, `du`, the estimate `s_est` with `du_est`, `ratio_est`, `angle_est_deg`, the model step next to the fixed step and the deviation between them in steps, freeze reason and note |
| Web UI | "estimate taps" (on by default when the case has tap changers) releases every tap changer of every in-service transformer except machine transformers, in the mode its changers give. Result table: [Web UI](webui_reference.md#State-estimation) |

**State and changer kinds.** The state of a tap changer is its electrical
step `s`, continuous and relative to the neutral step; the additional
voltage is `Delta U = s * du`. The changer data come from the tap changer
set of the branch: what the CGMES or PowSyBl import stored, the typed models
of a resolved winding, or the nameplate fields (MATPOWER
`mpc.sparlectra.tap_changers`, SCF). Each changer gives the factor `f(s)` of
its winding:

| Kind | Typical source | Factor |
|---|---|---|
| ratio | `RatioTapChanger`, a nameplate ratio step, the DTF Laengsregler | $f = 1 + s\,du$ |
| asymmetrical | `PhaseTapChangerAsymmetrical`, a nameplate Delta-u PST, the DTF Schraegregler | $f = 1 + s\,du\,e^{j\alpha}$ |
| asymmetrical at 90 degrees | quadrature booster, the DTF Querregler | $f = 1 + j\,s\,du$ |
| symmetrical | `PhaseTapChangerSymmetrical` | $f = e^{j\varphi}$ with $\varphi = 2\arctan(s\,du/2)$ |
| tabular | `RatioTapChangerTable`, `PhaseTapChangerTabular`, IIDM step tables | the table points, ratio and angle interpolated linearly between them |
| linear | `PhaseTapChangerLinear`, a nameplate phase step in degrees | $f = e^{j s\,\Delta\varphi}$ |

The factor multiplies the branch ratio (changer on the from side) or divides
it (to side), in the order of the import, starting from the ratio and shift
with every changer at neutral; a linear ratio grid (`ratio_grid` 1 in
MATPOWER, `grid: "linear"` in SCF) adds `s * step` instead. The direction
alpha, `du` and the tables are data from the import and never states. At
the start position the released model reproduces the stamped power flow to
machine precision. The Jacobian column of a step is a forward difference
along the analytic tangent `dt/ds` (with the impedance part when an
impedance source applies), with a difference step in `s` of `jac_eps` over
the size of one step, so `t` moves by about `jac_eps`. At most two tap
changers of one transformer can be told apart, because the complex ratio
has two real directions; a further changer, or one whose effect equals that
of another (two ratio changers, ratio and asymmetrical at alpha 0), is
frozen with `not_separable`. A transformer without tap changer has an empty
set and is never released.

**Fixation.** After convergence every estimated step is rounded to the
nearest physical step (ties away from zero; not a regulator with note
`relative_to_frozen`, see below) and clamped to the band of the
changer, read absolutely from its data (`lowStep` to `highStep`), not
relative to the import step. A phase table fixes on its nearest point, a
`RatioTapChangerTable` on its points inside the table and on whole steps of
its formula outside. Then one run with the steps fixed follows.
`out_of_range` marks an estimate more than half a step outside the band. A
frozen regulator keeps its start position exactly.
`SEResult.voltages`/`objectiveJ`/`dof` are those of the fixed run.
`tapFixation.delta_j_notice` gives `j_after - j_before` with `k` (the step
states the fixed run removed) and the 95 % quantile of chi2(k); it is a
notice, not part of any acceptance rule.

When the fixed run converged, the fixation raises J by more than the 95 % quantile of chi2(k) and
both runs gate the same Ia rows, the run logs the fixation line as a warning
instead of an info line and sets `tapFixation.j_rise_exceeds_quantile`; the
notice then names the rise and the quantile, with `k` of 2 or more also both J. With `k` of 2 or more
the steps are fixed once more one at a time, starting from the continuous
estimate: the step closest to its physical position first, the others
estimated again, and so on. The sequential result is kept when every solve
converged, it gates the same Ia rows and its J is smaller
(`tapFixation.sequential`; the J of the fixation at once is in
`tapFixation.j_simultaneous`). Strongly coupled regulators, such as two
ratio changers at one bus, can need this: rounded together they can land on
a pair far from the data although the true pair fits. The SE service
carries both flags in the result metadata (`se_tap_j_rise_exceeds_quantile`,
`se_tap_sequential_fixation`).

**Freeze reasons** (`frozen_reason`, with a `notice` text per row and a
`notice:` line in the run log of the SE service):

- `radial_no_voltage_pin`: a bridge transformer (its removal cuts the net)
  whose cut-off side carries no active voltage measurement; all its
  changers freeze.
- `not_separable`: see above.
- `no_angle_anchor`: a phase changer on a tree edge of the angle graph
  (its far side is tied to the reference by no Va measurement, no measured
  end of an Ia row and no meshed path) when the same transformer also has
  a released ratio changer. The ratio changer stays released.
- `symmetric_start`: a phase changer on such a tree edge without a
  released ratio changer, at a start where `|t|` does not move
  (asymmetrical at its turning point, for example 90 degrees with
  `s0 = 0`; symmetrical or linear; a table at a vertex or on a flat
  piece of its ratio, so a table with constant ratio at every start). A
  phase changer on such an edge that is not at a symmetric start stays
  released.
- `not_observable`: no measurement sees the step. Its Jacobian column at
  the start point is zero (rank test on the column scaled to norm 1) or
  zero within the rounding noise of its forward difference on every row
  (per row about `10 eps` times the size of the terms of the row, taken
  from its voltage derivatives, times the step size over `jac_eps`; a zero
  injection has terms of thousands of MW). The test is per row and does
  not read the weights, so neither a heavy row elsewhere nor a heavy row
  that sees only noise freezes a regulator that another row sees.
- `model_mismatch`: the branch ratio and angle do not lie on the curves of
  its changer data. The release is recorded, the regulator is frozen.
- `not_converged_released`: the estimate with the released taps did not
  converge (a set too thin for them, or a gross error that a regulator
  pinned only by its magnitude effect tries to absorb). Every regulator
  that was still a state is frozen again at its start position and the
  estimate is solved once more with the taps frozen, from the same start;
  `tapFixation.not_converged_released` is set, the run warns once, and the
  result is that frozen-tap estimate, never the diverged iterate.
- `no_changer`: a released transformer without a changer of the released
  role gets one report row (regulator 0), which the counts of the SE
  service leave out.

Notes (`note`) that do not freeze: `relative_only` for parallel members
(estimated only relative to the other member), `conditional` on the
ratio changer of a `no_angle_anchor` transformer when the frozen phase
changer moves `|t|` (asymmetrical, or a table whose ratio varies; not for
symmetrical or linear), `sign_from_start` for a released phase changer on
a tree edge without an anchor whose `|t|` has a turning point on its band
(it keeps the side of its start: a step that would cross the turning
point goes half the way and stays on the turning point once it is that
close), `relative_to_frozen` for a phase regulator whose loop closes only
over a frozen phase regulator (frozen for any reason, for example
`no_angle_anchor` or `model_mismatch`): it is only determined relative to
that frozen one, so it is not fixed (`fixed = false`), stays a continuous
state in the fixed run (its `s_fixed` is that continuous value, rounding it
would miss the measured loop angle), its notice names the frozen
transformer, `k_removed` does not count it, and `update_taps` keeps its
model position. Notices of the estimate: a phase changer on a tree edge without an
anchor whose estimate sits near the vertex (only the magnitude effect pins
it, its step is not reliable); a CGMES changer fixed at a step whose factor
leaves 0.5 to 2 (a re-import ignores it there); a released column more than
`1/jac_eps` times weaker than the strongest voltage column gets "weak
column, factor x below the strongest", also in the run log. A frozen
regulator solves bitwise identical to no release and still appears in
`tapEstimates` with its reason. Release and back-calculation are mutually
exclusive per transformer.

**Fallback of the SE service.** The Web UI's "estimate taps" first solves
once with the taps frozen and then estimates them from that state (from a
flat start many tap states crawl to the solution). The estimate is kept
only if its fixed positions fit the measurements at least as well as the
model positions of that first solve. If it is rejected, or the run does not
converge with released taps (`runse!` then already froze them again,
reason `not_converged_released`), every released tap is frozen back to its
model position and the estimation repeats once; the run log says so and the
result metadata carries
`se_tap_estimation_fallback = true` with
`se_tap_estimation_fallback_reason` `not_converged` or `worse_fit`: the tap
positions in such a result are model values, not estimates.

**Off-grid positions.** A band failure that appears only through the
fixation is not bad data: when the continuous fit passes or undershoots the
band and the fixed run fails it `:high`, `tapFixation.offgrid_residual` is
set (summary note `:offgrid_tap_residual`). The true position then sits
between steps, or the changer data (step size, neutral step) are wrong.

**Impedance per step.** `state_estimation.tap_estimation.impedance_correction`
(`auto`, `off`; default `auto`) decides whether a released regulator moves
the transformer impedance with its step. `auto` takes the per-step data the
model has: the r, x, g, b of a tap table row, the reactance curve of a phase
tap changer with `xMin`/`xMax`, the factors of a `TapStepImpedance`, and the
impedance correction of typed models under
`model.tap_changer_model = impedance_correction`; `off` keeps the impedance
of the start position. There is no `on`: without per-step data no curve
exists to force. With a correction the Jacobian column of the regulator
includes the part through the impedance.

**Fixation and dof.** The fixation raises `dof`: every estimated step is
one extra state, `dof_before = m - (voltage states + steps)`, and the
fixation removes them, `dof_after = m - voltage states = dof_before + k`.
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
current rows); the bridge test, the angle anchors and the per-column test
check the redundancy around the transformer (a meshed path or measurements
on both sides) directly. A frozen regulator keeps its exact model value.

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

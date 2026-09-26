State Estimation Extensions
============================

These estimator features build on the WLS core of
[State Estimation](state_estimation.md) and share its result structure
(`SEResult`). Write-back follows the
[Write-back policy](state_estimation.md#Write-back-policy) of the core.

## FACTS in state estimation (`se_view`)

State estimation is a snapshot: `runse!` never invokes the outer control
loop, every controller is frozen at its current operating point. The frozen
view is reported without mutating the net.

| Device (PF view) | SE view | estimable |
|---|---|---|
| OLTC (in-phase regulator) | transformer with fixed tap, or released tap state | yes, [tap estimation](#Transformer-tap-estimation) with `mode = :ratio` |
| PST (phase-shifting regulator) | transformer with fixed tap, or released tap state | yes, [tap estimation](#Transformer-tap-estimation) with `mode = :pst` (API) |
| Schraegregler | transformer with fixed taps, or released cascade | yes, [tap estimation](#Transformer-tap-estimation) with `mode = :both` (API) |
| SVC / STATCOM | shunt with Q injection | `B` per case A of [shunt estimation](state_estimation.md#Shunt-estimation-and-back-calculation-(ShuntQMeas)) |
| Compensation reactor / capacitor bank | shunt | `B` per case A/B |
| Series compensation | fixed branch impedance | no |

**Use**

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

**Use**

| | |
|---|---|
| Link flow measurement | `addPflowMeasurement!(net; linkNr = ...)` (positive from `link.fromBus` to `link.toBus`); constrains only the flow split |
| Allocation | `calcLinkFlowsSE!(net)` after `with_state_estimation_config(() -> runse!(...); update_net = true)`; the returned rows carry `source = :kcl` or `:measured_ls` plus the residual per link |
| Write-back | `updateNet = true` gives all cluster members the representative's estimated voltage |
| Diagnostics | `LINKAGG` rows appear with `measurement_index = 0` and are never eliminated |

**Notes** (measurement rules on a fused cluster)

- `VmMeas`/`VaMeas` on a member remap exactly (voltage equality).
- `PinjMeas`/`QinjMeas` aggregate only as a whole: when every member carries
  an active injection measurement of the kind, one `LINKAGG` row with the
  summed value and root-sum-square sigma replaces them; otherwise all member
  injections of that kind are excluded with a warning.
- `ShuntQMeas` and bus-referenced `ImagMeas` follow their device to the
  representative; branch measurements keep their branch.
- A branch shorted inside a cluster loses its flow measurements with a
  warning.

!!! details "Why it is built this way"
    A partial sum of member injections would bias the fused balance, hence
    whole-cluster aggregation or exclusion. `calcLinkFlowsSE!` extends the
    KCL allocation to a weighted least-squares split per link component:
    nodal balance rows keep weight 1, each link measurement adds a row with
    weight `1/sigma^2`, and without measurements the result equals
    `calcLinkFlowsKCL!` exactly. It runs after the write-back because the
    estimate fills the branch flows and shunt powers it consumes.

## Island-wise estimation

A net can split into several synchronous AC islands (multi-area MATPOWER
cases, CGMES deliveries with stub islands). Like the power flow, `runse!`
partitions the measurement set onto the islands and estimates every
measured island with its own angle reference; unmeasured islands are skipped
and reported.

**Use**

| | |
|---|---|
| Call | `runse!` partitions automatically; `evaluate_global_observability` judges each measured island on its own subnet (quality = worst measured island, note `:island_partition`, per-island results in `islands`) |
| Result field | `SEResult.islands`: one row per island (estimated/skipped, iterations, per-island `J`/`dof`, own Wilson-Hilferty verdict `band_reason`, `z_wh`); `SEResult.voltages` indexed by the original bus numbering; `objectiveJ`/`dof` summed over the estimated islands |
| Artifact | the island table in `se_diagnostics.md` (the printed diagnostics show the same) |
| Chain | `runpf_from_se!` solves multi-island nets island-wise automatically; the slack pickup sums over the island references |

**Notes**

- Closed links fuse islands: islands connected by a closed busbar coupler are
  one estimation group; an open link remains a real separation.
- Reference per island: an island with its own slack keeps it; one without
  gets the MATPOWER-like promotion the power flow uses
  (`detect_ac_islands`' reference bus). An unpowered stub island without any
  candidate still gets an angle pin on its first bus (only the angle is
  fixed, the magnitude stays a state).
- Unestimated islands keep their start values in `SEResult.voltages`.
- `validate_measurements` and the sequential elimination work across islands
  with the caller's measurement indices.

!!! details "Why it is built this way"
    Independent chi-squares add, so the summed `J`/`dof` is a valid overall
    statistic, while the per-island verdict keeps a local problem visible
    that the sum would average out.

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
| Write back to the model | `state_estimation.update_taps = true` (keyword `updateTaps`), off by default; see [Write-back policy](state_estimation.md#Write-back-policy) |
| Machine transformers | `calcMachineTrafoTapFromSE(net; trafo, v_machine_pu, p_mw, q_mvar)`: back-calculation after the run, reports the electrical and nearest mechanical step plus the reactive-power residual, writes nothing |
| Result | `SEResult.tapEstimates` (one row per released transformer: branch index and name, mRID for CGMES cases, branch and bus numbers for MATPOWER/DTF; model step, continuous electrical step, fixed mechanical step, change, range flag, freeze reason; on multi-island nets the island id), `SEResult.tapFixation` (`j_before`/`dof_before` with the taps still states, `j_after`/`dof_after` with the taps fixed, `offgrid_residual`) |
| Artifact | `se_tap_estimates.csv` (model step next to the fixed step and the change between them, so a wrong model position reads as "corrected") |
| Web UI | "estimate taps" releases every in-service transformer with a ratio tap changer (mode `:ratio`) except machine transformers (on by default); PST releases (`:pst`, `:both`) are an API-level choice. The result page reports per transformer the model step, the continuous electrical step and the fixed mechanical step, plus `J` before versus after the fixation: a small `J` before with a large `J` after means the true position sits between mechanical steps, both small means the fixed step explains the measurements |

**Notes**

- The regulator states start from the current branch position; a position
  contradicting the declared mode is an error.
- Fixation rounds regulator 1 on the fraction grid `tap_step` and regulator
  2 on the shift grid `phase_step_deg`; out-of-range positions are flagged
  and clamped. `SEResult.voltages`/`objectiveJ`/`dof` are those of the fixed
  run.
- A bridge transformer (its removal cuts the net) whose cut-off side carries
  no active voltage measurement is frozen at the current position
  (`radial_no_voltage_pin`); a remaining tap state the measurements cannot
  pin fails the per-column numerical test of the released shunt states and
  freezes with `not_observable` (a `:both` release can freeze partially). A
  transformer with all regulators frozen solves bitwise identical to no
  release and still appears in `tapEstimates` with its reason. Release and
  back-calculation are mutually exclusive per transformer.
- If the run does not converge with released taps (a set collectively too
  thin for the released taps, which the single-tap guards cannot see),
  every released tap is frozen back to its model position and the
  estimation repeats once; the run log says so and the result metadata
  carries `se_tap_estimation_fallback = true`, because the tap positions in
  such a result are model values, not estimates.
- A band failure that appears only through the fixation is not bad data:
  when the continuous fit passes or undershoots the band and the fixed run
  fails it `:high`, `tapFixation.offgrid_residual` is set (summary note
  `:offgrid_tap_residual`). The true position then sits between steps, or
  the step table (`tap_step`, neutral position) is wrong; the bad-data
  workflow is the wrong move there.

!!! details "Why it is built this way"
    **Cascade model.** The released tap enters as one or two WLS states
    with

    ```math
    t(r_1, r_2) = \frac{t_{base}}{(1 + r_1)\,(1 + r_2\, e^{j\alpha})}
    ```

    built from the `calcSkewAngleTap`/`calcAdmittance` convention of the
    importers. The stamped admittance terms leave the Ybus once and every
    prediction re-adds the branch at the cascade position, so at the
    initial position the released model reproduces the stamped power flow
    to machine precision.

    **Mandatory fixation and dof.** A tap changer sits on a mechanical
    grid, so the continuous estimate is rounded and re-solved. The fixation
    makes `dof` go up, not down: `dof = m - n` counts the redundancy of the
    residuals, not the model's freedom. A released tap is one extra state,
    `dof_before = m - (voltage states + taps)`; the fixation removes it and
    `dof_after = m - voltage states = dof_before + released taps`. `J`
    moves the other way, e.g. `J 58.5 (dof 41) -> J 62.5 (dof 42)`, both
    inside their band; both tiny means the fixed step explains the
    measurements, and only a fixation that jumps out of the band is a
    finding. On multi-island nets the sums aggregate like the island
    chi-squares.

    **Machine transformers.** A generator step-up transformer's machine
    terminal is AVR-held and not independently observed, so a released tap
    would only absorb the voltage control. The back-calculation uses the
    estimated network-side voltage, the AVR setpoint, the dispatch P and
    the measured machine reactive power (under AVR control Q is telemetry,
    not a schedule value). The Web UI mass release skips machine
    transformers (a generator bus hanging on the transformer alone); an
    explicit `setTapEstimation!` stays possible.

    **Guards instead of observability count.** Tap states are outside the
    observability count (like the gated current rows), so dedicated guards
    protect the release: estimating a tap needs redundancy around the
    transformer, a meshed path or measurements on both sides, which the
    bridge test and the per-column numerical test check directly. A frozen
    regulator keeps its exact model value inside a written cascade.

## Topology validation

A wrong service state (a breaker recorded closed that is open, or the
reverse) poisons an estimation differently from bad telemetry: the errors
cluster around one station and no elimination cures them. Three stages
validate the topology, all advisory: findings and recommendations only,
never a mutation of a status, a measurement, or the model.

**Use**

| | |
|---|---|
| Stage 1, pre-checks | `validate_topology(net, measurements)` runs linearly before any estimation, automatic at the start of `runse!` with `state_estimation.topology_precheck = true` (the default); findings land in `SEResult.topologyFindings` and warn, the result is logged on every service run |
| Stage 2, classification | `runse_diagnostics` reports `:topology_error_suspected_at_station` when the sequential elimination exhausts its budget while the band test still fails `:high` and the surviving suspects cluster at one station (at least `state_estimation.topology_cluster_min`) |
| Stage 3, hypothesis test | `test_topology_hypotheses(net, measurements; candidates = :auto, max_candidates)` toggles each candidate's service state on a working copy, re-runs the estimation, and reports J and the band z-score before versus after, ranked |
| Config keys | thresholds are sigma multiples (`state_estimation.topology_open_flow_k` and friends, see [configuration](state_estimation_configuration.md)) |
| Web UI | the stage-3 test sits behind a button on the SE result page and writes `topology_hypotheses.md` |

**Notes**

- Stage-1 findings: a measured flow over an open element
  (`:open_element_with_flow`), a closed branch reading dead at both measured
  ends while the model expects a clear flow (`:closed_element_without_flow`,
  low severity), voltage measurements disagreeing across a closed link
  (`:closed_link_voltage_mismatch`), and the node balance at completely
  measured nodes (`:kcl_violation`; shunt term from the measured voltage,
  partially measured nodes skipped).
- Out-of-service elements enter only the status-contradiction check.
- Stage 2: station = closed-link contraction cluster. Stage-1 findings at
  the same station are noted as `precheck_agreement`; with correlation
  columns on, a fully correlated cluster carries `correlated_group`.
- Stage 3: a hypothesis is supported when the `:high` failure disappears
  under the toggle (inside or below the band); several supported hypotheses
  are flagged ambiguous. Candidates come from the stage-2 stations and the
  stage-1 status contradictions (`candidates = :auto`, capped by
  `max_candidates`, with per-candidate timing), or from an explicit list.
- Nothing is switched automatically.

!!! details "Why it is built this way"
    A curable gross error never produces a topology finding, because the
    elimination succeeds first; the exhausted elimination against a `:high`
    band with suspects clustered at one station is the fingerprint of a
    wrong service state, which is why the classification waits for that
    state. The pre-checks are linear and cheap, so they run on every
    estimation; the hypothesis test re-estimates per candidate and is
    therefore explicit and capped.

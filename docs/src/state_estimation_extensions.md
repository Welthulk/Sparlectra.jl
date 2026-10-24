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
| OLTC (in-phase regulator) | transformer with fixed tap, or released tap state | yes, [tap estimation](@ref se-tap-estimation) with `mode = :ratio` |
| PST (phase-shifting regulator) | transformer with fixed tap, or released tap state | yes, [tap estimation](@ref se-tap-estimation) with `mode = :pst` (API) |
| Schraegregler | transformer with fixed taps, or released tap states | yes, [tap estimation](@ref se-tap-estimation) with `mode = :both` (API) |
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
separations. Link flow measurements on closed links are WLS rows; on an open
link they are excluded with a warning.

| | |
|---|---|
| Link flow measurement | `addPflowMeasurement!(net; linkNr = ...)` / `addQflowMeasurement!` (positive from `link.fromBus` to `link.toBus`), or a `P`/`Q` row on a link terminal of the measurement table; a WLS row with weight `1/sigma^2` |
| Ring circulations | one state per independent loop of a coupler ring that its measured links see, P and Q apart; `SEResult.linkCirculations` (`nothing` on a net with more than one AC island, the link flows carry them there) |
| Link flows | `calcLinkFlowsSE!(net)` after `with_state_estimation_config(() -> runse!(...); update_net = true)`: the balance split at the estimated state plus the estimated circulations; rows carry `source = :estimated` (measured link or estimated circulation) or `:kcl`, plus the residual `f - z` per active link measurement. The SE service calls it and writes `se_link_flows.csv` |
| Write-back | `updateNet = true` gives all cluster members the representative's estimated voltage |
| Diagnostics | link rows are ranked and eliminated like branch flows; `LINKAGG` rows appear with `measurement_index = 0` and are never eliminated |

**Method.** For a cluster with member buses $k$ and closed links $j$, $A$ is
the incidence matrix ($+1$ at the from bus, $-1$ at the to bus) and the link
flows $f$ solve $A f = b$ with

```math
b_k = P^{\mathrm{spec}}_k - P^{\mathrm{sh}}_k(x) - \sum_{\text{branches at } k} P^{\mathrm{br}}_{k}(x)
```

(generation minus load of the member, its shunt draw at the cluster
voltage, the flows of the branches leaving it, including a branch shorted
inside the cluster; the same for Q). These balances are equations. The
general solution is $f = A^{+} b(x) + N c$: $A^{+} b(x)$ is the part the
balances fix (on a tree the whole flow, on a ring the minimum-norm split),
$N$ spans the loops of the cluster. A link row is

```math
h_j(x, c) = \left(A^{+} b(x)\right)_j + (N c)_j
```

with weight $1/\sigma_j^2$. For each of P and Q the circulation states $c$
are the loop directions the measured links see (the row space of $N$
restricted to the measured links); unmeasured loop directions stay at the
minimum-norm split. A single coupler or a tree of couplers carries no
circulation state, so its row depends on the fused state alone: one wrong
meter pulls the state in proportion to its weight, shows a normalized
residual, and once the bad-data loop eliminates it the flow is the balance
value again. On a ring a single measured link pins the circulation and is
critical; with more measured links than loop directions the rows become
redundant and a gross error on one of them is localizable. After the
estimate `calcLinkFlowsSE!` writes $A^{+} b(\hat x) + N \hat c$ to the
links, with $b$ from the estimated branch flows and shunt powers of the
write-back. `examples/state_estimation/state_estimation_link_measurements.jl`
shows both cases.

On a fused cluster, `VmMeas`/`VaMeas` on a member remap exactly.
`PinjMeas`/`QinjMeas` aggregate only as a whole: when every member carries
an active injection measurement of the kind, one `LINKAGG` row with the
summed value and root-sum-square sigma replaces them, otherwise all member
injections of that kind are excluded with a warning (a partial sum would
bias the fused balance). `ShuntQMeas` and bus-referenced `ImagMeas` follow
their device to the representative; branch measurements keep their branch,
and a branch shorted inside a cluster loses its flow measurements with a
warning. `calcLinkFlowsSE!` runs after the write-back because it consumes
the estimated branch flows and shunt powers; a released shunt counts with
the estimated B only with `update_shunts = true`. It adds the estimated
circulations only while the node voltages are still those of the estimate.
On a net without a preceding estimate, or after a later power flow, it
splits by weighted least squares instead: each link measurement at
`1/sigma^2`, the nodal balances as equality constraints on the same KKT
system as `zero_injection = constraint`, so a measurement moves only the
circulations of a ring and an unmeasured circulation keeps the minimum-norm
split (`source = :measured_ls`); without measurements it equals
`calcLinkFlowsKCL!`.

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
| Exclusion | `state_estimation.tap_estimation.exclude` (names or mRIDs): released changers of these transformers stay out, `reason = :excluded`, not calculated. Machine transformers stay out unless `include_machine_transformers = true`, see [Tap estimation per transformer](@ref se-tap-exclude) |
| Machine transformers | `runse!` calculates the step of the first ratio changer of a machine transformer it kept out after a converged run with `update_net` (`status = :calculated`, `reason = :machine_transformer`; without a Q measurement on every bus of the machine node `:not_estimated` with `reason = :q_not_measured_in_group`, without a P measurement `note = :p_from_schedule`); `calcMachineTrafoTapFromSE(net; trafo, v_machine_pu, p_mw, q_mvar)` runs the same calculation for a transformer that is not released (`s_est`, the electrical and the nearest physical step, the reactive-power residual), writes nothing |
| Result | `SEResult.tapEstimates`, one row per released tap changer: branch, `regulator` (index in the changer set), transformer name and mRID, `changer` and `changer_mrid`, `kind`, `role`, `side`, `alpha_deg`, `du`, `neutral_step`, the estimate `s_est` with `du_est`, `ratio_est`, `angle_est_deg`, `s_fixed` (the position of the fixed run), `step_model`, `step_fixed`, `deviation_steps`, `sigma_s` (standard deviation of the step from the covariance) with `round_distance_sigma` (distance to the rounded step in units of it), `step_from_target` (see Reasons), `out_of_range`, `fixed`, `reactance_source` (`:table`, `:analytic`, `:nominal`, see Impedance per step), `status` (`:estimated`, `:out_of_range` for an estimate too far outside the band, see Fixation, `:not_estimated`, or `:calculated` for a machine transformer) with its `reason`, `note`, `step_assumed` (the import carried no step, `step_model` is the middle position), `notice`, and `island` on multi-island nets. `SEResult.tapFixation`: `fixed` (false when the steps were not fixed), `fallback_used` with `fallback_reason` (the run with tap estimation did not converge and the result is the run without it, see [Tap estimation switch and fallback](@ref se-tap-fallback)), `j_before`/`dof_before` with the steps as states, `j_after`/`dof_after` with the steps fixed, the iterations of both runs, `offgrid_residual`, `k_removed`, `delta_j_notice`, and `j_rise_exceeds_quantile`, `sequential` and `j_simultaneous` for the sequential fixation (see Fixation) |
| Artifact | `se_tap_estimates.csv`, one row per regulator, columns below; the metadata key `se_tap_out_of_range` counts the rows with status `out_of_range` |
| Web UI | "estimate taps" (`state_estimation.tap_estimation.enabled`, on by default) releases every declared tap changer of every in-service transformer, in the mode its changers give; the Case page sets per transformer whether it is estimated. Result table: [Web UI](webui_reference.md#State-estimation) |

**Columns of `se_tap_estimates.csv`.** Steps are the changer's own step
numbers. A state of the row is given in `status`, `reason` and `notice`,
never in a step column.

| Column | Meaning |
|---|---|
| `branch`, `regulator` | branch index and index of the changer in the changer set |
| `name`, `changer` | transformer name and changer name |
| `mrid`, `changer_mrid` | mRID of the transformer and of the changer (CGMES) |
| `step_transferred` | step of the input (SSH profile, case file, measurement file); the middle position where the import carried no step |
| `step_estimated` | step after rounding the estimate, the step of the fixed run; the transferred step where the row was not estimated or is `out_of_range` |
| `step_estimated_raw` | continuous estimate before rounding (neutral step plus `s`); the start position where the row was not estimated |
| `kind`, `side`, `alpha_deg`, `du` | changer data from the import |
| `du_est`, `ratio_est`, `angle_est_deg` | Delta U = `s du`, branch ratio and angle at the estimate |
| `deviation_steps` | `step_estimated` minus `step_transferred` |
| `sigma_s`, `round_distance_sigma` | standard deviation of the step, distance of the raw estimate from its nearest step in units of it |
| `step_from_target` | step from the control target of a regulator kept out automatically, not applied |
| `fixed` | the step was fixed for the final run |
| `reactance_source` | impedance source: `table`, `analytic`, `nominal` |
| `out_of_range` | estimate more than half a step outside the band; with status `estimated` it was fixed to the band end step |
| `status` | `estimated`, `out_of_range` (too far outside the band, transferred step kept), `not_estimated`, `calculated` (machine transformer) |
| `reason`, `note` | why a row was not estimated, and the note of an estimated row |
| `notice` | the text of the row, for example a mesh without P measurement or the band check |

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
machine precision. The start position is the transferred step; where the
import carried no step, the changer starts at its middle position (the
neutral step, the band middle without one) and its report row says "step
not transferred, middle position assumed", also in the fallback result
and on a machine transformer. The Jacobian column of a step is a forward difference
along the analytic tangent `dt/ds` (with the impedance part when an
impedance source applies), with a difference step in `s` of `jac_eps` over
the size of one step, so `t` moves by about `jac_eps`. At most two tap
changers of one transformer can be told apart, because the complex ratio
has two real directions; a further changer, or one whose effect equals that
of another (two ratio changers, ratio and asymmetrical at alpha 0), is
not estimated, reason `not_separable`. Only declared tap changers are
released: a transformer without a tap changer model (a MATPOWER branch
with a ratio only and no `mpc.sparlectra.tap_changers` row, a DTF
transformer without changer data, a CGMES or PowSyBl transformer without
`TapChanger`) has an empty set, is never released and has no tap row.

**Fixation.** After convergence every estimated step is rounded to the
nearest physical step (ties away from zero; not a regulator with note
`relative_to_frozen`, see below). The band of the changer is read
absolutely from its data (`lowStep` to `highStep`), not relative to the
import step; an estimate outside it is handled as below. A phase table fixes on its nearest point, a
`RatioTapChangerTable` on its points inside the table and on whole steps of
its formula outside. Then one run with the steps fixed follows.
`out_of_range` marks an estimate more than half a step outside the band.
Within 2 `sigma_s` of the band end step it is fixed to that end step and
stays `estimated`; farther out, or without a `sigma_s`, it is not taken
over: `status = out_of_range` (`reason = none`, `fixed = false`), the
transferred step stays. The notice names the band, the distance from the
band end in `sigma_s`, what the fixation did, and for an asymmetrical
changer whose `|t|` has a turning point the mirror step on the other side
of it; `runse!` logs one warning per such changer with the same text. The
same rule applies to a changer estimated again by the recheck below. A
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

**Observability per changer.** Before the first iteration a rank test on
the measurement Jacobian at the start point decides which steps are
estimated. Every column (voltage angles and magnitudes,
shunt and offset states, the steps) is scaled to norm 1; a step column
that adds no rank to the state columns and to the step columns before it
is not estimated: it keeps its transferred step (SSH profile or
measurement file) and its row says `status = not_estimated` with the
reason. The bound on the singular values is `1e-8` after the scaling, it
does not read the weights, so one heavy row (a pseudo measurement with a
small sigma) decides nothing. The test runs on the central difference of
the Jacobian (the forward difference has an error of about `jac_eps` times
the curvature, which hides an exact dependence). A current row (Imag, Ia)
counts only when its predicted current at the start passes three sigma of
its magnitude, the gate the estimator applies to the Ia rows: at zero
current the magnitude has a kink and the angle no value. At a flat start a
transformer at ratio 1 between equal voltages carries no current, so a
PMU current angle at its terminal does not anchor its phase changer there;
two Va rows of one time base do, and so does the Ia row from a warm start
(the staged start of the SE service). The step columns run in
branch order, then in the changer order of the transformer, the ones
without a structural cause (see Reasons) before the ones with. Of columns
that share one dependence the last one is not estimated: the changer a
structural rule names (the phase changer of a ratio and phase pair on a
bridge, the tree edge of parallel phase changers), else the last parallel
member. A column zero within the rounding noise of its forward
difference on every row is not estimated as well. No angle measurement
is needed: a phase changer in a mesh is estimated when the flows of the
mesh fix it; one on a bridge with a far-side voltage measurement is
estimated when it moves the magnitude at its start (an asymmetrical
changer at 60 degrees).

**Recheck at the estimate.** The rank test judges at the start point, and
from another start it can give another verdict. After the fixation the
noise and rank test run once more at the converged state with the fixed
steps, over every fixed changer and every changer the start test kept out.
A fixed changer that adds no rank there is withdrawn
(`not_observable_at_estimate`); a kept-out changer that adds rank is
estimated once more together with the other free changers of its
transformer, then fixed. A changer estimated this way is taken over only
when its `sigma_s` is below 0.5 steps; otherwise the changer keeps its
transferred step (`weak_at_estimate`), the final state is solved again
with that step, the row keeps the estimate in `s_est` and `sigma_s`, and
the re-estimate runs again without it while another newly determinable
changer is left. The bound applies to this
re-estimate only, not to the first estimate. A fixed phase
changer without an angle anchor that has a mirror step is solved once at
that step (`mirror_step_fits`). The angle anchors are judged at the
converged state as well: an Ia row whose current is zero at a flat start
anchors there. A run in which the start test keeps every changer out gets
the same recheck on its converged state. When the recheck withdraws a
changer, or the band verdict keeps one at its transferred step
(`out_of_range`), the other changers of its transformer still carry steps
from the joint estimate with it. They are estimated once more with it at
its transferred step, then the band verdict, the fixation and the mirror
check apply as before, and the final state is solved with the reported
steps. If that run does not converge, they keep their transferred step
(`not_converged_at_estimate`). Each withdrawal and re-estimate writes an
info line.

**Partner not determinable.** When the rank test or the recheck keeps one
changer of a transformer out (`radial_no_voltage_pin`, `not_separable`,
`no_angle_anchor`, `symmetric_start`, `not_observable`,
`not_observable_at_estimate`, `weak_at_estimate`, `mirror_step_fits`,
`not_converged_at_estimate`, see Reasons) or the band verdict does not
take its estimate over (`status = out_of_range`), and another changer of the
same transformer is estimated, the estimated row keeps
`status = estimated` and its notice names the partner, its reason and its
transferred step: the estimate holds only together with that step. A ratio
changer next to a phase changer without an angle anchor is the usual case.

**Reasons** of a changer with `status = not_estimated` (`reason`, with a
`notice` text per row and a `notice:` line in the run log of the SE
service). The structural rules name the cause of a dependence the rank
test found; on their own they freeze nothing:

- `radial_no_voltage_pin`: a changer that moves the magnitude of the ratio
  (ratio, asymmetrical, a table of varying ratio) on a bridge transformer
  (its removal cuts the net) whose cut-off side carries no active voltage
  measurement: the step and the voltage behind it move together. A pure
  angle changer there moves no magnitude; it gets `no_angle_anchor` or
  `symmetric_start`.
- `not_separable`: see above.
- `no_angle_anchor`: a phase changer on a tree edge of the angle graph
  (its far side is tied to the reference by no Va measurement, no measured
  end of a counted Ia row and no meshed path) on a transformer that also has a
  ratio changer, which covers what the phase changer moves of the
  magnitude.
- `symmetric_start`: a phase changer on such a tree edge without a ratio
  changer, at a start where `|t|` does not move (symmetrical or linear,
  an asymmetrical changer at its turning point, for example 90 degrees with
  `s0 = 0`, a table on a flat piece of its ratio). The notice says whether
  the changer moves the magnitude away from that step. At a vertex of a
  table ratio the step column carries the slope to the right of the start,
  so the rank test keeps the changer; it stays on the side of its start
  (`sign_from_start`) and the vertex notice watches it.
- `not_observable`: the rank test drops the column and no structural cause
  applies, or (notice "the measurements do not see this regulator") the
  column is zero within the rounding noise of its forward difference on
  every row (per row about `10 eps` times the size of the terms of the
  row, taken from its voltage derivatives, times the step size over
  `jac_eps`; a zero injection has terms of thousands of MW).
- `rank_budget`: the number of released changers puts the rank test at
  the start point above its work budget (`m k^2 + k^4/11` above `1e11` for
  `m` rows and `k` changers, about one minute); no changer of the run is
  estimated and the run is the run without tap states. It stays as a
  safeguard; it was hit on case6495rte while every transformer with a
  ratio was still released (8617 changers, the test would have run for
  days).
- `not_observable_at_estimate`: the changer was estimated and fixed, but
  the rank test at the converged state with the fixed steps drops its
  column (see the recheck above); the estimate is withdrawn and the
  transferred step stays.
- `weak_at_estimate`: the changer was kept out at the start point and
  became determinable at the converged state, but its re-estimate there
  has `sigma_s` of 0.5 steps or more. The estimate is not taken over, the
  transferred step stays (the state is solved with it), and `s_est`,
  `sigma_s` and the notice give the estimate and its `sigma_s`. The case this catches: a bridge transformer without a voltage
  measurement on its cut-off side, whose step gets rank only from the
  losses of its load.
- `mirror_step_fits`: a symmetrical or asymmetrical phase changer on a tree
  edge without an angle anchor whose mirror step (the sign of the step
  reversed about the neutral step, or about the turning point of `|t|`)
  fits the measurements as well as its fixed step: J of the fixed run at
  the mirror step differs by less than 9 (3 sigma of chi2 with one degree
  of freedom). The estimate is withdrawn, the transferred step stays, the
  notice gives both J.
- `not_converged_at_estimate`: a partner of a withdrawn changer was
  estimated again with that changer at its transferred step (see the
  recheck above), and that run did not converge; the transferred step
  stays.
- `model_mismatch`: the branch ratio and angle do not lie on the curves of
  its changer data. The release is recorded, the regulator is frozen.
- `disabled`: `state_estimation.tap_estimation.enabled = false`.
- `fallback_not_converged`: the run with tap estimation did not converge
  (a set too thin for the taps, or a gross error that a regulator pinned
  only by its magnitude effect tries to absorb); the result is the run
  without tap estimation, see
  [Tap estimation switch and fallback](@ref se-tap-fallback).
- `fallback_worse_fit`: the SE service rejected the fixed estimate because
  it fits worse than the staged pre-solve with the transferred steps
  (`worse_fit`); the result is the run without tap estimation.
- `no_changer`: a released transformer without a changer of the released
  role gets one report row (regulator 0), which the counts of the SE
  service leave out.
- `excluded`: the transformer is named in
  `state_estimation.tap_estimation.exclude`. The transferred step stays
  and is not calculated.
- `machine_transformer`: a machine transformer kept out by default
  (`include_machine_transformers = false`). After a converged run with
  `update_net = true` the row of its first ratio changer becomes
  `status = calculated` with this reason. Without `update_net`, when the
  calculation fails, and for any further ratio changer, the row stays
  `not_estimated`. The notice says why.
- `q_not_measured_in_group`: a machine transformer kept out by default
  whose machine node (the buses closed links tie together) has no Q
  measurement on every bus; its step is not calculated, the transferred
  step stays, and an info line names the buses without Q.

A changer the rank test kept out names in its notice the measured minus
estimated voltage at the terminal of its transformer with the larger
deviation (when a terminal has a Vm row): a wrong transferred step shows
there. For the reasons of the rank test and the two fallback reasons the
row also gets `step_from_target`, after a converged run with
`update_net = true`: the step that meets the control target of the
enabled tap controller of the transformer, rounded to the nearest step.
For a ratio changer with voltage control the voltage at the other
terminal and the power through the regulated terminal stay at the
estimate; for a phase changer with active-power control over the
transformer both terminal voltages stay, so the step is a lower bound of
the move. The step is not applied (the target is an assumption: a
controller at a limit, in hand, or with another target); when it differs
from the transferred step by more than one step the notice says so. The
column stays empty without a control target, for excluded and machine
transformers, and for `disabled`.

**Step sigma.** Every estimated step carries `sigma_s`, the square root of
its diagonal entry of the covariance of the continuous estimate
(`(H' W H)^-1` at the last iteration, with zero-injection constraints the
block of the KKT inverse), and `round_distance_sigma`, the distance of
`s_est` to the step it rounds to in units of `sigma_s`. Both are hints: a
large distance says the estimate lies far from a whole step for its
uncertainty (a true position between steps, or wrong step data), a
`sigma_s` of half a step or more says the rounding itself is uncertain.
The rounding does not read them. The SE service writes both to `se_tap_estimates.csv` and
to the regulator line in run.log.

Notes (`note`) that do not freeze: `p_from_schedule` on a calculated
machine transformer whose machine node has no P measurement on every bus
(the calculation took the scheduled P), `relative_only` for parallel members
(estimated only relative to the other member), `conditional` on the
ratio changer of a transformer whose phase changer is not estimated
(`no_angle_anchor`) while it moves `|t|` (asymmetrical, or a table whose
ratio varies; not for symmetrical or linear), `sign_from_start` for an
estimated phase changer on a tree edge without an anchor whose `|t|` has
a turning point on its band (it keeps the side of its start: a step that
would cross the turning point goes half the way and stays on the turning
point once it is that close), `relative_to_frozen` for a phase regulator
whose loop closes only over a phase regulator that is not estimated (for
any reason, for example `no_angle_anchor` or `model_mismatch`): it is only
determined relative to that one, so it is not fixed (`fixed = false`),
stays a continuous state in the fixed run (its `s_fixed` is that
continuous value, rounding it would miss the measured loop angle), its
notice names the other transformer, `k_removed` does not count it, and
`update_taps` keeps its model position. Notices of the estimate: a phase
changer on a tree edge without an anchor whose estimate sits near the
vertex (only the magnitude effect pins it, its step is not reliable); a
CGMES changer fixed at a step whose factor leaves 0.5 to 2 (a re-import
ignores it there); an estimated column more than `1/jac_eps` times weaker
than the strongest voltage column of the unweighted Jacobian gets "weak
column, factor x below the strongest", also in the run log. A changer
that is not estimated solves bitwise identical to no release and still
appears in `tapEstimates` with its reason. Release and back-calculation
are mutually exclusive per transformer.

**Staged start and fallback in the SE service.** The Web UI's "estimate
taps" first solves once with the taps frozen and then estimates them from
that state (from a flat start many tap states crawl to the solution). When
the run with tap estimation does not converge, or its fixed estimate fits
worse than that pre-solve (J above the pre-solve J by more than 1e-6
relative, `worse_fit`), the service repeats the estimation without tap
estimation, and that run is the result of the round: two lines in run.log,
`se_tap_fallback_used = true` in the metadata, reason
`fallback_not_converged` (or `fallback_worse_fit`) per changer. The repeat
starts from the start state and configuration of the round, not from the
pre-solve, so its result is the run with tap estimation switched off.
When the pre-solve itself does not converge (also not on the second try of
the start ladder from the case voltages), no released run is started: the
result is the run without tap estimation (`se_tap_fallback_reason =
staged_not_converged`, rows `fallback_not_converged`); with
`tap_estimation.fallback = false` the pre-solve is the result and the run
fails. The bad-data diagnostics then run on the end steps of that result; an
eliminated row starts the next round with released taps on the cleaned
set, see [Tap estimation switch and fallback](@ref se-tap-fallback).

**Off-grid positions.** A band failure that appears only through the
fixation is not bad data: when the continuous fit passes or undershoots the
band and the fixed run fails it `:high`, `tapFixation.offgrid_residual` is
set (summary note `:offgrid_tap_residual`). The true position then sits
between steps, or the changer data (step size, neutral step) are wrong.

**Impedance per step.** The rule is per transformer. Per-step data (the
r, x, g, b of a tap table row, the reactance curve of a phase tap changer
with `xMin`/`xMax`, the factors of a `TapStepImpedance`) always follow the
step, and the analytic correction is then off for that transformer. A
typed model (an SCF `tap_changer_models` entry, `addPIModelTrafo!`)
follows the power flow instead: its `x_min`/`x_max` curve or x table
follows the step under both values, and the correction factor on r that
the power flow applies next to the curve under
`model.tap_changer_model = impedance_correction` is taken under `auto`
only (source `table` under both).
`state_estimation.tap_estimation.impedance_correction` (`auto`, `off`;
default `auto`) decides for a transformer without per-step data: `auto`
applies the impedance correction of typed models under
`model.tap_changer_model = impedance_correction`, `off` keeps the
impedance of the start position. There is no `on`: without per-step data
or a typed correction no curve exists to force. Each row names the source
(`reactance_source`: `table`, `analytic`, `nominal`). With a source the
Jacobian column of the regulator includes the part through the impedance.

**Fixation and dof.** The fixation raises `dof`: every estimated step is
one extra state, `dof_before = m - (voltage states + steps)`, and the
fixation removes them, `dof_after = m - voltage states = dof_before + k`.
`J` moves the other way, e.g. `J 58.5 (dof 41) -> J 62.5 (dof 42)`, both
inside their band; only a fixation that jumps out of the band is a finding.
On multi-island nets the sums aggregate like the island chi-squares.

**Machine transformers.** The AVR holds the machine terminal voltage and
nothing observes it on its own, so a released tap would only absorb the
voltage control. `runse!` keeps such transformers out of the estimation,
also when released with `setTapEstimation!`, unless
`state_estimation.tap_estimation.include_machine_transformers` is set. A
machine transformer has a terminal with a generator and no other closed
branch. The buses that closed links tie together count as one terminal.
After the run the step is calculated from the estimated network-side
voltage, the AVR setpoint and the measured P and Q of the machine node.
Without a Q measurement on every bus of the node it is not calculated
(`q_not_measured_in_group`). Without a P measurement on every bus the
scheduled P is used (`note = p_from_schedule`).

**Guards.** Tap states are outside the observability count (like the gated
current rows); whether a step is observable decides the rank test of
"Observability per changer" above, the structural rules only name the
cause. A changer that is not estimated keeps its transferred step and its
exact model value.

**Known limit: mesh without active power.** A mesh without an active power
measurement can have a second converged solution with turned-over angle
differences. Q, |I| and Vm rows fix the size of an angle difference better
than its sign. On sp_casePST with no P flow, no |I| and no P injection in
the mesh 2-3-4-5-7 (the P row of bus 9 is dropped because its link partner
bus 2 is not measured), the WLS from the flat start converges with the
angle differences on br1 and br3 turned over, J 15.4 at dof 33 on exact
data, while the true state gives J near zero; br8 is then fixed one step
off. The topology precheck warns in this case
(`:mesh_without_p_measurement`, see Topology validation below). On
MicroGrid BE the same kind of set fires the warning as well: without tap
states the flat start still finds the true state there, with the four
changers released br8 and br12 end on wrong steps. A P row in the mesh
clears the warning but does not prove uniqueness: on sp_casePST a P
injection at bus 3 alone still gives a second solution (J 16.2), one at
bus 2 gives the true state. A different start value only changes which
solution is found, so the remedy is a P measurement in the mesh. A released
changer on a branch of such a mesh keeps its status and reason, its notice ends
with "lies in a mesh without P measurement: the step may not be unique"
(only with `state_estimation.topology_precheck = true`, the default).

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
status-contradiction check. One structural finding concerns uniqueness, not
a topology error: `:mesh_without_p_measurement` names the branches of a
mesh with no active power row on its branches or buses ("mesh without P
measurement: solution may not be unique"); `runse!` warns with its own log
line and the SE service lists it with the other findings. A P injection on
a link cluster counts only when every member is measured, as in the
estimator, and parallel branches count as one connection. See the known
limit under Transformer tap estimation.

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

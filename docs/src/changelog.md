# Version 0.9.16 - 2026-08-23

CGMES deliveries that arrive without their topology file, and a UPFC that
finally steers a line's active and reactive flow at the same time.

## Highlights

**Node-breaker CGMES imports without a TP profile.** EMS and substation
exports often ship only EQ and SSH, leaving the bus partition implicit in
the connectivity nodes and switch states. Sparlectra now works it out
itself: connectivity nodes merge across closed switches (an open switch in
the SSH wins over the equipment default, out-of-service counts as open),
retained switches stay as bus couplers, and boundary nodes hook into the
boundary set. Deliveries that do ship a TP are untouched, and on the
conformity sets the derived partition matches the sender's bus for bus. A
delivery with neither a TP nor connectivity nodes now stops with the import
analysis rather than a puzzling error further down. See
[CGMES Import](cgmes_import.md).

**The UPFC, both ways.** `addUpfcControl!` gives you the unified power flow
controller. The default `model = :quadrature` pairs an SSSC and a STATCOM as
one named device; `model = :full` (#326) lifts the quadrature restriction so
the series converter injects a voltage of arbitrary phase and the line holds
independent P and Q at once, with the active part balanced across the DC
link. Forcing the series phase back to quadrature collapses the full model
onto the composite. Honest about the first cut: the shunt follows a reactive
setpoint (closed-loop voltage regulation is still to come) and convergence
is dependable for realistic, moderate targets. See [FACTS Devices](facts.md).

**A CGMES workshop tour.** A new expert-level Colab tour walks the official
ENTSO-E conformity sets end to end: what a delivery is, reading a broken one
with the import analysis, bus-branch import with SV validation, node-breaker
with and without a TP, and the export round trip.

## The rest

- **The classical result knows about FACTS now.** Every outer-loop
  controller (STATCOM, SVC, TCSC/SSSC, HVDC pair, UPFC) shows up in the
  power-flow result: counted in the header, marked on its branch or bus, and
  summarised in the control footer. The FACTS tour prints these tables for
  each device.
- **Ferranti voltage in the bus table.** An isolated open-end bus shows the
  open-end voltage of its one-sided-open branch (flagged `open-end`) instead
  of a meaningless start value; the basic tour contrasts the correct pi-stub
  model against a full disconnect, which quietly drops the charging draw.
- **Smaller things.** `analyzeCGMES` now loads with `using Sparlectra`;
  `print(summarizeCGMES(...))` renders the readable multi-line report
  instead of a one-line dump; and every workshop example carries a unique
  number, a line diagram for each network it introduces, and a note on its
  AI-assisted origin.

# Version 0.9.15 - 2026-08-23

Fix and polish, and observability theory you can execute.

- **Dark states, named.** A not-observable set now reports
  `unobservable_state_columns`: the states the null space of the
  measurement Jacobian touches, i.e. exactly what no meter pins down. The
  column-restricted local check documents its limitation (necessary, not
  sufficient); the workshop shows the two disagreeing verdicts on the
  smallest possible counterexample.
- **The state-estimation workshop teaches the theory.** WLS foundation
  and normal equations, the null-space criterion with observable islands,
  the spanning-tree counting identity, criticality as a zero residual,
  J/dof as the after-run health number, matrices small enough to read,
  and a zero-injection section with the honest sigma trade-off.

- **The estimation problem shows its matrix.** New `measurement_jacobian(net)`
  returns the labeled Jacobian H (described measurement rows, named state
  columns); the state-estimation example suite writes a measurement-matrix
  page with a stability verdict (rank, redundancy, cond(H)) on every run.
  The workshop's observability deep dive now explains what global and
  local observability ARE and draws the meters into the line diagram.
- **Docs sweep.** Stale claims corrected across the pages (RealGrid has
  6209 buses, `control.controllers` is implemented, `hvdc_pair` is a
  supported type, chapter pointers follow the tour split, two release
  dates fixed); the drifted hand-written config key index now points at
  the canonical `configuration.yaml.example`; the transformer control
  theory lives once, on the Control Framework page.
- **Code sweep.** Dead code removed, including four dead public bindings
  (`addBusShuntPower!`, a deprecated no-op, plus `zero_row!`,
  `print_jacobian`, `recalc_trafo_model_data`) and six unused private
  helpers; the Y-bus assembly no longer rescans the isolated-bus list per
  branch, result tables build without quadratic string growth, and the
  reverse bus-name lookups are unified.

# Version 0.9.14 - 2026-08-22

Intentionally identical to 0.9.13: the registry needs the gapless number,
the project prefers not to stay on it. Install this one.

# Version 0.9.13 - 2026-08-22

Coordinated transformers, switched shunt banks, and observability you can
see.

- **Parallel transformers regulate as a group** (#322). Two independent
  tap controllers on one busbar fight each other through circulating
  reactive power, and CGMES deliveries produce exactly that setup whenever
  several tap changers share one `RegulatingControl`. Now the group is the
  model: the master runs the voltage loop, followers mirror its steps
  (`followers` on `addPowerTransformerControl!`), and the importer turns
  shared controls into one group instead of fighting controllers. Part 3
  of the workshop tour walks the whole story, from the circulating-current
  trap to the CGMES data.
- **MSC/MSR switched shunt banks** (#324). `step_mvar` turns the shunt
  voltage controller into a mechanically switched bank: whole blocks only,
  never overshooting the target, parking on the reached step instead of
  hunting between two.
- **Observability you can see.** The state-estimation workshop now shows
  what the observability numbers mean on a drawn network, where
  measurements must sit (and how a missing one darkens a corner), and how
  a PMU's time-base offset is estimated. Every workshop notebook states
  its level.

# Version 0.9.12 - 2026-08-22

The workshop splits in two, and every solver reads the same network.

## Behavior change

- **The DC power flow now reads the same specification as the AC solver**
  (issue #323). Both solvers build their injections from the prosumer
  objects; the node-level sums (`addBusLoadPower!`/`addBusGenPower!`) are
  a report layer that NO solver reads anymore. Before, a node-sum edit
  changed a DC solve but not an AC solve, so the two could silently
  disagree about what the network is. Affected: only workflows that
  edited node sums and then ran the DC solver; the fix is one line, add
  the load as a prosumer. Imported cases fill both layers identically and
  are not affected.

## Highlights

- **The workshop tour splits into basic and advanced.** The basic tour
  carries Newcomer to Advanced (first network to Q(U) control); the new
  [advanced tour](generated/workshop_tour_advanced.md) carries Expert and
  Beyond with room to breathe: remote voltage control, HVDC pairing,
  state estimation, the FACTS limit workshop (now with the in-range view
  and the generic element table), the N-1 workshop (now with worst-case
  ranking and CSV export), and the threaded sweeps.
- **Half-open branches show themselves in the classical result.** The
  summary prints one status line per one-sided open branch with the
  open-end voltage and its Ferranti rise over the feeding bus, the bus
  table marks the open-terminal bus (`open-end` in the Control column),
  and the branch row carries the open-end voltage next to its
  `open@to`/`open@from` marker.

## Improvements

- The Web UI start scripts launch Julia with `JULIA_NUM_THREADS=auto`
  (unless already set), so the threaded surfaces actually use the
  machine's cores; `runtime.parallel.*` was already on by default and
  sizes itself from the thread count.

# Version 0.9.11 - 2026-08-22

FACTS limit modes, one workshop notebook, and better discoverability.

## Highlights

- **STATCOM and SSSC limit modes.** A compensator shows its character at
  its limit, and that is usually during the voltage sag it was installed
  for. The machine voltage controller can now act as a STATCOM
  (`s_max_mva`): the bound is the converter current, so the deliverable Q
  follows the terminal voltage linearly instead of freezing at a fixed
  box, and instead of collapsing with V² like the SVC. The series
  controller can act as an SSSC (`v_inj_max_pu`): its reactance window is
  set by the injectable series voltage and narrows as the branch current
  grows. Both bounds follow the solved operating point from iteration to
  iteration. Why a UPFC is not on this list, and what to use instead,
  is explained on the new [FACTS Devices](facts.md) page. (#297)
- **The whole workshop in one notebook.** The introduction notebook and
  the workshop page are gone as separate things; the tour now carries
  everything in five tiers from Newcomer to Beyond, including two new
  Expert chapters on the FACTS limits and on N-1 contingency batches.

## Documentation

- Search engines and AI agents get proper signposts: `robots.txt`,
  `sitemap.xml`, `llms.txt`, meta description and preview image, and a
  page title that no longer collides with a bulk carrier of similar name.
- The feature matrix gives each analysis domain its own table. State
  estimation is no longer graded against power-flow features it never
  claimed to have, and every warning marker now says why it is there.

## Bugfixes

- The docs landing page showed broken image references where the README
  logo and the Web UI screenshots belong; the README embed now rewrites
  the image paths for the site.
- Re-running a configuration with a declared `series_reactance` controller
  crashed instead of recognizing the controller was already there.
- `addBusLoadPower!` and `addBusGenPower!` change sums the AC solver never
  reads, which made them silent no-ops for AC studies. Their docstrings
  now say so and point to the prosumer-based way; the underlying AC/DC
  inconsistency is #323.

# Version 0.9.10 - 2026-08-21

Multi-core execution, N-1 contingency analysis, and one-sided open branches.

## Highlights

- **Multi-core execution.** Independent work now runs on Julia threads
  (start with `julia --threads=auto`): AC islands
  (`power_flow.islands.mode: solve_parallel`), all-bus short-circuit
  sweeps, and contingency batches, gated by the new `runtime.parallel.*`
  keys (`enabled`, default true). Results stay identical to the serial
  runs, island for island and row for row. Measured on 16 cores: 4.8x on
  an 8000-bus fault sweep, 4.1x on a full N-1 over case1354pegase. The
  demos live in `examples/run_parallel_suite.jl`.
- **N-1 contingency analysis.** New `runContingencies!` batch API:
  branch-outage cases from `generateN1Branches` or imported FOR001 lists,
  solved warm-started on copies of the base case (never mutating it) and
  checked against the voltage band and `sn_MVA` loadings. Islanding and
  non-convergence are reported per case instead of thrown, parallel
  circuits get their own cases, and results print as a table or CSV. The
  solves run without the rescue ladder; `retry_flat_start = true` grants
  one flat retry per failed case. See
  [N-1 Contingency Analysis](contingency.md).
- **One-sided open branches.** Branches carry `from_status`/`to_status`
  next to the aggregate `status`: open at exactly one terminal, a branch
  stays in the model as its exact pi reduction, draws its full charging,
  and reports the open-end voltage (Ferranti rise) as a branch result.
  CGMES `Terminal.connected` maps to the flags on import and export
  (replacing the half-charging substitute shunt); MATPOWER export
  represents the state as the exact `Y_in` bus shunt with an
  `open_terminal=` marker.

## Performance

- **Takahashi sparse inverse** for all-bus fault sweeps, opt-in via
  `runShortCircuit!(net; sweep_method = :takahashi)`: the whole Thevenin
  diagonal of an island from one pass over the LU factors instead of one
  solve per bus, measured 34x to 264x over the serial sweep. Values agree
  with the default to machine precision; inapplicable islands fall back
  automatically. Theory in the
  [Short-Circuit Compendium](short_circuit.md).
- The solver status moved from global registries onto the `Net`, and
  per-island performance entries no longer overwrite each other. Two
  subtle serial-path details changed with the shared island solve: a
  failing island raises its error from the result handler (different
  backtrace origin, same error), and a nonfinite-voltage failure no
  longer counts its iterations into the aggregate.

## Notebooks

- New: distributed slack (where the participation weights come from and
  what happens without a valid participant) and transformer taps (OLTC vs
  PST device math, the phase-tap control loop, the 3WT star equivalent).
  Every notebook now opens with a topic warmup cell; the tour gained
  chapters on bus links, one-sided open branches, and parallel sweeps.

## Removed

- The `klu` linear-solver backend and the KLU.jl dependency: measured
  slower than UMFPACK on power-flow Jacobians, and a shared KLU
  factorization is unsafe under threads. Migration: replace
  `linear_solver: klu` with `linear_solver: umfpack_reuse`.

# Version 0.9.9 - 2026-08-20

HVDC link results and one-reference-per-island validation.

- HVDC links are first-class results: every link (Stage-0 or paired, from
  MATPOWER, CGMES, or the API) is tracked on the net and reported in a new
  `HVDC Link Flows` table, a `Link` row in the branch table (marked
  "HVDC, not a branch"), the result header, the generator-table terminal
  markers, a converter-loss line, `ACPFlowReport.hvdc_links`, and
  `hvdc_links.csv`.
- A synchronous island with more than one angle reference now fails with a
  message naming the island and its reference buses (both solver paths)
  instead of the generic unsupported-bus-type error; an `island_feed` pair
  whose grid-forming reference was demoted reports `invalid_topology`. See
  [HVDC Back-to-Back](hvdc_back_to_back.md), example
  `exp_hvdc_meshed_ac_tie.jl`, and the meshed-operation tour subsection.
- Branch `Type` column and report rows agree again: the printed label
  consults the import metadata first (PST heuristic pending the typed
  branch kind).

# Version 0.9.8 - 2026-08-19

**Steerable HVDC links.** `addHvdcPairControl!` couples the two converter
injections of a back-to-back link with exact loss accounting, per-terminal
Q or voltage targets, and honest `at_limit`. Both importers attach pairs
opt-in (`paired_control`); Stage-0 fixed injections stay the default. A
grid-forming `mode = :island_feed` feeds an island whose only source is
the receiving converter. See
[HVDC Back-to-Back](hvdc_back_to_back.md). (#297)

- The one-line install detects an older existing copy and offers the
  update (`SPARLECTRA_UPDATE=1/0`); the previous copy is kept as
  `Sparlectra.old`.
- Leaner, faster installs: two unused hard dependencies dropped (76 to 54
  packages), and the workshop notebooks now install into a fresh temporary
  environment from GitHub `main`, keeping Colab's preinstalled stack out
  of precompilation.
- The result header now counts every island reference and keeps the
  reference models apart (`Slack: n Source: m`) instead of always printing
  1; the workshop notebooks carry a network sketch per model.
- Each GitHub release now ships an SPDX SBOM (`Sparlectra-<tag>.spdx.json`)
  as a release asset, generated by CI via `tools/generate_sbom.jl`
  (PkgToSoftwareBOM, kept out of `Project.toml`).

# Version 0.9.7 - 2026-08-20

- Install scripts offer a desktop shortcut (Windows `.lnk`) or launcher
  (`.desktop` / symlink) for starting the Web UI;
  `SPARLECTRA_CREATE_SHORTCUT=1/0` for unattended installs.
- The Jacobian condition estimate is always part of the result log and the
  Web UI run overview (with verdict); the `output.condition_number` option
  is removed, a leftover YAML key is ignored.

# Version 0.9.6 - 2026-08-19
* #309 Fix Sysimage for windows v0.9.5
* #308 One-line Web UI install without a GitHub checkout; the install scripts now offer the optional fast-start sysimage build (y/N prompt, `SPARLECTRA_BUILD_SYSIMAGE=1/0` for unattended installs).
* Condition estimate now comes from the solver (lazy, computed only when reports ask); fixed `condestJacobian(net)` on nets with de-energized buses (DimensionMismatch), the estimate describes the active subsystem.

# Version 0.9.5 - 2026-08-18

## Highlights

**Fast start for the Web UI.** Optional PackageCompiler sysimage, built via
`tools/build_sysimage.jl` or the new Fast start page; the start scripts use
it automatically and skip the warm-up. See
[Fast Start (Sysimage)](fast_start.md).

**Controllers in configuration.** Outer-loop controllers can be declared
under `control.controllers`, one named entry per controller. See
[Control Framework](control_framework.md). (#305)

## Improvements

- Faster PowerFlow form renders (case-cache scan memoized) and an explicit
  list of saved case settings a newer configuration file replaced.
- MATPOWER case downloads retry transient failures; run paths tolerate a
  symlinked user root.


# Version 0.9.4 - 2026-08-14

## Highlights

**Series compensation (TCSC).** `addSeriesReactanceControl!` steers a line
branch onto an active-power target by moving `x_pu` within its range,
with honest `at_limit` when the target is out of reach. Shows up in the
result tables like the tap controllers. Theory in
[Series Compensation (TCSC)](series_compensation.md), example
`exp_tcsc_series_reactance_control.jl`, Colab notebook
[TCSC Flow Steering](generated/workshop_series_compensation.md). (#297)

## Improvements

- Distributed slack participation moved into the bus table (`dSl alpha`,
  `Pg eff MW`) instead of a separate block.


# Version 0.9.3 — 2026-08-12

## Improvements

**Distributed slack participation in the classical result print.** With
`power_flow.distributed_slack`, `printACPFlowResults` now lists each
participating generator with bus, alpha share, correction `dP`, and
scheduled vs. effective output. Also available in the solver status as
`distributed_slack_participation`, see
[Power-Flow Configuration](powerflow_configuration.md).

# Version 0.9.2 — 2026-08-11

## Highlights

**(Source) External grid element.** `addExternalGrid!` models the grid connection
as an IEC 60909-0 network feeder — ideal slack by default, non-ideal
source with `internal_impedance = true`. Its declared short-circuit power
feeds `runShortCircuit!` directly, so hand-built and MATPOWER networks
need no CGMES delivery for short-circuit studies. See
[Slack Bus and External Grid Sources](slack_vs_source.md). (#299)

## Improvements / Documentation

- **Try Sparlectra in your browser:** two Colab workshop notebooks:
  an introduction (since merged into the
  [workshop tour](generated/workshop_tour.md)) and
  [Slack Types and Short Circuit](generated/workshop_slack_short_circuit.md).
- New theory page:
  [Slack Bus and External Grid Sources](slack_vs_source.md). (#299)

# Version 0.9.1 — 2026-08-01

Fixes and refinements after the 0.9.0 release.

## Improvements

- Saved case settings can be reset.

- Rescue loader handles Q-limit-driven divergence.

- CGMES deliveries start from their own state by default.

- DC-line cases run again on old installations.


# Version 0.9.0 — 2026-07-31

IEC 60909 short-circuit stage (#277), CGMES export.

## Breaking Changes
 
- **Short-circuit calculation (IEC 60909-0).** `runShortCircuit!` computes Ik''max/min, Sk'' and i_p per fault bus from the CGMES short-circuit data. Where attributes were missing and defaults substituted (e.g. motors), the affected rows are flagged — a flagged Ik''max is a lower bound. Web UI gets a Short circuit button for CGMES cases that carry the data; writes `short_circuit_max.csv` / `short_circuit_min.csv`. (#277)
- **CGMES export.** `writeCGMESFiles` writes a complete 2.4.15 delivery (EQ/TP/SSH/SV, optionally zipped). Export → re-import solves to the same power flow and keeps the original mRIDs; new objects get deterministic uuid5 ids, renaming never changes an mRID. Export checkbox in the Web UI run form.
- **Auto slack** (`power_flow.auto_slack`, off by default): promotes the strongest injection to slack when the case registers none, logged.
- **Rescue ladder + DC fallback** (`power_flow.rescue`, `power_flow.dc.fallback`, both off by default): failed AC solves are retried with alternate start / autodamp / DC-seeded start; if all fails, the DC fallback at least leaves angles and branch P in the net. AC status stays non-converged, honestly.
- **Import analysis.** Deliveries that can't be imported now explain themselves: missing boundary sets are named by model id, unresolved references get a histogram, verdict in plain language. Runs automatically on Web UI upload (`<case>.import_analysis.txt`), also via `analyzeCGMES`.
- **BaseVoltage reconstruction** (`cgmes_import.infer_base_voltages`, off by default): deliveries missing their BaseVoltage catalog get levels seeded from SV voltages and transformer rated voltages, then propagated. One summary warning.

## Improved

- `compareWithSV` removes the median angle offset before judging deltas — an IGM cut from a CGM keeps the global angle reference (a real 50Hertz delivery: ~27° on every bus, meaning nothing). Offset is reported separately.
- Start scripts `start_webui.sh/.bat` and `install_webui.sh/.bat` in the repo root.
- SVC-style voltage control via `addShuntVoltageControl!`; at the MVAr limit the device clamps to constant B (Q then follows V²). `controllableElements(net)` lists all controllers in one vocabulary.
- PSTs update X(α) on tap moves when the winding carries reactance data. Fixed along the way: DTF importer dropped the phase-tap model, and the phase-probe compared stale flows against themselves.
- Config edits show up on form refresh; newest edit wins between YAML and saved case settings.
- `cgmes_import.placeholder_guards: strict` aborts on filler values instead of skipping them.
- Diagnose / self-check evaluates CGMES cases at their own SV state — all start-value machinery forced off, including a `flatstart: true` that previously wiped the SV start silently. New artifacts `self_check.log`, `self_check_residuals.csv`.
- CGMES runs: start values selectable (`flat`/`sv`), SV comparison now mandatory on every run (`sv_compare.csv`, `sv_compare_flows.csv`).
- Wrong-branch detection only judges the highest voltage level — healthy 45 kV feeders no longer flag SUSPECT while the 380 kV level is clean.

## Fixed

- `no slack bus found` now says what's actually wrong: all buses isolated (branch-less delivery), slack on an isolated bus, or genuinely no slack — and how to fix each.
- Failed multi-island runs no longer report `iterations=0 / before_nr` when the island did a full Newton solve, and never-attempted islands no longer copy the failed island's statistics into `ac_island_solver_summary.csv`.


# Version 0.8.16 — 2026-07-30

CGMES Import

## New Features

**CGMES import.** Sparlectra reads ENTSO-E CGMES 2.4.15 and 3.0 deliveries —
folders, ZIPs, nested ZIPs, boundary sets — and builds a bus-branch `Net`
from EQ+SSH+TP+SV: lines, 2- and 3-winding transformers, tap controllers,
switches, shunts, loads, machines and injections, with snapshot-faithful
handling of open terminals and out-of-service equipment. Each electrical
island gets its SV-declared angle reference; assembled multi-area deliveries
are supported, with HVDC border crossings kept as separate islands.
`summarizeCGMES` inspects a delivery, `compareWithSV` validates the solved
result against its SV profile. Configuration under `cgmes_import.*`; details
in [CGMES Import](cgmes_import.md).

**Distributed active-power slack.** An island's power imbalance can be shared
over the generators by participation factors (one `lambda_P` per island in
the Newton state) instead of landing entirely on the reference bus. Enable
with `power_flow.distributed_slack.enabled` (default `false`; disabled runs
stay bit-identical). Weight modes include the imported MATPOWER `APF` /
CGMES `normalPF` factors. See the [Solver Guide](solver.md).

**PMU angles in state estimation.** New `VaMeas` measurement type and
`addPmuPhasorMeasurement!` for the Vm+Va pair; the estimator solves the PMU
time-base offset as an extra state (`state_estimation.pmu_ref_offset`) and
reports it in `SEResult.vaRefOffsetDeg`.

**Tabular phase shifters land in the solved branch.** CGMES
`PhaseTapChangerTabular` was read but never applied. The table row at the
tap position now folds into the imported branch (ratio and angle; an end-2
angle enters negated per end-referral). On RealGrid this reduces the
deviation from the delivered SV state by roughly 40×.


**Remote voltage control by machines.** A machine can regulate the voltage
at a *different* bus via its reactive output — the counterpart of a CGMES
voltage `RegulatingControl` at a foreign terminal, previously held PV at the
machine's own bus. `addMachineVoltageControl!` programmatically, or
`cgmes_import.machine_control` (default `false`) from a delivery; parks
honestly `at_limit`. Theory in
[Remote Voltage Control](remote_voltage_control.md).

**Nonlinear shunt compensators are mapped.** Per-section points are summed
up to the active section count (an interpretation choice, documented — the
place to look if a delivery disagrees). Previously their reactive
contribution was missing from the solve.


## Improvements

- The rectangular Newton solver refreshes its Jacobian in place; a PV↔PQ
  switch triggers a structural rebuild automatically.

## Bugfixes

- **Active-set Q-limits now switch on the converged iterate.** The start
  gating could outlast convergence, rejecting runs with `remaining PV
  Q-limit violations` instead of enforcing them; a converged iterate now
  always counts as ready.


# Version 0.8.15 — 2026-07-22

DC power flow, faster linear solves, Web UI overhaul

## Breaking Changes

* **Standalone DC power flow** (`rundcpf!`, `power_flow.solver = :dc`) with
  per-island handling, optional DC-seeded AC starts, and a calculation-model
  selector in the Web UI.
* **Factorization-reuse linear solvers** for the rectangular Newton step
  (`power_flow.linear_solver`: `umfpack` | `klu` | `umfpack_reuse`, default
  unchanged). Symbolic analysis runs once per sparsity pattern, only numeric
  refactorization per iteration. On `case_SyntheticUSA`, `umfpack_reuse`
  cuts the linear-solve phase from 1.66s to 1.31s; KLU is honestly slower on
  power-flow Jacobians. New dependency: KLU.jl.
* **Diagnostic report** instead of a flat key/value dump, plus a
  fixed-reference self-check — both available via a dedicated "Diagnose"
  action in the Web UI.

## Improvements

* Cleaner PowerFlow form: collapsible option boxes, unified solver
  selection (AC / APSLF / DC), all texts in English, contextual help via
  docs-backed "?" links.
* Case selector is now a single editable combobox: type to filter, Enter
  downloads unknown cases, right-click deletes cached ones.
* Warm-up solve is deferred and the form render path pre-compiled — the
  first browser tab no longer stays blank.

## Bugfixes

* **Tap controller: probe solves now follow the configured start strategy.**
  `control_propose_update!` ran its ratio/phase probe power flows with the
  hard-coded default `opt_flatstart = true` instead of
  `pf_config.start_mode.flatstart`. Imported networks with off-nominal
  ratio branches (MATPOWER, DTF, CGMES) routinely fail from a flat start
  while converging from their stored voltage profile, so any outer-loop tap
  control on such a network aborted with `pf_failed` — with no hint that
  the probe start values were the cause.
* MATPOWER parser accepts newline-separated matrix rows (RTS-GMLC style);
  validated against the MATPOWER 8.0 RTS-GMLC reference.
* New Monte-Carlo examples: probabilistic power flow (case14) and a WLS
  state-estimation error study (7-bus workshop net).
* Fixed two silent Web UI form bugs where disabled submit elements dropped
  their values (Diagnose button and DC selection ran normal AC solves).
* Test suite: 1444 fast-profile and 1835 Web UI extended tests green.
  
# Version 0.8.14 — 2026-07-22

## Comment
* intentionally the same as 0.8.13

# Version 0.8.13 — 2026-07-22

## Breaking News
* AnalyticLoadFlow.jl (APSLF) is now integrated as an optional extension.
* **APSLF solver**: can run standalone, as the framework solver, or as a
  start-value generator ahead of the rectangular NR solve. Doesn't yet
  cover OLTC/PST/Q(U)/P(U), only simple Q-limit switching — details in
  `docs/src/external_solvers.md`, demo in `examples/powerflow/apslf_demo.jl`.
* **Merit-function line search** in the rectangular solver: optional
  Armijo line search inside the autodamp backtracking loop, off by
  default. Requires `autodamp = true`. Adds diagnostic fields and a
  `merit_linesearch.log`. Background in `docs/src/solver.md`.
* **Trust-region step control**: alternative to `autodamp`, caps the
  Newton step adaptively instead of relying on the mismatch criterion.
  Two modes (`scaled`/`dogleg`); `dogleg` helps when the Newton direction
  weakens as the radius shrinks, but it's not a fix for bad starting
  values. Mutually exclusive with `autodamp`.
* Web UI: solver dropdown including APSLF, status view shows which
  solver actually produced the result.

## Improvements
* Wrong-branch detection results are now visible everywhere (report, CSV,
  log, Web UI, API) instead of only showing up as console warnings. An
  automatic rescue-retry is deliberately out of scope — detection plus
  APSLF as a fallback are considered enough for the hard cases.

## Bug fixes
* The classical Q-limit outer loop was losing data per round because the
  inner solver treats every call as a standalone solve: the Q-limit log
  got reset each round (undercounting PV/Slack→PQ switches, often showing
  `0` even when switching actually happened), and diagnostic files from
  one round overwrote the previous round's. Fixed by keeping the log
  across rounds and giving each round its own file names.
* Independently solved AC islands were overwriting each other's
  diagnostic files (merit, trust-region, iteration-start logs) since they
  all shared one output directory. Each island now gets its own prefix.

# Version 0.8.12 — 2026-07-20

## New Features

* Added an `AbstractTapChangerModel` supertype and an explicit `PowerTransformerTaps.convention` field (currently only `:neutral_relative`), documenting the ratio-tap correction convention on the struct itself instead of leaving it implicit. See `docs/src/branchmodel.md` for the tap-changer/PST modeling layering (Issue #261 Stage 2).
* Added a CGMES-style `PhaseTapChangerModel <: AbstractTapChangerModel` (`:symmetrical`/`:asymmetrical`, `:tabular` staged for later) and three pure formula functions in `equicircuit.jl`: `calcPhaseTapFraction`, `calcPhaseTapAngleRatio` (effective ratio/shift/regulating vector per CGMES v2.4 ch. 4.2/6.2), and `calcPhaseTapReactance` (CGMES ch. 3 `X(α)` interpolation, standalone for now — not yet wired into branch `x_pu` or the outer control loop). `PowerTransformerWinding` gained a parallel `phase_taps` field (Issue #261 Stage 3).
* Added `kind = :tabular` support to `PhaseTapChangerModel`, backed by a new `TapTablePoint` struct (`step`, `ratio`, `angle_deg`, optional `x_pu`) and a `calcPhaseTapTable` exact lookup in `equicircuit.jl`. A table overrides formula-based reconstruction whenever present, integrated into `calcPhaseTapAngleRatio`/`calcPhaseTapFraction`/`calcPhaseTapReactance`. No interpolation between table steps; no importer produces tabular data yet (Issue #261 Stage 4).
* Added `phase_tap_side::Int` / `phase_taps::Union{Nothing,PhaseTapChangerModel}` keywords to `create3WTWindings!`, letting one winding of a three-winding transformer (MVA-method, star/AUX-bus model) carry a `PhaseTapChangerModel` — i.e. a Schrägregler on a single 3WT winding, optionally combined with a ratio tap on the same winding. Validated: `phase_tap_side ∈ 0:3`, and `phase_tap_side`/`phase_taps` must be set together. Resolving the attached model into an effective branch ratio/shift, and addressing a single 3WT winding from the outer-loop `PowerTransformerControl` framework, are analysed but intentionally not implemented — see `docs/dev/3wt_phase_tap_controller_addressing.md` (Issue #261).

### Improvements

* Consolidated the ratio-tap correction and tap-range formulas into `equicircuit.jl` (`calcRatioTapCorrection`, `calcRatioTapRange`), the single source of truth used by both `calcTransformerRatio` (`transformer.jl`) and the `Branch` tap-limit derivation (`branch.jl`). Removed the duplicated inline formulas from both call sites; no change in computed values.
* The native DTF importer's skew-angle tap computation (`_dtf_effective_transformer_tap`) now constructs a `PhaseTapChangerModel` and calls `calcPhaseTapAngleRatio` instead of computing `tap_fraction` and calling `calcSkewAngleTap` inline. Reproduces prior behavior exactly, including the pure-longitudinal case (Issue #261 Stage 3).

### Bug Fixes

* `create3WTWindings!` previously raised a `MethodError` for every call, including its own docstring example: the `PowerTransformerWinding(...)` positional call had drifted out of sync with the struct's field order (a `ratio` field was inserted ahead of `shift_degree` without updating this call site), so every argument from `shift_degree` onward landed one field too early. Fixed by inserting the missing `ratio` slot; the pre-existing `tap_side` side-selection logic itself (documented as 1-based `[1,2,3]`, `0` = no tap) was left unchanged — see `docs/dev/3wt_phase_tap_controller_addressing.md` for the remaining discrepancy between that documentation and the current implementation.

# Version 0.8.11 — 2026-07-17

## New Features

* `writeMatpowerCasefile` now has a `matpower_export.write_solution` option (default `true`). It writes the solved power flow back into the MATPOWER case — bus VM/VA get the solved values, branches gain result columns 14–17 (PF, QF, PT, QT), pulled from the existing report/loss path, no extra recomputation. A `mpc.sparlectra.solution_written = 1` marker flags that the case has a solution. Turned off, you get the old pure model file: 13 branch columns, VM = 1.0, VA = 0.0. Ask for a solution before solving, and it falls back to the 13-column output with a warning. Configurable via YAML, the config API, or the new Web UI option.

### Improvements

* The tap-impedance correction factor from `calcTapCorrectedRX` (MATPOWER and native DTF import) is now saved on the Branch/transformer metadata. Reimporting a Sparlectra-exported case with the `tap_changer_model = 'impedance_correction'` marker won't stack a second correction on top. Cases without the marker — including third-party ones — work as before.
* Shortened the "Tolerance" and "Export Solution" labels so the help icon isn't crammed against the text. The tolerance exponent stepper now looks and behaves like a real number spinbox.
* Removed the inline explanation text next to "Export Solution" and "Tap-changer model" — it was overlapping the help icon. Full explanations are still one click away; docs were expanded to cover it.
* Fixed the "Autodamping minimum" spinner step: was 1.0 (way too coarse for a 0–1 range), now 0.01.
* Removed the remaining inline hints across the PowerFlow form (Tolerance, CSV format, benchmarks, etc.) — help text belongs in the docs, not hardcoded in the views. Docs updated accordingly.
* Combined one-off Web UI messages (errors, import results, settings-loaded notices) into a single popup dialog instead of scattering them inline. The configuration notice stays inline since it has a link you need to see. Error history page unchanged.

### Bug Fixes

* Fixed wrong `Case` name in run.log/result files for MATPOWER cases with bus_name metadata — a loop variable was overwriting the case name with the last bus's name.
* Fixed misaligned PV→PQ header lines in the report (labels too long for the column).
* Fixed misaligned Solver time/Total time headers (3 characters too narrow).

# Version 0.8.10 — 2026-07-17

## New Features

* Added a configurable transformer tap-changer model (`transformer.tap_changer_model`, allowed values `ideal`/`impedance_correction`), applied centrally to all transformers of an imported case. `ideal` (default) preserves prior Sparlectra behavior where the tap changer only changes the complex winding ratio. `impedance_correction` re-refers the transformer series impedance through the tapped winding, scaling R and X with the squared magnitude of the regulating vector, `|1 + f·e^(jφ)|²`. The option is read by both the MATPOWER importer and the native DTF importer, but the correction math lives centrally in `calcTapCorrectedRX`/`calcTapImpedanceCorrectionFactor` (`src/equicircuit.jl`) so importers stay free of duplicated tap-impedance mathematics. Configurable via YAML, the GUI-editable configuration API, and a new Web UI expert option.

# Version 0.8.9 — 2026-07-16

## New Features

* Added multi-file import for MATPOWER `.m` and DTF `.DAT` cases in the Web UI. Imported runnable cases are added to the case selector.

### Improvements

* Enabled MATPOWER DC-line handling through terminal P/Q injections by default.
* Enabled independent AC-island solving and continued per-island diagnostics by default.
* Improved Web UI file-import validation, filename handling, conflict reporting, and per-file status output.

### Bug Fixes

* Improved test-run diagnostics, bounded failure output, and excluded generated files from repository-hygiene checks.
* Made DTF `.DAT` classification content-aware so only runnable cases and valid FOR002 references appear in the corresponding selectors.
* Improved the tolerance spinner behavior in the Web UI.
* Hidden unavailable commit information instead of displaying `commit unknown`.

# Version 0.8.8 — 2026-07-14

### New Features

* Added native DTF import with support for transformer ratio conventions, Schrägregler/skew-angle controls, trailing branch records, and FOR001 contingency metadata.
* Added an experimental DTF input path for the PowerFlow API and Web UI, including optional MATPOWER export, outage selection, import diagnostics, and explicit handling of unsupported DC-line data.
* Extended MATPOWER import and export with Sparlectra metadata for bus and branch names, branch types, FOR001 contingencies, and transformer loss data.
* Added support for the existing `mpc.sparlectra.transformer_losses` extension so transformer no-load conductance can be preserved across Sparlectra–MATPOWER roundtrips.
* Added optional MATPOWER DC-line handling through fixed terminal power injections while keeping rejection of active DC-line records as the default.
* Added user-facing documentation for the DTF format, transformer conventions, Schrägregler controls, trailing records, and FOR002 validation.

### Improvements

* Improved AC-island diagnostics and failure reporting, including per-island status, reference-bus selection, solver settings, mismatch information, Q-limit state, and downloadable diagnostic artifacts.
* Improved API failure handling so non-converged runs retain their generated artifacts and remain available for ZIP download.
* Improved Linux Web UI browser startup by using Chromium app windows where available and falling back cleanly to standard desktop browser launchers.
* Added validation workflows for FOR002 cases A–E, including transformer ratios, branch identity, no-load losses, shunt sensitivity, outage handling, and MATPOWER roundtrip checks.
* Simplified the DTF/FOR002 validation examples so normal use prints a concise summary while detailed diagnostics remain available on demand.

### Bug Fixes

* Preserved transformer conductance through the existing `PowerTransformerWinding.g` → `getTrafoRXBG` / `getTrafoRXBG_pu` → `Branch.g_pu` path without creating synthetic terminal bus shunts.
* Fixed MATPOWER reimport so `mpc.sparlectra.transformer_losses` restores transformer conductance correctly and does not duplicate loss contributions.
* Fixed transformer loss reporting so total branch losses can be separated into longitudinal copper losses and voltage-dependent no-load losses.
* Fixed the Web UI case selector so internal warm-up cases are no longer shown as selectable user cases.

## Version 0.8.7 – 2026-06-23

### New Features

* Added a guarded current-injection start pre-solve for large MATPOWER workflows, including YAML/API/Web UI configuration, sidecar persistence, and `current_iteration_start.log` diagnostics.
* Added conservative MATPOWER auto-profile handling with safer default recommendations, explicit apply/skip logging, and no silent solver-start or Q-limit overrides.
* Added explicit fail-fast detection for unsupported active MATPOWER `mpc.dcline` data in API/Web UI runs.
* Added classical Q-limit enforcement modes (`classic_simultaneous`, `classic_one_at_a_time`) for diagnostic and large-case analysis.
* Added Web UI case-sidecar settings and compact operation-log retention for local support workflows.

### Improvements

* Improved power-flow diagnostics with DC-start quality metrics, compact/full mismatch summaries, clearer final-mismatch status, and Q-limit validation output with MVAr units and non-converged validity labels.
* Improved Web UI/API transparency by recording runtime case metadata, MATPOWER import decisions, Q-limit settings, output options, and partial artifact status in logs and effective configuration artifacts.
* Improved the PowerFlow Web UI layout, Advanced/Expert option grouping, dismissible validation errors, and collapsed Last Errors handling.
* Improved detailed CSV handling for large and non-converged runs, including partial exports, structured skip reasons, and artifact rediscovery.
* Restored compact default test-run output while keeping verbose diagnostics available as an explicit opt-in.

### Experimental / Developer Tooling

* Moved the large-case Q-limit comparison utility to `examples/experimental/qlimit_large_case_comparison.jl` and removed it from the stable API and normal fast test profile.

### Bugfixes

* Fixed MATPOWER import option propagation for Web UI/API runs, including auto-profile, transformer ratio, phase shift, shunt, and voltage-reference options.
* Fixed Web UI case-sidecar save/load handling with type-safe YAML values and case-local persistence.
* Fixed Q-limit diagnostics so stale/intermediate violations are not reported as final violations when the final PV/REF check is OK.
* Fixed effective configuration artifacts so they include the actual resolved runtime casefile and case name.
* Fixed CSV artifact handling for non-converged runs with available network state.

## Version 0.8.6 – 2026-06-18
### Improvements
* Added configurable detailed CSV artifact writing with `auto`, `buffered`, and `streaming` modes so large output artifacts can stream directly to disk while normal cases keep buffered writes.
* Improved direct detailed CSV export with lower-allocation cell writing and per-file timing metadata for large Web UI/API runs.
* Removed the remaining per-bus prosumer scan from direct detailed CSV bus control labels, added operation-log CSV progress events, and bounded Q-limit run-log details while preserving full Q-limit CSV artifacts.

### Bugfixes
* Fixed Web UI Q-limit checkbox handling so an unchecked form submission explicitly disables `power_flow.qlimits.enabled` in the effective configuration.

## Version 0.8.5 – 2026-06-17
### Improvements
* Cleaned up single PowerFlow Web UI/API timing output so unavailable solver time is omitted and the measured run duration is shown as wall time.

### Bugfixes
* Fixed MATPOWER `.m` import for very large cases such as `case_ACTIVSg70k.m` by replacing regex-based matrix block extraction with deterministic string scanning.

## Version 0.8.4 – 2026-06-10
### New features
* A loopback‑only PowerFlow Web UI was added with run history, artifact downloads, logs, MATPOWER case management, help, and writable runtime directories. Optional artifacts now include additional logs, configuration output, and detailed CSV exports for voltages and branch flows.
### Improvements
* Added configurable detailed CSV formatting for technical, German Excel, and US Excel notation.

### Bugfixes
* Fixed rectangular PV/slack voltage initialization so replacing a voltage magnitude preserves the existing phasor angle.

## Version 0.8.3 – 2026-05-30

### Breaking changes

* Replaced the former keyword-heavy high-level runner surface with the configuration-driven `run_sparlectra` framework entry point.
* Kept `run_acpflow` only as a thin alias for `run_sparlectra`; it now accepts the same minimal framework arguments.
* High-level import, solver, control, benchmark, and output behavior is now controlled through `SparlectraConfig` or YAML configuration.
* Framework runs now consistently return `SparlectraRunResult` for both `casefile` and `net` workflows.

### New features

* Added config-driven multi-case MATPOWER execution via `matpower_import.cases` and the dedicated `run_sparlectra_cases` helper, while keeping `run_sparlectra` as a single-case workflow.

### Improvements

* Refactored the high-level ACP/MATPOWER workflow into a clearer framework path with separated import, execution, status, and output handling.
* Refined `SparlectraRunResult` status semantics so numerical convergence, solution availability, control-loop status, limit validation, and final framework acceptance are reported separately.
* Preserved detailed rectangular diagnostics in MATPOWER runner status rows, including Q-limit active-set information, final PV voltage residuals, and wrong-branch metrics.

### Bugfixes

* Fixed file-based MATPOWER start handling so projected/imported voltage and angle starts are actually honored by the rectangular solver instead of being overwritten by an effective flat start.
* Fixed framework and MATPOWER status handling for rejected numerical solutions, including wrong-branch rejection, active-set rejection, controlled-run PF failures, disabled controls, and synthetic benchmark summaries.
* Fixed file-based MATPOWER Q-limit lock handling so `power_flow.qlimits.lock_pv_to_pq_buses` is resolved from original MATPOWER bus IDs to internal Sparlectra bus indices.
* Fixed the public `ensure_casefile` binding and updated runnable documentation snippets so fresh checkouts and package installations can load example MATPOWER cases correctly.
* Fixed the tap-control demo classic-output toggle so the example no longer fails late with an undefined `show_classic` binding.

### Related issues

* #228 Remove the old compatibility surface from the high-level runner

## Version 0.8.2 – 2026-05-29
### New Features

* Added configurable wrong-branch diagnostics for rectangular power-flow results, including voltage, angle-spread, and branch-angle plausibility checks.

### Improvements

* Hardened `matpower_import.auto_profile` into a MATPOWER pre-run that logs recommendation evidence, preserves `recommend` mode without changing the active configuration, applies safe import-convention, comparison-reference, start-mode, and Q-limit guard changes in `apply` mode, and prints final effective options without rewriting YAML files.
* Refactored the rectangular complex-state power-flow implementation into focused modules under `src/powerflow_rectangular/`, with `runpf_rectangular!` as the network-integrated entry point and `run_complex_nr_rectangular` as the standalone array-level solver.

### Bugfixes

* Fixed `run_acpflow(casefile=...)` configuration forwarding so MATPOWER/file-based rectangular solves honor configured `power_flow.wrong_branch_*` options instead of falling back to default diagnostics.
* Aligned rectangular autodamping defaults so direct solver calls and configuration-driven runs use the same `autodamp_min = 0.05` default.


### Related Issues

* #193 Harden MATPOWER auto-profile recommendation and application
* #219 Detect wrong/false low-voltage branch convergence and retry safely
* #220 Mechanically split rectangular power-flow code into focused modules
  

## Version 0.8.1 – 2026-05-26

### Breaking Changes
* Removed the old `run_net_acpflow` public wrapper. Use `run_acpflow(net=...)` for already constructed networks and `run_acpflow(casefile=..., path=...)` for file-based workflows.

### Highlights
* Added a generic outer-loop control framework above `runpf!`.
  Transformer tap/phase control now uses this framework as the first concrete controller implementation.

### Improvements
* Added structured `ControlRunResult` output and `latest_control_result(net)` for inspecting controller status, outer iterations, PF solve count, controller rows, and trace rows.
* Added machine-readable control trace rows for transformer control.
* Added `run_acpflow(net=...)` as the preferred high-level entry point for already constructed networks, and made it the single public in-memory entry path.
* Simplified `examples/tap_control_demo_grid.jl` into a lightweight demo of the generic control framework.
* Documented the `control` configuration section, including that `control.controllers` is reserved for future YAML-based controller instantiation.


### Related
#179 Introduce a generic control framework above the existing power-flow solver. 

## Version 0.8.0 – 2026-05-25

### Breaking Changes

* The public AC power-flow path now supports only the sparse rectangular Newton-Raphson solver.
  Legacy polar/classic methods, dense PF matrices, and finite-difference PF Jacobian options are no longer supported as user-facing runtime choices.
* Power-flow configuration has moved to structured YAML sections and typed configuration objects.
  The old flat keyword-style configuration path is deprecated/removed for the cleaned rectangular workflow.
* Obsolete sparse switches such as `power_flow.sparse`, `opt_sparse`, and `state_estimation.sparse` are no longer valid configuration keys.
  Sparse matrix handling is mandatory for the production PF core.

### Highlights

* Added a central configuration workflow:
  * default template: `src/configuration.yaml.example`
  * optional user override: `examples/configuration.yaml`
  * typed config objects for power flow, MATPOWER import, state estimation, diagnostics, output, benchmarking, runtime, and performance profiling
  * early validation of unknown or obsolete keys
* Simplified the production power-flow path to the sparse rectangular AC solver with sparse Y-bus assembly, sparse analytic rectangular Jacobian, and sparse linear solves.
* Improved MATPOWER import and runner workflows:
  * central YAML-driven execution
  * configurable MATPOWER import options
  * cleaner compact summaries
  * better separation of numerical convergence, Q-limit validation, and solution availability
* Added configurable performance and timing output for MATPOWER and rectangular PF runs:
  * representative wall time
  * solver time
  * result-output time
  * timing coverage
  * optional allocation information
  * rectangular workspace metadata
* Improved large-case output handling:
  * configurable result-table row limits
  * summary/compact/full result-output modes
  * safer default behavior for large MATPOWER cases
* Added rectangular workspace reuse/preallocation controls:
  * `power_flow.rectangular_workspace_reuse`
  * `power_flow.rectangular_preallocate_workspace`
  * `power_flow.rectangular_workspace_min_buses`
* Reworked example scripts under top-level `examples/`:
  * `matpower_import.jl` is now a thin YAML-driven entry script
  * `tap_control_demo_grid.jl` uses central configuration helpers
  * `export_solution.jl` writes deterministic export artifacts under `examples/_out/export_solution/<case>_<timestamp>/`
* Added a clearer test profile structure:
  * default `fast` profile for normal development
  * `extended` profile for MATPOWER/output/documentation-heavy checks
  * `all` currently aliases `extended`

### Configuration Notes

* The main configuration template is now:

  ```text
  src/configuration.yaml.example
  ```

* A local user/example override can be placed at:

  ```text
  examples/configuration.yaml
  ```

* MATPOWER benchmarking moved to the top-level benchmark section:

  ```yaml
  benchmark:
    enabled: true
  ```

  The old `matpower_import.benchmark` key is rejected with a migration message.

### Output and Examples

* `examples/export_solution.jl` now produces files instead of only printing to the console.
  Typical output files are:

  ```text
  summary.txt
  internal_solution.csv
  external_solution.csv
  comparison.csv
  *_export.m
  ```

* `examples/matpower_import.jl` uses the central configuration and writes logs under:

  ```text
  examples/_out/
  ```

### Documentation

* Added/updated documentation for:
  * central configuration
  * power-flow configuration
  * MATPOWER import configuration
  * state-estimation configuration
  * performance profiling
  * test profiles
  * examples overview

### Related

* Issue #199: Central configuration, PF solver simplification, sparse-only PF core, and test-framework cleanup
* Issue #201: YAML redesign
  
## Version 0.7.8 – 2026-05-16
### Highlights
* Improved Q-limit handling for large MATPOWER imports, especially cases with many generators that have zero or very narrow reactive-power ranges.
* Added compact console reporting for large MATPOWER example runs while keeping the full diagnostics in the logfile.

### Fixes
* Fixed rectangular power-flow status caching to use weak network keys so repeated benchmark/example solves can release imported `Net` objects after callers consume the status.
* Fixed direct `run_acpflow(net=...)` and rectangular solver defaults so the Q-limit guard remains opt-in unless a caller or config explicitly enables it.
* Fixed direct `run_tap_controllers_outer!` defaults so the Q-limit guard remains opt-in when the exported tap-control API calls `runpf!`.
* Fixed MATPOWER example config forwarding so Q-limit guard options from YAML are preserved in the effective config and passed through all `run_acpflow` paths.
* Fixed MATPOWER example console row limiting so `console_max_rows` from YAML is forwarded to Q-limit event and final active-set row caps in the `run_acpflow` paths.
* Fixed rectangular NR status reporting: numerical convergence, Q-limit active-set consistency, final convergence, comparison status, and rejection reasons are now reported separately.
* Fixed Q-limit guard behavior so strongly violating active PV buses can be locked to PQ during eligible Q-limit checks, reducing final active-PV limit violations.

### Improvements
* Reduced console noise for large Q-limit active-set cases by replacing long PV→PQ and violation tables with compact summaries.
* Improved MATPOWER import example logging: the console now shows the essential run status, while detailed auto-profile evidence, diagnostics, and solver traces remain available in the logfile. 


## Version 0.7.7 – 2026-05-13

### Highlights

* Improved MATPOWER import diagnostics for large and difficult cases:
  * clearer YAML-backed logging and terminal summaries
  * VM/VA reference checks
  * branch-shift and transformer convention scans
  * bus-shunt and fixed-reference residual diagnostics

* Added configurable MATPOWER transformer and phase-shifter import conventions:
  * `matpower_ratio`
  * `matpower_shift_sign`
  * `matpower_shift_unit`

* Improved MATPOWER voltage reference handling:
  * configurable `BUS.VM` vs. `GEN.VG` handling
  * hybrid comparison mode
  * better diagnostics for PV/REF buses without online generators
  * correct handling of buses switched from PV to PQ by Q-limits

* Improved rectangular NR robustness for large MATPOWER cases:
  * preserved slack and PV setpoints during flat start
  * added DC-angle and blended-voltage start support
  * added wrong-branch diagnostics for suspicious low-voltage solutions

* Improved PV→PQ Q-limit handling:
  * clearer per-unit and MVAr diagnostics
  * hysteresis/deadband support
  * compact post-solve active-set consistency checks
  * separate handling of PV and REF/slack Q-limit violations

### Fixes

* Fixed MATPOWER slack/reference voltage import so VM/VA values are preserved unless explicitly overridden.
* Fixed nominal-tap transformer handling so explicit `TAP = 1` branches remain transformer models.
* Fixed repository-local MATPOWER `.jl` case loading and example regression-test handling.
* Fixed several Julia 1.12 / Revise world-age issues in MATPOWER diagnostic and example paths.
* Fixed flat-start setpoint extraction for networks with isolated buses.
* Improved warning/error capture in `matpower_import.jl` with compact logfile summaries.

### Diagnostics

* Added MATPOWER auto-profile pre-run mode (`matpower_import.auto_profile = recommend|apply`) to summarize or apply robust import, flat-start, PV/REF voltage-source, and Q-limit settings while preserving explicit YAML overrides.
* Added branch-neighborhood reports for selected high-residual buses.
* Added residual-cluster diagnostics for PEGASE-style mismatch regions.
* Added negative branch impedance scans while preserving signed MATPOWER `BR_R` / `BR_X` values.
* Added MATPOWER reference-data consistency diagnostics, including documentation of the known `case300.m` fixed-reference mismatch around `BUS_I 196 / 2040`.

### Related

* Issue #186: Singular Jacobian / world-age issue in rectangular NR solver.## Version 0.7.6 – 2026-05-11
### Bugfix
* Fixed rectangular MATPOWER runs so a singular Newton Jacobian is reported as non-convergence instead of aborting the example, and tightened Julia 1.12 / Revise entry-point calls to avoid world-age binding warnings.

## Version 0.7.5 – 2026-05-10
### New Features 
* Added PV→PQ Q-limit switching start controls (`qlimit_start_iter`, `qlimit_start_mode`, `qlimit_auto_q_delta_pu`) for rectangular power-flow runs.
* Added configurable bus-shunt modeling with classic Y-bus admittance stamping and a rectangular-solver voltage-dependent injection mode for keeping shunt effects in nonlinear mismatch terms.
* Added automatic rectangular Newton damping (`autodamp`) for difficult flat-start power-flow cases; the solver backtracks the Newton step from `damp` down to `autodamp_min` and accepts the first residual-reducing trial.
* Added start projection (`start_projection`) for difficult power-flow seeds, including optional DC-angle starts and raw/DC blend scans for both the internal rectangular solver and external-solver `PFModel` starts.

### Bugfix
* Replaced the singular sparse linear-solve fallback with a rank-revealing QR path before dense SVD fallback to avoid large `pinv`/LAPACK failures in ill-conditioned rectangular Newton steps.
* Wrapped MATPOWER voltage-angle comparison differences into the minimal ±180° range before applying angle tolerances.
* Synthetic tiled-grid benchmark example now falls back from `.yaml` to `.yaml.example` and reports when built-in defaults are used.

### Improvements
* Improved MATPOWER case parsing and network construction performance for large cases by reducing parser allocations, pre-sizing network containers, and avoiding repeated bus-name lookups while adding MATPOWER branches.
  
## Version 0.7.4 – 2026-05-08
### New Features 

### Bugfix
* Fixed YAML subset parser cleanup on parse errors so temporary config files are closed before `mktemp` cleanup on Windows.

### Improvements
* Improved large-network MATPOWER and rectangular-solver performance by aggregating prosumer-derived bus types and specified power injections in linear time. 
   
## Version 0.7.3 – 2026-05-04
### New Features
* Added transformer tap control in AC power flow with an outer-loop controller workflow.
* Added `addTapController!` support for:
  * voltage control via tap ratio (`mode = :voltage`)
  * branch active-power control via phase shift (`mode = :branch_active_power`)
  * coupled ratio + phase control (`mode = :voltage_and_branch_active_power`)
  * tap-control reporting in classic and structured ACP flow outputs

### Improvements
* Updated branch-model documentation with practical controller-direction probing guidance for phase-shift control.

### Tests
* Added regression tests for transformer tap controller behavior.

## Version 0.7.2 – 2026-04-15
### Bug Fixes
* Fixed the increased loading time for large test cases

## Version 0.7.1 – 2026-04-15
### Improvements
* Refactored decision logic for `qlimit_mode = :adjust_vset`:
  * consolidated Q-limit event handling via shared active-set flow
  * extracted voltage-step handling into dedicated helper logic
* Added typed `VoltageAdjustConfig` support for prosumers and integrated it into voltage-regulation detection.
* Simplified rectangular mismatch API by removing unused derivative keyword arguments from `mismatch_rectangular(...)`.

### Solver Robustness
* Rectangular solver now handles reduced Ybus matrices (caused by internal isolated buses) by expanding them back to full network dimension for mismatch/Jacobian processing.
* For rectangular runs with active-link merges and internal isolated buses, solver now uses a rectangular FD fallback path instead of switching to `:polar_full`.

## Version 0.7.0 – 2026-04-15
### New Features
* Added support for **P(U)** and **Q(U)** controller models in power flow calculations
* Characteristic curves can be defined via discrete points, with selectable interpolation methods:
  * linear
  * piecewise constant
  * cubic
  * polynomial
* Added support for controllers on **non-PV generators**
* Added documentation and examples for the new controller features
  
## Version 0.6.4 – 2026-04-12
### New Features
* Marked legacy Jacobian solvers as deprecated:
  * `runpf_full!` / `method = :polar_full`
  * `runpf_classic!` / `method = :classic`
  Both now emit a deprecation warning when used.
* Set rectangular complex Jacobian (`method = :rectangular`) as the default for:
  * `run_acpflow`
  * `run_acpflow(; net = net, ...)`
* Updated examples and user documentation to use `:rectangular` as the recommended/default solver method.
### Bugfixes
* For backwards compatibility, if vm is set, the coresponding bus type is set to PV

## Version 0.6.3 – 2026-04-11
### New Features
* Added `:adjust_vset` controller-based Q-limit handling at PV buses (adaptive Vset steps before optional PV→PQ fallback).
* Power-flow bus typing is now derived from attached prosumers (Slack > PV > PQ); `addBus!(busType=...)` is legacy-only and no longer defines operational PF type.

## Version 0.6.2 – 2026-04-02
### New Features
* adding sign validation and optional autocorrection of Q-limits before running power flows.
* Provide an option to lock selected PV buses from being switched to PQ 
* Added pre-run PV Q-limit preview logging in MVAr for easier diagnostics before the PF iteration loop.
## Version 0.6.1 – 2026-03-24
### New Features
* Bad Data Detection (BDD) and Statistical Diagnostics for State Estimation (SE)
  
### Bug Fixes
* Fix Issue 139

### Notes

* Bad-data detection and statistical diagnostics are available via
  `validate_measurements`, `runse_diagnostics`,
  `summarize_se_diagnostics`, and `print_se_diagnostics`.

## Version 0.6.0 – 2026-03-17
### New Features
* Added initial State Estimation (SE, WLS), documentation and examples.
* Added zero-injection-bus (ZIB) support in SE, documentation, and examples.
* Improved results reporting and logging.
* Refactor testsuite for clarity.

### Notes
* State Estimation is currently marked as **experimental**.
* Passive buses / ZIB are currently represented through zero-injection pseudo measurements in the WLS workflow.


## Version 0.5.0 – 2026-03-11
### Highlights
* Create and use machine-readable `ACPFlowReport` output
* Introduced bus-links +  documentation 
* Consolidated recent documentation

## Version 0.4.35 – 2026-02-02
### Fixes
* Stabilized MATPOWER case import and bus indexing in the presence of isolated buses.
* Fixed bus classification and ordering to ensure consistent solver input for flatstart and non-flatstart runs.
* Improved internal consistency checks for network topology before solver execution.

### Improvements
* Refined MATPOWER comparison tooling (Vm/Va diff diagnostics and angle alignment).
* Minor internal cleanups in MatpowerIO and network initialization to improve robustness.

## Version 0.4.34 – 2026-06-02
### Fixes
* Closed Issue 110:
* Corrected per‑unit conversion issues discovered during the refactoring of the MATPOWER import interface, including proper handling of line, transformer, and shunt parameters.
* Fixed PU conversion inconsistencies for shunts and aligned all shunt handling with the unified Y‑model.
* Harmonized Q‑limit checking across all solvers and ensured consistent PV/PQ switching behavior.

### New Features
* Updated shunt modeling: addShunts now exclusively supports the Y‑model; voltage‑dependent shunt power must be represented via loads or generators.
* Unified solver logic for reactive power limits and PV/PQ switching
* added new testcase matpower vs manual network

## Version 0.4.33 – 2026-02-02

### New Features
* Added on-demand MATPOWER case handling via `FetchMatpowerCase.ensure_casefile`
  (automatic download of `.m` cases and optional generation of `.jl` cases)
* Added option to control flat start behaviour in AC power flow (`flatstart`)

### Internal Changes
* Refactored MATPOWER case import logic and removed redundant example-based loaders

## Version 0.4.32 – 2026-30-01
## New Features
* Added function to export calculated network
* 
## Version 0.4.31 – 2026-28-01
### Title
* Change license from BSD-3-Clause to Apache License 2.0

## Version 0.4.30 – 2026-28-01
### Bug Fixes
* Importing Matpower files, wrong conversion to per unit system 
### New Features
* Option to choose flatstart 

## Version 0.4.29 – 2025-12-13
### New Features
* Support native 3-winding transformers in Sparlectra without AUX bus generation and allow non-PU parameterization #90
* Support 2-winding transformers as well in non-PU parameterization #90
  
## Version 0.4.28 – 2025-12-12
### Bug Fixes
* Close Issue #85: test totalBusPower vs. TotalLosses failed for solver != rectangular

## Version 0.4.27 – 2025-12-11
### New Features
* Added 3Bus Testcase with PV Generator
* Added createTest5BusNet with multiple generators at one bus and zero injection generator
* Issue #81 Fix Calculatation of Generators with Zero Injektion
* Mark PV->PQ switching in printACPFlowResults Log 
  
### Bug Fixes
* fixes per unit calculation for line shunt admittance
* fixes calculation of losses


## Version 0.4.26 – 2025-12-04
* Issue #74: Multiple Generators at One Bus Not Handled Correctly
* Documentation Updates 

## Version 0.4.25 – 2025-11-29

### Added
* Rectangular (Complex-State) Newton–Raphson Solver

## Version 0.4.24 – 2025-11-20

### Added
- Q-limit handling with automatic PV→PQ switching when a generator hits its reactive power limits.

## Version 0.4.23 (2025-11-11)

###  New Features
* Full-state Newton–Raphson solver

## Version 0.4.22 (2025-08-27)
### Bug Fixes
 - small fixes

## Version 0.4.21 (2025-03-14)
### New Features
 - adding functions to remove elements from a net

## Version 0.4.20 (2025-03-11)
 -  internal reorganization and small bugfixes 
 
## Version 0.4.19 (2024-10-14)
### Bug Fixes
 - pu calculation for transformer impedance

## Version 0.4.18 (2024-04-14)

### Bug Fixes
 - closes issue#48 "printout jacobian runs into error"

## Version 0.4.17 (2024-04-14)
### New Features
 - added testcase for importing Matpower files

### Bug Fixes
 - bugfix wrong function call in `createnet_powermat`

## Version 0.4.16 (2024-04-13)
### Bug Fixes
- bugfix shunt index for isolated buses, closes issue #38

## Version 0.4.15 (2024-04-12)
### New Features
- Implemented a function to detect isolated buses and incorporate them into the network calculation (Issue #38)

## Version 0.4.14 (2024-04-12)
### Bug Fixes
- bugfix addShuntPower, closes issue #36

## Version 0.4.13 (2024-04-12)
### New Features
- added attribute for Lineparameters based on length
- added update parameter function for network
- added workshop documentation

### Bug Fixes
- taking line length not (always) into account for line parameters
- parsing emtpy lines of Matpowerfiles
- documentation rendering

## Version 0.4.12 (2024-04-08)
### New Features
- added functions to facilitate the modification of networks.
- documentation available at https://welthulk.github.io/Sparlectra.jl/.

### Bug Fixes
- print prosumers

## Version 0.4.11 (2024-04-05)
### New Features
- make changes to imported Matpower networks after import.
- added functions to facilitate the creation of networks.

### Enhancements
- added documentation make file

### Bug Fixes
- import and parser for Matpower .m files

## Version 0.4.10 (2024-03-30)

### New Features
 - removed numerous redundant functions, partially restructured classes
 - removed support for CGMES due to the availability of numerous alternative libraries
 - removed support for the legacy custom JSON data format (potentially subject to reintroduction at a later stage)
 - added functions to facilitate the creation of networks
 - better performance
 
### Bug Fixes
- calculation of branch flow and losses
- branches out of service

## Version 0.4.8 (2024-03-26)
- first package release registered in the Julia registry

## Version 0.4.1 (2023-12-19)
- Initial release of Sparlectra

## Version 0.4.0 (2023-11-30)
- Initial public commit of Sparlectra

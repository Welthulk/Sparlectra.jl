# Feature Matrix (Quick Overview)

Sparlectra covers four analysis domains: power flow, balanced short
circuit, N-1 contingency analysis, and state estimation. Each domain has
its own table; the shared network model, the exchange formats, and the
workflow tooling follow below.

Legend:

* ✅ available
* ⚠️ available with limitations / specific workflow
* ❌ not available as dedicated feature

## Power flow

| Feature | Status | Notes |
|---|:---:|---|
| Rectangular NR solver (`runpf!`) | ✅ | Main PF entry point with the sparse rectangular complex Jacobian; polar/classic PF methods are not supported. |
| DC power flow (`rundcpf!`, `power_flow.solver = dc`) | ✅ | Standalone linear screening model (MATPOWER `rundcpf`/`makeBdc` equivalent): series-reactance-only `B'`, phase-shift injection vector, `Vm` implicitly 1.0 pu, lossless. Usable standalone or as the framework solver (including per-island handling). `rundcpf!(net; seed_ac_start=true)` optionally chains an AC Newton-Raphson solve seeded from the DC angles. No outer-loop controller support (rejected up front). |
| APSLF solver (`power_flow.solver = apslf`) | ⚠️ | Analytic power-series solver; usable standalone (`runpf_external!`), as the framework solver, or as a guarded start-value generator ahead of the rectangular NR solve (`power_flow.apslf_start`). The warning marker: it lives in the optional weak dependency AnalyticLoadFlow.jl (must be installed separately), has no selectable start voltage (always the canonical analytic germ), no OLTC/PST/Q(U)/P(U) controller support, and only simple internal PV→PQ Q-limit switching instead of the active-set/outer-loop modes. See [External Solver Interface](external_solvers.md#apslf-analyticloadflowjl). |
| PV/PQ reactive limit handling | ✅ | Default active-set Q-limit logic plus classical simultaneous and one-at-a-time outer-loop modes. |
| Automatic Newton damping (`autodamp`) | ✅ | The rectangular solver can backtrack the Newton step from `damp` down to `autodamp_min` for difficult flat starts. |
| Merit-function Armijo line search (`power_flow.merit`) | ✅ | Optional alternative step-acceptance criterion inside the autodamp backtracking loop (`f(x) = 1/2‖WF(x)‖²`, Armijo sufficient decrease); disabled by default and requires `autodamp = true`. |
| Trust-region step control (`power_flow.trust_region`) | ✅ | Optional alternative to `autodamp`: caps the Newton step norm at an adaptive radius, accepts/rejects by merit decrease, adapts the radius from the actual/predicted reduction ratio; `step_mode = :dogleg` blends toward a Cauchy step when the Newton direction degrades. Disabled by default, mutually exclusive with `autodamp = true`. |
| AC rescue ladder and DC fallback (`power_flow.rescue`, `power_flow.dc.fallback`) | ✅ | A non-converged AC solve is retried from the original start state with a fixed strategy ladder (alternate start, autodamp, DC-seeded projection); if nothing converges, the standalone DC power flow can leave usable angles and branch P flows. The AC status stays non-converged; strategy and fallback are logged. |
| Start projection and start modes | ✅ | DC-angle and blend-scan start projection with candidate measurement, imported-profile starts, and the optional guarded current-iteration pre-solve (`power_flow.start_current_iteration`). |
| Distributed active-power slack (`power_flow.distributed_slack`) | ✅ | The island's P imbalance is shared over participating generators by normalized participation factors instead of loading it entirely onto the reference bus. Weight modes `pg_weighted` (default), `pmax_weighted`, `headroom_weighted`, `imported` (MATPOWER `APF` / CGMES `normalPF`), `explicit`. Disabled runs are bit-identical to the classical solver. Theory in [Solver Guide](solver.md). |
| Automatic slack selection (`power_flow.auto_slack`) | ✅ | When a case registers no voltage reference, the solver promotes the strongest injection to slack instead of aborting; off by default so data errors stay visible. API: `ensureSlack!`. |
| External grid element (`addExternalGrid!`) | ✅ | Native IEC 60909-0 network feeder (issue #299): ideal slack by default, optional non-ideal load flow (`internal_impedance = true`) with the reference behind `z = Un²/Sk''`; the declared short-circuit data feeds `runShortCircuit!` on hand-built and MATPOWER nets. Theory on [Slack Bus and External Grid Sources](slack_vs_source.md). |
| Parallel island solving (`power_flow.islands.mode: solve_parallel`) | ✅ | Detected AC islands solve concurrently on Julia threads, gated by `runtime.parallel.*`; results are bitwise identical to the serial mode. Start Julia with `--threads=auto`. |
| Voltage-dependent prosumer control (`Q(U)`, `P(U)`) | ✅ | Controller-aware mismatch/Jacobian terms in the rectangular formulation. |
| Wrong-branch detection (`wrong_branch_detection`) | ✅ | Post-convergence plausibility guard for suspicious low-voltage/non-finite solutions (`off|warn|fail`); the result is surfaced in `ACPFlowReport.metadata`, the island diagnostics CSV, the console summary, and the Web UI run page. |
| Jacobian condition diagnostics (`condestJacobian`, `reportCondition`) | ✅ | Hager/Higham 1-norm condition estimate on the LU factorization with a digits-lost verdict, to tell an ill-conditioned Jacobian from a start/step-control problem. See [Solver](solver.md). |
| Narrative diagnostics and self-check (`diagnose.log`, `run_fixed_reference_self_check`) | ✅ | The diagnostics report names the worst-mismatch bus, the mismatch-history trend, autodamp health, and branch anomalies around the worst bus; the self-check evaluates a case's own stored operating point with every start-value machine forced off, separating model issues from solver issues. |
| External solver interface (`PFModel` / `PFSolution`) | ✅ | Integration point for solver backends outside the built-in rectangular NR. |

## Short circuit (IEC 60909-0)

Balanced short-circuit analysis stands on its own entry point
(`runShortCircuit!`) and its own result type; it needs no converged power
flow.

| Feature | Status | Notes |
|---|:---:|---|
| Initial symmetrical current per fault bus | ✅ | `Ik''` (max/min case), `Sk''`, and peak current `i_p`, positive sequence, series impedances only, IEC Table-1 voltage factors with `short_circuit.c_factor` as scalar expert override. See [Short-Circuit Analysis](short_circuit.md). |
| Data sources | ✅ | The always-on CGMES short-circuit harvest (machines, feeders, motors) and the native feeder records written by `addExternalGrid!`, field-identical contracts. |
| Safety-flag contract | ✅ | Substituted defaults and skipped contributions are flagged on the affected result rows, not only logged; a flagged maximum is a lower bound. |
| All-bus sweeps on Julia threads | ✅ | The fault-bus list fans out over task chunks (one factorization copy and reusable buffers per chunk, gated by `runtime.parallel.*`), row-identical to the serial sweep. |
| Takahashi sparse inverse (`sweep_method = :takahashi`) | ✅ | The whole Thevenin diagonal of an island from one selected-inverse pass over the LU factors instead of one triangular solve per bus; measured 34x to 264x over the serial sweep between 2000 and 16000 buses. Agrees with the default `:solves` to machine precision; inapplicable islands fall back automatically. Theory in the [Short-Circuit Compendium](short_circuit.md). |
| Web UI integration | ✅ | The "Short circuit" button runs both cases without a power-flow solve and writes `short_circuit_max.csv`/`short_circuit_min.csv`. |
| Unbalanced faults, `K_T`/`K_G` corrections | ❌ | Documented limitations; positive-sequence balanced faults only. |

## N-1 contingency analysis

| Feature | Status | Notes |
|---|:---:|---|
| Branch-outage batches (`runContingencies!`) | ✅ | Cases run warm-started on template copies of the solved base case (the base net is never mutated) and are checked against the voltage band and `sn_MVA` loadings. See [N-1 Contingency Analysis](contingency.md). |
| Case sources | ✅ | `generateN1Branches` (all in-service branches, transformer filter, parallel circuits disambiguated as `name#branchIdx`) and imported MATPOWER FOR001 lists. Screening filters `min_vn_kV` / `min_sn_MVA` / `name_pattern` keep the list to the outages worth simulating. |
| Generator outages | ✅ | `generateN1Generators` builds `kind = :gen` cases (one unit each, `min_pg_MW` / `name_pattern` filters). The slack absorbs the loss, or `distributed_slack_enabled = true` shares it; `auto_slack = true` promotes a survivor when the slack unit itself is the outage, and stranding an island with injection but no reference is reported as generation stranded. |
| Case weights | ✅ | Each case carries a `weight` (default 1.0) for a probability/severity-weighted ranking, carried into the result and CSV; `readContingencyWeightsCSV` + `applyContingencyWeights` attach per-branch outage rates. |
| Failure semantics | ✅ | Islanding without a promotable reference, non-convergence, and unresolvable elements are reported per case instead of thrown; a per-case start-value ladder (`contingency.rescue_ladder`) tries several starts, and a load-only island is reported as a quantified load-shed result. |
| Parallel execution | ✅ | The case list fans out over Julia threads (`runtime.parallel.*`), results identical to the serial run; full N-1 on case1354pegase measured 71.7 s serial vs 17.6 s on 16 threads. |
| Output | ✅ | `printContingencyResults` table (severity-ranked, failures first) and `writeContingencyResultsCSV`, both carrying the overload loadings with base/delta (`OverloadRecord`), shed load, weight, and severity per case. |
| Aggregate report | ✅ | `buildContingencyReport` folds a batch into a `ContingencyReport` (outcome counts, total/worst load shed, worst branch loading, worst weighted severity, most-overloaded branches); `printContingencyReport` prints it. |
| Web UI | ✅ | A "Contingency (N-1)" button with a branch/generator selector runs the batch through the shared service/run/case-cache path (a mode flag on `POST /powerflow/run`, no separate workflow), writes `contingency_n1.csv` + a `run.log` report, and shows an outcome summary; a slack-unit outage is named, not shown as a failure. |

## State estimation (WLS)

State estimation is its own workflow around `runse!`: it reconstructs the
network state from redundant, noisy measurements and brings its own
measurement model, observability analysis, and bad-data diagnostics. It
shares the network model and the importers with the power flow.

| Feature | Status | Notes |
|---|:---:|---|
| Nonlinear WLS estimator (`runse!`) | ✅ | Iterative WLS on the shared `Net` model; `updateNet = true` writes the estimated state back. |
| SCADA-style measurements | ✅ | `Vm`, `Pinj`, `Qinj`, `Pflow`, `Qflow` types with helper builders. |
| PMU voltage phasors | ✅ | Magnitude as tightly weighted `VmMeas`, angle as `VaMeas` (degrees) with an estimated common reference-angle offset between the PMU time base and the slack reference (`state_estimation.pmu_ref_offset`). Combined helper: `addPmuPhasorMeasurement!`. |
| Zero-injection buses | ⚠️ | Modeled as tightly weighted zero-injection pseudo-measurements. The warning marker: a pseudo-measurement enforces the exact P=Q=0 knowledge only up to its weight, and extreme weights degrade the gain-matrix conditioning; the textbook alternative (hard equality constraints, Hachtel/Lagrange formulation) is not implemented. |
| Observability analysis | ✅ | Structural (sparsity/matching), numerical (rank/SVD), global and local checks. |
| Bad-data diagnostics | ✅ | Global consistency, residual ranking, `runse_diagnostics`, `summarize_se_diagnostics`, `print_se_diagnostics` (`format = :markdown` or `:plain`), and a deactivate-and-rerun helper for the top suspicious measurement. |
| Synthetic measurements from a PF result | ✅ | `setMeasurementsFromPF!` builds SE test sets from a solved power flow, with or without noise. |
| Flat start control | ✅ | Same start discipline as the power flow. |
| PMU current phasors | ❌ | Not to be confused with the fully supported PMU VOLTAGE phasors above: branch current phasors (magnitude and angle of `I_ij`) have no dedicated measurement types yet; the `MeasurementType` enum covers `Vm`, `Va`, `Pinj`, `Qinj`, `Pflow`, `Qflow`. |
| Estimation-side controllers and Q-limit switching | ❌ | The measurement model estimates the state it is given; tap/voltage controllers and PV/PQ switching are power-flow concepts and stay there. |

## Controllers and FACTS (power-flow outer loop)

All controllers run in the generic outer loop above `runpf!`; results are
available through `ControlRunResult` / `latest_control_result(net)`, the
machine-readable trace rows on `ControlRunResult.trace`, and the uniform
`controllableElements` view (element, device, actuator with range, target,
live status). The device taxonomy and the limit-characteristic comparison
(constant-Q vs `V·S_max` vs `V²·B`) live on [FACTS Devices](facts.md).

| Feature | Status | Notes |
|---|:---:|---|
| Transformer regulation (OLTC / PST / combined) | ✅ | Voltage, branch-active-power, and combined modes, discrete steps with limits, split regulation with two independent controllers on one unit; the full per-capability breakdown is the [transformer support table](@ref transformer-support) below. |
| Remote voltage control by machines (`MachineVoltageControl`) | ✅ | A PQ machine regulates the voltage at a different bus via its reactive output (secant outer loop, `at_limit` on the reactive bounds); wired from CGMES remote `RegulatingControl`s with `cgmes_import.machine_control`. One controller per target bus. Theory in [Remote Voltage Control](remote_voltage_control.md). |
| STATCOM current-based limit mode (`s_max_mva`) | ✅ | The machine controller as a VSC shunt compensator: the reactive bound is the converter current, `Q_lim = V * S_max`, re-evaluated from the solved terminal voltage every outer iteration; at the limit the delivered Q tracks the voltage linearly (versus the SVC's quadratic `V² B` collapse). Device taxonomy and comparison on [FACTS Devices](facts.md). |
| SVC variable-shunt voltage control (`ShuntVoltageControl`) | ✅ | A continuous shunt-susceptance actuator holds the local bus voltage via secant iteration; at a limit the susceptance clamps and the reactive output follows V² through the Y-bus, the constant-B region of a real SVC. Imported CGMES SVCs still map as static PV injections. |
| MSC/MSR switched shunt banks (`step_mvar`) | ✅ | The shunt controller as a mechanically switched bank: the susceptance moves in whole blocks, truncated toward the target so the bank never overshoots (anti-hunting by construction), parking on the last step before crossing (`status = :parked`); the outermost block keeps the constant-B limit region. See [FACTS Devices](facts.md). |
| TCSC series-reactance flow control (`SeriesReactanceControl`) | ✅ | The series reactance of a line branch is the actuator, the branch active power the target; secant iteration with a bounded bootstrap probe, `at_limit` at a range end. Transformer branches are rejected (taps own transformer reactance). Theory in [Series Compensation (TCSC)](series_compensation.md). |
| SSSC injected-voltage limit mode (`v_inj_max_pu`) | ✅ | The series controller as a VSC series compensator: the admissible reactance deviation is bounded by the injectable series voltage, `\|x - x_base\| <= V_inj,max / \|I\|`, so the window shrinks with loading; at the limit the injected voltage sits at `V_inj,max` and tracks the branch current. |
| UPFC (combined shunt + series converter) | ⚠️ | Two models via `addUpfcControl!` (YAML type `upfc`). `model = :quadrature` (default): SSSC + STATCOM composite, one line quantity, no series active power. `model = :full` (#326): arbitrary-phase series injection steering line P AND Q independently, with the DC-link active balance (`P_se + P_sh = 0`) on the shunt. The warning marker on the full model: the shunt runs on a reactive setpoint, not closed-loop voltage (deferred to the sensitivity framework #217), no explicit series current limit, and convergence is reliable for feasible moderate targets but not for aggressive targets near the injection limit; stationary model, IPFC out of scope. Details on [FACTS Devices](facts.md). |
| Equipment impedance vs FACTS operating point (`r_base_pu`/`x_base_pu`) | ✅ | A series-FACTS control run (TCSC/SSSC/full UPFC) stamps its compensated impedance onto the live branch fields; the power flow uses that, but `runShortCircuit!` and the CGMES/MATPOWER exports read the physical base, so a fault study or export on the compensated net matches the equipment network with no manual reset (issue #329). `restoreBaseImpedances!` and `clearUpfcFullControllers!`/`clearSeriesReactanceControllers!` reset the live branch to the base. See [FACTS Devices](facts.md). |
| HVDC back-to-back pairing (`HvdcPairControl`) | ✅ | The two converter injections of one link stay coupled by the pairing invariant `P_to = P_transfer - loss`, per-terminal fixed Q or voltage-target secant, optional transfer rating with `at_limit`; no angle coupling, HVDC-joined areas stay separate islands. Opt-in from both importers (`matpower_dcline_mode` / `hvdc_mode` `= paired_control`); a grid-forming variant (`mode = :island_feed`) feeds an island whose only source is the receiving converter. Theory in [HVDC Back-to-Back](hvdc_back_to_back.md). |
| YAML controller instantiation | ✅ | Named controller definitions under `control.controllers`. See [Control Framework](control_framework.md). |

### [Transformer support](@id transformer-support)

Transformer types and regulation features by winding configuration. A 3-winding
transformer is modeled as a star equivalent with an internal AUX bus
(`add3WTPiModelTrafo!` / `create3WTWindings!`); tap and phase changers sit on a
chosen winding (`tap_side`, `phase_tap_side`) and are controlled per
star-equivalent leg branch. All tap/voltage control features act in the PF
outer loop only; SE uses the same transformer network model but has no
controller support.

| Type | 2-winding | 3-winding | Remarks |
|---|:---:|:---:|---|
| Fixed-tap transformer | ✅ | ✅ | Base network model, usable in PF and SE. |
| OLTC (voltage control on ratio tap) | ✅ | ✅ | `addTapController!` with `mode = :voltage`; discrete step operation with tap limits. Remote target-bus control via `target_bus` (⚠️ one controller channel with remote measurement). |
| PST (active-power control on phase tap) | ✅ | ✅ | `mode = :branch_active_power`; discrete step operation with phase limits. |
| Combined regulation, single controller | ✅ | ✅ | `mode = :voltage_and_branch_active_power`: one controller drives ratio and phase taps together. |
| Split combined regulation (Schrägregler) | ✅ | ✅ | Two independent controllers on one unit (voltage on the ratio tap plus active power on the phase tap), each with its own target, deadband, and convergence status; per-actuator exclusivity is enforced. Demo: `tap_control_schraeg_two_controllers.jl`; theory in `control_framework.md`. |
| Symmetrical PST model (CGMES) | ⚠️ | ⚠️ | Typed `PhaseTapChangerModel`. The warning marker: the typed model computes ratio and shift (the DTF importer uses it), but a persisted `phase_taps` model does not yet drive the solved branch on its own; details in the note below. |
| Asymmetrical PST model (CGMES) | ⚠️ | ⚠️ | Includes quadrature booster as ψ = 90°. Same restriction as the symmetrical model: modeling and import layer, no autonomous effect on the solved branch yet; details in the note below. |
| Tabular PST model (`:tabular`) | ✅ | ⚠️ | `TapTablePoint` lookup; tabular data overrides formula reconstruction. 2-winding: CGMES `PhaseTapChangerTabular` is wired into the imported branch (ratio and angle at the tap position; per-step r/x/g/b corrections reported, not applied) and validated against RealGrid's SV state. The 3-winding warning marker: the model can be defined on a star-equivalent winding (`create3WTWindings!` with `phase_taps`), but it shares the persisted-model restriction from the note below and has no comparable real-case validation yet. |
| Coordinated master/slave voltage control | ✅ | ✅ | Parallel transformers regulate as a GROUP (`followers` on `addPowerTransformerControl!`): the master runs the discrete voltage loop, followers mirror step-synchronously, which keeps the loop free of circulating reactive power. CGMES deliveries with several tap changers on one shared `TapChangerControl` import as one group instead of fighting controllers. No participation-factor allocation yet (pure position mirroring). See [Control Framework](control_framework.md). |
| Tap-changer impedance model (`tap_changer_model`) | ✅ | ✅ | `ideal` (default) keeps the tap changer free of series-impedance feedback; `impedance_correction` re-refers transformer R/X through the tapped winding (`|1 + f·e^(jφ)|²`). Applies to all transformers of an imported case, read by both the MATPOWER and the native DTF importer; implemented centrally in `calcTapCorrectedRX`/`calcTapImpedanceCorrectionFactor` (`src/equicircuit.jl`). SE reads the same imported `Net`, so the model choice is inherited but not independently configurable per SE run. |

!!! note "Status of the typed CGMES phase-tap-changer models"
    Phase-tap-changer models (`:symmetrical`, `:asymmetrical` (quadrature booster as ψ=90°), and `:tabular` with `TapTablePoint`) can be defined on transformer windings, directly (2WT) or via `create3WTWindings!`'s `phase_tap_side`/`phase_taps` keywords (3WT, see `examples/others/exp_3wt_phase_taps.jl`). The DTF importer uses these models to derive the branch ratio and shift. Limitations: a persisted `phase_taps` model does not yet change the solved branch on its own, and there is no per-transformer configuration surface yet.

## Network model and data exchange

| Feature | Status | Notes |
|---|:---:|---|
| Shared network model | ✅ | Buses, π-equivalent branches, transformers, generators, loads, shunts, and links in one `Net` used by every analysis domain. |
| Per-terminal branch status | ✅ | A branch open at exactly one terminal stays in the model as its exact pi reduction (full charging draw, open-end voltage as a branch result); CGMES `Terminal.connected` maps to the flags in both directions, MATPOWER export substitutes the exact `Y_in` bus shunt. |
| Topological bus links (`addLink!`) | ✅ | Impedance-less busbar couplers/sectionalizers: closed-link clusters contract onto one bus before the Y-bus is built, link flows are reconstructed after the solve (`calcLinkFlowsKCL!`, minimum-norm in zero-impedance loops). See [Links](links.md). |
| Configurable bus-shunt modeling | ⚠️ | `matpower_import.bus_shunt_model = "admittance"` (default) stamps bus shunts into the Y-bus; `"voltage_dependent_injection"` keeps them in the nonlinear mismatch terms instead. The warning marker: the injection variant is a MATPOWER-import option, and it refuses to run when active link merges are in play (both the rectangular and the APSLF path reject that combination with an error). |
| CGMES import (`importCGMES` / `createNetFromCGMES`) | ✅ | ENTSO-E CGMES 2.4.15 bus-branch import (EQ+SSH+TP+SV, boundary sets, folders/ZIP/ZIP-in-ZIP): lines incl. boundary lines across nominal-voltage steps, 2W/3W transformers with fixed tap positions or CGMES-defined outer-loop tap controllers, machines with remote voltage controllers, retained switches as bus links, always-on short-circuit data harvest, `summarizeCGMES` diagnostics and `compareWithSV` validation against the shipped SV profile. CGMES **3.0** deliveries are read as well (`dcat:Dataset` headers, per-border boundary files, SSH `Equipment.inService`); multi-area assemblies discard cancelling boundary equivalent pairs, and `VsConverter`/`CsConverter` map as Stage-0 fixed PCC injections. Placeholder guards keep completeness-set filler values out of the solve. Validated on MicroGrid, SmallGrid, FullGrid, and the 6209-bus RealGrid. |
| Node-breaker topology processor | ✅ | Imports node-breaker deliveries WITHOUT a TP profile: connectivity nodes aggregate across closed non-retained switches (SSH state overriding `normalOpen`, out-of-service counts as open), retained switches stay bus couplers, boundary nodes adopt the TP_BD set. Runs only when no non-boundary `TopologicalNode` exists; verified against the shipped TP on the MiniGrid/SmallGrid/FullGrid node-breaker conformity sets. See [CGMES Import](cgmes_import.md). |
| CGMES export (`writeCGMESFiles`) | ✅ | Complete CGMES 2.4.15 delivery, EQ + TP + SSH + SV (optionally one re-importable ZIP) with roundtrip-stable object identity: mRIDs recorded on import are reused, everything else is minted deterministically. An exported and re-imported net solves to the same power flow and reproduces the original short-circuit evaluation. See [CGMES Export](cgmes_export.md). |
| CGMES import analysis (`analyzeCGMES`) | ✅ | Explains a non-importable delivery: supplied models, declared prerequisites matched against the input, unresolved-reference histogram, plain-language verdict. |
| Base-voltage inference (`cgmes_import.infer_base_voltages`) | ✅ | Reconstructs missing nominal voltages from the SV state and transformer rated voltages when a delivery ships without its `BaseVoltage` catalog. |
| MATPOWER import / export | ✅ | Configurable SHIFT unit/sign and TAP ratio conventions, transformer-loss metadata round trips, auto-profile recommendations for robust large-case settings, `writeMatpowerCasefile` with optional solved-state columns. |
| Native DTF import | ✅ | Native `.dat` network cases incl. FOR001/FOR002 validation workflows. See [DTF Format](dtf_format.md). |
| Synthetic tiled-grid generator | ✅ | `build_synthetic_tiled_grid_net` creates artificial one-voltage-level benchmark networks. |

## Workflow, reporting, and tooling

| Feature | Status | Notes |
|---|:---:|---|
| Framework workflow (`run_sparlectra`) | ✅ | Configuration-driven import/control/solve/output orchestration, one `SparlectraRunResult` per run; `run_sparlectra_cases` executes configured MATPOWER batches in order. |
| Central typed configuration | ✅ | `SparlectraConfig` with cached YAML loading, typed validation, override precedence, and effective-configuration printing. |
| Parallel runtime (`runtime.parallel.*`) | ✅ | One switch set gates every threaded surface (islands, short-circuit sweeps, contingency batches): `enabled`, `max_tasks`, `min_work_items`. Serial fallbacks are the same functions. |
| Machine-readable report (`ACPFlowReport`) | ✅ | DataFrame-friendly rows for buses, branches, links, transformer controls, Q-limit events, and HVDC links. |
| GUI-ready programmatic run API | ✅ | `run_sparlectra_api` with stable run IDs, schema-versioned status, controlled configuration overrides, and explicit artifact discovery. |
| Local PowerFlow service boundary | ✅ | `start_powerflow_run`, persistent run indexing, restart recovery, result lookup, and safe artifact resolution, without HTTP dependencies. |
| Local browser Web UI | ⚠️ | PowerFlow forms with case management, contextual help, run history, and artifact viewing. The warning marker names the deliberate scope, not a defect: the server binds to loopback only, has no authentication and therefore no public deployment mode, and state estimation has no Web UI page yet. |

## Useful links

* [State Estimation](state_estimation.md)
* [Short-Circuit Compendium](short_circuit.md)
* [N-1 Contingency Analysis](contingency.md)
* [FACTS Devices](facts.md)
* [Links (bus couplers)](links.md)
* [External Solvers](external_solvers.md)
* [Network Reports](netreports.md)
* [Branch Model](branchmodel.md)
* [Workshop tour](generated/workshop_tour.md)
* [Changelog](changelog.md)

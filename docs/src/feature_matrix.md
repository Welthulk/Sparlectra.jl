# Feature Matrix (Quick Overview)

Sparlectra covers four analysis domains: power flow, balanced short
circuit, N-1 contingency analysis and state estimation. Each has its own
table; the shared network model, the exchange formats and the workflow
tooling follow.

Legend:

* ✅ available
* ◐ partial: not yet end to end, the note says what is missing
* ❌ not available: the note says what is offered instead

A row whose note starts with `Scope:` is available with a deliberate boundary.
Each table lists only concepts of its own domain; a missing row is not a gap.

## Power flow

| Feature | Status | Notes |
|---|:---:|---|
| Rectangular NR solver (`runpf!`) | ✅ | Main entry point, sparse rectangular complex Jacobian; no polar/classic PF method. [Solver](solver.md). |
| DC power flow (`rundcpf!`, `power_flow.solver = dc`) | ✅ | Linear screening model (MATPOWER `rundcpf`/`makeBdc` equivalent), standalone or as the framework solver per island; `seed_ac_start=true` chains an AC solve from the DC angles. No outer-loop controllers. |
| APSLF solver (`power_flow.solver = apslf`) | ✅ | Analytic power-series solver (AnalyticLoadFlow.jl): standalone (`runpf_external!`), framework solver, or guarded start-value generator (`power_flow.apslf_start`). Scope: plain power flow only; a network with controllers solves without them (one warning names them), simple PV to PQ switching. [External Solver Interface](external_solvers.md). |
| PV/PQ reactive limit handling | ✅ | Active-set logic by default, classical simultaneous and one-at-a-time outer-loop modes. [Power limits](powerlimits.md). |
| Automatic Newton damping (`autodamp`) | ✅ | Backtracks the Newton step from `damp` down to `autodamp_min`. |
| Merit-function Armijo line search (`power_flow.merit`) | ✅ | Alternative step acceptance inside the autodamp loop; off by default, requires `autodamp = true`. |
| Trust-region step control (`power_flow.trust_region`) | ✅ | Alternative to `autodamp` with an adaptive step-norm cap and `step_mode = :dogleg`; off by default, mutually exclusive with `autodamp = true`. |
| AC rescue ladder and DC fallback (`power_flow.rescue`, `power_flow.dc.fallback`) | ✅ | A non-converged AC solve is retried through a fixed strategy ladder; if nothing converges, the standalone DC power flow can leave usable angles and branch P flows, the AC status staying non-converged. |
| Start projection and start modes | ✅ | DC-angle, blend-scan and (on a flat start) transformer-ratio candidates, imported-profile starts, optional guarded current-iteration pre-solve (`power_flow.start_current_iteration`). [Start strategies](start_strategies.md). |
| Distributed active-power slack (`power_flow.distributed_slack`) | ✅ | The island's P imbalance is shared over participating generators; weight modes `pg_weighted` (default), `pmax_weighted`, `headroom_weighted`, `imported`, `explicit`. [Solver Guide](solver.md). |
| Automatic slack selection (`power_flow.auto_slack`) | ✅ | Promotes the best injection to slack when a case registers no reference; off by default. API: `ensureSlack!`. |
| Reference priority (`referencePriority`) | ✅ | Integer per unit (1 strongest, 0 none) ordering the units that may take over a reference; without priorities the strongest unit wins (`reference_candidate_rank`). Imported from MATPOWER type 3 and CGMES, carried by SCF and the CGMES export. [Reference priority](slack_vs_source.md#Reference-priority). |
| External grid element (`addExternalGrid!`) | ✅ | IEC 60909-0 network feeder: ideal slack by default, optional `internal_impedance = true`; its short-circuit data feeds `runShortCircuit!`. [Slack Bus and External Grid Sources](slack_vs_source.md). |
| Parallel island solving | ✅ | AC islands of at least `power_flow.islands.parallel_min_buses` buses solve concurrently on Julia threads (`runtime.parallel.enabled`), bitwise identical to serial. |
| Voltage-dependent prosumer control (`Q(U)`, `P(U)`) | ✅ | Controller-aware mismatch and Jacobian terms in the rectangular formulation. |
| Wrong-branch detection (`wrong_branch_detection`) | ✅ | Post-convergence plausibility guard (`off|warn|fail`) for low-voltage, non-finite and large-angle solutions, reported in `ACPFlowReport.metadata`, the console summary and the Web UI. |
| Jacobian condition diagnostics (`condestJacobian`, `reportCondition`) | ✅ | Hager/Higham 1-norm estimate on the LU factorization with the attainable accuracy. [Solver](solver.md). |
| Narrative diagnostics and self-check (`diagnose.log`, `run_fixed_reference_self_check`) | ✅ | Worst-mismatch bus, mismatch trend, autodamp health, branch anomalies; the self-check evaluates the case's stored operating point with every start machine off. |
| External solver interface (`PFModel` / `PFSolution`) | ✅ | Integration point for solver backends outside the built-in rectangular NR. |

## Short circuit (IEC 60909-0)

`runShortCircuit!` has its own result type and needs no converged power flow.

| Feature | Status | Notes |
|---|:---:|---|
| Initial symmetrical current per fault bus | ✅ | `Ik''` (max/min case), `Sk''` and peak current `i_p`; positive sequence, series impedances only, IEC Table-1 voltage factors with `short_circuit.c_factor` as override. [Short-Circuit Analysis](short_circuit.md). |
| Data sources | ✅ | The CGMES harvest (machines, feeders, motors), PowSyBl generators with `generatorShortCircuit`, the feeder records of `addExternalGrid!`. A machine without data enters on the default reactance and flags the rows it feeds; the service reports `succeeded` only when every source carries its data. |
| Safety-flag contract | ✅ | Substituted defaults and skipped contributions are flagged on the affected rows; a flagged maximum is a lower bound. |
| All-bus sweeps on Julia threads | ✅ | The fault-bus list fans out over task chunks (`runtime.parallel.*`), row-identical to the serial sweep. |
| Takahashi sparse inverse (`short_circuit.sweep_method`, default `auto`) | ✅ | The Thevenin diagonal of an island from one selected-inverse pass over the LU factors; `auto` applies it from `takahashi_min_buses` (default 50) buses, agreeing with `:solves` to machine precision. |
| Web UI integration | ✅ | **Short circuit** runs both cases without a power-flow solve and writes `short_circuit_max.csv`/`short_circuit_min.csv`. |
| Unbalanced faults, `K_T`/`K_G` corrections | ❌ | Balanced positive-sequence faults only, no impedance correction factors. Offered instead: the balanced study with c-factor selection, kappa and `i_p`, and per-row safety flags. |

## N-1 contingency analysis

| Feature | Status | Notes |
|---|:---:|---|
| Branch-outage batches (`runContingencies!`) | ✅ | Warm-started cases on template copies of the solved base case, checked against the voltage band and `sn_MVA` loadings. [N-1 Contingency Analysis](contingency.md). |
| Case sources | ✅ | `generateN1Branches` (in-service branches, transformer filter, parallel circuits as `name#branchIdx`) and imported MATPOWER FOR001 lists; filters `min_vn_kV` / `min_sn_MVA` / `name_pattern`. |
| Generator outages | ✅ | `generateN1Generators` (`kind = :gen`, filters `min_pg_MW` / `name_pattern`). The slack absorbs the loss or `distributed_slack_enabled = true` shares it; `auto_slack` (the default of the batch entries, the service and the Web UI) promotes a survivor when the slack unit is the outage, and an island that lost its reference takes its best unit by stated priority. |
| Case weights | ✅ | `weight` per case (default 1.0) for a severity-weighted ranking; `readContingencyWeightsCSV` + `applyContingencyWeights` attach outage rates. |
| Failure semantics | ✅ | Islanding without a promotable reference, non-convergence and unresolvable elements are reported per case, not thrown; `contingency.rescue_ladder` tries several starts; a load-only island is a quantified load-shed result. |
| Parallel execution | ✅ | The case list fans out over Julia threads (`runtime.parallel.*`), identical to the serial run. |
| Output | ✅ | `printContingencyResults` (severity-ranked, failures first) and `writeContingencyResultsCSV`, with `OverloadRecord` loadings, shed load, weight and severity per case. |
| Aggregate report | ✅ | `buildContingencyReport` folds a batch into a `ContingencyReport`; `printContingencyReport` prints it. |
| Web UI | ✅ | **Contingency (N-1)** with a branch/generator selector runs the batch through the shared service path, writes `contingency_n1.csv` and a `run.log` report and shows an outcome summary; a slack-unit outage is named, not shown as a failure. |

## State estimation (WLS)

`runse!` reconstructs the state from redundant, noisy measurements, sharing the network model and importers with the power flow.

| Feature | Status | Notes |
|---|:---:|---|
| Nonlinear WLS estimator (`runse!`) | ✅ | Iterative WLS on the shared `Net`; `state_estimation.update_net = true` writes the estimate back. [State Estimation](state_estimation.md). |
| SCADA-style measurements | ✅ | `Vm`, `Pinj`, `Qinj`, `Pflow`, `Qflow` with helper builders. |
| Branch current magnitudes (`ImagMeas`) | ✅ | Ampere-valued auxiliary measurements (`addImagMeasurement!`) at branch ends or shunt bays, excluded from observability and gated by `state_estimation.imag_activation_iteration`; they raise the localizability (`wii`) of the power measurements. |
| Shunt susceptance estimation (`setShuntEstimation!`) | ✅ | Case A: a released shunt's B becomes a state (`SEResult.shuntEstimates`, opt-in write-back `updateShunts`). Case B: `deriveShuntPseudoMeasurements!` derives protected `SHDERIV` pseudo-measurements from bay current and voltage. Injection-mode shunts are rejected. |
| PMU voltage phasors | ✅ | Magnitude as weighted `VmMeas`, angle as `VaMeas` with an estimated reference-angle offset (`state_estimation.pmu_ref_offset`); helper `addPmuPhasorMeasurement!`. |
| PMU current phasors (`IaMeas`) | ✅ | Branch-end current angles sharing the PMU offset, gated near zero current (`state_estimation.ia_current_floor_A`); helper `addCurrentPhasorMeasurement!`. |
| Transformer tap estimation (`setTapEstimation!`) | ✅ | A released tap (ratio, phase, or both) becomes a state; after convergence it is fixed to the nearest mechanical step and a final run reports J before versus after (`SEResult.tapEstimates`/`tapFixation`). Machine transformers are skipped, unobservable columns frozen. Opt-in write-back `state_estimation.update_taps`; machine taps are back-calculated (`calcMachineTrafoTapFromSE`). |
| Zero-injection buses | ✅ | Tightly weighted pseudo-measurements, protected from elimination and robust down-weighting. Hard equality constraints (Hachtel/Lagrange) are not implemented. |
| Observability analysis | ✅ | Two-stage global check (structural islands, then numeric rank with `state_estimation.rank_tol_factor`), structural matching, local column checks. [Observability](observability.md). |
| Critical and nearly critical measurements | ✅ | From the residual sensitivity diagonal in one selected-inverse pass (`state_estimation.criticality_method = omega`, default): `Omega_ii = 0` critical, `wii < 0.3` nearly critical. [Observability](observability.md). |
| Observable-island identification | ❌ | No decomposition into maximal observable islands. Offered instead: the global check reports rank deficiency and the structural islands, and `unobservable_state_columns` names every state the set does not pin down. |
| Observability restoration (automatic pseudo-measurements) | ❌ | No automatic pseudo-measurement selection. Offered instead: zero-injection pseudo-measurements on request and the dark-state list of the global check. |
| Optimal measurement / PMU placement | ❌ | No placement optimizer. Offered instead: local observability on a chosen state subset (`evaluate_local_observability`) and critical-measurement reporting. |
| Topology validation (`validate_topology`) | ✅ | Three advisory stages: linear pre-checks, the suspected-station classification when elimination exhausts against a `:high` band, and the status-toggle hypothesis test (`test_topology_hypotheses`). Nothing is switched automatically. |
| Bad-data diagnostics | ✅ | Wilson-Hilferty band test, residual ranking with `wii` and a localizable flag, sequential elimination with trace (`runse_diagnostics`), optional residual-correlation report, `summarize_se_diagnostics`, `print_se_diagnostics`. |
| Robust R modification (`state_estimation.robust`) | ✅ | Two-stage weight modification (tangential above 3 sigma, suppression above 6 sigma) with statistics on the original sigmas (`SEResult.robustRows`). |
| Synthetic measurements from a PF result | ✅ | `setMeasurementsFromPF!` builds test sets from a solved power flow, with or without noise. |
| Flat start control | ✅ | Same start discipline as the power flow. |

## Controllers and FACTS (power-flow outer loop)

All controllers run in the outer loop above `runpf!`, with results in `ControlRunResult` / `latest_control_result(net)` and the `controllableElements` view; device taxonomy and limit characteristics are on [FACTS Devices](facts.md).

| Feature | Status | Notes |
|---|:---:|---|
| Transformer regulation (OLTC / PST / combined) | ✅ | Voltage, branch-active-power and combined modes, discrete steps with limits, split regulation; per capability in the [transformer support table](@ref transformer-support). |
| Remote voltage control by machines (`MachineVoltageControl`) | ✅ | A PQ machine holds the voltage of another bus via its reactive output (secant outer loop, `at_limit` on the bounds); wired from CGMES `RegulatingControl`s with `cgmes_import.machine_control`. [Remote Voltage Control](remote_voltage_control.md). |
| STATCOM current-based limit mode (`s_max_mva`) | ✅ | The machine controller as VSC shunt compensator: `Q_lim = V * S_max`, re-evaluated every outer iteration. |
| SVC variable-shunt voltage control (`ShuntVoltageControl`) | ✅ | A continuous shunt susceptance holds the local voltage; at a limit it clamps and Q follows V². Imported CGMES SVCs still map as static PV injections. |
| MSC/MSR switched shunt banks (`step_mvar`) | ✅ | The shunt controller in whole blocks, truncated toward the target, parking on the last step before crossing (`status = :parked`). |
| TCSC series-reactance flow control (`SeriesReactanceControl`) | ✅ | The series reactance of a line is the actuator, the branch active power the target; `at_limit` at a range end. Transformer branches are rejected. [Series Compensation (TCSC)](series_compensation.md). |
| SSSC injected-voltage limit mode (`v_inj_max_pu`) | ✅ | The series controller as VSC compensator: the reactance window is bounded by the injectable series voltage and shrinks with loading. |
| UPFC (combined shunt + series converter) | ✅ | `addUpfcControl!` (YAML type `upfc`): `model = :quadrature` (default, SSSC + STATCOM composite) or `model = :full` (arbitrary-phase series injection steering line P and Q with the DC-link balance on the shunt; no explicit series current limit, converges for moderate targets). Scope: stationary model, no IPFC. |
| Equipment impedance vs FACTS operating point (`r_base_pu`/`x_base_pu`) | ✅ | A series-FACTS run stamps its compensated impedance onto the live branch; `runShortCircuit!` and the exports read the physical base. `restoreBaseImpedances!`, `clearUpfcFullControllers!`/`clearSeriesReactanceControllers!` reset the branch. |
| HVDC back-to-back pairing (`HvdcPairControl`) | ✅ | Two converter injections coupled by `P_to = P_transfer - loss`, per-terminal fixed Q or voltage-target secant, optional transfer rating; HVDC-joined areas stay separate islands. Opt-in from both importers (`matpower_dcline_mode` / `hvdc_mode` `= paired_control`); `mode = :island_feed` feeds an island fed only by the receiving converter. [HVDC Back-to-Back](hvdc_back_to_back.md). |
| YAML controller instantiation | ✅ | Named definitions under `control.controllers`. [Control Framework](control_framework.md). |

### [Transformer support](@id transformer-support)

A 3-winding transformer is a star equivalent with an internal AUX bus
(`add3WTPiModelTrafo!` / `create3WTWindings!`); tap and phase changers sit on
a chosen winding (`tap_side`, `phase_tap_side`) and are controlled per leg
branch. Tap and voltage control act in the PF outer loop only; SE uses the
same transformer model without controllers.

| Type | 2-winding | 3-winding | Remarks |
|---|:---:|:---:|---|
| Fixed-tap transformer | ✅ | ✅ | Base network model, PF and SE. |
| OLTC (voltage control on ratio tap) | ✅ | ✅ | `addTapController!` with `mode = :voltage`; discrete steps with tap limits; remote target bus via `target_bus`. |
| PST (active-power control on phase tap) | ✅ | ✅ | `mode = :branch_active_power`; discrete steps with phase limits. |
| Combined regulation, single controller | ✅ | ✅ | `mode = :voltage_and_branch_active_power`: one controller drives ratio and phase taps together. |
| Split combined regulation (Schrägregler) | ✅ | ✅ | Two independent controllers on one unit (voltage on the ratio tap, active power on the phase tap), each with its own target, deadband and status. Demo: `tap_control_schraeg_two_controllers.jl`. |
| Symmetrical PST model (CGMES) | ✅ | ✅ | A typed `PhaseTapChangerModel` on the winding drives the branch at its step ([Branch model](branchmodel.md)); a tap controller moves the step. Persisted in the SCF `tap_changer_models` block, exported as `PhaseTapChangerSymmetrical`. |
| Asymmetrical PST model (CGMES) | ✅ | ✅ | Includes the quadrature booster (ψ = 90°); same persistence and export (`PhaseTapChangerAsymmetrical`). The DTF regulators (Laengs-, Schraeg-, Querregler) are this kind (`tap_changer_kind`). |
| Tabular PST model (`:tabular`) | ✅ | ✅ | `TapTablePoint` lookup; CGMES `PhaseTapChangerTabular` and `PhaseTapChangerLinear` import through it (validated against RealGrid's SV state) and export as `PhaseTapChangerTabular`. 3-winding via `create3WTWindings!` with `phase_taps`; no real-case validation of a 3WT table yet. |
| Coordinated master/slave voltage control | ✅ | ✅ | Parallel transformers regulate as a group (`followers` on `addPowerTransformerControl!`): the master runs the discrete loop, followers mirror step-synchronously. CGMES deliveries with several tap changers on one `TapChangerControl` import as one group. No participation-factor allocation. [Control Framework](control_framework.md). |
| Tap-changer impedance model (`tap_changer_model`) | ✅ | ✅ | `ideal` (default) or `impedance_correction` (R/X re-referred through the tapped winding), for all transformers of an imported case (MATPOWER and native DTF importer), via `calcTapCorrectedRX`/`calcTapImpedanceCorrectionFactor`. SE inherits the choice of the imported `Net`. |

!!! note "The typed tap-changer models"
    Ratio (`PowerTransformerTaps`) and phase (`PhaseTapChangerModel`, kinds `:symmetrical`, `:asymmetrical`, `:tabular`) models sit on a transformer winding, directly (`addPIModelTrafo!` with `taps`/`phase_taps`) or via `create3WTWindings!`. A winding that carries a model drives its branch, controllers move the model step, the SCF file persists the model and the CGMES export writes its class. Per-transformer settings live in the SCF `sparlectra` block, not in the YAML configuration; the PowSyBl step tables are dumped, not yet built into a `:tabular` model.

## Network model and data exchange

| Feature | Status | Notes |
|---|:---:|---|
| Shared network model | ✅ | Buses, π-equivalent branches, transformers, generators, loads, shunts and links in one `Net` used by every analysis domain. |
| Asymmetric branch shunt | ✅ | A charging admittance per terminal (`g_from_pu`, `b_from_pu`, `g_to_pu`, `b_to_pu`); a transformer's magnetizing admittance stays on its own end, a PowSyBl line keeps `b1` and `b2`; SCF persists the split (`branch_shunt_split`), the MATPOWER export writes the excess as a named bus shunt. [Branch model](branchmodel.md). |
| Per-terminal branch status | ✅ | A branch open at one terminal stays in the model as its exact pi reduction; CGMES `Terminal.connected` maps both directions, the MATPOWER export substitutes the exact `Y_in` bus shunt. |
| Topological bus links (`addLink!`) | ✅ | Impedance-less couplers and sectionalizers: closed-link clusters contract onto one bus before the Y-bus is built, link flows are reconstructed after the solve (`calcLinkFlowsKCL!`). Imported from CGMES switches and the MATPOWER block `mpc.sparlectra.links`. [Links](links.md). |
| Links in state estimation | ✅ | SE runs on the contracted net; cluster injections aggregate only as a whole (`LINKAGG`), link flow measurements are allocation inputs (`calcLinkFlowsSE!`), never WLS rows. |
| FACTS SE view (`se_view`) | ✅ | Frozen-operating-point report: `runse!` never invokes outer-loop control; `se_view`/`print_se_view` list frozen controllers, shunt releases, link clusters and excluded measurements. |
| Configurable bus-shunt modeling | ✅ | `model.bus_shunt_model = "admittance"` (default) stamps bus shunts into the Y-bus; `"voltage_dependent_injection"` keeps them in the mismatch terms. Scope: the injection variant is a MATPOWER import option and refuses active link merges by name. |
| CGMES import (`importCGMES` / `createNetFromCGMES`) | ✅ | CGMES 2.4.15 bus-branch import (EQ+SSH+TP+SV, boundary sets, folders/ZIP/ZIP-in-ZIP) with transformers and their tap controllers, remote voltage controllers, retained switches as links, the short-circuit harvest, `summarizeCGMES` and `compareWithSV`. CGMES 3.0 deliveries are read as well; converters map as fixed PCC injections. Validated on MicroGrid, SmallGrid, FullGrid and the 6209-bus RealGrid. [CGMES Import](cgmes_import.md). |
| Node-breaker topology processor | ✅ | Imports node-breaker deliveries without a TP profile (connectivity nodes aggregate across closed non-retained switches, retained switches stay couplers); runs only when no non-boundary `TopologicalNode` exists; verified against the shipped TP on the MiniGrid/SmallGrid/FullGrid sets. |
| CGMES export (`writeCGMESFiles`) | ✅ | Complete CGMES 2.4.15 delivery (EQ + TP + SSH + SV, optionally one ZIP) with roundtrip-stable identity; a re-imported net solves to the same power flow and short circuit. [CGMES Export](cgmes_export.md). |
| CGMES import analysis (`analyzeCGMES`) | ✅ | Explains a non-importable delivery: supplied models, declared prerequisites, unresolved-reference histogram, verdict. |
| PowSyBl IIDM import (`PowsyblAdapter`, `read_iidm_tables`) | ✅ | XIIDM networks (schema 1.0 to 1.17; compressed files are not read) read in Julia: node-breaker topology, transformers at their current tap, three-winding star buses, HVDC stations as fixed injections, SVC, dangling and tie lines, operational and reactive limits; reproduces the OpenLoadFlow solution of five example networks. [PowSyBl Import](powsybl_import.md). |
| Base-voltage inference (`cgmes_import.infer_base_voltages`) | ✅ | Reconstructs missing nominal voltages from the SV state and transformer rated voltages. |
| MATPOWER import / export | ✅ | Configurable SHIFT and TAP conventions, transformer-loss metadata round trips, auto-profile recommendations, `writeMatpowerCasefile` with optional solved-state columns; generator costs are kept and written back unchanged ([Generator costs](@ref matpower_gencost)). [MATPOWER](matpower.md). |
| Native DTF import | ✅ | Native `.dat` network cases incl. FOR001/FOR002 validation workflows. [DTF Format](dtf_format.md). |
| Synthetic tiled-grid generator | ✅ | `build_synthetic_tiled_grid_net` creates artificial one-voltage-level benchmark networks. |

## Workflow, reporting, and tooling

| Feature | Status | Notes |
|---|:---:|---|
| Framework workflow (`run_sparlectra`) | ✅ | Configuration-driven import/control/solve/output, one `SparlectraRunResult` per run; `run_sparlectra_cases` runs configured MATPOWER batches. |
| Central typed configuration | ✅ | `SparlectraConfig` with cached YAML loading, typed validation and override precedence. [Configuration](configuration.md). |
| Parallel runtime (`runtime.parallel.*`) | ✅ | One switch set (`enabled`, `max_tasks`, `min_work_items`) gates every threaded surface; serial fallbacks are the same functions. [Parallel Execution](parallel_execution.md). |
| Machine-readable report (`ACPFlowReport`) | ✅ | DataFrame-friendly rows for buses, branches, links, transformer controls, Q-limit events and HVDC links. |
| GUI-ready programmatic run API | ✅ | `run_sparlectra_api` with stable run IDs, schema-versioned status, controlled overrides and artifact discovery. [Programmatic API](programmatic_api.md). |
| Local PowerFlow service boundary | ✅ | `start_powerflow_run`, persistent run index, restart recovery, result lookup and safe artifact resolution, no HTTP dependency. [Local PowerFlow Service](powerflow_service.md). |
| Local browser Web UI | ✅ | Case management, run forms with contextual help, run history, artifact viewing and the state-estimation section. Scope: loopback only, no authentication, no public deployment mode. [Web UI](webui.md). |

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

# Feature Matrix (Quick Overview)

Sparlectra covers four analysis domains: power flow (including
quasi-static time series), balanced short circuit, N-1 contingency
analysis and state estimation. Each has its own
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
| Start projection and start modes | ✅ | DC-angle, blend-scan and (on a flat start) transformer-ratio candidates, imported-profile starts, optional guarded current-iteration pre-solve (`power_flow.start_current_iteration`). A flat start puts every angle at the slack angle; `power_flow.start_mode.rotate_behind_pst` (default off) starts buses behind phase shifters at the shifted angle, and a solve without it that ends on a low-voltage root behind a shifter warns and names the key. [Start strategies](start_strategies.md). |
| Distributed active-power slack (`power_flow.distributed_slack`) | ✅ | The island's P imbalance is shared over participating generators; weight modes `pg_weighted` (default), `pmax_weighted`, `headroom_weighted`, `imported`, `explicit`. [Solver Guide](solver.md). |
| Automatic slack selection (`power_flow.auto_slack`) | ✅ | Promotes the best injection to slack when a case registers no reference; off by default. API: `ensureSlack!`. |
| Reference priority (`referencePriority`) | ✅ | Integer per unit (1 strongest, 0 none) ordering the units that may take over a reference; without priorities the strongest unit wins (`reference_candidate_rank`). Imported from MATPOWER type 3 and CGMES, carried by SCF and the CGMES export. [Reference priority](slack_vs_source.md#Reference-priority). |
| External grid element (`addExternalGrid!`) | ✅ | IEC 60909-0 network feeder: ideal slack by default, optional `internal_impedance = true`; its short-circuit data feeds `runShortCircuit!`. [Slack Bus and External Grid Sources](slack_vs_source.md). |
| Parallel island solving | ✅ | AC islands of at least `power_flow.islands.parallel_min_buses` buses solve concurrently on Julia threads (`runtime.parallel.enabled`), bitwise identical to serial. |
| Voltage-dependent prosumer control (`Q(U)`, `P(U)`) | ✅ | Controller-aware mismatch and Jacobian terms in the rectangular formulation. |
| Wrong-branch detection (`wrong_branch_detection`) | ✅ | Post-convergence plausibility guard (`off\|warn\|fail`) for low-voltage, non-finite and large-angle solutions, reported in `ACPFlowReport.metadata`, the console summary and the Web UI. |
| Jacobian condition diagnostics (`condestJacobian`, `reportCondition`) | ✅ | Hager/Higham 1-norm estimate on the LU factorization with the attainable accuracy. [Solver](solver.md). |
| Narrative diagnostics and self-check (`diagnose.log`, `run_fixed_reference_self_check`) | ✅ | Worst-mismatch bus, mismatch trend, autodamp health, branch anomalies; the self-check evaluates the case's stored operating point with every start machine off. |
| External solver interface (`PFModel` / `PFSolution`) | ✅ | Integration point for solver backends outside the built-in rectangular NR. |
| Quasi-static time series (`runQsts`, `qsts.mode`, default `off`) | ✅ | One power flow per step on an SCF time series or a series of CGMES snapshots: `snapshots` (independent steps, parallel) or `sequential` (SCF series, tap and bank states carried over), with event lists, switching counts and limit violations. Scope: no dynamics between steps. [QSTS](qsts.md). |

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
| Converter-connected units (IEC 60909-0:2016 current source) | ◐ | CGMES 3.0 `PowerElectronicsConnection` units (mapped with `cgmes_import.map_power_electronics`) feed the maximum case as current sources, `maxIFault` times the rated current from `ratedS` and `ratedU`, superposed on the equivalent voltage source. A unit without `maxIFault` or `ratedS` flags its island, so the maximum there reads as a lower bound. `maxIFault` is a CIM17 attribute in no CGMES 3.0 profile, so a conforming delivery has none; `short_circuit.converter_current_factor` (default off) then sets the current as a multiple of the rated current, otherwise every converter takes the flag. A synchronous machine on a wind or solar unit (CGMES `WindGeneratingUnit`/`SolarGeneratingUnit`, IIDM `energySource` `WIND`/`SOLAR`) takes the same path with the factor, rated current from `ratedS`, else from the maximum active power of the unit; without the factor or a size it flags its island. An island fed by converters alone reports the arithmetic sum of their currents, flagged as an upper bound. |

## N-1 contingency analysis

| Feature | Status | Notes |
|---|:---:|---|
| Branch-outage batches (`runContingencies!`) | ✅ | Warm-started cases on a working copy of the solved base case, checked against the voltage band and `sn_MVA` loadings. [N-1 Contingency Analysis](contingency.md). |
| Correlated events | ✅ | Double circuits, busbar sections, stations, couplers and file name groups as one case each (`generateCorrelatedEvents`), on the same path and screening as N-1; `contingency.events.*` in the service (SCF cases). |
| N-2 pairs | ✅ | `generateBranchPairs`: every pair of closed branches, with a cap (`contingency.pairs.max_pairs`); no pruning rule (the tested rules missed harmful pairs). |
| Display names | ✅ | Rows, elements and overloads carry the file names of an SCF case (`contingency_names`); the network keeps its generated names. |
| Case sources | ✅ | `generateN1Branches` (in-service branches, transformer filter, parallel circuits as `name#branchIdx`) and imported MATPOWER FOR001 lists; filters `min_vn_kV` / `min_sn_MVA` / `name_pattern`. |
| Generator outages | ✅ | `generateN1Generators` (`kind = :gen`, filters `min_pg_MW` / `name_pattern`). The slack absorbs the loss or `distributed_slack_enabled = true` shares it; `auto_slack` (the default of the batch entries, the service and the Web UI) promotes a survivor when the slack unit is the outage, and an island that lost its reference takes its best unit by stated priority. |
| Case weights | ✅ | `weight` per case (default 1.0) for a severity-weighted ranking; `readContingencyWeightsCSV` + `applyContingencyWeights` attach outage rates. |
| Failure semantics | ✅ | Islanding without a promotable reference, non-convergence and unresolvable elements are reported per case, not thrown; `contingency.rescue_ladder` tries several starts; a load-only island is a quantified load-shed result. |
| Parallel execution | ✅ | The case list fans out over Julia threads (`runtime.parallel.*`), identical to the serial run. |
| Output | ✅ | `printContingencyResults` (severity-ranked, failures first) and `writeContingencyResultsCSV`, with the case kind, the outaged elements, `OverloadRecord` loadings, shed load, weight and severity per case. |
| Aggregate report | ✅ | `buildContingencyReport` folds a batch into a `ContingencyReport` (counts by outcome and by kind); `printContingencyReport` prints it. |
| Web UI | ✅ | **Contingency (N-1)** with a branch/generator selector runs the batch through the shared service path, writes `contingency_n1.csv` and a `run.log` report and shows an outcome summary; a slack-unit outage is named, not shown as a failure. |

## State estimation (WLS)

`runse!` reconstructs the state from redundant, noisy measurements, sharing the network model and importers with the power flow.

| Feature | Status | Notes |
|---|:---:|---|
| Nonlinear WLS estimator (`runse!`) | ✅ | Iterative WLS on the shared `Net`; `state_estimation.update_net = true` writes the estimate back. [State Estimation](state_estimation.md). |
| SCADA-style measurements | ✅ | `Vm`, `Pinj`, `Qinj`, `Pflow`, `Qflow` with helper builders. |
| Measurement sigma by voltage level | ✅ | `state_estimation.sigma_by_voltage_level`: sigmas per level and measurement type in MW, MVAr, A and kV, optional share of the value; a row sigma wins, then the level, then the global default. Used by the measurement table import and its generators, editable on the State Estimation section of the Web UI. [Measurement sigma by voltage level](@ref se-sigma-by-level). |
| Branch current magnitudes (`ImagMeas`) | ✅ | Ampere-valued auxiliary measurements (`addImagMeasurement!`) at branch ends or shunt bays, excluded from observability and gated by `state_estimation.imag_activation_iteration`; they raise the localizability (`wii`) of the power measurements. |
| Shunt susceptance estimation (`setShuntEstimation!`) | ✅ | Case A: a released shunt's B becomes a state (`SEResult.shuntEstimates`, opt-in write-back `updateShunts`). Case B: `deriveShuntPseudoMeasurements!` derives protected `SHDERIV` pseudo-measurements from bay current and voltage. Injection-mode shunts are rejected. |
| PMU voltage phasors | ✅ | Magnitude as weighted `VmMeas`, angle as `VaMeas` with an estimated reference-angle offset (`state_estimation.pmu_ref_offset`); helper `addPmuPhasorMeasurement!`. |
| PMU current phasors (`IaMeas`) | ✅ | Branch-end current angles sharing the PMU offset, gated near zero current (`state_estimation.ia_current_floor_A`); helper `addCurrentPhasorMeasurement!`. |
| Transformer tap estimation (`setTapEstimation!`) | ✅ | Every released tap changer (ratio, asymmetrical incl. quadrature, symmetrical, tabular, linear; at most two per transformer) becomes one state, its step on the additional voltage, from the changer data of the CGMES, PowSyBl, SCF, MATPOWER or DTF import. After convergence each step is fixed to the nearest physical step on the absolute band and a final run reports J before versus after (`SEResult.tapEstimates`/`tapFixation`, one row per changer). Per-step impedance data (a CGMES tap table, `xMin`/`xMax`, `TapStepImpedance`, the typed `x_min`/`x_max` curve of an SCF file) always move the impedance with the step. `state_estimation.tap_estimation.impedance_correction` (`auto`, default, or `off`) only switches the analytic correction of a transformer without such data. Each row names the source of its reactance (`reactance_source`: `table`, `analytic`, `nominal`). Single transformers can be excluded (`state_estimation.tap_estimation.exclude`, by name or mRID; Web UI: the Case page), machine transformers are excluded by default and their step is calculated after the run (`include_machine_transformers`); a changer the measurements cannot determine (rank test at the start and again at the converged state, where a phase changer whose mirror step fits as well is withdrawn) is not estimated and keeps its transferred step, with its reason and a notice; a changer without a transferred step starts at its middle position. Switch `state_estimation.tap_estimation.enabled`; a run that does not converge with tap estimation, or whose fixed steps fit worse than the transferred steps (SE service, `worse_fit`), falls back to the run without it (`state_estimation.tap_estimation.fallback`). Opt-in write-back `state_estimation.update_taps`; machine taps are back-calculated (`runse!`, `calcMachineTrafoTapFromSE`). |
| Zero-injection buses | ✅ | Tightly weighted pseudo-measurements at passive buses (`state_estimation.zero_injection = on`, default, sigma `state_estimation.zero_injection_sigma_mva`), protected from elimination and from the `replacement` down-weighting and from the `staged` widening. |
| Zero injection as equality constraint | ✅ | `state_estimation.zero_injection = constraint`: the zero injections, also that of a fused coupler bus whose members are all passive (`LINKAGG`), are equality constraints of the WLS problem, solved as a scaled KKT system. They hold exactly, stay out of J, are never ranked as suspicious, and avoid the ill conditioning of a tiny sigma (`on` weights them with `zero_injection_sigma_mva`). The run reports the Lagrange multiplier per bus (`se_zero_injection_multipliers.csv`) and the condition number of the KKT or gain matrix. [Zero injections and coupler balances](@ref se_zero_injection). |
| Observability analysis | ✅ | Two-stage global check (structural islands, then numeric rank with `state_estimation.rank_tol_factor`), structural matching, local column checks. [Observability](observability.md). |
| Critical and nearly critical measurements | ✅ | From the residual sensitivity diagonal in one selected-inverse pass (`state_estimation.criticality_method = omega`, default): `Omega_ii = 0` critical, `wii < 0.3` nearly critical. [Observability](observability.md). |
| Topology validation (`validate_topology`) | ✅ | Three advisory stages: linear pre-checks, the suspected-station classification when elimination exhausts against a `:high` band, and the status-toggle hypothesis test (`test_topology_hypotheses`). Nothing is switched automatically. |
| Bad-data diagnostics | ✅ | Range check of the voltage measurements before the WLS (0 kV at a de-energized bus leaves the estimation, an implausible value at an energized bus becomes a pseudo measurement), Wilson-Hilferty band test, residual ranking with `wii` and a localizable flag, sequential elimination with trace (`runse_diagnostics`), optional residual-correlation report, `summarize_se_diagnostics`, `print_se_diagnostics`. |
| Robust R modification (`state_estimation.robust_mode`) | ✅ | `staged`: two-stage weight modification (tangential above `robust_k1`, suppression above `robust_k2`, default 3 and 6 sigma); `replacement`: rows above `k_suppress` solved at a fixed wide sigma, judged on the converged state. Statistics stay on the original sigmas (`SEResult.robustRows`). [Bad-data thresholds and robust modes](@ref se-bad-data). |
| Synthetic measurements from a PF result | ✅ | `generateMeasurementsFromPF` and `setMeasurementsFromPF!` build test sets from a solved power flow, with or without noise. |
| Measurement table | ✅ | The measurement format: CSV with one set per timestamp, keyed by CGMES mRID or by element name, carrying voltages and voltage angles, bus balances, flows, currents, current angles and shunt values; import report per row; `readMeasurementTable!`; generators from a solved power flow (`generateMeasurementTableFromPF`, with coverage) and from a CGMES SV profile; pseudo rows from the scheduled values for unobserved loads and generators (`generatePseudoMeasurementTable`). The Web UI generator writes it, the shipped sets use it. New format, may change before 1.0. [Measurement Table](state_estimation_measurement_table.md). |
| Measurements in the SCF case file | ✅ | PGM sensors in SI plus `sparlectra.measurements.rows` with the semantics of the measurement table (table types and units, node and terminal keys, shunt bay), since SCF format 1.2. The current format is 1.3; files from 1.0 on are read. [Measurements in SCF](scf.md#Measurements). |
| Flat start control (`state_estimation.flatstart`, default on) | ✅ | Flat start (every \|V\| 1 pu, every angle at the slack angle) or a warm start from the stored state, independent of the power-flow start (`power_flow.flatstart`, default off); results and the write-back keep the stored slack angle. `state_estimation.start_mode.rotate_behind_pst` turns buses behind phase shifters like the power flow flat start. [Start value](@ref se_start_value). |

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
| Controller priority (`priority`) | ✅ | An integer rank per controller: a controller with a larger number acts only when every controller with a smaller one is inside its deadband or at its limit, so a switched bank can stand behind the tap changer of its bus. Default 0 for all, the loop without ranks. [Control Framework](@ref controller_priority_ranks). |
| Tap control range (`tap_min`/`tap_max` per controller) | ✅ | A tap controller moves the ratio only inside its own range within the physical one and stands at its limit at either end, so the bank behind it takes over earlier; a ratio outside the range goes onto its nearer end in the first pass with a warning naming the transformer, its start and the limit. Set in `control.controllers` or on `addPowerTransformerControl!`; a range outside the physical one fails naming the entry. [Control Framework](control_framework.md). |

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
| Tap-changer impedance model (`model.tap_changer_model`) | ✅ | ✅ | `ideal` (default) or `impedance_correction` (R/X re-referred through the tapped winding, `calcTapCorrectedRX`/`calcTapImpedanceCorrectionFactor`), for all transformers of a MATPOWER or DTF import. An SCF file carries the choice per transformer in `tap_changer_models`. The correction is defined for the reciprocal convention; the SCF import warns for an asymmetrical changer with the direct convention. Per-step impedance data always take precedence. In the tap estimation, `state_estimation.tap_estimation.impedance_correction = off` switches the analytic correction off. |

!!! note "The typed tap-changer models"
    Ratio (`PowerTransformerTaps`) and phase (`PhaseTapChangerModel`, kinds `:symmetrical`, `:asymmetrical`, `:tabular`) models sit on a transformer winding, directly (`addPIModelTrafo!` with `taps`/`phase_taps`) or via `create3WTWindings!`. A winding that carries a model drives its branch, controllers move the model step, the SCF file persists the model and the CGMES export writes its class. Per-transformer settings live in the SCF `sparlectra` block, not in the YAML configuration. The PowSyBl step tables feed the tap estimation as table changers; power flow, controllers and export keep the ratio, shift and impedance of the current step.

## Network model and data exchange

| Feature | Status | Notes |
|---|:---:|---|
| Shared network model | ✅ | Buses, π-equivalent branches, transformers, generators, loads, shunts and links in one `Net` used by every analysis domain. |
| Line input in physical units | ✅ | `addACLine!` takes the line totals in ohm and siemens, `addACLinePerKm!` values per km with the length (b from the capacitance at the system frequency of the net), `addPIModelACLine!` per unit. [Component Types](component_types.md). |
| System frequency (50/60 Hz) | ✅ | `net.nominal_frequency_hz`, set from an SCF file's `meta.f_nom`, otherwise from `model.nominal_frequency_hz` (default 50 Hz, selectable in the Web UI). Enters only where a capacitance becomes a susceptance (`addACLinePerKm!`, SCF/PGM `c1`); power flow of cases in ohm, siemens or per unit and the short circuit do not depend on it. [System frequency](@ref system_frequency). |
| Asymmetric branch shunt | ✅ | A charging admittance per terminal (`g_from_pu`, `b_from_pu`, `g_to_pu`, `b_to_pu`); a transformer's magnetizing admittance stays on its own end in the PowSyBl import (the CGMES import puts half of the total on each terminal, as the ENTSO-E SV does), a PowSyBl line keeps `b1` and `b2`; SCF persists the split (`branch_shunt_split`), the MATPOWER and CGMES exports write the excess as a bus shunt. [Branch model](branchmodel.md). |
| Per-terminal branch status | ✅ | A branch open at one terminal stays in the model as its exact pi reduction; CGMES `Terminal.connected` maps both directions, the MATPOWER export substitutes the exact `Y_in` bus shunt. |
| Topological bus links (`addLink!`) | ✅ | Impedance-less couplers and sectionalizers: closed-link clusters contract onto one bus before the Y-bus is built, link flows are reconstructed after the solve (`calcLinkFlowsKCL!`). Imported from CGMES switches and the MATPOWER block `mpc.sparlectra.links`. [Links](links.md). |
| Links in state estimation | ✅ | SE runs on the contracted net; cluster injections aggregate only as a whole (`LINKAGG`). Link flow measurements are WLS rows (flow from the nodal balances of the cluster, circulation states on coupler rings), get normalized residuals and are eliminated like branch flows; `calcLinkFlowsSE!` gives the link flows after the estimate. The measurement table carries them as `P`/`Q` on the link terminals. [State Estimation Extensions](state_estimation_extensions.md). |
| FACTS SE view (`se_view`) | ✅ | Frozen-operating-point report: `runse!` never invokes outer-loop control; `se_view`/`print_se_view` list frozen controllers, shunt releases, link clusters and excluded measurements. |
| Configurable bus-shunt modeling | ✅ | `model.bus_shunt_model = "admittance"` (default) stamps bus shunts into the Y-bus; `"voltage_dependent_injection"` keeps them in the mismatch terms. Scope: the injection variant is a MATPOWER import option; it solves on nets with closed links like the admittance variant. |
| CGMES import (`importCGMES` / `createNetFromCGMES`) | ✅ | CGMES 2.4.15 bus-branch import (EQ+SSH+TP+SV, boundary sets, folders/ZIP/ZIP-in-ZIP) with transformers and their tap controllers, remote voltage controllers, retained switches as links, the short-circuit harvest, `summarizeCGMES` and `compareWithSV`. CGMES 3.0 deliveries are read as well; converters map as fixed PCC injections, `PowerElectronicsConnection` as PQ or voltage-regulating injection with `cgmes_import.map_power_electronics`. Validated on MicroGrid, SmallGrid, FullGrid and the 6209-bus RealGrid. [CGMES Import](cgmes_import.md). |
| Node-breaker topology processor | ✅ | Imports node-breaker deliveries without a TP profile (connectivity nodes aggregate across closed non-retained switches, retained switches stay couplers); runs only when no non-boundary `TopologicalNode` exists; verified against the shipped TP on the MiniGrid/SmallGrid/FullGrid sets. |
| CGMES export (`writeCGMESFiles`) | ✅ | Complete CGMES 2.4.15 delivery (EQ + TP + SSH + SV, optionally one ZIP) with roundtrip-stable identity; a re-imported net solves to the same power flow and short circuit. Converter units (`PowerElectronicsConnection`) go out as `EquivalentInjection` with their regulation and Q limits, without their fault current: a re-import flags their island in the short circuit as a lower bound. [CGMES Export](cgmes_export.md). |
| CGMES import analysis (`analyzeCGMES`) | ✅ | Explains a non-importable delivery: supplied models, declared prerequisites, unresolved-reference histogram, verdict. |
| PowSyBl IIDM import (`PowsyblAdapter`, `read_iidm_tables`) | ✅ | XIIDM networks (schema 1.0 to 1.17; compressed files are not read) read in Julia: node-breaker topology, transformers at their current tap, three-winding star buses, HVDC stations as fixed injections, SVC, dangling and tie lines, operational and reactive limits; reproduces the OpenLoadFlow solution of five example networks. [PowSyBl Import](powsybl_import.md). |
| Base-voltage inference (`cgmes_import.infer_base_voltages`) | ✅ | Reconstructs missing nominal voltages from the SV state and transformer rated voltages. |
| CGMES `EnergySource` (`cgmes_import.map_energy_source`) | ✅ | Off by default: the objects are counted in one import notice. On, each becomes a fixed PQ injection from the SSH values, without regulation and passed over by the import's slack choice; in an island fed only by EnergySource the importer makes the largest one the slack. [CGMES Import](@ref cgmes_energy_source). |
| CGMES `PowerElectronicsConnection` (`cgmes_import.map_power_electronics`) | ✅ | CGMES 3.0 converter-connected parks and batteries. Off by default: the objects are counted in one import notice. On, each becomes an injection from the SSH values, PV with an enabled voltage `RegulatingControl`, else PQ; passed over by the import's slack choice, the slack only of an island without other reference. [CGMES Import](@ref cgmes_power_electronics). |
| CGMES bank regulators (`cgmes_import.shunt_control`) | ✅ | A `LinearShuntCompensator` with `controlEnabled` and a voltage `RegulatingControl` at its own bus becomes a switched-bank voltage controller. Off by default: such banks stay fixed shunts. [Bank regulators](@ref cgmes_shunt_control). |
| Dynamic unit data (`DynamicData`) | ✅ | Inertia, damping, droop, governor lag, reserves and unit type on a generating unit (`ProSumer.dynamics`), persisted in SCF; the CGMES import reads them from the DY profile with `cgmes_import.read_dynamics` (default off). No calculation reads them; only the island slack choice of the CGMES import uses the unit type. [Dynamic Unit Data](@ref dynamic_data). |
| CGMES `EquivalentInjection` as source (`cgmes_import.equivalent_injection_as_source`) | ✅ | An `EquivalentInjection` on the SV angle reference or a boundary node becomes a voltage-regulating source at the SV voltage of its node, with the short-circuit power from its `r`/`x`; the one on the SV angle reference is the slack. Off by default: the injection stays as mapped. The SV angle reference takes the slack also when only an `EquivalentInjection` sits on it (an X-node, as in MicroGrid BE). [CGMES import configuration](cgmes_import.md#cgmes-import-config). |
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
| Local browser Web UI | ✅ | Case management, run forms with contextual help, run history, artifact viewing, the state-estimation section and QSTS runs on an SCF time series with a time series chart. Scope: loopback only, no authentication, no public deployment mode. [Web UI](webui.md). |

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

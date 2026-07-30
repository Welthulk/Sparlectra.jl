# Feature Matrix (Quick Overview)

This page gives a compact comparison of what is currently available in
Sparlectra for **Load Flow** and **State Estimation**.

Legend:

* ✅ available
* ⚠️ available with limitations / specific workflow
* ❌ not available as dedicated feature

## Network & modeling features

| Feature | Load Flow (`runpf!`) | State Estimation (`runse!`) | Notes |
|---|:---:|:---:|---|
| Framework workflow | ✅ | — | `run_sparlectra` is the preferred configuration-driven import/control/solve/output entry point and always returns one `SparlectraRunResult`; `run_acpflow` is its thin AC power-flow alias. `run_sparlectra_cases` executes configured MATPOWER batches sequentially and preserves case order. |
| Local browser Web UI | ⚠️ | ❌ | Local-only PowerFlow forms with automatic MATPOWER-case and example-configuration selection, optional standalone browser app-window launch, Markdown-backed contextual option help, an allowlisted documentation reader, result summaries, persistent run history, and safe artifact viewing/download. It binds to loopback only and intentionally has no State Estimation page or public deployment mode. |
| AC power flow (NR) | ✅ | — | Main PF entry point is `runpf!` with the sparse rectangular complex Jacobian; polar/classic PF methods are not supported. |
| AC state estimation (WLS) | — | ✅ | Main SE entry point is `runse!`. |
| Topological bus links (`addLink!`) | ✅ | ⚠️ | Links are fully integrated in PF workflow/reporting; in SE they are part of network topology context and should be used with care in measurement design. |
| Transformers | ✅ | ⚠️ | 2-winding ✅, 3-winding ✅, OLTC tap control ✅ (PF), phase shifter (PST) ✅ (PF), split combined regulation ✅ (PF), tap-changer impedance model ✅ — see the [transformer support table](@ref transformer-support) below for the full per-capability breakdown. |
| Generic outer-loop control framework | ✅ | ❌ | Reusable orchestration above `runpf!`; controller results are available through `ControlRunResult` / `latest_control_result(net)`. |
| Remote voltage control by machines (`MachineVoltageControl`) | ✅ | ❌ | A PQ machine regulates the voltage at a *different* bus via its reactive output (secant outer loop, honest `at_limit` on the reactive bounds); `addMachineVoltageControl!`, wired from CGMES remote `RegulatingControl`s with `cgmes_import.machine_control`. Theory in [Remote Voltage Control](remote_voltage_control.md). Coordinated Q-sharing of several machines on one target bus is not implemented (one controller per target). |
| Machine-readable control trace rows | ✅ | ❌ | Available through `ControlRunResult.trace`; avoids parsing console output. |
| YAML controller instantiation | ⚠️ | ❌ | `control.controllers` is reserved for future controller definitions; leave empty for current programmatic setup. |
| π-equivalent branch modeling | ✅ | ✅ | Common branch representation across PF/SE workflows. |
| Shunts / loads / generators in `Net` model | ✅ | ✅ | Shared physical network model and component handling. |
| Configurable bus-shunt modeling | ⚠️ | ❌ | `bus_shunt_model = "admittance"` is the default/classic Y-bus treatment; `"voltage_dependent_injection"` is available for rectangular PF formulations that keep shunt effects in nonlinear mismatch terms. |
| Voltage-dependent prosumer control (`Q(U)`, `P(U)`) | ✅ | ❌ | Implemented for PF with controller-aware mismatch/Jacobian terms in rectangular formulation; not part of SE model. |
| CGMES import (`importCGMES` / `createNetFromCGMES`) | ✅ | ⚠️ | ENTSO-E CGMES 2.4.15 bus-branch import (EQ+SSH+TP+SV, boundary sets, folders/ZIP/ZIP-in-ZIP): lines incl. boundary lines across nominal-voltage steps, 2W/3W transformers with fixed SvTapStep/SSH tap positions or CGMES-defined outer-loop tap controllers (`tap_control = true`), machines with a remote voltage `RegulatingControl` as outer-loop remote voltage controllers (`machine_control = true`), retained switches as bus links, loads/machines/equivalent injections with a documented slack-selection chain, always-on short-circuit data harvest, `summarizeCGMES` diagnostics and `compareWithSV` validation against the shipped SV profile (voltages and SvPowerFlow). **`Terminal.connected` uses real-snapshot semantics**: both ends open → branch out of service; exactly one end open → the branch carries no current but its half charging admittance is kept as a shunt at the closed bus. `StaticVarCompensator` maps MATPOWER-style to a P = 0 reactive injection with Q limits from the Ω ratings (`Q = vn²/X`), voltage-regulating per its control mode — the droop `slope` is deliberately ignored (notice), and a consistency guard demotes it to PQ with the SSH operating point when the target contradicts the delivered SV state. SE can use the imported `Net`; validated on MicroGrid, SmallGrid (bus-branch and node-breaker via TP), FullGrid and RealGrid (6249 buses). CGMES **3.0** deliveries are read as well (`dcat:Dataset` headers, per-border boundary files, SSH `Equipment.inService` honored for every mapped class); multi-area assemblies discard cancelling boundary equivalent pairs, split DC border crossings per side, and map `VsConverter`/`CsConverter` as Stage-0 fixed PCC injections (no angle coupling — DC-only-coupled areas stay separate islands); `AsynchronousMachine` maps as a Stage-0 fixed PQ operating point from SSH, and placeholder guards (shunt admittances > 10 × baseMVA, tap corrections outside 0.5…2.0) keep completeness-set filler values out of the solve; combined ReliCapGrid aliases (`relicapgrid_cgm`, `svedala_neighbours`) fetch several areas as one delivery. |
| Balanced short-circuit currents (`runShortCircuit!`, IEC 60909-0) | ✅ | — | Initial symmetrical short-circuit current `Ik''` (max/min), `Sk''`, and peak current `i_p` per fault bus from the CGMES short-circuit harvest: positive sequence, series impedances only, Z-bus column solve per island (one sparse LU per island/case), IEC Table-1 voltage factors with `short_circuit.c_factor` as scalar expert override. Safety-flag contract: substituted defaults and skipped contributions (motors without harvested attributes never contribute to `Ik''max` silently) are flagged on the affected result rows — a flagged maximum is a lower bound. Separate `ShortCircuitResult` + `printShortCircuitResult`; The Web UI "Short circuit" button runs both cases without a power-flow solve and writes `short_circuit_max.csv`/`short_circuit_min.csv`; unbalanced faults are not supported yet, and the `K_T`/`K_G` corrections are not applied (documented limitation). See [Short-Circuit Analysis](short_circuit.md), issue #277. |
| MATPOWER import / cases | ✅ | ✅ | Typical SE studies can start from imported PF-ready networks; PF import supports configurable SHIFT unit/sign and TAP ratio (`normal` or `reciprocal`) conventions, Sparlectra transformer-loss metadata round trips for FOR/DTF exports, plus example-workflow auto-profile recommendations for robust large-case settings. |
| Synthetic tiled-grid generator | ✅ | ⚠️ | `build_synthetic_tiled_grid_net` creates artificial one-voltage-level AC PF benchmark networks; SE can use the resulting `Net` as an artificial study case when measurements are supplied. |

### [Transformer support](@id transformer-support)

Transformer types and regulation features by winding configuration. A 3-winding
transformer is modeled as a star equivalent with an internal AUX bus
(`add3WTPiModelTrafo!` / `create3WTWindings!`); tap and phase changers sit on a
chosen winding (`tap_side`, `phase_tap_side`) and are controlled per
star-equivalent leg branch. All tap/voltage control features act in the PF
outer loop only — SE uses the same transformer network model but has no
controller support.

| Type | 2-winding | 3-winding | Remarks |
|---|:---:|:---:|---|
| Fixed-tap transformer | ✅ | ✅ | Base network model, usable in PF and SE. |
| OLTC (voltage control on ratio tap) | ✅ | ✅ | `addTapController!` with `mode = :voltage`; discrete step operation with tap limits. Remote target-bus control via `target_bus` (⚠️ one controller channel with remote measurement). |
| PST (active-power control on phase tap) | ✅ | ✅ | `mode = :branch_active_power`; discrete step operation with phase limits. |
| Combined regulation, single controller | ✅ | ✅ | `mode = :voltage_and_branch_active_power`: one controller drives ratio and phase taps together. |
| Split combined regulation (Schrägregler) | ✅ | ✅ | Two independent controllers on one unit — voltage on the ratio tap plus active power on the phase tap — each with its own target, deadband, and convergence status; per-actuator exclusivity is enforced. Demo: `tap_control_schraeg_two_controllers.jl`; theory in `control_framework.md`. |
| Symmetrical PST model (CGMES) | ⚠️ | ⚠️ | Typed `PhaseTapChangerModel` — modeling layer only, see note below. |
| Asymmetrical PST model (CGMES) | ⚠️ | ⚠️ | Includes quadrature booster as ψ = 90° — modeling layer only, see note below. |
| Tabular PST model (`:tabular`) | ✅ | ⚠️ | `TapTablePoint` lookup; tabular data overrides formula reconstruction. CGMES `PhaseTapChangerTabular` is wired into the imported branch (ratio and angle at the tap position; per-step r/x/g/b corrections reported, not applied) and validated against RealGrid's SV state. |
| Coordinated master/slave voltage control | ❌ | ❌ | Not yet implemented as dedicated multi-transformer coordination logic (no built-in participation-factor allocation/group dispatcher yet). |
| Tap-changer impedance model (`tap_changer_model`) | ✅ | ✅ | `ideal` (default) keeps the tap changer free of series-impedance feedback; `impedance_correction` re-refers transformer R/X through the tapped winding (`|1 + f·e^(jφ)|²`). Applies to all transformers of an imported case, read by both the MATPOWER and the native DTF importer; implemented centrally in `calcTapCorrectedRX`/`calcTapImpedanceCorrectionFactor` (`src/equicircuit.jl`). SE reads the same imported `Net`, so the model choice is inherited but not independently configurable per SE run. |

!!! note "Status of the typed CGMES phase-tap-changer models"
    Phase-tap-changer models (`:symmetrical`, `:asymmetrical` — quadrature booster as ψ=90° — and `:tabular` with `TapTablePoint`) can be defined on transformer windings, directly (2WT) or via `create3WTWindings!`'s `phase_tap_side`/`phase_taps` keywords (3WT, see `examples/others/exp_3wt_phase_taps.jl`). The DTF importer uses these models to derive the branch ratio and shift. Limitations: a persisted `phase_taps` model does not yet change the solved branch on its own, and there is no per-transformer configuration surface yet.

## Solvers, operations & limits

| Feature | Load Flow (`runpf!`) | State Estimation (`runse!`) | Notes |
|---|:---:|:---:|---|
| Polar full NR solver | ❌ | ⚠️ | Unsupported for PF; SE uses its own WLS iteration and Jacobian evaluation. |
| Rectangular NR solver | ✅ | ❌ | Available for PF, not as separate SE formulation. |
| Automatic rectangular Newton damping (`autodamp`) | ✅ | ❌ | PF rectangular solver can backtrack the Newton step from `damp` down to `autodamp_min` for difficult flat starts. |
| Merit-function Armijo line search (`power_flow.merit`) | ✅ | ❌ | Optional alternative step-acceptance criterion inside the autodamp backtracking loop (`f(x) = 1/2‖WF(x)‖²`, Armijo sufficient decrease); disabled by default and requires `autodamp = true`. Does not replace autodamp, the Newton solver, or candidate start-value ranking. |
| Distributed active-power slack (`power_flow.distributed_slack`) | ✅ | ❌ | The island's P imbalance (load + losses − schedule) is shared over participating generators by normalized participation factors instead of loading it entirely onto the reference bus; the reference keeps the angle reference while its P residual becomes an ordinary equation and one scalar `lambda_P` joins the state. Weight modes `pg_weighted` (default), `pmax_weighted`, `headroom_weighted`, `imported` (MATPOWER `APF` / CGMES `normalPF`), `explicit`. Disabled by default; disabled runs are bit-identical to the classical solver. Per-island `lambda_P` in island-wise runs; P limits are advisory (warn-only) in stage 1. Theory in [Solver Guide](solver.md), keys in [Power-Flow Configuration](powerflow_configuration.md). |
| Trust-region step control (`power_flow.trust_region`) | ✅ | ❌ | Optional alternative to `autodamp`: caps the Newton step norm at an adaptive radius, accepts/rejects by merit-function decrease (`rho`), and adapts the radius from the actual/predicted reduction ratio. Disabled by default and mutually exclusive with `autodamp = true`. `step_mode = :scaled` (default) rescales the full Newton direction to the radius; `step_mode = :dogleg` blends it with a steepest-descent (Cauchy) step along the dogleg path for graceful degradation when the Newton direction becomes a poor descent direction — not a rescue for infeasible cases or bad starts (Levenberg–Marquardt/Steihaug-CG remain out of scope). Reports `reason = :trust_region_collapsed` when the radius falls below `min_radius` without an accepted step (either mode). |
| Start projection (`start_projection`) | ✅ | ⚠️ | Internal PF and external-solver `PFModel` starts can use DC-angle and blend-scan projection; SE does not consume `PFModel`. |
| Guarded current-iteration start pre-solve | ⚠️ | ❌ | Optional PF-only start-value improver (`power_flow.start_current_iteration.enabled`) that runs after normal start modes and before the final rectangular NR solve; it accepts the prepared profile only when mismatch improves and guard checks pass. |
| Wrong-branch detection with full output visibility (`wrong_branch_detection`) | ✅ | ❌ | Post-convergence PF plausibility guard for suspicious low-voltage/non-finite solutions (`off|warn|fail|rescue`); a heuristic check, not a global-optimality proof or a replacement for start-value candidate selection. The result (not just the setting) is surfaced in `ACPFlowReport.metadata`, the AC island diagnostics CSV, a console summary line, the Web UI run result page, and `run_sparlectra_api` result metadata — see the wrong-branch section in [Configuration](configuration.md). The `rescue` mode is reserved (`wrong_branch_rescue_not_implemented`); the retry/rescue loop is intentionally not implemented — see below. |
| Narrative diagnostics report + fixed-reference self-check (`diagnose.log`, `run_fixed_reference_self_check`) | ✅ | ❌ | `diagnose.log` (`run_diagnostics = true`) is a diagnostic report, not a flat key/value dump: a "Diagnosis" section naming the worst-mismatch bus/equation, the mismatch-history trend (`monotonic`/`oscillatory`/`stagnant`/`diverging_to_nonfinite`), and autodamp health, a branch-anomaly scan of the branches incident to that bus (zero impedance, off-nominal tap ratio, large phase shift, unusually low X/R), and a "Recommendations" section. `run_fixed_reference_self_check` evaluates the mismatch at a case's own stored operating point (MATPOWER `VM`/`VA`, or the `SvVoltage` state of a CGMES delivery) with no corrective Newton step and every start-value machine forced off, to separate an imported-network-model issue from a solver start/step-control issue; it writes `self_check.log` (forced settings, start-state residual, CGMES SV coverage) and `self_check_residuals.csv` (full per-bus P/Q residuals with transformer/shunt adjacency for attribution). The Web UI PowerFlow form has a dedicated "Diagnose" action (distinct from "Start PowerFlow run") that runs this self-check and writes the same enriched `diagnose.log`. |
| Sparse PF matrices | ✅ | ⚠️ | PF requires sparse Y-bus and Jacobian matrices; SE internally builds Jacobians for WLS. |
| Flat start control | ✅ | ✅ | Available in both PF and SE workflows. |
| PV/PQ reactive limit handling | ✅ | ❌ | PF includes the default active-set Q-limit logic plus classical simultaneous and one-at-a-time outer-loop modes; SE currently does not expose PV/PQ switching logic. |
| `Q(U)` / `P(U)` controller solver support | ⚠️ | ❌ | Supported on the default rectangular PF path; legacy polar/classic PF modes are unsupported. |
| External solver interface | ✅ | ❌ | PF has external solver integration; SE is internal WLS. |
| APSLF solver (`power_flow.solver = apslf`) | ⚠️ | ❌ | Analytic power-series solver bridged from the optional weak dependency AnalyticLoadFlow.jl (`ApslfSolver`, `apslf_solver()`); usable standalone (`runpf_external!`), as the framework solver (`run_sparlectra`, including per-island handling), or as a guarded start-value generator ahead of the rectangular NR solve (`power_flow.apslf_start`). No selectable start voltage (always the canonical analytic germ), no OLTC/PST/Q(U)/P(U) controller support (rejected up front), and only simple internal PV→PQ Q-limit switching (no active-set guard/classical outer-loop modes). See [External Solver Interface](external_solvers.md#apslf-analyticloadflowjl). |
| DC power flow (`rundcpf!`, `power_flow.solver = dc`) | ✅ | ❌ | Standalone linear screening model (MATPOWER `rundcpf`/`makeBdc` equivalent): series-reactance-only `B'`, phase-shift injection vector, `Vm` implicitly 1.0 pu, lossless. Usable standalone (`rundcpf!`) or as the framework solver (`run_sparlectra`, including per-island handling). `rundcpf!(net; seed_ac_start=true)` optionally chains an AC Newton-Raphson solve seeded from the DC angles. A phase-shifting transformer's current, fixed angle is represented in the B′/injection math; no OLTC/PST/Q(U)/P(U) outer-loop controller support (rejected up front, mirrors `apslf`). |

## State-estimation measurements, observability & diagnostics

| Feature | State Estimation (`runse!`) | Notes |
|---|:---:|---|
| SCADA-style measurements (`Vm`, `Pinj`, `Qinj`, `Pflow`, `Qflow`) | ✅ | Public measurement types and helper builders are available. |
| PMU voltage-phasor measurements (`Vm` + `Va`) | ✅ | Magnitude as tightly weighted `VmMeas`, angle as `VaMeas` (degrees) with an estimated common reference-angle offset α between the PMU time base and the slack reference (`state_estimation.pmu_ref_offset`). Combined helper: `addPmuPhasorMeasurement!`. |
| PMU current-phasor measurement types | ❌ | No dedicated branch current phasor (`I_ij`) types yet. |
| Passive bus / zero-injection (ZIB) support | ⚠️ | Implemented via zero-injection pseudo-measurements (not separate hard-constraint block). |
| Global/local observability analysis | ✅ | Matrix-level and network-level observability helpers are available. |
| Structural observability checks | ✅ | Sparsity/matching-based checks are available. |
| Numerical observability checks | ✅ | Rank/SVD-based checks are available. |
| Local observability on selected state subset | ✅ | Dedicated local observability helpers are available. |
| Bad-data diagnostics (global consistency, residual ranking) | ✅ | `validate_measurements`, `runse_diagnostics`, `summarize_se_diagnostics`, `print_se_diagnostics`. |
| Deactivate-and-rerun helper | ✅ | Optional one-step rerun on top suspicious measurement; can improve objective but may still remain globally inconsistent. |
| Markdown/plain diagnostics output | ✅ | `print_se_diagnostics(...; format=:markdown|:plain)`. |

## Reporting, export & workflow helpers

| Feature | Load Flow (`runpf!`) | State Estimation (`runse!`) | Notes |
|---|:---:|:---:|---|
| Human-readable result printing | ✅ | ✅ | PF and SE both provide textual result output/reporting helpers. |
| Machine-readable report (`ACPFlowReport`) | ✅ | ❌ | Dedicated report container currently exists for PF workflow. |
| DataFrame-friendly report rows | ✅ | ❌ | PF report rows can be converted/used in tabular workflows. |
| Synthetic measurements from PF result | — | ✅ | PF + measurement generators support SE test-data workflows. |
| Central typed configuration | ✅ | ✅ | `SparlectraConfig` and module-specific config sections support cached YAML loading, typed validation, override precedence, and effective-configuration printing for application/example boundaries. |
| GUI-ready programmatic run API | ✅ | ❌ | `run_sparlectra_api` provides unique stable run IDs, schema-versioned structured status, controlled configuration overrides, effective configuration output, serialization, and explicit artifact discovery for MATPOWER power-flow runs. |
| Local PowerFlow service boundary | ✅ | ❌ | `start_powerflow_run`, persistent run indexing, restart recovery, result lookup, artifact listing, and safe artifact resolution provide a filesystem-backed boundary for a future local GUI without HTTP or Genie.jl dependencies. |
| Write-back solved states into `Net` | ✅ | ✅ | PF updates net states; SE supports `updateNet=true`. |

## Useful links

* [State Estimation](state_estimation.md)
* [Links (bus couplers)](links.md)
* [External Solvers](external_solvers.md)
* [Network Reports](netreports.md)
* [Branch Model](branchmodel.md)
* [Workshop](workshop.md)
* [Changelog](changelog.md)
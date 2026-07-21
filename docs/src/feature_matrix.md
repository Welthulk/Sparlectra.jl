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
| AC state estimation (WLS) | — | ✅ | Main SE entry point is `runse!` (experimental status). |
| Topological bus links (`addLink!`) | ✅ | ⚠️ | Links are fully integrated in PF workflow/reporting; in SE they are part of network topology context and should be used with care in measurement design. |
| 2-winding transformer | ✅ | ✅ | Supported in network model and usable in both workflows. |
| 3-winding transformer | ✅ | ✅ | Implemented via star-equivalent with AUX bus in network construction. |
| Generic outer-loop control framework | ✅ | ❌ | Reusable orchestration above `runpf!`; controller results are available through `ControlRunResult` / `latest_control_result(net)`. |
| Transformer tap/phase controller as outer controller | ✅ | ❌ | First concrete `AbstractOuterController`; supports ratio and phase updates outside the Newton system. |
| Machine-readable control trace rows | ✅ | ❌ | Available through `ControlRunResult.trace`; avoids parsing console output. |
| YAML controller instantiation | ⚠️ | ❌ | `control.controllers` is reserved for future controller definitions; leave empty for current programmatic setup. |
| Transformer tap control (`addTapController!`) | ✅ | ❌ | PF supports outer-loop tap control for ratio and/or phase (`:voltage`, `:branch_active_power`, `:voltage_and_branch_active_power`), including discrete step operation with tap/phase limits. |
| Remote target-bus voltage control (single-controller) | ⚠️ | ❌ | Supported in PF by setting `mode = :voltage` and `target_bus`; this is remote measurement with one controller channel. |
| Coordinated master/slave transformer voltage control | ❌ | ❌ | Not yet implemented as dedicated multi-transformer coordination logic (no built-in participation-factor allocation/group dispatcher yet). |
| π-equivalent branch modeling | ✅ | ✅ | Common branch representation across PF/SE workflows. |
| Shunts / loads / generators in `Net` model | ✅ | ✅ | Shared physical network model and component handling. |
| Configurable bus-shunt modeling | ⚠️ | ❌ | `bus_shunt_model = "admittance"` is the default/classic Y-bus treatment; `"voltage_dependent_injection"` is available for rectangular PF formulations that keep shunt effects in nonlinear mismatch terms. |
| Voltage-dependent prosumer control (`Q(U)`, `P(U)`) | ✅ | ❌ | Implemented for PF with controller-aware mismatch/Jacobian terms in rectangular formulation; not part of SE model. |
| MATPOWER import / cases | ✅ | ✅ | Typical SE studies can start from imported PF-ready networks; PF import supports configurable SHIFT unit/sign and TAP ratio (`normal` or `reciprocal`) conventions, Sparlectra transformer-loss metadata round trips for FOR/DTF exports, plus example-workflow auto-profile recommendations for robust large-case settings. |
| Tap-changer model (`transformer.tap_changer_model`) | ✅ | ⚠️ | `ideal` (default) keeps the tap changer free of series-impedance feedback; `impedance_correction` re-refers transformer R/X through the tapped winding (`|1 + f·e^(jφ)|²`). Applies to all transformers of an imported case, read by both the MATPOWER and the native DTF importer; implemented centrally in `calcTapCorrectedRX`/`calcTapImpedanceCorrectionFactor` (`src/equicircuit.jl`). SE reads the same imported `Net`, so the model choice is inherited but not independently configurable per SE run. |
| Typed phase-tap-changer models (CGMES PST) | ⚠️ | ❌ | `PhaseTapChangerModel` (`:symmetrical`/`:asymmetrical`, quadrature booster as ψ=90°) and `:tabular` models with `TapTablePoint` provide CGMES-oriented phase-shifter semantics; formula/lookup helpers (`calcPhaseTapAngleRatio`, `calcPhaseTapReactance`, `calcPhaseTapTable`) are centralized in `src/equicircuit.jl` (Issue #261). Currently a developer-facing modeling/equivalent-circuit layer: `PowerTransformerWinding.phase_taps` can be populated directly (2WT) or via `create3WTWindings!`'s `phase_tap_side`/`phase_taps` keywords (3WT, `examples/exp_3wt_phase_taps.jl`), and the DTF importer uses `calcPhaseTapAngleRatio` to derive branch `ratio`/`shift` (result unchanged) — but a persisted `phase_taps` model has no effect on the solved branch yet: `calcPhaseTapReactance` is not wired into the solver, `Branch.phase_min_deg/phase_max_deg/phase_step_deg` are still hard-coded constants regardless of `phase_taps`, and there is no config-driven per-transformer selection. Tabular data overrides formula reconstruction where present. |
| Synthetic tiled-grid generator | ✅ | ⚠️ | `build_synthetic_tiled_grid_net` creates artificial one-voltage-level AC PF benchmark networks; SE can use the resulting `Net` as an artificial study case when measurements are supplied. |

## Solvers, operations & limits

| Feature | Load Flow (`runpf!`) | State Estimation (`runse!`) | Notes |
|---|:---:|:---:|---|
| Polar full NR solver | ❌ | ⚠️ | Unsupported for PF; SE uses its own WLS iteration and Jacobian evaluation. |
| Rectangular NR solver | ✅ | ❌ | Available for PF, not as separate SE formulation. |
| Automatic rectangular Newton damping (`autodamp`) | ✅ | ❌ | PF rectangular solver can backtrack the Newton step from `damp` down to `autodamp_min` for difficult flat starts. |
| Merit-function Armijo line search (`power_flow.merit`) | ✅ | ❌ | Optional alternative step-acceptance criterion inside the autodamp backtracking loop (`f(x) = 1/2‖WF(x)‖²`, Armijo sufficient decrease); disabled by default and requires `autodamp = true`. Does not replace autodamp, the Newton solver, or candidate start-value ranking. |
| Trust-region step control (`power_flow.trust_region`) | ✅ | ❌ | Optional scaled-Newton alternative to `autodamp`: caps the Newton step norm at an adaptive radius, accepts/rejects by merit-function decrease (`rho`), and adapts the radius from the actual/predicted reduction ratio. Disabled by default and mutually exclusive with `autodamp = true`. Scaled-Newton only (no dogleg/Steihaug); reports `reason = :trust_region_collapsed` when the radius falls below `min_radius` without an accepted step. |
| Start projection (`start_projection`) | ✅ | ⚠️ | Internal PF and external-solver `PFModel` starts can use DC-angle and blend-scan projection; SE does not consume `PFModel`. |
| Guarded current-iteration start pre-solve | ⚠️ | ❌ | Optional PF-only start-value improver (`power_flow.start_current_iteration.enabled`) that runs after normal start modes and before the final rectangular NR solve; it accepts the prepared profile only when mismatch improves and guard checks pass. |
| Wrong-branch detection with full output visibility (`wrong_branch_detection`) | ✅ | ❌ | Post-convergence PF plausibility guard for suspicious low-voltage/non-finite solutions (`off|warn|fail|rescue`); a heuristic check, not a global-optimality proof or a replacement for start-value candidate selection. The result (not just the setting) is surfaced in `ACPFlowReport.metadata`, the AC island diagnostics CSV, a console summary line, the Web UI run result page, and `run_sparlectra_api` result metadata — see `docs/src/configuration.md#wrong-branch-detection-semantics-rectangular-pf`. The `rescue` mode is reserved (`wrong_branch_rescue_not_implemented`); the retry/rescue loop is intentionally not implemented — see below. |
| Sparse PF matrices | ✅ | ⚠️ | PF requires sparse Y-bus and Jacobian matrices; SE internally builds Jacobians for WLS. |
| Flat start control | ✅ | ✅ | Available in both PF and SE workflows. |
| PV/PQ reactive limit handling | ✅ | ❌ | PF includes the default active-set Q-limit logic plus classical simultaneous and one-at-a-time outer-loop modes; SE currently does not expose PV/PQ switching logic. |
| `Q(U)` / `P(U)` controller solver support | ⚠️ | ❌ | Supported on the default rectangular PF path; legacy polar/classic PF modes are unsupported. |
| External solver interface | ✅ | ❌ | PF has external solver integration; SE is internal WLS. |
| APSLF solver (`power_flow.solver = apslf`) | ⚠️ | ❌ | Analytic power-series solver bridged from the optional weak dependency AnalyticLoadFlow.jl (`ApslfSolver`, `apslf_solver()`); usable standalone (`runpf_external!`), as the framework solver (`run_sparlectra`, including per-island handling), or as a guarded start-value generator ahead of the rectangular NR solve (`power_flow.apslf_start`). No selectable start voltage (always the canonical analytic germ), no OLTC/PST/Q(U)/P(U) controller support (rejected up front), and only simple internal PV→PQ Q-limit switching (no active-set guard/classical outer-loop modes). See [External Solver Interface](external_solvers.md#apslf-analyticloadflowjl). |

## State-estimation measurements, observability & diagnostics

| Feature | State Estimation (`runse!`) | Notes |
|---|:---:|---|
| SCADA-style measurements (`Vm`, `Pinj`, `Qinj`, `Pflow`, `Qflow`) | ✅ | Public measurement types and helper builders are available. |
| PMU-specific native measurement model (e.g. direct phasor angle/current types) | ❌ | No dedicated PMU enum/types yet; current API is SCADA-style core WLS measurements. |
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
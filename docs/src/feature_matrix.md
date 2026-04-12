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
| AC power flow (NR) | ✅ | — | Main PF entry point is `runpf!` with rectangular complex Jacobian by default; `:polar_full` / `:classic` are deprecated. |
| AC state estimation (WLS) | — | ✅ | Main SE entry point is `runse!` (experimental status). |
| Topological bus links (`addLink!`) | ✅ | ⚠️ | Links are fully integrated in PF workflow/reporting; in SE they are part of network topology context and should be used with care in measurement design. |
| 2-winding transformer | ✅ | ✅ | Supported in network model and usable in both workflows. |
| 3-winding transformer | ✅ | ✅ | Implemented via star-equivalent with AUX bus in network construction. |
| π-equivalent branch modeling | ✅ | ✅ | Common branch representation across PF/SE workflows. |
| Shunts / loads / generators in `Net` model | ✅ | ✅ | Shared physical network model and component handling. |
| MATPOWER import / cases | ✅ | ✅ | Typical SE studies can start from imported PF-ready networks. |

## Solvers, operations & limits

| Feature | Load Flow (`runpf!`) | State Estimation (`runse!`) | Notes |
|---|:---:|:---:|---|
| Polar full NR solver (deprecated) | ✅ | ⚠️ | Legacy PF mode (deprecated); SE uses its own WLS iteration and Jacobian evaluation. |
| Rectangular NR solver | ✅ | ❌ | Available for PF, not as separate SE formulation. |
| Finite-difference Jacobian option (`opt_fd`) | ✅ | ❌ | In PF this toggles FD Jacobian vs analytic Jacobian (not fast-decoupled); not exposed as an SE user option. |
| Sparse Jacobian option (`opt_sparse`) | ✅ | ⚠️ | PF supports explicit sparse option; SE internally builds Jacobians for WLS. |
| Flat start control | ✅ | ✅ | Available in both PF and SE workflows. |
| PV/PQ reactive limit handling | ✅ | ❌ | PF includes Q-limit logic; SE currently does not expose PV/PQ switching logic. |
| External solver interface | ✅ | ❌ | PF has external solver integration; SE is internal WLS. |

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
| Write-back solved states into `Net` | ✅ | ✅ | PF updates net states; SE supports `updateNet=true`. |

## Useful links

* [State Estimation](state_estimation.md)
* [Links (bus couplers)](links.md)
* [External Solvers](external_solvers.md)
* [Network Reports](netreports.md)
* [Workshop](workshop.md)
* [Changelog](changelog.md)

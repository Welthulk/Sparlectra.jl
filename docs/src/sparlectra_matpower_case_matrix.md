# Sparlectra.jl MATPOWER Case Diagnostics Matrix

Date: 2026-07-14
Sparlectra versions covered: 0.7.7–0.8.9
Purpose: compact tracking table for MATPOWER import, reference-data consistency, and solver diagnostics.

> **Data and license note**  
> The MATPOWER case files referenced here are publicly available from the MATPOWER project or related public benchmark-data sources. This document does not redistribute those case files. It only records diagnostic observations from local Sparlectra.jl runs. MATPOWER itself is distributed under the 3-clause BSD license, but the case files included with MATPOWER may have their own origin, citation requests, or permission conditions. Before redistributing case data or derived data sets, check the header comments of the corresponding case file and the applicable MATPOWER documentation.
>
> **Non-binding diagnostic note**  
> The observations in this matrix are engineering diagnostics from specific Sparlectra.jl runs and configuration states. They are not legal advice, not a formal validation certificate, and not a legally binding assessment of the correctness, licensing status, or benchmark suitability of the referenced data sets.

This note is intended as a living diagnostic matrix. Add one row per MATPOWER case and keep the findings short.

The matrix distinguishes solver convergence failures, wrong Newton solution branches, Sparlectra import/Y-bus mismatches, and inconsistent stored MATPOWER reference data. When a large fixed-reference residual is already reproducible with a raw MATPOWER-style Y-bus using the stored `VM`/`VA` values and branch data, classify the case as a reference-data consistency diagnostic rather than as a clean solver benchmark failure.

## Latest addition: `case_SyntheticUSA.m` in v0.8.9

`case_SyntheticUSA.m` is documented from a successful v0.8.9 Web UI/API-path run. The network has 82,000 buses split across three disconnected AC islands: island 1 has 70,000 buses and 88,209 active AC branches, island 2 has 10,000 buses and 12,706 active AC branches, and island 3 has 2,000 buses and 3,206 active AC branches. The case contains nine active `mpc.dcline` rows. Each active DC-line row is imported as two fixed terminal injections; this steady-state PF injection approximation does not create an AC branch, dummy admittance, impedance bridge, or complete HVDC control/OPF model, and it does not electrically merge the AC Y-bus components.

Validated robust profile:

```yaml
power_flow:
  tol: 1.0e-5
  max_iter: 80
  autodamp: true
  autodamp_min: 0.01
  islands:
    enabled: true
    mode: solve_independent
    diagnostic_continue_after_failure: true
  start_mode:
    angle_mode: dc
    voltage_mode: profile_blend
    profile_source: matpower_reference
  qlimits:
    enabled: false
  start_current_iteration:
    enabled: false

matpower_import:
  matpower_dcline_mode: pf_injections
  auto_profile: apply
```

Verified result: Web UI/API status succeeded, all three islands converged, aggregate iteration count was 68, and final mismatch was approximately `5.597176064853215e-6`. Island 1 used reference bus 30902 and converged in 29 iterations with final mismatch approximately `5.597176064853215e-6`; island 2 promoted bus 70684 and converged in 32 iterations with final mismatch approximately `9.838130310413362e-12`; island 3 promoted bus 80004 and converged in 7 iterations with final mismatch approximately `9.086420504900161e-11`.

A diagnostic run with `autodamp=false`, Q-limits disabled, and current-iteration start disabled produced `nr_nonfinite` on islands 1 and 2 while island 3 converged. Use the validated autodamped profile above for reproducibility.

## Latest addition: `case_ACTIVSg25k.m`

`case_ACTIVSg25k.m` now converges successfully with the strengthened Q-limit guard profile. The final documented run reports `numerical_solution=OK`, `q_limit_active_set=OK`, `final_converged=true`, `iterations=8`, `pv2pq_events=2004`, `guarded narrow-Q PV buses=960`, and `compare=OK`. The remaining active PV/REF Q-limit violations are zero. This confirms that the earlier failure mode was primarily a Q-limit active-set stabilization issue, not a classical Newton divergence.

## Latest addition: `case_ACTIVSg10k.m` wrong-branch follow-up in v0.8.2

The historical ACTIVSg10k wrong-branch observation remains important, but it should now be classified as a **historical solver-path issue**, not as a currently reproducible natural wrong-branch case with the latest rectangular solver path.

In the current v0.8.2 checks, several default-near and degraded start configurations were tested without artificial branch-angle threshold tightening:

- `profile_blend` / MATPOWER reference voltage starts converge to a healthy solution.
- `classic` / `classic` voltage-angle starts no longer produce a formally converged wrong branch; they fail before numerical convergence with active-set instability.
- Delayed or disabled Q-limit variants either converge to the same healthy solution or fail Q-limit validation, but do not show wrong-branch quality warnings.
- No tested configuration with normal `wrong_branch_max_branch_angle_deg = 90.0`, normal `shift_sign = 1.0`, `shift_unit = deg`, and `ratio = normal` reproduced the old low-voltage wrong-branch solution.

The current observed healthy ACTIVSg10k branch-quality metrics are approximately:

```text
min_vm_pu                      ≈ 0.9234 … 0.9460
max_vm_pu                      ≈ 1.0889 … 1.1087
wrong_branch_angle_spread_deg  ≈ 107.1 … 107.7
wrong_branch_worst_branch_deg  ≈ 21.75 … 22.94
wrong_branch_status            = ok
```

The v0.8.2 wrong-branch configuration forwarding was separately validated with a **forced threshold test** on the same solved case:

- With `wrong_branch_max_branch_angle_deg = 90.0`, the status is `ok`.
- With `wrong_branch_max_branch_angle_deg = 20.0` and `wrong_branch_detection = warn`, the status becomes `warn`, reason `branch_angle_exceeded`, with 3 branch-angle violations.
- With the same `20.0` threshold and `wrong_branch_detection = fail`, the numerical solve still converges, but the final PF status becomes `wrong_branch_detected` and `final_converged = false`.

This forced-threshold test proves that the `wrong_branch_*` YAML values are now forwarded into the rectangular solver and that `warn`/`fail` final-status handling works. It is **not** evidence that the physical ACTIVSg10k case naturally lands on a wrong branch in v0.8.2.

## Case matrix

| Case | Size / type | Status in Sparlectra | Recommended YAML profile | Main findings | Case-specific anomalies / deviations | Wrong-branch risk | Open checks |
|---|---:|---|---|---|---|---|---|
| `case_SyntheticUSA.m` | 82,000 buses; three disconnected AC islands; nine active `mpc.dcline` rows | v0.8.9 Web UI/API path succeeds with independent island solving and PF-injection DC-line import: all islands converged, aggregate iterations 68, final mismatch `≈5.597176064853215e-6`. | Validated robust SyntheticUSA profile documented above | Active DC lines are fixed terminal P/Q injections only; no AC branch or Y-bus bridge is created. Reference buses: 30902, promoted 70684, promoted 80004. | Low with autodamping. Non-autodamped diagnostic produced `nr_nonfinite` on islands 1 and 2. | Do not redistribute the case file; keep MATPOWER/case licensing checks external. |
| `case_ACTIVSg10k.m` | 10,000 buses; large synthetic transmission case | Current v0.8.2 file-based MATPOWER runs converge to a healthy rectangular solution with profile-blend/MATPOWER-reference starts. Typical status: `status=converged`, `numerical_converged=true`, `final_converged=true`, `final_mismatch≈5e-8` or smaller, and `wrong_branch_status=ok`. | `large_flatstart_dc_blend` for robust production runs; `activsg10k_wrong_branch_historical_scan` only for manual diagnostics | MATPOWER import conventions confirmed: `shift_sign=1.0`, `shift_unit=deg`, `ratio=normal`, `bus_shunt_model=admittance`. Branch-shift scan strongly rejects opposite sign, radians, and reciprocal ratio. Bus shunts are required for reactive balance. Current healthy runs show approximately `min_vm_pu≈0.9234…0.9460`, `max_vm_pu≈1.0889…1.1087`, `angle_spread≈107°`, and worst active branch angle `≈22°`. | 273 PV/REF buses have no online generator; fallback to `BUS.VM` is required. Many PV buses have `Qmin=0` or `Qmax=0`. Large number of PV→PQ switches can occur when Q-limits are enforced. A delayed-Q-limit scan can produce `converged_limits_failed` due to remaining PV Q-limit violations while branch-quality metrics remain `ok`. | Historical, not currently reproduced. An earlier classic flat-start run converged formally to a low-voltage / wrong-branch solution with `max|dVm|≈1.0187 pu`, `max|dVa|≈90.73°`, `minVm=0.0 pu`, but current v0.8.2 scans either converge to a healthy branch or fail before numerical convergence. `wrong_branch_detection` itself is validated by forced-threshold tests (`90° => ok`, `20° warn => branch_angle_exceeded`, `20° fail => wrong_branch_detected`). | Preserve the historical wrong-branch note, but do not treat it as a current natural reproducer. If a true reproducer is still required, search old commits/configs rather than adding more YAML guesses. Keep compact Q-limit reporting, PV/REF-without-online-generator diagnostics, and forced-threshold wrong-branch regression coverage. |
| `case_ACTIVSg25k.m` | 25,000 buses; large ACTIVSg synthetic transmission case; 32,230 branches | Rectangular NR converges with strengthened Q-limit guard. Successful run: `numerical_solution=OK`, `q_limit_active_set=OK`, `final_converged=true`, `status=converged`, `iterations=8`, PF time `≈74.55 s`, total example runtime `≈85.11 s`, final mismatch `≈4.651e-08`. | `activsg25k_q_limit_guard` + `large_flatstart_dc_blend` | MATPOWER fixed-reference self-check is good for this size: `max|ΔP|≈7.118 MW`, `max|ΔQ|≈4.611 MVAr`. Auto-profile detects `2034/3779` online generators with zero or narrow Q range. Q-limit guard locks `960` narrow-range PV buses as PQ before rectangular NR. Active set converges cleanly: `PV violations=0`, `REF violations=0`, final `PV/REF Q-limit check=OK`. Switching statistics: `pv2pq_events=2004`, `pv2pq_buses=1996`, `oscillating_buses=0`. Compare against imported setpoints is `OK`: `max|dVm|≈0.04291 pu`, aligned `max|dVa|≈0.3683°`, slack Δ `≈0.0°`. | `2753/3235` PV/REF buses have online `GEN.VG` targets; `502` differ from `BUS.VM` by more than `1e-4 pu`. `482` PV/REF buses have no online generator and use `BUS.VM` fallback. Negative branch impedance is present and preserved: `BR_R<0` on `447` rows, `BR_X<0` on `503` rows, both negative on `447` rows. Compare reference kinds: `active_pv_imported_setpoint=1238`, `final_pq_after_qlimit=1996`, `pq_bus_vm=21765`, `ref_slack_imported_setpoint=1`. | Low with the successful profile. Previous failures were active-set stabilization failures caused by many zero/narrow-Q generators and infeasible PV setpoints, not classical Newton wrong-branch convergence. | Console output is still too verbose: keep full auto-profile, negative-branch diagnostics, PV voltage diagnostics, and per-bus PV→PQ events in the logfile; show only compact statistics on console. Consider whether `qlimit_guard_violation_mode=lock_pq` should become an auto-profile recommendation for large narrow-Q cases. |
| `case300.m` | 300 buses; medium test case | Rectangular NR converges in about 5 iterations with 3 PV→PQ switches. Final voltage magnitudes are close, but angle comparison fails mainly around the `BUS_I 196 / 2040` area. | `large_flatstart_dc_blend` + `diagnostic_import_scan` | The MATPOWER fixed-reference self-check is not power-balanced around `BUS_I 196 / 2040`. The dominant mismatch is already reproducible with a MATPOWER-style Y-bus using the stored `VM`/`VA` values and raw branch data. Main residual: about 926.9 MW at `BUS_I 2040` and `BUS_I 196`. | Not classified as a Sparlectra Newton wrong-branch issue. The dominant mismatch appears to be caused by the stored case reference angles together with a low-reactance active branch `196 -> 2040` (`x = 0.02`, `TAP = 1`, `SHIFT = 0`). Treat this case as a reference-data consistency diagnostic, not as a clean pass/fail solver benchmark. | Low for the documented run; angle-reference comparison is dominated by the stored fixed-reference inconsistency rather than a known wrong branch. | Keep the fixed-reference diagnostic classification visible when updating MATPOWER comparison tolerances or profiles. Keep the explicit `TAP = 1` preservation regression, but do not treat lost nominal TAP as the remaining root cause. |
| `case145.m` | 145 buses; medium/small test case | Rectangular NR converges: `converged=yes`, `iterations=6`, `pv2pq_events=1`, compare `status=OK`. | `large_flatstart_dc_blend` + `diagnostic_import_scan` | No branch-shift entries. Standard degree/radian/sign scans are irrelevant because `SHIFT=0`. Reciprocal ratio is rejected by diagnostics. Bus shunts are important: disabling them causes very large P/Q residuals. | One relevant PV→PQ switch: `BUS_I 104`, where `GEN.VG=1.045` differs from `BUS.VM≈1.0059`; after switching, comparison against imported setpoint shows `dVm≈-0.03909 pu`, but the hybrid compare is still within configured tolerances. Several PV buses have `Qmin=0`. | Low in the documented run. No evidence of a wrong Newton branch; angle deviation is very small (`max|dVa|≈0.0051°`). | Track whether `BUS_I 104` is a valid Q-limit-driven PV→PQ case. Keep zero-Q-limit and shunt sensitivity visible in future comparisons. |
| `case1354pegase.m` | 1,354 buses; PEGASE case | With standard MATPOWER-style shift interpretation (`sign=1.0`, `unit=deg`) the compare fails. With PEGASE-style interpretation `matpower_shift_sign=-1.0`, `matpower_shift_unit=rad`, rectangular NR converges in 3 iterations with `pv2pq_events=0`; compare `status=OK`, `max|dVm|≈0.00254 pu`, `max|dVa|≈0.47247°`. | `pegase_shift_rad_opposite` + `large_flatstart_dc_blend` + `diagnostic_import_scan` | Branch-shift scan identifies `opposite sign, radians` as the only plausible convention. It reduces the fixed-reference active-power residual from about 1299.8 MW to about 38.8 MW. Bus shunts remain important for Q balance. | Six non-zero shift entries. Remaining fixed-reference mismatch is dominated by `BUS_I 58` and `BUS_I 6153`. No PV/REF buses without online generators. Active PV/REF setpoint residual is essentially zero after solve. | Low in the documented successful run. The earlier failed compare was mainly an import-convention issue, not a demonstrated wrong Newton branch. | Add an automatic PEGASE shift-convention recommendation when branch-shift scan strongly prefers `sign=-1.0`, `unit=rad`. Keep residual classification visible because the fixed-reference check is not perfectly zero even with the preferred convention. |
| `case2869pegase.m` | 2,869 buses; PEGASE case | With PEGASE-style shift interpretation `matpower_shift_sign=-1.0`, `matpower_shift_unit=rad`, rectangular NR converges robustly: `converged=yes`, `iterations=3`, `pv2pq_events=0`, compare `status=OK`, `max|dVm|≈0.00129 pu`, `max|dVa|≈0.25599°`. | `pegase_shift_rad_opposite` + `large_flatstart_dc_blend` + `diagnostic_import_scan` | Branch-shift scan strongly rejects standard degree-based interpretations. The preferred `opposite sign, radians` convention gives the best fixed-reference residuals: `max|dP|≈20.55 MW`, `max|dQ|≈50.54 MVAr`; standard MATPOWER sign/degrees gives `max|dP|≈4202.63 MW`. Bus shunts are important for Q balance: disabling them raises `max|dQ|` to about 391.2 MVAr. | Twelve non-zero branch SHIFT entries. No PV→PQ switches occurred in the documented run despite many PV buses and Q-limit entries. No PV/REF buses without online generators. Active PV/REF setpoint residual is essentially zero (`≈9.6e-12 pu` over 510 buses). Top remaining compare deviations are small and concentrated around `BUS_I 1541`, `221`, `58`, and `6153`. | Low in the documented run. No evidence of a wrong Newton branch when the PEGASE shift convention is used. | Add automatic detection/reporting for PEGASE-like radian shift values with opposite sign. Keep this case as a positive regression for `sign=-1.0`, `unit=rad`, `ratio=normal`, `bus_shunt_model=admittance`. |
| `case9241pegase.m` | 9,241 buses; large PEGASE case | Latest run `run_case9241pegase.m_20260515_132616.log`: with `opt_flatstart=true`, DC-angle start, `profile_blend` with `profile_source: matpower_reference`, start projection enabled, PEGASE-style `matpower_shift_sign=-1.0`, `matpower_shift_unit=rad`, `matpower_ratio=normal`, and `ignore_q_limits=true`, rectangular NR converges: `converged=yes`, `iterations=4`, `pv2pq_events=0`, `pv2pq_buses=0`, PF time `≈10.43 s`. Strict compare still fails by angle tolerance: `max|dVm|≈0.02662 pu`, `max|dVa|≈1.48863°` with `tol_va=0.99°`; the later `status=OK` line is only the deliberately relaxed `1e9` tolerance check. | `pegase_shift_rad_opposite` + `large_flatstart_dc_blend` + `diagnostic_import_scan`; for this diagnostic run also `ignore_q_limits=true` | Solver convergence and PV setpoint tracking are good: post-solve active PV/REF residual `max|Vm_calc-imported_Vset|≈1.21e-13 pu` over 1445 buses and no PV/REF buses without online generators. The fixed-reference self-check remains inconsistent against the imported model: `max|ΔP|≈2.71965 pu` (`≈271.97 MW`) and `max|ΔQ|≈3.92020 pu` (`≈392.02 MVAr`). The branch-shift scan still selects `opposite sign + radians + normal ratio` as best global convention; disabling bus shunts does not improve P and worsens Q. | Remaining deviations are concentrated in local residual clusters and special branch regions. Residual clustering reports 31 clusters above 10 MW/MVAr; the dominant clusters include `5491/5673`, `7934/6621`, and `2309/8414`. Notable branches include rows `14747/14748` (`7934 -> 6621`, `TAP=0.976337`, `SHIFT=-0.045364`), row `14763` (`2309 -> 8414`, `TAP=0.946281`, `SHIFT=-0.010988`), row `14766` (`503 -> 8414`, `TAP=0.993682`, `SHIFT=+0.080858`), and negative-R high-residual branches `14764` (`2309 -> 8414`, `BR_R=-0.011591`) and `14767` (`503 -> 8414`, `BR_R=-0.010888`). The log reports 75 branches with `BR_R<0` and 16 with `BR_X<0`; these are preserved, not clipped. | Low in the documented run. No evidence of a low-voltage/wrong-branch solution; this is an import/reference-consistency diagnostic, not a clean Newton failure. | Inspect branch-local TAP/SHIFT conventions and residual attribution for the top clusters, especially duplicated transformer rows and adjacent negative-R branches. Add branch-group aggregation for ordinary lines, TAP-only, SHIFT-only, TAP+SHIFT, and negative R/X. Re-test with Q-limits enabled only after the fixed-reference mismatch is understood. |
| `case13659pegase.m` | 13,659 buses; large PEGASE case | Current Web UI/API work documents the import-convention scan before interpreting NR or Q-limit active-set failures. Full NR controls were not completed locally because the full auto-profile including shunt scoring was too expensive for this agent run, but the branch convention scan completed. | `pegase_shift_rad_opposite` + `large_flatstart_dc_blend` + `diagnostic_import_scan`; begin with Q-limits disabled | Branch convention residuals strongly prefer `shift_unit=rad`, `shift_sign=-1.0`, `ratio=normal`: `max|ΔP|≈70.228 MW`, `max|ΔQ|≈147.380 MVAr`, score `≈1.473797 pu`. With the Web UI/API baseline (`shift_unit=deg`, `shift_sign=1.0`, `ratio=normal`) the score is much worse: `max|ΔP|≈6299.572 MW`, `max|ΔQ|≈878.534 MVAr`, score `≈62.995724 pu`. Explicit reciprocal ratio is also rejected in the scan: best reciprocal rows have scores `≈61.175608 pu` or worse. This points to phase-shift unit/sign handling, not a normal-vs-reciprocal transformer ratio mistake. | The diagnostic scan covered `13658` active-power and `9567` reactive-power residual entries. The full shunt scan was interrupted after more than two minutes in local verification; finish the complete Web UI/API artifact run in CI or an environment sized for this case. | Do not classify the present non-convergence as wrong-branch evidence until the import convention is applied and a numerical solution exists. Wrong-branch checks are post-solve plausibility diagnostics, not a transformer-ratio convention test. | Run Web UI/API controls A/B/C with `run_import_profile_only` or equivalent profile-only workflow first, then run NR with Q-limits disabled under the winning convention, and only then repeat with Q-limits enabled to assess active-set stability. |
| `case1951rte.m` | 1,951 buses; large RTE transmission case | Latest documented rectangular run after: `converged=yes`, `iterations=6`, `pv2pq_events=15`, `pv2pq_buses=15`, PF time `≈2.23 s`. Final active PV Q-limit check is `OK`; REF/slack Q-limit diagnostic is `WARN` only. Strict MATPOWER compare still reports `FAIL` with `max|dVm|≈0.01025 pu`, `max|dVa|≈1.64992°`, but the later `status=OK` line is only the deliberately relaxed `1e9` tolerance check. | `standard_matpower_import` + `large_flatstart_dc_blend` + `diagnostic_import_scan`; keep Q-limit diagnostics enabled | Import conventions appear standard: `SHIFT sign=1.0`, `unit=deg`, `ratio=normal`; only 4 non-zero SHIFT entries. Fixed-reference self-check is nearly clean except one dominant active-power residual: `max|ΔP|≈0.144038 pu` (`≈14.404 MW`) at `BUS_I 46`; `max|ΔQ|≈0.000161 pu` (`≈0.016 MVAr`). Solver itself converges robustly, active PV setpoint residual after solve is essentially zero (`≈2.22e-16 pu` over 335 active PV/REF buses). | 15 PV buses switch to PQ due to Q-limits; therefore their original PV voltage setpoints are no longer enforced. The compare now classifies voltage references as `active_pv_imported_setpoint=354`, `final_pq_after_qlimit=15`, `pq_bus_vm=1581`, `ref_slack_imported_setpoint=1`. The largest Vm/angle deviations are dominated by final PQ-after-Q-limit buses and nearby PQ buses, especially around `BUS_I 46`, `817`, `1321`, `1429`, `58`, `1511`, and `908`. The remaining Q-limit violation is REF/slack-only: `BUS_I 1320`, `Qcalc≈1.325 MVAr`, `Qmax=1.000 MVAr`, violation `≈0.325 MVAr`; by design this is reported as WARN and does not trigger PV→PQ switching. There are 20 PV/REF buses without online generators, using `BUS.VM` fallback. 76 branches have `BR_X<0`; signed impedances are preserved. | Low in the documented run. No evidence of wrong Newton branch or low-voltage solution. Current issue is mainly comparison semantics after PV→PQ switching plus a small REF/slack Q-limit WARN, not a solver convergence failure. | Finish compare semantics: split strict compare into active-PV, final-PQ-after-Q-limit, PQ, and REF/slack categories. Consider `WARN` instead of hard `FAIL` when the only relevant voltage deviations are caused by buses that legitimately switched from PV to PQ. Keep REF/slack-only Q-limit violations as WARN by default. Investigate whether a distributed slack / distributed Q option is useful later, but do not mix it into the current cleanup. Reduce debug-log verbosity and consolidate tests to avoid test explosion. |

## YAML profiles

### `standard_matpower_import`

Use for normal MATPOWER reference comparison when imported voltage values are trusted as start values.

```yaml
opt_flatstart: false
matpower_shift_sign: 1.0
matpower_shift_unit: deg
matpower_ratio: normal
bus_shunt_model: admittance
matpower_pv_voltage_source: gen_vg
compare_voltage_reference: hybrid
```

### `large_flatstart_dc_blend`

Use for large cases where a true flat start is required, but a classic flat start may converge to an unstable or historically wrong Newton branch. In the current nested YAML schema this profile uses a DC-angle start with a blended voltage start and keeps wrong-branch detection active.

```yaml
power_flow:
  method: rectangular
  flatstart: true
  tol: 1.0e-5
  max_iter: 80
  autodamp: true
  autodamp_min: 0.05

  wrong_branch_detection: warn   # use fail for strict CI-style rejection
  # Maintainer decision (Issue #196): the wrong-branch rescue retry loop is
  # intentionally out of scope and will not be implemented as part of that
  # work. `wrong_branch_detection: rescue` stays a reserved mode that reports
  # `wrong_branch_rescue_not_implemented` rather than retrying. Detection
  # with full output visibility (report/CSV/console/Web UI/API metadata) plus
  # the APSLF solver as an alternative start/solve path (see
  # docs/src/external_solvers.md) are the supported mitigations for hard
  # flat-start cases.
  wrong_branch_rescue: false
  wrong_branch_min_vm_pu: 0.60
  wrong_branch_max_vm_pu: 1.30
  wrong_branch_max_angle_spread_deg: 180.0
  wrong_branch_max_branch_angle_deg: 90.0
  wrong_branch_min_low_vm_count: 1

  start_mode:
    angle_mode: dc
    voltage_mode: profile_blend
    profile_source: matpower_reference
    start_projection: false
    try_dc_start: true
    try_blend_scan: false
    branch_guard: true
    measure_candidates: true
    accept_unmeasured_dc_start: false
    reuse_import_data: true
    dc_angle_limit_deg: 60.0

matpower_import:
  shift_sign: 1.0
  shift_unit: deg
  ratio: normal
  bus_shunt_model: admittance
  pv_voltage_source: gen_vg
  compare_voltage_reference: hybrid
```

### `activsg10k_wrong_branch_historical_scan`

Use only for manual diagnostics when trying to reproduce the historical ACTIVSg10k wrong-branch report. This is **not** a production recommendation. Current v0.8.2 tests did not find a natural wrong-branch reproducer with normal shift/ratio conventions and normal `wrong_branch_*` thresholds.

```yaml
power_flow:
  method: rectangular
  flatstart: true
  tol: 1.0e-5
  max_iter: 90
  autodamp: false

  wrong_branch_detection: warn
  wrong_branch_min_vm_pu: 0.60
  wrong_branch_max_vm_pu: 1.30
  wrong_branch_max_angle_spread_deg: 180.0
  wrong_branch_max_branch_angle_deg: 90.0

  start_mode:
    angle_mode: classic
    voltage_mode: profile_blend
    profile_source: matpower_reference
    start_projection: false
    try_dc_start: false
    try_blend_scan: true
    branch_guard: true
    measure_candidates: true
    reuse_import_data: true

matpower_import:
  auto_profile: off
  pv_voltage_source: gen_vg
  compare_voltage_reference: hybrid
  bus_shunt_model: admittance
  shift_unit: deg
  shift_sign: 1.0
  ratio: normal
```

Observed v0.8.2 result with this default-near reproducer candidate:

```text
status                         = converged
final_converged                = true
wrong_branch_status            = ok
min_vm_pu                      ≈ 0.9235
max_vm_pu                      ≈ 1.1087
wrong_branch_angle_spread_deg  ≈ 107.1
wrong_branch_worst_branch_deg  ≈ 22.94
```

A more degraded `classic/classic` voltage-angle start does not currently reproduce a converged wrong branch either; it fails before the post-convergence wrong-branch check with `nr_mismatch_not_converged_active_set_unstable`.

### `activsg10k_forced_wrong_branch_detection_test`

Use only to validate that wrong-branch configuration forwarding and final-status handling work. This profile intentionally tightens the branch-angle threshold below the observed healthy worst-branch angle; it is not a physical wrong-branch reproducer.

```yaml
power_flow:
  wrong_branch_detection: warn   # change to fail for final rejection
  wrong_branch_max_branch_angle_deg: 20.0
```

Expected current result on the healthy ACTIVSg10k solution:

```text
warn mode:
  wrong_branch_status            = warn
  wrong_branch_reason            = branch_angle_exceeded
  wrong_branch_branch_violations = 3

fail mode:
  numerical_converged            = true
  final_converged                = false
  status                         = wrong_branch_detected
  reason                         = wrong_branch_detected
```

### `activsg25k_q_limit_guard`

Use for the documented successful `case_ACTIVSg25k.m` run. The key point is a robust DC-angle/blended-voltage flat start with start projection disabled, plus a strong Q-limit guard for zero/narrow-Q PV buses and PV buses whose calculated Q violates limits after the first eligible Q-limit check. Use `tol: 1.0e-5` for large benchmark validation unless strict validation is explicitly requested.

```yaml
opt_flatstart: true
flatstart_angle_mode: dc
flatstart_voltage_mode: profile_blend
profile_source: matpower_reference
start_projection: false
start_projection_try_blend_scan: false
start_projection_try_dc_start: true
start_projection_dc_angle_limit_deg: 70.0
start_projection_blend_lambdas:
  - 0.1
  - 0.25
  - 0.5
  - 0.75

autodamp: true
autodamp_min: 0.05
max_ite: 80
tol: 1.0e-5

ignore_q_limits: false
qlimit_start_iter: 3
qlimit_start_mode: iteration_or_auto
qlimit_auto_q_delta_pu: 0.0001
q_hyst_pu: 0.01
cooldown_iters: 1

qlimit_guard: true
qlimit_guard_min_q_range_pu: 0.02
qlimit_guard_zero_range_mode: lock_pq
qlimit_guard_narrow_range_mode: lock_pq
qlimit_guard_violation_mode: lock_pq
qlimit_guard_violation_threshold_pu: 0.0001
qlimit_guard_max_switches: 3
qlimit_guard_freeze_after_repeated_switching: true
qlimit_guard_accept_bounded_violations: false
qlimit_guard_max_remaining_violations: 0
qlimit_guard_log: true

matpower_shift_sign: 1.0
matpower_shift_unit: deg
matpower_ratio: normal
bus_shunt_model: admittance
matpower_pv_voltage_source: gen_vg
compare_voltage_reference: imported_setpoint
```

### `pegase_shift_rad_opposite`

Use for PEGASE cases where the branch-shift scan indicates that the stored shift values behave like radians with opposite sign.

```yaml
matpower_shift_sign: -1.0
matpower_shift_unit: rad
matpower_ratio: normal
bus_shunt_model: admittance
matpower_pv_voltage_source: gen_vg
compare_voltage_reference: hybrid
```

### `pegase9241_reference_diagnostic`

Use for the documented `case9241pegase.m` diagnostic runs while isolating import/reference consistency effects before enabling Q-limit behaviour again. The latest log used a DC-angle/blended flat start with start projection.

```yaml
opt_flatstart: true
flatstart_angle_mode: dc
flatstart_voltage_mode: profile_blend
profile_source: matpower_reference
start_projection: true
start_projection_try_blend_scan: false
start_projection_try_dc_start: true
ignore_q_limits: true
matpower_shift_sign: -1.0
matpower_shift_unit: rad
matpower_ratio: normal
bus_shunt_model: admittance
matpower_pv_voltage_source: gen_vg
compare_voltage_reference: hybrid
diagnose_matpower_reference: true
diagnose_branch_shift_conventions: true
diagnose_pv_voltage_references: true
show_diff: true
log_effective_config: true
```

### `diagnostic_import_scan`

Use when checking whether deviations originate from import conventions rather than the Newton solver.

```yaml
diagnose_matpower_reference: true
diagnose_branch_shift_conventions: true
diagnose_pv_voltage_references: true
show_diff: true
log_effective_config: true
```

## Diagnostic classification notes

Use the following short labels in the case matrix where useful:

| Label | Meaning |
|---|---|
| `ok` | Solver result and reference comparison are within the configured tolerances. |
| `wrong_branch_risk` | Newton may converge formally but to an implausible low-voltage or wrong-branch solution. |
| `reference_data_inconsistent` | The stored MATPOWER `VM`/`VA` reference is not power-balanced against the imported branch data in the fixed-reference self-check. |
| `import_convention_open` | SHIFT, TAP, ratio, shunt, or voltage-reference conventions still need verification. |
| `pv_limit_sensitive` | Result depends materially on PV→PQ switching, Q-limit hysteresis, or zero/near-zero Q limits. |
| `relaxed_compare_only` | A later comparison with artificially huge tolerances, such as `1e9`, is only a diagnostic sanity check and must not be counted as a real pass. |
| `pegase_shift_rad_opposite` | PEGASE-style shift interpretation appears necessary: `matpower_shift_sign=-1.0`, `matpower_shift_unit=rad`. |
| `ref_slack_q_warn` | Remaining Q-limit violation is on the REF/slack bus only and is reported as WARN by default. |
| `final_pq_after_qlimit` | Bus was originally PV but switched to PQ due to Q-limits; its original PV voltage setpoint is no longer enforced. |
| `qlimit_guard_stabilized` | Q-limit guard pre-locking and/or violation locking stabilizes many zero/narrow-Q PV buses and allows the active set to converge. |
| `compact_console_needed` | Run is technically successful, but console output is too verbose; detailed diagnostics should remain in the logfile. |
| `historical_wrong_branch_not_reproduced` | A previously observed wrong-branch or low-voltage branch remains documented, but current solver/config scans no longer reproduce it naturally. |
| `forced_wrong_branch_detection_test` | The threshold is intentionally tightened to validate wrong-branch forwarding/final-status handling; this is not a physical wrong-branch reproducer. |

## Notes for extending this table

Keep each row short. Prefer numeric evidence over prose:

- convergence status and iteration count,
- maximum voltage magnitude deviation,
- maximum angle deviation,
- number of PV→PQ switches,
- number of PV/REF buses without online generators,
- whether zero Q-limits occur,
- whether a classic flat start can enter a second Newton solution branch,
- whether the branch-shift scan indicates a non-standard convention.

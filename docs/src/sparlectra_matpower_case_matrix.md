# Sparlectra.jl MATPOWER Case Diagnostics Matrix

Date: 2026-05-16  
Sparlectra versions covered: 0.7.7–0.7.8  
Purpose: compact tracking table for MATPOWER import, reference-data consistency, and solver diagnostics.

> **Data and license note**  
> The MATPOWER case files referenced here are publicly available from the MATPOWER project or related public benchmark-data sources. This document does not redistribute those case files. It only records diagnostic observations from local Sparlectra.jl runs. MATPOWER itself is distributed under the 3-clause BSD license, but the case files included with MATPOWER may have their own origin, citation requests, or permission conditions. Before redistributing case data or derived data sets, check the header comments of the corresponding case file and the applicable MATPOWER documentation.
>
> **Non-binding diagnostic note**  
> The observations in this matrix are engineering diagnostics from specific Sparlectra.jl runs and configuration states. They are not legal advice, not a formal validation certificate, and not a legally binding assessment of the correctness, licensing status, or benchmark suitability of the referenced data sets.

This note is intended as a living diagnostic matrix. Add one row per MATPOWER case and keep the findings short.

The matrix distinguishes solver convergence failures, wrong Newton solution branches, Sparlectra import/Y-bus mismatches, and inconsistent stored MATPOWER reference data. When a large fixed-reference residual is already reproducible with a raw MATPOWER-style Y-bus using the stored `VM`/`VA` values and branch data, classify the case as a reference-data consistency diagnostic rather than as a clean solver benchmark failure.

## Latest addition: `case_ACTIVSg25k.m`

`case_ACTIVSg25k.m` now converges successfully with the strengthened Q-limit guard profile. The final documented run reports `numerical_solution=OK`, `q_limit_active_set=OK`, `final_converged=true`, `iterations=8`, `pv2pq_events=2004`, `guarded narrow-Q PV buses=960`, and `compare=OK`. The remaining active PV/REF Q-limit violations are zero. This confirms that the earlier failure mode was primarily a Q-limit active-set stabilization issue, not a classical Newton divergence.

## Case matrix

| Case | Size / type | Status in Sparlectra | Recommended YAML profile | Main findings | Case-specific anomalies / deviations | Wrong-branch risk | Open checks |
|---|---:|---|---|---|---|---|---|
| `case_ACTIVSg10k.m` | 10,000 buses; large synthetic transmission case | Converges with DC-angle flat start and blended voltage start. Latest documented run: `converged=yes`, `iterations=9`, `pv2pq_events=657`, `pv2pq_buses=640`, compare `status=OK`. | `large_flatstart_dc_blend` | MATPOWER import conventions confirmed: `shift_sign=1.0`, `shift_unit=deg`, `ratio=normal`, `bus_shunt_model=admittance`. Branch-shift scan strongly rejects opposite sign, radians, and reciprocal ratio. Bus shunts are required for reactive balance. | 273 PV/REF buses have no online generator; fallback to `BUS.VM` is required. Many PV buses have `Qmin=0` or `Qmax=0`. Large number of PV→PQ switches when Q-limits are enforced. Final active PV/REF setpoint residual is essentially zero for active PV buses, but hybrid comparison can show larger deviations where reference is `BUS.VM` or where buses switched to PQ. | High with classic flat start. A previous classic flat-start run converged formally to a low-voltage / wrong-branch solution with `max|dVm|≈1.0187 pu`, `max|dVa|≈90.73°`, `minVm=0.0 pu`. DC angle start removes this behaviour in the tested run. | Improve compact Q-limit reporting. Keep diagnostics for PV/REF buses without online generators. Track whether many zero Q-limits are data artefacts or valid MATPOWER case data. |
| `case_ACTIVSg25k.m` | 25,000 buses; large ACTIVSg synthetic transmission case; 32,230 branches | Rectangular NR converges with strengthened Q-limit guard. Successful run: `numerical_solution=OK`, `q_limit_active_set=OK`, `final_converged=true`, `status=converged`, `iterations=8`, PF time `≈74.55 s`, total example runtime `≈85.11 s`, final mismatch `≈4.651e-08`. | `activsg25k_q_limit_guard` + `large_flatstart_dc_blend` | MATPOWER fixed-reference self-check is good for this size: `max|ΔP|≈7.118 MW`, `max|ΔQ|≈4.611 MVAr`. Auto-profile detects `2034/3779` online generators with zero or narrow Q range. Q-limit guard locks `960` narrow-range PV buses as PQ before rectangular NR. Active set converges cleanly: `PV violations=0`, `REF violations=0`, final `PV/REF Q-limit check=OK`. Switching statistics: `pv2pq_events=2004`, `pv2pq_buses=1996`, `oscillating_buses=0`. Compare against imported setpoints is `OK`: `max|dVm|≈0.04291 pu`, aligned `max|dVa|≈0.3683°`, slack Δ `≈0.0°`. | `2753/3235` PV/REF buses have online `GEN.VG` targets; `502` differ from `BUS.VM` by more than `1e-4 pu`. `482` PV/REF buses have no online generator and use `BUS.VM` fallback. Negative branch impedance is present and preserved: `BR_R<0` on `447` rows, `BR_X<0` on `503` rows, both negative on `447` rows. Compare reference kinds: `active_pv_imported_setpoint=1238`, `final_pq_after_qlimit=1996`, `pq_bus_vm=21765`, `ref_slack_imported_setpoint=1`. | Low with the successful profile. Previous failures were active-set stabilization failures caused by many zero/narrow-Q generators and infeasible PV setpoints, not classical Newton wrong-branch convergence. | Console output is still too verbose: keep full auto-profile, negative-branch diagnostics, PV voltage diagnostics, and per-bus PV→PQ events in the logfile; show only compact statistics on console. Consider whether `qlimit_guard_violation_mode=lock_pq` should become an auto-profile recommendation for large narrow-Q cases. |
| `case300.m` | 300 buses; medium test case | Rectangular NR converges in about 5 iterations with 3 PV→PQ switches. Final voltage magnitudes are close, but angle comparison fails mainly around the `BUS_I 196 / 2040` area. | `large_flatstart_dc_blend` + `diagnostic_import_scan` | The MATPOWER fixed-reference self-check is not power-balanced around `BUS_I 196 / 2040`. The dominant mismatch is already reproducible with a MATPOWER-style Y-bus using the stored `VM`/`VA` values and raw branch data. Main residual: about 926.9 MW at `BUS_I 2040` and `BUS_I 196`. | Not classified as a Sparlectra Newton wrong-branch issue. The dominant mismatch appears to be caused by the stored case reference angles together with a low-reactance active branch `196 -> 2040` (`x = 0.02`, `TAP = 1`, `SHIFT = 0`). Treat this case as a reference-data consistency diagnostic, not as a clean pass/fail solver benchmark. | Low for the documented run; angle-reference comparison is dominated by the stored fixed-reference inconsistency rather than a known wrong branch. | Keep the fixed-reference diagnostic classification visible when updating MATPOWER comparison tolerances or profiles. Keep the explicit `TAP = 1` preservation regression, but do not treat lost nominal TAP as the remaining root cause. |
| `case145.m` | 145 buses; medium/small test case | Rectangular NR converges: `converged=yes`, `iterations=6`, `pv2pq_events=1`, compare `status=OK`. | `large_flatstart_dc_blend` + `diagnostic_import_scan` | No branch-shift entries. Standard degree/radian/sign scans are irrelevant because `SHIFT=0`. Reciprocal ratio is rejected by diagnostics. Bus shunts are important: disabling them causes very large P/Q residuals. | One relevant PV→PQ switch: `BUS_I 104`, where `GEN.VG=1.045` differs from `BUS.VM≈1.0059`; after switching, comparison against imported setpoint shows `dVm≈-0.03909 pu`, but the hybrid compare is still within configured tolerances. Several PV buses have `Qmin=0`. | Low in the documented run. No evidence of a wrong Newton branch; angle deviation is very small (`max|dVa|≈0.0051°`). | Track whether `BUS_I 104` is a valid Q-limit-driven PV→PQ case. Keep zero-Q-limit and shunt sensitivity visible in future comparisons. |
| `case1354pegase.m` | 1,354 buses; PEGASE case | With standard MATPOWER-style shift interpretation (`sign=1.0`, `unit=deg`) the compare fails. With PEGASE-style interpretation `matpower_shift_sign=-1.0`, `matpower_shift_unit=rad`, rectangular NR converges in 3 iterations with `pv2pq_events=0`; compare `status=OK`, `max|dVm|≈0.00254 pu`, `max|dVa|≈0.47247°`. | `pegase_shift_rad_opposite` + `large_flatstart_dc_blend` + `diagnostic_import_scan` | Branch-shift scan identifies `opposite sign, radians` as the only plausible convention. It reduces the fixed-reference active-power residual from about 1299.8 MW to about 38.8 MW. Bus shunts remain important for Q balance. | Six non-zero shift entries. Remaining fixed-reference mismatch is dominated by `BUS_I 58` and `BUS_I 6153`. No PV/REF buses without online generators. Active PV/REF setpoint residual is essentially zero after solve. | Low in the documented successful run. The earlier failed compare was mainly an import-convention issue, not a demonstrated wrong Newton branch. | Add an automatic PEGASE shift-convention recommendation when branch-shift scan strongly prefers `sign=-1.0`, `unit=rad`. Keep residual classification visible because the fixed-reference check is not perfectly zero even with the preferred convention. |
| `case2869pegase.m` | 2,869 buses; PEGASE case | With PEGASE-style shift interpretation `matpower_shift_sign=-1.0`, `matpower_shift_unit=rad`, rectangular NR converges robustly: `converged=yes`, `iterations=3`, `pv2pq_events=0`, compare `status=OK`, `max|dVm|≈0.00129 pu`, `max|dVa|≈0.25599°`. | `pegase_shift_rad_opposite` + `large_flatstart_dc_blend` + `diagnostic_import_scan` | Branch-shift scan strongly rejects standard degree-based interpretations. The preferred `opposite sign, radians` convention gives the best fixed-reference residuals: `max|dP|≈20.55 MW`, `max|dQ|≈50.54 MVAr`; standard MATPOWER sign/degrees gives `max|dP|≈4202.63 MW`. Bus shunts are important for Q balance: disabling them raises `max|dQ|` to about 391.2 MVAr. | Twelve non-zero branch SHIFT entries. No PV→PQ switches occurred in the documented run despite many PV buses and Q-limit entries. No PV/REF buses without online generators. Active PV/REF setpoint residual is essentially zero (`≈9.6e-12 pu` over 510 buses). Top remaining compare deviations are small and concentrated around `BUS_I 1541`, `221`, `58`, and `6153`. | Low in the documented run. No evidence of a wrong Newton branch when the PEGASE shift convention is used. | Add automatic detection/reporting for PEGASE-like radian shift values with opposite sign. Keep this case as a positive regression for `sign=-1.0`, `unit=rad`, `ratio=normal`, `bus_shunt_model=admittance`. |
| `case9241pegase.m` | 9,241 buses; large PEGASE case | Latest run `run_case9241pegase.m_20260515_132616.log`: with `opt_flatstart=true`, DC-angle start, `bus_vm_va_blend`, start projection enabled, PEGASE-style `matpower_shift_sign=-1.0`, `matpower_shift_unit=rad`, `matpower_ratio=normal`, and `ignore_q_limits=true`, rectangular NR converges: `converged=yes`, `iterations=4`, `pv2pq_events=0`, `pv2pq_buses=0`, PF time `≈10.43 s`. Strict compare still fails by angle tolerance: `max|dVm|≈0.02662 pu`, `max|dVa|≈1.48863°` with `tol_va=0.99°`; the later `status=OK` line is only the deliberately relaxed `1e9` tolerance check. | `pegase_shift_rad_opposite` + `large_flatstart_dc_blend` + `diagnostic_import_scan`; for this diagnostic run also `ignore_q_limits=true` | Solver convergence and PV setpoint tracking are good: post-solve active PV/REF residual `max|Vm_calc-imported_Vset|≈1.21e-13 pu` over 1445 buses and no PV/REF buses without online generators. The fixed-reference self-check remains inconsistent against the imported model: `max|ΔP|≈2.71965 pu` (`≈271.97 MW`) and `max|ΔQ|≈3.92020 pu` (`≈392.02 MVAr`). The branch-shift scan still selects `opposite sign + radians + normal ratio` as best global convention; disabling bus shunts does not improve P and worsens Q. | Remaining deviations are concentrated in local residual clusters and special branch regions. Residual clustering reports 31 clusters above 10 MW/MVAr; the dominant clusters include `5491/5673`, `7934/6621`, and `2309/8414`. Notable branches include rows `14747/14748` (`7934 -> 6621`, `TAP=0.976337`, `SHIFT=-0.045364`), row `14763` (`2309 -> 8414`, `TAP=0.946281`, `SHIFT=-0.010988`), row `14766` (`503 -> 8414`, `TAP=0.993682`, `SHIFT=+0.080858`), and negative-R high-residual branches `14764` (`2309 -> 8414`, `BR_R=-0.011591`) and `14767` (`503 -> 8414`, `BR_R=-0.010888`). The log reports 75 branches with `BR_R<0` and 16 with `BR_X<0`; these are preserved, not clipped. | Low in the documented run. No evidence of a low-voltage/wrong-branch solution; this is an import/reference-consistency diagnostic, not a clean Newton failure. | Inspect branch-local TAP/SHIFT conventions and residual attribution for the top clusters, especially duplicated transformer rows and adjacent negative-R branches. Add branch-group aggregation for ordinary lines, TAP-only, SHIFT-only, TAP+SHIFT, and negative R/X. Re-test with Q-limits enabled only after the fixed-reference mismatch is understood. |
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

Use for large cases where a true flat start is required, but a classic flat start may converge to a wrong branch.

```yaml
opt_flatstart: true
flatstart_angle_mode: dc
flatstart_voltage_mode: bus_vm_va_blend
start_projection: true
start_projection_try_blend_scan: false
start_projection_try_dc_start: true
autodamp: true
max_ite: 80
matpower_shift_sign: 1.0
matpower_shift_unit: deg
matpower_ratio: normal
bus_shunt_model: admittance
matpower_pv_voltage_source: gen_vg
compare_voltage_reference: hybrid
wrong_branch_detection: true
wrong_branch_rescue: true
```

### `activsg25k_q_limit_guard`

Use for the documented successful `case_ACTIVSg25k.m` run. The key point is not only a robust Newton start, but also a strong Q-limit guard for zero/narrow-Q PV buses and PV buses whose calculated Q violates limits after the first eligible Q-limit check.

```yaml
opt_flatstart: true
flatstart_angle_mode: dc
flatstart_voltage_mode: bus_vm_va_blend
start_projection: true
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
flatstart_voltage_mode: bus_vm_va_blend
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

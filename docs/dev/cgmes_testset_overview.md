# CGMES test-set overview — measured state

Measured 2026-07-30 on `dev/r0.8.16` with the full sweep in
`examples/run_cgmes_suite.jl` (import → `runpf!`, islands enabled, default
Q-limit handling, `max_iter = 60`, `tol = 1e-8` → `compareWithSV`; on
non-convergence one retry from a flat start on a fresh import, reported as
`ok (flat)`). Regenerate with:

```bash
julia --project=. examples/run_cgmes_suite.jl
```

Voltage/angle columns compare the solved state against the delivery's **own
SV profile** (isolated/de-energized buses excluded); flow columns compare
against `SvPowerFlow` per terminal. `iso` counts isolated buses after import.

## Result table

| Case | CGMES | Buses | iso | Solves? | it | max dvm (pu) | rms dvm | max dva | n flows | rms dp | max dp |
|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| MicroGrid BE | 2.4.15 | 12 | 0 | ✅ | 4 | 0.00007 | 0.00003 | 0.04° | 12 | 0.005 MW | 0.0 MW |
| MicroGrid NL | 2.4.15 | 10 | 0 | ✅ | 3 | 0.00012 | 0.00010 | 0.02° | 11 | 0.966 MW | 2.3 MW |
| MicroGrid Assembled | 2.4.15 | 17 | 0 | ✅ | 4 | 0.00009 | 0.00005 | 0.02° | 13 | 0.886 MW | 2.3 MW |
| SmallGrid BusBranch | 2.4.15 | 118 | 0 | ✅ | 3 | 0.00014 | 0.00005 | 0.01° | 139 | 0.009 MW | 0.1 MW |
| SmallGrid NodeBreaker | 2.4.15 | 118 | 0 | ✅ | 3 | 0.00014 | 0.00005 | 0.01° | 139 | 0.009 MW | 0.1 MW |
| MiniGrid BusBranch | 2.4.15 | 15 | 2 | ✅ | 2 | 0.00001 | 0.00001 | 0.00° | 3 | 0.000 MW | 0.0 MW |
| MiniGrid NodeBreaker | 2.4.15 | 15 | 2 | ✅ | 2 | 0.00001 | 0.00001 | 0.00° | 3 | 0.000 MW | 0.0 MW |
| FullGrid BusBranch | 2.4.15 | 25 | 1 | ✅ (flat) | 6 | (0.067) | (0.035) | (14.5°) | 15 | (101.8 MW) | (315.5 MW) |
| FullGrid NodeBreaker | 2.4.15 | 25 | 1 | ✅ (flat) | 6 | (0.067) | (0.035) | (14.5°) | 16 | (98.6 MW) | (315.5 MW) |
| **RealGrid** | 2.4.15 | 6209 | 158 | ✅ | 4 | 0.03118 | 0.00354 | 0.96° | 6861 | 0.026 MW | 2.1 MW |
| PST PTChLin PTE1 | 2.4.15 | 2 | 0 | ✅ | 2 | 0.00000 | 0.00000 | 0.00° | 2 | 0.000 MW | 0.0 MW |
| PST PTChLin PTE2 | 2.4.15 | 2 | 0 | ✅ | 4 | 0.00099 | 0.00070 | 0.29° | 2 | 0.559 MW | 0.8 MW |
| PST PTChTab PTE2 | 2.4.15 | 2 | 0 | ✅ | 4 | 0.00099 | 0.00070 | 0.29° | 2 | 0.559 MW | 0.8 MW |
| TransformerLineTest | 2.4.15 | 8 | 0 | ✅ | 2 | 0.00000 | 0.00000 | 0.00° | 8 | 0.000 MW | 0.0 MW |
| ExplicitLoadFlow | 2.4.15 | 2 | 0 | ✅ | 5 | 0.00227 | 0.00160 | 3.77° | 3 | 0.010 MW | 0.0 MW |
| **ReliCap Svedala** | 3.0 | 132 | 14 | ✅ | 7 | 0.08657 | 0.03613 | 3.36° | 418 | 1.603 MW | 11.0 MW |
| ReliCap Svedala+Neighbours | 3.0 | 273 | 14 | ❌ no conv | 60 | — | — | — | — | — | — |
| ReliCap Belgovia | 3.0 | 16 | 0 | ✅ | 4 | 0.01430 | 0.00792 | 0.22° | 46 | 0.235 MW | 0.7 MW |
| ReliCap Britheim | 3.0 | 6 | 0 | ✅ | 4 | 0.00000 | 0.00000 | 0.00° | 11 | 0.000 MW | 0.0 MW |
| ReliCap Espheim | 3.0 | 131 | 0 | ✅ | 13 | 0.16125 | 0.07050 | 24.25° | 549 | 8.661 MW | 89.9 MW |
| ReliCap Galia | 3.0 | 12 | 0 | ✅ | 4 | 0.00250 | 0.00187 | 0.37° | 33 | 0.034 MW | 0.1 MW |
| ReliCap Nordheim | 3.0 | 1 | 1 | ❌ no slack | — | — | — | — | — | — | — |
| ReliCap Portheim | 3.0 | 8 | 6 | ✅ | 6 | 0.76367 | 0.75687 | 1.48° | 6 | 0.008 MW | 0.0 MW |
| ReliCap CGM (all areas) | 3.0 | 294 | 21 | ❌ island 1 | — | — | — | — | — | — | — |

## Per-case notes

**Conformity sets (2.4.15).** MicroGrid (BE/NL/Assembled), SmallGrid (both
topology variants), the three PST/PSEI toys, TransformerLineTest and
ExplicitLoadFlow are SV-tight. The PTE2 rows deviate by ≈0.001 pu / 0.29° *by
design*: Sparlectra folds an end-2 tap angle with end-referral semantics
(θ_eff = θ1 − θ2, pinned by RealGrid), while these two conformity toys expect
the unflipped angle — documented in `docs/src/cgmes_import.md`. MicroGrid
NL/Assembled carry a ~1 MW flow deviation at a single machine terminal, a
data-level SSH↔SV inconsistency of the set.

**MiniGrid** is SV-tight (2 iterations, max dvm 1e-5). The residual reported
in the 2026-07-29 measurement (0.031 pu, flow rms 5.24 MW) was **not** a
data inconsistency as first noted — it was the then-unmapped
`AsynchronousMachine` load: MiniGrid carries three induction motors (M3,
M2a, M2b — 9 MW / 5 MVAr at bus 7), mapped since #294 point 6 as Stage-0 PQ
operating points. Two buses import isolated (sourceless stubs).

**FullGrid — solves from a flat start.** It is the import/export
*completeness* configuration: every CGMES class appears at least once (all
four PST types as four parallel transformers on one branch, VSC+CSC HVDC,
SVC, `AsynchronousMachine`, nonlinear shunts), with the systematic **X.99
placeholder scheme** filling otherwise-unused attributes (tabular PST table
row `ratio 9.99 / angle 0.99°` with `r = x = 9.99`, a
`NonlinearShuntCompensatorPoint` with `b = g = 0.99 S` = a 50-GW shunt at
225 kV, switch `ratedCurrent 999.99`, ENI limits 9.99/99.99 — independently
confirmed against pandapower's mirror of the v1 set). Two importer
plausibility guards (see `cgmes_import.md`) detect the tap-row and shunt
placeholders and keep them out of the model with warnings; the
`AsynchronousMachine` maps as a Stage-0 PQ point. With that, **the network
solves in 6 iterations from a flat start** (voltages 0.954–1.047 pu). The
shipped SV profile remains internally inconsistent — it declares a 14.5°
angle jump across a 0.3 Ω line (≈400 pu start residual), which is why the
default SV-based start fails and why the SV-deviation columns are shown in
parentheses: they measure the SV profile's own inconsistency, not solver
error. FullGrid's role stays reader completeness, short-circuit source data,
and — once #279 lands — the natural roundtrip-export fixture.

**RealGrid** (6209 buses, the only real-exchange-sized delivery): solves in 4
iterations from SV start. After the tabular-PST wiring (#294 point 2) the old
forensics verdict "local voltage collapse around the 63 kV level" is
**revised**: the collapse probe bus `1_38470111` sits at vm 1.0077 vs SV
1.0078, max dvm dropped 1.01 → 0.031, buses with |dvm| > 0.01 dropped
241 → 169. The remaining 0.01–0.03 band is genuine SSH↔SV residue. Q-limits
ON is bit-identical (active set consistent at SV); distributed slack is
unnecessary (λ_P ≈ +4.9 MW) and slightly worsens flows (rms dp 0.026 → 0.085
MW) — expected, not a regression: the SV snapshot was produced with a single
slack, so an SV comparison under `pg_weighted` measures the slack-distribution
difference (λ_P spread over 368 participants), not solver error; the
bit-identical voltages confirm only the P split moves. RealGrid contains **no** remote-regulating machines (528 machine
RegulatingControls, all with local terminals — the "10 controls" in the
package documentation refer to off-nominal target *values*), so
`machine_control = true` is a no-op here. Measurement script:
`examples/experimental/val_realgrid_remeasure.jl`.

*Fixed-reference self-check attribution (2026-07-30, `run_fixed_reference_self_check`
+ `self_check_residuals.csv`).* At the delivery's own SV state (full coverage:
0 of 6209 buses without `SvVoltage`) the residual is `inf`-norm **5.06 pu**,
totals sum|P| 48.4 / sum|Q| 21.6 pu. Attribution: 90 % of the Q residual
(19.5 pu) and 97 % of the P residual sit at the 1535 transformer-adjacent
buses (≈29× the per-bus level elsewhere); shunt buses are secondary (5.1 pu at
457 buses). The worst residual is one equal-and-opposite pair — buses
`Bus_2335_63`/`Bus_4799_63`, P ±5.06, Q ±1.65 pu across branch
`B_2WT_63_4799_2335` — which is exactly the tabular-PST winding of phase tap
`892946845` (step 12, shift 4.305°) whose **per-step r/x/g/b corrections are
read but not applied** (importer notice, 9 tables total). Interpretation: the
imported model is consistent at SV except for the known tabular-PST
correction gap; the divergence seen in WebUI runs is a start-state topic, not
a Q-side conversion error.

*Why WebUI runs diverge while the remeasure script solves (2026-07-30).* The
WebUI base `configuration.yaml` (and `src/configuration.yaml.example`) carries
the legacy top-level `power_flow.flatstart: true`; `PowerFlowConfig` merges
the power_flow section into the start_mode raw dict, so this silently forces a
**flat start** on CGMES runs — the imported SV voltages never reach the
solver, NR runs into the 435-PV Q-limit active-set cascade and diverges
(mismatch 2220 → 240 324, run `6412ae39`). From SV start
(`flatstart: false`) with the YAML default `autodamp: true` the solve stalls
at mismatch 4.5e-3; with `autodamp: false` it **converges** (58 iterations
with Q-limits ON, final mismatch 1.9e-7; the remeasure script's programmatic
config — plain Newton — needs 4). The fixed-reference self-check now forces
`flatstart=false` itself, so Diagnose runs measure the SV state regardless of
the base configuration.

**ReliCapGrid family (CGMES 3.0, second source).** Svedala solves with
multi-island references (14 isolated buses are parked out-of-service units and
stubs; the five implausible 0.001 kV regulation targets belong to units that
are out of service and are skipped before the setpoint check). The residual
0.087 pu / 11 MW reflects the delivery's own SSH↔SV gap. Belgovia, Britheim,
Galia are clean small areas. **Espheim** and **Portheim** solve but deviate
strongly (0.16 pu/24° resp. 0.76 pu on the two non-isolated buses) — both
sides of their shared border are modelled with single-sided equivalents whose
TN identities do not line up (the 3.0 split-TN border documented in
`cgmes_import.md`); standalone these areas lack the counterpart flows the SV
state was computed with. **Nordheim** standalone is a single HVDC-converter
bus with no AC source — no slack, honestly unsolvable. The **CGM** (all areas
assembled) fails in its main island under hull Q-limits: the unlimited
solution differs from SV (single-sided EIs, SM_NO1/G2 dispatch deltas), its
violations are not enforceable — the honest non-convergence analyzed in the
Q-limit sessions (task plan, 2026-07-30). **Svedala+Neighbours** shows the
same signature as the CGM main island.

## Open threads fed by this table

- Svedala+Neighbours / CGM island 1: same bracket as the Q-limit analysis —
  closing it further is #192 follow-up work (SD–BO deviation), not importer
  work.
- Espheim/Portheim standalone borders: expected artifact of single-sided
  3.0 deliveries; only the assembled CGM can resolve them, see above.

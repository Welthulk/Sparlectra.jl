# [State-Estimation Configuration](@id se-config)

The estimator reads every setting from the active configuration:
`runse!(net)`, `runse_diagnostics`, `validate_measurements` and
`validate_topology` take no settings.
`with_state_estimation_config(max_iter = 12) do runse!(net) end` installs
other values for one call and restores the previous configuration
afterwards. The service and the Web UI install the effective configuration
of the case (general file, case file, form values) the same way.

`flatstart`, `robust_mode`, `robust`, `k_eliminate`, `k_suppress`,
`max_eliminations`, `topology_precheck`, `report_residual_correlation`,
`zero_injection` and `sigma_by_voltage_level` are case-scope keys like `power_flow.solver`: a case configuration file may
carry them (Web UI: **Save settings for this case** on the State Estimation
page or the Settings page), and the State Estimation form shows the
effective value (case file, else the general file, else the default below).
The other keys are installation-wide: this file or the API.

Mechanisms: [Bad-data thresholds and robust modes](@ref se-bad-data),
[Sequential elimination budget](@ref se-max-eliminations),
[Residual correlations](@ref se-k-report), [Observability](observability.md)
(rank tolerance, rank source, critical measurements) and
[Topology validation](state_estimation_extensions.md#Topology-validation).

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `state_estimation.enabled` | Bool | `true` | `true`, `false` | Enables the SE stage. |
| `state_estimation.method` | Symbol/String | `wls` | `wls` | Estimation algorithm; other values fail validation. |
| `state_estimation.tol` | Float64 | `1e-6` | positive real | Convergence tolerance on the state step. The default is the finite-difference noise floor (`jac_eps`); a tighter value is unreachable and the run raises it with a log line. |
| `state_estimation.max_iter` | Int | `50` | positive integer | Iteration cap; the service and the Web UI form use the same value. A diverging run does not converge with a higher cap either. |
| `state_estimation.flatstart` | Bool | `true` | `true`, `false` | Start from a flat state instead of the model's current voltages. Switch off when a trusted start state exists. |
| `state_estimation.jac_eps` | Float64 | `1e-6` | positive real | Perturbation of the finite-difference Jacobian; sets the attainable `tol`. |
| `state_estimation.update_net` | Bool | `true` | `true`, `false` | Write the estimated state back to the network, so a following power flow starts from it. |
| `state_estimation.pmu_ref_offset` | Symbol/String | `auto` | `auto`, `off` | `auto` estimates the angle offset between the PMU time base and the slack reference as one extra state when `VaMeas` rows exist (the offset is near 0 when the bases coincide); `off` takes PMU angles as slack-referenced. Keyword `pmuRefOffset`. |
| `state_estimation.imag_activation_iteration` | Int | `2` | positive integer | First WLS iteration in which current-magnitude rows (`ImagMeas`) enter the update, a flat-start protection; no effect without them. Keyword `imagActivationIteration`. |
| `state_estimation.ia_current_floor_A` | Float64 | `10.0` | positive real | Activity gate of a current-angle row (`IaMeas`) without a paired magnitude: active while the predicted current exceeds 3 times this floor (ampere); with a paired `ImagMeas` its 3 sigma apply. Gating notes in the diagnostics. |
| `state_estimation.report_residual_correlation` | Bool | `true` | `true`, `false` | Residual-correlation columns in the bad-data report (max abs k per row, warning above 1/sqrt(2)); one dense m x m pass in the diagnostics, for simply redundant measurement groups. Keyword `reportResidualCorrelation`. |
| `state_estimation.update_shunts` | Bool | `true` | `true`, `false` | Write the estimated susceptances of released shunts (`setShuntEstimation!`) back into the model (`Shunt.y_pu_shunt`/`B_shunt`) after a converged run. Off keeps the model data authoritative. Keyword `updateShunts`; estimates in `SEResult.shuntEstimates`. |
| `state_estimation.update_taps` | Bool | `false` | `true`, `false` | Write the fixed tap positions (`tap_ratio`/`phase_shift_deg`) of released transformers (`setTapEstimation!`) back after a converged fixed run; frozen regulators keep their model value. Off: an estimation never overwrites tap data silently. Keyword `updateTaps`; results in `SEResult.tapEstimates`/`tapFixation`. |
| `state_estimation.measurement_table` | Bool | `true` | `true`, `false` | Let the SE service read a [measurement table](state_estimation_measurement_table.md) (CSV, one measurement set per timestamp, keys resolved by CGMES mRID for a CGMES case and by element name for every other format). Off: such a file is refused with this key named. |
| `state_estimation.zero_injection` | Symbol/String | `on` | `on`, `off` | Zero-injection rows (`ZI_PINJ_bus_<k>`/`ZI_QINJ_bus_<k>` at `ZERO_INJECTION_SIGMA`) at passive buses (`findPassiveBuses`: no generation, load or shunt; not the slack, not isolated) that carry no injection measurement. `on` adds them in every SE path: the measurement generator, a measurement CSV v1, the measurements of a case file and a [measurement table](state_estimation_measurement_table.md) (service, Web UI, and `applyMeasurementSet!` before `runse!`). `off` adds none, and the generator writes zero balance rows instead. Rows a set already contains are read as they are. The Web UI generator checkbox "passive as zero-injection constraint" shows and sets this key for the case. The key is case scope: a case with its own `<case>.config.yaml` takes the value from that file (written by the generator checkbox), not from the general configuration file. |
| `state_estimation.sigma_by_voltage_level` | mapping | `{}` (empty) | one row per level keyed by the nominal voltage in kV, block style | Measurement sigmas per voltage level in physical units, for rows without their own sigma and for the measurement table generators; see [Measurement sigma by voltage level](@ref se-sigma-by-level). |
| `state_estimation.topology_precheck` | Bool | `true` | `true`, `false` | Run the advisory stage-1 topology checks (`validate_topology`, linear and cheap) at the start of every estimation; findings warn and land in `SEResult.topologyFindings`, nothing is blocked. |
| `state_estimation.topology_open_flow_k` | Float64 | `4.0` | positive | Sigma multiple above which a measured flow over an open element is a status contradiction. |
| `state_estimation.topology_dead_flow_k` | Float64 | `3.0` | positive | Sigma multiple below which a closed branch's flow measurements read dead; fires only when the model state expects a flow above the same multiple. |
| `state_estimation.topology_voltage_k` | Float64 | `4.0` | positive | Sigma multiple of the voltage disagreement across a closed link. |
| `state_estimation.topology_kcl_k` | Float64 | `4.0` | positive | Sigma multiple of the node-balance check at completely measured nodes. |
| `state_estimation.topology_cluster_min` | Int | `3` | positive integer | Minimum surviving suspects at one station for the stage-2 classification `:topology_error_suspected_at_station` of `runse_diagnostics`; very small nets cluster trivially. |
| `state_estimation.robust_mode` | Symbol | `off` | `off`, `staged`, `replacement` | Solve-weight family: `off` is plain WLS; `staged` the two-stage R modification with the knees `robust_k1`/`robust_k2`; `replacement` pins rows whose converged normalized residual reaches `k_suppress` to `suppression_sigma` (an EMS-style suppression list, bounded re-solve rounds on the converged state). Statistics keep the original sigmas in every mode. The smoothing modes serve online runs with sporadic gross errors; identification uses the sequential elimination. |
| `state_estimation.robust` | Bool | `false` | `true`, `false` | Deprecated alias for `robust_mode = staged`, read only while `robust_mode` is `off`. Set `robust_mode` instead, never both. |
| `state_estimation.robust_start_iteration` | Int | `3` | positive integer | First WLS iteration with robust weight modification; earlier iterations solve unmodified, because residuals say nothing before the state settles. |
| `state_estimation.robust_k1` | Float64 | `3.0` | positive, `<= robust_k2` | Lower knee of `staged`, on the raw ratio `t = \|r\|/sigma` (not the normalized residual `rn` of `k_eliminate` and `k_suppress`): weights are unchanged up to `t = k1`, the tangential band starts there. |
| `state_estimation.robust_k2` | Float64 | `6.0` | `>= robust_k1` | Upper knee of `staged`, same raw ratio; beyond it the gradient contribution is suppressed (`sigma_mod = \|r\|/k1`). The defaults reproduce the classic 3/6 behavior bitwise. |
| `state_estimation.k_eliminate` | Float64 | `3.0` | positive | Elimination limit on the normalized residual `rn = r/sqrt(Omega_ii)` (keyword `normalizedThreshold`): from it a row is removed from the estimate, not just down-weighted. Rarely set below 3. The Web UI warns when `k_suppress < k_eliminate`. |
| `state_estimation.k_suppress` | Float64 | `4.0` | positive | Down-weight limit of `replacement`, on the same `rn` as `k_eliminate`, so the two are directly comparable: a row whose converged `rn` reaches it is solved with `suppression_sigma`. The row stays in every statistic. |
| `state_estimation.suppression_sigma` | Float64 | `2000.0` | positive | The sigma a suppressed row is solved with, in the measurement's unit; statistics keep the original sigma. It removes the row's pull without deleting it; values near the real sigma call for `staged` instead. |
| `state_estimation.max_eliminations` | Int | `3` | `>= 0` | Budget of the sequential elimination: at most this many rows per diagnostics run, a stop condition rather than a limit; `0` disables the elimination and leaves the diagnostics reporting only. Each increment costs one further estimation run. Past three removals the cause is usually a model or topology error, which stage 2 reports (stop reasons `:max_eliminations`, `:no_localizable_suspect`). |
| `state_estimation.rank_tol_factor` | Float64 | `10.0` | positive real | Observability rank tolerance `factor * jacEps * sigma_max` of the normalized Jacobian, matching the error floor of the finite-difference Jacobian. Lower it only for exact analytic Jacobians passed with an explicit `tol`. |
| `state_estimation.rank_method` | Symbol/String | `decomposition` | `decomposition`, `pivots` | Source of the observability rank: `decomposition` is the SVD of the Jacobian (sparse QR above 2000 states), exact with the deficit; `pivots` reads it from the LDLt factorization of the gain matrix that the criticality pass computes anyway, exact for full rank, and hands any deficit to the decomposition. The report names the method used in `rank_method`; a disagreement between the two is a finding, not a tolerance to tune. Use `pivots` on large networks once it has agreed with the decomposition on your cases. |
| `state_estimation.criticality_method` | Symbol/String | `omega` | `omega`, `rank` | How critical measurements (structurally zero residual) are found: `omega` reads `Omega_ii` for every row from one factorization and one selected-inverse pass, no size budget, and reports `wii` as the graded nearly-critical information; `rank` runs one rank test per row (budgeted at 300000 rows times states, `criticality_skipped` above) as a cross-check. |
| `state_estimation.takahashi_min_states` | Int | `200` | positive integer | State count from which the diagnostics compute `Omega_ii` with the sparse Takahashi selected inverse instead of the dense `pinv`; diagnostics only, `omega_path` in the result. K-matrix requests and guard failures use the dense path automatically. |
| `state_estimation.linear_solver` | Symbol/String | `umfpack` | `umfpack`, `klu` | Sparse LU of the normal equations `G dx = g` when the gain matrix is not exactly symmetric; an exactly symmetric one takes the Cholesky factorization on either setting. Both keep the symbolic analysis per sparsity pattern. `klu` needs the KLU package extension (`using KLU` next to Sparlectra) and falls back to `umfpack` with one warning without it; the estimates agree to rounding, not bitwise. `umfpack` stays the default when KLU is loaded for power mode (the power-flow counterpart is `power_flow.linear_solver`). |
| `state_estimation.symmetric_gain` | Bool | `false` | `true`, `false` | Build `G = H'WH` exactly symmetric (upper triangle mirrored), so the solve always takes the Cholesky route, several times faster per iteration on large networks; one extra copy of `G` per iteration. The estimates differ from the plain product at the 1e-10 pu level, so leave it off where results must match earlier runs. |
| `state_estimation.observability.enabled` | Bool | `true` | `true`, `false` | Run the observability checks. |

## [Measurement sigma by voltage level](@id se-sigma-by-level)

`measurementStdDevs` gives one sigma per measurement type for the whole
network. Instrument transformers and metering differ between voltage
levels, so `state_estimation.sigma_by_voltage_level` sets the sigmas per
level and measurement type in physical units. One row per level, keyed
by the nominal voltage of the level in kV as in the network (bus `Vn`,
matched to 1e-6 kV). The columns of a row:

| column | unit | meaning |
|---|---|---|
| `p_mw` | MW | sigma of P measurements (flows and injections) |
| `q_mvar` | MVAr | sigma of Q measurements (flows, injections, shunt bays) |
| `i_a` | A | sigma of current magnitudes |
| `u_kv` | kV | sigma of voltage magnitudes |
| `rel_pct` | % | optional, default 0: the sigma is `max(absolute, rel_pct/100 * \|value\|)` |

A missing column keeps the global default for that type, and voltage
angles (`VaMeas`, `IaMeas`) always do. The key must be a positive number
and name each level once (`380` and `0380` are the same level), every
given sigma positive, `rel_pct` not negative, and a row sets at least one
absolute sigma (`rel_pct` alone changes nothing). A bad row is an
`ArgumentError` that names the key, the row (its key in the file, its
position in a vector) and a reason: `bad_vn_kv`, `duplicate_vn_kv`,
`bad_value`, `no_sigma`, `unknown_column` (also for a `vn_kv` column,
the nominal voltage is the key), `not_a_row`, `flow_style`. For a table in a configuration file or a case configuration
file the message ends with the file name and the line of the offending
column or row, for example `(sp_case14.config.yaml line 9)`.

Order per measurement:

1. The sigma of the row itself (the `sigma` column of a measurement table,
   every row of a measurement CSV v1, a measurement set in code).
2. The row of the level of the measuring bus: the bus of the node, of the
   injection or shunt bay, or, at a branch end, the bus at that end. A
   measurement at the 110 kV side of a 380/110 kV transformer takes the
   110 kV row, one at the 380 kV side the 380 kV row.
3. The default (`measurementStdDevs`).

In practice the table acts on measurement tables with empty `sigma` cells
and on the generators: `generateMeasurementTableFromPF`,
`generateMeasurementTableFromSV`, `generatePseudoMeasurementTable` and the
measurement generator of the Web UI (a measurement table) write the level sigma
into each row and draw the noise with it, where the level of the row has a
value for the type; elsewhere their own sigma settings apply. The Web UI
generator reads the table of the case configuration and names the number
of such rows in the comments of the set. A pseudo row then carries the sigma of a meter at that level; pass
`sigma_by_voltage_level = Sparlectra.SESigmaLevel[]` to keep the larger pseudo
sigma. Each of the three takes the table as a keyword, default the active
configuration, like `importMeasurementTable`.

The estimator keeps P and Q in MW and MVAr, I in A and U in pu. The sigmas
enter at the import of a row, in the unit of the estimator:

```math
\sigma_{P,\mathrm{pu}} = \frac{\sigma_{P,\mathrm{MW}}}{S_\mathrm{base}}, \quad
\sigma_{Q,\mathrm{pu}} = \frac{\sigma_{Q,\mathrm{MVAr}}}{S_\mathrm{base}}, \quad
\sigma_{U,\mathrm{pu}} = \frac{\sigma_{U,\mathrm{kV}}}{V_n}, \quad
\sigma_{I,\mathrm{pu}} = \frac{\sigma_{I,\mathrm{A}}}{I_\mathrm{base}}, \quad
I_\mathrm{base} = \frac{1000\, S_\mathrm{base}}{\sqrt{3}\, V_n}\ \mathrm{A}
```

with $S_\mathrm{base}$ = `baseMVA` and $V_n$ the nominal voltage of the
measuring bus in kV. P, Q and I rows keep their sigma in MW, MVAr and A,
which gives the same normalized residual $(z - h)/\sigma$ as the pu form;
U rows get $\sigma_{U,\mathrm{pu}}$. A sigma below `ZERO_INJECTION_SIGMA`
in the unit of the estimator is raised to it and counted as
`sigma_lifted`, as for a sigma in the file. Example: `p_mw = 2` and
`rel_pct = 1` give 3 MW at a reading of 300 MW and 2 MW at 100 MW.

The table is written in block style, one row per level under its nominal
voltage, one column per line:

```yaml
state_estimation:
  sigma_by_voltage_level:
    380:
      p_mw: 2.0
      q_mvar: 2.0
      i_a: 20
      u_kv: 0.5
      rel_pct: 0.5
    220:
      p_mw: 1.5
      q_mvar: 1.5
      i_a: 25
      u_kv: 0.4
      rel_pct: 0.8
    110:
      p_mw: 1.0
      q_mvar: 1.0
      i_a: 30
      u_kv: 0.3
      rel_pct: 1.0
```

The numbers are placeholders; take them from your metering. The key of a
row is the nominal voltage in kV, with a decimal point where needed
(`10.5:`, `0.4:`); a key in quotes (`"10.5":`) reads the same. In code
the key may be any positive number or its text. A key written twice is
refused with `duplicate_vn_kv` and both lines, since the reader would keep
only the second row. A row written as a flow mapping (`380: {p_mw: 2.0}`)
is refused with the error `flow_style` and this block form, never read as
text; a line that starts with `- ` or `{` stops the YAML reader with the
same hint.

The key is case scope: it usually lives in the case configuration file
next to the case (`<case>.config.yaml`), which wins over the general file.
The State Estimation section of the Web UI shows the table with one row
per nominal voltage of the case (a configured level the case does not have
is shown and marked); an empty cell keeps the default. The cells show the
table of the case configuration and are empty without one, also when the
general file sets the key. An empty table that was empty on the page posts
nothing, so such a run resolves the key as an API run does; clearing the
cells of a stored table clears it. **Save settings for this case** writes
the table into the case configuration, keyed by the nominal voltage, and a
run takes the table of the form. The table is also shown before the case has a measurement set, with
its own save button, since the generator reads it. A stored table that cannot be read is shown with its error and
disabled cells; the form then posts no table, so a run stops with the same
error and a save keeps the stored value.

Each SE run reports, per nominal voltage of the case, the number of
measurements and how many of them carry their own sigma, the level sigma or
the default, and the configured sigmas in pu, in `run.log`,
`se_diagnostics.md`, on the result page and in `result.json`
(`metadata.se_sigma_by_voltage_level`). A configured level that no bus of
the case has is a note in the same places
(`metadata.se_sigma_levels_not_in_case`), not an error. Zero-injection rows
are constraints and are not counted.

# [State-Estimation Configuration](@id se-config)

The estimator reads every setting from the active configuration:
`runse!(net)`, `runse_diagnostics`, `validate_measurements` and
`validate_topology` take no settings.
`with_state_estimation_config(max_iter = 12) do runse!(net) end` installs
other values for one call and restores the previous configuration
afterwards. The service and the Web UI install the effective configuration
of the case (general file, case file, form values) the same way.

`flatstart`, `robust_mode`, `robust`, `k_eliminate`, `k_suppress`,
`max_eliminations`, `topology_precheck` and `report_residual_correlation`
are case-scope keys like `power_flow.solver`: a case configuration file may
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

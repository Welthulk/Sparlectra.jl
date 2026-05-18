# State-Estimation Configuration

| YAML path | Type | Default | Allowed values | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `state_estimation.enabled` | Bool | `true` | `true`, `false` | Enables SE stage. | Measurement-driven workflows. | No measurement model/data. | Extra nonlinear solve pass. | `update_net`, PF start source. |
| `state_estimation.method` | Symbol/String | `wls` | `wls` | State-estimation algorithm. | Standard weighted least squares. | Unsupported methods. | Fixed current implementation. | Validation rejects unsupported values. |
| `state_estimation.sparse` | Bool | `true` | `true`, `false` | Sparse linear algebra in SE. | Medium/large systems. | Rare tiny dense debug runs. | Better memory/runtime scaling. | Runtime thread/BLAS choices. |
| `state_estimation.tol` | Float64 | `1e-8` | positive real | SE convergence tolerance. | High-accuracy estimation. | Excessively tight noisy runs. | Tighter means more iterations. | `max_iter`. |
| `state_estimation.max_iter` | Int | `20` | positive integer | Iteration cap for SE solver. | Hard estimation scenarios. | Too low for challenging observability. | Runtime upper bound. | `tol`, observability quality. |
| `state_estimation.flatstart` | Bool | `true` | `true`, `false` | Start estimator from synthetic flat state. | No trusted historical/SCADA initial state. | When measured/historical start exists. | Low. | `update_net` and PF start chain. |
| `state_estimation.jac_eps` | Float64 | `1e-6` | positive real | Jacobian perturbation epsilon (numeric terms). | Numeric stability tuning. | Extreme epsilon values. | Low/medium sensitivity impact. | `tol` and convergence behavior. |
| `state_estimation.update_net` | Bool | `true` | `true`, `false` | Write estimated state back to network. | PF should consume SE result as start. | Isolated SE experiments only. | Low. | Links SE output to PF input. |
| `state_estimation.observability.enabled` | Bool | `true` | `true`, `false` | Run observability checks. | SCADA/SE validation workflows. | Fast-only runs where observability not needed. | Small extra analysis cost. | Diagnostics/output verbosity. |

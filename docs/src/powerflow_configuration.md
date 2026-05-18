# Power-Flow Configuration

## Start-source concepts

| Start source | Meaning | Requires previous solution? | Typical use | Risk |
|---|---|---:|---|---|
| Flat start | Synthetic `1.0 pu / 0°` style initialization | no | Small/simple cases | Can converge poorly on stressed grids |
| DC angle start | Synthetic DC-like angle estimate | no | Large transmission-like grids | Ignores full AC effects |
| MATPOWER imported reference | `BUS.VM/VA`, `GEN.VG`, imported setpoint logic | maybe | Validation against imported cases | Can resemble using known answer |
| Historical values | Previous trusted operating point | yes | Repeated operations | Stale operating point |
| SCADA/state-estimation values | Measured/estimated state | yes | Operations workflows | Measurement noise/outliers |
| Profile values | Explicit configured profile values | maybe | Reproducible studies | May not match topology changes |

## Solver core options

| YAML path | Type | Default | Allowed values | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `power_flow.method` | Symbol/String | `rectangular` | `rectangular` | AC solver formulation. | Always (current core). | N/A | Fixed single implementation. | Must match `benchmark.methods`. |
| `power_flow.sparse` | Bool | `true` | `true` | Sparse linear algebra mode. | Always (required). | N/A | Scales better for large systems. | Validated together with `method`. |
| `power_flow.flatstart` | Bool | `true` | `true`, `false` | Legacy flatstart toggle. | Synthetic starts. | When start projection and imported starts are used. | Low. | Combined into `start_mode`. |
| `power_flow.tol` | Float64 | `1.0e-5` | positive real | PF tolerance. | Accuracy-sensitive studies. | Overly tight on big batches. | Tighter means more iterations. | `max_iter`. |
| `power_flow.max_iter` | Int | `80` | positive integer | Iteration cap. | Hard cases. | Very low values. | Upper runtime bound. | `tol`, `qlimits`. |
| `power_flow.autodamp` | Bool | `true` | `true`, `false` | Adaptive damping. | Difficult convergence. | Strict algorithm comparison. | Small overhead, often fewer failures. | `autodamp_min`. |
| `power_flow.autodamp_min` | Float64 | `0.05` | positive real | Minimum damping factor. | Stabilizing hard cases. | Near-zero damping on easy grids. | Lower can increase iterations. | Active only with `autodamp=true`. |

## Start mode options

| YAML path | Type | Default | Allowed values | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `power_flow.start_mode.angle_mode` | Symbol/String | `dc` | `classic`, `dc`, `bus_va_blend`, `matpower_va` | Angle initialization source mode. | Transmission-like starts. | If trusted measured/historical state exists. | Can reduce iterations. | `try_dc_start`, `dc_angle_limit_deg`. |
| `power_flow.start_mode.voltage_mode` | Symbol/String | `bus_vm_va_blend` | `classic`, `pv_gen_vg`, `pv_bus_vm`, `all_bus_vm`, `bus_vm_va_blend` | Voltage-magnitude/angle blend strategy. | Imported-reference assisted starts. | Untrusted imported data. | Small startup overhead. | `blend_lambdas`, `reuse_import_data`. |
| `power_flow.start_mode.start_projection` | Bool | `true` | `true`, `false` | Enables start-projection workflow. | Robustness on hard cases. | Minimal-path microbench runs. | Extra startup pass. | Gates start-projection sub-options. |
| `power_flow.start_mode.try_dc_start` | Bool | `true` | `true`, `false` | Try DC candidate start. | Large transmission cases. | Highly resistive distribution cases. | Low overhead. | `dc_angle_limit_deg`. |
| `power_flow.start_mode.try_blend_scan` | Bool | `true` | `true`, `false` | Scan blend candidates. | Mixed/tricky start data. | Easy cases needing speed. | Startup cost ∝ lambda count. | `blend_lambdas`. |
| `power_flow.start_mode.branch_guard` | Bool | `true` | `true`, `false` | Branch sanity guard for candidate starts. | Stability-focused runs. | Rarely disabled. | Low. | Candidate measurement options. |
| `power_flow.start_mode.measure_candidates` | Bool | `true` | `true`, `false` | Score/select among candidates. | Multiple start candidates. | Fastest startup path only. | Low/medium startup overhead. | `try_dc_start`, `try_blend_scan`. |
| `power_flow.start_mode.accept_unmeasured_dc_start` | Bool | `false` | `true`, `false` | Allow DC start without measurement checks. | Synthetic studies. | Measurement-driven workflows. | Can avoid fallback retries. | `try_dc_start`. |
| `power_flow.start_mode.reuse_import_data` | Bool | `true` | `true`, `false` | Reuse imported MATPOWER references. | Trusted imports. | Uncertain conversion data. | Small reduction in recomputation. | `matpower_import.*` voltage reference keys. |
| `power_flow.start_mode.blend_lambdas` | Vector{Float64} | `[0.25,0.5,0.75]` | real vector (typ. 0..1) | Lambda candidates for blend scan. | Need robust candidate search. | Very large lambda sets. | Linear startup growth with vector size. | `try_blend_scan`. |
| `power_flow.start_mode.dc_angle_limit_deg` | Float64 | `60.0` | positive real | DC-start angle magnitude cap (deg). | Conservative angle starts. | Overly restrictive values. | Negligible. | `try_dc_start`. |

## Q-limit options and guard

| YAML path | Type | Default | Allowed values | Meaning | Use when | Avoid when | Performance impact | Interactions |
|---|---:|---:|---|---|---|---|---|---|
| `power_flow.qlimits.enabled` | Bool | `true` | `true`, `false` | Master switch for Q-limit enforcement. | Realistic PV/PQ behavior needed. | Pure unconstrained PF tests. | Can increase switching iterations. | Gates all `qlimits.*`. |
| `power_flow.qlimits.start_iter` | Int | `3` | integer | First iteration index for Q limits. | Delay switching noise early. | Very late switching on hard cases. | Affects convergence speed. | `start_mode`, `auto_q_delta_pu`. |
| `power_flow.qlimits.start_mode` | Symbol/String | `iteration_or_auto` | `iteration`, `auto`, `iteration_or_auto` | Activation policy. | Mixed robustness/perf runs. | Mismatched with expected policy. | Small control logic overhead. | `auto_q_delta_pu`. |
| `power_flow.qlimits.auto_q_delta_pu` | Float64 | `1e-4` | nonnegative real | Auto activation threshold. | Fine-tuning switch timing. | Extreme values. | Low. | `start_mode=auto` or `iteration_or_auto`. |
| `power_flow.qlimits.hysteresis_pu` | Float64 | `0.01` | nonnegative real | Hysteresis margin near Q limits. | Reduce switch chattering. | Too large if strict tracking needed. | Can reduce oscillatory iterations. | `cooldown_iters`, guard modes. |
| `power_flow.qlimits.cooldown_iters` | Int | `1` | nonnegative integer | Cooldown iterations after switching. | Reduce repeated toggling. | Too long cooldown on tight limits. | Affects convergence pace. | Hysteresis and freeze behavior. |
| `power_flow.qlimits.trace_buses` | Vector{Int} | `[]` | bus-id vector | Trace selected bus events. | Targeted diagnostics. | Large full-network trace. | Logging overhead if populated. | Output and diagnostics verbosity. |
| `power_flow.qlimits.lock_pv_to_pq_buses` | Vector{Int} | `[]` | bus-id vector | Force listed buses into PQ-lock behavior. | Known problematic buses. | Blindly on all buses. | Can simplify switching dynamics. | Guard modes. |
| `power_flow.qlimits.guard.enabled` | Bool | `true` | `true`, `false` | Enable guard subsystem. | Prevent unstable switching. | Pure baseline comparisons. | Small runtime overhead. | Guard fields below. |
| `power_flow.qlimits.guard.min_q_range_pu` | Float64 | `0.02` | nonnegative real | Range threshold for narrow/zero detection. | Robust zero/narrow range handling. | Too high threshold. | Low. | Narrow/zero modes. |
| `power_flow.qlimits.guard.narrow_range_mode` | Symbol/String | `lock_pq` | `prefer_pq`, `lock_pq` | Action for narrow Q range units. | Convergence protection. | If strict PV control required. | Can reduce oscillations. | Hysteresis/cooldown. |
| `power_flow.qlimits.guard.zero_range_mode` | Symbol/String | `lock_pq` | `lock_pq` | Action for zero Q range units. | Deterministic limit handling. | N/A | Low. | Lock lists and violation mode. |
| `power_flow.qlimits.guard.violation_mode` | Symbol/String | `lock_pq` | `delayed_switch`, `lock_pq` | Action on persistent violations. | Robustness under bad limits. | Aggressive switching studies. | Can add control logic. | Threshold and switch caps. |
| `power_flow.qlimits.guard.violation_threshold_pu` | Float64 | `1e-4` | nonnegative real | Violation threshold. | Tune sensitivity. | Extreme values. | Low. | `violation_mode`. |
| `power_flow.qlimits.guard.max_switches` | Int | `3` | nonnegative integer | Max switches before freeze logic. | Stop chattering. | Too low on valid dynamic cases. | Can reduce wasted iterations. | `freeze_after_repeated_switching`. |
| `power_flow.qlimits.guard.max_remaining_violations` | Int | `0` | nonnegative integer | Allowed violations at guarded exit. | Controlled tolerance policies. | Strict zero-violation policies. | Low. | `accept_bounded_violations`. |
| `power_flow.qlimits.guard.accept_bounded_violations` | Bool | `false` | `true`, `false` | Permit bounded residual violations. | Practical operations tradeoff. | Strict compliance studies. | May reduce retries. | `max_remaining_violations`. |
| `power_flow.qlimits.guard.freeze_after_repeated_switching` | Bool | `true` | `true`, `false` | Freeze after repeated switch cycling. | Anti-chatter behavior. | Cases requiring unrestricted switching. | Can stabilize solves. | `max_switches`. |
| `power_flow.qlimits.guard.log` | Bool | `true` | `true`, `false` | Emit guard logs. | Diagnostics/debugging. | Quiet batch runs. | I/O overhead when enabled. | `output.console_q_limit_events`. |

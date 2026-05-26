# Generic Control Framework

## Purpose

The inner numerical solver remains `runpf!`. The generic orchestration layer is `run_control!`, which runs controller-driven outer iterations around repeated PF solves.

`run_acpflow` automatically dispatches through `run_control!` when `collect_outer_controllers(net)` returns at least one controller.

## Architecture

```text
run_acpflow (public entry)
        |
        v
collect_outer_controllers(net)
        |
        v
run_control!
        |
        +--> runpf!
        +--> control_evaluate!
        +--> control_propose_update!
        +--> control_apply_update!
        +--> runpf! again
        |
        v
ControlRunResult stored on net.control_result
```

Preferred public entries:

```julia
run_acpflow(; net = net, ...)
run_acpflow(; casefile = "case14.m", path = "...", ...)
```

`run_acpflow` is the public high-level ACP runner for both in-memory and file-based workflows.

## Hook interface

- `AbstractOuterController`
- `AbstractControlState`
- `AbstractControlUpdate`
- `control_initialize!`
- `control_evaluate!`
- `control_propose_update!`
- `control_apply_update!`
- `control_is_converged`
- `control_is_blocked`
- `control_status`
- `control_report_rows`
- `control_trace_rows`
- `control_max_outer_iterations`

`control_max_outer_iterations` provides controller-specific outer-loop limits. The global outer-loop budget is `control.max_outer_iterations` and is combined with controller limits.

`control_max_outer_iterations` is currently treated as an internal extension hook (not exported). External custom controllers can still extend it via `Sparlectra.control_max_outer_iterations(::MyController)`.

## Result model

`ControlRunResult` contains:

- `status`
- `converged`
- `outer_iterations`
- `powerflow_solves`
- `last_pf_iterations`
- `last_pf_status`
- `controllers`
- `trace`

Terminal statuses:

- `:no_controllers`
- `:disabled`
- `:no_active_controllers`
- `:pf_failed`
- `:converged`
- `:blocked`
- `:max_outer_iterations`

## Legacy status boundary (`erg`)

At the legacy API boundary, `erg` reflects inner numerical PF success/failure only.

- `:pf_failed` maps to failure (`erg = 1`).
- Control-loop outcomes such as `:blocked` or `:max_outer_iterations` are not inner numerical PF failures.

Inspect control-loop outcome via `latest_control_result(net)` or `net.control_result`.

## Latest-result access

Use:

```julia
latest_control_result(net)
net.control_result
```

These expose the latest control run associated with the `Net` instance.

## Trace rows (transformer control)

Minimal row schema:

- `outer_iteration`
- `controller_name`
- `controller_type`
- `transformer_id`
- `mode`
- `status`
- `converged`
- `at_limit`
- `achieved_vm_pu`
- `target_vm_pu`
- `achieved_p_mw`
- `target_p_mw`
- `tap_ratio`
- `phase_shift_deg`

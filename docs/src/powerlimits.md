# Reactive Power Limits (Q-Limits) and PV→PQ Switching

A PV bus that violates its reactive limits becomes a PQ bus during the
Newton-Raphson iteration (PV→PQ switching); events are logged, and
hysteresis and a cooldown period allow switching back.

## Data Structures in `Net`

The Q-limit state lives on the network:

| field | meaning |
|---|---|
| `qmin_pu`, `qmax_pu` | per-bus limits in p.u., aggregated from the generators' `qMin`/`qMax` by `buildQLimits!` (or set with `setQLimits!`) before the power flow |
| `q_hyst_pu`, `cooldown_iters` | hysteresis band in p.u. and minimum iterations between PV/PQ state changes |
| `reenable_v_hyst_pu`, `final_q_accept_pu` | voltage margin of the PQ→PV release and size bound of the final check |
| `qLimitLog`, `qLimitEvents` | chronological `QLimitEvent` records and the last side (`:min` or `:max`) per bus |

## Pre-run Q-limit Preview in MVAr

`verbose > 0` prints the PV limits in MVAr (via `baseMVA`) before the
Newton loop, all rows with `verbose > 1`: a check that catches pu/MVAr
mixups early.

## Logging Q-Limit Hits

```julia
Base.@kwdef struct QLimitEvent
  iter::Int
  bus::Int
  side::Symbol         # :min | :max
  kind::Symbol = :hit  # :group_revert when a jointly released group returns to its limits
end
```

* `logQLimitHit!(net, iter, bus, side)`: appends a `QLimitEvent` to
  `net.qLimitLog` and sets `net.qLimitEvents[bus] = side`.
* `lastQLimitIter(net, bus)`: last iteration in which `bus` hit a limit.
* `resetQLimitLog!(net)`: clears both.
* `pv_hit_q_limit(net, ["B3", "B10"])`: `true` if any listed bus appears
  in `net.qLimitEvents` (for tests).
* `printQLimitLog(net; sort_by = :iter, io = stdout)`: the events as a
  table.

## Conceptual Algorithm

In every Newton-Raphson iteration of the rectangular solver:

1. **Solve one NR step**; obtain the voltages `V` and bus powers
   `S = P + jQ`.
2. **Compute the generator Q** per PV bus (net injection plus bus load).
3. **Check the limits**: bus `b` has hit a limit when its generator Q leaves
   `[Qmin_b, Qmax_b]` (`net.qmin_pu`, `net.qmax_pu`) by more than the
   hysteresis margin `q_hyst` (`power_flow.qlimits.hysteresis_pu`, default
   `0.01`; `0` is the plain limit test),

   ```math
   Q_b < Qmin_b - q_\mathrm{hyst} \quad \text{or} \quad Q_b > Qmax_b + q_\mathrm{hyst}
   ```

4. **Clip Q and switch the bus type**: set the generator's Q to the
   violated limit (`Q_b := Qmin_b` or `Q_b := Qmax_b`), then

   ```julia
   setBusType!(net, b, "PQ")
   logQLimitHit!(net, iter, b, :min)  # or :max
   ```

   From the next iteration on the solver's `bus_types` vector carries `:PQ`
   for the bus and its residual switches from `(ΔP, ΔV)` to `(ΔP, ΔQ)`
   ([Solver Guide](solver.md)).

5. **PV→PQ→PV re-enable**, only with `net.q_hyst_pu > 0` or
   `net.cooldown_iters > 0` (both at zero disables switching back). In the
   active-set step a bus clamped at Qmax is released when
   `Vm > Vset + reenable_v_hyst_pu`, one clamped at Qmin when
   `Vm < Vset - reenable_v_hyst_pu` (`power_flow.qlimits.reenable_v_hyst_pu`,
   default `1e-4` pu): a machine at Qmax whose voltage sits above its
   setpoint needs less than Qmax to hold it (the back-off rule of Sundaresh
   and Rao, Sadhana 40(4), 2015). Callers that pass no voltages fall back to
   the Q band test

   ```math
   Qmin_b + q_\mathrm{hyst} < Q_b < Qmax_b - q_\mathrm{hyst}
   ```

   In both cases enough iterations must have passed since the last hit:

   ```math
   \text{iter} - \text{lastQLimitIter}(b) \ge \text{cooldown\_iters}
   ```

   Cooldown and the one-retry guard apply: a machine flipping between the
   clamp and the voltage constraint is held after its first retry. At the
   converged point the machines held this way whose voltage is on the
   release side get one joint release; if one of them hits its limit again,
   the whole group goes back to its limits and stays there.
   `examples/powerflow/example_qlimit_reenable_voltage_rule.jl` shows the
   difference on the Zeng/Chiang 14-bus case.

6. **Repeat** with the updated bus types and Q values.

7. **Final check, the same for every enforcement mode.** After the solve
   the overshoot `d = Q - Qmax` (or `Qmin - Q`, in pu) of every PV bus is
   judged by its size, not by the number of buses:

   | class | condition | status | run |
   |---|---|---|---|
   | in order | `d <= hysteresis_pu` | `ok` (`within_hysteresis` when `d > tol`) | accepted |
   | bounded | `hysteresis_pu < d <= final_q_accept_pu` | `bounded_q_limit_violation` | accepted, warning with bus, `d` in pu and MVAr |
   | violation | `d > final_q_accept_pu` | `remaining_pv_q_limit_violations` | not accepted |

   `power_flow.qlimits.final_q_accept_pu` (default `auto`, twice
   `hysteresis_pu`; an explicit value must be `>= hysteresis_pu`) bounds
   what is still accepted; both at zero make the check strict. The
   count-based guard keys (`guard.accept_bounded_violations`,
   `guard.max_remaining_violations`) are not part of this classification.
   The check writes one line per bus (side, Q, limit, `d` in pu and MVAr,
   class) into the Q-limit block of the text report, `run.log`, the run
   metadata (`final_q_check_status`, with `qlimits_disabled` when Q-limit
   handling is off and `not_evaluated` without a converged solution,
   `final_q_check_buses`, `final_q_check_max_dev_pu`) and the Web UI result
   page, next to the voltage-side Q-V characteristic check.

## Alternative to immediate PV→PQ switching: `qlimit_mode = :adjust_vset`

```julia
runpf!(net, 40, 1e-8, 1; qlimit_mode = :adjust_vset)
```

When a PV bus reaches a limit, its voltage setpoint (`vm_pu`) is first
moved in steps of `vstep_pu` before the bus becomes PQ, which avoids a
switch where a small voltage correction keeps Q inside the limits. The
regulating generator prosumer at the bus (one controller definition per
bus) provides `vstep_pu` and optionally `tap_steps_down` / `tap_steps_up`
(maximum number of adjustments). Without a valid controller, or once the
step limits are exhausted, standard PV→PQ switching applies. Rectangular
solver only; the default is `:switch_to_pq`.

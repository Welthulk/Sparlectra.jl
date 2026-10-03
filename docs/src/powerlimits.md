# Reactive Power Limits (Q-Limits) and PV→PQ Switching

Sparlectra enforces per-bus reactive power limits for generators (PV
buses). A PV bus that violates its Q-limits becomes a PQ bus during the
Newton-Raphson iteration (PV→PQ switching). Events are logged; hysteresis
and a cooldown period optionally allow switching back.

## Data Structures in `Net`

```julia
struct Net
    name::String
    baseMVA::Float64
    slackVec::Vector{Int}
    vmin_pu::Float64
    vmax_pu::Float64
    nodeVec::Vector{Node}
    linesAC::Vector{ACLineSegment}
    trafos::Vector{PowerTransformer}
    branchVec::Vector{Branch}
    prosumpsVec::Vector{ProSumer}
    shuntVec::Vector{Shunt}
    busDict::Dict{String,Int}
    busOrigIdxDict::Dict{Int,Int}
    totalLosses::Vector{Tuple{Float64,Float64}}
    totalBusPower::Vector{Tuple{Float64,Float64}}
    _locked::Bool
    shuntDict::Dict{Int,Int}
    isoNodes::Vector{Int}

    qLimitLog::Vector{Any}          # chronological log of Q-limit events
    cooldown_iters::Int             # min. iterations between PV/PQ state changes
    q_hyst_pu::Float64              # hysteresis band in p.u. for Q-backtracking
    qmin_pu::Vector{Float64}        # per-bus Qmin (p.u.)
    qmax_pu::Vector{Float64}        # per-bus Qmax (p.u.)
    qLimitEvents::Dict{Int,Symbol}  # BusIdx -> :min | :max (last PV→PQ change)
end
```

The limits come from generators / prosumers (`qMin`, `qMax`), converted to
p.u. and aggregated per bus; `setQLimits!` / `buildQLimits!` pre-fill
`net.qmin_pu` and `net.qmax_pu` before the power flow.

## Pre-run Q-limit Preview in MVAr

With `verbose > 1` a Q-limit preview in MVAr at bus level (from the
per-unit limits via `baseMVA`) is printed before the Newton-Raphson loop,
a plausibility check that catches unit mixups (pu vs. MVAr) early.

## Logging Q-Limit Hits

Q-limit events are `QLimitEvent` records in `net.qLimitLog`:

```julia
Base.@kwdef struct QLimitEvent
    iter::Int
    bus::Int
    side::Symbol   # :min | :max
end
```

* `logQLimitHit!(net, iter, bus, side)`: appends a `QLimitEvent` to
  `net.qLimitLog` and sets `net.qLimitEvents[bus] = side`.
* `lastQLimitIter(net, bus)`: last iteration in which `bus` hit a limit.
* `resetQLimitLog!(net)`: clears `qLimitLog` and `qLimitEvents`.
* `pv_hit_q_limit(net, ["B3", "B10"])`: `true` if any listed PV bus
  appears in `net.qLimitEvents` (for tests).

```julia
printQLimitLog(net; sort_by = :iter, io = stdout)
```

```text
──────────────────────────────────
 Iteration │ Bus │ Side
──────────────────────────────────
         3 │   5 │ min
         5 │   3 │ max
──────────────────────────────────
Total events: 2
```

---

## Conceptual Algorithm

During a Newton-Raphson calculation (`calcNewtonRaphson_withPVIdentity!` or
the rectangular variant):

1. **Solve one NR step**; obtain bus voltages `V` and bus powers
   `S = P + jQ`.
2. **Compute actual generator Q** per PV bus.
3. **Check limits** per bus index `b` with `Q_b` the reactive injection
   and `[Qmin_b, Qmax_b]` from `net.qmin_pu` and `net.qmax_pu`; the bus has
   hit a limit if

   ```math
   Q_b < Qmin_b \quad \text{or} \quad Q_b > Qmax_b
   ```

4. **Clip Q and switch bus type**: set the generator's `Q` to the violated
   limit (`Q_b := Qmin_b` or `Q_b := Qmax_b`), then

   ```julia
   setBusType!(net, b, "PQ")
   logQLimitHit!(net, iter, b, :min)  # or :max
   ```

5. **Optional: PV→PQ→PV re-enable (hysteresis + cooldown).** If
   `net.q_hyst_pu > 0` and `net.cooldown_iters > 0`, a bus may switch back
   from PQ to PV only if its Q is inside the interval with margin,

   ```math
   Qmin_b + q_\mathrm{hyst} < Q_b < Qmax_b - q_\mathrm{hyst}
   ```

   and enough Newton iterations have passed since the last limit hit:

   ```math
   \text{iter} - \text{lastQLimitIter}(b) \ge \text{cooldown\_iters}
   ```

!!! note
    In many applications, switching back within the same power flow
    calculation is intentionally disabled (`q_hyst_pu = 0.0`,
    `cooldown_iters = 0`).

   **Voltage-side release (active set).** When the solver hands voltages
   and setpoints to the active-set step, a bus clamped at Qmax is released
   when `Vm > Vset + reenable_v_hyst_pu`, one clamped at Qmin when
   `Vm < Vset - reenable_v_hyst_pu` (`power_flow.qlimits.reenable_v_hyst_pu`,
   default `1e-4` pu): a machine at Qmax whose voltage sits above its
   setpoint needs less than Qmax to hold it (the "back off" rule of
   Sundaresh and Rao, Sadhana 40(4), 2015). The Q band test above is the
   fallback for callers that pass no voltages; cooldown and the one-retry
   guard still apply.
   `examples/powerflow/example_qlimit_reenable_voltage_rule.jl` shows the
   difference on the Zeng/Chiang 14-bus case.

6. **Repeat** with the updated bus types and Q-values in the next NR
   iteration.

7. **Final check, the same for every enforcement mode.** After the solve,
   every PV bus is checked against its reactive limits; the overshoot
   `d = Q - Qmax` (or `Qmin - Q`, in pu) is judged by its size, not by the
   number of buses:

   | class | condition | status | run |
   |---|---|---|---|
   | in order | `d <= hysteresis_pu` | `ok` (`within_hysteresis` when `d > tol`) | accepted |
   | bounded | `hysteresis_pu < d <= final_q_accept_pu` | `bounded_q_limit_violation` | accepted, warning with bus, `d` in pu and MVAr |
   | violation | `d > final_q_accept_pu` | `remaining_pv_q_limit_violations` | not accepted |

   `power_flow.qlimits.hysteresis_pu` is the switching hysteresis, which
   the final check tolerates as well; `power_flow.qlimits.final_q_accept_pu`
   (default `auto`, twice the hysteresis; an explicit value must be
   `>= hysteresis_pu`) bounds what is still accepted; both at zero make
   the check strict. The count-based guard keys
   (`guard.accept_bounded_violations`, `guard.max_remaining_violations`)
   are not part of this classification. The check writes one line per bus
   (side, Q, limit, `d` in pu and MVAr, class) into the Q-limit block of
   the text report, `run.log`, the run metadata (`final_q_check_status`,
   `qlimits_disabled` when Q-limit handling is off, `not_evaluated` without
   a converged solution,
   `final_q_check_buses`, `final_q_check_max_dev_pu`) and the Web UI
   result page, next to the separate, voltage-side Q-V characteristic
   check.

---

## Interaction with the Rectangular Solver

The rectangular solver keeps a `bus_types` vector with values `:PQ`, `:PV`,
`:Slack`. A PV bus that hits its Q-limit is treated as PQ from the next
iteration on, so the residual switches from `(ΔP, ΔV)` to `(ΔP, ΔQ)`:

```julia
bus_types[b] = :PQ
```

Solver details: [Solver Guide](solver.md).

---

## Alternative to immediate PV→PQ switching: `qlimit_mode = :adjust_vset`

Besides the default `:switch_to_pq`, the rectangular solver supports a
controller-based strategy:

```julia
runpf!(net, 40, 1e-8, 1; qlimit_mode = :adjust_vset)
```

When a PV bus reaches a reactive limit, the voltage setpoint (`vm_pu`) is
first changed in steps of `vstep_pu` before the bus becomes PQ, which
avoids switching where a small voltage correction keeps Q inside its
limits. One regulating generator prosumer at the PV bus (only one
controller definition per bus) provides `vstep_pu` and optionally
`tap_steps_down` / `tap_steps_up` (maximum number of downward/upward
adjustments). Without a valid controller, or once the step limits are
exhausted, standard PV→PQ switching applies. Supported on the rectangular
solver path.

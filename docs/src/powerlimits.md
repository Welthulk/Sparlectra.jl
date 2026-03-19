# Reactive Power Limits (Q-Limits) and PV→PQ Switching

## Overview

Sparlectra supports per-bus reactive power limits for generators (PV buses).
When a PV bus violates its Q-limits, it is automatically converted to a PQ bus
during Newton-Raphson iteration (PV→PQ switching). Events are logged, and
optionally hysteresis and a "cooldown" period can be used for later switching
back.

## Data Structures in `Net`

The `Net` type stores all Q-limit-related data:

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

Typical source of these limits:

* Q-limits come from generators / prosumers (e.g. `qMin`, `qMax`), are
  converted to p.u. and aggregated per bus.
* A helper like `setQLimits!` / `buildQLimits!` pre-fills `net.qmin_pu` and
  `net.qmax_pu` before the power flow.

## Logging Q-Limit Hits

Q-limit events are represented by `QLimitEvent` and are stored in
`net.qLimitLog`:

```julia
Base.@kwdef struct QLimitEvent
    iter::Int
    bus::Int
    side::Symbol   # :min | :max
end
```

Support functions:

* `logQLimitHit!(net, iter, bus, side)`

  * Appends a `QLimitEvent` to `net.qLimitLog`
  * Updates `net.qLimitEvents[bus] = side`
* `lastQLimitIter(net, bus)`

  * Returns the last iteration index in which `bus` hit a Q-limit.
* `resetQLimitLog!(net)`

  * Clears both `qLimitLog` and `qLimitEvents`.

Human-readable printout:

```julia
printQLimitLog(net; sort_by = :iter, io = stdout)
```

Example:

```text
──────────────────────────────────
 Iteration │ Bus │ Side
──────────────────────────────────
         3 │   5 │ min
         5 │   3 │ max
──────────────────────────────────
Total events: 2
```

A convenience function for tests:

```julia
pv_hit_q_limit(net, ["B3", "B10"])
```

returns `true` if any of the listed PV buses appears in `net.qLimitEvents`.

---

## Conceptual Algorithm

During a Newton-Raphson calculation (e.g., `calcNewtonRaphson_withPVIdentity!`
or the rectangular variant), the Q-limit logic works at its core as follows:

1. **Solve one NR step** and obtain updated bus voltages `V` and bus powers
   `S = P + jQ`.
2. **Compute actual generator Q** per PV bus (from `S` or from bus data).

3. **Check limits** per bus index `b`:

   * Let `Q_b` be the current reactive injection at bus `b`.
   * Let `[Qmin_b, Qmax_b]` be the limits from `net.qmin_pu` and
     `net.qmax_pu`.
   * If

     ```math
     Q_b < Qmin_b \quad \text{or} \quad Q_b > Qmax_b
     ```

     then the PV bus has hit a limit.

4. **Clip Q and switch bus type**:

   * Set the generator’s `Q` to the violated limit:
     `Q_b := Qmin_b` or `Q_b := Qmax_b`.
   * Change the bus type from PV to PQ:

     ```julia
     setBusType!(net, b, "PQ")
     ```
   * Log the event:

     ```julia
     logQLimitHit!(net, iter, b, :min)  # or :max
     ```

5. **Optional: PV→PQ→PV re-enable (hysteresis + cooldown)**

   If `net.q_hyst_pu > 0` and `net.cooldown_iters > 0`, a bus may be switched
   back from PQ to PV only if:

   * Its Q is now well inside the interval with margin:

     ```math
     Qmin_b + q_\mathrm{hyst} < Q_b < Qmax_b - q_\mathrm{hyst}
     ```

   * And sufficiently many Newton iterations have passed since the last limit
     hit:

     ```math
     \text{iter} - \text{lastQLimitIter}(b) \ge \text{cooldown\_iters}
     ```

!!! note
    In many applications, switching back within the same power flow calculation
    is intentionally **disabled** (`q_hyst_pu = 0.0`, `cooldown_iters = 0`).

6. **Repeat** with the updated bus types and Q-values in the next NR iteration.

---

## Interaction with the Rectangular Solver

The rectangular solver uses a separate `bus_types` vector with values
`:PQ`, `:PV`, `:Slack`. When a PV bus hits its Q-limit, the conceptual solver
action is the same as in the polar formulation: from the next iteration
onward the bus is treated as PQ, so the residual switches from
`(ΔP, ΔV)` to `(ΔP, ΔQ)`.

```julia
bus_types[b] = :PQ
```

This means the Q-limit logic changes the *equation set* used by the Newton
step, not only a post-processing label.

For the numerical solver details behind this equation system, see
[Solver Guide](solver.md).

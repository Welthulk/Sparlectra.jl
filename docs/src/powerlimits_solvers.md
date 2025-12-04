## Reactive Power Limits (Q-Limits) and PV→PQ Switching

### Overview

Sparlectra supports per-bus reactive power limits for generators (PV buses).
When a PV bus violates its Q-limits, it is automatically converted to a PQ bus during Newton-Raphson iteration (PV→PQ switching).
Events are logged, and optionally hysteresis and a "cooldown" period can be used for later switching back.

### Data Structures in `Net`

The `Net` type stores all Q-limit–related data:

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

* Q-limits come from generators / prosumers (e.g. `qMin`, `qMax`), are converted to p.u. and aggregated per bus.
* A helper like `setQLimits!` / `buildQLimits!` pre-fills `net.qmin_pu` and `net.qmax_pu` before the power flow.

### Logging Q-Limit Hits

Q-limit events are represented by `QLimitEvent` and are stored in `net.qLimitLog`:

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

### Conceptual Algorithm (per NR iteration)

During a Newton-Raphson calculation (e.g., `calcNewtonRaphson_withPVIdentity!` or the rectangular variant), the Q-limit logic works at its core as follows:

1. **Solve one NR step** and obtain updated bus voltages `V` and bus powers `S = P + jQ`.
2. **Compute actual generator Q** per PV bus (from `S` or from bus data).

3. **Check limits** per bus index `b`:

   * Let `Q_b` be the current reactive injection at bus `b`.
   * Let `[Qmin_b, Qmax_b]` be the limits from `net.qmin_pu` and `net.qmax_pu`.
   * If

     ```math
     Q_b < Qmin_b \quad \text{or} \quad Q_b > Qmax_b
     ```

     then the PV-bus has hit a limit.

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
  If `net.q_hyst_pu > 0` and `net.cooldown_iters > 0`, a bus may be switched back from PQ to PV only if:

    * Its Q is now well inside the interval with margin:

      ```math
      Qmin_b + q_\mathrm{hyst} < Q_b < Qmax_b - q_\mathrm{hyst}
      ```
    * And sufficiently many Newton iterations have passed since the last limit hit:

      ```math
      \text{iter} - \text{lastQLimitIter}(b) \ge \text{cooldown\_iters}
      ```

 !!!Note: In many applications, switching back within the same power flow calculation is intentionally **disabled** (`q_hyst_pu = 0.0`, `cooldown_iters = 0`).

6. **Repeat** with the updated bus types and Q-values in the next NR iteration.

---

## Rectangular Complex-State Newton–Raphson (`jacobian_complex.jl`)

### Motivation and State Vector

The rectangular complex-state NR uses the complex bus voltages

```math
V_i = V_{r,i} + j V_{i,i}
```

as state variables instead of magnitudes and angles. The state vector is

```math
x = 
\begin{bmatrix}
V_r(\text{non-slack}) \\
V_i(\text{non-slack})
\end{bmatrix}
\in \mathbb{R}^{2(n-1)}.
```



The complex bus powers are computed as

```math
I = Y_\mathrm{bus} V, \qquad
S = V \odot \overline{I} = V \odot \overline{Y_\mathrm{bus} V}
```

and the specified injections as

```math
S_\mathrm{spec} = P_\mathrm{spec} + j Q_\mathrm{spec}.
```

### Bus Types and Residual Definition

The function

```julia
mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
    -> Vector{Float64}
```



builds the real-valued residual vector `F(V)`:

* For each **PQ bus** $i \ne \text{slack}$:

  ```math
  \Delta P_i = \Re(S_{\text{calc},i}) - \Re(S_{\text{spec},i})
  \\
  \Delta Q_i = \Im(S_{\text{calc},i}) - \Im(S_{\text{spec},i})
  ```

* For each **PV bus** $i \ne \text{slack}$:

  ```math
  \Delta P_i = \Re(S_{\text{calc},i}) - \Re(S_{\text{spec},i})
  \\
  \Delta V_i = |V_i| - V_\text{set,i}
  ```

The stacked residual vector has size (2(n-1)):

```math
F = 
\begin{bmatrix}
\Delta P_2 \\
\Delta Q/\Delta V_2 \\
\vdots \\
\Delta P_n \\
\Delta Q/\Delta V_n
\end{bmatrix}.
```

* `bus_types::Vector{Symbol}` uses `:PQ`, `:PV`, `:Slack`.
* `Vset::Vector{Float64}` holds the PV voltage setpoints (unused for PQ).


### Analytic Rectangular Newton Step

`complex_newton_step_rectangular` performs one Newton step using an analytic Jacobian:

```julia
complex_newton_step_rectangular(Ybus, V, S;
    slack_idx::Int = 1,
    damp::Float64  = 1.0,
    bus_types::Vector{Symbol},
    Vset::Vector{Float64},
)
```

Key steps:

1. Compute currents and powers:

   ```julia
   I     = Ybus * V
   S_calc = V .* conj.(I)
   ΔS    = S_calc .- S
   ```

2. Build a complex Jacobian `Jc` such that

   ```math
   \delta S = J_c \, \delta x_\mathbb{C}
   ```

   by perturbing `dVr` and `dVi` basis vectors (Wirtinger style).

3. Convert to a real Jacobian `J ∈ ℝ^{2n×2n}` for

   ```math
   F = \begin{bmatrix} \Re(\Delta S) \\ \Im(\Delta S) \end{bmatrix}.
   ```

4. Eliminate slack bus variables and rows (only non-slack buses contribute).

5. Solve

   ```math
   J_\text{red} \, \delta x_\text{red} = -F_\text{red}
   ```

   with fallback to `pinv(J)` in case of singularity.

6. Apply damping and update:

   ```julia
   V_new = (Vr .+ δVr) .+ im .* (Vi .+ δVi)
   ```

### Sparse Analytic Jacobian (Rectangular)

For performance, a dedicated sparse Jacobian builder exists:

```julia
build_rectangular_jacobian_pq_pv_sparse(
    Ybus::SparseMatrixCSC{ComplexF64},
    V::Vector{ComplexF64},
    bus_types::Vector{Symbol},
    Vset::Vector{Float64},
    slack_idx::Int,
) -> SparseMatrixCSC{Float64}
```



It uses identities like

```math
S_i(V) = V_i \, \overline{(Y V)_i}
```

and

```math
\frac{\partial S}{\partial V} = \operatorname{diag}(\overline{I})
                              + \operatorname{diag}(V)\, \overline{Y},
```

then transforms these to derivatives w.r.t. (V_r) and (V_i), and finally to (\Delta P, \Delta Q, \Delta V). The sparsity pattern follows `Ybus`.

### High-Level Solver: `run_complex_nr_rectangular`

The main iteration routine is:

```julia
run_complex_nr_rectangular(
    Ybus,
    V0,
    S;
    slack_idx::Int = 1,
    maxiter::Int = 20,
    tol::Float64 = 1e-8,
    verbose::Bool = false,
    damp::Float64 = 1.0,
    bus_types::Vector{Symbol},
    Vset::Vector{Float64},
    use_fd::Bool = false,
    use_sparse::Bool = false,
) -> (V::Vector{ComplexF64}, converged::Bool, iters::Int, history::Vector{Float64})
```



* `use_fd = true` → use finite-difference Jacobian (`complex_newton_step_rectangular_fd`).
* `use_sparse = true` → use sparse analytic Jacobian where available.
* `history` stores the max mismatch per iteration.

Integration with `Net` happens über Hilfsfunktionen:

* `initial_Vrect_from_net(net) -> (V0, slack_idx)`
* `build_S_from_net(net) -> S`
* `update_net_voltages_from_complex!(net, V)`
* `run_complex_nr_rectangular_for_net!(net, ...)` (High-level wrapper, ruft intern `run_complex_nr_rectangular` auf und verbindet alles mit `runpf!`).

`runpf!` verwendet die Methode `:rectangular` und Optionen `opt_fd` / `opt_sparse`:

```julia
ite, status = runpf!(net, maxIte, tol, verbose;
    method    = :rectangular,
    opt_fd    = true,    # FD-Jacobian
    opt_sparse = true,   # sparse analytic J where applicable
)
```

---

## Finite-Difference Rectangular Jacobian (`jacobian_fd.jl`)

### Purpose

The file `jacobian_fd.jl` provides a reference / fallback implementation of the Newton step for the rectangular complex-state NR using finite differences. It is mainly useful for:

* Prototyping and validation of the analytic Jacobian
* Debugging difficult cases
* Small to medium test networks (the FD Jacobian is dense and relatively expensive)

### Function Signature

```julia
complex_newton_step_rectangular_fd(
    Ybus,
    V,
    S;
    slack_idx::Int       = 1,
    damp::Float64        = 1.0,
    h::Float64           = 1e-6,
    bus_types::Vector{Symbol},
    Vset::Vector{Float64},
) -> Vector{ComplexF64}
```



Arguments:

* `Ybus` — bus admittance matrix (n×n, complex)
* `V` — current complex voltages (length n)
* `S` — specified complex injections (length n)
* `slack_idx` — slack bus index
* `damp` — damping factor for the Newton step
* `h` — FD perturbation step (typically 1e-6)
* `bus_types` — `:PQ`, `:PV`, `:Slack`
* `Vset` — PV setpoints

### Internal Workflow

1. **Base mismatch**

   ```julia
   F0 = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
   m  = length(F0)      # expected = 2*(n-1)
   ```

2. **Non-slack variables**

   Only non-slack buses are state variables:

   ```julia
   non_slack = collect(1:n)
   deleteat!(non_slack, slack_idx)
   nvar = 2 * (n - 1)   # Vr(non-slack) and Vi(non-slack)
   ```

3. **Jacobian by finite differences**

   * First, perturb each real part (V_{r,k}):

     ```julia
     for (col_idx, bus) in enumerate(non_slack)
         V_pert = copy(V)
         V_pert[bus] = ComplexF64(real(V[bus]) + h, imag(V[bus]))
         Fp = mismatch_rectangular(...)
         J[:, col_idx] .= (Fp .- F0) ./ h
     end
     ```

   * Then, perturb each imaginary part (V_{i,k}):

     ```julia
     for (offset, bus) in enumerate(non_slack)
         col_idx = (n - 1) + offset
         V_pert = copy(V)
         V_pert[bus] = ComplexF64(real(V[bus]), imag(V[bus]) + h)
         Fp = mismatch_rectangular(...)
         J[:, col_idx] .= (Fp .- F0) ./ h
     end
     ```

4. **Solve linear system**

   ```julia
   δx = J \ (-F0)
   ```

   with fallback to `pinv(J)` if `J` is singular. Then apply damping

   ```julia
   δx .*= damp
   ```

   and map back to updates in `Vr` and `Vi` for all non-slack buses.

5. **Return updated voltages**

   ```julia
   V_new = ComplexF64.(Vr_new, Vi_new)
   ```

`run_complex_nr_rectangular` selects this FD step when `use_fd = true`.

---

## PV→PQ Switching in the Rectangular Solver

### Bus-Type Handling via `bus_types`

The rectangular solver mirrors the classical polar NR logic by encoding bus types in a separate vector:

```julia
bus_types[i] ∈ (:PQ, :PV, :Slack)
```

* `mismatch_rectangular` interprets each bus based on `bus_types[i]` and constructs either `(ΔP, ΔQ)` or `(ΔP, ΔV)` residuals (or none for the slack).
* The Jacobians (analytic and FD) implicitly assume this layout and generate matching rows.

### Integration with Q-Limits

When a PV bus hits its Q-limit, the following updates take place conceptually:

1. **Mark bus as PQ**
   In the rectangular formulation, the NR step uses

   ```julia
   bus_types[b] = :PQ
   ```

   instead of `:PV`. From the next iteration onward, the residual for this bus becomes `(ΔP, ΔQ)` instead of `(ΔP, ΔV)`.

2. **Fix generator Q at the limit**
   The specified injection vector `S` is updated such that the bus Q equals the active limit (`Qmin` or `Qmax`):

   ```math
   Q_{\text{spec},b} := Q_{\text{limit},b}.
   ```

   Thus `ΔQ_b` reflects deviations from this clipped reactive power.

3. **Log the event and enforce cooldown/hysteresis**
   The same `QLimitEvent` / `qLimitLog` / `qLimitEvents` mechanism applies as for the polar solver, allowing consistent diagnostics across formulations.

4. **Optional re-enable (PQ→PV)**
   If hysteresis and cooldown are enabled (`q_hyst_pu > 0`, `cooldown_iters > 0`), a bus that was switched to PQ can be re-enabled as PV:

   ```julia
   bus_types[b] = :PV
   ```

   once:

   * Its Q is sufficiently inside the allowed band (with hysteresis),
   * Enough iterations have passed since the last limit hit.

   
   !!! Hint: In many practical applications, however, a bus remains in the state it was switched to during a single power flow calculation (i.e., no PQ→PV return during the same calculation).

Gerne – hier ist **nur der fehlende Teil**, den du direkt *unten an deinen bestehenden workshop.md anhängen* kannst.
Er ergänzt:

1. **Ausgabe der Generatorleistungen** nach dem Power-Flow
2. **Waterfilling-Methode** zur fairen Q-Verteilung bei mehreren Generatoren pro Bus
3. **Einbettung in typische Sparlectra-Workflows** (rectangular & polar)

Formulierung ist bewusst **knapp, technisch, nachvollziehbar**, wie für einen Workshop.

---

# Generator Power Output and Waterfilling Allocation

## Overview

After a successful power-flow calculation (polar or rectangular), Sparlectra can output **per-generator active and reactive power**.
Since a bus may host *multiple* generators (ProSumers), and their combined injection must match the solved bus injection, a **post-processing step** allocates the solved values to individual units.

This allocation must handle:

* PV buses that stayed PV
* PV buses that were switched to PQ due to Q-limit violations
* Multiple generators per bus
* At least one generator possibly being Q-limited
* Fair distribution of remaining Q using a **waterfilling** scheme

The results are written into the `net.prosumpVec` entries and exported via `results.jl`.

---

## Generator Output After Power Flow

Once `V` and `S_calc` are known for all buses, Sparlectra computes:

```math
S_{\text{bus}} = P_{\text{bus}} + j Q_{\text{bus}}
```

and distributes this to all generators attached to that bus:

1. Read all ProSumers `g₁, g₂, …, g_k` at bus `b`.
2. Sum their nominal power references (or equal-share if no reference available).
3. Allocate:

```math
P_g = \alpha_g \, P_{\text{bus}}
```

equally or proportionally (depending on configuration).
4. Allocate Q using the Waterfilling algorithm described below.

Finally, each generator stores:

```julia
prosumer.P = Pg          # real power allocated
prosumer.Q = Qg          # reactive power allocated
prosumer.atQLimit = true | false
prosumer.qLimitSide = :min | :max | :none
```

These values appear in `results.jl` when printing generator results:

```
Generator results:
  Bus B5, Gen#1: P=…, Q=…, atQLimit=false
  Bus B5, Gen#2: P=…, Q=…, atQLimit=true (max)
```

---

## Waterfilling Reactive Power Distribution

If a bus hosts multiple generators *and* its total solved reactive power is:

```math
Q_{\text{bus}} = \sum_i Q_i
```

and at least one generator has a Q-limit, the Waterfilling method assigns Q fairly under constraints.

### Inputs

* Generators `g₁, …, g_k` with limits:

```math
Q_{i,\min} \le Q_i \le Q_{i,\max}
```

* Total bus reactive power `Q_bus` (from the NR solution)

### Goal

Find generator values:

```math
Q_1, Q_2, …, Q_k
```

such that:

1. All stay inside their limits
2. The sum matches the bus:

```math
\sum_i Q_i = Q_{\text{bus}}
```

3. The allocation is “fair”: all non-limited units are increased/decreased uniformly.

### Algorithm (conceptual)

1. Start with all generators marked **free** (not at limit).
2. Compute an equal provisional share:

```math
Q_i = Q_{\text{bus}} / k
```

3. Check limits for all generators:

   * If `Q_i > Q_{i,max}` → fix `Q_i = Q_{i,max}`, mark limited
   * If `Q_i < Q_{i,min}` → fix `Q_i = Q_{i,min}`, mark limited
4. Remove limited generators from the free pool.
5. Redistribute the remaining Q among the free generators:

```math
Q_{\text{remaining}} = Q_{\text{bus}} - \sum_{\text{limited}} Q_i
```

```math
Q_{\text{free}} = Q_{\text{remaining}} / k_{\text{free}}
```

6. Repeat steps 3–5 until:

   * All free generators respect their limits, and
   * The sum constraint is satisfied.

### Result

* Limited generators sit *exactly* at their Q-limit.
* Remaining required Q is spread uniformly among the others.
* The allocator guarantees feasibility when the bus Q is within the aggregate limits.

Each generator `g_i` receives flags:

```julia
g_i.atQLimit     = (Q_i == Q_i_min || Q_i == Q_i_max)
g_i.qLimitSide   = :min / :max / :none
```

---

## Integration Into The Power-Flow Workflow

### For Polar NR (`calcNewtonRaphson_withPVIdentity!`)

1. NR iterations may switch PV→PQ based on Q-limits
2. After convergence:

   * Compute bus powers
   * Apply waterfilling per bus
   * Write back generator values
   * Export via `results.jl`

### For Rectangular NR (`run_complex_nr_rectangular_for_net!`)

1. `bus_types` encodes PV/PQ
2. Q-limits handled identically (PV→PQ during iteration)
3. After convergence:

   * Compute `S_calc = V .* conj(Ybus*V)`
   * Allocate generator P and Q via waterfilling
   * Mark Q-limited generators
   * Export via `results.jl`

---

## Example Output Section in `workshop.md`

```
### Generator Outputs After Power Flow

After running the PF, Sparlectra prints a generator table such as:

Generator results:
─────────────────────────────────────────────
 Bus    Gen#   P (MW)     Q (MVAr)   Status
─────────────────────────────────────────────
 B5      1     12.40       3.10      ok
 B5      2     12.40       5.00      Q-max limit
 B9      1      8.20      -1.50      ok
─────────────────────────────────────────────
```

The Q-values follow the Waterfilling allocation. Generators at limits are marked accordingly.



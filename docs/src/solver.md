# Solver Guide

This page collects numerical notes for the supported internal power-flow path:
the sparse rectangular complex-state Newton–Raphson solver.

The rectangular power-flow implementation is split into focused layers:

* `src/powerflow_rectangular/rectangular_core_equations.jl` contains the core residual/equation helpers, including `mismatch_rectangular`.
* `src/powerflow_rectangular/rectangular_jacobian_builders.jl` contains the analytic rectangular Jacobian builders.
* `src/powerflow_rectangular/rectangular_newton_step.jl` contains Newton-step update and damping helpers.
* `src/powerflow_rectangular/rectangular_standalone_solver.jl` contains `run_complex_nr_rectangular`, the standalone array-level Newton driver.
* `src/powerflow_rectangular/rectangular_network_solver.jl` contains `runpf_rectangular!`, the network-integrated rectangular power-flow entry point and orchestration glue.

## Rectangular Complex-State Newton–Raphson (`powerflow_rectangular/rectangular_network_solver.jl`)

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

The stacked residual vector has size `2(n-1)` and matches the state dimension.

### Analytic Rectangular Newton Step

`complex_newton_step_rectangular` performs one Newton step using an analytic
Jacobian:

```julia
complex_newton_step_rectangular(Ybus, V, S;
    slack_idx::Int = 1,
    damp::Float64  = 1.0,
    autodamp::Bool = false,
    autodamp_min::Float64 = 0.05,
    bus_types::Vector{Symbol},
    Vset::Vector{Float64},
)
```

Key steps:

1. Compute currents and powers.
2. Form the residual `F(V)`.
3. Assemble the Jacobian `J = ∂F/∂x`.
4. Solve

   ```math
   J(x_k)\,\Delta x_k = -F(x_k)
   ```

5. Update the state with optional fixed damping or automatic damping.

This is the standard Newton linearization of a nonlinear root-finding problem.

### Automatic Rectangular Newton Damping

For difficult flat-start studies, enable residual-based backtracking:

```julia
runpf!(net, 60, 1e-8, 1;
    opt_flatstart = true,
    autodamp = true,
    autodamp_min = 0.05,
)
```

The rectangular solver first computes the ordinary Newton correction from
`J * Δx = -F`. With `autodamp = true`, it tests trial steps from `damp` down to
`autodamp_min` by halving the step length. The first trial that reduces the
maximum absolute mismatch is accepted. If no trial improves the mismatch, the
solver continues with the best finite conservative trial so that the next
iteration can rebuild the Jacobian and active-set state.

This is a line-search style damping strategy. It does not change the power-flow
model, bus equations, Q-limit logic, or Jacobian formulation. It only chooses the
scalar step length for the already computed rectangular Newton correction.

---

### Start Projection for Difficult Seeds

#### DC-angle flat-start background

A conventional AC flat start initializes most bus voltages near $1.0\,\mathrm{pu}$
and all non-slack voltage angles near $0^\circ$. That seed is simple and
reproducible, but it ignores the active-power flow pattern that is already
encoded in the network topology, branch reactances, transformer phase shifts,
and specified injections. In large or heavily phase-shifted MATPOWER cases, the
first Newton step may therefore start far away from the physically relevant
angle branch.

The DC-angle seed uses the active-power part of the power-flow model as a
linearized predictor for voltage angles. Under the usual high-voltage,
small-angle assumptions, voltage magnitudes are held near nominal values,
resistance and reactive-power coupling are neglected, and active-power transfer
across a branch is approximated by the angle difference divided by branch
reactance. This gives a sparse linear system of the form

```math
B'\theta = P
```

where $P$ is the net active-power injection vector and $B'$ is assembled from
the branch susceptance structure. The slack angle fixes the reference, and the
resulting angles are clipped by `start_projection_dc_angle_limit_deg` before
being used as a Newton seed. This is not a replacement for the AC solve; it is
only an initialization step that preserves the full AC equations, Q-limit logic,
and rectangular Newton formulation used by the main run.

Blended starts combine this DC-angle predictor with the stored MATPOWER `VM`/`VA`
data or the raw flat start. The projection scans the requested candidates and
keeps the one with the smallest rectangular mismatch, which helps avoid wrong
low-voltage or wrong-angle branches without changing the final convergence
criteria.

The rectangular solver can project the initial voltage before Newton iterations:

```julia
runpf!(net, 60, 1e-8, 1;
    opt_flatstart = true,
    start_projection = true,
    start_projection_try_dc_start = true,
    start_projection_try_blend_scan = true,
    start_projection_blend_lambdas = [0.25, 0.5, 0.75],
    start_projection_dc_angle_limit_deg = 60.0,
)
```

The projection is inspired by the other workflows. It sanitizes the raw seed,
optionally builds a DC-angle approximation from active-power injections and the
Y-bus off-diagonal susceptances, and optionally scans convex angle/magnitude
blends between the raw seed and the DC start. Sparlectra picks the candidate with
the lowest rectangular mismatch and then runs the ordinary solver.

The same projection options are also available through `buildPfModel` and
`runpf_external!`, so external solvers can receive the projected
`model.V0` without changing their own solver implementation.

---


## Solver-Specific Interaction with Power Limits

When a PV bus hits a Q-limit, the solver does **not** merely clamp a reported
reactive power value. It changes the equation type for that bus:

* from `(ΔP, ΔV)` at a PV bus
* to `(ΔP, ΔQ)` at a PQ bus.

In the rectangular solver this is represented through the `bus_types` vector.

The start of PV→PQ switching can be controlled for difficult cases:

```julia
runpf!(net, 60, 1e-8, 1;
    qlimit_start_iter = 4,
    qlimit_start_mode = :iteration,
)
```

Use `qlimit_start_mode = :auto_q_delta` to wait until the PV reactive-power
requests have stabilized. The threshold is `qlimit_auto_q_delta_pu` in p.u.
The hybrid `:iteration_or_auto` mode starts when either the configured iteration
or the reactive-power stabilization criterion is reached.

For the operational / data-model side of Q-limit handling, see
[Powerlimits Guide](powerlimits.md).

---

## Voltage-dependent P(U)/Q(U) controls

The rectangular solver also supports state-dependent specified injections via
`PUController` and `QUController` attached to prosumers. In that case, the
specified power vector is re-evaluated each Newton iteration as a function of
local `|V|`, and local chain-rule terms are added to the Jacobian.

For full derivation, formulas, API usage, and output semantics (`Type` vs
`Control`), see [Voltage-dependent Control](voltage_dependent_control.md).

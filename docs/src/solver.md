# Solver Guide

This page collects the numerical solver notes that were previously mixed with
power-limit documentation. It focuses on *why* the finite-difference (FD)
methods work and how they relate to the analytic Jacobian-based methods.

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

The stacked residual vector has size `2(n-1)` and matches the state dimension.

### Analytic Rectangular Newton Step

`complex_newton_step_rectangular` performs one Newton step using an analytic
Jacobian:

```julia
complex_newton_step_rectangular(Ybus, V, S;
    slack_idx::Int = 1,
    damp::Float64  = 1.0,
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

5. Update the state with optional damping.

This is the standard Newton linearization of a nonlinear root-finding problem.

---

## Finite-Difference Rectangular Jacobian (`jacobian_fd.jl`)

### Purpose

The file `jacobian_fd.jl` provides a reference / fallback implementation of the
Newton step for the rectangular complex-state NR using finite differences. It
is mainly useful for:

* prototyping and validation of the analytic Jacobian,
* debugging difficult cases,
* small to medium test networks.

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

### Why the FD Method Works

The FD method works because Newton's method only needs a **local linear model**
of the residual map near the current iterate. If

```math
F : \mathbb{R}^n \rightarrow \mathbb{R}^n
```

is differentiable, then for a small perturbation `δx`:

```math
F(x + \delta x) \approx F(x) + J(x)\,\delta x,
```

where `J(x)` is the Jacobian. Newton uses this first-order model and chooses
`δx` such that the linearized residual becomes zero:

```math
J(x)\,\delta x = -F(x).
```

If we do not assemble `J(x)` analytically, we can approximate each column by
perturbing one state variable at a time:

```math
\frac{\partial F}{\partial x_k}(x)
\approx
\frac{F(x + h e_k) - F(x)}{h},
```

with `e_k` the k-th unit vector and `h` a small step size. This is exactly what
the FD implementation does. The approximation is first-order accurate in `h`,
so when `h` is small enough, the FD Jacobian is a good surrogate for the true
Jacobian.

In short:

* the **physics** stays the same because `F` is still the same mismatch
  function,
* only the **derivative evaluation** changes,
* and Newton still solves the same linearized correction equation.

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

   * First, perturb each real part `V_r`:

     ```julia
     for (col_idx, bus) in enumerate(non_slack)
         V_pert = copy(V)
         V_pert[bus] = ComplexF64(real(V[bus]) + h, imag(V[bus]))
         Fp = mismatch_rectangular(...)
         J[:, col_idx] .= (Fp .- F0) ./ h
     end
     ```

   * Then, perturb each imaginary part `V_i`:

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

   with fallback to a pseudo-inverse path in the solver utilities if needed.
   Then apply damping:

   ```julia
   δx .*= damp
   ```

5. **Return updated voltages**

   ```julia
   V_new = ComplexF64.(Vr_new, Vi_new)
   ```

`run_complex_nr_rectangular` selects this FD step when `use_fd = true`.

### Practical Interpretation

The FD solver is therefore not a different load-flow model. It is the same NR
equation system with a numerically approximated Jacobian. This is why it is
useful as a validation path for the analytic Jacobian: both methods should
follow the same local Newton geometry, up to FD truncation and step-size
effects.

---

## Solver-Specific Interaction with Power Limits

When a PV bus hits a Q-limit, the solver does **not** merely clamp a reported
reactive power value. It changes the equation type for that bus:

* from `(ΔP, ΔV)` at a PV bus
* to `(ΔP, ΔQ)` at a PQ bus.

In the rectangular solver this is represented through the `bus_types` vector.

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
